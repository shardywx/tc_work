; Script to plot basic variables (vorticity, horizontal wind, vertical velocity)
; on model levels for a given tropical cyclone

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/rad_circ.ncl"

; ncl dat=\"02T12\" ens0=\"em11\" opt=\"x11\" dist=3.0 ts0=12 tf0=18 ls0=5 lf0=10 
; calc=1 full=1 clr=1 tc_basic_xy_cyl.ncl

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'ens0'  = ensemble simulation (em00-em11)
; 'opt'   = output file format ("pdf" or "x11")
; 'dist'  = size of box following storm (degrees)
; 'ts0'	  = start time index 
; 'tf0'   = end time index
; 'ls0'	  = start model level
; 'lf0'	  = end model level
; 'calc'  = use single storm centre (1) or calculate on each level (0)
; 'full'  = plot either 2 panels (1) or 1 (0)
; 'clr'	  = colour table --> "ncl_default" (1), "bdr_extra" (2) or "amwg256" (3)

begin

;==============================
; Prelash (reading in files)
;==============================

; Find 'pc' and 'pd' files (model height level data)
  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0
  fili_prefix = "20160702T1200Z_NPTK_4p4_L80_ra1t_"+ens0
  fili_c1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pc*.nc")
  fili_c      = diri+"/"+fili_c1
  fili_d1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pd*.nc")
  fili_d      = diri+"/"+fili_d1

; Calculate number of times in each file using built-in function ('nc_times')
  numINFO_c   = nc_times(fili_c)
  numINFO_d   = nc_times(fili_d)
  numINPUT_c  = numINFO_c[0]            ; Number of times in 'pc' stream
  time_c      = numINFO_c[1]            ; Array of times from 'pc' stream file
  numINPUT_d  = numINFO_d[0]            ; Number of times in 'pd' stream
  time_d      = numINFO_d[1]            ; Array of times from 'pd' stream file

; Calculate storm motion using built-in function ('st_rm')
; Argument 1  = ensemble simulation
; Argument 2  = path to input files
; Argument 3  = size of box following storm (degrees)

; Call external function 'st_rm'
  storm_rel = st_rm(dat, ens0, "$ar/text/", dist, 0)
  u_cyc     = storm_rel[0]                 ; Zonal wind
  v_cyc     = storm_rel[1]                 ; Meridional wind
  vel_cyc   = storm_rel[2]                 ; Vector wind
  numTIMES  = storm_rel[3]                 ; Number of times in file
  lat_arr   = storm_rel[4]                 ; Storm-following domain subset (lat)
  lon_arr   = storm_rel[5]                 ; Storm-following domain subset (lon)
  centre    = storm_rel[6]                 ; Storm centre at each hour (lat/lon)

;==========================================================================
; Create date/time string arrays for all times in file (use 'pc' stream)
;==========================================================================

  month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
                    "Jul","Aug","Sep","Oct","Nov","Dec"/)

  utc_date = cd_calendar(time_c, 0)

  year   = tointeger(utc_date(:,0))
  month  = tointeger(utc_date(:,1))
  day    = tointeger(utc_date(:,2))
  hour   = tointeger(utc_date(:,3))
  minute = tointeger(utc_date(:,4))
  second = utc_date(:,5)

; Correct for errors in the code (round up value of hour when minutes = 59)
  do ct = 0, numINPUT_c-1
   if (minute(ct).gt.30) then
     hour(ct) = hour(ct)+1
   end if
  end do
  delete(ct)

; Create arrays to hold strings for output plots
  date_str  = new(numINPUT_c,string)
  time_str  = new(numINPUT_c,string)
  time_arr  = new(numINPUT_c,string)
  title_arr = new(numINPUT_c,string)

  date_str = sprinti("%0.2iUTC ", hour) + \
             sprinti("%0.2i ", day) \
           + month_abbr(month)
  time_str = sprinti("%0.2i UTC ", hour) + \
             sprinti("%0.2i ", day) \
           + month_abbr(month)
  time_arr = sprinti("%0.2i", day) + \
             month_abbr(month) + \
             "_" + sprinti("%0.2iZ", hour)

;==================================================================================
;;;;;;;;;;;;;;;;;;;;;;;;;;; Start of main part of script ;;;;;;;;;;;;;;;;;;;;;;;;;;
;==================================================================================

  print_clock("Starting main part of script...")

;========================================
; Define cylindrical coordinate arrays
;========================================

  thetas     = new(73,float)
  pii        = 3.14159265
  radii      = fspan(0.0,1.5,31)     ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,73)   ; Array of azimuth angles (0, 5, 10, etc)
  thetas_rad = (pii/180)*thetas      ; Degrees to radians
  rad_size   = dimsizes(radii)-1     ; Index for later calculation
  num_l      = 59                    ; Number of model levels

; Define constants
  rd    = 287.0   ; Gas constant for dry air (J/kg/K^2)
  cp0   = 1004.0  ; Specific heat of dry air at constant pressure (J/kg/K)
  kp    = 0.286   ; For Exner function calculations (0.286)
  pref  = 1000.0  ; Reference pressure p0 (1000 hPa)

; Array to hold all 'xpos' and 'ypos' position arrays
  xpos_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  ypos_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)

; Define 3D arrays for mean quantities

; Mean zonal wind {u}
  mean_u0          = new((/num_l,numTIMES,dimsizes(radii)/),float)
  mean_u0!0        = "lev"
  mean_u0!1        = "time"
  mean_u0!2        = "rad"
  mean_u0@description = "Azimuthally averaged zonal wind"
  mean_u0@units    = "m s~S~-1~N~"

; Mean meridional wind {v}
  mean_v0          = new((/num_l,numTIMES,dimsizes(radii)/),float)
  mean_v0!0        = "lev"
  mean_v0!1        = "time"
  mean_v0!2        = "rad"
  mean_v0@description = "Azimuthally averaged meridional wind"
  mean_v0@units    = "m s~S~-1~N~"

; Mean tangential wind {vtan}
  mean_vtan0       = new((/num_l,numTIMES,dimsizes(radii)/),float)
  mean_vtan0!0     = "lev"
  mean_vtan0!1     = "time"
  mean_vtan0!2     = "rad"
  mean_vtan0@description = "Azimuthally averaged tangential wind"
  mean_vtan0@units = "m s~S~-1~N~"

; Mean radial wind {vrad}
  mean_vrad0       = new((/num_l,numTIMES,dimsizes(radii)/),float)
  mean_vrad0!0     = "lev"
  mean_vrad0!1     = "time"
  mean_vrad0!2     = "rad"
  mean_vrad0@description = "Azimuthally averaged radial wind"
  mean_vrad0@units = "m s~S~-1~N~"

; Mean vertical velocity {w}
  mean_vvel0       = new((/num_l,numTIMES,dimsizes(radii)/),float)
  mean_vvel0!0     = "lev"
  mean_vvel0!1     = "time"
  mean_vvel0!2     = "rad"
  mean_vvel0@description = "Azimuthally averaged vertical velocity"
  mean_vvel0@units = "m s~S~-1~N~"

; Mean relative vorticity {zeta}
  mean_vort0       = new((/num_l,numTIMES,dimsizes(radii)/),float)
  mean_vort0!0     = "lev"
  mean_vort0!1     = "time"
  mean_vort0!2     = "rad"
  mean_vort0@description = "Azimuthally averaged relative vorticity"
  mean_vort0@units = "s~S~-1~N~"

; Mean absolute vorticity {zeta + f0}
  mean_avo0        = new((/num_l,numTIMES,dimsizes(radii)/),float)
  mean_avo0!0      = "lev"
  mean_avo0!1      = "time"
  mean_avo0!2      = "rad"
  mean_avo0@description = "Azimuthally averaged absolute vorticity"
  mean_avo0@units  = "s~S~-1~N~"

;=========================
; Start multiple loops
;=========================

    print_clock("Working on ensemble member: "+ens0)

  ; 12Z 2 Jul [00], 18Z 2 Jul [01], 00Z 3 Jul [02], 06Z 3 Jul [03], 12Z 3 Jul [04],
  ; 18Z 3 Jul [05], 00Z 4 Jul [06], 06Z 4 Jul [07], 12Z 4 Jul [08], 18Z 4 Jul [09],
  ; 00Z 5 Jul [10], 06Z 5 Jul [11], 12Z 5 Jul [12], 18Z 5 Jul [13], 00Z 6 Jul [14],
  ; 06Z 6 Jul [15], 12Z 6 Jul [16], 18Z 6 Jul [17], 00Z 7 Jul [18], 06Z 7 Jul [19],
  ; 12Z 7 Jul [20]

  ; ct = 0   [00], ct = 5  [01], ct = 11  [02], ct = 17  [03], ct = 23  [04],
  ; ct = 29  [05], ct = 35 [06], ct = 41  [07], ct = 47  [08], ct = 53  [09],
  ; ct = 59  [10], ct = 65 [11], ct = 71  [12], ct = 77  [13], ct = 83  [14],
  ; ct = 89  [15], ct = 95 [16], ct = 101 [17], ct = 107 [18], ct = 113 [19],
  ; ct = 119 [20]

    numTIMES    = 120                       ; Total number of times in all files [36]
    ct          = ts0                       ; Counter variable
    llbox       = toint(dist*50)            ; Size of lat/lon arrays (see line below)
    dsize       = (/llbox,llbox/)           ; Size of lat/lon arrays (depends on 'dist')
    d0          = dsize(0)
    d1          = dsize(1)

    setvalues NhlGetWorkspaceObjectId
      "wsMaximumSize" : 1000000000
    end setvalues

;===============================
; Start loop over input files
;===============================

    c  = addfile(fili_c,"r")             ; Read in 'pc' stream [u,v,vort]
    d  = addfile(fili_d,"r")             ; Read in 'pd' stream [w,p,t,theta,kmh,tau]

;==================================
; Get the variables we will need
;==================================

   do it = ts0, tf0             ; Loop over times in file

    print("Working on time: "+time_str(it))
    itt = it + 1
    title_arr(it) = "Valid at "+time_str(it)+" (T+"+itt+")"

;=====================================================
; Read in variables from 'pc' stream (model levels)
;=====================================================

    t0    = lat_arr(cn,it,0)
    t1    = lat_arr(cn,it,1)
    n0    = lon_arr(cn,it,0)
    n1    = lon_arr(cn,it,1)
    print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

  ; Read in basic variables from 'pc' stream
    lon1C  = c->longitude_1({n0:n1})              ; longitude ['d1' grid points]
    lat1C  = c->latitude_1({t0:t1})               ; latitude  ['d0' grid points]
    lonC   = c->longitude({n0:n1})                ; longitude ['d1' grid points]
    latC   = c->latitude({t0:t1})                 ; latitude  ['d0' grid points]
    hybC   = c->hybrid_ht(:)                      ; 63 model ('theta') levels

    u_varname = "x-wind"
    v_varname = "y-wind"

  ; Read in horizontal wind components
  ; Both variables below -- [hybC | 63] * [latC | 150] * [lonC | 150]
    u  = c->$u_varname$(it,:,{t0:t1},{n0:n1}) ; Zonal wind
    v  = c->$v_varname$(it,:,{t0:t1},{n0:n1}) ; Meridional wind

  ; Calculate relative vorticity using centered finite differences
  ; Option '2' --> boundary points estimated using one-sided difference scheme
    vort = uv2vr_cfd(u, v, latC, lonC, 2)
    copy_VarCoords(u,vort)
    vort@units = "s~S~-1~N~"
    vort@name = "Relative vorticity on model levels"

  ; Coriolis parameter and absolute vorticity
    f1   = coriolis_param(latC)
    f0   = conform_dims(dimsizes(vort),f1,1)
    copy_VarMeta(vort,f0)
    copy_VarCoords(vort,f0)
    f0@description = "Coriolis parameter"
    f0@name        = "Coriolis parameter"

    avo  = vort + f0
    copy_VarMeta(vort,avo)
    copy_VarAtts(vort,avo)
    copy_VarCoords(vort,avo)
    avo@description = "Absolute vorticity"
    avo@name        = "Absolute vorticity on model levels"
    avo@long_name   = "Absolute vorticity"

;==================================================================
; Read in variables from 'pd' stream (model levels)
; B/c of differences in file size, access variables using 'it+2'
;==================================================================

  ; Read in basic variables from 'pd' stream
    lon1D  = d->longitude_1({n0:n1})                ; longitude [xxx grid points]
    lat1D  = d->latitude_1({t0:t1})                 ; latitude  [xxx grid points]
    lonD   = d->longitude({n0:n1})                  ; longitude [xxx grid points]
    latD   = d->latitude({t0:t1})                   ; latitude  [xxx grid points]
    hybD   = d->hybrid_ht(:)                        ; 63 model ('rho') levels

  ; Read in temperature, pressure, vertical velocity and geopotential height
  ; All variables below -- [hybD | 63] * [latD | 150] * [lon1D | 151]
    t      = d->temp(it-2,:,{t0:t1},{n0:n1})        ; Temperature (K)
    p      = d->p(it-2,:,{t0:t1},{n0:n1})           ; Pressure (Pa)
    w      = d->dz_dt(it-2,:,{t0:t1},{n0:n1})       ; Vertical velocity (m s-1)
    z      = d->ht(it-2,:,{t0:t1},{n0:n1})          ; Geopotential height (m)

  ; Convert units of pressure and geopotential height before continuing
    p  = p / 100
    z  = z / 10

  ; Add metadata
    p@units = "hPa"
    z@units = "dam"

    dx      = lonD(1) - lonD(0) ; Grid spacing (longitude)
    dy      = latD(1) - latD(0) ; Grid spacing (latitude)

;==========================
; Loop over model levels
;==========================

    do ilev = ls0, lf0 ; Loop over selected model levels

     hy = sprintf("%0.0f",hybD(ilev))
     print("Working on model level " + ilev+ " ("+hy+" m AGL)")

   ; Read in other variables on model levels
     u_plane                = u(ilev,:,:)          ; Zonal wind
     v_plane                = v(ilev,:,:)          ; Meridional wind
     w_plane                = w(ilev,:,:)          ; Vertical velocity
     vort_plane             = vort(ilev,:,:)       ; Relative vorticity
     avo_plane              = avo(ilev,:,:)        ; Absolute vorticity

   ; Calculate relative vorticity maximum on model level
     vort_max               = max(vort_plane)
     vort_smth              = smth9_Wrap(vort_plane, 0.5, 0.5, True)
     vort_max_smth          = max(vort_smth)

   ; Find domain index corresponding to vorticity maximum 
     dims                   = dimsizes(vort_plane)
     vort1d                 = ndtooned(vort_plane)
     inds                   = ind_resolve(maxind(vort1d),dims)
     vort1d_sm              = ndtooned(vort_smth)
     inds_sm                = ind_resolve(maxind(vort1d_sm),dims)

   ;==================================================================
   ; Either calculate the centre, or read in pre-calculated values
   ;==================================================================

   if (calc .eq. 0) then

   ; Find latitude and longitude of storm at each time and on each model level
     lat_max                = t0 + (dy * inds(0,0))
     lon_max                = n0 + (dx * inds(0,1))

   ; Difference between 'lat_max'/'lon_max' and values read in using NCL function
     dlat                   = abs(lat_max - centre(cn,it,0))
     dlon                   = abs(lon_max - centre(cn,it,1))

   ; If our calculations of latitude and longitude deviate too much from the values
   ; we read in using the NCL function ('st_rm'), then calculate a new centre
   ; using the minimum pressure. Constrain the new calculation using the values
   ; from 'st_rm', to prevent spurious centres being identified and analysed.

     if (dlat .gt. 0.75 .or. dlon .gt. 0.75) then

      lt1 = centre(cn,it,0) - 0.5
      lt2 = centre(cn,it,0) + 0.5
      ln1 = centre(cn,it,1) - 0.5
      ln2 = centre(cn,it,1) + 0.5

      prs_plane0 = prs_plane({lt1:lt2},{ln1:ln2})       ; Smaller grid (1.0 x 1.0 deg)
      prs_min0   = min(prs_plane0)                      ; Find minimum on smaller grid
      dims_p0    = dimsizes(prs_plane0)                 ; Size of smaller grid
      prs_1d0    = ndtooned(prs_plane0)                 ; Create 1-D array
      inds_p0    = ind_resolve(minind(prs_1d0),dims_p0) ; Find index of minimum

      lat_max    = lt1 + (dy * inds_p0(0,0))
      lon_max    = ln1 + (dx * inds_p0(0,1))
      print("New centre: "+lat_max+" degrees N, "+lon_max+" degrees E")

     end if

   else if (calc .eq. 1) then

      lat_max    = centre(cn,it,0)
      lon_max    = centre(cn,it,1)

   end if
   end if

   ;======================================================================
   ; Create lat/lon arrays for later switch to cylindrical coordinates
   ;======================================================================

   ;===============
   ; 'pc' stream
   ;===============

   ; Indices of lat/lon grid for 'u' and 'v' in 'pc' stream
     latC_plane  = ind(t0.le.latC.and.latC.le.t1)
     lonC_plane  = ind(n0.le.lonC.and.lonC.le.n1)

   ; Array sizes
     latC_size   = dimsizes(latC_plane)
     lonC_size   = dimsizes(lonC_plane)

   ; Array of latitude/longitude points corresponding to the indices above
     latC_sub    = latC(latC_plane(0):latC_plane(latC_size-1))
     lonC_sub    = lonC(lonC_plane(0):lonC_plane(lonC_size-1))

   ;===============
   ; 'pd' stream
   ;===============

   ; Indices of lat/lon grid for 'str_xz' in 'pd' stream
     latD_plane  = ind(t0.le.latD.and.latD.le.t1)
     lonD_plane  = ind(n0.le.lonD.and.lonD.le.n1)

   ; Array sizes
     latD_size   = dimsizes(latD_plane)
     lonD_size   = dimsizes(lonD_plane)

   ; Array of latitude/longitude points corresponding to the indices above
     latD_sub    = latD(latD_plane(0):latD_plane(latD_size-1))
     lonD_sub    = lonD(lonD_plane(0):lonD_plane(lonD_size-1))

   ; Indices of lat/lon grid for 'str_yz' in 'pd' stream
     lat1D_plane = ind(t0.le.lat1D.and.lat1D.le.t1)
     lon1D_plane = ind(n0.le.lon1D.and.lon1D.le.n1)

   ; Array sizes
     lat1D_size  = dimsizes(lat1D_plane)
     lon1D_size  = dimsizes(lon1D_plane)

   ; Array of latitude/longitude points corresponding to the indices above
     lat1D_sub   = lat1D(lat1D_plane(0):lat1D_plane(lat1D_size-1))
     lon1D_sub   = lon1D(lon1D_plane(0):lon1D_plane(lon1D_size-1))

     delete([/vort1d,inds,vort1d_sm,vort_smth,inds_sm/])

;===========================================
; Switch to cylindrical coordinate system
;===========================================

   ; Create new arrays describing cylindrical coordinate system (2D for now)
     xpos = new((/dimsizes(radii),dimsizes(thetas)/),float)
     ypos = new((/dimsizes(radii),dimsizes(thetas)/),float)

   ; New arrays containing variables we want to plot
     u_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Zonal wind
     v_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Meridional wind
     v_rad    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Radial wind
     v_tan    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Tangential wind
     vort_int = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Relative vorticity
     avo_int  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Absolute vorticity
     w_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vertical velocity

;=========================
; Add variable metadata
;=========================

     copy_VarAtts(u_plane,u_int)
     u_int@description = "Zonal velocity"
     u_int@units       = "m s~S~-1~N~"

     copy_VarAtts(v_plane,v_int)
     v_int@description = "Meridional velocity"
     v_int@units       = "m s~S~-1~N~"

     copy_VarAtts(w_plane,w_int)
     w_int@description = "Vertical velocity"
     w_int@units       = "m s~S~-1~N~"

     copy_VarAtts(vort_plane,vort_int)
     vort_int@description = "Relative vorticity"
     vort_int@units    = "s~S~-1~N~"

     copy_VarAtts(avo_plane,avo_int)
     avo_int@description = "Relative vorticity"
     avo_int@units    = "s~S~-1~N~"

     v_rad@description = "Radial velocity"
     v_rad@units       = "m s~S~-1~N~"

     v_tan@description = "Tangential velocity"
     v_tan@units       = "m s~S~-1~N~"


   ; Create cylindrincal coordinate array from original latitude/longitude grid
     do irad = 0, dimsizes(radii)-1
       r = radii(irad)
       do iang = 0, dimsizes(thetas)-1
         theta = thetas_rad(iang)
           xpos(irad,iang) = lon_max + r*cos(theta)
           ypos(irad,iang) = lat_max + r*sin(theta)
       end do
     end do

   ; Fill 'xpos_all' and 'ypos_all' with latitude & longitude values
     xpos_all(ilev,:,:,it) = xpos(:,:)
     ypos_all(ilev,:,:,it) = ypos(:,:)

   ; Fill new arrays with interpolated values on the cylindrical grid
   ; Which latitude/longitude grids are each of these variables defined on?
   ; u_plane     -- latC,  lonC
   ; v_plane     -- latC,  lonC
   ; vort_plane  -- latC,  lonC
   ; avo_plane   -- latC,  lonC
   ; w_plane     -- latD,  lon1D

     do iang = 0, dimsizes(thetas_rad)-1
      u_int(:,iang)    = linint2_points(lonC_sub,latC_sub,u_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      v_int(:,iang)    = linint2_points(lonC_sub,latC_sub,v_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      vort_int(:,iang) = linint2_points(lonC_sub,latC_sub,vort_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      avo_int(:,iang)  = linint2_points(lonC_sub,latC_sub,avo_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      w_int(:,iang)    = linint2_points(lon1D_sub,latD_sub,w_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
     end do

   ; Calculate radial/tangential windspeed and friction at each point (using unit vctrs)
      do iang = 0, dimsizes(thetas)-1
       theta = thetas_rad(iang)

      ; Tangential and radial wind
        v_rad(:,iang) = u_int(:,iang)*cos(theta) + v_int(:,iang)*sin(theta)
        v_tan(:,iang) = -u_int(:,iang)*sin(theta) + v_int(:,iang)*cos(theta)

      end do

   ; Add additional metadata to all arrays 
     xpos!0   = "rad"
     xpos!1   = "azi"
     xpos&rad = radii
     xpos&azi = thetas

     ypos!0   = "rad"
     ypos!1   = "azi"
     ypos&rad = radii
     ypos&azi = thetas

     u_int!0   = "rad"
     u_int!1   = "azi"
     u_int&rad = radii
     u_int&azi = thetas

     v_int!0   = "rad"
     v_int!1   = "azi"
     v_int&rad = radii
     v_int&azi = thetas

     v_tan!0 = "rad"
     v_tan!1 = "azi"
     v_tan&rad = radii
     v_tan&azi = thetas

     v_rad!0 = "rad"
     v_rad!1 = "azi"
     v_rad&rad = radii
     v_rad&azi = thetas

     vort_int!0   = "rad"
     vort_int!1   = "azi"
     vort_int&rad = radii
     vort_int&azi = thetas

     avo_int!0    = "rad"
     avo_int!1    = "azi"
     avo_int&rad  = radii
     avo_int&azi  = thetas

     w_int!0      = "rad"
     w_int!1      = "azi"
     w_int&rad    = radii
     w_int&azi    = thetas

    ;===================================================
    ; Calculate mean (azimuthally-averaged) variables
    ;===================================================

      do irad = 0, dimsizes(radii)-1

        mean_u0(ilev,ct,irad)	    = avg(u_int(irad,:))      ; [ {u} ]
	mean_v0(ilev,ct,irad)	    = avg(v_int(irad,:))      ; [ {v} ]
        mean_vtan0(ilev,ct,irad)    = avg(v_tan(irad,:))      ; [ {vtan} ]
        mean_vrad0(ilev,ct,irad)    = avg(v_rad(irad,:))      ; [ {vrad} ]
        mean_vvel0(ilev,ct,irad)    = avg(w_int(irad,:))      ; [ {w} ]
        mean_vort0(ilev,ct,irad)    = avg(vort_int(irad,:))   ; [ {zeta} ]
        mean_avo0(ilev,ct,irad)     = avg(avo_int(irad,:))    ; [ {zeta + f} ]

      end do

    ; Tidy up
      delete([/u_int,v_int,w_int,v_tan,v_rad,vort_int,avo_int/])
      delete([/u_plane,v_plane,w_plane,vort_plane,avo_plane,xpos,ypos/])
      delete([/latC_plane,lonC_plane,latD_plane,lonD_plane,lon1D_plane,lat1D_plane/])
      delete([/latC_sub,lonC_sub,latD_sub,lonD_sub,lat1D_sub,lon1D_sub/])

     end do                         ; End pressure level loop (do ilev = 0, num_l-1)

    ct  = ct + 1                    ; Counter variable (time)

  ; Tidy up before next iteration (all variables from 'pc' and 'pd' streams)
    delete([/lonC,latC,lon1C,lat1C,u,v,vort,avo,f0,f1/])
    delete([/lonD,latD,lon1D,lat1D,t,p,w,z/])

   end do     ; End time loop (do it = 0, times-1)

;==========================================
; Add metadata to arrays before plotting 
;==========================================

    ; Redefine model level array for vertical derivative calculations
      hybD   = d->hybrid_ht(:)         ; 63 model ('theta') levels
      mlevs  = dimsizes(hybD)          ; Size of model level array

      hybP   = hybD / 1000             ; Create array of levels in 'km'
      hybP@units = "km"                ; Model level units for plotting

    ; Add metadata to dimension 0 (model levels) of the arrays before plotting
      mean_vort0&lev = hybP
      mean_avo0&lev  = hybP
      mean_vtan0&lev = hybP
      mean_vrad0&lev = hybP
      mean_vvel0&lev = hybP
      mean_u0&lev    = hybP
      mean_v0&lev    = hybP

      mean_vort0&rad = radii
      mean_avo0&rad  = radii
      mean_vtan0&rad = radii
      mean_vrad0&rad = radii
      mean_vvel0&rad = radii
      mean_u0&rad    = radii
      mean_v0&rad    = radii

      if (calc .eq. 0) then

       mean_vort0 = smth9_Wrap(mean_vort0, 0.5, 0.5, True)
       mean_avo0  = smth9_Wrap(mean_avo0, 0.5,  0.5, True)
       mean_vtan0 = smth9_Wrap(mean_vtan0, 0.5, 0.5, True)
       mean_vrad0 = smth9_Wrap(mean_vrad0, 0.5, 0.5, True)
       mean_vvel0 = smth9_Wrap(mean_vvel0, 0.5, 0.5, True)
       mean_u0	  = smth9_Wrap(mean_u0, 0.5, 0.5, True)
       mean_v0	  = smth9_Wrap(mean_v0, 0.5, 0.5, True)

      end if

;============================================
; Define output file path and plot details
;============================================

 do it = ts0, tf0

  print("Working on time: "+it)

  do ilev = ls0, lf0

    hy = sprintf("%0.0f",hybP(ilev))
    print("Working on model level " + ilev+ " ("+hy+" m AGL)")

  ; Output file location and type
    output = "$sam/nepartak/images/cyl_coords/basic_xy_"+hy+"_"+ens0+"_"+time_arr(it)
    wks = gsn_open_wks(opt,output)

  ; Set up panel plot (plot either 1 or 2 panels)

    if (full .eq. 1) then
     panel = new(2,graphic)
    elseif (full .eq. 0) then
     panel = new(1,graphic)
    end if

  ; Read in desired colourmap
    if (clr .eq. 1) then
     gsn_define_colormap(wks,"ncl_default")
    elseif (clr .eq. 2) then
     gsn_define_colormap(wks,"bdr_extra")
    elseif (clr .eq. 3) then
     gsn_define_colormap(wks,"amwg256")
    end if

;==========================
; Options for plotting 
;==========================

  ; Relative vorticity
    opts_vort				  = True
    opts_vort@sfXArray			  = xpos_all(ilev,:,:,it)
    opts_vort@sfYArray			  = ypos_all(ilev,:,:,it)
    opts_vort@cnFillOn			  = True
    opts_vort@cnLineLabelInterval         = 2.0
    opts_vort@cnLineLabelFontHeightF      = 0.012
    opts_vort@cnLineLabelBackgroundColor  = "transparent"
    opts_vort@cnLineLabelPlacementMode    = "constant"
    opts_vort@cnLinesOn                   = False
    opts_vort@cnInfoLabelOn               = False
    opts_vort@cnLevelSelectionMode        = "ExplicitLevels"    

    if (clr .eq. 1) then
     opts_vort@cnFillColors               = (/66,78,90,102,120,0,0,\
                                              136,152,168,180,188,196,\
                                              204,216,228,240,252/)
     opts_vort@cnLevels                   = (/-5., -2., -1., -0.5, -0.2, \
                                               0, 0.2, 0.5, 1., 2., 3., \
                                               5., 10., 20., 50., 100., 150./)

    elseif (clr .eq. 2) then
     opts_vort@cnFillColors               = (/8,9,10,11,12,0,0,\
                                              15,16,17,18,19,20,21,22,23,24,25/)
     opts_vort@cnLevels                   = (/-5., -2., -1., -0.5, -0.2, \
                                               0, 0.2, 0.5, 1., 2., 3., \
                                               5., 10., 20., 50., 100., 150./)
    elseif (clr .eq. 3) then
     opts_vort@cnFillColors               = (/47,56,68,80,92,0,0,\
                                              150,158,168,180,188,196,\
                                              204,216,228,240,252/)
     opts_vort@cnLevels                   = (/-5., -2., -1., -0.5, -0.2, \
                                               0, 0.2, 0.5, 1., 2., 3., \
                                               5., 10., 20., 50., 100., 150./)
    end if

    opts_vort@gsnPaperOrientation          = "landscape"
    opts_vort@tiMainString                 = ""
    opts_vort@tiMainFontHeightF            = 0.0125
    opts_vort@gsnLeftString                = ""
    opts_vort@gsnRightString               = ""
    opts_vort@gsnMaximize                  = True
    opts_vort@lbLabelBarOn                 = False
    opts_vort@gsnDraw                      = False
    opts_vort@gsnFrame                     = False


    opts_vort@pmLabelBarWidthF	    	   = 0.6
    opts_vort@pmLabelBarHeightF	    	   = 0.15
    opts_vort@lbLabelFontHeightF	   = 0.01
    opts_vort@lbPerimOn		           = False
    opts_vort@lbLabelFont		   = "Helvetica"
    opts_vort@lbTitleOn		    	   = True
    opts_vort@lbTitleString = hy+" km relative vorticity (10~S~-4~N~ s~S~-1~N~)"
    opts_vort@lbTitleFontHeightF	   = 0.015
    opts_vort@lbTitleDirection	    	   = "Across"

  ; Additional plotting resources 
    opts_vort@mpDataBaseVersion	           = "Ncarg4_1" ; More recent database 
    opts_vort@mpDataSetName		   = "Earth..4" ; High resolution 
    opts_vort@mpOutlineBoundarySets 	   = "National" ; National borders
    opts_vort@mpGeophysicalLineColor	   = "black"    ; Colour borders black
    opts_vort@mpGeophysicalLineThicknessF  = 1.0        ; Border line thickness
    opts_vort@mpGridAndLimbOn	    	   = False      ; Turn on lat/lon lines
    opts_vort@pmTickMarkDisplayMode	   = "Always"   ; Turn on map tickmarks
    opts_vort@tmXBMajorLengthF	    	   = 0.005      ; Change tickmark length
    opts_vort@tmXTOn			   = "False"    ; No tickmarks on top x-axis
    opts_vort@tmYROn			   = "False"    ; No tickmarks on right y-axis
    opts_vort@gsnMaximize	 	   = True       ; Maximise plot size
    opts_vort@gsnAddCyclic                 = False

    radius				   = radii(rad_size) ; Size of plot (lat/lon)

    opts_vort@mpLimitMode		   = "Corners"
    opts_vort@mpLeftCornerLatF             = centre(cn,it,0)-radius
    opts_vort@mpLeftCornerLonF             = centre(cn,it,1)-radius
    opts_vort@mpRightCornerLatF            = centre(cn,it,0)+radius
    opts_vort@mpRightCornerLonF            = centre(cn,it,1)+radius
    opts_vort@gsnDraw		           = False      ; Do not draw the plot
    opts_vort@gsnFrame		           = False      ; Do not advance the frame

  ; Horizontal wind vector resources
    opts_vec				   = True

 
  ; Vertical velocity resources
    opts_vvel				   = True 

;=============
; Plot data 
;=============

  ; Plot relative vorticity 
    plot_vort  = gsn_csm_contour_map(wks,mean_vort0(ilev,it,:),opts_vort)

  ; Plot horizontal wind vectors
    plot_vec   = gsn_csm_vector(wks,u_plane,v_plane,opts_vec)

  ; Plot vertical velocity
    plot_vvel  = gsn_csm_contour(wks,w_plane,opts_vvel)

    overlay(plot_vort,plot_vvel)
    overlay(plot_vort,plot_vec)				       

;=========================================
; Add radial circles to plot (optional)
;=========================================

    if (rad .eq. 1) then 

      dtr = 0.017453292519943

    ; Get customizations for circles.
      radius  = radii(rad_size)
      xcenter = lon_max
      ycenter = lat_max
      outthck = 2
      inthck = 1
      spacing = 0.5
      dpattrn = 2
      outcolr = 1
      incolr = 1

    ;=============================================
    ; Construct and attach outer circle to plot
    ;=============================================

    ; Calculate arrays for outer circle (similar to what I have done)
      degrees = ispan(0,360,5)
      xcos    = cos(dtr * degrees)
      xsin    = sin(dtr * degrees)

    ; Array of points (x,y) representing the outer circle (radius = 25)
      xc      = xcenter + radius  * xcos
      yc      = ycenter + radius  * xsin

    ; Resources for outer circle
      lnres                   = True
      lnres@gsLineThicknessF  = outthck
      lnres@gsLineColor       = outcolr
      lnres@gsLineDashPattern = 0

    ; Attach outer circle to plot
      plot_outer     = gsn_add_polyline(wks,plot_vort,xc,yc,lnres)

      if(spacing.gt.radius) then
        print("add_radial_circles: spacing is > radius, can't draw inner circles.")
        return
      end if

    ;==================================================
    ; Now construct and attach inner circles to plot
    ;==================================================

    ; Draw inner circles if desired ('spacing' controls the number of inner circles)
      if(spacing.gt.0) then
        count = 0
              size = toint((radius-spacing)/spacing)
              plot_inner = new(size,graphic)
        do r = spacing,radius-spacing,spacing

        ; Calculate arrays for inner circle.
          xc = xcenter + (r * xcos)
          yc = ycenter + (r * xsin)

        ; Resources for inner circle
          delete(lnres@gsLineThicknessF)   ; Delete in case setting to
          delete(lnres@gsLineColor)        ; new type.
          delete(lnres@gsLineDashPattern)
          lnres@gsLineColor       = incolr
          lnres@gsLineThicknessF  = inthck
          lnres@gsLineDashPattern = dpattrn

	  ; Attach inner circle to plot
          plot_inner(count) = gsn_add_polyline(wks,plot_vort,xc,yc,lnres)
          count = count + 1
        end do
      end if

    ; Tidy up
      delete([/lnres, degrees, xcos, xsin, xc, yc/])

;=====================
; Add radial lines
;=====================

    ; Get customizations for radial lines (radius, xcenter, ycenter are already defined)
      angle   = 45
      lcolor  = 1
      dpattrn = 2
      thcknss = 1.

    ; Error checking
      if(angle.le.0.or.angle.ge.360) then
        print("add_radial_lines: angle must be between 0 and 360.")
        print("                  Can't draw radial lines.")
        return
      end if

    ; Resources for radial lines
      lnres                   = True
      lnres@gsLineDashPattern = dpattrn
      lnres@gsLineColor       = lcolor
      lnres@gsLineThicknessF  = thcknss

    ; Calculate arrays for lines
      xlines  = 360./angle
      nlines  = toint(ceil(xlines)) ; Smallest integer greater than/equal to xlines
      if(xlines.ne.nlines) then
        degrees = fspan(0,360,nlines)
      else
        degrees = fspan(0,360-angle,nlines) ; Contains angles (0, 45, 90, etc)
      end if

      nlines  = dimsizes(degrees)
      plot_lines = new(nlines,graphic)
      do i = 0,nlines-1
        xl = (/xcenter, xcenter + (radius  * cos(dtr * degrees(i)))/)
        yl = (/ycenter, ycenter + (radius  * sin(dtr * degrees(i)))/)

      ; Attach line to plot
        plot_lines(i) = gsn_add_polyline(wks,plot_vort,xl,yl,lnres)

      end do

    ; Tidy up
      delete(degrees)

;=============================================
; Add marker illustrating windspeed maximum
;=============================================

    ; Add circle indicating the radius of maximum windspeed
      rmw                     = inds_s(0,0)     ; Radial distance from TC centre
      thta                    = inds_s(0,1)     ; Azimuthal position (angle)
      xcentre                 = lon_max         ; TC centre (lon)
      ycentre                 = lat_max         ; TC centre (lat)

    ; Calculate circular array of points (every 5 deg)
      degrees                 = fspan(0.0,360.0,73)
      xcos                    = cos(pii/180 * degrees) ; cosine of angles in array (rad)
      xsin                    = sin(pii/180 * degrees) ; sine of angles in array (rad)

    ; First add marker corresponding to exact location of maximum
      marker = new(1,graphic)
      mres   = True
      mres@gsMarkerIndex = 5        ; Cross
      mres@gsMarkerThicknessF = 3.0 ; Marker thickness
      mres@gsMarkerSizeF = 16.0     ; Marker size
      mres@gsMarkerColor = "black"  ; Marker colour

    ; Calculate lat/lon of windspeed maximum on cylindrical grid
      cpos_x = xcentre + (radii(rmw) * xcos(thta))
      cpos_y = ycentre + (radii(rmw) * xsin(thta))

;      marker = gsn_add_polymarker(wks,plot_vort,cpos_x,cpos_y,mres)
;      delete(mres)

;================================================================
; Finally, draw circle to highlight the radius of maximum wind
;================================================================

    ; Additional options for circle representing maximum wind
      thck                    = 5
      color                   = 1

    ; Array of points (x,y) representing the circle
      xc                      = xcentre + radii(rmw) * xcos ; origin + [r * cos(theta)]
      yc                      = ycentre + radii(rmw) * xsin ; origin + [r * sin(theta)]

    ; Resources for circle indicating the radius of maximum wind
      lnres                   = True
      lnres@gsLineThicknessF  = thck
      lnres@gsLineColor       = color  ; Option only recognises numbers (color = 0)
      lnres@gsLineDashPattern = 0

    ; Attach circle to plot
      plot_rmw                = gsn_add_polyline(wks,plot_vort,xc,yc,lnres)
      delete([/lnres, xc, yc, degrees, xcos, xsin, rmw, xcentre, ycentre, thck, color/])

    end if 

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

      draw(plot_vort)
      frame(wks)

    ; Tidy up before next iteration
      delete([/vort_plane,u_plane,v_plane,t_plane,geo_plane,lat_plane,lon_plane/])
      delete([/u_int,v_int,v_tan,v_rad,speed/])

     end do	   ; End model level loop (do ilev = 0, plevs-1)

    end do     ; End time loop (do it = 0, 1)

  ; Tidy up before next iteration
    delete([/date_str, time_str, out_str/])
    delete([/time,times,u,v,z,t,w,vort,pres,plevs,div/])

end 