; Script to calculate the radial momentum budget as in Montgomery et al. (2018)

; Run using:

; ncl 'dat="02T12"' 'opt="x11"' 'ens0="em10"' 'ens1="em08"' 'ens2="em11"' 
; dist=3.0 ts0=96 tf0=96 calc=1 full=1 tc_basic_xz_new.ncl

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'opt'   = output file format ("pdf" or "x11")
; 'ens0'  = ensemble member #1 (em00 -- em11)
; 'ens1'  = ensemble member #2 (em00 -- em11)
; 'ens2'  = ensemble member #3 (em00 -- em11)
; 'dist'  = size of box following storm (degrees)
; 'ts0'   = starting time index (0–119)
; 'tf0'   = ending time index (0–119)
; 'calc'  = centre calculated on each level (0) or equal to 950 hPa centre (1)
; 'full'  = all 10 panels plotted (1) or only a subset of 4 selected panels (0)

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/nc_levs.ncl"

begin

;==============================
; Prelash (reading in files)
;==============================

; Create string array containing chosen ensemble members 
  ens_arr     = (/ens0, ens1, ens2/)
  en_size     = dimsizes(ens_arr)

  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens_arr(0)
  fili_c1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pc*.nc")
  fili_c      = diri+"/"+fili_c1

; Get dimension information from one of the files 
  numINFO_0   = nc_times(fili_c)        ; Retrieve number of times 
  numTIMES    = numINFO_0[0]            ; Number of times in 'pc' stream
  time_0      = numINFO_0[1]            ; Array of times from 'pc' stream file

  numINFO_1   = nc_levs(fili_c)		; Retrieve number of vertical levels
  num_l       = numINFO_1[1]            ; Number of vertical levels in 'pc' stream file
  delete([/fili_prefix,fili_c1,fili_c/])

;============================================================================
; Define arrays to hold output variables based upon dimension information 
;============================================================================

; Define cylindrical coordinate array
  thetas     = new(73,float)	     ; Array of azimuth angles (0, 5, 10, etc)
  pii        = 3.14159265	     
  radii      = fspan(0.0,1.5,31)     ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,73)   ; Fill azimuth angle array
  thetas_rad = (pii/180)*thetas      ; Degrees to radians
  rad_size   = dimsizes(radii)-1     ; Index for later calculation

; Define constants
  rd    = 287.0   ; Gas constant for dry air (J/kg/K^2)
  cp0   = 1004.0  ; Specific heat of dry air at constant pressure (J/kg/K)
  kp    = 0.286   ; For Exner function calculations (0.286)
  pref  = 1000.0  ; Reference pressure p0 (1000 hPa)

; Array to hold all 'xpos' and 'ypos' position arrays
  xpos_all  = new((/num_l,en_size,numTIMES,dimsizes(radii),dimsizes(thetas)/),float)
  ypos_all  = new((/num_l,en_size,numTIMES,dimsizes(radii),dimsizes(thetas)/),float)

;================================================
; Define 3D arrays for mean quantities
; Dimension 0 = vertical levels ("lev")
; Dimension 1 = times           ("time")
; Dimension 2 = radius          ("rad")
; This order is required for 'gsn_csm_contour'
;================================================

; Mean tangential wind {v}
  mean_vtan0       = new((/num_l,en_size,numTIMES,dimsizes(radii)/),float)
  mean_vtan0!0     = "lev"
  mean_vtan0!1	   = "ens"
  mean_vtan0!2     = "time"
  mean_vtan0!3     = "rad"
  mean_vtan0@description = "Azimuthally averaged tangential wind"
  mean_vtan0@units = "m s~S~-1~N~"

; Mean radial wind {u}
  mean_vrad0       = new((/num_l,en_size,numTIMES,dimsizes(radii)/),float)
  mean_vrad0!0     = "lev"
  mean_vrad0!1     = "ens"
  mean_vrad0!2     = "time"
  mean_vrad0!3     = "rad"
  mean_vrad0@description = "Azimuthally averaged radial wind"
  mean_vrad0@units = "m s~S~-1~N~"

; Mean vertical velocity {w}
  mean_vvel0       = new((/num_l,en_size,numTIMES,dimsizes(radii)/),float)
  mean_vvel0!0     = "lev"
  mean_vvel0!1     = "ens"
  mean_vvel0!2     = "time"
  mean_vvel0!3     = "rad"
  mean_vvel0@description = "Azimuthally averaged vertical velocity"
  mean_vvel0@units = "m s~S~-1~N~"

; Mean relative vorticity {zeta}
  mean_vort0       = new((/num_l,en_size,numTIMES,dimsizes(radii)/),float)
  mean_vort0!0     = "lev"
  mean_vort0!1     = "ens"
  mean_vort0!2     = "time"
  mean_vort0!3     = "rad"
  mean_vort0@description = "Azimuthally averaged relative vorticity"
  mean_vort0@units = "s~S~-1~N~"

; Mean absolute vorticity {zeta + f0}
  mean_avo0        = new((/num_l,en_size,numTIMES,dimsizes(radii)/),float)
  mean_avo0!0      = "lev"
  mean_avo0!1      = "ens"
  mean_avo0!2      = "time"
  mean_avo0!3      = "rad"
  mean_avo0@description = "Azimuthally averaged absolute vorticity"
  mean_avo0@units  = "s~S~-1~N~"

; Mean potential temperature {theta}
  mean_thta0       = new((/num_l,en_size,numTIMES,dimsizes(radii)/),float)
  mean_thta0!0     = "lev"
  mean_thta0!1     = "ens"
  mean_thta0!2     = "time"
  mean_thta0!3     = "rad"
  mean_thta0@description = "Azimuthally averaged potential temperature"
  mean_thta0@units = "K"

;=====================================
; Start loop over ensemble members 
;=====================================

do en = 0, en_size-1

  print_clock("Working on ensemble member: "+ens_arr(en) )

; Find 'pc' and 'pd' files (model height level data)
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens_arr(en)
  fili_c1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pc*.nc")
  fili_c      = diri+"/"+fili_c1
  fili_d1     =	systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pd*.nc")
  fili_d      = diri+"/"+fili_d1

; Calculate number of times in each file using built-in function ('nc_times')
  numINFO_c   = nc_times(fili_c)
  numINFO_d   = nc_times(fili_d)
  numINPUT_c  = numINFO_c[0]		; Number of times in 'pc' stream
  time_c      = numINFO_c[1]		; Array of times from 'pc' stream file

; Calculate storm motion using built-in function ('st_rm')
; Argument 1  = ensemble simulation
; Argument 2  = path to input files
; Argument 3  = size of box following storm (degrees)

; Call external function 'st_rm'
  storm_rel = st_rm(dat,ens_arr(en),"$ar/text/",3.0)
  u_cyc     = storm_rel[0]		   ; Zonal wind 
  v_cyc     = storm_rel[1]		   ; Meridional wind
  vel_cyc   = storm_rel[2]		   ; Vector wind 
  numTIMES  = storm_rel[3]                 ; Number of times in file
  lat_arr   = storm_rel[4]                 ; Storm-following domain subset (lat)
  lon_arr   = storm_rel[5]                 ; Storm-following domain subset (lon)
  centre    = storm_rel[6]		   ; Storm centre at each hour (lat/lon)

;==========================================================================
; Create date/time string arrays for all times in file (use 'pb' stream)
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

;=========================
; Start multiple loops
;=========================

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
    llbox	= toint(dist*50)	    ; Size of lat/lon arrays (see line below)
    dsize	= (/llbox,llbox/)           ; Size of lat/lon arrays (depends on 'dist')
    d0		= dsize(0)
    d1		= dsize(1)

    setvalues NhlGetWorkspaceObjectId
      "wsMaximumSize" : 1000000000
    end setvalues

;===============================
; Start loop over input files 
;===============================

    c  = addfile(fili_c,"r")         	 ; Read in 'pc' stream [u,v,vort]
    d  = addfile(fili_d,"r")         	 ; Read in 'pd' stream [w,p,t,theta,kmh,tau]

;==================================
; Get the variables we will need
;==================================

   do it = ts0, tf0		; Loop over times in file

    print("Working on time: "+time_str(it))
    itt = it + 1
    title_arr(it) = "Valid at "+time_str(it)+" (T+"+itt+")"

;=====================================================
; Read in variables from 'pc' stream (model levels)
;=====================================================

    t0    = lat_arr(it,0)
    t1    = lat_arr(it,1)
    n0    = lon_arr(it,0)
    n1    = lon_arr(it,1)
    print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

  ; Read in basic variables from 'pc' stream
    lon1C  = c->longitude_1({n0:n1})              ; longitude ['d1' grid points]
    lat1C  = c->latitude_1({t0:t1})               ; latitude  ['d0' grid points]
    lonC   = c->longitude({n0:n1})                ; longitude ['d1' grid points]
    latC   = c->latitude({t0:t1})                 ; latitude  ['d0' grid points]
    hybC   = c->hybrid_ht(:) 		          ; 63 model ('theta') levels

    u_varname = "x-wind"
    v_varname = "y-wind"

  ; Read in horizontal wind components 
  ; Both variables below -- [hybC | 63] * [latC | 150] * [lonC | 150]
    u  = c->$u_varname$(it,:,{t0:t1},{n0:n1}) ; Zonal wind	
    v  = c->$v_varname$(it,:,{t0:t1},{n0:n1}) ; Meridional wind 

    print(hybC)
    exit()

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
    hybD   = d->hybrid_ht(:)			    ; 63 model ('rho') levels

  ; Read in temperature, pressure, vertical velocity and geopotential height
  ; All variables below -- [hybD | 63] * [latD | 150] * [lon1D | 151]
    t      = d->temp(it-2,:,{t0:t1},{n0:n1})        ; Temperature (K)
    th	   = d->theta(it-2,:,{t0:t1},{n0:n1})	    ; Potential temperature (K)
    p      = d->p(it-2,:,{t0:t1},{n0:n1})           ; Pressure (Pa)
    w      = d->dz_dt(it-2,:,{t0:t1},{n0:n1})	    ; Vertical velocity (m s-1)
    z      = d->ht(it-2,:,{t0:t1},{n0:n1})    	    ; Geopotential height (m)

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

    do ilev = 0, num_l-1 ; Loop over all pressure levels

     hy = sprintf("%0.0f",hybD(ilev))
     print("Working on model level " + ilev+ " ("+hy+" m AGL)")

   ; Read in other variables on model levels
     u_plane 	            = u(ilev,:,:)	   ; Zonal wind 
     v_plane		    = v(ilev,:,:)	   ; Meridional wind 
     w_plane		    = w(ilev,:,:)	   ; Vertical velocity
     prs_plane		    = p(ilev,:,:)	   ; Pressure
     vort_plane	   	    = vort(ilev,:,:)	   ; Relative vorticity
     avo_plane		    = avo(ilev,:,:)	   ; Absolute vorticity
     thta_plane		    = th(ilev,:,:)	   ; Potential temperature

   ; Calculate maximum relative vorticity on model level
     vort_max      	    = max(vort_plane)
     vort_smth     	    = smth9_Wrap(vort_plane, 0.5, 0.5, True)
     vort_max_smth 	    = max(vort_smth)

   ; Find index of relative vorticity maximum
     dims          	    = dimsizes(vort_plane)
     vort1d        	    = ndtooned(vort_plane)
     inds          	    = ind_resolve(maxind(vort1d),dims)
     vort1d_sm     	    = ndtooned(vort_smth)
     inds_sm       	    = ind_resolve(maxind(vort1d_sm),dims)

   ;==================================================================
   ; Either calculate the centre, or read in pre-calculated values
   ;==================================================================

   if (calc .eq. 0) then

   ; Find latitude and longitude of storm at each time and on each model level
     lat_max 	            = t0 + (dy * inds(0,0))
     lon_max 	   	    = n0 + (dx * inds(0,1))

   ; Difference between 'lat_max'/'lon_max' and values read in using NCL function
     dlat 	   	    = abs(lat_max - centre(it,0))
     dlon	   	    = abs(lon_max - centre(it,1))

   ; If our calculations of latitude and longitude deviate too much from the values
   ; we read in using the NCL function ('st_rm'), then calculate a new centre
   ; using the minimum pressure. Constrain the new calculation using the values
   ; from 'st_rm', to prevent spurious centres being identified and analysed. 

     if (dlat .gt. 0.25 .or. dlon .gt. 0.25) then 

      lt1 = centre(it,0) - 0.5
      lt2 = centre(it,0) + 0.5
      ln1 = centre(it,1) - 0.5
      ln2 = centre(it,1) + 0.5

      prs_plane0 = prs_plane({lt1:lt2},{ln1:ln2})       ; Smaller grid (1.0 x 1.0 deg)
      prs_min0   = min(prs_plane0)                      ; Find minimum on smaller grid
      dims_p0    = dimsizes(prs_plane0)                 ; Size of smaller grid
      prs_1d0    = ndtooned(prs_plane0)                 ; Create 1-D array
      inds_p0    = ind_resolve(minind(prs_1d0),dims_p0) ; Find index of minimum

      lat_max 	 = lt1 + (dy * inds_p0(0,0))
      lon_max 	 = ln1 + (dx * inds_p0(0,1))
      print("New centre: "+lat_max+" degrees N, "+lon_max+" degrees E")

     end if 

   else if (calc .eq. 1) then

      lat_max    = centre(it,0)
      lon_max    = centre(it,1)

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
     
   ; Array of latitude/longitude points	corresponding to the indices above
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
     xpos     = new((/dimsizes(radii),dimsizes(thetas)/),float)
     ypos     = new((/dimsizes(radii),dimsizes(thetas)/),float)

   ; New arrays containing regular variables we want to plot
     u_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Zonal velocity
     v_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Meridional velocity
     w_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vertical velocity
     vort_int = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Relative vorticity
     avo_int  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Absolute vorticity
     thta_int = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Potential temperature
     v_rad    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Radial velocity
     v_tan    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Tangential velocity

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
 
     copy_VarAtts(thta_plane,thta_int)
     thta_int@description = "Potential temperature"
     thta_int@units    = "K"

   ; Also add metadata for radial and tangential wind components
     v_rad@description    = "Radial velocity"
     v_rad@units          = "m s~S~-1~N~"

     v_tan@description    = "Tangential velocity"
     v_tan@units          = "m s~S~-1~N~"

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
     xpos_all(ilev,en,it,:,:) = xpos(:,:)
     ypos_all(ilev,en,it,:,:) = ypos(:,:)

   ; Fill new arrays with interpolated values on the cylindrical grid 
   ; Which latitude/longitude grids are each of these variables defined on? 
   ; u_plane     -- latC,  lonC
   ; v_plane     -- latC,  lonC
   ; vort_plane	 -- latC,  lonC
   ; avo_plane	 -- latC,  lonC

   ; w_plane     -- latD,  lon1D
   ; thta_plane	 -- latD,  lon1D

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
      thta_int(:,iang) = linint2_points(lon1D_sub,latD_sub,thta_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
     end do

   ; Calculate radial/tangential windspeed and friction at each point (using unit vctrs)
      do iang = 0, dimsizes(thetas)-1
       theta = thetas_rad(iang)

      ; Tangential and radial wind 
        v_rad(:,iang) = u_int(:,iang)*cos(theta) + v_int(:,iang)*sin(theta)
        v_tan(:,iang) = -u_int(:,iang)*sin(theta) + v_int(:,iang)*cos(theta)

      end do

   ; Add metadata to new cylindrical coordinate arrays
     xpos!0       = "rad"
     xpos!1       = "azi"
     xpos&rad     = radii
     xpos&azi     = thetas

     ypos!0       = "rad"
     ypos!1       = "azi"
     ypos&rad     = radii
     ypos&azi     = thetas

     u_int!0      = "rad"
     u_int!1      = "azi"
     u_int&rad    = radii
     u_int&azi    = thetas

     v_int!0      = "rad"
     v_int!1      = "azi"
     v_int&rad    = radii
     v_int&azi    = thetas

     vort_int!0   = "rad"
     vort_int!1   = "azi"
     vort_int&rad = radii
     vort_int&azi = thetas

     avo_int!0    = "rad"
     avo_int!1    = "azi"
     avo_int&rad  = radii
     avo_int&azi  = thetas

     thta_int!0   = "rad"
     thta_int!1   = "azi"
     thta_int&rad = radii
     thta_int&azi = thetas

     w_int!0      = "rad"
     w_int!1      = "azi"
     w_int&rad    = radii
     w_int&azi    = thetas

     v_tan!0      = "rad"
     v_tan!1      = "azi"
     v_tan&rad    = radii
     v_tan&azi    = thetas

     v_rad!0      = "rad"
     v_rad!1      = "azi"
     v_rad&rad    = radii
     v_rad&azi    = thetas

    ;===================================================
    ; Calculate mean (azimuthally-averaged) variables 
    ;===================================================

      do irad = 0, dimsizes(radii)-1

        mean_vtan0(ilev,en,ct,irad)    = avg(v_tan(irad,:))      ; [ {v} ]
	mean_vrad0(ilev,en,ct,irad)    = avg(v_rad(irad,:))      ; [ {u} ]
	mean_vvel0(ilev,en,ct,irad)    = avg(w_int(irad,:))      ; [ {w} ]
	mean_vort0(ilev,en,ct,irad)    = avg(vort_int(irad,:))   ; [ {zeta} ]
	mean_avo0(ilev,en,ct,irad)     = avg(avo_int(irad,:))    ; [ {zeta + f} ]
	mean_thta0(ilev,en,ct,irad)    = avg(thta_int(irad,:))   ; [ {theta} ]

      end do

    ; Tidy up
      delete([/u_int,v_int,w_int,v_tan,v_rad,vort_int,avo_int/])
      delete([/u_plane,v_plane,w_plane,vort_plane,prs_plane,avo_plane,xpos,ypos/])
      delete([/latC_plane,lonC_plane,latD_plane,lonD_plane,lon1D_plane,lat1D_plane/])
      delete([/latC_sub,lonC_sub,latD_sub,lonD_sub,lat1D_sub,lon1D_sub/])

     end do			    ; End pressure level loop (do ilev = 0, num_l-1)

    ct  = ct + 1		    ; Counter variable (time)

  ; Tidy up before next iteration (all variables from 'pc' and 'pd' streams)
    delete([/lonC,latC,lon1C,lat1C,u,v,vort,avo,f0,f1/])
    delete([/lonD,latD,lon1D,lat1D,th,t,p,w,z/])

   end do     ; End time loop (do it = 0, times-1)

  end do      ; End ensemble simulation loop (do en = 0, en_size-1)

;================================================
; Calculate vertical derivatives of mean terms 
;================================================

    ; Redefine model level array for vertical derivative calculations
      hybC   = c->hybrid_ht(:)         ; 63 model ('theta') levels
      mlevs  = dimsizes(hybC)	       ; Size of model level array

      hybP   = hybC / 1000	       ; Create array of levels in 'km'
      hybP@units = "km"		       ; Model level units for plotting  

    ; Add metadata to dimension 0 (model levels) of the arrays before plotting 
      mean_thta0&lev = hybP
      mean_vort0&lev = hybP
      mean_avo0&lev  = hybP
      mean_vtan0&lev = hybP
      mean_vrad0&lev = hybP
      mean_vvel0&lev = hybP

      mean_thta0&ens = ens_arr
      mean_vort0&ens = ens_arr
      mean_avo0&ens  = ens_arr
      mean_vtan0&ens = ens_arr
      mean_vrad0&ens = ens_arr
      mean_vvel0&ens = ens_arr

      mean_thta0&rad = radii
      mean_vort0&rad = radii
      mean_avo0&rad  = radii
      mean_vtan0&rad = radii
      mean_vrad0&rad = radii
      mean_vvel0&rad = radii

      if (calc .eq. 0) then 

       mean_thta0 = smth9_Wrap(mean_thta0, 0.5, 0.5, True)
       mean_vort0 = smth9_Wrap(mean_vort0, 0.5,	0.5, True)
       mean_avo0  = smth9_Wrap(mean_avo0, 0.5,	0.5, True)
       mean_vtan0 = smth9_Wrap(mean_vtan0, 0.5,	0.5, True)
       mean_vrad0 = smth9_Wrap(mean_vrad0, 0.5,	0.5, True)
       mean_vvel0 = smth9_Wrap(mean_vvel0, 0.5,	0.5, True)
            
      end if 

;============================================
; Define output file path and plot details 
;============================================

 do it = ts0, tf0

    print("Working on time: "+it)

  ; Output file location and type
  
    if (full .eq. 1) then 
     output = "$sam/nepartak/images/cyl_coords/vrad_xz_"+ens0+"_"+ens1+"_"+time_arr(it)
    else
     output = "$sam/nepartak/images/cyl_coords/vrad_xz_"+ens0+"_"+ens1+\
     	      "_"+ens2+"_"+time_arr(it)
    end if 
   
    wks = gsn_open_wks(opt,output)

  ; Set up panel plot

    if (full .eq. 2) then 
     panel = new(3,graphic)
    else if (full .eq. 1) then 
     panel = new(2,graphic)
    else if (full .eq. 0) then 
     panel = new(1,graphic)
    end if 
    end if 
    end if 

  ; Read in desired colourmap
    gsn_define_colormap(wks,"BlueDarkRed18")

;==========================
; Options for plotting 
;==========================

  ; Radial velocity
    opts_vrad                             = True
    opts_vrad@cnFillOn                    = True
    opts_vrad@cnLineLabelInterval         = 2.0
    opts_vrad@cnLineLabelFontHeightF      = 0.012
    opts_vrad@cnLineLabelBackgroundColor  = "transparent"
    opts_vrad@cnLineLabelPlacementMode    = "constant"
    opts_vrad@cnLinesOn                   = False ; Contour lines off
    opts_vrad@cnInfoLabelOn               = False
    opts_vrad@cnLevelSelectionMode        = "ExplicitLevels"
    opts_vrad@cnLevels                    = (/-10., -5., -3., -2., -1., -0.5, -0.2, \
                                              0, 0.2, 0.5, 1., 2., 3., 5., 10./)
    opts_vrad@cnFillColors                = (/2,4,5,6,7,8,9,0,0,12,13,14,15,16,17,19/)
    opts_vrad@gsnPaperOrientation         = "landscape"
    opts_vrad@tiMainString                = ""
    opts_vrad@tiMainFontHeightF           = 0.0125
    opts_vrad@gsnLeftString               = ""
    opts_vrad@gsnRightString              = ""
    opts_vrad@gsnMaximize                 = True
    opts_vrad@lbLabelBarOn                = False
    opts_vrad@gsnDraw                     = False
    opts_vrad@gsnFrame                    = False

  ; Additional plotting resources
    opts_vrad@tiYAxisString                = "Height (km)"
    opts_vrad@tiXAxisString                = "Radial distance from cyclone centre (km)"
    opts_vrad@tiXAxisFontHeightF           = 0.025   ; X-axis title font height 

    opts_vrad@tmXTOn                       = "False" ; Turn off top x-axis TM
    opts_vrad@tmYROn                       = "False" ; Turn off right y-axis TM
    opts_vrad@gsnMaximize                  = True    ; Maximise plot size
    opts_vrad@gsnAddCyclic                 = False
    opts_vrad@tmXBMode                     = "Explicit"  ; Set tick marks explicitly
    opts_vrad@tmXBValues                   = (/0.0,0.25,0.5,0.75,1.0,1.25,1.5/)	; Tm pos
    opts_vrad@tmXBLabels                   = (/"0","25","50","75","100","125","150"/)
    opts_vrad@tmXBLabelFontHeightF	   = 0.02    ; Label font height 

    opts_vrad@tmYLMode			   = "Explicit"
    opts_vrad@tmYLValues		   = fspan(2.0, 18.0, 9)
    opts_vrad@tmYLLabels		   = (/"2.0","4.0","6.0","8.0","10.0",\
    					       "12.0","14.0","16.0","18.0"/)
    opts_vrad@tmYLLabelFontHeightF	   = 0.02
    opts_vrad@gsnYAxisIrregular2Linear	   = True	; Linearise the y-axis 

    opts_vrad@trYMinF			   = 0.0	; y-axis minimum value
    opts_vrad@trYMaxF			   = 18.0	; y-axis maximum value 

  ; Resources to overlay zero line
    opts_cont                             = True
    opts_cont@cnFillOn                    = False
    opts_cont@cnLineColor                 = "Black"
    opts_cont@cnLevelSelectionMode	  = "ExplicitLevels"
    opts_cont@cnLevels                    = (/0.0/)
    opts_cont@cnLineThicknessF            = 3.0
    opts_cont@cnInfoLabelOn               = False
    opts_cont@gsnDraw                     = False       ; Do not draw the plot
    opts_cont@gsnFrame                    = False       ; Do no advance the frame
    opts_cont@gsnContourZeroLineThicknessF = 3.0
    opts_cont@cnLineLabelsOn               = False      ; Turn off line labels
    opts_cont@tiMainString                = ""
    opts_cont@gsnLeftString               = ""
    opts_cont@gsnRightString              = ""

;=============================
; Plot the radial wind data 
;=============================

  ; Radial wind (LHS panel)
    panel0     = gsn_csm_contour(wks,mean_vrad0(:,0,it,:),opts_vrad) ; [vrad]
    cont0a     = gsn_csm_contour(wks,mean_vrad0(:,0,it,:),opts_cont) ; [zero contour]
    cont0      = ColorNegDashZeroPosContour(cont0a,"transparent","black","transparent")
    overlay(panel0, cont0a)

  ; Radial wind	(middle panel)
    panel1     = gsn_csm_contour(wks,mean_vrad0(:,1,it,:),opts_vrad) ; [vrad]
    cont1a     = gsn_csm_contour(wks,mean_vrad0(:,1,it,:),opts_cont) ; [zero contour]
    cont1      = ColorNegDashZeroPosContour(cont1a,"transparent","black","transparent")
    overlay(panel1, cont1a)

  ; Radial wind	(RHS panel)
    panel2     = gsn_csm_contour(wks,mean_vrad0(:,2,it,:),opts_vrad) ; [vrad]
    cont2a     = gsn_csm_contour(wks,mean_vrad0(:,2,it,:),opts_cont) ; [zero contour]
    cont2      = ColorNegDashZeroPosContour(cont2a,"transparent","black","transparent")
    overlay(panel2, cont2a)

  ; Read data into panels for plotting 
    if (full .eq. 2) then 
     panel(0) = panel0
     panel(1) = panel1
     panel(2) = panel2
    else if (full .eq. 1) then 
     panel(0) = panel0
     panel(1) = panel2
    else if (full .eq. 0) then 
     panel(0) = panel0
    end if 
    end if 
    end if

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    optsP                       = True			   ; Panel plot resources
    optsP@gsnFrame              = False			   ; Do not advance the frame
    optsP@gsnPanelLabelBar      = True       	           ; Panel labelbar
    optsP@pmLabelBarWidthF	= 0.6			   ; Labelbar width 
    optsP@pmLabelBarHeightF     = 0.15	     	           ; Labelbar height
    optsP@lbLabelFontHeightF	= 0.01	     	           ; Labelbar font height 
    optsP@lbLabelFont		= "Helvetica" 	           ; Choose font 
    optsP@lbPerimOn		= False	      	           ; Perimeter on label bar 

  ; Additional labelbar resources
    optsP@lbTitleString         = mean_vrad0@description+\
                                  " ("+mean_vrad0@units+")" ; Labelbar string
    optsP@lbTitlePosition	= "Bottom"     		    ; Title position
    optsP@lbTitleFontHeightF    = 0.0175                    ; Title font height
    optsP@lbBoxEndCapStyle	= "TriangleBothEnds"	    ; Labelbar style
    optsP@pmLabelBarOrthogonalPosF = -0.002		    ; Move lb up/down (h/l)
    optsP@lbTitleOffsetF        = 0		   	    ; x-axis title up/down (h/l)
    optsP@txString              = title_arr(it)		    ; Plot title string

    if (full .eq. 2) then 
     optsP@gsnPanelFigureStrings = (/"a) "+ens0+" ('GOOD')",\
                                     "b) "+ens1+" ('BAD')",\
                                     "c) "+ens2+" ('VC')"/)
     optsP@gsnPanelFigureStringsFontHeightF = 0.0125
    else if (full .eq. 1) then 
     optsP@gsnPanelFigureStrings = (/"a) "+ens0+" ('GOOD')",\
                                     "b) "+ens2+" ('VC')"/)
     optsP@gsnPanelFigureStringsFontHeightF = 0.0125
    else if (full .eq. 0) then 
     optsP@gsnPanelFigureStrings = (/"a) "+ens0/)
     optsP@gsnPanelFigureStringsFontHeightF = 0.015
    end if 
    end if
    end if 

    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.90
    optsP@gsnPanelBottom        = 0.10
    optsP@amJust                = "TopRight"

    if (full .eq. 2) then 
     gsn_panel(wks,panel,(/1,3/),optsP)             ; Draw as a single plot
     frame(wks)
    else if (full .eq. 1) then
     gsn_panel(wks,panel,(/1,2/),optsP)             ; Draw as a single plot
     frame(wks)
    else if (full .eq. 0) then 
     gsn_panel(wks,panel,(/1,1/),optsP)             ; Draw as a single plot
     frame(wks)
    end if 
    end if 
    end if 

  ; Tidy up
    delete([/panel0,panel1,panel2,cont0a,cont1a,cont2a,cont0,cont1,cont2/])

;================================
; Now plot the tangential wind 
;================================

  ; Output file location and type

    if (full .eq. 1) then
     output = "$sam/nepartak/images/cyl_coords/vtan_xz_"+ens0+"_"+ens1+"_"+time_arr(it)
    else
     output = "$sam/nepartak/images/cyl_coords/vtan_xz_"+ens0+"_"+ens1+\
              "_"+ens2+"_"+time_arr(it)
    end if

    wks = gsn_open_wks(opt,output)

  ; Set up panel plot

    if (full .eq. 2) then
     panel = new(3,graphic)
    else if (full .eq. 1) then 
     panel = new(2,graphic)
    else if (full .eq. 0) then
     panel = new(1,graphic)
    end if
    end if
    end if 

  ; Resources to overlay tangential wind as rough marker of eyewall (following M18)
    opts_vtan                              = opts_vrad
    delete(opts_vtan@cnLevels)
    delete(opts_vtan@cnFillColors)

    gsn_define_colormap(wks,"radar_new")
    cmap = read_colormap_file("radar_new")
    cmap(0,3) = 0.0

    opts_vtan@cnLineLabelInterval         = 2.0
    opts_vtan@cnLineLabelFontHeightF      = 0.012
    opts_vtan@cnLineLabelBackgroundColor  = "transparent"
    opts_vtan@cnLineLabelPlacementMode    = "constant"
    opts_vtan@cnLinesOn                   = False ; Contour lines off
    opts_vtan@cnInfoLabelOn               = False
    opts_vtan@cnFillPalette               = cmap
    opts_vtan@cnLevelSelectionMode        = "ExplicitLevels"
    opts_vtan@cnLevels                    = (/5.0, 10.0, 15.0, 20.0,\
                                              25.0, 30.0, 35.0, 40.0, \
                                              45.0, 50.0, 55.0, 60.0, \
                                              65.0/)
    opts_vtan@cnFillColors                = (/0,2,3,4,5,6,7,8,9,10,11,12,13,14/)

  ; Tangential wind (LHS panel)
    panel0     = gsn_csm_contour(wks,mean_vtan0(:,0,it,:),opts_vtan) ; [vtan]

  ; Tangential wind (middle panel)
    panel1     = gsn_csm_contour(wks,mean_vtan0(:,1,it,:),opts_vtan) ; [vtan]

  ; Tangential wind (RHS panel)
    panel2     = gsn_csm_contour(wks,mean_vtan0(:,2,it,:),opts_vtan) ; [vtan]

  ; Read data into panels for plotting
    if (full .eq. 2) then
     panel(0) = panel0
     panel(1) = panel1
     panel(2) = panel2
    else if (full .eq. 1) then 
     panel(0) = panel0
     panel(1) = panel2
    else if (full .eq. 0) then
     panel(0) = panel0
    end if
    end if
    end if 

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    optsP                       = True                     ; Panel plot resources
    optsP@gsnFrame              = False                    ; Do not advance the frame
    optsP@gsnPanelLabelBar      = True                     ; Panel labelbar
    optsP@pmLabelBarWidthF      = 0.6                      ; Labelbar width
    optsP@pmLabelBarHeightF     = 0.15                     ; Labelbar height
    optsP@lbLabelFontHeightF    = 0.01                     ; Labelbar font height
    optsP@lbLabelFont           = "Helvetica"              ; Choose font
    optsP@lbPerimOn             = False                    ; Perimeter on label bar

  ; Additional labelbar resources
    optsP@lbTitleString         = mean_vtan0@description+\
                                  " ("+mean_vtan0@units+")" ; Labelbar string
    optsP@lbTitlePosition       = "Bottom"                 ; Title position
    optsP@lbTitleFontHeightF    = 0.0175                    ; Title font height
    optsP@lbBoxEndCapStyle      = "TriangleBothEnds"       ; Labelbar style
    optsP@pmLabelBarOrthogonalPosF = -0.002                ; Move lb up/down (h/l)
    optsP@lbTitleOffsetF        = 0                        ; x-axis title up/down (h/l)
    optsP@txString              = title_arr(it)            ; Plot title string

    if (full .eq. 2) then
     optsP@gsnPanelFigureStrings = (/"a) "+ens0+" ('GOOD')",\
                                     "b) "+ens1+" ('BAD')",\
                                     "c) "+ens2+" ('VC')"/)
     optsP@gsnPanelFigureStringsFontHeightF = 0.0125
    else if (full .eq. 1) then
     optsP@gsnPanelFigureStrings = (/"a) "+ens0+" ('GOOD')",\
                                     "b) "+ens2+" ('VC')"/)
     optsP@gsnPanelFigureStringsFontHeightF = 0.0125
    else if (full .eq. 0) then
     optsP@gsnPanelFigureStrings = (/"a) "+ens0/)
     optsP@gsnPanelFigureStringsFontHeightF = 0.015
    end if
    end if
    end if

    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.90
    optsP@gsnPanelBottom        = 0.10
    optsP@amJust                = "TopRight"

    if (full .eq. 2) then
     gsn_panel(wks,panel,(/1,3/),optsP)             ; Draw as a single plot
     frame(wks)
    else if (full .eq. 1) then 
     gsn_panel(wks,panel,(/1,2/),optsP)
     frame(wks)
    else if (full .eq. 0) then
     gsn_panel(wks,panel,(/1,1/),optsP)             ; Draw as a single plot
     frame(wks)
    end if
    end if
    end if 

  ; Tidy up
    delete([/panel0,panel1,panel2/])

;===================================
; Now plot the vertical velocity
;===================================

  ; Output file location and type

    if (full .eq. 1) then
     output = "$sam/nepartak/images/cyl_coords/vvel_xz_"+ens0+"_"+ens1+"_"+time_arr(it)
    else
     output = "$sam/nepartak/images/cyl_coords/vvel_xz_"+ens0+"_"+ens1+\
              "_"+ens2+"_"+time_arr(it)
    end if

    wks = gsn_open_wks(opt,output)

  ; Set up panel plot

    if (full .eq. 2) then
     panel = new(3,graphic)
    else if (full .eq. 1) then 
     panel = new(2,graphic)
    else if (full .eq. 0) then
     panel = new(1,graphic)
    end if
    end if
    end if 

  ; Vertical velocity resources
    opts_vvel                             = opts_vtan
    delete(opts_vvel@cnLevels)
    delete(opts_vvel@cnFillColors)

    gsn_define_colormap(wks,"n11")

    opts_vvel@cnLineLabelInterval         = 2.0
    opts_vvel@cnLineLabelFontHeightF      = 0.012
    opts_vvel@cnLineLabelBackgroundColor  = "transparent"
    opts_vvel@cnLineLabelPlacementMode    = "constant"
    opts_vvel@cnLinesOn                   = False ; Contour lines off
    opts_vvel@cnInfoLabelOn               = False
    opts_vvel@cnLevelSelectionMode        = "ExplicitLevels"
    opts_vvel@cnLevels                    = (/0.1, 0.2, 0.3, 0.4, 0.5, 0.6, \
                                              0.7, 0.8, 1.0, 1.2, 1.5/)
    opts_vvel@cnFillColors                = (/0,2,3,4,5,6,7,8,9,10,11,12/)

  ; Vertical velocity (LHS panel)
    panel0     = gsn_csm_contour(wks,mean_vvel0(:,0,it,:),opts_vvel) ; [vvel]

  ; Vertical velocity (middle panel)
    panel1     = gsn_csm_contour(wks,mean_vvel0(:,1,it,:),opts_vvel) ; [vvel]

  ; Vertical velocity (RHS panel)
    panel2     = gsn_csm_contour(wks,mean_vvel0(:,2,it,:),opts_vvel) ; [vvel]

  ; Read data into panels for plotting
    if (full .eq. 2) then
     panel(0) = panel0
     panel(1) = panel1
     panel(2) = panel2
    else if (full .eq. 1) then 
     panel(0) = panel0
     panel(1) = panel2
    else if (full .eq. 0) then
     panel(0) = panel0
    end if
    end if
    end if 

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    optsP                       = True                     ; Panel plot resources
    optsP@gsnFrame              = False                    ; Do not advance the frame
    optsP@gsnPanelLabelBar      = True                     ; Panel labelbar
    optsP@pmLabelBarWidthF      = 0.6                      ; Labelbar width
    optsP@pmLabelBarHeightF     = 0.15                     ; Labelbar height
    optsP@lbLabelFontHeightF    = 0.01                     ; Labelbar font height
    optsP@lbLabelFont           = "Helvetica"              ; Choose font
    optsP@lbPerimOn             = False                    ; Perimeter on label bar

  ; Additional labelbar resources
    optsP@lbTitleString         = mean_vvel0@description+\
                                  " ("+mean_vvel0@units+")" ; Labelbar string
    optsP@lbTitlePosition       = "Bottom"                  ; Title position
    optsP@lbTitleFontHeightF    = 0.0175                    ; Title font height
    optsP@lbBoxEndCapStyle      = "TriangleBothEnds"        ; Labelbar style
    optsP@pmLabelBarOrthogonalPosF = -0.002                 ; Move lb up/down (h/l)
    optsP@lbTitleOffsetF        = 0                         ; x-axis title up/down (h/l)
    optsP@txString              = title_arr(it)             ; Plot title string

    if (full .eq. 2) then
     optsP@gsnPanelFigureStrings = (/"a) "+ens0+" ('GOOD')",\
                                     "b) "+ens1+" ('BAD')",\
                                     "c) "+ens2+" ('VC')"/)
     optsP@gsnPanelFigureStringsFontHeightF = 0.0125
    else if (full .eq. 1) then
     optsP@gsnPanelFigureStrings = (/"a) "+ens0+" ('GOOD')",\
                                     "b) "+ens2+" ('VC')"/)
     optsP@gsnPanelFigureStringsFontHeightF = 0.0125
    else if (full .eq. 0) then
     optsP@gsnPanelFigureStrings = (/"a) "+ens0/)
     optsP@gsnPanelFigureStringsFontHeightF = 0.015
    end if
    end if
    end if

    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.90
    optsP@gsnPanelBottom        = 0.10
    optsP@amJust                = "TopRight"

    if (full .eq. 2) then
     gsn_panel(wks,panel,(/1,3/),optsP)             ; Draw as a single plot
     frame(wks)
    else if (full .eq. 1) then 
     gsn_panel(wks,panel,(/1,2/),optsP)
     frame(wks)
    else if (full .eq. 0) then
     gsn_panel(wks,panel,(/1,1/),optsP)             ; Draw as a single plot
     frame(wks)
    end if
    end if
    end if

  ; Tidy up
    delete([/panel0,panel1,panel2/])

  end do     ; End time loop (do it = 0, ts, tf)

end 