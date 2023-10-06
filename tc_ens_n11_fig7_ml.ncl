; Script to reproduce Figure 6 from Nguyen et al. (2011)

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_rm.ncl"

; ncl 'dat="02T12"' 'opt="x11"' 'ens0="em11"' dist=3.0 ts=42 tf=86 mlev=13
; rclr=0 cnt=1 lay=0 mlev2=24 ar=1.0 nr=21 ar0=1.0 sr=1 tc_ens_n11_fig7_ml.ncl

; 'dat'    = initialisation time string ("02T12", "03T00", ...)
; 'ens0'   = ensemble simulation (em00-em11)
; 'dist'   = size of box following storm (degrees)
; 'opt'    = output file format ("pdf" or "x11")
; 'rclr'   = colour scheme for plotting radial gradient of relative vorticity (0-2)
; 'cnt'	   = use offline (1) or in-script (0) calculations of storm-centred grid 
; 'mlev0'  = model level for plotting relative velocity
; 'lay'    = calculate layer average (1) or single level (0)
; 'mlev1'  = upper model level for calculating layer average (only valid when lay=1)
; 'ar'	   = distance in degrees of outer radius (1.0, 2.0, etc)
; 'nr'	   = number of radial circles between r = 0 and r = 'ar' (21, 41, 61, etc)
; 'ar0'	   = outer radius (degrees) for final plot [must be <= 'ar']
; 'sr'	   = calculate storm-relative (1) or Earth-relative (0) vorticity

; 'mlev0'  = 13 and 'mlev1' = 24 give an approximate 1-3 km layer average
 
begin

;==============================
; Prelash (reading in files)
;==============================

; Find 'pb' stream data to read in and plot
  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0
  fili_c1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pc*.nc")
  fili_c      = diri+"/"+fili_c1
  fili_d1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pd*.nc")
  fili_d      = diri+"/"+fili_d1

; Calculate number of times in file using built-in function ('nc_times')
  numINFO_c   = nc_times(fili_c)
  numINFO_d   = nc_times(fili_d)
  numINPUT_c  = numINFO_c[0]            ; Number of times in 'pc' stream
  time_c      = numINFO_c[1]            ; Array of times from 'pc' stream file
  numINPUT_d  = numINFO_d[0]            ; Number of times in 'pd' stream
  time_d      = numINFO_d[1]            ; Array of times from 'pd' stream file

; Calculate storm motion using built-in function ('st_rm')
; 'dat'     = initialisation time (02T12, 03T00, etc)
; 'ens0'    = ensemble simulation (em00-em11)
; 'diri'    = path to input files (see above)
; 'dist'    = size of box following storm (degrees)
; 'mins'    = analyse 1-h (0) or 5-min (1) data
; 'dir'     = T+00 to T+53 (0), T+54 to T+66 (1), or T+67 to T+83 (2) [TEMPORARY]

  storm_rel = st_rm(dat, ens0, "$ar/text/", dist, 0, 0)
  u_cyc     = storm_rel[0]             ; Zonal wind
  v_cyc     = storm_rel[1]             ; Meridional wind
  vel_cyc   = storm_rel[2]             ; Vector wind
  numTIMES  = storm_rel[3]             ; Number of times in file
  lat_arr   = storm_rel[4]             ; Latitude array subset (following storm)
  lon_arr   = storm_rel[5]             ; Longitude array subset (following storm)
  centre    = storm_rel[6]             ; Storm track (position) array

; Create arrays for data and output information
  llbox     = toint(dist*50)                   ; Calculate domain size (grid points)
  dsize     = (/llbox,llbox/)                  ; Domain size (grid points)

; Extract string for each ensemble member
  str1    = str_split(fili_c(0),"_")
  ens_str = str1(5)
  delete(str1)

  print("Working on ensemble member: "+ens_str)

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 10000000000
  end setvalues

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

;========================================
; Define cylindrical coordinate arrays
;========================================

; ar = distance (degrees) of outer radial circle (1.0)
; nr = total number of radial circles (21)

  thetas     = new(73,float)
  pii        = 3.14159265
  radii      = fspan(0.0,ar,nr)      ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,73)   ; Array of azimuth angles (0, 5, 10, etc)
  thetas_rad = (pii/180)*thetas      ; Degrees to radians
  rad_size   = dimsizes(radii)-1     ; Index for later calculation

; Define arrays for later use in tangential velocity cross-section plots
  vtan_azi   = new((/dimsizes(radii)/),float) ; Tangential wind
  vrad_azi   = new((/dimsizes(radii)/),float) ; Radial wind
  spd_azi    = new((/dimsizes(radii)/),float) ; Vector wind
  pvrt_azi   = new((/dimsizes(radii)/),float) ; Potential vorticity
  vort_azi   = new((/dimsizes(radii)/),float) ; Relative vorticity
  vvel_azi   = new((/dimsizes(radii)/),float) ; Vertical velocity

  nTIMES     = (tf - ts) + 1		      ; Number of times in user-defined subset

; Define 2-dimensional array to hold all values from 'vtan_azi' at ALL times
  vtan_plot0   = new((/dimsizes(radii),numTIMES/),float)
  vtan_plot0!0 = "rad"
  vtan_plot0!1 = "time"

  vrad_plot0   = new((/dimsizes(radii),numTIMES/),float)
  vrad_plot0!0 = "rad"
  vrad_plot0!1 = "time"

  pvrt_plot0   = new((/dimsizes(radii),numTIMES/),float)
  pvrt_plot0!0 = "rad"
  pvrt_plot0!1 = "time"

  vort_plot0   = new((/dimsizes(radii),numTIMES/),float)
  vort_plot0!0 = "rad"
  vort_plot0!1 = "time"

  vvel_plot0   = new((/dimsizes(radii),numTIMES/),float)
  vvel_plot0!0 = "rad"
  vvel_plot0!1 = "time"

; Add metadata for azimuthally averaged variables
  vtan_azi!0           = "rad"
  vtan_azi&rad         = radii
  vtan_azi@description = "Azimuthally averaged tangential velocity"
  vtan_azi@units       = "m s~S~-1~N~"

  vrad_azi!0           = "rad"
  vrad_azi&rad         = radii
  vrad_azi@description = "Azimuthally averaged radial velocity"
  vrad_azi@units       = "m s~S~-1~N~"

  spd_azi!0            = "rad"
  spd_azi&rad          = radii
  spd_azi@description  = "Azimuthally averaged vector velocity"
  spd_azi@units        = "m s~S~-1~N~"

  pvrt_azi!0           = "rad"
  pvrt_azi&rad         = radii
  pvrt_azi@description = "Azimuthally averaged PV"
  pvrt_azi@units       = "10~S~-6~N~K kg~S~-1~N~ m~S~2~N~ s~S~-1~N~"

  vort_azi!0           = "rad"
  vort_azi&rad         = radii
  vort_azi@description = "Azimuthally averaged relative vorticity"
  vort_azi@units       = "s~S~-1~N~"

  vvel_azi!0            = "rad"
  vvel_azi&rad          = radii
  vvel_azi@description  = "Azimuthally averaged vertical velocity"
  vvel_azi@units        = "m s~S~-1~N~"

; Create time array for later plotting 
  time_int = new(numTIMES,"integer")
  time_int = ispan(0,119,1) ; ispan(ts,tf,1)
  time_int!0 = "time"

;===============================
; Start loop over input files
;===============================

  c  = addfile(fili_c,"r")             ; Read in 'pc' stream [u,v,vort]
  d  = addfile(fili_d,"r")             ; Read in 'pd' stream [w,p,t,theta,kmh,tau]

  ct = 0			       ; Counter variable (time)

;==================================
; Get the variables we will need
;==================================

  do it = ts, tf                   ; Loop over times in file

    print("Working on time: "+time_str(it))
    itt = it
    title_arr(it) = "Valid at "+time_str(it)+" (T+"+itt+")"

;=====================================================
; Read in variables from 'pc' stream (model levels)
;=====================================================

    t0    = lat_arr(it+1,0)
    t1    = lat_arr(it+1,1)
    n0    = lon_arr(it+1,0)
    n1    = lon_arr(it+1,1)
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
   if (lay .eq. 1) then
    u0  = c->$u_varname$(it,mlev0:mlev1,{t0:t1},{n0:n1})     ; Zonal wind (layer)
    v0  = c->$v_varname$(it,mlev0:mlev1,{t0:t1},{n0:n1})     ; Meridional wind (layer)
   else
    u0  = c->$u_varname$(it,mlev0,{t0:t1},{n0:n1})           ; Zonal wind (single level)
    v0  = c->$v_varname$(it,mlev0,{t0:t1},{n0:n1})           ; Meridional wind (single level)
   end if

;===================================
; Calculate storm-relative winds
;===================================

    u   = u0 - u_cyc(it)
    v   = v0 - v_cyc(it)

  ; Add metadata from original horizontal wind arrays
    copy_VarCoords(u0,u)
    copy_VarAtts(u0,u)

    copy_VarCoords(v0,v)
    copy_VarAtts(v0,v)

;==================================================================
; Calculate relative vorticity using centered finite differences
;==================================================================

  ; Option '2' --> boundary points estimated using one-sided difference scheme

    if (sr .eq. 1) then
     vort0 = uv2vr_cfd(u, v, latC, lonC, 2)             ; Storm-relative vorticity
     copy_VarCoords(u, vort0)
    elseif (sr .eq. 0) then
     vort0 = uv2vr_cfd(u0, v0, latC, lonC, 2)           ; Earth-relative vorticity
     copy_VarCoords(u0, vort0)
    end if

    vort0@units = "s~S~-1~N~"
    vort0@name  = "Relative vorticity on model levels"

    if (lay .eq. 1) then
     vort       = vort0(:,0:dsize(0)-1,0:dsize(1)-1)
    else
     vort       = vort0(0:dsize(0)-1,0:dsize(1)-1)
    end if

  ; Coriolis parameter and absolute vorticity
    f1   = coriolis_param(latC(0:dsize(0)-1))

    if (lay .eq. 0) then
     f0   = conform_dims(dimsizes(vort),f1,1)
    elseif (lay .eq. 1) then
     f0   = conform_dims(dimsizes(vort),f1,2)
    end if

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
    lon1D  = d->longitude_1({n0:n1})                   ; longitude [xxx grid points]
    lat1D  = d->latitude_1({t0:t1})                    ; latitude  [xxx grid points]
    lonD   = d->longitude({n0:n1})                     ; longitude [xxx grid points]
    latD   = d->latitude({t0:t1})                      ; latitude  [xxx grid points]
    hybD   = d->hybrid_ht(:)                           ; 63 model ('rho') levels

    hy     = hybD / 1000 ; Model level values in km
    hy0    = sprintf("%0.0f",hybD(mlev0))
    print("Working on model level: "+hy0+ " m")

    if (lay .eq. 1) then
     hy1    = sprintf("%0.0f",hybD(mlev1))
    end if

  ; Read in temperature, pressure, theta, vertical velocity and geopotential height
  ; All variables below -- [hybD | 63] * [latD | 150] * [lon1D | 151]

    if (lay .eq. 1) then 
     t      = d->temp(it-2,mlev0:mlev1,{t0:t1},{n0:n1})        ; Temperature (K)
     p      = d->p(it-2,mlev0:mlev1,{t0:t1},{n0:n1})           ; Pressure (Pa)
     th     = d->theta(it-2,mlev0:mlev1,{t0:t1},{n0:n1})       ; Potential temperature (K)
     w      = d->dz_dt(it-2,mlev0:mlev1,{t0:t1},{n0:n1})       ; Vertical velocity (m s-1)
     z      = d->ht(it-2,mlev0:mlev1,{t0:t1},{n0:n1})          ; Geopotential height (m)
     pv	   = d->field83(it-2,mlev0:mlev1,{t0:t1},{n0:n1})      ; Potential vorticity (PVU)
    elseif (lay .eq. 0) then
     t      = d->temp(it-2,mlev0,{t0:t1},{n0:n1})  	       ; Temperature (K)
     p      = d->p(it-2,mlev0,{t0:t1},{n0:n1})                 ; Pressure (Pa)
     th     = d->theta(it-2,mlev0,{t0:t1},{n0:n1})             ; Potential temperature (K)
     w      = d->dz_dt(it-2,mlev0,{t0:t1},{n0:n1})             ; Vertical velocity (m s-1)
     z      = d->ht(it-2,mlev0,{t0:t1},{n0:n1})                ; Geopotential height (m)
     pv     = d->field83(it-2,mlev0,{t0:t1},{n0:n1})           ; Potential vorticity (PVU)
    end if 

  ; Convert units of pressure and geopotential height before continuing
    p  = p / 100           ; Convert to 'hPa'
    z  = z / 10		   ; Convert to 'dam'
    pv = pv * (10 ^ 6)	   ; Convert to 'PVU' [10 ^ 6 K m2 s-1 kg-1]

  ; Add metadata
    p@units  = "hPa"
    z@units  = "dam"
    pv@units = "PVU"

    dx       = lonD(1) - lonD(0) ; Grid spacing (longitude)
    dy       = latD(1) - latD(0) ; Grid spacing (latitude)

;==============================================
; Analyse variables on chosen pressure level
;==============================================

    if (lay .eq. 1) then 
     u_plane       = dim_avg_n_Wrap(u(:,:,:),0)
     v_plane       = dim_avg_n_Wrap(v(:,:,:),0)
     vort_plane    = dim_avg_n_Wrap(vort(:,:,:),0)
     prs_plane     = dim_avg_n_Wrap(p(:,:,:),0)
     pv_plane	   = dim_avg_n_Wrap(pv(:,:,:),0)
     w_plane	   = dim_avg_n_Wrap(w(:,:,:),0)
    elseif (lay .eq. 0) then 
     u_plane       = u(:,:)
     v_plane       = v(:,:)
     vort_plane    = vort(:,:)
     prs_plane     = p(:,:)
     pv_plane      = pv(:,:)
     w_plane       = w(:,:)
    end if 

  ; Max and min pressure on height level
    pmin          = toint(min(prs_plane))
    pmax          = toint(max(prs_plane))

  ; Calculate maximum relative vorticity on pressure level
    vort_max      = max(vort_plane)
    vort_smth     = smth9_Wrap(vort_plane, 0.5, 0.5, True)
    vort_max_smth = max(vort_smth)

  ; Find index of relative vorticity maximum
    dims          = dimsizes(vort_plane)
    vort1d        = ndtooned(vort_plane)
    inds          = ind_resolve(maxind(vort1d),dims)
    vort1d_sm     = ndtooned(vort_smth)
    inds_sm       = ind_resolve(maxind(vort1d_sm),dims)

    lat_max 	  = t0 + (dy * inds(0,0))
    lon_max 	  = n0 + (dx * inds(0,1))

    print("ct = "+ct+" ; it = "+it)
    print("lat_max = "+centre(it-1,0)+" ; lon_max = "+centre(it-1,1))

  ; Difference between 'lat_max'/'lon_max' and values read in using NCL function
  ; 'it+1' results from difference in number of times in 'pb' and 'pc' stream arrays
    dlat          = abs(lat_max - centre(it-1,0))
    dlon          = abs(lon_max - centre(it-1,1))

  ; If our calculations of latitude and longitude deviate too much from the values
  ; we read in using the NCL function ('st_rm'), then calculate a new centre
  ; using the minimum pressure. Constrain the new calculation using the values
  ; from 'st_rm', to prevent spurious centres being identified and analysed.

    if (dlat .gt. 0.75 .or. dlon .gt. 0.75) then

   ; Should this be 'it-1' because pressure is in the 'pd' stream files?
     lt1 = centre(it-1,0) - 0.5
     lt2 = centre(it-1,0) + 0.5
     ln1 = centre(it-1,1) - 0.5
     ln2 = centre(it-1,1) + 0.5

     prs_plane0 = prs_plane({lt1:lt2},{ln1:ln2})       ; Smaller grid (1.0 x 1.0 deg)
     prs_min0   = min(prs_plane0)                      ; Find minimum on smaller grid
     dims_p0    = dimsizes(prs_plane0)                 ; Size of smaller grid
     prs_1d0    = ndtooned(prs_plane0)                 ; Create 1-D array
     inds_p0    = ind_resolve(minind(prs_1d0),dims_p0) ; Find index of minimum

     lat_max    = lt1 + (dy * inds_p0(0,0))
     lon_max    = ln1 + (dx * inds_p0(0,1))
     print("New centre: "+lat_max+" degrees N, "+lon_max+" degrees E")

    end if

   ;==============================================================
   ; Create arrays for later switch to cylindrical coordinates
   ;==============================================================

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
     spd_int  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vector wind (for RMW)
     pv_int   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Potential vorticity
     vort_int = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Relative vorticity
     w_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vertical velocity
     v_rad    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Radial wind
     v_tan    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Tangential wind

   ; Add variable metadata
     copy_VarAtts(u_plane,u_int)
     u_int@description = "Zonal wind"
     u_int@units       = "m s~S~-1~N~"

     copy_VarAtts(v_plane,v_int)
     v_int@description = "Meridional wind"
     v_int@units       = "m s~S~-1~N~"

     copy_VarAtts(w_plane,w_int)
     v_int@description = "Vertical velocity"
     v_int@units       = "m s~S~-1~N~"

     copy_VarAtts(pv_plane,pv_int)
     pv_int@description = "Potential vorticity"
     pv_int@units       = "10~S~-6~N~K kg~S~-1~N~ m~S~2~N~ s~S~-1~N~"

     copy_VarAtts(vort_plane,vort_int)
     vort_int@description = "Relative vorticity"
     vort_int@units       = "s~S~-1~N~"

     v_rad@description = "Radial wind"
     v_rad@units       = "m s~S~-1~N~"

     v_tan@description = "Tangential wind"
     v_tan@units       = "m s~S~-1~N~"

     xpos!0    = "rad"
     xpos!1    = "azi"
     xpos&rad  = radii
     xpos&azi  = thetas

     ypos!0    = "rad"
     ypos!1    = "azi"
     ypos&rad  = radii
     ypos&azi  = thetas

     u_int!0   = "rad"
     u_int!1   = "azi"
     u_int&rad = radii
     u_int&azi = thetas

     v_int!0   = "rad"
     v_int!1   = "azi"
     v_int&rad = radii
     v_int&azi = thetas

     spd_int!0   = "rad"
     spd_int!1   = "azi"
     spd_int&rad = radii
     spd_int&azi = thetas

     w_int!0   = "rad"
     w_int!1   = "azi"
     w_int&rad = radii
     w_int&azi = thetas

     pv_int!0   = "rad"
     pv_int!1   = "azi"
     pv_int&rad = radii
     pv_int&azi = thetas

     vort_int!0   = "rad"
     vort_int!1   = "azi"
     vort_int&rad = radii
     vort_int&azi = thetas

     v_tan!0   = "rad"
     v_tan!1   = "azi"
     v_tan&rad = radii
     v_tan&azi = thetas

     v_rad!0   = "rad"
     v_rad!1   = "azi"
     v_rad&rad = radii
     v_rad&azi = thetas

   ; Latitude/longitude (small grid)

     if (cnt .eq. 1) then

      do irad = 0, dimsizes(radii)-1
        r = radii(irad)
        do iang = 0, dimsizes(thetas)-1
          theta = thetas_rad(iang)
            xpos(irad,iang) = centre(it-1,1) + r*cos(theta)
            ypos(irad,iang) = centre(it-1,0) + r*sin(theta)
        end do
      end do

     elseif (cnt .eq. 0) then 

      do irad = 0, dimsizes(radii)-1
        r = radii(irad)
        do iang = 0, dimsizes(thetas)-1
          theta = thetas_rad(iang)
            xpos(irad,iang) = lon_max + r*cos(theta)
            ypos(irad,iang) = lat_max + r*sin(theta)
        end do
      end do

     end if 

   ; Fill new arrays with interpolated values on the cylindrical grid
   ; Which latitude/longitude grids are each of these variables defined on?
   ; u_plane     -- latC,  lonC
   ; v_plane     -- latC,  lonC
   ; vort_plane  -- latC,  lonC

   ; w_plane     -- latD,  lon1D
   ; pv_plane	 -- latD,  lon1D

  ; Loop over azimuth angles to fill entire circular array of data points
    do iang = 0, dimsizes(thetas_rad)-1
      u_int(:,iang)    = linint2_points(lonC_sub,latC_sub,u_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      v_int(:,iang)    = linint2_points(lonC_sub,latC_sub,v_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      w_int(:,iang)    = linint2_points(lon1D_sub,latD_sub,w_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      pv_int(:,iang)   = linint2_points(lon1D_sub,latD_sub,pv_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      vort_int(:,iang) = linint2_points(lonC_sub,latC_sub,vort_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
    end do
 
  ; Calculate vector windspeed (from 'u' and 'v')
    spd_int = sqrt((u_int ^ 2) + (v_int ^ 2))
    copy_VarCoords(u_int,spd_int)
    copy_VarMeta(u_int,spd_int)
    spd_int@units = "m s~S~-1~N~"

  ; Calculate radial and tangential windspeed at each point
    do iang = 0, dimsizes(thetas)-1
      theta = thetas_rad(iang)
        
	ra = u_int(:,iang)*cos(theta)
	rb = v_int(:,iang)*sin(theta)

	ta = -u_int(:,iang)*sin(theta)
	tb = v_int(:,iang)*cos(theta)

        v_rad(:,iang) = ra + rb 
        v_tan(:,iang) = ta + tb

    end do

  ;===================================================
  ; Average variables azimuthally around the storm 
  ;===================================================

    do irad = 0, dimsizes(radii)-1
      r = radii(irad)
      vtan_azi(irad) = avg(v_tan(irad,:))        ; Tangential wind
      vrad_azi(irad) = avg(v_rad(irad,:))	 ; Radial wind 
      spd_azi(irad)  = avg(spd_int(irad,:))	 ; Vector wind (for RMW)
      pvrt_azi(irad) = avg(pv_int(irad,:))	 ; Potential vorticity
      vort_azi(irad) = avg(vort_int(irad,:))	 ; Relative vorticity
      vvel_azi(irad) = avg(w_int(irad,:))	 ; Vertical velocity
    end do

  ; Input azimuthal averages into larger arrays before plotting  
    vtan_plot0(:,it-1)  = vtan_azi(:)            ; Tangential wind
    vrad_plot0(:,it-1)  = vrad_azi(:)	    	 ; Radial wind
    pvrt_plot0(:,it-1)  = pvrt_azi(:)	    	 ; PV
    vort_plot0(:,it-1)	= vort_azi(:)		 ; Relative vorticity
    vvel_plot0(:,it-1)  = vvel_azi(:)	    	 ; Vertical velocity

  ; Tidy up
    delete([/u_int,v_int,w_int,pv_int,vort_int,v_tan,v_rad/])
    delete([/u_plane,v_plane,w_plane,vort_plane,prs_plane,pv_plane/])

    ct  = ct + 1		    ; Counter variable (time)

  ; Tidy up (avoid dimension size errors)
    delete([/lonC,latC,lon1C,lat1C,u,v,vort/])
    delete([/lonD,latD,lon1D,lat1D,t,p,w,z,pv/])

   end do     ; End time loop (do it = ts, tf)

;====================================================
; Calculate radial gradient of relative vorticity
;====================================================

; Initialise new array
  vort_tend0	  = vort_plot0

; Add metadata
  vort_tend0!0    = "rad"
  vort_tend0&rad  = radii
  vort_tend0!1    = "time"
  vort_tend0&time = time_int

  print("Calculating radial gradient of relative vorticity!")

  do rd = 0, rad_size

   rp1  = min((/rd+1,rad_size/))
   rm1  = max((/rd-1,0/))

   vort_tend0(rd,:) = (vort_plot0(rp1,:) - vort_plot0(rm1,:)) / 2 ; Units [m s-1 h-1]

  end do     ; End loop over radii (do rd = 0, rad_size)

;======================================
; Reorder dimensions before plotting 
;======================================

 ; Arrays now ordered [TIME * RADIUS] rather than [RADIUS * TIME]
   vort_tend = vort_tend0(time|:,rad|:)
   vort_plot = vort_plot0(time|:,rad|:)

 ; Convert to units better suited to plotting 
   vort_tend = vort_tend * (10 ^ 4)

;=====================================================
; Create a panel plot of Figs. 6a and 6b from N11 
;=====================================================

  ; Output file location and type
    if (lay .eq. 1) then 
     output = "$sam/nepartak/images/cyl_coords/n11_fig7_"+dat+"_"+ens0+"_"+hy0+"_"+hy1
    elseif (lay .eq. 0) then 
     output = "$sam/nepartak/images/cyl_coords/n11_fig7_"+dat+"_"+ens0+"_"+hy0
    end if 
    wks = gsn_open_wks(opt,output)

  ; Set up panel plot
    panel = new(1,graphic)
    gsn_define_colormap(wks,"n11")

;==========================
; Options for plotting 
;==========================

  ; Radial gradient of relative vorticity
    opts_vort                              = True
    opts_vort@cnFillOn                     = True
    opts_vort@cnLineLabelInterval          = 2.0
    opts_vort@cnLineLabelFontHeightF       = 0.012
    opts_vort@cnLineLabelBackgroundColor   = "transparent"
    opts_vort@cnLineLabelPlacementMode     = "constant"
    opts_vort@cnLinesOn                    = False ; Contour lines off
    opts_vort@cnInfoLabelOn                = False
    opts_vort@cnLevelSelectionMode         = "ExplicitLevels"

    if (rclr .eq. 0) then
     opts_vort@cnLevels                     = (/0.0, -0.5, -1.0, -1.5, -2.0, \
                                                -2.5, -3.0, -4.0, -5.0/)
     opts_vort@cnFillColors                 = (/12,11,10,9,8,7,6,5,4,0/)
    elseif (rclr .eq. 1) then
     opts_vort@cnLevels                     = (/0.0, -0.5, -1.0, -1.5, -2.0, \
                                                -3.0, -4.0, -5.0, -6.0/)
     opts_vort@cnFillColors                 = (/12,11,10,9,8,7,6,5,4,0/)
    elseif (rclr .eq. 2) then
     opts_vort@cnLevels                     = (/0.0, -1.0, -2.0, -3.0, -4.0, \
                                                -5.0, -6.0, -7.0, -8.0/)
     opts_vort@cnFillColors                 = (/12,11,10,9,8,7,6,5,4,0/)
    elseif (rclr .eq. 3) then 
     opts_vort@cnLevels                     = (/0.0, -1.0, -2.0, -3.0, -4.0, \
                                                -5.0, -6.0, -8.0, -10.0/)
     opts_vort@cnFillColors                 = (/12,11,10,9,8,7,6,5,4,0/)
    end if

    opts_vort@gsnPaperOrientation          = "landscape"
    opts_vort@tiMainString                 = ""
    opts_vort@tiMainFontHeightF            = 0.0125
    opts_vort@gsnLeftString                = ""
    opts_vort@gsnRightString               = ""
    opts_vort@gsnMaximize                  = True
    opts_vort@lbLabelBarOn                 = True                  ; Labelbar on/off
    opts_vort@lbBoxEndCapStyle             = "TriangleBothEnds"    ; Labelbar end shape
    opts_vort@pmLabelBarWidthF             = 0.32                  ; Labelbar width
    opts_vort@pmLabelBarHeightF            = 0.08                  ; Labelbar height
    opts_vort@lbLabelFontHeightF           = 0.0125                ; Labelbar font size
    opts_vort@lbLabelFont                  = "Helvetica"           ; Labelbar font
    opts_vort@lbPerimOn                    = False                 ; Perimeter on/off
    opts_vort@gsnDraw                      = False                 ; Do not draw plot
    opts_vort@gsnFrame                     = False                 ; Do not advance frame
    opts_vort@vpWidthF                     = 0.30                  ; Match width in N11
    opts_vort@vpHeightF                    = 0.50                  ; Match height in N11

  ; Additional plotting resources
    opts_vort@tiYAxisString                = "Time (h)"
    opts_vort@tiXAxisString                = "Radius (km)"
    opts_vort@tiXAxisFontHeightF           = 0.018
    opts_vort@trYLog                       = False
    opts_vort@trYMinF                      = ts      ; Starting point (T+...) for plot
    opts_vort@trYMaxF                      = tf-2    ; Ending point (T+...) for plot
    opts_vort@tmXTOn                       = "False" ; Turn off top x-axis TM
    opts_vort@tmYROn                       = "False" ; Turn off right y-axis TM
    opts_vort@gsnMaximize                  = True    ; Maximise plot size
    opts_vort@gsnAddCyclic                 = False
    opts_vort@tmXBMode                     = "Explicit"                ; Choose tm mode
    opts_vort@tmXBValues                   = (/0.0,0.5,1.0,1.5,2.0/)   ; Tm positions
    opts_vort@tmXBLabels                   = (/"0","50","100",\
    					       "150","200"/)	       ; Tm values
    opts_vort@tmXBMinorOn                  = False

    opts_vort@pmLabelBarOrthogonalPosF     = 0.19    ; Move lb up/down (higher/lower)
    opts_vort@tiXAxisOffsetYF              = 0.17    ; Move x-axis title up/down (h/l)

    opts_vort@trXMinF                      = 0.0                       ; Set x-axis min
    opts_vort@trXMaxF                      = ar0                       ; Set x-axis max
    opts_vort@tmYLMode                     = "Manual"
    opts_vort@tmYLTickStartF               = opts_vort@trYMinF         ; See L730
    opts_vort@tmYLTickEndF                 = opts_vort@trYMaxF         ; See L731
    opts_vort@tmYLTickSpacingF             = 6
    opts_vort@tmYLMinorOn                  = False
 
;========================
; Construct panel plot
;========================

  ; Radial gradient of relative vorticity
    vort1 = gsn_csm_contour(wks,vort_tend,opts_vort)    ; Radial gradient of relative vorticity
    panel(0) = vort1
 
;=============================================
; Finally, put everything together and plot 
;=============================================

    optsP                       = True
    optsP@gsnFrame              = False    ; Do not advance the frame
    optsP@gsnPanelLabelBar      = False    ; Turn on panel labelbar
    optsP@txString              = "Simulation "+ens0+" ("+dat+")"

    optsP@gsnPanelFigureStrings = "" ; (/"a) Radial gradient of relative vorticity"/)

    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.90
    optsP@gsnPanelBottom        = 0.10
    optsP@amJust                = "TopLeft"
    optsP@gsnPanelFigureStringsFontHeightF = 0.0099 ; Reduce label size (default 0.01)

    gsn_panel(wks,panel,(/1,1/),optsP)             ; Draw as a single plot
    frame(wks)

end 