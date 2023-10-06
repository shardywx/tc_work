; Script to reproduce Figure 6 from Nguyen et al. (2011), using 5-minute data

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/diff_avg.ncl"

; ncl dat=\"02T12\" opt=\"x11\" ens0=\"em11\" dist=3.0 ts=815 tf=971 mlev=15
; tclr=0 rclr=0 wclr=0 lay=1 mlev2=22 ave=1 expl=1 ar=1.0 nr=21 nt=73 cn0=\"slp\" calc=1
; tc5_ens_n11_fig6_ml.ncl

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'ens0'  = ensemble simulation (em00-em11)
; 'dist'  = size of box following storm (degrees)
; 'opt'   = output file format ("pdf" or "x11")
; 'mlev'  = model level for plotting PV and vertical velocity (see L106-108)
; 'tclr'  = colour scheme for plotting tangential wind (0-2)
; 'rclr'  = colour scheme for plotting radial wind (0-2)
; 'wclr'  = colour for plotting vertical velocity (0-2)
; 'lay'   = calculate layer average (1) or single level (0)
; 'ave'	  = apply running average to 5-min data (1) or leave alone (0)
; 'mlev2' = upper model level for calculating layer average (only valid when lay=1)
; 'expl'  = draw user-defined labels on y-axis (1) or leave alone (0)
 
; Model levels
; (0) 
; (8) 
; (14) 

begin

;==============================
; Prelash (reading in files)
;==============================

; Find data streams to read in and plot
  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0+"/vc"
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0

  fili_c1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pc.nc")
  fili_d1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pd.nc")
  fili_j1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pj.nc")
  fili_k1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pk.nc")

  fili_c      = diri+"/"+fili_c1
  fili_d      = diri+"/"+fili_d1
  fili_j      = diri+"/"+fili_j1
  fili_k      = diri+"/"+fili_k1

; Calculate number of times in file using built-in function ('nc_times') - should be same here
  numINFO_c   = nc_times(fili_c)
  numINFO_d   = nc_times(fili_d)
  numINFO_j   = nc_times(fili_j)
  numINFO_k   = nc_times(fili_k)
  numINPUT_c  = numINFO_c[0]            ; Number of times in 'pc' stream
  time_c      = numINFO_c[1]            ; Array of times from 'pc' stream file
  numINPUT_d  = numINFO_d[0]            ; Number of times in 'pd' stream
  time_d      = numINFO_d[1]            ; Array of times from 'pd' stream file
  numINPUT_j  = numINFO_j[0]            ; Number of times in 'pj' stream
  time_j      = numINFO_j[1]            ; Array of times from 'pj' stream file
  numINPUT_k  = numINFO_k[0]            ; Number of times in 'pk' stream
  time_k      = numINFO_k[1]            ; Array of times from 'pk' stream file

; Calculate storm motion using built-in function ('st_rm')
; 'dat'     = initialisation time (02T12, 03T00, ...)
; 'ens0'    = ensemble simulation (em00, em01, ...)
; 'diri'    = path to input files (see above)
; 'dist'    = size of box following storm (degrees)
; 'mins'    = analyse 1-h (0) or 5-min (1) data

  storm_rel = st_rm(dat, ens0, "$ar/text/", dist, 1)
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

; Increase system memory
  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 10000000000
  end setvalues

;=======================================                                              
; Choose which storm track to read in                                                 
;=======================================                                              

  if (cn0 .eq. "slp") then
   cn = 0
  elseif (cn0 .eq. "slpf") then
   cn = 1
  end if

;===========================================================================
; Create date/time string arrays for all times with user-defined function
;===========================================================================

 ; Call external function 'times_arr'
 ; Returns 'time_str' (12 UTC 4 July), 'time_arr' (), and 'title_arr' ()
   all_times  = times_arr(dat, ens0)

 ; Return strings e.g. '0600 UTC 04 July'
   time_str   = all_times[0]

 ; Return strings e.g. '04Jul_0600Z'
   time_arr   = all_times[1]

 ; Return title strings e.g. 'Valid at 0600 UTC 04 July'
   title_arr  = all_times[2]

;==================================================================================
;;;;;;;;;;;;;;;;;;;;;;;;;;; Start of main part of script ;;;;;;;;;;;;;;;;;;;;;;;;;;
;==================================================================================

; Define cylindrical coordinate arrays
  pii        = 3.14159265
  radii      = fspan(0.0,ar,nr)      ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,nt)   ; Array of azimuth angles (0, 5, 10, etc)
  thetas_rad = (pii/180)*thetas      ; Degrees to radians
  rad_size   = dimsizes(radii)-1     ; Index for later calculation

; Define arrays for later use in tangential velocity cross-section plots
  vtan_azi   = new((/dimsizes(radii)/),float) ; Tangential wind
  vrad_azi   = new((/dimsizes(radii)/),float) ; Radial wind
  spd_azi    = new((/dimsizes(radii)/),float) ; Vector velocity
  pvrt_azi   = new((/dimsizes(radii)/),float) ; Potential vorticity
  vvel_azi   = new((/dimsizes(radii)/),float) ; Vertical velocity

  nTIMES     = (tf - ts) + 1		      ; Number of times in user-defined subset

; Define 2-dimensional array to hold all values from 'vtan_azi' at ALL times
  vtan_plot0   = new((/dimsizes(radii),nTIMES/),float)
  vtan_plot0!0 = "rad"
  vtan_plot0!1 = "time"

  vrad_plot0   = new((/dimsizes(radii),nTIMES/),float)
  vrad_plot0!0 = "rad"
  vrad_plot0!1 = "time"

  pvrt_plot0   = new((/dimsizes(radii),nTIMES/),float)
  pvrt_plot0!0 = "rad"
  pvrt_plot0!1 = "time"

  vvel_plot0   = new((/dimsizes(radii),nTIMES/),float)
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

  vvel_azi!0            = "rad"
  vvel_azi&rad          = radii
  vvel_azi@description  = "Azimuthally averaged vertical velocity"
  vvel_azi@units        = "m s~S~-1~N~"

; Create time array for later plotting 
  time_int = new(nTIMES,"integer")
  time_int = ispan(ts,tf,1)
  time_int!0 = "time"

;==================================
; Get the variables we will need
;==================================

  c  = addfile(fili_c,"r")         ; Read in 'pc' stream [u,v,vort]
  d  = addfile(fili_d,"r")         ; Read in 'pd' stream [w,p]
  j  = addfile(fili_j,"r")         ; Read in 'pj' stream [t,theta,mixing ratio]
  k  = addfile(fili_k,"r")         ; Read in 'pk' stream [tau,PV,w]
  ct = 0			   ; Counter variable (time)

  do it = ts, tf                   ; Loop over times in file (T+54 to T+66, every 5 min)

   time  = c->t(it)                ; Time in file
   times = dimsizes(time)

   print_clock("Working on "+time_str(it))

;=====================================================
; Read in variables from 'pc' stream (model levels)
;=====================================================

  ; Get domain boundaries at each time from storm track (calculated off-line)
    t0    = lat_arr(cn,it,0)
    t1    = lat_arr(cn,it,1)
    n0    = lon_arr(cn,it,0)
    n1    = lon_arr(cn,it,1)
    print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

 ; Create latitude/longitude array subsets to reduce computation time
   lonC  = c->longitude({n0:n1})    ; longitude [160 grid points]
   latC  = c->latitude({t0:t1})     ; latitude  [160 grid points]
   hybC  = c->hybrid_ht(:)          ; 63 model ('theta') levels

   u_varname = "x-wind"
   v_varname = "y-wind"

 ; Read in horizontal wind components
 ; Both variables below -- [hybC | 63] * [latC | 150] * [lonC | 150]
   if (lay .eq. 1) then
    u0  = c->$u_varname$(it,mlev:mlev2,{t0:t1},{n0:n1})     ; Zonal wind (layer)
    v0  = c->$v_varname$(it,mlev:mlev2,{t0:t1},{n0:n1})     ; Meridional wind (layer)
   else
    u0  = c->$u_varname$(it,mlev,{t0:t1},{n0:n1})           ; Zonal wind (single level)
    v0  = c->$v_varname$(it,mlev,{t0:t1},{n0:n1})           ; Meridional wind (single level)
   end if

;===================================                                                 
; Calculate storm-relative winds                                                     
;===================================                                                 

    u   = u0 - u_cyc(cn,it)
    v   = v0 - v_cyc(cn,it)

  ; Add metadata from original horizontal wind arrays                                
    copy_VarCoords(u0,u)
    copy_VarAtts(u0,u)

    copy_VarCoords(v0,v)
    copy_VarAtts(v0,v)

;==================================================================                  
; Calculate relative vorticity using centered finite differences                     
;================================================================== 

  ; Option '2' --> boundary points estimated using one-sided difference scheme
    vort = uv2vr_cfd(u, v, latC, lonC, 2)
    copy_VarCoords(u,vort)
    vort@units = "s~S~-1~N~"
    vort@name = "Relative vorticity on model levels"

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

;=====================================================
; Read in variables from 'pd' stream (model levels)
;=====================================================

  ; Read in basic variables from 'pd' stream
    lonD   = d->longitude({n0:n1})                  ; longitude [xxx grid points]
    latD   = d->latitude({t0:t1})                   ; latitude  [xxx grid points]
    hybD   = d->hybrid_ht(:)                        ; 63 model ('rho') levels

    hy     = hybD / 1000 ; Model level values in km
    hy0    = sprintf("%0.0f",hybD(mlev))
    print("Working on model level: "+hy0+ " m")

    if (lay .eq. 1) then 
     hy1    = sprintf("%0.0f",hybD(mlev2))
    end if 

  ; Read in temperature, pressure, vertical velocity and geopotential height
  ; All variables below -- [hybD | 63] * [latD | 150] * [lon1D | 151]

    if (lay .eq. 1) then
     print("Averaging PV between "+hy0+" and "+hy1+ " m")
     p     = d->p(it,mlev:mlev2,{t0:t1},{n0:n1})     ; Pressure (Pa)
     z     = d->ht(it,mlev:mlev2,{t0:t1},{n0:n1})    ; Geopotential height (m)
    else
     p     = d->p(it,mlev,{t0:t1},{n0:n1})           ; Pressure (Pa)
     z     = d->ht(it,mlev,{t0:t1},{n0:n1})          ; Geopotential height (m)
    end if 

  ; Convert units of pressure and geopotential height before continuing
    p      = p / 100          ; Convert to 'hPa'
    z      = z / 10           ; Convert to 'dam'

  ; Add metadata
    p@units  = "hPa"
    z@units  = "dam"

    dx0     = lonD(1) - lonD(0) ; Grid spacing (longitude)
    dy0     = latD(1) - latD(0) ; Grid spacing (latitude)

;======================================================
; Read in variables from 'pk' stream (model levels)
;======================================================

  ; Read in basic variables from 'pd' stream
    lon1K  = k->longitude_1({n0:n1})                ; longitude [xxx grid points]
    lat1K  = k->latitude_1({t0:t1})                 ; latitude  [xxx grid points]
    lonK   = k->longitude({n0:n1})                  ; longitude [xxx grid points]
    latK   = k->latitude({t0:t1})                   ; latitude  [xxx grid points]
    hybK   = k->hybrid_ht(:)                        ; 63 model ('rho') levels
    hyb1K  = k->hybrid_ht_1(:)                      ; 80 model ('rho') levels

    hyK    = hybK / 1000 ; Model level values in km
    hyk0   = sprintf("%0.0f",hybK(mlev))

    if (lay .eq. 1) then 
     hyk1  = sprintf("%0.0f",hybK(mlev2))
    end if 

    if (lay .eq. 1) then
     print("Averaging PV between "+hy0+" and "+hy1+ " m")
     pv    = k->field83(it,mlev:mlev2,{t0:t1},{n0:n1})   ; Potential vorticity (PVU)
     w     = k->dz_dt(it,mlev:mlev2,{t0:t1},{n0:n1})     ; Vertical velocity (m s-1)
    else
     pv    = k->field83(it,mlev,{t0:t1},{n0:n1})         ; Potential vorticity (PVU)
     w     = k->dz_dt(it,mlev,{t0:t1},{n0:n1})           ; Vertical velocity (m s-1)
    end if

  ; Convert units of pressure and geopotential height before continuing
    pv     = pv * (10 ^ 6)    ; Convert to 'PVU' [10 ^ 6 K m2 s-1 kg-1]

  ; Add metadata
    pv@units = "PVU"

;===========================================================
; Calculate rate of strain (following Nguyen et al. 2011)
;===========================================================

    dim1 = "latitude"
    dim2 = "longitude"

  ; Calculate horizontal derivatives of 'u' (du_dx, du_dy)
    du    = grad_latlon_cfd(u,u&$dim1$,u&$dim2$,False,False)
    du_dy = du[0]    ; Latitudinal derivative
    du_dx = du[1]    ; Longitudinal derivative
    delete(du)

  ; Calculate horizontal derivatives of 'v' (dv_dx, dv_dy)
    dv    = grad_latlon_cfd(v,v&$dim1$,v&$dim2$,False,False)
    dv_dy = dv[0]    ; Latitudinal derivative
    dv_dx = dv[1]    ; Longitudinal derivative
    delete(dv)

  ; Calculate strain --> sqrt[ (du_dx - dv_dy)^2 + (dv_dx + du_dy)^2 ]
    str1a   = du_dx - dv_dy   ; Stretching deformation
    str2a   = dv_dx + du_dy   ; Shearing deformation
    str1    = (str1a) ^ 2
    str2    = (str2a) ^ 2

    if (lay .eq. 1) then
     str_plane0  = sqrt( str1(:,:,:) + str2(:,:,:) )
    else
     str_plane0  = sqrt( str1(:,:) + str2(:,:) )
    end if
    delete([/str1a,str2a/])

;===========================================
; Analyse variables on chosen model level
;===========================================

    if (lay .eq. 1) then 
     u_plane       = dim_avg_n_Wrap(u(:,:,:),0)
     v_plane       = dim_avg_n_Wrap(v(:,:,:),0)
     vort_plane    = dim_avg_n_Wrap(vort(:,:,:),0)
     prs_plane     = dim_avg_n_Wrap(p(:,:,:),0)
     pv_plane	   = dim_avg_n_Wrap(pv(:,:,:),0)
     w_plane	   = dim_avg_n_Wrap(w(:,:,:),0)
     str_plane	   = dim_avg_n_Wrap(str_plane0(:,:,:),0)
     copy_VarCoords(u(0,:,:),str_plane)
    elseif (lay .eq. 0) then 
     u_plane       = u(:,:)
     v_plane       = v(:,:)
     vort_plane    = vort(:,:)
     prs_plane     = p(:,:)
     pv_plane      = pv(:,:)
     w_plane       = w(:,:)
     str_plane	   = str_plane0
     copy_VarCoords(u,str_plane)
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

  ;=================================================================
  ; Either calculate the centre, or read in pre-calculated values
  ;=================================================================

  if (calc .eq. 0) then 

  ; Find latitude and longitude of storm centre using vorticity calculation above 
    lat_max0 	  = t0 + (dy0 * inds(0,0))
    lon_max0	  = n0 + (dx0 * inds(0,1))

    print("Centre (vort): "+lat_max0+" degrees N, "+lon_max0+" degrees E")

  ; Create smaller grid around location identified using vorticity maximum               
    lt1 = lat_max0 - 0.15
    lt2 = lat_max0 + 0.15
    ln1 = lon_max0 - 0.15
    ln2 = lon_max0 + 0.15

  ; Calculate minimum pressure in region surrounding vorticity maximum                   
    prs_plane0 = prs_plane({lt1:lt2},{ln1:ln2})       ; Smaller grid (1.0 x 1.0 deg)     
    prs_min0   = min(prs_plane0)                      ; Find minimum on smaller grid     
    dims_p0    = dimsizes(prs_plane0)                 ; Size of smaller grid             
    prs_1d0    = ndtooned(prs_plane0)                 ; Create 1-D array                 
    inds_p0    = ind_resolve(minind(prs_1d0),dims_p0) ; Find index of minimum            

    lat_max    = lt1 + (dy0 * inds_p0(0,0))
    lon_max    = ln1 + (dx0 * inds_p0(0,1))
    print("Centre (vort + pres): "+lat_max+" degrees N, "+lon_max+" degrees E")
    delete([/prs_plane0,prs_min0,dims_p0,prs_1d0,inds_p0/])

  ; Difference between 'lat_max'/'lon_max' and values read in using NCL function
    dlat          = abs(lat_max - centre(cn,it,0))
    dlon          = abs(lon_max - centre(cn,it,1))

  ; If our calculations of latitude and longitude deviate too much from the values
  ; we read in using the NCL function ('st_rm'), then calculate a new centre
  ; using the minimum pressure. Constrain the new calculation using the values
  ; from 'st_rm', to prevent spurious centres being identified and analysed.

    if (dlat .gt. 1.00 .or. dlon .gt. 1.00) then

   ; Should this be 'it-1' because pressure is in the 'pd' stream files?
     lt1 = centre(cn,it,0) - 0.25
     lt2 = centre(cn,it,0) + 0.25
     ln1 = centre(cn,it,1) - 0.25
     ln2 = centre(cn,it,1) + 0.25

     prs_plane0 = prs_plane({lt1:lt2},{ln1:ln2})       ; Smaller grid (1.0 x 1.0 deg)
     prs_min0   = min(prs_plane0)                      ; Find minimum on smaller grid
     dims_p0    = dimsizes(prs_plane0)                 ; Size of smaller grid
     prs_1d0    = ndtooned(prs_plane0)                 ; Create 1-D array
     inds_p0    = ind_resolve(minind(prs_1d0),dims_p0) ; Find index of minimum

     lat_max    = lt1 + (dy0 * inds_p0(0,0))
     lon_max    = ln1 + (dx0 * inds_p0(0,1))
     print("New centre: "+lat_max+" degrees N, "+lon_max+" degrees E")

    end if

  elseif (calc .eq. 1) then 

    lat_max    = centre(cn,it,0)
    lon_max    = centre(cn,it,1)

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

   ;===============
   ; 'pk' stream
   ;===============

   ; Indices of lat/lon grid for 'str_xz' in 'pd' stream
     latK_plane  = ind(t0.le.latK.and.latK.le.t1)
     lonK_plane  = ind(n0.le.lonK.and.lonK.le.n1)

   ; Array sizes
     latK_size   = dimsizes(latK_plane)
     lonK_size   = dimsizes(lonK_plane)

   ; Array of latitude/longitude points corresponding to the indices above
     latK_sub    = latK(latK_plane(0):latK_plane(latK_size-1))
     lonK_sub    = lonK(lonK_plane(0):lonK_plane(lonK_size-1))

   ; Indices of lat/lon grid for 'str_yz' in 'pd' stream
     lat1K_plane = ind(t0.le.lat1K.and.lat1K.le.t1)
     lon1K_plane = ind(n0.le.lon1K.and.lon1K.le.n1)

   ; Array sizes
     lat1K_size  = dimsizes(lat1K_plane)
     lon1K_size  = dimsizes(lon1K_plane)

   ; Array of latitude/longitude points corresponding to the indices above
     lat1K_sub   = lat1K(lat1K_plane(0):lat1K_plane(lat1K_size-1))
     lon1K_sub   = lon1K(lon1K_plane(0):lon1K_plane(lon1K_size-1))

     delete([/vort1d,inds,vort1d_sm,vort_smth,inds_sm/])

;===========================================
; Switch to cylindrical coordinate system
;===========================================

   ; Create new arrays describing cylindrical coordinate system (2D for now)
     xpos = new((/dimsizes(radii),dimsizes(thetas)/),float)
     ypos = new((/dimsizes(radii),dimsizes(thetas)/),float)

   ; New arrays containing variables we want to plot
     u_int   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Zonal wind
     v_int   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Meridional wind
     spd_int = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vector wind (for RMW)
     pv_int  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Potential vorticity
     w_int   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vertical velocity
     v_rad   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Radial wind
     v_tan   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Tangential wind
     prs_int = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Pressure 

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

     copy_VarAtts(prs_plane,prs_int)
     prs_int@description = "Pressure"
     prs_int@units       = "hPa"

     v_rad@description = "Radial wind"
     v_rad@units       = "m s~S~-1~N~"

     v_tan@description = "Tangential wind"
     v_tan@units       = "m s~S~-1~N~"

     xpos!0    	 = "rad"
     xpos!1    	 = "azi"
     xpos&rad  	 = radii
     xpos&azi  	 = thetas

     ypos!0    	 = "rad"
     ypos!1    	 = "azi"
     ypos&rad  	 = radii
     ypos&azi  	 = thetas

     u_int!0   	 = "rad"
     u_int!1   	 = "azi"
     u_int&rad 	 = radii
     u_int&azi 	 = thetas

     v_int!0   	 = "rad"
     v_int!1     = "azi"
     v_int&rad 	 = radii
     v_int&azi 	 = thetas

     spd_int!0   = "rad"
     spd_int!1   = "azi"
     spd_int&rad = radii
     spd_int&azi = thetas

     w_int!0   	 = "rad"
     w_int!1   	 = "azi"
     w_int&rad 	 = radii
     w_int&azi 	 = thetas

     pv_int!0    = "rad"
     pv_int!1    = "azi"
     pv_int&rad  = radii
     pv_int&azi  = thetas

     prs_int!0   = "rad"
     prs_int!1   = "azi"
     prs_int&rad = radii
     prs_int&azi = thetas

     v_tan!0   	 = "rad"
     v_tan!1   	 = "azi"
     v_tan&rad 	 = radii
     v_tan&azi 	 = thetas

     v_rad!0   	 = "rad"
     v_rad!1   	 = "azi"
     v_rad&rad 	 = radii
     v_rad&azi 	 = thetas

   ; Latitude/longitude (small grid)
     do irad = 0, dimsizes(radii)-1
       r = radii(irad)
       do iang = 0, dimsizes(thetas)-1
         theta = thetas_rad(iang)
           xpos(irad,iang) = lon_max + r*cos(theta)
           ypos(irad,iang) = lat_max + r*sin(theta)
       end do
     end do

   ; Fill new arrays with interpolated values on the cylindrical grid
   ; Which latitude/longitude grids are each of these variables defined on?
   ; u_plane     -- latC,  lonC
   ; v_plane     -- latC,  lonC
   ; prs_plane	 -- latD,  lonD
   ; w_plane     -- latK,  lon1K
   ; pv_plane	 -- latK,  lon1K

  ; Loop over azimuth angles to fill entire circular array of data points
    do iang = 0, dimsizes(thetas_rad)-1
      u_int(:,iang)    = linint2_points(lonC_sub,latC_sub,u_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      v_int(:,iang)    = linint2_points(lonC_sub,latC_sub,v_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      prs_int(:,iang)  = linint2_points(lonD_sub,latD_sub,prs_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      w_int(:,iang)    = linint2_points(lon1K_sub,latK_sub,w_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      pv_int(:,iang)   = linint2_points(lon1K_sub,latK_sub,pv_plane,False,\
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
      vvel_azi(irad) = avg(w_int(irad,:))	 ; Vertical velocity
    end do

  ; Input azimuthal averages into larger arrays before plotting  
    vtan_plot0(:,ct)  = vtan_azi(:)              ; Tangential wind
    vrad_plot0(:,ct)  = vrad_azi(:)	    	 ; Radial wind
    pvrt_plot0(:,ct)  = pvrt_azi(:)	    	 ; PV
    vvel_plot0(:,ct)  = vvel_azi(:)	    	 ; Vertical velocity

  ; Tidy up
    delete([/u_int,v_int,w_int,pv_int,prs_int,v_tan,v_rad/])
    delete([/u_plane,v_plane,w_plane,vort_plane,prs_plane,pv_plane/])
    delete([/str_plane, str_plane0, str1, str2/])

    ct  = ct + 1		    ; Counter variable (time)

  ; Tidy up (avoid dimension size errors when iterating over times)
    delete([/du_dx,dv_dy,dv_dx,du_dy/])
    delete([/latC_sub,lonC_sub,latC_plane,lonC_plane,lonC,latC,u,v,vort,avo,f0,f1/])
    delete([/latD_sub,lonD_sub,latD_plane,lonD_plane,lonD,latD,p,z/])
    delete([/latK_sub,lonK_sub,lon1K_sub,lat1K_sub,latK_plane,lat1K_plane,lonK_plane,lon1K_plane/])
    delete([/lon1K,lat1K,lonK,latK,hybK,pv,w/])

   end do     ; End time loop (do it = ts, tf)

;=====================================================
; Calculate radial and tangential wind tendencies
;=====================================================

; Initialise new arrays
  vtan_tend0 	  = vtan_plot0
  vrad_tend0 	  = vrad_plot0

  vtan_tend0!0    = "rad"
  vtan_tend0&rad  = radii
  vtan_tend0!1    = "time"
  vtan_tend0&time = time_int

  vrad_tend0!0    = "rad"
  vrad_tend0&rad  = radii
  vrad_tend0!1    = "time"
  vrad_tend0&time = time_int

  print("Calculating radial and tangential wind tendencies...")

; EDIT -- check that this calculation is okay 
  do ct = 0, nTIMES-2

   tp1	= min((/ct+1,nTIMES/))
   tm1	= max((/ct-1,0/))

   vtan_tend0(:,ct) = (vtan_plot0(:,tp1) - vtan_plot0(:,tm1)) * 6 ; Units [m s-1 h-1]
   vrad_tend0(:,ct) = (vrad_plot0(:,tp1) - vrad_plot0(:,tm1)) * 6 ; Units [m s-1 h-1]

  end do     ; End time loop (do ct = nTIMES-1)

; Also calculate tendency for final element in array 
  vtan_tend0(:,nTIMES-1) = (vtan_plot0(:,nTIMES-1) - vtan_plot0(:,nTIMES-2) ) * 12 
  vrad_tend0(:,nTIMES-1) = (vrad_plot0(:,nTIMES-1) - vrad_plot0(:,nTIMES-2) ) * 12

;======================================
; Reorder dimensions before plotting 
;======================================

 ; Apply running average, if required (6 points in running ave; calculated on dimension 1)
   if (ave .eq. 1) then
    vtan_tend1 = runave_n_Wrap(vtan_tend0, 6, 0, 1)
    vrad_tend1 = runave_n_Wrap(vrad_tend0, 6, 0, 1)
    vvel_plot1 = runave_n_Wrap(vvel_plot0, 6, 0, 1)
    pvrt_plot1 = runave_n_Wrap(pvrt_plot0, 6, 0, 1)
   end if

 ; Arrays now ordered [TIME * RADIUS] rather than [RADIUS * TIME]
   vtan_tend = vtan_tend1(time|:,rad|:)
   vrad_tend = vrad_tend1(time|:,rad|:)

   vtan_plot = vtan_plot0(time|:,rad|:)
   vrad_plot = vrad_plot0(time|:,rad|:)
   vvel_plot = vvel_plot1(time|:,rad|:)
   pvrt_plot = pvrt_plot1(time|:,rad|:)

 ; Add additional metadata before plotting 
   vvel_plot!0    = "time"
   vvel_plot&time = time_int

;=====================================================
; Create a panel plot of Figs. 6a and 6b from N11 
;=====================================================

  ; Output file location and type
    if (lay .eq. 1) then 
     output = "$sam/nepartak/images/cyl_coords/n11_fig6_5min_"\
              +dat+"_"+ens0+"_"+ts+"_"+tf+"_"+hy0+"m_"+hy1"m"
    elseif (lay .eq. 0) then 
     output = "$sam/nepartak/images/cyl_coords/n11_fig6_5min_"\
              +dat+"_"+ens0+"_"+ts+"_"+tf+"_"+hy0"m"
    end if 
    wks = gsn_open_wks(opt,output)

  ; Set up panel plot
    panel = new(3,graphic)
    gsn_define_colormap(wks,"n11")

;==========================
; Options for plotting 
;==========================

  ; Tangential velocity
    opts_vtan                               = True
    opts_vtan@cnFillOn                      = True
    opts_vtan@cnLineLabelInterval           = 2.0
    opts_vtan@cnLineLabelFontHeightF        = 0.012
    opts_vtan@cnLineLabelBackgroundColor    = "transparent"
    opts_vtan@cnLineLabelPlacementMode      = "constant"
    opts_vtan@cnLinesOn                     = False ; Contour lines off
    opts_vtan@cnInfoLabelOn                 = False
    opts_vtan@cnLevelSelectionMode          = "ExplicitLevels"

    if (tclr .eq. 0) then 
     opts_vtan@cnLevels                     = (/1.0, 2.0, 2.5, 3.0, 3.5, 4.0, \
                                                4.5, 5.0, 6.0, 7.0/)
     opts_vtan@cnFillColors                 = (/0,2,3,4,5,6,7,8,9,10,11/)
    elseif (tclr .eq. 1) then 
     opts_vtan@cnLevels                     = (/0.5, 1.5, 2.0, 2.5, 3.0, 3.5, \
                                                4.0, 4.5, 5.0, 6.0/)
     opts_vtan@cnFillColors                 = (/0,2,3,4,5,6,7,8,9,10,11/)
    elseif (tclr .eq. 2) then		    
     opts_vtan@cnLevels                     = (/0.5, 1.0, 1.5, 2.0, 2.5, 3.0, \
                                                3.5, 4.0, 4.5, 5.0/)
     opts_vtan@cnFillColors                 = (/0,2,3,4,5,6,7,8,9,10,11/)
    end if 

    opts_vtan@gsnPaperOrientation           = "landscape"
    opts_vtan@tiMainString                  = ""
    opts_vtan@tiMainFontHeightF             = 0.0125
    opts_vtan@gsnLeftString                 = ""
    opts_vtan@gsnRightString                = ""
    opts_vtan@gsnMaximize                   = True
    opts_vtan@lbLabelBarOn                  = True                  ; Labelbar on/off
    opts_vtan@lbBoxEndCapStyle              = "TriangleBothEnds"    ; Labelbar end shape
    opts_vtan@pmLabelBarWidthF              = 0.32                  ; Labelbar width
    opts_vtan@pmLabelBarHeightF             = 0.08                  ; Labelbar height
    opts_vtan@lbLabelFontHeightF            = 0.0125                ; Labelbar font size
    opts_vtan@lbLabelFont                   = "Helvetica"           ; Labelbar font
    opts_vtan@lbPerimOn                     = False                 ; Perimeter on/off
    opts_vtan@gsnDraw                       = False                 ; Do not draw plot
    opts_vtan@gsnFrame                      = False                 ; Do not advance frame
    opts_vtan@vpWidthF                      = 0.30                  ; Match width in N11
    opts_vtan@vpHeightF                     = 0.50                  ; Match height in N11

  ; Additional plotting resources
    opts_vtan@tiYAxisString                 = "Time (h)"
    opts_vtan@tiXAxisString                 = "Radius (km)"
    opts_vtan@tiXAxisFontHeightF            = 0.018
    opts_vtan@trYLog                        = False
    opts_vtan@tmXTOn                        = "False" 		    ; Turn off top x-axis TM
    opts_vtan@tmYROn                        = "False" 		    ; Turn off right y-axis TM
    opts_vtan@gsnMaximize                   = True    		    ; Maximise plot size
    opts_vtan@gsnAddCyclic                  = False		    ; No cyclic point 
    opts_vtan@tmXBMode                      = "Explicit"                ; Choose tm mode
    opts_vtan@tmXBValues                    = (/0.0,0.5,1.0,1.5/)       ; Tm positions
    opts_vtan@tmXBLabels                    = (/"0","50","100","150"/)  ; Tm values
    opts_vtan@tmXBMinorOn                   = False

    opts_vtan@pmLabelBarOrthogonalPosF      = 0.10	       ; Move lb up/down (higher/lower)
    opts_vtan@tiXAxisOffsetYF               = 0.10    	       ; Move x-axis title up/down (h/l)

    opts_vtan@trXMinF                       = 0.0              ; Set x-axis min
    opts_vtan@trXMaxF                       = 1.0              ; Set x-axis max

  ; Y-axis labels
    if (expl .eq. 1) then
     opts_vtan@tmYLMode                     = "Explicit"
     opts_vtan@tmYLValues                   = (/12,24,36,48,60,72,84,96,108,120,132,144,156/)

     if(ts .eq. 563) then
      opts_vtan@tmYLLabels                   = (/"T+50","T+51","T+52","T+53","T+54","T+55","T+56",\
                                                 "T+57","T+58","T+59","T+60","T+61","T+62"/)
     elseif (ts .eq. 803) then
      opts_vtan@tmYLLabels                   = (/"T+69","T+70","T+71","T+72","T+73","T+74","T+75",\
                                                 "T+76","T+77","T+78","T+79","T+80","T+81"/)
     end if

    else
     opts_vtan@tmYLMode                     = "Manual"
     opts_vtan@tmYLTickStartF               = opts_vtan@trYMinF         ; See L730
     opts_vtan@tmYLTickEndF                 = opts_vtan@trYMaxF         ; See L731
     opts_vtan@tmYLTickSpacingF             = 6
    end if

  ; Y-axis limits (time)
    opts_vtan@trYMinF                       = ts + 12
    opts_vtan@trYMaxF                       = tf - 12
    opts_vtan@tmYLMinorOn                   = False

  ; Resources to overlay contours of mean tangential wind
    opts_cont				    = True
    opts_cont@cnFillOn                      = False       ; Filled contours on/off
    opts_cont@cnLineColor                   = "blue3"     ; Contour line colour
    opts_cont@cnLevelSelectionMode          = "ExplicitLevels"
    opts_cont@cnLevels                      = (/10,20,30,40,50,60,70/)
    opts_cont@cnLineThicknessF              = 2.0         ; Contour line thickness
    opts_cont@cnInfoLabelOn                 = False       ; Info label on/off
    opts_cont@gsnDraw                       = False       ; Do not draw the plot
    opts_cont@gsnFrame                      = False       ; Do not advance the frame
    opts_cont@cnLineLabelsOn                = True        ; Turn line labels on/off
    opts_cont@cnLineLabelFontColor          = "blue3"     ; Line label colour
    opts_cont@cnLineLabelFontHeightF        = 0.01        ; Default 0.006 for vpWidth=0.30
    opts_cont@cnLineLabelFontThicknessF     = 2.0         ; Line label thickness
    opts_cont@tiMainString                  = ""
    opts_cont@gsnLeftString                 = ""
    opts_cont@gsnRightString                = ""

  ; Resources to overlay zero line
    opts_zero                               = True
    opts_zero@cnFillOn                      = False
    opts_zero@cnLineColor                   = "gray82"
    opts_zero@cnLevels                      = (/-1.0, 0.0, 1.0/)
    opts_zero@cnLineThicknessF              = 2.0
    opts_zero@cnInfoLabelOn                 = False
    opts_zero@gsnDraw                       = False       ; Do not draw the plot
    opts_zero@gsnFrame                      = False       ; Do no advance the frame
    opts_zero@gsnContourZeroLineThicknessF  = 2.0
    opts_zero@cnLineLabelsOn                = False      ; Turn off line labels
    opts_zero@tiMainString                  = ""
    opts_zero@gsnLeftString                 = ""
    opts_zero@gsnRightString                = ""

  ; Radial velocity resources
    opts_vrad                               = opts_vtan
    delete(opts_vrad@cnLevels)
    delete(opts_vrad@cnFillColors)

    if (rclr .eq. 0) then 
     opts_vrad@cnLevels                     = (/0.0, -0.5, -1.0, -1.5, -2.0, \
                                                -2.5, -3.0, -4.0, -5.0/)
     opts_vrad@cnFillColors                 = (/12,11,10,9,8,7,6,5,4,0/)
    elseif (rclr .eq. 1) then 
     opts_vrad@cnLevels                     = (/0.0, -0.25, -0.5, -1.0, -1.5, \
                                                -2.0, -2.5, -3.0, -4.0/)
     opts_vrad@cnFillColors                 = (/12,11,10,9,8,7,6,5,4,0/)
    elseif (rclr .eq. 2) then 
     opts_vrad@cnLevels                     = (/0.0, -0.25, -0.5, -0.75, -1.0, \
                                                -1.5, -2.0, -2.5, -3.0/)
     opts_vrad@cnFillColors                 = (/12,11,10,9,8,7,6,5,4,0/)
    end if

  ; Vertical velocity resources
    opts_vvel                                = opts_vtan
    delete(opts_vvel@cnLevels)
    delete(opts_vvel@cnFillColors)

    if (wclr .eq. 0) then
     opts_vvel@cnLevels                      = (/0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, \
                                                 0.7, 0.8, 0.9, 1.0/)
     opts_vvel@cnFillColors                  = (/0,2,3,4,5,6,7,8,9,10,11,12/)
    elseif (wclr .eq. 1) then
     opts_vvel@cnLevels                      = (/0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, \
                                                 0.50, 0.60, 0.70, 0.80/)
     opts_vvel@cnFillColors                  = (/0,2,3,4,5,6,7,8,9,10,11,12/)    

    elseif (wclr .eq. 2) then
     opts_vvel@cnLevels                      = (/0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, \
                                                 0.40, 0.45, 0.50, 0.60/)
     opts_vvel@cnFillColors                  = (/0,2,3,4,5,6,7,8,9,10,11,12/)
    end if 

  ; Radius of maximum wind (RMW) resources
    opts_rmw                               = True
    opts_rmw@xyLineColors                  = "red3"
    opts_rmw@xyLineThicknesses             = 3.0
    opts_rmw@xyDashPatterns                = 0.0
    opts_rmw@gsnDraw                       = False       ; Do not draw the plot
    opts_rmw@gsnFrame                      = False       ; Do not advance the frame
    opts_rmw@tiMainString                  = ""
    opts_rmw@gsnLeftString                 = ""
    opts_rmw@gsnRightString                = ""

;==============================================
; Read in RMW data and construct panel plot
;==============================================

  ; Read in RMW data from text file (first 3 values equal to zero)
    spd_plot0 = asciiread("$ar/text/"+dat+"_"+ens0+"_rmw.txt",(120),"float")
    spd_plot  = new(numTIMES,"float")

  ; Read in data every hour and leave the other elements missing
    spd_plot(0)  = spd_plot0(53)

  ; Initialised counter variable
    nd = 54

;  ; Read in chosen values from 1-h RMW text file and leave the rest missing
;    do i = 11, tf, 12
;     spd_plot(i) = spd_plot0(nd)
;     nd = nd + 1
;    end do

;    spd_plot@_FillValue = 9.96921e+36

   ; Tangential wind (LHS panel)
     vtan0 = gsn_csm_contour(wks,vtan_tend,opts_vtan)    ; Tangential wind tendency (filled)
     zero0 = gsn_csm_contour(wks,vtan_tend,opts_zero)    ; Zero line (vtan tendency)
     zero  = ColorNegDashZeroPosContour(zero0,"transparent","gray69","transparent")
     cont0 = gsn_csm_contour(wks,vtan_plot,opts_cont)    ; Tangential wind
;     rmw0  = gsn_csm_xy(wks,spd_plot,time_int,opts_rmw)  ; RMW
     overlay(vtan0,zero)                                 ; Overlay vtan tendency zero line
;     overlay(vtan0,rmw0)				 ; Overlay RMW
     overlay(vtan0,cont0)                                ; Overlay tangential wind
     panel(0) = vtan0                                    ; Plot first panel
     delete([/zero0,zero/])                              ; Tidy up 

   ; Radial wind (middle panel)
     vrad0 = gsn_csm_contour(wks,vrad_tend,opts_vrad)    ; Radial wind tendency
     cont0 = gsn_csm_contour(wks,vtan_plot,opts_cont)    ; Tangential wind
;     rmw0  = gsn_csm_xy(wks,spd_plot,time_int,opts_rmw)  ; RMW
     overlay(vrad0,cont0)                                ; Overlay tangential wind
;     overlay(vrad0,rmw0)                                 ; Overlay RMW
     panel(1) = vrad0                                    ; Plot second panel 

   ; Vertical velocity (RHS panel)
     vvel0 = gsn_csm_contour(wks,vvel_plot,opts_vvel)    ; Vertical velocity
     zero0 = gsn_csm_contour(wks,vvel_plot,opts_zero)    ; Zero line (vertical velocity)
     zero  = ColorNegDashZeroPosContour(zero0,"transparent","gray69","transparent")
;     rmw0  = gsn_csm_xy(wks,spd_plot,time_int,opts_rmw)  ; RMW
     overlay(vvel0,zero)                                 ; Overlay vvel zero line
;     overlay(vvel0,rmw0)                                 ; Overlay RMW
     panel(2) = vvel0                                    ; Plot third panel
     delete([/zero0,zero/])                              ; Tidy up

;=============================================
; Finally, put everything together and plot 
;=============================================

    optsP                       = True
    optsP@gsnFrame              = False    ; Do not advance the frame
    optsP@gsnPanelLabelBar      = False    ; Turn on panel labelbar
    optsP@txString              = "RA1T ensemble simulation ("+ens0+")"
    optsP@gsnPanelFigureStrings = (/"a) vtan_tend","b) vrad_tend","c) vvel"/)
    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.90
    optsP@gsnPanelBottom        = 0.10
    optsP@amJust                = "TopLeft"
    optsP@gsnPanelFigureStringsFontHeightF = 0.0099 ; Reduce label size (default 0.01)

    optsP@gsnPanelRowSpec       = True 
    gsn_panel(wks,panel,(/3/),optsP)
    frame(wks)

end 