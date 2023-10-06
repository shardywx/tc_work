; Script to produce a Hovmoller plot of the 3D wind field, following Nguyen et al. (2011)

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_centre.ncl"
load "$sam/ncl_func/setup_cyl.ncl"
load "$sam/ncl_func/ring_mono.ncl"

; ncl opt=\"pdf\" dat=\"03T00\" ens0=\"em06\" dist=3.0 ts=60 tf=116
; rmw=1 calc=1 tclr=2 pclr=2 wclr=1 basic=0 lay=1 ar=1.5 nr=31 nt=73
; ar0=1.0 cn0=\"geo_sm\" sc=5 sm=0 vlev0=8 vlev1=14 wlev0=17 wlev1=28 tc_hov_mslp.ncl

; 'dat'    = initialisation time string ("02T12", "03T00", ...)
; 'ens0'   = ensemble simulation (em00-em11)
; 'dist'   = size of box following storm (degrees)
; 'opt'    = output file format ("pdf" or "x11")
; 'mlev'   = model level for plotting PV and vertical velocity (see L106-108)
; 'rmw'	   = read in RMW data on model (1) or pressure (0) levels
; 'calc'   = centre calculated on each level (0), calculated offline (1) or using 980 m vort (2)
; 'tclr'   = colour scheme for plotting tangential wind (0-2)
; 'pclr'   = colour scheme for plotting radial wind (0-2)
; 'wclr'   = colour for plotting vertical velocity (0-2)
; 'basic'  = plot u,v,w (1) or tendencies following N11 (0)
; 'lay'    = calculate layer average (1) or single level (0)
; 'mlev2'  = upper model level for calculating layer average (only valid when lay=1)
; 'ar'	   = distance in degrees of outer radius (1.0, 2.0, etc)
; 'nr'	   = number of radial circles between r = 0 and r = 'ar' (21, 41, 61, etc)
; 'ar0'	   = outer radius (degrees) for final plot [must be <= 'ar']
; 'cn0'    = method of storm tracking ["slp", "geo", "geo_sm", "vort"] - default is "geo_sm"

; 'mlev'   = 13 and 'mlev2' = 24 give an approximate 1-3 km layer average

; 'pc' stream --> 'it'
; 'centre'    --> 'it-1'
; 'pd' stream --> 'it-2'
 
begin

;==============================
; Prelash (reading in files)
;==============================

; Choose correct data stream(s) to read in and plot
  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0
  fili_a1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pa*.nc")
  fili_a      = diri+"/"+fili_a1
  fili_c1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pc*.nc")
  fili_c      = diri+"/"+fili_c1
  fili_d1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pd*.nc")
  fili_d      = diri+"/"+fili_d1

; How many input times in 'pa' stream? 
  num_fili_a  = dimsizes(fili_a)
  num_a       = num_fili_a * 12 

; Calculate number of times in file using built-in function ('nc_times')
  numINFO_c   = nc_times(fili_c)
  numINFO_d   = nc_times(fili_d)

; Number (numINPUT_x) and array (time_x) of times
  numINPUT_c  = numINFO_c[0]
  time_c      = numINFO_c[1]
  numINPUT_d  = numINFO_d[0]
  time_d      = numINFO_d[1]

; Calculate storm motion using built-in function ('st_rm')
; 'dat'     = initialisation time (02T12, 03T00, etc)
; 'ens0'    = ensemble simulation (em00-em11)
; 'diri'    = path to input files (see above)
; 'dist'    = size of box following storm (degrees)
; 'mins'    = analyse 1-h (0) or 5-min (1) data

  storm_rel = st_rm(dat, ens0, "$ar/text/", dist, 0)
  u_cyc     = storm_rel[0]             ; Zonal wind
  v_cyc     = storm_rel[1]             ; Meridional wind
  vel_cyc   = storm_rel[2]             ; Vector wind
  numTIMES  = storm_rel[3]             ; Number of times in file (120)
  lat_arr   = storm_rel[4]             ; Latitude array subset (following storm)
  lon_arr   = storm_rel[5]             ; Longitude array subset (following storm)
  centre    = storm_rel[6]             ; Storm track (position) array

; Create arrays for data and output information
  llbox     = toint(dist*50)                   ; Calculate domain size (grid points)
  dsize     = (/llbox,llbox/)                  ; Domain size (grid points)

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 10000000000
  end setvalues

;=======================================                                                
; Choose which storm track to read in                                                   
;=======================================                                                

  if (cn0 .eq. "slp") then
   cn = 0
  elseif (cn0 .eq. "vort") then
   cn = 1
  elseif (cn0 .eq. "geo") then
   cn = 2
  elseif (cn0 .eq. "geo_sm") then
   cn = 3
  end if

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

;========================================
; Define cylindrical coordinate arrays
;========================================

; ar = distance (degrees) of outer radial circle (1.0)
; nr = total number of radial circles (21)

  pii        = 3.14159265
  radii      = fspan(0.0,ar,nr)      ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,nt)   ; Array of azimuth angles (0, 5, 10, etc)
  thetas_rad = (pii/180)*thetas      ; Degrees to radians
  rad_size   = dimsizes(radii)-1     ; Index for later calculation

; Define arrays for later use in tangential velocity cross-section plots
  vtan_azi   = new((/dimsizes(radii)/),float) ; Tangential wind
  vrad_azi   = new((/dimsizes(radii)/),float) ; Radial wind
  spd_azi    = new((/dimsizes(radii)/),float) ; Vector wind
  pvrt_azi   = new((/dimsizes(radii)/),float) ; Potential vorticity
  vvel_azi   = new((/dimsizes(radii)/),float) ; Vertical velocity
  mslp_azi   = new((/dimsizes(radii)/),float) ; Mean sea level pressure 
  avel_azi   = new((/dimsizes(radii)/),float) ; Angular velocity 

  nTIMES     = (tf - ts) + 1		      ; Number of times in user-defined subset

; Define 2-dimensional array to hold all values from 'vtan_azi' at ALL times
  vtan_plot0   = new((/dimsizes(radii),numINPUT_c/),float)
  vtan_plot0!0 = "rad"
  vtan_plot0!1 = "time"

  vrad_plot0   = new((/dimsizes(radii),numINPUT_c/),float)
  vrad_plot0!0 = "rad"
  vrad_plot0!1 = "time"

  pvrt_plot0   = new((/dimsizes(radii),numINPUT_c/),float)
  pvrt_plot0!0 = "rad"
  pvrt_plot0!1 = "time"

  vvel_plot0   = new((/dimsizes(radii),numINPUT_c/),float)
  vvel_plot0!0 = "rad"
  vvel_plot0!1 = "time"

  mslp_plot0   = new((/dimsizes(radii),numINPUT_c/),float)
  mslp_plot0!0 = "rad"
  mslp_plot0!1 = "time"

  avel_plot0   = new((/dimsizes(radii),numINPUT_c/),float)
  avel_plot0!0 = "rad"
  avel_plot0!1 = "time"

; Add metadata for azimuthally averaged variables
  vtan_azi!0            = "rad"
  vtan_azi&rad          = radii
  vtan_azi@description  = "Azimuthally averaged tangential velocity"
  vtan_azi@units        = "m s~S~-1~N~"

  vrad_azi!0            = "rad"
  vrad_azi&rad          = radii
  vrad_azi@description  = "Azimuthally averaged radial velocity"
  vrad_azi@units        = "m s~S~-1~N~"

  spd_azi!0             = "rad"
  spd_azi&rad           = radii
  spd_azi@description   = "Azimuthally averaged vector velocity"
  spd_azi@units         = "m s~S~-1~N~"

  pvrt_azi!0            = "rad"
  pvrt_azi&rad          = radii
  pvrt_azi@description  = "Azimuthally averaged PV"
  pvrt_azi@units        = "10~S~-6~N~K kg~S~-1~N~ m~S~2~N~ s~S~-1~N~"

  vvel_azi!0            = "rad"
  vvel_azi&rad          = radii
  vvel_azi@description  = "Azimuthally averaged vertical velocity"
  vvel_azi@units        = "m s~S~-1~N~"

  mslp_azi!0            = "rad"
  mslp_azi&rad          = radii
  mslp_azi@description  = "Azimuthally averaged vertical velocity"
  mslp_azi@units        = "m s~S~-1~N~"

  avel_azi!0            = "rad"
  avel_azi&rad          = radii
  avel_azi@description  = "Azimuthally averaged vertical velocity"
  avel_azi@units        = "m s~S~-1~N~"

; Create time array for later plotting 
  time_int   = new(numINPUT_c,"integer")
  time_int   = ispan(0,120,1)
  time_int!0 = "time"

; Get MSLP array dimensions
  a1         = addfile(fili_a(0), "r")
  mslp1      = a1->p
  slp_size   = dimsizes(mslp1)

; Now create array to hold MSLP values from 'pa' stream loop
  mslp_arr   = new((/numTIMES,slp_size(2),slp_size(3)/),float)
  mslp_arr!0 = "time"
  mslp_arr!1 = "lat"
  mslp_arr!2 = "lon"

;===============================
; Start loop over input files
;===============================

; Create array to hold storm track positions using different methods
  num_l       = 59
  centre_new  = new((/7, num_l, numINPUT_c, 2/),float)

; Arrays to hold all 'xpos' and 'ypos' position arrays
  xpos_all   = new((/num_l,dimsizes(radii),dimsizes(thetas),numINPUT_c/),float)
  ypos_all   = new((/num_l,dimsizes(radii),dimsizes(thetas),numINPUT_c/),float)

; Arrays to hold lat/lon at centre of cylindrical grid
  xcen_all   = new((/num_l,dimsizes(radii),dimsizes(thetas),numINPUT_c/),float)
  ycen_all   = new((/num_l,dimsizes(radii),dimsizes(thetas),numINPUT_c/),float)

;=============================================
; Read in 'pa' stream data before main loop
;=============================================

  do nf      = 0, num_fili_a-1
   a         = addfile(fili_a(nf), "r")

;   print("Working on file number: "+nf)

 ; Times in file  
;   if (vc .eq. 1) then
;    time  = a->t_1
;    times = 12
;   elseif (vc .eq. 0) then
    time  = a->t
    times = dimsizes(time)
;   end if

 ; Time indices 
   ts0 = nf * 12 
   tf0 = ts0 + 11

;   print("start index = "+ts+" ; end index = "+tf)

 ; Read in mean sea level pressure (Pa) [t | z | lat | lon] 
 ; For 1-h data (vc = 0), ignore first time interval (nf = 0)
   if (nf .eq. 0) then ; .and. vc .eq. 0) then
    mslp  = a->p(1:times-1,:,:,:)
   else
    mslp  = a->p
   end if 

 ; Convert to hPa
   mslp  = mslp / 100 
   mslp@units = "hPa"

 ; Read into larger array
   mslp_arr(ts0:tf0,:,:) = (/mslp(:,0,:,:)/) 

 ; Tidy up
   delete([/mslp,time,times/])

  end do 

;==================================
; Get the variables we will need
;==================================

  c  = addfile(fili_c,"r")             ; Read in 'pc' stream [u,v,vort] 
  d  = addfile(fili_d,"r")             ; Read in 'pd' stream [w,p,t,theta,kmh,tau] 

  ct = 0                               ; Counter variable (time) 

  do it = ts, tf                   ; Loop over times in file

    print("Working on time: "+time_str(it))
    itt = it
    title_arr(it) = "Valid at "+time_str(it)+" (T+"+itt+")"

;=====================================================
; Read in variables from 'pc' stream (model levels)
;=====================================================

    t0    = lat_arr(cn,it-1,0)
    t1    = lat_arr(cn,it-1,1)
    n0    = lon_arr(cn,it-1,0)
    n1    = lon_arr(cn,it-1,1)
    print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

  ; Combine values above into array for input into external function
    ll_arr = (/t0,t1,n0,n1/)

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
    u0  = c->$u_varname$(it,vlev0:vlev1,{t0:t1},{n0:n1})     ; Zonal wind (layer)
    v0  = c->$v_varname$(it,vlev0:vlev1,{t0:t1},{n0:n1})     ; Meridional wind (layer)
   else
    u0  = c->$u_varname$(it,vlev0,{t0:t1},{n0:n1})           ; Zonal wind (single level)
    v0  = c->$v_varname$(it,vlev0,{t0:t1},{n0:n1})           ; Meridional wind (single level)
   end if

 ; Grid spacing 
   dy   = latC(1) - latC(0)
   dx   = lonC(1) - lonC(0)

;===================================
; Calculate storm-relative winds
;===================================

    u = u0 - u_cyc(cn,it-1)
    v = v0 - v_cyc(cn,it-1)

  ; Add metadata from original horizontal wind arrays
    copy_VarCoords(u0,u)
    copy_VarAtts(u0,u)

    copy_VarCoords(v0,v)
    copy_VarAtts(v0,v)

  ; Calculate relative vorticity using centered finite differences
    vort = uv2vr_cfd(u, v, latC, lonC, 2)
    copy_VarCoords(u, vort)
    vort@units = "s~S~-1~N~"
    vort@name  = "Relative vorticity on model levels"

;  ; EDIT -- convert to standard plotting units?
;    vort       = vort * (10 ^ 4)
;    vort@units = "10~S~-4~N~ s~S~-1~N~"

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

;================================================
; Create subset of 'pa' stream variable (MSLP)
;================================================

  ; Add latitude and longitude metadata to 'pa' stream arrays 
    ll = c->longitude_1
    lt = c->latitude_1
    mslp_arr&lat = lt
    mslp_arr&lon = ll
    delete([/ll,lt/])

  ; Now select grid subset to match other data streams 
    mslp = mslp_arr(it-1,{t0:t1},{n0:n1})
    mslp@units = "hPa"

;=====================================================
; Read in variables from 'pd' stream (model levels)
;=====================================================

  ; Read in basic variables from 'pd' stream
    lon1D  = d->longitude_1({n0:n1})                   ; longitude [xxx grid points]
    lat1D  = d->latitude_1({t0:t1})                    ; latitude  [xxx grid points]
    lonD   = d->longitude({n0:n1})                     ; longitude [xxx grid points]
    latD   = d->latitude({t0:t1})                      ; latitude  [xxx grid points]
    hybD   = d->hybrid_ht(:)                           ; 63 model ('rho') levels

  ; Model levels in km 
    hy     = hybD / 1000
    hy0    = sprintf("%0.0f",hybD(wlev0))

    if (lay .eq. 1) then
     hy1    = sprintf("%0.0f",hybD(wlev1))
    end if

  ; Read in temperature, pressure, theta, vertical velocity and geopotential height
  ; All variables below -- [hybD | 63] * [latD | 150] * [lon1D | 151]

    if (lay .eq. 1) then 
     t      = d->temp(it-2,wlev0:wlev1,{t0:t1},{n0:n1})        ; Temperature (K)
     p      = d->p(it-2,wlev0:wlev1,{t0:t1},{n0:n1})           ; Pressure (Pa)
     th     = d->theta(it-2,wlev0:wlev1,{t0:t1},{n0:n1})       ; Potential temperature (K)
     w      = d->dz_dt(it-2,wlev0:wlev1,{t0:t1},{n0:n1})       ; Vertical velocity (m s-1)
     z      = d->ht(it-2,wlev0:wlev1,{t0:t1},{n0:n1})          ; Geopotential height (m)
     pv	   = d->field83(it-2,wlev0:wlev1,{t0:t1},{n0:n1})      ; Potential vorticity (PVU)
    elseif (lay .eq. 0) then
     t      = d->temp(it-2,wlev0,{t0:t1},{n0:n1})	       ; Temperature (K)
     p      = d->p(it-2,wlev0,{t0:t1},{n0:n1})                 ; Pressure (Pa)
     th     = d->theta(it-2,wlev0,{t0:t1},{n0:n1})             ; Potential temperature (K)
     w      = d->dz_dt(it-2,wlev0,{t0:t1},{n0:n1})             ; Vertical velocity (m s-1)
     z      = d->ht(it-2,wlev0,{t0:t1},{n0:n1})                ; Geopotential height (m)
     pv     = d->field83(it-2,wlev0,{t0:t1},{n0:n1})           ; Potential vorticity (PVU)
    end if 

  ; Convert units of pressure and geopotential height before continuing
    p  = p / 100           ; Convert to 'hPa'
    z  = z / 10		   ; Convert to 'dam'
    pv = pv * (10 ^ 6)	   ; Convert to 'PVU' [10 ^ 6 K m2 s-1 kg-1]

  ; Add metadata
    p@units  = "hPa"
    z@units  = "dam"
    pv@units = "PVU"

    dx0      = lonD(1) - lonD(0) ; Grid spacing (longitude)
    dy0      = latD(1) - latD(0) ; Grid spacing (latitude)

;===========================================
; Analyse variables on chosen model level
;===========================================

    if (lay .eq. 1) then 

     print("Working on layer between " + hy0+ " and "+hy1+" m AGL")

   ; 'pc' stream 
     u_plane       = dim_avg_n_Wrap(u(:,:,:),0)
     v_plane       = dim_avg_n_Wrap(v(:,:,:),0)
     vort_plane    = dim_avg_n_Wrap(vort(:,:,:),0)
     f_plane       = dim_avg_n_Wrap(f0(:,:,:),0)

   ; 'pd' stream
     prs_plane     = dim_avg_n_Wrap(p(:,:,:),0)
     pv_plane	   = dim_avg_n_Wrap(pv(:,:,:),0)
     w_plane	   = dim_avg_n_Wrap(w(:,:,:),0)

    else

     print("Working on model level: "+hy0+ " m")

   ; 'pc' stream
     u_plane       = u(:,:)
     v_plane       = v(:,:)
     vort_plane    = vort(:,:)
     f_plane       = f0(:,:)

   ; 'pd' stream
     prs_plane     = p(:,:)
     pv_plane      = pv(:,:)
     w_plane       = w(:,:)

    end if 

  ; Mean sea level pressure (no layer-averaging)
    slp_plane      = mslp(:,:)

  ; Calculate vector windspeed and smooth
    spd_plane      = sqrt( (u_plane ^ 2) + (v_plane ^ 2) )
    copy_VarCoords(u_plane, spd_plane)
    vort_pl_smth   = smth9_Wrap(vort_plane, 0.5, 0.5, True)

  ; Combine all pc and pd stream variables into larger arrays 
  ; pc_plane needs to contain 7 variables; the final two arrays here are placeholders 
    pc_plane       = (/u_plane, v_plane, spd_plane, vort_plane, \
                       vort_plane, slp_plane, slp_plane, slp_plane/)
    pd_plane       = (/w_plane, prs_plane, pv_plane, pv_plane/)

;=============================================================
; Call external function to calculate storm centre position
;=============================================================

 ; Set radial distance (º) when looking for vort min in ringlike phase
   r0              = 0.10

   centre_arr      = st_centre(0, it, cn, r0, centre, vort_plane, prs_plane, spd_plane, \
                               centre_new, dy, dx, "single", it, "pd", "off", "1h")

 ; Updated storm track information contained in 'centre_new'   
   centre_new      = centre_arr[0]

 ; Changed 'it-2' to 'it-1' to reflect size of time dimension (120 not 119)
   print("Storm centre = "+centre_new(sc,0,it-1,0)+" ; "+centre_new(sc,0,it-1,1) )

;=====================================================================
; Call external function to switch to cylindrical coordinate system 
;=====================================================================

 ; Set 'ilev' option before function call below
   if (lay .eq. 1) then
    ilev = 0
   end if

 ; Call 'setup_cyl' and output pressure/vorticity on cylindrical grid 
 ; Last option = 0 --> basic variables output 
 ; Last option = 1 --> additional MBGT variables output 
   cyl_arr  = setup_cyl(it, dist, ll_arr, lonC, latC, lonD, latD, lon1D, lat1D, \
                        pc_plane, pd_plane, numTIMES, centre_new, \
                        sc, sm, ar, nr, nt, ilev, num_l, 0, "single", it, \
                        0, xpos_all, ypos_all, xcen_all, ycen_all)

   lat_max  = cyl_arr[4]
   lon_max  = cyl_arr[5]
   radii    = cyl_arr[6]
   rad_size = cyl_arr[7]
   thetas   = cyl_arr[19]

   u_int    = cyl_arr[8]
   v_int    = cyl_arr[9]
   spd_int  = cyl_arr[10]
   vort_int = cyl_arr[11]
   vr_int   = cyl_arr[13]
   vt_int   = cyl_arr[14]
   w_int    = cyl_arr[15]
   prs_int  = cyl_arr[16]
   pv_int   = cyl_arr[23]
   slp_int  = cyl_arr[25]

   xpos     = cyl_arr[17]
   ypos     = cyl_arr[18]

   print("Storm centre = "+lat_max+" ; "+lon_max)

;===================================================
; Average variables azimuthally around the storm 
;===================================================

   do irad = 0, dimsizes(radii)-1
     r = radii(irad)
     vtan_azi(irad) = avg(vt_int(irad,:))        ; Tangential wind
     vrad_azi(irad) = avg(vr_int(irad,:))	 ; Radial wind 
     spd_azi(irad)  = avg(spd_int(irad,:))	 ; Vector wind (for RMW)
     pvrt_azi(irad) = avg(pv_int(irad,:))	 ; Potential vorticity
     vvel_azi(irad) = avg(w_int(irad,:))	 ; Vertical velocity
     mslp_azi(irad) = avg(slp_int(irad,:))       ; Mean sea level pressure 
   end do

 ; Input azimuthal averages into larger arrays before plotting  
   vtan_plot0(:,it-1)  = vtan_azi(:)             ; Tangential wind
   vrad_plot0(:,it-1)  = vrad_azi(:)	    	 ; Radial wind
   pvrt_plot0(:,it-1)  = pvrt_azi(:)	    	 ; PV
   vvel_plot0(:,it-1)  = vvel_azi(:)	    	 ; Vertical velocity
   mslp_plot0(:,it-1)  = mslp_azi(:)             ; Mean sea level pressure 

 ; Tidy up
   delete([/u_int,v_int,w_int,pv_int,vt_int,vr_int,spd_int,slp_int/])
   delete([/u_plane,v_plane,w_plane,f_plane,vort_plane,prs_plane,spd_plane,pv_plane,slp_plane/])

   ct  = ct + 1		    ; Counter variable (time)

 ; Tidy up (avoid dimension size errors)
   delete([/lonC,latC,lon1C,lat1C,u,v,vort,mslp/])
   delete([/lonD,latD,lon1D,lat1D,t,p,w,z,pv/])

  end do     ; End time loop (do it = ts, tf)

;========================
; Calculate tendencies
;========================

; Initialise new arrays
  vtan_tend0 	  = vtan_plot0
  vrad_tend0 	  = vrad_plot0
  mslp_tend0      = mslp_plot0

  vtan_tend0!0    = "rad"
  vtan_tend0&rad  = radii
  vtan_tend0!1    = "time"
  vtan_tend0&time = time_int

  vrad_tend0!0    = "rad"
  vrad_tend0&rad  = radii
  vrad_tend0!1    = "time"
  vrad_tend0&time = time_int

  mslp_tend0!0    = "rad"
  mslp_tend0&rad  = radii
  mslp_tend0!1    = "time"
  mslp_tend0&time = time_int

  print("Calculating tendencies...")

  do it = ts, tf

   tp1  = min((/it+1,tf/))
   tm1  = max((/it-1,ts/))

   vtan_tend0(:,it) = (vtan_plot0(:,tp1) - vtan_plot0(:,tm1)) / 2 ; Units [m s-1 h-1]
   vrad_tend0(:,it) = (vrad_plot0(:,tp1) - vrad_plot0(:,tm1)) / 2 ; Units [m s-1 h-1]
   mslp_tend0(:,it) = (mslp_plot0(:,tp1) - mslp_plot0(:,tm1)) / 2 ; Units [hPa h-1]

  end do     ; End time loop (do it = ts, tf)

;======================================
; Reorder dimensions before plotting 
;======================================

 ; Arrays now ordered [TIME * RADIUS] rather than [RADIUS * TIME]
   vtan_tend = vtan_tend0(time|:,rad|:)
   vrad_tend = vrad_tend0(time|:,rad|:)
   mslp_tend = mslp_tend0(time|:,rad|:)

   vtan_plot = vtan_plot0(time|:,rad|:)
   vrad_plot = vrad_plot0(time|:,rad|:)
   vvel_plot = vvel_plot0(time|:,rad|:)
   pvrt_plot = pvrt_plot0(time|:,rad|:)
   mslp_plot = mslp_plot0(time|:,rad|:)

;=====================================================
; Create a panel plot of Figs. 6a and 6b from N11 
;=====================================================

  ; Manage output string depending on whether ERCs or VCs are analysed 
    if (erc .eq. 1) then 
     vstr="pv"
    else
     vstr="mslp"
    end if 

  ; Output file location and type
    if (lay .eq. 1) then 
     output = "$nep/nepartak/images/cyl_coords/hov_"+vstr+"_paper_"+dat+"_"+ens0+"_"+hy0+"_"+hy1+\
     	      "_sc"+sc
    elseif (lay .eq. 0) then 
     output = "$nep/nepartak/images/cyl_coords/hov_"+vstr+"_paper_"+dat+"_"+ens0+"_"+hy0+"_sc"+sc
    end if 
    wks = gsn_open_wks(opt,output)

  ; Set up panel plot
    panel = new(3,graphic)

    gsn_define_colormap(wks,"n11")

;==========================
; Options for plotting 
;==========================

  ; Tangential velocity
    opts_vtan                              = True
    opts_vtan@cnFillOn                     = True
    opts_vtan@cnLineLabelInterval          = 2.0
    opts_vtan@cnLineLabelFontHeightF       = 0.012
    opts_vtan@cnLineLabelBackgroundColor   = "transparent"
    opts_vtan@cnLineLabelPlacementMode     = "constant"
    opts_vtan@cnLinesOn                    = False
    opts_vtan@cnInfoLabelOn                = False
    opts_vtan@cnLevelSelectionMode         = "ExplicitLevels"

    if (basic .eq. 1) then 

     opts_vtan@cnLevels			   = (/5.0, 10.0, 15.0, 17.5, 20.0, \
                                                22.5, 25.0, 27.5, 30.0, \
                                                32.5, 35.0, 40.0/)
     opts_vtan@cnFillColors	           = (/-1,0,1,2,3,4,5,6,7,8,9,10,11/)

    else

     if (erc .eq. 0) then 
      if (tclr .eq. 0) then 
       opts_vtan@cnLevels                   = (/1.0, 2.0, 2.5, 3.0, 3.5, 4.0, \
                                                4.5, 5.0, 6.0, 7.0/)
       opts_vtan@cnFillColors               = (/0,2,3,4,5,6,7,8,9,10,11/)
      elseif (tclr .eq. 1) then 
       opts_vtan@cnLevels                   = (/0.5, 1.5, 2.0, 2.5, 3.0, 3.5, \
                                                4.0, 4.5, 5.0, 6.0/)
       opts_vtan@cnFillColors               = (/0,2,3,4,5,6,7,8,9,10,11/)
      elseif (tclr .eq. 2) then		    
       opts_vtan@cnLevels                   = (/0.5, 1.0, 1.5, 2.0, 2.5, 3.0, \
                                                3.5, 4.0, 4.5, 5.0/)
       opts_vtan@cnFillColors               = (/0,2,3,4,5,6,7,8,9,10,11/)
      end if 
     else
       opts_vtan@cnLevels                   = (/0.2, 0.3, 0.5, 0.7, 1.0, \
                                                1.2, 1.5, 2.0, 2.5, 3.0/)
       opts_vtan@cnFillColors               = (/0,2,3,4,5,6,7,8,9,10,11/)
     end if 

    end if 

    opts_vtan@gsnPaperOrientation          = "landscape"
    opts_vtan@tiMainString                 = ""
    opts_vtan@tiMainFontHeightF            = 0.0125
    opts_vtan@gsnLeftString                = ""
    opts_vtan@gsnRightString               = ""
    opts_vtan@gsnMaximize                  = True
    opts_vtan@lbLabelBarOn                 = True                  ; Labelbar on/off
    opts_vtan@lbBoxEndCapStyle             = "TriangleBothEnds"    ; Labelbar end shape
    opts_vtan@pmLabelBarWidthF             = 0.32                  ; Labelbar width
    opts_vtan@pmLabelBarHeightF            = 0.08                  ; Labelbar height
    opts_vtan@lbLabelFontHeightF           = 0.0125                ; Labelbar font size
    opts_vtan@lbLabelFont                  = "Helvetica"           ; Labelbar font
    opts_vtan@lbPerimOn                    = False                 ; Perimeter on/off
    opts_vtan@gsnDraw                      = False                 ; Do not draw plot
    opts_vtan@gsnFrame                     = False                 ; Do not advance frame
    opts_vtan@vpWidthF                     = 0.30                  ; Match width in N11
    opts_vtan@vpHeightF                    = 0.50                  ; Match height in N11

  ; Additional plotting resources
    opts_vtan@tiYAxisString                = "Time (h)"
    opts_vtan@tiXAxisString                = "Radius (km)"
    opts_vtan@tiXAxisFontHeightF           = 0.018
    opts_vtan@trYLog                       = False
    opts_vtan@trYMinF                      = ts      ; Starting point (T+...) for plot
    opts_vtan@trYMaxF                      = tf-2    ; Ending point (T+...) for plot
    opts_vtan@tmXTOn                       = "False" ; Turn off top x-axis TM
    opts_vtan@tmYROn                       = "False" ; Turn off right y-axis TM
    opts_vtan@gsnMaximize                  = True    ; Maximise plot size
    opts_vtan@gsnAddCyclic                 = False
    opts_vtan@tmXBMode                     = "Explicit"                ; Choose tm mode
    opts_vtan@tmXBValues                   = (/0.0,0.5,1.0,1.5,2.0/)   ; Tm positions
    opts_vtan@tmXBLabels                   = (/"0","50","100",\
    					       "150","200"/)	       ; Tm values
    opts_vtan@tmXBMinorOn                  = False

    opts_vtan@pmLabelBarOrthogonalPosF     = 0.10    ; Move lb up/down (higher/lower)
    opts_vtan@tiXAxisOffsetYF              = 0.10    ; Move x-axis title up/down (h/l)

    opts_vtan@trXMinF                      = 0.0                       ; Set x-axis min
    opts_vtan@trXMaxF                      = ar0                       ; Set x-axis max
    opts_vtan@tmYLMode                     = "Manual"
    opts_vtan@tmYLTickStartF               = opts_vtan@trYMinF         ; See L730
    opts_vtan@tmYLTickEndF                 = opts_vtan@trYMaxF         ; See L731
    opts_vtan@tmYLTickSpacingF             = 6
    opts_vtan@tmYLMinorOn                  = False

  ; Resources to overlay contours of mean tangential wind
    opts_cont				   = True

    opts_cont@vpWidthF                     = 0.30                  ; Match width in N11
    opts_cont@vpHeightF                    = 0.50                  ; Match height in N11

    opts_cont@cnFillOn                     = False       ; Filled contours on/off
    opts_cont@cnLineColor                  = "blue3"     ; Contour line colour
    opts_cont@cnLevelSelectionMode         = "ExplicitLevels"
    opts_cont@cnLevels                     = (/10,20,30,40,50,60,70/)

    opts_cont@tiYAxisString                = "Time (h)"
    opts_cont@tiXAxisString                = "Radius (km)"
    opts_cont@tiXAxisFontHeightF           = 0.018
    opts_cont@trYLog                       = False
    opts_cont@trYMinF                      = ts          ; Starting point (T+...) for plot
    opts_cont@trYMaxF                      = tf-2        ; Ending point (T+...) for plot
    opts_cont@tmXTOn                       = "False" 	 ; Turn off top x-axis TM
    opts_cont@tmYROn                       = "False" 	 ; Turn off right y-axis TM
    opts_cont@gsnMaximize                  = True    	 ; Maximise plot size
    opts_cont@gsnAddCyclic                 = False
    opts_cont@tmXBMode                     = "Explicit"                ; Choose tm mode
    opts_cont@tmXBValues                   = (/0.0,0.5,1.0,1.5/)       ; Tm positions
    opts_cont@tmXBLabels                   = (/"0","50","100","150"/)  ; Tm values
    opts_cont@tmXBMinorOn                  = False

    opts_cont@pmLabelBarOrthogonalPosF     = 0.10    ; Move lb up/down (higher/lower)
    opts_cont@tiXAxisOffsetYF              = 0.10    ; Move x-axis title up/down (h/l)

    opts_cont@trXMinF                      = 0.0                       ; Set x-axis min
    opts_cont@trXMaxF                      = 1.0                       ; Set x-axis max

    opts_cont@cnLineThicknessF             = 2.0         ; Contour line thickness
    opts_cont@cnInfoLabelOn                = False       ; Info label on/off
    opts_cont@gsnDraw                      = False       ; Do not draw the plot
    opts_cont@gsnFrame                     = False       ; Do not advance the frame
    opts_cont@cnLineLabelsOn               = True        ; Turn line labels on/off
    opts_cont@cnLineLabelFontColor         = "blue3"     ; Line label colour
    opts_cont@cnLineLabelFontHeightF       = 0.01        ; Default 0.006 for vpWidth=0.30
    opts_cont@cnLineLabelFontThicknessF    = 2.0         ; Line label thickness
    opts_cont@tiMainString                 = ""
    opts_cont@gsnLeftString                = ""
    opts_cont@gsnRightString               = ""

  ; Resources to overlay zero line
    opts_z                             = True
    opts_z@cnFillOn                    = False
    opts_z@cnLineColor                 = "gray82"
    opts_z@cnLevels                    = (/0.0/)
    opts_z@cnLineThicknessF            = 2.0
    opts_z@cnInfoLabelOn               = False
    opts_z@gsnDraw                     = False       ; Do not draw the plot
    opts_z@gsnFrame                    = False       ; Do no advance the frame
    opts_z@cnLineLabelsOn              = False       ; Turn off line labels
    opts_z@tiMainString                = ""
    opts_z@gsnLeftString               = ""
    opts_z@gsnRightString              = ""

  ; Mean sea level pressure tendency resources
    opts_mslp                          = opts_vtan
    delete([/opts_mslp@cnLevels,opts_mslp@cnFillColors/])

     if (pclr .eq. 0) then 
      opts_mslp@cnLevels                     = (/0.0, -0.5, -1.0, -1.5, -2.0, \
                                                 -2.5, -3.0, -4.0, -5.0, -6.0, -7.0/)
      opts_mslp@cnFillColors                 = (/12,11,10,9,8,7,6,5,4,3,2,0/)
     elseif (pclr .eq. 1) then 
      opts_mslp@cnLevels                     = (/0.0, -0.25, -0.5, -1.0, -1.5, \
                                                 -2.0, -2.5, -3.0, -4.0, -5.0, -6.0/)
      opts_mslp@cnFillColors                 = (/12,11,10,9,8,7,6,5,4,3,2,0/)
     elseif (pclr .eq. 2) then 
      opts_mslp@cnLevels                     = (/0.0, -0.25, -0.5, -0.75, -1.0, \
                                                 -1.5, -2.0, -2.5, -3.0, -4.0, -5.0/)
      opts_mslp@cnFillColors                 = (/12,11,10,9,8,7,6,5,4,3,2,0/)
     end if

  ; Potential vorticity resources 
    opts_pvrt                                = opts_vtan
    delete([/opts_pvrt@cnLevels, opts_pvrt@cnFillColors/])
    opts_pvrt@cnLevels                       = (/3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, \
                                                 15.0, 20.0, 25.0, 30.0, 35.0/)
    opts_pvrt@cnFillColors                   = (/0,2,3,4,5,6,7,8,9,10,11,12,13/)

  ; Vertical velocity resources
    opts_vvel                                = opts_vtan
    delete([/opts_vvel@cnLevels,opts_vvel@cnFillColors/])

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
    opts_rmw@xyLineColors                  = "black"
    opts_rmw@xyLineThicknesses             = 4.0
    opts_rmw@xyDashPatterns                = 0.0
    opts_rmw@gsnDraw                       = False       ; Do not draw the plot
    opts_rmw@gsnFrame                      = False       ; Do not advance the frame
    opts_rmw@tiMainString                  = ""
    opts_rmw@gsnLeftString                 = ""
    opts_rmw@gsnRightString                = ""

;========================
; Construct panel plot
;========================

  ; Read in RMW data from text file (first 3 values equal to zero)
    if (rmw .eq. 1) then 
     spd_plot = asciiread("$ar/text/mlev_"+dat+"_"+ens0+"_rmw_calc1.txt",(/numINPUT_c/),"float")
    elseif (rmw .eq. 0) then 
     spd_plot = asciiread("$ar/text/"+dat+"_"+ens0+"_rmw_calc1.txt",(/numINPUT_c/),"float")
    end if 

    ; Tangential wind (LHS panel)
      vtan0 = gsn_csm_contour(wks,vtan_tend,opts_vtan)    ; Tangential wind tendency (filled)
      zero0 = gsn_csm_contour(wks,vtan_tend,opts_z)    	  ; Zero line (vtan tendency)
      zero  = ColorNegDashZeroPosContour(zero0,"transparent","gray69","transparent")
      cont0 = gsn_csm_contour(wks,vtan_plot,opts_cont)    ; Tangential wind
      rmw0  = gsn_csm_xy(wks,spd_plot,time_int,opts_rmw)  ; RMW
      overlay(vtan0,zero)                                 ; Overlay vtan tendency zero line
      overlay(vtan0,rmw0)				  ; Overlay RMW
      overlay(vtan0,cont0)                                ; Overlay tangential wind
      panel(0) = vtan0                                    ; Plot first panel
      delete([/zero0,zero/])                              ; Tidy up 

    ; Potential vorticity/min sea level pressure (middle panel)
      if (erc .eq. 1) then 
       pvrt0 = gsn_csm_contour(wks,pvrt_plot,opts_pvrt)    ; Potential vorticity
       zero0 = gsn_csm_contour(wks,vtan_tend,opts_z)       ; Zero line (vtan tendency)          
       zero  = ColorNegDashZeroPosContour(zero0,"transparent","gray69","transparent")
       cont0 = gsn_csm_contour(wks,vtan_plot,opts_cont)    ; Tangential wind                    
       rmw0  = gsn_csm_xy(wks,spd_plot,time_int,opts_rmw)  ; RMW                                
       overlay(pvrt0,zero)
       overlay(pvrt0,cont0)                                ; Overlay tangential wind            
       overlay(pvrt0,rmw0)                                 ; Overlay RMW                        
       panel(1) = pvrt0                                    ; Plot second panel                  
       delete([/zero0,zero/])                              ; Tidy up 
      else
       mslp0 = gsn_csm_contour(wks,mslp_tend,opts_mslp)    ; MSLP tendency 
       zero0 = gsn_csm_contour(wks,vtan_tend,opts_z)       ; Zero line (vtan tendency)
       zero  = ColorNegDashZeroPosContour(zero0,"transparent","gray69","transparent")
       cont0 = gsn_csm_contour(wks,vtan_plot,opts_cont)    ; Tangential wind
       rmw0  = gsn_csm_xy(wks,spd_plot,time_int,opts_rmw)  ; RMW
       overlay(mslp0,zero)
       overlay(mslp0,cont0)                                ; Overlay tangential wind
       overlay(mslp0,rmw0)                                 ; Overlay RMW
       panel(1) = mslp0                                    ; Plot second panel 
       delete([/zero0,zero/])                              ; Tidy up 
      end if 

    ; Vertical velocity (RHS panel)
      vvel0 = gsn_csm_contour(wks,vvel_plot,opts_vvel)    ; Vertical velocity
      zero0 = gsn_csm_contour(wks,vvel_plot,opts_z)       ; Zero line (vertical velocity)
      zero  = ColorNegDashZeroPosContour(zero0,"transparent","gray69","transparent")
      rmw0  = gsn_csm_xy(wks,spd_plot,time_int,opts_rmw)  ; RMW
      overlay(vvel0,zero)                                 ; Overlay vvel zero line
      overlay(vvel0,rmw0)                                 ; Overlay RMW
      panel(2) = vvel0                                    ; Plot third panel
      delete([/zero0,zero/])                              ; Tidy up

;=============================================
; Finally, put everything together and plot 
;=============================================

    optsP                       = True
    optsP@gsnFrame              = False    ; Do not advance the frame
    optsP@gsnPanelLabelBar      = False    ; Turn on panel labelbar
    optsP@txString              = "Simulation "+ens0+" ("+dat+")"

    if (basic .eq. 0) then 
     print("No panel figure strings!")
    elseif (basic .eq. 1) then 
     optsP@gsnPanelFigureStrings = (/"a) vtan","b) mslp","c) vvel"/)
    end if 

    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.90
    optsP@gsnPanelBottom        = 0.10
    optsP@amJust                = "TopLeft"
    optsP@gsnPanelFigureStringsFontHeightF = 0.0099 ; Reduce label size (default 0.01)

    gsn_panel(wks,panel,(/1,3/),optsP)             ; Draw as a single plot
    frame(wks)

  ; Print output file path to screen
    print(output)

end 