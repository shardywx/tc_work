; Script calculates and writes max. standard deviation of vorticity
; and tangential wind tendency, to text files

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/setup_cyl.ncl"

; ncl opt=\"x11\" dat=\"02T12\" ens0=\"em11\" dist=3.0 calc=1
; ts=24 tf=96 tend=1 mlev0=13 lay=1 mlev1=25 cn0=\"geo_sm\" 
; sm=0 sc=3 dnt=1 ar=0.75 nr=16 tc_n11_fig2_step2a.ncl

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'ens0'  = ensemble simulation (em00-em11)
; 'dist'  = size of box following storm (degrees)
; 'calc'  = centre calculated on each level (0) or calculated offline (1)
; 'opt'   = output file format ("pdf" or "x11")
; 'ts'	  = start time index (0-119)
; 'tf'	  = end time index (0-119)
; 'ar'    = distance in degrees of outer radius (1.0, 2.0, etc)
; 'nr'    = number of radial circles between r = 0 and r = 'ar' (21, 41, 61, etc)
; 'tend'  = calculate wind tendency using 't2-t0 / 2' (2) or 't1-t0 / 1' (1)
; 'mlev0' = model level for plotting PV and vertical velocity (see L106-108)
; 'lay'   = calculate layer average (1)
; 'mlev1' = upper model level for calculating layer average (only valid when lay=1)
; 'cn0'   = method of storm tracking ["slp", "geo", "geo_sm", "vort"] 

; 'pc' stream --> 'it'
; 'centre'    --> 'it-1'
; 'pd' stream --> 'it-2'

begin

;==============================
; Prelash (reading in files)
;==============================

; Find 'pc' and 'pd' files (model height level data)
  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0
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
; 'dat'     = initialisation time (02T12, 03T00, etc)
; 'ens0'    = ensemble simulation (em00-em11)
; 'diri'    = path to input files (see above)
; 'dist'    = size of box following storm (degrees)
; 'mins'    = analyse 1-h (0) or 5-min (1) data

; Call external function 'st_rm'
  storm_rel = st_rm(dat, ens0, "$ar/text/", dist, 0)
  u_cyc     = storm_rel[0]                 ; Zonal wind
  v_cyc     = storm_rel[1]                 ; Meridional wind
  vel_cyc   = storm_rel[2]                 ; Vector wind
  numTIMES  = storm_rel[3]                 ; Number of times in file
  lat_arr   = storm_rel[4]                 ; Storm-following domain subset (lat)
  lon_arr   = storm_rel[5]                 ; Storm-following domain subset (lon)
  centre    = storm_rel[6]                 ; Storm centre at each hour (lat/lon)

; Create arrays for data and output information
  llbox     = toint(dist*50)                   ; Calculate domain size (grid points)
  dsize     = (/llbox,llbox/)                  ; Domain size (grid points)
  sym_arr   = new(numINPUT_d,float)            ; Values of symmetry parameter 

; Print out relevant information to user about storm tracking
  sc_arr      = (/"relative vorticity",\                          ; sc = 0 
                  "pressure",\                                    ; sc = 1 
                  "smoothed pressure",\                           ; sc = 2 
                  "extra smoothed pressure",\                     ; sc = 3 
                  "vort_min during symmetric phase",\             ; sc = 4 
                  "windspeed minimum",\                           ; sc = 5
                  "geopotential"/)                                ; sc = 6

  print("Using "+sc_arr(sc)+" to calculate storm centre position...") 

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

;==================================================================================
;;;;;;;;;;;;;;;;;;;;;;;;;;; Start of main part of script ;;;;;;;;;;;;;;;;;;;;;;;;;;
;==================================================================================

   print_clock("Starting main part of script...")

 ; ar = distance (degrees) of outer radial circle (1.0)
 ; nr = total number of radial circles (21)

   pii        = 3.14159265
   radii      = fspan(0.0,ar,nr)	; Size of radial circle array
   num_l      = 1			; Number of model levels

 ; Define constants
   rd    = 287.0   ; Gas constant for dry air (J/kg/K^2)
   cp0   = 1004.0  ; Specific heat of dry air at constant pressure (J/kg/K)
   kp    = 0.286   ; For Exner function calculations (0.286)
   pref  = 1000.0  ; Reference pressure p0 (1000 hPa)

 ; Create array for vorticity standard deviation and add metadata
   vort_sd1       = new((/dimsizes(radii),numTIMES/),float)
   vort_sd1!0     = "rad"
   vort_sd1!1     = "time"

;=============================
; Read in variables we want
;=============================

   print_clock("Working on ensemble member: "+ens0)
   print_clock("From simulation initialised at: "+dat)

   setvalues NhlGetWorkspaceObjectId
     "wsMaximumSize" : 1000000000
   end setvalues

 ; 1-dimensional array to hold max. 850 hPa tangential velocity for scatter plot
   vtan_max    = new((/numTIMES/),float)
   vtan_max!0  = "time"

 ; 1-dimensional array to hold tendency of 'vtan_max' for scatter plot
   vtan_tend   = new((/numTIMES/),float)
   vtan_tend!0 = "time"

   print("numTIMES = "+numTIMES)

;===============================
; Start loop over input files 
;===============================

 ; Create array to hold storm track positions using different methods 
   centre_new  = new((/7,num_l,numINPUT_d,2/),"float")

   c  = addfile(fili_c,"r")             ; Read in 'pc' stream [u,v,vort]
   d  = addfile(fili_d,"r")             ; Read in 'pd' stream [w,p,t,theta,kmh,tau]

   do it = ts, tf                   ; Loop over times in file

    print("Working on time: "+time_str(it)+" (T+"+it+")" )
    title_arr(it) = "Valid at "+time_str(it)+" (T+"+it+")"

  ; Get domain boundaries at each time from storm track (calculated off-line) 
    t0    = lat_arr(cn,it-1,0)
    t1    = lat_arr(cn,it-1,1)
    n0    = lon_arr(cn,it-1,0)
    n1    = lon_arr(cn,it-1,1)
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
     u0  = c->$u_varname$(it,13:25,{t0:t1},{n0:n1})     ; Zonal wind (layer)
     v0  = c->$v_varname$(it,13:25,{t0:t1},{n0:n1})     ; Meridional wind (layer)
    else
     u0  = c->$u_varname$(it,ilev,{t0:t1},{n0:n1})      ; Zonal wind (single level)
     v0  = c->$v_varname$(it,ilev,{t0:t1},{n0:n1})      ; Meridional wind (single level)
    end if

  ; Calculate storm-relative winds
    u = u0 - u_cyc(cn,it-1)
    v = v0 - v_cyc(cn,it-1)

  ; Add metadata from original horizontal wind arrays
    copy_VarCoords(u0,u)
    copy_VarAtts(u0,u)

    copy_VarCoords(v0,v)
    copy_VarAtts(v0,v)

  ; Calculate relative vorticity using centered finite differences
    vort0 = uv2vr_cfd(u, v, latC, lonC, 2)
    copy_VarCoords(u, vort0)
    vort0@units = "s~S~-1~N~"
    vort0@name = "Relative vorticity on model levels"

  ; Convert vorticity to standard units for plotting [10^-4 s-1]
    vort0 = vort0 * (10 ^ 4)

  ; EDIT -- is this part really required? 
    if (lay .eq. 1) then
     vort       = vort0(:,0:dsize(0)-1,0:dsize(1)-1)
    elseif (lay .eq. 0) then 
     vort       = vort0(0:dsize(0)-1,0:dsize(1)-1)
    end if

;======================================================
; Read in variables from 'pd' stream (model levels)
;======================================================

  ; Read in basic variables from 'pd' stream
    lon1D  = d->longitude_1({n0:n1})                   ; longitude [xxx grid points]
    lat1D  = d->latitude_1({t0:t1})                    ; latitude  [xxx grid points]
    lonD   = d->longitude({n0:n1})                     ; longitude [xxx grid points]
    latD   = d->latitude({t0:t1})                      ; latitude  [xxx grid points]
    hybD   = d->hybrid_ht(:)                           ; 63 model ('rho') levels

  ; Read in temperature, pressure, theta, vertical velocity and geopotential height
  ; All variables below -- [hybD | 63] * [latD | 150] * [lon1D | 151]

    if (lay .eq. 1) then 
     t      = d->temp(it-2,13:25,{t0:t1},{n0:n1})        ; Temperature (K)
     p      = d->p(it-2,13:25,{t0:t1},{n0:n1})           ; Pressure (Pa)
     z      = d->ht(it-2,13:25,{t0:t1},{n0:n1})          ; Geopotential height (m)
     pv     = d->field83(it-2,13:25,{t0:t1},{n0:n1})     ; Potential vorticity (PVU)
    elseif (lay .eq. 0) then 
     t      = d->temp(it-2,ilev,{t0:t1},{n0:n1})         ; Temperature (K)
     p      = d->p(it-2,ilev,{t0:t1},{n0:n1})            ; Pressure (Pa)
     z      = d->ht(it-2,ilev,{t0:t1},{n0:n1})           ; Geopotential height (m)
     pv     = d->field83(it-2,ilev,{t0:t1},{n0:n1})      ; Potential vorticity (PVU)
    end if 

  ; Convert units of pressure and geopotential height before continuing
    p  = p / 100           ; Convert to 'hPa'
    z  = z / 10            ; Convert to 'dam'
    pv = pv * (10 ^ 6)     ; Convert to 'PVU' [10 ^ 6 K m2 s-1 kg-1]

  ; Add metadata
    p@units  = "hPa"
    z@units  = "dam"
    pv@units = "PVU"

    dx       = lonD(1) - lonD(0) ; Grid spacing (longitude)
    dy       = latD(1) - latD(0) ; Grid spacing (latitude)

;===========================================
; Analyse variables on chosen model level
;===========================================

    if (lay .eq. 1) then 

   ; 'pc' stream
     u_plane       = dim_avg_n_Wrap(u(:,:,:), 0)
     v_plane       = dim_avg_n_Wrap(v(:,:,:), 0)
     vort_plane    = dim_avg_n_Wrap(vort(:,:,:), 0)

   ; 'pd' stream
     prs_plane     = dim_avg_n_Wrap(p(:,:,:), 0)
     pv_plane      = dim_avg_n_Wrap(pv(:,:,:), 0)

     hy0           = sprintf("%0.0f",hybD(13))
     hy1           = sprintf("%0.0f",hybD(25))
     print("Working on layer between "+hy0+ " m and "+hy1+" m")

    else

     hy2           = sprintf("%0.0f",hybD(ilev))
     print("Working on model level: "+hy2+ " m")

   ; 'pc' stream
     u_plane       = u(:,:)
     v_plane       = v(:,:)
     vort_plane    = vort(:,:)

   ; 'pd' stream
     prs_plane     = p(:,:)
     pv_plane      = pv(:,:)

    end if 

  ; Max and min pressure on height level
    pmin          = toint(min(prs_plane))
    pmax          = toint(max(prs_plane))

  ; Calculate vector windspeed and smooth pressure and vorticity
    spd_plane     = sqrt( (u_plane ^ 2) + (v_plane ^ 2) )
    copy_VarCoords(u_plane, spd_plane)
    vort_pl_smth  = smth9_Wrap(vort_plane, 0.5, 0.5, True)
    prs_pl_smth   = smth9_Wrap(prs_plane, 0.5, 0.5, True)
    prs_pl_smth2  = smth9_Wrap(prs_pl_smth, 0.5, 0.5, True)

;====================================================  
; Calculate storm centre using variety of methods
;====================================================

  ; Use already-calculated centre as a base for further calculations  
    lt1           = centre(cn,it-1,0) - 1.0
    lt2           = centre(cn,it-1,0) + 1.0
    ln1           = centre(cn,it-1,1) - 1.0
    ln2           = centre(cn,it-1,1) + 1.0

  ; RELATIVE VORTICITY MAXIMUM (sc=0)
    vort_plane0   = vort_plane({lt1:lt2},{ln1:ln2})
    vort_max0     = max(vort_plane0)
    dims_v0       = dimsizes(vort_plane0)
    vort_1d0      = ndtooned(vort_plane0)

  ; Find index of minimum --> (0,0) = latitude, (0,1) = longitude
    inds_v0       = ind_resolve(maxind(vort_1d0),dims_v0)

  ; Retrieve lat/lon information from grid subset
    ltN           = vort_plane0&latitude
    lnN           = vort_plane0&longitude

    lat_max_v0    = ltN(0) + (dy * inds_v0(0,0))
    lon_max_v0    = lnN(0) + (dx * inds_v0(0,1))
    print("Centre (vort): "+lat_max_v0+" degrees N, "+lon_max_v0+" degrees E")

    centre_new(0,0,it-2,0) = lat_max_v0
    centre_new(0,0,it-2,1) = lon_max_v0
    delete([/ltN, lnN/])

  ; PRESSURE MINIMUM (sc=1)
    prs_plane0     = prs_plane({lt1:lt2},{ln1:ln2})
    prs_min0       = min(prs_plane0)
    dims_p0        = dimsizes(prs_plane0)
    prs_1d0        = ndtooned(prs_plane0)

  ; Find index of minimum --> (0,0) = latitude, (0,1) = longitude 
    inds_p0        = ind_resolve(minind(prs_1d0),dims_p0)

  ; Retrieve lat/lon information from grid subset
    ltN            = prs_plane0&latitude
    lnN            = prs_plane0&longitude_1

    lat_max_p0     = ltN(0) + (dy * inds_p0(0,0))
    lon_max_p0     = lnN(0) + (dx * inds_p0(0,1))
    print("Centre (prs): "+lat_max_p0+" degrees N, "+lon_max_p0+" degrees E")

    centre_new(1,0,it-2,0) = lat_max_p0
    centre_new(1,0,it-2,1) = lon_max_p0
    delete([/ltN, lnN/])

  ; SMOOTHED PRESSURE MINIMUM (sc=2)
    prs_smth       = prs_pl_smth({lt1:lt2},{ln1:ln2})
    prs_min_smth   = min(prs_smth)
    dims_p1        = dimsizes(prs_smth)
    prs_1d1        = ndtooned(prs_smth)

  ; Find index of minimum --> (0,0) = latitude, (0,1) = longitude
    inds_p1        = ind_resolve(minind(prs_1d1),dims_p1)

  ; Retrieve lat/lon information from grid subset 
    ltN            = prs_smth&latitude
    lnN            = prs_smth&longitude_1

    lat_max_p1     = ltN(0) + (dy * inds_p1(0,0))
    lon_max_p1     = lnN(0) + (dx * inds_p1(0,1))
    print("Centre (prs_smth): "+lat_max_p1+" degrees N, "+lon_max_p1+" degrees E")

    centre_new(2,0,it-2,0) = lat_max_p1
    centre_new(2,0,it-2,1) = lon_max_p1
    delete([/ltN, lnN/])

  ; EXTRA SMOOTHED PRESSURE MINIMUM (sc=3) 
    prs_smth2      = prs_pl_smth2({lt1:lt2},{ln1:ln2})
    prs_min_smth2  = min(prs_smth2)
    dims_p2        = dimsizes(prs_smth2)
    prs_1d2        = ndtooned(prs_smth2)

  ; Find index of minimum --> (0,0) = latitude, (0,1) = longitude 
    inds_p2        = ind_resolve(minind(prs_1d2),dims_p2)

  ; Retrieve lat/lon information from grid subset 
    ltN            = prs_smth2&latitude
    lnN            = prs_smth2&longitude_1

    lat_max_p2     = ltN(0) + (dy * inds_p2(0,0))
    lon_max_p2     = lnN(0) + (dx * inds_p2(0,1))
    print("Centre (prs_smth_2): "+lat_max_p2+" degrees N, "+lon_max_p2+" degrees E")

    centre_new(3,0,it-2,0) = lat_max_p2
    centre_new(3,0,it-2,1) = lon_max_p2
    delete([/ltN, lnN/])

  ; DURING RINGLIKE PHASE, FIND THE CENTRE OF THE 'DOUGHNUT' (sc=4)
  ; Use centre calculated using smoothed pressure (above) as base
    r0             = 0.05
    lt3            = lat_max_p2 - r0
    lt4            = lat_max_p2 + r0
    ln3            = lon_max_p2 - r0
    ln4            = lon_max_p2 + r0

    vort_pl_new    = vort_plane({lt3:lt4},{ln3:ln4})
    vort_min_new   = min(vort_pl_new)
    dims_v3        = dimsizes(vort_pl_new)
    vort_1d3       = ndtooned(vort_pl_new)
 
  ; Find index of minimum --> (0,0) = latitude, (0,1) = longitude   
    inds_v3        = ind_resolve(minind(vort_1d3),dims_v3)

  ; Retrieve lat/lon information about grid subset 
    ltN            = vort_pl_new&latitude
    lnN            = vort_pl_new&longitude

    lat_min_v3     = ltN(0) + (dy * inds_v3(0,0))
    lon_min_v3     = lnN(0) + (dx * inds_v3(0,1))
    print("Centre (vort_min): "+lat_min_v3+" degrees N, "+lon_min_v3+" degrees E")

    centre_new(4,0,it-2,0) = lat_min_v3
    centre_new(4,0,it-2,1) = lon_min_v3
    delete([/ltN, lnN/])

  ; WINDSPEED MINIMUM (sc=5)  
    spd_pl_new     = spd_plane({lt3:lt4},{ln3:ln4})
    spd_min_new    = min(spd_pl_new)
    dims_s3        = dimsizes(spd_pl_new)
    spd_1d3        = ndtooned(spd_pl_new)

  ; Find index of minimum --> (0,0) = latitude, (0,1) = longitude 
    inds_s3        = ind_resolve(minind(spd_1d3),dims_s3)

  ; Retrieve lat/lon information from subset of grid 
    ltN            = spd_pl_new&latitude
    lnN            = spd_pl_new&longitude

    lat_min_s3     = ltN(0) + (dy * inds_s3(0,0))
    lon_min_s3     = lnN(0) + (dx * inds_s3(0,1))
    print("Centre (spd_min): "+lat_min_s3+" degrees N, "+lon_min_s3+" degrees E")

    centre_new(5,0,it-2,0) = lat_min_s3
    centre_new(5,0,it-2,1) = lon_min_s3
    delete([/ltN, lnN/])

  ; PRE-CALCULATED STORM POSITION (sc=6) 
    centre_new(6,0,it-2,0) = centre(cn,it-1,0)
    centre_new(6,0,it-2,1) = centre(cn,it-1,1)

;=====================================================================  
; Call external function to switch to cylindrical coordinate system
;=====================================================================  

  ; Call 'setup_cyl' and output pressure/vorticity on cylindrical grid
    cyl_arr  = setup_cyl(it, dist, t0, t1, n0, n1, lonC, latC, lonD, latD, \
                         "vort", vort_plane, prs_plane, centre_new, sc, sm, ar, nr)
    prs_int  = cyl_arr[0]
    vort_int = cyl_arr[1] 
    vort_azi = cyl_arr[2] 
    xpos     = cyl_arr[3]
    ypos     = cyl_arr[4]
    lat_max  = cyl_arr[5]
    lon_max  = cyl_arr[6]
    radii    = cyl_arr[7]
    rad_size = cyl_arr[8]

  ; Is our storm ringlike or monopolar? 
    sym_arr(it) = avg( vort_azi(0,0) + vort_azi(1,0) ) / avg ( vort_azi(2,0) + vort_azi(3,0) )
    if (sym_arr(it) .gt. 1) then
     print("Storm is in asymmetric phase (ratio = "+sym_arr(it)+")")
    else
     print("Storm in in symmetric phase (ratio = "+sym_arr(it)+")")
    end if    

  ; If storm is ringlike, recalculate the storm centre position
    if (dnt .eq. 1 .and. sym_arr(it) .lt. 1) then

     print("Re-calculating storm centre to fit ringlike vorticity profile...")

   ; Select the 'vort_min' method of storm tracking 
     sc = 4

   ; Call the external function again 
     cyl_arr  = setup_cyl(it, dist, t0, t1, n0, n1, lonC, latC, lonD, latD, \
                          "vort", vort_plane, prs_plane, centre_new, sc, sm, ar, nr)
     prs_int  = cyl_arr[0] 
     vort_int = cyl_arr[1]
     vort_azi = cyl_arr[2]
     xpos     = cyl_arr[3]
     ypos     = cyl_arr[4]
     lat_max  = cyl_arr[5]
     lon_max  = cyl_arr[6]
     radii    = cyl_arr[7]
     rad_size = cyl_arr[8] 

    end if 

  ; Also create array to hold relative vorticity standard deviation
    vort_sd0 = new((/dimsizes(radii)/),float)
    copy_VarAtts(vort_plane,vort_sd0)
    vort_sd0@description = "Relative vorticity standard deviation"
    vort_sd0@units       = ""
    vort_sd0!0           = "rad"
    vort_sd0&rad         = radii

  ; Azimuthally average PV around the storm and calculate standard deviation
    do irad = 0, dimsizes(radii)-1
     vort_sd0(irad)      = dim_stddev_n(vort_int(irad,:),0)
    end do

  ; Input into larger arrays
    vort_sd1(:,it)   = vort_sd0(:)     ; Standard deviation around radial circles 

  ; Tidy up
    delete([/vort_sd0, u_plane, v_plane, pv_plane, vort_plane, prs_plane/])
    delete([/vort_plane0, prs_plane0, prs_pl_smth, prs_pl_smth2, vort_pl_smth/])
    delete([/vort_1d0, prs_1d0, prs_1d1, prs_1d2, prs_smth, prs_smth2/])

  ; Tidy up (avoid dimension size errors)
    delete([/lonC,latC,lon1C,lat1C,u,v,vort,u0,v0,vort0/])
    delete([/lonD,latD,lon1D,lat1D,t,p,z,pv/])

  end do     ; End time loop (do it = ts, tf)

;==================================
; Write out values to text file 
;==================================

; Output file directory path (layer-averaged or single layer)
  diri_out  = "$ar/text/cp_sd_sc"+sc+"_"+dat+"_"+ens0+"_"

; Calculate maximum standard deviation
  vort_sd    = new(numTIMES,float)                ; Array to hold values below
  vort_sd(:) = dim_max_n(vort_sd1(:,:),0)         ; 1D array [numTIMES]

  if (lay .eq. 1) then
   print("Max. standard deviation of layer-averaged ("+\
         hy0+"-"+hy1+" m) relative vorticity = "+vort_sd)
  else
   print("Max. standard deviation of "+hy+" m relative vorticity = "+vort_sd)  
  end if 

; Write out relative vorticity standard deviation
  if (lay .eq. 1) then
   vort_out = diri_out+"vort_lay_"+hy0+"m_"+hy1+"m.txt"
   asciiwrite(vort_out,vort_sd)
  elseif (lay .eq. 0) then
   vort_out = diri_out+"vort_"+hy+"m.txt"
   asciiwrite(vort_out,vort_sd)
  end if

; Only write out tangential wind if user sets 'vt0' = 1
  if (vt0 .eq. 1) then 

  ; Read in the maximum tangential wind from text file (calculated on model levels)
    vtan_max    = asciiread("$ar/text/mlev_"+dat+"_"+ens0+"_vtan_max_"+cn0+".txt",\ 
  			   (/numTIMES/),"float")	

  ; Calculate the tendency of the maximum tangential wind
    do it = 0, numTIMES-1

     tp1 = min((/it+1,numTIMES-1/))
     tm1 = max((/it-1,0/))  

     if (tend .eq. 2) then 
      vtan_tend(it) = (vtan_max(tp1) - vtan_max(tm1)) / 2
     elseif (tend .eq. 1) then 
      vtan_tend(it) = (vtan_max(it) - vtan_max(tm1)) / 1
     end if 
  
    end do 

  ; Write out tangential wind tendency (calculated from text file)
    ten_out = diri_out+hy+"m_tend.txt"
    diri_out  = "$ar/text/cp_sd_sc"+sc+"_"+dat+"_"+ens0+"_"
    asciiwrite(ten_out,vtan_tend)

  end if 

end 