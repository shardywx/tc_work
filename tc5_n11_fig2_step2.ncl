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

; ncl opt=\"x11\" dat=\"02T12\" ens0=\"em11\" dist=3.0 calc=1
; ts=923 tf=995 ar=1.5 nr=31 vortO=1 tend=1 mlev0=13 lay=1 mlev1=25 cn0=\"slp\" 
; tc5_n11_fig2_step2.ncl

; T+12    = 131 (144)                                                               
; T+24    = 275 (288)                                                               
; T+36    = 419 (432)                                                               
; T+48    = 563 (576)                                                               
; T+60    = 707 (720)                                                               
; T+72    = 851 (864)                                                               
; T+84    = 995 (1008)

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'ens0'  = ensemble simulation (em00-em11)
; 'dist'  = size of box following storm (degrees)
; 'calc'  = centre calculated on each level (0) or calculated offline (1)
; 'opt'   = output file format ("pdf" or "x11")
; 'ts'	  = start time index (0-119)
; 'tf'	  = end time index (0-119)
; 'ar'    = distance in degrees of outer radius (1.0, 2.0, etc)
; 'nr'    = number of radial circles between r = 0 and r = 'ar' (21, 41, 61, etc)
; 'vortO' = choose relative (1) or potential (0) vorticity
; 'tend'  = calculate wind tendency using 't2-t0 / 2' (2) or 't1-t0 / 1' (1)
; 'mlev0' = model level for plotting PV and vertical velocity (see L106-108)
; 'lay'   = calculate layer average (1)
; 'mlev1' = upper model level for calculating layer average (only valid when lay=1)
; 'cn0'   = method of storm tracking ["slp" or "slpf"]

begin

;==============================
; Prelash (reading in files)
;==============================

; Find 'pc' and 'pd' stream data to read in and plot                                       
  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0+"/vc"
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0
  fili_c1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pc.nc")
  fili_d1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pd.nc")
  fili_k1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pk.nc")

  fili_c      = diri+"/"+fili_c1
  fili_d      = diri+"/"+fili_d1
  fili_k      = diri+"/"+fili_k1

  if ( .not. fileexists (fili_d) ) then
   print("Input file (pd stream) missing! Exiting script...")
   exit()
  elseif ( .not. fileexists (fili_c) ) then
   print("Input file (pc stream) missing! Exiting script...")
   exit()
  end if

; Calculate number of times in file using built-in function ('nc_times')          
  numINFO_c   = nc_times(fili_c)
  numINFO_d   = nc_times(fili_d)
  numINFO_k   = nc_times(fili_k)

  numINPUT_c  = numINFO_c[0]            ; Number of times in 'pc' stream          
  time_c      = numINFO_c[1]            ; Array of times from 'pc' stream file    

  numINPUT_d  = numINFO_d[0]            ; Number of times in 'pd' stream          
  time_d      = numINFO_d[1]            ; Array of times from 'pd' stream file    

  numINPUT_k  = numINFO_k[0]            ; Number of times in 'pk' stream              
  time_k      = numINFO_k[1]            ; Array of times from 'pk' stream file 

; Define number of vertical levels in input files                                 
  num_l       = 63

; Calculate storm motion using built-in function ('st_rm')
; 'dat'     = initialisation time (02T12, 03T00, etc)
; 'ens0'    = ensemble simulation (em00-em11)
; 'diri'    = path to input files (see above)
; 'dist'    = size of box following storm (degrees)
; 'mins'    = analyse 1-h (0) or 5-min (1) data

; Call external function 'st_rm'
  storm_rel = st_rm(dat, ens0, "$ar/text/", dist, 1)
  u_cyc     = storm_rel[0]                 ; Zonal wind
  v_cyc     = storm_rel[1]                 ; Meridional wind
  vel_cyc   = storm_rel[2]                 ; Vector wind
  lat_arr   = storm_rel[4]                 ; Storm-following domain subset (lat)
  lon_arr   = storm_rel[5]                 ; Storm-following domain subset (lon)
  centre    = storm_rel[6]                 ; Storm centre at each hour (lat/lon)

; Create arrays for data and output information
  llbox     = toint(dist*50)                   ; Calculate domain size (grid points)
  dsize     = (/llbox,llbox/)                  ; Domain size (grid points)
  d0        = dsize(0)
  d1        = dsize(1)

  print("Working on ensemble member: "+ens0+" ("+dat+")")

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 10000000000
  end setvalues

;==========================================================================       
; Create date/time string arrays for all times in file (use 'pc' stream)          
;==========================================================================       

  month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
                    "Jul","Aug","Sep","Oct","Nov","Dec"/)

  utc_date   = cd_calendar(time_c, 0)

  year       = tointeger(utc_date(:,0))
  month      = tointeger(utc_date(:,1))
  day        = tointeger(utc_date(:,2))
  hour       = tointeger(utc_date(:,3))
  minute     = tointeger(utc_date(:,4))
  second     = utc_date(:,5)

; Create array to hold finished date strings                                      
  hr_min = new(numINPUT_c,string)
  min0   = (/4,9,14,19,24,29,34,39,44,49,54,59/)      ; Minute array (before)     
  min1   = (/5,10,15,20,25,30,35,40,45,50,55,0/)      ; Minute array (after)      

; Before correcting, find the indices of the elements where 'minute' = 59         
  hr_ind = ind(minute .eq. 59)

; Correct for errors in the code (round up values of 'minute' to multiples of 5)  
  do i = 0, numINPUT_c-1
   do m = 0, dimsizes(min0)-1

    if (minute(i) .eq. min0(m) ) then
     minute(i) = min1(m)
    end if

   end do
  end do

; Edit all incorrect values of 'hour' (where 'minute' previously = 59)            
  do j = 0, dimsizes(hr_ind)-1
   hour(hr_ind(j)) = hour(hr_ind(j)) + 1
  end do
  delete(j)

; Create string array of all values in 'minute'                                   
  mins = tostring(minute)

; Create 'hhmm' strings for each time (for output and title strings)              
  do i = 0, numINPUT_c-1

 ; Add zeros in front of single-digit hours                                       
   if (hour(i) .lt. 10) then

    if (minute(i) .lt. 10) then
     hr_min(i) = "0"+hour(i)+"0"+mins(i)
    elseif (minute(i) .ge. 10) then
     hr_min(i) = "0"+hour(i)+mins(i)
    end if

   else

  ; Also add zeros in front of single-digit minutes                               
    if (minute(i) .lt. 10) then
     hr_min(i) = hour(i)+"0"+mins(i)
    elseif (minute(i) .ge. 10) then
     hr_min(i) = hour(i)+mins(i)
    end if

   end if

  end do
  delete(i)

  time_str  = new(numINPUT_c,string)
  time_arr  = new(numINPUT_c,string)
  title_arr = new(numINPUT_c,string)

  time_str  = hr_min + " UTC " + sprinti("%0.2i ", day) + month_abbr(month)
  time_arr  = sprinti("%0.2i", day) + month_abbr(month) + "_" + hr_min + "Z"
  title_arr = "Valid at "+time_str

; EDIT -- work out way to pause the script and start again if wrong times are chosen
  print_clock("Start time = "+time_arr(ts) )
  print_clock("End time = "+time_arr(tf) )

;=======================================                                          
; Choose which storm track to read in                                             
;=======================================                                          

  if (cn0 .eq. "slp") then
   cn = 0
  elseif (cn0 .eq. "slpf") then
   cn = 1
  end if

;==================================================================================
;;;;;;;;;;;;;;;;;;;;;;;;;;; Start of main part of script ;;;;;;;;;;;;;;;;;;;;;;;;;;
;==================================================================================

  print_clock("Starting main part of script...")

; ar = distance (degrees) of outer radial circle (1.0)
; nr = total number of radial circles (21)

  thetas     = new(73,float)
  pii        = 3.14159265
  radii      = fspan(0.0,ar,nr)	     ; Size of radial circle array
  thetas     = fspan(0.0,360.0,73)   ; Array of azimuth angles (0, 5, 10, etc)
  thetas_rad = (pii/180)*thetas      ; Degrees to radians
  rad_size   = dimsizes(radii)-1     ; Index for later calculation

; Define constants
  rd    = 287.0   ; Gas constant for dry air (J/kg/K^2)
  cp0   = 1004.0  ; Specific heat of dry air at constant pressure (J/kg/K)
  kp    = 0.286   ; For Exner function calculations (0.286)
  pref  = 1000.0  ; Reference pressure p0 (1000 hPa)

;===================================================
; Add metadata for azimuthally averaged variables 	
;===================================================

  if (vortO .eq. 0) then 

   pv_azi           = new((/dimsizes(radii)/),float)
   pv_azi!0         = "rad"
   pv_azi&rad       = radii
   pv_azi@description = "Azimuthally averaged PV"
   pv_azi@units     = "10~S~6~N~ K m~S~2~N~ s~S~-1~N~ kg~S~-1~N~"

  elseif (vortO .eq. 1) then 

   vort_azi         = new((/dimsizes(radii)/),float)
   vort_azi!0       = "rad"
   vort_azi&rad     = radii
   vort_azi@description = "Azimuthally averaged relative vorticity"
   vort_azi@units   = "s~S~-1~N~"

  end if 

;==============================================================
; Add metadata for larger arrays to hold values at all times
;==============================================================

  if (vortO .eq. 0) then 

  ; Potential vorticity
    pv_plot        = new((/dimsizes(radii),numINPUT_c/),float)
    pv_plot!0      = "rad"
    pv_plot!1      = "time"

  ; Potential vorticity standard deviation
    pv_sd1         = new((/dimsizes(radii),numINPUT_c/),float)
    pv_sd1!0       = "rad"
    pv_sd1!1       = "time"

  elseif (vortO .eq. 1) then 

  ; Relative vorticity
    vort_plot      = new((/dimsizes(radii),numINPUT_c/),float)
    vort_plot!0    = "rad"
    vort_plot!1    = "time"

  ; Relative vorticity standard deviation
    vort_sd1       = new((/dimsizes(radii),numINPUT_c/),float)
    vort_sd1!0     = "rad"
    vort_sd1!1     = "time"

  end if 

;=============================
; Read in variables we want
;=============================

 ; 1-dimensional array to hold max. 850 hPa tangential velocity for scatter plot
   vtan_max    = new((/numINPUT_c/),float)
   vtan_max!0  = "time"

 ; 1-dimensional array to hold tendency of 'vtan_max' for scatter plot
   vtan_tend   = new((/numINPUT_c/),float)
   vtan_tend!0 = "time"

;===============================
; Start loop over input files 
;===============================

   c  = addfile(fili_c,"r")             ; Read in 'pc' stream [u,v,vort]
   d  = addfile(fili_d,"r")             ; Read in 'pd' stream [p]
   k  = addfile(fili_k,"r")             ; Read in 'pk' stream [PV]


 ; Loop over times in file
   do it = ts, tf

    print("Working on time: "+time_str(it))
    itt = it
    title_arr(it) = "Valid at "+time_str(it)+" (T+"+itt+")"

  ; Get domain boundaries at each time from storm track (calculated offline)
    t0    = lat_arr(cn,it,0)
    t1    = lat_arr(cn,it,1)
    n0    = lon_arr(cn,it,0)
    n1    = lon_arr(cn,it,1)
    print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

;=====================================================
; Read in variables from 'pc' stream (model levels)
;=====================================================

  ; Latitude, longitude and model levels
    lonC   = c->longitude({n0:n1})                ; longitude ['d1' grid points]
    latC   = c->latitude({t0:t1})                 ; latitude  ['d0' grid points]
    hybC   = c->hybrid_ht(:)                      ; 63 model ('theta') levels

    u_varname = "x-wind"
    v_varname = "y-wind"

  ; Read in horizontal wind components
  ; Both variables below -- [hybC | 63] * [latC | 150] * [lonC | 150]

  ; Layer
    if (lay .eq. 1) then

     u0  = c->$u_varname$(it,mlev0:mlev1,{t0:t1},{n0:n1})
     v0  = c->$v_varname$(it,mlev0:mlev1,{t0:t1},{n0:n1})

  ; Single level
    else

     u0  = c->$u_varname$(it,mlev0,{t0:t1},{n0:n1})
     v0  = c->$v_varname$(it,mlev0,{t0:t1},{n0:n1})

    end if

;===================================
; Calculate storm-relative winds
;===================================

    u = u0 - u_cyc(cn,it)
    v = v0 - v_cyc(cn,it)

  ; Add metadata from original horizontal wind arrays
    copy_VarCoords(u0,u)
    copy_VarAtts(u0,u)

    copy_VarCoords(v0,v)
    copy_VarAtts(v0,v)

;==================================================================
; Calculate relative vorticity using centered finite differences
;==================================================================

  ; Option '2' --> boundary points estimated using one-sided difference scheme
    vort0 = uv2vr_cfd(u, v, latC, lonC, 2)
    copy_VarCoords(u, vort0)
    vort0@units = "s~S~-1~N~"
    vort0@name = "Relative vorticity on model levels"

    if (lay .eq. 1) then
     vort       = vort0(:, 0:d0-1, 0:d1-1)
    elseif (lay .eq. 0) then 
     vort       = vort0(0:d0-1, 0:d1-1)
    end if

;======================================================
; Read in variables from 'pd' stream (model levels)
;======================================================

  ; Latitude, longitude and model levels
    lonD   = d->longitude({n0:n1})                     ; longitude [xxx grid points]
    latD   = d->latitude({t0:t1})                      ; latitude  [xxx grid points]
    hybD   = d->hybrid_ht(:)                           ; 63 model ('rho') levels

    hy0    = sprintf("%0.0f",hybD(mlev0))

    if (lay .eq. 1) then 
     hy1    = sprintf("%0.0f",hybD(mlev1))
     print("Working on layer between "+hy0+ " m and "+hy1+" m")
    else
     print("Working on model level: "+hy0+ " m")
    end if 

  ; Read in pressure and geopotential height 
  ; All variables below -- [hybD | 63] * [latD | 150] * [lon1D | 151]

    if (lay .eq. 1) then 
     p       = d->p(it,mlev0:mlev1,{t0:t1},{n0:n1})           ; Pressure (Pa)
     z       = d->ht(it,mlev0:mlev1,{t0:t1},{n0:n1})          ; Geopotential height (m)
    elseif (lay .eq. 0) then 
     p       = d->p(it,mlev0,{t0:t1},{n0:n1})                 ; Pressure (Pa)
     z       = d->ht(it,mlev0,{t0:t1},{n0:n1})                ; Geopotential height (m)
    end if 

    dx0      = lonD(1) - lonD(0) ; Grid spacing (longitude)
    dy0      = latD(1) - latD(0) ; Grid spacing (latitude)

  ; EDIT -- NOT SURE 'z' IS ANY USE HERE ANYWAY (MODEL LEVELS)

;===============
; 'pk' stream 
;===============

  ; Read in potential vorticity
    if (lay .eq. 1) then 
     pv      = k->field83(it,mlev0:mlev1,{t0:t1},{n0:n1})
    elseif (lay .eq. 0) then 
     pv      = k->field83(it,mlev0,{t0:t1},{n0:n1}) 
    end if 

  ; Convert units before continuing                              
    p  = p / 100           ; Convert to 'hPa'
    z  = z / 10            ; Convert to 'dam'
    pv = pv * (10 ^ 6)     ; Convert to 'PVU' [10 ^ 6 K m2 s-1 kg-1]
    vort = vort * (10 ^ 4) ; Convert to standard units for plotting [10^-4 s-1]

  ; Add metadata
    p@units  = "hPa"
    z@units  = "dam"
    pv@units = "PVU"

;===========================================
; Analyse variables on chosen model level
;===========================================

    hy = sprintf("%0.0f",hybD(mlev0))
    print("Working on model level " + mlev0+ " ("+hy+" m AGL)")

    if (lay .eq. 1) then 
     u_plane       = dim_avg_n_Wrap(u(:,:,:),0)
     v_plane       = dim_avg_n_Wrap(v(:,:,:),0)
     vort_plane    = dim_avg_n_Wrap(vort(:,:,:),0)
     prs_plane     = dim_avg_n_Wrap(p(:,:,:),0)
     pv_plane      = dim_avg_n_Wrap(pv(:,:,:),0)
    elseif (lay .eq. 0) then 
     u_plane       = u(:,:)
     v_plane       = v(:,:)
     vort_plane    = vort(:,:)
     prs_plane     = p(:,:)
     pv_plane      = pv(:,:)
    end if 

  ; Max and min pressure on height level
    pmin           = toint(min(prs_plane))
    pmax           = toint(max(prs_plane))

  ; Calculate location of vorticity maximum (storm centre) on model level
    vort_max       = max(vort_plane)
    vort_smth      = smth9_Wrap(vort_plane, 0.5, 0.5, True)
    vort_max_smth  = max(vort_smth)

  ; Find index of relative vorticity maximum
    dims           = dimsizes(vort_plane)
    vort1d         = ndtooned(vort_plane)
    inds           = ind_resolve(maxind(vort1d),dims)
    vort1d_sm      = ndtooned(vort_smth)
    inds_sm        = ind_resolve(maxind(vort1d_sm),dims)

   ;==================================================================
   ; Either calculate the centre, or read in pre-calculated values
   ;==================================================================

   if (calc .eq. 0) then                        ; Find centre on each model level

   ; Find latitude and longitude of storm centre using vorticity calculation above
     lat_max       = t0 + (dy0 * inds(0,0))
     lon_max       = n0 + (dx0 * inds(0,1))

     print_clock("it = "+it)
     print_clock("lat_max = "+centre(cn,it,0)+" ; lon_max = "+centre(cn,it,1))

   ; Difference between 'lat_max'/'lon_max' and values read in using NCL function
     dlat          = abs(lat_max - centre(cn,it,0))
     dlon          = abs(lon_max - centre(cn,it,1))

   ; If our calculations of latitude and longitude deviate too much from the values
   ; we read in using the NCL function ('st_rm'), then calculate a new centre
   ; using the minimum pressure. Constrain the new calculation using the values
   ; from 'st_rm', to prevent spurious centres being identified and analysed.

     if (dlat .gt. 0.75 .or. dlon .gt. 0.75) then

    ; 'it' because 'pc' stream data is our base, and 'centre' has 1 element fewer
      lt1 = centre(cn,it,0) - 0.5
      lt2 = centre(cn,it,0) + 0.5
      ln1 = centre(cn,it,1) - 0.5
      ln2 = centre(cn,it,1) + 0.5

      prs_plane0 = prs_plane({lt1:lt2},{ln1:ln2})       ; Smaller grid (1.0 x 1.0 deg)
      prs_min0   = min(prs_plane0)                      ; Find minimum on smaller grid
      dims_p0    = dimsizes(prs_plane0)                 ; Size of smaller grid
      prs_1d0    = ndtooned(prs_plane0)                 ; Create 1-D array
      inds_p0    = ind_resolve(minind(prs_1d0),dims_p0) ; Find index of minimum

      lat_max    = lt1 + (dy0 * inds_p0(0,0))
      lon_max    = ln1 + (dx0 * inds_p0(0,1))
      print("New centre: "+lat_max+" degrees N, "+lon_max+" degrees E")

     end if

    elseif (calc .eq. 1) then				; Find centre using offline script

     lat_max     = centre(cn,it,0)
     lon_max     = centre(cn,it,1)

    end if 

;==============================================================
; Create arrays for later switch to cylindrical coordinates
;==============================================================

  ;===============
  ; 'pc' stream
  ;===============

  ; Indices of lat/lon grid for 'u' and 'v' in 'pc' stream
    latC_plane  = ind(t0 .le. latC .and. latC .le. t1)
    lonC_plane  = ind(n0 .le. lonC .and. lonC .le. n1)

  ; Array sizes
    latC_size   = dimsizes(latC_plane)
    lonC_size   = dimsizes(lonC_plane)

  ; Array of latitude/longitude points corresponding to the indices above
    latC_sub    = latC( latC_plane(0):latC_plane(d0-1) )
    lonC_sub    = lonC( lonC_plane(0):lonC_plane(d1-1) )

  ;===============
  ; 'pd' stream
  ;===============

  ; Indices of lat/lon grid for 'str_xz' in 'pd' stream
    latD_plane  = ind(t0 .le. latD .and. latD .le. t1)
    lonD_plane  = ind(n0 .le. lonD .and. lonD .le. n1)

  ; Array sizes
    latD_size   = dimsizes(latD_plane)
    lonD_size   = dimsizes(lonD_plane)

  ; Array of latitude/longitude points corresponding to the indices above
    latD_sub    = latD( latD_plane(0):latD_plane(d0-1) )
    lonD_sub    = lonD( lonD_plane(0):lonD_plane(d1-1) )

    delete([/vort1d,inds,vort1d_sm,vort_smth,inds_sm/])

;===========================================
; Switch to cylindrical coordinate system
;===========================================

  ; Create new arrays describing cylindrical coordinate system (2D for now)
    xpos = new((/dimsizes(radii),dimsizes(thetas)/),float)
    ypos = new((/dimsizes(radii),dimsizes(thetas)/),float)

  ; New arrays containing variables we want to plot
    pv_int   = new((/dimsizes(radii),dimsizes(thetas)/),float)  ; Potential vorticity (PV)
    pv_sd0   = new((/dimsizes(radii)/),float) 			; PV standard deviation (SD)

    vort_int = new((/dimsizes(radii),dimsizes(thetas)/),float)  ; Relative vorticity 
    vort_sd0 = new((/dimsizes(radii)/),float)                   ; Relative vorticity SD

  ; Add variable metadata
    copy_VarAtts(pv_plane,pv_int)
    pv_int@description = "Potential vorticity"
    pv_int@units       = "10~S~6~N~ K m~S~2~N~ s~S~-1~N~ kg~S~-1~N~"

    copy_VarAtts(pv_plane,pv_sd0)
    pv_sd0@description = "PV standard deviation"
    pv_sd0@units       = ""

    copy_VarAtts(vort_plane,vort_int)
    vort_int@description = "Relative vorticity"
    vort_int@units       = "s~S~-1~N~"

    copy_VarAtts(vort_plane,vort_sd0)
    vort_sd0@description = "Relative vorticity standard deviation"
    vort_sd0@units       = ""

    xpos!0    	 = "rad"
    xpos!1    	 = "azi"
    xpos&rad  	 = radii
    xpos&azi  	 = thetas

    ypos!0    	 = "rad"
    ypos!1    	 = "azi"
    ypos&rad  	 = radii
    ypos&azi     = thetas

    pv_int!0     = "rad"
    pv_int!1     = "azi"
    pv_int&rad   = radii
    pv_int&azi   = thetas

    pv_sd0!0     = "rad"
    pv_sd0&rad   = radii

    vort_int!0   = "rad"
    vort_int!1   = "azi"
    vort_int&rad = radii
    vort_int&azi = thetas

    vort_sd0!0   = "rad"
    vort_sd0&rad = radii

  ; Latitude/longitude (small grid)
    do irad = 0, dimsizes(radii)-1
      r = radii(irad)
      do iang = 0, dimsizes(thetas)-1
        theta = thetas_rad(iang)
          xpos(irad,iang) = lon_max + r*cos(theta)
          ypos(irad,iang) = lat_max + r*sin(theta)
      end do
    end do

  ; Loop over azimuth angles to fill entire circular array of data points
    do iang = 0, dimsizes(thetas_rad)-1
    ; Either read in potential vorticity
      if (vortO.eq. 0) then 
       pv_int(:,iang)   = linint2_points(lonD_sub,latD_sub,pv_plane,False,\
                                         xpos(:,iang),ypos(:,iang),0)
    ; Or relative vorticity
      elseif (vortO .eq. 1) then 
       vort_int(:,iang) = linint2_points(lonC_sub,latC_sub,vort_plane,False,\
                                         xpos(:,iang),ypos(:,iang),0)
      end if 
    end do

  ;=========================================================================
  ; Azimuthally average around the storm and calculate standard deviation
  ;=========================================================================

    do irad = 0, dimsizes(radii)-1

      r = radii(irad)
 
    ; Potential vorticity
      if (vortO .eq. 0) then 

       pv_sd0(irad)        = dim_stddev_n(pv_int(irad,:),0)   ; Standard deviation
       pv_azi(irad)        = avg(pv_int(irad,:))              ; Azimuthal average

    ; Relative vorticity
      elseif (vortO .eq. 1) then 

       vort_sd0(irad)      = dim_stddev_n(vort_int(irad,:),0) ; Standard deviation
       vort_azi(irad)      = avg(vort_int(irad,:))            ; Azimuthal average

      end if 

    end do

  ; Input into larger arrays

  ; Potential vorticity
    if (vortO .eq. 0) then 

     pv_plot(:,it)    = pv_azi(:)	; Input into larger arrays
     pv_sd1(:,it)     = pv_sd0(:)	; Input single value from each radial circle

  ; Relative vorticity
    elseif (vortO .eq. 1) then 

     vort_plot(:,it)  = vort_azi(:)     ; Input into larger arrays
     vort_sd1(:,it)   = vort_sd0(:)     ; Input single value from each radial circle

    end if 

  ; Tidy up
    delete([/pv_int, pv_sd0, vort_int, vort_sd0/])
    delete([/u_plane, v_plane, pv_plane, vort_plane, prs_plane/])
    delete([/latD_plane,lonD_plane,latC_plane,lonC_plane/])
    delete([/latD_sub,lonD_sub,latC_sub,lonC_sub/])

  ; Tidy up (avoid dimension size errors)
    delete([/lonC,latC,u,v,vort,u0,v0,vort0/])
    delete([/lonD,latD,p,z,pv/])

  end do     ; End time loop (do it = ts, tf)

;==================================
; Write out values to text file 
;==================================

; Output file directory path (layer-averaged or single layer)
  diri_out = "$ar/text/ml5_"+dat+"_"+ens0+"_"+cn0+"_"

; Calculate the maximum standard deviation

; Potential vorticity
  if (vortO .eq. 0) then

   pv_sd      = new(numINPUT_c,float)	           ; Array to hold values below 
   pv_sd(:)   = dim_max_n(pv_sd1(:,:),0)           ; 1D array [numINPUT_c]

   if (lay .eq. 1) then 
    print("Max. standard deviation of layer-averaged ("+\
          hy0+"-"+hy1+" m) PV = "+pv_sd)
   else
    print("Max. standard deviation of "+hy+" m PV = "+pv_sd)
   end if 

; Relative vorticity
  elseif (vortO .eq. 1) then

   vort_sd    = new(numINPUT_c,float)                ; Array to hold values below
   vort_sd(:) = dim_max_n(vort_sd1(:,:),0)         ; 1D array [numINPUT_c]

   if (lay .eq. 1) then
    print("Max. standard deviation of layer-averaged ("+\
          hy0+"-"+hy1+" m) relative vorticity = "+vort_sd)
   else
    print("Max. standard deviation of "+hy+" m relative vorticity = "+vort_sd)  
   end if 

  end if 

; Read in the maximum tangential wind from text file (calculated on model levels)
  vtan_max    = asciiread("$ar/text/ml5_"+dat+"_"+ens0+"_vtan_max_"+cn0+".txt",\
			  (/numINPUT_c/),"float")	

; Calculate the tendency of the maximum tangential wind
  do it = 0, numINPUT_c-1

   tp1 = min((/it+1,numINPUT_c-1/))
   tm1 = max((/it-1,0/))  

   if (tend .eq. 2) then 
    vtan_tend(it) = (vtan_max(tp1) - vtan_max(tm1)) / 2
   elseif (tend .eq. 1) then 
    vtan_tend(it) = (vtan_max(it) - vtan_max(tm1)) / 1
   end if 
  
  end do 

; Finally, write to text files

; Relative vorticity
  if (vortO .eq. 1) then 

   if (lay .eq. 1) then 
    vort_out = diri_out+"vort_lay_"+hy0+"m_"+hy1+"m.txt"
    asciiwrite(vort_out,vort_sd)
   elseif (lay .eq. 0) then 
    vort_out = diri_out+"vort_"+hy+"m.txt"
    asciiwrite(vort_out,vort_sd)
   end if 

; Potential vorticity
  elseif (vortO .eq. 0) then 

   if (lay .eq. 1) then 
    pv_out = diri_out+"pv_lay_"+hy0+"m_"+hy1+"m.txt"
    asciiwrite(pv_out,pv_sd)
   elseif (lay .eq. 0) then 
    pv_out = diri_out+"pv_"+hy+"m.txt"
    asciiwrite(pv_out,pv_sd)
   end if 

  end if 

; Write out tangential wind tendency (calculated from text file)
  ten_out = diri_out+hy+"m_tend.txt"
  asciiwrite(ten_out,vtan_tend)

end 