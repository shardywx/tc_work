; Produce a time series of vertical mass flux 

; Run using:

; ncl 'opt="x11"' 'dat="02T12"' 'ens0="em11"' dist=3.0 ts0=24 tf0=96 calc=1 write=0 
; vn=20 vx=80 fn=2.0 fx=14.0 multi=1 lev0=16 lev1=34 lev2=18 lev3=14 lev4=12 lev5=10 
; vtan=1 flux=1 ar=1.0 nr=21 cn0=\"slp\" tc_mflux_time_series.ncl

; 'opt'    = output file format ("pdf" or "x11")
; 'dat'	   = initialisation time string ("02T12", "03T00", ...)
; 'ens0'   = ensemble member (em00 -- em11)
; 'dist'   = size of box following storm (degrees)
; 'ts0'    = starting time index (0–119)
; 'tf0'    = ending time index (0–119)
; 'calc'   = centre calculated on each level (0) or equal to 950 hPa centre (1)
; 'write'  = write out mass flux values to text file (1) or leave alone (0)
; 'vn'	   = minimum tangential wind on final plot (m/s)
; 'vx'	   = maximum tangential wind on final plot (m/s)
; 'fn'	   = minimum mass flux on final plot
; 'fx'	   = maximum mass flux on final plot 
; 'vtan'   = overlay max vtan on model (2) or pressure (1) levels, or 10-m windspeed (0)
; 'multi'  = plot single (0) or multiple (1) vertical levels
; 'lev0'   = lower model level for mass flux plot (0-58)
; 'lev1'   = upper model level for mass flux plot (0-58)
; 'lev2'   = additional model level (upper or lower)
; 'lev3'   = additional model level (upper or lower)
; 'lev4'   = additional model level (upper or lower)
; 'lev5'   = additional model level (upper or lower)
; 'flux'   = plot mass flux (1) or vertical velocity (0)
; 'ar'	   = distance in degrees of outer radius (1.0, 2.0, etc)
; 'nr'	   = number of radial circles between r = 0 and r = 'ar' (21, 41, 61, etc)
; 'cn0'    = offline storm centre: "slp", "vort", "geo", "geo_sm"

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

  if (ts0 .lt. 3) then  
   print("Choose a starting index greater than 2! Quitting script...")
   exit()
  end if 

; Create string from user input ('ar'): '1.0' --> 100km, '1.5' --> 150km, etc
  dstr0 = ar * 100
  dstr = sprintf("%3.0fkm",dstr0)	; 3 characters, 0 decimal points, taking float as input

; Find 'pc' and 'pd' files (model height level data)
  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0
  fili_c1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pc.nc")
  fili_c      = diri+"/"+fili_c1
  fili_d1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pd.nc")
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

;==================================================
; Get model level heights from 'pd' stream file
;==================================================

  d      = addfile(fili_d,"r")          ; Read in 'pd' stream [w,p,t,theta,kmh,tau]	
  hybD   = d->hybrid_ht(:)              ; 63 model ('rho') levels

  hy = sprintf("%0.0f",hybD)
  lev_m0 = hy(lev0)
  lev_m1 = hy(lev1)

  if (multi .eq. 1) then 
   lev_m2 = hy(lev2)
   lev_m3 = hy(lev3)
   lev_m4 = hy(lev4)
   lev_m5 = hy(lev5)
  end if 

  delete([/d,hybD,hy/])

;======================================================================
; If text file doesn't already exist, then perform steps below...
; If it does exist, then skip to the plotting section of the script
;======================================================================

; File path for text file containing mass flux data
  fili_txt  = "$ar/text/mflux_"+dstr+"_"+dat+"_"+ens0+"_"

; File extensions for chosen model levels ('m0' and 'm1')
  fili_txt0 = fili_txt+lev_m0+"m_"+cn0+".txt"
  fili_txt1 = fili_txt+lev_m1+"m_"+cn0+".txt"

; Alternative file paths for text files containing vertical velocity
  fili_txt3 = "$ar/text/w_"+dstr+"_"+dat+"_"+ens0+"_"+lev_m0+"m_"+cn0+".txt"
  fili_txt4 = "$ar/text/w_"+dstr+"_"+dat+"_"+ens0+"_"+lev_m1+"m_"+cn0+".txt"

  if (vtan .eq. 2) then
   fili_txt2 = "$ar/text/mlev_"+dat+"_"+ens0+"_vtan_max_calc1.txt"
  elseif (vtan .eq. 1) then 
   fili_txt2 = "$ar/text/"+dat+"_"+ens0+"_vtan_max.txt"
  elseif (vtan .eq. 0) then 
   fili_txt2 = "$ar/text/cp_"+dat+"_"+ens0+"_10m.txt" 
  end if 

  if (flux .eq. 1) then 
   txt_file0 = fili_txt0
   txt_file1 = fili_txt1
  elseif (flux .eq. 0) then 
   txt_file0 = fili_txt3
   txt_file1 = fili_txt4
  end if 
 
  if ( .not. fileexists (txt_file0) .or. .not. fileexists (txt_file1) ) then
 
  print("Text file doesn't exist - calculating mass flux...")

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

; ar = distance (degrees) of outer radial circle (1.0)
; nr = total number of radial circles (21) 

  thetas     = new(73,float)
  pii        = 3.14159265
  radii      = fspan(0.0,ar,nr)      ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,73)   ; Array of azimuth angles (0, 5, 10, etc)
  thetas_rad = (pii/180)*thetas      ; Degrees to radians
  rad_size   = dimsizes(radii)-1     ; Index for later calculation
  num_l      = 59      	   	     ; Number of model levels (edit to evolve smoothly)

; Define constants
  rd	= 287.0	  ; Gas constant for dry air (J/kg/K^2)
  cp0	= 1004.0  ; Specific heat of dry air at constant pressure (J/kg/K)
  kp	= 0.286   ; For Exner function calculations (0.286)
  pref	= 1000.0  ; Reference pressure p0 (1000 hPa)

;===========================================================
; Array for distances between points in cylindrical grid 
;===========================================================

  rad_arr0  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)

;========================================
; Define 3D arrays for mean quantities 
;========================================

; Mean mass flux {flux}
  mean_flux0       = new((/num_l,dimsizes(radii),numTIMES/),float)
  mean_flux0!0     = "lev"
  mean_flux0!1     = "rad"
  mean_flux0!2     = "time"
  mean_flux0@description = "Azimuthally averaged mass flux"
  mean_flux0@units = "kg m~S~2 s~S~-1~N~"

; Mean vertical velocity {vvel}
  mean_w0          = new((/num_l,dimsizes(radii),numTIMES/),float)
  mean_w0!0        = "lev"
  mean_w0!1        = "rad"
  mean_w0!2        = "time"
  mean_w0@description = "Azimuthally averaged vertical velocity"
  mean_w0@units    = "m s~S~-1~N~"

; Mean pressure {p}
  mean_pres0       = new((/num_l,dimsizes(radii),numTIMES/),float)
  mean_pres0!0     = "lev"
  mean_pres0!1     = "rad"
  mean_pres0!2     = "time"
  mean_pres0@description = "Azimuthally averaged pressure"
  mean_pres0@units = "hPa"

; Create time array for later plotting 
  time_int = new(120,"integer")
  time_int = ispan(0,119,1)
  time_int!0 = "time"

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
   ct0	       = 0			    ; Counter variable (independent of position)
   llbox       = toint(dist*50)	    ; Size of lat/lon arrays (see line below)
   dsize       = (/llbox,llbox/)           ; Size of lat/lon arrays (depends on 'dist')
   d0	       = dsize(0)
   d1	       = dsize(1)

   setvalues NhlGetWorkspaceObjectId
     "wsMaximumSize" : 1000000000
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
    title_arr(it) = "Valid at "+time_str(it)+" (T+"+it+")"

  ; Index for 'pd' stream files (relative to 'pc')
    dt = it - 2	    	

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
    hybC   = c->hybrid_ht(:) 		          ; 63 model ('theta') levels

    u_varname = "x-wind"
    v_varname = "y-wind"

  ; Read in horizontal wind components 
  ; Both variables below -- [hybC | 63] * [latC | 150] * [lonC | 150]
    u0  = c->$u_varname$(it+1,:,{t0:t1},{n0:n1}) ; Zonal wind	
    v0  = c->$v_varname$(it+1,:,{t0:t1},{n0:n1}) ; Meridional wind 

  ; Calculate relative vorticity using centered finite differences 
  ; Option '2' --> boundary points estimated using one-sided difference scheme
    vort = uv2vr_cfd(u0, v0, latC, lonC, 2)
    copy_VarCoords(u0, vort)
    vort@units = "s~S~-1~N~"
    vort@name = "Relative vorticity on model levels"

  ; Calculate storm-relative horizontal wind
    u = u0 - u_cyc(cn,it)
    v = v0 - v_cyc(cn,it)

  ; Add metadata from original horizontal wind arrays
    copy_VarCoords(u0,u)
    copy_VarAtts(u0,u)

    copy_VarCoords(v0,v)
    copy_VarAtts(v0,v)

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

  ; Read in temperature, pressure, theta, vertical velocity and geopotential height
  ; All variables below -- [hybD | 63] * [latD | 150] * [lon1D | 151]
    t	   = d->temp(it-1,:,{t0:t1},{n0:n1})          ; Temperature (K)
    p      = d->p(it-1,:,{t0:t1},{n0:n1})             ; Pressure (Pa)
    w      = d->dz_dt(it-1,:,{t0:t1},{n0:n1})	    ; Vertical velocity (m s-1)
    z      = d->ht(it-1,:,{t0:t1},{n0:n1})    	    ; Geopotential height (m)

  ; Calculate dry air density from pressure and temperature
    rho	   = p / (rd * t)

  ; Add metadata
    copy_VarCoords(p, rho)
    rho@description = "Dry air density"
    rho@units = "kg m~S~-3~N~"

  ; Calculate mass flux 
    mflux = rho * w
    printVarSummary(rho)
    printVarSummary(w)
    copy_VarCoords(w,mflux)
    copy_VarAtts(w,mflux)
    mflux@units = "kg m~S~-2~N~ s~S~-1~N~"
    mflux@standard_name = "vertical_mass_flux"
    mflux@long_name     = "Vertical mass flux"
    mflux@title         = "Vertical mass flux"
    mflux@name          = "vert_mflux"

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

    do ilev = 0, num_l-1 ; Loop over all model levels

     hy = sprintf("%0.0f",hybD(ilev))
     print("Working on model level " + ilev+ " ("+hy+" m AGL)")

   ; Read in variables on model levels
     prs_plane              = p(ilev,:,:)	   ; Pressure
     vort_plane	   	    = vort(ilev,:,:)	   ; Relative vorticity
     flux_plane0            = mflux(ilev,:,:)	   ; Mass flux 
     w_plane0		    = w(ilev,:,:)	   ; Vertical velocity

   ; Remove extra array elements (keep all arrays same size)
     fsize = dimsizes(flux_plane0)
     wsize = dimsizes(w_plane0)

     if (fsize(0) .gt. d0 .or. fsize(1) .gt. d1) then 
      flux_plane = flux_plane0(0:d0-1,0:d1-1)
     else
      flux_plane = flux_plane0
     end if 

     if (wsize(0) .gt. d0 .or. wsize(1) .gt. d1) then
      w_plane = w_plane0(0:d0-1,0:d1-1)
     else
      w_plane = w_plane0
     end if

     vort_max      	    = max(vort_plane)
     vort_smth     	    = smth9_Wrap(vort_plane, 0.5, 0.5, True)
     vort_max_smth 	    = max(vort_smth)

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
     dlat 	   	    = abs(lat_max - centre(cn,it,0))
     dlon	   	    = abs(lon_max - centre(cn,it,1))

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

      lat_max 	 = lt1 + (dy * inds_p0(0,0))
      lon_max 	 = ln1 + (dx * inds_p0(0,1))
      print("New centre: "+lat_max+" degrees N, "+lon_max+" degrees E")

     end if 

   elseif (calc .eq. 1) then

      lat_max    = centre(cn,it,0)
      lon_max    = centre(cn,it,1)

   end if

   ;======================================================================
   ; Create lat/lon arrays for later switch to cylindrical coordinates
   ;======================================================================

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
     latD_sub0   = latD(latD_plane(0):latD_plane(latD_size-1))
     lonD_sub0   = lonD(lonD_plane(0):lonD_plane(lonD_size-1))

   ; Remove extra array	elements (keep all arrays same size)
     lsize0 = dimsizes(latD_sub0)
     lsize1 = dimsizes(lonD_sub0)

     if (lsize0 .gt. d0) then 
      latD_sub = latD_sub0(0:d0-1)
     else 
      latD_sub = latD_sub0
     end if 

     if (lsize1 .gt. d1) then
      lonD_sub = lonD_sub0(0:d0-1)
     else
      lonD_sub = lonD_sub0
     end if
     delete([/lsize0, lsize1, latD_sub0, lonD_sub0/])

   ; Indices of lat/lon grid for 'str_yz' in 'pd' stream
     lat1D_plane = ind(t0.le.lat1D.and.lat1D.le.t1)
     lon1D_plane = ind(n0.le.lon1D.and.lon1D.le.n1)

   ; Array sizes
     lat1D_size  = dimsizes(lat1D_plane)
     lon1D_size  = dimsizes(lon1D_plane)

   ; Array of latitude/longitude points corresponding to the indices above
     lat1D_sub0  = lat1D(lat1D_plane(0):lat1D_plane(lat1D_size-1))
     lon1D_sub0  = lon1D(lon1D_plane(0):lon1D_plane(lon1D_size-1))

   ; Remove extra array elements (keep all arrays same size)
     lsize0 = dimsizes(lat1D_sub0)
     lsize1 = dimsizes(lon1D_sub0)

     if (lsize0 .gt. d0) then
      lat1D_sub = lat1D_sub0(0:d0-1)
     else
      lat1D_sub = lat1D_sub0
     end if

     if (lsize1 .gt. d1) then
      lon1D_sub = lon1D_sub0(0:d0-1)
     else
      lon1D_sub = lon1D_sub0
     end if
     delete([/lsize0, lsize1, lon1D_sub0, lat1D_sub0/])

     delete([/vort1d,inds,vort1d_sm,vort_smth,inds_sm/])

;===========================================
; Switch to cylindrical coordinate system
;===========================================

   ; Create new arrays describing cylindrical coordinate system (2D for now)
     xpos      = new((/dimsizes(radii),dimsizes(thetas)/),float)
     ypos      = new((/dimsizes(radii),dimsizes(thetas)/),float)

   ; New arrays containing regular variables we want to plot
     flux_int0 = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Mass flux 
     w_int0    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vertical velocity

;=========================
; Add variable metadata 
;=========================

     copy_VarAtts(flux_plane,flux_int0)
     flux_int0@description = "Mass flux"
     flux_int0@units       = "kg m~S~2~N~ s~S~-1~N~"

     copy_VarAtts(w_plane,w_int0)
     w_int0@description = "Vertical velocity"
     w_int0@units       = "m s~S~-1~N~"

   ; Create cylindrincal coordinate array from original latitude/longitude grid
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
   ; flux_plane  -- latD,  lon1D
   ; w_plane	 -- latD,  lon1D

     do iang = 0, dimsizes(thetas_rad)-1
      flux_int0(:,iang) = linint2_points(lon1D_sub,latD_sub,flux_plane,False,\
                                         xpos(:,iang),ypos(:,iang),0)
      w_int0(:,iang)    = linint2_points(lon1D_sub,latD_sub,w_plane,False,\
                                         xpos(:,iang),ypos(:,iang),0)
     end do

     xpos!0       = "rad"
     xpos!1       = "azi"
     xpos&rad     = radii
     xpos&azi     = thetas

     ypos!0       = "rad"
     ypos!1       = "azi"
     ypos&rad     = radii
     ypos&azi     = thetas

    ;==========================================================
    ; Ignore all regions of negative mass flux (set to zero)
    ;==========================================================
    
     flux_int     = where( flux_int0 .gt. 0, flux_int0, 0 )
     w_int	  = where( w_int0 .gt. 0, flux_int0, 0)       ; ERROR IN THIS LINE???

    ;===================================================
    ; Calculate mean (azimuthally-averaged) variables 
    ;===================================================

      do irad = 0, dimsizes(radii)-1

	mean_flux0(ilev,irad,ct)    = avg(flux_int(irad,:))   ; [ {mflux} ]
	mean_w0(ilev,irad,ct)	    = avg(w_int(irad,:))      ; [ {vvel} ]

      end do

    ;======================================================================
    ; Calculate eddy terms involving radial derivatives [du'/dr, dp'/dr]
    ;======================================================================

       do irad = 1, dimsizes(radii)-1

        rp1 = min((/irad+1,dimsizes(radii)-1/))
        rm1 = max((/irad-1,0/))

        do th0 = 0, dimsizes(thetas)-1

         tp1 = min((/th0+1,dimsizes(thetas)-1/))
         tm1 = max((/th0-1,0/))

       ; Calculate the radial distance ('dr') between all points in array
       ; Arguments 1 and 2 represent the lat/lon of the first point, respectively 
       ; Arguments 3 and 4 represent the lat/lon of the second point, respectively 
       ; Argument 5 ('2') tells the script to avoid extra interpolated points
       ; Argument 6 ('3') returns the distance in metres

       ; Why were the latitude and longitude switched in the azimuthal derivative loop?
       ; Make sure the values I'm outputting here are correct, before continuing. 

       ; Radial distance between each point and its nearest neighbour
         dr  = gc_latlon(ypos(rm1,th0), xpos(rm1,th0), \
                         ypos(irad,th0), xpos(irad,th0), 2, 3)

         print("rp1 = "+rp1+" ; rm1 = "+rm1+" ; irad = "+irad)

       ; Input radial distances into large array for later calculation
       	 rad_arr0(ilev,irad,th0,ct) = dr 

        end do                        ; End loop over azimuth angles

       end do                         ; End loop over radii

    ; Tidy up
      delete([/flux_int,flux_int0,flux_plane,flux_plane0,vort_plane,prs_plane,w_plane0/])
      delete([/xpos,ypos/])
      delete([/latD_plane,lonD_plane,lon1D_plane,lat1D_plane/])
      delete([/latD_sub,lonD_sub,lat1D_sub,lon1D_sub/])

     end do			    ; End pressure level loop (do ilev = 0, num_l-1)

    ct  = ct + 1		    ; Counter variable (time)
    ct0 = ct0 + 1		    ; Independent counter variable (time)

  ; Tidy up before next iteration (all variables from 'pc' and 'pd' streams)
    delete([/lonC,latC,lon1C,lat1C,u,v,vort/])
    delete([/lonD,latD,lon1D,lat1D,t,p,w,z,rho,mflux/])

   end do     ; End time loop (do it = 0, times-1)

;===========================================
; Calculate radially integrated mass flux
;===========================================

   print(rad_arr0(0,:,:,0) )
   exit()
 
     ; Calculate average along azimuthal dimension of 'rad_arr0'
     ; Resulting array ('rad_arr') contains the azimuthally-averaged
     ; radial distance between one radial circle and the next
       rad_arr   = dim_avg_n_Wrap( rad_arr0(:,:,:,:), 2)  

     ; Calculate radial integral of mass flux
       flux_sum  = mean_flux0(:,:,:) * rad_arr(:,:,:)    ; Calculate [flux * dr]
       mean_flux = dim_sum_n_Wrap( flux_sum, 1)		 ; Integrate radially [sum along dim #1]

     ; Calculate radial integral of vertical velocity
       w_sum	 = mean_w0(:,:,:) * rad_arr(:,:,:)	 ; Calculate [w * dr]
       mean_w	 = dim_sum_n_Wrap( w_sum, 1)		 ; Integrate radially [sum along dim #1]

;=================================================
; Tidy up and produce final arrays for plotting 
;=================================================

     ; Add metadata
       hyb_plot       = hybC / 1000
       hyb_plot@units = "km"
 
       mean_flux!0    = "lev"
       mean_flux&lev  = hyb_plot
       mean_flux!1    = "time"
       mean_flux&time = time_int

       mean_w!0       = "lev"
       mean_w&lev     = hyb_plot
       mean_w!1       = "time"
       mean_w&time    = time_int

     ; Convert to usable units for ploting 
       mean_flux      = mean_flux / 10000
       mean_w	      = mean_w / 10000

     ; Extract data on chosen model levels 
       plot_flux0     = mean_flux(lev0,:)	 ; Default is 1.5 km (16)
       plot_flux1     = mean_flux(lev1,:)	 ; Default is 6 km (34)

       plot_w0	      = mean_w(lev0,:)		 ; Lower level
       plot_w1	      = mean_w(lev1,:)		 ; Upper level

     ; Combine into single array
       plot_flux      = (/plot_flux0, plot_flux1/)
       plot_w	      = (/plot_w0, plot_w1/)

     ; Write out data to text file for faster plotting 
       if (write .eq. 1) then 

        if (flux .eq. 0) then 
         w0_out     = fili_txt3
 	 asciiwrite(w0_out, plot_w0)
 	 w1_out     = fili_txt4 
	 asciiwrite(w1_out, plot_w1)
        elseif (flux .eq. 1) then 
         flux0_out     = fili_txt0
         asciiwrite(flux0_out, plot_flux0)
	 flux1_out     = fili_txt1
         asciiwrite(flux1_out, plot_flux1)
	end if 

       end if 

       exit()

 elseif ( fileexists (fili_txt0) ) then
  
  print("Text file already exists - reading in mass flux...")

; Create time array for plotting
  time_int = new(120,"integer")
  time_int = ispan(0,119,1)
  time_int!0 = "time"

; Read in mass flux data from text files 
  plot_flux0 = asciiread(fili_txt0,(/numTIMES/),"float")	; Lower model level	
  plot_flux1 = asciiread(fili_txt1,(/numTIMES/),"float")	; Higher model level
  plot_flux  = (/plot_flux0, plot_flux1/)

; Read in extra mass flux data for multi-line plots 
  if (multi .eq. 1) then 
   plot_flux2 = asciiread(fili_txt+lev_m2+"m.txt",(/numTIMES/),"float")
   plot_flux3 = asciiread(fili_txt+lev_m3+"m.txt",(/numTIMES/),"float")
   plot_flux4 = asciiread(fili_txt+lev_m4+"m.txt",(/numTIMES/),"float")
   plot_flux5 = asciiread(fili_txt+lev_m5+"m.txt",(/numTIMES/),"float")
  end if 

; Or, read in vertical velocity data from text files 
  if (flux .eq. 0) then 
   plot_w0    = asciiread(fili_txt3,(/numTIMES/),"float")        ; Lower model level
   plot_w1    = asciiread(fili_txt4,(/numTIMES/),"float")        ; Higher model level
   plot_w     = (/plot_w0, plot_w1/)
  end if 

; Also read in max tangential wind data from text file
  plot_vtan  = asciiread(fili_txt2,(/numTIMES/),"float")

 end if		; End of IF statement linked to existence of text file 
       
;============================================
; Define output file path and plot details 
;============================================

  ; Output file location and type
    if (flux .eq. 1) then  
   
     if (vtan .eq. 2) then 
      output = "$sam/nepartak/images/cyl_coords/mflux_vtanM_xz_"+dstr+"_"+dat+"_"+ens0+\
      	       "_"+lev_m0+"m_"+lev_m1+"m_"+cn0
     elseif (vtan .eq. 1) then
      output = "$sam/nepartak/images/cyl_coords/mflux_vtanP_xz_"+dstr+"_"+dat+"_"+ens0+\
               "_"+lev_m0+"m_"+lev_m1+"m_"+cn0
     elseif (vtan .eq. 0) then 
      output = "$sam/nepartak/images/cyl_coords/mflux_10m_xz_"+dstr+"_"+dat+"_"+ens0+\
      	       "_"+lev_m0+"m_"+lev_m1+"m_"+cn0
     end if 

    elseif (flux .eq. 0) then 

     if (vtan .eq. 2) then
      output = "$sam/nepartak/images/cyl_coords/vvel_vtanM_xz_"+dstr+"_"+dat+"_"+ens0+\
               "_"+lev_m0+"m_"+lev_m1+"m_"+cn0
     elseif (vtan .eq. 1) then
      output = "$sam/nepartak/images/cyl_coords/vvel_vtanP_xz_"+dstr+"_"+dat+"_"+ens0+\
               "_"+lev_m0+"m_"+lev_m1+"m_"+cn0
     elseif (vtan .eq. 0) then
      output = "$sam/nepartak/images/cyl_coords/vvel_10m_xz_"+dstr+"_"+dat+"_"+ens0+\
               "_"+lev_m0+"m_"+lev_m1+"m_"+cn0
     end if 

    end if 
    wks = gsn_open_wks(opt,output)

  ; Set up panel plot
    panel = new(1,graphic)

;==========================
; Options for plotting 
;==========================

  ; Mass flux/vertical velocity at 1.5 km
    opts_mflux0                            = True
    opts_mflux0@gsnMaximize		   = True

    opts_mflux0@vpWidthF             	   = 0.8
    opts_mflux0@vpHeightF            	   = 0.5
    opts_mflux0@tmYROn               	   = False		  ; Turn off tick marks
    opts_mflux0@tmXTOn               	   = False         	  ; Turn off tick marks
    opts_mflux0@tmXBLabelFontHeightF 	   = 0.01		  ; x-axis label font height
    opts_mflux0@trYMaxF              	   = fx 		  ; Set y-axis maximum
    opts_mflux0@trYMinF              	   = fn                   ; Set y-axis minimum

    opts_mflux0@trXMinF                    = ts0          	  ; Set x-axis minimum
    opts_mflux0@trXMaxF                    = tf0 		  ; Set x-axis maximum

    opts_mflux0@tmXBMode             	   = "Explicit"
    opts_mflux0@tmXBValues		   = ispan(0,120,12)
    opts_mflux0@tmXBLabels                 = (/"T+0","T+12","T+24","T+36","T+48","T+60",\
    					       "T+72","T+84","T+96","T+108","T+120"/)
    opts_mflux0@tiMainString         	   = "Mass flux time series: simulation "+ens0 
    opts_mflux0@tiYAxisString              = "Mass flux / kg s~S~-1~N~"
    opts_mflux0@tiXAxisString          	   = "Time"
    opts_mflux0@tmXBLabelFontHeightF       = 0.0150    ; x-axis	label font size
    opts_mflux0@tiXAxisFontHeightF	   = 0.0200    ; x-axis font size
    opts_mflux0@tiYAxisFontHeightF     	   = 0.0200    ; y-axis font size
    opts_mflux0@tiMainFontHeightF	   = 0.0225    ; title font size

    opts_mflux0@tmXMajorGridLineDashPattern  = 2    ; Dashed gridlines
    opts_mflux0@tmYMajorGridLineDashPattern  = 2
    opts_mflux0@gsnDraw                      = False
    opts_mflux0@gsnFrame                     = False
    opts_mflux0@gsnPaperOrientation          = "Portrait"

    opts_mflux0@xyMarkLineMode               = "MarkLines"
    opts_mflux0@xyLineColors                 = "red3"
    opts_mflux0@xyLineThicknesses            = 2.0
    opts_mflux0@xyMarkers                    = 16
    opts_mflux0@xyMarkerColors               = "red3"
    opts_mflux0@xyMarkerSizes                = 5.0
    opts_mflux0@xyMarkerThicknesses          = 2.0
    opts_mflux0@xyDashPatterns		     = 0

  ; Options for plotting tangential wind
    opts_vtan                                = True
    opts_vtan@tmXBLabelFontHeightF           = 0.01
    
    opts_vtan@trYMaxF                        = vx    ; Set maximum windspeed
    opts_vtan@trYMinF                        = vn    ; Set minimum windspeed

    opts_vtan@tiMainString                   = ""
    opts_vtan@tiYAxisString                  = "Maximum tangential windspeed (m s~S~-1~N~)"
    opts_vtan@tiYAxisFontHeightF             = 0.0125

    opts_vtan@xyLineColors                   = "royalblue1"
    opts_vtan@xyLineThicknesses              = 5.0
    opts_vtan@xyDashPatterns                 = 0.0
    opts_vtan@gsnDraw			     = False
    opts_vtan@gsnFrame			     = False

  ; Options for plotting 6-km mass flux 
    optsL_mflux		      	             = True
    optsL_mflux@gsLineColor		     = "black"
    optsL_mflux@gsLineThicknessF             = 2.0

    optsM_mflux                              = True
    optsM_mflux@gsMarkerIndex		     = 16
    optsM_mflux@gsMarkerColor                = "black"
    optsM_mflux@gsMarkerSizeF		     = 5.0
    optsM_mflux@gsMarkerThicknessF           = 2.0

  ; Options for additional model-level mass flux profiles
    multi_colours	   	       	     = (/"orange2","purple3","cyan","deeppink"/)

    optsL_mflux2	   	       	     = optsL_mflux
    optsL_mflux2@gsLineColor		     = multi_colours(0)
    optsM_mflux2			     = optsM_mflux
    optsM_mflux2@gsMarkerColor		     = multi_colours(0)

    optsL_mflux3           		     = optsL_mflux
    optsL_mflux3@gsLineColor                 = multi_colours(1)
    optsM_mflux3                             = optsM_mflux
    optsM_mflux3@gsMarkerColor               = multi_colours(1)

    optsL_mflux4           		     = optsL_mflux
    optsL_mflux4@gsLineColor                 = multi_colours(2)
    optsM_mflux4                             = optsM_mflux
    optsM_mflux4@gsMarkerColor               = multi_colours(2)

    optsL_mflux5           		     = optsL_mflux
    optsL_mflux5@gsLineColor                 = multi_colours(3)
    optsM_mflux5                             = optsM_mflux
    optsM_mflux5@gsMarkerColor               = multi_colours(3)

    if (flux .eq. 0) then 
     plot_w0 = plot_w0 / 10000
     plot_w1 = plot_w1 / 10000
    end if 

  ; Add conditions later in the script:
  ; if lev2,lev3,etc are greater than lev1 (upper), colour lev0 in bold and black
  ; likewise if lev2,lev3 are greater than lev0 but less than lev1, lev1 is bold and black

    if (multi .eq. 1 .and. lev_m2 .gt. lev_m1) then 

     opts_mflux0@xyLineColors	             = "black"
     opts_mflux0@xyMarkerColors		     = "black"
     opts_mflux0@xyLineThicknesses           = 4.0

     optsL_mflux@gsLineColor		     = "red3"
     optsL_mflux@gsLineThicknessF	     = 2.0
     optsM_mflux@gsMarkerColor		     = "red3"
     optsM_mflux@gsMarkerThicknessF	     = 2.0

    else 

     opts_mflux0@xyLineColors                = "red3"
     opts_mflux0@xyMarkerColors              = "red3"
     opts_mflux0@xyLineThicknesses           = 2.0

     optsL_mflux@gsLineColor                 = "black"
     optsL_mflux@gsLineThicknessF            = 4.0
     optsM_mflux@gsMarkerColor               = "black"
     optsM_mflux@gsMarkerThicknessF          = 4.0

    end if 

  ; Read data into panels for plotting (lower model level mass flux)
    if (flux .eq. 1) then 
     result   = gsn_csm_xy2(wks, time_int, plot_flux0, plot_vtan, opts_mflux0, opts_vtan)
    elseif (flux .eq. 0) then 
     result   = gsn_csm_xy2(wks, time_int, plot_w0, plot_vtan, opts_mflux0, opts_vtan)    
    end if 

  ; Now overlay 3rd line graph - line & markers separately - higher model level mass flux
    if (flux .eq. 1) then
     fluxL    = gsn_add_polyline(wks, result, time_int, plot_flux1, optsL_mflux)
     fluxM    = gsn_add_polymarker(wks, result, time_int, plot_flux1, optsM_mflux)
    elseif (flux .eq. 0) then
     fluxL    = gsn_add_polyline(wks, result, time_int, plot_w1, optsL_mflux)
     fluxM    = gsn_add_polymarker(wks, result, time_int, plot_w1, optsM_mflux)   
    end if 

  ; Also overlay additional mass flux plots
    if (flux .eq. 1 .and. multi .eq. 1) then 
     flux2L   = gsn_add_polyline(wks, result, time_int, plot_flux2, optsL_mflux2)	
     flux2M   = gsn_add_polymarker(wks, result, time_int, plot_flux2, optsM_mflux2)
     flux3L   = gsn_add_polyline(wks, result, time_int, plot_flux3, optsL_mflux3)
     flux3M   = gsn_add_polymarker(wks, result, time_int, plot_flux3, optsM_mflux3)
     flux4L   = gsn_add_polyline(wks, result, time_int, plot_flux4, optsL_mflux4)
     flux4M   = gsn_add_polymarker(wks, result, time_int, plot_flux4, optsM_mflux4)
     flux5L   = gsn_add_polyline(wks, result, time_int, plot_flux5, optsL_mflux5)
     flux5M   = gsn_add_polymarker(wks, result, time_int, plot_flux5, optsM_mflux5)
    end if 

;===============================
; Add and customise a legend
;===============================

    lg_opts                            = True
    lg_opts@pmLegendDisplayMode        = "Always"          ; Display a legend
    lg_opts@lgAutoManage               = False             ; Gives us control
    lg_opts@lgLabelFont                = "Helvetica"       ; Label font
    lg_opts@lgLabelFontHeightF         = 0.040             ; Height of labels

  ; Customise labels and line colours
    if (multi .eq. 1) then 
     lg_opts@lgLineColors               = (/optsL_mflux5@gsLineColor, \		; plot_flux5
     					    optsL_mflux4@gsLineColor, \		; plot_flux4
					    optsL_mflux3@gsLineColor, \		; plot_flux3
     					    opts_mflux0@xyLineColors, \		; Lower level
     					    optsL_mflux2@gsLineColor, \		; plot_flux2
     				       	    optsL_mflux@gsLineColor, \		; Higher level
 					    opts_vtan@xyLineColors/) 		; Windspeed
     lg_opts@lgDashIndexes              = (/0, 0, 0, 0, 0, 0, 0/)
     lsize			        = dimsizes(lg_opts@lgDashIndexes)
    else
     lg_opts@lgLineColors               = (/opts_mflux0@xyLineColors, \         ; Lower level
                                            optsL_mflux@gsLineColor, \           ; Higher level
                                            opts_vtan@xyLineColors/)             ; Windspeed
     lg_opts@lgDashIndexes              = (/0, 0, 0/)
     lsize                              = dimsizes(lg_opts@lgDashIndexes)
    end if 
    
    lg_opts@lgLabelPosition            = "Right"           ; Label position within box
    lg_opts@lgItemPlacement            = "ExplicitPlacement"

    if (multi .eq. 1) then 
     lg_opts@lgItemPositions            = fspan(0.02, 0.98, lsize)
     lg_opts@vpWidthF                   = 0.25
     lg_opts@vpHeightF                  = 0.12
    else
     lg_opts@lgItemPositions            = fspan(0.20, 0.80, lsize)
     lg_opts@vpWidthF                   = 0.22
     lg_opts@vpHeightF                  = 0.10
    end if 

    lg_opts@lgLineThicknessF           = 4.0
    lg_opts@lgLineDashSegLenF          = 0.20
    lg_opts@lgPerimColor               = "black"
    lg_opts@lgPerimThicknessF          = 3.0
    lg_opts@lgPerimFill                = 0           ; Solid fill
    lg_opts@lgPerimFillColor           = "white"

    if (flux .eq. 1 .and. multi .eq. 1) then 

      lab_arr                            = (/lev_m5+" m	mass flux", \
      					     lev_m4+" m	mass flux", \
					     lev_m3+" m	mass flux", \
      					     lev_m0+" m mass flux", \
      					     lev_m2+" m mass flux", \
                                             lev_m1+" m mass flux", \
                                             "Max. vtan"/)

    elseif (flux .eq. 1) then     

     if (vtan .eq. 1) then 
      lab_arr                            = (/lev_m0+" m mass flux", \
      					     lev_m1+" m mass flux", \
      				       	     "Max. vtan (pressure levels)"/)
     elseif (vtan .eq. 2) then
      lab_arr                            = (/lev_m0+" m mass flux", \
                                             lev_m1+" m mass flux", \
                                             "Max. vtan (model levels)"/)
     elseif (vtan .eq. 0) then 
      lab_arr                            = (/lev_m0+" m mass flux", \
      					     lev_m1+" m mass flux", \
                                             "10-m wind"/)    
     end if 

    elseif (flux .eq. 0) then 

     if (vtan .eq. 1 .or. vtan .eq. 2) then
      lab_arr                            = (/lev_m0+" m vertical velocity", \
                                             lev_m1+" m vertical velocity", \
                                             "Max. vtan"/)
     elseif (vtan .eq. 0) then
      lab_arr                            = (/lev_m0+" m vertical velocity", \
                                             lev_m1+" m vertical velocity", \
                                             "10-m wind"/)
     end if

    end if 

    lgd                                = gsn_create_legend(wks,lsize,lab_arr,lg_opts)

  ; Add legend to plot, and place it in top right hand corner
    am_opts                            = True
    am_opts@amJust                     = "TopLeft"
    am_opts@amParallelPosF             = -0.5
    am_opts@amOrthogonalPosF           = -0.5
    annotate_1                         = gsn_add_annotation(result,lgd,am_opts)

    draw(result)
    frame(wks)

end