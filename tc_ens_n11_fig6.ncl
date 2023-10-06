; Script to reproduce Figure 6 from Nguyen et al. (2011)

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_rm.ncl"

; ncl 'dat="02T12"' 'opt="x11"' 'ens0="em11"' dist=3.0 ts=39 tf=84 p_cp=3
; tclr=2 rclr=2 wclr=1 lay=1 p_cp2=7 cnt=0 tc_ens_n11_fig6.ncl

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'ens0'  = ensemble simulation (em00-em11)
; 'dist'  = size of box following storm (degrees)
; 'opt'   = output file format ("pdf" or "x11")
; 'p_cp'  = pressure level for plotting PV and vertical velocity (see L106-108)
; 'tclr'  = colour scheme for plotting tangential wind (0-2)
; 'rclr'  = colour scheme for plotting radial wind (0-2)
; 'wclr'  = colour for plotting vertical velocity (0-2)
; 'lay'	  = calculate layer average (1) or single level (0)
; 'p_cp2' = upper model level for calculating layer average (only valid when lay=1)
; 'cnt'   = use offline (1) or in-script (0) calculations of storm-centred grid	

begin

;==============================
; Prelash (reading in files)
;==============================

; Find 'pb' stream data to read in and plot
  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0
  fili_b1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pb*.nc")
  fili_b      = diri+"/"+fili_b1

; Calculate number of times in file using built-in function ('nc_times')
  numINFO     = nc_times(fili_b)
  numINPUT    = numINFO[0]

; Calculate storm motion using built-in function ('st_rm')
; 'ens0'    = ensemble simulation
; 'diri'    = path to input files (see above)
; 'dist'    = size of box following storm (degrees)

  storm_rel = st_rm(dat,ens0,"$ar/text/",dist,0)
  u_cyc     = storm_rel[0]             ; Zonal wind
  v_cyc     = storm_rel[1]             ; Meridional wind
  vel_cyc   = storm_rel[2]             ; Vector wind
  numTIMES  = storm_rel[3]             ; Number of times in file
  lat_arr   = storm_rel[4]             ; Latitude array subset (following storm)
  lon_arr   = storm_rel[5]             ; Longitude array subset (following storm)
  centre    = storm_rel[6]             ; Storm track (position) array

; Create arrays for data and output information
  date_str  = new(numINPUT,string)
  time_str  = new(numINPUT,string)
  time_arr  = new(numINPUT,string)
  title_arr = new(numINPUT,string)
  llbox     = toint(dist*50)                   ; Calculate domain size (grid points)
  dsize     = (/llbox,llbox/)                  ; Domain size (grid points)

; Extract string for each ensemble member
  str1    = str_split(fili_b(0),"_")
  ens_str = str1(5)
  delete(str1)

  print("Working on ensemble member: "+ens_str)

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 10000000000
  end setvalues

;==================================================================================
;;;;;;;;;;;;;;;;;;;;;;;;;;; Start of main part of script ;;;;;;;;;;;;;;;;;;;;;;;;;;
;==================================================================================

;========================================
; Define cylindrical coordinate arrays
;========================================

  thetas     = new(73,float)
  pii        = 3.14159265
  radii      = fspan(0.0,1.5,31)     ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,73)   ; Array of azimuth angles (0, 5, 10, etc)
  thetas_rad = (pii/180)*thetas      ; Degrees to radians
  rad_size   = dimsizes(radii)-1     ; Index for later calculation

; Define arrays for later use in tangential velocity cross-section plots
  vtan_azi = new((/dimsizes(radii)/),float) ; Tangential wind
  vrad_azi = new((/dimsizes(radii)/),float) ; Radial wind
  spd_azi  = new((/dimsizes(radii)/),float) ; Vector velocity
  vvel_azi = new((/dimsizes(radii)/),float) ; Vertical velocity
  the_azi  = new((/dimsizes(radii)/),float) ; Theta-e

; Define 2-dimensional array to hold all values from 'vtan_azi' at ALL times
  vtan_plot0   = new((/dimsizes(radii),numTIMES/),float)
  vtan_plot0!0 = "rad"
  vtan_plot0!1 = "time"

  vrad_plot0   = new((/dimsizes(radii),numTIMES/),float)
  vrad_plot0!0 = "rad"
  vrad_plot0!1 = "time"

  vvel_plot0   = new((/dimsizes(radii),numTIMES/),float)
  vvel_plot0!0 = "rad"
  vvel_plot0!1 = "time"

  the_plot0    = new((/dimsizes(radii),numTIMES/),float)
  the_plot0!0  = "rad"
  the_plot0!1  = "time"

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

  vvel_azi!0            = "rad"
  vvel_azi&rad          = radii
  vvel_azi@description  = "Azimuthally averaged vertical velocity"
  vvel_azi@units        = "m s~S~-1~N~"

  the_azi!0             = "rad"
  the_azi&rad           = radii
  the_azi@description   = "Azimuthally averaged equivalent potential temperature"
  the_azi@units         = "K"

; Create time array for later plotting 
  time_int = new(120,"integer")
  time_int = ispan(0,119,1)
  time_int!0 = "time"

;==================================
; Get the variables we will need
;==================================

  f = addfile(fili_b,"r")

  do it = ts, tf                   ; Loop over times in file

   time  = f->t(it)                ; Time in file
   times = dimsizes(time)

 ; Account for the first three missing times in the input file
   ct    = it + 3

 ; Get domain boundaries at each time from storm track (calculated off-line)
   t0    = lat_arr(ct,0)
   t1    = lat_arr(ct,1)
   n0    = lon_arr(ct,0)
   n1    = lon_arr(ct,1)

 ; Create latitude/longitude array subsets to reduce computation time
   lon10 = f->longitude_1({n0:n1})  ; longitude [160 grid points]
   lat10 = f->latitude_1({t0:t1})   ; latitude  [160 grid points]
   lon0  = f->longitude({n0:n1})    ; longitude [160 grid points]
   lat0  = f->latitude({t0:t1})     ; latitude  [160 grid points]

 ; Pressure levels (pres)
 ; (0) 1000, (1) 950, (2) 925, (3) 900, (4) 850, (5) 800, (6) 750, (7) 700
 ; (8) 650, (9) 600, (10) 500, (11) 400, (12) 300, (13) 250
 ; (14) 200, (15) 150, (16) 100

   if (lay .eq. 0) then 
    pres  = f->p(p_cp)              ; Subset of pressure level array
   elseif (lay .eq. 1) then 
    pres = f->p(p_cp:p_cp2)
   end if 

   plevs = dimsizes(pres)           ; Size of pressure level array
   pres@units = "hPa"               ; Pressure level units (for plotting)

 ; Arrays have () vertical levels, and horizontal extent depends on 'dist' option

   if (lay .eq. 1) then 
    u0     = f->u(it,p_cp:p_cp2,{t0:t1},{n0:n1})       ; Zonal velocity (m/s)
    v0     = f->v(it,p_cp:p_cp2,{t0:t1},{n0:n1})       ; Meridional velocity (m/s)
    vort0  = f->rvor(it,p_cp:p_cp2,{t0:t1},{n0:n1})    ; Relative vorticity (/s)
    t00    = f->temp(it,p_cp:p_cp2,{t0:t1},{n0:n1})    ; Temperature (K)
    z0     = f->ht(it,p_cp:p_cp2,{t0:t1},{n0:n1})      ; Geopotential height (m)
    w0     = f->dz_dt(it,p_cp:p_cp2,{t0:t1},{n0:n1})   ; Vertical velocity (m/s)
    thw0   = f->field22(it,p_cp:p_cp2,{t0:t1},{n0:n1}) ; Wet bulb potential temperature (K)
    rh0	  = f->rh(it,p_cp:p_cp2,{t0:t1},{n0:n1})       ; Relative humidity 
   elseif (lay .eq. 0)
    u0     = f->u(it,p_cp,{t0:t1},{n0:n1})	       ; Zonal velocity (m/s)
    v0     = f->v(it,p_cp,{t0:t1},{n0:n1})             ; Meridional velocity (m/s)
    vort0  = f->rvor(it,p_cp,{t0:t1},{n0:n1})          ; Relative vorticity (/s)
    t00    = f->temp(it,p_cp,{t0:t1},{n0:n1})          ; Temperature (K)
    z0     = f->ht(it,p_cp,{t0:t1},{n0:n1})            ; Geopotential height (m)
    w0     = f->dz_dt(it,p_cp,{t0:t1},{n0:n1})         ; Vertical velocity (m/s)
    thw0   = f->field22(it,p_cp,{t0:t1},{n0:n1})       ; Wet bulb potential temperature (K)
    rh0   = f->rh(it,p_cp,{t0:t1},{n0:n1})             ; Relative humidity
   end if 

 ; Grid information for all variables (for later switch to cylindrical coordinates)
 ; u0     = [lat , lon]
 ; v0	  = [lat , lon]
 ; vort0  = [lat , lon]
 ; pv0	  = [lat , lon]
 ; t00	  = [lat1, lon1]
 ; z0	  = [lat1, lon1]
 ; w0	  = [lat1, lon1]
 ; thw0	  = [lat1, lon1]
 ; rh0	  = [lat1, lon1]

 ;=========================================
 ; Calculate theta-e using NCL functions
 ;=========================================

; ; First calculate mixing ratio from pressure, temperature and RH 
; ; Option '1' returns the mixing ratio; '2' returns the specific humidity
;
; ; Create array with same dimensions as 't00' for mixing ratio calculation
;   tsize  = dimsizes(t00) 
;   pres0  = new((/tsize(0),tsize(1)/),float)
;   pres0  = pres
;   delete(tsize)
;
; ; Now calculate the mixing ratio from pressure, temperature and RH
;   mix0   = mixhum_ptrh(pres0,t00,rh0,1)
;
; ; Then calculate the lifting condensation level temperature using RH 
; ; Option '(/1,1/)' specifies units of input and output temperatures (K)
;   tlcl   = tlcl_rh_bolton(t00,rh0,(/1,1/))
;
; ; Finally calculate equivalent potential temperature
; ; Option '(/0,1,0,1/)' -- input p (hPa), temp (K), mixing ratio (kg/kg), theta-e (K)
;   the0	  = pot_temp_equiv_tlcl(pres,t0,tlcl,mix0,(/0,1,0,1/))
;

   vel0   = sqrt((u0^2) + (v0^2))             ; Vector wind (m/s)
   copy_VarCoords(v0,vel0)
   copy_VarMeta(v0,vel0)

 ; Avoid dimension size errors throughout the script
   usize    = dimsizes(u0)
   vsize    = dimsizes(v0)
   vrsize   = dimsizes(vort0)
   tsize    = dimsizes(t00)
   zsize    = dimsizes(z0)
   wsize    = dimsizes(w0)
   vcsize   = dimsizes(vel0)
   thsize   = dimsizes(thw0)
   latsize  = dimsizes(lat0)
   lonsize  = dimsizes(lon0)
   lat1size = dimsizes(lat10)
   lon1size = dimsizes(lon10)

   if ( any (usize .gt. dsize(0)) ) then
    if (lay .eq. 0) then 
     u = u0(0:dsize(0)-1,0:dsize(0)-1)
    elseif (lay .eq. 1) then 
     u = u0(:,0:dsize(0)-1,0:dsize(0)-1)
    end if 
   else
    u = u0
   end if

   if ( any (vsize .gt. dsize(0)) ) then
    if (lay .eq. 0) then
     v = v0(0:dsize(0)-1,0:dsize(0)-1)
    elseif (lay	.eq. 1) then
     v = v0(:,0:dsize(0)-1,0:dsize(0)-1)
    end if
   else
    v = v0
   end if

   if ( any (vrsize .gt. dsize(0)) ) then
    if (lay .eq. 0) then
     vort = vort0(0:dsize(0)-1,0:dsize(0)-1)
    elseif (lay .eq. 1) then
     vort = vort0(:,0:dsize(0)-1,0:dsize(0)-1)
    end if
   else
    vort = vort0
   end if

   if ( any (tsize .gt. dsize(0)) ) then
    if (lay .eq. 0) then
     t = t00(0:dsize(0)-1,0:dsize(0)-1)
    elseif (lay .eq. 1) then
     t = t00(:,0:dsize(0)-1,0:dsize(0)-1)
    end if
   else
    t = t00
   end if

   if ( any (zsize .gt. dsize(0)) ) then
    if (lay .eq. 0) then
     z = z0(0:dsize(0)-1,0:dsize(0)-1)
    elseif (lay .eq. 1) then
     z = z0(:,0:dsize(0)-1,0:dsize(0)-1)
    end if
   else
    z = z0
   end if

   if ( any (wsize .gt. dsize(0)) ) then
    if (lay .eq. 0) then
     w = w0(0:dsize(0)-1,0:dsize(0)-1)
    elseif (lay .eq. 1) then
     w = w0(:,0:dsize(0)-1,0:dsize(0)-1)
    end if
   else
    w = w0
   end if

   if ( any (vcsize .gt. dsize(0)) ) then
    if (lay .eq. 0) then
     vel = vel0(0:dsize(0)-1,0:dsize(0)-1)
    elseif (lay .eq. 1) then
     vel = vel0(:,0:dsize(0)-1,0:dsize(0)-1)
    end if
   else
    vel = vel0
   end if

   if ( any (thsize .gt. dsize(0)) ) then
    if (lay .eq. 0) then
     th_w = thw0(0:dsize(0)-1,0:dsize(0)-1)
    elseif (lay .eq. 1) then
     th_w = thw0(:,0:dsize(0)-1,0:dsize(0)-1)
    end if
   else
    th_w = thw0
   end if

 ;====================================================================
 ; Also follow the same procedure for latitude and longitude arrays 
 ; These arrays are only 1D, in contrast to the 2D arrays above 
 ;====================================================================

   if ( any (latsize .gt. dsize(0)) ) then
    lat = lat0(0:dsize(0)-1)
   else
    lat = lat0
   end if

   if ( any (lonsize .gt. dsize(0)) ) then
    lon = lon0(0:dsize(0)-1)
   else
    lon = lon0
   end if

   if ( any (lat1size .gt. dsize(0)) ) then
    lat1 = lat10(0:dsize(0)-1)
   else
    lat1 = lat10
   end if

   if ( any (lon1size .gt. dsize(0)) ) then
    lon1 = lon10(0:dsize(0)-1)
   else
    lon1 = lon10
   end if

 ; Convert to units we want
   z     = z/10             ; Convert to decametres
   z@units = "dam"          ; Add metadata
   vort  = vort * (10 ^ -6) ; Convert to /s

   dx    = lon(1) - lon(0) ; Grid spacing (longitude)
   dy    = lat(1) - lat(0) ; Grid spacing (latitude)

 ; Coriolis parameter and absolute vorticity
   f1   = coriolis_param(lat)

   if (lay .eq. 0) then 
    f0  = conform_dims(dimsizes(vort),f1,0)     ; Give 'f1' the dimensions of 'vort'
   elseif (lay .eq. 1) then 
    f0  = conform_dims(dimsizes(vort),f1,1)     ; Give 'f1' the dimensions of 'vort'
   end if 

   copy_VarMeta(vort,f0)
   copy_VarCoords(vort,f0)

   avo  = vort + f0
   copy_VarMeta(vort,avo)
   copy_VarAtts(vort,avo)
   copy_VarCoords(vort,avo)
   avo@description = "Absolute vorticity"
   avo@long_name   = "Absolute vorticity"

;===================================
; Calculate storm-relative winds
;===================================

   uf	= u - u_cyc(it-3)
   vf	= v - v_cyc(it-3)

 ; Add metadata from original horizontal wind arrays
   copy_VarCoords(u,uf)
   copy_VarAtts(u,uf)

   copy_VarCoords(v,vf)
   copy_VarAtts(v,vf)

;=====================================================
; Create correct date strings for each output time
;=====================================================

   month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
                     "Jul","Aug","Sep","Oct","Nov","Dec"/)

   utc_date = cd_calendar(time, 0)

   year   = tointeger(utc_date(:,0))
   month  = tointeger(utc_date(:,1))
   day    = tointeger(utc_date(:,2))
   hour   = tointeger(utc_date(:,3))
   minute = tointeger(utc_date(:,4))
   second = utc_date(:,5)

  ; Correct for errors in the code (round up value of hour when minutes = 59)
    if (.not.ismissing(minute).and.minute.gt.30) then
      hour = hour+1
    end if

   date_str(it)  = sprinti("%0.2iUTC ", hour) + \
                   sprinti("%0.2i ", day) \
                   + month_abbr(month)
   time_str(it)  = sprinti("%0.2i UTC ", hour) + \
                   sprinti("%0.2i ", day) \
                   + month_abbr(month)
   time_arr(it)  = sprinti("%0.2i", day) + \
                   month_abbr(month) + \
                   "_" + sprinti("%0.2iZ", hour)

   print("Working on time: "+time_str(it))
   print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")
   itt = ct + 1
   title_arr(it) = "valid at "+time_str(it)+" (T+"+itt+")"

;==============================================
; Analyse variables on chosen pressure level
;==============================================

     print("Working on: " + pres + " hPa")

     if (lay .eq. 1) then 
      u_plane       = dim_avg_n_Wrap(uf(:,:,:),0)
      v_plane       = dim_avg_n_Wrap(vf(:,:,:),0)
      vort_plane    = dim_avg_n_Wrap(vort(:,:,:),0)
      geo_plane     = dim_avg_n_Wrap(z(:,:,:),0)
      w_plane	    = dim_avg_n_Wrap(w(:,:,:),0)
     else 
      u_plane       = uf(:,:)
      v_plane       = vf(:,:)
      vort_plane    = vort(:,:)
      geo_plane     = z(:,:)
      w_plane       = w(:,:)
     end if 

   ; Max and min geopotential height on	pressure level
     zmin          = toint(min(geo_plane))
     zmax          = toint(max(geo_plane))

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

     lat_max = t0 + (dy * inds(0,0))
     lon_max = n0 + (dx * inds(0,1))

     print("ct = "+ct+" ; it = "+it)
     print("lat_max = "+centre(ct,0)+" ; lon_max = "+centre(ct,1))

   ; Difference between 'lat_max'/'lon_max' and values read in using NCL function
   ; Use 'ct' rather than 'it' because 'pb' stream starts 4 h into simulation
     dlat                   = abs(lat_max - centre(ct,0))
     dlon                   = abs(lon_max - centre(ct,1))

   ; If our calculations of latitude and longitude deviate too much from the values
   ; we read in using the NCL function ('st_rm'), then calculate a new centre
   ; using the minimum pressure. Constrain the new calculation using the values
   ; from 'st_rm', to prevent spurious centres being identified and analysed.

     if (dlat .gt. 0.75 .or. dlon .gt. 0.75) then

    ; EDITED FROM 'it' TO 'ct' -- CHECK THAT RESULTS ARE CORRECT WHEN RERUNNING
      lt1 = centre(ct,0) - 0.5
      lt2 = centre(ct,0) + 0.5
      ln1 = centre(ct,1) - 0.5
      ln2 = centre(ct,1) + 0.5

      geo_plane0 = geo_plane({lt1:lt2},{ln1:ln2})       ; Smaller grid (1.0 x 1.0 deg)
      geo_min0   = min(geo_plane0)                      ; Find minimum on smaller grid
      dims_p0    = dimsizes(geo_plane0)                 ; Size of smaller grid
      geo_1d0    = ndtooned(geo_plane0)                 ; Create 1-D array
      inds_p0    = ind_resolve(minind(geo_1d0),dims_p0) ; Find index of minimum

      lat_max    = lt1 + (dy * inds_p0(0,0))
      lon_max    = ln1 + (dx * inds_p0(0,1))
      print("New centre: "+lat_max+" degrees N, "+lon_max+" degrees E")

     end if

   ;==============================================================
   ; Create arrays for later switch to cylindrical coordinates
   ;==============================================================

   ; Indices of variables on 'lat'/'lon' grid for 'pb' stream
     lat_plane  = ind(t0.le.lat.and.lat.le.t1) ; Indices of latitude array
     lon_plane  = ind(n0.le.lon.and.lon.le.n1) ; Indices of longitude array

   ; Array sizes
     lat_size = dimsizes(lat_plane)
     lon_size = dimsizes(lon_plane)

   ; Array of latitude/longitude points corresponding to the indices above
     lat_sub = lat(lat_plane(0):lat_plane(lat_size-1))
     lon_sub = lon(lon_plane(0):lon_plane(lon_size-1))

   ; Indices of variables on 'lat1'/'lon1' grid for 'pb' stream
     lat1_plane  = ind(t0.le.lat1.and.lat1.le.t1) ; Indices of latitude array
     lon1_plane  = ind(n0.le.lon1.and.lon1.le.n1) ; Indices of longitude array

   ; Array sizes
     lat1_size = dimsizes(lat1_plane)
     lon1_size = dimsizes(lon1_plane)

   ; Array of latitude/longitude points corresponding to the indices above
     lat1_sub = lat1(lat1_plane(0):lat1_plane(lat1_size-1))
     lon1_sub = lon1(lon1_plane(0):lon1_plane(lon1_size-1))

     delete([/vort1d,inds,vort1d_sm,vort_smth,inds_sm/])

;===========================================
; Switch to cylindrical coordinate system
;===========================================

   ; Create new arrays describing cylindrical coordinate system (2D for now)
     xpos = new((/dimsizes(radii),dimsizes(thetas)/),float)
     ypos = new((/dimsizes(radii),dimsizes(thetas)/),float)

   ; New arrays containing variables we want to plot
     u_int   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Zonal velocity
     v_int   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Meridional velocity
     spd_int = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vector velocity (RMW)
     w_int   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vertical velocity
;     th_int  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Equiv. potential temp.
     v_rad   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Radial velocity
     v_tan   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Tangential velocity

   ; Add variable metadata
     copy_VarAtts(u_plane,u_int)
     u_int@description = "Zonal wind"
     u_int@units       = "m s~S~-1~N~"

     copy_VarAtts(v_plane,v_int)
     v_int@description = "Meridional wind"
     v_int@units       = "m s~S~-1~N~"

     copy_VarAtts(w_plane,w_int)
     w_int@description = "Vertical wind"
     w_int@units       = "m s~S~-1~N~"

;     copy_VarAtts(the_plane,th_int)
;     th_int@description = "Equivalent potential temperature"
;     th_int@units       = "K"

     v_rad@description = "Radial velocity"
     v_rad@units       = "m s~S~-1~N~"

     v_tan@description = "Tangential velocity"
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

;     th_int!0   = "rad"
;     th_int!1   = "azi"
;     th_int&rad = radii
;     th_int&azi = thetas

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
            xpos(irad,iang) = centre(ct,1) + r*cos(theta)
            ypos(irad,iang) = centre(ct,0) + r*sin(theta)
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
     
;     print("lon (offline) = "+centre(ct,1)+" ; lat (offline) = "+centre(ct,0))
;     print("lon (script) = "+lon_max+" ; lat (script) = "+lat_max)

  ; Loop over azimuth angles to fill entire circular array of data points
    do iang = 0, dimsizes(thetas_rad)-1
      u_int(:,iang)    = linint2_points(lon_sub,lat_sub,u_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      v_int(:,iang)    = linint2_points(lon_sub,lat_sub,v_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      w_int(:,iang)    = linint2_points(lon1_sub,lat1_sub,w_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
;      th_int(:,iang)   = linint2_points(lon1_sub,lat1_sub,the_plane,False,\
;                                        xpos(:,iang),ypos(:,iang),0)
    end do

;    printMinMax(u_int,False)
;    printMinMax(v_int,False)
 
  ; u0 	  =  [lat  , lon]
  ; v0     = [lat  , lon]
  ; vort0  = [lat  , lon]
  ; pv0    = [lat  , lon]
  ; t00    = [lat1 , lon1]
  ; z0     = [lat1 , lon1]
  ; w0     = [lat1 , lon1]
  ; thw0   = [lat1 , lon1]
  ; rh0    = [lat1 , lon1]

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
        vtan_azi(irad) = avg(v_tan(irad,:))      ; Tangential wind
	vrad_azi(irad) = avg(v_rad(irad,:))	 ; Radial wind 
	spd_azi(irad)  = avg(spd_int(irad,:))	 ; Vector wind (for RMW)
;	the_azi(irad)  = avg(th_int(irad,:))	 ; Theta-e
	vvel_azi(irad) = avg(w_int(irad,:))	 ; Vertical velocity
      end do

    ; Input azimuthal averages into larger arrays before plotting  
      vtan_plot0(:,ct)  = vtan_azi(:)            ; Tangential wind
      vrad_plot0(:,ct)  = vrad_azi(:)	    	 ; Radial wind
;      the_plot0(:,ct)	= the_azi(:)		 ; Theta-e
      vvel_plot0(:,ct)  = vvel_azi(:)	    	 ; Vertical velocity

    ; Tidy up
      delete([/u_int,v_int,w_int,v_tan,v_rad/])
      delete([/u_plane,v_plane,w_plane,vort_plane,geo_plane/])
;      delete([/pv_int,pv_plane/])

    ct  = ct + 1		    ; Counter variable (time)

  ; Tidy up (avoid dimension size errors)
    delete([/lat1,lon1,lat,lon,time,times,w,u,v,z,t,vel/])
    delete([/lat0,lon0,lat10,lon10,w0,u0,v0,z0,t00,vel0,vort0,rh0/])
    delete([/f0,f1,avo,vort/])
  ; delete([/mix0,tlcl,the0,pres0/])
    delete([/utc_date,year,month,day,hour,minute,second/])

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

  do it = ts, tf

   tp1  = min((/it+1,tf/))
   tm1  = max((/it-1,ts/))

   vtan_tend0(:,it) = (vtan_plot0(:,tp1) - vtan_plot0(:,tm1)) / 2 ; Units [m s-1 h-1]
   vrad_tend0(:,it) = (vrad_plot0(:,tp1) - vrad_plot0(:,tm1)) / 2 ; Units [m s-1 h-1]

  end do     ; End time loop (do it = ts, tf)

;======================================
; Reorder dimensions before plotting 
;======================================

 ; Arrays now ordered [TIME * RADIUS] rather than [RADIUS * TIME]
   vtan_tend = vtan_tend0(time|:,rad|:)
   vrad_tend = vrad_tend0(time|:,rad|:)

   vtan_plot = vtan_plot0(time|:,rad|:)
   vrad_plot = vrad_plot0(time|:,rad|:)
   vvel_plot = vvel_plot0(time|:,rad|:)
;   the_plot  = the_plot0(time|:,rad|:)

;=====================================================
; Create a panel plot of Figs. 6a and 6b from N11 
;=====================================================

  ; Get pressure-level strings for output file 
    delete(pres)
    pres = f->p
    p0 = pres(p_cp)
    p1 = pres(p_cp2)

  ; Output file location and type
    if (lay .eq. 1) then 
     output = "$sam/nepartak/images/cyl_coords/n11_fig6_"+dat+"_"+ens0+"_"+p0+"_"+p1
    elseif (lay .eq. 0) then 
     output = "$sam/nepartak/images/cyl_coords/n11_fig6_"+dat+"_"+ens0+"_"+p0
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
    opts_vtan@cnLinesOn                    = False ; Contour lines off
    opts_vtan@cnInfoLabelOn                = False
    opts_vtan@cnLevelSelectionMode         = "ExplicitLevels"

    if (tclr .eq. 0) then 
     opts_vtan@cnLevels                     = (/1.0, 2.0, 2.5, 3.0, 3.5, 4.0, \
                                                4.5, 5.0, 6.0, 7.0/)
     opts_vtan@cnFillColors                 = (/0,2,3,4,5,6,7,8,9,10,11/)
    else if (tclr .eq. 1) then 
     opts_vtan@cnLevels                     = (/0.5, 1.5, 2.0, 2.5, 3.0, 3.5, \
                                                4.0, 4.5, 5.0, 6.0/)
     opts_vtan@cnFillColors                 = (/0,2,3,4,5,6,7,8,9,10,11/)
    else if (tclr .eq. 2) then		    
     opts_vtan@cnLevels                     = (/0.5, 1.0, 1.5, 2.0, 2.5, 3.0, \
                                                3.5, 4.0, 4.5, 5.0/)
     opts_vtan@cnFillColors                 = (/0,2,3,4,5,6,7,8,9,10,11/)
    end if 
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
    opts_vtan@trYMinF                      = ts+3    ; Starting point (T+...) for plot
    opts_vtan@trYMaxF                      = tf      ; Ending point (T+...) for plot
    opts_vtan@tmXTOn                       = "False" ; Turn off top x-axis TM
    opts_vtan@tmYROn                       = "False" ; Turn off right y-axis TM
    opts_vtan@gsnMaximize                  = True    ; Maximise plot size
    opts_vtan@gsnAddCyclic                 = False
    opts_vtan@tmXBMode                     = "Explicit"                ; Choose tm mode
    opts_vtan@tmXBValues                   = (/0.0,0.5,1.0,1.5/)       ; Tm positions
    opts_vtan@tmXBLabels                   = (/"0","50","100","150"/)  ; Tm values
    opts_vtan@tmXBMinorOn                  = False

    opts_vtan@pmLabelBarOrthogonalPosF     = 0.10    ; Move lb up/down (higher/lower)
    opts_vtan@tiXAxisOffsetYF              = 0.10    ; Move x-axis title up/down (h/l)

    opts_vtan@trXMinF                      = 0.0                       ; Set x-axis min
    opts_vtan@trXMaxF                      = 1.0                       ; Set x-axis max
    opts_vtan@tmYLMode                     = "Manual"
    opts_vtan@tmYLTickStartF               = opts_vtan@trYMinF         ; See L730
    opts_vtan@tmYLTickEndF                 = opts_vtan@trYMaxF         ; See L731
    opts_vtan@tmYLTickSpacingF             = 6
    opts_vtan@tmYLMinorOn                  = False

  ; Resources to overlay contours of mean tangential wind
    opts_cont                              = True
    opts_cont@cnFillOn                     = False       ; Filled contours on/off
    opts_cont@cnLineColor                  = "blue3"     ; Contour line colour
    opts_cont@cnLevelSelectionMode         = "ExplicitLevels"
    opts_cont@cnLevels                     = (/10,20,30,40,50,60,70/)
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
    opts_zero                             = True
    opts_zero@cnFillOn                    = False
    opts_zero@cnLineColor                 = "gray82"
    opts_zero@cnLevels                    = (/-1.0, 0.0, 1.0/)
    opts_zero@cnLineThicknessF            = 2.0
    opts_zero@cnInfoLabelOn               = False
    opts_zero@gsnDraw                     = False       ; Do not draw the plot
    opts_zero@gsnFrame                    = False       ; Do no advance the frame
    opts_zero@gsnContourZeroLineThicknessF = 2.0
    opts_zero@cnLineLabelsOn               = False      ; Turn off line labels
    opts_zero@tiMainString                = ""
    opts_zero@gsnLeftString               = ""
    opts_zero@gsnRightString              = ""

  ; Radial velocity resources
    opts_vrad                              = opts_vtan
    delete(opts_vrad@cnLevels)
    delete(opts_vrad@cnFillColors)

    if (rclr .eq. 0) then 
     opts_vrad@cnLevels                     = (/0.0, -0.5, -1.0, -1.5, -2.0, \
                                                -2.5, -3.0, -4.0, -5.0/)
     opts_vrad@cnFillColors                 = (/12,11,10,9,8,7,6,5,4,0/)
    else if (rclr .eq. 1) then 
     opts_vrad@cnLevels                     = (/0.0, -0.25, -0.5, -1.0, -1.5, \
                                                -2.0, -2.5, -3.0, -4.0/)
     opts_vrad@cnFillColors                 = (/12,11,10,9,8,7,6,5,4,0/)
    else if (rclr .eq. 2) then 
     opts_vrad@cnLevels                     = (/0.0, -0.25, -0.5, -0.75, -1.0, \
                                                -1.5, -2.0, -2.5, -3.0/)
     opts_vrad@cnFillColors                 = (/12,11,10,9,8,7,6,5,4,0/)
    end if
    end if
    end if 

  ; Vertical velocity resources
    opts_vvel                              = opts_vtan
    delete(opts_vvel@cnLevels)
    delete(opts_vvel@cnFillColors)

    opts_vvel@cnLevels                     = (/0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, \
                                               0.7, 0.8, 0.9, 1.0/)
    opts_vvel@cnFillColors                 = (/0,2,3,4,5,6,7,8,9,10,11,12/)

    if (wclr .eq. 0) then
     opts_vvel@cnLevels                     = (/0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, \
                                                0.7, 0.8, 0.9, 1.0/)
     opts_vvel@cnFillColors                 = (/0,2,3,4,5,6,7,8,9,10,11,12/)
    else if (wclr .eq. 1) then
     opts_vvel@cnLevels                     = (/0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, \
                                                0.50, 0.60, 0.70, 0.80/)
     opts_vvel@cnFillColors                 = (/0,2,3,4,5,6,7,8,9,10,11,12/)

    else if (wclr .eq. 2) then
     opts_vvel@cnLevels                     = (/0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, \
                                                0.40, 0.45, 0.50, 0.60/)
     opts_vvel@cnFillColors                 = (/0,2,3,4,5,6,7,8,9,10,11,12/)

    end if
    end if
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


;========================
; Construct panel plot
;========================

  ; Read in RMW data from text file (first 3 values equal to zero)
    spd_plot = asciiread("$ar/text/"+dat+"_"+ens0+"_rmw.txt",(/numTIMES/),"float")

  ; Tangential wind (LHS panel)
    vtan0 = gsn_csm_contour(wks,vtan_tend,opts_vtan)    ; Tangential wind tendency
    zero0 = gsn_csm_contour(wks,vtan_tend,opts_zero)    ; Zero line (vtan tendency)
    zero  = ColorNegDashZeroPosContour(zero0,"transparent","gray69","transparent")
    cont0 = gsn_csm_contour(wks,vtan_plot,opts_cont)    ; Tangential wind
    rmw0  = gsn_csm_xy(wks,spd_plot,time_int,opts_rmw)  ; RMW
    overlay(vtan0,zero)                                 ; Overlay vtan tendency zero line
    overlay(vtan0,cont0)                                ; Overlay tangential wind
    overlay(vtan0,rmw0)                                 ; Overlay RMW
    panel(0) = vtan0                                    ; Plot first panel
    delete([/zero0,zero/])                              ; Tidy up

  ; Radial wind (middle panel)
    vrad0 = gsn_csm_contour(wks,vrad_tend,opts_vrad)    ; Radial wind tendency
    cont0 = gsn_csm_contour(wks,vtan_plot,opts_cont)    ; Tangential wind
    rmw0  = gsn_csm_xy(wks,spd_plot,time_int,opts_rmw)  ; RMW
    overlay(vrad0,cont0)                                ; Overlay tangential wind
    overlay(vrad0,rmw0)                                 ; Overlay RMW
    panel(1) = vrad0                                    ; Plot second panel

  ; Vertical velocity (RHS panel)
    vvel0 = gsn_csm_contour(wks,vvel_plot,opts_vvel)    ; Vertical velocity
    zero0 = gsn_csm_contour(wks,vvel_plot,opts_zero)    ; Zero line (vertical velocity)
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
    optsP@txString              = "RA1T ensemble simulation ("+ens0+")"
    optsP@gsnPanelFigureStrings = (/"a) vtan_tend","b) vrad_tend","c) vvel"/)
    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.90
    optsP@gsnPanelBottom        = 0.10
    optsP@amJust                = "TopLeft"
    optsP@gsnPanelFigureStringsFontHeightF = 0.0099 ; Reduce label size (default 0.01)

    gsn_panel(wks,panel,(/1,3/),optsP)             ; Draw as a single plot
    frame(wks)

end 