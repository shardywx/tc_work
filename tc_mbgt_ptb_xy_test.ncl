; Script to calculate the radial momentum budget as in Montgomery et al. (2018)

; Run using:

; ncl 'opt="x11"' 'ens0="em10"' 'dat="02T12"' ts0=36 tf0=39
; ls0=20 lf0=21 dist=3.0 calc=1 rad0=1 tc_mbgt_ptb_xy.ncl

; 'opt'   = output file format ("pdf" or "x11")
; 'ens0'  = 
; 'dat'	  = 
; 'dist'  = size of box following storm (degrees)
; 'calc'  = centre calculated on each level (0) or equal to 950 hPa centre (1)
; 'rad0'  = radius of final plots: 150 km (0) or 200 km (1)
; 'ts0'	  = 
; 'tf0'	  = 
; 'ls0'	  = 
; 'lf0'	  = 

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"

begin

;==============================
; Prelash (reading in files)
;==============================

; Find 'pc' and 'pd' files (model height level data)
  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0
  fili_c1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pc*.nc")
  fili_c      = diri+"/"+fili_c1
  fili_d1     =	systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pd*.nc")
  fili_d      = diri+"/"+fili_d1

; Calculate number of times in each file using built-in function ('nc_times')
  numINFO_c   = nc_times(fili_c)
  numINFO_d   = nc_times(fili_d)
  numINPUT_c  = numINFO_c[0]		; Number of times in 'pc' stream
  time_c      = numINFO_c[1]		; Array of times from 'pc' stream file
  numINPUT_d  = numINFO_d[0]		; Number of times in 'pd' stream
  time_d      = numINFO_d[1]		; Array of times from 'pd' stream file 

; Calculate storm motion using built-in function ('st_rm')
; Argument 1  = ensemble simulation
; Argument 2  = path to input files
; Argument 3  = size of box following storm (degrees)

; Call external function 'st_rm'
  storm_rel = st_rm(dat,ens0,"$ar/text/",3.0)
  u_cyc     = storm_rel[0]		   ; Zonal wind 
  v_cyc     = storm_rel[1]		   ; Meridional wind
  vel_cyc   = storm_rel[2]		   ; Vector wind 
  numTIMES  = storm_rel[3]                 ; Number of times in file
  lat_arr   = storm_rel[4]                 ; Storm-following domain subset (lat)
  lon_arr   = storm_rel[5]                 ; Storm-following domain subset (lon)
  centre    = storm_rel[6]		   ; Storm centre at each hour (lat/lon)

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

  nr	     = 41
  ar	     = 2.0

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

;===============================
; Define arrays for plotting 
;===============================

; Array to hold all 'xpos' and 'ypos' position arrays
  xpos_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  ypos_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)

; Arrays to hold lat/lon at centre of cylindrical grid
  xcen_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  ycen_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)

; Arrays for perturbation (eddy) variables (4D)
  ptb_vrad0 = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [u']
  ptb_vtan0 = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [v']
  ptb_vvel0 = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [w']
  ptb_pres0 = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [p']
  ptb_vort0 = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [zeta']
  ptb_thta0 = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [theta']

; Add metadata to perturbation arrays 
  ptb_vrad0@description = "Perturbation radial wind"
  ptb_vrad0@units = "m s~S~-1~N~"
  ptb_vrad0!0     = "lev"
  ptb_vrad0!1     = "rad"
  ptb_vrad0!2     = "azi"
  ptb_vrad0!3     = "time"

  ptb_vtan0@description = "Perturbation	tangential wind"
  ptb_vtan0@units = "m s~S~-1~N~"
  ptb_vtan0!0     = "lev"
  ptb_vtan0!1     = "rad"
  ptb_vtan0!2     = "azi"
  ptb_vtan0!3     = "time"

  ptb_vvel0@description = "Perturbation	vertical velocity"
  ptb_vvel0@units = "m s~S~-1~N~"
  ptb_vvel0!0     = "lev"
  ptb_vvel0!1     = "rad"
  ptb_vvel0!2     = "azi"
  ptb_vvel0!3     = "time"

  ptb_pres0@description = "Perturbation pressure"
  ptb_pres0@units = "hPa"
  ptb_pres0!0     = "lev"
  ptb_pres0!1     = "rad"
  ptb_pres0!2     = "azi"
  ptb_pres0!3     = "time"

  ptb_vort0@description = "Perturbation	relative vorticity"
  ptb_vort0@units = "s~S~-1~N~"
  ptb_vort0!0     = "lev"
  ptb_vort0!1     = "rad"
  ptb_vort0!2     = "azi"
  ptb_vort0!3     = "time"

  ptb_thta0@description = "Perturbation potential temperature"
  ptb_thta0@units = "K"
  ptb_thta0!0     = "lev"
  ptb_thta0!1     = "rad"
  ptb_thta0!2     = "azi"
  ptb_thta0!3     = "time"

;===========================================================
; Arrays to hold expanded 'v_int'/'vort_int' arrays (4D)
;===========================================================

  u_int0    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Zonal wind
  v_int0    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Meridional wind
  w_int0    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Vertical velocity
  vort_int0 = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Relative vorticity
  prs_int0  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Pressure 
  thta_int0 = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Potential temp.
  v_rad0    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Radial wind
  v_tan0    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Tangential wind

  u_int0@description = "Zonal velocity"
  u_int0@units       = "m s~S~-1~N~"
  u_int0!0           = "lev"
  u_int0!1           = "rad"
  u_int0!2           = "azi"
  u_int0!3           = "time"

  v_int0@description = "Meridional velocity"
  v_int0@units       = "m s~S~-1~N~"
  v_int0!0           = "lev"
  v_int0!1           = "rad"
  v_int0!2           = "azi"
  v_int0!3           = "time"

  w_int0@description = "Vertical velocity"
  w_int0@units       = "m s~S~-1~N~"
  w_int0!0           = "lev"
  w_int0!1           = "rad"
  w_int0!2           = "azi"
  w_int0!3           = "time"

  prs_int0@description = "Pressure"
  prs_int0@units       = "hPa"
  prs_int0!0           = "lev"
  prs_int0!1           = "rad"
  prs_int0!2           = "azi"
  prs_int0!3           = "time"

  vort_int0@description = "Relative vorticity"
  vort_int0@units       = "s~S~-1~N~"
  vort_int0!0           = "lev"
  vort_int0!1           = "rad"
  vort_int0!2           = "azi"
  vort_int0!3           = "time"

  thta_int0@description = "Potential temperature"
  thta_int0@units       = "K"
  thta_int0!0           = "lev"
  thta_int0!1           = "rad"
  thta_int0!2           = "azi"
  thta_int0!3           = "time"

  v_rad0@description   = "Radial velocity"
  v_rad0@units         = "m s~S~-1~N~"
  v_rad0!0             = "lev"
  v_rad0!1             = "rad"
  v_rad0!2             = "azi"
  v_rad0!3             = "time"

  v_tan0@description   = "Tangential velocity"
  v_tan0@units         = "m s~S~-1~N~"
  v_tan0!0             = "lev"
  v_tan0!1             = "rad"
  v_tan0!2             = "azi"
  v_tan0!3             = "time"

;========================================
; Define 4D arrays for mean quantities 
;========================================

; Mean tangential wind {v}
  mean_vtan0       = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  mean_vtan0!0     = "lev"
  mean_vtan0!1     = "rad"
  mean_vtan0!2     = "azi"
  mean_vtan0!3     = "time"
  mean_vtan0@description = "Azimuthally averaged tangential wind"
  mean_vtan0@units = "m s~S~-1~N~"

; Mean radial wind {u}
  mean_vrad0       = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  mean_vrad0!0     = "lev"
  mean_vrad0!1     = "rad"
  mean_vrad0!2     = "azi"
  mean_vrad0!3     = "time"
  mean_vrad0@description = "Azimuthally averaged radial wind"
  mean_vrad0@units = "m s~S~-1~N~"

; Mean vertical velocity {w}
  mean_vvel0       = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  mean_vvel0!0     = "lev"
  mean_vvel0!1     = "rad"
  mean_vvel0!2     = "azi"
  mean_vvel0!3     = "time"
  mean_vvel0@description = "Azimuthally averaged vertical velocity"
  mean_vvel0@units = "m s~S~-1~N~"

; Mean relative vorticity {zeta}
  mean_vort0       = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  mean_vort0!0     = "lev"
  mean_vort0!1     = "rad"
  mean_vort0!2     = "azi"
  mean_vort0!3     = "time"
  mean_vort0@description = "Azimuthally averaged relative vorticity"
  mean_vort0@units = "s~S~-1~N~"

; Mean potential temperature {theta}
  mean_thta0       = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  mean_thta0!0     = "lev"
  mean_thta0!1     = "rad"
  mean_thta0!2     = "azi"
  mean_thta0!3     = "time"
  mean_thta0@description = "Azimuthally averaged potential temperature"
  mean_thta0@units = "K"

; Mean pressure {p}
  mean_pres0       = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  mean_pres0!0     = "lev"
  mean_pres0!1     = "rad"
  mean_pres0!2	   = "azi"
  mean_pres0!3     = "time"
  mean_pres0@description = "Azimuthally averaged pressure"
  mean_pres0@units = "hPa"

;=========================
; Start multiple loops
;=========================

    print_clock("Working on ensemble member: "+ens0)

    numTIMES    = 120                       ; Total number of times in all files [36]
    ct          = ts0                       ; Counter variable
    ct0		= 0			    ; Counter variable (independent of position)
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
    d  = addfile(fili_d,"r")         	 ; Read in 'pd' stream [w,p,t,theta]

;==================================
; Get the variables we will need
;==================================

   do it = ts0, tf0		 ; Loop over times in file

    print("Working on time: "+time_str(it))
    itt = it
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

  ; Read in temperature, pressure, theta, vertical velocity and geopotential height
  ; All variables below -- [hybD | 63] * [latD | 150] * [lon1D | 151]
    t      = d->temp(it-2,:,{t0:t1},{n0:n1})        ; Temperature (K)
    p      = d->p(it-2,:,{t0:t1},{n0:n1})           ; Pressure (Pa)
    th     = d->theta(it-2,:,{t0:t1},{n0:n1})	    ; Potential temperature (K)
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
     print_clock("Working on model level " + ilev+ " ("+hy+" m AGL)")

   ; Read in other variables on model levels
     u_plane 	            = u(ilev,:,:)	   ; Zonal wind 
     v_plane		    = v(ilev,:,:)	   ; Meridional wind 
     w_plane		    = w(ilev,:,:)	   ; Vertical velocity
     prs_plane              = p(ilev,:,:)	   ; Pressure
     vort_plane	   	    = vort(ilev,:,:)	   ; Relative vorticity
     thta_plane		    = th(ilev,:,:)	   ; Potential temperature

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
     dlat 	   	    = abs(lat_max - centre(it,0))
     dlon	   	    = abs(lon_max - centre(it,1))

   ; If our calculations of latitude and longitude deviate too much from the values
   ; we read in using the NCL function ('st_rm'), then calculate a new centre
   ; using the minimum pressure. Constrain the new calculation using the values
   ; from 'st_rm', to prevent spurious centres being identified and analysed. 

     if (dlat .gt. 0.75 .or. dlon .gt. 0.75) then 

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

     xcen     = new((/dimsizes(radii),dimsizes(thetas)/),float)
     ycen     = new((/dimsizes(radii),dimsizes(thetas)/),float)

     xpos!0   = "rad"
     xpos!1   = "azi"
     xpos&rad = radii
     xpos&azi = thetas

     ypos!0   = "rad"
     ypos!1   = "azi"
     ypos&rad = radii
     ypos&azi = thetas

     xcen!0   = "rad"
     xcen!1   = "azi"
     xcen&rad = radii
     xcen&azi = thetas

     ycen!0   = "rad"
     ycen!1   = "azi"
     ycen&rad = radii
     ycen&azi = thetas

   ; New arrays containing regular variables we want to plot
     u_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Zonal velocity
     v_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Meridional velocity
     w_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vertical velocity
     vort_int = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Relative vorticity
     prs_int  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Pressure 
     thta_int = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Potential temperature
     v_rad    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Radial velocity
     v_tan    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Tangential velocity

   ; Create cylindrincal coordinate array from original latitude/longitude grid
     do irad = 0, dimsizes(radii)-1
       r = radii(irad)
       do iang = 0, dimsizes(thetas)-1
          theta = thetas_rad(iang)   

           xpos(irad,iang) = lon_max + r*cos(theta)
 	   ypos(irad,iang) = lat_max + r*sin(theta)
           xcen(irad,iang) = lon_max
           ycen(irad,iang) = lat_max

       end do
     end do

   ; Fill 'xpos_all' and 'ypos_all' with latitude & longitude values 
     xpos_all(ilev,:,:,it) = xpos(:,:)
     ypos_all(ilev,:,:,it) = ypos(:,:)
     xcen_all(ilev,:,:,it) = xcen(:,:)
     ycen_all(ilev,:,:,it) = ycen(:,:)

   ; Fill new arrays with interpolated values on the cylindrical grid 
   ; Which latitude/longitude grids are each of these variables defined on? 
   ; u_plane     -- latC,  lonC
   ; v_plane     -- latC,  lonC

   ; w_plane     -- latD,  lon1D
   ; prs_plane   -- latD,  lon1D

     do iang = 0, dimsizes(thetas_rad)-1
      u_int(:,iang)    = linint2_points(lonC_sub,latC_sub,u_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      v_int(:,iang)    = linint2_points(lonC_sub,latC_sub,v_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      vort_int(:,iang) = linint2_points(lonC_sub,latC_sub,vort_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      w_int(:,iang)    = linint2_points(lon1D_sub,latD_sub,w_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      prs_int(:,iang)  = linint2_points(lon1D_sub,latD_sub,prs_plane,False,\
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

    ;===================================================
    ; Calculate mean (azimuthally-averaged) variables 
    ;===================================================

      do irad = 0, dimsizes(radii)-1

        mean_vtan0(ilev,irad,:,ct)    = avg(v_tan(irad,:))      ; [ {v} ]
	mean_vrad0(ilev,irad,:,ct)    = avg(v_rad(irad,:))      ; [ {u} ]
	mean_vvel0(ilev,irad,:,ct)    = avg(w_int(irad,:))      ; [ {w} ]
	mean_pres0(ilev,irad,:,ct)    = avg(prs_int(irad,:))    ; [ {p} ]
	mean_thta0(ilev,irad,:,ct)    = avg(thta_int(irad,:))   ; [ {theta} ]
	mean_vort0(ilev,irad,:,ct)    = avg(vort_int(irad,:))   ; [ {zeta} ]

      end do

    ; Read 'v_int' etc into larger arrays
      u_int0(ilev,:,:,ct)    = (/u_int(:,:)/)
      v_int0(ilev,:,:,ct)    = (/v_int(:,:)/)
      w_int0(ilev,:,:,ct)    = (/w_int(:,:)/)
      prs_int0(ilev,:,:,ct)  = (/prs_int(:,:)/)
      vort_int0(ilev,:,:,ct) = (/vort_int(:,:)/)
      thta_int0(ilev,:,:,ct) = (/thta_int(:,:)/)
      v_rad0(ilev,:,:,ct)    = (/v_rad(:,:)/)
      v_tan0(ilev,:,:,ct)    = (/v_tan(:,:)/)

    ; Tidy up
      delete([/u_int,v_int,w_int,v_tan,v_rad/])
      delete([/prs_int,thta_int,vort_int/])
      delete([/u_plane,v_plane,w_plane,vort_plane,thta_plane,prs_plane/])
      delete([/xpos,ypos,xcen,ycen/])
      delete([/latC_plane,lonC_plane,latD_plane,lonD_plane,lon1D_plane,lat1D_plane/])
      delete([/latC_sub,lonC_sub,latD_sub,lonD_sub,lat1D_sub,lon1D_sub/])

     end do			    ; End pressure level loop (do ilev = 0, num_l-1)

    ct  = ct + 1		    ; Counter variable (time)
    ct0 = ct0 + 1		    ; Independent counter variable (time)

  ; Tidy up before next iteration (all variables from 'pc' and 'pd' streams)
    delete([/lonC,latC,lon1C,lat1C,u,v,vort/])
    delete([/lonD,latD,lon1D,lat1D,th,t,p,w,z/])

   end do     ; End time loop (do it = 0, times-1)

;====================================
; Calculate perturbation variables
;====================================

      ptb_vtan0(:,:,:,:) = v_tan0(:,:,:,:)    - mean_vtan0(:,:,:,:)    ; [v']
      ptb_vvel0(:,:,:,:) = w_int0(:,:,:,:)    - mean_vvel0(:,:,:,:)    ; [w']
      ptb_vrad0(:,:,:,:) = v_rad0(:,:,:,:)    - mean_vrad0(:,:,:,:)    ; [u']
      ptb_pres0(:,:,:,:) = prs_int0(:,:,:,:)  - mean_pres0(:,:,:,:)    ; [p']
      ptb_thta0(:,:,:,:) = thta_int0(:,:,:,:)  - mean_thta0(:,:,:,:)   ; [theta']
      ptb_vort0(:,:,:,:) = vort_int0(:,:,:,:)  - mean_vort0(:,:,:,:)   ; [zeta']

;=======================================================
; Remove unncessary 3rd dimension from mean variables
;=======================================================

      mean_vtan1 = mean_vtan0(:,:,0,:)
      mean_vrad1 = mean_vrad0(:,:,0,:)
      mean_vvel1 = mean_vvel0(:,:,0,:)
      mean_vort1 = mean_vort0(:,:,0,:)
      mean_pres1 = mean_pres0(:,:,0,:)
      mean_thta1 = mean_thta0(:,:,0,:)

;=======================================================
; Add metadata to perturbation arrays before plotting
;=======================================================

    ; Read in model level array again
      hybC       = c->hybrid_ht(:)
      hyb_plot   = hybC / 1000
      hyb_plot@units = "km"

    ; Radial wind
      ptb_vrad0&lev  = hyb_plot
      ptb_vrad0&rad  = radii
      ptb_vrad0&azi  = thetas

    ; Tangential wind
      ptb_vtan0&lev  = hyb_plot
      ptb_vtan0&rad  = radii
      ptb_vtan0&azi  = thetas

    ; Vertical velocity
      ptb_vvel0&lev  = hyb_plot
      ptb_vvel0&rad  = radii
      ptb_vvel0&azi  = thetas

    ; Relative vorticity
      ptb_vort0&lev  = hyb_plot
      ptb_vort0&rad  = radii
      ptb_vort0&azi  = thetas

    ; Pressure 
      ptb_pres0&lev  = hyb_plot
      ptb_pres0&rad  = radii
      ptb_pres0&azi  = thetas

    ; Potential temperature
      ptb_thta0&lev  = hyb_plot
      ptb_thta0&rad  = radii
      ptb_thta0&azi  = thetas

;==============================================
; Switch order of dimensions before plotting
;==============================================

    ; 'gsn_csm_contour' plots the leftmost (rightmost) dimension on the y-axis (x-axis)
    ; Switch 'rad' with 'time' so that 'rad' is the rightmost dimension
    ; New arrays: (lev * time * azi * rad)

    ; Radial wind 
      ptb_vrad        = (/ptb_vrad0(lev|:,time|:,azi|:,rad|:)/)

    ; Tangential wind
      ptb_vtan        = (/ptb_vtan0(lev|:,time|:,azi|:,rad|:)/)

    ; Vertical velocity
      ptb_vvel        = (/ptb_vvel0(lev|:,time|:,azi|:,rad|:)/)

    ; Relative vorticity
      ptb_vort        = (/ptb_vort0(lev|:,time|:,azi|:,rad|:)/)

    ; Pressure
      ptb_pres        = (/ptb_pres0(lev|:,time|:,azi|:,rad|:)/)

    ; Potential temperature
      ptb_thta        = (/ptb_thta0(lev|:,time|:,azi|:,rad|:)/)

;============================================
; Define output file path and plot details 
;============================================

  do it = ts0, tf0

   print("Working on time: "+it)

   do ilev = ls0, lf0

    hy = sprintf("%0.0f",hyb_plot(ilev))
    print("Working on model level " + ilev+ " ("+hy+" m AGL)")

  ; Output file location and type
    output = "$sam/nepartak/images/cyl_coords/mbgt_ptb_xy_"+ilev+"_"+it
    wks = gsn_open_wks(opt,output)

  ; Set up panel plot 
;    panel = new(6,graphic)
    gsn_define_colormap(wks,"BlueDarkRed18")

;==========================
; Options for plotting 
;==========================

  ; Momentum budget terms 
    opts_mbgt                              = True
    opts_mbgt@cnFillOn                     = True
    opts_mbgt@sfXArray			   = xpos_all(ilev,:,:,it)
    opts_mbgt@sfYArray			   = ypos_all(ilev,:,:,it)
    opts_mbgt@cnLineLabelInterval          = 2.0
    opts_mbgt@cnLineLabelFontHeightF       = 0.012
    opts_mbgt@cnLineLabelBackgroundColor   = "transparent"
    opts_mbgt@cnLineLabelPlacementMode     = "constant"
    opts_mbgt@cnLinesOn                    = False
    opts_mbgt@cnInfoLabelOn                = False
    opts_mbgt@cnLevelSelectionMode         = "ExplicitLevels"
    opts_mbgt@cnFillColors                 = (/2,4,5,6,7,8,9,0,0,12,13,14,15,16,17,18/)
    opts_mbgt@cnLevels  	           = (/-10., -5., -3., -2., -1., -0.5, -0.2, \
                                                0, 0.2, 0.5, 1., 2., 3., 5., 10./)
    opts_mbgt@gsnPaperOrientation          = "landscape"
    opts_mbgt@tiMainString                 = ""
    opts_mbgt@tiMainFontHeightF            = 0.0125
    opts_mbgt@gsnLeftString                = ""
    opts_mbgt@gsnRightString               = ""

  ; Additional plotting resources
    opts_mbgt@mpDataBaseVersion            = "Ncarg4_1" ; More recent database
    opts_mbgt@mpDataSetName                = "Earth..4" ; High resolution
    opts_mbgt@mpOutlineBoundarySets        = "National" ; National borders
    opts_mbgt@mpGeophysicalLineColor       = "black"    ; Colour borders black
    opts_mbgt@mpGeophysicalLineThicknessF  = 1.0        ; Border line thickness
    opts_mbgt@mpGridAndLimbOn              = False      ; Turn on lat/lon lines
    opts_mbgt@pmTickMarkDisplayMode        = "Never"    ; Turn on map tickmarks
    opts_mbgt@tmXBMajorLengthF             = 0.005      ; Change tickmark length
    opts_mbgt@tmXTOn                       = "False"    ; No tickmarks on top x-axis
    opts_mbgt@tmYROn                       = "False"    ; No tickmarks on right y-axis
    opts_mbgt@gsnMaximize                  = False      ; Maximise plot size
    opts_mbgt@gsnAddCyclic                 = False
    opts_mbgt@gsnDraw                      = False
    opts_mbgt@gsnFrame                     = False
    opts_mbgt@mpLimitMode                  = "Corners"

    radius                                 = radii(rad_size) ; Size of plot (lat/lon)

    opts_mbgt@mpLeftCornerLatF             = ypos_all(ilev,0,0,it)-radius
    opts_mbgt@mpLeftCornerLonF             = xpos_all(ilev,0,0,it)-radius
    opts_mbgt@mpRightCornerLatF            = ypos_all(ilev,0,0,it)+radius
    opts_mbgt@mpRightCornerLonF            = xpos_all(ilev,0,0,it)+radius

;===============================================================
; Resources to overlay panel labels without using 'gsn_panel'
;===============================================================

    txid			= new(6,graphic)
    amid	    	        = new(6,graphic)

  ; Label text resources 
    txres                       = True
    txres@txPerimOn             = True
    txres@txFontHeightF         = 0.010
    txres@txBackgroundFillColor = "White"

  ; Label position resources 
    amres 	   	        = True
    amres@amParallelPosF	= 0.5			; Right edge
    amres@amOrthogonalPosF	= -0.5			; Top edge
    amres@amJust		= "TopRight"

    panel_strings		= (/"a) vrad", "b) vtan", "c) vvel", \
                                    "d) vort", "e) pres", "f) theta"/)

;=================
; Plot the data 
;=================

;    panel_vrad  = gsn_csm_contour_map(wks,ptb_vrad0(ilev,:,:,it),opts_mbgt) ; [u]
;    panel_vtan  = gsn_csm_contour_map(wks,ptb_vtan0(ilev,:,:,it),opts_mbgt) ; [v]
;    panel_vvel  = gsn_csm_contour_map(wks,ptb_vvel0(ilev,:,:,it),opts_mbgt) ; [w]
;    panel_vort  = gsn_csm_contour_map(wks,ptb_vort0(ilev,:,:,it),opts_mbgt) ; [zeta]
;    panel_pres  = gsn_csm_contour_map(wks,ptb_pres0(ilev,:,:,it),opts_mbgt) ; [p]
;    panel_thta  = gsn_csm_contour_map(wks,ptb_thta0(ilev,:,:,it),opts_mbgt) ; [theta]

;    panel_vrad  = gsn_csm_contour(wks,ptb_vrad(ilev,:,:,it),opts_mbgt) ; [u]
;    panel_vtan  = gsn_csm_contour(wks,ptb_vtan(ilev,:,:,it),opts_mbgt) ; [v]
;    panel_vvel  = gsn_csm_contour(wks,ptb_vvel(ilev,:,:,it),opts_mbgt) ; [w]
;    panel_vort  = gsn_csm_contour(wks,ptb_vort(ilev,:,:,it),opts_mbgt) ; [zeta]
;    panel_pres  = gsn_csm_contour(wks,ptb_pres(ilev,:,:,it),opts_mbgt) ; [p]
;    panel_thta  = gsn_csm_contour(wks,ptb_thta(ilev,:,:,it),opts_mbgt) ; [theta]

;    panel(0)	= panel_vrad
;    panel(1)	= panel_vtan
;    panel(2)	= panel_vvel
;    panel(3)	= panel_vort
;    panel(4)	= panel_pres
;    panel(5)	= panel_thta

;===================================================
; Finally, draw the plot with everything overlaid
;===================================================

;    optsP                       = True
;    optsP@gsnFrame              = False    ; Do not advance the frame
;    optsP@gsnPanelLabelBar      = True     ; Panel labelbar

;    optsP@pmLabelBarWidthF      = 0.6      ; Label bar width
;    optsP@pmLabelBarHeightF     = 0.15     ; Label bar height
;    optsP@lbLabelFontHeightF    = 0.0125   ; Label bar font height
;    optsP@lbLabelFont           = "Helvetica"
;    optsP@lbPerimOn             = False
;    optsP@lbBoxEndCapStyle      = "TriangleBothEnds"

;    optsP@txString              = ""       ; Figure title
;    optsP@gsnPanelFigureStrings = (/"a) vtan", "b) vrad", "c) vvel", \
;                                    "d) vort", "e) pres", "f) theta"/)
;    optsP@gsnMaximize           = True
;    optsP@gsnPanelTop           = 0.90
;    optsP@gsnPanelBottom        = 0.10
;    optsP@amJust                = "TopLeft"
;    optsP@gsnPanelFigureStringsFontHeightF = 0.0125 ; Label size (default 0.01)

;    gsn_panel(wks,panel,(/2,3/),optsP)             ; Draw as a single plot
;    frame(wks)

;   end do                                                  ; End model level loop

;  end do                                                   ; End time loop

;end

;===================================================
;===================================================
;===================================================

   ; Radial wind
     opts_vrad2                             = opts_mbgt
     opts_vrad2@lbLabelBarOn                = True
     opts_vrad2@lbBoxEndCapStyle            = "TriangleBothEnds"
     opts_vrad2@lbOrientation               = "Vertical"
     opts_vrad2@vpXF                        = 0.14
     opts_vrad2@vpYF                        = 0.80
     opts_vrad2@vpWidthF                    = 0.20
     opts_vrad2@vpHeightF                   = 0.20 
     panel_vrad  = gsn_csm_contour_map(wks,ptb_vrad0(ilev,:,:,it),opts_vrad2) ; [u]
     txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
     amid(0)     = gsn_add_annotation(panel_vrad, txid(0), amres)

   ; Tangential wind
     opts_vtan2                             = opts_mbgt
     opts_vtan2@lbLabelBarOn                = True
     opts_vtan2@lbBoxEndCapStyle            = "TriangleBothEnds"
     opts_vtan2@lbOrientation               = "Vertical"
     opts_vtan2@vpXF                        = 0.44
     opts_vtan2@vpYF                        = 0.80
     opts_vtan2@vpWidthF                    = 0.20
     opts_vtan2@vpHeightF                   = 0.20
     panel_vtan  = gsn_csm_contour_map(wks,ptb_vtan0(ilev,:,:,it),opts_vtan2) ; [v]
     txid(1)     = gsn_create_text(wks, panel_strings(1), txres)
     amid(1)     = gsn_add_annotation(panel_vtan, txid(1), amres)

   ; Vertical velocity
     opts_vvel2                             = opts_mbgt
     opts_vvel2@lbLabelBarOn                = True
     opts_vvel2@lbBoxEndCapStyle            = "TriangleBothEnds"
     opts_vvel2@lbOrientation               = "Vertical"
     opts_vvel2@vpXF                        = 0.74
     opts_vvel2@vpYF                        = 0.80
     opts_vvel2@vpWidthF                    = 0.20
     opts_vvel2@vpHeightF                   = 0.20
     panel_vvel  = gsn_csm_contour_map(wks,ptb_vvel0(ilev,:,:,it),opts_vvel2) ; [w]
     txid(2)     = gsn_create_text(wks, panel_strings(2), txres)
     amid(2)     = gsn_add_annotation(panel_vvel, txid(2), amres)

   ; Relative vorticity
     opts_vort2                             = opts_mbgt
     opts_vort2@lbLabelBarOn                = True
     opts_vort2@lbBoxEndCapStyle            = "TriangleBothEnds"
     opts_vort2@lbOrientation               = "Vertical"
     opts_vort2@vpXF                        = 0.14
     opts_vort2@vpYF                        = 0.50
     opts_vort2@vpWidthF                    = 0.20
     opts_vort2@vpHeightF                   = 0.20
     panel_vort  = gsn_csm_contour_map(wks,ptb_vort0(ilev,:,:,it),opts_vort2) ; [zeta]
     txid(3)     = gsn_create_text(wks, panel_strings(3), txres)
     amid(3)     = gsn_add_annotation(panel_vort, txid(3), amres)

   ; Pressure
     opts_pres2                             = opts_mbgt
     opts_pres2@lbLabelBarOn                = True
     opts_pres2@lbBoxEndCapStyle            = "TriangleBothEnds"
     opts_pres2@lbOrientation               = "Vertical"
     opts_pres2@vpXF                        = 0.44
     opts_pres2@vpYF                        = 0.50
     opts_pres2@vpWidthF                    = 0.20
     opts_pres2@vpHeightF                   = 0.20
     panel_pres  = gsn_csm_contour_map(wks,ptb_pres0(ilev,:,:,it),opts_pres2) ; [p]
     txid(4)     = gsn_create_text(wks, panel_strings(4), txres)
     amid(4)     = gsn_add_annotation(panel_pres, txid(4), amres)

   ; Potential temperature
     opts_thta2                             = opts_mbgt
     opts_thta2@lbLabelBarOn                = True
     opts_thta2@lbBoxEndCapStyle            = "TriangleBothEnds"
     opts_thta2@lbOrientation               = "Vertical"
     opts_thta2@vpXF                        = 0.74
     opts_thta2@vpYF                        = 0.50
     opts_thta2@vpWidthF                    = 0.20
     opts_thta2@vpHeightF                   = 0.20
     panel_thta  = gsn_csm_contour_map(wks,ptb_thta0(ilev,:,:,it),opts_thta2) ; [theta]
     txid(5)     = gsn_create_text(wks, panel_strings(5), txres)
     amid(5)     = gsn_add_annotation(panel_thta, txid(5), amres)

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

     pres			= True 
     maximize_output(wks,pres)				   ; Calls 'draw' and 'frame'

    end do						   ; End model level loop

   end do						   ; End time loop

end 