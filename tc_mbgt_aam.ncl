; Script to calculate the radial momentum budget as in Montgomery et al. (2018)

; Run using:

; ncl 'opt="pdf"' 'ens0="em11"' dist=3.0 ts0=60 tf0=61 calc=1 rad0=1 thta=1 tc_mbgt_aam.ncl

; 'opt'   = output file format ("pdf" or "x11")
; 'ens0'  = ensemble member (em00 -- em11)
; 'dist'  = size of box following storm (degrees)
; 'ts0'   = starting time index (0–119)
; 'tf0'   = ending time index (0–119)
; 'calc'  = centre calculated on each level (0) or equal to 950 hPa centre (1)
; 'rad0'  = radius of final plots: 150 km (0) or 200 km (1)
; 'thta'  = plot either u,v,w,AAM (0) or AAM and theta-e (1)

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
  fili_prefix = "20160702T1200Z_NPTK_4p4_L80_ra1t_"+ens0
  fili_b1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pb*.nc")
  fili_b      = diri+"/"+fili_b1
  fili_c1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pc*.nc")
  fili_c      = diri+"/"+fili_c1
  fili_d1     =	systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pd*.nc")
  fili_d      = diri+"/"+fili_d1

; Calculate number of times in each file using built-in function ('nc_times')
  numINFO_b   = nc_times(fili_b)
  numINFO_c   = nc_times(fili_c)
  numINFO_d   = nc_times(fili_d)
  numINPUT_b  = numINFO_b[0]		; Number of times in 'pb' stream
  time_b      = numINFO_b[1]		; Array of times from 'pb' stream file
  numINPUT_c  = numINFO_c[0]		; Number of times in 'pc' stream
  time_c      = numINFO_c[1]		; Array of times from 'pc' stream file
  numINPUT_d  = numINFO_d[0]		; Number of times in 'pd' stream
  time_d      = numINFO_d[1]		; Array of times from 'pd' stream file 

; Calculate storm motion using built-in function ('st_rm')
; Argument 1  = ensemble simulation
; Argument 2  = path to input files
; Argument 3  = size of box following storm (degrees)

; Call external function 'st_rm'
  storm_rel = st_rm(ens0,"$ar/text/",3.0)
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
  num_p	     = 17		     ; Number of pressure levels 

; Define constants
  rd	= 287.0	  ; Gas constant for dry air (J/kg/K^2)
  cp0	= 1004.0  ; Specific heat of dry air at constant pressure (J/kg/K)
  kp	= 0.286   ; For Exner function calculations (0.286)
  pref	= 1000.0  ; Reference pressure p0 (1000 hPa)

;========================================================================================
;========================================================================================

; Array to hold all 'xpos' and 'ypos' position arrays
  xpos_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  ypos_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)

; Array for mean absolute angular momentum (3D)
; (1) {v} * r + (0.5 * ( {f} * r ^ 2 ) )
  mean_aam0  = new((/num_l,dimsizes(radii),numTIMES/),float) ; (1)

;========================================
; Define 3D arrays for mean quantities 
;========================================

; Mean tangential wind {v}
  mean_vtan0       = new((/num_l,dimsizes(radii),numTIMES/),float)
  mean_vtan0!0     = "lev"
  mean_vtan0!1     = "rad"
  mean_vtan0!2     = "time"
  mean_vtan0@description = "Azimuthally averaged tangential wind"
  mean_vtan0@units = "m s~S~-1~N~"

; Mean radial wind {u}
  mean_vrad0       = new((/num_l,dimsizes(radii),numTIMES/),float)
  mean_vrad0!0     = "lev"
  mean_vrad0!1     = "rad"
  mean_vrad0!2     = "time"
  mean_vrad0@description = "Azimuthally averaged radial wind"
  mean_vrad0@units = "m s~S~-1~N~"

; Mean vertical velocity {w}
  mean_vvel0       = new((/num_l,dimsizes(radii),numTIMES/),float)
  mean_vvel0!0     = "lev"
  mean_vvel0!1     = "rad"
  mean_vvel0!2     = "time"
  mean_vvel0@description = "Azimuthally averaged vertical velocity"
  mean_vvel0@units = "m s~S~-1~N~"

; Mean absolute angular momentum {w}
  mean_aam0        = new((/num_l,dimsizes(radii),numTIMES/),float)
  mean_aam0!0      = "lev"
  mean_aam0!1      = "rad"
  mean_aam0!2      = "time"
  mean_aam0@description = "Azimuthally averaged absolute angular momentum"
  mean_aam0@units  = "m s~S~-1~N~"

; Mean Coriolis force {f0}
  mean_cor0        = new((/num_l,dimsizes(radii),numTIMES/),float)
  mean_cor0!0      = "lev"
  mean_cor0!1      = "rad"
  mean_cor0!2      = "time"
  mean_cor0@description = "Azimuthally averaged Coriolis force"
  mean_cor0@units  = "s~S~-1~N~"

; Mean theta-e {the}
  mean_the0        = new((/num_p,dimsizes(radii),numTIMES/),float)
  mean_the0!0      = "lev"
  mean_the0!1      = "rad"
  mean_the0!2      = "time"
  mean_the0@description = "Azimuthally averaged equivalent potential temperature"
  mean_the0@units  = "K"

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

    b  = addfile(fili_b,"r")		 ; Read in 'pb' stream []
    c  = addfile(fili_c,"r")         	 ; Read in 'pc' stream [u,v,vort]
    d  = addfile(fili_d,"r")         	 ; Read in 'pd' stream [w,p,t,theta,kmh,tau]

;==================================
; Get the variables we will need
;==================================

   do it = ts0, tf0		; Loop over times in file

    print("Working on time: "+time_str(it))
    itt = it
    title_arr(it) = "Valid at "+time_str(it)+" (T+"+itt+")"

  ; Stream-specific time array indices ('pb', 'pc', 'pd' streams)
    bt = it - 4
    dt = it - 2

  ; Stream-specific indices for accessing storm track data 
    bt0 = it - 3 
    ct0	= it + 1
    dt0	= it - 1

;========================================================
; Read in variables from 'pb' stream (pressure levels)
;========================================================

    t0    = lat_arr(it,0)
    t1    = lat_arr(it,1)
    n0    = lon_arr(it,0)
    n1    = lon_arr(it,1)
    print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

 ; Create latitude/longitude array subsets to reduce computation time
   lon1B  = b->longitude_1({n0:n1})  ; longitude [160 grid points]
   lat1B  = b->latitude_1({t0:t1})   ; latitude  [160 grid points]
   lonB   = b->longitude({n0:n1})    ; longitude [160 grid points]
   latB   = b->latitude({t0:t1})     ; latitude  [160 grid points]

 ; Pressure levels (pres)
 ; (0) 1000, (1) 950, (2) 925, (3) 900, (4) 850, (5) 800, (6) 750, (7) 700
 ; (8) 650, (9) 600, (10) 500, (11) 400, (12) 300, (13) 250
 ; (14) 200, (15) 150, (16) 100

   pres  = b->p	  				; Subset of pressure level array
   pres@units = "hPa"                 	        ; Pressure level units (for plotting)

 ; Arrays have () vertical levels, and horizontal extent depends on 'dist' option
   t00    = b->temp(bt,:,{t0:t1},{n0:n1})    ; Temperature (K)
   rh0	  = b->rh(bt,:,{t0:t1},{n0:n1})      ; Relative humidity (%)

 ; Avoid dimension size errors throughout the script
   tsize  = dimsizes(t00)
   rhsize = dimsizes(rh0)

   if ( any (tsize .gt. dsize(0)) ) then
    t = t00(0:dsize(0)-1,0:dsize(0)-1)
   else
    t = t00
   end if

   if ( any (rhsize .gt. dsize(0)) ) then
    rh   = rh0(0:dsize(0)-1,0:dsize(0)-1)
   else
    rh   = rh0
   end if

; ; Horizontal grid spacing
;   dx    = lonB(1) - lonB(0) ; Grid spacing (longitude)
;   dy    = latB(1) - latB(0) ; Grid spacing (latitude)

 ;=========================================
 ; Calculate theta-e using NCL functions
 ;=========================================

 ; First calculate mixing ratio from pressure, temperature and RH
 ; Option '1' returns the mixing ratio; '2' returns the specific humidity

 ; Create pressure array with same dimensions as 't' for mixing ratio calculation
   tsize  = dimsizes(t)
   pres0  = conform_dims(tsize,pres,0)
   delete(tsize)

 ; Now calculate the mixing ratio from pressure, temperature and RH
   mixr   = mixhum_ptrh(pres0,t,rh,1)

 ; Then calculate the lifting condensation level temperature using RH
 ; Option '(/1,1/)' specifies units of input and output temperatures (K)
   tlcl   = tlcl_rh_bolton(t,rh,(/1,1/))

 ; Finally calculate equivalent potential temperature
 ; Option '(/0,1,0,1/)' -- input p (hPa), temp (K), mixing ratio (kg/kg), theta-e (K)
   th_e   = pot_temp_equiv_tlcl(pres0,t,tlcl,mixr,(/0,1,0,1/))

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
    td     = d->temp(dt,:,{t0:t1},{n0:n1})          ; Temperature (K)
    p      = d->p(dt,:,{t0:t1},{n0:n1})             ; Pressure (Pa)
    th     = d->theta(dt,:,{t0:t1},{n0:n1})	    ; Potential temperature (K)
    w      = d->dz_dt(dt,:,{t0:t1},{n0:n1})	    ; Vertical velocity (m s-1)
    z      = d->ht(dt,:,{t0:t1},{n0:n1})    	    ; Geopotential height (m)

  ; Calculate dry air density from pressure and temperature
    rho	   = p / (rd * td)

  ; Add metadata
    copy_VarCoords(p, rho)
    rho@description = "Dry air density"
    rho@units = "kg m~S~-3~N~"

  ; Also calculate specific volume (for	use in diffusive tendency source terms)
    rho_a  = 1 / rho
    copy_VarCoords(rho, rho_a)
    rho_a@description = "Specific volume"
    rho_a@units = "m~S~3~N~ kg~S~-1~N~"

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

   ; Read in other variables on model levels
     u_plane 	            = u(ilev,:,:)	   ; Zonal wind 
     v_plane		    = v(ilev,:,:)	   ; Meridional wind 
     w_plane		    = w(ilev,:,:)	   ; Vertical velocity
     f_plane                = f0(ilev,:,:)         ; Coriolis parameter
     prs_plane              = p(ilev,:,:)	   ; Pressure
     vort_plane	   	    = vort(ilev,:,:)	   ; Relative vorticity

   ; Read in theta-e on pressure levels
     if (ilev .le. num_p-1) then 
      print("Also working on pressure level: " + pres(ilev)+ " hPa")
      the_plane		    = th_e(ilev,:,:)	   ; Theta-e
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
   ; 'pb' stream
   ;===============

   ; Indices of lat/lon grid for 'str_yz' in 'pd' stream
     lat1B_plane = ind(t0.le.lat1B.and.lat1B.le.t1)
     lon1B_plane = ind(n0.le.lon1B.and.lon1B.le.n1)

   ; Array sizes
     lat1B_size  = dimsizes(lat1B_plane)
     lon1B_size  = dimsizes(lon1B_plane)

   ; Array of latitude/longitude points corresponding to the indices above
     lat1B_sub   = lat1B(lat1B_plane(0):lat1B_plane(lat1B_size-1))
     lon1B_sub   = lon1B(lon1B_plane(0):lon1B_plane(lon1B_size-1))

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
     f_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Coriolis force
     v_rad    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Radial velocity
     v_tan    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Tangential velocity

     if (ilev .le. num_p-1) then 
      the_int = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Theta-e
     end if 

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

     copy_VarAtts(f_plane,f_int)
     f_int@description = "Coriolis parameter"
     f_int@units       = "s~S~-1~N~"

     if (ilev .le. num_p-1) then 
      copy_VarAtts(the_plane,the_int)
      the_int@description = "Equivalent potential temperature"
      the_int@units       = "K"
     end if 

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
     xpos_all(ilev,:,:,it) = xpos(:,:)
     ypos_all(ilev,:,:,it) = ypos(:,:)

   ; Fill new arrays with interpolated values on the cylindrical grid 
   ; Which latitude/longitude grids are each of these variables defined on? 
   ; u_plane     -- latC,  lonC
   ; v_plane     -- latC,  lonC
   ; f_plane     -- latC,  lonC
   ; w_plane     -- latD,  lon1D
   ; the_plane	 -- lat1B, lon1B

     do iang = 0, dimsizes(thetas_rad)-1
      u_int(:,iang)    = linint2_points(lonC_sub,latC_sub,u_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      v_int(:,iang)    = linint2_points(lonC_sub,latC_sub,v_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      w_int(:,iang)    = linint2_points(lon1D_sub,latD_sub,w_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      f_int(:,iang)    = linint2_points(lonC_sub,latC_sub,f_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      if (ilev .le. num_p-1) then
       the_int(:,iang) = linint2_points(lon1B_sub,lat1B_sub,the_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      end if 

     end do

   ; Calculate radial/tangential windspeed and friction at each point (using unit vctrs)
      do iang = 0, dimsizes(thetas)-1
       theta = thetas_rad(iang)

      ; Tangential and radial wind 
        v_rad(:,iang) = u_int(:,iang)*cos(theta) + v_int(:,iang)*sin(theta)
        v_tan(:,iang) = -u_int(:,iang)*sin(theta) + v_int(:,iang)*cos(theta)

      end do

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

     w_int!0      = "rad"
     w_int!1      = "azi"
     w_int&rad    = radii
     w_int&azi    = thetas

     f_int!0      = "rad"
     f_int!1      = "azi"
     f_int&rad    = radii
     f_int&azi    = thetas

     if (ilev .le. num_p-1) then 
      the_int!0    = "rad"
      the_int!1    = "azi"
      the_int&rad  = radii
      the_int&azi  = thetas
     end if 

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

        mean_vtan0(ilev,irad,ct)    = avg(v_tan(irad,:))      ; [ {v} ]
	mean_vrad0(ilev,irad,ct)    = avg(v_rad(irad,:))      ; [ {u} ]
	mean_vvel0(ilev,irad,ct)    = avg(w_int(irad,:))      ; [ {w} ]
	mean_cor0(ilev,irad,ct)	    = avg(f_int(irad,:))      ; [ {f} ]

	if (ilev .le. num_p-1) then 
 	 mean_the0(ilev,irad,ct)    = avg(the_int(irad,:))    ; [ {the} ]
	end if 

      end do

    ;============================================
    ; Calculate mean absolute angular momentum 
    ;============================================
      
       do irad = 1, dimsizes(radii)-1

        rp1 = min((/irad+1,dimsizes(radii)-1/))
        rm1 = max((/irad-1,0/))

       ; Calculate the radial distance ('dr') between all points in array
       ; Arguments 1 and 2 represent the lat/lon of the first point, respectively
       ; Arguments 3 and 4 represent the lat/lon of the second point, respectively
       ; Argument 5 ('2') tells the script to avoid extra interpolated points
       ; Argument 6 ('3') returns the distance in metres

       ; Why were the latitude and longitude switched in the azimuthal derivative loop?
       ; Make sure the values I'm outputting here are correct, before continuing.

         th0 = (dimsizes(thetas) / 2)

       ; Radial distance between each point and its nearest neighbour
         dr  = gc_latlon(ypos(rm1,th0),xpos(rm1,th0),\
                         ypos(rp1,th0),xpos(rp1,th0),2,3)

       ; Radial distance between each point and the origin
         r1  = gc_latlon(ypos(0,th0),xpos(0,th0),\
                         ypos(irad,th0),xpos(irad,th0),2,3)

       ; Calculate increments
         vm   = mean_vtan0(ilev,irad,ct)
	 fm   = mean_cor0(ilev,irad,ct) * (r1 ^ 2)

       ; Also calculate absolute angular momentum components
	 aam0 = vm * r1 
	 aam1 = 0.5 * fm

       ; Calculate absolute angular momentum [ rv + 0.5 * (fr ^ 2) ]
         mean_aam0(ilev,irad,ct)  = aam0 + aam1

       end do                         ; End loop over radii

    ; Tidy up
      delete([/u_int,v_int,w_int,v_tan,v_rad,xpos,ypos/])
      delete([/f_plane,u_plane,v_plane,w_plane,vort_plane,prs_plane/])
      delete([/latC_plane,lonC_plane,latD_plane,lonD_plane,lon1D_plane,lat1D_plane/])
      delete([/latC_sub,lonC_sub,latD_sub,lonD_sub,lat1D_sub,lon1D_sub/])

      if (ilev .le. num_p-1) then 
       delete(the_plane)
      end if 

     end do			    ; End pressure level loop (do ilev = 0, num_l-1)

    ct  = ct + 1		    ; Counter variable (time)
    ct0 = ct0 + 1		    ; Independent counter variable (time)

  ; Tidy up before next iteration (all variables from 'pb', 'pc' and 'pd' streams)
    delete([/lat1B,lon1B,latB,lonB,t,rh,th_e,t00,rh0,pres0,mixr,tlcl/])
    delete([/lonC,latC,lon1C,lat1C,u,v,vort,avo,f0,f1/])
    delete([/lonD,latD,lon1D,lat1D,td,p,th,w,z,rho,rho_a/])

   end do     ; End time loop (do it = 0, times-1)

;=================================================
; Tidy up and produce final arrays for plotting 
;=================================================

     ; 'gsn_csm_contour' plots the leftmost (rightmost) dimension on the y-axis (x-axis)
     ; Switch 'rad' with 'time' so that 'rad' is the rightmost dimension

       hyb_plot   = hybC / 1000
       hyb_plot@units = "km"

     ; Mean radial wind -- also add metadata
       mean_vrad       = mean_vrad0(lev|:,time|:,rad|:)
       mean_vrad!0     = "lev"
       mean_vrad&lev   = hyb_plot
       mean_vrad!1     = "time"
       mean_vrad!2     = "rad"
       mean_vrad&rad   = radii
       mean_vrad@units = "m s~S~-1~N~"

     ; Mean tangential wind -- also add metadata
       mean_vtan       = mean_vtan0(lev|:,time|:,rad|:)
       mean_vtan!0     = "lev"
       mean_vtan&lev   = hyb_plot
       mean_vtan!1     = "time"
       mean_vtan!2     = "rad"
       mean_vtan&rad   = radii
       mean_vtan@units = "m s~S~-1~N~"

     ; Mean vertical velocity -- also add metadata
       mean_vvel       = mean_vvel0(lev|:,time|:,rad|:)
       mean_vvel!0     = "lev"
       mean_vvel&lev   = hyb_plot
       mean_vvel!1     = "time"
       mean_vvel!2     = "rad"
       mean_vvel&rad   = radii
       mean_vvel@units = "m s~S~-1~N~"

     ; Mean absolute angular momentum
       mean_aam        = mean_aam0(lev|:,time|:,rad|:)
       mean_aam!0      = "lev"
       mean_aam&lev    = hyb_plot
       mean_aam!1      = "time"
       mean_aam!2      = "rad"
       mean_aam&rad    = radii
       mean_aam	       = mean_aam / (10 ^ 6)
       mean_aam@units  = "10~S~6~N~ m~S~2~N~ s~S~-1~N~"

     ; Mean equivalent potential temperature (theta-e)
       mean_the	       = mean_the0(lev|:,time|:,rad|:)
       mean_the!0      = "lev"
       mean_the&lev    = pres
       mean_the!1      = "time"
       mean_the!2      = "rad"
       mean_the&rad    = radii

;============================================
; Define output file path and plot details 
;============================================

 do it = ts0, tf0

    print("Working on time: "+it)

  ; Output file location and type
    output = "$sam/nepartak/images/cyl_coords/aam_the_xz_"+ens0+"_"+time_arr(it)
    wks = gsn_open_wks(opt,output)

  ; Set up panel plot
    if (thta .eq. 0) then 
     panel = new(5,graphic)
    elseif (thta .eq. 1) then
     panel = new(1,graphic)
    end if

;==========================
; Options for plotting 
;==========================

  ; Momentum budget terms 
    opts_mbgt                              = True
    opts_mbgt@cnFillOn                     = True
    opts_mbgt@cnLineLabelInterval          = 2.0
    opts_mbgt@cnLineLabelFontHeightF       = 0.012
    opts_mbgt@cnLineLabelBackgroundColor   = "transparent"
    opts_mbgt@cnLineLabelPlacementMode     = "constant"
    opts_mbgt@cnLinesOn                    = False
    opts_mbgt@cnInfoLabelOn                = False
    opts_mbgt@cnLevelSelectionMode         = "ExplicitLevels"

    opts_mbgt@gsnPaperOrientation          = "landscape"
    opts_mbgt@tiMainString                 = ""
    opts_mbgt@tiMainFontHeightF            = 0.0125
    opts_mbgt@gsnLeftString                = ""
    opts_mbgt@gsnRightString               = ""
    opts_mbgt@gsnMaximize                  = True
    opts_mbgt@lbLabelBarOn                 = False
    opts_mbgt@gsnDraw                      = False
    opts_mbgt@gsnFrame                     = False

  ; Additional plotting resources
    opts_mbgt@tiYAxisString                = "Height (km)"
    opts_mbgt@tiXAxisString                = "Radial distance from cyclone centre (km)"
    opts_mbgt@tiXAxisFontHeightF           = 0.015   ; X-axis title font height 
;    opts_mbgt@tiXAxisOffsetYF      	   = 0.175   ; x-axis title up/down (h/l)

    opts_mbgt@tmXTOn                       = "False" ; Turn off top x-axis TM
    opts_mbgt@tmYROn                       = "False" ; Turn off right y-axis TM
    opts_mbgt@gsnMaximize                  = True    ; Maximise plot size
    opts_mbgt@gsnAddCyclic                 = False
    opts_mbgt@tmXBMode                     = "Explicit"  ; Set tick marks explicitly
    opts_mbgt@tmXBLabelFontHeightF	   = 0.015    ; Label font height 

    if (rad0 .eq. 0) then
     opts_mbgt@tmXBValues                   = (/0.0,0.25,0.5,0.75,1.0,1.25,1.5/)
     opts_mbgt@tmXBLabels                   = (/"0","25","50","75","100","125","150"/)
    else if (rad0 .eq. 1) then
     opts_mbgt@tmXBValues                   = (/0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0/)
     opts_mbgt@tmXBLabels                   = (/"0","25","50","75","100","125","150","175","200"/)
    end if
    end if

    opts_mbgt@tmYLMode			   = "Explicit"
    opts_mbgt@tmYLValues		   = fspan(2.0, 18.0, 9)
    opts_mbgt@tmYLLabels		   = (/"2.0","4.0","6.0","8.0","10.0",\
    					       "12.0","14.0","16.0","18.0"/)
    opts_mbgt@tmYLLabelFontHeightF	   = 0.015
    opts_mbgt@gsnYAxisIrregular2Linear	   = True	; Linearise the y-axis 

    opts_mbgt@trYMinF			   = 0.0	; y-axis minimum value
    opts_mbgt@trYMaxF			   = 17.6	; y-axis maximum value 

  ; Resources to overlay radial wind (following M18)
    opts_vrad				   = True
    opts_vrad@cnFillOn			   = False
    opts_vrad@cnLineColor		   = "black"
    opts_vrad@cnLevelSelectionMode	   = "ExplicitLevels"
    opts_vrad@cnLevels                     = (/-1.2, 1.2/)
    opts_vrad@gsnContourNegLineDashPattern = 16
    opts_vrad@gsnContourPosLineDashPattern = 0
    opts_vrad@cnLineThicknessF             = 2.5

    opts_vrad@cnInfoLabelOn                = False
    opts_vrad@gsnDraw			   = False
    opts_vrad@gsnFrame			   = False
    opts_vrad@cnLineLabelsOn		   = False
    opts_vrad@tiMainString                 = ""
    opts_vrad@gsnLeftString                = ""
    opts_vrad@gsnRightString               = ""

  ; Resources to overlay vertical velocity (following M18)
    opts_vvel 	 	 	  	   = True 
    opts_vvel@cnFillOn                     = False

  ; Set contour line colour to match M18 if we choose their colour scheme
    opts_vvel@cnLineColor                 = "grey45"

    opts_vvel@cnLevelSelectionMode         = "ExplicitLevels"
    opts_vvel@cnLevels                     = (/0.5/)
    opts_vvel@cnLineThicknessF             = 5.0

    opts_vvel@cnInfoLabelOn                = False
    opts_vvel@gsnDraw                      = False
    opts_vvel@gsnFrame                     = False
    opts_vvel@cnLineLabelsOn               = False
    opts_vvel@tiMainString                 = ""
    opts_vvel@gsnLeftString                = ""
    opts_vvel@gsnRightString               = ""

;==================================================
; Also add resources for mean quantities (u,v,w)
;==================================================

  ; Reading in colour table using 'cnFillPalette' removes the first 2 colours (black/white)
  ; The indices are thus all 2 values less than they'd normally be

  ; Radial wind 
    opts_vrad2                             = opts_mbgt
    opts_vrad2@cnFillPalette	           = "BlueDarkRed18"
    opts_vrad2@cnLevels			   = (/-10., -5., -3., -2., -1., -0.5, -0.2, \
                                                0, 0.2, 0.5, 1., 2., 3., 5., 10./)
    opts_vrad2@cnFillColors		   = (/0,2,3,4,5,6,7,-1,-1,10,11,12,13,14,15,17/)
    opts_vrad2@lbLabelBarOn                = True
    opts_vrad2@lbBoxEndCapStyle            = "TriangleBothEnds"
    opts_vrad2@lbOrientation		   = "Vertical"

  ; Tangential wind 
    opts_vtan2				   = opts_mbgt
    opts_vtan2@cnFillPalette		   = "radar_new"
    opts_vtan2@cnLevels			   = (/5.0, 10.0, 15.0, 20.0,\
                                               25.0, 30.0, 35.0, 40.0, \
                                               45.0, 50.0, 55.0, 60.0, 65.0/)

    opts_vtan2@cnFillColors		   = (/-1,0,1,2,3,4,5,6,7,8,9,10,11,12/)
    opts_vtan2@lbLabelBarOn                = True
    opts_vtan2@lbBoxEndCapStyle            = "TriangleBothEnds"
    opts_vtan2@lbOrientation               = "Vertical"

  ; Vertical velocity
    opts_vvel2                             = opts_mbgt
    opts_vvel2@cnFillPalette               = "n11"
    opts_vvel2@cnLevels                    = (/0.1, 0.2, 0.3, 0.4, 0.5, 0.6, \
                                               0.7, 0.8, 1.0, 1.2, 1.5, 2.0/)
    opts_vvel2@cnFillColors                = (/-1,0,1,2,3,4,5,6,7,8,9,10,11/)
    opts_vvel2@lbLabelBarOn                = True
    opts_vvel2@lbBoxEndCapStyle            = "TriangleBothEnds"
    opts_vvel2@lbOrientation               = "Vertical"

  ; Absolute angular momentum
    opts_aam				   = opts_vtan2
    delete(opts_aam@cnLevels)
    delete(opts_aam@cnFillColors)
    opts_aam@cnLevels                      = (/1.0, 1.5, 2.0,\
                                               2.5, 3.0, 3.5, 4.0, \
                                               4.5, 5.0, 5.5, 6.0, 6.5/)
    opts_aam@cnFillColors                  = (/-1,0,1,2,3,4,5,6,7,8,9,10,11/)
    opts_aam@lbLabelBarOn                  = True
    opts_aam@lbBoxEndCapStyle              = "TriangleBothEnds"
    opts_aam@lbOrientation                 = "Vertical"
;    opts_aam@lbLabelFontHeightF    	   = 0.01                     ; Label bar font height
;    opts_aam@lbLabelFont                   = "Helvetica"              ; Choose font
;    opts_aam@lbPerimOn                     = False                    ; Perimeter on label bar
;    opts_aam@lbTitleString           	   = \
;                                "Absolute angular momentum (10~S~6~N~ m~S~2~N~ s~S~-1~N~)"
;    opts_aam@lbTitleFontHeightF            = 0.0175                   ; Title font height
;    opts_aam@lbTitleOffsetF                = 0                        ; x-axis title up/down (h/l)
;    opts_aam@pmLabelBarWidthF      	   = 0.6                      ; Label bar width
;    opts_aam@pmLabelBarHeightF     	   = 0.15                     ; Label bar height
;    opts_aam@lbTitlePosition	   	   = "Bottom"		      ; Title position 
;    opts_aam@pmLabelBarOrthogonalPosF 	   = 0.09		      ; Move lb up/down (l/h)

  ; Equivalent potential temperature (theta-e)
    opts_the				   = opts_mbgt

    if (thta .eq. 1) then 
     opts_the@cnFillOn			   = False
     opts_the@cnLinesOn                    = True
     opts_the@cnLevelSelectionMode         = "ExplicitLevels"
     opts_the@cnLineColor		   = "black"
     opts_the@cnLevels			   = ispan(280,400,5)
     opts_the@cnLineThicknessF		   = 4.0
    end if 

  ; Additional plotting resources
    delete([/opts_the@tmYLMode,opts_the@tmYLValues,opts_the@tmYLLabels/])
    delete([/opts_the@trYMinF,opts_the@trYMaxF/])
    opts_the@tmYLMode                      = "Explicit"
    opts_the@tmYLValues                    = fspan(900.0, 100.0, 9)
    opts_the@tmYLLabels                    = (/"2.0","4.0","6.0","8.0","10.0",\
                                               "12.0","14.0","16.0","18.0"/)
    opts_the@tmYLLabelFontHeightF          = 0.02
    opts_the@gsnYAxisIrregular2Linear      = True        ; Linearise the y-axis
    opts_the@trYReverse			   = True	 ; Reverse the y-axis
    opts_the@tfDoNDCOverlay		   = True

    opts_the@trYMinF                       = 125.0       ; y-axis minimum value
    opts_the@trYMaxF                       = 1000.0      ; y-axis maximum value

    if (thta .eq. 0) then 

     opts_the@cnFillOn			   = True 
     opts_the@cnFillPalette                = "GMT_polar"
     opts_the@cnFillColors                 = (/3,4,5,6,7,8,9,-1,11,12,13,14,15,16,17,18/)
     opts_the@cnLevels                     = (/280,285,290,295,300,305,310,315,\
     					       320,325,330,335,340,345,350/)
     opts_the@cnLinesOn		   	   = False
     opts_the@vpWidthF			   = 0.28
     opts_the@vpHeightF			   = 0.28
     opts_the@vpXF			   = 0.20
     opts_the@vpYF			   = 0.20
;     opts_the@tmYLLabelFontHeightF         = 0.0075
     opts_the@lbBoxEndCapStyle             = "TriangleBothEnds"
     opts_the@lbOrientation                = "Vertical"
    
    end if 

;=================
; Plot the data 
;=================

  ; Radial wind (filled contours)
    panel_vrad  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad2) ; [u]

  ; Tangential wind (filled contours)
    panel_vtan  = gsn_csm_contour(wks,mean_vtan(:,it,:),opts_vtan2) ; [v]

  ; Vertical velocity (filled contours)
    panel_vvel  = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel2) ; [w]

  ; Absolute angular momentum (filled contours)
    panel_aam	= gsn_csm_contour(wks,mean_aam(:,it,:),opts_aam)    ; [AAM]
    panel_the	= gsn_csm_contour(wks,mean_the(:,it,:),opts_the)    ; [theta-e]

    if (thta .eq. 1) then 
     overlay(panel_aam, panel_the)
    end if 

  ; Read data into panels before plotting 
    if (thta .eq. 0) then
     panel(0)  = panel_vrad
     panel(1)  = panel_vtan
     panel(2)  = panel_vvel
     panel(3)  = panel_aam	
     panel(4)  = panel_the
    elseif (thta .eq. 1) then
     panel(0)  = panel_aam
    end if 

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    optsP                       = True			   ; Panel plot resources
    optsP@gsnFrame              = False			   ; Do not advance the frame
    optsP@gsnPanelLabelBar      = False       	           ; Panel labelbar

    optsP@pmLabelBarWidthF	= 0.6			   ; Label bar width 
    optsP@pmLabelBarHeightF     = 0.15	     	           ; Label bar height
    optsP@lbLabelFontHeightF	= 0.01	     	           ; Label bar font height 
    optsP@lbLabelFont		= "Helvetica" 	           ; Choose font 
    optsP@lbPerimOn		= False	      	           ; Perimeter on label bar 
    optsP@lbTitleString		= \
    				"Radial momentum budget tendency (m s~S~-1~N~ h~S~-1~N~)"
    optsP@lbTitlePosition	= "Bottom"     		   ; Title position
    optsP@lbTitleFontHeightF	= 0.015			   ; Title font height
    optsP@lbBoxEndCapStyle	= "TriangleBothEnds"	   ; Label bar style
    optsP@lbTitleOffsetF        = 0                        ; x-axis title up/down (h/l)
    optsP@pmLabelBarOrthogonalPosF = -0.002		   ; Move lb up/down (h/l)
    optsP@txString              = title_arr(it)		   ; Plot title string

    if (thta .eq. 0) then 
     optsP@gsnPanelFigureStrings = (/"a) vrad", "b) vtan", "c) vvel", \
     				     "d) AAM", "e) theta-e"/)
    elseif (thta .eq. 1) then 
     optsP@gsnPanelFigureStrings = "" ; (/"a) AAM/theta-e"/)
     optsP@gsnPanelFigureStringsFontHeightF = 0.0175
    end if 

    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.90
    optsP@gsnPanelBottom        = 0.10
    optsP@amJust                = "TopRight"

  ; Add an option to plot with uneven numbers of plots in each row 
    if (thta .eq. 0) then
     optsP@gsnPanelRowSpec = True    
     gsn_panel(wks,panel,(/2,3/),optsP)
     frame(wks)
    elseif (thta .eq. 1) then 
     gsn_panel(wks,panel,(/1,1/),optsP)
     frame(wks)
    end if

   end do     ; End time loop (do it = 0, ts, tf)

end