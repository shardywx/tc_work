; Script to calculate the tangential momentum budget as in Montgomery et al. (2018)

; Run using:

; ncl 'dat="02T12"' 'opt="x11"' 'ens0="em11"' dist=3.0 ts0=54 tf0=63 clr=5
; full=0 rad0=1 mlev0=13 mlev1=34 var0=0 tc_xy_pd.ncl

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'opt'   = output file format ("pdf" or "x11")
; 'ens0'  = ensemble member (em00 -- em11)
; 'dist'  = size of box following storm (degrees)
; 'ts0'   = starting time index (0-119)
; 'tf0'   = ending time index (0-119)
; 'clr'   = colour map for plots ("ncl_default" [1];  "bdr_extra" [2]; "amwg256" [3])
; 'full'  = plot vertical velocity on two levels on same (1) or separate (0) panels 
; 'rad0'  = radius of final plots: 150 km (0) or 200 km (1)
; 'mlev0' = model level for final plot (15, 23, 36, etc)
; 'mlev1' = additional level for final plot (if required) 
; 'var0'  = instantaneous (0) or perturbation (1) fields for final plot

; 'pd' stream --> 'it'
; 'centre'    --> 'it+1'
; 'pc' stream --> 'it+2'

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
  fili_d1     =	systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pd*.nc")
  fili_d      = diri+"/"+fili_d1

; Calculate number of times in each file using built-in function ('nc_times')
  numINFO_d   = nc_times(fili_d)
  numINPUT_d  = numINFO_d[0]		; Number of times in 'pd' stream
  time_d      = numINFO_d[1]		; Array of times from 'pd' stream file 

; Calculate storm motion using built-in function ('st_rm')
; 'dat'     = initialisation time (02T12, 03T00, ...)
; 'ens0'    = ensemble simulation (em00, em01, ...)
; 'diri'    = path to input files (see above)
; 'dist'    = size of box following storm (degrees)
; 'mins'    = analyse 1-h (0) or 5-min (1) data
; 'dir'     = T+00 to T+53 (0), T+54 to T+66 (1), or T+67 to T+83 (2) [ONLY IF 'mins' = 1]

; Call external function 'st_rm'
  storm_rel = st_rm(dat, ens0, "$ar/text/", dist, 0, 0)
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

  utc_date = cd_calendar(time_d, 0)

  year   = tointeger(utc_date(:,0))
  month  = tointeger(utc_date(:,1))
  day    = tointeger(utc_date(:,2))
  hour   = tointeger(utc_date(:,3))
  minute = tointeger(utc_date(:,4))
  second = utc_date(:,5)

; Correct for errors in the code (round up value of hour when minutes = 59)  do ct = 0, numINPUT_d-1
  do ct = 0, numINPUT_d-1
   if (minute(ct).gt.30) then
     hour(ct) = hour(ct)+1
   end if
  end do
  delete(ct)

; Create arrays to hold strings for output plots
  date_str  = new(numINPUT_d,string)
  time_str  = new(numINPUT_d,string)
  time_arr  = new(numINPUT_d,string)
  title_arr = new(numINPUT_d,string)

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
  num_l      = 59      	   	     ; Number of model levels

; Define constants
  rd	= 287.0	  ; Gas constant for dry air (J/kg/K^2)
  cp0	= 1004.0  ; Specific heat of dry air at constant pressure (J/kg/K)
  kp	= 0.286   ; For Exner function calculations (0.286)
  pref	= 1000.0  ; Reference pressure p0 (1000 hPa)

;========================================================================================
;========================================================================================
; Calculating the terms from Montgomery et al. (2018) - their Eq.(4)
; Ignore the perturbation PGF term because it is several orders	of magnitude smaller

; LHS (1) : d{v}/dt	     - Local tendency of the mean tangential wind

; RHS (1) : -{u} * {f+vort}  - Mean radial influx of absolute vertical vorticity
; RHS (2) : -{w} * (d{v}/dz) - Mean vertical advection of mean tangential momentum
; RHS (3) : -{u' * vort'}    - Eddy radial vorticity flux
; RHS (4) : -{w' * (dv'/dz)} - Vertical advection of eddy tangential momentum
; RHS (5) : {D_v}            - Combined diffusive and PBL tendency

;========================================================================================
;========================================================================================

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
  ptb_vort0 = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [vort']
  ptb_dvdz0 = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [dv'/dz]

;===========================================================
; Arrays to hold expanded 'v_int'/'vort_int' arrays (4D)
;===========================================================

  u_int0    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Zonal wind
  v_int0    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Meridional wind
  w_int0    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Vertical velocity
  f_int0    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Coriolis parameter
  vort_int0 = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Relative vorticity
  avo_int0  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Absolute vorticity
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

  vort_int0@description = "Relative vorticity"
  vort_int0@units    = "s~S~-1~N~"
  vort_int0!0           = "lev"
  vort_int0!1           = "rad"
  vort_int0!2           = "azi"
  vort_int0!3           = "time"

  avo_int0@description  = "Relative vorticity"
  avo_int0@units        = "s~S~-1~N~"
  avo_int0!0            = "lev"
  avo_int0!1            = "rad"
  avo_int0!2            = "azi"
  avo_int0!3            = "time"

  f_int0@description = "Coriolis parameter"
  f_int0@units       = "s~S~-1~N~"
  f_int0!0           = "lev"
  f_int0!1           = "rad"
  f_int0!2           = "azi"
  f_int0!3           = "time"

  v_rad0@description  = "Radial velocity"
  v_rad0@units        = "m s~S~-1~N~"
  v_rad0!0            = "lev"
  v_rad0!1            = "rad"
  v_rad0!2            = "azi"
  v_rad0!3            = "time"

  v_tan0@description  = "Tangential velocity"
  v_tan0@units        = "m s~S~-1~N~"
  v_tan0!0            = "lev"
  v_tan0!1            = "rad"
  v_tan0!2            = "azi"
  v_tan0!3            = "time"

;==============================================================
; Arrays for RHS source terms containing eddy variables (4D)
;==============================================================

; Eddy radial vorticity flux: RHS (3)
; {(u' * vort')}
; ptb_vez1  = (ptb_vrad0 * ptb_vort0)
  ptb_vez1  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [R3]

; Vertical advection of eddy tangential momentum: RHS (4)
; [-w' * dv'/dz]
; ptb_vev1  = - (ptb_vvel0 * ptb_dvdz0) 
  ptb_vev1  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [R4]

; Arrays for RHS source terms after azimuthally averaging (3D)
  ptb_vez0  = new((/num_l,dimsizes(radii),numTIMES/),float)		     ; [R3]
  ptb_vev0  = new((/num_l,dimsizes(radii),numTIMES/),float)  		     ; [R4]

;==================================================
; Arrays for horizontal diffusive tendency terms 
;==================================================

; Horizontal diffusivity (K_mh)
  k_mh0	    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)

; Arrays for fractional expressions in diffusive tendency source terms below (4D)
; (1) du/dy, (2) d(v/r)/dr, (3) u/r, (4) r, (5) d/dr, (6) d/dy
  dudy0	    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (1)
  dvdr0     = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (2)
  v_r0      = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (3)
  r_0       = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (4)
  ddr0	    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (5)
  ddy0	    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (6)

; Additional arrays for terms outside main loop
  du        = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ;
  dy        = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ;
  dr        = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ;
  r1        = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ;

; Stress tensor (tau_ry)
; { K_mh * ( (1/r * du/dlm) + d(v/r)/dr ) }
; tau_ry1   = k_mh0 * ( ( (1 / r_0) * dudy0) + ( v_r0 * ddr0) ) 
  tau_ry1   = new((/num_l,dimsizes(radii),numTIMES/),float)		     ; [3D array]

; Total horizontal diffusive tendency of tangential momentum (using 'tau_ry')
; ( (1/(r^2) * {rho}) * d( r ^ 2 * {rho} * {tau_ry} ) /dr)
; v_dh0	    = ( ( 1/( r_0 ^ 2 * mean_rho0) ) * 
; 	      (ddr0 * ( r_0 ^ 2 * mean_rho0 * tau_ry1) ) )
  v_dh0     = new((/num_l,dimsizes(radii),numTIMES/),float)		     ; [3D array]

; Total vertical diffusive tendency of tangential momentum 
; v_dz0	    = (1 /{rho0}) * ( (d/dz) * ({rho0} * {tau_yz}) )  
  v_dz0	    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [4D array]

;========================================
; Define 3D arrays for mean quantities 
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

; Mean absolute vorticity {zeta + f0}
  mean_avo0        = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  mean_avo0!0      = "lev"
  mean_avo0!1      = "rad"
  mean_avo0!2      = "azi"
  mean_avo0!3      = "time"
  mean_avo0@description = "Azimuthally averaged absolute vorticity"
  mean_avo0@units  = "s~S~-1~N~"

; Mean Coriolis force {f0}
  mean_cor0        = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  mean_cor0!0      = "lev"
  mean_cor0!1      = "rad"
  mean_cor0!2      = "azi"
  mean_cor0!3      = "time"
  mean_cor0@description = "Azimuthally averaged Coriolis force"
  mean_cor0@units  = "s~S~-1~N~"

; Vertical derivative of mean tangential wind [ d{v} / dz ]
  mean_dvdz0       = new((/num_l,dimsizes(radii),numTIMES/),float)
  mean_dvdz0!0     = "lev"
  mean_dvdz0!1     = "rad"
  mean_dvdz0!2     = "time"
  mean_dvdz0@description = "Azimuthally averaged 'dv_dz'"
  mean_dvdz0@units = "s~S~-1~N~"

;=========================
; Start multiple loops
;=========================

   print_clock("Working on ensemble member: "+ens0)

   numTIMES    = 120                       ; Total number of times in all files [36]
   ct          = ts0                       ; Counter variable
   llbox       = toint(dist*50)  	   ; Size of lat/lon arrays (see line below)
   dsize       = (/llbox,llbox/)           ; Size of lat/lon arrays (depends on 'dist')
   d0	       = dsize(0)
   d1	       = dsize(1)

   setvalues NhlGetWorkspaceObjectId
     "wsMaximumSize" : 1000000000
   end setvalues

;===============================
; Start loop over input files 
;===============================

   d  = addfile(fili_d,"r")         	 ; Read in 'pd' stream [w,p,t,theta,kmh,tau]

;==================================
; Get the variables we will need
;==================================

   do it = ts0, tf0		; Loop over times in file

    print("Working on time: "+time_str(it))
    itt = it + 1
    title_arr(it) = "Valid at "+time_str(it)+" (T+"+itt+")"

;======================================================
; Read in variables from 'pd' stream (model levels) 
;======================================================

    t0    = lat_arr(it+1,0)
    t1    = lat_arr(it+1,1)
    n0    = lon_arr(it+1,0)
    n1    = lon_arr(it+1,1)
    print_clock("Grid (pd stream) extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

  ; Read in basic variables from 'pd' stream 
    lon1D  = d->longitude_1({n0:n1})                ; longitude [xxx grid points]
    lat1D  = d->latitude_1({t0:t1})                 ; latitude  [xxx grid points]
    lonD   = d->longitude({n0:n1})                  ; longitude [xxx grid points]
    latD   = d->latitude({t0:t1})                   ; latitude  [xxx grid points]
    hybD   = d->hybrid_ht(:)			    ; 63 model ('rho') levels

    hy0    = sprintf("%0.0f",hybD(mlev0))
    hy1    = sprintf("%0.0f",hybD(mlev1))

  ; Read in temperature, pressure, vertical velocity and geopotential height
  ; All variables below -- [hybD | 63] * [latD | 150] * [lon1D | 151]
    t      = d->temp(it-2,:,{t0:t1},{n0:n1})        ; Temperature (K)
    p      = d->p(it-2,:,{t0:t1},{n0:n1})           ; Pressure (Pa)
    w      = d->dz_dt(it-2,:,{t0:t1},{n0:n1})	    ; Vertical velocity (m s-1)
    z      = d->ht(it-2,:,{t0:t1},{n0:n1})    	    ; Geopotential height (m)

  ; Convert units of pressure and geopotential height before continuing
    p  	   = p / 100
    z  	   = z / 10

  ; Add metadata
    p@units = "hPa"
    z@units = "dam"

    dx0     = lonD(1) - lonD(0) ; Grid spacing (longitude)
    dy0     = latD(1) - latD(0) ; Grid spacing (latitude)

;==========================
; Loop over model levels 
;==========================

    do ilev = 0, num_l-1 ; Loop over all pressure levels

     hy = sprintf("%0.0f",hybD(ilev))
     print("Working on model level " + ilev+ " ("+hy0+" m AGL)")

   ; Read in other variables on model levels
     w_plane		    = w(ilev,:,:)	   ; Vertical velocity
     prs_plane		    = p(ilev,:,:)	   ; Pressure 

   ;================================================
   ; Read in pre-calculated storm centre position
   ;================================================

     lat_max    = centre(it+1,0)
     lon_max    = centre(it+1,1)

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
     w_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vertical velocity

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
   ; w_plane     -- latD,  lon1D

     do iang = 0, dimsizes(thetas_rad)-1
      w_int(:,iang)    = linint2_points(lon1D_sub,latD_sub,w_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
     end do

    ;===================================================
    ; Calculate mean (azimuthally-averaged) variables 
    ;===================================================

      do irad = 0, dimsizes(radii)-1

        mean_vvel0(ilev,irad,:,ct)    = avg(w_int(irad,:))      ; [ {w} ]

      end do

    ; Read 'v_int' etc into larger arrays
      w_int0(ilev,:,:,ct)    = (/w_int(:,:)/)

    ; Tidy up
      delete([/w_int,w_plane/])
      delete([/latD_plane,lonD_plane,lon1D_plane,lat1D_plane/])
      delete([/latD_sub,lonD_sub,lat1D_sub,lon1D_sub/])
      delete([/xpos,ypos,xcen,ycen/])

     end do			    ; End model level loop (do ilev = 0, num_l-1)

    ct  = ct + 1		    ; Counter variable (time)

  ; Tidy up before next iteration (all variables from 'pd' stream)
    delete([/lonD,latD,lon1D,lat1D,t,p,w,z/])

   end do     ; End time loop (do it = ts0, tf0)

;====================================
; Calculate perturbation variables
;====================================

      ptb_vvel0(:,:,:,:) = w_int0(:,:,:,:)    - mean_vvel0(:,:,:,:)    ; [w']

;============================================
; Define output file path and plot details 
;============================================

  ; List of possible arrays for plotting 
  ; (1) w_int0 (2) ptb_vvel0

  ; Troubleshooting 
    printMinMax(w_int0(mlev0,:,:,tf0),False)

    do it = ts0, tf0             ; Loop over times in file

    ; Output file location and type 
      if (var0 .eq. 0) then
       output = "$sam/nepartak/images/cyl_coords/vvel_xy_"+dat+"_"+ens0+"_"\
       	         +hy0+"m_"+hy1+"m_inst_"+time_arr(it)
      elseif (var0 .eq. 1) then 
       output = "$sam/nepartak/images/cyl_coords/vvel_xy_"+dat+"_"+ens0+"_"\
                 +hy0+"m_"+hy1+"m_ptb_"+time_arr(it)
      end if 

    ; Set up workstation	     
      wks    = gsn_open_wks(opt,output)

    ; Create graphic for panel plot
      if (full .eq. 1) then 
       panel = new(2,graphic)
      elseif (full .eq. 0) then 
       panel = new(1,graphic)
      end if 

    ; Read in chosen colourmap
      if (var0 .eq. 0) then 
       gsn_define_colormap(wks,"n11")
      elseif (var0 .eq. 1) then 
       gsn_define_colormap(wks,"BlueDarkRed18")
      end if 
 
;==========================
; Options for plotting 
;==========================

    ; Vertical velocity
      opts_vvel                             = True
      opts_vvel@cnFillOn                    = True
      opts_vvel@sfXArray		    = xpos_all(mlev0,:,:,it)
      opts_vvel@sfYArray		    = ypos_all(mlev0,:,:,it)
      opts_vvel@cnLineLabelInterval         = 2.0
      opts_vvel@cnLineLabelFontHeightF      = 0.012
      opts_vvel@cnLineLabelBackgroundColor  = "transparent"
      opts_vvel@cnLineLabelPlacementMode    = "constant"
      opts_vvel@cnLinesOn                   = False
      opts_vvel@cnInfoLabelOn               = False
      opts_vvel@cnLevelSelectionMode        = "ExplicitLevels"

      opts_vvel@cnFillColors                = (/0,2,3,4,5,6,7,8,9,10,11,12/)
      opts_vvel@cnLevels                    = (/0.10, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, \
                                                1.2, 1.5, 2.0, 2.5/)

      opts_vvel@gsnPaperOrientation          = "landscape"
      opts_vvel@tiMainString                 = ""
      opts_vvel@tiMainFontHeightF            = 0.0125
      opts_vvel@gsnLeftString                = ""
      opts_vvel@gsnRightString               = ""
      opts_vvel@gsnMaximize                  = True
      opts_vvel@lbLabelBarOn                 = False
      opts_vvel@gsnDraw                      = False
      opts_vvel@gsnFrame                     = False

      opts_vvel@pmLabelBarWidthF             = 0.6
      opts_vvel@pmLabelBarHeightF            = 0.15
      opts_vvel@lbLabelFontHeightF           = 0.01
      opts_vvel@lbPerimOn                    = False
      opts_vvel@lbLabelFont                  = "Helvetica"
      opts_vvel@lbTitleOn                    = True
      opts_vvel@lbTitleString                = "Vertical velocity (m s~S~-1~N~)"
      opts_vvel@lbTitleFontHeightF           = 0.015
      opts_vvel@lbTitleDirection             = "Across"

    ; Additional plotting resources
      opts_vvel@mpDataBaseVersion            = "Ncarg4_1" ; More recent database
      opts_vvel@mpDataSetName                = "Earth..4" ; High resolution
      opts_vvel@mpOutlineBoundarySets        = "National" ; National borders
      opts_vvel@mpGeophysicalLineColor       = "black"    ; Colour borders black
      opts_vvel@mpGeophysicalLineThicknessF  = 1.0        ; Border line thickness
      opts_vvel@mpGridAndLimbOn              = False      ; Turn on lat/lon lines
      opts_vvel@pmTickMarkDisplayMode        = "Always"   ; Turn on map tickmarks
      opts_vvel@tmXBMajorLengthF             = 0.005      ; Change tickmark length
      opts_vvel@tmXTOn                       = "False"    ; No tickmarks on top x-axis
      opts_vvel@tmYROn                       = "False"    ; No tickmarks on right y-axis
      opts_vvel@gsnMaximize                  = True       ; Maximise plot size
      opts_vvel@gsnAddCyclic                 = False

      radius                                 = 1.0        ; Radius of plot (degrees)

      opts_vvel@mpLimitMode                  = "Corners"
      opts_vvel@mpLeftCornerLatF             = centre(it-1,0)-radius
      opts_vvel@mpLeftCornerLonF             = centre(it-1,1)-radius
      opts_vvel@mpRightCornerLatF            = centre(it-1,0)+radius
      opts_vvel@mpRightCornerLonF            = centre(it-1,1)+radius
      opts_vvel@gsnDraw                      = False      ; Do not draw the plot
      opts_vvel@gsnFrame                     = False      ; Do not advance the frame

    ; Perturbation plot resources 
      opts_ptb0				     = opts_vvel
      delete(opts_ptb0@cnFillColors)
      delete(opts_ptb0@cnLevels)
      opts_ptb0@cnFillColors                 = (/2,4,5,6,7,8,9,0,0,12,13,14,15,16,17,18/)
      opts_ptb0@cnLevels		     = (/-2.0, -1.5, -1.0, -0.5, -0.3, -0.2, -0.1, \
                                                 0, 0.1, 0.2, 0.3, 0.5, 1.0, 1.5, 2.0/)

;=============
; Plot data
;=============

    if (var0 .eq. 0) then				; Instantaneous variables

    ; Vertical velocity
      vvel_plane0      = w_int0(mlev0,:,:,it)
      vvel_plane1      = w_int0(mlev1,:,:,it)

    elseif (var0 .eq. 1) then				; Perturbation variables 

    ; Vertical velocity
      vvel_plane0      = ptb_vvel0(mlev0,:,:,it)
      vvel_plane1      = ptb_vvel0(mlev1,:,:,it)

    end if 

  ;=====================
  ; Lower model level
  ;=====================

  ; Plot filled contours of instantaneous or perturbation quantities
    if (var0 .eq. 0) then 
     plot_vvel0	= gsn_csm_contour_map(wks,vvel_plane0,opts_vvel)	; Vertical velocity (inst)
    elseif (var0 .eq. 1) then
     plot_vvel0 = gsn_csm_contour_map(wks,vvel_plane0,opts_ptb0)	; Vertical velocity (ptb)
    end if 

  ;=====================
  ; Upper model level
  ;=====================

  ; Plot instantaneous or perturbation quantities
    if (var0 .eq. 0) then 
     plot_vvel1 = gsn_csm_contour_map(wks,vvel_plane1,opts_vvel)	; Vertical velocity (inst)
    elseif (var0 .eq. 1) then 
     plot_vvel1 = gsn_csm_contour_map(wks,vvel_plane1,opts_ptb0)	; Vertical velocity (ptb)
    end if 

  ; Plot either one or two panel(s)

    if (full .eq. 0) then                  ; One panel

     panel(0) = plot_vvel0		   ; Panel #1

    elseif (full .eq. 1) then              ; Two panels

     panel(0) = plot_vvel0		   ; Panel #1 
     panel(1) = plot_vvel1		   ; Panel #2 

    end if

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    optsP                       = True                     ; Panel plot resources
    optsP@gsnFrame              = False                    ; Do not advance the frame
    optsP@gsnPanelLabelBar      = True                     ; Panel labelbar
    optsP@pmLabelBarWidthF      = 0.6                      ; Label bar width
    optsP@pmLabelBarHeightF     = 0.15                     ; Label bar height
    optsP@lbLabelFontHeightF    = 0.01                     ; Label bar font height
    optsP@lbLabelFont           = "Helvetica"              ; Choose font
    optsP@lbPerimOn             = False                    ; Perimeter on label bar
    optsP@lbTitleString         = "Vertical velocity (m s~S~-1~N~)"
    optsP@lbTitlePosition       = "Bottom"                 ; Title position
    optsP@lbTitleFontHeightF    = 0.015                    ; Title font height
    optsP@lbBoxEndCapStyle      = "TriangleBothEnds"       ; Label bar style
    optsP@pmLabelBarOrthogonalPosF = -0.002                ; Move lb up/down (h/l)
    optsP@pmLabelBarParallelPosF   = 0.04                  ; Move lb right/left (h/l)
    optsP@lbTitleOffsetF        = 0                        ; x-axis title up/down (h/l)
    optsP@txString              = title_arr(it)            ; Plot title string
    optsP@txPosXF               = 0.525                    ; Title pos (x-direction)

    optsP@gsnPanelFigureStrings = (/"a) Vertical velocity at "+hy0+" m"/)
    optsP@gsnPanelFigureStringsFontHeightF = 0.0125

    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.88
    optsP@gsnPanelBottom        = 0.12
    optsP@amJust                = "TopRight"

    if (full .eq. 1) then 
     gsn_panel(wks,panel,(/1,2/),optsP)		; Two panels
    elseif (full .eq. 0 .or. full .eq. 2) then 
     gsn_panel(wks,panel,(/1,1/),optsP)		; One panel
    end if
    frame(wks)

  ; Tidy up before next iteration
    delete([/vvel_plane0,vvel_plane1/])

   end do     ; End time loop (do it = ts0, tf0))

end