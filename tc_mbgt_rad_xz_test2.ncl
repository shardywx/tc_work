; Script to calculate the momentum budget as in Persing et al. (2013)

; Run using:

; ncl 'opt="x11"' 'ens0="em08"' clev=2 dist=3.0 nbox=1 ts0=30 tf0=36 'sy="lee"' 
; clr=1 gc=1 tc_mbgt_rad_xz_test.ncl

; 'opt'   = output file format ("pdf" or "x11")
; 'ens0'  = ensemble member (em00 -- em11)
; 'clev'  = pressure level for tangential wind panel (0 = 950, 1 = 900, 2 = 850, 3 = 800)
; 'dist'  = size of box following storm (degrees)
; 'ts0'   = starting time index (0–119)
; 'tf0'   = ending time index (0–119)
; 'sy'    = machine (ARCHER [ar] or LEEDS [lee]) for input/output options
; 'clr'	  = colour map for final plots ("BlueDarkRed18" [0] vs "bwr_new" [1])
; 'gc'	  = option for gc_latlon calculations (Instructions [0] vs previous script [1] )

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

; Find 'pb' stream data to read in and plot (pressure level data)
  diri        = "/nfs/a299/TCs/neptark/jul2018"
  fili_prefix = "20160702T1200Z_NPTK_4p4_L80_ra1t_"+ens0
  fili_b1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pb*.nc")
  fili_b      = diri+"/"+fili_b1

; Also find 'pc' and 'pd' streams (model height level data)
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
  numINPUT_d  = numINFO_d[0]		; Number of times in 'pd' stream

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
; Create date/time string arrays for all times in file (use 'pb' stream)
;==========================================================================

  month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
                    "Jul","Aug","Sep","Oct","Nov","Dec"/)

  utc_date = cd_calendar(time_b, 0)

  year   = tointeger(utc_date(:,0))
  month  = tointeger(utc_date(:,1))
  day    = tointeger(utc_date(:,2))
  hour   = tointeger(utc_date(:,3))
  minute = tointeger(utc_date(:,4))
  second = utc_date(:,5)

; Correct for errors in the code (round up value of hour when minutes = 59)
  do ct = 0, numINPUT_b-1
   if (minute(ct).gt.30) then
     hour(ct) = hour(ct)+1
   end if
  end do
  delete(ct)

; Create arrays to hold strings for output plots
  date_str  = new(numINPUT_b,string)
  time_str  = new(numINPUT_b,string)
  time_arr  = new(numINPUT_b,string)
  title_arr = new(numINPUT_b,string)

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

  thetas     = new(73,float)
  pii        = 3.14159265
  radii      = fspan(0.0,1.5,31)     ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,73)   ; Array of azimuth angles (0, 5, 10, etc)
  thetas_rad = (pii/180)*thetas      ; Degrees to radians
  rad_size   = dimsizes(radii)-1     ; Index for later calculation

; Pressure levels (pres)
; (0) 1000, (1) 950, (2) 925, (3) 900, (4) 850, (5) 800,
; (6) 750, (7) 700, (8) 650, (9) 600, (10) 500, (11) 400,
; (12) 300, (13) 250, (14) 200, (15) 150, (16) 100

  num_p = 17      ; Size of pressure level array above (defined on L296)

; Define constants
  rd	= 287.0	  ; Gas constant for dry air (J/kg/K^2)
  cp0	= 1004.0  ; Specific heat of dry air at constant pressure (J/kg/K)
  kp	= 0.286   ; For Exner function calculations (0.286)
  pref	= 1000.0  ; Reference pressure p0 (1000 hPa)

;========================================================================================
;========================================================================================
; Calculating the terms from Montgomery et al. (2018) - their Eq.(8)
; LHS (1) : d{u}/dt	     - Local tendency of the mean radial velocity
; LHS (2) : {u} * d{u}/dr    - Mean radial advection of	radial momentum
; LHS (3) : {(u' * du'/dr) + (v'/r * du'/dy)}  - Mean horiz. adv. of eddy radial momentum

; RHS (1) : -{w} * d{u}/dz   - Mean vertical advection of mean radial momentum
; RHS (2) : -{w' * d{u'}/dz} - Eddy vertical advection of eddy radial momentum
; RHS (3) : {v}^2/r + f{v} - 1/{rho} * d{p}/dr - Mean agradient force 
; RHS (4) : {v'^2/r - 1/rho * dp'/dr}	       - Eddy agradient force 
; RHS (5) : {F_r}	     - Mean radial/vertical diffusive tendency of radial momentum

; All terms except LHS (1) and RHS (5) are per unit mass (?)
;========================================================================================
;========================================================================================

; Arrays for perturbation (eddy) variables (5D)
  ptb_vrad0 = new((/nbox,num_p,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [u']
  ptb_vtan0 = new((/nbox,num_p,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [v']
  ptb_vvel0 = new((/nbox,num_p,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [w']
  ptb_pres0 = new((/nbox,num_p,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [p']
  ptb_rho0  = new((/nbox,num_p,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [ro']

; Arrays for eddy derivative terms (5D)
; (1) du'/dr, (2) du'/dy, (3) du'/dz (4) dp'/dr
  ptb_dudr0 = new((/nbox,num_p,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (1)
  ptb_dudy0 = new((/nbox,num_p,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (2)
  ptb_dudz0 = new((/nbox,num_p,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (3)
  ptb_dpdr0 = new((/nbox,num_p,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (4)

; Arrays for eddy 1/r terms (5D)
; (1) v'/r, (2) v'^2 / r
  ptb_vt_r0  = new((/nbox,num_p,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (1)
  ptb_vttr0  = new((/nbox,num_p,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (2)

; Array for mean 1/r term (4D)
; (1) {v}^2 / r
  mean_vttr0 = new((/nbox,num_p,dimsizes(radii),numTIMES/),float) ; (1) 

;=======================================================================
; Arrays for LHS and RHS source terms containing eddy variables (5D)
;=======================================================================

; Mean horizontal advection of eddy radial momentum: LHS (3)
; {(u' * du'/dr) + (v'/r * du'/dy)}
; ptb_ueh1  = (ptb_vrad0 * ptb_dudr0) + ( (ptb_vtan0 / r) * (ptb_dudy0) )
  ptb_ueh1  = new((/nbox,num_p,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [L1]

; Eddy vertical advection of eddy radial momentum: RHS (2)
; [-w' * du'/dz]
; ptb_uev1  = - (ptb_vvel0 * ptb_dudz0) 
  ptb_uev1  = new((/nbox,num_p,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [R2]

; Eddy agradient force: RHS (3)
; [ {v'^2/r - 1/rho * dp'/dr} ]
; ptb_uegf1 = ( ( (ptb_vtan0) ^ 2) / r) - ( (1/rho) * ptb_dpdr0)
  ptb_uegf1 = new((/nbox,num_p,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [R3]

; Arrays for LHS and RHS source terms after azimuthally averaging (4D)
  ptb_ueh0  = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)		  ; [L1]
  ptb_uev0  = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)  		  ; [R2]
  ptb_uegf0 = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)	          ; [R3]

;========================================
; Define 4D arrays for mean quantities 
;========================================

; Mean tangential wind {v}
  mean_vtan0       = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)
  mean_vtan0!1     = "lev"
  mean_vtan0!2     = "rad"
  mean_vtan0!3     = "time"
  mean_vtan0@description = "Azimuthally averaged tangential wind"
  mean_vtan0@units = "m s~S~-1~N~"

; Mean radial wind {u}
  mean_vrad0       = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)
  mean_vrad0!1     = "lev"
  mean_vrad0!2     = "rad"
  mean_vrad0!3     = "time"
  mean_vrad0@description = "Azimuthally averaged radial wind"
  mean_vrad0@units = "m s~S~-1~N~"

; Mean vertical velocity {w}
  mean_vvel0       = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)
  mean_vvel0!1     = "lev"
  mean_vvel0!2     = "rad"
  mean_vvel0!3     = "time"
  mean_vvel0@description = "Azimuthally averaged vertical velocity"
  mean_vvel0@units = "m s~S~-1~N~"

; Mean pressure {p}
  mean_pres0       = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)
  mean_pres0!1     = "lev"
  mean_pres0!2     = "rad"
  mean_pres0!3     = "time"
  mean_pres0@description = "Azimuthally averaged pressure"
  mean_pres0@units = "hPa"

; Mean density {rho}
  mean_rho0        = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)
  mean_rho0!1      = "lev"
  mean_rho0!2      = "rad"
  mean_rho0!3      = "time"
  mean_rho0@description = "Azimuthally averaged density"
  mean_rho0@units  = "kg m~S~-3~N~"

; Mean Coriolis force {f0}
  mean_cor0        = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)
  mean_cor0!1      = "lev"
  mean_cor0!2      = "rad"
  mean_cor0!3      = "time"
  mean_cor0@description = "Azimuthally averaged Coriolis force"
  mean_cor0@units  = "s~S~-1~N~"

; Mean vertical diffusion (x-component)
  mean_difx0       = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)
  mean_difx0!1     = "lev"
  mean_difx0!2     = "rad"
  mean_difx0!3     = "time"
  mean_difx0@description = "Azimuthally averaged vertical diffusion (x)"
  mean_difx0@units = "m s~S~-1~N~"

; Mean vertical diffusion (y-component)
  mean_dify0       = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)
  mean_dify0!1     = "lev"
  mean_dify0!2     = "rad"
  mean_dify0!3     = "time"
  mean_dify0@description = "Azimuthally averaged vertical diffusion (y)"
  mean_dify0@units = "m s~S~-1~N~"

; Mean vertical diffusion (total)
  mean_diff0       = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)
  mean_diff0!1     = "lev"
  mean_diff0!2     = "rad"
  mean_diff0!3     = "time"
  mean_diff0@description = "Azimuthally averaged vertical diffusion (tot)"
  mean_diff0@units = "m s~S~-1~N~"

; Radial derivative of mean radial wind [ d{u} / dr ]
  mean_dudr0       = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)
  mean_dudr0!1     = "lev"
  mean_dudr0!2     = "rad"
  mean_dudr0!3     = "time"
  mean_dudr0@description = "Azimuthally averaged 'du_dr'"
;  mean_dudr0@units = ""

; Vertical derivative of mean radial wind [ d{u} / dz ]
  mean_dudz0       = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)
  mean_dudz0!1     = "lev"
  mean_dudz0!2     = "rad"
  mean_dudz0!3     = "time"
  mean_dudz0@description = "Azimuthally averaged 'du_dz'"
;  mean_dudz0@units = ""

; Radial derivative of pressure [ d{p} / dr ] 
  mean_dpdr0       = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)
  mean_dpdr0!1     = "lev"
  mean_dpdr0!2     = "rad"
  mean_dpdr0!3     = "time"
  mean_dpdr0@description = "Azimuthally averaged 'dp_dr'"
;  mean_dpdr0@units = ""

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
    plevs	= 19			    ; Number of pressure levels
    llbox	= toint(dist*50)	    ; Size of lat/lon arrays (see line below)
    dsize	= (/llbox,llbox/)           ; Size of lat/lon arrays (depends on 'dist')
    d0		= dsize(0)
    d1		= dsize(1)

    setvalues NhlGetWorkspaceObjectId
      "wsMaximumSize" : 1000000000
    end setvalues

  ; Create array to hold values of frictional tendency
    dtaux_dz = new((/numTIMES,plevs,d0,d1/),float)
    dtauy_dz = new((/numTIMES,plevs,d0,d1/),float)

;===============================
; Start loop over input files 
;===============================

    b  = addfile(fili_b,"r")             ; Read in 'pb' stream [u,v,w,vort,div]
    c  = addfile(fili_c,"r")         	 ; Read in 'pc' stream []
    d  = addfile(fili_d,"r")         	 ; Read in 'pd' stream []

;==================================
; Get the variables we will need
;==================================

   do it = ts0, tf0		; Loop over times in file

    print("Working on time: "+time_str(it))
    itt = it + 1
    title_arr(it) = "Valid at "+time_str(it)+" (T+"+itt+")"

;========================================================
; Read in variables from 'pb' stream (pressure levels)
;========================================================

    t0    = lat_arr(it+3,0)
    t1    = lat_arr(it+3,1)
    n0    = lon_arr(it+3,0)
    n1    = lon_arr(it+3,1)
    print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

    lon1  = b->longitude_1({n0:n1})  ; longitude [xxx grid points]
    lat1  = b->latitude_1({t0:t1})   ; latitude  [xxx grid points]
    lon   = b->longitude({n0:n1})    ; longitude [xxx grid points]
    lat   = b->latitude({t0:t1})     ; latitude  [xxx grid points]

    pres  = b->p         	     ; Pressure (levels)
    plevs = dimsizes(pres)  	     ; Size of pressure level array
    pres@units = "hPa"

  ; Read in horizontal velocity, relative vorticity and geopotential height
    u       = b->u(it,:,{t0:t1},{n0:n1})        ; Zonal velocity (m/s)
    v       = b->v(it,:,{t0:t1},{n0:n1})        ; Meridional velocity (m/s)
    w       = b->dz_dt(it,:,{t0:t1},{n0:n1})    ; Vertical velocity (m/s)
    omega   = b->omega(it,:,{t0:t1},{n0:n1})	; Omega (Pa/s)
    vort    = b->rvor(it,:,{t0:t1},{n0:n1})     ; Relative vorticity (/s)
    pv	    = b->pv(it,:,{t0:t1},{n0:n1})	; Potential vorticity ()
    z       = b->ht(it,:,{t0:t1},{n0:n1})       ; Geopotential height (m)

  ; Convert to units we want
    z       = z/10
    z@units = "dam"			        ; Convert to dam
    vort    = vort * (10 ^ -6) 	 	        ; Convert to s-1
    pv	    = pv * (10 ^ 6)			; Convert to PVU [10 ^ 6 K m2 s-1 kg-1]

    dx      = lon(1) - lon(0) ; Grid spacing (longitude)
    dy      = lat(1) - lat(0) ; Grid spacing (latitude)

  ; Pressure levels (pres)
  ; (0) 1000, (1) 950, (2) 925, (3) 900, (4) 850, (5) 800,
  ; (6) 750, (7) 700, (8) 650, (9) 600, (10) 500, (11) 400,
  ; (12) 300, (13) 250, (14) 200, (15) 150, (16) 100

  ; Coriolis parameter and absolute vorticity
    f1   = coriolis_param(lat)
    f0   = conform_dims(dimsizes(vort),f1,1)
    copy_VarMeta(vort,f0)
    copy_VarCoords(vort,f0)

    avo  = vort + f0
    copy_VarMeta(vort,avo)
    copy_VarAtts(vort,avo)
    copy_VarCoords(vort,avo)
    avo@description = "Absolute vorticity"
    avo@long_name   = "Absolute vorticity"

;==================================================================
; Read in variables from 'pc' stream (model levels)
; B/c of differences in file size, access variables using 'it+4'
;==================================================================

  ; Read in basic variables from 'pc' stream
    lon1C  = c->longitude_1({n0:n1})              ; longitude [xxx grid points]
    lat1C  = c->latitude_1({t0:t1})               ; latitude  [xxx grid points]
    lonC   = c->longitude({n0:n1})                ; longitude [xxx grid points]
    latC   = c->latitude({t0:t1})                 ; latitude  [xxx grid points]
    hybC   = c->hybrid_ht(:) 		          ; 63 hybrid levels (model 'theta' levs)

    u_varname = "x-wind"
    v_varname = "y-wind"

  ; Read in horizontal wind components 
  ; Both variables below -- [hybC | 63] * [latC | 150] * [lonC | 150]
    u_hyb  = c->$u_varname$(it+4,:,{t0:t1},{n0:n1}) ; Zonal wind	
    v_hyb  = c->$v_varname$(it+4,:,{t0:t1},{n0:n1}) ; Meridional wind 

    dxC    = lonC(1) - lonC(0) ; Grid spacing (longitude)
    dyC    = latC(1) - latC(0) ; Grid spacing (latitude)

;==================================================================
; Read in variables from 'pd' stream (model levels) 
; B/c of differences in file size, access variables using 'it+2'
;==================================================================

  ; Read in basic variables from 'pd' stream 
    lon1D  = d->longitude_1({n0:n1})                ; longitude [xxx grid points]
    lat1D  = d->latitude_1({t0:t1})                 ; latitude  [xxx grid points]
    lonD   = d->longitude({n0:n1})                  ; longitude [xxx grid points]
    latD   = d->latitude({t0:t1})                   ; latitude  [xxx grid points]
    hybD   = d->hybrid_ht(:)			    ; 63 hybrid levels (model 'rho' levs)

  ; Read in diffusivity and wind stress
  ; 'str_xz' -- [hybD | 63] * [latD | 150]  * [lonD | 151]
  ; 'str_yz' -- [hybD | 63] * [lat1D | 150] * [lon1D | 151]
  ; 'kmh'    -- [hybD | 63] * [latD | 150]  * [lon1D | 151]
    str_xz = d->taux(it+2,:,{t0:t1},{n0:n1})        ; Wind stress, x-comp (N m/2)
    str_yz = d->tauy(it+2,:,{t0:t1},{n0:n1})        ; Wind stress, y-comp (N m/2)
    kmh	   = d->unspecified(it+2,:,{t0:t1},{n0:n1}) ; Horizontal diffusivity (m^2 s-1)

  ; Read in temperature, pressure, theta, vertical velocity and geopotential height
  ; All variables below -- [hybD | 63] * [latD | 150] * [lon1D | 151]
    t_hyb  = d->temp(it+2,:,{t0:t1},{n0:n1})        ; Temperature (K)
    p_hyb  = d->p(it+2,:,{t0:t1},{n0:n1})           ; Pressure (Pa)
    th_hyb = d->theta(it+2,:,{t0:t1},{n0:n1})	    ; Potential temperature (K)
    w_hyb  = d->dz_dt(it+2,:,{t0:t1},{n0:n1})	    ; Vertical velocity (m s-1)
    z_hyb  = d->ht(it+2,:,{t0:t1},{n0:n1})	    ; Geopotential height (m)

  ; Calculate dry air density from pressure and temperature
    rho	   = p_hyb / (rd * t_hyb)  	    	    ; Calculate dry air density
    copy_VarCoords(p_hyb, rho)
    rho@description = "Dry air density"		    ; Add metadata
    rho@units = "kg m~S~-3~N~"		            ; Add metadata

  ; Convert units of pressure and geopotential height before continuing
    p_hyb  = p_hyb / 100
    z_hyb  = z_hyb / 10

  ; Add metadata
    p_hyb@units = "hPa"
    z_hyb@units = "dam"

  ;;;;;;;;;;; EDIT TO CALCULATE VORT ON MODEL LEVELS ;;;;;;;;;;;
  ; Calculate relative vorticity from horizontal wind components 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;========================================
; Calculations of diffusive tendencies
;========================================

  ; Array of hybrid levels to approximately match pressure levels every 50 hPa
  ;;; EDIT TO MATCH NEW PRESSURE LEVELS (fewer and not always 50 hPa apart) ;;;

  ; Array [/1000, 950, 900, 850, 800, 750, etc, ..., 100/]
    lev_arr = (/1, 8, 13, 16, 19, 21, 23, 25, 27, 29, \
                31, 33, 35, 37, 39, 42, 47, 51, 56/)
    lev_arr@units = "m"

    l0    = 0  ; Subset of vertical levels (19 levels)
    l1    = 18 ; Level 56 of original array ~ 100 hPa (Int. Standard Atmosphere)

  ; Arrays to hold subset of 'str_xz' and 'str_yz', on levels in 'lev_arr'
  ; The dimensions of 'dsize' latitude (0) and longitude (1)
    fric_x_hgt  = new((/dimsizes(lev_arr),d0,d1/),float)
    fric_y_hgt  = new((/dimsizes(lev_arr),d0,d1/),float)

  ; Array of pressure levels with the same dimensions as 'lev_arr'
  ; Pressure levels (pres)
  ; (0) 1000, (1) 950, (2) 925, (3) 900, (4) 850, (5) 800,
  ; (6) 750, (7) 700, (8) 650, (9) 600, (10) 500, (11) 400,
  ; (12) 300, (13) 250, (14) 200, (15) 150, (16) 100

    pres       = (/1000,950,925,900,850,800,750,700,\
                   650,600,500,400,300,250,200,150,100/)
    pres@units = "hPa"
    plevs      = dimsizes(pres)

    cp         = 0    ; Counter variable

  ; Create array of pressure values with same dimensions as 'fric_x_hgt'
    pres_fr    = conform_dims(dimsizes(fric_x_hgt), pres, 0)
    pres_fr = pres_fr * 100
    pres_fr@units = "Pa"

  ;;;;;;;;;;;;;;; EDIT FROM HERE ;;;;;;;;;;;;;;;
  ;;;; CAN I SWITCH TO PURELY MODEL LEVELS? ;;;;

  ; TC latitude  -- centre(it,0)
  ; TC longitude -- centre(it,1)

;===============================
; Loop over pressure levels 
;===============================

    do ilev = 0, num_p-1 ; Loop over all pressure levels

   ;========================================
   ; Model levels ['pc' and 'pf' streams]
   ;========================================

     m0 = lev_arr(ilev)       ; Selected model levels, chosen to match pressure levels
     hy = sprintf("%0.0f",hyb(m0))
     print("Working on hybrid level " + m0+ " ("+hy+" m AGL)")

   ; Calculate components of friction term
     ip1 = min((/m0+1,lev_arr(l1)/)) ; Account for top and bottom levels
     im1 = max((/m0-1,lev_arr(l0)/))

   ; Current level [wind stress]
     tau_xz = str_xz(m0,0:d0-1,0:d1-1)
     tau_yz = str_yz(m0,0:d0-1,0:d1-1)

   ; Density ['pf' stream]
     rho    = rho1(m0,0:d0-1,0:d1-1)
     rho_a  = (1/rho) ; Specific volume (inverse density)

   ; Levels above/below for vertical derivatives [wind stress]
     tau_x1 = str_xz(ip1,0:d0-1,0:d1-1)
     tau_x0 = str_xz(im1,0:d0-1,0:d1-1)
     tau_y1 = str_yz(ip1,0:d0-1,0:d1-1)
     tau_y0 = str_yz(im1,0:d0-1,0:d1-1)

   ; Vertical derivative of 'str_xz'
     dtaux_dz(ct,ilev,:,:) = (tau_x1 - tau_x0) / (hyb(ip1) - hyb(im1))

   ; Vertical derivative of 'str_yz'
     dtauy_dz(ct,ilev,:,:) = (tau_y1 - tau_y0) / (hyb(ip1) - hyb(im1))

   ; X-component of vertical diffusion
     difv_xplane            = rho_a(:,:) * dtaux_dz(ct,ilev,:,:)
     copy_VarCoords(tau_xz,difv_xplane)
     copy_VarMeta(tau_xz,difv_xplane)
     difv_xplane@units       = "m s~S~-2~N~"
     difv_xplane@title       = "x-component of vertical diffusion"
     difv_xplane@long_name   = "x-component of vertical diffusion"
     difv_xplane@name        = "difv_x"

   ; Y-component of vertical diffusion
     difv_yplane             = rho_a(:,:) * dtauy_dz(ct,ilev,:,:)
     copy_VarCoords(tau_yz,difv_yplane)
     copy_VarMeta(tau_yz,difv_yplane)
     difv_yplane@units       = "m s~S~-2~N~"
     difv_yplane@title       = "y-component of vertical diffusion"
     difv_yplane@long_name   = "y-component of vertical diffusion"
     difv_yplane@name        = "difv_y"

   ; Agradient force (inverse density and pressure)
     rho_plane              = rho(:,:)
     prs_plane              = p_ht(m0,:,:)

   ;=================================
   ; Pressure levels ['pb' stream]
   ;=================================

     print("Working on: " + pres(ilev) + " hPa")
     p = pres(ilev)

     u_plane       = u(ilev,:,:)
     v_plane       = v(ilev,:,:)
     w_plane       = w(ilev,:,:)
     f_plane	   = f0(ilev,:,:)
     vort_plane	   = vort(ilev,:,:)
     avo_plane	   = avo(ilev,:,:)
     geo_plane     = z(ilev,:,:)

     vort_max      = max(vort_plane)
     vort_smth     = smth9_Wrap(vort_plane, 0.5, 0.5, True)
     vort_max_smth = max(vort_smth)

     zmin          = toint(min(geo_plane)) ; Minimum geopotential on pressure level
     zmax          = toint(max(geo_plane)) ; Maximum geopotential on pressure level

     dims          = dimsizes(vort_plane)
     vort1d        = ndtooned(vort_plane)
     inds          = ind_resolve(maxind(vort1d),dims)
     vort1d_sm     = ndtooned(vort_smth)
     inds_sm       = ind_resolve(maxind(vort1d_sm),dims)

   ; Find latitude and longitude of storm at each time (and pressure level)
     lat_max 	   = t0 + (dy * inds(0,0))
     lon_max 	   = n0 + (dx * inds(0,1))

   ; Difference between 'lat_max'/'lon_max' and values read in using NCL function
     dlat 	   = abs(lat_max - centre(it,0))
     dlon	   = abs(lon_max - centre(it,1))

   ; If our calculations of latitude and longitude deviate too much from the values
   ; we read in using the NCL function ('st_rm'), then calculate a new centre
   ; using the geopotential height. Constrain the new calculation using the values
   ; from 'st_rm', to prevent spurious centres being identified and analysed. 

     if (dlat .gt. 0.75 .or. dlon .gt. 0.75) then 

      lt1 = centre(it,0) - 0.5
      lt2 = centre(it,0) + 0.5
      ln1 = centre(it,1) - 0.5
      ln2 = centre(it,1) + 0.5

      geo_plane0 = geo_plane({lt1:lt2},{ln1:ln2})       ; Smaller grid (1.0 x 1.0 deg)
      geo_min0   = min(geo_plane0)                      ; Find minimum on smaller grid
      dims_h0    = dimsizes(geo_plane0)                 ; Size of smaller grid
      geo_1d0    = ndtooned(geo_plane0)                 ; Create 1-D array
      inds_h0    = ind_resolve(minind(geo_1d0),dims_h0) ; Find index of minimum

      lat_max 	 = lt1 + (dy * inds_h0(0,0))
      lon_max 	 = ln1 + (dx * inds_h0(0,1))
      print("New centre: "+lat_max+" degrees N, "+lon_max+" degrees E")

     end if 

   ;==============================================================
   ; Create arrays for later switch to cylindrical coordinates
   ;==============================================================

   ; u,v grid
     lat_plane  = ind(t0.le.lat.and.lat.le.t1) ; Indices of latitude array
     lon_plane  = ind(n0.le.lon.and.lon.le.n1) ; Indices of longitude array

     lat_size = dimsizes(lat_plane)		  	  	; Array sizes
     lon_size = dimsizes(lon_plane)

     lat_sub = lat(lat_plane(0):lat_plane(lat_size-1)) ; Create subset of points 
     lon_sub = lon(lon_plane(0):lon_plane(lon_size-1)) ; between values given above

   ; w grid
     lat1_plane  = ind(t0.le.lat1.and.lat1.le.t1) ; Indices of latitude array
     lon1_plane  = ind(n0.le.lon1.and.lon1.le.n1) ; Indices of longitude array

     lat1_size = dimsizes(lat1_plane)                             ; Array sizes
     lon1_size = dimsizes(lon1_plane)

     lat1_sub = lat1(lat1_plane(0):lat1_plane(lat1_size-1)) ; Create subset of points
     lon1_sub = lon1(lon1_plane(0):lon1_plane(lon1_size-1)) ; between values given above

     delete([/vort1d,inds,vort1d_sm,vort_smth,inds_sm/])

;=====================================================================================
; Create array of box positions either side of the original [Davis & Galarneau 2009]
;=====================================================================================

     if (nbox.eq.49) then
      lat_box = (/lat_max-(3*dy),lat_max-(2*dy),lat_max-(dy),lat_max,\
                  lat_max+(dy),lat_max+(2*dy),lat_max+(3*dy)/)
      lon_box = (/lon_max-(3*dx),lon_max-(2*dx),lon_max-(dx),lon_max,\
                  lon_max+(dx),lon_max+(2*dx),lon_max+(3*dx)/)
     else if (nbox.eq.25) then
      lat_box = (/lat_max-(2*dy),lat_max-(dy),lat_max,\
                  lat_max+(dy),lat_max+(2*dy)/)
      lon_box = (/lon_max-(2*dx),lon_max-(dx),lon_max,\
                  lon_max+(dx),lon_max+(2*dx)/)
     else if (nbox.eq.9) then
      lat_box = (/lat_max-(dy),lat_max,lat_max+(dy)/)
      lon_box = (/lon_max-(dx),lon_max,lon_max+(dx)/)
     else if (nbox.eq.1) then
      lat_box = (/lat_max/)
      lon_box = (/lon_max/)
     end if
     end if
     end if
     end if

     lsize   = dimsizes(lat_box) * dimsizes(lon_box)

;===========================================
; Switch to cylindrical coordinate system
;===========================================

; Create new arrays describing cylindrical coordinate system (2D for now)
  xpos = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float)
  ypos = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float)

; New arrays containing regular variables we want to plot
  u_int    = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float) ; Zonal velocity
  v_int    = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float) ; Meridional velocity
  w_int    = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float) ; Vertical velocity
  f_int	   = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float) ; Coriolis force 
  prs_int  = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float) ; Pressure 
  rho_int  = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float) ; Density
  v_rad    = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float) ; Radial velocity
  v_tan    = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float) ; Tangential velocity

; Additional arrays for horizontal and vertical diffusion
  difv_x   = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float) ; Vert. diff (x)
  difv_y   = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float) ; Vert. diff (y)
  difv_rad = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float) ; Vert. diff (rad)
  difv_tan = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float) ; Vert. diff (tan)
  difv_tot = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float) ; Vert. diff (total)

  difh_x   = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float) ; Horiz. diff (x)
  difh_y   = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float) ; Horiz. diff	(y)
  difh_rad = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float) ; Horiz. diff (rad)
  difh_tan = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float) ; Horiz. diff (tan)

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

;=======================
; Vertical diffusion
;=======================

  copy_VarAtts(difv_xplane,difv_x)
  difv_x@description   = "Vertical diffusion (x-component)"
  difv_x@units         = "m s~S~-2~N~"

  copy_VarAtts(difv_yplane,difv_y)
  difv_y@description   = "Vertical diffusion (y-component)"
  difv_y@units         = "m s~S~-2~N~"

  copy_VarAtts(difv_xplane,difv_rad)
  difv_rad@description = "Vertical diffusion (radial component)"
  difv_rad@units       = "m s~S~-2~N~"

  copy_VarAtts(difv_yplane,difv_tan)
  difv_tan@description = "Vertical diffusion (tangential component)"
  difv_tan@units       = "m s~S~-2~N~"

  copy_VarAtts(difv_yplane,difv_tot)
  difv_tot@description = "Vertical diffusion (total)"
  diff_int@units       = "m s~S~-2~N~"

;=========================
; Horizontal diffusion
;=========================

  copy_VarAtts(difv_xplane,difh_x)
  difv_x@description = "Horizontal diffusion (x-component)"
  difv_x@units       = "m s~S~-2~N~"

  copy_VarAtts(difv_yplane,difh_y)
  difv_y@description = "Horizontal diffusion (y-component)"
  difv_y@units       = "m s~S~-2~N~"

  copy_VarAtts(difv_xplane,difh_rad)
  difv_rad@description = "Horizontal diffusion (radial component)"
  difv_rad@units       = "m s~S~-2~N~"

  copy_VarAtts(difv_yplane,difh_tan)
  difv_tan@description = "Horizontal diffusion (tangential component)"
  difv_tan@units       = "m s~S~-2~N~"

  v_rad@description    = "Radial velocity"
  v_rad@units          = "m s~S~-1~N~"

  v_tan@description    = "Tangential velocity"
  v_tan@units          = "m s~S~-1~N~"

;========================================
; Loop over array of lat/lon locations
;========================================

   do i = 0, dimsizes(lon_box)-1	; Longitude
    do j = 0, dimsizes(lat_box)-1	; Latitude     

    a = (dimsizes(lon_box) * i) + j

  ; Create cylindrincal coordinate array from original latitude/longitude grid
    do irad = 0, dimsizes(radii)-1
      r = radii(irad)
      do iang = 0, dimsizes(thetas)-1
        theta = thetas_rad(iang)

          xpos(a,irad,iang) = (/lon_box(i)/) + r*cos(theta)
	  ypos(a,irad,iang) = (/lat_box(j)/) + r*sin(theta)

      end do
    end do

  ; Fill new arrays with interpolated values on the cylindrical grid 
  ; Check to make sure that diffusion arrays are matched correctly w/ 'lat/lon_sub'

    do iang = 0, dimsizes(thetas_rad)-1
      u_int(a,:,iang)    = linint2_points(lon_sub,lat_sub,u_plane,False,\
                                          xpos(a,:,iang),ypos(a,:,iang),0)
      v_int(a,:,iang)    = linint2_points(lon_sub,lat_sub,v_plane,False,\
                                          xpos(a,:,iang),ypos(a,:,iang),0)
      w_int(a,:,iang)    = linint2_points(lon1_sub,lat1_sub,w_plane,False,\
                                          xpos(a,:,iang),ypos(a,:,iang),0)
      f_int(a,:,iang)	 = linint2_points(lon_sub,lat_sub,f_plane,False,\
				          xpos(a,:,iang),ypos(a,:,iang),0)
      difv_x(a,:,iang)   = linint2_points(lon_sub,lat_sub,difv_xplane,False,\
					  xpos(a,:,iang),ypos(a,:,iang),0)
      difv_y(a,:,iang)   = linint2_points(lon_sub,lat_sub,difv_yplane,False,\
					  xpos(a,:,iang),ypos(a,:,iang),0)
      prs_int(a,:,iang)  = linint2_points(lon_sub,lat_sub,prs_plane,False,\
                                          xpos(a,:,iang),ypos(a,:,iang),0)
      rho_int(a,:,iang)  = linint2_points(lon_sub,lat_sub,rho_plane,False,\
					  xpos(a,:,iang),ypos(a,:,iang),0)
    end do

  ; Calculate radial/tangential windspeed and friction at each point (using unit vectors)
    do iang = 0, dimsizes(thetas)-1
      theta = thetas_rad(iang)

     ; Tangential and radial wind 
       v_rad(a,:,iang) = u_int(a,:,iang)*cos(theta) + v_int(a,:,iang)*sin(theta)
       v_tan(a,:,iang) = -u_int(a,:,iang)*sin(theta) + v_int(a,:,iang)*cos(theta)

     ; Total vertical diffusive	tendency of radial momentum
       difv_tot(a,:,iang) = difv_y(a,:,iang)*sin(theta) + difv_x(a,:,iang)*cos(theta)

;     ; Tangential and radial components of vertical diffusion
;       difv_rad(a,:,iang) = difv_x(a,:,iang)*cos(theta) + difv_y(a,:,iang)*sin(theta)
;       difv_tan(a,:,iang) = -difv_x(a,:,iang)*sin(theta) + difv_y(a,:,iang)*cos(theta)

;     ; Total vertical diffusive tendency of tangential momentum 
;       diff_int(a,:,iang) = difv_y(a,:,iang)*cos(theta) - difv_x(a,:,iang)*sin(theta)

    end do

     xpos!1       = "rad"
     xpos!2       = "azi"
     xpos&rad     = radii
     xpos&azi     = thetas

     ypos!1       = "rad"
     ypos!2       = "azi"
     ypos&rad     = radii
     ypos&azi     = thetas

     u_int!1      = "rad"
     u_int!2      = "azi"
     u_int&rad    = radii
     u_int&azi    = thetas

     v_int!1      = "rad"
     v_int!2      = "azi"
     v_int&rad    = radii
     v_int&azi    = thetas

     f_int!1      = "rad"
     f_int!2      = "azi"
     f_int&rad    = radii
     f_int&azi    = thetas

     w_int!1      = "rad"
     w_int!2      = "azi"
     w_int&rad    = radii
     w_int&azi    = thetas

     prs_int!1    = "rad"
     prs_int!2    = "azi"
     prs_int&rad  = radii
     prs_int&azi  = thetas

     rho_int!1    = "rad"
     rho_int!2    = "azi"
     rho_int&rad  = radii
     rho_int&azi  = thetas

     v_tan!1      = "rad"
     v_tan!2      = "azi"
     v_tan&rad    = radii
     v_tan&azi    = thetas

     v_rad!1      = "rad"
     v_rad!2      = "azi"
     v_rad&rad    = radii
     v_rad&azi    = thetas

     difv_x!1     = "rad"
     difv_x!2     = "azi"
     difv_x&rad   = radii
     difv_x&azi   = thetas

     difv_y!1     = "rad"
     difv_y!2     = "azi"
     difv_y&rad   = radii
     difv_y&azi   = thetas

     difv_rad!1   = "rad"
     difv_rad!2   = "azi"
     difv_rad&rad = radii
     difv_rad&azi = thetas

     difv_tan!1   = "rad"
     difv_tan!2   = "azi"
     difv_tan&rad = radii
     difv_tan&azi = thetas

     difv_tot!1   = "rad"
     difv_tot!2   = "azi"
     difv_tot&rad = radii
     difv_tot&azi = thetas

     difh_x!1     = "rad"
     difh_x!2     = "azi"
     difh_x&rad   = radii
     difh_x&azi   = thetas

     difh_y!1     = "rad"
     difh_y!2     = "azi"
     difh_y&rad   = radii
     difh_y&azi   = thetas

     difh_rad!1   = "rad"
     difh_rad!2   = "azi"
     difh_rad&rad = radii
     difh_rad&azi = thetas

     difh_tan!1   = "rad"
     difh_tan!2   = "azi"
     difh_tan&rad = radii
     difh_tan&azi = thetas

    ;===================================================
    ; Calculate mean (azimuthally-averaged) variables 
    ;===================================================

      do irad = 0, dimsizes(radii)-1

        mean_vtan0(a,ilev,irad,ct) = avg(v_tan(a,irad,:))      ; [ {v} ]
	mean_vrad0(a,ilev,irad,ct) = avg(v_rad(a,irad,:))      ; [ {u} ]
	mean_vvel0(a,ilev,irad,ct) = avg(w_int(a,irad,:))      ; [ {w} ]
	mean_pres0(a,ilev,irad,ct) = avg(prs_int(a,irad,:))    ; [ {p} ]
	mean_diff0(a,ilev,irad,ct) = avg(diff_int(a,irad,:))   ; [ {fric_xz + fric_yz} ]
	mean_rho0(a,ilev,irad,ct)  = avg(rho_int(a,irad,:))    ; [ {rho} ]
	mean_cor0(a,ilev,irad,ct)  = avg(f_int(a,irad,:))      ; [ {f} ]

      end do

      fcor = mean_cor0 * mean_vtan0

    ;====================================
    ; Calculate perturbation variables 
    ;====================================

      do irad = 0, dimsizes(radii)-1

        ptb_vtan0(a,ilev,irad,:,ct) = v_tan(a,irad,:) - mean_vtan0(a,ilev,irad,ct) ; [v']
	ptb_vvel0(a,ilev,irad,:,ct) = w_int(a,irad,:) - mean_vvel0(a,ilev,irad,ct) ; [w']
	ptb_vrad0(a,ilev,irad,:,ct) = v_rad(a,irad,:) - mean_vrad0(a,ilev,irad,ct) ; [u']
	ptb_pres0(a,ilev,irad,:,ct) = prs_int(a,irad,:)	- mean_pres0(a,ilev,irad,ct);[p']
	ptb_rho0(a,ilev,irad,:,ct)  = rho_int(a,irad,:)				  ; [rho]

      end do

    ;=================================================================
    ; Calculate eddy terms involving azimuthal derivatives [du'/dy]
    ;=================================================================

       do irad = 1, dimsizes(radii)-1

        do th0 = 0, dimsizes(thetas)-1
       
         tp1 = min((/th0+1,dimsizes(thetas)-1/))
       	 tm1 = max((/th0-1,0/))

       ; Calculate the azimuthal distance ('dy') between all points in array
       ; Arguments 1 and 2 represent the lat/lon of the first point, respectively
       ; Arguments 3 and 4 represent the lat/lon of the	second point, respectively
       ; Argument 5 ('2') tells the script to avoid extra interpolated points
       ; Argument 6 ('3') returns the distance in metres

         if (gc.eq.0) then 
  	  dl  = gc_latlon(xpos(a,irad,tm1),ypos(a,irad,tm1),\
	                  xpos(a,irad,tp1),ypos(a,irad,tp1),2,3) 
         else if (gc.eq.1) then 
          dl  = gc_latlon(ypos(a,irad,tm1),xpos(a,irad,tm1),\
                          ypos(a,irad,tp1),xpos(a,irad,tp1),2,3)
         end if 
	 end if 

       ; Calculate [du'/dy]
         ptb_dudy0(a,ilev,irad,th0,ct) = (ptb_vrad0(a,ilev,irad,tp1,ct) - \
 				          ptb_vrad0(a,ilev,irad,tm1,ct)) / dl

        end do			      ; End loop over azimuth angles
       
       end do			      ; End loop over radii 

    ;======================================================================
    ; Calculate eddy terms involving radial derivatives [du'/dr, dp'/dr]
    ;======================================================================

       ptb_pres0 = ptb_pres0 * 100

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

         if (gc.eq.0) then

        ; Radial distance between each point and its nearest neighbour
          dr  = gc_latlon(xpos(a,rm1,th0),ypos(a,rm1,th0),\
                          xpos(a,rp1,th0),ypos(a,rp1,th0),2,3)

        ; Radial distance between each point and the origin
          r1  = gc_latlon(xpos(a,0,th0),ypos(a,0,th0),\
                          xpos(a,irad,th0),ypos(a,irad,th0),2,3)

         else if (gc.eq.1) then 

        ; Radial distance between each point and its nearest neighbour
          dr  = gc_latlon(ypos(a,rm1,th0),xpos(a,rm1,th0),\
                          ypos(a,rp1,th0),xpos(a,rp1,th0),2,3)

        ; Radial distance between each point and the origin
          r1  = gc_latlon(ypos(a,0,th0),xpos(a,0,th0),\
                          ypos(a,irad,th0),xpos(a,irad,th0),2,3)

	 end if 
	 end if 		  

;	 print("Radius = "+radii(irad)+" degrees ; azimuth angle = "\
;	       +thetas(th0)+" degrees")
;	 print("Radial distance (r1) = "+r1+" metres")

       ; Calculate increments 
	 du   = ptb_vrad0(a,ilev,rp1,th0,ct) - ptb_vrad0(a,ilev,rm1,th0,ct)
	 dp   = ptb_pres0(a,ilev,rp1,th0,ct) - ptb_pres0(a,ilev,rm1,th0,ct)
	 vp   = ptb_vtan0(a,ilev,irad,th0,ct)
	 vpp  = (ptb_vtan0(a,ilev,irad,th0,ct) ^ 2)

       ; Calculate [du'/dr]
         ptb_dudr0(a,ilev,irad,th0,ct)  = du / dr

       ; Calculate [dp'/dr]
         ptb_dpdr0(a,ilev,irad,th0,ct)  = dp / dr

       ; Calculate [v'/r]
         ptb_vt_r0(a,ilev,irad,th0,ct)  = vp / r1 

       ; Calculate [v'^2/r]
       	 ptb_vttr0(a,ilev,irad,th0,ct)  = vpp / r1

        end do                        ; End loop over azimuth angles

       end do                         ; End loop over radii

       ptb_pres0 = ptb_pres0 / 100

    ;=========================================================================
    ; Calculate mean terms involving radial derivatives [ d{u}/dr, d{p}/dr ]
    ;=========================================================================
      
       mean_pres0 = mean_pres0 * 100

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

         if (gc.eq.0) then

        ; Radial distance between each point and its nearest neighbour
          dr  = gc_latlon(xpos(a,rm1,th0),ypos(a,rm1,th0),\
                          xpos(a,rp1,th0),ypos(a,rp1,th0),2,3)

        ; Radial distance between each point and the origin
          r1  = gc_latlon(xpos(a,0,th0),ypos(a,0,th0),\
                          xpos(a,irad,th0),ypos(a,irad,th0),2,3)

         else if (gc.eq.1) then

        ; Radial distance between each point and its nearest neighbour
          dr  = gc_latlon(ypos(a,rm1,th0),xpos(a,rm1,th0),\
                          ypos(a,rp1,th0),xpos(a,rp1,th0),2,3)

        ; Radial distance between each point and the origin
          r1  = gc_latlon(ypos(a,0,th0),xpos(a,0,th0),\
                          ypos(a,irad,th0),xpos(a,irad,th0),2,3)

         end if
         end if

       ; Calculate increments
         dmu  = mean_vrad0(a,ilev,rp1,ct) - mean_vrad0(a,ilev,rm1,ct)
         dmp  = mean_pres0(a,ilev,rp1,ct) - mean_pres0(a,ilev,rm1,ct)
         vm   = mean_vtan0(a,ilev,irad,ct)
	 vmm  = mean_vtan0(a,ilev,irad,ct) ^ 2

       ; Calculate [ d{u}/dr ]
         mean_dudr0(a,ilev,irad,ct) = dmu / dr

       ; Calculate [ d{p}/dr ]
         mean_dpdr0(a,ilev,irad,ct) = dmp / dr

       ; Calculate [ {v}^2/r ]
         mean_vttr0(a,ilev,irad,ct) = vmm / r1

       end do                         ; End loop over radii

       mean_pres0 = mean_pres0 / 100

    ;=========================================================
    ; End loop over latitude and longitude ('box' ensemble)
    ;=========================================================

      end do  		            ; End latitude loop (do j = 0, ...)
      end do			    ; End longitude loop (do i = 0, ...)
 
    ; Tidy up
      delete([/u_int,v_int,w_int,v_tan,v_rad/])
      delete([/prs_int,rho_int,difx_int,dify_int,diff_int,f_int/])
      delete([/u_plane,v_plane,w_plane,f_plane,vort_plane,avo_plane,geo_plane/])
      delete([/lat_plane,lat1_plane,lon_plane,lon1_plane,difv_xplane,difv_yplane/])
      delete([/lat_sub,lon_sub,lat1_sub,lon1_sub/])

     end do			    ; End pressure level loop (do ilev = 0, num_p-1)

    ct  = ct + 1		    ; Counter variable (time)
    ct0 = ct0 + 1		    ; Independent counter variable (time)

  ; Tidy up before next iteration
    delete([/year,month,day,hour,minute,second,utc_date/])
    delete([/lon,lat,lon1,lat1,time,times,u,v,z,w,vort,avo,f0,f1,pres,plevs/])
    delete([/lonF,latF,lon1F,lat1F,hyb,str_xz,str_yz,rho1,th_ht,p_ht,t_ht/])

   end do     ; End time loop (do it = 0, times-1)

;================================================
; Calculate vertical derivatives of mean terms 
;================================================

    ; Change units before calculations
      pres  = b->p                     ; Pressure (levels)
      plevs = dimsizes(pres)           ; Size of pressure level array
      pres@units = "hPa"

      pres = pres * 100
      pres@units = "Pa"

     do a = 0, nbox-1 

      do it = ts0, tf0

        do ilev = 0, plevs-1

        kp1 = min((/ilev+1,plevs-1/))
        km1 = max((/ilev-1,0/))
 
          do irad = 0, dimsizes(radii)-1

          ;=====================================================
	  ; Calculate vertical derivative of mean radial wind
	  ;=====================================================

          ; [ d{u}/dz ]
	    mean_dudz0(a,ilev,irad,it)    = (mean_vrad0(a,kp1,irad,it) - \
                                             mean_vrad0(a,km1,irad,it)) / \
                                            (pres(km1) - pres(kp1))

          end do		 ; End loop over radial circles (do irad = 0, ...)

        end do 		         ; End loop over pressure levels (do ilev = 0, num_p-1)

       end do			 ; End loop over times (do it = 0, numTIMES-1)

      end do			 ; End loop over box ensemble members (do a = 0, nbox-1)

     ; Return units to previous values
       pres = pres / 100
       pres@units = "hPa"

;================================================
; Calculate vertical derivatives of eddy terms 
;================================================

    ; Change units before calculations
      pres  = b->p                     ; Pressure (levels)
      plevs = dimsizes(pres)           ; Size of pressure level array
      pres@units = "hPa"

      pres = pres * 100
      pres@units = "Pa"

     do a = 0, nbox-1

      do it = ts0, tf0

        do ilev = 0, plevs-1

        kp1 = min((/ilev+1,plevs-1/))
        km1 = max((/ilev-1,0/))

          do irad = 0, dimsizes(radii)-1

           do th0 = 0, dimsizes(thetas)-1


           ;=====================================================
           ; Calculate vertical derivative of eddy radial wind
           ;=====================================================

           ; [ du'/dz ]
             ptb_dudz0(a,ilev,irad,th0,it) = (ptb_vrad0(a,kp1,irad,th0,it) - \
                                              ptb_vrad0(a,km1,irad,th0,it)) / \
                                             (pres(km1) - pres(kp1))

           end do	         ; End loop over azimuth angles (do th0 = 0, ...)

          end do                 ; End loop over radial circles (do irad = 0, ...)

        end do                   ; End loop over pressure levels (do ilev = 0, num_p-1)

       end do                    ; End loop over times (do it = 0, numTIMES-1)

      end do                     ; End loop over box ensemble members (do a = 0, nbox-1)

     ; Return units to previous values
       pres = pres / 100
       pres@units = "hPa"

;==========================================================
; Finally, calculate source terms and azimuthal averages 
;==========================================================

     ;===================================================================
     ; Mean radial advection of radial momentum per unit mass: LHS (2)
     ;===================================================================
       lhs2_st = mean_vrad0 * mean_dudr0

     ;=============================================================
     ; Mean horizontal advection of eddy radial momentum: LHS (3)
     ;=============================================================
       ptb_ueh1 = (ptb_vrad0 * ptb_dudr0) + (ptb_vt_r0 * ptb_dudy0)

     ; Calculate the azimuthal average of 'ptb_ueh1' and input into 'ptb_ueh0'
       do a = 0, nbox-1
        ptb_ueh0(a,:,:,:) = dim_avg_n(ptb_ueh1(a,:,:,:,:),2) 
       end do 
       
       lhs3_st = ptb_ueh0       

     ;=========================================================================
     ; Mean vertical advection of mean radial momentum per unit mass: RHS (1)
     ;=========================================================================
       rhs1_st = -(mean_vvel0 * mean_dudz0)

     ;=========================================================================
     ; Eddy vertical advection of eddy radial momentum per unit mass: RHS (2)
     ;=========================================================================
       ptb_uev1 = -(ptb_vvel0 * ptb_dudz0)

     ; Calculate the azimuthal average of 'ptb_uev1' and input into 'ptb_uev0'
       do a = 0, nbox-1
        ptb_uev0(a,:,:,:) = dim_avg_n(ptb_uev1(a,:,:,:,:),2)
       end do

       rhs2_st = ptb_uev0

     ;==============================================
     ; Mean agradient force per unit mass: RHS (3)
     ;==============================================

     ; Account for elements where 'mean_rho0' is equal to zero (should be none)
       mean_rho0 = where(mean_rho0.ne.0, mean_rho0, mean_rho0@_FillValue)

       rhs3_st = mean_vttr0 + (mean_cor0 * mean_vtan0) - ( (1 / mean_rho0) * mean_dpdr0 )

     ;==============================================
     ; Eddy agradient force per unit mass: RHS (4)
     ;==============================================

     ; Account for elements where 'ptb_rho0' is equal to zero
       ptb_rho0  = where(ptb_rho0.ne.0,ptb_rho0,ptb_rho0@_FillValue)
     
     ; Now calculate eddy agradient force 
       ptb_uegf1 = ptb_vttr0 - ( (1 / ptb_rho0) * ptb_dpdr0 )

     ; Calculate the azimuthal average of 'ptb_uegf1' and input into 'ptb_uegf0'
       do a = 0, nbox-1
        ptb_uegf0(a,:,:,:) = dim_avg_n(ptb_uegf1(a,:,:,:,:),2)
       end do

       rhs4_st = ptb_uegf0

     ;===================================================================================
     ; Combined mean radial and vertical diffusive tendency of radial momentum: RHS (5)
     ;===================================================================================
     ; Define this array once we have added diffusive tendencies into script
     ; rhs5_st = ...

;=========================================
; Average over all box ensemble members 
;=========================================

     ; Define new arrays before calculating averages 
       lhs1_plot = new((/num_p,dimsizes(radii),numTIMES/),float)
       lhs2_plot = new((/num_p,dimsizes(radii),numTIMES/),float)
       lhs3_plot = new((/num_p,dimsizes(radii),numTIMES/),float)
       rhs1_plot = new((/num_p,dimsizes(radii),numTIMES/),float)
       rhs2_plot = new((/num_p,dimsizes(radii),numTIMES/),float)
       rhs3_plot = new((/num_p,dimsizes(radii),numTIMES/),float)
       rhs4_plot = new((/num_p,dimsizes(radii),numTIMES/),float)
       rhs5_plot = new((/num_p,dimsizes(radii),numTIMES/),float)

     ; Average over box ensemble members 
       lhs1_plot  = dim_avg_n_Wrap(mean_vrad0(:,:,:,:),0) ; Mean vrad tendency
       vtan_plot  = dim_avg_n_Wrap(mean_vtan0(:,:,:,:),0) ; Mean vtan tendency
       lhs2_plot  = dim_avg_n_Wrap(lhs2_st(:,:,:,:),0)	  ; LHS (2)
       lhs3_plot  = dim_avg_n_Wrap(lhs3_st(:,:,:,:),0)	  ; LHS (3)
;       lhs4_plot  = dim_avg_n_Wrap(ptb_ueh0(:,:,:,:),0)  ; LHS (4)
       rhs1_plot  = dim_avg_n_Wrap(rhs1_st(:,:,:,:),0)	  ; RHS (1)
       rhs2_plot  = dim_avg_n_Wrap(rhs2_st(:,:,:,:),0)	  ; RHS (2)
       rhs3_plot  = dim_avg_n_Wrap(rhs3_st(:,:,:,:),0)	  ; RHS (3)
       rhs4_plot  = dim_avg_n_Wrap(rhs4_st(:,:,:,:),0)	  ; RHS (4)
;       rhs5_plot  = dim_avg_n_Wrap(rhs5_st(:,:,:,:),0)	  ; RHS (5)

       vtan_plot!0    = "lev"
       vtan_plot&lev  = pres

       lhs1_plot!0    = "lev"
       lhs1_plot&lev  = pres

       lhs2_plot!0    = "lev"
       lhs2_plot&lev  = pres

       lhs3_plot!0    = "lev"
       lhs3_plot&lev  = pres

       rhs1_plot!0    = "lev"
       rhs1_plot&lev  = pres

       rhs2_plot!0    = "lev"
       rhs2_plot&lev  = pres

       rhs3_plot!0    = "lev"
       rhs3_plot&lev  = pres

       rhs4_plot!0    = "lev"
       rhs4_plot&lev  = pres

;       rhs5_plot!0    = "lev"
;       rhs5_plot&lev  = pres

;=================================================
; Tidy up and produce final arrays for plotting 
;=================================================

; Calculating the terms from Montgomery et al. (2018) - their Eq.(8)
; LHS (1) : d{u}/dt          - Local tendency of the mean radial velocity
; LHS (2) : {u} * d{u}/dr    - Mean radial advection of radial momentum
; LHS (3) : {(u' * du'/dr) + (v'/r * du'/dy)}  - Mean horiz. adv. of eddy radial momentum

; RHS (1) : -{w} * d{u}/dz   - Mean vertical advection of mean radial momentum
; RHS (2) : -{w' * d{u'}/dz} - Eddy vertical advection of eddy radial momentum
; RHS (3) : {v}^2/r + f{v} - 1/{rho} * d{p}/dr - Mean agradient force
; RHS (4) : {v'^2/r - 1/rho * dp'/dr}          - Eddy agradient force
; RHS (5) : {F_r}            - Mean radial/vertical diffusive tendency of radial momentum

     ; Define new arrays for terms in momentum budget 
       lhs1A = new((/num_p,dimsizes(radii),numTIMES/),float)   ; Local vrad tendency
       lhs2A = new((/num_p,dimsizes(radii),numTIMES/),float)   ; Mean r-adv of radial mom
       lhs3A = new((/num_p,dimsizes(radii),numTIMES/),float)   ; Mean h-adv of eddy mom

       rhs1A = new((/num_p,dimsizes(radii),numTIMES/),float)   ; Mean v-adv of mean mom
       rhs2A = new((/num_p,dimsizes(radii),numTIMES/),float)   ; Eddy v-adv of eddy mom
       rhs3A = new((/num_p,dimsizes(radii),numTIMES/),float)   ; Mean agradient force 
       rhs4A = new((/num_p,dimsizes(radii),numTIMES/),float)   ; Eddy agradient force
       rhs5A = new((/num_p,dimsizes(radii),numTIMES/),float)   ; Mean rad/vert diffusion
       rhsTA = new((/num_p,dimsizes(radii),numTIMES/),float)   ; Sum of RHS terms
       resA  = new((/num_p,dimsizes(radii),numTIMES/),float)   ; Residual (LHS - RHS)

     ; Calculate mean radial wind tendency (LHS #1)
       do it = ts0, tf0

         cp1 = min((/it+1,tf0/))
	 cm1 = max((/it-1,ts0/))
 
         lhs1A(:,:,it) = (lhs1_plot(:,:,cp1) - lhs1_plot(:,:,cm1)) / 2

       end do	; End loop over times (do it = ts0, tf0)

     ; Add metadata and make units consistent for all terms
       lhs1A@units = "m s~S~-1~N~ h~S~-1~N~"

       lhs2A = lhs2_plot
       lhs2A = lhs2A * 3600 
       lhs2A@units = "m s~S~-1~N~ h~S~-1~N~"

       lhs3A = lhs3_plot
       lhs3A = lhs3A * 3600
       lhs3A@units = "m s~S~-1~N~ h~S~-1~N~"

;       lhs4A = lhs4_plot
;       lhs4A = lhs4A * 3600

       rhs1A = rhs1_plot
       rhs1A = rhs1A * 3600
       rhs1A@units = "m s~S~-1~N~ h~S~-1~N~"

       rhs2A = rhs2_plot
       rhs2A = rhs2A * 3600
       rhs2A@units = "m s~S~-1~N~ h~S~-1~N~"
       
       rhs3A = rhs3_plot
       rhs3A = rhs3A * 3600
       rhs3A@units = "m s~S~-1~N~ h~S~-1~N~"

       rhs4A = rhs4_plot
       rhs4A = rhs4A * 3600
       rhs4A@units = "m s~S~-1~N~ h~S~-1~N~"

;       rhs5A = rhs5_plot
;       rhs5A = rhs5A * 3600 
;       rhs5A@units = "s~S~-1~N~ h~S~-1~N~"

     ; Smooth data before plotting 
       lhs1 = smth9_Wrap(lhs1A, 0.5, 0.5, True)
       lhs2 = smth9_Wrap(lhs2A, 0.5, 0.5, True)
       lhs3 = smth9_Wrap(lhs3A, 0.5, 0.5, True)
;       lhs4 = smth9_Wrap(lhs4A, 0.5, 0.5, True)

       rhs1 = smth9_Wrap(rhs1A, 0.5, 0.5, True)
       rhs2 = smth9_Wrap(rhs2A, 0.5, 0.5, True)
       rhs3 = smth9_Wrap(rhs3A, 0.5, 0.5, True)
       rhs4 = smth9_Wrap(rhs4A, 0.5, 0.5, True)
;       rhs5 = smth9_Wrap(rhs5A, 0.5, 0.5, True)

     ; Calculate total RHS tendency
       rhsT = rhs1 + rhs2 + rhs3 + rhs4 ; + rhs5
       rhsT@units = "m s~S~-1~N~ h~S~-1~N~"

     ; Calculate total LHS tendency
       lhsT = lhs1 + lhs2 + lhs3
       lhsT@units = "m s~S~-1~N~ h~S~-1~N~"

;     ; Calculate residual (LHS terms minus RHS terms)
;       res  = (lhs1 + lhs2 + lhs3) - rhsT
;       res@units = "m s~S~-1~N~ h~S~-1~N~"

       lhs1!0   = "lev"
       lhs1&lev = pres
       lhs1!1   = "rad"
       lhs1&rad = radii
       lhs1!2   = "time"

       lhs2!0   = "lev"
       lhs2&lev = pres
       lhs2!1   = "rad"
       lhs2&rad = radii
       lhs2!2   = "time"

       lhs3!0   = "lev"
       lhs3&lev = pres
       lhs3!1   = "rad"
       lhs3&rad = radii
       lhs3!2   = "time"

;       lhs4!0   = "lev"
;       lhs4&lev = pres
;       lhs4!1   = "rad"
;       lhs4&rad = radii
;       lhs4!2   = "time"

       rhs1!0   = "lev"
       rhs1&lev = pres
       rhs1!1   = "rad"
       rhs1&rad = radii
       rhs1!2   = "time"

       rhs2!0   = "lev"
       rhs2&lev = pres
       rhs2!1   = "rad"
       rhs2&rad = radii
       rhs2!2   = "time"

       rhs3!0   = "lev"
       rhs3&lev = pres
       rhs3!1   = "rad"
       rhs3&rad = radii
       rhs3!2   = "time"

       rhs4!0   = "lev"
       rhs4&lev = pres
       rhs4!1   = "rad"
       rhs4&rad = radii
       rhs4!2   = "time"

;       rhs5!0   = "lev"
;       rhs5&lev = pres
;       rhs5!1   = "rad"
;       rhs5&rad = radii
;       rhs5!2   = "time"

       rhsT!0   = "lev"
       rhsT&lev = pres
       rhsT!1   = "rad"
       rhsT&rad = radii
       rhsT!2   = "time"

       lhsT!0   = "lev"
       lhsT&lev = pres
       lhsT!1   = "rad"
       lhsT&rad = radii
       lhsT!2   = "time"

;       res!0    = "lev"
;       res&lev  = pres
;       res!1    = "rad"
;       res&rad  = radii
;       res!2    = "time"

;============================================
; Define output file path and plot details 
;============================================

 do it = ts0, tf0

    print("Working on time: "+it)

  ; Output file location and type
    output = "$sam/nepartak/images/cyl_coords/mbgt_rad_xz_"+ens0+\
    	     "_box"+nbox+"_"+time_arr(it)
    wks = gsn_open_wks(opt,output)

  ; Set up panel plot (9 panels)
    panel = new(9,graphic)
    if (clr.eq.0) then    
     gsn_define_colormap(wks,"BlueDarkRed18")
    else if (clr.eq.1) then
     gsn_define_colormap(wks,"bwr_new")
    end if 
    end if

;==========================
; Options for plotting 
;==========================

  ; Momentum budget terms 
    opts_mbgt                             = True
    opts_mbgt@cnFillOn                    = True
    opts_mbgt@cnLineLabelInterval         = 2.0
    opts_mbgt@cnLineLabelFontHeightF      = 0.012
    opts_mbgt@cnLineLabelBackgroundColor  = "transparent"
    opts_mbgt@cnLineLabelPlacementMode    = "constant"
    opts_mbgt@cnLinesOn                   = False ; Contour lines off
    opts_mbgt@cnInfoLabelOn               = False
    opts_mbgt@cnLevelSelectionMode        = "ExplicitLevels"

    if (clr.eq.0) then					       
     opts_mbgt@cnFillColors                = (/2,4,5,6,7,8,9,0,0,12,13,14,15,16,17,19/)
     opts_mbgt@cnLevels                    = (/-10., -5., -3., -2., -1., -0.5, -0.2, \
                                                0, 0.2, 0.5, 1., 2., 3., 5., 10./)
    else if (clr.eq.1) then
     opts_mbgt@cnFillColors		   = (/5,6,7,8,11,12,13,14,15,-1,-1,\
     					       19,21,22,23,25,28,29,30,31/)
     opts_mbgt@cnLevels                    = (/-20., -15., -10., -5., \
     					       -3., -2., -1., -0.5, -0.2, \
                                                0, 0.2, 0.5, 1., 2., 3., \
						5., 10., 15., 20./)
    end if 
    end if

    opts_mbgt@gsnPaperOrientation         = "landscape"
    opts_mbgt@tiMainString                = ""
    opts_mbgt@tiMainFontHeightF           = 0.0125
    opts_mbgt@gsnLeftString               = ""
    opts_mbgt@gsnRightString              = ""
    opts_mbgt@gsnMaximize                 = True
    opts_mbgt@lbLabelBarOn                = False
    opts_mbgt@gsnDraw                     = False
    opts_mbgt@gsnFrame                    = False

  ; Additional plotting resources
    opts_mbgt@tiYAxisString                = "Pressure (hPa)"
    opts_mbgt@tiXAxisString                = "Radial distance from cyclone centre (km)"
    opts_mbgt@tiXAxisFontHeightF           = 0.025   ; X-axis title font height 
    opts_mbgt@trYLog                       = False   ; Log scale on/off

    opts_mbgt@tmXTOn                       = "False" ; Turn off top x-axis TM
    opts_mbgt@tmYROn                       = "False" ; Turn off right y-axis TM
    opts_mbgt@gsnMaximize                  = True    ; Maximise plot size
    opts_mbgt@gsnAddCyclic                 = False
    opts_mbgt@tmXBMode                     = "Explicit"  ; Set tick marks explicitly
    opts_mbgt@tmXBValues                   = (/0.0,0.25,0.5,0.75,1.0,1.25,1.5/)	; Tm pos
    opts_mbgt@tmXBLabels                   = (/"0","25","50","75","100","125","150"/)
    opts_mbgt@tmXBLabelFontHeightF	   = 0.02    ; Label font height 

  ; Resources to overlay zero line
    opts_cont                             = True
    opts_cont@cnFillOn                    = False
    opts_cont@cnLineColor                 = "Black"
    opts_cont@cnLevels                    = (/-1.0, 0.0, 1.0/)
    opts_cont@cnLineThicknessF            = 3.0
    opts_cont@cnInfoLabelOn               = False
    opts_cont@gsnDraw                     = False       ; Do not draw the plot
    opts_cont@gsnFrame                    = False       ; Do no advance the frame
    opts_cont@gsnContourZeroLineThicknessF = 3.0
    opts_cont@cnLineLabelsOn               = False      ; Turn off line labels
    opts_cont@tiMainString                = ""
    opts_cont@gsnLeftString               = ""
    opts_cont@gsnRightString              = ""

  ; Resources to overlay radial wind (following M18)
    opts_vrad				  = True
    opts_vrad@cnFillOn			  = False
    if (clr.eq.0) then
     opts_vrad@cnLineColor		  = "black"
    else if (clr.eq.1) then
     opts_vrad@cnLineColor		  = "black"
    end if 
    end if 
    opts_vrad@cnLevelSelectionMode	  = "ExplicitLevels"
    opts_vrad@cnLevels                    = (/-10., -5./)
    opts_vrad@cnInfoLabelOn               = False
    opts_vrad@cnLineThicknessF		  = 2.5
    opts_vrad@gsnDraw			  = False
    opts_vrad@gsnFrame			  = False
    opts_vrad@cnLineLabelsOn		  = False
    opts_vrad@tiMainString                = ""
    opts_vrad@gsnLeftString               = ""
    opts_vrad@gsnRightString              = ""

;=================
; Plot the data 
;=================

    lhs1_panel = gsn_csm_pres_hgt(wks,lhs1(:,:,it),opts_mbgt)     ; LHS (1)
    lhs1_plotF = gsn_csm_contour(wks,lhs1_plot(:,:,it),opts_vrad) ; Radial wind
    overlay(lhs1_panel,lhs1_plotF)				  ; Overlay
    delete(lhs1_plotF)

    lhs2_panel = gsn_csm_pres_hgt(wks,lhs2(:,:,it),opts_mbgt)     ; LHS (2)
    lhs1_plotF = gsn_csm_contour(wks,lhs1_plot(:,:,it),opts_vrad) ; Radial wind
    overlay(lhs2_panel,lhs1_plotF)                                ; Overlay
    delete(lhs1_plotF)

    lhs3_panel = gsn_csm_pres_hgt(wks,lhs3(:,:,it),opts_mbgt)     ; LHS (3)
    lhs1_plotF = gsn_csm_contour(wks,lhs1_plot(:,:,it),opts_vrad) ; Radial wind
    overlay(lhs3_panel,lhs1_plotF)                                ; Overlay
    delete(lhs1_plotF)

    rhsT_panel = gsn_csm_pres_hgt(wks,rhsT(:,:,it),opts_mbgt)     ; Total RHS
    lhs1_plotF = gsn_csm_contour(wks,lhs1_plot(:,:,it),opts_vrad) ; Radial wind
    overlay(rhsT_panel,lhs1_plotF)                                ; Overlay
    delete(lhs1_plotF)

    lhsT_panel = gsn_csm_pres_hgt(wks,lhsT(:,:,it),opts_mbgt)     ; Total LHS
    lhs1_plotF = gsn_csm_contour(wks,lhs1_plot(:,:,it),opts_vrad) ; Radial wind
    overlay(lhsT_panel,lhs1_plotF)                                ; Overlay
    delete(lhs1_plotF)

    rhs1_panel = gsn_csm_pres_hgt(wks,rhs1(:,:,it),opts_mbgt)     ; RHS (1)
    lhs1_plotF = gsn_csm_contour(wks,lhs1_plot(:,:,it),opts_vrad) ; Radial wind
    overlay(rhs1_panel,lhs1_plotF)                                ; Overlay
    delete(lhs1_plotF)

    rhs2_panel = gsn_csm_pres_hgt(wks,rhs2(:,:,it),opts_mbgt)     ; RHS (2)
    lhs1_plotF = gsn_csm_contour(wks,lhs1_plot(:,:,it),opts_vrad) ; Radial wind
    overlay(rhs2_panel,lhs1_plotF)                                ; Overlay
    delete(lhs1_plotF)

    rhs3_panel = gsn_csm_pres_hgt(wks,rhs3(:,:,it),opts_mbgt)     ; RHS (3)
    lhs1_plotF = gsn_csm_contour(wks,lhs1_plot(:,:,it),opts_vrad) ; Radial wind
    overlay(rhs3_panel,lhs1_plotF)                                ; Overlay
    delete(lhs1_plotF)

    rhs4_panel = gsn_csm_pres_hgt(wks,rhs4(:,:,it),opts_mbgt)     ; RHS (4)
    lhs1_plotF = gsn_csm_contour(wks,lhs1_plot(:,:,it),opts_vrad) ; Radial wind
    overlay(rhs4_panel,lhs1_plotF)                                ; Overlay
    delete(lhs1_plotF)

    panel(0) = lhs1_panel
    panel(1) = lhs2_panel
    panel(2) = lhs3_panel
    panel(3) = rhsT_panel
    panel(4) = lhsT_panel
    panel(5) = rhs1_panel
    panel(6) = rhs2_panel
    panel(7) = rhs3_panel
    panel(8) = rhs4_panel

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    optsP                       = True			   ; Panel plot resources
    optsP@gsnFrame              = False			   ; Do not advance the frame
    optsP@gsnPanelLabelBar      = True       	           ; Panel labelbar
    optsP@pmLabelBarWidthF	= 0.6			   ; Label bar width 
    optsP@pmLabelBarHeightF     = 0.15	     	           ; Label bar height
    optsP@lbLabelFontHeightF	= 0.01	     	           ; Label bar font height 
    optsP@lbLabelFont		= "Helvetica" 	           ; Choose font 
    optsP@lbPerimOn		= False	      	           ; Perimeter on label bar 
    optsP@lbTitleString		= "m s~S~-1~N~ h~S~-1~N~"  ; Title string
    optsP@lbTitlePosition	= "Bottom"     		   ; Title position
    optsP@lbTitleFontHeightF	= 0.015			   ; Title font height
    optsP@lbBoxEndCapStyle	= "TriangleBothEnds"	   ; Label bar style
    optsP@pmLabelBarOrthogonalPosF = -0.002		   ; Move lb up/down (h/l)
    optsP@lbTitleOffsetF        = 0		   	   ; x-axis title up/down (h/l)
    optsP@txString              = title_arr(it)		   ; Plot title string
    optsP@gsnPanelFigureStrings = (/"a) U_t", "b) U_mr", "c) U_eh", \
                                    "d) All u tend", "e) U_t + U_mr + U_eh", \
				    "f) U_mv", "g) U_ev", \
				    "h) U_magf", "i) U_eagf"/)
    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.90
    optsP@gsnPanelBottom        = 0.10
    optsP@amJust                = "TopRight"
    optsP@gsnPanelFigureStringsFontHeightF = 0.01 ; Label size (default 0.01)

    gsn_panel(wks,panel,(/3,3/),optsP)             ; Draw as a single plot
    frame(wks)

   end do     ; End time loop (do it = 0, ts, tf)

end