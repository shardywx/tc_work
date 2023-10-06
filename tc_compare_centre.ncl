; Script to calculate the tangential momentum budget as in Montgomery et al. (2018)

; Run using:

; ncl dat=\"02T12\" opt=\"x11\" ens0=\"em11\" dist=3.0 ts0=54 tf0=63 clr=5
; full=0 rad0=1 ave=0 wrt=0 new0=0 mlev0=13 mlev1=34 var0=0 lgd0=\"tl\"
; cn0=\"slp\" tc_compare_centre.ncl

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'opt'   = output file format ("pdf" or "x11")
; 'ens0'  = ensemble member (em00 -- em11)
; 'dist'  = size of box following storm (degrees)
; 'ts0'   = starting time index (0-119)
; 'tf0'   = ending time index (0-119)
; 'clr'   = colour map for plots ("ncl_default" [1];  "bdr_extra" [2]; "amwg256" [3])
; 'full'  = plot vertical velocity on two levels on same (1) or separate (0) panels 
; 'rad0'  = radius of final plots: 150 km (0) or 200 km (1)
; 'ave'   = calculate 3/6/12-h average (1) rather than analyse individual time interval (0)
; 'wrt'   = write out averages to text file (1) or leave alone (0 default)
; 'new0'  = reduced (1) or original (0) scale for plotting u,v,w,AAM
; 'pnum'  = number of panels in final plot
; 'mlev0' = model level for final plot (15, 23, 36, etc)
; 'mlev1' = additional level for final plot (if required) 
; 'var0'  = instantaneous (0) or perturbation (1) fields for final plot

; 'pc' stream --> 'it'
; 'centre'    --> 'it-1'
; 'pd' stream --> 'it-2'

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
; 'dat'     = initialisation time (02T12, 03T00, ...)
; 'ens0'    = ensemble simulation (em00, em01, ...)
; 'diri'    = path to input files (see above)
; 'dist'    = size of box following storm (degrees)
; 'mins'    = analyse 1-h (0) or 5-min (1) data

; Call external function 'st_rm'
  storm_rel = st_rm(dat, ens0, "$ar/text/", dist, 0)
  u_cyc     = storm_rel[0]		   ; Zonal wind 
  v_cyc     = storm_rel[1]		   ; Meridional wind
  vel_cyc   = storm_rel[2]		   ; Vector wind 
  numTIMES  = storm_rel[3]                 ; Number of times in file
  lat_arr   = storm_rel[4]                 ; Storm-following domain subset (lat)
  lon_arr   = storm_rel[5]                 ; Storm-following domain subset (lon)
  centre    = storm_rel[6]		   ; Storm centre at each hour (lat/lon)

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
; Create date/time string arrays for all times in file (use 'pb' stream)
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

  avo_int0@description = "Relative vorticity"
  avo_int0@units    = "s~S~-1~N~"
  avo_int0!0           = "lev"
  avo_int0!1           = "rad"
  avo_int0!2           = "azi"
  avo_int0!3           = "time"

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

;==========================================================================
; Also create array to hold all storm centre values for later annotation
;==========================================================================

; Dimension 0 [centre method]; dimension 1 [lat/lon]; dimension 2 [time]
  storm_centre     = new((/num_l,4,2,numTIMES/),float)
  
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

   c  = addfile(fili_c,"r")         	 ; Read in 'pc' stream [u,v,vort]
   d  = addfile(fili_d,"r")         	 ; Read in 'pd' stream [w,p,t,theta,kmh,tau]

;==================================
; Get the variables we will need
;==================================

   do it = ts0, tf0		; Loop over times in file

    print("Working on time: "+time_str(it))
    itt = it + 1
    title_arr(it) = "Valid at "+time_str(it)+" (T+"+itt+")"

;=====================================================
; Read in variables from 'pc' stream (model levels)
;=====================================================

    t0    = lat_arr(cn,it-1,0)
    t1    = lat_arr(cn,it-1,1)
    n0    = lon_arr(cn,it-1,0)
    n1    = lon_arr(cn,it-1,1)
    print_clock("Grid (pc stream) extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

  ; Read in basic variables from 'pc' stream
    lon1C  = c->longitude_1({n0:n1})              ; longitude ['d1' grid points]
    lat1C  = c->latitude_1({t0:t1})               ; latitude  ['d0' grid points]
    lonC   = c->longitude({n0:n1})                ; longitude ['d1' grid points]
    latC   = c->latitude({t0:t1})                 ; latitude  ['d0' grid points]
    hybC   = c->hybrid_ht(:) 		          ; 63 model ('theta') levels

    hy0	   = sprintf("%0.0f",hybC(mlev0))
    hy1	   = sprintf("%0.0f",hybC(mlev1))

    u_varname = "x-wind"
    v_varname = "y-wind"

  ; Read in horizontal wind components 
  ; Both variables below -- [hybC | 63] * [latC | 150] * [lonC | 150]
    u0  = c->$u_varname$(it,:,{t0:t1},{n0:n1}) ; Zonal wind	
    v0  = c->$v_varname$(it,:,{t0:t1},{n0:n1}) ; Meridional wind 

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

;==================================================================
; Calculate relative vorticity using centered finite differences
;==================================================================

  ; Option '2' --> boundary points estimated using one-sided difference scheme
    vort = uv2vr_cfd(u, v, latC, lonC, 2)
    copy_VarCoords(u, vort)
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

;======================================================
; Read in variables from 'pd' stream (model levels) 
;======================================================

  ; Read in basic variables from 'pd' stream 
    lon1D  = d->longitude_1({n0:n1})                ; longitude [xxx grid points]
    lat1D  = d->latitude_1({t0:t1})                 ; latitude  [xxx grid points]
    lonD   = d->longitude({n0:n1})                  ; longitude [xxx grid points]
    latD   = d->latitude({t0:t1})                   ; latitude  [xxx grid points]
    hybD   = d->hybrid_ht(:)			    ; 63 model ('rho') levels

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

;=======================================
; Calculate storm centre for 'calc=2'
;=======================================

  ; Model level closest to 1 km AGL (980 m)
    ilev = 13

  ; Read in relative vorticity on chosen model level
    vort_plane0 = vort(ilev,:,:)

  ; Calculate location of vorticity maximum (storm centre) on chosen model level
    vort_max0               = max(vort_plane0)
    dims0                   = dimsizes(vort_plane0)
    vort1d0                 = ndtooned(vort_plane0)
    inds0                   = ind_resolve(maxind(vort1d0),dims0)
    delete([/vort_plane0, vort_max0, dims0, vort1d0/])

;==========================
; Loop over model levels 
;==========================

    do ilev = 0, num_l-1 ; Loop over all pressure levels

     hy = sprintf("%0.0f",hybD(ilev))
     print("Working on model level " + ilev+ " ("+hy+" m AGL)")

   ; Read in other variables on model levels
     u_plane 	            = u(ilev,:,:)	   ; Zonal wind 
     v_plane		    = v(ilev,:,:)	   ; Meridional wind 
     w_plane		    = w(ilev,:,:)	   ; Vertical velocity
     f_plane		    = f0(ilev,:,:)	   ; Coriolis parameter
     prs_plane		    = p(ilev,:,:)	   ; Pressure 
     vort_plane	   	    = vort(ilev,:,:)	   ; Relative vorticity
     avo_plane		    = avo(ilev,:,:)	   ; Absolute vorticity

   ; Calculate location of vorticity maximum (storm centre) on model level
     vort_max      	    = max(vort_plane)
     vort_smth     	    = smth9_Wrap(vort_plane, 0.5, 0.5, True)
     vort_max_smth 	    = max(vort_smth)
   
   ; Find index of relative vorticity maximum
     dims          	    = dimsizes(vort_plane)
     vort1d        	    = ndtooned(vort_plane)
     inds          	    = ind_resolve(maxind(vort1d),dims)
     vort1d_sm     	    = ndtooned(vort_smth)
     inds_sm       	    = ind_resolve(maxind(vort1d_sm),dims)

   ;=========================================================
   ; Calculate the storm centre using a variety of methods
   ;=========================================================

   ;============
   ; METHOD 1
   ;============

   ; Find latitude and longitude of storm centre using vorticity calculation above
     lat_max0 	            = t0 + (dy0 * inds(0,0))
     lon_max0 	   	    = n0 + (dx0 * inds(0,1))

   ; Fill larger array with position of storm centre at each time 
     storm_centre(ilev,0,0,it)   = lat_max0
     storm_centre(ilev,0,1,it)   = lon_max0

     print("Centre (method 1): "+lat_max0+" degrees N, "+lon_max0+" degrees E")

   ; Access 'lat_max0' and 'lon_max0' later when overlaying storm centre

   ;============
   ; METHOD 2
   ;============

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

     lat_max1   = lt1 + (dy0 * inds_p0(0,0))
     lon_max1   = ln1 + (dx0 * inds_p0(0,1))
     print("Centre (method 2): "+lat_max1+" degrees N, "+lon_max1+" degrees E")
     delete([/prs_plane0,prs_min0,dims_p0,prs_1d0,inds_p0/])

   ; Fill larger array with position of storm centre at each time
     storm_centre(ilev,1,0,it)   = lat_max1
     storm_centre(ilev,1,1,it)   = lon_max1

   ;============
   ; METHOD 3
   ;============

   ; Distance between centre calculated above using relative vorticity, and that calculated offline
     dlat       = abs(lat_max1 - centre(cn,it-1,0))
     dlon	= abs(lon_max1 - centre(cn,it-1,1))

   ; If our storm centre calculations deviate too much (> 0.75º) from the values
   ; we read in using the NCL function 'st_rm', calculate a new centre using
   ; the minimum pressure. Constrain the new calculation using the values from
   ; 'st_rm', to prevent spurious centres being identified and analysed. 

     if (dlat .gt. 1.00 .or. dlon .gt. 1.00) then 

      lt1 = centre(cn,it-1,0) - 0.25
      lt2 = centre(cn,it-1,0) + 0.25
      ln1 = centre(cn,it-1,1) - 0.25
      ln2 = centre(cn,it-1,1) + 0.25

      prs_plane0 = prs_plane({lt1:lt2},{ln1:ln2})       ; Smaller grid (1.0 x 1.0 deg)
      prs_min0   = min(prs_plane0)                      ; Find minimum on smaller grid
      dims_p0    = dimsizes(prs_plane0)                 ; Size of smaller grid
      prs_1d0    = ndtooned(prs_plane0)                 ; Create 1-D array
      inds_p0    = ind_resolve(minind(prs_1d0),dims_p0) ; Find index of minimum

      lat_max2 	 = lt1 + (dy0 * inds_p0(0,0))
      lon_max2	 = ln1 + (dx0 * inds_p0(0,1))
      print("Centre (method 3): "+lat_max2+" degrees N, "+lon_max2+" degrees E")
      delete([/prs_plane0,prs_min0,dims_p0,prs_1d0,inds_p0/])

;    ; Fill larger array with position of storm centre at each time
;      storm_centre(2,0,it)   = lat_max2
;      storm_centre(2,1,it)   = lon_max2

     end if 

   ;============
   ; METHOD 4
   ;============

     lat_max3    = centre(cn,it-1,0)
     lon_max3    = centre(cn,it-1,1)
     print("Centre (method 4): "+lat_max3+" degrees N, "+lon_max3+" degrees E")

   ; Fill larger array with position of storm centre at each time
     storm_centre(ilev,2,0,it)   = lat_max3
     storm_centre(ilev,2,1,it)   = lon_max3

   ;============
   ; METHOD 5
   ;============

     lat_max4    = t0 + (dy0 * inds0(0,0))
     lon_max4    = n0 + (dx0 * inds0(0,1))
     print("Centre (method 5): "+lat_max4+" degrees N, "+lon_max4+" degrees E")      

   ; Fill larger array with position of storm centre at each time
     storm_centre(ilev,3,0,it)   = lat_max4
     storm_centre(ilev,3,1,it)   = lon_max4

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
     f_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Coriolis force 
     vort_int = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Relative vorticity
     avo_int  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Absolute vorticity
     v_rad    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Radial velocity
     v_tan    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Tangential velocity

   ; Create cylindrincal coordinate array from original latitude/longitude grid
   ; Use storm centre calculated offline as starting point for these calculations 

     do irad = 0, dimsizes(radii)-1
       r = radii(irad)
       do iang = 0, dimsizes(thetas)-1
          theta = thetas_rad(iang)   

           xpos(irad,iang) = lon_max3 + r*cos(theta)
 	   ypos(irad,iang) = lat_max3 + r*sin(theta)
           xcen(irad,iang) = lon_max3 
           ycen(irad,iang) = lat_max3
	   
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
   ; vort_plane	 -- latC,  lonC
   ; avo_plane	 -- latC,  lonC
   ; f_plane     -- latC,  lonC

   ; w_plane     -- latD,  lon1D
   ; prs_plane   -- latD,  lon1D

     do iang = 0, dimsizes(thetas_rad)-1
      u_int(:,iang)    = linint2_points(lonC_sub,latC_sub,u_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      v_int(:,iang)    = linint2_points(lonC_sub,latC_sub,v_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      vort_int(:,iang) = linint2_points(lonC_sub,latC_sub,vort_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      avo_int(:,iang)  = linint2_points(lonC_sub,latC_sub,avo_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      f_int(:,iang)    = linint2_points(lonC_sub,latC_sub,f_plane,False,\
				          xpos(:,iang),ypos(:,iang),0)
      w_int(:,iang)    = linint2_points(lon1D_sub,latD_sub,w_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
     end do

   ; Calculate radial/tangential windspeed and friction at each point (using unit vctrs)
      do iang = 0, dimsizes(thetas)-1
       theta = thetas_rad(iang)

      ; Tangential and radial wind 
        v_rad(:,iang) = u_int(:,iang) * cos(theta) + v_int(:,iang) * sin(theta)
        v_tan(:,iang) = -u_int(:,iang) * sin(theta) + v_int(:,iang) * cos(theta)

      end do

    ;===================================================
    ; Calculate mean (azimuthally-averaged) variables 
    ;===================================================

      do irad = 0, dimsizes(radii)-1

        mean_vtan0(ilev,irad,:,ct)    = avg(v_tan(irad,:))      ; [ {v} ]
        mean_vrad0(ilev,irad,:,ct)    = avg(v_rad(irad,:))      ; [ {u} ]
        mean_vvel0(ilev,irad,:,ct)    = avg(w_int(irad,:))      ; [ {w} ]
	mean_vort0(ilev,irad,:,ct)    = avg(vort_int(irad,:))   ; [ {zeta} ]
        mean_avo0(ilev,irad,:,ct)     = avg(avo_int(irad,:))    ; [ {zeta + f} ]
        mean_cor0(ilev,irad,:,ct)     = avg(f_int(irad,:))      ; [ {f} ]

      end do

      fcor = mean_cor0 * mean_vtan0

    ; Read 'v_int' etc into larger arrays
      u_int0(ilev,:,:,ct)    = (/u_int(:,:)/)
      v_int0(ilev,:,:,ct)    = (/v_int(:,:)/)
      w_int0(ilev,:,:,ct)    = (/w_int(:,:)/)
      vort_int0(ilev,:,:,ct) = (/vort_int(:,:)/)
      avo_int0(ilev,:,:,ct)  = (/avo_int(:,:)/)
      f_int0(ilev,:,:,ct)    = (/f_int(:,:)/)
      v_rad0(ilev,:,:,ct)    = (/v_rad(:,:)/)
      v_tan0(ilev,:,:,ct)    = (/v_tan(:,:)/)

    ; Tidy up
      delete([/u_int,v_int,w_int,v_tan,v_rad,f_int,vort_int,avo_int/])
      delete([/u_plane,v_plane,w_plane,f_plane,vort_plane,avo_plane/])
      delete([/latC_plane,lonC_plane,latD_plane,lonD_plane,lon1D_plane,lat1D_plane/])
      delete([/latC_sub,lonC_sub,latD_sub,lonD_sub,lat1D_sub,lon1D_sub/])
      delete([/xpos,ypos,xcen,ycen/])

     end do			    ; End model level loop (do ilev = 0, num_l-1)

    ct  = ct + 1		    ; Counter variable (time)

  ; Tidy up before next iteration (all variables from 'pc' and 'pd' streams)
    delete([/lonC,latC,lon1C,lat1C,u,v,vort,avo,f0,f1/])
    delete([/lonD,latD,lon1D,lat1D,t,p,w,z/])

   end do     ; End time loop (do it = 0, times-1)

;====================================
; Calculate perturbation variables
;====================================

     ptb_vtan0(:,:,:,:) = v_tan0(:,:,:,:)    - mean_vtan0(:,:,:,:)    ; [v']
     ptb_vvel0(:,:,:,:) = w_int0(:,:,:,:)    - mean_vvel0(:,:,:,:)    ; [w']
     ptb_vrad0(:,:,:,:) = v_rad0(:,:,:,:)    - mean_vrad0(:,:,:,:)    ; [u']
     ptb_vort0(:,:,:,:) = vort_int0(:,:,:,:) - mean_vort0(:,:,:,:)    ; [zeta']

;============================================
; Define output file path and plot details 
;============================================

    do it = ts0, tf0
     print("Working on time: "+it)
     print(storm_centre(mlev0,:,:,it))
    end do 

  ; Before new loop, multiply vorticity ready for plotting 

  ; Relative vorticity
    vort_int0	     = vort_int0 * (10 ^ 4)
    vort_int0@units  = "10~S~-4~N~ s~S~-1~N~"

  ; Perturbation relative vorticity
    ptb_vort0 	     = ptb_vort0 * (10 ^ 4)
    ptb_vort0@units  = "10~S~-4~N~ s~S~-1~N~"

    do it = ts0, tf0             ; Loop over times in file

    ; Output file location and type 
      if (var0 .eq. 0) then
       output = "$sam/nepartak/images/cyl_coords/centre_xy_"+dat+"_"+ens0+"_"\
       	         +hy0+"m_"+hy1+"m_inst_"+time_arr(it)
      elseif (var0 .eq. 1) then 
       output = "$sam/nepartak/images/cyl_coords/centre_xy_"+dat+"_"+ens0+"_"\
                 +hy0+"m_"+hy1+"m_ptb_"+time_arr(it)
      end if 

    ; Set up workstation	     
      wks    = gsn_open_wks(opt,output)

    ; Create graphic for panel plot
      if (full .eq. 1) then 
       panel = new(6,graphic)
      elseif (full .eq. 0) then 
       panel = new(1,graphic)
      end if 

    ; Read in chosen colourmap
      if (var0 .eq. 0) then 
       gsn_define_colormap(wks,"prcp_new")
      elseif (var0 .eq. 1) then 
       gsn_define_colormap(wks,"BlueDarkRed18")
      end if 
 
;==========================
; Options for plotting 
;==========================

    ; Relative vorticity
      opts_vort                             = True
      opts_vort@cnFillOn                    = True
      opts_vort@sfXArray		    = xpos_all(mlev0,:,:,it)
      opts_vort@sfYArray		    = ypos_all(mlev0,:,:,it)
      opts_vort@cnLineLabelInterval         = 2.0
      opts_vort@cnLineLabelFontHeightF      = 0.012
      opts_vort@cnLineLabelBackgroundColor  = "transparent"
      opts_vort@cnLineLabelPlacementMode    = "constant"
      opts_vort@cnLinesOn                   = False
      opts_vort@cnInfoLabelOn               = False
      opts_vort@cnLevelSelectionMode        = "ExplicitLevels"

      opts_vort@cnLevels                    = (/5., 10., 20., 30., 40., 50., \
                                                60., 80., 100., 120., 150./)
      opts_vort@cnFillColors                = (/0,4,5,6,7,8,9,10,11,12,13,14/)

      opts_vort@gsnPaperOrientation          = "landscape"
      opts_vort@tiMainString                 = ""
      opts_vort@tiMainFontHeightF            = 0.0125
      opts_vort@gsnLeftString                = ""
      opts_vort@gsnRightString               = ""
      opts_vort@gsnMaximize                  = True
      opts_vort@lbLabelBarOn                 = False
      opts_vort@gsnDraw                      = False
      opts_vort@gsnFrame                     = False

      opts_vort@pmLabelBarWidthF             = 0.6
      opts_vort@pmLabelBarHeightF            = 0.15
      opts_vort@lbLabelFontHeightF           = 0.01
      opts_vort@lbPerimOn                    = False
      opts_vort@lbLabelFont                  = "Helvetica"
      opts_vort@lbTitleOn                    = True
      opts_vort@lbTitleString                = "Relative vorticity (10~S~-4~N~ s~S~-1~N~)"
      opts_vort@lbTitleFontHeightF           = 0.015
      opts_vort@lbTitleDirection             = "Across"

    ; Additional plotting resources
      opts_vort@mpDataBaseVersion            = "Ncarg4_1" ; More recent database
      opts_vort@mpDataSetName                = "Earth..4" ; High resolution
      opts_vort@mpOutlineBoundarySets        = "National" ; National borders
      opts_vort@mpGeophysicalLineColor       = "black"    ; Colour borders black
      opts_vort@mpGeophysicalLineThicknessF  = 1.0        ; Border line thickness
      opts_vort@mpGridAndLimbOn              = False      ; Turn on lat/lon lines
      opts_vort@pmTickMarkDisplayMode        = "Always"   ; Turn on map tickmarks
      opts_vort@tmXBMajorLengthF             = 0.005      ; Change tickmark length
      opts_vort@tmXTOn                       = "False"    ; No tickmarks on top x-axis
      opts_vort@tmYROn                       = "False"    ; No tickmarks on right y-axis
      opts_vort@gsnMaximize                  = True       ; Maximise plot size
      opts_vort@gsnAddCyclic                 = False

      radius                                 = 1.0        ; Radius of plot (degrees)

      opts_vort@mpLimitMode                  = "Corners"
      opts_vort@mpLeftCornerLatF             = centre(cn,it-1,0)-radius
      opts_vort@mpLeftCornerLonF             = centre(cn,it-1,1)-radius
      opts_vort@mpRightCornerLatF            = centre(cn,it-1,0)+radius
      opts_vort@mpRightCornerLonF            = centre(cn,it-1,1)+radius
      opts_vort@gsnDraw                      = False      ; Do not draw the plot
      opts_vort@gsnFrame                     = False      ; Do not advance the frame

    ; Perturbation plot resources 
      opts_ptb0				     = opts_vort
      delete(opts_ptb0@cnFillColors)
      delete(opts_ptb0@cnLevels)
      opts_ptb0@cnFillColors                 = (/2,4,5,6,7,8,9,0,0,12,13,14,15,16,17,18/)
      opts_ptb0@cnLevels		     = (/-10., -5., -3., -2., -1., -0.5, -0.2, \
                                                 0, 0.2, 0.5, 1., 2., 3., 5., 10./)

    ; Horizontal wind vector resources
      opts_vec                               = True
      opts_vec@gsnDraw                       = False
      opts_vec@gsnFrame                      = False
      opts_vec@gsnAddCyclic                  = False
 
      opts_vec@vcGlyphStyle                  = "LineArrow"   ; draw windbarbs
      opts_vec@vcRefMagnitudeF               = 15.0          ; ref vec = max vector
      opts_vec@vcRefLengthF                  = 0.04          ; length of reference vector
      opts_vec@vcMinDistanceF                = 0.08          ; distance between wind barbs
      opts_vec@vcLineArrowThicknessF         = 2.0           ; windbarb line thickness
      opts_vec@vcLineArrowColor              = "black"       ; windbarb colour
      opts_vec@gsnMaximize                   = True          ; Maximise plot size
      opts_vec@gsnRightString                = ""
      opts_vec@gsnLeftString                 = ""
      opts_vec@tiXAxisString                 = ""
      opts_vec@vcRefAnnoOrthogonalPosF       = 0.02          ; Y-location of ref vec

    ; Vertical velocity resources
      opts_vvel                              = True
      opts_vvel@cnFillOn                     = False
      opts_vvel@cnInfoLabelOn                = False
      opts_vvel@gsnDraw                      = False
      opts_vvel@gsnFrame                     = False
      opts_vvel@gsnAddCyclic                 = False
      opts_vvel@tiMainString                 = ""
      opts_vvel@gsnLeftString                = ""
      opts_vvel@gsnRightString               = ""
      opts_vvel@cnLineThicknessF             = 4.0
      opts_vvel@cnLineLabelsOn               = False
      opts_vvel@cnLevelSelectionMode         = "ExplicitLevels"
      opts_vvel@cnLevels                     = 1.0

    ; Yellow vertical velocity contours if we choose Roger Smith's colour map
      if (clr .eq. 4) then
       opts_vvel@cnLineColor                 = "yellow"
      elseif (clr .eq. 5) then
       opts_vvel@cnLineColor                 = "black"
      else
       opts_vvel@cnLineColor                 = "gray45"
      end if

    ; Resources for mid-tropospheric vertical velocity
      opts_vvel1                             = opts_vvel
      opts_vvel1@cnLineColor                 = "royalblue1"

    ; Resources for perturbation vertical velocity 
      opts_ptb1	    		 	     = opts_vvel
      delete(opts_ptb1@cnLevels)
      opts_ptb1@cnLevels		     = (/-0.20, -0.10, -0.05, 0.05, 0.10, 0.20/)

      opts_ptb2				     = opts_vvel1
      delete(opts_ptb2@cnLevels)
      opts_ptb2@cnLevels		     = (/-0.20, -0.10, -0.05, 0.05, 0.10, 0.20/)

;=============
; Plot data
;=============

    if (var0 .eq. 0) then				; Instantaneous variables

    ; Relative vorticity
      vort_plane0      = vort_int0(mlev0,:,:,it)
      vort_plane1      = vort_int0(mlev1,:,:,it)
    ; Vertical velocity
      vvel_plane0      = w_int0(mlev0,:,:,it)
      vvel_plane1      = w_int0(mlev1,:,:,it)
    ; Tangential wind
      vtan_plane0      = v_tan0(mlev0,:,:,it)
      vtan_plane1      = v_tan0(mlev1,:,:,it)
    ; Radial wind 
      vrad_plane0      = v_rad0(mlev0,:,:,it)
      vrad_plane1      = v_rad0(mlev1,:,:,it)

    elseif (var0 .eq. 1) then				; Perturbation variables 

    ; Relative vorticity
      vort_plane0      = ptb_vort0(mlev0,:,:,it)
      vort_plane1      = ptb_vort0(mlev1,:,:,it)
    ; Vertical velocity
      vvel_plane0      = ptb_vvel0(mlev0,:,:,it)
      vvel_plane1      = ptb_vvel0(mlev1,:,:,it)
    ; Tangential wind
      vtan_plane0      = ptb_vtan0(mlev0,:,:,it)
      vtan_plane1      = ptb_vtan0(mlev1,:,:,it)
    ; Radial wind 
      vrad_plane0      = ptb_vrad0(mlev0,:,:,it)
      vrad_plane1      = ptb_vrad0(mlev1,:,:,it)

    end if 

;=====================
; Lower model level
;=====================

  ; Plot instantaneous or perturbation quantities
    if (var0 .eq. 0) then 
     plot_vort0 = gsn_csm_contour_map(wks,vort_plane0,opts_vort)	; Relative vorticity (inst)
     plot_vvel0 = gsn_csm_contour(wks,vvel_plane0,opts_vvel)		; Vertical velocity (inst)
    elseif (var0 .eq. 1) then
     plot_vort0 = gsn_csm_contour_map(wks,vort_plane0,opts_ptb0)	; Relative vorticity (ptb)
     plot_vvel0 = gsn_csm_contour(wks,vvel_plane0,opts_ptb1)		; Vertical velocity (ptb)
    end if 

    overlay(plot_vort0,plot_vvel0)

;=====================
; Upper model level
;=====================

  ; Plot instantaneous or perturbation quantities
    if (var0 .eq. 0) then 
     plot_vort1 = gsn_csm_contour_map(wks,vort_plane1,opts_vort)	; Relative vorticity (inst)
     plot_vvel1 = gsn_csm_contour(wks,vvel_plane1,opts_vvel)		; Vertical velocity (inst)
    elseif (var0 .eq. 1) then 
     plot_vort1 = gsn_csm_contour_map(wks,vort_plane1,opts_ptb0)	; Relative vorticity (ptb)
     plot_vvel1 = gsn_csm_contour(wks,vvel_plane1,opts_ptb2)		; Vertical velocity (ptb)
    end if 

  ; Plot either one or two panel(s)

    if (full .eq. 0) then                  ; One panel

     overlay(plot_vort0, plot_vvel1)	   ; Also overlay upper-level vertical velocity
     panel(0) = plot_vort0

    elseif (full .eq. 1) then              ; Two panels

     panel(0) = plot_vort0		   ; Panel #1 
     overlay(plot_vort1, plot_vvel1)
     panel(1) = plot_vort1		   ; Panel #2

    end if

;=========================================================
; Overlay storm centre positions from different methods
;=========================================================

  ; Method 1 -- Variable level relative vorticity 
  ; Method 2 -- Variable level relative vorticity and pressure
  ; Method 3 -- Offline
  ; Method 4 -- Fixed-level relative vorticity (980 m)

    m_colours  = (/"grey55", "royalblue1", "black", "purple3"/)
    m_indices  = (/6, 16, 12, 9/)
    m_sizes    = (/10.0, 10.0, 10.0, 10.0/)
    m_thick    = (/5.0, 5.0, 5.0, 5.0/)

    markers    = new(dimsizes(m_colours),graphic)

    do i = 0, dimsizes(m_colours)-1
     mres0                    = True
     mres0@gsMarkerIndex      = m_indices(i)
     mres0@gsMarkerSizeF      = m_sizes(i)
     mres0@gsMarkerColor      = m_colours(i)
     mres0@gsMarkerThicknessF = m_thick(i)
     markers(i)               = gsn_add_polymarker(wks,panel(0),storm_centre(mlev0,i,1,it),\
                                                   storm_centre(mlev0,i,0,it),mres0)
    end do

;===============================
; Add and customise a legend
;===============================

    lg_opts                            = True
    lg_opts@pmLegendDisplayMode        = "Always"          ; Display a legend
    lg_opts@lgAutoManage               = False             ; Gives us control
    lg_opts@lgLabelFont                = "Helvetica"       ; Label font
    lg_opts@lgLabelFontHeightF         = 0.050             ; Height of labels

  ; Customise labels and line colours
    lg_opts@lgItemType		       = "Markers"
    lg_opts@lgMarkerColors	       = m_colours
    lg_opts@lgMarkerIndexes	       = m_indices
    lg_opts@lgMarkerSizeF	       = 10.0
    lg_opts@lgMarkerThicknessF	       = 5.0
    lsize                              = dimsizes(lg_opts@lgMarkerColors)

    lg_opts@lgLabelPosition            = "Right"           ; Label position within box
    lg_opts@lgItemPlacement            = "ExplicitPlacement"
    lg_opts@lgItemPositions            = fspan(0.20, 0.80, lsize)

    lg_opts@vpWidthF                   = 0.250
    lg_opts@vpHeightF                  = 0.150
    lg_opts@lgPerimColor               = "black"
    lg_opts@lgPerimThicknessF          = 3.0
    lg_opts@lgPerimFill                = "SolidFill"
    lg_opts@lgPerimFillColor           = "white"

    lab_arr                            = (/"vort ("+hy0+" m)", "vort + pres ("+hy0+" m)",\
                                           "mslp", "vort (980 m)"/)
    lgd                                = gsn_create_legend(wks,lsize,lab_arr,lg_opts)

  ; Add legend to plot (top left or right hand corner)
    am_opts                            = True

    if (lgd0 .eq. "tr") then
     am_opts@amJust                    = "TopRight"
     am_opts@amParallelPosF            = 0.5
     am_opts@amOrthogonalPosF          = -0.5
    elseif (lgd0 .eq. "tl") then
     am_opts@amJust                    = "TopLeft"
     am_opts@amParallelPosF            = -0.5
     am_opts@amOrthogonalPosF          = -0.5
    elseif (lgd0 .eq. "br") then
     am_opts@amJust                    = "BottomRight"
     am_opts@amParallelPosF            = 0.5
     am_opts@amOrthogonalPosF          = 0.5
    elseif (lgd0 .eq. "bl") then
     am_opts@amJust                    = "BottomLeft"
     am_opts@amParallelPosF            = -0.5
     am_opts@amOrthogonalPosF          = 0.5
    end if

    annotate_0                         = gsn_add_annotation(panel(0),lgd,am_opts)
    annotate_1                         = gsn_add_annotation(panel(1),lgd,am_opts)

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
    optsP@lbTitleString         = "Relative vorticity (10~S~-4~N~ s~S~-1~N~)"
    optsP@lbTitlePosition       = "Bottom"                 ; Title position
    optsP@lbTitleFontHeightF    = 0.015                    ; Title font height
    optsP@lbBoxEndCapStyle      = "TriangleBothEnds"       ; Label bar style
    optsP@pmLabelBarOrthogonalPosF = -0.002                ; Move lb up/down (h/l)
    optsP@pmLabelBarParallelPosF   = 0.04                  ; Move lb right/left (h/l)
    optsP@lbTitleOffsetF        = 0                        ; x-axis title up/down (h/l)
    optsP@txString              = title_arr(it)            ; Plot title string
    optsP@txPosXF               = 0.525                    ; Title pos (x-direction)

    optsP@gsnPanelFigureStrings = (/"a) Relative vorticity at "+hy0+" m"/)
    optsP@gsnPanelFigureStringsFontHeightF = 0.0125

    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.88
    optsP@gsnPanelBottom        = 0.12
    optsP@amJust                = "TopRight"

    if (full .eq. 1) then 
     gsn_panel(wks,panel,(/1,2/),optsP)		; Two panels
    elseif (full .eq. 0) then 
     gsn_panel(wks,panel,(/1,1/),optsP)		; One panel
    end if
    frame(wks)

  ; Tidy up before next iteration
    delete([/vort_plane0,vvel_plane0/])
    delete([/vort_plane1,vvel_plane1/])

   end do     ; End time loop (do it = ts0, tf0))

end