; Script to calculate the tangential momentum budget as in Montgomery et al. (2018)

; Run using:

; ncl dat=\"02T12\" opt=\"x11\" ens0=\"em11\" dist=3.0 ts0=54 tf0=63 clr=5
; full=1 rad0=1 ave=0 wrt=0 new0=0 cn0=\"geo_sm\" res0=0 sc=5 sm=0
; ar=0.5 nr=11 nt=73 lgd0=\"tr\" dnt=1 w0=0.3 plt=1 tc_mbgt_tan_xz.ncl

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'opt'   = output file format ("pdf" or "x11")
; 'ens0'  = ensemble member (em00 -- em11)
; 'dist'  = size of box following storm (degrees)
; 'ts0'   = starting time index (0-119)
; 'tf0'   = ending time index (0-119)
; 'clr'   = colour map for plots ("ncl_default" [1];  "bdr_extra" [2]; "amwg256" [3])
; 'full'  = plot 13 panels inc. u,v,w (2), 10 panels w/mbgt terms only (1) or 4 panels (0)
; 'rad0'  = radius of final plots: 150 km (0) or 200 km (1)
; 'ave'   = calculate 3/6/12-h average (1) rather than analyse individual time interval (0)
; 'wrt'   = write out averages to text file (1) or leave alone (0 – default)
; 'new0'  = reduced (1) or original (0) scale for plotting u,v,w,AAM
; 'cn0'   = method of storm tracking ["slp", "geo", "geo_sm", "vort"]
; 'res'   = plot residual (1) or no residual (0) 
; 'w0'    = vertical velocity threshold for re-calculating storm centre position (0.1, 0.3, etc)

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
load "$sam/ncl_func/st_centre.ncl"
load "$sam/ncl_func/setup_cyl.ncl"
load "$sam/ncl_func/ring_mono.ncl"

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

  date_str  = sprinti("%0.2iUTC ", hour) + \
              sprinti("%0.2i ", day) \
            + month_abbr(month)
  time_str  = sprinti("%0.2i UTC ", hour) + \
              sprinti("%0.2i ", day) \
            + month_abbr(month)
  time_arr  = sprinti("%0.2i", day) + \
              month_abbr(month) + \
              "_" + sprinti("%0.2iZ", hour)

;==================================================================================
;;;;;;;;;;;;;;;;;;;;;;;;;;; Start of main part of script ;;;;;;;;;;;;;;;;;;;;;;;;;;
;==================================================================================

  print_clock("Starting main part of script...")

;========================================
; Define cylindrical coordinate arrays
;========================================

  pii        = 3.14159265
  radii      = fspan(0.0,ar,nr)      ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,nt)   ; Array of azimuth angles (0, 5, 10, etc)
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
  kmh_int0  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Horiz. diffusivity
  vort_int0 = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Relative vorticity
  avo_int0  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Absolute vorticity
  rho_int0  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Density
  v_rad0    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Radial wind
  v_tan0    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Tangential wind

  difv_x0   = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Vert. diffusion (x)
  difv_y0   = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Vert. diffusion (y)
  difv_tot0 = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Vert. diff. (tot)

  difh_x0   = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Horiz. diffusion (x)
  difh_y0   = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Horiz. diffusion (y)

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

  kmh_int0@description = "Horizontal diffusivity"
  kmh_int0@units       = "m~S~2~N~ s~S~-1~N~"
  kmh_int0!0           = "lev"
  kmh_int0!1           = "rad"
  kmh_int0!2           = "azi"
  kmh_int0!3           = "time"

  difv_y0@description   = "Vertical diffusion (y-component)"
  difv_y0@units         = "m s~S~-2~N~"
  difv_y0!0           = "lev"
  difv_y0!1           = "rad"
  difv_y0!2           = "azi"
  difv_y0!3           = "time"

  difv_x0@description   = "Vertical diffusion (x-component)"
  difv_x0@units         = "m s~S~-2~N~"
  difv_x0!0           = "lev"
  difv_x0!1           = "rad"
  difv_x0!2           = "azi"
  difv_x0!3           = "time"

  difh_x0@description = "Horizontal diffusion (x-component)"
  difh_x0@units       = "m s~S~-2~N~"
  difh_x0!0           = "lev"
  difh_x0!1           = "rad"
  difh_x0!2           = "azi"
  difh_x0!3           = "time"

  difh_y0@description = "Horizontal diffusion (y-component)"
  difh_y0@units       = "m s~S~-2~N~"
  difh_y0!0           = "lev"
  difh_y0!1           = "rad"
  difh_y0!2           = "azi"
  difh_y0!3           = "time"

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

; Additional arrays for terms outside main loop (grid spacing; radius)
  du        = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  dy        = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  dr        = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  r1        = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)

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

; Mean absolute angular momentum {aam}
  mean_aam0        = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  mean_aam0!0      = "lev"
  mean_aam0!1      = "rad"
  mean_aam0!2      = "azi"
  mean_aam0!3      = "time"
  mean_aam0@description = "Azimuthally averaged absolute angular momentum"
  mean_aam0@units  = "m s~S~-1~N~"

; Mean horizontal diffusivity {k_mh}
  mean_kmh0        = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  mean_kmh0!0      = "lev"
  mean_kmh0!1      = "rad"
  mean_kmh0!2      = "azi"
  mean_kmh0!3      = "time"
  mean_kmh0@description = "Azimuthally averaged horizontal diffusivity"
  mean_kmh0@units  = "m~S~2~N~ s~S~-1~N~"

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

; Mean density {rho}
  mean_rho0        = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  mean_rho0!0      = "lev"
  mean_rho0!1      = "rad"
  mean_rho0!2      = "azi"
  mean_rho0!3      = "time"
  mean_rho0@description = "Azimuthally averaged density"
  mean_rho0@units  = "kg m~S~-3~N~"

; Mean vertical diffusion (total)
  mean_difv0          = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  mean_difv0!0        = "lev"
  mean_difv0!1        = "rad"
  mean_difv0!2        = "azi"
  mean_difv0!3        = "time"
  mean_difv0@description = "Azimuthally averaged vertical diffusion (tot)"
  mean_difv0@units    = "m s~S~-1~N~"

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

;=========================================================
; Read in required variables from 'pc' and 'pd' streams
;=========================================================

 ; Create array to hold storm track positions using different methods
   centre_new  = new((/7, num_l, numINPUT_d, 2/),float)
   icp_arr     = new((/numINPUT_d/),float)

   c  = addfile(fili_c,"r")         	 ; Read in 'pc' stream [u,v,vort]
   d  = addfile(fili_d,"r")         	 ; Read in 'pd' stream [w,p,t,theta,kmh,tau]

;==================================
; Get the variables we will need
;==================================

   do it = ts0, tf0		; Loop over times in file

    print("it = "+it+" ; ct = "+ct)

  ; 'time_c' array (121) chosen, so 'it' time index used 
    print_clock("Working on time: "+time_str(it)+" (T+"+it+")" )
    title_arr(it) = "Valid at "+time_str(it)+" (T+"+it+")"

;=====================================================
; Read in variables from 'pc' stream (model levels)
;=====================================================

    t0    = lat_arr(cn,it-1,0)
    t1    = lat_arr(cn,it-1,1)
    n0    = lon_arr(cn,it-1,0)
    n1    = lon_arr(cn,it-1,1)
    print_clock("Grid (pc stream) extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

  ; Combine values above into array for input into external function
    ll_arr = (/t0,t1,n0,n1/)

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

;===========================================================
; Calculate rate of strain (following Nguyen et al. 2011)
;=========================================================== 

    dim1   = "latitude"
    dim2   = "longitude"

  ; Calculate horizontal derivatives of 'u' (du_dx, du_dy) 
    du0    = grad_latlon_cfd(u,u&$dim1$,u&$dim2$,False,False)
    du_dy  = du0[0]
    du_dx  = du0[1]
    delete(du0)

  ; Calculate horizontal derivatives of 'v' (dv_dx, dv_dy) 
    dv0    = grad_latlon_cfd(v,v&$dim1$,v&$dim2$,False,False)
    dv_dy  = dv0[0]
    dv_dx  = dv0[1]
    delete(dv0)

  ; Calculate strain rate --> sqrt[ (du_dx - dv_dy)^2 + (dv_dx + du_dy)^2 ] 

  ; Stretching deformation
    str1a   = du_dx - dv_dy
  ; Shearing deformation
    str2a   = dv_dx + du_dy
    str1    = (str1a)^2
    str2    = (str2a)^2

  ; Strain rate (3D array)
    str0    = sqrt( str1(:,:,:) + str2(:,:,:) )

;======================================================
; Read in variables from 'pd' stream (model levels) 
;======================================================

  ; Read in basic variables from 'pd' stream 
    lon1D  = d->longitude_1({n0:n1})                ; longitude [xxx grid points]
    lat1D  = d->latitude_1({t0:t1})                 ; latitude  [xxx grid points]
    lonD   = d->longitude({n0:n1})                  ; longitude [xxx grid points]
    latD   = d->latitude({t0:t1})                   ; latitude  [xxx grid points]
    hybD   = d->hybrid_ht(:)			    ; 63 model ('rho') levels

  ; Read in diffusivity and wind stress
  ; 'str_xz' -- [hybD | 63] * [latD | 150]  * [lonD | 151]
  ; 'str_yz' -- [hybD | 63] * [lat1D | 150] * [lon1D | 151]
  ; 'kmh'    -- [hybD | 63] * [latD | 150]  * [lon1D | 151]
    str_xz = d->taux(it-2,:,{t0:t1},{n0:n1})        ; Wind stress, x-comp (N m/2)
    str_yz = d->tauy(it-2,:,{t0:t1},{n0:n1})        ; Wind stress, y-comp (N m/2)
    kmh	   = d->unspecified(it-2,:,{t0:t1},{n0:n1}) ; Horizontal diffusivity (m^2 s-1)

  ; Read in temperature, pressure, vertical velocity and geopotential height
  ; All variables below -- [hybD | 63] * [latD | 150] * [lon1D | 151]
    t      = d->temp(it-2,:,{t0:t1},{n0:n1})        ; Temperature (K)
    p      = d->p(it-2,:,{t0:t1},{n0:n1})           ; Pressure (Pa)
    w      = d->dz_dt(it-2,:,{t0:t1},{n0:n1})	    ; Vertical velocity (m s-1)
    z      = d->ht(it-2,:,{t0:t1},{n0:n1})    	    ; Geopotential height (m)

  ; Calculate dry air density from pressure and temperature
    rho	   = p / (rd * t)

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

    do ilev = 0, num_l-1 ; Loop over chosen model levels

     hy = sprintf("%0.0f",hybD(ilev))
     print_clock("Working on model level " + ilev+ " ("+hy+" m AGL)")

   ; Calculate components of vertical diffusion term
   ; accounting for top and bottom levels
     ip1 = min((/ilev+1,num_l-1/))
     im1 = max((/ilev-1,0/))

   ; Wind stress components on current model level
     tau_xz = str_xz(ilev,0:d0-1,0:d1-1)
     tau_yz = str_yz(ilev,0:d0-1,0:d1-1)

   ; Wind stress components on model levels above/below current level
     tau_x1 = str_xz(ip1,0:d0-1,0:d1-1)
     tau_x0 = str_xz(im1,0:d0-1,0:d1-1)
     tau_y1 = str_yz(ip1,0:d0-1,0:d1-1)
     tau_y0 = str_yz(im1,0:d0-1,0:d1-1)
     dz	    = hybD(ip1) - hybD(im1)

   ; Vertical derivative of x-component of wind stress
     dtaux_dz = (tau_x1 - tau_x0) / dz

   ; Vertical derivative of y-component of wind stress
     dtauy_dz = (tau_y1 - tau_y0) / dz

   ; X-component of vertical diffusive tendency
     difv_xplane             = rho_a(ilev,:,:) * dtaux_dz
     copy_VarCoords(tau_xz,difv_xplane)
     copy_VarMeta(tau_xz,difv_xplane)
     difv_xplane@units       = "m s~S~-2~N~"
     difv_xplane@title       = "vertical diffusive tendency (x-comp)"
     difv_xplane@long_name   = "vertical diffusive tendency (x-comp)"
     difv_xplane@name        = "difv_x"

   ; Y-component of vertical diffusive tendency
     difv_yplane             = rho_a(ilev,:,:) * dtauy_dz
     copy_VarCoords(tau_yz,difv_yplane)
     copy_VarMeta(tau_yz,difv_yplane)
     difv_yplane@units       = "m s~S~-2~N~"
     difv_yplane@title       = "vertical diffusive tendency (y-comp)"
     difv_yplane@long_name   = "vertical diffusive tendency (y-comp)"
     difv_yplane@name        = "difv_y"

   ; Read in other variables on model levels (u,v,w,f,p,str,rho,vort,avo,kmh)
     u_plane 	            = u(ilev,:,:)
     v_plane		    = v(ilev,:,:)
     w_plane		    = w(ilev,:,:)
     f_plane		    = f0(ilev,:,:)
     prs_plane		    = p(ilev,:,:)
     str_plane              = str0(ilev,:,:)
     rho_plane              = rho(ilev,:,:)
     vort_plane	   	    = vort(ilev,:,:)
     avo_plane		    = avo(ilev,:,:)
     kmh_plane		    = kmh(ilev,:,:)

   ; Calculate vector windspeed and smooth
     spd_plane              = sqrt( (u_plane ^ 2) + (v_plane ^ 2) )
     copy_VarCoords(u_plane, spd_plane)
     vort_pl_smth           = smth9_Wrap(vort_plane, 0.5, 0.5, True)

   ; Combine all arrays above into single, larger array for input into 'setup_cyl'
     pc_plane               = (/u_plane, v_plane, spd_plane, vort_plane, \
                                avo_plane, f_plane, str_plane, str_plane/)
     
   ; NOTE: second 'prs_plane' is a placeholder for PV, which we haven't calculated 
     pd_plane               = (/w_plane, prs_plane, prs_plane, prs_plane, kmh_plane, \
                                difv_xplane, difv_yplane, rho_plane/)

;=======================================================
; Calculate storm centre position (external function)  
;======================================================= 

   ; Set radial distance (º) when looking for speed/vort min (in ringlike phase)
     r0            = 0.12

     print_clock("Calculating storm centre position...")
     centre_arr    = st_centre(ilev, it, cn, r0, centre, vort_plane, prs_plane, spd_plane, \
                               centre_new, dy0, dx0, "single", ct, "pd", "off", "1h")

   ; Updated storm track information contained in 'centre_new' 
   ; centre_new(0,0,it-2,:) --> vort_max 
   ; centre_new(1,0,it-2,:) --> prs_min     
   ; centre_new(2,0,it-2,:) --> prs_min_sm
   ; centre_new(3,0,it-2,:) --> prs_min_sm2   
   ; centre_new(4,0,it-2,:) --> vort_min 
   ; centre_new(5,0,it-2,:) --> spd_min  
   ; centre_new(6,0,it-2,:) --> offline 
     centre_new     = centre_arr[0]

     print("Storm centre = "+centre_new(sc,ilev,it-2,0)+" ; "+centre_new(sc,ilev,it-2,1) )

;=====================================================================
; Make sure that the storm centre is not within the eyewall updraft 
;=====================================================================

   ; Retrieve coordinates of storm centre
     cen_loc0 = centre_new(sc,ilev,it-2,0)
     cen_loc1 = centre_new(sc,ilev,it-2,1)

   ; Coordinates of grid points either side  
     cen_0n = cen_loc0 - dy0
     cen_0x = cen_loc0 + dy0
     cen_1n = cen_loc1 - dy0
     cen_1x = cen_loc1 + dy0

   ; Calculate vertical velocity at these grid points 
     cen_w    = w_plane( {cen_0n:cen_0x}, {cen_1n:cen_1x} )

   ; If average over several grid points is > threshold value, search for centre again
     ave_w    = avg(cen_w)
     ave_out  = sprintf("%0.1f",ave_w)

   ; Coordinates of grid points either side (bigger grid than above)  
     delete([/cen_0n, cen_0x, cen_1n, cen_1x/])
     cen_0n = cen_loc0 - (2 * dy0)
     cen_0x = cen_loc0 + (2 * dy0)
     cen_1n = cen_loc1 - (2 * dy0)
     cen_1x = cen_loc1 + (2 * dy0)

     if (ave_w .gt. w0) then

      print("ave. vertical velocity = "+ave_w)

    ; Create smaller grid around previous centre 
      w_pl_new  = w_plane({cen_0n:cen_0x},{cen_1n:cen_1x})

    ; Find minimum vertical velocity on this smaller grid 
      w_min_new = min(w_pl_new)

    ; Reshape to 1D array, and find index of minimum vertical velocity  
      dims_w    = dimsizes(w_pl_new)
      w_1d      = ndtooned(w_pl_new)
      inds_w    = ind_resolve(minind(w_1d), dims_w)

    ; Retrieve lat/lon information from subset of grid 
      ltN       = w_pl_new&latitude
      lnN       = w_pl_new&longitude_1

      lat_min_w = ltN(0) + (dy0 * inds_w(0,0) )
      lon_min_w = lnN(0) + (dx0 * inds_w(0,1) )
      delete([/ltN, lnN, w_1d, w_pl_new/])

    ; Print new storm centre to screen     
      print("Centre (w_min): "+lat_min_w+" degrees N, "+lon_min_w+" degrees E")

    ; Replace values in 'centre_new' before calling 'setup_cyl' below  
      centre_new(sc,ilev,it-2,0) = lat_min_w
      centre_new(sc,ilev,it-2,1) = lon_min_w

     end if 

   ; Tidy up 
     delete([/cen_w, cen_0n, cen_0x, cen_1n, cen_1x/])

;========================================================
; Switch to cylindrical coordinates (external function) 
;========================================================
 
   ; Call 'setup_cyl' and output pressure/vorticity on cylindrical grid 
   ; Second to last option = 0 --> basic variables output 
   ; Second to last option = 1 --> additional MBGT variables output
   ; Last option --> composite ("comp") or single simulation ("single"); same as 'st_centre' call
   ; Currently, when second to last option = "single", last option is redundant
     print_clock("Translating to cylindrical grid...")
     cyl_arr  = setup_cyl(it, dist, ll_arr, lonC, latC, lonD, latD, lon1D, lat1D, \
                          pc_plane, pd_plane, numTIMES, centre_new, \
                          sc, sm, ar, nr, nt, ilev, num_l, 1, "single", it-2, \
                          0, xpos_all, ypos_all, xcen_all, ycen_all)

     xpos_all = cyl_arr[0]
     ypos_all = cyl_arr[1]
     xcen_all = cyl_arr[2]
     ycen_all = cyl_arr[3]
     lat_max  = cyl_arr[4]
     lon_max  = cyl_arr[5]
     radii    = cyl_arr[6]
     rad_size = cyl_arr[7]

     u_int    = cyl_arr[8]
     v_int    = cyl_arr[9]
     spd_int  = cyl_arr[10]
     vort_int = cyl_arr[11]
     avo_int  = cyl_arr[12]
     f_int    = cyl_arr[13]
     vr_int   = cyl_arr[14]
     vt_int   = cyl_arr[15]
     w_int    = cyl_arr[16]
     prs_int  = cyl_arr[17]
     kmh_int  = cyl_arr[18]
     dvx_int  = cyl_arr[19]
     dvy_int  = cyl_arr[20]
     rho_int  = cyl_arr[21]
     dvt_int  = cyl_arr[22]
     str_int  = cyl_arr[27]

     print("Storm centre = "+lat_max+" ; "+lon_max)

   ; Call external function to determine whether inner core is ringlike or monopolar
     innr_arr = ring_mono(it, radii, vort_int, spd_int, numINPUT_d, icp_arr)
     icp_arr  = innr_arr[0]

   ; If we're in the symmetric phase, put the storm centre in the middle of the vorticity ring
     if (dnt .eq. 1 ) then ; .and. icp_arr(it) .lt. 0.90) then

      print("Re-calculating storm centre to fit ringlike vorticity profile...")

    ; Select the 'vort_min' method of storm tracking 
      sc0 = 4

    ; Call the external function again
      cyl_arr  = setup_cyl(it, dist, ll_arr, lonC, latC, lonD, latD, lon1D, lat1D, \
                           pc_plane, pd_plane, numTIMES, centre_new, \
                           sc0, sm, ar, nr, nt, ilev, num_l, 1, "single", it-2, \
                           0, xpos_all, ypos_all, xcen_all, ycen_all)

      xpos_all = cyl_arr[0]
      ypos_all = cyl_arr[1]
      xcen_all = cyl_arr[2]
      ycen_all = cyl_arr[3]
      lat_max  = cyl_arr[4]
      lon_max  = cyl_arr[5]
      radii    = cyl_arr[6]
      rad_size = cyl_arr[7]

      u_int    = cyl_arr[8]
      v_int    = cyl_arr[9] 
      spd_int  = cyl_arr[10]
      vort_int = cyl_arr[11]
      avo_int  = cyl_arr[12]
      f_int    = cyl_arr[13]
      vr_int   = cyl_arr[14]
      vt_int   = cyl_arr[15]
      w_int    = cyl_arr[16]
      prs_int  = cyl_arr[17]
      kmh_int  = cyl_arr[18]
      dvx_int  = cyl_arr[19]
      dvy_int  = cyl_arr[20]
      rho_int  = cyl_arr[21]
      dvt_int  = cyl_arr[22]
      str_int  = cyl_arr[27]

      print("New storm centre = "+lat_max+" ; "+lon_max)

     end if

;=====================================
; START MID-SCRIPT PLOTTING ROUTINE 
;=====================================

    if (plt .eq. 1) then 

     if (ilev .eq. 20) then 

   ; Plot the data as we go
     output_loop = "$nep/nepartak/images/cyl_coords/xy_4panel_"+dat+"_"+ens0+\
                   "_"+hy+"m_sc"+sc+"_"+time_arr(it)
     wks = gsn_open_wks(opt, output_loop)

     gsn_define_colormap(wks,"prcp_new")

   ; Relative vorticity resources
     opts_rv                               = True
     opts_rv@cnFillOn                      = True
     opts_rv@sfXArray                      = xpos_all(ilev,:,:,it-2)
     opts_rv@sfYArray                      = ypos_all(ilev,:,:,it-2)
     opts_rv@cnLineLabelInterval           = 2.0
     opts_rv@cnLineLabelFontHeightF        = 0.012
     opts_rv@cnLineLabelBackgroundColor    = "transparent"
     opts_rv@cnLineLabelPlacementMode      = "constant"
     opts_rv@cnLinesOn                     = False
     opts_rv@cnInfoLabelOn                 = False
     opts_rv@cnLevelSelectionMode          = "ExplicitLevels"
     opts_rv@cnLevels                      = (/2., 3., 5., 10., 15., 20., 25., 30., \
                                               40., 50., 60., 70., 80./)
     opts_rv@cnFillColors                  = (/0,3,4,5,6,7,8,9,10,11,12,13,14,15/)
     opts_rv@gsnMaximize                   = False

   ; No title or additional strings
     opts_rv@tiMainString                  = ""
     opts_rv@gsnLeftString                 = ""
     opts_rv@gsnRightString                = ""

   ; Turn individual label bar on/off
     opts_rv@lbLabelBarOn                   = True
     opts_rv@lbBoxEndCapStyle               = "TriangleBothEnds"
     opts_rv@lbOrientation                  = "Vertical"

     opts_rv@trYMaxF                        = lat_max + ar
     opts_rv@trYMinF                        = lat_max - ar
     opts_rv@trXMaxF                        = lon_max + ar
     opts_rv@trXMinF                        = lon_max - ar

   ; Do not draw the plot or advance the frame                                        
     opts_rv@gsnDraw                       = False
     opts_rv@gsnFrame                      = False
     opts_rv@gsnAddCyclic                  = False

   ; Background plot resources                                                        
     opts_xy                               = True

   ; No tickmarks on any axis                                                         
     opts_xy@tmYROn                        = False
     opts_xy@tmYLOn                        = True
     opts_xy@tmXTOn                        = False
     opts_xy@tmXBOn                        = True

     opts_xy@xyLineColor                   = -1
     opts_xy@gsnDraw                       = False
     opts_xy@gsnFrame                      = False

     opts_xy@trYMaxF                       = lat_max + ar
     opts_xy@trYMinF                       = lat_max - ar
     opts_xy@trXMaxF                       = lon_max + ar
     opts_xy@trXMinF                       = lon_max - ar
     opts_xy@gsnMaximize                   = False

   ; Resources for windspeed
     opts_wind                             = opts_rv
     delete([/opts_wind@cnLevels, opts_wind@cnFillColors/])
     opts_wind@cnFillOn                    = True
     opts_wind@cnLinesOn                   = False
     opts_wind@cnLineLabelsOn              = False
     opts_wind@cnInfoLabelOn               = False
     opts_wind@cnLevelSelectionMode        = "ExplicitLevels"
     opts_wind@cnFillPalette               = "radar_new"
     opts_wind@cnLevels                    = (/5.0, 10.0, 15.0, 20.0,\
                                               25.0, 30.0, 35.0, 40.0, \
                                               45.0, 50.0, 55.0, 60.0, 65.0/)
     opts_wind@cnFillColors                = (/-1,0,1,2,3,4,5,6,7,8,9,10,11,12/)
     opts_wind@gsnMaximize                 = False

   ; Resources for vertical velocity                                                             
     opts_vvel                             = opts_wind
     delete([/opts_vvel@cnLevels, opts_vvel@cnFillColors, opts_vvel@cnFillPalette/])

     opts_vvel@cnFillPalette               = "BlueDarkRed18"
     opts_vvel@cnLevels                    = (/-3.0, -2.0, -1.5, -1.0, -0.5, -0.2, -0.1, -0.05,\
                                                0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0/)
     opts_vvel@cnFillColors                = (/1,2,3,4,5,6,7,8,-1,9,10,11,12,13,14,15,16/)

   ; Resources for strain rate                                                                   
     opts_strn                             = opts_vvel
     delete([/opts_strn@cnLevels, opts_strn@cnFillColors, opts_strn@cnFillPalette/])
     opts_strn@cnFillPalette               = "circ_budget"
     opts_strn@cnLevels                    = (/10,12,14,16,18,20,22,25,30/)
     opts_strn@cnFillColors                = (/-1,11,13,15,16,17,18,19,20,21/)

   ; Resources for divergence                                                                    
     opts_divg                             = opts_vvel
     delete([/opts_divg@cnLevels, opts_divg@cnFillColors, opts_divg@cnFillPalette/])
     opts_divg@cnFillPalette               = "BlueDarkRed18"
     opts_divg@cnLevels                    = (/-100, -50, -30, -20, -10, -5, \
                                               5, 10, 20, 30, 50, 100/)
     opts_divg@cnFillColors                = (/1,2,3,4,6,8,-1,9,11,13,14,15,16/)

   ; Resources for radial wind                                                                   
     opts_vrad                             = opts_vvel
     delete([/opts_vrad@cnLevels, opts_vrad@cnFillColors, opts_vrad@cnFillPalette/])
     opts_vrad@cnFillPalette               = "BlueDarkRed18"
     opts_vrad@cnLevels                    = (/-20., -15., -10., -5., -3., -2., -1., -0.5, \
                                                0, 0.5, 1., 2., 3., 5., 10., 15., 20./)
     opts_vrad@cnFillColors                = (/0,2,3,4,5,6,7,8,-1,-1,10,11,12,13,14,15,16,17/)

   ; Panel plot label resources                                                                  
     txid                                  = new(4,graphic)
     amid                                  = new(4,graphic)

     panel_strings                         = (/ "a) Relative vorticity",\
                                                "b) Vertical velocity",\
                                                "c) Tangential wind",\
                                                "d) Radial wind"/)

     txres                                  = True
     txres@txPerimOn                        = True
     txres@txFontHeightF                    = 0.0125
     txres@txBackgroundFillColor            = "White"

     amres                                  = True
     amres@amParallelPosF                   = -0.48
     amres@amOrthogonalPosF                 = -0.48
     amres@amJust                           = "TopLeft"

   ; Coordinates of panels                                                              
     xf = (/0.05, 0.57, 0.05, 0.57/)
     yf = (/0.90, 0.90, 0.52, 0.52/)
     wf = 0.35
     hf = 0.32

   ; Panel 1 (relative vorticity)                                                       
     opts_rv@vpXF       = xf(0)
     opts_rv@vpYF       = yf(0)
     opts_rv@vpWidthF   = wf
     opts_rv@vpHeightF  = hf

     opts_xy0           = opts_xy
     opts_xy0@vpXF      = xf(0)
     opts_xy0@vpYF      = yf(0)
     opts_xy0@vpWidthF  = wf
     opts_xy0@vpHeightF = hf

     vort_int_plot      = vort_int * (10 ^ 4)

     plot_tl            = gsn_csm_xy(wks, xpos_all(ilev,:,:,it-2), \
                                     ypos_all(ilev,:,:,it-2), opts_xy0)
     plot_tl0           = gsn_csm_contour(wks, vort_int_plot(:,:), opts_rv)
     overlay(plot_tl, plot_tl0)

   ; Overlay labels 
     txid(0)      = gsn_create_text(wks, panel_strings(0), txres)
     amid(0)      = gsn_add_annotation(plot_tl, txid(0), amres)

   ; Panel 2 (vertical velocity)
     opts_vvel@vpXF       = xf(1)
     opts_vvel@vpYF       = yf(1)
     opts_vvel@vpWidthF   = wf
     opts_vvel@vpHeightF  = hf

     opts_xy1             = opts_xy
     opts_xy1@vpXF        = xf(1)
     opts_xy1@vpYF        = yf(1)
     opts_xy1@vpWidthF    = wf
     opts_xy1@vpHeightF   = hf

     plot_tr              = gsn_csm_xy(wks, xpos_all(ilev,:,:,it-2), \
                                       ypos_all(ilev,:,:,it-2), opts_xy1)
     plot_tr0             = gsn_csm_contour(wks, w_int(:,:), opts_vvel)
     overlay(plot_tr, plot_tr0)

   ; Overlay labels 
     txid(1)      = gsn_create_text(wks, panel_strings(1), txres)
     amid(1)      = gsn_add_annotation(plot_tr, txid(1), amres)

   ; Panel 3 (tangential wind)
     opts_wind@vpXF       = xf(2)
     opts_wind@vpYF       = yf(2)
     opts_wind@vpWidthF   = wf
     opts_wind@vpHeightF  = hf

     opts_xy2             = opts_xy
     opts_xy2@vpXF        = xf(2)
     opts_xy2@vpYF        = yf(2)
     opts_xy2@vpWidthF    = wf
     opts_xy2@vpHeightF   = hf

     plot_bl              = gsn_csm_xy(wks, xpos_all(ilev,:,:,it-2), \
                                       ypos_all(ilev,:,:,it-2), opts_xy2)
     plot_bl0             = gsn_csm_contour(wks, vt_int(:,:), opts_wind)
     overlay(plot_bl, plot_bl0)

   ; Overlay labels 
     txid(2)      = gsn_create_text(wks, panel_strings(2), txres)
     amid(2)      = gsn_add_annotation(plot_bl, txid(2), amres)

   ; Panel 4 (radial wind)                                                              
     opts_vrad0           = opts_vrad
     opts_vrad0@vpXF      = xf(3)
     opts_vrad0@vpYF      = yf(3)
     opts_vrad0@vpWidthF  = wf
     opts_vrad0@vpHeightF = hf

     opts_xy3             = opts_xy
     opts_xy3@vpXF        = xf(3)
     opts_xy3@vpYF        = yf(3)
     opts_xy3@vpWidthF    = wf
     opts_xy3@vpHeightF   = hf

     plot_br              = gsn_csm_xy(wks, xpos_all(ilev,:,:,it-2), \
                                       ypos_all(ilev,:,:,it-2), opts_xy3)
     plot_br0             = gsn_csm_contour(wks, vr_int(:,:), opts_vrad0)
     overlay(plot_br, plot_br0)

   ; Overlay labels                                                                     
     txid(3)      = gsn_create_text(wks, panel_strings(3), txres)
     amid(3)      = gsn_add_annotation(plot_br, txid(3), amres)

   ; Tidy up before final plots later in script 
     delete([/opts_vvel, opts_wind, opts_rv, opts_strn, opts_divg, opts_vrad, \
              opts_xy, opts_xy0, opts_xy1, opts_xy2, opts_xy3,\
              txid, amid, panel_strings/])

;==================================================
; Overlay markers for each storm centre estimate
;================================================== 

; (0) vort, (1) pressure, (2) pressure_smth, (3), pressure_smth2, (4) vort_min,
; (5) speed_min, (6) offline 
  m_colours  = (/"royalblue1", "orange4", "orange4", "orange4", "royalblue1", \
                 "black", "red4"/)
  m_indices  = (/12, 12, 16, 8, 8, 8, 12/)
  m_sizes    = (/10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0/)
  m_thick    = (/5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0/)

  markers    = new( (/4, dimsizes(m_colours)/), graphic)

  do i = 0, dimsizes(m_colours)-1
    mres0                    = True
    mres0@gsMarkerIndex      = m_indices(i)
    mres0@gsMarkerSizeF      = m_sizes(i)
    mres0@gsMarkerColor      = m_colours(i)
    mres0@gsMarkerThicknessF = m_thick(i)
    markers(0,i)             = gsn_add_polymarker(wks, plot_tl, centre_new(i,ilev,it-2,1),\
                                                  centre_new(i,ilev,it-2,0), mres0)
    markers(1,i)             = gsn_add_polymarker(wks, plot_tr, centre_new(i,ilev,it-2,1),\
                                                  centre_new(i,ilev,it-2,0), mres0)
    markers(2,i)             = gsn_add_polymarker(wks, plot_bl, centre_new(i,ilev,it-2,1),\
                                                  centre_new(i,ilev,it-2,0), mres0)
    markers(3,i)             = gsn_add_polymarker(wks, plot_br, centre_new(i,ilev,it-2,1),\
                                                  centre_new(i,ilev,it-2,0), mres0)
  end do 

;===================================
; Add a legend (only on one panel)
;===================================

; Display legend; give us control; label font; label height 
  lg_opts                            = True
  lg_opts@pmLegendDisplayMode        = "Always"
  lg_opts@lgAutoManage               = False
  lg_opts@lgLabelFont                = "Helvetica"
  lg_opts@lgLabelFontHeightF         = 0.060

; Customise labels and line colours 
  lg_opts@lgItemType                 = "Markers"
  lg_opts@lgMarkerColors             = m_colours
  lg_opts@lgMarkerIndexes            = m_indices
  lg_opts@lgMarkerSizeF              = mres0@gsMarkerSizeF
  lg_opts@lgMarkerThicknessF         = mres0@gsMarkerThicknessF
  lsize                              = dimsizes(lg_opts@lgMarkerColors)

  lg_opts@lgLabelPosition            = "Right"
  lg_opts@lgItemPlacement            = "ExplicitPlacement"
  lg_opts@lgItemPositions            = fspan(0.05, 0.95, lsize)
  lg_opts@lgLeftMarginF              = 0.01
  lg_opts@lgRightMarginF             = 0.15

  lg_opts@vpWidthF                   = 0.150
  lg_opts@vpHeightF                  = 0.090
  lg_opts@lgPerimColor               = "black"
  lg_opts@lgPerimThicknessF          = 3.0
  lg_opts@lgPerimFill                = "SolidFill"
  lg_opts@lgPerimFillColor           = "white"

  lab_arr                            = (/"vort", "pres", \
                                         "pres_sm", "pres_sm2", \
                                         "vort_min", "windspeed_min", \
                                         "offline ("+cn0+")"/)
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

  annotate_0                         = gsn_add_annotation(plot_bl, lgd, am_opts)

;==============================
; Add radial circles to plot
;==============================

      dtr = 0.017453292519943

    ; Get customizations for circles
      radius   = radii(rad_size)
      xcenter  = lon_max
      ycenter  = lat_max
      out_thck = 2
      in_thck  = 1
      spacing  = 0.05
      dpattrn  = 2
      out_clr  = 1
      in_clr   = 1

;============================================= 
; Construct and attach outer circle to plot
;============================================= 

    ; Calculate arrays for outer circle
      degrees  = ispan(0, 360, 5)
      xcos     = cos(dtr * degrees)
      xsin     = sin(dtr * degrees)

    ; Array of points (x,y) representing the outer circle
      xc       = xcenter + (radius * xcos)
      yc       = ycenter + (radius * xsin)

    ; Resources for outer circle
      lnres    = True
      lnres@gsLineThicknessF  = out_thck
      lnres@gsLineColor       = out_clr
      lnres@gsLineDashPattern = 0

    ; Attach outer circle to each panel plot
      plot_outer0 = gsn_add_polyline(wks, plot_tl, xc, yc, lnres)
      plot_outer1 = gsn_add_polyline(wks, plot_tr, xc, yc, lnres)
      plot_outer2 = gsn_add_polyline(wks, plot_bl, xc, yc, lnres)
      plot_outer3 = gsn_add_polyline(wks, plot_br, xc, yc, lnres)

      if (spacing .gt. radius) then
       print("add_radial_circles: spacing is > radius, can't draw inner circles.")
       return
      end if

;================================================== 
; Now construct and attach inner circles to plot 
;================================================== 

    ; Draw inner circles if desired
      if (spacing .gt. 0) then

       count = 0
       size  = toint( (radius-spacing) / spacing )
       plot_inner = new( (/4,size/), graphic)

       do r = spacing, radius, spacing

       ; Calculate arrays for inner circles
         xc = xcenter + (r * xcos)
         yc = ycenter + (r * xsin)

       ; Resources for inner circles
         delete([/lnres@gsLineThicknessF,lnres@gsLineColor,lnres@gsLineDashPattern/])
         lnres@gsLineColor       = in_clr
         lnres@gsLineThicknessF  = in_thck
         lnres@gsLineDashPattern = dpattrn

       ; Attach inner circles to plot
         plot_inner(0,count) = gsn_add_polyline(wks, plot_tl, xc, yc, lnres)
         plot_inner(1,count) = gsn_add_polyline(wks, plot_tr, xc, yc, lnres)
         plot_inner(2,count) = gsn_add_polyline(wks, plot_bl, xc, yc, lnres)
         plot_inner(3,count) = gsn_add_polyline(wks, plot_br, xc, yc, lnres)

         count = count + 1

       end do

      end if

    ; Tidy up
      delete([/lnres, degrees, xcos, xsin, xc, yc/])

    ; Create a main title at the top of the page
      title_res                   = True
      title_str                   = "Panel plot at "+hy+" m AGL: "+ens0+", "+dat+" (T+"+it+")"
      title_res@txFontHeightF     = 0.02
      title_res@txFont            = "Helvetica"
      xpos_title                  = 0.50
      ypos_title                  = 0.94

      gsn_text_ndc(wks, title_str, xpos_title, ypos_title, title_res)

    ; Call 'draw' and 'frame'
      maximize_output(wks, True)

    ; Only loop over selected model levels
      end if 

  ; END MID-SCRIPT PLOTTING ROUTINE 
    end if   

   ; Fill 'u_dz0' array with values from 'difv_tot' (vertical diffusive tendency)
     v_dz0(ilev,:,:,ct) = dvt_int(:,:)

;===================================================
; Calculate mean (azimuthally-averaged) variables 
;===================================================

     do irad = 0, dimsizes(radii)-1
      mean_vtan0(ilev,irad,:,ct)    = avg(vt_int(irad,:))      ; [ {v} ]
      mean_vrad0(ilev,irad,:,ct)    = avg(vr_int(irad,:))      ; [ {u} ]
      mean_vvel0(ilev,irad,:,ct)    = avg(w_int(irad,:))      ; [ {w} ]
      mean_vort0(ilev,irad,:,ct)    = avg(vort_int(irad,:))   ; [ {zeta} ]
      mean_avo0(ilev,irad,:,ct)     = avg(avo_int(irad,:))    ; [ {zeta + f} ]
      mean_difv0(ilev,irad,:,ct)    = avg(dvt_int(irad,:))    ; [ {diff_xz + diff_yz} ]
      mean_rho0(ilev,irad,:,ct)     = avg(rho_int(irad,:))    ; [ {rho} ]
      mean_kmh0(ilev,irad,:,ct)     = avg(kmh_int(irad,:))    ; [ {k_mh} ]
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
     kmh_int0(ilev,:,:,ct)  = (/kmh_int(:,:)/)
     rho_int0(ilev,:,:,ct)  = (/rho_int(:,:)/)
     difv_x0(ilev,:,:,ct)   = (/dvx_int(:,:)/)
     difv_y0(ilev,:,:,ct)   = (/dvy_int(:,:)/)
     difv_tot0(ilev,:,:,ct) = (/dvt_int(:,:)/)
     v_rad0(ilev,:,:,ct)    = (/vr_int(:,:)/)
     v_tan0(ilev,:,:,ct)    = (/vt_int(:,:)/)

   ; Tidy up
     delete([/u_int,v_int,w_int,vt_int,vr_int/])
     delete([/rho_int,dvx_int,dvy_int,f_int,kmh_int,vort_int,avo_int/])
     delete([/u_plane,v_plane,w_plane,f_plane,vort_plane,avo_plane/])
     delete([/difv_xplane,difv_yplane,kmh_plane/])

    end do     	            ; End model level loop (do ilev = 0, num_l-1)

   ct  = ct + 1		    ; Counter variable (time)

 ; Tidy up before next iteration (all variables from 'pc' and 'pd' streams)
   delete([/lonC,latC,lon1C,lat1C,u,v,vort,avo,f0,f1/])
   delete([/lonD,latD,lon1D,lat1D,str_xz,str_yz,kmh,t,p,w,z,rho,rho_a/])

 end do     ; End time loop (do it = 0, times-1)

;====================================
; Calculate perturbation variables
;====================================

      ptb_vtan0(:,:,:,:) = v_tan0(:,:,:,:)    - mean_vtan0(:,:,:,:)    ; [v']
      ptb_vvel0(:,:,:,:) = w_int0(:,:,:,:)    - mean_vvel0(:,:,:,:)    ; [w']
      ptb_vrad0(:,:,:,:) = v_rad0(:,:,:,:)    - mean_vrad0(:,:,:,:)    ; [u']
      ptb_vort0(:,:,:,:) = vort_int0(:,:,:,:) - mean_vort0(:,:,:,:)    ; [zeta']

      k_mh0(:,:,:,:)     = kmh_int0(:,:,:,:)                           ; [kmh]

;===========================================================
; Calculate term involving azimuthal derivatives [du/dy]
;===========================================================

    ; Calculate the azimuthal distance ('dy') between all points in array
    ; Arguments 1 and 2 represent the lat/lon of the first point, respectively
    ; Arguments 3 and 4 represent the lat/lon of the second point, respectively
    ; Argument 5 ('2') tells the script to avoid extra interpolated points
    ; Argument 6 ('3') returns the distance in metres

      dth  = dimsizes(thetas)

      dy(:,:,0:dth-2,:)  = gc_latlon( ypos_all(:,:,0:dth-2,:),xpos_all(:,:,0:dth-2,:),\
                                      ypos_all(:,:,1:dth-1,:),xpos_all(:,:,1:dth-1,:), 2, 3)

    ; Replace zeros with positive values
      dy  = where( dy .ne. 0, dy, 2000)

    ; Calculate radial wind increment in azimuthal direction
      du(:,:,0:dth-2,:)  = v_rad0(:,:,1:dth-1,:) - v_rad0(:,:,0:dth-2,:)

    ; Calculate [du/dy -- for horizontal diffusive tendency]
      dudy0(:,:,:,:) = du / dy

;===========================================================
; Calculate term involving radial derivatives [d(v/r)/dr]
;===========================================================

    ; Calculate the radial distance ('dr') between all points in array
    ; Arguments 1 and 2 represent the lat/lon of the first point, respectively
    ; Arguments 3 and 4 represent the lat/lon of the second point, respectively
    ; Argument 5 ('2') tells the script to avoid extra interpolated points
    ; Argument 6 ('3') returns the distance in metres

      drad = dimsizes(radii)

      dr(:,0:drad-2,:,:)  = gc_latlon( ypos_all(:,0:drad-2,:,:),xpos_all(:,0:drad-2,:,:),\
                                       ypos_all(:,1:drad-1,:,:),xpos_all(:,1:drad-1,:,:), 2, 3)

    ; Replace zeros with positive values
      dr = where( dr .ne. 0, dr, 2000)

      r1(:,0:drad-2,:,:)  = gc_latlon( ycen_all(:,1:drad-1,:,:),xcen_all(:,1:drad-1,:,:),\
                                       ypos_all(:,1:drad-1,:,:),xpos_all(:,1:drad-1,:,:), 2, 3)

    ; Replace zeros with positive values
      r1 = where( r1 .ne. 0, r1, 2000)

    ; Calculate components of horizontal diffusive tendency source term

    ; Calculate [v/r]
      v_r0(:,:,:,:)               = v_tan0(:,:,:,:) / r1

    ; Calculate [d(v/r)/dr]
      dvdr0(:,:,:,:)              = v_r0(:,:,:,:) / dr

    ; Calculate [r]
      r_0(:,:,:,:)                = r1

    ; Calculate [d/dr]
      ddr0(:,:,:,:)               = 1 / dr

;=============================================
; Calculate absolute angular momentum (AAM)
;=============================================

;      printMinMax(mean_vtan0(:,:,:,:), False)
;      printMinMax(mean_cor0(:,:,:,:), False)
;      printMinMax(r1, False)
;      exit()

    ; Calculate increments
      vm   = mean_vtan0(:,:,:,:)
      fm   = mean_cor0(:,:,:,:) * (r1 ^ 2)

    ; Calculate absolute angular momentum (AAM)
      aam0 = vm * r1
      aam1 = 0.5 * fm

    ; Calculate absolute angular momentum [ rv + 0.5 * (fr ^ 2) ]
      mean_aam0(:,:,:,:)  = aam0 + aam1

;      printMinMax(aam0, False)
;      printMinMax(aam1, False)

;=======================================================
; Remove unncessary 3rd dimension from mean variables
;=======================================================

      mean_vtan1 = mean_vtan0(:,:,0,:)
      mean_vrad1 = mean_vrad0(:,:,0,:)
      mean_vvel1 = mean_vvel0(:,:,0,:)
      mean_vort1 = mean_vort0(:,:,0,:)
      mean_avo1  = mean_avo0(:,:,0,:)
      mean_rho1  = mean_rho0(:,:,0,:)
      mean_aam1  = mean_aam0(:,:,0,:)

;================================================
; Calculate vertical derivatives of mean terms 
;================================================

    ; Redefine model level array for vertical derivative calculations
      hybC   = c->hybrid_ht(:)         ; 63 model ('theta') levels
      mlevs  = dimsizes(hybC)	       ; Size of model level array
      hybC@units = "m"		       ; Model level units 

      do it = ts0, tf0

	do ilev = 0, num_l-1		; Loop over all model levels
	
        mp1 = min((/ilev+1,num_l-1/))
        mm1 = max((/ilev-1,0/))
 
          do irad = 0, dimsizes(radii)-1

          ;=========================================================
	  ; Calculate vertical derivative of mean tangential wind
	  ;=========================================================

          ; [ d{v}/dz ]
            mean_dvdz0(ilev,irad,it)      = (mean_vtan1(mp1,irad,it) - \
                                             mean_vtan1(mm1,irad,it)) / \
                                            (hybC(mp1) - hybC(mm1))

          end do                 ; End loop over radial circles (do irad = 0, ...)

        end do 		         ; End loop over model levels (do ilev = 0, num_l-1)

       end do			 ; End loop over times (do it = 0, numTIMES-1)

;================================================
; Calculate vertical derivatives of eddy terms 
;================================================

      do it = ts0, tf0

        do ilev = 0, num_l-1

        mp1 = min((/ilev+1,num_l-1/))
        mm1 = max((/ilev-1,0/))

          do irad = 0, dimsizes(radii)-1

           do th0 = 0, dimsizes(thetas)-1

           ;=========================================================
           ; Calculate vertical derivative of eddy tangential wind
           ;=========================================================

           ; [ dv'/dz ]
             ptb_dvdz0(ilev,irad,th0,it) = (ptb_vtan0(mp1,irad,th0,it) - \
                                            ptb_vtan0(mm1,irad,th0,it)) / \
                                           (hybC(mp1) - hybC(mm1))

           end do	         ; End loop over azimuth angles (do th0 = 0, ...)

          end do                 ; End loop over radial circles (do irad = 0, ...)

        end do                   ; End loop over model levels (do ilev = 0, num_l-1)

       end do                    ; End loop over times (do it = ts0, tf0)

;==========================================================
; Finally, calculate source terms and azimuthal averages 
;==========================================================

       print_clock("Calculating source terms and azimuthal averages...")

     ;=============================================================
     ; Mean radial influx of absolute vertical vorticity: RHS (1)
     ;=============================================================

       rhs1_st = -( mean_vrad1 * mean_avo1 )

     ;===============================================================
     ; Mean vertical advection of mean tangential momentum: RHS (2)
     ;===============================================================

       rhs2_st = -(mean_vvel1 * mean_dvdz0)
     
     ;=======================================
     ; Eddy radial vorticity flux: RHS (3) 
     ;=======================================

       ptb_vez1 = -(ptb_vrad0 * ptb_vort0)

     ; Calculate the azimuthal average of 'ptb_vez1' and input into 'rhs3_st'
       rhs3_st	= dim_avg_n(ptb_vez1(:,:,:,:), 2)

     ;===========================================================
     ; Vertical advection of eddy tangential momentum: RHS (4)
     ;===========================================================

       ptb_vev1 = -(ptb_vvel0 * ptb_dvdz0)

     ; Calculate the azimuthal average of 'ptb_vev1' and input into 'rhs4_st'
       rhs4_st  = dim_avg_n(ptb_vev1(:,:,:,:), 2)

     ;===================================================================================
     ; Combined mean tangential/vert diffusive tendency of tangential momentum: RHS (5)
     ;===================================================================================

     ; Stress tensor
       tau_ry1	 = dim_avg_n( ( k_mh0 * ( ( (1 / r_0) * dudy0) + dvdr0 ) ), 2)

     ;=============================================================================
     ; Use stress tensor to calculate horizontal diffusive tendency source term
     ;=============================================================================

     ; Azimuthally average the 'r' and 'd/dr' expressions
       r_1   = dim_avg_n(r_0(:,:,:,:),2)
       ddr1  = dim_avg_n(ddr0(:,:,:,:),2)

       v_dh0 = ( ( 1/( r_1 ^ 2 * mean_rho1) ) * \
               (ddr1 * ( r_1 ^ 2 * mean_rho1 * tau_ry1) ) )

       rhs5a_st = v_dh0	       

     ; Vertical diffusive tendency source term has already been calculated (L827)
     ; Now calculate the azimuthal average using 'dim_avg_n'
       rhs5b_st = dim_avg_n( v_dz0(:,:,:,:),2 )

;=================================================
; Tidy up and produce final arrays for plotting 
;=================================================

     ; Calculating the terms from Montgomery et al. (2018) - their Eq.(4)
     ; Ignore the perturbation PGF term because it is several orders of magnitude smaller

     ; LHS (1) : d{v}/dt          - Mean tangential momentum tendency

     ; RHS (1) : -{u} * {f+vort}  - Mean radial influx of absolute vertical vorticity
     ; RHS (2) : -{w} * (d{v}/dz) - Mean vertical advection of mean tangential momentum
     ; RHS (3) : -{u' * vort'}    - Eddy radial vorticity flux
     ; RHS (4) : -{w' * (dv'/dz)} - Vertical advection of eddy tangential momentum
     ; RHS (5) : {D_v}            - Combined diffusive and PBL tendency

     ; Define new arrays for terms in momentum budget 
       lhs1A  = new((/num_l,dimsizes(radii),numTIMES/),float)   ; d{v}/dt

       rhs1A  = new((/num_l,dimsizes(radii),numTIMES/),float)   ; -{u} * {f+vort}
       rhs2A  = new((/num_l,dimsizes(radii),numTIMES/),float)   ; -{w} * (d{v}/dz)
       rhs3A  = new((/num_l,dimsizes(radii),numTIMES/),float)   ; -{u' * vort'}
       rhs4A  = new((/num_l,dimsizes(radii),numTIMES/),float)   ; -{w' * (dv'/dz)}
       rhs5aA = new((/num_l,dimsizes(radii),numTIMES/),float)   ; {D_r}
       rhs5bA = new((/num_l,dimsizes(radii),numTIMES/),float)   ; {D_z}
       rhsTA  = new((/num_l,dimsizes(radii),numTIMES/),float)   ; Sum of RHS terms

     ; Define additional arrays if calculating averages

       if (ave .eq. 1) then

         lhs1_av0  = new((/num_l,dimsizes(radii)/),float)

         rhs1_av0  = new((/num_l,dimsizes(radii)/),float)
         rhs2_av0  = new((/num_l,dimsizes(radii)/),float)
         rhs3_av0  = new((/num_l,dimsizes(radii)/),float)
         rhs4_av0  = new((/num_l,dimsizes(radii)/),float)
         rhs5a_av0 = new((/num_l,dimsizes(radii)/),float)
         rhs5b_av0 = new((/num_l,dimsizes(radii)/),float)
         rhsT_av0  = new((/num_l,dimsizes(radii)/),float)

       end if

     ;=================================================
     ; Calculate mean tangential wind tendency (LHS)
     ;=================================================

       if (ave .eq. 1) then 

         dt0 = tf0 - ts0 
	
 	 lhs1_av0(:,:) = (mean_vtan1(:,:,tf0) - mean_vtan1(:,:,ts0)) / dt0

       elseif (ave .eq. 0) then 

         do it = ts0, tf0
 
           cp1 = min((/it+1,tf0/))
  	   cm1 = max((/it-1,ts0/))
 
           lhs1A(:,:,it) = (mean_vtan1(:,:,cp1) - mean_vtan1(:,:,cm1)) / 2

         end do		 ; End loop over times (do it = ts0, tf0)

       end if 		 ; End IF statement 

     ; Add metadata and make units consistent for all terms
       lhs1A@units = "m s~S~-1~N~ h~S~-1~N~"

       rhs1A = (/rhs1_st/)
       rhs1A = rhs1A * 3600
       rhs1A@units = "m s~S~-1~N~ h~S~-1~N~"

       rhs2A = (/rhs2_st/)
       rhs2A = rhs2A * 3600
       rhs2A@units = "m s~S~-1~N~ h~S~-1~N~"
       
       rhs3A = (/rhs3_st/)
       rhs3A = rhs3A * 3600
       rhs3A@units = "m s~S~-1~N~ h~S~-1~N~"

       rhs4A = (/rhs4_st/)
       rhs4A = rhs4A * 3600
       rhs4A@units = "m s~S~-1~N~ h~S~-1~N~"

       rhs5aA = (/rhs5a_st/)
       rhs5aA = rhs5aA * 3600 
       rhs5aA@units = "m s~S~-1~N~ h~S~-1~N~"

       rhs5bA = (/rhs5b_st/)
       rhs5bA = rhs5bA * 3600
       rhs5bA@units = "m s~S~-1~N~ h~S~-1~N~"

       hyb_plot   = hybC / 1000
       hyb_plot@units = "km"

     ; Add metadata for the time-averaged tangential wind tendency
       if (ave .eq. 1) then
        lhs1_av0!0    = "lev"
        lhs1_av0&lev  = hyb_plot
        lhs1_av0!1    = "rad"
        lhs1_av0&rad  = radii
       end if

     ; Add metadata for all the other arrays (instantaneous)
       lhs1A!0    = "lev"
       lhs1A&lev  = hyb_plot
       lhs1A!1    = "rad"
       lhs1A&rad  = radii
       lhs1A!2    = "time"

       rhs1A!0    = "lev"
       rhs1A&lev  = hyb_plot
       rhs1A!1    = "rad"
       rhs1A&rad  = radii
       rhs1A!2    = "time"

       rhs2A!0    = "lev"
       rhs2A&lev  = hyb_plot
       rhs2A!1    = "rad"
       rhs2A&rad  = radii
       rhs2A!2    = "time"

       rhs3A!0    = "lev"
       rhs3A&lev  = hyb_plot
       rhs3A!1    = "rad"
       rhs3A&rad  = radii
       rhs3A!2    = "time"

       rhs4A!0    = "lev"
       rhs4A&lev  = hyb_plot
       rhs4A!1    = "rad"
       rhs4A&rad  = radii
       rhs4A!2    = "time"

       rhs5aA!0   = "lev"
       rhs5aA&lev = hyb_plot
       rhs5aA!1   = "rad"
       rhs5aA&rad = radii
       rhs5aA!2   = "time"

       rhs5bA!0   = "lev"
       rhs5bA&lev = hyb_plot
       rhs5bA!1   = "rad"
       rhs5bA&rad = radii
       rhs5bA!2   = "time"

       rhsTA!0    = "lev"
       rhsTA&lev  = hyb_plot
       rhsTA!1    = "rad"
       rhsTA&rad  = radii
       rhsTA!2    = "time"

     ;==============================================
     ; Calculate time averages of all terms above
     ;==============================================

       if (ave .eq. 1) then

        rhs1_av0   = dim_avg_n_Wrap(rhs1A(:,:,ts0:tf0),2)
        rhs2_av0   = dim_avg_n_Wrap(rhs2A(:,:,ts0:tf0),2)
        rhs3_av0   = dim_avg_n_Wrap(rhs3A(:,:,ts0:tf0),2)
        rhs4_av0   = dim_avg_n_Wrap(rhs4A(:,:,ts0:tf0),2)
        rhs5a_av0  = dim_avg_n_Wrap(rhs5aA(:,:,ts0:tf0),2)
        rhs5b_av0  = dim_avg_n_Wrap(rhs5bA(:,:,ts0:tf0),2)
        rhsT_av0   = dim_avg_n_Wrap(rhsTA(:,:,ts0:tf0),2)

       end if

;==============================================
; Switch order of dimensions before plotting
;==============================================

     ; 'gsn_csm_contour' plots the leftmost (rightmost) dimension on the y-axis (x-axis)
     ; Switch 'rad' with 'time' so that 'rad' is the rightmost dimension

     ; Mean tangential wind -- also add metadata
       mean_vtan       = mean_vtan1(lev|:,time|:,rad|:)
       mean_vtan!0     = "lev"
       mean_vtan&lev   = hyb_plot
       mean_vtan!1     = "time"
       mean_vtan!2     = "rad"
       mean_vtan&rad   = radii
       mean_vtan@units = "m s~S~-1~N~"

     ; Mean radial wind -- also add metadata
       mean_vrad       = mean_vrad1(lev|:,time|:,rad|:)
       mean_vrad!0     = "lev"
       mean_vrad&lev   = hyb_plot
       mean_vrad!1     = "time"
       mean_vrad!2     = "rad"
       mean_vrad&rad   = radii
       mean_vrad@units = "m s~S~-1~N~"

     ; Mean vertical velocity -- also add metadata
       mean_vvel       = mean_vvel1(lev|:,time|:,rad|:)
       mean_vvel!0     = "lev"
       mean_vvel&lev   = hyb_plot
       mean_vvel!1     = "time"
       mean_vvel!2     = "rad"
       mean_vvel&rad   = radii
       mean_vvel@units = "m s~S~-1~N~"

     ; Mean absolute angular momentum
       mean_aam        = mean_aam1(lev|:,time|:,rad|:)
       mean_aam!0      = "lev"
       mean_aam&lev    = hyb_plot
       mean_aam!1      = "time"
       mean_aam!2      = "rad"
       mean_aam&rad    = radii
       mean_aam        = mean_aam / (10 ^ 6)
       mean_aam@units  = "10~S~6~N~ m~S~2~N~ s~S~-1~N~"

     ;===================================================
     ; Also calculate time averages of mean variables
     ;===================================================

       if (ave .eq. 1) then

        vrad_av0   = dim_avg_n_Wrap(mean_vrad(:,ts0:tf0,:),1)     ; Radial wind 
        vtan_av0   = dim_avg_n_Wrap(mean_vtan(:,ts0:tf0,:),1)     ; Tangential wind
        vvel_av0   = dim_avg_n_Wrap(mean_vvel(:,ts0:tf0,:),1)     ; Vertical velocity
	aam_av0    = dim_avg_n_Wrap(mean_aam(:,ts0:tf0,:),1)      ; Absolute angular momentum

       end if

     ;==================================
     ; Reorder arrays before plotting
     ;==================================

       lhs1  	  = lhs1A(lev|:,time|:,rad|:)

       rhs1  	  = rhs1A(lev|:,time|:,rad|:)
       rhs2  	  = rhs2A(lev|:,time|:,rad|:)
       rhs3  	  = rhs3A(lev|:,time|:,rad|:)
       rhs4  	  = rhs4A(lev|:,time|:,rad|:)
       rhs5a 	  = rhs5aA(lev|:,time|:,rad|:)
       rhs5b 	  = rhs5bA(lev|:,time|:,rad|:)

       if (ave .eq. 1) then

         lhs1_av  = lhs1_av0(lev|:,rad|:)
         rhs1_av  = rhs1_av0(lev|:,rad|:)
         rhs2_av  = rhs2_av0(lev|:,rad|:)
         rhs3_av  = rhs3_av0(lev|:,rad|:)
         rhs4_av  = rhs4_av0(lev|:,rad|:)
         rhs5a_av = rhs5a_av0(lev|:,rad|:)
         rhs5b_av = rhs5b_av0(lev|:,rad|:)

       end if

;=========================================
; Finally, smooth data before plotting
;=========================================

       lhs1	  = smth9_Wrap(lhs1, 0.5, 0.5, True)
       rhs1  	  = smth9_Wrap(rhs1, 0.5, 0.5, True)
       rhs2  	  = smth9_Wrap(rhs2, 0.5, 0.5, True)
       rhs3  	  = smth9_Wrap(rhs3, 0.5, 0.5, True)
       rhs4  	  = smth9_Wrap(rhs4, 0.5, 0.5, True)
       rhs5a 	  = smth9_Wrap(rhs5a, 0.5, 0.5, True)
       rhs5b 	  = smth9_Wrap(rhs5b, 0.5, 0.5, True)

     ; Smooth again
       lhs1       = smth9_Wrap(lhs1, 0.5, 0.5, True)
       rhs1       = smth9_Wrap(rhs1, 0.5, 0.5, True)
       rhs2       = smth9_Wrap(rhs2, 0.5, 0.5, True)
       rhs3       = smth9_Wrap(rhs3, 0.5, 0.5, True)
       rhs4       = smth9_Wrap(rhs4, 0.5, 0.5, True)
       rhs5a      = smth9_Wrap(rhs5a, 0.5, 0.5, True)
       rhs5b      = smth9_Wrap(rhs5b, 0.5, 0.5, True)

     ; Calculate total RHS tendency
       rhsT       = rhs1 + rhs2 + rhs3 + rhs4 + rhs5a + rhs5b
       rhsT!0     = "lev"
       rhsT&lev   = hyb_plot
       rhsT!1     = "time"
       rhsT!2     = "rad"
       rhsT&rad   = radii
       rhsT@units = "m s~S~-1~N~ h~S~-1~N~"

     ; Calculate residual (difference between LHS and RHS)
       res        = rhsT - lhs1
       res!0      = "lev"
       res&lev    = hyb_plot
       res!1      = "time"
       res!2      = "rad"
       res&rad    = radii
       res@units  = "m s~S~-1~N~ h~S~-1~N~"

     ; Also calculate total diffusive tendency                                          
       rhs5       = rhs5a + rhs5b
       rhs5!0     = "lev"
       rhs5&lev   = hyb_plot
       rhs5!1     = "time"
       rhs5!2     = "rad"
       rhs5&rad   = radii
       rhs5@units = "m s~S~-1~N~ h~S~-1~N~"

     ; Calculate sum of all mean terms (vadv + hadv + diff)
       mean       = rhs1 + rhs2 + rhs5a + rhs5b
       mean!0     = "lev"
       mean&lev   = hyb_plot
       mean!1     = "time"
       mean!2     = "rad"
       mean&rad   = radii
       mean@units = "m s~S~-1~N~ h~S~-1~N~"

     ; Calculate sum of all eddy terms (vadv + hadv)
       eddy       = rhs3 + rhs4
       eddy!0     = "lev"
       eddy&lev   = hyb_plot
       eddy!1     = "time"
       eddy!2     = "rad"
       eddy&rad   = radii
       eddy@units = "m s~S~-1~N~ h~S~-1~N~"

     ; Smooth one final time 
       mean       = smth9_Wrap(mean, 0.5, 0.5, True)
       eddy       = smth9_Wrap(eddy, 0.5, 0.5, True)
       rhsT       = smth9_Wrap(rhsT, 0.5, 0.5, True)
       rhs5       = smth9_Wrap(rhs5, 0.5, 0.5, True)
       res        = smth9_Wrap(res, 0.5, 0.5, True)
       
      if (ave .eq. 1) then

      ; Total RHS
        rhsT_av       = rhs1_av + rhs2_av + rhs3_av + rhs4_av + rhs5a_av + rhs5b_av
        rhsT_av!0     = "lev"
        rhsT_av&lev   = hyb_plot
        rhsT_av!1     = "rad"
        rhsT_av&rad   = radii
        rhsT_av@units = "m s~S~-1~N~ h~S~-1~N~"

      ; Residual                                                                        
        res_av        = rhsT_av - lhs1_av
        res_av!0      = "lev"
        res_av&lev    = hyb_plot
        res_av!1      = "rad"
        res_av&rad    = radii
        res_av@units  = "m s~S~-1~N~ h~S~-1~N~"

      ; Total diffusive tendency                                                        
        rhs5_av       = rhs5a_av + rhs5b_av
        rhs5_av!0     = "lev"
        rhs5_av&lev   = hyb_plot
        rhs5_av!1     = "rad"
        rhs5_av&rad   = radii
        rhs5_av@units = "m s~S~-1~N~ h~S~-1~N~"

      ; All mean terms
        mean_av       = rhs1_av + rhs2_av + rhs5a_av + rhs5b_av
        mean_av!0     = "lev"
        mean_av&lev   = hyb_plot
        mean_av!1     = "rad"
        mean_av&rad   = radii
        mean_av@units = "m s~S~-1~N~ h~S~-1~N~"

      ; All eddy terms
        eddy_av       = rhs3_av + rhs4_av
        eddy_av!0     = "lev"
        eddy_av&lev   = hyb_plot
        eddy_av!1     = "rad"
        eddy_av&rad   = radii
        eddy_av@units = "m s~S~-1~N~ h~S~-1~N~"

      ; Smooth data 
        lhs1       = smth9_Wrap(lhs1, 0.5, 0.5, True)

        rhs1       = smth9_Wrap(rhs1, 0.5, 0.5, True)
        rhs2       = smth9_Wrap(rhs2, 0.5, 0.5, True)
        rhs3       = smth9_Wrap(rhs3, 0.5, 0.5, True)
        rhs4       = smth9_Wrap(rhs4, 0.5, 0.5, True)
        rhs5a      = smth9_Wrap(rhs5a, 0.5, 0.5, True)
        rhs5b      = smth9_Wrap(rhs5b, 0.5, 0.5, True)

      end if

;==================================================================
; Before plotting, write out arrays to text file (write_matrix)
;==================================================================

      if (wrt .eq. 1) then

      ; Title resources for 'write_matrix' procedure
        opts_mat       = True
        opts_mat@row   = True
        opts_mat@title = ""
        opt_arr        = (/"lhs1", "rhs1", "rhs2", \
                           "rhs3", "rhs4", "rhs5a", "rhs5b"/)
        mbgt_arr       = (/lhs1_av, rhs1_av, rhs2_av, \
                           rhs3_av, rhs4_av, rhs5a_av, rhs5b_av/)

      ; Additional metadata for procedure
        ncol      = dimsizes(hyb_plot)          ; No. of cols (model levels)
        nrow      = dimsizes(radii)             ; No. of rows (radial circles)
        fmtf      = nrow + "f10.3"              ; String format (Fortran style - 10 digits; 3 dpts)

        print_clock("Writing 2D arrays to matrices for further analysis...")

      ; Call procedure and write time-averaged arrays
        do im = 0, 6
          opts_mat@fout = "./mbgt_tan_"+opt_arr(im)+"_"+dat+"_"+ens0+"_"+ts0+"_"+tf0+".txt"
          print_clock("Working on budget term: "+opt_arr(im))
          write_matrix( mbgt_arr(im,:,:), fmtf, opts_mat )
        end do

      ; Now exit script before plotting, to save time
        exit()

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

    if (clr .eq. 1) then
     opts_mbgt@cnFillPalette               = "ncl_default"
     opts_mbgt@cnFillColors                = (/2,8,16,28,40,52,64,76,88,100,118,-1,-1,\
                                               134,150,166,178,186,194,\
                                               202,214,226,238,250/)
     opts_mbgt@cnLevels                    = (/-150., -100., -50., -20., -10., -5., \
                                                -3., -2., -1., -0.5, -0.2, \
                                                0, 0.2, 0.5, 1., 2., 3., \
                                                5., 10., 20., 50., 100., 150./)
    elseif (clr .eq. 2) then
     opts_mbgt@cnFillPalette               = "bdr_extra"
     opts_mbgt@cnFillColors                = (/0,1,2,3,4,5,6,7,8,9,10,-1,-1,\
                                               13,14,15,16,17,18,19,20,21,22,23/)
     opts_mbgt@cnLevels                    = (/-150., -100., -50., -20., -10., -5., \
                                                -3., -2., -1., -0.5, -0.2, \
                                                0, 0.2, 0.5, 1., 2., 3., \
                                                5., 10., 20., 50., 100., 150./)
    elseif (clr .eq. 3) then
     opts_mbgt@cnFillPalette               = "amwg256"
     opts_mbgt@cnFillColors                = (/0,6,12,18,27,36,45,54,66,78,90,-1,-1,\
                                               148,156,166,178,186,194,\
                                               202,214,226,238,250/)
     opts_mbgt@cnLevels                    = (/-150., -100., -50., -20., -10., -5., \
                                                -3., -2., -1., -0.5, -0.2, \
                                                0, 0.2, 0.5, 1., 2., 3., \
                                                5., 10., 20., 50., 100., 150./)
    elseif (clr .eq. 4) then
     opts_mbgt@cnFillPalette               = "GMT_polar"
     opts_mbgt@cnFillColors                = (/3,5,7,-1,11,14,18/)
     opts_mbgt@cnLevels                    = (/-20,-5,0,0.5,5,20/)
    elseif (clr .eq. 5) then 
     opts_mbgt@cnFillPalette               = "BlueDarkRed18"
     if (ave .eq. 0) then 
      opts_mbgt@cnLevels                    = (/-20.0, -15.0, -10.0, -5.0, -2.0, -1.0, -0.5, \
                                                 0, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0/)
      opts_mbgt@cnFillColors                = (/0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17/)
     elseif (ave .eq. 1) then 
      opts_mbgt@cnLevels                    = (/-10.0, -5.0, -3.0, -2.0, -1.0, -0.5, -0.2, \
                                                 0, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0/)
      opts_mbgt@cnFillColors                = (/0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17/)
     end if 
    end if

    opts_mbgt@gsnPaperOrientation          = "landscape"
    opts_mbgt@tiMainString                 = ""
    opts_mbgt@tiMainFontHeightF            = 0.0125
    opts_mbgt@gsnLeftString                = ""
    opts_mbgt@gsnRightString               = ""
    opts_mbgt@gsnMaximize                  = False
    opts_mbgt@lbLabelBarOn                 = False
    opts_mbgt@gsnDraw                      = False
    opts_mbgt@gsnFrame                     = False

  ; Additional plotting resources
    opts_mbgt@tiYAxisString                = "Height (km)"
    opts_mbgt@tiXAxisString                = "Radial distance from cyclone centre (km)"

  ; Turn off tickmarks on top and right axes
    opts_mbgt@tmXTOn                       = "False"
    opts_mbgt@tmYROn                       = "False"
    opts_mbgt@gsnMaximize                  = False
    opts_mbgt@gsnAddCyclic                 = False
    opts_mbgt@tmXBMode                     = "Explicit"

    if (rad0 .eq. 0) then 
     opts_mbgt@tmXBValues                  = (/0.0,0.25,0.5,0.75,1.0,1.25,1.5/)
     opts_mbgt@tmXBLabels                  = (/"0","25","50","75","100","125","150"/)
    elseif (rad0 .eq. 1) then 
     opts_mbgt@tmXBValues                  = (/0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0/)
     opts_mbgt@tmXBLabels                  = (/"0","25","50","75","100","125","150","175","200"/)
    end if 

    opts_mbgt@tmYLMode			   = "Explicit"
    opts_mbgt@tmYLValues		   = fspan(2.0, 18.0, 9)
    opts_mbgt@tmYLLabels		   = (/"2.0","4.0","6.0","8.0","10.0",\
    					       "12.0","14.0","16.0","18.0"/)

  ; Linearise y-axis; set min and max
    opts_mbgt@gsnYAxisIrregular2Linear	   = True
    opts_mbgt@trYMinF			   = 0.0
    opts_mbgt@trYMaxF			   = 17.6

  ; Resources to overlay radial wind (following M18)
    opts_vrad                              = True
    opts_vrad@cnFillOn                     = False
    opts_vrad@cnLineColor                  = "black"
    opts_vrad@cnLevelSelectionMode         = "ExplicitLevels"
    opts_vrad@cnLevels                     = (/-1.2, 1.2/)
    opts_vrad@gsnContourNegLineDashPattern = 16
    opts_vrad@gsnContourPosLineDashPattern = 0
    opts_vrad@cnLineThicknessF             = 2.5

    opts_vrad@cnInfoLabelOn                = False
    opts_vrad@gsnDraw                      = False
    opts_vrad@gsnFrame                     = False
    opts_vrad@cnLineLabelsOn               = False
    opts_vrad@tiMainString                 = ""
    opts_vrad@gsnLeftString                = ""
    opts_vrad@gsnRightString               = ""

  ; Resources to overlay vertical velocity (following M18)
    opts_vvel                              = True
    opts_vvel@cnFillOn                     = False

  ; Set contour line colour to match M18 if we choose their colour scheme
    if (clr .eq. 4 .or. clr .eq. 5) then
     opts_vvel@cnLineColor                 = "yellow"
    else
     opts_vvel@cnLineColor                 = "grey45"
    end if

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

  ; Resources to overlay zero line
    opts_zero                             = True
    opts_zero@cnFillOn                    = False
    opts_zero@cnLineColor                 = "grey50"
    opts_zero@cnLineDashPattern		  = 3
    opts_zero@cnLevelSelectionMode        = "ExplicitLevels"
    opts_zero@cnLevels                    = 0.0
    opts_zero@cnLineThicknessF            = 3.0
    opts_zero@cnInfoLabelOn               = False
    opts_zero@gsnDraw                     = False       ; Do not draw the plot
    opts_zero@gsnFrame                    = False       ; Do no advance the frame
    opts_zero@gsnContourZeroLineThicknessF = 2.0
    opts_zero@cnLineLabelsOn               = False      ; Turn off line labels
    opts_zero@tiMainString                = ""
    opts_zero@gsnLeftString               = ""
    opts_zero@gsnRightString              = ""

;==================================================================
; Also add resources for mean quantities (u,v,w) when 'full' = 2
;==================================================================

  ; Reading in colour table using 'cnFillPalette' removes the first 2 colours (black/white)
  ; This means that the indices are all 2 values lower than they'd normally be

  ; Radial wind
    opts_vrad2                             = opts_mbgt
    delete(opts_vrad2@cnLevels)
    delete(opts_vrad2@cnFillColors)
    opts_vrad2@cnFillPalette               = "BlueDarkRed18"
    opts_vrad2@cnLevels                    = (/-10., -5., -3., -2., -1., -0.5, -0.2, \
                                                0, 0.2, 0.5, 1., 2., 3., 5., 10./)
    opts_vrad2@cnFillColors                = (/0,2,3,4,5,6,7,-1,-1,10,11,12,13,14,15,17/)
    opts_vrad2@lbLabelBarOn                = True
    opts_vrad2@lbBoxEndCapStyle            = "TriangleBothEnds"
    opts_vrad2@lbOrientation               = "Vertical"
    opts_vrad2@vpXF                        = 0.07
    opts_vrad2@vpYF                        = 0.90
    opts_vrad2@vpWidthF                    = 0.15
    opts_vrad2@vpHeightF                   = 0.15

  ; Tangential wind
    opts_vtan2                             = opts_mbgt
    delete(opts_vtan2@cnLevels)
    delete(opts_vtan2@cnFillColors)
    opts_vtan2@cnFillPalette               = "radar_new"

    if (new0 .eq. 1) then 
      opts_vtan2@cnLevels                    = (/5.0, 10.0, 15.0, 17.5, 20.0, \
                                                22.5, 25.0, 27.5, 30.0, \
                                                32.5, 35.0, 37.5, 40.0/)
      opts_vtan2@cnFillColors                = (/-1,0,1,2,3,4,5,6,7,8,9,10,11,12/)
    else 
      opts_vtan2@cnLevels                    = (/5.0, 10.0, 15.0, 20.0,\
                                                25.0, 30.0, 35.0, 40.0, \
                                                45.0, 50.0, 55.0, 60.0, \
                                                65.0/)
      opts_vtan2@cnFillColors                = (/-1,0,1,2,3,4,5,6,7,8,9,10,11,12/)
    end if 

    opts_vtan2@lbLabelBarOn                = True
    opts_vtan2@lbBoxEndCapStyle            = "TriangleBothEnds"
    opts_vtan2@lbOrientation               = "Vertical"
    opts_vtan2@vpXF                        = 0.31
    opts_vtan2@vpYF                        = 0.90
    opts_vtan2@vpWidthF                    = 0.15
    opts_vtan2@vpHeightF                   = 0.15

  ; Vertical velocity
    opts_vvel2                             = opts_mbgt
    delete(opts_vvel2@cnLevels)
    delete(opts_vvel2@cnFillColors)
    opts_vvel2@cnFillPalette               = "n11"
    opts_vvel2@cnLevels                    = (/0.1, 0.2, 0.3, 0.4, 0.5, 0.6, \
                                               0.7, 0.8, 1.0, 1.2, 1.5/)
    opts_vvel2@cnFillColors                = (/-1,0,1,2,3,4,5,6,7,8,9,10/)
    opts_vvel2@lbLabelBarOn                = True
    opts_vvel2@lbBoxEndCapStyle            = "TriangleBothEnds"
    opts_vvel2@lbOrientation               = "Vertical"
    opts_vvel2@vpXF                        = 0.55
    opts_vvel2@vpYF                        = 0.90
    opts_vvel2@vpWidthF                    = 0.15
    opts_vvel2@vpHeightF                   = 0.15

  ; Absolute angular momentum
    opts_aam                               = opts_vtan2
    delete(opts_aam@cnLevels)
    delete(opts_aam@cnFillColors)

    if (new0 .eq. 2) then
      opts_aam@cnLevels                      = (/0.5, 1.0, 1.25, 1.5, 1.75, 2.0, \
                                                 2.25, 2.5, 2.75, 3.0, 3.25, 3.5/)
      opts_aam@cnFillColors                  = (/-1,0,1,2,3,4,5,6,7,8,9,10,11/)
    elseif (new0 .eq. 1) then
      opts_aam@cnLevels                      = (/1.0, 1.5, 1.75, 2.0, \
                                                 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 4.0, 4.5/)
      opts_aam@cnFillColors                  = (/-1,0,1,2,3,4,5,6,7,8,9,10,11/)
    else
      opts_aam@cnLevels                      = (/1.0, 1.5, 2.0, 2.5, 3.0, \
                                                 3.5, 4.0, 4.5, 5.0, \
                                                 5.5, 6.0, 6.5/)
      opts_aam@cnFillColors                  = (/-1,0,1,2,3,4,5,6,7,8,9,10,11/)
    end if

    opts_aam@lbLabelBarOn                  = True
    opts_aam@lbBoxEndCapStyle              = "TriangleBothEnds"
    opts_aam@lbOrientation                 = "Vertical"
    opts_aam@vpXF                          = 0.79
    opts_aam@vpYF                          = 0.90
    opts_aam@vpWidthF                      = 0.15
    opts_aam@vpHeightF                     = 0.15

;===============================================================
; Resources to overlay panel labels without using 'gsn_panel'
;===============================================================

    if (full .eq. 1) then 
     txid                       = new(12,graphic)
     amid                       = new(12,graphic)
    elseif (full .eq. 0) then 
     txid  	      	        = new(2,graphic)
     amid			= new(2,graphic)
    end if 

  ; Label text resources
    txres                       = True
    txres@txPerimOn             = True
    txres@txFontHeightF         = 0.010
    txres@txBackgroundFillColor = "White"

  ; Label position resources
    amres                       = True
    amres@amParallelPosF        = 0.5                   ; Right edge
    amres@amOrthogonalPosF      = -0.5                  ; Top edge
    amres@amJust                = "TopRight"

    if (full .eq. 1 .and. res0 .eq. 0) then 
     panel_strings              = (/"a) vrad", "b) vtan", "c) vvel", "d) aam", \
                                    "e) V_mz", "f) V_mv", "g) V_ez", \
                                    "h) V_ev", "i) V_dh", "j) V_dz", \
                                    "k) All V tend", "l) V_t"/)
    elseif (full .eq. 1 .and. res0 .eq. 1) then 
     panel_strings              = (/"a) vrad", "b) vtan", "c) vvel", "d) aam", \
                                    "e) V_mz", "f) V_mv", "g) V_ez", \
                                    "h) V_ev", "i) V_d", "j) All V tend", \
                                    "k) V_t", "l) Residual"/)
    elseif (full .eq. 0) then 
     panel_strings    	        = (/"a) Mean terms","b) Eddy terms"/)
    end if 

;=================
; Plot the data 
;=================

  ; Calculating the terms from Montgomery et al. (2018) - their Eq.(4)
  ; Ignore the perturbation PGF term because it is several orders of magnitude smaller

  ; LHS (1) : d{v}/dt          - Mean tangential momentum tendency

  ; RHS (1) : -{u} * {f+vort}  - Mean radial influx of absolute vertical vorticity
  ; RHS (2) : -{w} * (d{v}/dz) - Mean vertical advection of mean tangential momentum
  ; RHS (3) : -{u' * vort'}    - Eddy radial vorticity flux
  ; RHS (4) : -{w' * (dv'/dz)} - Vertical advection of eddy tangential momentum
  ; RHS (5) : {D_v}            - Combined diffusive and PBL tendency

  if (ave .eq. 0) then 

   do it = ts0, tf0
   
    print("Working on time: "+it)

  ; Output file location and type 

    if (res0 .eq. 0) then 

     output = "$nep/nepartak/images/cyl_coords/mbgt_tan_xz_"+dat+"_"+ens0+\
              "_sc"+sc+"_"+time_arr(it)

    elseif (res0 .eq. 1) then 

     output = "$nep/nepartak/images/cyl_coords/mbgt_tan_xz_res_"+dat+"_"+ens0+\
              "_sc"+sc+"_"+time_arr(it)

    end if 
    wks     = gsn_open_wks(opt,output)

  if (full .eq. 1) then 

  ; Radial wind (filled contours)
    panel_vrad  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad2) ; [u]
    panel_zero  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_zero)  ; Zero line contour
    zero0       = ColorNegDashZeroPosContour(panel_zero,"transparent","black","transparent")
    overlay(panel_vrad, zero0)
    txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
    amid(0)     = gsn_add_annotation(panel_vrad, txid(0), amres)

  ; Tangential wind (filled contours)
    panel_vtan  = gsn_csm_contour(wks,mean_vtan(:,it,:),opts_vtan2) ; [v]
    txid(1)     = gsn_create_text(wks, panel_strings(1), txres)
    amid(1)     = gsn_add_annotation(panel_vtan, txid(1), amres)

  ; Vertical velocity (filled contours)
    panel_vvel  = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel2) ; [w]
    txid(2)     = gsn_create_text(wks, panel_strings(2), txres)
    amid(2)     = gsn_add_annotation(panel_vvel, txid(2), amres)

  ; Absolute angular momentum (AAM; filled contours)
    panel_aam   = gsn_csm_contour(wks,mean_aam(:,it,:),opts_aam)   ; [AAM]
    txid(3)     = gsn_create_text(wks, panel_strings(3), txres)
    amid(3)     = gsn_add_annotation(panel_aam, txid(3), amres)

  ; Mean radial influx of absolute vertical vorticity [V_mzeta]
    opts_mbgt@vpXF       = 0.07
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel0      = gsn_csm_contour(wks,rhs1(:,it,:),opts_mbgt)	    ; [V_mzeta]
    cont0	= gsn_csm_contour(wks,rhs1(:,it,:),opts_zero)	    ; Draw contours

    panel0_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad)  ; [u]
    panel0_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel)  ; [w]

    overlay(panel0, cont0)                                          ; Overlay contours
    overlay(panel0, panel0_ovr)
    overlay(panel0, panel0_ovr2)
    txid(4)     = gsn_create_text(wks, panel_strings(4), txres)
    amid(4)     = gsn_add_annotation(panel0, txid(4), amres)
    delete([/cont0/])

  ; Mean vertical advection of mean tangential momentum [V_mv]
    opts_mbgt@vpXF       = 0.31
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel1     = gsn_csm_contour(wks,rhs2(:,it,:),opts_mbgt)        ; [V_mv]
    cont0      = gsn_csm_contour(wks,rhs2(:,it,:),opts_zero)        ; Draw contours

    panel1_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad)  ; [u]
    panel1_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel)  ; [w]

    overlay(panel1, cont0)                                          ; Overlay contours
    overlay(panel1, panel1_ovr)
    overlay(panel1, panel1_ovr2)
    txid(5)     = gsn_create_text(wks, panel_strings(5), txres)
    amid(5)     = gsn_add_annotation(panel1, txid(5), amres)
    delete([/cont0/])

  ; Eddy radial vorticity flux [V_ezeta]
    opts_mbgt@vpXF       = 0.55
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel2     = gsn_csm_contour(wks,rhs3(:,it,:),opts_mbgt)        ; [V_ezeta]
    cont0      = gsn_csm_contour(wks,rhs3(:,it,:),opts_zero)        ; Draw contours

    panel2_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad)  ; [u]
    panel2_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel)  ; [w]

    overlay(panel2, cont0)                                          ; Overlay contours
    overlay(panel2, panel2_ovr)
    overlay(panel2, panel2_ovr2)
    txid(6)     = gsn_create_text(wks, panel_strings(6), txres)
    amid(6)     = gsn_add_annotation(panel2, txid(6), amres)
    delete([/cont0/])

  ; Vertical advection of eddy tangential momentum [V_ev]
    opts_mbgt@vpXF       = 0.79
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel3      = gsn_csm_contour(wks,rhs4(:,it,:),opts_mbgt)       ; [V_ev]
    cont0       = gsn_csm_contour(wks,rhs4(:,it,:),opts_zero)       ; Draw contours

    panel3_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad)  ; [u]
    panel3_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel)  ; [w]

    overlay(panel3, cont0)                                          ; Overlay contours
    overlay(panel3, panel3_ovr)
    overlay(panel3, panel3_ovr2)
    txid(7)     = gsn_create_text(wks, panel_strings(7), txres)
    amid(7)     = gsn_add_annotation(panel3, txid(7), amres)
    delete([/cont0/])

    if (res0 .eq. 0) then 

    ; Horizontal diffusive tendency of tangential momentum [V_dh]
      opts_mbgt@vpXF       = 0.07
      opts_mbgt@vpYF       = 0.50
      opts_mbgt@vpWidthF   = 0.15
      opts_mbgt@vpHeightF  = 0.15
      panel4      = gsn_csm_contour(wks,rhs5a(:,it,:),opts_mbgt)      ; [V_dh]
      cont0       = gsn_csm_contour(wks,rhs5a(:,it,:),opts_zero)      ; Draw contours    
 
      panel4_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad)  ; [u]
      panel4_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel)  ; [w]
 
      overlay(panel4, cont0)                                          ; Overlay contours
      overlay(panel4, panel4_ovr)
      overlay(panel4, panel4_ovr2)
      txid(8)     = gsn_create_text(wks, panel_strings(8), txres)
      amid(8)     = gsn_add_annotation(panel4, txid(8), amres)
      delete([/cont0/])

    ; Vertical diffusive tendency of tangential momentum [V_dz]
      opts_mbgt@vpXF       = 0.31
      opts_mbgt@vpYF       = 0.50
      opts_mbgt@vpWidthF   = 0.15
      opts_mbgt@vpHeightF  = 0.15
      panel5     = gsn_csm_contour(wks,rhs5b(:,it,:),opts_mbgt)       ; [V_dz]
      cont0      = gsn_csm_contour(wks,rhs5b(:,it,:),opts_zero)       ; Draw contours
     
      panel5_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad)  ; [u]
      panel5_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel)  ; [w] 

      overlay(panel5, cont0)                                          ; Overlay contours
      overlay(panel5, panel5_ovr)
      overlay(panel5, panel5_ovr2)
      txid(9)     = gsn_create_text(wks, panel_strings(9), txres)
      amid(9)     = gsn_add_annotation(panel5, txid(9), amres)
      delete([/cont0/])

    ; Total RHS [All V tend]
      opts_mbgt@vpXF       = 0.55
      opts_mbgt@vpYF       = 0.50
      opts_mbgt@vpWidthF   = 0.15
      opts_mbgt@vpHeightF  = 0.15
      panel6      = gsn_csm_contour(wks,rhsT(:,it,:),opts_mbgt)       ; [All V tend]
      cont0       = gsn_csm_contour(wks,rhsT(:,it,:),opts_zero)       ; Draw contours 

      panel6_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad)  ; [u]
      panel6_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel)  ; [w]

      overlay(panel6, cont0)                                          ; Overlay contours
      overlay(panel6, panel6_ovr)
      overlay(panel6, panel6_ovr2)
      txid(10)    = gsn_create_text(wks, panel_strings(10), txres)
      amid(10)    = gsn_add_annotation(panel6, txid(10), amres)
      delete([/cont0/])

    ; Local tendency of the mean tangential wind [V_t]
      opts_mbgt@vpXF       = 0.79
      opts_mbgt@vpYF       = 0.50
      opts_mbgt@vpWidthF   = 0.15
      opts_mbgt@vpHeightF  = 0.15
      panel7      = gsn_csm_contour(wks,lhs1(:,it,:),opts_mbgt)       ; [V_t]
      cont0       = gsn_csm_contour(wks,lhs1(:,it,:),opts_zero)       ; Draw contours

      panel7_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad)  ; [u]
      panel7_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel)  ; [w]

      overlay(panel7, cont0)                                          ; Overlay contours
      overlay(panel7, panel7_ovr)
      overlay(panel7, panel7_ovr2)
      txid(11)    = gsn_create_text(wks, panel_strings(11), txres)
      amid(11)    = gsn_add_annotation(panel7, txid(11), amres)
      delete([/cont0/])

    elseif (res0 .eq. 1) then 

    ; Total diffusive tendency of tangential momentum [V_d]                               
      opts_mbgt@vpXF       = 0.07
      opts_mbgt@vpYF       = 0.50
      opts_mbgt@vpWidthF   = 0.15
      opts_mbgt@vpHeightF  = 0.15
      panel4      = gsn_csm_contour(wks,rhs5(:,it,:),opts_mbgt)       ; [V_d]             
      cont0       = gsn_csm_contour(wks,rhs5(:,it,:),opts_zero)       ; Draw contours     

      panel4_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad)  ; [u]               
      panel4_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel)  ; [w]               

      overlay(panel4, cont0)                                          ; Overlay contours  
      overlay(panel4, panel4_ovr)
      overlay(panel4, panel4_ovr2)
      txid(8)     = gsn_create_text(wks, panel_strings(8), txres)
      amid(8)     = gsn_add_annotation(panel4, txid(8), amres)
      delete([/cont0/])

    ; Total RHS [All V tend]                                                              
      opts_mbgt@vpXF       = 0.31
      opts_mbgt@vpYF       = 0.50
      opts_mbgt@vpWidthF   = 0.15
      opts_mbgt@vpHeightF  = 0.15
      panel5     = gsn_csm_contour(wks,rhsT(:,it,:),opts_mbgt)        ; [All V tend]      
      cont0      = gsn_csm_contour(wks,rhsT(:,it,:),opts_zero)        ; Draw contours     

      panel5_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad)  ; [u]               
      panel5_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel)  ; [w]               

      overlay(panel5, cont0)                                          ; Overlay contours  
      overlay(panel5, panel5_ovr)
      overlay(panel5, panel5_ovr2)
      txid(9)     = gsn_create_text(wks, panel_strings(9), txres)
      amid(9)     = gsn_add_annotation(panel5, txid(9), amres)
      delete([/cont0/])

    ; Local tendency of the mean tangential wind [V_t]                                    
      opts_mbgt@vpXF       = 0.55
      opts_mbgt@vpYF       = 0.50
      opts_mbgt@vpWidthF   = 0.15
      opts_mbgt@vpHeightF  = 0.15
      panel6      = gsn_csm_contour(wks,lhs1(:,it,:),opts_mbgt)       ; [V_t]             
      cont0       = gsn_csm_contour(wks,lhs1(:,it,:),opts_zero)       ; Draw contours     

      panel6_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad)  ; [u]               
      panel6_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel)  ; [w]               

      overlay(panel6, cont0)                                          ; Overlay contours  
      overlay(panel6, panel6_ovr)
      overlay(panel6, panel6_ovr2)
      txid(10)    = gsn_create_text(wks, panel_strings(10), txres)
      amid(10)    = gsn_add_annotation(panel6, txid(10), amres)
      delete([/cont0/])

    ; Residual (LHS - RHS)                                                                
      opts_mbgt@vpXF       = 0.79
      opts_mbgt@vpYF       = 0.50
      opts_mbgt@vpWidthF   = 0.15
      opts_mbgt@vpHeightF  = 0.15
      panel7      = gsn_csm_contour(wks,res(:,it,:),opts_mbgt)        ; [Residual]        
      cont0       = gsn_csm_contour(wks,res(:,it,:),opts_zero)        ; Draw contours     

      panel7_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad)  ; [u]               
      panel7_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel)  ; [w]               

      overlay(panel7, cont0)                                          ; Overlay contours  
      overlay(panel7, panel7_ovr)
      overlay(panel7, panel7_ovr2)
      txid(11)    = gsn_create_text(wks, panel_strings(11), txres)
      amid(11)    = gsn_add_annotation(panel7, txid(11), amres)
      delete([/cont0/])

    end if 

  elseif (full .eq. 0) then

  ; All eddy terms
    opts_mbgt@vpXF       = 0.05
    opts_mbgt@vpYF       = 0.75
    opts_mbgt@vpWidthF   = 0.42
    opts_mbgt@vpHeightF  = 0.42
    panel_eddy      = gsn_csm_contour(wks,eddy(:,it,:),opts_mbgt)      ; [All U tend]
    cont0      	    = gsn_csm_contour(wks,eddy(:,it,:),opts_zero)      ; Draw contours

    panel_eddy_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad) ; [U]
    panel_eddy_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel) ; [w]

    overlay(panel_eddy, cont0)                                         ; Overlay contours
    overlay(panel_eddy, panel_eddy_ovr)
    overlay(panel_eddy, panel_eddy_ovr2)
    txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
    amid(0)     = gsn_add_annotation(panel_eddy, txid(0), amres)
    delete([/cont0/])

  ; All mean terms
    opts_mbgt@vpXF       = 0.55
    opts_mbgt@vpYF       = 0.75
    opts_mbgt@vpWidthF   = 0.42
    opts_mbgt@vpHeightF  = 0.42
    panel_mean      = gsn_csm_contour(wks,mean(:,it,:),opts_mbgt)      ; [All U tend]
    cont0      	    = gsn_csm_contour(wks,mean(:,it,:),opts_zero)      ; Draw contours

    panel_mean_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad) ; [U]
    panel_mean_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel) ; [w]

    overlay(panel_mean, cont0)                                         ; Overlay contours
    overlay(panel_mean, panel_mean_ovr)
    overlay(panel_mean, panel_mean_ovr2)
    txid(1)     = gsn_create_text(wks, panel_strings(1), txres)
    amid(1)     = gsn_add_annotation(panel_mean, txid(1), amres)
    delete([/cont0/])

  end if

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

  ; Create a main title at top of frame                                                  
    title_res                   = True
    title_str                   = "Tangential momentum budget: "+time_str(it)+" (T+"+it+")"
    title_res@txFontHeightF     = 0.02
    title_res@txFont            = "Helvetica"
    xpos_title                  = 0.50
    ypos_title                  = 0.96
    gsn_text_ndc(wks, title_str, xpos_title, ypos_title, title_res)

  ; Also create labelbar at bottom of page
    lopt                        = True
    lopt@vpWidthF               = 0.7
    lopt@vpHeightF              = 0.1
    lopt@lbAutoManage           = False
    lopt@lbOrientation          = "horizontal"
    lopt@lbLabelAlignment       = "InteriorEdges"
    lopt@lbFillPattern          = "SolidFill"
    lopt@lbMonoFillPattern      = True
    lopt@lbBoxEndCapStyle       = "TriangleBothEnds"

  ; Define colourmap and set labelbar colours
    gsn_define_colormap(wks,"BlueDarkRed18")
    lopt@lbFillColors           = (/2,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19/)

  ; Customise labelbar (number of boxes, strings, size, etc)
    label_str                  = (/"-20.0","-15.0","-10.0","-5.0","-2.0","-1.0","-0.5","0.0",\
                                    "0.5", "1.0", "2.0", "5.0", "10.0", "15.0", "20.0", "25.0"/)

    nlabels                     = dimsizes(label_str)
    lopt@lbLabelFontHeightF     = 0.010
    lopt@txFont                 = "Helvetica"
    lopt@lbTitleOn              = True 
    lopt@lbTitleString          = "Momentum budget ("+rhs1A@units+")"
    lopt@lbTitleFontHeightF     = 0.0125
    lopt@lbTitlePosition        = "Bottom"
    xpos_lbar                   = 0.16
    ypos_lbar                   = 0.30

  ; Use routine to overlay labelbar
    gsn_labelbar_ndc(wks, nlabels, label_str, xpos_lbar, ypos_lbar, lopt)

    pres			= True		; Resource to call function below
    maximize_output(wks,pres)			; Calls 'draw' and 'frame'

   end do     ; End time loop (do it = 0, ts, tf)

  elseif (ave .eq. 1) then 

  ; Define output file path and type
    if (full .eq. 0) then

     output = "$nep/nepartak/images/cyl_coords/mbgt_tan_ave_xz_em_"+dat+"_"+ens0+\
              "_"+ts0+"_"+tf0+"_sc"+sc

    elseif (full .eq. 1 .and. res0 .eq. 1) then
  
     output = "$nep/nepartak/images/cyl_coords/mbgt_tan_ave_xz_res_"+dat+"_"+ens0+\
              "_"+ts0+"_"+tf0+"_sc"+sc   

    elseif (full .eq. 1 .and. res0 .eq. 0) then 

     output = "$nep/nepartak/images/cyl_coords/mbgt_tan_ave_xz_"+dat+"_"+ens0+\
              "_"+ts0+"_"+tf0+"_sc"+sc

    end if

    wks = gsn_open_wks(opt,output)

   if (full .eq. 1) then 

  ; Radial wind (filled contours)
    panel_vrad  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad2) ; [u]
    panel_zero  = gsn_csm_contour(wks,vrad_av0(:,:),opts_zero)  ; Zero line contour
    zero0       = ColorNegDashZeroPosContour(panel_zero,"transparent","black","transparent")
    overlay(panel_vrad, zero0)
    txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
    amid(0)     = gsn_add_annotation(panel_vrad, txid(0), amres)

  ; Tangential wind (filled contours)
    panel_vtan  = gsn_csm_contour(wks,vtan_av0(:,:),opts_vtan2) ; [v]
    txid(1)     = gsn_create_text(wks, panel_strings(1), txres)
    amid(1)     = gsn_add_annotation(panel_vtan, txid(1), amres)

  ; Vertical velocity (filled contours)
    panel_vvel  = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel2) ; [w]
    txid(2)     = gsn_create_text(wks, panel_strings(2), txres)
    amid(2)     = gsn_add_annotation(panel_vvel, txid(2), amres)

  ; Absolute angular momentum (filled contours)
    panel_aam   = gsn_csm_contour(wks,aam_av0(:,:),opts_aam)    ; [AAM]
    txid(3)     = gsn_create_text(wks, panel_strings(3), txres)
    amid(3)     = gsn_add_annotation(panel_aam, txid(3), amres)


  ; Mean radial influx of absolute vertical vorticity [V_mzeta]
    opts_mbgt@vpXF       = 0.07
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel0      = gsn_csm_contour(wks,rhs1_av(:,:),opts_mbgt)   ; [V_mzeta]
    cont0       = gsn_csm_contour(wks,rhs1_av(:,:),opts_zero)   ; Draw contours

    panel0_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad)  ; [u]
    panel0_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel)  ; [w]

    overlay(panel0, cont0)                                      ; Overlay contours
    overlay(panel0, panel0_ovr)
    overlay(panel0, panel0_ovr2)
    txid(4)     = gsn_create_text(wks, panel_strings(4), txres)
    amid(4)     = gsn_add_annotation(panel0, txid(4), amres)
    delete([/cont0/])

  ; Mean vertical advection of mean tangential momentum [V_mv]
    opts_mbgt@vpXF       = 0.31
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel1     = gsn_csm_contour(wks,rhs2_av(:,:),opts_mbgt)    ; [V_mv]
    cont0      = gsn_csm_contour(wks,rhs2_av(:,:),opts_zero)    ; Draw contours

    panel1_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad)  ; [u]
    panel1_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel)  ; [w]

    overlay(panel1, cont0)                                          ; Overlay contours
    overlay(panel1, panel1_ovr)
    overlay(panel1, panel1_ovr2)
    txid(5)     = gsn_create_text(wks, panel_strings(5), txres)
    amid(5)     = gsn_add_annotation(panel1, txid(5), amres)
    delete([/cont0/])

  ; Eddy radial vorticity flux [V_ezeta]
    opts_mbgt@vpXF       = 0.55
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel2     = gsn_csm_contour(wks,rhs3_av(:,:),opts_mbgt)    ; [V_ezeta]
    cont0      = gsn_csm_contour(wks,rhs3_av(:,:),opts_zero)    ; Draw contours

    panel2_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad)  ; [u]
    panel2_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel)  ; [w]

    overlay(panel2, cont0)                                      ; Overlay contours
    overlay(panel2, panel2_ovr)
    overlay(panel2, panel2_ovr2)
    txid(6)     = gsn_create_text(wks, panel_strings(6), txres)
    amid(6)     = gsn_add_annotation(panel2, txid(6), amres)
    delete([/cont0/])

  ; Vertical advection of eddy tangential momentum [V_ev]
    opts_mbgt@vpXF       = 0.79
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel3     = gsn_csm_contour(wks,rhs4_av(:,:),opts_mbgt)    ; [V_ev]
    cont0      = gsn_csm_contour(wks,rhs4_av(:,:),opts_zero)    ; Draw contours

    panel3_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad)  ; [u]
    panel3_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel)  ; [w]

    overlay(panel3, cont0)                                      ; Overlay contours
    overlay(panel3, panel3_ovr)
    overlay(panel3, panel3_ovr2)
    txid(7)     = gsn_create_text(wks, panel_strings(7), txres)
    amid(7)     = gsn_add_annotation(panel3, txid(7), amres)
    delete([/cont0/])

    if (res0 .eq. 0) then 

    ; Horizontal diffusive tendency of tangential momentum [V_dh]
      opts_mbgt@vpXF       = 0.07
      opts_mbgt@vpYF       = 0.50
      opts_mbgt@vpWidthF   = 0.15
      opts_mbgt@vpHeightF  = 0.15
      panel4     = gsn_csm_contour(wks,rhs5a_av(:,:),opts_mbgt)   ; [V_dh]
      cont0      = gsn_csm_contour(wks,rhs5a_av(:,:),opts_zero)   ; Draw contours

      panel4_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad)  ; [u]
      panel4_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel)  ; [w]

      overlay(panel4, cont0)                                      ; Overlay contours
      overlay(panel4, panel4_ovr)
      overlay(panel4, panel4_ovr2)
      txid(8)     = gsn_create_text(wks, panel_strings(8), txres)
      amid(8)     = gsn_add_annotation(panel4, txid(8), amres)
      delete([/cont0/])

    ; Vertical diffusive tendency of tangential momentum [V_dz]
      opts_mbgt@vpXF       = 0.31
      opts_mbgt@vpYF       = 0.50
      opts_mbgt@vpWidthF   = 0.15
      opts_mbgt@vpHeightF  = 0.15
      panel5     = gsn_csm_contour(wks,rhs5b_av(:,:),opts_mbgt)   ; [V_dz]
      cont0      = gsn_csm_contour(wks,rhs5b_av(:,:),opts_zero)   ; Draw contours

      panel5_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad)  ; [u]
      panel5_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel)  ; [w]

      overlay(panel5, cont0)                                      ; Overlay contours
      overlay(panel5, panel5_ovr)
      overlay(panel5, panel5_ovr2)
      txid(9)     = gsn_create_text(wks, panel_strings(9), txres)
      amid(9)     = gsn_add_annotation(panel5, txid(9), amres)
      delete([/cont0/])

    ; Total RHS [All V tend]
      opts_mbgt@vpXF       = 0.55
      opts_mbgt@vpYF       = 0.50
      opts_mbgt@vpWidthF   = 0.15
      opts_mbgt@vpHeightF  = 0.15
      panel6      = gsn_csm_contour(wks,rhsT_av(:,:),opts_mbgt)   ; [All V tend]
      cont0       = gsn_csm_contour(wks,rhsT_av(:,:),opts_zero)   ; Draw contours

      panel6_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad)  ; [u]
      panel6_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel)  ; [w]

      overlay(panel6, cont0)                                      ; Overlay contours
      overlay(panel6, panel6_ovr)
      overlay(panel6, panel6_ovr2)
      txid(10)    = gsn_create_text(wks, panel_strings(10), txres)
      amid(10)    = gsn_add_annotation(panel6, txid(10), amres)
      delete([/cont0/])

    ; Local tendency of the mean tangential wind [V_t]
      opts_mbgt@vpXF       = 0.79
      opts_mbgt@vpYF       = 0.50
      opts_mbgt@vpWidthF   = 0.15
      opts_mbgt@vpHeightF  = 0.15
      panel7     = gsn_csm_contour(wks,lhs1_av(:,:),opts_mbgt)    ; [V_t]
      cont0      = gsn_csm_contour(wks,lhs1_av(:,:),opts_zero)    ; Draw contours

      panel7_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad)  ; [u]
      panel7_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel)  ; [w]
  
      overlay(panel7, cont0)                                      ; Overlay contours
      overlay(panel7, panel7_ovr)
      overlay(panel7, panel7_ovr2)
      txid(11)    = gsn_create_text(wks, panel_strings(11), txres)
      amid(11)    = gsn_add_annotation(panel7, txid(11), amres)
      delete([/cont0/])

    elseif (res0 .eq. 1) then 

    ; Total diffusive tendency of tangential momentum [V_d]                                       
      opts_mbgt@vpXF       = 0.07
      opts_mbgt@vpYF       = 0.50
      opts_mbgt@vpWidthF   = 0.15
      opts_mbgt@vpHeightF  = 0.15
      panel4     = gsn_csm_contour(wks,rhs5_av(:,:),opts_mbgt)    ; [V_d]                         
      cont0      = gsn_csm_contour(wks,rhs5_av(:,:),opts_zero)    ; Draw contours                 

      panel4_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad)  ; [u]                           
      panel4_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel)  ; [w]                           

      overlay(panel4, cont0)                                      ; Overlay contours              
      overlay(panel4, panel4_ovr)
      overlay(panel4, panel4_ovr2)
      txid(8)     = gsn_create_text(wks, panel_strings(8), txres)
      amid(8)     = gsn_add_annotation(panel4, txid(8), amres)
      delete([/cont0/])

    ; Total RHS [All V tend]                                                                      
      opts_mbgt@vpXF       = 0.31
      opts_mbgt@vpYF       = 0.50
      opts_mbgt@vpWidthF   = 0.15
      opts_mbgt@vpHeightF  = 0.15
      panel5     = gsn_csm_contour(wks,rhsT_av(:,:),opts_mbgt)    ; [All V tend]                  
      cont0      = gsn_csm_contour(wks,rhsT_av(:,:),opts_zero)    ; Draw contours                 

      panel5_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad)  ; [u]                           
      panel5_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel)  ; [w]                           

      overlay(panel5, cont0)                                      ; Overlay contours              
      overlay(panel5, panel5_ovr)
      overlay(panel5, panel5_ovr2)
      txid(9)     = gsn_create_text(wks, panel_strings(9), txres)
      amid(9)     = gsn_add_annotation(panel5, txid(9), amres)
      delete([/cont0/])

    ; Local tendency of the mean tangential wind [V_t]                                            
      opts_mbgt@vpXF       = 0.55
      opts_mbgt@vpYF       = 0.50
      opts_mbgt@vpWidthF   = 0.15
      opts_mbgt@vpHeightF  = 0.15
      panel6      = gsn_csm_contour(wks,lhs1_av(:,:),opts_mbgt)   ; [All V tend]                  
      cont0       = gsn_csm_contour(wks,lhs1_av(:,:),opts_zero)   ; Draw contours                 

      panel6_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad)  ; [u]                           
      panel6_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel)  ; [w]                           

      overlay(panel6, cont0)                                      ; Overlay contours              
      overlay(panel6, panel6_ovr)
      overlay(panel6, panel6_ovr2)
      txid(10)    = gsn_create_text(wks, panel_strings(10), txres)
      amid(10)    = gsn_add_annotation(panel6, txid(10), amres)
      delete([/cont0/])

    ; Residual (LHS - RHS)                                                                        
      opts_mbgt@vpXF       = 0.79
      opts_mbgt@vpYF       = 0.50
      opts_mbgt@vpWidthF   = 0.15
      opts_mbgt@vpHeightF  = 0.15
      panel7     = gsn_csm_contour(wks,res_av(:,:),opts_mbgt)     ; [Residual]                    
      cont0      = gsn_csm_contour(wks,res_av(:,:),opts_zero)     ; Draw contours                 

      panel7_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad)  ; [u]                           
      panel7_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel)  ; [w]                           

      overlay(panel7, cont0)                                      ; Overlay contours              
      overlay(panel7, panel7_ovr)
      overlay(panel7, panel7_ovr2)
      txid(11)    = gsn_create_text(wks, panel_strings(11), txres)
      amid(11)    = gsn_add_annotation(panel7, txid(11), amres)
      delete([/cont0/])

    end if 

   elseif (full .eq. 0) then

  ; All eddy terms
    opts_mbgt@vpXF       = 0.05
    opts_mbgt@vpYF       = 0.75
    opts_mbgt@vpWidthF   = 0.42
    opts_mbgt@vpHeightF  = 0.42
    panel_eddy      = gsn_csm_contour(wks,eddy_av(:,:),opts_mbgt)  ; [All U tend]
    cont0      	    = gsn_csm_contour(wks,eddy_av(:,:),opts_zero)  ; Draw contours

    panel_eddy_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad) ; [U]
    panel_eddy_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel) ; [w]

    overlay(panel_eddy, cont0)                                     ; Overlay contours
    overlay(panel_eddy, panel_eddy_ovr)
    overlay(panel_eddy, panel_eddy_ovr2)
    txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
    amid(0)     = gsn_add_annotation(panel_eddy, txid(0), amres)
    delete([/cont0/])

  ; All mean terms
    opts_mbgt@vpXF       = 0.55
    opts_mbgt@vpYF       = 0.75
    opts_mbgt@vpWidthF   = 0.42
    opts_mbgt@vpHeightF  = 0.42
    panel_mean      = gsn_csm_contour(wks,mean_av(:,:),opts_mbgt)  ; [All U tend]
    cont0      	    = gsn_csm_contour(wks,mean_av(:,:),opts_zero)  ; Draw contours

    panel_mean_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad) ; [U]
    panel_mean_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel) ; [w]

    overlay(panel_mean, cont0)                                     ; Overlay contours
    overlay(panel_mean, panel_mean_ovr)
    overlay(panel_mean, panel_mean_ovr2)
    txid(1)     = gsn_create_text(wks, panel_strings(1), txres)
    amid(1)     = gsn_add_annotation(panel_mean, txid(1), amres)
    delete([/cont0/])

   end if

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

  ; Create a main title at top of frame
    title_res                   = True 
    title_str                   = "Tangential momentum budget: "+time_str(ts0)+" to "\
                                  +time_str(tf0)
    title_res@txFontHeightF     = 0.02
    title_res@txFont            = "Helvetica"
    xpos_title                  = 0.50
    ypos_title                  = 0.96
    gsn_text_ndc(wks, title_str, xpos_title, ypos_title, title_res)

  ; Also create labelbar at bottom of page
    lopt                        = True
    lopt@vpWidthF               = 0.7
    lopt@vpHeightF              = 0.1
    lopt@lbAutoManage           = False
    lopt@lbOrientation          = "horizontal"
    lopt@lbLabelAlignment       = "InteriorEdges"
    lopt@lbFillPattern          = "SolidFill"
    lopt@lbMonoFillPattern      = True
    lopt@lbBoxEndCapStyle       = "TriangleBothEnds"

  ; Define colourmap and set labelbar colours
    gsn_define_colormap(wks,"BlueDarkRed18")
    lopt@lbFillColors           = (/2,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19/)

  ; Customise labelbar (number of boxes, strings, size, etc)
    label_str                  = (/"-10.0","-5.0","-3.0","-2.0","-1.0","-0.5","-0.2", "0.0",\
                                    "0.2", "0.5", "1.0", "2.0", "3.0", "5.0", "10.0", "15.0"/)

    nlabels                     = dimsizes(label_str)
    lopt@lbLabelFontHeightF     = 0.010
    lopt@txFont                 = "Helvetica"
    lopt@lbTitleOn              = True
    lopt@lbTitleString          = "Momentum budget ("+rhs1A@units+")"
    lopt@lbTitleFontHeightF     = 0.0125
    lopt@lbTitlePosition        = "Bottom"
    xpos_lbar                   = 0.16
    ypos_lbar                   = 0.30

  ; Use routine to overlay labelbar                                                               
    gsn_labelbar_ndc(wks, nlabels, label_str, xpos_lbar, ypos_lbar, lopt)

    pres                        = True
    maximize_output(wks,pres)                              ; Calls 'draw' and 'frame'

  end if      ; End IF statement (if (ave .eq. 0) then ...)

end