; Script to calculate the tangential momentum budget as in Montgomery et al. (2018)

; Run using:

; T+12    = 131 (144)         ; T+18  = 203 (216)
; T+24    = 275 (288)         ; T+30  = 347 (360) 
; T+36    = 419 (432)         ; T+42  = 491 (504)
; T+48    = 563 (576)         ; T+54  = 635 (648)
; T+60    = 707 (720)         ; T+66  = 779 (792)
; T+72    = 851 (864)         ; T+78  = 923 (936)
; T+84    = 995 (1008)             

; ncl dat=\"02T12\" opt=\"pdf\" ens0=\"em11\" dist=1.5 ts0=539 tf0=767
; calc=1 rad0=0 ave=0 cn0=\"slp\" diff=0 test=1 nz=73 tc_mbgt5_tan_h18.ncl

; ts0     = 539 (T+46) --> tf0 = 767 (T+65)  [FIRST VC -- 48, 51, 54, 57, 60, 63]
; ts0     = 791 (T+67) --> tf0 = 971 (T+82)  [SECOND VC -- 68, 71, 74, 77, 80]

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'opt'   = output file format ("pdf" or "x11")
; 'ens0'  = ensemble member (em00 -- em11)
; 'dist'  = size of box following storm (degrees)
; 'ts0'   = starting time index (0-119)
; 'tf0'   = ending time index (0-119)
; 'calc'  = centre calculated on each level (0), calculated offline (1) or using 980 m vort (2)
; 'rad0'  = radius of final plots: 100 km (0), 150 km (1) or 200 km (2)
; 'ave'   = calculate 3/6/12-h average (1) rather than analyse individual time interval (0)
; 'diff'  = calculate 1h-averaged diffusivity/wind stress (1) or every 5 mins (0) 
; 'test'  = calculate full time averages (0) or just a small subset (1)
; 'nz'    = number of azimuth angles in 360º (73 = 5º spacing) 

; Load main NCL functions 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"

; Also load user-defined functions 
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/diff_avg.ncl"
load "$sam/ncl_func/times_arr.ncl"

begin

;=================
; Read in files 
;=================

; Find 'pc' and 'pd' files (model height level data)
  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0+"/vc"
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0

  fili_c1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pc.nc")
  fili_d1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pd.nc")
  fili_j1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pj.nc")
  fili_k1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pk.nc")

  fili_c      = diri+"/"+fili_c1
  fili_d      = diri+"/"+fili_d1
  fili_j      = diri+"/"+fili_j1
  fili_k      = diri+"/"+fili_k1

; Calculate number of times in each file using built-in function ('nc_times')
  numINFO_c   = nc_times(fili_c)
  numINFO_d   = nc_times(fili_d)
  numINFO_j   = nc_times(fili_j)
  numINFO_k   = nc_times(fili_k)
  numINPUT_c  = numINFO_c[0]            ; Number of times in 'pc' stream
  time_c      = numINFO_c[1]            ; Array of times from 'pc' stream file
  numINPUT_d  = numINFO_d[0]            ; Number of times in 'pd' stream
  time_d      = numINFO_d[1]            ; Array of times from 'pd' stream file
  numINPUT_j  = numINFO_j[0]            ; Number of times in 'pj' stream
  time_j      = numINFO_j[1]            ; Array of times from 'pj' stream file
  numINPUT_k  = numINFO_k[0]            ; Number of times in 'pk' stream
  time_k      = numINFO_k[1]            ; Array of times from 'pk' stream file

; Calculate storm motion using built-in function ('st_rm')
; 'dat'     = initialisation time (02T12, 03T00, ...)
; 'ens0'    = ensemble simulation (em00, em01, ...)
; 'diri'    = path to input files (see above)
; 'dist'    = size of box following storm (degrees)
; 'mins'    = analyse 1-h (0) or 5-min (1) data

; Call external function 'st_rm'
  storm_rel = st_rm(dat, ens0, "$ar/text/", dist, 1)
  u_cyc     = storm_rel[0]		   ; Zonal wind 
  v_cyc     = storm_rel[1]		   ; Meridional wind
  vel_cyc   = storm_rel[2]		   ; Vector wind 
  numTIMES  = storm_rel[3]                 ; Number of times in file
  lat_arr   = storm_rel[4]                 ; Storm-following domain subset (lat)
  lon_arr   = storm_rel[5]                 ; Storm-following domain subset (lon)
  centre    = storm_rel[6]		   ; Storm centre at each hour (lat/lon)

; Create arrays for data and output information
  llbox     = toint(dist*50)                   ; Calculate domain size (grid points)
  dsize     = (/llbox,llbox/)                  ; Domain size (grid points)
  d0        = dsize(0)
  d1        = dsize(1)
  ct        = ts0                              ; Counter variable
  nh	    = (tf0 - ts0) / 12		       ; Number of hours analysed (for averaging later)
  num_l     = 63                    	       ; Number of model levels
  num_t	    = numINPUT_k		       ; Number of times in 5-min data files 

  print_clock("Working on ensemble member: "+ens0)

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 1000000000
  end setvalues

;=======================================
; Choose which storm track to read in
;=======================================

  if (cn0 .eq. "slp") then
   cn = 0
  elseif (cn0 .eq. "slpf") then
   cn = 1
  end if 

;===========================================================================            
; Create date/time string arrays for all times with user-defined function               
;===========================================================================            

 ; Call external function 'times_arr'                                                   
 ; Returns 'time_str' (12 UTC 4 July), 'time_arr' (), and 'title_arr' ()                
   all_times  = times_arr(dat, ens0)

 ; Return strings e.g. '0600 UTC 04 July'                                               
   time_str   = all_times[0]

 ; Return strings e.g. '04Jul_0600Z'                                                    
   time_arr   = all_times[1]

 ; Return title strings e.g. 'Valid at 0600 UTC 04 July'                                
   title_arr  = all_times[2]

;================================================================
; Calculate hourly averages from 5-min wind stress diagnostics
;================================================================

  if (diff .eq. 1) then 

  ; Call user-defined function ('diff_avg')
  ; 'dat'     = initialisation time (02T12, 03T00, ...)
  ; 'ens0'    = ensemble simulation (em00, em01, ...)
  ; 'dist'    = size of box following storm (degrees)
  ; 'num_l'   = number of vertical levels
  ; 'ts0'     = start time index 
  ; 'tf0'     = finish time index
  ; 'cn'      = use original (0) or finescale (1) MSLP to calculate storm centre

  ; Call user-defined function 'diff_avg' to output wind stress and horizontal diffusivity
    diff_smth = diff_avg(dat, ens0, dist, num_l, ts0, tf0, cn)
    taux_arr  = diff_smth[0]                 ; Wind stress (x-component)
    tauy_arr  = diff_smth[1]                 ; Wind stress (y-component)
    kmh_arr   = diff_smth[2]                 ; Horizontal diffusivity

  end if 

;==================================================================================
;;;;;;;;;;;;;;;;;;;;;;;;;;; Start of main part of script ;;;;;;;;;;;;;;;;;;;;;;;;;;
;==================================================================================

  print_clock("Starting main part of script...")

;========================================
; Define cylindrical coordinate arrays
;========================================

  if (rad0 .eq. 0) then
   nr        = 11
   ar        = 0.5
   dstr      = "50km"
  elseif (rad0 .eq. 1) then
   nr        = 16
   ar        = 0.75
   dstr     = "75km"
  elseif (rad0 .eq. 2) then
   nr        = 21
   ar        = 1.0
   dstr     = "100km"
  elseif (rad0 .eq. 3) then
   nr        = 31
   ar        = 1.5
   dstr     = "150km"
  elseif (rad0 .eq. 4) then
   nr        = 41
   ar        = 2.0
   dstr     = "200km"
  end if

  pii        = 3.14159265
  radii      = fspan(0.0,ar,nr)      ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,nz)   ; Array of azimuth angles (0, 5, 10, etc)
  thetas_rad = (pii/180)*thetas      ; Degrees to radians
  rad_size   = dimsizes(radii)-1     ; Index for later calculation

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

;==================================
; Get the variables we will need
;==================================

   c  = addfile(fili_c,"r")              ; Read in 'pc' stream [u,v,vort]
   d  = addfile(fili_d,"r")              ; Read in 'pd' stream [w,p]
   j  = addfile(fili_j,"r")              ; Read in 'pj' stream [t,theta,mixing ratio]
   k  = addfile(fili_k,"r")              ; Read in 'pk' stream [tau,PV,w]

   do it = ts0, tf0			 ; Loop over times in file

    print("Working on time: "+time_str(it))

;=====================================================
; Read in variables from 'pc' stream (model levels)
;=====================================================

    t0    = lat_arr(cn,it,0)
    t1    = lat_arr(cn,it,1)
    n0    = lon_arr(cn,it,0)
    n1    = lon_arr(cn,it,1)
    print_clock("Grid (pc stream) extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

  ; Standard latitude/longitude grid
    lon    = c->longitude({n0:n1})                ; longitude ['d1' grid points]
    lat    = c->latitude({t0:t1})                 ; latitude  ['d0' grid points]

    dx0    = lon(1) - lon(0) ; Grid spacing (longitude)
    dy0    = lat(1) - lat(0) ; Grid spacing (latitude)

  ; Model levels (vertical)
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

    u   = u0 - u_cyc(cn,it)
    v   = v0 - v_cyc(cn,it)

  ; Add metadata from original horizontal wind arrays
    copy_VarCoords(u0,u)
    copy_VarAtts(u0,u)

    copy_VarCoords(v0,v)
    copy_VarAtts(v0,v)

;==================================================================
; Calculate relative vorticity using centered finite differences
;==================================================================

  ; Option '2' --> boundary points estimated using one-sided difference scheme
    vort = uv2vr_cfd(u, v, lat, lon, 2)
    copy_VarCoords(u, vort)
    vort@units = "s~S~-1~N~"
    vort@name = "Relative vorticity on model levels"

  ; Coriolis parameter and absolute vorticity
    f1   = coriolis_param(lat)
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
    hybD   = d->hybrid_ht(:)			    ; 63 model ('rho') levels

  ; Read in temperature, pressure, vertical velocity and geopotential height
  ; All variables below -- [hybD | 63] * [latD | 150] * [lon1D | 151]
    p      = d->p(it,:,{t0:t1},{n0:n1})             ; Pressure (Pa)
    z      = d->ht(it,:,{t0:t1},{n0:n1})    	    ; Geopotential height (m)

;======================================
; Read in variables from 'pj' stream
;======================================

  ; Read in basic variables from 'pj' stream
    hybJ   = j->hybrid_ht(:)                        ; 63 model ('rho') levels

  ; Read in temperature, potential temperature and mixing ratio
  ; All variables below  -- [hybJ | 63] * [latJ | 150]  * [lonJ | 150]
    t      = j->temp(it,:,{t0:t1},{n0:n1})          ; Temperature
    th     = j->theta(it,:,{t0:t1},{n0:n1})         ; Potential temperature
    mix    = j->mix(it,:,{t0:t1},{n0:n1}) 	    ; Mixing ratio

;======================================
; Read in variables from 'pk' stream
;======================================

  ; Read in basic variables from 'pk' stream
    hybK   = k->hybrid_ht(:)                        ; 63 model ('rho') levels

  ; Alternative latitude/longitude grid
    lon1   = k->longitude_1({n0:n1})                ; longitude
    lat1   = k->latitude_1({t0:t1})                 ; latitude

  ; Read in potential vorticity and vertical velocity
    w      = k->dz_dt(it,:,{t0:t1},{n0:n1})         ; Vertical velocity (m s-1)
    pv	   = k->field83(it,0:62,{t0:t1},{n0:n1})    ; Potential vorticity

  ; Read in horizontal wind stress components and diffusivity
    if (diff .eq. 0) then
     str_xz = k->taux(it,:,{t0:t1},{n0:n1})
     str_yz = k->tauy(it,:,{t0:t1},{n0:n1})
     kmh    = k->unspecified(it,:,{t0:t1},{n0:n1})
    end if

;================================================================
; Also calculate dry air density from pressure and temperature
;================================================================

    rho    = p / (rd * t)

  ; Add metadata
    copy_VarCoords(p, rho)
    rho@description = "Dry air density"
    rho@units = "kg m~S~-3~N~"

  ; Also calculate specific volume (for use in diffusive tendency source terms)
    rho_a  = 1 / rho
    copy_VarCoords(rho, rho_a)
    rho_a@description = "Specific volume"
    rho_a@units = "m~S~3~N~ kg~S~-1~N~"

  ; Now convert units of pressure and geopotential height before continuing
    p      = p / 100
    z      = z / 10

  ; Add metadata
    p@units = "hPa"
    z@units = "dam"

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

    do ilev = 0, num_l-1 ; Loop over chosen model levels

     hy = sprintf("%0.0f",hybD(ilev))
     print("Working on model level " + ilev+ " ("+hy+" m AGL)")

   ; Calculate components of vertical diffusion term
   ; accounting for top and bottom levels
     ip1 = min((/ilev+1,num_l-1/))
     im1 = max((/ilev-1,0/))

   ; Read in data as normal (diff=0) or from the pre-calculated arrays (diff=1)
     if (diff .eq. 1) then 
 
     ; Wind stress components on current model level
       tau_xz = taux_arr(it,ilev,0:d0-1,0:d1-1)
       tau_yz = tauy_arr(it,ilev,0:d0-1,0:d1-1)

     ; Wind stress components on model levels above/below current level
       tau_x1 = taux_arr(it,ip1,0:d0-1,0:d1-1)
       tau_x0 = taux_arr(it,im1,0:d0-1,0:d1-1)
       tau_y1 = tauy_arr(it,ip1,0:d0-1,0:d1-1)
       tau_y0 = tauy_arr(it,im1,0:d0-1,0:d1-1)
       dz     = hybD(ip1) - hybD(im1)

     elseif (diff .eq. 0) then 

     ; Wind stress components on current model level
       tau_xz = str_xz(ilev,0:d0-1,0:d1-1)
       tau_yz = str_yz(ilev,0:d0-1,0:d1-1)

     ; Wind stress components on model levels above/below current level 
       tau_x1 = str_xz(ip1,0:d0-1,0:d1-1)
       tau_x0 = str_xz(im1,0:d0-1,0:d1-1)
       tau_y1 = str_yz(ip1,0:d0-1,0:d1-1)
       tau_y0 = str_yz(im1,0:d0-1,0:d1-1)
       dz     = hybD(ip1) - hybD(im1)       

     end if 

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

   ; Read in other variables on model levels
     u_plane 	            = u(ilev,:,:)	   ; Zonal wind 
     v_plane		    = v(ilev,:,:)	   ; Meridional wind 
     w_plane		    = w(ilev,:,:)	   ; Vertical velocity
     f_plane		    = f0(ilev,:,:)	   ; Coriolis parameter
     prs_plane		    = p(ilev,:,:)	   ; Pressure 
     rho_plane              = rho(ilev,:,:)	   ; Density
     vort_plane	   	    = vort(ilev,:,:)	   ; Relative vorticity
     avo_plane		    = avo(ilev,:,:)	   ; Absolute vorticity

   ; Read in horizontal diffusivity as normal (diff=0) or from pre-calculated array (diff=1)
     if (diff .eq. 1) then 
      kmh_plane		    = kmh_arr(it,ilev,:,:)
     elseif (diff .eq. 0) then 
      kmh_plane             = kmh(ilev,:,:)    
     end if 

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

   ;==================================================================
   ; Either calculate the centre, or read in pre-calculated values
   ;==================================================================

   if (calc .eq. 0) then

   ; Find latitude and longitude of storm centre using vorticity calculation above
     lat_max0 	            = t0 + (dy0 * inds(0,0))
     lon_max0 	   	    = n0 + (dx0 * inds(0,1))

     print("Centre (vort): "+lat_max0+" degrees N, "+lon_max0+" degrees E")

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

     lat_max    = lt1 + (dy0 * inds_p0(0,0))
     lon_max    = ln1 + (dx0 * inds_p0(0,1))
     print("Centre (vort + pres): "+lat_max+" degrees N, "+lon_max+" degrees E")
     delete([/prs_plane0,prs_min0,dims_p0,prs_1d0,inds_p0/])

   ; Distance between centre calculated above using relative vorticity, and that calculated offline
     dlat       = abs(lat_max - centre(cn,it,0))
     dlon	= abs(lon_max - centre(cn,it,1))

   ; If our storm centre calculations deviate too much (> 0.75º) from the values
   ; we read in using the NCL function 'st_rm', calculate a new centre using
   ; the minimum pressure. Constrain the new calculation using the values from
   ; 'st_rm', to prevent spurious centres being identified and analysed. 

     if (dlat .gt. 1.00 .or. dlon .gt. 1.00) then 

      lt1 = centre(cn,it,0) - 0.25
      lt2 = centre(cn,it,0) + 0.25
      ln1 = centre(cn,it,1) - 0.25
      ln2 = centre(cn,it,1) + 0.25

      prs_plane0 = prs_plane({lt1:lt2},{ln1:ln2})       ; Smaller grid (1.0 x 1.0 deg)
      prs_min0   = min(prs_plane0)                      ; Find minimum on smaller grid
      dims_p0    = dimsizes(prs_plane0)                 ; Size of smaller grid
      prs_1d0    = ndtooned(prs_plane0)                 ; Create 1-D array
      inds_p0    = ind_resolve(minind(prs_1d0),dims_p0) ; Find index of minimum

      lat_max 	 = lt1 + (dy0 * inds_p0(0,0))
      lon_max 	 = ln1 + (dx0 * inds_p0(0,1))
      print("New centre: "+lat_max+" degrees N, "+lon_max+" degrees E")
      delete([/prs_plane0,prs_min0,dims_p0,prs_1d0,inds_p0/])

     end if 

   elseif (calc .eq. 1) then

      lat_max    = centre(cn,it,0)		; Does this relate to 'pc' or 'pd' stream???
      lon_max    = centre(cn,it,1)

   elseif (calc .eq. 2) then 

      lat_max    = t0 + (dy0 * inds0(0,0))
      lon_max	 = n0 + (dx0 * inds0(0,1))
      print("Centre (980 m): "+lat_max+" degrees N, "+lon_max+" degrees E")      

   end if

   delete([/vort1d,inds,vort1d_sm,vort_smth,inds_sm/])

   ;======================================================================
   ; Create lat/lon arrays for later switch to cylindrical coordinates
   ;======================================================================

   ; Indices of 'lat'/'lon' grids
     lat_plane   = ind(t0.le.lat.and.lat.le.t1)
     lon_plane   = ind(n0.le.lon.and.lon.le.n1)

   ; Array sizes
     lat_size    = dimsizes(lat_plane)
     lon_size    = dimsizes(lon_plane)
     
   ; Array of latitude/longitude points	corresponding to the indices above
     lat_sub     = lat(lat_plane(0):d0-1)
     lon_sub     = lon(lon_plane(0):d1-1)

   ; Indices of 'lat1'/'lon1' grids
     lat1_plane  = ind(t0.le.lat1.and.lat1.le.t1)
     lon1_plane  = ind(n0.le.lon1.and.lon1.le.n1)

   ; Array sizes
     lat1_size   = dimsizes(lat1_plane)
     lon1_size   = dimsizes(lon1_plane)

   ; Array of latitude/longitude points corresponding to the indices above
     lat1_sub    = lat1(lat1_plane(0):d0-1)
     lon1_sub    = lon1(lon1_plane(0):d1-1)

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
     kmh_int  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Horizontal diffusivity
     vort_int = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Relative vorticity
     avo_int  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Absolute vorticity
     rho_int  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Density
     v_rad    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Radial velocity
     v_tan    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Tangential velocity

   ; Additional arrays for horizontal and vertical diffusion
     difv_x   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vert. diff (x)
     difv_y   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vert. diff (y)
     difv_tot = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vert. diff (total)

     difh_x   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Horiz. diff (x)
     difh_y   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Horiz. diff (y)

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
   ; vort_plane	 -- latC,  lonC
   ; avo_plane	 -- latC,  lonC
   ; f_plane     -- latC,  lonC

   ; w_plane     -- latD,  lon1D
   ; prs_plane   -- latD,  lon1D
   ; rho_plane	 -- latD,  lon1D
   ; kmh_plane	 -- latD,  lon1D

   ; difv_xplane -- latD,  lonD
   ; difv_yplane -- lat1D, lon1D 

     do iang = 0, dimsizes(thetas_rad)-1
      u_int(:,iang)    = linint2_points(lon_sub,lat_sub,u_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      v_int(:,iang)    = linint2_points(lon_sub,lat_sub,v_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      vort_int(:,iang) = linint2_points(lon_sub,lat_sub,vort_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      avo_int(:,iang)  = linint2_points(lon_sub,lat_sub,avo_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      f_int(:,iang)    = linint2_points(lon_sub,lat_sub,f_plane,False,\
				        xpos(:,iang),ypos(:,iang),0)
      difv_x(:,iang)   = linint2_points(lon_sub,lat_sub,difv_xplane,False,\
					xpos(:,iang),ypos(:,iang),0)
      difv_y(:,iang)   = linint2_points(lon1_sub,lat1_sub,difv_yplane,False,\
					xpos(:,iang),ypos(:,iang),0)
      w_int(:,iang)    = linint2_points(lon1_sub,lat_sub,w_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      rho_int(:,iang)  = linint2_points(lon1_sub,lat_sub,rho_plane,False,\
					xpos(:,iang),ypos(:,iang),0)
      kmh_int(:,iang)  = linint2_points(lon1_sub,lat_sub,kmh_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
     end do

   ; Calculate radial/tangential windspeed and friction at each point (using unit vctrs)
      do iang = 0, dimsizes(thetas)-1
       theta = thetas_rad(iang)

      ; Tangential and radial wind 
        v_rad(:,iang) = u_int(:,iang)*cos(theta) + v_int(:,iang)*sin(theta)
        v_tan(:,iang) = -u_int(:,iang)*sin(theta) + v_int(:,iang)*cos(theta)

      ; Total vertical diffusive tendency of tangential momentum
        difv_tot(:,iang) = difv_y(:,iang)*cos(theta) - difv_x(:,iang)*sin(theta)

      end do

   ; Fill 'u_dz0' array with values from 'difv_tot' (vertical diffusive tendency)
     v_dz0(ilev,:,:,ct)   = difv_tot(:,:)    

    ;===================================================
    ; Calculate mean (azimuthally-averaged) variables 
    ;===================================================

      do irad = 0, dimsizes(radii)-1

        mean_vtan0(ilev,irad,:,ct)    = avg(v_tan(irad,:))      ; [ {v} ]
        mean_vrad0(ilev,irad,:,ct)    = avg(v_rad(irad,:))      ; [ {u} ]
        mean_vvel0(ilev,irad,:,ct)    = avg(w_int(irad,:))      ; [ {w} ]
	mean_vort0(ilev,irad,:,ct)    = avg(vort_int(irad,:))   ; [ {zeta} ]
        mean_avo0(ilev,irad,:,ct)     = avg(avo_int(irad,:))    ; [ {zeta + f} ]
        mean_difv0(ilev,irad,:,ct)    = avg(difv_tot(irad,:))   ; [ {diff_xz + diff_yz} ]
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
      difv_x0(ilev,:,:,ct)   = (/difv_x(:,:)/)
      difv_y0(ilev,:,:,ct)   = (/difv_y(:,:)/)
      difv_tot0(ilev,:,:,ct) = (/difv_tot(:,:)/)
      v_rad0(ilev,:,:,ct)    = (/v_rad(:,:)/)
      v_tan0(ilev,:,:,ct)    = (/v_tan(:,:)/)

    ; Tidy up
      delete([/u_int,v_int,w_int,v_tan,v_rad/])
      delete([/rho_int,difv_x,difv_y,f_int,kmh_int,vort_int,avo_int/])
      delete([/u_plane,v_plane,w_plane,f_plane,vort_plane,avo_plane/])
      delete([/difv_xplane,difv_yplane,kmh_plane,xpos,ypos,xcen,ycen/])
      delete([/lat_plane,lon_plane,lon1_plane,lat1_plane/])
      delete([/lat_sub,lon_sub,lat1_sub,lon1_sub/])

     end do			    ; End model level loop (do ilev = 0, num_l-1)

    ct  = ct + 1		    ; Counter variable (time)

  ; Tidy up before next iteration (all variables from 'pc' and 'pd' streams)
    delete([/lon,lat,lon1,lat1,u,v,u0,v0,vort,avo,f0,f1/])
    delete([/p,z/])
    delete([/t,th,mix/])
    delete([/w,pv,rho,rho_a/])

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

 ; Calculate increments
   vm   = mean_vtan0(:,:,:,:)
   fm   = mean_cor0(:,:,:,:) * (r1 ^ 2)

 ; Calculate absolute angular momentum (AAM)
   aam0 = vm * r1
   aam1 = 0.5 * fm

 ; Calculate absolute angular momentum [ rv + 0.5 * (fr ^ 2) ]
   mean_aam0(:,:,:,:)  = aam0 + aam1

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

          end do		 ; End loop over radial circles (do irad = 0, ...)

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

 ; LHS (1) : d{v}/dt          - Local tangential wind tendency

 ; RHS (1) : -{u} * {f+vort}  - Mean radial advection
 ; RHS (2) : -{w} * (d{v}/dz) - Mean vertical advection
 ; RHS (3) : -{u' * vort'}    - Eddy radial advection
 ; RHS (4) : -{w' * (dv'/dz)} - Eddy vertical advection
 ; RHS (5) : {D_v}            - Combined diffusive and PBL tendency

 ; Create array to hold tangential wind tendency at all times
   lhs1F_t = new((/numTIMES,num_l,dimsizes(radii)/),float)    ; Local tangential wind tendency (all)

 ; Now create additional arrays for final plots (6 time averages)
   if (ts0 .eq. 539) then 
    num_plots = 6
   elseif (ts0 .eq. 791) then 
    num_plots = 5
   end if 

   lhs1F   = new((/num_plots,num_l,dimsizes(radii)/),float)   ; Local tangential wind tendency

   rhs1F   = new((/num_plots,num_l,dimsizes(radii)/),float)   ; Mean radial advection
   rhs2F   = new((/num_plots,num_l,dimsizes(radii)/),float)   ; Mean vertical advection
   rhs3F   = new((/num_plots,num_l,dimsizes(radii)/),float)   ; Eddy radial advection
   rhs4F   = new((/num_plots,num_l,dimsizes(radii)/),float)   ; Eddy vertical advection
   rhs5aF  = new((/num_plots,num_l,dimsizes(radii)/),float)   ; Horizontal diffusion
   rhs5bF  = new((/num_plots,num_l,dimsizes(radii)/),float)   ; Vertical diffusion

   hyb_plot   = hybC / 1000
   hyb_plot@units = "km"

 ; Add time metadata
   lhs1F!0     = "time"
   lhs1F!1     = "lev"
   lhs1F&lev   = hyb_plot

   rhs1F!0     = "time"
   rhs1F!1     = "lev"
   rhs1F&lev   = hyb_plot

   rhs2F!0     = "time"
   rhs2F!1     = "lev"
   rhs2F&lev   = hyb_plot

   rhs3F!0     = "time"
   rhs3F!1     = "lev"
   rhs3F&lev   = hyb_plot

   rhs4F!0     = "time"
   rhs4F!1     = "lev"
   rhs4F&lev   = hyb_plot

   rhs5aF!0    = "time"
   rhs5aF!1    = "lev"
   rhs5aF&lev  = hyb_plot

   rhs5bF!0    = "time"
   rhs5bF!1    = "lev"
   rhs5bF&lev  = hyb_plot

 ; Same as above arrays but averaged in the radial direction ('radii')
   lhs1_f  = new((/num_plots,num_l/),float)   ; Local tangential wind tendency
   rhs1_f  = new((/num_plots,num_l/),float)   ; Mean radial advection
   rhs2_f  = new((/num_plots,num_l/),float)   ; Mean vertical advection
   rhs3_f  = new((/num_plots,num_l/),float)   ; Eddy radial advection
   rhs4_f  = new((/num_plots,num_l/),float)   ; Eddy vertical advection
   rhs5a_f = new((/num_plots,num_l/),float)   ; Horizontal diffusion
   rhs5b_f = new((/num_plots,num_l/),float)   ; Vertical diffusion

 ; Add metadata
   lhs1_f!0    = "time"
   lhs1_f!1    = "lev"
   lhs1_f&lev  = hyb_plot

   rhs1_f!0    = "time"
   rhs1_f!1    = "lev"
   rhs1_f&lev  = hyb_plot

   rhs2_f!0    = "time"
   rhs2_f!1    = "lev"
   rhs2_f&lev  = hyb_plot

   rhs3_f!0    = "time"
   rhs3_f!1    = "lev"
   rhs3_f&lev  = hyb_plot

   rhs4_f!0    = "time"
   rhs4_f!1    = "lev"
   rhs4_f&lev  = hyb_plot

   rhs5a_f!0   = "time"
   rhs5a_f!1   = "lev"
   rhs5a_f&lev = hyb_plot

   rhs5b_f!0   = "time"
   rhs5b_f!1   = "lev"
   rhs5b_f&lev = hyb_plot

 ; Add metadata and make units consistent for all terms
   rhs1_st = rhs1_st * 3600
   rhs1_st@units = "m s~S~-1~N~ h~S~-1~N~"

   rhs2_st = rhs2_st * 3600
   rhs2_st@units = "m s~S~-1~N~ h~S~-1~N~"
       
   rhs3_st = rhs3_st * 3600
   rhs3_st@units = "m s~S~-1~N~ h~S~-1~N~"

   rhs4_st = rhs4_st * 3600
   rhs4_st@units = "m s~S~-1~N~ h~S~-1~N~"

   rhs5a_st = rhs5a_st * 3600 
   rhs5a_st@units = "m s~S~-1~N~ h~S~-1~N~"

   rhs5b_st = rhs5b_st * 3600
   rhs5b_st@units = "m s~S~-1~N~ h~S~-1~N~"

 ; Add metadata for all the other arrays (instantaneous)
   rhs1_st!0    = "lev"
   rhs1_st&lev  = hyb_plot
   rhs1_st!2    = "time"
   rhs1_st&time = time_c
   rhs1_st!1    = "rad"
   rhs1_st&rad  = radii

   rhs2_st!0    = "lev"
   rhs2_st&lev  = hyb_plot
   rhs2_st!2    = "time"
   rhs2_st&time = time_c
   rhs2_st!1    = "rad"
   rhs2_st&rad  = radii

   rhs3_st!0    = "lev"
   rhs3_st&lev  = hyb_plot
   rhs3_st!2    = "time"
   rhs3_st&time = time_c
   rhs3_st!1    = "rad"
   rhs3_st&rad  = radii

   rhs4_st!0    = "lev"
   rhs4_st&lev  = hyb_plot
   rhs4_st!2    = "time"
   rhs4_st&time = time_c
   rhs4_st!1    = "rad"
   rhs4_st&rad  = radii

   rhs5a_st!0   = "lev"
   rhs5a_st&lev = hyb_plot
   rhs5a_st!2   = "time"
   rhs5a_st&time = time_c
   rhs5a_st!1   = "rad"
   rhs5a_st&rad = radii

   rhs5b_st!0   = "lev"
   rhs5b_st&lev = hyb_plot
   rhs5b_st!2   = "time"
   rhs5b_st&time = time_c
   rhs5b_st!1   = "rad"
   rhs5b_st&rad = radii

 ; Loop over times to calculate tangential wind tendency at all times 
   do it = ts0, tf0

    cp1 = min((/it+1,tf0/))
    cm1 = max((/it-1,ts0/))

  ; Multiply by 6 to maintain units of [m s-1 h-1] -- RERUNNING NOW TO CHECK
    lhs1F_t(it,:,:) = (mean_vtan1(:,:,cp1) - mean_vtan1(:,:,cm1)) / 2 ; * 6

   end do                ; End loop over times (do it = ts0, tf0) 

 ; Add metadata to 'lhs1F_t'
   lhs1F_t@units = "m s~S~-1~N~ h~S~-1~N~"

   lhs1F_t!0     = "time"
   lhs1F_t&time  = time_c
   lhs1F_t!1     = "lev"
   lhs1F_t&lev   = hyb_plot
   lhs1F_t!2     = "rad"
   lhs1F_t&rad   = radii

 ; Calculate 2-h averages of all terms above at 6 specified times 
   if (test .eq. 1) then 
    ti0 = (/0, 2, 4, 6, 8, 10/)
    ti1 = (/2, 4, 6, 8, 10, 12/)
   else 
    if (ts0 .eq. 539) then 
     ti0 = (/12, 48, 84, 120, 156, 192/)
     ti1 = (/36, 72, 108, 144, 180, 216/)
    elseif (ts0 .eq. 791) then 
     ti0 = (/12, 48, 84, 120, 156/)
     ti1 = (/36, 72, 108, 144, 180/)
    end if 
   end if 

   do tt = 0, num_plots-1

  ; Local tendency 
    lhs1F(tt,:,:)   = dim_sum_n_Wrap( lhs1F_t(ts0+ti0(tt):ts0+ti1(tt),:,:), 0)

  ; Source terms 
    rhs1F(tt,:,:)   = dim_avg_n_Wrap( rhs1_st(:,:,ts0+ti0(tt):ts0+ti1(tt) ), 2)
    rhs2F(tt,:,:)   = dim_avg_n_Wrap( rhs2_st(:,:,ts0+ti0(tt):ts0+ti1(tt) ), 2)
    rhs3F(tt,:,:)   = dim_avg_n_Wrap( rhs3_st(:,:,ts0+ti0(tt):ts0+ti1(tt) ), 2)
    rhs4F(tt,:,:)   = dim_avg_n_Wrap( rhs4_st(:,:,ts0+ti0(tt):ts0+ti1(tt) ), 2)
    rhs5aF(tt,:,:)  = dim_avg_n_Wrap( rhs5a_st(:,:,ts0+ti0(tt):ts0+ti1(tt) ), 2)
    rhs5bF(tt,:,:)  = dim_avg_n_Wrap( rhs5b_st(:,:,ts0+ti0(tt):ts0+ti1(tt) ), 2)

   end do 

 ; Also average over radial band to create 2D arrays 
   lhs1_f(:,:)      = dim_avg_n_Wrap( lhs1F(:,:,:), 2)
   rhs1_f(:,:)      = dim_avg_n_Wrap( rhs1F(:,:,:), 2)
   rhs2_f(:,:)      = dim_avg_n_Wrap( rhs2F(:,:,:), 2)
   rhs3_f(:,:)      = dim_avg_n_Wrap( rhs3F(:,:,:), 2)
   rhs4_f(:,:)      = dim_avg_n_Wrap( rhs4F(:,:,:), 2)
   rhs5a_f(:,:)     = dim_avg_n_Wrap( rhs5aF(:,:,:), 2)
   rhs5b_f(:,:)     = dim_avg_n_Wrap( rhs5bF(:,:,:), 2)

;=========================================
; Finally, smooth data before plotting
;=========================================

 ; Local tendency 
   lhs1           = smth9_Wrap(lhs1_f, 0.5, 0.5, True)
  
 ; Source terms 
   rhs1  	  = smth9_Wrap(rhs1_f, 0.5, 0.5, True)      ; Mean radial advection
   rhs2  	  = smth9_Wrap(rhs2_f, 0.5, 0.5, True)      ; Mean vertical advection
   rhs3  	  = smth9_Wrap(rhs3_f, 0.5, 0.5, True)      ; Eddy radial advection
   rhs4  	  = smth9_Wrap(rhs4_f, 0.5, 0.5, True)      ; Eddy vertical advection
   rhs5a 	  = smth9_Wrap(rhs5a_f, 0.5, 0.5, True)     ; Horizontal diffusion
   rhs5b 	  = smth9_Wrap(rhs5b_f, 0.5, 0.5, True)     ; Vertical diffusion

 ; Calculate total RHS tendency
   rhsT       = rhs1 + rhs2 + rhs3 + rhs4 + rhs5a + rhs5b
   rhsT!1     = "lev"
   rhsT&lev   = hyb_plot
   rhsT!0     = "time"
   rhsT@units = "m s~S~-1~N~ h~S~-1~N~"

 ; Calculate sum of all mean terms (vadv + hadv + diff)
   mean       = rhs1 + rhs2 + rhs5a + rhs5b
   mean!1     = "lev"
   mean&lev   = hyb_plot
   mean!0     = "time"
   mean@units = "m s~S~-1~N~ h~S~-1~N~"

 ; Calculate sum of all eddy terms (vadv + hadv)
   eddy       = rhs3 + rhs4
   eddy!1     = "lev"
   eddy&lev   = hyb_plot
   eddy!0     = "time"
   eddy@units = "m s~S~-1~N~ h~S~-1~N~"

 ; Total diffusive tendency 
   rhs5       = rhs5a + rhs5b
   rhs5!1     = "lev"
   rhs5&lev   = hyb_plot
   rhs5!0     = "time"
   rhs5@units = "m s~S~-1~N~ h~S~-1~N~"

 ; Total horizontal advection
   hadv       = rhs1 + rhs3 
   hadv!1     = "lev"
   hadv&lev   = hyb_plot
   hadv!0     = "time"
   hadv@units = "m s~S~-1~N~ h~S~-1~N~"

 ; Total vertical advection
   vadv       = rhs2 + rhs4
   vadv!1     = "lev"
   vadv&lev   = hyb_plot
   vadv!0     = "time"
   vadv@units = "m s~S~-1~N~ h~S~-1~N~"

;========================
; Options for plotting 
;========================

 ; Set up output file type and path
   output     = "$sam/nepartak/images/cyl_coords/mbgt5_tan_h18_"+dat+"_"+ens0+\
                "_"+ts0+"_"+tf0+"_"+cn0+"_"+dstr
   wks        = gsn_open_wks(opt,output)

 ; Set up panel plot 
   panel      = new(7,graphic)

 ; Y-axis resources (set range later in script)                                           
   opts_mbgt                              = True 
   opts_mbgt@tiYAxisString                = "Height (km)"        ; Axis title 
   opts_mbgt@tiYAxisFontHeightF           = 0.035                ; Axis title font size 
   opts_mbgt@tmYROn                       = False                ; Tick marks off (RHS)   
   opts_mbgt@tmYLLabelFontHeightF         = 0.025                ; Axis label font size 
   opts_mbgt@tmYRMode                     = "Explicit"
   opts_mbgt@trYMinF                      = 0.0                  ; Minimum value 
   opts_mbgt@trYMaxF                      = 16.0                 ; Maximum value 

 ; Set y-axis labels
   dy_axes                                = toint(opts_mbgt@trYMaxF - opts_mbgt@trYMinF)
   opts_mbgt@tmYLValues                   = fspan(opts_mbgt@trYMinF, opts_mbgt@trYMaxF, dy_axes)
   opts_mbgt@tmYLLabels                   = sprintf("%0.0f",opts_mbgt@tmYLValues)

 ; X-axis resources                                                         
 ; EDIT TO ADD LEGEND
   opts_mbgt@trXMinF                      = -9                   ; Minimum value
   opts_mbgt@trXMaxF                      = 9                    ; Maximum value
   opts_mbgt@tiXAxisOn                    = False                ; Turn off x-axis title
   opts_mbgt@tiXAxisString                = "Momentum budget ("+eddy@units+")"
   opts_mbgt@tiXAxisFontHeightF           = 0.035                ; Axis title font size 
   opts_mbgt@tmXTOn                       = False                ; Tick marks off (top)   
   opts_mbgt@tmXBLabelFontHeightF         = 0.025                ; Axis label font size
   opts_mbgt@tmXBMode                     = "Explicit"

 ; Produce labels for x-axis based on range 
   dx                                     = opts_mbgt@trXMaxF - opts_mbgt@trXMinF  ; Range 
   nl                                     = (dx / 3) + 1                           ; No. of labels
   opts_mbgt@tmXBValues                   = fspan(opts_mbgt@trXMinF, \
                                                  opts_mbgt@trXMaxF, nl)
   opts_mbgt@tmXBLabels                   = sprintf("%0.0f",opts_mbgt@tmXBValues)

 ; Additional ploting resources
   opts_mbgt@gsnDraw                      = False
   opts_mbgt@gsnFrame                     = False
   opts_mbgt@tmXMajorGridLineDashPattern  = 2
   opts_mbgt@tmYMajorGridLineDashPattern  = 2

 ; Line plot resources
   opts_mbgt@xyLineColors                 = (/"royalblue1","cyan","yellow2",\
                                              "orange2","red3","purple3"/)
   opts_mbgt@xyDashPatterns               = (/0, 0, 0, 0, 0, 0/)
   opts_mbgt@xyLineThicknessF             = 2.0            

 ; Plot each panel 
   panel(0)   = gsn_csm_xy(wks, lhs1(:,:), hyb_plot, opts_mbgt)     ; Total LHS (local tendency)
   panel(1)   = gsn_csm_xy(wks, rhsT(:,:), hyb_plot, opts_mbgt)     ; Total RHS 
   panel(2)   = gsn_csm_xy(wks, rhs1(:,:), hyb_plot, opts_mbgt)     ; Mean radial advection
   panel(3)   = gsn_csm_xy(wks, rhs2(:,:), hyb_plot, opts_mbgt)     ; Mean vertical advection

   panel(4)   = gsn_csm_xy(wks, rhs3(:,:), hyb_plot, opts_mbgt)     ; Eddy radial avection 
   panel(5)   = gsn_csm_xy(wks, rhs4(:,:), hyb_plot, opts_mbgt)     ; Eddy vertical advection
   panel(6)   = gsn_csm_xy(wks, rhs5(:,:), hyb_plot, opts_mbgt)     ; Diffusive tendency (h + v)

 ;================================                                                           
 ; Overlay vertical line at x=0                                                              
 ;================================                                                           

 ; New graphical object to hold lines for each panel                                         
   lines = new( dimsizes(panel), graphic)

 ; Line resources                                                                            
 ; See https://www.ncl.ucar.edu/Document/Graphics/Images/dashpatterns.png for dash patterns  
   opts_l = True
   opts_l@gsLineThicknessF  = 2.0       ; Line thickness                                     
   opts_l@gsLineDashPattern = 1         ; Dash pattern                                       
   opts_l@gsLineColor       = "black"   ; Line colour                                        

 ; Overlay vertical line onto each panel (opts_mbgt@trYMaxF = 16.0)
   do n = 0, dimsizes(panel)-1
    lines(n) = gsn_add_polyline(wks,panel(n),(/0.0,0.0/),(/0.0,opts_mbgt@trYMaxF/),opts_l)
   end do

 ;=========================                                                                  
 ; Final panel resources                                                                     
 ;=========================                                                                  
   optsP                       = True
   optsP@gsnFrame              = False

 ; Turn on panel labelbar and customise size of final plot                               
   optsP@gsnPanelLabelBar      = False
   optsP@amJust                = "TopLeft"
   optsP@gsnPanelTop           = 0.98
   optsP@gsnPanelBottom        = 0.02
   optsP@gsnMaximize           = True

 ; Leave title blank                                                                     
   optsP@txString              = ""

 ; Panel figure labels                                                                   
   optsP@gsnPanelFigureStrings = (/"a) Local tendency (LHS)", "b) Source terms (RHS)", \
                                   "c) Mean_h", "d) Mean_v", "e) Eddy_h", \
                                   "f) Eddy_v", "h) Diff"/)

 ; Reduce size of labels (default 0.01)
   optsP@gsnPanelFigureStringsFontHeightF = 0.0099                                   

 ; Draw panel plot 
   optsP@gsnPanelRowSpec = True            ; Option 2 in call below is no. of plots on each row
   gsn_panel(wks, panel, (/4,3/), optsP)
   frame(wks)

end