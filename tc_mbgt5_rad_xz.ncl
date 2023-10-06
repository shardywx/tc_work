; Script to calculate the radial momentum budget as in Montgomery et al. (2018)

; Run using:

; ncl dat=\"02T12\" opt=\"pdf\" ens0=\"em11\" dist=3.0 ts0=815 tf0=839 clr=5
; calc=1 full=1 rad0=1 ave=0 wrt=0 ueh=0 new0=0 cn0=\"slp\" diff0=0 tc_mbgt5_rad_xz.ncl

; 'dat'   = initialisation time string ("02T12", "03T00", ...) 
; 'opt'   = output file format ("pdf" or "x11")
; 'ens0'  = ensemble member (em00 -- em11)
; 'dist'  = size of box following storm (degrees)
; 'ts0'   = starting time index (0–119)
; 'tf0'   = ending time index (0–119)
; 'clr'	  = colour map for plots (option 4 is used in Roger Smith's papers):
; 	    "ncl_default" [1],  "bdr_extra" [2], "amwg256" [3], "GMT_polar" [4]
; 'calc'  = centre calculated on each level (0) or equal to 950 hPa centre (1)
; 'full'  = plot 14 panels inc. u,v,w,mom (1), or 4 panels (0) 
; 'rad0'  = radius of final plots: 150 km (0) or 200 km (1)
; 'ave'	  = calculate 3/6/12-h average (1) rather than analyse individual time interval (0)
; 'wrt'	  = write out averages to text file (1) or leave alone (0 – default)
; 'ueh'	  = plot eddy horizontal advection of radial momentum (1) or leave out, as in M18 (0)
; 'new0'  = reduced colour scale for vtan (1) or original (0)
; 'diff0' = calculate 1h-averaged diffusivity/wind stress (1) or every 5 mins (0)

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/diff_avg.ncl"

begin

;==============================
; Prelash (reading in files)
;==============================

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
  u_cyc     = storm_rel[0]                 ; Zonal wind
  v_cyc     = storm_rel[1]                 ; Meridional wind
  vel_cyc   = storm_rel[2]                 ; Vector wind
  numTIMES  = storm_rel[3]                 ; Number of times in file
  lat_arr   = storm_rel[4]                 ; Storm-following domain subset (lat)
  lon_arr   = storm_rel[5]                 ; Storm-following domain subset (lon)
  centre    = storm_rel[6]                 ; Storm centre at each hour (lat/lon)

; Create arrays for data and output information
  llbox     = toint(dist*50)                   ; Calculate domain size (grid points)
  dsize     = (/llbox,llbox/)                  ; Domain size (grid points)
  d0        = dsize(0)
  d1        = dsize(1)
  ct        = ts0                              ; Counter variable
  nh        = (tf0 - ts0) / 12                 ; Number of hours analysed (for averaging later)
  num_l     = 63                               ; Number of model levels
  num_t     = numINPUT_k                       ; Number of times in 5-min data files

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

   if (hour(i) .lt. 10) then                    ; Add zeros in front of single-digit hours

    if (minute(i) .lt. 10) then
     hr_min(i) = "0"+hour(i)+"0"+mins(i)
    elseif (minute(i) .ge. 10) then
     hr_min(i) = "0"+hour(i)+mins(i)
    end if

   else

    if (minute(i) .lt. 10) then
     hr_min(i) = hour(i)+"0"+mins(i)            ; Also add zeros in front of single-digit minutes
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

;================================================================
; Calculate hourly averages from 5-min wind stress diagnostics
;================================================================

  if (diff0 .eq. 1) then 

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

  nr	     = 41
  ar	     = 2.0

  thetas     = new(73,float)
  pii        = 3.14159265
  radii      = fspan(0.0,ar,nr)      ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,73)   ; Array of azimuth angles (0, 5, 10, etc)
  thetas_rad = (pii/180)*thetas      ; Degrees to radians
  rad_size   = dimsizes(radii)-1     ; Index for later calculation

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

; Additional variables (same dimensions but not eddy quantities)
  rho0      = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [ro]

; Arrays for eddy derivative terms (4D)
; (1) du'/dr, (2) du'/dy, (3) du'/dz (4) dp'/dr
  ptb_dudr0 = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (1)
  ptb_dudy0 = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (2)
  ptb_dudz0 = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (3)
  ptb_dpdr0 = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (4)

; Arrays for eddy 1/r terms (4D)
; (1) v'/r, (2) v'^2 / r
  ptb_vt_r0  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (1)
  ptb_vttr0  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (2)

; Array for mean 1/r term (3D)
; (1) {v}^2 / r
  mean_vttr1 = new((/num_l,dimsizes(radii),numTIMES/),float) ; (1)

; Vertical derivative of mean radial wind [ d{u} / dz ]
  mean_dudz1          = new((/num_l,dimsizes(radii),numTIMES/),float)
  mean_dudz1!0        = "lev"
  mean_dudz1!1        = "rad"
  mean_dudz1!2        = "time"
  mean_dudz1@description = "Azimuthally averaged 'du_dz'"
  mean_dudz1@units = "s~S~-1~N~"

;===========================================================
; Arrays to hold expanded 'v_int'/'vort_int' arrays (4D)
;===========================================================

  u_int0    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Zonal wind
  v_int0    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Meridional wind
  w_int0    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Vertical velocity
  f_int0    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Coriolis parameter
  kmh_int0  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Horiz. diffusivity
  prs_int0  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Pressure
  rho_int0  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Density
  v_rad0    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Radial wind
  v_tan0    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Tangential wind

  difv_x0   = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Vert. diffusion (x)
  difv_y0   = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Vert. diffusion (y)
  difv_tot0 = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Vert. diff. (tot)

  difh_x0   = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Horiz. diffusion (x)
  difh_y0   = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)  ; Horiz. diffusion (y)

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

  prs_int0@description = "Pressure"
  prs_int0@units       = "hPa"
  prs_int0!0           = "lev"
  prs_int0!1           = "rad"
  prs_int0!2           = "azi"
  prs_int0!3           = "time"

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

  difv_y0@description  = "Vertical diffusion (y-component)"
  difv_y0@units        = "m s~S~-2~N~"
  difv_y0!0            = "lev"
  difv_y0!1            = "rad"
  difv_y0!2            = "azi"
  difv_y0!3            = "time"

  difv_x0@description  = "Vertical diffusion (x-component)"
  difv_x0@units        = "m s~S~-2~N~"
  difv_x0!0            = "lev"
  difv_x0!1            = "rad"
  difv_x0!2            = "azi"
  difv_x0!3            = "time"

  difh_x0@description  = "Horizontal diffusion (x-component)"
  difh_x0@units        = "m s~S~-2~N~"
  difh_x0!0            = "lev"
  difh_x0!1            = "rad"
  difh_x0!2            = "azi"
  difh_x0!3            = "time"

  difh_y0@description  = "Horizontal diffusion (y-component)"
  difh_y0@units        = "m s~S~-2~N~"
  difh_y0!0            = "lev"
  difh_y0!1            = "rad"
  difh_y0!2            = "azi"
  difh_y0!3            = "time"

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

;=======================================================================
; Arrays for LHS and RHS source terms containing eddy variables (4D)
;=======================================================================

; Mean horizontal advection of eddy radial momentum: LHS (3)
; {(u' * du'/dr) + (v'/r * du'/dy)}
; ptb_ueh1  = (ptb_vrad0 * ptb_dudr0) + ( (ptb_vtan0 / r) * (ptb_dudy0) )
  ptb_ueh1  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [L1]

; Eddy vertical advection of eddy radial momentum: RHS (2)
; [-w' * du'/dz]
; ptb_uev1  = - (ptb_vvel0 * ptb_dudz0) 
  ptb_uev1  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [R2]

; Eddy agradient force: RHS (3)
; [ {v'^2/r - 1/rho * dp'/dr} ]
; ptb_uegf1 = ( ( (ptb_vtan0) ^ 2) / r) - ( (1/rho) * ptb_dpdr0)
  ptb_uegf1 = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [R3]

; Arrays for LHS and RHS source terms after azimuthally averaging (4D)
  ptb_ueh0  = new((/num_l,dimsizes(radii),numTIMES/),float)		  ; [L1]
  ptb_uev0  = new((/num_l,dimsizes(radii),numTIMES/),float)  		  ; [R2]
  ptb_uegf0 = new((/num_l,dimsizes(radii),numTIMES/),float)	          ; [R3]

;==================================================
; Arrays for horizontal diffusive tendency terms 
;==================================================

; Horizontal diffusivity (K_mh)
  k_mh0	    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)

; Arrays for fractional expressions in diffusive tendency source terms below (4D)
; (1) du/dr, (2) dv/dy, (3) u/r, (4) r, (5) d/dr
  dudr0     = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (1)
  dvdy0	    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (2)
  u_r0      = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (3)
  r_0       = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (4)
  ddr0	    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (5)

; Additional arrays for terms outside main loop
  du        = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ;
  dv        = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ;
  dl        = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ;
  dr        = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ;
  r1        = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ;

  du0       = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; Radial 'u' increment
  du1       = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; Radial 'u'' increment
  dp1       = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; Radial 'p'' increment

  dmp       = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; Radial 'p' increment
  dmu       = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; Radial 'u' increment

; Radial stress tensor: RHS (5a)
; { (2 * K_mh) * (du/dr) }
; tau_rr1   = (2 * k_mh0) * (dudr0)
  tau_rr1   = new((/num_l,dimsizes(radii),numTIMES/),float)		     ; [3D array]

; Azimuthal stress tensor: RHS (5b)
; { (2 * K_mh) * ( (1/r * dv/dy) + u/r ) }
; tau_yy1   = (2 * k_mh0) * ( ( (1 / r_0) * dvdy0) + u_r0 ) 
  tau_yy1   = new((/num_l,dimsizes(radii),numTIMES/),float)		     ; [3D array]

; Total horizontal diffusive tendency of radial momentum 
; ( (1/r*{rho}) * d( r * {rho} * {tau_rr} )/dr) - tau_yy/r
; u_dh0	    = ((1/(r_0*mean_rho0)) * (ddr0*(r_0*mean_rho0*tau_rr1)))-((1/r_0)*(tau_yy1))
  u_dh0     = new((/num_l,dimsizes(radii),numTIMES/),float)		     ; [3D array]

; Total vertical diffusive tendency of radial momentum 
; u_dz0	    = (1 /{rho0}) * ( (d/dz) * ({rho0} * {tau_rz}) )  
  u_dz0	    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [4D array]

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

; Mean pressure {p}
  mean_pres0       = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  mean_pres0!0     = "lev"
  mean_pres0!1     = "rad"
  mean_pres0!2     = "azi"
  mean_pres0!3     = "time"
  mean_pres0@description = "Azimuthally averaged pressure"
  mean_pres0@units = "hPa"

; Mean density {rho}
  mean_rho0        = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  mean_rho0!0      = "lev"
  mean_rho0!1      = "rad"
  mean_rho0!2      = "azi"
  mean_rho0!3      = "time"
  mean_rho0@description = "Azimuthally averaged density"
  mean_rho0@units  = "kg m~S~-3~N~"

; Mean Coriolis force {f0}
  mean_cor0        = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  mean_cor0!0      = "lev"
  mean_cor0!1      = "rad"
  mean_cor0!2      = "azi"
  mean_cor0!3      = "time"
  mean_cor0@description = "Azimuthally averaged Coriolis force"
  mean_cor0@units  = "s~S~-1~N~"

; Mean vertical diffusion (total)
  mean_difv0          = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  mean_difv0!0        = "lev"
  mean_difv0!1        = "rad"
  mean_difv0!2        = "azi"
  mean_difv0!3        = "time"
  mean_difv0@description = "Azimuthally averaged vertical diffusion (tot)"
  mean_difv0@units    = "m s~S~-1~N~"

; Radial derivative of mean radial wind [ d{u} / dr ]
  mean_dudr0          = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  mean_dudr0!0        = "lev"
  mean_dudr0!1        = "rad"
  mean_dudr0!2        = "azi"
  mean_dudr0!3        = "time"
  mean_dudr0@description = "Azimuthally averaged 'du_dr'"

; Radial derivative of pressure [ d{p} / dr ]
  mean_dpdr0          = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  mean_dpdr0!0        = "lev"
  mean_dpdr0!1        = "rad"
  mean_dpdr0!2        = "azi"
  mean_dpdr0!3        = "time"
  mean_dpdr0@description = "Azimuthally averaged 'dp_dr'"

;==================================
; Get the variables we will need
;==================================

   c  = addfile(fili_c,"r")              ; Read in 'pc' stream [u,v,vort]
   d  = addfile(fili_d,"r")              ; Read in 'pd' stream [w,p]
   j  = addfile(fili_j,"r")              ; Read in 'pj' stream [t,theta,mixing ratio]
   k  = addfile(fili_k,"r")              ; Read in 'pk' stream [tau,PV,w]

   do it = ts0, tf0                      ; Loop over times in file

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
    hybC   = c->hybrid_ht(:)                      ; 63 model ('theta') levels

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

;======================================
; Read in variables from 'pd' stream
;======================================

  ; Read in basic variables from 'pd' stream
    hybD   = d->hybrid_ht(:)                        ; 63 model ('rho') levels

  ; Read in temperature, pressure, vertical velocity and geopotential height
  ; All variables below -- [hybD | 63] * [latD | 150] * [lon1D | 151]
    p      = d->p(it,:,{t0:t1},{n0:n1})             ; Pressure (Pa)
    z      = d->ht(it,:,{t0:t1},{n0:n1})            ; Geopotential height (m)

;======================================
; Read in variables from 'pj' stream
;======================================

  ; Read in basic variables from 'pj' stream
    hybJ   = j->hybrid_ht(:)                        ; 63 model ('rho') levels

  ; Read in temperature, potential temperature and mixing ratio
  ; All variables below  -- [hybJ | 63] * [latJ | 150]  * [lonJ | 150]
    t      = j->temp(it,:,{t0:t1},{n0:n1})          ; Temperature
    th     = j->theta(it,:,{t0:t1},{n0:n1})         ; Potential temperature
    mix    = j->mix(it,:,{t0:t1},{n0:n1})           ; Mixing ratio

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
    pv     = k->field83(it,0:62,{t0:t1},{n0:n1})    ; Potential vorticity

  ; Read in horizontal wind stress components and diffusivity                     
    if (diff0 .eq. 0) then
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

  ; Convert units of pressure and geopotential height before continuing
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
     if (diff0 .eq. 1) then 

     ; Wind stress components on current model level
       tau_xz = taux_arr(it,ilev,0:d0-1,0:d1-1)
       tau_yz = tauy_arr(it,ilev,0:d0-1,0:d1-1)

     ; Wind stress components on model levels above/below current level
       tau_x1 = taux_arr(it,ip1,0:d0-1,0:d1-1)
       tau_x0 = taux_arr(it,im1,0:d0-1,0:d1-1)
       tau_y1 = tauy_arr(it,ip1,0:d0-1,0:d1-1)
       tau_y0 = tauy_arr(it,im1,0:d0-1,0:d1-1)
       dz     = hybD(ip1) - hybD(im1)

     elseif (diff0 .eq. 0) then 

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
     rho_plane              = rho(ilev,:,:)	   ; Density
     prs_plane              = p(ilev,:,:)	   ; Pressure
     vort_plane	   	    = vort(ilev,:,:)	   ; Relative vorticity

   ; Read in horizontal diffusivity as normal (diff=0) or from pre-calculated array (diff=1)  
     if (diff0 .eq. 1) then
      kmh_plane             = kmh_arr(it,ilev,:,:)
     elseif (diff0 .eq. 0) then
      kmh_plane             = kmh(ilev,:,:)
     end if

   ; Calculate location of vorticity maximum (storm centre) on model level
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

   ; Difference between 'lat_max'/'lon_max' and values read in using NCL function
     dlat 	   	    = abs(lat_max - centre(cn,it-1,0))
     dlon	   	    = abs(lon_max - centre(cn,it-1,1))

   ; If our calculations of latitude and longitude deviate too much from the values
   ; we read in using the NCL function ('st_rm'), then calculate a new centre
   ; using the minimum pressure. Constrain the new calculation using the values
   ; from 'st_rm', to prevent spurious centres being identified and analysed. 

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

      lat_max 	 = lt1 + (dy0 * inds_p0(0,0))
      lon_max 	 = ln1 + (dx0 * inds_p0(0,1))
      print("New centre: "+lat_max+" degrees N, "+lon_max+" degrees E")
      delete([/prs_plane0,prs_min0,dims_p0,prs_1d0,inds_p0/])      

     end if 

   elseif (calc .eq. 1) then

      lat_max    = centre(cn,it-1,0)
      lon_max    = centre(cn,it-1,1)

   elseif (calc .eq. 2) then 

      lat_max    = t0 + (dy0 * inds0(0,0))
      lon_max    = n0 + (dx0 * inds0(0,1))
      print("Centre (980 m): "+lat_max+" degrees N, "+lon_max+" degrees E")

   end if

   ;======================================================================
   ; Create lat/lon arrays for later switch to cylindrical coordinates
   ;======================================================================

   ; Indices of 'lat'/'lon' grids
     lat_plane   = ind(t0.le.lat.and.lat.le.t1)
     lon_plane   = ind(n0.le.lon.and.lon.le.n1)

   ; Array sizes
     lat_size    = dimsizes(lat_plane)
     lon_size    = dimsizes(lon_plane)

   ; Array of latitude/longitude points corresponding to the indices above
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
     prs_int  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Pressure 
     rho_int  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Density
     v_rad    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Radial velocity
     v_tan    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Tangential velocity

   ; Additional arrays for vertical diffusion (calculate horizontal diffusion later)
     difv_x   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vert. diff (x)
     difv_y   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vert. diff (y)
     difv_tot = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vert. diff (total)

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
      f_int(:,iang)    = linint2_points(lon_sub,lat_sub,f_plane,False,\
				          xpos(:,iang),ypos(:,iang),0)
      difv_x(:,iang)   = linint2_points(lon_sub,lat_sub,difv_xplane,False,\
					  xpos(:,iang),ypos(:,iang),0)
      difv_y(:,iang)   = linint2_points(lon1_sub,lat1_sub,difv_yplane,False,\
					  xpos(:,iang),ypos(:,iang),0)
      w_int(:,iang)    = linint2_points(lon1_sub,lat_sub,w_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      prs_int(:,iang)  = linint2_points(lon1_sub,lat_sub,prs_plane,False,\
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

      ; Total vertical diffusive tendency of radial momentum
        difv_tot(:,iang) = difv_y(:,iang)*sin(theta) + difv_x(:,iang)*cos(theta)

      end do

   ; Fill 'u_dz0' array with values from 'difv_tot' (vertical diffusive tendency)
     u_dz0(ilev,:,:,ct)   = difv_tot(:,:)    

    ;===================================================
    ; Calculate mean (azimuthally-averaged) variables 
    ;===================================================

      do irad = 0, dimsizes(radii)-1

        mean_vtan0(ilev,irad,:,ct)    = avg(v_tan(irad,:))      ; [ {v} ]
        mean_vrad0(ilev,irad,:,ct)    = avg(v_rad(irad,:))      ; [ {u} ]
        mean_vvel0(ilev,irad,:,ct)    = avg(w_int(irad,:))      ; [ {w} ]
        mean_pres0(ilev,irad,:,ct)    = avg(prs_int(irad,:))    ; [ {p} ]
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
      prs_int0(ilev,:,:,ct)  = (/prs_int(:,:)/)
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
      delete([/rho_int,difv_x,difv_y,f_int,kmh_int,prs_int/])
      delete([/u_plane,v_plane,w_plane,f_plane,vort_plane,prs_plane/])
      delete([/difv_xplane,difv_yplane,kmh_plane,xpos,ypos,xcen,ycen/])
      delete([/lat_plane,lon_plane,lon1_plane,lat1_plane/])
      delete([/lat_sub,lon_sub,lat1_sub,lon1_sub/])

     end do			    ; End pressure level loop (do ilev = 0, num_l-1)

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
      ptb_pres0(:,:,:,:) = prs_int0(:,:,:,:)  - mean_pres0(:,:,:,:)    ; [p']

      k_mh0(:,:,:,:)     = kmh_int0(:,:,:,:)                           ; [kmh]
      rho0(:,:,:,:)      = rho_int0(:,:,:,:)                           ; [rho]

;=================================================================
; Calculate eddy terms involving azimuthal derivatives [du'/dy]
;=================================================================

    ; Calculate the azimuthal distance ('dy') between all points in array
    ; Arguments 1 and 2 represent the lat/lon of the first point, respectively
    ; Arguments 3 and 4 represent the lat/lon of the second point, respectively
    ; Argument 5 ('2') tells the script to avoid extra interpolated points
    ; Argument 6 ('3') returns the distance in metres

      dth  = dimsizes(thetas)
      dtt  = numTIMES

      dl(:,:,0:dth-2,:)  = gc_latlon( ypos_all(:,:,0:dth-2,:),xpos_all(:,:,0:dth-2,:),\
                                      ypos_all(:,:,1:dth-1,:),xpos_all(:,:,1:dth-1,:), 2, 3)

    ; Replace zeros with positive values
      dl  = where( dl .ne. 0, dl, 2000)

    ; Calculate increments
      du(:,:,0:dth-2,:)    = ptb_vrad0(:,:,1:dth-1,:) - ptb_vrad0(:,:,0:dth-2,:)
      dv(:,:,:,0:dtt-2)    = v_tan0(:,:,:,1:dtt-1) - v_tan0(:,:,:,0:dtt-2)

    ; Calculate [du'/dl]
      ptb_dudy0(:,:,:,:)   = du / dl

    ; Calculate [dv/dy -- for horizontal diffusive tendency]
      dvdy0(:,:,:,:)       = dv / dl

;======================================================================
; Calculate eddy terms involving radial derivatives [du'/dr, dp'/dr]
;======================================================================

      ptb_pres0 = ptb_pres0 * 100

    ; Calculate the radial distance ('dr') between all points in array
    ; Arguments 1 and 2 represent the lat/lon of the first point, respectively
    ; Arguments 3 and 4 represent the lat/lon of the second point, respectively
    ; Argument 5 ('2') tells the script to avoid extra interpolated points
    ; Argument 6 ('3') returns the distance in metres

      drad = dimsizes(radii)

      dr(:,0:drad-2,:,:)  = gc_latlon( ypos_all(:,0:drad-2,:,:),xpos_all(:,0:drad-2,:,:),\
                                       ypos_all(:,1:drad-1,:,:),xpos_all(:,1:drad-1,:,:), 2, 3)

    ; Replace zeros with positive values
      dr = where( dr .ne. 0, dr, 10000)

      r1(:,0:drad-2,:,:)  = gc_latlon( ycen_all(:,1:drad-1,:,:),xcen_all(:,1:drad-1,:,:),\
                                       ypos_all(:,1:drad-1,:,:),xpos_all(:,1:drad-1,:,:), 2, 3)

    ; Replace zeros with positive values
      r1 = where( r1 .ne. 0, r1, 10000)

    ; Calculate increments ('u', 'u_dash' and 'p_dash') in radial direction
      du0(:,0:drad-2,:,:)   = v_rad0(:,1:drad-1,:,:) - v_rad0(:,0:drad-2,:,:)
      du1(:,0:drad-2,:,:)   = ptb_vrad0(:,1:drad-1,:,:) - ptb_vrad0(:,0:drad-2,:,:)

      dp1(:,0:drad-2,:,:)   = ptb_pres0(:,1:drad-1,:,:) - ptb_pres0(:,0:drad-2,:,:)

      vp   = ptb_vtan0(:,:,:,:)
      vpp  = vp ^ 2

    ; Calculate [du'/dr]
      ptb_dudr0(:,:,:,:)    = du1 / dr

    ; Calculate [du/dr -- for horizontal diffusive tendency]
      dudr0(:,:,:,:)        = du0 / dr

    ; Calculate [dp'/dr]
      ptb_dpdr0(:,:,:,:)    = dp1 / dr

    ; Calculate [u/r -- for horizontal diffusive tendency]
      u_r0(:,:,:,:)         = v_rad0(:,:,:,:) / r1

    ; Calculate [r -- for horizontal diffusive tendency]
      r_0(:,:,:,:)          = r1

    ; Calculate [d/dr -- for horizontal diffusive tendency]
      ddr0(:,:,:,:)         = 1 / dr

    ; Calculate [v'/r]
      ptb_vt_r0(:,:,:,:)    = vp / r1

    ; Calculate [v'^2/r]
      ptb_vttr0(:,:,:,:)    = vpp / r1

      ptb_pres0 = ptb_pres0 / 100

;=========================================================================
; Calculate mean terms involving radial derivatives [ d{u}/dr, d{p}/dr ]
;=========================================================================

    ; Calculate the radial distance ('dr') between all points in array
    ; Arguments 1 and 2 represent the lat/lon of the first point, respectively
    ; Arguments 3 and 4 represent the lat/lon of the second point, respectively
    ; Argument 5 ('2') tells the script to avoid extra interpolated points
    ; Argument 6 ('3') returns the distance in metres

    ; 'r1' and 'dr' have already been defined

;=============================================
; Calculate absolute angular momentum (AAM)
;=============================================

    ; Temporarily change units of pressure (to Pa)
      mean_pres0 = mean_pres0 * 100

    ; Calculate increments ('u' and 'p') in radial direction
      dmu(:,0:drad-2,:,:)  = mean_vrad0(:,1:drad-1,:,:) - mean_vrad0(:,0:drad-2,:,:)
      dmp(:,0:drad-2,:,:)  = mean_pres0(:,1:drad-1,:,:) - mean_pres0(:,0:drad-2,:,:)

      vm   = mean_vtan0(:,:,:,:)
      vmm  = vm ^ 2
      fm   = mean_cor0(:,:,:,:) * (r1 ^ 2)

    ; Change units of pressure back (to hPa)
      mean_pres0 = mean_pres0 / 100

    ; Calculate absolute angular momentum components
      aam0 = vm * r1
      aam1 = 0.5 * fm

    ; Finally, calculate absolute angular momentum [ rv + 0.5 * (fr ^ 2) ]
      mean_aam0(:,:,:,:)  = aam0 + aam1

;====================================================
; Calculate additional components including radius
;====================================================

    ; Calculate [ d{u}/dr ]
      mean_dudr0          = dmu / dr

    ; Calculate [ d{p}/dr ]
      mean_dpdr0          = dmp / dr

    ; Calculate [ {v}^2/r ]
      mean_vttr0          = vmm / r1

;=======================================================
; Remove unncessary 3rd dimension from mean variables
;=======================================================

      mean_vtan1 = mean_vtan0(:,:,0,:)
      mean_vrad1 = mean_vrad0(:,:,0,:)
      mean_vvel1 = mean_vvel0(:,:,0,:)
      mean_rho1  = mean_rho0(:,:,0,:)
      mean_aam1  = mean_aam0(:,:,0,:)
      mean_cor1  = mean_cor0(:,:,0,:)

      mean_dudr1 = mean_dudr0(:,:,0,:)
      mean_dpdr1 = mean_dpdr0(:,:,0,:)
      mean_vttr1 = mean_vttr0(:,:,0,:)

;================================================
; Calculate vertical derivatives of mean terms 
;================================================

    ; Redefine model level array for vertical derivative calculations
      mlevs      = dimsizes(hybC)      ; Size of model level array
      hybC@units = "m"		       ; Model level units 

      do it = ts0, tf0

	do ilev = 0, num_l-1		; Loop over all model levels
	
        mp1 = min((/ilev+1,num_l-1/))
        mm1 = max((/ilev-1,0/))
 
          do irad = 0, dimsizes(radii)-1

          ;=====================================================
	  ; Calculate vertical derivative of mean radial wind
	  ;=====================================================

          ; [ d{u}/dz ]
            mean_dudz1(ilev,irad,it)    = (mean_vrad1(mp1,irad,it) - \
                                           mean_vrad1(mm1,irad,it)) / \
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

           ;=====================================================
           ; Calculate vertical derivative of eddy radial wind
           ;=====================================================

           ; [ du'/dz ]
             ptb_dudz0(ilev,irad,th0,it) = (ptb_vrad0(mp1,irad,th0,it) - \
                                            ptb_vrad0(mm1,irad,th0,it)) / \
                                           (hybC(mp1) - hybC(mm1))

           end do	         ; End loop over azimuth angles (do th0 = 0, ...)

          end do                 ; End loop over radial circles (do irad = 0, ...)

        end do                   ; End loop over model levels (do ilev = 0, num_l-1)

       end do                    ; End loop over times (do it = ts0, tf0)

;==========================================================
; Finally, calculate source terms and azimuthal averages 
;==========================================================
	
       print_clock("Calculating source terms and azimuthal averages...")

     ;===================================================================
     ; Mean radial advection of radial momentum per unit mass: LHS (2)
     ;===================================================================
       lhs2_st = mean_vrad1 * mean_dudr1

     ;=============================================================
     ; Mean horizontal advection of eddy radial momentum: LHS (3)
     ;=============================================================
       ptb_ueh1 = (ptb_vrad0 * ptb_dudr0) + (ptb_vt_r0 * ptb_dudy0)

     ; Calculate the azimuthal average of 'ptb_ueh1' and input into 'lhs3_st
       lhs3_st = dim_avg_n_Wrap(ptb_ueh1(:,:,:,:),2)

     ;=========================================================================
     ; Mean vertical advection of mean radial momentum per unit mass: RHS (1)
     ;=========================================================================
       rhs1_st = -(mean_vvel1 * mean_dudz1)

     ;=========================================================================
     ; Eddy vertical advection of eddy radial momentum per unit mass: RHS (2)
     ;=========================================================================
       ptb_uev1 = -(ptb_vvel0 * ptb_dudz0)

     ; Calculate the azimuthal average of 'ptb_uev1' and input into 'rhs2_st'
       rhs2_st  = dim_avg_n_Wrap(ptb_uev1(:,:,:,:),2)

     ;==============================================
     ; Mean agradient force per unit mass: RHS (3)
     ;==============================================

     ; Account for elements where 'mean_rho0' is equal to zero (should be none)
       mean_rho1 = where(mean_rho1.ne.0, mean_rho1, mean_rho1@_FillValue)

       rhs3_st = mean_vttr1 + (mean_cor1 * mean_vtan1) - ( (1 / mean_rho1) * mean_dpdr1 )

     ;==============================================
     ; Eddy agradient force per unit mass: RHS (4)
     ;==============================================

     ; Account for elements where 'rho0' is equal to zero
       rho0  = where(rho0 .ne. 0, rho0, rho0@_FillValue)
     
     ; Now calculate eddy agradient force 
       ptb_uegf1 = ptb_vttr0 - ( (1 / rho0) * ptb_dpdr0 )

       printMinMax(ptb_uegf1, False)

       printMinMax(ptb_vttr0, False)
       printMinMax(ptb_dpdr0, False)
       exit()

     ; CONTINUE EDITING HERE 

     ; Calculate the azimuthal average of 'ptb_uegf1' and input into 'rhs4_st'
       rhs4_st   = dim_avg_n_Wrap(ptb_uegf1(:,:,:,:),2)

     ;===================================================================================
     ; Combined mean radial and vertical diffusive tendency of radial momentum: RHS (5)
     ;===================================================================================

     ; Stress tensor (radial component) 
       tau_rr1   = dim_avg_n_Wrap( ((2 * k_mh0) * (dudr0)), 2)
     
     ; Stress tensor (azimuthal component)
       tau_yy1	 = dim_avg_n_Wrap( ((2 * k_mh0) * ( ( (1 / r_0) * dvdy0) + u_r0 )), 2)

     ;=============================================================================
     ; Use stress tensors to calculate horizontal diffusive tendency source term
     ;=============================================================================

     ; Azimuthally average the 'r' and 'd/dr' expressions
       r_1   = dim_avg_n_Wrap(r_0(:,:,:,:),2)
       ddr1  = dim_avg_n_Wrap(ddr0(:,:,:,:),2)

       u_dh0 = ( ( 1/(r_1 * mean_rho1) ) * (ddr1 * (r_1 * mean_rho1 * tau_rr1)) ) - \
               ( (1/r_1) * (tau_yy1) )

       rhs5a_st = u_dh0	       

     ; Vertical diffusive tendency source term has already been calculated (L827)
     ; Now calculate the azimuthal average using 'dim_avg_n_Wrap'
       rhs5b_st = dim_avg_n_Wrap( u_dz0(:,:,:,:),2 )

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
       lhs1A    = new((/num_l,dimsizes(radii),numTIMES/),float)   ; Local vrad tendency
       lhs2A    = new((/num_l,dimsizes(radii),numTIMES/),float)   ; Mean r-adv of rad mom
       lhs3A    = new((/num_l,dimsizes(radii),numTIMES/),float)   ; Mean h-adv of eddy mom
       lhsTA    = new((/num_l,dimsizes(radii),numTIMES/),float)   ; Sum of LHS terms

       rhs1A    = new((/num_l,dimsizes(radii),numTIMES/),float)   ; Mean v-adv of mean mom
       rhs2A    = new((/num_l,dimsizes(radii),numTIMES/),float)   ; Eddy v-adv of eddy mom
       rhs3A    = new((/num_l,dimsizes(radii),numTIMES/),float)   ; Mean agradient force 
       rhs4A    = new((/num_l,dimsizes(radii),numTIMES/),float)   ; Eddy agradient force
       rhs5aA   = new((/num_l,dimsizes(radii),numTIMES/),float)   ; Mean horiz diffusion
       rhs5bA   = new((/num_l,dimsizes(radii),numTIMES/),float)   ; Mean vert diffusion
       rhsTA    = new((/num_l,dimsizes(radii),numTIMES/),float)   ; Sum of RHS terms

     ; Define additional arrays if calculating averages 

       if (ave .eq. 1) then 

         lhs1_av0  = new((/num_l,dimsizes(radii)/),float)
         lhs2_av0  = new((/num_l,dimsizes(radii)/),float)
         lhs3_av0  = new((/num_l,dimsizes(radii)/),float)
         lhsT_av0  = new((/num_l,dimsizes(radii)/),float)
 
         rhs1_av0  = new((/num_l,dimsizes(radii)/),float)
         rhs2_av0  = new((/num_l,dimsizes(radii)/),float)
         rhs3_av0  = new((/num_l,dimsizes(radii)/),float)
         rhs4_av0  = new((/num_l,dimsizes(radii)/),float)
         rhs5a_av0 = new((/num_l,dimsizes(radii)/),float)
         rhs5b_av0 = new((/num_l,dimsizes(radii)/),float)
         rhsT_av0  = new((/num_l,dimsizes(radii)/),float)

       end if 

     ;=================================================
     ; Calculate the mean radial wind tendency (LHS)
     ;=================================================

       if (ave .eq. 1) then 

       ; 'dt' calculated in hours
         dt0 = tofloat( (tf0 - ts0) / 12)

	 lhs1_av0(:,:) = (mean_vrad1(:,:,tf0) - mean_vrad1(:,:,ts0)) / dt0

       elseif (ave .eq. 0) then 

         do it = ts0, tf0

           cp1 = min((/it+1,tf0/))
  	   cm1 = max((/it-1,ts0/))
 
           lhs1A(:,:,it) = (mean_vrad1(:,:,cp1) - mean_vrad1(:,:,cm1)) * 6     ; Units [m s-1 h-1]

         end do	           ; End loop over times (do it = ts0, tf0)

       end if 

     ; Add metadata and make units consistent for all terms
       lhs1A@units = "m s~S~-1~N~ h~S~-1~N~"

       lhs2A = (/lhs2_st/)
       lhs2A = lhs2A * 3600
       lhs2A@units = "m s~S~-1~N~ h~S~-1~N~"

       lhs3A = (/lhs3_st/)
       lhs3A = lhs3A * 3600
       lhs3A@units = "m s~S~-1~N~ h~S~-1~N~"

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

     ; Add metadata for the time-averaged radial wind tendency
       if (ave .eq. 1) then
        lhs1_av0!0    = "lev"
        lhs1_av0&lev  = hyb_plot
        lhs1_av0!1    = "rad"
        lhs1_av0&rad  = radii
       end if 

     ; Add metadata for all the other important arrays (instantaneous)
       lhs1A!0    = "lev"
       lhs1A&lev  = hyb_plot
       lhs1A!1    = "rad"
       lhs1A&rad  = radii
       lhs1A!2    = "time"

       lhs2A!0    = "lev"
       lhs2A&lev  = hyb_plot
       lhs2A!1    = "rad"
       lhs2A&rad  = radii
       lhs2A!2    = "time"

       lhs3A!0    = "lev"
       lhs3A&lev  = hyb_plot
       lhs3A!1    = "rad"
       lhs3A&rad  = radii
       lhs3A!2    = "time"

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

       lhsTA!0    = "lev"
       lhsTA&lev  = hyb_plot
       lhsTA!1    = "rad"
       lhsTA&rad  = radii
       lhsTA!2    = "time"

     ;==============================================
     ; Calculate time averages of all terms above
     ;==============================================

       if (ave .eq. 1) then 

        lhs2_av0   = dim_avg_n_Wrap(lhs2A(:,:,ts0:tf0),2)	  ; Mean r-adv of rad mom
        lhs3_av0   = dim_avg_n_Wrap(lhs3A(:,:,ts0:tf0),2)	  ; Mean h-adv of eddy mom
        lhsT_av0   = dim_avg_n_Wrap(lhsTA(:,:,ts0:tf0),2)	  ; Sum of LHS terms
 
        rhs1_av0   = dim_avg_n_Wrap(rhs1A(:,:,ts0:tf0),2)	  ; Mean v-adv of mean mom
        rhs2_av0   = dim_avg_n_Wrap(rhs2A(:,:,ts0:tf0),2)	  ; Eddy v-adv of eddy mom
        rhs3_av0   = dim_avg_n_Wrap(rhs3A(:,:,ts0:tf0),2)	  ; Mean agradient force
        rhs4_av0   = dim_avg_n_Wrap(rhs4A(:,:,ts0:tf0),2)	  ; Eddy agradient force
        rhs5a_av0  = dim_avg_n_Wrap(rhs5aA(:,:,ts0:tf0),2)	  ; Mean horiz diffusion
        rhs5b_av0  = dim_avg_n_Wrap(rhs5bA(:,:,ts0:tf0),2)	  ; Mean vert diffusion
        rhsT_av0   = dim_avg_n_Wrap(rhsTA(:,:,ts0:tf0),2)	  ; Sum of RHS terms

       end if 

;==============================================
; Switch order of dimensions before plotting
;==============================================

     ; 'gsn_csm_contour' plots the leftmost (rightmost) dimension on the y-axis (x-axis)
     ; Switch 'rad' with 'time' so that 'rad' is the rightmost dimension

     ; Mean radial wind -- also add metadata
       mean_vrad       = mean_vrad1(lev|:,time|:,rad|:)
       mean_vrad!0     = "lev"
       mean_vrad&lev   = hyb_plot
       mean_vrad!1     = "time"
       mean_vrad!2     = "rad"
       mean_vrad&rad   = radii
       mean_vrad@units = "m s~S~-1~N~"

     ; Mean tangential wind -- also add metadata
       mean_vtan       = mean_vtan1(lev|:,time|:,rad|:)
       mean_vtan!0     = "lev"
       mean_vtan&lev   = hyb_plot
       mean_vtan!1     = "time"
       mean_vtan!2     = "rad"
       mean_vtan&rad   = radii
       mean_vtan@units = "m s~S~-1~N~"

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
       mean_aam	       = mean_aam / (10 ^ 6)
       mean_aam@units  = "10~S~6~N~ m~S~2~N~ s~S~-1~N~"

     ;==================================================
     ; Also calculate time averages of mean u,v,w,aam
     ;==================================================

       if (ave .eq. 1) then

        vrad_av0   = dim_avg_n_Wrap(mean_vrad(:,ts0:tf0,:),1)     ; Mean r-adv of rad mom
        vtan_av0   = dim_avg_n_Wrap(mean_vtan(:,ts0:tf0,:),1)     ; Mean h-adv of eddy mom
        vvel_av0   = dim_avg_n_Wrap(mean_vvel(:,ts0:tf0,:),1)     ; Sum of LHS terms
        aam_av0    = dim_avg_n_Wrap(mean_aam(:,ts0:tf0,:),1)      ; Mean v-adv of mean mom

       end if

     ;==================================
     ; Reorder arrays before plotting
     ;==================================
       
       lhs1  	  = lhs1A(lev|:,time|:,rad|:)
       lhs2  	  = lhs2A(lev|:,time|:,rad|:)
       lhs3  	  = lhs3A(lev|:,time|:,rad|:)

       rhs1  	  = rhs1A(lev|:,time|:,rad|:)
       rhs2  	  = rhs2A(lev|:,time|:,rad|:)
       rhs3  	  = rhs3A(lev|:,time|:,rad|:)
       rhs4  	  = rhs4A(lev|:,time|:,rad|:)
       rhs5a 	  = rhs5aA(lev|:,time|:,rad|:)
       rhs5b 	  = rhs5bA(lev|:,time|:,rad|:)

       if (ave .eq. 1) then 

         lhs1_av  = lhs1_av0(lev|:,rad|:)
	 lhs2_av  = lhs2_av0(lev|:,rad|:)
	 lhs3_av  = lhs3_av0(lev|:,rad|:)

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
       lhs2  	  = smth9_Wrap(lhs2, 0.5, 0.5, True)
       lhs3  	  = smth9_Wrap(lhs3, 0.5, 0.5, True)

       rhs1  	  = smth9_Wrap(rhs1, 0.5, 0.5, True)
       rhs2  	  = smth9_Wrap(rhs2, 0.5, 0.5, True)
       rhs3  	  = smth9_Wrap(rhs3, 0.5, 0.5, True)
       rhs4  	  = smth9_Wrap(rhs4, 0.5, 0.5, True)
       rhs5a 	  = smth9_Wrap(rhs5a, 0.5, 0.5, True)
       rhs5b 	  = smth9_Wrap(rhs5b, 0.5, 0.5, True)

     ; Calculate total RHS tendency
       rhsT       = rhs1 + rhs2 + rhs3 + rhs4 + rhs5a + rhs5b
       rhsT!0     = "lev"
       rhsT&lev   = hyb_plot
       rhsT!1     = "time"
       rhsT!2     = "rad"
       rhsT&rad   = radii
       rhsT@units = "m s~S~-1~N~ h~S~-1~N~"

     ; Calculate total LHS tendency
       lhsT       = lhs1 + lhs2 + lhs3
       lhsT!0     = "lev"
       lhsT&lev   = hyb_plot
       lhsT!1     = "time"
       lhsT!2     = "rad"
       lhsT&rad   = radii
       lhsT@units = "m s~S~-1~N~ h~S~-1~N~"

     ; Calculate sum of all mean terms (vadv + hadv + diff)
       mean       = lhs2 + rhs1 + rhs5a + rhs5b
       mean!0     = "lev"
       mean&lev   = hyb_plot
       mean!1     = "time"
       mean!2     = "rad"
       mean&rad   = radii
       mean@units = "m s~S~-1~N~ h~S~-1~N~"

     ; Calculate sum of all eddy terms (vadv + hadv)
       eddy       = rhs2 + lhs3
       eddy!0     = "lev"
       eddy&lev   = hyb_plot
       eddy!1     = "time"
       eddy!2     = "rad"
       eddy&rad   = radii
       eddy@units = "m s~S~-1~N~ h~S~-1~N~"

     ; Calculate sum of both diffusion terms
       diff       = rhs5a + rhs5b
       diff!0     = "lev"
       diff&lev   = hyb_plot
       diff!1     = "time"
       diff!2     = "rad"
       diff&rad   = radii
       diff@units = "m s~S~-1~N~ h~S~-1~N~"

       if (ave .eq. 1) then 

      ; Total LHS        
        lhsT_av       = lhs1_av + lhs2_av + lhs3_av
        lhsT_av!0     = "lev"
        lhsT_av&lev   = hyb_plot
        lhsT_av!1     = "rad"
        lhsT_av&rad   = radii
        lhsT_av@units = "m s~S~-1~N~ h~S~-1~N~"

      ; Total RHS 
	rhsT_av	      = rhs1_av + rhs2_av + rhs3_av + rhs4_av + rhs5a_av + rhs5b_av
        rhsT_av!0     = "lev"
        rhsT_av&lev   = hyb_plot
        rhsT_av!1     = "rad"
        rhsT_av&rad   = radii
        rhsT_av@units = "m s~S~-1~N~ h~S~-1~N~"

      ; All mean terms 
        mean_av       = lhs2_av + rhs1_av + rhs5a_av + rhs5b_av
        mean_av!0     = "lev"
        mean_av&lev   = hyb_plot
        mean_av!1     = "rad"
        mean_av&rad   = radii
        mean_av@units = "m s~S~-1~N~ h~S~-1~N~"
 
      ; All eddy terms 
        eddy_av       = rhs2_av + lhs3_av
        eddy_av!0     = "lev"
        eddy_av&lev   = hyb_plot
        eddy_av!1     = "rad"
        eddy_av&rad   = radii
        eddy_av@units = "m s~S~-1~N~ h~S~-1~N~"

      ; All diffusion terms
        diff_av       = rhs5a_av + rhs5b_av
	diff_av!0     = "lev"
        diff_av&lev   = hyb_plot
        diff_av!1     = "rad"
        diff_av&rad   = radii
        diff_av@units = "m s~S~-1~N~ h~S~-1~N~"

       end if 

;==================================================================
; Before plotting, write out arrays to text file (write_matrix)
;==================================================================

;       if (wrt .eq. 1) then 

;      ; Title resources for 'write_matrix' procedure
;        opts_mat       = True 
;	opts_mat@row   = True
;	opts_mat@title = ""
;	opt_arr        = (/"lhs1", "lhs2", "lhs3", "rhs1", "rhs2", \
;		           "rhs3", "rhs4", "rhs5a", "rhs5b", \	   
;			   "vrad", "vtan", "vvel", "aam"/)
;        mbgt_arr       = (/lhs1_av, lhs2_av, lhs3_av, rhs1_av, rhs2_av, \
;		       	   rhs3_av, rhs4_av, rhs5a_av, rhs5b_av, vrad_av0, \
;			   vtan_av0, vvel_av0, aam_av0/)

;      ; Fill composite array with values from terms in 'mbgt_arr'
;        do mb = 0, ndiag-1
;         plot_arr(:,st,mb,:)  = mbgt_arr(mb,:,:)
;        end do

;      end if

;============================================
; Define output file path and plot details 
;============================================

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
     opts_mbgt@cnFillPalette		   = "ncl_default"
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

    opts_mbgt@tmXTOn                       = "False" ; Turn off top x-axis TM
    opts_mbgt@tmYROn                       = "False" ; Turn off right y-axis TM
    opts_mbgt@gsnMaximize                  = False   ; Maximise plot size
    opts_mbgt@gsnAddCyclic                 = False
    opts_mbgt@tmXBMode                     = "Explicit"  ; Set tick marks explicitly

    if (rad0 .eq. 0) then
     opts_mbgt@tmXBValues                   = (/0.0,0.25,0.5,0.75,1.0,1.25,1.5/)
     opts_mbgt@tmXBLabels                   = (/"0","25","50","75","100","125","150"/)
    elseif (rad0 .eq. 1) then
     opts_mbgt@tmXBValues                   = (/0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0/)
     opts_mbgt@tmXBLabels                   = (/"0","25","50","75","100","125","150","175","200"/)
    end if

    opts_mbgt@tmYLMode			   = "Explicit"
    opts_mbgt@gsnYAxisIrregular2Linear     = True       ; Linearise the y-axis
    opts_mbgt@trYMinF                      = 0.0        ; y-axis minimum value (0.0)
    opts_mbgt@trYMaxF                      = 6.0        ; y-axis maximum value (17.6)

  ; Set y-axis labels
    opts_mbgt@tmYLValues                   = fspan(opts_mbgt@trYMinF+1, opts_mbgt@trYMaxF, 6)
    opts_mbgt@tmYLLabels                   = tostring(opts_mbgt@tmYLValues)

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
    if (clr .eq. 4 .or. clr .eq. 5) then
     opts_vvel@cnLineColor		   = "yellow"
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
    opts_zero@cnLineDashPattern           = 3
    opts_zero@cnLevelSelectionMode        = "ExplicitLevels"
    opts_zero@cnLevels			  = 0.0
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
    opts_vrad2@cnFillPalette	           = "BlueDarkRed18"
    opts_vrad2@cnLevels			   = (/-10., -5., -3., -2., -1., -0.5, -0.2, \
                                                0, 0.2, 0.5, 1., 2., 3., 5., 10./)
    opts_vrad2@cnFillColors		   = (/0,2,3,4,5,6,7,-1,-1,10,11,12,13,14,15,17/)
    opts_vrad2@lbLabelBarOn                = True
    opts_vrad2@lbBoxEndCapStyle            = "TriangleBothEnds"
    opts_vrad2@lbOrientation		   = "Vertical"
    opts_vrad2@vpXF			   = 0.07
    opts_vrad2@vpYF			   = 0.90
    opts_vrad2@vpWidthF                    = 0.15
    opts_vrad2@vpHeightF		   = 0.15

  ; Tangential wind 
    opts_vtan2				   = opts_mbgt
    delete(opts_vtan2@cnLevels)
    delete(opts_vtan2@cnFillColors)
    opts_vtan2@cnFillPalette		   = "radar_new"

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
    opts_aam				   = opts_vtan2
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
     txid			= new(14,graphic)
     amid	    	        = new(14,graphic)
    elseif (full .eq. 0) then 
     txid                       = new(2,graphic)
     amid                       = new(2,graphic)
    end if

  ; Label text resources 
    txres                       = True
    txres@txPerimOn             = True
    txres@txFontHeightF         = 0.010
    txres@txBackgroundFillColor = "White"

  ; Label position resources 
    amres 	   	        = True
    amres@amParallelPosF	= 0.5			; Right edge (-0.5 for left edge)
    amres@amOrthogonalPosF	= -0.5			; Top edge
    amres@amJust		= "TopRight"		; Switch to "TopLeft" if required

    if (full .eq. 1) then
     if (ueh .eq. 0) then 
      panel_strings		 = (/"a) vrad", "b) vtan", "c) vvel", "d) aam", \
                                     "e) U_magf", "f) U_mr", "g) U_mv", \
                                     "h) U_eagf", "i) U_dh", \
                                     "j) U_dz", "k) U_ev", "l) U_t", \
                                     "m) U_t + U_eh + U_mr", "n) All U tend"/)
     elseif (ueh .eq. 1) then 
      panel_strings              = (/"a) vrad", "b) vtan", "c) vvel", "d) aam", \
                                     "e) U_magf", "f) U_mr", "g) U_mv", \
                                     "h) U_eagf", "i) U_eh", \
                                     "j) U_dz + U_dh", "k) U_ev", "l) U_t", \
                                     "m) U_t + U_eh + U_mr", "n) All U tend"/)
     end if 
    elseif (full .eq. 0) then 
     panel_strings              = (/"a) Mean terms","b) Eddy terms"/)
    end if 

;=================
; Plot the data 
;=================

  ; Calculating the terms from Montgomery et al. (2018) - their Eq.(8)
  ; LHS (1) : d{u}/dt          - Local tendency of the mean radial wind
  ; LHS (2) : {u} * d{u}/dr    - Mean radial advection of rad mom
  ; LHS (3) : {(u' * du'/dr) + (v'/r * du'/dy)}  - Mean horiz. adv. of eddy rad mom
 
  ; RHS (1) : -{w} * d{u}/dz   - Mean vertical advection of mean rad mom
  ; RHS (2) : -{w' * d{u'}/dz} - Eddy vertical advection of eddy rad mom
  ; RHS (3) : {v}^2/r + f{v} - 1/{rho} * d{p}/dr - Mean agradient force
  ; RHS (4) : {v'^2/r - 1/rho * dp'/dr}          - Eddy agradient force
  ; RHS (5) : {F_r}            - Mean radial/vertical diffusive tendency of rad mom

  if (ave .eq. 0) then 

   do it = ts0, tf0

    print("Working on time: "+it)

  ; Output file location and type
    output = "$sam/nepartak/images/cyl_coords/mbgt5_rad_xz_"+dat+"_"+ens0+\
             "_"+cn0+"_"+time_arr(it)
    wks = gsn_open_wks(opt,output)

  if (full .eq. 1) then 

  ; Radial wind (filled contours)
    panel_vrad  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad2) ; [u]
    panel_zero	= gsn_csm_contour(wks,mean_vrad(:,it,:),opts_zero)  ; Zero line contour
    zero0	= ColorNegDashZeroPosContour(panel_zero,"transparent","black","transparent")
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

  ; Absolute angular momentum (filled contours)
    panel_aam	= gsn_csm_contour(wks,mean_aam(:,it,:),opts_aam)   ; [AAM]
    txid(3)     = gsn_create_text(wks, panel_strings(3), txres)
    amid(3)     = gsn_add_annotation(panel_aam, txid(3), amres)

  ; Mean agradient force [U_magf]
    opts_mbgt@vpXF       = 0.05
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel0      = gsn_csm_contour(wks,rhs3(:,it,:),opts_mbgt)	   ; [U_magf]
    cont0	= gsn_csm_contour(wks,rhs3(:,it,:),opts_zero)	   ; Draw contours

    panel0_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad) ; [u]
    panel0_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel) ; [w] 

    overlay(panel0, cont0)					   ; Overlay contours
    overlay(panel0, panel0_ovr)
    overlay(panel0, panel0_ovr2)
    txid(4)     = gsn_create_text(wks, panel_strings(4), txres)
    amid(4)     = gsn_add_annotation(panel0, txid(4), amres)
    delete([/cont0/])

  ; Mean radial advection of mean radial momentum [U_mr]
    opts_mbgt@vpXF       = 0.24
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel1      = gsn_csm_contour(wks,lhs2(:,it,:),opts_mbgt)      ; [U_mr]
    cont0	= gsn_csm_contour(wks,lhs2(:,it,:),opts_zero)      ; Draw contours

    panel1_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad) ; [U]
    panel1_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel) ; [w]

    overlay(panel1, cont0)                                         ; Overlay contours
    overlay(panel1, panel1_ovr)
    overlay(panel1, panel1_ovr2)
    txid(5)     = gsn_create_text(wks, panel_strings(5), txres)
    amid(5)     = gsn_add_annotation(panel1, txid(5), amres)
    delete([/cont0/])

  ; Minus the mean vertical advection of mean radial momentum [U_mv]
    opts_mbgt@vpXF       = 0.44
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel2      = gsn_csm_contour(wks,rhs1(:,it,:),opts_mbgt)      ; [U_mv]
    cont0       = gsn_csm_contour(wks,rhs1(:,it,:),opts_zero)      ; Draw contours

    panel2_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad) ; [U]
    panel2_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel) ; [w]

    overlay(panel2, cont0)                                         ; Overlay contours
    overlay(panel2, panel2_ovr)
    overlay(panel2, panel2_ovr2)
    txid(6)     = gsn_create_text(wks, panel_strings(6), txres)
    amid(6)     = gsn_add_annotation(panel2, txid(6), amres)
    delete([/cont0/])

  ; Eddy agradient force [U_eagf]
    opts_mbgt@vpXF       = 0.64
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel3      = gsn_csm_contour(wks,rhs4(:,it,:),opts_mbgt)      ; [U_eagf]
    cont0       = gsn_csm_contour(wks,rhs4(:,it,:),opts_zero)      ; Draw contours

    panel3_ovr = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad)  ; [U]
    panel3_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel) ; [w]

    overlay(panel3, cont0)                                         ; Overlay contours
    overlay(panel3, panel3_ovr)
    overlay(panel3, panel3_ovr2)
    txid(7)     = gsn_create_text(wks, panel_strings(7), txres)
    amid(7)     = gsn_add_annotation(panel3, txid(7), amres)
    delete([/cont0/])

  ; Radial diffusive tendency of radial momentum [U_dh] OR 
  ; Eddy horizontal advection of eddy radial momentum [U_eh]
    opts_mbgt@vpXF       = 0.83
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15

    if (ueh .eq. 0) then 

     panel4      = gsn_csm_contour(wks,rhs5a(:,it,:),opts_mbgt)     ; [U_dh]
     cont0       = gsn_csm_contour(wks,rhs5a(:,it,:),opts_zero)     ; Draw contours
    
    elseif (ueh .eq. 1) then 

     panel4      = gsn_csm_contour(wks,lhs3(:,it,:),opts_mbgt)      ; [U_eh]
     cont0       = gsn_csm_contour(wks,lhs3(:,it,:),opts_zero)      ; Draw contours

    end if 

    panel4_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad) ; [U]
    panel4_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel) ; [w]

    overlay(panel4, cont0)                                         ; Overlay contours
    overlay(panel4, panel4_ovr)
    overlay(panel4, panel4_ovr2)
    txid(8)     = gsn_create_text(wks, panel_strings(8), txres)
    amid(8)     = gsn_add_annotation(panel4, txid(8), amres)
    delete([/cont0/])

  ; Vertical diffusive tendency of radial momentum [U_dz] (total tendency if 'ueh' = 1)
    opts_mbgt@vpXF       = 0.05
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15

    if (ueh .eq. 0) then

     panel5      = gsn_csm_contour(wks,rhs5b(:,it,:),opts_mbgt)     ; [U_dz]
     cont0       = gsn_csm_contour(wks,rhs5b(:,it,:),opts_zero)     ; Draw contours

    elseif (ueh .eq. 1) then

     panel5      = gsn_csm_contour(wks,diff(:,it,:),opts_mbgt)      ; [U_dz]
     cont0       = gsn_csm_contour(wks,diff(:,it,:),opts_zero)      ; Draw contours

    end if 

    panel5_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad)  ; [U]
    panel5_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel)  ; [w]

    overlay(panel5, cont0)                                          ; Overlay contours
    overlay(panel5, panel5_ovr)
    overlay(panel5, panel5_ovr2)
    txid(9)     = gsn_create_text(wks, panel_strings(9), txres)
    amid(9)     = gsn_add_annotation(panel5, txid(9), amres)
    delete([/cont0/])

  ; Minus the eddy vertical advection of eddy radial momentum [U_ev]
    opts_mbgt@vpXF       = 0.24
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel6      = gsn_csm_contour(wks,rhs2(:,it,:),opts_mbgt)      ; [U_ev]
    cont0       = gsn_csm_contour(wks,rhs2(:,it,:),opts_zero)      ; Draw contours

    panel6_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad) ; [U]
    panel6_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel) ; [w]

    overlay(panel6, cont0)                                         ; Overlay contours
    overlay(panel6, panel6_ovr)
    overlay(panel6, panel6_ovr2)
    txid(10)    = gsn_create_text(wks, panel_strings(10), txres)
    amid(10)    = gsn_add_annotation(panel6, txid(10), amres)
    delete([/cont0/])

  ; Local tendency of the mean radial wind [U_t]
    opts_mbgt@vpXF       = 0.44
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel7      = gsn_csm_contour(wks,lhs1(:,it,:),opts_mbgt)      ; [U_t]
    cont0       = gsn_csm_contour(wks,lhs1(:,it,:),opts_zero)      ; Draw contours

    panel7_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad) ; [U]
    panel7_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel) ; [w]

    overlay(panel7, cont0)                                         ; Overlay contours
    overlay(panel7, panel7_ovr)
    overlay(panel7, panel7_ovr2)
    txid(11)    = gsn_create_text(wks, panel_strings(11), txres)
    amid(11)    = gsn_add_annotation(panel7, txid(11), amres)
    delete([/cont0/])

  ; Total LHS [U_t + U_eh + U_mr]
    opts_mbgt@vpXF       = 0.64
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel8      = gsn_csm_contour(wks,lhsT(:,it,:),opts_mbgt)      ; [U_t + U_eh + U_mr]
    cont0       = gsn_csm_contour(wks,lhsT(:,it,:),opts_zero)      ; Draw contours

    panel8_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad) ; [U]
    panel8_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel) ; [w]

    overlay(panel8, cont0)                                         ; Overlay contours
    overlay(panel8, panel8_ovr)
    overlay(panel8, panel8_ovr2)
    txid(12)    = gsn_create_text(wks, panel_strings(12), txres)
    amid(12)    = gsn_add_annotation(panel8, txid(12), amres)
    delete([/cont0/])

  ; Total RHS [All U tend]
    opts_mbgt@vpXF       = 0.83
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel9      = gsn_csm_contour(wks,rhsT(:,it,:),opts_mbgt)      ; [All U tend]
    cont0       = gsn_csm_contour(wks,rhsT(:,it,:),opts_zero)      ; Draw contours

    panel9_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad) ; [U]
    panel9_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel) ; [w]

    overlay(panel9, cont0)                                         ; Overlay contours
    overlay(panel9, panel9_ovr)
    overlay(panel9, panel9_ovr2)
    txid(13)    = gsn_create_text(wks, panel_strings(13), txres)
    amid(13)    = gsn_add_annotation(panel9, txid(13), amres)
    delete([/cont0/])

  elseif (full .eq. 0) then 

  ; All eddy terms
    opts_mbgt@vpXF       = 0.05
    opts_mbgt@vpYF       = 0.75
    opts_mbgt@vpWidthF   = 0.42
    opts_mbgt@vpHeightF  = 0.42
    panel_eddy  = gsn_csm_contour(wks,eddy(:,it,:),opts_mbgt)          ; [All U tend]
    cont0       = gsn_csm_contour(wks,eddy(:,it,:),opts_zero)          ; Draw contours

    panel_eddy_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad) ; [U]
    panel_eddy_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel) ; [w]

    overlay(panel_eddy, cont0)                                         ; Overlay contours
    overlay(panel_eddy, panel_eddy_ovr)
    overlay(panel_eddy, panel_eddy_ovr2)
    delete([/cont0/])

  ; All mean terms 
    opts_mbgt@vpXF       = 0.55
    opts_mbgt@vpYF       = 0.75
    opts_mbgt@vpWidthF   = 0.42
    opts_mbgt@vpHeightF  = 0.42
    panel_mean  = gsn_csm_contour(wks,mean(:,it,:),opts_mbgt)          ; [All U tend]
    cont0       = gsn_csm_contour(wks,mean(:,it,:),opts_zero)          ; Draw contours

    panel_mean_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad) ; [U]
    panel_mean_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel) ; [w]

    overlay(panel_mean, cont0)                                         ; Overlay contours
    overlay(panel_mean, panel_mean_ovr)
    overlay(panel_mean, panel_mean_ovr2)
    delete([/cont0/])

  end if 

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

  ; Create a main title at top of frame                                                 
    title_res                   = True
    title_str                   = "Radial momentum budget: "+time_str(it)
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
    label_str                   = (/"-20.0","-15.0","-10.0","-5.0","-2.0","-1.0","-0.5","0",\
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

    pres			= True 
    maximize_output(wks,pres)				   ; Calls 'draw' and 'frame'

   end do     ; End time loop (do it = 0, ts, tf)

  elseif (ave .eq. 1) then 

  ; Define output file path and type
    if (full .eq. 0) then
     output = "$sam/nepartak/images/cyl_coords/mbgt5_rad_ave_xz_em_"+dat+\
     	      "_"+ens0+"_"+ts0+"_"+tf0+"_"+cn0
    else
     output = "$sam/nepartak/images/cyl_coords/mbgt5_rad_ave_xz_"+dat+\
     	      "_"+ens0+"_"+ts0+"_"+tf0+"_"+cn0
    end if

    wks = gsn_open_wks(opt,output)

   if (full .eq. 1) then 

  ; Radial wind (filled contours)
    panel_vrad  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad2) ; [u]
    panel_zero	= gsn_csm_contour(wks,vrad_av0(:,:),opts_zero)  ; Zero line contour
    zero0	= ColorNegDashZeroPosContour(panel_zero,"transparent","black","transparent")
    overlay(panel_vrad,	zero0)
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

  ; Mean agradient force [U_magf]
    opts_mbgt@vpXF       = 0.05
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel0      = gsn_csm_contour(wks,rhs3_av(:,:),opts_mbgt)  ; [U_magf]
    cont0       = gsn_csm_contour(wks,rhs3_av(:,:),opts_zero)  ; Draw contours

    panel0_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad) ; [u]
    panel0_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel) ; [w]

    overlay(panel0, cont0)                                     ; Overlay contours
    overlay(panel0, panel0_ovr)
    overlay(panel0, panel0_ovr2)
    txid(4)     = gsn_create_text(wks, panel_strings(4), txres)
    amid(4)     = gsn_add_annotation(panel0, txid(4), amres)
    delete([/cont0/])

  ; Mean radial advection of mean radial momentum [U_mr]
    opts_mbgt@vpXF       = 0.24
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel1      = gsn_csm_contour(wks,lhs2_av(:,:),opts_mbgt)  ; [U_mr]
    cont0       = gsn_csm_contour(wks,lhs2_av(:,:),opts_zero)  ; Draw contours

    panel1_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad) ; [U]
    panel1_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel) ; [w]

    overlay(panel1, cont0)                                     ; Overlay contours
    overlay(panel1, panel1_ovr)
    overlay(panel1, panel1_ovr2)
    txid(5)     = gsn_create_text(wks, panel_strings(5), txres)
    amid(5)     = gsn_add_annotation(panel1, txid(5), amres)
    delete([/cont0/])

  ; Minus the mean vertical advection of mean radial momentum [U_mv]
    opts_mbgt@vpXF       = 0.44
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel2      = gsn_csm_contour(wks,rhs1_av(:,:),opts_mbgt)      ; [U_mv]
    cont0       = gsn_csm_contour(wks,rhs1_av(:,:),opts_zero)      ; Draw contours

    panel2_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad) 	   ; [U]
    panel2_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel) 	   ; [w]

    overlay(panel2, cont0)                                         ; Overlay contours
    overlay(panel2, panel2_ovr)
    overlay(panel2, panel2_ovr2)
    txid(6)     = gsn_create_text(wks, panel_strings(6), txres)
    amid(6)     = gsn_add_annotation(panel2, txid(6), amres)
    delete([/cont0/])

  ; Eddy agradient force [U_eagf]
    opts_mbgt@vpXF       = 0.64
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel3      = gsn_csm_contour(wks,rhs4_av(:,:),opts_mbgt)  ; [U_eagf]
    cont0       = gsn_csm_contour(wks,rhs4_av(:,:),opts_zero)  ; Draw contours

    panel3_ovr = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad)  ; [U]
    panel3_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel) ; [w]

    overlay(panel3, cont0)                                     ; Overlay contours
    overlay(panel3, panel3_ovr)
    overlay(panel3, panel3_ovr2)
    txid(7)     = gsn_create_text(wks, panel_strings(7), txres)
    amid(7)     = gsn_add_annotation(panel3, txid(7), amres)
    delete([/cont0/])

  ; Radial diffusive tendency of radial momentum [U_dh] OR 
  ; Eddy horizontal advection of eddy radial momentum [U_eh]
    opts_mbgt@vpXF       = 0.83
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15

    if (ueh .eq. 0) then 

     panel4      = gsn_csm_contour(wks,rhs5a_av(:,:),opts_mbgt) ; [U_dh]
     cont0       = gsn_csm_contour(wks,rhs5a_av(:,:),opts_zero) ; Draw contours

    elseif (ueh .eq. 1) then 

     panel4      = gsn_csm_contour(wks,lhs3_av(:,:),opts_mbgt)  ; [U_eh]
     cont0       = gsn_csm_contour(wks,lhs3_av(:,:),opts_zero)  ; Draw contours

    end if 

    panel4_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad)  ; [U]
    panel4_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel)  ; [w]

    overlay(panel4, cont0)                                      ; Overlay contours
    overlay(panel4, panel4_ovr)
    overlay(panel4, panel4_ovr2)
    txid(8)     = gsn_create_text(wks, panel_strings(8), txres)
    amid(8)     = gsn_add_annotation(panel4, txid(8), amres)
    delete([/cont0/])

  ; Vertical diffusive tendency of radial momentum [U_dz]
    opts_mbgt@vpXF       = 0.05
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15

    if (ueh .eq. 0) then

     panel5      = gsn_csm_contour(wks,rhs5b_av(:,:),opts_mbgt) ; [U_dz]
     cont0       = gsn_csm_contour(wks,rhs5b_av(:,:),opts_zero) ; Draw contours

    elseif (ueh .eq. 1) then

     panel5      = gsn_csm_contour(wks,diff_av(:,:),opts_mbgt)  ; [U_dz]
     cont0       = gsn_csm_contour(wks,diff_av(:,:),opts_zero)  ; Draw contours

    end if

    panel5_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad)  ; [U]
    panel5_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel)  ; [w]

    overlay(panel5, cont0)                                      ; Overlay contours
    overlay(panel5, panel5_ovr)
    overlay(panel5, panel5_ovr2)
    txid(9)     = gsn_create_text(wks, panel_strings(9), txres)
    amid(9)     = gsn_add_annotation(panel5, txid(9), amres)
    delete([/cont0/])

  ; Minus the eddy vertical advection of eddy radial momentum [U_ev]
    opts_mbgt@vpXF       = 0.24
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel6      = gsn_csm_contour(wks,rhs2_av(:,:),opts_mbgt)      ; [U_ev]
    cont0       = gsn_csm_contour(wks,rhs2_av(:,:),opts_zero)      ; Draw contours

    panel6_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad) 	   ; [U]
    panel6_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel)     ; [w]
    
    overlay(panel6, cont0)                                         ; Overlay contours
    overlay(panel6, panel6_ovr)
    overlay(panel6, panel6_ovr2)
    txid(10)    = gsn_create_text(wks, panel_strings(10), txres)
    amid(10)    = gsn_add_annotation(panel6, txid(10), amres)
    delete([/cont0/])

  ; Local tendency of the mean radial wind [U_t]
    opts_mbgt@vpXF       = 0.44
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel7      = gsn_csm_contour(wks,lhs1_av(:,:),opts_mbgt)      ; [U_t]
    cont0       = gsn_csm_contour(wks,lhs1_av(:,:),opts_zero)      ; Draw contours

    panel7_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad) 	   ; [U]
    panel7_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel) 	   ; [w]

    overlay(panel7, cont0)                                         ; Overlay contours
    overlay(panel7, panel7_ovr)
    overlay(panel7, panel7_ovr2)
    txid(11)    = gsn_create_text(wks, panel_strings(11), txres)
    amid(11)    = gsn_add_annotation(panel7, txid(11), amres)
    delete([/cont0/])

  ; Total LHS [U_t + U_eh + U_mr]
    opts_mbgt@vpXF       = 0.64
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel8      = gsn_csm_contour(wks,lhsT_av(:,:),opts_mbgt)      ; [U_t + U_eh + U_mr]
    cont0       = gsn_csm_contour(wks,lhsT_av(:,:),opts_zero)      ; Draw contours

    panel8_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad) 	   ; [U]
    panel8_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel)     ; [w]

    overlay(panel8, cont0)                                         ; Overlay contours
    overlay(panel8, panel8_ovr)
    overlay(panel8, panel8_ovr2)
    txid(12)    = gsn_create_text(wks, panel_strings(12), txres)
    amid(12)    = gsn_add_annotation(panel8, txid(12), amres)
    delete([/cont0/])

  ; Total RHS [All U tend]
    opts_mbgt@vpXF       = 0.83
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel9      = gsn_csm_contour(wks,rhsT_av(:,:),opts_mbgt)      ; [All U tend]
    cont0       = gsn_csm_contour(wks,rhsT_av(:,:),opts_zero)      ; Draw contours

    panel9_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad) 	   ; [U]
    panel9_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel) 	   ; [w]

    overlay(panel9, cont0)                                         ; Overlay contours
    overlay(panel9, panel9_ovr)
    overlay(panel9, panel9_ovr2)
    txid(13)    = gsn_create_text(wks, panel_strings(13), txres)
    amid(13)    = gsn_add_annotation(panel9, txid(13), amres)
    delete([/cont0/])

  elseif (full .eq. 0) then 

  ; All eddy terms 
    opts_mbgt@vpXF       = 0.05
    opts_mbgt@vpYF       = 0.75
    opts_mbgt@vpWidthF   = 0.42
    opts_mbgt@vpHeightF  = 0.42
    panel_eddy  = gsn_csm_contour(wks,eddy_av(:,:),opts_mbgt)      ; [All eddy terms]
    cont0       = gsn_csm_contour(wks,eddy_av(:,:),opts_zero)      ; Draw contours

    panel_eddy_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad) ; [U]
    panel_eddy_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel) ; [w]

    overlay(panel_eddy, cont0)					   ; Overlay contours
    overlay(panel_eddy, panel_eddy_ovr)
    overlay(panel_eddy, panel_eddy_ovr2)
    txid(0)         = gsn_create_text(wks, panel_strings(0), txres)
    amid(0)         = gsn_add_annotation(panel_eddy, txid(0), amres)
    delete([/cont0/])

  ; All mean terms 
    opts_mbgt@vpXF       = 0.55
    opts_mbgt@vpYF       = 0.75
    opts_mbgt@vpWidthF   = 0.42
    opts_mbgt@vpHeightF  = 0.42
    panel_mean      = gsn_csm_contour(wks,mean_av(:,:),opts_mbgt)  ; [All mean terms]
    cont0           = gsn_csm_contour(wks,mean_av(:,:),opts_zero)  ; Draw contours

    panel_mean_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad) ; [U]
    panel_mean_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel) ; [w]

    overlay(panel_mean, cont0)                                     ; Overlay contours
    overlay(panel_mean, panel_mean_ovr)
    overlay(panel_mean, panel_mean_ovr2)
    txid(1)         = gsn_create_text(wks, panel_strings(1), txres)
    amid(1)         = gsn_add_annotation(panel_mean, txid(1), amres)
    delete([/cont0/])

   end if 	      		; End IF statement (full .eq. 1)

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

  ; Create a main title at top of frame                                                   
    title_res                   = True
    title_str                   = "Radial momentum budget: "+time_str(ts0)+" to "\
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
    label_str                   = (/"-10.0","-5.0","-3.0","-2.0","-1.0","-0.5","-0.2", "0.0",\
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

  end if 			; End IF statement (if (ave .eq. 1) then ...)

end