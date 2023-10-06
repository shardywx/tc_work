; Script to calculate the radial momentum budget as in Montgomery et al. (2018)

; Run using:

; ncl opt=\"x11\" dist=1.04 rad0=0 ave=1 nt=73 nr=11 ar=0.5 lgd0=\"tr\" 
; wrt=1 new0=0 cn0=\"geo_sm\" sm=0 sc=5 full=9 var0=\"div1\" typ=\"sym\" tc_plt_xy.ncl

; 'opt'   = output file format ("pdf" or "x11")
; 'dist'  = size of box following storm (degrees)
; 'rad0'  = radius of final plots: 150 km (0) or 200 km (1)
; 'ave'	  = calculate 3/6/12-h average (1) rather than analyse individual time interval (0)
; 'wrt'	  = write out averages to text file (1) or leave alone (0 – default)
; 'typ'	  = VC phase ("sym","asym","stoa","atos")
; 'full'  = top 18 (1); top 8 (-5); top 6 (-6)
; 'ar'    = distance in degrees of outer radius (1.0, 2.0, etc)
; 'nr'    = number of radial circles between r = 0 and r = 'ar' (21, 41, 61, etc)
; 'nt'    = number of azimuth angles in cylindrical grid
; 'new0'  = reduced scale for plotting tangential wind and AAM (1)
; 'cn0'   = offline storm centre: "slp", "vort", "geo", "geo_sm"

; 500 m [10], 1 km [14], 1.5 km [17], 2 km [20], 3 km [24], 4 km [28], 5 km [31], 6 km [34]  
; 7 km [37], 8 km [39], 9 km [42], 10 km [45], 11 km [46], 12 km [48], 13 km [50], 14 km [52]
; 15 km [54], 16 km [56]

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

; Dimensions of text files to read data in from (S-phase, S-A, A-phase, A-S)    
; 'vc_arr_all(:,0)' = Group 1 [17, 16, 16, 9]                                                       
; 'vc_arr_all(:,1)' = Group 2 [25, 14, 29, 15]
; 'vc_arr_all(:,2)' = Group 3 [31, 30, 31, 31]
; 'vc_arr_all(:,3)' = Group 4 [29, 24, 38, 24]
  vc_arr_all = (/ (/17,25,31,29,24/), (/16,14,30,24,24/), (/16,29,31,38,24/), (/9,15,31,24,24/) /)
  vc_arr_nw  = (/ (/14,19,31/), (/12,12,30/), (/13,24,31/), (/9,12,31/)  /)

; EDIT -- are the time intervals correct (it-2, etc)?
;  ct = 0 

; Correct description of VC phase 
  if (typ.eq. "sym") then
   typ0 = "ring"
  elseif (typ .eq. "stoa") then
   typ0 = "r2m"
  elseif (typ .eq. "asym") then
   typ0 = "mono"
  elseif (typ .eq. "atos") then
   typ0 = "m2r"
  end if

; For now, explicitly tell NCL which dates and simulations to analyse
  if (typ .eq. "sym") then

   if (full .eq. 9) then

    dat  = asciiread("$ar/text/vc_sym_dat_group4a.txt",(/vc_arr_all(0,4)/),"string")
    ens0 = asciiread("$ar/text/vc_sym_sim_group4a.txt",(/vc_arr_all(0,4)/),"string")
    ts0  = asciiread("$ar/text/vc_sym_ts_group4a.txt",(/vc_arr_all(0,4)/),"integer")
    tf0  = ts0 + 1

   elseif (full .eq. 10) then

    dat  = asciiread("$ar/text/vc_sym_dat_group4b.txt",(/vc_arr_all(0,4)/),"string")
    ens0 = asciiread("$ar/text/vc_sym_sim_group4b.txt",(/vc_arr_all(0,4)/),"string")
    ts0  = asciiread("$ar/text/vc_sym_ts_group4b.txt",(/vc_arr_all(0,4)/),"integer")
    tf0  = ts0 + 1

   elseif (full .eq. 11) then

    dat  = asciiread("$ar/text/vc_sym_dat_group4c.txt",(/vc_arr_all(0,4)/),"string")
    ens0 = asciiread("$ar/text/vc_sym_sim_group4c.txt",(/vc_arr_all(0,4)/),"string")
    ts0  = asciiread("$ar/text/vc_sym_ts_group4c.txt",(/vc_arr_all(0,4)/),"integer")
    tf0  = ts0 + 1

   end if

  elseif (typ .eq. "stoa") then

   if (full .eq. 9) then

    dat  = asciiread("$ar/text/vc_s2a_dat_group4a.txt",(/vc_arr_all(1,4)/),"string")
    ens0 = asciiread("$ar/text/vc_s2a_sim_group4a.txt",(/vc_arr_all(1,4)/),"string")
    ts0  = asciiread("$ar/text/vc_s2a_ts_group4a.txt",(/vc_arr_all(1,4)/),"integer")
    tf0  = ts0 + 1

   elseif (full .eq. 10) then

    dat  = asciiread("$ar/text/vc_s2a_dat_group4b.txt",(/vc_arr_all(1,4)/),"string")
    ens0 = asciiread("$ar/text/vc_s2a_sim_group4b.txt",(/vc_arr_all(1,4)/),"string")
    ts0  = asciiread("$ar/text/vc_s2a_ts_group4b.txt",(/vc_arr_all(1,4)/),"integer")
    tf0  = ts0 + 1

   elseif (full .eq. 11) then

    dat  = asciiread("$ar/text/vc_s2a_dat_group4c.txt",(/vc_arr_all(1,4)/),"string")
    ens0 = asciiread("$ar/text/vc_s2a_sim_group4c.txt",(/vc_arr_all(1,4)/),"string")
    ts0  = asciiread("$ar/text/vc_s2a_ts_group4c.txt",(/vc_arr_all(1,4)/),"integer")
    tf0  = ts0 + 1

   end if 

  elseif (typ .eq. "asym") then

   if (full .eq. 9) then

    dat  = asciiread("$ar/text/vc_asym_dat_group4a.txt",(/vc_arr_all(2,4)/),"string")
    ens0 = asciiread("$ar/text/vc_asym_sim_group4a.txt",(/vc_arr_all(2,4)/),"string")
    ts0  = asciiread("$ar/text/vc_asym_ts_group4a.txt",(/vc_arr_all(2,4)/),"integer")
    tf0  = ts0 + 1

   elseif (full .eq. 10) then

    dat  = asciiread("$ar/text/vc_asym_dat_group4b.txt",(/vc_arr_all(2,4)/),"string")
    ens0 = asciiread("$ar/text/vc_asym_sim_group4b.txt",(/vc_arr_all(2,4)/),"string")
    ts0  = asciiread("$ar/text/vc_asym_ts_group4b.txt",(/vc_arr_all(2,4)/),"integer")
    tf0  = ts0 + 1

   elseif (full .eq. 11) then

    dat  = asciiread("$ar/text/vc_asym_dat_group4c.txt",(/vc_arr_all(2,4)/),"string")
    ens0 = asciiread("$ar/text/vc_asym_sim_group4c.txt",(/vc_arr_all(2,4)/),"string")
    ts0  = asciiread("$ar/text/vc_asym_ts_group4c.txt",(/vc_arr_all(2,4)/),"integer")
    tf0  = ts0 + 1

   end if 

  elseif (typ .eq. "atos") then

   if (full .eq. 9) then

    dat  = asciiread("$ar/text/vc_a2s_dat_group4a.txt",(/vc_arr_all(3,4)/),"string")
    ens0 = asciiread("$ar/text/vc_a2s_sim_group4a.txt",(/vc_arr_all(3,4)/),"string")
    ts0  = asciiread("$ar/text/vc_a2s_ts_group4a.txt",(/vc_arr_all(3,4)/),"integer")
    tf0  = ts0 + 1

   elseif (full .eq. 10) then

    dat  = asciiread("$ar/text/vc_a2s_dat_group4b.txt",(/vc_arr_all(3,4)/),"string")
    ens0 = asciiread("$ar/text/vc_a2s_sim_group4b.txt",(/vc_arr_all(3,4)/),"string")
    ts0  = asciiread("$ar/text/vc_a2s_ts_group4b.txt",(/vc_arr_all(3,4)/),"integer")
    tf0  = ts0 + 1

   elseif (full .eq. 11) then

    dat  = asciiread("$ar/text/vc_a2s_dat_group4c.txt",(/vc_arr_all(3,4)/),"string")
    ens0 = asciiread("$ar/text/vc_a2s_sim_group4c.txt",(/vc_arr_all(3,4)/),"string")
    ts0  = asciiread("$ar/text/vc_a2s_ts_group4c.txt",(/vc_arr_all(3,4)/),"integer")
    tf0  = ts0 + 1

   end if

  end if

;==========================================================
; Now start loop over these chosen input files and times
;==========================================================

 do st = 0, dimsizes(dat)-1 

  dt0 = tf0(st) - ts0(st)

; Create array to hold time-averaged values for composites
  if (wrt .eq. 1) then 
   print("Creating composite chart array...")

   if (st .eq. 0) then

 ; Dimensions of 4D array (levels; time averages; diagnostics; radial circles)
   nlev	       = 59		
   ntimes      = dimsizes(dat)	
   ndiag       = 11		
   nrad	       = nr		

 ; Create array to hold time averages 
   plot_arr    = new( (/nlev, ntimes, ndiag, nrad/), "float")

   end if 

  end if 

; Find 'pc' and 'pd' files (model height level data)
  diri        = "/nfs/a319/earshar/"+dat(st)+"/"+ens0(st)
  fili_prefix = "201607"+dat(st)+"00Z_NPTK_4p4_L80_ra1t_"+ens0(st)
  fili_c1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pc*.nc")
  fili_c      = diri+"/"+fili_c1
  fili_d1     =	systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pd*.nc")
  fili_d      = diri+"/"+fili_d1

; Calculate number of times in each file using built-in function ('nc_times')
  numINFO_c   = nc_times(fili_c)
  numINFO_d   = nc_times(fili_d)

; Number/array of times from each file
  numINPUT_c  = numINFO_c[0]
  time_c      = numINFO_c[1]
  numINPUT_d  = numINFO_d[0]
  time_d      = numINFO_d[1]

;============================================================
; Calculate storm motion using built-in function ('st_rm')
;============================================================

; 'dat'       = initialisation time (02T12, 03T00, ...)
; 'ens0'      = ensemble simulation (em00, em01, ...)
; 'diri'      = path to input files (see above)
; 'dist'      = size of box following storm (degrees)
; 'mins'      = analyse 1-h (0) or 5-min (1) data

; Zonal; meridional; vector wind; times; lat; lon; centre
  storm_rel = st_rm(dat(st), ens0(st), "$ar/text/", dist, 0)
  u_cyc     = storm_rel[0]
  v_cyc     = storm_rel[1]
  vel_cyc   = storm_rel[2]
  numTIMES  = storm_rel[3]
  lat_arr   = storm_rel[4]
  lon_arr   = storm_rel[5]
  centre    = storm_rel[6]

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
  do ct0 = 0, numINPUT_c-1
   if (minute(ct0).gt.30) then
     hour(ct0) = hour(ct0)+1
   end if
  end do
  delete(ct0)

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

  pii        = 3.14159265
  radii      = fspan(0.0,ar,nr)      ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,nt)   ; Array of azimuth angles (0, 5, 10, etc)
  thetas_rad = (pii/180)*thetas      ; Degrees to radians
  rad_size   = dimsizes(radii)-1     ; Index for later calculation
  num_l      = 59      	   	     ; Number of model levels (edit to evolve smoothly)

; Define constants
  rd	     = 287.0	             ; Gas constant for dry air (J/kg/K^2)
  cp0	     = 1004.0                ; Specific heat of dry air at constant pressure (J/kg/K)
  kp	     = 0.286                 ; For Exner function calculations (0.286)
  pref	     = 1000.0                ; Reference pressure p0 (1000 hPa)

;========================================================================================
;========================================================================================
; Calculating the terms from Montgomery et al. (2018) - their Eq.(4)
; Ignore the perturbation PGF term because it is several orders of magnitude smaller

; LHS (1) : d{v}/dt          - Local tendency of the mean tangential wind

; RHS (1) : -{u} * {f+vort}  - Mean radial influx of absolute vertical vorticity
; RHS (2) : -{w} * (d{v}/dz) - Mean vertical advection of mean tangential momentum
; RHS (3) : -{u' * vort'}    - Eddy radial vorticity flux
; RHS (4) : -{w' * (dv'/dz)} - Vertical advection of eddy tangential momentum
; RHS (5) : {D_v}            - Combined diffusive and PBL tendency

;========================================================================================
;========================================================================================

; Arrays to hold all 'xpos' and 'ypos' position arrays
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
  avo_int0@units       = "s~S~-1~N~"
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
  ptb_vez0  = new((/num_l,dimsizes(radii),numTIMES/),float)                  ; [R3]
  ptb_vev0  = new((/num_l,dimsizes(radii),numTIMES/),float)                  ; [R4]

;==================================================
; Arrays for horizontal diffusive tendency terms
;==================================================

; Horizontal diffusivity (K_mh)
  k_mh0     = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)

; Arrays for fractional expressions in diffusive tendency source terms below (4D)
; (1) du/dy, (2) d(v/r)/dr, (3) u/r, (4) r, (5) d/dr, (6) d/dy
  dudy0     = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (1)
  dvdr0     = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (2)
  v_r0      = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (3)
  r_0       = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (4)
  ddr0      = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (5)
  ddy0      = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (6)

; Additional arrays for terms outside main loop
  du	    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ;
  dy	    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ;
  dr	    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ;
  r1	    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ;

; Stress tensor (tau_ry)
; { K_mh * ( (1/r * du/dy) + d(v/r)/dr ) }
; tau_ry1   = k_mh0 * ( ( (1 / r_0) * dudy0) + ( v_r0 * ddr0) )
  tau_ry1   = new((/num_l,dimsizes(radii),numTIMES/),float)                  ; [3D array]

; Total horizontal diffusive tendency of tangential momentum (using 'tau_ry')
; ( (1/(r^2) * {rho}) * d( r ^ 2 * {rho} * {tau_ry} ) /dr)
; v_dh0     = ( ( 1/( r_0 ^ 2 * mean_rho0) ) *
;             (ddr0 * ( r_0 ^ 2 * mean_rho0 * tau_ry1) ) )
  v_dh0     = new((/num_l,dimsizes(radii),numTIMES/),float)                  ; [3D array]

; Total vertical diffusive tendency of tangential momentum
; v_dz0     = (1 /{rho0}) * ( (d/dz) * ({rho0} * {tau_yz}) )
  v_dz0     = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [4D array]

;========================================
; Define 4D arrays for mean quantities
;========================================

; Mean tangential wind {v}
  mean_vtan0       = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  mean_vtan0!0     = "lev"
  mean_vtan0!1     = "rad"
  mean_vtan0!2	   = "azi"
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

  print_clock("Working on ensemble member: "+ens0(st))

; Set counter variable and size of arrays 
  ct          = ts0(st)
  llbox	      = toint(dist*50)
  dsize	      = (/llbox,llbox/)
  d0          = dsize(0)
  d1	      = dsize(1)

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

;==========================================================  
; Read in required variables from 'pc' and 'pd' streams
;========================================================== 

; Create array to hold storm track positions using different methods  
  centre_new  = new((/7,num_l,numTIMES,2/),"float")

  c  = addfile(fili_c,"r")         	 ; Read in 'pc' stream [u,v,vort]
  d  = addfile(fili_d,"r")         	 ; Read in 'pd' stream [w,p,t,theta,kmh,tau]

;==================================
; Get the variables we will need
;==================================

  do it = ts0(st), tf0(st)		 ; Loop over times in file

    print("Working on time: "+time_str(it))
    title_arr(it) = "Valid at "+time_str(it)+" (T+"+it+")"

;=====================================================
; Read in variables from 'pc' stream (model levels)
;=====================================================

    t0    = lat_arr(cn,it-1,0)
    t1    = lat_arr(cn,it-1,1)
    n0    = lon_arr(cn,it-1,0)
    n1    = lon_arr(cn,it-1,1)
    print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

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

;==========================================================
; Calculate divergence using one-sided difference scheme
;==========================================================

    div = uv2dv_cfd(u, v, latC, lonC, 2)
    copy_VarCoords(u, div)
    div@units = "s~S~-1~N~"
    div@name = "Divergence on pressure levels"

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

  ; Read in temperature, pressure, theta, vertical velocity and geopotential height
  ; All variables below -- [hybD | 63] * [latD | 150] * [lon1D | 151]
    t      = d->temp(it-2,:,{t0:t1},{n0:n1})        ; Temperature (K)
    p      = d->p(it-2,:,{t0:t1},{n0:n1})           ; Pressure (Pa)
    th     = d->theta(it-2,:,{t0:t1},{n0:n1})	    ; Potential temperature (K)
    w      = d->dz_dt(it-2,:,{t0:t1},{n0:n1})	    ; Vertical velocity (m s-1)
    z      = d->ht(it-2,:,{t0:t1},{n0:n1})    	    ; Geopotential height (m)
    pv     = d->field83(it-2,:,{t0:t1},{n0:n1})     ; Potential vorticity (PVU)

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
    p  = p / 100
    z  = z / 10
    pv = pv * (10 * 6)

  ; Add metadata
    p@units = "hPa"
    z@units = "dam"
    pv@units= "PVU"

    dx0     = lonD(1) - lonD(0) ; Grid spacing (longitude)
    dy0     = latD(1) - latD(0) ; Grid spacing (latitude)

;==========================
; Loop over model levels 
;==========================

  ; Loop over selected model levels
    do ilev = 40, 42, 2 ; num_l-7, 2

   ; Create output strings (model level info)
     hy  = sprintf("%0.0f",hybD(ilev) )
     hy1 = toint(hy)
     hy0 = sprintf("%05g",hy1) 
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

   ; Read in other variables on model levels
     u_plane 	            = u(ilev,:,:)	   ; Zonal wind 
     v_plane		    = v(ilev,:,:)	   ; Meridional wind 
     w_plane		    = w(ilev,:,:)	   ; Vertical velocity
     f_plane		    = f0(ilev,:,:)	   ; Coriolis parameter
     rho_plane              = rho(ilev,:,:)	   ; Density
     prs_plane              = p(ilev,:,:)	   ; Pressure
     str_plane              = str0(ilev,:,:)       ; Strain rate
     vort_plane	   	    = vort(ilev,:,:)	   ; Relative vorticity
     avo_plane		    = avo(ilev,:,:)	   ; Absolute vorticity
     kmh_plane		    = kmh(ilev,:,:)	   ; Horizontal diffusivity (K_mh)
     div_plane              = div(ilev,:,:)        ; Divergence 
     pv_plane               = pv(ilev,:,:)         ; Potential vorticity 

   ; Calculate vector windspeed and smooth 
     spd_plane              = sqrt( (u_plane ^ 2) + (v_plane ^ 2) )
     copy_VarCoords(u_plane, spd_plane)
     vort_pl_smth           = smth9_Wrap(vort_plane, 0.5, 0.5, True)

   ; Combine all arrays above into single, larger arrays for input into 'setup_cyl'

   ; 'pc' stream 
   ; full = 0 [u, v, spd, vort, avo, str, slp, div]
   ; full = 1 [u, v, spd, vort, avo, f, str, div]
     pc_plane               = (/u_plane, v_plane, spd_plane, vort_plane, \
                                avo_plane, str_plane, str_plane, div_plane/)

   ; 'pd' stream 
   ; full = 0 [w_plane, prs_plane, pv_plane, flux_plane]
   ; full = 1 [w_plane, prs_plane, ..., ..., kmh_plane, difv_xplane, difv_yplane, rho_plane]
     pd_plane               = (/w_plane, prs_plane, pv_plane, w_plane/)
     
;=======================================================
; Calculate storm centre position (external function)
;=======================================================

   ; Set radial distance (º) when looking for speed/vort min (in ringlike phase) 
     r0            = 0.12

     print_clock("Calculating storm centre position...")
     centre_arr    = st_centre(ilev, it, cn, r0, centre, vort_plane, prs_plane, spd_plane, \
                               centre_new, dy0, dx0, "comp", ct, "pd", "off")

   ; Updated storm track information contained in 'centre_new'
     centre_new    = centre_arr[0] 

;======================================================================
; Make sure that the storm 'centre' is not within the eyewall updraft 
;======================================================================

   ; If vertical velocity is above a threshold value at the storm centre, try again

   ; Retrieve coordinates of storm centre 
     cen_loc0 = centre_new(sc,ilev,ct,0)
     cen_loc1 = centre_new(sc,ilev,ct,1)

   ; Coordinates of grid points either side 
     cen_0n = cen_loc0 - dy0
     cen_0x = cen_loc0 + dy0
     cen_1n = cen_loc1 - dy0
     cen_1x = cen_loc1 + dy0

   ; Calculate vertical velocity and divergence at these grid points 
     cen_w    = w_plane( {cen_0n:cen_0x}, {cen_1n:cen_1x} )
     cen_d    = div_plane( {cen_0n:cen_0x}, {cen_1n:cen_1x} )

   ; If average over several grid points is above a threshold value, continue searching for centre
     ave_w    = avg(cen_w)
     ave_d    = avg(cen_d)
     print("ave. vertical velocity = "+ave_w+" ; ave. divergence = "+ave_d)

   ; Coordinates of grid points either side (bigger grid than above)
     delete([/cen_0n, cen_0x, cen_1n, cen_1x/])
     cen_0n = cen_loc0 - (2 * dy0)
     cen_0x = cen_loc0 + (2 * dy0)
     cen_1n = cen_loc1 - (2 * dy0)
     cen_1x = cen_loc1 + (2 * dy0)

     if (ave_w .gt. 0.1) then 

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
      delete([/ltN, lnN/])

    ; Print new storm centre to screen
      print("Centre (w_min): "+lat_min_w+" degrees N, "+lon_min_w+" degrees E")

    ; Replace values in 'centre_new' before calling 'setup_cyl' below 
      centre_new(sc,ilev,ct,0) = lat_min_w
      centre_new(sc,ilev,ct,1) = lon_min_w

     end if 

   ; Tidy up
     delete([/cen_w, cen_d/])

;========================================================
; Switch to cylindrical coordinates (external function)
;========================================================

   ; Last option = 0 --> basic variables output 
   ; Last option = 1 --> additional MBGT variables output  
     print_clock("Translating to cylindrical grid...")
     cyl_arr  = setup_cyl(it, dist, ll_arr, lonC, latC, lonD, latD, lon1D, lat1D, \
                          pc_plane, pd_plane, numTIMES, centre_new, \
                          sc, sm, ar, nr, ilev, num_l, 0, "comp", ct, \
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

     vr_int   = cyl_arr[13]
     vt_int   = cyl_arr[14]
     w_int    = cyl_arr[15]
     prs_int  = cyl_arr[16]

     xpos     = cyl_arr[17]
     ypos     = cyl_arr[18]
     thetas   = cyl_arr[19]
     xcen     = cyl_arr[20]
     ycen     = cyl_arr[21]

     str_int  = cyl_arr[22]
     pv_int   = cyl_arr[23]
     flux_int = cyl_arr[24]
     slp_int  = cyl_arr[25]
     div_int  = cyl_arr[26]

     print("Storm centre = "+lat_max+" ; "+lon_max) 

   ; Plot the data as we go  
     output   = "$sam/nepartak/images/cyl_coords/plt_xy_"+var0+"_"+dat(st)+"_"+ens0(st)+\
                "_"+typ0+"_"+hy0+"m_sc"+sc+"_"+time_arr(it)
     wks = gsn_open_wks(opt, output)

     gsn_define_colormap(wks,"prcp_new")

   ; Relative vorticity resources  
     opts_rv                               = True
     opts_rv@cnFillOn                      = True
     opts_rv@sfXArray                      = xpos_all(ilev,:,:,ct)
     opts_rv@sfYArray                      = ypos_all(ilev,:,:,ct)
     opts_rv@cnLineLabelInterval           = 2.0
     opts_rv@cnLineLabelFontHeightF        = 0.012
     opts_rv@cnLineLabelBackgroundColor    = "transparent"
     opts_rv@cnLineLabelPlacementMode      = "constant"
     opts_rv@cnLinesOn                     = False
     opts_rv@cnInfoLabelOn                 = False
     opts_rv@cnLevelSelectionMode          = "ExplicitLevels"
     opts_rv@cnLevels                      = (/5., 10., 15., 20., 25., 30., \
                                                40., 50., 60., 70., 80./)
     opts_rv@cnFillColors                  = (/0,4,5,6,7,8,9,10,11,12,13,14/)
     opts_rv@gsnMaximize                   = False

   ; No title or additional strings 
     opts_rv@tiMainString                  = ""
     opts_rv@gsnLeftString                 = ""
     opts_rv@gsnRightString                = ""

   ; Turn individual label bar on/off
     opts_rv@lbLabelBarOn                  = True    
     opts_rv@lbBoxEndCapStyle              = "TriangleBothEnds"
     opts_rv@lbOrientation                 = "Vertical"

     opts_rv@trYMaxF                       = lat_max + ar
     opts_rv@trYMinF                       = lat_max - ar 
     opts_rv@trXMaxF                       = lon_max + ar 
     opts_rv@trXMinF                       = lon_max - ar

   ; Do not draw the plot or advance the frame 
     opts_rv@gsnDraw                       = False
     opts_rv@gsnFrame                      = False
     opts_rv@gsnAddCyclic                  = False

   ; Background plot resources 
     opts_xy                                = True

   ; No tickmarks on any axis
     opts_xy@tmYROn                        = False
     opts_xy@tmYLOn                        = False
     opts_xy@tmXTOn                        = False
     opts_xy@tmXBOn                        = False

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
     opts_vvel@cnLevels                    = (/-1.5, -1.0, -0.5, -0.2, -0.1, -0.05, \
                                                0.05, 0.1, 0.2, 0.5, 1.0, 1.5/)
     opts_vvel@cnFillColors                = (/1,2,3,4,6,8,-1,9,11,13,14,15,16/)

   ; Resources for strain rate  
     opts_strn                              = opts_vvel
     delete([/opts_strn@cnLevels, opts_strn@cnFillColors, opts_strn@cnFillPalette/])
     opts_strn@cnFillPalette                = "circ_budget"
     opts_strn@cnLevels                     = (/10,12,14,16,18,20,22,25,30/)
     opts_strn@cnFillColors                 = (/-1,11,13,15,16,17,18,19,20,21/)

   ; Resources for divergence 
     opts_divg                              = opts_vvel
     delete([/opts_divg@cnLevels, opts_divg@cnFillColors, opts_divg@cnFillPalette/])
     if (var0 .eq. "div0") then 
      opts_divg@cnFillPalette               = "prcp_new"
      opts_divg@cnLevels                    = (/5,10,15,20,30,40,60,80,100,150/)
      opts_divg@cnFillColors                = (/-1,2,3,4,5,6,7,8,9,10,11,12,13/)
     elseif (var0 .eq. "div1") then 
      opts_divg@cnFillPalette               = "BlueDarkRed18"
      opts_divg@cnLevels                    = (/-100, -50, -30, -20, -10, -5, \
                                                5, 10, 20, 30, 50, 100/)
      opts_divg@cnFillColors                = (/1,2,3,4,6,8,-1,9,11,13,14,15,16/)
     end if 

   ; Panel plot label resources 
     txid                                   = new(4,graphic)
     amid                                   = new(4,graphic)

     if (var0 .eq. "str") then 
      panel_strings                         = (/ "a) Relative vorticity",\
                                                 "b) Vertical velocity",\
                                                 "c) Tangential wind",\
                                                 "d) Strain rate"/)
     else
      panel_strings                         = (/ "a) Relative vorticity",\
                                                 "b) Vertical velocity",\
                                                 "c) Tangential wind",\
                                                 "d) Divergence"/) 
     end if 

     txres                                  = True
     txres@txPerimOn                        = True
     txres@txFontHeightF                    = 0.0125
     txres@txBackgroundFillColor            = "White"

     amres                                  = True
     amres@amParallelPosF                   = -0.48
     amres@amOrthogonalPosF                 = -0.48
     amres@amJust                           = "TopLeft"

   ; Coordinates of panels 
      xf = (/0.05, 0.52, 0.05, 0.52/)
      yf = (/0.90, 0.90, 0.49, 0.49/)
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

     plot_tl            = gsn_csm_xy(wks, xpos_all(ilev,:,:,ct), ypos_all(ilev,:,:,ct), opts_xy0)
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

     plot_tr              = gsn_csm_xy(wks, xpos_all(ilev,:,:,ct), ypos_all(ilev,:,:,ct), opts_xy1)
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

     plot_bl              = gsn_csm_xy(wks, xpos_all(ilev,:,:,ct), ypos_all(ilev,:,:,ct), opts_xy2)
     plot_bl0             = gsn_csm_contour(wks, vt_int(:,:), opts_wind)
     overlay(plot_bl, plot_bl0)

   ; Overlay labels  
     txid(2)      = gsn_create_text(wks, panel_strings(2), txres)
     amid(2)      = gsn_add_annotation(plot_bl, txid(2), amres)

   ; Panel 4 (strain rate / divergence)
     if (var0 .eq. "str") then 
      opts_strn0          = opts_strn
     else
      opts_strn0          = opts_divg
     end if 

     opts_strn0@vpXF      = xf(3)
     opts_strn0@vpYF      = yf(3)
     opts_strn0@vpWidthF  = wf
     opts_strn0@vpHeightF = hf

     opts_xy3             = opts_xy
     opts_xy3@vpXF        = xf(3)
     opts_xy3@vpYF        = yf(3)
     opts_xy3@vpWidthF    = wf
     opts_xy3@vpHeightF   = hf

   ; Change units before plotting 
     str_int_plot         = str_int * (10 ^ 4)
     div_int_plot         = div_int * (10 ^ 5)

   ; Plot strain rate/divergence and overlay cylindrical grid 
     plot_br              = gsn_csm_xy(wks, xpos_all(ilev,:,:,ct), ypos_all(ilev,:,:,ct), opts_xy3)
     if (var0 .eq. "str") then 
      plot_br0            = gsn_csm_contour(wks, str_int_plot(:,:), opts_strn0)
     else 
      plot_br0            = gsn_csm_contour(wks, div_int_plot(:,:), opts_strn0)
     end if 
     overlay(plot_br, plot_br0)

   ; Overlay labels
     txid(3)      = gsn_create_text(wks, panel_strings(3), txres)
     amid(3)      = gsn_add_annotation(plot_br, txid(3), amres)

   ; Tidy up before final plots later in script
     delete([/opts_vvel, opts_wind, opts_rv, opts_strn, opts_divg, \
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
    markers(0,i)             = gsn_add_polymarker(wks, plot_tl, centre_new(i,ilev,ct,1),\
                                                  centre_new(i,ilev,ct,0), mres0)
    markers(1,i)             = gsn_add_polymarker(wks, plot_tr, centre_new(i,ilev,ct,1),\
                                                  centre_new(i,ilev,ct,0), mres0)
    markers(2,i)             = gsn_add_polymarker(wks, plot_bl, centre_new(i,ilev,ct,1),\
                                                  centre_new(i,ilev,ct,0), mres0)
    markers(3,i)             = gsn_add_polymarker(wks, plot_br, centre_new(i,ilev,ct,1),\
                                                  centre_new(i,ilev,ct,0), mres0)
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
      title_str                   = "Panel plot at "+hy+" m AGL: "+ens0(st)+\
                                    ", "+dat(st)+" (T+"+it+")"
      title_res@txFontHeightF     = 0.02
      title_res@txFont            = "Helvetica"
      xpos_title                  = 0.50
      ypos_title                  = 0.94

      gsn_text_ndc(wks, title_str, xpos_title, ypos_title, title_res)

    ; Call 'draw' and 'frame'
      maximize_output(wks, True)

  ; Tidy up
    delete([/u_int,v_int,spd_int,w_int,vt_int,vr_int,div_int/])
    delete([/prs_int,str_int,pv_int,flux_int,vort_int,avo_int,slp_int/])
    delete([/u_plane,v_plane,w_plane,f_plane,vort_plane,avo_plane/])
    delete([/difv_xplane,difv_yplane,kmh_plane,pv_plane/])

   end do                           ; End model level loop (do ilev = 0, num_l-1)
   
   ct = ct + 1 

 ; Tidy up before next iteration (all variables from 'pc' and 'pd' streams)
   delete([/lonC,latC,lon1C,lat1C,u,v,vort,avo,f0,f1,div/])
   delete([/lonD,latD,lon1D,lat1D,str_xz,str_yz,kmh,t,p,w,z,rho,rho_a/])

  end do     ; End loop over times (do it = ts0(st), tf0(st) )

 end do      ; End loop over sets of dates (do st = 0, dimsizes(dat)-1)

end 