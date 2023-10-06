; Calculate composite plots of relative vorticity and strain in different VC phases

; Run using:

; ncl opt=\"x11\" dist=2.0 nt=73 nr=11 ar=0.5 typ=\"sym\" full=9 cn0=\"geo_sm\"
; mlev0=0 mlev1=28 sc=5 sm=0 rad0=0.5 diff=1 azi=1 var=\"w\" tc_mbgt_scatter.ncl

; 'opt'   = output file format ("pdf" or "x11")
; 'dist'  = size of box following storm (degrees)
; 'typ'	  = VC phase ("sym","asym","stoa","atos")
; 'full'  = top 18 (1); top 8 (-5); top 6 (-6)
; 'ar'    = distance in degrees of outer radius (1.0, 2.0, etc)
; 'nr'    = number of radial circles between r = 0 and r = 'ar' (21, 41, 61, etc)
; 'nt'    = number of azimuth angles in cylindrical grid
; 'cn0'   = offline storm centre: "slp", "vort", "geo", "geo_sm"
; 'rad0'  = radius (degrees) for final plots 
; 'diff'  = calculate difference between two levels (1) 
; 'azi'   = calculate azimuthal average (1) or plot all points (0)

; 500 m [10], 1 km [14], 1.5 km [17], 2 km [20], 3 km [24], 4 km [28], 5 km [31], 6 km [34]
; 7 km [37], 8 km [39], 9 km [42], 10 km [45], 11 km [46], 12 km [48], 13 km [50], 14 km [52]  
; 15 km [54], 16 km [56]

; 'pc' stream --> 'it'
; 'centre'    --> 'it-1'
; 'pd' stream --> 'it-2'

; Load main NCL functions and procedures
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"

; Also load user-defined functions and procedures
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_centre.ncl"
load "$sam/ncl_func/setup_cyl.ncl"
load "$sam/ncl_func/draw_cyl.ncl"

begin

;==============================
; Prelash (reading in files)
;==============================

; Dimensions of text files to read data in from (S-phase, S-A, A-phase, A-S) 
; EDIT -- make clear how these groups correspond to those in 'cxz' scripts
; 'vc_arr(:,0)' = Group 1 
; 'vc_arr(:,1)' = Group 2 
; 'vc_arr(:,2)' = Group 3 
  vc_arr_all = (/ (/17,25,31,29,24/), (/16,14,30,24,24/), (/16,29,31,38,24/), (/9,15,31,24,24/) /)

; String containing info on width of radial averaging band  
  ar0  = ar * 100
  rstr = sprintf("%0.0fkm", ar0)

; String for output file 
; EDIT -- remember discussion with Montgomery and Kilroy about terminology (avoid 'sym'/'asym')
  if (typ .eq. "sym") then 
   out_str = "ringlike"
  elseif (typ .eq. "s2a") then 
   out_str = "ring->mono transition"
  elseif (typ .eq. "asym") then 
   out_str = "monopole"
  elseif (typ .eq. "a2s") then 
   out_str = "mono->ring transition"
  end if 

  if (var .eq. "w") then 
   var_str = "Vertical velocity"
  elseif (var .eq. "vr") then 
   var_str = "Radial wind"
  elseif (var .eq. "vort") then
   var_str = "Relative vorticity"
  elseif (var .eq. "flux") then 
   var_str = "Eddy radial vorticity flux"
  end if 

; Troubleshooting
  if (full .lt. 9) then 
   print_clock("Not using new method of calculation - sure you want to continue?")  
  end if 

; For now, explicitly tell NCL which dates and simulations to analyse
  if (typ .eq. "sym") then

   if (full .eq. 8) then

    dat  = asciiread("$ar/text/vc_sym_dat_group4.txt",(/vc_arr_all(0,3)/),"string")
    ens0 = asciiread("$ar/text/vc_sym_sim_group4.txt",(/vc_arr_all(0,3)/),"string")
    ts0  = asciiread("$ar/text/vc_sym_ts_group4.txt",(/vc_arr_all(0,3)/),"integer")
    tf0  = ts0 + 1

   elseif (full .eq. 9) then

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

  elseif (typ .eq. "s2a") then

   if (full .eq. 8) then

    dat  = asciiread("$ar/text/vc_s2a_dat_group4.txt",(/vc_arr_all(1,3)/),"string")
    ens0 = asciiread("$ar/text/vc_s2a_sim_group4.txt",(/vc_arr_all(1,3)/),"string")
    ts0  = asciiread("$ar/text/vc_s2a_ts_group4.txt",(/vc_arr_all(1,3)/),"integer")
    tf0  = ts0 + 1

   elseif (full .eq. 9) then

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

   if (full .eq. 8) then

    dat  = asciiread("$ar/text/vc_asym_dat_group4.txt",(/vc_arr_all(2,3)/),"string")
    ens0 = asciiread("$ar/text/vc_asym_sim_group4.txt",(/vc_arr_all(2,3)/),"string")
    ts0  = asciiread("$ar/text/vc_asym_ts_group4.txt",(/vc_arr_all(2,3)/),"integer")
    tf0  = ts0 + 1

   elseif (full .eq. 9) then

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

  elseif (typ .eq. "a2s") then

   if (full .eq. 8) then

    dat  = asciiread("$ar/text/vc_a2s_dat_group4.txt",(/vc_arr_all(3,3)/),"string")
    ens0 = asciiread("$ar/text/vc_a2s_sim_group4.txt",(/vc_arr_all(3,3)/),"string")
    ts0  = asciiread("$ar/text/vc_a2s_ts_group4.txt",(/vc_arr_all(3,3)/),"integer")
    tf0  = ts0 + 1

   elseif (full .eq. 9) then

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

;=====================================================
; Define cylindrical coordinate arrays for plotting 
;=====================================================

; Calculate total number of time intervals analysed
  num_int = sum( (tf0 - ts0) + 1)
  print("Total number of times in sample: "+num_int)

; Create new array to hold storm centre positions [7 * 59 * num_int * 2]
  num_l        = (mlev1 - mlev0) + 1
  centre_new   = new( (/7, num_l, num_int, 2/), float)

  pii        = 3.14159265
  radii      = fspan(0.0,ar,nr)      ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,nt)   ; Array of azimuth angles (0, 5, 10, etc)
  thetas_rad = (pii/180)*thetas      ; Degrees to radians
  rad_size   = dimsizes(radii)-1     ; Index for later calculation

; Define constants (gas constant for dry air; specific heat of dry air at constant p)
  rd    = 287.0 
  cp0   = 1004.0

; 4D arrays to hold 'xpos' and 'ypos' values for all analysis times and model levels
; Even if running for all 4 phases, don't need extra dimension (values will overwrite)
  xpos_all    = new((/num_l,dimsizes(radii),dimsizes(thetas),num_int/),float)
  ypos_all    = new((/num_l,dimsizes(radii),dimsizes(thetas),num_int/),float)
  xcen_all    = new((/num_l,dimsizes(radii),dimsizes(thetas),num_int/),float)
  ycen_all    = new((/num_l,dimsizes(radii),dimsizes(thetas),num_int/),float)

; 4D arrays for all variables (u,v,w,vort,spd,vt,vr)
  w_int0      = new((/num_l,dimsizes(radii),dimsizes(thetas),num_int/),float)
  vort_int0   = new((/num_l,dimsizes(radii),dimsizes(thetas),num_int/),float)
  spd_int0    = new((/num_l,dimsizes(radii),dimsizes(thetas),num_int/),float)
  vt_int0     = new((/num_l,dimsizes(radii),dimsizes(thetas),num_int/),float)
  vr_int0     = new((/num_l,dimsizes(radii),dimsizes(thetas),num_int/),float)

; 4D arrays for mean and perturbation variables, for all phases (w,u,vort)
  wp_int0     = new((/num_l,dimsizes(radii),dimsizes(thetas),num_int/),float)
  wm_int0     = new((/num_l,dimsizes(radii),dimsizes(thetas),num_int/),float)
  up_int0     = new((/num_l,dimsizes(radii),dimsizes(thetas),num_int/),float)
  um_int0     = new((/num_l,dimsizes(radii),dimsizes(thetas),num_int/),float)
  zp_int0     = new((/num_l,dimsizes(radii),dimsizes(thetas),num_int/),float)
  zm_int0     = new((/num_l,dimsizes(radii),dimsizes(thetas),num_int/),float)
  rf_int0     = new((/num_l,dimsizes(radii),dimsizes(thetas),num_int/),float)

  wp_int0@description  = "Perturbation vertical velocity"
  wp_int0@units        = "m s~S~-1~N~"
  wp_int0!0            = "lev"
  wp_int0!1            = "rad"
  wp_int0!2            = "azi"
  wp_int0!3            = "time"

  wm_int0@description  = "Mean vertical velocity"
  wm_int0@units        = "m s~S~-1~N~"
  wm_int0!0            = "lev"
  wm_int0!1            = "rad"
  wm_int0!2            = "azi"
  wm_int0!3            = "time"

  up_int0@description  = "Perturbation radial wind"
  up_int0@units        = "m s~S~-1~N~"
  up_int0!0            = "lev"
  up_int0!1            = "rad"
  up_int0!2            = "azi"
  up_int0!3            = "time"

  um_int0@description  = "Mean radial wind"
  um_int0@units        = "m s~S~-1~N~"
  um_int0!0            = "lev"
  um_int0!1            = "rad"
  um_int0!2            = "azi"
  um_int0!3            = "time"

  zp_int0@description  = "Perturbation relative vorticity"
  zp_int0@units        = "s~S~-1~N~"
  zp_int0!0            = "lev"
  zp_int0!1            = "rad"
  zp_int0!2            = "azi"
  zp_int0!3            = "time"

  zm_int0@description  = "Mean relative vorticity"
  zm_int0@units        = "s~S~-1~N~"
  zm_int0!0            = "lev"
  zm_int0!1            = "rad"
  zm_int0!2            = "azi"
  zm_int0!3            = "time"

  rf_int0@description  = "Eddy radial vorticity flux"
  rf_int0@units        = "s~S~-1~N~"
  rf_int0!0            = "lev"
  rf_int0!1            = "rad"
  rf_int0!2            = "azi"
  rf_int0!3            = "time"

  w_int0@description   = "Vertical velocity"
  w_int0@units         = "m s~S~-1~N~"
  w_int0!0             = "rad"
  w_int0!1             = "azi"
  w_int0!2             = "time"

  vort_int0@description = "Relative vorticity"
  vort_int0@units       = "s~S~-1~N~"
  vort_int0!0           = "rad"
  vort_int0!1           = "azi"
  vort_int0!2           = "time"

  spd_int0@description  = "Vector wind"
  spd_int0@units        = "m s~S~-1~N~"
  spd_int0!0            = "rad"
  spd_int0!1            = "azi"
  spd_int0!2            = "time"

  vt_int0@description   = "Tangential wind"
  vt_int0@units         = "m s~S~-1~N~"
  vt_int0!0             = "rad"
  vt_int0!1             = "azi"
  vt_int0!2             = "time"

; Initialise counter variable 
  ct                    = 0

; Print out relevant information about storm tracking (sc=0,1,2,3,4,5,6)
  sc_arr                = (/"relative vorticity",\
                            "pressure",\          
                            "smoothed pressure",\ 
                            "extra smoothed pressure",\
                            "vort_min during symmetric phase",\
                            "windspeed minimum",\
                            "geopotential"/)

  print_clock("Using "+sc_arr(sc)+" to calculate storm centre position...")

;==========================================================
; Now start loop over these chosen input files and times
;==========================================================

 do st = 0, dimsizes(dat)-1 

  dt0 = tf0(st) - ts0(st)

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
  numINPUT_c  = numINFO_c[0]		; Number of times in 'pc' stream
  time_c      = numINFO_c[1]		; Array of times from 'pc' stream file
  numINPUT_d  = numINFO_d[0]		; Number of times in 'pd' stream
  time_d      = numINFO_d[1]		; Array of times from 'pd' stream file 

;============================================================
; Calculate storm motion using built-in function ('st_rm')
;============================================================

; 'dat'      = initialisation time (02T12, 03T00, ...)
; 'ens0'     = ensemble simulation (em00, em01, ...)
; 'diri'     = path to input files (see above)
; 'dist'     = size of box following storm (degrees)
; 'mins'     = analyse 1-h (0) or 5-min (1) data

  storm_rel  = st_rm(dat(st), ens0(st), "$ar/text/", dist, 0)
  u_cyc      = storm_rel[0]		   ; Zonal wind 
  v_cyc      = storm_rel[1]		   ; Meridional wind
  vel_cyc    = storm_rel[2]		   ; Vector wind 
  lat_arr    = storm_rel[4]                ; Storm-following domain subset (lat)
  lon_arr    = storm_rel[5]                ; Storm-following domain subset (lon)
  centre     = storm_rel[6]		   ; Storm centre at each hour (lat/lon)

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

;=========================
; Start multiple loops
;=========================

   print_clock("Working on ensemble member: "+ens0(st)+" ("+dat(st)+")" )

   llbox       = toint(dist*50)	           ; Size of lat/lon arrays (see line below)
   dsize       = (/llbox,llbox/)           ; Size of lat/lon arrays (depends on 'dist')
   d0	       = dsize(0)
   d1	       = dsize(1)

 ; Array to hold values of inner-core parameter
   sym_arr     = new(numINPUT_d,float)

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

;===============================
; Start loop over input files 
;===============================

 ; Read in 'pc' [u,v,vort] and 'pd' streams [w,p,t,theta,kmh,tau]
   c  = addfile(fili_c,"r")
   d  = addfile(fili_d,"r")

;==================================
; Get the variables we will need
;==================================

 ; Loop over times in file
   do it = ts0(st), tf0(st)

    print("Working on time: "+time_str(it)+" (ct = "+ct+")" )
    itt = it
    title_arr(it) = "Valid at "+time_str(it)+" (T+"+itt+")"

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
    hyb_m  = sprintf("%0.0f",hybC(:))

  ; Read in horizontal wind components 
    u_varname = "x-wind"
    v_varname = "y-wind"
    u0        = c->$u_varname$(it,mlev0:mlev1,{t0:t1},{n0:n1})
    v0        = c->$v_varname$(it,mlev0:mlev1,{t0:t1},{n0:n1})

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

  ; Calculate relative vorticity using centered finite differences
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

  ; Calculate horizontal [lat,lon] derivatives of 'u' (du_dx, du_dy)
    du0    = grad_latlon_cfd(u,u&$dim1$,u&$dim2$,False,False)
    du_dy  = du0[0]
    du_dx  = du0[1]
    delete(du0)

  ; Calculate horizontal derivatives of 'v' (dv_dx, dv_dy) 
    dv0    = grad_latlon_cfd(v,v&$dim1$,v&$dim2$,False,False)
    dv_dy  = dv0[0]
    dv_dx  = dv0[1]
    delete(dv0)

  ; Calculate strain rate components --> sqrt[ (str)^2 + (shr)^2 ]
  ; Stretching (du_dx - dv_dy) and shearing (dv_dx + du_dy) deformation
    str1a   = du_dx - dv_dy
    str2a   = dv_dx + du_dy
    str1    = str1a ^ 2
    str2    = str2a ^ 2

  ; Strain rate
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

  ; Read in temperature (K), pressure (Pa), vertical velocity (m/s) and geopotential height (m)
    t      = d->temp(it-2,mlev0:mlev1,{t0:t1},{n0:n1})
    p      = d->p(it-2,mlev0:mlev1,{t0:t1},{n0:n1}) 
    w      = d->dz_dt(it-2,mlev0:mlev1,{t0:t1},{n0:n1}) 
    z      = d->ht(it-2,mlev0:mlev1,{t0:t1},{n0:n1})  

  ; Calculate dry air density from pressure and temperature
    rho    = p / (rd * t)

  ; Add metadata
    copy_VarCoords(p, rho)
    rho@description = "Dry air density"
    rho@units = "kg m~S~-3~N~"

  ; Convert units of pressure and geopotential height before continuing 
    p      = p / 100
    z      = z / 10

  ; Add metadata
    p@units = "hPa"
    z@units = "dam"

  ; Grid spacing (lon, lat)
    dx0    = lonD(1) - lonD(0)
    dy0    = latD(1) - latD(0)

;==========================
; Loop over model levels 
;==========================

    do ilev = 0, num_l-1 

     hy = sprintf("%0.0f",hybD(ilev))
     print_clock("Working on model level " + ilev+ " ("+hy+" m AGL)")

   ; Read in other variables on model levels (u,v,w,f,p,str,rho,vort,avo,kmh)
     u_plane                = u(ilev,:,:)
     v_plane                = v(ilev,:,:)
     w_plane                = w(ilev,:,:)
     prs_plane              = p(ilev,:,:)
     str_plane              = str0(ilev,:,:)
     vort_plane             = vort(ilev,:,:)
     avo_plane              = avo(ilev,:,:)

   ; Calculate vector windspeed and smooth 
     spd_plane              = sqrt( (u_plane ^ 2) + (v_plane ^ 2) )
     copy_VarCoords(u_plane, spd_plane)
     vort_pl_smth           = smth9_Wrap(vort_plane, 0.5, 0.5, True)

   ; Combine all arrays above into single, larger array for input into 'setup_cyl' 
     pc_plane               = (/u_plane, v_plane, spd_plane, \
                                vort_plane, avo_plane, str_plane/)

   ; 2nd/3rd 'prs_plane' are both placeholders for PV/mflux, which we haven't calculated  
     pd_plane               = (/w_plane, prs_plane, prs_plane, prs_plane/)

;=============================================================  
; Call external function to calculate storm centre position
;============================================================= 

    ; Set radial distance (º) when looking for vort min in ringlike phase   
      r0            = 0.10 

    ; Lower level/layer
      centre_arr    = st_centre(0, it, cn, r0, centre, vort_plane, prs_plane, spd_plane, \
                                centre_new, dy0, dx0, "comp", ct, "pd", "off")

    ; Updated storm track information contained in 'centre_new'
      centre_new  = centre_arr[0]

;=====================================================================  
; Call external function to switch to cylindrical coordinate system  
;=====================================================================  

    ; Call 'setup_cyl' and output variables on cylindrical grid 
    ; Option before "comp" --> basic ('0') or MBGT ('1') variables output    
    ; Option before 'num_l' --> 'ilev=0' 
    ; Option before 'xpos_all' etc --> 'lay=0/1' (layer or individual levels)
      cyl_arr    = setup_cyl(ct, dist, ll_arr, lonC, latC, lonD, latD, lon1D, lat1D, \
                             pc_plane, pd_plane, num_int, centre_new, \
                             sc, sm, ar, nr, 0, num_l, 0, "comp", ct, \
                             0, xpos_all, ypos_all, xcen_all, ycen_all)

      xpos_all   = cyl_arr[0]
      ypos_all   = cyl_arr[1]
      xcen_all   = cyl_arr[2]
      ycen_all   = cyl_arr[3] 

      lat_max    = cyl_arr[4]
      lon_max    = cyl_arr[5]
      radii      = cyl_arr[6]
      rad_size   = cyl_arr[7]
      thetas     = cyl_arr[19] 

      spd_int    = cyl_arr[10]
      vort_int   = cyl_arr[11]
      vr_int     = cyl_arr[13]
      w_int      = cyl_arr[15]

    ; Calculate azimuthal averages
      do irad = 0, dimsizes(radii)-1
       wm_int0(ilev,irad,:,ct) = avg(w_int(irad,:) )
       um_int0(ilev,irad,:,ct) = avg(vr_int(irad,:) )
       zm_int0(ilev,irad,:,ct) = avg(vort_int(irad,:) )
      end do 

    ; Fill larger arrays with key variables (vr,w,vort)
      vr_int0(ilev,:,:,ct)     = (/vr_int/) 
      w_int0(ilev,:,:,ct)      = (/w_int/) 
      vort_int0(ilev,:,:,ct)   = (/vort_int/)

  ; Tidy up
    delete([/spd_int, vort_int, vr_int, w_int/])
    delete([/w_plane, vort_plane, u_plane, v_plane, spd_plane, str_plane, avo_plane, prs_plane/])

    end do                          ; End model level loop (do ilev = 0, num_l-1)

    ct  = ct + 1		    ; Counter variable (time)

 ; Tidy up before next iteration (all variables from 'pc' and 'pd' streams)
   delete([/lonC, latC, lon1C, lat1C, u0, v0, u, v, vort/])
   delete([/du_dy, du_dx, dv_dy, dv_dx, str0/])
   delete([/lonD, latD, lon1D, lat1D, w, p, t, rho/])

  end do     ; End time loop (do it = 0, times-1)

 end do      ; End loop over input file times (do st = 0, dimsizes(dat)-1

;========================
; Calculate eddy terms
;========================

; Vertical velocity, radial wind, relative vorticity and radial vorticity flux
  wp_int0    = w_int0(:,:,:,:) - wm_int0(:,:,:,:)
  up_int0    = vr_int0(:,:,:,:) - um_int0(:,:,:,:)
  zp_int0    = vort_int0(:,:,:,:) - zm_int0(:,:,:,:)
  rf_int0    = up_int0 * zp_int0

;====================================
; Now produce scatter or box plots
;====================================
 
 ; Plotting resources for perturbation variables 
   pt_opts                      = True 
   pt_opts@gsnDraw              = False
   pt_opts@gsnFrame             = False
   pt_opts@gsnMaximize          = False 

 ; Y-axis
   pt_opts@tiMainString         = var_str+" ("+out_str+")"
   pt_opts@tiYAxisString        = ""
   pt_opts@tiYAxisFontHeightF   = 0.0150

 ; EDIT -- how to reduce whitespace between 'columns'?

 ; X-axis (define our own labels)
   pt_opts@tiXAxisString        = "Radius (km)"
   pt_opts@tiXAxisFontHeightF   = 0.0150
   pt_opts@tmXBLabelFontHeightF = 0.0125
   pt_opts@tmXBMode             = "Explicit"
   pt_opts@tmXBValues           = (/0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45/)
   pt_opts@tmXBLabels           = (/"0", "5", "10", "15", "20", "25", "30", "35", "40", "45"/)

;   pt_opts@vpWidthF             = 0.5
;   pt_opts@vpHeightF            = 0.7

 ; Additional resources 
   pt_opts@xyMarkLineMode       = "Markers"
   pt_opts@xyMonoMarkerColor    = True
   pt_opts@xyMarkerColor        = "RoyalBlue1"
   pt_opts@xyMonoMarkerSize     = True 
   pt_opts@xyMonoMarkerThickness= True 
   pt_opts@xyMarkerSizeF        = 0.002
   pt_opts@xyMarkerThicknessF   = 0.02

 ; Turn off tickmarks
   pt_opts@tmXTOn               = "False"
   pt_opts@tmYROn               = "False"

 ; Additional resources (x-axis, etc)
   pt_opts@trXMaxF              = radii(dimsizes(radii)-2)
   pt_opts@trXMinF              = 0.05

 ; Open workstation and define output file path
   output = "$sam/nepartak/images/vc/mbgt_scatter_"+var+"_full"+full+"_sc"+sc+"_"+typ
   wks    = gsn_open_wks(opt, output) 

 ; Switch dimension order before plotting (don't want all vertical levels)
   vvel_p_plot = wp_int0(time|:, lev|:, azi|:, rad|:)
   vrad_p_plot = up_int0(time|:, lev|:, azi|:, rad|:)
   vort_p_plot = zp_int0(time|:, lev|:, azi|:, rad|:)
   flux_p_plot = rf_int0(time|:, lev|:, azi|:, rad|:)

   if (var .eq. "w") then 
    plot_final0 = vvel_p_plot
   elseif (var .eq. "vr") then 
    plot_final0 = vrad_p_plot
   elseif (var .eq. "vort") then 
    plot_final0 = vort_p_plot
   elseif (var .eq. "flux") then 
    plot_final0 = flux_p_plot
   end if 

 ; Create object to hold all data for plotting 
   plot  = new( num_int, graphic)

 ; Reshape array for plotting (first make 1D, then 2D)
   plot_final1 = ndtooned(plot_final0)
   plot_final  = onedtond(plot_final1, (/101616, 11/) )

 ; Add metadat
   plot_final!1   = "rad"
   plot_final&rad = radii

;==========================
; Finally, plot the data
;==========================

 ; Plot first time as base
   plot = gsn_csm_xy(wks, radii, plot_final(:,:), pt_opts)

   draw(plot)
   frame(wks)

; ; Some additional plotting resources 
;   pt2_opts                     = pt_opts
;   pt2_opts@tiMainString        = ""
;   pt2_opts@tiYAxisString       = ""
;   pt2_opts@tiXAxisString       = ""
;   pt2_opts@tmXBLabels          = (/"", "", "", "", "", "", "", "", "", ""/)

;   pt2_opts@tmYLMode            = "Explicit"
;   pt2_opts@tmYLValues          = (/-4.0, -2.0, 0.0, 2.0, 4.0, 6.0/)
;   pt2_opts@tmYLLabels          = (/"", "", "", "", "", ""/)

;   do ct = 1, num_int-1
;    plot(ct) = gsn_csm_xy(wks, radii, plot_final(:,:), pt2_opts)
;   end do

end 