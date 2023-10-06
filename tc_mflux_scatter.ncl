; Produce a scatter plot of vertical mass flux for chosen VC phase 

; Run using:

; ncl opt=\"x11\" dist=2.0 nt=73 nr=11 ar=0.5 typ=\"sym\" full=9 cn0=\"geo_sm\"
; lay=1 mlev0=17 mlev1=34 int=2 sc=5 sm=0 rad0=0.5 diff=1 azi=1 tc_mflux_scatter.ncl

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
   out_str = "Symmetric phase"
  elseif (typ .eq. "s2a") then 
   out_str = "Sym to asym transition"
  elseif (typ .eq. "asym") then 
   out_str = "Asymmetric phase"
  elseif (typ .eq. "a2s") then 
   out_str = "Asym to sym transition"
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
  num_l        = 59 
  centre_new_l = new( (/7, 1, num_int, 2/), float)
  centre_new_u = new( (/7, 1, num_int, 2/), float)

  pii        = 3.14159265
  radii      = fspan(0.0,ar,nr)      ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,nt)   ; Array of azimuth angles (0, 5, 10, etc)
  thetas_rad = (pii/180)*thetas      ; Degrees to radians
  rad_size   = dimsizes(radii)-1     ; Index for later calculation

; Define constants (gas constant for dry air; specific heat of dry air at constant p)
  rd    = 287.0 
  cp0   = 1004.0

; 3D arrays to hold 'xpos' and 'ypos' values for all analysis times 
  xpos_u_all  = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)
  ypos_u_all  = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)
  xcen_u_all  = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)
  ycen_u_all  = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)

  xpos_l_all  = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)
  ypos_l_all  = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)
  xcen_l_all  = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)
  ycen_l_all  = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)

; 3D arrays for all variables
  u_int0      = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)  ; Zonal wind         
  v_int0      = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)  ; Meridional wind    
  w_int0      = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)  ; Vertical velocity  
  vort_int0   = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)  ; Relative vorticity

  spd_int0    = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)  ; Vector windspeed
  vt_int0     = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)  ; Tangential wind 

  rad_arr0    = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)  ; Radius (for v/r)
  flux_u_int0 = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)  ; Mass flux (upper)
  flux_l_int0 = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)  ; Mass flux (lower)
  flux_d_int0 = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)  ; Mass flux (difference)

  flux_u_int0@description = "Vertical mass flux"
  flux_u_int0@units       = "kg m~S~-2~N~ s~S~-1~N~"
  flux_u_int0!0           = "rad"
  flux_u_int0!1           = "azi"
  flux_u_int0!2           = "time"

  flux_l_int0@description = "Vertical mass flux"
  flux_l_int0@units       = "kg m~S~-2~N~ s~S~-1~N~"
  flux_l_int0!0           = "rad"
  flux_l_int0!1           = "azi"
  flux_l_int0!2           = "time"

  flux_d_int0@description = "Vertical mass flux"
  flux_d_int0@units       = "kg m~S~-2~N~ s~S~-1~N~"
  flux_d_int0!0           = "rad"
  flux_d_int0!1           = "azi"
  flux_d_int0!2           = "time"

  u_int0@description   = "Zonal wind"
  u_int0@units         = "m s~S~-1~N~"
  u_int0!0             = "rad"
  u_int0!1             = "azi"
  u_int0!2             = "time"

  v_int0@description   = "Meridional wind"
  v_int0@units         = "m s~S~-1~N~"
  v_int0!0             = "rad"
  v_int0!1             = "azi"
  v_int0!2             = "time"

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

  spd_int0@description  = "Vector windspeed"
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

  ; Calculate model level height information
    hy0    = sprintf("%0.0f",hybC(mlev0))
    hy1    = sprintf("%0.0f",hybC(mlev1))

  ; Read in horizontal wind components 
    u_varname = "x-wind"
    v_varname = "y-wind"

  ; Read in values on two single levels (lay=0) or over two layers (lay=1)
    if (lay .eq. 1) then

     ul0 = c->$u_varname$(it,mlev0-int:mlev0+int,{t0:t1},{n0:n1})
     vl0 = c->$v_varname$(it,mlev0-int:mlev0+int,{t0:t1},{n0:n1})

     uu0 = c->$u_varname$(it,mlev1-int:mlev1+int,{t0:t1},{n0:n1})
     vu0 = c->$v_varname$(it,mlev1-int:mlev1+int,{t0:t1},{n0:n1})

    else

     ul0 = c->$u_varname$(it,mlev0,{t0:t1},{n0:n1})
     vl0 = c->$v_varname$(it,mlev0,{t0:t1},{n0:n1})

     uu0 = c->$u_varname$(it,mlev1,{t0:t1},{n0:n1})
     vu0 = c->$v_varname$(it,mlev1,{t0:t1},{n0:n1})

    end if

;===================================
; Calculate storm-relative winds
;===================================

  ; Upper level
    ul = ul0 - u_cyc(cn,it-1)
    vl = vl0 - v_cyc(cn,it-1)

  ; Lower level
    uu = uu0 - u_cyc(cn,it-1)
    vu = vu0 - v_cyc(cn,it-1)

  ; Add metadata from original horizontal wind arrays
    copy_VarCoords(ul0, ul)
    copy_VarAtts(ul0, ul)

    copy_VarCoords(vl0, vl)
    copy_VarAtts(vl0, vl)

    copy_VarCoords(uu0, uu)
    copy_VarAtts(uu0, uu)

    copy_VarCoords(vu0, vu)
    copy_VarAtts(vu0, vu)

  ; Calculate relative vorticity using centered finite differences
  ; Option '2' --> boundary points estimated using one-sided difference scheme
    vort_l = uv2vr_cfd(ul, vl, latC, lonC, 2)
    copy_VarCoords(ul, vort_l)
    vort_l@units = "s~S~-1~N~"
    vort_l@name = "Relative vorticity on model levels"

    vort_u = uv2vr_cfd(uu, vu, latC, lonC, 2)
    copy_VarCoords(uu, vort_u)
    vort_u@units = "s~S~-1~N~"
    vort_u@name = "Relative vorticity on model levels"

  ; Coriolis parameter and absolute vorticity
    f1   = coriolis_param(latC)
    f0   = conform_dims(dimsizes(vort_l),f1,1)
    copy_VarMeta(vort_l,f0)
    copy_VarCoords(vort_l,f0)
    f0@description = "Coriolis parameter"
    f0@name        = "Coriolis parameter"

  ; Don't need two 'avo' variables, right?
    avo  = vort_l + f0
    copy_VarMeta(vort_l,avo)
    copy_VarAtts(vort_l,avo)
    copy_VarCoords(vort_l,avo)
    avo@description = "Absolute vorticity"
    avo@name        = "Absolute vorticity on model levels"
    avo@long_name   = "Absolute vorticity"

;===========================================================
; Calculate rate of strain (following Nguyen et al. 2011)
;===========================================================

    dim1   = "latitude"
    dim2   = "longitude"

  ; Calculate horizontal [lat,lon] derivatives of 'u' (du_dx, du_dy)
    dul    = grad_latlon_cfd(ul, ul&$dim1$, ul&$dim2$, False, False)
    dul_dy = dul[0]
    dul_dx = dul[1]
    delete(dul)

    duu    = grad_latlon_cfd(uu, uu&$dim1$, uu&$dim2$, False, False)
    duu_dy = duu[0]
    duu_dx = duu[1]
    delete(duu)

  ; Calculate horizontal [lat,lon] derivatives of 'v' (dv_dx, dv_dy)
    dvl     = grad_latlon_cfd(vl, vl&$dim1$, vl&$dim2$, False, False)
    dvl_dy  = dvl[0]
    dvl_dx  = dvl[1]
    delete(dvl)

    dvu     = grad_latlon_cfd(vu, vu&$dim1$, vu&$dim2$, False, False)
    dvu_dy  = dvu[0]
    dvu_dx  = dvu[1]
    delete(dvu)

  ; Calculate strain rate components --> sqrt[ (str)^2 + (shr)^2 ]
  ; Stretching (du_dx - dv_dy) and shearing (dv_dx + du_dy) deformation
    str1a_l = dul_dx - dvl_dy
    str2a_l = dvl_dx + dul_dy
    str1_l  = str1a_l ^ 2
    str2_l  = str2a_l ^ 2

    str1a_u = duu_dx - dvu_dy
    str2a_u = dvu_dx + duu_dy
    str1_u  = str1a_u ^ 2
    str2_u  = str2a_u ^ 2

  ; Finally, calculate strain rate 
    if (lay .eq. 1) then 
     str_l  = sqrt( str1_l(:,:,:) + str2_l(:,:,:) )
     str_u  = sqrt( str1_u(:,:,:) + str2_u(:,:,:) )
    else
     str_l  = sqrt( str1_l(:,:) + str2_l(:,:) )
     str_u  = sqrt( str1_u(:,:) + str2_u(:,:) )
    end if 

;======================================================
; Read in variables from 'pd' stream (model levels)
;======================================================

  ; Read in basic variables from 'pd' stream 
    lon1D  = d->longitude_1({n0:n1})                ; longitude [xxx grid points]
    lat1D  = d->latitude_1({t0:t1})                 ; latitude  [xxx grid points]
    lonD   = d->longitude({n0:n1})                  ; longitude [xxx grid points]
    latD   = d->latitude({t0:t1})                   ; latitude  [xxx grid points]
    hybD   = d->hybrid_ht(:)			    ; 63 model ('rho') levels

  ; Read in pressure and vertical velocity
    if (lay .eq. 1) then 

     p_l   = d->p(it-2,mlev0-int:mlev0+int,{t0:t1},{n0:n1})
     w_l   = d->dz_dt(it-2,mlev0-int:mlev0+int,{t0:t1},{n0:n1})
     th_l  = d->theta(it-2,mlev0-int:mlev0+int,{t0:t1},{n0:n1})
     t_l   = d->temp(it-2,mlev0-int:mlev0+int,{t0:t1},{n0:n1})

     p_u   = d->p(it-2,mlev1-int:mlev1+int,{t0:t1},{n0:n1})
     w_u   = d->dz_dt(it-2,mlev1-int:mlev1+int,{t0:t1},{n0:n1})
     th_u  = d->theta(it-2,mlev1-int:mlev1+int,{t0:t1},{n0:n1})
     t_u   = d->temp(it-2,mlev1-int:mlev1+int,{t0:t1},{n0:n1})

    else

     p_l   = d->p(it-2,mlev0,{t0:t1},{n0:n1})
     w_l   = d->dz_dt(it-2,mlev0,{t0:t1},{n0:n1})
     th_l  = d->theta(it-2,mlev0,{t0:t1},{n0:n1})
     t_l   = d->temp(it-2,mlev0,{t0:t1},{n0:n1})

     p_u   = d->p(it-2,mlev1,{t0:t1},{n0:n1})
     w_u   = d->dz_dt(it-2,mlev1,{t0:t1},{n0:n1})
     th_u  = d->theta(it-2,mlev1,{t0:t1},{n0:n1})
     t_u   = d->temp(it-2,mlev1,{t0:t1},{n0:n1})

    end if 

  ; Calculate dry air density from pressure and temperature
    rho_l  = p_l / (rd * t_l)
    rho_u  = p_u / (rd * t_u)

  ; Add metadata
    copy_VarCoords(p_l, rho_l)
    rho_l@description = "Dry air density"
    rho_l@units = "kg m~S~-3~N~"

    copy_VarCoords(p_u, rho_u)
    rho_u@description = "Dry air density"
    rho_u@units = "kg m~S~-3~N~"

  ; Calculate mass flux (lower level/layer)
    mflux_l = rho_l * w_l
    copy_VarCoords(w_l, mflux_l)
    copy_VarAtts(w_l, mflux_l)
    mflux_l@units = "kg m~S~-2~N~ s~S~-1~N~"
    mflux_l@standard_name = "vertical_mass_flux"
    mflux_l@long_name     = "Vertical mass flux"
    mflux_l@title         = "Vertical mass flux"
    mflux_l@name          = "vert_mflux"

  ; Calculate mass flux (upper level/layer)
    mflux_u = rho_u * w_u
    copy_VarCoords(w_u, mflux_u)
    copy_VarAtts(w_u, mflux_u)
    mflux_u@units = "kg m~S~-2~N~ s~S~-1~N~"
    mflux_u@standard_name = "vertical_mass_flux"
    mflux_u@long_name     = "Vertical mass flux"
    mflux_u@title         = "Vertical mass flux"
    mflux_u@name          = "vert_mflux"

  ; Grid spacing (lon, lat)
    dx0    = lonD(1) - lonD(0)
    dy0    = latD(1) - latD(0)

;==========================
; Loop over model levels 
;==========================

    hy0   = sprintf("%0.0f",hybD(mlev0))
    hy0_l = sprintf("%0.0f",hybD(mlev0-int))
    hy1_l = sprintf("%0.0f",hybD(mlev0+int))

    hy1   = sprintf("%0.0f",hybD(mlev1))
    hy0_u = sprintf("%0.0f",hybD(mlev1-int))
    hy1_u = sprintf("%0.0f",hybD(mlev1+int))

  ; Read in variables on model levels
    if (lay .eq. 1) then

     print("Working on lower layer between " + hy0_l+ " and "+hy1_l+" m AGL")
     print("Working on upper layer between " + hy0_u+ " and "+hy1_u+" m AGL")

   ; 'pc' stream (upper level/layer)
     uu_plane     = dim_avg_n_Wrap( uu(:,:,:),0)
     vu_plane     = dim_avg_n_Wrap( vu(:,:,:),0)
     vrt_u_plane  = dim_avg_n_Wrap( vort_u(:,:,:),0)

     str_u_plane  = dim_avg_n_Wrap( str_u(:,:,:),0)
     copy_VarCoords(uu(0,:,:),str_u_plane)

   ; 'pc' stream (lower level/layer)
     ul_plane     = dim_avg_n_Wrap( ul(:,:,:),0)
     vl_plane     = dim_avg_n_Wrap( vl(:,:,:),0)
     vrt_l_plane  = dim_avg_n_Wrap( vort_l(:,:,:),0)

     str_l_plane  = dim_avg_n_Wrap( str_l(:,:,:),0)
     copy_VarCoords(ul(0,:,:),str_l_plane)

   ; 'pd' stream (upper level/layer)
     wu_plane     = dim_avg_n_Wrap( w_u(:,:,:),0)
     prs_u_plane  = dim_avg_n_Wrap( p_u(:,:,:), 0)
     th_u_plane   = dim_avg_n_Wrap( th_u(:,:,:),0)

     flux_u_plane = dim_avg_n_Wrap( mflux_u(:,:,:),0)
     copy_VarCoords(uu(0,:,:),flux_u_plane)

   ; 'pd' stream (lower level/layer)
     wl_plane     = dim_avg_n_Wrap( w_l(:,:,:),0)
     prs_l_plane  = dim_avg_n_Wrap( p_l(:,:,:), 0)
     th_l_plane   = dim_avg_n_Wrap( th_l(:,:,:),0)

     flux_l_plane = dim_avg_n_Wrap( mflux_l(:,:,:),0)
     copy_VarCoords(ul(0,:,:),flux_l_plane)

    else

      print("Working on " +hy0+ " and " +hy1+ " m AGL")

    ; 'pc' stream (upper level)
      uu_plane     = uu(:,:)
      vu_plane     = vu(:,:)
      vrt_u_plane  = vort_u(:,:)

      str_u_plane  = str_u(:,:)
      copy_VarCoords(uu, str_u_plane)

    ; 'pc' stream (lower level)
      ul_plane     = ul(:,:)
      vl_plane     = vl(:,:)
      vrt_l_plane  = vort_l(:,:)

      str_l_plane  = str_l(:,:)
      copy_VarCoords(ul, str_l_plane)

    ; 'pd' stream (upper level)
      wu_plane     = w_u(:,:)
      prs_u_plane  = p_u(:,:)
      th_u_plane   = th_u(:,:)
      flux_u_plane = mflux_u(:,:)

    ; 'pd' stream (lower level)
      wl_plane     = w_l(:,:)
      prs_l_plane  = p_l(:,:)
      th_l_plane   = th_l(:,:)
      flux_l_plane = mflux_l(:,:)

    end if

  ; Calculate vector windspeed and smooth
    spd_u_plane    = sqrt( (uu_plane ^ 2) + (vu_plane ^ 2) )
    copy_VarCoords(uu_plane, spd_u_plane)
    vort_u_smth    = smth9_Wrap(vrt_u_plane, 0.5, 0.5, True)

    spd_l_plane    = sqrt( (ul_plane ^ 2) + (vl_plane ^ 2) )
    copy_VarCoords(ul_plane, spd_l_plane)
    vort_l_smth    = smth9_Wrap(vrt_l_plane, 0.5, 0.5, True)

  ; Combine all pc and pd stream variables into larger arrays
    pc_u_plane     = (/uu_plane, vu_plane, spd_u_plane, \
                       vrt_u_plane, vrt_u_plane, str_u_plane/)

    pc_l_plane     = (/ul_plane, vl_plane, spd_l_plane, \
                       vrt_l_plane, vrt_l_plane, str_l_plane/)

  ; Second 'prs_plane' is a placeholder
    pd_u_plane     = (/wu_plane, prs_u_plane, prs_u_plane, flux_u_plane/)
    pd_l_plane     = (/wl_plane, prs_l_plane, prs_l_plane, flux_l_plane/)

;=============================================================  
; Call external function to calculate storm centre position
;============================================================= 

    ; Set radial distance (º) when looking for vort min in ringlike phase   
      r0            = 0.10 

    ; Lower level/layer
      centre_arr_l  = st_centre(0, it, cn, r0, centre, vrt_l_plane, prs_l_plane, spd_l_plane, \
                                centre_new_l, dy0, dx0, "comp", ct, "pd", "off")

    ; Upper level/layer
      centre_arr_u  = st_centre(0, it, cn, r0, centre, vrt_u_plane, prs_u_plane, spd_u_plane, \
                                centre_new_u, dy0, dx0, "comp", ct, "pd", "off")

    ; Updated storm track information contained in 'centre_new'
      centre_new_l  = centre_arr_l[0]
      centre_new_u  = centre_arr_u[0]

;=====================================================================  
; Call external function to switch to cylindrical coordinate system  
;=====================================================================  

    ; Call 'setup_cyl' and output pressure/vorticity on cylindrical grid 
    ; Option before "comp" --> basic ('0') or MBGT ('1') variables output    
    ; Option before 'num_l' --> 'ilev=0' 
    ; Run 'setup_cyl' for upper level/layer
      cyl_arr_u  = setup_cyl(ct, dist, ll_arr, lonC, latC, lonD, latD, lon1D, lat1D, \
                             pc_u_plane, pd_u_plane, num_int, centre_new_u, \
                             sc, sm, ar, nr, 0, num_l, 0, "comp", ct, 1, \
                             xpos_u_all, ypos_u_all, xcen_u_all, ycen_u_all)

      xpos_u_all = cyl_arr_u[0]
      ypos_u_all = cyl_arr_u[1]
      xcen_u_all = cyl_arr_u[2]
      ycen_u_all = cyl_arr_u[3] 

      lat_u_max  = cyl_arr_u[4]
      lon_u_max  = cyl_arr_u[5]
      radii_u    = cyl_arr_u[6]
      rad_u_size = cyl_arr_u[7]
      thetas_u   = cyl_arr_u[19] 

      uu_int     = cyl_arr_u[8]
      vu_int     = cyl_arr_u[9]
      spd_u_int  = cyl_arr_u[10]
      vort_u_int = cyl_arr_u[11]
      vt_u_int   = cyl_arr_u[14]
      w_u_int    = cyl_arr_u[15]
      flux_u_int = cyl_arr_u[24]

    ; Also run 'setup_cyl' for lower level/layer
      cyl_arr_l  = setup_cyl(ct, dist, ll_arr, lonC, latC, lonD, latD, lon1D, lat1D, \
                             pc_l_plane, pd_l_plane, num_int, centre_new_l, \
                             sc, sm, ar, nr, 0, num_l, 0, "comp", ct, 1, \
                             xpos_l_all, ypos_l_all, xcen_l_all, ycen_l_all)

      xpos_l_all = cyl_arr_l[0]
      ypos_l_all = cyl_arr_l[1]
      xcen_l_all = cyl_arr_l[2]
      ycen_l_all = cyl_arr_l[3]

      lat_l_max  = cyl_arr_l[4]
      lon_l_max  = cyl_arr_l[5]
      radii_l    = cyl_arr_l[6]
      rad_l_size = cyl_arr_l[7]
      thetas_l   = cyl_arr_l[19]

      ul_int     = cyl_arr_l[8]
      vl_int     = cyl_arr_l[9]
      spd_l_int  = cyl_arr_l[10]
      vort_l_int = cyl_arr_l[11]
      vt_l_int   = cyl_arr_l[14]
      w_l_int    = cyl_arr_l[15]
      flux_l_int = cyl_arr_l[24]

;===============================================
; Tidy up and input values into larger arrays 
;===============================================

  ; Define new arrays to hold azimuthally-averaged variables
    flux_u_azi = new( (/dimsizes(radii)/), float)
    flux_l_azi = new( (/dimsizes(radii)/), float)

  ; Now calculate the azimuthal average of relative vorticity
    do irad = 0, dimsizes(radii)-1
     flux_u_azi(irad) = avg( flux_u_int(irad,:) )
     flux_l_azi(irad) = avg( flux_l_int(irad,:) )
    end do

  ; Fill larger arrays with calculated vertical mass flux 
    flux_u_int0(:,:,ct) = (/flux_u_int(:,:)/)
    flux_l_int0(:,:,ct) = (/flux_l_int(:,:)/)

  ; Tidy up
    delete([/ul_int, vl_int, spd_l_int, vort_l_int, vt_l_int, w_l_int, flux_l_int/])
    delete([/uu_int, vu_int, spd_u_int, vort_u_int, vt_u_int, w_u_int, flux_u_int/])
    delete([/wu_plane, vrt_u_plane, uu_plane, vu_plane, spd_u_plane, str_u_plane/])
    delete([/wl_plane, vrt_l_plane, ul_plane, vl_plane, spd_l_plane, str_l_plane/])
    delete([/prs_u_plane, th_u_plane, flux_u_plane, prs_l_plane, th_l_plane, flux_l_plane/])

    ct  = ct + 1		    ; Counter variable (time)

 ; Tidy up before next iteration (all variables from 'pc' and 'pd' streams)
   delete([/lonC,latC,lon1C,lat1C,ul0,vl0,uu0,vu0,ul,vl,uu,vu,vort_l,vort_u/])
   delete([/dul_dy, dul_dx, duu_dy, duu_dx, dvl_dy, dvl_dx, dvu_dy, dvu_dx, str_l, str_u/])
   delete([/lonD, latD, lon1D, lat1D, w_l, p_l, th_l, t_l, p_u, w_u, th_u, t_u/])
   delete([/rho_l, rho_u, mflux_l, mflux_u/])

  end do     ; End time loop (do it = 0, times-1)

 end do      ; End loop over input file times (do st = 0, dimsizes(dat)-1

;===============================================
; Produce scatter plots of vertical mass flux
;=============================================== 

 ; Calculate difference in mass flux between two levels/layers 
   flux_d_int0 = flux_u_int0(:,:,:) - flux_l_int0(:,:,:)

 ; Also calculate azimuthally-averaged fluxes
   flux_u_avg0 = dim_avg_n_Wrap(flux_u_int0(:,:,:), 1)
   flux_l_avg0 = dim_avg_n_Wrap(flux_l_int0(:,:,:), 1)
   flux_d_avg0 = dim_avg_n_Wrap(flux_d_int0(:,:,:), 1)

 ; Vertical mass flux plotting resources 
   mf_opts                      = True 
   mf_opts@gsnDraw              = False
   mf_opts@gsnFrame             = False
   mf_opts@gsnMaximize          = False 

 ; Y-axis
   if (diff .eq. 1) then 
    mf_opts@tiMainString        = "Vertical mass flux ("+hy1+" m - "+hy0+" m)"
   else
    mf_opts@tiMainString        = "Vertical mass flux at "+hy1+" m "
   end if 

   mf_opts@tiYAxisString        = "Vertical mass flux (kg m~S~-2~N~ s~S~-1~N~)"
   mf_opts@tiYAxisFontHeightF   = 0.0150

   if (azi .eq. 1) then 
    mf_opts@trYMaxF             = 1.5
    mf_opts@trYMinF             = -1.5
   else 
    mf_opts@trYMaxF             = 6.0
    mf_opts@trYMinF             = -4.0
   end if 

 ; X-axis (define our own labels)
   mf_opts@tiXAxisString        = "Radius (km)"
   mf_opts@tiXAxisFontHeightF   = 0.0150
   mf_opts@tmXBLabelFontHeightF = 0.0125
   mf_opts@tmXBMode             = "Explicit"
   mf_opts@tmXBValues           = (/0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45/)
   mf_opts@tmXBLabels           = (/"0", "5", "10", "15", "20", "25", "30", "35", "40", "45"/)

 ; EDIT -- how to reduce whitespace between 'columns'?
;   mf_opts@vpWidthF             = 0.5
;   mf_opts@vpHeightF            = 0.7

 ; Additional resources 
   mf_opts@xyMarkLineMode       = "Markers"
   mf_opts@xyMonoMarkerColor    = True
   mf_opts@xyMarkerColor        = "RoyalBlue1"
   mf_opts@xyMonoMarkerSize     = True 
   mf_opts@xyMonoMarkerThickness= True 
 
   if (azi .eq. 1) then 
    mf_opts@xyMarkerSizeF       = 0.03
    mf_opts@xyMarkerThicknessF  = 0.5
   else
    mf_opts@xyMarkerSizeF       = 0.002
    mf_opts@xyMarkerThicknessF  = 0.02
   end if 

 ; Turn off tickmarks
   mf_opts@tmXTOn               = "False"
   mf_opts@tmYROn               = "False"

 ; Additional resources (x-axis, etc)
   mf_opts@trXMaxF              = radii(dimsizes(radii)-2)
   mf_opts@trXMinF              = 0.05

 ; Open workstation and define output file path
   output = "$sam/nepartak/images/vc/mflux_scatter_"+typ+"_full"+full+"_sc"+sc+"_lay"+lay+"_"\
            +hy0+"_"+hy1+"_int"+int
   wks    = gsn_open_wks(opt, output) 

 ; Switch dimension order before plotting
   if (azi .eq. 1) then 
    flux_u_plot = flux_u_avg0(time|:, rad|:)
    flux_l_plot = flux_l_avg0(time|:, rad|:)
    flux_d_plot = flux_d_avg0(time|:, rad|:)
   else
    flux_u_plot = flux_u_int0(time|:, azi|:, rad|:)
    flux_l_plot = flux_l_int0(time|:, azi|:, rad|:)
    flux_d_plot = flux_d_int0(time|:, azi|:, rad|:)
   end if 

 ; Create object to hold all data for plotting 
   plot_u  = new( num_int, graphic)
   plot_l  = new( num_int, graphic)
   plot_d  = new( num_int, graphic)

   psize = dimsizes(plot_u) / 2

;==========================
; Finally, plot all data
;==========================

 if (diff .eq. 0) then 

 ; Coordinates of panels 
   xf = (/0.08, 0.58/)
   yf = (/0.70, 0.70/)
   wf = 0.40
   hf = 0.40

   mf_opts@vpXF      = xf(0)
   mf_opts@vpYF      = yf(0)
   mf_opts@vpWidthF  = wf
   mf_opts@vpHeightF = hf

 ; Plot size and position 
   mf2_opts                     = mf_opts
   mf2_opts@tiMainString        = ""
   mf2_opts@tiYAxisString       = ""
   mf2_opts@tiXAxisString       = ""
   mf2_opts@tmXBLabels          = (/"", "", "", "", "", "", "", "", "", ""/)

   mf2_opts@tmYLMode            = "Explicit"
   mf2_opts@tmYLValues          = (/-4.0, -2.0, 0.0, 2.0, 4.0, 6.0/)
   mf2_opts@tmYLLabels          = (/"", "", "", "", "", ""/)

 ; Upper level/layer
   plot_u(0) = gsn_csm_xy(wks, radii, flux_u_plot(0,:,:), mf_opts)
   do ct = 1, psize-1
    plot_u(ct) = gsn_csm_xy(wks, radii, flux_u_plot(ct,:,:), mf2_opts)
   end do 

   mfl_opts                     = mf_opts
   mfl_opts@tiMainString        = "Vertical mass flux at "+hy0+" m "
   mf_opts@xyMarkerColor        = "red4"
   mfl_opts@vpXF                = xf(1)
   mfl_opts@vpYF                = yf(1)
   mfl_opts@vpWidthF            = wf
   mfl_opts@vpHeightF           = hf

   mf3_opts                     = mfl_opts
   mf3_opts                     = mf2_opts
   mf3_opts@xyMarkerColor       = "red4"

 ; Lower level/layer 
   plot_l(0) = gsn_csm_xy(wks, radii, flux_l_plot(0,:,:), mfl_opts)
   do ct = 1, psize-1
    plot_l(ct) = gsn_csm_xy(wks, radii, flux_l_plot(ct,:,:), mf3_opts)    
   end do 

   draw(plot_u)
   draw(plot_l)
   frame(wks)

 else

 ; Difference between two levels/layers 
   if (azi .eq. 1) then 
    plot_d(0) = gsn_csm_xy(wks, radii, flux_d_plot(0,:), mf_opts)
   else
    plot_d(0) = gsn_csm_xy(wks, radii, flux_d_plot(0,:,:), mf_opts)
   end if 

 ; Some additional plotting resources 
   mf2_opts                     = mf_opts
   mf2_opts@tiMainString        = ""
   mf2_opts@tiYAxisString       = ""
   mf2_opts@tiXAxisString       = ""
   mf2_opts@tmXBLabels          = (/"", "", "", "", "", "", "", "", "", ""/)

   mf2_opts@tmYLMode            = "Explicit"
   mf2_opts@tmYLValues          = (/-4.0, -2.0, 0.0, 2.0, 4.0, 6.0/)
   mf2_opts@tmYLLabels          = (/"", "", "", "", "", ""/)

   if (azi .eq. 1) then 
    do ct = 1, psize-1
     plot_d(ct) = gsn_csm_xy(wks, radii, flux_d_plot(ct,:), mf2_opts)
    end do
   else
    do ct = 1, psize-1
     plot_d(ct) = gsn_csm_xy(wks, radii, flux_d_plot(ct,:,:), mf2_opts)
    end do
   end if 

   draw(plot_d)
   frame(wks)

 end if 

end 