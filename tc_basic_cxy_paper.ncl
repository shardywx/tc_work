; Script to calculate and save azimuthally averaged variables

; Run using:

; ncl dist=1.08 calc=1 nt=73 nr=21 ar=1.0 typ=\"stoa\" full=6 cn0=\"slp\" 
; sc=5 sm=0 tc_basic_cxz_paper.ncl

; 'dist'  = size of box following storm (degrees)
; 'calc'  = centre calculated on each level (0), calculated offline (1) or using 980 m vort (2)
; 'typ'	  = VC phase ("sym","asym","stoa","atos")
; 'full'  = choice of input data
; 'ar'    = distance in degrees of outer radius (1.0, 2.0, etc)
; 'nr'    = number of radial circles between r = 0 and r = 'ar' (21, 41, 61, etc)
; 'nt'    = number of azimuth angles in cylindrical grid
; 'cn0'   = offline storm centre: "slp", "vort", "geo", "geo_sm"

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

begin

;==============================
; Prelash (reading in files)
;==============================

; Dimensions of text files to read data in from (S-phase, S-A, A-phase, A-S)    
; 'vc_arr_all(:,0)' = Group 1   ; 'vc_arr_nw(:,0)' = Group 1a (no weak VCs)
; 'vc_arr_all(:,1)' = Group 2   ; 'vc_arr_nw(:,1)' = Group 2a (no weak VCs)
; 'vc_arr_all(:,2)' = Group 3
; 'vc_arr_all(:,3)' = Group 4
  vc_arr_all = (/ (/17,25,31,29/), (/16,14,30,24/), (/16,29,31,38/), (/9,15,31,24/)  /)
  vc_arr_nw  = (/ (/14,19,31/), (/12,12,30/), (/13,24,31/), (/9,12,31/)  /)

; String containing info on width of radial averaging band
  ar0  = ar * 100
  rstr = sprintf("%0.0fkm", ar0)

; For now, explicitly tell NCL which dates and simulations to analyse
  if (typ .eq. "sym") then

   if (full .eq. 3) then         ; Group 2 (12 UTC 5 -- 6 July)

    dat  = asciiread("$ar/text/vc_sym_dat_group2.txt",(/vc_arr_all(0,1)/),"string")
    ens0 = asciiread("$ar/text/vc_sym_sim_group2.txt",(/vc_arr_all(0,1)/),"string")
    ts0  = asciiread("$ar/text/vc_sym_ts_group2.txt",(/vc_arr_all(0,1)/),"integer")
    tf0  = ts0 + 1

   elseif (full .eq. 4) then         ; Group 1 (12 UTC 4 -- 5 July)

    dat  = asciiread("$ar/text/vc_sym_dat_group1.txt",(/vc_arr_all(0,0)/),"string")
    ens0 = asciiread("$ar/text/vc_sym_sim_group1.txt",(/vc_arr_all(0,0)/),"string")
    ts0 = asciiread("$ar/text/vc_sym_ts_group1.txt",(/vc_arr_all(0,0)/),"integer")
    tf0 = ts0 + 1

   elseif (full .eq. 6) then         ; Group 1 (no 'weak' cases)

    dat  = asciiread("$ar/text/vc_sym_dat_group1a.txt",(/vc_arr_nw(0,0)/),"string")
    ens0 = asciiread("$ar/text/vc_sym_sim_group1a.txt",(/vc_arr_nw(0,0)/),"string")
    ts0 = asciiread("$ar/text/vc_sym_ts_group1a.txt",(/vc_arr_nw(0,0)/),"integer")
    tf0 = ts0 + 1

   elseif (full .eq. 5) then         ; Group 2 (no 'weak' cases)

    dat  = asciiread("$ar/text/vc_sym_dat_group2a.txt",(/vc_arr_nw(0,1)/),"string")
    ens0 = asciiread("$ar/text/vc_sym_sim_group2a.txt",(/vc_arr_nw(0,1)/),"string")
    ts0 = asciiread("$ar/text/vc_sym_ts_group2a.txt",(/vc_arr_nw(0,1)/),"integer")
    tf0 = ts0 + 1

   elseif (full .eq. 7) then 

    dat  = asciiread("$ar/text/vc_sym_dat_group3.txt",(/vc_arr_all(0,2)/),"string")
    ens0 = asciiread("$ar/text/vc_sym_sim_group3.txt",(/vc_arr_all(0,2)/),"string")
    ts0  = asciiread("$ar/text/vc_sym_ts_group3.txt",(/vc_arr_all(0,2)/),"integer")
    tf0  = ts0 + 1

   elseif (full .eq. 8) then

    dat  = asciiread("$ar/text/vc_sym_dat_group4.txt",(/vc_arr_all(0,3)/),"string")
    ens0 = asciiread("$ar/text/vc_sym_sim_group4.txt",(/vc_arr_all(0,3)/),"string")
    ts0  = asciiread("$ar/text/vc_sym_ts_group4.txt",(/vc_arr_all(0,3)/),"integer")
    tf0  = ts0 + 1

   end if

  elseif (typ .eq. "stoa") then

   if (full .eq. 3) then         ; Group 2 (12 UTC 5 -- 6 July)

    dat  = asciiread("$ar/text/vc_stoa_dat_group2.txt",(/vc_arr_all(1,1)/),"string")
    ens0 = asciiread("$ar/text/vc_stoa_sim_group2.txt",(/vc_arr_all(1,1)/),"string")
    ts0 = asciiread("$ar/text/vc_stoa_ts_group2.txt",(/vc_arr_all(1,1)/),"integer")
    tf0 = ts0 + 1

  elseif (full.eq. 4) then       ; Group 1 (12 UTC 4 -- 5 July)                                  

    dat  = asciiread("$ar/text/vc_stoa_dat_group1.txt",(/vc_arr_all(1,0)/),"string")
    ens0 = asciiread("$ar/text/vc_stoa_sim_group1.txt",(/vc_arr_all(1,0)/),"string")
    ts0 = asciiread("$ar/text/vc_stoa_ts_group1.txt",(/vc_arr_all(1,0)/),"integer")
    tf0 = ts0 + 1

   elseif (full .eq. 6) then         ; Group 1 (no 'weak' cases)                       

    dat  = asciiread("$ar/text/vc_stoa_dat_group1a.txt",(/vc_arr_nw(1,0)/),"string")
    ens0 = asciiread("$ar/text/vc_stoa_sim_group1a.txt",(/vc_arr_nw(1,0)/),"string")
    ts0 = asciiread("$ar/text/vc_stoa_ts_group1a.txt",(/vc_arr_nw(1,0)/),"integer")
    tf0 = ts0 + 1

   elseif (full .eq. 5) then         ; Group 2 (no 'weak' cases)                       

    dat  = asciiread("$ar/text/vc_stoa_dat_group2a.txt",(/vc_arr_nw(1,1)/),"string")
    ens0 = asciiread("$ar/text/vc_stoa_sim_group2a.txt",(/vc_arr_nw(1,1)/),"string")
    ts0 = asciiread("$ar/text/vc_stoa_ts_group2a.txt",(/vc_arr_nw(1,1)/),"integer")
    tf0 = ts0 + 1

   elseif (full .eq. 7) then 
  
    dat  = asciiread("$ar/text/vc_s2a_dat_group3.txt",(/vc_arr_all(1,2)/),"string")
    ens0 = asciiread("$ar/text/vc_s2a_sim_group3.txt",(/vc_arr_all(1,2)/),"string")
    ts0  = asciiread("$ar/text/vc_s2a_ts_group3.txt",(/vc_arr_all(1,2)/),"integer")
    tf0  = ts0 + 1

   elseif (full .eq. 8) then

    dat  = asciiread("$ar/text/vc_s2a_dat_group4.txt",(/vc_arr_all(1,3)/),"string")
    ens0 = asciiread("$ar/text/vc_s2a_sim_group4.txt",(/vc_arr_all(1,3)/),"string")
    ts0  = asciiread("$ar/text/vc_s2a_ts_group4.txt",(/vc_arr_all(1,3)/),"integer")
    tf0  = ts0 + 1

   end if 

  elseif (typ .eq. "asym") then

   if (full .eq. 3) then         ; Group 2 (12 UTC 5 -- 6 July)

    dat  = asciiread("$ar/text/vc_asym_dat_group2.txt",(/vc_arr_all(2,1)/),"string")
    ens0 = asciiread("$ar/text/vc_asym_sim_group2.txt",(/vc_arr_all(2,1)/),"string")
    ts0 = asciiread("$ar/text/vc_asym_ts_group2.txt",(/vc_arr_all(2,1)/),"integer")
    tf0 = ts0 + 1

  elseif (full.eq. 4) then       ; Group 1 (12 UTC 4 -- 5 July)                                  

    dat  = asciiread("$ar/text/vc_asym_dat_group1.txt",(/vc_arr_all(2,0)/),"string")
    ens0 = asciiread("$ar/text/vc_asym_sim_group1.txt",(/vc_arr_all(2,0)/),"string")
    ts0 = asciiread("$ar/text/vc_asym_ts_group1.txt",(/vc_arr_all(2,0)/),"integer")
    tf0 = ts0 + 1

   elseif (full .eq. 6) then         ; Group 1 (no 'weak' cases)                       

    dat  = asciiread("$ar/text/vc_asym_dat_group1a.txt",(/vc_arr_nw(2,0)/),"string")
    ens0 = asciiread("$ar/text/vc_asym_sim_group1a.txt",(/vc_arr_nw(2,0)/),"string")
    ts0 = asciiread("$ar/text/vc_asym_ts_group1a.txt",(/vc_arr_nw(2,0)/),"integer")
    tf0 = ts0 + 1

   elseif (full .eq. 5) then         ; Group 2 (no 'weak' cases)                       

    dat  = asciiread("$ar/text/vc_asym_dat_group2a.txt",(/vc_arr_nw(2,1)/),"string")
    ens0 = asciiread("$ar/text/vc_asym_sim_group2a.txt",(/vc_arr_nw(2,1)/),"string")
    ts0 = asciiread("$ar/text/vc_asym_ts_group2a.txt",(/vc_arr_nw(2,1)/),"integer")
    tf0 = ts0 + 1

   elseif (full .eq. 7) then 

    dat  = asciiread("$ar/text/vc_asym_dat_group3.txt",(/vc_arr_all(2,2)/),"string")
    ens0 = asciiread("$ar/text/vc_asym_sim_group3.txt",(/vc_arr_all(2,2)/),"string")
    ts0  = asciiread("$ar/text/vc_asym_ts_group3.txt",(/vc_arr_all(2,2)/),"integer")
    tf0  = ts0 + 1

   elseif (full .eq. 8) then

    dat  = asciiread("$ar/text/vc_asym_dat_group4.txt",(/vc_arr_all(2,3)/),"string")
    ens0 = asciiread("$ar/text/vc_asym_sim_group4.txt",(/vc_arr_all(2,3)/),"string")
    ts0  = asciiread("$ar/text/vc_asym_ts_group4.txt",(/vc_arr_all(2,3)/),"integer")
    tf0  = ts0 + 1

   end if 

  elseif (typ .eq. "atos") then

   if (full .eq. 3) then         ; Group 2 (12 UTC 5 -- 6 July)

    dat  = asciiread("$ar/text/vc_atos_dat_group2.txt",(/vc_arr_all(3,1)/),"string")
    ens0 = asciiread("$ar/text/vc_atos_sim_group2.txt",(/vc_arr_all(3,1)/),"string")
    ts0 = asciiread("$ar/text/vc_atos_ts_group2.txt",(/vc_arr_all(3,1)/),"integer")
    tf0 = ts0 + 1

  elseif (full.eq. 4) then       ; Group 1 (12 UTC 4 -- 5 July)                                  

    dat  = asciiread("$ar/text/vc_atos_dat_group1.txt",(/vc_arr_all(3,0)/),"string")
    ens0 = asciiread("$ar/text/vc_atos_sim_group1.txt",(/vc_arr_all(3,0)/),"string")
    ts0 = asciiread("$ar/text/vc_atos_ts_group1.txt",(/vc_arr_all(3,0)/),"integer")
    tf0 = ts0 + 1

   elseif (full .eq. 6) then         ; Group 1 (no 'weak' cases)                       

    dat  = asciiread("$ar/text/vc_atos_dat_group1a.txt",(/vc_arr_nw(3,0)/),"string")
    ens0 = asciiread("$ar/text/vc_atos_sim_group1a.txt",(/vc_arr_nw(3,0)/),"string")
    ts0 = asciiread("$ar/text/vc_atos_ts_group1a.txt",(/vc_arr_nw(3,0)/),"integer")
    tf0 = ts0 + 1

   elseif (full .eq. 5) then         ; Group 2 (no 'weak' cases)                       

    dat  = asciiread("$ar/text/vc_atos_dat_group2a.txt",(/vc_arr_nw(3,1)/),"string")
    ens0 = asciiread("$ar/text/vc_atos_sim_group2a.txt",(/vc_arr_nw(3,1)/),"string")
    ts0 = asciiread("$ar/text/vc_atos_ts_group2a.txt",(/vc_arr_nw(3,1)/),"integer")
    tf0 = ts0 + 1

   elseif (full .eq. 7) then 

    dat  = asciiread("$ar/text/vc_a2s_dat_group3.txt",(/vc_arr_all(3,2)/),"string")
    ens0 = asciiread("$ar/text/vc_a2s_sim_group3.txt",(/vc_arr_all(3,2)/),"string")
    ts0  = asciiread("$ar/text/vc_a2s_ts_group3.txt",(/vc_arr_all(3,2)/),"integer")
    tf0  = ts0 + 1

   elseif (full .eq. 8) then

    dat  = asciiread("$ar/text/vc_a2s_dat_group4.txt",(/vc_arr_all(3,3)/),"string")
    ens0 = asciiread("$ar/text/vc_a2s_sim_group4.txt",(/vc_arr_all(3,3)/),"string")
    ts0  = asciiread("$ar/text/vc_a2s_ts_group4.txt",(/vc_arr_all(3,3)/),"integer")
    tf0  = ts0 + 1
 
   end if

  end if

;=====================================================================================
; Before starting the multiple loops, define the arrays to hold the averaged values
;=====================================================================================

; Define array dimension
  numTIMES = (dimsizes(dat) * 2 )
  print("There are "+numTIMES+" times in this composite.")

; 59 model levels in UM files 
  num_l = 59

; Mean tangential wind {v}
  mean_vtan0       = new((/num_l,numTIMES/),float)
  mean_vtan0!0     = "lev"
  mean_vtan0!1     = "time"
  mean_vtan0@description = "Azimuthally averaged tangential wind"
  mean_vtan0@units = "m s~S~-1~N~"

; Mean radial wind {u}
  mean_vrad0       = new((/num_l,numTIMES/),float)
  mean_vrad0!0     = "lev"
  mean_vrad0!1     = "time"
  mean_vrad0@description = "Azimuthally averaged radial wind"
  mean_vrad0@units = "m s~S~-1~N~"

; Mean vertical velocity {w}
  mean_vvel0       = new((/num_l,numTIMES/),float)
  mean_vvel0!0     = "lev"
  mean_vvel0!1     = "time"
  mean_vvel0@description = "Azimuthally averaged vertical velocity"
  mean_vvel0@units = "m s~S~-1~N~"

; Mean relative vorticity {zeta}
  mean_vort0       = new((/num_l,numTIMES/),float)
  mean_vort0!0     = "lev"
  mean_vort0!1     = "time"
  mean_vort0@description = "Azimuthally averaged relative vorticity"
  mean_vort0@units = "s~S~-1~N~"  

; Mean absolute vorticity {zeta + f0}
  mean_avo0        = new((/num_l,numTIMES/),float)
  mean_avo0!0      = "lev"
  mean_avo0!1      = "time"
  mean_avo0@description = "Azimuthally averaged absolute vorticity"
  mean_avo0@units  = "s~S~-1~N~"

; Print out relevant information about storm tracking 
  sc_arr      = (/"relative vorticity",\                          ; sc = 0  
                  "pressure",\                                    ; sc = 1
                  "smoothed pressure",\                           ; sc = 2 
                  "extra smoothed pressure",\                     ; sc = 3  
                  "vort_min during symmetric phase",\             ; sc = 4  
                  "windspeed minimum",\                           ; sc = 5 
                  "geopotential"/)                                ; sc = 6  

  print("Using "+sc_arr(sc)+" to calculate storm centre position...")

;==========================================================
; Now start loop over these chosen input files and times
;==========================================================

 ct = 0      ; Initialise counter variable (time)

 do st = 0, dimsizes(dat)-1 

  dt0 = tf0(st) - ts0(st)

; Find 'pc' and 'pd' files (model height level data)
  diri        = "/nfs/a319/earshar/"+dat(st)+"/"+ens0(st)
  fili_prefix = "201607"+dat(st)+"00Z_NPTK_4p4_L80_ra1t_"+ens0(st)
  fili_c1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pc.nc")
  fili_c      = diri+"/"+fili_c1
  fili_d1     =	systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pd.nc")
  fili_d      = diri+"/"+fili_d1

; Calculate number of times in each file using built-in function ('nc_times')
  numINFO_c   = nc_times(fili_c)
  numINPUT_c  = numINFO_c[0]		; Number of times in 'pc' stream
  time_c      = numINFO_c[1]		; Array of times from 'pc' stream file

  numINFO_d   = nc_times(fili_d)
  numINPUT_d  = numINFO_d[0]		; Number of times in 'pd' stream
  time_d      = numINFO_d[1]		; Array of times from 'pd' stream file 

;============================================================
; Calculate storm motion using built-in function ('st_rm')
;============================================================

; 'dat'       = initialisation time (02T12, 03T00, ...)
; 'ens0'      = ensemble simulation (em00, em01, ...)
; 'diri'      = path to input files (see above)
; 'dist'      = size of box following storm (degrees)
; 'mins'      = analyse 1-h (0) or 5-min (1) data

  storm_rel   = st_rm(dat(st), ens0(st), "$ar/text/", dist, 0)
  u_cyc       = storm_rel[0]		   ; Zonal wind 
  v_cyc       = storm_rel[1]		   ; Meridional wind
  vel_cyc     = storm_rel[2]		   ; Vector wind 
  numTIMES    = storm_rel[3]               ; Number of times in file
  lat_arr     = storm_rel[4]               ; Storm-following domain subset (lat)
  lon_arr     = storm_rel[5]               ; Storm-following domain subset (lon)
  centre      = storm_rel[6]		   ; Storm centre at each hour (lat/lon)

  print_clock("Working on ensemble member: "+ens0(st))

  llbox       = toint(dist*50)            ; Size of lat/lon arrays (see line below)
  dsize       = (/llbox,llbox/)           ; Size of lat/lon arrays (depends on 'dist')
  d0          = dsize(0)
  d1          = dsize(1)

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 1000000000
  end setvalues

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
  do dt = 0, numINPUT_c-1
   if (minute(dt).gt.30) then
     hour(dt) = hour(dt)+1
   end if
  end do
  delete(dt)

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

; Arrays to hold all 'xpos' and 'ypos' position arrays
  xpos_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  ypos_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)

; Arrays to hold lat/lon at centre of cylindrical grid 
  xcen_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  ycen_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)

; Choose which storm track to read in
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

  ; Create new array to hold centre positions
    centre_new = new( (/7, num_l, numINPUT_d, 2/), float)

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
    if (lay .eq. 1) then
     u0  = c->$u_varname$(it,13:25,{t0:t1},{n0:n1})    ; Zonal wind
     v0  = c->$v_varname$(it,13:25,{t0:t1},{n0:n1})    ; Meridional wind
    else
     u0  = c->$u_varname$(it,ilev,{t0:t1},{n0:n1})     ; Zonal wind
     v0  = c->$v_varname$(it,ilev,{t0:t1},{n0:n1})     ; Meridional wind
    end if

  ; Calculate storm-relative winds
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
    vort@name  = "Relative vorticity on model levels"

  ; Convert to standard units for plotting 
    vort       = vort * (10 ^ 4)
    vort@units = "10~S~-4~N~ s~S~-1~N~"

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

  ; Read in temperature, pressure, theta, vertical velocity and geopotential height
    if (lay .eq. 1) then
     p     = d->p(it-2,13:25,{t0:t1},{n0:n1})      ; Pressure (Pa) 
     w     = d->dz_dt(it-2,13:25,{t0:t1},{n0:n1})  ; Vertical velocity (m s-1)
     z     = d->ht(it-2,13:25,{t0:t1},{n0:n1})     ; Geopotential height (m)
     th    = d->theta(it-2,13:25,{t0:t1},{n0:n1})  ; Potential temperature (K)
    else
     p     = d->p(it-2,ilev,{t0:t1},{n0:n1})       ; Pressure (Pa)
     w     = d->dz_dt(it-2,ilev,{t0:t1},{n0:n1})   ; Vertical velocity (m s-1)
     z     = d->ht(it-2,ilev,{t0:t1},{n0:n1})      ; Geopotential height (m)
     th    = d->theta(it-2,ilev,{t0:t1},{n0:n1})   ; Potential temperature (K)
    end if

  ; Convert units of pressure and geopotential height before continuing
    p  = p / 100
    z  = z / 10

  ; Add metadata
    p@units = "hPa"
    z@units = "dam"

    dx0      = lonD(1) - lonD(0) ; Grid spacing (longitude)
    dy0      = latD(1) - latD(0) ; Grid spacing (latitude)

;==========================
; Loop over model levels 
;==========================

    hy0 = sprintf("%0.0f",hybD(13))
    hy1 = sprintf("%0.0f",hybD(25))

    do ilev = 0, num_l-1 ; Loop over all pressure levels

     hy = sprintf("%0.0f",hybD(ilev))
     print_clock("Working on model level " + ilev+ " ("+hy+" m AGL)")

     if (lay .eq. 1) then 

       print_clock("Working on layer between " + hy0+ " and "+hy1+" m AGL")

     ; pc stream
       u_plane      = dim_avg_n_Wrap( u(:,:,:), 0)
       v_plane      = dim_avg_n_Wrap( v(:,:,:), 0)
       vort_plane   = dim_avg_n_Wrap( vort(:,:,:), 0)
       avo_plane    = dim_avg_n_Wrap( avo(:,:,:), 0)

     ; pd stream
       w_plane      = dim_avg_n_Wrap( w(:,:,:), 0)
       prs_plane    = dim_avg_n_Wrap( p(:,:,:), 0)
       th_plane     = dim_avg_n_Wrap( th(:,:,:),0)

     else 

       hy2 = sprintf("%0.0f",hybD(ilev))
       print_clock("Working on " + hy2+ " m AGL")

     ; pc stream
       u_plane      = u(:,:)
       v_plane      = v(:,:)
       vort_plane   = vort(:,:)
       avo_plane    = avo(:,:)    

     ; pd stream
       w_plane      = w(:,:)
       prs_plane    = p(:,:)
       th_plane     = th(:,:)

     end if 

   ; Calculate vector windspeed and smooth  
     spd_plane              = sqrt( (u_plane ^ 2) + (v_plane ^ 2) )
     copy_VarCoords(u_plane, spd_plane)
     vort_pl_smth           = smth9_Wrap(vort_plane, 0.5, 0.5, True)

   ; Combine all 'pc' and 'pd' stream variables into larger arrays
     pc_plane               = (/u_plane, v_plane, spd_plane, \
                                vort_plane, avo_plane/)
     pd_plane               = (/w_plane, prs_plane/)

;=======================================================  
; Calculate storm centre position (external function)
;=======================================================  

   ; Set radial distance (º) when looking for vort min in ringlike phase 
     r0            = 0.08

     print_clock("Calculating storm centre position...")
     centre_arr    = st_centre(ilev, it, cn, r0, centre, vort_plane, prs_plane, spd_plane, \
                               centre_new, dy0, dx0)

   ; Updated storm track information contained in 'centre_new' 
   ; centre_new(0,0,it-2,:) --> vort_max  
   ; centre_new(1,0,it-2,:) --> prs_min
   ; centre_new(2,0,it-2,:) --> prs_min_sm 
   ; centre_new(3,0,it-2,:) --> prs_min_sm2   
   ; centre_new(4,0,it-2,:) --> vort_min   
   ; centre_new(5,0,it-2,:) --> spd_min
   ; centre_new(6,0,it-2,:) --> offline 
     centre_new     = centre_arr[0]

;=====================================================================
; Call external function to switch to cylindrical coordinate system
;=====================================================================

   ; Set 'ilev' option before function call below 
     if (lay .eq. 1) then
      ilev = 0
     end if

   ; Last option = 0 --> basic variables output 
   ; Last option = 1 --> additional MBGT variables output 
     print_clock("Translating to cylindrical grid...")
     cyl_arr  = setup_cyl(it, dist, ll_arr, lonC, latC, lonD, latD, lon1D, lat1D, \
                          pc_plane, pd_plane, numTIMES, centre_new, \
                          sc, sm, ar, nr, ilev, num_l, 0)

;     xpos_all = cyl_arr[0]
;     ypos_all = cyl_arr[1]
;     xcen_all = cyl_arr[2]
;     ycen_all = cyl_arr[3]

     lat_max  = cyl_arr[4]
     lon_max  = cyl_arr[5]
     radii    = cyl_arr[6]
     rad_size = cyl_arr[7]

     xpos     = cyl_arr[17]
     ypos     = cyl_arr[18]
     thetas   = cyl_arr[19]

     u_int    = cyl_arr[8]
     v_int    = cyl_arr[9]
     spd_int  = cyl_arr[10]
     vort_int = cyl_arr[11]
     avo_int  = cyl_arr[12]
     vr_int   = cyl_arr[13]
     vt_int   = cyl_arr[14]
     w_int    = cyl_arr[15]
     prs_int  = cyl_arr[16]

   ; Calculate mean (azimuthally-averaged) variables 
     do irad = 0, dimsizes(radii)-1
      mean_vtan = avg(vt_int(irad,:) )
      mean_vrad = avg(vr_int(irad,:) )
      mean_vvel = avg(w_int(irad,:) )
      mean_vort = avg(vort_int(irad,:) )
      mean_avo  = avg(avo_int(irad,:) )
     end do

   ; Now average radially and input into 2D arrays
     mean_vtan0(ilev,ct)    = avg(mean_vtan)
     mean_vrad0(ilev,ct)    = avg(mean_vrad)
     mean_vvel0(ilev,ct)    = avg(mean_vvel)
     mean_vort0(ilev,ct)    = avg(mean_vort)
     mean_avo0(ilev,ct)     = avg(mean_avo)

   ; Tidy up
     delete([/u_int,v_int,spd_int,w_int,vt_int,vr_int,vort_int,avo_int/])
     delete([/u_plane,v_plane,spd_plane,w_plane,vort_plane,avo_plane/])

     end do			    ; End model level loop (do ilev = 0, num_l-1)

    ct  = ct + 1		    ; Counter variable (time)

  ; Tidy up before next iteration (all variables from 'pc' and 'pd' streams)
    delete([/lonC,latC,lon1C,lat1C,u,v,vort,avo,f0,f1/])
    delete([/lonD,latD,lon1D,lat1D,str_xz,str_yz,kmh,t,p,w,z/])

   end do                   ; End time loop (do it = 0, times-1)

 ; Tidy up before next iteration
   delete([/utc_date,year,month,day,hour,minute,second,date_str,time_str,time_arr,title_arr/])
   delete([/numINFO_c,numINFO_d,numINPUT_c,time_c,numINPUT_d,time_d/])
   delete([/storm_rel,u_cyc,v_cyc,vel_cyc,lat_arr,lon_arr,centre/])

  end do                    ; End loop over times (do it = ts0, tf0

;=============================================================
; Average over vertical levels before writing to text files
;=============================================================

; 980 -- 3380 m layer-averaged relative vorticity (lower tropospheric)
  mean_vortf = dim_avg_n_Wrap(mean_vort0(13:25,:), 0)

; 3920 -- 15360 m layer-averaged vertical velocity (mid to upper-tropospheric)
  mean_vvelf = dim_avg_n_Wrap(mean_vvel0(27:54,:), 0)

;=======================
; Write to text files 
;=======================

; Add string to output file to denote chosen group of input files 
  if (full .eq. 3) then
   dstr = "group2"
  elseif (full .eq. 4) then
   dstr = "group1"
  elseif (full .eq. 5) then
   dstr = "group2a"
  elseif (full .eq. 6) then
   dstr = "group1a"
  elseif (full .eq. 7) then 
   dstr = "group3"
  elseif (full .eq. 8) then 
   dstr = "group4"
  end if

; Define output file paths
  diri_out = "$ar/text/comp_"+typ+"_"+dstr+"_"+rstr+"_sc"+sc+"_"
  vort_out = diri_out+"mean_vort.txt"
  vvel_out = diri_out+"mean_vvel.txt"

; Write out to text file using 'asciiwrite'
  asciiwrite(vort_out, mean_vortf)
  asciiwrite(vvel_out, mean_vvelf)

end 