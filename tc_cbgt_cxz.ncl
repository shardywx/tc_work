; Script to calculate circulation budget following a tropical storm

; Run using: 

; ncl fct=5 op=1 opt=\"pdf\" rad=0.76 dist=2.2 nbox=1 plt=0 \
; cn0=\"geo_sm\" tr=0 int=0 sc=5 typ=\"sym\" full=9 wrt=0 tc_cbgt_cxz.ncl

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'opt'   = output file format ("pdf" or "x11")
; 'ens0'  = ensemble member (em00 -- em11)
; 'dist'  = size of box following storm (degrees)
; 'ts0'   = starting time index (0-119)
; 'tf0'   = ending time index (0-119)
; 'nbox'  = number of ensemble box positions for circ. budget calculations (49,25,9,1)
; 'rad'   = radius of box around TC for budget calculations (integer multiple of dx,dy)
; 'fct'	  = 
; 'op'	  = 
; 'cn0'	  = offline storm centre: "slp", "vort", "geo", "geo_sm"
; 'sc'    = centre-finding method (replacing 'calc' in old scripts)

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"

load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/nc_levs.ncl"
load "$sam/ncl_func/st_centre.ncl"
load "$sam/ncl_func/setup_cyl.ncl"
load "$sam/ncl_func/ring_mono.ncl"

begin

;==============================
; Prelash (reading in files)
;==============================

; Dimensions of text files to read data in from (SYM, S2A, ASYM, A2S)
  vc_arr_all = (/ (/17,25,31,29,24/), (/16,14,30,24,24/), (/16,29,31,38,24/), (/9,15,31,24,24/) /)
  vc_arr_nw  = (/ (/14,19,31/), (/12,12,30/), (/13,24,31/), (/9,12,31/)  /)

  if (typ .eq. "sym") then 

; Newest methods only 
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

;=========================================================
; Now start loop over these chosen input files and times 
;=========================================================

; Counter variable (time)
  ct = 0

; Find 'pb' and 'pd' stream data (pressure and model height level data)
  diri        = "/nfs/a319/earshar/"+dat(0)+"/"+ens0(0)
  fili_prefix = "201607"+dat(0)+"00Z_NPTK_4p4_L80_ra1t_"+ens0(0)

  fili_b1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pb*.nc")
  fili_b      = diri+"/"+fili_b1

; 'pb' stream [number of times; array of times]
  numINFO_b   = nc_times(fili_b)
  numINPUT_b  = numINFO_b[0]
  time_b      = numINFO_b[1]

; 'pb' stream [number of levels; array of levels]
  nLEVS_b     = nc_levs(fili_b)
  num_l       = nLEVS_b[1]
  levs_b      = nLEVS_b[0]

; CHECK -- define diagnostic arrays outside of all main loops? 
; Create arrays to hold circulation tendency values 
; [box ensemble (nbox) * time (numINPUT_b) * levs (num_l)]                                          
  circ_inst   = new((/nbox,numINPUT_b,num_l/),float)
  hadv_ten    = new((/nbox,numINPUT_b,num_l/),float)
  str_ten     = new((/nbox,numINPUT_b,num_l/),float)
  flux_ten    = new((/nbox,numINPUT_b,num_l/),float)
  tilt_ten    = new((/nbox,numINPUT_b,num_l/),float)
  fric_ten    = new((/nbox,numINPUT_b,num_l/),float)

; Create array to hold storm track positions using different methods
  centre_new  = new((/7, num_l, numINPUT_b, 2/),"float")

;=============================
; Now loop over simulations
;=============================

 do st = 0, dimsizes(dat)-1

  dt0 = tf0(st) - ts0(st)

  print("ct = "+ct)

; Create array to hold time-averaged values for composites
  if (wrt .eq. 1) then
   print("Creating composite chart array...")

   if (st .eq. 0) then

 ; Dimensions of 4D array (levels; time averages; diagnostics; radial circles) 
   nlev        = 59
   ntimes      = dimsizes(dat)
   ndiag       = 11
   nrad        = nr

 ; Create array to hold time averages 
   plot_arr    = new( (/nlev, ntimes, ndiag, nrad/), "float")

   end if

  end if

; Find 'pb' and 'pd' stream data (pressure and model height level data)
  diri        = "/nfs/a319/earshar/"+dat(st)+"/"+ens0(st)
  fili_prefix = "201607"+dat(st)+"00Z_NPTK_4p4_L80_ra1t_"+ens0(st)

  fili_b1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pb*.nc")
  fili_b      = diri+"/"+fili_b1

  fili_d1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pd*.nc")
  fili_d      = diri+"/"+fili_d1

; 'pb' stream [number of times; array of times]
  numINFO_b   = nc_times(fili_b)
  numINPUT_b  = numINFO_b[0]
  time_b      = numINFO_b[1]

; 'pd' stream [number of times; array of times]
  numINFO_d   = nc_times(fili_d)
  numINPUT_d  = numINFO_d[0]
  time_d      = numINFO_d[1]

; 'pb' stream [number of levels; array of levels]
  nLEVS_b     = nc_levs(fili_b)
  numLEVS_b   = nLEVS_b[1]
  levs_b      = nLEVS_b[0]

;==================================================
; Calculate storm motion using built-in function
;==================================================

; 'dat'       = initialisation time (02T12, 03T00, ...)
; 'ens0'      = ensemble simulation (em00, em01, ...)
; 'diri'      = path to input files (see above)
; 'dist'      = size of box following storm (degrees)
; 'mins'      = analyse 1-h (0) or 5-min (1) data

; Zonal; meridional; vector wind; lat; lon; centre
  storm_rel   = st_rm(dat(st), ens0(st), "$ar/text/", dist, 0)
  u_cyc       = storm_rel[0]
  v_cyc       = storm_rel[1]
  vel_cyc     = storm_rel[2]
  lat_arr     = storm_rel[4]
  lon_arr     = storm_rel[5]
  centre      = storm_rel[6]

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 10000000000
  end setvalues

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
  do ct0 = 0, numINPUT_b-1
   if (minute(ct0).gt.30) then
     hour(ct0) = hour(ct0)+1
   end if
  end do
  delete(ct0)

; Create arrays to hold strings for output plots
  date_str  = new(numINPUT_b,string)
  time_str  = new(numINPUT_b,string)
  time_arr  = new(numINPUT_b,string)
  title_arr = new(numINPUT_b,string)

; 117 values in these arrays, vs 121 for 'time_c'
  date_str  = sprinti("%0.2iUTC ", hour) + \
              sprinti("%0.2i ", day) \
            + month_abbr(month)
  time_str  = sprinti("%0.2i UTC ", hour) + \
              sprinti("%0.2i ", day) \
            + month_abbr(month)
  time_arr  = sprinti("%0.2i", day) + \
              month_abbr(month) + \
              "_" + sprinti("%0.2iZ", hour)

; NOTE -- this shouldn't matter than much, because we're not looking at a single case 

;==================================================================================
;;;;;;;;;;;;;;;;;;;;;;;;;;; Start of main part of script ;;;;;;;;;;;;;;;;;;;;;;;;;;
;==================================================================================

  print_clock("Starting main part of script...")

; Define constants
  r_earth     = 6.37e6  ; Earth's radius (constant)
  rd          = 287.0   ; Gas constant for dry air (J/kg/K^2)
  cp0         = 1004.0  ; Specific heat of dry air at constant pressure (J/kg/K)
  kp          = 0.286   ; For Exner function calculations (0.286)
  pref        = 1000.0  ; Reference pressure p0 (1000 hPa)

; Retrieve size of arrays for plotting 
  llbox       = toint(dist*50)
  dsize       = (/llbox,llbox/)
  d0          = dsize(0)
  d1          = dsize(1)

; Choose which offline storm track to read in
  if (cn0 .eq. "slp") then
   cn = 0
  elseif (cn0 .eq. "vort") then
   cn = 1
  elseif (cn0 .eq. "geo") then
   cn = 2
  elseif (cn0 .eq. "geo_sm") then
   cn = 3
  end if

;========================================================================
; Read in required variables from 'pb' stream (add 'pd' stream shortly)
;========================================================================

; pb [u,v,w,vort,t,z,rh], pd [w,p,t,theta,kmh,tau]
  b	    = addfile(fili_b,"r")
  d  	    = addfile(fili_d,"r")

;==================================
; Get the variables we will need
;==================================

; Loop over times in each simulation
  do it  = ts0(st), tf0(st)	

 ; Print out time information for user 
   print("Working on time: "+time_str(it-4))
   title_arr(it-4) = "Valid at "+time_str(it-4)+" (T+"+it+")"

 ; Use this index as base (120 times)
   t0    = lat_arr(cn,it-1,0)
   t1    = lat_arr(cn,it-1,1)
   n0    = lon_arr(cn,it-1,0)
   n1    = lon_arr(cn,it-1,1)
   print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

 ; Combine values above into array for input into external function
   ll_arr = (/t0,t1,n0,n1/)

 ; Create latitude/longitude array subsets to reduce computation time
   lon1  = b->longitude_1({n0:n1})  ; longitude ['d0' grid points]
   lat1  = b->latitude_1({t0:t1})   ; latitude  ['d1' grid points]
   lon   = b->longitude({n0:n1})    ; longitude ['d0' grid points]
   lat   = b->latitude({t0:t1})     ; latitude  ['d1' grid points]

 ; Pressure levels (pres)
 ; (0) 1000, (1) 950, (2) 925, (3) 900, (4) 850, (5) 800, (6) 750, (7) 700, (8) 650, 
 ; (9) 600, (10) 500, (11) 400, (12) 300, (13) 250, (14) 200, (15) 150, (16) 100
   pres  = b->p                             ; Pressure (hPa)
   num_l = dimsizes(pres)
   pres@units = "hPa"

   u0    = b->u(it-4,:,{t0:t1},{n0:n1})     ; Zonal velocity (m/s)
   v0    = b->v(it-4,:,{t0:t1},{n0:n1})     ; Meridional velocity (m/s)
   vort  = b->rvor(it-4,:,{t0:t1},{n0:n1})  ; Relative vorticity (/s)

   t     = b->temp(it-4,:,{t0:t1},{n0:n1})  ; Temperature (K)
   z     = b->ht(it-4,:,{t0:t1},{n0:n1})    ; Geopotential height (m)
   w     = b->dz_dt(it-4,:,{t0:t1},{n0:n1}) ; Vertical velocity (m/s)

 ; Convert to units we want (decametres; /s)
   z     = z/10
   z@units = "dam"
   vort  = vort * (10 ^ -6)

   dx0   = lon(1) - lon(0) ; Grid spacing (longitude)
   dy0   = lat(1) - lat(0) ; Grid spacing (latitude)

 ; Calculate storm-relative winds
   u = u0 - u_cyc(cn,it-1)
   v = v0 - v_cyc(cn,it-1)

 ; Add metadata from original horizontal wind arrays
   copy_VarCoords(u0,u)
   copy_VarAtts(u0,u)

   copy_VarCoords(v0,v)
   copy_VarAtts(v0,v)

   delete([/u0,v0/])

 ; Coriolis parameter and absolute vorticity
   f1    = coriolis_param(lat)
   f0    = conform_dims(dimsizes(vort),f1,1)
   copy_VarMeta(vort,f0)
   copy_VarCoords(vort,f0)
   f0@description  = "Coriolis parameter"
   f0@name         = "Coriolis parameter"

   avo   = vort + f0
   copy_VarMeta(vort,avo)
   copy_VarAtts(vort,avo)
   copy_VarCoords(vort,avo)
   avo@description = "Absolute vorticity"
   avo@long_name   = "Vorticity (absolute)"

;==============================================================
; Make sure all arrays are same size (for later calculations)
;==============================================================

 ; 2-D array (6 * 3)
   size  = (/dimsizes(u),dimsizes(v),dimsizes(vort),dimsizes(t),dimsizes(z),dimsizes(w)/)

 ; Reduce size of any dimension larger than the rest (u,v,w,t,etc)
   do sz  = 0, 2
    nind = min( size(:,sz) )
    if ( any (size(:,sz) .gt. nind) ) then 
     size(:,sz) = nind
    end if
   end do  

 ; Resize arrays based upon results of calculations above
   u2    = u(0:size(0,0)-1,0:size(0,1)-1,0:size(0,2)-1)
   v2    = v(0:size(1,0)-1,0:size(1,1)-1,0:size(1,2)-1)
   vort2 = vort(0:size(2,0)-1,0:size(2,1)-1,0:size(2,2)-1)
   t2 	 = t(0:size(3,0)-1,0:size(3,1)-1,0:size(3,2)-1)
   z2 	 = z(0:size(4,0)-1,0:size(4,1)-1,0:size(4,2)-1)
   w2    = w(0:size(5,0)-1,0:size(5,1)-1,0:size(5,2)-1)
   lat2	 = lat(0:nind-1)
   lon2	 = lon(0:nind-1)
   delete([/u,v,vort,t,z,w,lat,lon/])

 ; Reset names to avoid further code changes 
   u     = u2
   v 	 = v2
   vort	 = vort2
   t	 = t2
   z	 = z2
   w	 = w2
   lat	 = lat2
   lon	 = lon2
   delete([/u2,v2,vort2,t2,z2,w2,lat2,lon2/])

;==========================================================
; Calculate divergence using one-sided difference scheme
;==========================================================

   div = uv2dv_cfd(u, v, lat, lon, 2)
   copy_VarCoords(u, div)
   div@units = "s~S~-1~N~"
   div@name = "Divergence on pressure levels"

;==============================================
; Get the variables we will need [pd stream]
;==============================================

 ; Read in basic variables from 'pd' stream
   lon1D  = d->longitude_1({n0:n1})                ; longitude [xxx grid points]
   lat1D  = d->latitude_1({t0:t1})                 ; latitude  [xxx grid points]
   lonD   = d->longitude({n0:n1})                  ; longitude [xxx grid points]
   latD   = d->latitude({t0:t1})                   ; latitude  [xxx grid points]
   hybD   = d->hybrid_ht(:)                        ; 63 model ('rho') levels

 ; Read in diffusivity and wind stress
 ; 'str_xz' -- [hybD | 63] * [latD | 150]  * [lonD | 151]
 ; 'str_yz' -- [hybD | 63] * [lat1D | 150] * [lon1D | 151]
 ; 'kmh'    -- [hybD | 63] * [latD | 150]  * [lon1D | 151]
   str_xz = d->taux(it-2,:,{t0:t1},{n0:n1})        ; Wind stress, x-comp (N m/2)
   str_yz = d->tauy(it-2,:,{t0:t1},{n0:n1})        ; Wind stress, y-comp (N m/2)
   kmh    = d->unspecified(it-2,:,{t0:t1},{n0:n1}) ; Horizontal diffusivity (m^2 s-1)

 ; Read in temperature, pressure, vertical velocity and geopotential height
 ; All variables below -- [hybD | 63] * [latD | 150] * [lon1D | 151]
   td     = d->temp(it-2,:,{t0:t1},{n0:n1})        ; Temperature (K)
   pd     = d->p(it-2,:,{t0:t1},{n0:n1})           ; Pressure (Pa)
   wd     = d->dz_dt(it-2,:,{t0:t1},{n0:n1})       ; Vertical velocity (m s-1)

 ; Calculate dry air density from pressure and temperature
   rho    = pd / (rd * td)

 ; Add metadata
   copy_VarCoords(pd, rho)
   rho@description = "Dry air density"
   rho@units = "kg m~S~-3~N~"

 ; Also calculate specific volume (for use in diffusive tendency source terms)
   rho_a  = 1 / rho
   copy_VarCoords(rho, rho_a)
   rho_a@description = "Specific volume"
   rho_a@units = "m~S~3~N~ kg~S~-1~N~"

   pd      = pd / 100 	     ; Convert units of pressure before continuing
   pd@units = "hPa"	     ; Add metadata

;==========================================================================
; Define new arrays for terms in vorticity equation (no friction for now)
;==========================================================================

   pii        = 3.14159265           ; Constant
   lat_rad    = lat * (pii/180)      ; Convert to radians for calculation below
   lon_rad    = lon * (pii/180)

 ; Troubleshooting (look for abnormally large values and exit the script if found)
   if( any (u .gt. 200) ) then
    print("Quitting because of storm-relative velocity error!!!")
   end if

   d_vort      = dimsizes(vort)
   d_1	       = d_vort(1)
   d_2	       = d_vort(2)

 ; Create new arrays for horizontal components of tilting term
   if (d_1 .eq. d0 .and. d_2 .eq. d1) then
    tilt_x     = new((/num_l,d0,d1/),float)
    tilt_y     = new((/num_l,d0,d1/),float)
   else
    tilt_x     = vort
    tilt_y     = vort
   end if 

 ; Add basic metadata (modify details later)
   copy_VarMeta(vort,tilt_x)
   copy_VarAtts(vort,tilt_x)
   copy_VarCoords(vort,tilt_x)

   copy_VarMeta(vort,tilt_y)
   copy_VarAtts(vort,tilt_y)
   copy_VarCoords(vort,tilt_y)

;==========================================================
; Define latitude and longitude for finite differencing
;==========================================================

 ; Create lat/lon grid with same dimensions as all input variables
   if (d_1 .eq. d0 .and. d_2 .eq. d1) then
    y2d      = new((/d0,d1/),float)	
   else 
    y2d	     = vort(0,:,:)   
   end if 

   y2d!0     = "lat"
   y2d!1     = "lon"
   y2d&lat   =  lat
   y2d&lon   =  lon

   y2d       = conform_dims(dimsizes(y2d),lat,0) ; Array of latitude values
   coslat    = cos(y2d*pii/180)                  ; Same array in radians

 ; Determine horizontal resolution of data
   delta_lon = y2d
   delta_lat = y2d

   delta_lon = lon(1) - lon(0) ; Zonal grid spacing
   delta_lat = lat(1) - lat(0) ; Meridional grid spacing

   delta_lon = sqrt(delta_lon ^ 2)
   delta_lat = sqrt(delta_lat ^ 2)

   lons      = lon(1) - lon(0)
   lats      = lat(1) - lat(0)
   lons      = sqrt(lons^2)
   lats      = sqrt(lats^2)

 ; Calculate dx and dy (changes with longitude but not with latitude)
   ddx 	     = y2d
   ddx 	     = r_earth*(pii/180.0)*delta_lon*cos(y2d*pii/180)
   ddy 	     = y2d
   ddy 	     = r_earth*(pii/180)*delta_lat

 ; Give 'coslat' the dimensions of 'u' (3-D) and hold values in new array ('coslat3d')
   coslat3d  = conform_dims(dimsizes(u),coslat,(/1,2/))

;=================================
; Friction-related calculations 
;=================================

 ; Array of hybrid levels to approximately match pressure levels in 'pres'
 ; Array [/1000, 950, 925, 900, 850, 800, 750, 700, \
 ;         650, 600, 500, 400, 300, 250, 200, 150, 100/]

   lev_arr = (/1, 8, 13, 16, 19, 21, 23, 25, 27, 29, \
               31, 33, 35, 37, 39, 42, 47, 51, 56/)
   lev_arr@units = "m"

   l0    = 0  ; Subset of vertical levels (19 levels)
   l1    = 18 ; Level 56 of original array ~ 100 hPa (Int. Standard Atmosphere)

 ; EDIT (???) -- see L385-389 in 'tc_cp_ens_xz_cbgt3_core.ncl'

;===============================================================
; Tilting term (w * dv_dp, -w * du_dp)
;===============================================================

   if (d_1 .eq. d0 .and. d_2 .eq. d1) then
    du_dp      = new((/num_l,d0,d1/),float)
    dv_dp      = new((/num_l,d0,d1/),float)
   else
    du_dp      = vort
    dv_dp      = vort
   end if

 ; Convert pressure values to Pa for calculation below
   pres0 = pres * 100
   pres0@units = "Pa"

 ; Fill arrays with vertical derivatives at each pressure level
   do pn = 0, num_l-1

  ; Account for top and bottom levels 
    kp1 = min((/pn+1,num_l-1/))
    km1 = max((/pn-1,0/))

    du_dp(pn,:,:) = (u(kp1,:,:) - u(km1,:,:)) / ( (pres0(km1))-(pres0(kp1)) )
    dv_dp(pn,:,:) = (v(kp1,:,:) - v(km1,:,:)) / ( (pres0(km1))-(pres0(kp1)) )

   end do
   delete(pres0)

 ; Calculate [dv_dp * w] -- x-component of tilting term
   tilt_y = dv_dp * w

 ; Calculate [-du_dp * w] -- y-component of tilting term
   tilt_x = -du_dp * w

 ; EDIT -- check dimensions of these arrays

 ; Metadata
   tilt_x@units = "m s~S~-2~N~"
   tilt_x@standard_name = "tilting_x"
   tilt_x@long_name = "Tilting (x-component)"
   tilt_x@title = "Tilting (x-component)"
   tilt_x@name  = "tilting_x"

   tilt_y@units = "m s~S~-2~N~"
   tilt_y@standard_name = "tilting_y"
   tilt_y@long_name = "Tilting (y-component)"
   tilt_y@title = "Tilting (y-component)"
   tilt_y@name  = "tilting_y"

;============================================================
; Loop over pressure levels to calculate vertical profile
;============================================================

 ; Counter variable (pressure)
   cp = 0		
 
 ; Loop over pressure levels
   do ilev = 0, num_l-1	

    p = pres(ilev)
    print_clock("Working on: " + p + " hPa")

;=================================
; Pressure levels ['pb' stream]
;=================================

  ; EDIT -- replace '0:d0-1,0:d1-1' with ':,:'???

  ; Read in storm-relative wind on pressure levels
    u_plane       = u(ilev,0:d0-1,0:d1-1)
    v_plane       = v(ilev,0:d0-1,0:d1-1)

  ; Read in other variables on pressure levels
    vort_plane    = vort(ilev,0:d0-1,0:d1-1)
    avo_plane	  = avo(ilev,0:d0-1,0:d1-1)
    div_plane	  = div(ilev,0:d0-1,0:d1-1)
    geo_plane     = z(ilev,0:d0-1,0:d1-1)
    t_plane       = t(ilev,0:d0-1,0:d1-1)
    w_plane       = w(ilev,0:d0-1,0:d1-1)

  ; EDIT FROM HERE -- what additional arrays do I need to run 'st_centre' successfully?
    spd_plane     = sqrt( (u_plane ^ 2) + (v_plane ^ 2) )
    copy_VarCoords(u_plane, spd_plane)
    vort_pl_smth  = smth9_Wrap(vort_plane, 0.5, 0.5, True)

;======================================================= 
; Calculate storm centre position (external function) 
;=======================================================

  ; Set radial distance (º) when looking for vort min in ringlike phase
    r0            = 0.05

  ; Call 'st_centre', replacing 'prs_plane' with 'geo_plane'
    print_clock("Calculating storm centre position...")
    centre_arr    = st_centre(ilev, it, cn, r0, centre, vort_plane, geo_plane, spd_plane, \
                              centre_new, dy0, dx0, "comp", ct, "pb", "off")    

  ; Updated storm track information contained in 'centre_new'
    centre_new    = centre_arr[0]

  ; Retrieve storm centre position at current time 
    lat_max0      = centre_new(sc,ilev,ct,0)
    lon_max0      = centre_new(sc,ilev,ct,1)    

;===================================== 
; START MID-SCRIPT PLOTTING ROUTINE 
;===================================== 

   if (plt .eq. 1) then 

   ; Plot the data as we go                                                        
     output_loop = "$sam/nepartak/images/cart_coords/xy_cbgt_"+dat+"_"+ens0+\
                   "_"+p+"_sc"+sc+"_"+time_arr(it-4)
     wks = gsn_open_wks(opt, output_loop)

     gsn_define_colormap(wks,"prcp_new")

   ; Relative vorticity resources                                                  
      opts_rv                               = True
      opts_rv@cnFillOn                      = True
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
      opts_rv@lbLabelBarOn                   = True

   ; Label bar resources (only active if lb switched on)                           
      opts_rv@lbPerimOn                     = False
      opts_rv@lbLabelFont                   = "Helvetica"
      opts_rv@lbBoxEndCapStyle              = "TriangleBothEnds"
      opts_rv@lbOrientation                 = "Vertical"

      rad0 = (rad + 0.25)

      opts_rv@mpLimitMode                   = "Corners"
      opts_rv@mpRightCornerLatF             = lat_max0 + rad0
      opts_rv@mpLeftCornerLatF              = lat_max0 - rad0
      opts_rv@mpRightCornerLonF             = lon_max0 + rad0
      opts_rv@mpLeftCornerLonF              = lon_max0 - rad0

   ; Do not draw the plot or advance the frame                                          
      opts_rv@gsnDraw                       = False
      opts_rv@gsnFrame                      = False
      opts_rv@gsnAddCyclic                  = False
      opts_rv@tmXTOn                        = False
      opts_rv@tmYROn                        = False

   ; Additional plotting resources (database; resolution; borders; colour; thickness; lat/lon) 
      opts_rv@mpDataBaseVersion             = "Ncarg4_1"
      opts_rv@mpDataSetName                 = "Earth..4"
      opts_rv@mpOutlineBoundarySets         = "National"
      opts_rv@mpGeophysicalLineColor        = "black"
      opts_rv@mpGeophysicalLineThicknessF   = 1.0
      opts_rv@mpGridAndLimbOn               = False
      opts_rv@pmTickMarkDisplayMode         = "Always"
      opts_rv@tmXBMajorLengthF              = 0.005

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
      opts_vvel@cnFillPalette               = "n11"
      opts_vvel@cnLevels                    = (/0.1, 0.2, 0.4, 0.6, 0.8, \
                                                1.0, 1.2, 1.5, 2.0, 3.0/)
      opts_vvel@cnFillColors                = (/0,2,3,4,5,6,7,8,9,10,11/)

   ; Resources for divergence                                                          
      opts_divg                              = opts_rv
      delete([/opts_divg@cnLevels, opts_divg@cnFillColors/])
      opts_divg@cnLevels                     = (/5,10,15,20,30,40,60,80,100,150/)
      opts_divg@cnFillColors                 = (/0,4,5,6,7,8,9,10,11,12,13,14,15/)

   ; Wind vector resources
      opts_vec                            = True
      opts_vec@gsnDraw                    = False
      opts_vec@gsnFrame                   = False
      opts_vec@vcGlyphStyle               = "LineArrow"
      opts_vec@vcRefAnnoOn                = False
      opts_vec@vcRefMagnitudeF            = 20
      opts_vec@vcRefLengthF               = 0.04
      opts_vec@vcMinDistanceF             = 0.05
      opts_vec@vcLineArrowThicknessF      = 1.0
      opts_vec@vcLineArrowColor           = "black"
      opts_vec@gsnMaximize                = True
      opts_vec@gsnRightString             = ""
      opts_vec@gsnLeftString              = ""
      opts_vec@tiXAxisString              = ""
      opts_vec@vcRefAnnoOrthogonalPosF    = 0.05

   ; Panel plot label resources                                                         
      txid                                   = new(4,graphic)
      amid                                   = new(4,graphic)
      panel                                  = new(4,graphic)

      panel_strings                          = (/ "a) Relative vorticity",\
                                                  "b) Vertical velocity",\
                                                  "c) Vector wind",\
                                                  "d) Divergence"/)

      txres                                  = True
      txres@txPerimOn                        = True
      txres@txFontHeightF                    = 0.0125
      txres@txBackgroundFillColor            = "White"

      amres                                  = True
      amres@amParallelPosF                   = -0.48
      amres@amOrthogonalPosF                 = -0.48
      amres@amJust                           = "TopLeft"

   ; Coordinates of panels                                                        
      xf = (/0.05, 0.55, 0.05, 0.55/)
      yf = (/0.90, 0.90, 0.49, 0.49/)
      wf = 0.40
      hf = 0.35

;==============================
; Overlay box onto x-y plots 
;==============================

    ; Box coordinates 
      lat_d0 = lat_max0 - rad
      lon_d0 = lon_max0 - rad
      lat_u0 = lat_max0 + rad
      lon_u0 = lon_max0 + rad

      lat_box = new(5,float)
      lon_box = new(5,float)      

      lat_box  = (/lat_d0, lat_d0, lat_u0, lat_u0, lat_d0/)
      lon_box  = (/lon_d0, lon_u0, lon_u0, lon_d0, lon_d0/)
   
   ; Panel 1 (relative vorticity)                                                 
      opts_rv@vpXF       = xf(0)
      opts_rv@vpYF       = yf(0)
      opts_rv@vpWidthF   = wf
      opts_rv@vpHeightF  = hf

    ; Change units before plotting
      vort_plane0        = vort_plane * (10 ^ 4)
      copy_VarCoords(vort_plane, vort_plane0)

      vec_0              = gsn_csm_vector(wks, u_plane, v_plane, opts_vec)
      vort_0             = gsn_csm_contour_map(wks, vort_plane0, opts_rv)
      overlay(vort_0, vec_0)
      panel(0)           = vort_0

   ; Overlay labels                                                               
      txid(0)      = gsn_create_text(wks, panel_strings(0), txres)
      amid(0)      = gsn_add_annotation(panel(0), txid(0), amres)

   ; Panel 2 (vertical velocity)                                                  
      opts_vvel@vpXF       = xf(1)
      opts_vvel@vpYF       = yf(1)
      opts_vvel@vpWidthF   = wf
      opts_vvel@vpHeightF  = hf

      vec_1                = gsn_csm_vector(wks, u_plane,v_plane, opts_vec)
      vvel_1               = gsn_csm_contour_map(wks, w_plane(:,:), opts_vvel)
      overlay(vvel_1, vec_1)
      panel(1)             = vvel_1

   ; Overlay labels                                                               
      txid(1)      = gsn_create_text(wks, panel_strings(1), txres)
      amid(1)      = gsn_add_annotation(panel(1), txid(1), amres)

   ; Panel 3 (total wind)                                                         
      opts_wind@vpXF       = xf(2)
      opts_wind@vpYF       = yf(2)
      opts_wind@vpWidthF   = wf
      opts_wind@vpHeightF  = hf

      vec_2                = gsn_csm_vector(wks, u_plane,v_plane, opts_vec)
      wind_2               = gsn_csm_contour_map(wks, spd_plane(:,:), opts_wind)
      overlay(wind_2, vec_2)
      panel(2)             = wind_2

   ; Overlay labels                                                               
      txid(2)      = gsn_create_text(wks, panel_strings(2), txres)
      amid(2)      = gsn_add_annotation(panel(2), txid(2), amres)

   ; Panel 4 (divergence)                                                        
      opts_divg0           = opts_divg
      opts_divg0@vpXF      = xf(3)
      opts_divg0@vpYF      = yf(3)
      opts_divg0@vpWidthF  = wf
      opts_divg0@vpHeightF = hf

    ; Change units before plotting 
      div_plane0           = div_plane * (10 ^ 5)
      copy_VarCoords(div_plane, div_plane0)

      vec_3                = gsn_csm_vector(wks, u_plane,v_plane, opts_vec)
      div_3                = gsn_csm_contour_map(wks, div_plane0(:,:), opts_divg0)
      overlay(div_3, vec_3)
      panel(3)             = div_3

   ; Overlay labels                                                               
      txid(3)      = gsn_create_text(wks, panel_strings(3), txres)
      amid(3)      = gsn_add_annotation(panel(3), txid(3), amres)

;====================================
; Overlay box onto all panel plots 
;====================================

      dres     = True
      dres@gsLineThicknessF  = 3.0
      dres@gsLineDashPattern = 1
      dres@gsLineColor       = "black"

      lines    = new(4,graphic)

      do n = 0, dimsizes(panel)-1
       lines(n) = gsn_add_polyline(wks,panel(n),lon_box,lat_box,dres)
      end do 

;==================================================
; Overlay markers for each storm centre estimate 
;==================================================

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
    markers(0,i)             = gsn_add_polymarker(wks, panel(0), centre_new(i,ilev,it-4,1),\
                                                  centre_new(i,ilev,it-4,0), mres0)
    markers(1,i)             = gsn_add_polymarker(wks, panel(1), centre_new(i,ilev,it-4,1),\
                                                  centre_new(i,ilev,it-4,0), mres0)
    markers(2,i)             = gsn_add_polymarker(wks, panel(2), centre_new(i,ilev,it-4,1),\
                                                  centre_new(i,ilev,it-4,0), mres0)
    markers(3,i)             = gsn_add_polymarker(wks, panel(3), centre_new(i,ilev,it-4,1),\
                                                  centre_new(i,ilev,it-4,0), mres0)
  end do

   ; Tidy up before final plots later in script
     delete([/opts_vvel, opts_wind, opts_rv, opts_divg, panel, \
              txid, amid, panel_strings, lat_box, lon_box, lat_d0, lon_d0, lat_u0, lon_u0/])

;===============================================
; Finally, create a title and plot everything
;===============================================

    ; Create a main title at the top of the page                                               
      title_res                   = True
      title_str                   = "Panel plot at "+p+" hPa: "+ens0+", "+dat+" (T+"+it+")"
      title_res@txFontHeightF     = 0.02
      title_res@txFont            = "Helvetica"
      xpos_title                  = 0.50
      ypos_title                  = 0.94

      gsn_text_ndc(wks, title_str, xpos_title, ypos_title, title_res)

    ; Call 'draw' and 'frame'                                                                  
      maximize_output(wks, True)

    end if 

;====================================  
; END MID-SCRIPT PLOTTING ROUTINE 
;====================================  

  ; Create temporary variables to make circulation budget calculation less messy
  ; 'tilt_x' and 'tilt_y' contain values on all pressure levels for each time
  ; 'tilt_x_plane' and 'tilt_y_plane' contain values on single pressure levels
    tilt_x_plane  = tilt_x(ilev,:,:)
    tilt_y_plane  = tilt_y(ilev,:,:)

;=====================================================================================
; Create array of box positions either side of the original [Davis & Galarneau 2009]
;===================================================================================== 

  ; Get information about box size for later ensemble calculations
    box       = ddx( {lat_max0-rad:lat_max0+rad}, {lon_max0-rad:lon_max0+rad} )
    box_size  = dimsizes(box)
    size0     = box_size(0)
    size1     = box_size(1)

    s0        = size0-1
    s1        = size1-1

   if (nbox .eq. 49) then
    lat_box   = (/lat_max0-(3*dy0),lat_max0-(2*dy0),lat_max0-(dy0),lat_max0,\
                  lat_max0+(dy0),lat_max0+(2*dy0),lat_max0+(3*dy0)/)
    lon_box   = (/lon_max0-(3*dx0),lon_max0-(2*dx0),lon_max0-(dx0),lon_max0,\
                  lon_max0+(dx0),lon_max0+(2*dx0),lon_max0+(3*dx0)/)
   elseif (nbox .eq. 25) then
    lat_box   = (/lat_max0-(2*dy0),lat_max0-(dy0),lat_max0,\
   	          lat_max0+(dy0),lat_max0+(2*dy0)/)
    lon_box   = (/lon_max0-(2*dx0),lon_max0-(dx0),lon_max0,\
                  lon_max0+(dx0),lon_max0+(2*dx0)/)
   elseif (nbox .eq. 9) then
    lat_box   = (/lat_max0-(dy0),lat_max0,lat_max0+(dy0)/)
    lon_box   = (/lon_max0-(dx0),lon_max0,lon_max0+(dx0)/)
   elseif (nbox .eq. 1) then
    lat_box   = (/lat_max0/)
    lon_box   = (/lon_max0/)
   end if

   lsize      = dimsizes(lat_box) * dimsizes(lon_box)

;====================================================
; Define arrays to hold ensemble of box positions
;====================================================

 ; EDIT --> NEATER WAY TO INITIALISE MULTIPLE ARRAYS (???)
 ; Ensemble of boxes (all grid points, not just box centres)
   ddx_box    = new((/lsize,size0,size1/),float)
   ddy_box    = new((/lsize,size0,size1/),float)
   tilt_x_box = new((/lsize,size0,size1/),float)
   tilt_y_box = new((/lsize,size0,size1/),float)
   vort_box   = new((/lsize,size0,size1/),float)
   avo_box    = new((/lsize,size0,size1/),float)
   div_box    = new((/lsize,size0,size1/),float)
   u_box      = new((/lsize,size0,size1/),float)
   v_box      = new((/lsize,size0,size1/),float)

   lat_d      = new(lsize,float)
   lon_d      = new(lsize,float)
   lat_u      = new(lsize,float)
   lon_u      = new(lsize,float)

   d_ave      = new(lsize,float)
   dz0        = new(lsize,float)
   dz1        = new(lsize,float)
   div_ave    = new(lsize,float)
   vort_ave   = new(lsize,float)

 ; Initialise array to hold perimeter averages (for each side of box)
   vort_avg   = new((/lsize,5/),float)
   avo_avg    = new((/lsize,5/),float)
   u_avg      = new((/lsize,5/),float)
   v_avg      = new((/lsize,5/),float)

;=============================================
; Calculate terms in the circulation budget
;=============================================

 ; Radius of box for calculation (degrees)
   r          = rad
   r1 	      = sprintf("%0.1fdeg",2*r) ; Create string to add to output file

;===================================================================
; Loop over perturbed box positions, either side of the TC centre
;===================================================================

   do i = 0, dimsizes(lon_box)-1    ; Longitude
    do j = 0, dimsizes(lat_box)-1   ; Latitude

     a	  = (dimsizes(lon_box) * i) + j
     print_clock("Working on box number "+a)

     ddx_box(a,:,:)    = ddx({lat_box(j)-r:lat_box(j)+r},{lon_box(i)-r:lon_box(i)+r})
     ddy_box(a,:,:)    = ddy({lat_box(j)-r:lat_box(j)+r},{lon_box(i)-r:lon_box(i)+r})
     tilt_x_box(a,:,:) = tilt_x_plane({lat_box(j)-r:lat_box(j)+r},\
     		      	              {lon_box(i)-r:lon_box(i)+r})
     tilt_y_box(a,:,:) = tilt_y_plane({lat_box(j)-r:lat_box(j)+r},\
				      {lon_box(i)-r:lon_box(i)+r})

     vort_box(a,:,:)   = vort_plane({lat_box(j)-r:lat_box(j)+r},\
				    {lon_box(i)-r:lon_box(i)+r})
     avo_box(a,:,:)    = avo_plane({lat_box(j)-r:lat_box(j)+r},\
      		      	           {lon_box(i)-r:lon_box(i)+r})
     div_box(a,:,:)    = div_plane({lat_box(j)-r:lat_box(j)+r},\
				   {lon_box(i)-r:lon_box(i)+r})
     u_box(a,:,:)      = u_plane({lat_box(j)-r:lat_box(j)+r},\
     		      	         {lon_box(i)-r:lon_box(i)+r})
     v_box(a,:,:)      = v_plane({lat_box(j)-r:lat_box(j)+r},\
			 	 {lon_box(i)-r:lon_box(i)+r})

   ; EDIT -- this is where additional code may be required
   ; Box edge coordinates (used earlier to overlay box onto panel plots)
     lat_d(a) = lat_box(j)-r
     lon_d(a) = lon_box(i)-r
     lat_u(a) = lat_box(j)+r
     lon_u(a) = lon_box(i)+r

;=====================================================
; Calculate distances around each perimeter (m)
;=====================================================

   ; Meridional boundaries (same both sides)
     dm0 = r_earth * ((pii/180) * (2*r))
     dm0@units = "m"

   ; 'Bottom' boundary
     dz0(a) = r_earth * ((pii/180) * (2*r)) * cos((pii/180)*lat_d(a))
     dz0@units = "m"

   ; 'Top' boundary
     dz1(a) = r_earth * ((pii/180) * (2*r)) * cos((pii/180)*lat_u(a))
     dz1@units = "m"

   ; Use approximation to calculate the area of the box
     d_ave(a) = ( (dz0(a) + dz1(a) ) / 2) * dm0
     d_ave@units = "m~S~-2~N~"

   ; Also calculate area of box using 'gc_latlon' function
     ddz0 = gc_latlon( lat_d(a), lon_d(a), lat_d(a), lon_u(a), 2, 3)	; Bottom
     ddz1 = gc_latlon( lat_d(a), lon_u(a), lat_u(a), lon_u(a), 2, 3)	; Right
     ddz2 = gc_latlon( lat_u(a), lon_u(a), lat_u(a), lon_d(a), 2, 3)	; Top	
     ddz3 = gc_latlon( lat_u(a), lon_d(a), lat_d(a), lon_d(a), 2, 3)	; Left
     aaz0 = ( (ddz0 + ddz2) / 2) * ddz1	   	     	       	        ; Area

   ; Troubleshooting 
     if (tr .eq. 1) then 
      print("Right boundary length = "+dm0+" "+dm0@units+" ("+ddz1+" m)")     
      print("Left boundary length = "+dm0+" "+dm0@units+" ("+ddz3+" m)")
      print("Bottom boundary length = "+dz0+" "+dz0@units+" ("+ddz0+" m)")
      print("Top boundary length = "+dz1+" "+dz1@units+" ("+ddz2+" m)")
      print("Area of box = "+d_ave+" "+d_ave@units+" ("+aaz0+" m)")
     end if 

   ; Calculate divergence averaged over the box
     div_ave(a) = wgt_areaave(div_box(a,:,:), 1.0, 1.0, 0) ; Divergence
     div_ave@units = "s~S~-1~N~"

   ; Calculate area-averaged vorticity 
     vort_ave(a) = wgt_areaave(vort_box(a,:,:), 1.0, 1.0, 0) ; Relative vorticity
     vort_ave@units = "s~S~-1~N~"

   ; Number of grid points along each side of box surrounding cyclone
     numPOINTS_x0 = dimsizes(ddx_box(a,:,:))
     numPOINTS_x  = numPOINTS_x0(1)   ; Longitude points
     numPOINTS_y  = numPOINTS_x0(0)   ; Latitude points

   ; Troubleshootin
     if (tr .eq. 1) then
      print("Right boundary length = "+dm0+" "+dm0@units+" ("+ddz1+" m)")
      print("Left boundary length = "+dm0+" "+dm0@units+" ("+ddz3+" m)")
      print("Bottom boundary length = "+dz0+" "+dz0@units+" ("+ddz0+" m)")
      print("Top boundary length = "+dz1+" "+dz1@units+" ("+ddz2+" m)")
      print("Area of box = "+d_ave+" "+d_ave@units+" ("+aaz0+" m)")
      print("Area-averaged divergence = "+div_ave(a)+" ("+div_ave@units+")")
      print("Area-averaged vorticity = "+vort_ave(a)+" ("+vort_ave@units+")")
     end if

;===========================================
; Calculate averages along each perimeter
;===========================================

   ; Now calculate terms along the four edges of the box, for each ensemble member
   ; Hence four arrays required for each term in the budget

   ; All arrays have dimensions [lsize * 13] -> ALL members on a SINGLE pressure level

   ; South 
     vort_avg(a,0) = avg(vort_box(a,0,:))
     avo_avg(a,0)  = avg(avo_box(a,0,:))
     u_avg(a,0)	   = avg(u_box(a,0,:))
     v_avg(a,0)	   = avg(v_box(a,0,:))

   ; East
     vort_avg(a,1) = avg(vort_box(a,:,numPOINTS_x-1))
     avo_avg(a,1)  = avg(avo_box(a,:,numPOINTS_x-1))
     u_avg(a,1)    = avg(u_box(a,:,numPOINTS_x-1))
     v_avg(a,1)    = avg(v_box(a,:,numPOINTS_x-1))

   ; North
     vort_avg(a,2) = avg(vort_box(a,numPOINTS_y-1,:))
     avo_avg(a,2)  = avg(avo_box(a,numPOINTS_y-1,:))
     u_avg(a,2)    = avg(u_box(a,numPOINTS_y-1,:))
     v_avg(a,2)    = avg(v_box(a,numPOINTS_y-1,:))

   ; West 
     vort_avg(a,3) = avg(vort_box(a,:,0))
     avo_avg(a,3)  = avg(avo_box(a,:,0))
     u_avg(a,3)    = avg(u_box(a,:,0))
     v_avg(a,3)    = avg(v_box(a,:,0))

     if (tr .eq. 1) then 
      print("Average vort: "+vort_avg(a,0)+" (S), "+vort_avg(a,1)+\
      		           " (N), "+vort_avg(a,2)+" (E), "+vort_avg(a,3)+" (W)")
      print("Average avo: "+avo_avg(a,0)+" (S), "+avo_avg(a,1)+\
      		    	  " (N), "+avo_avg(a,2)+" (E), "+avo_avg(a,3)+" (W)")	
      print("Average u: "+u_avg(a,0)+" (S), "+u_avg(a,1)+\
      		        " (N), "+u_avg(a,2)+" (E), "+u_avg(a,3)+" (W)")
      print("Average v: "+v_avg(a,0)+" (S), "+v_avg(a,1)+\
      		        " (N), "+v_avg(a,2)+" (E), "+v_avg(a,3)+" (W)")
     end if

   ; Also calculate average around the whole box
     vort_avg(a,4) = avg(vort_avg(a,0:3) )
     avo_avg(a,4)  = avg(avo_avg(a,0:3) )
     u_avg(a,4)    = avg(u_avg(a,0:3) )
     v_avg(a,4)    = avg(v_avg(a,0:3) )

;====================================================================================
; Second round of calculations (Horizontal advection, relative vorticity, tilting)
;====================================================================================

   ; Calculate terms on northern/southern boundaries

     if (int .eq. 1) then 

    ; Horizontal advection (into the box)
      hadv_s      = dim_sum_n_Wrap((-(avo_box(a,0,:) - avo_avg(a,0)) * \
                                   (v_box(a,0,:) - v_avg(a,0)) * \
                                    ddx_box(a,0,:)), 0)
      hadv_n	  = dim_sum_n_Wrap((-(avo_box(a,numPOINTS_y-1,:) - avo_avg(a,2)) * \
                                   (v_box(a,numPOINTS_y-1,:) - v_avg(a,2)) * \
                                    ddx_box(a,numPOINTS_y-1,:)), 0)

    ; Tilting
      tilt_x0     = dim_sum_n_Wrap((tilt_x_box(a,:,:) * ddx_box(a,:,:)),1)
      tilt_s      = tilt_x0(0)              ; Tilting (south)
      tilt_n      = tilt_x0(numPOINTS_y-1)  ; Tilting (north)

    ; Circulation
      circ_x0     = dim_sum_n_Wrap((u_box(a,:,:) * ddx_box(a,:,:)),1)
      circ_s      = circ_x0(0)              ; Circulation (south)
      circ_n      = -circ_x0(numPOINTS_y-1) ; Circulation (north)

     elseif (int .eq. 0) then 

    ; Tilting 
      tilt_x0     = dim_sum_n_Wrap((tilt_x_box(a,:,:) * ddx_box(a,:,:)),1)

    ; Horizontal advection (into the box, across the perimeter)
      hadv_x0     = dim_sum_n_Wrap((-(avo_box(a,:,:) - avo_avg(a,4)) * \
                                   (v_box(a,:,:) - v_avg(a,4)) * \
                                    ddx_box(a,:,:)),1)
      circ_x0     = dim_sum_n_Wrap((u_box(a,:,:) * ddx_box(a,:,:)),1)
      tilt_s      = tilt_x0(0)              ; Tilting (south)
      tilt_n      = tilt_x0(numPOINTS_y-1)  ; Tilting (north)
      hadv_s      = hadv_x0(0)              ; Horizontal advection (south)
      hadv_n      = -hadv_x0(numPOINTS_y-1) ; Horizontal advection (north)
      circ_s      = circ_x0(0)              ; Circulation (south)
      circ_n      = -circ_x0(numPOINTS_y-1) ; Circulation (north)	
      delete(hadv_x0)

     end if 

   ; Calculate terms on eastern/western boundaries

     if (int .eq. 1) then

    ; Horizontal advection (into the box)
      hadv_w      = dim_sum_n_Wrap((-(avo_box(a,:,0) - avo_avg(a,3)) * \
                                   (v_box(a,:,0) - v_avg(a,3)) * \
                                    ddx_box(a,:,0)), 0)
      hadv_e      = dim_sum_n_Wrap((-(avo_box(a,:,numPOINTS_y-1) - avo_avg(a,1)) * \
                                   (v_box(a,:,numPOINTS_y-1) - v_avg(a,1)) * \
                                    ddx_box(a,:,numPOINTS_y-1)), 0)
   
    ; Tilting 
      tilt_y0     = dim_sum_n_Wrap((tilt_y_box(a,:,:) * ddy_box(a,:,:)),0)
      tilt_w      = tilt_y0(0)              ; Tilting (west)
      tilt_e      = tilt_y0(numPOINTS_x-1)  ; Tilting (east)

    ; Circulation
      circ_y0     = dim_sum_n_Wrap((v_box(a,:,:) * ddy_box(a,:,:)),0)
      circ_w      = -circ_y0(0)             ; Circulation (west)
      circ_e      = circ_y0(numPOINTS_x-1)  ; Circulation (east)

     elseif (int .eq. 0) then 

      tilt_y0     = dim_sum_n_Wrap((tilt_y_box(a,:,:) * ddy_box(a,:,:)),0)
      hadv_y0     = dim_sum_n_Wrap((-(avo_box(a,:,:) - avo_avg(a,4)) * \
                                   (u_box(a,:,:) - u_avg(a,4)) * \
                                    ddy_box(a,:,:)),0)
      circ_y0     = dim_sum_n_Wrap((v_box(a,:,:) * ddy_box(a,:,:)),0)
      tilt_w      = tilt_y0(0)              ; Tilting (west)
      tilt_e      = tilt_y0(numPOINTS_x-1)  ; Tilting (east)
      hadv_w      = hadv_y0(0)              ; Horizontal advection (west)
      hadv_e      = -hadv_y0(numPOINTS_x-1) ; Horizontal advection (east)
      circ_w      = -circ_y0(0)             ; Circulation (west)
      circ_e      = circ_y0(numPOINTS_x-1)  ; Circulation (east)
      delete(hadv_y0)

     end if 

;========================
; Further calculations
;========================

   ; Calculate stretching term using averaged quantities
     str_tenT   = -(div_ave(a) * avo_avg(a,4))

   ; Calculate horizontal advection from line integrals (4 sides of box)
     hadv_tenT  = hadv_s + hadv_n + hadv_w + hadv_e
     hadv_tenT  = hadv_tenT / d_ave(a) ; Divide by area of box [Rios-Berrios et al. 2016]

   ; Calculate tilting tendency from line integrals (4 sides of box)
     tilt_tenT  = tilt_s + tilt_n + tilt_w + tilt_e
     tilt_tenT  = tilt_tenT / d_ave(a) ; Divide by area of box [Rios-Berrios et al. 2016]

   ; Calculate circulation from line integrals (4 sides of box)
     circ_instT0 = circ_s + circ_n + circ_w + circ_e	     ; Sum 4 line integrals
     circ_instT0 = circ_instT0 / d_ave(a)	             ; Divide by area of box 
     circ_instT  = vort_ave(a) 				     ; Area-averaged vorticity

     if (tr .eq. 1) then 
      print("Circulation = "+circ_instT0+" (v); "+circ_instT+" (vort)")
     end if

   ; Enter values into predefined arrays
     hadv_ten(a,ct,cp)      = hadv_tenT		     ; Horizontal advection
     str_ten(a,ct,cp)       = str_tenT		     ; Stretching
     tilt_ten(a,ct,cp)      = tilt_tenT		     ; Tilting
     circ_inst(a,ct,cp)     = circ_instT             ; Circulation (LHS)

   ; Tidy up (delete variables from within lat/lon loops before next iteration)
     delete([/hadv_tenT, str_tenT, tilt_tenT, circ_instT, circ_instT0/])
     delete([/numPOINTS_x0, numPOINTS_x, numPOINTS_y/])
     delete([/tilt_y0,tilt_x0,circ_y0,circ_x0/])

     end do     ; End latitude loop (do j = 0, dimsizes(lat_box)-1)
    end do	; End longitude loop (do i = 0, dimsizes(lon_box)-1)

    cp                   = cp + 1       ; Counter variable (pressure)

  ; Tidy up (delete variables from within pressure level loop)
    delete([/u_plane, v_plane, tilt_x_plane, tilt_y_plane/])
    delete([/vort_plane, avo_plane, div_plane/])
    delete([/ddx_box, ddy_box, u_box, v_box, box/])
    delete([/tilt_x_box, tilt_y_box, vort_box, avo_box, div_box, lat_box, lon_box/])
    delete([/div_ave,vort_ave,vort_avg,avo_avg,u_avg,v_avg/])

   end do ; End pressure level loop (do ilev = 0, num_l-1)

 ; Tidy up (delete variables from within time loop before next iteration)
   delete([/lon_rad,lat_rad,lon,lat,lon1,lat1,u,v,vort,div,f0,f1,avo,t,z,w/])
   delete([/lon1D,lat1D,lonD,latD,hybD,str_xz,str_yz,kmh,td,p,pd,wd,rho,rho_a/])

   ct = ct + 1

 ; End time loop (do it = 0, times-1)
   end do

   delete([/time_b, time_d/])
   delete([/utc_date,year,month,day,hour,minute,second,date_str,time_str,time_arr,title_arr/])

; End loop over simulations (do st = 0, ...)
  end do

;=====================================================
; Calculate term averages over all ensemble members
;=====================================================

   print_clock("Finished ensemble calculations! Now averaging over all members...")

   hadv_final = new((/numINPUT_b,num_l/),"float")
   str_final  = new((/numINPUT_b,num_l/),"float")
   tilt_final = new((/numINPUT_b,num_l/),"float")
   circ_final = new((/numINPUT_b,num_l/),"float")
   circ_ten   = new((/numINPUT_b,num_l/),"float")

 ; Average over box calculations (no effect if 'box=1')
   hadv_final(:,:) = dim_avg_n_Wrap(hadv_ten(:,:,:),0)
   str_final(:,:)  = dim_avg_n_Wrap(str_ten(:,:,:),0)
   tilt_final(:,:) = dim_avg_n_Wrap(tilt_ten(:,:,:),0)
   circ_final(:,:) = dim_avg_n_Wrap(circ_inst(:,:,:),0)

;=========================================
; Integrate over Nepartak's life cycle
;=========================================

 ; Each calculation produces an array of [num_l] values [pressure levels]
 ; We do this calculation [numINPUT_b] times
 ; Create arrays with [numINPUT_b] * [num_l] elements

   hadv_intA = new((/numINPUT_b,num_l/),float)
   str_intA  = new((/numINPUT_b,num_l/),float)
   tilt_intA = new((/numINPUT_b,num_l/),float)
   tot_intA  = new((/numINPUT_b,num_l/),float)
   circ_intA = new((/numINPUT_b,num_l/),float)

 ; Create array of time intervals
   time_int = new(numINPUT_b,"integer")
   time_int = ispan(0,numINPUT_b-1,1)
   time_int!0 = "time"

   print_clock("Starting plots...")

   do ic = 0, numINPUT_b-1  ; 1 h integrations over the whole of Nepartak's life cycle

    cp1 = min((/ic+1,numINPUT_b-1/))
    cm1 = max((/ic-1,0/))

     f_str = sprinti("%0.3i",ic)+"h"    ; Formatted string for output file

     hadv_intA(ic,:) = dim_sum_n_Wrap((hadv_final(ic:cp1,:) * 1), 0)
     str_intA(ic,:)  = dim_sum_n_Wrap((str_final(ic:cp1,:) * 1), 0)
     tilt_intA(ic,:) = dim_sum_n_Wrap((tilt_final(ic:cp1,:) * 1), 0)
     circ_intA(ic,:) = (circ_final(ic,:) - circ_final(cm1,:)) / 3600

     hadv_intA(ic,:) = hadv_intA(ic,:) * (10 ^ fct)
     str_intA(ic,:)  = str_intA(ic,:) * (10 ^ fct)
     tilt_intA(ic,:) = tilt_intA(ic,:) * (10 ^ fct)
     circ_intA(ic,:) = circ_intA(ic,:) * (10 ^ fct)

   end do               ; End time loop (do ic = 0, numINPUT_b-1)

 ; Add metadata after calculating tendencies 
   hadv_intA@units = "10~S~-"+fct+"~N~ s~S~-2~N~"
   str_intA@units  = "10~S~-"+fct+"~N~ s~S~-2~N~"
   tilt_intA@units = "10~S~-"+fct+"~N~ s~S~-2~N~"
   tot_intA@units  = "10~S~-"+fct+"~N~ s~S~-2~N~"
   circ_intA@units = "10~S~-"+fct+"~N~ s~S~-2~N~"

   tot_intA  = str_intA + hadv_intA + tilt_intA
   
 ; Add metadata before time-height plotting
   circ_intA!0    = "time"
   circ_intA&time = time_int
   circ_intA!1    = "lev"
   circ_intA&lev  = pres

   tot_intA!0     = "time"
   tot_intA&time  = time_int
   tot_intA!1     = "lev"
   tot_intA&lev   = pres

   hadv_intA!0    = "time"
   hadv_intA&time = time_int
   hadv_intA!1    = "lev"
   hadv_intA&lev  = pres

   str_intA!0     = "time"
   str_intA&time  = time_int
   str_intA!1     = "lev"
   str_intA&lev   = pres

   tilt_intA!0    = "time"
   tilt_intA&time = time_int
   tilt_intA!1    = "lev"
   tilt_intA&lev  = pres

;======================================
; Reorder dimensions before plotting
;======================================

 ; Define new arrays for final plotting 
   circ_int = new((/num_l,numINPUT_b/),float)
   hadv_int = new((/num_l,numINPUT_b/),float)
   str_int  = new((/num_l,numINPUT_b/),float)
   tilt_int = new((/num_l,numINPUT_b/),float)
   tot_int  = new((/num_l,numINPUT_b/),float)

 ; Reorder from [time * pressure] to [pressure * time]
   circ_int(:,:) = circ_intA(lev|:,time|:)
   hadv_int(:,:) = hadv_intA(lev|:,time|:)
   str_int(:,:)  = str_intA(lev|:,time|:)
   tilt_int(:,:) = tilt_intA(lev|:,time|:)
   tot_int(:,:)  = tot_intA(lev|:,time|:)

 ; Calculate running average for RHS 
   hadv_int(:,:) = runave_n_Wrap(hadv_int, 3, 0, 1)
   str_int(:,:)  = runave_n_Wrap(str_int, 3, 0, 1)
   tilt_int(:,:) = runave_n_Wrap(tilt_int, 3, 0, 1)
   tot_int(:,:)  = runave_n_Wrap(tot_int, 3, 0, 1)
   
;========================
; Options for plotting 
;========================

     output = "$sam/nepartak/images/cart_coords/cbgt_xz_paper_"+dat+"_"+ens0+"_"+r1+"_"+cn0
     wks    = gsn_open_wks(opt,output)

     panel = new(5,graphic)
     gsn_define_colormap(wks,"BlueDarkRed18")

   ; Circulation budget terms
     opts_bdgt                             = True
     opts_bdgt@cnFillOn                    = True
     opts_bdgt@cnLineLabelInterval         = 2.0
     opts_bdgt@cnLineLabelFontHeightF      = 0.012
     opts_bdgt@cnLineLabelBackgroundColor  = "transparent"
     opts_bdgt@cnLineLabelPlacementMode    = "constant"
     opts_bdgt@cnLinesOn                   = False
     opts_bdgt@cnInfoLabelOn               = False
     opts_bdgt@cnLevelSelectionMode        = "ExplicitLevels"

     if (fct .eq. 8 .and. r .eq. 1.0) then
      opts_bdgt@cnLevels                    = (/-2.0, -1.5, -1.0, -0.5, \
                                                -0.3, -0.2, -0.1, -0.05, 0.05, \
                                                 0.1, 0.2, 0.3, 0.5, 1.0, 1.5, 2.0/)
     elseif (fct .eq. 8 .and. r .eq. 0.76) then
      opts_bdgt@cnLevels                    = (/-3.0, -2.0, -1.5, -1.0, \
                                                -0.5, -0.2, -0.1, -0.05, 0.05, \
                                                 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0/)
     elseif (fct .eq. 8 .and. r .eq. 0.48) then
      opts_bdgt@cnLevels                    = (/-3.0, -2.0, -1.5, -1.0, \
                                                -0.5, -0.2, -0.1, -0.05, 0.05, \
                                                 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0/)
     elseif (fct .eq. 8 .and. r .eq. 0.24) then
      opts_bdgt@cnLevels                    = (/-4.0, -3.0, -2.0, -1.5, \
                                                -1.0, -0.5, -0.2, -0.1, 0.1, \
                                                 0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0/)
     elseif (fct .eq. 9 .and. op .eq. 0) then
      opts_bdgt@cnLevels                    = (/-30.0, -25.0, -20.0, -15.0, -10.0, \
                                                -5.0, -2.0, -1.0, 1.0, 2.0, \
                                                 5.0, 10.0, 15.0, 20.0, 25.0, 30.0/)
     elseif (fct .eq. 9 .and. op .eq. 1) then
      opts_bdgt@cnLevels                    = (/-25.0, -20.0, -15.0, -10.0, -5.0, \
                                                -2.0, -1.0, -0.5, 0.5, 1.0, \
                                                 2.0, 5.0, 10.0, 15.0, 20.0, 25.0/)
     end if

     opts_bdgt@cnFillColors                = (/3,4,5,6,7,8,9,10,0,11,12,\
                                               13,14,15,16,18,19/)

     opts_bdgt@gsnPaperOrientation         = "landscape"
     opts_bdgt@tiMainFontHeightF           = 0.0125
     opts_bdgt@gsnLeftString               = ""
     opts_bdgt@gsnRightString              = ""

     opts_bdgt@gsnMaximize                 = True
     opts_bdgt@tmXTOn                      = "False" ; Turn off top x-axis TM
     opts_bdgt@tmYROn                      = "False" ; Turn off right y-axis TM
     opts_bdgt@tiYAxisString               = "Pressure (hPa)"
     opts_bdgt@tiYAxisFontHeightF          = 0.030   ; Larger y-axis title
     opts_bdgt@tiXAxisString               = "Time"  ; X-axis title
     opts_bdgt@tiXAxisFontHeightF          = 0.030   ; Larger x-axis title
     opts_bdgt@vpWidthF                    = 0.6     ; Plot width
     opts_bdgt@vpHeightF                   = 0.5     ; Plot height
     opts_bdgt@tmXBMode                    = "Explicit"
     opts_bdgt@tmXBLabelFontHeightF        = 0.025

   ; Create labels based on length of integration
   ; EDIT -- so length of integration can be smaller than currently 
     tt0 = tf0 - ts0
     if (tt0 .eq. 6) then 
      opts_bdgt@tmXBValues                  = (/ts0, ts0+3, tf0/)
     elseif (tt0 .eq. 12) then 
      opts_bdgt@tmXBValues                  = (/ts0, ts0+6, tf0/)
     elseif (tt0 .eq. 24) then 
      opts_bdgt@tmXBValues                  = (/ts0, ts0+12, tf0/)
     elseif (tt0 .eq. 36) then 
      opts_bdgt@tmXBValues                  = (/ts0, ts0+12, ts0+24, tf0/)
     elseif (tt0 .eq. 48) 
      opts_bdgt@tmXBValues                  = (/ts0, ts0+12, ts0+24, ts0+36, tf0/)
     elseif (tt0 .eq. 60)
      opts_bdgt@tmXBValues                  = (/ts0, ts0+12, ts0+24, ts0+36, ts0+48, tf0/)
     elseif (tt0 .eq. 72)
      opts_bdgt@tmXBValues                  = (/ts0, ts0+12, ts0+24, ts0+36, \
      					      	ts0+48, ts0+60, tf0/)
     elseif (tt0 .eq. 84)
      opts_bdgt@tmXBValues                  = (/ts0, ts0+12, ts0+24, ts0+36, ts0+48, \
      					        ts0+60, ts0+72, tf0/)
     elseif (tt0 .eq. 96)
      opts_bdgt@tmXBValues                  = (/ts0, ts0+12, ts0+24, ts0+36, ts0+48, \
      					        ts0+60, ts0+72, ts0+84, tf0/)
     end if

   ; Retrieve correct indices (from 'time_str') corresponding to the times above
     lab_ind = opts_bdgt@tmXBValues(:) - 4

   ; Create x-axis labels from values in 'time_str'
     opts_bdgt@tmXBLabels                  = time_str(lab_ind)

     opts_bdgt@trXMaxF			   = tf0	    ; X-axis maximum (12Z 8 Jul)
     opts_bdgt@trXMinF			   = ts0	    ; X-axis minimum (00Z 5 Jul)

     opts_bdgt@trYMaxF                     = 1000           ; Ignore BL
     opts_bdgt@trYMinF                     = 150            ; Include upper troposphere

     opts_bdgt@tiMainString                = ""
     opts_bdgt@gsnDraw                     = False   ; Do not draw the plot
     opts_bdgt@gsnFrame                    = False   ; Do not advance the frame
     opts_bdgt@lbLabelBarOn                = False   ; Keep label bar on

   ; Additional resources for black contours of circulation (LHS)
     opts_cont                             = True           ; New resources
     opts_cont@cnFillOn                    = False          ; Turn off filled contours
     opts_cont@cnLineColor                 = "black"        ; Line colour
     opts_cont@gsnDraw                     = False          ; Do not draw the plot
     opts_cont@gsnFrame                    = False          ; Do not advance the frame
     opts_cont@cnLevelSelectionMode        = "ExplicitLevels"

     if (r .eq. 0.76) then
      opts_cont@cnLevels                   = (/10.0, 14.0, 18.0, 22.0, 26.0/)
     elseif (r .eq. 0.48) then
      opts_cont@cnLevels                   = (/6.0, 8.0, 10.0, 14.0, 18.0/)
     elseif (r .eq. 0.24) then
      opts_cont@cnLevels                   = (/2.0, 4.0, 6.0, 8.0, 10.0/)
     end if

     opts_cont@cnLineThicknessF             = 1.5           ; Contour line thickness
     opts_cont@cnLineLabelsOn               = False         ; Turn off contour labels
     opts_cont@gsnContourNegLineDashPattern = 2             ; Negative contours dashed
     opts_cont@gsnContourPosLineDashPattern = 0             ; Positive contours solid
     opts_cont@cnInfoLabelOn                = False         ; Turn off info label
     circ_int@units                         = ""            ; Turn off string above plot

;=================
; Plot the data
;=================

     printMinMax(circ_int(:,:), False)
     printMinMax(tot_int(:,:), False)
     printMinMax(str_int(:,:), False)
     printMinMax(hadv_int(:,:), False)
     printMinMax(tilt_int(:,:), False)

     plot_circ0   = gsn_csm_pres_hgt(wks,circ_int(:,:),opts_bdgt)
     plot_tot0    = gsn_csm_pres_hgt(wks,tot_int(:,:),opts_bdgt)
     plot_str0    = gsn_csm_pres_hgt(wks,str_int(:,:),opts_bdgt)	
     plot_hadv0	  = gsn_csm_pres_hgt(wks,hadv_int(:,:),opts_bdgt)
     plot_tilt0   = gsn_csm_pres_hgt(wks,tilt_int(:,:),opts_bdgt)
;     plot_fric0   = gsn_csm_pres_hgt(wks,fric_int(:,:),opts_bdgt)

   ; First row
     panel(0)  = plot_circ0 
     panel(1)  = plot_tot0
     panel(2)  = plot_str0
     panel(3)  = plot_hadv0
     panel(4)  = plot_tilt0

;     panel(5)  = plot_fric0

     optsP                       = True
     optsP@gsnFrame              = False    ; Do not advance the frame
     optsP@gsnPanelLabelBar      = True     ; Turn on panel labelbar
     optsP@pmLabelBarWidthF      = 0.6
     optsP@pmLabelBarHeightF     = 0.15
     optsP@lbLabelFontHeightF    = 0.01
     optsP@lbPerimOn             = False
     optsP@lbLabelFont           = "Helvetica"
     optsP@lbBoxEndCapStyle      = "TriangleBothEnds"    ; Labelbar end shape

     if (fct .eq. 4) then
      optsP@lbTitleString = "Circulation budget tendency terms (10~S~-4~N~ s~S~-2~N~)"
     elseif (fct .eq. 5) then
      optsP@lbTitleString = "Circulation budget tendency terms (10~S~-5~N~ s~S~-2~N~)"
     end if 

     optsP@lbTitleFontHeightF    = 0.015
     optsP@lbTitleDirection      = "Across"
     optsP@lbTitlePosition       = "Bottom"
     optsP@txString              = ""
     optsP@gsnPanelFigureStrings = (/"a) Total LHS", \
                                     "b) Total RHS", \
                                     "c) Stretching", \
				     "d) Horizontal advection", \
                                     "e) Tilting"/)
     optsP@gsnMaximize           = True
     optsP@gsnPanelTop           = 0.98
     optsP@gsnPanelBottom        = 0.02
     optsP@gsnPaperOrientation	 = "landscape"
     optsP@amJust                = "TopLeft"
     optsP@gsnPanelFigureStringsFontHeightF = 0.0099 ; Label size (default 0.01)
     
     optsP@gsnPanelRowSpec       = True 
     gsn_panel(wks,panel,(/2,3/),optsP)              ; Draw as a single plot
     frame(wks)

end