; Script to calculate the radial momentum budget as in Montgomery et al. (2018)

; Run using:

; ncl opt=\"pdf\" dist=3.0 clr=5 calc=1 rad0=1 ave=1 nt=73 nr=41 ar=2.0
; wrt=1 grid=0 typ=\"stoa\" full=1 new0=0 mean0=0 cn0=\"slp\" tc_mbgt_tan_xz_comp.ncl

; 'opt'   = output file format ("pdf" or "x11")
; 'dist'  = size of box following storm (degrees)
; 'clr'	  = colour map for plots (option 4 is used in Roger Smith's papers):
; 	    "ncl_default" [1],  "bdr_extra" [2], "amwg256" [3], "GMT_polar" [4]
; 'calc'  = centre calculated on each level (0), calculated offline (1) or using 980 m vort (2)
; 'rad0'  = radius of final plots: 150 km (0) or 200 km (1)
; 'ave'	  = calculate 3/6/12-h average (1) rather than analyse individual time interval (0)
; 'wrt'	  = write out averages to text file (1) or leave alone (0 – default)
; 'grid'  = overlay NDC grid (1) or leave alone (0)
; 'typ'	  = VC phase ("sym","asym","stoa","atos")
; 'full'  = top 18 (1); top 8 (-5); top 6 (-6)
; 'ar'    = distance in degrees of outer radius (1.0, 2.0, etc)
; 'nr'    = number of radial circles between r = 0 and r = 'ar' (21, 41, 61, etc)
; 'nt'    = number of azimuth angles in cylindrical grid
; 'new0'  = reduced scale for plotting tangential wind and AAM (1)
; 'mean0' = plot eddy and mean terms (1) instead of diffusive tendencies (0)
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

begin

;==============================
; Prelash (reading in files)
;==============================

; For now, explicitly tell NCL which dates and simulations to analyse
  if (typ .eq. "stoa") then

   if (full .eq. 0) then             ; Top 5 simulations (previous choice of intervals)

    dat  = (/"02T12","02T12","02T12","02T12",\
             "03T00","03T00",\
             "03T12","03T12","03T12","03T12"/)

    ens0 = (/"em11", "em11", "em06", "em06",\
             "em09", "em09",\
             "em01", "em01", "em07", "em07"/)

    ts0  = (/58    , 73    , 51    , 69, \
             65    , 78    , \
             35    , 58    , 54    , 65/)

    tf0  = (/59    , 74    , 52    , 70, \
             66    , 79    , \
             36    , 59    , 55    , 66/)

   elseif (full .eq. -1) then	     ; Top 5 simulations (1-h intervals)

    dat  = (/"02T12","02T12","03T00","03T12","03T12"/)

    ens0 = (/"em11","em06","em09","em01","em07"/)

    ts0  = (/57    , 70   , 79   , 58   , 81/)

    tf0  = (/58    , 71   , 80   , 59   , 82/)

   elseif (full .eq. -5) then        ; Top 8 simulations (1-h intervals)

    dat  = (/"02T12","02T12","03T00","03T12","03T12","03T00","04T00","04T00"/)

    ens0 = (/"em11","em06","em09","em01","em07","em08","em01","em03"/)

    ts0  = (/57, 70, 79, 58, 81, 74, 58, 57/)

    tf0  = (/58, 71, 80, 59, 82, 75, 59, 58/)

   elseif (full .eq. -6) then        ; Top 6 simulations (1-h intervals)

    dat  = (/"02T12","02T12","03T00","03T12","03T12","03T00"/)

    ens0 = (/"em11","em06","em09","em01","em07","em08"/)

    ts0  = (/57, 70, 79, 58, 81, 74/)

    tf0  = (/58, 71, 80, 59, 82, 75/)

   elseif (full .eq. -2) then        ; Top 5 simulations (2-h intervals)

    dat  = (/"02T12","02T12","02T12","02T12","03T00","03T00","03T12","03T12","03T12","03T12"/)

    ens0 = (/"em11","em11","em06","em06","em09","em09","em01","em01","em07","em07"/)

    ts0  = (/57, 58, 69, 70, 79, 64, 58, 59, 81, 65/)

    tf0  = (/58, 59, 70, 71, 80, 65, 59, 60, 82, 66/)

   elseif (full .eq. -3) then	     ; All 11 'best-defined' VC simulations 

    dat =(/"02T12","02T12","02T12","03T00","03T00","03T12","03T12","03T12","04T00","04T00","04T00"/)

    ens0 = (/"em11","em06","em00","em09","em08","em00","em01","em07","em01","em02","em03"/)

    ts0  = (/57, 69, 65, 79, 73, 46, 58, 81, 51, 64, 57/)

    tf0  = (/58, 70, 66, 80, 74, 47, 59, 82, 52, 65, 58/)

   elseif (full .eq. -4) then	     ; Same 'best-defined' VC sims. (minor changes vs full = -3)

    dat =(/"02T12","02T12","02T12","03T00","03T00","03T12","03T12","03T12","04T00","04T00","04T00"/)

    ens0 = (/"em11","em06","em00","em09","em08","em00","em01","em07","em01","em02","em03"/)

    ts0  = (/57, 69, 65, 79, 73, 46, 58, 81, 51, 64, 57/)

    tf0  = (/58, 70, 66, 80, 74, 47, 59, 82, 52, 65, 58/)

   end if

  elseif (typ .eq. "sym") then

   if (full .eq. 1) then			; Top 18 simulations (objective method)

    dat  = (/"02T12","02T12","02T12",\
    	     "03T00","03T00","03T00","03T00","03T00",\
    	     "03T12","03T12","03T12","03T12",\
	     "04T00","04T00","04T00","04T00","04T00","04T00"/)

    ens0 = (/"em00","em06","em11",\
    	     "em02","em06","em08","em09","em10",\
    	     "em00","em01","em04","em07",\
	     "em01","em02","em03","em04","em06","em11"/)

    ts0  = (/74, 58, 50, \					; em00, em06, em11
    	     76, 69, 84, 74, 84, \				; em02, em06, em08, em09, em10
	     74, 70, 56, 58, \	 				; em00, em01, em04, em07
	     64, 79, 84, 33, 77, 45/)	     			; em01, em02, em03, em04, em06, em11

    tf0  = (/75, 59, 51, \   	 				; em00,	em06, em11
    	     77, 70, 85, 75, 85, \				; em02,	em06, em08, em09, em10
	     75, 71, 57, 59, \	 				; em00,	em01, em04, em07
	     65, 80, 85, 34, 78, 46/)	     	     		; em01,	em02, em03, em04, em06,	em11

   elseif (full .eq. 0) then 	     	        ; Top 5 simulations (previous choice of intervals)

    dat  = (/"02T12","02T12","02T12","02T12",\
             "03T00","03T00",\
             "03T12","03T12","03T12","03T12"/)

    ens0 = (/"em11", "em11", "em06", "em06",\
             "em09", "em09",\
             "em01", "em01", "em07", "em07"/)

    ts0  = (/56    , 71    , 49    , 67, \
             63    , 76    , \
             33    , 57    , 52    , 64/)

    tf0  = (/57    , 72    , 50    , 68, \
             64    , 77    , \
             34    , 58    , 53    , 65/)

   elseif (full .eq. -1) then          ; Top 5 simulations (1-h intervals)

    dat  = (/"02T12","02T12","03T00","03T12","03T12"/)

    ens0 = (/"em11","em06","em09","em01","em07"/)

    ts0  = (/55    , 48   , 77   , 54   , 72/)

    tf0  = (/56    , 49   , 78   , 55   , 73/)

   elseif (full .eq. -5) then        ; Top 8 simulations (1-h intervals)

    dat  = (/"02T12","02T12","03T00","03T12","03T12","03T00","04T00","04T00"/)

    ens0 = (/"em11","em06","em09","em01","em07","em08","em01","em03"/)

    ts0  = (/55, 48, 77, 54, 72, 71, 47, 51/)

    tf0  = (/56, 49, 78, 55, 73, 72, 48, 52/)

   elseif (full .eq. -6) then        ; Top 6 simulations (1-h intervals)

    dat  = (/"02T12","02T12","03T00","03T12","03T12","03T00"/)

    ens0 = (/"em11","em06","em09","em01","em07","em08"/)

    ts0  = (/55, 48, 77, 54, 72, 71/)

    tf0  = (/56, 49, 78, 55, 73, 72/)

   elseif (full .eq. -2) then	       ; Top 5 simulations (2-h intervals)

    dat  = (/"02T12","02T12","02T12","02T12","03T00","03T00","03T12","03T12","03T12","03T12"/)

    ens0 = (/"em11","em11","em06","em06","em09","em09","em01","em01","em07","em07"/)

    ts0  = (/55, 54, 48, 49, 77, 63, 54, 76, 72, 73/)

    tf0  = (/56, 55, 49, 50, 78, 64, 55, 77, 73, 74/)

   elseif (full .eq. -3) then	       ; All 11 'best-defined' VC simulations

    dat =(/"02T12","02T12","02T12","03T00","03T00","03T12","03T12","03T12","04T00","04T00","04T00"/)

    ens0 = (/"em11","em06","em00","em09","em08","em00","em01","em07","em01","em02","em03"/)

    ts0  = (/55, 48, 63, 77, 67, 55, 54, 72, 48, 40, 51/)

    tf0  = (/56, 49, 64, 78, 68, 56, 55, 73, 49, 41, 52/)

   elseif (full .eq. -4) then	       ; Same 'best-defined' VC sims. (minor changes vs full = -3)

    dat =(/"02T12","02T12","02T12","03T00","03T00","03T12","03T12","03T12","04T00","04T00","04T00"/)

    ens0 = (/"em11","em06","em00","em09","em08","em00","em01","em07","em01","em02","em03"/)

    ts0  = (/55, 48, 63, 77, 67, 78, 54, 72, 49, 46, 51/)

    tf0  = (/56, 49, 64, 78, 68, 79, 55, 73, 50, 47, 52/)

   elseif (full .eq. 3) then

    dat  = (/"02T12","03T00","03T12"/)

    ens0 = (/"em11","em09","em01"/)

    ts0  = (/56, 63, 33/)

    tf0  = (/57, 64, 34/)

   end if

  elseif (typ .eq. "asym") then

   if (full .eq. 1) then		; Top 18 simulations (objective	method)

    dat  = (/"02T12","02T12","02T12",\
             "03T00","03T00","03T00","03T00","03T00",\
             "03T12","03T12","03T12","03T12",\
             "04T00","04T00","04T00","04T00","04T00","04T00"/)

    ens0 = (/"em00","em06","em11",\
             "em02","em06","em08","em09","em10",\
             "em00","em01","em04","em07",\
             "em01","em02","em03","em04","em06","em11"/)

    ts0  = (/80, 49, 57, \				 ; em00, em06, em11
    	     61, 77, 89, 64, 77, \			 ; em02, em06, em08, em09, em10
             78, 83, 63, 64, \	 			 ; em00, em01, em04, em07
             71, 65, 73, 38, 83, 40/)			 ; em01, em02, em03, em04, em06, em11

    tf0  = (/81, 50, 58, \   	 			 ; em00, em06, em11
    	     62, 78, 90, 65, 78, \			 ; em02, em06, em08, em09, em10
             79, 84, 64, 65, \	 			 ; em00, em01, em04, em07
             72, 66, 74, 39, 84, 41/)			 ; em01, em02, em03, em04, em06, em11

   elseif (full .eq. 0) then

    dat  = (/"02T12","02T12","02T12","02T12",\
             "03T00","03T00",\
             "03T12","03T12","03T12","03T12"/)

    ens0 = (/"em11", "em11", "em06", "em06",\
             "em09", "em09",\
             "em01", "em01", "em07", "em07"/)

    ts0  = (/63    , 77    , 56    , 71, \
             70    , 80    , \
             36    , 61    , 55    , 66/)

    tf0  = (/64    , 78    , 57    , 72, \
             71    , 81    , \
             37    , 62    , 56    , 67/)

   elseif (full .eq. -1) then          ; Top 5 simulations (1-h intervals)

    dat  = (/"02T12","02T12","03T00","03T12","03T12"/)

    ens0 = (/"em11","em06","em09","em01","em07"/)

    ts0  = (/63    , 73   , 69   , 62   , 57/)

    tf0  = (/64    , 74   , 70   , 63   , 58/)

   elseif (full .eq. -5) then        ; Top 8 simulations (1-h intervals)

    dat  = (/"02T12","02T12","03T00","03T12","03T12","03T00","04T00","04T00"/)

    ens0 = (/"em11","em06","em09","em01","em07","em08","em01","em03"/)

    ts0  = (/63, 73, 69, 62, 57, 77, 61, 61/)

    tf0  = (/64, 74, 70, 63, 58, 78, 62, 62/)

   elseif (full .eq. -6) then        ; Top 6 simulations (1-h intervals)

    dat  = (/"02T12","02T12","03T00","03T12","03T12","03T00"/)

    ens0 = (/"em11","em06","em09","em01","em07","em08"/)

    ts0  = (/63, 73, 69, 62, 57, 77/)

    tf0  = (/64, 74, 70, 63, 58, 78/)

   elseif (full .eq. -2) then	       ; Top 5 simulations (2-h intervals)

    dat  = (/"02T12","02T12","02T12","02T12","03T00","03T00","03T12","03T12","03T12","03T12"/)

    ens0 = (/"em11","em11","em06","em06","em09","em09","em01","em01","em07","em07"/)

    ts0  = (/63, 64, 72, 73, 68, 69, 61, 62, 56, 57/)

    tf0  = (/64, 65, 73, 74, 69, 70, 62, 63, 57, 58/)

   elseif (full .eq. -3) then	       ; All 11 'best-defined' VC simulations

    dat =(/"02T12","02T12","02T12","03T00","03T00","03T12","03T12","03T12","04T00","04T00","04T00"/)

    ens0 = (/"em11","em06","em00","em09","em08","em00","em01","em07","em01","em02","em03"/)

    ts0  = (/63, 73, 67, 69, 76, 50, 62, 57, 60, 49, 59/)

    tf0  = (/64, 74, 68, 70, 77, 51, 63, 58, 61, 50, 60/)

   elseif (full .eq. -4) then	       ; Same 'best-defined' VC sims. (minor changes vs full = -3)

    dat =(/"02T12","02T12","02T12","03T00","03T00","03T12","03T12","03T12","04T00","04T00","04T00"/)

    ens0 = (/"em11","em06","em00","em09","em08","em00","em01","em07","em01","em02","em03"/)

    ts0  = (/63, 73, 67, 69, 76, 50, 62, 57, 60, 49, 59/)

    tf0  = (/64, 74, 68, 70, 77, 51, 63, 58, 61, 50, 60/)

   elseif (full .eq. 3) then

    dat  = (/"02T12","03T00","03T12"/)

    ens0 = (/"em11","em09","em01"/)

    ts0  = (/63, 70, 36/)

    tf0  = (/64, 71, 37/)

   end if

  elseif (typ .eq. "atos") then

   if (full .eq. 1) then 

   ;;; EDIT WHEN THESE ARE DEFINED ;;;

   elseif (full .eq. 0) then

    dat  = (/"02T12","02T12","03T00","03T12","03T12"/)

    ens0 = (/"em11" ,"em06" ,"em09" ,"em01" ,"em07"/)

    ts0  = (/66     , 66    , 72    , 40    , 58   /)

    tf0  = (/67     , 67    , 73    , 41    , 59   /)

   elseif (full .eq. -1) then          ; Top 5 simulations (1-h intervals)

    dat  = (/"02T12","02T12","03T00","03T12","03T12"/)

    ens0 = (/"em11","em06","em09","em01","em07"/)

    ts0  = (/69    , 66   , 75   , 74   , 60/)

    tf0  = (/70    , 67   , 76   , 75   , 61/)

   elseif (full .eq. -5) then        ; Top 8 simulations (1-h intervals)

    dat  = (/"02T12","02T12","03T00","03T12","03T12","03T00","04T00","04T00"/)

    ens0 = (/"em11","em06","em09","em01","em07","em08","em01","em03"/)

    ts0  = (/69, 66, 75, 74, 60, 87, 43, 69/)

    tf0  = (/70, 67, 76, 75, 61, 88, 44, 70/)

   elseif (full .eq. -6) then        ; Top 6 simulations (1-h intervals)

    dat  = (/"02T12","02T12","03T00","03T12","03T12","03T00"/)

    ens0 = (/"em11","em06","em09","em01","em07","em08"/)

    ts0  = (/69, 66, 75, 74, 60, 87/)

    tf0  = (/70, 67, 76, 75, 61, 88/)

   elseif (full .eq. -2) then          ; Top 5 simulations (2-h intervals)

    dat  = (/"02T12","02T12","02T12","02T12","03T00","03T00","03T12","03T12","03T12","03T12"/)

    ens0 = (/"em11","em11","em06","em06","em09","em09","em01","em01","em07","em07"/)

    ts0  = (/69, 70, 66, 78, 75, 76, 74, 75, 60, 61/)

    tf0  = (/70, 71, 67, 79, 76, 77, 75, 76, 61, 62/)

   elseif (full .eq. -3) then	       ; All 11 'best-defined' VC simulations

    dat =(/"02T12","02T12","02T12","03T00","03T00","03T12","03T12","03T12","04T00","04T00","04T00"/)

    ens0 = (/"em11","em06","em00","em09","em08","em00","em01","em07","em01","em02","em03"/)

    ts0  = (/69, 66, 76, 75, 87, 76, 74, 60, 43, 38, 73/)

    tf0  = (/70, 67, 77, 76, 88, 77, 75, 61, 44, 39, 74/)

   elseif (full .eq. -4) then	       ; Same 'best-defined' VC sims. (minor changes vs full = -3)

    dat =(/"02T12","02T12","02T12","03T00","03T00","03T12","03T12","03T12","04T00","04T00","04T00"/)

    ens0 = (/"em11","em06","em00","em09","em08","em00","em01","em07","em01","em02","em03"/)

    ts0  = (/69, 66, 77, 75, 87, 76, 74, 60, 43, 38, 72/)

    tf0  = (/70, 67, 78, 76, 88, 77, 75, 61, 44, 39, 73/)

   elseif (full .eq. 3) then

    dat  = (/"02T12","03T00","03T12"/)

    ens0 = (/"em11" ,"em09" ,"em01"/)

    ts0  = (/66     , 72    , 58   /)

    tf0  = (/67     , 73    , 59   /)

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

 ; Dimensions of 4D array
   nlev	       = 59			 ; Number of model levels
   ntimes      = dimsizes(dat)	         ; Number of time averages
   ndiag       = 11			 ; Number of diagnostics (MBGT terms)
   nrad	       = 41			 ; Number of radial circles for x-z plots 

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
  numINPUT_c  = numINFO_c[0]		; Number of times in 'pc' stream
  time_c      = numINFO_c[1]		; Array of times from 'pc' stream file
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

  storm_rel = st_rm(dat(st), ens0(st), "$ar/text/", 3.0, 0)
  u_cyc     = storm_rel[0]		   ; Zonal wind 
  v_cyc     = storm_rel[1]		   ; Meridional wind
  vel_cyc   = storm_rel[2]		   ; Vector wind 
  numTIMES  = storm_rel[3]                 ; Number of times in file
  lat_arr   = storm_rel[4]                 ; Storm-following domain subset (lat)
  lon_arr   = storm_rel[5]                 ; Storm-following domain subset (lon)
  centre    = storm_rel[6]		   ; Storm centre at each hour (lat/lon)

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

  thetas     = new(nt,float)
  pii        = 3.14159265
  radii      = fspan(0.0,ar,nr)      ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,nt)   ; Array of azimuth angles (0, 5, 10, etc)
  thetas_rad = (pii/180)*thetas      ; Degrees to radians
  rad_size   = dimsizes(radii)-1     ; Index for later calculation
  num_l      = 59      	   	     ; Number of model levels (edit to evolve smoothly)

; Define constants
  rd	= 287.0	  ; Gas constant for dry air (J/kg/K^2)
  cp0	= 1004.0  ; Specific heat of dry air at constant pressure (J/kg/K)
  kp	= 0.286   ; For Exner function calculations (0.286)
  pref	= 1000.0  ; Reference pressure p0 (1000 hPa)

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

    numTIMES    = 120                       ; Total number of times in all files [36]
    ct          = ts0(st)                   ; Counter variable
    llbox	= toint(dist*50)	    ; Size of lat/lon arrays (see line below)
    dsize	= (/llbox,llbox/)           ; Size of lat/lon arrays (depends on 'dist')
    d0		= dsize(0)
    d1		= dsize(1)

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

    c  = addfile(fili_c,"r")         	 ; Read in 'pc' stream [u,v,vort]
    d  = addfile(fili_d,"r")         	 ; Read in 'pd' stream [w,p,t,theta,kmh,tau]

;==================================
; Get the variables we will need
;==================================

   do it = ts0(st), tf0(st)		 ; Loop over times in file

    print("Working on time: "+time_str(it))
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

  ; Add metadata
    p@units = "hPa"
    z@units = "dam"

    dx0      = lonD(1) - lonD(0) ; Grid spacing (longitude)
    dy0      = latD(1) - latD(0) ; Grid spacing (latitude)

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
     vort_plane	   	    = vort(ilev,:,:)	   ; Relative vorticity
     avo_plane		    = avo(ilev,:,:)	   ; Absolute vorticity
     kmh_plane		    = kmh(ilev,:,:)	   ; Horizontal diffusivity (K_mh)

   ; Calculate location of vorticity maximum (storm centre) on model level
     vort_max      	    = max(vort_plane)
     vort_smth     	    = smth9_Wrap(vort_plane, 0.5, 0.5, True)
     vort_max_smth 	    = max(vort_smth)

   ; Find grid index of relative vorticity maximum
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
     dlat 	= abs(lat_max - centre(cn,it-1,0))
     dlon	= abs(lon_max - centre(cn,it-1,1))

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

     end if 

   elseif (calc .eq. 1) then

      lat_max    = centre(cn,it-1,0)               ; Does this relate to 'pc' or 'pd' stream???
      lon_max    = centre(cn,it-1,1)

   elseif (calc .eq. 2) then 

      lat_max    = t0 + (dy0 * inds0(0,0))
      lon_max    = n0 + (dx0 * inds0(0,1))
      print("Centre (980 m): "+lat_max+" degrees N, "+lon_max+" degrees E")

   end if

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
   ; vort_plane  -- latC,  lonC
   ; avo_plane   -- latC,  lonC
   ; f_plane     -- latC,  lonC

   ; w_plane     -- latD,  lon1D
   ; rho_plane	 -- latD,  lon1D
   ; kmh_plane	 -- latD,  lon1D

   ; difv_xplane -- latD,  lonD
   ; difv_yplane -- lat1D, lon1D 

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
      difv_x(:,iang)   = linint2_points(lonD_sub,latD_sub,difv_xplane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      difv_y(:,iang)   = linint2_points(lon1D_sub,lat1D_sub,difv_yplane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      w_int(:,iang)    = linint2_points(lon1D_sub,latD_sub,w_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      rho_int(:,iang)  = linint2_points(lon1D_sub,latD_sub,rho_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      kmh_int(:,iang)  = linint2_points(lon1D_sub,latD_sub,kmh_plane,False,\
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
      delete([/latC_plane,lonC_plane,latD_plane,lonD_plane,lon1D_plane,lat1D_plane/])
      delete([/latC_sub,lonC_sub,latD_sub,lonD_sub,lat1D_sub,lon1D_sub/])

     end do			    ; End model level loop (do ilev = 0, num_l-1)

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
      v_r0(:,:,:,:)		  = v_tan0(:,:,:,:) / r1

    ; Calculate [d(v/r)/dr]
      dvdr0(:,:,:,:)		  = v_r0(:,:,:,:) / dr

    ; Calculate [r]
      r_0(:,:,:,:)		  = r1

    ; Calculate [d/dr]
      ddr0(:,:,:,:)		  = 1 / dr

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
      mean_aam1	 = mean_aam0(:,:,0,:)

;================================================
; Calculate vertical derivatives of mean terms 
;================================================

    ; Redefine model level array for vertical derivative calculations
      hybC   = c->hybrid_ht(:)         ; 63 model ('theta') levels
      mlevs  = dimsizes(hybC)	       ; Size of model level array
      hybC@units = "m"		       ; Model level units 

      print_clock("Calculating vertical derivatives of mean terms...")
      
      do it = ts0(st), tf0(st)

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

      print_clock("Calculating vertical derivatives of eddy terms...")

      do it = ts0(st), tf0(st)

        do ilev = 0, num_l-1

        mp1 = min((/ilev+1,num_l-1/))
        mm1 = max((/ilev-1,0/))

          do irad = 0, dimsizes(radii)-1

           do th0 = 0, dimsizes(thetas)-1

           ;==========================================================
           ; Calculate vertical derivative of eddy tangential wind
           ;==========================================================

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

       rhs2_st = -( mean_vvel1 * mean_dvdz0 )

     ;=======================================
     ; Eddy radial vorticity flux: RHS (3)
     ;=======================================

       ptb_vez1 = -(ptb_vrad0 * ptb_vort0)

     ; Calculate the azimuthal average of 'ptb_vez1' and input into 'rhs3_st'
       rhs3_st  = dim_avg_n(ptb_vez1(:,:,:,:), 2)

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
       tau_ry1   = dim_avg_n( ( k_mh0 * ( ( (1 / r_0) * dudy0) + dvdr0 ) ), 2)

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

         dt0 = tf0(st) - ts0(st)

	 lhs1_av0(:,:) = (mean_vtan1(:,:,tf0(st)) - mean_vtan1(:,:,ts0(st))) / dt0

       elseif (ave .eq. 0) then 

         do it = ts0(st), tf0(st)

           cp1 = min((/it+1,tf0(st)/))
  	   cm1 = max((/it-1,ts0(st)/))
 
           lhs1A(:,:,it) = (mean_vtan1(:,:,cp1) - mean_vtan1(:,:,cm1)) / 2

         end do	           ; End loop over times (do it = ts0(st), tf0(st))

       end if 		   ; End IF statement 

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

     ; Add metadata for all the other important arrays (instantaneous)
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

        rhs1_av0   = dim_avg_n_Wrap(rhs1A(:,:,ts0(st):tf0(st)),2)
        rhs2_av0   = dim_avg_n_Wrap(rhs2A(:,:,ts0(st):tf0(st)),2)
        rhs3_av0   = dim_avg_n_Wrap(rhs3A(:,:,ts0(st):tf0(st)),2)
        rhs4_av0   = dim_avg_n_Wrap(rhs4A(:,:,ts0(st):tf0(st)),2)
        rhs5a_av0  = dim_avg_n_Wrap(rhs5aA(:,:,ts0(st):tf0(st)),2)
        rhs5b_av0  = dim_avg_n_Wrap(rhs5bA(:,:,ts0(st):tf0(st)),2)
        rhsT_av0   = dim_avg_n_Wrap(rhsTA(:,:,ts0(st):tf0(st)),2)

       end if 

;==============================================
; Switch order of dimensions before plotting
;==============================================

     ; 'gsn_csm_contour' plots the leftmost (rightmost) dimension on the y-axis (x-axis)
     ; Switch 'rad' with 'time' so that 'rad' is the rightmost dimension

     ; Mean radial wind -- also add metadata
       mean_vrad       = (/mean_vrad1(lev|:,time|:,rad|:)/)
       mean_vrad!0     = "lev"
       mean_vrad&lev   = hyb_plot
       mean_vrad!1     = "time"
       mean_vrad!2     = "rad"
       mean_vrad&rad   = radii
       mean_vrad@units = "m s~S~-1~N~"

     ; Mean tangential wind -- also add metadata
       mean_vtan       = (/mean_vtan1(lev|:,time|:,rad|:)/)
       mean_vtan!0     = "lev"
       mean_vtan&lev   = hyb_plot
       mean_vtan!1     = "time"
       mean_vtan!2     = "rad"
       mean_vtan&rad   = radii
       mean_vtan@units = "m s~S~-1~N~"

     ; Mean vertical velocity -- also add metadata
       mean_vvel       = (/mean_vvel1(lev|:,time|:,rad|:)/)
       mean_vvel!0     = "lev"
       mean_vvel&lev   = hyb_plot
       mean_vvel!1     = "time"
       mean_vvel!2     = "rad"
       mean_vvel&rad   = radii
       mean_vvel@units = "m s~S~-1~N~"

     ; Mean absolute angular momentum
       mean_aam        = (/mean_aam1(lev|:,time|:,rad|:)/)
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

        vrad_av0   = dim_avg_n_Wrap(mean_vrad(:,ts0(st):tf0(st),:),1)
        vtan_av0   = dim_avg_n_Wrap(mean_vtan(:,ts0(st):tf0(st),:),1)
        vvel_av0   = dim_avg_n_Wrap(mean_vvel(:,ts0(st):tf0(st),:),1)
        aam_av0    = dim_avg_n_Wrap(mean_aam(:,ts0(st):tf0(st),:),1) 

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

     ; Calculate total RHS tendency
       rhsT       = rhs1 + rhs2 + rhs3 + rhs4 + rhs5a + rhs5b
       rhsT!0     = "lev"
       rhsT&lev   = hyb_plot
       rhsT!1     = "time"
       rhsT!2     = "rad"
       rhsT&rad   = radii
       rhsT@units = "m s~S~-1~N~ h~S~-1~N~"

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

       if (ave .eq. 1) then 

      ; Total RHS 
	rhsT_av	      = rhs1_av + rhs2_av + rhs3_av + rhs4_av + rhs5a_av + rhs5b_av
        rhsT_av!0     = "lev"
        rhsT_av&lev   = hyb_plot
        rhsT_av!1     = "rad"
        rhsT_av&rad   = radii
        rhsT_av@units = "m s~S~-1~N~ h~S~-1~N~"

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

       end if 

;===================================================
; Fill 'plot_arr' with values from time averages
;===================================================

      if (wrt .eq. 1) then 

      ; Title resources for 'write_matrix' procedure
        opts_mat       = True 
	opts_mat@row   = True
	opts_mat@title = ""
	opt_arr        = (/"lhs1", "rhs1", "rhs2", \
		           "rhs3", "rhs4", "rhs5a", "rhs5b", \	   
			   "vrad", "vtan", "vvel", "aam"/)
        mbgt_arr       = (/lhs1_av, rhs1_av, rhs2_av, \
		       	   rhs3_av, rhs4_av, rhs5a_av, rhs5b_av, \
			   vrad_av0, vtan_av0, vvel_av0, aam_av0/)

      ; Fill composite array with values from terms in 'mbgt_arr'
        do mb = 0, ndiag-1 
         plot_arr(:,st,mb,:)  = mbgt_arr(mb,:,:)
	end do 

      end if 

     ; Tidy up before next iteration
       delete([/utc_date,year,month,day,hour,minute,second,date_str,time_str,time_arr,title_arr/])
       delete([/numINFO_c,numINFO_d,numINPUT_c,time_c,numINPUT_d,time_d/])
       delete([/storm_rel,u_cyc,v_cyc,vel_cyc,numTIMES,lat_arr,lon_arr,centre/])

     end do		; End loop over sets of dates (do st = 0, dimsizes(dat)-1)

     print_clock("Finished looping over input times! Now calculating averages...")

;=========================================================
; Before plotting, average over all chosen time periods
;=========================================================

    lhs1_c0  = plot_arr(:,:,0,:)

    rhs1_c0  = plot_arr(:,:,1,:)
    rhs2_c0  = plot_arr(:,:,2,:)
    rhs3_c0  = plot_arr(:,:,3,:)
    rhs4_c0  = plot_arr(:,:,4,:)
    rhs5a_c0 = plot_arr(:,:,5,:)
    rhs5b_c0 = plot_arr(:,:,6,:)

    vrad_c0  = plot_arr(:,:,7,:)
    vtan_c0  = plot_arr(:,:,8,:)
    vvel_c0  = plot_arr(:,:,9,:)
    aam_c0   = plot_arr(:,:,10,:)
 
    lhs1_c0@name  = "Local tendency of mean tangential wind"

    rhs1_c0@name  = "Mean radial vorticity flux"
    rhs2_c0@name  = "Mean vertical advection of mean tangential momentum"
    rhs3_c0@name  = "Eddy radial vorticity flux"
    rhs4_c0@name  = "Eddy vertical advection of eddy tangential momentum"
    rhs5a_c0@name = "Horizontal diffusive tendency of tangential momentum"
    rhs5b_c0@name = "Vertical diffusive tendency of tangential momentum"

    vrad_c0@name  = "Radial wind"
    vtan_c0@name  = "Tangential wind"
    vvel_c0@name  = "Vertical velocity"
    aam_c0@name   = "Absolute angular momentum"

;===============================================
; Calculate additional terms before averaging
;===============================================

  ; Total RHS
    rhsT_c0    = rhs1_c0 + rhs2_c0 + rhs3_c0 + rhs4_c0 + rhs5a_c0 + rhs5b_c0
    copy_VarCoords(rhs1_c0,rhsT_c0)
    copy_VarMeta(rhs1_c0,rhsT_c0)
    rhsT_c0@name = "Total RHS"

  ; All mean terms
    mean_c0    = rhs1_c0 + rhs2_c0 + rhs5a_c0 + rhs5b_c0
    copy_VarCoords(rhs1_c0,mean_c0)
    copy_VarMeta(rhs1_c0,mean_c0)
    mean_c0@name = "Sum of mean terms"

  ; All eddy terms
    eddy_c0    = rhs3_c0 + rhs4_c0
    copy_VarCoords(rhs2_c0,eddy_c0)
    copy_VarMeta(rhs2_c0,eddy_c0)
    eddy_c0@name = "Sum of eddy terms"

  ; Radial vorticity flux (eddy + mean)
    vflux_c0   = rhs1_c0 + rhs3_c0
    copy_VarCoords(rhs1_c0,vflux_c0)
    copy_VarMeta(rhs1_c0,vflux_c0)
    vflux_c0@name = "Total radial vorticity flux"

  ; Vertical advection (eddy + mean)
    vadv_c0   = rhs2_c0 + rhs4_c0
    copy_VarCoords(rhs1_c0,vadv_c0)
    copy_VarMeta(rhs1_c0,vadv_c0)
    vflux_c0@name = "Total vertical advection"
 
;==================================================================
; Average over time dimension to calculate composite diagnostics
;================================================================== 

    lhs1_c  = dim_avg_n_Wrap(lhs1_c0(:,:,:),1)

    rhs1_c  = dim_avg_n_Wrap(rhs1_c0(:,:,:),1)
    rhs2_c  = dim_avg_n_Wrap(rhs2_c0(:,:,:),1)
    rhs3_c  = dim_avg_n_Wrap(rhs3_c0(:,:,:),1)
    rhs4_c  = dim_avg_n_Wrap(rhs4_c0(:,:,:),1)
    rhs5a_c = dim_avg_n_Wrap(rhs5a_c0(:,:,:),1)
    rhs5b_c = dim_avg_n_Wrap(rhs5b_c0(:,:,:),1)

    rhsT_c  = dim_avg_n_Wrap(rhsT_c0(:,:,:),1)
    mean_c  = dim_avg_n_Wrap(mean_c0(:,:,:),1)
    eddy_c  = dim_avg_n_Wrap(eddy_c0(:,:,:),1)

    vflux_c = dim_avg_n_Wrap(vflux_c0(:,:,:),1)
    vadv_c  = dim_avg_n_Wrap(vadv_c0(:,:,:),1)

    vrad_c  = dim_avg_n_Wrap(vrad_c0(:,:,:),1)
    vtan_c  = dim_avg_n_Wrap(vtan_c0(:,:,:),1)
    vvel_c  = dim_avg_n_Wrap(vvel_c0(:,:,:),1)
    aam_c   = dim_avg_n_Wrap(aam_c0(:,:,:),1)
  
;===============================
; Add metadata before ploting 
;===============================

    lhs1_c!0    = "lev"
    lhs1_c&lev  = hyb_plot
    lhs1_c!1	= "rad"
    lhs1_c&rad	= radii

    rhs1_c!0    = "lev"
    rhs1_c&lev  = hyb_plot
    rhs1_c!1	= "rad"
    rhs1_c&rad	= radii

    rhs2_c!0    = "lev"
    rhs2_c&lev  = hyb_plot
    rhs2_c!1	= "rad"
    rhs2_c&rad	= radii

    rhs3_c!0    = "lev"
    rhs3_c&lev  = hyb_plot
    rhs3_c!1	= "rad"
    rhs3_c&rad	= radii

    rhs4_c!0    = "lev"
    rhs4_c&lev  = hyb_plot
    rhs4_c!1	= "rad"
    rhs4_c&rad	= radii

    rhs5a_c!0    = "lev"
    rhs5a_c&lev  = hyb_plot
    rhs5a_c!1	= "rad"
    rhs5a_c&rad	= radii

    rhs5b_c!0    = "lev"
    rhs5b_c&lev  = hyb_plot
    rhs5b_c!1	= "rad"
    rhs5b_c&rad	= radii

    rhsT_c!0    = "lev"
    rhsT_c&lev  = hyb_plot
    rhsT_c!1	= "rad"
    rhsT_c&rad	= radii

    mean_c!0    = "lev"
    mean_c&lev  = hyb_plot
    mean_c!1	= "rad"
    mean_c&rad	= radii

    eddy_c!0    = "lev"
    eddy_c&lev  = hyb_plot
    eddy_c!1	= "rad"
    eddy_c&rad	= radii

    vflux_c!0   = "lev"
    vflux_c&lev = hyb_plot
    vflux_c!1   = "rad"
    vflux_c&rad = radii

    vadv_c!0    = "lev"
    vadv_c&lev  = hyb_plot
    vadv_c!1    = "rad"
    vadv_c&rad  = radii

    vrad_c!0    = "lev"
    vrad_c&lev  = hyb_plot
    vrad_c!1	= "rad"
    vrad_c&rad	= radii

    vtan_c!0    = "lev"
    vtan_c&lev  = hyb_plot
    vtan_c!1	= "rad"
    vtan_c&rad	= radii

    vvel_c!0    = "lev"
    vvel_c&lev  = hyb_plot
    vvel_c!1	= "rad"
    vvel_c&rad	= radii

    aam_c!0    = "lev"
    aam_c&lev  = hyb_plot
    aam_c!1	= "rad"
    aam_c&rad	= radii

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
     opts_mbgt@cnFillColors               = (/3,5,7,-1,11,14,18/)
     opts_mbgt@cnLevels                   = (/-20,-5,0,0.5,5,20/)
    elseif (clr .eq. 5) then 
     opts_mbgt@cnFillPalette               = "BlueDarkRed18"
     opts_mbgt@cnLevels                    = (/-10.0, -5.0, -3.0, -2.0, -1.0, -0.5, -0.2, \
                                                0, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0/)
     opts_mbgt@cnFillColors                = (/0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17/)
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
;    opts_mbgt@tiXAxisFontHeightF           = 0.025   ; X-axis title font height 

    opts_mbgt@tmXTOn                       = "False" ; Turn off top x-axis TM
    opts_mbgt@tmYROn                       = "False" ; Turn off right y-axis TM
    opts_mbgt@gsnMaximize                  = False   ; Maximise plot size
    opts_mbgt@gsnAddCyclic                 = False
    opts_mbgt@tmXBMode                     = "Explicit"  ; Set tick marks explicitly
;    opts_mbgt@tmXBLabelFontHeightF	   = 0.02    ; Label font height 

    if (rad0 .eq. 0) then
     opts_mbgt@tmXBValues                   = (/0.0,0.25,0.5,0.75,1.0,1.25,1.5/)
     opts_mbgt@tmXBLabels                   = (/"0","25","50","75","100","125","150"/)
    else if (rad0 .eq. 1) then
     opts_mbgt@tmXBValues                   = (/0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0/)
     opts_mbgt@tmXBLabels                   = (/"0","25","50","75","100","125","150","175","200"/)
    end if
    end if

    opts_mbgt@tmYLMode			   = "Explicit"
    opts_mbgt@tmYLValues		   = fspan(2.0, 18.0, 9)
    opts_mbgt@tmYLLabels		   = (/"2.0","4.0","6.0","8.0","10.0",\
    					       "12.0","14.0","16.0","18.0"/)
;    opts_mbgt@tmYLLabelFontHeightF	   = 0.02
    opts_mbgt@gsnYAxisIrregular2Linear	   = True	; Linearise the y-axis 

    opts_mbgt@trYMinF			   = 0.0	; y-axis minimum value
    opts_mbgt@trYMaxF			   = 17.6	; y-axis maximum value 

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
                                               0.7, 0.8, 1.0, 1.2, 1.5, 2.0/)
    opts_vvel2@cnFillColors                = (/-1,0,1,2,3,4,5,6,7,8,9,10,11/)
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

    txid			= new(12,graphic)
    amid	    	        = new(12,graphic)

  ; Label text resources 
    txres                       = True
    txres@txPerimOn             = True
    txres@txFontHeightF         = 0.010
    txres@txBackgroundFillColor = "White"

  ; Label position resources 
    amres 	   	        = True
    amres@amParallelPosF	= 0.5			; Right edge (-0.5 for left edge)
    amres@amOrthogonalPosF	= -0.5			; Top edge
    amres@amJust		= "TopRight"

    if (mean0 .eq. 1) then 
     panel_strings              = (/"a) vrad", "b) vtan", "c) vvel", "d) aam", \
                                    "e) V_mz", "f) V_mv", "g) V_ez", \
                                    "h) V_ev", "i) Mean", "j) Eddy", \
                                    "k) All V tend", "l) V_t"/)
    elseif (mean0 .eq. 0) then 
     panel_strings              = (/"a) vrad", "b) vtan", "c) vvel", "d) aam", \
                                    "e) V_mz", "f) V_mv", "g) V_ez", \
                                    "h) V_ev", "i) V_dh", "j) V_dz", \
                                    "k) All V tend", "l) V_t"/)    
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

  ; Output file location and type
    output = "$sam/nepartak/images/cyl_coords/mbgt_tan_xz_comp_"+typ+"_full"+full+"_"+cn0
    wks = gsn_open_wks(opt,output)

  ; Radial wind (filled contours)
    panel_vrad  = gsn_csm_contour(wks,vrad_c(:,:),opts_vrad2) ; [u]
    panel_zero  = gsn_csm_contour(wks,vrad_c(:,:),opts_zero)  ; Zero line contour
    zero0       = ColorNegDashZeroPosContour(panel_zero,"transparent","black","transparent")
    overlay(panel_vrad, zero0)
    txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
    amid(0)     = gsn_add_annotation(panel_vrad, txid(0), amres)

  ; Tangential wind (filled contours)
    panel_vtan  = gsn_csm_contour(wks,vtan_c(:,:),opts_vtan2) ; [v]
    txid(1)     = gsn_create_text(wks, panel_strings(1), txres)
    amid(1)     = gsn_add_annotation(panel_vtan, txid(1), amres)

  ; Vertical velocity (filled contours)
    panel_vvel  = gsn_csm_contour(wks,vvel_c(:,:),opts_vvel2) ; [w]
    txid(2)     = gsn_create_text(wks, panel_strings(2), txres)
    amid(2)     = gsn_add_annotation(panel_vvel, txid(2), amres)

  ; Absolute angular momentum (filled contours)
    panel_aam   = gsn_csm_contour(wks,aam_c(:,:),opts_aam)   ; [AAM]
    txid(3)     = gsn_create_text(wks, panel_strings(3), txres)
    amid(3)     = gsn_add_annotation(panel_aam, txid(3), amres)

  ; Mean radial influx of absolute vertical vorticity [V_mzeta]
    opts_mbgt@vpXF       = 0.07
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF	 = 0.15
    panel0      = gsn_csm_contour(wks,rhs1_c(:,:),opts_mbgt) ; [V_mzeta]
    cont0       = gsn_csm_contour(wks,rhs1_c(:,:),opts_zero) ; Draw contours

    panel0_ovr  = gsn_csm_contour(wks,vrad_c(:,:),opts_vrad) ; [u]
    panel0_ovr2 = gsn_csm_contour(wks,vvel_c(:,:),opts_vvel) ; [w]

    overlay(panel0, cont0)				     ; Overlay contours
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
    panel1      = gsn_csm_contour(wks,rhs2_c(:,:),opts_mbgt) ; [V_mv]
    cont0       = gsn_csm_contour(wks,rhs2_c(:,:),opts_zero) ; Draw contours

    panel1_ovr  = gsn_csm_contour(wks,vrad_c(:,:),opts_vrad) ; [U]
    panel1_ovr2 = gsn_csm_contour(wks,vvel_c(:,:),opts_vvel) ; [w]

    overlay(panel1, cont0)                                   ; Overlay contours    
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
    panel2      = gsn_csm_contour(wks,rhs3_c(:,:),opts_mbgt) ; [V_ezeta]
    cont0       = gsn_csm_contour(wks,rhs3_c(:,:),opts_zero) ; Draw contours

    panel2_ovr  = gsn_csm_contour(wks,vrad_c(:,:),opts_vrad) ; [U]
    panel2_ovr2 = gsn_csm_contour(wks,vvel_c(:,:),opts_vvel) ; [w]

    overlay(panel2, cont0)                                   ; Overlay contours
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
    panel3      = gsn_csm_contour(wks,rhs4_c(:,:),opts_mbgt) ; [V_ev]
    cont0       = gsn_csm_contour(wks,rhs4_c(:,:),opts_zero) ; Draw contours

    panel3_ovr  = gsn_csm_contour(wks,vrad_c(:,:),opts_vrad) ; [U]
    panel3_ovr2 = gsn_csm_contour(wks,vvel_c(:,:),opts_vvel) ; [w]

    overlay(panel3, cont0)                                   ; Overlay contours
    overlay(panel3, panel3_ovr)
    overlay(panel3, panel3_ovr2)
    txid(7)     = gsn_create_text(wks, panel_strings(7), txres)
    amid(7)     = gsn_add_annotation(panel3, txid(7), amres)
    delete([/cont0/])

  ; Horizontal diffusive tendency of tangential momentum [V_dh]
    opts_mbgt@vpXF       = 0.07
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15

    if (mean0 .eq. 1) then 
     panel4     = gsn_csm_contour(wks,mean_c(:,:),opts_mbgt) ; [V_dh]
     cont0      = gsn_csm_contour(wks,mean_c(:,:),opts_zero) ; Draw contours
    elseif (mean0 .eq. 0) then  
     panel4     = gsn_csm_contour(wks,rhs5a_c(:,:),opts_mbgt) ; [V_dh]
     cont0      = gsn_csm_contour(wks,rhs5a_c(:,:),opts_zero) ; Draw contours
    end if 

    panel4_ovr  = gsn_csm_contour(wks,vrad_c(:,:),opts_vrad) ; [U]
    panel4_ovr2 = gsn_csm_contour(wks,vvel_c(:,:),opts_vvel) ; [w]

    overlay(panel4, cont0)                                   ; Overlay contours
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

    if (mean0 .eq. 1) then 
     panel5     = gsn_csm_contour(wks,eddy_c(:,:),opts_mbgt)  ; [V_dz]
     cont0      = gsn_csm_contour(wks,eddy_c(:,:),opts_zero)  ; Draw contours
    elseif (mean0 .eq. 0) then 
     panel5     = gsn_csm_contour(wks,rhs5b_c(:,:),opts_mbgt) ; [V_dz]
     cont0      = gsn_csm_contour(wks,rhs5b_c(:,:),opts_zero) ; Draw contours
    end if 

    panel5_ovr  = gsn_csm_contour(wks,vrad_c(:,:),opts_vrad)  ; [U]
    panel5_ovr2 = gsn_csm_contour(wks,vvel_c(:,:),opts_vvel)  ; [w]

    overlay(panel5, cont0)                                    ; Overlay contours
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
    panel6      = gsn_csm_contour(wks,rhsT_c(:,:),opts_mbgt) ; [All V tend]
    cont0       = gsn_csm_contour(wks,rhsT_c(:,:),opts_zero) ; Draw contours

    panel6_ovr  = gsn_csm_contour(wks,vrad_c(:,:),opts_vrad) ; [U]
    panel6_ovr2 = gsn_csm_contour(wks,vvel_c(:,:),opts_vvel) ; [w]

    overlay(panel6, cont0)                                   ; Overlay contours
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
    panel7      = gsn_csm_contour(wks,lhs1_c(:,:),opts_mbgt) ; [V_t]
    cont0       = gsn_csm_contour(wks,lhs1_c(:,:),opts_zero) ; Draw contours

    panel7_ovr  = gsn_csm_contour(wks,vrad_c(:,:),opts_vrad) ; [U]
    panel7_ovr2 = gsn_csm_contour(wks,vvel_c(:,:),opts_vvel) ; [w]

    overlay(panel7, cont0)                                   ; Overlay contours
    overlay(panel7, panel7_ovr)
    overlay(panel7, panel7_ovr2)
    txid(11)    = gsn_create_text(wks, panel_strings(11), txres)
    amid(11)    = gsn_add_annotation(panel7, txid(11), amres)
    delete([/cont0/])

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    pres			= True 
    maximize_output(wks,pres)				   ; Calls 'draw' and 'frame'

end 