; Script to compute the storm track for a TC using MSLP and relative vorticity
; Write out max 10-m windspeed and min. MSLP following a TC (storm tracker) 

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_rm.ncl"

; ncl dat=\"02T12\" int=1 opt=\"x11\" ens0=\"em11\" rad0=3.0 dist=3.0 vc=0 
; lay=0 p_cp=3 zoom=1 lgd0=\"tr\" tc_fc_vort_mslp.ncl

; 'dat'  = initialisation time string ("02T12", "03T00", ...)
; 'ens0' = ensemble simulation (em00-em11)
; 'int'	 = 
; 'opt'  = output file format ("pdf" or "x11")
; 'rad0' = 
; 'dist' = size of box following storm (degrees)
; 'vc'   = 
; 'lay'  = calculate layer average (1) or single level (0)
; 'p_cp' = pressure level for calculating storm centre
; 'lgd0' = legend position ("tr","tl","br","bl")

begin

;==============================================
; Find 'pa' stream data to read in and plot
;==============================================

  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0
  fili_a1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pa*.nc")
  fili_a      = diri+"/"+fili_a1

  numFILES    = dimsizes(fili_a)		; 10 input files (T+000 to T+108)
  numINPUT    = numFILES * 12 			; 120 input time intervals (12 per input file)

  print("Working on ensemble member "+ens0+" (initialised: "+dat+")")

  setvalues NhlGetWorkspaceObjectId
  "wsMaximumSize" : 10000000000
  end setvalues

; Create arrays to hold values of storm position, minimum MSLP and max 10-m windspeed
  slp_box   = new(numINPUT,float)  		  ; Minimum MSLP
  vec_box   = new(numINPUT,float)		  ; Maximum 10-m windspeed
  centre    = new((/numINPUT,2/),"float")	  ; Storm centre (lat/lon)
  	      					  ; [120 * 2]

; Set counter variable (for writing above diagnostics to text files)
  ct        = 0

; Create arrays to hold date strings 
  time_str  = new(numINPUT,string)
  time_arr  = new(numINPUT,string)
  title_arr = new(numINPUT,string)

; Create initialisation string
  if (dat .eq. "02T12") then 
   date_str = "12 UTC 2 July"
  elseif (dat .eq. "03T00") then 
   date_str = "00 UTC 3 July"
  elseif (dat .eq. "03T12") then 
   date_str = "12 UTC 3 July"
  elseif (dat .eq. "04T00") then 
   date_str = "00 UTC 4 July"
  end if  

;=================================================
; Read in best track information from text file
;=================================================

; Number of times in IBTrACS dataset for Nepartak
  dsn 	    = 28

; Read in all latitude/longitude/date information
  lat_bt0   = asciiread("$sam/um/cp/nepartak/lat_nep.txt",(/dsn/),"float")
  lon_bt0   = asciiread("$sam/um/cp/nepartak/lon_nep.txt",(/dsn/),"float")
  dat_bt0   = asciiread("$sam/um/cp/nepartak/dat_nep.txt",(/dsn/),"string")

; From these arrays, extract storm position at start of VC analysis period
; Use these data to constrain the initial position of the storm in the UM analysis
  if (dat .eq. "02T12") then
   lat0     = lat_bt0(1)
   lon0     = lon_bt0(1)
  elseif (dat .eq. "03T00") then
   lat0     = lat_bt0(3)
   lon0     = lon_bt0(3)
  elseif (dat .eq. "03T12") then
   lat0     = lat_bt0(5)
   lon0     = lon_bt0(5)
  elseif (dat .eq. "04T00") then
   lat0     = lat_bt0(7)
   lon0     = lon_bt0(7)
  end if

;==========================================================
; Get the variables we will need (loop over input files)
;==========================================================

  do nf     = 0, numFILES-1
  f  	    = addfile(fili_a(nf),"r")

;================================================
; Number of times in different variable arrays
;================================================

; Times in file 
  if (vc .eq. 1) then 
   time  = f->t_1
   times = 12
  elseif (vc .eq. 0) then	
   time  = f->t
   times = dimsizes(time)
  end if 

;=======================
; All other variables 
;=======================

  lon   = f->longitude    ; longitude (1098 points --> 109.04 to 152.92 degrees E)
  lat   = f->latitude     ; latitude  (810 points --> 1.8 to 34.16 degrees N)
  lon1  = f->longitude_1  ; longitude (1098 points --> 109.02 to 152.90 degrees E)
  lat1  = f->latitude_1   ; latitude  (811 points --> 1.78 to 34.18 degrees N)

  u_varname = "x-wind"
  v_varname = "y-wind"

  u     = f->$u_varname$  ; 10-m zonal velocity (m/s) [t_3 | z | lat_1 | lon_1]
  v     = f->$v_varname$  ; 10-m meridional velocity (m/s) [t_3 | z | lat_1 | lon_1]

; Read in mean sea level pressure (Pa) [t | z | lat | lon]
; For 1-h data (vc = 0), ignore first time interval (nf = 0)
  if (nf .eq. 0 .and. vc .eq. 0) then
   mslp  = f->p(1:times-1,:,:,:)
  else
   mslp  = f->p
  end if

  mslp  = mslp / 100      ; Convert to hPa
  mslp@units = "hPa"

  dx    = lon(1) - lon(0) ; Grid spacing (longitude)
  dy    = lat(1) - lat(0) ; Grid spacing (latitude)

;===========================
; Loop over times in file
;===========================

  if (nf .eq. 0 .and. vc .eq. 0) then		; For 1-h data, ignore first time interval
   tf = times-2
  else
   tf = times-1
  end if

  do it = 0, tf, int        ; Loop over times in each file

   ctt = ct + 1 
   print("Working on forecast hour: "+ctt)
   title = ctt+" h forecast, initialised at "+date_str

;=====================================================
; Create correct date strings for each output time
;=====================================================

; CAN RUN WITHOUT EDITING, BUT EDIT ONCE RUNNING SMOOTHLY

  month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
	         "Jul","Aug","Sep","Oct","Nov","Dec"/)

  if (nf .eq. 0 .and. vc .eq. 0) then
   utc_date = cd_calendar(time(1:times-1), 0)
  else 
   utc_date = cd_calendar(time, 0)
  end if

  year   = tointeger(utc_date(:,0))
  month  = tointeger(utc_date(:,1))
  day    = tointeger(utc_date(:,2))
  hour   = tointeger(utc_date(:,3))
  minute = tointeger(utc_date(:,4))
  second = utc_date(:,5)

; Correct for errors in the code (round up value of hour when minutes = 59)
  if (minute(it) .eq. 9) then
   minute(it) = 10
  elseif (minute(it) .eq. 19) then 
   minute(it) = 20
  elseif (minute(it) .eq. 29) then
   minute(it) = 30
  elseif (minute(it) .eq. 39) then
   minute(it) = 40
  elseif (minute(it) .eq. 54) then
   minute(it) = 55
  end if

  min_arr      = (/":00",":05",":10",":15",":20",":25",\
		   ":30",":35",":40",":45",":50",":55"/)

  hr_min       = hour(it)+min_arr(it)

  time_arr(ct) = sprinti("%0.2i", day(it)) + month_abbr(month(it)) + \
	       	         "_" + hr_min + "Z"
  time_str(ct) = hr_min + " UTC " + sprinti("%0.2i ", day(it)) + month_abbr(month(it))

;  print("Working on time: "+time_str(ct))
  title_arr(ct) = "Valid at "+time_str(ct)+" (T+"+ctt+")"

;===============================================
; Find cyclone centre on RA1T grid (use MSLP)
;===============================================

; Read in mean sea level pressure and 10-m windspeed data 
; At the initial time, focus on the region close to the best track location

  if (ct .eq. 0) then 

 ; Define small grid subset centred on early best track location 	
   latS	     = lat0-1.0
   latF	     = lat0+1.0
   lonS	     = lon0-1.0
   lonF	     = lon0+1.0

 ; Focus on smaller region, constrained by the best track location
   u_plane       = u(it,0,{latS:latF},{lonS:lonF})
   v_plane       = v(it,0,{latS:latF},{lonS:lonF})
   slp_plane     = mslp(it,0,{latS:latF},{lonS:lonF})

  else

   u_plane       = u(it,0,:,:)
   v_plane       = v(it,0,:,:)
   slp_plane     = mslp(it,0,:,:)

  end if 

; Also calculate 10-m vector windspeed
  speed_plane = sqrt(u_plane^2.0 + v_plane^2.0)
  copy_VarCoords(u_plane, speed_plane)
  copy_VarAtts(u_plane, speed_plane)
  copy_VarMeta(u_plane, speed_plane)
  speed_plane@units = "m s~S~-1~N~"

; MSLP
  slp_min       = min(slp_plane)
  slp_smth      = smth9_Wrap(slp_plane, 0.5, 0.5, True)
  slp_max_smth  = max(slp_smth)

  dims          = dimsizes(slp_plane)
  slp1d         = ndtooned(slp_plane)
  inds          = ind_resolve(minind(slp1d),dims)
  slp1d_sm      = ndtooned(slp_smth)
  inds_sm       = ind_resolve(maxind(slp1d_sm),dims)

  if (ct .eq. 0) then 
   lat_minslp   = latS + (dy * inds(0,0))
   lon_minslp   = lonS + (dx * inds(0,1))
  else 
   lat_minslp   = lat(0) + (dy * inds(0,0))
   lon_minslp   = lon(0) + (dx * inds(0,1))
  end if 

;========================================================
; EMPLOY A SAFETY NET IF WE IDENTIFY THE WRONG CYCLONE
;========================================================

  if (ct .eq. 0) then

   centre(ct,0) = lat_minslp
   centre(ct,1) = lon_minslp

  else

   lt0  = centre(ct-1,0)             ; Previous TC latitude
   ln0  = centre(ct-1,1)             ; Previous TC longitude
   d_lat = abs(lat_minslp-lt0)       ; Change in latitude (t1 - t0)
   d_lon = abs(lon_minslp-ln0)       ; Change in longitude

   if (d_lat .gt. 0.75 .or. d_lon .gt. 0.75) then

    print("Previous storm latitude ("+lt0+") ; longitude ("+ln0+")")
    print("Change in lat = "+d_lat+" deg")
    print("Change in lon = "+d_lon+" deg")

    print_clock("Now finding correct centre based upon previous storm position...")

  ; Create smaller grid [0.5 degrees]
    lt1 = lt0 - 0.5
    lt2 = lt0 + 0.5
    ln1 = ln0 - 0.5
    ln2 = ln0 + 0.5

    slp_plane0 = slp_plane({lt1:lt2},{ln1:ln2})	      ; Smaller grid (1.0 x 1.0 deg)
    slp_min0   = min(slp_plane0)                      ; Find minimum on smaller grid
    dims_p0    = dimsizes(slp_plane0)                 ; Size of smaller grid
    slp_1d0    = ndtooned(slp_plane0)                 ; Create 1-D array
    inds_p0    = ind_resolve(minind(slp_1d0),dims_p0) ; Find index of minimum

    lat_minslp = lt1 + (dy * inds_p0(0,0))
    lon_minslp = ln1 + (dy * inds_p0(0,1))
    print("New centre: "+lat_minslp+" deg N, "+lon_minslp+" deg E")
    centre(ct,0) = lat_minslp
    centre(ct,1) = lon_minslp

    delete([/slp_plane0, slp_min0, dims_p0, slp_1d0, inds_p0/])

   else

    print("Cyclone centre: "+lat_minslp+" deg N, "+lon_minslp+" deg E")
    centre(ct,0) = lat_minslp
    centre(ct,1) = lon_minslp

   end if

  end if

; Tidy up
  delete([/slp1d,inds,slp1d_sm,slp_smth,inds_sm/])

;=====================================================================================
; Once we've correctly identified the storm, write MSLP and 10-m windspeed to files
;=====================================================================================

  rad = 1.0

  slp_box(ct) = min(slp_plane({centre(ct,0)-rad:centre(ct,0)+rad},\
			      {centre(ct,1)-rad:centre(ct,1)+rad}))
  vec_box(ct) = max(speed_plane({centre(ct,0)-rad:centre(ct,0)+rad},\
			      {centre(ct,1)-rad:centre(ct,1)+rad}))

; Tidy up 
  delete([/slp_plane,u_plane,v_plane,speed_plane/])
  delete([/utc_date,year,month,day,hour,minute,second/])

  ct = ct + 1 ; Counter variable (time)

  end do     ; End time loop (do it = 0, times-1, int)

; Tidy up before next iteration 
  delete([/time,times,u,v,mslp/])

 end do	      ; End input file loop (do nf = 0, numFILES-1)

;=================================================
; Now find 'pb' stream data to read in and plot
;=================================================

; Then plot MSLP and wind/vorticity with both centres overlaid for comparison
; Also write out the new storm centre positions for comparison

  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0
  fili_b1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pb.nc")
  fili_b      = diri+"/"+fili_b1 

; Get number of times in file using built-in function ('nc_times')
  numINFO     = nc_times(fili_b)
  numINPUT_b  = numINFO[0]              ; Number of time intervals
  time_b      = numINFO[1]              ; Array containing timestamp info at each time
  print("Total number of time intervals = "+numINPUT_b)

; Create array to hold updated storm centre positions (relative vorticity and geopotential height)
  centre_new  = new((/3,numINPUT_b,2/),"float")         ; [3 * 117 * 2] --> same size as 'time_b'

; Remove any old arrays 
  delete([/lat0, lon0, lat1, lon1, lat, lon/])
  delete([/date_str, time_str, time_arr, title_arr/])

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
  do ct = 0, numINPUT_b-1
   if (minute(ct).gt.30) then
     hour(ct) = hour(ct)+1
   end if
  end do
  delete(ct)

; Create arrays to hold strings for output plots
  date_str  = new(numINPUT_b,string)
  time_str  = new(numINPUT_b,string)
  time_arr  = new(numINPUT_b,string)
  title_arr = new(numINPUT_b,string)

  date_str  = sprinti("%0.2iUTC ", hour) + \
              sprinti("%0.2i ", day) \
            + month_abbr(month)
  time_str  = sprinti("%0.2i UTC ", hour) + \
              sprinti("%0.2i ", day) \
            + month_abbr(month)
  time_arr  = sprinti("%0.2i", day) + \
              month_abbr(month) + \
              "_" + sprinti("%0.2iZ", hour)

;========================================================================
; Mimic the features of 'st_rm' without calling the external function
;========================================================================

  lat_arr      = new((/numINPUT,2/),float)
  lon_arr      = new((/numINPUT,2/),float)

; Create domain subset arrays following the storm
  lat_arr(:,0) = centre(:,0) - dist
  lat_arr(:,1) = centre(:,0) + dist
  lon_arr(:,0) = centre(:,1) - dist
  lon_arr(:,1) = centre(:,1) + dist

; Create arrays for data and output information
  llbox        = toint(dist*50)           ; Calculate domain size (grid points)
  dsize        = (/llbox,llbox/)          ; Domain size (grid points)
  d0	       = dsize(0)		  ; Dimension 0
  d1	       = dsize(1)		  ; Dimension 1

;==========================================================
; Read in 'pb' stream data and extract variables we need
;==========================================================

  b	       = addfile(fili_b,"r")

  do it	       = 0, numINPUT_b-4	  ; In 'pb' stream, it = 0 --> 16 UTC 2 July

 ; pb stream data starts 3 h after pa stream (see README.txt)
   bt 	       = it + 3 

   print("Working on time: "+time_str(it))
   itt = bt + 1
   title_arr(it) = "Valid at "+time_str(it)+" (T+"+itt+")"

 ; Get domain boundaries at each time from storm track (calculated off-line)
   t0          = lat_arr(bt,0)
   t1          = lat_arr(bt,1)
   n0          = lon_arr(bt,0)
   n1          = lon_arr(bt,1)

 ; Create latitude/longitude array subsets to reduce computation time
   lon10       = b->longitude_1({n0:n1})  ; longitude [d1 grid points]
   lat10       = b->latitude_1({t0:t1})   ; latitude  [d0 grid points]
   lon0        = b->longitude({n0:n1})    ; longitude [d1 grid points]
   lat0        = b->latitude({t0:t1})     ; latitude  [d0 grid points]

 ; Pressure levels (pres)
 ; (0) 1000, (1) 950, (2) 925, (3) 900, (4) 850, (5) 800, (6) 750, (7) 700
 ; (8) 650, (9) 600, (10) 500, (11) 400, (12) 300, (13) 250
 ; (14) 200, (15) 150, (16) 100 

 ; Either choose single pressure level or multiple levels (for layer average)
   if (lay .eq. 0) then
    pres       = b->p(p_cp)               ; Single pressure level
   elseif (lay .eq. 1) then
    pres       = b->p(p_cp:p_cp2)	  ; Multiple pressure levels
   end if 

   plevs       = dimsizes(pres)           ; Size of pressure level array
   pres@units  = "hPa"                    ; Pressure level units (for plotting)

   if (lay .eq. 1) then
    vort0      = b->rvor(it,p_cp:p_cp2,{t0:t1},{n0:n1})    ; Relative vorticity (/s)
    z0         = b->ht(it,p_cp:p_cp2,{t0:t1},{n0:n1})      ; Geopotential height (m)
   elseif (lay .eq. 0)
    vort0      = b->rvor(it,p_cp,{t0:t1},{n0:n1})          ; Relative vorticity (/s)
    z0         = b->ht(it,p_cp,{t0:t1},{n0:n1})            ; Geopotential height (m)
   end if

   vort    = vort0
   z       = z0
   lat1    = lat10
   lon1    = lon10
   lat     = lat0
   lon     = lon0

 ; Convert to preferred units 
   z 	   = z / 10	      ; Convert to decametres
   z@units = "dam"	      ;	Add metadata
   vort	   = vort * (10 ^ -6) ; Convert to s-1

;==============================================
; Analyse variables on chosen pressure level
;==============================================

   print("Working on: " + pres + " hPa")

   if (lay .eq. 1) then
    vort_plane    = dim_avg_n_Wrap(vort(:,:,:),0)
    geo_plane     = dim_avg_n_Wrap(z(:,:,:),0)
   elseif (lay .eq. 0) then 
    vort_plane    = vort(:,:)
    geo_plane     = z(:,:)
   end if

;====================================================
; Use centre already calculated using MSLP as base
;====================================================

   lt1		  = centre(bt,0) - 0.5
   lt2		  = centre(bt,0) + 0.5
   ln1		  = centre(bt,1) - 0.5
   ln2		  = centre(bt,1) + 0.5

 ; Find relative vorticity maximum close to previously-calculated storm centre
   vort_plane0    = vort_plane({lt1:lt2},{ln1:ln2})
   vort_max0	  = max(vort_plane0)
   dims_v0        = dimsizes(vort_plane0)
   vort_1d0       = ndtooned(vort_plane0)
   inds_v0        = ind_resolve(maxind(vort_1d0),dims_v0)

   lat_max_v0 	  = lt1 + (dy * inds_v0(0,0))
   lon_max_v0 	  = ln1 + (dx * inds_v0(0,1))
   print("Centre (vort): "+lat_max_v0+" degrees N, "+lon_max_v0+" degrees E")

   centre_new(0,it,0) = lat_max_v0
   centre_new(0,it,1) = lon_max_v0

 ; Also find geopotential height minimum close to previously-calculated storm centre
   geo_plane0     = geo_plane({lt1:lt2},{ln1:ln2})
   geo_min0       = min(geo_plane0)
   dims_p0        = dimsizes(geo_plane0)
   geo_1d0        = ndtooned(geo_plane0)
   inds_p0        = ind_resolve(minind(geo_1d0),dims_p0)

   lat_max_z0     = lt1 + (dy * inds_p0(0,0))
   lon_max_z0     = ln1 + (dx * inds_p0(0,1))   
   print("Centre (geo): "+lat_max_z0+" degrees N, "+lon_max_z0+" degrees E")

   centre_new(1,it,0) = lat_max_z0
   centre_new(1,it,1) = lon_max_z0

 ; Smooth geopotential height field and search again (ignore small-scale detail)
   geo_smth	  = smth9_Wrap(geo_plane0, 0.5, 0.5, True)
   geo_min_smth	  = min(geo_smth)
   dims_p1	  = dimsizes(geo_smth)
   geo_1d1	  = ndtooned(geo_smth)
   inds_p1	  = ind_resolve(minind(geo_1d1),dims_p1)

   lat_max_z1     = lt1 + (dy * inds_p1(0,0))
   lon_max_z1     = ln1 + (dx * inds_p1(0,1))
   print("Centre (geo_smth): "+lat_max_z1+" degrees N, "+lon_max_z1+" degrees E")

   centre_new(2,it,0) = lat_max_z1
   centre_new(2,it,1) = lon_max_z1   

 ; Calculate geopotential height max/min values for later contour plotting
   z_max0	  = toint(max(geo_plane))
   z_min0	  = toint(geo_min0)
   z_avg0	  = ( (z_max0 + z_min0) / 2)

;==========================================
; Now overlay centres onto contour plots
;==========================================

; Output file location and type
  output = "$sam/nepartak/images/basic/compare_centre_"+dat+"_"+ens0+"_"+time_arr(it)
  wks    = gsn_open_wks(opt,output)

;==============================================================
; Load NIMROD colour table for precip duration data
;==============================================================

; Load colour table 
  gsn_define_colormap(wks,"prcp_new")

; Create panel plot (MSLP and 10-m windspeed)
  panel = new(1,graphic)

;==========================
; Options for plotting
;==========================

; Relative vorticity
  opts_rvort                              = True
  opts_rvort@cnFillOn                     = True
  opts_rvort@cnLineLabelInterval          = 2.0
  opts_rvort@cnLineLabelFontHeightF       = 0.012
  opts_rvort@cnLineLabelBackgroundColor   = "transparent"
  opts_rvort@cnLineLabelPlacementMode     = "constant"
  opts_rvort@cnLinesOn                    = False ; Contour lines off
  opts_rvort@cnInfoLabelOn                = False ; Contour labels off
  opts_rvort@cnLevelSelectionMode         = "ExplicitLevels"
  opts_rvort@cnLevels                     = (/5., 10., 20., 30., 40., 50., \
                                              60., 80., 100., 120., 150./)
  opts_rvort@cnFillColors                 = (/0,4,5,6,7,8,9,10,11,12,13,14/)

  opts_rvort@gsnPaperOrientation          = "landscape"
  opts_rvort@tiMainString                 = ""
  opts_rvort@tiMainFontHeightF            = 0.0125
  opts_rvort@gsnLeftString                = ""
  opts_rvort@gsnRightString               = ""

  opts_rvort@lbLabelBarOn                 = True
  opts_rvort@lbBoxEndCapStyle             = "TriangleBothEnds"    ; Labelbar end shape
  opts_rvort@lbLabelFontHeightF           = 0.0125                ; Labelbar font size
  opts_rvort@lbLabelFont                  = "Helvetica"           ; Labelbar font
  opts_rvort@lbTitleString                = pres+" hPa relative vorticity (10~S~-4~N~ s~S~-1~N~)"
  opts_rvort@lbTitlePosition              = "Bottom"              ; Title position
  opts_rvort@lbTitleFontHeightF           = 0.015                 ; Title font height
  opts_rvort@pmLabelBarOrthogonalPosF     = 0.1                   ; Labelbar position
  opts_rvort@lbTitleOffsetF               = 0                     ; Lb title up/down
  opts_rvort@lbPerimOn                    = False                 ; Perimeter on/off

  opts_rvort@gsnDraw                      = False                 ; Do not draw plot
  opts_rvort@gsnFrame                     = False                 ; Do not adv. frame
  opts_rvort@gsnAddCyclic                 = False

  opts_rvort@mpLimitMode                  = "LatLon"

  if (zoom .eq. 1) then
   opts_rvort@mpMinLatF                   = centre(bt,0) - rad0
   opts_rvort@mpMinLonF                   = centre(bt,1) - rad0
   opts_rvort@mpMaxLatF                   = centre(bt,0) + rad0
   opts_rvort@mpMaxLonF                   = centre(bt,1) + rad0
  else
   opts_rvort@mpMinLatF                   = lat(0)
   opts_rvort@mpMinLonF                   = lon(0)
   opts_rvort@mpMaxLatF                   = lat(dimsizes(lat)-1)
   opts_rvort@mpMaxLonF                   = lon(dimsizes(lon)-1)
  end if

  opts_rvort@tmXTOn                       = "False"    ; No tickmarks on top x-axis
  opts_rvort@tmYROn                       = "False"    ; No tickmarks on right y-axis

  if (rad0 .le. 2.0) then 
   opts_rvort@tmXBMode			   = "Explicit"
   opts_rvort@tmXBValues	           = fspan(100.0, 150.0, 51)
   xb_labels0				   = sprintf("%0.1f~S~o~N~E",opts_rvort@tmXBValues)
   opts_rvort@tmXBLabels		   = xb_labels0

   opts_rvort@tmYLMode                     = "Explicit"
   opts_rvort@tmYLValues                   = fspan(0.0, 30.0, 31)
   yl_labels0                              = sprintf("%0.1f~S~o~N~N",opts_rvort@tmYLValues)
   opts_rvort@tmYLLabels                   = yl_labels0
  elseif (rad0 .gt. 2.0) then 
   opts_rvort@tmXBMode                     = "Explicit"
   opts_rvort@tmXBValues                   = fspan(100.0, 150.0, 26)
   xb_labels0                              = sprintf("%0.1f~S~o~N~E",opts_rvort@tmXBValues)
   opts_rvort@tmXBLabels                   = xb_labels0

   opts_rvort@tmYLMode                     = "Explicit"
   opts_rvort@tmYLValues                   = fspan(0.0, 30.0, 16)
   yl_labels0                              = sprintf("%0.1f~S~o~N~N",opts_rvort@tmYLValues)
   opts_rvort@tmYLLabels                   = yl_labels0
  end if 

  opts_rvort@mpGridAndLimbOn              = True       ; Grid lines on/off
  opts_rvort@mpGridLineDashPattern        = 2          ; Dash pattern
  opts_rvort@mpGridLatSpacingF            = 2.0        ; Spacing (latitude)
  opts_rvort@mpGridLonSpacingF            = 2.0        ; Spacing (longitude)

; Geopotential height
  opts_geo                                = True
  opts_geo@cnFillOn			  = False
  opts_geo@cnLineColor			  = "black"
  opts_geo@cnInfoLabelOn		  = False
  opts_geo@cnLineLabelsOn		  = True
  opts_geo@cnLineLabelInterval		  = 2.0
  opts_geo@cnLevelSelectionMode		  = "ExplicitLevels"
  opts_geo@cnLevels                       = ispan(z_min0, z_max0, 4)
  opts_geo@cnLineLabelPlacementMode	  = "constant"
  opts_geo@cnLineLabelPerimOn		  = False
  opts_geo@gsnContourLineThicknessesScale = 3.0
  opts_geo@tiMainString                   = ""
  opts_geo@gsnLeftString                  = ""
  opts_geo@gsnRightString                 = ""
  opts_geo@gsnDraw                        = False      ; Do not draw the plot
  opts_geo@gsnFrame                       = False      ; Do not advance the frame
  opts_geo@pmTickMarkDisplayMode          = "Always"   ; Turn on map tickmarks
  opts_geo@tmXBMajorLengthF               = 0.005      ; Change tickmark length
  opts_geo@tmXTOn                         = "False"    ; No tickmarks on top x-axis
  opts_geo@tmYROn                         = "False"    ; No tickmarks on right y-axis
  opts_geo@gsnMaximize                    = True       ; Maximise plot size

;=================
; Plot the data
;=================

  vort_plane = vort_plane * (10 ^ 4)
  vort_plane@units = "10~S~-4~N~ s~S~-1~N~"

  vort_plot = gsn_csm_contour_map(wks,vort_plane,opts_rvort)  ; Relative vorticity
  geo_plot  = gsn_csm_contour(wks,geo_plane,opts_geo)         ; Geopotential height
  overlay(vort_plot,geo_plot)
  panel(0)  = vort_plot

;==================================================
; Overlay markers for each storm centre estimate
;==================================================

; (1) Vorticity, (2) Geopotential height
  m_colours  = (/"purple3", "royalblue1", "orange4"/)
  m_indices  = (/12, 16, 9/)
  m_sizes    = (/10.0, 10.0, 10.0/)
  m_thick    = (/5.0, 5.0, 5.0/)

  markers    = new(dimsizes(m_colours),graphic)

  do i = 0, dimsizes(m_colours)-1
   mres0                    = True
   mres0@gsMarkerIndex      = m_indices(i)
   mres0@gsMarkerSizeF      = m_sizes(i)
   mres0@gsMarkerColor      = m_colours(i)
   mres0@gsMarkerThicknessF = m_thick(i)
   markers(i)               = gsn_add_polymarker(wks,vort_plot,centre_new(i,it,1),\
                                                 centre_new(i,it,0),mres0)
  end do

;===============================================
; Also overlay existing storm centre estimate
;===============================================

  marker = new(1,graphic)

  tc0 = NhlNewMarker(wks, "p", 137, 0.0, 0.0, 1.0, 1.0, 0.0)

  mres0               = True
  mres0@gsMarkerIndex = tc0         ; TC symbol
  mres0@gsMarkerSizeF = 10.0        ; Marker size
  mres0@gsMarkerColor = "red3"      ; Marker colour
  mres0@gsMarkerThicknessF = 5.0    ; Marker thickness (increase)

  marker    = gsn_add_polymarker(wks,vort_plot,centre(bt,1),centre(bt,0),mres0)

; Create single resources for legend 
  m_colour_all	      = (/m_colours(0),m_colours(1),m_colours(2),mres0@gsMarkerColor/)
  m_index_all	      = (/m_indices(0),m_indices(1),m_indices(2),mres0@gsMarkerIndex/)

;===============================
; Add and customise a legend
;===============================

  lg_opts                            = True
  lg_opts@pmLegendDisplayMode        = "Always"          ; Display a legend
  lg_opts@lgAutoManage               = False             ; Gives us control
  lg_opts@lgLabelFont                = "Helvetica"       ; Label font
  lg_opts@lgLabelFontHeightF         = 0.060             ; Height of labels

; Customise labels and line colours
  lg_opts@lgItemType                 = "Markers"
  lg_opts@lgMarkerColors             = m_colour_all
  lg_opts@lgMarkerIndexes            = m_index_all
  lg_opts@lgMarkerSizeF              = mres0@gsMarkerSizeF
  lg_opts@lgMarkerThicknessF         = mres0@gsMarkerThicknessF
  lsize                              = dimsizes(lg_opts@lgMarkerColors)

  lg_opts@lgLabelPosition            = "Right"           ; Label position within box
  lg_opts@lgItemPlacement            = "ExplicitPlacement"
  lg_opts@lgItemPositions            = fspan(0.12, 0.88, lsize)
  lg_opts@lgLeftMarginF		     = 0.01
  lg_opts@lgRightMarginF	     = 0.15

  lg_opts@vpWidthF                   = 0.125
  lg_opts@vpHeightF                  = 0.075
  lg_opts@lgPerimColor               = "black"
  lg_opts@lgPerimThicknessF          = 3.0
  lg_opts@lgPerimFill                = "SolidFill"
  lg_opts@lgPerimFillColor           = "white"

  lab_arr                            = (/"vort ("+pres+" hPa)", "hgt ("+pres+" hPa)", \
  				         "hgt_smth ("+pres+" hPa)", "mslp"/)
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

  annotate_1                         = gsn_add_annotation(vort_plot,lgd,am_opts)

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

  optsP                       = True
  optsP@gsnFrame              = False             ; Do not advance the frame
  optsP@gsnPanelLabelBar      = False             ; Turn off panel labelbar
  optsP@txString              = "Simulation "+ens0+" ("+dat+"): "+title_arr(it)
  optsP@gsnPanelFigureStrings = ""
  optsP@gsnMaximize           = True
  optsP@gsnPanelTop           = 0.90
  optsP@gsnPanelBottom        = 0.10
  optsP@amJust                = "TopLeft"
  optsP@gsnPanelFigureStringsFontHeightF = 0.0125 ; Label size (default 0.01)

  gsn_panel(wks,panel,(/1,1/),optsP)              ; Draw as a single plot
  frame(wks)

; Tidy up before next iteration
  delete([/vort,z,vort0,z0,lat10,lon10,lat1,lon1,lat0,lon0,lat,lon/])
  delete([/vort_plane,geo_plane,vort_plane0,geo_plane0,geo_smth,vort_1d0,geo_1d0,geo_1d1/])
  delete([/opts_geo@cnLevels/])

 end do                ; End loop over input times (do it = 0, numINPUT_b-1)

;===================================
; Write out values to text files
;===================================

  diri_out = "$ar/text/cp_new_"

;===========================
; MSLP and 10-m windspeed
;===========================
 
; MSLP
  slp_out = diri_out+dat+"_"+ens0+"_slp.txt"
  asciiwrite(slp_out,slp_box)

; Windspeed
  vec_out  = diri_out+dat+"_"+ens0+"_10m.txt"
  asciiwrite(vec_out,vec_box)

;================================================
; Storm centre position using multiple methods 
;================================================

; (1) Original method 

; Latitude 
  lat_out  = diri_out+dat+"_"+ens0+"_lat.txt"
  asciiwrite(lat_out,centre(:,0))

; Longitude
  lon_out  = diri_out+dat+"_"+ens0+"_lon.txt"
  asciiwrite(lon_out,centre(:,1))

; (2) Relative vorticity

; Latitude
  lat_out  = diri_out+dat+"_"+ens0+"_lat_rv.txt"
  asciiwrite(lat_out,centre_new(0,:,0))

; Longitude
  lon_out  = diri_out+dat+"_"+ens0+"_lon_rv.txt"
  asciiwrite(lon_out,centre_new(0,:,1))

; (3) Geopotential height

; Latitude
  lat_out  = diri_out+dat+"_"+ens0+"_lat_z0.txt"
  asciiwrite(lat_out,centre_new(1,:,0))

; Longitude
  lon_out  = diri_out+dat+"_"+ens0+"_lon_z0.txt"
  asciiwrite(lon_out,centre_new(1,:,1))

; (4) Geopotential height (smoothed)

; Latitude
  lat_out  = diri_out+dat+"_"+ens0+"_lat_z1.txt"
  asciiwrite(lat_out,centre_new(2,:,0))

; Longitude
  lon_out  = diri_out+dat+"_"+ens0+"_lon_z1.txt"
  asciiwrite(lon_out,centre_new(2,:,1))

end 