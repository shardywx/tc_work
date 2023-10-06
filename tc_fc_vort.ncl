; Script to compute the storm track for a TC using relative vorticity

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_rm.ncl"

; ncl dat=\"02T12\" int=1 plt=0 opt=\"x11\" ens0=\"em11\" rad=2.5 trb=0 p0=3 dist=3.0 tc_fc_vort.ncl

begin

;==============================================
; Find 'pb' stream data to read in and plot 
;==============================================

  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0+"/vc"
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0
  fili_b1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pb.nc")
  fili_b      = diri+"/"+fili_b1

; Get number of times in file using built-in function ('nc_times')
  numINFO     = nc_times(fili_b)
  numTIMES    = numINFO[0]		; Number of times
  time	      = numINFO[1]		; Array containing timestamp info at each time]
  print("Total number of time intervals = "+numTIMES)

;============================================================
; Calculate storm motion using built-in function ('st_rm')
;============================================================

; 'ens0'    = ensemble simulation
; 'diri'    = path to input files (see above)
; 'dist'    = size of box following storm (degrees)
; 'mins'    = analyse 1-h (0) or 5-min (1) data

  storm_rel = st_rm(dat, ens0, "$ar/text/", dist, 1)
  u_cyc     = storm_rel[0]             ; Zonal wind
  v_cyc     = storm_rel[1]             ; Meridional wind
  vel_cyc   = storm_rel[2]             ; Vector wind
  lat_arr   = storm_rel[4]             ; Latitude array subset (following storm)
  lon_arr   = storm_rel[5]             ; Longitude array subset (following storm)
  centre    = storm_rel[6]             ; Storm track (position) array

; Create arrays for data and output information
  llbox     = toint(dist*50)           ; Calculate domain size (grid points)
  dsize     = (/llbox,llbox/)          ; Domain size (grid points)
 
; Extract string for each ensemble member
  str1      = str_split(fili_b(0),"_")
  ens_str   = str1(5)

  print("Working on ensemble member "+ens_str)

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 10000000000
  end setvalues

;=====================================================
; Create correct date strings for each output time
;=====================================================

  month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
                    "Jul","Aug","Sep","Oct","Nov","Dec"/)

  utc_date = cd_calendar(time, 0)

  year   = tointeger(utc_date(:,0))
  month  = tointeger(utc_date(:,1))
  day    = tointeger(utc_date(:,2))
  hour   = tointeger(utc_date(:,3))
  minute = tointeger(utc_date(:,4))
  second = utc_date(:,5)

; Create array to hold finished date strings
  hr_min = new(numTIMES,string)

  min0   = (/4,9,14,19,24,29,34,39,44,49,54,59/)      ; Minute array (before)
  min1   = (/5,10,15,20,25,30,35,40,45,50,55,0/)      ; Minute array (after)

; Before correcting, find the indices of the elements where 'minute' = 59
  hr_ind = ind(minute .eq. 59)

; Correct for errors in the code (round up values of 'minute' to multiples of 5)
  ind_f = numTIMES-1

  do i = 0, ind_f
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

; Create string array of all values in 'minute'
  mins = tostring(minute)

; Create 'hhmm' strings for each time (for output and title strings)
  do i = 0, ind_f

   if (hour(i) .lt. 10) then                  ; Add zeros in front of single-digit hours

    if (minute(i) .lt. 10) then
     hr_min(i) = "0"+hour(i)+"0"+mins(i)
    elseif (minute(i) .ge. 10) then
     hr_min(i) = "0"+hour(i)+mins(i)
    end if

   else

    if (minute(i) .lt. 10) then
     hr_min(i) = hour(i)+"0"+mins(i)          ; Also add zeros in front of single-digit minutes
    elseif (minute(i) .ge. 10) then
     hr_min(i) = hour(i)+mins(i)
    end if

   end if

  end do

  time_str   = new(numTIMES,string)
  time_arr   = new(numTIMES,string)
  title_arr  = new(numTIMES,string)

  time_str   = hr_min + " UTC " + sprinti("%0.2i ", day) + month_abbr(month)
  time_arr   = sprinti("%0.2i", day) + month_abbr(month) + "_" + hr_min + "Z"
  title_arr  = "Valid at "+time_str

;===================
; More prelashing
;===================

; Create arrays to hold values of maximum vorticity and storm centre position
  vort_box = new(numTIMES,float)  		  ; Maximum 900 hPa vorticity
  centre   = new((/numTIMES,2/),"float")	  ; Storm centre (lat/lon)

;=================================================
; Read in best track information from text file
;=================================================

; Number of times in IBTrACS dataset for Nepartak
  dsn = 28

; Read in all latitude/longitude/date information
  lat_bt0  = asciiread("$sam/um/cp/nepartak/lat_nep.txt",(/dsn/),"float")
  lon_bt0  = asciiread("$sam/um/cp/nepartak/lon_nep.txt",(/dsn/),"float")
  dat_bt0  = asciiread("$sam/um/cp/nepartak/dat_nep.txt",(/dsn/),"string")

; From these arrays, extract storm position at start of VC analysis period
; Use these data to constrain the initial position of the storm in the UM analysis
  if (dat .eq. "02T12") then  	  	    ; T+1 (13 UTC 2 July)
   lat0     = lat_bt0(1)
   lon0     = lon_bt0(1)
  elseif (dat .eq. "03T00") then	    ; T+1 (01 UTC 3 July)
   lat0     = lat_bt0(3)
   lon0     = lon_bt0(3)
  end if 

;======================
; Read in input file 
;======================

  f = addfile(fili_b,"r")
  
;==================================
; Get the variables we will need
;==================================

  do it = 0, numTIMES-1

   print("Working on time: "+time_str(it))

;=======================
; All other variables 
;=======================

 ; Get domain boundaries at each time from storm track (calculated off-line using MSLP)
   t0    = lat_arr(it,0)
   t1    = lat_arr(it,1)
   n0    = lon_arr(it,0)
   n1    = lon_arr(it,1)

 ; Create latitude/longitude array subsets (dim = dsize) to reduce computation time
   lon   = f->longitude({n0:n1})    ; longitude
   lat   = f->latitude({t0:t1})     ; latitude
   lon1  = f->longitude_1({n0:n1})  ; longitude
   lat1  = f->latitude_1({t0:t1})   ; latitude

 ; Pressure levels (pres)
 ; (0) 1000, (1) 950, (2) 925, (3) 900, (4) 850, (5) 800, (6) 750, (7) 700
 ; (8) 650, (9) 600, (10) 500, (11) 400, (12) 300, (13) 250
 ; (14) 200, (15) 150, (16) 100

 ; Read in pressure levels
   pres  = f->p
   pres@units = "hPa"
   plevs = dimsizes(pres)

 ; Read in relative vorticity on single pressure level
   vort  = f->rvor(it,p0,{t0:t1},{n0:n1})
   print("Working on pressure level: "+pres(p0)+ " hPa")

   dx    = lon(1) - lon(0)  ; Grid spacing (longitude)
   dy    = lat(1) - lat(0)  ; Grid spacing (latitude)

 ; Convert units 
   vort  = vort * (10 ^ -6) ; Convert to s-1

;====================================================
; Find cyclone centre on RA1T grid (use vorticity)
;====================================================

 ; Read in mean sea level pressure and 10-m windspeed data 
 ; At the initial time, focus on the region close to the best track location

   if (it .eq. 0) then 

  ; Define small grid subset centred on early best track location 	
    latS  	     = lat0-1.0
    latF	     = lat0+1.0
    lonS	     = lon0-1.0
    lonF	     = lon0+1.0
    vort_plane       = vort({latS:latF},{lonS:lonF})	; Smaller region
   else
    vort_plane       = vort(:,:)			; Whole grid 
   end if 

   dims          = dimsizes(vort_plane)
   vort1d        = ndtooned(vort_plane)
   inds          = ind_resolve(maxind(vort1d),dims)

   if (it .eq. 0) then 
    lat_max_rv = latS + (dy * inds(0,0))
    lon_max_rv = lonS + (dx * inds(0,1))
   else 
    lat_max_rv = lat(0) + (dy * inds(0,0))
    lon_max_rv = lon(0) + (dx * inds(0,1))
   end if 

   delete([/dims,vort1d,inds/])

;========================================================
; EMPLOY A SAFETY NET IF WE IDENTIFY THE WRONG CYCLONE
;========================================================

   if (it .eq. 0) then
    centre(it,0) = lat_max_rv
    centre(it,1) = lon_max_rv

   else
    lt0  = centre(it-1,0)             ; Previous TC latitude
    ln0  = centre(it-1,1)             ; Previous TC longitudem
    d_lat = abs(lat_max_rv-lt0)       ; Change in latitude (t1 - t0)
    d_lon = abs(lon_max_rv-ln0)       ; Change in longitude

    if (d_lat .gt. 0.75 .or. d_lon .gt. 0.75) then

     if (trb .eq. 1) then 
      print("Previous storm latitude ("+lt0+") ; longitude ("+ln0+")")
      print("Change in latitude = "+d_lat+" deg")
      print("Change in longitude = "+d_lon+" deg")
     end if 

      print_clock("Calculating new storm centre...")

    ; Create smaller grid [0.5 degrees]
      lt1 = lt0 - 0.5
      lt2 = lt0 + 0.5
      ln1 = ln0 - 0.5
      ln2 = ln0 + 0.5

      rv_plane0  = vort_plane({lt1:lt2},{ln1:ln2})      ; Smaller grid (1.0 x 1.0 deg)
      rv_max0    = max(rv_plane0)                       ; Find minimum on smaller grid
      dims_v0    = dimsizes(rv_plane0)                  ; Size of smaller grid
      rv_1d0     = ndtooned(rv_plane0)                  ; Create 1-D array
      inds_v0    = ind_resolve(maxind(rv_1d0),dims_v0)  ; Find index of minimum

      lat_max_rv = lt1 + (dy * inds_v0(0,0))
      lon_max_rv = ln1 + (dy * inds_v0(0,1))
      print("New centre: "+lat_max_rv+" deg N, "+lon_max_rv+" deg E")
      centre(it,0) = lat_max_rv
      centre(it,1) = lon_max_rv

      delete([/rv_plane0, rv_max0, dims_v0, rv_1d0, inds_v0/])

    else
      print("Cyclone centre: "+lat_max_rv+" deg N, "+lon_max_rv+" deg E")
      centre(it,0) = lat_max_rv
      centre(it,1) = lon_max_rv
    end if

  end if

;==================================================
; Plot data over whole domain (skip if required)
;==================================================

   if (plt .eq. 1) then 

  ; Output file location and type
    output = "$sam/nepartak/images/basic/track_vort_5min_"+dat+"_"+ens0+"_"+time_arr(it)
    wks    = gsn_open_wks(opt,output)

;==============================================================
; Load NIMROD colour table for precip duration data
;==============================================================

  ; Load colour table
    gsn_define_colormap(wks,"prcp_new")

  ; Create panel plot (relative vorticity)
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
    opts_rvort@lbTitleString                = pres(p0)+ \
    					      " hPa relative vorticity (10~S~-4~N~ s~S~-1~N~)"
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
     opts_rvort@mpMinLatF                   = centre(it,0) - rad
     opts_rvort@mpMinLonF                   = centre(it,1) - rad
     opts_rvort@mpMaxLatF                   = centre(it,0) + rad
     opts_rvort@mpMaxLonF                   = centre(it,1) + rad
    else
     opts_rvort@mpMinLatF                   = lat(0)
     opts_rvort@mpMinLonF                   = lon(0)
     opts_rvort@mpMaxLatF                   = lat(dimsizes(lat)-1)
     opts_rvort@mpMaxLonF                   = lon(dimsizes(lon)-1)
    end if 

    opts_rvort@tmXTOn                       = "False"    ; No tickmarks on top x-axis
    opts_rvort@tmYROn                       = "False"    ; No tickmarks on right y-axis
    opts_rvort@mpGridAndLimbOn              = True       ; Grid lines on/off
    opts_rvort@mpGridLineDashPattern        = 2          ; Dash pattern
    opts_rvort@mpGridLatSpacingF            = 2.0        ; Spacing (latitude)
    opts_rvort@mpGridLonSpacingF            = 2.0        ; Spacing (longitude)

;=================
; Plot the data
;=================

    vort_plane = vort_plane * (10 ^ 4)
    vort_plane@units = "10~S~-4~N~ s~S~-1~N~"
   
    rv_plot   = gsn_csm_contour_map(wks,vort_plane,opts_rvort) ; Relative vorticity
    panel(0)  = rv_plot

    marker    = new(1,graphic)

  ; Add NH tropical cyclone symbol to index
    tc0 = NhlNewMarker(wks, "p", 137, 0.0, 0.0, 1.0, 1.0, 0.0)

    mres0               = True
    mres0@gsMarkerIndex = tc0         ; TC symbol
    mres0@gsMarkerSizeF = 20.0        ; Marker size
    mres0@gsMarkerColor = "red3"      ; Marker colour
    mres0@gsMarkerThicknessF = 5.0    ; Marker thickness (increase)

    marker    = gsn_add_polymarker(wks,rv_plot,centre(it,1),centre(it,0),mres0)

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    optsP                       = True
    optsP@gsnFrame              = False             ; Do not advance the frame
    optsP@gsnPanelLabelBar      = False             ; Turn off panel labelbar
    optsP@txString              = "Typhoon Nepartak (2016): "+title_arr(it)
    optsP@gsnPanelFigureStrings = ""
    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.90
    optsP@gsnPanelBottom        = 0.10
    optsP@amJust                = "TopLeft"
    optsP@gsnPanelFigureStringsFontHeightF = 0.0125 ; Label size (default 0.01)

    gsn_panel(wks,panel,(/1,1/),optsP)              ; Draw as a single plot
    frame(wks)

   end if		; Only plot images if required 

  ; Tidy up 
    delete([/vort_plane,vort,lat,lon,lat1,lon1/])

   end do     ; End time loop (do it = 0, times-1, int)

;=================================================
; Write out storm centre position to text files
;=================================================

    diri_out = "$ar/text/cp_5min_"

  ; 'ens_str' is the number of the ensemble member (000 -- 011)

  ; Latitude 
    lat_out  = diri_out+dat+"_"+ens_str+"_lat_rv.txt"
    asciiwrite(lat_out,centre(:,0))

  ; Longitude
    lon_out  = diri_out+dat+"_"+ens_str+"_lon_rv.txt"
    asciiwrite(lon_out,centre(:,1))

end 