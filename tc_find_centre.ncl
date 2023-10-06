; Script to compute the storm track for a TC using mean sea level pressure

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_rm.ncl"

; ncl 'dat="02T12"' int=1 plt=0 'opt="x11"' 'ens0="em11"' rad=3.0 tc_find_centre.ncl

begin

;==============================================
; Find 'pa' stream data to read in and plot
;==============================================

  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0+"/vc"
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0
  fili_a1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pa*all.nc")
  fili_a      = diri+"/"+fili_a1

  numFILES    = dimsizes(fili_a)		; 13 input files (T+54 to T+66)
  numINPUT    = numFILES * 12 			; 156 input time intervals (12 per input file)

; Extract string for each ensemble member
  str1    = str_split(fili_a(0),"_")
  ens_str = str1(5)

  print("Working on ensemble member "+ens_str)

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 10000000000
  end setvalues

; Create arrays to hold values of storm position, minimum MSLP and max 10-m windspeed
  slp_box  = new(numINPUT,float)  		  ; Minimum MSLP
  vec_box  = new(numINPUT,float)		  ; Maximum 10-m windspeed
  centre   = new((/numINPUT,2/),"float")	  ; Storm centre (lat/lon)

; Set counter variable (for writing above diagnostics to text files)
  ct       = 0

; Create arrays to hold date strings 
  date_str  = new(numINPUT,string)
  time_str  = new(numINPUT,string)
  time_arr  = new(numINPUT,string)
  title_arr = new(numINPUT,string)

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
  if (dat .eq. "02T12") then      	    ; T+54
   lat0     = lat_bt0(9)
   lon0     = lon_bt0(9)
  elseif (dat .eq. "03T00") then	    ; T+??
   lat0     = lat_bt0(3)
   lon0     = lon_bt0(3)
  end if 

;==========================================================
; Get the variables we will need (loop over input files)
;==========================================================

  do nf = 0, numFILES-1
  f = addfile(fili_a(nf),"r")
  
 ;================================================
 ; Number of times in different variable arrays
 ;================================================

   time  = f->t_1           ; Times in file [12]
   times = 12

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
   mslp  = f->p

   mslp  = mslp / 100      ; Convert to hPa
   mslp@units = "hPa"

   dx    = lon(1) - lon(0) ; Grid spacing (longitude)
   dy    = lat(1) - lat(0) ; Grid spacing (latitude)

;===========================
; Loop over times in file
;===========================

   tf = times-1

   do it = 0, tf, int        ; Loop over times in each file

    ctt = ct + 1 
    print("Working on forecast hour: "+ctt)
    title = ctt+" h forecast, initialised at 12 UTC 2 July"

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

   hr_min	= hour(it)+min_arr(it)

   time_arr(ct) = sprinti("%0.2i", day(it)) + month_abbr(month(it)) + \
                  "_" + hr_min + "Z"
   time_str(ct) = hr_min + " UTC " + sprinti("%0.2i ", day(it)) + month_abbr(month(it))

   print("Working on time: "+time_str(ct))
   title_arr(ct) = "Valid at "+time_str(ct)+" (T+"+ctt+")"

;===============================================
; Find cyclone centre on RA1T grid (use MSLP)
;===============================================

 ; Read in mean sea level pressure and 10-m windspeed data 
 ; At the initial time, focus on the region close to the best track location

   if (ct .eq. 0) then 

  ; Define small grid subset centred on early best track location 	
    latS  	     = lat0-1.0
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
    lat_minslp = latS + (dy * inds(0,0))
    lon_minslp = lonS + (dx * inds(0,1))
   else 
    lat_minslp = lat(0) + (dy * inds(0,0))
    lon_minslp = lon(0) + (dx * inds(0,1))
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
      print("Change in latitude = "+d_lat+" deg")
      print("Change in longitude = "+d_lon+" deg")

      print_clock("Now finding correct centre based upon previous storm position...")

    ; Create smaller grid [0.5 degrees]
      lt1 = lt0 - 0.5
      lt2 = lt0 + 0.5
      ln1 = ln0 - 0.5
      ln2 = ln0 + 0.5

      slp_plane0 = slp_plane({lt1:lt2},{ln1:ln2})	 ; Smaller grid (1.0 x 1.0 deg)
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

   slp_box(ct) = min(slp_plane({lat_minslp-rad:lat_minslp+rad},\
			       {lon_minslp-rad:lon_minslp+rad}))
   vec_box(ct) = max(speed_plane({lat_minslp-rad:lat_minslp+rad},\
    		                 {lon_minslp-rad:lon_minslp+rad}))

;==================================================
; Plot data over whole domain (skip if required)
;==================================================

   if (plt .eq. 1) then 

  ; Output file location and type
    output = "$sam/nepartak/images/basic/track_5min_"+dat+"_"+ens0+"_"+time_arr(ct)
    wks    = gsn_open_wks(opt,output)

;==============================================================
; Load NIMROD colour table for precip duration data
;==============================================================

  ; Load colour table
    gsn_define_colormap(wks,"radar_new")

  ; Change the colour map around & make first colour transparent
    cmap = read_colormap_file("radar_new")
    cmap(0,3) = 0.0

  ; Create panel plot (MSLP and 10-m windspeed)
    panel = new(1,graphic)

;==========================
; Options for plotting
;==========================

  ; 10-m windspeed 
    opts_speed                              = True
    opts_speed@cnFillOn                     = True
    opts_speed@cnLineLabelInterval          = 2.0
    opts_speed@cnLineLabelFontHeightF       = 0.012
    opts_speed@cnLineLabelBackgroundColor   = "transparent"
    opts_speed@cnLineLabelPlacementMode     = "constant"
    opts_speed@cnLinesOn                    = False ; Contour lines off
    opts_speed@cnInfoLabelOn                = False ; Contour labels off
    opts_speed@cnLevelSelectionMode         = "ExplicitLevels"
    opts_speed@cnFillPalette                = cmap
    opts_speed@cnLevels                     = (/15.0, 20.0,\
                                                25.0, 30.0, 35.0, 40.0, \
                                                45.0, 50.0, 55.0, 60.0, \
                                                65.0/)
    opts_speed@cnFillColors                 = (/0,4,5,6,7,8,9,10,11,12,13,14/)

    opts_speed@gsnPaperOrientation          = "landscape"
    opts_speed@tiMainString                 = ""
    opts_speed@tiMainFontHeightF            = 0.0125
    opts_speed@gsnLeftString                = ""
    opts_speed@gsnRightString               = ""

    opts_speed@lbLabelBarOn                 = True
    opts_speed@lbBoxEndCapStyle             = "TriangleBothEnds"    ; Labelbar end shape
    opts_speed@lbLabelFontHeightF           = 0.0125                ; Labelbar font size
    opts_speed@lbLabelFont                  = "Helvetica"           ; Labelbar font
    opts_speed@lbTitleString                = "Wind speed (m s~S~-1~N~)"
    opts_speed@lbTitlePosition              = "Bottom"              ; Title position
    opts_speed@lbTitleFontHeightF           = 0.015                 ; Title font height
    opts_speed@pmLabelBarOrthogonalPosF     = 0.1                   ; Labelbar position
    opts_speed@lbTitleOffsetF               = 0                     ; Lb title up/down
    opts_speed@lbPerimOn                    = False                 ; Perimeter on/off
    opts_speed@gsnDraw                      = False                 ; Do not draw plot
    opts_speed@gsnFrame                     = False                 ; Do not adv. frame
    opts_speed@gsnAddCyclic                 = False

    opts_speed@mpLimitMode                  = "LatLon"

    if (zoom .eq. 1) then 
     opts_speed@mpMinLatF                    = centre(ct,0) - rad
     opts_speed@mpMinLonF                    = centre(ct,1) - rad
     opts_speed@mpMaxLatF                    = centre(ct,0) + rad
     opts_speed@mpMaxLonF                    = centre(ct,1) + rad
    else
     opts_speed@mpMinLatF                   = lat(0)
     opts_speed@mpMinLonF                   = lon(0)
     opts_speed@mpMaxLatF                   = lat(dimsizes(lat)-1)
     opts_speed@mpMaxLonF                   = lon(dimsizes(lon)-1)
    end if 

    opts_speed@tmXTOn                       = "False"    ; No tickmarks on top x-axis
    opts_speed@tmYROn                       = "False"    ; No tickmarks on right y-axis
    opts_speed@mpGridAndLimbOn              = True       ; Grid lines on/off
    opts_speed@mpGridLineDashPattern        = 2          ; Dash pattern
    opts_speed@mpGridLatSpacingF            = 2.0        ; Spacing (latitude)
    opts_speed@mpGridLonSpacingF            = 2.0        ; Spacing (longitude)

  ; Mean sea level pressure resources
    opts_mslp                               = True
    opts_mslp@cnFillOn                      = False
    opts_mslp@cnLineColor                   = "black"
    opts_mslp@cnLevelSelectionMode          = "ExplicitLevels"
    opts_mslp@cnLevels                      = ispan(880,1020,4)
    opts_mslp@cnInfoLabelOn                 = False
    opts_mslp@cnLineThicknessF              = 3.0
    opts_mslp@cnLineLabelsOn                = False
    opts_mslp@gsnDraw                       = False
    opts_mslp@gsnFrame                      = False
    opts_mslp@gsnAddCyclic                  = False
    opts_mslp@gsnLeftString                 = ""
    opts_mslp@gsnRightString                = ""

;=================
; Plot the data
;=================

    spd_plot  = gsn_csm_contour_map(wks,speed_plane,opts_speed) ; 10-m wind speed
    slp_plot  = gsn_csm_contour(wks,slp_plane,opts_mslp)      	; Mean sea level pressure
    overlay(spd_plot,slp_plot)
    panel(0)  = spd_plot

    marker    = new(1,graphic)

  ; Add NH tropical cyclone symbol to index
    tc0 = NhlNewMarker(wks, "p", 137, 0.0, 0.0, 1.0, 1.0, 0.0)

    mres0               = True
    mres0@gsMarkerIndex = tc0         ; TC symbol
    mres0@gsMarkerSizeF = 20.0        ; Marker size
    mres0@gsMarkerColor = "red3"      ; Marker colour
    mres0@gsMarkerThicknessF = 5.0    ; Marker thickness (increase)

    marker    = gsn_add_polymarker(wks,spd_plot,lon_minslp,lat_minslp,mres0)

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    optsP                       = True
    optsP@gsnFrame              = False             ; Do not advance the frame
    optsP@gsnPanelLabelBar      = False             ; Turn off panel labelbar
    optsP@txString              = "Nepartak (2-9 July 2016): "+title_arr(ct)
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
    delete([/slp_plane,u_plane,v_plane,speed_plane/])
    delete([/utc_date,year,month,day,hour,minute,second/])

    ct = ct + 1 ; Counter variable (time)

   end do     ; End time loop (do it = 0, times-1, int)

 ; Tidy up before next iteration 
   delete([/time,times,u,v,mslp/])

 end do	      ; End input file loop (do nf = 0, numFILES-1)

;===================================
; Write out values to text files
;===================================

    diri_out = "$ar/text/cp_5min_"

  ; 'ens_str' is the number of the ensemble member (000 -- 011)

  ; MSLP
    slp_out = diri_out+dat+"_"+ens_str+"_slp.txt"
    asciiwrite(slp_out,slp_box)

  ; Windspeed
    vec_out  = diri_out+dat+"_"+ens_str+"_10m.txt"
    asciiwrite(vec_out,vec_box)

  ; Latitude 
    lat_out  = diri_out+dat+"_"+ens_str+"_lat.txt"
    asciiwrite(lat_out,centre(:,0))

  ; Longitude
    lon_out  = diri_out+dat+"_"+ens_str+"_lon.txt"
    asciiwrite(lon_out,centre(:,1))

end 