; Script to compute the storm track for a TC using MSLP and relative vorticity

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_rm.ncl"

; ncl dat=\"02T12\" int=1 opt=\"x11\" rad0=3.0 dist=3.0 vc=0 
; lay=0 p_cp=3 zoom=1 lgd0=\"tr\" tc_fc_haima.ncl

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

  diri        = "/nfs/a37/earshar/um/cp/haima/"
  fili_prefix = "201610"+dat+"00Z_HAIMA_4p4_L80_ra1t_"

; Find 'pa' stream data (horizontal windspeed)
  fili_a1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"pvera???_new.nc")
  fili_a      = diri+"/"+fili_a1

; Find 'pb' stream (MSLP)
  fili_b1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"pverb???_new.nc")
  fili_b      = diri+"/"+fili_b1

  numFILES    = dimsizes(fili_a)		; 16 input files (T+006 to T+096)
  numINPUT    = numFILES * 6     		; 96 input time intervals (16 * 6)

  print("Working on initialisation time "+dat)

  setvalues NhlGetWorkspaceObjectId
  "wsMaximumSize" : 10000000000
  end setvalues

; Create arrays to hold values of storm position, minimum MSLP and max 10-m windspeed
  slp_box   = new(numINPUT,float)  		  ; Minimum MSLP
  vec_box   = new(numINPUT,float)		  ; Maximum 10-m windspeed
  centre    = new((/numINPUT,2/),"float")	  ; Storm centre [numINPUT * 2]

; Set counter variable (for writing above diagnostics to text files)
  ct        = 0

; Create arrays to hold date strings 
  time_str  = new(numINPUT,string)
  time_arr  = new(numINPUT,string)
  title_arr = new(numINPUT,string)

; Create initialisation string
  if (dat .eq. "15T00") then 
   date_str = "00 UTC 15 Oct"
  elseif (dat .eq. "15T12") then 
   date_str = "12 UTC 15 Oct"
  elseif (dat .eq. "16T00") then 
   date_str = "00 UTC 16 Oct"
  elseif (dat .eq. "16T12") then 
   date_str = "12 UTC 16 Oct"
  elseif (dat .eq. "17T00") then
   date_str = "00 UTC 17 Oct"
  elseif (dat .eq. "17T12") then
   date_str = "12 UTC 17 Oct"
  elseif (dat .eq. "18T00") then
   date_str = "00 UTC 18 Oct"
  elseif (dat .eq. "18T12") then
   date_str = "12 UTC 18 Oct"
  end if  

;=================================================
; Read in best track information from text file
;=================================================

; Number of times in IBTrACS dataset for Haima (28 for Nepartak)
  dsn 	    = 30

; Read in all latitude/longitude/date information
  lat_bt0   = asciiread("$sam/um/cp/haima/lat_hai.txt",(/dsn/),"float")
  lon_bt0   = asciiread("$sam/um/cp/haima/lon_hai.txt",(/dsn/),"float")
  dat_bt0   = asciiread("$sam/um/cp/haima/dat_hai.txt",(/dsn/),"string")

; From these arrays, extract storm position at start of VC analysis period
; Use these data to constrain the initial position of the storm in the UM analysis
  if (dat .eq. "15T00") then
   lat0     = lat_bt0(1)
   lon0     = lon_bt0(1)
  elseif (dat .eq. "15T12") then
   lat0     = lat_bt0(3)
   lon0     = lon_bt0(3)
  elseif (dat .eq. "16T00") then
   lat0     = lat_bt0(5)
   lon0     = lon_bt0(5)
  elseif (dat .eq. "16T12") then
   lat0     = lat_bt0(7)
   lon0     = lon_bt0(7)
  elseif (dat .eq. "17T00") then
   lat0     = lat_bt0(9)
   lon0     = lon_bt0(9)
  elseif (dat .eq. "17T12") then
   lat0     = lat_bt0(11)
   lon0     = lon_bt0(11)
  elseif (dat .eq. "18T00") then
   lat0     = lat_bt0(13)
   lon0     = lon_bt0(13)
  elseif (dat .eq. "18T12") then
   lat0     = lat_bt0(15)
   lon0     = lon_bt0(15)
  end if

;==========================================================
; Get the variables we will need (loop over input files)
;==========================================================

  do nf     = 0, numFILES-2
  a  	    = addfile(fili_a(nf),"r")
  b         = addfile(fili_b(nf),"r")

;=======================
; All other variables 
;=======================

  lon   = a->longitude    ; longitude (1098 points --> 109.04 to 152.92 degrees E)
  lat   = a->latitude     ; latitude  (810 points --> 1.8 to 34.16 degrees N)
  lon1  = a->longitude_1  ; longitude (1098 points --> 109.02 to 152.90 degrees E)
  lat1  = a->latitude_1   ; latitude  (811 points --> 1.78 to 34.18 degrees N)

  u_varname = "x-wind"
  v_varname = "y-wind"

  u0    = a->$u_varname$  ; 10-m zonal velocity (m/s) [t_3 | z | lat_1 | lon_1]
  v0    = a->$v_varname$  ; 10-m meridional velocity (m/s) [t_3 | z | lat_1 | lon_1]

; Remove extra record dimensions we added with 'ncecat'
  u     = u0(0,:,0,:,:)
  v     = v0(0,:,0,:,:)

; Read in mean sea level pressure (Pa) [t | z | lat | lon]
; For 1-h data (vc = 0), ignore first time interval (nf = 0)
  mslp0 = b->p

; Remove record dimension
  mslp  = mslp0(0,:,0,:,:)

  mslp  = mslp / 100      ; Convert to hPa
  mslp@units = "hPa"
 
  dx    = lon(1) - lon(0) ; Grid spacing (longitude)
  dy    = lat(1) - lat(0) ; Grid spacing (latitude)

; Also read in time arrays 
  time_a = a->t_1         ; 'pa' stream
  time_b = b->t           ; 'pb' stream

; Number of times in file 
  numTIMES = dimsizes(time_a)

; Tidy up before continuing
  delete([/u0,v0,mslp0/])

;===========================
; Loop over times in file
;===========================

  do it = 0, numTIMES-1                     ; Loop over all input times 

   ctt = ct + 7 
   print("Working on forecast hour: "+ctt)
   title = ctt+" h forecast, initialised at "+date_str

;=====================================================
; Create correct date strings for each output time
;=====================================================

  month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
	         "Jul","Aug","Sep","Oct","Nov","Dec"/)

  utc_date   = cd_calendar(time_a, 0)
  year       = tointeger(utc_date(:,0))
  month      = tointeger(utc_date(:,1))
  day        = tointeger(utc_date(:,2))
  hour       = tointeger(utc_date(:,3))
  minute     = tointeger(utc_date(:,4))
  second     = utc_date(:,5)

; Correct for errors in the code (round up value of hour when minutes = 59)
  do tt = 0, numTIMES-1
   if (minute(tt) .gt. 30) then
     hour(tt) = hour(tt)+1
   end if
  end do
  delete(tt)

; Correct for errors in 'minutes'
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
   u_plane       = u(it,{latS:latF},{lonS:lonF})
   v_plane       = v(it,{latS:latF},{lonS:lonF})
   slp_plane     = mslp(it,{latS:latF},{lonS:lonF})

  else

   u_plane       = u(it,:,:)
   v_plane       = v(it,:,:)
   slp_plane     = mslp(it,:,:)

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

    slp_plane0 = slp_plane({lt1:lt2},{ln1:ln2})       ; Smaller grid (1.0 x 1.0 deg)
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

  r = 1.0

  slp_box(ct) = min(slp_plane({centre(ct,0)-r:centre(ct,0)+r},\
			    {centre(ct,1)-r:centre(ct,1)+r}))
  vec_box(ct) = max(speed_plane({centre(ct,0)-r:centre(ct,0)+r},\
			      {centre(ct,1)-r:centre(ct,1)+r}))

; Output file location and type
  output = "$sam/haima/images/track_new_"+dat+"_"+time_arr(it)
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
    opts_speed@tmXTOn                       = "False"    ; No tickmarks on top x-axis
    opts_speed@tmYROn                       = "False"    ; No tickmarks on right y-axis

    if (zoom .eq. 1) then
     opts_speed@mpMinLatF                    = centre(ct,0) - rad0
     opts_speed@mpMinLonF                    = centre(ct,1) - rad0
     opts_speed@mpMaxLatF                    = centre(ct,0) + rad0
     opts_speed@mpMaxLonF                    = centre(ct,1) + rad0
    else
     opts_speed@mpMinLatF                   = lat(0)
     opts_speed@mpMinLonF                   = lon(0)
     opts_speed@mpMaxLatF                   = lat(dimsizes(lat)-1)
     opts_speed@mpMaxLonF                   = lon(dimsizes(lon)-1)
    end if

    if (rad0 .le. 2.0) then 
     opts_speed@tmXBMode	            = "Explicit"
     opts_speed@tmXBValues	            = fspan(100.0, 150.0, 51)
     xb_labels0				    = sprintf("%0.1f~S~o~N~E",opts_speed@tmXBValues)
     opts_speed@tmXBLabels		    = xb_labels0

     opts_speed@tmYLMode                     = "Explicit"
     opts_speed@tmYLValues                   = fspan(0.0, 30.0, 31)
     yl_labels0                              = sprintf("%0.1f~S~o~N~N",opts_speed@tmYLValues)
     opts_speed@tmYLLabels                   = yl_labels0

     opts_speed@mpGridAndLimbOn              = True       ; Grid lines on/off
     opts_speed@mpGridLineDashPattern        = 2          ; Dash pattern      
     opts_speed@mpGridLatSpacingF            = 1.0        ; Spacing (latitude)
     opts_speed@mpGridLonSpacingF            = 1.0        ; Spacing (longitude)

    elseif (rad0 .gt. 2.0) then 

     opts_speed@tmXBMode                     = "Explicit"
     opts_speed@tmXBValues                   = fspan(100.0, 150.0, 26)
     xb_labels0                              = sprintf("%0.1f~S~o~N~E",opts_speed@tmXBValues)
     opts_speed@tmXBLabels                   = xb_labels0 

     opts_speed@tmYLMode                     = "Explicit"
     opts_speed@tmYLValues                   = fspan(0.0, 30.0, 16)
     yl_labels0                              = sprintf("%0.1f~S~o~N~N",opts_speed@tmYLValues)
     opts_speed@tmYLLabels                   = yl_labels0

     opts_speed@mpGridAndLimbOn              = True       ; Grid lines on/off
     opts_speed@mpGridLineDashPattern        = 2          ; Dash pattern
     opts_speed@mpGridLatSpacingF            = 2.0        ; Spacing (latitude)
     opts_speed@mpGridLonSpacingF            = 2.0        ; Spacing (longitude)

    end if 

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
    opts_mslp@gsnMaximize                   = True
    opts_mslp@pmTickMarkDisplayMode         = "Always"   ; Turn on map tickmarks
    opts_mslp@tmXBMajorLengthF              = 0.005      ; Change tickmark length
    opts_mslp@tmXTOn                        = "False"    ; No tickmarks on top x-axis
    opts_mslp@tmYROn                        = "False"    ; No tickmarks on right y-axis

;=================
; Plot the data
;=================

    spd_plot  = gsn_csm_contour_map(wks,speed_plane,opts_speed) ; 10-m wind speed           
    slp_plot  = gsn_csm_contour(wks,slp_plane,opts_mslp)        ; Mean sea level pressure   
    overlay(spd_plot,slp_plot)
    panel(0)  = spd_plot

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

    marker    = gsn_add_polymarker(wks,spd_plot,centre(ct,1),centre(ct,0),mres0)

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    optsP                       = True
    optsP@gsnFrame              = False             ; Do not advance the frame
    optsP@gsnPanelLabelBar      = False             ; Turn off panel labelbar
    optsP@txString              = "Haima ("+dat+"): (T+"+ctt+")"
    optsP@gsnPanelFigureStrings = ""
    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.90
    optsP@gsnPanelBottom        = 0.10
    optsP@amJust                = "TopLeft"
    optsP@gsnPanelFigureStringsFontHeightF = 0.0125 ; Label size (default 0.01)

   gsn_panel(wks,panel,(/1,1/),optsP)              ; Draw as a single plot
   frame(wks)

 ; Tidy up
   delete([/slp_plane,u_plane,v_plane,speed_plane/])
   delete([/utc_date,year,month,day,hour,minute,second/])

   ct = ct + 1 ; Counter variable (time)

  end do     ; End time loop (do it = 0, times-1, int)

; Tidy up before next iteration
  delete([/time_a,time_b,u,v,mslp/])

 end do       ; End input file loop (do nf = 0, numFILES-1)  

;===================================
; Write out values to text files
;===================================

  diri_out = "$ar/text/cp_hai_"

; MSLP
  slp_out = diri_out+dat+"_slp.txt"
  asciiwrite(slp_out,slp_box)

; Windspeed
  vec_out  = diri_out+dat+"_10m.txt"
  asciiwrite(vec_out,vec_box)

; Latitude 
  lat_out  = diri_out+dat+"_lat.txt"
  asciiwrite(lat_out,centre(:,0))

; Longitude
  lon_out  = diri_out+dat+"_lon.txt"
  asciiwrite(lon_out,centre(:,1))

end 