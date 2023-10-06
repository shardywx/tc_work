; Script to calculate and write out storm track and intensity for basic plotting 

; Run using:

; ncl 'opt="x11"' 'tc="mer"' 'dat_str="09T0000Z"' dist=3.0 tc_overview_total.ncl

; For HAIMA ("hai"); dat_str   = 15T0000Z --> 21T0000Z
; For MERANTI ("mer"); dat_str = 09T0000Z --> 14T0000Z

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"
load "$sam/ncl_func/st_rm.ncl"

begin

;==================================================================================
;;;;;;;;;;;;;;;;;;;;;;;;;;; Start of main part of script ;;;;;;;;;;;;;;;;;;;;;;;;;;
;==================================================================================

  print_clock("Starting main part of script...")

;=========================
; Start multiple loops
;=========================

  ; List all files ('pa' and 'pb' streams) to read in and analyse
    if (tc.eq."hai") then
     diri        = "$sam/um/cp/haima/201610"+dat_str+"_HAIMA_4p4_L80_ra1t_"
     print_clock("Working on Haima – initialisation time: "+dat_str)
    else if (tc.eq."mer") then
     diri        = "$sam/um/cp/meranti/201609"+dat_str+"_MRNTI_4p4_L80_ra1t_"
     print_clock("Working on Meranti – initialisation time: "+dat_str)
    end if
    end if 

    fili_pa     = diri+"pvera"
    fili_pb	= diri+"pverb"
;    fili_pc	= diri+"pverc"
    input_pa	= systemfunc("ls "+fili_pa+"*.nc")
    input_pb	= systemfunc("ls "+fili_pb+"*.nc")
;    input_pc	= systemfunc("ls "+fili_pc+"*.nc")

    numINPUT0   = dimsizes(input_pb)        ; Number of input files
    numTIMES	= (numINPUT0 * 6)	    ; Number of input times 
    ct       	= 6                         ; Counter variable
    ct0		= 0			    ; Counter variable (independent of position)
    lgc		= False			    ; Boolean operator (for use in date array)
;    llbox	= toint(dist*50)	    ; Size of lat/lon arrays (see line below)
;    dsize	= (/llbox,llbox/)           ; Size of lat/lon arrays (depends on 'dist')

    setvalues NhlGetWorkspaceObjectId
      "wsMaximumSize" : 1000000000
    end setvalues

    time_str  = new(numTIMES,string)	    ; 18 UTC 10 Sep
    dat_um    = new(numTIMES,string)	    ; 10Sep_18utc
  ; dat_bt0 (defined below)                 ; 10Sep_18utc (same	as 'dat_um')

    time_arr  = new(numTIMES,string)	    ; 10Sep_18Z
    title_arr = new(numTIMES,string)	    ; Valid at 18 UTC 10 Sep (T+6)

    slp_box   = new(numTIMES,float)	    ; Minimum MSLP values 
    vec_box   = new(numTIMES,float)	    ; Maximum windspeed values 
    lat_box   = new(numTIMES,float)	    ; Storm position (latitude)
    lon_box   = new(numTIMES,float)	    ; Storm position (longitude)

;=================================================
; Read in best track information from text file 
;=================================================

  ; Number of times in each set of arrays
    dsh = 30  	    ; Haima
    dsm = 28	    ; Meranti

    if (tc.eq."hai") then
     lat_bt0  = asciiread("$sam/um/cp/haima/lat_hai.txt",(/dsh/),"float")
     lon_bt0  = asciiread("$sam/um/cp/haima/lon_hai.txt",(/dsh/),"float")
     dat_bt0  = asciiread("$sam/um/cp/haima/dat_hai.txt",(/dsh/),"string")
    else if (tc.eq."mer") then
     lat_bt0  = asciiread("$sam/um/cp/meranti/lat_mer.txt",(/dsm/),"float")
     lon_bt0  = asciiread("$sam/um/cp/meranti/lon_mer.txt",(/dsm/),"float")
     dat_bt0  = asciiread("$sam/um/cp/meranti/dat_mer.txt",(/dsm/),"string")
    end if 
    end if 

;===========================================================
; Create enlarged arrays of hourly storm position (basic)
;===========================================================

  ; Number of elements in the expanded best-track arrays
    btTIMES   = dimsizes(lat_bt0) * 6

  ; Estimate hourly storm position from 6-h best track data
    lat_bt1   = new(btTIMES,float)
    lon_bt1   = new(btTIMES,float)
    lat_bt    = new(btTIMES,float)
    lon_bt    = new(btTIMES,float)
    dat_bt    = new(btTIMES,string)
    dlat      = new(btTIMES,float)
    dlon      = new(btTIMES,float)

    do il = 0, dimsizes(lat_bt0)-1

     il0 = il * 6

     lp1 = min((/il+1,dimsizes(lat_bt0)-1/))
     lm1 = max((/il-1,0/))

   ; Estimated 1-h change in storm position 
     dlat(il0) 	         = 0
     dlon(il0)	         = 0

     dlat(il0+1)         = 1 * ((lat_bt0(lp1) - lat_bt0(il)) / 6)
     dlon(il0+1)         = 1 * ((lon_bt0(lp1) - lon_bt0(il)) / 6)

     dlat(il0+2)         = 2 * ((lat_bt0(lp1) - lat_bt0(il)) / 6)
     dlon(il0+2)         = 2 * ((lon_bt0(lp1) - lon_bt0(il)) / 6)

     dlat(il0+3)         = 3 * ((lat_bt0(lp1) - lat_bt0(il)) / 6)
     dlon(il0+3)         = 3 * ((lon_bt0(lp1) - lon_bt0(il)) / 6)

     dlat(il0+4)         = 4 * ((lat_bt0(lp1) - lat_bt0(il)) / 6)
     dlon(il0+4)         = 4 * ((lon_bt0(lp1) - lon_bt0(il)) / 6)

     dlat(il0+5)         = 5 * ((lat_bt0(lp1) - lat_bt0(il)) / 6)
     dlon(il0+5)         = 5 * ((lon_bt0(lp1) - lon_bt0(il)) / 6)

     lat_bt(il0:il0+5)   = lat_bt0(il)
     lon_bt(il0:il0+5)   = lon_bt0(il)
     dat_bt(il0:il0+5)   = dat_bt0(il)

    end do

    lat_bt = lat_bt + dlat
    lon_bt = lon_bt + dlon

;===============================
; Start loop over input files 
;===============================

  do nf = 0,  numINPUT0-1		; Read in files one by one 
  a = addfile(input_pa(nf),"r")         ; Read in 'pa' stream [10-m wind]
 
  tm0 = a->t
  tm1 = dimsizes(tm0)

;==================================
; Get the variables we will need
;==================================

  do it = 0, tm1-1		; Loop over times in file (6) 

   time  = a->t(it)            ; Read in all times in file

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
   if (minute.gt.30) then
     hour = hour+1
   end if

   time_str(ct0) = sprinti("%0.2i UTC ", hour) + \
                   sprinti("%0.2i ", day) \
                   + month_abbr(month)
   time_arr(ct0) = sprinti("%0.2i", day) + \
                   month_abbr(month) + \
                   "_" + sprinti("%0.2iZ", hour)

   if (hour.gt.9) then 
    dat_um(ct0)	 = day+month_abbr(month)+"_"+hour+"utc"
   else if (hour.le.9) then
    dat_um(ct0)	 = day+month_abbr(month)+"_0"+hour+"utc"
   end if 
   end if

   ct1 = ct0 + 6
   title_arr(ct0) = "Valid at "+time_str(ct0)+" (T+"+ct1+")"

   ct0 = ct0 + 1 
   ct  = ct + 1

  end do	 ; End time loop (do it = 0, tm1-1)

 end do		 ; End input file loop (do nf = 0,  numINPUT0-1)

;====================================
; Start main loop over input files
;====================================

; Redefine counter variables 
  ct          = 6
  ct0         = 0                       ; Independent of position

; Set starting index for storm-following domain
 if (tc.eq."mer") then
  if (dat_str.eq."09T0000Z") then
   index      = 12
  elseif (dat_str.eq."09T1200Z") then 
   index      = 24
  elseif (dat_str.eq."10T0000Z") then
   index      = 36
  elseif (dat_str.eq."10T1200Z") then
   index      =	48
  elseif (dat_str.eq."11T0000Z") then
   index      =	60
  elseif (dat_str.eq."11T1200Z") then
   index      =	72
  elseif (dat_str.eq."12T0000Z") then
   index      =	84
  elseif (dat_str.eq."12T1200Z") then
   index      =	96
  elseif (dat_str.eq."13T0000Z") then
   index      =	108
  elseif (dat_str.eq."13T1200Z") then
   index      =	120
  elseif (dat_str.eq."14T0000Z") then
   index      =	132
  end if 

 elseif (tc.eq."hai") then 
  if (dat_str.eq."15T0000Z") then
   index      = 12
  elseif (dat_str.eq."15T1200Z") then
   index      = 24
  elseif (dat_str.eq."16T0000Z") then
   index      = 36
  elseif (dat_str.eq."16T1200Z") then
   index      = 48
  elseif (dat_str.eq."17T0000Z") then
   index      = 60
  elseif (dat_str.eq."17T1200Z") then
   index      = 72
  elseif (dat_str.eq."18T0000Z") then
   index      = 84
  elseif (dat_str.eq."18T1200Z") then
   index      = 96
  elseif (dat_str.eq."19T0000Z") then
   index      = 108
  elseif (dat_str.eq."19T1200Z") then
   index      = 120
  elseif (dat_str.eq."20T0000Z") then
   index      = 132
  elseif (dat_str.eq."20T1200Z") then
   index      = 144
  elseif (dat_str.eq."21T0000Z") then
   index      = 156
  end if
 end if

  do nf = 0,  numINPUT0-1               ; Read in files one by one
  a = addfile(input_pa(nf),"r")         ; Read in 'pa' stream [10-m wind]
  b = addfile(input_pb(nf),"r")         ; Read in 'pb' stream [PMSL, 1-km reflectivity]
;  c = addfile(input_pc(nf),"r")         ; Read in 'pc' stream [u,v,theta-w on p-levs]

  tm0 = a->t
  tm1 = dimsizes(tm0)
  FirstTime = True 

;==================================
; Get the variables we will need
;==================================

  do it = 0, tm1-1              ; Loop over times in file (6)

    time  = a->t(it)            ; Read in all times in file

    print("Working on time: "+ct0)

;=================================================
; Continue reading in variables from input file
;=================================================

    lat_storm = lat_bt(index)
    lon_storm = lon_bt(index)

    t0        = lat_storm - dist 
    t1        = lat_storm + dist 
    n0        = lon_storm - dist
    n1        = lon_storm + dist

    if ( tc .eq. "hai" .and. lon_storm .lt. (115.6) ) then
     print("lon_storm = "+lon_storm)
     print_clock("Cyclone has reached the edge of the domain. Exiting...")
     break
    end if

    if ( tc .eq. "mer" .and. lon_storm .lt. (115.6) ) then
     print("lon_storm =	     "+lon_storm)
     print_clock("Cyclone has reached the edge of the domain. Exiting...")
     break
    end if

    ct1       = ct0 + 6

    print("Best track storm position at "+time_str(ct0)+" (T+"+ct1+"): "+lat_storm+\
    	  " deg N, "+lon_storm+" deg E")

    lon1  = a->longitude_1({n0:n1})  ; longitude [xxx grid points]
    lat1  = a->latitude_1({t0:t1})   ; latitude  [xxx grid points]
    lon   = a->longitude({n0:n1})    ; longitude [xxx grid points]
    lat   = a->latitude({t0:t1})     ; latitude  [xxx grid points]

;===============
; 'pa' stream
;===============

  ; Read in horizontal wind 
    u_varname = "x-wind"
    v_varname = "y-wind"
    u     = a->$u_varname$(it,0,{t0:t1},{n0:n1})  ; 10-m zonal wind (m/s)
    v     = a->$v_varname$(it,0,{t0:t1},{n0:n1})  ; 10-m meridional wind (m/s)

    dx    = lon(1) - lon(0) 	      ; Grid spacing (longitude)
    dy    = lat(1) - lat(0) 	      ; Grid spacing (latitude)

  ; Calculate vector wind 
    speed = sqrt(u ^ 2.0 + v ^ 2.0)
    copy_VarCoords(u,speed)
    copy_VarAtts(u,speed)
    copy_VarMeta(u,speed)
    speed@units = "m s~S~-1~N~"

  ; Write maximum vector wind to array
    vec_box(ct0) = max(speed)

;===============
; 'pb' stream
;===============

  ; Read in mean sea level pressure and 1-km radar reflectivity
    mslp  = b->p(it,0,{t0:t1},{n0:n1})	          ; Mean sea level pressure
    dbz	  = b->unspecified(it,0,{t0:t1},{n0:n1})  ; Radar reflectivity at 1 km

    mslp  = mslp / 100		      		  ; Convert from Pa to hPa
    mslp@units = "hPa"

  ; Calculate 10-m wind speed
    speed = sqrt( (u ^ 2.0) + (v ^ 2.0) ) 
    copy_VarCoords(u, speed)
    copy_VarAtts(u, speed)
    copy_VarMeta(u, speed)
    speed@units = "m s~S~-1~N~"

  ; Write minimum MSLP to array and calculate storm centre
    slp_box(ct0) = min(mslp)
    dims_p0      = dimsizes(mslp)
    slp_1d0      = ndtooned(mslp)
    inds_p0      = ind_resolve(minind(slp_1d0),dims_p0)  ; Find index of minimum

    lat_slp      = lat(0) + (dy * inds_p0(0,0))
    lon_slp      = lon(0) + (dx * inds_p0(0,1))
    print("Storm centre: "+lat_slp+" deg N, "+lon_slp+" deg E")

  ; If difference exceeds a threshold, use the previous storm location as our base
    if (.not.FirstTime) then 

     dlat0 = abs(lat_slp - lat0)	; Difference between current and previous centre
     dlon0 = abs(lon_slp - lon0)	; Lat/lon at previous time = [lat0/lon0]

     if (dlat0 .gt. 1.0 .or. dlon0 .gt. 1.0) then

      print("Change in lat = "+dlat0+" ; change in lon = "+dlon0)
      inc = 0.75
      lt1 = lat0 - inc
      lt2 = lat0 + inc
      ln1 = lon0 - inc
      ln2 = lon0 + inc

      slp_plane0  = mslp({lt1:lt2},{ln1:ln2})
      dims_p00	  = dimsizes(mslp({lt1:lt2},{ln1:ln2}))
      slp_1d00	  = ndtooned(mslp({lt1:lt2},{ln1:ln2}))
      inds_p00	  = ind_resolve(minind(slp_1d00),dims_p00)

      lat_slp	  = lt1 + (dy * inds_p00(0,0))
      lon_slp	  = ln1 + (dy * inds_p00(0,1))
      print("New storm centre: "+lat_slp+" deg N, "+lon_slp+" deg E")

      lat_box(ct0) = lat_slp
      lon_box(ct0) = lon_slp

     else

      lat_box(ct0) = lat_slp
      lon_box(ct0) = lon_slp

     end if 

    else if (FirstTime) then 
    
      lat_box(ct0) = lat_slp
      lon_box(ct0) = lon_slp

    end if 
    end if

    if (tc.eq."hai") then 
     if (lon_box(ct0) .le. (115.6) ) then 
      print_clock("Cyclone has reached the edge of the domain. Exiting...")
      break
     end if 
    else if (tc.eq."mer") then 
     if (lon_box(ct0) .le. (117.8) ) then 
      print_clock("Cyclone has reached the edge of the domain. Exiting...")
      break
     end if 
    end if 
    end if 

    delete([/dims_p0,slp_1d0,inds_p0/])

;===============
; 'pc' stream
;===============

;  ; Read in 850 hPa relative vorticity 
;    p	= c->p
;    t	= c->t
;    u_p = c->u
;    v_p = c->v

;===============================
; Plot data over whole domain
;===============================

  ; Output file location and type
    if (tc.eq."hai") then
     output = "$sam/haima/images/track_"+dat_str+"_"+time_arr(ct0)
    else if (tc.eq."mer") then
     output = "$sam/meranti/images/track_"+dat_str+"_"+time_arr(ct0)
    end if
    end if 

    wks = gsn_open_wks(opt,output)

;==============================================================
; Load NIMROD colour table for precip duration data
;==============================================================

  ; Load colour table
    gsn_define_colormap(wks,"radar_new")

  ; Change the colour map around & make first colour transparent
    cmap = read_colormap_file("radar_new")
    cmap(0,3) = 0.0

  ; Create panel plot (MSLP/wind and reflectivity)
    panel = new(1,graphic)

;==========================
; Options for plotting
;==========================

  ; Rainfall rate
    opts_prcp                              = True
    opts_prcp@cnFillOn                     = True
    opts_prcp@cnLineLabelInterval          = 2.0
    opts_prcp@cnLineLabelFontHeightF       = 0.012
    opts_prcp@cnLineLabelBackgroundColor   = "transparent"
    opts_prcp@cnLineLabelPlacementMode     = "constant"
    opts_prcp@cnLinesOn                    = False ; Contour lines off
    opts_prcp@cnInfoLabelOn                = False ; Contour labels off
    opts_prcp@cnLevelSelectionMode         = "ExplicitLevels"
    opts_prcp@cnFillPalette                = cmap
    opts_prcp@cnLevels                     = (/1,2,5,10,20,30,40,50,60,80,100/)
    opts_prcp@cnFillColors                 = (/0,14,2,3,5,7,8,9,10,11,12,13/)

    opts_prcp@gsnPaperOrientation          = "landscape"
    opts_prcp@tiMainString                 = ""
    opts_prcp@tiMainFontHeightF            = 0.0125
    opts_prcp@gsnLeftString                = ""
    opts_prcp@gsnRightString               = ""

    opts_prcp@lbLabelBarOn                 = True
    opts_prcp@lbBoxEndCapStyle             = "TriangleBothEnds"    ; Labelbar end shape
    opts_prcp@lbLabelFontHeightF           = 0.0125                ; Labelbar font size
    opts_prcp@lbLabelFont                  = "Helvetica"           ; Labelbar font
    opts_prcp@lbTitleString		   = "Wind speed (m s~S~-1~N~)"
    opts_prcp@lbTitlePosition		   = "Bottom"	           ; Title position
    opts_prcp@lbTitleFontHeightF	   = 0.015		   ; Title font height
    opts_prcp@pmLabelBarOrthogonalPosF	   = 0.1		   ; Labelbar position
    opts_prcp@lbTitleOffsetF		   = 0			   ; Lb title up/down
    opts_prcp@lbPerimOn                    = False                 ; Perimeter on/off
    opts_prcp@gsnDraw                      = False                 ; Do not draw plot
    opts_prcp@gsnFrame                     = False                 ; Do not advance frame
    opts_prcp@gsnAddCyclic                 = False

    opts_prcp@mpLimitMode                  = "LatLon"
    opts_prcp@mpMinLatF                    = lat(0)
    opts_prcp@mpMinLonF                    = lon(0)
    opts_prcp@mpMaxLatF                    = lat(dimsizes(lat)-1)
    opts_prcp@mpMaxLonF                    = lon(dimsizes(lon)-1)

    opts_prcp@tmXTOn                       = "False"    ; No tickmarks on top x-axis
    opts_prcp@tmYROn                       = "False"    ; No tickmarks on right y-axis
    opts_prcp@tmXBMode			   = "Explicit"	; Tickmark labels

    opts_prcp@tmXBValues                   = fspan(116.0,148.0,9)  ; Labels every 4 deg
    opts_prcp@tmYLValues                   = fspan(0.0,30.0,16)    ; Label every 2 deg
    opts_prcp@tmXBLabels                   = (/"116.0","102.0","124.0",\
                                               "128.0","132.0","136.0",\
                                               "140.0","144.0","148.0"/)
    opts_prcp@tmYLLabels                   = (/"0.0","2.0","4.0","6.0","8.0","10.0",\
                                               "12.0","14.0","16.0","18.0","20.0",\
                                               "22.0","24.0","26.0","28.0","30.0"/)

    opts_prcp@mpGridAndLimbOn		   = True	; Grid lines on/off
    opts_prcp@mpGridLineDashPattern	   = 2		; Dash pattern
    opts_prcp@mpGridLatSpacingF		   = 2.0	; Spacing (latitude)
    opts_prcp@mpGridLonSpacingF		   = 2.0	; Spacing (longitude)
    
  ; Wind speed resources
    opts_speed                             = opts_prcp
    delete(opts_speed@cnLevels)
    delete(opts_speed@cnFillColors)
    opts_speed@cnLevels                    = (/15.0, 20.0,\
                                               25.0, 30.0, 35.0, 40.0, \
                                               45.0, 50.0, 55.0, 60.0, \
                                               65.0/)
    opts_speed@cnFillColors                = (/0,4,5,6,7,8,9,10,11,12,13,14/)

  ; Mean sea level pressure resources 
    opts_mslp				   = True
    opts_mslp@cnFillOn			   = False
    opts_mslp@cnLineColor		   = "black"
    opts_mslp@cnLevelSelectionMode         = "ExplicitLevels"
    opts_mslp@cnLevels			   = ispan(880,1020,4)
    opts_mslp@cnInfoLabelOn		   = False
    opts_mslp@cnLineThicknessF		   = 3.0
    opts_mslp@cnLineLabelsOn		   = False
    opts_mslp@gsnDraw			   = False
    opts_mslp@gsnFrame			   = False
    opts_mslp@gsnAddCyclic                 = False
    opts_mslp@gsnLeftString                = ""
    opts_mslp@gsnRightString               = ""

;=================
; Plot the data
;=================

;    prcp_plot = gsn_csm_contour_map(wks,dbz,opts_prcp)   ; Radar reflectivity (1-km)
;    panel(0)  = prcp_plot

    spd_plot  = gsn_csm_contour_map(wks,speed,opts_speed) ; 10-m wind speed
    slp_plot  = gsn_csm_contour(wks,mslp,opts_mslp)       ; Mean sea level pressure
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

    marker    = gsn_add_polymarker(wks,spd_plot,lon_box(ct0),lat_box(ct0),mres0)

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    optsP                       = True
    optsP@gsnFrame              = False             ; Do not advance the frame
    optsP@gsnPanelLabelBar      = False             ; Turn off panel labelbar
    if (tc.eq."hai") then
     optsP@txString             = "Haima (2016): "+title_arr(ct0)
    else if (tc.eq."mer") then
     optsP@txString             = "Meranti (2016): "+title_arr(ct0)    
    end if
    end if
    optsP@gsnPanelFigureStrings = ""
    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.90
    optsP@gsnPanelBottom        = 0.10
    optsP@amJust                = "TopLeft"
    optsP@gsnPanelFigureStringsFontHeightF = 0.0125 ; Label size (default 0.01)

    gsn_panel(wks,panel,(/1,1/),optsP)              ; Draw as a single plot
    frame(wks)

  ; Tidy up
    delete([/lat,lon,lat1,lon1,dx,dy,u,v,mslp,dbz,speed/])

  ; Current lat/lon values --> if dlat/dlon at any time is too large, exit the loop
    lat0 = lat_box(ct0)
    lon0 = lon_box(ct0)

    ct  = ct + 1		    ; Counter variable (time)
    ct0 = ct0 + 1		    ; Independent counter variable (time)

    index = index + 1		    ; Counter variable (best-track lat/lon array)

    FirstTime = False

   end do     ; End time loop (do it = 0, tm1-1)

   delete([/tm0,tm1/])

 end do     ; End input file loop (do nf = 0, numINPUT0-1)

;==================================
; Write out values to text files 
;==================================

  ; Minimum mean sea level pressure 
    if (tc.eq."hai") then 
     slp_out = "$sam/haima/"+dat_str+"_slp.txt"
     asciiwrite(slp_out,slp_box)
    else if (tc.eq."mer") then 
     slp_out = "$sam/meranti/"+dat_str+"_slp.txt"
     asciiwrite(slp_out,slp_box)
    end if 
    end if 

  ; Maximum vector wind 
    if (tc.eq."hai") then
     vec_out = "$sam/haima/"+dat_str+"_vec.txt"
     asciiwrite(vec_out,vec_box)
    else if (tc.eq."mer") then	
     vec_out = "$sam/meranti/"+dat_str+"_vec.txt"
     asciiwrite(vec_out,vec_box)
    end	if
    end	if

  ; Storm track (latitude / longitude)
    if (tc.eq."hai") then
     lat_out = "$sam/haima/"+dat_str+"_lat.txt"
     asciiwrite(lat_out,lat_box)
    else if (tc.eq."mer") then
     lat_out = "$sam/meranti/"+dat_str+"_lat.txt"
     asciiwrite(lat_out,lat_box)
    end if
    end if
    
    if (tc.eq."hai") then
     lon_out = "$sam/haima/"+dat_str+"_lon.txt"
     asciiwrite(lon_out,lon_box)
    else if (tc.eq."mer") then
     lon_out = "$sam/meranti/"+dat_str+"_lon.txt"
     asciiwrite(lon_out,lon_box)
    end if
    end if

end