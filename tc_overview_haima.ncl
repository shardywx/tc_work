; Script to calculate the momentum budget as in Persing et al. (2013)

; Run using:

; ncl opt=\"x11\" dat_str=\"15T00\" dist=3.0 tc_overview_haima.ncl

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

    print_clock("Working on initialisation time: "+dat_str)

  ; 00Z 4 Jul [00], 06Z 4 Jul [01], 12Z 4 Jul [02], 18Z 4 Jul [03], 00Z 5 Jul [04],
  ; 06Z 5 Jul [05], 12Z 5 Jul [06], 18Z 5 Jul [07], 00Z 6 Jul [08], 06Z 6 Jul [09],
  ; 12Z 6 Jul [10], 18Z 6 Jul [11], 00Z 7 Jul [12], 06Z 7 Jul [13], 12Z 7 Jul [14],
  ; 18Z 7 Jul [15], 00Z 8 Jul [16], 06Z 8 Jul [17], 12Z 8 Jul [18], 18Z 8 Jul [19],
  ; 00Z 9 Jul [20]

  ; ct = 0   [00], ct = 5  [01], ct = 11  [02], ct = 17  [03], ct = 23  [04],
  ; ct = 29  [05], ct = 35 [06], ct = 41  [07], ct = 47  [08], ct = 53  [09],
  ; ct = 59  [10], ct = 65 [11], ct = 71  [12], ct = 77  [13], ct = 83  [14],
  ; ct = 89  [15], ct = 95 [16], ct = 101 [17], ct = 107 [18], ct = 113 [19],
  ; ct = 119 [20]

  ; List all files ('pa' and 'pb' streams) to read in and analyse
    diri        = "$sam/um/cp/haima/201610"+dat_str+"00Z_HAIMA_4p4_L80_ra1t_"
    fili_pa     = diri+"pvera"
    fili_pb	= diri+"pverb"
    input_pa	= systemfunc("ls "+fili_pa+"*.nc")
    input_pb	= systemfunc("ls "+fili_pb+"*.nc")

    numINPUT0   = dimsizes(input_pb)        ; Number of input files
    numTIMES	= (numINPUT0 * 6)	    ; Number of input times 
    ct          = 0                         ; Counter variable
    ct0		= 0			    ; Counter variable (independent of position)
;    llbox	= toint(dist*50)	    ; Size of lat/lon arrays (see line below)
;    dsize	= (/llbox,llbox/)           ; Size of lat/lon arrays (depends on 'dist')

    setvalues NhlGetWorkspaceObjectId
      "wsMaximumSize" : 1000000000
    end setvalues

    date_str  = new(numTIMES,string)
    time_str  = new(numTIMES,string)
    out_str   = new(numTIMES,string)

    time_arr  = new((/numTIMES/),string)
    title_arr = new((/numTIMES/),string)

;===============================
; Start loop over input files 
;===============================

  do nf = 0,  numINPUT0-1		; Read in files one by one 
  a = addfile(input_pa(nf),"r")         ; Read in 'pa' stream [10-m wind]
  b = addfile(input_pb(nf),"r")         ; Read in 'pb' stream [PMSL, 1-km reflectivity]

  na = getfilevarnames(a)
  nb = getfilevarnames(b)

  print(na)
  print(nb)


  tm0 = a->t
  tm1	= dimsizes(tm0)

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

   date_str(ct0) = sprinti("%0.2iUTC ", hour) + \
                   sprinti("%0.2i ", day) \
                   + month_abbr(month)
   time_str(ct0) = sprinti("%0.2i UTC ", hour) + \
                   sprinti("%0.2i ", day) \
                   + month_abbr(month)
   out_str(ct0)  = sprinti("%0.2i", day) + \
                   month_abbr(month) + \
                   "_" + sprinti("%0.2iZ", hour)

   time_arr(ct)  = out_str(ct0)

   print("Working on time: "+time_str(ct0))
   title_arr(ct) = "Valid at "+time_str(ct0)+" (T+"+hour+")"

;=================================================
; Continue reading in variables from input file
;=================================================

;    t0    = lat_arr(ct,0)
;    t1    = lat_arr(ct,1)
;    n0    = lon_arr(ct,0)
;    n1    = lon_arr(ct,1)
;    print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

;    lon1  = a->longitude_1({n0:n1})  ; longitude [xxx grid points]
;    lat1  = a->latitude_1({t0:t1})   ; latitude  [xxx grid points]
;    lon   = a->longitude({n0:n1})    ; longitude [xxx grid points]
;    lat   = a->latitude({t0:t1})     ; latitude  [xxx grid points]

;===============
; 'pa' stream
;===============

  ; Read in latitude/longitude
    lon1  = a->longitude_1  	      ; longitude [950 grid points]
    lat1  = a->latitude_1     	      ; latitude  [551 grid points]
    lon   = a->longitude      	      ; longitude [950 grid points]
    lat   = a->latitude               ; latitude  [550 grid points]

  ; Read in horizontal velocity
    u_varname = "x-wind"
    v_varname = "y-wind"
    u     = a->$u_varname$(it,0,:,:)  ; 10-m zonal wind (m/s)
    v     = a->$v_varname$(it,0,:,:)  ; 10-m meridional wind (m/s)

    dx    = lon(1) - lon(0) 	      ; Grid spacing (longitude)
    dy    = lat(1) - lat(0) 	      ; Grid spacing (latitude)

;===============
; 'pb' stream
;===============

  ; Read in mean sea level pressure and 1-km radar reflectivity
    mslp  = b->p(it,0,:,:)	      ; Mean sea level pressure
    dbz	  = b->unspecified(it,0,:,:)  ; Radar reflectivity at 1 km

    mslp  = mslp / 100		      ; Convert from Pa to hPa
    mslp@units = "hPa"

  ; Calculate 10-m wind speed
    speed = sqrt( (u ^ 2.0) + (v ^ 2.0) ) 
    copy_VarCoords(u, speed)
    copy_VarAtts(u, speed)
    copy_VarMeta(u, speed)
    speed@units = "m s~S~-1~N~"

;===============================
; Plot data over whole domain
;===============================

  ; Output file location and type
    output = "$sam/haima/images/pmsl_uv_"+dat_str+"_"+time_arr(ct)
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

;    opts_prcp@mpLimitMode		   = "Corners"
;    opts_prcp@mpLeftCornerLatF		   = lat(0)
;    opts_prcp@mpLeftCornerLonF		   = lon(0)
;    opts_prcp@mpRightCornerLatF	   = lat(dimsizes(lat)-1)
;    opts_prcp@mpRightCornerLonF           = lon(dimsizes(lon)-1)

    opts_prcp@mpLimitMode                  = "LatLon"
    opts_prcp@mpMinLatF                    = lat(0)
    opts_prcp@mpMinLonF                    = lon(0)
    opts_prcp@mpMaxLatF                    = lat(dimsizes(lat)-1)
    opts_prcp@mpMaxLonF                    = lon(dimsizes(lon)-1)

    opts_prcp@tmXTOn                       = "False"    ; No tickmarks on top x-axis
    opts_prcp@tmYROn                       = "False"    ; No tickmarks on right y-axis

  ; Wind speed resources
    opts_speed                             = opts_prcp
    delete(opts_speed@cnLevels)
    delete(opts_speed@cnFillColors)
    opts_speed@cnLevels                    = (/15.0, 20.0,\
                                               25.0, 30.0, 35.0, 40.0, \
                                               45.0, 50.0, 55.0, 60.0, \
                                               65.0/)
    opts_speed@cnFillColors                = (/0,4,5,6,7,8,9,10,11,12,13,14/)
;    opts_speed@mpLimitMode                 = "LatLon"
;    opts_speed@mpMinLatF                   = lat(0)
;    opts_speed@mpMinLonF                   = lon(0)
;    opts_speed@mpMaxLatF                   = lat(dimsizes(lat)-1)
;    opts_speed@mpMaxLonF                   = lon(dimsizes(lon)-1)

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
;    opts_mslp@mpLimitMode                  = "LatLon"
;    opts_mslp@mpMinLatF             	   = lat(0)
;    opts_mslp@mpMinLon			   = lon(0)
;    opts_mslp@mpMaxLatF            	   = lat(dimsizes(lat)-1)
;    opts_mslp@mpMaxLonF		   = lon(dimsizes(lon)-1)
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

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    optsP                       = True
    optsP@gsnFrame              = False             ; Do not advance the frame
    optsP@gsnPanelLabelBar      = False             ; Turn off panel labelbar
    optsP@txString              = "Haima (2016): "+title_arr(ct)
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

    ct  = ct + 1		    ; Counter variable (time)
    ct0 = ct0 + 1		    ; Independent counter variable (time)

   end do     ; End time loop (do it = 0, tm1-1)

   delete([/tm0,tm1/])

 end do     ; End input file loop (do nf = 0, numINPUT0-1)

end 