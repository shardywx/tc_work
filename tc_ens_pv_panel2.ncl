; Script to plot PV and vertical velocity in the TC inner core
; Run using:

; ncl 'ens="em08"' dist=4.0 p_cp=2 'opt="x11"' ts=12 tf=18 r0=0.75 tc_ens_pv_panel2.ncl

; 'ens'	  = ensemble simulation (em00-em11)
; 'opt'   = output file format ("pdf" or "x11")
; 'dist'  = size of box following storm (degrees)
; 'p_cp'  = pressure level for plotting PV and vertical velocity (see L106-108)
; 'ts'	  = starting time (0-120) 
; 'tf'	  = ending time (0-120)
; 'r0'	  = radius for final plot (degrees)
 

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"

begin

;================================================
; Read in cyclone position data from text file
;================================================

   diri        = "/nfs/a37/earshar/ncl_scripts/nepartak/archer_march2018/"
   numTIMES    = 120
   centre      = new((/numTIMES,2/),float) ; Array to hold all [120] lat/lon values

 ; Read in values from chosen simulation
   centre(:,0)  = asciiread(diri+"cp_"+ens+"_lat.txt",(/numTIMES,1/),"float")
   centre(:,1)  = asciiread(diri+"cp_"+ens+"_lon.txt",(/numTIMES,1/),"float")

;================================================================================
; Create lat/lon arrays containing grid subsets at each time, following the TC
;================================================================================
  lat_arr = new((/numTIMES,2/),float)
  lon_arr = new((/numTIMES,2/),float)

; Lat/lon arrays for 'em02'
  lat_arr(:,0) = centre(:,0)-dist
  lat_arr(:,1) = centre(:,0)+dist
  lon_arr(:,0) = centre(:,1)-dist
  lon_arr(:,1) = centre(:,1)+dist

;==================================================================================
;;;;;;;;;;;;;;;;;;;;;;;;;;; Start of main part of script ;;;;;;;;;;;;;;;;;;;;;;;;;;
;==================================================================================

; Prelash (define arrays for loops below)
  ens_arr     = ens
  numTIMES    = 120

  time_arr    = new(numTIMES,string)			    ; Time array
  title_arr   = new(numTIMES,string)			    ; Title array
  llbox       = toint(dist*50)            		    ; Calculate box size (user)
  dsize       = (/llbox,llbox/)           		    ; Box size
  pv_arr      = new((/numTIMES,dsize(0),dsize(1)/),"float") ; PV array
  str_arr     = new((/numTIMES,dsize(0),dsize(1)/),"float") ; Strain rate array
  geo_arr     = new((/numTIMES,dsize(0),dsize(1)/),"float") ; Height array
  vv_arr      = new((/numTIMES,dsize(0),dsize(1)/),"float") ; Omega array
  vel_arr     = new((/numTIMES,dsize(0),dsize(1)/),"float") ; Vector wind array

  date_str    = new(numTIMES,string)
  time_str    = new(numTIMES,string)
  out_str     = new(numTIMES,string)

  print_clock("Working on PV calculations!")

;====================================
; Start loop over ensemble members
;====================================

; List all files to read in and analyse
  diri        = "$ens/reg_netcdf/20160704T0000Z_ra1t_"
  input_pb    = diri+ens_arr+"_pb.nc"

  f  = addfile(input_pb,"r")
  print("Working on ensemble simulation: "+ens_arr)

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 1000000000
  end setvalues

;==================================
; Get the variables we will need
;==================================

  do it = ts, tf   		   ; Loop over subset of times

   time  = f->t(it)        	   ; Time in file
   times = dimsizes(time)	   

 ; Get domain boundaries at each time from previously-calcuated storm track ('centre')
   t0    = lat_arr(it,0)
   t1    = lat_arr(it,1)
   n0    = lon_arr(it,0)
   n1    = lon_arr(it,1)

 ; Create latitude/longitude array subsets to reduce computation time
   lon1  = f->longitude_1({n0:n1})  ; longitude [160 grid points]
   lat1  = f->latitude_1({t0:t1})   ; latitude  [160 grid points]
   lon   = f->longitude({n0:n1})    ; longitude [160 grid points]
   lat   = f->latitude({t0:t1})     ; latitude  [160 grid points]

 ; Pressure levels (pres)
 ; (0) 1000, (1) 950, (2) 900, (3) 850, (4) 800, (5) 750, (6) 700, (7) 650
 ; (8) 600, (9) 550, (10) 500, (11) 450, (12) 400, (13) 350
 ; (14) 300, (15) 250, (16) 200, (17) 150, (18) 100

   pres  = f->p(p_cp)               ; Subset of pressure level array
   plevs = dimsizes(pres)           ; Size of pressure level array
   pres@units = "hPa"		    ; Pressure level units (for plotting)

 ; Arrays have () vertical levels, and horizontal extent depends on 'dist' option
   u     = f->u(it,p_cp,{t0:t1},{n0:n1})     ; Zonal velocity (m/s)
   v     = f->v(it,p_cp,{t0:t1},{n0:n1})     ; Meridional velocity (m/s)
   vort  = f->rvor(it,p_cp,{t0:t1},{n0:n1})  ; Relative vorticity (/s)
   div   = f->div(it,p_cp,{t0:t1},{n0:n1})   ; Divergence (/s)
   pv    = f->pv(it,p_cp,{t0:t1},{n0:n1})    ; PV (PVU)

   t     = f->temp(it,p_cp,{t0:t1},{n0:n1})  ; Temperature (K)
   z     = f->ht(it,p_cp,{t0:t1},{n0:n1})    ; Geopotential height (m)
   w     = f->dz_dt(it,p_cp,{t0:t1},{n0:n1}) ; Vertical velocity (m/s)

   vel	 = sqrt((u^2) + (v^2))		     ; Vector velocity (m/s)
   copy_VarCoords(v,vel)
   copy_VarMeta(v,vel)

   w_size = dimsizes(w)
   ws1	  = w_size(0)-1
   ws2	  = w_size(1)-1

 ; Convert to units we want
   z     = z/10	            ; Convert to decametres
   z@units = "dam"	    ; Add metadata
   vort  = vort * (10 ^ -6) ; Convert to /s
   div   = div  * (10 ^ -6) ; Convert to /s
   pv    = pv * (10 ^ 6)    ; Convert to PVU [10 ^ 6 K m2 s-1 kg-1]

   dx    = lon(1) - lon(0) ; Grid spacing (longitude)
   dy    = lat(1) - lat(0) ; Grid spacing (latitude)

 ; Coriolis parameter and absolute vorticity
   f1   = coriolis_param(lat)
   f0   = conform_dims(dimsizes(vort),f1,1)
   copy_VarMeta(vort,f0)
   copy_VarCoords(vort,f0)

   avo  = vort + f0
   copy_VarMeta(vort,avo)
   copy_VarAtts(vort,avo)
   copy_VarCoords(vort,avo)
   avo@description = "Absolute vorticity"
   avo@long_name   = "Vorticity (absolute)"

;=====================================================
; Create correct date strings for each output time
;=====================================================

    month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
                      "Jul","Aug","Sep","Oct","Nov","Dec"/)
 
    times = dimsizes(time) ; Files are not all same size
    utc_date = cd_calendar(time, 0)

    year   = tointeger(utc_date(:,0))
    month  = tointeger(utc_date(:,1))
    day    = tointeger(utc_date(:,2))
    hour   = tointeger(utc_date(:,3))
    minute = tointeger(utc_date(:,4))
    second = utc_date(:,5)

  ; Correct for errors in the code (round up value of hour when minutes = 59)
    if (.not.ismissing(minute).and.minute.gt.30) then
      hour = hour+1
    end if  

    date_str(it) = sprinti("%0.2iUTC ", hour) + \
                   sprinti("%0.2i ", day) \
                   + month_abbr(month)
    time_str(it) = sprinti("%0.2i UTC ", hour) + \
                   sprinti("%0.2i ", day) \
                   + month_abbr(month)
    out_str(it)  = sprinti("%0.2i", day) + \
                   month_abbr(month) + \
                   "_" + sprinti("%0.2iZ", hour)

    print_clock("Working on time: "+time_str(it))
    print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")
    time_arr(it)  = out_str(it)
    itt = it + 1
    title_arr(it) = "Valid at "+time_str(it)+" (T+"+itt+")"

;====================================================================
; Find cyclone centre using both geopotential height and vorticity
;====================================================================

  ; TC latitude  -- centre(it,0)
  ; TC longitude -- centre(it,1)

;===========================================================
; Calculate rate of strain (following Nguyen et al. 2011)
;===========================================================

     dim1 = "latitude"
     dim2 = "longitude"

   ; Calculate horizontal derivatives of 'u' (du_dx, du_dy)
     du    = grad_latlon_cfd(u,u&$dim1$,u&$dim2$,False,False)
     du_dy = du[0]    ; Latitudinal derivative
     du_dx = du[1]    ; Longitudinal derivative
     delete(du)

   ; Calculate horizontal derivatives of 'v' (dv_dx, dv_dy)
     dv    = grad_latlon_cfd(v,v&$dim1$,v&$dim2$,False,False)
     dv_dy = dv[0]    ; Latitudinal derivative
     dv_dx = dv[1]    ; Longitudinal derivative
     delete(dv)

   ; Calculate strain --> sqrt[ (du_dx - dv_dy)^2 + (dv_dx + du_dy)^2 ]
     str1a   = du_dx - dv_dy   ; Stretching deformation
     str2a   = dv_dx + du_dy   ; Shearing deformation
     str1    = (str1a)^2
     str2    = (str2a)^2
     str_plane0  = sqrt( str1(:,:) + str2(:,:) )

;===================================================================
; Also read in lower-tropospheric PV [proxy for convective cells]
;===================================================================

     d0 = dsize(0)-1
     d1 = dsize(1)-1

   ; Read in PV, geopotential height and vertical velocity
     pv_plane   = pv(0:d0,0:d1)
     geo_plane  = z(0:d0,0:d1)
     w_plane	= w(0:d0,0:d1)
     str_plane  = str_plane0(0:d0,0:d1)
     copy_VarCoords(u,str_plane)

   ; Calculate max, min and mean geopotential height (for contours)
     zmax0	 = toint(max(geo_plane))
     zmin0	 = toint(min(geo_plane))
     zavg0       = ( (zmax0 + zmin0) / 2)

     zmax1	 = toint(max(geo_plane))
     zmin1	 = toint(min(geo_plane))
     zavg1	 = ( (zmax1 + zmin1) / 2)

     panel 	 = new(2,graphic)	; Panel plot with xxx images

;=================================
; Output file type and location
;=================================

     output = "$sam/nepartak/images/cart_coords/"+ens+"_str_pv"+pres+"_"+time_arr(it)
     wks = gsn_open_wks(opt,output)

   ; Load colour table
     gsn_define_colormap(wks,"prcp_new")

;==========================
; Options for plotting
;==========================

   ; Potential vorticity
     opts_pv                               = True
     opts_pv@cnFillOn                      = True
     opts_pv@cnLineLabelInterval           = 2.0
     opts_pv@cnLineLabelFontHeightF        = 0.012
     opts_pv@cnLineLabelBackgroundColor    = "transparent"
     opts_pv@cnLineLabelPlacementMode      = "constant"
     opts_pv@cnLinesOn                     = False ; Contour lines off
     opts_pv@cnInfoLabelOn                 = False ; Contour labels off
     opts_pv@cnLevelSelectionMode          = "ExplicitLevels"

     opts_pv@cnLevels                      = (/9.0, 10.0 ,12.0, 14.0, \
                                               16.0, 18.0, 20.0, 22.0, 25.0, \
                                               30.0, 35.0, 40.0, 45.0, 50.0/)
     opts_pv@cnFillColors                  = (/0,2,3,4,5,\
                                               6,7,8,9,10,\
                                               11,12,13,14,15/)

     opts_pv@gsnPaperOrientation           = "landscape"
     opts_pv@tiMainString                  = ""
     opts_pv@tiMainFontHeightF             = 0.0125
     opts_pv@gsnLeftString                 = ""
     opts_pv@gsnRightString                = ""
     opts_pv@lbLabelBarOn		   = True       ; Individual lb on/off

     opts_pv@pmLabelBarWidthF      	   = 0.52
     opts_pv@pmLabelBarHeightF     	   = 0.13
     opts_pv@lbLabelFontHeightF    	   = 0.0125
     opts_pv@lbPerimOn             	   = False
     opts_pv@lbLabelFont           	   = "Helvetica"
     opts_pv@lbTitleString         	   = pres+" hPa potential vorticity (PVU)"
     opts_pv@lbBoxEndCapStyle              = "TriangleBothEnds"    ; Labelbar end shape
     opts_pv@lbTitleFontHeightF    	   = 0.02
     opts_pv@lbTitleDirection      	   = "Across"
     opts_pv@lbTitlePosition       	   = "Bottom"

   ; Additional plotting resources
     opts_pv@mpDataBaseVersion             = "Ncarg4_1" ; More recent database
     opts_pv@mpDataSetName                 = "Earth..4" ; High resolution
     opts_pv@mpOutlineBoundarySets         = "National" ; National borders
     opts_pv@mpGeophysicalLineColor        = "black"    ; Colour borders black
     opts_pv@mpGeophysicalLineThicknessF   = 1.0        ; Border line thickness
     opts_pv@mpGridAndLimbOn               = False      ; Turn on lat/lon lines
     opts_pv@pmTickMarkDisplayMode         = "Always"   ; Turn on map tickmarks
     opts_pv@tmXBMajorLengthF              = 0.005      ; Change tickmark length
     opts_pv@tmXTOn                        = "False"    ; No tickmarks on top x-axis
     opts_pv@tmYROn                        = "False"    ; No tickmarks on right y-axis
     opts_pv@gsnMaximize                   = True       ; Maximise plot size 

     r					   = r0		; Radius of plot (user input)

     opts_pv@gsnAddCyclic                  = False
     opts_pv@mpLimitMode		   = "Corners"
     opts_pv@mpLeftCornerLatF              = centre(it,0)-r
     opts_pv@mpLeftCornerLonF              = centre(it,1)-r
     opts_pv@mpRightCornerLatF             = centre(it,0)+r
     opts_pv@mpRightCornerLonF             = centre(it,1)+r
     opts_pv@gsnDraw		           = False      ; Do not draw the plot
     opts_pv@gsnFrame		           = False      ; Do not advance the frame

   ; Vertical velocity
     opts_vv                               = True
     opts_vv@cnFillOn                      = False
     opts_vv@cnLineColor                   = "black"
     opts_vv@cnInfoLabelOn                 = False
     opts_vv@cnLineLabelsOn                = True
     opts_vv@cnLineLabelInterval           = 2.0
     opts_vv@cnLevelSelectionMode          = "ExplicitLevels"
     opts_vv@cnLevels                      = (/1.0,2.0/)
     opts_vv@cnLineLabelPlacementMode      = "constant"
     opts_vv@cnLineLabelPerimOn            = False
     opts_vv@gsnContourLineThicknessesScale = 3.0
     opts_vv@tiMainString                  = ""
     opts_vv@gsnLeftString                 = ""
     opts_vv@gsnRightString                = ""
     opts_vv@gsnDraw                       = False      ; Do not draw the plot
     opts_vv@gsnFrame                      = False      ; Do not advance the frame
     opts_vv@pmTickMarkDisplayMode         = "Always"   ; Turn on map tickmarks
     opts_vv@tmXBMajorLengthF              = 0.005      ; Change tickmark length
     opts_vv@tmXTOn                        = "False"    ; No tickmarks on top x-axis
     opts_vv@tmYROn                        = "False"    ; No tickmarks on right y-axis
     opts_vv@gsnMaximize                   = True       ; Maximise plot size

   ; Geopotential height
     opts_geo				   = opts_vv
     delete(opts_geo@cnLevels)
     opts_geo@cnLevels		           = ispan(toint(zmin0),toint(zmax0),4)
     opts_geo@gsnContourLineThicknessesScale = 3.0

;=========================================================
; Now plot rate of strain using different colour table
;=========================================================

   ; 850 hPa rate of strain
     opts_str                              = opts_pv
     opts_str@cnFillOn                     = True
     opts_str@cnLineLabelInterval          = 2.0
     opts_str@cnLineLabelFontHeightF       = 0.012
     opts_str@cnLineLabelBackgroundColor   = "transparent"
     opts_str@cnLineLabelPlacementMode     = "constant"
     opts_str@cnLinesOn                    = False ; Contour lines off
     opts_str@cnInfoLabelOn                = False ; Contour labels off
     opts_str@cnFillPalette		   = "circ_budget"
     opts_str@cnLevelSelectionMode         = "ExplicitLevels"
     delete(opts_str@cnLevels)
     delete(opts_str@cnFillColors)
     opts_str@cnLevels                     = (/10,12,14,16,18,20,22,25,30/)
     opts_str@cnFillColors                 = (/-1,11,13,15,16,17,18,19,20,21/)
     opts_str@gsnPaperOrientation          = "landscape"
     opts_str@tiMainString                 = ""
     opts_str@tiMainFontHeightF            = 0.0125
     opts_str@gsnLeftString                = ""
     opts_str@gsnRightString               = ""
     opts_str@lbLabelBarOn                 = True	      ; Individual labelbar 

     opts_str@pmLabelBarWidthF             = 0.52
     opts_str@pmLabelBarHeightF            = 0.13
     opts_str@lbLabelFontHeightF           = 0.0125
     opts_str@lbPerimOn                    = False
     opts_str@lbLabelFont                  = "Helvetica"
     opts_str@lbTitleString = pres+" hPa strain rate (10~S~-4~N~ s~S~-2~N~)"
     opts_str@lbBoxEndCapStyle             = "TriangleBothEnds"    ; Labelbar end shape
     opts_str@lbTitleFontHeightF           = 0.02
     opts_str@lbTitleDirection             = "Across"
     opts_str@lbTitlePosition              = "Bottom"

;=============
; Plot data 
;=============

     str_plane   = str_plane * (10 ^ 4)
     str_plane@units = "10~S~-4~N~ s~S~-2~N~"

     plot_pv	 = gsn_csm_contour_map(wks,pv_plane,opts_pv)  ; Potential vorticity
     plot_vv     = gsn_csm_contour(wks,w_plane,opts_vv)       ; Vertical velocity
     overlay(plot_pv, plot_vv)
     panel(0)   = plot_pv

     plot_str	= gsn_csm_contour_map(wks,str_plane,opts_str) ; Strain rate 
     panel(1)	= plot_str

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    optsP                       = True
    optsP@gsnFrame              = False    ; Do not advance the frame
    optsP@gsnPanelLabelBar      = False    ; Turn off panel labelbar
    optsP@txString              = title_arr(it)
    optsP@gsnPanelFigureStrings = ""
    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.90
    optsP@gsnPanelBottom        = 0.10
    optsP@amJust                = "TopLeft"
    optsP@gsnPanelFigureStringsFontHeightF = 0.0125 ; Label size (default 0.01)

    gsn_panel(wks,panel,(/1,2/),optsP)              ; Draw as a single plot
    frame(wks)

  ; Tidy up before next iteration
    delete([/pv_plane, geo_plane, w_plane, str_plane/])

  ; Tidy up (avoid dimension size errors)
    delete([/lat1,lon1,lat,lon,time,times,w,pv,u,v,vort,div,z,t,vel/])
    delete([/du_dx,dv_dy,dv_dx,du_dy/])
    delete([/utc_date,year,month,day,hour,minute,second/])

  end do      ; End time loop (do it = 0, times-1)

end 