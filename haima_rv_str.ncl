; Script to identify vacillation cycles (VCs) in UM ensemble simulations

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_rm.ncl"

; ncl dat=\"15T00\" opt=\"x11\" cp0=3 dist=3.0 r0=1.0 \
; tc=0 lay=1 cp1=6 haima_rv_str.ncl

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'ens0'  = ensemble simulation (em00-em11)
; 'dist'  = size of box following storm (degrees)
; 'opt'   = output file format ("pdf" or "x11")
; 'cp0'   = pressure level for plotting PV and vertical velocity (see L106-108)
; 'r0'    = radius for final plot (degrees)
; 'tc'	  = annotate images with TC symbol (1) or leave as they were (0)
; 'lay'	  = calculate layer average (1) 
; 'cp1'   = upper pressure level for calculating layer average (only valid when lay=1)

begin

;==============================
; Prelash (reading in files)
;==============================

  diri         = "/nfs/a37/earshar/um/cp/haima/"
  fili_prefix  = "201610"+dat+"00Z_HAIMA_4p4_L80_ra1t_"

; Find 'pc' stream data (u,v on pressure levels)
  fili_c1      = systemfunc("cd "+diri+" ; ls "+fili_prefix+"pverc???_new.nc")
  fili_c       = diri+"/"+fili_c1

  numFILES     = dimsizes(fili_c)                ; 16 input files (T+006 to T+096)        
  numTEXT      = 96                              ; No. of values in 1-h text files (96)
  numINPUT     = (numTEXT / 3) - 1               ; No. of values in 3-h arrays in this script (31)

  print("Working on initialisation time "+dat)

; Create arrays to hold storm centre positions
  centre       = new((/numINPUT,2/),float)
  lat_arr      = new((/numINPUT,2/),float)
  lon_arr      = new((/numINPUT,2/),float)

; Read in storm centre positions from text files 
  lat1         = asciiread("$ar/text/cp_hai_"+dat+"_lat.txt",(/numTEXT/),"float")
  lon1         = asciiread("$ar/text/cp_hai_"+dat+"_lon.txt",(/numTEXT/),"float")

; Select every 3rd value to match up 3-h time intervals (31 values)
  centre(:,0)  = lat1(3::3)
  centre(:,1)  = lon1(3::3)
  delete([/lat1, lon1/])
  
; Create domain subset arrays following the storm 
  lat_arr(:,0) = centre(:,0) - dist
  lat_arr(:,1) = centre(:,0) + dist
  lon_arr(:,0) = centre(:,1) - dist
  lon_arr(:,1) = centre(:,1) + dist

; Calculate domain size (subset)
  llbox	       = toint(dist*50)	
  dsize	       = (/llbox,llbox/)
  d0           = dsize(0)
  d1           = dsize(1)

; Counter variable (time)
  ct           = 0

; Create arrays to hold date strings                                                     
  time_str  = new(numINPUT,string)
  time_arr  = new(numINPUT,string)

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
  end if

;==================================
; Get the variables we will need 
;==================================

; Loop over input files 
  do nf        = 0, numFILES-1
  c            = addfile(fili_c(nf),"r")       ; Read in pressure-level data

  time_c       = c->t

 ; Loop over input times
   if (nf .eq. numFILES-1) then 
    tf = 0
   else 
    tf = 1
   end if 

    do it       = 0, tf

     ctt  = ct + 7
     print("Working on forecast hour: "+ctt)
     title = ctt+" h forecast, initialised at "+date_str

   ; Get domain boundaries from offline calculations 
     t0   = lat_arr(ct,0)
     t1   = lat_arr(ct,1)
     n0   = lon_arr(ct,0)
     n1   = lon_arr(ct,1)

   ; Read in all variables 
     lat  = c->latitude({t0:t1})
     lon  = c->longitude({n0:n1})
     lat1 = c->latitude_1({t0:t1})
     lon1 = c->longitude_1({n0:n1})

   ; Pressure level array
     pres = c->p
     print("Working on "+pres(cp0)+" hPa")

   ; EDIT TO ADD OPTIONS FOR LAYER CALCULATIONS

   ; Horizontal windspeed (remove record dimension we added with 'ncecat')
   ; This gives us windspeed on single pressure level (2D lat/lon array)
     u    = c->u(0,it,cp0,{t0:t1},{n0:n1})
     v    = c->v(0,it,cp0,{t0:t1},{n0:n1})

;=============================================
; Calculate storm-relative winds (add later)
;=============================================

;    u = u0 - u_cyc(it-1)
;    v = v0 - v_cyc(it-1)

;  ; Add metadata from original horizontal wind arrays
;    copy_VarCoords(u0,u)
;    copy_VarAtts(u0,u)

;    copy_VarCoords(v0,v)
;    copy_VarAtts(v0,v)

;==================================================================
; Calculate relative vorticity using centered finite differences
;==================================================================

  ; Option '2' --> boundary points estimated using one-sided difference scheme
    vort = uv2vr_cfd(u, v, lat, lon, 2)
    copy_VarCoords(u,vort)
    vort@units = "s~S~-1~N~"
    vort@name = "Relative vorticity on model levels"

  ; Coriolis parameter and absolute vorticity
    f1   = coriolis_param(lat)
    f0   = conform_dims(dimsizes(vort),f1,0)
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
     
     if (lay .eq. 1) then 
      str0   = sqrt( str1(:,:,:) + str2(:,:,:) )
     else 
      str0   = sqrt( str1(:,:) + str2(:,:) )
     end if

;===============================
; Consolidate these variables 
;===============================

   ; Read in relative vorticity and rate of strain
     if (lay .eq. 1) then 
      vort_plane = dim_avg_n_Wrap(vort(:,:,:),0)
      str_plane	 = dim_avg_n_Wrap(str0(:,:,:),0)
      copy_VarCoords(u(0,:,:),str)
     else
      vort_plane = vort
      str_plane	 = str0
      copy_VarCoords(u,str_plane)
     end if 

   ; EDIT - ALSO ADD CORRECT UNITS
     str0@description = "Rate of strain"
     str0@name        = "Rate of strain"

   ; Choose size of panel plot 
     panel       = new(2,graphic)

;=====================================================                         
; Create correct date strings for each output time                             
;=====================================================                         

  month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
                 "Jul","Aug","Sep","Oct","Nov","Dec"/)

  utc_date   = cd_calendar(time_c, 0)
  year       = tointeger(utc_date(:,0))
  month      = tointeger(utc_date(:,1))
  day        = tointeger(utc_date(:,2))
  hour       = tointeger(utc_date(:,3))
  minute     = tointeger(utc_date(:,4))
  second     = utc_date(:,5)

; Correct for errors in the code (round up value of hour when minutes = 59)    
  do tt = 0, dimsizes(time_c)-1
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

  time_str(ct)  = sprinti("%0.2i UTC ", hour(it) ) + \
                  sprinti("%0.2i ", day(it) ) \
                + month_abbr(month(it) )
  time_arr(ct)  = sprinti("%0.2i", day(it) ) + \
                  month_abbr(month(it) ) + \
                  "_" + sprinti("%0.2iZ", hour(it) )

;=================================
; Output file type and location
;=================================

     if (lay .eq. 1) then 
      output = "$sam/haima/images/"+dat+"_str_rv_"+pres(cp0)+"_"+pres(cp1)+"hpa_"+time_arr(ct)
     elseif (lay .eq. 0) then 
      output = "$sam/haima/images/"+dat+"_str_rv_"+pres(cp0)+"hpa_"+time_arr(ct)
     end if 

     wks = gsn_open_wks(opt,output)

   ; Load colour table
     gsn_define_colormap(wks,"prcp_new")

;==========================
; Options for plotting
;==========================

   ; Relative vorticity
     opts_rv                               = True
     opts_rv@cnFillOn                      = True
     opts_rv@cnLineLabelInterval           = 2.0
     opts_rv@cnLineLabelFontHeightF        = 0.012
     opts_rv@cnLineLabelBackgroundColor    = "transparent"
     opts_rv@cnLineLabelPlacementMode      = "constant"
     opts_rv@cnLinesOn                     = False ; Contour lines off
     opts_rv@cnInfoLabelOn                 = False ; Contour labels off
     opts_rv@cnLevelSelectionMode          = "ExplicitLevels"

     opts_rv@cnLevels                     = (/5., 10., 15., 20., 25., 30., \
                                              40., 50., 60., 70., 80./)
     opts_rv@cnFillColors                 = (/0,4,5,6,7,8,9,10,11,12,13,14/)

     opts_rv@gsnPaperOrientation           = "landscape"
     opts_rv@tiMainString                  = ""
     opts_rv@tiMainFontHeightF             = 0.0125
     opts_rv@gsnLeftString                 = ""
     opts_rv@gsnRightString                = ""
     opts_rv@lbLabelBarOn                  = True       ; Individual lb on/off

     opts_rv@pmLabelBarWidthF              = 0.52
     opts_rv@pmLabelBarHeightF             = 0.13
     opts_rv@lbLabelFontHeightF            = 0.0125
     opts_rv@lbPerimOn                     = False
     opts_rv@lbLabelFont                   = "Helvetica"

     if (lay .eq. 0) then 

      opts_rv@lbTitleString                 = pres(cp0)+\
       			        	        " hPa relative vorticity (10~S~-4~N~ s~S~-1~N~)"  

     elseif (lay .eq. 1) then 

      opts_rv@lbTitleString             = "Layer-averaged ("+pres(cp0)+"-"+pres(cp1)+\
      					  " hPa) ~C~     relative vorticity (10~S~-4~N~ s~S~-1~N~)"

     end if 

     opts_rv@lbBoxEndCapStyle              = "TriangleBothEnds"    ; Labelbar end shape
     opts_rv@lbTitleFontHeightF            = 0.02
     opts_rv@lbTitleDirection              = "Across"
     opts_rv@lbTitlePosition               = "Bottom"

   ; Additional plotting resources
     opts_rv@mpDataBaseVersion             = "Ncarg4_1" ; More recent database
     opts_rv@mpDataSetName                 = "Earth..4" ; High resolution
     opts_rv@mpOutlineBoundarySets         = "National" ; National borders
     opts_rv@mpGeophysicalLineColor        = "black"    ; Colour borders black
     opts_rv@mpGeophysicalLineThicknessF   = 1.0        ; Border line thickness
     opts_rv@mpGridAndLimbOn               = False      ; Turn on lat/lon lines
     opts_rv@pmTickMarkDisplayMode         = "Always"   ; Turn on map tickmarks
     opts_rv@tmXBMajorLengthF              = 0.005      ; Change tickmark length
     opts_rv@tmXTOn                        = "False"    ; No tickmarks on top x-axis
     opts_rv@tmYROn                        = "False"    ; No tickmarks on right y-axis
     opts_rv@gsnMaximize                   = True       ; Maximise plot size

     r                                     = r0         ; Radius of plot (user input)

     opts_rv@gsnAddCyclic                  = False
     opts_rv@gsnDraw                       = False      ; Do not draw the plot
     opts_rv@gsnFrame                      = False      ; Do not advance the frame

   ; When reading in 'pd' stream data with 119 time intervals (PV)
   ; use the index one less than that in 'centre'
     opts_rv@mpLimitMode                   = "Corners"
     opts_rv@mpLeftCornerLatF              = centre(ct,0)-r
     opts_rv@mpLeftCornerLonF              = centre(ct,1)-r
     opts_rv@mpRightCornerLatF             = centre(ct,0)+r
     opts_rv@mpRightCornerLonF             = centre(ct,1)+r

;=========================================================
; Now plot rate of strain using different colour table
;=========================================================

   ; Rate of strain 
     opts_str                              = opts_rv
     opts_str@cnFillOn                     = True
     opts_str@cnLineLabelInterval          = 2.0
     opts_str@cnLineLabelFontHeightF       = 0.012
     opts_str@cnLineLabelBackgroundColor   = "transparent"
     opts_str@cnLineLabelPlacementMode     = "constant"
     opts_str@cnLinesOn                    = False ; Contour lines off
     opts_str@cnInfoLabelOn                = False ; Contour labels off
     opts_str@cnFillPalette                = "circ_budget"
     opts_str@cnLevelSelectionMode         = "ExplicitLevels"
     delete([/opts_str@cnLevels, opts_str@cnFillColors/])
     opts_str@cnLevels                     = (/10,12,14,16,18,20,22,25,30/)
     opts_str@cnFillColors                 = (/-1,11,13,15,16,17,18,19,20,21/)
     opts_str@gsnPaperOrientation          = "landscape"
     opts_str@tiMainString                 = ""
     opts_str@tiMainFontHeightF            = 0.0125
     opts_str@gsnLeftString                = ""
     opts_str@gsnRightString               = ""
     opts_str@lbLabelBarOn                 = True             ; Individual labelbar

     opts_str@pmLabelBarWidthF             = 0.52
     opts_str@pmLabelBarHeightF            = 0.13
     opts_str@lbLabelFontHeightF           = 0.0125
     opts_str@lbPerimOn                    = False
     opts_str@lbLabelFont                  = "Helvetica"

     if (lay .eq. 0) then 
      opts_str@lbTitleString = pres(cp0)+" hPa strain rate (10~S~-4~N~ s~S~-2~N~)"
     elseif (lay .eq. 1) then 
      opts_str@lbTitleString = "Layer-averaged ("+pres(cp0)+"-"+pres(cp1)+\
      			       " hPa) ~C~         strain rate (10~S~-4~N~ s~S~-2~N~)"
     end if 

     opts_str@lbBoxEndCapStyle             = "TriangleBothEnds"    ; Labelbar end shape
     opts_str@lbTitleFontHeightF           = 0.02
     opts_str@lbTitleDirection             = "Across"
     opts_str@lbTitlePosition              = "Bottom"

     opts_str@mpLeftCornerLatF             = centre(ct,0)-r
     opts_str@mpLeftCornerLonF             = centre(ct,1)-r
     opts_str@mpRightCornerLatF            = centre(ct,0)+r
     opts_str@mpRightCornerLonF            = centre(ct,1)+r

;=============
; Plot data
;=============

     str_plane   = str_plane * (10 ^ 4)
     str_plane@units = "10~S~-4~N~ s~S~-2~N~"

     vort_plane  = vort_plane * (10 ^ 4)
     vort_plane@units = "10~S~-4~N~ s~S~-1~N~"

     plot_rv     = gsn_csm_contour_map(wks,vort_plane,opts_rv) ; Relative vorticity

     panel(0)    = plot_rv

     plot_str    = gsn_csm_contour_map(wks,str_plane,opts_str) ; Strain rate
     panel(1)    = plot_str

;==================================================
; Add NH tropical cyclone symbol to panel plots
;==================================================

  if (tc .eq. 1) then 

    markers    = new((/dimsizes(panel)/),graphic)

  ; Add NH tropical cyclone symbol to index
    tc0 = NhlNewMarker(wks, "p", 137, 0.0, 0.0, 1.0, 1.0, 0.0)

    mres0               = True
    mres0@gsMarkerIndex = tc0         ; TC symbol
    mres0@gsMarkerSizeF = 20.0        ; Marker size
    mres0@gsMarkerColor = "grey55"    ; Marker colour
    mres0@gsMarkerThicknessF = 5.0    ; Marker thickness (increase)

    markers(0) = gsn_add_polymarker(wks,panel(0),centre(it-1,1),centre(it-1,0),mres0)
    markers(1) = gsn_add_polymarker(wks,panel(1),centre(it-1,1),centre(it-1,0),mres0)

  end if 

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    optsP                       = True
    optsP@gsnFrame              = False    ; Do not advance the frame
    optsP@gsnPanelLabelBar      = False    ; Turn off panel labelbar
    optsP@txString              = "Haima (2016): Initialised at "+dat+", "\
    				  +time_str(ct)
    optsP@gsnPanelFigureStrings = ""
    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.90
    optsP@gsnPanelBottom        = 0.10
    optsP@amJust                = "TopLeft"
    optsP@gsnPanelFigureStringsFontHeightF = 0.0125 ; Label size (default 0.01)

    gsn_panel(wks,panel,(/1,2/),optsP)
    frame(wks)

  ; Tidy up before next iteration
    delete([/vort_plane, str_plane, str0, str1, str2, str1a, str2a/])
    delete([/du_dx,dv_dy,dv_dx,du_dy/])
    delete([/lon,lat,lon1,lat1,u,v,vort,avo,f0,f1/])

  ; Counter variable (time)
    ct = ct + 1
 
  end do      ; End time loop (do it = 0, numINPUT-1)

 end do

end