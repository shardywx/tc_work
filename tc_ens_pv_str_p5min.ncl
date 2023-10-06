; Script to identify vacillation cycles (VCs) in UM ensemble simulations

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_rm.ncl"

; ncl 'dat="02T12"' 'opt="x11"' p_cp=3 dist=3.0 r0=1.0 'ens0="em11"' \
; ts=0 tf=155 var0=0 var1=1 full=0 tc=0 lay=0 nep_cp_ens_vc5min_p.ncl 

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'ens0'  = ensemble simulation (em00-em11)
; 'dist'  = size of box following storm (degrees)
; 'opt'   = output file format ("pdf" or "x11")
; 'p_cp'  = pressure level for plotting PV and vertical velocity (see L106-108)
; 'r0'    = radius for final plot (degrees)
; 'var0'  = plot either PV (0) or relative vorticity (1)
; 'var1'  = plot either theta-w (0) or theta-e (1)
; 'ts'	  = starting time index (0-116)
; 'tf'	  = ending time index (0-116)
; 'full'  = produce 3-panel (1) or 2-panel (0) plot
; 'tc'    = annotate images with TC symbol (1) or leave as they were (0)
; 'lay'   = calculate layer average (1)

begin

;==============================
; Prelash (reading in files)
;==============================

; Troubleshooting
  if (r0 .gt. dist) then 
   print("Choose a bigger domain ('dist') or smaller plot size ('r0')!")
   exit()
  end if  

; Find 'pb' stream data to read in and plot 
  diri        = "/nfs/a319/earshar/02T12/em11/vc"
  fili_prefix = "20160702T1200Z_NPTK_4p4_L80_ra1t_em11"
  fili_b1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pb*.nc")
  fili_b      = diri+"/"+fili_b1

; Calculate number of times in file using built-in function ('nc_times')
  numINFO_b   = nc_times(fili_b)
  numINPUT_b  = numINFO_b[0]
  time_b      = numINFO_b[1]

; Calculate storm motion using built-in function ('st_rm')
; 'ens0'    = ensemble simulation
; 'diri'    = path to input files (see above)
; 'dist'    = size of box following storm (degrees)
; 'mins'    = analyse 1-h (0) or 5-min (1) data

  storm_rel = st_rm(dat,ens0,"$ar/text/",dist,1)
  u_cyc     = storm_rel[0]             ; Zonal wind
  v_cyc     = storm_rel[1]             ; Meridional wind
  vel_cyc   = storm_rel[2]             ; Vector wind
  numTIMES  = storm_rel[3]             ; Number of times in file
  lat_arr   = storm_rel[4]             ; Latitude array subset (following storm)
  lon_arr   = storm_rel[5]             ; Longitude array subset (following storm)
  centre    = storm_rel[6]             ; Storm track (position) array 

; Create arrays for data and output information
  llbox	    = toint(dist*50)		       ; Calculate domain size (grid points)
  dsize	    = (/llbox,llbox/)		       ; Domain size (grid points)

; Extract string for each ensemble member
  str1    = str_split(fili_b(0),"_")
  ens_str = str1(5)
  delete(str1)

  print("Working on ensemble member "+ens_str)

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 10000000000
  end setvalues

;=====================================================
; Create correct date strings for each output time
;=====================================================

    month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
                      "Jul","Aug","Sep","Oct","Nov","Dec"/)

    utc_date = cd_calendar(time_b, 0)

    year   = tointeger(utc_date(:,0))
    month  = tointeger(utc_date(:,1))
    day    = tointeger(utc_date(:,2))
    hour   = tointeger(utc_date(:,3))
    minute = tointeger(utc_date(:,4))
    second = utc_date(:,5)

  ; Create array to hold finished date strings
    hr_min = new(numINPUT_b,string)
    min0   = (/4,9,14,19,24,29,34,39,44,49,54,59/)      ; Minute array (before)
    min1   = (/5,10,15,20,25,30,35,40,45,50,55,0/)      ; Minute array (after)

  ; Correct for errors in the code (round up values of 'minute' to multiples of 5)
    do i = 0, numINPUT_b-1
     do m = 0, dimsizes(min0)-1

      if (minute(i) .eq. min0(m) ) then
       minute(i) = min1(m)
      end if

     end do
    end do

  ; Correct for instances where 'minute' = 59 and 'hour' is incorrectly set to the previous hour
    do i = 23, 131, 36
     hour(i) = hour(i) + 1
    end do

  ; Create string array of all values in 'minute'
    mins = tostring(minute)

  ; Create 'hhmm' strings for each time (for output and title strings)
    do i = 0, numINPUT_b-1

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

    time_str  = new(numINPUT_b,string)
    time_arr  = new(numINPUT_b,string)
    title_arr = new(numINPUT_b,string)

    time_str  = hr_min + " UTC " + sprinti("%0.2i ", day) + month_abbr(month)
    time_arr  = sprinti("%0.2i", day) + month_abbr(month) + "_" + hr_min + "Z"
    title_arr = "Valid at "+time_str

;==================================
; Get the variables we will need 
;==================================

  f = addfile(fili_b,"r")
  
  do it = ts, tf		   ; Loop over times in file (T+54 to T+66, every 5 min)

   time  = f->t(it)                ; Time in file
   times = dimsizes(time)
   print_clock("Working on time: "+time_str(it))

 ; Get domain boundaries at each time from storm track (calculated off-line)
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
 ; (0) 1000, (1) 950, (2) 925, (3) 900, (4) 850, (5) 800, (6) 750, (7) 700
 ; (8) 650, (9) 600, (10) 500, (11) 400, (12) 300, (13) 250
 ; (14) 200, (15) 150, (16) 100

   pres  = f->p(p_cp)               ; Subset of pressure level array
   plevs = dimsizes(pres)           ; Size of pressure level array
   pres@units = "hPa"               ; Pressure level units (for plotting)

 ; Arrays have () vertical levels, and horizontal extent depends on 'dist' option
   u0     = f->u(it,p_cp,{t0:t1},{n0:n1})       ; Zonal velocity (m/s)
   v0     = f->v(it,p_cp,{t0:t1},{n0:n1})       ; Meridional velocity (m/s)
   vort0  = f->rvor(it,p_cp,{t0:t1},{n0:n1})    ; Relative vorticity (/s)
   pv0    = f->pv(it,p_cp,{t0:t1},{n0:n1})      ; PV (PVU)

   t00    = f->temp(it,p_cp,{t0:t1},{n0:n1})    ; Temperature (K)
   w0     = f->dz_dt(it,p_cp,{t0:t1},{n0:n1})   ; Vertical velocity (m/s)
   rh0    = f->rh(it,p_cp,{t0:t1},{n0:n1})      ; Relative humidity (%)

   vel0   = sqrt((u0^2) + (v0^2))             ; Vector velocity (m/s)
   copy_VarCoords(v0,vel0)
   copy_VarMeta(v0,vel0)

 ; Avoid dimension size errors throughout the script 
   usize  = dimsizes(u0)
   vsize  = dimsizes(v0)
   vrsize = dimsizes(vort0)
   pvsize = dimsizes(pv0)
   tsize  = dimsizes(t00)
   wsize  = dimsizes(w0)
   vcsize = dimsizes(vel0)
   rhsize = dimsizes(rh0)

   if ( any (usize .gt. dsize(0)) ) then 
    u = u0(0:dsize(0)-1,0:dsize(0)-1)   
   else 
    u = u0
   end if 

   if ( any (vsize .gt. dsize(0)) ) then
    v = v0(0:dsize(0)-1,0:dsize(0)-1)
   else	
    v =	v0
   end if

   if ( any (vrsize .gt. dsize(0)) ) then
    vort = vort0(0:dsize(0)-1,0:dsize(0)-1)
   else	 
    vort = vort0
   end if

   if ( any (pvsize .gt. dsize(0)) ) then
    pv = pv0(0:dsize(0)-1,0:dsize(0)-1)
   else 
    pv = pv0
   end if

   if ( any (tsize .gt. dsize(0)) ) then
    t = t00(0:dsize(0)-1,0:dsize(0)-1)
   else 
    t = t00
   end if

   if ( any (wsize .gt. dsize(0)) ) then
    w = w0(0:dsize(0)-1,0:dsize(0)-1)
   else 
    w = w0
   end if

   if ( any (vcsize .gt. dsize(0)) ) then
    vel = vel0(0:dsize(0)-1,0:dsize(0)-1)
   else 
    vel = vel0
   end if

   if ( any (rhsize .gt. dsize(0)) ) then
    rh   = rh0(0:dsize(0)-1,0:dsize(0)-1)
   else
    rh   = rh0
   end if

 ; Convert to units we want
   vort  = vort * (10 ^ -6) ; Convert to /s
   pv    = pv * (10 ^ 6)    ; Convert to PVU [10 ^ 6 K m2 s-1 kg-1]

   dx    = lon(1) - lon(0) ; Grid spacing (longitude)
   dy    = lat(1) - lat(0) ; Grid spacing (latitude)

 ;===============================================================================
 ; Calculate Coriolis parameter and absolute vorticity from relative vorticity
 ;===============================================================================

   f1   = coriolis_param(lat(0:dsize(0)-1))
   f0   = conform_dims(dimsizes(vort),f1,0)	; Give 'f1' the dimensions of 'vort' 
   copy_VarMeta(vort,f0)
   copy_VarCoords(vort,f0)

   avo  = vort + f0
   copy_VarMeta(vort,avo)
   copy_VarAtts(vort,avo)
   copy_VarCoords(vort,avo)
   avo@description = "Absolute vorticity"
   avo@long_name   = "Absolute vorticity"

 ; Calculate storm-relative winds
   uf = u - u_cyc(it)
   vf = v - v_cyc(it)

 ; Add metadata from original horizontal wind arrays
   copy_VarCoords(u,uf)
   copy_VarAtts(u,uf)

   copy_VarCoords(v,vf)
   copy_VarAtts(v,vf)

;===========================================================
; Calculate rate of strain (following Nguyen et al. 2011)
;===========================================================

     dim1 = "latitude"
     dim2 = "longitude"

   ; Calculate horizontal derivatives of 'u' (du_dx, du_dy)
     du    = grad_latlon_cfd(uf,uf&$dim1$,uf&$dim2$,False,False)
     du_dy = du[0]    ; Latitudinal derivative
     du_dx = du[1]    ; Longitudinal derivative
     delete(du)

   ; Calculate horizontal derivatives of 'v' (dv_dx, dv_dy)
     dv    = grad_latlon_cfd(vf,vf&$dim1$,vf&$dim2$,False,False)
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

   ; Read in PV, geopotential height and vertical velocity
     vort_plane = vort
     pv_plane   = pv
     w_plane    = w

     str_plane  = str_plane0
     copy_VarCoords(u,str_plane)

     if (full .eq. 1) then 
      panel      = new(3,graphic)       ; Panel plot with 3 images
     elseif (full .eq. 0) then 
      panel      = new(2,graphic)       ; Panel plot with 2 images
     end if 

;=================================
; Output file type and location
;=================================

     output = "$sam/nepartak/images/vc/"+ens0+"_"+dat+"_5min_str_pv_"+pres+"_"+time_arr(it)
     wks = gsn_open_wks(opt,output)

   ; Load colour table
     if (var0 .eq. 0) then
      gsn_define_colormap(wks,"prcp_new")		; Potential vorticity
     elseif (var0 .eq. 1) then
      gsn_define_colormap(wks,"bwr_new")		; Relative vorticity
     end if

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

     if (var0 .eq. 0) then 
      opts_pv@cnLevels                     = (/9.0, 10.0 ,12.0, 14.0, \
                                               16.0, 18.0, 20.0, 22.0, 25.0, \
                                               30.0, 35.0, 40.0, 45.0, 50.0/)
      opts_pv@cnFillColors                 = (/0,2,3,4,5,\
                                               6,7,8,9,10,\
                                               11,12,13,14,15/) 
     elseif (var0 .eq. 1) then 
      opts_pv@cnLevels                     = (/-5., -2., -1., 1., 2., 5., 10., 20., \
                                                30., 40., 60., 80., 100., 150./)
      opts_pv@cnFillColors                 = (/13,14,15,0,\
                                               18,19,21,23,24,26,27,28,30,31,32/)
     end if 

     opts_pv@gsnPaperOrientation           = "landscape"
     opts_pv@tiMainString                  = ""
     opts_pv@tiMainFontHeightF             = 0.0125
     opts_pv@gsnLeftString                 = ""
     opts_pv@gsnRightString                = ""
     opts_pv@lbLabelBarOn                  = True       ; Individual lb on/off

     opts_pv@pmLabelBarWidthF              = 0.52
     opts_pv@pmLabelBarHeightF             = 0.13
     opts_pv@lbLabelFontHeightF            = 0.0125
     opts_pv@lbPerimOn                     = False
     opts_pv@lbLabelFont                   = "Helvetica"

     if (var0 .eq. 0) then 
      opts_pv@lbTitleString                = pres+" hPa potential vorticity (PVU)"
     elseif (var0 .eq. 1) then 
      opts_pv@lbTitleString                = pres+" hPa relative vorticity (10~S~-4~N~ s~S~-1~N~)"
     end if      

     opts_pv@lbBoxEndCapStyle              = "TriangleBothEnds"    ; Labelbar end shape
     opts_pv@lbTitleFontHeightF            = 0.02
     opts_pv@lbTitleDirection              = "Across"
     opts_pv@lbTitlePosition               = "Bottom"

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

     r                                     = r0         ; Radius of plot (user input)

     opts_pv@gsnAddCyclic                  = False
     opts_pv@mpLimitMode                   = "Corners"
     opts_pv@mpLeftCornerLatF              = centre(it,0)-r
     opts_pv@mpLeftCornerLonF              = centre(it,1)-r
     opts_pv@mpRightCornerLatF             = centre(it,0)+r
     opts_pv@mpRightCornerLonF             = centre(it,1)+r
     opts_pv@gsnDraw                       = False      ; Do not draw the plot
     opts_pv@gsnFrame                      = False      ; Do not advance the frame

   ; Vertical velocity
     opts_vv                               = True
     opts_vv@cnFillOn                      = False
     opts_vv@cnLineColor                   = "black"
     opts_vv@cnInfoLabelOn                 = False
     opts_vv@cnLineLabelsOn                = True
     opts_vv@cnLineLabelInterval           = 2.0
     opts_vv@cnLevelSelectionMode          = "ExplicitLevels"
     opts_vv@cnLevels                      = (/1.5/)
     opts_vv@cnLineLabelPlacementMode      = "constant"
     opts_vv@cnLineLabelPerimOn            = False
     opts_vv@gsnContourLineThicknessesScale = 4.0
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
     opts_str@cnFillPalette                = "circ_budget"
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
     opts_str@lbLabelBarOn                 = True             ; Individual labelbar

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

;================================================
; Finally plot wet bulb potential temperature 
;================================================

   if (full .eq. 1) then 

   ; Read in "radar_new" colour map
     cmap = read_colormap_file("radar_new")
   
   ; Make the first 2 colours (black and white) transparent
     cmap(0:1,3) = 0.0

   ; Wet bulb potential temperature
     opts_the                              = opts_pv
     opts_the@cnFillOn                     = True
     opts_the@cnLineLabelInterval          = 2.0
     opts_the@cnLineLabelFontHeightF       = 0.012
     opts_the@cnLineLabelBackgroundColor   = "transparent"
     opts_the@cnLineLabelPlacementMode     = "constant"
     opts_the@cnLinesOn                    = False ; Contour lines off
     opts_the@cnInfoLabelOn                = False ; Contour labels off
     opts_the@cnFillPalette                = "radar_new"
     opts_the@cnLevelSelectionMode         = "ExplicitLevels"
     delete(opts_the@cnLevels)
     delete(opts_the@cnFillColors)
     opts_the@cnLevels			   = ispan(tavg0-5,tavg0+6,1)
     opts_the@cnFillColors		   = (/2,3,4,5,6,7,8,9,10,11,12,13,14/)
     opts_the@gsnPaperOrientation          = "landscape"
     opts_the@tiMainString                 = ""
     opts_the@tiMainFontHeightF            = 0.0125
     opts_the@gsnLeftString                = ""
     opts_the@gsnRightString               = ""
     opts_the@lbLabelBarOn                 = True             ; Individual labelbar

     opts_the@pmLabelBarWidthF             = 0.52
     opts_the@pmLabelBarHeightF            = 0.13
     opts_the@lbLabelFontHeightF           = 0.0125
     opts_the@lbPerimOn                    = False
     opts_the@lbLabelFont                  = "Helvetica"

     opts_the@lbTitleString = pres+" hPa theta-e (K)"

     opts_the@lbBoxEndCapStyle             = "TriangleBothEnds"    ; Labelbar end shape
     opts_the@lbTitleFontHeightF           = 0.02
     opts_the@lbTitleDirection             = "Across"
     opts_the@lbTitlePosition              = "Bottom"

   end if 

;=============
; Plot data
;=============

     str_plane   = str_plane * (10 ^ 4)
     str_plane@units = "10~S~-4~N~ s~S~-2~N~"

     vort_plane  = vort_plane * (10 ^ 4)
     vort_plane@units = "10~S~-4~N~ s~S~-1~N~"

     if (var0 .eq. 0) then 
      plot_pv    = gsn_csm_contour_map(wks,pv_plane,opts_pv)   ; Potential vorticity
     elseif (var0 .eq. 1) then 
      plot_pv    = gsn_csm_contour_map(wks,vort_plane,opts_pv) ; Relative vorticity     
     end if 

     plot_vv     = gsn_csm_contour(wks,w_plane,opts_vv)       ; Vertical velocity
     overlay(plot_pv, plot_vv)
     panel(0)   = plot_pv

     plot_str   = gsn_csm_contour_map(wks,str_plane,opts_str) ; Strain rate
     panel(1)   = plot_str

;     if (full .eq. 1) then 
;
;      plot_the  = gsn_csm_contour_map(wks,the_plane,opts_the) ; Equivalent potential temp
;      panel(2)	= plot_the
;
;     end if 

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

    do m = 0, dimsizes(panel)-1
     markers(m) = gsn_add_polymarker(wks,panel(m),centre(it,1),centre(it,0),mres0)
    end do      

  end if 

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    optsP                       = True
    optsP@gsnFrame              = False    ; Do not advance the frame
    optsP@gsnPanelLabelBar      = False    ; Turn off panel labelbar
    optsP@txString              = "Nepartak (2016): Simulation "+ens0+", "\
                                  +title_arr(it)
    optsP@gsnPanelFigureStrings = ""
    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.90
    optsP@gsnPanelBottom        = 0.10
    optsP@amJust                = "TopLeft"
    optsP@gsnPanelFigureStringsFontHeightF = 0.0125 ; Label size (default 0.01)

    if (full .eq. 1) then
     gsn_panel(wks,panel,(/1,3/),optsP)
     frame(wks)
    elseif (full .eq. 0) then
     gsn_panel(wks,panel,(/1,2/),optsP)
     frame(wks)
    end if

  ; Tidy up before next iteration
    delete([/pv_plane, w_plane/])
    delete([/str_plane, str_plane0, str1, str2/])

  ; Tidy up (avoid dimension size errors)
    delete([/lat1,lon1,lat,lon,time,times,w,pv,u,v,t,vel,rh/])
    delete([/w0,pv0,u0,v0,t00,vel0,vort0,rh0/])
    delete([/f0,f1,avo,vort/])
    delete([/du_dx,dv_dy,dv_dx,du_dy/])

  end do      ; End time loop (do it = 0, numINPUT-1)

end 