; Script to plot vorticity and strain rate on model height levels, following a TC 

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/times_arr.ncl"

; ncl dat=\"02T12\" opt=\"x11\" dist=3.0 r0=1.0 ens0=\"em11\" sr=1
; ts=563 tf=671 ti=6 var0=\"rv\" tc=0 cn0=\"slp\" plt=1 zoom=1 rad0=3.0 lg0=\"tl\" tc_fc5.ncl

; 563 --> 671 = T+48 --> T+57
; 131 --> 923 = T+12 --> T+78

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'ens0'  = ensemble simulation (em00-em11)
; 'dist'  = size of box following storm (degrees)
; 'opt'   = output file format ("pdf" or "x11")
; 'mlev'  = model level for plotting PV and vertical velocity (see L106-108)
; 'r0'    = radius for final plot (degrees)
; 'var0'  = plot either PV ("pv") or relative vorticity ("rv")
; 'ts'	  = starting time index (0-116)
; 'tf'	  = ending time index (0-116)
; 'full'  = produce 3-panel (1) or 2-panel (0) plot
; 'tc'    = annotate images with TC symbol (1) or leave as they were (0)
; 'lay'   = calculate layer average (1)
; 'mlev2' = upper model level for calculating layer average (only valid when lay=1)
; 'sr'	  = calculate storm-relative vorticity (1) or leave alone (0) 
; 'cn0'   = storm tracking method ("slp" or "slpf" --> using finescale tracking method)

begin

;==============================
; Prelash (reading in files)
;==============================

; Troubleshooting
  if (r0 .gt. dist) then 
   print("Choose a bigger domain ('dist') or smaller plot size ('r0')!")
   exit()
  end if  

; Find data streams to read in and plot 
  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0+"/vc"
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0

  fili_c1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pc.nc")
  fili_d1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pd.nc")
  fili_c      = diri+"/"+fili_c1
  fili_d      = diri+"/"+fili_d1

; Calculate number of times in file using built-in function ('nc_times') - should be same here
; 'pc' stream --> u,v
; 'pd' stream --> z,p
  numINFO_c   = nc_times(fili_c)
  numINFO_d   = nc_times(fili_d)
  numINPUT_c  = numINFO_c[0]            ; Number of times in 'pc' stream
  time_c      = numINFO_c[1]            ; Array of times from 'pc' stream file
  numINPUT_d  = numINFO_d[0]            ; Number of times in 'pd' stream
  time_d      = numINFO_d[1]            ; Array of times from 'pd' stream file
  num_l       = 63                      ; Number of model levels 

; Calculate storm motion using built-in function ('st_rm')
; 'dat'     = initialisation time (02T12, 03T00, etc)
; 'ens0'    = ensemble simulation (em00-em11)
; 'diri'    = path to input files (see above)
; 'dist'    = size of box following storm (degrees)
; 'mins'    = analyse 1-h (0) or 5-min (1) data

  storm_rel = st_rm(dat, ens0, "$ar/text/", dist, 1)
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
  d0        = dsize(0)                         ; Dim 0 of 'dsize'
  d1        = dsize(1)                         ; Dim 1 of 'dsize'

; Increase system memory
  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 10000000000
  end setvalues

; Choose which storm track to read in                                                   
  cn = 0

;===========================================================================            
; Create date/time string arrays for all times with user-defined function               
;===========================================================================            

 ; Call external function 'times_arr'                                                   
 ; Returns 'time_str' (12 UTC 4 July), 'time_arr' (), and 'title_arr' ()                
   all_times  = times_arr(dat, ens0)

 ; Return strings e.g. '0600 UTC 04 July'                                               
   time_str   = all_times[0]

 ; Return strings e.g. '04Jul_0600Z'                                                    
   time_arr   = all_times[1]

 ; Return title strings e.g. 'Valid at 0600 UTC 04 July'                                
   title_arr  = all_times[2]

;================================================================
; Create arrays for relative vorticity and strain at each time 
;================================================================
  
 ; Number of time intervals to loop over 
   nt_all   = (tf - ts) + 1

 ; Array for storm track from different centre-finding methods 
   centre_new = new( (/6, num_l, numINPUT_c, 2/), "float")

;==================================
; Get the variables we will need 
;==================================

  c = addfile(fili_c,"r")          ; Read in 'pc' stream [u,v,vort]
  d = addfile(fili_d,"r")          ; Read in 'pd' stream [w,p]

  do it = ts, tf, ti		   ; Loop over times in file (T+54 to T+66, every 5 min)

   time  = c->t(it)                ; Time in file
   times = dimsizes(time)

   print_clock("Working on "+time_str(it))

 ; Get domain boundaries at each time from storm track (calculated off-line)
   t0    = lat_arr(cn,it,0)
   t1    = lat_arr(cn,it,1)
   n0    = lon_arr(cn,it,0)
   n1    = lon_arr(cn,it,1)
   print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

 ; Create latitude/longitude array subsets to reduce computation time
   lonC  = c->longitude({n0:n1})    ; longitude [160 grid points]
   latC  = c->latitude({t0:t1})     ; latitude  [160 grid points]
   hybC  = c->hybrid_ht(:)          ; 63 model ('theta') levels

   u_varname = "x-wind"
   v_varname = "y-wind"

 ; Read in horizontal wind components
 ; Both variables below -- [hybC | 63] * [latC | 150] * [lonC | 150]
   u0  = c->$u_varname$(it,:,{t0:t1},{n0:n1})           ; Zonal wind (single level)
   v0  = c->$v_varname$(it,:,{t0:t1},{n0:n1})           ; Meridional wind (single level)

 ; Calculate storm-relative winds
   u   = u0 - u_cyc(cn,it)
   v   = v0 - v_cyc(cn,it)

 ; Add metadata from original horizontal wind arrays
   copy_VarCoords(u0,u)
   copy_VarAtts(u0,u)

   copy_VarCoords(v0,v)
   copy_VarAtts(v0,v)

 ; Calculate relative vorticity using centered finite differences
   vort0 = uv2vr_cfd(u, v, latC, lonC, 2)
   copy_VarCoords(u, vort0)
   vort0@units = "s~S~-1~N~"
   vort0@name  = "Relative vorticity on model levels"

 ; Anticipate any dimension errors later in the script 
   vort        = vort0(:, 0:dsize(0)-1, 0:dsize(1)-1)

;=====================================================
; Read in variables from 'pd' stream (model levels)
;=====================================================

 ; Read in basic variables from 'pd' stream
   lonD   = d->longitude({n0:n1})                  ; longitude [xxx grid points]
   latD   = d->latitude({t0:t1})                   ; latitude  [xxx grid points]
   hybD   = d->hybrid_ht(:)                        ; 63 model ('rho') levels

;   hy     = hybD / 1000 ; Model level values in km

 ; Read in temperature, pressure, vertical velocity and geopotential height
 ; All variables below -- [hybD | 63] * [latD | 150] * [lon1D | 151]
   p      = d->p(it,:,{t0:t1},{n0:n1})             ; Pressure (Pa)

 ; Convert units of pressure before continuing
   p      = p / 100

 ; Add metadata
   p@units  = "hPa"

   dx     = lonD(1) - lonD(0) ; Grid spacing (longitude)
   dy     = latD(1) - latD(0) ; Grid spacing (latitude)

;==========================
; Loop over model levels 
;========================== 

   do ilev = 0, num_l-1

    hy = sprintf("%0.0f",hybD(ilev))
    print("Working on model level " + ilev+ " ("+hy+" m AGL)")
  
  ; Read in variables on model levels 
    u_plane   = u(ilev,:,:)
    v_plane   = v(ilev,:,:)
    spd_plane = sqrt( (u_plane ^ 2) + (v_plane ^ 2) )
    copy_VarCoords(u_plane, spd_plane)

    vort_plane   = vort(ilev,:,:)
    vort_pl_smth = smth9_Wrap(vort_plane, 0.5, 0.5, True)

    prs_plane    = p(ilev,:,:)
    prs_pl_smth  = smth9_Wrap(prs_plane, 0.5, 0.5, True)
    prs_pl_smth2 = smth9_Wrap(prs_pl_smth, 0.5, 0.5, True)

  ; Use centre already calculated using MSLP as base
    lt1          = centre(cn,it-1,0) - 1.0
    lt2          = centre(cn,it-1,0) + 1.0
    ln1          = centre(cn,it-1,1) - 1.0
    ln2          = centre(cn,it-1,1) + 1.0    

  ; RELATIVE VORTICITY (maximum)
    vort_plane0  = vort_plane({lt1:lt2},{ln1:ln2})
    vort_max0    = max(vort_plane0)
    dims_v0      = dimsizes(vort_plane0)
    vort_1d0     = ndtooned(vort_plane0)
    inds_v0      = ind_resolve(maxind(vort_1d0),dims_v0)

    lat_max_v0   = lt1 + (dy * inds_v0(0,0))
    lon_max_v0   = ln1 + (dx * inds_v0(0,1))
    print("Centre (vort): "+lat_max_v0+" degrees N, "+lon_max_v0+" degrees E")

    centre_new(0,ilev,it,0) = lat_max_v0
    centre_new(0,ilev,it,1) = lon_max_v0
 
  ; SMOOTHED RELATIVE VORTICITY (maximum) 
    vort_smth      = vort_pl_smth({lt1:lt2},{ln1:ln2})
    vort_max_smth  = max(vort_smth)
    dims_v1        = dimsizes(vort_smth)
    vort_1d1       = ndtooned(vort_smth)
    inds_v1        = ind_resolve(minind(vort_1d1),dims_v1)

    lat_max_v1     = lt1 + (dy * inds_v1(0,0))
    lon_max_v1     = ln1 + (dx * inds_v1(0,1))
    print("Centre (vort_smth): "+lat_max_v1+" degrees N, "+lon_max_v1+" degrees E")

    centre_new(1,ilev,it,0) = lat_max_v1
    centre_new(1,ilev,it,1) = lon_max_v1

  ; PRESSURE (minimum)
    prs_plane0     = prs_plane({lt1:lt2},{ln1:ln2})
    prs_min0       = min(prs_plane0)
    dims_p0        = dimsizes(prs_plane0)
    prs_1d0        = ndtooned(prs_plane0)
    inds_p0        = ind_resolve(minind(prs_1d0),dims_p0)

    lat_max_p0     = lt1 + (dy * inds_p0(0,0))
    lon_max_p0     = ln1 + (dx * inds_p0(0,1))
    print("Centre (prs): "+lat_max_p0+" degrees N, "+lon_max_p0+" degrees E")

    centre_new(2,ilev,it,0) = lat_max_p0
    centre_new(2,ilev,it,1) = lon_max_p0

  ; SMOOTHED PRESSURE (minimum) 
    prs_smth       = prs_pl_smth({lt1:lt2},{ln1:ln2})
    prs_min_smth   = min(prs_smth)
    dims_p1        = dimsizes(prs_smth)
    prs_1d1        = ndtooned(prs_smth)
    inds_p1        = ind_resolve(minind(prs_1d1),dims_p1)

    lat_max_p1     = lt1 + (dy * inds_p1(0,0))
    lon_max_p1     = ln1 + (dx * inds_p1(0,1))
    print("Centre (prs_smth): "+lat_max_p1+" degrees N, "+lon_max_p1+" degrees E")

    centre_new(3,ilev,it,0) = lat_max_p1
    centre_new(3,ilev,it,1) = lon_max_p1

  ; EXTRA SMOOTHED PRESSURE (minimum) 
    prs_smth2      = prs_pl_smth2({lt1:lt2},{ln1:ln2})
    prs_min_smth2  = min(prs_smth2)
    dims_p2        = dimsizes(prs_smth2)
    prs_1d2        = ndtooned(prs_smth2)
    inds_p2        = ind_resolve(minind(prs_1d2),dims_p2)

    lat_max_p2     = lt1 + (dy * inds_p2(0,0))
    lon_max_p2     = ln1 + (dx * inds_p2(0,1))
    print("Centre (prs_smth_2): "+lat_max_p2+" degrees N, "+lon_max_p2+" degrees E")

    centre_new(4,ilev,it,0) = lat_max_p2
    centre_new(4,ilev,it,1) = lon_max_p2

  ; PRE-CALCULATED STORM POSITION (MSLP)
    centre_new(5,ilev,it,0) = centre(cn,it,0)
    centre_new(5,ilev,it,1) = centre(cn,it,1)

  ; Calculate geopotential height max/min values for later contour plotting
    p_max0         = toint(max(prs_plane))
    p_min0         = toint(prs_min0)
    p_avg0         = ( (p_max0 + p_min0) / 2)

;=======================
; Now plot if desired 
;=======================

    if (plt .eq. 1) then 

   ; Output file location and type
     output = "$sam/nepartak/images/basic/fc5_"+var0+"_"+dat+"_"+ens0+"_"+hy+"m_"+time_arr(it)
     wks    = gsn_open_wks(opt,output)
  
   ; Load colour table 
     gsn_define_colormap(wks,"prcp_new")

   ; Create panel plot for desired variable (vector wind, relative vorticity, etc)
     panel  = new(1,graphic)

;========================
; Options for plotting 
;========================
 
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
     opts_rvort@cnLevels                     = (/5., 10., 15., 20., 25., 30., \
                                              40., 50., 60., 70., 80./)
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
     opts_rvort@lbTitleString                = \
                                             "Relative vorticity at "+hy+" m (10~S~-4~N~ s~S~-1~N~)"
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
      opts_rvort@mpMinLatF                   = centre(cn,it-1,0) - 1
      opts_rvort@mpMinLonF                   = centre(cn,it-1,1) - 1
      opts_rvort@mpMaxLatF                   = centre(cn,it-1,0) + 1
      opts_rvort@mpMaxLonF                   = centre(cn,it-1,1) + 1
     else
      opts_rvort@mpMinLatF                   = lat(0)
      opts_rvort@mpMinLonF                   = lon(0)
      opts_rvort@mpMaxLatF                   = lat(dimsizes(lat)-1)
      opts_rvort@mpMaxLonF                   = lon(dimsizes(lon)-1)
     end if

     opts_rvort@tmXTOn                       = "False"    ; No tickmarks on top x-axis
     opts_rvort@tmYROn                       = "False"    ; No tickmarks on right y-axis

     if (rad0 .le. 2.0) then
      opts_rvort@tmXBMode                     = "Explicit"
      opts_rvort@tmXBValues                   = fspan(100.0, 150.0, 51)
      xb_labels0                              = sprintf("%0.1f~S~o~N~E",   opts_rvort@tmXBValues)
      opts_rvort@tmXBLabels                   = xb_labels0

      opts_rvort@tmYLMode                     = "Explicit"
      opts_rvort@tmYLValues                   = fspan(0.0, 30.0, 31)
      yl_labels0                              = sprintf("%0.1f~S~o~N~N",   opts_rvort@tmYLValues)
      opts_rvort@tmYLLabels                   = yl_labels0
     elseif (rad0 .gt. 2.0) then
      opts_rvort@tmXBMode                     = "Explicit"
      opts_rvort@tmXBValues                   = fspan(100.0, 150.0, 26)
      xb_labels0                              = sprintf("%0.1f~S~o~N~E",   opts_rvort@tmXBValues)
      opts_rvort@tmXBLabels                   = xb_labels0

      opts_rvort@tmYLMode                     = "Explicit"
      opts_rvort@tmYLValues                   = fspan(0.0, 30.0, 16)
      yl_labels0                              = sprintf("%0.1f~S~o~N~N",   opts_rvort@tmYLValues)
      opts_rvort@tmYLLabels                   = yl_labels0
     end if

     opts_rvort@mpGridAndLimbOn              = True       ; Grid lines on/off
     opts_rvort@mpGridLineDashPattern        = 2          ; Dash pattern 
     opts_rvort@mpGridLatSpacingF            = 2.0        ; Spacing (latitude) 
     opts_rvort@mpGridLonSpacingF            = 2.0        ; Spacing (longitude)

   ; Windspeed
     opts_spd                                = opts_rvort
     delete([/opts_spd@cnLevels, opts_spd@cnFillColors, opts_spd@lbTitleString/])
     opts_spd@cnLevels                       = (/10,15,20,25,30,35,40,45,50,55/)
     opts_spd@cnFillColors                   = (/0,4,5,6,7,8,9,10,11,12,13,14,15/)
     opts_spd@lbTitleString                  = "Horizontal windspeed at "+hy+" m (m s~S~-1~N~)"

   ; Pressure 
     opts_prs                                = True
     opts_prs@cnFillOn                       = False
     opts_prs@cnLineColor                    = "black"
     opts_prs@cnInfoLabelOn                  = False
     opts_prs@cnLineLabelsOn                 = True
     opts_prs@cnLineLabelInterval            = 2.0
     opts_prs@cnLevelSelectionMode           = "ExplicitLevels"
     opts_prs@cnLevels                       = ispan(p_min0, p_max0, 4)
     opts_prs@cnLineLabelPlacementMode       = "constant"
     opts_prs@cnLineLabelPerimOn             = False
     opts_prs@gsnContourLineThicknessesScale = 3.0
     opts_prs@tiMainString                   = ""
     opts_prs@gsnLeftString                  = ""
     opts_prs@gsnRightString                 = ""
     opts_prs@gsnDraw                        = False      ; Do not draw the plot
     opts_prs@gsnFrame                       = False      ; Do not advance the frame
     opts_prs@pmTickMarkDisplayMode          = "Always"   ; Turn on map tickmarks 
     opts_prs@tmXBMajorLengthF               = 0.005      ; Change tickmark length
     opts_prs@tmXTOn                         = "False"    ; No tickmarks on top x-axis
     opts_prs@tmYROn                         = "False"    ; No tickmarks on right y-axis
     opts_prs@gsnMaximize                    = True       ; Maximise plot size

;=================
; Plot the data
;=================

     vort_plane = vort_plane * (10 ^ 4)
     vort_plane@units = "10~S~-4~N~ s~S~-1~N~"

   ; Decide which variable to contour                                                               
     if (var0 .eq. "rv") then
      main_plot = gsn_csm_contour_map(wks,vort_plane,opts_rvort)   ; Relative vorticity
     elseif (var0 .eq. "uv") then
      main_plot = gsn_csm_contour_map(wks,spd_plane,opts_spd)      ; Vector wind
     end if

     prs_plot   = gsn_csm_contour(wks,prs_plane,opts_prs)          ; Pressure
     overlay(main_plot, prs_plot)
     panel      = main_plot

;==================================================
; Overlay markers for each storm centre estimate
;================================================== 

   ; (1) vort, (2) vort_smth, (3) pressure, (4) pressure_smth, (5), pressure_smth2
     m_colours  = (/"royalblue1", "royalblue1", "orange4", "orange4", "orange4", "black"/)
     m_indices  = (/12, 16, 12, 16, 8, 4/)
     m_sizes    = (/10.0, 10.0, 10.0, 10.0, 10.0, 10.0/)
     m_thick    = (/5.0, 5.0, 5.0, 5.0, 5.0, 5.0/)
     markers    = new(dimsizes(m_colours),graphic)

     do i = 0, dimsizes(m_colours)-1
      mres0                    = True
      mres0@gsMarkerIndex      = m_indices(i)
      mres0@gsMarkerSizeF      = m_sizes(i)
      mres0@gsMarkerColor      = m_colours(i)
      mres0@gsMarkerThicknessF = m_thick(i)
      markers(i)               = gsn_add_polymarker(wks,panel,centre_new(i,ilev,it,1),\
                                                    centre_new(i,ilev,it,0),mres0)
     end do

;==============================
; Add and customise a legend
;==============================

     lg_opts                            = True
     lg_opts@pmLegendDisplayMode        = "Always"          ; Display a legend
     lg_opts@lgAutoManage               = False             ; Gives us control
     lg_opts@lgLabelFont                = "Helvetica"       ; Label font
     lg_opts@lgLabelFontHeightF         = 0.060             ; Height of labels

   ; Customise labels and line colours
     lg_opts@lgItemType                 = "Markers"
     lg_opts@lgMarkerColors             = m_colours
     lg_opts@lgMarkerIndexes            = m_indices
     lg_opts@lgMarkerSizeF              = mres0@gsMarkerSizeF
     lg_opts@lgMarkerThicknessF         = mres0@gsMarkerThicknessF
     lsize                              = dimsizes(lg_opts@lgMarkerColors)

     lg_opts@lgLabelPosition            = "Right"           ; Label position within box
     lg_opts@lgItemPlacement            = "ExplicitPlacement"
     lg_opts@lgItemPositions            = fspan(0.08, 0.92, lsize)
     lg_opts@lgLeftMarginF              = 0.01
     lg_opts@lgRightMarginF             = 0.15

     lg_opts@vpWidthF                   = 0.125
     lg_opts@vpHeightF                  = 0.090
     lg_opts@lgPerimColor               = "black"
     lg_opts@lgPerimThicknessF          = 3.0
     lg_opts@lgPerimFill                = "SolidFill"
     lg_opts@lgPerimFillColor           = "white"

     lab_arr                            = (/"vort ("+hy+" m)", "vort_sm ("+hy+" m)", \
                                            "pres ("+hy+" m)", "pres_sm ("+hy+" m)", \
                                            "pres_sm2 ("+hy+" m)", "offline (MSLP)"/)
     lgd                                = gsn_create_legend(wks,lsize,lab_arr,lg_opts)

   ; Add legend to plot (user-controlled)
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

     annotate_1                         = gsn_add_annotation(panel,lgd,am_opts)

;===================================================
; Finally, draw the plot with everything overlaid 
;===================================================

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

     delete([/opts_prs@cnLevels/])

    end if 
  
    delete([/u_plane,v_plane,spd_plane,vort_plane,vort_plane0,vort_pl_smth/])
    delete([/prs_plane,prs_plane0,prs_pl_smth,prs_pl_smth2/])
    delete([/vort_1d0,vort_1d1,prs_1d0,prs_1d1,prs_1d2/])
    delete([/vort_smth,prs_smth,prs_smth2/])

   end do     ; End model level loop 

   delete([/vort,u,v,u0,v0,latC,lonC,p,latD,lonD/])

  end do      ; End time loop (do it = 0, numINPUT_d-1) 

;===================================
; Write out values to text files 
;===================================

  diri_out = "$ar/text/levs5/cp_"

; Loop over times and write values on all model levels (59) to single text file 
  do it = ts, tf 

   print("Working on time: "+time_arr(it) )

 ; (1) Relative vorticity
   lat_out  = diri_out+dat+"_"+ens0+"_"+time_arr(it)+"_lat_v0.txt"
   asciiwrite(lat_out,centre_new(0,:,it,0) )

   lon_out  = diri_out+dat+"_"+ens0+"_"+time_arr(it)+"_lon_v0.txt"
   asciiwrite(lon_out,centre_new(0,:,it,1))

 ; (2) Smoothed relative vorticity
   lat_out  = diri_out+dat+"_"+ens0+"_"+time_arr(it)+"_lat_v1.txt"
   asciiwrite(lat_out,centre_new(1,:,it,0))

   lon_out  = diri_out+dat+"_"+ens0+"_"+time_arr(it)+"_lon_v1.txt"
   asciiwrite(lon_out,centre_new(1,:,it,1))

 ; (3) Pressure
   lat_out  = diri_out+dat+"_"+ens0+"_"+time_arr(it)+"_lat_p0.txt"
   asciiwrite(lat_out,centre_new(1,:,it,0))

   lon_out  = diri_out+dat+"_"+ens0+"_"+time_arr(it)+"_lon_p0.txt"
   asciiwrite(lon_out,centre_new(2,:,it,1))

 ; (4) Smoothed pressure
   lat_out  = diri_out+dat+"_"+ens0+"_"+time_arr(it)+"_lat_p1.txt"
   asciiwrite(lat_out,centre_new(3,:,it,0))

   lon_out  = diri_out+dat+"_"+ens0+"_"+time_arr(it)+"_lon_p1.txt"
   asciiwrite(lon_out,centre_new(3,:,it,1))

 ; (5) Doubly smoothed pressure
   lat_out  = diri_out+dat+"_"+ens0+"_"+time_arr(it)+"_lat_p2.txt"
   asciiwrite(lat_out,centre_new(4,:,it,0))

   lon_out  = diri_out+dat+"_"+ens0+"_"+time_arr(it)+"_lon_p2.txt"
   asciiwrite(lon_out,centre_new(4,:,it,1))

  end do

end 