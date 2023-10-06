; Compare storm-tracking methods and overlay on relative vorticity or windspeed 

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_rm.ncl"

; ncl dat=\"02T12\" opt=\"x11\" ens0=\"em11\" rad0=2.0 dist=3.0 vc=0 
; zoom=1 lgd0=\"tr\" wrt=0 plt=1 cn0=\"geo_sm\" ts0=51 tf0=52 tc_fc_vort_sd.ncl

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

;=====================================
; Read in 'pc' and 'pd' stream data
;=====================================

  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0

; Horizontal windspeed on model levels 
  fili_c1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pc.nc")
  fili_c      = diri+"/"+fili_c1

; Pressure on model levels 
  fili_d1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pd.nc")
  fili_d      = diri+"/"+fili_d1 

; Get number of times in file using built-in function ('nc_times')
; [0] = number of time intervals
; [1] = information corresponding to time intervals
  numINFO_c   = nc_times(fili_c)
  numINPUT_c  = numINFO_c[0]
  time_c      = numINFO_c[1]

  numINFO_d   = nc_times(fili_d)
  numINPUT_d  = numINFO_d[0]
  time_d      = numINFO_d[1]
  num_l       = 59 

; Calculate storm motion using built-in function ('st_rm')
; 'dat'     = initialisation time (02T12, 03T00, etc) 
; 'ens0'    = ensemble simulation (em00-em11) 
; 'diri'    = path to input files (see above)
; 'dist'    = size of box following storm (degrees)
; 'mins'    = analyse 1-h (0) or 5-min (1) data

  storm_rel = st_rm(dat, ens0, "$ar/text/", dist, 0)
  u_cyc     = storm_rel[0]             ; Zonal wind 
  v_cyc     = storm_rel[1]             ; Meridional wind 
  vel_cyc   = storm_rel[2]             ; Vector wind
  numTIMES  = storm_rel[3]             ; Number of times in file 
  lat_arr   = storm_rel[4]             ; Latitude array subset (following storm) 
  lon_arr   = storm_rel[5]             ; Longitude array subset (following storm) 
  centre    = storm_rel[6]             ; Storm track (position) array 

; Create arrays for data and output information
  llbox     = toint(dist*50)                   ; Calculate domain size (grid points)
  dsize     = (/llbox,llbox/)                  ; Domain size (grid points) 

; Choose which storm track to read in 
  if (cn0 .eq. "slp") then
   cn = 0
  elseif (cn0 .eq. "vort") then
   cn = 1
  elseif (cn0 .eq. "geo") then
   cn = 2
  elseif (cn0 .eq. "geo_sm") then
   cn = 3
  end if

;========================================================================== 
; Create date/time string arrays for all times in file (use 'pc' stream)
;========================================================================== 

  month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
                    "Jul","Aug","Sep","Oct","Nov","Dec"/)

  utc_date = cd_calendar(time_d, 0)

  year   = tointeger(utc_date(:,0))
  month  = tointeger(utc_date(:,1))
  day    = tointeger(utc_date(:,2))
  hour   = tointeger(utc_date(:,3))
  minute = tointeger(utc_date(:,4))
  second = utc_date(:,5)

; Correct for errors in the code (round up value of hour when minutes = 59) 
  do ct = 0, numINPUT_d-1
   if (minute(ct).gt.30) then
     hour(ct) = hour(ct)+1
   end if
  end do
  delete(ct)

; Create arrays to hold strings for output plots
  date_str  = new(numINPUT_d,string)
  time_str  = new(numINPUT_d,string)
  time_arr  = new(numINPUT_d,string)
  title_arr = new(numINPUT_d,string)

  date_str  = sprinti("%0.2iUTC ", hour) + \
              sprinti("%0.2i ", day) \
            + month_abbr(month)
  time_str  = sprinti("%0.2i UTC ", hour) + \
              sprinti("%0.2i ", day) \
            + month_abbr(month)
  time_arr  = sprinti("%0.2i", day) + \
              month_abbr(month) + \
              "_" + sprinti("%0.2iZ", hour)

  print("Working on ensemble member "+ens0+" ("+dat+")")

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 10000000000
  end setvalues

;==========================================================
; Read in required variables from 'pc' and 'pd' streams 
;==========================================================

; Create array to hold storm track positions using different methods 
  centre_new  = new((/6,num_l,numINPUT_d,2/),"float")

  c  = addfile(fili_c,"r")              ; Read in 'pc' stream [u,v]
  d  = addfile(fili_d,"r")              ; Read in 'pd' stream [w,p,t,theta,kmh,tau]

  do it	       = ts0, tf0               ; Loop over times in file 

   print("Working on time: "+time_str(it))
   itt = it + 1
   title_arr(it) = "Valid at "+time_str(it)+" (T+"+itt+")"

 ; Get domain boundaries at each time from storm track (calculated off-line)
   t0          = lat_arr(cn,it-1,0)
   t1          = lat_arr(cn,it-1,1)
   n0          = lon_arr(cn,it-1,0)
   n1          = lon_arr(cn,it-1,1)

 ; Read in basic variables from 'pc' stream
   lon1C  = c->longitude_1({n0:n1})              ; longitude ['d1' grid points]
   lat1C  = c->latitude_1({t0:t1})               ; latitude  ['d0' grid points]
   lonC   = c->longitude({n0:n1})                ; longitude ['d1' grid points]
   latC   = c->latitude({t0:t1})                 ; latitude  ['d0' grid points]
   hybC   = c->hybrid_ht(:)                      ; 63 model ('theta') levels

   dy     = latC(1) - latC(0)
   dx     = lonC(1) - lonC(0)

   u_varname = "x-wind"
   v_varname = "y-wind"

 ; Read in horizontal wind components 
   if (lay .eq. 1) then 
    u0  = c->$u_varname$(it,13:25,{t0:t1},{n0:n1}) ; Zonal wind
    v0  = c->$v_varname$(it,13:25,{t0:t1},{n0:n1}) ; Meridional wind
   else
    u0  = c->$u_varname$(it,:,{t0:t1},{n0:n1})     ; Zonal wind 
    v0  = c->$v_varname$(it,:,{t0:t1},{n0:n1})     ; Meridional wind
   end if 

 ; Storm-relative wind calculations
   u = u0 - u_cyc(cn,it-1)
   v = v0 - v_cyc(cn,it-1)

 ; Add metadata from original horizontal wind arrays 
   copy_VarCoords(u0,u)
   copy_VarAtts(u0,u)

   copy_VarCoords(v0,v)
   copy_VarAtts(v0,v)

 ; Calculate relative vorticity using centered finite differences 
   vort = uv2vr_cfd(u, v, latC, lonC, 2)
   copy_VarCoords(u, vort)
   vort@units = "s~S~-1~N~"
   vort@name  = "Relative vorticity on model levels"

;======================================
; Read in variables from 'pd' stream
;======================================

 ; Read in basic variables 
   lon1D  = d->longitude_1({n0:n1})                ; longitude [xxx grid points]
   lat1D  = d->latitude_1({t0:t1})                 ; latitude  [xxx grid points]
   lonD   = d->longitude({n0:n1})                  ; longitude [xxx grid points]
   latD   = d->latitude({t0:t1})                   ; latitude  [xxx grid points]
   hybD   = d->hybrid_ht(:)                        ; 63 model ('rho') levels 

 ; Read in pressure and vertical velocity 
   if (lay .eq. 1) then  
    p      = d->p(it-2,13:25,{t0:t1},{n0:n1})      ; Pressure (Pa)
    w      = d->dz_dt(it-2,13:25,{t0:t1},{n0:n1})  ; Vertical velocity (m s-1)
    z      = d->ht(it-2,13:25,{t0:t1},{n0:n1})     ; Geopotential height (m)
   else
    p      = d->p(it-2,:,{t0:t1},{n0:n1})          ; Pressure (Pa)
    w      = d->dz_dt(it-2,:,{t0:t1},{n0:n1})      ; Vertical velocity (m s-1) 
    z      = d->ht(it-2,:,{t0:t1},{n0:n1})         ; Geopotential height (m) 
   end if 

 ; Convert units of pressure and geopotential height before continuing 
   p      = p / 100
   z      = z / 10

 ; Add metadata
   p@units = "hPa"
   z@units = "dam"

;==========================
; Loop over model levels 
;==========================

;  do ilev = 0, num_l-1 
  do ilev = 10, 46, 12

   hy = sprintf("%0.0f",hybD(ilev))
   print("Working on model level " + ilev+ " ("+hy+" m AGL)")

 ; Read in variables on model levels 
   u_plane      = u(ilev,:,:)
   v_plane      = v(ilev,:,:)
   
   spd_plane    = sqrt( (u_plane ^ 2) + (v_plane ^ 2) ) 
   copy_VarCoords(u_plane, spd_plane)

   vvel_plane   = w(ilev,:,:)

   vort_plane   = vort(ilev,:,:)
   vort_pl_smth = smth9_Wrap(vort_plane, 0.5, 0.5, True) 

   prs_plane    = p(ilev,:,:) 
   prs_pl_smth  = smth9_Wrap(prs_plane, 0.5, 0.5, True)
   prs_pl_smth2 = smth9_Wrap(prs_pl_smth, 0.5, 0.5, True)

;====================================================
; Use centre already calculated using MSLP as base
;====================================================

   lt1		  = centre(cn,it-1,0) - 1.0
   lt2		  = centre(cn,it-1,0) + 1.0
   ln1		  = centre(cn,it-1,1) - 1.0
   ln2		  = centre(cn,it-1,1) + 1.0

 ; RELATIVE VORTICITY (maximum)
   vort_plane0    = vort_plane({lt1:lt2},{ln1:ln2})
   vort_max0	  = max(vort_plane0)
   dims_v0        = dimsizes(vort_plane0)
   vort_1d0       = ndtooned(vort_plane0)
   inds_v0        = ind_resolve(maxind(vort_1d0),dims_v0)

   lat_max_v0 	  = lt1 + (dy * inds_v0(0,0))
   lon_max_v0 	  = ln1 + (dx * inds_v0(0,1))
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
   prs_min_smth	  = min(prs_smth)
   dims_p1	  = dimsizes(prs_smth)
   prs_1d1	  = ndtooned(prs_smth)
   inds_p1	  = ind_resolve(minind(prs_1d1),dims_p1)

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

 ; PRE-CALCULATED STORM POSITION 
   centre_new(5,ilev,it,0) = centre(cn,it-1,0)
   centre_new(5,ilev,it,1) = centre(cn,it-1,1)

 ; Calculate geopotential height max/min values for later contour plotting
   p_max0	  = toint(max(prs_plane))
   p_min0	  = toint(prs_min0)
   p_avg0	  = ( (p_max0 + p_min0) / 2)

; Output file location and type
  output = "$sam/nepartak/images/basic/fc_all_"+dat+"_"+ens0+"_"+hy+"m_"+time_arr(it)
  wks    = gsn_open_wks(opt,output)

; Load colour table 
  gsn_define_colormap(wks,"prcp_new")

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
  opts_rvort@lbTitleString                = "Relative vorticity at "+hy+" m (10~S~-4~N~ s~S~-1~N~)"
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
   opts_rvort@mpMinLatF                   = centre(cn,it-1,0) - 0.5
   opts_rvort@mpMinLonF                   = centre(cn,it-1,1) - 0.5
   opts_rvort@mpMaxLatF                   = centre(cn,it-1,0) + 0.5
   opts_rvort@mpMaxLonF                   = centre(cn,it-1,1) + 0.5
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
   opts_rvort@tmXBValues	           = fspan(100.0, 150.0, 101)
   xb_labels0				   = sprintf("%0.1f~S~o~N~E",opts_rvort@tmXBValues)
   opts_rvort@tmXBLabels		   = xb_labels0

   opts_rvort@tmYLMode                     = "Explicit"
   opts_rvort@tmYLValues                   = fspan(0.0, 30.0, 61)
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
  opts_rvort@gsnMaximize                  = True 

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

; Single panel 
  vort_plot  = gsn_csm_contour_map(wks,vort_plane,opts_rvort)   ; Relative vorticity
  prs_plot   = gsn_csm_contour(wks,prs_plane,opts_prs)          ; Pressure
;  overlay(vort_plot, prs_plot)

;==================================================
; Overlay markers for each storm centre estimate
;==================================================

; (1) vort, (2) vort_smth, (3) pressure, (4) pressure_smth, (5), pressure_smth2, (6) offline
  m_colours  = (/"royalblue1", "royalblue1", "orange4", "orange4", "orange4", "red4"/)
  m_indices  = (/12, 16, 12, 16, 8, 12/)
  m_sizes    = (/10.0, 10.0, 10.0, 10.0, 10.0, 10.0/)
  m_thick    = (/5.0, 5.0, 5.0, 5.0, 5.0, 5.0/)

  markers    = new( (/dimsizes(m_colours)/), graphic)

  do i = 0, dimsizes(m_colours)-1
    mres0                    = True
    mres0@gsMarkerIndex      = m_indices(i)
    mres0@gsMarkerSizeF      = m_sizes(i)
    mres0@gsMarkerColor      = m_colours(i)
    mres0@gsMarkerThicknessF = m_thick(i)
    markers(i)               = gsn_add_polymarker(wks,vort_plot,centre_new(i,ilev,it,1),\
                                                  centre_new(i,ilev,it,0),mres0)
  end do

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
  lg_opts@lgMarkerColors             = m_colours
  lg_opts@lgMarkerIndexes            = m_indices
  lg_opts@lgMarkerSizeF              = mres0@gsMarkerSizeF
  lg_opts@lgMarkerThicknessF         = mres0@gsMarkerThicknessF
  lsize                              = dimsizes(lg_opts@lgMarkerColors)

  lg_opts@lgLabelPosition            = "Right"           ; Label position within box
  lg_opts@lgItemPlacement            = "ExplicitPlacement"
  lg_opts@lgItemPositions            = fspan(0.12, 0.88, lsize)
  lg_opts@lgLeftMarginF		     = 0.01
  lg_opts@lgRightMarginF	     = 0.15

  lg_opts@vpWidthF                   = 0.300
  lg_opts@vpHeightF                  = 0.180
  lg_opts@lgPerimColor               = "black"
  lg_opts@lgPerimThicknessF          = 3.0
  lg_opts@lgPerimFill                = "SolidFill"
  lg_opts@lgPerimFillColor           = "white"

  lab_arr                            = (/"vort ("+hy+" m)", "vort_sm ("+hy+" m)", \
  				         "pres ("+hy+" m)", "pres_sm ("+hy+" m)", \
                                         "pres_sm2 ("+hy+" m)", "offline (900 hPa)"/)
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

;==============================
; Add radial circles to plot 
;==============================

  dtr = 0.017453292519943

; Set up properties of radial circles to overlay
  thetas   = new(73,float)
  pii      = 3.14159265
  radii    = fspan(0.0, 0.5, 11)
  rad_size = dimsizes(radii)-1

; Get customizations for circles
  radius   = radii(rad_size)
  xcenter  = centre(cn,it-1,1)
  ycenter  = centre(cn,it-1,0)
  out_thck = 2
  in_thck  = 1
  spacing  = 0.05
  dpattrn  = 2
  out_clr  = 1
  in_clr   = 1

;=============================================
; Construct and attach outer circle to plot 
;=============================================

; Calculate arrays for outer circle
  degrees  = ispan(0, 360, 5)
  xcos     = cos(dtr * degrees)
  xsin     = sin(dtr * degrees)

; Array of points (x,y) representing the outer circle
  xc       = xcenter + (radius * xcos)
  yc       = ycenter + (radius * xsin)

; Resources for outer circle 
  lnres    = True 
  lnres@gsLineThicknessF  = out_thck
  lnres@gsLineColor       = out_clr
  lnres@gsLineDashPattern = 0

; Attach outer circle to plot
  plot_outer = gsn_add_polyline(wks, vort_plot, xc, yc, lnres)

  if (spacing .gt. radius) then
   print("add_radial_circles: spacing is > radius, can't draw inner circles.")
   return
  end if

;==================================================
; Now construct and attach inner circles to plot 
;==================================================

; Draw inner circles if desired 
  if (spacing .gt. 0) then 

   count = 0
   size  = toint( (radius-spacing) / spacing )
   plot_inner = new(size,graphic)

   do r = spacing, radius, spacing 

   ; Calculate arrays for inner circles
     xc = xcenter + (r * xcos)
     yc = ycenter + (r * xsin)

   ; Resources for inner circles
     delete([/lnres@gsLineThicknessF,lnres@gsLineColor,lnres@gsLineDashPattern/])
     lnres@gsLineColor       = in_clr
     lnres@gsLineThicknessF  = in_thck
     lnres@gsLineDashPattern = dpattrn

   ; Attach inner circles to plot 
     plot_inner(count) = gsn_add_polyline(wks, vort_plot, xc, yc, lnres)

     count = count + 1 

   end do 

  end if 

; Tidy up 
  delete([/lnres, degrees, xcos, xsin, xc, yc/])

  draw(vort_plot)
  frame(wks)

  delete([/opts_prs@cnLevels/])

  delete([/vvel_plane,u_plane,v_plane,spd_plane/])
  delete([/vort_plane,prs_plane,vort_plane0,prs_plane0,prs_pl_smth,prs_pl_smth2,vort_pl_smth/])
  delete([/vort_1d0,vort_1d1,prs_1d0,prs_1d1,prs_1d2/])
  delete([/vort_smth,prs_smth,prs_smth2/]) 

 end do                                          ; End loop over model levels 

 delete([/vort,u,v,u0,v0,latC,lonC,lat1C,lon1C/]) ; Tidy up 
 delete([/p,w,z,latD,lonD,lat1D,lon1D/])

end do                ; End loop over input times (do it = 0, numINPUT_d-1)

end 