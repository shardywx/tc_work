; Script to calculate area-averaged quantities around a developing TC 

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
   numTIMES    = 21
   centre      = new((/numTIMES,2/),float) ; Array to hold all [120] lat/lon values

 ; Read in latitude and longitude of cyclone track from text file
   centre(:,0) = asciiread(diri+"tc_ecm_lat.txt",(/numTIMES/),"float")
   centre(:,1) = asciiread(diri+"tc_ecm_lon.txt",(/numTIMES/),"float")

;============================================================================
; Calculate distance between two points --> then calculate system velocity
;============================================================================

   print_clock("Starting system velocity calculations...")

 ; Define new arrays (cyclone position, speed)
   x_cyc   = new(numTIMES,float)
   y_cyc   = new(numTIMES,float)
   u_cyc   = new(numTIMES,float)
   v_cyc   = new(numTIMES,float)
   vel_cyc = new(numTIMES,float)

 ; Loop to calculate system velocity from 6-hourly cyclone position [lat,lon]
   do iang = 1, numTIMES-1
    x_cyc(iang)   = gc_latlon(centre(iang-1,0),centre(iang-1,1),centre(iang-1,0),\
                              centre(iang,1),2,3)
    y_cyc(iang)   = gc_latlon(centre(iang-1,0),centre(iang-1,1),centre(iang,0),\
                              centre(iang-1,1),2,3)
    u_cyc(iang)   = x_cyc(iang)/(21600)                 ; 6 h = 21600 seconds
    v_cyc(iang)   = y_cyc(iang)/(21600)
    vel_cyc(iang) = sqrt(u_cyc(iang)^2.0 + v_cyc(iang)^2.0)
   end do

 ; Set initial system velocity [0] equal to that at the second time [1]
   u_cyc(0)   = u_cyc(1)
   v_cyc(0)   = v_cyc(1)
   vel_cyc(0) = vel_cyc(1)

 ; Add metadata to new arrays
   u_cyc@units         = "m s~S~-1~N~"
   u_cyc@description   = "Cyclone zonal velocity"

   v_cyc@units         = "m s~S~-1~N~"
   v_cyc@description   = "Cyclone meridional velocity"

   vel_cyc@units       = "m s~S~-1~N~"
   vel_cyc@description = "Cyclone speed"

   u_cyc = -(u_cyc) ; Easterly winds are negative

;================================================================================
; Create lat/lon arrays containing grid subsets at each time, following the TC
;================================================================================
  lat_arr = new((/numTIMES,2/),float)
  lon_arr = new((/numTIMES,2/),float)

; Lat/lon arrays from ECMWF analysis data
  lat_arr(:,0) = centre(:,0)-dist
  lat_arr(:,1) = centre(:,0)+dist
  lon_arr(:,0) = centre(:,1)-dist
  lon_arr(:,1) = centre(:,1)+dist

;==================================================================================
;;;;;;;;;;;;;;;;;;;;;;;;;;; Start of main part of script ;;;;;;;;;;;;;;;;;;;;;;;;;;
;==================================================================================

  print_clock("Starting main part of script...")

  ecm_input = (/"$sam/ecm/ecm_new_20160704.nc","$sam/ecm/ecm_new_20160705.nc",\
                "$sam/ecm/ecm_new_20160706.nc","$sam/ecm/ecm_new_20160707.nc",\
                "$sam/ecm/ecm_new_20160708.nc"/)

;===================================
; Prelash (before start of loops)
;===================================

   numINPUT  = dimsizes(ecm_input)    ; Number of input files
   FirstTime = True                   ; Skip first time when calculating vort tendencies
   r_earth   = 6.37e6                 ; Earth's radius (constant)
   numTIMES  = numINPUT * 4           ; Number of time intervals across ALL input files
   ct        = 0                      ; Counter variable
   plevs     = 14		      ; Number of presssure levels (subset)

   time_arr  = new((/numTIMES/),string)  ; Array to hold strings for second loop
   title_arr = new((/numTIMES/),string)  ; Array to hold strings for plot titles
   pos_arr   = new((/numTIMES,2/),float) ; Array for cyclone position [lat/lon]

   date_str  = new(numTIMES,string)
   time_str  = new(numTIMES,string)
   out_str   = new(numTIMES,string)

   setvalues NhlGetWorkspaceObjectId
     "wsMaximumSize" : 1000000000
   end setvalues

;================================
; Start loop over input files
;================================

  do nf = 0, numINPUT-1
  f = addfile(ecm_input(nf),"r")

;==================================
; Get the variables we will need
;==================================

  do it = 0, 3		   ; Loop over 00, 06, 12, 18 UTC in each file 

; 14 (vertical) x 241 (latitude) x 481 (longitude) 

  print("ct = "+ct)

  t0    = lat_arr(ct,0)
  t1    = lat_arr(ct,1)
  n0    = lon_arr(ct,0)
  n1    = lon_arr(ct,1)
  print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

  lon   = f->longitude({n0:n1})        ; Longitude (grid subset)
  lat   = f->latitude({t0:t1})         ; Latitude  (grid subset)

 ; Pressure levels
 ; (0)  150, (1) 200, (2) 250, (3) 300, (4) 400, (5) 500, (6) 600
 ; (7)  700, (8) 800, (9) 850, (10) 900, (11) 925, (12) 950, (13) 1000

  p_ec	= (/12,9,7,5,1/)	       ; 5 levels (950, 850, 700, 500, 200)
;  p_ec  = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13/)
  pres  = f->level(p_ec)               ; Pressure
  pres@units = "hPa"
  plevs = dimsizes(pres)

  time	= f->time(it)		     	  ; Time 	        [1 time]
  pvA	= f->pv(it,p_ec,{t0:t1},{n0:n1})  ; PV 	     	        [PVU]
  zA	= f->z(it,p_ec,{t0:t1},{n0:n1})   ; Geopotential        [m**2 / s**2]
  tA    = f->t(it,p_ec,{t0:t1},{n0:n1})   ; Temperature 	[K]     
  omgA  = f->w(it,p_ec,{t0:t1},{n0:n1})   ; Vertical velocity   [Pa/s]
  qA	= f->q(it,p_ec,{t0:t1},{n0:n1})   ; Specific humidity   [kg/kg]
  rhA	= f->r(it,p_ec,{t0:t1},{n0:n1})   ; Relative humidity   [%]
  uA    = f->u(it,p_ec,{t0:t1},{n0:n1})   ; Zonal velocity      [m/s]
  vA    = f->v(it,p_ec,{t0:t1},{n0:n1})   ; Meridional velocity [m/s]
  vortA = f->vo(it,p_ec,{t0:t1},{n0:n1})  ; Relative vorticity  [/s]
  divA  = f->d(it,p_ec,{t0:t1},{n0:n1})   ; Divergence	        [/s] 

  dx    = lon(1) - lon(0) ; Grid spacing (longitude)
  dy    = lat(1) - lat(0) ; Grid spacing (latitude)  

;==========================================================
; Unpack the data using scale factors/offset and tidy up
;==========================================================

  pv    = tofloat((pvA * pvA@scale_factor) + pvA@add_offset)
  copy_VarCoords(pvA, pv)
  copy_VarAtts(pvA, pv)
  copy_VarMeta(pvA, pv)
  delete(pv@missing_value_original)
  delete(pv@_FillValue_original)
  delete(pvA)

  q    = tofloat((qA * qA@scale_factor) + qA@add_offset)
  copy_VarCoords(qA, q)
  copy_VarAtts(qA, q)
  copy_VarMeta(qA, q)
  delete(q@missing_value_original)
  delete(q@_FillValue_original)
  delete(qA)

  rh   = tofloat((rhA * rhA@scale_factor) + rhA@add_offset)
  copy_VarCoords(rhA, rh)
  copy_VarAtts(rhA, rh)
  copy_VarMeta(rhA, rh)
  delete(rh@missing_value_original)
  delete(rh@_FillValue_original)
  delete(rhA)

  z    = tofloat((zA * zA@scale_factor) + zA@add_offset)
  copy_VarCoords(zA, z)
  copy_VarAtts(zA, z)
  copy_VarMeta(zA, z)
  delete(z@missing_value_original)
  delete(z@_FillValue_original)
  delete(zA)

  gr   = 9.80665
  z    = z / (10 * gr)
  z@units = "dam"

  t    = tofloat((tA * tA@scale_factor) + tA@add_offset)
  copy_VarCoords(tA, t)
  copy_VarAtts(tA, t)
  copy_VarMeta(tA, t)
  delete(t@missing_value_original)
  delete(t@_FillValue_original)
  delete(tA)

  u    = tofloat((uA * uA@scale_factor) + uA@add_offset)
  copy_VarCoords(uA, u)
  copy_VarAtts(uA, u)
  copy_VarMeta(uA, u)
  delete(u@missing_value_original)
  delete(u@_FillValue_original)
  delete(uA)

  v    = tofloat((vA * vA@scale_factor) + vA@add_offset)
  copy_VarCoords(vA, v)
  copy_VarAtts(vA, v)
  copy_VarMeta(vA, v)
  delete(v@missing_value_original)
  delete(v@_FillValue_original)
  delete(vA)

  omg    = tofloat((omgA * omgA@scale_factor) + omgA@add_offset)
  copy_VarCoords(omgA, omg)
  copy_VarAtts(omgA, omg)
  copy_VarMeta(omgA, omg)
  delete(omg@missing_value_original)
  delete(omg@_FillValue_original)
  delete(omgA)

  vort = tofloat((vortA * vortA@scale_factor) + vortA@add_offset)
  copy_VarCoords(vortA, vort)
  copy_VarAtts(vortA, vort)
  copy_VarMeta(vortA, vort)
  delete(vort@missing_value_original)
  delete(vort@_FillValue_original)
  delete(vortA)

  div  = tofloat((divA * divA@scale_factor) + divA@add_offset)
  copy_VarCoords(divA, div)
  copy_VarAtts(divA, div)
  copy_VarMeta(divA, div)
  delete(div@missing_value_original)
  delete(div@_FillValue_original)
  delete(divA)

  f1     = coriolis_param(lat)                       ; Coriolis parameter
  f0     = conform_dims(dimsizes(vort),f1,1)         ; Create same-sized array
  avo    = vort + f0                                 ; Absolute vorticity
  copy_VarCoords(vort,avo)
  copy_VarAtts(vort,avo)
  copy_VarMeta(vort,avo)

;=============================================
; Convert omega to vertical velocity in m/s
;=============================================

  prs = conform(omg,pres,0)      ; Create array with same dimensions
  prs = prs * 100                ; Convert to Pa
  prs@units = "Pa"               ; Add metadata

  w   = omega_to_w(omg, prs, t)  ; Convert from Pa/s to m/s
  w@units = "m s~S~-1~N~"        ; Add metadata

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
   if (minute.gt.30) then
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

   time_arr(ct) = out_str(it)

   print("Working on time: "+time_str(it))
   title_arr(ct) = "Valid at: "+time_str(it)

;========================================
; Define new arrays for advective flux
;========================================

; ; Storm-relative velocity
;   u_rel      = u ; - u_cyc(ct)
;   v_rel      = v ; - v_cyc(ct)
;   printVarSummary(u_rel)
;   printVarSummary(v_rel)

   d          = dimsizes(vort)
   d0         = d(0)-1
   d1         = d(1)-1
   d2         = d(2)-1
   pii        = 3.14159265           ; Constant
   lat_rad    = lat * (pii/180)      ; Convert to radians for calculation below
   lon_rad    = lon * (pii/180)

;==========================================================
; Define latitude and longitude for finite differencing
;==========================================================

   y2d     = new((/d(1),d(2)/),float)
   y2d!0   = "lat"
   y2d!1   = "lon"
   y2d&lat =  lat
   y2d&lon =  lon

   y2d     = conform_dims(dimsizes(y2d),lat,0) ; Array of latitude values
   coslat  = cos(y2d*pii/180)                      ; Same array in radians

 ; Determine horizontal resolution of data
   delta_lon = y2d
   delta_lat = y2d
   delta_lon = lon(1) - lon(0) ; Zonal grid spacing
   delta_lat = lat(1) - lat(0) ; Meridional grid spacing
   delta_lon = sqrt(delta_lon^2)
   delta_lat = sqrt(delta_lat^2)

   lons = lon(1) - lon(0)
   lats = lat(1) - lat(0)
   lons = sqrt(lons^2)
   lats = sqrt(lats^2)

 ; Calculate dx and dy (changes with longitude but not with latitude)
   ddx = y2d
   ddx = r_earth*(pii/180.0)*delta_lon*cos(y2d*pii/180)
   ddy = y2d
   ddy = r_earth*(pii/180)*delta_lat

 ; 3-D latitude array
   coslat3d = conform_dims(dimsizes(u),coslat,(/1,2/))

;========================
; Find cyclone centre
;========================

  ; New pressure levels
  ; (0)  950, (1) 850, (2) 700, (3) 500, (4) 200

  ; Find cyclone centre using 950 hPa vorticity and geopotential height 
    ilev = 0

  ; RELATIVE VORTICITY
    vort_plane    = vort(ilev,:,:)
    vort_max      = max(vort_plane)
    vort_smth     = smth9_Wrap(vort_plane, 0.5, 0.5, True)
    vort_max_smth = max(vort_smth)

    dims          = dimsizes(vort_plane)
    vort1d        = ndtooned(vort_plane)
    inds          = ind_resolve(maxind(vort1d),dims)
    vort1d_sm     = ndtooned(vort_smth)
    inds_sm       = ind_resolve(maxind(vort1d_sm),dims)

    lat_max950 = t0 + (dy * inds(0,0)) ; Latitude of max. vorticity
    lon_max950 = n0 + (dx * inds(0,1)) ; Longitude of max. vorticity
    print("Cyclone centre (vort): "+lat_max950+" degrees N, "+lon_max950+" degrees E")

  ; GEOPOTENTIAL HEIGHT
    geo_plane	   = z(ilev,:,:)
    geo_min        = min(geo_plane)
    geo_smth       = smth9_Wrap(geo_plane, 0.5, 0.5, True)
    geo_min_smth   = min(geo_smth)

    dims_h         = dimsizes(geo_plane)
    geo1d          = ndtooned(geo_plane)
    inds_h         = ind_resolve(minind(geo1d),dims_h)
    geo1d_sm       = ndtooned(geo_smth)
    inds_h_sm      = ind_resolve(minind(geo1d_sm),dims_h)

    lat_min950 = t0 + (dy * inds_h(0,0))
    lon_min950 = n0 + (dx * inds_h(0,1))
    print("Cyclone centre (geo): "+lat_min950+" degrees N, "+lon_min950+" degrees E")  

   ;========================================================
   ; EMPLOY A SAFETY NET IF WE IDENTIFY THE WRONG CYCLONE
   ;========================================================
     if (ct.eq.0) then
      centre(ct,0) = lat_max950
      centre(ct,1) = lon_max950

     else
      lt0   = centre(ct-1,0)            ; Previous TC latitude
      ln0   = centre(ct-1,1)            ; Previous TC longitude
      print("lt0 = "+lt0+" ; ln0 = "+ln0)
      d_lat = abs(lat_max950-lt0)       ; Change in latitude (t1 - t0)
      d_lon = abs(lon_max950-ln0)       ; Change in longitude (note extra minus sign)
      print("d_lat = "+d_lat+" ; d_lon = "+d_lon)

      if (d_lat.gt.5.or.d_lon.gt.5) then        ; Recalculate TC centre if incorrect

     ; Create smaller grid [0.5 degrees]
       inc = 0.5
       lt1 = lt0 - inc
       lt2 = lt0 + inc
       ln1 = ln0 - inc
       ln2 = ln0 + inc

       vort_plane0 = vort(ilev,{lt1:lt2},{ln1:ln2})        ; Smaller grid (1.0 x 1.0 deg)
       vort_max0   = max(vort_plane0)                      ; Find minimum on smaller grid
       dims_v0     = dimsizes(vort_plane0)                 ; Size of smaller grid
       vort_1d0    = ndtooned(vort_plane0)                 ; Create 1-D array
       inds_v0     = ind_resolve(minind(vort_1d0),dims_v0) ; Find index of minimum

       lat_max950 = lt1 + (dy * inds_v0(0,0))
       lon_max950 = ln1 + (dy * inds_v0(0,1))
       print("New centre (vort): "+lat_max950+" degrees N, "+lon_max950+" degrees E")
       centre(ct,0) = lat_max950
       centre(ct,1) = lon_max950

       delete([/vort_plane0, vort_max0, dims_v0, vort_1d0, inds_v0/])

      else
       print("Cyclone centre (vort): "+lat_max950+" degrees N, "+lon_max950+" degrees E")
       centre(ct,0) = lat_max950
       centre(ct,1) = lon_max950
      end if

     end if

   ; Tidy up before looping over pressure levels
     delete([/vort_plane, vort1d, inds, vort1d_sm, vort_smth, inds_sm/])

;==============================
; Loop over pressure levels 
;==============================

     ilev  = clev
     ilev2 = clev2
     p     = pres(ilev)
     p2	   = pres(ilev2)
     print("Working on: " + p + " and "+p2+ " hPa")

   ; Calculate variables on single pressure level, ready for plotting 
     vort_plane    = vort(ilev,:,:)
     avo_plane	   = avo(ilev,:,:)
     geo_plane     = z(ilev,:,:)
     u_plane	   = u(ilev,:,:)
     v_plane	   = v(ilev,:,:)
     t_plane	   = t(ilev,:,:)
     w_plane	   = w(ilev,:,:)

   ; Calculate variables on second level for panel plot
     vort_plane1    = vort(ilev2,:,:)
     avo_plane1     = avo(ilev2,:,:)
     geo_plane1     = z(ilev2,:,:)
     w_plane1       = w(ilev2,:,:)

   ; Max/min geopotential on pressure level
     zmin          = toint(min(geo_plane))
     zmax          = toint(max(geo_plane))

     zmin1          = toint(min(geo_plane1))
     zmax1          = toint(max(geo_plane1))

     copy_VarCoords(avo_plane,u_plane)
     copy_VarCoords(avo_plane,v_plane)
     
     spd_plane	   = sqrt((u_plane^2)+(v_plane^2))
     copy_VarCoords(u_plane,spd_plane)
     copy_VarMeta(u_plane,spd_plane)

;=============================================================
; Open workstation for plotting (time-height cross section)
;=============================================================

;     output = "$sam/nepartak/images/cart_coords/ecm_overview_"+p+"_"+time_arr(ct)
     output = "$sam/nepartak/images/cart_coords/ecm_overview_"+time_arr(ct)
     wks = gsn_open_wks(opt,output)

;=========================
; Options for plotting 
;=========================

    panel = new(2,graphic)
    gsn_define_colormap(wks,"vort_panel")

  ; Absolute vorticity
    opts_vort                            = True
    opts_vort@cnFillOn                   = True
    opts_vort@cnLinesOn                  = False
    opts_vort@cnLineLabelInterval        = 2.0
    opts_vort@cnLineLabelFontHeightF     = 0.012
    opts_vort@cnLineLabelBackgroundColor = "transparent"
    opts_vort@cnLineLabelPlacementMode   = "constant"
    opts_vort@cnLinesOn                  = False ; Contour lines off
    opts_vort@cnInfoLabelOn              = False
    opts_vort@cnLevelSelectionMode       = "ExplicitLevels"
;    opts_vort@cnLevels                   = (/-1.0, 0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, \
;                                              12.0, 15.0, 20.0, 25.0, 30.0, \
;                                              35.0, 40.0, 45.0, 50.0/)
    opts_vort@cnLevels                   = (/-1.0, 0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, \
                                              12.0, 14.0, 16.0, 18.0, 20.0, \
                                              22.0, 24.0, 26.0, 28.0/)
    opts_vort@cnFillColors               = (/6,8,0,12,13,15,17,18,19,\
                                             20,21,22,23,24,\
                                             25,26,27,28/)
    opts_vort@gsnPaperOrientation        = "landscape"
    opts_vort@tiMainString               = ""
    opts_vort@tiMainFontHeightF          = 0.0125
    opts_vort@gsnLeftString              = ""
    opts_vort@gsnRightString             = ""
    opts_vort@lbLabelBarOn               = False

  ; Additional plotting resources
    opts_vort@mpDataBaseVersion          = "Ncarg4_1"     ; More recent database
    opts_vort@mpDataSetName              = "Earth..4"     ; High resolution
    opts_vort@mpOutlineBoundarySets      = "National"     ; National borders
    opts_vort@mpGeophysicalLineColor     = "black"        ; Colour borders black
    opts_vort@mpGeophysicalLineThicknessF = 1.0           ; Border line thickness
    opts_vort@mpGridAndLimbOn            = True           ; Turn on lat/lon lines
    opts_vort@mpGridLineDashPattern      = 2              ; Dash pattern
    opts_vort@mpGridSpacingF             = 2.0            ; Spacing between grid lines
    opts_vort@pmTickMarkDisplayMode      = "Always"       ; Turn on map tickmarks
    opts_vort@tmXBMajorLengthF           = 0.005          ; Change tickmark length
    opts_vort@tmXBMinorOn                = False          ; Turn off minor tm (x-axis)
    opts_vort@tmYLMinorOn                = False          ; Turn off minor tm (y-axis)
    opts_vort@tmXTOn                     = "False"        ; Turn off tm (top x-axis)
    opts_vort@tmYROn                     = "False"        ; Turn off tm (right y-axis)
    opts_vort@gsnMaximize                = False          ; Maximise plot size
    opts_vort@gsnPaperOrientation        = "landscape"    ; change orientation
    opts_vort@gsnAddCyclic               = False
    opts_vort@gsnDraw                    = False          ; Do not draw the plot
    opts_vort@gsnFrame                   = False          ; Do not advance the frame
    opts_vort@mpLimitMode                = "Corners"

    radius                               = rad            ; Size of plot
    opts_vort@mpLeftCornerLatF           = lat_min950-radius
    opts_vort@mpLeftCornerLonF           = lon_min950-radius
    opts_vort@mpRightCornerLatF          = lat_min950+radius
    opts_vort@mpRightCornerLonF          = lon_min950+radius

  ; Geopotential height
    opts_geo                              = True
    opts_geo@cnLineColor                  = "black"
    opts_geo@cnInfoLabelOn                = False
    opts_geo@cnLevelSelectionMode         = "ExplicitLevels"
    opts_geo@cnLevels                     = ispan(zmin,zmax,4) ; Depends on min/max
    opts_geo@cnLineLabelInterval          = 2.0
    opts_geo@cnLineLabelPlacementMode     = "constant"
    opts_geo@gsnContourLineThicknessesScale = 2.0
    opts_geo@cnLineLabelFontHeightF       = 0.01
    opts_geo@gsnMaximize                  = True       ; Maximise plot size
    opts_geo@tiMainString                 = ""
    opts_geo@gsnLeftString                = ""
    opts_geo@gsnRightString               = ""
    opts_geo@gsnDraw                      = False      ; Do not draw the plot
    opts_geo@gsnFrame                     = False      ; Do not advance the frame

    opts_geo1				  = opts_geo
    delete(opts_geo1@cnLevels)
    opts_geo1@cnLevels			  = ispan(zmin1,zmax1,4)

  ; Wind vectors
    opts_vec                            = True
    opts_vec@gsnDraw                    = False         ; Do not draw the plot
    opts_vec@gsnFrame                   = False         ; Do not advance the frame
    opts_vec@vcGlyphStyle               = "LineArrow"   ; draw windbarbs
    opts_vec@vcRefMagnitudeF            = 10            ; magnitude of reference vector
    opts_vec@vcRefLengthF               = 0.04          ; length of reference vector
    opts_vec@vcMinDistanceF             = 0.05          ; distance between wind barbs
    opts_vec@vcLineArrowThicknessF      = 1.0           ; windbarb line thickness
    opts_vec@vcLineArrowColor           = "black"       ; windbarb colour
    opts_vec@gsnMaximize                = True          ; Maximise plot size
    opts_vec@gsnRightString             = ""
    opts_vec@gsnLeftString              = ""
    opts_vec@tiXAxisString              = ""
    opts_vec@vcRefAnnoOrthogonalPosF    = 0.05

   ; Vertical velocity
     opts_vvel                               = True
     opts_vvel@cnFillOn                      = False
     opts_vvel@cnLineColor                   = "navyblue"
     opts_vvel@cnInfoLabelOn                 = False
     opts_vvel@cnLineLabelsOn                = True
     opts_vvel@cnLineLabelInterval           = 2.0
     opts_vvel@cnLevelSelectionMode          = "ExplicitLevels"
     opts_vvel@cnLevels                      = (/0.40,0.60,0.80/)
     opts_vvel@cnLineLabelPlacementMode      = "constant"
     opts_vvel@cnLineLabelPerimOn            = False
     opts_vvel@gsnContourLineThicknessesScale = 3.0
     opts_vvel@tiMainString                  = ""
     opts_vvel@gsnLeftString                 = ""
     opts_vvel@gsnRightString                = ""
     opts_vvel@gsnDraw                       = False      ; Do not draw the plot
     opts_vvel@gsnFrame                      = False      ; Do not advance the frame
     opts_vvel@pmTickMarkDisplayMode         = "Always"   ; Turn on map tickmarks
     opts_vvel@tmXBMajorLengthF              = 0.005      ; Change tickmark length
     opts_vvel@tmXTOn                        = "False"    ; No tickmarks on top x-axis
     opts_vvel@tmYROn                        = "False"    ; No tickmarks on right y-axis
     opts_vvel@gsnMaximize                   = True       ; Maximise plot size

;=============
; Plot data
;=============

    avo_plane = avo_plane * (10 ^ 4)
    avo_plane@units = "10~S~-4~N~ s~S~-1~N~"	

    avo_plane1 = avo_plane1 * (10 ^ 4)
    avo_plane1@units = "10~S~-4~N~ s~S~-1~N~"

    vort0    = gsn_csm_contour_map(wks,avo_plane,opts_vort)
    geo0     = gsn_csm_contour(wks,geo_plane,opts_geo)
    vvel0    = gsn_csm_contour(wks,w_plane,opts_vvel)
    overlay(vort0,geo0)
    overlay(vort0,vvel0)
    panel(0) = vort0

    vort1    = gsn_csm_contour_map(wks,avo_plane1,opts_vort)
    geo1     = gsn_csm_contour(wks,geo_plane1,opts_geo1)
    vvel1    = gsn_csm_contour(wks,w_plane1,opts_vvel)
    overlay(vort1,geo1)
    overlay(vort1,vvel1)
    panel(1) = vort1

    optsP                       = True
    optsP@gsnFrame              = False	   ; Do not advance the frame
    optsP@gsnPanelLabelBar      = True	   ; Turn off panel labelbar
;    optsP@lbTitleString         = p+" hPa absolute vorticity (10~S~-4~N~ s~S~-1~N~)"
    optsP@lbTitleString         = "Absolute vorticity (10~S~-4~N~ s~S~-1~N~)"
    optsP@lbTitlePosition	= "Bottom"
    optsP@lbTitleFontHeightF	= 0.015
    optsP@pmLabelBarWidthF      = 0.6
    optsP@pmLabelBarHeightF     = 0.15
    optsP@lbLabelFontHeightF    = 0.0125
    optsP@txString              = ""
    optsP@gsnPanelFigureStrings = ""
;    optsP@gsnPanelFigureStrings = (/"a) "+p+" hPa vorticity and geopotential height",\
;    				    "b) "+p2+" hPa vorticity and geopotential height"/)
    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.90
    optsP@gsnPanelBottom        = 0.10
    optsP@amJust                = "TopLeft"
    optsP@gsnPanelFigureStringsFontHeightF = 0.0125 ; Reduce size of labels (0.01)

    gsn_panel(wks,panel,(/1,2/),optsP) ; Draw as a single plot
    frame(wks)

     ct = ct + 1	; Counter variable (time)
  
   ; Tidy up before next iteration
     delete([/lat,lon,vort_plane,avo_plane,geo_plane,u_plane,v_plane,t_plane/])
     delete([/opts_geo@cnLevels,opts_geo1@cnLevels/])

   end do      ; End time loop (do it = 0, 3)

  end do        ; End ensemble member loop (do nf = 0, numINPUT-1)

end 