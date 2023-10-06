; Script to calculate area-averaged quantities around a developing TC 

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

 ecm_input = (/"$sam/ecm/ecm_new_20160704.nc",\
               "$sam/ecm/ecm_new_20160705.nc",\
               "$sam/ecm/ecm_new_20160706.nc",\
               "$sam/ecm/ecm_new_20160707.nc",\
               "$sam/ecm/ecm_new_20160708.nc"/)

;===================================
; Prelash (before start of loops)
;===================================

   numINPUT  = dimsizes(ecm_input)    ; Number of input files
   FirstTime = True                   ; Skip first time when calculating vort tendencies
   r_earth   = 6.37122e6              ; Earth's radius (constant)
   numTIMES  = numINPUT * 4           ; Number of time intervals across ALL input files
   ct        = 0                      ; Counter variable

   circ_inst = new((/14,numTIMES/),float) ; Array to hold circulation values (14 * 28)
   rh_ave    = new((/14,numTIMES/),float) ; Relative humidity (14 * 28)
   t_pert    = new((/14,numTIMES/),float) ; Temperature perturbation (14 * 28)
   t_init    = new((/14/),float)	  ; Initial temperature
   time_arr  = new(numTIMES,string)       ; Array to hold strings for second loop
   title_arr = new(numTIMES,string)       ; Second array to hold strings for plot titles
   pos_arr   = new((/numTIMES,2/),float)  ; Array for lat/lon values of cyclone centre

;================================
; Start loop over input files
;================================

 do ifil = 0, numINPUT-1
 f = addfile(ecm_input(ifil),"r")

 setvalues NhlGetWorkspaceObjectId
   "wsMaximumSize" : 300000000
 end setvalues

;==================================
; Get the variables we will need
;==================================

; 14 (vertical) x 241 (latitude) x 481 (longitude) 

  lon   = f->longitude     ; longitude   	 [100 - 160 degrees E]
  lat   = f->latitude      ; latitude 	         [0   - 30  degrees S]
  lat_size = dimsizes(lat)  
  lon_size = dimsizes(lon)

  pres  = f->level         ; Pressure            [14 levels --> 150 - 1000 hPa]
  time	= f->time	   ; Time 	         [4 times --> 00, 06, 12, 18] 
  pv1	= f->pv		   ; Potential vorticity [PVU]
  z1	= f->z		   ; Geopotential        [m**2 / s**2]
  t1    = f->t             ; Temperature 	 [K]     
  w1    = f->w             ; Vertical velocity   [Pa/s]
  q1	= f->q		   ; Specific humidity   [kg / kg]
  rh1	= f->r		   ; Relative humidity   [%]
  u1    = f->u             ; Zonal velocity      [m / s]
  v1    = f->v             ; Meridional velocity [m / s]
  vort1 = f->vo            ; Relative vorticity  [ / s]
  div1  = f->d             ; Divergence		 [ / s] 

  lat   = lat(::-1)       ; Make latitude increase monotonically
  pres@units = "hPa"
  plevs = dimsizes(pres)

  dx    = lon(1) - lon(0) ; Grid spacing (longitude)
  dy    = lat(1) - lat(0) ; Grid spacing (latitude)  

  circ_inst!0    = "lev"
  circ_inst&lev  = pres

  rh_ave!0	 = "lev"
  rh_ave&lev	 = pres

  t_pert!0	 = "lev"
  t_pert&lev	 = pres

;==========================================================
; Unpack the data using scale factors/offset and tidy up
;==========================================================

 pv    = tofloat((pv1 * pv1@scale_factor) + pv1@add_offset)
 copy_VarCoords(pv1, pv)
 copy_VarAtts(pv1, pv)
 copy_VarMeta(pv1, pv)
 delete(pv@missing_value_original)
 delete(pv@_FillValue_original)
 delete(pv1)

 q    = tofloat((q1 * q1@scale_factor) + q1@add_offset)
 copy_VarCoords(q1, q)
 copy_VarAtts(q1, q)
 copy_VarMeta(q1, q)
 delete(q@missing_value_original)
 delete(q@_FillValue_original)
 delete(q1)

 rh   = tofloat((rh1 * rh1@scale_factor) + rh1@add_offset)
 copy_VarCoords(rh1, rh)
 copy_VarAtts(rh1, rh)
 copy_VarMeta(rh1, rh)
 delete(rh@missing_value_original)
 delete(rh@_FillValue_original)
 delete(rh1)

 z    = tofloat((z1 * z1@scale_factor) + z1@add_offset)
 copy_VarCoords(z1, z)
 copy_VarAtts(z1, z)
 copy_VarMeta(z1, z)
 delete(z@missing_value_original)
 delete(z@_FillValue_original)
 delete(z1)

 gr   = 9.80665
 z    = z / (10 * gr)
 z@units = "dam"

 t    = tofloat((t1 * t1@scale_factor) + t1@add_offset)
 copy_VarCoords(t1, t)
 copy_VarAtts(t1, t)
 copy_VarMeta(t1, t)
 delete(t@missing_value_original)
 delete(t@_FillValue_original)
 delete(t1)

 u    = tofloat((u1 * u1@scale_factor) + u1@add_offset)
 copy_VarCoords(u1, u)
 copy_VarAtts(u1, u)
 copy_VarMeta(u1, u)
 delete(u@missing_value_original)
 delete(u@_FillValue_original)
 delete(u1)

 v    = tofloat((v1 * v1@scale_factor) + v1@add_offset)
 copy_VarCoords(v1, v)
 copy_VarAtts(v1, v)
 copy_VarMeta(v1, v)
 delete(v@missing_value_original)
 delete(v@_FillValue_original)
 delete(v1)

 w    = tofloat((w1 * w1@scale_factor) + w1@add_offset)
 copy_VarCoords(w1, w)
 copy_VarAtts(w1, w)
 copy_VarMeta(w1, w)
 delete(w@missing_value_original)
 delete(w@_FillValue_original)
 delete(w1)

 vort = tofloat((vort1 * vort1@scale_factor) + vort1@add_offset)
 copy_VarCoords(vort1, vort)
 copy_VarAtts(vort1, vort)
 copy_VarMeta(vort1, vort)
 delete(vort@missing_value_original)
 delete(vort@_FillValue_original)
 delete(vort1)

 div  = tofloat((div1 * div1@scale_factor) + div1@add_offset)
 copy_VarCoords(div1, div)
 copy_VarAtts(div1, div)
 copy_VarMeta(div1, div)
 delete(div@missing_value_original)
 delete(div@_FillValue_original)
 delete(div1)

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
   do it = 0, times-1
     if (minute(it).gt.30) then
       hour(it) = hour(it)+1
     end if
   end do

   date_str = new(times,string)
   time_str = new(times,string)
   out_str  = new(times,string)

;======================================================================
; Make latitude arrays increase monotonically (for grad_latlon_cfd)
;======================================================================

  z    = z(:,:,::-1,:)
  t    = t(:,:,::-1,:)
  u    = u(:,:,::-1,:)
  v    = v(:,:,::-1,:)
  w    = w(:,:,::-1,:)
  vort = vort(:,:,::-1,:)
  div  = div(:,:,::-1,:)
  pv   = pv(:,:,::-1,:)
  q    = q(:,:,::-1,:)
  rh   = rh(:,:,::-1,:)

;====================================================
; Pressure levels to loop over later in the script
;====================================================

 ; Pressure levels
 ; (0)  150
 ; (1)  200
 ; (2)  250
 ; (3)  300
 ; (4)  400
 ; (5)  500
 ; (6)  600
 ; (7)  700
 ; (8)  800
 ; (9)  850
 ; (10) 900
 ; (11) 925
 ; (12) 950
 ; (13) 1000

 ; Calculate absolute vorticity
   d = dimsizes(vort)
   f0  = new((/d(0),d(1),d(2),d(3)/),float)  ; Earth's vorticity
   avo = new((/d(0),d(1),d(2),d(3)/),float)  ; Absolute vorticity 

 ; Add metadata
   copy_VarMeta(vort,avo)
   copy_VarAtts(vort,avo)
   copy_VarCoords(vort,avo)
  
 ; Degrees to radians (for calculation below)
   pii = 3.14159265
   lat_rad = (pii/180)*lat
   lon_rad = (pii/180)*lon     

 ; Calculate the Coriolis parameter at all grid points 
   do ilon = 0, d(3) - 1
      do ilat = 0, d(2) - 1
         f0(:,:,ilat,ilon) = 2 * (7.27 * 10 ^ -5) * sin(lat_rad(ilat))
      end do 
   end do  

   avo = vort + f0   ; Absolute vorticity 
   copy_VarMeta(vort,avo)

   avo@long_name = "Vorticity (planetary)"
   avo@standard_name = "atmosphere_planetary_vorticity"
   avo@title     = "Absolute vorticity"
   avo@name      = "avo"

;===========================
; Loop over times in file
;===========================

   do it = 0, times-1   ; Loop over all times in the file

     date_str(it) = sprinti("%0.2iUTC ", hour(it)) + \
                    sprinti("%0.2i ", day(it)) \
                    + month_abbr(month(it))
     time_str(it) = sprinti("%0.2i UTC ", hour(it)) + \
                    sprinti("%0.2i ", day(it)) \
                    + month_abbr(month(it))
     out_str(it) = sprinti("%0.2i", day(it)) + \
                   month_abbr(month(it)) + \
                   "_" + sprinti("%0.2iZ", hour(it))

     time_arr(ct) = out_str(it)

     print("Working on time: "+time_str(it))
     title_arr(ct) = "Valid at: "+time_str(it)

;========================================
; Define new arrays for advective flux
;========================================

  ; Create array subset (not in use atm, but potentially later)
    lat_0 = 0   ; -15 ; -30
    lat_1 = 30  ; -5  ;  0
    lon_0 = 115 ; 115 ; 100
    lon_1 = 150 ; 135 ; 160

    u_sub     = u(it,:,{lat_0:lat_1},{lon_0:lon_1})
    v_sub     = v(it,:,{lat_0:lat_1},{lon_0:lon_1})
    omega_sub = w(it,:,{lat_0:lat_1},{lon_0:lon_1})
    vort_sub  = vort(it,:,{lat_0:lat_1},{lon_0:lon_1})
    avo_sub   = avo(it,:,{lat_0:lat_1},{lon_0:lon_1})
    pv_sub    = pv(it,:,{lat_0:lat_1},{lon_0:lon_1})
    div_sub   = div(it,:,{lat_0:lat_1},{lon_0:lon_1})
    geo_sub   = z(it,:,{lat_0:lat_1},{lon_0:lon_1})
    t_sub     = t(it,:,{lat_0:lat_1},{lon_0:lon_1})
    init_sub  = t(0,:,{lat_0:lat_1},{lon_0:lon_1})
    q_sub     = q(it,:,{lat_0:lat_1},{lon_0:lon_1})
    rh_sub    = rh(it,:,{lat_0:lat_1},{lon_0:lon_1})

;    avo_sub   = -(avo_sub)	   ; -ve absolute vorticity corresponds to cyclonic flow

    d_sub     = dimsizes(u_sub)
    lat_sub   = lat({lat_0:lat_1})    ; Subset of lat array
    lon_sub   = lon({lon_0:lon_1})    ; Subset of lon array
    lat_subR  = lat_sub * (pii/180)   ; Convert to radians for calculation below
    lon_subR  = lon_sub * (pii/180)

;==============================================================
; Define latitude and longitude (following Juliane's script)
;==============================================================

     y2d     = new((/d_sub(1),d_sub(2)/),float)
     y2d!0   = "lat"
     y2d!1   = "lon"
     y2d&lat =  lat_sub
     y2d&lon =  lon_sub

     y2d     = conform_dims(dimsizes(y2d),lat_sub,0) ; Array of latitude values
     coslat  = cos(y2d*pii/180)                      ; Same array in radians

   ; Determine horizontal resolution of data
     delta_lon = y2d
     delta_lat = y2d
     delta_lon = lon_sub(1) - lon_sub(0) ; Zonal grid spacing      (0.25 degrees)
     delta_lat = lat_sub(1) - lat_sub(0) ; Meridional grid spacing (0.25 degrees)
     delta_lon = sqrt(delta_lon^2)       ; Still 0.25 degrees (array)
     delta_lat = sqrt(delta_lat^2)       ; Still 0.25 degrees (array)

     lons = lon_sub(1) - lon_sub(0)      ; 0.25 degrees (single value)
     lats = lat_sub(1) - lat_sub(0)      ; 0.25 degrees (single value)
     lons = sqrt(lons^2)
     lats = sqrt(lats^2)

   ; Calculate dx and dy (dy is constant; dx changes with latitude)
     ddx = y2d
     ddx = r_earth * ((pii/180.0) * delta_lon) * cos(y2d*pii/180)
     ddy = y2d
     ddy = r_earth * ((pii/180) * delta_lat)

   ; 3-D latitude array
     coslat3d = conform_dims(dimsizes(u_sub),coslat,(/1,2/))

;========================
; Find cyclone centre
;========================

  ; Find cyclone centre using 850 hPa geopotential height
    ilev = 9

    geo_plane    = geo_sub(ilev,:,:)
    geo_min      = min(geo_plane)
    geo_smth     = smth9_Wrap(geo_plane, 0.5, 0.5, True)
    geo_min_smth = min(geo_smth)

    dims         = dimsizes(geo_plane)
    geo1d        = ndtooned(geo_plane)
    inds         = ind_resolve(minind(geo1d),dims)
    geo1d_sm     = ndtooned(geo_smth)
    inds_sm      = ind_resolve(minind(geo1d_sm),dims)

    lat_min = lat_0 + (dy * inds(0,0)) ; Subset starts at 5 degrees N
    lon_min = lon_0 + (dx * inds(0,1)) ; Subset starts at 105 degrees E
    print("Cyclone centre: "+lat_min+" degrees N, "+lon_min+" degrees E")

    delete([/geo1d, inds, geo1d_sm, geo_smth, inds_sm, geo_plane/])

;==============================
; Loop over pressure levels 
;==============================

   cp = 0	         ; Counter variable (pressure)

   do ilev = 0, plevs-1  ; 800, 850, 900 hPa

     p = pres(ilev)
     print("Working on: " + p + "hPa")

   ; Calculate vort max. on each pressure level
     vort_plane    = vort_sub(ilev,:,:)
     avo_plane	   = avo_sub(ilev,:,:)
     geo_plane     = geo_sub(ilev,:,:)
     u_plane	   = u_sub(ilev,:,:)
     v_plane	   = v_sub(ilev,:,:)
     rh_plane	   = rh_sub(ilev,:,:)
     init_plane	   = init_sub(ilev,:,:)
     t_plane	   = t_sub(ilev,:,:)

     zmin          = toint(min(geo_plane)) ; Minimum geopotential on pressure level
     zmax          = toint(max(geo_plane)) ; Maximum geopotential on pressure level

;==========================================================
; Calculate variables around a 2 degree by 2 degree box 
;==========================================================

   ; Radius of box for calculation (degrees)
     r = 1.0
     r1 = sprintf("%0.1fdeg",2*r) ; Create string to add to output file

   ; Design box around the centre of the TC 
     ddx_box    = ddx({lat_min-r:lat_min+r},{lon_min-r:lon_min+r})
     ddy_box    = ddy({lat_min-r:lat_min+r},{lon_min-r:lon_min+r})
     u_box      = u_plane({lat_min-r:lat_min+r},{lon_min-r:lon_min+r})
     v_box      = v_plane({lat_min-r:lat_min+r},{lon_min-r:lon_min+r})
     rh_box	= rh_plane({lat_min-r:lat_min+r},{lon_min-r:lon_min+r})
     init_box	= init_plane({lat_min-r:lat_min+r},{lon_min-r:lon_min+r})
     t_box	= t_plane({lat_min-r:lat_min+r},{lon_min-r:lon_min+r})

   ; Convert velocity to km/s (for circulation budget calculations)
     u_box       = u_box / 1000
     u_box@units = "km s~S~-1~N~"
     v_box       = v_box / 1000
     v_box@units = "km s~S~-1~N~"

   ; Convert grid spacing to kilometres
     ddx_box = ddx_box / 1000
     ddy_box = ddy_box / 1000
     ddx_box@units = "km"
     ddy_box@units = "km"

   ; Calculate distance between four corner points of the 2 x 2 degree box
     lat_d = lat_min-r
     lon_d = lon_min-r
     lat_u = lat_min+r
     lon_u = lon_min+r

   ;=====================================================
   ; Calculate distances around each perimeter in km
   ;=====================================================

   ; Meridional boundaries (same both sides)
     dm0 = r_earth * ((pii/180) * (2*r))
     dm0 = dm0 / 1000
     dm0@units = "km"

   ; 'Bottom' boundary
     dz0 = r_earth * ((pii/180) * (2*r)) * cos((pii/180)*lat_d)
     dz0 = dz0 / 1000
     dz0@units = "km"

   ; 'Top' boundary
     dz1 = r_earth * ((pii/180) * (2*r)) * cos((pii/180)*lat_u)
     dz1 = dz1 / 1000
     dz1@units = "km"

   ; Use approximation to calculate the area of the box
     d_ave = ((dz0 + dz1) / 2) * dm0
     d_ave@units = "km~S~-2~N~"

   ; Calculate relative humidity averaged over the box
     rh_aveT = wgt_areaave(rh_box, 1.0, 1.0, 0) ; Relative humidity
     rh_aveT@units = "%"

     initT   = wgt_areaave(init_box, 1.0, 1.0, 0) ; Initial temperature
     tempT   = wgt_areaave(t_box, 1.0, 1.0, 0)    ; Temperature
     t_pertT = initT - tempT	       	    	  ; Perturbation temperature
     t_pertT@units = "K"

   ;=============================================================
   ; Now calculate circulation around the perimeter of the box
   ;=============================================================

   ; Number of grid points along each side of box surrounding cyclone
     numPOINTS_x0 = dimsizes(ddx_box)
     numPOINTS_x  = numPOINTS_x0(1)   ; Longitude points
     numPOINTS_y  = numPOINTS_x0(0)   ; Latitude points

   ; New array to hold integral sum (circulation)
     circ_x_final0 = new(numPOINTS_x,float)
     circ_x_final1 = new(numPOINTS_x,float)
     circ_y_final0 = new(numPOINTS_y,float)
     circ_y_final1 = new(numPOINTS_y,float)

   ; Bottom side
     do ilon = 0, numPOINTS_x-1
       circ_x_final0(ilon) = u_box(0,ilon) * ddx_box(0,ilon) ; Circulation
     end do

   ; Right side
     do ilat = 0, numPOINTS_y-1
       circ_y_final0(ilat) = v_box(ilat,numPOINTS_x-1) * \
                             ddy_box(ilat,numPOINTS_x-1)   ; Circulation
     end do

   ; Top side
     do ilon = 0, numPOINTS_x-1
       circ_x_final1(ilon) = u_box(numPOINTS_y-1,ilon) * \
       			     ddx_box(numPOINTS_y-1,ilon)   ; Circulation
     end do
     circ_x_final1 = -(circ_x_final1) ; Easterly winds contribute to +ve circulation

   ; Left side
     do ilat = 0, numPOINTS_y-1
       circ_y_final1(ilat) = v_box(ilat,0) * ddy_box(ilat,0) ; Circulation
     end do
     circ_y_final1 = -(circ_y_final1) ; Northerly winds contribute to +ve circulation

   ; Calculate circulation (LHS)
     circ_instT = sum(circ_x_final0) + sum(circ_x_final1) + \
                  sum(circ_y_final0) + sum(circ_y_final1) 

     circ_inst(cp,ct)     = circ_instT   ; Circulation (instantaneous)
     rh_ave(cp,ct)	  = rh_aveT	 ; Relative humidity
     t_pert(cp,ct)	  = t_pertT	 ; Perturbation temperature
     cp                   = cp + 1       ; Counter variable (pressure)  

     delete([/numPOINTS_x0, numPOINTS_x, numPOINTS_y/])
     delete([/circ_x_final0, circ_x_final1, circ_y_final0, circ_y_final1/])
     delete([/ddx_box, ddy_box, u_box, v_box, rh_box, init_box, t_box/])

    end do       ; End pressure level loop (do ilev = 0, plevs-1)

    ct = ct + 1	 ; Counter variable (time)

   end do        ; End time loop (do it = 0, times-1)

 ; Tidy up
   delete([/lon,lat,pres,time,pv,z,t,w,q,r,u,v,vort,div/])
   delete([/time_str, date_str, out_str/])

   end do	 ; End input file loop (do ifil = 0, numINPUT-1)

   ; Change sign so that a positive change is cyclonic (see Narsey et al. 2017)
;     circ_inst = -(circ_inst)

;=============================================================
; Open workstation for plotting (time-height cross section)
;=============================================================

     output = "$sam/nepartak/images/cart_coords/fig8"
     wks = gsn_open_wks(opt, output)

;=========================
; Options for plotting 
;=========================

    panel = new(3,graphic)

    gsn_define_colormap(wks,"GMT_polar")

  ; Circulation
    opts_circ                             = True
    opts_circ@cnFillOn                    = True
    opts_circ@cnLineLabelInterval         = 2.0
    opts_circ@cnLineLabelFontHeightF      = 0.012
    opts_circ@cnLineLabelBackgroundColor  = "transparent"
    opts_circ@cnLineLabelPlacementMode    = "constant"
    opts_circ@cnLinesOn                   = False
    opts_circ@cnInfoLabelOn               = False
    opts_circ@cnLevelSelectionMode        = "ExplicitLevels"
    opts_circ@cnLevels                    = (/-1.0, 0.0, 3.0, 6.0, 9.0/)
    opts_circ@cnFillColors                = (/9,11,0,13,16,20/)
    opts_circ@gsnPaperOrientation         = "landscape"
    opts_circ@tiMainFontHeightF           = 0.0125
    opts_circ@gsnLeftString               = ""
    opts_circ@gsnRightString              = ""

    opts_circ@tmXTOn                      = "False" ; Turn off top x-axis TM
    opts_circ@tmYROn                      = "False" ; Turn off right y-axis TM
    opts_circ@vpWidthF			  = 0.8     ; Plot width
    opts_circ@vpHeightF			  = 0.4	    ; Plot height
    opts_circ@gsnMaximize                 = True    ; Maximise plot size
    opts_circ@trXMaxF			  = 24	    ; Number of 6 h time intervals 
    opts_circ@gsnDraw			  = False   ; Do not draw the plot
    opts_circ@gsnFrame			  = False   ; Do not advance the frame	
    opts_circ@lbLabelBarOn		  = True    ; Keep labelbar on

  ; Perturbation temperature
    opts_temp				  = opts_circ
    opts_temp@cnLevels			  = (/-1.0, -0.5, 0.0, 0.5, 1.0/)

  ; Relative humidity
    opts_rh				  = opts_circ
    opts_rh@cnLevels			  = (/20, 40, 60, 80, 85/)

;=============
; Plot data
;=============

    panel(0) = gsn_csm_pres_hgt(wks,t_pert,opts_temp)     ; Perturbation temperature
    panel(1) = gsn_csm_pres_hgt(wks,rh_ave,opts_rh)       ; Relative humidity
    panel(2) = gsn_csm_pres_hgt(wks,circ_inst,opts_circ)  ; Circulation

    optsP                       = True
    optsP@gsnFrame              = False	   ; Do not advance the frame
    optsP@gsnPanelLabelBar      = False	   ; Turn off panel labelbar
    optsP@txString              = "Tropical cyclone Narelle: area-averaged quantities"
    optsP@gsnPanelBottom        = 0.05
    optsP@gsnPanelFigureStrings = (/"a) Perturbation temperature", \
                                    "b) Relative humidity", \
                                    "c) Circulation"/)
    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.9
    optsP@gsnPanelBottom        = 0.1
    optsP@amJust                = "TopLeft"
    optsP@gsnPanelFigureStringsFontHeightF = 0.008 ; Reduce size of labels (0.01)

    gsn_panel(wks,panel,(/3,1/),optsP) ; Draw as a single plot
    frame(wks)

end 