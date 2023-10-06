; Script to compute the flux form of the vorticity equation 
; and to calculate changes in circulation due to its different terms  

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

 pre = "$sam/um/cp/20160704T0000Z_PHI_4p4_L80_ra1t_pb"
 cp_input = (/pre+"012.nc",pre+"018.nc",pre+"024.nc",\
 	      pre+"030.nc",pre+"036.nc",pre+"042.nc",\
	      pre+"048.nc",pre+"054.nc",pre+"060.nc",\
	      pre+"066.nc",pre+"072.nc",pre+"078.nc",\
	      pre+"084.nc",pre+"090.nc",pre+"096.nc"/)

 pre1 = "$sam/um/gl/20160704T0000Z_glm_pb" ; 54 h forecast missing
 gl_input = (/pre1+"012.nc",pre1+"018.nc",pre1+"024.nc",\
              pre1+"030.nc",pre1+"036.nc",pre1+"042.nc",\
              pre1+"048.nc",pre1+"056.nc",pre1+"060.nc",\
              pre1+"066.nc",pre1+"072.nc",pre1+"078.nc",\
              pre1+"084.nc",pre1+"090.nc",pre1+"096.nc"/)

 pre2 = "$sam/ecm/ecmwf_2016070"
 ecm_input = (/pre2+"412.grib",pre2+"418.grib",pre2+"500.grib",\
               pre2+"506.grib",pre2+"512.grib",pre2+"518.grib",\
               pre2+"600.grib",pre2+"606.grib",pre2+"612.grib",\
               pre2+"618.grib",pre2+"700.grib",pre2+"706.grib",\
               pre2+"712.grib",pre2+"718.grib",pre2+"800.grib"/)

 pre3 = "$sam/era5/era5_2016070"
 era_input = (/pre3+"4_12.nc",pre3+"4_18.nc",pre3+"5_00.nc",\
               pre3+"5_06.nc",pre3+"5_12.nc",pre3+"5_18.nc",\
               pre3+"6_00.nc",pre3+"6_06.nc",pre3+"6_12.nc",\
               pre3+"6_18.nc",pre3+"7_00.nc",pre3+"7_06.nc",\
               pre3+"7_12.nc",pre3+"7_18.nc",pre3+"8_00.nc"/)

 numINPUT  = dimsizes(era_input)       ; Number of input files
 time_arr  = new(numINPUT,string)      ; Array to hold strings for second loop
 title_arr = new(numINPUT,string)      ; Second array to hold strings for plot titles

;================================================
; Create colourmap based on Brewer pattern
;================================================

 CBR_wet = (/ "(/1.0000,    1.0000,    1.0000   /)", \
              "(/0.0000,    0.0000,    0.0000   /)", \
              "(/0.7250,    0.7250,    0.7250   /)", \
              "(/0.8644,    0.8644,    0.8644   /)", \
              "(/0.9686,    0.9882,    0.9412   /)", \
              "(/0.8784,    0.9529,    0.8588   /)", \
              "(/0.8000,    0.9216,    0.7725   /)", \
              "(/0.6588,    0.8667,    0.7098   /)", \
              "(/0.4824,    0.8000,    0.7686   /)", \
              "(/0.3059,    0.7020,    0.8275   /)", \
              "(/0.1686,    0.5490,    0.7451   /)", \
              "(/0.0314,    0.4078,    0.6745   /)", \
              "(/0.0314,    0.2510,    0.5059   /)", \
              "(/0.6140,    0.5986,    0.7817   /)", \
              "(/0.4272,    0.3407,    0.6512   /)", \
              "(/0.2703,    0.0432,    0.5101   /)" /)

;===============================
; Start loop over input files 
;===============================

   do ifil = 0, numINPUT-1
   c = addfile(cp_input(ifil),"r")
   d = addfile(gl_input(ifil),"r")
   e = addfile(ecm_input(ifil),"r")
   f = addfile(era_input(ifil),"r")

   setvalues NhlGetWorkspaceObjectId
     "wsMaximumSize" : 300000000
   end setvalues

;=======================================================
; Get the variables we will need, from all 4 datasets
;=======================================================

;=============================
; Convection-permitting UM 
;=============================

; 2 (time) x 19	(vertical) x 811 (latitude) x 1098 (longitude)

  time_cp  = c->t            ; Times in file (2 - every hour)
  times_cp = dimsizes(time_cp)

  lon1_cp  = c->longitude_1  ; longitude (1098 points --> 109.04 to 152.92 degrees E)
  lat1_cp  = c->latitude_1   ; latitude  (810 points --> 1.8 to 34.16 degrees N)
  lon_cp   = c->longitude    ; longitude (1098 points --> 109.02 to 152.90 degrees E)
  lat_cp   = c->latitude     ; latitude (811 points --> 1.78 to 34.18 degrees N)

  pres_cp  = c->p            ; Pressure (levels)
  plevs_cp = dimsizes(pres_cp) ; Size of pressure level array

  u_cp     = c->u            ; Zonal velocity (m/s) [time | pres | lat | lon]
  v_cp     = c->v            ; Meridional velocity (m/s) [time | pres | lat | lon]
  vort_cp  = c->rvor         ; Relative vorticity (/s) [time | pres | lat | lon]

  t_cp     = c->temp         ; Temperature (K) [time | pres | lat1 | lon1]
  z_cp     = c->ht           ; Geopotential height (m) [time | pres | lat1 | lon1]

; Convert to units we want
  z_cp     = z_cp/10
  z_cp@units = "dam" ; Convert to decametres
  vort_cp  = vort_cp * (10 ^ -6)   ; Convert to /s
  vort_cp  = vort_cp * (10 ^ 5)    ; Convert to preferred units for plotting
  dx_cp    = lon_cp(1) - lon_cp(0) ; Grid spacing (longitude)
  dy_cp    = lat_cp(1) - lat_cp(0) ; Grid spacing (latitude)
  p_cp     = (/0,1,2,3,4,6,8,10,12,14,15,16,17/) ; Elements of 'pres_cp' to access later

;==============
; Global UM 
;==============

; 2 (time) x 19 (vertical) x 1153 (latitude) x 1536 (longitude)

  time_gl  = d->t            ; Times in file (2 - every hour)
  times_gl = dimsizes(time_gl)

  lon1_gl  = d->longitude_1  ; longitude (1536 points --> 0.12 to 359.88 degrees E)
  lat1_gl  = d->latitude_1   ; latitude  (1152 points --> -89.92 to 89.92 degrees N)
  lon_gl   = d->longitude    ; longitude (1536 points --> 0 to 359.77 degrees E)
  lat_gl   = d->latitude     ; latitude (1153 points --> -90 to 90 degrees N)

  pres_gl  = d->p            ; Pressure (levels)
  plevs_gl = dimsizes(pres_gl) ; Size of pressure level array

  u_gl     = d->u            ; Zonal velocity (m/s) [time | pres | lat | lon]
  v_gl     = d->v            ; Meridional velocity (m/s) [time | pres | lat | lon]
  vort_gl  = d->rvor         ; Relative vorticity (/s) [time | pres | lat | lon]

  t_gl     = d->temp         ; Temperature (K) [time | pres | lat1 | lon1]
  z_gl     = d->ht           ; Geopotential height (m) [time | pres | lat1 | lon1]

; Convert to units we want
  z_gl     = z_gl/10
  z_gl@units = "dam" ; Convert to decametres
  vort_gl  = vort_gl * (10 ^ -6) ; Convert to /s
  vort_gl  = vort_gl * (10 ^ 5)	 ; Convert to preferred units for plotting
  dx_gl    = lon_gl(1) - lon_gl(0) ; Grid spacing (longitude)
  dy_gl    = lat_gl(1) - lat_gl(0) ; Grid spacing (latitude)
  p_gl	   = (/0,1,2,3,4,6,8,10,12,14,15,16,17/) ; Elements of "pres_gl" to access later

;===================
; ECMWF analysis
;===================

 ; 25 (vertical) x 1441 (latitude) x 2800 (longitude)

   lon_ec   = e->g0_lon_2     ; longitude
   lat_ec   = e->g0_lat_1     ; latitude
   pres_ec  = e->lv_ISBL0     ; Pressure
   tk_ec    = e->T_GDS0_ISBL  ; Temperature (K)
   u_ec     = e->U_GDS0_ISBL  ; Zonal velocity (m/s)
   v_ec     = e->V_GDS0_ISBL  ; Meridional velocity (m/s)
   vort_ec  = e->VO_GDS0_ISBL ; Relative vorticity (/s)
   vort_ec  = vort_ec * (10 ^ 5) ; Convert to preferred units for plotting
   vort_ec  = vort_ec(::-1,:,:)	 ; Switch order of pressure levels in array

   lat_ec   = lat_ec(::-1)       ; Make latitude increase monotonically
   pres_ec  = pres_ec(::-1)	 ; Switch order of pressure levels (1000,975,950,...)
   plevs_ec = dimsizes(pres_ec)

   dy_ec    = lat_ec(1) - lat_ec(0) ; Grid spacing (lat)
   dx_ec    = lon_ec(1) - lon_ec(0) ; Grid spacing (lon)

   p_ec     = (/0,1,3,4,5,6,7,8,9,10,11,12,13/) ; Elements of "pres_ec" to access later

;===================
; ERA5 reanalysis
;===================

; 6 (time) x 25 (vertical) x 121 (latitude) x 161 (longitude)

   lon_era   = f->longitude      ; longitude (110 - 150 degrees E)
   lat_era   = f->latitude       ; latitude (5 - 35 degrees N)
   time_era  = f->time           ; time (hours since 1900-01-01 00:00:0.0)
   pres_era  = f->level          ; Pressure (hPa)
   z1_era    = f->z              ; Geopotential (m**2 s**-2)
   t1_era    = f->t              ; Temperature (K)
   u1_era    = f->u              ; Zonal velocity (m/s)
   v1_era    = f->v              ; Meridional velocity (m/s)
   vo1_era   = f->vo             ; Relative vorticity (/s)

   lat_era   = lat_era(::-1)     ; Make latitude increase monotonically
   pres_era  = pres_era(::-1)	 ; Switch order of pressure levels (1000,975,950,...)

   dy_era    = lat_era(1) - lat_era(0)   ; Grid spacing (lat)
   dx_era    = lon_era(1) - lon_era(0)   ; Grid spacing (lon)

   p_era     = (/0,2,4,6,8,11,13,15,17,19,20,22,24/) ; Elements of "pres_era" to access

;========================================================================
; Create correct date strings for each output time (using C-P UM data)
;========================================================================

   month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
                     "Jul","Aug","Sep","Oct","Nov","Dec"/)

   times = dimsizes(time_cp) ; Files are not all same size
   utc_date = cd_calendar(time_cp, 0)

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

;==========================================================
; Unpack the data using scale factors/offset and tidy up
;==========================================================

   z_era = tofloat((z1_era * z1_era@scale_factor) + z1_era@add_offset)
   copy_VarCoords(z1_era, z_era)
   copy_VarAtts(z1_era, z_era)
   copy_VarMeta(z1_era, z_era)
   delete(z_era@missing_value_original)
   delete(z_era@_FillValue_original)
   delete(z1_era)

   t_era = tofloat((t1_era * t1_era@scale_factor) + t1_era@add_offset)
   copy_VarCoords(t1_era, t_era)
   copy_VarAtts(t1_era, t_era)
   copy_VarMeta(t1_era, t_era)
   delete(t_era@missing_value_original)
   delete(t_era@_FillValue_original)
   delete(t1_era)
 
   u_era = tofloat((u1_era * u1_era@scale_factor) + u1_era@add_offset)
   copy_VarCoords(u1_era, u_era)
   copy_VarAtts(u1_era, u_era)
   copy_VarMeta(u1_era, u_era)
   delete(u_era@missing_value_original)
   delete(u_era@_FillValue_original)
   delete(u1_era)

   v_era    = tofloat((v1_era * v1_era@scale_factor) + v1_era@add_offset)
   copy_VarCoords(v1_era, v_era)
   copy_VarAtts(v1_era, v_era)
   copy_VarMeta(v1_era, v_era)
   delete(v_era@missing_value_original)
   delete(v_era@_FillValue_original)
   delete(v1_era)

   vort_era = tofloat((vo1_era * vo1_era@scale_factor) + vo1_era@add_offset)
   copy_VarCoords(vo1_era, vort_era)
   copy_VarAtts(vo1_era, vort_era)
   copy_VarMeta(vo1_era, vort_era)
   delete(vort_era@missing_value_original)
   delete(vort_era@_FillValue_original)
   delete(vo1_era)

   vort_era = vort_era * (10 ^ 5)  ; Convert to preferred units for plotting
   vort_era = vort_era(:,::-1,:,:) ; Switch order of pressure levels in array 

;===========================
; Additional calculations
;===========================

  gr = 9.80665  ; Gravitational acceleration
  gr@units = "m s~S~-2~N~"

  z_era  = z_era / (10 * gr) ; Calculate height from geopotential, convert to decametres
  z_era@units  = "dam"

  pii = 3.14159265 ; Degrees to radians

;=========================================
; Only interested in first time in file
;=========================================

     it = 0

     date_str(it) = sprinti("%0.2iUTC ", hour(it)) + \
                    sprinti("%0.2i ", day(it)) \
                    + month_abbr(month(it))
     time_str(it) = sprinti("%0.2i UTC ", hour(it)) + \
                    sprinti("%0.2i ", day(it)) \
                    + month_abbr(month(it))
     out_str(it) = sprinti("%0.2i", day(it)) + \
                   month_abbr(month(it)) + \
                   "_" + sprinti("%0.2iZ", hour(it))

     time_arr(ifil) = out_str(it)		   

     print("Working on time: "+time_str(it))
     title_arr(ifil) = "Valid at: "+time_str(it)

;==========================================================================
; Define new arrays for terms in vorticity equation (no friction for now)
;==========================================================================

  ; Create array subsets to speed up calculations
    lat_0 = 5
    lat_1 = 30
    lon_0 = 115
    lon_1 = 148

  ; Find cyclone centre using 950 hPa relative vorticity
    ilev_era = 2
    ilev_um  = 1 

    vort_pl_cp  = vort_cp(it,ilev_um,{lat_0:lat_1},{lon_0:lon_1})   ; C-P UM
    vort_pl_gl  = vort_gl(it,ilev_um,{lat_0:lat_1},{lon_0:lon_1})   ; Global UM
    vort_pl_ec  = vort_ec(ilev_um,{lat_0:lat_1},{lon_0:lon_1})      ; ECMWF analysis (3D)
    vort_pl_era = vort_era(it,ilev_era,{lat_0:lat_1},{lon_0:lon_1}) ; ERA5 reanalysis

    vort_max_cp  = max(vort_pl_cp)  ; Convection-permitting UM
    vort_max_gl	 = max(vort_pl_gl)  ; Global UM
    vort_max_ec	 = max(vort_pl_ec)  ; ECMWF analysis
    vort_max_er	 = max(vort_pl_era) ; ERA5 reanalysis

  ;========================================================
  ; Additionally, find centre using geopotential height
  ;========================================================

    geo_pl_cp	 = z_cp(it,ilev_um,{lat_0:lat_1},{lon_0:lon_1})
    geo_min_cp	 = min(geo_pl_cp)

  ;==============================================================
  ; Now calculate cyclone centre for each of the four datasets
  ;==============================================================

  ; Convection-permitting UM 
    dims_cp       = dimsizes(vort_pl_cp)
    vort1d_cp     = ndtooned(vort_pl_cp)
    inds_cp       = ind_resolve(maxind(vort1d_cp),dims_cp)
    lat_max950_cp = lat_0 + (dy_cp * inds_cp(0,0)) ; Max latitude
    lon_max950_cp = lon_0 + (dx_cp * inds_cp(0,1)) ; Max longitude

    dims_z_cp	  = dimsizes(geo_pl_cp)
    geo1d_cp	  = ndtooned(geo_pl_cp)
    inds_h_cp	  = ind_resolve(minind(geo1d_cp),dims_z_cp)
    lat_min950_cp = lat_0 + (dy_cp * inds_h_cp(0,0))
    lon_min950_cp = lon_0 + (dy_cp * inds_h_cp(0,1))

    print("Cyclone centre (C-P UM): "+lat_max950_cp+\
    	  " degrees N, "+lon_max950_cp+" degrees E")

    print("Cyclone centre (C-P UM - height): "+lat_min950_cp+\
    	  " degrees N, "+lon_min950_cp+" degrees E")	  

  ; Global UM
    dims_gl       = dimsizes(vort_pl_gl)
    vort1d_gl     = ndtooned(vort_pl_gl)
    inds_gl       = ind_resolve(maxind(vort1d_gl),dims_gl)
    lat_max950_gl = lat_0 + (dy_gl * inds_gl(0,0)) ; Max latitude
    lon_max950_gl = lon_0 + (dx_gl * inds_gl(0,1)) ; Max longitude
    print("Cyclone centre (Global UM): "+lat_max950_gl+\
          " degrees N, "+lon_max950_gl+" degrees E")

  ; ECMWF analysis
    dims_ec       = dimsizes(vort_pl_ec)
    vort1d_ec     = ndtooned(vort_pl_ec)
    inds_ec       = ind_resolve(maxind(vort1d_ec),dims_ec)
    lat_max950_ec = lat_0 + (dy_ec * inds_ec(0,0)) ; Max latitude
    lon_max950_ec = lon_0 + (dx_ec * inds_ec(0,1)) ; Max longitude
    print("Cyclone centre (ECMWF): "+lat_max950_ec+\
          " degrees N, "+lon_max950_ec+" degrees E")

  ; ERA5 reanalysis
    dims_era       = dimsizes(vort_pl_era)
    vort1d_era     = ndtooned(vort_pl_era)
    inds_era       = ind_resolve(maxind(vort1d_era),dims_era)
    lat_max950_era = lat_0 + (dy_era * inds_era(0,0)) ; Max latitude
    lon_max950_era = lon_0 + (dx_era * inds_era(0,1)) ; Max longitude
    print("Cyclone centre (ERA5): "+lat_max950_era+\
          " degrees N, "+lon_max950_era+" degrees E")

  ; Tidy up 
    delete([/vort_pl_cp,vort1d_cp,inds_cp,dims_cp/])
    delete([/vort_pl_gl,vort1d_gl,inds_gl,dims_gl/])
    delete([/vort_pl_ec,vort1d_ec,inds_ec,dims_ec/])
    delete([/vort_pl_era,vort1d_era,inds_era,dims_era/])

  ; Loop over levels (p_cp, p_gl, p_ec, p_era)

  do ilev = 0, 12   ; 13 pressure levels that match up in all 4 datasets

     p0 = p_cp(ilev)
     p1 = p_gl(ilev)
     p2 = p_ec(ilev)
     p3 = p_era(ilev)

     print("Working on: " + pres_cp(p0) + "hPa")

     vort_pl_cp   = vort_cp(it,p0,{lat_0:lat_1},{lon_0:lon_1})   ; C-P UM
     vort_pl_gl   = vort_gl(it,p1,{lat_0:lat_1},{lon_0:lon_1})   ; Global UM
     vort_pl_ec   = vort_ec(p2,{lat_0:lat_1},{lon_0:lon_1})      ; ECMWF
     vort_pl_era  = vort_era(it,p3,{lat_0:lat_1},{lon_0:lon_1})  ; ERA5

     vort_max_cp  = max(vort_pl_cp)  ; Convection-permitting UM
     vort_max_gl  = max(vort_pl_gl)  ; Global UM
     vort_max_ec  = max(vort_pl_ec)  ; ECMWF analysis
     vort_max_er  = max(vort_pl_era) ; ERA5 reanalysis

     ;============================
     ; Convection-permitting UM 
     ;============================

     dims_cp       = dimsizes(vort_pl_cp)
     vort1d_cp     = ndtooned(vort_pl_cp)
     inds_cp       = ind_resolve(maxind(vort1d_cp),dims_cp)
     lat_max_cp    = lat_0 + (dy_cp * inds_cp(0,0)) ; Max latitude
     lon_max_cp    = lon_0 + (dx_cp * inds_cp(0,1)) ; Max longitude
     print("Cyclone centre (C-P UM): "+lat_max_cp+\
           " degrees N, "+lon_max_cp+" degrees E")

     lat_pl_cp = ind(lat_0.le.lat_cp.and.lat_cp.le.lat_1) ; Indices where lat is in range
     lon_pl_cp = ind(lon_0.le.lon_cp.and.lon_cp.le.lon_1) ; Indices where lon is in range

     lat_sz_cp  = dimsizes(lat_pl_cp)
     lon_sz_cp  = dimsizes(lon_pl_cp)
     lat_sub_cp = lat_cp(lat_pl_cp(0):lat_pl_cp(lat_sz_cp-1))
     lon_sub_cp	= lon_cp(lon_pl_cp(0):lon_pl_cp(lon_sz_cp-1))

     delete([/vort1d_cp,inds_cp,dims_cp/]) ; Don't delete 'vort_pl_cp' this time

   ; Account for cases in which another maximum is chosen in error
     lat1_cp1 = abs(lat_max_cp-lat_min950_cp) ; abs(lat_max_cp-lat_max950_cp)
     lon1_cp1 = abs(lon_max_cp-lon_min950_cp) ; abs(lon_max_cp-lon_max950_cp)

     if (lon1_cp1.gt.2.0.or.lat1_cp1.gt.2.0) then
        print("Fixed an error in the tracking method!")
        lat_max_cp = lat_min950_cp ; lat_max950_cp
        lon_max_cp = lon_min950_cp ; lon_max950_cp
     end if
     print("Cyclone centre (C-P UM): "+lat_max_cp+\
           " degrees N, "+lon_max_cp+" degrees E")

     ;===============
     ; Global UM
     ;===============

     dims_gl       = dimsizes(vort_pl_gl)
     vort1d_gl     = ndtooned(vort_pl_gl)
     inds_gl       = ind_resolve(maxind(vort1d_gl),dims_gl)
     lat_max_gl    = lat_0 + (dy_gl * inds_gl(0,0)) ; Max latitude
     lon_max_gl    = lon_0 + (dx_gl * inds_gl(0,1)) ; Max longitude
     print("Cyclone centre (Global UM): "+lat_max_gl+\
           " degrees N, "+lon_max_gl+" degrees E")

     lat_pl_gl = ind(lat_0.le.lat_gl.and.lat_gl.le.lat_1) ; Indices where lat is in range
     lon_pl_gl = ind(lon_0.le.lon_gl.and.lon_gl.le.lon_1) ; Indices where lon is in range

     lat_sz_gl  = dimsizes(lat_pl_gl)
     lon_sz_gl  = dimsizes(lon_pl_gl)
     lat_sub_gl = lat_gl(lat_pl_gl(0):lat_pl_gl(lat_sz_gl-1))
     lon_sub_gl = lon_gl(lon_pl_gl(0):lon_pl_gl(lon_sz_gl-1))

     delete([/vort1d_gl,inds_gl,dims_gl/]) ; Don't delete 'vort_pl_gl' this time

   ; Account for cases in which another maximum is chosen in error
     lat1_gl1 = abs(lat_max_gl-lat_max950_gl)
     lon1_gl1 = abs(lon_max_gl-lon_max950_gl)

     if (lon1_gl1.gt.2.0.or.lat1_gl1.gt.2.0) then
        print("Fixed an error in the tracking method!")
        lat_max_gl = lat_max950_gl
        lon_max_gl = lon_max950_gl
     end if
     print("Cyclone centre (Global UM): "+lat_max_gl+\
     	   " degrees N, "+lon_max_gl+" degrees E")

     ;===================
     ; ECMWF analysis
     ;===================

     dims_ec       = dimsizes(vort_pl_ec)
     vort1d_ec     = ndtooned(vort_pl_ec)
     inds_ec       = ind_resolve(maxind(vort1d_ec),dims_ec)
     lat_max_ec    = lat_0 + (dy_ec * inds_ec(0,0)) ; Max latitude
     lon_max_ec    = lon_0 + (dx_ec * inds_ec(0,1)) ; Max longitude
     print("Cyclone centre (ECMWF): "+lat_max_ec+\
           " degrees N, "+lon_max_ec+" degrees E")

     lat_pl_ec = ind(lat_0.le.lat_ec.and.lat_ec.le.lat_1) ; Indices where lat is in range
     lon_pl_ec = ind(lon_0.le.lon_ec.and.lon_ec.le.lon_1) ; Indices where lon is in range

     lat_sz_ec  = dimsizes(lat_pl_ec)
     lon_sz_ec  = dimsizes(lon_pl_ec)
     lat_sub_ec = lat_ec(lat_pl_ec(0):lat_pl_ec(lat_sz_ec-1))
     lon_sub_ec = lon_ec(lon_pl_ec(0):lon_pl_ec(lon_sz_ec-1))

     delete([/vort1d_ec,inds_ec,dims_ec/]) ; Don't delete 'vort_pl_ec' this time

   ; Account for cases in which another maximum is chosen in error
     lat1_ec1 = abs(lat_max_ec-lat_max950_ec)
     lon1_ec1 = abs(lon_max_ec-lon_max950_ec)

     if (lon1_ec1.gt.2.0.or.lat1_ec1.gt.2.0) then
        print("Fixed an error in the tracking method!")
        lat_max_ec = lat_max950_ec
        lon_max_ec = lon_max950_ec
     end if
     print("Cyclone centre (ECMWF): "+lat_max_ec+\
           " degrees N, "+lon_max_ec+" degrees E")

     ;=====================
     ; ERA5 reanalysis
     ;=====================

     dims_era       = dimsizes(vort_pl_era)
     vort1d_era     = ndtooned(vort_pl_era)
     inds_era       = ind_resolve(maxind(vort1d_era),dims_era)
     lat_max_era    = lat_0 + (dy_era * inds_era(0,0)) ; Max latitude
     lon_max_era    = lon_0 + (dx_era * inds_era(0,1)) ; Max longitude
     print("Cyclone centre (ERA5): "+lat_max_era+\
           " degrees N, "+lon_max_era+" degrees E")

     lat_pl_era = ind(lat_0.le.lat_era.and.lat_era.le.lat_1)
     lon_pl_era = ind(lon_0.le.lon_era.and.lon_era.le.lon_1)

     lat_sz_era  = dimsizes(lat_pl_era)
     lon_sz_era  = dimsizes(lon_pl_era)
     lat_sub_era = lat_era(lat_pl_era(0):lat_pl_era(lat_sz_era-1))
     lon_sub_era = lon_era(lon_pl_era(0):lon_pl_era(lon_sz_era-1))

     delete([/vort1d_era,inds_era,dims_era/]) ; Don't delete 'vort_pl_era' this time

   ; Account for cases in which another maximum is chosen in error
     lat1_era1 = abs(lat_max_era-lat_max950_era)
     lon1_era1 = abs(lon_max_era-lon_max950_era)

     if (lon1_era1.gt.2.0.or.lat1_era1.gt.2.0) then
        print("Fixed an error in the tracking method!")
        lat_max_era = lat_max950_era
        lon_max_era = lon_max950_era
     end if
     print("Cyclone centre (ERA5): "+lat_max_era+\
     	   " degrees N, "+lon_max_era+" degrees E")

;=================================
; Output file type and location
;=================================

     diri = "$sam/nepartak/images/cyl_coords"
     output = diri+"/vort_panel_cyl_"+pres_cp(p0)+"_"+time_arr(ifil)
     wks = gsn_open_wks(opt,output)

;===========================================
; Switch to cylindrical coordinate system
;===========================================

     radii  = new(17,float)   ; New arrays
     thetas = new(73,float)

     radii      = fspan(0.0,2,17)       ; Array of radial distances (3 deg from centre)
     thetas     = fspan(0.0,360.0,73)   ; Array of azimuth angles (0, 5, 10, ..., 360)
     thetas_rad = (pii/180)*thetas      ; Degrees to radians
     rad_size   = dimsizes(radii)-1     ; Index for later calculation

   ; Create new arrays describing cylindrical coordinate system (2D for now)
     xpos_cp  = new((/dimsizes(radii),dimsizes(thetas)/),float)
     ypos_cp  = new((/dimsizes(radii),dimsizes(thetas)/),float)
     xpos_gl  = new((/dimsizes(radii),dimsizes(thetas)/),float)
     ypos_gl  = new((/dimsizes(radii),dimsizes(thetas)/),float)
     xpos_ec  = new((/dimsizes(radii),dimsizes(thetas)/),float)
     ypos_ec  = new((/dimsizes(radii),dimsizes(thetas)/),float)
     xpos_era = new((/dimsizes(radii),dimsizes(thetas)/),float)
     ypos_era = new((/dimsizes(radii),dimsizes(thetas)/),float)

   ; New arrays containing variables we want to plot
     vort_int_cp  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; C-P UM
     vort_int_gl  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Global UM
     vort_int_ec  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; ECMWF 
     vort_int_era = new((/dimsizes(radii),dimsizes(thetas)/),float) ; ERA5
 
   ; Latitude/longitude (small grid)
     do irad = 0, dimsizes(radii)-1
       r = radii(irad)
       do iang = 0, dimsizes(thetas)-1
         theta = thetas_rad(iang)
           xpos_cp(irad,iang)  = lon_max_cp + r*cos(theta)
           ypos_cp(irad,iang)  = lat_max_cp + r*sin(theta)
           xpos_gl(irad,iang)  = lon_max_gl + r*cos(theta)
           ypos_gl(irad,iang)  = lat_max_gl + r*sin(theta)
           xpos_ec(irad,iang)  = lon_max_ec + r*cos(theta)
           ypos_ec(irad,iang)  = lat_max_ec + r*sin(theta)
           xpos_era(irad,iang) = lon_max_era + r*cos(theta)
           ypos_era(irad,iang) = lat_max_era + r*sin(theta)
       end do
     end do

  ; Loop over azimuth angles to fill entire circular array of data points
    do iang = 0, dimsizes(thetas_rad)-1
      vort_int_cp(:,iang) = linint2_points(lon_sub_cp,lat_sub_cp,vort_pl_cp,False,\
                                           xpos_cp(:,iang),ypos_cp(:,iang),0)
      vort_int_gl(:,iang) = linint2_points(lon_sub_gl,lat_sub_gl,vort_pl_gl,False,\
                                           xpos_gl(:,iang),ypos_gl(:,iang),0)
      vort_int_ec(:,iang) = linint2_points(lon_sub_ec,lat_sub_ec,vort_pl_ec,False,\
                                           xpos_ec(:,iang),ypos_ec(:,iang),0)
      vort_int_era(:,iang) = linint2_points(lon_sub_era,lat_sub_era,vort_pl_era,False,\
                                           xpos_era(:,iang),ypos_era(:,iang),0)
    end do

   ; Add variable metadata
     copy_VarAtts(vort_pl_cp,vort_int_cp)
     vort_int_cp@description = "Relative vorticity"
     copy_VarAtts(vort_pl_gl,vort_int_gl)
     vort_int_gl@description = "Relative vorticity"
     copy_VarAtts(vort_pl_ec,vort_int_ec)
     vort_int_ec@description = "Relative vorticity"
     copy_VarAtts(vort_pl_era,vort_int_era)
     vort_int_era@description = "Relative vorticity"

   ; Position arrays (C-P UM)
     xpos_cp!0   = "rad"
     xpos_cp!1   = "azi"
     xpos_cp&rad = radii
     xpos_cp&azi = thetas

     ypos_cp!0   = "rad"
     ypos_cp!1   = "azi"
     ypos_cp&rad = radii
     ypos_cp&azi = thetas

   ; Position arrays (Global UM)
     xpos_gl!0   = "rad"
     xpos_gl!1   = "azi"
     xpos_gl&rad = radii
     xpos_gl&azi = thetas

     ypos_gl!0   = "rad"
     ypos_gl!1   = "azi"
     ypos_gl&rad = radii
     ypos_gl&azi = thetas

   ; Position arrays (ECMWF)
     xpos_ec!0   = "rad"
     xpos_ec!1   = "azi"
     xpos_ec&rad = radii
     xpos_ec&azi = thetas

     ypos_ec!0   = "rad"
     ypos_ec!1   = "azi"
     ypos_ec&rad = radii
     ypos_ec&azi = thetas

   ; Position arrays (ERA5)
     xpos_era!0   = "rad"
     xpos_era!1   = "azi"
     xpos_era&rad = radii
     xpos_era&azi = thetas

     ypos_era!0   = "rad"
     ypos_era!1   = "azi"
     ypos_era&rad = radii
     ypos_era&azi = thetas

   ; Vorticity interpolated to new arrays
     vort_int_cp!0   = "rad"
     vort_int_cp!1   = "azi"
     vort_int_cp&rad = radii
     vort_int_cp&azi = thetas

     vort_int_gl!0   = "rad"
     vort_int_gl!1   = "azi"
     vort_int_gl&rad = radii
     vort_int_gl&azi = thetas

     vort_int_ec!0   = "rad"
     vort_int_ec!1   = "azi"
     vort_int_ec&rad = radii
     vort_int_ec&azi = thetas

     vort_int_era!0   = "rad"
     vort_int_era!1   = "azi"
     vort_int_era&rad = radii
     vort_int_era&azi = thetas

;==========================
; Options for plotting 
;==========================

    panel = new(4,graphic)
    gsn_define_colormap(wks,CBR_wet)

  ; Relative vorticity
    opts_vort_cp       	 	     	    = True
    opts_vort_cp@sfXArray                   = xpos_cp
    opts_vort_cp@sfYArray                   = ypos_cp    
    opts_vort_cp@cnFillOn		    = True
    opts_vort_cp@cnLinesOn                  = False
    opts_vort_cp@cnLineLabelInterval	    = 2.0
    opts_vort_cp@cnLineLabelFontHeightF     = 0.012
    opts_vort_cp@cnLineLabelBackgroundColor = "transparent"
    opts_vort_cp@cnLineLabelPlacementMode   = "constant"
    opts_vort_cp@cnLinesOn                  = False ; Contour lines off
    opts_vort_cp@cnInfoLabelOn              = False
    opts_vort_cp@cnLevelSelectionMode       = "ExplicitLevels"
    opts_vort_cp@cnLevels                   = (/-2., 5., 10., 20., 30., \
                                                 40., 60., 80., 100., 150., \
                                                 200., 300./)
    opts_vort_cp@cnFillColors               = (/2,3,4,5,7,8,9,10,11,12,13,14,15/)
    opts_vort_cp@gsnPaperOrientation	    = "landscape"
    opts_vort_cp@tiMainString		    = ""
    opts_vort_cp@tiMainFontHeightF	    = 0.0125
    opts_vort_cp@gsnLeftString		    = ""
    opts_vort_cp@gsnRightString		    = ""
    opts_vort_cp@lbLabelBarOn		    = False

  ; Additional plotting resources 
    opts_vort_cp@mpDataBaseVersion	    = "Ncarg4_1"     ; More recent database 
    opts_vort_cp@mpDataSetName		    = "Earth..4"     ; High resolution 
    opts_vort_cp@mpOutlineBoundarySets	    = "National"     ; National borders
    opts_vort_cp@mpGeophysicalLineColor	    = "black"        ; Colour borders black
    opts_vort_cp@mpGeophysicalLineThicknessF = 1.0	     ; Border line thickness
    opts_vort_cp@mpGridAndLimbOn	    = False	     ; Turn on lat/lon lines
    opts_vort_cp@mpGridSpacingF		    = 2.0	     ; Spacing between grid lines
    opts_vort_cp@pmTickMarkDisplayMode	    = "Always"       ; Turn on map tickmarks
    opts_vort_cp@tmXBMajorLengthF	    = 0.005	     ; Change tickmark length
    opts_vort_cp@tmXBMinorOn		    = False	     ; Turn off minor tm (x-axis)
    opts_vort_cp@tmYLMinorOn		    = False	     ; Turn off minor tm (y-axis)
    opts_vort_cp@tmXTOn			    = "False"        ; Turn off tm (top x-axis)
    opts_vort_cp@tmYROn			    = "False"        ; Turn off tm (right y-axis)
    opts_vort_cp@gsnMaximize		    = False          ; Maximise plot size
    opts_vort_cp@gsnAddCyclic		    = False
    opts_vort_cp@gsnDraw                    = False	     ; Do not draw the plot
    opts_vort_cp@gsnFrame                   = False          ; Do not advance the frame
    opts_vort_cp@mpLimitMode                = "Corners"

    radius				    = radii(rad_size) ; Size of plot
    opts_vort_cp@mpLeftCornerLatF           = lat_max_cp-radius
    opts_vort_cp@mpLeftCornerLonF           = lon_max_cp-radius
    opts_vort_cp@mpRightCornerLatF          = lat_max_cp+radius
    opts_vort_cp@mpRightCornerLonF          = lon_max_cp+radius

  ; Global UM
    opts_vort_gl			    = opts_vort_cp
    delete(opts_vort_gl@sfXArray)
    delete(opts_vort_gl@sfYArray)
    opts_vort_gl@sfXArray                   = xpos_gl
    opts_vort_gl@sfYArray                   = ypos_gl
    opts_vort_gl@mpLeftCornerLatF           = lat_max_gl-radius
    opts_vort_gl@mpLeftCornerLonF           = lon_max_gl-radius
    opts_vort_gl@mpRightCornerLatF          = lat_max_gl+radius
    opts_vort_gl@mpRightCornerLonF          = lon_max_gl+radius

  ; ECMWF
    opts_vort_ec                            = opts_vort_gl
    delete(opts_vort_ec@sfXArray)
    delete(opts_vort_ec@sfYArray)
    opts_vort_ec@sfXArray                   = xpos_ec
    opts_vort_ec@sfYArray                   = ypos_ec
    opts_vort_ec@mpLeftCornerLatF           = lat_max_ec-radius
    opts_vort_ec@mpLeftCornerLonF           = lon_max_ec-radius
    opts_vort_ec@mpRightCornerLatF          = lat_max_ec+radius
    opts_vort_ec@mpRightCornerLonF          = lon_max_ec+radius

  ; ERA5
    opts_vort_era                           = opts_vort_ec
    delete(opts_vort_era@sfXArray)
    delete(opts_vort_era@sfYArray)
    opts_vort_era@sfXArray                  = xpos_era
    opts_vort_era@sfYArray                  = ypos_era
    opts_vort_era@mpLeftCornerLatF          = lat_max_era-radius
    opts_vort_era@mpLeftCornerLonF          = lon_max_era-radius
    opts_vort_era@mpRightCornerLatF         = lat_max_era+radius
    opts_vort_era@mpRightCornerLonF         = lon_max_era+radius

;=============
; Plot data 
;=============

    panel(0) = gsn_csm_contour_map(wks,vort_int_cp,opts_vort_cp)   ; C-P UM 
    panel(1) = gsn_csm_contour_map(wks,vort_int_gl,opts_vort_gl)   ; Global UM
    panel(2) = gsn_csm_contour_map(wks,vort_int_ec,opts_vort_ec)   ; ECMWF 
    panel(3) = gsn_csm_contour_map(wks,vort_int_era,opts_vort_era) ; ERA5

    opts_vortP                       = True
    opts_vortP@gsnFrame              = False
    opts_vortP@gsnPanelLabelBar      = True
    opts_vortP@lbTitleString         = "(10~S~-5~N~ s~S~-1~N~)"
    opts_vortP@lbTitlePosition       = "Bottom"
    opts_vortP@lbTitleFontHeightF    = 0.01
    opts_vortP@pmLabelBarWidthF      = 0.6
    opts_vortP@pmLabelBarHeightF     = 0.15
    opts_vortP@lbLabelFontHeightF    = 0.01
    opts_vortP@txString = title_arr(ifil)
    opts_vortP@gsnPanelFigureStrings = (/"a) Convection-permitting UM", \
    				         "b) Global UM", \
    				         "c) ECMWF analysis", \
                                         "d) ERA5 reanalysis"/)
    opts_vortP@gsnMaximize           = True
    opts_vortP@gsnPanelTop           = 0.95
    opts_vortP@gsnPanelBottom        = 0.05
    opts_vortP@amJust                = "TopLeft"
    opts_vortP@gsnPanelFigureStringsFontHeightF = 0.01 ; Reduce label size (def 0.01)

;==============================
; Add radial circles to plot
;==============================

      DEGTORAD = 0.017453292519943

    ; Create array to access and run following code multiple times 
      lon_arr = (/lon_max_cp, lon_max_gl, lon_max_ec, lon_max_era/)
      lat_arr = (/lat_max_cp, lat_max_gl, lat_max_ec, lat_max_era/)

    ; Create 4 * 73 element arrays to hold x,y values for each panel
      xc      = new((/dimsizes(lon_arr),dimsizes(thetas)/),float)
      yc      = new((/dimsizes(lon_arr),dimsizes(thetas)/),float)

    ; Create array to hold outer circle grid points
      plot_outer = new(dimsizes(lon_arr),graphic)

    ; Create array to hold inner circle grid points
      radius     = radii(rad_size)
      spacing    = 0.5
      size       = toint((radius-spacing)/spacing) 
      plot_inner = new((/dimsizes(lon_arr),size/),graphic)

    ; Also create array to hold radial lines 
      angle      = 45
      xlines  	 = 360./angle
      nlines  	 = toint(ceil(xlines))
      plot_lines = new((/dimsizes(lon_arr),nlines/),graphic)

   do j = 0, dimsizes(lon_arr)-1

    ; Get customizations for circles.
      radius  = radii(rad_size)
      xcenter = lon_arr(j)
      ycenter = lat_arr(j)
      outthck = 2
      inthck = 1
      spacing = 0.5
      dpattrn = 2
      outcolr = 1
      incolr = 1

    ; Calculate arrays for outer circle (similar to what I have done)
      degrees = ispan(0,360,5)
      xcos    = cos(DEGTORAD * degrees)
      xsin    = sin(DEGTORAD * degrees)

    ; Array of points (x,y) representing the outer circle (radius = 25)
      xc(j,:)    = xcenter + (radius * xcos)
      yc(j,:)    = ycenter + (radius * xsin)

    ; Resources for outer circle
      lnres                   = True
      lnres@gsLineThicknessF  = outthck
      lnres@gsLineColor       = outcolr
      lnres@gsLineDashPattern = 0

    ; Attach circle to plot
      plot_outer(j) = gsn_add_polyline(wks,panel(j),xc(j,:),yc(j,:),lnres)

      if(spacing.gt.radius) then
        print("add_radial_circles: spacing is > radius, can't draw inner circles.")
        return
      end if

    ; Draw inner circles if desired ('spacing' controls the number of inner circles)
      if(spacing.gt.0) then
        count = 0
        do r = spacing,radius-spacing,spacing

	; Calculate arrays for inner circle.
          x1 = xcenter + (r * xcos)
          y1 = ycenter + (r * xsin)

	  ; Resources for inner circle
          delete(lnres@gsLineThicknessF)   ; Delete in case setting to
          delete(lnres@gsLineColor)        ; new type.
          delete(lnres@gsLineDashPattern)

          lnres@gsLineColor       = incolr
          lnres@gsLineThicknessF  = inthck
          lnres@gsLineDashPattern = dpattrn

          ; Attach inner circle to plot
          plot_inner(j,count) = gsn_add_polyline(wks,panel(j),x1,y1,lnres)
          count = count + 1
        end do
        delete(lnres)
        delete(degrees)
      end if

;=====================
; Add radial lines
;=====================

    ; Get customizations for radial lines (radius, xcenter, ycenter are already defined)
      angle   = 45
      lcolor  = 1
      dpattrn = 2
      thcknss = 1.

    ; Error checking
      if(angle.le.0.or.angle.ge.360) then
        print("add_radial_lines: angle must be between 0 and 360.")
        print("                  Can't draw radial lines.")
        return
      end if

    ; Resources for radial lines
      lnres                   = True
      lnres@gsLineDashPattern = dpattrn
      lnres@gsLineColor       = lcolor
      lnres@gsLineThicknessF  = thcknss

    ; Calculate arrays for lines
      if(xlines.ne.nlines) then
        degreesL = fspan(0,360,nlines)
      else
        degreesL = fspan(0,360-angle,nlines) ; Contains angles (0, 45, 90, etc)
      end if

      nlines  = dimsizes(degreesL)
      do i = 0,nlines-1
        xl = (/xcenter, xcenter + (radius  * cos(DEGTORAD * degreesL(i)))/)
        yl = (/ycenter, ycenter + (radius  * sin(DEGTORAD * degreesL(i)))/)

      ; Attach line to plot
        plot_lines(j,i) = gsn_add_polyline(wks,panel(j),xl,yl,lnres)

      end do

      delete([/xcenter, ycenter/]) ; Tidy up after each iteration
    end do	; End loop over array objects (do i = 0, dimsizes(lon_arr)-1)

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    gsn_panel(wks,panel,(/2,2/),opts_vortP) ; Draw as a single plot
    frame(wks)

  ; Tidy up (delete variables from within pressure level loop before next iteration)
    delete([/vort_pl_cp, vort_pl_gl, vort_pl_ec, vort_pl_era/])
    delete([/vort_int_cp, vort_int_gl, vort_int_ec, vort_int_era/])
    delete([/xc, yc, plot_outer, plot_inner, plot_lines/])

    end do ; End pressure level loop (do ilev = 0, 12)

  ; Tidy up before next iteration
    delete([/out_str,time_str,date_str/])

  end do     ; End input file loop (do time = 0, size-2)

end 