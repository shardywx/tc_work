; Script to compute the flux form of the vorticity equation 
; and to calculate changes in circulation due to its different terms  

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"

begin

;=========================================================
; Calculate storm motion [for system-relative velocity]
;=========================================================

  ; List all files to read in and analyse [00Z 4-9 July]
    diri        = "$sam/era5/cdf/"
    fili_prefix = "era5_201607"
    fili_a      = systemfunc("cd "+diri+" ; ls "+fili_prefix+"*.nc")
    fili        = diri+fili_a
    numINPUT    = dimsizes(fili)            ; Number of input files [20]
    numTIMES    = numINPUT * 6              ; Total number of times in all files [120]

    centre	= new((/numTIMES,2/),float) ; Cyclone centre at each time interval
    lat_arr	= new((/numTIMES/),float)   ; Latitude array
    lon_arr	= new((/numTIMES/),float)   ; Longitude array
    ct		= 0			    ; Counter variable

    time_arr    = new(numTIMES,string)      ; Array to hold strings for second loop
    title_arr   = new(numTIMES,string)      ; Array to hold strings for plot titles
    pos_arr     = new((/numTIMES,2/),float) ; Array for lat/lon values of cyclone centre
    r_earth     = 6.37e6                    ; Earth's radius (constant)

    setvalues NhlGetWorkspaceObjectId
      "wsMaximumSize" : 1000000000
    end setvalues

;==========================
; Loop over input files
;==========================

  do nf = 0,  numINPUT-1     ; 00Z 4 -- 9 July
  f = addfile(fili(nf),"r")

;==================================
; Get the variables we will need
;==================================

 ; 6 (time) x 25 (vertical) x 121 (latitude) x 181 (longitude)
   lon   = f->longitude      ; longitude (110 - 150 degrees E)
   lat   = f->latitude       ; latitude (5 - 35 degrees N)
   time  = f->time           ; time (hours since 1900-01-01 00:00:0.0)
   pres  = f->level          ; Pressure (hPa)
   plevs = dimsizes(pres)
   pres@units = "hPa"
   z1    = f->z              ; Geopotential (m**2 s**-2)
   t1    = f->t              ; Temperature (K)
   u1    = f->u              ; Zonal velocity (m/s)
   v1    = f->v              ; Meridional velocity (m/s)
   w1    = f->w              ; Vertical velocity (Pa/s)
   vo1   = f->vo             ; Relative vorticity (/s)
   div1  = f->d              ; Divergence (/s)
   lat   = lat(::-1)         ; Make latitude increase monotonically
   dy    = lat(1) - lat(0)   ; Grid spacing (lat)
   dx    = lon(1) - lon(0)   ; Grid spacing (lon)

; Pressure levels (pres)
; (0) 150, (1) 175, (2) 200, (3) 225, (4) 250, (5) 300, (6) 350, (7) 400
; (8) 450, (9) 500, (10) 550, (11) 600, (12) 650, (13) 700
; (14) 750, (15) 775, (16) 800, (17) 825, (18) 850
; (19) 875, (20) 900, (21) 925, (22) 950, (23) 975, (24) 1000

;==========================================================
; Unpack the data using scale factors/offset and tidy up
;==========================================================

   z    = tofloat((z1 * z1@scale_factor) + z1@add_offset)
   copy_VarCoords(z1, z)
   copy_VarAtts(z1, z)
   copy_VarMeta(z1, z)
   delete([/z@missing_value_original, z@_FillValue_original, z1/])

   t    = tofloat((t1 * t1@scale_factor) + t1@add_offset)
   copy_VarCoords(t1, t)
   copy_VarAtts(t1, t)
   copy_VarMeta(t1, t)
   delete([/t@missing_value_original, t@_FillValue_original, t1/])

   u    = tofloat((u1 * u1@scale_factor) + u1@add_offset)
   copy_VarCoords(u1, u)
   copy_VarAtts(u1, u)
   copy_VarMeta(u1, u)
   delete([/u@missing_value_original, u@_FillValue_original, u1/])

   v    = tofloat((v1 * v1@scale_factor) + v1@add_offset)
   copy_VarCoords(v1, v)
   copy_VarAtts(v1, v)
   copy_VarMeta(v1, v)
   delete([/v@missing_value_original, v@_FillValue_original, v1/])

   w    = tofloat((w1 * w1@scale_factor) + w1@add_offset)
   copy_VarCoords(w1, w)
   copy_VarAtts(w1, w)
   copy_VarMeta(w1, w)
   delete([/w@missing_value_original, w@_FillValue_original, w1/])

   vort = tofloat((vo1 * vo1@scale_factor) + vo1@add_offset)
   copy_VarCoords(vo1, vort)
   copy_VarAtts(vo1, vort)
   copy_VarMeta(vo1, vort)
   delete([/vort@missing_value_original, vort@_FillValue_original, vo1/])

   div  = tofloat((div1 * div1@scale_factor) + div1@add_offset)
   copy_VarCoords(div1, div)
   copy_VarAtts(div1, div)
   copy_VarMeta(div1, div)
   delete([/div@missing_value_original, div@_FillValue_original, div1/])

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

;===========================
; Additional calculations
;===========================

   gr = 9.80665  ; Gravitational acceleration
   gr@units = "m s~S~-2~N~"

   z  = z / (10 * gr)    ; Calculate height from geopotential and convert to decametres
   z@units  = "dam"

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

;===========================
; Loop over times in file
;===========================

   do it = 0, times-1

     date_str(it) = sprinti("%0.2iUTC ", hour(it)) + \
                    sprinti("%0.2i ", day(it)) \
                    + month_abbr(month(it))
     time_str(it) = sprinti("%0.2i UTC ", hour(it)) + \
                    sprinti("%0.2i ", day(it)) \
                    + month_abbr(month(it))
     out_str(it)  = sprinti("%0.2i", day(it)) + \
                    month_abbr(month(it)) + \
                    "_" + sprinti("%0.2iZ", hour(it))

     time_arr(ct) = out_str(it)

     print_clock("Working on time: "+time_str(it))
     title_arr(ct) = "Valid at: "+time_str(it)

;================================================================================
; Find cyclone centre on subset of global grid (use both height and vorticity)
;================================================================================

   ; Grid subset to include all entire cylcone track [81 * 82 grid points]
     lat_0 = 0.00
     lat_1 = 31.00
     lon_0 = 110.00
     lon_1 = 150.00

   ; Choose a level in the lower troposphere (950 hPa)
     ilev = 22

   ; First, concentrate on a subset of the global grid
     u_plane    = u(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
     v_plane    = v(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
     geo_plane  = z(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
     t_plane	= t(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
     vort_plane = vort(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
     dsize	= dimsizes(vort_plane)

     vort_max      = max(vort_plane)
     vort_smth     = smth9_Wrap(vort_plane, 0.5, 0.5, True)
     vort_max_smth = max(vort_smth)

     dims          = dimsizes(vort_plane)
     vort1d        = ndtooned(vort_plane)
     inds          = ind_resolve(maxind(vort1d),dims)
     vort1d_sm     = ndtooned(vort_smth)
     inds_sm       = ind_resolve(maxind(vort1d_sm),dims)

     lat_max950 = lat_0 + (dy * inds(0,0)) ; Latitude of max. vorticity
     lon_max950 = lon_0 + (dx * inds(0,1)) ; Longitude of max. vorticity
     print("Cyclone centre (vort): "+lat_max950+" degrees N, "+lon_max950+" degrees E")

   ; GEOPOTENTIAL HEIGHT
     geo_min        = min(geo_plane)
     geo_smth       = smth9_Wrap(geo_plane, 0.5, 0.5, True)
     geo_min_smth   = min(geo_smth)

     dims_h         = dimsizes(geo_plane)
     geo1d          = ndtooned(geo_plane)
     inds_h         = ind_resolve(minind(geo1d),dims_h)
     geo1d_sm       = ndtooned(geo_smth)
     inds_h_sm      = ind_resolve(minind(geo1d_sm),dims_h)

     lat_min950 = lat_0 + (dy * inds_h(0,0))
     lon_min950 = lon_0 + (dx * inds_h(0,1))
     print("Cyclone centre (hgt): "+lat_min950+" degrees N, "+lon_min950+" degrees E")

   ;========================================================
   ; EMPLOY A SAFETY NET IF WE IDENTIFY THE WRONG CYCLONE
   ;========================================================
     if (ct.eq.0) then
;      centre(ct,0) = lat_min950
;      centre(ct,1) = lon_min950
      centre(ct,0) = lat_max950
      centre(ct,1) = lon_max950

     else
      lt0   = centre(ct-1,0)            ; Previous TC latitude
      ln0   = centre(ct-1,1)            ; Previous TC longitude
      print("lt0 = "+lt0+" ; ln0 = "+ln0)
;      d_lat = abs(lat_min950-lt0)       ; Change in latitude (t1 - t0)
;      d_lon = abs(lon_min950-ln0)       ; Change in longitude (note extra minus sign)
      d_lat = abs(lat_max950-lt0)       ; Change in latitude (t1 - t0)
      d_lon = abs(lon_max950-ln0)       ; Change in longitude (note extra minus sign)
      print("d_lat = "+d_lat+" ; d_lon = "+d_lon)

      if (d_lat.gt.1.or.d_lon.gt.1) then        ; Recalculate TC centre if incorrect

     ; Create smaller grid [0.5 degrees]
       inc = 0.5
       lt1 = lt0 - inc
       lt2 = lt0 + inc
       ln1 = ln0 - inc
       ln2 = ln0 + inc

       geo_plane0 = z(it,ilev,{lt1:lt2},{ln1:ln2})       ; Smaller grid (1.0 x 1.0 deg)
       geo_min0   = min(geo_plane0)                      ; Find minimum on smaller grid
       dims_h0    = dimsizes(geo_plane0)                 ; Size of smaller grid
       geo_1d0    = ndtooned(geo_plane0)                 ; Create 1-D array
       inds_h0    = ind_resolve(minind(geo_1d0),dims_h0) ; Find index of minimum

       vort_plane0 = vort(it,ilev,{lt1:lt2},{ln1:ln2})     ; Smaller grid (1.0 x 1.0 deg)
       vort_max0   = max(vort_plane0)                      ; Find minimum on smaller grid
       dims_v0     = dimsizes(vort_plane0)                 ; Size of smaller grid
       vort_1d0    = ndtooned(vort_plane0)                 ; Create 1-D array
       inds_v0     = ind_resolve(minind(vort_1d0),dims_v0) ; Find index of minimum

       lat_min950 = lt1 + (dy * inds_h0(0,0))
       lon_min950 = ln1 + (dy * inds_h0(0,1))
       lat_max950 = lt1 + (dy * inds_v0(0,0))
       lon_max950 = ln1 + (dy * inds_v0(0,1))
       print("New centre (geo): "+lat_min950+" degrees N, "+lon_min950+" degrees E")
       print("New centre (vort): "+lat_max950+" degrees N, "+lon_max950+" degrees E")
       centre(ct,0) = lat_max950 ; lat_min950
       centre(ct,1) = lon_max950 ; lon_min950

       delete([/geo_plane0, geo_min0, dims_h0, geo_1d0, inds_h0/])
       delete([/vort_plane0, vort_max0, dims_v0, vort_1d0, inds_v0/])

      else
       print("Cyclone centre (geo): "+lat_min950+" degrees N, "+lon_min950+" degrees E")
       print("Cyclone centre (vort): "+lat_max950+" degrees N, "+lon_max950+" degrees E")
       centre(ct,0) = lat_max950 ; lat_min950
       centre(ct,1) = lon_max950 ; lon_min950
      end if

     end if

    ; Tidy up before looping over pressure levels
      delete([/vort1d, inds, vort1d_sm, vort_smth, inds_sm/])
      delete([/geo1d, inds_h, geo1d_sm, geo_smth, inds_h_sm/])
      delete([/u_plane, v_plane, vort_plane, geo_plane/])

      lat_arr(ct) = centre(ct,0)
      lon_arr(ct) = centre(ct,1)

      ct = ct + 1

    end do    ; End of time loop (do it = 0, 1)

  ; Tidy up to avoid array dimension errors [but leave alone on final iteration]
    if (nf.lt.numINPUT-1)
      delete([/time,times,pres,z,u,v,w,t,vort,div/])
      delete([/utc_date,year,month,day,hour,\
               minute,second,date_str,time_str,out_str/])
    end if

   end do     ; End of input file loop (do nf = 0, numINPUT-1) 

   diri = "$sam/ncl_scripts/nepartak/archer_march2018/tc_era5_"
   lat_out = diri+"lat.txt"
   lon_out = diri+"lon.txt"
   asciiwrite(lat_out,lat_arr)
   asciiwrite(lon_out,lon_arr)

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

 ; Loop to calculate system velocity from hourly cyclone position [lat,lon]
   do iang = 1, numTIMES-1
    x_cyc(iang)   = gc_latlon(centre(iang-1,0),centre(iang-1,1),centre(iang-1,0),\
                              centre(iang,1),2,3)
    y_cyc(iang)   = gc_latlon(centre(iang-1,0),centre(iang-1,1),centre(iang,0),\
                              centre(iang-1,1),2,3)
    u_cyc(iang)   = x_cyc(iang)/(3600)
    v_cyc(iang)   = y_cyc(iang)/(3600)
    vel_cyc(iang) = sqrt(u_cyc(iang)^2.0 + v_cyc(iang)^2.0)
   end do

 ; Set initial system velocity [0] equal to that at the second time [1]
   u_cyc(0)   = u_cyc(1)
   v_cyc(0)   = v_cyc(1)
   vel_cyc(0) = vel_cyc(1)

 ; Add metadata to new arrays
   copy_VarAtts(u,u_cyc)
   u_cyc@units         = "m s~S~-1~N~"
   u_cyc@description   = "Cyclone zonal velocity"

   copy_VarAtts(v,v_cyc)
   v_cyc@units         = "m s~S~-1~N~"
   v_cyc@description   = "Cyclone meridional velocity"

   vel_cyc@units       = "m s~S~-1~N~"
   vel_cyc@description = "Cyclone speed"

   u_cyc = -(u_cyc) ; Easterly winds are negative

 ; Tidy up
   delete([/lat,lon,pres,u,v,w,z,t,vort,div/])
   delete([/time,times,utc_date,year,month,day,hour,\
            minute,second,date_str,time_str,out_str/])

;==================================================================================
;;;;;;;;;;;;;;;;;;;;;;;;;;; Start of main part of script ;;;;;;;;;;;;;;;;;;;;;;;;;;
;==================================================================================

    print_clock("Starting main part of script...")

  ; List all files to read in and analyse
  ; Simulation runs from 00Z 4 to 9 July [120 h]
    diri        = "$sam/era5/cdf/"
    fili_prefix = "era5_201607"
    fili_a      = systemfunc("cd "+diri+" ; ls "+fili_prefix+"*.nc")
    fili        = diri+fili_a
    numINPUT    = dimsizes(fili)            ; Number of input files [20]
    numTIMES    = numINPUT * 6              ; Total number of times in all files [120]

    ct	        = 0		            ; Counter variable
    time_arr    = new(numTIMES,string)	    ; Array to hold strings for second loop
    title_arr   = new(numTIMES,string)	    ; Array to hold strings for plot titles
    pos_arr     = new((/numTIMES,2/),float) ; Array for lat/lon values of cyclone centre
    r_earth     = 6.37e6                    ; Earth's radius (constant)

  ; Create array to hold all values of relative vorticity, to compute tendencies
    circ_size   = (/plevs,dsize(0),dsize(1)/) ; [plevs * lat_dims * lon_dims]
    circ_arr    = new((/numTIMES,circ_size(0),circ_size(1),circ_size(2)/),float)
    circ_arr!0  = "time"         ; Add coordinate variable to dimension 0

    setvalues NhlGetWorkspaceObjectId
      "wsMaximumSize" : 1000000000
    end setvalues

  ; Create arrays to hold circulation tendency values [ens (25) * time (120) * levs (25)]
    circ_inst = new((/25,numTIMES,plevs/),float)

;===============================
; Start loop over input files 
;===============================

  do nf = 0,  numINPUT-1     ; 00Z 4 -- 9 July
  f = addfile(fili(nf),"r")

;==================================
; Get the variables we will need
;==================================

 ; 6 (time) x 25 (vertical) x 121 (latitude) x 181 (longitude)
   lon   = f->longitude      ; longitude (110 - 150 degrees E)
   lat   = f->latitude       ; latitude (5 - 35 degrees N)
   time  = f->time           ; time (hours since 1900-01-01 00:00:0.0)
   pres  = f->level          ; Pressure (hPa)
   z1    = f->z              ; Geopotential (m**2 s**-2)
   t1    = f->t              ; Temperature (K)
   u1    = f->u              ; Zonal velocity (m/s)
   v1    = f->v              ; Meridional velocity (m/s)
   w1    = f->w              ; Vertical velocity (Pa/s)
   vo1   = f->vo             ; Relative vorticity (/s)
   div1  = f->d              ; Divergence (/s)
   lat   = lat(::-1)         ; Make latitude increase monotonically
   dy    = lat(1) - lat(0)   ; Grid spacing (lat)
   dx    = lon(1) - lon(0)   ; Grid spacing (lon)

; Pressure levels (pres)
; (0) 150, (1) 175, (2) 200, (3) 225, (4) 250, (5) 300, (6) 350, (7) 400
; (8) 450, (9) 500, (10) 550, (11) 600, (12) 650, (13) 700
; (14) 750, (15) 775, (16) 800, (17) 825, (18) 850
; (19) 875, (20) 900, (21) 925, (22) 950, (23) 975, (24) 1000

;==========================================================
; Unpack the data using scale factors/offset and tidy up
;==========================================================

   z    = tofloat((z1 * z1@scale_factor) + z1@add_offset)
   copy_VarCoords(z1, z)
   copy_VarAtts(z1, z)
   copy_VarMeta(z1, z)
   delete([/z@missing_value_original, z@_FillValue_original, z1/])

   t    = tofloat((t1 * t1@scale_factor) + t1@add_offset)
   copy_VarCoords(t1, t)
   copy_VarAtts(t1, t)
   copy_VarMeta(t1, t)
   delete([/t@missing_value_original, t@_FillValue_original, t1/])

   u    = tofloat((u1 * u1@scale_factor) + u1@add_offset)
   copy_VarCoords(u1, u)
   copy_VarAtts(u1, u)
   copy_VarMeta(u1, u)
   delete([/u@missing_value_original, u@_FillValue_original, u1/])

   v    = tofloat((v1 * v1@scale_factor) + v1@add_offset)
   copy_VarCoords(v1, v)
   copy_VarAtts(v1, v)
   copy_VarMeta(v1, v)
   delete([/v@missing_value_original, v@_FillValue_original, v1/])

   w    = tofloat((w1 * w1@scale_factor) + w1@add_offset)
   copy_VarCoords(w1, w)
   copy_VarAtts(w1, w)
   copy_VarMeta(w1, w)
   delete([/w@missing_value_original, w@_FillValue_original, w1/])

   vort = tofloat((vo1 * vo1@scale_factor) + vo1@add_offset)
   copy_VarCoords(vo1, vort)
   copy_VarAtts(vo1, vort)
   copy_VarMeta(vo1, vort)
   delete([/vort@missing_value_original, vort@_FillValue_original, vo1/])

   div  = tofloat((div1 * div1@scale_factor) + div1@add_offset)
   copy_VarCoords(div1, div)
   copy_VarAtts(div1, div)
   copy_VarMeta(div1, div)
   delete([/div@missing_value_original, div@_FillValue_original, div1/])

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

;===========================
; Additional calculations
;===========================

   gr = 9.80665  ; Gravitational acceleration
   gr@units = "m s~S~-2~N~"

   z  = z / (10 * gr)    ; Calculate height from geopotential and convert to decametres
   z@units  = "dam"

;========================================
; Also calculation absolute vorticity
;========================================

; Degrees to radians (for calculation below)
  pii = 3.14159265
  lat_rad = (pii/180)*lat    ; u,v grid
  lon_rad = (pii/180)*lon     

; Calculate absolute vorticity
  d = dimsizes(vort)
  f0  = new((/d(0),d(1),d(2),d(3)/),float)  ; Earth's vorticity
  avo = new((/d(0),d(1),d(2),d(3)/),float)  ; Absolute vorticity

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
  copy_VarMeta(vort,f0)

  avo           = vort + f0 ; Absolute vorticity
  copy_VarMeta(vort,avo)

  avo@long_name = "Vorticity (planetary)"
  avo@standard_name = "atmosphere_planetary_vorticity"
  avo@title     = "Absolute vorticity"
  avo@name      = "avo"

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

;===========================
; Loop over times in file
;===========================

   do it = 0, times-1

     date_str(it) = sprinti("%0.2iUTC ", hour(it)) + \
                    sprinti("%0.2i ", day(it)) \
                    + month_abbr(month(it))
     time_str(it) = sprinti("%0.2i UTC ", hour(it)) + \
                    sprinti("%0.2i ", day(it)) \
                    + month_abbr(month(it))
     out_str(it)  = sprinti("%0.2i", day(it)) + \
                    month_abbr(month(it)) + \
                    "_" + sprinti("%0.2iZ", hour(it))

     time_arr(ct) = out_str(it)

     print_clock("Working on time: "+time_str(it))
     title_arr(ct) = "Valid at: "+time_str(it)

;==========================================================================
; Define new arrays for terms in vorticity equation (no friction for now)
;==========================================================================

  ; Grid subset to include all entire cylcone track [dsize(0) * dsize(1) grid points]
    lat_0 = 0.00
    lat_1 = 31.00
    lon_0 = 110.00
    lon_1 = 150.00

    geo_sub    = z(it,:,{lat_0:lat_1},{lon_0:lon_1})
    tk_sub     = t(it,:,{lat_0:lat_1},{lon_0:lon_1})
    f0_sub     = f0(it,:,{lat_0:lat_1},{lon_0:lon_1})
    omega_sub  = w(it,:,{lat_0:lat_1},{lon_0:lon_1})
    vort_sub   = vort(it,:,{lat_0:lat_1},{lon_0:lon_1})
    avo_sub    = avo(it,:,{lat_0:lat_1},{lon_0:lon_1})
    div_sub    = div(it,:,{lat_0:lat_1},{lon_0:lon_1})
    u_sub      = u(it,:,{lat_0:lat_1},{lon_0:lon_1})
    v_sub      = v(it,:,{lat_0:lat_1},{lon_0:lon_1})

  ; Storm-relative velocity
    u_rel      = u_sub - u_cyc(ct)
    v_rel      = v_sub - v_cyc(ct)

  ; Add metadata
    avo_sub@units  = "s~S~-1~N~"
    f0_sub@units   = "s~S~-1~N~"

    d_sub     = dimsizes(vort_sub)    ; Size of arrays
    lat_sub   = lat({lat_0:lat_1})    ; Subset of lat array
    lon_sub   = lon({lon_0:lon_1})    ; Subset of lon array
    lat_subR  = lat_sub * (pii/180)   ; Convert to radians for calculation below
    lon_subR  = lon_sub * (pii/180)

   ; Troubleshooting (look for abnormally large values and exit the script if found)
     do p1 = 0, d_sub(0)-1    ; pressure
       do l1 = 0, d_sub(1)-1   ; latitude
         do l2 = 0, d_sub(2)-1 ; longitude
           if (u_rel(p1,l1,l2) .gt. 200)
            print("Quitting because of storm-relative velocity error!!!")
             exit()
           end if
         end do
       end do
     end do

    copy_VarMeta(u_sub,u_rel)
    copy_VarAtts(u_sub,u_rel)
    copy_VarCoords(u_sub,u_rel)

    copy_VarMeta(v_sub,v_rel)
    copy_VarAtts(v_sub,v_rel)
    copy_VarCoords(v_sub,v_rel)

;==========================================================
; Define latitude and longitude for finite differencing
;==========================================================

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
    delta_lon = lon_sub(1) - lon_sub(0) ; Zonal grid spacing
    delta_lat = lat_sub(1) - lat_sub(0) ; Meridional grid spacing
    delta_lon = sqrt(delta_lon^2)
    delta_lat = sqrt(delta_lat^2)

    lons = lon_sub(1) - lon_sub(0)
    lats = lat_sub(1) - lat_sub(0)
    lons = sqrt(lons^2)
    lats = sqrt(lats^2)

  ; Calculate dx and dy (changes with longitude but not with latitude)
    ddx = y2d
    ddx = r_earth*(pii/180.0)*delta_lon*cos(y2d*pii/180)
    ddy = y2d
    ddy = r_earth*(pii/180)*delta_lat

  ; 3-D latitude array
    coslat3d = conform_dims(dimsizes(u_sub),coslat,(/1,2/))

;================================================================================
; Find cyclone centre on subset of global grid (use both height and vorticity)
;================================================================================

   ; Choose a level in the lower troposphere (950 hPa)
     ilev = 22

   ; First, concentrate on a subset of the global grid 
     geo_plane     = geo_sub(ilev,{lat_0:lat_1},{lon_0:lon_1})
     vort_plane    = vort_sub(ilev,{lat_0:lat_1},{lon_0:lon_1})

   ; RELATIVE VORTICITY
     vort_max      = max(vort_plane)
     vort_smth     = smth9_Wrap(vort_plane, 0.5, 0.5, True)
     vort_max_smth = max(vort_smth)

     dims          = dimsizes(vort_plane)
     vort1d        = ndtooned(vort_plane)
     inds          = ind_resolve(maxind(vort1d),dims)
     vort1d_sm     = ndtooned(vort_smth)
     inds_sm       = ind_resolve(maxind(vort1d_sm),dims)

     lat_max950 = lat_0 + (dy * inds(0,0)) ; Latitude of max. vorticity
     lon_max950 = lon_0 + (dx * inds(0,1)) ; Longitude of max. vorticity
     print("Cyclone centre (vort): "+lat_max950+" degrees N, "+lon_max950+" degrees E")

   ; GEOPOTENTIAL HEIGHT 
     geo_min	    = min(geo_plane)
     geo_smth	    = smth9_Wrap(geo_plane, 0.5, 0.5, True)
     geo_min_smth   = min(geo_smth)
      
     dims_h	    = dimsizes(geo_plane)
     geo1d	    = ndtooned(geo_plane)
     inds_h	    = ind_resolve(minind(geo1d),dims_h)
     geo1d_sm	    = ndtooned(geo_smth)
     inds_h_sm	    = ind_resolve(minind(geo1d_sm),dims_h)            

     lat_min950 = lat_0 + (dy * inds_h(0,0))
     lon_min950 = lon_0 + (dx * inds_h(0,1))
     print("Cyclone centre (hgt): "+lat_min950+" degrees N, "+lon_min950+" degrees E")

   ;========================================================
   ; EMPLOY A SAFETY NET IF WE IDENTIFY THE WRONG CYCLONE
   ;========================================================
     if (ct.eq.0) then
;      centre(ct,0) = lat_min950
;      centre(ct,1) = lon_min950
      centre(ct,0) = lat_max950
      centre(ct,1) = lon_max950

     else
      lt0   = centre(ct-1,0)            ; Previous TC latitude
      ln0   = centre(ct-1,1)            ; Previous TC longitude
      print("lt0 = "+lt0+" ; ln0 = "+ln0)
;      d_lat = abs(lat_min950-lt0)       ; Change in latitude (t1 - t0)
;      d_lon = abs(lon_min950-ln0)       ; Change in longitude (note extra minus sign)
      d_lat = abs(lat_max950-lt0)       ; Change in latitude (t1 - t0)
      d_lon = abs(lon_max950-ln0)       ; Change in longitude (note extra minus sign)
      print("d_lat = "+d_lat+" ; d_lon = "+d_lon)

      if (d_lat.gt.1.or.d_lon.gt.1) then        ; Recalculate TC centre if incorrect

     ; Create smaller grid [0.5 degrees]
       inc = 0.5
       lt1 = lt0 - inc
       lt2 = lt0 + inc
       ln1 = ln0 - inc
       ln2 = ln0 + inc

       geo_plane0 = z(it,ilev,{lt1:lt2},{ln1:ln2})       ; Smaller grid (1.0 x 1.0 deg)
       geo_min0   = min(geo_plane0)                      ; Find minimum on smaller grid
       dims_h0    = dimsizes(geo_plane0)                 ; Size of smaller grid
       geo_1d0    = ndtooned(geo_plane0)                 ; Create 1-D array
       inds_h0    = ind_resolve(minind(geo_1d0),dims_h0) ; Find index of minimum

       vort_plane0 = vort(it,ilev,{lt1:lt2},{ln1:ln2})     ; Smaller grid (1.0 x 1.0 deg)
       vort_max0   = max(vort_plane0)                      ; Find minimum on smaller grid
       dims_v0     = dimsizes(vort_plane0)                 ; Size of smaller grid
       vort_1d0    = ndtooned(vort_plane0)                 ; Create 1-D array
       inds_v0     = ind_resolve(minind(vort_1d0),dims_v0) ; Find index of minimum

       lat_min950 = lt1 + (dy * inds_h0(0,0))
       lon_min950 = ln1 + (dy * inds_h0(0,1))
       lat_max950 = lt1 + (dy * inds_v0(0,0))
       lon_max950 = ln1 + (dy * inds_v0(0,1))
       print("New centre (geo): "+lat_min950+" degrees N, "+lon_min950+" degrees E")
       print("New centre (vort): "+lat_max950+" degrees N, "+lon_max950+" degrees E")
       centre(ct,0) = lat_max950 ; lat_min950
       centre(ct,1) = lon_max950 ; lon_min950

       delete([/geo_plane0, geo_min0, dims_h0, geo_1d0, inds_h0/])
       delete([/vort_plane0, vort_max0, dims_v0, vort_1d0, inds_v0/])

      else
       print("Cyclone centre (geo): "+lat_min950+" degrees N, "+lon_min950+" degrees E")
       print("Cyclone centre (vort): "+lat_max950+" degrees N, "+lon_max950+" degrees E")
       centre(ct,0) = lat_max950 ; lat_min950
       centre(ct,1) = lon_max950 ; lon_min950
      end if

     end if

   ; Get information about box size for later ensemble calculations (use VORT not GEO)
     r = rad   ; Box radius (degrees)
     box = ddx({lat_max950-r:lat_max950+r},{lon_max950-r:lon_max950+r})
     box_size = dimsizes(box)
     size0    = box_size(0)
     size1    = box_size(1)

   ; Tidy up before looping over pressure levels
     delete([/vort1d, inds, vort1d_sm, vort_smth, inds_sm/])
     delete([/geo1d, inds_h, geo1d_sm, geo_smth, inds_h_sm/])
     delete([/vort_plane, geo_plane/])

;============================================================
; Loop over pressure levels to calculate vertical profile
;============================================================

    cp = 0               ; Counter variable

    do ilev = 0, plevs-1 ; 150 - 250 (25 hPa), 250 - 750 (50 hPa), 750 - 1000 (25 hPa)

     p = pres(ilev)
     print_clock("Working on: " + p + " hPa")

   ; Storm-relative velocity
     u_plane       = u_rel(ilev,{lat_0:lat_1},{lon_0:lon_1})
     v_plane       = v_rel(ilev,{lat_0:lat_1},{lon_0:lon_1})

   ; Other variables
     vort_plane    = vort_sub(ilev,{lat_0:lat_1},{lon_0:lon_1})
     avo_plane	   = avo_sub(ilev,{lat_0:lat_1},{lon_0:lon_1})
     div_plane	   = div_sub(ilev,{lat_0:lat_1},{lon_0:lon_1})
     geo_plane     = geo_sub(ilev,{lat_0:lat_1},{lon_0:lon_1})
     t_plane       = tk_sub(ilev,{lat_0:lat_1},{lon_0:lon_1})

     vort_max      = max(vort_plane)
     vort_smth     = smth9_Wrap(vort_plane, 0.5, 0.5, True)
     vort_max_smth = max(vort_smth)

     zmin          = toint(min(geo_plane)) ; Minimum geopotential on pressure level
     zmax          = toint(max(geo_plane)) ; Maximum geopotential on pressure level

     dims          = dimsizes(vort_plane)
     vort1d        = ndtooned(vort_plane)
     inds          = ind_resolve(maxind(vort1d),dims)
     vort1d_sm     = ndtooned(vort_smth)
     inds_sm       = ind_resolve(maxind(vort1d_sm),dims)

     lat_max = lat_0 + (dy * inds(0,0))
     lon_max = lon_0 + (dx * inds(0,1))

     lat_plane = ind(lat_0.le.lat.and.lat.le.lat_1) ; Indices of latitude array
     lon_plane = ind(lon_0.le.lon.and.lon.le.lon_1) ; Indices of longitude array

   ; Modify latitude/longitude arrays
     lat_size = dimsizes(lat_plane)		  	  	; Array sizes
     lon_size = dimsizes(lon_plane)

     lat_sub = lat(lat_plane(0):lat_plane(lat_size-1)) ; Create subset of points 
     lon_sub = lon(lon_plane(0):lon_plane(lon_size-1)) ; between values given above
     delete([/vort1d, inds, vort1d_sm, vort_smth, inds_sm/])

   ; Account for cases in which another maximum is chosen in error 
     lat_diff = abs(lat_max-lat_max950)
     lon_diff = abs(lon_max-lon_max950)

     if (lon_diff.gt.2.0.or.lat_diff.gt.2.0) then
        lat_max = lat_max950
        lon_max = lon_max950
     end if
     print("Cyclone centre: "+lat_max+" degrees N, "+lon_max+" degrees E")

   ; Add lat/lon values to array for later calculations
     pos_arr(ct,0) = lat_max
     pos_arr(ct,1) = lon_max

;=====================================================================================
; Create array of box positions either side of the original [Davis & Galarneau 2009]
;=====================================================================================

     lat_box = (/lat_max-(2*dy),lat_max-(dy),lat_max,\
     	         lat_max+(dy),lat_max+(2*dy)/)
     lon_box = (/lon_max-(2*dx),lon_max-(dx),lon_max,\
                 lon_max+(dx),lon_max+(2*dx)/)

     lsize   = dimsizes(lat_box) * dimsizes(lon_box)
   
;====================================================
; Define arrays to hold ensemble of box positions
;====================================================

     ddx_box    = new((/lsize,size0,size1/),float)
     ddy_box    = new((/lsize,size0,size1/),float)
     vort_box   = new((/lsize,size0,size1/),float)
     avo_box    = new((/lsize,size0,size1/),float)
     div_box    = new((/lsize,size0,size1/),float)
     u_box      = new((/lsize,size0,size1/),float)
     v_box      = new((/lsize,size0,size1/),float)

     lat_d      = new(lsize,float)
     lon_d      = new(lsize,float)
     lat_u      = new(lsize,float)
     lon_u      = new(lsize,float)

     d_ave      = new(lsize,float)
     dz0        = new(lsize,float)
     dz1        = new(lsize,float)
     div_ave    = new(lsize,float)
     vort_ave   = new(lsize,float)
     vort_avg   = new(lsize,float)
     avo_avg    = new(lsize,float)
     u_avg      = new(lsize,float)
     v_avg      = new(lsize,float)

     vort_avg1  = new(lsize,float)
     avo_avg1   = new(lsize,float)
     u_avg1     = new(lsize,float)
     v_avg1     = new(lsize,float)
 
;=============================================
; Calculate terms in the circulation budget
;=============================================

   ; Radius of box for calculation (degrees)
     r = rad
     r1 = sprintf("%0.1fdeg",2*r) ; Create string to add to output file

;===================================================================
; Loop over perturbed box positions, either side of the TC centre
;===================================================================

   do i = 0, dimsizes(lon_box)-1    ; Longitude
     do j = 0, dimsizes(lat_box)-1  ; Latitude

      a	  = (dimsizes(lon_box) * i) + j

      ddx_box(a,:,:)    = ddx({lat_box(j)-r:lat_box(j)+r},{lon_box(i)-r:lon_box(i)+r})
      ddy_box(a,:,:)    = ddy({lat_box(j)-r:lat_box(j)+r},{lon_box(i)-r:lon_box(i)+r})
      vort_box(a,:,:)   = vort_plane({lat_box(j)-r:lat_box(j)+r},\
				     {lon_box(i)-r:lon_box(i)+r})
      avo_box(a,:,:)    = avo_plane({lat_box(j)-r:lat_box(j)+r},\
      		      	            {lon_box(i)-r:lon_box(i)+r})
      div_box(a,:,:)    = div_plane({lat_box(j)-r:lat_box(j)+r},\
				    {lon_box(i)-r:lon_box(i)+r})
      u_box(a,:,:)      = u_plane({lat_box(j)-r:lat_box(j)+r},\
      		      	          {lon_box(i)-r:lon_box(i)+r})
      v_box(a,:,:)      = v_plane({lat_box(j)-r:lat_box(j)+r},\
			 	  {lon_box(i)-r:lon_box(i)+r})

    ; Change units 
      ddx_box(a,:,:) = ddx_box(a,:,:) / 1000
      ddy_box(a,:,:) = ddy_box(a,:,:) / 1000
      u_box(a,:,:)   = u_box(a,:,:)   / 1000
      v_box(a,:,:)   = v_box(a,:,:)   / 1000

    ; Calculate distance between four corner points of the 2 x 2 degree box
      lat_d(a) = lat_box(j)-r
      lon_d(a) = lon_box(i)-r
      lat_u(a) = lat_box(j)+r
      lon_u(a) = lon_box(i)+r

    ;=====================================================
    ; Calculate distances around each perimeter (m)
    ;=====================================================

    ; Meridional boundaries (same both sides)
      dm0 = r_earth * ((pii/180) * (2*r))
      dm0@units = "m"

    ; 'Bottom' boundary
      dz0(a) = r_earth * ((pii/180) * (2*r)) * cos((pii/180)*lat_d(a))
      dz0@units = "m"

    ; 'Top' boundary
      dz1(a) = r_earth * ((pii/180) * (2*r)) * cos((pii/180)*lat_u(a))
      dz1@units = "m"

    ; Use approximation to calculate the area of the box
      d_ave(a) = ((dz0(a) + dz1(a)) / 2) * dm0
      d_ave@units = "m~S~-2~N~"

    ; Calculate divergence averaged over the box
      div_ave(a) = wgt_areaave(div_box(a,:,:), 1.0, 1.0, 0) ; Divergence
      div_ave@units = "s~S~-1~N~"

    ; Calculate area-averaged vorticity 
      vort_ave(a) = wgt_areaave(vort_box(a,:,:), 1.0, 1.0, 0) ; Relative vorticity
      vort_ave@units = "s~S~-1~N~"

    ; Number of grid points along each side of box surrounding cyclone
      numPOINTS_x0 = dimsizes(ddx_box(a,:,:))
      numPOINTS_x  = numPOINTS_x0(1)   ; Longitude points
      numPOINTS_y  = numPOINTS_x0(0)   ; Latitude points

    ;==========================================================
    ; New arrays to hold integral sums for each side of box
    ;==========================================================

    ; Now calculate terms along the four edges of the box, for each ensemble member
    ; Hence four arrays required for each term in the budget

    ; All arrays have dimensions [lsize * 13] -> ALL members on a SINGLE pressure level

    ; Relative vorticity (to compute perimeter average)
      vort_x_final0 = new((/lsize,numPOINTS_x/),float)
      vort_x_final1 = new((/lsize,numPOINTS_x/),float)
      vort_y_final0 = new((/lsize,numPOINTS_y/),float)
      vort_y_final1 = new((/lsize,numPOINTS_y/),float)

    ; Absolute vorticity (to compute perimeter average)
      avo_x_final0 = new((/lsize,numPOINTS_x/),float)
      avo_x_final1 = new((/lsize,numPOINTS_x/),float)
      avo_y_final0 = new((/lsize,numPOINTS_y/),float)
      avo_y_final1 = new((/lsize,numPOINTS_y/),float) 

    ; Zonal velocity (to compute perimeter average)
      ux_final0 = new((/lsize,numPOINTS_x/),float)
      ux_final1 = new((/lsize,numPOINTS_x/),float)
      uy_final0 = new((/lsize,numPOINTS_y/),float)
      uy_final1 = new((/lsize,numPOINTS_y/),float)

    ; Meridional velocity (to compute perimeter average)
      vx_final0 = new((/lsize,numPOINTS_x/),float)
      vx_final1 = new((/lsize,numPOINTS_x/),float)
      vy_final0 = new((/lsize,numPOINTS_y/),float)
      vy_final1 = new((/lsize,numPOINTS_y/),float)

    ;==========================================================================
    ; Integrate over each perimeter to calculate average relative vorticity
    ;==========================================================================

    ; Bottom
      do ilon = 0, numPOINTS_x-1
        vort_x_final0(a,ilon) = vort_box(a,0,ilon)
        avo_x_final0(a,ilon)  = avo_box(a,0,ilon)
        ux_final0(a,ilon)     = u_box(a,0,ilon)
        vx_final0(a,ilon)     = v_box(a,0,ilon)
      end do 

      rvx0 = avg(vort_x_final0(a,:))
      avx0 = avg(avo_x_final0(a,:))
      ux0  = avg(ux_final0(a,:))
      vx0  = avg(vx_final0(a,:))

    ; Right
      do ilat = 0, numPOINTS_y-1
        vort_y_final0(a,ilat) = vort_box(a,ilat,numPOINTS_x-1)
        avo_y_final0(a,ilat)  = avo_box(a,ilat,numPOINTS_x-1)
        uy_final0(a,ilat)     = u_box(a,ilat,numPOINTS_x-1)
        vy_final0(a,ilat)     = v_box(a,ilat,numPOINTS_x-1)
      end do

      rvy0 = avg(vort_y_final0(a,:))
      avy0 = avg(avo_y_final0(a,:))
      uy0  = avg(uy_final0(a,:))
      vy0  = avg(vy_final0(a,:))

    ; Top
      do ilon = 0, numPOINTS_x-1
        vort_x_final1(a,ilon) = vort_box(a,numPOINTS_y-1,ilon)
        avo_x_final1(a,ilon)  = avo_box(a,numPOINTS_y-1,ilon)
        ux_final1(a,ilon)     = u_box(a,numPOINTS_y-1,ilon)
        vx_final1(a,ilon)     = v_box(a,numPOINTS_y-1,ilon)
      end do

      rvx1 = avg(vort_x_final1(a,:))
      avx1 = avg(avo_x_final1(a,:))
      ux1 = avg(ux_final1(a,:))
      vx1 = avg(vx_final1(a,:))

    ; Left
      do ilat = 0, numPOINTS_y-1
        vort_y_final1(a,ilat) = vort_box(a,ilat,0)
        avo_y_final1(a,ilat)  = avo_box(a,ilat,0)
        uy_final1(a,ilat)     = u_box(a,ilat,0)
        vy_final1(a,ilat)     = v_box(a,ilat,0)
      end do

      rvy1 = avg(vort_y_final1(a,:))
      avy1 = avg(avo_y_final1(a,:))
      uy1 = avg(uy_final1(a,:))
      vy1 = avg(vy_final1(a,:))

    ; Averages using method above
      vort_avg(a) = (rvx0 + rvx1 + rvy0 + rvy1) / 4
      avo_avg(a)  = (avx0 + avx1 + avy0 + avy1) / 4
      u_avg(a)    = (ux0 + ux1 + uy0 + uy1) / 4
      v_avg(a)    = (vx0 + vx1 + vy0 + vy1) / 4

      if (u_avg(a).gt.100)
        print("Should be quitting because of box-average velocity error!!!")
        print("u [bottom] = "+ux0)
        print("u [top] = "+ux1)
        print("u [right] = "+uy0)
        print("u [left] = "+uy1)
        print("u_avg(a) = "+u_avg(a))
        exit()
      end if

;====================================================================================
; Second round of calculations (Horizontal advection, relative vorticity, tilting)
;====================================================================================

   ; Calculate terms on northern/southern boundaries
     circ_x0     = dim_sum_n_Wrap((u_box(a,:,:) * ddx_box(a,:,:)),1)
     circ_s      = circ_x0(0)              ; Circulation (south)
     circ_n      = -circ_x0(numPOINTS_y-1) ; Circulation (north)

   ; Calculate terms on eastern/western boundaries
     circ_y0     = dim_sum_n_Wrap((v_box(a,:,:) * ddy_box(a,:,:)),0)
     circ_w      = -circ_y0(0)             ; Circulation (west)
     circ_e      = circ_y0(numPOINTS_x-1)  ; Circulation (east)

;========================
; Further calculations
;========================

   ; Calculate circulation from line integrals (4 sides of box)
     circ_instT = circ_s + circ_n + circ_w + circ_e

   ; Enter values into predefined arrays
     circ_inst(a,ct,cp)     = circ_instT   	             ; Circulation (LHS)

   ; Tidy up (delete variables from within lat/lon loops before next iteration)
     delete([/vort_x_final0, vort_x_final1, vort_y_final0, vort_y_final1/])
     delete([/avo_x_final0,  avo_x_final1,  avo_y_final0,  avo_y_final1/])
     delete([/uy_final0,     uy_final1,     ux_final0,     ux_final1/])
     delete([/vy_final0,     vy_final1,     vx_final0,     vx_final1/])
     delete([/circ_instT, numPOINTS_x0, numPOINTS_x, numPOINTS_y/])

     end do     ; End latitude loop (do j = 0, dimsizes(lat_box)-1)
    end do	; End longitude loop (do i = 0, dimsizes(lon_box)-1)

    cp                   = cp + 1       ; Counter variable (pressure)

  ; Tidy up (delete variables from within pressure level loop)
    delete([/u_plane, v_plane/])
    delete([/lat_plane, lon_plane, vort_plane, avo_plane, div_plane/])
    delete([/ddx_box, ddy_box, u_box, v_box, vort_box, avo_box, div_box/])
    delete([/div_ave,vort_ave,vort_avg,avo_avg,u_avg,v_avg/])

   end do ; End pressure level loop (do ilev = 0, plevs-1)

   ct                  = ct + 1  ; Counter variable (time)

  end do                ; End time loop (do it = 0, times-1)

; Tidy up (delete variables from within input file loop before next iteration)
  delete([/lon,lat,u,v,w,t,z,vort,div,f0,avo,d/])
  delete([/time,times,utc_date,year,month,day,hour,minute,\
           second,date_str,time_str,out_str/])

 end do     ; End input file loop (do nf = 0, numINPUT-1)

;=====================================================
; Calculate term averages over all ensemble members
;=====================================================

   print_clock("Finished ensemble calculations! Now averaging over all members...")
   circ_instA =	dim_avg_n_Wrap(circ_inst,0)

;=========================================
; Integrate over Nepartak's life cycle
;=========================================

 ; Each calculation produces an array of 19 values [pressure levels]
 ; We do this calculation [numTIMES] times
 ; Create array with [numTIMES] * [plevs] elements

   arr      = new((/numTIMES,plevs/),float)
   csize    = dimsizes(arr)
   csize0   = csize(0)

   time_int = new(csize0,"integer")
   time_int = ispan(0,csize0-1,1)
   time_int!0 = "time"

 ; Add metadata before time-height plotting [shouldn't need lines below]
   circ_instA!0    = "time"
   circ_instA&time = time_int
   circ_instA!1    = "lev"
   circ_instA&lev  = pres

 ; Reorder from [time * pressure] to [pressure * time]
   circ_instF      = circ_instA(lev|:,time|:)

;====================================================
; Open workstation for plotting (time-height plots)
;====================================================

     output = "$sam/nepartak/images/cart_coords/era5_"+r1+"_xz_circ"
     wks = gsn_open_wks(opt,output)

;========================
; Options for plotting 
;========================

     gsn_define_colormap(wks,"vort_panel")

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

;==========================================
; Set levels depending on radius of box
;==========================================

     if (r.ge.0.70.and.r.le.0.80) then
     opts_circ@cnLevels                    = (/-1.0, 0.0, 1.0, 2.0, \
                                                4.0, 6.0, 8.0, 10.0, 12.0, 14.0, \
                                                16.0, 18.0, 20.0, 22.0, 24.0, 26.0, \
                                                28.0, 30.0/)
     opts_circ@cnFillColors                = (/6,8,0,12,13,14,15,\
                                               16,18,19,20,21,22,\
                                               23,24,25,26,27,28/)
     else if (r.ge.0.40.and.r.le.0.50) then
     opts_circ@cnLevels                    = (/-1.0, 0.0, 1.0, 2.0, \
                                                3.0, 4.0, 5.0, 6.0, 8.0, 10.0, \
                                                12.0, 14.0, 16.0, 18.0, 20.0, 22.0, \
                                                24.0, 26.0/)
     opts_circ@cnFillColors                = (/6,8,0,12,13,14,15,\
                                               16,18,19,20,21,22,\
                                               23,24,25,26,27,28/)
     else if (r.ge.0.25.and.r.le.0.35) then
     opts_circ@cnLevels                    = (/-0.5, 0.0, 0.5, 1.0, \
                                                1.5, 2.0, 3.0, 4.0, 5.0, 6.0, \
                                                7.0, 8.0, 9.0, 10.0, 11.0, 12.0, \
                                                14.0, 16.0/)
     opts_circ@cnFillColors                = (/6,8,0,12,13,14,15,\
                                               16,18,19,20,21,22,\
                                               23,24,25,26,27,28/)
     end if
     end if
     end if

     opts_circ@gsnPaperOrientation         = "landscape"
     opts_circ@tiMainFontHeightF           = 0.0125
     opts_circ@gsnLeftString               = ""
     opts_circ@gsnRightString              = ""

     opts_circ@pmLabelBarWidthF            = 0.6
     opts_circ@pmLabelBarHeightF           = 0.15
     opts_circ@lbLabelFontHeightF          = 0.01
     opts_circ@lbPerimOn                   = False
     opts_circ@lbLabelFont                 = "Helvetica"
     opts_circ@lbTitleOn                   = True
     opts_circ@lbTitleString = "Normalised circulation (s~S~-1~N~)"
     opts_circ@lbTitleFontHeightF          = 0.015
     opts_circ@lbTitleDirection            = "Across"

     opts_circ@gsnMaximize                 = True
     opts_circ@tmXTOn                      = "False" ; Turn off top x-axis TM
     opts_circ@tmYROn                      = "False" ; Turn off right y-axis TM
     opts_circ@tiYAxisString 		   = "Pressure (hPa)"
     opts_circ@tiXAxisString 		   = ""	     ; No label beneath x-axis
     opts_circ@vpWidthF                    = 0.6     ; Plot width
     opts_circ@vpHeightF                   = 0.5     ; Plot height
     opts_circ@tmXBMode                    = "Explicit"
     opts_circ@tmXBValues                  = (/0,24,48,72,96,120,144/)
     opts_circ@tmXBLabels                  = (/"00Z 3 Jul","00Z 4 Jul",\
                                               "00Z 5 Jul","00Z 6 Jul",\
                                               "00Z 7 Jul","00Z 8 Jul",\
					       "00Z 9 Jul"/)
     opts_circ@tmXBLabelFontHeightF        = 0.01
     opts_circ@tmXBMinorValues             = ispan(0,numTIMES,6)      ; Minor tickmarks
     opts_circ@trYMaxF                     = 1000    ; Include boundary layer
     opts_circ@trYMinF                     = 150     ; Include upper troposphere

     opts_circ@trXMaxF			   = 144     ; Finish at 00Z 9 July
     opts_circ@trXMinF			   = 24	     ; Start at 00Z 4 July 

     opts_circ@tiMainString                = ""
     opts_circ@tiXAxisFontHeightF          = 0.0125
     opts_circ@tiYAxisFontHeightF          = 0.0125
     opts_circ@gsnDraw                     = False   ; Do not draw the plot
     opts_circ@gsnFrame                    = False   ; Do not advance the frame
     opts_circ@lbLabelBarOn                = True    ; Keep label bar on

;=================
; Plot the data
;=================

     plot_circ = gsn_csm_pres_hgt(wks,circ_instF,opts_circ)
     draw(plot_circ)
     frame(wks)

end