; Script to compute the flux form of the vorticity equation 
; and to calculate changes in circulation due to its different terms  

; Run using [ncl 'opt="x11"' rad=0.25 dist=2.0 nbox=25 tc_era5_xz_flux.ncl]

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
   numTIMES    = 157
   centre      = new((/numTIMES,2/),float) ; Array to hold all [120] lat/lon values

 ; Read in values from text files 
   centre(:,0) = asciiread(diri+"tc_era5_lat.txt",(/numTIMES/),"float")
   centre(:,1) = asciiread(diri+"tc_era5_lon.txt",(/numTIMES/),"float")

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

; Lat/lon arrays from ERA5 reanalysis data
  lat_arr(:,0) = centre(:,0)-dist
  lat_arr(:,1) = centre(:,0)+dist
  lon_arr(:,0) = centre(:,1)-dist
  lon_arr(:,1) = centre(:,1)+dist

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
    numTIMES    = 157              	    ; Total number of times in file [157]
    plevs	= 25			    ; Pressure levels

    time_arr    = new(numTIMES,string)	    ; Array to hold strings for second loop
    title_arr   = new(numTIMES,string)	    ; Array to hold strings for plot titles
    pos_arr     = new((/numTIMES,2/),float) ; Array for lat/lon values of cyclone centre
    r_earth     = 6.37e6                    ; Earth's radius (constant)

    ct		= 0			    ; Counter variable (time)

    date_str 	= new(numTIMES,string)
    time_str 	= new(numTIMES,string)
    out_str  	= new(numTIMES,string)

  ; Create arrays to hold circulation tendency values
  ; [box ensemble (nbox) * time (numTIMES) * levs (plevs)]
    circ_inst 	= new((/nbox,numTIMES,plevs/),float)
    hadv_ten  	= new((/nbox,numTIMES,plevs/),float)
    str_ten   	= new((/nbox,numTIMES,plevs/),float)
    flux_ten  	= new((/nbox,numTIMES,plevs/),float)
    tilt_ten  	= new((/nbox,numTIMES,plevs/),float)
    tot_ten   	= new((/nbox,numTIMES,plevs/),float)

;===============================
; Start loop over input file 
;===============================

  do nf = 0, 25			; Leave out last file (one time only) 
  f = addfile(fili(nf),"r")

;==================================
; Get the variables we will need
;==================================

   do it = 0, 5     ; Loop over ALL times in each file

 ; Original grid --> 157 (time) x 25 (vertical) x 121 (latitude) x 161 (longitude)

   print_clock("Working on time: "+ct)

   t0    = lat_arr(ct,0)
   t1    = lat_arr(ct,1)
   n0    = lon_arr(ct,0)
   n1    = lon_arr(ct,1)
   print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

   lon   = f->longitude({n0:n1})	       ; Longitude
   lat   = f->latitude({t0:t1})                ; Latitude

   time  = f->time(it)                         ; Single time
   pres  = f->level(24:0)                      ; Pressure (hPa) - reverse pressure levels
   z1    = f->z(it,24:0,{t0:t1},{n0:n1})       ; Geopotential (m**2 s**-2)
   t2    = f->t(it,24:0,{t0:t1},{n0:n1})       ; Temperature (K)
   u1    = f->u(it,24:0,{t0:t1},{n0:n1})       ; Zonal velocity (m/s)
   v1    = f->v(it,24:0,{t0:t1},{n0:n1})       ; Meridional velocity (m/s)
   omg1  = f->w(it,24:0,{t0:t1},{n0:n1})       ; Vertical velocity (Pa/s)
   vo1   = f->vo(it,24:0,{t0:t1},{n0:n1})      ; Relative vorticity (/s)
   div1	 = f->d(it,24:0,{t0:t1},{n0:n1})       ; Divergence (/s) 
   dy    = lat(1) - lat(0)		       ; Grid spacing (lat)
   dx    = lon(1) - lon(0)   	               ; Grid spacing (lon)

; Pressure levels (pres)
; (0) 1000, (1) 975, (2) 950, (3) 925, (4) 900, (5) 875, (6) 850, (7) 825
; (8) 800, (9) 775, (10) 750, (11) 700, (12) 650, (13) 600
; (14) 550, (15) 500, (16) 450, (17) 400, (18) 350
; (19) 300, (20) 250, (21) 225, (22) 200, (23) 175, (24) 150

;==========================================================
; Unpack the data using scale factors/offset and tidy up
;==========================================================

   z    = tofloat((z1 * z1@scale_factor) + z1@add_offset)
   copy_VarCoords(z1, z)
   copy_VarAtts(z1, z)
   copy_VarMeta(z1, z)
   delete([/z@missing_value_original, z@_FillValue_original, z1/])

   t    = tofloat((t2 * t2@scale_factor) + t2@add_offset)
   copy_VarCoords(t2, t)
   copy_VarAtts(t2, t)
   copy_VarMeta(t2, t)
   delete([/t@missing_value_original, t@_FillValue_original, t2/])

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

   omg  = tofloat((omg1 * omg1@scale_factor) + omg1@add_offset)
   copy_VarCoords(omg1, omg)
   copy_VarAtts(omg1, omg)
   copy_VarMeta(omg1, omg)
   delete([/omg@missing_value_original, omg@_FillValue_original, omg1/])

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

   f1    = coriolis_param(lat)                           ; Coriolis parameter
   f0    = conform_dims(dimsizes(vort),f1,1)             ; Create same-sized array
   avo   = vort + f0                                     ; Absolute vorticity
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

;===========================
; Additional calculations
;===========================

   gr = 9.80665  ; Gravitational acceleration
   gr@units = "m s~S~-2~N~"

   z  = z / (10 * gr)    ; Calculate height from geopotential and convert to decametres
   z@units  = "dam"

   pii = 3.14159625	 ; Degrees to radians

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

   date_str(ct) = sprinti("%0.2iUTC ", hour) + \
                  sprinti("%0.2i ", day) \
                  + month_abbr(month)
   time_str(ct) = sprinti("%0.2i UTC ", hour) + \
                  sprinti("%0.2i ", day) \
                  + month_abbr(month)
   out_str(ct)  = sprinti("%0.2i", day) + \
                  month_abbr(month) + \
                  "_" + sprinti("%0.2iZ", hour)

   time_arr(ct) = out_str(ct)

   print_clock("Working on time: "+time_str(ct))
   title_arr(ct) = "Valid at: "+time_str(ct)

;==========================================================================
; Define new arrays for terms in vorticity equation (no friction for now)
;==========================================================================

  ; Storm-relative velocity
    u_rel      = u - u_cyc(ct)
    v_rel      = v - v_cyc(ct)

    d          = dimsizes(vort)       ; Size of all arrays [19 * 180 * 180]
    d0         = d(0)-1
    d1         = d(1)-1
    d2         = d(2)-1
    lat_rad    = lat * (pii/180)      ; Convert to radians for calculation below
    lon_rad    = lon * (pii/180)
 
  ; Troubleshooting (look for abnormally large values and exit the script if found)
    do p1 = 0, d0       ; pressure
      do l1 = 0, d1     ; latitude
        do l2 = 0, d2   ; longitude
          if (u_rel(p1,l1,l2) .gt. 200)
           print("Quitting because of storm-relative velocity error!!!")
            exit()
          end if
        end do
      end do
    end do

    tilt_x   = new(d,float) ; Tilting (x-component)
    tilt_y   = new(d,float) ; Tilting (y-component)

  ; Add basic metadata (modify details later)
    copy_VarMeta(vort,tilt_x)
    copy_VarAtts(vort,tilt_x)
    copy_VarCoords(vort,tilt_x)

    copy_VarMeta(vort,tilt_y)
    copy_VarAtts(vort,tilt_y)
    copy_VarCoords(vort,tilt_y)

    copy_VarMeta(u,u_rel)
    copy_VarAtts(u,u_rel)
    copy_VarCoords(u,u_rel)

    copy_VarMeta(v,v_rel)
    copy_VarAtts(v,v_rel)
    copy_VarCoords(v,v_rel)

;==========================================================
; Define latitude and longitude for finite differencing
;==========================================================

    y2d     = new((/d(1),d(2)/),float)
    y2d!0   = "lat"
    y2d!1   = "lon"
    y2d&lat =  lat
    y2d&lon =  lon

    y2d     = conform_dims(dimsizes(y2d),lat,0) ; Array of latitude values
    coslat  = cos(y2d*pii/180)                  ; Same array in radians

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

;===============================================================
; Tilting term (w * dv_dp, -w * du_dp)
;===============================================================

  ; Calculate vertical derivatives of horizontal [storm-relative] velocity
    du_dp = u_rel
    dv_dp = v_rel

  ; Convert pressure values to Pa for calculation below
    pres = pres * 100
    pres@units = "Pa"

  ; Fill arrays with vertical derivatives at each pressure level
  ; (0) 1000, (1) 950, ..., (17) 150, (18) 100
    do plevsN = 0, plevs-1

       kp1 = min((/plevsN+1,plevs-1/)) ; Account for top and bottom levels
       km1 = max((/plevsN-1,0/))

       du_dp(plevsN,:,:) = (u_rel(kp1,:,:) - \
                            u_rel(km1,:,:))\
                            /((pres(km1))-(pres(kp1)))
       dv_dp(plevsN,:,:) = (v_rel(kp1,:,:) - \
                            v_rel(km1,:,:))\
                            /((pres(km1))-(pres(kp1)))

    end do

  ; Convert pressure values back to hPa
    pres = pres / 100
    pres@units = "hPa"

  ; Calculate [dv_dp * w]
    tilt_y = dv_dp * w    ; X-component of relative vorticity

  ; Calculate [-du_dp * w]
    tilt_x = -du_dp * w   ; Y-component of relative vorticity

  ; Metadata
    tilt_x@units = "m s~S~-2~N~"
    tilt_x@standard_name = "tilting_x"
    tilt_x@long_name = "Tilting (x-component)"
    tilt_x@title = "Tilting (x-component)"
    tilt_x@name = "tilting_x"

    tilt_y@units = "m s~S~-2~N~"
    tilt_y@standard_name = "tilting_y"
    tilt_y@long_name = "Tilting (y-component)"
    tilt_y@title = "Tilting (y-component)"
    tilt_y@name = "tilting_y"

;================================================================================
; Find cyclone centre on subset of global grid (use both height and vorticity)
;================================================================================

   ; Choose a level in the lower troposphere (950 hPa)
     ilev = 2

   ; First, concentrate on a subset of the global grid 
     geo_plane     = z(ilev,:,:)
     vort_plane    = vort(ilev,:,:)

   ; RELATIVE VORTICITY
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
     geo_min	    = min(geo_plane)
     geo_smth	    = smth9_Wrap(geo_plane, 0.5, 0.5, True)
     geo_min_smth   = min(geo_smth)
      
     dims_h	    = dimsizes(geo_plane)
     geo1d	    = ndtooned(geo_plane)
     inds_h	    = ind_resolve(minind(geo1d),dims_h)
     geo1d_sm	    = ndtooned(geo_smth)
     inds_h_sm	    = ind_resolve(minind(geo1d_sm),dims_h)            

     lat_min950 = t0 + (dy * inds_h(0,0))
     lon_min950 = n0 + (dx * inds_h(0,1))
     print("Cyclone centre (hgt): "+lat_min950+" degrees N, "+lon_min950+" degrees E")

   ;========================================================
   ; EMPLOY A SAFETY NET IF WE IDENTIFY THE WRONG CYCLONE
   ;========================================================
     if (ct.eq.0) then
;      centre(ct,0) = lat_max950
;      centre(ct,1) = lon_max950
      centre(ct,0) = lat_min950
      centre(ct,1) = lon_min950

     else
      lt0   = centre(ct-1,0)            ; Previous TC latitude
      ln0   = centre(ct-1,1)            ; Previous TC longitude
      print("lt0 = "+lt0+" ; ln0 = "+ln0)
;      d_lat = abs(lat_max950-lt0)       ; Change in latitude (t1 - t0)
;      d_lon = abs(lon_max950-ln0)       ; Change in longitude (n1 - n0)
      d_lat = abs(lat_min950-lt0)       ; Change in latitude (t1 - t0)
      d_lon = abs(lon_min950-ln0)       ; Change in longitude (n1 - n0)
      print("d_lat = "+d_lat+" ; d_lon = "+d_lon)

      if (d_lat.gt.1.0.or.d_lon.gt.1.0) then      ; Recalculate TC centre if incorrect

     ; Create smaller grid [0.5 degrees]
       inc = 0.5
       lt1 = lt0 - inc
       lt2 = lt0 + inc
       ln1 = ln0 - inc
       ln2 = ln0 + inc

       geo_plane0 = z(ilev,{lt1:lt2},{ln1:ln2})          ; Smaller grid (1.0 x 1.0 deg)
       geo_min0   = min(geo_plane0)                      ; Find minimum on smaller grid
       dims_h0    = dimsizes(geo_plane0)                 ; Size of smaller grid
       geo_1d0    = ndtooned(geo_plane0)                 ; Create 1-D array
       inds_h0    = ind_resolve(minind(geo_1d0),dims_h0) ; Find index of minimum

       vort_plane0 = vort(ilev,{lt1:lt2},{ln1:ln2})        ; Smaller grid (1.0 x 1.0 deg)
       vort_max0   = max(vort_plane0)                      ; Find minimum on smaller grid
       dims_v0     = dimsizes(vort_plane0)                 ; Size of smaller grid
       vort_1d0    = ndtooned(vort_plane0)                 ; Create 1-D array
       inds_v0     = ind_resolve(minind(vort_1d0),dims_v0) ; Find index of minimum

       lat_max950  = lt1 + (dy * inds_v0(0,0))
       lon_max950  = ln1 + (dy * inds_v0(0,1))
       print("New centre (vort): "+lat_max950+" degrees N, "+lon_max950+" degrees E")

       lat_min950 = lt1 + (dy * inds_h0(0,0))
       lon_min950 = ln1 + (dy * inds_h0(0,1))
       print("New centre (geo): "+lat_min950+" degrees N, "+lon_min950+" degrees E")

;       centre(ct,0) = lat_max950
;       centre(ct,1) = lon_max950
       centre(ct,0) = lat_min950
       centre(ct,1) = lon_min950

       delete([/geo_plane0, geo_min0, dims_h0, geo_1d0, inds_h0/])
       delete([/vort_plane0, vort_max0, dims_v0, vort_1d0, inds_v0/])

      else
       print("Cyclone centre (geo): "+lat_min950+" degrees N, "+lon_min950+" degrees E")
       print("Cyclone centre (vort): "+lat_max950+" degrees N, "+lon_max950+" degrees E")
;       centre(ct,0) = lat_max950
;       centre(ct,1) = lon_max950
       centre(ct,0) = lat_min950
       centre(ct,1) = lon_min950
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

    do ilev = 0, plevs-1 ; 1000 - 150 hPa 

     p = pres(ilev)
     print_clock("Working on: " + p + " hPa")

   ; Storm-relative velocity
     u_plane       = u_rel(ilev,:,:)
     v_plane       = v_rel(ilev,:,:)

   ; Other variables
     vort_plane    = vort(ilev,:,:)
     avo_plane	   = avo(ilev,:,:)
     div_plane	   = div(ilev,:,:)
     geo_plane     = z(ilev,:,:)
     t_plane       = t(ilev,:,:)

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

     lat_max 	   = t0 + (dy * inds(0,0))
     lon_max 	   = n0 + (dx * inds(0,1))

     delete([/vort1d, inds, vort1d_sm, vort_smth, inds_sm/])

   ; Account for cases in which another maximum is chosen in error 
     lat_diff = abs(lat_max-lat_max950)
     lon_diff = abs(lon_max-lon_max950)

     if (lon_diff.gt.1.0.or.lat_diff.gt.1.0) then
;        lat_max = lat_max950
;        lon_max = lon_max950
        lat_max = lat_min950
        lon_max = lon_min950
     end if
     print("Cyclone centre: "+lat_max+" degrees N, "+lon_max+" degrees E")

   ; Add lat/lon values to array for later calculations
     pos_arr(ct,0) = lat_max
     pos_arr(ct,1) = lon_max

   ; Terms in the circulation budget
   ; 'tilt_x' and 'tilt_y' contain values on all pressure levels for each time
   ; 'tilt_x_plane' and 'tilt_y_plane' contain values on single pressure levels
     tilt_x_plane               = tilt_x(ilev,:,:)       ; Tilting (x-component)
     tilt_y_plane               = tilt_y(ilev,:,:)       ; Tilting (y-component)

;=====================================================================================
; Create array of box positions either side of the original [Davis & Galarneau 2009]
;=====================================================================================

     if (nbox.eq.49) then
      lat_box = (/lat_max-(3*dy),lat_max-(2*dy),lat_max-(dy),lat_max,\
                  lat_max+(dy),lat_max+(2*dy),lat_max+(3*dy)/)
      lon_box = (/lon_max-(3*dx),lon_max-(2*dx),lon_max-(dx),lon_max,\
                  lon_max+(dx),lon_max+(2*dx),lon_max+(3*dx)/)
     else if (nbox.eq.25) then
      lat_box = (/lat_max-(2*dy),lat_max-(dy),lat_max,\
                  lat_max+(dy),lat_max+(2*dy)/)
      lon_box = (/lon_max-(2*dx),lon_max-(dx),lon_max,\
                  lon_max+(dx),lon_max+(2*dx)/)
     else if (nbox.eq.9) then
      lat_box = (/lat_max-(dy),lat_max,lat_max+(dy)/)
      lon_box = (/lon_max-(dx),lon_max,lon_max+(dx)/)
     else if (nbox.eq.1) then
      lat_box = (/lat_max/)
      lon_box = (/lon_max/)
     end if
     end if
     end if
     end if

     lsize   = dimsizes(lat_box) * dimsizes(lon_box)
   
;====================================================
; Define arrays to hold ensemble of box positions
;====================================================

     ddx_box    = new((/lsize,size0,size1/),float)
     ddy_box    = new((/lsize,size0,size1/),float)
     tilt_x_box = new((/lsize,size0,size1/),float)
     tilt_y_box = new((/lsize,size0,size1/),float)
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
      print_clock("Working on box number "+a)

      ddx_box(a,:,:)    = ddx({lat_box(j)-r:lat_box(j)+r},{lon_box(i)-r:lon_box(i)+r})
      ddy_box(a,:,:)    = ddy({lat_box(j)-r:lat_box(j)+r},{lon_box(i)-r:lon_box(i)+r})
      tilt_x_box(a,:,:) = tilt_x_plane({lat_box(j)-r:lat_box(j)+r},\
                                       {lon_box(i)-r:lon_box(i)+r})
      tilt_y_box(a,:,:) = tilt_y_plane({lat_box(j)-r:lat_box(j)+r},\
                                       {lon_box(i)-r:lon_box(i)+r})
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
     tilt_x0     = dim_sum_n_Wrap((tilt_x_box(a,:,:) * ddx_box(a,:,:)),1)
     flux_x0     = dim_sum_n_Wrap((-avo_box(a,:,:) * v_box(a,:,:) * \
                                   ddx_box(a,:,:)),1)
     hadv_x0     = dim_sum_n_Wrap((-(avo_box(a,:,:) - avo_avg(a)) * \
                                  (v_box(a,:,:) - v_avg(a)) * \
                                   ddx_box(a,:,:)),1)
     circ_x0     = dim_sum_n_Wrap((u_box(a,:,:) * ddx_box(a,:,:)),1)
     tilt_s      = tilt_x0(0)              ; Tilting (south)
     tilt_n      = tilt_x0(numPOINTS_y-1)  ; Tilting (north)
     flux_s      = flux_x0(0)              ; Advective flux divergence (south)
     flux_n      = -flux_x0(numPOINTS_y-1) ; Advective flux divergence (north)
     hadv_s      = hadv_x0(0)              ; Horizontal advection (south)
     hadv_n      = -hadv_x0(numPOINTS_y-1) ; Horizontal advection (north)
     circ_s      = circ_x0(0)              ; Circulation (south)
     circ_n      = -circ_x0(numPOINTS_y-1) ; Circulation (north)

   ; Calculate terms on eastern/western boundaries
     tilt_y0    = dim_sum_n_Wrap((tilt_y_box(a,:,:) * ddy_box(a,:,:)),0)
     flux_y0    = dim_sum_n_Wrap((-avo_box(a,:,:) * u_box(a,:,:) * \
                                  ddy_box(a,:,:)),0)
     hadv_y0    = dim_sum_n_Wrap((-(avo_box(a,:,:) - avo_avg(a)) * \
                                 (u_box(a,:,:) - u_avg(a)) * \
                                  ddy_box(a,:,:)),0)
     circ_y0    = dim_sum_n_Wrap((v_box(a,:,:) * ddy_box(a,:,:)),0)
     tilt_w     = tilt_y0(0)              ; Tilting (west)
     tilt_e     = tilt_y0(numPOINTS_x-1)  ; Tilting (east)
     flux_w     = flux_y0(0)              ; Advective flux divergence (west)
     flux_e     = -flux_y0(numPOINTS_x-1) ; Advective flux divergence (east)
     hadv_w     = hadv_y0(0)              ; Horizontal advection (west)
     hadv_e     = -hadv_y0(numPOINTS_x-1) ; Horizontal advection (east)
     circ_w     = -circ_y0(0)             ; Circulation (west)
     circ_e     = circ_y0(numPOINTS_x-1)  ; Circulation (east)

;========================
; Further calculations
;========================

   ; Calculate stretching term using averaged quantities
     str_tenT   = -(div_ave(a) * avo_avg(a))

   ; Calculate horizontal advection from line integrals (4 sides of box)
     hadv_tenT  = hadv_s + hadv_n + hadv_w + hadv_e
     hadv_tenT  = hadv_tenT / d_ave(a) ; Divide by area of box [Rios-Berrios et al. 2016]

;   ; Calculate advective flux divergence from line integrals (4 sides of box)
;     flux_tenT  = flux_s + flux_n + flux_w + flux_e
;     flux_tenT  = flux_tenT / d_ave(a) ; Divide by area of box[Rios-Berrios et al. 2016]

   ; Calculate tilting tendency from line integrals (4 sides of box)
     tilt_tenT  = tilt_s + tilt_n + tilt_w + tilt_e
     tilt_tenT  = tilt_tenT / d_ave(a) ; Divide by area of box [Rios-Berrios et al. 2016]

   ; Calculate circulation from line integrals (4 sides of box)
     circ_instT = circ_s + circ_n + circ_w + circ_e

   ; Enter values into predefined arrays
   ; Enter values into predefined arrays
     hadv_ten(a,ct,cp)      = hadv_tenT                   ; Horizontal advection
     str_ten(a,ct,cp)       = str_tenT                    ; Stretching
;     flux_ten(a,ct,cp)      = flux_tenT                  ; Advective flux
     tilt_ten(a,ct,cp)      = tilt_tenT                   ; Tilting
     circ_inst(a,ct,cp)     = circ_instT   	          ; Circulation (LHS)

   ; Tidy up (delete variables from within lat/lon loops before next iteration)
;     delete([/flux_tenT/])
     delete([/hadv_tenT, str_tenT, tilt_tenT, circ_instT/])
     delete([/vort_x_final0, vort_x_final1, vort_y_final0, vort_y_final1/])
     delete([/avo_x_final0,  avo_x_final1,  avo_y_final0,  avo_y_final1/])
     delete([/uy_final0,     uy_final1,     ux_final0,     ux_final1/])
     delete([/vy_final0,     vy_final1,     vx_final0,     vx_final1/])
     delete([/numPOINTS_x0, numPOINTS_x, numPOINTS_y/])
     delete([/tilt_y0,tilt_x0,circ_y0,circ_x0,flux_x0,flux_y0,hadv_x0,hadv_y0/])

     end do     ; End latitude loop (do j = 0, dimsizes(lat_box)-1)
    end do	; End longitude loop (do i = 0, dimsizes(lon_box)-1)

    cp                   = cp + 1       ; Counter variable (pressure)

  ; Tidy up (delete variables from within pressure level loop)
    delete([/u_plane, v_plane, vort_plane, avo_plane, div_plane/])
    delete([/tilt_x_plane, tilt_y_plane/])
    delete([/ddx_box, ddy_box, u_box, v_box, vort_box, avo_box, div_box/])
    delete([/div_ave,vort_ave,vort_avg,avo_avg,u_avg,v_avg/])

   end do ; End pressure level loop (do ilev = 0, plevs-1)

 ; Tidy up (delete variables from within input file loop before next iteration)
   delete([/lon,lat,u,v,w,t,z,vort,div,f0,avo,d,box/])
   delete([/time,times,utc_date,year,month,day,hour,minute,second/])

   ct = ct + 1	        ; Counter variable (time)

  end do                ; End time loop (do it = 0, 5)

 end do			; End input file loop (do nf = 0, 25)

;=====================================================
; Calculate term averages over all ensemble members
;=====================================================

   print_clock("Finished ensemble calculations! Now averaging over all members...")

   hadv_final = new((/numTIMES,plevs/),"float")
   str_final  = new((/numTIMES,plevs/),"float")
;   flux_final = new((/numTIMES,plevs/),"float")
   tilt_final = new((/numTIMES,plevs/),"float")
   circ_final = new((/numTIMES,plevs/),"float")

   circ_inst  = circ_inst / 1000000 ; Divide by 1000000 [following 'tc_ens_xz_circ']

  ; Average over 0th dimension (ensemble of box positions)
    hadv_final(:,:) = dim_avg_n_Wrap(hadv_ten(:,:,:),0)
    str_final(:,:)  = dim_avg_n_Wrap(str_ten(:,:,:),0)
;    flux_final(:,:) = dim_avg_n_Wrap(flux_ten(:,:,:),0)
    tilt_final(:,:) = dim_avg_n_Wrap(tilt_ten(:,:,:),0)
    circ_final(:,:) = dim_avg_n_Wrap(circ_inst(:,:,:),0)

;=========================================
; Integrate over Nepartak's life cycle
;=========================================

 ; Each calculation produces an array of [plevs] values (pressure levels)
 ; We do this calculation [numTIMES] times
 ; Create array with [numTIMES] * [plevs] elements

   hadv_intA = new((/numTIMES,plevs/),float)
   str_intA  = new((/numTIMES,plevs/),float)
   flux_intA = new((/numTIMES,plevs/),float)
   tilt_intA = new((/numTIMES,plevs/),float)
   tot_intA  = new((/numTIMES,plevs/),float)

   hadv_intA@units = "10~S~-4~N~ s~S~-1~N~"
   str_intA@units  = "10~S~-4~N~ s~S~-1~N~"
   flux_intA@units = "10~S~-4~N~ s~S~-1~N~"
   tilt_intA@units = "10~S~-4~N~ s~S~-1~N~"
   tot_intA@units  = "10~S~-4~N~ s~S~-1~N~"

   csize    = dimsizes(flux_intA)
   csize0   = csize(0)

   time_int = new(csize0,"integer")
   time_int = ispan(0,csize0-1,1)
   time_int!0 = "time"

   print_clock("Starting plots...")

;   do ic = 24, 144, 1  ; 1 h intervals between 00Z 4 and 9 July

   do ic = 0, csize0-2, 1 ; 1 h intervals (skip final time)

     if (ic.lt.csize0-1) then
      cs = ic                 ; Start index
      cf = cs + 1             ; Finish index
     else
      cs = ic
      cf = ic
     end if

      f_str = sprinti("%0.3i",cs)+"h"    ; Formatted string for output file

      hadv_intA(ic,:) = dim_sum_n_Wrap((hadv_final(cs:cf,:) * 3600), 0)
      str_intA(ic,:)  = dim_sum_n_Wrap((str_final(cs:cf,:) * 3600), 0)
;      flux_intA(ic,:) = dim_sum_n_Wrap((flux_final(cs:cf,:) * 3600), 0)
      tilt_intA(ic,:) = dim_sum_n_Wrap((tilt_final(cs:cf,:) * 3600), 0)

      if (fct.eq.4) then
       hadv_intA(ic,:) = hadv_intA(ic,:) * (10 ^ 4)
       str_intA(ic,:)  = str_intA(ic,:) * (10 ^ 4)
       tilt_intA(ic,:) = tilt_intA(ic,:) * (10 ^ 4)
      else if (fct.eq.5) then
       hadv_intA(ic,:) = hadv_intA(ic,:) * (10 ^ 5)
       str_intA(ic,:)  = str_intA(ic,:) * (10 ^ 5)
       tilt_intA(ic,:) = tilt_intA(ic,:) * (10 ^ 5)
      end if
      end if

    end do               ; End time loop (do ic = 0, csize0-1)

    flux_intA = hadv_intA + str_intA     ; Advective flux divergence
    tot_intA = flux_intA + tilt_intA     ; Calculate total circulation tendency (RHS)

    printMinMax(flux_intA,False)
    printMinMax(tilt_intA,False)
    printMinMax(tot_intA,False)

  ; Add metadata before time-height plotting
    tot_intA!0     = "time"
    tot_intA&time  = time_int
    tot_intA!1     = "lev"
    tot_intA&lev   = pres

    hadv_intA!0    = "time"
    hadv_intA&time = time_int
    hadv_intA!1    = "lev"
    hadv_intA&lev  = pres

    str_intA!0    = "time"
    str_intA&time = time_int
    str_intA!1    = "lev"
    str_intA&lev  = pres

    flux_intA!0    = "time"
    flux_intA&time = time_int
    flux_intA!1    = "lev"
    flux_intA&lev  = pres

    tilt_intA!0    = "time"
    tilt_intA&time = time_int
    tilt_intA!1    = "lev"
    tilt_intA&lev  = pres

;======================================
; Reorder dimensions before plotting
;======================================

 ; Define new arrays for final plotting
   circ_int = new((/plevs,numTIMES/),float)
   hadv_int = new((/plevs,numTIMES/),float)
   str_int  = new((/plevs,numTIMES/),float)
   flux_int = new((/plevs,numTIMES/),float)
   tilt_int = new((/plevs,numTIMES/),float)
   tot_int  = new((/plevs,numTIMES/),float)

 ; Add metadata to circulation array
   circ_final!0     = "time"
   circ_final&time  = time_int
   circ_final!1     = "lev"
   circ_final&lev   = pres

 ; Reorder from [time * pressure] to [pressure * time]
   circ_int(:,:) = circ_final(lev|:,time|:)
   hadv_int(:,:) = hadv_intA(lev|:,time|:)
   str_int(:,:)  = str_intA(lev|:,time|:)
   flux_int(:,:) = flux_intA(lev|:,time|:)
   tilt_int(:,:) = tilt_intA(lev|:,time|:)
   tot_int(:,:)  = tot_intA(lev|:,time|:)

;====================================================
; Open workstation for plotting (time-height plots)
;====================================================

     output = "$sam/nepartak/images/cart_coords/era5_"+r1+"_xz_flux"
     wks = gsn_open_wks(opt,output)

;========================
; Options for plotting 
;========================

     panel = new(1,graphic)
     gsn_define_colormap(wks,"BlueDarkRed18")

   ; Circulation budget terms
     opts_bdgt                             = True
     opts_bdgt@cnFillOn                    = True
     opts_bdgt@cnLineLabelInterval         = 2.0
     opts_bdgt@cnLineLabelFontHeightF      = 0.012
     opts_bdgt@cnLineLabelBackgroundColor  = "transparent"
     opts_bdgt@cnLineLabelPlacementMode    = "constant"
     opts_bdgt@cnLinesOn                   = False
     opts_bdgt@cnInfoLabelOn               = False
     opts_bdgt@cnLevelSelectionMode        = "ExplicitLevels"

     if (fct.eq.4) then
      opts_bdgt@cnLevels                    = (/-15.0, -10.0, -5.0, -2.0, \
                                                -1.5, -1.0, -0.5, -0.2, 0.2, \
                                                 0.5, 1.0, 1.5, 2.0, 5.0, 10.0, 15.0/)
     else if (fct.eq.5) then
      opts_bdgt@cnLevels                    = (/-30.0, -25.0, -20.0, -15.0, -10.0, \
                                                -5.0, -2.0, -1.0, 1.0, 2.0, \
                                                 5.0, 10.0, 15.0, 20.0, 25.0, 30.0/)
     end if
     end if

     opts_bdgt@cnFillColors                = (/3,4,5,6,7,8,9,10,0,11,12,\
                                               13,14,15,16,18,19/)

     opts_bdgt@gsnPaperOrientation         = "landscape"
     opts_bdgt@tiMainFontHeightF           = 0.0125
     opts_bdgt@gsnLeftString               = ""
     opts_bdgt@gsnRightString              = ""

     opts_bdgt@gsnMaximize                 = True
     opts_bdgt@tmXTOn                      = "False" ; Turn off top x-axis TM
     opts_bdgt@tmYROn                      = "False" ; Turn off right y-axis TM
     opts_bdgt@tiYAxisString 		   = "Pressure (hPa)"
     opts_bdgt@tiYAxisFontHeightF          = 0.030   ; Larger y-axis title
     opts_bdgt@tiXAxisString               = "Time"  ; X-axis title
     opts_bdgt@tiXAxisFontHeightF          = 0.030   ; Larger x-axis title
     opts_bdgt@vpWidthF                    = 0.6     ; Plot width
     opts_bdgt@vpHeightF                   = 0.5     ; Plot height
     opts_bdgt@tmXBMode                    = "Explicit"
;     opts_bdgt@tmXBValues                  = (/24,48,72,96,120,144/)
;     opts_bdgt@tmXBLabels                  = (/"00Z 4 Jul","00Z 5 Jul",\
;                                               "00Z 6 Jul","00Z 7 Jul",\
;                                               "00Z 8 Jul","00Z 9 Jul"/)

     opts_bdgt@tmXBValues                  = (/48,72,96,120/)
     opts_bdgt@tmXBLabels                  = (/"00Z 5 Jul","00Z 6 Jul",\
                                               "00Z 7 Jul","00Z 8 Jul"/)

     opts_bdgt@tmXBLabelFontHeightF        = 0.025
     opts_bdgt@trYMaxF                     = 1000    ; Include boundary layer
     opts_bdgt@trYMinF                     = 150     ; Include upper troposphere
     opts_bdgt@trXMaxF			   = 132     ; Finish at 12Z 8 July
     opts_bdgt@trXMinF			   = 48	     ; Start at 00Z 5 July 

     opts_bdgt@tmXBMinorValues             = ispan(opts_bdgt@trXMinF,\
     					     opts_bdgt@trXMaxF,6)	; Minor tickmarks

     opts_bdgt@tiMainString                = ""
     opts_bdgt@gsnDraw                     = False   ; Do not draw the plot
     opts_bdgt@gsnFrame                    = False   ; Do not advance the frame
     opts_bdgt@lbLabelBarOn                = False   ; Keep label bar on

   ; Additional resources for black contours of circulation (LHS)
     opts_cont                             = True           ; New resources
     opts_cont@cnFillOn                    = False          ; Turn off filled contours
     opts_cont@cnLineColor                 = "black"        ; Line colour
     opts_cont@gsnDraw                     = False          ; Do not draw the plot
     opts_cont@gsnFrame                    = False          ; Do not advance the frame
     opts_cont@cnLevelSelectionMode        = "ExplicitLevels"

     if (r.ge.0.70.and.r.le.0.80) then
      opts_cont@cnLevels                     = (/10.0, 14.0, 18.0, 22.0, 26.0/)
     else if (r.ge.0.40.and.r.le.0.50) then
      opts_cont@cnLevels                     = (/6.0, 8.0, 10.0, 14.0, 18.0/)
     else if (r.ge.0.25.and.r.le.0.35) then
      opts_cont@cnLevels                     = (/2.0, 4.0, 6.0, 8.0, 10.0/)
     end if
     end if
     end if

     opts_cont@cnLineThicknessF             = 5.0           ; Contour line thickness
     opts_cont@cnLineLabelsOn               = False         ; Turn off contour labels
     opts_cont@gsnContourNegLineDashPattern = 2             ; Negative contours dashed
     opts_cont@gsnContourPosLineDashPattern = 0             ; Positive contours solid
     opts_cont@cnInfoLabelOn                = False         ; Turn off info label
     circ_int@units                         = ""            ; Turn off string above plot

;=================
; Plot the data
;=================

     plot_flux06  = gsn_csm_pres_hgt(wks,flux_int(:,:),opts_bdgt)
     plot_cont06a = gsn_csm_contour(wks,circ_int(:,:),opts_cont)

   ; Overlay contours of circulation onto filled contour plots
     overlay(plot_flux06,plot_cont06a)

     panel(0)  = plot_flux06

     optsP                       = True
     optsP@gsnFrame              = False    ; Do not advance the frame
     optsP@gsnPanelLabelBar      = False    ; Turn on panel labelbar
     optsP@pmLabelBarWidthF      = 0.6
     optsP@pmLabelBarHeightF     = 0.15
     optsP@lbLabelFontHeightF    = 0.01
     optsP@lbPerimOn             = False
     optsP@lbLabelFont           = "Helvetica"
     optsP@lbTitleString = "Circulation budget tendency terms (10~S~-4~N~ s~S~-1~N~)"
     optsP@lbTitleFontHeightF    = 0.015
     optsP@lbTitleDirection      = "Across"
     optsP@lbTitlePosition       = "Bottom"
     optsP@txString              = ""
     optsP@gsnPanelFigureStrings = "" ; (/"d) ERA5 reanalysis (31 km)"/)
     optsP@gsnMaximize           = True
     optsP@gsnPanelTop           = 0.90
     optsP@gsnPanelBottom        = 0.10
     optsP@gsnPaperOrientation   = "landscape"
     optsP@amJust                = "TopLeft"
     optsP@gsnPanelFigureStringsFontHeightF = 0.0125 ; Label size (default 0.01)

     gsn_panel(wks,panel,(/1,1/),optsP)            ; Draw as a single plot
     frame(wks)

end