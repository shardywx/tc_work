; Script to illustrate the representation of vorticity aggregation in different datasets

; ncl 'opt="x11"' 'ens0="em07"' 'ens1="em07"' rad=1.5 nd=0.75 \
; dist=4.5 ir=0.5 tc_vort_panel3_paper.ncl

; Where 'opt' = output file format (x11 or PDF)
; 'ens0' and 'ens1' = CP and global MetUM ensemble simulation extensions (00-11)
; 'rad' = preferred radius for the final plots (degrees)
; 'nd' = radius from the TC centre when finding the geopotential height minimum
; 'dist' = extent of gridded data either side of cyclone (subset - speeds up script)
; 'ir' = radius of box (circulation budget) to overlay onto panel plots 

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"

begin

;=================================================
; Read in cyclone position data from text files
;=================================================

   diri        = "/nfs/a37/earshar/ncl_scripts/nepartak/archer_march2018/"
   numTIMES    = 120
   centre      = new((/4,numTIMES,2/),float) ; Array to hold positions for all 4 datasets

 ;===============================
 ; Convection-permitting MetUM 
 ;===============================
 ; Read in values from 'em07'
   centre(0,:,0) = asciiread(diri+"cp_lat_"+ens0+".txt",(/numTIMES/),"float")
   centre(0,:,1) = asciiread(diri+"cp_lon_"+ens0+".txt",(/numTIMES/),"float")

 ;=================
 ; Global MetUM
 ;=================
 ; Read in values from 'em02'
   centre(1,:,0) = asciiread(diri+"cp_lat_"+ens1+".txt",(/numTIMES/),"float")
   centre(1,:,1) = asciiread(diri+"cp_lon_"+ens1+".txt",(/numTIMES/),"float")

 ;===================
 ; ERA5 reanalysis
 ;===================
 ; Read in values directly from ERA5 reanalysis
   centre(2,:,0) = asciiread(diri+"tc_era5_lat.txt",(/numTIMES/),"float")
   centre(2,:,1) = asciiread(diri+"tc_era5_lon.txt",(/numTIMES/),"float")

;============================================================================
; Calculate distance between two points --> then calculate system velocity
;============================================================================

  print_clock("Starting system velocity calculations...")

; Define new arrays (cyclone position, speed)
  x_cyc   = new((/3,numTIMES/),float)
  y_cyc   = new((/3,numTIMES/),float)
  u_cyc   = new((/3,numTIMES/),float)
  v_cyc   = new((/3,numTIMES/),float)
  vel_cyc = new((/3,numTIMES/),float)

 do pt = 0, 2

 ; Loop to calculate system velocity from hourly cyclone position [lat,lon]
   do iang = 1, numTIMES-1
    x_cyc(pt,iang)   = gc_latlon(centre(pt,iang-1,0),centre(pt,iang-1,1),\
                       centre(pt,iang-1,0),centre(pt,iang,1),2,3)
    y_cyc(pt,iang)   = gc_latlon(centre(pt,iang-1,0),centre(pt,iang-1,1),\
                       centre(pt,iang,0),centre(pt,iang-1,1),2,3)
    u_cyc(pt,iang)   = x_cyc(pt,iang)/(3600)
    v_cyc(pt,iang)   = y_cyc(pt,iang)/(3600)
    vel_cyc(pt,iang) = sqrt(u_cyc(pt,iang)^2.0 + v_cyc(pt,iang)^2.0)
   end do

 ; Set initial system velocity [0] equal to that at the second time [1]
   x_cyc(pt,0)   = 0
   y_cyc(pt,0)   = 0
   u_cyc(pt,0)   = u_cyc(pt,1)
   v_cyc(pt,0)   = v_cyc(pt,1)
   vel_cyc(pt,0) = vel_cyc(pt,1)

 ; Add metadata to new arrays
   u_cyc@units         = "m s~S~-1~N~"
   u_cyc@description   = "Cyclone zonal velocity"

   vel_cyc@units       = "m s~S~-1~N~"
   vel_cyc@description = "Cyclone speed"

   u_cyc(pt,:) = -(u_cyc(pt,:)) ; Easterly motion is negative

  end do      ; End of loop over datasets (do pt = 0, 2)

;================================================================================
; Create lat/lon arrays containing grid subsets at each time, following the TC
;================================================================================

  lat_arr = new((/3,numTIMES,2/),float)
  lon_arr = new((/3,numTIMES,2/),float)

; Lat/lon arrays from CP MetUM data
  lat_arr(0,:,0) = centre(0,:,0)-dist
  lat_arr(0,:,1) = centre(0,:,0)+dist
  lon_arr(0,:,0) = centre(0,:,1)-dist
  lon_arr(0,:,1) = centre(0,:,1)+dist

; Lat/lon arrays from global MetUM data
  lat_arr(1,:,0) = centre(1,:,0)-dist
  lat_arr(1,:,1) = centre(1,:,0)+dist
  lon_arr(1,:,0) = centre(1,:,1)-dist
  lon_arr(1,:,1) = centre(1,:,1)+dist

; Lat/lon arrays from ERA5 reanalysis data
  lat_arr(2,:,0) = centre(2,:,0)-dist
  lat_arr(2,:,1) = centre(2,:,0)+dist
  lon_arr(2,:,0) = centre(2,:,1)-dist
  lon_arr(2,:,1) = centre(2,:,1)+dist

;==================================================================================
;;;;;;;;;;;;;;;;;;;;;;;;;;; Start of main part of script ;;;;;;;;;;;;;;;;;;;;;;;;;;
;==================================================================================

  print_clock("Starting main part of script...")

    time_arr    = new(numTIMES,string)      ; Array to hold strings for second loop
    title_arr   = new(numTIMES,string)      ; Array to hold strings for plot titles
    pos_arr     = new((/numTIMES,2/),float) ; Array for lat/lon values of cyclone centre
    r_earth     = 6.37e6                    ; Earth's radius (constant)

    date_str = new(numTIMES,string)
    time_str = new(numTIMES,string)
    out_str  = new(numTIMES,string)

    setvalues NhlGetWorkspaceObjectId
      "wsMaximumSize" : 1000000000
    end setvalues

    cp_input  = "$ens/reg_netcdf/20160704T0000Z_ra1t_em07_pb.nc"
    gl_input  = "$ens/global_netcdf/20160704T0000Z_engl_em02_pb.nc"
    era_input = "$sam/era5/cdf/era5_nepartak_all.nc"
    
    f	     = addfile(cp_input,"r")
    g	     = addfile(gl_input,"r")
    i	     = addfile(era_input,"r")

;================================================================
; Start loop over input files: 00Z 4 July to 8 July, every 1 h
;================================================================

   do it = 0, 96

   setvalues NhlGetWorkspaceObjectId
     "wsMaximumSize" : 300000000
   end setvalues

   ct = it + 25		; Counter variable (ERA5 reanalysis)

;=======================================================
; Get the variables we will need, from all 4 datasets
;=======================================================

;=============
; CP MetUM 
;=============

  time_cp  = f->t(it)
  times_cp = dimsizes(time_cp)

  t0c = lat_arr(0,it,0)
  t1c = lat_arr(0,it,1)
  n0c = lon_arr(0,it,0)
  n1c = lon_arr(0,it,1)
  print_clock("Grid extends from "+t0c+" to "+t1c+" deg N, and "+n0c+" to "+n1c+" deg E")

; Pressure levels (pres)
; (0) 1000, (1) 950, (2) 900, (3) 850, (4) 800, (5) 750, (6) 700, (7) 650
; (8) 600, (9) 550, (10) 500, (11) 450, (12) 400, (13) 350
; (14) 300, (15) 250, (16) 200, (17) 150, (18) 100

  p_cp     = (/1,3,6,10,16/)          ; Pressure levels to access (950,850,700,500,200)

  lon1_cp  = f->longitude_1({n0c:n1c})  ; longitude
  lat1_cp  = f->latitude_1({t0c:t1c})   ; latitude
  lon_cp   = f->longitude({n0c:n1c})    ; longitude
  lat_cp   = f->latitude({t0c:t1c})     ; latitude

  pres_cp  = f->p(p_cp)        	        ; Pressure (levels)
  plevs_cp = dimsizes(pres_cp) 	        ; Size of pressure level array	
  pres_cp@units = "hPa"

; Arrays have 5 vertical levels, and horizontal extent depends on 'dist' option
  u_cp     = f->u(it,p_cp,{t0c:t1c},{n0c:n1c})      ; Zonal velocity
  v_cp     = f->v(it,p_cp,{t0c:t1c},{n0c:n1c})      ; Meridional velocity
  vort_cp  = f->rvor(it,p_cp,{t0c:t1c},{n0c:n1c})   ; Relative vorticity

  t_cp     = f->temp(it,p_cp,{t0c:t1c},{n0c:n1c})   ; Temperature
  z_cp     = f->ht(it,p_cp,{t0c:t1c},{n0c:n1c})     ; Geopotential height

; Convert to units we want
  z_cp     = z_cp/10
  z_cp@units = "dam" ; Convert to decametres
  vort_cp  = vort_cp * (10 ^ -6)   ; Convert to /s

  f_cp0	   = coriolis_param(lat_cp)			   ; Coriolis parameter 
  f_cp     = conform_dims(dimsizes(vort_cp),f_cp0,1)	   ; Create same-sized array
  avo_cp   = vort_cp + f_cp				   ; Absolute vorticity
  copy_VarCoords(vort_cp,avo_cp)
  copy_VarAtts(vort_cp,avo_cp)
  copy_VarMeta(vort_cp,avo_cp)
  avo_cp   = avo_cp * (10 ^ 4)				   ; Convert units for plotting

  vort_cp  = vort_cp * (10 ^ 4)    ; Convert units for plotting
  dx_cp    = lon_cp(1) - lon_cp(0) ; Grid spacing (longitude)
  dy_cp    = lat_cp(1) - lat_cp(0) ; Grid spacing (latitude)

;==============
; Global UM 
;==============

; 6 [time] x 19 [vertical] x 600 [latitude] x 800 [longitude]

  time_gl  = g->t(it)
  times_gl = dimsizes(time_gl)

  t0g = lat_arr(1,it,0)
  t1g = lat_arr(1,it,1)
  n0g = lon_arr(1,it,0)
  n1g = lon_arr(1,it,1)
  print_clock("Grid extends from "+t0g+" to "+t1g+" deg N, and "+n0g+" to "+n1g+" deg E")

; Pressure levels (pres)
; (0) 1000, (1) 950, (2) 900, (3) 850, (4) 800, (5) 750, (6) 700, (7) 650
; (8) 600, (9) 550, (10) 500, (11) 450, (12) 400, (13) 350
; (14) 300, (15) 250, (16) 200, (17) 150, (18) 100

  p_gl     = (/1,3,6,10,16/)          ; Pressure levels to access (950,850,700,500,200)

  lon1_gl  = g->longitude_1({n0g:n1g})  ; longitude
  lat1_gl  = g->latitude_1({t0g:t1g})   ; latitude
  lon_gl   = g->longitude({n0g:n1g})    ; longitude
  lat_gl   = g->latitude({t0g:t1g})     ; latitude

  pres_gl  = g->p(p_gl)               ; Pressure (levels)
  plevs_gl = dimsizes(pres_gl)        ; Size of pressure level array
  pres_gl@units = "hPa"

; Arrays have 5 vertical levels, and horizontal extent depends on 'dist' option
  u_gl     = g->u(it,p_gl,{t0g:t1g},{n0g:n1g})      ; Zonal velocity
  v_gl     = g->v(it,p_gl,{t0g:t1g},{n0g:n1g})      ; Meridional velocity
  vort_gl  = g->rvor(it,p_gl,{t0g:t1g},{n0g:n1g})   ; Relative vorticity

  t_gl     = g->temp(it,p_gl,{t0g:t1g},{n0g:n1g})   ; Temperature
  z_gl     = g->ht(it,p_gl,{t0g:t1g},{n0g:n1g})     ; Geopotential height

; Convert to units we want
  z_gl     = z_gl/10
  z_gl@units = "dam" ; Convert to decametres
  vort_gl  = vort_gl * (10 ^ -6) ; Convert to /s

  f_gl0    = coriolis_param(lat_gl)                        ; Coriolis parameter
  f_gl	   = conform_dims(dimsizes(vort_gl),f_gl0,1)       ; Create same-sized array
  avo_gl   = vort_gl + f_gl                                ; Absolute vorticity
  copy_VarCoords(vort_gl,avo_gl)
  copy_VarAtts(vort_gl,avo_gl)
  copy_VarMeta(vort_gl,avo_gl)
  avo_gl   = avo_gl * (10 ^ 4)				   ; Convert units for plotting

  vort_gl  = vort_gl * (10 ^ 4)	 ; Convert to preferred units for plotting
  dx_gl    = lon_gl(1) - lon_gl(0) ; Grid spacing (longitude)
  dy_gl    = lat_gl(1) - lat_gl(0) ; Grid spacing (latitude)

;===================
; ERA5 reanalysis
;===================

; Original grid --> 6 (time) x 25 (vertical) x 121 (latitude) x 161 (longitude)

  t0a = lat_arr(2,it,0)
  t1a = lat_arr(2,it,1)
  n0a = lon_arr(2,it,0)
  n1a = lon_arr(2,it,1)
  print_clock("Grid extends from "+t0a+" to "+t1a+" deg N, and "+n0a+" to "+n1a+" deg E")

  lon_era   = i->longitude({n0a:n1a})      ; longitude
  lat_era   = i->latitude({t0a:t1a})       ; latitude

  p_era     = (/22,18,13,9,2/)           ; [950, 850, 700, 500, 200 hPa] for arrays below

  time_era  = i->time(ct)                         ; Single time
  pres_era  = i->level(p_era)                     ; Pressure (hPa)
  z1_era    = i->z(ct,p_era,{t0a:t1a},{n0a:n1a})  ; Geopotential (m**2 s**-2)
  t1_era    = i->t(ct,p_era,{t0a:t1a},{n0a:n1a})  ; Temperature (K)
  u1_era    = i->u(ct,p_era,{t0a:t1a},{n0a:n1a})  ; Zonal velocity (m/s)
  v1_era    = i->v(ct,p_era,{t0a:t1a},{n0a:n1a})  ; Meridional velocity (m/s)
  vo1_era   = i->vo(ct,p_era,{t0a:t1a},{n0a:n1a}) ; Relative vorticity (/s)

  dy_era    = lat_era(1) - lat_era(0)   ; Grid spacing (lat)
  dx_era    = lon_era(1) - lon_era(0)   ; Grid spacing (lon)

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

   time_arr(it) = out_str(it)

   print_clock("Working on time: "+time_str(it))
   title_arr(it) = "Valid at: "+time_str(it)
 
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

   f_era0    = coriolis_param(lat_era)                       ; Coriolis parameter
   f_era     = conform_dims(dimsizes(vort_era),f_era0,1)     ; Create same-sized array
   avo_era   = vort_era + f_era                              ; Absolute vorticity
   copy_VarCoords(vort_era,avo_era)
   copy_VarAtts(vort_era,avo_era)
   copy_VarMeta(vort_era,avo_era)
   avo_era   = avo_era * (10 ^ 4)                            ; Convert units for plotting
   vort_era = vort_era * (10 ^ 4)  			     ; Convert units for plotting

;===========================
; Additional calculations
;===========================

   gr = 9.80665  ; Gravitational acceleration
   gr@units = "m s~S~-2~N~"

   z_era  = z_era / (10 * gr) ; Calculate height from geopotential, convert to decametres
   z_era@units  = "dam"

   pii = 3.14159265 ; Degrees to radians

;==========================================================================
; Define new arrays for terms in vorticity equation (no friction for now)
;==========================================================================

 ; Find cyclone centre using relative vorticity or geopotential height (950 hPa)
   ilev = 0

   vort_pl_cp   = vort_cp(ilev,:,:)			    ; CP MetUM
   vort_pl_gl   = vort_gl(ilev,:,:)			    ; Global MetUM
   vort_pl_era  = vort_era(ilev,:,:)			    ; ERA5 reanalysis

   vort_max_cp  = max(vort_pl_cp)  			    ; CP MetUM
   vort_max_gl  = max(vort_pl_gl)  			    ; Global UM
   vort_max_era = max(vort_pl_era) 			    ; ERA5 reanalysis

  ;==============================================================
  ; Now calculate cyclone centre for each of the four datasets
  ;==============================================================

  ;============================
  ; Convection-permitting UM 
  ;============================
    dims_cp       = dimsizes(vort_pl_cp)
    vort1d_cp     = ndtooned(vort_pl_cp)
    inds_cp       = ind_resolve(maxind(vort1d_cp),dims_cp)
    lat_max950_cp = t0c + (dy_cp * inds_cp(0,0)) ; Max latitude
    lon_max950_cp = n0c + (dx_cp * inds_cp(0,1)) ; Max longitude

    print("Pre-script centre = "+centre(0,it,0)+" deg N, "+centre(0,it,1)+" deg E")
    print("In-script centre = "+lat_max950_cp+" deg N, "+lon_max950_cp+" deg E")

  ; Now search for minimum geopotential height surrounding vorticity maximum
    lt1 = lat_max950_cp - nd
    lt2 = lat_max950_cp + nd
    ln1 = lon_max950_cp - nd
    ln2 = lon_max950_cp + nd

    geo_pl_cp0 = z_cp(ilev,{lt1:lt2},{ln1:ln2})
    geo_min0   = min(geo_pl_cp0)
    dims_h_cp0 = dimsizes(geo_pl_cp0)			    ; Size of smaller grid
    geo_1d_cp0 = ndtooned(geo_pl_cp0)                 	    ; Create 1-D array
    inds_h_cp0 = ind_resolve(minind(geo_1d_cp0),dims_h_cp0) ; Find index of minimum    

    lat_min950_cp = lt1 + (dy_cp * inds_h_cp0(0,0))
    lon_min950_cp = ln1 + (dx_cp * inds_h_cp0(0,1))

    print("Cyclone centre (C-P UM): "+lat_max950_cp+\
    	  " degrees N, "+lon_max950_cp+" degrees E")

    print("Cyclone centre (C-P UM - final): "+lat_min950_cp+\
    	  " degrees N, "+lon_min950_cp+" degrees E")	  

   ;========================================================
   ; EMPLOY A SAFETY NET IF WE IDENTIFY THE WRONG CYCLONE
   ;========================================================
    if (it.eq.0) then
     centre(0,it,0) = lat_min950_cp
     centre(0,it,1) = lon_min950_cp
    else
     lt0	    = centre(0,it-1,0)
     ln0	    = centre(0,it-1,1)
     print("lt0 = "+lt0+" ; ln0 = "+ln0)
     d_lat = abs(lat_min950_cp-lt0)       ; Change in latitude (t1 - t0)
     d_lon = abs(lon_min950_cp-ln0)       ; Change in longitude (note extra minus sign)
     print("d_lat = "+d_lat+" ; d_lon = "+d_lon)     

     if (d_lat.gt.1.or.d_lon.gt.1) then        ; Recalculate TC centre if incorrect

     ; Create smaller grid [0.5 degrees]
       inc = 0.5
       lt3 = lt0 - inc
       lt4 = lt0 + inc
       ln3 = ln0 - inc
       ln4 = ln0 + inc

       gp0        = z_cp(ilev,{lt3:lt4},{ln3:ln4})       ; Smaller grid (1.0 x 1.0 deg)
       gm0        = min(gp0)                             ; Find minimum on smaller grid
       dh0        = dimsizes(gp0)                 	 ; Size of smaller grid
       gd0        = ndtooned(gp0)                 	 ; Create 1-D array
       ih0        = ind_resolve(minind(gd0),dh0)         ; Find index of minimum

       lat_min950_cp  = lt3 + (dy_cp * ih0(0,0))
       lon_min950_cp  = ln3 + (dy_cp * ih0(0,1))
       print("FINAL centre: "+lat_min950_cp+" degrees N, "+lon_min950_cp+" degrees E")
       centre(0,it,0) = lat_min950_cp
       centre(0,it,1) = lon_min950_cp

       delete([/gp0, gm0, dh0, gd0, ih0/])

     else
       print("Cyclone centre (C-P UM - final): "+lat_min950_cp+\
          " degrees N, "+lon_min950_cp+" degrees E")
       centre(0,it,0) = lat_min950_cp
       centre(0,it,1) = lon_min950_cp
     end if
    end if

  ;=============
  ; Global UM
  ;=============
    dims_gl       = dimsizes(vort_pl_gl)
    vort1d_gl     = ndtooned(vort_pl_gl)
    inds_gl       = ind_resolve(maxind(vort1d_gl),dims_gl)
    lat_max950_gl = t0g + (dy_gl * inds_gl(0,0)) ; Max latitude
    lon_max950_gl = n0g + (dx_gl * inds_gl(0,1)) ; Max longitude

    print("Pre-script centre = "+centre(1,it,0)+" deg N, "+centre(1,it,1)+" deg E")
    print("In-script centre = "+lat_max950_gl+" deg N, "+lon_max950_gl+" deg E")

  ; Now search for minimum geopotential height surrounding vorticity maximum
    lt1 = lat_max950_gl - nd
    lt2 = lat_max950_gl + nd
    ln1 = lon_max950_gl - nd
    ln2 = lon_max950_gl + nd

    geo_pl_gl0 = z_gl(ilev,{lt1:lt2},{ln1:ln2})
    geo_min0   = min(geo_pl_gl0)
    dims_h_gl0 = dimsizes(geo_pl_gl0)                       ; Size of smaller grid
    geo_1d_gl0 = ndtooned(geo_pl_gl0)                       ; Create 1-D array
    inds_h_gl0 = ind_resolve(minind(geo_1d_gl0),dims_h_gl0) ; Find index of minimum

    lat_min950_gl = lt1 + (dy_gl * inds_h_gl0(0,0))
    lon_min950_gl = ln1 + (dx_gl * inds_h_gl0(0,1))

    print("Cyclone centre (Global UM): "+lat_max950_gl+\
          " degrees N, "+lon_max950_gl+" degrees E")

    print("Cyclone centre (Global UM - final): "+lat_min950_gl+\
          " degrees N, "+lon_min950_gl+" degrees E")

   ;========================================================
   ; EMPLOY A SAFETY NET IF WE IDENTIFY THE WRONG CYCLONE
   ;========================================================
    if (it.eq.0) then
     centre(1,it,0) = lat_min950_gl
     centre(1,it,1) = lon_min950_gl
    else
     lt0            = centre(1,it-1,0)
     ln0            = centre(1,it-1,1)
     print("lt0 = "+lt0+" ; ln0 = "+ln0)
     d_lat = abs(lat_min950_gl-lt0)       ; Change in latitude (t1 - t0)
     d_lon = abs(lon_min950_gl-ln0)       ; Change in longitude (note extra minus sign)
     print("d_lat = "+d_lat+" ; d_lon = "+d_lon)

     if (d_lat.gt.1.or.d_lon.gt.1) then        ; Recalculate TC centre if incorrect

     ; Create smaller grid [0.5 degrees]
       inc = 0.5
       lt3 = lt0 - inc
       lt4 = lt0 + inc
       ln3 = ln0 - inc
       ln4 = ln0 + inc

       gp0        = z_gl(ilev,{lt3:lt4},{ln3:ln4})       ; Smaller grid (1.0 x 1.0 deg)
       gm0        = min(gp0)                             ; Find minimum on smaller grid
       dh0        = dimsizes(gp0)                        ; Size of smaller grid
       gd0        = ndtooned(gp0)                        ; Create 1-D array
       ih0        = ind_resolve(minind(gd0),dh0)         ; Find index of minimum

       lat_min950_gl  = lt3 + (dy_gl * ih0(0,0))
       lon_min950_gl  = ln3 + (dy_gl * ih0(0,1))
       print("FINAL centre: "+lat_min950_gl+" degrees N, "+lon_min950_gl+" degrees E")
       centre(1,it,0) = lat_min950_gl
       centre(1,it,1) = lon_min950_gl

       delete([/gp0, gm0, dh0, gd0, ih0/])

     else
       print("Cyclone centre (Global UM - final): "+lat_min950_gl+\
          " degrees N, "+lon_min950_gl+" degrees E")
       centre(1,it,0) = lat_min950_gl
       centre(1,it,1) = lon_min950_gl
     end if
    end if

  ;===================
  ; ERA5 reanalysis
  ;===================
    dims_era       = dimsizes(vort_pl_era)
    vort1d_era     = ndtooned(vort_pl_era)
    inds_era       = ind_resolve(maxind(vort1d_era),dims_era)
    lat_max950_era = t0a + (dy_era * inds_era(0,0)) ; Max latitude
    lon_max950_era = n0a + (dx_era * inds_era(0,1)) ; Max longitude

  ; Now search for minimum geopotential height surrounding vorticity maximum
    lt1 = lat_max950_era - nd
    lt2 = lat_max950_era + nd
    ln1 = lon_max950_era - nd
    ln2 = lon_max950_era + nd

    geo_pl_era0 = z_era(ilev,{lt1:lt2},{ln1:ln2})
    geo_min0   = min(geo_pl_era0)
    dims_h_era0 = dimsizes(geo_pl_era0)                       ; Size of smaller grid
    geo_1d_era0 = ndtooned(geo_pl_era0)                       ; Create 1-D array
    inds_h_era0 = ind_resolve(minind(geo_1d_era0),dims_h_era0) ; Find index of minimum

    lat_min950_era = lt1 + (dy_era * inds_h_era0(0,0))
    lon_min950_era = ln1 + (dx_era * inds_h_era0(0,1))

    print("Cyclone centre (ERA5): "+lat_max950_era+\
          " degrees N, "+lon_max950_era+" degrees E")

    print("Cyclone centre (ERA5 - final): "+lat_min950_era+\
          " degrees N, "+lon_min950_era+" degrees E")

   ;========================================================
   ; EMPLOY A SAFETY NET IF WE IDENTIFY THE WRONG CYCLONE
   ;========================================================
    if (it.eq.0) then
     centre(2,it,0) = lat_min950_era
     centre(2,it,1) = lon_min950_era
    else
     lt0            = centre(2,it-1,0)
     ln0            = centre(2,it-1,1)
     print("lt0 = "+lt0+" ; ln0 = "+ln0)
     d_lat = abs(lat_min950_era-lt0)       ; Change in latitude (t1 - t0)
     d_lon = abs(lon_min950_era-ln0)       ; Change in longitude (note extra minus sign)
     print("d_lat = "+d_lat+" ; d_lon = "+d_lon)

     if (d_lat.gt.1.or.d_lon.gt.1) then        ; Recalculate TC centre if incorrect

     ; Create smaller grid [0.5 degrees]
       inc = 0.5
       lt3 = lt0 - inc
       lt4 = lt0 + inc
       ln3 = ln0 - inc
       ln4 = ln0 + inc

       gp0        = z_era(ilev,{lt3:lt4},{ln3:ln4})      ; Smaller grid (1.0 x 1.0 deg)
       gm0        = min(gp0)                             ; Find minimum on smaller grid
       dh0        = dimsizes(gp0)                        ; Size of smaller grid
       gd0        = ndtooned(gp0)                        ; Create 1-D array
       ih0        = ind_resolve(minind(gd0),dh0)         ; Find index of minimum

       lat_min950_era  = lt3 + (dy_era * ih0(0,0))
       lon_min950_era  = ln3 + (dy_era * ih0(0,1))
       print("FINAL centre: "+lat_min950_era+" degrees N, "+lon_min950_era+" degrees E")
       centre(2,it,0) = lat_min950_era
       centre(2,it,1) = lon_min950_era

       delete([/gp0, gm0, dh0, gd0, ih0/])

     else
       print("Cyclone centre (ERA5 - final): "+lat_min950_era+\
          " degrees N, "+lon_min950_era+" degrees E")
       centre(2,it,0) = lat_min950_era
       centre(2,it,1) = lon_min950_era
     end if
    end if

  ;===============================================
  ; Tidy up before looping over pressure levels
  ;===============================================
    delete([/vort_pl_cp,vort1d_cp,inds_cp,dims_cp/])
    delete([/vort_pl_gl,vort1d_gl,inds_gl,dims_gl/])
    delete([/vort_pl_era,vort1d_era,inds_era,dims_era/])
    delete([/geo_pl_era0,dims_h_era0,geo_1d_era0,inds_h_era0/])
    delete([/geo_pl_gl0,dims_h_gl0,geo_1d_gl0,inds_h_gl0/])
    delete([/geo_pl_cp0,dims_h_cp0,geo_1d_cp0,inds_h_cp0/])

  ; Loop over levels (p_cp, p_gl, p_ec, p_era)

;  do ilev = 0, dimsizes(p_cp)-1    ; 950, 850, 700, 500, 200 hPa
  do ilev = 1, 1

     print("Working on: " + pres_cp(ilev) + "hPa")

     lon_box = new((/3,5/),float)				 ; Longitude array
     lat_box = new((/3,5/),float)				 ; Latitude array

     vort_pl_cp   = avo_cp(ilev,:,:)				 ; CP MetUM
     vort_pl_gl   = avo_gl(ilev,:,:)				 ; Global MetUM
     vort_pl_era  = avo_era(ilev,:,:)			         ; ERA5

     u_pl_cp   	  = u_cp(ilev,:,:)                               ; CP MetUM
     u_pl_gl   	  = u_gl(ilev,:,:)                            	 ; Global MetUM
     u_pl_era  	  = u_era(ilev,:,:)			         ; ERA5

     v_pl_cp   	  = v_cp(ilev,:,:)                            	 ; CP MetUM
     v_pl_gl   	  = v_gl(ilev,:,:)                            	 ; Global MetUM
     v_pl_era  	  = v_era(ilev,:,:)				 ; ERA5

     vort_max_cp  = max(vort_pl_cp)  ; CP MetUM
     vort_max_gl  = max(vort_pl_gl)  ; Global MetUM
     vort_max_er  = max(vort_pl_era) ; ERA5 reanalysis

     ;============================
     ; Convection-permitting UM 
     ;============================
     dims_cp       = dimsizes(vort_pl_cp)
     vort1d_cp     = ndtooned(vort_pl_cp)
     inds_cp       = ind_resolve(maxind(vort1d_cp),dims_cp)
     lat_max_cp    = t0c + (dy_cp * inds_cp(0,0)) ; Latitude of maximum vorticity
     lon_max_cp    = n0c + (dx_cp * inds_cp(0,1)) ; Longitude of maximum vorticity
     print("Cyclone centre (C-P UM): "+lat_max_cp+\
           " degrees N, "+lon_max_cp+" degrees E")

     delete([/vort1d_cp,inds_cp,dims_cp/]) ; Don't delete 'vort_pl_cp' this time

     ;============================================================================
     ; Now search for minimum geopotential height surrounding vorticity maximum
     ;============================================================================
     lt1 = lat_max_cp - nd
     lt2 = lat_max_cp + nd
     ln1 = lon_max_cp - nd
     ln2 = lon_max_cp + nd

     geo_pl_cp0 = z_cp(ilev,{lt1:lt2},{ln1:ln2})
     geo_min0   = min(geo_pl_cp0)
     dims_h_cp0 = dimsizes(geo_pl_cp0)                       ; Size of smaller grid
     geo_1d_cp0 = ndtooned(geo_pl_cp0)                       ; Create 1-D array
     inds_h_cp0 = ind_resolve(minind(geo_1d_cp0),dims_h_cp0) ; Find index of minimum

     lat_min_cp = lt1 + (dy_cp * inds_h_cp0(0,0))
     lon_min_cp = ln1 + (dx_cp * inds_h_cp0(0,1))

     ;========================================================================
     ; Account for cases in which another cyclone centre is chosen in error
     ;========================================================================
     lat1_cp1 = abs(lat_min_cp-lat_min950_cp)
     lon1_cp1 = abs(lon_min_cp-lon_min950_cp)

     if (lon1_cp1.gt.0.5.or.lat1_cp1.gt.0.5) then
        print("Fixed an error in the tracking method!")
        lat_min_cp = lat_min950_cp
        lon_min_cp = lon_min950_cp
     end if
     print("Cyclone centre (C-P UM): "+lat_min_cp+\
           " degrees N, "+lon_min_cp+" degrees E")

     ;===============
     ; Global UM
     ;===============
     dims_gl       = dimsizes(vort_pl_gl)
     vort1d_gl     = ndtooned(vort_pl_gl)
     inds_gl       = ind_resolve(maxind(vort1d_gl),dims_gl)
     lat_max_gl    = t0g + (dy_gl * inds_gl(0,0)) ; Max latitude
     lon_max_gl    = n0g + (dx_gl * inds_gl(0,1)) ; Max longitude
     print("Cyclone centre (Global UM): "+lat_max_gl+\
           " degrees N, "+lon_max_gl+" degrees E")

     delete([/vort1d_gl,inds_gl,dims_gl/]) ; Don't delete 'vort_pl_gl' this time

     ;============================================================================
     ; Now search for minimum geopotential height surrounding vorticity maximum
     ;============================================================================
     lt1 = lat_max_gl - nd
     lt2 = lat_max_gl + nd
     ln1 = lon_max_gl - nd
     ln2 = lon_max_gl + nd

     geo_pl_gl0 = z_gl(ilev,{lt1:lt2},{ln1:ln2})
     geo_min0   = min(geo_pl_gl0)
     dims_h_gl0 = dimsizes(geo_pl_gl0)                       ; Size of smaller grid
     geo_1d_gl0 = ndtooned(geo_pl_gl0)                       ; Create 1-D array
     inds_h_gl0 = ind_resolve(minind(geo_1d_gl0),dims_h_gl0) ; Find index of minimum

     lat_min_gl = lt1 + (dy_gl * inds_h_gl0(0,0))
     lon_min_gl = ln1 + (dx_gl * inds_h_gl0(0,1))

     ;========================================================================
     ; Account for cases in which another cyclone centre is chosen in error
     ;========================================================================
     lat1_gl1 = abs(lat_min_gl-lat_min950_gl)
     lon1_gl1 = abs(lon_min_gl-lon_min950_gl)

     if (lon1_gl1.gt.0.9.or.lat1_gl1.gt.0.9) then
        print("Fixed an error in the tracking method!")
        lat_min_gl = lat_min950_gl
        lon_min_gl = lon_min950_gl
     end if
     print("Cyclone centre (Global UM): "+lat_min_gl+\
     	   " degrees N, "+lon_min_gl+" degrees E")

     ;=====================
     ; ERA5 reanalysis
     ;=====================
     dims_era       = dimsizes(vort_pl_era)
     vort1d_era     = ndtooned(vort_pl_era)
     inds_era       = ind_resolve(maxind(vort1d_era),dims_era)
     lat_max_era    = t0a + (dy_era * inds_era(0,0)) ; Max latitude
     lon_max_era    = n0a + (dx_era * inds_era(0,1)) ; Max longitude
     print("Cyclone centre (ERA5): "+lat_max_era+\
           " degrees N, "+lon_max_era+" degrees E")

     delete([/vort1d_era,inds_era,dims_era/]) ; Don't delete 'vort_pl_era' this time

     ;============================================================================
     ; Now search for minimum geopotential height surrounding vorticity maximum
     ;============================================================================
     lt1 = lat_max_era - nd
     lt2 = lat_max_era + nd
     ln1 = lon_max_era - nd
     ln2 = lon_max_era + nd

     geo_pl_era0 = z_era(ilev,{lt1:lt2},{ln1:ln2})
     geo_min0   = min(geo_pl_era0)
     dims_h_era0 = dimsizes(geo_pl_era0)                       ; Size of smaller grid
     geo_1d_era0 = ndtooned(geo_pl_era0)                       ; Create 1-D array
     inds_h_era0 = ind_resolve(minind(geo_1d_era0),dims_h_era0) ; Find index of minimum

     lat_min_era = lt1 + (dy_era * inds_h_era0(0,0))
     lon_min_era = ln1 + (dx_era * inds_h_era0(0,1))

     ;========================================================================
     ; Account for cases in which another cyclone centre is chosen in error
     ;========================================================================
     lat1_era1 = abs(lat_min_era-lat_min950_era)
     lon1_era1 = abs(lon_min_era-lon_min950_era)

     if (lon1_era1.gt.0.75.or.lat1_era1.gt.0.75) then
        print("Fixed an error in the tracking method!")
        lat_min_era = lat_min950_era
        lon_min_era = lon_min950_era
     end if
     print("Cyclone centre (ERA5): "+lat_min_era+\
     	   " degrees N, "+lon_min_era+" degrees E")

;=====================================
; Design boxes to draw around cyclone
;=====================================

   ; User defines 'ir' (radius of integration box around cyclone)
     lat_box(0,:) = (/lat_min_cp-ir,lat_min_cp+ir,lat_min_cp+ir,\
                      lat_min_cp-ir,lat_min_cp-ir/)
     lon_box(0,:) = (/lon_min_cp-ir,lon_min_cp-ir,lon_min_cp+ir,\
                      lon_min_cp+ir,lon_min_cp-ir/)

     lat_box(1,:) = (/lat_min_gl-ir,lat_min_gl+ir,lat_min_gl+ir,\
                      lat_min_gl-ir,lat_min_gl-ir/)
     lon_box(1,:) = (/lon_min_gl-ir,lon_min_gl-ir,lon_min_gl+ir,\
                      lon_min_gl+ir,lon_min_gl-ir/)

     lat_box(2,:) = (/lat_min_era-ir,lat_min_era+ir,lat_min_era+ir,\
                      lat_min_era-ir,lat_min_era-ir/)
     lon_box(2,:) = (/lon_min_era-ir,lon_min_era-ir,lon_min_era+ir,\
                      lon_min_era+ir,lon_min_era-ir/)

;=================================
; Output file type and location
;=================================

     diri = "$sam/nepartak/images/cart_coords/"
     output = diri+\
     "vort_panel3_cart_"+ens0+"_"+ens1+"_"+ir+"box_"+pres_cp(ilev)+"_"+time_arr(it)
     wks = gsn_open_wks(opt,output)

;==========================
; Options for plotting 
;==========================

    panel = new(3,graphic)
    gsn_define_colormap(wks,"vort_panel")

  ; Absolute vorticity
    opts_vort_cp       	 	     	    = True
    opts_vort_cp@cnFillOn		    = True
    opts_vort_cp@cnLinesOn                  = False
    opts_vort_cp@cnLineLabelInterval	    = 2.0
    opts_vort_cp@cnLineLabelFontHeightF     = 0.012
    opts_vort_cp@cnLineLabelBackgroundColor = "transparent"
    opts_vort_cp@cnLineLabelPlacementMode   = "constant"
    opts_vort_cp@cnLinesOn                  = False ; Contour lines off
    opts_vort_cp@cnInfoLabelOn              = False
    opts_vort_cp@cnLevelSelectionMode       = "ExplicitLevels"
    opts_vort_cp@cnLevels                   = (/1.0, 2.0, 4.0, 6.0, 8.0, 10.0, \
						12.0, 15.0, 20.0, 25.0, 30.0, \
						35.0, 40.0, 45.0, 50.0/)
    opts_vort_cp@cnFillColors               = (/0,12,13,15,17,18,19,\
                                                20,21,22,23,24,\
						25,26,27,28/)
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
    opts_vort_cp@gsnPaperOrientation        = "landscape"    ; change orientation
    opts_vort_cp@gsnAddCyclic		    = False
    opts_vort_cp@gsnDraw                    = False	     ; Do not draw the plot
    opts_vort_cp@gsnFrame                   = False          ; Do not advance the frame
    opts_vort_cp@mpLimitMode                = "Corners"

    radius				    = rad	     ; Size of plot
    opts_vort_cp@mpLeftCornerLatF           = lat_min_cp-radius
    opts_vort_cp@mpLeftCornerLonF           = lon_min_cp-radius
    opts_vort_cp@mpRightCornerLatF          = lat_min_cp+radius
    opts_vort_cp@mpRightCornerLonF          = lon_min_cp+radius

  ; Global UM
    opts_vort_gl			    = opts_vort_cp
    opts_vort_gl@mpLeftCornerLatF           = lat_min_gl-radius
    opts_vort_gl@mpLeftCornerLonF           = lon_min_gl-radius
    opts_vort_gl@mpRightCornerLatF          = lat_min_gl+radius
    opts_vort_gl@mpRightCornerLonF          = lon_min_gl+radius

  ; ERA5
    opts_vort_era                           = opts_vort_gl
    opts_vort_era@mpLeftCornerLatF          = lat_min_era-radius
    opts_vort_era@mpLeftCornerLonF          = lon_min_era-radius
    opts_vort_era@mpRightCornerLatF         = lat_min_era+radius
    opts_vort_era@mpRightCornerLonF         = lon_min_era+radius

  ; Wind vectors
    opts_vec                            = True
    opts_vec@gsnDraw                    = False         ; Do not draw the plot
    opts_vec@gsnFrame                   = False         ; Do not advance the frame
    opts_vec@vcGlyphStyle               = "LineArrow"   ; draw windbarbs
    opts_vec@vcRefMagnitudeF            = 20            ; magnitude of reference vector
    opts_vec@vcRefLengthF               = 0.04          ; length of reference vector
    opts_vec@vcMinDistanceF             = 0.05          ; distance between wind barbs
    opts_vec@vcLineArrowThicknessF      = 1.0           ; windbarb line thickness
    opts_vec@vcLineArrowColor           = "black"       ; windbarb colour
    opts_vec@gsnMaximize                = True          ; Maximise plot size
    opts_vec@gsnRightString             = ""
    opts_vec@gsnLeftString              = ""
    opts_vec@tiXAxisString              = ""
    opts_vec@vcRefAnnoOrthogonalPosF    = 0.05

;=============
; Plot data 
;=============

  ; Convection-permitting MetUM
    vort0    = gsn_csm_contour_map(wks,vort_pl_cp,opts_vort_cp)
    vec0     = gsn_csm_vector(wks,u_pl_cp,v_pl_cp,opts_vec)
    overlay(vort0,vec0)
    panel(0) = vort0
    
  ; Global MetUM
    vort1    = gsn_csm_contour_map(wks,vort_pl_gl,opts_vort_gl)
    vec1     = gsn_csm_vector(wks,u_pl_gl,v_pl_gl,opts_vec)
    overlay(vort1,vec1)
    panel(1) = vort1

  ; ERA5 reanalysis
    vort3    = gsn_csm_contour_map(wks,vort_pl_era,opts_vort_era)
    vec3     = gsn_csm_vector(wks,u_pl_era,v_pl_era,opts_vec)
    overlay(vort3,vec3)  
    panel(2) = vort3

;================================
; Overlay boxes onto the plots
;================================

    dres                   = True
    dres@gsLineThicknessF  = 3.0
    dres@gsLineDashPattern = 1
    dres@gsLineColor       = "black"

    lines = new(dimsizes(panel),graphic)

    do n = 0, dimsizes(panel)-1
      lines(n) = gsn_add_polyline(wks,panel(n),lon_box(n,:),lat_box(n,:),dres)
    end do

;=================================
; Additional panel plot options
;=================================

    p	     = pres_cp(ilev)
    
    opts_vortP                       = True
    opts_vortP@gsnFrame              = False
    opts_vortP@gsnPanelLabelBar      = True
    opts_vortP@lbTitleString         = p+" hPa absolute vorticity (10~S~-4~N~ s~S~-1~N~)"
    opts_vortP@lbTitlePosition       = "Bottom"
    opts_vortP@lbTitleFontHeightF    = 0.015
    opts_vortP@pmLabelBarWidthF      = 0.6
    opts_vortP@pmLabelBarHeightF     = 0.15
    opts_vortP@lbLabelFontHeightF    = 0.0125
    opts_vortP@txString 	     = ""
;    opts_vortP@gsnPanelFigureStrings = (/"a) CP MetUM: "+ens0+" (4.4 km)", \
;    				         "b) Global MetUM: "+ens1+" (20 km)", \
;                                         "d) ERA5 reanalysis (31 km)"/)
    opts_vortP@gsnPanelFigureStrings = (/"a)", "b)", "c)"/)
    opts_vortP@gsnMaximize           = True
    opts_vortP@gsnPanelTop           = 0.90
    opts_vortP@gsnPanelBottom        = 0.10
    opts_vortP@amJust                = "TopLeft"
    opts_vortP@gsnPanelFigureStringsFontHeightF = 0.0125 ; Label size (default 0.01)

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    gsn_panel(wks,panel,(/1,3/),opts_vortP) ; Draw as a single plot
    frame(wks)

  ; Tidy up (delete variables from within pressure level loop before next iteration)
    delete([/vort_pl_cp, vort_pl_gl, vort_pl_era/])
    delete([/u_pl_cp, u_pl_gl, u_pl_era, v_pl_cp, v_pl_gl, v_pl_era/])
    delete([/geo_pl_era0,dims_h_era0,geo_1d_era0,inds_h_era0/])
    delete([/geo_pl_gl0,dims_h_gl0,geo_1d_gl0,inds_h_gl0/])
    delete([/geo_pl_cp0,dims_h_cp0,geo_1d_cp0,inds_h_cp0/])
    delete([/lon_box,lat_box/])

    end do ; End pressure level loop (do ilev = 0, dimsizes(p_cp)-1)

  ; Size of grid may change from one time to next, so delete after each iteration
    delete([/lon1_cp,lat1_cp,lon_cp,lat_cp,lon1_gl,lat1_gl,lon_gl,lat_gl/])
    delete([/u_cp,v_cp,vort_cp,t_cp,z_cp,f_cp0,f_cp,avo_cp/])
    delete([/u_gl,v_gl,vort_gl,t_gl,z_gl,f_gl0,f_gl,avo_gl/])
    delete([/u_era,v_era,vort_era,t_era,z_era,f_era0,f_era,avo_era/])

  end do     ; End input file loop (do it = 0, 96)

end 