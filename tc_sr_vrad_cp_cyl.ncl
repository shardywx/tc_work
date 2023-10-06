; Script to compute the flux form of the vorticity equation 
; and to calculate changes in circulation due to its different terms  

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

; List all files to read in and analyse
  diri        = "$sam/um/cp"
  fili_prefix = "20160704T0000Z_PHI_4p4_L80_ra1t_pb"
  fili_a      = systemfunc("cd "+diri+" ; ls "+fili_prefix+"*.nc")
  fili        = diri+"/"+fili_a
  numINPUT    = dimsizes(fili)
  numTIMES    = numINPUT * 2

  ct          = 0
  centre      = new((/numTIMES,2/),float)
  x_cyc       = new((/numTIMES/),float)   ; Cyclone longitude
  y_cyc       = new((/numTIMES/),float)   ; Cyclone latitude
  u_cyc       = new((/numTIMES/),float)   ; Zonal cyclone velocity
  v_cyc       = new((/numTIMES/),float)   ; Meridional cyclone velocity
  vel_cyc     = new((/numTIMES/),float)   ; Cyclone velocity (vector)  

;===============================
; Start loop over input files
;===============================

  do nf = 0, numINPUT-1
  f = addfile(fili(nf),"r")

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 300000000
  end setvalues

;==================================
; Get the variables we will need
;==================================

  time  = f->t            ; Times in file (2 - every hour)
  times = dimsizes(time)

  lon1  = f->longitude_1  ; longitude (1098 points --> 109.04 to 152.92 degrees E)
  lat1  = f->latitude_1   ; latitude  (810 points --> 1.8 to 34.16 degrees N)
  lon   = f->longitude    ; longitude (1098 points --> 109.02 to 152.90 degrees E)
  lat   = f->latitude     ; latitude (811 points --> 1.78 to 34.18 degrees N)

  pres  = f->p            ; Pressure (levels)
  plevs = dimsizes(pres)  ; Size of pressure level array

  u     = f->u            ; Zonal velocity (m/s) [time | pres | lat | lon]
  v     = f->v            ; Meridional velocity (m/s) [time | pres | lat | lon]
  vort  = f->rvor         ; Relative vorticity (/s) [time | pres | lat | lon]
  z     = f->ht           ; Geopotential height (m) [time | pres | lat1 | lon1]

  dx    = lon(1) - lon(0) ; Grid spacing (longitude)
  dy    = lat(1) - lat(0) ; Grid spacing (latitude)

; Convert to units we want
  z     = z/10
  z@units = "dam" ; Convert to decametres

  vort  = vort * (10 ^ -6) ; Convert to /s

 ; Degrees to radians (for calculation below)
   pii = 3.14159265
   lat_rad = (pii/180)*lat    ; u,v grid
   lon_rad = (pii/180)*lon

   lat1_rad = (pii/180)*lat1  ; height, temperature grid
   lon1_rad = (pii/180)*lon1

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

   do it = 0, 1

     date_str(it) = sprinti("%0.2iUTC ", hour(it)) + \
                    sprinti("%0.2i ", day(it)) \
                    + month_abbr(month(it))
     time_str(it) = sprinti("%0.2i UTC ", hour(it)) + \
                    sprinti("%0.2i ", day(it)) \
                    + month_abbr(month(it))
     out_str(it)  = sprinti("%0.2i", day(it)) + \
                    month_abbr(month(it)) + \
                    "_" + sprinti("%0.2iZ", hour(it))

     print("Working on time: "+time_str(it))

;================================================================================
; Find cyclone centre on subset of global grid (use both height and vorticity)
;================================================================================

    ; Choose a level in the lower troposphere (950 hPa)
      ilev = 1

      lat_0 = 5
      lat_1 = 30
      lon_0 = 115
      lon_1 = 148

    ; First, concentrate on a subset of the global grid
      u_plane    = u(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
      v_plane    = v(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
      geo_plane  = z(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
      vort_plane = vort(it,ilev,{lat_0:lat_1},{lon_0:lon_1})

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
      geo_min       = min(geo_plane)
      geo_smth      = smth9_Wrap(geo_plane, 0.5, 0.5, True)
      geo_min_smth  = min(geo_smth)

      dims_h        = dimsizes(geo_plane)
      geo1d         = ndtooned(geo_plane)
      inds_h        = ind_resolve(minind(geo1d),dims_h)
      geo1d_sm      = ndtooned(geo_smth)
      inds_h_sm     = ind_resolve(minind(geo1d_sm),dims_h)

      lat_min950 = lat_0 + (dy * inds_h(0,0))
      lon_min950 = lon_0 + (dx * inds_h(0,1))
      print("Cyclone centre (hgt): "+lat_min950+" degrees N, "+lon_min950+" degrees E")

   ;========================================================
   ; EMPLOY A SAFETY NET IF WE IDENTIFY THE WRONG CYCLONE
   ;========================================================
     if (ct.eq.0) then
      centre(ct,0) = lat_min950
      centre(ct,1) = lon_min950

     else
      lt0   = centre(ct-1,0)         ; Previous TC latitude
      ln0   = centre(ct-1,1)         ; Previous TC longitude
      print("lt0 = "+lt0+" ; ln0 = "+ln0)
      d_lat = abs(lat_min950-lt0)       ; Change in latitude (t1 - t0)
      d_lon = abs(lon_min950-ln0)       ; Change in longitude
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

       lat_min950 = lt1 + (dy * inds_h0(0,0))
       lon_min950 = ln1 + (dy * inds_h0(0,1))
       print("New centre (hgt): "+lat_min950+" degrees N, "+lon_min950+" degrees E")
       centre(ct,0) = lat_min950
       centre(ct,1) = lon_min950

       delete([/geo_plane0, geo_min0, dims_h0, geo_1d0, inds_h0/])

      else
       print("Cyclone centre (hgt): "+lat_min950+" degrees N, "+lon_min950+" degrees E")
       centre(ct,0) = lat_min950
       centre(ct,1) = lon_min950
      end if

     end if

    ; Tidy up before looping over pressure levels
      delete([/vort1d, inds, vort1d_sm, vort_smth, inds_sm/])
      delete([/u_plane, v_plane, vort_plane, geo_plane/])

      ct = ct + 1

    end do    ; End of time loop (do it = 0, times-1)

  ; Tidy up to avoid array dimension errors [but leave alone on final iteration]
    if (nf.lt.numINPUT-1)
      delete([/time,times,z,u,v,vort/])
      delete([/utc_date,year,month,day,hour,\
              minute,second,date_str,time_str,out_str/])
    end if

   end do     ; End of input file loop (do nf = 0, numINPUT-1)

;============================================================================
; Calculate distance between two points --> then calculate system velocity
;============================================================================

   print("Starting system velocity calculations...")

 ; Loop to calculate system velocity from hourly cyclone position [lat,lon]
   do iang = 1, numTIMES-1
    x_cyc(iang)   = gc_latlon(centre(iang-1,0),centre(iang-1,1),\
                       centre(iang-1,0),centre(iang,1),2,3)
    y_cyc(iang)   = gc_latlon(centre(iang-1,0),centre(iang-1,1),\
                       centre(iang,0),centre(iang-1,1),2,3)
    u_cyc(iang)   = x_cyc(iang)/(3600)
    v_cyc(iang)   = y_cyc(iang)/(3600)
    vel_cyc(iang) = sqrt(u_cyc(iang)^2.0 + v_cyc(iang)^2.0)
   end do

 ; Set initial system velocity [0] equal to that at the second time [1]
   u_cyc(0)   = u_cyc(1)
   v_cyc(0)   = v_cyc(1)
   vel_cyc(0) = vel_cyc(1)

 ; Add metadata to new arrays
   copy_VarAtts(u,u_cyc(:))
   u_cyc@units         = "m s~S~-1~N~"
   u_cyc@description   = "Cyclone zonal velocity"

   copy_VarAtts(v,v_cyc(:))
   v_cyc@units         = "m s~S~-1~N~"
   v_cyc@description   = "Cyclone meridional velocity"

   vel_cyc@units       = "m s~S~-1~N~"
   vel_cyc@description = "Cyclone speed"

   u_cyc(:) = -(u_cyc(:)) ; Easterly motion is negative [comment if not needed]

 ; Tidy up
   delete([/lat,lon,lat1,lon1,pres,plevs,u,v,z,vort/])
   delete([/time,times,utc_date,year,month,day,hour,\
            minute,second,date_str,time_str,out_str/])

;=================================
; Now start main part of script
;=================================

; List all files to read in and analyse
  diri	      = "$sam/um/cp"
  fili_prefix = "20160704T0000Z_PHI_4p4_L80_ra1t_pb"
  fili_a      = systemfunc("cd "+diri+" ; ls "+fili_prefix+"*.nc")
  fili	      = diri+"/"+fili_a
  numINPUT    = dimsizes(fili)
  numTIMES    = numINPUT * 2 

  ct	      = 0
  time_arr    = new(numTIMES,string)
  title_arr   = new(numTIMES,string)

;===============================
; Start loop over input files 
;===============================

  do nf = 0, numINPUT-1
  f = addfile(fili(nf),"r")

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 300000000
  end setvalues

;==================================
; Get the variables we will need
;==================================

  time  = f->t            ; Times in file (2 - every hour)
  times = dimsizes(time)

  lon1  = f->longitude_1  ; longitude (1098 points --> 109.04 to 152.92 degrees E)
  lat1  = f->latitude_1   ; latitude  (810 points --> 1.8 to 34.16 degrees N)
  lon   = f->longitude    ; longitude (1098 points --> 109.02 to 152.90 degrees E)
  lat   = f->latitude     ; latitude (811 points --> 1.78 to 34.18 degrees N)

  pres  = f->p		  ; Pressure (levels)
  plevs	= dimsizes(pres)  ; Size of pressure level array

  u     = f->u            ; Zonal velocity (m/s) [time | pres | lat | lon]
  v     = f->v            ; Meridional velocity (m/s) [time | pres | lat | lon]
  vort  = f->rvor         ; Relative vorticity (/s) [time | pres | lat | lon]
  div   = f->div          ; Divergence (/s) [time | pres | lat | lon]

  t     = f->temp	  ; Temperature (K) [time | pres | lat1 | lon1]
  z     = f->ht   	  ; Geopotential height (m) [time | pres | lat1 | lon1]
  w     = f->dz_dt	  ; Vertical velocity (m/s) [time | pres | lat1 | lon1]

  dx	= lon(1) - lon(0) ; Grid spacing (longitude)
  dy	= lat(1) - lat(0) ; Grid spacing (latitude)

 ; Convert to units we want 
   z     = z/10
   z@units = "dam" ; Convert to decametres

   vort  = vort * (10 ^ -6) ; Convert to /s
   div   = div  * (10 ^ -6) ; Convert to /s

 ; Pressure levels (pres)
 ; (0)   1000	  	
 ; (1)   950		
 ; (2)   900		
 ; (3)   850
 ; (4)	 800
 ; (5)	 750		
 ; (6)   700		
 ; (7)   650 		
 ; (8)   600		
 ; (9)   550		
 ; (10)  500 		
 ; (11)  450		
 ; (12)  400		
 ; (13)  350		
 ; (14)  300		
 ; (15)  250		
 ; (16)  200		
 ; (17)  150		
 ; (18)  100

 ; Degrees to radians (for calculation below)
   pii = 3.14159265
   lat_rad = (pii/180)*lat    ; u,v grid
   lon_rad = (pii/180)*lon     

   lat1_rad = (pii/180)*lat1  ; height, temperature grid
   lon1_rad = (pii/180)*lon1

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

   do it = 0, 1

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

     print("Working on time: "+time_str(it))
     title_arr(ct) = "Valid at: "+time_str(it)

;================================================================================
; Find cyclone centre on subset of global grid (use both height and vorticity)
;================================================================================

    ; Choose a level in the lower troposphere (950 hPa)
      ilev = 1

      lat_0 = 5
      lat_1 = 30
      lon_0 = 115
      lon_1 = 148

    ; First, concentrate on a subset of the global grid 
      u_plane    = u(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
      v_plane    = v(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
      geo_plane  = z(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
      t_plane    = t(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
      vort_plane = vort(it,ilev,{lat_0:lat_1},{lon_0:lon_1})	

    ; Storm-relative velocity
      u_plane    = u_plane - u_cyc(ct)
      v_plane    = v_plane - v_cyc(ct)

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
      geo_min_smth  = min(geo_smth)
      
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
      centre(ct,0) = lat_min950
      centre(ct,1) = lon_min950

     else
      lt0   = centre(ct-1,0)         ; Previous TC latitude
      ln0   = centre(ct-1,1)         ; Previous TC longitude
      print("lt0 = "+lt0+" ; ln0 = "+ln0)
      d_lat = abs(lat_min950-lt0)       ; Change in latitude (t1 - t0)
      d_lon = abs(lon_min950-ln0)       ; Change in longitude (note extra minus sign)
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

       lat_min950 = lt1 + (dy * inds_h0(0,0))
       lon_min950 = ln1 + (dy * inds_h0(0,1))
       print("New centre (hgt): "+lat_min950+" degrees N, "+lon_min950+" degrees E")
       centre(ct,0) = lat_min950
       centre(ct,1) = lon_min950

       delete([/geo_plane0, geo_min0, dims_h0, geo_1d0, inds_h0/])

      else
       print("Cyclone centre (hgt): "+lat_min950+" degrees N, "+lon_min950+" degrees E")
       centre(ct,0) = lat_min950
       centre(ct,1) = lon_min950
      end if

     end if

    ; Tidy up before looping over pressure levels
      delete([/vort1d,inds,vort1d_sm,vort_smth,inds_sm/])
      delete([/u_plane,v_plane,vort_plane,t_plane,geo_plane/])

    do ilev = 0, 1 ; plevs-1 ; 1000,950,900,850,800,750,700,650,600,550,...,100

     print("Working on: " + pres(ilev) + " hPa")
     p = pres(ilev)

     u_plane       = u(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
     v_plane       = v(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
     vort_plane	   = vort(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
     geo_plane     = z(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
     t_plane  	   = t(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
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

     delete(vort1d)
     delete(inds)
     delete(vort1d_sm)
     delete(vort_smth)
     delete(inds_sm)

   ; Account for cases in which another maximum is chosen in error 
     lat_diff = abs(lat_max-lat_min950)
     lon_diff = abs(lon_max-lon_min950)

     if (lon_diff.gt.2.0.or.lat_diff.gt.2.0) then
        lat_max = lat_min950
        lon_max = lon_min950
     end if
     print("Cyclone centre: "+lat_max+" degrees N, "+lon_max+" degrees E")

;=================================
; Output file type and location
;=================================

    output = "$sam/nepartak/images/cyl_coords/cp_rmw_vrad"+p+"_"+out_str(it)
    wks = gsn_open_wks(opt,output)

;===========================================
; Switch to cylindrical coordinate system
;===========================================

     radii  = new(25,float)   ; New arrays
     thetas = new(73,float)

     radii      = fspan(0.0,3,25)       ; Array of radial distances (3 deg from centre)
     thetas     = fspan(0.0,360.0,73)   ; Array of azimuth angles (0, 5, 10, ..., 360)
     thetas_rad = (pii/180)*thetas      ; Degrees to radians
     rad_size   = dimsizes(radii)-1     ; Index for later calculation

   ; Create new arrays describing cylindrical coordinate system (2D for now)
     xpos = new((/dimsizes(radii),dimsizes(thetas)/),float)
     ypos = new((/dimsizes(radii),dimsizes(thetas)/),float)

   ; New arrays containing variables we want to plot
     u_int      = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Zonal wind
     v_int      = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Meridional wind
     v_rad   	= new((/dimsizes(radii),dimsizes(thetas)/),float) ; Radial wind
     v_tan   	= new((/dimsizes(radii),dimsizes(thetas)/),float) ; Tangential wind

   ; Latitude/longitude (small grid)
     do irad = 0, dimsizes(radii)-1
       r = radii(irad)
       do iang = 0, dimsizes(thetas)-1
         theta = thetas_rad(iang)
           xpos(irad,iang) = lon_max + r*cos(theta)
           ypos(irad,iang) = lat_max + r*sin(theta)
       end do
     end do

  ; Loop over azimuth angles to fill entire circular array of data points
    do iang = 0, dimsizes(thetas_rad)-1
      u_int(:,iang)    = linint2_points(lon_sub,lat_sub,u_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      v_int(:,iang)    = linint2_points(lon_sub,lat_sub,v_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
    end do

  ; Calculate radial and tangential windspeed at each point
    do iang = 0, dimsizes(thetas)-1
      theta = thetas_rad(iang)
        v_rad(:,iang) = u_int(:,iang)*cos(theta) + v_int(:,iang)*sin(theta)
        v_tan(:,iang) = -u_int(:,iang)*sin(theta) + v_int(:,iang)*cos(theta)
    end do

   ; Add variable metadata
     copy_VarAtts(u_plane,u_int)
     u_int@description = "Zonal velocity"
     u_int@units       = "m s~S~-1~N~"

     copy_VarAtts(v_plane,v_int)
     v_int@description = "Meridional velocity"
     v_int@units       = "m s~S~-1~N~"

     v_rad@description = "Radial velocity"
     v_rad@units       = "m s~S~-1~N~"

     v_tan@description = "Tangential velocity"
     v_tan@units       = "m s~S~-1~N~"

     xpos!0   = "rad"
     xpos!1   = "azi"
     xpos&rad = radii
     xpos&azi = thetas

     ypos!0   = "rad"
     ypos!1   = "azi"
     ypos&rad = radii
     ypos&azi = thetas

     u_int!0   = "rad"
     u_int!1   = "azi"
     u_int&rad = radii
     u_int&azi = thetas

     v_int!0   = "rad"
     v_int!1   = "azi"
     v_int&rad = radii
     v_int&azi = thetas

     v_tan!0 = "rad"
     v_tan!1 = "azi"
     v_tan&rad = radii
     v_tan&azi = thetas

     v_rad!0 = "rad"
     v_rad!1 = "azi"
     v_rad&rad = radii
     v_rad&azi = thetas

   ; Windspeed (zonal/meridional)
     speed = sqrt(u_int^2.0 + v_int^2.0)
     copy_VarCoords(u_int,speed)
     copy_VarAtts(u_int,speed)
     copy_VarMeta(u_int,speed)
     speed@units = "m s~S~-1~N~"
     speed@standard_name = "horizontal_velocity"
     speed@long_name     = "Horizontal velocity"
     speed@title         = "Horizontal velocity"
     speed@name          = "wind velocity"
     speed@description   = "Horizontal velocity"

;===============================================================
; Calculate location of windspeed maximum on cylindrical grid
;===============================================================

     speed_max  = max(speed)                           ; Calculate maximum windspeed
     dims_s     = dimsizes(speed)                      ; Size of windspeed array
     speed_1d   = ndtooned(speed)                      ; Reshape windspeed array to 1D
     inds_s     = ind_resolve(maxind(speed_1d),dims_s) ; Calculate location of maximum

;==========================
; Options for plotting 
;==========================

    gsn_define_colormap(wks,"BlueDarkRed18")

  ; Wind speed
    opts_vel                              = True
    opts_vel@sfXArray                     = xpos
    opts_vel@sfYArray                     = ypos
    opts_vel@cnFillOn                     = True
    opts_vel@cnLineLabelInterval          = 2.0
    opts_vel@cnLineLabelFontHeightF       = 0.012
    opts_vel@cnLineLabelBackgroundColor   = "transparent"
    opts_vel@cnLineLabelPlacementMode     = "constant"
    opts_vel@cnLinesOn                    = False ; Contour lines off
    opts_vel@cnInfoLabelOn                = False
    opts_vel@cnLevelSelectionMode         = "ExplicitLevels"
    opts_vel@cnLevels                     = (/-20., -15., -10., -5., -2., -1., \
                                                1., 2., 5., 10., 15., 20./)
    opts_vel@cnFillColors                 = (/2,4,6,7,8,9,0,12,13,14,15,17,19/)

;    opts_vel@cnLevels                     = (/-10., -5., -3., -2., -1., -0.5, -0.2, \
;                                               0, 0.2, 0.5, 1., 2., 3., 5., 10./)
;    opts_vel@cnFillColors                 = (/2,4,5,6,7,8,9,0,0,12,13,14,15,16,17,19/)

    opts_vel@gsnPaperOrientation          = "landscape"
    opts_vel@tiMainString                 = title_arr(ct)
    opts_vel@tiMainFontHeightF		  = 0.0125
    opts_vel@gsnLeftString                = ""
    opts_vel@gsnRightString               = ""
    opts_vel@gsnContourZeroLineThicknessF = 3.0 ; Double the thickness of the zero line

    opts_vel@pmLabelBarWidthF	    	  = 0.6
    opts_vel@pmLabelBarHeightF	    	  = 0.15
    opts_vel@lbLabelFontHeightF	  	  = 0.01
    opts_vel@lbPerimOn		          = False
    opts_vel@lbLabelFont		  = "Helvetica"
    opts_vel@lbTitleOn		    	  = True
    opts_vel@lbTitleString = p+" hPa radial windspeed (m s~S~-1~N~)"
    opts_vel@lbTitleFontHeightF	     	  = 0.015
    opts_vel@lbTitleDirection	    	  = "Across"

  ; Additional plotting resources 
    opts_vel@mpDataBaseVersion	          = "Ncarg4_1" ; More recent database 
    opts_vel@mpDataSetName		  = "Earth..4" ; High resolution 
    opts_vel@mpOutlineBoundarySets 	  = "National" ; National borders
    opts_vel@mpGeophysicalLineColor	  = "black"    ; Colour borders black
    opts_vel@mpGeophysicalLineThicknessF  = 1.0        ; Border line thickness
    opts_vel@mpGridAndLimbOn	    	  = False      ; Turn on lat/lon lines
    opts_vel@pmTickMarkDisplayMode	  = "Always"   ; Turn on map tickmarks
    opts_vel@tmXBMajorLengthF	    	  = 0.005      ; Change tickmark length
    opts_vel@tmXTOn			  = "False"    ; No tickmarks on top x-axis
    opts_vel@tmYROn			  = "False"    ; No tickmarks on right y-axis
    opts_vel@gsnMaximize	 	  = True       ; Maximise plot size
    opts_vel@gsnAddCyclic                 = False

    radius				  = radii(rad_size) ; Size of plot (lat/lon)

    opts_vel@mpLimitMode		  = "Corners"
    opts_vel@mpLeftCornerLatF             = lat_max-radius
    opts_vel@mpLeftCornerLonF             = lon_max-radius
    opts_vel@mpRightCornerLatF            = lat_max+radius
    opts_vel@mpRightCornerLonF            = lon_max+radius
    opts_vel@gsnDraw		          = False      ; Do not draw the plot
    opts_vel@gsnFrame		          = False      ; Do not advance the frame

;=============
; Plot data 
;=============

    plot_vort  = gsn_csm_contour_map(wks,vort_plane,opts_vort) ; Relative vorticity
    plot_vec   = gsn_csm_contour(wks,u_plane,v_plane,opts_vec) ; Horizontal wind vectors
    plot_vvel  = gsn_csm_contour(wks,w_plane,opts_vvel)	       ; Vertical velocity
    overlay(plot_vort,plot_vec)				       
    overlay(plot_vort,plot_vvel)

    plot_speed = gsn_csm_contour_map(wks,v_rad,opts_vel) ; Windspeed

;==============================
; Add radial circles to plot
;==============================

      DEGTORAD = 0.017453292519943

    ; Get customizations for circles.
      radius  = radii(rad_size)
      xcenter = lon_max
      ycenter = lat_max
      outthck = 2
      inthck = 1
      spacing = 0.5
      dpattrn = 2
      outcolr = 1
      incolr = 1

    ;=============================================
    ; Construct and attach outer circle to plot
    ;=============================================

    ; Calculate arrays for outer circle (similar to what I have done)
      degrees = ispan(0,360,5)
      xcos    = cos(DEGTORAD * degrees)
      xsin    = sin(DEGTORAD * degrees)

    ; Array of points (x,y) representing the outer circle (radius = 25)
      xc      = xcenter + radius  * xcos
      yc      = ycenter + radius  * xsin

    ; Resources for outer circle
      lnres                   = True
      lnres@gsLineThicknessF  = outthck
      lnres@gsLineColor       = outcolr
      lnres@gsLineDashPattern = 0

    ; Attach outer circle to plot
      plot_outer     = gsn_add_polyline(wks,plot_speed,xc,yc,lnres)

      if(spacing.gt.radius) then
        print("add_radial_circles: spacing is > radius, can't draw inner circles.")
        return
      end if

    ;==================================================
    ; Now construct and attach inner circles to plot
    ;==================================================

    ; Draw inner circles if desired ('spacing' controls the number of inner circles)
      if(spacing.gt.0) then
        count = 0
              size = toint((radius-spacing)/spacing)
              plot_inner = new(size,graphic)
        do r = spacing,radius-spacing,spacing

        ; Calculate arrays for inner circle.
          xc = xcenter + (r * xcos)
          yc = ycenter + (r * xsin)

        ; Resources for inner circle
          delete(lnres@gsLineThicknessF)   ; Delete in case setting to
          delete(lnres@gsLineColor)        ; new type.
                delete(lnres@gsLineDashPattern)
          lnres@gsLineColor       = incolr
          lnres@gsLineThicknessF  = inthck
          lnres@gsLineDashPattern = dpattrn

	  ; Attach inner circle to plot
          plot_inner(count) = gsn_add_polyline(wks,plot_speed,xc,yc,lnres)
          count = count + 1
        end do
      end if

    ; Tidy up
      delete([/lnres, degrees, xcos, xsin, xc, yc/])

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
      xlines  = 360./angle
      nlines  = toint(ceil(xlines)) ; Smallest integer greater than/equal to xlines
      if(xlines.ne.nlines) then
        degrees = fspan(0,360,nlines)
      else
        degrees = fspan(0,360-angle,nlines) ; Contains angles (0, 45, 90, etc)
      end if

      nlines  = dimsizes(degrees)
      plot_lines = new(nlines,graphic)
      do i = 0,nlines-1
        xl = (/xcenter, xcenter + (radius  * cos(DEGTORAD * degrees(i)))/)
        yl = (/ycenter, ycenter + (radius  * sin(DEGTORAD * degrees(i)))/)

      ; Attach line to plot
        plot_lines(i) = gsn_add_polyline(wks,plot_speed,xl,yl,lnres)

      end do

    ; Tidy up
      delete(degrees)

;=============================================
; Add marker illustrating windspeed maximum
;=============================================

    ; Add circle indicating the radius of maximum windspeed
      rmw                     = inds_s(0,0)     ; Radial distance from TC centre
      thta                    = inds_s(0,1)     ; Azimuthal position (angle)
      xcentre                 = lon_max         ; TC centre (lon)
      ycentre                 = lat_max         ; TC centre (lat)

    ; Calculate circular array of points (every 5 deg)
      degrees                 = fspan(0.0,360.0,73)
      xcos                    = cos(pii/180 * degrees) ; cosine of angles in array (rad)
      xsin                    = sin(pii/180 * degrees) ; sine of angles in array (rad)

    ; First add marker corresponding to exact location of maximum
      marker = new(1,graphic)
      mres   = True
      mres@gsMarkerIndex = 5        ; Cross
      mres@gsMarkerThicknessF = 3.0 ; Marker thickness
      mres@gsMarkerSizeF = 16.0     ; Marker size
      mres@gsMarkerColor = "black"  ; Marker colour

    ; Calculate lat/lon of windspeed maximum on cylindrical grid
      cpos_x = xcentre + (radii(rmw) * xcos(thta))
      cpos_y = ycentre + (radii(rmw) * xsin(thta))

;      marker = gsn_add_polymarker(wks,plot_speed,cpos_x,cpos_y,mres)
;      delete(mres)

;================================================================
; Finally, draw circle to highlight the radius of maximum wind
;================================================================

    ; Additional options for circle representing maximum wind
      thck                    = 5
      color                   = 1

    ; Array of points (x,y) representing the circle
      xc                      = xcentre + radii(rmw) * xcos ; origin + [r * cos(theta)]
      yc                      = ycentre + radii(rmw) * xsin ; origin + [r * sin(theta)]

    ; Resources for circle indicating the radius of maximum wind
      lnres                   = True
      lnres@gsLineThicknessF  = thck
      lnres@gsLineColor       = color  ; Option only recognises numbers (color = 0)
      lnres@gsLineDashPattern = 0

    ; Attach circle to plot
      plot_rmw                = gsn_add_polyline(wks,plot_speed,xc,yc,lnres)
      delete([/lnres, xc, yc, degrees, xcos, xsin, rmw, xcentre, ycentre, thck, color/])

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    draw(plot_speed)
    frame(wks)

    ; Tidy up before next iteration
      delete([/vort_plane,u_plane,v_plane,t_plane,geo_plane,lat_plane,lon_plane/])
      delete([/u_int,v_int,v_tan,v_rad,speed/])

     end do	   ; End pressure level loop (do ilev = 0, plevs-1)

    end do     ; End time loop (do it = 0, 1)

  ; Tidy up before next iteration
    delete([/date_str, time_str, out_str/])
    delete([/time,times,u,v,z,t,w,vort,pres,plevs,div/])

  end do     ; End input file loop

end 