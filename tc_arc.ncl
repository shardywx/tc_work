; Produce Hovmoller plots of momentum budget terms 

; Run using:

; ncl 'opt="x11"' 'ens0="em02"' clev=2 dist=3.0 'sy="lee"' ts0=30 tf0=72 tc_arc.ncl

; 'opt'   = output file format ("pdf" or "x11")
; 'ens0'  = ensemble member (em00 -- em11)
; 'clev'  = pressure level for tangential wind panel (0 = 950, 1 = 900, 2 = 850, 3 = 800)
; 'dist'  = size of box following storm (degrees)
; 'ts0'	  = starting time index (0–119)
; 'tf0'   = ending time index (0–119)
; 'sy'	  = machine (ARCHER [ar] or LEEDS [lee]) for input/output options

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"

begin

;=====================================
; Calculate storm-relative velocity
;=====================================

   if (sy.eq."lee") then
     diri = "/nfs/a37/earshar/ncl_scripts/nepartak/archer_march2018/"
   else if (sy.eq."ar") then
     diri = "$sam/nepartak/images/cart_coords/cp_um/ens/"
   end if
   end if 
 
   numSIMS     = 12                                ; Number of simulations
   numTIMES    = 120                               ; Number of times in each simulation
   centre      = new((/numTIMES,2/),float)         ; Cyclone track in all simulations

 ; Read in values from ensemble members

   centre(:,0)  = asciiread(diri+"cp_"+ens0+"_lat.txt",(/numTIMES,1/),"float")
   centre(:,1)  = asciiread(diri+"cp_"+ens0+"_lon.txt",(/numTIMES,1/),"float")

;============================================================================
; Calculate distance between two points --> then calculate system velocity
;============================================================================

   print_clock("Starting system velocity calculations...")

 ; Define new arrays (cyclone position, speed)
   x_cyc   = new((/numTIMES/),float)
   y_cyc   = new((/numTIMES/),float)
   u_cyc   = new((/numTIMES/),float)
   v_cyc   = new((/numTIMES/),float)
   vel_cyc = new((/numTIMES/),float)

 ; Loop to calculate system velocity from 6-hourly cyclone position [lat,lon]
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
   x_cyc(0)   = 0
   y_cyc(0)   = 0
   u_cyc(0)   = u_cyc(1)
   v_cyc(0)   = v_cyc(1)
   vel_cyc(0) = vel_cyc(1)

 ; Add metadata to new arrays
   u_cyc@units         = "m s~S~-1~N~"
   u_cyc@description   = "Cyclone zonal velocity"

   vel_cyc@units       = "m s~S~-1~N~"
   vel_cyc@description = "Cyclone speed"

   u_cyc(:) = -(u_cyc(:)) ; Easterly motion is negative

;================================================================================
; Create lat/lon arrays containing grid subsets at each time, following the TC
;================================================================================
  lat_arr = new((/numTIMES,2/),float)
  lon_arr = new((/numTIMES,2/),float)

; Lat/lon arrays for ensemble simulations
  lat_arr(:,0) = centre(:,0)-dist
  lat_arr(:,1) = centre(:,0)+dist
  lon_arr(:,0) = centre(:,1)-dist
  lon_arr(:,1) = centre(:,1)+dist

;==================================================================================
;;;;;;;;;;;;;;;;;;;;;;;;;;; Start of main part of script ;;;;;;;;;;;;;;;;;;;;;;;;;;
;==================================================================================

  print_clock("Starting main part of script...")

; Prelash
  time_arr    = new((/numTIMES/),string)
  title_arr   = new((/numTIMES/),string)

;========================================
; Define cylindrical coordinate arrays
;========================================

  thetas     = new(73,float)
  pii        = 3.14159265
  r_e	     = 6.371e6		     ; Earth's radius (metres)
  radii      = fspan(0.0,1.5,31)     ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,73)   ; Array of azimuth angles (0, 5, 10, etc)
  thetas_rad = (pii/180)*thetas      ; Degrees to radians
  rad_size   = dimsizes(radii)-1     ; Index for later calculation
  num_p	     = 19		     ; Size of pressure level subset array (L344-345)

;=====================================================
; Define arrays for later use in radius-time plots
;=====================================================

  vtan_azi = new((/num_p,dimsizes(radii)/),float)  ; Tangential wind      [ {v} ]
  vrad_azi = new((/num_p,dimsizes(radii)/),float)  ; Radial wind          [ {u} ]
  vort_azi = new((/num_p,dimsizes(radii)/),float)  ; Rel. vorticity       [ {vort} ]
  avo_azi  = new((/num_p,dimsizes(radii)/),float)  ; Abs. vorticity       [ {f + vort} ]

;==================================================================
; Define 3D arrays (hold values from above arrays at ALL times)
;==================================================================

; Mean tangential wind 
  vtan_plot0      = new((/num_p,dimsizes(radii),numTIMES/),float)
  vtan_plot0!0	  = "lev"
  vtan_plot0!1    = "rad"
  vtan_plot0&rad  = radii
  vtan_plot0!2    = "time"

; Mean radial wind 
  vrad_plot0      = new((/num_p,dimsizes(radii),numTIMES/),float)
  vrad_plot0!0	  = "lev"
  vrad_plot0!1    = "rad"
  vrad_plot0&rad  = radii
  vrad_plot0!2    = "time"

; Mean vertical velocity
  vvel_plot0      = new((/num_p,dimsizes(radii),numTIMES/),float)
  vvel_plot0!0    = "lev"
  vvel_plot0!1    = "rad"
  vvel_plot0&rad  = radii
  vvel_plot0!2    = "time"

; Mean relative vorticity
  vort_plot0 	  = new((/num_p,dimsizes(radii),numTIMES/),float)
  vort_plot0!0    = "lev"
  vort_plot0!1    = "rad"
  vort_plot0&rad  = radii
  vort_plot0!2    = "time"

; Mean absolute vorticity
  avo_plot0       = new((/num_p,dimsizes(radii),numTIMES/),float)
  avo_plot0!0	  = "lev"
  avo_plot0!1     = "rad"
  avo_plot0&rad   = radii
  avo_plot0!2     = "time"
;===================================================
; Add metadata for azimuthally averaged variables
;===================================================

  vtan_azi!0   	       = "lev"
  vtan_azi!1           = "rad"
  vtan_azi&rad         = radii
  vtan_azi@description = "Azimuthally averaged tangential velocity"
  vtan_azi@units       = "m s~S~-1~N~"

  vrad_azi!0           = "lev"
  vrad_azi!1           = "rad"
  vrad_azi&rad         = radii
  vrad_azi@description = "Azimuthally averaged radial velocity"
  vrad_azi@units       = "m s~S~-1~N~"

  avo_azi!0            = "lev"
  avo_azi!1            = "rad"
  avo_azi&rad          = radii
  avo_azi@description  = "Azimuthally averaged absolute vorticity"
  avo_azi@units        = "s~S~-1~N~"

  vort_azi!0           = "lev"
  vort_azi!1           = "rad"
  vort_azi&rad         = radii
  vort_azi@description = "Azimuthally averaged relative vorticity"
  vort_azi@units       = "s~S~-1~N~"

; Create time array for later plotting 
  time_int = new(120,"integer")
  time_int = ispan(0,119,1)
  time_int!0 = "time"

;=========================
; Start multiple loops
;=========================

    print_clock("Working on ensemble member: "+ens0)

  ; 00Z 4 Jul [00], 06Z 4 Jul [01], 12Z 4 Jul [02], 18Z 4 Jul [03], 00Z 5 Jul [04],
  ; 06Z 5 Jul [05], 12Z 5 Jul [06], 18Z 5 Jul [07], 00Z 6 Jul [08], 06Z 6 Jul [09],
  ; 12Z 6 Jul [10], 18Z 6 Jul [11], 00Z 7 Jul [12], 06Z 7 Jul [13], 12Z 7 Jul [14],
  ; 18Z 7 Jul [15], 00Z 8 Jul [16], 06Z 8 Jul [17], 12Z 8 Jul [18], 18Z 8 Jul [19],
  ; 00Z 9 Jul [20]

  ; ct = 0   [00], ct = 5  [01], ct = 11  [02], ct = 17  [03], ct = 23  [04],
  ; ct = 29  [05], ct = 35 [06], ct = 41  [07], ct = 47  [08], ct = 53  [09],
  ; ct = 59  [10], ct = 65 [11], ct = 71  [12], ct = 77  [13], ct = 83  [14],
  ; ct = 89  [15], ct = 95 [16], ct = 101 [17], ct = 107 [18], ct = 113 [19],
  ; ct = 119 [20]

  ; List all files to read in and analyse
  ; Analyse from 12Z 5 -- 6 July [36 h]
    diri        = "$ens/reg_netcdf/20160704T0000Z_ra1t_"

    input_pb    = diri+ens0+"_pb.nc"

    numINPUT    = dimsizes(input_pb)        ; Number of input files [4]
    numTIMES    = 120                       ; Total number of times in all files [36]
    ct          = ts0                       ; Counter variable
    ct0		= 0			    ; Counter variable (independent of position)
    plevs	= 19			    ; Number of pressure levels
    llbox       = toint(dist*50)            ; Size of lat/lon arrays (see line below)
    dsize       = (/llbox,llbox/)           ; Size of lat/lon arrays (depends on 'dist')
    d0          = dsize(0)
    d1          = dsize(1)

    setvalues NhlGetWorkspaceObjectId
      "wsMaximumSize" : 1000000000
    end setvalues

    date_str = new(numTIMES,string)
    time_str = new(numTIMES,string)
    out_str  = new(numTIMES,string)

;===============================
; Start loop over input files 
;===============================

    nf = 0
    f  = addfile(input_pb(nf),"r")         ; Read in 'pb' stream [u,v,w,vort,div]

;==================================
; Get the variables we will need
;==================================

   do it = ts0, tf0		; Loop over times in file

    time  = f->t(it)            ; Read in one time interval
    times = dimsizes(time)

;=====================================================
; Create correct date strings for each output time
;=====================================================

   month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
                     "Jul","Aug","Sep","Oct","Nov","Dec"/)

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

   date_str(ct0) = sprinti("%0.2iUTC ", hour) + \
                   sprinti("%0.2i ", day) \
                   + month_abbr(month)
   time_str(ct0) = sprinti("%0.2i UTC ", hour) + \
                   sprinti("%0.2i ", day) \
                   + month_abbr(month)
   out_str(ct0)  = sprinti("%0.2i", day) + \
                   month_abbr(month) + \
                   "_" + sprinti("%0.2iZ", hour)

   time_arr(ct)  = out_str(ct0)

   print_clock("Working on time: "+time_str(ct0))
   print_clock("ct = "+ct+" ; ct0 = "+ct0)
   ctt = ct + 1
   title_arr(ct) = "Valid at: "+time_str(ct0)+" (T+"+ctt+")"

;=================================================
; Continue reading in variables from input file
;=================================================

    t0    = lat_arr(ct,0)
    t1    = lat_arr(ct,1)
    n0    = lon_arr(ct,0)
    n1    = lon_arr(ct,1)
    print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

    lon1  = f->longitude_1({n0:n1})  ; longitude [xxx grid points]
    lat1  = f->latitude_1({t0:t1})   ; latitude  [xxx grid points]
    lon   = f->longitude({n0:n1})    ; longitude [xxx grid points]
    lat   = f->latitude({t0:t1})     ; latitude  [xxx grid points]

    pres  = f->p		     ; Pressure levels
    plevs = dimsizes(pres)  	     ; Size of pressure level array
    pres@units = "hPa"

  ; Read in horizontal velocity, relative vorticity and geopotential height
    u       = f->u(it,:,{t0:t1},{n0:n1})     ; Zonal velocity (m/s)
    v       = f->v(it,:,{t0:t1},{n0:n1})     ; Meridional velocity (m/s)
    pv	    = f->pv(it,:,{t0:t1},{n0:n1})	; Potential vorticity 
    vort    = f->rvor(it,:,{t0:t1},{n0:n1})  ; Relative vorticity (/s)
    z       = f->ht(it,:,{t0:t1},{n0:n1})    ; Geopotential height (m)
    w	    = f->dz_dt(it,:,{t0:t1},{n0:n1}) ; Vertical velocity (m/s)
    t	    = f->temp(it,:,{t0:t1},{n0:n1})	; Temperature (K)

  ; Convert to units we want
    z       = z/10
    z@units = "dam"			        ; Convert to decametres 
    vort    = vort * (10 ^ -6) 	 	        ; Convert to /s
    pv      = pv * (10 ^ 6)                     ; Convert to PVU [10 ^ 6 K m2 s-1 kg-1]

    dx      = lon(1) - lon(0) ; Grid spacing (longitude)
    dy      = lat(1) - lat(0) ; Grid spacing (latitude)

  ; Pressure levels (pres)
  ; (0) 1000, (1) 950, (2) 900, (3) 850, (4) 800, (5) 750, (6) 700, (7) 650
  ; (8) 600, (9) 550, (10) 500, (11) 450, (12) 400, (13) 350
  ; (14) 300, (15) 250, (16) 200, (17) 150, (18) 100

  ; Coriolis parameter and absolute vorticity
    f1   = coriolis_param(lat)
    f0   = conform_dims(dimsizes(vort),f1,1)
    copy_VarMeta(vort,f0)
    copy_VarCoords(vort,f0)

    avo  = vort + f0
    copy_VarMeta(vort,avo)
    copy_VarAtts(vort,avo)
    copy_VarCoords(vort,avo)
    avo@description = "Absolute vorticity"
    avo@long_name   = "Absolute vorticity"

;===================================================
; Add pressure level metadata for selected arrays 
;===================================================

; 3D arrays holding values for final plotting  
  vtan_plot0&lev = pres
  vrad_plot0&lev = pres
  vort_plot0&lev = pres
  avo_plot0&lev  = pres

; Arrays holding azimuthally-averaged variables
  vtan_azi&lev   = pres
  vrad_azi&lev   = pres
  avo_azi&lev    = pres
  vort_azi&lev   = pres

;===================================================================================
; Find cyclone centre on subset of global grid (use 950 hPa height and vorticity)
;===================================================================================

    ; Choose a level in the lower troposphere (950 hPa)
      ilev = 0

    ; First, concentrate on a subset of the global grid 
      u_plane    = u(ilev,:,:)
      v_plane    = v(ilev,:,:)
      geo_plane  = z(ilev,:,:)
      vort_plane = vort(ilev,:,:)

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
      dims_h	    = dimsizes(geo_plane)
      geo1d	    = ndtooned(geo_plane)
      inds_h	    = ind_resolve(minind(geo1d),dims_h)

      lat_min950 = t0 + (dy * inds_h(0,0))
      lon_min950 = n0 + (dx * inds_h(0,1))
      print("Cyclone centre (hgt): "+lat_min950+" degrees N, "+lon_min950+" degrees E")

      delete([/geo1d, inds_h/])

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

       geo_plane0 = z(ilev,{lt1:lt2},{ln1:ln2})          ; Smaller grid (1.0 x 1.0 deg)
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

;===============================
; Loop over pressure levels 
;===============================

    ilev = clev		   ; Pressure level for tangential velocity

     p = pres(ilev)
     print("Working on: " + p + " hPa")

     u_plane0      = u(ilev,:,:)
     v_plane0      = v(ilev,:,:)
     vort_plane	   = vort(ilev,:,:)
     avo_plane	   = avo(ilev,:,:)
     geo_plane     = z(ilev,:,:)

     u_plane	   = u_plane0 - u_cyc(ct)
     v_plane	   = v_plane0 - v_cyc(ct)

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

     lat_max = t0 + (dy * inds(0,0))
     lon_max = n0 + (dx * inds(0,1))

   ;==============================================================
   ; Create arrays for later switch to cylindrical coordinates
   ;==============================================================

   ; u,v grid
     lat_plane  = ind(t0.le.lat.and.lat.le.t1) ; Indices of latitude array
     lon_plane  = ind(n0.le.lon.and.lon.le.n1) ; Indices of longitude array

     lat_size = dimsizes(lat_plane)		  	  	; Array sizes
     lon_size = dimsizes(lon_plane)

     lat_sub = lat(lat_plane(0):lat_plane(lat_size-1)) ; Create subset of points 
     lon_sub = lon(lon_plane(0):lon_plane(lon_size-1)) ; between values given above

   ; w grid
     lat1_plane  = ind(t0.le.lat1.and.lat1.le.t1) ; Indices of latitude array
     lon1_plane  = ind(n0.le.lon1.and.lon1.le.n1) ; Indices of longitude array

     lat1_size = dimsizes(lat1_plane)                             ; Array sizes
     lon1_size = dimsizes(lon1_plane)

     lat1_sub = lat1(lat1_plane(0):lat1_plane(lat1_size-1)) ; Create subset of points
     lon1_sub = lon1(lon1_plane(0):lon1_plane(lon1_size-1)) ; between values given above

     delete([/vort1d,inds,vort1d_sm,vort_smth,inds_sm/])

   ; Account for cases in which another maximum is chosen in error 
     lat_diff = abs(lat_max-lat_min950)
     lon_diff = abs(lon_max-lon_min950)

     if (lon_diff.gt.1.0.or.lat_diff.gt.1.0) then
        lat_max = lat_min950
        lon_max = lon_min950
     end if
     print("Cyclone centre: "+lat_max+" degrees N, "+lon_max+" degrees E")

;===========================================
; Switch to cylindrical coordinate system
;===========================================

   ; Create new arrays describing cylindrical coordinate system (2D for now)
     xpos = new((/dimsizes(radii),dimsizes(thetas)/),float)
     ypos = new((/dimsizes(radii),dimsizes(thetas)/),float)
     arc0 = new((/dimsizes(radii),dimsizes(thetas)/),float)
     arc1 = new((/dimsizes(radii),dimsizes(thetas)/),float)

   ; New arrays containing variables we want to plot
     u_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Zonal velocity
     v_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Meridional velocity
     v_rad    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Radial velocity
     v_tan    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Tangential velocity
     vort_int = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Relative vorticity
     avo_int  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Absolute vorticity

;=========================
; Add variable metadata
;=========================

     copy_VarAtts(u_plane,u_int)
     u_int@description = "Zonal velocity"
     u_int@units       = "m s~S~-1~N~"

     copy_VarAtts(v_plane,v_int)
     v_int@description = "Meridional velocity"
     v_int@units       = "m s~S~-1~N~"

     copy_VarAtts(vort_plane,vort_int)
     vort_int@description = "Relative vorticity"
     vort_int@units       = "s~S~-1~N~"

     copy_VarAtts(avo_plane,avo_int)
     avo_int@description = "Absolute vorticity"
     avo_int@units       = "s~S~-1~N~"

     xpos!0    = "rad"
     xpos!1    = "azi"
     xpos&rad  = radii
     xpos&azi  = thetas

     ypos!0    = "rad"
     ypos!1    = "azi"
     ypos&rad  = radii
     ypos&azi  = thetas

     arc0!0    = "rad"
     arc0!1    = "azi"
     arc0&rad  = radii
     arc0&azi  = thetas

     arc1!0    = "rad"
     arc1!1    = "azi"
     arc1&rad  = radii
     arc1&azi  = thetas

     u_int!0   = "rad"
     u_int!1   = "azi"
     u_int&rad = radii
     u_int&azi = thetas

     v_int!0   = "rad"
     v_int!1   = "azi"
     v_int&rad = radii
     v_int&azi = thetas

     vort_int!0   = "rad"
     vort_int!1   = "azi"
     vort_int&rad = radii
     vort_int&azi = thetas

     avo_int!0    = "rad"
     avo_int!1    = "azi"
     avo_int&rad  = radii
     avo_int&azi  = thetas

     v_tan!0      = "rad"
     v_tan!1      = "azi"
     v_tan&rad    = radii
     v_tan&azi    = thetas

     v_rad!0      = "rad"
     v_rad!1      = "azi"
     v_rad&rad    = radii
     v_rad&azi    = thetas

   ;=========================================================
   ; Loop over azimuth and radius to calculate arc lengths
   ;=========================================================

     do irad = 0, dimsizes(radii)-1
       rd = radii(irad)			; Radius (degrees)
       rr = rd * (pii/180)		; Radius (radians)

       do iang = 0, dimsizes(thetas)-1

         tp1 = min((/iang+1,dimsizes(thetas)-1/))
	 tm1 = max((/iang-1,0/))

         ; Calculate arc length (radians)
           theta = thetas_rad(iang)                                    ; theta (rad)
           dthta = theta - thetas_rad(tm1)                             ; d-theta (rad)

           xpos(irad,iang) = centre(ct,1) + rd*cos(theta)
           ypos(irad,iang) = centre(ct,0) + rd*sin(theta)	

         ; Calculate arc length (metres)
	   d_e   = r_e * (rr * dthta)				       ; Distance (m)
	   arc0(irad,iang)  = d_e				       ; Input into array

       end do

     end do

   ;================================================
   ; Repeat above loop for gc_latlon NCL function
   ;================================================

     do irad = 0, dimsizes(radii)-1
       rd = radii(irad)                 ; Radius (degrees)
       rr = rd * (pii/180)              ; Radius (radians)

       do iang = 0, dimsizes(thetas)-1

         tp1 = min((/iang+1,dimsizes(thetas)-1/))
         tm1 = max((/iang-1,0/))

         ; Calculate distance using notes
           theta = thetas_rad(iang)                                    ; theta (rad)
           dthta = theta - thetas_rad(tm1)                             ; d-theta (rad)

         ; Calculate distance using great circle formula (built-in)
           dist0 = gc_latlon(ypos(irad,iang),xpos(irad,iang),ypos(irad,tp1),\
                             xpos(irad,tp1),2,3)

           arc1(irad,iang)  = dist0                                    ; Input into array

       end do

     end do

     print("manual calculations = "+arc0(20,5:10))
     print("gc_latlon calculations = "+arc1(20,5:10))

  ; Loop over azimuth angles to fill entire circular array of data points
    do iang = 0, dimsizes(thetas_rad)-1
      u_int(:,iang)    = linint2_points(lon_sub,lat_sub,u_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      v_int(:,iang)    = linint2_points(lon_sub,lat_sub,v_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      vort_int(:,iang) = linint2_points(lon_sub,lat_sub,vort_plane,False,\
					xpos(:,iang),ypos(:,iang),0)
      avo_int(:,iang)  = linint2_points(lon_sub,lat_sub,avo_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
    end do

  ; Calculate radial and tangential windspeed at each point
    do iang = 0, dimsizes(thetas)-1
      theta = thetas_rad(iang)
        v_rad(:,iang) = u_int(:,iang)*cos(theta) + v_int(:,iang)*sin(theta)
        v_tan(:,iang) = -u_int(:,iang)*sin(theta) + v_int(:,iang)*cos(theta)
    end do

    ;===================================================
    ; Average variables azimuthally around the storm 
    ;===================================================

      do irad = 0, dimsizes(radii)-1
        r = radii(irad)
        vtan_azi(ilev,irad) = avg(v_tan(irad,:))         ; Tangential wind
	vrad_azi(ilev,irad) = avg(v_rad(irad,:))	 ; Radial wind 
	vort_azi(ilev,irad) = avg(vort_int(irad,:))	 ; Relative vorticity
	avo_azi(ilev,irad)  = avg(avo_int(irad,:))	 ; Absolute vorticity
      end do

    ;===================================
    ; Input values into larger arrays
    ;===================================

    ; Azimuthally-averaged variables
      vtan_plot0(ilev,:,ct)  = (/vtan_azi(ilev,:)/)      ; Tangential wind
      vrad_plot0(ilev,:,ct)  = (/vrad_azi(ilev,:)/)	 ; Radial wind
      vort_plot0(ilev,:,ct)  = (/vort_azi(ilev,:)/)	 ; Relative vorticity
      avo_plot0(ilev,:,ct)   = (/avo_azi(ilev,:)/)	 ; Absolute vorticity

    ; Tidy up
      delete([/xpos,ypos,arc0,arc1/])
      delete([/u_int,v_int,v_tan,v_rad/])
      delete([/u_plane,v_plane,vort_plane,geo_plane/])
      delete([/lat_plane,lon_plane,lat1_plane,lon1_plane/])
      delete([/lat_sub,lat1_sub,lon_sub,lon1_sub/])

      ct  = ct + 1		    ; Counter variable (time)
      ct0 = ct0 + 1		    ; Independent counter variable (time)

  ; Tidy up before next iteration
    delete([/year,month,day,hour,minute,second,utc_date/])
    delete([/lon,lat,lon1,lat1,time,times,u,v,z,t,w,vort,avo,f0,f1,pres,plevs/])

   end do     ; End time loop (do it = 0, times-1)

end 