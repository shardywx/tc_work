; Script to reproduce Figure 6 from Nguyen et al. (2011) for a given ensemble member
; Run using:

; ncl 'opt="x11"' 'ens0="em02"' clev=2 dist=3.0 tc_ens_n11_fig7.ncl

load "$LIB/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

;=====================================
; Calculate storm-relative velocity
;=====================================

   diri        = "$sam/nepartak/images/cart_coords/cp_um/ens/"
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
  radii      = fspan(0.0,1.5,31)     ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,73)   ; Array of azimuth angles (0, 5, 10, etc)
  thetas_rad = (pii/180)*thetas      ; Degrees to radians
  rad_size   = dimsizes(radii)-1     ; Index for later calculation

; Define arrays for later use in tangential velocity cross-section plots
  vtan_azi = new((/dimsizes(radii)/),float) ; Tangential wind
  vrad_azi = new((/dimsizes(radii)/),float) ; Radial wind
  spd_azi  = new((/dimsizes(radii)/),float) ; Vector velocity
  pvrt_azi = new((/dimsizes(radii)/),float) ; Potential vorticity
  vvel_azi = new((/dimsizes(radii)/),float) ; Vertical velocity
  the_azi  = new((/dimsizes(radii)/),float) ; Theta-e

; Define 2-dimensional array to hold all values from 'vtan_azi' at ALL times
  vtan_plot0   = new((/dimsizes(radii),numTIMES/),float)
  vtan_plot0!0 = "rad"
  vtan_plot0!1 = "time"

  vrad_plot0   = new((/dimsizes(radii),numTIMES/),float)
  vrad_plot0!0 = "rad"
  vrad_plot0!1 = "time"

  pvrt_plot0   = new((/dimsizes(radii),numTIMES/),float)
  pvrt_plot0!0 = "rad"
  pvrt_plot0!1 = "time"

  vvel_plot0   = new((/dimsizes(radii),numTIMES/),float)
  vvel_plot0!0 = "rad"
  vvel_plot0!1 = "time"

  the_plot0    = new((/dimsizes(radii),numTIMES/),float)
  the_plot0!0  = "rad"
  the_plot0!1  = "time"

  spd_plot     = new(numTIMES,float)
  spd_plot!0   = "time"

; Add metadata for azimuthally averaged variables
  vtan_azi!0           = "rad"
  vtan_azi&rad         = radii
  vtan_azi@description = "Azimuthally averaged tangential velocity"
  vtan_azi@units       = "m s~S~-1~N~"

  vrad_azi!0           = "rad"
  vrad_azi&rad         = radii
  vrad_azi@description = "Azimuthally averaged radial velocity"
  vrad_azi@units       = "m s~S~-1~N~"

  spd_azi!0            = "rad"
  spd_azi&rad          = radii
  spd_azi@description  = "Azimuthally averaged vector velocity"
  spd_azi@units        = "m s~S~-1~N~"

  pvrt_azi!0           = "rad"
  pvrt_azi&rad         = radii
  pvrt_azi@description = "Azimuthally averaged PV"
  pvrt_azi@units       = "10~S~-6~N~K kg~S~-1~N~ m~S~2~N~ s~S~-1~N~"

  vvel_azi!0            = "rad"
  vvel_azi&rad          = radii
  vvel_azi@description  = "Azimuthally averaged vertical velocity"
  vvel_azi@units        = "m s~S~-1~N~"

  the_azi!0             = "rad"
  the_azi&rad           = radii
  the_azi@description   = "Azimuthally averaged equivalent potential temperature"
  the_azi@units         = "K"

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
    diri        = "$sam/um/cp/ens/20160704T0000Z_ra1t_"
    fili_pb     = diri+ens0+"_pb_"
    input_pb    = (/fili_pb+"00.nc",fili_pb+"01.nc",fili_pb+"02.nc",fili_pb+"03.nc",\
    		    fili_pb+"04.nc",fili_pb+"05.nc",fili_pb+"06.nc",fili_pb+"07.nc",\
    		    fili_pb+"08.nc",fili_pb+"09.nc",fili_pb+"10.nc",fili_pb+"11.nc",\
		    fili_pb+"12.nc",fili_pb+"13.nc",fili_pb+"14.nc",fili_pb+"15.nc",\
		    fili_pb+"16.nc",fili_pb+"17.nc",fili_pb+"18.nc",fili_pb+"19.nc",\
		    fili_pb+"20.nc"/)
    numINPUT    = dimsizes(input_pb)        ; Number of input files [4]
    numTIMES    = 120                       ; Total number of times in all files [36]
    ct          = 23                        ; Counter variable
    ct0		= 0			    ; Counter variable (independent of position)

    setvalues NhlGetWorkspaceObjectId
      "wsMaximumSize" : 1000000000
    end setvalues

    date_str = new(numTIMES,string)
    time_str = new(numTIMES,string)
    out_str  = new(numTIMES,string)

;===============================
; Start loop over input files 
;===============================

;    do nf = 0,  numINPUT-1     ; 00Z 4 -- 9 July
    do nf = 4,  14	       ; 00Z 5 -- 12Z 7 July
    f = addfile(input_pb(nf),"r")

    tm0 = f->t
    tm1	= dimsizes(tm0)

;==================================
; Get the variables we will need
;==================================

   do it = 0, tm1-1

    time  = f->t(it)            ; Times in file (6 per file - hourly)
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

   print("Working on time: "+time_str(ct0))
   title_arr(ct) = "Valid at: "+time_str(ct0)

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

    p_cp  = (/1,2,3/)          	     ; Pressure level(s) to access (950,900,850)

    pres  = f->p(p_cp)         	     ; Pressure (levels)
    plevs = dimsizes(pres)  	     ; Size of pressure level array
    pres@units = "hPa"

  ; Read in horizontal velocity, relative vorticity and geopotential height
    u       = f->u(it,p_cp,{t0:t1},{n0:n1})     ; Zonal velocity (m/s)
    v       = f->v(it,p_cp,{t0:t1},{n0:n1})     ; Meridional velocity (m/s)
    pv	    = f->pv(it,p_cp,{t0:t1},{n0:n1})	; Potential vorticity 
    vort    = f->rvor(it,p_cp,{t0:t1},{n0:n1})  ; Relative vorticity (/s)
    z       = f->ht(it,p_cp,{t0:t1},{n0:n1})    ; Geopotential height (m)
    w	    = f->dz_dt(it,p_cp,{t0:t1},{n0:n1}) ; Vertical velocity (m/s)
    t	    = f->temp(it,p_cp,{t0:t1},{n0:n1})	; Temperature (K)
    q	    = f->q(it,p_cp,{t0:t1},{n0:n1})	; Specific humidity (kg/kg)

  ; Calculate theta-e using NCL function
    the   = pot_temp_equiv(pres,t,q,0,"q")

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
    avo@long_name   = "Vorticity (absolute)"

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

    ; Storm-relative velocity 
      u_rel	 = u_plane - u_cyc(ct)
      v_rel      = v_plane - v_cyc(ct)

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

     ilev = clev	; Work on both 900 and 850 hPa

     print("Working on: " + pres(ilev) + " hPa")
     p = pres(ilev)

     u_plane       = u(ilev,:,:)
     v_plane       = v(ilev,:,:)
     vort_plane	   = vort(ilev,:,:)
     geo_plane     = z(ilev,:,:)
     pv_plane	   = pv(ilev,:,:)
     w_plane	   = w(ilev,:,:)
     the_plane	   = the(ilev,:,:)
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

   ; LAT/LON grid
     lat_plane  = ind(t0.le.lat.and.lat.le.t1) ; Indices of latitude array
     lon_plane  = ind(n0.le.lon.and.lon.le.n1) ; Indices of longitude array

     lat_size = dimsizes(lat_plane)		  	  	; Array sizes
     lon_size = dimsizes(lon_plane)

     lat_sub = lat(lat_plane(0):lat_plane(lat_size-1)) ; Create subset of points 
     lon_sub = lon(lon_plane(0):lon_plane(lon_size-1)) ; between values given above

   ; LAT1/LON1 grid
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

   ; New arrays containing variables we want to plot
     u_int   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Zonal velocity
     v_int   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Meridional velocity
     spd_int = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vector velocity (RMW)
     pv_int  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Potential vorticity
     w_int   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vertical velocity
     th_int  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Equiv. potential temp.
     v_rad   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Radial velocity
     v_tan   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Tangential velocity

   ; Add variable metadata
     copy_VarAtts(u_plane,u_int)
     u_int@description = "Zonal velocity"
     u_int@units       = "m s~S~-1~N~"

     copy_VarAtts(v_plane,v_int)
     v_int@description = "Meridional velocity"
     v_int@units       = "m s~S~-1~N~"

     copy_VarAtts(w_plane,w_int)
     v_int@description = "Vertical velocity"
     v_int@units       = "m s~S~-1~N~"

     copy_VarAtts(the_plane,th_int)
     v_int@description = "Equivalent potential temperature"
     v_int@units       = "K"

     copy_VarAtts(pv_plane,pv_int)
     v_int@description = "Potential vorticity"
     v_int@units       = "10~S~-6~N~K kg~S~-1~N~ m~S~2~N~ s~S~-1~N~"

     v_rad@description = "Radial velocity"
     v_rad@units       = "m s~S~-1~N~"

     v_tan@description = "Tangential velocity"
     v_tan@units       = "m s~S~-1~N~"

     xpos!0    = "rad"
     xpos!1    = "azi"
     xpos&rad  = radii
     xpos&azi  = thetas

     ypos!0    = "rad"
     ypos!1    = "azi"
     ypos&rad  = radii
     ypos&azi  = thetas

     u_int!0   = "rad"
     u_int!1   = "azi"
     u_int&rad = radii
     u_int&azi = thetas

     v_int!0   = "rad"
     v_int!1   = "azi"
     v_int&rad = radii
     v_int&azi = thetas

     spd_int!0   = "rad"
     spd_int!1   = "azi"
     spd_int&rad = radii
     spd_int&azi = thetas

     w_int!0   = "rad"
     w_int!1   = "azi"
     w_int&rad = radii
     w_int&azi = thetas

     th_int!0   = "rad"
     th_int!1   = "azi"
     th_int&rad = radii
     th_int&azi = thetas

     pv_int!0   = "rad"
     pv_int!1   = "azi"
     pv_int&rad = radii
     pv_int&azi = thetas

     v_tan!0   = "rad"
     v_tan!1   = "azi"
     v_tan&rad = radii
     v_tan&azi = thetas

     v_rad!0   = "rad"
     v_rad!1   = "azi"
     v_rad&rad = radii
     v_rad&azi = thetas

   ; Latitude/longitude (small grid)
     do irad = 0, dimsizes(radii)-1
       r = radii(irad)
       do iang = 0, dimsizes(thetas)-1
         theta = thetas_rad(iang)
           xpos(irad,iang) = centre(ct,1) + r*cos(theta)
           ypos(irad,iang) = centre(ct,0) + r*sin(theta)
       end do
     end do

  ; Loop over azimuth angles to fill entire circular array of data points
    do iang = 0, dimsizes(thetas_rad)-1
      u_int(:,iang)    = linint2_points(lon_sub,lat_sub,u_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      v_int(:,iang)    = linint2_points(lon_sub,lat_sub,v_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      w_int(:,iang)    = linint2_points(lon1_sub,lat1_sub,w_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      th_int(:,iang)   = linint2_points(lon1_sub,lat1_sub,the_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      pv_int(:,iang)   = linint2_points(lon_sub,lat_sub,pv_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
    end do

    spd_int = sqrt((u_int ^ 2) + (v_int ^ 2))
    copy_VarCoords(u_int,spd_int)
    copy_VarMeta(u_int,spd_int)
    spd_int@units = "m s~S~-1~N~"

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
        vtan_azi(irad) = avg(v_tan(irad,:))      ; Tangential wind
	vrad_azi(irad) = avg(v_rad(irad,:))	 ; Radial wind 
	spd_azi(irad)  = avg(spd_int(irad,:))	 ; Vector velocity (for RMW)
	pvrt_azi(irad) = avg(pv_int(irad,:))	 ; PV 
	the_azi(irad)  = avg(th_int(irad,:))	 ; Theta-e
	vvel_azi(irad) = avg(w_int(irad,:))	 ; Vertical velocity
      end do
 
    ;=============================
    ; Calculate location of RMW
    ;=============================

      spd_max = max(spd_azi)
      dims_s  = dimsizes(spd_azi)
      spd_1d  = ndtooned(spd_azi)
      inds_s  = ind_resolve(maxind(spd_1d),dims_s)
      rmw     = radii(inds_s(0,0))

      spd_plot(ct)      = rmw			 ; Radius of maximum wind       

      vtan_plot0(:,ct)  = vtan_azi(:)            ; Tangential wind
      vrad_plot0(:,ct)  = vrad_azi(:)	    	 ; Radial wind
      pvrt_plot0(:,ct)  = pvrt_azi(:)	    	 ; PV
      the_plot0(:,ct)	= the_azi(:)		 ; Theta-e
      vvel_plot0(:,ct)  = vvel_azi(:)	    	 ; Vertical velocity

    ; Tidy up
      delete([/u_int,v_int,w_int,pv_int,v_tan,v_rad/])
      delete([/u_plane,v_plane,w_plane,the_plane,pv_plane,vort_plane,geo_plane/])

    ct  = ct + 1		    ; Counter variable (time)
    ct0 = ct0 + 1		    ; Independent counter variable (time)

  ; Tidy up before next iteration
    delete([/year,month,day,hour,minute,second,utc_date/])
    delete([/lon,lat,lon1,lat1,time,times,u,v,z,t,w,vort,pres,plevs/])

   end do     ; End time loop (do it = 0, times-1)

   delete([/tm0,tm1/])

 end do     ; End input file loop (do nf = 0, numINPUT-1)

;=====================================================
; Calculate radial and tangential wind tendencies
;=====================================================

; Initialise new arrays 
  pvrt_tend0 = pvrt_plot0
 
  pvrt_tend0!0    = "rad"
  pvrt_tend0&rad  = radii
  pvrt_tend0!1    = "time"
  pvrt_tend0&time = time_int

;;;;;;;;;;;;;;;;;;; EDIT FROM HERE ;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;; MAY NEED TO ADD EXTRA DIM TO LOOP ;;;;;;;;;;

  do rd = 0, rad_size

   rp1 = min((/rd+1,rad_size/)) 
   rm1	= max((/rd-1,0/))

   pvrt_tend0(rd,:) = (pvrt_plot0(rp1,:) - pvrt_plot0(rm1,:))/2 ; Units [m s-1 h-1] 

  end do    ; End radius loop (do rd = 0, rad_size)

;======================================
; Reorder dimensions before plotting 
;======================================

 ; Arrays now ordered [TIME * RADIUS] rather than [RADIUS * TIME]
   pvrt_tend = pvrt_tend0(time|:,rad|:)

   vtan_plot = vtan_plot0(time|:,rad|:)
   vrad_plot = vrad_plot0(time|:,rad|:)
   vvel_plot = vvel_plot0(time|:,rad|:)
   the_plot  = the_plot0(time|:,rad|:)
   pvrt_plot = pvrt_plot0(time|:,rad|:)

;=====================================================
; Create a panel plot of Figs. 6a and 6b from N11 
;=====================================================

  ; Output file location and type
    output = "$sam/nepartak/images/cyl_coords/n11_fig7_"+ens0
    wks = gsn_open_wks(opt,output)

  ; Set up panel plot
    panel = new(3,graphic)
    gsn_define_colormap(wks,"n11")

;==========================
; Options for plotting 
;==========================

  ; Potential vorticity
    opts_pvrt                              = True
    opts_pvrt@cnFillOn                     = True
    opts_pvrt@cnLineLabelInterval          = 2.0
    opts_pvrt@cnLineLabelFontHeightF       = 0.012
    opts_pvrt@cnLineLabelBackgroundColor   = "transparent"
    opts_pvrt@cnLineLabelPlacementMode     = "constant"
    opts_pvrt@cnLinesOn                    = False ; Contour lines off
    opts_pvrt@cnInfoLabelOn                = False
    opts_pvrt@cnLevelSelectionMode         = "ExplicitLevels"
    opts_pvrt@cnLevels                     = (/0.0, -0.5, -1.0, -1.5, -2.0, \
                                               -3.0, -4.0, -5.0, -6.0, -8.0, -10.0/)
    opts_pvrt@cnFillColors		   = (/13,12,11,10,9,8,7,6,5,4,3,0/)
    opts_pvrt@gsnPaperOrientation          = "landscape"
    opts_pvrt@tiMainString                 = ""
    opts_pvrt@tiMainFontHeightF            = 0.0125
    opts_pvrt@gsnLeftString                = ""
    opts_pvrt@gsnRightString               = ""
    opts_pvrt@gsnMaximize		   = True
    opts_pvrt@lbLabelBarOn		   = True		   ; Labelbar on/off
    opts_pvrt@lbBoxEndCapStyle		   = "TriangleBothEnds"	   ; Labelbar end shape
    opts_pvrt@pmLabelBarWidthF		   = 0.32		   ; Labelbar width
    opts_pvrt@pmLabelBarHeightF		   = 0.08		   ; Labelbar height
    opts_pvrt@lbLabelFontHeightF	   = 0.0125		   ; Labelbar font size
    opts_pvrt@lbLabelFont		   = "Helvetica"	   ; Labelbar font
    opts_pvrt@lbPerimOn			   = False		   ; Perimeter on/off
    opts_pvrt@gsnDraw			   = False		   ; Do not draw plot
    opts_pvrt@gsnFrame			   = False		   ; Do not advance frame

    opts_pvrt@vpWidthF            	   = 0.30		   ; Match width in N11
    opts_pvrt@vpHeightF                    = 0.50		   ; Match height in N11

  ; Additional plotting resources
    opts_pvrt@tiYAxisString                = "Time (hour)"
    opts_pvrt@tiXAxisString		   = "Radius (km)"
    opts_pvrt@tiXAxisFontHeightF	   = 0.018
    opts_pvrt@trYLog                       = False
    opts_pvrt@trYMinF			   = 30	     ; Starting point (T+...) for plot
    opts_pvrt@trYMaxF			   = 72	     ; Ending point (T+...) for plot

    opts_pvrt@trXMinF                      = 0.0                       ; Set x-axis min
    opts_pvrt@trXMaxF                      = 1.0                       ; Set x-axis max
    opts_pvrt@tmXTOn			   = "False" ; Turn off top x-axis TM
    opts_pvrt@tmYROn                       = "False" ; Turn off right y-axis TM
    opts_pvrt@gsnMaximize                  = True    ; Maximise plot size
    opts_pvrt@gsnAddCyclic                 = False
    opts_pvrt@tmXBMode			   = "Explicit"	; Set tick marks explicitly
    opts_pvrt@tmXBValues		   = (/0,0.5,1.0,1.5/)	       ; Tm positions
    opts_pvrt@tmXBLabels                   = (/"0","50","100","150"/)  ; Tm values
    opts_pvrt@tmXBMinorOn		   = False		       ; No minor tm
    opts_pvrt@tmYLMode			   = "Manual"		       
    opts_pvrt@tmYLTickStartF		   = opts_pvrt@trYMinF	       ; See L730
    opts_pvrt@tmYLTickEndF		   = opts_pvrt@trYMaxF	       ; See L731
    opts_pvrt@tmYLTickSpacingF		   = 6			       ; Tick spacing
    opts_pvrt@tmYLMinorOn		   = False		       ; Minor TMs off

    opts_pvrt@tiXAxisOffsetYF              = 0.10        ; Move x-axis title up/dn (h/l)
    opts_pvrt@pmLabelBarOrthogonalPosF     = 0.10        ; Move lb up/down (higher/lower)

  ; Resources to overlay contours of mean tangential wind
    opts_cont 	 	 	           = True    
    opts_cont@cnFillOn			   = False	 ; Filled contours on/off
    opts_cont@cnLineColor		   = "blue3"	 ; Contour line colour
    opts_cont@cnLevelSelectionMode         = "ExplicitLevels"
    opts_cont@cnLevels			   = (/10,20,30,40,50,60,70/)
    opts_cont@cnLineThicknessF		   = 2.0	 ; Contour line thickness
    opts_cont@cnInfoLabelOn                = False	 ; Info label on/off
    opts_cont@gsnDraw                      = False       ; Do not draw the plot
    opts_cont@gsnFrame                     = False       ; Do not advance the frame
    opts_cont@cnLineLabelsOn               = True        ; Turn line labels on/off
    opts_cont@cnLineLabelFontColor	   = "blue3"	 ; Line label colour
    opts_cont@cnLineLabelFontHeightF	   = 0.01	 ; Default 0.006 for vpWidth=0.30
    opts_cont@cnLineLabelFontThicknessF	   = 2.0	 ; Line label thickness
    opts_cont@tiMainString                 = ""
    opts_cont@gsnLeftString                = ""
    opts_cont@gsnRightString               = ""

  ; Theta-e resources
    opts_th				   = opts_pvrt
    delete(opts_th@cnLevels)
    delete(opts_th@cnFillColors)
    opts_th@cnLevels			   = (/351.0, 354.0, 357.0, 360.0, 363.0, \
                                               366.0, 369.0, 372.0, 376.0, 380.0/)
    opts_th@cnFillColors                   = (/0,2,3,4,5,6,7,8,9,10,11/)

  ; Vertical velocity resources
    opts_vvel				   = opts_pvrt
    delete(opts_vvel@cnLevels)
    delete(opts_vvel@cnFillColors)
    opts_vvel@cnLevels			   = (/0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, \
    					       0.7, 0.8, 0.9, 1.0/)
    opts_vvel@cnFillColors		   = (/0,2,3,4,5,6,7,8,9,10,11,12/)

  ; Resources to overlay zero line
    opts_zero                             = True
    opts_zero@cnFillOn                    = False
    opts_zero@cnLineColor                 = "gray69"
    opts_zero@cnLevels                    = (/-1.0, 0.0, 1.0/)
    opts_zero@cnLineThicknessF            = 2.0
    opts_zero@cnInfoLabelOn               = False
    opts_zero@gsnDraw                     = False       ; Do not draw the plot
    opts_zero@gsnFrame                    = False       ; Do no advance the frame
    opts_zero@gsnContourZeroLineThicknessF = 2.0
    opts_zero@cnLineLabelsOn               = False      ; Turn off line labels
    opts_zero@tiMainString                = ""
    opts_zero@gsnLeftString               = ""
    opts_zero@gsnRightString              = ""

  ; Radius of maximum wind (RMW) resources
    opts_rmw                               = True
    opts_rmw@xyLineColors                  = "black"
    opts_rmw@xyLineThicknesses             = 3.0
    opts_rmw@xyDashPatterns                = 0.0
    opts_rmw@gsnDraw                       = False       ; Do not draw the plot
    opts_rmw@gsnFrame                      = False       ; Do not advance the frame
    opts_rmw@tiMainString                  = ""
    opts_rmw@gsnLeftString                 = ""
    opts_rmw@gsnRightString                = ""

;==============
; Panel plot 
;==============

    pvrt0 = gsn_csm_contour(wks,pvrt_tend,opts_pvrt)	; Potential vorticity
    cont0 = gsn_csm_contour(wks,vtan_plot,opts_cont)	; Tangential wind 
    rmw0  = gsn_csm_xy(wks,spd_plot,time_int,opts_rmw)  ; RMW	     
    overlay(pvrt0,cont0)                                ; Overlay tangential wind
    overlay(pvrt0,rmw0)					; Overlay RMW 
    panel(0) = pvrt0					; Plot first panel

    the0  = gsn_csm_contour(wks,the_plot,opts_th)	; Theta-e
    cont0 = gsn_csm_contour(wks,vtan_plot,opts_cont)	; Tangential wind
    rmw0  = gsn_csm_xy(wks,spd_plot,time_int,opts_rmw)  ; RMW
    overlay(the0,cont0)    				; Overlay tangential wind
    overlay(pvrt0,rmw0)                                 ; Overlay RMW
    panel(1) = the0					; Plot second panel

    vvel0 = gsn_csm_contour(wks,vvel_plot,opts_vvel)	; Vertical velocity
    zero0 = gsn_csm_contour(wks,vvel_plot,opts_zero)	; Zero line (vertical velocity)
    zero  = ColorNegDashZeroPosContour(zero0,"transparent","grey69","transparent")
    rmw0  = gsn_csm_xy(wks,spd_plot,time_int,opts_rmw)  ; RMW
    overlay(vvel0,zero)					; Overlay vvel zero line
    overlay(pvrt0,rmw0)                                 ; Overlay RMW
    panel(2) = vvel0					; Plot third panel
    delete([/zero0,zero/])				; Tidy up 

    optsP                       = True
    optsP@gsnFrame              = False    ; Do not advance the frame
    optsP@gsnPanelLabelBar      = False    ; Turn on panel labelbar
    optsP@txString              = ""       ; title_arr(0,it)
    optsP@gsnPanelFigureStrings = ""       ; (/"a) "+ens0, "b) "+ens1/)
    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.98
    optsP@gsnPanelBottom        = 0.02
    optsP@amJust                = "TopLeft"
    optsP@gsnPanelFigureStringsFontHeightF = 0.0099 ; Reduce label size (default 0.01)

    gsn_panel(wks,panel,(/1,3/),optsP)             ; Draw as a single plot
    frame(wks)

end 