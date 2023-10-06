; Script to calculate azimuthally averaged vertical velocity following a TC 
; [ncl int=1 dist=1 'opt="x11"' 'ens0="em02"' 'ens1="em08"' tc_pv_ens_csec_paper.ncl]

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$LIB/gsn_csm.ncl"

begin

;=========================================================
; Calculate storm motion [for system-relative velocity]
;=========================================================

  ; Prelash (define arrays for loops below)
    ens_arr = (/ens0,ens1/)            ; Ensemble member string array
    centre  = new((/2,18,78,2/),float) ; TC centre at each time interval/pressure level
    x_cyc   = new((/2,78/),float)      ; Cyclone longitude
    y_cyc   = new((/2,78/),float)      ; Cyclone latitude
    u_cyc   = new((/2,78/),float)      ; Zonal cyclone velocity
    v_cyc   = new((/2,78/),float)      ; Meridional cyclone velocity
    vel_cyc = new((/2,78/),float)      ; Cyclone velocity (vector)

  ; 00Z 4 Jul [00], 06Z 4 Jul [01], 12Z 4 Jul [02], 18Z 4 Jul [03], 00Z 5 Jul [04],
  ; 06Z 5 Jul [05], 12Z 5 Jul [06], 18Z 5 Jul [07], 00Z 6 Jul [08], 06Z 6 Jul [09],
  ; 12Z 6 Jul [10], 18Z 6 Jul [11], 00Z 7 Jul [12], 06Z 7 Jul [13], 12Z 7 Jul [14],
  ; 18Z 7 Jul [15], 00Z 8 Jul [16], 06Z 8 Jul [17], 12Z 8 Jul [18], 18Z 8 Jul [19],
  ; 00Z 9 Jul [20]

 do en = 0, 1

    print_clock("Working on ensemble member: "+ens_arr(en))  

  ; List all files to read in and analyse
  ; Analyse from 12Z 5 -- 6 July [78 h]
    diri        = "$sam/um/cp/ens/20160704T0000Z_ra1t_"
    fili_pb     = diri+ens_arr(en)+"_pb_"
    input_pb    = (/fili_pb+"04.nc",fili_pb+"05.nc",fili_pb+"06.nc",fili_pb+"07.nc",\
                    fili_pb+"08.nc",fili_pb+"09.nc",fili_pb+"10.nc",\
                    fili_pb+"11.nc",fili_pb+"12.nc",fili_pb+"13.nc",\
                    fili_pb+"14.nc",fili_pb+"15.nc",fili_pb+"16.nc"/)
    numINPUT    = dimsizes(input_pb)        ; Number of input files
    numTIMES    = numINPUT * 6              ; Total number of times in all files
    ct          = 0                         ; Counter variable

    setvalues NhlGetWorkspaceObjectId
      "wsMaximumSize" : 1000000000
    end setvalues

;==========================
; Loop over input files
;==========================

  do nf = 0, numINPUT-1     ; Loop over all times in file [00Z - 12Z 5 July]
  f = addfile(input_pb(nf),"r")

;==================================
; Get the variables we will need
;==================================

  time  = f->t            ; Times in file (6 per file - hourly)

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

  z     = z/10
  z@units = "dam"         ; Convert to decametres
  vort  = vort * (10 ^ -6) ; Convert to /s

  dx    = lon(1) - lon(0) ; Grid spacing (longitude)
  dy    = lat(1) - lat(0) ; Grid spacing (latitude)

  ; Pressure levels (pres)
  ; (0) 1000, (1) 950, (2) 900, (3) 850, (4) 800, (5) 750, (6) 700, (7) 650
  ; (8) 600, (9) 550, (10) 500, (11) 450, (12) 400, (13) 350
  ; (14) 300, (15) 250, (16) 200, (17) 150, (18) 100

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

     print_clock("Working on time: "+time_str(it))

;================================================================================
; Find cyclone centre on subset of global grid (use both height and vorticity)
;================================================================================

   ; Create array subsets to speed up calculations
     lat_0 = 6.92   ; 5
     lat_1 = 31.02  ; 26.02
     lon_0 = 109.48 ; 114
     lon_1 = 147.02 ; 152.02

   ; Choose a level in the lower troposphere (950 hPa)
     ilev = 1

   ; First, concentrate on a subset of the global grid
     u_plane    = u(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
     v_plane    = v(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
     geo_plane  = z(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
     vort_plane = vort(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
     dsize      = dimsizes(u_plane)

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
     print("Cyclone centre (geo): "+lat_min950+" degrees N, "+lon_min950+" degrees E")

   ;========================================================
   ; EMPLOY A SAFETY NET IF WE IDENTIFY THE WRONG CYCLONE
   ;========================================================
     if (ct.eq.0) then
      centre(en,ilev,ct,0) = lat_min950
      centre(en,ilev,ct,1) = lon_min950

     else
      lt0   = centre(en,ilev,ct-1,0)         ; Previous TC latitude
      ln0   = centre(en,ilev,ct-1,1)         ; Previous TC longitude
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
       print("New centre (geo): "+lat_min950+" degrees N, "+lon_min950+" degrees E")
       centre(en,ilev,ct,0) = lat_min950
       centre(en,ilev,ct,1) = lon_min950

       delete([/geo_plane0, geo_min0, dims_h0, geo_1d0, inds_h0/])

      else
       print("Cyclone centre (geo): "+lat_min950+" degrees N, "+lon_min950+" degrees E")
       centre(en,ilev,ct,0) = lat_min950
       centre(en,ilev,ct,1) = lon_min950
      end if

     end if

    ; Tidy up before looping over pressure levels
      delete([/vort1d, inds, vort1d_sm, vort_smth, inds_sm/])
      delete([/u_plane, v_plane, vort_plane, geo_plane/])

      ct = ct + 1

    end do    ; End of time loop (do it = 0, times-1)

  ; Tidy up to avoid array dimension errors [but leave alone on final iteration]
     if (en.eq.0) then
      delete([/time,times,z,u,v,vort/])
      delete([/utc_date,year,month,day,hour,\
               minute,second,date_str,time_str,out_str/])
     else
      if (nf.lt.numINPUT-1)
        delete([/time,times,z,u,v,vort/])
        delete([/utc_date,year,month,day,hour,\
                minute,second,date_str,time_str,out_str/])
      end if
     end if

   end do     ; End of input file loop (do nf = 0, numINPUT-1)

  end do      ; End of ensemble loop (do en = 0, 1)

;============================================================================
; Calculate distance between two points --> then calculate system velocity
;============================================================================

   print_clock("Starting system velocity calculations...")

 do en = 0, 1

 ; Loop to calculate system velocity from hourly cyclone position [lat,lon]
   do iang = 1, numTIMES-1
    x_cyc(en,iang)   = gc_latlon(centre(en,ilev,iang-1,0),centre(en,ilev,iang-1,1),\
                       centre(en,ilev,iang-1,0),centre(en,ilev,iang,1),2,3)
    y_cyc(en,iang)   = gc_latlon(centre(en,ilev,iang-1,0),centre(en,ilev,iang-1,1),\
                       centre(en,ilev,iang,0),centre(en,ilev,iang-1,1),2,3)
    u_cyc(en,iang)   = x_cyc(en,iang)/(3600)
    v_cyc(en,iang)   = y_cyc(en,iang)/(3600)
    vel_cyc(en,iang) = sqrt(u_cyc(en,iang)^2.0 + v_cyc(en,iang)^2.0)
   end do

 ; Set initial system velocity [0] equal to that at the second time [1]
   u_cyc(en,0)   = u_cyc(en,1)
   v_cyc(en,0)   = v_cyc(en,1)
   vel_cyc(en,0) = vel_cyc(en,1)

 ; Add metadata to new arrays
   copy_VarAtts(u,u_cyc(en,:))
   u_cyc@units         = "m s~S~-1~N~"
   u_cyc@description   = "Cyclone zonal velocity"

   copy_VarAtts(v,v_cyc(en,:))
   v_cyc@units         = "m s~S~-1~N~"
   v_cyc@description   = "Cyclone meridional velocity"

   vel_cyc@units       = "m s~S~-1~N~"
   vel_cyc@description = "Cyclone speed"

   u_cyc(en,:) = -(u_cyc(en,:)) ; Easterly motion is negative

  end do      ; End of ensemble loop (do en = 0, 1)

 ; Tidy up
   delete([/lat,lon,lat1,lon1,pres,plevs,u,v,z,vort/])
   delete([/time,times,utc_date,year,month,day,hour,\
            minute,second,date_str,time_str,out_str/])

;==================================================================================
;;;;;;;;;;;;;;;;;;;;;;;;;;; Start of main part of script ;;;;;;;;;;;;;;;;;;;;;;;;;;
;==================================================================================

    print_clock("Starting main part of script...")

   ; Prelash
     time_arr    = new((/2,numTIMES/),string)
     title_arr   = new((/2,numTIMES/),string)

   ;========================================
   ; Define cylindrical coordinate arrays
   ;========================================

   ; Loop over pressure levels (1000 - 100 hPa, every 50 hPa)
     p_start = 0   ; 1000 hPa
     p_end   = 17  ; 150 hPa
     p_int   = 1

     count_p = 0                   ; Counter for loop
     num_p   = (p_end - p_start)+1 ; Vertical array size

     thetas = new(73,float)

   ; Define radial distances [r1 = outer distance; r2 = spacing] 
   ; r1 = 0.50, r2 = 11 draws radial circles every 0.05 degrees
   ; r1 = 0.50, r2 = 5 draws radial circles every 0.125 degrees
   ; r1 = 3.00, r2 = 25 draws radial circles every 0.125 degrees

     pii        = 3.14159265
     radii      = fspan(0.0,0.75,16)    ; Size of radial circle array
     thetas     = fspan(0.0,360.0,73)   ; Array of azimuth angles (0, 5, 10, etc)
     thetas_rad = (pii/180)*thetas      ; Degrees to radians
     rad_size   = dimsizes(radii)-1     ; Index for later calculation

   ; Define arrays for later use in cross-section plots
     pv_azi = new((/2,num_p,dimsizes(radii)/),float) ; 0   - 360 (Whole cyclone)

   ; Define 4-dimensional array to hold all values from 'pv_azi' at ALL times
     pv_plot = new((/num_p,dimsizes(radii),numTIMES,2/),float)
     pv_plot!0 = "lev"
     pv_plot!1 = "rad"
     pv_plot!2 = "time"
     pv_plot!3 = "ens"

   ; Add metadata for azimuthally averaged potential vorticity
     pv_azi!0     = "ens"
     pv_azi!1     = "pres"
     pv_azi!2     = "rad"
     pv_azi&rad   = radii
     pv_azi@description = "Azimuthally averaged potential vorticity (Whole Cyclone)"
     pv_azi@units = "PVU"

;=========================
; Start multiple loops
;=========================

 do en = 0, 1

    print_clock("Working on ensemble member: "+ens_arr(en))

  ; List all files to read in and analyse
  ; Analyse from 12Z 5 -- 6 July [78 h]
    diri        = "$sam/um/cp/ens/20160704T0000Z_ra1t_"
    fili_pb     = diri+ens_arr(en)+"_pb_"
    input_pb    = (/fili_pb+"04.nc",fili_pb+"05.nc",fili_pb+"06.nc",fili_pb+"07.nc",\
                    fili_pb+"08.nc",fili_pb+"09.nc",fili_pb+"10.nc",\
                    fili_pb+"11.nc",fili_pb+"12.nc",fili_pb+"13.nc",\
                    fili_pb+"14.nc",fili_pb+"15.nc",fili_pb+"16.nc"/)
    numINPUT    = dimsizes(input_pb)        ; Number of input files
    numTIMES    = numINPUT * 6              ; Total number of times in all files
    ct          = 0                         ; Counter variable

    setvalues NhlGetWorkspaceObjectId
      "wsMaximumSize" : 1000000000
    end setvalues

;===============================
; Start loop over input files 
;===============================

    do nf = 0,  numINPUT-1     ; 12Z 5 -- 6 July
    f = addfile(input_pb(nf),"r")

;==================================
; Get the variables we will need
;==================================

    time  = f->t            ; Times in file (6 per file - hourly)
    times = dimsizes(time)

    lon1  = f->longitude_1  ; longitude (1098 points --> 109.04 to 152.92 degrees E)
    lat1  = f->latitude_1   ; latitude  (810 points --> 1.8 to 34.16 degrees N)
    lon   = f->longitude    ; longitude (1098 points --> 109.02 to 152.90 degrees E)
    lat   = f->latitude     ; latitude (811 points --> 1.78 to 34.18 degrees N)

    pres  = f->p 	    ; Pressure (levels)
    plevs = dimsizes(pres)  ; Size of pressure level array
    pres@units = "hPa"

    u     = f->u            ; Zonal velocity (m/s) [time | pres | lat | lon]
    v     = f->v            ; Meridional velocity (m/s) [time | pres | lat | lon]
    vort  = f->rvor         ; Relative vorticity (/s) [time | pres | lat | lon]

    t     = f->temp	    ; Temperature (K) [time | pres | lat1 | lon1]
    z     = f->ht   	    ; Geopotential height (m) [time | pres | lat1 | lon1]
    w     = f->dz_dt	    ; Vertical velocity (m/s) [time | pres | lat1 | lon1]

  ; Read in potential vorticity (PV)
    pv    = f->pv           ; Read in PV
    pv    = pv * (10 ^ 6)   ; Convert to PVU [10^6 K m2 s-1 kg-1]

  ; Convert to units we want
    z     = z/10
    z@units = "dam" ; Convert to decametres 

    vort  = vort * (10 ^ -6) ; Convert to /s

    dx    = lon(1) - lon(0) ; Grid spacing (longitude)
    dy    = lat(1) - lat(0) ; Grid spacing (latitude)

  ; Pressure levels (pres)
  ; (0) 1000, (1) 950, (2) 900, (3) 850, (4) 800, (5) 750, (6) 700, (7) 650
  ; (8) 600, (9) 550, (10) 500, (11) 450, (12) 400, (13) 350
  ; (14) 300, (15) 250, (16) 200, (17) 150, (18) 100

 ; Degrees to radians (for calculation below)
   pii = 3.14159265
   lat_rad = (pii/180)*lat
   lon_rad = (pii/180)*lon

 ; Calculate absolute vorticity
   d  = dimsizes(vort)
   f0 = new((/d(0),d(1),d(2),d(3)/),float)  ; Earth's vorticity

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

 ; Absolute vorticity
   avo = vort + f0
   copy_VarMeta(vort,avo)
   avo@long_name = "Vorticity (planetary)"
   avo@title     = "Absolute vorticity"
   avo@name      = "avo"
   avo@standard_name = "atmosphere_planetary_vorticity"
   avo@units     = "s~S~-1~N~"

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

     time_arr(en,ct) = out_str(it)

     print_clock("Working on time: "+time_str(it))
     title_arr(en,ct) = "Valid at: "+time_str(it)

;================================================================================
; Find cyclone centre on subset of global grid (use both height and vorticity)
;================================================================================

    ; Choose a level in the lower troposphere (950 hPa)
      ilev = 1

    ; Create array subsets to speed up calculations
      lat_0 = 6.92
      lat_1 = 31.02
      lon_0 = 109.48
      lon_1 = 147.02

    ; First, concentrate on a subset of the global grid 
      u_plane    = u(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
      v_plane    = v(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
      geo_plane  = z(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
      vort_plane = vort(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
      avo_plane	 = avo(it,ilev,{lat_0:lat_1},{lon_0:lon_1})

    ; Storm-relative velocity 
      u_rel	 = u_plane - u_cyc(en,ct)
      v_rel      = v_plane - v_cyc(en,ct)

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
      print("Cyclone centre (geo): "+lat_min950+" degrees N, "+lon_min950+" degrees E")

   ;========================================================
   ; EMPLOY A SAFETY NET IF WE IDENTIFY THE WRONG CYCLONE
   ;========================================================
     if (ct.eq.0) then
      centre(en,ilev,ct,0) = lat_min950
      centre(en,ilev,ct,1) = lon_min950

     else
      lt0   = centre(en,ilev,ct-1,0)         ; Previous TC latitude
      ln0   = centre(en,ilev,ct-1,1)         ; Previous TC longitude
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
       print("New centre (geo): "+lat_min950+" degrees N, "+lon_min950+" degrees E")
       centre(en,ilev,ct,0) = lat_min950
       centre(en,ilev,ct,1) = lon_min950

       delete([/geo_plane0, geo_min0, dims_h0, geo_1d0, inds_h0/])

      else
       print("Cyclone centre (geo): "+lat_min950+" degrees N, "+lon_min950+" degrees E")
       centre(en,ilev,ct,0) = lat_min950
       centre(en,ilev,ct,1) = lon_min950
      end if

     end if

    ; Tidy up before looping over pressure levels
      delete([/vort1d, inds, vort1d_sm, vort_smth, inds_sm/])
      delete([/u_plane, v_plane, vort_plane, geo_plane/])

;===============================
; Loop over pressure levels 
;===============================

    count_p = 0                   ; Counter variable (pressure)	 

    do ilev = 0, p_end, p_int

     print("Working on: " + pres(ilev) + " hPa")
     p = pres(ilev)

     u_plane       = u(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
     v_plane       = v(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
     w_plane       = w(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
     vort_plane	   = vort(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
     pv_plane	   = pv(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
     avo_plane     = avo(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
     geo_plane     = z(it,ilev,{lat_0:lat_1},{lon_0:lon_1})

     r		   = 1.5

     geo_plane1	   = z(it,ilev,{lat_min950-r:lat_min950+r},\
     		       {lon_min950-r:lon_min950+r})
     vort_plane1   = vort(it,ilev,{lat_min950-r:lat_min950+r},\
     		          {lon_min950-r:lon_min950+r})

   ;==============================================
   ; Find cyclone centre on each pressure level
   ;==============================================

   ; GEOPOTENTIAL HEIGHT 
     geo_min1       = toint(min(geo_plane1)) ; Minimum geopotential on pressure level
     geo_max1       = toint(max(geo_plane1)) ; Maximum geopotential on pressure level
     geo_smth1      = smth9_Wrap(geo_plane1, 0.5, 0.5, True)
     geo_min_smth1  = min(geo_smth1)     

     dims_h1        = dimsizes(geo_plane1)
     geo1d_1        = ndtooned(geo_plane1)
     inds_h1        = ind_resolve(minind(geo1d_1),dims_h1)
     geo1d_sm1      = ndtooned(geo_smth1)
     inds_h_sm1     = ind_resolve(minind(geo1d_sm1),dims_h1)     

     lat_min	    = lat_min950-r + (dy * inds_h1(0,0))
     lon_min 	    = lon_min950-r + (dy * inds_h1(0,1))
     print("Cyclone centre (geo): "+lat_min+" degrees N, "+lon_min+" degrees E")

   ; RELATIVE VORTICITY
     vort_max1      = max(vort_plane1)
     vort_smth1     = smth9_Wrap(vort_plane1, 0.5, 0.5, True)
     vort_max_smth1 = max(vort_smth1)

     dims1          = dimsizes(vort_plane1)
     vort1d_1       = ndtooned(vort_plane1)
     inds1          = ind_resolve(maxind(vort1d_1),dims1)
     vort1d_sm1     = ndtooned(vort_smth1)
     inds_sm1       = ind_resolve(maxind(vort1d_sm1),dims1)

     lat_max        = lat_max950-r + (dy * inds1(0,0))
     lon_max        = lon_max950-r + (dy * inds1(0,1))
     print("Cyclone centre (vort): "+lat_max+" degrees N, "+lon_max+" degrees E")

     lat_plane = ind(lat_0.le.lat.and.lat.le.lat_1) ; Indices of latitude array
     lon_plane = ind(lon_0.le.lon.and.lon.le.lon_1) ; Indices of longitude array

   ; Modify latitude/longitude arrays
     lat_size = dimsizes(lat_plane)		  	  	; Array sizes
     lon_size = dimsizes(lon_plane)

     lat_sub = lat(lat_plane(0):lat_plane(lat_size-1)) ; Create subset of points 
     lon_sub = lon(lon_plane(0):lon_plane(lon_size-1)) ; between values given above

     delete([/vort1d_1,inds1,vort1d_sm1,vort_smth1,inds_sm1/])
     delete([/geo1d_1,inds_h1,geo1d_sm1,geo_smth1,inds_h_sm1/])

   ;===========================================================================
   ; If centre of cyclone on pressure level is too far from 950 hPa position,
   ; then use the 950 hPa position (calculated using geopotential height)
   ;===========================================================================

   ; Calculate difference between centres
     lat_diff = abs(lat_max-lat_min950)
     lon_diff = abs(lon_max-lon_min950)

   ; TC centre position on each pressure level in loop (1000,950,..)
     centre(en,ilev,ct,0) = lat_min
     centre(en,ilev,ct,1) = lon_min

;===========================================
; Switch to cylindrical coordinate system
;===========================================

   ; Create new arrays describing cylindrical coordinate system (2D for now)
     xpos = new((/dimsizes(radii),dimsizes(thetas)/),float)
     ypos = new((/dimsizes(radii),dimsizes(thetas)/),float)

   ; New arrays containing variables we want to plot
     u_int   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Zonal velocity
     v_int   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Meridional velocity
     pv_int  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Potential vorticity

   ; Add variable metadata
     copy_VarAtts(u_plane,u_int)
     u_int@description = "Zonal velocity"
     u_int@units       = "m s~S~-1~N~"

     copy_VarAtts(v_plane,v_int)
     v_int@description = "Meridional velocity"
     v_int@units       = "m s~S~-1~N~"

     copy_VarAtts(w_plane,pv_int)
     pv_int@description = "Potential vorticity"
     pv_int@units       = "PVU"

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

     pv_int!0   = "rad"
     pv_int!1   = "azi"
     pv_int&rad = radii
     pv_int&azi = thetas

   ; Latitude/longitude (small grid)
     do irad = 0, dimsizes(radii)-1
       r = radii(irad)
       do iang = 0, dimsizes(thetas)-1
         theta = thetas_rad(iang)
           xpos(irad,iang) = centre(en,ilev,ct,1) + r*cos(theta)
           ypos(irad,iang) = centre(en,ilev,ct,0) + r*sin(theta)
       end do
     end do

  ; Loop over azimuth angles to fill entire circular array of data points
    do iang = 0, dimsizes(thetas_rad)-1
      u_int(:,iang)    = linint2_points(lon_sub,lat_sub,u_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      v_int(:,iang)    = linint2_points(lon_sub,lat_sub,v_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      pv_int(:,iang)   = linint2_points(lon_sub,lat_sub,pv_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
    end do

    ;========================================================
    ; Average potential vorticity over storm quadrants
    ;========================================================

      do irad = 0, dimsizes(radii)-1
        r = radii(irad)
        pv_azi(en,count_p,irad) = avg(pv_int(irad,:))             ; Whole Cyclone
      end do

      printMinMax(pv_azi(en,count_p,:),False)		  ; Print variable information
      pv_plot(count_p,:,ct,en) = pv_azi(en,count_p,:)     ; Input into larger array

      count_p = count_p + 1  ; Counter for loop

    ; Tidy up
      delete([/u_int,v_int,pv_int/])
      delete([/w_plane,u_plane,v_plane/])
      delete([/avo_plane,vort_plane,geo_plane,pv_plane/])
      delete([/geo_plane1,vort_plane1/])

    end do                  ; End pressure level loop

    ct = ct + 1		    ; Counter variable (time)

   end do     ; End time loop (do it = 0, 1)

 ; Add correct metadata to array before plotting
   pv_plot!0   = "lev"
   pv_plot&lev = pres(p_start:p_end)

 ; Tidy up before next iteration
   delete([/year,month,day,hour,minute,second,utc_date,date_str,time_str,out_str/])
   delete([/time,times,u,v,z,t,w,vort,pres,plevs/])

 end do     ; End input file loop (do nf = 0, numINPUT-1)

end do	    ; End ensemble member loop (do en = 0, 1)

;=====================================================
; Create a panel plot of multiple ensemble members
;=====================================================

  do it = 0, numTIMES-1

  ; Output file location and type
    output = \
   "$sam/nepartak/images/cyl_coords/panel_pv_csec_"+ens0+"_"+ens1+"_"+time_arr(0,it)
    wks = gsn_open_wks(opt,output)

  ; Set up panel plot
    panel = new(2,graphic)
    gsn_define_colormap(wks,"vort_panel")

;==========================
; Options for plotting 
;==========================

  ; Potential vorticity
    opts_em00                             = True
    opts_em00@cnFillOn                    = True
    opts_em00@cnLineLabelInterval         = 2.0
    opts_em00@cnLineLabelFontHeightF      = 0.012
    opts_em00@cnLineLabelBackgroundColor  = "transparent"
    opts_em00@cnLineLabelPlacementMode    = "constant"
    opts_em00@cnLinesOn                   = False ; Contour lines off
    opts_em00@cnInfoLabelOn               = False
    opts_em00@cnLevelSelectionMode        = "ExplicitLevels"

    opts_em00@cnLevels                    = (/4.0, 5.0, 6.0, 8.0, \
                                              10.0, 12.0, 14.0, 16.0, 18.0, \
                                              20.0, 25.0, 30.0, 35.0/)
    opts_em00@cnFillColors                = (/0,13,15,17,19,\
                                              20,21,22,23,24,\
                                              25,26,27,28/)

    opts_em00@gsnPaperOrientation         = "landscape"
    opts_em00@tiMainString                = ""
    opts_em00@tiMainFontHeightF           = 0.0125
    opts_em00@gsnLeftString               = ""
    opts_em00@gsnRightString              = ""
    opts_em00@gsnMaximize		  = True
    opts_em00@lbLabelBarOn		  = False
    opts_em00@gsnDraw			  = False
    opts_em00@gsnFrame			  = False

  ; Additional plotting resources
    opts_em00@tiYAxisString               = "Pressure (hPa)"
    opts_em00@tiXAxisString		  = "Radial distance from cyclone centre (km)"
    opts_em00@tiXAxisFontHeightF	  = 0.018
    opts_em00@trYLog                      = False

    opts_em00@tmXTOn			  = "False" ; Turn off top x-axis TM
    opts_em00@tmYROn                      = "False" ; Turn off right y-axis TM
    opts_em00@gsnMaximize                 = True    ; Maximise plot size
    opts_em00@gsnAddCyclic                = False
    opts_em00@tmXBMode			  = "Explicit"	; Set tick marks explicitly
    opts_em00@tmXBValues                  = (/0,1,2,3/) ; Tick mark positions
    opts_em00@tmXBLabels                  = (/"0","25","50","75"/)   ; Tick mark values
    opts_em00@tmXBMinorValues             = radii		     ; Minor tick marks

;==============
; Panel plot 
;==============

    panel(0) = gsn_csm_pres_hgt(wks,pv_plot(:,:,it,0),opts_em00)
    panel(1) = gsn_csm_pres_hgt(wks,pv_plot(:,:,it,1),opts_em00)

    optsP                       = True
    optsP@gsnFrame              = False    ; Do not advance the frame
    optsP@gsnPanelLabelBar      = True     ; Turn on panel labelbar
    optsP@pmLabelBarWidthF      = 0.6
    optsP@pmLabelBarHeightF     = 0.15
    optsP@lbLabelFontHeightF    = 0.01
    optsP@lbPerimOn             = False
    optsP@lbLabelFont           = "Helvetica"
    optsP@lbTitleString = "Azimuthally averaged potential vorticity (PVU)"
    optsP@lbTitleFontHeightF    = 0.015
    optsP@lbTitleDirection      = "Across"
    optsP@lbTitlePosition       = "Bottom"
    optsP@txString              = title_arr(0,it)
    optsP@gsnPanelFigureStrings = (/"a) "+ens0, "b) "+ens1/)
    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.98
    optsP@gsnPanelBottom        = 0.02
    optsP@amJust                = "TopLeft"
    optsP@gsnPanelFigureStringsFontHeightF = 0.0125 ; Label size (default 0.01)

    gsn_panel(wks,panel,(/1,2/),optsP)              ; Draw as a single plot
    frame(wks)

  end do     ; End time loop (do it = 0, numTIMES-1)

end 