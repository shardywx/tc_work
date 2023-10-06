; Script to calculate tangential and radial velocity following a TC 
; and to plot radius-height cross sections of both variables in a panel plot 

; Run using: 

; ncl dat=\"02T12\" opt=\"pdf\" ens0=\"em07\" dist=3.0 calc=1 cn0=\"slp\" tc_vboth_xz.ncl

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'opt'   = output file format ("pdf" or "x11")
; 'ens0'  = ensemble member (em00 -- em11)
; 'dist'  = size of box following storm (degrees)
; 'calc'  = method of finding storm centre
; 'cn0'   = offline storm centre: "slp", "vort", "geo", "geo_sm"

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/nc_levs.ncl"

begin

;==============================
; Prelash (reading in files)
;==============================

; Find 'pb' stream data to read in and plot
  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0
  fili_b1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pb*.nc")
  fili_b      = diri+"/"+fili_b1

;=============================================================
; Calculate number of times in file using built-in function
;=============================================================

  numINFO_b   = nc_times(fili_b)
  numINPUT_b  = numINFO_b[0]            ; Number of times
  time_b      = numINFO_b[1]            ; Array of times

;===============================================
; Calculate number of vertical levels in file
;===============================================

  nLEVS_b     = nc_levs(fili_b)
  numLEVS_b   = nLEVS_b[1]              ; Number of vertical levels
  levs_b      = nLEVS_b[0]              ; Array of vertical levels
  levs_b@units = "hPa"

;============================================================
; Calculate storm motion using built-in function ('st_rm')
;============================================================

; 'dat'       = initialisation time (02T12, 03T00, ...)
; 'ens0'      = ensemble simulation (em00, em01, ...)
; 'diri'      = path to input files (see above)
; 'dist'      = size of box following storm (degrees)
; 'mins'      = analyse 1-h (0) or 5-min (1) data

  storm_rel   = st_rm(dat, ens0, "$ar/text/", dist, 0)
  u_cyc       = storm_rel[0]             ; Zonal wind
  v_cyc       = storm_rel[1]             ; Meridional wind
  vel_cyc     = storm_rel[2]             ; Vector wind
  numINPUT_a  = storm_rel[3]             ; Number of times in file
  lat_arr     = storm_rel[4]             ; Latitude array subset (following storm)
  lon_arr     = storm_rel[5]             ; Longitude array subset (following storm)
  centre      = storm_rel[6]             ; Storm track (position) array

; Create arrays for data and output information
  llbox       = toint(dist*50)                   ; Calculate domain size (grid points)
  dsize       = (/llbox,llbox/)                  ; Domain size (grid points)
  d0          = dsize(0)                         ; Size of dimension 0
  d1          = dsize(1)                         ; Size of dimension 1

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 10000000000
  end setvalues

;=======================================
; Choose which storm track to read in
;=======================================

  if (cn0 .eq. "slp") then
   cn = 0
  elseif (cn0 .eq. "vort") then
   cn = 1
  elseif (cn0 .eq. "geo") then
   cn = 2
  elseif (cn0 .eq. "geo_sm") then
   cn = 3
  end if

;==========================================================================
; Create date/time string arrays for all times in file (use 'pb' stream)
;==========================================================================

  month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
                    "Jul","Aug","Sep","Oct","Nov","Dec"/)

  utc_date = cd_calendar(time_b, 0)

  year   = tointeger(utc_date(:,0))
  month  = tointeger(utc_date(:,1))
  day    = tointeger(utc_date(:,2))
  hour   = tointeger(utc_date(:,3))
  minute = tointeger(utc_date(:,4))
  second = utc_date(:,5)

; Correct for errors in the code (round up value of hour when minutes = 59)
  do ct = 0, numINPUT_b-1
   if (minute(ct).gt.30) then
     hour(ct) = hour(ct)+1
   end if
  end do
  delete(ct)

; Create arrays to hold strings for output plots
  date_str  = new(numINPUT_b,string)
  time_str  = new(numINPUT_b,string)
  time_arr  = new(numINPUT_b,string)
  title_arr = new(numINPUT_b,string)

  date_str = sprinti("%0.2iUTC ", hour) + \
             sprinti("%0.2i ", day) \
           + month_abbr(month)
  time_str = sprinti("%0.2i UTC ", hour) + \
             sprinti("%0.2i ", day) \
           + month_abbr(month)
  time_arr = sprinti("%0.2i", day) + \
             month_abbr(month) + \
             "_" + sprinti("%0.2iZ", hour)

;==================================================================================
;;;;;;;;;;;;;;;;;;;;;;;;;;; Start of main part of script ;;;;;;;;;;;;;;;;;;;;;;;;;;
;==================================================================================

  print_clock("Starting main part of script...")

;========================================
; Define cylindrical coordinate arrays
;========================================

  thetas     = new(73,float)
  pii        = 3.14159265
  radii      = fspan(0.0,1.5,31)     ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,73)   ; Array of azimuth angles (0, 5, 10, etc)
  thetas_rad = (pii/180)*thetas      ; Degrees to radians
  rad_size   = dimsizes(radii)-1     ; Index for later calculation

;================
; Radial wind 
;================

; Define 3D array to hold all values from 'vrad_av4' at ALL times
  vrad_plot 	       = new((/numLEVS_b,dimsizes(radii),numINPUT_a/),float)
  vrad_plot!0	       = "lev"
  vrad_plot&lev	       = levs_b
  vrad_plot!1	       = "rad"
  vrad_plot!2 	       = "time"

; Define array to hold azimuthally averaged radial wind
  vrad_av4 	       = new((/numLEVS_b,dimsizes(radii)/),float)
  vrad_av4!0           = "lev" 				  ; Change to 'lev' ???
  vrad_av4&lev	       = levs_b
  vrad_av4!1           = "rad"
  vrad_av4&rad         = radii
  vrad_av4@description = "Azimuthally averaged radial wind"
  vrad_av4@units       = "m s~S~-1~N~"

;===================
; Tangential wind 
;===================

; Define 3D array to hold all values from 'vtan_av4' at ALL times
  vtan_plot            = new((/numLEVS_b,dimsizes(radii),numINPUT_a/),float)
  vtan_plot!0          = "lev"
  vtan_plot&lev	       = levs_b
  vtan_plot!1          = "rad"
  vtan_plot!2          = "time"

; Define array to hold azimuthally averaged tangential wind
  vtan_av4             = new((/numLEVS_b,dimsizes(radii)/),float) ; 0 - 360 (Whole cyclone)
  vtan_av4!0           = "lev"
  vtan_av4&lev	       = levs_b
  vtan_av4!1           = "rad"
  vtan_av4&rad         = radii
  vtan_av4@description = "Azimuthally averaged tangential wind"
  vtan_av4@units       = "m s~S~-1~N~"

;=========================
; Start multiple loops
;=========================

  b = addfile(fili_b,"r")

;==================================
; Get the variables we will need
;==================================

  do it = ts0, tf0                              ; Loop over desired times

 ; MAYBE EDIT THESE LINES (NOT SURE)
   print_clock("Working on time: "+time_str(it))
   itt = it + 1
   title_arr(it) = "Valid at "+time_str(it)+" (T+"+itt+")"

   time  = b->t(it)        ; Time in file
   times = dimsizes(time)

   t0    = lat_arr(cn,it,0)
   t1    = lat_arr(cn,it,1)
   n0    = lon_arr(cn,it,0)
   n1    = lon_arr(cn,it,1)
   print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

 ; Create latitude/longitude array subsets to reduce computation time
   lon1  = b->longitude_1({n0:n1})  ; longitude ['d0' grid points]
   lat1  = b->latitude_1({t0:t1})   ; latitude  ['d1' grid points]
   lon   = b->longitude({n0:n1})    ; longitude ['d0' grid points]
   lat   = b->latitude({t0:t1})     ; latitude  ['d1' grid points]

 ; Pressure levels (pres)
 ; (0) 1000, (1) 950, (2) 925, (3) 900, (4) 850, (5) 800, (6) 750, (7) 700, (8) 650,
 ; (9) 600, (10) 500, (11) 400, (12) 300, (13) 250, (14) 200, (15) 150, (16) 100

   u0    = b->u(it-3,:,{t0:t1},{n0:n1})     ; Zonal wind (m/s)
   v0    = b->v(it-3,:,{t0:t1},{n0:n1})     ; Meridional wind (m/s)
   vort  = b->rvor(it-3,:,{t0:t1},{n0:n1})  ; Relative vorticity (/s)

   t     = b->temp(it-3,:,{t0:t1},{n0:n1})  ; Temperature (K)
   z     = b->ht(it-3,:,{t0:t1},{n0:n1})    ; Geopotential height (m)
   w     = b->dz_dt(it-3,:,{t0:t1},{n0:n1}) ; Vertical wind (m/s)

   pres  = b->p                             ; Pressure (hPa)
   plevs = dimsizes(pres)
   pres@units = "hPa"

 ; Convert to units we want
   z     = z/10
   z@units = "dam" ; Convert to decametres
   vort  = vort * (10 ^ -6) ; Convert to /s

   dx0   = lon(1) - lon(0) ; Grid spacing (longitude)
   dy0   = lat(1) - lat(0) ; Grid spacing (latitude)

 ; Coriolis parameter and absolute vorticity
   f1    = coriolis_param(lat)
   f0    = conform_dims(dimsizes(vort),f1,1)
   copy_VarMeta(vort,f0)
   copy_VarCoords(vort,f0)
   f0@description  = "Coriolis parameter"
   f0@name         = "Coriolis parameter"

   avo   = vort + f0
   copy_VarMeta(vort,avo)
   copy_VarAtts(vort,avo)
   copy_VarCoords(vort,avo)
   avo@description = "Absolute vorticity"
   avo@long_name   = "Vorticity (absolute)"

;==============================================================
; Make sure all arrays are same size (for later calculations)
;==============================================================

 ; 2-D array (6 * 3)
   size  = (/dimsizes(u0),dimsizes(v0),dimsizes(vort),dimsizes(t),dimsizes(z),dimsizes(w)/)

 ; Reduce size of any dimension larger than the rest (u,v,w,t,etc)
   do sz  = 0, 2
    nind = min( size(:,sz) )
    if ( any (size(:,sz) .gt. nind) ) then
     size(:,sz) = nind
    end if
   end do

 ; Resize arrays based upon results of calculations above
   u2    = u0(0:size(0,0)-1,0:size(0,1)-1,0:size(0,2)-1)
   v2    = v0(0:size(1,0)-1,0:size(1,1)-1,0:size(1,2)-1)
   vort2 = vort(0:size(2,0)-1,0:size(2,1)-1,0:size(2,2)-1)
   t2    = t(0:size(3,0)-1,0:size(3,1)-1,0:size(3,2)-1)
   z2    = z(0:size(4,0)-1,0:size(4,1)-1,0:size(4,2)-1)
   w2    = w(0:size(5,0)-1,0:size(5,1)-1,0:size(5,2)-1)
   lat2  = lat(0:nind-1)
   lon2  = lon(0:nind-1)
   delete([/u0,v0,vort,t,z,w,lat,lon/])

 ; Reset names to avoid further code changes
   u0    = u2
   v0    = v2
   vort  = vort2
   t     = t2
   z     = z2
   w     = w2
   lat   = lat2
   lon   = lon2
   delete([/u2,v2,vort2,t2,z2,w2,lat2,lon2/])

;===================================
; Calculate storm-relative winds
;===================================

   u = u0 - u_cyc(cn,it)
   v = v0 - v_cyc(cn,it)

 ; Add metadata from original horizontal wind arrays
   copy_VarCoords(u0,u)
   copy_VarAtts(u0,u)

   copy_VarCoords(v0,v)
   copy_VarAtts(v0,v)

;==========================================================
; Calculate divergence using one-sided difference scheme
;==========================================================

   div = uv2dv_cfd(u, v, lat, lon, 2)
   copy_VarCoords(u, div)
   div@units = "s~S~-1~N~"
   div@name = "Divergence on pressure levels"

;===============================
; Loop over pressure levels 
;===============================

   do ilev = 0, plevs-1

    p = pres(ilev)
    print_clock("Working on: " + p + " hPa")

  ; Read in storm-relative wind on pressure levels
    u_plane       = u(ilev,0:d0-1,0:d1-1)
    v_plane       = v(ilev,0:d0-1,0:d1-1)

  ; Read in other variables on pressure levels
    vort_plane    = vort(ilev,0:d0-1,0:d1-1)
    avo_plane     = avo(ilev,0:d0-1,0:d1-1)
    div_plane     = div(ilev,0:d0-1,0:d1-1)
    geo_plane     = z(ilev,0:d0-1,0:d1-1)
    t_plane       = t(ilev,0:d0-1,0:d1-1)

  ; Calculate minimum and maximum geopotential height on pressure level
    zmin          = toint(min(geo_plane))
    zmax          = toint(max(geo_plane))

  ; Locate relative vorticity maximum (storm centre) on pressure level
    vort_max      = max(vort_plane)
    vort_smth     = smth9_Wrap(vort_plane, 0.5, 0.5, True)
    vort_max_smth = max(vort_smth)

  ; Find index of relative vorticity maximum
    dims          = dimsizes(vort_plane)
    vort1d        = ndtooned(vort_plane)
    inds          = ind_resolve(maxind(vort1d),dims)
    vort1d_sm     = ndtooned(vort_smth)
    inds_sm       = ind_resolve(maxind(vort1d_sm),dims)

   ;==================================================================
   ; Either calculate the centre, or read in pre-calculated values
   ;==================================================================

   if (calc .eq. 0) then

   ; Find latitude and longitude of storm centre using vorticity calculation above
     lat_max0        = t0 + (dy0 * inds(0,0))
     lon_max0        = n0 + (dx0 * inds(0,1))

     print("Centre (vort): "+lat_max0+" degrees N, "+lon_max0+" degrees E")

   ; Create smaller grid around location identified using vorticity maximum
     lt1 = lat_max0 - 0.15
     lt2 = lat_max0 + 0.15
     ln1 = lon_max0 - 0.15
     ln2 = lon_max0 + 0.15

   ; Calculate minimum geopotential height in region immediately surrounding 'storm centre'
     geo_plane0 = geo_plane({lt1:lt2},{ln1:ln2})       ; Smaller grid (1.0 x 1.0 deg)
     geo_min0   = min(geo_plane0)                      ; Find minimum on smaller grid
     dims_z0    = dimsizes(geo_plane0)                 ; Size of smaller grid
     geo_1d0    = ndtooned(geo_plane0)                 ; Create 1-D array
     inds_z0    = ind_resolve(minind(geo_1d0),dims_z0) ; Find index of minimum

     lat_max0   = lt1 + (dy0 * inds_z0(0,0))
     lon_max0   = ln1 + (dx0 * inds_z0(0,1))

     print("Centre (vort + geo): "+lat_max0+" degrees N, "+lon_max0+" degrees E")
     delete([/geo_plane0,geo_min0,dims_z0,geo_1d0,inds_z0/])

   ; Distance between centre calculated above using relative vorticity, and that calculated offline
     dlat       = abs(lat_max0 - centre(cn,it,0))
     dlon       = abs(lon_max0 - centre(cn,it,1))

   ; If our storm centre calculations deviate too much (> 0.75º) from the values
   ; we read in using the NCL function 'st_rm', calculate a new centre using
   ; the minimum geopotential height. Constrain the new calculation using the values
   ; from 'st_rm', to prevent spurious centres being identified and analysed.

     if (dlat .gt. 1.00 .or. dlon .gt. 1.00) then

      lt1 = centre(cn,it,0) - 0.25
      lt2 = centre(cn,it,0) + 0.25
      ln1 = centre(cn,it,1) - 0.25
      ln2 = centre(cn,it,1) + 0.25

      geo_plane0 = geo_plane({lt1:lt2},{ln1:ln2})       ; Smaller grid (1.0 x 1.0 deg)
      geo_min0   = min(geo_plane0)                      ; Find minimum on smaller grid
      dims_z0    = dimsizes(geo_plane0)                 ; Size of smaller grid
      geo_1d0    = ndtooned(geo_plane0)                 ; Create 1-D array
      inds_z0    = ind_resolve(minind(geo_1d0),dims_z0) ; Find index of minimum

      lat_max0   = lt1 + (dy0 * inds_z0(0,0))
      lon_max0   = ln1 + (dx0 * inds_z0(0,1))
      print("New centre: "+lat_max0+" degrees N, "+lon_max0+" degrees E")
      delete([/geo_plane0,geo_min0,dims_z0,geo_1d0,inds_z0/])

     end if

   elseif (calc .eq. 1) then

      lat_max0   = centre(cn,it,0)
      lon_max0   = centre(cn,it,1)

   end if

   ;==============================================================
   ; Create arrays for later switch to cylindrical coordinates
   ;==============================================================

   ; Indices of variables on 'lat'/'lon' grid for 'pb' stream
     lat_plane  = ind(t0.le.lat.and.lat.le.t1) ; Indices of latitude array
     lon_plane  = ind(n0.le.lon.and.lon.le.n1) ; Indices of longitude array

   ; Array sizes
     lat_size = dimsizes(lat_plane)
     lon_size = dimsizes(lon_plane)

   ; Array of latitude/longitude points corresponding to the indices above
     lat_sub = lat(lat_plane(0):lat_plane(lat_size-1))
     lon_sub = lon(lon_plane(0):lon_plane(lon_size-1))

   ; Indices of variables on 'lat1'/'lon1' grid for 'pb' stream
     lat1_plane  = ind(t0.le.lat1.and.lat1.le.t1) ; Indices of latitude array
     lon1_plane  = ind(n0.le.lon1.and.lon1.le.n1) ; Indices of longitude array

   ; Array sizes
     lat1_size = dimsizes(lat1_plane)
     lon1_size = dimsizes(lon1_plane)

   ; Array of latitude/longitude points corresponding to the indices above
     lat1_sub = lat1(lat1_plane(0):lat1_plane(lat1_size-1))
     lon1_sub = lon1(lon1_plane(0):lon1_plane(lon1_size-1))

     delete([/vort1d,inds,vort1d_sm,vort_smth,inds_sm/])

;===========================================
; Switch to cylindrical coordinate system
;===========================================

 ; Create new arrays describing cylindrical coordinate system (2D for now)
   xpos = new((/dimsizes(radii),dimsizes(thetas)/),float)
   ypos = new((/dimsizes(radii),dimsizes(thetas)/),float)

 ; New arrays containing variables we want to plot
   u_int   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Zonal wind
   v_int   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Meridional wind
   v_rad   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Radial wind
   v_tan   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Tangential wind

 ; Add variable metadata
   copy_VarAtts(u_plane,u_int)
   u_int@description = "Zonal wind"
   u_int@units       = "m s~S~-1~N~"

   copy_VarAtts(v_plane,v_int)
   v_int@description = "Meridional wind"
   v_int@units       = "m s~S~-1~N~"

   v_rad@description = "Radial wind"
   v_rad@units       = "m s~S~-1~N~"

   v_tan@description = "Tangential wind"
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
         xpos(irad,iang) = centre(cn,it,1) + r*cos(theta)
         ypos(irad,iang) = centre(cn,it,0) + r*sin(theta)
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

;==================================================================
; Average the radial & tangential wind over storm quadrants
;==================================================================

   do irad = 0, dimsizes(radii)-1
    r = radii(irad)

     vrad_av4(ilev,irad) = avg(v_rad(irad,:))             ; Whole Cyclone
     vtan_av4(ilev,irad) = avg(v_tan(irad,:))             ; Whole Cyclone

   end do

 ; Input into larger array
   vrad_plot(ilev,:,it) = vrad_av4(ilev,:)
   vtan_plot(ilev,:,it) = vtan_av4(ilev,:)

 ; Tidy up
   delete([/u_int,v_int,v_tan,v_rad,u_plane,v_plane,vort_plane,geo_plane/])

  end do                  ; End pressure level loop

; Tidy up before next iteration
  delete([/lon1,lat1,lon,lat,u0,v0,u,v,vort,t,z,w,f0,f1,avo,div/])

 end do     ; End time loop (do it = ts0, tf0)

;=====================================================
; Create a panel plot of multiple ensemble members
;=====================================================

 do it = ts0, tf0

  ; Output file location and type
    output = "$sam/nepartak/images/cyl_coords/vboth_xz_"+dat+"_"+ens0+"_"+time_arr(it)
    wks    = gsn_open_wks(opt,output)

  ; Set up panel plot
    panel  = new(2,graphic)
    gsn_define_colormap(wks,"BlueDarkRed18")

;==========================
; Options for plotting 
;==========================

  ; Radial wind
    opts_rad                             = True
    opts_rad@cnFillOn                    = True
    opts_rad@cnLineLabelInterval         = 2.0
    opts_rad@cnLineLabelFontHeightF      = 0.012
    opts_rad@cnLineLabelBackgroundColor  = "transparent"
    opts_rad@cnLineLabelPlacementMode    = "constant"
    opts_rad@cnLinesOn                   = False ; Contour lines off
    opts_rad@cnInfoLabelOn               = False
    opts_rad@cnLevelSelectionMode        = "ExplicitLevels"
    opts_rad@cnLevels                    = (/-10., -5., -3., -2., -1., -0.5, -0.2, \
                                              0, 0.2, 0.5, 1., 2., 3., 5., 10./)
    opts_rad@cnFillColors                = (/2,4,5,6,7,8,9,0,0,12,13,14,15,16,17,19/)
    opts_rad@gsnPaperOrientation         = "landscape"
    opts_rad@tiMainString                = ""
    opts_rad@tiMainFontHeightF           = 0.0125
    opts_rad@gsnLeftString               = ""
    opts_rad@gsnRightString              = ""
    opts_rad@gsnMaximize		 = True
    opts_rad@lbLabelBarOn		 = True
    opts_rad@lbBoxEndCapStyle            = "TriangleBothEnds"    ; Labelbar end shape
    opts_rad@lbLabelFontHeightF          = 0.0125                ; Labelbar font size
    opts_rad@lbLabelFont                 = "Helvetica"           ; Labelbar font
    opts_rad@lbPerimOn                   = False                 ; Perimeter on/off
    opts_rad@gsnDraw			 = False
    opts_rad@gsnFrame			 = False

  ; Additional plotting resources
    opts_rad@tiYAxisString                = "Pressure (hPa)"
    opts_rad@tiXAxisString		  = "Radial distance from cyclone centre (km)"
    opts_rad@tiXAxisFontHeightF		  = 0.020
    opts_rad@trYLog                       = False

    opts_rad@pmLabelBarOrthogonalPosF     = 0.10    ; Move lb up/down (higher/lower)
    opts_rad@tiXAxisOffsetYF              = 0.12    ; Move x-axis title up/down (h/l)

    opts_rad@tmXTOn			  = "False" ; Turn off top x-axis TM
    opts_rad@tmYROn                       = "False" ; Turn off right y-axis TM
    opts_rad@gsnMaximize                  = True    ; Maximise plot size
    opts_rad@gsnAddCyclic                 = False
    opts_rad@tmXBMode			  = "Explicit"	; Set tick marks explicitly
    opts_rad@tmXBValues                   = (/0,0.5,1.0,1.5/) ; Tick mark positions
    opts_rad@tmXBLabels                   = (/"0","50","100","150"/) ; Tick mark values
    opts_rad@tmXBMinorValues              = fspan(0.0, 3.0, 16)      ; Minor tick marks

  ; Resources to overlay zero line 
    opts_cont 	 	      	          = True
    opts_cont@cnFillOn			  = False
    opts_cont@cnLineColor		  = "Black"
    opts_cont@cnLevelSelectionMode	  = "ExplicitLevels"
    opts_cont@cnLevels			  = (/0.0/)
    opts_cont@cnLineThicknessF		  = 3.0
    opts_cont@cnInfoLabelOn		  = False
    opts_cont@gsnDraw			  = False	; Do not draw the plot
    opts_cont@gsnFrame			  = False	; Do no advance the frame
    opts_cont@gsnContourZeroLineThicknessF = 3.0
    opts_cont@cnLineLabelsOn		   = False	; Turn off line labels
    opts_cont@tiMainString                = ""
    opts_cont@gsnLeftString               = ""
    opts_cont@gsnRightString              = ""

    vrad0    = gsn_csm_pres_hgt(wks,vrad_plot(:,:,it),opts_rad)
    cont0a   = gsn_csm_contour(wks,vrad_plot(:,:,it),opts_cont)
    cont0    = ColorNegDashZeroPosContour(cont0a,"transparent","black","transparent")
    overlay(vrad0,cont0)
    panel(0) = vrad0

  ;=================================
  ; Now plot tangential wind
  ;=================================

    opts_tan                             = opts_rad
    delete([/opts_tan@cnLevels, opts_tan@cnFillColors/])

    opts_tan@cnLineLabelInterval         = 2.0
    opts_tan@cnLineLabelFontHeightF      = 0.012
    opts_tan@cnLineLabelBackgroundColor  = "transparent"
    opts_tan@cnLineLabelPlacementMode    = "constant"
    opts_tan@cnLinesOn                   = False ; Contour lines off
    opts_tan@cnInfoLabelOn               = False
    opts_tan@cnFillPalette               = "radar_new"
    opts_tan@cnLevelSelectionMode        = "ExplicitLevels"
    opts_tan@cnLevels                    = (/5.0, 10.0, 15.0, 20.0,\
                                             25.0, 30.0, 35.0, 40.0, \
                                             45.0, 50.0, 55.0, 60.0, 65.0/)
    opts_tan@cnFillColors 		 = (/0,2,3,4,5,6,7,8,9,10,11,12,13,14/)
    opts_tan@gsnPaperOrientation         = "landscape"
    opts_tan@tiMainString                = ""
    opts_tan@tiMainFontHeightF           = 0.0125
    opts_tan@gsnLeftString               = ""
    opts_tan@gsnRightString              = ""
    opts_tan@gsnMaximize                 = True
    opts_tan@gsnDraw                     = False
    opts_tan@gsnFrame                    = False

    panel(1) = gsn_csm_pres_hgt(wks,vtan_plot(:,:,it),opts_tan)

;==============
; Panel plot 
;==============

    optsP                       = True
    optsP@gsnFrame              = False    ; Do not advance the frame
    optsP@gsnPanelLabelBar      = False    ; Turn panel labelbar on/off
    optsP@txString              = ""       ; title_arr(0,it)
    optsP@gsnPanelFigureStrings = ""       ; (/"a) 'Strong'", "b) 'Weak'"/)
    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.98
    optsP@gsnPanelBottom        = 0.51
    optsP@amJust                = "TopLeft"
    optsP@gsnPanelFigureStringsFontHeightF = 0.0098 ; Reduce label size (default 0.01)

    gsn_panel(wks,panel,(/1,2/),optsP)             ; Draw 2 panels on 1 row 
    frame(wks)

  end do     ; End time loop (do it = ts0, tf0)

end 