; Script to calculate circulation following a TC, and produce a time-height plot 

; Run using:

; ncl dat=\"02T12\" opt=\"pdf\" ens0=\"em07\" rad=0.76 dist=2.2 nbox=25 \
; ts0=24 tf0=48 calc=1 cn0=\"slp\" tc_circ_new.ncl

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'opt'   = output file format ("pdf" or "x11")
; 'ens0'  = ensemble member (em00 -- em11)
; 'dist'  = size of box following storm (degrees)
; 'ts0'   = starting time index (0-119)
; 'tf0'   = ending time index (0-119)
; 'nbox'  = number of ensemble box positions for circ. budget calculations (49,25,9,1)
; 'rad'   = radius of box around TC for budget calculations (integer multiple of dx,dy)
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

; Exit script if missing times are chosen
  if (ts0 .lt. 3) then
   print("Starting time index must be 3 or greater (data missing)! Please try again...")
   exit()
  end if

; Find 'pb' stream data (pressure level data)
  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0
  fili_b1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pb*.nc")
  fili_b      = diri+"/"+fili_b1

;=============================================================
; Calculate number of times in file using built-in function
;=============================================================

; 'pb' stream
  numINFO_b   = nc_times(fili_b)
  numINPUT_b  = numINFO_b[0]            ; Number of times
  time_b      = numINFO_b[1]            ; Array of times

;===============================================
; Calculate number of vertical levels in file
;===============================================

; 'pb' stream
  nLEVS_b     = nc_levs(fili_b)
  numLEVS_b   = nLEVS_b[1]              ; Number of vertical levels
  levs_b      = nLEVS_b[0]              ; Array of vertical levels

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

; Prelash [define arrays for use in loops below]

  pos_arr     = new((/numINPUT_a,2/),float) ; Array for cyclone position [lat/lon]
  r_earth     = 6.37e6                      ; Earth's radius (constant)
  plevs       = numLEVS_b                   ; Pressure levels (17)

; Define constants
  rd    = 287.0   ; Gas constant for dry air (J/kg/K^2)
  cp0   = 1004.0  ; Specific heat of dry air at constant pressure (J/kg/K)
  kp    = 0.286   ; For Exner function calculations (0.286)
  pref  = 1000.0  ; Reference pressure p0 (1000 hPa)
  pii	= 3.14159265

; Array to hold circulation tendency values [ens (nbox) * time (120) * levs (19)]
  circ_inst0 = new((/nbox,numINPUT_a,plevs/),float)
  circ_inst1 = new((/nbox,numINPUT_a,plevs/),float)

;===============================
; Start loop over input files 
;===============================

  b         = addfile(fili_b,"r")               ; Read in 'pb' stream [u,v,w,vort,t,z,rh]

;==============================================
; Get the variables we will need [pb stream]
;==============================================

  do it = ts0, tf0                              ; Loop over desired times

 ; MAYBE EDIT THESE LINES (NOT SURE)
   print("Working on time: "+time_str(it))
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

   u0    = b->u(it-3,:,{t0:t1},{n0:n1})     ; Zonal velocity (m/s)
   v0    = b->v(it-3,:,{t0:t1},{n0:n1})     ; Meridional velocity (m/s)
   vort  = b->rvor(it-3,:,{t0:t1},{n0:n1})  ; Relative vorticity (/s)

   t     = b->temp(it-3,:,{t0:t1},{n0:n1})  ; Temperature (K)
   z     = b->ht(it-3,:,{t0:t1},{n0:n1})    ; Geopotential height (m)
   w     = b->dz_dt(it-3,:,{t0:t1},{n0:n1}) ; Vertical velocity (m/s)

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



;==========================================================
; Define latitude and longitude for finite differencing
;==========================================================

 ; Get dimensions of vorticity array (???)
   d_vort     = dimsizes(vort)
   d_1        = d_vort(1)
   d_2        = d_vort(2)

 ; Create lat/lon grid with same dimensions as all input variables
   if (d_1 .eq. d0 .and. d_2 .eq. d1) then
    y2d      = new((/d0,d1/),float)
   else
    y2d      = vort(0,:,:)
   end if

   y2d!0     = "lat"
   y2d!1     = "lon"
   y2d&lat   =  lat
   y2d&lon   =  lon

   y2d       = conform_dims(dimsizes(y2d),lat,0) ; Array of latitude values
   coslat    = cos(y2d*pii/180)                  ; Same array in radians

 ; Determine horizontal resolution of data
   delta_lon = y2d
   delta_lat = y2d

   delta_lon = lon(1) - lon(0) ; Zonal grid spacing
   delta_lat = lat(1) - lat(0) ; Meridional grid spacing

   delta_lon = sqrt(delta_lon ^ 2)
   delta_lat = sqrt(delta_lat ^ 2)

   lons      = lon(1) - lon(0)
   lats      = lat(1) - lat(0)
   lons      = sqrt(lons^2)
   lats      = sqrt(lats^2)

 ; Calculate dx and dy (changes with longitude but not with latitude)
   ddx       = y2d
   ddx       = r_earth*(pii/180.0)*delta_lon*cos(y2d*pii/180)
   ddy       = y2d
   ddy       = r_earth*(pii/180)*delta_lat

 ; Give 'coslat' the dimensions of 'u' (3-D) and hold values in new array ('coslat3d')
   coslat3d  = conform_dims(dimsizes(u),coslat,(/1,2/))

 ; Get information about box size for later ensemble calculations
   r = rad   ; Box radius (degrees)
   box = ddx({centre(cn,it,0)-r:centre(cn,it,0)+r},{centre(cn,it,1)-r:centre(cn,it,1)+r})
   box_size = dimsizes(box)
   size0    = box_size(0)
   size1    = box_size(1)

   s0       = size0-1
   s1       = size1-1

;============================================================
; Loop over pressure levels to calculate vertical profile
;============================================================

    cp = 0               ; Counter variable

    do ilev = 0, plevs-1 ; 1000,950,900,850,800,750,700,650,600,550,...,100

     p = pres(ilev)
     print_clock("Working on: " + p + " hPa")

;=================================
; Pressure levels ['pb' stream]
;=================================

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

 ; Add lat/lon values to array for later calculations
   pos_arr(it,0) = lat_max0
   pos_arr(it,1) = lon_max0

;=====================================================================================
; Create array of box positions either side of the original [Davis & Galarneau 2009]
;=====================================================================================

   if (nbox .eq. 49) then
    lat_box = (/lat_max0-(3*dy0),lat_max0-(2*dy0),lat_max0-(dy0),lat_max0,\
                lat_max0+(dy0),lat_max0+(2*dy0),lat_max0+(3*dy0)/)
    lon_box = (/lon_max0-(3*dx0),lon_max0-(2*dx0),lon_max0-(dx0),lon_max0,\
                lon_max0+(dx0),lon_max0+(2*dx0),lon_max0+(3*dx0)/)
   elseif (nbox .eq. 25) then
    lat_box = (/lat_max0-(2*dy0),lat_max0-(dy0),lat_max0,\
                lat_max0+(dy0),lat_max0+(2*dy0)/)
    lon_box = (/lon_max0-(2*dx0),lon_max0-(dx0),lon_max0,\
                lon_max0+(dx0),lon_max0+(2*dx0)/)
   elseif (nbox .eq. 9) then
    lat_box = (/lat_max0-(dy0),lat_max0,lat_max0+(dy0)/)
    lon_box = (/lon_max0-(dx0),lon_max0,lon_max0+(dx0)/)
   elseif (nbox .eq. 1) then
    lat_box = (/lat_max0/)
    lon_box = (/lon_max0/)
   end if

   lsize = dimsizes(lat_box) * dimsizes(lon_box)
   
;====================================================
; Define arrays to hold ensemble of box positions
;====================================================

 ; EDIT --> NEATER WAY TO INITIALISE MULTIPLE ARRAYS
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

;    ; Change units 
;      ddx_box(a,:,:) = ddx_box(a,:,:) / 1000
;      ddy_box(a,:,:) = ddy_box(a,:,:) / 1000
;      u_box(a,:,:)   = u_box(a,:,:)   / 1000
;      v_box(a,:,:)   = v_box(a,:,:)   / 1000

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
      print("Area of box = "+d_ave(a))
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
     circ_instT0 = circ_s + circ_n + circ_w + circ_e
     circ_instT0 = circ_instT0 / d_ave(a)		     ; Divide by area of box 
     print("Circulation (wind) = "+circ_instT0)

   ; Enter values into predefined arrays
     circ_inst0(a,it,cp)     = circ_instT0   	             ; Circulation (winds)

   ; Tidy up (delete variables from within lat/lon loops before next iteration)
     delete([/circ_instT0/])
     delete([/vort_x_final0, vort_x_final1, vort_y_final0, vort_y_final1/])
     delete([/avo_x_final0,  avo_x_final1,  avo_y_final0,  avo_y_final1/])
     delete([/uy_final0,     uy_final1,     ux_final0,     ux_final1/])
     delete([/vy_final0,     vy_final1,     vx_final0,     vx_final1/])
     delete([/numPOINTS_x0, numPOINTS_x, numPOINTS_y/])

     end do     ; End latitude loop (do j = 0, dimsizes(lat_box)-1)
    end do	; End longitude loop (do i = 0, dimsizes(lon_box)-1)

    cp                   = cp + 1       ; Counter variable (pressure)

  ; Tidy up (delete variables from within pressure level loop)
    delete([/u_plane, v_plane/])
    delete([/vort_plane, avo_plane, div_plane, geo_plane, t_plane/])
    delete([/ddx_box, ddy_box, u_box, v_box/])
    delete([/vort_box, avo_box, div_box/])
    delete([/div_ave,vort_ave,vort_avg,avo_avg,u_avg,v_avg/])

   end do ; End pressure level loop (do ilev = 0, plevs-1)

 ; Tidy up (delete variables from within time loop before next iteration)
   delete([/lon,lat,lon1,lat1,u0,v0,u,v,vort,div,f0,f1,avo,t,z,w,box/])

  end do                ; End time loop (do it = ts0, tf0)

;=====================================================
; Calculate term averages over all ensemble members
;=====================================================

   print_clock("Finished ensemble calculations! Now averaging over all members...")
   circ_instA0 = dim_avg_n_Wrap(circ_inst0,0)
   circ_instA0 = circ_instA0 * (10 ^ 4)
   circ_tendA0 = circ_instA0

;=========================================
; Integrate over Nepartak's life cycle
;=========================================

 ; Each calculation produces an array of 19 values [pressure levels]
 ; We do this calculation [numINPUT_a] times
 ; Create array with [numINPUT_a] * 19 elements

   time_int = new(numINPUT_a,"integer")
   time_int = ispan(0,numINPUT_a-1,1)
   time_int!0 = "time"

 ; Add metadata before time-height plotting [shouldn't need lines below]

 ; Circulation using wind speed
   circ_instA0!0    = "time"
   circ_instA0&time = time_int
   circ_instA0!1    = "lev"
   circ_instA0&lev  = pres

 ; Circulation tendency 
   circ_tendA0!0    = "time"
   circ_tendA0&time = time_int
   circ_tendA0!1    = "lev"
   circ_tendA0&lev  = pres

   do it = 0, numINPUT_a-1 

    tp1  = min((/it+1,numINPUT_a-1/))
    tm1  = max((/it-1,0/))  

    circ_tendA0(it,:) = (circ_instA0(tp1,:) - circ_instA0(tm1,:)) / 7200 ; Units [s-2]

   end do 

   circ_tendA0	    = circ_tendA0 * (10 ^ 4)

 ; Reorder from [time * pressure] to [pressure * time]
   circ_tendF0	    = circ_tendA0(lev|:,time|:)
   circ_instF0      = circ_instA0(lev|:,time|:)

   circ_tendF1	    = runave_n_Wrap(circ_tendF0, 3, 0, 1)
   circ_tendF2	    = runave_n_Wrap(circ_tendF0, 6, 0, 1)

;====================================================
; Open workstation for plotting (time-height plots)
;====================================================

     output = "$sam/nepartak/images/cart_coords/circ_xz_"+dat+"_"+ens0+"_"+r1
     wks = gsn_open_wks(opt,output)

;========================
; Options for plotting 
;========================

     panel = new(1,graphic)
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

     if (r .eq. 0.76) then
     opts_circ@cnLevels                    = (/-0.5, 0.0, 0.5, 1.0, \
                                                1.5, 2.0, 3.0, 4.0, 5.0, 6.0, \
                                                7.0, 8.0, 9.0, 10.0, 11.0, 12.0, \
                                                13.0, 14.0/)
     opts_circ@cnFillColors                = (/6,8,0,12,13,14,15,\
                                               16,18,19,20,21,22,\
					       23,24,25,26,27,28/)
     elseif (r .eq. 0.48) then 
     opts_circ@cnLevels                    = (/-1.0, 0.0, 1.0, 2.0, \
                                                3.0, 4.0, 5.0, 6.0, 8.0, 10.0, \
                                                12.0, 14.0, 16.0, 18.0, 20.0, 22.0, \
                                                24.0, 26.0/)
     opts_circ@cnFillColors                = (/6,8,0,12,13,14,15,\
                                               16,18,19,20,21,22,\
                                               23,24,25,26,27,28/)     
     elseif (r .eq. 0.24) then
     opts_circ@cnLevels                    = (/-0.5, 0.0, 0.5, 1.0, \
                                                1.5, 2.0, 3.0, 4.0, 5.0, 6.0, \
                                                7.0, 8.0, 9.0, 10.0, 11.0, 12.0, \
                                                14.0, 16.0/)
     opts_circ@cnFillColors                = (/6,8,0,12,13,14,15,\
                                               16,18,19,20,21,22,\
                                               23,24,25,26,27,28/)
     end if 

     opts_circ@gsnPaperOrientation         = "landscape"
     opts_circ@tiMainFontHeightF           = 0.0125
     opts_circ@gsnLeftString               = ""
     opts_circ@gsnRightString              = ""

     opts_circ@lbLabelBarOn                = True    ; Keep label bar on
     opts_circ@pmLabelBarWidthF            = 0.6
     opts_circ@pmLabelBarHeightF           = 0.15
     opts_circ@lbLabelFontHeightF          = 0.01
     opts_circ@lbPerimOn                   = False
     opts_circ@lbLabelFont                 = "Helvetica"
     opts_circ@lbTitleOn                   = True
     opts_circ@lbTitleString 		   = "Normalised circulation (s~S~-1~N~)"
     opts_circ@lbTitleFontHeightF          = 0.015
     opts_circ@lbTitleDirection            = "Across"
     opts_circ@lbBoxEndCapStyle             = "TriangleBothEnds"    ; Labelbar end shape

     opts_circ@gsnMaximize                 = True
     opts_circ@tmXTOn                      = "False" ; Turn off top x-axis TM
     opts_circ@tmYROn                      = "False" ; Turn off right y-axis TM
     opts_circ@tiYAxisString 		   = "Pressure (hPa)"
     opts_circ@tiXAxisString 		   = ""	     ; No label beneath x-axis
     opts_circ@vpWidthF                    = 0.6     ; Plot width
     opts_circ@vpHeightF                   = 0.5     ; Plot height
     opts_circ@tmXBMode                    = "Explicit"
     opts_circ@tmXBLabelFontHeightF        = 0.01

   ; Create labels based on length of integration
     tt0 = tf0 - ts0
     if (tt0 .eq. 24) then
      opts_circ@tmXBValues                  = (/ts0, ts0+12, tf0/)
     elseif (tt0 .eq. 36) then
      opts_circ@tmXBValues                  = (/ts0, ts0+12, ts0+24, tf0/)
     elseif (tt0 .eq. 48)
      opts_circ@tmXBValues                  = (/ts0, ts0+12, ts0+24, ts0+36, tf0/)
     elseif (tt0 .eq. 60)
      opts_circ@tmXBValues                  = (/ts0, ts0+12, ts0+24, ts0+36, ts0+48, tf0/)
     elseif (tt0 .eq. 72)
      opts_circ@tmXBValues                  = (/ts0, ts0+24, ts0+48, tf0/)
     elseif (tt0 .eq. 84)
      opts_circ@tmXBValues                  = (/ts0, ts0+24, ts0+48, ts0+72/)
     elseif (tt0 .eq. 96)
      opts_circ@tmXBValues                  = (/ts0, ts0+24, ts0+48, ts0+72, tf0/)
     end if

   ; Retrieve correct indices (from 'time_str') corresponding to the times above
     lab_ind = opts_circ@tmXBValues(:) - 3

   ; Create x-axis labels from values in 'time_str'
     opts_circ@tmXBLabels                  = time_str(lab_ind)

     opts_circ@trXMaxF                     = tf0            ; X-axis maximum (12Z 8 Jul)
     opts_circ@trXMinF                     = ts0            ; X-axis minimum (00Z 5 Jul)

     opts_circ@trYMaxF                     = 1000    ; Include boundary	layer
     opts_circ@trYMinF                     = 150     ; Include upper troposphere

     opts_circ@tiMainString                = ""
     opts_circ@tiXAxisFontHeightF          = 0.0125
     opts_circ@tiYAxisFontHeightF          = 0.0125
     opts_circ@gsnDraw                     = False   ; Do not draw the plot
     opts_circ@gsnFrame                    = False   ; Do not advance the frame

     opts_ten				   = opts_circ
     delete(opts_ten@cnFillColors)
     delete(opts_ten@cnLevels)
     opts_ten@cnLevels                     = (/-4.0, -2.0, -1.0, -0.5, -0.2, 0.0, \
     					        0.2, 0.5, 1.0, 2.0, 4.0/)
     opts_ten@cnFillColors                 = (/2,4,6,7,8,0,12,14,15,16,18,19/)

;=================
; Plot the data
;=================

     panel(0) = gsn_csm_pres_hgt(wks,circ_instF0,opts_circ)

;=================================
; Additional panel plot options
;=================================

     opts_vortP                       = True
     opts_vortP@gsnFrame              = False
     opts_vortP@gsnPanelLabelBar      = False
     opts_vortP@txString              = ""
     opts_vortP@gsnPanelFigureStrings = ""
     opts_vortP@gsnMaximize           = True
     opts_vortP@gsnPanelTop           = 0.98
     opts_vortP@gsnPanelBottom        = 0.02
     opts_vortP@amJust                = "TopLeft"
     opts_vortP@gsnPanelFigureStringsFontHeightF = 0.01 ; Reduce label size (0.01)

     gsn_panel(wks,panel,(/1,1/),opts_vortP) ; Draw as a single plot
     frame(wks)

end