; Script to calculate circulation budget following a tropical storm

; Run using: 

; ncl dat=\"02T12\" opt=\"x11\" ens0=\"em11\" rad=0.76 dist=2.2 nbox=25 \
; fct=5 op=1 ts0=24 tf0=48 calc=1 cn0=\"slp\" tc_cbgt_new_5min.ncl

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'opt'   = output file format ("pdf" or "x11")
; 'ens0'  = ensemble member (em00 -- em11)
; 'dist'  = size of box following storm (degrees)
; 'ts0'   = starting time index (0-119)
; 'tf0'   = ending time index (0-119)
; 'nbox'  = number of ensemble box positions for circ. budget calculations (49,25,9,1)
; 'rad'   = radius of box around TC for budget calculations (integer multiple of dx,dy)
; 'fct'	  = 
; 'op'	  = 
; 'calc'  = method of finding storm centre
; 'cn0'	  = offline storm centre: "slp", "vort", "geo", "geo_sm"

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

  fili_d1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pd*.nc")
  fili_d      = diri+"/"+fili_d1

;=============================================================
; Calculate number of times in file using built-in function
;=============================================================

; 'pb' stream
  numINFO_b   = nc_times(fili_b)
  numINPUT_b  = numINFO_b[0]		; Number of times
  time_b      = numINFO_b[1]		; Array of times

; 'pd' stream
  numINFO_d   = nc_times(fili_d)
  numINPUT_d  = numINFO_d[0]            ; Number of times
  time_d      = numINFO_d[1]            ; Array of times

;===============================================
; Calculate number of vertical levels in file 
;===============================================

; 'pb' stream
  nLEVS_b     = nc_levs(fili_b)
  numLEVS_b   = nLEVS_b[1]              ; Number of vertical levels
  levs_b      = nLEVS_b[0]		; Array of vertical levels

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
  d0	      = dsize(0)			 ; Size of dimension 0
  d1	      = dsize(1)			 ; Size of dimension 1

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

  date_str  = sprinti("%0.2iUTC ", hour) + \
              sprinti("%0.2i ", day) \
            + month_abbr(month)
  time_str  = sprinti("%0.2i UTC ", hour) + \
              sprinti("%0.2i ", day) \
            + month_abbr(month)
  time_arr  = sprinti("%0.2i", day) + \
              month_abbr(month) + \
              "_" + sprinti("%0.2iZ", hour)

;==================================================================================
;;;;;;;;;;;;;;;;;;;;;;;;;;; Start of main part of script ;;;;;;;;;;;;;;;;;;;;;;;;;;
;==================================================================================

  print_clock("Starting main part of script...")

; Prelash [define arrays for use in loops below]

  pos_arr     = new((/numINPUT_a,2/),float) ; Array for cyclone position [lat/lon]
  r_earth     = 6.37e6                      ; Earth's radius (constant)
  plevs	      = numLEVS_b	            ; Pressure levels (17)

; Define constants
  rd    = 287.0   ; Gas constant for dry air (J/kg/K^2)
  cp0   = 1004.0  ; Specific heat of dry air at constant pressure (J/kg/K)
  kp    = 0.286   ; For Exner function calculations (0.286)
  pref  = 1000.0  ; Reference pressure p0 (1000 hPa)

; Create arrays to hold circulation tendency values 
; [ensemble member (nsim) * box ensemble (nbox) * time (numINPUT_a) * levs (plevs)]
  circ_inst = new((/nbox,numINPUT_a,plevs/),float)
  hadv_ten  = new((/nbox,numINPUT_a,plevs/),float)
  str_ten   = new((/nbox,numINPUT_a,plevs/),float)
  flux_ten  = new((/nbox,numINPUT_a,plevs/),float)
  tilt_ten  = new((/nbox,numINPUT_a,plevs/),float)
  fric_ten  = new((/nbox,numINPUT_a,plevs/),float)

;===============================
; Start loop over input files 
;===============================

  b	    = addfile(fili_b,"r")		; Read in 'pb' stream [u,v,w,vort,t,z,rh]
  d  	    = addfile(fili_d,"r")		; Read in 'pd' stream [w,p,t,theta,kmh,tau]

;==============================================
; Get the variables we will need [pb stream]
;==============================================

  do it = ts0, tf0				; Loop over desired times 

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

   pres	 = b->p				    ; Pressure (hPa)
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
   t2 	 = t(0:size(3,0)-1,0:size(3,1)-1,0:size(3,2)-1)
   z2 	 = z(0:size(4,0)-1,0:size(4,1)-1,0:size(4,2)-1)
   w2    = w(0:size(5,0)-1,0:size(5,1)-1,0:size(5,2)-1)
   lat2	 = lat(0:nind-1)
   lon2	 = lon(0:nind-1)
   delete([/u0,v0,vort,t,z,w,lat,lon/])

 ; Reset names to avoid further code changes 
   u0    = u2
   v0	 = v2
   vort	 = vort2
   t	 = t2
   z	 = z2
   w	 = w2
   lat	 = lat2
   lon	 = lon2
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

;==============================================
; Get the variables we will need [pd stream]
;==============================================

 ; Read in basic variables from 'pd' stream
   lon1D  = d->longitude_1({n0:n1})                ; longitude [xxx grid points]
   lat1D  = d->latitude_1({t0:t1})                 ; latitude  [xxx grid points]
   lonD   = d->longitude({n0:n1})                  ; longitude [xxx grid points]
   latD   = d->latitude({t0:t1})                   ; latitude  [xxx grid points]
   hybD   = d->hybrid_ht(:)                        ; 63 model ('rho') levels

 ; Read in diffusivity and wind stress
 ; 'str_xz' -- [hybD | 63] * [latD | 150]  * [lonD | 151]
 ; 'str_yz' -- [hybD | 63] * [lat1D | 150] * [lon1D | 151]
 ; 'kmh'    -- [hybD | 63] * [latD | 150]  * [lon1D | 151]
   str_xz = d->taux(it-1,:,{t0:t1},{n0:n1})        ; Wind stress, x-comp (N m/2)
   str_yz = d->tauy(it-1,:,{t0:t1},{n0:n1})        ; Wind stress, y-comp (N m/2)
   kmh    = d->unspecified(it-1,:,{t0:t1},{n0:n1}) ; Horizontal diffusivity (m^2 s-1)

 ; Read in temperature, pressure, vertical velocity and geopotential height
 ; All variables below -- [hybD | 63] * [latD | 150] * [lon1D | 151]
   td     = d->temp(it-1,:,{t0:t1},{n0:n1})        ; Temperature (K)
   pd     = d->p(it-1,:,{t0:t1},{n0:n1})           ; Pressure (Pa)
   wd     = d->dz_dt(it-1,:,{t0:t1},{n0:n1})       ; Vertical velocity (m s-1)

 ; Calculate dry air density from pressure and temperature
   rho    = pd / (rd * td)

 ; Add metadata
   copy_VarCoords(pd, rho)
   rho@description = "Dry air density"
   rho@units = "kg m~S~-3~N~"

 ; Also calculate specific volume (for use in diffusive tendency source terms)
   rho_a  = 1 / rho
   copy_VarCoords(rho, rho_a)
   rho_a@description = "Specific volume"
   rho_a@units = "m~S~3~N~ kg~S~-1~N~"

   pd      = pd / 100 	     ; Convert units of pressure before continuing
   pd@units = "hPa"	     ; Add metadata

;==========================================================================
; Define new arrays for terms in vorticity equation (no friction for now)
;==========================================================================

   pii        = 3.14159265           ; Constant
   lat_rad    = lat * (pii/180)      ; Convert to radians for calculation below
   lon_rad    = lon * (pii/180)

 ; Troubleshooting (look for abnormally large values and exit the script if found)
   if( any (u .gt. 200) ) then
    print("Quitting because of storm-relative velocity error!!!")
   end if

   d_vort     = dimsizes(vort)
   d_1	      = d_vort(1)
   d_2	      = d_vort(2)

 ; Create new arrays for horizontal components of tilting term
   if (d_1 .eq. d0 .and. d_2 .eq. d1) then
    tilt_x     = new((/plevs,d0,d1/),float)
    tilt_y     = new((/plevs,d0,d1/),float)
   else
    tilt_x     = vort
    tilt_y     = vort
   end if 

 ; Add basic metadata (modify details later)
   copy_VarMeta(vort,tilt_x)
   copy_VarAtts(vort,tilt_x)
   copy_VarCoords(vort,tilt_x)

   copy_VarMeta(vort,tilt_y)
   copy_VarAtts(vort,tilt_y)
   copy_VarCoords(vort,tilt_y)

;==========================================================
; Define latitude and longitude for finite differencing
;==========================================================

 ; Create lat/lon grid with same dimensions as all input variables
   if (d_1 .eq. d0 .and. d_2 .eq. d1) then
    y2d      = new((/d0,d1/),float)	
   else 
    y2d	     = vort(0,:,:)   
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
   ddx 	     = y2d
   ddx 	     = r_earth*(pii/180.0)*delta_lon*cos(y2d*pii/180)
   ddy 	     = y2d
   ddy 	     = r_earth*(pii/180)*delta_lat

 ; Give 'coslat' the dimensions of 'u' (3-D) and hold values in new array ('coslat3d')
   coslat3d  = conform_dims(dimsizes(u),coslat,(/1,2/))

;=================================
; Friction-related calculations 
;=================================

 ; Array of hybrid levels to approximately match pressure levels in 'pres'
 ; Array [/1000, 950, 925, 900, 850, 800, 750, 700, \
 ;         650, 600, 500, 400, 300, 250, 200, 150, 100/]

   lev_arr = (/1, 8, 13, 16, 19, 21, 23, 25, 27, 29, \
               31, 33, 35, 37, 39, 42, 47, 51, 56/)
   lev_arr@units = "m"

   l0    = 0  ; Subset of vertical levels (19 levels)
   l1    = 18 ; Level 56 of original array ~ 100 hPa (Int. Standard Atmosphere)

;===============================================================
; Tilting term (w * dv_dp, -w * du_dp)
;===============================================================

   if (d_1 .eq. d0 .and. d_2 .eq. d1) then
    du_dp      = new((/plevs,d0,d1/),float)
    dv_dp      = new((/plevs,d0,d1/),float)
   else
    du_dp      = vort
    dv_dp      = vort
   end if

 ; Convert pressure values to Pa for calculation below
   pres = pres * 100
   pres@units = "Pa"

 ; Fill arrays with vertical derivatives at each pressure level
   do pn = 0, plevs-1

    kp1 = min((/pn+1,plevs-1/)) ; Account for top and bottom levels
    km1 = max((/pn-1,0/))

    du_dp(pn,:,:) = (u(kp1,:,:) - u(km1,:,:)) / ((pres(km1))-(pres(kp1)))
    dv_dp(pn,:,:) = (v(kp1,:,:) - v(km1,:,:)) / ((pres(km1))-(pres(kp1)))

   end do

 ; Convert pressure values back to hPa
   pres = pres / 100
   pres@units = "hPa"

 ; Calculate [dv_dp * w] -- x-component of tilting term
   tilt_y = dv_dp * w

 ; Calculate [-du_dp * w] -- y-component of tilting term
   tilt_x = -du_dp * w

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

 ; Get information about box size for later ensemble calculations
   r = rad   ; Box radius (degrees)
   box = ddx({centre(cn,it,0)-r:centre(cn,it,0)+r},{centre(cn,it,1)-r:centre(cn,it,1)+r})
   box_size = dimsizes(box)
   size0    = box_size(0)
   size1    = box_size(1)

   s0       = size0-1
   s1	    = size1-1

;============================================================
; Loop over pressure levels to calculate vertical profile
;============================================================

   cp = 0			     ; Counter variable
 
   do ilev = 0, plevs-1		     ; Loop over chosen pressure levels

    p = pres(ilev)
    print("Working on: " + p + " hPa")

;=================================
; Pressure levels ['pb' stream]
;=================================

  ; Read in storm-relative wind on pressure levels
    u_plane       = u(ilev,0:d0-1,0:d1-1)
    v_plane       = v(ilev,0:d0-1,0:d1-1)

  ; Read in other variables on pressure levels
    vort_plane    = vort(ilev,0:d0-1,0:d1-1)
    avo_plane	  = avo(ilev,0:d0-1,0:d1-1)
    div_plane	  = div(ilev,0:d0-1,0:d1-1)
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
     lon_max0 	     = n0 + (dx0 * inds(0,1))

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

 ; Terms in the circulation budget
 ; 'tilt_x' and 'tilt_y' contain values on all pressure levels for each time
 ; 'tilt_x_plane' and 'tilt_y_plane' contain values on single pressure levels
   tilt_x_plane  = tilt_x(ilev,:,:)       ; Tilting (x-component)
   tilt_y_plane  = tilt_y(ilev,:,:)       ; Tilting (y-component)

;=====================================================================================
; Create array of box positions either side of the original [Davis & Galarneau 2009]
;=====================================================================================

   if (nbox .eq. 49) then
    lat_box   = (/lat_max0-(3*dy0),lat_max0-(2*dy0),lat_max0-(dy0),lat_max0,\
                  lat_max0+(dy0),lat_max0+(2*dy0),lat_max0+(3*dy0)/)
    lon_box   = (/lon_max0-(3*dx0),lon_max0-(2*dx0),lon_max0-(dx0),lon_max0,\
                  lon_max0+(dx0),lon_max0+(2*dx0),lon_max0+(3*dx0)/)
   elseif (nbox .eq. 25) then
    lat_box   = (/lat_max0-(2*dy0),lat_max0-(dy0),lat_max0,\
   	          lat_max0+(dy0),lat_max0+(2*dy0)/)
    lon_box   = (/lon_max0-(2*dx0),lon_max0-(dx0),lon_max0,\
                  lon_max0+(dx0),lon_max0+(2*dx0)/)
   elseif (nbox .eq. 9) then
    lat_box   = (/lat_max0-(dy0),lat_max0,lat_max0+(dy0)/)
    lon_box   = (/lon_max0-(dx0),lon_max0,lon_max0+(dx0)/)
   elseif (nbox .eq. 1) then
    lat_box   = (/lat_max0/)
    lon_box   = (/lon_max0/)
   end if

   lsize      = dimsizes(lat_box) * dimsizes(lon_box)

;====================================================
; Define arrays to hold ensemble of box positions
;====================================================

 ; EDIT --> NEATER WAY TO INITIALISE MULTIPLE ARRAYS (???)
 ; Ensemble of boxes (all grid points, not just box centres)
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

 ; Initialise array to hold perimeter averages (for each side of box)
   vort_avg   = new((/lsize,5/),float)
   avo_avg    = new((/lsize,5/),float)
   u_avg      = new((/lsize,5/),float)
   v_avg      = new((/lsize,5/),float)

;=============================================
; Calculate terms in the circulation budget
;=============================================

 ; Radius of box for calculation (degrees)
   r          = rad
   r1 	      = sprintf("%0.1fdeg",2*r) ; Create string to add to output file

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

   ; Read box edge coordinates into 'lat_d' and 'lon_d' arrays
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
     d_ave(a) = ( (dz0(a) + dz1(a) ) / 2) * dm0
     d_ave@units = "m~S~-2~N~"

   ; Also calculate area of box using 'gc_latlon' function
     ddz0 = gc_latlon( lat_d(a), lon_d(a), lat_d(a), lon_u(a), 2, 3)	; Bottom
     ddz1 = gc_latlon( lat_d(a), lon_u(a), lat_u(a), lon_u(a), 2, 3)	; Right
     ddz2 = gc_latlon( lat_u(a), lon_u(a), lat_u(a), lon_d(a), 2, 3)	; Top	
     ddz3 = gc_latlon( lat_u(a), lon_d(a), lat_d(a), lon_d(a), 2, 3)	; Left
     aaz0 = ( (ddz0 + ddz2) / 2) * ddz1	   	     	       	        ; Area

   ; Troubleshooting 
     if (tr .eq. 1) then 
      print("Right boundary length = "+dm0+" "+dm0@units+" ("+ddz1+" m)")     
      print("Left boundary length = "+dm0+" "+dm0@units+" ("+ddz3+" m)")
      print("Bottom boundary length = "+dz0+" "+dz0@units+" ("+ddz0+" m)")
      print("Top boundary length = "+dz1+" "+dz1@units+" ("+ddz2+" m)")
      print("Area of box = "+d_ave+" "+d_ave@units+" ("+aaz0+" m)")
     end if 

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

;===========================================
; Calculate averages along each perimeter
;===========================================

   ; Now calculate terms along the four edges of the box, for each ensemble member
   ; Hence four arrays required for each term in the budget

   ; All arrays have dimensions [lsize * 13] -> ALL members on a SINGLE pressure level

   ; South 
     vort_avg(a,0) = avg(vort_box(a,0,:))
     avo_avg(a,0)  = avg(avo_box(a,0,:))
     u_avg(a,0)	   = avg(u_box(a,0,:))
     v_avg(a,0)	   = avg(v_box(a,0,:))

   ; East
     vort_avg(a,1) = avg(vort_box(a,:,numPOINTS_x-1))
     avo_avg(a,1)  = avg(avo_box(a,:,numPOINTS_x-1))
     u_avg(a,1)    = avg(u_box(a,:,numPOINTS_x-1))
     v_avg(a,1)    = avg(v_box(a,:,numPOINTS_x-1))

   ; North
     vort_avg(a,2) = avg(vort_box(a,numPOINTS_y-1,:))
     avo_avg(a,2)  = avg(avo_box(a,numPOINTS_y-1,:))
     u_avg(a,2)    = avg(u_box(a,numPOINTS_y-1,:))
     v_avg(a,2)    = avg(v_box(a,numPOINTS_y-1,:))

   ; West 
     vort_avg(a,3) = avg(vort_box(a,:,0))
     avo_avg(a,3)  = avg(avo_box(a,:,0))
     u_avg(a,3)    = avg(u_box(a,:,0))
     v_avg(a,3)    = avg(v_box(a,:,0))

     if (tr .eq. 1) then 
      print("Average vort: "+vort_avg(a,0)+" (S), "+vort_avg(a,1)+\
      		           " (N), "+vort_avg(a,2)+" (E), "+vort_avg(a,3)+" (W)")
      print("Average avo: "+avo_avg(a,0)+" (S), "+avo_avg(a,1)+\
      		    	  " (N), "+avo_avg(a,2)+" (E), "+avo_avg(a,3)+" (W)")	
      print("Average u: "+u_avg(a,0)+" (S), "+u_avg(a,1)+\
      		        " (N), "+u_avg(a,2)+" (E), "+u_avg(a,3)+" (W)")
      print("Average v: "+v_avg(a,0)+" (S), "+v_avg(a,1)+\
      		        " (N), "+v_avg(a,2)+" (E), "+v_avg(a,3)+" (W)")
     end if

   ; Also calculate average around the whole box
     vort_avg(a,4) = avg(vort_avg(a,0:3) )
     avo_avg(a,4)  = avg(avo_avg(a,0:3) )
     u_avg(a,4)    = avg(u_avg(a,0:3) )
     v_avg(a,4)    = avg(v_avg(a,0:3) )

;====================================================================================
; Second round of calculations (Horizontal advection, relative vorticity, tilting)
;====================================================================================

   ; Calculate terms on northern/southern boundaries

     if (int .eq. 1) then 

    ; Horizontal advection (into the box)
      hadv_s      = dim_sum_n_Wrap((-(avo_box(a,0,:) - avo_avg(a,0)) * \
                                   (v_box(a,0,:) - v_avg(a,0)) * \
                                    ddx_box(a,0,:)), 0)
      hadv_n	  = dim_sum_n_Wrap((-(avo_box(a,numPOINTS_y-1,:) - avo_avg(a,2)) * \
                                   (v_box(a,numPOINTS_y-1,:) - v_avg(a,2)) * \
                                    ddx_box(a,numPOINTS_y-1,:)), 0)

    ; Tilting
      tilt_x0     = dim_sum_n_Wrap((tilt_x_box(a,:,:) * ddx_box(a,:,:)),1)
      tilt_s      = tilt_x0(0)              ; Tilting (south)
      tilt_n      = tilt_x0(numPOINTS_y-1)  ; Tilting (north)

    ; Circulation
      circ_x0     = dim_sum_n_Wrap((u_box(a,:,:) * ddx_box(a,:,:)),1)
      circ_s      = circ_x0(0)              ; Circulation (south)
      circ_n      = -circ_x0(numPOINTS_y-1) ; Circulation (north)

     elseif (int .eq. 0) then 

    ; Tilting 
      tilt_x0     = dim_sum_n_Wrap((tilt_x_box(a,:,:) * ddx_box(a,:,:)),1)

    ; Horizontal advection (into the box, across the perimeter)
      hadv_x0     = dim_sum_n_Wrap((-(avo_box(a,:,:) - avo_avg(a,4)) * \
                                   (v_box(a,:,:) - v_avg(a,4)) * \
                                    ddx_box(a,:,:)),1)
      circ_x0     = dim_sum_n_Wrap((u_box(a,:,:) * ddx_box(a,:,:)),1)
      tilt_s      = tilt_x0(0)              ; Tilting (south)
      tilt_n      = tilt_x0(numPOINTS_y-1)  ; Tilting (north)
      hadv_s      = hadv_x0(0)              ; Horizontal advection (south)
      hadv_n      = -hadv_x0(numPOINTS_y-1) ; Horizontal advection (north)
      circ_s      = circ_x0(0)              ; Circulation (south)
      circ_n      = -circ_x0(numPOINTS_y-1) ; Circulation (north)	
      delete(hadv_x0)

     end if 

   ; Calculate terms on eastern/western boundaries

     if (int .eq. 1) then

    ; Horizontal advection (into the box)
      hadv_w      = dim_sum_n_Wrap((-(avo_box(a,:,0) - avo_avg(a,3)) * \
                                   (v_box(a,:,0) - v_avg(a,3)) * \
                                    ddx_box(a,:,0)), 0)
      hadv_e      = dim_sum_n_Wrap((-(avo_box(a,:,numPOINTS_y-1) - avo_avg(a,1)) * \
                                   (v_box(a,:,numPOINTS_y-1) - v_avg(a,1)) * \
                                    ddx_box(a,:,numPOINTS_y-1)), 0)
   
    ; Tilting 
      tilt_y0     = dim_sum_n_Wrap((tilt_y_box(a,:,:) * ddy_box(a,:,:)),0)
      tilt_w      = tilt_y0(0)              ; Tilting (west)
      tilt_e      = tilt_y0(numPOINTS_x-1)  ; Tilting (east)

    ; Circulation
      circ_y0     = dim_sum_n_Wrap((v_box(a,:,:) * ddy_box(a,:,:)),0)
      circ_w      = -circ_y0(0)             ; Circulation (west)
      circ_e      = circ_y0(numPOINTS_x-1)  ; Circulation (east)

     elseif (int .eq. 0) then 

      tilt_y0     = dim_sum_n_Wrap((tilt_y_box(a,:,:) * ddy_box(a,:,:)),0)
      hadv_y0     = dim_sum_n_Wrap((-(avo_box(a,:,:) - avo_avg(a,4)) * \
                                   (u_box(a,:,:) - u_avg(a,4)) * \
                                    ddy_box(a,:,:)),0)
      circ_y0     = dim_sum_n_Wrap((v_box(a,:,:) * ddy_box(a,:,:)),0)
      tilt_w      = tilt_y0(0)              ; Tilting (west)
      tilt_e      = tilt_y0(numPOINTS_x-1)  ; Tilting (east)
      hadv_w      = hadv_y0(0)              ; Horizontal advection (west)
      hadv_e      = -hadv_y0(numPOINTS_x-1) ; Horizontal advection (east)
      circ_w      = -circ_y0(0)             ; Circulation (west)
      circ_e      = circ_y0(numPOINTS_x-1)  ; Circulation (east)
      delete(hadv_y0)

     end if 

;========================
; Further calculations
;========================

   ; Calculate stretching term using averaged quantities
     str_tenT   = -(div_ave(a) * avo_avg(a,4))

   ; Calculate horizontal advection from line integrals (4 sides of box)
     hadv_tenT  = hadv_s + hadv_n + hadv_w + hadv_e
     hadv_tenT  = hadv_tenT / d_ave(a) ; Divide by area of box [Rios-Berrios et al. 2016]

   ; Calculate tilting tendency from line integrals (4 sides of box)
     tilt_tenT  = tilt_s + tilt_n + tilt_w + tilt_e
     tilt_tenT  = tilt_tenT / d_ave(a) ; Divide by area of box [Rios-Berrios et al. 2016]

   ; Calculate circulation from line integrals (4 sides of box)
     circ_instT0 = circ_s + circ_n + circ_w + circ_e	     ; Sum 4 line integrals
     circ_instT0 = circ_instT0 / d_ave(a)	             ; Divide by area of box 
     circ_instT  = vort_ave(a) 				     ; Area-averaged vorticity

     if (tr .eq. 1) then 
      print("Circulation = "+circ_instT0+" (v); "+circ_instT+" (vort)")
     end if

   ; Enter values into predefined arrays
     hadv_ten(a,it,cp)      = hadv_tenT		     ; Horizontal advection
     str_ten(a,it,cp)       = str_tenT		     ; Stretching
     tilt_ten(a,it,cp)      = tilt_tenT		     ; Tilting
     circ_inst(a,it,cp)     = circ_instT             ; Circulation (LHS)

   ; Tidy up (delete variables from within lat/lon loops before next iteration)
     delete([/hadv_tenT, str_tenT, tilt_tenT, circ_instT, circ_instT0/])
     delete([/numPOINTS_x0, numPOINTS_x, numPOINTS_y/])
     delete([/tilt_y0,tilt_x0,circ_y0,circ_x0/])

     end do     ; End latitude loop (do j = 0, dimsizes(lat_box)-1)
    end do	; End longitude loop (do i = 0, dimsizes(lon_box)-1)

    cp                   = cp + 1       ; Counter variable (pressure)

  ; Tidy up (delete variables from within pressure level loop)
    delete([/u_plane, v_plane, tilt_x_plane, tilt_y_plane/])
    delete([/vort_plane, avo_plane, div_plane/])
    delete([/ddx_box, ddy_box, u_box, v_box/])
    delete([/tilt_x_box, tilt_y_box, vort_box, avo_box, div_box, lat_box, lon_box/])
    delete([/div_ave,vort_ave,vort_avg,avo_avg,u_avg,v_avg/])

   end do ; End pressure level loop (do ilev = 0, plevs-1)

 ; Tidy up (delete variables from within time loop before next iteration)
   delete([/lon_rad,lat_rad,lon,lat,lon1,lat1,u0,v0,u,v,vort,div,f0,f1,avo,t,z,w,box/])
   delete([/lon1D,lat1D,lonD,latD,hybD,str_xz,str_yz,kmh,td,p,pd,wd,rho,rho_a/])

  end do                ; End time loop (do it = 0, times-1)

;=====================================================
; Calculate term averages over all ensemble members
;=====================================================

   print_clock("Finished ensemble calculations! Now averaging over all members...")

   hadv_final = new((/numINPUT_a,plevs/),"float")
   str_final  = new((/numINPUT_a,plevs/),"float")
   tilt_final = new((/numINPUT_a,plevs/),"float")
   circ_final = new((/numINPUT_a,plevs/),"float")
   circ_ten   = new((/numINPUT_a,plevs/),"float")

 ; Average over ensemble members (do we still need this step?)
   hadv_final(:,:) = dim_avg_n_Wrap(hadv_ten(:,:,:),0)
   str_final(:,:)  = dim_avg_n_Wrap(str_ten(:,:,:),0)
   tilt_final(:,:) = dim_avg_n_Wrap(tilt_ten(:,:,:),0)
   circ_final(:,:) = dim_avg_n_Wrap(circ_inst(:,:,:),0)

 ; Get consistent units before next set of calculations 
   if (fct .eq. 4) then 
    circ_final = circ_final * (10 ^ 4)
   elseif (fct .eq. 5) then 
    circ_final = circ_final * (10 ^ 5)
   end if 

;=========================================
; Integrate over Nepartak's life cycle
;=========================================

 ; Each calculation produces an array of [plevs] values [pressure levels]
 ; We do this calculation [numINPUT_a] times
 ; Create arrays with [numINPUT_a] * [plevs] elements

   hadv_intA = new((/numINPUT_a,plevs/),float)
   str_intA  = new((/numINPUT_a,plevs/),float)
   tilt_intA = new((/numINPUT_a,plevs/),float)
   tot_intA  = new((/numINPUT_a,plevs/),float)
   circ_intA = new((/numINPUT_a,plevs/),float)

 ; Create array of time intervals
   time_int = new(numINPUT_a,"integer")
   time_int = ispan(0,numINPUT_a-1,1)
   time_int!0 = "time"

   print_clock("Starting plots...")

   do ic = 0, numINPUT_a-1  ; 1 h integrations over the whole of Nepartak's life cycle

    cp1 = min((/ic+1,numINPUT_a-1/))
    cm1 = max((/ic-1,0/))

     f_str = sprinti("%0.3i",ic)+"h"    ; Formatted string for output file

     hadv_intA(ic,:) = dim_sum_n_Wrap((hadv_final(ic:cp1,:) * 3600), 0)
     str_intA(ic,:)  = dim_sum_n_Wrap((str_final(ic:cp1,:) * 3600), 0)
     tilt_intA(ic,:) = dim_sum_n_Wrap((tilt_final(ic:cp1,:) * 3600), 0)
     circ_intA(ic,:) = (circ_final(cp1,:) - circ_final(cm1,:)) / 7200

     hadv_intA(ic,:) = hadv_intA(ic,:) * (10 ^ fct)
     str_intA(ic,:)  = str_intA(ic,:) * (10 ^ fct)
     tilt_intA(ic,:) = tilt_intA(ic,:) * (10 ^ fct)
     circ_intA(ic,:) = circ_intA(ic,:) * (10 ^ fct)

   end do               ; End time loop (do ic = 0, numINPUT_a-1)

 ; Add metadata after calculating tendencies 
   hadv_intA@units = "10~S~-"+fct+"~N~ s~S~-2~N~"
   str_intA@units  = "10~S~-"+fct+"~N~ s~S~-2~N~"
   tilt_intA@units = "10~S~-"+fct+"~N~ s~S~-2~N~"
   tot_intA@units  = "10~S~-"+fct+"~N~ s~S~-2~N~"
   circ_intA@units = "10~S~-"+fct+"~N~ s~S~-2~N~"

   tot_intA  = str_intA + hadv_intA + tilt_intA
   
 ; Add metadata before time-height plotting
   circ_intA!0    = "time"
   circ_intA&time = time_int
   circ_intA!1    = "lev"
   circ_intA&lev  = pres

   tot_intA!0     = "time"
   tot_intA&time  = time_int
   tot_intA!1     = "lev"
   tot_intA&lev   = pres

   hadv_intA!0    = "time"
   hadv_intA&time = time_int
   hadv_intA!1    = "lev"
   hadv_intA&lev  = pres

   str_intA!0     = "time"
   str_intA&time  = time_int
   str_intA!1     = "lev"
   str_intA&lev   = pres

   tilt_intA!0    = "time"
   tilt_intA&time = time_int
   tilt_intA!1    = "lev"
   tilt_intA&lev  = pres

;======================================
; Reorder dimensions before plotting
;======================================

 ; Define new arrays for final plotting 
   circ_int = new((/plevs,numINPUT_a/),float)
   hadv_int = new((/plevs,numINPUT_a/),float)
   str_int  = new((/plevs,numINPUT_a/),float)
   tilt_int = new((/plevs,numINPUT_a/),float)
   tot_int  = new((/plevs,numINPUT_a/),float)

 ; Reorder from [time * pressure] to [pressure * time]
   circ_int(:,:) = circ_intA(lev|:,time|:)
   hadv_int(:,:) = hadv_intA(lev|:,time|:)
   str_int(:,:)  = str_intA(lev|:,time|:)
   tilt_int(:,:) = tilt_intA(lev|:,time|:)
   tot_int(:,:)  = tot_intA(lev|:,time|:)

 ; Calculate running average for LHS 
   circ_int1     = runave_n_Wrap(circ_int, 3, 0, 1)

;======================================
; Create an x-y time series plot
;======================================

     output = "$sam/nepartak/images/cart_coords/cbgt_xz_"+dat+"_"+ens0+"_"+r1+"_"+cn0
     wks    = gsn_open_wks(opt,output)

     panel = new(5,graphic)
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

     if (fct .eq. 4 .and. r .eq. 1.0) then
      opts_bdgt@cnLevels                    = (/-2.0, -1.5, -1.0, -0.5, \
                                                -0.3, -0.2, -0.1, -0.05, 0.05, \
                                                 0.1, 0.2, 0.3, 0.5, 1.0, 1.5, 2.0/)
     elseif (fct .eq. 4 .and. r .eq. 0.76) then
      opts_bdgt@cnLevels                    = (/-3.0, -2.0, -1.5, -1.0, \
                                                -0.5, -0.2, -0.1, -0.05, 0.05, \
                                                 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0/)
     elseif (fct .eq. 4 .and. r .eq. 0.48) then
      opts_bdgt@cnLevels                    = (/-3.0, -2.0, -1.5, -1.0, \
                                                -0.5, -0.2, -0.1, -0.05, 0.05, \
                                                 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0/)
     elseif (fct .eq. 4 .and. r .eq. 0.24) then
      opts_bdgt@cnLevels                    = (/-4.0, -3.0, -2.0, -1.5, \
                                                -1.0, -0.5, -0.2, -0.1, 0.1, \
                                                 0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0/)
     elseif (fct .eq. 5 .and. op .eq. 0) then
      opts_bdgt@cnLevels                    = (/-30.0, -25.0, -20.0, -15.0, -10.0, \
                                                -5.0, -2.0, -1.0, 1.0, 2.0, \
                                                 5.0, 10.0, 15.0, 20.0, 25.0, 30.0/)
     elseif (fct .eq. 5 .and. op .eq. 1) then
      opts_bdgt@cnLevels                    = (/-25.0, -20.0, -15.0, -10.0, -5.0, \
                                                -2.0, -1.0, -0.5, 0.5, 1.0, \
                                                 2.0, 5.0, 10.0, 15.0, 20.0, 25.0/)
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
     opts_bdgt@tiYAxisString               = "Pressure (hPa)"
     opts_bdgt@tiYAxisFontHeightF          = 0.030   ; Larger y-axis title
     opts_bdgt@tiXAxisString               = "Time"  ; X-axis title
     opts_bdgt@tiXAxisFontHeightF          = 0.030   ; Larger x-axis title
     opts_bdgt@vpWidthF                    = 0.6     ; Plot width
     opts_bdgt@vpHeightF                   = 0.5     ; Plot height
     opts_bdgt@tmXBMode                    = "Explicit"
     opts_bdgt@tmXBLabelFontHeightF        = 0.025

   ; Create labels based on length of integration
     tt0 = tf0 - ts0
     if (tt0 .eq. 24) then 
      opts_bdgt@tmXBValues                  = (/ts0, ts0+12, tf0/)
     elseif (tt0 .eq. 36) then 
      opts_bdgt@tmXBValues                  = (/ts0, ts0+12, ts0+24, tf0/)
     elseif (tt0 .eq. 48) 
      opts_bdgt@tmXBValues                  = (/ts0, ts0+12, ts0+24, ts0+36, tf0/)
     elseif (tt0 .eq. 60)
      opts_bdgt@tmXBValues                  = (/ts0, ts0+12, ts0+24, ts0+36, ts0+48, tf0/)
     elseif (tt0 .eq. 72)
      opts_bdgt@tmXBValues                  = (/ts0, ts0+12, ts0+24, ts0+36, \
      					      	ts0+48, ts0+60, tf0/)
     elseif (tt0 .eq. 84)
      opts_bdgt@tmXBValues                  = (/ts0, ts0+12, ts0+24, ts0+36, ts0+48, \
      					        ts0+60, ts0+72, tf0/)
     elseif (tt0 .eq. 96)
      opts_bdgt@tmXBValues                  = (/ts0, ts0+12, ts0+24, ts0+36, ts0+48, \
      					        ts0+60, ts0+72, ts0+84, tf0/)
     end if

   ; Retrieve correct indices (from 'time_str') corresponding to the times above
     lab_ind = opts_bdgt@tmXBValues(:) - 3

   ; Create x-axis labels from values in 'time_str'
     opts_bdgt@tmXBLabels                  = time_str(lab_ind)

     opts_bdgt@trXMaxF			   = tf0	    ; X-axis maximum (12Z 8 Jul)
     opts_bdgt@trXMinF			   = ts0	    ; X-axis minimum (00Z 5 Jul)

     opts_bdgt@trYMaxF                     = 1000           ; Ignore BL
     opts_bdgt@trYMinF                     = 150            ; Include upper troposphere

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

     if (r .eq. 0.76) then
      opts_cont@cnLevels                   = (/10.0, 14.0, 18.0, 22.0, 26.0/)
     elseif (r .eq. 0.48) then
      opts_cont@cnLevels                   = (/6.0, 8.0, 10.0, 14.0, 18.0/)
     elseif (r .eq. 0.24) then
      opts_cont@cnLevels                   = (/2.0, 4.0, 6.0, 8.0, 10.0/)
     end if

     opts_cont@cnLineThicknessF             = 1.5           ; Contour line thickness
     opts_cont@cnLineLabelsOn               = False         ; Turn off contour labels
     opts_cont@gsnContourNegLineDashPattern = 2             ; Negative contours dashed
     opts_cont@gsnContourPosLineDashPattern = 0             ; Positive contours solid
     opts_cont@cnInfoLabelOn                = False         ; Turn off info label
     circ_int@units                         = ""            ; Turn off string above plot

;=================
; Plot the data
;=================

     plot_circ0  = gsn_csm_pres_hgt(wks,circ_int(:,:),opts_bdgt)
     plot_tot0    = gsn_csm_pres_hgt(wks,tot_int(:,:),opts_bdgt)
     plot_str0    = gsn_csm_pres_hgt(wks,str_int(:,:),opts_bdgt)	
     plot_hadv0	  = gsn_csm_pres_hgt(wks,hadv_int(:,:),opts_bdgt)
     plot_tilt0   = gsn_csm_pres_hgt(wks,tilt_int(:,:),opts_bdgt)
;     plot_fric0   = gsn_csm_pres_hgt(wks,fric_int(:,:),opts_bdgt)

   ; First row
     panel(0)  = plot_circ0 
     panel(1)  = plot_tot0
     panel(2)  = plot_str0
     panel(3)  = plot_hadv0
     panel(4)  = plot_tilt0

;     panel(5)  = plot_fric0

     optsP                       = True
     optsP@gsnFrame              = False    ; Do not advance the frame
     optsP@gsnPanelLabelBar      = True     ; Turn on panel labelbar
     optsP@pmLabelBarWidthF      = 0.6
     optsP@pmLabelBarHeightF     = 0.15
     optsP@lbLabelFontHeightF    = 0.01
     optsP@lbPerimOn             = False
     optsP@lbLabelFont           = "Helvetica"
     optsP@lbBoxEndCapStyle      = "TriangleBothEnds"    ; Labelbar end shape

     if (fct .eq. 4) then
      optsP@lbTitleString = "Circulation budget tendency terms (10~S~-4~N~ s~S~-2~N~)"
     elseif (fct .eq. 5) then
      optsP@lbTitleString = "Circulation budget tendency terms (10~S~-5~N~ s~S~-2~N~)"
     end if 

     optsP@lbTitleFontHeightF    = 0.015
     optsP@lbTitleDirection      = "Across"
     optsP@lbTitlePosition       = "Bottom"
     optsP@txString              = ""
;     optsP@gsnPanelFigureStrings = (/"a) Total RHS", \
;                                     "b) Stretching", \
;				     "c) Horizontal advection", \
;                                     "d) Tilting"/)
     optsP@gsnMaximize           = True
     optsP@gsnPanelTop           = 0.98
     optsP@gsnPanelBottom        = 0.02
     optsP@gsnPaperOrientation	 = "landscape"
     optsP@amJust                = "TopLeft"
     optsP@gsnPanelFigureStringsFontHeightF = 0.0099 ; Label size (default 0.01)

     gsn_panel(wks,panel,(/1,5/),optsP)              ; Draw as a single plot
     frame(wks)

end