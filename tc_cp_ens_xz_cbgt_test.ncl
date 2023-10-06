; Script to calculate circulation budget terms every hour for CP ensemble simulations
; To run script, use command line options ['ens', 'opt', rad, run, nave]
; e.g. [ncl p1=1 p2=2 int=1 'opt="pdf"' rad=3 run=1 tc_ens_xz_cbgt.ncl]

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "$LIB/gsn_csm.ncl"

begin

;================================================
; Read in cyclone position data from text file
;================================================

   diri        = "/nfs/a37/earshar/ncl_scripts/nepartak/archer_march2018/"
   numTIMES    = 120
   centre      = new((/2,numTIMES,2/),float) ; Array to hold all [120] lat/lon values

 ; Read in values from 'em07'
   centre(0,:,0) = asciiread(diri+"tc_lat_em07.txt",(/numTIMES,1/),"float")
   centre(0,:,1) = asciiread(diri+"tc_lon_em07.txt",(/numTIMES,1/),"float")

 ; Read in values from 'em05'
   centre(1,:,0) = asciiread(diri+"tc_lat_em05.txt",(/numTIMES,1/),"float")
   centre(1,:,1) = asciiread(diri+"tc_lon_em05.txt",(/numTIMES,1/),"float")

;============================================================================
; Calculate distance between two points --> then calculate system velocity
;============================================================================

  print_clock("Starting system velocity calculations...")

; Define new arrays (cyclone position, speed)
  x_cyc   = new((/2,numTIMES/),float)
  y_cyc   = new((/2,numTIMES/),float)
  u_cyc   = new((/2,numTIMES/),float)
  v_cyc   = new((/2,numTIMES/),float)
  vel_cyc = new((/2,numTIMES/),float)

 do pt = 0, 1

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
   y_cyc(pt,0)	 = 0
   u_cyc(pt,0)   = u_cyc(pt,1)
   v_cyc(pt,0)   = v_cyc(pt,1)
   vel_cyc(pt,0) = vel_cyc(pt,1)

 ; Add metadata to new arrays
   u_cyc@units         = "m s~S~-1~N~"
   u_cyc@description   = "Cyclone zonal velocity"

   v_cyc@units         = "m s~S~-1~N~"
   v_cyc@description   = "Cyclone meridional velocity"

   vel_cyc@units       = "m s~S~-1~N~"
   vel_cyc@description = "Cyclone speed"

   u_cyc(pt,:) = -(u_cyc(pt,:)) ; Easterly motion is negative

  end do      ; End of ensemble loop (do pt = 0, 1)

;================================================================================
; Create lat/lon arrays containing grid subsets at each time, following the TC
;================================================================================
  lat_arr = new((/2,numTIMES,2/),float)
  lon_arr = new((/2,numTIMES,2/),float)

;  dist    = 3.6        ; Distance in degrees (whole number of grid points)

; Lat/lon arrays for 'em07'
  lat_arr(0,:,0) = centre(0,:,0)-dist
  lat_arr(0,:,1) = centre(0,:,0)+dist
  lon_arr(0,:,0) = centre(0,:,1)-dist
  lon_arr(0,:,1) = centre(0,:,1)+dist

; Lat/lon arrays for 'em05'
  lat_arr(1,:,0) = centre(1,:,0)-dist
  lat_arr(1,:,1) = centre(1,:,0)+dist
  lon_arr(1,:,0) = centre(1,:,1)-dist
  lon_arr(1,:,1) = centre(1,:,1)+dist

;==================================================================================
;;;;;;;;;;;;;;;;;;;;;;;;;;; Start of main part of script ;;;;;;;;;;;;;;;;;;;;;;;;;;
;==================================================================================

  print_clock("Starting main part of script...")

; Prelash [define arrays for use in loops below]

  time_arr    = new((/2,numTIMES/),string)      ; Array to hold strings for second loop
  title_arr   = new((/2,numTIMES/),string)      ; Array to hold strings for plot titles
  pos_arr     = new((/2,numTIMES,2/),float)     ; Array for cyclone position [lat/lon]
  r_earth     = 6.37e6                          ; Earth's radius (constant)
  pt	      = 0				; Counter variable (ensemble member)
  plevs	      = 19				; Pressure levels

  date_str = new(numTIMES,string)
  time_str = new(numTIMES,string)
  out_str  = new(numTIMES,string)

; Create arrays to hold circulation tendency values 
; [ensemble member (2) * box ensemble (25) * time (numTIMES) * levs (plevs)]
  circ_inst = new((/2,9,numTIMES,plevs/),float)
  hadv_ten  = new((/2,9,numTIMES,plevs/),float)
  str_ten   = new((/2,9,numTIMES,plevs/),float)
  flux_ten  = new((/2,9,numTIMES,plevs/),float)
  tilt_ten  = new((/2,9,numTIMES,plevs/),float)
  fric_ten  = new((/2,9,numTIMES,plevs/),float)

;==================================================
; Start loop over ensemble members [em01, em05]
;==================================================

  ens_arr = (/"em07","em05"/)

  do en = 0, 1

   print("Working on ensemble member: "+ens_arr(en))

 ; List all files to read in and analyse
   diri        = "$ens/reg_netcdf/20160704T0000Z_ra1t_"
   input_pb    = diri+ens_arr(en)+"_pb.nc"
   input_pc    = diri+ens_arr(en)+"_pc.nc"
   input_pf    = diri+ens_arr(en)+"_pf.nc"

   setvalues NhlGetWorkspaceObjectId
     "wsMaximumSize" : 1000000000
   end setvalues

;===============================
; Start loop over input files 
;===============================

   nf = 0
   f  = addfile(input_pb(nf),"r")	; Read in 'pb' stream [u,v,w,vort,div]
   g  = addfile(input_pc(nf),"r")	; Read in 'pc' stream [friction]
   h  = addfile(input_pf(nf),"r")	; Read in 'pf' stream [density on model levs]

;==============================================
; Get the variables we will need [pb stream]
;==============================================

 do it = 0, numTIMES-1   ; Loop over all times in file

   time  = f->t(it)        ; Time in file
   times = dimsizes(time)

   t0    = lat_arr(en,it,0)
   t1    = lat_arr(en,it,1)
   n0    = lon_arr(en,it,0)
   n1    = lon_arr(en,it,1)
   print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

   lon1  = f->longitude_1({n0:n1})  ; longitude [160 grid points]
   lat1  = f->latitude_1({t0:t1})   ; latitude  [160 grid points]
   lon   = f->longitude({n0:n1})    ; longitude [160 grid points]
   lat   = f->latitude({t0:t1})     ; latitude  [160 grid points]

   pres  = f->p            ; Pressure (levels)
   plevs = dimsizes(pres)  ; Size of pressure level array
   pres@units = "hPa"

   u     = f->u(it,:,{t0:t1},{n0:n1})     ; Zonal velocity (m/s)
   v     = f->v(it,:,{t0:t1},{n0:n1})     ; Meridional velocity (m/s)
   vort  = f->rvor(it,:,{t0:t1},{n0:n1})  ; Relative vorticity (/s)
   div   = f->div(it,:,{t0:t1},{n0:n1})   ; Divergence (/s)

   t     = f->temp(it,:,{t0:t1},{n0:n1})  ; Temperature (K)
   z     = f->ht(it,:,{t0:t1},{n0:n1})    ; Geopotential height (m)
   w     = f->dz_dt(it,:,{t0:t1},{n0:n1}) ; Vertical velocity (m/s)

 ; Convert to units we want
   z     = z/10
   z@units = "dam" ; Convert to decametres

   vort  = vort * (10 ^ -6) ; Convert to /s
   div   = div  * (10 ^ -6) ; Convert to /s

   dx    = lon(1) - lon(0) ; Grid spacing (longitude)
   dy    = lat(1) - lat(0) ; Grid spacing (latitude)

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

 ; Pressure levels (pres)
 ; (0) 1000, (1) 950, (2) 900, (3) 850, (4) 800, (5) 750, (6) 700, (7) 650
 ; (8) 600, (9) 550, (10) 500, (11) 450, (12) 400, (13) 350
 ; (14) 300, (15) 250, (16) 200, (17) 150, (18) 100

;==============================================
; Get the variables we will need [pc stream]
;==============================================

 ; Get additional variables from the model-level file
   lon1F = g->longitude_1({n0:n1})       ; longitude [160 grid points]
   lat1F = g->latitude_1({t0:t1})        ; latitude  [160 grid points]
   lonF  = g->longitude({n0:n1})         ; longitude [160 grid points]
   latF  = g->latitude({t0:t1})          ; latitude  [160 grid points]
   hyb   = g->hybrid_ht                  ; Array of 63 levels, from 5 m up to ~ 20 km
   strx  = g->taux(it,:,{t0:t1},{n0:n1}) ; Wind stress, x-comp (N m/2)
   stry  = g->tauy(it,:,{t0:t1},{n0:n1}) ; Wind stress, y-comp (N m/2)

   dxF   = lonF(1) - lonF(0) ; Grid spacing (longitude)
   dyF   = latF(1) - latF(0) ; Grid spacing (latitude)   

;=========================================================================
; Get density from pf stream [same size as 'strx', smaller than 'stry']
;=========================================================================

   rho1  = h->field27(it,:,{t0:t1},{n0:n1}) ; Dry air density (kg/m-3)

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

   date_str(it) = sprinti("%0.2iUTC ", hour) + \
                  sprinti("%0.2i ", day) \
                  + month_abbr(month)
   time_str(it) = sprinti("%0.2i UTC ", hour) + \
                  sprinti("%0.2i ", day) \
                  + month_abbr(month)
   out_str(it)  = sprinti("%0.2i", day) + \
                  month_abbr(month) + \
                  "_" + sprinti("%0.2iZ", hour)

   time_arr(en,it) = out_str(it)

   print_clock("Working on time: "+time_str(it))
   title_arr(en,it) = "Valid at: "+time_str(it)

;==========================================================================
; Define new arrays for terms in vorticity equation (no friction for now)
;==========================================================================

 ; Storm-relative velocity
   u_rel      = u - u_cyc(en,it)
   v_rel      = v - v_cyc(en,it)

   d          = dimsizes(vort)       ; Size of all arrays [19 * 180 * 180]
   d0	      = d(0)-1
   d1	      = d(1)-1
   d2	      = d(2)-1
   pii        = 3.14159265           ; Constant
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
    coslat  = cos(y2d*pii/180)                      ; Same array in radians

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

;=================================
; Friction-related calculations 
;=================================

  ; Array of hybrid levels to approximately match pressure levels every 50 hPa
  ; Array [/1000, 950, 900, 850, 800, 750, etc, ..., 100/]
    lev_arr = (/1, 8, 13, 16, 19, 21, 23, 25, 27, 29, \
                31, 33, 35, 37, 39, 42, 47, 51, 56/)
    lev_arr@units = "m"

    l0    = 0  ; Subset of vertical levels (19 levels)
    l1    = 18 ; Level 56 of original array ~ 100 hPa (Int. Standard Atmosphere)

  ; Array of pressure levels with the same dimensions as 'lev_arr'
    pres       = (/1000,950,900,850,800,750,700,650,600,550,\
                   500,450,400,350,300,250,200,150,100/)
    pres@units = "hPa"
    plevs      = dimsizes(pres)

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
    tilt_y = dv_dp(:,0:d1,0:d2) * w(:,0:d1,0:d2)    ; X-component of relative vorticity

  ; Calculate [-du_dp * w]
    tilt_x = -du_dp(:,0:d1,0:d2) * w(:,0:d1,0:d2)   ; Y-component of relative vorticity

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
     ilev = 1

   ; First, concentrate on a subset of the global grid 
     geo_plane     = z(ilev,0:d1,0:d2)
     vort_plane    = vort(ilev,0:d1,0:d2)

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
     if (it.eq.0) then
      centre(en,it,0) = lat_min950
      centre(en,it,1) = lon_min950

     else
      lt0   = centre(en,it-1,0)         ; Previous TC latitude
      ln0   = centre(en,it-1,1)         ; Previous TC longitude
      print("lt0 = "+lt0+" ; ln0 = "+ln0)
      d_lat = abs(lat_min950-lt0)       ; Change in latitude (t1 - t0)
      d_lon = abs(lon_min950-ln0)       ; Change in longitude (note extra minus sign)
      print("d_lat = "+d_lat+" ; d_lon = "+d_lon)

      if (d_lat.gt.1.or.d_lon.gt.1) then        ; Recalculate TC centre if incorrect

     ; Create smaller grid [0.5 degrees]
       lt1 = lt0 - 0.5
       lt2 = lt0 + 0.5
       ln1 = ln0 - 0.5
       ln2 = ln0 + 0.5

       geo_plane0 = z(ilev,{lt1:lt2},{ln1:ln2})          ; Smaller grid (1.0 x 1.0 deg)
       geo_min0   = min(geo_plane0)                      ; Find minimum on smaller grid
       dims_h0    = dimsizes(geo_plane0)                 ; Size of smaller grid
       geo_1d0    = ndtooned(geo_plane0)                 ; Create 1-D array
       inds_h0    = ind_resolve(minind(geo_1d0),dims_h0) ; Find index of minimum

       lat_min950 = lt1 + (dy * inds_h0(0,0))
       lon_min950 = ln1 + (dy * inds_h0(0,1))
       print("New centre (hgt): "+lat_min950+" degrees N, "+lon_min950+" degrees E")
       centre(en,it,0) = lat_min950
       centre(en,it,1) = lon_min950

       delete([/geo_plane0, geo_min0, dims_h0, geo_1d0, inds_h0/])

      else
       print("Cyclone centre (hgt): "+lat_min950+" degrees N, "+lon_min950+" degrees E")
       centre(en,it,0) = lat_min950
       centre(en,it,1) = lon_min950
      end if

     end if

   ; Get information about box size for later ensemble calculations
     r = rad   ; Box radius (degrees)
     box = ddx({lat_min950-r:lat_min950+r},{lon_min950-r:lon_min950+r})
     box_size = dimsizes(box)
     size0    = box_size(0)
     size1    = box_size(1)

     s0	      = size0-1
     s1	      = size1-1

   ; Tidy up before looping over pressure levels
     delete([/vort1d, inds, vort1d_sm, vort_smth, inds_sm/])
     delete([/geo1d, inds_h, geo1d_sm, geo_smth, inds_h_sm/])
     delete([/vort_plane, geo_plane/])

;============================================================
; Loop over pressure levels to calculate vertical profile
;============================================================

    cp = 0               ; Counter variable

    do ilev = 0, plevs-1 ; 1000,950,900,850,800,750,700,650,600,550,...,100

   ; Pressure levels ['pb' stream]
     p = pres(ilev)
     print("Working on: " + p + " hPa")

   ;========================================
   ; Model levels ['pc' and 'pf' streams]
   ;========================================

     hy = sprintf("%0.0f",hyb(ilev))
     print("Working on hybrid level " + ilev+ " ("+hy+" m AGL)")
     m0 = lev_arr(ilev)	      ; Selected model levels, chosen to match pressure levels

   ; Calculate components of friction term
     ip1 = min((/m0+1,lev_arr(l1)/)) ; Account for top and bottom levels
     im1 = max((/m0-1,lev_arr(l0)/))

   ; Current level [wind stress]
     taux   = strx(m0,0:d1,0:d2)
     tauy   = stry(m0,0:d1,0:d2)

   ; Density ['pf' stream]
     rho    = rho1(m0,0:d1,0:d2)
     rho_a = (1/rho) ; Specific volume (inverse density)

   ; Levels above/below for vertical derivatives [wind stress]
     tau_x1 = strx(ip1,0:d1,0:d2)
     tau_x0 = strx(im1,0:d1,0:d2)
     tau_y1 = stry(ip1,0:d1,0:d2)
     tau_y0 = stry(im1,0:d1,0:d2)

   ; Vertical derivative of 'strx'
     dtaux_dz = (tau_x1 - tau_x0) / (hyb(ip1) - hyb(im1))

   ; Vertical derivative of 'stry'
     dtauy_dz = (tau_y1 - tau_y0) / (hyb(ip1) - hyb(im1))

   ; X-component of friction
     fric_x_plane           = rho_a(:,:) * dtaux_dz
     copy_VarCoords(taux,fric_x_plane)
     copy_VarMeta(taux,fric_x_plane)
     fric_x_plane@units     = "m s~S~-2~N~"
     fric_x_plane@title     = "x-component of friction"
     fric_x_plane@long_name = "x-component of friction"
     fric_x_plane@name      = "fric_x"     

   ; Y-component of friction
     fric_y_plane           = rho_a(:,:) * dtauy_dz
     copy_VarCoords(tauy,fric_y_plane)
     copy_VarMeta(tauy,fric_y_plane)
     fric_y_plane@units     = "m s~S~-2~N~"
     fric_y_plane@title     = "y-component of friction"
     fric_y_plane@long_name = "y-component of friction"
     fric_y_plane@name      = "fric_y"

   ;=================================
   ; Pressure levels ['pb' stream]
   ;=================================

   ; Storm-relative velocity
     u_plane       = u_rel(ilev,0:d1,0:d2)
     v_plane       = v_rel(ilev,0:d1,0:d2)

   ; Other variables
     vort_plane    = vort(ilev,0:d1,0:d2)
     avo_plane	   = avo(ilev,0:d1,0:d2)
     div_plane	   = div(ilev,0:d1,0:d2)
     geo_plane     = z(ilev,0:d1,0:d2)
     t_plane       = t(ilev,0:d1,0:d2)

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

     delete([/vort1d, inds, vort1d_sm, vort_smth, inds_sm/])

   ; Account for cases in which another maximum is chosen in error 
     lat_diff = abs(lat_max-lat_min950)
     lon_diff = abs(lon_max-lon_min950)

     if (lon_diff.gt.0.5.or.lat_diff.gt.0.5) then
        lat_max = lat_min950
        lon_max = lon_min950
     end if
     print("Cyclone centre: "+lat_max+" degrees N, "+lon_max+" degrees E")

   ; Add lat/lon values to array for later calculations
     pos_arr(en,it,0) = lat_max
     pos_arr(en,it,1) = lon_max

   ; Terms in the circulation budget
   ; 'tilt_x' and 'tilt_y' contain values on all pressure levels for each time
   ; 'tilt_x_plane' and 'tilt_y_plane' contain values on single pressure levels
     tilt_x_plane               = tilt_x(ilev,:,:)       ; Tilting (x-component)
     tilt_y_plane               = tilt_y(ilev,:,:)       ; Tilting (y-component)

;=====================================================================================
; Create array of box positions either side of the original [Davis & Galarneau 2009]
;=====================================================================================

;     lat_box = (/lat_max-(2*dy),lat_max-(dy),lat_max,\
;     	         lat_max+(dy),lat_max+(2*dy)/)
;     lon_box = (/lon_max-(2*dx),lon_max-(dx),lon_max,\
;                 lon_max+(dx),lon_max+(2*dx)/)

     lat_box = (/lat_max-(dy),lat_max,lat_max+(dy)/)
     lon_box = (/lon_max-(dx),lon_max,lon_max+(dx)/)

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
     fric_x_box = new((/lsize,size0,size1/),float)
     fric_y_box = new((/lsize,size0,size1/),float)

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

    ; Size of original array (x-component)
      dim1 = dimsizes(fric_x_plane({lat_box(j)-r:lat_box(j)+r},\
                                       {lon_box(i)-r:lon_box(i)+r}))

    ; Size of original array (y-component)
      dim2 = dimsizes(fric_y_plane({lat_box(j)-r:lat_box(j)+r},\
                                       {lon_box(i)-r:lon_box(i)+r}))

    ; Size of new, larger array (x-component)
      dim3 = dimsizes(fric_x_plane({lat_box(j)-r-dy:lat_box(j)+r},\
                                         {lon_box(i)-r:lon_box(i)+r}))

    ; Size of new, smaller array (x-component)
      dim4 = dimsizes(fric_x_plane({lat_box(j)-r+dy:lat_box(j)+r},\
                                         {lon_box(i)-r:lon_box(i)+r}))

      printVarSummary(fric_x_box(a,:,:))
      printVarSummary((fric_x_plane({lat_box(j)-r:lat_box(j)+r},\
                                       {lon_box(i)-r:lon_box(i)+r})))
      printVarSummary(fric_y_box(a,:,:))
      printVarSummary((fric_y_plane({lat_box(j)-r:lat_box(j)+r},\
                                       {lon_box(i)-r:lon_box(i)+r})))

      print("box_size = "+box_size)
      print(dimsizes(fric_x_plane({lat_box(j)-r:lat_box(j)+r},\
                                         {lon_box(i)-r:lon_box(i)+r})))
      exit()
      test = conform_dims(box_size,fric_x_plane({lat_box(j)-r:lat_box(j)+r},\
                                         {lon_box(i)-r:lon_box(i)+r}),1)
      printVarSummary(test)
      exit()

;      print("Size of LHS array: "+dimsizes(fric_x_box(a,:,:)))
;      print("Size of RHS array: "+dim1)
;      print("Size of larger RHS array: "+dim3)
;      print("Size of smaller RHS array: "+dim4)	
;      exit()

    ; Make sure all arrays are the same size
      if (size0.gt.dim1(0)) then
        print("LHS latitude array is bigger than RHS!")
        fric_x_box(a,:,:) = fric_x_plane({lat_box(j)-r-dxy:lat_box(j)+r+dxy},\
 				         {lon_box(i)-r:lon_box(i)+r})       
      else if (size0.lt.dim1(0)) then
        print("LHS latitude array is smaller than RHS!")
        fric_x_box(a,:,:) = fric_x_plane({lat_box(j)-r+dxy:lat_box(j)+r-dxy},\
                                         {lon_box(i)-r:lon_box(i)+r})
      else if (size1.gt.dim1(1)) then
        print("LHS longitude array is bigger than RHS!")
        fric_x_box(a,:,:) = fric_x_plane({lat_box(j)-r:lat_box(j)+r},\
                                         {lon_box(i)-r-dxy:lon_box(i)+r+dxy})      
      else if (size1.lt.dim1(1)) then
        print("LHS longitude array is smaller than RHS!")
        fric_x_box(a,:,:) = fric_x_plane({lat_box(j)-r:lat_box(j)+r},\
                                         {lon_box(i)-r+dxy:lon_box(i)+r-dxy})
      else
        fric_x_box(a,:,:) = fric_x_plane({lat_box(j)-r:lat_box(j)+r},\
                                         {lon_box(i)-r:lon_box(i)+r})
      end if
      end if
      end if
      end if

      if (size0.gt.dim2(0)) then
        fric_y_box(a,:,:) = fric_y_plane({lat_box(j)-r-dxy:lat_box(j)+r+dxy},\
                                         {lon_box(i)-r:lon_box(i)+r})
      else if (size0.lt.dim2(0)) then
        fric_y_box(a,:,:) = fric_y_plane({lat_box(j)-r+dxy:lat_box(j)+r-dxy},\
                                         {lon_box(i)-r:lon_box(i)+r})
      else if (size1.gt.dim2(1)) then
        fric_y_box(a,:,:) = fric_y_plane({lat_box(j)-r:lat_box(j)+r},\
                                         {lon_box(i)-r-dxy:lon_box(i)+r+dxy})      
      else if (size1.lt.dim2(1)) then
        fric_y_box(a,:,:) = fric_y_plane({lat_box(j)-r:lat_box(j)+r},\
                                         {lon_box(i)-r+dxy:lon_box(i)+r-dxy})
      else
        fric_y_box(a,:,:) = fric_y_plane({lat_box(j)-r:lat_box(j)+r},\
	 				 {lon_box(i)-r:lon_box(i)+r})
      end if
      end if
      end if
      end if

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
     fric_x0     = dim_sum_n_Wrap((fric_x_box(a,:,:) * ddx_box(a,:,:)),1)
     tilt_s      = tilt_x0(0)              ; Tilting (south)
     tilt_n      = tilt_x0(numPOINTS_y-1)  ; Tilting (north)
     flux_s      = flux_x0(0)              ; Horizontal advection (south)
     flux_n      = -flux_x0(numPOINTS_y-1) ; Horizontal advection (north)
     hadv_s      = hadv_x0(0)              ; Horizontal advection (south)
     hadv_n      = -hadv_x0(numPOINTS_y-1) ; Horizontal advection (north)
     circ_s      = circ_x0(0)              ; Circulation (south)
     circ_n      = -circ_x0(numPOINTS_y-1) ; Circulation (north)
     fric_s      = fric_x0(0)              ; Friction (south)
     fric_n      = -fric_x0(numPOINTS_y-1) ; Friction (north)

   ; Calculate terms on eastern/western boundaries
     tilt_y0    = dim_sum_n_Wrap((tilt_y_box(a,:,:) * ddy_box(a,:,:)),0)
     flux_y0    = dim_sum_n_Wrap((-avo_box(a,:,:) * u_box(a,:,:) * \
                                  ddy_box(a,:,:)),0)
     hadv_y0    = dim_sum_n_Wrap((-(avo_box(a,:,:) - avo_avg(a)) * \
     		  		 (u_box(a,:,:) - u_avg(a)) * \
				  ddy_box(a,:,:)),0)
     circ_y0    = dim_sum_n_Wrap((v_box(a,:,:) * ddy_box(a,:,:)),0)
     fric_y0    = dim_sum_n_Wrap((fric_y_box(a,:,:) * ddy_box(a,:,:)),0)
     tilt_w     = tilt_y0(0)              ; Tilting (west)
     tilt_e     = tilt_y0(numPOINTS_x-1)  ; Tilting (east)
     flux_w     = flux_y0(0)              ; Advective flux divergence (west)
     flux_e     = -flux_y0(numPOINTS_x-1) ; Advective flux divergence (east)
     hadv_w     = hadv_y0(0)              ; Horizontal advection (west)
     hadv_e     = -hadv_y0(numPOINTS_x-1) ; Horizontal advection (east)
     circ_w     = -circ_y0(0)             ; Circulation (west)
     circ_e     = circ_y0(numPOINTS_x-1)  ; Circulation (east)
     fric_w     = -fric_y0(0)             ; Friction (west)
     fric_e     = fric_y0(numPOINTS_x-1)  ; Friction (east)

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
     circ_instT = circ_s + circ_n + circ_w + circ_e	     ; Sum 4 line integrals
     circ_instT = circ_instT				     ; Do not divide by area

   ; Calculate horizontal advection from line integrals (4 sides of box)
     fric_tenT  = fric_s + fric_n + fric_w + fric_e
     fric_tenT  = fric_tenT / d_ave(a) ; Divide by box area [Rios-Berrios et al. 2016]

   ; Enter values into predefined arrays
     hadv_ten(en,a,it,cp)      = hadv_tenT		     ; Horizontal advection
     str_ten(en,a,it,cp)       = str_tenT		     ; Stretching
;     flux_ten(en,a,it,cp)      = flux_tenT		     ; Advective flux
     tilt_ten(en,a,it,cp)      = tilt_tenT		     ; Tilting
     circ_inst(en,a,it,cp)     = circ_instT                  ; Circulation (LHS)
     fric_ten(en,a,it,cp)      = fric_tenT		     ; Friction 

   ; Tidy up (delete variables from within lat/lon loops before next iteration)
;     delete([/flux_tenT/])
     delete([/hadv_tenT, str_tenT, tilt_tenT, circ_instT, fric_tenT/])
     delete([/vort_x_final0, vort_x_final1, vort_y_final0, vort_y_final1/])
     delete([/avo_x_final0,  avo_x_final1,  avo_y_final0,  avo_y_final1/])
     delete([/uy_final0,     uy_final1,     ux_final0,     ux_final1/])
     delete([/vy_final0,     vy_final1,     vx_final0,     vx_final1/])
     delete([/numPOINTS_x0, numPOINTS_x, numPOINTS_y/])

     end do     ; End latitude loop (do j = 0, dimsizes(lat_box)-1)
    end do	; End longitude loop (do i = 0, dimsizes(lon_box)-1)

    cp                   = cp + 1       ; Counter variable (pressure)

  ; Tidy up (delete variables from within pressure level loop)
    delete([/u_plane, v_plane, tilt_x_plane, tilt_y_plane/])
    delete([/fric_x_plane, fric_y_plane, dtaux_dz, dtauy_dz/])
    delete([/rho,rho_a/])
    delete([/vort_plane, avo_plane, div_plane/])
    delete([/ddx_box, ddy_box, u_box, v_box, fric_x_box, fric_y_box/])
    delete([/tilt_x_box, tilt_y_box, vort_box, avo_box, div_box, lat_box, lon_box/])
    delete([/div_ave,vort_ave,vort_avg,avo_avg,u_avg,v_avg/])

   end do ; End pressure level loop (do ilev = 0, plevs-1)

  end do                ; End time loop (do it = 0, times-1)

; Tidy up (delete variables from within input file loop before next iteration)
  delete([/lon,lat,lon1,lat1,u,v,vort,div,f0,avo,t,z,w/])
  delete([/time,times,utc_date,year,month,day,hour,minute,second/])
  delete([/lon1F,lat1F,lonF,latF,hyb,strx,stry,dxF,dyF,rho1/])

 end do	      ; End ensemble member loop (do en = 0, 1)

;=====================================================
; Calculate term averages over all ensemble members
;=====================================================

   print_clock("Finished ensemble calculations! Now averaging over all members...")
   exit()

   hadv_final = new((/2,numTIMES,plevs/),"float")
   str_final  = new((/2,numTIMES,plevs/),"float")
;   flux_final = new((/2,numTIMES,plevs/),"float")
   tilt_final = new((/2,numTIMES,plevs/),"float")
   circ_final = new((/2,numTIMES,plevs/),"float")
   fric_final = new((/2,numTIMES,plevs/),"float")

   circ_inst  = circ_inst / 1000000 ; Divide by 1000000 [following 'tc_ens_xz_circ']

   do en = 0, 1    ; Loop over ensemble members
 
    hadv_final(en,:,:) = dim_avg_n_Wrap(hadv_ten(en,:,:,:),0)
    str_final(en,:,:)  = dim_avg_n_Wrap(str_ten(en,:,:,:),0)
;    flux_final(en,:,:) = dim_avg_n_Wrap(flux_ten(en,:,:,:),0)
    tilt_final(en,:,:) = dim_avg_n_Wrap(tilt_ten(en,:,:,:),0)
    circ_final(en,:,:) = dim_avg_n_Wrap(circ_inst(en,:,:,:),0)
    fric_final(en,:,:) = dim_avg_n_Wrap(fric_ten(en,:,:,:),0)

   end do	   ; End loop over ensemble members

;=========================================
; Integrate over Nepartak's life cycle
;=========================================

 ; Each calculation produces an array of 19 values [pressure levels]
 ; We do this calculation [numTIMES] times
 ; Create arrays with [numTIMES] * 19 elements

   hadv_intA = new((/2,numTIMES,plevs/),float)
   str_intA  = new((/2,numTIMES,plevs/),float)
   flux_intA = new((/2,numTIMES,plevs/),float)
   tilt_intA = new((/2,numTIMES,plevs/),float)
   tot_intA  = new((/2,numTIMES,plevs/),float)
   fric_intA = new((/2,numTIMES,plevs/),float)

   csize     = dimsizes(hadv_intA)
   csize0    = csize(1)

   time_int = new(csize0,"integer")
   time_int = ispan(0,csize0-1,1)
   time_int!0 = "time"

   hadv_intA@units = "10~S~-5~N~ s~S~-1~N~"
   str_intA@units  = "10~S~-5~N~ s~S~-1~N~"
   flux_intA@units = "10~S~-5~N~ s~S~-1~N~"
   tilt_intA@units = "10~S~-5~N~ s~S~-1~N~"
   tot_intA@units  = "10~S~-5~N~ s~S~-1~N~"
   fric_intA@units = "10~S~-5~N~ s~S~-1~N~"

   print_clock("Starting plots...")

   do en = 0, 1	         ; Loop over ensemble members [em01, em05]
    do ic = 0, csize0-1  ; 1 h integrations over the whole of Nepartak's life cycle

     if (ic.lt.csize0-1) then
      cs = ic            ; Start index
      cf = cs + 1        ; Finish index
     else
      cs = ic 
      cf = ic
     end if

      f_str = sprinti("%0.3i",cs)+"h"    ; Formatted string for output file

      hadv_intA(en,ic,:) = dim_sum_n_Wrap((hadv_final(en,cs:cf,:) * 3600), 0)
      str_intA(en,ic,:)  = dim_sum_n_Wrap((str_final(en,cs:cf,:) * 3600), 0)
;      flux_intA(en,ic,:) = dim_sum_n_Wrap((flux_final(en,cs:cf,:) * 3600), 0)
      tilt_intA(en,ic,:) = dim_sum_n_Wrap((tilt_final(en,cs:cf,:) * 3600), 0)
      fric_intA(en,ic,:) = dim_sum_n_Wrap((fric_final(en,cs:cf,:) * 3600), 0)

      hadv_intA(en,ic,:) = hadv_intA(en,ic,:) * (10 ^ 5)
      str_intA(en,ic,:)  = str_intA(en,ic,:) * (10 ^ 5)
;      flux_intA(en,ic,:) = flux_intA(en,ic,:) * (10 ^ 5)
      tilt_intA(en,ic,:) = tilt_intA(en,ic,:) * (10 ^ 5)
      fric_intA(en,ic,:) = fric_intA(en,ic,:) * (10 ^ 5)

    end do               ; End time loop (do ic = 0, csize0-1)
   end do		 ; End ensemble member loop (do en = 0, 1)

   flux_intA = hadv_intA + str_intA
   tot_intA  = flux_intA + tilt_intA + fric_intA

 ; Add metadata before time-height plotting
   tot_intA!0	  = "ens"
   tot_intA&ens   = ens_arr(p1:p2:int)
   tot_intA!1     = "time"
   tot_intA&time  = time_int
   tot_intA!2     = "lev"
   tot_intA&lev   = pres

   hadv_intA!0	  = "ens"
   hadv_intA&ens  = ens_arr(p1:p2:int)
   hadv_intA!1    = "time"
   hadv_intA&time = time_int
   hadv_intA!2    = "lev"
   hadv_intA&lev  = pres

   str_intA!0	  = "ens"
   str_intA&ens   = ens_arr(p1:p2:int)
   str_intA!1     = "time"
   str_intA&time  = time_int
   str_intA!2     = "lev"
   str_intA&lev   = pres

   flux_intA!0    = "ens"
   flux_intA&ens  = ens_arr(p1:p2:int)
   flux_intA!1    = "time"
   flux_intA&time = time_int
   flux_intA!2    = "lev"
   flux_intA&lev  = pres

   tilt_intA!0	  = "ens"
   tilt_intA&ens  = ens_arr(p1:p2:int)
   tilt_intA!1    = "time"
   tilt_intA&time = time_int
   tilt_intA!2    = "lev"
   tilt_intA&lev  = pres

   fric_intA!0    = "ens"
   fric_intA&ens  = ens_arr(p1:p2:int)
   fric_intA!1    = "time"
   fric_intA&time = time_int
   fric_intA!2    = "lev"
   fric_intA&lev  = pres	;;; CHANGE TO MODEL LEVELS (?) ;;;

;======================================
; Reorder dimensions before plotting
;======================================

 ; Define new arrays for final plotting 
   circ_int = new((/plevs,2,numTIMES/),float)
   hadv_int = new((/plevs,2,numTIMES/),float)
   str_int  = new((/plevs,2,numTIMES/),float)
   flux_int = new((/plevs,2,numTIMES/),float)
   tilt_int = new((/plevs,2,numTIMES/),float)
   tot_int  = new((/plevs,2,numTIMES/),float)
   fric_int = new((/plevs,2,numTIMES/),float)

 ; Add metadata to circulation array
   circ_final!0     = "ens"
   circ_final&ens   = ens_arr(p1:p2:int)
   circ_final!1     = "time"
   circ_final&time  = time_int
   circ_final!2     = "lev"
   circ_final&lev   = pres

 ; Reorder from [time * pressure] to [pressure * time]
   do en = 0, 1   ; Loop over ensemble members

    circ_int(:,en,:) = circ_final(lev|:,ens|en,time|:)
    hadv_int(:,en,:) = hadv_intA(lev|:,ens|en,time|:)
    str_int(:,en,:)  = str_intA(lev|:,ens|en,time|:)
    flux_int(:,en,:) = flux_intA(lev|:,ens|en,time|:)
    tilt_int(:,en,:) = tilt_intA(lev|:,ens|en,time|:)
    tot_int(:,en,:)  = tot_intA(lev|:,ens|en,time|:)
    fric_int(:,en,:) = fric_intA(lev|:,ens|en,time|:)

   end do	  ; End ensemble member loop (do en = 0, 1)

;==============================================
; Calculate running average for budget terms
;==============================================

 ; Running average over time dimension (2)
 ; User input: only calculates average when "run" = 0 when calling script
 ; User input: "nave" [number of points used to calculate running average]

   if (run.eq.0) then
    tot_int  = runave_n_Wrap(tot_int, nave, 0, 2)
    hadv_int = runave_n_Wrap(hadv_int, nave, 0, 2)
    str_int  = runave_n_Wrap(str_int, nave, 0, 2)
    tilt_int = runave_n_Wrap(tilt_int, nave, 0, 2)
    fric_int = runave_n_Wrap(fric_int, nave, 0, 2)
    output = "$sam/nepartak/images/cart_coords/ens9_"+r1+"_xz_cbgtF_paperF_"+nave+"ave"
   else
    output = "$sam/nepartak/images/cart_coords/ens9_"+r1+"_xzF_cbgtF_paperF_intF_"+\
    	       ens_arr(p1)+"_"+ens_arr(p2)
   end if
 
;====================================================
; Open workstation for plotting (time-height plots)
;====================================================

     wks = gsn_open_wks(opt,output)

;======================================
; Create an x-y time series plot
;======================================

     panel = new(8,graphic)
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
     opts_bdgt@cnLevels                    = (/-1.0, -0.5, -0.2, -0.1, -0.05, \
                                                0.05, 0.1, 0.2, 0.5, 1.0/)
     opts_bdgt@cnFillColors                = (/3,4,6,8,10,0,12,14,16,18,19/)
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
     opts_bdgt@tmXBValues                  = (/0,23,47,71/)
     opts_bdgt@tmXBLabels                  = (/"00Z 5 Jul","00Z 6 Jul",\
     					       "00Z 7 Jul","00Z 8 Jul"/)
     opts_bdgt@tmXBLabelFontHeightF        = 0.025
     opts_bdgt@tmXBMinorValues             = ispan(2,78,3) ; Minor tickmarks (every 3h)
     opts_bdgt@trYMaxF                     = 1000          ; Ignore BL
     opts_bdgt@trYMinF                     = 150           ; Include upper troposphere

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
     opts_cont@cnLevels                     = (/24.0, 28.0, 32.0, 36.0, 40.0/)
;    opts_cont@cnLevels			    = (/14.0, 16.0, 18.0, 20.0/)
     opts_cont@cnLineThicknessF             = 2             ; Contour line thickness
     opts_cont@cnLineLabelsOn               = False         ; Turn off contour labels
     opts_cont@gsnContourNegLineDashPattern = 2             ; Negative contours dashed
     opts_cont@gsnContourPosLineDashPattern = 0             ; Positive contours solid
     opts_cont@cnInfoLabelOn                = False         ; Turn off info label
     circ_int@units                         = ""            ; Turn off string above plot

;=================
; Plot the data
;=================

     plot_tot06   = gsn_csm_pres_hgt(wks,tot_int(:,0,:),opts_bdgt)
     plot_flux06  = gsn_csm_pres_hgt(wks,flux_int(:,0,:),opts_bdgt)
     plot_tilt06  = gsn_csm_pres_hgt(wks,tilt_int(:,0,:),opts_bdgt)
     plot_fric06  = gsn_csm_pres_hgt(wks,fric_int(:,0,:),opts_bdgt)

     plot_tot07	  = gsn_csm_pres_hgt(wks,tot_int(:,1,:),opts_bdgt)
     plot_flux07  = gsn_csm_pres_hgt(wks,flux_int(:,1,:),opts_bdgt)
     plot_tilt07  = gsn_csm_pres_hgt(wks,tilt_int(:,1,:),opts_bdgt)
     plot_fric07  = gsn_csm_pres_hgt(wks,fric_int(:,1,:),opts_bdgt)

     plot_cont06a = gsn_csm_contour(wks,circ_int(:,0,:),opts_cont)
     plot_cont06b = gsn_csm_contour(wks,circ_int(:,0,:),opts_cont)
     plot_cont06c = gsn_csm_contour(wks,circ_int(:,0,:),opts_cont)
     plot_cont06d = gsn_csm_contour(wks,circ_int(:,0,:),opts_cont)
     plot_cont07a = gsn_csm_contour(wks,circ_int(:,1,:),opts_cont)
     plot_cont07b = gsn_csm_contour(wks,circ_int(:,1,:),opts_cont)
     plot_cont07c = gsn_csm_contour(wks,circ_int(:,1,:),opts_cont)
     plot_cont07d = gsn_csm_contour(wks,circ_int(:,1,:),opts_cont)

   ; Overlay contours of circulation onto filled contour plots
     overlay(plot_tot06,plot_cont06a)
     overlay(plot_flux06,plot_cont06b)
     overlay(plot_tilt06,plot_cont06c)
     overlay(plot_fric06,plot_cont06d)
     overlay(plot_tot07,plot_cont07a)
     overlay(plot_flux07,plot_cont07b)
     overlay(plot_tilt07,plot_cont07c)
     overlay(plot_fric07,plot_cont07d)

     panel(0)  = plot_tot06
     panel(1)  = plot_flux06
     panel(2)  = plot_tilt06
     panel(3)  = plot_fric06
     panel(4)  = plot_tot07
     panel(5)  = plot_flux07
     panel(6)  = plot_tilt07
     panel(7)  = plot_fric07

     optsP                       = True
     optsP@gsnFrame              = False    ; Do not advance the frame
     optsP@gsnPanelLabelBar      = True     ; Turn on panel labelbar
     optsP@pmLabelBarWidthF      = 0.6
     optsP@pmLabelBarHeightF     = 0.15
     optsP@lbLabelFontHeightF    = 0.01
     optsP@lbPerimOn             = False
     optsP@lbLabelFont           = "Helvetica"
     optsP@lbTitleString = "Circulation budget tendency terms (10~S~-5~N~ s~S~-1~N~)"
     optsP@lbTitleFontHeightF    = 0.015
     optsP@lbTitleDirection      = "Across"
     optsP@lbTitlePosition       = "Bottom"
     optsP@txString              = "Nepartak: circulation budget terms"
     optsP@gsnPanelFigureStrings = (/"a) Sum of RHS terms", \
                                     "b) Advective flux divergence", \
                                     "c) Tilting", \
				     "d) Friction", \
                                     "e) Sum of RHS terms", \
                                     "f) Advective flux divergence", \
                                     "g) Tilting", \
				     "h) Friction"/)
     optsP@gsnMaximize           = True
     optsP@gsnPanelTop           = 0.98
     optsP@gsnPanelBottom        = 0.02
     optsP@gsnPaperOrientation	 = "landscape"
     optsP@amJust                = "TopLeft"
     optsP@gsnPanelFigureStringsFontHeightF = 0.01 ; Label size (default 0.01)

     gsn_panel(wks,panel,(/2,4/),optsP)              ; Draw as a single plot
     frame(wks)
end