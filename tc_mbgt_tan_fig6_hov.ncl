; Script to reproduce Figure 6 from Nguyen et al. (2011) for a given ensemble member
; Run using:

; ncl 'opt="x11"' 'ens0="em02"' clev=2 clev0=2 dist=3.0 'sy="lee"' 
; ts0=30 tf0=72 mb=3 cont=0 tc_mbgt_fig6_hov.ncl

; 'opt'   = output file format ("pdf" or "x11")
; 'ens0'  = ensemble member (em00 -- em11)
; 'clev'  = pressure level for tangential wind panel (0 = 950, 1 = 900, 2 = 850, 3 = 800)
; 'clev0' = pressure level for vertical velocity panel
; 'dist'  = size of box following storm (degrees)
; 'sy'    = machine (ARCHER [ar] or LEEDS [lee]) for input/output options
; 'mb'	  = 
; 'cont'  = 

load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_func/st_rm.ncl"

begin

;===========
; Prelash
;===========

   if (sy.eq."lee") then
     diri = "/nfs/a37/earshar/ncl_scripts/nepartak/archer_march2018/"
   else if (sy.eq."ar") then
     diri = "$sam/nepartak/images/cart_coords/cp_um/ens/"
   end if
   end if

 ; Calculate storm motion using built-in function ('st_rm')
 ; 'ens0'    = ensemble simulation
 ; 'diri'    = path to input files (see above)
 ; 'dist'    = size of box following storm (degrees)

   storm_rel = st_rm(ens0,diri,3.0)
   u_cyc     = storm_rel[0]             ; Zonal wind
   v_cyc     = storm_rel[1]             ; Meridional wind
   vel_cyc   = storm_rel[2]             ; Vector wind
   numTIMES  = storm_rel[3]             ; Number of times in file
   lat_arr   = storm_rel[4]             ; Latitude array subset (following storm)
   lon_arr   = storm_rel[5]             ; Longitude array subset (following storm)

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
  num_p      = 19                    ; Size of pressure level subset array (L344-345)

;=====================================================
; Define arrays for later use in radius-time plots
;=====================================================

  spd_azi  = new((/dimsizes(radii)/),float)  ; Vector wind          [ {V} ]
  pvrt_azi = new((/dimsizes(radii)/),float)  ; Potential vorticity  [ {PV} ]
  vtan_azi = new((/dimsizes(radii)/),float)  ; Tangential wind      [ {v} ]
  vrad_azi = new((/dimsizes(radii)/),float)  ; Radial wind          [ {u} ]
  vvel_azi = new((/dimsizes(radii)/),float)  ; Vertical velocity    [ {w} ]
  vort_azi = new((/dimsizes(radii)/),float)  ; Rel. vorticity       [ {vort} ]
  avo_azi  = new((/dimsizes(radii)/),float)  ; Abs. vorticity       [ {f + vort} ]
  tht_azi  = new((/dimsizes(radii)/),float)  ; Theta                [ {theta} ]
  frcx_azi = new((/dimsizes(radii)/),float)  ; Vert. diffusion (x)  [ {frcx} ]
  frcy_azi = new((/dimsizes(radii)/),float)  ; Vert. diffusion (y)  [ {frcy} ]
  fric_azi = new((/dimsizes(radii)/),float)  ; Vert. diffusion (tot)[ {fric} ]

; Add metadata for azimuthally averaged variables
  vtan_azi!0            = "rad"
  vtan_azi&rad          = radii
  vtan_azi@description  = "Azimuthally averaged tangential velocity"
  vtan_azi@units        = "m s~S~-1~N~"

  vrad_azi!0            = "rad"
  vrad_azi&rad          = radii
  vrad_azi@description  = "Azimuthally averaged radial velocity"
  vrad_azi@units        = "m s~S~-1~N~"

  spd_azi!0             = "rad"
  spd_azi&rad           = radii
  spd_azi@description   = "Azimuthally averaged vector velocity"
  spd_azi@units         = "m s~S~-1~N~"

  pvrt_azi!0            = "rad"
  pvrt_azi&rad          = radii
  pvrt_azi@description  = "Azimuthally averaged PV"
  pvrt_azi@units        = "10~S~-6~N~K kg~S~-1~N~ m~S~2~N~ s~S~-1~N~"

  vvel_azi!0            = "rad"
  vvel_azi&rad          = radii
  vvel_azi@description  = "Azimuthally averaged vertical velocity"
  vvel_azi@units        = "m s~S~-1~N~"

  avo_azi!0             = "rad"
  avo_azi&rad           = radii
  avo_azi@description   = "Azimuthally averaged absolute vorticity"
  avo_azi@units         = "s~S~-1~N~"

  vort_azi!0            = "rad"
  vort_azi&rad          = radii
  vort_azi@description  = "Azimuthally averaged relative vorticity"
  vort_azi@units        = "s~S~-1~N~"

  tht_azi!0             = "rad"
  tht_azi&rad           = radii
  tht_azi@description   = "Azimuthally averaged potential temperature"
  tht_azi@units         = "K"

  frcx_azi!0            = "rad"
  frcx_azi&rad          = radii
  frcx_azi@description  = "Azimuthally averaged vertical diffusion (x)"
  frcx_azi@units        = "m s~S~-1~N~"

  frcy_azi!0            = "rad"
  frcy_azi&rad          = radii
  frcy_azi@description  = "Azimuthally averaged vertical diffusion (y)"
  frcy_azi@units        = "m s~S~-1~N~"

  fric_azi!0            = "rad"
  fric_azi&rad          = radii
  fric_azi@description  = "Azimuthally averaged vertical diffusion (tot)"
  fric_azi@units        = "m s~S~-1~N~"

;==================================================================
; Define 3D arrays (hold values from above arrays at ALL times)
;==================================================================

; Mean tangential wind 
  vtan_plot0 	   = new((/dimsizes(radii),numTIMES/),float)
  vtan_plot0!0     = "rad"
  vtan_plot0!1     = "time"

; Mean radial wind 
  vrad_plot0       = new((/dimsizes(radii),numTIMES/),float)
  vrad_plot0!0     = "rad"
  vrad_plot0!1     = "time"

; Mean potential vorticity
  pvrt_plot0       = new((/dimsizes(radii),numTIMES/),float)
  pvrt_plot0!0     = "rad"
  pvrt_plot0!1     = "time"

; Mean vertical velocity
  vvel_plot0 	   = new((/dimsizes(radii),numTIMES/),float)
  vvel_plot0!0	   = "rad"
  vvel_plot0!1     = "time"

; Mean relative vorticity
  vort_plot0       = new((/dimsizes(radii),numTIMES/),float)
  vort_plot0!0     = "rad"
  vort_plot0!1     = "time"

; Mean absolute vorticity
  avo_plot0        = new((/dimsizes(radii),numTIMES/),float)
  avo_plot0!0      = "rad"
  avo_plot0!1      = "time"

; Mean potential temperature
  tht_plot0        = new((/dimsizes(radii),numTIMES/),float)
  tht_plot0!0      = "rad"
  tht_plot0!1      = "time"

; Mean vertical diffusion
  fric_plot0       = new((/dimsizes(radii),numTIMES/),float)
  fric_plot0!0     = "rad"
  fric_plot0!1     = "time"

; Radius of maximum wind
  spd_plot0    	   = new(numTIMES,float)
  spd_plot0!0 	   = "time"

; Maximum vertical velocity
  maxw_plot0       = new((/numTIMES/),float)
  maxw_plot0!0     = "time"

; Perturbation tangential wind
  pvtan_plot0     = new((/dimsizes(radii),numTIMES/),float)
  pvtan_plot0!0   = "rad"
  pvtan_plot0!1   = "time"

; Perturbation radial wind
  pvrad_plot0     = new((/dimsizes(radii),numTIMES/),float)
  pvrad_plot0!0   = "rad"
  pvrad_plot0!1   = "time"

; Perturbation vertical velocity
  pvvel_plot0     = new((/dimsizes(radii),numTIMES/),float)
  pvvel_plot0!0   = "rad"
  pvvel_plot0!1   = "time"

; Perturbation relative vorticity
  pvort_plot0     = new((/dimsizes(radii),numTIMES/),float)
  pvort_plot0!0   = "rad"
  pvort_plot0!1   = "time"

; Vertical derivative of mean tangential wind
  mvtdz_plot0     = new((/dimsizes(radii),numTIMES/),float)
  mvtdz_plot0!0   = "rad"
  mvtdz_plot0!1   = "time"

; Vertical derivative of perturbation tangential wind
  pvtdz_plot0     = new((/dimsizes(radii),numTIMES/),float)
  pvtdz_plot0!0   = "rad"
  pvtdz_plot0!1   = "time"

; Eddy radial vorticity flux (RHS #3)
  pvflx_plot0     = new((/dimsizes(radii),numTIMES/),float)
  pvflx_plot0!0   = "rad"
  pvflx_plot0!1   = "time"

; Vertical advection of eddy tangential momentum (RHS #4)
  pvadv_plot0     = new((/dimsizes(radii),numTIMES/),float)
  pvadv_plot0!0   = "rad"
  pvadv_plot0!1   = "time"

; Create time array for later plotting 
  time_int = new(120,"integer")
  time_int = ispan(0,119,1)
  time_int!0 = "time"

;=========================================================
; Additional 3D arrays for momentum budget source terms
;=========================================================

; Arrays for perturbation (eddy) variables (2D)
  ptb_vrad0 = new((/dimsizes(radii),dimsizes(thetas)/),float) ; [u']
  ptb_vort0 = new((/dimsizes(radii),dimsizes(thetas)/),float) ; [vort']
  ptb_vflx0 = new((/dimsizes(radii),dimsizes(thetas)/),float) ; -[u' * vort']
  ptb_vadv0 = new((/dimsizes(radii),dimsizes(thetas)/),float) ; -[w' * dv'/dz]

; 3D arrays for vertical derivative of perturbation tangential wind [-w' * dv'/dz]
  ptb_vvel1 = new((/dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [w']
  ptb_vtan1 = new((/dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [v']
  ptb_vadv1 = new((/dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [R4]

; Arrays for azimuthally averaged eddy variables (1D)
  ptb_vrad  = new((/dimsizes(radii)/),float)       ; [ {u'} ]
  ptb_vtan  = new((/dimsizes(radii)/),float)       ; [ {v'} ]
  ptb_vvel  = new((/dimsizes(radii)/),float)       ; [ {w'} ]
  ptb_vort  = new((/dimsizes(radii)/),float)       ; [ {vort'} ]

  mvtan_dz  = new((/dimsizes(radii)/),float)       ; d{v}/dz
  ptb_vflx  = new((/dimsizes(radii)/),float)       ; -{u' * vort'}    [R3]
  ptb_vadv  = new((/dimsizes(radii)/),float)       ; -{w' * (dv'/dz)} [R4]

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

;=============================
; Read in and analyse files   
;=============================

;    diri        = "$sam/um/cp/ens/20160704T0000Z_ra1t_"
    diri        = "$ens/reg_netcdf/20160704T0000Z_ra1t_"

    input_pb    = diri+ens0+"_pb.nc"
    input_pc    = diri+ens0+"_pc.nc"
    input_pf    = diri+ens0+"_pf.nc"

    numINPUT    = dimsizes(input_pb)        ; Number of input files [4]
    numTIMES    = 120                       ; Total number of times in all files [36]
    ct          = ts0                       ; Counter variable
    ct0		= 0			    ; Counter variable (independent of position)
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
    g  = addfile(input_pc(nf),"r")         ; Read in 'pc' stream [friction]
    h  = addfile(input_pf(nf),"r")         ; Read in 'pf' stream [density on model levs]

;==================================
; Get the variables we will need
;==================================

   do it = ts0, tf0             ; Loop over times in file

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

   print("Working on time: "+time_str(ct0))
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

    pres  = f->p                     ; Pressure levels 
    plevs = dimsizes(pres)  	     ; Size of pressure level array
    pres@units = "hPa"

  ; Read in horizontal velocity, relative vorticity and geopotential height
    u       = f->u(it,:,{t0:t1},{n0:n1})        ; Zonal velocity (m/s)
    v       = f->v(it,:,{t0:t1},{n0:n1})        ; Meridional velocity (m/s)
    pv	    = f->pv(it,:,{t0:t1},{n0:n1})	; Potential vorticity 
    vort    = f->rvor(it,:,{t0:t1},{n0:n1})  	; Relative vorticity (/s)
    z       = f->ht(it,:,{t0:t1},{n0:n1})       ; Geopotential height (m)
    w	    = f->dz_dt(it,:,{t0:t1},{n0:n1})    ; Vertical velocity (m/s)
    t	    = f->temp(it,:,{t0:t1},{n0:n1})	; Temperature (K)

  ; Calculate potential temperature (theta) using NCL function
    th      = pot_temp(pres,t,0,False)

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

     ilev = clev	; Work on both 900 and 850 hPa
     ilev0 = clev0	; Option of a different level for vertical velocity

     p = pres(ilev)
     print("Working on: " + p + " hPa")

     u_plane0      = u(ilev,:,:)
     v_plane0      = v(ilev,:,:)
     vort_plane	   = vort(ilev,:,:)
     avo_plane     = avo(ilev,:,:)
     tht_plane     = th(ilev,:,:)
     geo_plane     = z(ilev,:,:)
     pv_plane	   = pv(ilev,:,:)
     w_plane	   = w(ilev0,:,:)

     u_plane       = u_plane0 - u_cyc(ct)
     v_plane       = v_plane0 - v_cyc(ct)

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

   ; New arrays containing variables we want to plot
     u_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Zonal velocity
     v_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Meridional velocity
     pv_int   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Potential vorticity
     w_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vertical velocity
     v_rad    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Radial velocity
     v_tan    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Tangential velocity
     vort_int = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Relative vorticity
     avo_int  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Absolute vorticity
     tht_int  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Potential temperature

;=========================
; Add variable metadata
;=========================

     copy_VarAtts(u_plane,u_int)
     u_int@description = "Zonal velocity"
     u_int@units       = "m s~S~-1~N~"

     copy_VarAtts(v_plane,v_int)
     v_int@description = "Meridional velocity"
     v_int@units       = "m s~S~-1~N~"

     copy_VarAtts(w_plane,w_int)
     v_int@description = "Vertical velocity"
     v_int@units       = "m s~S~-1~N~"

     copy_VarAtts(pv_plane,pv_int)
     v_int@description = "Potential vorticity"
     v_int@units       = "10~S~-6~N~K kg~S~-1~N~ m~S~2~N~ s~S~-1~N~"

     copy_VarAtts(vort_plane,vort_int)
     vort_int@description = "Relative vorticity"
     vort_int@units       = "s~S~-1~N~"

     copy_VarAtts(avo_plane,avo_int)
     avo_int@description = "Absolute vorticity"
     avo_int@units       = "s~S~-1~N~"

     copy_VarAtts(tht_plane,tht_int)
     tht_int@description = "Potential temperature"
     tht_int@units       = "K"

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

     w_int!0   = "rad"
     w_int!1   = "azi"
     w_int&rad = radii
     w_int&azi = thetas

     pv_int!0   = "rad"
     pv_int!1   = "azi"
     pv_int&rad = radii
     pv_int&azi = thetas

     vort_int!0   = "rad"
     vort_int!1   = "azi"
     vort_int&rad = radii
     vort_int&azi = thetas

     avo_int!0    = "rad"
     avo_int!1    = "azi"
     avo_int&rad  = radii
     avo_int&azi  = thetas

     tht_int!0     = "rad"
     tht_int!1     = "azi"
     tht_int&rad   = radii
     tht_int&azi   = thetas

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

   ; 'xpos' = longitude grid (31 * 73)
   ; 'ypos' = latitude grid (31 * 73)

   ; Loop over azimuth angles to fill entire circular array of data points
     do iang = 0, dimsizes(thetas_rad)-1
       u_int(:,iang)    = linint2_points(lon_sub,lat_sub,u_plane,False,\
                                         xpos(:,iang),ypos(:,iang),0)
       v_int(:,iang)    = linint2_points(lon_sub,lat_sub,v_plane,False,\
                                         xpos(:,iang),ypos(:,iang),0)
       w_int(:,iang)    = linint2_points(lon1_sub,lat1_sub,w_plane,False,\
                                         xpos(:,iang),ypos(:,iang),0)
       pv_int(:,iang)   = linint2_points(lon_sub,lat_sub,pv_plane,False,\
                                         xpos(:,iang),ypos(:,iang),0)
       tht_int(:,iang)   = linint2_points(lon1_sub,lat1_sub,tht_plane,False,\
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
        vtan_azi(irad) = avg(v_tan(irad,:))      ; Tangential wind
	vrad_azi(irad) = avg(v_rad(irad,:))	 ; Radial wind 
	pvrt_azi(irad) = avg(pv_int(irad,:))	 ; Potential vorticity
	vvel_azi(irad) = avg(w_int(irad,:))	 ; Vertical velocity
        vort_azi(irad) = avg(vort_int(irad,:))   ; Relative vorticity
        avo_azi(irad)  = avg(avo_int(irad,:))    ; Absolute vorticity
        tht_azi(irad)  = avg(tht_int(irad,:))    ; Potential temperature
      end do

    ;====================================
    ; Calculate perturbation variables
    ;====================================

      do irad = 0, dimsizes(radii)-1

        ptb_vtan1(irad,:,ct) = v_tan(irad,:)    - vtan_azi(irad) ; [v']
        ptb_vvel1(irad,:,ct) = w_int(irad,:)    - vvel_azi(irad) ; [w']
        ptb_vrad0(irad,:)    = v_rad(irad,:)    - vrad_azi(irad) ; [u']
        ptb_vort0(irad,:)    = vort_int(irad,:) - vort_azi(irad) ; [vort']
        ptb_vflx0(irad,:)    = ptb_vrad0(irad,:) * \        	 ; [u' * vort']
                               ptb_vort0(irad,:)

      end do

    ;====================================================================================
    ; Azimuthally average the perturbation variables ( {} = brackets in Persing et al.)
    ;====================================================================================

      do irad = 0, dimsizes(radii)-1

        ptb_vflx(irad) = avg(ptb_vflx0(irad,:))			 ; [ {u' * vort'} ]
        ptb_vort(irad) = avg(ptb_vort0(irad,:))              	 ; [ {vort'} ]
        ptb_vrad(irad) = avg(ptb_vrad0(irad,:))              	 ; [ {u'} ]

      end do
 
    ; Azimuthally-averaged variables
      vtan_plot0(:,ct)  = (/vtan_azi(:)/)        ; Tangential wind
      vrad_plot0(:,ct)  = (/vrad_azi(:)/)    	 ; Radial wind
      pvrt_plot0(:,ct)  = (/pvrt_azi(:)/)    	 ; PV
      vvel_plot0(:,ct)  = (/vvel_azi(:)/)    	 ; Vertical velocity
      vort_plot0(:,ct)  = (/vort_azi(:)/)        ; Relative vorticity
      avo_plot0(:,ct)   = (/avo_azi(:)/)         ; Absolute vorticity
      tht_plot0(:,ct)   = (/tht_azi(:)/)         ; Potential temperature

    ; Perturbation variables
      pvflx_plot0(:,ct) = (/ptb_vflx(:)/)        ; Eddy vorticity flux (RHS #3)
      pvort_plot0(:,ct) = (/ptb_vort(:)/)        ; Relative vorticity
      pvrad_plot0(:,ct) = (/ptb_vrad(:)/)        ; Radial wind

    ; Tidy up
      delete([/u_int,v_int,w_int,pv_int,v_tan,v_rad,tht_int,avo_int,vort_int/])
      delete([/u_plane0,v_plane0,u_plane,v_plane,w_plane,pv_plane,vort_plane,geo_plane/])
      delete([/lat_plane,lon_plane,lat1_plane,lon1_plane,avo_plane,tht_plane/])
      delete([/lat_sub,lat1_sub,lon_sub,lon1_sub/])

    ct  = ct + 1		    ; Counter variable (time)
    ct0 = ct0 + 1		    ; Independent counter variable (time)

  ; Tidy up before next iteration
    delete([/year,month,day,hour,minute,second,utc_date/])
    delete([/lon,lat,lon1,lat1,time,times,u,v,z,t,w,vort,pres,plevs/])

   end do     ; End time loop (do it = 0, times-1)

;=====================================================
; Calculate radial and tangential wind tendencies
;=====================================================

; Initialise new arrays 
  vtan_tend0 	  = vtan_plot0
  vrad_tend0 	  = vrad_plot0
 
  vtan_tend0!0    = "rad"
  vtan_tend0&rad  = radii
  vtan_tend0!1    = "time"
  vtan_tend0&time = time_int

  vrad_tend0!0    = "rad"
  vrad_tend0&rad  = radii
  vrad_tend0!1    = "time"
  vrad_tend0&time = time_int

  do it = ts0, tf0

   tp1  = min((/it+1,tf0/)) 
   tm1	= max((/it-1,ts0/))

   vtan_tend0(:,it) = (vtan_plot0(:,tp1) - vtan_plot0(:,tm1)) / 2 ; Units [m s-1 h-1] 
   vrad_tend0(:,it) = (vrad_plot0(:,tp1) - vrad_plot0(:,tm1)) / 2 ; Units [m s-1 h-1]

  end do     ; End time loop (do it = ts0, tf0)

; Calculate terms from Persing et al. (2013) - their Eq.(12)
; LHS     : d{v}/dt          - Rate of change of mean tangential wind
; RHS #1  : -{u} * {f+vort}  - Mean radial influx of absolute vertical vorticity
; RHS #2  : -{w} * (d{v}/dz) - Mean vert. advection of mean tangential momentum
; RHS #3  : -{u' * vort'}    - Eddy radial vorticity flux
; RHS #4  : -{w' * (dv'/dz)} - Vertical advection of eddy tangential momentum
; RHS #5  : c_p{(theta'/r) * (dpi'/dlambda)} - Azimuthal pert. PGF per unit mass
; RHS #6  : {D_v1}           - Combined diffusive and PBL tendency (vertical)
; RHS #7  ; {D_v2}           - Combined diffusive and PBL tendency (horizontal)

; Radial influx of absolute (vertical) vorticity
  rhs1A      = -vrad_plot0 * avo_plot0
  rhs1A      = rhs1A * 3600
  rhs1A@units = "m s~S~-1~N~ h~S~-1~N~"

; Mean vert. adv. of mean tangential momentum
  rhs2A     = -vvel_plot0 * mvtdz_plot0
  rhs2A      = rhs2A * 3600
  rhs2A@units = "m s~S~-1~N~ h~S~-1~N~"

; Eddy radial vorticity flux
  rhs3A      = -pvflx_plot0
  rhs3A      = rhs3A * 3600
  rhs3A@units = "m s~S~-1~N~ h~S~-1~N~"

; Vertical advection of eddy tangential momentum
  rhs4A = -pvadv_plot0
  rhs4A = rhs4A * 3600
  rhs4A@units = "m s~S~-1~N~ h~S~-1~N~"

  rhs1A!0         = "rad"
  rhs1A&rad       = radii
  rhs1A!1         = "time"
  rhs1A&time      = time_int

  rhs2A!0         = "rad"
  rhs2A&rad       = radii
  rhs2A!1         = "time"
  rhs2A&time      = time_int

  rhs3A!0         = "rad"
  rhs3A&rad       = radii
  rhs3A!1         = "time"
  rhs3A&time      = time_int

  rhs4A!0         = "rad"
  rhs4A&rad       = radii
  rhs4A!1         = "time"
  rhs4A&time      = time_int

  pres  = f->p                     ; Pressure levels
  plevs = dimsizes(pres)           ; Size of pressure level array
  pres@units = "hPa"

;======================================
; Reorder dimensions before plotting 
;======================================

 ; Arrays now ordered [TIME * RADIUS] rather than [RADIUS * TIME]
   vtan_tend = vtan_tend0(time|:,rad|:)
   vrad_tend = vrad_tend0(time|:,rad|:)

   vtan_plot = vtan_plot0(time|:,rad|:)
   vrad_plot = vrad_plot0(time|:,rad|:)
   vvel_plot = vvel_plot0(time|:,rad|:)

 ; Add metadata before plotting 
   vvel_plot!0     = "time"
   vvel_plot&time  = time_int
   vvel_plot!1     = "rad"
   vvel_plot&rad   = radii

   vtan_plot!0     = "time"
   vtan_plot&time  = time_int
   vtan_plot!1     = "rad"
   vtan_plot&rad   = radii

   vrad_plot!0     = "time"
   vrad_plot&time  = time_int
   vrad_plot!1     = "rad"
   vrad_plot&rad   = radii

 ; Momentum budget source terms
   rhs1      = rhs1A(time|:,rad|:)
   rhs2      = rhs2A(time|:,rad|:)
   rhs3      = rhs3A(time|:,rad|:)
   rhs4      = rhs4A(time|:,rad|:)

;=====================================================
; Create a panel plot of Figs. 6a and 6b from N11 
;=====================================================

  ; Output file location and type
    output = "$sam/nepartak/images/cyl_coords/mbgt_hov_all_"+p+"_"+ens0
    wks = gsn_open_wks(opt,output)

  ; Set up panel plot
    panel = new(3,graphic)
    gsn_define_colormap(wks,"n11")

;==========================
; Options for plotting 
;==========================

  ; Tangential velocity
    opts_vtan                              = True
    opts_vtan@cnFillOn                     = True
    opts_vtan@cnLineLabelInterval          = 2.0
    opts_vtan@cnLineLabelFontHeightF       = 0.012
    opts_vtan@cnLineLabelBackgroundColor   = "transparent"
    opts_vtan@cnLineLabelPlacementMode     = "constant"
    opts_vtan@cnLinesOn                    = False ; Contour lines off
    opts_vtan@cnInfoLabelOn                = False
    opts_vtan@cnLevelSelectionMode         = "ExplicitLevels"
    opts_vtan@cnLevels                     = (/1.0, 2.0, 2.5, 3.0, 3.5, 4.0, \
                                               4.5, 5.0, 6.0, 7.0/)
    opts_vtan@cnFillColors                 = (/0,2,3,4,5,6,7,8,9,10,11/)
    opts_vtan@gsnPaperOrientation          = "landscape"
    opts_vtan@tiMainString                 = ""
    opts_vtan@tiMainFontHeightF            = 0.0125
    opts_vtan@gsnLeftString                = ""
    opts_vtan@gsnRightString               = ""
    opts_vtan@gsnMaximize		   = True
    opts_vtan@lbLabelBarOn		   = True		   ; Labelbar on/off
    opts_vtan@lbBoxEndCapStyle		   = "TriangleBothEnds"	   ; Labelbar end shape
    opts_vtan@pmLabelBarWidthF		   = 0.32		   ; Labelbar width
    opts_vtan@pmLabelBarHeightF		   = 0.08		   ; Labelbar height
    opts_vtan@lbLabelFontHeightF	   = 0.0125		   ; Labelbar font size
    opts_vtan@lbLabelFont		   = "Helvetica"	   ; Labelbar font
    opts_vtan@lbPerimOn			   = False		   ; Perimeter on/off
    opts_vtan@gsnDraw			   = False		   ; Do not draw plot
    opts_vtan@gsnFrame			   = False		   ; Do not advance frame
    opts_vtan@vpWidthF            	   = 0.30		   ; Match width in N11
    opts_vtan@vpHeightF                    = 0.50		   ; Match height in N11

  ; Additional plotting resources
    opts_vtan@tiYAxisString                = "Time (hour)"
    opts_vtan@tiXAxisString		   = "Radius (km)"
    opts_vtan@tiXAxisFontHeightF           = 0.018
    opts_vtan@trYLog                       = False
    opts_vtan@trYMinF                      = 30      ; Starting point (T+...) for plot
    opts_vtan@trYMaxF                      = 72      ; Ending point (T+...) for plot
    opts_vtan@tmXTOn                       = "False" ; Turn off top x-axis TM
    opts_vtan@tmYROn                       = "False" ; Turn off right y-axis TM
    opts_vtan@gsnMaximize                  = True    ; Maximise plot size
    opts_vtan@gsnAddCyclic                 = False
    opts_vtan@tmXBMode                     = "Explicit"                ; Choose tm mode
    opts_vtan@tmXBValues                   = (/0.0,0.5,1.0,1.5/)       ; Tm positions
    opts_vtan@tmXBLabels                   = (/"0","50","100","150"/)  ; Tm values
    opts_vtan@tmXBMinorOn                  = False

    opts_vtan@pmLabelBarOrthogonalPosF     = 0.10    ; Move lb up/down (higher/lower)
    opts_vtan@tiXAxisOffsetYF		   = 0.10    ; Move x-axis title up/down (h/l)

    opts_vtan@trXMinF			   = 0.0		       ; Set x-axis min
    opts_vtan@trXMaxF			   = 1.0		       ; Set x-axis max
    opts_vtan@tmYLMode			   = "Manual"
    opts_vtan@tmYLTickStartF		   = opts_vtan@trYMinF	       ; See L730
    opts_vtan@tmYLTickEndF		   = opts_vtan@trYMaxF	       ; See L731
    opts_vtan@tmYLTickSpacingF		   = 6
    opts_vtan@tmYLMinorOn		   = False

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

  ; Resources to overlay zero line
    opts_zero                             = True
    opts_zero@cnFillOn                    = False
    opts_zero@cnLineColor                 = "gray82"
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

  ; Radial velocity resources
    opts_vrad				   = opts_vtan
    delete(opts_vrad@cnLevels)
    delete(opts_vrad@cnFillColors)
    opts_vrad@cnLevels			   = (/0.0, -0.5, -1.0, -1.5, -2.0, \
                                               -2.5, -3.0, -4.0, -5.0/)
    opts_vrad@cnFillColors                 = (/12,11,10,9,8,7,6,5,4,0/)

  ; Vertical velocity resources
    opts_vvel				   = opts_vtan
    delete(opts_vvel@cnLevels)
    delete(opts_vvel@cnFillColors)
    opts_vvel@cnLevels			   = (/0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, \
    					       0.7, 0.8, 0.9, 1.0/)
    opts_vvel@cnFillColors		   = (/0,2,3,4,5,6,7,8,9,10,11,12/)

  ; Radius of maximum wind (RMW) resources 
    opts_rmw	     	       	           = True
    opts_rmw@xyLineColors      		   = "red3"
    opts_rmw@xyLineThicknesses 		   = 3.0
    opts_rmw@xyDashPatterns    		   = 0.0
    opts_rmw@gsnDraw                       = False       ; Do not draw the plot
    opts_rmw@gsnFrame                      = False       ; Do not advance the frame
    opts_rmw@tiMainString                  = ""
    opts_rmw@gsnLeftString                 = ""
    opts_rmw@gsnRightString                = ""

  ; Momentum budget term resources (negative dashed, positive solid)
    opts_mbgt				   = opts_cont
    opts_mbgt@cnLineColor		   = "black"
    delete(opts_mbgt@cnLevels)
    opts_mbgt@cnLevels			   = (/-30., -15., 15., 30./)
    opts_mbgt@cnLineLabelsOn	           = False
    opts_mbgt@cnLineThicknessF		   = 2.5
    opts_mbgt@gsnContourNegLineDashPattern = 14

;==============
; Panel plot 
;==============

  ; Read in RMW data from text file
    spd_plot0 = asciiread("./"+ens0+"_rmw.txt",(/numTIMES/),"float")

  ;===============================
  ; Tangential wind (LHS panel)
  ;===============================
    vtan0 = gsn_csm_contour(wks,vtan_tend,opts_vtan)	; Tangential wind tendency 
    zero0 = gsn_csm_contour(wks,vtan_tend,opts_zero)	; Zero line (vtan tendency)
    zero  = ColorNegDashZeroPosContour(zero0,"transparent","gray69","transparent")
    cont0 = gsn_csm_contour(wks,vtan_plot,opts_cont)	; Tangential wind 
    rmw0  = gsn_csm_xy(wks,spd_plot0,time_int,opts_rmw)	; RMW

    if(mb.eq.1) then
     mbgt1 = gsn_csm_contour(wks,rhs1(:,:),opts_mbgt)    ; MBGT term #1
    else if (mb.eq.3) then
     mbgt1 = gsn_csm_contour(wks,rhs3(:,:),opts_mbgt)    ; MBGT term #3
    end if 
    end if

    overlay(vtan0,zero)					; Overlay vtan tendency zero line
    overlay(vtan0,cont0)				; Overlay tangential wind
    overlay(vtan0,rmw0)					; Overlay RMW

    if (cont.eq.1) then
     overlay(vtan0,mbgt1)                               ; Overlay MBGT source term
     panel(0) = vtan0					; Plot first panel
     delete([/zero0,zero/])				; Tidy up 
     delete(mbgt1)
    else if (cont.eq.0) then
     panel(0) = vtan0                                   ; Plot first panel
     delete([/zero0,zero/])                             ; Tidy up
     delete(mbgt1)
    end if	
    end if 

  ;===============================
  ; Radial wind (middle panel)
  ;===============================
    vrad0 = gsn_csm_contour(wks,vrad_tend,opts_vrad)	; Radial wind tendency
    cont0 = gsn_csm_contour(wks,vtan_plot,opts_cont)	; Tangential wind
    rmw0  = gsn_csm_xy(wks,spd_plot0,time_int,opts_rmw)	; RMW 

    if(mb.eq.1) then
     mbgt1 = gsn_csm_contour(wks,rhs1(:,:),opts_mbgt)    ; MBGT term #1
    else if (mb.eq.3) then
     mbgt1 = gsn_csm_contour(wks,rhs3(:,:),opts_mbgt)    ; MBGT term #3
    end if
    end if

    overlay(vrad0,cont0)    				; Overlay tangential wind
    overlay(vrad0,rmw0)					; Overlay RMW

    if (cont.eq.1) then
     overlay(vrad0,mbgt1)				; Overlay MBGT term #1
     panel(1) = vrad0					; Plot second panel
     delete(mbgt1)
    else if (cont.eq.0) then
     panel(1) = vrad0                                   ; Plot second panel
     delete(mbgt1)
    end if 
    end if

  ;===============================
  ; Vertical velocity (RHS panel)
  ;===============================
    vvel0 = gsn_csm_contour(wks,vvel_plot,opts_vvel)	; Vertical velocity
    zero0 = gsn_csm_contour(wks,vvel_plot,opts_zero)	; Zero line (vertical velocity)
    zero  = ColorNegDashZeroPosContour(zero0,"transparent","gray69","transparent")
    rmw0  = gsn_csm_xy(wks,spd_plot0,time_int,opts_rmw)	; RMW

    if(mb.eq.1) then
     mbgt1 = gsn_csm_contour(wks,rhs1(:,:),opts_mbgt)    ; MBGT term #1
    else if (mb.eq.3) then
     mbgt1 = gsn_csm_contour(wks,rhs3(:,:),opts_mbgt)    ; MBGT term #3
    end if
    end if

    overlay(vvel0,zero)					; Overlay vvel zero line
    overlay(vvel0,rmw0)					; Overlay RMW

    if (cont.eq.1) then
     overlay(vvel0,mbgt1)				; Overlay MBGT term #1 
     panel(2) = vvel0					; Plot third panel
     delete([/zero0,zero,mbgt1/])			; Tidy up
    else if (cont.eq.0) then
     panel(2) = vvel0                                   ; Plot third panel
     delete([/zero0,zero,mbgt1/])                       ; Tidy up
    end if
    end if 

  ;===============================
  ; Final panel plot resources 
  ;===============================
    optsP                       = True
    optsP@gsnFrame              = False    ; Do not advance the frame
    optsP@gsnPanelLabelBar      = False    ; Turn on panel labelbar
    optsP@txString              = "Figure 6 from Nguyen et al (2011): simulation "\
    				   +ens0+" ("+pres(clev)+" hPa)"
    optsP@gsnPanelFigureStrings = (/"a) vtan_tend ","b) vrad_tend ", \
    				    "c) vvel"/)
    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.90
    optsP@gsnPanelBottom        = 0.10
    optsP@amJust                = "TopLeft"
    optsP@gsnPanelFigureStringsFontHeightF = 0.0099 ; Reduce label size (default 0.01)

    gsn_panel(wks,panel,(/1,3/),optsP)             ; Draw as a single plot
    frame(wks)

end 