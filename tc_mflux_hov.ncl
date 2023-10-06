; Produce Hovmoller plots of momentum budget terms 

; Run using:

; ncl 'opt="x11"' 'ens0="em02"' clev=2 dist=3.0 'clr="bwr_new"' 'sy="lee"'
; ptb=1 ts0=30 tf0=72 lb=1 tc_mbgt_hov.ncl

; 'opt'   = output file format ("pdf" or "x11")
; 'ens0'  = ensemble member (em00 -- em11)
; 'clev'  = pressure level for tangential wind panel (0 = 950, 1 = 900, 2 = 850, 3 = 800)
; 'dist'  = size of box following storm (degrees)
; 'ptb'	  = type of plot (ptb=1 --> perturbation; ptb=0 --> mean)
; 'ts0'	  = starting time index (0–119)
; 'tf0'   = ending time index (0–119)
; 'sy'	  = machine (ARCHER [ar] or LEEDS [lee]) for input/output options
; 'lb'    = individual labelbars (lb=1) or a single panel labelbar (lb=0)

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"
load "$sam/ncl_func/st_rm.ncl"

;if (sy.eq."lee") then
;  path = "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"
;else if (sy.eq."ar") then
;  path = "$LIB/gsn_csm.ncl"
;end if
;end if

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
   u_cyc     = storm_rel[0]		; Zonal	wind
   v_cyc     = storm_rel[1]    	        ; Meridional wind
   vel_cyc   = storm_rel[2]	        ; Vector wind
   numTIMES  = storm_rel[3] 	        ; Number of times in file
   lat_arr   = storm_rel[4]	      	; Latitude array subset (following storm)
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
  num_p	     = 19		     ; Size of pressure level subset array (L344-345)

;=====================================================
; Define arrays for later use in radius-time plots
;=====================================================

  spd_azi  = new((/num_p,dimsizes(radii)/),float)  ; Vector wind	  [ {V} ]
  pvrt_azi = new((/num_p,dimsizes(radii)/),float)  ; Potential vorticity  [ {PV} ]
  vtan_azi = new((/num_p,dimsizes(radii)/),float)  ; Tangential wind      [ {v} ]
  vrad_azi = new((/num_p,dimsizes(radii)/),float)  ; Radial wind          [ {u} ]
  vvel_azi = new((/num_p,dimsizes(radii)/),float)  ; Vertical velocity    [ {w} ]
  vort_azi = new((/num_p,dimsizes(radii)/),float)  ; Rel. vorticity       [ {vort} ]
  avo_azi  = new((/num_p,dimsizes(radii)/),float)  ; Abs. vorticity       [ {f + vort} ]
  flux_azi = new((/num_p,dimsizes(radii)/),float)  ; Mass flux            [{mflux} ]
  tht_azi  = new((/num_p,dimsizes(radii)/),float)  ; Theta                [ {theta} ]
  frcx_azi = new((/num_p,dimsizes(radii)/),float)  ; Vert. diffusion (x)  [ {frcx} ]
  frcy_azi = new((/num_p,dimsizes(radii)/),float)  ; Vert. diffusion (y)  [ {frcy} ]
  fric_azi = new((/num_p,dimsizes(radii)/),float)  ; Vert. diffusion (tot)[ {fric} ]

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

; Mean potential vorticity
  pvrt_plot0      = new((/num_p,dimsizes(radii),numTIMES/),float)
  pvrt_plot0!0    = "lev"
  pvrt_plot0!1    = "rad"
  pvrt_plot0&rad  = radii
  pvrt_plot0!2    = "time"

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

; Mean mass flux
  flux_plot0      = new((/num_p,dimsizes(radii),numTIMES/),float)
  flux_plot0!0    = "lev"
  flux_plot0!1    = "rad"
  flux_plot0&rad  = radii
  flux_plot0!2    = "time"

; Mean potential temperature
  tht_plot0 	  = new((/num_p,dimsizes(radii),numTIMES/),float)
  tht_plot0!0	  = "lev"
  tht_plot0!1 	  = "rad"
  tht_plot0&rad   = radii
  tht_plot0!2 	  = "time"

; Mean vertical diffusion
  fric_plot0       = new((/num_p,dimsizes(radii),numTIMES/),float)
  fric_plot0!0     = "lev"
  fric_plot0!1     = "rad"
  fric_plot0&rad   = radii
  fric_plot0!2     = "time"

; Radius of maximum wind
  spd_plot0   	  = new((/numTIMES/),float)
  spd_plot0!0 	  = "time"

; Maximum vertical velocity
  maxw_plot0	   = new((/num_p,numTIMES/),float)
  maxw_plot0!0	   = "lev"
  maxw_plot0!1	   = "time"

; Perturbation tangential wind
  pvtan_plot0     = new((/num_p,dimsizes(radii),numTIMES/),float)
  pvtan_plot0!0   = "lev"
  pvtan_plot0!1   = "rad"
  pvtan_plot0&rad = radii
  pvtan_plot0!2   = "time"

; Perturbation radial wind
  pvrad_plot0     = new((/num_p,dimsizes(radii),numTIMES/),float)
  pvrad_plot0!0   = "lev"
  pvrad_plot0!1   = "rad"
  pvrad_plot0&rad = radii
  pvrad_plot0!2   = "time"

; Perturbation vertical velocity
  pvvel_plot0     = new((/num_p,dimsizes(radii),numTIMES/),float)
  pvvel_plot0!0   = "lev"
  pvvel_plot0!1   = "rad"
  pvvel_plot0&rad = radii
  pvvel_plot0!2   = "time"

; Perturbation relative vorticity
  pvort_plot0     = new((/num_p,dimsizes(radii),numTIMES/),float)
  pvort_plot0!0   = "lev"
  pvort_plot0!1   = "rad"
  pvort_plot0&rad = radii
  pvort_plot0!2   = "time"

; Vertical derivative of mean tangential wind
  mvtdz_plot0     = new((/num_p,dimsizes(radii),numTIMES/),float)
  mvtdz_plot0!0   = "lev"
  mvtdz_plot0!1   = "rad"
  mvtdz_plot0&rad = radii
  mvtdz_plot0!2   = "time"

; Vertical derivative of perturbation tangential wind
  pvtdz_plot0     = new((/num_p,dimsizes(radii),numTIMES/),float)
  pvtdz_plot0!0   = "lev"
  pvtdz_plot0!1   = "rad"
  pvtdz_plot0&rad = radii
  pvtdz_plot0!2   = "time"

; Eddy radial vorticity flux (RHS #3)
  pvflx_plot0     = new((/num_p,dimsizes(radii),numTIMES/),float)
  pvflx_plot0!0   = "lev"
  pvflx_plot0!1   = "rad"
  pvflx_plot0&rad = radii
  pvflx_plot0!2   = "time"

; Vertical advection of eddy tangential momentum (RHS #4)
  pvadv_plot0     = new((/num_p,dimsizes(radii),numTIMES/),float)
  pvadv_plot0!0   = "lev"
  pvadv_plot0!1   = "rad"
  pvadv_plot0&rad = radii
  pvadv_plot0!2   = "time"

;=========================================================
; Additional 3D arrays for momentum budget source terms 
;=========================================================

; Arrays for perturbation (eddy) variables (3D)
  ptb_vrad0 = new((/num_p,dimsizes(radii),dimsizes(thetas)/),float) ; [u']
  ptb_vort0 = new((/num_p,dimsizes(radii),dimsizes(thetas)/),float) ; [vort']
  ptb_vflx0 = new((/num_p,dimsizes(radii),dimsizes(thetas)/),float) ; -[u' * vort']
  ptb_vadv0 = new((/num_p,dimsizes(radii),dimsizes(thetas)/),float) ; -[w' * dv'/dz]

; 4D arrays for vertical derivative of perturbation tangential wind [-w' * dv'/dz]
  ptb_vvel1 = new((/num_p,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [w']
  ptb_vtan1 = new((/num_p,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [v']
  ptb_vadv1 = new((/num_p,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [R4]

; Arrays for azimuthally averaged eddy variables (2D)
  ptb_vrad  = new((/num_p,dimsizes(radii)/),float)       ; [ {u'} ]
  ptb_vtan  = new((/num_p,dimsizes(radii)/),float)       ; [ {v'} ]
  ptb_vvel  = new((/num_p,dimsizes(radii)/),float)       ; [ {w'} ]
  ptb_vort  = new((/num_p,dimsizes(radii)/),float)       ; [ {vort'} ]

  mvtan_dz  = new((/num_p,dimsizes(radii)/),float)       ; d{v}/dz
  ptb_vflx  = new((/num_p,dimsizes(radii)/),float)       ; -{u' * vort'}    [R3]
  ptb_vadv  = new((/num_p,dimsizes(radii)/),float)       ; -{w' * (dv'/dz)} [R4]

; Add metadata for azimuthally averaged variables
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

  spd_azi!0            = "lev"
  spd_azi!1            = "rad"
  spd_azi&rad          = radii
  spd_azi@description  = "Azimuthally averaged vector velocity"
  spd_azi@units        = "m s~S~-1~N~"

  pvrt_azi!0           = "lev"
  pvrt_azi!1           = "rad"
  pvrt_azi&rad         = radii
  pvrt_azi@description = "Azimuthally averaged PV"
  pvrt_azi@units       = "10~S~-6~N~K kg~S~-1~N~ m~S~2~N~ s~S~-1~N~"

  vvel_azi!0           = "lev"
  vvel_azi!1           = "rad"
  vvel_azi&rad         = radii
  vvel_azi@description = "Azimuthally averaged vertical velocity"
  vvel_azi@units       = "m s~S~-1~N~"

  avo_azi!0            = "lev"
  avo_azi!1            = "rad"
  avo_azi&rad          = radii
  avo_azi@description  = "Azimuthally averaged absolute vorticity"
  avo_azi@units        = "s~S~-1~N~"

  flux_azi!0            = "lev"
  flux_azi!1            = "rad"
  flux_azi&rad          = radii
  flux_azi@description  = "Azimuthally averaged mass flux"
  flux_azi@units        = "s~S~-1~N~"

  vort_azi!0           = "lev"
  vort_azi!1           = "rad"
  vort_azi&rad         = radii
  vort_azi@description = "Azimuthally averaged relative vorticity"
  vort_azi@units       = "s~S~-1~N~"

  tht_azi!0             = "lev"
  tht_azi!1             = "rad"
  tht_azi&rad           = radii
  tht_azi@description   = "Azimuthally averaged potential temperature"
  tht_azi@units         = "K"

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
;    diri        = "$sam/um/cp/ens/20160704T0000Z_ra1t_"
    diri        = "$ens/reg_netcdf/20160704T0000Z_ra1t_"

    input_pb    = diri+ens0+"_pb.nc"
    input_pc	= diri+ens0+"_pc.nc"
    input_pf	= diri+ens0+"_pf.nc"

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

  ; Create array to hold values of frictional tendency
    dtaux_dz = new((/numTIMES,plevs,d0,d1/),float)
    dtauy_dz = new((/numTIMES,plevs,d0,d1/),float)

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

;=======================
; Calculate mass flux
;=======================

  ; First calculate dry air density
    pres_fr = conform_dims(dimsizes(t),pres,0) ; Create p array [plevs * d(0) * d(1)]
    pres_fr = pres_fr * 100
    pres_fr@units = "Pa"

    rho   = pres_fr(:,:,:) / (287 * t(:,:,:))
    rho@description = "Dry air density"
    rho@units = "kg m~S~-3~N~"

  ; Now calculate mass flux
    mflux = rho * w
    copy_VarCoords(w,mflux)
    copy_VarAtts(w,mflux)
    mflux@units = "kg m~S~-2~N~ s~S~-1~N~"
    mflux@standard_name = "vertical_mass_flux"
    mflux@long_name     = "Vertical mass flux"
    mflux@title         = "Vertical mass flux"
    mflux@name          = "vert_mflux"

;==============================================
; Get the variables we will need [pc stream]
;==============================================

  ; Get additional variables from the model-level file
    lon1F = g->longitude_1({n0:n1})            ; longitude [xxx grid points]
    lat1F = g->latitude_1({t0:t1})             ; latitude  [xxx grid points]
    lonF  = g->longitude({n0:n1})              ; longitude [xxx grid points]
    latF  = g->latitude({t0:t1})               ; latitude  [xxx grid points]
    hyb   = g->hybrid_ht(:)                    ; 63 hybrid levels (5 m up to ~20 km)
    strx  = g->taux(it,:,{t0:t1},{n0:n1})      ; Wind stress, x-comp (N m/2)
    stry  = g->tauy(it,:,{t0:t1},{n0:n1})      ; Wind stress, y-comp (N m/2)

    dxF   = lonF(1) - lonF(0) ; Grid spacing (longitude)
    dyF   = latF(1) - latF(0) ; Grid spacing (latitude)

;=========================================================================
; Get density from pf stream [same size as 'strx', smaller than 'stry']
;=========================================================================

    rho1  = h->field27(it,:,{t0:t1},{n0:n1})  ; Dry air density (kg/m-3)

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

  ; Arrays to hold subset of 'strx' and 'stry', on levels in 'lev_arr'
  ; The dimensions of 'dsize' latitude (0) and longitude (1)
    fric_x_hgt  = new((/dimsizes(lev_arr),d0,d1/),float)
    fric_y_hgt  = new((/dimsizes(lev_arr),d0,d1/),float)

  ; Array of pressure levels with the same dimensions as 'lev_arr'
    pres       = (/1000,950,900,850,800,750,700,650,600,550,\
                   500,450,400,350,300,250,200,150,100/)
    pres@units = "hPa"
    plevs      = dimsizes(pres)

  ; Create array of pressure values with same dimensions as 'fric_x_hgt'
    pres_fr    = conform_dims(dimsizes(fric_x_hgt), pres, 0)
    pres_fr = pres_fr * 100
    pres_fr@units = "Pa"

  ; TC latitude  -- centre(it,0)
  ; TC longitude -- centre(it,1)

;===================================================
; Add pressure level metadata for selected arrays 
;===================================================

; 3D arrays holding values for final plotting  
  vtan_plot0&lev = pres
  vrad_plot0&lev = pres
  pvrt_plot0&lev = pres
  vvel_plot0&lev = pres
  vort_plot0&lev = pres
  avo_plot0&lev  = pres
  flux_plot0&lev = pres
  tht_plot0&lev  = pres
  fric_plot0&lev = pres
  maxw_plot0&lev = pres

; Arrays holding azimuthally-averaged variables
  vtan_azi&lev   = pres
  vrad_azi&lev   = pres
  spd_azi&lev    = pres
  pvrt_azi&lev   = pres
  vvel_azi&lev   = pres
  avo_azi&lev    = pres
  flux_azi&lev	 = pres
  vort_azi&lev   = pres
  tht_azi&lev    = pres  

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

    do ilev = 0, plevs-1   ; All pressure levels in vertical level array subset

   ;========================================
   ; Model levels ['pc' and 'pf' streams]
   ;========================================

     hy = sprintf("%0.0f",hyb(ilev))
     print("Working on hybrid level " + ilev+ " ("+hy+" m AGL)")
     m0 = lev_arr(ilev)       ; Selected model levels, chosen to match pressure levels

   ; Calculate components of friction term
     ip1 = min((/m0+1,lev_arr(l1)/)) ; Account for top and bottom levels
     im1 = max((/m0-1,lev_arr(l0)/))

   ; Current level [wind stress]
     taux   = strx(m0,0:d0-1,0:d1-1)
     tauy   = stry(m0,0:d0-1,0:d1-1)

   ; Density ['pf' stream]
     rho0   = rho1(m0,0:d0-1,0:d1-1)
     rho_a  = (1/rho0) ; Specific volume (inverse density)

   ; Levels above/below for vertical derivatives [wind stress]
     tau_x1 = strx(ip1,0:d0-1,0:d1-1)
     tau_x0 = strx(im1,0:d0-1,0:d1-1)
     tau_y1 = stry(ip1,0:d0-1,0:d1-1)
     tau_y0 = stry(im1,0:d0-1,0:d1-1)

   ; Vertical derivative of 'strx'
     dtaux_dz(ct,ilev,:,:) = (tau_x1 - tau_x0) / (hyb(ip1) - hyb(im1))

   ; Vertical derivative of 'stry'
     dtauy_dz(ct,ilev,:,:) = (tau_y1 - tau_y0) / (hyb(ip1) - hyb(im1))

   ; X-component of friction
     frcx_plane             = rho_a(:,:) * dtaux_dz(ct,ilev,:,:)
     copy_VarCoords(taux,frcx_plane)
     copy_VarMeta(taux,frcx_plane)
     frcx_plane@units       = "m s~S~-2~N~"
     frcx_plane@title       = "x-component of friction"
     frcx_plane@long_name   = "x-component of friction"
     frcx_plane@name        = "fric_x"

   ; Y-component of friction
     frcy_plane             = rho_a(:,:) * dtauy_dz(ct,ilev,:,:)
     copy_VarCoords(tauy,frcy_plane)
     copy_VarMeta(tauy,frcy_plane)
     frcy_plane@units       = "m s~S~-2~N~"
     frcy_plane@title       = "y-component of friction"
     frcy_plane@long_name   = "y-component of friction"
     frcy_plane@name        = "fric_y"

   ;=================================
   ; Pressure levels ['pb' stream]
   ;=================================

     p = pres(ilev)
     print("Working on: " + p + " hPa")

     u_plane0      = u(ilev,:,:)
     v_plane0      = v(ilev,:,:)
     vort_plane	   = vort(ilev,:,:)
     avo_plane	   = avo(ilev,:,:)
     flux_plane	   = mflux(ilev,:,:)
     tht_plane	   = th(ilev,:,:)	
     geo_plane     = z(ilev,:,:)
     pv_plane	   = pv(ilev,:,:)
     w_plane	   = w(ilev,:,:)

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

   ; New arrays containing variables we want to plot
     u_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Zonal velocity
     v_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Meridional velocity
     spd_int  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vector velocity (RMW)
     pv_int   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Potential vorticity
     w_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vertical velocity
     v_rad    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Radial velocity
     v_tan    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Tangential velocity
     vort_int = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Relative vorticity
     avo_int  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Absolute vorticity
     flux_int = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Mass flux
     tht_int  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Potential temperature
     frcx_int = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Diffusion (x)
     frcy_int = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Diffusion (y)
     fric_int = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Diffusion (total)

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
     w_int@description = "Vertical velocity"
     w_int@units       = "m s~S~-1~N~"

     copy_VarAtts(pv_plane,pv_int)
     pv_int@description = "Potential vorticity"
     pv_int@units       = "10~S~-6~N~K kg~S~-1~N~ m~S~2~N~ s~S~-1~N~"

     copy_VarAtts(vort_plane,vort_int)
     vort_int@description = "Relative vorticity"
     vort_int@units       = "s~S~-1~N~"

     copy_VarAtts(avo_plane,avo_int)
     avo_int@description = "Absolute vorticity"
     avo_int@units       = "s~S~-1~N~"

     copy_VarAtts(flux_plane,flux_int)
     flux_int@description = "Mass flux"
     flux_int@units       = "kg m~S~-2~N~ s~S~-1~N~" 

     copy_VarAtts(tht_plane,tht_int)
     tht_int@description = "Potential temperature"
     tht_int@units       = "K"

     copy_VarAtts(frcx_plane,frcx_int)
     frcx_int@description = "Diffusion (x-comp.)"
     frcx_int@units          = "m s~S~-1~N~"

     copy_VarAtts(frcy_plane,frcy_int)
     frcy_int@description = "Diffusion (y-comp.)"
     frcy_int@units       = "m s~S~-1~N~" 

     fric_int@description = "Diffusion (total)"
     fric_int@units       = "m s~S~-1~N~"

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

     flux_int!0    = "rad"
     flux_int!1    = "azi"
     flux_int&rad  = radii
     flux_int&azi  = thetas

     tht_int!0     = "rad"
     tht_int!1     = "azi"
     tht_int&rad   = radii
     tht_int&azi   = thetas

     v_tan!0      = "rad"
     v_tan!1      = "azi"
     v_tan&rad    = radii
     v_tan&azi    = thetas

     v_rad!0      = "rad"
     v_rad!1      = "azi"
     v_rad&rad    = radii
     v_rad&azi    = thetas

     fric_int!0     = "rad"
     fric_int!1     = "azi"
     fric_int&rad   = radii
     fric_int&azi   = thetas

     frcx_int!0     = "rad"
     frcx_int!1     = "azi"
     frcx_int&rad   = radii
     frcx_int&azi   = thetas

     frcy_int!0     = "rad"
     frcy_int!1     = "azi"
     frcy_int&rad   = radii
     frcy_int&azi   = thetas

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
      pv_int(:,iang)   = linint2_points(lon_sub,lat_sub,pv_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      tht_int(:,iang)   = linint2_points(lon1_sub,lat1_sub,tht_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      vort_int(:,iang) = linint2_points(lon_sub,lat_sub,vort_plane,False,\
					xpos(:,iang),ypos(:,iang),0)
      avo_int(:,iang)  = linint2_points(lon_sub,lat_sub,avo_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      flux_int(:,iang) = linint2_points(lon1_sub,lat1_sub,flux_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      frcx_int(:,iang) = linint2_points(lon_sub,lat_sub,frcx_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      frcy_int(:,iang) = linint2_points(lon_sub,lat_sub,frcy_plane,False,\
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
        fric_int(:,iang) = frcy_int(:,iang)*cos(theta) - frcx_int(:,iang)*sin(theta)
    end do

    ;===================================================
    ; Average variables azimuthally around the storm 
    ;===================================================

      do irad = 0, dimsizes(radii)-1
        r = radii(irad)
        vtan_azi(ilev,irad) = avg(v_tan(irad,:))         ; Tangential wind
	vrad_azi(ilev,irad) = avg(v_rad(irad,:))	 ; Radial wind 
	spd_azi(ilev,irad)  = avg(spd_int(irad,:))	 ; Vector velocity (for RMW)
	pvrt_azi(ilev,irad) = avg(pv_int(irad,:))	 ; PV 
	vvel_azi(ilev,irad) = avg(w_int(irad,:))	 ; Vertical velocity
	vort_azi(ilev,irad) = avg(vort_int(irad,:))	 ; Relative vorticity
	avo_azi(ilev,irad)  = avg(avo_int(irad,:))	 ; Absolute vorticity
	flux_azi(ilev,irad) = avg(flux_int(irad,:))	 ; Mass flux 
	tht_azi(ilev,irad)  = avg(tht_int(irad,:))	 ; Potential temperature
	fric_azi(ilev,irad) = avg(fric_int(irad,:))	 ; Vertical diffusion
      end do

    ;====================================
    ; Calculate perturbation variables 
    ;====================================

      do irad = 0, dimsizes(radii)-1

        ptb_vtan1(ilev,irad,:,ct) = v_tan(irad,:)    - vtan_azi(ilev,irad) ; [v']
        ptb_vvel1(ilev,irad,:,ct) = w_int(irad,:)    - vvel_azi(ilev,irad) ; [w']
        ptb_vrad0(ilev,irad,:)    = v_rad(irad,:)    - vrad_azi(ilev,irad) ; [u']
	ptb_vort0(ilev,irad,:)    = vort_int(irad,:) - vort_azi(ilev,irad) ; [vort']
        ptb_vflx0(ilev,irad,:)    = ptb_vrad0(ilev,irad,:) * \             ; [u' * vort']
                                    ptb_vort0(ilev,irad,:)        

      end do 

    ;====================================================================================
    ; Azimuthally average the perturbation variables ( {} = brackets in Persing et al.)
    ;====================================================================================

      do irad = 0, dimsizes(radii)-1

        ptb_vflx(ilev,irad) = avg(ptb_vflx0(ilev,irad,:))	       ; [ {u' * vort'} ]
	ptb_vort(ilev,irad) = avg(ptb_vort0(ilev,irad,:))	       ; [ {vort'} ]
	ptb_vrad(ilev,irad) = avg(ptb_vrad0(ilev,irad,:))	       ; [ {u'} ]

	test		    = sum(ptb_vrad0(ilev,irad,:)) / 73	       

      end do

    ;======================================================
    ; Also calculate radius of maximum vertical velocity
    ;======================================================
    
      omg_max = max(vvel_azi(ilev,:))
      dims_w  = dimsizes(vvel_azi(ilev,:))
      omg_1d  = ndtooned(vvel_azi(ilev,:))
      inds_w  = ind_resolve(maxind(omg_1d),dims_w)
      rmvv    = radii(inds_w(0,0))

      maxw_plot0(ilev,ct) = rmvv

    ;===================================
    ; Input values into larger arrays
    ;===================================

    ; Azimuthally-averaged variables
      vtan_plot0(ilev,:,ct)  = (/vtan_azi(ilev,:)/)      ; Tangential wind
      vrad_plot0(ilev,:,ct)  = (/vrad_azi(ilev,:)/)	 ; Radial wind
      pvrt_plot0(ilev,:,ct)  = (/pvrt_azi(ilev,:)/)	 ; PV
      vvel_plot0(ilev,:,ct)  = (/vvel_azi(ilev,:)/)	 ; Vertical velocity
      vort_plot0(ilev,:,ct)  = (/vort_azi(ilev,:)/)	 ; Relative vorticity
      avo_plot0(ilev,:,ct)   = (/avo_azi(ilev,:)/)	 ; Absolute vorticity
      flux_plot0(ilev,:,ct)  = (/flux_azi(ilev,:)/)	 ; Mass flux 
      tht_plot0(ilev,:,ct)   = (/tht_azi(ilev,:)/)       ; Potential temperature
      fric_plot0(ilev,:,ct)  = (/fric_azi(ilev,:)/)	 ; Vertical diffusion

    ; Perturbation variables 
      pvflx_plot0(ilev,:,ct) = (/ptb_vflx(ilev,:)/)	 ; Eddy vorticity flux (RHS #3)
      pvort_plot0(ilev,:,ct) = (/ptb_vort(ilev,:)/)	 ; Relative vorticity
      pvrad_plot0(ilev,:,ct) = (/ptb_vrad(ilev,:)/)	 ; Radial wind

    ; Tidy up
      delete([/u_int,v_int,w_int,pv_int,v_tan,v_rad/])
      delete([/flux_int,tht_int,fric_int,frcx_int,frcy_int/])
      delete([/flux_plane,u_plane,v_plane,w_plane,pv_plane,vort_plane,geo_plane/])
      delete([/lat_plane,lon_plane,lat1_plane,lon1_plane,frcx_plane,frcy_plane/])
      delete([/lat_sub,lat1_sub,lon_sub,lon1_sub/])

      end do			    ; End pressure level loop (do ilev = clev-1, clev+1)

    ct  = ct + 1		    ; Counter variable (time)
    ct0 = ct0 + 1		    ; Independent counter variable (time)

  ; Tidy up before next iteration
    delete([/year,month,day,hour,minute,second,utc_date/])
    delete([/lon,lat,lon1,lat1,time,times,u,v,z,t,th,w,vort,avo,f0,f1,pres,plevs/])
    delete([/lonF,latF,lon1F,lat1F,hyb,strx,stry,rho1,rho0,mflux,rho/])

   end do     ; End time loop (do it = 0, times-1)

;======================================================
; Calculate vertical derivatives and remaining terms
;======================================================

    ; Change units before calculations
      pres   = f->p
      plevs = dimsizes(pres)                ; Size of pressure level array
      pres@units = "hPa"                    ; Pressure level array units

      pres = pres * 100
      pres@units = "Pa"

      do it = ts0, tf0

      do ilev = 0, plevs-1

        plv = pres(ilev) / 100

        print_clock("Working on level: "+plv+" hPa")

        kp1 = min((/ilev+1,plevs-1/))
        km1 = max((/ilev-1,0/))

          do irad = 0, dimsizes(radii)-1

          ;=======================
          ; Calculate mean term
          ;=======================

            mvtdz_plot0(ilev,irad,it) = (vtan_plot0(kp1,irad,it) - \
                                         vtan_plot0(km1,irad,it)) / \
                                        (pres(km1) - pres(kp1))

          ;======================
          ; Calcuate eddy term
          ;======================

          ; Vertical derivative of eddy tangential wind
            pvtan_dz            = (ptb_vtan1(kp1,irad,:,it) - \
                                   ptb_vtan1(km1,irad,:,it)) / \
                                  (pres(km1) - pres(kp1))

          ; Vertical advection of eddy tangential momentum [-{w' * dv'/dz}]
            ptb_vadv1(ilev,irad,:,it) = ptb_vvel1(ilev,irad,:,it) * pvtan_dz

          end do                 ; End loop over radial circles (do irad = 0, ...)

          end do                  ; End loop over pressure levels (do ilev = 0, num_p-1)

       end do                    ; End loop over times (do it = 0, numTIMES-1)

     ; Return units to previous values
       pres = pres / 100
       pres@units = "hPa"

     ; Also calculate the azimuthal average of [-w' * dv'/dz]
       pvadv_plot0(:,:,:) = dim_avg_n(ptb_vadv1(:,:,:,:), 2)      ; [ {w' * dv'/dz} ]

;======================================
; Initialise new arrays for plotting 
;======================================

; Tangential and radial wind tendency
  vtan_tend0 = vtan_plot0		 ; Tangential wind tendency
  vrad_tend0 = vrad_plot0		 ; Radial wind tendency

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
  rhs1A	     = -vrad_plot0 * avo_plot0
  rhs1A	     = rhs1A * 3600
  rhs1A@units = "m s~S~-1~N~ h~S~-1~N~"

; Mean vert. adv. of mean tangential momentum
  rhs2A     = -vvel_plot0 * mvtdz_plot0
  rhs2A      = rhs2A * 3600
  rhs2A@units = "m s~S~-1~N~ h~S~-1~N~"

; Eddy radial vorticity flux
  rhs3A	     = -pvflx_plot0
  rhs3A	     = rhs3A * 3600
  rhs3A@units = "m s~S~-1~N~ h~S~-1~N~"

; Vertical advection of eddy tangential momentum
  rhs4A = -pvadv_plot0
  rhs4A = rhs4A * 3600
  rhs4A@units = "m s~S~-1~N~ h~S~-1~N~"

; Vertical diffusion 
  rhs6A = fric_plot0
  rhs6a = rhs6A * 3600
  rhs6A@units = "m s~S~-1~N~ h~S~-1~N~"

  vtan_tend0!0	  = "lev"
  vtan_tend0&lev  = pres 
  vtan_tend0!1    = "rad"
  vtan_tend0&rad  = radii
  vtan_tend0!2    = "time"
  vtan_tend0&time = time_int

  vrad_tend0!0	  = "lev"
  vrad_tend0&lev  = pres
  vrad_tend0!1    = "rad"
  vrad_tend0&rad  = radii
  vrad_tend0!2    = "time"
  vrad_tend0&time = time_int

  rhs1A!0	  = "lev"
  rhs1A&lev	  = pres
  rhs1A!1	  = "rad"
  rhs1A&rad	  = radii
  rhs1A!2	  = "time"
  rhs1A&time	  = time_int

  rhs2A!0         = "lev"
  rhs2A&lev       = pres
  rhs2A!1         = "rad"
  rhs2A&rad       = radii
  rhs2A!2         = "time"
  rhs2A&time      = time_int

  rhs3A!0	  = "lev"
  rhs3A&lev	  = pres
  rhs3A!1         = "rad"
  rhs3A&rad       = radii
  rhs3A!2         = "time"
  rhs3A&time      = time_int

  rhs4A!0         = "lev"
  rhs4A&lev       = pres
  rhs4A!1         = "rad"
  rhs4A&rad       = radii
  rhs4A!2         = "time"
  rhs4A&time      = time_int

  rhs6A!0         = "lev"
  rhs6A&lev       = pres
  rhs6A!1         = "rad"
  rhs6A&rad       = radii
  rhs6A!2         = "time"
  rhs6A&time      = time_int

  do it = ts0, tf0

   tp1  = min((/it+1,tf0/)) 
   tm1	= max((/it-1,ts0/))

   vtan_tend0(:,:,it) = (vtan_plot0(:,:,tp1) - \
   		      	 vtan_plot0(:,:,tm1)) / 2 ; Units [m s-1 h-1] 
   vrad_tend0(:,:,it) = (vrad_plot0(:,:,tp1) - \
   			 vrad_plot0(:,:,tm1)) / 2 ; Units [m s-1 h-1]

  end do     ; End time loop (do it = ts0, tf0)

;======================================
; Reorder dimensions before plotting 
;======================================

 ; Arrays now ordered [PRES * TIME * RADIUS] rather than [PRES * RADIUS * TIME]

 ; Tangential and radial wind tendencies
   vtan_tend = vtan_tend0(lev|:,time|:,rad|:)
   vrad_tend = vrad_tend0(lev|:,time|:,rad|:)

 ; Mean quantities
   mvtan_plot = vtan_plot0(lev|:,time|:,rad|:)
   mvrad_plot = vrad_plot0(lev|:,time|:,rad|:)
   mvvel_plot = vvel_plot0(lev|:,time|:,rad|:)
   mvort_plot = vort_plot0(lev|:,time|:,rad|:)
   mavo_plot  = avo_plot0(lev|:,time|:,rad|:)
   mflux_plot = flux_plot0(lev|:,time|:,rad|:)
   mfric_plot = fric_plot0(lev|:,time|:,rad|:)

 ; Perturbation quantities
   pvort_plot = pvort_plot0(lev|:,time|:,rad|:)
   pvrad_plot = pvrad_plot0(lev|:,time|:,rad|:)

 ; Momentum budget source terms
   rhs1	     = rhs1A(lev|:,time|:,rad|:)
   rhs2	     = rhs2A(lev|:,time|:,rad|:)
   rhs3	     = rhs3A(lev|:,time|:,rad|:)
   rhs4	     = rhs4A(lev|:,time|:,rad|:)
   rhs6	     = rhs6A(lev|:,time|:,rad|:)

;=========================================================
; Create a panel plot of terms from the momentum budget 
;=========================================================

  ; Output file location and type
    output = "$sam/nepartak/images/cyl_coords/mflux_hov_"+pres(clev)+"_"+ens0
    wks = gsn_open_wks(opt,output)

  ; Set up panel plot
    panel = new(3,graphic)
    gsn_define_colormap(wks,"bwr_new")

;==========================
; Options for plotting 
;==========================

  ; Mean terms
    opts_mean                              = True
    opts_mean@cnFillOn                     = True
    opts_mean@cnLineLabelInterval          = 2.0
    opts_mean@cnLineLabelFontHeightF       = 0.012
    opts_mean@cnLineLabelBackgroundColor   = "transparent"
    opts_mean@cnLineLabelPlacementMode     = "constant"
    opts_mean@cnLinesOn                    = False ; Contour lines off
    opts_mean@cnInfoLabelOn                = False
    opts_mean@cnLevelSelectionMode         = "ExplicitLevels"

    opts_mean@cnLevels                     = (/-2., -1., 1., 2., 5., 10., 20., \
    					        30., 40., 60., 80., 100./)
    opts_mean@cnFillColors                 = (/13,15,-1,18,19,21,23,24,25,28,29,30,31/)

    opts_mean@gsnPaperOrientation          = "landscape"
    opts_mean@tiMainString                 = ""
    opts_mean@tiMainFontHeightF            = 0.0125
    opts_mean@gsnLeftString                = ""
    opts_mean@gsnRightString               = ""

    if (lb.eq.0) then
     opts_mean@lbLabelBarOn                = False
     opts_mean@gsnDraw                     = False
     opts_mean@gsnFrame                    = False
    else if (lb.eq.1) then
     opts_mean@lbLabelBarOn                = True
     opts_mean@lbBoxEndCapStyle            = "TriangleBothEnds"    ; Labelbar end shape
     opts_mean@pmLabelBarWidthF		   = 0.32		   ; Labelbar width
     opts_mean@pmLabelBarHeightF	   = 0.08		   ; Labelbar height
     opts_mean@lbLabelFontHeightF          = 0.0175                ; Labelbar font size
     opts_mean@lbLabelFont                 = "Helvetica"           ; Labelbar font
     opts_mean@lbPerimOn                   = False                 ; Perimeter on/off
     opts_mean@pmLabelBarOrthogonalPosF    = 0.1                   ; Move lb up/down
     opts_mean@lbTitleString               = "m s~S~-1~N~"         ; Title string
     opts_mean@lbTitlePosition             = "Bottom"              ; Title position
     opts_mean@lbTitleFontHeightF          = 0.025                 ; Title font height
     opts_mean@gsnDraw                     = False                 ; Do not draw plot
     opts_mean@gsnFrame                    = False                 ; Do not advance frame
    end if
    end if

  ; Additional plotting resources
    opts_mean@tiYAxisString                = "Time (hour)"
    opts_mean@tiXAxisString		   = "Radius (km)"
    opts_mean@tiXAxisFontHeightF           = 0.018
    opts_mean@trYLog                       = False
    opts_mean@trYMinF                      = 30      ; Starting point (T+...) for plot
    opts_mean@trYMaxF                      = 72      ; Ending point (T+...) for plot
    opts_mean@tmXTOn                       = "False" ; Turn off top x-axis TM
    opts_mean@tmYROn                       = "False" ; Turn off right y-axis TM
    opts_mean@gsnMaximize                  = True    ; Maximise plot size
    opts_mean@gsnAddCyclic                 = False
    opts_mean@tmXBMode                     = "Explicit"                   ; Tm mode
    opts_mean@tmXBValues                   = (/0.0,0.25,0.50,0.75,1.0/)   ; Tm positions
    opts_mean@tmXBLabels                   = (/"0","25","50","75","100"/) ; Tm values
    opts_mean@tmXBMinorOn                  = False
    opts_mean@vpWidthF			   = 0.30
    opts_mean@vpHeightF			   = 0.50

;    opts_mean@pmLabelBarOrthogonalPosF     = 0.10    ; Move lb up/down (higher/lower)
;    opts_mean@tiXAxisOffsetYF		   = 0.095   ; Move x-axis title up/down (h/l)

    opts_mean@trXMinF			   = 0.0		       ; Set x-axis min
    opts_mean@trXMaxF			   = 1.0		       ; Set x-axis max
    opts_mean@tmYLMode			   = "Manual"
    opts_mean@tmYLTickStartF		   = opts_mean@trYMinF	       ; See L730
    opts_mean@tmYLTickEndF		   = opts_mean@trYMaxF	       ; See L731
    opts_mean@tmYLTickSpacingF		   = 6
    opts_mean@tmYLMinorOn		   = False

  ; Tangential wind resources 
    opts_vtan			           = opts_mean
    opts_vtan@cnFillPalette		   = "radar_new"
    delete(opts_vtan@cnLevels)
    delete(opts_vtan@cnFillColors)
    opts_vtan@cnLevels                     = (/5.0, 10.0, 15.0, 20.0,\
                                               25.0, 30.0, 35.0, 40.0, \
                                               45.0, 50.0, 55.0, 60.0, \
                                               65.0/)
    opts_vtan@cnFillColors 		   = (/0,2,3,4,5,6,7,8,9,10,11,12,13,14/)

  ; Vertical velocity resources
    opts_vvel                              = opts_vtan
    delete(opts_vvel@cnLevels)
    delete(opts_vvel@cnFillColors)
    opts_vvel@cnFillPalette		   = "n11"
    opts_vvel@cnLevels                     = (/0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, \
                                               0.7, 0.8, 1.0/)
    opts_vvel@cnFillColors                 = (/0,2,3,4,5,6,7,8,9,10,11/)
 
  ; Radial wind resources
    opts_vrad			           = opts_vvel
    delete(opts_vrad@cnLevels)
    delete(opts_vrad@cnFillColors)
    opts_vrad@cnFillPalette		   = "BlueDarkRed18"
    opts_vrad@cnLevels                     = (/-20., -15., -10., -5., -3., -2., -1., \
                                                 1., 2., 3., 5., 10./)
    opts_vrad@cnFillColors                 = (/3,4,5,6,7,8,9,10,11,13,14,15,16,17/)

  ; Momentum budget resources
    opts_mbgt				   = opts_vrad
    opts_mbgt@cnFillPalette		   = "bwr_new"
    delete(opts_mbgt@cnLevels)
    delete(opts_mbgt@cnFillColors)
    delete(opts_mbgt@lbTitleString)
    opts_mbgt@lbTitleString                = "m s~S~-1~N~ h~S~-1~N~"
    opts_mbgt@cnLevels                     = (/-40., -30., -25., -20., -15., \
    					       -10., -5., -3., -2., -1., \
                                                1., 2., 3., 5., 10., \
						15., 20., 25., 30., 40./)
    opts_mbgt@cnFillColors                 = (/3,4,5,6,9,10,11,12,13,14,-1,\
                                               17,19,20,21,22,24,26,27,28,29/)

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
    opts_zero                              = True
    opts_zero@cnFillOn                     = False
    opts_zero@cnLineColor                  = "gray82"
    opts_zero@cnLevels                     = (/-1.0, 0.0, 1.0/)
    opts_zero@cnLineThicknessF             = 2.0
    opts_zero@cnInfoLabelOn                = False
    opts_zero@gsnDraw                      = False       ; Do not draw the plot
    opts_zero@gsnFrame                     = False       ; Do no advance the frame
    opts_zero@gsnContourZeroLineThicknessF = 2.0
    opts_zero@cnLineLabelsOn               = False      ; Turn off line labels
    opts_zero@tiMainString                 = ""
    opts_zero@gsnLeftString                = ""
    opts_zero@gsnRightString               = ""

  ; Radius of maximum wind (RMW) resources 
    opts_rmw	     	       	           = True
    opts_rmw@xyLineColors      		   = "gray55"
    opts_rmw@xyLineThicknesses 		   = 3.0
    opts_rmw@xyDashPatterns    		   = 0.0
    opts_rmw@gsnDraw                       = False       ; Do not draw the plot
    opts_rmw@gsnFrame                      = False       ; Do not advance the frame
    opts_rmw@tiMainString                  = ""
    opts_rmw@gsnLeftString                 = ""
    opts_rmw@gsnRightString                = ""

  ; Tangential wind tendency resources 
    opts_tend                              = True
    opts_tend@cnLineColor                  = "black"
    opts_tend@cnFillOn			   = False
    opts_tend@cnLevelSelectionMode	   = "ExplicitLevels"
    opts_tend@cnLevels                     = (/2.0/)
    opts_tend@cnLineLabelsOn               = False
    opts_tend@cnInfoLabelOn		   = False
    opts_tend@cnLineThicknessF             = 3.5
    opts_tend@gsnContourNegLineDashPattern = 14
    opts_tend@gsnDraw			   = False
    opts_tend@gsnFrame			   = False
    opts_tend@tiMainString                 = ""
    opts_tend@gsnLeftString                = ""
    opts_tend@gsnRightString               = ""

  ; Mass flux resources
    opts_flux                             = opts_mean
    opts_flux@cnLinesOn                   = False ; Contour lines off
    opts_flux@cnInfoLabelOn               = False
    opts_flux@cnLevelSelectionMode        = "ExplicitLevels"
    opts_flux@cnFillPalette               = "radar_new"
    delete(opts_flux@cnLevels)
    delete(opts_flux@cnFillColors)
    delete(opts_flux@lbTitleString)
    opts_flux@lbTitleString               = "kg m~S~-2~N~ s~S~-1~N~"
    opts_flux@cnLevels                    = (/0.05,0.10,0.20,0.40,0.60,\
                                              0.80,1.00,1.50,2.00,3.00/)
    opts_flux@cnFillColors                = (/-1,5,6,7,8,9,10,11,12,13,14/)

;==============
; Panel plot 
;==============

    ilev = clev

  ; Change units for absolute/relative vorticity
    mavo_plot  = mavo_plot * (10 ^ 4)
    mvort_plot = mvort_plot * (10 ^ 4)

  ; Read in RMW data from text file
    spd_plot0 = asciiread("./"+ens0+"_rmw.txt",(/numTIMES/),"float")

  ; Panel 1 
    if (ptb.eq.0) then
     vrad0 = gsn_csm_contour(wks,mvrad_plot(ilev,:,:),opts_vrad) ; Mean vrad
    else if (ptb.eq.1) then
     vrad0 = gsn_csm_contour(wks,rhs1(ilev,:,:),opts_mbgt)	 ; RHS #1 
    end if 
    end if 
    rmw0  = gsn_csm_xy(wks,spd_plot0,time_int,opts_rmw)		 ; RMW 
    tend0 = gsn_csm_contour(wks,vtan_tend(ilev,:,:),opts_tend)   ; Vtan tendency
    overlay(vrad0,rmw0)						 ; Overlay RMW
    overlay(vrad0,tend0)					 ; Overlay vtan_tend
    panel(0) = vrad0						 ; Plot second panel
    delete(tend0)

  ; Panel 2 
    flux1 = gsn_csm_contour(wks,mflux_plot(ilev,:,:),opts_flux)  ; Mass flux
    rmw1  = gsn_csm_xy(wks,spd_plot0,time_int,opts_rmw)		 ; RMW
    overlay(flux1,rmw1)                                          ; Overlay RMW
    panel(1) = flux1                                             ; Plot second panel
    delete(flux1)

  ; Panel 3 
    if (ptb.eq.0) then
     vvel0 = gsn_csm_contour(wks,rhs1(ilev,:,:),opts_mean)	 ; RHS term 1
    else if (ptb.eq.1) then
     vvel0 = gsn_csm_contour(wks,rhs3(ilev,:,:),opts_mbgt)       ; RHS #3
    end if 
    end if 
    rmw2  = gsn_csm_xy(wks,spd_plot0,time_int,opts_rmw)		 ; RMW
    tend0 = gsn_csm_contour(wks,vtan_tend(ilev,:,:),opts_tend)   ; Vtan tendency
    overlay(vvel0,rmw2)						 ; Overlay RMW
    overlay(vvel0,tend0)					 ; Overlay vtan_tend
    panel(2) = vvel0						 ; Plot third panel
    delete(tend0)

    optsP                       = True
    optsP@gsnFrame              = False    ; Do not advance the frame

    if (lb.eq.0) then 
     optsP@gsnPanelLabelBar     = True                     ; Panel labelbar
    else if (lb.eq.1) then 
     optsP@gsnPanelLabelBar     = False
    end if 
    end if 

    optsP@pmLabelBarWidthF      = 0.6                      ; Label bar width
    optsP@pmLabelBarHeightF     = 0.15                     ; Label bar height
    optsP@lbLabelFontHeightF    = 0.01                     ; Label bar font height
    optsP@lbLabelFont           = "Helvetica"              ; Choose font
    optsP@lbPerimOn             = False                    ; Perimeter on label bar
    optsP@lbTitleString         = "m s~S~-1~N~ h~S~-1~N~"  ; Title string
    optsP@lbTitlePosition       = "Bottom"                 ; Title position
    optsP@lbTitleFontHeightF	= 0.015			   ; Title font height
    optsP@lbBoxEndCapStyle      = "TriangleBothEnds"       ; Label bar style

    optsP@txString              = "Simulation: "+ens0+" ("+pres(clev)+" hPa)"
    optsP@gsnPanelFigureStrings = (/"a) Mean vort. flux ","b) Mass flux",\
    				    "c) Eddy vort. flux "/)
    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.90
    optsP@gsnPanelBottom        = 0.10
    optsP@amJust                = "TopRight"
    optsP@gsnPanelFigureStringsFontHeightF = 0.0125 ; Reduce label size (default 0.01)

    gsn_panel(wks,panel,(/1,3/),optsP)             ; Draw as a single plot
    frame(wks)

end 