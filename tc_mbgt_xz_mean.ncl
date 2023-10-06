; Script to calculate the momentum budget as in Persing et al. (2013)

; Run using:

; ncl 'opt="x11"' 'ens0="em08"' clev=2 dist=3.0 nbox=1 ts0=0 tf0=72 'sys="lee"' 
; tc_mbgt_xz_mean.ncl

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"

;if (sy.eq."lee") then
;  path = "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"
;else if (sy.eq."ar") then
;  path = "$LIB/gsn_csm.ncl"
;end if 
;end if

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

 ; Read in values from chosen simulation
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

; Pressure levels (pres)
; (0) 1000, (1) 950, (2) 900, (3) 850, (4) 800, (5) 750, (6) 700, (7) 650
; (8) 600, (9) 550, (10) 500, (11) 450, (12) 400, (13) 350
; (14) 300, (15) 250, (16) 200, (17) 150, (18) 100

  num_p = 19      ; Size of pressure level array above (defined on L296)

;==================================================
; Define all other cylindrical coordinate arrays 
;==================================================

; Define 3D arrays for later use in tangential velocity cross-section plots
  vtan_azi = new((/nbox,num_p,dimsizes(radii)/),float)  ; Tangential wind      [ {v} ]
  vrad_azi = new((/nbox,num_p,dimsizes(radii)/),float)  ; Radial wind          [ {u} ]
  vvel_azi = new((/nbox,num_p,dimsizes(radii)/),float)  ; Vertical velocity    [ {w} ]
  vort_azi = new((/nbox,num_p,dimsizes(radii)/),float)  ; Rel. vorticity    [ {vort} ]
  avo_azi  = new((/nbox,num_p,dimsizes(radii)/),float)  ; Abs. vorticity    [{f + vort} ]
  tht_azi  = new((/nbox,num_p,dimsizes(radii)/),float)  ; Theta  	    [ {theta} ]
  frcx_azi = new((/nbox,num_p,dimsizes(radii)/),float)  ; Vert. diffusion (x)  [ {frcx} ]
  frcy_azi = new((/nbox,num_p,dimsizes(radii)/),float)  ; Vert. diffusion (y)  [ {frcy} ]
  fric_azi = new((/nbox,num_p,dimsizes(radii)/),float)  ; Vert. diffusion (tot)[ {fric} ]

;========================================================================================
;========================================================================================
; Calculating the terms from Persing et al. (2013) - their Eq.(12)
; LHS     : d{v}/dt	     - Mean tangential momentum tendency
; RHS (1) : -{u} * {f+vort}  - Mean radial influx of absolute vertical vorticity
; RHS (2) : -{w} * (d{v}/dz) - Mean vertical advection of mean tangential momentum
; RHS (3) : -{u' * vort'}    - Eddy radial vorticity flux 
; RHS (4) : -{w' * (dv'/dz)} - Vertical advection of eddy tangential momentum
; RHS (5) : c_p{(theta'/r) * (dpi'/dlambda)} - Azimuthal perturbation PGF per unit mass
; RHS (6) : {D_v}	     - Combined diffusive and PBL tendency
;========================================================================================
;========================================================================================

; Arrays for perturbation (eddy) variables (4D)
  ptb_vrad0 = new((/nbox,num_p,dimsizes(radii),dimsizes(thetas)/),float) ; [u']
  ptb_vort0 = new((/nbox,num_p,dimsizes(radii),dimsizes(thetas)/),float) ; [vort']
  ptb_vflx0 = new((/nbox,num_p,dimsizes(radii),dimsizes(thetas)/),float) ; -[u' * vort']
  ptb_vadv0 = new((/nbox,num_p,dimsizes(radii),dimsizes(thetas)/),float) ; -[w' * dv'/dz]

; 5D arrays for vertical derivative of perturbation tangential wind [-w' * dv'/dz]
  ptb_vvel1 = new((/nbox,num_p,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [w']
  ptb_vtan1 = new((/nbox,num_p,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [v']
  ptb_vadv1 = new((/nbox,num_p,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [R4]

; Arrays for azimuthally averaged eddy variables (3D)
  ptb_vrad  = new((/nbox,num_p,dimsizes(radii)/),float)       ; [ {u'} ]
  ptb_vtan  = new((/nbox,num_p,dimsizes(radii)/),float)       ; [ {v'} ]
  ptb_vvel  = new((/nbox,num_p,dimsizes(radii)/),float)       ; [ {w'} ]
  ptb_vort  = new((/nbox,num_p,dimsizes(radii)/),float)       ; [ {vort'} ]

  mvtan_dz  = new((/nbox,num_p,dimsizes(radii)/),float)       ; d{v}/dz  
  ptb_vflx  = new((/nbox,num_p,dimsizes(radii)/),float)       ; -{u' * vort'}    [R3]
  ptb_vadv  = new((/nbox,num_p,dimsizes(radii)/),float)       ; -{w' * (dv'/dz)} [R4]  

;================================================================
; Define 4D arrays (hold values from above arrays at ALL times)
;================================================================

; Mean tangential wind 
  vtan_plot0   = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)
  vtan_plot0!1 = "lev"
  vtan_plot0!2 = "rad"
  vtan_plot0!3 = "time"

; Mean radial wind 
  vrad_plot0   = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)
  vrad_plot0!1 = "lev"
  vrad_plot0!2 = "rad"
  vrad_plot0!3 = "time"

; Mean vertical velocity
  vvel_plot0   = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)
  vvel_plot0!1 = "lev"
  vvel_plot0!2 = "rad"
  vvel_plot0!3 = "time"

; Mean relative vorticity
  vort_plot0   = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)
  vort_plot0!1 = "lev"
  vort_plot0!2 = "rad"
  vort_plot0!3 = "time"

; Mean absolute vorticity
  avo_plot0    = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)
  avo_plot0!1  = "lev"
  avo_plot0!2  = "rad"
  avo_plot0!3  = "time"

; Mean potential temperature
  tht_plot0    = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)
  tht_plot0!1  = "lev"
  tht_plot0!2  = "rad"
  tht_plot0!3  = "time"

; Mean vertical diffusion (x-component)
  frcx_plot0   = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)
  frcx_plot0!1 = "lev"
  frcx_plot0!2 = "rad"
  frcx_plot0!3 = "time"

; Mean vertical diffusion (y-component)
  frcy_plot0   = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)
  frcy_plot0!1 = "lev"
  frcy_plot0!2 = "rad"
  frcy_plot0!3 = "time"

; Mean vertical diffusion (total)
  fric_plot0   = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)
  fric_plot0!1 = "lev"
  fric_plot0!2 = "rad"
  fric_plot0!3 = "time"

; Perturbation tangential wind
  pvtan_plot0   = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)
  pvtan_plot0!1 = "lev"
  pvtan_plot0!2 = "rad"
  pvtan_plot0!3 = "time"

; Perturbation radial wind
  pvrad_plot0   = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)
  pvrad_plot0!1 = "lev"
  pvrad_plot0!2 = "rad"
  pvrad_plot0!3 = "time"

; Perturbation vertical velocity
  pvvel_plot0   = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)
  pvvel_plot0!1 = "lev"
  pvvel_plot0!2 = "rad"
  pvvel_plot0!3 = "time"

; Perturbation relative vorticity
  pvort_plot0   = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)
  pvort_plot0!1 = "lev"
  pvort_plot0!2 = "rad"
  pvort_plot0!3 = "time"

; Vertical derivative of mean tangential wind
  mvtdz_plot0   = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)
  mvtdz_plot0!1 = "lev"
  mvtdz_plot0!2 = "rad"
  mvtdz_plot0!3 = "time"

; Vertical derivative of perturbation tangential wind
  pvtdz_plot0   = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)
  pvtdz_plot0!1 = "lev"
  pvtdz_plot0!2 = "rad"
  pvtdz_plot0!3 = "time"

; Eddy radial vorticity flux (RHS #3)
  pvflx_plot0   = new((/nbox,num_p,dimsizes(radii),numTIMES/),float)
  pvflx_plot0!1 = "lev"
  pvflx_plot0!2	= "rad"
  pvflx_plot0!3	= "time"

; Vertical advection of eddy tangential momentum (RHS #4)
  pvadv_plot0 	= new((/nbox,num_p,dimsizes(radii),numTIMES/),float)
  pvadv_plot0!1	= "lev"
  pvadv_plot0!2	= "rad"
  pvadv_plot0!3 = "time"

;===================================================
; Add metadata for azimuthally averaged variables	
;===================================================

  vtan_azi!1   	       = "lev"
  vtan_azi!2           = "rad"
  vtan_azi&rad         = radii
  vtan_azi@description = "Azimuthally averaged tangential velocity"
  vtan_azi@units       = "m s~S~-1~N~"

  vrad_azi!1	       = "lev"
  vrad_azi!2           = "rad"
  vrad_azi&rad         = radii
  vrad_azi@description = "Azimuthally averaged radial velocity"
  vrad_azi@units       = "m s~S~-1~N~"

  vort_azi!1	       = "lev"
  vort_azi!2           = "rad"
  vort_azi&rad         = radii
  vort_azi@description = "Azimuthally averaged relative vorticity"
  vort_azi@units       = "s~S~-1~N~"

  avo_azi!1	        = "lev"
  avo_azi!2             = "rad"
  avo_azi&rad           = radii
  avo_azi@description   = "Azimuthally averaged relative vorticity"
  avo_azi@units         = "s~S~-1~N~"

  vvel_azi!1		= "lev"
  vvel_azi!2            = "rad"
  vvel_azi&rad          = radii
  vvel_azi@description  = "Azimuthally averaged vertical velocity"
  vvel_azi@units        = "m s~S~-1~N~"

  tht_azi!1		= "lev"
  tht_azi!2             = "rad"
  tht_azi&rad           = radii
  tht_azi@description   = "Azimuthally averaged equivalent potential temperature"
  tht_azi@units         = "K"

  frcx_azi!1            = "lev"
  frcx_azi!2            = "rad"
  frcx_azi&rad          = radii
  frcx_azi@description  = "Azimuthally averaged friction (x-component)"
  frcx_azi@units        = "m s~S~-1~N~"

  frcy_azi!1            = "lev"
  frcy_azi!2            = "rad"
  frcy_azi&rad          = radii
  frcy_azi@description  = "Azimuthally averaged friction (y-component)"
  frcy_azi@units        = "m s~S~-1~N~"

  fric_azi!1            = "lev"
  fric_azi!2            = "rad"
  fric_azi&rad          = radii
  fric_azi@description  = "Azimuthally averaged friction (y-component)"
  fric_azi@units        = "m s~S~-1~N~"

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

    numTIMES    = 120                       ; Total number of times in all files [36]
    ct          = ts0                       ; Counter variable
    ct0		= 0			    ; Counter variable (independent of position)
    plevs	= 19			    ; Number of pressure levels
    llbox	= toint(dist*50)	    ; Size of lat/lon arrays (see line below)
    dsize	= (/llbox,llbox/)           ; Size of lat/lon arrays (depends on 'dist')
    d0		= dsize(0)
    d1		= dsize(1)

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

    time  = f->t(it)            ; Read in all times in file
    times = dimsizes(time)	; Define number of times in file

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
;   title_arr(ct) = "Valid at: "+time_str(ct0)
   itt = it + 1
   title_arr(it) = "Valid at "+time_str(ct0)+" (T+"+itt+")"

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

    pres  = f->p         	     ; Pressure (levels)
    plevs = dimsizes(pres)  	     ; Size of pressure level array
    pres@units = "hPa"

  ; Read in horizontal velocity, relative vorticity and geopotential height
    u       = f->u(it,:,{t0:t1},{n0:n1})        ; Zonal velocity (m/s)
    v       = f->v(it,:,{t0:t1},{n0:n1})        ; Meridional velocity (m/s)
    w       = f->dz_dt(it,:,{t0:t1},{n0:n1})    ; Vertical velocity (m/s)
    vort    = f->rvor(it,:,{t0:t1},{n0:n1})     ; Relative vorticity (/s)
    z       = f->ht(it,:,{t0:t1},{n0:n1})       ; Geopotential height (m)
    t	    = f->temp(it,:,{t0:t1},{n0:n1})	; Temperature (K)

  ; Calculate potential temperature (theta) using NCL function
    th      = pot_temp(pres,t,0,False)

  ; Convert to units we want
    z       = z/10
    z@units = "dam"			        ; Convert to decametres 
    vort    = vort * (10 ^ -6) 	 	        ; Convert to /s

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

;==============================================
; Get the variables we will need [pc stream]
;==============================================

  ; Get additional variables from the model-level file
    lon1F = g->longitude_1({n0:n1})            ; longitude [xxx grid points]
    lat1F = g->latitude_1({t0:t1})             ; latitude  [xxx grid points]
    lonF  = g->longitude({n0:n1})              ; longitude [xxx grid points]
    latF  = g->latitude({t0:t1})               ; latitude  [xxx grid points]
    hyb   = g->hybrid_ht(:) 		       ; 63 hybrid levels (5 m up to ~20 km)
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

    cp         = 0    ; Counter variable

  ; Create array of pressure values with same dimensions as 'fric_x_hgt'
    pres_fr    = conform_dims(dimsizes(fric_x_hgt), pres, 0)
    pres_fr = pres_fr * 100
    pres_fr@units = "Pa"

  ; TC latitude  -- centre(it,0)
  ; TC longitude -- centre(it,1)

;===============================
; Loop over pressure levels 
;===============================

    do ilev = 0, num_p-1 ; Loop over all pressure levels

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
     rho    = rho1(m0,0:d0-1,0:d1-1)
     rho_a = (1/rho) ; Specific volume (inverse density)

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

     print("Working on: " + pres(ilev) + " hPa")
     p = pres(ilev)

     u_plane       = u(ilev,:,:)
     v_plane       = v(ilev,:,:)
     w_plane       = w(ilev,:,:)
     vort_plane	   = vort(ilev,:,:)
     avo_plane	   = avo(ilev,:,:)
     geo_plane     = z(ilev,:,:)
     th_plane	   = th(ilev,:,:)

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

     lat_max 	   = t0 + (dy * inds(0,0))
     lon_max 	   = n0 + (dx * inds(0,1))

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

;=====================================================================================
; Create array of box positions either side of the original [Davis & Galarneau 2009]
;=====================================================================================

     if (nbox.eq.49) then
      lat_box = (/lat_max-(3*dy),lat_max-(2*dy),lat_max-(dy),lat_max,\
                  lat_max+(dy),lat_max+(2*dy),lat_max+(3*dy)/)
      lon_box = (/lon_max-(3*dx),lon_max-(2*dx),lon_max-(dx),lon_max,\
                  lon_max+(dx),lon_max+(2*dx),lon_max+(3*dx)/)
     else if (nbox.eq.25) then
      lat_box = (/lat_max-(2*dy),lat_max-(dy),lat_max,\
                  lat_max+(dy),lat_max+(2*dy)/)
      lon_box = (/lon_max-(2*dx),lon_max-(dx),lon_max,\
                  lon_max+(dx),lon_max+(2*dx)/)
     else if (nbox.eq.9) then
      lat_box = (/lat_max-(dy),lat_max,lat_max+(dy)/)
      lon_box = (/lon_max-(dx),lon_max,lon_max+(dx)/)
     else if (nbox.eq.1) then
      lat_box = (/lat_max/)
      lon_box = (/lon_max/)
     end if
     end if
     end if
     end if

     lsize   = dimsizes(lat_box) * dimsizes(lon_box)

;===========================================
; Switch to cylindrical coordinate system
;===========================================

; Create new arrays describing cylindrical coordinate system (2D for now)
  xpos = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float)
  ypos = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float)

; New arrays containing variables we want to plot
  u_int    = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float) ; Zonal velocity
  v_int    = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float) ; Meridional velocity
  w_int    = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float) ; Vertical velocity
  vort_int = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float) ; Relative vorticity
  avo_int  = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float) ; Absolute vorticity
  tht_int  = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float) ; Potential temperature
  v_rad    = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float) ; Radial velocity
  v_tan    = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float) ; Tangential velocity
  frcx_int = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float) ; Diffusion (x-comp.)
  frcy_int = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float) ; Diffusion (y-comp.)
  fric_int = new((/nbox,dimsizes(radii),dimsizes(thetas)/),float) ; Diffusion (total)

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

  copy_VarAtts(th_plane,tht_int)
  tht_int@description = "Potential temperature"
  tht_int@units       = "K"

  copy_VarAtts(vort_plane,vort_int)
  vort_int@description = "Relative vorticity"
  vort_int@units       = "s~S~-1~N~"

  copy_VarAtts(avo_plane,avo_int)
  avo_int@description = "Absolute vorticity"
  avo_int@units       = "s~S~-1~N~"

  copy_VarAtts(frcx_plane,frcx_int)
  frcx_int@description = "Diffusion (x-comp.)"
  frcx_int@units	  = "m s~S~-1~N~"

  copy_VarAtts(frcy_plane,frcy_int)
  frcy_int@description = "Diffusion (y-comp.)"
  frcy_int@units       = "m s~S~-1~N~"

  fric_int@description = "Diffusion (total)"
  fric_int@units       = "m s~S~-1~N~"

  v_rad@description = "Radial velocity"
  v_rad@units       = "m s~S~-1~N~"

  v_tan@description = "Tangential velocity"
  v_tan@units       = "m s~S~-1~N~"

;========================================
; Loop over array of lat/lon locations
;========================================

   do i = 0, dimsizes(lon_box)-1	; Longitude
    do j = 0, dimsizes(lat_box)-1	; Latitude     

    a = (dimsizes(lon_box) * i) + j
    print_clock("Working on box number "+a)

  ; Create cylindrincal coordinate array from original latitude/longitude grid
    do irad = 0, dimsizes(radii)-1
      r = radii(irad)
      do iang = 0, dimsizes(thetas)-1
        theta = thetas_rad(iang)

          xpos(a,irad,iang) = (/lon_box(i)/) + r*cos(theta)
	  ypos(a,irad,iang) = (/lat_box(j)/) + r*sin(theta)

      end do
    end do

  ; Fill new arrays with interpolated values on the cylindrical grid 
  ; Check to make sure that diffusion arrays are matched correctly w/ 'lat/lon_sub'

    do iang = 0, dimsizes(thetas_rad)-1
      u_int(a,:,iang)    = linint2_points(lon_sub,lat_sub,u_plane,False,\
                                          xpos(a,:,iang),ypos(a,:,iang),0)
      v_int(a,:,iang)    = linint2_points(lon_sub,lat_sub,v_plane,False,\
                                          xpos(a,:,iang),ypos(a,:,iang),0)
      w_int(a,:,iang)    = linint2_points(lon1_sub,lat1_sub,w_plane,False,\
                                          xpos(a,:,iang),ypos(a,:,iang),0)
      tht_int(a,:,iang)  = linint2_points(lon1_sub,lat1_sub,th_plane,False,\
                                          xpos(a,:,iang),ypos(a,:,iang),0)
      vort_int(a,:,iang) = linint2_points(lon_sub,lat_sub,vort_plane,False,\
                                          xpos(a,:,iang),ypos(a,:,iang),0)
      avo_int(a,:,iang)  = linint2_points(lon_sub,lat_sub,avo_plane,False,\
                                          xpos(a,:,iang),ypos(a,:,iang),0)
      frcx_int(a,:,iang) = linint2_points(lon_sub,lat_sub,frcx_plane,False,\
					  xpos(a,:,iang),ypos(a,:,iang),0)
      frcy_int(a,:,iang) = linint2_points(lon_sub,lat_sub,frcy_plane,False,\
					  xpos(a,:,iang),ypos(a,:,iang),0)
    end do

  ; Calculate radial/tangential windspeed and friction at each point
    do iang = 0, dimsizes(thetas)-1
      theta = thetas_rad(iang)
       v_rad(a,:,iang) = u_int(a,:,iang)*cos(theta) + v_int(a,:,iang)*sin(theta)
       v_tan(a,:,iang) = -u_int(a,:,iang)*sin(theta) + v_int(a,:,iang)*cos(theta)
       fric_int(a,:,iang) = frcy_int(a,:,iang)*cos(theta) - frcx_int(a,:,iang)*sin(theta)
    end do

     xpos!1       = "rad"
     xpos!2       = "azi"
     xpos&rad     = radii
     xpos&azi     = thetas

     ypos!1       = "rad"
     ypos!2       = "azi"
     ypos&rad     = radii
     ypos&azi     = thetas

     u_int!1      = "rad"
     u_int!2      = "azi"
     u_int&rad    = radii
     u_int&azi    = thetas

     v_int!1      = "rad"
     v_int!2      = "azi"
     v_int&rad    = radii
     v_int&azi    = thetas

     w_int!1      = "rad"
     w_int!2      = "azi"
     w_int&rad    = radii
     w_int&azi    = thetas

     tht_int!1    = "rad"
     tht_int!2    = "azi"
     tht_int&rad  = radii
     tht_int&azi  = thetas

     vort_int!1   = "rad"
     vort_int!2   = "azi"
     vort_int&rad = radii
     vort_int&azi = thetas

     avo_int!1    = "rad"
     avo_int!2    = "azi"
     avo_int&rad  = radii
     avo_int&azi  = thetas

     v_tan!1      = "rad"
     v_tan!2      = "azi"
     v_tan&rad    = radii
     v_tan&azi    = thetas

     v_rad!1      = "rad"
     v_rad!2      = "azi"
     v_rad&rad    = radii
     v_rad&azi    = thetas

     frcx_int!1   = "rad"
     frcx_int!2   = "azi"
     frcx_int&rad = radii
     frcx_int&azi = thetas

     frcy_int!1   = "rad"
     frcy_int!2   = "azi"
     frcy_int&rad = radii
     frcy_int&azi = thetas

     fric_int!1   = "rad"
     fric_int!2   = "azi"
     fric_int&rad = radii
     fric_int&azi = thetas

    ;============================================
    ; Calculate azimuthally-averaged variables 
    ;============================================

    ;;;;;;;;;; ADD METADATA TO '...AZI' VARIABLES ;;;;;;;;;;

      do irad = 0, dimsizes(radii)-1
        r = radii(irad)

        vtan_azi(a,ilev,irad) = avg(v_tan(a,irad,:))      ; [ {v} ]
	vrad_azi(a,ilev,irad) = avg(v_rad(a,irad,:))      ; [ {u} ]
	vvel_azi(a,ilev,irad) = avg(w_int(a,irad,:))      ; [ {w} ]
	vort_azi(a,ilev,irad) = avg(vort_int(a,irad,:))   ; [ {vort} ]
	avo_azi(a,ilev,irad)  = avg(avo_int(a,irad,:))    ; [ {f + vort} ]
	tht_azi(a,ilev,irad)  = avg(tht_int(a,irad,:))    ; [ {theta} ]
	fric_azi(a,ilev,irad) = avg(fric_int(a,irad,:))   ; [ {fric_xz + fric_yz} ]

      end do

;        vtan_all(a,ilev,irad,:) = v_tan(a,irad,:)
;        vrad_all(a,ilev,irad,:) = v_rad(a,irad,:)
;        vvel_all(a,ilev,irad,:) = w_int(a,irad,:)
;        vort_all(a,ilev,irad,:) = vort_int(a,irad,:)
;        avo_all(a,ilev,irad,:)  = avo_int(a,irad,:)
;        tht_all(a,ilev,irad,:)  = tht_int(a,irad,:)
;        fric_all(a,ilev,irad,:) = fric_int(a,irad,:)

    ;;;;; EDIT FROM HERE ;;;;;
    ;;;;; 'fric_int' arrays have 3 dimensions - (1) box * (2) radius * (3) azimuth

    ;====================================
    ; Calculate perturbation variables 
    ;====================================

      do irad = 0, dimsizes(radii)-1

        ptb_vtan1(a,ilev,irad,:,ct) = v_tan(a,irad,:) - vtan_azi(a,ilev,irad)    ; [v']
	ptb_vvel1(a,ilev,irad,:,ct) = w_int(a,irad,:) - vvel_azi(a,ilev,irad)    ; [w']	
	ptb_vrad0(a,ilev,irad,:)    = v_rad(a,irad,:) - vrad_azi(a,ilev,irad)    ; [u']
	ptb_vort0(a,ilev,irad,:)    = vort_int(a,irad,:) - vort_azi(a,ilev,irad) ; [vrt']

	ptb_vflx0(a,ilev,irad,:)    = ptb_vrad0(a,ilev,irad,:) * \
			       	      ptb_vort0(a,ilev,irad,:)              ; [u' * vrt']

      end do

    ;====================================================================================
    ; Azimuthally average the perturbation variables ( {} = brackets in Persing et al.)
    ;====================================================================================

      do irad = 0, dimsizes(radii)-1 

	ptb_vflx(a,ilev,irad) = avg(ptb_vflx0(a,ilev,irad,:))	      ; [ {u' * vort'} ]

      end do 
 
    ;===================================
    ; Input values into larger arrays
    ;===================================

    ; Azimuthally-averaged variables
      vtan_plot0(a,ilev,:,ct)  = vtan_azi(a,ilev,:)     ; Tangential wind       [v]
      vrad_plot0(a,ilev,:,ct)  = vrad_azi(a,ilev,:)     ; Radial wind 	        [u]
      vvel_plot0(a,ilev,:,ct)  = vvel_azi(a,ilev,:)     ; Vertical velocity     [w]
      vort_plot0(a,ilev,:,ct)  = vort_azi(a,ilev,:)     ; Rel. vorticity        [vort]
      avo_plot0(a,ilev,:,ct)   = avo_azi(a,ilev,:)      ; Abs. vorticity    [f + vort]
      tht_plot0(a,ilev,:,ct)   = tht_azi(a,ilev,:)      ; Potential temperature [theta] 
      fric_plot0(a,ilev,:,ct)  = fric_azi(a,ilev,:)     ; Vertical diffusion    [fx+fy]

    ; Perturbation variables 
      pvflx_plot0(a,ilev,:,ct) = (/ptb_vflx(a,ilev,:)/) ; Eddy vort. flux  [{u'*vort'}]

      end do  		            ; End latitude loop (do j = 0, ...)
      end do			    ; End longitude loop (do i = 0, ...)

    ; Tidy up
      delete([/u_int,v_int,w_int,vort_int,avo_int,tht_int,v_tan,v_rad/])
      delete([/u_plane,v_plane,w_plane,th_plane,vort_plane,avo_plane,geo_plane/])
      delete([/lat_plane,lat1_plane,lon_plane,lon1_plane,frcx_plane,frcy_plane/])
      delete([/lat_sub,lon_sub,lat1_sub,lon1_sub/])

     end do			    ; End pressure level loop (do ilev = 0, num_p-1)

    ct  = ct + 1		    ; Counter variable (time)
    ct0 = ct0 + 1		    ; Independent counter variable (time)

  ; Tidy up before next iteration
    delete([/year,month,day,hour,minute,second,utc_date/])
    delete([/lon,lat,lon1,lat1,time,times,u,v,z,t,th,w,vort,avo,f0,f1,pres,plevs/])
    delete([/lonF,latF,lon1F,lat1F,hyb,strx,stry,rho1/])

   end do     ; End time loop (do it = 0, times-1)

;======================================================
; Calculate vertical derivatives and remaining terms
;======================================================

    ; Change units before calculations
      pres  = f->p                     ; Pressure (levels)
      plevs = dimsizes(pres)           ; Size of pressure level array
      pres@units = "hPa"

      pres = pres * 100
      pres@units = "Pa"

     do a = 0, nbox-1 

      do it = ts0, tf0

        do ilev = 0, plevs-1

        kp1 = min((/ilev+1,plevs-1/))
        km1 = max((/ilev-1,0/))
 
          do irad = 0, dimsizes(radii)-1

          ;=======================
	  ; Calculate mean term
	  ;=======================

            mvtdz_plot0(a,ilev,irad,it) = (vtan_plot0(a,kp1,irad,it) - \
	    	    		  	   vtan_plot0(a,km1,irad,it)) / \
					  (pres(km1) - pres(kp1))  

          ;======================
	  ; Calcuate eddy term
	  ;======================

          ; Vertical derivative of eddy tangential wind 
            pvtan_dz            = (ptb_vtan1(a,kp1,irad,:,it) - \
	    			   ptb_vtan1(a,km1,irad,:,it)) / \
                                  (pres(km1) - pres(kp1))

          ; Vertical advection of eddy tangential momentum [-{w' * dv'/dz}] 
            ptb_vadv1(a,ilev,irad,:,it) = ptb_vvel1(a,ilev,irad,:,it) * pvtan_dz

          end do		 ; End loop over radial circles (do irad = 0, ...)

        end do 		         ; End loop over pressure levels (do ilev = 0, num_p-1)

       end do			 ; End loop over times (do it = 0, numTIMES-1)

      end do			 ; End loop over box ensemble members (do a = 0, nbox-1)

     ; Return units to previous values
       pres = pres / 100
       pres@units = "hPa"

     ; Also calculate the azimuthal average of [-w' * dv'/dz]  
       do a = 0, nbox-1
        pvadv_plot0(a,:,:,:) = dim_avg_n(ptb_vadv1(a,:,:,:,:), 2) ; [ {w' * dv'/dz} ]
       end do

;=========================================
; Average over all box ensemble members 
;=========================================

    ; Define new arrays before calculating averages 
      vrad_plot  = new((/num_p,dimsizes(radii),numTIMES/),float)
      vtan_plot  = new((/num_p,dimsizes(radii),numTIMES/),float)
      vvel_plot  = new((/num_p,dimsizes(radii),numTIMES/),float)
      vort_plot  = new((/num_p,dimsizes(radii),numTIMES/),float)
      avo_plot   = new((/num_p,dimsizes(radii),numTIMES/),float)
      tht_plot	 = new((/num_p,dimsizes(radii),numTIMES/),float)
      mvtdz_plot = new((/num_p,dimsizes(radii),numTIMES/),float)
      pvflx_plot = new((/num_p,dimsizes(radii),numTIMES/),float)
      pvadv_plot = new((/num_p,dimsizes(radii),numTIMES/),float)
      fric_plot	 = new((/num_p,dimsizes(radii),numTIMES/),float)

    ; Average over box ensemble members 
      vrad_plot  = dim_avg_n_Wrap(vrad_plot0(:,:,:,:),0)
      vtan_plot  = dim_avg_n_Wrap(vtan_plot0(:,:,:,:),0)
      vvel_plot  = dim_avg_n_Wrap(vvel_plot0(:,:,:,:),0)
      vort_plot  = dim_avg_n_Wrap(vort_plot0(:,:,:,:),0)
      avo_plot   = dim_avg_n_Wrap(avo_plot0(:,:,:,:),0)
      tht_plot	 = dim_avg_n_Wrap(tht_plot0(:,:,:,:),0)
      mvtdz_plot = dim_avg_n_Wrap(mvtdz_plot0(:,:,:,:),0)
      pvflx_plot = dim_avg_n_Wrap(pvflx_plot0(:,:,:,:),0)
      pvadv_plot = dim_avg_n_Wrap(pvadv_plot0(:,:,:,:),0)
      fric_plot  = dim_avg_n_Wrap(fric_plot0(:,:,:,:),0)

      vrad_plot!0    = "lev"
      vrad_plot&lev  = pres

      vtan_plot!0    = "lev"
      vtan_plot&lev  = pres

      vvel_plot!0    = "lev"
      vvel_plot&lev  = pres

      vort_plot!0    = "lev"
      vort_plot&lev  = pres

      avo_plot!0     = "lev"
      avo_plot&lev   = pres

      tht_plot!0     = "lev"
      tht_plot&lev   = pres

      mvtdz_plot!0   = "lev"
      mvtdz_plot&lev = pres

      pvflx_plot!0   = "lev"
      pvflx_plot&lev = pres

      pvadv_plot!0   = "lev"
      pvadv_plot&lev = pres

      fric_plot!0    = "lev"
      fric_plot&lev  = pres

;=================================================
; Tidy up and produce final arrays for plotting 
;=================================================

     ; Calculate terms from Persing et al. (2013) - their Eq.(12)
     ; LHS     : d{v}/dt          - Rate of change of mean tangential wind 
     ; RHS #1  : -{u} * {f+vort}  - Mean radial influx of absolute vertical vorticity
     ; RHS #2  : -{w} * (d{v}/dz) - Mean vert. advection of mean tangential momentum
     ; RHS #3  : -{u' * vort'}    - Eddy radial vorticity flux
     ; RHS #4  : -{w' * (dv'/dz)} - Vertical advection of eddy tangential momentum
     ; RHS #5  : c_p{(theta'/r) * (dpi'/dlambda)} - Azimuthal pert. PGF per unit mass
     ; RHS #6  : {D_v1}           - Combined diffusive and PBL tendency (vertical)
     ; RHS #7  ; {D_v2}		  - Combined diffusive and PBL tendency (horizontal)

     ; Define new arrays for terms in momentum budget 
       lhs  = new((/num_p,dimsizes(radii),numTIMES/),float)   ; d{v}/dt
       rhs1 = new((/num_p,dimsizes(radii),numTIMES/),float)   ; -{u} * {f+vort}
       rhs2 = new((/num_p,dimsizes(radii),numTIMES/),float)   ; -{w} * (d{v}/dz)
       rhs3 = new((/num_p,dimsizes(radii),numTIMES/),float)   ; -{u' * vort'} 
       rhs4 = new((/num_p,dimsizes(radii),numTIMES/),float)   ; -{w' * (dv'/dz)}
       rhs5 = new((/num_p,dimsizes(radii),numTIMES/),float)   ; PGF term
       rhs6 = new((/num_p,dimsizes(radii),numTIMES/),float)   ; {D_v1}
       rhs7 = new((/num_p,dimsizes(radii),numTIMES/),float)   ; {D_v2}
       rhsT = new((/num_p,dimsizes(radii),numTIMES/),float)   ; Sum of rhs1-7
       res  = new((/num_p,dimsizes(radii),numTIMES/),float)   ; Residual (LHS - RHS)

       if (ts0.eq.0) then
         ts = 0
         tf = tf0-1
       else if (ts0.ge.1) then
         ts = ts0+1
         tf = tf0-1
       end if
       end if

     ; Calculate mean tangential momentum tendency (LHS)
       do it = ts, tf

         cp1 = min((/it+1,tf/))
	 cm1 = max((/it-1,ts/))
 
         lhs(:,:,it) = (vtan_plot(:,:,cp1) - vtan_plot(:,:,cm1)) / 2

       end do	; End loop over times (do it = 0, numTIMES-1)

     ; Add metadata and make units consistent for all terms
       lhs@units = "m s~S~-1~N~ h~S~-1~N~"

       rhs1 = -vrad_plot0 * avo_plot0   ; Mean radial influx of abs. vert. vorticity
       rhs1 = rhs1 * 3600
       rhs1@units = "m s~S~-1~N~ h~S~-1~N~"

       rhs2 = -vvel_plot0 * mvtdz_plot0 ; Mean vert. adv. of mean tangential momentum
       rhs2 = rhs2 * 3600
       rhs2@units = "m s~S~-1~N~ h~S~-1~N~"
       
       rhs3 = -pvflx_plot0  	        ; Eddy radial vorticity flux 
       rhs3 = rhs3 * 3600 
       rhs3@units = "m s~S~-1~N~ h~S~-1~N~"

       rhs4 = -pvadv_plot0	        ; Vertical advection of eddy tangential momentum
       rhs4 = rhs4 * 3600
       rhs4@units = "m s~S~-1~N~ h~S~-1~N~"

       rhs6 = fric_plot0		; Vertical diffusion
       rhs6 = rhs6 * 3600
       rhs6@units = "m s~S~-1~N~ h~S~-1~N~"

       rhsT = rhs1 + rhs2 + rhs3 + rhs4 + rhs6   ; Sum of RHS terms (add rhs5 & rhs7)
       rhsT@units = "m s~S~-1~N~ h~S~-1~N~"

       res  = lhs - rhsT			 ; Residual (LHS - RHS)
       res@units = "m s~S~-1~N~ h~S~-1~N~"

       lhs!0   = "lev"
       lhs&lev = pres
       lhs!1   = "rad"
       lhs&rad = radii
       lhs!2   = "time"

       rhs1!0   = "lev"
       rhs1&lev = pres
       rhs1!1   = "rad"
       rhs1&rad = radii
       rhs1!2   = "time"

       rhs2!0   = "lev"
       rhs2&lev = pres
       rhs2!1   = "rad"
       rhs2&rad = radii
       rhs2!2   = "time"

       rhs3!0   = "lev"
       rhs3&lev = pres
       rhs3!1   = "rad"
       rhs3&rad = radii
       rhs3!2   = "time"

       rhs4!0   = "lev"
       rhs4&lev = pres
       rhs4!1   = "rad"
       rhs4&rad = radii
       rhs4!2   = "time"

       rhs6!0   = "lev"
       rhs6&lev = pres
       rhs6!1   = "rad"
       rhs6&rad = radii
       rhs6!2   = "time"

       rhsT!0   = "lev"
       rhsT&lev = pres
       rhsT!1   = "rad"
       rhsT&rad = radii
       rhsT!2   = "time"

       res!0    = "lev"
       res&lev  = pres
       res!1    = "rad"
       res&rad  = radii
       res!2    = "time"

;============================================
; Define output file path and plot details 
;============================================

 if (ts0.eq.0) then
  ts = 0
  tf = tf0-1
 else if (ts0.ge.1) then
  ts = ts0+1
  tf = tf0-1
 end if 
 end if 

 do it = ts, tf

    print("Working on time: "+it)

  ; Output file location and type
    output = "$sam/nepartak/images/cyl_coords/mbgt_xz_mean_"+ens0+\
    	     "_box"+nbox+"_"+time_arr(it)
    wks = gsn_open_wks(opt,output)

  ; Set up panel plot (6 panels)
    panel = new(6,graphic)
    gsn_define_colormap(wks,"BlueDarkRed18")

;==========================
; Options for plotting 
;==========================

  ; Momentum budget terms 
    opts_mbgt                             = True
    opts_mbgt@cnFillOn                    = True
    opts_mbgt@cnLineLabelInterval         = 2.0
    opts_mbgt@cnLineLabelFontHeightF      = 0.012
    opts_mbgt@cnLineLabelBackgroundColor  = "transparent"
    opts_mbgt@cnLineLabelPlacementMode    = "constant"
    opts_mbgt@cnLinesOn                   = False ; Contour lines off
    opts_mbgt@cnInfoLabelOn               = False
    opts_mbgt@cnLevelSelectionMode        = "ExplicitLevels"
    opts_mbgt@cnLevels                     = (/-15., -10., -5., -3., -2., \
    					       -1., -0.5, -0.2, \
                                                0, 0.2, 0.5, 1., \
						2., 3., 5., 10., 15./)
    opts_mbgt@cnFillColors                 = (/2,3,4,5,6,7,8,9,0,0,\
    					       11,12,13,14,15,16,17,19/)
    opts_mbgt@gsnPaperOrientation         = "landscape"
    opts_mbgt@tiMainString                = ""
    opts_mbgt@tiMainFontHeightF           = 0.0125
    opts_mbgt@gsnLeftString               = ""
    opts_mbgt@gsnRightString              = ""
    opts_mbgt@gsnMaximize                 = True
    opts_mbgt@lbLabelBarOn                = True
    opts_mbgt@lbBoxEndCapStyle            = "TriangleBothEnds"    ; Labelbar end shape
;    opts_mbgt@pmLabelBarWidthF            = 0.32                  ; Labelbar width
;    opts_mbgt@pmLabelBarHeightF           = 0.08                  ; Labelbar height
    opts_mbgt@lbLabelFontHeightF          = 0.0125                ; Labelbar font size
    opts_mbgt@lbLabelFont                 = "Helvetica"           ; Labelbar font
    opts_mbgt@lbPerimOn                   = False                 ; Perimeter on/off
    opts_mbgt@gsnDraw                     = False                 ; Do not draw plot
    opts_mbgt@gsnFrame                    = False                 ; Do not advance frame
;    opts_mbgt@vpWidthF                    = 0.30                  ; Match width in N11
;    opts_mbgt@vpHeightF                   = 0.50                  ; Match height in N11
    opts_mbgt@gsnDraw                     = False
    opts_mbgt@gsnFrame                    = False

  ; Additional plotting resources
    opts_mbgt@tiYAxisString                = "Pressure (hPa)"
    opts_mbgt@tiXAxisString                = "Radial distance from cyclone centre (km)"
    opts_mbgt@tiXAxisFontHeightF           = 0.018
    opts_mbgt@trYLog                       = False

    opts_mbgt@tmXTOn                       = "False" ; Turn off top x-axis TM
    opts_mbgt@tmYROn                       = "False" ; Turn off right y-axis TM
    opts_mbgt@gsnMaximize                  = True    ; Maximise plot size
    opts_mbgt@gsnAddCyclic                 = False
    opts_mbgt@tmXBMode                     = "Explicit"  ; Set tick marks explicitly
    opts_mbgt@tmXBValues                   = (/0,0.5,1.0,1.5/) ; Tick mark positions
    opts_mbgt@tmXBLabels                   = (/"0","50","100","150"/) ; Tick mark values
    opts_mbgt@tmXBMinorValues              = fspan(0.0, 3.0, 16)      ; Minor tick marks

  ; Resources to overlay tangential wind (as rough marker of eyewall)
    opts_cont				  = True
    opts_cont@cnFillOn			  = False
    opts_cont@cnLineColor		  = "black"
    opts_cont@cnLevelSelectionMode	  = "ExplicitLevels"
    opts_cont@cnLevels			  = (/20,25,30,35,40,45,50,55/)
    opts_cont@cnInfoLabelOn		  = False
    opts_cont@cnLineThicknessF		  = 2.0
    opts_cont@gsnDraw			  = False
    opts_cont@gsnFrame			  = False
    opts_cont@cnLineLabelsOn		  = False
    opts_cont@tiMainString                = ""
    opts_cont@gsnLeftString               = ""
    opts_cont@gsnRightString              = ""

  ; Tangential wind resources
    opts_vtan                             = opts_mbgt
    opts_vtan@cnFillPalette               = "radar_new"
    delete(opts_vtan@cnLevels)
    delete(opts_vtan@cnFillColors)
    opts_vtan@cnLevels                    = (/5.0, 10.0, 15.0, 20.0,\
                                              25.0, 30.0, 35.0, 40.0, \
                                              45.0, 50.0, 55.0, 60.0, \
                                              65.0/)
    opts_vtan@cnFillColors                = (/0,2,3,4,5,6,7,8,9,10,11,12,13,14/)

  ; Radial wind resources
    opts_vrad                             = opts_mbgt
    delete(opts_vrad@cnLevels)
    delete(opts_vrad@cnFillColors)
    opts_vrad@cnFillPalette               = "BlueDarkRed18"
    opts_vrad@cnLevels                    = (/-20., -15., -10., -5., -2., -1., \
                                               1., 2., 5., 10., 15., 20./)
    opts_vrad@cnFillColors                = (/3,4,6,7,8,9,-1,11,12,13,14,15,16/)

  ; Vertical velocity resources
    opts_vvel                             = opts_mbgt
    delete(opts_vvel@cnLevels)
    delete(opts_vvel@cnFillColors)
    opts_vvel@cnFillPalette               = "n11"
    opts_vvel@cnLevels                    = (/0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, \
                                              0.7, 0.8, 1.0/)
    opts_vvel@cnFillColors                = (/0,2,3,4,5,6,7,8,9,10,11/)

  ; Vorticity resources
    opts_vort			          = opts_mbgt 
    delete(opts_vort@cnLevels)
    delete(opts_vort@cnFillColors)
    opts_vort@cnFillPalette	          = "circ_budget"
    opts_vort@cnLevels                    = (/-2., 2., 5., 10., 20., 30., \
                                               40., 60., 80., 100., 150., 200./)
    opts_vort@cnFillColors                = (/6,-1,10,11,13,14,15,16,17,18,19,20,21/)

  ; Potential temperature resources
    opts_thta                              = opts_vtan

    ilev      = clev

  ; Determine contour levels for theta using crude method
    th_mean0  = (max(tht_plot(:,:,it)) + min(tht_plot(:,:,it))) / 2
    th_mean   = toint(th_mean0)
    cnLevs    = ispan(th_mean-24,th_mean+24,4)

    delete(opts_thta@cnLevels)
    delete(opts_thta@cnFillColors)
    opts_thta@cnLevels                     = cnLevs
    opts_thta@cnFillColors                 = (/0,2,3,4,5,6,7,8,9,10,11,12,13,14/)

;=================
; Plot the data 
;=================

  ; Change units for vorticity and vertical diffusion
    vort_plot = vort_plot * (10 ^ 4)
    avo_plot  = avo_plot * (10 ^ 4)
    fric_plot = fric_plot * (10 ^ 4)

    plot0      = gsn_csm_pres_hgt(wks,vtan_plot(:,:,it),opts_vtan)  ; Tangential wind

    plot1      = gsn_csm_pres_hgt(wks,vrad_plot(:,:,it),opts_mbgt)  ; Radial wind 
    vtan_plotF = gsn_csm_contour(wks,vtan_plot(:,:,it),opts_cont)   ; Tangential wind
    overlay(plot1,vtan_plotF)                                       ; Overlay
    delete(vtan_plotF)

    plot2      = gsn_csm_pres_hgt(wks,vvel_plot(:,:,it),opts_vvel)  ; Vertical velocity
    vtan_plotF = gsn_csm_contour(wks,vtan_plot(:,:,it),opts_cont)   ; Tangential wind
    overlay(plot2,vtan_plotF)                                       ; Overlay
    delete(vtan_plotF)

    plot3      = gsn_csm_pres_hgt(wks,vort_plot(:,:,it),opts_vort)  ; Relative vorticity
    vtan_plotF = gsn_csm_contour(wks,vtan_plot(:,:,it),opts_cont)   ; Tangential wind
    overlay(plot3,vtan_plotF)                                       ; Overlay
    delete(vtan_plotF)

    plot4      = gsn_csm_pres_hgt(wks,avo_plot(:,:,it),opts_vort)   ; Absolute vorticity
    vtan_plotF = gsn_csm_contour(wks,vtan_plot(:,:,it),opts_cont)   ; Tangential wind
    overlay(plot4,vtan_plotF)                                       ; Overlay
    delete(vtan_plotF)

    plot5      = gsn_csm_pres_hgt(wks,tht_plot(:,:,it),opts_thta)   ; Potential temp
    vtan_plotF = gsn_csm_contour(wks,vtan_plot(:,:,it),opts_cont)   ; Tangential wind
    overlay(plot5,vtan_plotF)                                       ; Overlay
    delete(vtan_plotF)

    plot6      = gsn_csm_pres_hgt(wks,fric_plot(:,:,it),opts_mbgt)  ; Vertical diffusion
    vtan_plotF = gsn_csm_contour(wks,vtan_plot(:,:,it),opts_cont)   ; Tangential wind
    overlay(plot6,vtan_plotF)                                  	    ; Overlay
    delete(vtan_plotF)

    panel(0) = plot0
    panel(1) = plot1
    panel(2) = plot2
    panel(3) = plot4
    panel(4) = plot5
    panel(5) = plot6

  ; Change units back before next iteration
    vort_plot = vort_plot / (10 ^ 4)
    avo_plot  = avo_plot / (10 ^ 4)
    fric_plot = fric_plot / (10 ^ 4)

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    optsP                       = True
    optsP@gsnFrame              = False    ; Do not advance the frame
    optsP@gsnPanelLabelBar      = False    ; Panel labelbar

;    optsP@pmLabelBarWidthF	= 0.6	   ; Label bar width 
;    optsP@pmLabelBarHeightF     = 0.15	   ; Label bar height
;    optsP@lbLabelFontHeightF	= 0.0125   ; Label bar font height 
;    optsP@lbLabelFont		= "Helvetica"
;    optsP@lbPerimOn		= False
;    optsP@lbBoxEndCapStyle	= "TriangleBothEnds"

    optsP@txString              = title_arr(it)
    optsP@gsnPanelFigureStrings = (/"a) vtan", "b) vrad", "c) vvel", \
                                    "d) avo", "e) theta", "f) diff "/)
    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.90
    optsP@gsnPanelBottom        = 0.10
    optsP@amJust                = "TopLeft"
    optsP@gsnPanelFigureStringsFontHeightF = 0.0125 ; Label size (default 0.01)

    gsn_panel(wks,panel,(/2,3/),optsP)             ; Draw as a single plot
    frame(wks)

   end do     ; End time loop (do it = 0, ts, tf)

end 