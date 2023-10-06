; Script to reproduce Figure 6 from Nguyen et al. (2011) for a given ensemble member
; Run using:

; ncl 'opt="x11"' 'ens0="em02"' clev=2 clev0=2 dist=3.0 tc_ens_mbgt_fig6_test.ncl

; 'opt' = output file format; 'ens0' = ensemble member; 'clev' = pressure level (vtan)
; 'clev0' = pressure level (omega); 'dist' = size of box following storm (degrees)

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
  radii      = fspan(0.0,1.0,21)     ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,73)   ; Array of azimuth angles (0, 5, 10, etc)
  thetas_rad = (pii/180)*thetas      ; Degrees to radians
  rad_size   = dimsizes(radii)-1     ; Index for later calculation

; Pressure levels (pres)
; (0) 1000, (1) 950, (2) 900, (3) 850, (4) 800, (5) 750, (6) 700, (7) 650
; (8) 600, (9) 550, (10) 500, (11) 450, (12) 400, (13) 350
; (14) 300, (15) 250, (16) 200, (17) 150, (18) 100

  num_p = 19      ; Size of pressure level array above (defined on L296)

;=====================================================
; Define arrays for later use in radius-time plots 
;=====================================================

  vtan_azi = new((/num_p,dimsizes(radii)/),float) ; Tangential wind
  vrad_azi = new((/num_p,dimsizes(radii)/),float) ; Radial wind
  spd_azi  = new((/num_p,dimsizes(radii)/),float) ; Vector wind
  vvel_azi = new((/num_p,dimsizes(radii)/),float) ; Vertical velocity
  vort_azi = new((/num_p,dimsizes(radii)/),float) ; Relative vorticity
  avo_azi  = new((/num_p,dimsizes(radii)/),float) ; Absolute vorticity
  th_azi   = new((/num_p,dimsizes(radii)/),float) ; Theta
  frcx_azi = new((/num_p,dimsizes(radii)/),float) ; Vert. diffusion (x)
  frcy_azi = new((/num_p,dimsizes(radii)/),float) ; Vert. diffusion (y)
  fric_azi = new((/num_p,dimsizes(radii)/),float) ; Vert. diffusion (tot)

;========================================================================================
;========================================================================================
; Calculating the terms from Persing et al. (2013) - their Eq.(12)
; LHS     : d{v}/dt          - Mean tangential momentum tendency
; RHS (1) : -{u} * {f+vort}  - Mean radial influx of absolute vertical vorticity
; RHS (2) : -{w} * (d{v}/dz) - Mean vertical advection of mean tangential momentum
; RHS (3) : -{u' * vort'}    - Eddy radial vorticity flux
; RHS (4) : -{w' * (dv'/dz)} - Vertical advection of eddy tangential momentum
; RHS (5) : c_p{(theta'/r) * (dpi'/dlambda)} - Azimuthal perturbation PGF per unit mass
; RHS (6) : {D_v}            - Combined diffusive and PBL tendency
;========================================================================================
;========================================================================================

; Arrays for perturbation (eddy) variables (3D)
  ptb_vrad0 = new((/num_p,dimsizes(radii),dimsizes(thetas)/),float) ; [u']
  ptb_vort0 = new((/num_p,dimsizes(radii),dimsizes(thetas)/),float) ; [vort']
  ptb_vflx0 = new((/num_p,dimsizes(radii),dimsizes(thetas)/),float) ; -[u' * vort']
  ptb_vadv0 = new((/num_p,dimsizes(radii),dimsizes(thetas)/),float) ; -[w' * dv'/dz]

; 4D arrays for vertical derivative of perturbation tangential wind [-w' * dv'/dz]
  ptb_vvel1 = new((/num_p,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [w']
  ptb_vtan1 = new((/num_p,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [v']
  ptb_vadv1 = new((/num_p,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [R4]

; Arrays for azimuthally averaged eddy variables (1D)
  ptb_vrad  = new((/num_p,dimsizes(radii)/),float)       ; [ {u'} ]
  ptb_vtan  = new((/num_p,dimsizes(radii)/),float)       ; [ {v'} ]
  ptb_vvel  = new((/num_p,dimsizes(radii)/),float)       ; [ {w'} ]
  ptb_vort  = new((/num_p,dimsizes(radii)/),float)       ; [ {vort'} ]

  mvtan_dz  = new((/num_p,dimsizes(radii)/),float)       ; d{v}/dz
  ptb_vflx  = new((/num_p,dimsizes(radii)/),float)       ; -{u' * vort'}    [R3]
  ptb_vadv  = new((/num_p,dimsizes(radii)/),float)       ; -{w' * (dv'/dz)} [R4]

;==================================================================
; Define 3D arrays (hold values from above arrays at ALL times)
;==================================================================

; Mean tangential wind 
  vtan_plot0   = new((/num_p,dimsizes(radii),numTIMES/),float)
  vtan_plot0!1 = "rad"
  vtan_plot0!2 = "time"

; Mean radial wind 
  vrad_plot0   = new((/num_p,dimsizes(radii),numTIMES/),float)
  vrad_plot0!1 = "rad"
  vrad_plot0!2 = "time"

; Mean vector wind 
  spd_plot     = new((/num_p,numTIMES/),float)
  spd_plot!1   = "time"

; Mean vertical velocity
  vvel_plot0   = new((/num_p,dimsizes(radii),numTIMES/),float)
  vvel_plot0!1 = "rad"
  vvel_plot0!2 = "time"

; Mean relative vorticity
  vort_plot0   = new((/num_p,dimsizes(radii),numTIMES/),float)
  vort_plot0!1 = "rad"
  vort_plot0!2 = "time"

; Mean absolute vorticity
  avo_plot0    = new((/num_p,dimsizes(radii),numTIMES/),float)
  avo_plot0!1  = "rad"
  avo_plot0!2  = "time"

; Mean potential temperature
  th_plot0     = new((/num_p,dimsizes(radii),numTIMES/),float)
  th_plot0!1   = "rad"
  th_plot0!2   = "time"

; Mean vertical diffusion (x-component)
  frcx_plot0   = new((/num_p,dimsizes(radii),numTIMES/),float)
  frcx_plot0!1 = "rad"
  frcx_plot0!2 = "time"

; Mean vertical diffusion (y-component)
  frcy_plot0   = new((/num_p,dimsizes(radii),numTIMES/),float)
  frcy_plot0!1 = "rad"
  frcy_plot0!2 = "time"

; Mean vertical diffusion (total)
  fric_plot0   = new((/num_p,dimsizes(radii),numTIMES/),float)
  fric_plot0!1 = "rad"
  fric_plot0!2 = "time"

; Perturbation tangential wind
  pvtan_plot0   = new((/num_p,dimsizes(radii),numTIMES/),float)
  pvtan_plot0!1 = "rad"
  pvtan_plot0!2 = "time"

; Perturbation radial wind
  pvrad_plot0   = new((/num_p,dimsizes(radii),numTIMES/),float)
  pvrad_plot0!1 = "rad"
  pvrad_plot0!2 = "time"

; Perturbation vertical velocity
  pvvel_plot0   = new((/num_p,dimsizes(radii),numTIMES/),float)
  pvvel_plot0!1 = "rad"
  pvvel_plot0!2 = "time"

; Perturbation relative vorticity
  pvort_plot0   = new((/num_p,dimsizes(radii),numTIMES/),float)
  pvort_plot0!1 = "rad"
  pvort_plot0!2 = "time"

; Vertical derivative of mean tangential wind
  mvtdz_plot0   = new((/num_p,dimsizes(radii),numTIMES/),float)
  mvtdz_plot0!1 = "rad"
  mvtdz_plot0!2 = "time"

; Vertical derivative of perturbation tangential wind
  pvtdz_plot0   = new((/num_p,dimsizes(radii),numTIMES/),float)
  pvtdz_plot0!1 = "rad"
  pvtdz_plot0!2 = "time"

; Eddy radial vorticity flux (RHS #3)
  pvflx_plot0   = new((/num_p,dimsizes(radii),numTIMES/),float)
  pvflx_plot0!1 = "rad"
  pvflx_plot0!2 = "time"

; Vertical advection of eddy tangential momentum (RHS #4)
  pvadv_plot0   = new((/num_p,dimsizes(radii),numTIMES/),float)
  pvadv_plot0!1 = "rad"
  pvadv_plot0!2 = "time"

;===================================================
; Add metadata for azimuthally averaged variables
;===================================================

  vtan_azi!0	       = "lev"
  vtan_azi!1           = "rad"
  vtan_azi&rad         = radii
  vtan_azi@description = "Azimuthally averaged tangential wind"
  vtan_azi@units       = "m s~S~-1~N~"

  vrad_azi!0	       = "lev"
  vrad_azi!1           = "rad"
  vrad_azi&rad         = radii
  vrad_azi@description = "Azimuthally averaged radial wind"
  vrad_azi@units       = "m s~S~-1~N~"

  spd_azi!0	       = "lev"
  spd_azi!1            = "rad"
  spd_azi&rad          = radii
  spd_azi@description  = "Azimuthally averaged vector wind"
  spd_azi@units        = "m s~S~-1~N~"

  vvel_azi!0	       = "lev"
  vvel_azi!1           = "rad"
  vvel_azi&rad         = radii
  vvel_azi@description = "Azimuthally averaged vertical velocity"
  vvel_azi@units       = "m s~S~-1~N~"

  vort_azi!0	       = "lev"
  vort_azi!1           = "rad"
  vort_azi&rad         = radii
  vort_azi@description = "Azimuthally averaged relative vorticity"
  vort_azi@units       = "s~S~-1~N~"

  avo_azi!0	       = "lev"
  avo_azi!1            = "rad"
  avo_azi&rad          = radii
  avo_azi@description  = "Azimuthally averaged absolute vorticity"
  avo_azi@units        = "s~S~-1~N~"

  th_azi!0	       = "lev"
  th_azi!1             = "rad"
  th_azi&rad           = radii
  th_azi@description   = "Azimuthally averaged potential temperature"
  th_azi@units         = "K"

  frcx_azi!0	       = "lev"
  frcx_azi!1           = "rad"
  frcx_azi&rad         = radii
  frcx_azi@description = "Azimuthally averaged vertical diffusion (x)"
  frcx_azi@units       = "m s~S~-1~N~"

  frcy_azi!0	       = "lev"
  frcy_azi!1           = "rad"
  frcy_azi&rad         = radii
  frcy_azi@description = "Azimuthally averaged vertical diffusion (y)"
  frcy_azi@units       = "m s~S~-1~N~"

  fric_azi!0	       = "lev"
  fric_azi!1           = "rad"
  fric_azi&rad         = radii
  fric_azi@description = "Azimuthally averaged vertical diffusion (tot)"
  fric_azi@units       = "m s~S~-1~N~"

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
    fili_pc     = diri+ens0+"_pc_"
    fili_pf     = diri+ens0+"_pf_"
    input_pb    = (/fili_pb+"05.nc",fili_pb+"06.nc"/)
    input_pc    = (/fili_pc+"05.nc",fili_pc+"06.nc"/)
    input_pf    = (/fili_pf+"05.nc",fili_pf+"06.nc"/)
    numINPUT0   = dimsizes(input_pb)        ; Number of input files
    numINPUT	= (numINPUT0 * 6) + 1	    ; Number of input files
    numTIMES    = 120                       ; Total number of times in all files
    ct          = 29                        ; Counter variable
    ct0		= 0			    ; Counter variable (independent of position)
    llbox	= toint(dist*50)	    ; Take user input to calculate line below
    dsize	= (/llbox,llbox/)	    ; No. of lat/lon grid points in domain subset

    setvalues NhlGetWorkspaceObjectId
      "wsMaximumSize" : 1000000000
    end setvalues

    date_str = new(numTIMES,string)
    time_str = new(numTIMES,string)
    out_str  = new(numTIMES,string)

  ; Create array to hold values of frictional tendency
    dtaux_dz = new((/numTIMES,num_p,dsize(0),dsize(1)/),float)
    dtauy_dz = new((/numTIMES,num_p,dsize(0),dsize(1)/),float)

;===============================
; Start loop over input files 
;===============================

    do nf = 0,  numINPUT0-1		  ; Loop over input times 
    f = addfile(input_pb(nf),"r")         ; Read in 'pb' stream [u,v,w,vort,div]
    g = addfile(input_pc(nf),"r")         ; Read in 'pc' stream [friction]
    h = addfile(input_pf(nf),"r")         ; Read in 'pf' stream [density on model levs]

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

   print_clock("Working on time: "+time_str(ct0))
   print_clock("ct = "+ct+" ; ct0 = "+ct0)
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

;    p_cp = (/clev-1,clev,clev+1/)    ; User-defined levels

    pres  = f->p         	     ; Pressure (levels)
;    pres  = f->p(p_cp)               ; Pressure (levels)
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

    dxF   = lonF(1) - lonF(0) 		       ; Grid spacing (longitude)
    dyF   = latF(1) - latF(0) 		       ; Grid spacing (latitude)

;=========================================================================
; Get density from pf stream [same size as 'strx', smaller than 'stry']
;=========================================================================

    rho1  = h->field27(it,:,{t0:t1},{n0:n1})   ; Dry air density (kg/m-3)

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
    fric_x_hgt  = new((/dimsizes(lev_arr),dsize(0),dsize(1)/),float)
    fric_y_hgt  = new((/dimsizes(lev_arr),dsize(0),dsize(1)/),float)

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
     ilev0 = clev0	; Option of a different level for vertical velocity

   do ilev = clev-1, clev+1	      		; Single level (+ level above/below)

   ;========================================
   ; Model levels ['pc' and 'pf' streams]
   ;========================================

;     hy = sprintf("%0.0f",hyb(ilev))
     hy = sprintf("%0.0f",lev_arr(ilev))
     print("Working on hybrid level " + ilev+ " ("+hy+" m AGL)")
     m0 = lev_arr(ilev)       ; Selected model levels, chosen to match pressure levels

   ; Calculate components of friction term
     ip1 = min((/m0+1,lev_arr(l1)/)) ; Account for top and bottom levels
     im1 = max((/m0-1,lev_arr(l0)/))

   ; Current level [wind stress]
     taux   = strx(m0,:,:)
     tauy   = stry(m0,:,:)

   ; Density ['pf' stream]
     rho    = rho1(m0,:,:)
     rho_a = (1/rho) ; Specific volume (inverse density)

   ; Levels above/below for vertical derivatives [wind stress]
     tau_x1 = strx(ip1,:,:)
     tau_x0 = strx(im1,:,:)
     tau_y1 = stry(ip1,:,:)
     tau_y0 = stry(im1,:,:)

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
     u_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Zonal wind
     v_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Meridional wind
     spd_int  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vector wind 
     w_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vertical velocity
     vort_int = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Relative vorticity
     avo_int  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Absolute vorticity
     th_int   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Potential temp.
     v_rad    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Radial wind
     v_tan    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Tangential wind
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
     v_int@description = "Vertical velocity"
     v_int@units       = "m s~S~-1~N~"

     copy_VarAtts(vort_plane,vort_int)
     vort_int@description = "Relative vorticity"
     vort_int@units       = "s~S~-1~N~"

     copy_VarAtts(avo_plane,avo_int)
     avo_int@description  = "Relative vorticity"
     avo_int@units        = "s~S~-1~N~"

     copy_VarAtts(th_plane,th_int)
     th_int@description   = "Potential temperature"
     th_int@units         = "K"

     copy_VarAtts(frcx_plane,frcx_int)
     frcx_int@description = "Vertical diffusion (x-component)"
     frcx_int@units       = "m s~S~-1~N~"

     copy_VarAtts(frcy_plane,frcy_int)
     frcy_int@description = "Vertical diffusion (y-component)"
     frcy_int@units       = "m s~S~-1~N~"

     copy_VarAtts(frcx_plane,fric_int)
     fric_int@description = "Vertical diffusion (total)"
     fric_int@units       = "m s~S~-1~N~"

     v_rad@description = "Radial velocity"
     v_rad@units       = "m s~S~-1~N~"

     v_tan@description = "Tangential velocity"
     v_tan@units       = "m s~S~-1~N~"

     xpos!0       = "rad"
     xpos!1       = "azi"
     xpos&rad     = radii
     xpos&azi     = thetas

     ypos!0       = "rad"
     ypos!1       = "azi"
     ypos&rad     = radii
     ypos&azi     = thetas

     u_int!0      = "rad"
     u_int!1      = "azi"
     u_int&rad    = radii
     u_int&azi    = thetas

     v_int!0      = "rad"
     v_int!1      = "azi"
     v_int&rad    = radii
     v_int&azi    = thetas

     spd_int!0   = "rad"
     spd_int!1   = "azi"
     spd_int&rad = radii
     spd_int&azi = thetas

     w_int!0      = "rad"
     w_int!1      = "azi"
     w_int&rad    = radii
     w_int&azi    = thetas

     vort_int!0   = "rad"
     vort_int!1   = "azi"
     vort_int&rad = radii
     vort_int&azi = thetas

     avo_int!0    = "rad"
     avo_int!1    = "azi"
     avo_int&rad  = radii
     avo_int&azi  = thetas

     th_int!0     = "rad"
     th_int!1     = "azi"
     th_int&rad   = radii
     th_int&azi   = thetas

     frcx_int!0   = "rad"
     frcx_int!1   = "azi"
     frcx_int&rad = radii
     frcx_int&azi = thetas

     frcy_int!0   = "rad"
     frcy_int!1   = "azi"
     frcy_int&rad = radii
     frcy_int&azi = thetas

     fric_int!0   = "rad"
     fric_int!1   = "azi"
     fric_int&rad = radii
     fric_int&azi = thetas

     v_tan!0      = "rad"
     v_tan!1      = "azi"
     v_tan&rad    = radii
     v_tan&azi    = thetas

     v_rad!0      = "rad"
     v_rad!1      = "azi"
     v_rad&rad    = radii
     v_rad&azi    = thetas

   ; Latitude/longitude (small grid)
     do irad = 0, dimsizes(radii)-1
       r = radii(irad)
       do iang = 0, dimsizes(thetas)-1
         theta = thetas_rad(iang)
           xpos(irad,iang) = centre(ct,1) + r*cos(theta)
           ypos(irad,iang) = centre(ct,0) + r*sin(theta)
       end do
     end do

     n2 = dimsizes(lon1_sub)
     t2	= dimsizes(lat1_sub)

  ; Loop over azimuth angles to fill entire circular array of data points
    do iang = 0, dimsizes(thetas_rad)-1
      u_int(:,iang)    = linint2_points(lon_sub,lat_sub,u_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      v_int(:,iang)    = linint2_points(lon_sub,lat_sub,v_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      w_int(:,iang)    = linint2_points(lon1_sub,lat1_sub,w_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      th_int(:,iang)   = linint2_points(lon1_sub,lat1_sub,th_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      vort_int(:,iang) = linint2_points(lon_sub,lat_sub,vort_plane,False,\
					xpos(:,iang),ypos(:,iang),0)
      avo_int(:,iang)  = linint2_points(lon_sub,lat_sub,avo_plane,False,\
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

  ; Calculate vertical diffusion and radial/tangential wind at each point
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
        vtan_azi(ilev,irad) = avg(v_tan(irad,:))      ; Tangential wind
	vrad_azi(ilev,irad) = avg(v_rad(irad,:))      ; Radial wind 
	spd_azi(ilev,irad)  = avg(spd_int(irad,:))    ; Vector wind (for RMW)
	vvel_azi(ilev,irad) = avg(w_int(irad,:))      ; Vertical velocity
	vort_azi(ilev,irad) = avg(vort_int(irad,:))   ; Relative vorticity
	avo_azi(ilev,irad)  = avg(avo_int(irad,:))    ; Absolute vorticity
	th_azi(ilev,irad)   = avg(th_int(irad,:))     ; Potential temperature
	fric_azi(ilev,irad) = avg(fric_int(irad,:))   ; Vertical diffusion (total)
      end do

    ;====================================
    ; Calculate perturbation variables
    ;====================================

      do irad = 0, dimsizes(radii)-1

        ptb_vtan1(ilev,irad,:,ct) = v_tan(irad,:) - vtan_azi(ilev,irad)    ; [v']
        ptb_vvel1(ilev,irad,:,ct) = w_int(irad,:) - vvel_azi(ilev,irad)    ; [w']
        ptb_vrad0(ilev,irad,:)    = v_rad(irad,:) - vrad_azi(ilev,irad)    ; [u']
        ptb_vort0(ilev,irad,:)    = vort_int(irad,:) - vort_azi(ilev,irad) ; [vrt']

	ptb_vflx0(ilev,irad,:)    = ptb_vrad0(ilev,irad,:) * \		   ; [u' * vrt']
			       	    ptb_vort0(ilev,irad,:)

      end do

    ;====================================================================================
    ; Azimuthally average the perturbation variables ( {} = brackets in Persing et al.)
    ;====================================================================================

      do irad = 0, dimsizes(radii)-1

        ptb_vflx(ilev,irad) = avg(ptb_vflx0(ilev,irad,:))         ; [ {u' * vort'} ]

      end do
 
;    ;=============================
;    ; Calculate location of RMW
;    ;=============================

;      spd_max = max(spd_azi)
;      dims_s  = dimsizes(spd_azi)
;      spd_1d  = ndtooned(spd_azi)
;      inds_s  = ind_resolve(maxind(spd_1d),dims_s)
;      rmw     = radii(inds_s(0,0))

;    ; DO THE ARRAYS BELOW NEED TO BE REFERENCED BY 'ct0' RATHER THAN 'ct'? 
;      spd_plot(ct)      = rmw			 ; Radius of maximum wind       

    ;===================================
    ; Input values into larger arrays 
    ;===================================

    ; Azimuthally-averaged variables 
      vtan_plot0(ilev,:,ct)  = vtan_azi(ilev,:)        ; Tangential wind       [v]
      vrad_plot0(ilev,:,ct)  = vrad_azi(ilev,:)        ; Radial wind           [u]
      vvel_plot0(ilev,:,ct)  = vvel_azi(ilev,:)        ; Vertical velocity     [w]
      vort_plot0(ilev,:,ct)  = vort_azi(ilev,:)        ; Rel. vorticity        [vort]
      avo_plot0(ilev,:,ct)   = avo_azi(ilev,:)         ; Abs. vorticity        [f + vort]
      th_plot0(ilev,:,ct)    = th_azi(ilev,:)          ; Potential temperature [theta]
      fric_plot0(ilev,:,ct)  = fric_azi(ilev,:)        ; Vertical diffusion    [fx+fy]

    ; Perturbation variable(s)
      pvflx_plot0(ilev,:,ct) = (/ptb_vflx(ilev,:)/)    ; Eddy vorticity flux   [fx + fy]

    ; Tidy up
      delete([/u_int,v_int,w_int,spd_int,vort_int,avo_int,th_int,v_tan,v_rad/])
      delete([/frcx_int,frcy_int,fric_int,frcx_plane,frcy_plane/])
      delete([/u_plane,v_plane,w_plane,vort_plane,avo_plane,geo_plane,th_plane/])
      delete([/lat_plane,lon_plane,lat1_plane,lon1_plane/])
      delete([/lat_sub,lat1_sub,lon_sub,lon1_sub/])

      end do			    ; End pressure level loop (do ilev = clev-1, clev+1)

    ct  = ct + 1		    ; Counter variable (time)
    ct0 = ct0 + 1		    ; Independent counter variable (time)

  ; Tidy up before next iteration
    delete([/year,month,day,hour,minute,second,utc_date/])
    delete([/lon,lat,lon1,lat1,time,times,u,v,z,t,w,vort,pres,plevs/])
    delete([/lonF,latF,lon1F,lat1F,hyb,strx,stry/])

   end do     ; End time loop (do it = 0, times-1)

   delete([/tm0,tm1/])

 end do     ; End input file loop (do nf = 0, numINPUT0-1)

;======================================================
; Calculate vertical derivatives and remaining terms
;======================================================

    ; Change units before calculations
      pres  = f->p                     ; Pressure (levels)
      plevs = dimsizes(pres)           ; Size of pressure level array
      pres@units = "hPa"

      pres = pres * 100
      pres@units = "Pa"

      do it = 0, numTIMES-1

      do ilev = 0, plevs-1

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

	  end do                   ; End loop over pressure levels (do ilev = 0, num_p-1)

       end do                    ; End loop over times (do it = 0, numTIMES-1)

     ; Return units to previous values
       pres = pres / 100
       pres@units = "hPa"

     ; Also calculate the azimuthal average of [-w' * dv'/dz]
     ;;;;; MAKE SURE THIS IS CORRECT ;;;;;

     pvadv_plot0(:,:,:) = dim_avg_n(ptb_vadv1(:,:,:,:), 2)      ; [ {w' * dv'/dz} ]

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
     ; RHS #7  ; {D_v2}           - Combined diffusive and PBL tendency (horizontal)

     ; Define new arrays for terms in momentum budget
       lhsA  = new((/num_p,dimsizes(radii),numTIMES/),float)   ; d{v}/dt
       rhs1A = new((/num_p,dimsizes(radii),numTIMES/),float)   ; -{u} * {f+vort}
       rhs2A = new((/num_p,dimsizes(radii),numTIMES/),float)   ; -{w} * (d{v}/dz)
       rhs3A = new((/num_p,dimsizes(radii),numTIMES/),float)   ; -{u' * vort'}
       rhs4A = new((/num_p,dimsizes(radii),numTIMES/),float)   ; -{w' * (dv'/dz)}
       rhs5A = new((/num_p,dimsizes(radii),numTIMES/),float)   ; PGF term
       rhs6A = new((/num_p,dimsizes(radii),numTIMES/),float)   ; {D_v1}
       rhs7A = new((/num_p,dimsizes(radii),numTIMES/),float)   ; {D_v2}
       rhsTA = new((/num_p,dimsizes(radii),numTIMES/),float)   ; Sum of rhs1-7
       resA  = new((/num_p,dimsizes(radii),numTIMES/),float)   ; Residual (LHS - RHS)

       ts   = 30
       tf   = 36

     ; Calculate mean tangential momentum tendency (LHS)
       do it = ts, tf

         cp1 = min((/it+1,tf/))
         cm1 = max((/it-1,ts/))

         lhsA(:,:,it) = (vtan_plot0(:,:,cp1) - vtan_plot0(:,:,cm1)) / 2

       end do  ; End loop over times (do it = ts, tf)

     ; Add metadata and make units consistent for all terms
       lhsA@units = "m s~S~-1~N~ h~S~-1~N~"

       rhs1A = -vrad_plot0 * avo_plot0   ; Mean radial influx of abs. vert. vorticity
       rhs1A = rhs1A * 3600
       rhs1A@units = "m s~S~-1~N~ h~S~-1~N~"

       rhs2A = -vvel_plot0 * mvtdz_plot0 ; Mean vert. adv. of mean tangential momentum
       rhs2A = rhs2A * 3600
       rhs2A@units = "m s~S~-1~N~ h~S~-1~N~"

       rhs3A = -pvflx_plot0              ; Eddy radial vorticity flux
       rhs3A = rhs3A * 3600
       rhs3A@units = "m s~S~-1~N~ h~S~-1~N~"

       rhs4A = -pvadv_plot0              ; Vertical advection of eddy tangential momentum
       rhs4A = rhs4A * 3600
       rhs4A@units = "m s~S~-1~N~ h~S~-1~N~"

       rhs6A = fric_plot0                ; Vertical diffusion
       rhs6A = rhs6A * 3600
       rhs6A@units = "m s~S~-1~N~ h~S~-1~N~"

       rhsTA = rhs1A + rhs2A + rhs3A + rhs4A + rhs6A ; Sum of RHS terms (add rhs5 & rhs7)
       rhsTA@units = "m s~S~-1~N~ h~S~-1~N~"

       resA  = lhsA - rhsTA                         ; Residual (LHS - RHS)
       resA@units = "m s~S~-1~N~ h~S~-1~N~"

       lhsA!0   = "lev"
       lhsA&lev = pres
       lhsA!1   = "rad"
       lhsA&rad = radii
       lhsA!2   = "time"

       rhs1A!0   = "lev"
       rhs1A&lev = pres
       rhs1A!1   = "rad"
       rhs1A&rad = radii
       rhs1A!2   = "time"

       rhs2A!0   = "lev"
       rhs2A&lev = pres
       rhs2A!1   = "rad"
       rhs2A&rad = radii
       rhs2A!2   = "time"

       rhs3A!0   = "lev"
       rhs3A&lev = pres
       rhs3A!1   = "rad"
       rhs3A&rad = radii
       rhs3A!2   = "time"

       rhs4A!0   = "lev"
       rhs4A&lev = pres
       rhs4A!1   = "rad"
       rhs4A&rad = radii
       rhs4A!2   = "time"

       rhs6A!0   = "lev"
       rhs6A&lev = pres
       rhs6A!1   = "rad"
       rhs6A&rad = radii
       rhs6A!2   = "time"

       rhsTA!0   = "lev"
       rhsTA&lev = pres
       rhsTA!1   = "rad"
       rhsTA&rad = radii
       rhsTA!2   = "time"

       resA!0    = "lev"
       resA&lev  = pres
       resA!1    = "rad"
       resA&rad  = radii
       resA!2    = "time"

;=======================================================
; Reorder dimensions before ploting (gsn_csm_contour)
;=======================================================

       lhs  = lhsA(lev|:,time|:,rad|:)
       rhs1 = rhs1A(lev|:,time|:,rad|:)
       rhs2 = rhs2A(lev|:,time|:,rad|:)
       rhs3 = rhs3A(lev|:,time|:,rad|:)
       rhs4 = rhs4A(lev|:,time|:,rad|:)
       rhs6 = rhs6A(lev|:,time|:,rad|:)
       rhsT = rhsTA(lev|:,time|:,rad|:)
       res  = resA(lev|:,time|:,rad|:)

;==========================================================
; Create a Hovmoller panel plot of momentum budget terms
;==========================================================

  ; Output file location and type
    output = "$sam/nepartak/images/cyl_coords/mbgt_fig6_"+ens0
    wks = gsn_open_wks(opt,output)

  ; Set up panel plot
    panel = new(3,graphic)
    gsn_define_colormap(wks,"BlueDarkRed18")

;==========================
; Options for plotting 
;==========================

  ; Momentum budget terms
    opts_mbgt                              = True
    opts_mbgt@cnFillOn                     = True
    opts_mbgt@cnLineLabelInterval          = 2.0
    opts_mbgt@cnLineLabelFontHeightF       = 0.012
    opts_mbgt@cnLineLabelBackgroundColor   = "transparent"
    opts_mbgt@cnLineLabelPlacementMode     = "constant"
    opts_mbgt@cnLinesOn                    = False ; Contour lines off
    opts_mbgt@cnInfoLabelOn                = False
    opts_mbgt@cnLevelSelectionMode         = "ExplicitLevels"
    opts_mbgt@cnLevels                     = (/-10., -5., -3., -2., -1., -0.5, -0.2, \
                                               0, 0.2, 0.5, 1., 2., 3., 5., 10./)
    opts_mbgt@cnFillColors                 = (/2,4,5,6,7,8,9,0,0,12,13,14,15,16,17,19/)
    opts_mbgt@gsnPaperOrientation          = "landscape"
    opts_mbgt@tiMainString                 = ""
    opts_mbgt@tiMainFontHeightF            = 0.0125
    opts_mbgt@gsnLeftString                = ""
    opts_mbgt@gsnRightString               = ""
    opts_mbgt@gsnMaximize		   = True
    opts_mbgt@lbLabelBarOn		   = False		   ; Labelbar on/off
    opts_mbgt@gsnDraw			   = False		   ; Do not draw plot
    opts_mbgt@gsnFrame			   = False		   ; Do not advance frame
    opts_mbgt@vpWidthF            	   = 0.30		   ; Match width in N11
    opts_mbgt@vpHeightF                    = 0.50		   ; Match height in N11

  ; Additional plotting resources
    opts_mbgt@tiYAxisString                = "Time (hour)"
    opts_mbgt@tiXAxisString		   = "Radius (km)"
    opts_mbgt@tiXAxisFontHeightF           = 0.018
    opts_mbgt@trYLog                       = False
    ys = 29
    yf = 40
    opts_mbgt@trYMinF                      = ys      ; Set y-axis min (starting time)
    opts_mbgt@trYMaxF                      = yf      ; Set y-axis max (ending time)
    opts_mbgt@trXMinF                      = 0.0     ; Set x-axis min (centre of TC)
    opts_mbgt@trXMaxF                      = 1.0     ; Set x-axis max (outer radius)
    opts_mbgt@tmXTOn                       = "False" ; Turn off top x-axis TM
    opts_mbgt@tmYROn                       = "False" ; Turn off right y-axis TM
    opts_mbgt@gsnMaximize                  = True    ; Maximise plot size
    opts_mbgt@gsnAddCyclic                 = False
    opts_mbgt@tmXBMode                     = "Explicit"                ; Choose tm mode
    opts_mbgt@tmXBValues                   = (/0.0,0.5,1.0,1.5/)       ; Tm positions
    opts_mbgt@tmXBLabels                   = (/"0","50","100","150"/)  ; Tm values
    opts_mbgt@tmXBMinorOn                  = False

    opts_mbgt@pmLabelBarOrthogonalPosF     = 0.10    ; Move lb up/down (higher/lower)
    opts_mbgt@tiXAxisOffsetYF		   = 0.10    ; Move x-axis title up/down (h/l)

    opts_mbgt@tmYLMode			   = "Manual"
    opts_mbgt@tmYLTickStartF		   = opts_mbgt@trYMinF	       ; See L730
    opts_mbgt@tmYLTickEndF		   = opts_mbgt@trYMaxF	       ; See L731
    opts_mbgt@tmYLTickSpacingF		   = 6
    opts_mbgt@tmYLMinorOn		   = False

;==============
; Panel plot 
;==============

    print("clev = "+clev)
    printVarSummary(rhs1(clev,:,:))
    printVarSummary(rhs2(clev,:,:))
    printVarSummary(rhs3(clev,:,:))

    print("RHS term #1 : "+rhs1(clev,ys:yf,:))
    print("RHS term #2 : "+rhs2(clev,ys:yf,:))
    print("RHS term #3 : "+rhs3(clev,ys:yf,:))

  ; LHS panel
    rhs_plot1 = gsn_csm_contour(wks,rhs1(clev,ys:yf,:),opts_mbgt) ; RHS #1 
    panel(0)  = rhs_plot1

  ; Middle panel
    rhs_plot2 = gsn_csm_contour(wks,rhs2(clev,ys:yf,:),opts_mbgt) ; RHS #2 
    panel(1)  = rhs_plot2

  ; RHS panel
    rhs_plot3 = gsn_csm_contour(wks,rhs3(clev,ys:yf,:),opts_mbgt) ; RHS #3
    panel(2)  = rhs_plot3

    optsP                       = True			; Panel plot resources
    optsP@gsnFrame              = False			; Do not advance the frame
    optsP@gsnPanelLabelBar      = True       	        ; Turn on panel labelbar
    optsP@lbBoxEndCapStyle      = "TriangleBothEnds"    ; Labelbar end shape
    optsP@pmLabelBarWidthF      = 0.32                  ; Labelbar width
    optsP@pmLabelBarHeightF     = 0.08                  ; Labelbar height
    optsP@lbLabelFontHeightF    = 0.0125                ; Labelbar font size
    optsP@lbLabelFont           = "Helvetica"           ; Labelbar font
    optsP@lbPerimOn             = False                 ; Perimeter on/off
    optsP@gsnDraw               = False                 ; Do not draw plot
    optsP@gsnFrame              = False                 ; Do not advance frame
    optsP@txString              = ""
    optsP@gsnPanelFigureStrings = ""
    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.98
    optsP@gsnPanelBottom        = 0.02
    optsP@amJust                = "TopLeft"
    optsP@gsnPanelFigureStringsFontHeightF = 0.0099 ; Reduce label size (default 0.01)

    gsn_panel(wks,panel,(/1,3/),optsP)             ; Draw as a single plot
    frame(wks)

end 