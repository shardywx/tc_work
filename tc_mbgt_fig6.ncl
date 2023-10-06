; Script to reproduce Figure 6 from Nguyen et al. (2011) for a given ensemble member
; Run using:

; ncl 'opt="x11"' 'ens0="em02"' clev=2 clev0=2 dist=3.0 'clr="bwr_new"' tc_mbgt_fig6.ncl

; 'opt'   = output file format ("pdf" or "x11")
; 'ens0'  = ensemble member (em00 -- em11)
; 'clev'  = pressure level for tangential wind panel (0 = 950, 1 = 900, 2 = 850, 3 = 800)
; 'clev0' = pressure level for vertical velocity panel
; 'dist'  = size of box following storm (degrees)

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
  num_p	     = 7		     ; Size of pressure level subset array (L344-345)

;=====================================================
; Define arrays for later use in radius-time plots
;=====================================================

  vtan_azi = new((/num_p,dimsizes(radii)/),float) ; Tangential wind
  vrad_azi = new((/num_p,dimsizes(radii)/),float) ; Radial wind
  spd_azi  = new((/num_p,dimsizes(radii)/),float) ; Vector velocity
  pvrt_azi = new((/num_p,dimsizes(radii)/),float) ; Potential vorticity
  vort_azi = new((/num_p,dimsizes(radii)/),float) ; Relative vorticity
  avo_azi  = new((/num_p,dimsizes(radii)/),float) ; Absolute vorticity
  th_azi   = new((/num_p,dimsizes(radii)/),float) ; Potential temperature
  vvel_azi = new((/num_p,dimsizes(radii)/),float) ; Vertical velocity

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

; Mean potential temperature
  th_plot0 	  = new((/num_p,dimsizes(radii),numTIMES/),float)
  th_plot0!0	  = "lev"
  th_plot0!1 	  = "rad"
  th_plot0&rad    = radii
  th_plot0!2 	  = "time"

; Mean vector wind 
  spd_plot0   	  = new((/num_p,numTIMES/),float)
  spd_plot0!0	  = "lev"
  spd_plot0!1 	  = "time"

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

  vort_azi!0           = "lev"
  vort_azi!1           = "rad"
  vort_azi&rad         = radii
  vort_azi@description = "Azimuthally averaged relative vorticity"
  vort_azi@units       = "s~S~-1~N~"

  th_azi!0             = "lev"
  th_azi!1             = "rad"
  th_azi&rad           = radii
  th_azi@description   = "Azimuthally averaged potential temperature"
  th_azi@units         = "K"

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

    p_cp  = (/1,2,3,4,6,8,10/) 	     ; Chosen pressure levels to access
    pres  = f->p(p_cp)         	     ; Pressure levels (950,900,850,800,700,600,500)
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
    avo@long_name   = "Vorticity (absolute)"

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
  th_plot0&lev   = pres
  spd_plot0&lev  = pres
  maxw_plot0&lev = pres

; Arrays holding azimuthally-averaged variables
  vtan_azi&lev   = pres
  vrad_azi&lev   = pres
  spd_azi&lev    = pres
  pvrt_azi&lev   = pres
  vvel_azi&lev   = pres
  avo_azi&lev    = pres
  vort_azi&lev   = pres
  th_azi&lev     = pres  

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

    ilev = clev		   ; Pressure level for tangential velocity
    ilev0 = clev0	   ; Pressure level for radial velocity

    do ilev = 0, plevs-1   ; All pressure levels in vertical level array subset

     p = pres(ilev)
     print("Working on: " + p + " hPa")

     u_plane       = u_rel(:,:)
     v_plane       = v_rel(:,:)
     vort_plane	   = vort(ilev,:,:)
     avo_plane	   = avo(ilev,:,:)
     th_plane	   = th(ilev,:,:)	
     geo_plane     = z(ilev,:,:)
     pv_plane	   = pv(ilev,:,:)
     w_plane	   = w(ilev,:,:)

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
     th_int   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Potential temperature

   ; Add variable metadata
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

     copy_VarAtts(th_plane,th_int)
     th_int@description = "Potential temperature"
     th_int@units       = "K"

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

     th_int!0     = "rad"
     th_int!1     = "azi"
     th_int&rad   = radii
     th_int&azi   = thetas

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
      pv_int(:,iang)   = linint2_points(lon_sub,lat_sub,pv_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      th_int(:,iang)   = linint2_points(lon1_sub,lat1_sub,th_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      vort_int(:,iang) = linint2_points(lon_sub,lat_sub,vort_plane,False,\
					xpos(:,iang),ypos(:,iang),0)
      avo_int(:,iang)  = linint2_points(lon_sub,lat_sub,avo_plane,False,\
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
        vtan_azi(ilev,irad) = avg(v_tan(irad,:))         ; Tangential wind
	vrad_azi(ilev,irad) = avg(v_rad(irad,:))	 ; Radial wind 
	spd_azi(ilev,irad)  = avg(spd_int(irad,:))	 ; Vector velocity (for RMW)
	pvrt_azi(ilev,irad) = avg(pv_int(irad,:))	 ; PV 
	vvel_azi(ilev,irad) = avg(w_int(irad,:))	 ; Vertical velocity
	vort_azi(ilev,irad) = avg(vort_int(irad,:))	 ; Relative vorticity
	avo_azi(ilev,irad)  = avg(avo_int(irad,:))	 ; Absolute vorticity
	th_azi(ilev,irad)   = avg(th_int(irad,:))	 ; Potential temperature
      end do

    ;====================================
    ; Calculate perturbation variables 
    ;====================================

      do irad = 0, dimsizes(radii)-1

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

      end do
 
    ;=============================
    ; Calculate location of RMW
    ;=============================

      spd_max = max(spd_azi(ilev,:))
      dims_s  = dimsizes(spd_azi(ilev,:))
      spd_1d  = ndtooned(spd_azi(ilev,:))
      inds_s  = ind_resolve(maxind(spd_1d),dims_s)
      rmw     = radii(inds_s(0,0))

      spd_plot0(ilev,ct)  = rmw			         ; Radius of maximum wind       

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
      th_plot0(ilev,:,ct)    = (/th_azi(ilev,:)/)        ; Potential temperature

    ; Perturbation variables 
      pvflx_plot0(ilev,:,ct) = (/ptb_vflx(ilev,:)/)	 ; Eddy vorticity flux (RHS #3)

    ; Tidy up
      delete([/u_int,v_int,w_int,pv_int,v_tan,v_rad/])
      delete([/u_plane,v_plane,w_plane,pv_plane,vort_plane,geo_plane/])
      delete([/lat_plane,lon_plane,lat1_plane,lon1_plane/])
      delete([/lat_sub,lat1_sub,lon_sub,lon1_sub/])

      end do			    ; End pressure level loop (do ilev = clev-1, clev+1)

    ct  = ct + 1		    ; Counter variable (time)
    ct0 = ct0 + 1		    ; Independent counter variable (time)

  ; Tidy up before next iteration
    delete([/year,month,day,hour,minute,second,utc_date/])
    delete([/lon,lat,lon1,lat1,time,times,u,v,z,t,w,vort,pres,plevs/])

   end do     ; End time loop (do it = 0, times-1)

   delete([/tm0,tm1/])

 end do     ; End input file loop (do nf = 0, numINPUT-1)

;======================================================
; Calculate vertical derivatives and remaining terms
;======================================================

    ; Change units before calculations
      p_cp  = (/1,2,3,4,6,8,10/)            ; Chosen pressure levels to access
      pres  = f->p(p_cp)                    ; Pressure levs (950,900,850,800,700,600,500)
      plevs = dimsizes(pres)                ; Size of pressure level array
      pres@units = "hPa"                    ; Pressure level array units

      pres = pres * 100
      pres@units = "Pa"

      do it = 23, 88

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

  do it = 23, 88

   tp1  = min((/it+1,88/)) 
   tm1	= max((/it-1,23/))

   vtan_tend0(:,:,it) = (vtan_plot0(:,:,tp1) - \
   			 vtan_plot0(:,:,tm1)) / 2 ; Units [m s-1 h-1] 
   vrad_tend0(:,:,it) = (vrad_plot0(:,:,tp1) - \
   			 vrad_plot0(:,:,tm1)) / 2 ; Units [m s-1 h-1]

  end do     ; End time loop (do it = 23, 88)

;======================================
; Reorder dimensions before plotting 
;======================================

 ; Arrays now ordered [PRES * TIME * RADIUS] rather than [PRES * RADIUS * TIME]

 ; Tangential and radial wind tendencies
   vtan_tend = vtan_tend0(lev|:,time|:,rad|:)
   vrad_tend = vrad_tend0(lev|:,time|:,rad|:)

 ; Tangential/radial/vertical wind 
   vtan_plot = vtan_plot0(lev|:,time|:,rad|:)
   vrad_plot = vrad_plot0(lev|:,time|:,rad|:)
   vvel_plot = vvel_plot0(lev|:,time|:,rad|:)

 ; Momentum budget source terms
   rhs1	     = rhs1A(lev|:,time|:,rad|:)
   rhs2	     = rhs2A(lev|:,time|:,rad|:)
   rhs3	     = rhs3A(lev|:,time|:,rad|:)
   rhs4	     = rhs4A(lev|:,time|:,rad|:)

;=====================================================
; Create a panel plot of Figs. 6a and 6b from N11 
;=====================================================

  ; Output file location and type
    output = "$sam/nepartak/images/cyl_coords/mbgt_fig6_"+pres(clev)+"_"+ens0
    wks = gsn_open_wks(opt,output)

  ; Set up panel plot
    panel = new(3,graphic)
    gsn_define_colormap(wks,"bwr_new")

;==========================
; Options for plotting 
;==========================

  ; Momentum budget source terms
    opts_mbgt                              = True
    opts_mbgt@cnFillOn                     = True
    opts_mbgt@cnLineLabelInterval          = 2.0
    opts_mbgt@cnLineLabelFontHeightF       = 0.012
    opts_mbgt@cnLineLabelBackgroundColor   = "transparent"
    opts_mbgt@cnLineLabelPlacementMode     = "constant"
    opts_mbgt@cnLinesOn                    = False ; Contour lines off
    opts_mbgt@cnInfoLabelOn                = False
    opts_mbgt@cnLevelSelectionMode         = "ExplicitLevels"

    if (clr.eq."BlueDarkRed18") then 
      opts_mbgt@cnLevels                     = (/-10., -5., -3., -2., -1., \
                                                 1., 2., 3., 5., 10./)
      opts_mbgt@cnFillColors                 = (/4,5,6,7,8,9,11,13,14,15,16,17/)
    else if (clr.eq."bwr_new") then
      opts_mbgt@cnLevels                     = (/-20., -15., -10., -5., -3., -2., -1., \
                                                 1., 2., 3., 5., 10., 15., 20./)
      opts_mbgt@cnFillColors                 = (/3,5,7,9,10,11,13,0,\
      					         19,21,23,25,27,28,29/)
    end if 
    end if 

;    opts_mbgt@cnLevels                     = (/-20., -15., -10., -5., -3., -2., -1., \
;                                               1., 2., 3., 5., 10., 15., 20./)
;    opts_mbgt@cnFillColors                 = (/3,5,7,9,10,11,13,0,\
;                                               19,21,23,25,27,28,29/)

    opts_mbgt@gsnPaperOrientation          = "landscape"
    opts_mbgt@tiMainString                 = ""
    opts_mbgt@tiMainFontHeightF            = 0.0125
    opts_mbgt@gsnLeftString                = ""
    opts_mbgt@gsnRightString               = ""
    opts_mbgt@gsnMaximize		   = True
    opts_mbgt@lbLabelBarOn		   = True		   ; Labelbar on/off
    opts_mbgt@lbBoxEndCapStyle		   = "TriangleBothEnds"	   ; Labelbar end shape
    opts_mbgt@pmLabelBarWidthF		   = 0.32		   ; Labelbar width
    opts_mbgt@pmLabelBarHeightF		   = 0.08		   ; Labelbar height
    opts_mbgt@lbLabelFontHeightF	   = 0.0125		   ; Labelbar font size
    opts_mbgt@lbLabelFont		   = "Helvetica"	   ; Labelbar font
    opts_mbgt@lbPerimOn			   = False		   ; Perimeter on/off
    opts_mbgt@gsnDraw			   = False		   ; Do not draw plot
    opts_mbgt@gsnFrame			   = False		   ; Do not advance frame
    opts_mbgt@vpWidthF            	   = 0.30		   ; Match width in N11
    opts_mbgt@vpHeightF                    = 0.50		   ; Match height in N11

  ; Additional plotting resources
    opts_mbgt@tiYAxisString                = "Time (hour)"
    opts_mbgt@tiXAxisString		   = "Radius (km)"
    opts_mbgt@tiXAxisFontHeightF           = 0.018
    opts_mbgt@trYLog                       = False
    opts_mbgt@trYMinF                      = 30      ; Starting point (T+...) for plot
    opts_mbgt@trYMaxF                      = 72      ; Ending point (T+...) for plot
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

    opts_mbgt@trXMinF			   = 0.0		       ; Set x-axis min
    opts_mbgt@trXMaxF			   = 1.0		       ; Set x-axis max
    opts_mbgt@tmYLMode			   = "Manual"
    opts_mbgt@tmYLTickStartF		   = opts_mbgt@trYMinF	       ; See L730
    opts_mbgt@tmYLTickEndF		   = opts_mbgt@trYMaxF	       ; See L731
    opts_mbgt@tmYLTickSpacingF		   = 6
    opts_mbgt@tmYLMinorOn		   = False

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
    opts_rmw@xyLineColors      		   = "black"
    opts_rmw@xyLineThicknesses 		   = 3.0
    opts_rmw@xyDashPatterns    		   = 0.0
    opts_rmw@gsnDraw                       = False       ; Do not draw the plot
    opts_rmw@gsnFrame                      = False       ; Do not advance the frame
    opts_rmw@tiMainString                  = ""
    opts_rmw@gsnLeftString                 = ""
    opts_rmw@gsnRightString                = ""

;==============
; Panel plot 
;==============

    ilev = clev

  ; Mean radial influx of absolute vertical vorticity (RHS #1)
    vrad0 = gsn_csm_contour(wks,rhs1(ilev,:,:),opts_mbgt)       ; MBGT - source term #1
    rmw0  = gsn_csm_xy(wks,spd_plot0(ilev,:),time_int,opts_rmw)	; RMW 
    overlay(vrad0,rmw0)						; Overlay RMW
    panel(0) = vrad0						; Plot second panel

  ; Mean vertical advection of mean tangential momentum (RHS #2)
    vrad1 = gsn_csm_contour(wks,rhs2(ilev,:,:),opts_mbgt)       ; MBGT - source term 2
    rmw1  = gsn_csm_xy(wks,spd_plot0(ilev,:),time_int,opts_rmw) ; RMW
    overlay(vrad1,rmw1)                                         ; Overlay RMW
    panel(1) = vrad1                                            ; Plot second panel

  ; Eddy radial vorticity flux (RHS #3)
    vvel0 = gsn_csm_contour(wks,rhs3(ilev,:,:),opts_mbgt)	; MBGT - source term #3
    rmw2  = gsn_csm_xy(wks,spd_plot0(ilev,:),time_int,opts_rmw)	; RMW
    overlay(vvel0,rmw2)						; Overlay RMW
    panel(2) = vvel0						; Plot third panel

    optsP                       = True
    optsP@gsnFrame              = False    ; Do not advance the frame
    optsP@gsnPanelLabelBar      = False    ; Turn on panel labelbar
    optsP@txString              = ""
    optsP@gsnPanelFigureStrings = ""
    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.90
    optsP@gsnPanelBottom        = 0.10
    optsP@amJust                = "TopLeft"
    optsP@gsnPanelFigureStringsFontHeightF = 0.0099 ; Reduce label size (default 0.01)

    gsn_panel(wks,panel,(/1,3/),optsP)             ; Draw as a single plot
    frame(wks)

end 