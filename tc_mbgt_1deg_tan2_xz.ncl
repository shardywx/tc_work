; Script to calculate the tangential momentum budget as in Montgomery et al. (2018)

; Run using:

; ncl 'dat="02T12"' 'opt="pdf"' 'ens0="em10"' dist=3.0 ts0=72 tf0=90 clr=4
; calc=1 full=1 rad0=1 ave=1 wrt=0 new0=0 tc_mbgt_tan2_xz.ncl

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'opt'   = output file format ("pdf" or "x11")
; 'ens0'  = ensemble member (em00 -- em11)
; 'dist'  = size of box following storm (degrees)
; 'ts0'   = starting time index (0-119)
; 'tf0'   = ending time index (0-119)
; 'clr'   = colour map for plots ("ncl_default" [1];  "bdr_extra" [2]; "amwg256" [3])
; 'calc'  = centre calculated on each level (0) or equal to 950 hPa centre (1)
; 'full'  = plot 13 panels inc. u,v,w (2), 10 panels w/mbgt terms only (1) or 4 panels (0)
; 'rad0'  = radius of final plots: 150 km (0) or 200 km (1)
; 'ave'   = calculate 3/6/12-h average (1) rather than analyse individual time interval (0)
; 'wrt'   = write out averages to text file (1) or leave alone (0 – default)

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"

begin

;==============================
; Prelash (reading in files)
;==============================

; Find 'pc' and 'pd' files (model height level data)
  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0
  fili_c1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pc*.nc")
  fili_c      = diri+"/"+fili_c1
  fili_d1     =	systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pd*.nc")
  fili_d      = diri+"/"+fili_d1

; Calculate number of times in each file using built-in function ('nc_times')
  numINFO_c   = nc_times(fili_c)
  numINFO_d   = nc_times(fili_d)
  numINPUT_c  = numINFO_c[0]		; Number of times in 'pc' stream
  time_c      = numINFO_c[1]		; Array of times from 'pc' stream file
  numINPUT_d  = numINFO_d[0]		; Number of times in 'pd' stream
  time_d      = numINFO_d[1]		; Array of times from 'pd' stream file 

; Calculate storm motion using built-in function ('st_rm')
; Argument 1  = ensemble simulation
; Argument 2  = path to input files
; Argument 3  = size of box following storm (degrees)

; Call external function 'st_rm'
  storm_rel = st_rm(dat,ens0,"$ar/text/",3.0)
  u_cyc     = storm_rel[0]		   ; Zonal wind 
  v_cyc     = storm_rel[1]		   ; Meridional wind
  vel_cyc   = storm_rel[2]		   ; Vector wind 
  numTIMES  = storm_rel[3]                 ; Number of times in file
  lat_arr   = storm_rel[4]                 ; Storm-following domain subset (lat)
  lon_arr   = storm_rel[5]                 ; Storm-following domain subset (lon)
  centre    = storm_rel[6]		   ; Storm centre at each hour (lat/lon)

;==========================================================================
; Create date/time string arrays for all times in file (use 'pb' stream)
;==========================================================================

  month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
                    "Jul","Aug","Sep","Oct","Nov","Dec"/)

  utc_date = cd_calendar(time_c, 0)

  year   = tointeger(utc_date(:,0))
  month  = tointeger(utc_date(:,1))
  day    = tointeger(utc_date(:,2))
  hour   = tointeger(utc_date(:,3))
  minute = tointeger(utc_date(:,4))
  second = utc_date(:,5)

; Correct for errors in the code (round up value of hour when minutes = 59)
  do ct = 0, numINPUT_c-1
   if (minute(ct).gt.30) then
     hour(ct) = hour(ct)+1
   end if
  end do
  delete(ct)

; Create arrays to hold strings for output plots
  date_str  = new(numINPUT_c,string)
  time_str  = new(numINPUT_c,string)
  time_arr  = new(numINPUT_c,string)
  title_arr = new(numINPUT_c,string)

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

  nr	     = 41
  ar	     = 2.0

  thetas     = new(121,float)
  pii        = 3.14159265
  radii      = fspan(0.0,ar,nr)      ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,121)  ; Array of azimuth angles (0, 5, 10, etc)
  thetas_rad = (pii/180)*thetas      ; Degrees to radians
  rad_size   = dimsizes(radii)-1     ; Index for later calculation
  num_l      = 59      	   	     ; Number of model levels

; Define constants
  rd	= 287.0	  ; Gas constant for dry air (J/kg/K^2)
  cp0	= 1004.0  ; Specific heat of dry air at constant pressure (J/kg/K)
  kp	= 0.286   ; For Exner function calculations (0.286)
  pref	= 1000.0  ; Reference pressure p0 (1000 hPa)

;========================================================================================
;========================================================================================
; Calculating the terms from Montgomery et al. (2018) - their Eq.(4)
; Ignore the perturbation PGF term because it is several orders	of magnitude smaller

; LHS (1) : d{v}/dt	     - Local tendency of the mean tangential wind

; RHS (1) : -{u} * {f+vort}  - Mean radial influx of absolute vertical vorticity
; RHS (2) : -{w} * (d{v}/dz) - Mean vertical advection of mean tangential momentum
; RHS (3) : -{u' * vort'}    - Eddy radial vorticity flux
; RHS (4) : -{w' * (dv'/dz)} - Vertical advection of eddy tangential momentum
; RHS (5) : {D_v}            - Combined diffusive and PBL tendency

;========================================================================================
;========================================================================================

; Array to hold all 'xpos' and 'ypos' position arrays
  xpos_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  ypos_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)

; Arrays for perturbation (eddy) variables (4D)
  ptb_vrad0 = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [u']
  ptb_vtan0 = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [v']
  ptb_vvel0 = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [w']
  ptb_vort0 = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [vort']
  ptb_dvdz0 = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [dv'/dz]

;==============================================================
; Arrays for RHS source terms containing eddy variables (4D)
;==============================================================

; Eddy radial vorticity flux: RHS (3)
; {(u' * vort')}
; ptb_vez1  = (ptb_vrad0 * ptb_vort0)
  ptb_vez1  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [R3]

; Vertical advection of eddy tangential momentum: RHS (4)
; [-w' * dv'/dz]
; ptb_vev1  = - (ptb_vvel0 * ptb_dvdz0) 
  ptb_vev1  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [R4]

; Arrays for RHS source terms after azimuthally averaging (3D)
  ptb_vez0  = new((/num_l,dimsizes(radii),numTIMES/),float)		     ; [R3]
  ptb_vev0  = new((/num_l,dimsizes(radii),numTIMES/),float)  		     ; [R4]

;==================================================
; Arrays for horizontal diffusive tendency terms 
;==================================================

; Horizontal diffusivity (K_mh)
  k_mh0	    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)

; Arrays for fractional expressions in diffusive tendency source terms below (4D)
; (1) du/dy, (2) d(v/r)/dr, (3) u/r, (4) r, (5) d/dr, (6) d/dy
  dudy0	    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (1)
  dvdr0     = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (2)
  v_r0      = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (3)
  r_0       = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (4)
  ddr0	    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (5)
  ddy0	    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; (6)

; Stress tensor (tau_ry)
; { K_mh * ( (1/r * du/dy) + d(v/r)/dr ) }
; tau_ry1   = k_mh0 * ( ( (1 / r_0) * dudy0) + ( v_r0 * ddr0) ) 
  tau_ry1   = new((/num_l,dimsizes(radii),numTIMES/),float)		     ; [3D array]

; Total horizontal diffusive tendency of tangential momentum (using 'tau_ry')
; ( (1/(r^2) * {rho}) * d( r ^ 2 * {rho} * {tau_ry} ) /dr)
; v_dh0	    = ( ( 1/( r_0 ^ 2 * mean_rho0) ) * 
; 	      (ddr0 * ( r_0 ^ 2 * mean_rho0 * tau_ry1) ) )
  v_dh0     = new((/num_l,dimsizes(radii),numTIMES/),float)		     ; [3D array]

; Total vertical diffusive tendency of tangential momentum 
; v_dz0	    = (1 /{rho0}) * ( (d/dz) * ({rho0} * {tau_yz}) )  
  v_dz0	    = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float) ; [4D array]

;========================================
; Define 3D arrays for mean quantities 
;========================================

; Mean tangential wind {v}
  mean_vtan0       = new((/num_l,dimsizes(radii),numTIMES/),float)
  mean_vtan0!0     = "lev"
  mean_vtan0!1     = "rad"
  mean_vtan0!2     = "time"
  mean_vtan0@description = "Azimuthally averaged tangential wind"
  mean_vtan0@units = "m s~S~-1~N~"

; Mean radial wind {u}
  mean_vrad0       = new((/num_l,dimsizes(radii),numTIMES/),float)
  mean_vrad0!0     = "lev"
  mean_vrad0!1     = "rad"
  mean_vrad0!2     = "time"
  mean_vrad0@description = "Azimuthally averaged radial wind"
  mean_vrad0@units = "m s~S~-1~N~"

; Mean vertical velocity {w}
  mean_vvel0       = new((/num_l,dimsizes(radii),numTIMES/),float)
  mean_vvel0!0     = "lev"
  mean_vvel0!1     = "rad"
  mean_vvel0!2     = "time"
  mean_vvel0@description = "Azimuthally averaged vertical velocity"
  mean_vvel0@units = "m s~S~-1~N~"

; Mean absolute angular momentum {aam}
  mean_aam0        = new((/num_l,dimsizes(radii),numTIMES/),float)
  mean_aam0!0      = "lev"
  mean_aam0!1      = "rad"
  mean_aam0!2      = "time"
  mean_aam0@description = "Azimuthally averaged absolute angular momentum"
  mean_aam0@units  = "m s~S~-1~N~"

; Mean horizontal diffusivity {k_mh}
  mean_kmh0        = new((/num_l,dimsizes(radii),numTIMES/),float)
  mean_kmh0!0      = "lev"
  mean_kmh0!1      = "rad"
  mean_kmh0!2      = "time"
  mean_kmh0@description = "Azimuthally averaged horizontal diffusivity"
  mean_kmh0@units  = "m~S~2~N~ s~S~-1~N~"

; Mean relative vorticity {zeta}
  mean_vort0       = new((/num_l,dimsizes(radii),numTIMES/),float)
  mean_vort0!0     = "lev"
  mean_vort0!1     = "rad"
  mean_vort0!2     = "time"
  mean_vort0@description = "Azimuthally averaged relative vorticity"
  mean_vort0@units = "s~S~-1~N~"

; Mean absolute vorticity {zeta + f0}
  mean_avo0        = new((/num_l,dimsizes(radii),numTIMES/),float)
  mean_avo0!0      = "lev"
  mean_avo0!1      = "rad"
  mean_avo0!2      = "time"
  mean_avo0@description = "Azimuthally averaged absolute vorticity"
  mean_avo0@units  = "s~S~-1~N~"

; Mean Coriolis force {f0}
  mean_cor0        = new((/num_l,dimsizes(radii),numTIMES/),float)
  mean_cor0!0      = "lev"
  mean_cor0!1      = "rad"
  mean_cor0!2      = "time"
  mean_cor0@description = "Azimuthally averaged Coriolis force"
  mean_cor0@units  = "s~S~-1~N~"

; Mean density {rho}
  mean_rho0        = new((/num_l,dimsizes(radii),numTIMES/),float)
  mean_rho0!0      = "lev"
  mean_rho0!1      = "rad"
  mean_rho0!2      = "time"
  mean_rho0@description = "Azimuthally averaged density"
  mean_rho0@units  = "kg m~S~-3~N~"

; Mean vertical diffusion (total)
  mean_difv0          = new((/num_l,dimsizes(radii),numTIMES/),float)
  mean_difv0!0        = "lev"
  mean_difv0!1        = "rad"
  mean_difv0!2        = "time"
  mean_difv0@description = "Azimuthally averaged vertical diffusion (tot)"
  mean_difv0@units    = "m s~S~-1~N~"

; Vertical derivative of mean tangential wind [ d{v} / dz ]
  mean_dvdz0       = new((/num_l,dimsizes(radii),numTIMES/),float)
  mean_dvdz0!0     = "lev"
  mean_dvdz0!1     = "rad"
  mean_dvdz0!2     = "time"
  mean_dvdz0@description = "Azimuthally averaged 'dv_dz'"
  mean_dvdz0@units = "s~S~-1~N~"

;=========================
; Start multiple loops
;=========================

    print_clock("Working on ensemble member: "+ens0)

  ; 12Z 2 Jul [00], 18Z 2 Jul [01], 00Z 3 Jul [02], 06Z 3 Jul [03], 12Z 3 Jul [04],
  ; 18Z 3 Jul [05], 00Z 4 Jul [06], 06Z 4 Jul [07], 12Z 4 Jul [08], 18Z 4 Jul [09],
  ; 00Z 5 Jul [10], 06Z 5 Jul [11], 12Z 5 Jul [12], 18Z 5 Jul [13], 00Z 6 Jul [14],
  ; 06Z 6 Jul [15], 12Z 6 Jul [16], 18Z 6 Jul [17], 00Z 7 Jul [18], 06Z 7 Jul [19],
  ; 12Z 7 Jul [20]

  ; ct = 0   [00], ct = 5  [01], ct = 11  [02], ct = 17  [03], ct = 23  [04],
  ; ct = 29  [05], ct = 35 [06], ct = 41  [07], ct = 47  [08], ct = 53  [09],
  ; ct = 59  [10], ct = 65 [11], ct = 71  [12], ct = 77  [13], ct = 83  [14],
  ; ct = 89  [15], ct = 95 [16], ct = 101 [17], ct = 107 [18], ct = 113 [19],
  ; ct = 119 [20]

    numTIMES    = 120                       ; Total number of times in all files [36]
    ct          = ts0                       ; Counter variable
    llbox	= toint(dist*50)	    ; Size of lat/lon arrays (see line below)
    dsize	= (/llbox,llbox/)           ; Size of lat/lon arrays (depends on 'dist')
    d0		= dsize(0)
    d1		= dsize(1)

    setvalues NhlGetWorkspaceObjectId
      "wsMaximumSize" : 1000000000
    end setvalues

;===============================
; Start loop over input files 
;===============================

    c  = addfile(fili_c,"r")         	 ; Read in 'pc' stream [u,v,vort]
    d  = addfile(fili_d,"r")         	 ; Read in 'pd' stream [w,p,t,theta,kmh,tau]

;==================================
; Get the variables we will need
;==================================

   do it = ts0, tf0		; Loop over times in file

    print("Working on time: "+time_str(it))
    itt = it + 1
    title_arr(it) = "Valid at "+time_str(it)+" (T+"+itt+")"

;=====================================================
; Read in variables from 'pc' stream (model levels)
;=====================================================

    t0    = lat_arr(it,0)
    t1    = lat_arr(it,1)
    n0    = lon_arr(it,0)
    n1    = lon_arr(it,1)
    print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

  ; Read in basic variables from 'pc' stream
    lon1C  = c->longitude_1({n0:n1})              ; longitude ['d1' grid points]
    lat1C  = c->latitude_1({t0:t1})               ; latitude  ['d0' grid points]
    lonC   = c->longitude({n0:n1})                ; longitude ['d1' grid points]
    latC   = c->latitude({t0:t1})                 ; latitude  ['d0' grid points]
    hybC   = c->hybrid_ht(:) 		          ; 63 model ('theta') levels

    u_varname = "x-wind"
    v_varname = "y-wind"

  ; Read in horizontal wind components 
  ; Both variables below -- [hybC | 63] * [latC | 150] * [lonC | 150]
    u0  = c->$u_varname$(it,:,{t0:t1},{n0:n1}) ; Zonal wind	
    v0  = c->$v_varname$(it,:,{t0:t1},{n0:n1}) ; Meridional wind 

  ; Calculate relative vorticity using centered finite differences 
  ; Option '2' --> boundary points estimated using one-sided difference scheme
    vort = uv2vr_cfd(u0, v0, latC, lonC, 2)
    copy_VarCoords(u0, vort)
    vort@units = "s~S~-1~N~"
    vort@name = "Relative vorticity on model levels"

  ; Coriolis parameter and absolute vorticity
    f1   = coriolis_param(latC)
    f0   = conform_dims(dimsizes(vort),f1,1)
    copy_VarMeta(vort,f0)
    copy_VarCoords(vort,f0)
    f0@description = "Coriolis parameter"
    f0@name        = "Coriolis parameter"

    avo  = vort + f0
    copy_VarMeta(vort,avo)
    copy_VarAtts(vort,avo)
    copy_VarCoords(vort,avo)
    avo@description = "Absolute vorticity"
    avo@name        = "Absolute vorticity on model levels"
    avo@long_name   = "Absolute vorticity"

  ; Calculate storm-relative horizontal wind
    u = u0 - u_cyc(it+1)
    v = v0 - v_cyc(it+1)

  ; Add metadata from original horizontal wind arrays
    copy_VarCoords(u0,u)
    copy_VarAtts(u0,u)

    copy_VarCoords(v0,v)
    copy_VarAtts(v0,v)

;==================================================================
; Read in variables from 'pd' stream (model levels) 
; B/c of differences in file size, access variables using 'it+2'
;==================================================================

  ; Read in basic variables from 'pd' stream 
    lon1D  = d->longitude_1({n0:n1})                ; longitude [xxx grid points]
    lat1D  = d->latitude_1({t0:t1})                 ; latitude  [xxx grid points]
    lonD   = d->longitude({n0:n1})                  ; longitude [xxx grid points]
    latD   = d->latitude({t0:t1})                   ; latitude  [xxx grid points]
    hybD   = d->hybrid_ht(:)			    ; 63 model ('rho') levels

  ; Read in diffusivity and wind stress
  ; 'str_xz' -- [hybD | 63] * [latD | 150]  * [lonD | 151]
  ; 'str_yz' -- [hybD | 63] * [lat1D | 150] * [lon1D | 151]
  ; 'kmh'    -- [hybD | 63] * [latD | 150]  * [lon1D | 151]
    str_xz = d->taux(it-2,:,{t0:t1},{n0:n1})        ; Wind stress, x-comp (N m/2)
    str_yz = d->tauy(it-2,:,{t0:t1},{n0:n1})        ; Wind stress, y-comp (N m/2)
    kmh	   = d->unspecified(it-2,:,{t0:t1},{n0:n1}) ; Horizontal diffusivity (m^2 s-1)

  ; Read in temperature, pressure, vertical velocity and geopotential height
  ; All variables below -- [hybD | 63] * [latD | 150] * [lon1D | 151]
    t      = d->temp(it-2,:,{t0:t1},{n0:n1})        ; Temperature (K)
    p      = d->p(it-2,:,{t0:t1},{n0:n1})           ; Pressure (Pa)
    w      = d->dz_dt(it-2,:,{t0:t1},{n0:n1})	    ; Vertical velocity (m s-1)
    z      = d->ht(it-2,:,{t0:t1},{n0:n1})    	    ; Geopotential height (m)

  ; Calculate dry air density from pressure and temperature
    rho	   = p / (rd * t)

  ; Add metadata
    copy_VarCoords(p, rho)
    rho@description = "Dry air density"
    rho@units = "kg m~S~-3~N~"

  ; Also calculate specific volume (for	use in diffusive tendency source terms)
    rho_a  = 1 / rho
    copy_VarCoords(rho, rho_a)
    rho_a@description = "Specific volume"
    rho_a@units = "m~S~3~N~ kg~S~-1~N~"

  ; Convert units of pressure and geopotential height before continuing
    p  = p / 100
    z  = z / 10

  ; Add metadata
    p@units = "hPa"
    z@units = "dam"

    dx      = lonD(1) - lonD(0) ; Grid spacing (longitude)
    dy      = latD(1) - latD(0) ; Grid spacing (latitude)

;==========================
; Loop over model levels 
;==========================

    do ilev = 0, num_l-1 ; Loop over all pressure levels

     hy = sprintf("%0.0f",hybD(ilev))
     print("Working on model level " + ilev+ " ("+hy+" m AGL)")

   ; Calculate components of vertical diffusion term
   ; accounting for top and bottom levels
     ip1 = min((/ilev+1,num_l-1/))
     im1 = max((/ilev-1,0/))

   ; Wind stress components on current model level
     tau_xz = str_xz(ilev,0:d0-1,0:d1-1)
     tau_yz = str_yz(ilev,0:d0-1,0:d1-1)

   ; Wind stress components on model levels above/below current level
     tau_x1 = str_xz(ip1,0:d0-1,0:d1-1)
     tau_x0 = str_xz(im1,0:d0-1,0:d1-1)
     tau_y1 = str_yz(ip1,0:d0-1,0:d1-1)
     tau_y0 = str_yz(im1,0:d0-1,0:d1-1)
     dz	    = hybD(ip1) - hybD(im1)

   ; Vertical derivative of x-component of wind stress
     dtaux_dz = (tau_x1 - tau_x0) / dz

   ; Vertical derivative of y-component of wind stress
     dtauy_dz = (tau_y1 - tau_y0) / dz

   ; X-component of vertical diffusive tendency
     difv_xplane             = rho_a(ilev,:,:) * dtaux_dz
     copy_VarCoords(tau_xz,difv_xplane)
     copy_VarMeta(tau_xz,difv_xplane)
     difv_xplane@units       = "m s~S~-2~N~"
     difv_xplane@title       = "vertical diffusive tendency (x-comp)"
     difv_xplane@long_name   = "vertical diffusive tendency (x-comp)"
     difv_xplane@name        = "difv_x"

   ; Y-component of vertical diffusive tendency
     difv_yplane             = rho_a(ilev,:,:) * dtauy_dz
     copy_VarCoords(tau_yz,difv_yplane)
     copy_VarMeta(tau_yz,difv_yplane)
     difv_yplane@units       = "m s~S~-2~N~"
     difv_yplane@title       = "vertical diffusive tendency (y-comp)"
     difv_yplane@long_name   = "vertical diffusive tendency (y-comp)"
     difv_yplane@name        = "difv_y"

   ; Read in other variables on model levels
     u_plane 	            = u(ilev,:,:)	   ; Zonal wind 
     v_plane		    = v(ilev,:,:)	   ; Meridional wind 
     w_plane		    = w(ilev,:,:)	   ; Vertical velocity
     f_plane		    = f0(ilev,:,:)	   ; Coriolis parameter
     rho_plane              = rho(ilev,:,:)	   ; Density
     prs_plane		    = p(ilev,:,:)	   ; Pressure 
     vort_plane	   	    = vort(ilev,:,:)	   ; Relative vorticity
     avo_plane		    = avo(ilev,:,:)	   ; Absolute vorticity
     kmh_plane		    = kmh(ilev,:,:)	   ; Horizontal diffusivity (K_mh)

     vort_max      	    = max(vort_plane)
     vort_smth     	    = smth9_Wrap(vort_plane, 0.5, 0.5, True)
     vort_max_smth 	    = max(vort_smth)

     dims          	    = dimsizes(vort_plane)
     vort1d        	    = ndtooned(vort_plane)
     inds          	    = ind_resolve(maxind(vort1d),dims)
     vort1d_sm     	    = ndtooned(vort_smth)
     inds_sm       	    = ind_resolve(maxind(vort1d_sm),dims)

   ;==================================================================
   ; Either calculate the centre, or read in pre-calculated values
   ;==================================================================

   if (calc .eq. 0) then

   ; Find latitude and longitude of storm at each time and on each model level
     lat_max 	            = t0 + (dy * inds(0,0))
     lon_max 	   	    = n0 + (dx * inds(0,1))

   ; Difference between 'lat_max'/'lon_max' and values read in using NCL function
     dlat 	   	    = abs(lat_max - centre(it,0))
     dlon	   	    = abs(lon_max - centre(it,1))

   ; If our calculations of latitude and longitude deviate too much from the values
   ; we read in using the NCL function ('st_rm'), then calculate a new centre
   ; using the minimum pressure. Constrain the new calculation using the values
   ; from 'st_rm', to prevent spurious centres being identified and analysed. 

     if (dlat .gt. 0.75 .or. dlon .gt. 0.75) then 

      lt1 = centre(it,0) - 0.5
      lt2 = centre(it,0) + 0.5
      ln1 = centre(it,1) - 0.5
      ln2 = centre(it,1) + 0.5

      prs_plane0 = prs_plane({lt1:lt2},{ln1:ln2})       ; Smaller grid (1.0 x 1.0 deg)
      prs_min0   = min(prs_plane0)                      ; Find minimum on smaller grid
      dims_p0    = dimsizes(prs_plane0)                 ; Size of smaller grid
      prs_1d0    = ndtooned(prs_plane0)                 ; Create 1-D array
      inds_p0    = ind_resolve(minind(prs_1d0),dims_p0) ; Find index of minimum

      lat_max 	 = lt1 + (dy * inds_p0(0,0))
      lon_max 	 = ln1 + (dx * inds_p0(0,1))
      print("New centre: "+lat_max+" degrees N, "+lon_max+" degrees E")

     end if 

   else if (calc .eq. 1) then

      lat_max    = centre(it,0)
      lon_max    = centre(it,1)

   end if
   end if

   ;======================================================================
   ; Create lat/lon arrays for later switch to cylindrical coordinates
   ;======================================================================

   ;===============
   ; 'pc' stream
   ;===============
  
   ; Indices of lat/lon grid for 'u' and 'v' in 'pc' stream
     latC_plane  = ind(t0.le.latC.and.latC.le.t1)
     lonC_plane  = ind(n0.le.lonC.and.lonC.le.n1)

   ; Array sizes
     latC_size   = dimsizes(latC_plane)
     lonC_size   = dimsizes(lonC_plane)
     
   ; Array of latitude/longitude points	corresponding to the indices above
     latC_sub    = latC(latC_plane(0):latC_plane(latC_size-1))
     lonC_sub    = lonC(lonC_plane(0):lonC_plane(lonC_size-1))

   ;===============
   ; 'pd' stream
   ;===============

   ; Indices of lat/lon grid for 'str_xz' in 'pd' stream
     latD_plane  = ind(t0.le.latD.and.latD.le.t1)
     lonD_plane  = ind(n0.le.lonD.and.lonD.le.n1)

   ; Array sizes
     latD_size   = dimsizes(latD_plane)
     lonD_size   = dimsizes(lonD_plane)

   ; Array of latitude/longitude points corresponding to the indices above
     latD_sub    = latD(latD_plane(0):latD_plane(latD_size-1))
     lonD_sub    = lonD(lonD_plane(0):lonD_plane(lonD_size-1))

   ; Indices of lat/lon grid for 'str_yz' in 'pd' stream
     lat1D_plane = ind(t0.le.lat1D.and.lat1D.le.t1)
     lon1D_plane = ind(n0.le.lon1D.and.lon1D.le.n1)

   ; Array sizes
     lat1D_size  = dimsizes(lat1D_plane)
     lon1D_size  = dimsizes(lon1D_plane)

   ; Array of latitude/longitude points corresponding to the indices above
     lat1D_sub   = lat1D(lat1D_plane(0):lat1D_plane(lat1D_size-1))
     lon1D_sub   = lon1D(lon1D_plane(0):lon1D_plane(lon1D_size-1))

     delete([/vort1d,inds,vort1d_sm,vort_smth,inds_sm/])

;===========================================
; Switch to cylindrical coordinate system
;===========================================

   ; Create new arrays describing cylindrical coordinate system (2D for now)
     xpos     = new((/dimsizes(radii),dimsizes(thetas)/),float)
     ypos     = new((/dimsizes(radii),dimsizes(thetas)/),float)

   ; New arrays containing regular variables we want to plot
     u_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Zonal velocity
     v_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Meridional velocity
     w_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vertical velocity
     f_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Coriolis force 
     kmh_int  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Horizontal diffusivity
     vort_int = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Relative vorticity
     avo_int  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Absolute vorticity
     rho_int  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Density
     v_rad    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Radial velocity
     v_tan    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Tangential velocity

   ; Additional arrays for horizontal and vertical diffusion
     difv_x   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vert. diff (x)
     difv_y   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vert. diff (y)
     difv_tot = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vert. diff (total)

     difh_x   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Horiz. diff (x)
     difh_y   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Horiz. diff (y)

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

     copy_VarAtts(vort_plane,vort_int)
     vort_int@description = "Relative vorticity"
     vort_int@units    = "s~S~-1~N~"

     copy_VarAtts(avo_plane,avo_int)
     avo_int@description = "Relative vorticity"
     avo_int@units    = "s~S~-1~N~"

     copy_VarAtts(f_plane,f_int)
     f_int@description = "Coriolis parameter"
     f_int@units       = "s~S~-1~N~"

     copy_VarAtts(kmh_plane,kmh_int)
     kmh_int@description = "Horizontal diffusivity"
     kmh_int@units       = "m~S~2~N~ s~S~-1~N~"

;=======================
; Vertical diffusion
;=======================

     copy_VarAtts(difv_xplane,difv_x)
     difv_x@description   = "Vertical diffusion (x-component)"
     difv_x@units         = "m s~S~-2~N~"

     copy_VarAtts(difv_yplane,difv_y)
     difv_y@description   = "Vertical diffusion (y-component)"
     difv_y@units         = "m s~S~-2~N~"

     copy_VarAtts(difv_yplane,difv_tot)
     difv_tot@description = "Vertical diffusion (total)"
     difv_tot@units       = "m s~S~-2~N~"

;=========================
; Horizontal diffusion
;=========================

     copy_VarAtts(difv_xplane,difh_x)
     difh_x@description = "Horizontal diffusion (x-component)"
     difh_x@units       = "m s~S~-2~N~"

     copy_VarAtts(difv_yplane,difh_y)
     difh_y@description = "Horizontal diffusion (y-component)"
     difh_y@units       = "m s~S~-2~N~"

     v_rad@description    = "Radial velocity"
     v_rad@units          = "m s~S~-1~N~"

     v_tan@description    = "Tangential velocity"
     v_tan@units          = "m s~S~-1~N~"

   ; Create cylindrincal coordinate array from original latitude/longitude grid
     do irad = 0, dimsizes(radii)-1
       r = radii(irad)
       do iang = 0, dimsizes(thetas)-1
          theta = thetas_rad(iang)   

           xpos(irad,iang) = lon_max + r*cos(theta)
 	   ypos(irad,iang) = lat_max + r*sin(theta)

       end do
     end do

   ; Fill 'xpos_all' and 'ypos_all' with latitude & longitude values 
     xpos_all(ilev,:,:,it) = xpos(:,:)
     ypos_all(ilev,:,:,it) = ypos(:,:)

   ; Fill new arrays with interpolated values on the cylindrical grid 
   ; Which latitude/longitude grids are each of these variables defined on? 
   ; u_plane     -- latC,  lonC
   ; v_plane     -- latC,  lonC
   ; vort_plane	 -- latC,  lonC
   ; avo_plane	 -- latC,  lonC
   ; f_plane     -- latC,  lonC

   ; w_plane     -- latD,  lon1D
   ; prs_plane   -- latD,  lon1D
   ; rho_plane	 -- latD,  lon1D
   ; kmh_plane	 -- latD,  lon1D

   ; difv_xplane -- latD,  lonD
   ; difv_yplane -- lat1D, lon1D 

     do iang = 0, dimsizes(thetas_rad)-1
      u_int(:,iang)    = linint2_points(lonC_sub,latC_sub,u_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      v_int(:,iang)    = linint2_points(lonC_sub,latC_sub,v_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      vort_int(:,iang) = linint2_points(lonC_sub,latC_sub,vort_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      avo_int(:,iang)  = linint2_points(lonC_sub,latC_sub,avo_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      f_int(:,iang)    = linint2_points(lonC_sub,latC_sub,f_plane,False,\
				          xpos(:,iang),ypos(:,iang),0)
      difv_x(:,iang)   = linint2_points(lonD_sub,latD_sub,difv_xplane,False,\
					  xpos(:,iang),ypos(:,iang),0)
      difv_y(:,iang)   = linint2_points(lon1D_sub,lat1D_sub,difv_yplane,False,\
					  xpos(:,iang),ypos(:,iang),0)
      w_int(:,iang)    = linint2_points(lon1D_sub,latD_sub,w_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      rho_int(:,iang)  = linint2_points(lon1D_sub,latD_sub,rho_plane,False,\
					  xpos(:,iang),ypos(:,iang),0)
      kmh_int(:,iang)  = linint2_points(lon1D_sub,latD_sub,kmh_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
     end do

   ; Calculate radial/tangential windspeed and friction at each point (using unit vctrs)
      do iang = 0, dimsizes(thetas)-1
       theta = thetas_rad(iang)

      ; Tangential and radial wind 
        v_rad(:,iang) = u_int(:,iang)*cos(theta) + v_int(:,iang)*sin(theta)
        v_tan(:,iang) = -u_int(:,iang)*sin(theta) + v_int(:,iang)*cos(theta)

      ; Total vertical diffusive tendency of tangential momentum
        difv_tot(:,iang) = difv_y(:,iang)*cos(theta) - difv_x(:,iang)*sin(theta)

      end do

   ; Fill 'u_dz0' array with values from 'difv_tot' (vertical diffusive tendency)
     v_dz0(ilev,:,:,ct)   = difv_tot(:,:)    

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

     vort_int!0   = "rad"
     vort_int!1   = "azi"
     vort_int&rad = radii
     vort_int&azi = thetas

     avo_int!0    = "rad"
     avo_int!1    = "azi"
     avo_int&rad  = radii
     avo_int&azi  = thetas

     f_int!0      = "rad"
     f_int!1      = "azi"
     f_int&rad    = radii
     f_int&azi    = thetas

     w_int!0      = "rad"
     w_int!1      = "azi"
     w_int&rad    = radii
     w_int&azi    = thetas

     rho_int!0    = "rad"
     rho_int!1    = "azi"
     rho_int&rad  = radii
     rho_int&azi  = thetas

     kmh_int!0    = "rad"
     kmh_int!1    = "azi"
     kmh_int&rad  = radii
     kmh_int&azi  = thetas

     v_tan!0      = "rad"
     v_tan!1      = "azi"
     v_tan&rad    = radii
     v_tan&azi    = thetas

     v_rad!0      = "rad"
     v_rad!1      = "azi"
     v_rad&rad    = radii
     v_rad&azi    = thetas

     difv_x!0     = "rad"
     difv_x!1     = "azi"
     difv_x&rad   = radii
     difv_x&azi   = thetas

     difv_y!0     = "rad"
     difv_y!1     = "azi"
     difv_y&rad   = radii
     difv_y&azi   = thetas

     difv_tot!0   = "rad"
     difv_tot!1   = "azi"
     difv_tot&rad = radii
     difv_tot&azi = thetas

     difh_x!0     = "rad"
     difh_x!1     = "azi"
     difh_x&rad   = radii
     difh_x&azi   = thetas

     difh_y!0     = "rad"
     difh_y!1     = "azi"
     difh_y&rad   = radii
     difh_y&azi   = thetas

    ;===================================================
    ; Calculate mean (azimuthally-averaged) variables 
    ;===================================================

      do irad = 0, dimsizes(radii)-1

        mean_vtan0(ilev,irad,ct)    = avg(v_tan(irad,:))      ; [ {v} ]
	mean_vrad0(ilev,irad,ct)    = avg(v_rad(irad,:))      ; [ {u} ]
	mean_vvel0(ilev,irad,ct)    = avg(w_int(irad,:))      ; [ {w} ]
	mean_vort0(ilev,irad,ct)    = avg(vort_int(irad,:))   ; [ {zeta} ]
	mean_avo0(ilev,irad,ct)	    = avg(avo_int(irad,:))    ; [ {zeta + f} ]
	mean_difv0(ilev,irad,ct)    = avg(difv_tot(irad,:))   ; [ {diff_xz + diff_yz} ]
	mean_rho0(ilev,irad,ct)     = avg(rho_int(irad,:))    ; [ {rho} ]
	mean_kmh0(ilev,irad,ct)	    = avg(kmh_int(irad,:))    ; [ {k_mh} ]
	mean_cor0(ilev,irad,ct)     = avg(f_int(irad,:))      ; [ {f} ]

      end do

      fcor = mean_cor0 * mean_vtan0

    ;====================================
    ; Calculate perturbation variables 
    ;====================================

      do irad = 0, dimsizes(radii)-1

        ptb_vtan0(ilev,irad,:,ct) = v_tan(irad,:) - mean_vtan0(ilev,irad,ct)    ; [v']
	ptb_vvel0(ilev,irad,:,ct) = w_int(irad,:) - mean_vvel0(ilev,irad,ct)    ; [w']
	ptb_vrad0(ilev,irad,:,ct) = v_rad(irad,:) - mean_vrad0(ilev,irad,ct)    ; [u']
	ptb_vort0(ilev,irad,:,ct) = vort_int(irad,:) - mean_vort0(ilev,irad,ct) ; [zeta']

	k_mh0(ilev,irad,:,ct)     = kmh_int(irad,:)				; [kmh]

      end do

    ;=========================================================
    ; Calculate term involving azimuthal derivative [du/dy]
    ;=========================================================

       do irad = 1, dimsizes(radii)-1

        do th0 = 0, dimsizes(thetas)-1
       
         tp1 = min((/th0+1,dimsizes(thetas)-1/))
       	 tm1 = max((/th0-1,0/))

       ; Calculate the azimuthal distance ('dy') between all points in array
       ; Arguments 1 and 2 represent the lat/lon of the first point, respectively
       ; Arguments 3 and 4 represent the lat/lon of the	second point, respectively
       ; Argument 5 ('2') tells the script to avoid extra interpolated points
       ; Argument 6 ('3') returns the distance in metres

         dl  = gc_latlon(ypos(irad,tm1),xpos(irad,tm1),\
                         ypos(irad,tp1),xpos(irad,tp1),2,3)

;         print("lat0 = "+ypos(irad,tm1)+ " ; lon0 = "+xpos(irad,tm1))
;	 print("lat1 = "+ypos(irad,tp1)+ " ; lon1 = "+xpos(irad,tp1))
;	 print("Distance between points  = "+dl+ " (m)")

       ; Calculate radial wind increment
         du   = v_rad(irad,tp1) - v_rad(irad,tm1)

       ; Calculate [du/dy -- for horizontal diffusive tendency]
       	 dudy0(ilev,irad,th0,ct)     = du / dl 

        end do			      ; End loop over azimuth angles
       
       end do			      ; End loop over radii 

    ;===========================================================
    ; Calculate term involving radial derivatives [d(v/r)/dr]
    ;===========================================================

       do irad = 1, dimsizes(radii)-1

        rp1 = min((/irad+1,dimsizes(radii)-1/))
        rm1 = max((/irad-1,0/))

        do th0 = 0, dimsizes(thetas)-1

         tp1 = min((/th0+1,dimsizes(thetas)-1/))
         tm1 = max((/th0-1,0/))

       ; Calculate the radial distance ('dr') between all points in array
       ; Arguments 1 and 2 represent the lat/lon of the first point, respectively 
       ; Arguments 3 and 4 represent the lat/lon of the second point, respectively 
       ; Argument 5 ('2') tells the script to avoid extra interpolated points
       ; Argument 6 ('3') returns the distance in metres

       ; Why were the latitude and longitude switched in the azimuthal derivative loop?
       ; Make sure the values I'm outputting here are correct, before continuing. 

       ; Radial distance between each point and its nearest neighbour
         dr  = gc_latlon(ypos(rm1,th0),xpos(rm1,th0),\
                         ypos(irad,th0),xpos(irad,th0),2,3)

       ; Radial distance between each point and the origin
         r1  = gc_latlon(ypos(0,th0),xpos(0,th0),\
                         ypos(irad,th0),xpos(irad,th0),2,3)

       ; Calculate components of horizontal diffusive tendency source term

       ; Calculate [v/r]
       	 v_r0(ilev,irad,th0,ct)	      = v_tan(irad,th0) / r1 

       ; Calculate [d(v/r)/dr]
         dvdr0(ilev,irad,th0,ct)      = v_r0(ilev,irad,th0,ct) / dr

       ; Calculate [r]
	 r_0(ilev,irad,th0,ct)        = r1

       ; Calculate [d/dr]
       	 ddr0(ilev,irad,th0,ct)	      = 1 / dr

        end do                        ; End loop over azimuth angles

       end do                         ; End loop over radii

    ;=============================================
    ; Calculate absolute angular momentum (AAM)
    ;=============================================

       do irad = 1, dimsizes(radii)-1
       
        rp1 = min((/irad+1,dimsizes(radii)-1/))
	rm1 = max((/irad-1,0/))

	th0 = (dimsizes(thetas) / 2)

       ; Radial distance between each point and its nearest neighbour
         dr  = gc_latlon(ypos(rm1,th0),xpos(rm1,th0),\
                         ypos(irad,th0),xpos(irad,th0),2,3)

       ; Radial distance between each point and the origin
         r1  = gc_latlon(ypos(0,th0),xpos(0,th0),\
                         ypos(irad,th0),xpos(irad,th0),2,3)

       ; Calculate increments
         vm   = mean_vtan0(ilev,irad,ct)
         fm   = mean_cor0(ilev,irad,ct) * (r1 ^ 2)

       ; Calculate absolute angular momentum (AAM)
         aam0 = vm * r1
         aam1 = 0.5 * fm

       ; Calculate absolute angular momentum [ rv + 0.5 * (fr ^ 2) ]
         mean_aam0(ilev,irad,ct)  = aam0 + aam1

       end do                         ; End loop over radii

    ; Tidy up
      delete([/u_int,v_int,w_int,v_tan,v_rad/])
      delete([/rho_int,difv_x,difv_y,f_int,kmh_int,vort_int,avo_int/])
      delete([/u_plane,v_plane,w_plane,f_plane,vort_plane,avo_plane/])
      delete([/difv_xplane,difv_yplane,kmh_plane,prs_plane,xpos,ypos/])
      delete([/latC_plane,lonC_plane,latD_plane,lonD_plane,lon1D_plane,lat1D_plane/])
      delete([/latC_sub,lonC_sub,latD_sub,lonD_sub,lat1D_sub,lon1D_sub/])

     end do			    ; End model level loop (do ilev = 0, num_l-1)

    ct  = ct + 1		    ; Counter variable (time)

  ; Tidy up before next iteration (all variables from 'pc' and 'pd' streams)
    delete([/lonC,latC,lon1C,lat1C,u,v,vort,avo,f0,f1/])
    delete([/lonD,latD,lon1D,lat1D,str_xz,str_yz,kmh,t,p,w,z,rho,rho_a/])

   end do     ; End time loop (do it = 0, times-1)

;================================================
; Calculate vertical derivatives of mean terms 
;================================================

    ; Redefine model level array for vertical derivative calculations
      hybC   = c->hybrid_ht(:)         ; 63 model ('theta') levels
      mlevs  = dimsizes(hybC)	       ; Size of model level array
      hybC@units = "m"		       ; Model level units 

      do it = ts0, tf0

	do ilev = 0, num_l-1		; Loop over all model levels
	
        mp1 = min((/ilev+1,num_l-1/))
        mm1 = max((/ilev-1,0/))
 
          do irad = 0, dimsizes(radii)-1

          ;=========================================================
	  ; Calculate vertical derivative of mean tangential wind
	  ;=========================================================

          ; [ d{v}/dz ]
	    mean_dvdz0(ilev,irad,it)    = (mean_vtan0(mp1,irad,it) - \
                                           mean_vtan0(mm1,irad,it)) / \
                                          (hybC(mp1) - hybC(mm1))

          end do		 ; End loop over radial circles (do irad = 0, ...)

        end do 		         ; End loop over model levels (do ilev = 0, num_l-1)

       end do			 ; End loop over times (do it = 0, numTIMES-1)

;================================================
; Calculate vertical derivatives of eddy terms 
;================================================

      do it = ts0, tf0

        do ilev = 0, num_l-1

        mp1 = min((/ilev+1,num_l-1/))
        mm1 = max((/ilev-1,0/))

          do irad = 0, dimsizes(radii)-1

           do th0 = 0, dimsizes(thetas)-1

           ;=========================================================
           ; Calculate vertical derivative of eddy tangential wind
           ;=========================================================

           ; [ dv'/dz ]
             ptb_dvdz0(ilev,irad,th0,it) = (ptb_vtan0(mp1,irad,th0,it) - \
                                            ptb_vtan0(mm1,irad,th0,it)) / \
                                           (hybC(mp1) - hybC(mm1))

           end do	         ; End loop over azimuth angles (do th0 = 0, ...)

          end do                 ; End loop over radial circles (do irad = 0, ...)

        end do                   ; End loop over model levels (do ilev = 0, num_l-1)

       end do                    ; End loop over times (do it = ts0, tf0)

;==========================================================
; Finally, calculate source terms and azimuthal averages 
;==========================================================

     ;=============================================================
     ; Mean radial influx of absolute vertical vorticity: RHS (1)
     ;=============================================================

       rhs1_st = -(mean_vrad0 * mean_avo0)

     ;===============================================================
     ; Mean vertical advection of mean tangential momentum: RHS (2)
     ;===============================================================

       rhs2_st = -(mean_vvel0 * mean_dvdz0)
     
     ;=======================================
     ; Eddy radial vorticity flux: RHS (3) 
     ;=======================================

       ptb_vez1 = -(ptb_vrad0 * ptb_vort0)

     ; Calculate the azimuthal average of 'ptb_vez1' and input into 'rhs3_st'
       rhs3_st	= dim_avg_n(ptb_vez1(:,:,:,:), 2)

     ;===========================================================
     ; Vertical advection of eddy tangential momentum: RHS (4)
     ;===========================================================

       ptb_vev1 = -(ptb_vvel0 * ptb_dvdz0)

     ; Calculate the azimuthal average of 'ptb_vev1' and input into 'rhs4_st'
       rhs4_st  = dim_avg_n(ptb_vev1(:,:,:,:), 2)

     ;===================================================================================
     ; Combined mean tangential/vert diffusive tendency of tangential momentum: RHS (5)
     ;===================================================================================

     ; Stress tensor
       tau_ry1	 = dim_avg_n( ( k_mh0 * ( ( (1 / r_0) * dudy0) + dvdr0 ) ), 2)

     ;=============================================================================
     ; Use stress tensor to calculate horizontal diffusive tendency source term
     ;=============================================================================

     ; Azimuthally average the 'r' and 'd/dr' expressions
       r_1   = dim_avg_n(r_0(:,:,:,:),2)
       ddr1  = dim_avg_n(ddr0(:,:,:,:),2)

       v_dh0 = ( ( 1/( r_1 ^ 2 * mean_rho0) ) * \
       	       (ddr1 * ( r_1 ^ 2 * mean_rho0 * tau_ry1) ) )

       rhs5a_st = v_dh0	       

     ; Vertical diffusive tendency source term has already been calculated (L827)
     ; Now calculate the azimuthal average using 'dim_avg_n'
       rhs5b_st = dim_avg_n( v_dz0(:,:,:,:),2 )

;=================================================
; Tidy up and produce final arrays for plotting 
;=================================================

     ; Calculating the terms from Montgomery et al. (2018) - their Eq.(4)
     ; Ignore the perturbation PGF term because it is several orders of magnitude smaller

     ; LHS (1) : d{v}/dt          - Mean tangential momentum tendency

     ; RHS (1) : -{u} * {f+vort}  - Mean radial influx of absolute vertical vorticity
     ; RHS (2) : -{w} * (d{v}/dz) - Mean vertical advection of mean tangential momentum
     ; RHS (3) : -{u' * vort'}    - Eddy radial vorticity flux
     ; RHS (4) : -{w' * (dv'/dz)} - Vertical advection of eddy tangential momentum
     ; RHS (5) : {D_v}            - Combined diffusive and PBL tendency

     ; Define new arrays for terms in momentum budget 
       lhs1A  = new((/num_l,dimsizes(radii),numTIMES/),float)   ; d{v}/dt

       rhs1A  = new((/num_l,dimsizes(radii),numTIMES/),float)   ; -{u} * {f+vort}
       rhs2A  = new((/num_l,dimsizes(radii),numTIMES/),float)   ; -{w} * (d{v}/dz)
       rhs3A  = new((/num_l,dimsizes(radii),numTIMES/),float)   ; -{u' * vort'}
       rhs4A  = new((/num_l,dimsizes(radii),numTIMES/),float)   ; -{w' * (dv'/dz)}
       rhs5aA = new((/num_l,dimsizes(radii),numTIMES/),float)   ; {D_r}
       rhs5bA = new((/num_l,dimsizes(radii),numTIMES/),float)   ; {D_z}
       rhsTA  = new((/num_l,dimsizes(radii),numTIMES/),float)   ; Sum of RHS terms

     ; Define additional arrays if calculating averages

       if (ave .eq. 1) then

         lhs1_av0  = new((/num_l,dimsizes(radii)/),float)

         rhs1_av0  = new((/num_l,dimsizes(radii)/),float)
         rhs2_av0  = new((/num_l,dimsizes(radii)/),float)
         rhs3_av0  = new((/num_l,dimsizes(radii)/),float)
         rhs4_av0  = new((/num_l,dimsizes(radii)/),float)
         rhs5a_av0 = new((/num_l,dimsizes(radii)/),float)
         rhs5b_av0 = new((/num_l,dimsizes(radii)/),float)
         rhsT_av0  = new((/num_l,dimsizes(radii)/),float)

       end if

     ;=================================================
     ; Calculate mean tangential wind tendency (LHS)
     ;=================================================

       if (ave .eq. 1) then 

         dt0 = tf0 - ts0 
	
 	 lhs1_av0(:,:) = (mean_vtan0(:,:,tf0) - mean_vtan0(:,:,ts0)) / dt0

       elseif (ave .eq. 0) then 

         do it = ts0, tf0
 
           cp1 = min((/it+1,tf0/))
  	   cm1 = max((/it-1,ts0/))
 
           lhs1A(:,:,it) = (mean_vtan0(:,:,cp1) - mean_vtan0(:,:,cm1)) / 2

         end do		 ; End loop over times (do it = ts0, tf0)

       end if 		 ; End IF statement 


     ; Add metadata and make units consistent for all terms
       lhs1A@units = "m s~S~-1~N~ h~S~-1~N~"

       rhs1A = rhs1_st
       rhs1A = rhs1A * 3600
       rhs1A@units = "m s~S~-1~N~ h~S~-1~N~"

       rhs2A = rhs2_st
       rhs2A = rhs2A * 3600
       rhs2A@units = "m s~S~-1~N~ h~S~-1~N~"
       
       rhs3A = rhs3_st
       rhs3A = rhs3A * 3600
       rhs3A@units = "m s~S~-1~N~ h~S~-1~N~"

       rhs4A = rhs4_st
       rhs4A = rhs4A * 3600
       rhs4A@units = "m s~S~-1~N~ h~S~-1~N~"

       rhs5aA = rhs5a_st
       rhs5aA = rhs5aA * 3600 
       rhs5aA@units = "m s~S~-1~N~ h~S~-1~N~"

       rhs5bA = rhs5b_st
       rhs5bA = rhs5bA * 3600
       rhs5bA@units = "m s~S~-1~N~ h~S~-1~N~"

       hyb_plot   = hybC / 1000
       hyb_plot@units = "km"

     ; Add metadata for the time-averaged tangential wind tendency
       if (ave .eq. 1) then
        lhs1_av0!0    = "lev"
        lhs1_av0&lev  = hyb_plot
        lhs1_av0!1    = "rad"
        lhs1_av0&rad  = radii
       end if

     ; Add metadata for all the other arrays (instantaneous)
       lhs1A!0    = "lev"
       lhs1A&lev  = hyb_plot
       lhs1A!1    = "rad"
       lhs1A&rad  = radii
       lhs1A!2    = "time"

       rhs1A!0    = "lev"
       rhs1A&lev  = hyb_plot
       rhs1A!1    = "rad"
       rhs1A&rad  = radii
       rhs1A!2    = "time"

       rhs2A!0    = "lev"
       rhs2A&lev  = hyb_plot
       rhs2A!1    = "rad"
       rhs2A&rad  = radii
       rhs2A!2    = "time"

       rhs3A!0    = "lev"
       rhs3A&lev  = hyb_plot
       rhs3A!1    = "rad"
       rhs3A&rad  = radii
       rhs3A!2    = "time"

       rhs4A!0    = "lev"
       rhs4A&lev  = hyb_plot
       rhs4A!1    = "rad"
       rhs4A&rad  = radii
       rhs4A!2    = "time"

       rhs5aA!0   = "lev"
       rhs5aA&lev = hyb_plot
       rhs5aA!1   = "rad"
       rhs5aA&rad = radii
       rhs5aA!2   = "time"

       rhs5bA!0   = "lev"
       rhs5bA&lev = hyb_plot
       rhs5bA!1   = "rad"
       rhs5bA&rad = radii
       rhs5bA!2   = "time"

       rhsTA!0    = "lev"
       rhsTA&lev  = hyb_plot
       rhsTA!1    = "rad"
       rhsTA&rad  = radii
       rhsTA!2    = "time"

     ;==============================================
     ; Calculate time averages of all terms above
     ;==============================================

       if (ave .eq. 1) then

        rhs1_av0   = dim_avg_n_Wrap(rhs1A(:,:,ts0:tf0),2)
        rhs2_av0   = dim_avg_n_Wrap(rhs2A(:,:,ts0:tf0),2)
        rhs3_av0   = dim_avg_n_Wrap(rhs3A(:,:,ts0:tf0),2)
        rhs4_av0   = dim_avg_n_Wrap(rhs4A(:,:,ts0:tf0),2)
        rhs5a_av0  = dim_avg_n_Wrap(rhs5aA(:,:,ts0:tf0),2)
        rhs5b_av0  = dim_avg_n_Wrap(rhs5bA(:,:,ts0:tf0),2)
        rhsT_av0   = dim_avg_n_Wrap(rhsTA(:,:,ts0:tf0),2)

       end if

;==============================================
; Switch order of dimensions before plotting
;==============================================

     ; 'gsn_csm_contour' plots the leftmost (rightmost) dimension on the y-axis (x-axis)
     ; Switch 'rad' with 'time' so that 'rad' is the rightmost dimension

     ; Mean tangential wind -- also add metadata
       mean_vtan       = mean_vtan0(lev|:,time|:,rad|:)
       mean_vtan!0     = "lev"
       mean_vtan&lev   = hyb_plot
       mean_vtan!1     = "time"
       mean_vtan!2     = "rad"
       mean_vtan&rad   = radii
       mean_vtan@units = "m s~S~-1~N~"

     ; Mean radial wind -- also add metadata
       mean_vrad       = mean_vrad0(lev|:,time|:,rad|:)
       mean_vrad!0     = "lev"
       mean_vrad&lev   = hyb_plot
       mean_vrad!1     = "time"
       mean_vrad!2     = "rad"
       mean_vrad&rad   = radii
       mean_vrad@units = "m s~S~-1~N~"

     ; Mean vertical velocity -- also add metadata
       mean_vvel       = mean_vvel0(lev|:,time|:,rad|:)
       mean_vvel!0     = "lev"
       mean_vvel&lev   = hyb_plot
       mean_vvel!1     = "time"
       mean_vvel!2     = "rad"
       mean_vvel&rad   = radii
       mean_vvel@units = "m s~S~-1~N~"

     ; Mean absolute angular momentum
       mean_aam        = mean_aam0(lev|:,time|:,rad|:)
       mean_aam!0      = "lev"
       mean_aam&lev    = hyb_plot
       mean_aam!1      = "time"
       mean_aam!2      = "rad"
       mean_aam&rad    = radii
       mean_aam        = mean_aam / (10 ^ 6)
       mean_aam@units  = "10~S~6~N~ m~S~2~N~ s~S~-1~N~"

     ;===================================================
     ; Also calculate time averages of mean variables
     ;===================================================

       if (ave .eq. 1) then

        vrad_av0   = dim_avg_n_Wrap(mean_vrad(:,ts0:tf0,:),1)     ; Mean r-adv of rad mom
        vtan_av0   = dim_avg_n_Wrap(mean_vtan(:,ts0:tf0,:),1)     ; Mean h-adv of eddy mom
        vvel_av0   = dim_avg_n_Wrap(mean_vvel(:,ts0:tf0,:),1)     ; Sum of LHS terms
	aam_av0    = dim_avg_n_Wrap(mean_aam(:,ts0:tf0,:),1)      ; Mean v-adv of mean mom

       end if

     ;==================================
     ; Reorder arrays before plotting
     ;==================================

       lhs1  	  = lhs1A(lev|:,time|:,rad|:)

       rhs1  	  = rhs1A(lev|:,time|:,rad|:)
       rhs2  	  = rhs2A(lev|:,time|:,rad|:)
       rhs3  	  = rhs3A(lev|:,time|:,rad|:)
       rhs4  	  = rhs4A(lev|:,time|:,rad|:)
       rhs5a 	  = rhs5aA(lev|:,time|:,rad|:)
       rhs5b 	  = rhs5bA(lev|:,time|:,rad|:)

       if (ave .eq. 1) then

         lhs1_av  = lhs1_av0(lev|:,rad|:)
         rhs1_av  = rhs1_av0(lev|:,rad|:)
         rhs2_av  = rhs2_av0(lev|:,rad|:)
         rhs3_av  = rhs3_av0(lev|:,rad|:)
         rhs4_av  = rhs4_av0(lev|:,rad|:)
         rhs5a_av = rhs5a_av0(lev|:,rad|:)
         rhs5b_av = rhs5b_av0(lev|:,rad|:)

       end if

;=========================================
; Finally, smooth data before plotting
;=========================================

       lhs1	  = smth9_Wrap(lhs1, 0.5, 0.5, True)
 
       rhs1  	  = smth9_Wrap(rhs1, 0.5, 0.5, True)
       rhs2  	  = smth9_Wrap(rhs2, 0.5, 0.5, True)
       rhs3  	  = smth9_Wrap(rhs3, 0.5, 0.5, True)
       rhs4  	  = smth9_Wrap(rhs4, 0.5, 0.5, True)
       rhs5a 	  = smth9_Wrap(rhs5a, 0.5, 0.5, True)
       rhs5b 	  = smth9_Wrap(rhs5b, 0.5, 0.5, True)

     ; Calculate total RHS tendency
       rhsT       = rhs1 + rhs2 + rhs3 + rhs4 + rhs5a + rhs5b
       rhsT!0     = "lev"
       rhsT&lev   = hyb_plot
       rhsT!1     = "time"
       rhsT!2     = "rad"
       rhsT&rad   = radii
       rhsT@units = "m s~S~-1~N~ h~S~-1~N~"

     ; Calculate sum of all mean terms (vadv + hadv + diff)
       mean       = rhs1 + rhs2 + rhs5a + rhs5b
       mean!0     = "lev"
       mean&lev   = hyb_plot
       mean!1     = "time"
       mean!2     = "rad"
       mean&rad   = radii
       mean@units = "m s~S~-1~N~ h~S~-1~N~"

     ; Calculate sum of all eddy terms (vadv + hadv)
       eddy       = rhs3 + rhs4
       eddy!0     = "lev"
       eddy&lev   = hyb_plot
       eddy!1     = "time"
       eddy!2     = "rad"
       eddy&rad   = radii
       eddy@units = "m s~S~-1~N~ h~S~-1~N~"

      if (ave .eq. 1) then

      ; Total RHS
        rhsT_av       = rhs1_av + rhs2_av + rhs3_av + rhs4_av + rhs5a_av + rhs5b_av
        rhsT_av!0     = "lev"
        rhsT_av&lev   = hyb_plot
        rhsT_av!1     = "rad"
        rhsT_av&rad   = radii
        rhsT_av@units = "m s~S~-1~N~ h~S~-1~N~"

      ; All mean terms
        mean_av       = rhs1_av + rhs2_av + rhs5a_av + rhs5b_av
        mean_av!0     = "lev"
        mean_av&lev   = hyb_plot
        mean_av!1     = "rad"
        mean_av&rad   = radii
        mean_av@units = "m s~S~-1~N~ h~S~-1~N~"

      ; All eddy terms
        eddy_av       = rhs3_av + rhs4_av
        eddy_av!0     = "lev"
        eddy_av&lev   = hyb_plot
        eddy_av!1     = "rad"
        eddy_av&rad   = radii
        eddy_av@units = "m s~S~-1~N~ h~S~-1~N~"

      end if

;==================================================================
; Before plotting, write out arrays to text file (write_matrix)
;==================================================================

      if (wrt .eq. 1) then

      ; Title resources for 'write_matrix' procedure
        opts_mat       = True
        opts_mat@row   = True
        opts_mat@title = ""
        opt_arr        = (/"lhs1", "rhs1", "rhs2", \
                           "rhs3", "rhs4", "rhs5a", "rhs5b"/)
        mbgt_arr       = (/lhs1_av, rhs1_av, rhs2_av, \
                           rhs3_av, rhs4_av, rhs5a_av, rhs5b_av/)

      ; Additional metadata for procedure
        ncol      = dimsizes(hyb_plot)          ; No. of cols (model levels)
        nrow      = dimsizes(radii)             ; No. of rows (radial circles)
        fmtf      = nrow + "f10.3"              ; String format (Fortran style - 10 digits; 3 dpts)

        print_clock("Writing 2D arrays to matrices for further analysis...")

      ; Call procedure and write time-averaged arrays
        do im = 0, 6
          opts_mat@fout = "./mbgt_tan_"+opt_arr(im)+"_"+dat+"_"+ens0+"_"+ts0+"_"+tf0+".txt"
          print_clock("Working on budget term: "+opt_arr(im))
          write_matrix( mbgt_arr(im,:,:), fmtf, opts_mat )
        end do

      ; Now exit script before plotting, to save time
        exit()

      end if 

;============================================
; Define output file path and plot details 
;============================================

;    if (full .eq. 1) then 
;      markers  = new(12,graphic)
;    elseif (full .eq. 0) then 
;      markers  = new(2,graphic)
;    end if

;  ; Add array for markers denoting RMW 
;    mres      	  	    = True
;    mres@gsMarkerIndex	    = 5		 ; Marker symbol
;    mres@gsMarkerSizeF	    = 20.0	 ; Marker size
;    mres@gsMarkerColor	    = "black"    ; Marker colour
;    mres@gsMarkerThicknessF = 5.0 	 ; Marker thickness

;  ; Read in RMW data from text file (first 3 values equal to zero)
;    spd_plot      = new((/2,numTIMES/),float)
;    spd_plot(0,:) = asciiread("$ar/text/"+ens0+"_rmw.txt",(/numTIMES/),"float")  ; x-coords
;    spd_plot(1,:) = 1.0                                                          ; y-coords

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
    opts_mbgt@cnLinesOn                    = False
    opts_mbgt@cnInfoLabelOn                = False
    opts_mbgt@cnLevelSelectionMode         = "ExplicitLevels"

    if (clr.eq.1) then
     opts_mbgt@cnFillPalette               = "ncl_default"
     opts_mbgt@cnFillColors                = (/2,8,16,28,40,52,64,76,88,100,118,-1,-1,\
                                               134,150,166,178,186,194,\
                                               202,214,226,238,250/)
     opts_mbgt@cnLevels                    = (/-150., -100., -50., -20., -10., -5., \
                                                -3., -2., -1., -0.5, -0.2, \
                                                0, 0.2, 0.5, 1., 2., 3., \
                                                5., 10., 20., 50., 100., 150./)
    elseif (clr.eq.2) then
     opts_mbgt@cnFillPalette               = "bdr_extra"
     opts_mbgt@cnFillColors                = (/0,1,2,3,4,5,6,7,8,9,10,-1,-1,\
                                               13,14,15,16,17,18,19,20,21,22,23/)
     opts_mbgt@cnLevels                    = (/-150., -100., -50., -20., -10., -5., \
                                                -3., -2., -1., -0.5, -0.2, \
                                                0, 0.2, 0.5, 1., 2., 3., \
                                                5., 10., 20., 50., 100., 150./)
    elseif (clr.eq.3) then
     opts_mbgt@cnFillPalette               = "amwg256"
     opts_mbgt@cnFillColors                = (/0,6,12,18,27,36,45,54,66,78,90,-1,-1,\
                                               148,156,166,178,186,194,\
                                               202,214,226,238,250/)
     opts_mbgt@cnLevels                    = (/-150., -100., -50., -20., -10., -5., \
                                                -3., -2., -1., -0.5, -0.2, \
                                                0, 0.2, 0.5, 1., 2., 3., \
                                                5., 10., 20., 50., 100., 150./)
    elseif (clr.eq.4) then
     opts_mbgt@cnFillPalette               = "GMT_polar"
     opts_mbgt@cnFillColors               = (/3,5,7,-1,11,14,18/)
     opts_mbgt@cnLevels                   = (/-20,-5,0,0.5,5,20/)
    end if

    opts_mbgt@gsnPaperOrientation          = "landscape"
    opts_mbgt@tiMainString                 = ""
    opts_mbgt@tiMainFontHeightF            = 0.0125
    opts_mbgt@gsnLeftString                = ""
    opts_mbgt@gsnRightString               = ""
    opts_mbgt@gsnMaximize                  = False
    opts_mbgt@lbLabelBarOn                 = False
    opts_mbgt@gsnDraw                      = False
    opts_mbgt@gsnFrame                     = False

  ; Additional plotting resources
    opts_mbgt@tiYAxisString                = "Height (km)"
    opts_mbgt@tiXAxisString                = "Radial distance from cyclone centre (km)"
;    opts_mbgt@tiXAxisFontHeightF           = 0.025   ; X-axis title font height 

    opts_mbgt@tmXTOn                       = "False" ; Turn off top x-axis TM
    opts_mbgt@tmYROn                       = "False" ; Turn off right y-axis TM
    opts_mbgt@gsnMaximize                  = False   ; Maximise plot size (TURN OFF)
    opts_mbgt@gsnAddCyclic                 = False
    opts_mbgt@tmXBMode                     = "Explicit"  ; Set tick marks explicitly
;    opts_mbgt@tmXBLabelFontHeightF         = 0.02    ; Label font height

    if (rad0 .eq. 0) then 
     opts_mbgt@tmXBValues                  = (/0.0,0.25,0.5,0.75,1.0,1.25,1.5/)
     opts_mbgt@tmXBLabels                  = (/"0","25","50","75","100","125","150"/)
    else if (rad0 .eq. 1) then 
     opts_mbgt@tmXBValues                  = (/0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0/)
     opts_mbgt@tmXBLabels                  = (/"0","25","50","75","100","125","150","175","200"/)
    end if 
    end if 

    opts_mbgt@tmYLMode			   = "Explicit"
    opts_mbgt@tmYLValues		   = fspan(2.0, 18.0, 9)
    opts_mbgt@tmYLLabels		   = (/"2.0","4.0","6.0","8.0","10.0",\
    					       "12.0","14.0","16.0","18.0"/)
;    opts_mbgt@tmYLLabelFontHeightF	   = 0.02
    opts_mbgt@gsnYAxisIrregular2Linear	   = True	; Linearise the y-axis 

    opts_mbgt@trYMinF			   = 0.0	; y-axis minimum value
    opts_mbgt@trYMaxF			   = 17.6	; y-axis maximum value 

  ; Resources to overlay radial wind (following M18)
    opts_vrad                              = True
    opts_vrad@cnFillOn                     = False
    opts_vrad@cnLineColor                  = "black"
    opts_vrad@cnLevelSelectionMode         = "ExplicitLevels"
    opts_vrad@cnLevels                     = (/-1.2, 1.2/)
    opts_vrad@gsnContourNegLineDashPattern = 16
    opts_vrad@gsnContourPosLineDashPattern = 0
    opts_vrad@cnLineThicknessF             = 2.5

    opts_vrad@cnInfoLabelOn                = False
    opts_vrad@gsnDraw                      = False
    opts_vrad@gsnFrame                     = False
    opts_vrad@cnLineLabelsOn               = False
    opts_vrad@tiMainString                 = ""
    opts_vrad@gsnLeftString                = ""
    opts_vrad@gsnRightString               = ""

  ; Resources to overlay vertical velocity (following M18)
    opts_vvel                              = True
    opts_vvel@cnFillOn                     = False

  ; Set contour line colour to match M18 if we choose their colour scheme
    if (clr .eq. 4) then
     opts_vvel@cnLineColor                 = "yellow"
    else
     opts_vvel@cnLineColor                 = "grey45"
    end if

    opts_vvel@cnLevelSelectionMode         = "ExplicitLevels"
    opts_vvel@cnLevels                     = (/0.5/)
    opts_vvel@cnLineThicknessF             = 5.0

    opts_vvel@cnInfoLabelOn                = False
    opts_vvel@gsnDraw                      = False
    opts_vvel@gsnFrame                     = False
    opts_vvel@cnLineLabelsOn               = False
    opts_vvel@tiMainString                 = ""
    opts_vvel@gsnLeftString                = ""
    opts_vvel@gsnRightString               = ""

  ; Resources to overlay zero line
    opts_zero                             = True
    opts_zero@cnFillOn                    = False
    opts_zero@cnLineColor                 = "royalblue1"
    opts_zero@cnLevelSelectionMode        = "ExplicitLevels"
    opts_zero@cnLevels                    = (/-5.0, 0.0, 5.0/)
    opts_zero@cnLineThicknessF            = 2.0
    opts_zero@cnInfoLabelOn               = False
    opts_zero@gsnDraw                     = False       ; Do not draw the plot
    opts_zero@gsnFrame                    = False       ; Do no advance the frame
    opts_zero@gsnContourZeroLineThicknessF = 2.0
    opts_zero@cnLineLabelsOn               = False      ; Turn off line labels
    opts_zero@tiMainString                = ""
    opts_zero@gsnLeftString               = ""
    opts_zero@gsnRightString              = ""

;==================================================================
; Also add resources for mean quantities (u,v,w) when 'full' = 2
;==================================================================

  ; Reading in colour table using 'cnFillPalette' removes the first 2 colours (black/white)
  ; This means that the indices are all 2 values lower than they'd normally be

  ; Radial wind
    opts_vrad2                             = opts_mbgt
    delete(opts_vrad2@cnLevels)
    delete(opts_vrad2@cnFillColors)
    opts_vrad2@cnFillPalette               = "BlueDarkRed18"
    opts_vrad2@cnLevels                    = (/-10., -5., -3., -2., -1., -0.5, -0.2, \
                                                0, 0.2, 0.5, 1., 2., 3., 5., 10./)
    opts_vrad2@cnFillColors                = (/0,2,3,4,5,6,7,-1,-1,10,11,12,13,14,15,17/)
    opts_vrad2@lbLabelBarOn                = True
    opts_vrad2@lbBoxEndCapStyle            = "TriangleBothEnds"
    opts_vrad2@lbOrientation               = "Vertical"
    opts_vrad2@vpXF                        = 0.07
    opts_vrad2@vpYF                        = 0.90
    opts_vrad2@vpWidthF                    = 0.15
    opts_vrad2@vpHeightF                   = 0.15

  ; Tangential wind
    opts_vtan2                             = opts_mbgt
    delete(opts_vtan2@cnLevels)
    delete(opts_vtan2@cnFillColors)
    opts_vtan2@cnFillPalette               = "radar_new"

    if (new0 .eq. 1) then 
      opts_vtan2@cnLevels                    = (/5.0, 10.0, 15.0, 20.0,\
                                                25.0, 30.0, 32.5, 35.0, \
                                                37.5, 40.0, 42.5, 45.0, \
                                                50.0/)
      opts_vtan2@cnFillColors                = (/-1,0,1,2,3,4,5,6,7,8,9,10,11,12/)
    else 
      opts_vtan2@cnLevels                    = (/5.0, 10.0, 15.0, 20.0,\
                                                25.0, 30.0, 35.0, 40.0, \
                                                45.0, 50.0, 55.0, 60.0, \
                                                65.0/)
      opts_vtan2@cnFillColors                = (/-1,0,1,2,3,4,5,6,7,8,9,10,11,12/)
    end if 

    opts_vtan2@lbLabelBarOn                = True
    opts_vtan2@lbBoxEndCapStyle            = "TriangleBothEnds"
    opts_vtan2@lbOrientation               = "Vertical"
    opts_vtan2@vpXF                        = 0.31
    opts_vtan2@vpYF                        = 0.90
    opts_vtan2@vpWidthF                    = 0.15
    opts_vtan2@vpHeightF                   = 0.15

  ; Vertical velocity
    opts_vvel2                             = opts_mbgt
    delete(opts_vvel2@cnLevels)
    delete(opts_vvel2@cnFillColors)
    opts_vvel2@cnFillPalette               = "n11"
    opts_vvel2@cnLevels                    = (/0.1, 0.2, 0.3, 0.4, 0.5, 0.6, \
                                               0.7, 0.8, 1.0, 1.2, 1.5/)
    opts_vvel2@cnFillColors                = (/-1,0,1,2,3,4,5,6,7,8,9,10/)
    opts_vvel2@lbLabelBarOn                = True
    opts_vvel2@lbBoxEndCapStyle            = "TriangleBothEnds"
    opts_vvel2@lbOrientation               = "Vertical"
    opts_vvel2@vpXF                        = 0.55
    opts_vvel2@vpYF                        = 0.90
    opts_vvel2@vpWidthF                    = 0.15
    opts_vvel2@vpHeightF                   = 0.15

  ; Absolute angular momentum
    opts_aam                               = opts_vtan2
    delete(opts_aam@cnLevels)
    delete(opts_aam@cnFillColors)

    if (new0 .eq. 1) then
      opts_aam@cnLevels                      = (/1.0, 1.5, 2.0, 2.5, 3.0, \
                                                 3.5, 4.0, 4.5, 5.0/)
      opts_aam@cnFillColors                  = (/-1,3,4,5,6,7,8,9,10,11/)
    else
      opts_aam@cnLevels                      = (/1.0, 1.5, 2.0, 2.5, 3.0, \
                                                 3.5, 4.0, 4.5, 5.0, \
                                                 5.5, 6.0, 6.5/)
      opts_aam@cnFillColors                  = (/-1,0,1,2,3,4,5,6,7,8,9,10,11/)
    end if 

    opts_aam@lbLabelBarOn                  = True
    opts_aam@lbBoxEndCapStyle              = "TriangleBothEnds"
    opts_aam@lbOrientation                 = "Vertical"
    opts_aam@vpXF                          = 0.79
    opts_aam@vpYF                          = 0.90
    opts_aam@vpWidthF                      = 0.15
    opts_aam@vpHeightF                     = 0.15

;===============================================================
; Resources to overlay panel labels without using 'gsn_panel'
;===============================================================

    if (full .eq. 1) then 
     txid                       = new(12,graphic)
     amid                       = new(12,graphic)
    elseif (full .eq. 0) then 
     txid  	      	        = new(2,graphic)
     amid			= new(2,graphic)
    end if 

  ; Label text resources
    txres                       = True
    txres@txPerimOn             = True
    txres@txFontHeightF         = 0.010
    txres@txBackgroundFillColor = "White"

  ; Label position resources
    amres                       = True
    amres@amParallelPosF        = 0.5                   ; Right edge
    amres@amOrthogonalPosF      = -0.5                  ; Top edge
    amres@amJust                = "TopRight"

    if (full .eq. 1) then 
     panel_strings              = (/"a) vrad", "b) vtan", "c) vvel", "d) aam", \
                                    "e) V_mz", "f) V_mv", "g) V_ez", \
                                    "h) V_ev", "i) V_dh", "j) V_dz", \
                                    "k) All V tend", "l) V_t"/)
    elseif (full .eq. 0) then 
     panel_strings    	        = (/"a) Mean terms","b) Eddy terms"/)
    end if 

;=================
; Plot the data 
;=================

  ; Calculating the terms from Montgomery et al. (2018) - their Eq.(4)
  ; Ignore the perturbation PGF term because it is several orders of magnitude smaller

  ; LHS (1) : d{v}/dt          - Mean tangential momentum tendency

  ; RHS (1) : -{u} * {f+vort}  - Mean radial influx of absolute vertical vorticity
  ; RHS (2) : -{w} * (d{v}/dz) - Mean vertical advection of mean tangential momentum
  ; RHS (3) : -{u' * vort'}    - Eddy radial vorticity flux
  ; RHS (4) : -{w' * (dv'/dz)} - Vertical advection of eddy tangential momentum
  ; RHS (5) : {D_v}            - Combined diffusive and PBL tendency

  if (ave .eq. 0) then 

   do it = ts0, tf0
   
    print("Working on time: "+it)

  ; Output file location and type 
    output = "$sam/nepartak/images/cyl_coords/mbgt_1deg_tan2_xz_"+dat+"_"+ens0+"_"+time_arr(it)
    wks    = gsn_open_wks(opt,output)

  if (full .eq. 1) then 

  ; Radial wind (filled contours)
    panel_vrad  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad2) ; [u]
    panel_zero  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_zero)  ; Zero line contour
    zero0       = ColorNegDashZeroPosContour(panel_zero,"transparent","black","transparent")
    overlay(panel_vrad, zero0)
    txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
    amid(0)     = gsn_add_annotation(panel_vrad, txid(0), amres)

  ; Tangential wind (filled contours)
    panel_vtan  = gsn_csm_contour(wks,mean_vtan(:,it,:),opts_vtan2) ; [v]
    txid(1)     = gsn_create_text(wks, panel_strings(1), txres)
    amid(1)     = gsn_add_annotation(panel_vtan, txid(1), amres)

  ; Vertical velocity (filled contours)
    panel_vvel  = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel2) ; [w]
    txid(2)     = gsn_create_text(wks, panel_strings(2), txres)
    amid(2)     = gsn_add_annotation(panel_vvel, txid(2), amres)

  ; Absolute angular momentum (filled contours)
    panel_aam   = gsn_csm_contour(wks,mean_aam(:,it,:),opts_aam)   ; [AAM]
    txid(3)     = gsn_create_text(wks, panel_strings(3), txres)
    amid(3)     = gsn_add_annotation(panel_aam, txid(3), amres)


  ; Mean radial influx of absolute vertical vorticity [V_mzeta]
    opts_mbgt@vpXF       = 0.07
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel0      = gsn_csm_contour(wks,rhs1(:,it,:),opts_mbgt)	    ; [V_mzeta]
    cont0	= gsn_csm_contour(wks,rhs1(:,it,:),opts_zero)	    ; Draw contours
    cont        = ColorNegDashZeroPosContour(cont0,"royalblue1","grey50","IndianRed")

    panel0_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad)  ; [u]
    panel0_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel)  ; [w]

    overlay(panel0, cont)                                          ; Overlay contours
    overlay(panel0, panel0_ovr)
    overlay(panel0, panel0_ovr2)
    txid(4)     = gsn_create_text(wks, panel_strings(4), txres)
    amid(4)     = gsn_add_annotation(panel0, txid(4), amres)
    delete([/cont0,cont/])

  ; Mean vertical advection of mean tangential momentum [V_mv]
    opts_mbgt@vpXF       = 0.31
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel1     = gsn_csm_contour(wks,rhs2(:,it,:),opts_mbgt)        ; [V_mv]
    cont0      = gsn_csm_contour(wks,rhs2(:,it,:),opts_zero)       ; Draw contours
    cont        = ColorNegDashZeroPosContour(cont0,"royalblue1","grey50","IndianRed")

    panel1_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad)  ; [u]
    panel1_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel)  ; [w]

    overlay(panel1, cont)                                          ; Overlay contours
    overlay(panel1, panel1_ovr)
    overlay(panel1, panel1_ovr2)
    txid(5)     = gsn_create_text(wks, panel_strings(5), txres)
    amid(5)     = gsn_add_annotation(panel1, txid(5), amres)
    delete([/cont0,cont/])

  ; Eddy radial vorticity flux [V_ezeta]
    opts_mbgt@vpXF       = 0.55
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel2     = gsn_csm_contour(wks,rhs3(:,it,:),opts_mbgt)        ; [V_ezeta]
    cont0      = gsn_csm_contour(wks,rhs3(:,it,:),opts_zero)        ; Draw contours
    cont        = ColorNegDashZeroPosContour(cont0,"royalblue1","grey50","IndianRed")

    panel2_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad)  ; [u]
    panel2_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel)  ; [w]

    overlay(panel2, cont)                                          ; Overlay contours
    overlay(panel2, panel2_ovr)
    overlay(panel2, panel2_ovr2)
    txid(6)     = gsn_create_text(wks, panel_strings(6), txres)
    amid(6)     = gsn_add_annotation(panel2, txid(6), amres)
    delete([/cont0,cont/])

  ; Vertical advection of eddy tangential momentum [V_ev]
    opts_mbgt@vpXF       = 0.79
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel3     = gsn_csm_contour(wks,rhs4(:,it,:),opts_mbgt)        ; [V_ev]
    cont0      = gsn_csm_contour(wks,rhs4(:,it,:),opts_zero)        ; Draw contours
    cont        = ColorNegDashZeroPosContour(cont0,"royalblue1","grey50","IndianRed")

    panel3_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad)  ; [u]
    panel3_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel)  ; [w]

    overlay(panel3, cont)                                          ; Overlay contours
    overlay(panel3, panel3_ovr)
    overlay(panel3, panel3_ovr2)
    txid(7)     = gsn_create_text(wks, panel_strings(7), txres)
    amid(7)     = gsn_add_annotation(panel3, txid(7), amres)
    delete([/cont0,cont/])

  ; Horizontal diffusive tendency of tangential momentum [V_dh]
    opts_mbgt@vpXF       = 0.07
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel4     = gsn_csm_contour(wks,rhs5a(:,it,:),opts_mbgt)       ; [V_dh]
    cont0      = gsn_csm_contour(wks,rhs5a(:,it,:),opts_zero)       ; Draw contours
    cont        = ColorNegDashZeroPosContour(cont0,"royalblue1","grey50","IndianRed")

    panel4_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad)  ; [u]
    panel4_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel)  ; [w]

    overlay(panel4, cont)                                          ; Overlay contours
    overlay(panel4, panel4_ovr)
    overlay(panel4, panel4_ovr2)
    txid(8)     = gsn_create_text(wks, panel_strings(8), txres)
    amid(8)     = gsn_add_annotation(panel4, txid(8), amres)
    delete([/cont0,cont/])

  ; Vertical diffusive tendency of tangential momentum [V_dz]
    opts_mbgt@vpXF       = 0.31
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel5     = gsn_csm_contour(wks,rhs5b(:,it,:),opts_mbgt)       ; [V_dz]
    cont0      = gsn_csm_contour(wks,rhs5b(:,it,:),opts_zero)       ; Draw contours
    cont        = ColorNegDashZeroPosContour(cont0,"royalblue1","grey50","IndianRed")
    
    panel5_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad)  ; [u]
    panel5_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel)  ; [w]

    overlay(panel5, cont)                                          ; Overlay contours
    overlay(panel5, panel5_ovr)
    overlay(panel5, panel5_ovr2)
    txid(9)     = gsn_create_text(wks, panel_strings(9), txres)
    amid(9)     = gsn_add_annotation(panel5, txid(9), amres)
    delete([/cont0,cont/])

  ; Total RHS [All V tend]
    opts_mbgt@vpXF       = 0.55
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel6      = gsn_csm_contour(wks,rhsT(:,it,:),opts_mbgt)        ; [All V tend]
    cont0       = gsn_csm_contour(wks,rhsT(:,it,:),opts_zero)        ; Draw contours
    cont        = ColorNegDashZeroPosContour(cont0,"royalblue1","grey50","IndianRed")

    panel6_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad)  ; [u]
    panel6_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel)  ; [w]

    overlay(panel6, cont)                                          ; Overlay contours
    overlay(panel6, panel6_ovr)
    overlay(panel6, panel6_ovr2)
    txid(10)    = gsn_create_text(wks, panel_strings(10), txres)
    amid(10)    = gsn_add_annotation(panel6, txid(10), amres)
    delete([/cont0,cont/])

  ; Local tendency of the mean tangential wind [V_t]
    opts_mbgt@vpXF       = 0.79
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel7     = gsn_csm_contour(wks,lhs1(:,it,:),opts_mbgt)        ; [V_t]
    cont0      = gsn_csm_contour(wks,lhs1(:,it,:),opts_zero)        ; Draw contours
    cont        = ColorNegDashZeroPosContour(cont0,"royalblue1","grey50","IndianRed")

    panel7_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad)  ; [u]
    panel7_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel)  ; [w]


    overlay(panel7, cont)                                          ; Overlay contours
    overlay(panel7, panel7_ovr)
    overlay(panel7, panel7_ovr2)
    txid(11)    = gsn_create_text(wks, panel_strings(11), txres)
    amid(11)    = gsn_add_annotation(panel7, txid(11), amres)
    delete([/cont0,cont/])

  elseif (full .eq. 0) then

  ; All eddy terms
    opts_mbgt@vpXF       = 0.05
    opts_mbgt@vpYF       = 0.75
    opts_mbgt@vpWidthF   = 0.42
    opts_mbgt@vpHeightF  = 0.42
    panel_eddy      = gsn_csm_contour(wks,eddy(:,it,:),opts_mbgt)      ; [All U tend]
    cont0      	    = gsn_csm_contour(wks,eddy(:,it,:),opts_zero)        ; Draw contours
    cont        = ColorNegDashZeroPosContour(cont0,"royalblue1","grey50","IndianRed")

    panel_eddy_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad) ; [U]
    panel_eddy_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel) ; [w]

    overlay(panel_eddy, cont)                                          ; Overlay contours
    overlay(panel_eddy, panel_eddy_ovr)
    overlay(panel_eddy, panel_eddy_ovr2)
    txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
    amid(0)     = gsn_add_annotation(panel_eddy, txid(0), amres)
    delete([/cont0,cont/])

  ; All mean terms
    opts_mbgt@vpXF       = 0.55
    opts_mbgt@vpYF       = 0.75
    opts_mbgt@vpWidthF   = 0.42
    opts_mbgt@vpHeightF  = 0.42
    panel_mean      = gsn_csm_contour(wks,mean(:,it,:),opts_mbgt)      ; [All U tend]
    cont0      	    = gsn_csm_contour(wks,mean(:,it,:),opts_zero)        ; Draw contours
    cont        = ColorNegDashZeroPosContour(cont0,"royalblue1","grey50","IndianRed")

    panel_mean_ovr  = gsn_csm_contour(wks,mean_vrad(:,it,:),opts_vrad) ; [U]
    panel_mean_ovr2 = gsn_csm_contour(wks,mean_vvel(:,it,:),opts_vvel) ; [w]

    overlay(panel_mean, cont)                                          ; Overlay contours
    overlay(panel_mean, panel_mean_ovr)
    overlay(panel_mean, panel_mean_ovr2)
    txid(1)     = gsn_create_text(wks, panel_strings(1), txres)
    amid(1)     = gsn_add_annotation(panel_mean, txid(1), amres)
    delete([/cont0,cont/])

  end if

;=======================================
; Overlay RMW markers onto panel plot
;=======================================

;  ; First define an array containing all of the plots above 
;    if (full .eq. 1) then
;     panel = (/panel_vrad, panel_vtan, panel_vvel, panel_aam, \
;     	       panel0, panel1, panel2, panel3, panel4, panel5, panel6, panel7/)
;    elseif (full .eq. 0) then
;     panel = (/panel_mean, panel_eddy/)
;    end if 

;  ; Overlay markers on each panel corresponding to RMW
;    if (full .eq. 1) then
;      do m = 0, 11
;        markers(m) = gsn_add_polymarker(wks,panel(m),spd_plot(0,it),spd_plot(1,it),mres)
;      end do
;    elseif (full .eq. 0) then
;      do m = 0,  1
;        markers(m) = gsn_add_polymarker(wks,panel(m),spd_plot(0,it),spd_plot(1,it),mres)
;      end do
;    end if

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    pres			= True		; Resource to call function below
    maximize_output(wks,pres)			; Calls 'draw' and 'frame'

;    delete([/panel,markers/])

   end do     ; End time loop (do it = 0, ts, tf)

  elseif (ave .eq. 1) then 

  ; Define output file path and type
    if (full .eq. 0) then
     output = "$sam/nepartak/images/cyl_coords/mbgt_1deg_tan_eddy_mean_xz_"+dat+\
     	      "_"+ens0+"_"+ts0+"_"+tf0
    elseif (full .eq. 1) then 
     output = "$sam/nepartak/images/cyl_coords/mbgt_1deg_tan_ave_xz_"+dat+"_"+ens0+"_"+ts0+"_"+tf0
    end if

    wks = gsn_open_wks(opt,output)

   if (full .eq. 1) then 

  ; Radial wind (filled contours)
    panel_vrad  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad2) ; [u]
    panel_zero  = gsn_csm_contour(wks,vrad_av0(:,:),opts_zero)  ; Zero line contour
    zero0       = ColorNegDashZeroPosContour(panel_zero,"transparent","black","transparent")
    overlay(panel_vrad, zero0)
    txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
    amid(0)     = gsn_add_annotation(panel_vrad, txid(0), amres)

  ; Tangential wind (filled contours)
    panel_vtan  = gsn_csm_contour(wks,vtan_av0(:,:),opts_vtan2) ; [v]
    txid(1)     = gsn_create_text(wks, panel_strings(1), txres)
    amid(1)     = gsn_add_annotation(panel_vtan, txid(1), amres)

  ; Vertical velocity (filled contours)
    panel_vvel  = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel2) ; [w]
    txid(2)     = gsn_create_text(wks, panel_strings(2), txres)
    amid(2)     = gsn_add_annotation(panel_vvel, txid(2), amres)

  ; Absolute angular momentum (filled contours)
    panel_aam   = gsn_csm_contour(wks,aam_av0(:,:),opts_aam)    ; [AAM]
    txid(3)     = gsn_create_text(wks, panel_strings(3), txres)
    amid(3)     = gsn_add_annotation(panel_aam, txid(3), amres)


  ; Mean radial influx of absolute vertical vorticity [V_mzeta]
    opts_mbgt@vpXF       = 0.07
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel0      = gsn_csm_contour(wks,rhs1_av(:,:),opts_mbgt)   ; [V_mzeta]
    cont0       = gsn_csm_contour(wks,rhs1_av(:,:),opts_zero)   ; Draw contours
    cont        = ColorNegDashZeroPosContour(cont0,"royalblue1","grey50","IndianRed")

    panel0_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad)  ; [u]
    panel0_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel)  ; [w]

    overlay(panel0, cont)                                          ; Overlay contours
    overlay(panel0, panel0_ovr)
    overlay(panel0, panel0_ovr2)
    txid(4)     = gsn_create_text(wks, panel_strings(4), txres)
    amid(4)     = gsn_add_annotation(panel0, txid(4), amres)
    delete([/cont0,cont/])

  ; Mean vertical advection of mean tangential momentum [V_mv]
    opts_mbgt@vpXF       = 0.31
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel1     = gsn_csm_contour(wks,rhs2_av(:,:),opts_mbgt)    ; [V_mv]
    cont0      = gsn_csm_contour(wks,rhs2_av(:,:),opts_zero)    ; Draw contours
    cont        = ColorNegDashZeroPosContour(cont0,"royalblue1","grey50","IndianRed")

    panel1_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad)  ; [u]
    panel1_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel)  ; [w]

    overlay(panel1, cont)                                          ; Overlay contours
    overlay(panel1, panel1_ovr)
    overlay(panel1, panel1_ovr2)
    txid(5)     = gsn_create_text(wks, panel_strings(5), txres)
    amid(5)     = gsn_add_annotation(panel1, txid(5), amres)
    delete([/cont0,cont/])

  ; Eddy radial vorticity flux [V_ezeta]
    opts_mbgt@vpXF       = 0.55
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel2     = gsn_csm_contour(wks,rhs3_av(:,:),opts_mbgt)    ; [V_ezeta]
    cont0      = gsn_csm_contour(wks,rhs3_av(:,:),opts_zero)    ; Draw contours
    cont        = ColorNegDashZeroPosContour(cont0,"royalblue1","grey50","IndianRed")

    panel2_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad)  ; [u]
    panel2_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel)  ; [w]

    overlay(panel2, cont)                                          ; Overlay contours
    overlay(panel2, panel2_ovr)
    overlay(panel2, panel2_ovr2)
    txid(6)     = gsn_create_text(wks, panel_strings(6), txres)
    amid(6)     = gsn_add_annotation(panel2, txid(6), amres)
    delete([/cont0,cont/])

  ; Vertical advection of eddy tangential momentum [V_ev]
    opts_mbgt@vpXF       = 0.79
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel3     = gsn_csm_contour(wks,rhs4_av(:,:),opts_mbgt)    ; [V_ev]
    cont0      = gsn_csm_contour(wks,rhs4_av(:,:),opts_zero)    ; Draw contours
    cont        = ColorNegDashZeroPosContour(cont0,"royalblue1","grey50","IndianRed")

    panel3_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad)  ; [u]
    panel3_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel)  ; [w]

    overlay(panel3, cont)                                          ; Overlay contours
    overlay(panel3, panel3_ovr)
    overlay(panel3, panel3_ovr2)
    txid(7)     = gsn_create_text(wks, panel_strings(7), txres)
    amid(7)     = gsn_add_annotation(panel3, txid(7), amres)
    delete([/cont0,cont/])

  ; Horizontal diffusive tendency of tangential momentum [V_dh]
    opts_mbgt@vpXF       = 0.07
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel4     = gsn_csm_contour(wks,rhs5a_av(:,:),opts_mbgt)   ; [V_dh]
    cont0      = gsn_csm_contour(wks,rhs5a_av(:,:),opts_zero)   ; Draw contours
    cont        = ColorNegDashZeroPosContour(cont0,"royalblue1","grey50","IndianRed")

    panel4_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad)  ; [u]
    panel4_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel)  ; [w]

    overlay(panel4, cont)                                          ; Overlay contours
    overlay(panel4, panel4_ovr)
    overlay(panel4, panel4_ovr2)
    txid(8)     = gsn_create_text(wks, panel_strings(8), txres)
    amid(8)     = gsn_add_annotation(panel4, txid(8), amres)
    delete([/cont0,cont/])

  ; Vertical diffusive tendency of tangential momentum [V_dz]
    opts_mbgt@vpXF       = 0.31
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel5     = gsn_csm_contour(wks,rhs5b_av(:,:),opts_mbgt)   ; [V_dz]
    cont0      = gsn_csm_contour(wks,rhs5b_av(:,:),opts_zero)   ; Draw contours
    cont        = ColorNegDashZeroPosContour(cont0,"royalblue1","grey50","IndianRed")

    panel5_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad)  ; [u]
    panel5_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel)  ; [w]

    overlay(panel5, cont)                                          ; Overlay contours
    overlay(panel5, panel5_ovr)
    overlay(panel5, panel5_ovr2)
    txid(9)     = gsn_create_text(wks, panel_strings(9), txres)
    amid(9)     = gsn_add_annotation(panel5, txid(9), amres)
    delete([/cont0,cont/])

  ; Total RHS [All V tend]
    opts_mbgt@vpXF       = 0.55
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel6      = gsn_csm_contour(wks,rhsT_av(:,:),opts_mbgt)   ; [All V tend]
    cont0       = gsn_csm_contour(wks,rhsT_av(:,:),opts_zero)   ; Draw contours
    cont        = ColorNegDashZeroPosContour(cont0,"royalblue1","grey50","IndianRed")

    panel6_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad)  ; [u]
    panel6_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel)  ; [w]

    overlay(panel6, cont)                                          ; Overlay contours
    overlay(panel6, panel6_ovr)
    overlay(panel6, panel6_ovr2)
    txid(10)    = gsn_create_text(wks, panel_strings(10), txres)
    amid(10)    = gsn_add_annotation(panel6, txid(10), amres)
    delete([/cont0,cont/])

  ; Local tendency of the mean tangential wind [V_t]
    opts_mbgt@vpXF       = 0.79
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel7     = gsn_csm_contour(wks,lhs1_av(:,:),opts_mbgt)    ; [V_t]
    cont0      = gsn_csm_contour(wks,lhs1_av(:,:),opts_zero)    ; Draw contours
    cont        = ColorNegDashZeroPosContour(cont0,"royalblue1","grey50","IndianRed")

    panel7_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad)  ; [u]
    panel7_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel)  ; [w]

    overlay(panel7, cont)                                          ; Overlay contours
    overlay(panel7, panel7_ovr)
    overlay(panel7, panel7_ovr2)
    txid(11)    = gsn_create_text(wks, panel_strings(11), txres)
    amid(11)    = gsn_add_annotation(panel7, txid(11), amres)
    delete([/cont0,cont/])

   elseif (full .eq. 0) then

  ; All eddy terms
    opts_mbgt@vpXF       = 0.05
    opts_mbgt@vpYF       = 0.75
    opts_mbgt@vpWidthF   = 0.42
    opts_mbgt@vpHeightF  = 0.42
    panel_eddy      = gsn_csm_contour(wks,eddy_av(:,:),opts_mbgt)  ; [All U tend]
    cont0      	    = gsn_csm_contour(wks,eddy_av(:,:),opts_zero)  ; Draw contours
    cont        = ColorNegDashZeroPosContour(cont0,"royalblue1","grey50","IndianRed")

    panel_eddy_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad) ; [U]
    panel_eddy_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel) ; [w]

    overlay(panel_eddy, cont)                                          ; Overlay contours
    overlay(panel_eddy, panel_eddy_ovr)
    overlay(panel_eddy, panel_eddy_ovr2)
    txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
    amid(0)     = gsn_add_annotation(panel_eddy, txid(0), amres)
    delete([/cont0,cont/])

  ; All mean terms
    opts_mbgt@vpXF       = 0.55
    opts_mbgt@vpYF       = 0.75
    opts_mbgt@vpWidthF   = 0.42
    opts_mbgt@vpHeightF  = 0.42
    panel_mean      = gsn_csm_contour(wks,mean_av(:,:),opts_mbgt)  ; [All U tend]
    cont0      	    = gsn_csm_contour(wks,mean_av(:,:),opts_zero)  ; Draw contours
    cont        = ColorNegDashZeroPosContour(cont0,"royalblue1","grey50","IndianRed")

    panel_mean_ovr  = gsn_csm_contour(wks,vrad_av0(:,:),opts_vrad) ; [U]
    panel_mean_ovr2 = gsn_csm_contour(wks,vvel_av0(:,:),opts_vvel) ; [w]

    overlay(panel_mean, cont)                                          ; Overlay contours
    overlay(panel_mean, panel_mean_ovr)
    overlay(panel_mean, panel_mean_ovr2)
    txid(1)     = gsn_create_text(wks, panel_strings(1), txres)
    amid(1)     = gsn_add_annotation(panel_mean, txid(1), amres)
    delete([/cont0,cont/])

   end if

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    pres                        = True
    maximize_output(wks,pres)                              ; Calls 'draw' and 'frame'

  end if      ; End IF statement (if (ave .eq. 0) then ...)

end