; Calculate composite plots of relative vorticity and strain in different VC phases

; Run using:

; ncl opt=\"x11\" dist=2.0 nt=73 nr=11 ar=0.5 typ=\"sym\" cn0=\"geo_sm\"
; lay=1 mlev0=13 mlev1=25 sc=5 sm=0 rad0=0.5 dnt=0 plt=0 grp=\"7\" var=\"va\" tc_ang_vel.ncl

; 'opt'   = output file format ("pdf" or "x11")
; 'dist'  = size of box following storm (degrees)
; 'typ'	  = VC phase ("sym","asym","stoa","atos")
; 'full'  = top 18 (1); top 8 (-5); top 6 (-6)
; 'ar'    = distance in degrees of outer radius (1.0, 2.0, etc)
; 'nr'    = number of radial circles between r = 0 and r = 'ar' (21, 41, 61, etc)
; 'nt'    = number of azimuth angles in cylindrical grid
; 'cn0'   = offline storm centre: "slp", "vort", "geo", "geo_sm"
; 'rad0'  = radius (degrees) for final plots 

; 'pc' stream --> 'it'
; 'centre'    --> 'it-1'
; 'pd' stream --> 'it-2'

; Load main NCL functions and procedures
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"

; Also load user-defined functions and procedures
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_centre.ncl"
load "$sam/ncl_func/setup_cyl.ncl"
load "$sam/ncl_func/draw_cyl.ncl"

begin

;==============================
; Prelash (reading in files)
;==============================

; Group 5; Group 6; Group 7; Group 8; Group 9 (28/10/2019)
  vc_arr_all = (/ (/56,38,51,30,22/), (/55,35,49,27,23/), \
                  (/69,39,64,37,27/), (/52,33,43,26,19/) /)

; Link the chosen group to the way we read in data below
  if (grp .eq. "5") then
   vc_ind = 0
  elseif (grp .eq. "6") then
   vc_ind = 1
  elseif (grp .eq. "7") then
   vc_ind = 2
  elseif (grp .eq. "8") then
   vc_ind = 3
  elseif (grp .eq. "9") then
   vc_ind = 4
  end if

; Read in data from text files according to VC phase  
  if (typ .eq. "sym") then
   typ0    = "ring"
   gr_ind  = 0
   out_str = "Ringlike phase"
  elseif (typ .eq. "s2a") then
   typ0    = "r2m"
   gr_ind  = 1
   out_str = "Ringlike to monopole transition"
  elseif (typ .eq. "asym") then
   typ0    = "mono"
   gr_ind  = 2
   out_str = "Monopole phase"
  elseif (typ .eq. "a2s") then
   typ0    = "m2r"
   gr_ind  = 3
   out_str = "Monopole to ringlike transition"
  elseif (typ .eq. "novc") then
   typ0 = "novc"
   out_str = "No vacillation cycles"
  end if

; String containing info on width of radial averaging band  
  ar0  = ar * 100
  rstr = sprintf("%0.0fkm", ar0)

; Which set of simulations do we want to analyse?

; No VCs
  if (typ .eq. "novc") then
    novc_size = 16
    dat  = asciiread("$ar/text/novc_dat.txt",(/novc_size/),"string")
    ens0 = asciiread("$ar/text/novc_sim.txt",(/novc_size/),"string")
    ts0  = asciiread("$ar/text/novc_ts.txt",(/novc_size/),"integer")
    tf0  = asciiread("$ar/text/novc_tf.txt",(/novc_size/),"integer")
; Any other VC phase
  else
   dat =asciiread("$ar/text/vc_"+typ+"_dat_group"+grp+".txt",(/vc_arr_all(gr_ind,vc_ind)/),"string")
   ens0=asciiread("$ar/text/vc_"+typ+"_sim_group"+grp+".txt",(/vc_arr_all(gr_ind,vc_ind)/),"string")
   ts0 =asciiread("$ar/text/vc_"+typ+"_ts_group"+grp+".txt",(/vc_arr_all(gr_ind,vc_ind)/),"integer")
   tf0 =asciiread("$ar/text/vc_"+typ+"_tf_group"+grp+".txt",(/vc_arr_all(gr_ind,vc_ind)/),"integer")
  end if

;=====================================================
; Define cylindrical coordinate arrays for plotting 
;=====================================================

; Calculate total number of time intervals analysed
  num_int = sum( (tf0 - ts0) + 1)
  print("Total number of times in sample: "+num_int)

; Create new array to hold storm centre positions [7 * 59 * num_int * 2]
  num_l      = 59 
  centre_new = new( (/7, num_l, num_int, 2/), float)

  pii        = 3.14159265
  radii      = fspan(0.0,ar,nr)      ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,nt)   ; Array of azimuth angles (0, 5, 10, etc)
  thetas_rad = (pii/180)*thetas      ; Degrees to radians
  rad_size   = dimsizes(radii)-1     ; Index for later calculation

; 3D arrays to hold 'xpos' and 'ypos' values for all analysis times 
  xpos_all  = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)
  ypos_all  = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)
  xcen_all  = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)
  ycen_all  = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)

; 3D arrays for all variables
  u_int0    = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)  ; Zonal wind         
  v_int0    = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)  ; Meridional wind    
  w_int0    = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)  ; Vertical velocity  
  str_int0  = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)  ; Strain rate        
  vort_int0 = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)  ; Relative vorticity

  spd_int0  = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)  ; Vector windspeed
  vt_int0   = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)  ; Tangential wind 

  rad_arr0  = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)  ; Radius (for v/r)
  rdr_arr0  = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)  ; Radial distance
  va_int0   = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)  ; Angular velocity

  u_int0@description   = "Zonal wind"
  u_int0@units         = "m s~S~-1~N~"
  u_int0!0             = "rad"
  u_int0!1             = "azi"
  u_int0!2             = "time"

  v_int0@description   = "Meridional wind"
  v_int0@units         = "m s~S~-1~N~"
  v_int0!0             = "rad"
  v_int0!1             = "azi"
  v_int0!2             = "time"

  w_int0@description   = "Vertical velocity"
  w_int0@units         = "m s~S~-1~N~"
  w_int0!0             = "rad"
  w_int0!1             = "azi"
  w_int0!2             = "time"

  str_int0@description = "Strain rate"
  str_int0@units       = "s~S~-2~N~"
  str_int0!0           = "rad"
  str_int0!1           = "azi"
  str_int0!2           = "time"

  vort_int0@description = "Relative vorticity"
  vort_int0@units       = "s~S~-1~N~"
  vort_int0!0           = "rad"
  vort_int0!1           = "azi"
  vort_int0!2           = "time"

  spd_int0@description  = "Vector windspeed"
  spd_int0@units        = "m s~S~-1~N~"
  spd_int0!0            = "rad"
  spd_int0!1            = "azi"
  spd_int0!2            = "time"

  vt_int0@description   = "Tangential wind"
  vt_int0@units         = "m s~S~-1~N~"
  vt_int0!0             = "rad"
  vt_int0!1             = "azi"
  vt_int0!2             = "time"

  va_int0@description   = "Angular velocity"
  va_int0@units         = "s~S~-1~N~"
  va_int0!0             = "rad"
  va_int0!1             = "azi"
  va_int0!2             = "time"

; Initialise counter variable 
  ct                    = 0

; Print out relevant information about storm tracking
  sc_arr      = (/"relative vorticity",\                          ; sc = 0    
                  "pressure",\                                    ; sc = 1   
                  "smoothed pressure",\                           ; sc = 2
                  "extra smoothed pressure",\                     ; sc = 3  
                  "vort_min during symmetric phase",\             ; sc = 4
                  "windspeed minimum",\                           ; sc = 5
                  "geopotential"/)                                ; sc = 6

  print_clock("Using "+sc_arr(sc)+" to calculate storm centre position...")

;==========================================================
; Now start loop over these chosen input files and times
;==========================================================

 do st = 0, dimsizes(dat)-1 

  dt0 = tf0(st) - ts0(st)

; Find 'pc' and 'pd' files (model height level data)
  diri        = "/nfs/a319/earshar/"+dat(st)+"/"+ens0(st)
  fili_prefix = "201607"+dat(st)+"00Z_NPTK_4p4_L80_ra1t_"+ens0(st)
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

;============================================================
; Calculate storm motion using built-in function ('st_rm')
;============================================================

; 'dat'      = initialisation time (02T12, 03T00, ...)
; 'ens0'     = ensemble simulation (em00, em01, ...)
; 'diri'     = path to input files (see above)
; 'dist'     = size of box following storm (degrees)
; 'mins'     = analyse 1-h (0) or 5-min (1) data

  storm_rel  = st_rm(dat(st), ens0(st), "$ar/text/", dist, 0)
  u_cyc      = storm_rel[0]		   ; Zonal wind 
  v_cyc      = storm_rel[1]		   ; Meridional wind
  vel_cyc    = storm_rel[2]		   ; Vector wind 
  lat_arr    = storm_rel[4]                ; Storm-following domain subset (lat)
  lon_arr    = storm_rel[5]                ; Storm-following domain subset (lon)
  centre     = storm_rel[6]		   ; Storm centre at each hour (lat/lon)

;==========================================================================
; Create date/time string arrays for all times in file (use 'pc' stream)
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
  do ct0 = 0, numINPUT_c-1
   if (minute(ct0).gt.30) then
     hour(ct0) = hour(ct0)+1
   end if
  end do
  delete(ct0)

; Create arrays to hold strings for output plots
  date_str  = new(numINPUT_c,string)
  time_str  = new(numINPUT_c,string)
  time_arr  = new(numINPUT_c,string)
  title_arr = new(numINPUT_c,string)

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

   print("Starting main part of script...")
   print_clock("Working on ensemble member: "+ens0(st)+" ("+dat(st)+")" )

   llbox       = toint(dist*50)	           ; Size of lat/lon arrays (see line below)
   dsize       = (/llbox,llbox/)           ; Size of lat/lon arrays (depends on 'dist')
   d0	       = dsize(0)
   d1	       = dsize(1)

 ; Array to hold values of inner-core parameter
   sym_arr     = new(numINPUT_d,float)

   setvalues NhlGetWorkspaceObjectId
     "wsMaximumSize" : 1000000000
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

;===============================
; Start loop over input files 
;===============================

   c  = addfile(fili_c,"r")         	 ; Read in 'pc' stream [u,v,vort]
   d  = addfile(fili_d,"r")         	 ; Read in 'pd' stream [w,p,t,theta,kmh,tau]

;==================================
; Get the variables we will need
;==================================

   do it = ts0(st), tf0(st)		 ; Loop over times in file

    print("Working on time: "+time_str(it)+" (ct = "+ct+")" )
    itt = it
    title_arr(it) = "Valid at "+time_str(it)+" (T+"+itt+")"

;=====================================================
; Read in variables from 'pc' stream (model levels)
;=====================================================

    t0    = lat_arr(cn,it-1,0)
    t1    = lat_arr(cn,it-1,1)
    n0    = lon_arr(cn,it-1,0)
    n1    = lon_arr(cn,it-1,1)
    print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

  ; Combine values above into array for input into external function 
    ll_arr = (/t0,t1,n0,n1/)

  ; Read in basic variables from 'pc' stream
    lon1C  = c->longitude_1({n0:n1})              ; longitude ['d1' grid points]
    lat1C  = c->latitude_1({t0:t1})               ; latitude  ['d0' grid points]
    lonC   = c->longitude({n0:n1})                ; longitude ['d1' grid points]
    latC   = c->latitude({t0:t1})                 ; latitude  ['d0' grid points]
    hybC   = c->hybrid_ht(:) 		          ; 63 model ('theta') levels

  ; Calculate model level height information
    hy0    = sprintf("%0.0f",hybC(mlev0))
    hy1    = sprintf("%0.0f",hybC(mlev1))

  ; Read in horizontal wind components 
    u_varname = "x-wind"
    v_varname = "y-wind"

  ; Either read in multiple values over a layer or on a single level
    if (lay .eq. 1) then
     u0 = c->$u_varname$(it,mlev0:mlev1,{t0:t1},{n0:n1})
     v0 = c->$v_varname$(it,mlev0:mlev1,{t0:t1},{n0:n1})
    else
     u0 = c->$u_varname$(it,mlev0,{t0:t1},{n0:n1})
     v0 = c->$v_varname$(it,mlev0,{t0:t1},{n0:n1})
    end if

;===================================
; Calculate storm-relative winds
;===================================

    u = u0 - u_cyc(cn,it-1)
    v = v0 - v_cyc(cn,it-1)

  ; Add metadata from original horizontal wind arrays
    copy_VarCoords(u0,u)
    copy_VarAtts(u0,u)

    copy_VarCoords(v0,v)
    copy_VarAtts(v0,v)

  ; Calculate relative vorticity using centered finite differences
  ; Option '2' --> boundary points estimated using one-sided difference scheme
    vort = uv2vr_cfd(u, v, latC, lonC, 2)
    copy_VarCoords(u, vort)
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

;===========================================================
; Calculate rate of strain (following Nguyen et al. 2011)
;===========================================================

    dim1   = "latitude"
    dim2   = "longitude"

  ; Calculate horizontal derivatives of 'u' (du_dx, du_dy)
    du     = grad_latlon_cfd(u,u&$dim1$,u&$dim2$,False,False)
    du_dy  = du[0]    ; Latitudinal derivative
    du_dx  = du[1]    ; Longitudinal derivative
    delete(du)

  ; Calculate horizontal derivatives of 'v' (dv_dx, dv_dy)
    dv     = grad_latlon_cfd(v,v&$dim1$,v&$dim2$,False,False)
    dv_dy  = dv[0]    ; Latitudinal derivative
    dv_dx  = dv[1]    ; Longitudinal derivative
    delete(dv)

  ; Calculate strain --> sqrt[ (du_dx - dv_dy)^2 + (dv_dx + du_dy)^2 ]
    str1a   = du_dx - dv_dy   ; Stretching deformation
    str2a   = dv_dx + du_dy   ; Shearing deformation
    str1    = (str1a)^2
    str2    = (str2a)^2

    if (lay .eq. 1) then
     str0   = sqrt( str1(:,:,:) + str2(:,:,:) )
    else
     str0   = sqrt( str1(:,:) + str2(:,:) )
    end if

;======================================================
; Read in variables from 'pd' stream (model levels)
;======================================================

  ; Read in basic variables from 'pd' stream 
    lon1D  = d->longitude_1({n0:n1})                ; longitude [xxx grid points]
    lat1D  = d->latitude_1({t0:t1})                 ; latitude  [xxx grid points]
    lonD   = d->longitude({n0:n1})                  ; longitude [xxx grid points]
    latD   = d->latitude({t0:t1})                   ; latitude  [xxx grid points]
    hybD   = d->hybrid_ht(:)			    ; 63 model ('rho') levels

  ; Read in pressure and vertical velocity
    if (lay .eq. 1) then 
     p     = d->p(it-2,mlev0:mlev1,{t0:t1},{n0:n1})
     w     = d->dz_dt(it-2,mlev0:mlev1,{t0:t1},{n0:n1})
     th    = d->theta(it-2,mlev0:mlev1,{t0:t1},{n0:n1})
    else
     p     = d->p(it-2,mlev0,{t0:t1},{n0:n1})
     w     = d->dz_dt(it-2,mlev0,{t0:t1},{n0:n1})
     th    = d->theta(it-2,mlev0,{t0:t1},{n0:n1})
    end if 

    dx0    = lonD(1) - lonD(0) ; Grid spacing (longitude)
    dy0    = latD(1) - latD(0) ; Grid spacing (latitude)

;==========================
; Loop over model levels 
;==========================

    hy0 = sprintf("%0.0f",hybD(mlev0))
    hy1 = sprintf("%0.0f",hybD(mlev1))

  ; Read in variables on model levels
    if (lay .eq. 1) then

      print("Working on layer between " + hy0+ " and "+hy1+" m AGL")

    ; pc stream
      u_plane    = dim_avg_n_Wrap( u(:,:,:),0)
      v_plane    = dim_avg_n_Wrap( v(:,:,:),0)
      vort_plane = dim_avg_n_Wrap( vort(:,:,:),0)
      avo_plane  = dim_avg_n_Wrap( avo(:,:,:),0)
      str_plane  = dim_avg_n_Wrap( str0(:,:,:),0)
      copy_VarCoords(u(0,:,:),str_plane)

    ; pd stream
      w_plane    = dim_avg_n_Wrap( w(:,:,:),0)
      prs_plane  = dim_avg_n_Wrap( p(:,:,:), 0)
      th_plane   = dim_avg_n_Wrap( th(:,:,:),0)
      
    else

      print("Working on " + hybD(mlev0)+ " m AGL")

    ; pc stream
      u_plane    = u(:,:)
      v_plane    = v(:,:)
      vort_plane = vort(:,:)
      avo_plane  = avo(:,:)
      str_plane  = str0(:,:)
      copy_VarCoords(u,str_plane)

    ; pd stream
      w_plane    = w(:,:)
      prs_plane  = w(:,:)
      th_plane   = w(:,:)

    end if

    ; Calculate vector windspeed and smooth
      spd_plane     = sqrt( (u_plane ^ 2) + (v_plane ^ 2) )
      copy_VarCoords(u_plane, spd_plane)
      vort_pl_smth  = smth9_Wrap(vort_plane, 0.5, 0.5, True)

    ; Combine all pc and pd stream variables into larger arrays
    ; 28/10/2019 --> tidy up this function so that we don't have empty options
      pc_plane      = (/u_plane, v_plane, spd_plane, vort_plane, \
                        avo_plane, str_plane, str_plane, str_plane/)

    ; Second 'prs_plane' is a placeholder
      pd_plane      = (/w_plane, prs_plane, prs_plane, prs_plane/)

;=============================================================  
; Call external function to calculate storm centre position
;============================================================= 

    ; Set radial distance (º) when looking for vort min in ringlike phase   
      r0            = 0.10 

    ; Composite ("comp") or single ("single") simulation
    ; EDIT -- check that we don't need to redo our composite plot calculations 
      centre_arr    = st_centre(0, it, cn, r0, centre, vort_plane, prs_plane, spd_plane, \
                                centre_new, dy0, dx0, "comp", ct, "pd", "off", "1h")

    ; Updated storm track information contained in 'centre_new'
      centre_new    = centre_arr[0]

;=====================================================================  
; Call external function to switch to cylindrical coordinate system  
;=====================================================================  

    ; Set 'ilev' option before function call below    
      if (lay .eq. 1) then
       ilev = 0
      end if

    ; Call 'setup_cyl' and output pressure/vorticity on cylindrical grid 
    ; Second to last option = 0 --> basic variables output    
    ; Second to last option = 1 --> additional MBGT variables output  
    ; Last option --> composite ("comp") or single simulation ("single"); same as 'st_centre' call

      print("mlev0 = "+mlev0)

    ; EDIT -- editing 'setup_cyl' to define xpos_all, ypos_all, etc outside function
    ; Option after 'nr' MUST be called 'ilev', or script will produce missing values 
      cyl_arr  = setup_cyl(ct, dist, ll_arr, lonC, latC, lonD, latD, lon1D, lat1D, \
                           pc_plane, pd_plane, num_int, centre_new, \
                           sc, sm, ar, nr, nt, ilev, num_l, 0, "comp", ct, 1, \
                           xpos_all, ypos_all, xcen_all, ycen_all)

      xpos_all = cyl_arr[0]
      ypos_all = cyl_arr[1]
      xcen_all = cyl_arr[2]
      ycen_all = cyl_arr[3] 

      lat_max  = cyl_arr[4]
      lon_max  = cyl_arr[5]
      radii    = cyl_arr[6]
      rad_size = cyl_arr[7]
      thetas   = cyl_arr[19] 

      u_int    = cyl_arr[8]
      v_int    = cyl_arr[9]
      spd_int  = cyl_arr[10]
      vort_int = cyl_arr[11]
      vt_int   = cyl_arr[14]
      w_int    = cyl_arr[15]
      str_int  = cyl_arr[22]

;===================================================================================
; Follow method in 'tc_ring_mono.ncl' (find doughnut centre during ringlike phase)
;===================================================================================

 ; Define new arrays to hold azimuthally-averaged variables
   vort_azi = new( (/dimsizes(radii)/), float)
   spd_azi  = new( (/dimsizes(radii)/), float)

 ; Now calculate the azimuthal average of relative vorticity
   do irad  = 0, dimsizes(radii)-1
    vort_azi(irad) = avg( vort_int(irad,:) )
    spd_azi(irad)  = avg( spd_int(irad,:) )
   end do

 ; Use more sophisticated method to calculate whether vortex is ringlike or monopolar 
 ; This code will go into an NCL function 

 ; (1) Find radius (index) corresponding to relative vorticity maximum 
   ind_vmax    = ind_resolve(maxind(vort_azi), dimsizes(vort_azi) )

 ; Calculate vorticity at ind_vmax ('vort_max'), and at the vortex centre ('vort_centre')
   vort_max    = vort_azi(ind_vmax(0,0) )
   vort_centre = vort_azi(0)

 ; (2) Calculate the ratio of 'ind_max' to 'ind_centre' 
   sym_arr(it) = vort_centre / vort_max

 ; If 'ind_vmax' = 'ind_centre' (i.e. monopole), then 'sym_arr' will be equal to 1 
   if (sym_arr(it) .gt. 0.95) then

    print("Inner-core structure is monopolar (ratio = "+sym_arr(it)+")")
    sym_str = "Monopole: "+sprintf("%0.2f", sym_arr(it) )

   elseif (sym_arr(it) .lt. 0.80) then

    print("Inner-core structure is ringlike (ratio = "+sym_arr(it)+")")
    sym_str = "Ringlike: "+sprintf("%0.2f", sym_arr(it) )

   else
    print("Inner-core structure is not clearly monopolar or ringlike (ratio = "+sym_arr(it)+")")
    sym_str = "Not clear: "+sprintf("%0.2f", sym_arr(it) )

   end if

 ; If we're in the symmetric phase, put the storm centre in the middle of the vorticity ring 
   if (dnt .eq. 1 .and. typ .eq. "sym") then

    print("Re-calculating storm centre to fit ringlike vorticity profile...")

  ; Select the 'vort_min' method of storm tracking
    sc0 = 4

  ; Call the external function again 
    delete(cyl_arr)
    cyl_arr  = setup_cyl(ct, dist, ll_arr, lonC, latC, lonD, latD, lon1D, lat1D, \
                         pc_plane, pd_plane, num_int, centre_new, \
                         sc0, sm, ar, nr, nt, mlev0, num_l, 0, "comp", ct, 1, \
                         xpos_all, ypos_all, xcen_all, ycen_all)

    xpos_all = cyl_arr[0]
    ypos_all = cyl_arr[1]
    xcen_all = cyl_arr[2]
    ycen_all = cyl_arr[3]

    lat_max  = cyl_arr[4]
    lon_max  = cyl_arr[5]
    radii    = cyl_arr[6]
    rad_size = cyl_arr[7]
    thetas   = cyl_arr[19]

    u_int    = cyl_arr[8]
    v_int    = cyl_arr[9]
    spd_int  = cyl_arr[10]
    vort_int = cyl_arr[11]
    vt_int   = cyl_arr[14]
    w_int    = cyl_arr[15]
    prs_int  = cyl_arr[16]
    str_int  = cyl_arr[22]

  ; Recalculate azimuthal average of relative vorticity and windspeed 
    do irad = 0, dimsizes(radii)-1
     vort_azi(irad) = avg( vort_int(irad,:) )
     spd_azi(irad)  = avg( spd_int(irad,:) )
    end do

   end if

;=================================
; Continue with original script
;=================================

  ; Fill larger arrays with calculated values (vvel; vort; speed; vtan; strain)
    w_int0(:,:,ct)    = (/w_int(:,:)/)
    vort_int0(:,:,ct) = (/vort_int(:,:)/)
    spd_int0(:,:,ct)  = (/spd_int(:,:)/)
    vt_int0(:,:,ct)   = (/vt_int(:,:)/)
    str_int0(:,:,ct)  = (/str_int(:,:)/)

    if (plt .eq. 1) then 

    ; Plot the data as we go 
      output_loop = "$sam/nepartak/images/cyl_coords/xy_comp_"+typ+"_group"+grp+"_sc"+sc+"_t"+ct
      wks = gsn_open_wks(opt2,output_loop)

      gsn_define_colormap(wks,"prcp_new")

    ; Relative vorticity resources 
      opts_rv                               = True
      opts_rv@cnFillOn                      = True
      opts_rv@sfXArray                      = xpos(ilev,:,:,ct)
      opts_rv@sfYArray                      = ypos(ilev,:,:,ct)
      opts_rv@cnLineLabelInterval           = 2.0
      opts_rv@cnLineLabelFontHeightF        = 0.012
      opts_rv@cnLineLabelBackgroundColor    = "transparent"
      opts_rv@cnLineLabelPlacementMode      = "constant"
      opts_rv@cnLinesOn                     = False
      opts_rv@cnInfoLabelOn                 = False
      opts_rv@cnLevelSelectionMode          = "ExplicitLevels"
      opts_rv@cnLevels                      = (/5., 10., 15., 20., 25., 30., \
                                                40., 50., 60., 70., 80./)
      opts_rv@cnFillColors                  = (/0,4,5,6,7,8,9,10,11,12,13,14/)
      opts_rv@gsnMaximize                   = False      

      opts_rv@tiMainString                  = "Simulation "+ens0(st)+" ("+dat(st)+"): T+"+it
      opts_rv@gsnLeftString                 = ""
      opts_rv@gsnRightString                = ""
      opts_rv@lbLabelBarOn                  = False

      opts_rv@trYMaxF                       = ycen + rad0
      opts_rv@trYMinF                       = ycen - rad0
      opts_rv@trXMaxF                       = xcen + rad0
      opts_rv@trXMinF                       = xcen - rad0
      opts_rv@gsnDraw                       = False
      opts_rv@gsnFrame                      = False
      opts_rv@gsnAddCyclic                  = False

    ; Background plot resources 
      opts_xy                               = True

    ; No tickmarks on any axis 
      opts_xy@tmYROn                        = False
      opts_xy@tmYLOn                        = False
      opts_xy@tmXTOn                        = False
      opts_xy@tmXBOn                        = False

      opts_xy@xyLineColor                   = -1
      opts_xy@gsnDraw                       = False
      opts_xy@gsnFrame                      = False

      opts_xy@trYMaxF                       = ycen + rad0
      opts_xy@trYMinF                       = ycen - rad0
      opts_xy@trXMaxF                       = xcen + rad0
      opts_xy@trXMinF                       = xcen - rad0
      opts_xy@gsnMaximize                   = False

    ; Panel position
      xf = 0.12
      yf = 0.90
      wf = 0.80
      hf = 0.75

    ; Produce panel plot 
      opts_rv@vpXF      = xf
      opts_rv@vpYF      = yf
      opts_rv@vpWidthF  = wf
      opts_rv@vpHeightF = hf

      opts_xy@vpXF      = xf
      opts_xy@vpYF      = yf
      opts_xy@vpWidthF  = wf
      opts_xy@vpHeightF = hf

    ; Changed units for plotting 
      vort_int(:,:)     = vort_int(:,:) * (10 ^ 4)

      plot_bck          = gsn_csm_xy(wks, xpos(ilev,:,:,ct), ypos(ilev,:,:,ct), opts_xy)
      plot_cont         = gsn_csm_contour(wks, vort_int(:,:), opts_rv)
      overlay(plot_bck, plot_cont)

;==============================
; Add radial circles to plot 
;==============================

      dtr = 0.017453292519943

    ; Get customizations for circles
      radius   = radii(rad_size)
      xcenter  = lon_max
      ycenter  = lat_max
      out_thck = 2
      in_thck  = 1
      spacing  = 0.05
      dpattrn  = 2
      out_clr  = 1
      in_clr   = 1

;=============================================
; Construct and attach outer circle to plot 
;=============================================

    ; Calculate arrays for outer circle
      degrees  = ispan(0, 360, 5)
      xcos     = cos(dtr * degrees)
      xsin     = sin(dtr * degrees)      

    ; Array of points (x,y) representing the outer circle
      xc       = xcenter + (radius * xcos)
      yc       = ycenter + (radius * xsin)

    ; Resources for outer circle 
      lnres    = True
      lnres@gsLineThicknessF  = out_thck
      lnres@gsLineColor       = out_clr
      lnres@gsLineDashPattern = 0      

    ; Attach outer circle to plot 
      plot_outer = gsn_add_polyline(wks, plot_bck, xc, yc, lnres)

      if (spacing .gt. radius) then
       print("add_radial_circles: spacing is > radius, can't draw inner circles.")
       return
      end if    

;==================================================
; Now construct and attach inner circles to plot
;==================================================

    ; Draw inner circles if desired
      if (spacing .gt. 0) then

       count = 0
       size  = toint( (radius-spacing) / spacing )
       plot_inner = new(size,graphic)

       do r = spacing, radius-spacing, spacing

       ; Calculate arrays for inner circles 
         xc = xcenter + (r * xcos)
         yc = ycenter + (r * xsin)

       ; Resources for inner circles  
         delete([/lnres@gsLineThicknessF,lnres@gsLineColor,lnres@gsLineDashPattern/])
         lnres@gsLineColor       = in_clr
         lnres@gsLineThicknessF  = in_thck
         lnres@gsLineDashPattern = dpattrn

       ; Attach inner circles to plot
         plot_inner(count) = gsn_add_polyline(wks, plot_bck, xc, yc, lnres)

         count = count + 1

       end do

      end if

    ; Tidy up
      delete([/lnres, degrees, xcos, xsin, xc, yc/])
      delete([/opts_rv, opts_xy, plot_bck, plot_cont, xf, yf, wf, hf, plot_inner/])

    ; Call 'draw' and 'frame'
      maximize_output(wks,True)

  ; End IF statement for plotting as we go
    end if 

    ; Tidy up
      delete([/w_int,vort_int,u_int,v_int,spd_int,vt_int/])
      delete([/w_plane, vort_plane, avo_plane, u_plane, v_plane, spd_plane, str_plane/])

     ct  = ct + 1		    ; Counter variable (time)

  ; Tidy up before next iteration (all variables from 'pc' and 'pd' streams)
    delete([/lonC,latC,lon1C,lat1C,u,v,vort/])
    delete([/lonD,latD,lon1D,lat1D,w/])

   end do     ; End time loop (do it = 0, times-1)

  end do      ; End loop over input file times (do st = 0, dimsizes(dat)-1)

;======================================================
; Calculate radius of maximum wind (tangential wind)
;======================================================

; CODE GOES HERE (see tc_write_rmw_vtan.ncl, ~L540)

;===================================================
; Calculate angular velocity from tangential wind 
;===================================================

  drad = dimsizes(radii)

; Calculate radial distance from storm centre at all grid points 
; rad_arr0 dims: [dimsizes(radii), dimsizes(thetas), num_int]
  do st = 0, num_int-1
   rad_arr0(1:drad-1,:,st)  = gc_latlon( ycen_all(1:drad-1,:,st),xcen_all(1:drad-1,:,st),\
                                         ypos_all(1:drad-1,:,st),xpos_all(1:drad-1,:,st), 2, 3)
  end do 

; Replace zeros with small, positive values (avoid divide-by-zero errors)
  rad_arr0(0,:,:) = 1

; Also calculate radial distance between all grid points 
  do st = 0, num_int-1 
   rdr_arr0(0:drad-2,:,st)  = gc_latlon( ypos_all(0:drad-2,:,st),xpos_all(0:drad-2,:,st),\
                                         ypos_all(1:drad-1,:,st),xpos_all(1:drad-1,:,st), 2, 3)
  end do 

; Calculate angular velocity and change units 
  va_int0   = vt_int0(:,:,:) / rad_arr0(:,:,:)
  va_int0   = va_int0 * (10 ^ 4)
  va_int0@units   = "10~S~-4~N~ s~S~-1~N~"

; Also change units of relative vorticity 
  vort_int0 = vort_int0 * (10 ^ 4)
  vort_int0@units = "10~S~-4~N~ s~S~-1~N~"

;=================================================
; Produce time-averaged arrays of key variables
;=================================================

 ; Calculate azimuthal average of angular velocity 
   va_avg0       = dim_avg_n_Wrap(va_int0(:,:,:), 1)
   va_avg0!0     = "rad"
   va_avg0&rad   = radii
   va_avg0!1     = "time"

 ; Tangential wind 
   vt_avg0       = dim_avg_n_Wrap(vt_int0(:,:,:), 1)
   vt_avg0!0     = "rad"
   vt_avg0&rad   = radii
   vt_avg0!1     = "time"

 ; Relative vorticity 
   vort_avg0     = dim_avg_n_Wrap(vort_int0(:,:,:), 1)
   vort_avg0!0   = "rad"
   vort_avg0&rad = radii
   vort_avg0!1   = "time"

;===============================================
; Produce radial profiles of angular velocity 
;=============================================== 

 ; Do I need to create new object to hold values? 
 ; 'va_plot' should hold each radial profile (num_int = 48)
   va_plot   = new( (/num_int, dimsizes(radii)/), float)
   vt_plot   = new( (/num_int, dimsizes(radii)/), float)
   vort_plot = new( (/num_int, dimsizes(radii)/), float)

   do ct = 0, num_int-1
    va_plot(ct,:)   = va_avg0(:,ct)
    vt_plot(ct,:)   = vt_avg0(:,ct)
    vort_plot(ct,:) = vort_avg0(:,ct) 
   end do 

 ; Angular velocity plotting options 
   va_opts                      = True 
   va_opts@gsnDraw              = False
   va_opts@gsnFrame             = False

 ; Y-axis
   if (var .eq. "va") then 
    va_opts@tiMainString        = "" ; Radial profile of angular velocity ("+out_str+")"
    va_opts@tiYAxisString       = "Angular velocity (10~S~-4~N~ s~S~-1~N~)"
    va_opts@trYMaxF             = 50.0
    va_opts@trYMinF             = 0.0
   elseif (var .eq. "vt") then 
    va_opts@tiMainString        = "" ; Radial profile of tangential wind ("+out_str+")"
    va_opts@tiYAxisString       = "Tangential wind (m s~S~-1~N~)"
    va_opts@trYMaxF             = 65.0
    va_opts@trYMinF             = 0.0
   else
    va_opts@tiMainString        = "" ; Radial profile of relative vorticity ("+out_str+")"
    va_opts@tiYAxisString       = "Relative vorticity (10~S~-4~N~ s~S~-1~N~)"
    va_opts@trYMaxF             = 120.0
    va_opts@trYMinF             = 0.0
   end if 
   va_opts@tiYAxisFontHeightF   = 0.0150

 ; X-axis (define our own labels)
   va_opts@tiXAxisString        = "Radius (km)"
   va_opts@tiXAxisFontHeightF   = 0.0150
   va_opts@tmXBLabelFontHeightF = 0.0125
   va_opts@tmXBMode             = "Explicit"
   va_opts@tmXBValues           = (/0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45/)
   va_opts@tmXBLabels           = (/"0", "5", "10", "15", "20", "25", "30", "35", "40", "45"/)

 ; Additional resources 
   va_opts@xyLineColor          = "black"
   va_opts@xyLineThicknessF     = 2.0
   va_opts@xyDashPatterns       = 0.0
   va_opts@xyMarkLineMode       = "Lines"

 ; Turn off tickmarks
   va_opts@tmXTOn               = "False"
   va_opts@tmYROn               = "False"

 ; Additional resources (x-axis, etc)
   va_opts@trXMaxF              = radii(dimsizes(radii)-2)
   va_opts@trXMinF              = 0.00

 ; Open workstation and define output file path
   output = "$sam/nepartak/images/vc/"+var+"_"+typ+"_group"+grp+"_sc"+sc 
   wks    = gsn_open_wks(opt, output) 

 ; Finally, produce radial plot of chosen variable
   if (var .eq. "va") then 
    plot  = gsn_csm_xy(wks, radii, va_plot(:,:), va_opts)
   elseif (var .eq. "vt") then 
    plot  = gsn_csm_xy(wks, radii, vt_plot(:,:), va_opts)
   else
    plot  = gsn_csm_xy(wks, radii, vort_plot(:,:), va_opts)
   end if 

   draw(plot)
   frame(wks)

   exit()

;=================================================
; Now average all profiles and produce new plot
;=================================================

 ; Now average along time dimension 
   va_arr    = dim_avg_n_Wrap(va_avg0(:,:), 1)

 ; Add metadata before plotting 
   va_arr&rad     = radii

 ; Write values to text file 
   va_out = "$ar/text/"+var+"_"+typ+"_sc"+sc+".txt"
   asciiwrite(va_out, va_arr)

   exit()

 ; Output file location and type
   output = "$sam/nepartak/images/cyl_coords/xy_comp_"+typ+"_group"+grp+"_sc"+sc
   wks = gsn_open_wks(opt,output)
   
 ; Load colour table 
   gsn_define_colormap(wks,"prcp_new")

;==========================
; Options for plotting 
;==========================

 ; Relative vorticity                                                             
   opts_rv                               = True
   opts_rv@cnFillOn                      = True
   opts_rv@sfXArray                      = xpos_all(ilev,:,:,ct)
   opts_rv@sfYArray                      = ypos_all(ilev,:,:,ct)
   opts_rv@cnLineLabelInterval           = 2.0
   opts_rv@cnLineLabelFontHeightF        = 0.012
   opts_rv@cnLineLabelBackgroundColor    = "transparent"
   opts_rv@cnLineLabelPlacementMode      = "constant"

 ; Turn contour lines and labels off
   opts_rv@cnLinesOn                     = False
   opts_rv@cnInfoLabelOn                 = False
   opts_rv@cnLevelSelectionMode          = "ExplicitLevels"
   opts_rv@cnLevels                      = (/5., 10., 15., 20., 25., 30., \
                                             40., 50., 60., 70., 80./)
   opts_rv@cnFillColors                  = (/0,4,5,6,7,8,9,10,11,12,13,14/)
   opts_rv@gsnMaximize                   = False

 ; No title or additional strings 
   opts_rv@tiMainString                  = ""
   opts_rv@gsnLeftString                 = ""
   opts_rv@gsnRightString                = ""

 ; Turn individual labelbar on/off
   opts_rv@lbLabelBarOn                  = False

 ; Labelbar resources (only active if lb switched on)
   opts_rv@pmLabelBarWidthF              = 0.52
   opts_rv@pmLabelBarHeightF             = 0.13
   opts_rv@lbLabelFontHeightF            = 0.0125
   opts_rv@lbPerimOn                     = False
   opts_rv@lbLabelFont                   = "Helvetica"
   opts_rv@lbBoxEndCapStyle              = "TriangleBothEnds"
   opts_rv@lbTitleFontHeightF            = 0.02
   opts_rv@lbTitleDirection              = "Across"
   opts_rv@lbTitlePosition               = "Bottom"

   opts_rv@trYMaxF                       = ycen + rad0
   opts_rv@trYMinF                       = ycen - rad0
   opts_rv@trXMaxF                       = xcen + rad0
   opts_rv@trXMinF                       = xcen - rad0

   if (lay .eq. 0) then

    opts_rv@lbTitleString                = hy0+\
                                              " m relative vorticity (10~S~-4~N~ s~S~-1~N~)"

   else

    opts_rv@lbTitleString                = "Layer-averaged ("+hy0+"-"+hy1+\
                                           " hPa) ~C~     relative vorticity (10~S~-4~N~ s~S~-1~N~)"

   end if

 ; Do not draw the plot or advance the frame 
   opts_rv@gsnDraw                       = False
   opts_rv@gsnFrame                      = False
   opts_rv@gsnAddCyclic                  = False

 ; Resources for background plot (cylindrical grid)
 ; EDIT -- make these lines transparent 
   opts_xy                               = True 
   opts_xy@tmYROn                        = "False"
   opts_xy@tmXTOn                        = "False"
   opts_xy@xyLineColor                   = -1
   opts_xy@gsnDraw                       = False
   opts_xy@gsnFrame                      = False

   opts_xy@trYMaxF                       = ycen + rad0
   opts_xy@trYMinF                       = ycen - rad0
   opts_xy@trXMaxF                       = xcen + rad0
   opts_xy@trXMinF                       = xcen - rad0
   opts_xy@gsnMaximize                   = False

 ; Resources for windspeed
   opts_wind                             = opts_rv
   delete([/opts_wind@cnLevels, opts_wind@cnFillColors/])
   opts_wind@cnFillOn                    = True
   opts_wind@cnLinesOn                   = False 
   opts_wind@cnLineLabelsOn              = False 
   opts_wind@cnInfoLabelOn               = False  
   opts_wind@cnLevelSelectionMode        = "ExplicitLevels" 
   opts_wind@cnFillPalette               = "radar_new"
   opts_wind@cnLevels                    = (/5.0, 10.0, 15.0, 20.0,\
                                             25.0, 30.0, 35.0, 40.0, \
                                             45.0, 50.0, 55.0, 60.0, 65.0/)
   opts_wind@cnFillColors                = (/-1,0,1,2,3,4,5,6,7,8,9,10,11,12/)
   opts_wind@lbLabelBarOn                = False  
   opts_wind@gsnMaximize                 = False

 ; Resources for vertical velocity 
   opts_vvel                             = opts_wind
   delete([/opts_vvel@cnLevels, opts_vvel@cnFillColors, opts_vvel@cnFillPalette/])
   opts_vvel@cnFillPalette               = "n11"
   opts_vvel@cnLevels                    = (/0.05, 0.1, 0.2, 0.3, 0.4, \
                                             0.5, 0.6, 0.7, 0.8, 1.0/)
   opts_vvel@cnFillColors                = (/0,2,3,4,5,6,7,8,9,10,11/)

 ; Resources for strain rate 
   opts_strn                              = opts_vvel
   delete([/opts_strn@cnLevels, opts_strn@cnFillColors, opts_strn@cnFillPalette/])
   opts_strn@cnFillPalette                = "circ_budget"
   opts_strn@cnLevels                     = (/10,12,14,16,18,20,22,25,30/)
   opts_strn@cnFillColors                 = (/-1,11,13,15,16,17,18,19,20,21/)

 ; Panel plot label resources 
   txid                                   = new(4,graphic)
   amid                                   = new(4,graphic)

   panel_strings                          = (/ "a) Relative vorticity",\
                                               "b) Vertical velocity",\
                                               "c) Vector wind",\
                                               "d) Strain rate"/)

   txres                                  = True
   txres@txPerimOn                        = True
   txres@txFontHeightF                    = 0.0125
   txres@txBackgroundFillColor            = "White"

   amres                                  = True
   amres@amParallelPosF                   = -0.48
   amres@amOrthogonalPosF                 = -0.48
   amres@amJust                           = "TopLeft"

;=============                                                                              
; Plot data                                                                                 
;=============                                                                              

 ; Coordinates of panels
   xf = (/0.08, 0.58/)
   yf = (/0.70, 0.70/)
   wf = 0.40
   hf = 0.35

 ; Panel 1 (tangential wind)
   opts_rv@vpXF       = xf(0)
   opts_rv@vpYF       = yf(0)
   opts_rv@vpWidthF   = wf
   opts_rv@vpHeightF  = hf   

   opts_xy0           = opts_xy
   opts_xy0@vpXF      = xf(0)
   opts_xy0@vpYF      = yf(0)
   opts_xy0@vpWidthF  = wf
   opts_xy0@vpHeightF = hf

   plot_tl            = gsn_csm_xy(wks, xpos_all(ilev,:,:,ct), ypos_all(ilev,:,:,ct), opts_xy0)
   plot_tl0           = gsn_csm_contour(wks, vtan_arr, opts_rv)
   overlay(plot_tl, plot_tl0)

 ; Overlay labels 
   txid(0)      = gsn_create_text(wks, panel_strings(0), txres)
   amid(0)      = gsn_add_annotation(plot_tl, txid(0), amres)

 ; Panel 2 (angular velocity)
   opts_vvel@vpXF       = xf(1)
   opts_vvel@vpYF       = yf(1)
   opts_vvel@vpWidthF   = wf
   opts_vvel@vpHeightF  = hf

   opts_xy1             = opts_xy
   opts_xy1@vpXF        = xf(1)
   opts_xy1@vpYF        = yf(1)
   opts_xy1@vpWidthF    = wf
   opts_xy1@vpHeightF   = hf

   plot_tr              = gsn_csm_xy(wks, xpos_all(ilev,:,:,ct), ypos_all(ilev,:,:,ct), opts_xy1)
   plot_tr0             = gsn_csm_contour(wks, vang_arr, opts_vvel)
   overlay(plot_tr, plot_tr0)

 ; Overlay labels
   txid(1)      = gsn_create_text(wks, panel_strings(1), txres)
   amid(1)      = gsn_add_annotation(plot_tr, txid(1), amres)

;==============================
; Add radial circles to plot 
;==============================

   dtr = 0.017453292519943

 ; Get customizations for circles
   radius   = radii(rad_size)
   xcenter  = lon_max
   ycenter  = lat_max
   out_thck = 2
   in_thck  = 1
   spacing  = 0.05
   dpattrn  = 2
   out_clr  = 1
   in_clr   = 1

;============================================= 
; Construct and attach outer circle to plot
;============================================= 

 ; Calculate arrays for outer circle
   degrees  = ispan(0, 360, 5)
   xcos     = cos(dtr * degrees)
   xsin     = sin(dtr * degrees)

 ; Array of points (x,y) representing the outer circle
   xc       = xcenter + (radius * xcos)
   yc       = ycenter + (radius * xsin)

 ; Resources for outer circle 
   lnres    = True
   lnres@gsLineThicknessF  = out_thck
   lnres@gsLineColor       = out_clr
   lnres@gsLineDashPattern = 0

 ; Attach outer circle to each panel plot  
   plot_outer0 = gsn_add_polyline(wks, plot_tl, xc, yc, lnres)
   plot_outer1 = gsn_add_polyline(wks, plot_tr, xc, yc, lnres)
   plot_outer2 = gsn_add_polyline(wks, plot_bl, xc, yc, lnres)
   plot_outer3 = gsn_add_polyline(wks, plot_br, xc, yc, lnres)

   if (spacing .gt. radius) then
    print("add_radial_circles: spacing is > radius, can't draw inner circles.")
    return
   end if

;================================================== 
; Now construct and attach inner circles to plot 
;================================================== 

 ; Draw inner circles if desired 
   if (spacing .gt. 0) then

    count = 0
    size  = toint( (radius-spacing) / spacing )
    plot_inner = new( (/4,size/), graphic)

    do r = spacing, radius, spacing

    ; Calculate arrays for inner circles
      xc = xcenter + (r * xcos)
      yc = ycenter + (r * xsin)

    ; Resources for inner circles 
      delete([/lnres@gsLineThicknessF,lnres@gsLineColor,lnres@gsLineDashPattern/])
      lnres@gsLineColor       = in_clr
      lnres@gsLineThicknessF  = in_thck
      lnres@gsLineDashPattern = dpattrn 

    ; Attach inner circles to plot  
      plot_inner(0,count) = gsn_add_polyline(wks, plot_tl, xc, yc, lnres)
      plot_inner(1,count) = gsn_add_polyline(wks, plot_tr, xc, yc, lnres)
      plot_inner(2,count) = gsn_add_polyline(wks, plot_bl, xc, yc, lnres)
      plot_inner(3,count) = gsn_add_polyline(wks, plot_br, xc, yc, lnres)

      count = count + 1

    end do

   end if

 ; Tidy up 
   delete([/lnres, degrees, xcos, xsin, xc, yc/])

;====================================================                                       
; Finally, draw the plot with everything overlaid                                           
;====================================================                                       

 ; Create a main title at the top of the page 
   title_res                   = True
   title_str                   = out_str+": composite plots"
   title_res@txFontHeightF     = 0.02
   title_res@txFont            = "Helvetica"
   xpos_title                  = 0.50
   ypos_title                  = 0.94

   gsn_text_ndc(wks, title_str, xpos_title, ypos_title, title_res)

 ; Call 'draw' and 'frame'
   maximize_output(wks,True)

end 