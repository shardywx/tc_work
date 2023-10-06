; Produce a boxplot of the difference in vertical mass flux between two levels/layers 

; Run using:

; ncl dist=0.6 nt=73 nr=9 ar=0.40 sr=1 wt=0.5 cn0=\"geo_sm\" lay=1 mlev0=17 mlev1=46 int=2 sc=5
; sm=0 diff=1 grp=\"6\" trb=0 typ=\"ring\" w_check=1 w0=0.3 tc_mflux_write.ncl

; 'dist'  = size of box following storm (degrees)
; 'typ'	  = VC phase ("sym","asym","stoa","atos")
; 'full'  = beginning points (9); middle points (10); end points (11)
; 'ar'    = distance in degrees of outer radius (1.0, 2.0, etc)
; 'sr'    = radial circle number (1 = 5 km, 2 = 10 km, etc) to start flux calculations
; 'nr'    = number of radial circles between r = 0 and r = 'ar' (21, 41, 61, etc)
; 'nt'    = number of azimuth angles in cylindrical grid
; 'cn0'   = offline storm centre: "slp", "vort", "geo", "geo_sm"
; 'diff'  = write out upper and lower layers separately (0) or write out only the difference (1)
; 'azi'   = (2) calculate azimuthal average and integrate radially
;           (1) calculate azimuthal average alone
;           (0) leave each grid point value alone
; 'w_check' = recalculate storm centre when sampling the eyewall updraft (w > 0)
; 'w0'    = 
; 'int'   = no. of levels either side of 'mlev0'/'mlev1' to calculate layer average
; 'wt'    = vertical velocity threshold (e.g. 0.5 m s-1) when integrating radially 

; 500 m [10], 1 km [14], 1.5 km [17], 2 km [20], 3 km [24], 4 km [28], 5 km [31], 6 km [34]
; 7 km [37], 8 km [39], 9 km [42], 10 km [45], 11 km [46], 12 km [48], 13 km [50], 14 km [52]  
; 15 km [54], 16 km [56]

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

; Group 5 --> Group 9 
  vc_arr_all = (/ (/56,38,51,30,22/), (/55,35,49,27,23/), \
                  (/69,39,64,37,27/), (/52,33,43,26,19/) /)

; String containing info on width of radial averaging band  
  ar0  = ar * 100
  rf = sprintf("%0.0fkm", ar0)

  sr0  = sr * 5
  rs = sprintf("%0.0fkm", sr0)

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
   typ0   = "ring"
   gr_ind = 0
  elseif (typ .eq. "s2a") then
   typ0   = "r2m"
   gr_ind = 1
  elseif (typ .eq. "asym") then
   typ0   = "mono"
   gr_ind = 2
  elseif (typ .eq. "a2s") then
   typ0   = "m2r"
   gr_ind = 3
  elseif (typ .eq. "novc") then
   typ0 = "novc"
  end if

; Which set of simulations do we want to analyse?

; No VCs   
  if (typ .eq. "novc") then
    novc_size = 16
    dat  = asciiread("$ar/text/novc_dat.txt",(/novc_size/),"string")
    ens0 = asciiread("$ar/text/novc_sim.txt",(/novc_size/),"string")
    ts0  = asciiread("$ar/text/novc_ts.txt",(/novc_size/),"integer")
    tf0  = asciiread("$ar/text/novc_tf.txt",(/novc_size/),"integer")
; Wind speed or MSLP tendency 
  elseif (typ .eq. "v" .or. typ .eq. "slp") then
    dat  = asciiread("$ar/text/vc_scatter_sc"+sc+"_"+typ+"_"+thr+"_dat.txt",(/arr_size/),"string")
    ens0 = asciiread("$ar/text/vc_scatter_sc"+sc+"_"+typ+"_"+thr+"_ens.txt",(/arr_size/),"string")
    ts0  = asciiread("$ar/text/vc_scatter_sc"+sc+"_"+typ+"_"+thr+"_time.txt",(/arr_size/),"integer")
    tf0  = ts0
; Any other VC phase 
  else
   dat =asciiread("$ar/text/vc_"+typ+"_dat_group"+grp+".txt",(/vc_arr_all(gr_ind,vc_ind)/),"string")
   ens0=asciiread("$ar/text/vc_"+typ+"_sim_group"+grp+".txt",(/vc_arr_all(gr_ind,vc_ind)/),"string")
   ts0 =asciiread("$ar/text/vc_"+typ+"_ts_group"+grp+".txt",(/vc_arr_all(gr_ind,vc_ind)/),"integer")
   tf0 =asciiread("$ar/text/vc_"+typ+"_tf_group"+grp+".txt",(/vc_arr_all(gr_ind,vc_ind)/),"integer")
  end if

;========================================================================
; Before main loops, calculate max number of times in any one VC phase 
;========================================================================

; End minus start time indices for selected VC phase 
  nts  = max( (tf0 - ts0) + 1)
  nsim = dimsizes(dat)

; Total number of elements in time dimension (sims * times)
  ntot  = nts * nsim
  print("Upper limit on times analysed: "+ntot)

; Define constants and important array sizes 
  num_l      = 59 
  pii        = 3.14159265
  radii      = fspan(0.0,ar,nr)      ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,nt)   ; Array of azimuth angles (0, 5, 10, etc)
  thetas_rad = (pii/180)*thetas      ; Degrees to radians
  rad_size   = dimsizes(radii)-1     ; Index for later calculation

; Define constants (gas constant for dry air; specific heat of dry air at constant p)
  rd    = 287.0 
  cp0   = 1004.0

; 3D arrays to hold 'xpos' and 'ypos' values for all analysis times 
; 14/10/2019 --> for compatability with 'setup_cyl' and 'st_centre', do not split time dimension
  xpos_u_all  = new((/dimsizes(radii),dimsizes(thetas),ntot/),float)
  ypos_u_all  = new((/dimsizes(radii),dimsizes(thetas),ntot/),float)
  xcen_u_all  = new((/dimsizes(radii),dimsizes(thetas),ntot/),float)
  ycen_u_all  = new((/dimsizes(radii),dimsizes(thetas),ntot/),float)

  xpos_l_all  = new((/dimsizes(radii),dimsizes(thetas),ntot/),float)
  ypos_l_all  = new((/dimsizes(radii),dimsizes(thetas),ntot/),float)
  xcen_l_all  = new((/dimsizes(radii),dimsizes(thetas),ntot/),float)
  ycen_l_all  = new((/dimsizes(radii),dimsizes(thetas),ntot/),float)

; Mass flux (all grid points --> upper; lower; difference)
  flux_u_int0 = new((/dimsizes(radii),dimsizes(thetas),nsim,nts/),float)
  flux_l_int0 = new((/dimsizes(radii),dimsizes(thetas),nsim,nts/),float)
  flux_d_int0 = new((/dimsizes(radii),dimsizes(thetas),nsim,nts/),float)

; Mass flux (azimuthally averaged --> upper; lower; difference)
  flux_u_azi0 = new((/dimsizes(radii),nsim,nts/),float)
  flux_l_azi0 = new((/dimsizes(radii),nsim,nts/),float)
  flux_d_azi0 = new((/dimsizes(radii),nsim,nts/),float)

; Troubleshooting --> vertical velocity (14/10/2019)
  vvel_u_azi0 = new((/dimsizes(radii),nsim,nts/),float)
  vvel_l_azi0 = new((/dimsizes(radii),nsim,nts/),float)
  vvel_d_azi0 = new((/dimsizes(radii),nsim,nts/),float)

; Mass flux (all grid points)
  flux_u_int0@description = "Vertical mass flux"
  flux_u_int0@units       = "kg m~S~-2~N~ s~S~-1~N~"
  flux_u_int0!0           = "rad"
  flux_u_int0!1           = "azi"
  flux_u_int0!2           = "sim"
  flux_u_int0!3           = "time"

  flux_l_int0@description = "Vertical mass flux"
  flux_l_int0@units       = "kg m~S~-2~N~ s~S~-1~N~"
  flux_l_int0!0           = "rad"
  flux_l_int0!1           = "azi"
  flux_l_int0!2           = "sim"
  flux_l_int0!3           = "time"

  flux_d_int0@description = "Vertical mass flux"
  flux_d_int0@units       = "kg m~S~-2~N~ s~S~-1~N~"
  flux_d_int0!0           = "rad"
  flux_d_int0!1           = "azi"
  flux_d_int0!2           = "sim"
  flux_d_int0!3           = "time"

; Mass flux (azimuthally-averaged)
  flux_u_azi0@description = "Vertical mass flux"
  flux_u_azi0@units       = "kg m~S~-2~N~ s~S~-1~N~"
  flux_u_azi0!0           = "rad"
  flux_u_azi0!1           = "sim"
  flux_u_azi0!2           = "time"

  flux_l_azi0@description = "Vertical mass flux"
  flux_l_azi0@units       = "kg m~S~-2~N~ s~S~-1~N~"
  flux_l_azi0!0           = "rad"
  flux_l_azi0!1           = "sim"  
  flux_l_azi0!2           = "time"

  flux_d_azi0@description = "Vertical mass flux"
  flux_d_azi0@units       = "kg m~S~-2~N~ s~S~-1~N~"
  flux_d_azi0!0           = "rad"
  flux_d_azi0!1           = "sim"
  flux_d_azi0!2           = "time"

  vvel_u_azi0@description = "Vertical velocity"
  vvel_u_azi0@units       = "m s~S~-1~N~"
  vvel_u_azi0!0           = "rad"
  vvel_u_azi0!1           = "sim"
  vvel_u_azi0!2           = "time"

  vvel_l_azi0@description = "Vertical velocity"
  vvel_l_azi0@units       = "m s~S~-1~N~"
  vvel_l_azi0!0           = "rad"
  vvel_l_azi0!1           = "sim"
  vvel_l_azi0!2           = "time"

  vvel_d_azi0@description = "Vertical velocity"
  vvel_d_azi0@units       = "m s~S~-1~N~"
  vvel_d_azi0!0           = "rad"
  vvel_d_azi0!1           = "sim"
  vvel_d_azi0!2           = "time"

; Print out relevant information about storm tracking (sc=0,1,2,3,4,5,6)
  sc_arr                = (/"relative vorticity",\
                            "pressure",\          
                            "smoothed pressure",\ 
                            "extra smoothed pressure",\
                            "vort_min during symmetric phase",\
                            "windspeed minimum",\
                            "geopotential"/)

  print_clock("Using "+sc_arr(sc)+" to calculate storm centre position...")

;==========================================================
; Now start loop over these chosen input files and times
;==========================================================

; 14/10/2019 --> for compatability with 'setup_cyl' and 'st_centre', do not split time dimension
  centre_new_l = new( (/7, 1, ntot, 2/), float)
  centre_new_u = new( (/7, 1, ntot, 2/), float)

; Initialise counter variable (time)
  ct = 0

; Initialise counter variable (missing mass flux values)
  ms = 0 

 do st = 0, dimsizes(dat)-1 

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

   print_clock("Starting main part of script...")
   print("Working on ensemble member: "+ens0(st)+" ("+dat(st)+")" )

 ; Size of lat/lon arrays 
   llbox       = toint(dist*50)
   dsize       = (/llbox,llbox/)
   d0	       = dsize(0)
   d1	       = dsize(1)

   setvalues NhlGetWorkspaceObjectId
     "wsMaximumSize" : 1000000000
   end setvalues

 ; Choose which storm track to read in
   if (cn0 .eq. "slp") then
    cn = 0
   elseif (cn0 .eq. "vort") then
    cn = 1
   elseif (cn0 .eq. "geo") then
    cn = 2
   elseif (cn0 .eq. "geo_sm") then
    cn = 3
   end if

 ; Abbreviate layer information for reading in data 
   if (lay .eq. 1) then 
    lr0    = mlev0-int
    lr1    = mlev0+int
    ur0    = mlev1-int
    ur1    = mlev1+int
   end if 

 ; Read in 'pc' [u,v,vort] and 'pd' streams [w,p,t,theta,kmh,tau]
   c  = addfile(fili_c,"r")
   d  = addfile(fili_d,"r")

 ; Initialise additional counter variable (reset to zero for each simulation) 
   dt = 0 

;==================================
; Start loop over input files 
;==================================

 ; Loop over times in file
   do it = ts0(st), tf0(st)

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
    hyb_m  = sprintf("%0.0f",hybC(:))

  ; Calculate model level height information
    hy0    = sprintf("%0.0f",hybC(mlev0))
    hy1    = sprintf("%0.0f",hybC(mlev1))

  ; Read in horizontal wind components 
    u_varname = "x-wind"
    v_varname = "y-wind"

  ; Read in values on two single levels (lay=0) or over two layers (lay=1)
    if (lay .eq. 1) then
     ul0 = c->$u_varname$(it,lr0:lr1,{t0:t1},{n0:n1})
     vl0 = c->$v_varname$(it,lr0:lr1,{t0:t1},{n0:n1})

     uu0 = c->$u_varname$(it,ur0:ur1,{t0:t1},{n0:n1})
     vu0 = c->$v_varname$(it,ur0:ur1,{t0:t1},{n0:n1})
    else
     ul0 = c->$u_varname$(it,mlev0,{t0:t1},{n0:n1})
     vl0 = c->$v_varname$(it,mlev0,{t0:t1},{n0:n1})

     uu0 = c->$u_varname$(it,mlev1,{t0:t1},{n0:n1})
     vu0 = c->$v_varname$(it,mlev1,{t0:t1},{n0:n1})
    end if

;===================================
; Calculate storm-relative winds
;===================================

  ; Upper level
    ul = ul0 - u_cyc(cn,it-1)
    vl = vl0 - v_cyc(cn,it-1)

  ; Lower level
    uu = uu0 - u_cyc(cn,it-1)
    vu = vu0 - v_cyc(cn,it-1)

  ; Add metadata from original horizontal wind arrays
    copy_VarCoords(ul0, ul)
    copy_VarAtts(ul0, ul)

    copy_VarCoords(vl0, vl)
    copy_VarAtts(vl0, vl)

    copy_VarCoords(uu0, uu)
    copy_VarAtts(uu0, uu)

    copy_VarCoords(vu0, vu)
    copy_VarAtts(vu0, vu)

  ; Calculate relative vorticity using centered finite differences
  ; Option '2' --> boundary points estimated using one-sided difference scheme
    vort_l = uv2vr_cfd(ul, vl, latC, lonC, 2)
    copy_VarCoords(ul, vort_l)
    vort_l@units = "s~S~-1~N~"
    vort_l@name = "Relative vorticity on model levels"

    vort_u = uv2vr_cfd(uu, vu, latC, lonC, 2)
    copy_VarCoords(uu, vort_u)
    vort_u@units = "s~S~-1~N~"
    vort_u@name = "Relative vorticity on model levels"

  ; Coriolis parameter and absolute vorticity
    f1   = coriolis_param(latC)
    f0   = conform_dims(dimsizes(vort_l),f1,1)
    copy_VarMeta(vort_l,f0)
    copy_VarCoords(vort_l,f0)
    f0@description = "Coriolis parameter"
    f0@name        = "Coriolis parameter"

  ; Don't need two 'avo' variables, right?
    avo  = vort_l + f0
    copy_VarMeta(vort_l,avo)
    copy_VarAtts(vort_l,avo)
    copy_VarCoords(vort_l,avo)
    avo@description = "Absolute vorticity"
    avo@name        = "Absolute vorticity on model levels"
    avo@long_name   = "Absolute vorticity"

;===========================================================
; Calculate rate of strain (following Nguyen et al. 2011)
;===========================================================

    dim1   = "latitude"
    dim2   = "longitude"

  ; Calculate horizontal [lat,lon] derivatives of 'u' (du_dx, du_dy)
    dul    = grad_latlon_cfd(ul, ul&$dim1$, ul&$dim2$, False, False)
    dul_dy = dul[0]
    dul_dx = dul[1]
    delete(dul)

    duu    = grad_latlon_cfd(uu, uu&$dim1$, uu&$dim2$, False, False)
    duu_dy = duu[0]
    duu_dx = duu[1]
    delete(duu)

  ; Calculate horizontal [lat,lon] derivatives of 'v' (dv_dx, dv_dy)
    dvl     = grad_latlon_cfd(vl, vl&$dim1$, vl&$dim2$, False, False)
    dvl_dy  = dvl[0]
    dvl_dx  = dvl[1]
    delete(dvl)

    dvu     = grad_latlon_cfd(vu, vu&$dim1$, vu&$dim2$, False, False)
    dvu_dy  = dvu[0]
    dvu_dx  = dvu[1]
    delete(dvu)

  ; Calculate strain rate components --> sqrt[ (str)^2 + (shr)^2 ]
  ; Stretching (du_dx - dv_dy) and shearing (dv_dx + du_dy) deformation
    str1a_l = dul_dx - dvl_dy
    str2a_l = dvl_dx + dul_dy
    str1_l  = str1a_l ^ 2
    str2_l  = str2a_l ^ 2

    str1a_u = duu_dx - dvu_dy
    str2a_u = dvu_dx + duu_dy
    str1_u  = str1a_u ^ 2
    str2_u  = str2a_u ^ 2

  ; Finally, calculate strain rate 
    if (lay .eq. 1) then 
     str_l  = sqrt( str1_l(:,:,:) + str2_l(:,:,:) )
     str_u  = sqrt( str1_u(:,:,:) + str2_u(:,:,:) )
    else
     str_l  = sqrt( str1_l(:,:) + str2_l(:,:) )
     str_u  = sqrt( str1_u(:,:) + str2_u(:,:) )
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

     p_l   = d->p(it-2,lr0:lr1,{t0:t1},{n0:n1})
     w_l   = d->dz_dt(it-2,lr0:lr1,{t0:t1},{n0:n1})
     th_l  = d->theta(it-2,lr0:lr1,{t0:t1},{n0:n1})
     t_l   = d->temp(it-2,lr0:lr1,{t0:t1},{n0:n1})

     p_u   = d->p(it-2,ur0:ur1,{t0:t1},{n0:n1})
     w_u   = d->dz_dt(it-2,ur0:ur1,{t0:t1},{n0:n1})
     th_u  = d->theta(it-2,ur0:ur1,{t0:t1},{n0:n1})
     t_u   = d->temp(it-2,ur0:ur1,{t0:t1},{n0:n1})

    else

     p_l   = d->p(it-2,mlev0,{t0:t1},{n0:n1})
     w_l   = d->dz_dt(it-2,mlev0,{t0:t1},{n0:n1})
     th_l  = d->theta(it-2,mlev0,{t0:t1},{n0:n1})
     t_l   = d->temp(it-2,mlev0,{t0:t1},{n0:n1})

     p_u   = d->p(it-2,mlev1,{t0:t1},{n0:n1})
     w_u   = d->dz_dt(it-2,mlev1,{t0:t1},{n0:n1})
     th_u  = d->theta(it-2,mlev1,{t0:t1},{n0:n1})
     t_u   = d->temp(it-2,mlev1,{t0:t1},{n0:n1})

    end if 

  ; Calculate dry air density from pressure and temperature
    rho_l  = p_l / (rd * t_l)
    rho_u  = p_u / (rd * t_u)

  ; Add metadata
    copy_VarCoords(p_l, rho_l)
    rho_l@description = "Dry air density"
    rho_l@units = "kg m~S~-3~N~"

    copy_VarCoords(p_u, rho_u)
    rho_u@description = "Dry air density"
    rho_u@units = "kg m~S~-3~N~"

  ; Calculate mass flux (lower level/layer)
    mflux_l = rho_l * w_l
    copy_VarCoords(w_l, mflux_l)
    copy_VarAtts(w_l, mflux_l)
    mflux_l@units = "kg m~S~-2~N~ s~S~-1~N~"
    mflux_l@standard_name = "vertical_mass_flux"
    mflux_l@long_name     = "Vertical mass flux"
    mflux_l@title         = "Vertical mass flux"
    mflux_l@name          = "vert_mflux"

  ; Calculate mass flux (upper level/layer)
    mflux_u = rho_u * w_u
    copy_VarCoords(w_u, mflux_u)
    copy_VarAtts(w_u, mflux_u)
    mflux_u@units = "kg m~S~-2~N~ s~S~-1~N~"
    mflux_u@standard_name = "vertical_mass_flux"
    mflux_u@long_name     = "Vertical mass flux"
    mflux_u@title         = "Vertical mass flux"
    mflux_u@name          = "vert_mflux"

  ; Grid spacing (lon, lat)
    dx0    = lonD(1) - lonD(0)
    dy0    = latD(1) - latD(0)

;==========================
; Loop over model levels 
;==========================

    hy0   = sprintf("%0.0f",hybD(mlev0))
    hy0_l = sprintf("%0.0f",hybD(mlev0-int))
    hy0_u = sprintf("%0.0f",hybD(mlev1-int))

    hy1   = sprintf("%0.0f",hybD(mlev1))
    hy1_l = sprintf("%0.0f",hybD(mlev0+int))
    hy1_u = sprintf("%0.0f",hybD(mlev1+int))

  ; Read in variables on model levels
    if (lay .eq. 1) then

     print("Working on lower layer between " + hy0_l+ " and "+hy1_l+" m AGL")
     print("Working on upper layer between " + hy0_u+ " and "+hy1_u+" m AGL")

   ; 'pc' stream (upper level/layer)
     uu_plane     = dim_avg_n_Wrap( uu(:,:,:),0)
     vu_plane     = dim_avg_n_Wrap( vu(:,:,:),0)
     vrt_u_plane  = dim_avg_n_Wrap( vort_u(:,:,:),0)

     str_u_plane  = dim_avg_n_Wrap( str_u(:,:,:),0)
     copy_VarCoords(uu(0,:,:),str_u_plane)

   ; 'pc' stream (lower level/layer)
     ul_plane     = dim_avg_n_Wrap( ul(:,:,:),0)
     vl_plane     = dim_avg_n_Wrap( vl(:,:,:),0)
     vrt_l_plane  = dim_avg_n_Wrap( vort_l(:,:,:),0)

     str_l_plane  = dim_avg_n_Wrap( str_l(:,:,:),0)
     copy_VarCoords(ul(0,:,:),str_l_plane)

   ; 'pd' stream (upper level/layer)
     wu_plane     = dim_avg_n_Wrap( w_u(:,:,:),0)
     prs_u_plane  = dim_avg_n_Wrap( p_u(:,:,:), 0)
     th_u_plane   = dim_avg_n_Wrap( th_u(:,:,:),0)

     flux_u_plane = dim_avg_n_Wrap( mflux_u(:,:,:),0)
     copy_VarCoords(uu(0,:,:),flux_u_plane)

   ; 'pd' stream (lower level/layer)
     wl_plane     = dim_avg_n_Wrap( w_l(:,:,:),0)
     prs_l_plane  = dim_avg_n_Wrap( p_l(:,:,:), 0)
     th_l_plane   = dim_avg_n_Wrap( th_l(:,:,:),0)

     flux_l_plane = dim_avg_n_Wrap( mflux_l(:,:,:),0)
     copy_VarCoords(ul(0,:,:),flux_l_plane)

   ; 25/09/2019: calculate both layer averages on Cartesian grid ('flux_l_plane', 'flux_u_plane')

    else

      print("Working on " +hy0+ " and " +hy1+ " m AGL")

    ; 'pc' stream (upper level)
      uu_plane     = uu(:,:)
      vu_plane     = vu(:,:)
      vrt_u_plane  = vort_u(:,:)

      str_u_plane  = str_u(:,:)
      copy_VarCoords(uu, str_u_plane)

    ; 'pc' stream (lower level)
      ul_plane     = ul(:,:)
      vl_plane     = vl(:,:)
      vrt_l_plane  = vort_l(:,:)

      str_l_plane  = str_l(:,:)
      copy_VarCoords(ul, str_l_plane)

    ; 'pd' stream (upper level)
      wu_plane     = w_u(:,:)
      prs_u_plane  = p_u(:,:)
      th_u_plane   = th_u(:,:)
      flux_u_plane = mflux_u(:,:)

    ; 'pd' stream (lower level)
      wl_plane     = w_l(:,:)
      prs_l_plane  = p_l(:,:)
      th_l_plane   = th_l(:,:)
      flux_l_plane = mflux_l(:,:)

    end if

  ; Calculate vector windspeed and smooth
    spd_u_plane    = sqrt( (uu_plane ^ 2) + (vu_plane ^ 2) )
    copy_VarCoords(uu_plane, spd_u_plane)
    vort_u_smth    = smth9_Wrap(vrt_u_plane, 0.5, 0.5, True)

    spd_l_plane    = sqrt( (ul_plane ^ 2) + (vl_plane ^ 2) )
    copy_VarCoords(ul_plane, spd_l_plane)
    vort_l_smth    = smth9_Wrap(vrt_l_plane, 0.5, 0.5, True)

  ; Combine all pc and pd stream variables into larger arrays
    pc_u_plane     = (/uu_plane, vu_plane, spd_u_plane, \
                       vrt_u_plane, vrt_u_plane, str_u_plane, str_u_plane, str_u_plane/)

    pc_l_plane     = (/ul_plane, vl_plane, spd_l_plane, \
                       vrt_l_plane, vrt_l_plane, str_l_plane, str_l_plane, str_l_plane/)

  ; Second 'prs_x_plane' is a placeholder
    pd_u_plane     = (/wu_plane, prs_u_plane, prs_u_plane, flux_u_plane/)
    pd_l_plane     = (/wl_plane, prs_l_plane, prs_l_plane, flux_l_plane/)

;=============================================================  
; Call external function to calculate storm centre position
;============================================================= 

    ; Set radial distance (º) when looking for vort min in ringlike phase   
      r0            = 0.10 

    ; Lower level/layer
      centre_arr_l  = st_centre(0, it, cn, r0, centre, vrt_l_plane, prs_l_plane, spd_l_plane, \
                                centre_new_l, dy0, dx0, "comp", ct, "pd", "off")

    ; Upper level/layer
      centre_arr_u  = st_centre(0, it, cn, r0, centre, vrt_u_plane, prs_u_plane, spd_u_plane, \
                                centre_new_u, dy0, dx0, "comp", ct, "pd", "off")

    ; Updated storm track information contained in 'centre_new'
      centre_new_l  = centre_arr_l[0]
      centre_new_u  = centre_arr_u[0]

    ; When accessing the arrays above, remember there's only one vertical level/layer 
      ilev          = 0 

;====================================================================== 
; Make sure that the storm 'centre' is not within the eyewall updraft
;====================================================================== 

    if (w_check .eq. 1) then

     print("Original method sampled eyewall updraft. Recalculating...")

   ; Retrieve coordinates of storm centre  
     cen_loc0_l = centre_new_l(sc,ilev,ct,0)
     cen_loc1_l = centre_new_l(sc,ilev,ct,1)
     cen_loc0_u = centre_new_u(sc,ilev,ct,0)
     cen_loc1_u = centre_new_u(sc,ilev,ct,1)

   ; Coordinates of grid points either side
   ; EDIT -- lines below could be edited to increase the area (± 2/3dy0, etc)
     cen_0n_l = cen_loc0_l - dy0
     cen_0x_l = cen_loc0_l + dy0
     cen_1n_l = cen_loc1_l - dy0
     cen_1x_l = cen_loc1_l + dy0

     cen_0n_u = cen_loc0_u - dy0
     cen_0x_u = cen_loc0_u + dy0
     cen_1n_u = cen_loc1_u - dy0
     cen_1x_u = cen_loc1_u + dy0

   ; Calculate vertical velocity over small subset of grid points either side of original centre
   ; 13/10/2019 --> add divergence 
     cen_w_l  = wl_plane( {cen_0n_l:cen_0x_l}, {cen_1n_l:cen_1x_l} )
     cen_w_u  = wu_plane( {cen_0n_u:cen_0x_u}, {cen_1n_u:cen_1x_u} )

   ; If average over several grid points is above a threshold value, continue searching for centre 
   ; 14/10/2019 --> calculate averaged vertical velocity over small box...
     ave_w_l  = avg(cen_w_l)
     ave_w_u  = avg(cen_w_u)
     ave_out_l= sprintf("%0.1f",ave_w_l)
     ave_out_u= sprintf("%0.1f",ave_w_u)

   ; Coordinates of grid points either side (bigger grid than above) 
     delete([/cen_0n_l, cen_0x_l, cen_1n_l, cen_1x_l/])
     delete([/cen_0n_u, cen_0x_u, cen_1n_u, cen_1x_u/])

   ; 14/10/2019 --> ... but use larger box to constrain your search for new storm centre (L944)
   ; Lower level
     cen_0n_l = cen_loc0_l - (2 * dy0)
     cen_0x_l = cen_loc0_l + (2 * dy0)
     cen_1n_l = cen_loc1_l - (2 * dy0)
     cen_1x_l = cen_loc1_l + (2 * dy0)

   ; Higher level 
     cen_0n_u = cen_loc0_u - (2 * dy0)
     cen_0x_u = cen_loc0_u + (2 * dy0)
     cen_1n_u = cen_loc1_u - (2 * dy0)
     cen_1x_u = cen_loc1_u + (2 * dy0)

   ; Implement modified policy if there is even a hint of enhanced vertical velocity 
     if (ave_w_l .gt. w0 .or. ave_w_u .gt. w0) then

      print("ave. vertical velocity (lower) = "+ave_w_l+" ; ave. vvel (upper) = "+ave_w_u)

    ;===================================================
    ; First calculate new centre on lower level/layer
    ;===================================================

    ; Create smaller grid around previous centre (L923-926)
      w_pl_new  = wl_plane({cen_0n_l:cen_0x_l},{cen_1n_l:cen_1x_l})

    ; Find minimum vertical velocity on this smaller grid 
      w_min_new = min(w_pl_new)

    ; Reshape to 1D array, and find index of minimum vertical velocity 
      dims_w    = dimsizes(w_pl_new)
      w_1d      = ndtooned(w_pl_new)
      inds_w    = ind_resolve(minind(w_1d), dims_w)

    ; Retrieve lat/lon information from subset of grid 
      ltN       = w_pl_new&latitude
      lnN       = w_pl_new&longitude_1

      lat_min_w = ltN(0) + (dy0 * inds_w(0,0) )
      lon_min_w = lnN(0) + (dx0 * inds_w(0,1) )
      delete([/ltN, lnN, w_1d, w_pl_new/])

    ; Print new storm centre to screen 
      print("Centre (w_min): "+lat_min_w+" degrees N, "+lon_min_w+" degrees E")

    ; Replace values in 'centre_new' before calling 'setup_cyl' below 
      centre_new_l(sc,ilev,ct,0) = lat_min_w
      centre_new_l(sc,ilev,ct,1) = lon_min_w
      delete([/lat_min_w, lon_min_w/])

    ;===================================================
    ; Now calculate new centre on higher level/layer
    ;===================================================

    ; Create smaller grid around previous centre (L929-932)
      w_pl_new  = wu_plane({cen_0n_u:cen_0x_u},{cen_1n_u:cen_1x_u})

    ; Find minimum vertical velocity on this smaller grid
      w_min_new = min(w_pl_new)

    ; Reshape to 1D array, and find index of minimum vertical velocity
      dims_w    = dimsizes(w_pl_new)
      w_1d      = ndtooned(w_pl_new)
      inds_w    = ind_resolve(minind(w_1d), dims_w)

    ; Retrieve lat/lon information from subset of grid
      ltN       = w_pl_new&latitude
      lnN       = w_pl_new&longitude_1

      lat_min_w = ltN(0) + (dy0 * inds_w(0,0) )
      lon_min_w = lnN(0) + (dx0 * inds_w(0,1) )
      delete([/ltN, lnN, w_1d, w_pl_new/])

    ; Print new storm centre to screen
      print("Centre (w_min): "+lat_min_w+" degrees N, "+lon_min_w+" degrees E")

    ; Replace values in 'centre_new' before calling 'setup_cyl' below
      centre_new_u(sc,ilev,ct,0) = lat_min_w
      centre_new_u(sc,ilev,ct,1) = lon_min_w
      delete([/lat_min_w, lon_min_w/])

     end if

   ; Tidy up
     delete([/cen_w_l, cen_0n_l, cen_0x_l, cen_1n_l, cen_1x_l/])
     delete([/cen_w_u, cen_0n_u, cen_0x_u, cen_1n_u, cen_1x_u/])

   end if 

;=====================================================================  
; Call external function to switch to cylindrical coordinate system  
;=====================================================================  

    ; Call 'setup_cyl' and output pressure/vorticity on cylindrical grid 
    ; Option before "comp" --> basic ('0') or MBGT ('1') variables output    
    ; Option before 'num_l' --> 'ilev=0' 
    ; Run 'setup_cyl' for upper level/layer
      cyl_arr_u  = setup_cyl(ct, dist, ll_arr, lonC, latC, lonD, latD, lon1D, lat1D, \
                             pc_u_plane, pd_u_plane, ntot, centre_new_u, \
                             sc, sm, ar, nr, 0, num_l, 0, "comp", ct, 1, \
                             xpos_u_all, ypos_u_all, xcen_u_all, ycen_u_all)

      xpos_u_all = cyl_arr_u[0]
      ypos_u_all = cyl_arr_u[1]
      xcen_u_all = cyl_arr_u[2]
      ycen_u_all = cyl_arr_u[3] 

      lat_u_max  = cyl_arr_u[4]
      lon_u_max  = cyl_arr_u[5]
      radii_u    = cyl_arr_u[6]
      rad_u_size = cyl_arr_u[7]
      thetas_u   = cyl_arr_u[19] 

      vvel_u_int = cyl_arr_u[15]
      flux_u_int = cyl_arr_u[24]

    ; Also run 'setup_cyl' for lower level/layer
      cyl_arr_l  = setup_cyl(ct, dist, ll_arr, lonC, latC, lonD, latD, lon1D, lat1D, \
                             pc_l_plane, pd_l_plane, ntot, centre_new_l, \
                             sc, sm, ar, nr, 0, num_l, 0, "comp", ct, 1, \
                             xpos_l_all, ypos_l_all, xcen_l_all, ycen_l_all)

      xpos_l_all = cyl_arr_l[0]
      ypos_l_all = cyl_arr_l[1]
      xcen_l_all = cyl_arr_l[2]
      ycen_l_all = cyl_arr_l[3]

      lat_l_max  = cyl_arr_l[4]
      lon_l_max  = cyl_arr_l[5]
      radii_l    = cyl_arr_l[6]
      rad_l_size = cyl_arr_l[7]
      thetas_l   = cyl_arr_l[19]

      vvel_l_int = cyl_arr_l[15]
      flux_l_int = cyl_arr_l[24]

;===============================================
; Tidy up and input values into larger arrays 
;===============================================

  ; Define new arrays to hold azimuthally-averaged variables
    flux_u_azi = new( (/dimsizes(radii)/), float)
    flux_l_azi = new( (/dimsizes(radii)/), float)

  ; Also define arrays to hold azimuthally-averaged vertical velocity 
    vvel_u_azi = new( (/dimsizes(radii)/), float)
    vvel_l_azi = new( (/dimsizes(radii)/), float)

  ; Loop over radii and calculate the azimuthally-averaged mass flux
    do irad = sr, rad_size
     flux_u_azi(irad)   = avg( flux_u_int(irad,:) )
     flux_l_azi(irad)   = avg( flux_l_int(irad,:) )
     vvel_u_azi(irad)   = avg( vvel_u_int(irad,:) )
     vvel_l_azi(irad)   = avg( vvel_l_int(irad,:) )
    end do

  ; Find all indices where vertical velocity exceeds a threshold (15/10/2019)
    vvel_u_ind = ind(vvel_u_azi .lt. wt)
    vvel_l_ind = ind(vvel_l_azi .lt. wt)
    u_size     = dimsizes(vvel_u_ind)
    l_size     = dimsizes(vvel_l_ind)

  ; Ignore data points corresponding to vertical velocities less than our threshold (15/10/2019)
    if (.not. all(ismissing(vvel_u_ind) ) ) then 
     do ut = 0, u_size-1
      flux_u_azi(vvel_u_ind(ut) ) = flux_u_azi@_FillValue
     end do 
    end if 

    if (.not. all(ismissing(vvel_l_ind) ) ) then 
     do lt = 0, l_size-1
      flux_l_azi(vvel_l_ind(lt) ) = flux_l_azi@_FillValue
     end do 
    end if 

  ; Fill larger arrays with azimuthally-averaged mass flux
  ; 15/10/2019 --> skip time interval if all values are missing for either layer 
    if ( all(ismissing(flux_u_azi) .or. all(ismissing(flux_l_azi) ) ) ) then 
     print("Skip this time! Don't have data on both layers...")
     ms = ms + 1 
    else
     flux_u_azi0(:,st,dt)  = (/flux_u_azi(:)/)
     flux_l_azi0(:,st,dt)  = (/flux_l_azi(:)/)
    end if

  ; 14/10/2019 --> fill arrays below with azimuthally-averaged vertical velocity
    if (trb .eq. 1) then 
     vvel_u_azi0(:,st,dt)  = (/vvel_u_azi(:)/)
     vvel_l_azi0(:,st,dt)  = (/vvel_l_azi(:)/)
    end if 

  ; Tidy up
    delete([/flux_l_int, vvel_l_int, vvel_l_ind/])
    delete([/flux_u_int, vvel_u_int, vvel_u_ind/])
    delete([/wu_plane, vrt_u_plane, uu_plane, vu_plane, spd_u_plane, str_u_plane/])
    delete([/wl_plane, vrt_l_plane, ul_plane, vl_plane, spd_l_plane, str_l_plane/])
    delete([/prs_u_plane, th_u_plane, flux_u_plane, prs_l_plane, th_l_plane, flux_l_plane/])

    ct  = ct + 1		    ; Counter variable (time)

 ; Tidy up before next iteration (all variables from 'pc' and 'pd' streams)
   delete([/lonC,latC,lon1C,lat1C,ul0,vl0,uu0,vu0,ul,vl,uu,vu,vort_l,vort_u/])
   delete([/dul_dy, dul_dx, duu_dy, duu_dx, dvl_dy, dvl_dx, dvu_dy, dvu_dx, str_l, str_u/])
   delete([/lonD, latD, lon1D, lat1D, w_l, p_l, th_l, t_l, p_u, w_u, th_u, t_u/])
   delete([/rho_l, rho_u, mflux_l, mflux_u/])

  end do     ; End time loop (do it = 0, times-1)

 end do      ; End loop over input file times (do st = 0, dimsizes(dat)-1

;================================================
; Integrate over radius and write to text file
;================================================

; Summarise results from main loops above 
  print("Ignoring "+ms+" times from mass flux calculations (missing values on either layer)")

; 14/10/2019 --> Troubleshooting: do the same as above for vertical velocity 
  if (trb .eq. 1) then 

   vvel_u_all0 = dim_sum_n_Wrap(vvel_u_azi0(:,:,:), 0)
   vvel_l_all0 = dim_sum_n_Wrap(vvel_l_azi0(:,:,:), 0)

   vvel_d_all0 = vvel_u_all0 - vvel_l_all0

 ; 28/10/2019 --> Convert to 1D arrays before writing to text files 
   vvel_u_out0 = ndtooned(vvel_u_all0)
   vvel_l_out0 = ndtooned(vvel_l_all0)
   vvel_d_out0 = ndtooned(vvel_d_all0)

   if (diff .eq. 0) then 
    diri_u_out = "$ar/text/group"+grp+"_vvel_"+typ0+"_sc"+sc+"_lay"+lay+"_"\
                  +hy1+"_w"+w0+"_"+rs+"_"+rf+"_wt"+wt+".txt"
    asciiwrite(diri_u_out, vvel_u_out0)
    diri_l_out = "$ar/text/group"+grp+"_vvel_"+typ0+"_sc"+sc+"_lay"+lay+"_"\
                  +hy0+"_w"+w0+"_"+rs+"_"+rf+"_wt"+wt+".txt"
    asciiwrite(diri_l_out, vvel_l_out0)
   else
    diri_out   = "$ar/text/group"+grp+"_vvel_"+typ0+"_sc"+sc+"_lay"+lay+"_"\
                  +hy0+"_"+hy1+"_w"+w0+"_"+rs+"_"+rf+"_wt"+wt+".txt"
    asciiwrite(diri_out, vvel_d_out0)
   end if 

  else

 ; Calculate radial integral of azimuthally-averaged data (1D arrays --> "time")        
   flux_u_all0 = dim_sum_n_Wrap(flux_u_azi0(:,:,:), 0)
   flux_l_all0 = dim_sum_n_Wrap(flux_l_azi0(:,:,:), 0)

 ; Calculate difference between upper and lower layer                                   
   flux_d_all0 = flux_u_all0 - flux_l_all0

 ; Convert to 1D array before writing to text file                                      
   flux_d_out0 = ndtooned(flux_d_all0)

 ; 28/10/2019 --> also convert 'flux_u' and 'flux_l' to 1D arrays 
   flux_u_out0 = ndtooned(flux_u_all0)
   flux_l_out0 = ndtooned(flux_l_all0)

 ; Write out values to text file using 'asciiwrite'
   if (diff .eq. 0) then 
    diri_u_out = "$ar/text/group"+grp+"_"+typ0+"_sc"+sc+"_lay"+lay+"_"\
                 +hy1+"_w"+w0+"_"+rs+"_"+rf+"_wt"+wt+".txt"
    asciiwrite(diri_u_out, flux_u_out0)
    diri_l_out = "$ar/text/group"+grp+"_"+typ0+"_sc"+sc+"_lay"+lay+"_"\
                  +hy0+"_w"+w0+"_"+rs+"_"+rf+"_wt"+wt+".txt"
    asciiwrite(diri_l_out, flux_l_out0)
   else
    diri_out   = "$ar/text/group"+grp+"_"+typ0+"_sc"+sc+"_lay"+lay+"_"\
                  +hy0+"_"+hy1+"_w"+w0+"_"+rs+"_"+rf+"_wt"+wt+".txt"
    asciiwrite(diri_out, flux_d_out0)
   end if 

  end if

end 