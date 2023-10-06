; Calculate the vertical mass flux over two layers, and write out values to text files 

; Run using:

; ncl dist=0.6 nt=73 nr=11 ar=0.50 sr=0 wt=0.0 cn0=\"slp\" lay=1 mlev1=45
; int=2 sc=5 sm=0 trb=0 w_check=1 w0=0.3 run0=0 run=12 typ=\"all\" vc_mflux5_write.ncl

; 'dist'  = size of box following storm (degrees)
; 'full'  = beginning points (9); middle points (10); end points (11)
; 'ar'    = distance in degrees of outer radius (1.0, 2.0, etc)
; 'sr'    = radial circle number (1 = 5 km, 2 = 10 km, etc) to start flux calculations
; 'nr'    = number of radial circles between r = 0 and r = 'ar' (21, 41, 61, etc)
; 'nt'    = number of azimuth angles in cylindrical grid
; 'cn0'   = offline storm centre: "slp", "vort", "geo", "geo_sm"
; 'azi'   = (2) calculate azimuthal average and integrate radially
;           (1) calculate azimuthal average alone
;           (0) leave each grid point value alone
; 'w_check' = recalculate storm centre when sampling the eyewall updraft (w > 0)
; 'w0'    = 
; 'int'   = no. of levels either side of 'mlev1' to calculate layer average
; 'wt'    = vertical velocity threshold (e.g. 0.5 m s-1) when integrating radially 
; 'run0'  = calculate running average inline (1), or leave step out (0)

; 500 m [10], 1 km [14], 1.5 km [17], 2 km [20], 3 km [24], 4 km [28], 5 km [31], 6 km [34]
; 7 km [37], 8 km [39], 9 km [42], 10 km [45], 11 km [46], 12 km [48], 13 km [50], 14 km [52]  
; 15 km [54], 16 km [56]

; T+12    = 131 (144)         ; T+18  = 203 (216)          ; T+24    = 275 (288)  
; T+30    = 347 (360)         ; T+36  = 419 (432)          ; T+42    = 491 (504)   
; T+48    = 563 (576)         ; T+54  = 635 (648)          ; T+60    = 707 (720) 
; T+66    = 779 (792)         ; T+72  = 851 (864)          ; T+78    = 923 (936) 
; T+84    = 995 (1008)        ; T+90  = 1067 (1080)

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

; Which set of simulations do we want to analyse?                                            
  dat  = "02T12"
  ens0 = "em11"

; Ring-like phase (T+54 to T+56)                                                             
  if (typ .eq. "ring") then
   ts0  = 643
   tf0  = 674
; Ring-like to monopole transition (T+58 to T+60)                                            
  elseif (typ .eq. "r2m") then
   ts0  = 679
   tf0  = 710
; Monopole phase (T+63 to T+66)                                                              
  elseif (typ .eq. "mono") then
   ts0  = 739
   tf0  = 782
; Monopole to ring-like transition (T+71 to T+73)                                            
  elseif (typ .eq. "m2r") then
   ts0  = 835
   tf0  = 866
  elseif (typ .eq. "all") then 
   ts0  = 347
   tf0  = 995
  end if

; Starting minus ending time indices of all VC phases                                        
  ntot  = (tf0 - ts0) + 1

;======================================================================================      
; Create array to hold all values for composite (previous used with 'wrt=1' option)          
;======================================================================================      

; Dimensions of 5D array (levs; times; radial circles; diagnostics; sims)                    
  nlev      = 63
  ndiag     = 10

; Not sure if we need this array
  plot_arr0 = new( (/nlev, ntot, nr, ndiag/), "float")

;============================================================                                 
; Calculate storm motion using built-in function ('st_rm')                                    
;============================================================                                 

; 'dat'       = initialisation time (02T12, 03T00, ...)                                       
; 'ens0'      = ensemble simulation (em00, em01, ...)                                         
; 'diri'      = path to input files (see above)                                               
; 'dist'      = size of box following storm (degrees)                                         
; 'mins'      = analyse 1-h (0) or 5-min (1) data                                             

; Zonal; meridional; vector wind; times; lat; lon; centre                                     
  storm_rel = st_rm(dat, ens0, "$ar/text/", dist, 1)
  u_cyc     = storm_rel[0]
  v_cyc     = storm_rel[1]
  vel_cyc   = storm_rel[2]
  lat_arr   = storm_rel[4]
  lon_arr   = storm_rel[5]
  centre    = storm_rel[6]

;==========================================================                            
; Find 'pc' and 'pd' files (model height level data)                                   
;==========================================================                            

  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0+"/vc"
  fili_p      = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0

  fili_c1     = systemfunc("cd "+diri+" ; ls "+fili_p+"_pc.nc")
  fili_d1     = systemfunc("cd "+diri+" ; ls "+fili_p+"_pd.nc")
  fili_j1     = systemfunc("cd "+diri+" ; ls "+fili_p+"_pj.nc")
  fili_k1     = systemfunc("cd "+diri+" ; ls "+fili_p+"_pk.nc")

  fili_c      = diri+"/"+fili_c1
  fili_d      = diri+"/"+fili_d1
  fili_j      = diri+"/"+fili_j1
  fili_k      = diri+"/"+fili_k1

; Calculate number of times in each file using built-in function ('nc_times')          
  numINFO_c   = nc_times(fili_c)
  numINFO_d   = nc_times(fili_d)
  numINFO_j   = nc_times(fili_j)
  numINFO_k   = nc_times(fili_k)

; Number of times in each file                                                         
  numINPUT_c  = numINFO_c[0]
  numINPUT_d  = numINFO_d[0]
  numINPUT_j  = numINFO_j[0]
  numINPUT_k  = numINFO_k[0]

; Array of times from each file                                                        
  time_c      = numINFO_c[1]
  time_d      = numINFO_d[1]
  time_j      = numINFO_j[1]
  time_k      = numINFO_k[1]

; Create arrays for data and output information                                        
  llbox     = toint(dist*50)
  dsize     = (/llbox,llbox/)
  d0        = dsize(0)
  d1        = dsize(1)

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

; Create arrays to hold finished date strings                                       
  hr_min = new(numINPUT_c,string)
  min0   = (/4,9,14,19,24,29,34,39,44,49,54,59/)      ; Minute array (before)       
  min1   = (/5,10,15,20,25,30,35,40,45,50,55,0/)      ; Minute array (after)        

; Before correcting, find the indices of the elements where 'minute' = 59           
  hr_ind = ind(minute .eq. 59)

; Correct for errors in the code (round up values of 'minute' to multiples of 5)    
  do i = 0, numINPUT_c-1
   do m = 0, dimsizes(min0)-1

    if (minute(i) .eq. min0(m) ) then
     minute(i) = min1(m)
    end if

   end do
  end do

; Edit all incorrect values of 'hour' (where 'minute' previously = 59)              
  do j = 0, dimsizes(hr_ind)-1
   hour(hr_ind(j)) = hour(hr_ind(j)) + 1
  end do
  delete(j)

; Create string array of all values in 'minute'                                     
  mins = tostring(minute)

; Create 'hhmm' strings for each time (for output and title strings)                
  do i = 0, numINPUT_c-1

 ; Add zeros in front of single-digit hours                                         
   if (hour(i) .lt. 10) then

    if (minute(i) .lt. 10) then
     hr_min(i) = "0"+hour(i)+"0"+mins(i)
    elseif (minute(i) .ge. 10) then
     hr_min(i) = "0"+hour(i)+mins(i)
    end if

   else

  ; Also add zeros in front of single-digit minutes                                 
    if (minute(i) .lt. 10) then
     hr_min(i) = hour(i)+"0"+mins(i)
    elseif (minute(i) .ge. 10) then
     hr_min(i) = hour(i)+mins(i)
    end if

   end if

  end do
  delete(i)

; 18/12/2019 --> EDIT for consistency with 'tc5_ring_mono.ncl'
  time_str  = new(numINPUT_c,string)
  time_arr  = new(numINPUT_c,string)
  title_arr = new(numINPUT_c,string)

  time_str  = hr_min + " UTC " + sprinti("%0.2i ", day) + month_abbr(month)
  time_arr  = sprinti("%0.2i", day) + month_abbr(month) + "_" + hr_min + "Z"
  title_arr = "Valid at "+time_str

;=============================================================
; Define constants and create arrays before reading in data
;=============================================================

; String containing info on width of radial averaging band  
  ar0  = ar * 100
  rf = sprintf("%0.0fkm", ar0)

  sr0  = sr * 5
  rs = sprintf("%0.0fkm", sr0)

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

; Mass flux (azimuthally averaged --> upper; lower; difference)
  flux_u_azi0 = new((/dimsizes(radii),ntot/),float)

; Troubleshooting --> vertical velocity (14/10/2019)
  vvel_u_azi0 = new((/dimsizes(radii),ntot/),float)

; Mass flux (azimuthally-averaged)
  flux_u_azi0@description = "Vertical mass flux"
  flux_u_azi0@units       = "kg m~S~-2~N~ s~S~-1~N~"
  flux_u_azi0!0           = "rad"
  flux_u_azi0!1           = "time"

  vvel_u_azi0@description = "Vertical velocity"
  vvel_u_azi0@units       = "m s~S~-1~N~"
  vvel_u_azi0!0           = "rad"
  vvel_u_azi0!1           = "time"

; Print out relevant information about storm tracking (sc=0,1,2,3,4,5,6)
  sc_arr                  = (/"relative vorticity",\
                              "pressure",\          
                              "smoothed pressure",\ 
                              "extra smoothed pressure",\
                              "vort_min during symmetric phase",\
                              "windspeed minimum",\
                              "geopotential"/)

  print_clock("Using "+sc_arr(sc)+" to calculate storm centre position...")

 ; 14/10/2019 --> for compatability with 'setup_cyl' and 'st_centre', do not split time dimension
   centre_new_u = new( (/7, 1, ntot, 2/), float)

 ; Initialise counter variable (missing mass flux values)
   ms = 0 

   setvalues NhlGetWorkspaceObjectId
     "wsMaximumSize" : 1000000000
   end setvalues

 ; Choose which storm track to read in
   if (cn0 .eq. "slp") then
    cn = 0
   elseif (cn0 .eq. "slpf") then
    cn = 1
   else
    print("Wrong storm tracking method used! Exiting...")
    exit()
   end if

 ; Abbreviate layer information for reading in data 
   if (lay .eq. 1) then 
    ur0    = mlev1-int
    ur1    = mlev1+int
   end if 

 ; Read in 'pc' [u,v,vort] and 'pd' streams [w,p,t,theta,kmh,tau]
   c  = addfile(fili_c,"r")              ; Read in 'pc' stream [u,v,vort] 
   d  = addfile(fili_d,"r")              ; Read in 'pd' stream [w,p] 
   j  = addfile(fili_j,"r")              ; Read in 'pj' stream [t,theta,mixing ratio] 
   k  = addfile(fili_k,"r")              ; Read in 'pk' stream [tau,PV,w]  

 ; Initialise additional counter variable (reset to zero for each simulation) 
   ct = 0 

 ; Loop over chosen times 
   do it = ts0, tf0

    print("Working on time: "+time_str(it)+" (ct = "+ct+")" )
    title_arr(it) = "Valid at "+time_str(it)+" (T+"+it+")"

;=====================================================
; Read in variables from 'pc' stream (model levels)
;=====================================================

    t0    = lat_arr(cn,it,0)
    t1    = lat_arr(cn,it,1)
    n0    = lon_arr(cn,it,0)
    n1    = lon_arr(cn,it,1)
    print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

  ; Combine values above into array for input into external function 
    ll_arr = (/t0,t1,n0,n1/)

  ; Read in basic variables from 'pc' stream
    lonC   = c->longitude({n0:n1})                ; longitude ['d1' grid points]
    latC   = c->latitude({t0:t1})                 ; latitude  ['d0' grid points]
    hybC   = c->hybrid_ht(:) 		          ; 63 model ('theta') levels
    hyb_plot = hybC / 1000 
    hyb_plot@units = "km"
    hyb_m  = sprintf("%0.0f",hybC(:))

  ; Calculate model level height information
    hy1    = sprintf("%0.0f",hybC(mlev1))

  ; Read in horizontal wind components 
    u_varname = "x-wind"
    v_varname = "y-wind"

  ; Read in values on two single levels (lay=0) or over two layers (lay=1)
    if (lay .eq. 1) then
     uu0 = c->$u_varname$(it,ur0:ur1,{t0:t1},{n0:n1})
     vu0 = c->$v_varname$(it,ur0:ur1,{t0:t1},{n0:n1})
    else
     uu0 = c->$u_varname$(it,mlev1,{t0:t1},{n0:n1})
     vu0 = c->$v_varname$(it,mlev1,{t0:t1},{n0:n1})
    end if

;===================================
; Calculate storm-relative winds
;===================================

  ; Upper level
    uu = uu0 - u_cyc(cn,it)
    vu = vu0 - v_cyc(cn,it)

  ; Add metadata from original horizontal wind arrays
    copy_VarCoords(uu0, uu)
    copy_VarAtts(uu0, uu)

    copy_VarCoords(vu0, vu)
    copy_VarAtts(vu0, vu)

  ; Calculate relative vorticity using centered finite differences
  ; Option '2' --> boundary points estimated using one-sided difference scheme
    vort_u = uv2vr_cfd(uu, vu, latC, lonC, 2)
    copy_VarCoords(uu, vort_u)
    vort_u@units = "s~S~-1~N~"
    vort_u@name = "Relative vorticity on model levels"

  ; Coriolis parameter and absolute vorticity
    f1   = coriolis_param(latC)
    f0   = conform_dims(dimsizes(vort_u),f1,1)
    copy_VarMeta(vort_u,f0)
    copy_VarCoords(vort_u,f0)
    f0@description = "Coriolis parameter"
    f0@name        = "Coriolis parameter"

  ; Don't need two 'avo' variables, right?
    avo  = vort_u + f0
    copy_VarMeta(vort_u,avo)
    copy_VarAtts(vort_u,avo)
    copy_VarCoords(vort_u,avo)
    avo@description = "Absolute vorticity"
    avo@name        = "Absolute vorticity on model levels"
    avo@long_name   = "Absolute vorticity"

;===========================================================
; Calculate rate of strain (following Nguyen et al. 2011)
;===========================================================

    dim1   = "latitude"
    dim2   = "longitude"

  ; Calculate horizontal [lat,lon] derivatives of 'u' (du_dx, du_dy)
    duu    = grad_latlon_cfd(uu, uu&$dim1$, uu&$dim2$, False, False)
    duu_dy = duu[0]
    duu_dx = duu[1]
    delete(duu)

  ; Calculate horizontal [lat,lon] derivatives of 'v' (dv_dx, dv_dy)
    dvu     = grad_latlon_cfd(vu, vu&$dim1$, vu&$dim2$, False, False)
    dvu_dy  = dvu[0]
    dvu_dx  = dvu[1]
    delete(dvu)

  ; Calculate strain rate components --> sqrt[ (str)^2 + (shr)^2 ]
  ; Stretching (du_dx - dv_dy) and shearing (dv_dx + du_dy) deformation
    str1a_u = duu_dx - dvu_dy
    str2a_u = dvu_dx + duu_dy
    str1_u  = str1a_u ^ 2
    str2_u  = str2a_u ^ 2

  ; Finally, calculate strain rate 
    if (lay .eq. 1) then 
     str_u  = sqrt( str1_u(:,:,:) + str2_u(:,:,:) )
    else
     str_u  = sqrt( str1_u(:,:) + str2_u(:,:) )
    end if 

;====================================================== 
; Read in variables from 'pd' stream (model levels)
;======================================================

  ; Read in basic variables from 'pd' stream 
    hybD   = d->hybrid_ht(:)			    ; 63 model ('rho') levels

  ; Read in pressure and vertical velocity
    if (lay .eq. 1) then 

     p_u   = d->p(it,ur0:ur1,{t0:t1},{n0:n1})
     z_u   = d->ht(it,ur0:ur1,{t0:t1},{n0:n1})

    else

     p_u   = d->p(it,mlev1,{t0:t1},{n0:n1})
     z_u   = d->ht(it,mlev1,{t0:t1},{n0:n1})

    end if 

;======================================
; Read in variables from 'pj' stream
;======================================

  ; Read in temperature and potential temperature
    if (lay .eq. 1) then 
     t_u   = j->temp(it,ur0:ur1,{t0:t1},{n0:n1})
     th_u  = j->theta(it,ur0:ur1,{t0:t1},{n0:n1})
    else
     t_u   = j->temp(it,mlev1,{t0:t1},{n0:n1})
     th_u  = j->theta(it,mlev1,{t0:t1},{n0:n1})
    end if 

;======================================
; Read in variables from 'pk' stream
;======================================

  ; Read in basic variables 
    hybK   = k->hybrid_ht(:) 

  ; Alternative lat/lon grid 
    lon1K  = k->longitude_1({n0:n1})
    lat1K  = k->latitude_1({t0:t1})
    lonK   = k->longitude({n0:n1})
    latK   = k->latitude({t0:t1})

  ; Read in vertical velocity 
    if (lay .eq. 1) then 
     w_u   = k->dz_dt(it,ur0:ur1,{t0:t1},{n0:n1})
    else
     w_u   = k->dz_dt(it,mlev1,{t0:t1},{n0:n1})
    end if 

  ; Calculate dry air density from pressure and temperature
    rho_u  = p_u / (rd * t_u)

  ; Add metadata
    copy_VarCoords(p_u, rho_u)
    rho_u@description = "Dry air density"
    rho_u@units = "kg m~S~-3~N~"

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
    dx0    = lonC(1) - lonC(0)
    dy0    = latC(1) - latC(0)

;==========================
; Loop over model levels 
;==========================

    hy0_u = sprintf("%0.0f",hybD(mlev1-int))
    hy1   = sprintf("%0.0f",hybD(mlev1))
    hy1_u = sprintf("%0.0f",hybD(mlev1+int))

  ; Read in variables on model levels
    if (lay .eq. 1) then

     print("Working on layer between " + hy0_u+ " and "+hy1_u+" m AGL")

   ; 'pc' stream (upper level/layer)
     uu_plane     = dim_avg_n_Wrap( uu(:,:,:),0)
     vu_plane     = dim_avg_n_Wrap( vu(:,:,:),0)
     vrt_u_plane  = dim_avg_n_Wrap( vort_u(:,:,:),0)

     str_u_plane  = dim_avg_n_Wrap( str_u(:,:,:),0)
     copy_VarCoords(uu(0,:,:),str_u_plane)

   ; 'pd' stream (upper level/layer)
     wu_plane     = dim_avg_n_Wrap( w_u(:,:,:),0)
     prs_u_plane  = dim_avg_n_Wrap( p_u(:,:,:), 0)
     th_u_plane   = dim_avg_n_Wrap( th_u(:,:,:),0)

     flux_u_plane = dim_avg_n_Wrap( mflux_u(:,:,:),0)
     copy_VarCoords(uu(0,:,:),flux_u_plane)

    else

      print("Working on " +hy1+ " m AGL")

    ; 'pc' stream (upper level)
      uu_plane     = uu(:,:)
      vu_plane     = vu(:,:)
      vrt_u_plane  = vort_u(:,:)

      str_u_plane  = str_u(:,:)
      copy_VarCoords(uu, str_u_plane)

    ; 'pd' stream (upper level)
      wu_plane     = w_u(:,:)
      prs_u_plane  = p_u(:,:)
      th_u_plane   = th_u(:,:)
      flux_u_plane = mflux_u(:,:)

    end if

  ; Calculate vector windspeed and smooth
    spd_u_plane    = sqrt( (uu_plane ^ 2) + (vu_plane ^ 2) )
    copy_VarCoords(uu_plane, spd_u_plane)
    vort_u_smth    = smth9_Wrap(vrt_u_plane, 0.5, 0.5, True)

  ; Combine all pc and pd stream variables into larger arrays
    pc_u_plane     = (/uu_plane, vu_plane, spd_u_plane, vrt_u_plane, \
                       vrt_u_plane, str_u_plane, str_u_plane, str_u_plane/)

  ; Second 'prs_x_plane' is a placeholder
    pd_u_plane     = (/wu_plane, prs_u_plane, prs_u_plane, flux_u_plane/)

;=============================================================  
; Call external function to calculate storm centre position
;============================================================= 

    ; Set radial distance (º) when looking for vort min in ringlike phase   
      r0            = 0.10 

    ; Upper level/layer
      centre_arr_u  = st_centre(0, it, cn, r0, centre, vrt_u_plane, prs_u_plane, spd_u_plane, \
                                centre_new_u, dy0, dx0, "comp", ct, "pd", "off", "5min")

    ; Updated storm track information contained in 'centre_new'
      centre_new_u  = centre_arr_u[0]

;====================================================================== 
; Make sure that the storm 'centre' is not within the eyewall updraft
;====================================================================== 

    if (w_check .eq. 1) then

     print("Original method sampled eyewall updraft. Recalculating...")

   ; Retrieve coordinates of storm centre  
   ; 18/12/2019 --> argument 1 is '0', because 'centre_new_u' only has one vertical level
   ; Same as for call to 'setup_cyl' --> we've already vertically averaged over levels
   ; So the normal notation is redundant throughout this whole script 
     cen_loc0_u = centre_new_u(sc,0,ct,0)
     cen_loc1_u = centre_new_u(sc,0,ct,1)

   ; Coordinates of grid points either side
   ; EDIT -- lines below could be edited to increase the area (± 2/3dy0, etc)
     cen_0n_u = cen_loc0_u - dy0
     cen_0x_u = cen_loc0_u + dy0
     cen_1n_u = cen_loc1_u - dy0
     cen_1x_u = cen_loc1_u + dy0

   ; Calculate vertical velocity over small subset of grid points either side of original centre
   ; 13/10/2019 --> add divergence 
     cen_w_u  = wu_plane( {cen_0n_u:cen_0x_u}, {cen_1n_u:cen_1x_u} )

   ; If average over several grid points is above a threshold value, continue searching for centre 
   ; 14/10/2019 --> calculate averaged vertical velocity over small box...
     ave_w_u  = avg(cen_w_u)
     ave_out_u= sprintf("%0.1f",ave_w_u)

   ; Coordinates of grid points either side (bigger grid than above) 
     delete([/cen_0n_u, cen_0x_u, cen_1n_u, cen_1x_u/])

   ; 14/10/2019 --> ... but use larger box to constrain your search for new storm centre (L944)
   ; Higher level 
     cen_0n_u = cen_loc0_u - (2 * dy0)
     cen_0x_u = cen_loc0_u + (2 * dy0)
     cen_1n_u = cen_loc1_u - (2 * dy0)
     cen_1x_u = cen_loc1_u + (2 * dy0)

   ; Implement modified policy if there is even a hint of enhanced vertical velocity 
     if (ave_w_u .gt. w0) then

      print("ave. vertical velocity = "+ave_w_u)

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
    ; 18/12/2019 --> 'mlev1' corresponds to the level at the centre of the vertical layer
    ; ... previously this argument had been set to '0'
      centre_new_u(sc,0,ct,0) = lat_min_w
      centre_new_u(sc,0,ct,1) = lon_min_w
      delete([/lat_min_w, lon_min_w/])

     end if

   ; Tidy up
     delete([/cen_w_u, cen_0n_u, cen_0x_u, cen_1n_u, cen_1x_u/])

   end if 

;=====================================================================  
; Call external function to switch to cylindrical coordinate system  
;=====================================================================  

    ; Call 'setup_cyl' and output pressure/vorticity on cylindrical grid 
    ; Option before "comp" --> basic ('0') or MBGT ('1') variables output    
    ; 18/12/2019 --> 'mlev1' corresponds to the level at the centre of the vertical layer
    ; Run 'setup_cyl' for upper level/layer
      cyl_arr_u  = setup_cyl(it, dist, ll_arr, lonC, latC, lonK, latK, lon1K, lat1K, \
                             pc_u_plane, pd_u_plane, ntot, centre_new_u, \
                             sc, sm, ar, nr, nt, 0, num_l, 0, "comp", ct, 1, \
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

      if (trb .eq. 1) then 
       print("Upper level")
       printMinMax(flux_u_plane,False)
       printMinMax(flux_u_int,False)
       exit()
      end if 

;===============================================
; Tidy up and input values into larger arrays 
;===============================================

  ; Define new arrays to hold azimuthally-averaged variables
    flux_u_azi = new( (/dimsizes(radii)/), float)

  ; Also define arrays to hold azimuthally-averaged vertical velocity 
    vvel_u_azi = new( (/dimsizes(radii)/), float)

  ; Loop over radii and calculate the azimuthally-averaged mass flux
    do irad = sr, rad_size
     flux_u_azi(irad)   = avg( flux_u_int(irad,:) )
     vvel_u_azi(irad)   = avg( vvel_u_int(irad,:) )
    end do

  ; Find all indices where vertical velocity exceeds a threshold (15/10/2019)
    vvel_u_ind = ind(vvel_u_azi .lt. wt)
    u_size     = dimsizes(vvel_u_ind)

  ; Ignore data points corresponding to vertical velocities less than our threshold (15/10/2019)
  ; 28/10/2019 --> to bypass this partf of script, set 'wt' to a large negative value (e.g. -5)
    if (.not. all(ismissing(vvel_u_ind) ) ) then 
     do ut = 0, u_size-1
      flux_u_azi(vvel_u_ind(ut) ) = flux_u_azi@_FillValue
     end do 
    end if 

  ; Fill larger arrays with azimuthally-averaged mass flux
  ; 15/10/2019 --> skip time interval if all values are missing for either layer 
  ; 28/10/2019 --> have replaced 'dt' with 'it' [fill arrays with correct indices from 0-->119]
    if ( all(ismissing(flux_u_azi) ) ) then 
     print("Skip this time! Don't have data on both layers...")
     ms = ms + 1 
    else
     flux_u_azi0(:,ct)  = (/flux_u_azi(:)/)
    end if

  ; 14/10/2019 --> fill arrays below with azimuthally-averaged vertical velocity
    if (trb .eq. 1) then 
     vvel_u_azi0(:,ct)  = (/vvel_u_azi(:)/)
    end if 

  ; Tidy up
    delete([/flux_u_int, vvel_u_int, vvel_u_ind/])
    delete([/wu_plane, vrt_u_plane, uu_plane, vu_plane, spd_u_plane, str_u_plane/])
    delete([/prs_u_plane, th_u_plane, flux_u_plane/])

    ct  = ct + 1		    ; Counter variable (time)

 ; Tidy up before next iteration (all variables from 'pc' and 'pd' streams)
 ; Switched 'lonD' to 'lonK' (17/12/2019)
   delete([/lonC,latC,uu0,vu0,uu,vu,vort_u/])
   delete([/duu_dy, duu_dx, dvu_dy, dvu_dx, str_u/])
   delete([/lonK, latK, lon1K, lat1K, p_u, w_u, th_u, t_u/])
   delete([/rho_u, mflux_u/])

  end do     ; End time loop (do it = 0, times-1)

;=======================================================
; Integrate over radius and calculate running average
;=======================================================

; Summarise results from main loops above 
  print("Ignoring "+ms+" times from mass flux calculations (missing values on either layer)")

; Calculate radial integral of azimuthally-averaged data --> 2D arrays
  flux_u_out = dim_sum_n_Wrap(flux_u_azi0(:,:), 0)

;; Calculate running average for each simulation (on dimension 1 --> time) 
;  flux_u_out  = runave_n_Wrap(flux_u_all0(:), run, 0, 0)
  
; Write out to multiple text files, one for each simulation
  print("Writing to text file for: "+dat+" ("+ens0+")")

  if (run0 .eq. 1) then 
   diri_u_out = "$ar/text/flux5_"+typ+"_sc"+sc+"_lay"+lay+"_"\
                +hy1+"_w"+w0+"_"+rs+"_"+rf+"_wt"+wt+"_run"+run+".txt"
  else
   diri_u_out = "$ar/text/flux5_"+typ+"_sc"+sc+"_lay"+lay+"_"\
                +hy1+"_w"+w0+"_"+rs+"_"+rf+"_wt"+wt+".txt"
  end if 
  asciiwrite(diri_u_out, flux_u_out(:) )

end 