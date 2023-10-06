; Calculate the radius of maximum wind (model levels), and write to text file
; Written on 30/10/2019

; ncl dist=1.0 nt=73 ar=1.0 nr=21 cn0=\"geo_sm\" w_check=1 w0=0.3 sc=5
; sm=0 cn0=\"geo_sm\" typ=\"sym\" grp=\"7\" tc_write_ml_rmw_vtan.ncl

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'opt'   = output file format ("pdf" or "x11")
; 'ens0'  = ensemble member (em00 -- em11)
; 'dist'  = size of box following storm (degrees)
; 'ts0'	  = starting time (0–120)
; 'tf0'   = ending time (0–120)
; 'tan0'  = write out tangential wind (1) or RMW (0)

; Pressure levels (pres)
; (0) 1000, (1) 950, (2) 925, (3) 900, (4) 850, (5) 800, (6) 750, (7) 700
; (8) 650, (9) 600, (10) 500, (11) 400, (12) 300, (13) 250
; (14) 200, (15) 150, (16) 100

load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/st_centre.ncl"
load "$sam/ncl_func/setup_cyl.ncl"
load "$sam/ncl_func/ring_mono.ncl"

begin

;==============================
; Prelash (reading in files)
;==============================

; Define possible group dimension (30/10/2019)
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
  elseif (typ .eq. "test") then 
    dat  = (/"02T12","03T00"/)
    ens0 = (/"em00","em06"/)
    ts0  = (/75,40/)
    tf0  = (/76,41/)
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

; Ending time indices of all VC phases                                               
  tx   = tf0 + 1

; Starting time indices of all VC phases                                             
  tn   = ts0 - 1

; Starting minus ending time indices of all VC phases                                
  td   = (tx - tn) + 1

; Maximum number of times in any one VC phase                                        
  nts  = max(td)

; Number of individual simulations to analyse for composite                          
  nsim = dimsizes(dat)

; Total number of elements in time-dimension                                         
  ntot = sum( (tf0 - ts0) + 1)
  print("Total number of times in sample: "+ntot)

;======================================================================================  
; Create array to hold all values for composite (previous used with 'wrt=1' option)      
;======================================================================================

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

;===================================================
; Define arrays to hold tangential and total wind 
;===================================================

; Mean tangential wind {v}                                                            
  mean_vtan0       = new((/num_l,dimsizes(radii),ntot/),float)
  mean_vtan0!0     = "lev"
  mean_vtan0!1     = "rad"
  mean_vtan0!2     = "time"
  mean_vtan0@description = "Azimuthally averaged tangential wind"
  mean_vtan0@units = "m s~S~-1~N~"

; Mean total wind {v}
  mean_vtot0       = new((/num_l,dimsizes(radii),ntot/),float)
  mean_vtot0!0     = "lev"
  mean_vtot0!1     = "rad"
  mean_vtot0!2     = "time"
  mean_vtot0@description = "Azimuthally averaged total wind"
  mean_vtot0@units = "m s~S~-1~N~"

; Array to hold RMW information
  rmw_arr   = new( (/2,ntot/), float)
  max_arr   = new( (/2,ntot/), float)
  z_arr     = new( (/2,ntot/), float)

; Arrays to hold all 'xpos' and 'ypos' position arrays                                      
  xpos_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),ntot/),float)
  ypos_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),ntot/),float)

; Arrays to hold lat/lon at centre of cylindrical grid                                      
  xcen_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),ntot/),float)
  ycen_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),ntot/),float)

; Set counter variable before the start of any loops 
  ct        = 0 

;==========================================================                               
; Now start loop over these chosen input files and times                                  
;==========================================================  

 do st = 0, nsim-1

  dt0 = tf0(st) - ts0(st)

; Find 'pc' and 'pd' files (model height level data)                                      
  diri        = "/nfs/a319/earshar/"+dat(st)+"/"+ens0(st)
  fili_p      = "201607"+dat(st)+"00Z_NPTK_4p4_L80_ra1t_"+ens0(st)
  fili_c1     = systemfunc("cd "+diri+" ; ls "+fili_p+"_pc*.nc")
  fili_c      = diri+"/"+fili_c1
  fili_d1     = systemfunc("cd "+diri+" ; ls "+fili_p+"_pd*.nc")
  fili_d      = diri+"/"+fili_d1

; Calculate number of times in each file using built-in function ('nc_times')             
  numINFO_c   = nc_times(fili_c)
  numINFO_d   = nc_times(fili_d)

; Number/array of times from each file                                                    
  numINPUT_c  = numINFO_c[0]
  time_c      = numINFO_c[1]
  numINPUT_d  = numINFO_d[0]
  time_d      = numINFO_d[1]

;============================================================                             
; Calculate storm motion using built-in function ('st_rm')                                
;============================================================                             

; 'dat'       = initialisation time (02T12, 03T00, ...)                                   
; 'ens0'      = ensemble simulation (em00, em01, ...)                                     
; 'diri'      = path to input files (see above)                                           
; 'dist'      = size of box following storm (degrees)                                     
; 'mins'      = analyse 1-h (0) or 5-min (1) data                                         

; Zonal; meridional; vector wind; times; lat; lon; centre                                 
  storm_rel = st_rm(dat(st), ens0(st), "$ar/text/", dist, 0)
  u_cyc     = storm_rel[0]
  v_cyc     = storm_rel[1]
  vel_cyc   = storm_rel[2]
  numTIMES  = storm_rel[3]
  lat_arr   = storm_rel[4]
  lon_arr   = storm_rel[5]
  centre    = storm_rel[6]

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

  date_str = sprinti("%0.2iUTC ", hour) + \
             sprinti("%0.2i ", day) \
           + month_abbr(month)
  time_str = sprinti("%0.2i UTC ", hour) + \
             sprinti("%0.2i ", day) \
           + month_abbr(month)
  time_arr = sprinti("%0.2i", day) + \
             month_abbr(month) + \
             "_" + sprinti("%0.2iZ", hour)

  print_clock("Starting main part of script...")

;========================================                                                 
; Define cylindrical coordinate arrays                                                    
;========================================                                                 

  pii        = 3.14159265
  radii      = fspan(0.0,ar,nr)      ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,nt)   ; Array of azimuth angles (0, 5, 10, etc)            
  thetas_rad = (pii/180)*thetas      ; Degrees to radians                                 
  rad_size   = dimsizes(radii)-1     ; Index for later calculation                        
  num_l      = 59                    ; Number of model levels (edit to evolve smoothly) 

  print_clock("Working on ensemble member: "+ens0(st))

; Set size of arrays                                                               
  llbox      = toint(dist*50)
  dsize      = (/llbox,llbox/)
  d0         = dsize(0)
  d1         = dsize(1)

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

;==========================================================                              
; Read in required variables from 'pc' and 'pd' streams                                  
;==========================================================                              

; Create array to hold storm track positions using different methods                     
  centre_new  = new((/7, num_l, ntot, 2/),"float")

  c  = addfile(fili_c,"r")               ; Read in 'pc' stream [u,v,vort]                
  d  = addfile(fili_d,"r")               ; Read in 'pd' stream [w,p,t,theta,kmh,tau] 

;==================================                                                        
; Get the variables we will need                                                           
;================================== 

  do it = ts0(st), tf0(st)               ; Loop over times in file    

    print_clock("Looping between T+"+ts0(st)+" and T+"+tf0(st)+" ("+ens0(st)+", "+dat(st)+")")

    print("Working on time: "+time_str(it)+" (T+"+it+")" )
    title_arr(it) = "Valid at "+time_str(it)+" (T+"+it+")"

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
    hybC   = c->hybrid_ht(:)                      ; 63 model ('theta') levels              

    u_varname = "x-wind"
    v_varname = "y-wind"

  ; Read in horizontal wind components                                                     
  ; Both variables below -- [hybC | 63] * [latC | 150] * [lonC | 150]                      
    u0  = c->$u_varname$(it,:,{t0:t1},{n0:n1}) ; Zonal wind                                
    v0  = c->$v_varname$(it,:,{t0:t1},{n0:n1}) ; Meridional wind 

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

;==================================================================                        
; Calculate relative vorticity using centered finite differences                           
;==================================================================                        

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
    du0    = grad_latlon_cfd(u,u&$dim1$,u&$dim2$,False,False)
    du_dy  = du0[0]
    du_dx  = du0[1]
    delete(du0)

  ; Calculate horizontal derivatives of 'v' (dv_dx, dv_dy)                                 
    dv0    = grad_latlon_cfd(v,v&$dim1$,v&$dim2$,False,False)
    dv_dy  = dv0[0]
    dv_dx  = dv0[1]
    delete(dv0)

  ; Calculate strain rate --> sqrt[ (du_dx - dv_dy)^2 + (dv_dx + du_dy)^2 ]                

  ; Stretching deformation                                                                 
    str1a   = du_dx - dv_dy
  ; Shearing deformation                                                                   
    str2a   = dv_dx + du_dy
    str1    = (str1a)^2
    str2    = (str2a)^2

  ; Strain rate (3D array)                                                                 
    str0    = sqrt( str1(:,:,:) + str2(:,:,:) )

;==========================================================                                
; Calculate divergence using one-sided difference scheme                                   
;==========================================================                                

    div = uv2dv_cfd(u, v, latC, lonC, 2)
    copy_VarCoords(u, div)
    div@units = "s~S~-1~N~"
    div@name = "Divergence on pressure levels"

;======================================================                                    
; Read in variables from 'pd' stream (model levels)                                        
;======================================================                                    

  ; Read in basic variables from 'pd' stream                                               
    lon1D  = d->longitude_1({n0:n1})                ; longitude [xxx grid points]          
    lat1D  = d->latitude_1({t0:t1})                 ; latitude  [xxx grid points]          
    lonD   = d->longitude({n0:n1})                  ; longitude [xxx grid points]          
    latD   = d->latitude({t0:t1})                   ; latitude  [xxx grid points]          
    hybD   = d->hybrid_ht(:)                        ; 63 model ('rho') levels

  ; Read in temperature, pressure, theta, vertical velocity and geopotential height        
  ; All variables below -- [hybD | 63] * [latD | 150] * [lon1D | 151]                      
    t      = d->temp(it-2,:,{t0:t1},{n0:n1})        ; Temperature (K)                      
    p      = d->p(it-2,:,{t0:t1},{n0:n1})           ; Pressure (Pa)                        
    th     = d->theta(it-2,:,{t0:t1},{n0:n1})       ; Potential temperature (K)            
    w      = d->dz_dt(it-2,:,{t0:t1},{n0:n1})       ; Vertical velocity (m s-1)            
    z      = d->ht(it-2,:,{t0:t1},{n0:n1})          ; Geopotential height (m)  

  ; Calculate dry air density from pressure and temperature                                
    rho    = p / (rd * t)

  ; Add metadata                                                                           
    copy_VarCoords(p, rho)
    rho@description = "Dry air density"
    rho@units = "kg m~S~-3~N~"

  ; Also calculate specific volume (for use in diffusive tendency source terms)            
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

    dx0     = lonD(1) - lonD(0) ; Grid spacing (longitude)                                 
    dy0     = latD(1) - latD(0) ; Grid spacing (latitude) 

;==========================                                                         
; Loop over model levels                                                            
;==========================                                                         

  ; Loop over chosen model levels                                                   
    do ilev = 0, num_l-1

   ; Create output strings (model level info)                                       
     hy  = sprintf("%0.0f",hybD(ilev) )
     hy1 = toint(hy)
     hy0 = sprintf("%05g",hy1)
     print_clock("Working on model level " + ilev+ " ("+hy+" m AGL)")

   ; Read in other variables on model levels                                            
     u_plane                = u(ilev,:,:)          ; Zonal wind                         
     v_plane                = v(ilev,:,:)          ; Meridional wind                    
     w_plane                = w(ilev,:,:)          ; Vertical velocity                  
     f_plane                = f0(ilev,:,:)         ; Coriolis parameter                 
     rho_plane              = rho(ilev,:,:)        ; Density                            
     prs_plane              = p(ilev,:,:)          ; Pressure                           
     str_plane              = str0(ilev,:,:)       ; Strain rate                        
     vort_plane             = vort(ilev,:,:)       ; Relative vorticity                 
     avo_plane              = avo(ilev,:,:)        ; Absolute vorticity                 
     div_plane              = div(ilev,:,:)        ; Divergence  

   ; Calculate vector windspeed and smooth                                              
     spd_plane              = sqrt( (u_plane ^ 2) + (v_plane ^ 2) )
     copy_VarCoords(u_plane, spd_plane)
     vort_pl_smth           = smth9_Wrap(vort_plane, 0.5, 0.5, True)

   ; Combine all arrays above into single, larger array for input into 'setup_cyl'      

   ; 'pc' stream                                                                        
   ; full = 0 [u, v, spd, vort, avo, str, slp, div]                                     
   ; full = 1 [u, v, spd, vort, avo, f, str, div]                                       
     pc_plane               = (/u_plane, v_plane, spd_plane, vort_plane, \
                                avo_plane, f_plane, str_plane, div_plane/)

   ; 'pd' stream                                                                             
   ; full = 0 [w_plane, prs_plane, pv_plane, flux_plane]                                     
   ; full = 1 [w_plane, prs_plane, ..., ..., km_phlane, difv_xplane, difv_yplane, rho_plane]
     pd_plane               = (/w_plane, prs_plane, prs_plane, prs_plane, prs_plane, \
                                rho_plane, rho_plane, rho_plane/)

;=======================================================                                     
; Calculate storm centre position (external function)                                        
;=======================================================                                     

   ; Set radial distance (º) when looking for speed/vort min (in ringlike phase)             
     r0            = 0.12

   ; 08/10/2019 --> 'ct' = 53, 54, etc (counter variable doesn't start at 0...)              
     print_clock("Calculating storm centre position...")
     centre_arr    = st_centre(ilev, it, cn, r0, centre, vort_plane, prs_plane, spd_plane, \
                               centre_new, dy0, dx0, "comp", ct, "pd", "off")

   ; Updated storm track information contained in 'centre_new'                               
     centre_new    = centre_arr[0]

;======================================================================                      
; Make sure that the storm 'centre' is not within the eyewall updraft                        
;======================================================================                      

   if (w_check .eq. 1) then

   ; If vertical velocity is above a threshold value at the storm centre, try again          

   ; Retrieve coordinates of storm centre                                                    
     cen_loc0 = centre_new(sc,ilev,ct,0)
     cen_loc1 = centre_new(sc,ilev,ct,1)

   ; Coordinates of grid points either side                                                  
     cen_0n = cen_loc0 - dy0
     cen_0x = cen_loc0 + dy0
     cen_1n = cen_loc1 - dy0
     cen_1x = cen_loc1 + dy0

   ; Calculate vertical velocity and divergence at these grid points                         
     cen_w    = w_plane( {cen_0n:cen_0x}, {cen_1n:cen_1x} )
     cen_d    = div_plane( {cen_0n:cen_0x}, {cen_1n:cen_1x} )

   ; If average over several grid points is above a threshold value, continue searching for centre
     ave_w    = avg(cen_w)
     ave_d    = avg(cen_d)
     ave_out  = sprintf("%0.1f",ave_w)

   ; Coordinates of grid points either side (bigger grid than above)                 
     delete([/cen_0n, cen_0x, cen_1n, cen_1x/])
     cen_0n = cen_loc0 - (2 * dy0)
     cen_0x = cen_loc0 + (2 * dy0)
     cen_1n = cen_loc1 - (2 * dy0)
     cen_1x = cen_loc1 + (2 * dy0)

     if (ilev .gt. 0 .and. ave_w .gt. w0) then

      print("ave. vertical velocity = "+ave_w+" ; ave. divergence = "+ave_d)

    ; Create smaller grid around previous centre                                     
      w_pl_new  = w_plane({cen_0n:cen_0x},{cen_1n:cen_1x})

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
      centre_new(sc,ilev,ct,0) = lat_min_w
      centre_new(sc,ilev,ct,1) = lon_min_w

     end if

   ; Tidy up                                                                         
     delete([/cen_w, cen_d, cen_0n, cen_0x, cen_1n, cen_1x/])

   end if

;========================================================                            
; Switch to cylindrical coordinates (external function)                              
;========================================================  

   ; Last option = 0 --> basic variables output                                             
   ; Last option = 1 --> additional MBGT variables output                                   
   ; 08/10/2019 --> 'nts' option changed from 'numTIMES'. Currently does not impact output  
     print_clock("Translating to cylindrical grid...")
     cyl_arr  = setup_cyl(it, dist, ll_arr, lonC, latC, lonD, latD, lon1D, lat1D, \
                          pc_plane, pd_plane, nts, centre_new, \
                          sc, sm, ar, nr, ilev, num_l, 1, "comp", ct, \
                          0, xpos_all, ypos_all, xcen_all, ycen_all)

     xpos_all = cyl_arr[0]
     ypos_all = cyl_arr[1]
     xcen_all = cyl_arr[2]
     ycen_all = cyl_arr[3]
     lat_max  = cyl_arr[4]
     lon_max  = cyl_arr[5]
     radii    = cyl_arr[6]
     rad_size = cyl_arr[7]
     spd_int  = cyl_arr[10]
     vort_int = cyl_arr[11]

     vr_int   = cyl_arr[14]
     vt_int   = cyl_arr[15]

     print("Storm centre = "+lat_max+" ; "+lon_max)

;===================================================                                     
; Calculate mean (azimuthally-averaged) variables                                        
;===================================================                                     

    do irad = 0, dimsizes(radii)-1
     mean_vtan0(ilev,irad,ct) = avg(vt_int(irad,:))     ; [ {vtan} ]                  
     mean_vtot0(ilev,irad,ct) = avg(spd_int(irad,:))    ; [ {vtot} ]                  
    end do

  ; Tidy up (EDIT)
    delete([/spd_int,vr_int,vt_int,vort_int/])
    delete([/u_plane,v_plane,w_plane,f_plane,rho_plane/])
    delete([/prs_plane,str_plane,vort_plane,avo_plane,div_plane/])

  ; End loop over model levels
    end do

  ; Advance counter variable (time)
    ct  = ct + 1

  ; Tidy up (EDIT)
    delete([/u0,v0,vort,f0,f1,avo,du_dx,du_dy,dv_dx,dv_dy,str0,str1,str2,str1a,str2a,div/])
    delete([/t,p,th,w,z,rho,rho_a/])

 ; End loop over times
   end do

 ; Tidy up 
   delete([/numINPUT_c,numINPUT_d,time_c,time_d,u_cyc,v_cyc,vel_cyc/])
   delete([/utc_date,year,month,day,hour,minute,second,date_str,time_str,time_arr,title_arr/])

; End loop over simulations 
  end do

;===================================
; Now calculate RMW for all times
;===================================

; For each time, calculate the maximum tangential and total wind 
  do it = 0, ntot-1

   print("Working on time: "+it)

 ; Maximum tangential (0) and total (1) wind 
   max_arr(0,it) = max(mean_vtan0(:,:,it) ) 
   max_arr(1,it) = max(mean_vtot0(:,:,it) )

 ; Reshape 'mean_vtan0' and 'mean_vtot0' arrays to 1D
   dims_vtan     = dimsizes(mean_vtan0(:,:,it) )
   vtan_1d       = ndtooned(mean_vtan0(:,:,it) )
   dims_vtot     = dimsizes(mean_vtot0(:,:,it) )
   vtot_1d       = ndtooned(mean_vtot0(:,:,it) )

 ; Now find indices corresponding to wind maxima
   inds_vtan     = ind_resolve(maxind(vtan_1d), dims_vtan)
   inds_vtot     = ind_resolve(maxind(vtot_1d), dims_vtot)

 ; Retrieve corresponding height (model level) of wind maxima
   z_arr(0,it)   = hybD(inds_vtan(0,0) )
   z_arr(1,it)   = hybD(inds_vtot(0,0) )

 ; Also retrieve radius of maximum wind in both cases 
   rmw_arr(0,it) = radii(inds_vtan(0,1) )
   rmw_arr(1,it) = radii(inds_vtot(0,1) )

  end do 

;=================================
; Write out values to text file 
;=================================

  print_clock("Writing out variables to text files...")
  diri_out = "$ar/text/"

; Write out RMW
  rmw_out0 = diri_out+"rmw_vtan_"+typ0+"_group"+grp+"_sc"+sc+".txt"
  rmw_out1 = diri_out+"rmw_vtot_"+typ0+"_group"+grp+"_sc"+sc+".txt"
  asciiwrite(rmw_out0, rmw_arr(0,:) )
  asciiwrite(rmw_out1, rmw_arr(1,:) )

; Write out corresponding model level information
  z0_out0  = diri_out+"rmw_z_vtan_"+typ0+"_group"+grp+"_sc"+sc+".txt"
  z0_out1  = diri_out+"rmw_z_vtot_"+typ0+"_group"+grp+"_sc"+sc+".txt"
  asciiwrite(z0_out0, z_arr(0,:) )
  asciiwrite(z0_out1, z_arr(1,:) )

; Finally, write out actual wind speed values 
  max_out0 = diri_out+"rmw_v_vtan_"+typ0+"_group"+grp+"_sc"+sc+".txt"
  max_out1 = diri_out+"rmw_v_vtot_"+typ0+"_group"+grp+"_sc"+sc+".txt"
  asciiwrite(max_out0, max_arr(0,:) )
  asciiwrite(max_out1, max_arr(1,:) )
 
end 