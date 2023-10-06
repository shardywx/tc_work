; Produce panel plots of tangential and radial wind for different centre-finding methods

; Run using:

; ncl opt=\"png\" w0=0.5 z0=0 dist=1.4 clr=5 calc=1 rad0=0 nt=73 nr=22 ar=1.05 lgd0=\"tr\"
; cn0=\"slp\" sm=0 sc=5 w_check=1 typ=\"ring\" run=1 ri=12 tc_mbgt5_centre.ncl

; 'opt'   = output file format ("pdf" or "x11")
; 'dist'  = size of box following storm (degrees)
; 'clr'	  = colour map for plots (option 4 is used in Roger Smith's papers):
; 	    "ncl_default" [1],  "bdr_extra" [2], "amwg256" [3], "GMT_polar" [4]
; 'calc'  = centre calculated on each level (0), calculated offline (1) or using 980 m vort (2)
; 'rad0'  = radius of final plots: 150 km (0) or 200 km (1)
; 'typ'	  = VC phase ("sym","asym","stoa","atos")
; 'ar'    = distance in degrees of outer radius (1.0, 2.0, etc)
; 'nr'    = number of radial circles between r = 0 and r = 'ar' (21, 41, 61, etc)
; 'nt'    = number of azimuth angles in cylindrical grid
; 'cn0'   = offline storm centre: "slp", "vort", "geo", "geo_sm"

; 500 m [10], 1 km [14], 1.5 km [17], 2 km [20], 3 km [24], 4 km [28], 5 km [31], 6 km [34]  
; 7 km [37], 8 km [39], 9 km [42], 10 km [45], 11 km [46], 12 km [48], 13 km [50], 14 km [52]
; 15 km [54], 16 km [56]

; T+12    = 131 (144)         ; T+18  = 203 (216)          ; T+24    = 275 (288)
; T+30    = 347 (360)         ; T+36  = 419 (432)          ; T+42    = 491 (504)
; T+48    = 563 (576)         ; T+54  = 635 (648)          ; T+60    = 707 (720)
; T+66    = 779 (792)         ; T+72  = 851 (864)          ; T+78    = 923 (936)
; T+84    = 995 (1008)        ; T+90  = 1067 (1080)

; 'pc' stream --> 'it'
; 'centre'    --> 'it-1'
; 'pd' stream --> 'it-2'

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_centre.ncl"
load "$sam/ncl_func/setup_cyl.ncl"
load "$sam/ncl_func/setup_cyl_quick.ncl"
load "$sam/ncl_func/ring_mono.ncl"

begin

;==============================
; Prelash (reading in files)
;==============================

; Which set of simulations do we want to analyse? 
  dat  = "02T12"
  ens0 = "em11"

; Ring-like phase (T+54 to T+56)
  if (typ .eq. "ring") then 
   ts0  = 647 ; 635
   tf0  = 670 ; 658
; Ring-like to monopole transition (T+58 to T+60)
  elseif (typ .eq. "r2m") then 
   ts0  = 683
   tf0  = 706
; Monopole phase (T+63 to T+66)
  elseif (typ .eq. "mono") then 
   ts0  = 743
   tf0  = 778
; Monopole to ring-like transition (T+71 to T+73)
  elseif (typ .eq. "m2r") then 
   ts0  = 839
   tf0  = 862
  elseif (typ .eq. "test") then 
   ts0  = 647
   tf0  = 649
  end if 

; Starting minus ending time indices of all VC phases
  ntot  = (tf0 - ts0) + 1

;======================================================================================
; Create array to hold all values for composite (previous used with 'wrt=1' option)
;======================================================================================

; Dimensions of 5D array (levs; times; radial circles; diagnostics; sims)
  nlev      = 63
  ndiag     = 10

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

  if (run .eq. 1) then
   u_cyc    = runave_Wrap(storm_rel[0],ri,0)
   v_cyc    = runave_Wrap(storm_rel[1],ri,0)
   rstr     = "run"+ri
  else
   u_cyc    = storm_rel[0]
   v_cyc    = storm_rel[1]
  end if

  numTIMES  = storm_rel[3]
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
  llbox     = toint(dist*50)                   ; Calculate domain size (grid points)
  dsize     = (/llbox,llbox/)                  ; Domain size (grid points)
  d0        = dsize(0)
  d1        = dsize(1)

; Counter variable (equal to zero or 'ts0'?)
  ct        = ts0

  nh        = (tf0 - ts0) / 12                 ; Number of hours analysed (for averaging later)
  num_l     = 63                               ; Number of model levels
  num_t     = numINPUT_k                       ; Number of times in 5-min data files 

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

; Create array to hold finished date strings                                                
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

  time_str  = new(numINPUT_c,string)
  time_arr  = new(numINPUT_c,string)
  title_arr = new(numINPUT_c,string)

  time_str  = hr_min + " UTC " + sprinti("%0.2i ", day) + month_abbr(month)
  time_arr  = sprinti("%0.2i", day) + month_abbr(month) + "_" + hr_min + "Z"
  title_arr = "Valid at "+time_str

;==================================================================================
;;;;;;;;;;;;;;;;;;;;;;;;;;; Start of main part of script ;;;;;;;;;;;;;;;;;;;;;;;;;;
;==================================================================================

  print_clock("Starting main part of script...")

;========================================
; Define cylindrical coordinate arrays
;========================================

  pii        = 3.14159265
  radii      = fspan(0.0,ar,nr)      ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,nt)   ; Array of azimuth angles (0, 5, 10, etc)
  thetas_rad = (pii/180)*thetas      ; Degrees to radians
  rad_size   = dimsizes(radii)-1     ; Index for later calculation

; Define constants
  rd	     = 287.0	             ; Gas constant for dry air (J/kg/K^2)
  cp0	     = 1004.0                ; Specific heat of dry air at constant pressure (J/kg/K)
  kp	     = 0.286                 ; For Exner function calculations (0.286)
  pref	     = 1000.0                ; Reference pressure p0 (1000 hPa)

; Arrays to hold all 'xpos' and 'ypos' position arrays
  xpos_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),ntot/),float)
  ypos_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),ntot/),float)

; Arrays to hold lat/lon at centre of cylindrical grid 
  xcen_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),ntot/),float)
  ycen_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),ntot/),float)

  xpos_all!0          = "lev"
  xpos_all!3          = "time"
  
  ypos_all!0          = "lev"
  ypos_all!3          = "time"

;=========================
; Start multiple loops
;=========================

; Also set time-dependent counter variable (08/10/2019 EDIT)
; 'ct' refreshes to zero (0) at the end of each simulation loop (do st = ...)
  ct          = 0

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 1000000000
  end setvalues

;=======================================
; Choose which storm track to read in
;=======================================

  if (cn0 .eq. "slp") then
   cn = 0
  elseif (cn0 .eq. "slpf") then
   cn = 1
  end if

;==========================================================  
; Read in required variables from 'pc' and 'pd' streams
;========================================================== 

; Create array to hold storm track positions using different methods  
  centre_new  = new((/10, num_l, ntot, 2/),"float")

  c  = addfile(fili_c,"r")              ; Read in 'pc' stream [u,v,vort]                    
  d  = addfile(fili_d,"r")              ; Read in 'pd' stream [w,p]                         
  j  = addfile(fili_j,"r")              ; Read in 'pj' stream [t,theta,mixing ratio]        
  k  = addfile(fili_k,"r")              ; Read in 'pk' stream [tau,PV,w] 

; Loop over times in file
  do it = ts0, tf0

    print_clock("Looping between T+"+ts0+" and T+"+tf0+" ("+ens0+", "+dat+")")

    print("Working on time: "+time_str(it)+" (T+"+it+")" )
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
    lon    = c->longitude({n0:n1})                ; longitude ['d1' grid points]
    lat    = c->latitude({t0:t1})                 ; latitude  ['d0' grid points]  

  ; Horizontal grid spacing
    dy0    = lat(1) - lat(0)
    dx0    = lon(1) - lon(0)

  ; Model levels
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

    u = u0 - u_cyc(cn,it)
    v = v0 - v_cyc(cn,it)

  ; Add metadata from original horizontal wind arrays
    copy_VarCoords(u0,u)
    copy_VarAtts(u0,u)

    copy_VarCoords(v0,v)
    copy_VarAtts(v0,v)

;==================================================================
; Calculate relative vorticity using centered finite differences
;==================================================================

  ; Option '2' --> boundary points estimated using one-sided difference scheme
    vort = uv2vr_cfd(u, v, lat, lon, 2)
    copy_VarCoords(u, vort)
    vort@units = "s~S~-1~N~"
    vort@name = "Relative vorticity on model levels"

  ; Coriolis parameter and absolute vorticity
    f1   = coriolis_param(lat)
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
    copy_VarCoords(u, str0)
    str0@name = "Strain rate on model levels"

;==========================================================
; Calculate divergence using one-sided difference scheme
;==========================================================

    div = uv2dv_cfd(u, v, lat, lon, 2)
    copy_VarCoords(u, div)
    div@units = "s~S~-1~N~"
    div@name = "Divergence on model levels"

;======================================================                                 
; Read in variables from 'pd' stream (model levels)                                     
;======================================================                                 

  ; Read in basic variables from 'pd' stream                                            
    hybD   = d->hybrid_ht(:)                        ; 63 model ('rho') levels           

  ; Read in temperature, pressure, vertical velocity and geopotential height            
  ; All variables below -- [hybD | 63] * [latD | 150] * [lon1D | 151]                   
    p      = d->p(it,:,{t0:t1},{n0:n1})             ; Pressure (Pa)                     
    z      = d->ht(it,:,{t0:t1},{n0:n1})            ; Geopotential height (m)           

;======================================                                                 
; Read in variables from 'pj' stream                                                    
;======================================                                                 

  ; Read in basic variables from 'pj' stream                                            
    hybJ   = j->hybrid_ht(:)                        ; 63 model ('rho') levels           

  ; Read in temperature, potential temperature and mixing ratio                         
  ; All variables below  -- [hybJ | 63] * [latJ | 150]  * [lonJ | 150]                  
    t      = j->temp(it,:,{t0:t1},{n0:n1})          ; Temperature                       
    th     = j->theta(it,:,{t0:t1},{n0:n1})         ; Potential temperature             
    mix    = j->mix(it,:,{t0:t1},{n0:n1})           ; Mixing ratio 

;======================================                                                 
; Read in variables from 'pk' stream                                                    
;======================================                                                 

  ; Read in basic variables from 'pk' stream                                            
    hybK   = k->hybrid_ht(:)                        ; 63 model ('rho') levels           

  ; Alternative latitude/longitude grid                                                 
    lon1K  = k->longitude_1({n0:n1})
    lat1K  = k->latitude_1({t0:t1}) 
    lonK   = k->longitude({n0:n1})
    latK   = k->latitude({t0:t1})
    
   
  ; Read in potential vorticity and vertical velocity                                   
    w      = k->dz_dt(it,:,{t0:t1},{n0:n1})         ; Vertical velocity (m s-1)         
    pv     = k->field83(it,0:62,{t0:t1},{n0:n1})    ; Potential vorticity               

  ; Read in horizontal wind stress components and diffusivity                           
    str_xz = k->taux(it,:,{t0:t1},{n0:n1})
    str_yz = k->tauy(it,:,{t0:t1},{n0:n1})
    kmh    = k->unspecified(it,:,{t0:t1},{n0:n1})

;================================================================                       
; Also calculate dry air density from pressure and temperature                          
;================================================================                       

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

  ; Now convert units of pressure and geopotential height before continuing             
    p      = p / 100
    z      = z / 10

  ; Add metadata                                                                        
    p@units = "hPa"
    z@units = "dam"

;==========================
; Loop over model levels 
;==========================

  ; Loop over chosen model levels
    do ilev = 0, num_l-1

   ; Create output strings (model level info)
     hy  = sprintf("%0.0f",hybD(ilev) )
     hy1 = toint(hy)
     hy0 = sprintf("%05g",hy1) 

     hy2 = sprintf("%0.0f",hybD(z0) )
     hy3 = tofloat(hy2) / 1000
     print_clock("Working on model level " + ilev+ " ("+hy+" m AGL)")

   ; Read in other variables on model levels
     u_plane 	            = u(ilev,:,:)	   ; Zonal wind 
     v_plane		    = v(ilev,:,:)	   ; Meridional wind 
     w_plane		    = w(ilev,:,:)	   ; Vertical velocity
     f_plane		    = f0(ilev,:,:)	   ; Coriolis parameter
     rho_plane              = rho(ilev,:,:)	   ; Density
     prs_plane              = p(ilev,:,:)	   ; Pressure
     str_plane              = str0(ilev,:,:)       ; Strain rate
     vort_plane	   	    = vort(ilev,:,:)	   ; Relative vorticity
     avo_plane		    = avo(ilev,:,:)	   ; Absolute vorticity
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
                                prs_plane, prs_plane, rho_plane/)

;=====================================================                                         
; Interpolate onto finer grid before finding centre                                            
;=====================================================                                         

   ; get size information from 'pc' and 'pd' stream variables                                  
     u_size = dimsizes(u_plane)
     u_t    = u_size(0)
     u_n    = u_size(1)
     w_size = dimsizes(w_plane)
     w_t    = w_size(0)
     w_n    = w_size(1)

   ; create latitude and longitude arrays with finer grid spacing                              
     lat_pc = fspan(u_plane&latitude(0), u_plane&latitude(u_t-1), (u_t*4)-1)
     lon_pc = fspan(u_plane&longitude(0), u_plane&longitude(u_n-1), (u_n*4)-1)

   ; add metadata                                                                              
     lat_pc!0 = "latitude"
     lat_pc@units = "degrees_north"
     lon_pc!0 = "longitude"
     lon_pc@units = "degrees_east"

   ; interpolate the 'pc' stream variables from original to finer grid using 'linint2'
     u_pl_fine    = linint2_Wrap(u_plane&longitude, u_plane&latitude, u_plane, False, \
                                 lon_pc, lat_pc, 0)
     v_pl_fine    = linint2_Wrap(v_plane&longitude, v_plane&latitude, v_plane, False, \
                                 lon_pc, lat_pc, 0)
     spd_pl_fine  = linint2_Wrap(spd_plane&longitude, spd_plane&latitude, spd_plane, False, \
                                 lon_pc, lat_pc, 0)
     vort_pl_fine = linint2_Wrap(vort_plane&longitude, vort_plane&latitude, vort_plane, False, \
                                 lon_pc, lat_pc, 0)
     avo_pl_fine  = linint2_Wrap(avo_plane&longitude, avo_plane&latitude, avo_plane, False, \
                                 lon_pc, lat_pc, 0)
     str_pl_fine  = linint2_Wrap(str_plane&longitude, str_plane&latitude, str_plane, False, \
                                 lon_pc, lat_pc, 0)
     div_pl_fine  = linint2_Wrap(div_plane&longitude, div_plane&latitude, div_plane, False, \
                                 lon_pc, lat_pc, 0)
     f_pl_fine    = linint2_Wrap(f_plane&longitude, f_plane&latitude, f_plane, False, \
                                 lon_pc, lat_pc, 0)

   ; also interpolate the 'pd' stream variables onto the same grid
     w_pl_fine    = linint2_Wrap(w_plane&longitude_1, w_plane&latitude, w_plane, False, \
                                 lon_pc, lat_pc, 0)
     prs_pl_fine  = linint2_Wrap(prs_plane&longitude, prs_plane&latitude, prs_plane, False, \
                                 lon_pc, lat_pc, 0)
     rho_pl_fine  = linint2_Wrap(rho_plane&longitude, rho_plane&latitude, rho_plane, False, \
                                 lon_pc, lat_pc, 0)

   ; new grid spacing
     dy0          = abs(lat_pc(1) - lat_pc(0))
     dx0          = abs(lon_pc(1) - lon_pc(0))
    
;=======================================================
; Calculate storm centre position (external function)
;=======================================================

   ; Set radial distance (º) when looking for speed/vort min (in ringlike phase) 
     r0            = 0.05

   ; 08/10/2019 --> 'ct' = 53, 54, etc (counter variable doesn't start at 0...)
     print_clock("Calculating storm centre position...")
     centre_arr    = st_centre(ilev, it, cn, r0, centre, vort_pl_fine, prs_pl_fine, spd_pl_fine, \
                               centre_new, dy0, dx0, "comp", ct, "pd", "off", "5min")

   ; Updated storm track information contained in 'centre_new'
     centre_new = centre_arr[0] 

   ; Extract storm-centre calculations using (1) horizontal wind and (2) Laplacian of (1)
     centre_spd = centre_new(5,ilev,ct,:)
     centre_s2  = centre_new(6,ilev,ct,:)

     centre_psm = centre_new(7,ilev,ct,:)
     centre_ps2 = centre_new(8,ilev,ct,:)

;================================================================================
; Apply an algorithm to prevent the centre from moving too much between levels
;================================================================================

     if (ilev .ge. 1) then 
    ; calculate change in latitude and longitude between model levels
      d_lat = centre_new(sc,ilev,ct,0) - centre_new(sc,ilev-1,ct,0)
      d_lon = centre_new(sc,ilev,ct,1) - centre_new(sc,ilev-1,ct,1)
    ; print to screen if needed 
      print("change in lat = "+d_lat+", change in lon = "+d_lon)
    ; if storm centre has moved more than a threshold distance, recalculate 
;      if ( abs(d_lat) .gt. 0.02 .or. abs(d_lon) .gt. 0.02) then 
;       print("Change in storm centre position too large - using position from level below...")
;       centre_new(sc,ilev,ct,0) = centre_new(sc,ilev-1,ct,0)
;       centre_new(sc,ilev,ct,1) = centre_new(sc,ilev-1,ct,1)
;       ;centre_spd(0) = centre_new(sc,ilev-1,ct,0)
;       ;centre_spd(1) = centre_new(sc,ilev-1,ct,1)
;      end if 
     end if 

;===================================== 
; START MID-SCRIPT PLOTTING ROUTINE 
;===================================== 

   ; Print out x-y plots of u,v,w at selected heights (2,4,6,10,12,14 km)  
     ;if (ilev .eq. 28) then   

   ; Plot the data as we go  
     output = "$nep/nepartak/images/cart_coords/cxy_4panel_"+dat+\
              "_"+ens0+"_"+typ+"_"+hy0+"m_sc"+sc+"_T"+it
     wks    = gsn_open_wks(opt, output)

     print("ct = "+ct)

   ; Windspeed minimum resources 
     opts_spd                               = True
     opts_spd@cnFillOn                      = True
     opts_spd@cnLineLabelInterval           = 2.0
     opts_spd@cnLineLabelFontHeightF        = 0.012
     opts_spd@cnLineLabelBackgroundColor    = "transparent"
     opts_spd@cnLineLabelPlacementMode      = "constant"
     opts_spd@cnLinesOn                     = False
     opts_spd@cnInfoLabelOn                 = False
     opts_spd@cnLevelSelectionMode          = "ExplicitLevels"
     opts_spd@cnFillPalette                 = "radar_new"
;     opts_spd@cnLevels                      = (/5.0, 10.0, 15.0, 20.0,\
;                                                25.0, 30.0, 35.0, 40.0, \
;                                                45.0, 50.0, 55.0, 60.0, 65.0/)
     opts_spd@cnLevels                      = (/1.0, 2.0, 3.0, 5.0,\
                                                10.0, 15.0, 20.0, 30.0, \
                                                40.0, 50.0, 55.0, 60.0, 65.0/)
     opts_spd@cnFillColors                  = (/-1,0,1,2,3,4,5,6,7,8,9,10,11,12/)
     if (single .eq. 1) then 
      opts_spd@gsnMaximize                  = True 
     else
      opts_spd@gsnMaximize                  = False
     end if 

   ; No title or additional strings 
     opts_spd@tiMainString                  = ""
     opts_spd@gsnLeftString                 = ""
     opts_spd@gsnRightString                = ""

   ; Turn individual label bar on/off
     opts_spd@lbLabelBarOn                  = True
     opts_spd@lbBoxEndCapStyle              = "TriangleBothEnds"
     opts_spd@lbOrientation                 = "Vertical"

     opts_spd@trYMaxF                       = centre_new(5,ilev,ct,0) + ar
     opts_spd@trYMinF                       = centre_new(5,ilev,ct,0) - ar 
     opts_spd@trXMaxF                       = centre_new(5,ilev,ct,1) + ar 
     opts_spd@trXMinF                       = centre_new(5,ilev,ct,1) - ar

   ; Do not draw the plot or advance the frame 
     opts_spd@gsnDraw                       = False
     opts_spd@gsnFrame                      = False
     opts_spd@gsnAddCyclic                  = False

     opts_spd@tmYROn                        = False
     opts_spd@tmYLOn                        = True
     opts_spd@tmXTOn                        = False
     opts_spd@tmXBOn                        = True

   ; Resources for smoothed wind speed
     opts_ssm                               = opts_spd
     opts_ssm@trYMaxF                       = centre_new(6,ilev,ct,0) + ar
     opts_ssm@trYMinF                       = centre_new(6,ilev,ct,0) - ar
     opts_ssm@trXMaxF                       = centre_new(6,ilev,ct,1) + ar
     opts_ssm@trXMinF                       = centre_new(6,ilev,ct,1) - ar

     opts_sm2                               = opts_ssm
     opts_sm2@trYMaxF                       = centre_new(7,ilev,ct,0) + ar
     opts_sm2@trYMinF                       = centre_new(7,ilev,ct,0) - ar
     opts_sm2@trXMaxF                       = centre_new(7,ilev,ct,1) + ar
     opts_sm2@trXMinF                       = centre_new(7,ilev,ct,1) - ar

     opts_sm3                               = opts_sm2
     opts_sm3@trYMaxF                       = centre_new(8,ilev,ct,0) + ar
     opts_sm3@trYMinF                       = centre_new(8,ilev,ct,0) - ar
     opts_sm3@trXMaxF                       = centre_new(8,ilev,ct,1) + ar
     opts_sm3@trXMinF                       = centre_new(8,ilev,ct,1) - ar

   ; Resources for pressure on model height levels (filled contours)
     opts_prs                               = opts_spd 
     pmin                                   = toint(min(prs_pl_fine))
     pmax                                   = toint(max(prs_pl_fine))
     delete([/opts_prs@cnLevels, opts_prs@cnFillColors/])
     opts_prs@cnLevels                      = ispan(pmin,pmax,4)
     opts_prs@cnFillColors                  = (/-1,3,4,5,6,7,8,9,10,11,12/)

   ; Resources for Laplacian of horizontal wind speed (line contours)
     opts_lap                               = True 
     opts_lap@cnFillOn                      = False
     opts_lap@cnLineColor                   = "black"
     opts_lap@cnLevelSelectionMode          = "ExplicitLevels"
     opts_lap@cnLevels                      = (/2.0, 5.0/)
     opts_lap@gsnContourNegLineDashPattern  = 16
     opts_lap@gsnContourPosLineDashPattern  = 0
     opts_lap@cnLineThicknessF              = 3.0

     opts_lap@cnInfoLabelOn                 = False
     opts_lap@gsnDraw                       = False
     opts_lap@gsnFrame                      = False
     opts_lap@cnLineLabelsOn                = False
     opts_lap@tiMainString                  = ""
     opts_lap@gsnLeftString                 = ""
     opts_lap@gsnRightString                = ""

   ; Panel plot label resources 
     if (single .eq. 1) then 
      txid                                  = new(1,graphic)
      amid                                  = new(1,graphic)
      panel_strings                         = (/ "a) wind speed"/)
     else
    ; create objects to hold plot captions 
      txid                                  = new(4,graphic)
      amid                                  = new(4,graphic)
    ; plot captions 
;      panel_strings                         = (/ "a) wind speed", \
;                                                 "b) smoothed wind speed", \
;                                                 "c) pressure", \
;                                                 "d) smoothed pressure"/)

      panel_strings                         = (/ "a) wind speed", \
                                                 "b) smth1 wind speed", \
                                                 "c) smth2 wind speed", \
                                                 "d) smth3 wind speed"/)

    ; panel coordinates
      xf = (/0.05, 0.57, 0.05, 0.57/)
      yf = (/0.90, 0.90, 0.52, 0.52/)
      wf = 0.35
      hf = 0.32
     end if 

     txres                                  = True
     txres@txPerimOn                        = True
     txres@txFontHeightF                    = 0.0125
     txres@txBackgroundFillColor            = "White"

     amres                                  = True
     amres@amParallelPosF                   = -0.48
     amres@amOrthogonalPosF                 = -0.48
     amres@amJust                           = "TopLeft"

   ; produce single plot, or a 4-panel plot 

     if (single .eq. 1) then 
    ; plot wind speed and overlay Laplacian 
      plot_spd            = gsn_csm_contour_map(wks, spd_pl_fine(:,:), opts_spd)
    ; overlay figure labels
      txid(0)             = gsn_create_text(wks, panel_strings(0), txres)
      amid(0)             = gsn_add_annotation(plot_spd, txid(0), amres)
     else
    ; panel 1 (wind speed)
      opts_spd@vpXF       = xf(0)
      opts_spd@vpYF       = yf(0)
      opts_spd@vpWidthF   = wf
      opts_spd@vpHeightF  = hf
      plot_spd            = gsn_csm_contour(wks, spd_pl_fine(:,:), opts_spd)
      txid(0)             = gsn_create_text(wks, panel_strings(0), txres)
      amid(0)             = gsn_add_annotation(plot_spd, txid(0), amres)
    ; panel 2 (smoothed wind speed)
      opts_ssm@vpXF       = xf(1)
      opts_ssm@vpYF       = yf(1)
      opts_ssm@vpWidthF   = wf
      opts_ssm@vpHeightF  = hf
      spd_pl_sm           = smth9_Wrap(spd_pl_fine, 0.5, 0.5, True)
      plot_ssm            = gsn_csm_contour(wks, spd_pl_sm(:,:), opts_ssm)
      txid(1)             = gsn_create_text(wks, panel_strings(1), txres)
      amid(1)             = gsn_add_annotation(plot_ssm, txid(1), amres)
    ; panel 3 (geopotential height)
      opts_sm2@vpXF       = xf(2)
      opts_sm2@vpYF       = yf(2)
      opts_sm2@vpWidthF   = wf
      opts_sm2@vpHeightF  = hf
      spd_pl_sm2          = smth9_Wrap(spd_pl_sm, 0.5, 0.5, True)
      plot_prs            = gsn_csm_contour(wks, spd_pl_sm2(:,:), opts_sm2)
      txid(2)             = gsn_create_text(wks, panel_strings(2), txres)
      amid(2)             = gsn_add_annotation(plot_prs, txid(2), amres)
    ; panel 4 (smoothed geopotential height)
      ;opts_psm            = opts_prs
      opts_sm3@vpXF       = xf(3)
      opts_sm3@vpYF       = yf(3)
      opts_sm3@vpWidthF   = wf
      opts_sm3@vpHeightF  = hf
      spd_pl_sm3          = smth9_Wrap(spd_pl_sm2, 0.5, 0.5, True)
      ;prs_pl_sm           = smth9_Wrap(prs_pl_sm2, 0.5, 0.5, True)
      plot_psm            = gsn_csm_contour(wks, spd_pl_sm3(:,:), opts_sm3)
      txid(3)             = gsn_create_text(wks, panel_strings(3), txres)
      amid(3)             = gsn_add_annotation(plot_psm, txid(3), amres)
     end if 

   ; Overlay position of max Laplacian (relative to wind speed field)
     mres                 = True
     m_colours            = (/"black","red4","purple4", "grey45"/)
     m_indices            = (/12, 12, 12, 12/)
     m_sizes              = (/10.0, 10.0, 10.0, 10.0/)
     m_thick              = (/3.0, 3.0, 3.0, 3.0/)
     markers              = new( (/4,dimsizes(m_colours)/),graphic)

   ; Storm centre position from each method 
     centres              = (/centre_spd, centre_s2, centre_psm, centre_ps2/)

   ; Array of panels 
     plots                = (/plot_spd, plot_ssm, plot_prs, plot_psm/)

   ; Overlay all these onto plot 
     do i = 0, dimsizes(m_colours)-1
      mres                = True 
      mres@gsMarkerIndex  = m_indices(i)
      mres@gsMarkerSizeF  = m_sizes(i)
      mres@gsMarkerColor  = m_colours(i)
      mres@gsMarkerThicknessF = m_thick(i)
      markers(0,i)        = gsn_add_polymarker(wks, plot_spd, centres(i,1), centres(i,0), mres)
      markers(1,i)        = gsn_add_polymarker(wks, plot_ssm, centres(i,1), centres(i,0), mres)
      markers(2,i)        = gsn_add_polymarker(wks, plot_prs, centres(i,1), centres(i,0), mres)
      markers(3,i)        = gsn_add_polymarker(wks, plot_psm, centres(i,1), centres(i,0), mres)
     end do 

   ; Tidy up before final plots later in script
     delete([/opts_spd, opts_lap, txid, amid, panel_strings/])

;=================================== 
; Add a legend (only on one panel) 
;===================================

; Display legend; give us control; label font; label height
  lg_opts                            = True
  lg_opts@pmLegendDisplayMode        = "Always"
  lg_opts@lgAutoManage               = False
  lg_opts@lgLabelFont                = "Helvetica"
  lg_opts@lgLabelFontHeightF         = 0.060

; Customise labels and line colours 
  lg_opts@lgItemType                 = "Markers"
  lg_opts@lgMarkerColors             = m_colours
  lg_opts@lgMarkerIndexes            = m_indices
  lg_opts@lgMarkerSizeF              = mres@gsMarkerSizeF
  lg_opts@lgMarkerThicknessF         = mres@gsMarkerThicknessF
  lsize                              = dimsizes(lg_opts@lgMarkerColors)

  lg_opts@lgLabelPosition            = "Right"
  lg_opts@lgItemPlacement            = "ExplicitPlacement"
  lg_opts@lgItemPositions            = fspan(0.05, 0.95, lsize)
  lg_opts@lgLeftMarginF              = 0.01
  lg_opts@lgRightMarginF             = 0.15

  lg_opts@vpWidthF                   = 0.150
  lg_opts@vpHeightF                  = 0.090
  lg_opts@lgPerimColor               = "black"
  lg_opts@lgPerimThicknessF          = 3.0
  lg_opts@lgPerimFill                = "SolidFill"
  lg_opts@lgPerimFillColor           = "white"

  lab_arr                            = (/"spd_min", "spd_smth_min", "prs_smth", "prs_smth2"/)
  lgd                                = gsn_create_legend(wks,lsize,lab_arr,lg_opts)

; Add legend to plot (top left or right hand corner)  
  am_opts                            = True

  if (lgd0 .eq. "tr") then
   am_opts@amJust                    = "TopRight"
   am_opts@amParallelPosF            = 0.5
   am_opts@amOrthogonalPosF          = -0.5
  elseif (lgd0 .eq. "tl") then
   am_opts@amJust                    = "TopLeft"
   am_opts@amParallelPosF            = -0.5
   am_opts@amOrthogonalPosF          = -0.5
  elseif (lgd0 .eq. "br") then
   am_opts@amJust                    = "BottomRight"
   am_opts@amParallelPosF            = 0.5
   am_opts@amOrthogonalPosF          = 0.5
  elseif (lgd0 .eq. "bl") then
   am_opts@amJust                    = "BottomLeft"
   am_opts@amParallelPosF            = -0.5
   am_opts@amOrthogonalPosF          = 0.5
  end if

  annotate_0                         = gsn_add_annotation(plot_spd, lgd, am_opts)

;========================================================
; Create a main title at the top of the page, and plot 
;======================================================== 

  title_res                   = True
  title_str                   = "Panel plot at "+hy+" m AGL: "+ens0+\
                                ", "+dat+" (T+"+it+")"
  title_res@txFontHeightF     = 0.02
  title_res@txFont            = "Helvetica"
  xpos_title                  = 0.50
  ypos_title                  = 0.94

  gsn_text_ndc(wks, title_str, xpos_title, ypos_title, title_res)

; Call 'draw' and 'frame'
  maximize_output(wks, True)

; Tidy up
  delete([/u_plane,v_plane,w_plane,f_plane,vort_plane,avo_plane,spd_plane/])

 end do			    ; End model level loop (do ilev = 0, num_l-1)

 ct  = ct + 1		    ; Counter variable (time)

; Tidy up before next iteration (all variables from 'pc' and 'pd' streams)
  delete([/lon,lat,u,v,vort,avo,f0,f1/])
  delete([/lonK,latK,lon1K,lat1K,str_xz,str_yz,kmh,t,p,w,z,rho,rho_a/])

 end do     ; End time loop (do it = 0, times-1)

end 