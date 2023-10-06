; Produce panel plots of tangential and radial wind for different centre-finding methods

; Run using:

; ncl opt=\"png\" w0=0.5 z0=0 dist=1.4 clr=5 calc=1 rad0=0 nt=121 nr=22 ar=1.05 lgd0=\"tr\"
; cn0=\"slp\" sm=0 sc=5 w_check=1 typ=\"ring\" run=1 ri=12 tc_mbgt5_both_centre.ncl

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
   tf0  = 647
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
; 07/04/20 --> this information feeds into 'setup_cyl' (interpolating to cylindrical grids)
; dx = 0.04 (dist*50); dx = 0.02 (dist*100); dx = 0.01 (dist*200)
  llbox     = toint(dist*200)-1                ; Calculate domain size (grid points)
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
     lat_pd = fspan(w_plane&latitude(0), w_plane&latitude(w_t-1), (w_t*4)-1)
     lon_pd = fspan(w_plane&longitude_1(0), w_plane&longitude_1(w_n-1), (w_n*4)-1)

   ; add metadata 
     lat_pc!0 = "latitude"
     lat_pc@units = "degrees_north"
     lon_pc!0 = "longitude"
     lon_pc@units = "degrees_east"
     lat_pd!0 = "latitude"
     lat_pd@units = "degrees_north"
     lon_pd!0 = "longitude"
     lon_pd@units = "degrees_east"
   
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
     r0            = 0.10

   ; 08/10/2019 --> 'ct' = 53, 54, etc (counter variable doesn't start at 0...)
     print_clock("Calculating storm centre position...")
     centre_arr    = st_centre(ilev, it, cn, r0, centre, vort_pl_fine, prs_pl_fine, spd_pl_fine, \
                               centre_new, dy0, dx0, "comp", ct, "pd", "on", "5min")

   ; Updated storm track information contained in 'centre_new'
     centre_new = centre_arr[0] 

     print(centre(cn,it,0)+", "+centre(cn,it,1))

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
     cen_w    = w_pl_fine( {cen_0n:cen_0x}, {cen_1n:cen_1x} )
     cen_d    = div_pl_fine( {cen_0n:cen_0x}, {cen_1n:cen_1x} )

   ; If average over several grid points is > threshold value, search for centre again
     ave_w    = avg(cen_w)
     ave_d    = avg(cen_d)
     ave_out  = sprintf("%0.1f",ave_w)
     print("ave_w = "+ave_w)

   ; Coordinates of grid points either side (bigger grid than above)
     delete([/cen_0n, cen_0x, cen_1n, cen_1x/])
     cen_0n = cen_loc0 - (2 * dy0)
     cen_0x = cen_loc0 + (2 * dy0)
     cen_1n = cen_loc1 - (2 * dy0)
     cen_1x = cen_loc1 + (2 * dy0)

     if (ilev .gt. z0 .and. ave_w .gt. w0) then

      print("ave. vertical velocity = "+ave_w+" ; ave. divergence = "+ave_d)

    ; Create smaller grid around previous centre  
      w_pl_new  = w_pl_fine({cen_0n:cen_0x},{cen_1n:cen_1x})

    ; Find minimum vertical velocity on this smaller grid  
      w_min_new = min(w_pl_new)

    ; Reshape to 1D array, and find index of minimum vertical velocity 
      dims_w    = dimsizes(w_pl_new)
      w_1d      = ndtooned(w_pl_new)
      inds_w    = ind_resolve(minind(w_1d), dims_w)

    ; Retrieve lat/lon information from subset of grid
    ; 07/04/20 --> changed 'longitude_1' to 'longitude' after creating finescale grid 
      ltN       = w_pl_new&latitude
      lnN       = w_pl_new&longitude

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

   ; temporary solution
     pc_pl_fine = (/u_pl_fine, v_pl_fine, spd_pl_fine, vort_pl_fine, \
                    avo_pl_fine, f_pl_fine, str_pl_fine, div_pl_fine/)

     pd_pl_fine = (/w_pl_fine, prs_pl_fine, prs_pl_fine, prs_pl_fine, \
                    prs_pl_fine, prs_pl_fine, rho_pl_fine, rho_pl_fine/)

   ; Last option = 0 --> basic variables output 
   ; Last option = 1 --> additional MBGT variables output  
     print_clock("Translating to cylindrical grid...")
     cyl_arr_0  = setup_cyl(it, dist, ll_arr, lon_pc, lat_pc, lon_pd, lat_pd, lon_pd, lat_pd, \
                            pc_pl_fine, pd_pl_fine, ntot, centre_new, \
                            5, sm, ar, nr, nt, ilev, num_l, 1, "comp", ct, \
                            0, xpos_all, ypos_all, xcen_all, ycen_all, llbox)

;     cyl_arr_0  = setup_cyl(it, dist, ll_arr, lon, lat, lonK, latK, lon1K, lat1K, \
;                          pc_plane, pd_plane, ntot, centre_new, \
;                          5, sm, ar, nr, nt, ilev, num_l, 1, "comp", ct, \
;                          0, xpos_all, ypos_all, xcen_all, ycen_all)

     xpos_all_0 = cyl_arr_0[0]
     ypos_all_0 = cyl_arr_0[1]
     lat_max_0  = cyl_arr_0[4]
     lon_max_0  = cyl_arr_0[5]
     radii_0    = cyl_arr_0[6]
     rad_size_0 = cyl_arr_0[7]

     vr_int_0   = cyl_arr_0[14]
     vt_int_0   = cyl_arr_0[15]

;==========================================================================
; Repeat for multiple centre-finding methods to output cylindrical grids
;==========================================================================

   ; sc = 5
   ; wind speed minimum (not suitable for the development stage) 
     cyl_arr_1  = setup_cyl(it, dist, ll_arr, lon_pc, lat_pc, lon_pd, lat_pd, lon_pd, lat_pd, \
                          pc_pl_fine, pd_pl_fine, ntot, centre_new, \
                          6, sm, ar, nr, nt, ilev, num_l, 1, "comp", ct, \
                          0, xpos_all, ypos_all, xcen_all, ycen_all, llbox)

;     cyl_arr_1  = setup_cyl(it, dist, ll_arr, lon, lat, lonK, latK, lon1K, lat1K, \
;                          pc_plane, pd_plane, ntot, centre_new, \
;                          6, sm, ar, nr, nt, ilev, num_l, 1, "comp", ct, \
;                          0, xpos_all, ypos_all, xcen_all, ycen_all)

     xpos_all_1 = cyl_arr_1[0]
     ypos_all_1 = cyl_arr_1[1]
     lat_max_1  = cyl_arr_1[4]
     lon_max_1  = cyl_arr_1[5]
     radii_1    = cyl_arr_1[6]
     rad_size_1 = cyl_arr_1[7]

     vr_int_1   = cyl_arr_1[14]
     vt_int_1   = cyl_arr_1[15]

;===================================== 
; START MID-SCRIPT PLOTTING ROUTINE 
;===================================== 

   ; Print out x-y plots of u,v,w at selected heights (2,4,6,10,12,14 km)  
     ;if (ilev .eq. 28) then   

   ; Plot the data as we go  
     output = "$nep/nepartak/images/cyl_coords/cxy_4panel_"+dat+\
              "_"+ens0+"_both_"+typ+"_"+hy0+"m_sc"+sc+"_T"+it
     wks    = gsn_open_wks(opt, output)

     print("ct = "+ct)

   ; Tangential wind resources (sc = 3)
     opts_vt0                               = True
     opts_vt0@cnFillOn                      = True
     opts_vt0@sfXArray                      = xpos_all_0(ilev,:,:,ct)
     opts_vt0@sfYArray                      = ypos_all_0(ilev,:,:,ct)
     opts_vt0@cnLineLabelInterval           = 2.0
     opts_vt0@cnLineLabelFontHeightF        = 0.012
     opts_vt0@cnLineLabelBackgroundColor    = "transparent"
     opts_vt0@cnLineLabelPlacementMode      = "constant"
     opts_vt0@cnLinesOn                     = False
     opts_vt0@cnInfoLabelOn                 = False
     opts_vt0@cnLevelSelectionMode          = "ExplicitLevels"
     opts_vt0@cnFillPalette                 = "radar_new"
     opts_vt0@cnLevels                      = (/5.0, 10.0, 15.0, 20.0,\
                                               25.0, 30.0, 35.0, 40.0, \
                                               45.0, 50.0, 55.0, 60.0, 65.0/)
     opts_vt0@cnFillColors                  = (/-1,0,1,2,3,4,5,6,7,8,9,10,11,12/)
     opts_vt0@gsnMaximize                   = False

   ; No title or additional strings 
     opts_vt0@tiMainString                  = ""
     opts_vt0@gsnLeftString                 = ""
     opts_vt0@gsnRightString                = ""

   ; Turn individual label bar on/off
     opts_vt0@lbLabelBarOn                  = False
     opts_vt0@lbBoxEndCapStyle              = "TriangleBothEnds"
     opts_vt0@lbOrientation                 = "Vertical"

     opts_vt0@trYMaxF                       = lat_max_0 + ar
     opts_vt0@trYMinF                       = lat_max_0 - ar 
     opts_vt0@trXMaxF                       = lon_max_0 + ar 
     opts_vt0@trXMinF                       = lon_max_0 - ar

   ; Do not draw the plot or advance the frame 
     opts_vt0@gsnDraw                       = False
     opts_vt0@gsnFrame                      = False
     opts_vt0@gsnAddCyclic                  = False

   ; Resources for other centre-finding methods 
     opts_vt1                               = opts_vt0
     opts_vt1@sfXArray                      = xpos_all_1(ilev,:,:,ct)
     opts_vt1@sfYArray                      = ypos_all_1(ilev,:,:,ct)
     opts_vt1@trYMaxF                       = lat_max_1 + ar
     opts_vt1@trYMinF                       = lat_max_1 - ar
     opts_vt1@trXMaxF                       = lon_max_1 + ar
     opts_vt1@trXMinF                       = lon_max_1 - ar

   ; Background plot resources 
     opts_xy0                               = True
     opts_xy0@tmYROn                        = False
     opts_xy0@tmYLOn                        = True
     opts_xy0@tmXTOn                        = False
     opts_xy0@tmXBOn                        = True
     opts_xy0@xyLineColor                   = -1
     opts_xy0@gsnDraw                       = False
     opts_xy0@gsnFrame                      = False
     opts_xy0@trYMaxF                       = lat_max_0 + ar
     opts_xy0@trYMinF                       = lat_max_0 - ar
     opts_xy0@trXMaxF                       = lon_max_0 + ar
     opts_xy0@trXMinF                       = lon_max_0 - ar

     opts_xy1                               = opts_xy0
     opts_xy1@trYMaxF                       = lat_max_1 + ar
     opts_xy1@trYMinF                       = lat_max_1 - ar
     opts_xy1@trXMaxF                       = lon_max_1 + ar
     opts_xy1@trXMinF                       = lon_max_1 - ar

   ; Resources for radial wind 
     opts_vr0                               = opts_vt0
     delete([/opts_vr0@cnLevels, opts_vr0@cnFillColors, opts_vr0@cnFillPalette/])
     opts_vr0@cnFillPalette                 = "BlueDarkRed18"
     opts_vr0@cnLevels                      = (/-20., -15., -10., -5., -3., -2., -1., -0.5, \
                                                0, 0.5, 1., 2., 3., 5., 10., 15., 20./)
     opts_vr0@cnFillColors                  = (/0,2,3,4,5,6,7,8,-1,-1,10,11,12,13,14,15,16,17/)

     opts_vr1                               = opts_vt1
     delete([/opts_vr1@cnLevels, opts_vr1@cnFillColors, opts_vr1@cnFillPalette/])
     opts_vr1@cnFillPalette                 = "BlueDarkRed18"
     opts_vr1@cnLevels                      = (/-20., -15., -10., -5., -3., -2., -1., -0.5, \
                                                0, 0.5, 1., 2., 3., 5., 10., 15., 20./)
     opts_vr1@cnFillColors                  = (/0,2,3,4,5,6,7,8,-1,-1,10,11,12,13,14,15,16,17/)

   ; Resources for wind speed on Cartesian grid 
     opts_sp0                               = opts_vt0
     opts_sp0@tmYROn                        = False
     opts_sp0@tmYLOn                        = True
     opts_sp0@tmXTOn                        = False
     opts_sp0@tmXBOn                        = True
     delete([/opts_sp0@sfXArray, opts_sp0@sfYArray/])

     opts_sp1                               = opts_vt1
     opts_sp1@tmYROn                        = False
     opts_sp1@tmYLOn                        = True
     opts_sp1@tmXTOn                        = False
     opts_sp1@tmXBOn                        = True
     delete([/opts_sp1@sfXArray, opts_sp1@sfYArray/])

   ; Panel plot label resources 
     txid                                   = new(6,graphic)
     amid                                   = new(6,graphic)

     panel_strings                          = (/ "a) speed (sc = 5)",\
                                                 "b) vtan (sc = 5)",\
                                                 "c) vrad (sc = 5)",\
                                                 "d) speed (sc = 6)",\
                                                 "e) vtan (sc = 6)",\
                                                 "f) vrad (sc = 6)"/)

     txres                                  = True
     txres@txPerimOn                        = True
     txres@txFontHeightF                    = 0.0125
     txres@txBackgroundFillColor            = "White"

     amres                                  = True
     amres@amParallelPosF                   = -0.48
     amres@amOrthogonalPosF                 = -0.48
     amres@amJust                           = "TopLeft"

   ; Coordinates of panels 
     xf = (/0.05, 0.38, 0.71, 0.05, 0.38, 0.71/)
     yf = (/0.82, 0.82, 0.82, 0.48, 0.48, 0.48/)
     wf = 0.24
     hf = 0.24

   ; Panel 1 (wind speed on Cartesian grid)
     opts_sp0@vpXF        = xf(0)
     opts_sp0@vpYF        = yf(0)
     opts_sp0@vpWidthF    = wf
     opts_sp0@vpHeightF   = hf
     plot_sp0             = gsn_csm_contour(wks, spd_plane(:,:), opts_sp0)
     txid(0)              = gsn_create_text(wks, panel_strings(0), txres)
     amid(0)              = gsn_add_annotation(plot_sp0, txid(0), amres)

   ; Panel 2 (tangential wind --> first centre-finding method)
     opts_vt0@vpXF        = xf(1)
     opts_vt0@vpYF        = yf(1)
     opts_vt0@vpWidthF    = wf
     opts_vt0@vpHeightF   = hf

     opts_xy_0            = opts_xy0
     opts_xy_0@vpXF       = xf(1)
     opts_xy_0@vpYF       = yf(1)
     opts_xy_0@vpWidthF   = wf
     opts_xy_0@vpHeightF  = hf

     plot_vt              = gsn_csm_xy(wks, xpos_all_0(ilev,:,:,ct), \
                                       ypos_all_0(ilev,:,:,ct), opts_xy_0)
     plot_vt0             = gsn_csm_contour(wks, vt_int_0(:,:), opts_vt0)
     overlay(plot_vt, plot_vt0)

   ; Overlay labels
     txid(1)              = gsn_create_text(wks, panel_strings(1), txres)
     amid(1)              = gsn_add_annotation(plot_vt, txid(1), amres)

   ; Panel 3 (radial wind --> first centre-finding method)
     opts_vr0@vpXF        = xf(2)
     opts_vr0@vpYF        = yf(2)
     opts_vr0@vpWidthF    = wf
     opts_vr0@vpHeightF   = hf

     opts_xy_1            = opts_xy0
     opts_xy_1@vpXF       = xf(2)
     opts_xy_1@vpYF       = yf(2)
     opts_xy_1@vpWidthF   = wf
     opts_xy_1@vpHeightF  = hf

     plot_vr              = gsn_csm_xy(wks, xpos_all_0(ilev,:,:,ct), \
                                       ypos_all_0(ilev,:,:,ct), opts_xy_1)
     plot_vr0             = gsn_csm_contour(wks, vr_int_0(:,:), opts_vr0)
     overlay(plot_vr, plot_vr0)

   ; Overlay labels                                                                           
     txid(2)              = gsn_create_text(wks, panel_strings(2), txres)
     amid(2)              = gsn_add_annotation(plot_vr, txid(2), amres)

   ; Panel 4 (wind speed on Cartesian grid)                                                 
     opts_sp1@vpXF        = xf(3)
     opts_sp1@vpYF        = yf(3)
     opts_sp1@vpWidthF    = wf
     opts_sp1@vpHeightF   = hf
     spd_pl_sm            = smth9_Wrap(spd_plane, 0.5, 0.5, True)
     plot_sp1             = gsn_csm_contour(wks, spd_pl_sm(:,:), opts_sp1)
     txid(3)              = gsn_create_text(wks, panel_strings(3), txres)
     amid(3)              = gsn_add_annotation(plot_sp1, txid(3), amres)

   ; Panel 5 (tangential wind --> second centre-finding method)
     opts_vt1@vpXF        = xf(4)
     opts_vt1@vpYF        = yf(4)
     opts_vt1@vpWidthF    = wf
     opts_vt1@vpHeightF   = hf

     opts_xy_2            = opts_xy1
     opts_xy_2@vpXF       = xf(4)
     opts_xy_2@vpYF       = yf(4)
     opts_xy_2@vpWidthF   = wf
     opts_xy_2@vpHeightF  = hf

     plot_vtt             = gsn_csm_xy(wks, xpos_all_1(ilev,:,:,ct), \
                                       ypos_all_1(ilev,:,:,ct), opts_xy_2)
     plot_vtt0            = gsn_csm_contour(wks, vt_int_1(:,:), opts_vt1)
     overlay(plot_vtt, plot_vtt0)

   ; Overlay labels  
     txid(4)              = gsn_create_text(wks, panel_strings(4), txres)
     amid(4)              = gsn_add_annotation(plot_vtt, txid(4), amres)

   ; Panel 6 (radial wind --> second centre-finding method)
     opts_vr1@vpXF        = xf(5)
     opts_vr1@vpYF        = yf(5)
     opts_vr1@vpWidthF    = wf
     opts_vr1@vpHeightF   = hf

     opts_xy_3            = opts_xy1
     opts_xy_3@vpXF       = xf(5)
     opts_xy_3@vpYF       = yf(5)
     opts_xy_3@vpWidthF   = wf
     opts_xy_3@vpHeightF  = hf

     plot_vrr             = gsn_csm_xy(wks, xpos_all_1(ilev,:,:,ct), \
                                       ypos_all_1(ilev,:,:,ct), opts_xy_3)
     plot_vrr0            = gsn_csm_contour(wks, vr_int_1(:,:), opts_vr1)
     overlay(plot_vrr, plot_vrr0)

   ; Overlay labels                                                                   
     txid(5)              = gsn_create_text(wks, panel_strings(5), txres)
     amid(5)              = gsn_add_annotation(plot_vrr, txid(5), amres)

;==================================================
; Overlay markers for each storm centre estimate 
;==================================================

; (0) speed_min, (1) speed_smth_min
  m_colours  = (/"royalblue1", "black"/)
  m_indices  = (/8, 8/)
  m_sizes    = (/10.0, 10.0/)
  m_thick    = (/5.0, 5.0/)
  m_st       = (/5,6/)

  markers    = new( (/6/), graphic)

  mres0                    = True
  mres0@gsMarkerIndex      = m_indices(0)
  mres0@gsMarkerSizeF      = m_sizes(0)
  mres0@gsMarkerColor      = m_colours(0)
  mres0@gsMarkerThicknessF = m_thick(0)
  markers(0)               = gsn_add_polymarker(wks, plot_sp0, centre_new(m_st(0),ilev,ct,1),\
                                                centre_new(m_st(0),ilev,ct,0), mres0)
  markers(1)               = gsn_add_polymarker(wks, plot_vt0, centre_new(m_st(0),ilev,ct,1),\
                                                centre_new(m_st(0),ilev,ct,0), mres0)
  markers(2)               = gsn_add_polymarker(wks, plot_vr0, centre_new(m_st(0),ilev,ct,1),\
                                                centre_new(m_st(0),ilev,ct,0), mres0)

  mres1                    = True
  mres1@gsMarkerIndex      = m_indices(1)
  mres1@gsMarkerSizeF      = m_sizes(1)
  mres1@gsMarkerColor      = m_colours(1)
  mres1@gsMarkerThicknessF = m_thick(1)
  markers(3)               = gsn_add_polymarker(wks, plot_sp1, centre_new(m_st(1),ilev,ct,1),\
                                                centre_new(m_st(1),ilev,ct,0), mres1)
  markers(4)               = gsn_add_polymarker(wks, plot_vtt0, centre_new(m_st(1),ilev,ct,1),\
                                                centre_new(m_st(1),ilev,ct,0), mres1)
  markers(5)               = gsn_add_polymarker(wks, plot_vrr0, centre_new(m_st(1),ilev,ct,1),\
                                                centre_new(m_st(1),ilev,ct,0), mres1)

; Also overlay original centre position (from sea level pressure) for comparison
  mark_slp                 = new(2,graphic)
  slp_res                  = mres1
  slp_res@gsMarkerIndex    = 12
  slp_res@gsMarkerColor    = "red4"
  mark_slp(0)              = gsn_add_polymarker(wks, plot_sp0, centre(cn,it,1), \
                                                centre(cn,it,0), slp_res)
  mark_slp(1)              = gsn_add_polymarker(wks, plot_sp1, centre(cn,it,1),\
                                                centre(cn,it,0), slp_res)

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
;  lg_opts@lgItemType                 = "Markers"
;  lg_opts@lgMarkerColors             = m_colours
;  lg_opts@lgMarkerIndexes            = m_indices
;  lg_opts@lgMarkerSizeF              = mres0@gsMarkerSizeF
;  lg_opts@lgMarkerThicknessF         = mres0@gsMarkerThicknessF
;  lsize                              = dimsizes(lg_opts@lgMarkerColors)

;  lg_opts@lgLabelPosition            = "Right"
;  lg_opts@lgItemPlacement            = "ExplicitPlacement"
;  lg_opts@lgItemPositions            = fspan(0.05, 0.95, lsize)
;  lg_opts@lgLeftMarginF              = 0.01
;  lg_opts@lgRightMarginF             = 0.15

;  lg_opts@vpWidthF                   = 0.150
;  lg_opts@vpHeightF                  = 0.090
;  lg_opts@lgPerimColor               = "black"
;  lg_opts@lgPerimThicknessF          = 3.0
;  lg_opts@lgPerimFill                = "SolidFill"
;  lg_opts@lgPerimFillColor           = "white"

;  lab_arr                            = (/"pres_sm2", "vort_min", "windspeed_min"/)
;  lgd                                = gsn_create_legend(wks,lsize,lab_arr,lg_opts)

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

;  annotate_0                         = gsn_add_annotation(plot_1a, lgd, am_opts)

;============================== 
; Add radial circles to plot
;============================== 

      dtr = 0.017453292519943

    ; Get customizations for circles
      radius   = radii(rad_size)
      xct      = (/lon_max_0,lon_max_0,lon_max_0,lon_max_1,lon_max_1,lon_max_1/)
      yct      = (/lat_max_0,lat_max_0,lat_max_0,lat_max_1,lat_max_1,lat_max_1/)
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
      xc       = (/xct(0) + (radius * xcos), xct(1) + (radius * xcos), xct(2) + (radius * xcos), \
                   xct(3) + (radius * xcos), xct(4) + (radius * xcos), xct(5) + (radius * xcos)/)
      yc       = (/yct(0) + (radius * xsin), yct(1) + (radius * xsin), yct(2) + (radius * xsin), \
                   yct(3) + (radius * xsin), yct(4) + (radius * xsin), yct(5) + (radius * xsin)/)

    ; Resources for outer circle
      lnres    = True
      lnres@gsLineThicknessF  = out_thck
      lnres@gsLineColor       = out_clr
      lnres@gsLineDashPattern = 0

    ; Attach outer circle to each panel plot
      plot_outer0 = gsn_add_polyline(wks, plot_sp0, xc(0,:), yc(0,:), lnres)
      plot_outer1 = gsn_add_polyline(wks, plot_vt0, xc(1,:), yc(1,:), lnres)
      plot_outer2 = gsn_add_polyline(wks, plot_vr0, xc(2,:), yc(2,:), lnres)
      plot_outer3 = gsn_add_polyline(wks, plot_sp1, xc(3,:), yc(3,:), lnres)
      plot_outer4 = gsn_add_polyline(wks, plot_vtt0, xc(4,:), yc(4,:), lnres)
      plot_outer5 = gsn_add_polyline(wks, plot_vrr0, xc(5,:), yc(5,:), lnres)

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
       plot_inner = new( (/6,size/), graphic)

       do r = spacing, (radius-spacing), spacing

       ; Calculate arrays for inner circles
         xc  = (/xct(0) + (r * xcos), xct(1) + (r * xcos), xct(2) + (r * xcos), \
                 xct(3) + (r * xcos), xct(4) + (r * xcos), xct(5) + (r * xcos)/)
         yc  = (/yct(0) + (r * xsin), yct(1) + (r * xsin), yct(2) + (r * xsin), \
                 yct(3) + (r * xsin), yct(4) + (r * xsin), yct(5) + (r * xsin)/)

       ; Resources for inner circles
         delete([/lnres@gsLineThicknessF,lnres@gsLineColor,lnres@gsLineDashPattern/])
         lnres@gsLineColor       = in_clr
         lnres@gsLineThicknessF  = in_thck
         lnres@gsLineDashPattern = dpattrn

       ; Attach inner circles to plot
         plot_inner(0,count) = gsn_add_polyline(wks, plot_sp0, xc(0,:), yc(0,:), lnres)
         plot_inner(1,count) = gsn_add_polyline(wks, plot_vt0, xc(1,:), yc(1,:), lnres)
         plot_inner(2,count) = gsn_add_polyline(wks, plot_vr0, xc(2,:), yc(2,:), lnres)
         plot_inner(3,count) = gsn_add_polyline(wks, plot_sp1, xc(3,:), yc(3,:), lnres)
         plot_inner(4,count) = gsn_add_polyline(wks, plot_vtt0, xc(4,:), yc(4,:), lnres)
         plot_inner(5,count) = gsn_add_polyline(wks, plot_vrr0, xc(5,:), yc(5,:), lnres)

         count = count + 1

       end do

      end if

    ; Tidy up
      delete([/lnres, degrees, xcos, xsin, xc, yc/])

    ; Create a main title at the top of the page
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

    ; Add labelbar at bottom of page
      lopt                        = True 
      lopt@lbAutoManage           = False
      lopt@lbOrientation          = "horizontal"
      lopt@lbLabelAlignment       = "InteriorEdges"
      lopt@lbFillPattern          = "SolidFill"
      lopt@lbMonoFillPattern      = True
      lopt@lbBoxEndCapStyle       = "TriangleBothEnds"      

    ; Tangential wind 
      gsn_define_colormap(wks,"radar_new")
      lopt@lbFillColors           = (/0,2,3,4,5,6,7,8,9,10,11,12,13,14/)
      label_str                   = (/"5.0","10.0","15.0","20.0","25.0","30.0","35.0",\
                                      "40.0","45.0","50.0", "55.0", "60.0", "65.0"/)
      nlabels                     = dimsizes(label_str)
      lopt@lbLabelFontHeightF     = 0.010
      lopt@txFont                 = "Helvetica"
      lopt@lbTitleOn              = True
      lopt@lbTitleString          = "Tangential wind (m s~S~-1~N~)"
      lopt@lbTitleFontHeightF     = 0.0125
      lopt@lbTitlePosition        = "Bottom"
      lopt@vpWidthF               = 0.40
      lopt@vpHeightF              = 0.08 
      xpos_lbar                   = 0.08
      ypos_lbar                   = 0.15

    ; Use routine to overlay labelbar (tangential wind)
      gsn_labelbar_ndc(wks, nlabels, label_str, xpos_lbar, ypos_lbar, lopt)

    ; Also overlay labelbar for radial wind 
      lopt2                       = lopt
      gsn_define_colormap(wks,"BlueDarkRed18")
      delete([/lopt2@lbFillColors,lopt2@vpWidthF,lopt2@vpHeightF/])
      lopt2@lbFillColors          = (/3,4,5,6,7,8,9,10,0,0,12,13,14,15,16,17,18,19/)
      label_str2                  = (/"-20.0","-15.0","-10.0","-5.0","-3.0","-2.0","-1.0",\
                                      "-0.5","0.0","0.5", "1.0", "2.0", "3.0", \
                                      "5.0", "10.0", "15.0", "20.0", "25.0"/)
      nlabels2                    = dimsizes(label_str2)
      lopt2@vpWidthF              = 0.40
      lopt2@vpHeightF             = 0.08
      lopt2@lbTitleString         = "Radial wind (m s~S~-1~N~)"
      xpos_lbar2                  = 0.53
      ypos_lbar2                  = 0.15

    ; Use routine to overlay labelbar (radial wind)
      gsn_labelbar_ndc(wks, nlabels2, label_str2, xpos_lbar2, ypos_lbar2, lopt2)

    ; Draw with everything overlaid
      maximize_output(wks,True)

  ; Tidy up
    delete([/vt_int_0,vt_int_1,vr_int_0,vr_int_1/])
    delete([/u_plane,v_plane,w_plane,f_plane,vort_plane,avo_plane/])

  ; Tidy up before final plots later in script                                              
    delete([/opts_vt0, opts_vt1, opts_vr0, opts_vr1, opts_xy0, opts_xy1, opts_xy_0, \
             opts_xy_1, opts_xy_2, opts_xy_3, txid, amid, panel_strings/])

   end do			    ; End model level loop (do ilev = 0, num_l-1)

  ct  = ct + 1		    ; Counter variable (time)

; Tidy up before next iteration (all variables from 'pc' and 'pd' streams)
  delete([/lon,lat,u,v,vort,avo,f0,f1/])
  delete([/lonK,latK,lon1K,lat1K,str_xz,str_yz,kmh,t,p,w,z,rho,rho_a/])

 end do     ; End time loop (do it = 0, times-1)

end 