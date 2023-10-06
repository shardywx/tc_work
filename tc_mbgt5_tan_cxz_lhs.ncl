; Script to calculate the radial momentum budget (LHS) following Montgomery et al. (2018)

; Run using:

; ncl opt=\"x11\" w0=0.5 z0=0 dist=1.4 clr=5 calc=1 rad0=0 nt=73 nr=22 ar=1.05 lgd0=\"tr\" plt=0
; new0=0 mean0=0 cn0=\"slp\" sm=0 sc=5 typ=\"ring\" thr=2 max0=0 grp=\"8\" 
; run=1 ri=12 tc_mbgt5_tan_cxz_lhs.ncl

; 'opt'   = output file format ("pdf" or "x11")
; 'dist'  = size of box following storm (degrees)
; 'clr'	  = colour map for plots (option 4 is used in Roger Smith's papers):
; 	    "ncl_default" [1],  "bdr_extra" [2], "amwg256" [3], "GMT_polar" [4]
; 'calc'  = centre calculated on each level (0), calculated offline (1) or using 980 m vort (2)
; 'rad0'  = radius of final plots: 150 km (0) or 200 km (1)
; 'typ'	  = VC phase ("sym","asym","stoa","atos")
; 'thr'   = if plotting "v" or "slp" data, threshold above/below which we're plotting (e.g. 2/-2)
; 'grp'   = which set of times to analyse: "6", "5", "4c", "4b", "4a", "4"
; 'ar'    = distance in degrees of outer radius (1.0, 2.0, etc)
; 'nr'    = number of radial circles between r = 0 and r = 'ar' (21, 41, 61, etc)
; 'nt'    = number of azimuth angles in cylindrical grid
; 'new0'  = reduced scale for plotting tangential wind and AAM (1)
; 'mean0' = plot eddy and mean terms (1) instead of diffusive tendencies (0)
; 'cn0'   = offline storm centre: "slp", "vort", "geo", "geo_sm"
; 'nts'   = start time for no-VC composite 
; 'ntf'   = finish time for no-VC composite

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
   if (max0 .eq. 1) then 
    ts0  = 645
    tf0  = 672
   else
    ts0  = 646
    tf0  = 671 
   end if 
; Ring-like to monopole transition (T+58 to T+60)                                            
  elseif (typ .eq. "r2m") then
   if (max0 .eq. 1) then 
    ts0  = 681
    tf0  = 708
   else
    ts0  = 682
    tf0  = 707
   end if 
; Monopole phase (T+63 to T+66)                                                              
  elseif (typ .eq. "mono") then
   if (max0 .eq. 1) then 
    ts0  = 741
    tf0  = 780
   else
    ts0  = 742
    tf0  = 779
   end if 
; Monopole to ring-like transition (T+71 to T+73)                                            
  elseif (typ .eq. "m2r") then
   if (max0 .eq. 1) then 
    ts0  = 837
    tf0  = 864
   else
    ts0  = 838
    tf0  = 863
   end if 
; For now, this covers the period from T+48 to T+84
  elseif (typ .eq. "all") then 
    ts0  = 563
    tf0  = 995
; Troubleshooting 
  else
    ts0  = 647
    tf0  = 649
  end if

; Starting minus ending time indices of all VC phases                                        
  ntot  = (tf0 - ts0) + 1

; Additional times on each end, to calculate tendencies
  tx    = tf0 + 1
  tn    = ts0 - 1
  td    = (tx - tn) + 1

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

  time_str  = new(numINPUT_c,string)
  time_arr  = new(numINPUT_c,string)
  title_arr = new(numINPUT_c,string)

  time_str  = hr_min + " UTC " + sprinti("%0.2i ", day) + month_abbr(month)
  time_arr  = sprinti("%0.2i", day) + month_abbr(month) + "_" + hr_min + "Z"
  title_arr = "Valid at "+time_str


;========================================                                                       
; Define cylindrical coordinate arrays                                                          
;========================================                                                       

  pii        = 3.14159265
  radii      = fspan(0.0,ar,nr)      ; Array of radial distances ('Start','End','Number')       
  thetas     = fspan(0.0,360.0,nt)   ; Array of azimuth angles (0, 5, 10, etc)                  
  thetas_rad = (pii/180)*thetas      ; Degrees to radians                                       
  rad_size   = dimsizes(radii)-1     ; Index for later calculation                              

; Define constants                                                                              
  rd         = 287.0                 ; Gas constant for dry air (J/kg/K^2)                      
  cp0        = 1004.0                ; Specific heat of dry air at constant pressure (J/kg/K)   
  kp         = 0.286                 ; For Exner function calculations (0.286)                  
  pref       = 1000.0                ; Reference pressure p0 (1000 hPa)  

;============================================================================
; Now calculate tangential wind tendency using centred finite differencing
;============================================================================

; Tangential wind 
  mean_vtan_all       = new( (/nlev, dimsizes(radii), ntot/), "float")
  mean_vtan_all!0     = "lev"
  mean_vtan_all!1     = "rad"
  mean_vtan_all!2     = "time"
  mean_vtan_all@description = "Azimuthally averaged tangential wind"
  mean_vtan_all@units = "m s~S~-1~N~"

; Radial wind 
  mean_vrad_all       = new( (/nlev, dimsizes(radii), ntot/), "float")
  mean_vrad_all!0     = "lev"
  mean_vrad_all!1     = "rad"
  mean_vrad_all!2     = "time"
  mean_vrad_all@description = "Azimuthally averaged radial wind"
  mean_vrad_all@units = "m s~S~-1~N~"

; Coriolis parameter
  mean_cor_all        = new( (/nlev, dimsizes(radii), ntot/), "float")
  mean_cor_all!0      = "lev"
  mean_cor_all!1      = "rad"
  mean_cor_all!2      = "time"
  mean_cor_all@description = "Azimuthally averaged Coriolis parameter"
  mean_cor_all@units  = "s~S~-1~N~"

; Total wind 
  mean_vtot_all       = new( (/nlev, dimsizes(radii), ntot/), "float")
  mean_vtot_all!0     = "lev"
  mean_vtot_all!1     = "rad"
  mean_vtot_all!2     = "time"
  mean_vtot_all@description = "Azimuthally averaged total wind"
  mean_vtot_all@units = "m s~S~-1~N~"

  if (max0 .eq. 0) then 
; Tangential wind tendency
   mean_vtan_ten      = new( (/nlev, dimsizes(radii), ntot/), "float")
   mean_vtan_ten!0    = "lev"
   mean_vtan_ten!1    = "rad"
   mean_vtan_ten!2    = "time"
   mean_vtan_ten@description = "Azimuthally averaged tangential wind tendency"
   mean_vtan_ten@units = "m s~S~-1~N~"
  else
   mean_vtan_ten      = new( (/ntot/), "float")
   mean_vtan_ten!0    = "time"
   mean_vtan_ten@description = "Tendency of the azimuthally-averaged maximum tangential wind"
   mean_vtan_ten@units = "m s~S~-1~N~"
  end if 

; Maximum tangential wind (1D)
  mean_vtan_max       = new( (/ntot/), "float")
  mean_vtan_max!0     = "time"
  mean_vtan_max@description = "Azimuthally-averaged maximum tangential wind"
  mean_vtan_max@units = "m s~S~-1~N~"

; Radial wind at vmax
  vrad_vmax           = new( (/ntot/), "float")
  vrad_vmax!0         = "time"
  vrad_vmax@description = "Radial wind at vmax"
  vrad_vmax@units     = "m s~S~-1~N~"

; Coriolis parameter at vmax (for AAM calculation)
  cor_vmax            = new( (/ntot/), "float")
  cor_vmax!0          = "time"
  cor_vmax@description = "Azimuthally-averaged Coriolis parameter"
  cor_vmax@units      = "s~S~-1~N~"


; Choose which storm track to read in  
  if (cn0 .eq. "slp") then
   cn = 0
  elseif (cn0 .eq. "slpf") then
   cn = 1
  end if 

; Create array to hold storm track positions (L727)
  centre_new = new( (/10, nlev, ntot, 2/), "float")

; Arrays to hold all 'xpos' and 'ypos' position arrays                                       
  xpos_all  = new((/nlev,dimsizes(radii),dimsizes(thetas),ntot/),float)
  ypos_all  = new((/nlev,dimsizes(radii),dimsizes(thetas),ntot/),float)

; Arrays to hold lat/lon at centre of cylindrical grid                                       
  xcen_all  = new((/nlev,dimsizes(radii),dimsizes(thetas),ntot/),float)
  ycen_all  = new((/nlev,dimsizes(radii),dimsizes(thetas),ntot/),float)

; Initialise time-varying counter variable
  ct = 0

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
    hyb_plot  = hybC / 1000
    hyb_plot@units = "km"

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

;==========================================================                                    
; Calculate divergence using one-sided difference scheme                                       
;==========================================================                                    

    div = uv2dv_cfd(u, v, lat, lon, 2)
    copy_VarCoords(u, div)
    div@units = "s~S~-1~N~"
    div@name = "Divergence on pressure levels"

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

  ; Convert units of pressure and geopotential height before continuing                    
    p      = p / 100
    z      = z / 10

  ; Add metadata                                                                               
    p@units = "hPa"
    z@units = "dam"

;==========================                                                                   
; Loop over model levels                                                                      
;==========================                                                                   

  ; Loop over chosen levels (L907)                                                            
    do ilev = 0, nlev-1

   ; Create output strings (model level info)                                                 
     hy  = sprintf("%0.0f",hybD(ilev) )
     hy1 = toint(hy)
     hy0 = sprintf("%05g",hy1)

     hy2 = sprintf("%0.0f",hybD(z0) )
     hy3 = tofloat(hy2) / 1000
     print_clock("Working on model level " + ilev+ " ("+hy+" m AGL)")

   ; Read in variables on model levels (L973)                                                 
     u_plane                = u(ilev,:,:)          ; Zonal wind                               
     v_plane                = v(ilev,:,:)          ; Meridional wind                          
     vort_plane             = vort(ilev,:,:)       ; Relative vorticity                       
     w_plane                = w(ilev,:,:)          ; Vertical velocity
     prs_plane              = p(ilev,:,:)          ; Pressure 

     avo_plane              = avo(ilev,:,:)        ; Absolute vorticity
     f_plane                = f0(ilev,:,:)         ; Coriolis parameter

   ; Calculate vector windspeed and smooth                                                    
     spd_plane              = sqrt( (u_plane ^ 2) + (v_plane ^ 2) )
     copy_VarCoords(u_plane, spd_plane)
     vort_pl_smth           = smth9_Wrap(vort_plane, 0.5, 0.5, True)

   ; Combine all arrays above into single, larger array for input into 'setup_cyl'            

   ; 'pc' stream                                                                              
   ; 07/10/2019 --> second line arrays are placeholders                                       
     pc_plane               = (/u_plane, v_plane, spd_plane, vort_plane, \
                                avo_plane, f_plane, vort_plane, vort_plane/)

   ; 'pd' stream                                                                              
   ; 08/10/2019 --> after first two sub-arrays, whole array is a placeholder
     pd_plane               = (/w_plane, prs_plane, prs_plane, prs_plane, prs_plane, \
                                w_plane, w_plane, w_plane/)

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
     lat_pc = fspan(u_plane&latitude(0), u_plane&latitude(u_t-1), (u_t*fc)-1)
     lon_pc = fspan(u_plane&longitude(0), u_plane&longitude(u_n-1), (u_n*fc)-1)

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
     f_pl_fine    = linint2_Wrap(f_plane&longitude, f_plane&latitude, f_plane, False, \
                                 lon_pc, lat_pc, 0)

   ; also interpolate the 'pd' stream variables onto the same grid
     w_pl_fine    = linint2_Wrap(w_plane&longitude_1, w_plane&latitude, w_plane, False, \
                                 lon_pc, lat_pc, 0)
     prs_pl_fine  = linint2_Wrap(prs_plane&longitude, prs_plane&latitude, prs_plane, False, \
                                 lon_pc, lat_pc, 0)

   ; new grid spacing
     dy0          = abs(lat_pc(1) - lat_pc(0))
     dx0          = abs(lon_pc(1) - lon_pc(0))

;==========================================                                                   
; Calculate storm-centre position (L996)                                                      
;==========================================                                                   

   ; Set radial distance (º) when looking for speed/vort min (in ringlike phase)              
     r0            = 0.10

     print_clock("Calculating storm centre position...")
;     centre_arr    = st_centre(ilev, it, cn, r0, centre, vort_plane, prs_plane, spd_plane, \
;                               centre_new, dy0, dx0, "comp", ct, "pd", "off", "5min")

     centre_arr    = st_centre(ilev, it, cn, r0, centre, vort_pl_fine, prs_pl_fine, spd_pl_fine, \
                               centre_new, dy0, dx0, "comp", ct, "pd", "on", "5min")

   ; Updated storm track information contained in 'centre_new'
   ; [5] wind speed min; [6] smoothed wind speed min; [8] Laplacian of [5]; [9] Laplacian of [6]
     centre_new  = centre_arr[0]

;===============================================================================              
; Make sure that the storm 'centre' is not within the eyewall updraft (L1006)                 
;===============================================================================              

   if (w_check .eq. 1) then 

   ; Retrieve coordinates of storm centre                                                     
     cen_loc0 = centre_new(sc,ilev,ct,0)
     cen_loc1 = centre_new(sc,ilev,ct,1)

   ; Coordinates of grid points either side                                                   
     cen_0n = cen_loc0 - dy0
     cen_0x = cen_loc0 + dy0
     cen_1n = cen_loc1 - dy0
     cen_1x = cen_loc1 + dy0

   ; Calculate vertical velocity at these grid points                                         
     cen_w    = w_pl_fine( {cen_0n:cen_0x}, {cen_1n:cen_1x} )

   ; If average over several grid points is above a threshold value, continue searching for centre
     ave_w    = avg(cen_w)
     ave_out  = sprintf("%0.1f",ave_w)

   ; Coordinates of grid points either side (bigger grid than above)                               
     delete([/cen_0n, cen_0x, cen_1n, cen_1x/])
     cen_0n = cen_loc0 - (2 * dy0)
     cen_0x = cen_loc0 + (2 * dy0)
     cen_1n = cen_loc1 - (2 * dy0)
     cen_1x = cen_loc1 + (2 * dy0)

     if (ave_w .gt. w0) then

      print("ave. vertical velocity = "+ave_w)

    ; Create smaller grid around previous centre                                                   
      w_pl_new  = w_pl_fine({cen_0n:cen_0x},{cen_1n:cen_1x})

    ; Find minimum vertical velocity on this smaller grid                                          
      w_min_new = min(w_pl_new)

    ; Reshape to 1D array, and find index of minimum vertical velocity                             
      dims_w    = dimsizes(w_pl_new)
      w_1d      = ndtooned(w_pl_new)
      inds_w    = ind_resolve(minind(w_1d), dims_w)

    ; Retrieve lat/lon information from subset of grid                                             
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
     delete([/cen_w, cen_0n, cen_0x, cen_1n, cen_1x/])

   end if 

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
      if ( abs(d_lat) .gt. 0.02 .or. abs(d_lon) .gt. 0.02) then
       print("Change in storm centre position too large - using position from level below...")
       centre_new(sc,ilev,ct,0) = centre_new(sc,ilev-1,ct,0)
       centre_new(sc,ilev,ct,1) = centre_new(sc,ilev-1,ct,1)
      end if
     end if

;=============================================                                           
; Switch to cylindrical coordinates (L1077)                                              
;=============================================                                           

   ; temporary solution
     pc_pl_fine             = (/u_pl_fine, v_pl_fine, spd_pl_fine, vort_pl_fine, \
                                avo_pl_fine, f_pl_fine, vort_pl_fine, vort_pl_fine/)

     pd_pl_fine             = (/w_pl_fine, prs_pl_fine, prs_pl_fine, prs_pl_fine, prs_pl_fine, \
                                w_pl_fine, w_pl_fine, w_pl_fine/)

   ; Call user-defined function 'setup_cyl'
     print_clock("Translating to cylindrical grid...")
     llbox = toint( dist * (fc*50) )-1
;     cyl_arr  = setup_cyl(it, dist, ll_arr, lon, lat, lonK, latK, lon1K, lat1K, \
;                          pc_plane, pc_plane, ntot, centre_new, \
;                          sc, sm, ar, nr, nt, ilev, nlev, 1, "comp", ct, \
;                          0, xpos_all, ypos_all, xcen_all, ycen_all)

     cyl_arr  = setup_cyl(it, dist, ll_arr, lon_pc, lat_pc, lon_pc, lat_pc, lon_pc, lat_pc, \
                          pc_pl_fine, pd_pl_fine, ntot, centre_new, \
                          sc, sm, ar, nr, nt, ilev, num_l, 1, "comp", ct, \
                          0, xpos_all, ypos_all, xcen_all, ycen_all, llbox)

   ; Output the tangential wind on a cylindrical grid 
     vt_int   = cyl_arr[15]
     spd_int  = cyl_arr[10]
     vr_int   = cyl_arr[14]
     f_int    = cyl_arr[13]

   ; Calculate azimuthally-averaged tangential wind                                      
     do irad  = 0, dimsizes(radii)-1
      mean_vtan_all(ilev,irad,ct) = avg(vt_int(irad,:) )
      mean_vtot_all(ilev,irad,ct) = avg(spd_int(irad,:) )
      mean_vrad_all(ilev,irad,ct) = avg(vr_int(irad,:) )
     end do

   ; EDIT --> possibly add 'delete' command 

  ; End loop over model levels (do ilev = 0, nlev-1)                                    
    end do

  ; Increase values of time-varying counter variables 
    ct  = ct  + 1 

  ; EDIT --> possibly add 'delete' command 

 ; End loop over times (do it = ts0(st)-1, tf0(st)+1)
   end do

;==================================================
; Calculate maximum tangential wind at each time 
;==================================================

 ; Array to hold indices of maximum mean tangential wind 
   vtan_all_ind = new( (/2, ntot/), "integer")

 ; Loop over all required times 
   do ct = 0, ntot-1

   ; Reshape subset of 'mean_vtan_all' into 1D array
     vtan_1d = ndtooned(mean_vtan_all(:,:,ct) )

   ; Read maximum value from 1D array into 'mean_vtan_max'
     mean_vtan_max(ct) = max(vtan_1d)

   ; Retrieve the indices of the maximum value at each time 
     vtan_ind_temp     = ind_resolve(maxind(vtan_1d), dimsizes(mean_vtan_all(:,:,ct) ) )

   ; Index of height of maximum value  
     vtan_all_ind(0,ct) = vtan_ind_temp(0,0)

   ; Index of radius of maximum value
     vtan_all_ind(1,ct) = vtan_ind_temp(0,1)

   ; Retrieve values of 'vrad' and 'f0' at the location of max tangential wind 
     cor_vmax(ct)       = (/mean_cor_all(vtan_ind_temp(0,0), vtan_ind_temp(0,1), ct)/)
     vrad_vmax(ct)      = (/mean_vrad_all(vtan_ind_temp(0,0), vtan_ind_temp(0,1), ct)/)

 ; End loop over times 
   end do

;============================================================
; Calculate tendency of maximum tangential wind (and exit)
;============================================================

 ; Exit the script if we're calculating the maximum tangential wind data 
   if (max0 .eq. 1) then 

  ; Loop over time indices (not including first and last)                                    
    do ct = 1, ntot-2
     cp1  = ct+1
     cm1  = ct-1
     mean_vtan_ten(ct) = (mean_vtan_max(cp1) - mean_vtan_max(cm1) ) / 2
    end do
  ; End loop over times (do ct = 0, ntot-1)

  ; Write all variables out to text files 
    diri_out = "$ar/text/vtan5_max_"+typ+"_group"+grp+"_sc"+sc+"_w"+w0+"_"
    asciiwrite(diri_out+"hgt.txt", vtan_all_ind(0,:) )
    asciiwrite(diri_out+"rad.txt", vtan_all_ind(1,:) )
    asciiwrite(diri_out+"max.txt", mean_vtan_max)
    asciiwrite(diri_out+"ten.txt", mean_vtan_ten)
    asciiwrite(diri_out+"vrad.txt", vrad_vmax)
    asciiwrite(diri_out+"f0.txt", cor_vmax)

    exit()

   end if 

;=============================================================================
; Calculate tangential wind tendency (LHS) using centred finite differences
;=============================================================================

 ; Loop over times (not including first and last)
   do it = 1, ntot-2

    cp1 = it+1
    cm1 = it-1

  ; Multiply by 6 (divide by 2, multiply by 12)
    mean_vtan_ten(:,:,it) = (mean_vtan_all(:,:,cp1) - mean_vtan_all(:,:,cm1) ) * 6 

 ; End loop over time indices (do it = 1, ntot(st)-2)
   end do

;====================================================================================
; Read all values of 'mean_vtan_ten' (not including missing values) into new array  
;====================================================================================

; Arrays to hold maximum d(vtan)/dt and maximum vtan
  vtan_tend_ind  = new ( (/2, td/), "integer")
  vtan_ind       = new ( (/2, td/), "integer")

; Initialise counter variable to calculate number of non-missing values 
  mt = 0 

  print_clock("Working on maximum tangential wind tendency")
  do tt = 0, ntot-1

 ; Reshape to 1D arrays for use with 'ind_resolve'
   vtan_tend_1d        = ndtooned(mean_vtan_ten(:,:,tt) )

 ; Retrieve the indices of the maximum values at each time 
   ten_ind_temp        = ind_resolve(maxind(vtan_tend_1d), dimsizes(mean_vtan_ten(:,:,tt) ) )

 ; Index of height of maximum value 
   z_ind_ten           = ten_ind_temp(0,0)

 ; Index of radius of maximum value 
   r_ind_ten           = ten_ind_temp(0,1)

 ; Read values into larger array 
   vtan_tend_ind(0,tt) = z_ind_ten
   vtan_tend_ind(1,tt) = r_ind_ten

 ; Count number of non-missing values 
   if (.not. ismissing(z_ind_ten) ) then 
    mt = mt + 1 
   end if 
   
  end do 

; Create new array without missing values                                                   
  print("Array contains "+mt+" pairs of non-missing values.")
  vtan_tend_out  = new( (/2, mt/), "integer")

  mt = 0
  do tt = 0, ntot-1
   if (.not. ismissing(vtan_tend_ind(0,tt) ) ) then
    vtan_tend_out(0,mt) = vtan_tend_ind(0,tt)
    vtan_tend_out(1,mt) = vtan_tend_ind(1,tt)
    mt = mt + 1
   end if
  end do

;==========================================
; Now repeat for maximum tangential wind 
;==========================================

; Re-initialise counter variable to calculate number of non-missing values 
  mt = 0

  print_clock("Working on maximum tangential wind")
  do tt = 0, ntot-1

 ; Reshape to 1D arrays for use with 'ind_resolve'
   vtan_1d             = ndtooned(mean_vtan_all(:,:,tt) )

 ; Retrieve the indices of the maximum values at each time
   vtan_ind_temp       = ind_resolve(maxind(vtan_1d), dimsizes(mean_vtan_all(:,:,tt) ) )

 ; Index of height of maximum value
   z_ind_vtan          = vtan_ind_temp(0,0)

 ; Index of radius of maximum value
   r_ind_vtan          = vtan_ind_temp(0,1)

 ; Read values into larger array
   vtan_ind(0,tt)      = z_ind_vtan
   vtan_ind(1,tt)      = r_ind_vtan

 ; Count number of non-missing values
   if (.not. ismissing(z_ind_vtan) ) then
    mt = mt + 1
   end if

  end do

; Create new array without missing values 
  print("Array contains "+mt+" pairs of non-missing values.")
  vtan_all_out   = new( (/2, mt/), "integer")

  mt = 0 
  do tt = 0, ntot-1
   if (.not. ismissing(vtan_ind(0,tt) ) ) then 
    vtan_all_out(0,mt) = vtan_ind(0,tt) 
    vtan_all_out(1,mt) = vtan_ind(1,tt)
    mt = mt + 1  
   end if 
  end do 

;=================================================================
; Write out values to text files for overlaying onto MBGT plots 
;=================================================================

  if (sct .ne. 0) then 
   diri_ten  = "$ar/text/5min_vtan_ten_scatter_"+typ+"_group"+grp+"_sc"+sc+"_w"+w0+"_"
   diri_vtan = "$ar/text/5min_vtan_all_scatter_"+typ+"_group"+grp+"_sc"+sc+"_w"+w0+"_"
   asciiwrite(diri_ten+"hgt.txt", vtan_tend_out(0,:) )
   asciiwrite(diri_ten+"rad.txt", vtan_tend_out(1,:) )
   asciiwrite(diri_vtan+"hgt.txt", vtan_all_out(0,:) )
   asciiwrite(diri_vtan+"rad.txt", vtan_all_out(1,:) )
  end if 

;========================================================================
; Write values to output matrices (text files) using in-built function
;========================================================================

; Additional metadata for 'write_matrix' procedure [21 * 59]
  ncol           = nr-1
  nrow           = nlev
  fmt            = nrow + "f10.3"

  print("Number of columns = "+ncol+" ; number of rows = "+nrow)

; Resources for writing out 
  opts_mat       = True
  opts_mat@row   = True
  opts_mat@title = ""

;================================================================
; Either average over all times, or output values at each time 
;================================================================ 

; Either write out data at each time 
  if (each .eq. 1) then

 ; Set counter variable (time)
   ct = 0

   do it = ts0, tf0

    print("Working on time: T+"+it)

  ; Extract data for a single time 
    lhs1_c = mean_vtan_ten(:,:,ct)

  ; Call procedure and write to text files 
    print_clock("Writing 2D array to matrix for further analysis...")
 
  ; Set file path
    t_dir = "./5min_"
    ;t_out = "lhs1_"+rstr+"_sc"+sc+"_group"+grp+"_"+w0+"w_"+hy3+"km_"+typ+"_T"+it+".txt"
    t_out = "lhs1_"+rstr+"_sc"+sc+"_"+w0+"w_"+hy3+"km_"+typ+"_T"+it+".txt"

    t_path = t_dir + t_out 
    print_clock("Writing to file: "+t_path)
    opts_mat@fout = t_path

  ; Check for existing files before writing out 
    system("/bin/rm -f " + t_path)
    write_matrix( lhs1_c(:,0:nr-2), fmt, opts_mat)

    ct = ct + 1 
  ; Advance counter variable (time)

   end do 
 ; End loop over times

 ; Tidy up text files by deleting indices in first column
   system("sed -i -r 's/.{10}//' ./5min_lhs1*.txt")

 ; Move text files to correct directory 
   system("mv ./5min_lhs1*.txt $ar/text/")

 ; Exit script 
   exit()

; OR, continue as before and write out time-averaged data to text files (01/04/20)
  end if 

; Before plotting, average over all chosen time periods
  lhs1_c         = dim_avg_n_Wrap(mean_vtan_ten, 2)
  lhs1_c@name    = "Local tangential wind tendency"

; Add metadata before plotting
  lhs1_c!0       = "lev"
  lhs1_c&lev     = hyb_plot
  lhs1_c!1       = "rad"
  lhs1_c&rad     = radii

; Name all variables we want to write 
  opt_arr = (/"lhs1"/)

; Array of all composite terms
  mbgt_arr = (/lhs1_c/)

; Call procedure and write time-averaged
  print_clock("Writing 2D array to matrices for further analysis...")
  do mb = 0, dimsizes(opt_arr)-1

    print("Working on "+opt_arr(mb) )

  ; Set file path, etc
    t_dir  = "./5min_"
    if (run .eq. 1) then 
     t_out  = opt_arr(mb)+"_"+rstr+"_sc"+sc+"_group"+grp+"_"+w0+"w_"+hy3+"km_"+typ+".txt" 
    else
     t_out  = opt_arr(mb)+"_sc"+sc+"_group"+grp+"_"+w0+"w_"+hy3+"km_"+typ+".txt"
    end if 
    t_path = t_dir + t_out
    print_clock("Writing to file: "+t_path)
    opts_mat@fout = t_path

  ; Check for existing files and then write 
    system("/bin/rm -f " + t_path)
    write_matrix( mbgt_arr(:,0:nr-2), fmt, opts_mat)

  end do

; Exit the script, or continue plotting                                                    
  if (plt .ne. 1) then

   ; Tidy up text files by deleting indices in first column ('0,1,2,3,...')                
     system("sed -i -r 's/.{10}//' ./5min_lhs1*"+typ+"*.txt")

   ; Move text files to correct directory                                                  
     system("mv ./*"+typ+"*.txt $ar/text/")

   ; Now, exit script                                                                      
     exit()

; Produce a panel plot of key variables (plt=1)...
  else

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
    opts_mbgt@cnlevelSelectionMode         = "ExplicitLevels"

    if (clr .eq. 1) then					       
     opts_mbgt@cnFillPalette		   = "ncl_default"
     opts_mbgt@cnFillColors                = (/2,8,16,28,40,52,64,76,88,100,118,-1,-1,\
                                               134,150,166,178,186,194,\
                                               202,214,226,238,250/)
     opts_mbgt@cnLevels                    = (/-150., -100., -50., -20., -10., -5., \
                                                -3., -2., -1., -0.5, -0.2, \
                                                0, 0.2, 0.5, 1., 2., 3., \
                                                5., 10., 20., 50., 100., 150./)

    elseif (clr .eq. 2) then 
     opts_mbgt@cnFillPalette               = "bdr_extra"
     opts_mbgt@cnFillColors                = (/0,1,2,3,4,5,6,7,8,9,10,-1,-1,\
                                               13,14,15,16,17,18,19,20,21,22,23/)
     opts_mbgt@cnLevels                    = (/-150., -100., -50., -20., -10., -5., \
                                                -3., -2., -1., -0.5, -0.2, \
                                                0, 0.2, 0.5, 1., 2., 3., \
                                                5., 10., 20., 50., 100., 150./)

    elseif (clr .eq. 3) then
     opts_mbgt@cnFillPalette               = "amwg256"
     opts_mbgt@cnFillColors                = (/0,6,12,18,27,36,45,54,66,78,90,-1,-1,\
                                               148,156,166,178,186,194,\
                                               202,214,226,238,250/)
     opts_mbgt@cnLevels                    = (/-150., -100., -50., -20., -10., -5., \
                                                -3., -2., -1., -0.5, -0.2, \
                                                0, 0.2, 0.5, 1., 2., 3., \
                                                5., 10., 20., 50., 100., 150./)
    elseif (clr .eq. 4) then
     opts_mbgt@cnFillPalette               = "GMT_polar"
     opts_mbgt@cnFillColors               = (/3,5,7,-1,11,14,18/)
     opts_mbgt@cnLevels                   = (/-20,-5,0,0.5,5,20/)
    elseif (clr .eq. 5) then 
     opts_mbgt@cnFillPalette               = "BlueDarkRed18"
     opts_mbgt@cnLevels                    = (/-10.0, -5.0, -3.0, -2.0, -1.0, -0.5, -0.2, \
                                                0, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0/)
     opts_mbgt@cnFillColors                = (/0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17/)
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
    opts_mbgt@gsnMaximize                  = False   ; Maximise plot size
    opts_mbgt@gsnAddCyclic                 = False
    opts_mbgt@tmXBMode                     = "Explicit"  ; Set tick marks explicitly
;    opts_mbgt@tmXBLabelFontHeightF	   = 0.02    ; Label font height 

    if (rad0 .eq. 0) then
     opts_mbgt@tmXBValues                   = (/0.0,0.25,0.5,0.75,1.0,1.25,1.5/)
     opts_mbgt@tmXBLabels                   = (/"0","25","50","75","100","125","150"/)
    else if (rad0 .eq. 1) then
     opts_mbgt@tmXBValues                   = (/0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0/)
     opts_mbgt@tmXBLabels                   = (/"0","25","50","75","100","125","150","175","200"/)
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

  ; Resources to overlay zero line
    opts_zero                             = True
    opts_zero@cnFillOn                    = False
    opts_zero@cnLineColor                 = "grey50"
    opts_zero@cnLineDashPattern           = 3
    opts_zero@cnLevelSelectionMode        = "ExplicitLevels"
    opts_zero@cnLevels                    = 0.0
    opts_zero@cnLineThicknessF            = 3.0
    opts_zero@cnInfoLabelOn               = False
    opts_zero@gsnDraw                     = False       ; Do not draw the plot
    opts_zero@gsnFrame                    = False       ; Do no advance the frame
    opts_zero@gsnContourZeroLineThicknessF = 2.0
    opts_zero@cnLineLabelsOn               = False      ; Turn off line labels
    opts_zero@tiMainString                = ""
    opts_zero@gsnLeftString               = ""
    opts_zero@gsnRightString              = ""

;===============================================================
; Resources to overlay panel labels without using 'gsn_panel'
;===============================================================

    txid			= new(1,graphic)
    amid	    	        = new(1,graphic)

  ; Label text resources 
    txres                       = True
    txres@txPerimOn             = True
    txres@txFontHeightF         = 0.010
    txres@txBackgroundFillColor = "White"

  ; Label position resources 
    amres 	   	        = True
    amres@amParallelPosF	= 0.5			; Right edge (-0.5 for left edge)
    amres@amOrthogonalPosF	= -0.5			; Top edge
    amres@amJust		= "TopRight"
    panel_strings               = (/"Local tangential wind tendency"/)

;=================
; Plot the data 
;=================

  ; Plotting the LHS of the tangential momentum budget (Huang et al. 2018)
  ; Ignore the perturbation PGF term because it is several orders of magnitude smaller

  ; LHS (1): d{v}/dt - Local tangential wind tendency
  ; Output file location and type (including w and z threshold information)
    if (typ .eq. "v" .or. typ .eq. "slp") then
     output = "$sam/nepartak/images/cyl_coords/mbgt5_tan_cxz_"+typ+thr+"_group"+grp+\
              "_sc"+sc+"_"+hy3+"km_"+w0+"w_mean"+mean0
    else
     if (typ .eq. "novc") then 
      output = "$sam/nepartak/images/cyl_coords/mbgt5_tan_cxz_"+typ+"_"+nts+"_"+ntf+"_group"+grp+\
               "_sc"+sc+"_"+hy3+"km_"+w0+"w_mean"+mean0
     else
      output = "$sam/nepartak/images/cyl_coords/mbgt5_tan_cxz_"+typ+"_group"+grp+\
               "_sc"+sc+"_"+hy3+"km_"+w0+"w_mean"+mean0
     end if 
    end if 
    wks    = gsn_open_wks(opt,output)

  ; Local tendency of the mean tangential wind [V_t]
    opts_mbgt@vpXF       = 0.10
    opts_mbgt@vpYF       = 0.90
    opts_mbgt@vpWidthF   = 0.70
    opts_mbgt@vpHeightF  = 0.70
    panel7      = gsn_csm_contour(wks,lhs1_c(:,:),opts_mbgt) ; [V_t]
    cont0       = gsn_csm_contour(wks,lhs1_c(:,:),opts_zero) ; Draw contours

    overlay(panel7, cont0)                                   ; Overlay contours
    txid(0)    = gsn_create_text(wks, panel_strings(0), txres)
    amid(0)    = gsn_add_annotation(panel7, txid(0), amres)
    delete([/cont0/])

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    pres			= True 
    maximize_output(wks,pres)				   ; Calls 'draw' and 'frame'

; End IF statement (choose whether to plot in addition to writing out to text files) 
  end if 

end 