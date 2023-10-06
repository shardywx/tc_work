; Calculate inner-core structure (ring-like vs monopole) for 5-min data and write to text files

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/setup_cyl.ncl"
load "$sam/ncl_func/st_centre.ncl"

; ncl opt=\"x11\" rad0=2.0 dist=1.0 sm=0 zoom=1 lgd0=\"tr\" cn0=\"slp\" lay=1 cyl=1 dnt=0 prs=0
; trb=0 ar=0.75 nr=16 nt=73 rad1=500 w_check=1 w0=0.3 var0=\"vort\" sc=5 tc5_ring_mono.ncl

; 'dat'  = initialisation time string ("02T12", "03T00", ...)
; 'ens0' = ensemble simulation (em00-em11)
; 'int'	 = 
; 'opt'  = output file format ("pdf" or "x11")
; 'rad0' = 
; 'dist' = size of box following storm (degrees)
; 'lay'  = calculate layer average (1) or single level (0)
; 'p_cp' = pressure level for calculating storm centre
; 'lgd0' = legend position ("tr","tl","br","bl")
; 'sc'   = method of choosing storm centre position
; 'ke'   = plot angular velocity ("01") or relative vorticity ("...")

; T+12    = 131 (144)         ; T+18  = 203 (216)          ; T+24    = 275 (288)
; T+30    = 347 (360)         ; T+36  = 419 (432)          ; T+42    = 491 (504)
; T+48    = 563 (576)         ; T+54  = 635 (648)          ; T+60    = 707 (720)
; T+66    = 779 (792)         ; T+72  = 851 (864)          ; T+78    = 923 (936) 
; T+84    = 995 (1008)        ; T+90  = 1067 (1080) 

begin

;=====================================
; Read in 'pc' and 'pd' stream data
;=====================================

; Which set of simulations do we want to analyse?
  dat  = "02T12"
  ens0 = "em11"

; T+48 to T+84 (T+30 = 347) 
  ts0  = 563
  tf0  = 995

; Total number of times to analyse
  ntot = (tf0 - ts0) + 1 

; Define file path 
  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0+"/vc"
  fili_p      = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0

; Other important dimensions
  nlev        = 63

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
  llbox     = toint(dist*50)                   ; Calculate domain size (grid points)            
  dsize     = (/llbox,llbox/)                  ; Domain size (grid points)                      
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

; Create arrays to hold time information for ALL times 
  time_str0  = new(numINPUT_c,string)
  time_arr0  = new(numINPUT_c,string)
  title_arr0 = new(numINPUT_c,string)

  time_str0  = hr_min + " UTC " + sprinti("%0.2i ", day) + month_abbr(month)
  time_arr0  = hr_min+ "Z " + sprinti("%0.2i", day) + " " + month_abbr(month)
  title_arr0 = "Valid at "+time_str0

; Now get the subset of dates we're interested in  
; 18 UTC 3 July (it = 347) --> 00 UTC 6 July (it = 995)
  time_str   = time_str0(ts0:tf0)
  time_arr   = time_arr0(ts0:tf0)
  title_arr  = title_arr0(ts0:tf0)

; Print out relevant information to user about storm tracking 
  sc_arr      = (/"relative vorticity",\                          ; sc = 0
                  "pressure",\                                    ; sc = 1
                  "smoothed pressure",\                           ; sc = 2
                  "extra smoothed pressure",\                     ; sc = 3 
                  "vort_min during symmetric phase",\             ; sc = 4 
                  "windspeed minimum",\                           ; sc = 5
                  "geopotential"/)                                ; sc = 6

  print("Using "+sc_arr(sc)+" to calculate storm centre position...")  

;========================================
; Create arrays before reading in data
;======================================== 

; Array to hold values of relative vorticity ratio
  vr_arr    = new(ntot,"float")

; Array to hold values of maximum standard deviation of relative vorticity
  std_arr   = new(ntot,"float")

; Choose which storm track to read in 
  if (cn0 .eq. "slp") then
   cn = 0
  elseif (cn0 .eq. "slpf") then
   cn = 1
  else
   print("Wrong storm tracking method used! Exiting...")
   exit()
  end if

; Increase system memory before reading in data
  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 10000000000
  end setvalues

;==========================================================
; Read in required variables from 'pc' and 'pd' streams 
;==========================================================

; Create array to hold storm track positions using different methods 
; 18/12/2019 --> change second argument from 'nlev' to '1'
  centre_new  = new((/7,1,ntot,2/),"float")

  pii         = 3.14159265
  radii       = fspan(0.0,ar,nr)
  thetas      = fspan(0.0,360.0,nt) 
  thetas_rad  = (pii/180)*thetas   
  rad_size    = dimsizes(radii)-1

; 3D arrays to hold 'xpos' and 'ypos' values for all analysis times 
; 18/12/2019 --> every array defined in this script, with 'ntot' times, must be referenced with 'ct'
  xpos_all    = new((/dimsizes(radii),dimsizes(thetas),ntot/),float)
  ypos_all    = new((/dimsizes(radii),dimsizes(thetas),ntot/),float)
  xcen_all    = new((/dimsizes(radii),dimsizes(thetas),ntot/),float)
  ycen_all    = new((/dimsizes(radii),dimsizes(thetas),ntot/),float)

; Read in all data streams
  c  = addfile(fili_c,"r")              ; Read in 'pc' stream [u,v,vort]
  d  = addfile(fili_d,"r")              ; Read in 'pd' stream [w,p]                              
  j  = addfile(fili_j,"r")              ; Read in 'pj' stream [t,theta,mixing ratio]             
  k  = addfile(fili_k,"r")              ; Read in 'pk' stream [tau,PV,w] 

; Initialise counter variable
  ct = 0

; Loop over chosen times 
  do it	       = ts0, tf0

   print("Working on time: "+time_str(ct)+" (T+"+it+")" )
   title_arr(ct) = "Valid at "+time_str(ct)+" (T+"+it+")"

;=====================================================  
; Read in variables from 'pc' stream (model levels)
;=====================================================  

 ; Get domain boundaries at each time from storm track (calculated off-line)
   t0          = lat_arr(cn,it,0)
   t1          = lat_arr(cn,it,1)
   n0          = lon_arr(cn,it,0)
   n1          = lon_arr(cn,it,1)
   print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

 ; Combine values above into array for input into external function 
   ll_arr = (/t0,t1,n0,n1/)

 ; Read in basic variables from 'pc' stream
   lonC   = c->longitude({n0:n1})                ; longitude ['d1' grid points]          
   latC   = c->latitude({t0:t1})                 ; latitude  ['d0' grid points]  
   hybC   = c->hybrid_ht(:)                      ; 63 model ('theta') levels
   hyb_plot = hybC / 1000
   hyb_plot@units = "km"
   hyb_m  = sprintf("%0.0f",hybC(:))

   dy0    = latC(1) - latC(0)
   dx0    = lonC(1) - lonC(0)

   u_varname = "x-wind"
   v_varname = "y-wind"

 ; Read in horizontal wind components 
   if (lay .eq. 1) then 
    u0  = c->$u_varname$(it,13:25,{t0:t1},{n0:n1})    ; Zonal wind
    v0  = c->$v_varname$(it,13:25,{t0:t1},{n0:n1})    ; Meridional wind
   else
    u0  = c->$u_varname$(it,ilev,{t0:t1},{n0:n1})     ; Zonal wind 
    v0  = c->$v_varname$(it,ilev,{t0:t1},{n0:n1})     ; Meridional wind
   end if 

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

 ; Calculate relative vorticity using centered finite differences 
   vort = uv2vr_cfd(u, v, latC, lonC, 2)
   copy_VarCoords(u, vort)
   vort@units = "s~S~-1~N~"
   vort@name  = "Relative vorticity on model levels"

 ; Convert to standard units for plotting 
   vort       = vort * (10 ^ 4) 
   vort@units = "10~S~-4~N~ s~S~-1~N~"

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

;======================================
; Read in variables from 'pd' stream
;======================================

 ; Read in basic variables 
   hybD   = d->hybrid_ht(:)                        ; 63 model ('rho') levels 

 ; Read in pressure and vertical velocity 
   if (lay .eq. 1) then  
    p     = d->p(it,13:25,{t0:t1},{n0:n1})      ; Pressure (Pa)
    z     = d->ht(it,13:25,{t0:t1},{n0:n1})     ; Geopotential height (m)
   else
    p     = d->p(it,ilev,{t0:t1},{n0:n1})       ; Pressure (Pa)
    z     = d->ht(it,ilev,{t0:t1},{n0:n1})      ; Geopotential height (m) 
   end if 

 ; Convert units of pressure and geopotential height before continuing 
   p      = p / 100
   z      = z / 10

 ; Add metadata
   p@units = "hPa"
   z@units = "dam"

;===============
; 'pj' stream
;===============

 ; Read in temperature and potential temperature 
   if (lay .eq. 1) then
    t     = j->temp(it,13:25,{t0:t1},{n0:n1})
    th    = j->theta(it,13:25,{t0:t1},{n0:n1})
   else
    t     = j->temp(it,ilev,{t0:t1},{n0:n1})
    th    = j->theta(it,ilev,{t0:t1},{n0:n1})
   end if

;===============
; 'pk' stream
;===============

 ; Read in basic variables                                                               
   hybK   = k->hybrid_ht(:)

 ; Alternative lat/lon grid                                                              
   lon1K  = k->longitude_1({n0:n1})
   lat1K  = k->latitude_1({t0:t1})
   lonK   = k->longitude({n0:n1})
   latK   = k->latitude({t0:t1})

 ; Read in vertical velocity                                                             
   if (lay .eq. 1) then
    w     = k->dz_dt(it,13:25,{t0:t1},{n0:n1})
   else
    w     = k->dz_dt(it,ilev,{t0:t1},{n0:n1})
   end if

;==========================
; Loop over model levels 
;==========================

   hy0 = sprintf("%0.0f",hybD(13) )
   hy1 = sprintf("%0.0f",hybD(25) )

 ; Read in variables on model levels 
   if (lay .eq. 1) then 

    print("Working on layer between " + hy0+ " and "+hy1+" m AGL")

  ; 'pc' stream
    u_plane      = dim_avg_n_Wrap( u(:,:,:), 0)
    v_plane      = dim_avg_n_Wrap( v(:,:,:), 0)
    vort_plane   = dim_avg_n_Wrap( vort(:,:,:), 0)
    avo_plane    = dim_avg_n_Wrap( avo(:,:,:), 0)

  ; 'pd' stream 
    w_plane      = dim_avg_n_Wrap( w(:,:,:), 0)
    prs_plane    = dim_avg_n_Wrap( p(:,:,:), 0)
    th_plane     = dim_avg_n_Wrap( th(:,:,:),0)

   else

    hy2 = sprintf("%0.0f",hybD(ilev))
    print("Working on " + hy2+ " m AGL")

  ; 'pc' stream
    u_plane      = u(:,:)
    v_plane      = v(:,:)
    vort_plane   = vort(:,:)
    avo_plane    = avo(:,:)

  ; 'pd' stream
    w_plane      = w(:,:)
    prs_plane    = p(:,:)
    th_plane     = th(:,:)

   end if 

 ; Calculate vector windspeed and smooth 
   spd_plane     = sqrt( (u_plane ^ 2) + (v_plane ^ 2) )
   copy_VarCoords(u_plane, spd_plane)
   vort_pl_smth  = smth9_Wrap(vort_plane, 0.5, 0.5, True)

 ; Final two 'avo_plane' options are placeholders 
   pc_plane      = (/u_plane, v_plane, spd_plane, vort_plane, \
                     avo_plane, avo_plane, avo_plane, avo_plane/)

 ; Final two 'prs_plane' options are placeholders 
   pd_plane      = (/w_plane, prs_plane, prs_plane, prs_plane/)

;=============================================================
; Call external function to calculate storm centre position
;=============================================================

 ; Set radial distance (º) when looking for vort min in ringlike phase 
   r0            = 0.10

 ; 18/12/2019 --> "single" vs "comp" argument does nothing when res="5min"
   centre_arr    = st_centre(0, it, cn, r0, centre, vort_plane, prs_plane, spd_plane, \
                             centre_new, dy0, dx0, "single", ct, "pd", "off", "5min")

 ; Updated storm track information contained in 'centre_new'
 ; centre_new([0,1,2,3,...]) --> vort_max, prs_min, prs_min_sm, prs_min_sm2, vort_min, spd_min 
   centre_new     = centre_arr[0]

   print("Storm centre = "+centre_new(sc,0,ct,0)+" ; "+centre_new(sc,0,ct,1) )

;======================================================================                           
; Make sure that the storm 'centre' is not within the eyewall updraft                             
;======================================================================      

    if (w_check .eq. 1) then

     print("Original method sampled eyewall updraft. Recalculating...")

   ; Retrieve coordinates of storm centre                                                         
   ; 18/12/2019 --> argument 1 is '0', because 'centre_new_u' only has one vertical level         
   ; Same as for call to 'setup_cyl' --> we've already vertically averaged over levels            
   ; So the normal notation is redundant throughout this whole script                             
     cen_loc0 = centre_new(sc,0,ct,0)
     cen_loc1 = centre_new(sc,0,ct,1)

   ; Coordinates of grid points either side                                                       
   ; Could edit lines below to increase the area (± 2/3dy0, etc)                     
     cen_0n = cen_loc0 - dy0
     cen_0x = cen_loc0 + dy0
     cen_1n = cen_loc1 - dy0
     cen_1x = cen_loc1 + dy0

   ; Calculate vertical velocity over small subset of grid points either side of original centre    
   ; 13/10/2019 --> add divergence                                                                  
     cen_w  = w_plane( {cen_0n:cen_0x}, {cen_1n:cen_1x} )

   ; If average over several grid points is above a threshold value, continue searching for centre 
   ; 14/10/2019 --> calculate averaged vertical velocity over small box...                          
     ave_w  = avg(cen_w)
     ave_out= sprintf("%0.1f",ave_w)

   ; Coordinates of grid points either side (bigger grid than above)                                
     delete([/cen_0n, cen_0x, cen_1n, cen_1x/])

   ; 14/10/2019 --> ... but use larger box to constrain your search for new storm centre (L944)     
   ; Higher level                                                                                   
     cen_0n = cen_loc0 - (2 * dy0)
     cen_0x = cen_loc0 + (2 * dy0)
     cen_1n = cen_loc1 - (2 * dy0)
     cen_1x = cen_loc1 + (2 * dy0)

   ; Implement modified policy if there is even a hint of enhanced vertical velocity                
     if (ave_w .gt. w0) then

      print("ave. vertical velocity = "+ave_w)

    ;===================================================                                       
    ; Now calculate new centre on higher level/layer                                           
    ;===================================================                                       

    ; Create smaller grid around previous centre (L929-932)                                    
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
    ; 18/12/2019 --> 'mlev1' corresponds to the level at the centre of the vertical layer      
    ; ... previously this argument had been set to '0'                                         
      centre_new(sc,0,ct,0) = lat_min_w
      centre_new(sc,0,ct,1) = lon_min_w
      delete([/lat_min_w, lon_min_w/])

     end if

   ; Tidy up                                                                                   
     delete([/cen_w, cen_0n, cen_0x, cen_1n, cen_1x/])

   end if

;=====================================================================
; Call external function to switch to cylindrical coordinate system
;=====================================================================

 ; Set 'ilev' option before function call below 
   if (lay .eq. 1) then
    ilev = 0
   end if

 ; Call 'setup_cyl' and output pressure/vorticity on cylindrical grid 
 ; Last option = 0 --> basic variables output  
 ; Last option = 1 --> additional MBGT variables output 
 ; 18/12/2019 --> changed "single" to "comp" and 'it' to 'ct' immediately after
   cyl_arr  = setup_cyl(it, dist, ll_arr, lonC, latC, lonK, latK, lon1K, lat1K, \
                        pc_plane, pd_plane, ntot, centre_new, \
                        sc, sm, ar, nr, nt, ilev, nlev, 0, "comp", ct, 1, \
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
   prs_int  = cyl_arr[16]
   
   xpos     = cyl_arr[17]
   ypos     = cyl_arr[18]
   xcen     = cyl_arr[20]
   ycen     = cyl_arr[21]

   print("Storm centre = "+lat_max+" ; "+lon_max)

 ; Define new arrays to hold azimuthally-averaged variables 
   vort_azi = new( (/dimsizes(radii)/), float)
   spd_azi  = new( (/dimsizes(radii)/), float)
   vt_azi   = new( (/dimsizes(radii)/), float)
   va_azi   = new( (/dimsizes(radii)/), float)

 ; Calculate radial distance from storm centre at all grid points (see tc_ang_vel.ncl - L1115)
   drad     = dimsizes(radii)
   rad_arr  = new( (/dimsizes(radii), dimsizes(thetas)/), float)
   rad_arr(1:drad-1,:) = gc_latlon( ycen(1:drad-1,:), xcen(1:drad-1,:), \
                                    ypos(1:drad-1,:), xpos(1:drad-1,:), 2, 3)

 ; Replace zeros at r=0 with small, positive values (avoid divide by zero errors)
   rad_arr(0,:) = rad1

 ; Calculate angular velocity from tangential wind 
   va_int   = vt_int / rad_arr
   va_int   = va_int * (10 ^ 4)
   va_int@units = "10~S~-4~N~ s~S~-1~N~"
   va_int@_FillValue = default_fillvalue("float")

 ; Convert the radius values into 'km'
   rad_arr  = rad_arr / 1000
   rad_arr@units = "km"

 ; Now calculate the azimuthal average of relative vorticity & angular velocity
   do irad  = 0, dimsizes(radii)-1
    vort_azi(irad) = avg( vort_int(irad,:) )
    spd_azi(irad)  = avg( spd_int(irad,:) )
    vt_azi(irad)   = avg( vt_int(irad,:) )
    va_azi(irad)   = avg( va_int(irad,:) )
   end do 

 ; Use more sophisticated method to calculate whether vortex is ringlike or monopolar
 ; This code will go into an NCL function 

 ; (1) Find radius (index) corresponding to relative vorticity maximum
   if (var0 .eq. "ang") then 
    ind_vmax   = ind_resolve(maxind(va_azi), dimsizes(va_azi) )
   else
    ind_vmax   = ind_resolve(maxind(vort_azi), dimsizes(vort_azi) )
   end if 

 ; Calculate vort/ang-vel at ind_vmax ('vort_max'), and at the vortex centre ('vort_centre')
   if (var0 .eq. "ang") then 
    v_max      = va_azi(ind_vmax(0,0) )
    v_centre   = va_azi(0)
   else 
    v_max      = vort_azi(ind_vmax(0,0) )
    v_centre   = vort_azi(0)
   end if 

 ; (2) Calculate the ratio of 'ind_max' to 'ind_centre' 
   vr_arr(ct) = v_centre / v_max
   print("v_centre = "+v_centre+"; v_max = "+v_max)

 ; If 'ind_vmax' = 'ind_centre' (i.e. monopole), then 'vr_arr' will be equal to 1  
   if (vr_arr(ct) .gt. 0.95) then 

    print("Inner-core structure is monopolar (ratio = "+vr_arr(ct)+")")
    sym_str = "Monopole: "+sprintf("%0.2f", vr_arr(ct) )

   elseif (vr_arr(ct) .lt. 0.80) then 

    print("Inner-core structure is ringlike (ratio = "+vr_arr(ct)+")")
    sym_str = "Ringlike: "+sprintf("%0.2f", vr_arr(ct) )

   else 
    print("Inner-core structure is not clearly monopolar or ringlike (ratio = "+vr_arr(ct)+")")
    sym_str = "Not clear: "+sprintf("%0.2f", vr_arr(ct) )

   end if 

 ; Before recalculating and potentially changing value of 'sc', define output file path and type 
   if (lay .eq. 1) then 

    if (dnt .eq. 0) then
     output="$nep/nepartak/images/basic/vr5_sc"+sc+"_"+var0+"_"+time_arr(ct)
    else 
     output="$nep/nepartak/images/basic/vr5_sc"+sc+"a_"+var0+"_"+time_arr(ct)
    end if

   else

    if (dnt .eq. 0) then
     output="$nep/nepartak/images/basic/vr5_sc"+sc+"_"+var0+\
            "_"+time_arr(ct)+"_"+ilev
    else
     output="$nep/nepartak/images/basic/fc_sc"+sc+"a_"+var0+\
            "_"+time_arr(ct)+"_"+ilev
    end if

   end if 

   if (trb .eq. 1) then 
    print("Value of 'sc' = "+sc)
    print("File path = "+output)
   end if 

   wks = gsn_open_wks(opt,output)

 ; Load colour table                                                                               
   gsn_define_colormap(wks,"prcp_new")

 ; If we're in the symmetric phase, put the storm centre in the middle of the vorticity ring 
   if (dnt .eq. 1 .and. vr_arr(ct) .lt. 0.80) then 

    print("Re-calculating storm centre to fit ringlike vorticity profile...")

  ; Select the 'vort_min' method of storm tracking 
    sc0 = 4

  ; Call the external function again
    delete(cyl_arr)
    cyl_arr  = setup_cyl(it, dist, ll_arr, lonC, latC, lonK, latK, lon1K, lat1K, \
                         pc_plane, pd_plane, ntot, centre_new, \
                         sc0, sm, ar, nr, nt, ilev, nlev, 0, "single", ct, 1, \
                         xpos_all, ypos_all, xcen_all, ycen_all)

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
    prs_int  = cyl_arr[16]

    xpos     = cyl_arr[17]
    ypos     = cyl_arr[18]
    xcen     = cyl_arr[20]
    ycen     = cyl_arr[21]

  ; Recalculate azimuthal average of relative vorticity and windspeed
    do irad = 0, dimsizes(radii)-1
     vort_azi(irad) = avg( vort_int(irad,:) )
     spd_azi(irad)  = avg( spd_int(irad,:) )
     vt_azi(irad)   = avg( vt_int(irad,:) )
     va_azi(irad)   = avg( va_int(irad,:) )
    end do

   end if 

;================================================================================
; Throwback --> calculate the maximum standard deviation of relative vorticity
;================================================================================

; Calculate perturbation relative vorticity
  vort_ptb = new( (/dimsizes(radii),dimsizes(thetas)/), "float")
  do irad  = 0, dimsizes(radii)-1
   vort_ptb(irad,:) = vort_azi(irad) - vort_int(irad,:)
  end do

; Calculate standard deviation of relative vorticity 
  std_lev  = new(dimsizes(radii),"float")
  do irad  = 0, dimsizes(radii)-1
   std_lev(irad) = dim_stddev_n_Wrap(vort_ptb(irad,:),0)
  end do

; Calculate maximum value of array above, and read into 'std_arr'
  std_arr(ct) = max(std_lev)

; Print value out as troubleshooting if required 
  if (std .eq. 1) then 
   print("Max. standard deviation = "+std_arr(ct) )
  end if 

; Tidy up
  delete([/w_plane,u_plane,v_plane,spd_plane/])
  delete([/avo_plane,vort_plane,prs_plane,vort_pl_smth/])
  delete([/vort,u,v,u0,v0,latC,lonC/])
  delete([/p,w,z,th,latK,lonK,lat1K,lon1K/])

  ct = ct + 1  ; Advance counter variable (time)

 end do       ; End loop over input times (do it = 0, ntot-1)

;===================================================
; Write out relative vorticity ratio to text file
;===================================================

; Skip this part if we're only looking at a few times (e.g. troubleshooting)
  if (ntot .gt. 48) then 
 
   print("Working through all times; write to text file!")

   if (std .eq. 0) then 

    if (dnt .eq. 0) then 
     vr_out   = "$ar/text/cp5_vr_sc"+sc+".txt"
    else 
     vr_out   = "$ar/text/cp5_vr_sc"+sc+"a.txt"
    end if 
    asciiwrite(vr_out, vr_arr)

   else

    std_out  = "$ar/text/cp5_std_sc"+sc+"_.txt"
    asciiwrite(std_out, std_arr)

   end if 

  end if 

end 