; Calculate the radius of gale force wind following a TC, and write to 1D text file 

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/setup_cyl_test.ncl"
load "$sam/ncl_func/st_centre.ncl"

; ncl dist=3.2 sm=0 cn0=\"geo_sm\" ar=3.0 nr=61 nt=73
; w0=0.3 sc=5 trb=1 thr=18.0 grp=\"nvc\" ts=30 tf=36 ti=6 vc_size_paper.ncl

; 'dist' = size of box following storm (degrees)
; 'sc'   = method of choosing storm centre position (5 = min wind speed)
; 'w0'   = vertical velocity threshold for choosing storm centre position (0.3)
; 'ar'   = outer radius for integration (degrees)
; 'nr'   = number of radial bands between centre and outer radius
; 'grp'  = group ("nvc", "vc", "nvcs", "vcs")
; 'trb'  = turn troubleshooting on (1) or off (0)
; 'ts'   = start time index for non-VC composite (30)
; 'tf'   = end time index for non-VC composite (36)
; 'ti'   = number of hours to include in VC composite, before the start of each VC (6)

begin

;============
; Prelash 
;============

; Calculations over a layer, not multiple vertical levels
  num_l      = 1

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

; Increase memory before reading in data
  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 10000000000
  end setvalues

;==========================================================================
; (1) Read in all data for selected 6-h periods for all non-VC forecasts
;==========================================================================

; Read in non-VC data after getting file size information
  size_novc0 = systemfunc("wc -l < $ar/text/novc_dat.txt" )
  size_novc  = toint(size_novc0)

; Initialisation time and simulation information
  dat_novc   = asciiread("$ar/text/novc_dat.txt",size_novc,"string")
  ens0_novc  = asciiread("$ar/text/novc_sim.txt",size_novc,"string")

; Start and end times for the basic VC and non-VC datasets (add option for user interaction)    
  ts0_novc   = new(size_novc, integer)
  tf0_novc   = new(size_novc, integer)
  ts0_novc   = ts
  tf0_novc   = tf

; Number of points in array 
  nv_size    = dimsizes(dat_novc) * ( (tf0_novc(0) - ts0_novc(0))+1 )

;================================================================
; (2) Also read in data for the VC subset of the VC forecasts
;================================================================

; Get size of file containing all VC forecast information
  size_vc0   = systemfunc("wc -l < $ar/text/vc_dat.txt" )
  size_vc    = toint(size_vc0)

; Initialisation time and simulation information
  dat_vc     = asciiread("$ar/text/vc_dat.txt",size_vc,"string")
  ens0_vc    = asciiread("$ar/text/vc_sim.txt",size_vc,"string")

; 04/12/2019 --> array containing the start times of each VC event 
; For this analysis, we want to look at the pre-VC environment, and so set this array as 'tf0'
  tf0_vc     = (/46, 67, 53, 39, 66, 57, 34, 62, 44,\
                 51, 58, 62, 41,\
                 27, 64, 48, 28, 40, 55, 33, 36,\
                 59, 59, 29, 37, 52, 40, 48, 35/)

; Look at the 6-h period before the start of VC development, for each simulation
  ts0_vc     = tf0_vc - ti

; Calculate how many times we need                                                             
  nsim       = (tf0_vc - ts0_vc) + 1
  v_size     = sum(nsim)

; Sum of each of the values in 'nsim' for use when filling array in loop below                 
  ns         = new(dimsizes(nsim),"integer")

; Calculate sum of each of the values in 'nsim' for use when filling array in loop below       
  do nt0     = 0, dimsizes(nsim)-1
   ns(nt0)   = sum(nsim(0:nt0) ) - 1
  end do

; Add one to the values above                                                                  
  ns1        = ns + 1

;==================================================================
; Read in desired dataset (VC or non-VC; subset or full dataset)
;==================================================================

  pii          = 3.14159265
  radii        = fspan(0.0,ar,nr)
  thetas       = fspan(0.0,360.0,nt) 
  thetas_rad   = (pii/180)*thetas   
  rad_size     = dimsizes(radii)-1

; 3D arrays to hold 'xpos' and 'ypos' values for all analysis times 
  if (grp .eq. "nvc") then 
   size = nv_size
  elseif (grp .eq. "vc") then 
   size = v_size
  end if 

; Create array to hold storm track positions using different methods 
  centre_all  = new( (/7,num_l,size,2/), "float")

; Select group depending on user input argument ("grp")
  xpos_all    = new((/dimsizes(radii),dimsizes(thetas),size/),float)
  ypos_all    = new((/dimsizes(radii),dimsizes(thetas),size/),float)
  xcen_all    = new((/dimsizes(radii),dimsizes(thetas),size/),float)
  ycen_all    = new((/dimsizes(radii),dimsizes(thetas),size/),float)

; 1D array to hold the radius of gale force winds for each dataset 
  arr_all     = new( (/size/), "float")

; Add metadata                                                                                  
  arr_all!0   = "time"

; Set other arrays now to save space in loops below 
  if (grp .eq. "nvc") then
    dat  = dat_novc
    ens0 = ens0_novc
    ts0  = ts0_novc
    tf0  = tf0_novc
  elseif (grp .eq. "vc") then
    dat  = dat_vc
    ens0 = ens0_vc
    ts0  = ts0_vc
    tf0  = tf0_vc
  end if

;======================================
; (1) Read in data for chosen group 
;======================================

  if (grp .eq. "nvc") then 
   print("Reading in non-VC forecast data...")
  elseif (grp .eq. "vc") then 
   print("Reading in VC forecast data...")
  end if 

; Initialise counter variable 
  ct    = 0

; Loop over all forecasts                                                                  
  do st = 0, dimsizes(dat)-1

 ; Print out useful information to the screen                                                 
   print_clock("Working on simulation "+ens0(st)+" ("+dat(st)+")")

 ; Define path for model-level data
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

 ; Calculate storm motion using built-in function ('st_rm')  
 ; Zonal; meridional; vector wind; times; lat; lon; centre                         
   storm_rel = st_rm(dat(st), ens0(st), "$ar/text/", dist, 0)
   u_cyc     = storm_rel[0]
   v_cyc     = storm_rel[1]
   vel_cyc   = storm_rel[2]
   numTIMES  = storm_rel[3]
   lat_arr   = storm_rel[4]
   lon_arr   = storm_rel[5]
   centre    = storm_rel[6]

 ; Size of arrays below 
   dsize     = toint(dist * 50) - 1

 ; Read in 'pc' (121) and 'pd' (119) stream data
   c  = addfile(fili_c,"r")
   d  = addfile(fili_d,"r")

 ; Loop over chosen times 
   do it     = ts0(st), tf0(st)

    print("it = "+it+" ; ct = "+ct)

  ; Get domain boundaries at each time from storm track (calculated off-line)
    t0       = lat_arr(cn,it-1,0)
    t1       = lat_arr(cn,it-1,1)
    n0       = lon_arr(cn,it-1,0)
    n1       = lon_arr(cn,it-1,1)

  ; Combine values above into array for input into external function 
    ll_arr = (/t0,t1,n0,n1/)

  ; Read in basic variables from 'pc' stream (lon, lat, model levs)
    lon1C  = c->longitude_1({n0:n1})
    lat1C  = c->latitude_1({t0:t1}) 
    lonC   = c->longitude({n0:n1})  
    latC   = c->latitude({t0:t1})   
    hybC   = c->hybrid_ht(:)        

    dy     = latC(1) - latC(0)
    dx     = lonC(1) - lonC(0)

    u_varname = "x-wind"
    v_varname = "y-wind"

  ; Read in horizontal wind components (zonal, meridional) 
    u0  = c->$u_varname$(it,13:25,{t0:t1},{n0:n1})
    v0  = c->$v_varname$(it,13:25,{t0:t1},{n0:n1})

  ; Storm-relative wind calculations
    u = u0 - u_cyc(cn,it-1)
    v = v0 - v_cyc(cn,it-1)

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
    lon1D  = d->longitude_1({n0:n1})                ; longitude [xxx grid points]
    lat1D  = d->latitude_1({t0:t1})                 ; latitude  [xxx grid points]
    lonD   = d->longitude({n0:n1})                  ; longitude [xxx grid points]
    latD   = d->latitude({t0:t1})                   ; latitude  [xxx grid points]
    hybD   = d->hybrid_ht(:)                        ; 63 model ('rho') levels 

  ; Read in pressure, vertical velocity, geopotential height and potential temperature 
    p      = d->p(it-2,13:25,{t0:t1},{n0:n1})      ; Pressure (Pa)
    w      = d->dz_dt(it-2,13:25,{t0:t1},{n0:n1})  ; Vertical velocity (m s-1)
    z      = d->ht(it-2,13:25,{t0:t1},{n0:n1})     ; Geopotential height (m)

  ; Convert units of pressure and geopotential height before continuing 
    p      = p / 100
    z      = z / 10

  ; Add metadata
    p@units = "hPa"
    z@units = "dam"

;=============================
; Calculate layer averages
;=============================

    hy0 = sprintf("%0.0f",hybD(13) )
    hy1 = sprintf("%0.0f",hybD(25) )

  ; Read in variables on model levels 
    print("Working on layer between " + hy0+ " and "+hy1+" m AGL")
    
  ; 'pc' stream
    if (grp .eq. "nvcs") then 
     u_plane     = dim_avg_n_Wrap( u(:,:dsize,:dsize), 0)
     v_plane     = dim_avg_n_Wrap( v(:,:dsize,:dsize), 0)
     vort_plane  = dim_avg_n_Wrap( vort(:,:dsize,:dsize), 0)
     avo_plane   = dim_avg_n_Wrap( avo(:,:dsize,:dsize), 0)
    else
     u_plane     = dim_avg_n_Wrap( u(:,:dsize,:dsize), 0)
     v_plane     = dim_avg_n_Wrap( v(:,:dsize,:dsize), 0)
     vort_plane  = dim_avg_n_Wrap( vort(:,:dsize,:dsize), 0)
     avo_plane   = dim_avg_n_Wrap( avo(:,:dsize,:dsize), 0)
    end if 

  ; 'pd' stream 
    if (grp .eq. "nvcs") then 
     w_plane     = dim_avg_n_Wrap( w(:,:dsize,:dsize), 0)
     prs_plane   = dim_avg_n_Wrap( p(:,:dsize,:dsize), 0)
    else
     w_plane     = dim_avg_n_Wrap( w(:,:dsize,:dsize), 0)
     prs_plane   = dim_avg_n_Wrap( p(:,:dsize,:dsize), 0)
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
    centre_arr    = st_centre(0, it, cn, r0, centre, vort_plane, prs_plane, spd_plane, \
                              centre_all, dy, dx, "comp", ct, "pd", "off")

  ; Set 'ilev' to zero (vertical levels are ignored here)
    ilev = 0                               

  ; Updated storm track information contained in 'centre_all'
  ; centre_all([0,1,2,3,...]) --> vort_max, prs_min, prs_min_sm, prs_min_sm2, vort_min, spd_min 
    centre_all     = centre_arr[0]

;====================================================================
; Make sure that storm 'centre' is not within the eyewall updraft
;====================================================================

  ; Retrieve coordinates of storm centre (should match what we used for 'pd' stream above)
    cen_loc0 = centre_all(sc,ilev,ct,0)
    cen_loc1 = centre_all(sc,ilev,ct,1)

  ; Coordinates of grid points either side                                               
    cen_0n = cen_loc0 - dy
    cen_0x = cen_loc0 + dy
    cen_1n = cen_loc1 - dy
    cen_1x = cen_loc1 + dy

  ; Calculate vertical velocity and divergence at these grid points                      
    cen_w    = w_plane( {cen_0n:cen_0x}, {cen_1n:cen_1x} )

  ; If average over several grid points is > threshold value, search for centre again    
    ave_w    = avg(cen_w)
    ave_out  = sprintf("%0.1f",ave_w)

  ; Coordinates of grid points either side (bigger grid than above)                      
    delete([/cen_0n, cen_0x, cen_1n, cen_1x/])
    cen_0n = cen_loc0 - (2 * dy)
    cen_0x = cen_loc0 + (2 * dy)
    cen_1n = cen_loc1 - (2 * dy)
    cen_1x = cen_loc1 + (2 * dy)

    if (ilev .ge. 0 .and. ave_w .gt. w0) then
     print("ave. vertical velocity = "+ave_w)

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

     lat_min_w = ltN(0) + (dy * inds_w(0,0) )
     lon_min_w = lnN(0) + (dx * inds_w(0,1) )
     delete([/ltN, lnN, w_1d, w_pl_new/])

   ; Print new storm centre to screen                                                    
     print("Centre (w_min): "+lat_min_w+" degrees N, "+lon_min_w+" degrees E")

   ; Replace values in 'centre_all' before calling 'setup_cyl_test' below                     
     centre_all(sc,ilev,ct,0) = lat_min_w
     centre_all(sc,ilev,ct,1) = lon_min_w

    end if

  ; Tidy up                                                                              
    delete([/cen_w, cen_0n, cen_0x, cen_1n, cen_1x/])

;=====================================================================
; Call external function to switch to cylindrical coordinate system
;=====================================================================

  ; Call 'setup_cyl_test' and output pressure/vorticity on cylindrical grid
  ; '0' before "single" --> basic variables output
  ; '1' on last row --> layer calculations rather than multiple levels
  ; 'numTIMES' doesn't have any effect on the output of this function 
    cyl_arr  = setup_cyl_test(it, dist, ll_arr, lonC, latC, lonD, latD, lon1D, lat1D, \
                              pc_plane, pd_plane, numTIMES, centre_all, \
                              sc, sm, ar, nr, ilev, num_l, 0, "comp", ct, \
                              1, xpos_all, ypos_all, xcen_all, ycen_all)

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
    vt_int   = cyl_arr[14]
   
    xpos     = cyl_arr[17]
    ypos     = cyl_arr[18]
    xcen     = cyl_arr[20]
    ycen     = cyl_arr[21]

    print("Storm centre = "+lat_max+" ; "+lon_max)

  ; Define new arrays to hold azimuthally-averaged variables 
    spd_azi  = new( (/dimsizes(radii)/), float)
    vt_azi   = new( (/dimsizes(radii)/), float)

  ; Now calculate the azimuthal average of tangential and total wind 
    do irad  = 0, dimsizes(radii)-1
     spd_azi(irad)  = avg( spd_int(irad,:) )
     vt_azi(irad)   = avg( vt_int(irad,:) )
    end do 

;=======================================================================
; Find the radius at which wind speed falls below gale force (18 m/s)
;=======================================================================

    do irad  = 1, dimsizes(radii)-1
   ; Calculate the radial gradient at each radius
     dr = spd_azi(irad) - spd_azi(irad-1)
   ; Troubleshooting
     if (trb .eq. 1) then 
      print("dr = "+dr+" ; spd_azi(irad) = "+spd_azi(irad) )
      print("ct = "+ct)
     end if 
   ; If the gradient is negative (decreasing outward) and speed is below gale force
     if (dr .lt. 0 .and. spd_azi(irad) .lt. thr) then
    ; Then print the radius to screen, and add the value to the large 1D array
      print("Radius of gale force winds is "+radii(irad)+" at T+"+it)
      arr_all(ct) = (/radii(irad)/)
    ; Advance counter variable (27/11/2019 --> can't have this line outside the IF statement)
    ; b/c 'ct' would increase in value each time "nvcs" should be ignored (L509)
      ct = ct + 1 
    ; Make sure that the counter variable has increased in value by 1 
      if (trb .eq. 1) then 
       print("New value of ct = "+ct)
      end if 
    ; Now abort the loop and move onto the next time interval
      break
     end if
    end do   

 ; Tidy up 
   delete([/u_plane,v_plane,vort_plane,avo_plane,w_plane,prs_plane,spd_plane,vort_pl_smth/])
   delete([/pc_plane,pd_plane/])
   delete([/lon1C,lat1C,lonC,latC,u0,v0,u,v,vort,avo,f0,f1/])
   delete([/lon1D,lat1D,lonD,latD,p,w,z/])

  end do ; End loop over times (do it = ts0(st), ...)

; Tidy up 
  delete([/time_c,time_d/])

 end do  ; End loop over simulations (do st = ...)

;=============================
; Write data to text files 
;=============================

; Output file path 
  if (grp .eq. "novc") then 
   diri_out = "$ar/text/"+grp+"_"+thr+"_"+ts+"_"+tf+".txt"
  else
   diri_out = "$ar/text/"+grp+"_"+thr+"_tm"+ti+".txt"
  end if 

; Write the 1D array to this text file using asciiwrite 
  asciiwrite(diri_out, arr_all)

end 