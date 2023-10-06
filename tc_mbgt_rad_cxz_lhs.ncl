; Script to calculate the local radial wind tendency

; Run using:

; ncl opt=\"png\" w0=0.3 z0=0 dist=1.04 clr=5 calc=1 rad0=0 nt=73
; nr=22 ar=1.05 lgd0=\"tr\" plt=0 new0=0 mean0=0 cn0=\"geo_sm\" sm=0 sc=5 w_check=1
; nts=36 ntf=54 max0=0 thr=2 grp=\"8\" typ=\"erc1\" tc_mbgt_rad_cxz_lhs.ncl

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
; 'max'   = 1 (calculate maximum mean tangential wind tendency, and quit)

; 500 m [10], 1 km [14], 1.5 km [17], 2 km [20], 3 km [24], 4 km [28], 5 km [31], 6 km [34]  
; 7 km [37], 8 km [39], 9 km [42], 10 km [45], 11 km [46], 12 km [48], 13 km [50], 14 km [52]
; 15 km [54], 16 km [56]

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

; Read in data from text files according to VC phase
  if (typ .eq. "sym") then
   typ0   = "ring"
  elseif (typ .eq. "s2a") then
   typ0   = "r2m"
  elseif (typ .eq. "asym") then
   typ0   = "mono"
  elseif (typ .eq. "a2s") then
   typ0   = "m2r"
  elseif (typ .eq. "novc") then
   typ0 = "novc"
  elseif (typ .eq. "test") then 
   typ0 = "test"
  elseif (typ .eq. "erc1") then
   typ0 = "erc1"
  elseif (typ .eq. "erc2") then
   typ0 = "erc2"
  elseif (typ .eq. "erc3") then
   typ0 = "erc3"
  end if

; Retrieve file size information using 'systemfunc'
  if (typ .ne. "erc1" .and. typ .ne. "erc2" .and. typ .ne. "erc3") then
   fili  = "$ar/text/vc_"+typ+"_dat_group"+grp+".txt"
   size0 = systemfunc("wc -l < "+fili)
   size  = toint(size0)
  end if 

; Which set of simulations do we want to analyse? 

; No VCs
  if (typ .eq. "novc") then
   novc_size = 16
   dat  = asciiread("$ar/text/novc_dat.txt",(/novc_size/),"string")
   ens0 = asciiread("$ar/text/novc_sim.txt",(/novc_size/),"string")
   ts0  = new(novc_size,integer)
   ts0  = nts0
   tf0  = new(novc_size,integer)
   tf0  = ntf0
; First ERC case (em04, 03T00)
  elseif (typ .eq. "erc1") then
   dat  = (/"03T00"/)
   ens0 = (/"em04"/)
   ts0  = (/78/)
   tf0  = (/108/)
; Second ERC case (em11, 03T00)
  elseif (typ .eq. "erc2") then
   dat  = (/"03T00"/)
   ens0 = (/"em11"/)
   ts0  = (/72/)
   tf0  = (/102/)
; Test mode                                                                               
  elseif (typ .eq. "test") then
   dat  = (/"03T00"/)
   ens0 = (/"em04"/)
   ts0  = (/88/)
   tf0  = (/89/)
; Any other VC phase 
  else
   dat  = asciiread("$ar/text/vc_"+typ+"_dat_group"+grp+".txt",size,"string")
   ens0 = asciiread("$ar/text/vc_"+typ+"_sim_group"+grp+".txt",size,"string")
   ts0  = asciiread("$ar/text/vc_"+typ+"_ts_group"+grp+".txt",size,"integer")
   tf0  = asciiread("$ar/text/vc_"+typ+"_tf_group"+grp+".txt",size,"integer")
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

; Total number of elements in time-dimension (sims * times per sim)
  ntot = nts * nsim 

;======================================================================================
; Create array to hold all values for composite (previous used with 'wrt=1' option)
;======================================================================================

; Dimensions of 5D array (levs; radial circles; diagnostics; sims; times)
  nlev       = 59
  ndiag      = 10
  nrad       = nr

  pii        = 3.14159265
  radii      = fspan(0.0,ar,nr)      ; Array of radial distances ('Start','End','Number') 
  thetas     = fspan(0.0,360.0,nt)   ; Array of azimuth angles (0, 5, 10, etc)
  thetas_rad = (pii/180)*thetas      ; Degrees to radians 
  rad_size   = dimsizes(radii)-1     ; Index for later calculation 

  plot_arr0  = new( (/nlev, nts, nrad, ndiag, nsim/), "float")

;============================================================================
; Now calculate tangential wind tendency using centred finite differencing
;============================================================================

; Maximum number of times we could analyse (will never reach this number)
  ntimes_all = nsim * nts

; Now reshape 'plot_arr0' into a 4D array
; 09/10/2019 --> test out 'reshape' function until it works; currently introduces errors
  plot_arr   = reshape(plot_arr0, (/nlev, dimsizes(radii), ndiag, ntimes_all/) )

; Tangential wind (4D)
  mean_vtan_all       = new( (/nlev, dimsizes(radii), nsim, nts/), "float")
  mean_vtan_all!0     = "lev"
  mean_vtan_all!1     = "rad"
  mean_vtan_all!2     = "sim"
  mean_vtan_all!3     = "time"
  mean_vtan_all@description = "Azimuthally averaged tangential wind"
  mean_vtan_all@units = "m s~S~-1~N~"

; Radial wind (4D)
  mean_vrad_all       = new( (/nlev, dimsizes(radii), nsim, nts/), "float")
  mean_vrad_all!0     = "lev"
  mean_vrad_all!1     = "rad"
  mean_vrad_all!2     = "sim"
  mean_vrad_all!3     = "time"
  mean_vrad_all@description = "Azimuthally averaged radial wind"
  mean_vrad_all@units = "m s~S~-1~N~"

; Coriolis parameter (4D)
  mean_cor_all        = new( (/nlev, dimsizes(radii), nsim, nts/), "float")
  mean_cor_all!0      = "lev"
  mean_cor_all!1      = "rad"
  mean_cor_all!2      = "sim"
  mean_cor_all!3      = "time"
  mean_cor_all@description = "Azimuthally averaged Coriolis parameter"
  mean_cor_all@units  = "s~S~-1~N~"

; Maximum tangential wind (2D)
  mean_vtan_max       = new( (/nsim, nts/), "float")
  mean_vtan_max!0     = "sim"
  mean_vtan_max!1     = "time"
  mean_vtan_max@description = "Azimuthally averaged tangential wind"
  mean_vtan_max@units = "m s~S~-1~N~"

; Tendency of maximum tangential wind (2D)
  max_vtan_ten        = new( (/nsim, nts/), "float")
  max_vtan_ten!0      = "sim"
  max_vtan_ten!1      = "time"
  max_vtan_ten@description = "Azimuthally averaged tangential wind"
  max_vtan_ten@units  = "m s~S~-1~N~"

; Radial wind at vmax
  vrad_vmax           = new( (/nsim, nts/), "float")
  vrad_vmax!0         = "sim"
  vrad_vmax!1         = "time"
  vrad_vmax@description = "Azimuthally averaged tangential wind"
  vrad_vmax@units     = "m s~S~-1~N~"

; Coriolis parameter at vmax
  f0_vmax             = new( (/nsim, nts/), "float")
  f0_vmax!0           = "sim"
  f0_vmax!1           = "time"
  f0_vmax@description = "Azimuthally averaged tangential wind"
  f0_vmax@units       = "m s~S~-1~N~"

; Radial wind tendency
  mean_vrad_ten       = new( (/nlev, dimsizes(radii), nsim, nts/), "float")
  mean_vrad_ten!0     = "lev"
  mean_vrad_ten!1     = "rad"
  mean_vrad_ten!2     = "sim"
  mean_vrad_ten!3     = "time"
  mean_vrad_ten@description = "Azimuthally averaged radial wind tendency"
  mean_vrad_ten@units = "m s~S~-1~N~"

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

; Create array to hold storm track positions (L727)
  centre_newer = new( (/7, nlev, ntimes_all, 2/), "float")

; Redefine arrays to hold all 'xpos' and 'ypos' position arrays
  xpos_newer   = new((/nlev, dimsizes(radii), dimsizes(thetas), ntimes_all/),float)
  ypos_newer   = new((/nlev, dimsizes(radii), dimsizes(thetas), ntimes_all/),float)

; Redefine arrays to hold lat/lon at centre of cylindrical grid
  xcen_newer   = new((/nlev, dimsizes(radii), dimsizes(thetas), ntimes_all/),float)
  ycen_newer   = new((/nlev, dimsizes(radii), dimsizes(thetas), ntimes_all/),float)

; Initialise time-varying counter variable
  ct = 0

; Read in data corresponding to times in text files (see L266)                            
  do st    = 0, nsim-1

 ; Print information to screen 
   print_clock("Working on simulation: "+ens0(st)+" ("+dat(st)+")" )

 ; Define file path and data stream (L290)                                                
   diri    = "/nfs/a319/earshar/"+dat(st)+"/"+ens0(st)
   fili_p  = "201607"+dat(st)+"00Z_NPTK_4p4_L80_ra1t_"+ens0(st)
   fili_c1 = systemfunc("cd "+diri+" ; ls "+fili_p+"_pc*.nc")
   fili_c  = diri+"/"+fili_c1
   fili_d1 = systemfunc("cd "+diri+" ; ls "+fili_p+"_pd*.nc")
   fili_d  = diri+"/"+fili_d1

 ; Retrieve important file information (L298)                                             
   numINFO_c  = nc_times(fili_c)
   numINPUT_c = numINFO_c[0]
   time_c     = numINFO_c[1]

   numINFO_d  = nc_times(fili_d)
   numINPUT_d = numINFO_d[0]
   time_d     = numINFO_d[1]

 ; Calculate storm motion using built-in function (L318)                                  
   storm_rel  = st_rm(dat(st), ens0(st), "$ar/text/", dist, 0)
   u_cyc      = storm_rel[0]
   v_cyc      = storm_rel[1]
   vel_cyc    = storm_rel[2]
   numTIMES   = storm_rel[3]
   lat_arr    = storm_rel[4]
   lon_arr    = storm_rel[5]
   centre     = storm_rel[6]

 ; Time and domain-size information                                                              
   llbox = toint(dist*50)
   dsize = (/llbox, llbox/)
   d0    = dsize(0)
   d1    = dsize(1)

 ; Loop over additional times (same as on L736, + immediately before and after)                  
   c     = addfile(fili_c, "r")
   d     = addfile(fili_d, "r")

 ; Initialise another time-varying counter variable (resets after each simulation)
   dt0   = 0

   do it = ts0(st)-1, tf0(st)+1

  ; Print information to screen
    print_clock("Looping between T+"+ts0(st)+" and T+"+tf0(st) )
    print("ct = "+ct+" ; dt0 = "+dt0)

;==============================================                                               
; Read in variables from 'pc' stream (L747)                                                   
;==============================================                                               

    t0    = lat_arr(cn,it-1,0)
    t1    = lat_arr(cn,it-1,1)
    n0    = lon_arr(cn,it-1,0)
    n1    = lon_arr(cn,it-1,1)

  ; Combine values above into array for input into external function                          
    ll_arr = (/t0,t1,n0,n1/)

  ; Read in basic variables from 'pc' stream                                                  
    lon1C  = c->longitude_1({n0:n1})
    lat1C  = c->latitude_1({t0:t1})
    lonC   = c->longitude({n0:n1})
    latC   = c->latitude({t0:t1})

    hybC   = c->hybrid_ht(:)
    hybC@units = "m"
    hyb_plot   = hybC / 1000 
    hyb_plot@units = "km"

    u_varname  = "x-wind"
    v_varname  = "y-wind"

  ; Read in horizontal wind components                                                        
  ; Both variables below -- [hybC | 63] * [latC | 150] * [lonC | 150]                         
    u0  = c->$u_varname$(it,:,{t0:t1},{n0:n1}) ; Zonal wind                                   
    v0  = c->$v_varname$(it,:,{t0:t1},{n0:n1}) ; Meridional wind                              

  ; Calculate storm-relative winds                                                            
    u = u0 - u_cyc(cn,it-1)
    v = v0 - v_cyc(cn,it-1)

  ; Add metadata from original horizontal wind arrays                                         
    copy_VarCoords(u0,u)
    copy_VarAtts(u0,u)

    copy_VarCoords(v0,v)
    copy_VarAtts(v0,v)

  ; Calculate relative vorticity using centred finite differences                             
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

    dx0     = lonC(1) - lonC(0) ; Grid spacing (longitude)
    dy0     = latC(1) - latC(0) ; Grid spacing (latitude)

;=============================================
; Read in variables from 'pd' stream (L894)
;=============================================

  ; Model height levels
    hybD   = d->hybrid_ht(:)
    hybD@units = "m"

  ; Pressure (Pa), vertical velocity (m s-1) and geopotential height (m) 
    p       = d->p(it-2,:,{t0:t1},{n0:n1}) 
    w       = d->dz_dt(it-2,:,{t0:t1},{n0:n1})

  ; Convert units of pressure and add metadata
    p       = p / 100
    p@units = "hPa"

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
     avo_plane              = avo(ilev,:,:)        ; Absolute vorticity
     f_plane                = f0(ilev,:,:)         ; Coriolis parameter
     w_plane                = w(ilev,:,:)          ; Vertical velocity
     prs_plane              = p(ilev,:,:)          ; Pressure 

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

;==========================================                                                   
; Calculate storm-centre position (L996)                                                      
;==========================================                                                   

   ; Set radial distance (º) when looking for speed/vort min (in ringlike phase)              
     r0            = 0.12

     print_clock("Calculating storm centre position...")
     centre_arr    = st_centre(ilev, it, cn, r0, centre, vort_plane, prs_plane, spd_plane, \
                               centre_newer, dy0, dx0, "comp", ct, "pd", "off", "1h")

   ; Updated storm track information contained in 'centre_new'                                
     centre_newer  = centre_arr[0]

;===============================================================================              
; Make sure that the storm 'centre' is not within the eyewall updraft (L1006)                 
;===============================================================================              

   ; Retrieve coordinates of storm centre                                                     
     cen_loc0 = centre_newer(sc,ilev,ct,0)
     cen_loc1 = centre_newer(sc,ilev,ct,1)

   ; Coordinates of grid points either side                                                   
     cen_0n = cen_loc0 - dy0
     cen_0x = cen_loc0 + dy0
     cen_1n = cen_loc1 - dy0
     cen_1x = cen_loc1 + dy0

   ; Calculate vertical velocity at these grid points                                         
     cen_w    = w_plane( {cen_0n:cen_0x}, {cen_1n:cen_1x} )

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

    ; Replace values in 'centre_newer' before calling 'setup_cyl' below
      centre_newer(sc,ilev,ct,0) = lat_min_w
      centre_newer(sc,ilev,ct,1) = lon_min_w

     end if

   ; Tidy up
     delete([/cen_w, cen_0n, cen_0x, cen_1n, cen_1x/])

;=============================================                                           
; Switch to cylindrical coordinates (L1077)                                              
;=============================================                                           

   ; Call user-defined function 'setup_cyl'
     print_clock("Translating to cylindrical grid...")
     cyl_arr  = setup_cyl(it, dist, ll_arr, lonC, latC, lonC, latC, lon1C, lat1C, \
                          pc_plane, pc_plane, ntimes_all, centre_newer, \
                          sc, sm, ar, nr, nt, ilev, nlev, 1, "comp", ct, \
                          0, xpos_newer, ypos_newer, xcen_newer, ycen_newer)

   ; Output the tangential wind on a cylindrical grid 
     vt_int   = cyl_arr[15]
     vr_int   = cyl_arr[14]
     f_int    = cyl_arr[13]

   ; Calculate azimuthally-averaged tangential wind (+ ...)
     do irad  = 0, dimsizes(radii)-1
      mean_vtan_all(ilev,irad,st,dt0) = avg(vt_int(irad,:) )
      mean_vrad_all(ilev,irad,st,dt0) = avg(vr_int(irad,:) )
      mean_cor_all(ilev,irad,st,dt0)  = avg(f_int(irad,:) )
     end do

  ; End loop over model levels (do ilev = 0, nlev-1)                                    
    end do

  ; Increase values of time-varying counter variables 
    dt0 = dt0 + 1
    ct  = ct  + 1 

 ; End loop over times (do it = ts0(st)-1, tf0(st)+1)
   end do

 ; Tidy up before next iteration (21/10/2019)
   delete([/time_c, time_d, storm_rel, u_cyc,v_cyc,vel_cyc,lat_arr,lon_arr,centre/])

; End loop over simulations (do st = 0, nsim-1)                                 
  end do

;=============================================================================
; Calculate tangential wind tendency (LHS) using centred finite differences
;=============================================================================

; Loop over all simulations 
  do st = 0, nsim-1

 ; Loop over time indices (not including first and last)
   do it = 1, td(st)-2

    cp1 = it+1
    cm1 = it-1

    mean_vrad_ten(:,:,st,it) = (mean_vrad_all(:,:,st,cp1) - mean_vrad_all(:,:,st,cm1) ) / 2

 ; End loop over time indices (do it = 1, td(st)-2)
   end do

; End loop over simulations (do st = 0, nsim-1)
  end do

;====================================================================================
; Read all values of 'mean_vrad_ten' (not including missing values) into new array  
;====================================================================================

; Reshape 'mean_vrad_ten' into 3D array 
; 08/10/2019 --> [levs * radii * times]
; 08/10/2019 --> Edit to discard missing values 
  vrad_tend_out0 = reshape(mean_vrad_ten, (/nlev, dimsizes(radii), ntimes_all/) )

;====================================================================
; Write out values at individual times, if required (ERC analysis)
;====================================================================

; If we're analysing ERCs rather than VCs, write out data at each time individually
  if (typ .eq. "erc1" .or. typ .eq. "erc2" .or. typ .eq. "erc3" .or. typ .eq. "test") then

 ; Get dimension sizes
   vsize      = dimsizes(vrad_tend_out0)

 ; Create new array to hold only the non-missing values 
   lhs1_c     = new ( (/vsize(0),vsize(1),vsize(2)-2/), "float")

 ; Rename array before writing out 
   dt = 0
   do ct = 1, vsize(2)-2
    lhs1_c(:,:,dt) = vrad_tend_out0(:,:,ct)
    dt = dt + 1
   end do 

 ; Add metadata
   lhs1_c!0   = "lev"
   lhs1_c&lev = hyb_plot
   lhs1_c!1   = "rad"
   lhs1_c&rad = radii

 ; Set counter variable (time)
   ct = 0

 ; Loop over chosen times 
   do it = ts0(0), tf0(0)

    print("Working on time: T+"+it)

  ; Additional metadata for 'write_matrix' procedure [21 * 59]                              
    ncol           = nr
    nrow           = nlev
    fmt            = nrow + "f10.3"

  ; Resources for writing out                                                               
    opts_mat       = True
    opts_mat@row   = True
    opts_mat@title = ""

  ; Name all variables we want to write                                                     
    opt_arr = (/"lhs1"/)

  ; Array of all composite terms                                                            
    mbgt_arr = (/lhs1_c/)

  ; Call procedure and write time-averaged                                                  
    print_clock("Writing 2D arrays to matrices for further analysis...")
    do mb = 0, dimsizes(opt_arr)-1

     print("Working on "+opt_arr(mb) )

   ; Set file path, etc
     t_dir  = "./rad_mbgt_"
     t_out  = opt_arr(mb)+"_sc"+sc+"_"+w0+"w_"+hy3+"km_"+typ0+"_T"+it+".txt"

     t_path = t_dir + t_out
     print_clock("Writing to file: "+t_path)
     opts_mat@fout = t_path

   ; Check for existing files and then write                                                
     system("/bin/rm -f " + t_path)
     write_matrix( lhs1_c(:,0:nr-2,ct), fmt, opts_mat)

    end do
  ; End loop over diagnostics                                                               

    ct = ct + 1
  ; Advance counter variable (time) 

   end do 
 ; End loop over chosen times 

 ; Either tidy up the text files and exit the script...   
   if (plt .ne. 1) then

   ; Tidy up text files by deleting indices in first column ('0,1,2,3,...')   
     system("sed -i -r 's/.{10}//' ./*"+typ+"_T*.txt")

   ; Move text files to correct directory   
     system("mv ./*"+typ+"_T*.txt $ar/text/rad/")

   ; Now, exit script 
     exit()

 ; ...or, continue to plotting (delete existing arrays before continuing) 
   else
     delete([/lhs1_c,mbgt_arr/])
   end if

  end if 

; Before plotting, average over all chosen time periods
  lhs1_c         = dim_avg_n_Wrap(vrad_tend_out0, 2)
  lhs1_c@name    = "Local radial wind tendency"

; Add metadata before plotting 
  lhs1_c!0       = "lev"
  lhs1_c&lev     = hyb_plot
  lhs1_c!1	 = "rad"
  lhs1_c&rad	 = radii

;========================================================================
; Write values to output matrices (text files) using in-built function
;========================================================================

; Additional metadata for 'write_matrix' procedure [21 * 59]
  ncol           = nrad
  nrow           = nlev
  fmt            = nrow + "f10.3"

  print("Number of columns = "+ncol+" ; number of rows = "+nrow)

; Resources for writing out 
  opts_mat       = True
  opts_mat@row   = True
  opts_mat@title = ""

; Name all variables we want to write 
  opt_arr = (/"lhs1"/)

; Array of all composite terms
  mbgt_arr = (/lhs1_c/)

; Call procedure and write time-averaged
  print_clock("Writing 2D array to matrices for further analysis...")
  do mb = 0, dimsizes(opt_arr)-1

    print("Working on "+opt_arr(mb) )

  ; Set file path, etc
    t_dir  = "./rad_mbgt_"
    t_out  = opt_arr(mb)+"_sc"+sc+"_"+w0+"w_"+hy3+"km_"+typ0+".txt" 
    t_path = t_dir + t_out
    print_clock("Writing to file: "+t_path)
    opts_mat@fout = t_path

  ; Check for existing files and then write 
    system("/bin/rm -f " + t_path)
    write_matrix( mbgt_arr(:,:), fmt, opts_mat)

  end do

;==========================
; Options for plotting 
;==========================

  if (plt .eq. 1) then 

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

  ; Linearise the y-axis
    opts_mbgt@gsnYAxisIrregular2Linear	   = True

  ; y-axis limits 
    opts_mbgt@trYMinF			   = 0.0
    opts_mbgt@trYMaxF			   = 17.6

  ; x-axis limits
    opts_mbgt@trXMinF                      = 0.0
    opts_mbgt@trXMaxF                      = 1.0

  ; Resources to overlay zero line
    opts_zero                              = True
    opts_zero@cnFillOn                     = False
    opts_zero@cnLineColor                  = "grey50"
    opts_zero@cnLineDashPattern            = 3
    opts_zero@cnLevelSelectionMode         = "ExplicitLevels"
    opts_zero@cnLevels                     = 0.0
    opts_zero@cnLineThicknessF             = 3.0
    opts_zero@cnInfoLabelOn                = False
    opts_zero@gsnDraw                      = False       ; Do not draw the plot
    opts_zero@gsnFrame                     = False       ; Do no advance the frame
    opts_zero@gsnContourZeroLineThicknessF = 2.0
    opts_zero@cnLineLabelsOn               = False      ; Turn off line labels
    opts_zero@tiMainString                 = ""
    opts_zero@gsnLeftString                = ""
    opts_zero@gsnRightString               = ""

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
    panel_strings               = (/"Local radial wind tendency"/)

;=================
; Plot the data 
;=================

  ; Plotting the LHS of the radial momentum budget (Huang et al. 2018)
  ; LHS (1): d{u}/dt - Local radial wind tendency
    output = "$nep/nepartak/images/cyl_coords/mbgt_rad_cxz_"+typ0+\
             "_sc"+sc+"_"+hy3+"km_"+w0+"w_mean"+mean0
    wks    = gsn_open_wks(opt,output)

  ; Local tendency of the mean radial wind [U_t]
    opts_mbgt@vpXF       = 0.10
    opts_mbgt@vpYF       = 0.90
    opts_mbgt@vpWidthF   = 0.70
    opts_mbgt@vpHeightF  = 0.70
    panel7      = gsn_csm_contour(wks,lhs1_c(:,:),opts_mbgt) ; [U_t]
    cont0       = gsn_csm_contour(wks,lhs1_c(:,:),opts_zero) ; Draw contours

    overlay(panel7, cont0)                                   ; Overlay contours
    txid(0)    = gsn_create_text(wks, panel_strings(0), txres)
    amid(0)    = gsn_add_annotation(panel7, txid(0), amres)
    delete([/cont0/])

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    pres = True 
    maximize_output(wks,pres)
  ; Call 'draw' and 'frame'

  end if
; End IF statement (choose whether to plot in addition to writing out to text files)

end 