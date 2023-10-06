; Produce 4-panel plots of layer-averaged diagnostics following a TC 

; Run using:

; ncl dat=\"03T00\" opt=\"x11\" ens0=\"em04\" dist=1.2 ts0=78 tf0=102 cn0=\"geo_sm\"
; sc=5 sm=0 ar=1.0 nr=21 nt=73 vlev0=8 vlev1=14 wlev0=17 wlev1=28 w0=0.3 plt=1 tc_erc_xy.ncl

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'opt'   = output file format ("pdf" or "x11")
; 'ens0'  = ensemble member (em00 -- em11)
; 'dist'  = size of box following storm (degrees)
; 'ts0'   = starting time index (0-119)
; 'tf0'   = ending time index (0-119)
; 'clr'   = colour map for plots ("ncl_default" [1];  "bdr_extra" [2]; "amwg256" [3])
; 'cn0'   = method of storm tracking ["slp", "geo", "geo_sm", "vort"]
; 'res'   = plot residual (1) or no residual (0) 
; 'w0'    = vertical velocity threshold for re-calculating storm centre position (0.1, 0.3, etc)
; 'vlev0' = lower level for tangential wind calculation
; 'vlev1' = upper level for tangential wind calculation
; 'wlev0' = lower level for vertical velocity calculation
; 'wlev1' = upper level for vertical velocity calculation 

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

; Find 'pc' and 'pd' files (model height level data)
  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0
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

; Calculate storm motion using built-in function ('st_rm')
; 'dat'     = initialisation time (02T12, 03T00, ...)
; 'ens0'    = ensemble simulation (em00, em01, ...)
; 'diri'    = path to input files (see above)
; 'dist'    = size of box following storm (degrees)
; 'mins'    = analyse 1-h (0) or 5-min (1) data

; Call external function 'st_rm'
  storm_rel = st_rm(dat, ens0, "$ar/text/", dist, 0)
  u_cyc     = storm_rel[0]		   ; Zonal wind 
  v_cyc     = storm_rel[1]		   ; Meridional wind
  vel_cyc   = storm_rel[2]		   ; Vector wind 
  numTIMES  = storm_rel[3]                 ; Number of times in file
  lat_arr   = storm_rel[4]                 ; Storm-following domain subset (lat)
  lon_arr   = storm_rel[5]                 ; Storm-following domain subset (lon)
  centre    = storm_rel[6]		   ; Storm centre at each hour (lat/lon)

;==========================================================================
; Create date/time string arrays for all times in file (use 'pb' stream)
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
  do ct = 0, numINPUT_c-1
   if (minute(ct).gt.30) then
     hour(ct) = hour(ct)+1
   end if
  end do
  delete(ct)

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

;========================================
; Define cylindrical coordinate arrays
;========================================

  pii        = 3.14159265
  radii      = fspan(0.0,ar,nr)      ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,nt)   ; Array of azimuth angles (0, 5, 10, etc)
  thetas_rad = (pii/180)*thetas      ; Degrees to radians
  rad_size   = dimsizes(radii)-1     ; Index for later calculation
  num_l      = 59      	   	     ; Number of model levels

; Array to hold all 'xpos' and 'ypos' position arrays
  xpos_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  ypos_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)

; Arrays to hold lat/lon at centre of cylindrical grid
  xcen_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)
  ycen_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),numTIMES/),float)

;=========================
; Start multiple loops
;=========================

   print_clock("Working on ensemble member: "+ens0)

   numTIMES    = 120                       ; Total number of times in all files [36]
   ct          = ts0                       ; Counter variable
   llbox       = toint(dist*50)  	   ; Size of lat/lon arrays (see line below)
   dsize       = (/llbox,llbox/)           ; Size of lat/lon arrays (depends on 'dist')
   d0	       = dsize(0)
   d1	       = dsize(1)

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

;=========================================================
; Read in required variables from 'pc' and 'pd' streams
;=========================================================

 ; Create array to hold storm track positions using different methods
   centre_new  = new((/7, num_l, numINPUT_d, 2/),float)
   icp_arr     = new((/numINPUT_d/),float)

   c  = addfile(fili_c,"r")         	 ; Read in 'pc' stream [u,v,vort]
   d  = addfile(fili_d,"r")         	 ; Read in 'pd' stream [w,p,t,theta]

;==================================
; Get the variables we will need
;==================================

   do it = ts0, tf0		; Loop over times in file

    print("it = "+it+" ; ct = "+ct)

  ; 'time_c' array (121) chosen, so 'it' time index used 
    print_clock("Working on time: "+time_str(it)+" (T+"+it+")" )
    title_arr(it) = "Valid at "+time_str(it)+" (T+"+it+")"

;=====================================================
; Read in variables from 'pc' stream (model levels)
;=====================================================

    t0    = lat_arr(cn,it-1,0)
    t1    = lat_arr(cn,it-1,1)
    n0    = lon_arr(cn,it-1,0)
    n1    = lon_arr(cn,it-1,1)
    print_clock("Grid (pc stream) extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

  ; Combine values above into array for input into external function
    ll_arr = (/t0,t1,n0,n1/)

  ; Read in basic variables from 'pc' stream
    lon1C  = c->longitude_1({n0:n1})              ; longitude ['d1' grid points]
    lat1C  = c->latitude_1({t0:t1})               ; latitude  ['d0' grid points]
    lonC   = c->longitude({n0:n1})                ; longitude ['d1' grid points]
    latC   = c->latitude({t0:t1})                 ; latitude  ['d0' grid points]
    hybC   = c->hybrid_ht(:) 		          ; 63 model ('theta') levels

    u_varname = "x-wind"
    v_varname = "y-wind"

  ; Read in horizontal wind components 
  ; Both variables below -- [hybC | 63] * [latC | 150] * [lonC | 150]
    u0  = c->$u_varname$(it,vlev0:vlev1,{t0:t1},{n0:n1}) ; Zonal wind	
    v0  = c->$v_varname$(it,vlev0:vlev1,{t0:t1},{n0:n1}) ; Meridional wind 

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

;======================================================
; Read in variables from 'pd' stream (model levels) 
;======================================================

  ; Read in basic variables from 'pd' stream 
    lon1D  = d->longitude_1({n0:n1})                ; longitude [xxx grid points]
    lat1D  = d->latitude_1({t0:t1})                 ; latitude  [xxx grid points]
    lonD   = d->longitude({n0:n1})                  ; longitude [xxx grid points]
    latD   = d->latitude({t0:t1})                   ; latitude  [xxx grid points]
    hybD   = d->hybrid_ht(:)			    ; 63 model ('rho') levels

  ; Model levels in km   
    hy     = hybD / 1000
    hy0    = sprintf("%0.0f",hybD(wlev0))
    hy1    = sprintf("%0.0f",hybD(wlev1))

  ; Read in temperature, pressure, vertical velocity and geopotential height
  ; All variables below -- [hybD | 63] * [latD | 150] * [lon1D | 151]
    t      = d->temp(it-2,wlev0:wlev1,{t0:t1},{n0:n1})        ; Temperature (K)
    p      = d->p(it-2,wlev0:wlev1,{t0:t1},{n0:n1})           ; Pressure (Pa)
    w      = d->dz_dt(it-2,wlev0:wlev1,{t0:t1},{n0:n1})       ; Vertical velocity (m s-1)
    z      = d->ht(it-2,wlev0:wlev1,{t0:t1},{n0:n1})          ; Geopotential height (m)
    pv     = d->field83(it-2,wlev0:wlev1,{t0:t1},{n0:n1})     ; Potential vorticity (PVU)

  ; Convert units before continuing 
    p  	   = p / 100           ; Convert to 'hPa' 
    z  	   = z / 10            ; Convert to 'dam'
    pv     = pv * (10 ^ 6)     ; Convert to 'PVU' [10 ^ 6 K m2 s-1 kg-1] 

  ; Add metadata
    p@units  = "hPa"
    z@units  = "dam"
    pv@units = "PVU"

    dx0      = lonD(1) - lonD(0) ; Grid spacing (longitude)
    dy0      = latD(1) - latD(0) ; Grid spacing (latitude)

;=========================================
; Calculate layer-averaged diagnostics
;=========================================

     print_clock("Working on layer between " + hy0+ " and "+hy1+"m AGL")

   ; Read in other variables on model levels (u,v,w,f,p,str,rho,vort,avo)
     u_plane 	            = dim_avg_n_Wrap(u(:,:,:),0)
     v_plane		    = dim_avg_n_Wrap(v(:,:,:),0)
     w_plane		    = dim_avg_n_Wrap(w(:,:,:),0)
     f_plane		    = dim_avg_n_Wrap(f0(:,:,:),0)
     prs_plane		    = dim_avg_n_Wrap(p(:,:,:),0)
     str_plane              = dim_avg_n_Wrap(str0(:,:,:),0)
     vort_plane	   	    = dim_avg_n_Wrap(vort(:,:,:),0)
     avo_plane		    = dim_avg_n_Wrap(avo(:,:,:),0)
     pv_plane               = dim_avg_n_Wrap(pv(:,:,:),0)

   ; Calculate vector windspeed and smooth
     spd_plane              = sqrt( (u_plane ^ 2) + (v_plane ^ 2) )
     copy_VarCoords(u_plane, spd_plane)
     vort_pl_smth           = smth9_Wrap(vort_plane, 0.5, 0.5, True)

   ; Combine all arrays above into single, larger array for input into 'setup_cyl'
     pc_plane               = (/u_plane, v_plane, spd_plane, vort_plane, \
                                avo_plane, f_plane, str_plane, str_plane/)
     
   ; NOTE: second 'pv_plane' is a placeholder 
     pd_plane               = (/w_plane, prs_plane, pv_plane, pv_plane/)

;=======================================================
; Calculate storm centre position (external function)  
;======================================================= 

   ; Set radial distance (º) when looking for speed/vort min (in ringlike phase)
     r0            = 0.12

   ; Set 'ilev' equal to zero ahead of function call
     ilev          = 0 

     print_clock("Calculating storm centre position...")
     centre_arr    = st_centre(ilev, it, cn, r0, centre, vort_plane, prs_plane, spd_plane, \
                               centre_new, dy0, dx0, "single", ct, "pd", "off", "1h")

   ; Updated storm track information contained in 'centre_new' 
     centre_new     = centre_arr[0]

     print("Storm centre = "+centre_new(sc,ilev,it-2,0)+" ; "+centre_new(sc,ilev,it-2,1) )

;=====================================================================
; Make sure that the storm centre is not within the eyewall updraft 
;=====================================================================

   ; Retrieve coordinates of storm centre
     cen_loc0 = centre_new(sc,ilev,it-2,0)
     cen_loc1 = centre_new(sc,ilev,it-2,1)

   ; Coordinates of grid points either side  
     cen_0n = cen_loc0 - dy0
     cen_0x = cen_loc0 + dy0
     cen_1n = cen_loc1 - dy0
     cen_1x = cen_loc1 + dy0

   ; Calculate vertical velocity at these grid points 
     cen_w    = w_plane( {cen_0n:cen_0x}, {cen_1n:cen_1x} )

   ; If average over several grid points is > threshold value, search for centre again
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

    ; Replace values in 'centre_new' before calling 'setup_cyl' below  
      centre_new(sc,ilev,it-2,0) = lat_min_w
      centre_new(sc,ilev,it-2,1) = lon_min_w

     end if 

   ; Tidy up 
     delete([/cen_w, cen_0n, cen_0x, cen_1n, cen_1x/])

;========================================================
; Switch to cylindrical coordinates (external function) 
;========================================================

   ; Set 'ilev' option before function call 
     ilev = 0
 
   ; Call 'setup_cyl' and output pressure/vorticity on cylindrical grid 
   ; Second to last option = 0 --> basic variables output 
   ; Second to last option = 1 --> additional MBGT variables output
   ; Last option --> composite ("comp") or single simulation ("single"); same as 'st_centre' call
   ; 13/01/2020 --> when typ = "single", next option is redundant in this call (ct = it-2)
     print_clock("Translating to cylindrical grid...")
     cyl_arr  = setup_cyl(it, dist, ll_arr, lonC, latC, lonD, latD, lon1D, lat1D, \
                          pc_plane, pd_plane, numTIMES, centre_new, \
                          sc, sm, ar, nr, nt, ilev, num_l, 0, "single", it-2, \
                          0, xpos_all, ypos_all, xcen_all, ycen_all)

     xpos_all = cyl_arr[0]
     ypos_all = cyl_arr[1]
     xcen_all = cyl_arr[2]
     ycen_all = cyl_arr[3]
     lat_max  = cyl_arr[4]
     lon_max  = cyl_arr[5]
     radii    = cyl_arr[6]
     rad_size = cyl_arr[7]

     u_int    = cyl_arr[8]
     v_int    = cyl_arr[9]
     spd_int  = cyl_arr[10]
     vort_int = cyl_arr[11]
     avo_int  = cyl_arr[12]
     vr_int   = cyl_arr[13]
     vt_int   = cyl_arr[14]
     w_int    = cyl_arr[15]
     prs_int  = cyl_arr[16]
     pv_int   = cyl_arr[23]
     slp_int  = cyl_arr[25]

     print("Storm centre = "+lat_max+" ; "+lon_max)

;==========================================
; Produce 4-panel plot of key variables
;==========================================

   ; Plot the data as we go
     output_loop = "$nep/nepartak/images/cyl_coords/xy_4panel_"+dat+"_"+ens0+\
                   "_"+hy0+"m_"+hy1+"m_sc"+sc+"_"+time_arr(it)
     wks = gsn_open_wks(opt, output_loop)

     gsn_define_colormap(wks,"prcp_new")

   ; Relative vorticity resources
     opts_rv                               = True
     opts_rv@cnFillOn                      = True
     opts_rv@sfXArray                      = xpos_all(ilev,:,:,it-2)
     opts_rv@sfYArray                      = ypos_all(ilev,:,:,it-2)
     opts_rv@cnLineLabelInterval           = 2.0
     opts_rv@cnLineLabelFontHeightF        = 0.012
     opts_rv@cnLineLabelBackgroundColor    = "transparent"
     opts_rv@cnLineLabelPlacementMode      = "constant"
     opts_rv@cnLinesOn                     = False
     opts_rv@cnInfoLabelOn                 = False
     opts_rv@cnLineLabelsOn                = False
     opts_rv@cnLevelSelectionMode          = "ExplicitLevels"
     opts_rv@cnLevels                      = (/2., 3., 5., 10., 15., 20., 25., 30., \
                                               40., 50., 60., 70., 80./)
     opts_rv@cnFillColors                  = (/0,3,4,5,6,7,8,9,10,11,12,13,14,15/)
     opts_rv@gsnMaximize                   = False

   ; No title or additional strings
     opts_rv@tiMainString                  = ""
     opts_rv@gsnLeftString                 = ""
     opts_rv@gsnRightString                = ""

   ; Turn individual label bar on/off
     opts_rv@lbLabelBarOn                   = True
     opts_rv@lbBoxEndCapStyle               = "TriangleBothEnds"
     opts_rv@lbOrientation                  = "Vertical"

     opts_rv@trYMaxF                        = lat_max + ar
     opts_rv@trYMinF                        = lat_max - ar
     opts_rv@trXMaxF                        = lon_max + ar
     opts_rv@trXMinF                        = lon_max - ar

   ; Do not draw the plot or advance the frame                                        
     opts_rv@gsnDraw                       = False
     opts_rv@gsnFrame                      = False
     opts_rv@gsnAddCyclic                  = False

   ; Background plot resources                                                        
     opts_xy                               = True

   ; No tickmarks on any axis                                                         
     opts_xy@tmYROn                        = False
     opts_xy@tmYLOn                        = True
     opts_xy@tmXTOn                        = False
     opts_xy@tmXBOn                        = True

     opts_xy@xyLineColor                   = -1
     opts_xy@gsnDraw                       = False
     opts_xy@gsnFrame                      = False

     opts_xy@trYMaxF                       = lat_max + ar
     opts_xy@trYMinF                       = lat_max - ar
     opts_xy@trXMaxF                       = lon_max + ar
     opts_xy@trXMinF                       = lon_max - ar
     opts_xy@gsnMaximize                   = False

   ; Resources for windspeed
     opts_wind                             = opts_rv
     delete([/opts_wind@cnLevels, opts_wind@cnFillColors/])
     opts_wind@cnFillOn                    = True
     opts_wind@cnLevelSelectionMode        = "ExplicitLevels"
     opts_wind@cnFillPalette               = "radar_new"
     opts_wind@cnLevels                    = (/5.0, 10.0, 15.0, 20.0,\
                                               25.0, 30.0, 35.0, 40.0, \
                                               45.0, 50.0, 55.0, 60.0, 65.0/)
     opts_wind@cnFillColors                = (/-1,0,1,2,3,4,5,6,7,8,9,10,11,12/)

   ; Resources for vertical velocity                                                             
     opts_vvel                             = opts_wind
     delete([/opts_vvel@cnLevels, opts_vvel@cnFillColors, opts_vvel@cnFillPalette/])

     opts_vvel@cnFillPalette               = "BlueDarkRed18"
     opts_vvel@cnLevels                    = (/-3.0, -2.0, -1.5, -1.0, -0.5, -0.2, -0.1, -0.05,\
                                                0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0/)
     opts_vvel@cnFillColors                = (/1,2,3,4,5,6,7,8,-1,9,10,11,12,13,14,15,16/)

   ; Resources for strain rate                                                                   
     opts_strn                             = opts_vvel
     delete([/opts_strn@cnLevels, opts_strn@cnFillColors, opts_strn@cnFillPalette/])
     opts_strn@cnFillPalette               = "circ_budget"
     opts_strn@cnLevels                    = (/10,12,14,16,18,20,22,25,30/)
     opts_strn@cnFillColors                = (/-1,11,13,15,16,17,18,19,20,21/)

   ; Resources for divergence                                                                    
     opts_divg                             = opts_vvel
     delete([/opts_divg@cnLevels, opts_divg@cnFillColors, opts_divg@cnFillPalette/])
     opts_divg@cnFillPalette               = "BlueDarkRed18"
     opts_divg@cnLevels                    = (/-100, -50, -30, -20, -10, -5, \
                                               5, 10, 20, 30, 50, 100/)
     opts_divg@cnFillColors                = (/1,2,3,4,6,8,-1,9,11,13,14,15,16/)

   ; Resources for potential vorticity 
     opts_pvrt                             = opts_rv
     delete([/opts_pvrt@cnLevels, opts_pvrt@cnFillColors/])
     opts_pvrt@cnLevels                    = (/9.0, 10.0 ,12.0, 14.0, 16.0, 18.0, 20.0, 22.0, \
                                               25.0, 30.0, 35.0, 40.0, 45.0, 50.0/)
     opts_pvrt@cnFillColors                = (/0,2,3,4,5,6,7,8,9,10,\
                                               10,11,12,13,14,15/)

   ; Resources for radial wind                                                                   
     opts_vrad                             = opts_vvel
     delete([/opts_vrad@cnLevels, opts_vrad@cnFillColors, opts_vrad@cnFillPalette/])
     opts_vrad@cnFillPalette               = "BlueDarkRed18"
     opts_vrad@cnLevels                    = (/-20., -15., -10., -5., -3., -2., -1., -0.5, \
                                                0, 0.5, 1., 2., 3., 5., 10., 15., 20./)
     opts_vrad@cnFillColors                = (/0,2,3,4,5,6,7,8,-1,-1,10,11,12,13,14,15,16,17/)

   ; Panel plot label resources                                                                  
     txid                                  = new(4,graphic)
     amid                                  = new(4,graphic)

     panel_strings                         = (/ "a) Relative vorticity",\
                                                "b) Vertical velocity",\
                                                "c) Tangential wind",\
                                                "d) Radial wind"/)

     txres                                  = True
     txres@txPerimOn                        = True
     txres@txFontHeightF                    = 0.0125
     txres@txBackgroundFillColor            = "White"

     amres                                  = True
     amres@amParallelPosF                   = -0.48
     amres@amOrthogonalPosF                 = -0.48
     amres@amJust                           = "TopLeft"

   ; Coordinates of panels                                                              
     xf = (/0.05, 0.57, 0.05, 0.57/)
     yf = (/0.90, 0.90, 0.52, 0.52/)
     wf = 0.35
     hf = 0.32

   ; Panel 1 (relative vorticity)                                                       
     opts_rv@vpXF       = xf(0)
     opts_rv@vpYF       = yf(0)
     opts_rv@vpWidthF   = wf
     opts_rv@vpHeightF  = hf

     opts_xy0           = opts_xy
     opts_xy0@vpXF      = xf(0)
     opts_xy0@vpYF      = yf(0)
     opts_xy0@vpWidthF  = wf
     opts_xy0@vpHeightF = hf

     vort_int_plot      = vort_int * (10 ^ 4)

     plot_tl            = gsn_csm_xy(wks, xpos_all(ilev,:,:,it-2), \
                                     ypos_all(ilev,:,:,it-2), opts_xy0)
     plot_tl0           = gsn_csm_contour(wks, vort_int_plot(:,:), opts_rv)
     overlay(plot_tl, plot_tl0)

   ; Overlay labels 
     txid(0)      = gsn_create_text(wks, panel_strings(0), txres)
     amid(0)      = gsn_add_annotation(plot_tl, txid(0), amres)

   ; Panel 2 (vertical velocity)
     opts_vvel@vpXF       = xf(1)
     opts_vvel@vpYF       = yf(1)
     opts_vvel@vpWidthF   = wf
     opts_vvel@vpHeightF  = hf

     opts_xy1             = opts_xy
     opts_xy1@vpXF        = xf(1)
     opts_xy1@vpYF        = yf(1)
     opts_xy1@vpWidthF    = wf
     opts_xy1@vpHeightF   = hf

     plot_tr              = gsn_csm_xy(wks, xpos_all(ilev,:,:,it-2), \
                                       ypos_all(ilev,:,:,it-2), opts_xy1)
     plot_tr0             = gsn_csm_contour(wks, w_int(:,:), opts_vvel)
     overlay(plot_tr, plot_tr0)

   ; Overlay labels 
     txid(1)      = gsn_create_text(wks, panel_strings(1), txres)
     amid(1)      = gsn_add_annotation(plot_tr, txid(1), amres)

   ; Panel 3 (tangential wind)
     opts_wind@vpXF       = xf(2)
     opts_wind@vpYF       = yf(2)
     opts_wind@vpWidthF   = wf
     opts_wind@vpHeightF  = hf

     opts_xy2             = opts_xy
     opts_xy2@vpXF        = xf(2)
     opts_xy2@vpYF        = yf(2)
     opts_xy2@vpWidthF    = wf
     opts_xy2@vpHeightF   = hf

     plot_bl              = gsn_csm_xy(wks, xpos_all(ilev,:,:,it-2), \
                                       ypos_all(ilev,:,:,it-2), opts_xy2)
     plot_bl0             = gsn_csm_contour(wks, vt_int(:,:), opts_wind)
     overlay(plot_bl, plot_bl0)

   ; Overlay labels 
     txid(2)      = gsn_create_text(wks, panel_strings(2), txres)
     amid(2)      = gsn_add_annotation(plot_bl, txid(2), amres)

   ; Panel 4 (radial wind)                                                              
     opts_vrad0           = opts_vrad
     opts_vrad0@vpXF      = xf(3)
     opts_vrad0@vpYF      = yf(3)
     opts_vrad0@vpWidthF  = wf
     opts_vrad0@vpHeightF = hf

     opts_xy3             = opts_xy
     opts_xy3@vpXF        = xf(3)
     opts_xy3@vpYF        = yf(3)
     opts_xy3@vpWidthF    = wf
     opts_xy3@vpHeightF   = hf

     plot_br              = gsn_csm_xy(wks, xpos_all(ilev,:,:,it-2), \
                                       ypos_all(ilev,:,:,it-2), opts_xy3)
     plot_br0             = gsn_csm_contour(wks, vr_int(:,:), opts_vrad0)
     overlay(plot_br, plot_br0)

   ; Overlay labels                                                                     
     txid(3)      = gsn_create_text(wks, panel_strings(3), txres)
     amid(3)      = gsn_add_annotation(plot_br, txid(3), amres)

   ; Tidy up before final plots later in script 
     delete([/opts_vvel, opts_wind, opts_rv, opts_strn, opts_divg, opts_vrad, \
              opts_xy, opts_xy0, opts_xy1, opts_xy2, opts_xy3,\
              txid, amid, panel_strings/])

;==================================================
; Overlay markers for each storm centre estimate
;================================================== 

; (0) vort, (1) pressure, (2) pressure_smth, (3), pressure_smth2, (4) vort_min,
; (5) speed_min, (6) offline 
  m_colours  = (/"royalblue1", "orange4", "orange4", "orange4", "royalblue1", \
                 "black", "red4"/)
  m_indices  = (/12, 12, 16, 8, 8, 8, 12/)
  m_sizes    = (/10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0/)
  m_thick    = (/5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0/)

  markers    = new( (/4, dimsizes(m_colours)/), graphic)

  do i = 0, dimsizes(m_colours)-1
    mres0                    = True
    mres0@gsMarkerIndex      = m_indices(i)
    mres0@gsMarkerSizeF      = m_sizes(i)
    mres0@gsMarkerColor      = m_colours(i)
    mres0@gsMarkerThicknessF = m_thick(i)
    markers(0,i)             = gsn_add_polymarker(wks, plot_tl, centre_new(i,ilev,it-2,1),\
                                                  centre_new(i,ilev,it-2,0), mres0)
    markers(1,i)             = gsn_add_polymarker(wks, plot_tr, centre_new(i,ilev,it-2,1),\
                                                  centre_new(i,ilev,it-2,0), mres0)
    markers(2,i)             = gsn_add_polymarker(wks, plot_bl, centre_new(i,ilev,it-2,1),\
                                                  centre_new(i,ilev,it-2,0), mres0)
    markers(3,i)             = gsn_add_polymarker(wks, plot_br, centre_new(i,ilev,it-2,1),\
                                                  centre_new(i,ilev,it-2,0), mres0)
  end do 

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

    ; Create a main title at the top of the page
      title_res                   = True
      title_str                   = "Panel plot of "+hy0+" to "+hy1+" m AGL: "+ens0+\
                                    ", "+dat+" (T+"+it+")"
      title_res@txFontHeightF     = 0.02
      title_res@txFont            = "Helvetica"
      xpos_title                  = 0.50
      ypos_title                  = 0.94

      gsn_text_ndc(wks, title_str, xpos_title, ypos_title, title_res)

    ; Call 'draw' and 'frame'
      maximize_output(wks, True)

   ; Tidy up
     delete([/u_int,v_int,w_int,vt_int,vr_int/])
     delete([/vort_int,avo_int/])
     delete([/u_plane,v_plane,w_plane,f_plane,vort_plane,avo_plane/])

   ct  = ct + 1		    ; Counter variable (time)

 ; Tidy up before next iteration (all variables from 'pc' and 'pd' streams)
   delete([/lonC,latC,lon1C,lat1C,u,v,vort,avo,f0,f1/])
   delete([/lonD,latD,lon1D,lat1D,t,p,w,z/])

 end do     ; End time loop (do it = 0, times-1)

end