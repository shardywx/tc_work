; Script to identify vacillation cycles (VCs) in UM ensemble simulations

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_centre.ncl"
load "$sam/ncl_func/setup_cyl.ncl"
load "$sam/ncl_func/ring_mono.ncl"

; ncl opt=\"x11\" dat=\"03T00\" ens0=\"em04\" mlev0=13 dist=1.0 r0=0.50 cn0=\"geo_sm\" ts=86
; tf=91 ti=1 var0=\"rv\" full=0 tc=0 lay=1 mlev1=25 sc=5 sm=0 paper=1 diff=1 tc_vort_str_paper.ncl

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'ens0'  = ensemble simulation (em00-em11)
; 'dist'  = size of box following storm (degrees)
; 'opt'   = output file format ("pdf" or "x11")
; 'mlev0' = model level for plotting PV and vertical velocity (see L106-108)
; 'r0'    = radius for final plot (degrees)
; 'var0'  = plot either PV ("pv") or relative vorticity ("rv")
; 'ts'	  = starting time index (0-116)
; 'tf'	  = ending time index (0-116)
; 'ti'	  = time interval (1,2,3,4,etc)
; 'full'  = produce 3-panel (1) or 2-panel (0) plot 
; 'tc'	  = annotate images with TC symbol (1) or leave as they were (0)
; 'lay'	  = calculate layer average (1) 
; 'mlev1' = upper model level for calculating layer average (only valid when lay=1)
; 'paper' = settings to produce plot for December 2019 version of VC paper (1)

begin

;==============================
; Prelash (reading in files)
;==============================

; Troubleshooting
  if (r0 .gt. dist) then 
   print("Choose a bigger domain ('dist') or smaller plot size ('r0')!")
   exit()
  end if  

; Find 'pc' and 'pd' files (model height level data)
  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0
  fili_c1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pc*.nc")
  fili_c      = diri+"/"+fili_c1
  fili_d1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pd*.nc")
  fili_d      = diri+"/"+fili_d1

; Calculate number of times in each file using built-in function ('nc_times')
  numINFO_c   = nc_times(fili_c)
  numINFO_d   = nc_times(fili_d)
  numINPUT_c  = numINFO_c[0]            ; Number of times in 'pc' stream
  print("pc stream files contain "+numINFO_c[0]+" times")
  time_c      = numINFO_c[1]            ; Array of times from 'pc' stream file
  numINPUT_d  = numINFO_d[0]            ; Number of times in 'pd' stream
  time_d      = numINFO_d[1]            ; Array of times from 'pd' stream file
  print("pd stream files contain "+numINFO_d[0]+" times")

; Calculate storm motion using built-in function ('st_rm')
; 'dat'	    = initialisation time (02T12, 03T00, etc)
; 'ens0'    = ensemble simulation (em00-em11)
; 'diri'    = path to input files (see above)
; 'dist'    = size of box following storm (degrees)
; 'mins'    = analyse 1-h (0) or 5-min (1) data

  storm_rel = st_rm(dat, ens0, "$ar/text/", dist, 0)
  u_cyc     = storm_rel[0]             ; Zonal wind
  v_cyc     = storm_rel[1]             ; Meridional wind
  vel_cyc   = storm_rel[2]             ; Vector wind
  numTIMES  = storm_rel[3]             ; Number of times in file
  lat_arr   = storm_rel[4]             ; Latitude array subset (following storm)
  lon_arr   = storm_rel[5]             ; Longitude array subset (following storm)
  centre    = storm_rel[6]             ; Storm track (position) array 

; Create arrays for data and output information
  llbox	    = toint(dist*50)		       ; Calculate domain size (grid points)
  dsize	    = (/llbox,llbox/)		       ; Domain size (grid points)

; Extract string for each ensemble member
  str1    = str_split(fili_c(0),"_")
  ens_str = str1(5)
  delete(str1)

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

  date_str = sprinti("%0.2iUTC ", hour) + \
             sprinti("%0.2i ", day) \
           + month_abbr(month)
  time_str = sprinti("%0.2i UTC ", hour) + \
             sprinti("%0.2i ", day) \
           + month_abbr(month)
  time_arr = sprinti("%0.2i", day) + \
             month_abbr(month) + \
             "_" + sprinti("%0.2iZ", hour)

  print("Working on ensemble member "+ens_str)

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 10000000000
  end setvalues

;==================================
; Get the variables we will need 
;==================================

; Create array to hold storm track positions using different methods 
  num_l       = 59
  centre_new  = new((/8, num_l, numINPUT_c, 2/),float)

  c = addfile(fili_c,"r")	   ; Read in 'pc' stream [u,v,vort]
  d = addfile(fili_d,"r")  	   ; Read in 'pd' stream [w,p,t,theta,PV]

  do it = ts, tf, ti		   ; Loop over times in file 

    time  = c->t(it)                ; Time in file
    times = dimsizes(time)

    title_arr(it) = "valid at "+time_str(it)+" (T+"+it+")"
    print("Working on "+time_str(it)+" (T+"+it+")")

  ; Get domain boundaries at each time from storm track (calculated offline)
    t0    = lat_arr(cn,it-1,0)
    t1    = lat_arr(cn,it-1,1)
    n0    = lon_arr(cn,it-1,0)
    n1    = lon_arr(cn,it-1,1)

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
    if (lay .eq. 1) then 
     u0 = c->$u_varname$(it,mlev0:mlev1,{t0:t1},{n0:n1})     ; Zonal wind (layer)
     v0 = c->$v_varname$(it,mlev0:mlev1,{t0:t1},{n0:n1})     ; Meridional wind (layer)
    else 
     u0 = c->$u_varname$(it,mlev0,{t0:t1},{n0:n1})	    ; Zonal wind (single level)
     v0 = c->$v_varname$(it,mlev0,{t0:t1},{n0:n1})     	    ; Meridional wind (single level)
    end if 

 ; Grid spacing  
   dy   = latC(1) - latC(0)
   dx   = lonC(1) - lonC(0)

 ; Calculate storm-relative winds
   u = u0 - u_cyc(cn,it-1)
   v = v0 - v_cyc(cn,it-1)

 ; Add metadata from original horizontal wind arrays
   copy_VarCoords(u0,u)
   copy_VarAtts(u0,u)

   copy_VarCoords(v0,v)
   copy_VarAtts(v0,v)

 ; Calculate relative vorticity using centered finite differences
   vort = uv2vr_cfd(u, v, latC, lonC, 2)
   copy_VarCoords(u,vort)
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

;==================================================================
; Read in variables from 'pd' stream (model levels)
; B/c of differences in file size, access variables using 'it+2'
;==================================================================

  ; Read in basic variables from 'pd' stream
    lon1D  = d->longitude_1({n0:n1})                ; longitude [xxx grid points]
    lat1D  = d->latitude_1({t0:t1})                 ; latitude  [xxx grid points]
    lonD   = d->longitude({n0:n1})                  ; longitude [xxx grid points]
    latD   = d->latitude({t0:t1})                   ; latitude  [xxx grid points]
    hybD   = d->hybrid_ht(:)                        ; 63 model ('rho') levels

    hy     = hybD / 1000 ; Model level values in km
    hy0	   = sprintf("%0.0f",hybD(mlev0))
    hy1	   = sprintf("%0.0f",hybD(mlev1))
    print("Working on model level: "+hy0+ " m")

  ; Read in desired variables 
    if (lay .eq. 1) then 

     print("Averaging PV between "+hy0+" and "+hy1+ " m")
     pv	   = d->field83(it-2,mlev0-5:mlev0+5,{t0:t1},{n0:n1}) ; Potential vorticity (PVU)
     w	   = d->dz_dt(it-2,mlev0-5:mlev0+5,{t0:t1},{n0:n1})   ; Vertical velocity (m s-1; lower)
     w2    = d->dz_dt(it-2,mlev1-5:mlev1+5,{t0:t1},{n0:n1})   ; Vertical velocity (m s-1; upper)
     t     = d->temp(it-2,mlev0-5:mlev0+5,{t0:t1},{n0:n1})    ; Temperature (K) 
     p     = d->p(it-2,mlev0-5:mlev0+5,{t0:t1},{n0:n1})       ; Pressure (Pa) 
     z     = d->ht(it-2,mlev0-5:mlev0+5,{t0:t1},{n0:n1})      ; Geopotential height (m)

    else

     pv    = d->field83(it-2,mlev0,{t0:t1},{n0:n1})           ; Potential vorticity (PVU)
     w     = d->dz_dt(it-2,mlev0,{t0:t1},{n0:n1})             ; Vertical velocity (m s-1; lower)
     w2	   = d->dz_dt(it-2,mlev1,{t0:t1},{n0:n1})             ; Vertical velocity (m s-1; upper)
     t     = d->temp(it-2,mlev0,{t0:t1},{n0:n1})              ; Temperature (K) 
     p     = d->p(it-2,mlev0,{t0:t1},{n0:n1})                 ; Pressure (Pa) 
     z     = d->ht(it-2,mlev0,{t0:t1},{n0:n1})                ; Geopotential height (m)

    end if 

  ; Convert units of pressure and geopotential height before continuing
    p      = p / 100          ; Convert to 'hPa'
    z      = z / 10	      ; Convert to 'dam'
    pv     = pv * (10 ^ 6)    ; Convert to 'PVU' [10 ^ 6 K m2 s-1 kg-1]

  ; Add metadata
    p@units  = "hPa"
    z@units  = "dam"
    pv@units = "PVU"

;===========================================================
; Calculate rate of strain (following Nguyen et al. 2011)
;===========================================================

     dim1 = "latitude"
     dim2 = "longitude"

   ; Calculate horizontal derivatives of 'u' (du_dx, du_dy)
     du    = grad_latlon_cfd(u,u&$dim1$,u&$dim2$,False,False)
     du_dy = du[0]    ; Latitudinal derivative
     du_dx = du[1]    ; Longitudinal derivative
     delete(du)

   ; Calculate horizontal derivatives of 'v' (dv_dx, dv_dy)
     dv    = grad_latlon_cfd(v,v&$dim1$,v&$dim2$,False,False)
     dv_dy = dv[0]    ; Latitudinal derivative
     dv_dx = dv[1]    ; Longitudinal derivative
     delete(dv)

   ; Calculate strain --> sqrt[ (du_dx - dv_dy)^2 + (dv_dx + du_dy)^2 ]
     str1a   = du_dx - dv_dy   ; Stretching deformation
     str2a   = dv_dx + du_dy   ; Shearing deformation
     str1    = (str1a)^2
     str2    = (str2a)^2
     
     if (lay .eq. 1) then 
      str0  = sqrt( str1(:,:,:) + str2(:,:,:) )
      strT  = str1(:,:,:) + str2(:,:,:)
     else 
      str0  = sqrt( str1(:,:) + str2(:,:) )
      strT  = str1(:,:) + str2(:,:)
     end if

;===========================================
; Analyse variables on chosen model level
;===========================================

   ; Calculate single-level or layer-averaged variables 
     
     if (lay .eq. 1) then 

      print("Working on layer between " + hy0+ " and "+hy1+" m AGL")

    ; 'pc' stream
      u_plane    = dim_avg_n_Wrap(u(:,:,:),0)
      v_plane    = dim_avg_n_Wrap(v(:,:,:),0)
      vort_plane = dim_avg_n_Wrap(vort(:,:,:),0)
      f_plane    = dim_avg_n_Wrap(f0(:,:,:),0)

    ; 'pd' stream
      pv_plane   = dim_avg_n_Wrap(pv(:,:,:),0)
      w_plane	 = dim_avg_n_Wrap(w(:,:,:),0)
      w_plane2	 = dim_avg_n_Wrap(w2(:,:,:),0)
      prs_plane  = dim_avg_n_Wrap(p(:,:,:),0)
      str_plane	 = dim_avg_n_Wrap(str0(:,:,:),0)
      strT_plane = dim_avg_n_Wrap(strT(:,:,:),0)
      copy_VarCoords(u(0,:,:),str_plane)
      copy_VarCoords(u(0,:,:),strT_plane)

     else

      print("Working on model level: "+hy0+ " m")

    ; 'pc' stream
      u_plane       = u(:,:)
      v_plane       = v(:,:)
      vort_plane    = vort(:,:)
      f_plane       = f0(:,:)

    ; 'pd' stream
      pv_plane	    = pv(:,:)
      w_plane	    = w(:,:)
      w_plane2	    = w2(:,:)
      prs_plane     = p(:,:)
      str_plane	    = str0(:,:)
      strT_plane    = strT(:,:)
      copy_VarCoords(u,str_plane)
      copy_VarCoords(u,strT_plane)

     end if 

   ; Calculate vector windspeed and smooth  
     spd_plane     = sqrt( (u_plane ^ 2) + (v_plane ^ 2) )
     copy_VarCoords(u_plane, spd_plane)
     vort_pl_smth  = smth9_Wrap(vort_plane, 0.5, 0.5, True)

   ; Also calculate Laplacian of the horizontal wind speed 
     spd_grad = grad_latlon_cfd(spd_plane, latC, lonC, False, False)
     spd_grad_lat = spd_grad[0]
     spd_grad_lon = spd_grad[1]

     lap_plane  = uv2dv_cfd(spd_grad_lon, spd_grad_lat, latC, lonC, 2)
     copy_VarCoords(spd_plane, lap_plane)
     lap_plane  = lap_plane * (10 ^ 7)
     lap_plane@units = "10~S~-7~N~ m~S~-1~N~ s~S~-1~N~"
     lap_plane@name  = "Laplacian of horizontal wind speed"

   ; Combine all pc and pd stream variables into larger arrays
   ; First 'f_plane' is a placeholder (pc_plane)
   ; 'w_plane2' is a placeholder (pd_plane)
     pc_plane      = (/u_plane, v_plane, spd_plane, \
                       vort_plane, f_plane, f_plane, str_plane/)
     pd_plane      = (/w_plane, prs_plane, pv_plane/)

;=============================================================
; Call external function to calculate storm centre position
;=============================================================

 ; Set radial distance (º) when looking for vort min in ringlike phase 
   r1            = 0.10 
   centre_arr    = st_centre(0, it, cn, r1, centre, vort_plane, prs_plane, spd_plane, \
                             centre_new, dy, dx, "single", it, "pd", "off", "1h", lap_plane)

 ; Output centre positions 
   centre_new    = centre_arr[0]

 ; Tidy up storm centre for later plotting 
   lat_max       = centre_new(sc,0,it-2,0)
   lon_max       = centre_new(sc,0,it-2,1)

   print("Storm centre = "+lat_max+" ; "+lon_max)

;=================
; Plot the data 
;=================

   ; Choose size of panel plot 
     if (full .eq. 1) then
      panel = new(3,graphic)
     else
      panel = new(2,graphic)
     end if

   ; Set output file type and location
     if (lay .eq. 1) then 
      output = "$nep/nepartak/images/vc/"+ens0+"_"+dat+"_str_"+var0+\
               "_"+hy0+"m_"+hy1+"m_"+time_arr(it)
     elseif (lay .eq. 0) then 
      output = "$nep/nepartak/images/vc/"+ens0+"_"+dat+"_str_"+var0+\
      	       "_"+hy0+"m_"+time_arr(it)
     end if 

     wks = gsn_open_wks(opt,output)

   ; Load colour table
     gsn_define_colormap(wks,"prcp_new")

;==========================
; Options for plotting
;==========================

   ; Vorticity
     opts_pv                               = True
     opts_pv@cnFillOn                      = True
     opts_pv@cnLineLabelInterval           = 2.0
     opts_pv@cnLineLabelFontHeightF        = 0.012
     opts_pv@cnLineLabelBackgroundColor    = "transparent"
     opts_pv@cnLineLabelPlacementMode      = "constant"
     opts_pv@cnLinesOn                     = False
     opts_pv@cnInfoLabelOn                 = False
     opts_pv@cnLevelSelectionMode          = "ExplicitLevels"

     if (var0 .eq. "pv") then 

      opts_pv@cnLevels                     = (/9.0, 10.0 ,12.0, 14.0, \
                                               16.0, 18.0, 20.0, 22.0, 25.0, \
                                               30.0, 35.0, 40.0, 45.0, 50.0/)
      opts_pv@cnFillColors                 = (/0,2,3,4,5,\
                                               6,7,8,9,10,\
                                               11,12,13,14,15/)
     elseif (var0 .eq. "rv") then 

      opts_pv@cnLevels                     = (/5., 10., 15., 20., 25., 30., \
                                               40., 50., 60., 70., 80./)
      opts_pv@cnFillColors                 = (/0,4,5,6,7,8,9,10,11,12,13,14/)

     end if 

     opts_pv@tiMainString                  = ""
     opts_pv@tiMainFontHeightF             = 0.0125
     opts_pv@gsnLeftString                 = ""
     opts_pv@gsnRightString                = ""

   ; Labelbar on/off
     opts_pv@lbLabelBarOn                  = True
     opts_pv@lbOrientation                 = "vertical"

   ; Labelbar width; height; orthogonal position; (left/right); parallel position (moves up/down)
     opts_pv@pmLabelBarWidthF              = 0.08
     opts_pv@pmLabelBarHeightF             = 0.48
     opts_pv@pmLabelBarOrthogonalPosF      = 0.10
     opts_pv@pmLabelBarParallelPosF        = 0.60

   ; Additional lb resources (perimeter; label font; box end cap style)
     opts_pv@lbPerimOn                     = False
     opts_pv@lbLabelFont                   = "Helvetica"
     opts_pv@lbBoxEndCapStyle              = "TriangleBothEnds"

     if (lay .eq. 0) then 

      if (var0 .eq. "pv") then 
       opts_pv@lbTitleString               = hy0+" m potential vorticity (PVU)"
      elseif (var0 .eq. "rv") then 
       opts_pv@lbTitleString               = hy0+\
       			        	        " m relative vorticity (10~S~-4~N~ s~S~-1~N~)"  
      end if 
    
     elseif (lay .eq. 1) then 
      if (paper .eq. 1) then 
       if (var0 .eq. "pv") then
        opts_pv@lbTitleString               = ""
       elseif (var0 .eq. "rv") then
        opts_pv@lbTitleString               = ""
       end if 
      else
       if (var0 .eq. "pv") then
        opts_pv@lbTitleString               = "Layer-averaged ("+hy0+"-"+hy1+\
                                               " m) ~C~  potential vorticity (PVU)"
       elseif (var0 .eq. "rv") then
        opts_pv@lbTitleString               = "Layer-averaged ("+hy0+"-"+hy1+\
                                              " m) ~C~   relative vorticity (10~S~-4~N~ s~S~-1~N~)"
       end if 
      end if 
     end if 

   ; Additional plotting resources
     opts_pv@mpDataBaseVersion             = "Ncarg4_1" ; More recent database
     opts_pv@mpDataSetName                 = "Earth..4" ; High resolution
     opts_pv@mpOutlineBoundarySets         = "National" ; National borders
     opts_pv@mpGeophysicalLineColor        = "black"    ; Colour borders black
     opts_pv@mpGeophysicalLineThicknessF   = 1.0        ; Border line thickness
     opts_pv@mpGridAndLimbOn               = False      ; Turn on lat/lon lines
     opts_pv@pmTickMarkDisplayMode         = "Always"   ; Turn on map tickmarks
     opts_pv@tmXBMajorLengthF              = 0.005      ; Change tickmark length
     opts_pv@tmXTOn                        = "False"    ; No tickmarks on top x-axis
     opts_pv@tmYROn                        = "False"    ; No tickmarks on right y-axis
     opts_pv@gsnMaximize                   = False      ; Maximise plot size

   ; Make lat/lon labels bigger (01/12/2019)
     opts_pv@tmXBLabelFontHeightF          = 0.02
     opts_pv@tmYLLabelFontHeightF          = 0.02

     opts_pv@gsnAddCyclic                  = False
     opts_pv@gsnDraw                       = False      ; Do not draw the plot
     opts_pv@gsnFrame                      = False      ; Do not advance the frame

   ; When reading in 'pd' stream data with 119 time intervals (PV)
   ; use the index one less than that in 'centre'
     opts_pv@mpLimitMode                   = "Corners"
     opts_pv@mpLeftCornerLatF              = lat_max - r0 ; centre(cn,it-1,0)-r
     opts_pv@mpLeftCornerLonF              = lon_max - r0 ; centre(cn,it-1,1)-r
     opts_pv@mpRightCornerLatF             = lat_max + r0 ; centre(cn,it-1,0)+r
     opts_pv@mpRightCornerLonF             = lon_max + r0 ; centre(cn,it-1,1)+r

   ; Position panel on the page 
     opts_pv@vpXF                          = 0.28
     opts_pv@vpYF                          = 0.95
     opts_pv@vpWidthF                      = 0.40
     opts_pv@vpHeightF                     = 0.40

   ; Vertical velocity (lower)
     opts_vv                               = True
     opts_vv@cnFillOn                      = False
     opts_vv@cnLineColor                   = "black"
     opts_vv@cnInfoLabelOn                 = False
     opts_vv@cnLineLabelsOn                = False
     opts_vv@cnLineLabelInterval           = 2.0
     opts_vv@cnLevelSelectionMode          = "ExplicitLevels"
     opts_vv@cnLevels                      = (/1.5/)
     opts_vv@cnLineLabelPlacementMode      = "constant"
     opts_vv@cnLineLabelPerimOn            = False
     opts_vv@gsnContourLineThicknessesScale = 4.0
     opts_vv@tiMainString                  = ""
     opts_vv@gsnLeftString                 = ""
     opts_vv@gsnRightString                = ""
     opts_vv@gsnDraw                       = False      ; Do not draw the plot
     opts_vv@gsnFrame                      = False      ; Do not advance the frame
     opts_vv@pmTickMarkDisplayMode         = "Always"   ; Turn on map tickmarks
     opts_vv@tmXBMajorLengthF              = 0.005      ; Change tickmark length
     opts_vv@tmXTOn                        = "False"    ; No tickmarks on top x-axis
     opts_vv@tmYROn                        = "False"    ; No tickmarks on right y-axis
     opts_vv@gsnMaximize                   = False      ; Maximise plot size

   ; Vertical velocity (upper)
     opts_vv2				   = opts_vv
     opts_vv2@cnLineColor		   = "purple3"

   ; Change units before plotting 
     vort_plane  = vort_plane * (10 ^ 4)
     vort_plane@units = "10~S~-4~N~ s~S~-1~N~"

     copy_VarCoords(u_plane, vort_plane)

   ; Potential vorticity; relative vorticity
     if (var0 .eq. "pv") then
      plot_pv     = gsn_csm_contour_map(wks,pv_plane,opts_pv)
     elseif (var0 .eq. "rv") then
      plot_pv     = gsn_csm_contour_map(wks,vort_plane,opts_pv)
     end if

   ; Vertical velocity (lower; upper)
     plot_vv      = gsn_csm_contour(wks,w_plane,opts_vv) 
     overlay(plot_pv, plot_vv)
     plot_vv2     = gsn_csm_contour(wks,w_plane2,opts_vv2)
     overlay(plot_pv, plot_vv2)

;=========================================================
; Now plot rate of strain using different colour table
;=========================================================

   ; 850 hPa rate of strain
     opts_str                              = opts_pv
     opts_str@cnFillOn                     = True
     opts_str@cnLineLabelInterval          = 2.0
     opts_str@cnLineLabelFontHeightF       = 0.012
     opts_str@cnLineLabelBackgroundColor   = "transparent"
     opts_str@cnLineLabelPlacementMode     = "constant"

   ; Contour lines off; contour labels off 
     opts_str@cnLinesOn                    = False
     opts_str@cnInfoLabelOn                = False

   ; Colour table; contour levels
     opts_str@cnLevelSelectionMode         = "ExplicitLevels"
     delete(opts_str@cnLevels)
     delete(opts_str@cnFillColors)
     if (diff .eq. 0) then 
      opts_str@cnFillPalette               = "circ_budget"
      opts_str@cnLevels                    = (/10,12,14,16,18,20,22,25,30/)
      opts_str@cnFillColors                = (/-1,11,13,15,16,17,18,19,20,21/)
     elseif (diff .eq. 1) then 
      opts_str@cnFillPalette               = "BlueDarkRed18"
      opts_str@cnLevels                    = (/-3., -2., -1., -0.5, -0.3, -0.2, -0.1, \
                                               0, 0.1, 0.2, 0.3, 0.5, 1., 2., 3./)
      opts_str@cnFillColors                = (/0,2,3,4,5,6,7,-1,-1,10,11,12,13,14,15,17/)
     else
      opts_str@cnFillPalette               = "BlueDarkRed18"
      opts_str@cnLevels                    = (/-10., -5., -3., -2, -1, -0.5, -0.2, \
                                               0, 0.2, 0.5, 1, 2, 3., 5., 10./)
      opts_str@cnFillColors                = (/0,2,3,4,5,6,7,-1,-1,10,11,12,13,14,15,17/)
     end if 
     opts_str@gsnPaperOrientation          = "landscape"
     opts_str@tiMainString                 = ""
     opts_str@tiMainFontHeightF            = 0.0125
     opts_str@gsnLeftString                = ""
     opts_str@gsnRightString               = ""

     if (paper .eq. 1) then 
      opts_str@lbTitleString  = ""
     else
      if (lay .eq. 0) then
       opts_str@lbTitleString = hy0+" m strain rate (10~S~-4~N~ s~S~-2~N~)"
      else
       opts_str@lbTitleString = "Layer-averaged ("+hy0+"-"+hy1+\
                               " m) ~C~         strain rate (10~S~-4~N~ s~S~-2~N~)"
      end if
     end if 

   ; When reading in 'pc' stream data with 121 time intervals (u,v)
   ; use the index one higher than that in 'centre'
     opts_str@mpLeftCornerLatF             = lat_max - r0 ; centre(cn,it-1,0)-r
     opts_str@mpLeftCornerLonF             = lon_max - r0 ; centre(cn,it-1,1)-r
     opts_str@mpRightCornerLatF            = lat_max + r0 ; centre(cn,it-1,0)+r
     opts_str@mpRightCornerLonF            = lon_max + r0 ; centre(cn,it-1,1)+r

   ; Now position strain rate panel on page 
     opts_str@vpXF                         = 0.28
     opts_str@vpYF                         = 0.48
     opts_str@vpWidthF                     = 0.40
     opts_str@vpHeightF                    = 0.40

;================================================================================== 
; Calculate whether flow is dominated by rotation or strain (Rozoff et al. 2006)
;==================================================================================

   ; Calculate square of the relative vorticity 
     vort_plane0 = vort_plane / (10 ^ 4)
     vort_plane1 = vort_plane0 ^ 2
     copy_VarCoords(vort_plane, vort_plane1)
     copy_VarAtts(vort_plane, vort_plane1)
     copy_VarMeta(vort_plane, vort_plane1)

   ; Calculate difference between 'vort_plane1' and 'strT_plane'
     str_diff    = strT_plane - vort_plane1
     copy_VarCoords(strT_plane, str_diff)

   ; Change units before plotting 
     str_plane   = str_plane * (10 ^ 4)
     str_plane@units = "10~S~-4~N~ s~S~-2~N~"

     str_diff    = str_diff * (10 ^ 5)
     str_diff@units  = "10~S~-5~N~ s~S~-2~N~"

   ; Plot strain rate using resources above
     if (diff .eq. 1) then 
      plot_str   = gsn_csm_contour_map(wks,str_diff,opts_str)
     elseif (diff .eq. 0) then 
      plot_str   = gsn_csm_contour_map(wks,str_plane,opts_str)
     else 
      plot_str   = gsn_csm_contour_map(wks,lap_plane,opts_str)
     end if 
     panel(1)    = plot_str

   ; Finally, draw the plot with everything overlaid (call 'draw' and 'frame')
     pres        = True 
     maximize_output(wks,pres)

  ; Tidy up before next iteration
    delete([/pv_plane, w_plane, w_plane2, prs_plane, lap_plane/])
    delete([/str_plane, str0, str1, str2/])
    delete([/strT_plane, str_diff/])

  ; Tidy up (avoid dimension size errors)
    delete([/du_dx,dv_dy,dv_dx,du_dy/])
    delete([/lonC,latC,lon1C,lat1C,u,v,vort,avo,f0,f1/])
    delete([/lonD,latD,lon1D,lat1D,t,p,w,z,pv/])

    if (full .eq. 1) then 
     delete([/th_w,th_e/])
    end if 

  end do      ; End time loop (do it = 0, numINPUT-1)

end