; Script to plot vorticity and strain rate on model height levels, following a TC 

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/times_arr.ncl"

; ncl dat=\"02T12\" opt=\"x11\" mlev=13 dist=3.0 r0=1.0 ens0=\"em11\" sr=1 \
; ts=563 tf=671 ti=6 var0=\"rv\" full=0 tc=0 lay=1 mlev2=25 cn0=\"slp\" tc_ens_pv_str_ml5min.ncl

; 563 - 671 = T+48 --> T+57

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'ens0'  = ensemble simulation (em00-em11)
; 'dist'  = size of box following storm (degrees)
; 'opt'   = output file format ("pdf" or "x11")
; 'mlev'  = model level for plotting PV and vertical velocity (see L106-108)
; 'r0'    = radius for final plot (degrees)
; 'var0'  = plot either PV ("pv") or relative vorticity ("rv")
; 'ts'	  = starting time index (0-116)
; 'tf'	  = ending time index (0-116)
; 'full'  = produce 3-panel (1) or 2-panel (0) plot
; 'tc'    = annotate images with TC symbol (1) or leave as they were (0)
; 'lay'   = calculate layer average (1)
; 'mlev2' = upper model level for calculating layer average (only valid when lay=1)
; 'sr'	  = calculate storm-relative vorticity (1) or leave alone (0) 
; 'cn0'   = storm tracking method ("slp" or "slpf" --> using finescale tracking method)

begin

;==============================
; Prelash (reading in files)
;==============================

; Troubleshooting
  if (r0 .gt. dist) then 
   print("Choose a bigger domain ('dist') or smaller plot size ('r0')!")
   exit()
  end if  

; Find data streams to read in and plot 
  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0+"/vc"                 ; T+54 to T+66
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0

  fili_c1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pc.nc")
  fili_d1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pd.nc")
  fili_j1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pj.nc")	
  fili_k1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pk.nc")

  fili_c      = diri+"/"+fili_c1
  fili_d      = diri+"/"+fili_d1
  fili_j      = diri+"/"+fili_j1
  fili_k      = diri+"/"+fili_k1

; Calculate number of times in file using built-in function ('nc_times') - should be same here
  numINFO_c   = nc_times(fili_c)
  numINFO_d   = nc_times(fili_d)
  numINFO_j   = nc_times(fili_j)
  numINFO_k   = nc_times(fili_k)
  numINPUT_c  = numINFO_c[0]            ; Number of times in 'pc' stream
  time_c      = numINFO_c[1]            ; Array of times from 'pc' stream file
  numINPUT_d  = numINFO_d[0]            ; Number of times in 'pd' stream
  time_d      = numINFO_d[1]            ; Array of times from 'pd' stream file
  numINPUT_j  = numINFO_j[0]            ; Number of times in 'pj' stream
  time_j      = numINFO_j[1]            ; Array of times from 'pj' stream file
  numINPUT_k  = numINFO_k[0]            ; Number of times in 'pk' stream
  time_k      = numINFO_k[1]            ; Array of times from 'pk' stream file

; Calculate storm motion using built-in function ('st_rm')
; 'dat'     = initialisation time (02T12, 03T00, etc)
; 'ens0'    = ensemble simulation (em00-em11)
; 'diri'    = path to input files (see above)
; 'dist'    = size of box following storm (degrees)
; 'mins'    = analyse 1-h (0) or 5-min (1) data

  storm_rel = st_rm(dat, ens0, "$ar/text/", dist, 1)
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
  d0        = dsize(0)                         ; Dim 0 of 'dsize'
  d1        = dsize(1)                         ; Dim 1 of 'dsize'

; Increase system memory
  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 10000000000
  end setvalues

;=======================================                                                
; Choose which storm track to read in                                                   
;=======================================                                                

  if (cn0 .eq. "slp") then
   cn = 0
  elseif (cn0 .eq. "slpf") then
   cn = 1
  end if

;===========================================================================            
; Create date/time string arrays for all times with user-defined function               
;===========================================================================            

 ; Call external function 'times_arr'                                                   
 ; Returns 'time_str' (12 UTC 4 July), 'time_arr' (), and 'title_arr' ()                
   all_times  = times_arr(dat, ens0)

 ; Return strings e.g. '0600 UTC 04 July'                                               
   time_str   = all_times[0]

 ; Return strings e.g. '04Jul_0600Z'                                                    
   time_arr   = all_times[1]

 ; Return title strings e.g. 'Valid at 0600 UTC 04 July'                                
   title_arr  = all_times[2]

;================================================================
; Create arrays for relative vorticity and strain at each time 
;================================================================
  
 ; Number of time intervals to loop over 
   nt_all   = (tf - ts) + 1

 ; Output on one model level in either case, either single or layer-averaged
   vort_arr = new( (/nt_all, d0, d1/), "float")
   str_arr  = new( (/nt_all, d0, d1/), "float")

;==================================
; Get the variables we will need 
;==================================

  c = addfile(fili_c,"r")          ; Read in 'pc' stream [u,v,vort]
  d = addfile(fili_d,"r")          ; Read in 'pd' stream [w,p,]
  j = addfile(fili_j,"r")	   ; Read in 'pj' stream [t,theta,mixing ratio]
  k = addfile(fili_k,"r")	   ; Read in 'pk' stream [tau,PV,w]

  do it = ts, tf, ti		   ; Loop over times in file (T+54 to T+66, every 5 min)

   time  = c->t(it)                ; Time in file
   times = dimsizes(time)

   print_clock("Working on "+time_str(it))

 ; Get domain boundaries at each time from storm track (calculated off-line)
   t0    = lat_arr(cn,it,0)
   t1    = lat_arr(cn,it,1)
   n0    = lon_arr(cn,it,0)
   n1    = lon_arr(cn,it,1)
   print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

 ; Create latitude/longitude array subsets to reduce computation time
   lonC  = c->longitude({n0:n1})    ; longitude [160 grid points]
   latC  = c->latitude({t0:t1})     ; latitude  [160 grid points]
   hybC  = c->hybrid_ht(:)          ; 63 model ('theta') levels

   u_varname = "x-wind"
   v_varname = "y-wind"

 ; Read in horizontal wind components
 ; Both variables below -- [hybC | 63] * [latC | 150] * [lonC | 150]
   if (lay .eq. 1) then
    u0  = c->$u_varname$(it,mlev:mlev2,{t0:t1},{n0:n1})     ; Zonal wind (layer)
    v0  = c->$v_varname$(it,mlev:mlev2,{t0:t1},{n0:n1})     ; Meridional wind (layer)
   else
    u0  = c->$u_varname$(it,mlev,{t0:t1},{n0:n1})           ; Zonal wind (single level)
    v0  = c->$v_varname$(it,mlev,{t0:t1},{n0:n1})           ; Meridional wind (single level)
   end if

;===================================
; Calculate storm-relative winds
;===================================

    u   = u0 - u_cyc(cn,it)
    v   = v0 - v_cyc(cn,it)

  ; Add metadata from original horizontal wind arrays
    copy_VarCoords(u0,u)
    copy_VarAtts(u0,u)

    copy_VarCoords(v0,v)
    copy_VarAtts(v0,v)

;==================================================================
; Calculate relative vorticity using centered finite differences
;==================================================================

  ; Option '2' --> boundary points estimated using one-sided difference scheme

    if (sr .eq. 1) then 
     vort0 = uv2vr_cfd(u, v, latC, lonC, 2)		; Storm-relative vorticity
     copy_VarCoords(u, vort0)
    elseif (sr .eq. 0) then 
     vort0 = uv2vr_cfd(u0, v0, latC, lonC, 2)		; Earth-relative vorticity
     copy_VarCoords(u0, vort0)
    end if 

    vort0@units = "s~S~-1~N~"
    vort0@name  = "Relative vorticity on model levels"
    
    if (lay .eq. 1) then 
     vort       = vort0(:,0:dsize(0)-1,0:dsize(1)-1)
    else 
     vort       = vort0(0:dsize(0)-1,0:dsize(1)-1)
    end if 

  ; Coriolis parameter and absolute vorticity
    f1   = coriolis_param(latC(0:dsize(0)-1))
    
    if (lay .eq. 0) then 
     f0   = conform_dims(dimsizes(vort),f1,1)
    elseif (lay .eq. 1) then 
     f0   = conform_dims(dimsizes(vort),f1,2)
    end if 

    copy_VarMeta(vort,f0)
    copy_VarCoords(vort,f0)
    f0@description  = "Coriolis parameter"
    f0@name         = "Coriolis parameter"

    avo  	    = vort + f0
    copy_VarMeta(vort,avo)
    copy_VarAtts(vort,avo)
    copy_VarCoords(vort,avo)
    avo@description = "Absolute vorticity"
    avo@name        = "Absolute vorticity on model levels"
    avo@long_name   = "Absolute vorticity"

;=====================================================
; Read in variables from 'pd' stream (model levels)
;=====================================================

  ; Read in basic variables from 'pd' stream
    lonD   = d->longitude({n0:n1})                  ; longitude [xxx grid points]
    latD   = d->latitude({t0:t1})                   ; latitude  [xxx grid points]
    hybD   = d->hybrid_ht(:)                        ; 63 model ('rho') levels

    hy     = hybD / 1000 ; Model level values in km
    hy0    = sprintf("%0.0f",hybD(mlev))
    hy1    = sprintf("%0.0f",hybD(mlev2))
    print("Working on model level: "+hy0+ " m")

  ; Read in temperature, pressure, vertical velocity and geopotential height
  ; All variables below -- [hybD | 63] * [latD | 150] * [lon1D | 151]
    p      = d->p(it,mlev,{t0:t1},{n0:n1})           ; Pressure (Pa)

  ; Convert units of pressure before continuing
    p      = p / 100          ; Convert to 'hPa'

  ; Add metadata
    p@units  = "hPa"

    dx     = lonD(1) - lonD(0) ; Grid spacing (longitude)
    dy     = latD(1) - latD(0) ; Grid spacing (latitude)

;======================================================
; Read in variables from 'pk' stream (model levels)
;======================================================

  ; Read in basic variables from 'pd' stream
    lon1K  = k->longitude_1({n0:n1})                ; longitude [xxx grid points]
    lat1K  = k->latitude_1({t0:t1})                 ; latitude  [xxx grid points]
    lonK   = k->longitude({n0:n1})                  ; longitude [xxx grid points]
    latK   = k->latitude({t0:t1})                   ; latitude  [xxx grid points]
    hybK   = k->hybrid_ht(:)                        ; 63 model ('rho') levels
    hyb1K  = k->hybrid_ht_1(:)			    ; 80 model ('rho') levels

    hyK    = hybK / 1000 ; Model level values in km
    hyk0   = sprintf("%0.0f",hybK(mlev))
    hyk1   = sprintf("%0.0f",hybK(mlev2))

    if (lay .eq. 1) then
     print("Averaging PV between "+hy0+" and "+hy1+ " m")
     pv    = k->field83(it,mlev:mlev2,{t0:t1},{n0:n1})   ; Potential vorticity (PVU)
     w     = k->dz_dt(it,mlev:mlev2,{t0:t1},{n0:n1})     ; Vertical velocity (m s-1)
    else
     pv    = k->field83(it,mlev,{t0:t1},{n0:n1})         ; Potential vorticity (PVU)
     w     = k->dz_dt(it,mlev,{t0:t1},{n0:n1})           ; Vertical velocity (m s-1)
    end if

  ; Convert units of pressure and geopotential height before continuing
    pv     = pv * (10 ^ 6)    ; Convert to 'PVU' [10 ^ 6 K m2 s-1 kg-1]

  ; Add metadata
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
     str1    = (str1a) ^ 2
     str2    = (str2a) ^ 2

     if (lay .eq. 1) then
      str_plane0  = sqrt( str1(:,:,:) + str2(:,:,:) )
     else
      str_plane0  = sqrt( str1(:,:) + str2(:,:) )
     end if
     delete([/str1a,str2a/])

;===================================================================
; Also read in lower-tropospheric PV [proxy for convective cells]
;===================================================================

   ; Read in PV, geopotential height and vertical velocity

     if (lay .eq. 1) then
      pv_plane   = dim_avg_n_Wrap(pv(:,:,:),0)
      vort_plane = dim_avg_n_Wrap(vort(:,:,:),0)
      w_plane    = dim_avg_n_Wrap(w(:,:,:),0)
      str_plane  = dim_avg_n_Wrap(str_plane0(:,:,:),0)
      copy_VarCoords(u(0,:,:),str_plane)
     else
      pv_plane   = pv
      vort_plane = vort
      w_plane    = w
      str_plane  = str_plane0
      copy_VarCoords(u,str_plane)
     end if

   ; Read single-level values into 3D array before plotting 

   vort_arr = new( (/nt_all, d0, d1/), "float")
   str_arr  = new( (/nt_all, d0, d1/), "float")

   ; Choose size of panel plot
     if (full .eq. 1) then
      panel       = new(3,graphic)
     elseif (full .eq. 0) then
      panel       = new(2,graphic)
     end if

;=================================
; Output file type and location
;=================================

     if (lay .eq. 1) then
      output = "$sam/nepartak/images/vc/"+ens0+"_"+dat+"_5min_str_"+var0+\
	       "_"+hy0+"m_"+hy1+"m_"+time_arr(it)
     else
      output = "$sam/nepartak/images/vc/"+ens0+"_"+dat+"_5min_str_"+var0+\
      	       "_"+hy0+"m_"+time_arr(it)
     end if
     wks = gsn_open_wks(opt,output)

   ; Load colour table
     if (var0 .eq. "pv") then
      gsn_define_colormap(wks,"prcp_new")		; Potential vorticity
     elseif (var0 .eq. "rv") then
      gsn_define_colormap(wks,"prcp_new")		; Relative vorticity
     end if

;==========================
; Options for plotting
;==========================

   ; Potential vorticity
     opts_pv                               = True
     opts_pv@cnFillOn                      = True
     opts_pv@cnLineLabelInterval           = 2.0
     opts_pv@cnLineLabelFontHeightF        = 0.012
     opts_pv@cnLineLabelBackgroundColor    = "transparent"
     opts_pv@cnLineLabelPlacementMode      = "constant"
     opts_pv@cnLinesOn                     = False ; Contour lines off
     opts_pv@cnInfoLabelOn                 = False ; Contour labels off
     opts_pv@cnLevelSelectionMode          = "ExplicitLevels"

     if (var0 .eq. "pv") then 
      opts_pv@cnLevels                     = (/9.0, 10.0 ,12.0, 14.0, \
                                               16.0, 18.0, 20.0, 22.0, 25.0, \
                                               30.0, 35.0, 40.0, 45.0, 50.0/)
      opts_pv@cnFillColors                 = (/0,2,3,4,5,\
                                               6,7,8,9,10,\
                                               11,12,13,14,15/) 
     elseif (var0 .eq. "rv") then 
      opts_pv@cnLevels                     = (/5., 10., 20., 30., 40., 50., \
                                               60., 80., 100., 120., 150./)
      opts_pv@cnFillColors                 = (/0,4,5,6,7,8,9,10,11,12,13,14/)
     end if 

     opts_pv@gsnPaperOrientation           = "landscape"
     opts_pv@tiMainString                  = ""
     opts_pv@tiMainFontHeightF             = 0.0125
     opts_pv@gsnLeftString                 = ""
     opts_pv@gsnRightString                = ""
     opts_pv@lbLabelBarOn                  = True       ; Individual lb on/off

     opts_pv@pmLabelBarWidthF              = 0.52
     opts_pv@pmLabelBarHeightF             = 0.13
     opts_pv@lbLabelFontHeightF            = 0.0125
     opts_pv@lbPerimOn                     = False
     opts_pv@lbLabelFont                   = "Helvetica"

     if (lay .eq. 0) then

      if (var0 .eq. "pv") then
       opts_pv@lbTitleString                 = hy0+" m potential vorticity (PVU)"
      elseif (var0 .eq. "rv") then
       opts_pv@lbTitleString                 = hy0+\
						" m relative vorticity (10~S~-4~N~ s~S~-1~N~)"
      end if

     elseif (lay .eq. 1) then

      if (var0 .eq. "pv") then
       opts_pv@lbTitleString             = "Layer-averaged ("+\
                                           hy0+"-"+hy1+" m) ~C~     potential vorticity (PVU)"
      elseif (var0 .eq. "rv") then
       opts_pv@lbTitleString             = "Layer-averaged ("+hy0+"-"+hy1+\
                                           " m) ~C~     relative vorticity (10~S~-4~N~ s~S~-1~N~)"

      end if

     end if

     opts_pv@lbBoxEndCapStyle              = "TriangleBothEnds"    ; Labelbar end shape
     opts_pv@lbTitleFontHeightF            = 0.02
     opts_pv@lbTitleDirection              = "Across"
     opts_pv@lbTitlePosition               = "Bottom"

   ; Additional plotting resources (database; resolution; borders; colour; thickness; lat/lon)
     opts_pv@mpDataBaseVersion             = "Ncarg4_1"
     opts_pv@mpDataSetName                 = "Earth..4"
     opts_pv@mpOutlineBoundarySets         = "National"
     opts_pv@mpGeophysicalLineColor        = "black"   
     opts_pv@mpGeophysicalLineThicknessF   = 1.0       
     opts_pv@mpGridAndLimbOn               = False     

   ; Tickmarks; length; maximise?
     opts_pv@pmTickMarkDisplayMode         = "Always"
     opts_pv@tmXBMajorLengthF              = 0.005   
     opts_pv@tmXTOn                        = "False" 
     opts_pv@tmYROn                        = "False" 
     opts_pv@gsnMaximize                   = True    

     r                                     = r0         ; Radius of plot (user input)

     opts_pv@gsnAddCyclic                  = False
     opts_pv@mpLimitMode                   = "Corners"
     opts_pv@mpLeftCornerLatF              = centre(cn,it,0)-r
     opts_pv@mpLeftCornerLonF              = centre(cn,it,1)-r
     opts_pv@mpRightCornerLatF             = centre(cn,it,0)+r
     opts_pv@mpRightCornerLonF             = centre(cn,it,1)+r
     opts_pv@gsnDraw                       = False      ; Do not draw the plot
     opts_pv@gsnFrame                      = False      ; Do not advance the frame

   ; Vertical velocity
     opts_vv                               = True
     opts_vv@cnFillOn                      = False
     opts_vv@cnLineColor                   = "black"
     opts_vv@cnInfoLabelOn                 = False
     opts_vv@cnLineLabelsOn                = True
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
     opts_vv@gsnMaximize                   = True       ; Maximise plot size

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
     opts_str@cnLinesOn                    = False ; Contour lines off
     opts_str@cnInfoLabelOn                = False ; Contour labels off
     opts_str@cnFillPalette                = "circ_budget"
     opts_str@cnLevelSelectionMode         = "ExplicitLevels"
     delete(opts_str@cnLevels)
     delete(opts_str@cnFillColors)
     opts_str@cnLevels                     = (/10,12,14,16,18,20,22,25,30/)
     opts_str@cnFillColors                 = (/-1,11,13,15,16,17,18,19,20,21/)
     opts_str@gsnPaperOrientation          = "landscape"
     opts_str@tiMainString                 = ""
     opts_str@tiMainFontHeightF            = 0.0125
     opts_str@gsnLeftString                = ""
     opts_str@gsnRightString               = ""
     opts_str@lbLabelBarOn                 = True             ; Individual labelbar

     opts_str@pmLabelBarWidthF             = 0.52
     opts_str@pmLabelBarHeightF            = 0.13
     opts_str@lbLabelFontHeightF           = 0.0125
     opts_str@lbPerimOn                    = False
     opts_str@lbLabelFont                  = "Helvetica"

     if (lay .eq. 0) then
      opts_str@lbTitleString = hy0+" m strain rate (10~S~-4~N~ s~S~-2~N~)"
     elseif (lay .eq. 1) then
      opts_str@lbTitleString = "Layer-averaged ("+hy0+"-"+hy1+\
                               " m) ~C~         strain rate (10~S~-4~N~ s~S~-2~N~)"
     end if

     opts_str@lbBoxEndCapStyle             = "TriangleBothEnds"    ; Labelbar end shape
     opts_str@lbTitleFontHeightF           = 0.02
     opts_str@lbTitleDirection             = "Across"
     opts_str@lbTitlePosition              = "Bottom"

;================================================
; Finally plot wet bulb potential temperature 
;================================================

   if (full .eq. 1) then 

   ; Read in "radar_new" colour map
     cmap = read_colormap_file("radar_new")
   
   ; Make the first 2 colours (black and white) transparent
     cmap(0:1,3) = 0.0

   ; Wet bulb potential temperature
     opts_the                              = opts_pv
     opts_the@cnFillOn                     = True
     opts_the@cnLineLabelInterval          = 2.0
     opts_the@cnLineLabelFontHeightF       = 0.012
     opts_the@cnLineLabelBackgroundColor   = "transparent"
     opts_the@cnLineLabelPlacementMode     = "constant"
     opts_the@cnLinesOn                    = False ; Contour lines off
     opts_the@cnInfoLabelOn                = False ; Contour labels off
     opts_the@cnFillPalette                = "radar_new"
     opts_the@cnLevelSelectionMode         = "ExplicitLevels"
     delete(opts_the@cnLevels)
     delete(opts_the@cnFillColors)
     opts_the@cnLevels			   = ispan(tavg0-5,tavg0+6,1)
     opts_the@cnFillColors		   = (/2,3,4,5,6,7,8,9,10,11,12,13,14/)
     opts_the@gsnPaperOrientation          = "landscape"
     opts_the@tiMainString                 = ""
     opts_the@tiMainFontHeightF            = 0.0125
     opts_the@gsnLeftString                = ""
     opts_the@gsnRightString               = ""
     opts_the@lbLabelBarOn                 = True             ; Individual labelbar

     opts_the@pmLabelBarWidthF             = 0.52
     opts_the@pmLabelBarHeightF            = 0.13
     opts_the@lbLabelFontHeightF           = 0.0125
     opts_the@lbPerimOn                    = False
     opts_the@lbLabelFont                  = "Helvetica"

     opts_the@lbTitleString = pres+" hPa theta-e (K)"

     opts_the@lbBoxEndCapStyle             = "TriangleBothEnds"    ; Labelbar end shape
     opts_the@lbTitleFontHeightF           = 0.02
     opts_the@lbTitleDirection             = "Across"
     opts_the@lbTitlePosition              = "Bottom"

   end if 

;=============
; Plot data
;=============

     str_plane   = str_plane * (10 ^ 4)
     str_plane@units = "10~S~-4~N~ s~S~-2~N~"

     vort_plane  = vort_plane * (10 ^ 4)
     vort_plane@units = "10~S~-4~N~ s~S~-1~N~"

     if (var0 .eq. "pv") then 
      plot_pv    = gsn_csm_contour_map(wks,pv_plane,opts_pv)   ; Potential vorticity
     elseif (var0 .eq. "rv") then 
      plot_pv    = gsn_csm_contour_map(wks,vort_plane,opts_pv) ; Relative vorticity     
     end if 

     plot_vv     = gsn_csm_contour(wks,w_plane,opts_vv)       ; Vertical velocity
     overlay(plot_pv, plot_vv)
     panel(0)   = plot_pv

     plot_str   = gsn_csm_contour_map(wks,str_plane,opts_str) ; Strain rate
     panel(1)   = plot_str

;     if (full .eq. 1) then 
;
;      plot_the  = gsn_csm_contour_map(wks,the_plane,opts_the) ; Equivalent potential temp
;      panel(2)	= plot_the
;
;     end if 

;==================================================
; Add NH tropical cyclone symbol to panel plots
;==================================================

  if (tc .eq. 1) then 

    markers    = new((/dimsizes(panel)/),graphic)

  ; Add NH tropical cyclone symbol to index
    tc0 = NhlNewMarker(wks, "p", 137, 0.0, 0.0, 1.0, 1.0, 0.0)

    mres0               = True
    mres0@gsMarkerIndex = tc0         ; TC symbol
    mres0@gsMarkerSizeF = 20.0        ; Marker size
    mres0@gsMarkerColor = "grey55"    ; Marker colour
    mres0@gsMarkerThicknessF = 5.0    ; Marker thickness (increase)

    do m = 0, dimsizes(panel)-1
     markers(m) = gsn_add_polymarker(wks,panel(m),centre(it,1),centre(it,0),mres0)
    end do      

  end if 

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    optsP                       = True
    optsP@gsnFrame              = False    ; Do not advance the frame
    optsP@gsnPanelLabelBar      = False    ; Turn off panel labelbar
    optsP@txString              = "Simulation "+ens0+" : "+title_arr(it)
    optsP@gsnPanelFigureStrings = ""
    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.90
    optsP@gsnPanelBottom        = 0.10
    optsP@amJust                = "TopLeft"
    optsP@gsnPanelFigureStringsFontHeightF = 0.0125 ; Label size (default 0.01)

    if (full .eq. 1) then
     gsn_panel(wks,panel,(/1,3/),optsP)
     frame(wks)
    elseif (full .eq. 0) then
     gsn_panel(wks,panel,(/1,2/),optsP)
     frame(wks)
    end if

  ; Tidy up before next iteration
    delete([/vort_plane, pv_plane, w_plane/])
    delete([/str_plane, str_plane0, str1, str2/])

  ; Tidy up (avoid dimension size errors)
    delete([/du_dx,dv_dy,dv_dx,du_dy/])
    delete([/lonC,latC,u,v,u0,v0,vort,vort0,avo,f0,f1/])
    delete([/lonD,latD,p/])
    delete([/lon1K,lat1K,lonK,latK,hybK,pv,w/])

  end do      ; End time loop (do it = 0, numINPUT-1)

end 