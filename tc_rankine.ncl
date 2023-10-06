; Compare storm-tracking methods and overlay on relative vorticity or windspeed 

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/setup_cyl.ncl"
load "$sam/ncl_func/st_centre.ncl"

; ncl dat=\"02T12\" opt=\"x11\" ens0=\"em11\" rad0=2.0 dist=3.0 vc=0 sm=0 zoom=1 plt=1
; lgd0=\"tr\" cn0=\"geo_sm\" ts0=51 tf0=52 sc=4 lay=1 cyl=1 dnt=0 prs=0 
; ar=0.75 nr=16 var0=\"vort\" tc_rankine.ncl

; 'dat'  = initialisation time string ("02T12", "03T00", ...)
; 'ens0' = ensemble simulation (em00-em11)
; 'int'	 = 
; 'opt'  = output file format ("pdf" or "x11")
; 'rad0' = 
; 'dist' = size of box following storm (degrees)
; 'vc'   = 
; 'lay'  = calculate layer average (1) or single level (0)
; 'p_cp' = pressure level for calculating storm centre
; 'lgd0' = legend position ("tr","tl","br","bl")
; 'sc'   = method of choosing storm centre position

begin

;=====================================
; Read in 'pc' and 'pd' stream data
;=====================================

  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0

; Print out relevant information to user about storm tracking 
  sc_arr      = (/"relative vorticity",\                          ; sc = 0
                  "pressure",\                                    ; sc = 1
                  "smoothed pressure",\                           ; sc = 2
                  "extra smoothed pressure",\                     ; sc = 3 
                  "vort_min during symmetric phase",\             ; sc = 4 
                  "windspeed minimum",\                           ; sc = 5
                  "geopotential"/)                                ; sc = 6

  print("Using "+sc_arr(sc)+" to calculate storm centre position...")  

; Horizontal windspeed on model levels 
  fili_c1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pc.nc")
  fili_c      = diri+"/"+fili_c1

; Pressure on model levels 
  fili_d1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pd.nc")
  fili_d      = diri+"/"+fili_d1 

; Get number of times in file using built-in function ('nc_times')
; [0] = number of time intervals
; [1] = information corresponding to time intervals
  numINFO_c   = nc_times(fili_c)
  numINPUT_c  = numINFO_c[0]
  time_c      = numINFO_c[1]       ; 121 times 

  numINFO_d   = nc_times(fili_d)
  numINPUT_d  = numINFO_d[0]
  time_d      = numINFO_d[1]       ; 119 times 
  num_l       = 1

; Calculate storm motion using built-in function ('st_rm')
; 'dat'     = initialisation time (02T12, 03T00, etc) 
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
  llbox     = toint(dist*50)                   ; Calculate domain size (grid points)
  dsize     = (/llbox,llbox/)                  ; Domain size (grid points) 

; Array to hold values of symmetry parameter to write out 
  sym_arr   = new(numINPUT_d,float)

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

  utc_date = cd_calendar(time_d, 0)

  year   = tointeger(utc_date(:,0))
  month  = tointeger(utc_date(:,1))
  day    = tointeger(utc_date(:,2))
  hour   = tointeger(utc_date(:,3))
  minute = tointeger(utc_date(:,4))
  second = utc_date(:,5)

; Correct for errors in the code (round up value of hour when minutes = 59) 
  do ct = 0, numINPUT_d-1
   if (minute(ct).gt.30) then
     hour(ct) = hour(ct)+1
   end if
  end do
  delete(ct)

; Create arrays to hold strings for output plots
  date_str  = new(numINPUT_d,string)
  time_str  = new(numINPUT_d,string)
  time_arr  = new(numINPUT_d,string)
  title_arr = new(numINPUT_d,string)

  date_str  = sprinti("%0.2iUTC ", hour) + \
              sprinti("%0.2i ", day) \
            + month_abbr(month)
  time_str  = sprinti("%0.2i UTC ", hour) + \
              sprinti("%0.2i ", day) \
            + month_abbr(month)
  time_arr  = sprinti("%0.2i", day) + \
              month_abbr(month) + \
              "_" + sprinti("%0.2iZ", hour)

  print("Working on ensemble member "+ens0+" ("+dat+")")

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 10000000000
  end setvalues

;==========================================================
; Read in required variables from 'pc' and 'pd' streams 
;==========================================================

; Create array to hold storm track positions using different methods 
  centre_new  = new((/7,num_l,numINPUT_d,2/),"float")

; We are using 'time_d' as our foundation, which starts at T+2 
; So if ts0=48 and tf0=72, we're actually choosing T+50 and T+74
  c  = addfile(fili_c,"r")              ; Read in 'pc' stream [u,v] - 121 times
  d  = addfile(fili_d,"r")              ; Read in 'pd' stream [w,p,t,theta,kmh,tau] - 119 times 

  do it	       = ts0, tf0               ; Loop over times in file 

   print("Working on time: "+time_str(it-2)+" (T+"+it+")" )
   title_arr(it-2) = "Valid at "+time_str(it-2)+" (T+"+it+")"

 ; Get domain boundaries at each time from storm track (calculated off-line)
   t0          = lat_arr(cn,it-1,0)
   t1          = lat_arr(cn,it-1,1)
   n0          = lon_arr(cn,it-1,0)
   n1          = lon_arr(cn,it-1,1)

 ; Combine values above into array for input into external function 
   ll_arr = (/t0,t1,n0,n1/)

 ; Read in basic variables from 'pc' stream
   lon1C  = c->longitude_1({n0:n1})              ; longitude ['d1' grid points]
   lat1C  = c->latitude_1({t0:t1})               ; latitude  ['d0' grid points]
   lonC   = c->longitude({n0:n1})                ; longitude ['d1' grid points]
   latC   = c->latitude({t0:t1})                 ; latitude  ['d0' grid points]
   hybC   = c->hybrid_ht(:)                      ; 63 model ('theta') levels

   dy     = latC(1) - latC(0)
   dx     = lonC(1) - lonC(0)

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

 ; Read in pressure and vertical velocity 
   if (lay .eq. 1) then  
    p     = d->p(it-2,13:25,{t0:t1},{n0:n1})      ; Pressure (Pa)
    w     = d->dz_dt(it-2,13:25,{t0:t1},{n0:n1})  ; Vertical velocity (m s-1)
    z     = d->ht(it-2,13:25,{t0:t1},{n0:n1})     ; Geopotential height (m)
    th    = d->theta(it-2,13:25,{t0:t1},{n0:n1})  ; Potential temperature (K)
   else
    p     = d->p(it-2,ilev,{t0:t1},{n0:n1})       ; Pressure (Pa)
    w     = d->dz_dt(it-2,ilev,{t0:t1},{n0:n1})   ; Vertical velocity (m s-1) 
    z     = d->ht(it-2,ilev,{t0:t1},{n0:n1})      ; Geopotential height (m) 
    th    = d->theta(it-2,ilev,{t0:t1},{n0:n1})   ; Potential temperature (K)
   end if 

 ; Convert units of pressure and geopotential height before continuing 
   p      = p / 100
   z      = z / 10

 ; Add metadata
   p@units = "hPa"
   z@units = "dam"

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

 ; Combine all pc and pd stream variables into larger arrays
 ; pc_plane needs to contain 6 arrays; the final array is just a placeholder here
   pc_plane      = (/u_plane, v_plane, spd_plane, \
                     vort_plane, avo_plane, avo_plane/)
   pd_plane      = (/w_plane, prs_plane/)

;=============================================================
; Call external function to calculate storm centre position
;=============================================================

 ; Set radial distance (º) when looking for vort min in ringlike phase 
   r0            = 0.10

   centre_arr    = st_centre(0, it, cn, r0, centre, vort_plane, prs_plane, spd_plane, \
                             centre_new, dy, dx, "single", it, "pd")

 ; Updated storm track information contained in 'centre_new'
 ; centre_new(0,0,it-2,:) --> vort_max 
 ; centre_new(1,0,it-2,:) --> prs_min
 ; centre_new(2,0,it-2,:) --> prs_min_sm
 ; centre_new(3,0,it-2,:) --> prs_min_sm2
 ; centre_new(4,0,it-2,:) --> vort_min
 ; centre_new(5,0,it-2,:) --> spd_min
 ; centre_new(6,0,it-2,:) --> offline
   centre_new     = centre_arr[0]

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
   cyl_arr  = setup_cyl(it, dist, ll_arr, lonC, latC, lonD, latD, lon1D, lat1D, \
                        pc_plane, pd_plane, numTIMES, centre_new, \
                        sc, sm, ar, nr, ilev, num_l, 0, "single", it)

   lat_max  = cyl_arr[4]
   lon_max  = cyl_arr[5]
   radii    = cyl_arr[6]
   rad_size = cyl_arr[7]
   thetas   = cyl_arr[19]

   u_int    = cyl_arr[8]
   v_int    = cyl_arr[9]
   spd_int  = cyl_arr[10]
   vort_int = cyl_arr[11]
   prs_int  = cyl_arr[16]

   xpos     = cyl_arr[17]
   ypos     = cyl_arr[18]

   vort_azi = new( (/dimsizes(radii)/), float)
   spd_azi  = new( (/dimsizes(radii)/), float)

 ; Now calculate the azimuthal average of relative vorticity 
   do irad = 0, dimsizes(radii)-1
    vort_azi(irad) = avg( vort_int(irad,:) )
    spd_azi(irad)  = avg( spd_int(irad,:) )
   end do 

 ; Use more sophisticated method to calculate whether vortex is ringlike or monopolar
 ; This code will go into an NCL function 

 ; (1) Find radius (index) corresponding to relative vorticity maximum
   ind_vmax    = ind_resolve(maxind(vort_azi), dimsizes(vort_azi) )

 ; Calculate vorticity at ind_vmax ('vort_max'), and at the vortex centre ('vort_centre')
   vort_max    = vort_azi(ind_vmax(0,0) )
   vort_centre = vort_azi(0)

 ; (2) Calculate the ratio of 'ind_max' to 'ind_centre' 
   sym_arr(it) = vort_centre / vort_max

 ; If 'ind_vmax' = 'ind_centre' (i.e. monopole), then 'sym_arr' will be equal to 1  
   if (sym_arr(it) .gt. 0.95) then 

    print("Inner-core structure is monopolar (ratio = "+sym_arr(it)+")")
    sym_str = "Monopole: "+sprintf("%0.2f", sym_arr(it) )

   elseif (sym_arr(it) .lt. 0.85) then 

    print("Inner-core structure is ringlike (ratio = "+sym_arr(it)+")")
    sym_str = "Ringlike: "+sprintf("%0.2f", sym_arr(it) )

   else 
    print("Inner-core structure is not clearly monopolar or ringlike (ratio = "+sym_arr(it)+")")
    sym_str = "Not clear: "+sprintf("%0.2f", sym_arr(it) )

   end if 

 ; If we're in the symmetric phase, put the storm centre in the middle of the vorticity ring 
   if (dnt .eq. 1 .and. sym_arr(it) .lt. 0.85) then 

    print("Re-calculating storm centre to fit ringlike vorticity profile...")

  ; Select the 'vort_min' method of storm tracking 
    sc = 4

  ; Call the external function again
    delete(cyl_arr)
    cyl_arr  = setup_cyl(it, dist, ll_arr, lonC, latC, lonD, latD, lon1D, lat1D, \
                         pc_plane, pd_plane, numTIMES, centre_new, \
                         sc, sm, ar, nr, ilev, num_l, 0, "single", it)

    lat_max  = cyl_arr[4]
    lon_max  = cyl_arr[5]
    radii    = cyl_arr[6]
    rad_size = cyl_arr[7]
    vort_int = cyl_arr[11]
    prs_int  = cyl_arr[16]
    xpos     = cyl_arr[17]
    ypos     = cyl_arr[18]

  ; Recalculate azimuthal average of relative vorticity
    do irad = 0, dimsizes(radii)-1
     vort_azi(irad) = avg( vort_int(irad,:) )
     spd_azi(irad)  = avg( spd_int(irad,:) )
    end do

   end if 

 ; Reset value of 'sc'
   sc = sc 

 ; EDIT -- customise further
 ; Plot raw data alongside the azimuthal mean ?

 ; Output file location and type
   if (dnt .eq. 0) then 
    output = "$sam/nepartak/images/basic/vc_vort_rad_"+sc+"_"+dat+"_"+ens0+"_"+time_arr(it-2)
   else
    output = "$sam/nepartak/images/basic/vc_vort_rad_"+sc+"a_"+dat+"_"+ens0+"_"+time_arr(it-2)
   end if 

   wks    = gsn_open_wks(opt,output)

   print(vort_azi)
   exit()

 ; EDIT FROM HERE -- output plots of radius vs relative vorticity 

; ; Load colour table 
;   gsn_define_colormap(wks,"prcp_new")

;==========================
; Options for plotting
;==========================

; Relative vorticity
  opts_rvort                              = True
  opts_rvort@cnFillOn                     = True

  if (cyl .eq. 1 .or. cyl .eq. 2) then 
   opts_rvort@sfXArray                    = xpos(:,:)
   opts_rvort@sfYArray                    = ypos(:,:)
  end if 

  opts_rvort@cnLineLabelInterval          = 2.0
  opts_rvort@cnLineLabelFontHeightF       = 0.012
  opts_rvort@cnLineLabelBackgroundColor   = "transparent"
  opts_rvort@cnLineLabelPlacementMode     = "constant"
  opts_rvort@cnLinesOn                    = False ; Contour lines off
  opts_rvort@cnInfoLabelOn                = False ; Contour labels off
  opts_rvort@cnLevelSelectionMode         = "ExplicitLevels"

  if (var0 .eq. "vort") then 
   opts_rvort@cnLevels                    = (/5., 10., 15., 20., 25., 30., \
                                              40., 50., 60., 70., 80./)
   opts_rvort@cnFillColors                = (/0,4,5,6,7,8,9,10,11,12,13,14/)
  elseif (var0 .eq. "th") then
   mean0                                  = avg(th_int(:,:) )
   mean                                   = toint(mean0)
   levs                                   = ispan(mean-8, mean+8, 1)
   opts_rvort@cnFillPalette               = "BlAqGrYeOrReVi200"
   opts_rvort@cnLevels                    = levs
   opts_rvort@cnFillColors                = ispan(18,188,10)
  end if 

  opts_rvort@gsnPaperOrientation          = "landscape"
  opts_rvort@tiMainString                 = ""
  opts_rvort@tiMainFontHeightF            = 0.0125
  opts_rvort@gsnLeftString                = ""
  opts_rvort@gsnRightString               = ""

  opts_rvort@lbLabelBarOn                 = True
  opts_rvort@lbBoxEndCapStyle             = "TriangleBothEnds"    ; Labelbar end shape
  opts_rvort@lbLabelFontHeightF           = 0.0125                ; Labelbar font size
  opts_rvort@lbLabelFont                  = "Helvetica"           ; Labelbar font

  if (var0 .eq. "vort") then 
   if (lay .eq. 1) then 
    opts_rvort@lbTitleString              = "Layer-averaged ("+hy0+"-"+hy1+\
                                            " m) ~C~     relative vorticity (10~S~-4~N~ s~S~-1~N~)"
   else
    opts_rvort@lbTitleString              = hy2+" m relative vorticity (10~S~-4~N~ s~S~-1~N~)"
   end if 
  elseif (var0 .eq. "th") then 
   if (lay .eq. 1) then
    opts_rvort@lbTitleString              = "Layer-averaged ("+hy0+"-"+hy1+\
                                            " m) ~C~     potential temperature (K)"
   else
    opts_rvort@lbTitleString              = hy2+" m potential temperature (K)"
   end if
  end if 

  opts_rvort@lbTitlePosition              = "Bottom"              ; Title position
  opts_rvort@lbTitleFontHeightF           = 0.015                 ; Title font height
  opts_rvort@pmLabelBarOrthogonalPosF     = 0.1                   ; Labelbar position
  opts_rvort@lbTitleOffsetF               = 0                     ; Lb title up/down
  opts_rvort@lbPerimOn                    = False                 ; Perimeter on/off
  opts_rvort@gsnDraw                      = False                 ; Do not draw plot
  opts_rvort@gsnFrame                     = False                 ; Do not adv. frame
  opts_rvort@gsnAddCyclic                 = False

  opts_rvort@mpLimitMode                  = "LatLon"

  if (zoom .eq. 1) then
   opts_rvort@mpMinLatF                   = lat_max - ar ; centre(cn,it-1,0) - 0.5
   opts_rvort@mpMinLonF                   = lon_max - ar ; centre(cn,it-1,1) - 0.5
   opts_rvort@mpMaxLatF                   = lat_max + ar ; centre(cn,it-1,0) + 0.5
   opts_rvort@mpMaxLonF                   = lon_max + ar ; centre(cn,it-1,1) + 0.5
  else
   opts_rvort@mpMinLatF                   = lat(0)
   opts_rvort@mpMinLonF                   = lon(0)
   opts_rvort@mpMaxLatF                   = lat(dimsizes(lat)-1)
   opts_rvort@mpMaxLonF                   = lon(dimsizes(lon)-1)
  end if

  opts_rvort@tmXTOn                       = "False"    ; No tickmarks on top x-axis
  opts_rvort@tmYROn                       = "False"    ; No tickmarks on right y-axis

  if (rad0 .le. 2.0) then 
   opts_rvort@tmXBMode			   = "Explicit"
   opts_rvort@tmXBValues	           = fspan(100.0, 150.0, 101)
   xb_labels0				   = sprintf("%0.1f~S~o~N~E",opts_rvort@tmXBValues)
   opts_rvort@tmXBLabels		   = xb_labels0

   opts_rvort@tmYLMode                     = "Explicit"
   opts_rvort@tmYLValues                   = fspan(0.0, 30.0, 61)
   yl_labels0                              = sprintf("%0.1f~S~o~N~N",opts_rvort@tmYLValues)
   opts_rvort@tmYLLabels                   = yl_labels0
  elseif (rad0 .gt. 2.0) then 
   opts_rvort@tmXBMode                     = "Explicit"
   opts_rvort@tmXBValues                   = fspan(100.0, 150.0, 26)
   xb_labels0                              = sprintf("%0.1f~S~o~N~E",opts_rvort@tmXBValues)
   opts_rvort@tmXBLabels                   = xb_labels0

   opts_rvort@tmYLMode                     = "Explicit"
   opts_rvort@tmYLValues                   = fspan(0.0, 30.0, 16)
   yl_labels0                              = sprintf("%0.1f~S~o~N~N",opts_rvort@tmYLValues)
   opts_rvort@tmYLLabels                   = yl_labels0
  end if 

  opts_rvort@mpGridAndLimbOn               = True       ; Grid lines on/off
  opts_rvort@mpGridLineDashPattern         = 2          ; Dash pattern
  opts_rvort@mpGridLatSpacingF             = 2.0        ; Spacing (latitude)
  opts_rvort@mpGridLonSpacingF             = 2.0        ; Spacing (longitude)
  opts_rvort@gsnMaximize                   = True 

; Pressure 
  opts_prs                                = True

  if (cyl .eq. 1 .or. cyl .eq. 2) then
   opts_prs@sfXArray                      = xpos(:,:)
   opts_prs@sfYArray                      = ypos(:,:)
  end if

  opts_prs@cnFillOn                       = False
  opts_prs@cnLineColor                    = "black"
  opts_prs@cnInfoLabelOn                  = False
  opts_prs@cnLineLabelsOn                 = True
  opts_prs@cnLineLabelInterval            = 2.0
  opts_prs@cnLevelSelectionMode           = "ExplicitLevels"
  opts_prs@cnLevels                       = ispan(200, 1040, 4)
  opts_prs@cnLineLabelPlacementMode       = "constant"
  opts_prs@cnLineLabelPerimOn             = False
  opts_prs@gsnContourLineThicknessesScale = 3.0
  opts_prs@tiMainString                   = ""
  opts_prs@gsnLeftString                  = ""
  opts_prs@gsnRightString                 = ""
  opts_prs@gsnDraw                        = False      ; Do not draw the plot
  opts_prs@gsnFrame                       = False      ; Do not advance the frame
  opts_prs@gsnMaximize                    = True       ; Maximise plot size 

;=================
; Plot the data
;=================

  if (var0 .eq. "vort") then 
   plane0 = vort_int
   plane1 = vort_azi
   plane2 = vort_plane
  elseif (var0 .eq. "spd") then 
   plane0 = spd_int
   plane1 = spd_azi
   plane2 = spd_plane
  elseif (var0 .eq. "th") then 
   plane0 = th_int
   plane1 = th_azi
   plane2 = th_plane
  end if 

; Single panel 
  if (cyl .eq. 1) then

   vort_plot  = gsn_csm_contour_map(wks,plane0,opts_rvort)          ; Relative vorticity (CYL)
   prs_plot   = gsn_csm_contour(wks,prs_int,opts_prs)               ; Pressure (CYL)
   if (prs .eq. 1) then 
    overlay(vort_plot, prs_plot)
   end if 

  elseif (cyl .eq. 2) then 

   vort_plot  = gsn_csm_contour_map(wks,plane1,opts_rvort)          ; Azimuthal average (CYL)

  else

   vort_plot  = gsn_csm_contour_map(wks,plane2,opts_rvort)          ; Relative vorticity (CART)
   prs_plot   = gsn_csm_contour(wks,prs_plane,opts_prs)             ; Pressure (CART)
   if (prs .eq. 1) then 
    overlay(vort_plot, prs_plot)
   end if 

  end if 

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

  markers    = new( (/dimsizes(m_colours)/), graphic)

  do i = 0, dimsizes(m_colours)-1
    mres0                    = True
    mres0@gsMarkerIndex      = m_indices(i)
    mres0@gsMarkerSizeF      = m_sizes(i)
    mres0@gsMarkerColor      = m_colours(i)
    mres0@gsMarkerThicknessF = m_thick(i)
    markers(i)               = gsn_add_polymarker(wks, vort_plot, centre_new(i,0,it-2,1),\
                                                  centre_new(i,0,it-2,0), mres0)
  end do

;===============================
; Add and customise a legend
;===============================

  lg_opts                            = True
  lg_opts@pmLegendDisplayMode        = "Always"          ; Display a legend
  lg_opts@lgAutoManage               = False             ; Gives us control
  lg_opts@lgLabelFont                = "Helvetica"       ; Label font
  lg_opts@lgLabelFontHeightF         = 0.060             ; Height of labels

; Customise labels and line colours
  lg_opts@lgItemType                 = "Markers"
  lg_opts@lgMarkerColors             = m_colours
  lg_opts@lgMarkerIndexes            = m_indices
  lg_opts@lgMarkerSizeF              = mres0@gsMarkerSizeF
  lg_opts@lgMarkerThicknessF         = mres0@gsMarkerThicknessF
  lsize                              = dimsizes(lg_opts@lgMarkerColors)

  lg_opts@lgLabelPosition            = "Right"           ; Label position within box
  lg_opts@lgItemPlacement            = "ExplicitPlacement"
  lg_opts@lgItemPositions            = fspan(0.05, 0.95, lsize)
  lg_opts@lgLeftMarginF		     = 0.01
  lg_opts@lgRightMarginF	     = 0.15

  lg_opts@vpWidthF                   = 0.200
  lg_opts@vpHeightF                  = 0.120
  lg_opts@lgPerimColor               = "black"
  lg_opts@lgPerimThicknessF          = 3.0
  lg_opts@lgPerimFill                = "SolidFill"
  lg_opts@lgPerimFillColor           = "white"

  lab_arr                            = (/"vort", "pres", \
                                         "pres_sm", "pres_sm2", \
                                         "vort_min", "windspeed_min", \
                                         "offline ("+cn0+")"/)
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

  annotate_1                         = gsn_add_annotation(vort_plot,lgd,am_opts)

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

; Attach outer circle to plot
  plot_outer = gsn_add_polyline(wks, vort_plot, xc, yc, lnres)

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
   plot_inner = new(size,graphic)

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
     plot_inner(count) = gsn_add_polyline(wks, vort_plot, xc, yc, lnres)

     count = count + 1 

   end do 

  end if 

; Tidy up 
  delete([/lnres, degrees, xcos, xsin, xc, yc/])

; Add information on the inner-core structure using 'gsn_text_ndc'
  title_res                   = True
  title_str                   = sym_str+" (T+"+it+")"
  title_res@txFontHeightF     = 0.015
  title_res@txFont            = "Helvetica"
  xpos_title                  = 0.12
  ypos_title                  = 0.96

  gsn_text_ndc(wks, title_str, xpos_title, ypos_title, title_res)

; Draw and frame the final plot 
  draw(vort_plot)
  frame(wks)

  delete([/opts_prs@cnLevels/])
  delete([/w_plane,u_plane,v_plane,spd_plane/])
  delete([/avo_plane,vort_plane,prs_plane,vort_pl_smth/])

 delete([/vort,u,v,u0,v0,latC,lonC,lat1C,lon1C/]) ; Tidy up 
 delete([/p,w,z,latD,lonD,lat1D,lon1D/])

 end do                ; End loop over input times (do it = 0, numINPUT_d-1)

; Stop script from overwriting text file unless all times are chosen
  t_int = tf0 - ts0

; Write out symmetry parameter values to text file 
  if (t_int .gt. 48) then 
 
   print("Working through all times; write to text file!")

   if (dnt .eq. 0) then 
    diri_out = "$ar/text/cp_sym_sc"+sc+"_"
    sym_out  = diri_out+dat+"_"+ens0+".txt"
   else 
    diri_out = "$ar/text/cp_sym_sc"+sc+"a_"
    sym_out  = diri_out+dat+"_"+ens0+".txt"
   end if 
   
   asciiwrite(sym_out, sym_arr)

  end if 

end 