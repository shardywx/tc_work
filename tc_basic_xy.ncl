; Script to plot basic variables (vorticity, horizontal wind, vertical velocity)
; on model levels for a given tropical cyclone

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/rad_circ.ncl"

; Run using:

; ncl 'opt="x11"' 'ens0="em11"' ts0=36 tf0=42 ls0=13 lf0=34 dist=3.0 
; 'dat="02T12"' calc=1 full=1 clr=1 tc=0 vvel0=0.8 tc_basic_xy.ncl

; 'opt'   = output file format ("pdf" or "x11")
; 'ens0'  = ensemble member (em00 -- em11)
; 'dist'  = size of box following storm (degrees)
; 'ts0'   = starting time index (0–119)
; 'tf0'   = ending time index (0–119)
; 'ls0'   = starting model level index (0–62)
; 'lf0'   = ending model level index (0–62)
; 'calc'  = centre calculated on each level (0) or equal to 950 hPa centre (1)
; 'full'  = all 10 panels plotted (1) or only a subset of 4 selected panels (0)
; 'clr'   = colour map for plots ("ncl_default" [1];  "bdr_extra" [2]; "amwg256" [3])
; 'tc'	  = overlay TC symbol on storm centre location (1) or leave blank (0)
; 'dat'	  = initialisation time ("02T12", "03T00", "03T12", "04T00")
; 'vvel0' = single vertical velocity contour in final plot (0.8, 1.5, etc)

begin

;==============================
; Prelash (reading in files)
;==============================

; Find 'pc' and 'pd' files (model height level data)
  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0
  fili_prefix = "20160702T1200Z_NPTK_4p4_L80_ra1t_"+ens0
  fili_c1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pc*.nc")
  fili_c      = diri+"/"+fili_c1
  fili_d1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pd*.nc")
  fili_d      = diri+"/"+fili_d1

; Calculate number of times in each file using built-in function ('nc_times')
  numINFO_c   = nc_times(fili_c)
  numINFO_d   = nc_times(fili_d)
  numINPUT_c  = numINFO_c[0]            ; Number of times in 'pc' stream
  time_c      = numINFO_c[1]            ; Array of times from 'pc' stream file
  numINPUT_d  = numINFO_d[0]            ; Number of times in 'pd' stream
  time_d      = numINFO_d[1]            ; Array of times from 'pd' stream file

; Calculate storm motion using built-in function ('st_rm')
; Argument 1  = ensemble simulation
; Argument 2  = path to input files
; Argument 3  = size of box following storm (degrees)

; Call external function 'st_rm'
  storm_rel = st_rm(dat, ens0, "$ar/text/", 3.0, 0, 0)
  u_cyc     = storm_rel[0]                 ; Zonal wind
  v_cyc     = storm_rel[1]                 ; Meridional wind
  vel_cyc   = storm_rel[2]                 ; Vector wind
  numTIMES  = storm_rel[3]                 ; Number of times in file
  lat_arr   = storm_rel[4]                 ; Storm-following domain subset (lat)
  lon_arr   = storm_rel[5]                 ; Storm-following domain subset (lon)
  centre    = storm_rel[6]                 ; Storm centre at each hour (lat/lon)

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

;==================================================================================
;;;;;;;;;;;;;;;;;;;;;;;;;;; Start of main part of script ;;;;;;;;;;;;;;;;;;;;;;;;;;
;==================================================================================

  print_clock("Starting main part of script...")

;=========================
; Start multiple loops
;=========================

    print_clock("Working on ensemble member: "+ens0)

    numTIMES    = 120                       ; Total number of times in all files
    ct          = 0                         ; Counter variable (time)
    llbox       = toint(dist*50)            ; Size of lat/lon arrays (see line below)
    dsize       = (/llbox,llbox/)           ; Size of lat/lon arrays (depends on 'dist')
    d0          = dsize(0)
    d1          = dsize(1)
    tsize	= (tf0 - ts0) + 1	    ; Number of times (depends on user input)

    centre_calc = new((/tsize,2,2/),float)  ; Storm centre (calculated in script) 

    setvalues NhlGetWorkspaceObjectId
      "wsMaximumSize" : 1000000000
    end setvalues

;==========================================================================
; Create arrays to hold values of key variables on multiple model levels
;==========================================================================

    u_arr       = new((/2,tsize,d0,d1/),float)	; Zonal wind	
    v_arr	= new((/2,tsize,d0,d1/),float)	; Meridional wind
    w_arr       = new((/2,tsize,d0,d1/),float)	; Vertical velocity
    vort_arr	= new((/2,tsize,d0,d1/),float)	; Relative vorticity
    avo_arr	= new((/2,tsize,d0,d1/),float)	; Absolute vorticity

;===============================
; Start loop over input files
;===============================

    c  = addfile(fili_c,"r")             ; Read in 'pc' stream [u,v,vort]
    d  = addfile(fili_d,"r")             ; Read in 'pd' stream [w,p,t,theta,kmh,tau]

;==================================
; Get the variables we will need
;==================================

   do it = ts0, tf0             ; Loop over times in file

    print("Working on time: "+time_str(it))
    itt = it + 1
    title_arr(it) = "Valid at "+time_str(it)+" (T+"+it+")"

    print("Time index = "+ct)

;=====================================================
; Read in variables from 'pc' stream (model levels)
;=====================================================

    t0    = lat_arr(it-1,0)
    t1    = lat_arr(it-1,1)
    n0    = lon_arr(it-1,0)
    n1    = lon_arr(it-1,1)
    print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

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
    u0  = c->$u_varname$(it,:,{t0:t1},{n0:n1}) ; Zonal wind
    v0  = c->$v_varname$(it,:,{t0:t1},{n0:n1}) ; Meridional wind

;===================================
; Calculate storm-relative winds
;===================================

    u = u0 - u_cyc(it-1)
    v = v0 - v_cyc(it-1)

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

;======================================================
; Read in variables from 'pd' stream (model levels)
;======================================================

  ; Read in basic variables from 'pd' stream
    lon1D  = d->longitude_1({n0:n1})                ; longitude [xxx grid points]
    lat1D  = d->latitude_1({t0:t1})                 ; latitude  [xxx grid points]
    lonD   = d->longitude({n0:n1})                  ; longitude [xxx grid points]
    latD   = d->latitude({t0:t1})                   ; latitude  [xxx grid points]
    hybD   = d->hybrid_ht(:)                        ; 63 model ('rho') levels

    hy     = hybD / 1000 ; Model level vales in km

  ; Read in temperature, pressure, vertical velocity and geopotential height
  ; All variables below -- [hybD | 63] * [latD | 150] * [lon1D | 151]
    t      = d->temp(it-2,:,{t0:t1},{n0:n1})        ; Temperature (K)
    p      = d->p(it-2,:,{t0:t1},{n0:n1})           ; Pressure (Pa)
    w      = d->dz_dt(it-2,:,{t0:t1},{n0:n1})       ; Vertical velocity (m s-1)
    z      = d->ht(it-2,:,{t0:t1},{n0:n1})          ; Geopotential height (m)
    pv     = d->field83(it-2,:,{t0:t1},{n0:n1})     ; Potential vorticity

  ; Convert units of pressure and geopotential height before continuing
    p      = p / 100          ; Convert to 'hPa'
    z      = z / 10           ; Convert to 'dam'
    pv     = pv * (10 ^ 6)    ; Convert to 'PVU' [10 ^ 6 K m2 s-1 kg-1]

  ; Add metadata
    p@units  = "hPa"
    z@units  = "dam"
    pv@units = "PVU"

    dx0      = lonD(1) - lonD(0) ; Grid spacing (longitude)
    dy0      = latD(1) - latD(0) ; Grid spacing (latitude)

;==========================
; Loop over model levels
;==========================

     m	     = 0			; Counter variable (model levels)

    do ilev  = ls0, lf0, (lf0 - ls0)	; Loop over selected model levels

     hy0     = sprintf("%0.0f",hy(ilev))
     hy1     = sprintf("%0.0f",hy(ls0))
     hy2     = sprintf("%0.0f",hy(lf0))
     print("Working on model level " + ilev+ " ("+hy0+" km AGL)")

   ; Read in other variables on model levels
     u_plane                = u(ilev,:,:)          ; Zonal wind
     v_plane                = v(ilev,:,:)          ; Meridional wind
     w_plane                = w(ilev,:,:)          ; Vertical velocity
     vort_plane             = vort(ilev,:,:)       ; Relative vorticity
     avo_plane              = avo(ilev,:,:)        ; Absolute vorticity

   ; Enter values into pre-existing, larger arrays for later plotting 
     u_arr(m,ct,:,:)        = u_plane
     v_arr(m,ct,:,:)	    = v_plane
     w_arr(m,ct,:,:)	    = w_plane
     vort_arr(m,ct,:,:)	    = vort_plane
     avo_arr(m,ct,:,:)	    = avo_plane

   ; Calculate relative vorticity maximum on model level
     vort_max               = max(vort_plane)
     vort_smth              = smth9_Wrap(vort_plane, 0.5, 0.5, True)
     vort_max_smth          = max(vort_smth)

   ; Find grid index corresponding to vorticity maximum 
     dims                   = dimsizes(vort_plane)
     vort1d                 = ndtooned(vort_plane)
     inds                   = ind_resolve(maxind(vort1d),dims)
     vort1d_sm              = ndtooned(vort_smth)
     inds_sm                = ind_resolve(maxind(vort1d_sm),dims)

   ;==================================================================
   ; Either calculate the centre, or read in pre-calculated values
   ;==================================================================

   if (calc .eq. 0) then

   ; Find latitude and longitude of storm at each time and on each model level
     lat_max                = t0 + (dy0 * inds(0,0))
     lon_max                = n0 + (dx0 * inds(0,1))

   ; Difference between 'lat_max'/'lon_max' and values read in using NCL function
     dlat                   = abs(lat_max - centre(it-1,0))
     dlon                   = abs(lon_max - centre(it-1,1))

   ; If our calculations of latitude and longitude deviate too much from the values
   ; we read in using the NCL function ('st_rm'), then calculate a new centre
   ; using the minimum pressure. Constrain the new calculation using the values
   ; from 'st_rm', to prevent spurious centres being identified and analysed.

     if (dlat .gt. 0.75 .or. dlon .gt. 0.75) then

      lt1 = centre(it-1,0) - 0.5
      lt2 = centre(it-1,0) + 0.5
      ln1 = centre(it-1,1) - 0.5
      ln2 = centre(it-1,1) + 0.5

      prs_plane0 = prs_plane({lt1:lt2},{ln1:ln2})       ; Smaller grid (1.0 x 1.0 deg)
      prs_min0   = min(prs_plane0)                      ; Find minimum on smaller grid
      dims_p0    = dimsizes(prs_plane0)                 ; Size of smaller grid
      prs_1d0    = ndtooned(prs_plane0)                 ; Create 1-D array
      inds_p0    = ind_resolve(minind(prs_1d0),dims_p0) ; Find index of minimum

      lat_max    = lt1 + (dy0 * inds_p0(0,0))
      lon_max    = ln1 + (dx0 * inds_p0(0,1))
      print("New centre: "+lat_max+" degrees N, "+lon_max+" degrees E")

    ; Input calculated values into array for later plotting 
      centre_calc(ct,m,0) = lat_max
      centre_calc(ct,m,1) = lon_max

     end if

   elseif (calc .eq. 1) then

      lat_max    = centre(it-1,0)
      lon_max    = centre(it-1,1)

    ; Input calculated values into array for later plotting
      centre_calc(ct,m,0) = lat_max
      centre_calc(ct,m,1) = lon_max

   end if

 ; Tidy up before next iteration
   delete([/vort_plane,u_plane,v_plane,w_plane,avo_plane/])

   m = m + 1	 ; Counter variable (model level)

   end do        ; End model level loop (do ilev = ls0, lf0)

   ct = ct + 1   ; Counter variable (time)

  end do     ; End time loop (do it = ts0, tf0))

; Tidy up before next iteration (all variables from 'pc' and 'pd' streams)
  delete([/lonC,latC,lon1C,lat1C,u,v,vort,avo,f0,f1/])
  delete([/lonD,latD,lon1D,lat1D,t,p,w,z/])

;================================================
; Begin second loop to plot and write out data
;================================================

  print_clock("Producing panel plots...")

  ct    = 0		                ; Counter variable (time)

  do it = ts0, tf0			; Loop over all times

  print("Working on time: "+time_str(it))
  print("ct = "+ct+" ; it = "+it)

; Read in latitude and longitude information again (metadata)
  t0    = lat_arr(it-1,0)
  t1    = lat_arr(it-1,1)
  n0    = lon_arr(it-1,0)
  n1    = lon_arr(it-1,1)
  print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

; Read in latitude and longitude from 'pc' and 'pd' streams
  lonC   = c->longitude({n0:n1})                ; longitude ['d1' grid points]
  latC   = c->latitude({t0:t1})                 ; latitude  ['d0' grid points]

  lon1D  = d->longitude_1({n0:n1})              ; longitude [xxx grid points]
  latD   = d->latitude({t0:t1})                 ; latitude  [xxx grid points]

; Read in data on lower model level
  vort_plane0 = vort_arr(0,ct,:,:)
  u_plane0    = u_arr(0,ct,:,:)
  v_plane0    = v_arr(0,ct,:,:)
  w_plane0    = w_arr(0,ct,:,:)

; Read in data on higher model level
  vort_plane1 = vort_arr(1,ct,:,:)
  u_plane1    = u_arr(1,ct,:,:)
  v_plane1    = v_arr(1,ct,:,:)
  w_plane1    = w_arr(1,ct,:,:)

; Add metadata to arrays so we plot the correct region
  vort_plane0&latitude  = latC
  vort_plane0&longitude = lonC
  v_plane0&latitude     = latC
  v_plane0&longitude    = lonC
  u_plane0&latitude     = latC
  u_plane0&longitude    = lonC
  w_plane0&latitude     = latD
  w_plane0&longitude_1  = lon1D  

  vort_plane1&latitude  = latC
  vort_plane1&longitude = lonC
  v_plane1&latitude     = latC
  v_plane1&longitude    = lonC
  u_plane1&latitude     = latC
  u_plane1&longitude    = lonC
  w_plane0&latitude     = latD
  w_plane0&longitude_1  = lon1D

;============================================
; Define output file path and plot details
;============================================

; Output file location and type
  if (clr .eq. 4) then 
   output = "$sam/nepartak/images/cyl_coords/m18_xy_"+hy1+"km_"+hy2+"km_"+ens0+"_"+time_arr(it)  
  else
   output = "$sam/nepartak/images/cyl_coords/basic_xy_"+hy1+"km_"+hy2+"km_"+ens0+"_"+time_arr(it)
  end if 
  wks = gsn_open_wks(opt,output)

; Set up panel plot (plot either 1 or 2 panels)

  if (full .eq. 1) then
   panel = new(2,graphic)
  elseif (full .eq. 0) then
   panel = new(1,graphic)
  end if

  ; Read in desired colourmap
    if (clr.eq.1) then
     gsn_define_colormap(wks,"ncl_default")
    elseif (clr.eq.2) then
     gsn_define_colormap(wks,"bdr_extra")
    elseif (clr.eq.3) then
     gsn_define_colormap(wks,"amwg256")
    elseif (clr.eq.4) then 
     gsn_define_colormap(wks,"GMT_polar")
    elseif (clr.eq.5) then 
     gsn_define_colormap(wks,"prcp_new")
    end if 

;==========================
; Options for plotting
;==========================

  ; Relative vorticity
    opts_vort                             = True
    opts_vort@cnFillOn                    = True
    opts_vort@cnLineLabelInterval         = 2.0
    opts_vort@cnLineLabelFontHeightF      = 0.012
    opts_vort@cnLineLabelBackgroundColor  = "transparent"
    opts_vort@cnLineLabelPlacementMode    = "constant"
    opts_vort@cnLinesOn                   = False
    opts_vort@cnInfoLabelOn               = False
    opts_vort@cnLevelSelectionMode        = "ExplicitLevels"

    if (clr.eq.1) then
     opts_vort@cnFillColors               = (/66,78,90,102,120,0,0,\
                                              136,152,168,180,188,196,\
                                              204,216,228,240,252/)
     opts_vort@cnLevels                   = (/-3., -2., -1., -0.5, -0.2, \
                                               0, 0.2, 0.5, 1., 2., 3., \
                                               5., 10., 20., 30., 40., 50./)

    elseif (clr.eq.2) then
     opts_vort@cnFillColors               = (/8,9,10,11,12,0,0,\
                                              15,16,17,18,19,20,21,22,23,24,25/)
     opts_vort@cnLevels                   = (/-3., -2., -1., -0.5, -0.2, \
                                               0, 0.2, 0.5, 1., 2., 3., \
                                               5., 10., 20., 30., 40., 50./)
    elseif (clr.eq.3) then
     opts_vort@cnFillColors               = (/47,56,68,80,92,0,0,\
                                              150,158,168,180,188,196,\
                                              204,216,228,240,252/)
     opts_vort@cnLevels                   = (/-3., -2., -1., -0.5, -0.2, \
                                               0, 0.2, 0.5, 1., 2., 3., \
                                               5., 10., 20., 30., 40., 50./)
    elseif (clr.eq.4) then
     opts_vort@cnFillColors               = (/9,0,13,16,20/)
     opts_vort@cnLevels                   = (/-5,5,10,20/)
    elseif (clr .eq. 5) then 
     opts_vort@cnLevels                   = (/5., 10., 20., 30., 40., 50., \
                                              60., 80., 100., 120., 150./)
     opts_vort@cnFillColors               = (/0,4,5,6,7,8,9,10,11,12,13,14/)
    end if

    opts_vort@gsnPaperOrientation          = "landscape"
    opts_vort@tiMainString                 = ""
    opts_vort@tiMainFontHeightF            = 0.0125
    opts_vort@gsnLeftString                = ""
    opts_vort@gsnRightString               = ""
    opts_vort@gsnMaximize                  = True
    opts_vort@lbLabelBarOn                 = False
    opts_vort@gsnDraw                      = False
    opts_vort@gsnFrame                     = False

    opts_vort@pmLabelBarWidthF             = 0.6
    opts_vort@pmLabelBarHeightF            = 0.15
    opts_vort@lbLabelFontHeightF           = 0.01
    opts_vort@lbPerimOn                    = False
    opts_vort@lbLabelFont                  = "Helvetica"
    opts_vort@lbTitleOn                    = True
    opts_vort@lbTitleString 		   = "Relative vorticity (10~S~-4~N~ s~S~-1~N~)"
    opts_vort@lbTitleFontHeightF           = 0.015
    opts_vort@lbTitleDirection             = "Across"

  ; Additional plotting resources
    opts_vort@mpDataBaseVersion            = "Ncarg4_1" ; More recent database
    opts_vort@mpDataSetName                = "Earth..4" ; High resolution
    opts_vort@mpOutlineBoundarySets        = "National" ; National borders
    opts_vort@mpGeophysicalLineColor       = "black"    ; Colour borders black
    opts_vort@mpGeophysicalLineThicknessF  = 1.0        ; Border line thickness
    opts_vort@mpGridAndLimbOn              = False      ; Turn on lat/lon lines
    opts_vort@pmTickMarkDisplayMode        = "Always"   ; Turn on map tickmarks
    opts_vort@tmXBMajorLengthF             = 0.005      ; Change tickmark length
    opts_vort@tmXTOn                       = "False"    ; No tickmarks on top x-axis
    opts_vort@tmYROn                       = "False"    ; No tickmarks on right y-axis
    opts_vort@gsnMaximize                  = True       ; Maximise plot size
    opts_vort@gsnAddCyclic                 = False

    radius                                 = 1.0	; Radius of plot (degrees)

    opts_vort@mpLimitMode                  = "Corners"
    opts_vort@mpLeftCornerLatF             = centre(it-1,0)-radius
    opts_vort@mpLeftCornerLonF             = centre(it-1,1)-radius
    opts_vort@mpRightCornerLatF            = centre(it-1,0)+radius
    opts_vort@mpRightCornerLonF            = centre(it-1,1)+radius
    opts_vort@gsnDraw                      = False      ; Do not draw the plot
    opts_vort@gsnFrame                     = False      ; Do not advance the frame

  ; Horizontal wind vector resources
    opts_vec                               = True
    opts_vec@gsnDraw			   = False
    opts_vec@gsnFrame			   = False
    opts_vec@gsnAddCyclic		   = False

    opts_vec@vcGlyphStyle                  = "LineArrow"   ; draw windbarbs
    opts_vec@vcRefMagnitudeF               = 15.0          ; ref vec = max vector
    opts_vec@vcRefLengthF                  = 0.04          ; length of reference vector
    opts_vec@vcMinDistanceF                = 0.08          ; distance between wind barbs
    opts_vec@vcLineArrowThicknessF         = 2.0           ; windbarb line thickness
    opts_vec@vcLineArrowColor              = "black"       ; windbarb colour
    opts_vec@gsnMaximize                   = True          ; Maximise plot size
    opts_vec@gsnRightString                = ""
    opts_vec@gsnLeftString                 = ""
    opts_vec@tiXAxisString                 = ""
    opts_vec@vcRefAnnoOrthogonalPosF       = 0.02	   ; Y-location of ref vec

  ; Vertical velocity resources
    opts_vvel                              = True
    opts_vvel@cnFillOn			   = False
    opts_vvel@cnInfoLabelOn		   = False
    opts_vvel@gsnDraw			   = False
    opts_vvel@gsnFrame			   = False
    opts_vvel@gsnAddCyclic		   = False
    opts_vvel@tiMainString		   = ""
    opts_vvel@gsnLeftString                = ""
    opts_vvel@gsnRightString               = ""
    opts_vvel@cnLineThicknessF		   = 4.0
    opts_vvel@cnLineLabelsOn		   = False
    opts_vvel@cnLevelSelectionMode	   = "ExplicitLevels"
    opts_vvel@cnLevels			   = vvel0

  ; Yellow vertical velocity contours if we choose Roger Smith's colour map
    if (clr .eq. 4) then 
     opts_vvel@cnLineColor		   = "yellow"
    elseif (clr .eq. 5) then 
     opts_vvel@cnLineColor		   = "black"
    else
     opts_vvel@cnLineColor                 = "gray45"
    end if 

  ; Resources for mid-tropospheric vertical velocity
    opts_vvel1	  		   	   = opts_vvel
    opts_vvel1@cnLineColor		   = "royalblue1"

;=============
; Plot data 
;=============

    vort_plane0 = vort_plane0 * (10 ^ 4)
    vort_plane0@units = "10~S~-4~N~ s~S~-1~N~"

    vort_plane1 = vort_plane1 * (10 ^ 4)
    vort_plane1@units = "10~S~-4~N~ s~S~-1~N~"

  ;=====================
  ; Lower model level
  ;=====================

  ; Plot relative vorticity
    plot_vort0  = gsn_csm_contour_map(wks,vort_plane0,opts_vort)

  ; Plot horizontal wind vectors
    plot_vec0   = gsn_csm_vector(wks,u_plane0,v_plane0,opts_vec)

  ; Plot vertical velocity
    plot_vvel0  = gsn_csm_contour(wks,w_plane0,opts_vvel)

    overlay(plot_vort0,plot_vec0)
    overlay(plot_vort0,plot_vvel0)

  ;=====================
  ; Upper model level
  ;=====================

  ; Plot relative vorticity
    plot_vort1  = gsn_csm_contour_map(wks,vort_plane1,opts_vort)

  ; Plot horizontal wind vectors
    plot_vec1   = gsn_csm_vector(wks,u_plane1,v_plane1,opts_vec)

  ; Plot vertical velocity
    plot_vvel1  = gsn_csm_contour(wks,w_plane1,opts_vvel)

    overlay(plot_vort1,plot_vec1)
    overlay(plot_vort1,plot_vvel1)

  ; Plot either one or two panel(s)

    if (full .eq. 0) then    	    	   ; One panel

     plot_vvel1 = gsn_csm_contour(wks,w_plane1,opts_vvel1)
     overlay(plot_vort0, plot_vvel1) 
     panel(0)   = plot_vort0   

    elseif (full .eq. 1) then		   ; Two panels 

     panel(0) = plot_vort0
     panel(1) = plot_vort1

    end if 

;==============================================
; Add the NH tropical cyclone symbol to plot
;==============================================

    marker    = new(dimsizes(panel),graphic)

  ; Add NH tropical cyclone symbol to index
    tc0 = NhlNewMarker(wks, "p", 137, 0.0, 0.0, 1.0, 1.0, 0.0)

    mres0               = True
    mres0@gsMarkerIndex = tc0         ; TC symbol
    mres0@gsMarkerSizeF = 20.0        ; Marker size
    mres0@gsMarkerColor = "red3"      ; Marker colour
    mres0@gsMarkerThicknessF = 5.0    ; Marker thickness (increase)

  ; Add annotation if flag is set to '1'
    if (tc .eq. 1) then 
     do m = 0, dimsizes(panel)-1
      marker(m) = gsn_add_polymarker(wks,panel(m),centre(it-1,1),centre(it-1,0),mres0)
     end do
    end if 

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    optsP                       = True                     ; Panel plot resources
    optsP@gsnFrame              = False                    ; Do not advance the frame
    optsP@gsnPanelLabelBar      = True                     ; Panel labelbar
    optsP@pmLabelBarWidthF      = 0.6                      ; Label bar width
    optsP@pmLabelBarHeightF     = 0.15                     ; Label bar height
    optsP@lbLabelFontHeightF    = 0.01                     ; Label bar font height
    optsP@lbLabelFont           = "Helvetica"              ; Choose font
    optsP@lbPerimOn             = False                    ; Perimeter on label bar
    optsP@lbTitleString         = \
				"Relative vorticity (10~S~-4~N~ s~S~-1~N~)"
    optsP@lbTitlePosition       = "Bottom"                 ; Title position
    optsP@lbTitleFontHeightF    = 0.015                    ; Title font height
    optsP@lbBoxEndCapStyle      = "TriangleBothEnds"       ; Label bar style
    optsP@pmLabelBarOrthogonalPosF = -0.002                ; Move lb up/down (h/l)
    optsP@pmLabelBarParallelPosF   = 0.04		   ; Move lb right/left (h/l)
    optsP@lbTitleOffsetF        = 0                        ; x-axis title up/down (h/l)
    optsP@txString              = title_arr(it)            ; Plot title string
    optsP@txPosXF		= 0.525			   ; Title pos (x-direction)

    if (full .eq. 1) then
     optsP@gsnPanelFigureStrings = (/"a) Relative vorticity at "+hy1+" km", \
     				     "b) Relative vorticity at "+hy2+" km"/)
     optsP@gsnPanelFigureStringsFontHeightF = 0.0125
    else if (full .eq. 0) then
     optsP@gsnPanelFigureStrings = (/"a) Relative vorticity at "+hy1+" km"/)
     optsP@gsnPanelFigureStringsFontHeightF = 0.0125
    end if
    end if

    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.88
    optsP@gsnPanelBottom        = 0.12
    optsP@amJust                = "TopRight"

    if (full .eq. 1) then
     gsn_panel(wks,panel,(/1,2/),optsP)             ; Draw as a single plot
     frame(wks)
    else if (full .eq. 0) then
     gsn_panel(wks,panel,(/1,1/),optsP)             ; Draw as a single plot
     frame(wks)
    end if
    end if

  ; Tidy up before next iteration
    delete([/vort_plane0,u_plane0,v_plane0,w_plane0/])
    delete([/vort_plane1,u_plane1,v_plane1,w_plane1/])
    delete([/lonC,latC,lon1D,latD/])

    ct = ct + 1	  ; Counter variable (time)

   end do     ; End time loop (do it = ts0, tf0))

end