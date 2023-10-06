; Compare storm-tracking methods and overlay on relative vorticity or windspeed 

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/setup_cyl.ncl"

; ncl dat=\"02T12\" opt=\"x11\" ens0=\"em11\" rad0=2.0 dist=3.0 vc=0 sm=0
; zoom=1 lgd0=\"tr\" plt=1 cn0=\"geo_sm\" ts0=51 tf0=52 sc=4 lay=1 cyl=1 dnt=0 tc_fc_uvmin.ncl

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
                  "potential temperature maximum",\               ; sc = 6
                  "geopotential"/)                                ; sc = 7

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

;   printVarSummary(u0)

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
    p      = d->p(it-2,13:25,{t0:t1},{n0:n1})      ; Pressure (Pa)
    w      = d->dz_dt(it-2,13:25,{t0:t1},{n0:n1})  ; Vertical velocity (m s-1)
    z      = d->ht(it-2,13:25,{t0:t1},{n0:n1})     ; Geopotential height (m)
   else
    p      = d->p(it-2,ilev,{t0:t1},{n0:n1})       ; Pressure (Pa)
    w      = d->dz_dt(it-2,ilev,{t0:t1},{n0:n1})   ; Vertical velocity (m s-1) 
    z      = d->ht(it-2,ilev,{t0:t1},{n0:n1})      ; Geopotential height (m) 
   end if 

;   printVarSummary(p)

 ; Convert units of pressure and geopotential height before continuing 
   p      = p / 100
   z      = z / 10

 ; Add metadata
   p@units = "hPa"
   z@units = "dam"

;==========================
; Loop over model levels 
;==========================

   hy0 = sprintf("%0.0f",hybD(13))
   hy1 = sprintf("%0.0f",hybD(25))
   print("Working on layer between " + hy0+ " and "+hy1+" m AGL")

 ; Read in variables on model levels 
   if (lay .eq. 1) then 

  ; 'pc' stream
    u_plane      = dim_avg_n_Wrap( u(:,:,:), 0)
    v_plane      = dim_avg_n_Wrap( v(:,:,:), 0)
    vort_plane   = dim_avg_n_Wrap( vort(:,:,:), 0)

  ; 'pd' stream 
    vvel_plane   = dim_avg_n_Wrap( w(:,:,:), 0)
    prs_plane    = dim_avg_n_Wrap( p(:,:,:), 0)
   else

  ; 'pc' stream
    u_plane      = u(ilev,:,:)
    v_plane      = v(ilev,:,:)
    vort_plane   = vort(ilev,:,:)

  ; 'pd' stream
    vvel_plane   = w(ilev,:,:)
    prs_plane    = p(ilev,:,:)

   end if 

 ; Calculate vector windspeed and smooth 
   spd_plane    = sqrt( (u_plane ^ 2) + (v_plane ^ 2) )
   copy_VarCoords(u_plane, spd_plane)
   vort_pl_smth = smth9_Wrap(vort_plane, 0.5, 0.5, True)
   prs_pl_smth  = smth9_Wrap(prs_plane, 0.5, 0.5, True)
   prs_pl_smth2 = smth9_Wrap(prs_pl_smth, 0.5, 0.5, True)

;====================================================
; Calculate storm centre using variety of methods
;====================================================

 ; Use already-calculated centre as a base for further calculations
   lt1		  = centre(cn,it-1,0) - 1.0
   lt2		  = centre(cn,it-1,0) + 1.0
   ln1		  = centre(cn,it-1,1) - 1.0
   ln2		  = centre(cn,it-1,1) + 1.0

 ; RELATIVE VORTICITY (maximum)
   vort_plane0    = vort_plane({lt1:lt2},{ln1:ln2})
   vort_max0	  = max(vort_plane0)
   dims_v0        = dimsizes(vort_plane0)
   vort_1d0       = ndtooned(vort_plane0)

 ; Find index of minimum --> (0,0) = latitude, (0,1) = longitude
   inds_v0        = ind_resolve(maxind(vort_1d0),dims_v0)

 ; Retrieve lat/lon information from grid subset
   ltN = vort_plane0&latitude
   lnN = vort_plane0&longitude

   lat_max_v0 	  = ltN(0) + (dy * inds_v0(0,0))
   lon_max_v0 	  = lnN(0) + (dx * inds_v0(0,1))
   print("Centre (vort): "+lat_max_v0+" degrees N, "+lon_max_v0+" degrees E")

   centre_new(0,0,it-2,0) = lat_max_v0
   centre_new(0,0,it-2,1) = lon_max_v0
   delete([/ltN, lnN/])

 ; PRESSURE (minimum)
   prs_plane0     = prs_plane({lt1:lt2},{ln1:ln2})
   prs_min0       = min(prs_plane0)
   dims_p0        = dimsizes(prs_plane0)
   prs_1d0        = ndtooned(prs_plane0)

 ; Find index of minimum --> (0,0) = latitude, (0,1) = longitude
   inds_p0        = ind_resolve(minind(prs_1d0),dims_p0)

 ; Retrieve lat/lon information from grid subset
   ltN            = prs_plane0&latitude
   lnN            = prs_plane0&longitude_1

   lat_max_p0     = ltN(0) + (dy * inds_p0(0,0))
   lon_max_p0     = lnN(0) + (dx * inds_p0(0,1))   
   print("Centre (prs): "+lat_max_p0+" degrees N, "+lon_max_p0+" degrees E")

   centre_new(1,0,it-2,0) = lat_max_p0
   centre_new(1,0,it-2,1) = lon_max_p0
   delete([/ltN, lnN/])

 ; SMOOTHED PRESSURE (minimum)
   prs_smth       = prs_pl_smth({lt1:lt2},{ln1:ln2})
   prs_min_smth	  = min(prs_smth)
   dims_p1	  = dimsizes(prs_smth)
   prs_1d1	  = ndtooned(prs_smth)

 ; Find index of minimum --> (0,0) = latitude, (0,1) = longitude
   inds_p1	  = ind_resolve(minind(prs_1d1),dims_p1)

 ; Retrieve lat/lon information from grid subset
   ltN            = prs_smth&latitude
   lnN            = prs_smth&longitude_1

   lat_max_p1     = ltN(0) + (dy * inds_p1(0,0))
   lon_max_p1     = lnN(0) + (dx * inds_p1(0,1))
   print("Centre (prs_smth): "+lat_max_p1+" degrees N, "+lon_max_p1+" degrees E")

   centre_new(2,0,it-2,0) = lat_max_p1
   centre_new(2,0,it-2,1) = lon_max_p1   
   delete([/ltN, lnN/])

 ; EXTRA SMOOTHED PRESSURE (minimum)
   prs_smth2      = prs_pl_smth2({lt1:lt2},{ln1:ln2})
   prs_min_smth2  = min(prs_smth2)
   dims_p2        = dimsizes(prs_smth2)
   prs_1d2        = ndtooned(prs_smth2)

 ; Find index of minimum --> (0,0) = latitude, (0,1) = longitude
   inds_p2        = ind_resolve(minind(prs_1d2),dims_p2)

 ; Retrieve lat/lon information from grid subset
   ltN            = prs_smth2&latitude
   lnN            = prs_smth2&longitude_1

   lat_max_p2     = ltN(0) + (dy * inds_p2(0,0))
   lon_max_p2     = lnN(0) + (dx * inds_p2(0,1))
   print("Centre (prs_smth_2): "+lat_max_p2+" degrees N, "+lon_max_p2+" degrees E")

   centre_new(3,0,it-2,0) = lat_max_p2
   centre_new(3,0,it-2,1) = lon_max_p2 
   delete([/ltN, lnN/])

 ; DURING RINGLIKE PHASE, FIND THE CENTRE OF THE 'DOUGHNUT'
 ; Use centre calculated using smoothed pressure (above) as base
   r0             = 0.10
   lt3            = lat_max_p2 - r0
   lt4            = lat_max_p2 + r0
   ln3            = lon_max_p2 - r0
   ln4            = lon_max_p2 + r0
   
   vort_pl_new    = vort_plane({lt3:lt4},{ln3:ln4})
   vort_min_new   = min(vort_pl_new)
   dims_v3        = dimsizes(vort_pl_new)
   vort_1d3       = ndtooned(vort_pl_new)

 ; Find index of minimum --> (0,0) = latitude, (0,1) = longitude
   inds_v3        = ind_resolve(minind(vort_1d3),dims_v3)

 ; Retrieve lat/lon information about grid subset
   ltN            = vort_pl_new&latitude
   lnN            = vort_pl_new&longitude

   lat_min_v3     = ltN(0) + (dy * inds_v3(0,0))
   lon_min_v3     = lnN(0) + (dx * inds_v3(0,1))
   print("Centre (vort_min): "+lat_min_v3+" degrees N, "+lon_min_v3+" degrees E")

   centre_new(4,0,it-2,0) = lat_min_v3
   centre_new(4,0,it-2,1) = lon_min_v3
   delete([/ltN, lnN/])

 ; WINDSPEED MINIMUM 
   spd_pl_new     = spd_plane({lt3:lt4},{ln3:ln4})          ; Smaller grid 
   spd_min_new    = min(spd_pl_new)                         ; Find minimum on smaller grid
   dims_s3        = dimsizes(spd_pl_new)                    ; Size of smaller grid
   spd_1d3        = ndtooned(spd_pl_new)                    ; Create 1D array

 ; Find index of minimum --> (0,0) = latitude, (0,1) = longitude
   inds_s3        = ind_resolve(minind(spd_1d3),dims_s3)    ; Find index of minimum

 ; Retrieve lat/lon information from subset of grid 
   ltN = spd_pl_new&latitude
   lnN = spd_pl_new&longitude

 ; Find lat/lon corresponding to windspeed minimum on this grid subset
   lat_min_s3     = ltN(0) + (dy * inds_s3(0,0))
   lon_min_s3     = lnN(0) + (dx * inds_s3(0,1))
   print("Centre (spd_min): "+lat_min_s3+" degrees N, "+lon_min_s3+" degrees E")

   centre_new(5,0,it-2,0) = lat_min_s3
   centre_new(5,0,it-2,1) = lon_min_s3
   delete([/ltN, lnN/])

 ; PRE-CALCULATED STORM POSITION 
   centre_new(6,0,it-2,0) = centre(cn,it-1,0)
   centre_new(6,0,it-2,1) = centre(cn,it-1,1)

 ; Calculate geopotential height max/min values for later contour plotting
   p_max0	  = toint(max(prs_plane))
   p_min0	  = toint(prs_min0)
   p_avg0	  = ( (p_max0 + p_min0) / 2)

;=====================================================================
; Call external function to switch to cylindrical coordinate system
;=====================================================================

 ; Call 'setup_cyl' and output pressure/vorticity on cylindrical grid 
   cyl_arr  = setup_cyl(it, dist, t0, t1, n0, n1, lonC, latC, lonD, latD, \
                        "spd", spd_plane, prs_plane, centre_new, sc, sm, ar, nr)
   prs_int  = cyl_arr[0]
   spd_int  = cyl_arr[1]
   spd_azi  = cyl_arr[2]
   xpos     = cyl_arr[3]
   ypos     = cyl_arr[4]
   lat_max  = cyl_arr[5]
   lon_max  = cyl_arr[6]
   radii    = cyl_arr[7]
   rad_size = cyl_arr[8]

 ; How symmetric is our storm? 
 ; EDIT -- customise this part of script; how does our choice of radii affect the results?
   sym_arr(it) = avg( spd_azi(0,0) + spd_azi(1,0) ) / avg ( spd_azi(2,0) + spd_azi(3,0) )
   if (sym_arr(it) .gt. 1) then 
    print("Storm is in asymmetric phase (ratio = "+sym_arr(it)+")")   
   else
    print("Storm in in symmetric phase (ratio = "+sym_arr(it)+")")
   end if 

 ; If we're in the symmetric phase, put the storm centre in the middle of the vorticity ring 
   if (dnt .eq. 1 .and. sym_arr(it) .lt. 1) then 

    print("Re-calculating storm centre to fit ringlike vorticity profile...")

  ; Select the 'vort_min' method of storm tracking 
    sc = 4

  ; Call the external function again
    cyl_arr  = setup_cyl(it, dist, t0, t1, n0, n1, lonC, latC, lonD, latD, \
                         "spd", spd_plane, prs_plane, centre_new, sc, sm, ar, nr)
    prs_int  = cyl_arr[0]
    spd_int  = cyl_arr[1]
    spd_azi  = cyl_arr[2]    
    xpos     = cyl_arr[3]
    ypos     = cyl_arr[4]
    lat_max  = cyl_arr[5]
    lon_max  = cyl_arr[6]
    radii    = cyl_arr[7]
    rad_size = cyl_arr[8]

   end if 

 ; EDIT -- customise further
 ; (2) maybe add SYM/ASYM information on panel itself 
 ; (3) plot raw data alongside the azimuthal mean 

 ; Output file location and type
   output = "$sam/nepartak/images/basic/fc_uvmin_sc"+sc+"_"+dat+"_"+ens0+"_"+time_arr(it-2)
   wks    = gsn_open_wks(opt,output)

 ; Load colour table 
   gsn_define_colormap(wks,"prcp_new")

;==========================
; Options for plotting
;==========================

; Horizontal windspeed
  opts_uv                              = True
  opts_uv@cnFillOn                     = True

  if (cyl .eq. 1 .or. cyl .eq. 2) then 
   opts_uv@sfXArray                    = xpos(:,:)
   opts_uv@sfYArray                    = ypos(:,:)
  end if 

  opts_uv@cnLineLabelInterval          = 2.0
  opts_uv@cnLineLabelFontHeightF       = 0.012
  opts_uv@cnLineLabelBackgroundColor   = "transparent"
  opts_uv@cnLineLabelPlacementMode     = "constant"
  opts_uv@cnLinesOn                    = False ; Contour lines off
  opts_uv@cnInfoLabelOn                = False ; Contour labels off
  opts_uv@cnLevelSelectionMode         = "ExplicitLevels"
  opts_uv@cnLevels                     = (/15,20,25,30,35,40,45,50,55,60/)
  opts_uv@cnFillColors                 = (/0,4,5,6,7,8,9,10,11,12,13,14,15/)

  opts_uv@gsnPaperOrientation          = "landscape"
  opts_uv@tiMainString                 = ""
  opts_uv@tiMainFontHeightF            = 0.0125
  opts_uv@gsnLeftString                = ""
  opts_uv@gsnRightString               = ""

  opts_uv@lbLabelBarOn                 = True
  opts_uv@lbBoxEndCapStyle             = "TriangleBothEnds"    ; Labelbar end shape
  opts_uv@lbLabelFontHeightF           = 0.0125                ; Labelbar font size
  opts_uv@lbLabelFont                  = "Helvetica"           ; Labelbar font
  opts_uv@lbTitleString                = "Layer-averaged ("+hy0+"-"+hy1+\
                                           " m) ~C~     horizontal windspeed (m s~S~-1~N~)"
  opts_uv@lbTitlePosition              = "Bottom"              ; Title position
  opts_uv@lbTitleFontHeightF           = 0.015                 ; Title font height
  opts_uv@pmLabelBarOrthogonalPosF     = 0.1                   ; Labelbar position
  opts_uv@lbTitleOffsetF               = 0                     ; Lb title up/down
  opts_uv@lbPerimOn                    = False                 ; Perimeter on/off
  opts_uv@gsnDraw                      = False                 ; Do not draw plot
  opts_uv@gsnFrame                     = False                 ; Do not adv. frame
  opts_uv@gsnAddCyclic                 = False

  opts_uv@mpLimitMode                  = "LatLon"

  if (zoom .eq. 1) then
   opts_uv@mpMinLatF                   = lat_max - 0.5 ; centre(cn,it-1,0) - 0.5
   opts_uv@mpMinLonF                   = lon_max - 0.5 ; centre(cn,it-1,1) - 0.5
   opts_uv@mpMaxLatF                   = lat_max + 0.5 ; centre(cn,it-1,0) + 0.5
   opts_uv@mpMaxLonF                   = lon_max + 0.5 ; centre(cn,it-1,1) + 0.5
  else
   opts_uv@mpMinLatF                   = lat(0)
   opts_uv@mpMinLonF                   = lon(0)
   opts_uv@mpMaxLatF                   = lat(dimsizes(lat)-1)
   opts_uv@mpMaxLonF                   = lon(dimsizes(lon)-1)
  end if

  opts_uv@tmXTOn                       = "False"    ; No tickmarks on top x-axis
  opts_uv@tmYROn                       = "False"    ; No tickmarks on right y-axis

  if (rad0 .le. 2.0) then 

   opts_uv@tmXBMode			   = "Explicit"
   opts_uv@tmXBValues    	           = fspan(100.0, 150.0, 101)
   xb_labels0				   = sprintf("%0.1f~S~o~N~E",opts_uv@tmXBValues)
   opts_uv@tmXBLabels	        	   = xb_labels0

   opts_uv@tmYLMode                        = "Explicit"
   opts_uv@tmYLValues                      = fspan(0.0, 30.0, 61)
   yl_labels0                              = sprintf("%0.1f~S~o~N~N",opts_uv@tmYLValues)
   opts_uv@tmYLLabels                      = yl_labels0

  elseif (rad0 .gt. 2.0) then 

   opts_uv@tmXBMode                        = "Explicit"
   opts_uv@tmXBValues                      = fspan(100.0, 150.0, 26)
   xb_labels0                              = sprintf("%0.1f~S~o~N~E",opts_uv@tmXBValues)
   opts_uv@tmXBLabels                      = xb_labels0

   opts_uv@tmYLMode                        = "Explicit"
   opts_uv@tmYLValues                      = fspan(0.0, 30.0, 16)
   yl_labels0                              = sprintf("%0.1f~S~o~N~N",opts_uv@tmYLValues)
   opts_uv@tmYLLabels                      = yl_labels0

  end if 

  opts_uv@mpGridAndLimbOn               = True       ; Grid lines on/off
  opts_uv@mpGridLineDashPattern         = 2          ; Dash pattern
  opts_uv@mpGridLatSpacingF             = 2.0        ; Spacing (latitude)
  opts_uv@mpGridLonSpacingF             = 2.0        ; Spacing (longitude)
  opts_uv@gsnMaximize                   = True 

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
  opts_prs@cnLevels                       = ispan(p_min0, p_max0, 4)
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

; Single panel 
  if (cyl .eq. 1) then
   spd_plot   = gsn_csm_contour_map(wks,spd_int,opts_uv)         ; Horizontal windspeed (CYL)
   prs_plot   = gsn_csm_contour(wks,prs_int,opts_prs)            ; Pressure (CYL)
   if (prs .eq. 1) then 
    overlay(spd_plot, prs_plot)
   end if 
  elseif (cyl .eq. 2) then 
   spd_plot   = gsn_csm_contour_map(wks,spd_azi,opts_uv)         ; Azimuthal average (CYL)
  else
   spd_plot   = gsn_csm_contour_map(wks,spd_plane,opts_uv)       ; Horizontal windspeed (CART)
   prs_plot   = gsn_csm_contour(wks,prs_plane,opts_prs)          ; Pressure (CART)
   if (prs .eq. 1) then 
    overlay(spd_plot, prs_plot)
   end if 
  end if 

;==================================================
; Overlay markers for each storm centre estimate
;==================================================

; (0) vort, (1) pressure, (2) pressure_smth, (3), pressure_smth2, 
; (4) vort_min, (5) windspeed_min, (6) offline
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
    markers(i)               = gsn_add_polymarker(wks,spd_plot,centre_new(i,0,it-2,1),\
                                                  centre_new(i,0,it-2,0),mres0)
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

  annotate_1                         = gsn_add_annotation(spd_plot,lgd,am_opts)

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
  plot_outer = gsn_add_polyline(wks, spd_plot, xc, yc, lnres)

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
     plot_inner(count) = gsn_add_polyline(wks, spd_plot, xc, yc, lnres)

     count = count + 1 

   end do 

  end if 

; Tidy up 
  delete([/lnres, degrees, xcos, xsin, xc, yc/])

  draw(spd_plot)
  frame(wks)

  delete([/opts_prs@cnLevels/])

  delete([/vvel_plane,u_plane,v_plane,spd_plane/])
  delete([/vort_plane,prs_plane,vort_plane0,prs_plane0,prs_pl_smth,prs_pl_smth2,vort_pl_smth/])
  delete([/vort_1d0,prs_1d0,prs_1d1,prs_1d2/])
  delete([/prs_smth,prs_smth2/]) 

 delete([/vort,u,v,u0,v0,latC,lonC,lat1C,lon1C/]) ; Tidy up 
 delete([/p,w,z,latD,lonD,lat1D,lon1D/])

 end do                ; End loop over input times (do it = 0, numINPUT_d-1)

; Write out symmetry parameter values to text file 
  diri_out = "$ar/text/cp_sym_"
  sym_out  = diri_out+dat+"_"+ens0+".txt"
  asciiwrite(sym_out,sym_arr)

end 