; Script to identify vacillation cycles (VCs) in UM ensemble simulations

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_rm.ncl"

; ncl 'opt="x11"' p_cp=3 dist=3.0 r0=1.0 'ens0="em11"' \
; ts=60 tf=66 tc=0 nep_cp_ens_thta.ncl 

; 'ens0'  = ensemble simulation (em00-em11)
; 'dist'  = size of box following storm (degrees)
; 'opt'   = output file format ("pdf" or "x11")
; 'p_cp'  = pressure level for plotting PV and vertical velocity (see L106-108)
; 'r0'    = radius for final plot (degrees)
; 'ts'	  = starting time index (0-116)
; 'tf'	  = ending time index (0-116)
; 'tc'	  = overlay TC symbol on storm centre location (1) or leave blank (0)

begin

;==============================
; Prelash (reading in files)
;==============================

; Troubleshooting
  if (r0 .gt. dist) then 
   print("Choose a bigger domain ('dist') or smaller plot size ('r0')!")
   exit()
  end if  

; Find 'pb' stream data to read in and plot 
  diri	      = "/nfs/a299/TCs/neptark/jul2018"
  fili_prefix = "20160702T1200Z_NPTK_4p4_L80_ra1t_"+ens0
  fili_b1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pb*.nc")
  fili_b      = diri+"/"+fili_b1

; Calculate number of times in file using built-in function ('nc_times')
  numINFO     = nc_times(fili_b)
  numINPUT    = numINFO[0]

; Calculate storm motion using built-in function ('st_rm')
; 'ens0'    = ensemble simulation
; 'diri'    = path to input files (see above)
; 'dist'    = size of box following storm (degrees)

  storm_rel = st_rm(ens0,"$ar/text/",dist)
  u_cyc     = storm_rel[0]             ; Zonal wind
  v_cyc     = storm_rel[1]             ; Meridional wind
  vel_cyc   = storm_rel[2]             ; Vector wind
  numTIMES  = storm_rel[3]             ; Number of times in file
  lat_arr   = storm_rel[4]             ; Latitude array subset (following storm)
  lon_arr   = storm_rel[5]             ; Longitude array subset (following storm)
  centre    = storm_rel[6]             ; Storm track (position) array 

; Create arrays for data and output information
  date_str  = new(numINPUT,string)
  time_str  = new(numINPUT,string)
  time_arr  = new(numINPUT,string)
  title_arr = new(numINPUT,string)
  llbox	    = toint(dist*50)		       ; Calculate domain size (grid points)
  dsize	    = (/llbox,llbox/)		       ; Domain size (grid points)

; Extract string for each ensemble member
  str1    = str_split(fili_b(0),"_")
  ens_str = str1(5)
  delete(str1)

  print("Working on ensemble member "+ens_str)

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 10000000000
  end setvalues

;==================================
; Get the variables we will need 
;==================================

  f = addfile(fili_b,"r")
  
  do it = ts, tf		   ; Loop over times in file 

   time  = f->t(it)                ; Time in file
   times = dimsizes(time)

 ; Account for the first three missing times in the input file 
   ct	 = it + 3

 ; Get domain boundaries at each time from storm track (calculated off-line)
   t0    = lat_arr(ct,0)
   t1    = lat_arr(ct,1)
   n0    = lon_arr(ct,0)
   n1    = lon_arr(ct,1)

 ; Create latitude/longitude array subsets to reduce computation time
   lon1  = f->longitude_1({n0:n1})  ; longitude [160 grid points]
   lat1  = f->latitude_1({t0:t1})   ; latitude  [160 grid points]
   lon   = f->longitude({n0:n1})    ; longitude [160 grid points]
   lat   = f->latitude({t0:t1})     ; latitude  [160 grid points]

 ; Pressure levels (pres)
 ; (0) 1000, (1) 950, (2) 925, (3) 900, (4) 850, (5) 800, (6) 750, (7) 700
 ; (8) 650, (9) 600, (10) 500, (11) 400, (12) 300, (13) 250
 ; (14) 200, (15) 150, (16) 100

   pres  = f->p(p_cp)               ; Subset of pressure level array
   plevs = dimsizes(pres)           ; Size of pressure level array
   pres@units = "hPa"               ; Pressure level units (for plotting)

 ; Arrays have () vertical levels, and horizontal extent depends on 'dist' option
   u0     = f->u(it,p_cp,{t0:t1},{n0:n1})       ; Zonal wind (m/s)
   v0     = f->v(it,p_cp,{t0:t1},{n0:n1})       ; Meridional wind (m/s)
   vort0  = f->rvor(it,p_cp,{t0:t1},{n0:n1})    ; Relative vorticity (/s)

   t00    = f->temp(it,p_cp,{t0:t1},{n0:n1})    ; Temperature (K)
   z0     = f->ht(it,p_cp,{t0:t1},{n0:n1})      ; Geopotential height (m)
   w0     = f->dz_dt(it,p_cp,{t0:t1},{n0:n1})   ; Vertical velocity (m/s)
   rh0    = f->rh(it,p_cp,{t0:t1},{n0:n1})      ; Relative humidity (%)

 ; Avoid dimension size errors throughout the script 
   usize  = dimsizes(u0)
   vsize  = dimsizes(v0)
   vrsize = dimsizes(vort0)
   tsize  = dimsizes(t00)
   zsize  = dimsizes(z0)
   wsize  = dimsizes(w0)
   rhsize = dimsizes(rh0)

   if ( any (usize .gt. dsize(0)) ) then 
    u = u0(0:dsize(0)-1,0:dsize(0)-1)   
   else 
    u = u0
   end if 

   if ( any (vsize .gt. dsize(0)) ) then
    v = v0(0:dsize(0)-1,0:dsize(0)-1)
   else	
    v =	v0
   end if

   if ( any (vrsize .gt. dsize(0)) ) then
    vort = vort0(0:dsize(0)-1,0:dsize(0)-1)
   else	 
    vort = vort0
   end if

   if ( any (tsize .gt. dsize(0)) ) then
    t = t00(0:dsize(0)-1,0:dsize(0)-1)
   else 
    t = t00
   end if

   if ( any (zsize .gt. dsize(0)) ) then
    z = z0(0:dsize(0)-1,0:dsize(0)-1)
   else
    z = z0
   end if

   if ( any (wsize .gt. dsize(0)) ) then
    w = w0(0:dsize(0)-1,0:dsize(0)-1)
   else 
    w = w0
   end if

   if ( any (rhsize .gt. dsize(0)) ) then
    rh   = rh0(0:dsize(0)-1,0:dsize(0)-1)
   else
    rh   = rh0
   end if

 ; Convert to units we want
   z     = z/10             ; Convert to decametres
   z@units = "dam"          ; Add metadata
   vort  = vort * (10 ^ -6) ; Convert to /s

   dx    = lon(1) - lon(0) ; Grid spacing (longitude)
   dy    = lat(1) - lat(0) ; Grid spacing (latitude)

 ; Calculate vector wind (m/s) from 'u' and 'v'
   vel   = sqrt((u ^ 2) + (v ^ 2))
   copy_VarCoords(v,vel)
   copy_VarMeta(v,vel)

 ;===============================================================================
 ; Calculate Coriolis parameter and absolute vorticity from relative vorticity
 ;===============================================================================

   f1   = coriolis_param(lat)
   f0   = conform_dims(dimsizes(vort),f1,0)	; Give 'f1' the dimensions of 'vort' 
   copy_VarMeta(vort,f0)
   copy_VarCoords(vort,f0)

   avo  = vort + f0
   copy_VarMeta(vort,avo)
   copy_VarAtts(vort,avo)
   copy_VarCoords(vort,avo)
   avo@description = "Absolute vorticity"
   avo@long_name   = "Absolute vorticity"

 ;=========================================
 ; Calculate theta-e using NCL functions
 ;=========================================

 ; First calculate mixing ratio from pressure, temperature and RH
 ; Option '1' returns the mixing ratio; '2' returns the specific humidity

 ; Create array with same dimensions as 't' for mixing ratio calculation
   tsize  = dimsizes(t)
   pres0  = new((/tsize(0),tsize(1)/),float)
   pres0  = pres
   delete(tsize)

 ; Now calculate the mixing ratio from pressure, temperature and RH
   mixr   = mixhum_ptrh(pres0,t,rh,1)

 ; Then calculate the lifting condensation level temperature using RH
 ; Option '(/1,1/)' specifies units of input and output temperatures (K)
   tlcl   = tlcl_rh_bolton(t,rh,(/1,1/))

 ; Finally calculate equivalent potential temperature
 ; Option '(/0,1,0,1/)' -- input p (hPa), temp (K), mixing ratio (kg/kg), theta-e (K)
   th_e   = pot_temp_equiv_tlcl(pres0,t,tlcl,mixr,(/0,1,0,1/))

;=====================================================
; Create correct date strings for each output time
;=====================================================

    month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
                      "Jul","Aug","Sep","Oct","Nov","Dec"/)

    times = dimsizes(time) ; Files are not all same size
    utc_date = cd_calendar(time, 0)

    year   = tointeger(utc_date(:,0))
    month  = tointeger(utc_date(:,1))
    day    = tointeger(utc_date(:,2))
    hour   = tointeger(utc_date(:,3))
    minute = tointeger(utc_date(:,4))
    second = utc_date(:,5)

  ; Correct for errors in the code (round up value of hour when minutes = 59)
    if (.not.ismissing(minute).and.minute.gt.30) then
      hour = hour+1
    end if

    date_str(it) = sprinti("%0.2iUTC ", hour) + \
                   sprinti("%0.2i ", day) \
                   + month_abbr(month)
    time_str(it) = sprinti("%0.2i UTC ", hour) + \
                   sprinti("%0.2i ", day) \
                   + month_abbr(month)
    time_arr(it)  = sprinti("%0.2i", day) + \
                   month_abbr(month) + \
                   "_" + sprinti("%0.2iZ", hour)

    print_clock("Working on time: "+time_str(it))
    print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")
    itt = ct + 1
    title_arr(it) = "valid at "+time_str(it)+" (T+"+itt+")"

;==========================================================
; Now read in variables of interest and calculate maxima
;==========================================================

     u_plane	= u
     v_plane	= v

   ; Calculate max, min and mean geopotential height (for contours)
     zmax0       = toint(max(z))
     zmin0       = toint(min(z))
     zavg0       = ( (zmax0 + zmin0) / 2)

   ; Also calculate same quantities for theta-e (uncomment to use these lines)
     tmax0       = toint(max(th_e))
     tmin0       = toint(min(th_e))
     tavg0       = ( (tmax0 + tmin0) / 2)

     panel       = new(1,graphic)       ; Panel plot with 1 image

;=================================
; Output file type and location
;=================================

     output = "$sam/nepartak/images/vc/"+ens0+"_the"+pres+"_"+time_arr(it)
     wks = gsn_open_wks(opt,output)

   ; Load colour table
     gsn_define_colormap(wks,"radar_new")

;==========================
; Options for plotting
;==========================

   ; Equivalent potential temperature
     opts_th                               = True
     opts_th@cnFillOn                      = True
     opts_th@cnLineLabelInterval           = 2.0
     opts_th@cnLineLabelFontHeightF        = 0.012
     opts_th@cnLineLabelBackgroundColor    = "transparent"
     opts_th@cnLineLabelPlacementMode      = "constant"
     opts_th@cnLinesOn                     = False ; Contour lines off
     opts_th@cnInfoLabelOn                 = False ; Contour labels off
     opts_th@cnLevelSelectionMode          = "ExplicitLevels"

   ; Pressure levels (pres)
   ; (0) 1000, (1) 950, (2) 925, (3) 900, (4) 850, (5) 800, (6) 750, (7) 700
   ; (8) 650, (9) 600, (10) 500, (11) 400, (12) 300, (13) 250
   ; (14) 200, (15) 150, (16) 100

   ; Set contour levels depending on pressure level
     if (p_cp .le. 3) then 
      opts_th@cnLevels                      = ispan(302,315,1)
      opts_th@cnFillColors                  = (/2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/)
     else if (p_cp .gt. 3 .and. p_cp .le. 5) then
      opts_th@cnLevels                      = ispan(305,318,1)
      opts_th@cnFillColors                  = (/2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/)     
     else if (p_cp .eq. 6) then 
      opts_th@cnLevels                      = ispan(310,323,1)
      opts_th@cnFillColors                  = (/2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/)
     else if (p_cp .ge. 7 .and. p_cp .le. 8) then 
      opts_th@cnLevels                      = ispan(315,328,1)
      opts_th@cnFillColors                  = (/2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/)
     else if (p_cp .eq. 9) then
      opts_th@cnLevels                      = ispan(322,335,1)
      opts_th@cnFillColors                  = (/2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/)
     else if (p_cp .eq. 10) then
      opts_th@cnLevels                      = ispan(327,340,1)
      opts_th@cnFillColors                  = (/2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/)
     end if 
     end if 
     end if  
     end if 
     end if 
     end if 

     opts_th@gsnPaperOrientation           = "landscape"
     opts_th@tiMainString                  = ""
     opts_th@tiMainFontHeightF             = 0.0125
     opts_th@gsnLeftString                 = ""
     opts_th@gsnRightString                = ""
     opts_th@lbLabelBarOn                  = True       ; Individual lb on/off

     opts_th@pmLabelBarWidthF              = 0.52
     opts_th@pmLabelBarHeightF             = 0.13
     opts_th@lbLabelFontHeightF            = 0.0125
     opts_th@lbPerimOn                     = False
     opts_th@lbLabelFont                   = "Helvetica"

     opts_th@lbTitleString                 = pres+" hPa theta-e (K)"

     opts_th@lbBoxEndCapStyle              = "TriangleBothEnds"    ; Labelbar end shape
     opts_th@lbTitleFontHeightF            = 0.02
     opts_th@lbTitleDirection              = "Across"
     opts_th@lbTitlePosition               = "Bottom"

   ; Additional plotting resources
     opts_th@mpDataBaseVersion             = "Ncarg4_1" ; More recent database
     opts_th@mpDataSetName                 = "Earth..4" ; High resolution
     opts_th@mpOutlineBoundarySets         = "National" ; National borders
     opts_th@mpGeophysicalLineColor        = "black"    ; Colour borders black
     opts_th@mpGeophysicalLineThicknessF   = 1.0        ; Border line thickness
     opts_th@mpGridAndLimbOn               = False      ; Turn on lat/lon lines
     opts_th@pmTickMarkDisplayMode         = "Always"   ; Turn on map tickmarks
     opts_th@tmXBMajorLengthF              = 0.005      ; Change tickmark length
     opts_th@tmXTOn                        = "False"    ; No tickmarks on top x-axis
     opts_th@tmYROn                        = "False"    ; No tickmarks on right y-axis
     opts_th@gsnMaximize                   = True       ; Maximise plot size

     r                                     = r0         ; Radius of plot (user input)

     opts_th@gsnAddCyclic                  = False
     opts_th@mpLimitMode                   = "Corners"
     opts_th@mpLeftCornerLatF              = centre(ct,0)-r
     opts_th@mpLeftCornerLonF              = centre(ct,1)-r
     opts_th@mpRightCornerLatF             = centre(ct,0)+r
     opts_th@mpRightCornerLonF             = centre(ct,1)+r
     opts_th@gsnDraw                       = False      ; Do not draw the plot
     opts_th@gsnFrame                      = False      ; Do not advance the frame

   ; Vertical velocity
     opts_vv                               = True
     opts_vv@cnFillOn                      = False
     opts_vv@cnLineColor                   = "white"
     opts_vv@cnInfoLabelOn                 = False
     opts_vv@cnLineLabelsOn                = False
     opts_vv@cnLevelSelectionMode          = "ExplicitLevels"
     opts_vv@cnLevels                      = (/1.0/)
     opts_vv@gsnContourLineThicknessesScale = 4.0
     opts_vv@tiMainString                  = ""
     opts_vv@gsnLeftString                 = ""
     opts_vv@gsnRightString                = ""
     opts_vv@gsnDraw                       = False      ; Do not draw the plot
     opts_vv@gsnFrame                      = False      ; Do not advance the frame

   ; Horizontal winds 
     opts_vec				   = True
     opts_vec@gsnDraw                      = False
     opts_vec@gsnFrame                     = False
     opts_vec@gsnAddCyclic                 = False

     opts_vec@vcGlyphStyle                 = "LineArrow"   ; draw windbarbs
     opts_vec@vcRefMagnitudeF              = 15.0          ; ref vec = max vector
     opts_vec@vcRefLengthF                 = 0.04          ; length of reference vector
     opts_vec@vcMinDistanceF               = 0.08          ; distance between wind barbs
     opts_vec@vcLineArrowThicknessF        = 2.5           ; windbarb line thickness
     opts_vec@vcLineArrowColor             = "black"       ; windbarb colour
     opts_vec@gsnMaximize                  = True          ; Maximise plot size
     opts_vec@gsnRightString               = ""
     opts_vec@gsnLeftString                = ""
     opts_vec@tiXAxisString                = ""
     opts_vec@vcRefAnnoOrthogonalPosF      = 0.02          ; Y-location of ref vec


;==============================================
; Plot theta-e and overlay vertical velocity
;==============================================

     plot_the  = gsn_csm_contour_map(wks,th_e,opts_th)   ; Equivalent potential temp
     plot_vv   = gsn_csm_contour(wks,w,opts_vv)    	 ; Vertical velocity
     plot_vec  = gsn_csm_vector(wks,u,v,opts_vec)	 ; Horizontal winds

     overlay(plot_the, plot_vv)
     overlay(plot_the, plot_vec)

     panel(0)  = plot_the

;==================================================
; Add NH tropical cyclone symbol to panel plots
;==================================================

    markers    = new((/dimsizes(panel)/),graphic)

  ; Add NH tropical cyclone symbol to index
    tc0 = NhlNewMarker(wks, "p", 137, 0.0, 0.0, 1.0, 1.0, 0.0)

    mres0               = True
    mres0@gsMarkerIndex = tc0         ; TC symbol
    mres0@gsMarkerSizeF = 20.0        ; Marker size
    mres0@gsMarkerColor = "black"     ; Marker colour
    mres0@gsMarkerThicknessF = 5.0    ; Marker thickness (increase)

  ; Add annotation if flag is set to '1'
    if (tc .eq. 1) then
     do m = 0, dimsizes(panel)-1
      markers(m) = gsn_add_polymarker(wks,panel(m),centre(ct,1),centre(ct,0),mres0)
     end do      
    end if 

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    optsP                       = True
    optsP@gsnFrame              = False    ; Do not advance the frame
    optsP@gsnPanelLabelBar      = False    ; Turn off panel labelbar
    optsP@txString              = "Nepartak (2016): Ensemble simulation "+ens0+", "\
    				  +title_arr(it)
    optsP@gsnPanelFigureStrings = ""
    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.90
    optsP@gsnPanelBottom        = 0.10
    optsP@amJust                = "TopLeft"
    optsP@gsnPanelFigureStringsFontHeightF = 0.0125 ; Label size (default 0.01)

    gsn_panel(wks,panel,(/1,1/),optsP)              ; Draw as a single plot
    frame(wks)

  ; Tidy up (avoid dimension size errors)
    delete([/lat1,lon1,lat,lon,time,times,w,u,v,z,t,vel,th_e,rh/])
    delete([/w0,u0,v0,z0,t00,vort0,rh0/])
    delete([/f0,f1,avo,vort/])
    delete([/utc_date,year,month,day,hour,minute,second/])

  end do      ; End time loop (do it = 0, numINPUT-1)

end 