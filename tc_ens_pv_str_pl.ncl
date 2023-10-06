; Script to identify vacillation cycles (VCs) in UM ensemble simulations

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_rm.ncl"

; ncl 'dat="03T00"' 'opt="x11"' p_cp=3 dist=3.0 r0=1.0 'ens0="em00"' \
; ts=60 tf=66 var1=1 vc=0 nep_cp_ens_vc.ncl 

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'ens0'  = ensemble simulation (em00-em11)
; 'dist'  = size of box following storm (degrees)
; 'opt'   = output file format ("pdf" or "x11")
; 'p_cp'  = pressure level for plotting PV and vertical velocity (see L106-108)
; 'r0'    = radius for final plot (degrees)
; 'var1'  = plot either theta-w (0) or theta-e (1)
; 'ts'	  = starting time index (0-116)
; 'tf'	  = ending time index (0-116)
; 'vc'	  = plot 1-h (0) or 5-min (1) data

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
  if (vc .eq. 0) then 
   diri        = "/nfs/a319/earshar/"+dat+"/"+ens0
   fili_prefix = "20160702T1200Z_NPTK_4p4_L80_ra1t_"+ens0
   fili_b1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pb*.nc")
   fili_b      = diri+"/"+fili_b1
  elseif (vc .eq. 1) then 
   diri        = "/nfs/a319/earshar/02T12/em11/vc
   fili_prefix = "20160702T1200Z_NPTK_4p4_L80_ra1t_em11"
   fili_b1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pb*.nc")
   fili_b      = diri+"/"+fili_b1
  end if 

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
   u0     = f->u(it,p_cp,{t0:t1},{n0:n1})       ; Zonal velocity (m/s)
   v0     = f->v(it,p_cp,{t0:t1},{n0:n1})       ; Meridional velocity (m/s)
   vort0  = f->rvor(it,p_cp,{t0:t1},{n0:n1})    ; Relative vorticity (/s)
   div0   = f->div(it,p_cp,{t0:t1},{n0:n1})     ; Divergence (/s)
   pv0    = f->pv(it,p_cp,{t0:t1},{n0:n1})      ; PV (PVU)

   t00    = f->temp(it,p_cp,{t0:t1},{n0:n1})    ; Temperature (K)
   z0     = f->ht(it,p_cp,{t0:t1},{n0:n1})      ; Geopotential height (m)
   w0     = f->dz_dt(it,p_cp,{t0:t1},{n0:n1})   ; Vertical velocity (m/s)
   thw0	  = f->field22(it,p_cp,{t0:t1},{n0:n1}) ; Wet bulb potential temperature (K)
   rh0    = f->rh(it,p_cp,{t0:t1},{n0:n1})      ; Relative humidity (%)

   vel0   = sqrt((u0^2) + (v0^2))             ; Vector velocity (m/s)
   copy_VarCoords(v0,vel0)
   copy_VarMeta(v0,vel0)

 ; Avoid dimension size errors throughout the script 
   usize  = dimsizes(u0)
   vsize  = dimsizes(v0)
   vrsize = dimsizes(vort0)
   dvsize = dimsizes(div0)
   pvsize = dimsizes(pv0)
   tsize  = dimsizes(t00)
   zsize  = dimsizes(z0)
   wsize  = dimsizes(w0)
   vcsize = dimsizes(vel0)
   thsize = dimsizes(thw0)
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

   if ( any (dvsize .gt. dsize(0)) ) then
    div = div0(0:dsize(0)-1,0:dsize(0)-1)
   else 
    div = div0
   end if

   if ( any (pvsize .gt. dsize(0)) ) then
    pv = pv0(0:dsize(0)-1,0:dsize(0)-1)
   else 
    pv = pv0
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

   if ( any (vcsize .gt. dsize(0)) ) then
    vel = vel0(0:dsize(0)-1,0:dsize(0)-1)
   else 
    vel = vel0
   end if

   if ( any (thsize .gt. dsize(0)) ) then
    th_w = thw0(0:dsize(0)-1,0:dsize(0)-1)
   else
    th_w = thw0
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
   div   = div  * (10 ^ -6) ; Convert to /s
   pv    = pv * (10 ^ 6)    ; Convert to PVU [10 ^ 6 K m2 s-1 kg-1]

   dx    = lon(1) - lon(0) ; Grid spacing (longitude)
   dy    = lat(1) - lat(0) ; Grid spacing (latitude)

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
     str_plane0  = sqrt( str1(:,:) + str2(:,:) )

;===================================================================
; Also read in lower-tropospheric PV [proxy for convective cells]
;===================================================================

   ; Read in PV, geopotential height and vertical velocity
     vort_plane = vort
     pv_plane   = pv
     geo_plane  = z
     w_plane    = w
     thw_plane	= th_w
     the_plane	= th_e

     str_plane  = str_plane0
     copy_VarCoords(u,str_plane)

   ; Calculate max, min and mean geopotential height (for contours)
     zmax0       = toint(max(geo_plane))
     zmin0       = toint(min(geo_plane))
     zavg0       = ( (zmax0 + zmin0) / 2)

   ; Also calculate same quantities for theta-e (uncomment to use these lines)
     tmax0       = toint(max(the_plane))
     tmin0       = toint(min(the_plane))
     tavg0       = ( (tmax0 + tmin0) / 2)

     panel       = new(3,graphic)       ; Panel plot with xxx images

;=================================
; Output file type and location
;=================================

     output = "$sam/nepartak/images/vc/"+ens0+"_"+dat+"_str_pv_"+pres+"_"+time_arr(it)
     wks = gsn_open_wks(opt,output)

   ; Load colour table
     gsn_define_colormap(wks,"prcp_new")

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
     opts_pv@cnLevels                     = (/9.0, 10.0 ,12.0, 14.0, \
                                              16.0, 18.0, 20.0, 22.0, 25.0, \
                                              30.0, 35.0, 40.0, 45.0, 50.0/)
     opts_pv@cnFillColors                 = (/0,2,3,4,5,\
                                              6,7,8,9,10,\
                                              11,12,13,14,15/)
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
     opts_pv@lbTitleString                 = pres+" hPa potential vorticity (PVU)"
     opts_pv@lbBoxEndCapStyle              = "TriangleBothEnds"    ; Labelbar end shape
     opts_pv@lbTitleFontHeightF            = 0.02
     opts_pv@lbTitleDirection              = "Across"
     opts_pv@lbTitlePosition               = "Bottom"

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
     opts_pv@gsnMaximize                   = True       ; Maximise plot size

     r                                     = r0         ; Radius of plot (user input)

     opts_pv@gsnAddCyclic                  = False
     opts_pv@mpLimitMode                   = "Corners"
     opts_pv@mpLeftCornerLatF              = centre(ct,0)-r
     opts_pv@mpLeftCornerLonF              = centre(ct,1)-r
     opts_pv@mpRightCornerLatF             = centre(ct,0)+r
     opts_pv@mpRightCornerLonF             = centre(ct,1)+r
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
     opts_vv@cnLevels                      = (/1.0,2.0/)
     opts_vv@cnLineLabelPlacementMode      = "constant"
     opts_vv@cnLineLabelPerimOn            = False
     opts_vv@gsnContourLineThicknessesScale = 3.0
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

   ; Geopotential height
     opts_geo                              = opts_vv
     delete(opts_geo@cnLevels)
     opts_geo@cnLevels                     = ispan(toint(zmin0),toint(zmax0),4)
     opts_geo@gsnContourLineThicknessesScale = 3.0

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
     opts_str@lbTitleString = pres+" hPa strain rate (10~S~-4~N~ s~S~-2~N~)"
     opts_str@lbBoxEndCapStyle             = "TriangleBothEnds"    ; Labelbar end shape
     opts_str@lbTitleFontHeightF           = 0.02
     opts_str@lbTitleDirection             = "Across"
     opts_str@lbTitlePosition              = "Bottom"

;================================================
; Finally plot wet bulb potential temperature 
;================================================

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
;     opts_the@cnFillPalette                = cmap; (2:,:)        ; Skip first 2 colours 
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

;=============
; Plot data
;=============

     str_plane   = str_plane * (10 ^ 4)
     str_plane@units = "10~S~-4~N~ s~S~-2~N~"

     vort_plane  = vort_plane * (10 ^ 4)
     vort_plane@units = "10~S~-4~N~ s~S~-1~N~"

     plot_pv     = gsn_csm_contour_map(wks,pv_plane,opts_pv)   ; Potential vorticity

     plot_vv     = gsn_csm_contour(wks,w_plane,opts_vv)       ; Vertical velocity
     overlay(plot_pv, plot_vv)
     panel(0)   = plot_pv

     plot_str   = gsn_csm_contour_map(wks,str_plane,opts_str) ; Strain rate
     panel(1)   = plot_str

     plot_the  = gsn_csm_contour_map(wks,the_plane,opts_the) ; Equivalent potential temp
     panel(2)	= plot_the

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

    do m = 0, dimsizes(panel)-1
     markers(m) = gsn_add_polymarker(wks,panel(m),centre(ct,1),centre(ct,0),mres0)
    end do      

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    optsP                       = True
    optsP@gsnFrame              = False    ; Do not advance the frame
    optsP@gsnPanelLabelBar      = False    ; Turn off panel labelbar
;    optsP@txString              = "Nepartak (2-9 July 2016): "+title_arr(it)
    optsP@txString              = "Nepartak (2016): Ensemble simulation "+ens0+", "\
    				  +title_arr(it)
    optsP@gsnPanelFigureStrings = ""
    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.90
    optsP@gsnPanelBottom        = 0.10
    optsP@amJust                = "TopLeft"
    optsP@gsnPanelFigureStringsFontHeightF = 0.0125 ; Label size (default 0.01)

    gsn_panel(wks,panel,(/1,3/),optsP)              ; Draw as a single plot
    frame(wks)

  ; Tidy up before next iteration
    delete([/pv_plane, geo_plane, thw_plane, w_plane/])
    delete([/str_plane, str_plane0, str1, str2/])

  ; Tidy up (avoid dimension size errors)
    delete([/lat1,lon1,lat,lon,time,times,w,pv,u,v,z,t,vel,th_w,th_e,rh/])
    delete([/w0,pv0,u0,v0,z0,t00,vel0,vort0,div0,thw0,rh0/])
    delete([/f0,f1,avo,vort,div/])
    delete([/du_dx,dv_dy,dv_dx,du_dy/])
    delete([/utc_date,year,month,day,hour,minute,second/])

  end do      ; End time loop (do it = 0, numINPUT-1)

end 