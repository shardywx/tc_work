; Script to plot PV and vertical velocity in the TC inner core

; Run using e.g. [ncl dist=4.0 clev= 1 'opt="x11"' tc_ens_pv_panel2.ncl]
; Where 'int' is the interval in the time loop, 'dist' is the threshold for the cyclone
; tracking code, and 'opt' is the output file format

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"

begin

;=================================================
; Read in cyclone position data from text files
;=================================================

   diri        = "/nfs/a37/earshar/ncl_scripts/nepartak/archer_march2018/"
   numTIMES    = 120
   centre      = new((/2,numTIMES,2/),float) ; Array to hold all [120] lat/lon values

 ; Read in values from ensemble members

   centre(0,:,0)  = asciiread(diri+"cp_lat_em00.txt",(/numTIMES,1/),"float")
   centre(0,:,1)  = asciiread(diri+"cp_lon_em00.txt",(/numTIMES,1/),"float")

   centre(1,:,0)  = asciiread(diri+"cp_lat_em01.txt",(/numTIMES,1/),"float")
   centre(1,:,1)  = asciiread(diri+"cp_lon_em01.txt",(/numTIMES,1/),"float")

   centre(2,:,0)  = asciiread(diri+"cp_lat_em02.txt",(/numTIMES,1/),"float")
   centre(2,:,1)  = asciiread(diri+"cp_lon_em02.txt",(/numTIMES,1/),"float")

   centre(3,:,0)  = asciiread(diri+"cp_lat_em03.txt",(/numTIMES,1/),"float")
   centre(3,:,1)  = asciiread(diri+"cp_lon_em03.txt",(/numTIMES,1/),"float")

   centre(4,:,0)  = asciiread(diri+"cp_lat_em04.txt",(/numTIMES,1/),"float")
   centre(4,:,1)  = asciiread(diri+"cp_lon_em04.txt",(/numTIMES,1/),"float")

   centre(5,:,0)  = asciiread(diri+"cp_lat_em05.txt",(/numTIMES,1/),"float")
   centre(5,:,1)  = asciiread(diri+"cp_lon_em05.txt",(/numTIMES,1/),"float")

   centre(6,:,0)  = asciiread(diri+"cp_lat_em06.txt",(/numTIMES,1/),"float")
   centre(6,:,1)  = asciiread(diri+"cp_lon_em06.txt",(/numTIMES,1/),"float")

   centre(7,:,0)  = asciiread(diri+"cp_lat_em07.txt",(/numTIMES,1/),"float")
   centre(7,:,1)  = asciiread(diri+"cp_lon_em07.txt",(/numTIMES,1/),"float")

   centre(8,:,0)  = asciiread(diri+"cp_lat_em08.txt",(/numTIMES,1/),"float")
   centre(8,:,1)  = asciiread(diri+"cp_lon_em08.txt",(/numTIMES,1/),"float")

   centre(9,:,0)  = asciiread(diri+"cp_lat_em09.txt",(/numTIMES,1/),"float")
   centre(9,:,1)  = asciiread(diri+"cp_lon_em09.txt",(/numTIMES,1/),"float")

   centre(10,:,0) = asciiread(diri+"cp_lat_em10.txt",(/numTIMES,1/),"float")
   centre(10,:,1) = asciiread(diri+"cp_lon_em10.txt",(/numTIMES,1/),"float")

   centre(11,:,0) = asciiread(diri+"cp_lat_em11.txt",(/numTIMES,1/),"float")
   centre(11,:,1) = asciiread(diri+"cp_lon_em11.txt",(/numTIMES,1/),"float")

;================================================================================
; Create lat/lon arrays containing grid subsets at each time, following the TC
;================================================================================
  lat_arr = new((/12,numTIMES,2/),float)
  lon_arr = new((/12,numTIMES,2/),float)

; Lat/lon arrays for ensemble simulations
  lat_arr(0,:,0) = centre(0,:,0)-dist
  lat_arr(0,:,1) = centre(0,:,0)+dist
  lon_arr(0,:,0) = centre(0,:,1)-dist
  lon_arr(0,:,1) = centre(0,:,1)+dist

  lat_arr(1,:,0) = centre(1,:,0)-dist
  lat_arr(1,:,1) = centre(1,:,0)+dist
  lon_arr(1,:,0) = centre(1,:,1)-dist
  lon_arr(1,:,1) = centre(1,:,1)+dist

  lat_arr(2,:,0) = centre(2,:,0)-dist
  lat_arr(2,:,1) = centre(2,:,0)+dist
  lon_arr(2,:,0) = centre(2,:,1)-dist
  lon_arr(2,:,1) = centre(2,:,1)+dist

  lat_arr(3,:,0) = centre(3,:,0)-dist
  lat_arr(3,:,1) = centre(3,:,0)+dist
  lon_arr(3,:,0) = centre(3,:,1)-dist
  lon_arr(3,:,1) = centre(3,:,1)+dist

  lat_arr(4,:,0) = centre(4,:,0)-dist
  lat_arr(4,:,1) = centre(4,:,0)+dist
  lon_arr(4,:,0) = centre(4,:,1)-dist
  lon_arr(4,:,1) = centre(4,:,1)+dist

  lat_arr(5,:,0) = centre(5,:,0)-dist
  lat_arr(5,:,1) = centre(5,:,0)+dist
  lon_arr(5,:,0) = centre(5,:,1)-dist
  lon_arr(5,:,1) = centre(5,:,1)+dist

  lat_arr(6,:,0) = centre(6,:,0)-dist
  lat_arr(6,:,1) = centre(6,:,0)+dist
  lon_arr(6,:,0) = centre(6,:,1)-dist
  lon_arr(6,:,1) = centre(6,:,1)+dist

  lat_arr(7,:,0) = centre(7,:,0)-dist
  lat_arr(7,:,1) = centre(7,:,0)+dist
  lon_arr(7,:,0) = centre(7,:,1)-dist
  lon_arr(7,:,1) = centre(7,:,1)+dist

  lat_arr(8,:,0) = centre(8,:,0)-dist
  lat_arr(8,:,1) = centre(8,:,0)+dist
  lon_arr(8,:,0) = centre(8,:,1)-dist
  lon_arr(8,:,1) = centre(8,:,1)+dist

  lat_arr(9,:,0) = centre(9,:,0)-dist
  lat_arr(9,:,1) = centre(9,:,0)+dist
  lon_arr(9,:,0) = centre(9,:,1)-dist
  lon_arr(9,:,1) = centre(9,:,1)+dist

  lat_arr(10,:,0) = centre(10,:,0)-dist
  lat_arr(10,:,1) = centre(10,:,0)+dist
  lon_arr(10,:,0) = centre(10,:,1)-dist
  lon_arr(10,:,1) = centre(10,:,1)+dist

  lat_arr(11,:,0) = centre(11,:,0)-dist
  lat_arr(11,:,1) = centre(11,:,0)+dist
  lon_arr(11,:,0) = centre(11,:,1)-dist
  lon_arr(11,:,1) = centre(11,:,1)+dist

;==================================================================================
;;;;;;;;;;;;;;;;;;;;;;;;;;; Start of main part of script ;;;;;;;;;;;;;;;;;;;;;;;;;;
;==================================================================================

; Prelash (define arrays for loops below)
  ens_arr     = (/"em00","em01","em02","em03","em04","em05",\
  	      	  "em06","em07","em08","em09","em10","em11"/)
  numSIMS     = dimsizes(ens_arr)
  numTIMES    = 120

  time_arr    = new(numTIMES,string)			  ; Time equal for both plots
  title_arr   = new(numTIMES,string)			  ; Title equal for both plots
  vel_arr     = new((/numSIMS,numTIMES,150,150/),"float") ; Velocity array

  date_str    = new(numTIMES,string)
  time_str    = new(numTIMES,string)
  out_str     = new(numTIMES,string)

  print_clock("Working on windspeed calculations!")

;====================================
; Start loop over ensemble members
;====================================

  do en = 0, numSIMS-1

; List all files to read in and analyse
  diri        = "$ens/reg_netcdf/20160704T0000Z_ra1t_"
  input_pb    = diri+ens_arr(en)+"_pb.nc"

  f  = addfile(input_pb,"r")
  print("Working on ensemble simulation: "+ens_arr(en))

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 1000000000
  end setvalues

;==================================
; Get the variables we will need
;==================================

;  do it = 0, numTIMES-1	   ; Loop over all times in file
  do it = 47, 53	   ; Loop over subset of times

   time  = f->t(it)        ; Time in file
   times = dimsizes(time)

   t0    = lat_arr(en,it,0)
   t1    = lat_arr(en,it,1)
   n0    = lon_arr(en,it,0)
   n1    = lon_arr(en,it,1)
   print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

   lon1  = f->longitude_1({n0:n1})  ; longitude [160 grid points]
   lat1  = f->latitude_1({t0:t1})   ; latitude  [160 grid points]
   lon   = f->longitude({n0:n1})    ; longitude [160 grid points]
   lat   = f->latitude({t0:t1})     ; latitude  [160 grid points]

 ; Pressure levels (pres)
 ; (0) 1000, (1) 950, (2) 900, (3) 850, (4) 800, (5) 750, (6) 700, (7) 650
 ; (8) 600, (9) 550, (10) 500, (11) 450, (12) 400, (13) 350
 ; (14) 300, (15) 250, (16) 200, (17) 150, (18) 100

   p_cp  = (/1,3,6,10,16/)          ; Pressure levels to access (950,850,700,500,200)

   pres  = f->p(p_cp)               ; Pressure (levels)
   plevs = dimsizes(pres)           ; Size of pressure level array
   pres@units = "hPa"

 ; Arrays have 5 vertical levels, and horizontal extent depends on 'dist' option
   u     = f->u(it,p_cp,{t0:t1},{n0:n1})     ; Zonal velocity (m/s)
   v     = f->v(it,p_cp,{t0:t1},{n0:n1})     ; Meridional velocity (m/s)
   vort  = f->rvor(it,p_cp,{t0:t1},{n0:n1})  ; Relative vorticity (/s)

   t     = f->temp(it,p_cp,{t0:t1},{n0:n1})  ; Temperature (K)
   z     = f->ht(it,p_cp,{t0:t1},{n0:n1})    ; Geopotential height (m)

   vel	 = sqrt((u^2) + (v^2))		     ; Vector velocity (m/s)
   copy_VarCoords(v,vel)
   copy_VarMeta(v,vel)

 ; Convert to units we want
   z     = z/10
   z@units = "dam"	    ; Convert to decametres
   vort  = vort * (10 ^ -6) ; Convert to /s

   dx    = lon(1) - lon(0) ; Grid spacing (longitude)
   dy    = lat(1) - lat(0) ; Grid spacing (latitude)

 ; Coriolis parameter and absolute vorticity
   f1   = coriolis_param(lat)
   f0   = conform_dims(dimsizes(vort),f1,1)
   copy_VarMeta(vort,f0)
   copy_VarCoords(vort,f0)

   avo  = vort + f0
   copy_VarMeta(vort,avo)
   copy_VarAtts(vort,avo)
   copy_VarCoords(vort,avo)
   avo@description = "Absolute vorticity"
   avo@long_name   = "Vorticity (absolute)"

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
    out_str(it)  = sprinti("%0.2i", day) + \
                   month_abbr(month) + \
                   "_" + sprinti("%0.2iZ", hour)

    print_clock("Working on time: "+time_str(it))
    time_arr(it)  = out_str(it)
    title_arr(it) = "Valid at "+time_str(it)

;====================================================================
; Find cyclone centre using both geopotential height and vorticity
;====================================================================

  ; Look at 950 hPa
    ilev	    = 0			
    vort_plane      = vort(ilev,:,:)
    geo_plane	    = z(ilev,:,:)

    vort_max      = max(vort_plane)
    vort_smth     = smth9_Wrap(vort_plane, 0.5, 0.5, True)
    vort_max_smth = max(vort_smth)

    dims          = dimsizes(vort_plane)
    vort1d        = ndtooned(vort_plane)
    inds          = ind_resolve(maxind(vort1d),dims)
    vort1d_sm     = ndtooned(vort_smth)
    inds_sm       = ind_resolve(maxind(vort1d_sm),dims)

    lat_max950 = t0 + (dy * inds(0,0)) ; Latitude of max. vorticity
    lon_max950 = n0 + (dx * inds(0,1)) ; Longitude of max. vorticity
    print("Cyclone centre (pre-script): "+centre(en,it,0)+" deg N, "+centre(en,it,1)+\
          " deg E")
    print("Cyclone centre (script - VORT): "+lat_max950+" deg N, "+lon_max950+\
          " deg E")

  ; GEOPOTENTIAL HEIGHT
    geo_min        = min(geo_plane)
    geo_smth       = smth9_Wrap(geo_plane, 0.5, 0.5, True)
    geo_min_smth   = min(geo_smth)
 
    dims_h         = dimsizes(geo_plane)
    geo1d          = ndtooned(geo_plane)
    inds_h         = ind_resolve(minind(geo1d),dims_h)
    geo1d_sm       = ndtooned(geo_smth)
    inds_h_sm      = ind_resolve(minind(geo1d_sm),dims_h)

    lat_min950 = t0 + (dy * inds_h(0,0))
    lon_min950 = n0 + (dx * inds_h(0,1))
    print("Cyclone centre (script - HGT): "+lat_min950+" deg N, "+lon_min950+\
          " deg E")

    delete([/geo_smth, dims_h, geo1d, inds_h, geo1d_sm, inds_h_sm/])

   ;========================================================
   ; EMPLOY A SAFETY NET IF WE IDENTIFY THE WRONG CYCLONE
   ;========================================================
     if (it.eq.0) then
;      centre(en,it,0) = lat_min950
;      centre(en,it,1) = lon_min950
      lt0	      = centre(en,it,0)
      ln0	      = centre(en,it,1)

     else
      lt0   = centre(en,it-1,0)         ; Previous TC latitude
      ln0   = centre(en,it-1,1)         ; Previous TC longitude
      print("lt0 = "+lt0+" ; ln0 = "+ln0)
      d_lat = abs(lat_min950-lt0)       ; Change in latitude (t1 - t0)
      d_lon = abs(lon_min950-ln0)       ; Change in longitude (note extra minus sign)
      print("d_lat = "+d_lat+" ; d_lon = "+d_lon)

      if (d_lat.gt.0.75.or.d_lon.gt.0.75) then   ; Recalculate TC centre if incorrect

     ; Create smaller grid [0.5 degrees]
       inc = 0.5
       lt1 = lt0 - inc
       lt2 = lt0 + inc
       ln1 = ln0 - inc
       ln2 = ln0 + inc

       geo_plane0 = z(ilev,{lt1:lt2},{ln1:ln2})          ; Smaller grid (1.0 x 1.0 deg)
       geo_min0   = min(geo_plane0)                      ; Find minimum on smaller grid
       dims_h0    = dimsizes(geo_plane0)                 ; Size of smaller grid
       geo_1d0    = ndtooned(geo_plane0)                 ; Create 1-D array
       inds_h0    = ind_resolve(minind(geo_1d0),dims_h0) ; Find index of minimum

       lat_min950 = lt1 + (dy * inds_h0(0,0))
       lon_min950 = ln1 + (dy * inds_h0(0,1))
       print("New centre (hgt): "+lat_min950+" degrees N, "+lon_min950+" degrees E")
;       centre(en,it,0) = lat_min950
;       centre(en,it,1) = lon_min950

       delete([/geo_plane0, geo_min0, dims_h0, geo_1d0, inds_h0/])

      else
       print("Cyclone centre (hgt): "+lat_min950+" degrees N, "+lon_min950+" degrees E")
;       centre(en,it,0) = lat_min950
;       centre(en,it,1) = lon_min950
      end if

     end if

   ; Tidy up before looping over pressure levels
     delete([/vort1d, inds, dims, vort1d_sm, vort_smth, inds_sm/])
     delete([/vort_plane, geo_plane/])

;===================================================================
; Also read in lower-tropospheric PV [proxy for convective cells]
;===================================================================

      ilev = clev				 ; Choose pressure level
      p	   = pres(ilev)				 ; String for output file
      vel_arr(en,it,:,:) = vel(ilev,:,:)         ; Input values into velocity array

      print("en = "+en)
      print("it = "+it)
      if(any(ismissing(vel(ilev,:,:)))) then
       print("Missing values in PV array...")
      end if

    ; Tidy up (avoid dimension size errors)
      delete([/lat1,lon1,lat,lon,time,times,w,pv,u,v,vort,div,z,t,vel/])
      delete([/utc_date,year,month,day,hour,minute,second/])

    end do	; End time loop (do it = 0, times-1)

  end do	; End ensemble member loop (do en = 0, dimsizes(ens_arr)-1)

;===============================================
; Produce panel plots for each time interval
;===============================================

   print("Working on panel plots!")
   times = 120   		  

;   do it = 23, times-1, 1
   do it = 47, 53, 1

   print("Working on time: "+it)

   ; Read 850 hPa PV values
     pv_plane00 = pv_arr(0,it,:,:)
     pv_plane01 = pv_arr(1,it,:,:)

   ; Read vertical velocity
     vv_plane00 = vv_arr(0,it,:,:)
     vv_plane01 = vv_arr(1,it,:,:)

     panel = new(2,graphic)	; Panel plot with 6 images

;=================================
; Output file type and location
;=================================

     output = "$sam/nepartak/images/cart_coords/panel_pv"+p+"_paper_"+time_arr(it)
     wks = gsn_open_wks(opt,output)

   ; Load colour table
     gsn_define_colormap(wks,"prcp_new")

;==========================
; Options for plotting
;==========================

   ; Potential vorticity
     opts_em00                              = True
     opts_em00@cnFillOn                     = True
     opts_em00@cnLineLabelInterval          = 2.0
     opts_em00@cnLineLabelFontHeightF       = 0.012
     opts_em00@cnLineLabelBackgroundColor   = "transparent"
     opts_em00@cnLineLabelPlacementMode     = "constant"
     opts_em00@cnLinesOn                    = False ; Contour lines off
     opts_em00@cnInfoLabelOn                = False ; Contour labels off
     opts_em00@cnLevelSelectionMode         = "ExplicitLevels"

;;     opts_em00@cnLevels                     = (/9,12,15,20,25,30,40,50/)
;     opts_em00@cnLevels                     = (/15,20,25,30,35,40,45,50/)
;;     opts_em00@cnLevels                     = (/20,25,30,35,40,45,50,60/) 
;     opts_em00@cnFillColors                 = (/0,13,15,17,19,20,21,22,23/)

     opts_em00@cnLevels                     = (/9.0, 10.0 ,12.0, 14.0, \
                                                16.0, 18.0, 20.0, 22.0, 25.0, \
                                                30.0, 35.0, 40.0, 45.0, 50.0/)
     opts_em00@cnFillColors                 = (/0,2,3,4,5,\
                                                6,7,8,9,10,\
                                                11,12,13,14,15/)

     opts_em00@gsnPaperOrientation          = "landscape"
     opts_em00@tiMainString                 = ""
     opts_em00@tiMainFontHeightF            = 0.0125
     opts_em00@gsnLeftString                = ""
     opts_em00@gsnRightString               = ""
     opts_em00@lbLabelBarOn		    = False	 ; Turn off individual labelbars

   ; Additional plotting resources
     opts_em00@mpDataBaseVersion            = "Ncarg4_1" ; More recent database
     opts_em00@mpDataSetName                = "Earth..4" ; High resolution
     opts_em00@mpOutlineBoundarySets        = "National" ; National borders
     opts_em00@mpGeophysicalLineColor       = "black"    ; Colour borders black
     opts_em00@mpGeophysicalLineThicknessF  = 1.0        ; Border line thickness
     opts_em00@mpGridAndLimbOn              = False      ; Turn on lat/lon lines
     opts_em00@pmTickMarkDisplayMode        = "Always"   ; Turn on map tickmarks
     opts_em00@tmXBMajorLengthF             = 0.005      ; Change tickmark length
;     opts_em00@tmXBLabelFontHeightF         = 0.021
;     opts_em00@tmYLLabelFontHeightF         = 0.021
     opts_em00@tmXTOn                       = "False"    ; No tickmarks on top x-axis
     opts_em00@tmYROn                       = "False"    ; No tickmarks on right y-axis
     opts_em00@gsnMaximize                  = True       ; Maximise plot size 

     r					    = 1.00	 ; Radius of plot

     opts_em00@gsnAddCyclic                 = False
     opts_em00@mpLimitMode		    = "Corners"
     opts_em00@mpLeftCornerLatF             = centre(0,it,0)-r
     opts_em00@mpLeftCornerLonF             = centre(0,it,1)-r
     opts_em00@mpRightCornerLatF            = centre(0,it,0)+r
     opts_em00@mpRightCornerLonF            = centre(0,it,1)+r
     opts_em00@gsnDraw		            = False      ; Do not draw the plot
     opts_em00@gsnFrame		            = False      ; Do not advance the frame

   ; Vertical velocity
     opts_vv00                               = True
     opts_vv00@cnFillOn                      = False
     opts_vv00@cnLineColor                   = "black"
     opts_vv00@cnInfoLabelOn                 = False
     opts_vv00@cnLineLabelsOn                = True
     opts_vv00@cnLineLabelInterval           = 2.0
     opts_vv00@cnLevelSelectionMode          = "ExplicitLevels"
     opts_vv00@cnLevels                      = (/0.5,1.0,1.5/)
     opts_vv00@cnLineLabelPlacementMode      = "constant"
     opts_vv00@cnLineLabelPerimOn            = False
     opts_vv00@gsnContourLineThicknessesScale = 2.0
     opts_vv00@tiMainString                  = ""
     opts_vv00@gsnLeftString                 = ""
     opts_vv00@gsnRightString                = ""
     opts_vv00@gsnDraw                       = False      ; Do not draw the plot
     opts_vv00@gsnFrame                      = False      ; Do not advance the frame
     opts_vv00@pmTickMarkDisplayMode         = "Always"   ; Turn on map tickmarks
     opts_vv00@tmXBMajorLengthF              = 0.005      ; Change tickmark length
     opts_vv00@tmXTOn                        = "False"    ; No tickmarks on top x-axis
     opts_vv00@tmYROn                        = "False"    ; No tickmarks on right y-axis
     opts_vv00@gsnMaximize                   = True       ; Maximise plot size

     opts_em01                              = opts_em00
     opts_em01@mpLeftCornerLatF             = centre(1,it,0)-r
     opts_em01@mpLeftCornerLonF             = centre(1,it,1)-r
     opts_em01@mpRightCornerLatF            = centre(1,it,0)+r
     opts_em01@mpRightCornerLonF            = centre(1,it,1)+r

;=============
; Plot data 
;=============

     plot_pv0     = gsn_csm_contour_map(wks,pv_plane00,opts_em00)      ; Panel 1
;     plot_vv00    = gsn_csm_contour(wks,vv_plane00,opts_vv00)          ; Omega
;     overlay(plot_pv0, plot_vv00)

     plot_pv1     = gsn_csm_contour_map(wks,pv_plane01,opts_em01)      ; Panel 2
;     plot_vv01    = gsn_csm_contour(wks,vv_plane01,opts_vv00)          ; Omega
;     overlay(plot_pv1, plot_vv01)

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    panel(0)   = plot_pv0
    panel(1)   = plot_pv1

    optsP                       = True
    optsP@gsnFrame              = False    ; Do not advance the frame
    optsP@gsnPanelLabelBar      = True     ; Turn off panel labelbar
    optsP@pmLabelBarWidthF      = 0.6
    optsP@pmLabelBarHeightF     = 0.15
    optsP@lbLabelFontHeightF    = 0.01
    optsP@lbPerimOn             = False
    optsP@lbLabelFont           = "Helvetica"
    optsP@lbTitleString         = p+" hPa potential vorticity (PVU)"
    optsP@lbTitleFontHeightF    = 0.015
    optsP@lbTitleDirection      = "Across"
    optsP@lbTitlePosition       = "Bottom"
    optsP@txString              = p+" hPa PV and vertical velocity: "+title_arr(it)
    optsP@gsnPanelFigureStrings = (/"a) "+ens_arr(0), "b) "+ens_arr(1)/)
    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.98
    optsP@gsnPanelBottom        = 0.02
    optsP@amJust                = "TopLeft"
    optsP@gsnPanelFigureStringsFontHeightF = 0.0125 ; Label size (default 0.01)

    gsn_panel(wks,panel,(/1,2/),optsP)             ; Draw as a single plot
    frame(wks)

  ; Tidy up before next iteration
    delete([/pv_plane00, pv_plane01, vv_plane00, vv_plane01/])

   end do     ; End time loop

end 