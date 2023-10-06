; Script to plot MSLP, 10-m winds and precipitation following a TC 

; [ncl 'ens0="em02"' 'ens1="em08"' int=1 'opt="pdf"' clev=3 tc_ens_vv850_panel2.ncl]

; Where 'int' is the interval in the time loop, 'dist' is the threshold for the cyclone
; tracking code, and 'opt' is the output file format

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$LIB/gsn_csm.ncl"

begin

;======================================================================
; Read in PV and vertical velocity; start loop over ensemble members 
;======================================================================

; Prelash (define arrays for loops below)
  ens_arr     = (/ens0,ens1/)
  numSIMS     = dimsizes(ens_arr)
  numTIMES    = 120
  time_arr    = new(numTIMES,string)
  title_arr   = new(numTIMES,string)
  vv_arr      = new((/numSIMS,numTIMES,602,926/),"float") ; Vertical velocity array
  vvave_arr   = new((/numSIMS,numTIMES,602,926/),"float") ; Layer-averaged vvel array
  centre      = new((/numSIMS,numTIMES,2/),"float")       ; TC track in each simulation

  print_clock("Working on vertical velocity calculations!")

;====================================
; Start loop over ensemble members
;====================================

  do en = 0, dimsizes(ens_arr)-1

; List all files to read in and analyse
  diri        = "$sam/um/cp/ens/"
  fili_pr     = diri+"20160704T0000Z_ra1t_"+ens_arr(en)
  fili_b      = systemfunc("cd "+diri+" ; ls "+fili_pr+"*pb*.nc")
  numFILES_b  = dimsizes(fili_b)    ; Number of files (simulations)
  ct	      = 0		    ; Counter variable (time)

;===============================
; Start loop over input files
;===============================

   do nf = 0, numFILES_b-1
   f = addfile(fili_b(nf),"r")

   setvalues NhlGetWorkspaceObjectId
     "wsMaximumSize" : 1000000000
   end setvalues

;==================================
; Get the variables we will need
;==================================

   time  = f->t            ; Times in file (2 - every hour)
   times = dimsizes(time)

   lon1  = f->longitude_1  ; longitude (1098 points --> 109.04 to 152.92 degrees E)
   lat1  = f->latitude_1   ; latitude  (810 points --> 1.8 to 34.16 degrees N)
   lon   = f->longitude    ; longitude (1098 points --> 109.02 to 152.90 degrees E)
   lat   = f->latitude     ; latitude (811 points --> 1.78 to 34.18 degrees N)

   pres  = f->p            ; Pressure (levels)
   plevs = dimsizes(pres)  ; Size of pressure level array
   pres@units = "hPa"

 ; Grid subset to include all entire cylcone track
   lat_0 = 6.98
   lat_1 = 31.04
   lon_0 = 110.02
   lon_1 = 147.04

 ; Pressure levels (pres)
 ; (0) 1000, (1) 950, (2) 900, (3) 850, (4) 800, (5) 750, (6) 700, (7) 650
 ; (8) 600, (9) 550, (10) 500, (11) 450, (12) 400, (13) 350
 ; (14) 300, (15) 250, (16) 200, (17) 150, (18) 100

  ; Read in potential vorticity (PV)
    pv 	 = f->pv(:,:,{lat_0:lat_1},{lon_0:lon_1}) ; Read in PV 
    pv 	 = pv * (10 ^ 6)   			  ; Convert to PVU [10 ^ 6 K m2 s-1 kg-1]

  ; Read in vertical velocity (m/s) and concentrate on 700-300 hPa over a subset
    w  = f->dz_dt(:,:,{lat_0:lat_1},{lon_0:lon_1})
    wa = f->dz_dt(:,6:14,{lat_0:lat_1},{lon_0:lon_1}) ; Dimensions [5 * 9 * 602 * 926]

  ; Read in horizontal velocity, relative vorticity and geopotential height
    u     = f->u(:,:,{lat_0:lat_1},{lon_0:lon_1})     ; Zonal velocity
    v     = f->v(:,:,{lat_0:lat_1},{lon_0:lon_1})     ; Meridional velocity
    vort  = f->rvor(:,:,{lat_0:lat_1},{lon_0:lon_1})  ; Relative vorticity
    z     = f->ht(:,:,{lat_0:lat_1},{lon_0:lon_1})    ; Geopotential height

    z     = z/10
    z@units = "dam"         ; Convert to decametres
    vort    = vort * (10 ^ -6) ; Convert to /s

    dy = lat(1) - lat(0)    ; Grid spacing (latitude)
    dx = lon(1) - lon(0)    ; Grid spacing (longitude)

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
    do it = 0, times-1
      if (.not.ismissing(minute(it)).and.minute(it).gt.30) then
       hour(it) = hour(it)+1
      end if  
    end do

    date_str = new(times,string)
    time_str = new(times,string)
    out_str  = new(times,string)

;===========================
; Loop over times in file
;===========================

    do it = 0, times-1
 
      time_str(it) = sprinti("%0.2i UTC ", hour(it)) + \
                     sprinti("%0.2i ", day(it)) \
                     + month_abbr(month(it))
      out_str(it)  = sprinti("%0.2i", day(it)) + \
                     month_abbr(month(it)) + \
                     "_" + sprinti("%0.2iZ", hour(it))
      date_str(it) = sprinti("%0.2iUTC ", hour(it)) + \
                     sprinti("%0.2i ", day(it)) \
                     + month_abbr(month(it))

      print("Working on time: "+time_str(it))
      print("ct = "+ct)
      time_arr(ct)  = out_str(it)
      title_arr(ct) = "Valid at "+time_str(it)

;===========================
; Calculate layer average 
;===========================

      w_box = wa(it,:,:,:)			; Look at each time individually
      w_ave = dim_avg_n_Wrap(w_box,0)		; Average over levels (700-300 hPa)

      vvave_arr(en,ct,:,:) = w_ave		; Input values into larger array

;====================================================================
; Find cyclone centre using both geopotential height and vorticity
;====================================================================

      ilev = 1
      vort_plane = vort(it,ilev,:,:)
      geo_plane  = z(it,ilev,:,:)

      vort_max      = max(vort_plane)
      vort_smth     = smth9_Wrap(vort_plane, 0.5, 0.5, True)
      vort_max_smth = max(vort_smth)

      dims          = dimsizes(vort_plane)
      vort1d        = ndtooned(vort_plane)
      inds          = ind_resolve(maxind(vort1d),dims)
      vort1d_sm     = ndtooned(vort_smth)
      inds_sm       = ind_resolve(maxind(vort1d_sm),dims)

      lat_max950 = lat_0 + (dy * inds(0,0)) ; Latitude of max. vorticity
      lon_max950 = lon_0 + (dx * inds(0,1)) ; Longitude of max. vorticity
      print("Cyclone centre (vort): "+lat_max950+" degrees N, "+lon_max950+" degrees E")

    ; GEOPOTENTIAL HEIGHT
      geo_min        = min(geo_plane)
      geo_smth       = smth9_Wrap(geo_plane, 0.5, 0.5, True)
      geo_min_smth   = min(geo_smth)

      dims_h         = dimsizes(geo_plane)
      geo1d          = ndtooned(geo_plane)
      inds_h         = ind_resolve(minind(geo1d),dims_h)
      geo1d_sm       = ndtooned(geo_smth)
      inds_h_sm      = ind_resolve(minind(geo1d_sm),dims_h)

      lat_min950 = lat_0 + (dy * inds_h(0,0))
      lon_min950 = lon_0 + (dx * inds_h(0,1))
      print("Cyclone centre (hgt): "+lat_min950+" degrees N, "+lon_min950+" degrees E")

   ;========================================================
   ; EMPLOY A SAFETY NET IF WE IDENTIFY THE WRONG CYCLONE
   ;========================================================
     if (ct.eq.0) then
      centre(en,ct,0) = lat_min950
      centre(en,ct,1) = lon_min950

     else
      lt0   = centre(en,ct-1,0)         ; Previous TC latitude
      ln0   = centre(en,ct-1,1)         ; Previous TC longitude
      print("lt0 = "+lt0+" ; ln0 = "+ln0)
      d_lat = abs(lat_min950-lt0)       ; Change in latitude (t1 - t0)
      d_lon = abs(lon_min950-ln0)       ; Change in longitude (note extra minus sign)
      print("d_lat = "+d_lat+" ; d_lon = "+d_lon)

      if (d_lat.gt.1.or.d_lon.gt.1) then        ; Recalculate TC centre if incorrect

     ; Create smaller grid [0.5 degrees]
       inc = 0.5
       lt1 = lt0 - inc
       lt2 = lt0 + inc
       ln1 = ln0 - inc
       ln2 = ln0 + inc

       geo_plane0 = z(it,ilev,{lt1:lt2},{ln1:ln2})       ; Smaller grid (1.0 x 1.0 deg)
       geo_min0   = min(geo_plane0)                      ; Find minimum on smaller grid
       dims_h0    = dimsizes(geo_plane0)                 ; Size of smaller grid
       geo_1d0    = ndtooned(geo_plane0)                 ; Create 1-D array
       inds_h0    = ind_resolve(minind(geo_1d0),dims_h0) ; Find index of minimum

       lat_min950 = lt1 + (dy * inds_h0(0,0))
       lon_min950 = ln1 + (dy * inds_h0(0,1))
       print("New centre (hgt): "+lat_min950+" degrees N, "+lon_min950+" degrees E")
       centre(en,ct,0) = lat_min950
       centre(en,ct,1) = lon_min950

       delete([/geo_plane0, geo_min0, dims_h0, geo_1d0, inds_h0/])

      else
       print("Cyclone centre (hgt): "+lat_min950+" degrees N, "+lon_min950+" degrees E")
       centre(en,ct,0) = lat_min950
       centre(en,ct,1) = lon_min950
      end if

     end if

   ; Tidy up before looping over pressure levels
     delete([/vort1d, inds, vort1d_sm, vort_smth, inds_sm/])
     delete([/vort_plane, geo_plane/])

;===================================================================
; Also read in vertical velocity [proxy for convective cells]
;===================================================================

      ilev = clev
      p = pres(ilev)
      pstr = p+"hPa"

    ; Plot single-level vertical velocity
      vv_arr(en,ct,:,:) = w(it,ilev,:,:)

      print("en = "+en)
      print("ct = "+ct)
      if(any(ismissing(w(it,ilev,:,:)))) then
       print("Missing values in 'w(it,ilev,:,:)'...")
      end if

      ct = ct + 1				; Counter variable (time)

    end do	; End time loop (do it = 0, times-1)

  ; Tidy up (avoid dimension size errors)
    delete([/time,times,w,wa,pv,u,v,vort,z/])
    delete([/utc_date,year,month,day,hour,minute,second/])
    delete([/date_str,time_str,out_str/])
 
   end do	; End file loop (do nf = 0, numFILES_b-1)

  end do	; End ensemble member loop (do en = 0, dimsizes(ens_arr)-1)

;===============================================
; Produce panel plots for each time interval
;===============================================

   print_clock("Working on panel plots!")
   times = 120   		  

   do it = 23, times-1, int

   ; Read 850 hPa vertical velocity values
     vv_plane00 = vv_arr(0,it,:,:)
     vv_plane01 = vv_arr(1,it,:,:)

     panel = new(2,graphic)	; Panel plot with 6 images

;=================================
; Output file type and location
;=================================

     output = "$sam/nepartak/images/cart_coords/panel_vvel"+pstr+"_paper_"+time_arr(it)
     wks = gsn_open_wks(opt,output)

;==============================================================
; Load NIMROD colour table for precip duration data
;==============================================================

   ; Load colour table
     gsn_define_colormap(wks,"radar_new")

   ; Change the colour map around & make first colour transparent
     cmap = read_colormap_file("radar_new")
     cmap(0,3) = 0.0

;==========================
; Options for plotting
;==========================

   ; Vertical velocity on single pressure level
     opts_em00                              = True
     opts_em00@cnFillOn                     = True
     opts_em00@cnLineLabelInterval          = 2.0
     opts_em00@cnLineLabelFontHeightF       = 0.012
     opts_em00@cnLineLabelBackgroundColor   = "transparent"
     opts_em00@cnLineLabelPlacementMode     = "constant"
     opts_em00@cnLinesOn                    = False ; Contour lines off
     opts_em00@cnInfoLabelOn                = False ; Contour labels off
     opts_em00@cnLevelSelectionMode         = "ExplicitLevels"
     opts_em00@cnLevels                     = (/0.5,0.6,0.8,1.0,1.2,1.4,\
     					        1.6,1.8,2.0,2.2,2.5,3.0/)
     opts_em00@cnFillColors                 = (/0,7,8,9,10,11,12,13,14,15,16,17,18/)

     opts_em00@gsnPaperOrientation          = "landscape"
     opts_em00@tiMainString                 = ""
     opts_em00@tiMainFontHeightF            = 0.0125
     opts_em00@gsnLeftString                = ""
     opts_em00@gsnRightString               = ""
     opts_em00@lbLabelBarOn		    = False	; Turn off individual labelbars

   ; Additional plotting resources
     opts_em00@mpDataBaseVersion            = "Ncarg4_1" ; More recent database
     opts_em00@mpDataSetName                = "Earth..4" ; High resolution
     opts_em00@mpOutlineBoundarySets        = "National" ; National borders
     opts_em00@mpGeophysicalLineColor       = "black"    ; Colour borders black
     opts_em00@mpGeophysicalLineThicknessF  = 1.0        ; Border line thickness
     opts_em00@mpGridAndLimbOn              = False      ; Turn on lat/lon lines
     opts_em00@pmTickMarkDisplayMode        = "Always"   ; Turn on map tickmarks
     opts_em00@tmXBMajorLengthF             = 0.005      ; Change tickmark length
     opts_em00@tmXBLabelFontHeightF         = 0.021
     opts_em00@tmYLLabelFontHeightF         = 0.021
     opts_em00@tmXTOn                       = "False"    ; No tickmarks on top x-axis
     opts_em00@tmYROn                       = "False"    ; No tickmarks on right y-axis
     opts_em00@gsnMaximize                  = True       ; Maximise plot size 

     r					    = 0.75	 ; Radius of plot

     opts_em00@gsnAddCyclic                 = False
     opts_em00@mpLimitMode		    = "Corners"
     opts_em00@mpLeftCornerLatF             = centre(0,it,0)-r
     opts_em00@mpLeftCornerLonF             = centre(0,it,1)-r
     opts_em00@mpRightCornerLatF            = centre(0,it,0)+r
     opts_em00@mpRightCornerLonF            = centre(0,it,1)+r
     opts_em00@gsnDraw		            = False      ; Do not draw the plot
     opts_em00@gsnFrame		            = False      ; Do not advance the frame

     opts_em01                              = opts_em00
     opts_em01@mpLeftCornerLatF             = centre(1,it,0)-r
     opts_em01@mpLeftCornerLonF             = centre(1,it,1)-r
     opts_em01@mpRightCornerLatF            = centre(1,it,0)+r
     opts_em01@mpRightCornerLonF            = centre(1,it,1)+r

;=============
; Plot data 
;=============

     plot_vv0     = gsn_csm_contour_map(wks,vv_plane00,opts_em00)      ; Panel 1
     plot_vv1     = gsn_csm_contour_map(wks,vv_plane01,opts_em01)      ; Panel 2

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    panel(0)   = plot_vv0
    panel(1)   = plot_vv1

    optsP                       = True
    optsP@gsnFrame              = False    ; Do not advance the frame
    optsP@gsnPanelLabelBar      = True     ; Turn off panel labelbar
    optsP@pmLabelBarWidthF      = 0.6
    optsP@pmLabelBarHeightF     = 0.15
    optsP@lbLabelFontHeightF    = 0.01
    optsP@lbPerimOn             = False
    optsP@lbLabelFont           = "Helvetica"
    optsP@lbTitleString         = p+" hPa vertical velocity (m s~S~-1~N~)"
    optsP@lbTitleFontHeightF    = 0.015
    optsP@lbTitleDirection      = "Across"
    optsP@lbTitlePosition       = "Bottom"
    optsP@txString              = p+" hPa vertical velocity: "+title_arr(it)
    optsP@gsnPanelFigureStrings = (/"a) "+ens0, "b) "+ens1/)
    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.98
    optsP@gsnPanelBottom        = 0.02
    optsP@amJust                = "TopLeft"
    optsP@gsnPanelFigureStringsFontHeightF = 0.0125 ; Label size (default 0.01)

    gsn_panel(wks,panel,(/1,2/),optsP)            ; Draw as a single plot
    frame(wks)

   end do     ; End time loop

end 