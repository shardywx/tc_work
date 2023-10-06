; Script to plot MSLP, 10-m winds and precipitation following a TC 
; Run using e.g. [ncl int=1 dist=1 'opt="x11"' tc_ens_pv_panel_vp.ncl]
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
  ens_arr     = (/"em00","em01","em02","em03","em04","em05",\
                  "em06","em07","em08","em09","em10","em11"/)
  numSIMS     = dimsizes(ens_arr)
  numTIMES    = 120
  time_arr    = new(numTIMES,string)
  title_arr   = new(numTIMES,string)
  pv_arr      = new((/numSIMS,numTIMES,602,926/),"float") ; PV array
  vv_arr      = new((/numSIMS,numTIMES,602,926/),"float") ; Vertical velocity array
  centre      = new((/numSIMS,numTIMES,2/),"float")       ; TC track in each simulation

  print_clock("Working on PV calculations!")

;====================================
; Start loop over ensemble members
;====================================

  do en = 0, dimsizes(ens_arr)-1

; List all files to read in and analyse
  diri        = "$sam/um/cp/ens/"
  fili_pr     = diri+"20160704T0000Z_ra1t_"+ens_arr(en)
  fili_b      = systemfunc("cd "+diri+" ; ls "+fili_pr+"*pb*.nc")
  numFILES_b  = dimsizes(fili_b)    ; Number of files (simulations)
  numTIMES_b  = 120                 ; Number of times in each simulation  
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
   wa    = f->dz_dt(:,6:14,{lat_0:lat_1},{lon_0:lon_1}) ; Dimensions [5 * 9 * 602 * 926]
   w	 = f->dz_dt(:,:,{lat_0:lat_1},{lon_0:lon_1})	; All pressure levels

 ; Read in horizontal velocity, relative vorticity and geopotential height
   u     = f->u(:,:,{lat_0:lat_1},{lon_0:lon_1})     ; Zonal velocity
   v     = f->v(:,:,{lat_0:lat_1},{lon_0:lon_1})     ; Meridional velocity
   vort  = f->rvor(:,:,{lat_0:lat_1},{lon_0:lon_1})  ; Relative vorticity
   z	 = f->ht(:,:,{lat_0:lat_1},{lon_0:lon_1})    ; Geopotential height 

   z     = z/10
   z@units = "dam"         ; Convert to decametres
   vort    = vort * (10 ^ -6) ; Convert to /s

   dy = lat(1) - lat(0)	   ; Grid spacing (latitude)
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

;====================================================================
; Find cyclone centre using both geopotential height and vorticity
;====================================================================

      ilev = 1
      vort_plane = vort(it,ilev,:,:)
      geo_plane	 = z(it,ilev,:,:)

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
; Also read in lower-tropospheric PV [proxy for convective cells]
;===================================================================

      ilev = 3					; Look at 850 hPa

    ; Plot either layer-averaged or single-level 'w'
      if (lay.eq.1) then
       w_ave = dim_avg_n_Wrap(wa(it,:,:,:),0)    ; Average over levels (700-300 hPa)
       vv_arr(en,ct,:,:) = w_ave                 ; Input values into vvel array
       pv_ave = dim_avg_n_Wrap(pv(it,1:4,:,:),0) ; Average PV (950-800 hPa)
       pv_arr(en,ct,:,:) = pv_ave
      else
       pv_arr(en,ct,:,:) = pv(it,ilev,:,:)
       vv_arr(en,ct,:,:) = w(it,ilev,:,:)        ; Input values into vvel array
      end if 

      print("en = "+en)
      print("ct = "+ct)
      if(any(ismissing(pv(it,ilev,:,:)))) then
       print("Missing values in PV array...")
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

   print("Working on panel plots!")
   times = 120   		  

   do it = 23, times-1, int

   ; Read 850 hPa PV values
     pv_plane00 = pv_arr(0,it,:,:)
     pv_plane01 = pv_arr(1,it,:,:)
     pv_plane02 = pv_arr(2,it,:,:)
     pv_plane03 = pv_arr(3,it,:,:)
     pv_plane04 = pv_arr(4,it,:,:)
     pv_plane05 = pv_arr(5,it,:,:)
     pv_plane06 = pv_arr(6,it,:,:)
     pv_plane07 = pv_arr(7,it,:,:)
     pv_plane08 = pv_arr(8,it,:,:)
     pv_plane09 = pv_arr(9,it,:,:)
     pv_plane10 = pv_arr(10,it,:,:)
     pv_plane11 = pv_arr(11,it,:,:)
  
   ; Read vertical velocity
     vv_plane00 = vv_arr(0,it,:,:)
     vv_plane01 = vv_arr(1,it,:,:)
     vv_plane02 = vv_arr(2,it,:,:)
     vv_plane03 = vv_arr(3,it,:,:)
     vv_plane04 = vv_arr(4,it,:,:)
     vv_plane05 = vv_arr(5,it,:,:)
     vv_plane06 = vv_arr(6,it,:,:)
     vv_plane07 = vv_arr(7,it,:,:)
     vv_plane08 = vv_arr(8,it,:,:)
     vv_plane09 = vv_arr(9,it,:,:)
     vv_plane10 = vv_arr(10,it,:,:)
     vv_plane11 = vv_arr(11,it,:,:)

     panel = new(12,graphic)	; Panel plot with 6 images

;=================================
; Output file type and location
;=================================

     if (lay.eq.1) then
      output = "$sam/nepartak/images/cart_coords/panel_pv850_lay_paper_"+time_arr(it)
      wks = gsn_open_wks(opt,output)     
     else
      output = "$sam/nepartak/images/cart_coords/panel_pv850_paper_"+time_arr(it)
      wks = gsn_open_wks(opt,output)
     end if

   ; Load colour table
     gsn_define_colormap(wks,"circ_budget")

;==========================
; Options for plotting
;==========================

   ; 850 hPa potential vorticity
     opts_em00                              = True
     opts_em00@cnFillOn                     = True
     opts_em00@cnLineLabelInterval          = 2.0
     opts_em00@cnLineLabelFontHeightF       = 0.012
     opts_em00@cnLineLabelBackgroundColor   = "transparent"
     opts_em00@cnLineLabelPlacementMode     = "constant"
     opts_em00@cnLinesOn                    = False ; Contour lines off
     opts_em00@cnInfoLabelOn                = False ; Contour labels off
     opts_em00@cnLevelSelectionMode         = "ExplicitLevels"
;     opts_em00@cnLevels                     = (/9,12,15,20,25,30,40,50/)
     opts_em00@cnLevels                     = (/15,20,25,30,35,40,45,50/)
;     opts_em00@cnLevels                     = (/20,25,30,35,40,45,50,60/) 
     opts_em00@cnFillColors                 = (/0,13,15,17,19,20,21,22,23/)
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
     opts_em00@tmXBLabelFontHeightF         = 0.021
     opts_em00@tmYLLabelFontHeightF         = 0.021
     opts_em00@tmXTOn                       = "False"    ; No tickmarks on top x-axis
     opts_em00@tmYROn                       = "False"    ; No tickmarks on right y-axis
     opts_em00@gsnMaximize                  = True       ; Maximise plot size 

     r					    = 0.50	 ; Radius of plot

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
     opts_vv00@cnLevels                      = (/0.6,0.8,1.0,1.2/)
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

;===========================================
; Additional options for ensemble members
;===========================================

     opts_em01                              = opts_em00
     opts_em01@mpLeftCornerLatF             = centre(1,it,0)-r
     opts_em01@mpLeftCornerLonF             = centre(1,it,1)-r
     opts_em01@mpRightCornerLatF            = centre(1,it,0)+r
     opts_em01@mpRightCornerLonF            = centre(1,it,1)+r

     opts_em02                              = opts_em00
     opts_em02@mpLeftCornerLatF             = centre(2,it,0)-r
     opts_em02@mpLeftCornerLonF             = centre(2,it,1)-r
     opts_em02@mpRightCornerLatF            = centre(2,it,0)+r
     opts_em02@mpRightCornerLonF            = centre(2,it,1)+r

     opts_em03                              = opts_em00
     opts_em03@mpLeftCornerLatF             = centre(3,it,0)-r
     opts_em03@mpLeftCornerLonF             = centre(3,it,1)-r
     opts_em03@mpRightCornerLatF            = centre(3,it,0)+r
     opts_em03@mpRightCornerLonF            = centre(3,it,1)+r

     opts_em04                              = opts_em00
     opts_em04@mpLeftCornerLatF             = centre(4,it,0)-r
     opts_em04@mpLeftCornerLonF             = centre(4,it,1)-r
     opts_em04@mpRightCornerLatF            = centre(4,it,0)+r
     opts_em04@mpRightCornerLonF            = centre(4,it,1)+r

     opts_em05                              = opts_em00
     opts_em05@mpLeftCornerLatF             = centre(5,it,0)-r
     opts_em05@mpLeftCornerLonF             = centre(5,it,1)-r
     opts_em05@mpRightCornerLatF            = centre(5,it,0)+r
     opts_em05@mpRightCornerLonF            = centre(5,it,1)+r

     opts_em06                              = opts_em00
     opts_em06@mpLeftCornerLatF             = centre(6,it,0)-r
     opts_em06@mpLeftCornerLonF             = centre(6,it,1)-r
     opts_em06@mpRightCornerLatF            = centre(6,it,0)+r
     opts_em06@mpRightCornerLonF            = centre(6,it,1)+r

     opts_em07                              = opts_em00
     opts_em07@mpLeftCornerLatF             = centre(7,it,0)-r
     opts_em07@mpLeftCornerLonF             = centre(7,it,1)-r
     opts_em07@mpRightCornerLatF            = centre(7,it,0)+r
     opts_em07@mpRightCornerLonF            = centre(7,it,1)+r

     opts_em08                              = opts_em00
     opts_em08@mpLeftCornerLatF             = centre(8,it,0)-r
     opts_em08@mpLeftCornerLonF             = centre(8,it,1)-r
     opts_em08@mpRightCornerLatF            = centre(8,it,0)+r
     opts_em08@mpRightCornerLonF            = centre(8,it,1)+r

     opts_em09                              = opts_em00
     opts_em09@mpLeftCornerLatF             = centre(9,it,0)-r
     opts_em09@mpLeftCornerLonF             = centre(9,it,1)-r
     opts_em09@mpRightCornerLatF            = centre(9,it,0)+r
     opts_em09@mpRightCornerLonF            = centre(9,it,1)+r

     opts_em10                              = opts_em00
     opts_em10@mpLeftCornerLatF             = centre(10,it,0)-r
     opts_em10@mpLeftCornerLonF             = centre(10,it,1)-r
     opts_em10@mpRightCornerLatF            = centre(10,it,0)+r
     opts_em10@mpRightCornerLonF            = centre(10,it,1)+r

     opts_em11                              = opts_em00
     opts_em11@mpLeftCornerLatF             = centre(11,it,0)-r
     opts_em11@mpLeftCornerLonF             = centre(11,it,1)-r
     opts_em11@mpRightCornerLatF            = centre(11,it,0)+r
     opts_em11@mpRightCornerLonF            = centre(11,it,1)+r

;=============
; Plot data 
;=============

     plot_pv0     = gsn_csm_contour_map(wks,pv_plane00,opts_em00)      ; Panel 1
     plot_vv00    = gsn_csm_contour(wks,vv_plane00,opts_vv00)          ; Omega
     overlay(plot_pv0, plot_vv00)

     plot_pv1     = gsn_csm_contour_map(wks,pv_plane01,opts_em01)      ; Panel 2
     plot_vv01    = gsn_csm_contour(wks,vv_plane01,opts_vv00)          ; Omega
     overlay(plot_pv1, plot_vv01)

     plot_pv2     = gsn_csm_contour_map(wks,pv_plane02,opts_em02)      ; Panel 3
     plot_vv02    = gsn_csm_contour(wks,vv_plane02,opts_vv00)          ; Omega
     overlay(plot_pv2, plot_vv02)

     plot_pv3     = gsn_csm_contour_map(wks,pv_plane03,opts_em03)      ; Panel 4
     plot_vv03    = gsn_csm_contour(wks,vv_plane03,opts_vv00)          ; Omega
     overlay(plot_pv3, plot_vv03)

     plot_pv4     = gsn_csm_contour_map(wks,pv_plane04,opts_em04)      ; Panel 5	
     plot_vv04    = gsn_csm_contour(wks,vv_plane04,opts_vv00)          ; Omega
     overlay(plot_pv4, plot_vv04)

     plot_pv5     = gsn_csm_contour_map(wks,pv_plane05,opts_em05)      ; Panel 6
     plot_vv05    = gsn_csm_contour(wks,vv_plane05,opts_vv00)          ; Omega
     overlay(plot_pv5, plot_vv05)

     plot_pv6     = gsn_csm_contour_map(wks,pv_plane06,opts_em06)      ; Panel 7
     plot_vv06    = gsn_csm_contour(wks,vv_plane06,opts_vv00)          ; Omega
     overlay(plot_pv6, plot_vv06)

     plot_pv7     = gsn_csm_contour_map(wks,pv_plane07,opts_em07)      ; Panel 8
     plot_vv07    = gsn_csm_contour(wks,vv_plane07,opts_vv00)          ; Omega
     overlay(plot_pv7, plot_vv07)

     plot_pv8     = gsn_csm_contour_map(wks,pv_plane08,opts_em08)      ; Panel 9
     plot_vv08    = gsn_csm_contour(wks,vv_plane08,opts_vv00)          ; Omega
     overlay(plot_pv8, plot_vv08)

     plot_pv9     = gsn_csm_contour_map(wks,pv_plane09,opts_em09)      ; Panel 10
     plot_vv09    = gsn_csm_contour(wks,vv_plane09,opts_vv00)          ; Omega
     overlay(plot_pv9, plot_vv09)

     plot_pv10    = gsn_csm_contour_map(wks,pv_plane10,opts_em10)      ; Panel 11
     plot_vv10    = gsn_csm_contour(wks,vv_plane10,opts_vv00)          ; Omega
     overlay(plot_pv10, plot_vv10)

     plot_pv11    = gsn_csm_contour_map(wks,pv_plane11,opts_em11)      ; Panel 12
     plot_vv11    = gsn_csm_contour(wks,vv_plane11,opts_vv00)          ; Omega
     overlay(plot_pv11, plot_vv11)

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    panel(0)   = plot_pv0
    panel(1)   = plot_pv1
    panel(2)   = plot_pv2
    panel(3)   = plot_pv3
    panel(4)   = plot_pv4
    panel(5)   = plot_pv5
    panel(6)   = plot_pv6
    panel(7)   = plot_pv7
    panel(8)   = plot_pv8
    panel(9)   = plot_pv9
    panel(10)  = plot_pv10
    panel(11)  = plot_pv11

    optsP                       = True
    optsP@gsnFrame              = False    ; Do not advance the frame
    optsP@gsnPanelLabelBar      = True     ; Turn off panel labelbar
    optsP@pmLabelBarWidthF      = 0.6
    optsP@pmLabelBarHeightF     = 0.15
    optsP@lbLabelFontHeightF    = 0.0125
    optsP@lbPerimOn             = False
    optsP@lbLabelFont           = "Helvetica"
    optsP@lbTitleString         = "850 hPa potential vorticity (PVU)"
    optsP@lbTitleFontHeightF    = 0.015
    optsP@lbTitleDirection      = "Across"
    optsP@lbTitlePosition       = "Bottom"
    optsP@txString              = "850 hPa PV and vertical velocity: "+title_arr(it)
    optsP@gsnPanelFigureStrings = (/"a) em00", "b) em01", "c) em02", "d) em03", \
                                    "e) em04", "f) em05", "g) em06", "h) em07", \
                                    "i) em08", "j) em09", "k) em10", "l) em11"/)
    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.98
    optsP@gsnPanelBottom        = 0.02
    optsP@amJust                = "TopLeft"
    optsP@gsnPanelFigureStringsFontHeightF = 0.0125 ; Reduce size of labels (0.01)

    gsn_panel(wks,panel,(/3,4/),optsP)             ; Draw as a single plot
    frame(wks)

   end do     ; End time loop

end 