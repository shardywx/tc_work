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

; Prelash (define arrays for loops below)
  ens_arr  = (/"em01","em02","em05",\
               "em07","em08","em09"/)

  numSIMS     = dimsizes(ens_arr)
  numTIMES    = 120

  time_arr    = new(numTIMES,string)
  title_arr   = new(numTIMES,string)
  centre      = new((/numSIMS,numTIMES,2/),"float")         ; TC track in each simulation

  do en = 0, dimsizes(ens_arr)-1

; List all files to read in and analyse
  diri        = "$sam/um/cp/ens/"				   ; Path to directory
  fili_a      = diri+"20160704T0000Z_ra1t_"+ens_arr(en)+"_pa.nc"   ; Only one file
  numTIMES    = 120		    				   ; No. of times in file

;===============================
; Start loop over input files
;===============================

  f = addfile(fili_a,"r")		; Read in data for chosen simulation (en=0,1,etc)

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 1000000000
  end setvalues

;==================================
; Get the variables we will need
;==================================

 ;================================================
 ; Number of times in different variable arrays
 ;================================================

   time  = f->t            ; Times in file [119]
   times = dimsizes(time)

   time1 = f->t_1          ; Times in file [120]
   t1    = dimsizes(time1)

   time2 = f->t_2          ; Times in file [120]
   t2    = dimsizes(time2)

   time3 = f->t_3          ; Times in file [120]
   t3    = dimsizes(time3)

   time4 = f->t_4          ; Times in file [121]
   t4     = dimsizes(time4)

 ;=======================
 ; All other variables
 ;=======================

   lon   = f->longitude    ; longitude (1098 points --> 109.04 to 152.92 degrees E)
   lat   = f->latitude     ; latitude  (810 points --> 1.8 to 34.16 degrees N)
   lon1  = f->longitude_1  ; longitude (1098 points --> 109.02 to 152.90 degrees E)
   lat1  = f->latitude_1   ; latitude (811 points --> 1.78 to 34.18 degrees N)

   mslp  = f->p            ; Mean sea level pressure (Pa) [time2 | z | lat | lon]
   mslp  = mslp / 100      ; Convert to hPa
   mslp@units = "hPa"

   dx    = lon(1) - lon(0) ; Grid spacing (longitude)
   dy    = lat(1) - lat(0) ; Grid spacing (latitude)

;=====================================================
; Create correct date strings for each output time
;=====================================================

   month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
                     "Jul","Aug","Sep","Oct","Nov","Dec"/)

   utc_date = cd_calendar(time2, 0)

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

   date_str = new(t2,string)
   time_str = new(t2,string)
   out_str  = new(t2,string)

   ct	    = 0		        ; Counter variable (time)

;===========================
; Loop over times in file
;===========================

   do it = 0, t2-1, int		; Loop over times in file

     date_str(it) = sprinti("%0.2iUTC ", hour(it)) + \
                    sprinti("%0.2i ", day(it)) \
                    + month_abbr(month(it))
     time_str(it) = sprinti("%0.2i UTC ", hour(it)) + \
                    sprinti("%0.2i ", day(it)) \
                    + month_abbr(month(it))
     out_str(it)  = sprinti("%0.2i", day(it)) + \
                    month_abbr(month(it)) + \
                    "_" + sprinti("%0.2iZ", hour(it))

     print("ct = "+ct)		    
     time_arr(ct) = out_str(it)

     print("Working on time: "+time_str(it))
     title_arr(ct) = "Valid at "+time_str(it)

;===========================================================
; Find cyclone centre on subset of global grid (use MSLP)
;===========================================================

    ; Choose grid subset
      lat_0 = 6.98
      lat_1 = 31.04
      lon_0 = 110.02
      lon_1 = 147.04

    ; Account for different sizes of arrays [u,v one element larger than MSLP]
      if (.not.ismissing(it).and.it.eq.0) then
        slp_plane  = mslp(it,0,{lat_0:lat_1},{lon_0:lon_1})
      else if (.not.ismissing(it).and.it.gt.0) then
        slp_plane  = mslp(it-1,0,{lat_0:lat_1},{lon_0:lon_1}) ; MSLP (hPa)
      end if
      end if

    ; Minimum sea level pressure
      slp_min	    = min(slp_plane)
      slp_smth	    = smth9_Wrap(slp_plane, 0.5, 0.5, True)
      slp_min_smth  = min(slp_smth)
      
      dims_h	    = dimsizes(slp_plane)
      slp1d	    = ndtooned(slp_plane)
      inds_h	    = ind_resolve(minind(slp1d),dims_h)
      slp1d_sm	    = ndtooned(slp_smth)
      inds_h_sm	    = ind_resolve(minind(slp1d_sm),dims_h)            

      lat_minslp = lat_0 + (dy * inds_h(0,0))
      lon_minslp = lon_0 + (dx * inds_h(0,1))
      print("Cyclone centre (MSLP): "+lat_minslp+" degrees N, "+lon_minslp+" degrees E")

    ;=========================================================
    ; EMPLOY A SAFETY NET IF WE IDENTIFY THE WRONG CYCLONE
    ;=========================================================
      if (.not.ismissing(ct).and.ct.eq.0) then
       centre(en,ct,0) = lat_minslp
       centre(en,ct,1) = lon_minslp

      else if (.not.ismissing(ct).and.ct.gt.0) then
       lt0  = centre(en,ct-1,0)             ; Previous TC latitude
       ln0  = centre(en,ct-1,1)             ; Previous TC longitude
       print("Previous TC latitude = "+lt0+" ; Previous TC longitude = "+ln0)
       d_lat = abs(lat_minslp-lt0)       ; Change in latitude (t1 - t0)
       d_lon = abs(lon_minslp-ln0)       ; Change in longitude
       print("Change in latitude = "+d_lat+" ; Change in longitude = "+d_lon)

     ; ; Recalculate TC centre if incorrect
       if (.not.ismissing(d_lat).and.d_lat.gt.dist.or.d_lon.gt.dist) then

      ; Create smaller grid [0.5 degrees]
        lt1 = lt0 - 0.5
        lt2 = lt0 + 0.5
        ln1 = ln0 - 0.5
        ln2 = ln0 + 0.5

        slp_plane0 = mslp(it-1,0,{lt1:lt2},{ln1:ln2})     ; Smaller grid (1.0 x 1.0 deg)
        slp_min0   = min(slp_plane0)                      ; Find minimum on smaller grid
        dims_p0    = dimsizes(slp_plane0)                 ; Size of smaller grid
        slp_1d0    = ndtooned(slp_plane0)                 ; Create 1-D array
        inds_p0    = ind_resolve(minind(slp_1d0),dims_p0) ; Find index of minimum

        lat_minslp = lt1 + (dy * inds_p0(0,0))
        lon_minslp = ln1 + (dy * inds_p0(0,1))
        print("New centre: "+lat_minslp+" degrees N, "+lon_minslp+" degrees E")
        centre(en,ct,0) = lat_minslp
        centre(en,ct,1) = lon_minslp

        delete([/slp_plane0, slp_min0, dims_p0, slp_1d0, inds_p0/])

       else if (.not.ismissing(d_lat).and.d_lat.le.dist.or.d_lon.le.dist) then
        print("Cyclone centre: "+lat_minslp+" degrees N, "+lon_minslp+" degrees E")
        centre(en,ct,0) = lat_minslp
        centre(en,ct,1) = lon_minslp
       end if
       end if

      end if
      end if

      delete([/dims_h,slp1d,inds_h,slp1d_sm,inds_h_sm/])

      ct = ct + 1	; Counter variable (time)

    end do		; End of time loop (do it = 0, t2-1, int)
 
  end do		; End of input file loop (do en = 0, dimsizes(ens_arr)-1)

;=======================================
; Tidy up before reading in new files 
;=======================================

  delete([/time,times,time1,t1,time2,t2,time3,t3,time4,t4/])
  delete([/lon,lat,lon1,lat1,mslp,dx,dy/])
  delete([/utc_date,year,month,day,hour,minute,second/])
  delete([/date_str,time_str,out_str/])

;==========================================================================
; Now read in PV and vertical velocity; start loop over ensemble members 
;==========================================================================

; Prelash (define arrays for loops below)
  ens_arr     = (/"em01","em02","em05",\
                  "em07","em08","em09"/)
  pv_arr      = new((/numSIMS,numTIMES,602,926/),"float") ; PV array
  vvel_arr    = new((/numSIMS,numTIMES,602,926/),"float") ; Vertical velocity array

  print("Working on PV calculations!")

;================================================
; Create colourmap based on Brewer pattern
;================================================

 CBR_wet = (/ "(/1.0000,    1.0000,    1.0000   /)", \
              "(/0.0000,    0.0000,    0.0000   /)", \
              "(/0.6588,    0.8667,    0.7098   /)", \ ;         c < -15.0 (Green 1)
              "(/0.2286,    0.8000,    0.3170   /)", \ ; -15.0 < c < -10.0 (Green 2)
              "(/0.0000,    0.6000,    0.1500   /)", \ ; -10.0 < c < -5.0  (Green 3)
	      "(/1.0000,    1.0000,    1.0000   /)", \ ; -5.0  < c < 5.0   (White)
              "(/0.3059,    0.7020,    0.8275   /)", \ ; 5.0   < c < 10.0  (Blue 1)
              "(/0.1686,    0.5490,    0.7451   /)", \ ; 10.0  < c < 15.0  (Blue 2)
              "(/0.0314,    0.4078,    0.6745   /)", \ ; 15.0  < c < 20.0  (Blue 3)
              "(/0.0314,    0.2510,    0.5059   /)", \ ; 20.0  < c < 25.0  (Blue 4)
              "(/0.6140,    0.5986,    0.7817   /)", \ ; 25.0  < c < 30.0  (Purple 1)
              "(/0.4272,    0.3407,    0.6512   /)", \ ; 30.0  < c < 35.0  (Purple 2)
              "(/0.4001,    0.2000,    0.6500   /)" /) ;         c > 35.0  (Purple 3)

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
   w = f->dz_dt(:,6:14,{lat_0:lat_1},{lon_0:lon_1})     ; Dimensions [5 * 9 * 602 * 926]

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

      print("Working on time: "+time_str(it))

;===========================
; Calculate layer average 
;===========================

      w_box = w(it,:,:,:)			; Look at each time individually
      w_ave = dim_avg_n_Wrap(w_box,0)		; Average over levels (700-300 hPa)

      vvel_arr(en,ct,:,:) = w_ave		; Input values into vvel array

;===================================================================
; Also read in lower-tropospheric PV [proxy for convective cells]
;===================================================================

      pv_plane = pv(it,3,:,:)			; 850 hPa PV 
      pv_arr(en,ct,:,:) = pv_plane		; Input values into PV array      

      print("en = "+en)
      print("ct = "+ct)
      if(any(ismissing(pv_plane))) then
       print("Missing values in 'pv_plane'...")
      end if

      ct = ct + 1				; Counter variable (time)

    end do	; End time loop (do it = 0, times-1)

    delete([/time,times,w,pv/])        	        ; Tidy up (avoid dimension size errors)
    delete([/utc_date,year,month,day,hour,minute,second/])
    delete([/date_str,time_str,out_str/])
 
   end do	; End file loop (do nf = 0, numFILES_b-1)

  end do	; End ensemble member loop (do en = 0, dimsizes(ens_arr)-1)

;===============================================
; Produce panel plots for each time interval
;===============================================

   print("Working on panel plots!")
   times = 120   		  

   do it = 24, times-1, int

   ; Read 850 hPa PV values
     pv_plane00 = pv_arr(0,it,:,:)
     pv_plane01 = pv_arr(1,it,:,:)
     pv_plane02 = pv_arr(2,it,:,:)
     pv_plane03 = pv_arr(3,it,:,:)
     pv_plane04 = pv_arr(4,it,:,:)
     pv_plane05 = pv_arr(5,it,:,:)

     if(any(ismissing(pv_plane00))) then
      print("Some missing values in 'pv_plane00'...")
     end if
     if(any(ismissing(pv_plane01))) then
      print("Some missing values in 'pv_plane01'...")
     end if
     if(any(ismissing(pv_plane02))) then
      print("Some missing values in 'pv_plane02'...")
     end if
     if(any(ismissing(pv_plane03))) then
      print("Some missing values in 'pv_plane03'...")
     end if
     if(any(ismissing(pv_plane04))) then
      print("Some missing values in 'pv_plane04'...")
     end if
     if(any(ismissing(pv_plane05))) then
      print("Some missing values in 'pv_plane05'...")
     end if

;   ; Read vertical velocity values 
;     vv_plane00 = vvel_arr(0,it,:,:)
;     vv_plane01 = vvel_arr(1,it,:,:)
;     vv_plane02 = vvel_arr(2,it,:,:)
;     vv_plane03 = vvel_arr(3,it,:,:)
;     vv_plane04 = vvel_arr(4,it,:,:)
;     vv_plane05 = vvel_arr(5,it,:,:)

;     panel = new(6,graphic)	; Panel plot with 6 images

;=================================
; Output file type and location
;=================================

     output = "$sam/nepartak/images/cart_coords/panel_pv_vvel_"+time_arr(it)
     wks = gsn_open_wks(opt,output)

;==============================================================
; Load NIMROD colour table for precip duration data
;==============================================================

   ; Load colour table
;     gsn_define_colormap(wks,CBR_wet)
     gsn_define_colormap(wks,"circ_budget")

;==========================
; Options for plotting
;==========================

;;;;;;;;;; EDIT OPTIONS ;;;;;;;;;; [turn off "Ertel PV..." label, labelbar, resize...]

   ; Rainfall rate
     opts_em00                              = True
     opts_em00@cnFillOn                     = True
     opts_em00@cnLineLabelInterval          = 2.0
     opts_em00@cnLineLabelFontHeightF       = 0.012
     opts_em00@cnLineLabelBackgroundColor   = "transparent"
     opts_em00@cnLineLabelPlacementMode     = "constant"
     opts_em00@cnLinesOn                    = False ; Contour lines off
     opts_em00@cnInfoLabelOn                = False ; Contour labels off
     opts_em00@cnLevelSelectionMode         = "ExplicitLevels"
     opts_em00@cnLevels			    = (/-8,-4,-2,2,4,8,12,16,20,25,30,40/)
;     opts_em00@cnFillColors		    = (/2,3,4,5,6,7,8,9,10,11,12/)
     opts_em00@cnFillColors                 = (/4,6,8,0,12,13,15,17,19,20,21,22,23/)

;     opts_em00@gsnPaperOrientation          = "landscape"
;     opts_em00@tiMainString                 = ""
;     opts_em00@tiMainFontHeightF            = 0.0125
;     opts_em00@gsnLeftString                = ""
;     opts_em00@gsnRightString               = ""

     opts_em00@lbLabelBarOn		    = False	; Turn off individual labelbars
;     opts_em00@pmLabelBarWidthF             = 0.6
;     opts_em00@pmLabelBarHeightF            = 0.15
;     opts_em00@lbLabelFontHeightF           = 0.01
;     opts_em00@lbPerimOn                    = False
;     opts_em00@lbLabelFont                  = "Helvetica"
;     opts_em00@lbTitleOn                    = True
;     opts_em00@lbTitleString                = "Potential vorticity (PVU)"
;     opts_em00@lbTitleFontHeightF           = 0.015
;     opts_em00@lbTitleDirection             = "Across"

;   ; Additional plotting resources
;     opts_em00@mpDataBaseVersion            = "Ncarg4_1" ; More recent database
;     opts_em00@mpDataSetName                = "Earth..4" ; High resolution
;     opts_em00@mpOutlineBoundarySets        = "National" ; National borders
;     opts_em00@mpGeophysicalLineColor       = "black"    ; Colour borders black
;     opts_em00@mpGeophysicalLineThicknessF  = 1.0        ; Border line thickness
;     opts_em00@mpGridAndLimbOn              = False      ; Turn on lat/lon lines
;     opts_em00@pmTickMarkDisplayMode        = "Always"   ; Turn on map tickmarks
;     opts_em00@tmXBMajorLengthF             = 0.005      ; Change tickmark length
;     opts_em00@tmXTOn                       = "False"    ; No tickmarks on top x-axis
;     opts_em00@tmYROn                       = "False"    ; No tickmarks on right y-axis
;     opts_em00@gsnMaximize                  = True       ; Maximise plot size 

     r					    = 0.50	 ; Radius of plot

     opts_em00@gsnAddCyclic                 = False
     opts_em00@mpLimitMode		    = "Corners"
     opts_em00@mpLeftCornerLatF             = centre(0,it,0)-r
     opts_em00@mpLeftCornerLonF             = centre(0,it,1)-r
     opts_em00@mpRightCornerLatF            = centre(0,it,0)+r
     opts_em00@mpRightCornerLonF            = centre(0,it,1)+r
     opts_em00@gsnDraw		            = True       ; Do not draw the plot
     opts_em00@gsnFrame		            = False      ; Do not advance the frame

     opts_em00@vpXF			    = 0.08
     opts_em00@vpYF			    = 0.94
     opts_em00@vpWidthF			    = 0.25
     opts_em00@vpHeightF		    = 0.25

;     plot_pv00 = gsn_csm_contour_map(wks,pv_plane00,opts_em00)     ; 850 hPa PV
     plot0     = gsn_csm_contour_map(wks,pv_plane00,opts_em00)      ; Panel 1

;===========================================
; Additional options for ensemble members
;===========================================

     opts_em00@mpLeftCornerLatF             = centre(1,it,0)-r
     opts_em00@mpLeftCornerLonF             = centre(1,it,1)-r
     opts_em00@mpRightCornerLatF            = centre(1,it,0)+r
     opts_em00@mpRightCornerLonF            = centre(1,it,1)+r     	
     opts_em00@vpXF                         = 0.38
     opts_em00@vpYF                         = 0.94
;     plot_pv01 = gsn_csm_contour_map(wks,pv_plane01,opts_em00)     ; 850 hPa PV
     plot1     = gsn_csm_contour_map(wks,pv_plane01,opts_em00)      ; Panel 2

     opts_em00@mpLeftCornerLatF             = centre(2,it,0)-r
     opts_em00@mpLeftCornerLonF             = centre(2,it,1)-r
     opts_em00@mpRightCornerLatF            = centre(2,it,0)+r
     opts_em00@mpRightCornerLonF            = centre(2,it,1)+r
     opts_em00@vpXF                         = 0.68
     opts_em00@vpYF                         = 0.94
;     plot_pv02 = gsn_csm_contour_map(wks,pv_plane02,opts_em00)     ; 850 hPa PV
     plot2     = gsn_csm_contour_map(wks,pv_plane02,opts_em00)      ; Panel 3

     opts_em00@mpLeftCornerLatF             = centre(3,it,0)-r
     opts_em00@mpLeftCornerLonF             = centre(3,it,1)-r
     opts_em00@mpRightCornerLatF            = centre(3,it,0)+r
     opts_em00@mpRightCornerLonF            = centre(3,it,1)+r
     opts_em00@vpXF                         = 0.08
     opts_em00@vpYF                         = 0.56
;     plot_pv03 = gsn_csm_contour_map(wks,pv_plane03,opts_em00)     ; 850 hPa PV
     plot3     = gsn_csm_contour_map(wks,pv_plane03,opts_em00)      ; Panel 4

     opts_em00@mpLeftCornerLatF             = centre(4,it,0)-r
     opts_em00@mpLeftCornerLonF             = centre(4,it,1)-r
     opts_em00@mpRightCornerLatF            = centre(4,it,0)+r
     opts_em00@mpRightCornerLonF            = centre(4,it,1)+r
     opts_em00@vpXF                         = 0.38
     opts_em00@vpYF                         = 0.56
;     plot_pv04 = gsn_csm_contour_map(wks,pv_plane04,opts_em00)     ; 850 hPa PV
     plot4     = gsn_csm_contour_map(wks,pv_plane04,opts_em00)      ; Panel 5

     opts_em00@mpLeftCornerLatF             = centre(5,it,0)-r
     opts_em00@mpLeftCornerLonF             = centre(5,it,1)-r
     opts_em00@mpRightCornerLatF            = centre(5,it,0)+r
     opts_em00@mpRightCornerLonF            = centre(5,it,1)+r
     opts_em00@vpXF                         = 0.68
     opts_em00@vpYF                         = 0.56
;     plot_pv05 = gsn_csm_contour_map(wks,pv_plane05,opts_em00)     ; 850 hPa PV
     plot5     = gsn_csm_contour_map(wks,pv_plane05,opts_em00)      ; Panel 6

;     drawNDCGrid(wks)
     frame(wks)

    end do     ; End time loop

end 