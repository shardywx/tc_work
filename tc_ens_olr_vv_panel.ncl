; Script to plot MSLP, 10-m winds and precipitation following a TC 
; Run using e.g. [ncl int=1 dist=1 'opt="x11"' tc_ens_olr_panel.ncl]
; Where 'int' is the interval in the time loop, 'dist' is the threshold for the cyclone
; tracking code, and 'opt' is the output file format

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"

begin

;======================================================================
; Read in vertical velocity; start loop over ensemble members
;======================================================================

; Prelash (define arrays for loops below)
  ens_arr     = (/"em00","em01","em02","em03","em04","em05",\
                  "em06","em07","em08","em09","em10","em11"/)
  numSIMS     = dimsizes(ens_arr)
  numTIMES    = 120
  time_arr    = new(numTIMES,string)
  title_arr   = new(numTIMES,string)
  vv_arr      = new((/numSIMS,numTIMES,602,926/),"float") ; Vertical velocity array
  centre      = new((/numSIMS,numTIMES,2/),"float")       ; TC track in each simulation
  markers     = new((/numSIMS,numTIMES/),graphic)	  ; Add TC marker to plots

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
  ct          = 0                   ; Counter variable (time)

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

  ; Read in vertical velocity (m/s) and concentrate on 700-300 hPa over a subset
    w  = f->dz_dt(:,:,{lat_0:lat_1},{lon_0:lon_1})

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

;=============================
; Read in vertical velocity
;=============================

      ilev = 16		; 200 hPa vertical velocity
      p = pres(ilev)
      pstr = p+"hPa"

      w_plane = w(it,ilev,:,:)                  ; Vertical velocity on chosen level
      vv_arr(en,ct,:,:) = w_plane               ; Input values into larger array

      print("en = "+en)
      print("ct = "+ct)
      if(any(ismissing(w_plane))) then
       print("Missing values in 'w_plane'...")
      end if

      ct = ct + 1                               ; Counter variable (time)

    end do      ; End time loop (do it = 0, times-1)

  ; Tidy up (avoid dimension size errors)
    delete([/time,times,w,u,v,vort,z,lat,lon,lat1,lon1/])
    delete([/utc_date,year,month,day,hour,minute,second/])
    delete([/date_str,time_str,out_str/])

   end do       ; End file loop (do nf = 0, numFILES_b-1)

  end do        ; End ensemble member loop (do en = 0, dimsizes(ens_arr)-1)

;=================================================
; Read in OLR; start loop over ensemble members
;=================================================

; List all files to read in and analyse
  diri        = "$sam/um/cp/ens/"
  fili_prefix = "20160704T0000Z_ra1t_em"
  fili_a      = systemfunc("cd "+diri+" ; ls "+fili_prefix+"*pa.nc")
  fili        = diri+fili_a

  numSIMS     = dimsizes(fili)      ; Number of files (simulations)
  numTIMES    = 120                 ; Number of times in each simulation
  time_arr    = new(numTIMES,string)
  title_arr   = new(numTIMES,string)
  olr_arr     = new((/numSIMS,numTIMES,602,926/),"float") ; OLR array

  print_clock("Working on OLR calculations!")

;===============================
; Start loop over input files
;===============================

  do nf = 0, numSIMS-1
  f = addfile(fili(nf),"r")

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

   lsrn  = f->lsrain       ; Rainfall rate (kg/m**2/s) [t1 | _ | lat | lon]
   lsrn  = lsrn * 3600     ; Convert to mm / hr
   lsrn@units = "mm hr~S~-1~N~"

   olr   = f->olr          ; Outgoing longwave radiation [time | TOA | lat | lon]

   u_varname = "x-wind"
   v_varname = "y-wind"

   u     = f->$u_varname$  ; 10-m zonal velocity (m/s) [time4 | z | lat_1 | lon_1]
   v     = f->$v_varname$  ; 10-m meridional velocity (m/s) [time4 | z | lat_1 | lon_1]
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

 ; Local time arrays [to add to final plots]
   day_loc  = new(120,"integer")
   hour_loc = new(120,"integer")

   day_loc(0:111) = day(8:119)
   day_loc(112:119) = 9

   hour_loc(0:111) = hour(8:119)
   hour_loc(112:119) = ispan(1,8,1)

 ; Correct for errors in the code (round up value of hour when minutes = 59)
   do it = 0, times-1
     if (minute(it).gt.30) then
       hour(it) = hour(it)+1
     end if
   end do

   date_str = new(t2,string)
   time_str = new(t2,string)
   loc_str  = new(t2,string)
   out_str  = new(t2,string)

   ct       = 0                 ; Counter variable (time)

;===========================
; Loop over times in file
;===========================

   do it = 0, t2-1, int         ; Loop over times in file

     date_str(it) = sprinti("%0.2iUTC ", hour(it)) + \
                    sprinti("%0.2i ", day(it)) \
                    + month_abbr(month(it))
     time_str(it) = sprinti("%0.2iZ ", hour(it)) + \
                    sprinti("%0.2i ", day(it)) \
                    + month_abbr(month(it))
     loc_str(it)  = sprinti("%0.2iZ ", hour_loc(it)) + \
                    sprinti("%0.2i ", day_loc(it)) \
                    + month_abbr(month(it))
     out_str(it)  = sprinti("%0.2i", day(it)) + \
                    month_abbr(month(it)) + \
                    "_" + sprinti("%0.2iZ", hour(it))

     print("ct = "+ct)
     time_arr(ct) = out_str(it)

     print("Working on time: "+time_str(it))
     title_arr(ct) = "Valid at "+time_str(it)+" (local time "+loc_str(it)+")"

;===========================================================
; Find cyclone centre on subset of global grid (use MSLP)
;===========================================================

    ; Choose grid subset
      lat_0 = 6.98
      lat_1 = 31.04
      lon_0 = 110.02
      lon_1 = 147.04

    ; Account for different sizes of arrays [u,v one element larger than MSLP]
      if (it.eq.0) then
        u_plane    = u(it,0,{lat_0:lat_1},{lon_0:lon_1})
        v_plane    = v(it,0,{lat_0:lat_1},{lon_0:lon_1})
	slp_plane  = mslp(it,0,{lat_0:lat_1},{lon_0:lon_1})
        olr_plane  = olr(it,0,{lat_0:lat_1},{lon_0:lon_1})    ; Conv. rain rate (mm/h)
      else
        u_plane    = u(it,0,{lat_0:lat_1},{lon_0:lon_1})      ; Zonal velocity (m/s)
        v_plane    = v(it,0,{lat_0:lat_1},{lon_0:lon_1})      ; Meridional velocity (m/s)
        slp_plane  = mslp(it-1,0,{lat_0:lat_1},{lon_0:lon_1}) ; MSLP (hPa)
        olr_plane  = olr(it-1,0,{lat_0:lat_1},{lon_0:lon_1})  ; Conv. rain rate (mm/h)
      end if

    ; Calculate 10-m wind speed
      speed      = sqrt(u_plane ^ 2.0 + v_plane ^ 2.0)    ; 10m windspeed (m/s)
      copy_VarCoords(u_plane, speed)
      copy_VarAtts(u_plane, speed)
      copy_VarMeta(u_plane, speed)
      speed@units = "m s~S~-1~N~"

    ; Input values into larger arrays
      olr_arr(nf,ct,:,:) = olr_plane

      ct = ct + 1       ; Counter variable (time)

    end do              ; End of time loop (do it = 0, t2-1, int)

  end do                ; End of input file loop (do nf = 0, numSIMS-1)

;===============================================
; Produce panel plots for each time interval
;===============================================

   print_clock("Working on panel plots!")

   do it = 60, 78, 1

     olr_plane00 = olr_arr(0,it,:,:)
     olr_plane01 = olr_arr(1,it,:,:)
     olr_plane02 = olr_arr(2,it,:,:)
     olr_plane03 = olr_arr(3,it,:,:)
     olr_plane04 = olr_arr(4,it,:,:)
     olr_plane05 = olr_arr(5,it,:,:)
     olr_plane06 = olr_arr(6,it,:,:)
     olr_plane07 = olr_arr(7,it,:,:)
     olr_plane08 = olr_arr(8,it,:,:)
     olr_plane09 = olr_arr(9,it,:,:)
     olr_plane10 = olr_arr(10,it,:,:)
     olr_plane11 = olr_arr(11,it,:,:)

   ; Read 200 hPa vertical velocity values
     vv_plane00 = vv_arr(0,it-1,:,:)
     vv_plane01 = vv_arr(1,it-1,:,:)
     vv_plane02 = vv_arr(2,it-1,:,:)
     vv_plane03 = vv_arr(3,it-1,:,:)
     vv_plane04 = vv_arr(4,it-1,:,:)
     vv_plane05 = vv_arr(5,it-1,:,:)
     vv_plane06 = vv_arr(6,it-1,:,:)
     vv_plane07 = vv_arr(7,it-1,:,:)
     vv_plane08 = vv_arr(8,it-1,:,:)
     vv_plane09 = vv_arr(9,it-1,:,:)
     vv_plane10 = vv_arr(10,it-1,:,:)
     vv_plane11 = vv_arr(11,it-1,:,:)

     panel = new(12,graphic)	; Panel plot with 12 images
     print_clock("Working on time "+it)

;=================================
; Output file type and location
;=================================

     output = "$sam/nepartak/images/cart_coords/panel_olr_vv_"+time_arr(it)
     wks = gsn_open_wks(opt,output)

;==============================================================
; Load NIMROD colour table for precip duration data
;==============================================================

;    ; Load colour table (Greyscale)
;      gsn_define_colormap(wks,"MPL_gist_yarg")

;   ; Select colours from the "MPL_gist_yarg" colour table [27 colours]
;     colours                               = (/2, 18, 34, 50, 66, \
;                                               74, 90, 106, 120/)
;     opts_em00@cnLevels                     = (/75,80,85,90,95,100,105,110/)

;    ; Load colour table (Greyscale + purple)
;      gsn_define_colormap(wks,"olr_new")

;    ; Select colours from the "MPL_gist_yarg" colour table [31 colours]
;      colours                              = (/137, 136, 135, 134, 133, \
;                                               2, 8, 14, 19, 24, 29, 33, 37, 41, 45, \
;                                               50, 56, 62, 68, 74, 80, 85, 90,\
;                                               95, 100, 106, 112/)

    ; Load colour table (Greyscale + purple)
      gsn_define_colormap(wks,"olr_new")

    ; Select colours from the "MPL_gist_yarg" colour table [31 colours]
      colours = (/2, 10, 18, 26, 36, 47, 58, 66, 74, 82, 90, 98, 106/)

;==========================
; Options for plotting
;==========================

   ; OLR
     opts_em00                              = True
     opts_em00@cnFillOn                     = True
     opts_em00@cnLineLabelInterval          = 2.0
     opts_em00@cnLineLabelFontHeightF       = 0.012
     opts_em00@cnLineLabelBackgroundColor   = "transparent"
     opts_em00@cnLineLabelPlacementMode     = "constant"
     opts_em00@cnLinesOn                    = False ; Contour lines off
     opts_em00@cnInfoLabelOn                = False ; Contour labels off
     opts_em00@cnLevelSelectionMode         = "ExplicitLevels"
     opts_em00@cnLevels                     = (/76,78,80,82,\
                                                84,86,88,90,\
                                                95,100,110/) ; 12 levels

;     opts_em00@cnLevels                     = (/80,82,84,86,\
;                                                88,90,95,100,\
;                                                105,110,120/) ; 12 levels

;     opts_em00@cnLevels                     = (/80,85,90,\
;                                                100,110,120,130,\
;                                                140,150,160,170,180/) ; 12 levels

;     opts_em00@cnLevels                     = (/65,70,75,80,85,90,95,\
;                                                100,105,110,115,120,125,130,\
;                                                140,150,160,170,180,190,200,\
;                                                210,220,230,240,250/) ; 26 levels

     opts_em00@cnFillColors                 = colours
     opts_em00@gsnPaperOrientation          = "landscape"
     opts_em00@tiMainString                 = ""
     opts_em00@tiMainFontHeightF            = 0.0125
     opts_em00@gsnLeftString                = ""
     opts_em00@gsnRightString               = ""

     opts_em00@lbLabelBarOn		    = False	; Turn off individual labelbars
     opts_em00@pmLabelBarWidthF             = 0.6
     opts_em00@pmLabelBarHeightF            = 0.15
     opts_em00@lbLabelFontHeightF           = 0.01
     opts_em00@lbPerimOn                    = False
     opts_em00@lbLabelFont                  = "Helvetica"
     opts_em00@lbTitleOn                    = True
     opts_em00@lbTitleString                = "OLR (W m~S~-2~N~)"
     opts_em00@lbTitleFontHeightF           = 0.015
     opts_em00@lbTitleDirection             = "Across"

   ; Additional plotting resources
     opts_em00@mpDataBaseVersion            = "Ncarg4_1" ; More recent database
     opts_em00@mpDataSetName                = "Earth..4" ; High resolution
     opts_em00@mpOutlineBoundarySets        = "National" ; National borders
     opts_em00@mpGeophysicalLineColor       = "black"    ; Colour borders black
     opts_em00@mpGeophysicalLineThicknessF  = 1.0        ; Border line thickness
     opts_em00@mpGridAndLimbOn              = False      ; Turn on lat/lon lines
     opts_em00@pmTickMarkDisplayMode        = "Always"   ; Turn on map tickmarks
     opts_em00@tmXBMajorLengthF             = 0.007      ; Change tickmark length
     opts_em00@tmXBLabelStride		    = 1		 ; Skip every other label (x)
     opts_em00@tmYLLabelStride		    = 1		 ; Skip every other label (y)
     opts_em00@tmXTOn                       = "False"    ; No tickmarks on top x-axis
     opts_em00@tmYROn                       = "False"    ; No tickmarks on right y-axis
     opts_em00@gsnMaximize                  = True       ; Maximise plot size 

     r					    = 2.50	 ; Determine radius of plot

     opts_em00@gsnAddCyclic                 = False
     opts_em00@mpLimitMode		    = "Corners"
     opts_em00@mpLeftCornerLatF             = centre(0,it-1,0)-r
     opts_em00@mpLeftCornerLonF             = centre(0,it-1,1)-r
     opts_em00@mpRightCornerLatF            = centre(0,it-1,0)+r
     opts_em00@mpRightCornerLonF            = centre(0,it-1,1)+r
     opts_em00@gsnDraw		            = False      ; Do not draw the plot
     opts_em00@gsnFrame		            = False      ; Do not advance the frame

   ; Vertical velocity
     opts_vv00                               = True
     opts_vv00@cnFillOn                      = False
     opts_vv00@cnLineColor                   = "dodgerblue"
     opts_vv00@cnInfoLabelOn                 = False
     opts_vv00@cnLineLabelsOn                = True
     opts_vv00@cnLineLabelInterval           = 2.0
     opts_vv00@cnLevelSelectionMode          = "ExplicitLevels"
     opts_vv00@cnLevels                      = (/1.0,2.0,3.0/)
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

     opts_em01				    = opts_em00
     opts_em01@mpLeftCornerLatF             = centre(1,it-1,0)-r
     opts_em01@mpLeftCornerLonF             = centre(1,it-1,1)-r
     opts_em01@mpRightCornerLatF            = centre(1,it-1,0)+r
     opts_em01@mpRightCornerLonF            = centre(1,it-1,1)+r     	

     opts_em02				    = opts_em00
     opts_em02@mpLeftCornerLatF             = centre(2,it-1,0)-r
     opts_em02@mpLeftCornerLonF             = centre(2,it-1,1)-r
     opts_em02@mpRightCornerLatF            = centre(2,it-1,0)+r
     opts_em02@mpRightCornerLonF            = centre(2,it-1,1)+r

     opts_em03				    = opts_em00
     opts_em03@mpLeftCornerLatF             = centre(3,it-1,0)-r
     opts_em03@mpLeftCornerLonF             = centre(3,it-1,1)-r
     opts_em03@mpRightCornerLatF            = centre(3,it-1,0)+r
     opts_em03@mpRightCornerLonF            = centre(3,it-1,1)+r

     opts_em04				    = opts_em00
     opts_em04@mpLeftCornerLatF             = centre(4,it-1,0)-r
     opts_em04@mpLeftCornerLonF             = centre(4,it-1,1)-r
     opts_em04@mpRightCornerLatF            = centre(4,it-1,0)+r
     opts_em04@mpRightCornerLonF            = centre(4,it-1,1)+r

     opts_em05				    = opts_em00
     opts_em05@mpLeftCornerLatF             = centre(5,it-1,0)-r
     opts_em05@mpLeftCornerLonF             = centre(5,it-1,1)-r
     opts_em05@mpRightCornerLatF            = centre(5,it-1,0)+r
     opts_em05@mpRightCornerLonF            = centre(5,it-1,1)+r

     opts_em06				    = opts_em00
     opts_em06@mpLeftCornerLatF             = centre(6,it-1,0)-r
     opts_em06@mpLeftCornerLonF             = centre(6,it-1,1)-r
     opts_em06@mpRightCornerLatF            = centre(6,it-1,0)+r
     opts_em06@mpRightCornerLonF            = centre(6,it-1,1)+r

     opts_em07				    = opts_em00
     opts_em07@mpLeftCornerLatF             = centre(7,it-1,0)-r
     opts_em07@mpLeftCornerLonF             = centre(7,it-1,1)-r
     opts_em07@mpRightCornerLatF            = centre(7,it-1,0)+r
     opts_em07@mpRightCornerLonF            = centre(7,it-1,1)+r

     opts_em08				    = opts_em00
     opts_em08@mpLeftCornerLatF             = centre(8,it-1,0)-r
     opts_em08@mpLeftCornerLonF             = centre(8,it-1,1)-r
     opts_em08@mpRightCornerLatF            = centre(8,it-1,0)+r
     opts_em08@mpRightCornerLonF            = centre(8,it-1,1)+r

     opts_em09				    = opts_em00
     opts_em09@mpLeftCornerLatF             = centre(9,it-1,0)-r
     opts_em09@mpLeftCornerLonF             = centre(9,it-1,1)-r
     opts_em09@mpRightCornerLatF            = centre(9,it-1,0)+r
     opts_em09@mpRightCornerLonF            = centre(9,it-1,1)+r

     opts_em10				    = opts_em00
     opts_em10@mpLeftCornerLatF             = centre(10,it-1,0)-r
     opts_em10@mpLeftCornerLonF             = centre(10,it-1,1)-r
     opts_em10@mpRightCornerLatF            = centre(10,it-1,0)+r
     opts_em10@mpRightCornerLonF            = centre(10,it-1,1)+r

     opts_em11				    = opts_em00
     opts_em11@mpLeftCornerLatF             = centre(11,it-1,0)-r
     opts_em11@mpLeftCornerLonF             = centre(11,it-1,1)-r
     opts_em11@mpRightCornerLatF            = centre(11,it-1,0)+r
     opts_em11@mpRightCornerLonF            = centre(11,it-1,1)+r

;=============
; Plot data 
;=============

    plot_olr00 = gsn_csm_contour_map(wks,olr_plane00,opts_em00)
;    plot_vv00  = gsn_csm_contour(wks,vv_plane00,opts_vv00)
;    overlay(plot_olr00, plot_vv00)

    plot_olr01 = gsn_csm_contour_map(wks,olr_plane01,opts_em01)
;    plot_vv01  = gsn_csm_contour(wks,vv_plane01,opts_vv00)
;    overlay(plot_olr01,	plot_vv01)

    plot_olr02 = gsn_csm_contour_map(wks,olr_plane02,opts_em02)
;    plot_vv02  = gsn_csm_contour(wks,vv_plane02,opts_vv00)
;    overlay(plot_olr02,	plot_vv02)

    plot_olr03 = gsn_csm_contour_map(wks,olr_plane03,opts_em03)
;    plot_vv03  = gsn_csm_contour(wks,vv_plane03,opts_vv00)
;    overlay(plot_olr03,	plot_vv03)

    plot_olr04 = gsn_csm_contour_map(wks,olr_plane04,opts_em04)
;    plot_vv04  = gsn_csm_contour(wks,vv_plane04,opts_vv00)
;    overlay(plot_olr04,	plot_vv04)

    plot_olr05 = gsn_csm_contour_map(wks,olr_plane05,opts_em05)
;    plot_vv05  = gsn_csm_contour(wks,vv_plane05,opts_vv00)
;    overlay(plot_olr05,	plot_vv05)

    plot_olr06 = gsn_csm_contour_map(wks,olr_plane06,opts_em06)
;    plot_vv06  = gsn_csm_contour(wks,vv_plane06,opts_vv00)
;    overlay(plot_olr06,	plot_vv06)

    plot_olr07 = gsn_csm_contour_map(wks,olr_plane07,opts_em07)
;    plot_vv07  = gsn_csm_contour(wks,vv_plane07,opts_vv00)
;    overlay(plot_olr07,	plot_vv07)

    plot_olr08 = gsn_csm_contour_map(wks,olr_plane08,opts_em08)
;    plot_vv08  = gsn_csm_contour(wks,vv_plane08,opts_vv00)
;    overlay(plot_olr08,	plot_vv08)

    plot_olr09 = gsn_csm_contour_map(wks,olr_plane09,opts_em09)
;    plot_vv09  = gsn_csm_contour(wks,vv_plane09,opts_vv00)
;    overlay(plot_olr09,	plot_vv09)

    plot_olr10 = gsn_csm_contour_map(wks,olr_plane10,opts_em10)
;    plot_vv10  = gsn_csm_contour(wks,vv_plane10,opts_vv00)
;    overlay(plot_olr10,	plot_vv10)

    plot_olr11 = gsn_csm_contour_map(wks,olr_plane11,opts_em11)
;    plot_vv11  = gsn_csm_contour(wks,vv_plane11,opts_vv00)
;    overlay(plot_olr11,	plot_vv11)

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    panel(0)  = plot_olr00
    panel(1)  = plot_olr01
    panel(2)  = plot_olr02
    panel(3)  = plot_olr03
    panel(4)  = plot_olr04
    panel(5)  = plot_olr05
    panel(6)  = plot_olr06
    panel(7)  = plot_olr07
    panel(8)  = plot_olr08
    panel(9)  = plot_olr09
    panel(10) = plot_olr10
    panel(11) = plot_olr11

;=============================================
; Add markers denoting TC position to plots
;=============================================

   ; Add NH tropical cyclone symbol to index
     tc = NhlNewMarker(wks, "p", 137, 0.0, 0.0, 1.0, 1.0, 0.0)

     mres0               = True
     mres0@gsMarkerIndex = tc          ; TC symbol
     mres0@gsMarkerSizeF = 20.0        ; Marker size
     mres0@gsMarkerColor = "red3"      ; Marker colour
     mres0@gsMarkerThicknessF = 5.0    ; Marker thickness (increase)

     markers(0,it)  = gsn_add_polymarker(wks,panel(0),\
                      centre(0,it-1,1),centre(0,it-1,0),mres0)
     markers(1,it)  = gsn_add_polymarker(wks,panel(1),\
                      centre(1,it-1,1),centre(1,it-1,0),mres0)
     markers(2,it)  = gsn_add_polymarker(wks,panel(2),\
                      centre(2,it-1,1),centre(2,it-1,0),mres0)
     markers(3,it)  = gsn_add_polymarker(wks,panel(3),\
                      centre(3,it-1,1),centre(3,it-1,0),mres0)
     markers(4,it)  = gsn_add_polymarker(wks,panel(4),\
                      centre(4,it-1,1),centre(4,it-1,0),mres0)
     markers(5,it)  = gsn_add_polymarker(wks,panel(5),\
                      centre(5,it-1,1),centre(5,it-1,0),mres0)
     markers(6,it)  = gsn_add_polymarker(wks,panel(6),\
                      centre(6,it-1,1),centre(6,it-1,0),mres0)
     markers(7,it)  = gsn_add_polymarker(wks,panel(7),\
                      centre(7,it-1,1),centre(7,it-1,0),mres0)
     markers(8,it)  = gsn_add_polymarker(wks,panel(8),\
                      centre(8,it-1,1),centre(8,it-1,0),mres0)
     markers(9,it)  = gsn_add_polymarker(wks,panel(9),\
                      centre(9,it-1,1),centre(9,it-1,0),mres0)
     markers(10,it) = gsn_add_polymarker(wks,panel(10),\
                      centre(10,it-1,1),centre(10,it-1,0),mres0)
     markers(11,it) = gsn_add_polymarker(wks,panel(11),\
                      centre(11,it-1,1),centre(11,it-1,0),mres0)
     delete(mres0)

    optsP                       = True
    optsP@gsnFrame              = False    ; Do not advance the frame
    optsP@gsnPanelLabelBar      = True     ; Turn off panel labelbar
    optsP@pmLabelBarWidthF      = 0.6
    optsP@pmLabelBarHeightF     = 0.15
    optsP@lbLabelFontHeightF    = 0.01
    optsP@lbPerimOn             = False
    optsP@lbLabelFont           = "Helvetica"
    optsP@lbTitleString         = "OLR (W m~S~-2~N~)"
    optsP@lbTitleFontHeightF    = 0.015
    optsP@lbTitleDirection      = "Across"
    optsP@lbTitlePosition	= "Bottom"
    optsP@txString              = "OLR: "+title_arr(it)
    optsP@gsnPanelFigureStrings = (/"a) em00", "b) em01", "c) em02", \
                                    "d) em03", "e) em04", "f) em05", \
			            "g) em06", "h) em07", "i) em08", \
				    "j) em09", "k) em10", "l) em11"/)
    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.98
    optsP@gsnPanelBottom        = 0.02
    optsP@amJust                = "TopLeft"
    optsP@gsnPanelFigureStringsFontHeightF = 0.012 ; Label size (default 0.01)

    gsn_panel(wks,panel,(/3,4/),optsP)             ; Draw as a single plot
    frame(wks)

    end do     ; End time loop

end 