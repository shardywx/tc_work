; Script to plot windspeed following a cyclone for multiple ensemble members
; Run using [ncl int=1 dist=4.0 clev=1 'opt="pdf"' tc_ens_speed_panel_paper.ncl]
; Where 'int' is the interval in the time loop, 'dist' is the threshold for the cyclone
; 'opt' is the output file format ("pdf","x11",etc)
; 'clev' is the chosen pressure level for analysis [0 (950); 1 (850); 2 (700); 3 (500)]

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"

begin

;=====================================
; Calculate storm-relative velocity
;=====================================

   diri        = "/nfs/a37/earshar/ncl_scripts/nepartak/archer_march2018/"
   numSIMS     = 2				   ; Number of simulations
   numTIMES    = 120				   ; Number of times in each simulation
   centre      = new((/numSIMS,numTIMES,2/),float) ; Cyclone track in all simulations

 ; Read in values from ensemble members

   centre(0,:,0)  = asciiread(diri+"cp_lat_"+ens0+".txt",(/numTIMES,1/),"float")
   centre(0,:,1)  = asciiread(diri+"cp_lon_"+ens0+".txt",(/numTIMES,1/),"float")

   centre(1,:,0)  = asciiread(diri+"cp_lat_"+ens1+".txt",(/numTIMES,1/),"float")
   centre(1,:,1)  = asciiread(diri+"cp_lon_"+ens1+".txt",(/numTIMES,1/),"float")

;============================================================================
; Calculate distance between two points --> then calculate system velocity
;============================================================================

  print_clock("Starting system velocity calculations...")

; Define new arrays (cyclone position, speed)
  x_cyc   = new((/numSIMS,numTIMES/),float)
  y_cyc   = new((/numSIMS,numTIMES/),float)
  u_cyc   = new((/numSIMS,numTIMES/),float)
  v_cyc   = new((/numSIMS,numTIMES/),float)
  vel_cyc = new((/numSIMS,numTIMES/),float)

 do pt = 0, numSIMS-1

 ; Loop to calculate system velocity from 6-hourly cyclone position [lat,lon]
   do iang = 1, numTIMES-1
    x_cyc(pt,iang)   = gc_latlon(centre(pt,iang-1,0),centre(pt,iang-1,1),\
                       centre(pt,iang-1,0),centre(pt,iang,1),2,3)
    y_cyc(pt,iang)   = gc_latlon(centre(pt,iang-1,0),centre(pt,iang-1,1),\
                       centre(pt,iang,0),centre(pt,iang-1,1),2,3)
    u_cyc(pt,iang)   = x_cyc(pt,iang)/(3600)
    v_cyc(pt,iang)   = y_cyc(pt,iang)/(3600)
    vel_cyc(pt,iang) = sqrt(u_cyc(pt,iang)^2.0 + v_cyc(pt,iang)^2.0)
   end do

 ; Set initial system velocity [0] equal to that at the second time [1]
   x_cyc(pt,0)   = 0
   y_cyc(pt,0)   = 0
   u_cyc(pt,0)   = u_cyc(pt,1)
   v_cyc(pt,0)   = v_cyc(pt,1)
   vel_cyc(pt,0) = vel_cyc(pt,1)

 ; Add metadata to new arrays
   u_cyc@units         = "m s~S~-1~N~"
   u_cyc@description   = "Cyclone zonal velocity"

   vel_cyc@units       = "m s~S~-1~N~"
   vel_cyc@description = "Cyclone speed"

   u_cyc(pt,:) = -(u_cyc(pt,:)) ; Easterly motion is negative

  end do      ; End of loop over datasets (do pt = 0, numSIMS-1)

;================================================================================
; Create lat/lon arrays containing grid subsets at each time, following the TC
;================================================================================
  lat_arr = new((/numSIMS,numTIMES,2/),float)
  lon_arr = new((/numSIMS,numTIMES,2/),float)

; Lat/lon arrays for ensemble simulations
  lat_arr(0,:,0) = centre(0,:,0)-dist
  lat_arr(0,:,1) = centre(0,:,0)+dist
  lon_arr(0,:,0) = centre(0,:,1)-dist
  lon_arr(0,:,1) = centre(0,:,1)+dist

  lat_arr(1,:,0) = centre(1,:,0)-dist
  lat_arr(1,:,1) = centre(1,:,0)+dist
  lon_arr(1,:,0) = centre(1,:,1)-dist
  lon_arr(1,:,1) = centre(1,:,1)+dist

;======================================================================
; Read in horizontal windspeed; start loop over ensemble members 
;======================================================================

; Prelash (define arrays for loops below)
  ens_arr     = (/ens0,ens1/)
  numSIMS     = dimsizes(ens_arr)
  numTIMES    = 120
  time_arr    = new(numTIMES,string)
  title_arr   = new(numTIMES,string)
  vel_arr     = new((/numSIMS,numTIMES,200,200/),"float") ; Windspeed array

  date_str = new(numTIMES,string)
  time_str = new(numTIMES,string)
  out_str  = new(numTIMES,string)

  print_clock("Working on windspeed calculations!")

;====================================
; Start loop over ensemble members
;====================================

  do en = 0, dimsizes(ens_arr)-1

; List all files to read in and analyse
  diri        = "$ens/reg_netcdf/20160704T0000Z_ra1t_"
  input_pb    = diri+ens_arr(en)+"_pb.nc"

; Create extra space for large memory variables
  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 1000000000
  end setvalues

; Read in data
  f = addfile(input_pb,"r")

;==================================
; Get the variables we will need
;==================================

  do it = 23, 35         ; 00Z 5 July to 12Z 7 July (60 h)
;  do it = 0, numTIMES-1  ; ALL times 

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

   p_cp  = (/1,3,6,10,16/)          ; Pressure levels to access (950,850,700,500,200)

   pres  = f->p            ; Pressure (levels)
   plevs = dimsizes(pres)  ; Size of pressure level array
   pres@units = "hPa"

 ; Read in horizontal velocity, relative vorticity and geopotential height
   u       = f->u(it,p_cp,{t0:t1},{n0:n1})     ; Zonal velocity
   v       = f->v(it,p_cp,{t0:t1},{n0:n1})     ; Meridional velocity
   vort    = f->rvor(it,p_cp,{t0:t1},{n0:n1})  ; Relative vorticity
   z       = f->ht(it,p_cp,{t0:t1},{n0:n1})    ; Geopotential height

   z       = z/10
   z@units = "dam"              ; Convert to decametres
   vort    = vort * (10 ^ -6)   ; Convert to /s

   dy 	   = lat(1) - lat(0)    ; Grid spacing (latitude)
   dx 	   = lon(1) - lon(0)    ; Grid spacing (longitude)

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
   if (minute.gt.30) then
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
   time_arr(it) = out_str(it)
   title_arr(it) = "Valid at: "+time_str(it)

;======================================
; Calculate storm-relative velocity 
;======================================

 ; Storm-relative velocity
   u_rel      = u - u_cyc(en,it)
   v_rel      = v - v_cyc(en,it)

   d          = dimsizes(vort)       ; Size of all arrays [19 * 180 * 180]
   d0         = d(0)-1
   d1         = d(1)-1
   d2         = d(2)-1
   pii        = 3.14159265           ; Constant
   lat_rad    = lat * (pii/180)      ; Convert to radians for calculation below
   lon_rad    = lon * (pii/180)

 ; Add basic metadata
   copy_VarMeta(u,u_rel)
   copy_VarAtts(u,u_rel)
   copy_VarCoords(u,u_rel)

   copy_VarMeta(v,v_rel)
   copy_VarAtts(v,v_rel)
   copy_VarCoords(v,v_rel)

;====================================================================
; Find cyclone centre using both geopotential height and vorticity
;====================================================================

   ilev = clev
   vort_plane    = vort(ilev,:,:)
   geo_plane     = z(ilev,:,:)

 ; RELATIVE VORTICITY
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
   print("Cyclone centre (vort): "+lat_max950+" degrees N, "+lon_max950+" degrees E")

 ; GEOPOTENTIAL HEIGHT
   geo_min        = min(geo_plane)
   dims_h         = dimsizes(geo_plane)
   geo1d          = ndtooned(geo_plane)
   inds_h         = ind_resolve(minind(geo1d),dims_h)

   lat_min950 = t0 + (dy * inds_h(0,0))
   lon_min950 = n0 + (dx * inds_h(0,1))
   print("Cyclone centre (hgt): "+lat_min950+" degrees N, "+lon_min950+" degrees E")

   delete([/geo1d, inds_h/])

   ;========================================================
   ; EMPLOY A SAFETY NET IF WE IDENTIFY THE WRONG CYCLONE
   ;========================================================
     if (it.eq.0) then
      centre(en,it,0) = lat_min950
      centre(en,it,1) = lon_min950
     else
      lt0   = centre(en,it-1,0)         ; Previous TC latitude
      ln0   = centre(en,it-1,1)         ; Previous TC longitude
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

       geo_plane0 = z(ilev,{lt1:lt2},{ln1:ln2})          ; Smaller grid (1.0 x 1.0 deg)
       geo_min0   = min(geo_plane0)                      ; Find minimum on smaller grid
       dims_h0    = dimsizes(geo_plane0)                 ; Size of smaller grid
       geo_1d0    = ndtooned(geo_plane0)                 ; Create 1-D array
       inds_h0    = ind_resolve(minind(geo_1d0),dims_h0) ; Find index of minimum

       lat_min950 = lt1 + (dy * inds_h0(0,0))
       lon_min950 = ln1 + (dy * inds_h0(0,1))
       print("New centre (hgt): "+lat_min950+" degrees N, "+lon_min950+" degrees E")
       centre(en,it,0) = lat_min950
       centre(en,it,1) = lon_min950

       delete([/geo_plane0, geo_min0, dims_h0, geo_1d0, inds_h0/])

      else
       print("Cyclone centre (hgt): "+lat_min950+" degrees N, "+lon_min950+" degrees E")
       centre(en,it,0) = lat_min950
       centre(en,it,1) = lon_min950
      end if

     end if

   ; Tidy up before looping over pressure levels
     delete([/vort1d, inds, vort1d_sm, vort_smth, inds_sm/])
     delete([/vort_plane, geo_plane/])

;===================================================================
; Also read in vertical velocity [proxy for convective cells]
;===================================================================

     ilev = 1		; (0) 950 hPa, (1) 850 hPa, (2) 700 hPa, (3) 500 hPa, (4) 200 hPa
     p = pres(ilev)
     pstr = p+"hPa"

     u_plane = u_rel(ilev,:,:)
     v_plane = v_rel(ilev,:,:)
     spd_plane = sqrt((u_plane ^ 2) + (v_plane ^ 2))
     copy_VarCoords(v(ilev,:,:),spd_plane)
     copy_VarMeta(v(ilev,:,:),spd_plane)

     vel_arr(en,it,:,:) = spd_plane		; Input values into larger array      

     print("en = "+en)
     print("it = "+it)
     if(any(ismissing(spd_plane))) then
      print("Missing values in 'spd_plane'...")
     end if

  ; Tidy up (avoid dimension size errors)
    delete([/time,times,u,v,vort,z,lat,lon,lat1,lon1/])
    delete([/utc_date,year,month,day,hour,minute,second/])

  end do  ; End time loop (do it = 23, 107)
 
 end do	; End ensemble member loop (do en = 0, dimsizes(ens_arr)-1)

;===============================================
; Produce panel plots for each time interval
;===============================================

   print("Working on panel plots!")
   times = 120   		  

   do it = 24, times-1, int

   ; Read 850 hPa horizontal windspeed
     vel_plane00 = vel_arr(0,it,:,:)
     vel_plane01 = vel_arr(1,it,:,:)

     panel 	 = new(2,graphic)	; Panel plot with 6 images

;=================================
; Output file type and location
;=================================

     output = "$sam/nepartak/images/cart_coords/panel_speed"+p+"_paper_"+time_arr(it)
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

   ; 850 hPa horizontal windspeed
     opts_em00                              = True
     opts_em00@cnFillOn                     = True
     opts_em00@cnLineLabelInterval          = 2.0
     opts_em00@cnLineLabelFontHeightF       = 0.012
     opts_em00@cnLineLabelBackgroundColor   = "transparent"
     opts_em00@cnLineLabelPlacementMode     = "constant"
     opts_em00@cnLinesOn                    = False ; Contour lines off
     opts_em00@cnInfoLabelOn                = False ; Contour labels off
     opts_em00@cnLevelSelectionMode         = "ExplicitLevels"
     opts_em00@cnFillPalette		    = cmap
     opts_em00@cnLevels                     = (/5.0, 10.0, 15.0, 20.0,\
                                                25.0, 30.0, 35.0, 40.0, \
                                                45.0, 50.0, 55.0, 60.0, \
                                                65.0/)
     opts_em00@cnFillColors 		    = (/0,2,3,4,5,6,7,8,9,10,11,12,13,14/)

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

     r					    = 0.75	 ; Radius of plot

     opts_em00@gsnAddCyclic                 = False
     opts_em00@mpLimitMode		    = "Corners"
     opts_em00@mpLeftCornerLatF             = centre(0,it,0)-r
     opts_em00@mpLeftCornerLonF             = centre(0,it,1)-r
     opts_em00@mpRightCornerLatF            = centre(0,it,0)+r
     opts_em00@mpRightCornerLonF            = centre(0,it,1)+r

;     print_clock("Printing information for each ensemble member...")

   ; Try and customise the tickmarks to avoid size mismatching
     opts_em00@tmXBMode                     = "Manual"			    ; Set mode
     opts_em00@tmXBTickStartF               = toint(centre(0,it,1)-r) 	    ; LHS
     opts_em00@tmXBTickEndF                 = toint(centre(0,it,1)+r)-0.333 ; RHS
     opts_em00@tmXBTickSpacingF             = 0.5                    	    ; Spacing

     print_clock("Start longitude = "+opts_em00@tmXBTickStartF)
     print_clock("End longitude = "+opts_em00@tmXBTickEndF)

     opts_em00@tmYLMode                     = "Manual"               ; Set tickmark mode
     opts_em00@tmYLTickStartF               = centre(0,it,0)-r       ; Leftmost tickmark
     opts_em00@tmYLTickEndF                 = centre(0,it,0)-r-0.333 ; Rightmost tickmark
     opts_em00@tmYLTickSpacingF             = 0.5                    ; Tickmark spacing

     opts_em00@tmXBMajorLengthF             = 0.005      ; Change tickmark length
;     opts_em00@tmXBLabelFontHeightF         = 0.021
;     opts_em00@tmYLLabelFontHeightF         = 0.021
     opts_em00@tmXTOn                       = "False"    ; No tickmarks on top x-axis
     opts_em00@tmYROn                       = "False"    ; No tickmarks on right y-axis
     opts_em00@gsnMaximize                  = True       ; Maximise plot size

     opts_em00@gsnDraw		            = False      ; Do not draw the plot
     opts_em00@gsnFrame		            = False      ; Do not advance the frame

     opts_em01                              = opts_em00
     opts_em01@mpLeftCornerLatF             = centre(1,it,0)-r
     opts_em01@mpLeftCornerLonF             = centre(1,it,1)-r
     opts_em01@mpRightCornerLatF            = centre(1,it,0)+r
     opts_em01@mpRightCornerLonF            = centre(1,it,1)+r

   ; Try and customise the tickmarks to avoid size mismatching
     opts_em00@tmXBMode                     = "Manual"               ; Set tickmark mode
     opts_em00@tmXBTickStartF               = centre(1,it,1)-r       ; Leftmost tickmark
     opts_em00@tmXBTickEndF                 = centre(1,it,1)+r-0.333 ; Rightmost tickmark
     opts_em00@tmXBTickSpacingF             = 0.5                    ; Tickmark spacing

     opts_em00@tmYLMode                     = "Manual"               ; Set tickmark mode
     opts_em00@tmYLTickStartF               = centre(1,it,0)-r       ; Leftmost tickmark
     opts_em00@tmYLTickEndF                 = centre(1,it,0)-r-0.333 ; Rightmost tickmark
     opts_em00@tmYLTickSpacingF             = 0.5                    ; Tickmark spacing

;=============
; Plot data 
;=============

     plot_vel0     = gsn_csm_contour_map(wks,vel_plane00,opts_em00)      ; Panel 1
     plot_vel1     = gsn_csm_contour_map(wks,vel_plane01,opts_em01)      ; Panel 2

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    panel(0)   = plot_vel0
    panel(1)   = plot_vel1

    optsP                       = True
    optsP@gsnFrame              = False    ; Do not advance the frame
    optsP@gsnPanelLabelBar      = True     ; Turn off panel labelbar
    optsP@pmLabelBarWidthF      = 0.6
    optsP@pmLabelBarHeightF     = 0.15
    optsP@lbLabelFontHeightF    = 0.01
    optsP@lbPerimOn             = False
    optsP@lbLabelFont           = "Helvetica"
    optsP@lbTitleString         = p+" hPa horizontal windspeed (m s~S~-1~N~)"
    optsP@lbTitleFontHeightF    = 0.015
    optsP@lbTitleDirection      = "Across"
    optsP@lbTitlePosition       = "Bottom"
    optsP@txString              = ""
    optsP@gsnPanelFigureStrings = (/"a) "+ens0, "b) "+ens1/)
    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.98
    optsP@gsnPanelBottom        = 0.02
    optsP@amJust                = "TopLeft"
    optsP@gsnPanelFigureStringsFontHeightF = 0.0099 ; Reduce size of labels (0.01)

    gsn_panel(wks,panel,(/1,2/),optsP)             ; Draw as a single plot
    frame(wks)

   end do     ; End time loop

end 