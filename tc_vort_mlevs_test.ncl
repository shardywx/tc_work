; Script to calculate vorticity and geopotential height on model levels

; Run using:

; ncl 'opt="x11"' 'ens0="em10"' dist=5.0 ts0=30 tf0=36 r0=4.0 tc_vort_mlevs_test.ncl

; 'opt'   = output file format ("pdf" or "x11")
; 'ens0'  = ensemble member (em00 -- em11)
; 'dist'  = radius of box following storm (degrees)
; 'r0'	  = radius of final plot (degrees)
; 'ts0'   = starting time index (0–119)
; 'tf0'   = ending time index (0–119)

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"

begin

;==============================
; Prelash (reading in files)
;==============================

; Find 'pc' and 'pd' streams (model height level data)
;  diri        = "/nfs/a299/TCs/neptark/jul2018"
  diri        = "/nfs/a319/earshar"
  fili_prefix = "20160702T1200Z_NPTK_4p4_L80_ra1t_"+ens0

  fili_c1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pc*.nc")
  fili_c      = diri+"/"+fili_c1
  fili_d1     =	systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pd*.nc")
  fili_d      = diri+"/"+fili_d1

; Calculate number of times in each file using built-in function ('nc_times')
  numINFO_c   = nc_times(fili_c)
  numINFO_d   = nc_times(fili_d)
  numINPUT_c  = numINFO_c[0]		; Number of times in 'pc' stream
  time_c      = numINFO_c[1]		; Array of times from 'pc' stream file 
  numINPUT_d  = numINFO_d[0]		; Number of times in 'pd' stream

; Calculate storm motion using built-in function ('st_rm')
; Argument 1  = ensemble simulation
; Argument 2  = path to input files
; Argument 3  = size of box following storm (degrees)

; Call external function 'st_rm'
; Values extend from 12Z 2 July to 11Z 7 July (120 h)
  storm_rel = st_rm(ens0,"$ar/text/",dist)
  u_cyc     = storm_rel[0]		   ; Zonal wind 
  v_cyc     = storm_rel[1]		   ; Meridional wind
  vel_cyc   = storm_rel[2]		   ; Vector wind 
  numTIMES  = storm_rel[3]                 ; Number of times in file
  lat_arr   = storm_rel[4]                 ; Storm-following domain subset (lat)
  lon_arr   = storm_rel[5]                 ; Storm-following domain subset (lon)
  centre    = storm_rel[6]		   ; Storm centre at each hour (lat/lon)

;==========================================================================
; Create date/time string arrays for all times in file (use 'pb' stream)
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

; Fill arrays with date strings at each time interval 
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

  ; 12Z 2 Jul [00], 18Z 2 Jul [01], 00Z 3 Jul [02], 06Z 3 Jul [03], 12Z 3 Jul [04],
  ; 18Z 3 Jul [05], 00Z 4 Jul [06], 06Z 4 Jul [07], 12Z 4 Jul [08], 18Z 4 Jul [09],
  ; 00Z 5 Jul [10], 06Z 5 Jul [11], 12Z 5 Jul [12], 18Z 5 Jul [13], 00Z 6 Jul [14],
  ; 06Z 6 Jul [15], 12Z 6 Jul [16], 18Z 6 Jul [17], 00Z 7 Jul [18], 06Z 7 Jul [19],
  ; 12Z 7 Jul [20]

  ; ct = 0   [00], ct = 5  [01], ct = 11  [02], ct = 17  [03], ct = 23  [04],
  ; ct = 29  [05], ct = 35 [06], ct = 41  [07], ct = 47  [08], ct = 53  [09],
  ; ct = 59  [10], ct = 65 [11], ct = 71  [12], ct = 77  [13], ct = 83  [14],
  ; ct = 89  [15], ct = 95 [16], ct = 101 [17], ct = 107 [18], ct = 113 [19],
  ; ct = 119 [20]

    numTIMES    = 120                       ; Total number of times in all files [36]
    ct          = ts0                       ; Counter variable
    ct0		= 0			    ; Counter variable (independent of position)
    plevs	= 19			    ; Number of pressure levels
    llbox	= toint(dist*50)	    ; Size of lat/lon arrays (see line below)
    dsize	= (/llbox,llbox/)           ; Size of lat/lon arrays (depends on 'dist')
    d0		= dsize(0)
    d1		= dsize(1)

    setvalues NhlGetWorkspaceObjectId
      "wsMaximumSize" : 1000000000
    end setvalues

;===============================
; Start loop over input files 
;===============================

    c  = addfile(fili_c,"r")         	 ; Read in 'pc' stream [u,v]
    d  = addfile(fili_d,"r")         	 ; Read in 'pd' stream [w,p,t,theta,kmh,tau]

;==================================
; Get the variables we will need
;==================================

   do it = ts0, tf0		; Loop over times in file

    print("Working on time: "+time_str(it))
    itt = it + 1
    title_arr(it) = "Valid at "+time_str(it)+" (T+"+itt+")"

;=====================================================
; Read in variables from 'pc' stream (model levels)
;=====================================================

    t0    = lat_arr(it,0)
    t1    = lat_arr(it,1)
    n0    = lon_arr(it,0)
    n1    = lon_arr(it,1)
    print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

  ; Read in basic variables from 'pc' stream
    lon1C  = c->longitude_1({n0:n1})              ; longitude [xxx grid points]
    lat1C  = c->latitude_1({t0:t1})               ; latitude  [xxx grid points]
    lonC   = c->longitude({n0:n1})                ; longitude [xxx grid points]
    latC   = c->latitude({t0:t1})                 ; latitude  [xxx grid points]
    hybC   = c->hybrid_ht(:) 		          ; 63 hybrid levels (model 'theta' levs)
    num_m  = dimsizes(hybC)			  ; Number of model levels

    u_varname = "x-wind"
    v_varname = "y-wind"

  ; Read in horizontal wind components 
  ; Both variables below -- [hybC | 63] * [latC | 150] * [lonC | 150]
    u  = c->$u_varname$(it,:,{t0:t1},{n0:n1}) ; Zonal wind	
    v  = c->$v_varname$(it,:,{t0:t1},{n0:n1}) ; Meridional wind 

  ; Calculate relative vorticity using centered finite differences 
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
    f0@name	   = "Coriolis parameter"

    avo  = vort + f0
    copy_VarMeta(vort,avo)
    copy_VarAtts(vort,avo)
    copy_VarCoords(vort,avo)
    avo@description = "Absolute vorticity"
    avo@name	    = "Absolute vorticity on model levels"
    avo@long_name   = "Absolute vorticity"

;==================================================================
; Read in variables from 'pd' stream (model levels) 
; B/c of differences in file size, access variables using 'it+2'
;==================================================================

  ; Read in basic variables from 'pd' stream 
    lon1D   = d->longitude_1({n0:n1})                ; longitude [xxx grid points]
    lat1D   = d->latitude_1({t0:t1})                 ; latitude  [xxx grid points]
    lonD    = d->longitude({n0:n1})                  ; longitude [xxx grid points]
    latD    = d->latitude({t0:t1})                   ; latitude  [xxx grid points]
    hybD    = d->hybrid_ht(:)			     ; 63 hybrid levs (model 'rho' levs)

  ; Read in geopotential height -- [hybD | 63] * [latD | 150] * [lon1D | 151]
    z       = d->ht(it-2,:,{t0:t1},{n0:n1})	     ; Geopotential height (m)
    z       = z / 10				     ; Convert units before continuing
    z@units = "dam"				     ; Add metadata

    dx      = lonD(1) - lonD(0) ; Grid spacing (longitude)
    dy      = latD(1) - latD(0) ; Grid spacing (latitude)

;==========================
; Loop over model levels 
;==========================

    do ilev = 0, num_m-1, 3 ; Loop over all model levels

     hy = sprintf("%0.0f",hybC(ilev))
     print("Working on hybrid level " + ilev+ " ("+hy+" m AGL)")

     u_plane	   = u(ilev,:,:)
     v_plane	   = v(ilev,:,:)
     vort_plane    = vort(ilev,:,:)

     vort_max      = max(vort_plane)
     vort_smth     = smth9_Wrap(vort_plane, 0.5, 0.5, True)
     vort_max_smth = max(vort_smth)

     dims          = dimsizes(vort_plane)
     vort1d        = ndtooned(vort_plane)
     inds          = ind_resolve(maxind(vort1d),dims)
     vort1d_sm     = ndtooned(vort_smth)
     inds_sm       = ind_resolve(maxind(vort1d_sm),dims)

   ; Find latitude and longitude of storm at each time (and pressure level)
     lat_max 	   = t0 + (dy * inds(0,0))
     lon_max 	   = n0 + (dx * inds(0,1))

   ; Difference between 'lat_max'/'lon_max' and values read in using NCL function
     dlat 	   = abs(lat_max - centre(it,0))
     dlon	   = abs(lon_max - centre(it,1))

   ; If our calculations of latitude and longitude deviate too much from the values
   ; we read in using the NCL function ('st_rm'), then calculate a new centre
   ; using the geopotential height. Constrain the new calculation using the values
   ; from 'st_rm', to prevent spurious centres being identified and analysed. 

     if (dlat .gt. 0.75 .or. dlon .gt. 0.75) then 

      lt1 = centre(it,0) - 0.5
      lt2 = centre(it,0) + 0.5
      ln1 = centre(it,1) - 0.5
      ln2 = centre(it,1) + 0.5

      vort_plane0 = vort_plane({lt1:lt2},{ln1:ln2})       ; Smaller grid (1.0 x 1.0 deg)
      vort_min0   = min(vort_plane0)                      ; Find maximum on smaller grid
      dims_h0     = dimsizes(vort_plane0)                 ; Size of smaller grid
      vort_1d0    = ndtooned(vort_plane0)                 ; Create 1-D array
      inds_h0     = ind_resolve(maxind(vort_1d0),dims_h0) ; Find index of minimum

      lat_max 	 = lt1 + (dy * inds_h0(0,0))
      lon_max 	 = ln1 + (dx * inds_h0(0,1))
      print("New centre: "+lat_max+" degrees N, "+lon_max+" degrees E")

     end if 

  ; Plot data
    output = "vort_mlevs_test"
    wks = gsn_open_wks(opt,output)

    panel = new(1,graphic)
    gsn_define_colormap(wks,"bwr_new")

  ; Options for plotting vorticity and geopotential height 
    opts_vort                              = True
    opts_vort@cnFillOn                     = True
    opts_vort@cnLineLabelInterval          = 2.0
    opts_vort@cnLineLabelFontHeightF       = 0.012
    opts_vort@cnLineLabelBackgroundColor   = "transparent"
    opts_vort@cnLineLabelPlacementMode     = "constant"
    opts_vort@cnLinesOn                    = False ; Contour lines off
    opts_vort@cnInfoLabelOn                = False
    opts_vort@cnLevelSelectionMode         = "ExplicitLevels"

    opts_vort@cnLevels                     = (/-2., 2., 5., 10., 20., \
                                                30., 40., 60., 80., 100./)
    opts_vort@cnFillColors                 = (/13,0,19,21,23,25,24,26,27,28,29/)

    opts_vort@gsnPaperOrientation          = "landscape"
    opts_vort@tiMainString                 = ""
    opts_vort@tiMainFontHeightF            = 0.0125
    opts_vort@gsnLeftString                = ""
    opts_vort@gsnRightString               = ""
    opts_vort@tmXTOn                       = "False"
    opts_vort@tmYROn                       = "False"
    opts_vort@gsnMaximize                  = True
    opts_vort@lbLabelBarOn                 = False
    opts_vort@gsnDraw                      = False
    opts_vort@gsnFrame                     = False
    opts_vort@gsnAddCyclic		   = False

    opts_vort@mpLimitMode		   = "Corners"
    opts_vort@mpLeftCornerLatF		   = centre(it,0) - r0
    opts_vort@mpLeftCornerLonF		   = centre(it,1) - r0
    opts_vort@mpRightCornerLatF		   = centre(it,0) + r0
    opts_vort@mpRightCornerLonF		   = centre(it,1) + r0

  ; Winds
    opts_vec                               = True
    opts_vec@gsnDraw                       = False
    opts_vec@gsnFrame                      = False
    opts_vec@vcGlyphStyle                  = "LineArrow"   ; draw windbarbs
    opts_vec@vcRefMagnitudeF               = 20.0          ; ref vec = max vector
    opts_vec@vcRefLengthF                  = 0.04          ; length of reference vector
    opts_vec@vcMinDistanceF                = 0.05          ; distance between wind barbs
    opts_vec@vcLineArrowThicknessF         = 1.2           ; windbarb line thickness
    opts_vec@vcLineArrowColor              = "black"       ; windbarb colour
    opts_vec@gsnMaximize                   = True          ; Maximise plot size
    opts_vec@gsnRightString                = ""
    opts_vec@gsnLeftString                 = ""
    opts_vec@tiXAxisString                 = ""
    opts_vec@vcRefAnnoOrthogonalPosF       = -0.08


  ; Plot vorticity and winds 
    vort_plane = vort_plane * (10 ^ 4)
    plot_vort = gsn_csm_contour_map(wks,vort_plane,opts_vort)
    plot_vec  = gsn_csm_vector(wks,u_plane,v_plane,opts_vec)
    overlay(plot_vort,plot_vec)
    panel(0)  = plot_vort

;==================================================
; Add NH tropical cyclone symbol to panel plots
;==================================================

    markers    = new(1,graphic)

  ; Add NH tropical cyclone symbol to index
    tc0 = NhlNewMarker(wks, "p", 137, 0.0, 0.0, 1.0, 1.0, 0.0)

    mres0               = True
    mres0@gsMarkerIndex = tc0         ; TC symbol
    mres0@gsMarkerSizeF = 20.0        ; Marker size
    mres0@gsMarkerColor = "black"     ; Marker colour
    mres0@gsMarkerThicknessF = 5.0    ; Marker thickness (increase)

    markers = gsn_add_polymarker(wks,panel,centre(it,1),centre(it,0),mres0)

  ; Panel plot resources
    optsP                       = True                     ; Panel plot resources
    optsP@gsnFrame              = False                    ; Do not advance the frame
    optsP@gsnPanelLabelBar      = True                     ; Panel labelbar
    optsP@pmLabelBarWidthF      = 0.6                      ; Label bar width
    optsP@pmLabelBarHeightF     = 0.15                     ; Label bar height
    optsP@lbLabelFontHeightF    = 0.01                     ; Label bar font height
    optsP@lbLabelFont           = "Helvetica"              ; Choose font
    optsP@lbPerimOn             = False                    ; Perimeter on label bar
    optsP@lbTitleString         = "10 ~S~-4~N~ s~S~-1~N~"  ; Title string
    optsP@lbTitlePosition       = "Bottom"                 ; Title position
    optsP@lbTitleFontHeightF    = 0.015                    ; Title font height
    optsP@lbBoxEndCapStyle      = "TriangleBothEnds"       ; Label bar style
    optsP@pmLabelBarOrthogonalPosF = -0.002                ; Move lb up/down (h/l)
    optsP@lbTitleOffsetF        = 0                        ; x-axis title up/down (h/l)

    gsn_panel(wks,panel,(/1,1/),optsP)
    frame(wks)

  ; Tidy up
    delete([/vort_plane,u_plane,v_plane/])

    end do			    ; End pressure level loop (do ilev = 0, num_p-1)

  ; Tidy up before next iteration
    delete([/year,month,day,hour,minute,second,utc_date/])
    delete([/lonD,latD,lon1D,lat1D,time,times,z,vort,avo,f0,f1/])
    delete([/lonC,latC,lon1C,lat1C/])

   end do     ; End time loop (do it = 0, times-1)

end