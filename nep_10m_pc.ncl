; Script to compute the storm track for given TCs

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_rm.ncl"

; ncl 'dat="03T00"' int=1 plt=0 'opt="pdf"' 'ens0="em00"' dist=3.0 nep_10m_pc.ncl

begin

;========================
; Prelash (file input)
;========================

; Find 'pc' stream data to read in and plot 
  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0
  fili_c1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pc*.nc")
  fili_c      = diri+"/"+fili_c1

; Get dimension information from one of the files
  numINFO_0   = nc_times(fili_c)        ; Retrieve number of times
  numINPUT    = numINFO_0[0]            ; Number of times in 'pc' stream
  print("pc stream files contain "+numINPUT+" times")
  time_c      = numINFO_0[1]            ; Array of times from 'pc' stream file

; Calculate storm motion using built-in function ('st_rm')
; 'ens0'    = ensemble simulation
; 'diri'    = path to input files (see above)
; 'dist'    = size of box following storm (degrees)

  storm_rel = st_rm(dat,ens0,"$ar/text/",dist)
  u_cyc     = storm_rel[0]             ; Zonal wind
  v_cyc     = storm_rel[1]             ; Meridional wind
  vel_cyc   = storm_rel[2]             ; Vector wind
  lat_arr   = storm_rel[4]             ; Latitude array subset (following storm)
  lon_arr   = storm_rel[5]             ; Longitude array subset (following storm)
  centre    = storm_rel[6]             ; Storm track (position) array (lat,lon)

; Create arrays for data and output information
  llbox     = toint(dist*50)                   ; Calculate domain size (grid points)
  dsize     = (/llbox,llbox/)                  ; Domain size (grid points)

; Extract string for each ensemble member
  str1    = str_split(fili_c(0),"_")
  ens_str = str1(5)

; Create array to hold values of maximum 10-m (12.5 m) windspeed
  vec_box  = new(numINPUT,float)	       ; Maximum 10-m windspeed

; Set counter variable (for writing above diagnostics to text files)
  ct       = 0

;=================================================
; Read in best track information from text file
;=================================================

 ; Number of times in IBTrACS dataset for Nepartak
   dsn = 28

 ; Read in all latitude/longitude/data information
   lat_bt0  = asciiread("$sam/um/cp/nepartak/lat_nep.txt",(/dsn/),"float")
   lon_bt0  = asciiread("$sam/um/cp/nepartak/lon_nep.txt",(/dsn/),"float")
   dat_bt0  = asciiread("$sam/um/cp/nepartak/dat_nep.txt",(/dsn/),"string")

 ; From these arrays, extract earliest storm position
 ; and use this to constrain the initial position of the storm in the UM data
   if (dat .eq. "02T12") then  
    lat0     = lat_bt0(1)
    lon0     = lon_bt0(1)
   elseif (dat .eq. "03T00") then 
    lat0     = lat_bt0(3)
    lon0     = lon_bt0(3)
   elseif (dat	.eq. "03T12") then
    lat0     = lat_bt0(5)
    lon0     = lon_bt0(5)
   elseif (dat .eq. "04T00") then
    lat0     = lat_bt0(7)
    lon0     = lon_bt0(7)
   end if 

;=====================================================
; Create correct date strings for each output time
;=====================================================

   month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
                     "Jul","Aug","Sep","Oct","Nov","Dec"/)

   utc_date   = cd_calendar(time_c, 0)

   year       = tointeger(utc_date(:,0))
   month      = tointeger(utc_date(:,1))
   day        = tointeger(utc_date(:,2))
   hour       = tointeger(utc_date(:,3))
   minute     = tointeger(utc_date(:,4))
   second     = utc_date(:,5)

 ; Correct for errors in the code (round up value of hour when minutes = 59)
   do ct = 0, numINPUT-1
    if (minute(ct).gt.30) then
      hour(ct) = hour(ct)+1
    end if
   end do
   delete(ct)

 ; Create arrays to hold date strings
   date_str  = new(numINPUT,string)
   time_str  = new(numINPUT,string)
   time_arr  = new(numINPUT,string)
   title_arr = new(numINPUT,string)

   date_str = sprinti("%0.2iUTC ", hour) + \
              sprinti("%0.2i ", day) \
            + month_abbr(month)
   time_str = sprinti("%0.2i UTC ", hour) + \
              sprinti("%0.2i ", day) \
            + month_abbr(month)
   time_arr = sprinti("%0.2i", day) + \
              month_abbr(month) + \
              "_" + sprinti("%0.2iZ", hour)

   print("Working on ensemble member "+ens_str)

   setvalues NhlGetWorkspaceObjectId
     "wsMaximumSize" : 10000000000
   end setvalues

;==================================
; Get the variables we will need 
;==================================

    c   = addfile(fili_c,"r")          ; Read in 'pc' stream [u,v,vort]
    ct  = 0			       ; Counter variable
   
    if (dat .eq. "03T00" .and. ens0 .eq. "em07") then 
    numINPUT = numINPUT-1
    end if 

  do it = 1, numINPUT-1              ; Loop over times in file

    time  = c->t(it)                 ; Time in file
    times = dimsizes(time)

    title_arr(it) = "valid at "+time_str(it)+" (T+"+it+")"
    print("Working on "+time_str(it)+" (T+"+it+")")

  ; Get domain boundaries at each time from storm track (calculated off-line)
    t0    = lat_arr(ct,0)
    t1    = lat_arr(ct,1)
    n0    = lon_arr(ct,0)
    n1    = lon_arr(ct,1)

  ; Read in basic variables from 'pc' stream
    lon1C  = c->longitude_1({n0:n1})              ; longitude ['d1' grid points]
    lat1C  = c->latitude_1({t0:t1})               ; latitude  ['d0' grid points]
    lonC   = c->longitude({n0:n1})                ; longitude ['d1' grid points]
    latC   = c->latitude({t0:t1})                 ; latitude  ['d0' grid points]
    hybC   = c->hybrid_ht(:)                      ; 63 model ('theta') levels

    u_varname = "x-wind"
    v_varname = "y-wind"

  ; Read in horizontal wind components on model level 1 (12.5 m)
  ; Both variables below -- [latC | 150] * [lonC | 150]
    u  = c->$u_varname$(it,1,{t0:t1},{n0:n1})     ; Zonal wind
    v  = c->$v_varname$(it,1,{t0:t1},{n0:n1})     ; Meridional wind

    vec = sqrt( (u ^ 2) + (v ^ 2) )
    copy_VarCoords(u,vec)

;=====================================================================================
; Once we've correctly identified the storm, write MSLP and 10-m windspeed to files
;=====================================================================================

    rad = 1.0

    vec_box(it) = max(vec({centre(ct,0)-rad:centre(ct,0)+rad},\
    		          {centre(ct,1)-rad:centre(ct,1)+rad}))

;==================================================
; Plot data over whole domain (skip if required)
;==================================================

   if (plt .eq. 1) then 

  ; Output file location and type
    output = "$sam/nepartak/images/basic/track_pc_"+dat+"_"+ens_str+"_"+time_arr(it)
    wks    = gsn_open_wks(opt,output)

;==============================================================
; Load NIMROD colour table for precip duration data
;==============================================================

  ; Load colour table
    gsn_define_colormap(wks,"radar_new")

  ; Change the colour map around & make first colour transparent
    cmap = read_colormap_file("radar_new")
    cmap(0,3) = 0.0

  ; Create panel plot (MSLP and 10-m windspeed)
    panel = new(1,graphic)

;==========================
; Options for plotting
;==========================

  ; 10-m windspeed 
    opts_speed                              = True
    opts_speed@cnFillOn                     = True
    opts_speed@cnLineLabelInterval          = 2.0
    opts_speed@cnLineLabelFontHeightF       = 0.012
    opts_speed@cnLineLabelBackgroundColor   = "transparent"
    opts_speed@cnLineLabelPlacementMode     = "constant"
    opts_speed@cnLinesOn                    = False ; Contour lines off
    opts_speed@cnInfoLabelOn                = False ; Contour labels off
    opts_speed@cnLevelSelectionMode         = "ExplicitLevels"
    opts_speed@cnFillPalette                = cmap
    opts_speed@cnLevels                     = (/15.0, 20.0,\
                                                25.0, 30.0, 35.0, 40.0, \
                                                45.0, 50.0, 55.0, 60.0, \
                                                65.0/)
    opts_speed@cnFillColors                 = (/0,4,5,6,7,8,9,10,11,12,13,14/)

    opts_speed@gsnPaperOrientation          = "landscape"
    opts_speed@tiMainString                 = ""
    opts_speed@tiMainFontHeightF            = 0.0125
    opts_speed@gsnLeftString                = ""
    opts_speed@gsnRightString               = ""

    opts_speed@lbLabelBarOn                 = True
    opts_speed@lbBoxEndCapStyle             = "TriangleBothEnds"    ; Labelbar end shape
    opts_speed@lbLabelFontHeightF           = 0.0125                ; Labelbar font size
    opts_speed@lbLabelFont                  = "Helvetica"           ; Labelbar font
    opts_speed@lbTitleString                = "Wind speed (m s~S~-1~N~)"
    opts_speed@lbTitlePosition              = "Bottom"              ; Title position
    opts_speed@lbTitleFontHeightF           = 0.015                 ; Title font height
    opts_speed@pmLabelBarOrthogonalPosF     = 0.1                   ; Labelbar position
    opts_speed@lbTitleOffsetF               = 0                     ; Lb title up/down
    opts_speed@lbPerimOn                    = False                 ; Perimeter on/off
    opts_speed@gsnDraw                      = False                 ; Do not draw plot
    opts_speed@gsnFrame                     = False                 ; Do not adv. frame
    opts_speed@gsnAddCyclic                 = False

    opts_speed@mpLimitMode                  = "LatLon"
    opts_speed@mpMinLatF                    = lat(0)
    opts_speed@mpMinLonF                    = lon(0)
    opts_speed@mpMaxLatF                    = lat(dimsizes(lat)-1)
    opts_speed@mpMaxLonF                    = lon(dimsizes(lon)-1)

    opts_speed@tmXTOn                       = "False"    ; No tickmarks on top x-axis
    opts_speed@tmYROn                       = "False"    ; No tickmarks on right y-axis
    opts_speed@mpGridAndLimbOn              = True       ; Grid lines on/off
    opts_speed@mpGridLineDashPattern        = 2          ; Dash pattern
    opts_speed@mpGridLatSpacingF            = 2.0        ; Spacing (latitude)
    opts_speed@mpGridLonSpacingF            = 2.0        ; Spacing (longitude)

;=================
; Plot the data
;=================

    spd_plot  = gsn_csm_contour_map(wks,vec,opts_speed)		; 10-m wind speed
    panel(0)  = spd_plot

    marker    = new(1,graphic)

  ; Add NH tropical cyclone symbol to index
    tc0 = NhlNewMarker(wks, "p", 137, 0.0, 0.0, 1.0, 1.0, 0.0)

    mres0               = True
    mres0@gsMarkerIndex = tc0         ; TC symbol
    mres0@gsMarkerSizeF = 20.0        ; Marker size
    mres0@gsMarkerColor = "red3"      ; Marker colour
    mres0@gsMarkerThicknessF = 5.0    ; Marker thickness (increase)

    marker    = gsn_add_polymarker(wks,spd_plot,centre(ct,1),centre(ct,0),mres0)

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    optsP                       = True
    optsP@gsnFrame              = False             ; Do not advance the frame
    optsP@gsnPanelLabelBar      = False             ; Turn off panel labelbar
    optsP@txString              = "Nepartak (2-9 July 2016): "+title_arr(ct)
    optsP@gsnPanelFigureStrings = ""
    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.90
    optsP@gsnPanelBottom        = 0.10
    optsP@amJust                = "TopLeft"
    optsP@gsnPanelFigureStringsFontHeightF = 0.0125 ; Label size (default 0.01)

    gsn_panel(wks,panel,(/1,1/),optsP)              ; Draw as a single plot
    frame(wks)

   end if		; Only plot images if required 

  ; Tidy up 
    delete([/u,v,vec,lon1C,lat1C,lonC,latC,hybC/])	

    ct = ct + 1 	; Update counter variable (time)

   end do     ; End time loop (do it = 0, times-1, int)

;===================================
; Write out values to text files
;===================================

    diri_out = "$ar/text/cp_pc_"

  ; 'ens_str' is the number of the ensemble member (000 -- 011)

  ; Windspeed
    vec_out  = diri_out+dat+"_"+ens_str+"_10m.txt"
    asciiwrite(vec_out,vec_box)

end 