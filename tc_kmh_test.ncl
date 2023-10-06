; Script to plot wind stress and diffusivity (test for spurious patterns)

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"

; ncl dat=\"02T12\" ens0=\"em11\" opt=\"x11\" ts0=60 tf0=61 dist=3.0 cn0=\"slp\" tc_kmh_test.ncl

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'opt'   = output file format ("pdf" or "x11")
; 'ens0'  = ensemble member (em00 -- em11)
; 'dist'  = size of box following storm (degrees)
; 'ts0'   = starting time index (0-119)
; 'tf0'   = ending time index (0-119)
; 'cn0'   = offline storm centre: "slp", "vort", "geo", "geo_sm"

;==============================
; Prelash (reading in files)
;==============================

; Find model level files containing wind stress and diffusivity ('pd' stream)
  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0
  fili_c1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pc*.nc")
  fili_c      = diri+"/"+fili_c1
  fili_d1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pd*.nc")
  fili_d      = diri+"/"+fili_d1

; Calculate number of times in each file using built-in function ('nc_times')
  numINFO_d   = nc_times(fili_d)
  numINPUT_d  = numINFO_d[0]            ; Number of times in 'pd' stream
  time_d      = numINFO_d[1]            ; Array of times from 'pd' stream

  numINFO_c   = nc_times(fili_c)
  numINPUT_c  = numINFO_c[0]            ; Number of times in 'pd' stream
  time_c      = numINFO_c[1]            ; Array of times from 'pd' stream

; Calculate storm motion using built-in function ('st_rm')
; 'dat'     = initialisation time (02T12, 03T00, ...)
; 'ens0'    = ensemble simulation (em00, em01, ...)
; 'diri'    = path to input files (see above)
; 'dist'    = size of box following storm (degrees)
; 'mins'    = analyse 1-h (0) or 5-min (1) data

; Call external function 'st_rm'
  storm_rel = st_rm(dat, ens0, "$ar/text/", dist, 0)
  u_cyc     = storm_rel[0]                 ; Zonal wind
  v_cyc     = storm_rel[1]                 ; Meridional wind
  vel_cyc   = storm_rel[2]                 ; Vector wind
  numTIMES  = storm_rel[3]                 ; Number of times in file
  lat_arr   = storm_rel[4]                 ; Storm-following domain subset (lat)
  lon_arr   = storm_rel[5]                 ; Storm-following domain subset (lon)
  centre    = storm_rel[6]                 ; Storm centre at each hour (lat/lon)

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
  time_str  = new(numINPUT_c,string)
  time_arr  = new(numINPUT_c,string)

  time_str  = sprinti("%0.2i UTC ", hour) + \
              sprinti("%0.2i ", day) \
            + month_abbr(month)
  time_arr  = sprinti("%0.2i", day) + \
              month_abbr(month) + \
              "_" + sprinti("%0.2iZ", hour)

; Get dimensions of box following storm (dependent on user input 'dist')
  llbox     = toint(dist*50)      
  dsize     = (/llbox,llbox/)     
  d0        = dsize(0)
  d1        = dsize(1)

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 10000000000
  end setvalues

;=======================================
; Choose which storm track to read in
;=======================================

  if (cn0 .eq. "slp") then
   cn = 0
  elseif (cn0 .eq. "vort") then
   cn = 1
  elseif (cn0 .eq. "geo") then
   cn = 2
  elseif (cn0 .eq. "geo_sm") then
   cn = 3
  end if

;===================================
; Get the variables we will need
;===================================

  d  = addfile(fili_d,"r")              ; Read in 'pd' stream [w,p,t,theta,kmh,tau]

  do it = ts0, tf0			; Loop over desired times 

   print("Working on time: "+time_str(it))

   t0    = lat_arr(cn,it,0)
   t1    = lat_arr(cn,it,1)
   n0    = lon_arr(cn,it,0)
   n1    = lon_arr(cn,it,1)   
   print_clock("Domain subset: "+t0+" to "+t1+" deg N; "+n0+" to "+n1+" deg E.")

 ; Read in basic variables from 'pd' stream
   lon1D  = d->longitude_1({n0:n1})                ; longitude [xxx grid points]
   lat1D  = d->latitude_1({t0:t1})                 ; latitude  [xxx grid points]
   lonD   = d->longitude({n0:n1})                  ; longitude [xxx grid points]
   latD   = d->latitude({t0:t1})                   ; latitude  [xxx grid points]
   hybD   = d->hybrid_ht(:)                        ; 63 model ('rho') levels

 ; Read in diffusivity and wind stress
 ; 'str_xz' -- [hybD | 15] * [latD | 150]  * [lonD | 151]
 ; 'str_yz' -- [hybD | 15] * [lat1D | 150] * [lon1D | 151]
 ; 'kmh'    -- [hybD | 15] * [latD | 150]  * [lon1D | 151]
   str_xz = d->taux(it-1,0:12,{t0:t1},{n0:n1})        ; Wind stress, x-comp (N m/2)
   str_yz = d->tauy(it-1,0:12,{t0:t1},{n0:n1})        ; Wind stress, y-comp (N m/2)
   kmh    = d->unspecified(it-1,0:12,{t0:t1},{n0:n1}) ; Horizontal diffusivity (m^2 s-1)
   dstr	  = dimsizes(str_xz)

 ; Loop over model height levels 
   do ilev = 2, dstr(0)-1, 4

    hy = sprintf("%0.0f",hybD(ilev))
    print("Working on model level " + ilev+ " ("+hy+" m AGL)")

    ip1 = min((/ilev+1,dstr(0)-1/))
    im1 = max((/ilev-1,0/))
 
    kmh_plane  = kmh(ilev,:,:)
    strx_plane = str_xz(ilev,:,:)
    stry_plane = str_yz(ilev,:,:)

  ; Calculate levels above/below for derivatives

  ; x-component of wind stress
    str_xz0    = str_xz(im1,:,:)
    str_xz1    = str_xz(ip1,:,:)

  ; y-component of wind stress
    str_yz0    = str_yz(im1,:,:)
    str_yz1    = str_yz(ip1,:,:)

  ; Calculate derivatives 
    dtaux_dz   = (str_xz1 - str_xz0) / ( hybD(ip1) - hybD(im1) )
    dtauy_dz   = (str_yz1 - str_yz0) / ( hybD(ip1) - hybD(im1) )

  ; Plot data
    output     = "$sam/nepartak/images/cart_coords/str_xy_test_"+hy+"m_"+cn0
    wks        = gsn_open_wks(opt,output)
   
    panel      = new(3,graphic)
    gsn_define_colormap(wks,"BlueDarkRed18")

  ; Diffusivity resources 
    opts_mbgt                             = True
    opts_mbgt@cnFillOn                    = True
    opts_mbgt@cnLineLabelInterval         = 2.0
    opts_mbgt@cnLineLabelFontHeightF      = 0.012
    opts_mbgt@cnLineLabelBackgroundColor  = "transparent"
    opts_mbgt@cnLineLabelPlacementMode    = "constant"
    opts_mbgt@cnLinesOn                   = False ; Contour lines off
    opts_mbgt@cnInfoLabelOn               = False
    opts_mbgt@cnLevelSelectionMode        = "ExplicitLevels"

    opts_mbgt@cnFillColors		 = (/0,12,13,14,15,16,17,19/)
    opts_mbgt@cnLevels			 = (/0., 5., 10., 20., 30., 50., 100./)

    opts_mbgt@gsnPaperOrientation         = "landscape"
    opts_mbgt@tiMainString                = ""
    opts_mbgt@tiMainFontHeightF           = 0.0125
    opts_mbgt@gsnLeftString               = ""
    opts_mbgt@gsnRightString              = ""
    opts_mbgt@gsnMaximize                 = True
    opts_mbgt@lbLabelBarOn                = True
    opts_mbgt@lbLabelFontHeightF          = 0.015
    opts_mbgt@lbPerimOn             	  = False
    opts_mbgt@lbLabelFont           	  = "Helvetica"
    opts_mbgt@lbTitleString               = hy+" m diffusivity ("+kmh@units+")"
    opts_mbgt@lbTitlePosition             = "Bottom"              ; Title position
    opts_mbgt@lbTitleFontHeightF          = 0.025                 ; Title font height
    opts_mbgt@pmLabelBarOrthogonalPosF    = 0.1                   ; Labelbar position
    opts_mbgt@lbTitleOffsetF              = 0                     ; Lb title up/down

    opts_mbgt@lbBoxEndCapStyle      	  = "TriangleBothEnds"    ; Labelbar end shape 
    opts_mbgt@pmTickMarkDisplayMode       = "Always"   ; Turn on map tickmarks
    opts_mbgt@tmXBMajorLengthF            = 0.005      ; Change tickmark length
    opts_mbgt@tmXBLabelFontHeightF        = 0.021
    opts_mbgt@tmYLLabelFontHeightF        = 0.021
    opts_mbgt@tmXTOn                      = "False" ; Turn off top x-axis TM
    opts_mbgt@tmYROn                      = "False" ; Turn off right y-axis TM

  ; Set up labels every 1º in latitude and longitude
    opts_mbgt@tmXBMode			  = "Explicit"
    opts_mbgt@tmXBValues		  = fspan(100.0, 150.0, 26)
    xb_labels0                            = sprintf("%0.1f~S~o~N~E",opts_mbgt@tmXBValues)
    opts_mbgt@tmXBLabels                  = xb_labels0

    opts_mbgt@tmYLValues                  = fspan(0.0, 30.0, 16)
    yl_labels0                            = sprintf("%0.1f~S~o~N~N",opts_mbgt@tmYLValues)
    opts_mbgt@tmYLLabels                  = yl_labels0

    opts_mbgt@gsnDraw                     = False
    opts_mbgt@gsnFrame                    = False
    opts_mbgt@gsnAddCyclic		  = False

  ; Wind stress resources
    opts_taux				  = opts_mbgt
    delete([/opts_taux@cnFillColors,opts_taux@cnLevels/])
    opts_taux@cnFillColors		  = (/5,6,7,8,9,10,0,0,12,13,14,15,16,17/)
    opts_taux@cnLevels			  = (/-2.0, -1.0, -0.5, -0.3, -0.2, -0.1, \
   					      0.0, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0/)
    delete(opts_taux@lbTitleString)
    opts_taux@lbTitleString               = hy+" m wind stress x-comp ("+str_xz@units+")"

    opts_tauy				  = opts_taux
    delete(opts_tauy@lbTitleString)
    opts_tauy@lbTitleString               = hy+" m wind stress y-comp ("+str_yz@units+")"

  ; Plot panels side by side
    panel(0) = gsn_csm_contour(wks,kmh_plane,opts_mbgt)
    panel(1) = gsn_csm_contour(wks,strx_plane,opts_taux)
    panel(2) = gsn_csm_contour(wks,stry_plane,opts_tauy)

  ; Panel plot resources
    optsP                       = True                     ; Panel plot resources
    optsP@gsnFrame              = False                    ; Do not advance the frame
    optsP@gsnPanelLabelBar      = False                    ; Panel labelbar
    optsP@pmLabelBarWidthF      = 0.6                      ; Label bar width
    optsP@pmLabelBarHeightF     = 0.15                     ; Label bar height
    optsP@lbLabelFontHeightF    = 0.01                     ; Label bar font height
    optsP@lbLabelFont           = "Helvetica"              ; Choose font
    optsP@lbPerimOn             = False                    ; Perimeter on label bar
    optsP@lbTitleString         = "m s~S~-2~N~"            ; Title string
    optsP@lbTitlePosition       = "Bottom"                 ; Title position
    optsP@lbTitleFontHeightF    = 0.015                    ; Title font height
    optsP@lbBoxEndCapStyle      = "TriangleBothEnds"       ; Label bar style
    optsP@pmLabelBarOrthogonalPosF = -0.002                ; Move lb up/down (h/l)
    optsP@lbTitleOffsetF        = 0                        ; x-axis title up/down (h/l)
    optsP@gsnPanelFigureStrings = (/"a) K_mh", \
                                    "b) tau_x", \
                                    "c) tau_y"/)
 
    gsn_panel(wks,panel,(/1,3/),optsP)
    frame(wks)

   end do		; End model level loop (do ilev = ...)

  end do		; End time loop (do it = ...)