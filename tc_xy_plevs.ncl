; Script to identify ERCs using Hövmoller plots of the three-dimensional wind field  

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_centre.ncl"
load "$sam/ncl_func/setup_cyl.ncl"
load "$sam/ncl_func/ring_mono.ncl"

; ncl opt=\"pdf\" dat=\"01T00\" ens0=\"em00\" dist=3.0 ts=0 tf=35 rmw=1 lay=1 ar=1.5 nr=31 nt=73
; ar0=1.0 vlev0=2 vlev1=5 wlev0=5 wlev1=7 tc_xy_plevs.ncl

; 'dat'    = initialisation time string ("02T12", "03T00", ...)
; 'ens0'   = ensemble simulation (em00-em11)
; 'dist'   = size of box following storm (degrees)
; 'opt'    = output file format ("pdf" or "x11")
; 'mlev'   = model level for plotting PV and vertical velocity (see L106-108)
; 'rmw'	   = read in RMW data on model (1) or pressure (0) levels
; 'calc'   = centre calculated on each level (0), calculated offline (1) or using 980 m vort (2)
; 'tclr'   = colour scheme for plotting tangential wind (0-2)
; 'pclr'   = colour scheme for plotting radial wind (0-2)
; 'wclr'   = colour for plotting vertical velocity (0-2)
; 'basic'  = plot u,v,w (1) or tendencies following N11 (0)
; 'lay'    = calculate layer average (1) or single level (0)
; 'mlev2'  = upper model level for calculating layer average (only valid when lay=1)
; 'ar'	   = distance in degrees of outer radius (1.0, 2.0, etc)
; 'nr'	   = number of radial circles between r = 0 and r = 'ar' (21, 41, 61, etc)
; 'ar0'	   = outer radius (degrees) for final plot [must be <= 'ar']
; 'cn0'    = method of storm tracking ["slp", "geo", "geo_sm", "vort"] - default is "geo_sm"

; 'mlev'   = 13 and 'mlev2' = 24 give an approximate 1-3 km layer average

; 'pc' stream --> 'it'
; 'centre'    --> 'it-1'
; 'pd' stream --> 'it-2'
 
begin

;==============================
; Prelash (reading in files)
;==============================

; Choose correct data stream(s) to read in and plot
; dat = '01T00', '01T12', '02T00', '02T12', '03T00', '04T00'
  diri        = "/nfs/a319/earshar/halong"
  fili_prefix = "201911"+dat+"00Z_SEA5_km4p4_ra2t_"+ens0
  fili_0      = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pb*.nc")
  fili        = diri+"/"+fili_0

; 10/02/2020 --> modify 'nc_times' to deal with multiple input files
  numINFO   = nc_times(fili)
  numINPUT  = numINFO[0]
  time      = numINFO[1]

; 10/02/2020 --> modify 'st_rm' or similar to calculate storm centre for other cases

; Create arrays for data and output information
  llbox     = toint(dist*50)                   ; Calculate domain size (grid points)
  dsize     = (/llbox,llbox/)                  ; Domain size (grid points)

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 10000000000
  end setvalues

;==========================================================================
; Create date/time string arrays for all times in file (use 'pc' stream)
;==========================================================================

  month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
                    "Jul","Aug","Sep","Oct","Nov","Dec"/)

  utc_date = cd_calendar(time, 0)

  year   = tointeger(utc_date(:,0))
  month  = tointeger(utc_date(:,1))
  day    = tointeger(utc_date(:,2))
  hour   = tointeger(utc_date(:,3))
  minute = tointeger(utc_date(:,4))
  second = utc_date(:,5)

; Correct for errors in the code (round up value of hour when minutes = 59)
  do ct = 0, numINPUT-1
   if (minute(ct).gt.30) then
     hour(ct) = hour(ct)+1
   end if
  end do
  delete(ct)

; Create arrays to hold strings for output plots
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

; Number of vertical (pressure) levels
  num_l      = 12

; Read in 'pb' stream data [u,v,w]
  b  = addfile(fili,"r")

; Counter variable (time)
  ct = 0                 

;============================
; Loop over times in file
;============================

  do it = ts, tf

    print("Working on time: "+time_str(it))
    title_arr(it) = "Valid at "+time_str(it)+" (T+"+it+")"

;========================================
; Read in variables on pressure levels
;========================================

  ; Grid subset (reduce computation time)
    t0  = 8.00
    n0  = 135.00

  ; Read in latitude and longitude arrays
    lon = b->longitude({n0::})
    lat = b->latitude({t0::})
 
  ; Pressure level array 
    p   = b->p

  ; Grid spacing 
    dy  = lat(1) - lat(0)
    dx  = lon(1) - lon(0)

 ; Read in horizontal wind components and vertical velocity (pressure coordinates)
 ; Both variables below -- [plevs | 12] * [latC | 150] * [lonC | 150]
   if (lay .eq. 1) then
    u0  = b->u(it,vlev0:vlev1,{t0::},{n0::})
    v0  = b->v(it,vlev0:vlev1,{t0::},{n0::})
    w0  = b->omega(it,wlev0:wlev1,{t0::},{n0::})
   else
    u0  = b->u(it,vlev0,{t0::},{n0::})      
    v0  = b->v(it,vlev0,{t0::},{n0::})      
    w0  = b->omega(it,wlev0,{t0::},{n0::})
   end if

 ; Calculate relative vorticity using centered finite differences  
   vort = uv2vr_cfd(u0, v0, lat, lon, 2)
   copy_VarCoords(u0, vort)
   vort@units = "s~S~-1~N~"
   vort@name  = "Relative vorticity on pressure levels"

;===========================================
; Analyse variables on chosen model level
;===========================================

    if (lay .eq. 1) then 
     print("Working on layer between " + p(vlev0)+ " and "+p(vlev1)+" hPa")
     u_plane       = dim_avg_n_Wrap(u0(:,:,:),0)
     v_plane       = dim_avg_n_Wrap(v0(:,:,:),0)
     w_plane       = dim_avg_n_Wrap(w0(:,:,:),0)
     vort_plane    = dim_avg_n_Wrap(vort(:,:,:),0)
    else
     print("Working on pressure level: "+vlev0+ " hPa")
     u_plane       = u0(:,:)
     v_plane       = v0(:,:)
     vort_plane    = vort(:,:)
     w_plane       = w0(:,:)
    end if 

  ; Horizontal wind vector
    spd_plane  = sqrt( u_plane ^ 2 + v_plane ^ 2)
    copy_VarCoords(u_plane, spd_plane)
    copy_VarMeta(u_plane, spd_plane)

    copy_VarCoords(u_plane, vort_plane)
    copy_VarMeta(u_plane, vort_plane)

;===============================
; Mid-script plotting routine 
;===============================

    out_loop = "$nep/halong/images/cxy_4panel_"+dat+"_"+ens0+"_"+time_arr(it)
    wks      = gsn_open_wks(opt, out_loop)

  ; Relative vorticity resources
    opts_rv                               = True
    opts_rv@cnFillOn                      = True
    opts_rv@cnLineLabelInterval           = 2.0
    opts_rv@cnLineLabelFontHeightF        = 0.012
    opts_rv@cnLineLabelBackgroundColor    = "transparent"
    opts_rv@cnLineLabelPlacementMode      = "constant"
    opts_rv@cnLinesOn                     = False
    opts_rv@cnInfoLabelOn                 = False
    opts_rv@cnLevelSelectionMode          = "ExplicitLevels"

  ; Load colour table and set contour levels
    opts_rv@cnFillPalette                 = "prcp_new"
    opts_rv@cnLevels                      = (/5., 10., 15., 20., 25., 30., \
                                                  40., 50., 60., 70., 80./)
    opts_rv@cnFillColors                  = (/-1,2,3,4,5,6,7,8,9,10,11,12/)

  ; No title or additional strings                                                     
    opts_rv@tiMainString                  = ""
    opts_rv@gsnLeftString                 = ""
    opts_rv@gsnRightString                = ""
    opts_rv@tiMainFontHeightF             = 0.0125

  ; Turn individual label bar on/off                                                   
    opts_rv@lbLabelBarOn                  = True
    opts_rv@lbBoxEndCapStyle              = "TriangleBothEnds"
    opts_rv@lbOrientation                 = "Vertical"

  ; Additional plotting resources                                                             
    opts_rv@mpDataBaseVersion             = "Ncarg4_1"     ; More recent database           
    opts_rv@mpDataSetName                 = "Earth..4"     ; High resolution                
    opts_rv@mpOutlineBoundarySets         = "National"     ; National borders               
    opts_rv@mpGeophysicalLineColor        = "black"        ; Colour borders black           
    opts_rv@mpGeophysicalLineThicknessF   = 1.0            ; Border line thickness          
    opts_rv@mpGridAndLimbOn               = True           ; Turn on lat/lon lines          
    opts_rv@mpGridSpacingF                = 2.0            ; Spacing between grid lines     
    opts_rv@pmTickMarkDisplayMode         = "Always"       ; Turn on map tickmarks          
    opts_rv@tmXBMajorLengthF              = 0.005          ; Change tickmark length         
    opts_rv@tmXBMinorOn                   = "False"        ; Turn off minor tm (x-axis)     
    opts_rv@tmYLMinorOn                   = "False"        ; Turn off minor tm (y-axis)     
    opts_rv@tmXTOn                        = "False"        ; Turn off tm (top x-axis)       
    opts_rv@tmYROn                        = "False"        ; Turn off tm (right y-axis)     
    opts_rv@gsnMaximize                   = "False"        ; Maximise plot size             
    opts_rv@gsnPaperOrientation           = "landscape"    ; change orientation             

    opts_rv@mpLimitMode                   = "Corners"
    opts_rv@mpLeftCornerLatF              = lat(0) + 0.5
    opts_rv@mpLeftCornerLonF              = lon(0) + 0.5
    opts_rv@mpRightCornerLatF             = lat(dimsizes(lat)-1) - 0.5
    opts_rv@mpRightCornerLonF             = lon(dimsizes(lon)-1) - 0.5

  ; Do not draw the plot or advance the frame                                          
    opts_rv@gsnDraw                       = False
    opts_rv@gsnFrame                      = False
    opts_rv@gsnAddCyclic                  = False

  ; Position panel on the page 
    opts_rv@vpXF                          = 0.05
    opts_rv@vpYF                          = 0.90
    opts_rv@vpWidthF                      = 0.35
    opts_rv@vpHeightF                     = 0.32

  ; Change units before plotting 
    vort_plane         = vort_plane * (10 ^ 4)
    vort_plane@units   = "10~S~-4~N~ s~S~-1~N~"

  ; Plot relative vorticity 
    plot_vort           = gsn_csm_contour_map(wks, spd_plane, opts_rv)

  ; Finally, draw the plot with everything overlaid (call 'draw' and 'frame')
    maximize_output(wks, True)

  ; Tidy up
    delete([/u_plane,v_plane,w_plane,vort_plane,spd_plane/])

  ; Advance counter variable (time)
    ct  = ct + 1

 ; Tidy up (avoid dimension size errors)
   delete([/lon,lat,u0,v0,w0,vort/])

  end do     ; End time loop (do it = ts, tf)

end 