; Script to identify vacillation cycles (VCs) in UM ensemble simulations

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

; ncl opt=\"x11\" km=\"8p8\" dat=\"03T00\" ens0=\"em00\" dist=1.0 r0=0.50 ts=0
; tf=3 tc=0 lay=1 vlev0=2 vlev1=5 wlev0=5 wlev1=7 paper=1 tc_xy_halong.ncl

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'ens0'  = ensemble simulation (em00-em11)
; 'dist'  = size of box following storm (degrees)
; 'opt'   = output file format ("pdf" or "x11")
; 'mlev0' = model level for plotting PV and vertical velocity (see L106-108)
; 'r0'    = radius for final plot (degrees)
; 'ts'	  = starting time index (0-116)
; 'tf'	  = ending time index (0-116)
; 'ti'	  = time interval (1,2,3,4,etc)
; 'full'  = produce 3-panel (1) or 2-panel (0) plot 
; 'tc'	  = annotate images with TC symbol (1) or leave as they were (0)
; 'lay'	  = calculate layer average (1) 
; 'mlev1' = upper model level for calculating layer average (only valid when lay=1)
; 'paper' = settings to produce plot for December 2019 version of VC paper (1)

begin

;==============================
; Prelash (reading in files)
;==============================

; Troubleshooting
  if (r0 .gt. dist) then 
   print("Choose a bigger domain ('dist') or smaller plot size ('r0')!")
   exit()
  end if  

; Choose correct data stream(s) to read in and plot  
; dat = '01T00', '01T12', '02T00', '02T12', '03T00', '04T00' 
  diri        = "/nfs/a319/earshar/halong"
  fili_prefix = "201911"+dat+"00Z_SEA5_km"+km+"_ra2t_"+ens0
  fili_0      = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pb*.nc")
  fili        = diri+"/"+fili_0

; 10/02/2020 --> modify 'nc_times' to deal with multiple input files 
  numINFO   = nc_times(fili)
  numINPUT  = numINFO[0]
  time      = numINFO[1]

; Create arrays for data and output information
  llbox	    = toint(dist*50)		       ; Calculate domain size (grid points)
  dsize	    = (/llbox,llbox/)		       ; Domain size (grid points)

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

  print("Working on ensemble member "+ens0)

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 10000000000
  end setvalues

; Number of vertical levels (pressure)
  num_l       = 12

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

  ; Read in variables on pressure levels

  ; Grid subset (reduce computation time)                                                
    t0  = 10.00
    t1  = 28.00
    n0  = 140.00

  ; Read in latitude and longitude arrays                                                
    lon = b->longitude({n0::})
    lat = b->latitude({t0:t1})

  ; Pressure level array                                                                 
    p   = b->p

  ; Grid spacing                                                                         
    dy  = lat(1) - lat(0)
    dx  = lon(1) - lon(0)

 ; Read in horizontal wind components and vertical velocity (pressure coordinates)       
 ; Both variables below -- [plevs | 12] * [latC | 150] * [lonC | 150]                    
   if (lay .eq. 1) then
    u0  = b->u(it,vlev0:vlev1,{t0:t1},{n0::})
    v0  = b->v(it,vlev0:vlev1,{t0:t1},{n0::})
    w0  = b->omega(it,wlev0:wlev1,{t0:t1},{n0::})
   else
    u0  = b->u(it,vlev0,{t0:t1},{n0::})
    v0  = b->v(it,vlev0,{t0:t1},{n0::})
    w0  = b->omega(it,wlev0,{t0:t1},{n0::})
   end if

 ; Calculate relative vorticity using centered finite differences                        
   vort = uv2vr_cfd(u0, v0, lat, lon, 2)
   copy_VarCoords(u0, vort)
   vort@units = "s~S~-1~N~"
   vort@name  = "Relative vorticity on pressure levels"

;==============================================
; Analyse variables on chosen pressure level
;==============================================

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

;=================
; Plot the data 
;=================

   ; Choose size of panel plot 
     panel = new(3,graphic)

   ; Set output file type and location
     out_loop = "$nep/halong/images/cxy_4panel_"+dat+"_"+ens0+"_"+km+"_"+time_arr(it)
     wks = gsn_open_wks(opt,out_loop)

   ; Load colour table
     gsn_define_colormap(wks,"prcp_new")

;==========================
; Options for plotting
;==========================

   ; Vorticity
     opts_rv                               = True
     opts_rv@cnFillOn                      = True
     opts_rv@cnLineLabelInterval           = 2.0
     opts_rv@cnLineLabelFontHeightF        = 0.012
     opts_rv@cnLineLabelBackgroundColor    = "transparent"
     opts_rv@cnLineLabelPlacementMode      = "constant"
     opts_rv@cnLinesOn                     = False
     opts_rv@cnInfoLabelOn                 = False
     opts_rv@cnLevelSelectionMode          = "ExplicitLevels"

     opts_rv@cnLevels                      = (/5., 10., 15., 20., 25., 30., \
                                               40., 50., 60., 70., 80./)
     opts_rv@cnFillColors                  = (/0,4,5,6,7,8,9,10,11,12,13,14/)

     opts_rv@tiMainString                  = ""
     opts_rv@tiMainFontHeightF             = 0.0125
     opts_rv@gsnLeftString                 = ""
     opts_rv@gsnRightString                = ""

   ; Labelbar on/off
     opts_rv@lbLabelBarOn                  = True
     opts_rv@lbOrientation                 = "horizontal"

   ; Additional lb resources (perimeter; label font; box end cap style)
     opts_rv@lbPerimOn                     = False
     opts_rv@lbLabelFont                   = "Helvetica"
     opts_rv@lbBoxEndCapStyle              = "TriangleBothEnds"

   ; Additional plotting resources
     opts_rv@mpDataBaseVersion             = "Ncarg4_1" ; More recent database
     opts_rv@mpDataSetName                 = "Earth..4" ; High resolution
     opts_rv@mpOutlineBoundarySets         = "National" ; National borders
     opts_rv@mpGeophysicalLineColor        = "black"    ; Colour borders black
     opts_rv@mpGeophysicalLineThicknessF   = 1.0        ; Border line thickness
     opts_rv@mpGridAndLimbOn               = False      ; Turn on lat/lon lines
     opts_rv@pmTickMarkDisplayMode         = "Always"   ; Turn on map tickmarks
     opts_rv@tmXBMajorLengthF              = 0.005      ; Change tickmark length
     opts_rv@tmXTOn                        = "False"    ; No tickmarks on top x-axis
     opts_rv@tmYROn                        = "False"    ; No tickmarks on right y-axis
     opts_rv@gsnMaximize                   = False      ; Maximise plot size

   ; Make lat/lon labels bigger (01/12/2019)
;     opts_rv@tmXBLabelFontHeightF          = 0.02
;     opts_rv@tmYLLabelFontHeightF          = 0.02

     opts_rv@gsnAddCyclic                  = False
     opts_rv@gsnDraw                       = False      ; Do not draw the plot
     opts_rv@gsnFrame                      = False      ; Do not advance the frame

   ; When reading in 'pd' stream data with 119 time intervals (PV)
   ; use the index one less than that in 'centre'
     opts_rv@mpLimitMode                   = "Corners"
     opts_rv@mpLeftCornerLatF              = lat(0) + 0.5
     opts_rv@mpLeftCornerLonF              = lon(0) + 0.5
     opts_rv@mpRightCornerLatF             = lat(dimsizes(lat)-1) - 0.5
     opts_rv@mpRightCornerLonF             = lon(dimsizes(lon)-1) - 0.5

   ; Position panel on the page 
     opts_rv@vpXF                          = 0.05
     opts_rv@vpYF                          = 0.80
     opts_rv@vpWidthF                      = 0.28
     opts_rv@vpHeightF                     = 0.28

   ; Change units before plotting 
     vort_plane  = vort_plane * (10 ^ 4)
     vort_plane@units = "10~S~-4~N~ s~S~-1~N~"

     copy_VarCoords(u_plane, vort_plane)

   ; Plot relative vorticity 
     plot_rv      = gsn_csm_contour_map(wks,vort_plane,opts_rv)

;==============================
; Also plot horizontal wind 
;==============================

   ; Horizontal wind resources 
     opts_vel                             = opts_rv
     delete([/opts_vel@cnLevels,opts_vel@cnFillColors/])
     opts_vel@cnFillPalette               = "radar_new"
     opts_vel@cnLevels                    = (/5.0, 10.0, 15.0, 20.0,\
                                                25.0, 30.0, 35.0, 40.0, \
                                                45.0, 50.0, 55.0, 60.0, \
                                                65.0/)
     opts_vel@cnFillColors                = (/-1,0,1,2,3,4,5,6,7,8,9,10,11,12/)

   ; Position panel on page 
     opts_vel@vpXF                        = 0.35
     opts_vel@vpYF                        = 0.80
     opts_vel@vpWidthF                    = 0.28
     opts_vel@vpHeightF                   = 0.28

   ; Smooth the wind field several times 
     spd_smth   = smth9_Wrap(spd_plane, 0.5, 0.5, False)
     spd_smth   = smth9_Wrap(spd_smth, 0.5, 0.5, False)
     spd_smth   = smth9_Wrap(spd_smth, 0.5, 0.5, False)     

   ; Plot horizontal wind 
     if (smth .eq. 1) then 
      plot_vel                            = gsn_csm_contour_map(wks,spd_smth,opts_vel)
     else
      plot_vel                            = gsn_csm_contour_map(wks,spd_plane,opts_vel)
     end if 

;============================================================
; Now plot vertical velocity using different colour table
;============================================================

   ; Vertical velocity resources
     opts_omg                              = opts_rv
     opts_omg@cnFillOn                     = True
     opts_omg@cnLineLabelInterval          = 2.0
     opts_omg@cnLineLabelFontHeightF       = 0.012
     opts_omg@cnLineLabelBackgroundColor   = "transparent"
     opts_omg@cnLineLabelPlacementMode     = "constant"

   ; Contour lines off; contour labels off 
     opts_omg@cnLinesOn                    = False
     opts_omg@cnInfoLabelOn                = False

   ; Colour table; contour levels
     opts_omg@cnLevelSelectionMode         = "ExplicitLevels"
     delete([/opts_omg@cnLevels,opts_omg@cnFillColors/])
     opts_omg@cnFillPalette               = "BlueDarkRed18"
     opts_omg@cnLevels                    = (/-3., -2., -1., -0.5, \
                                               0, 0.5, 1., 2., 3./)
     opts_omg@cnFillColors                = (/0,3,5,7,-1,-1,10,12,14,16/)

     opts_omg@tiMainString                 = ""
     opts_omg@tiMainFontHeightF            = 0.0125
     opts_omg@gsnLeftString                = ""
     opts_omg@gsnRightString               = ""

     opts_omg@lbTitleString                = ""

   ; Now position strain rate panel on page 
     opts_omg@vpXF                         = 0.65
     opts_omg@vpYF                         = 0.80
     opts_omg@vpWidthF                     = 0.28
     opts_omg@vpHeightF                    = 0.28

   ; Plot vertical velocity
     w_plane     = smth9_Wrap(w_plane, 0.5, 0.5, False)
     plot_omg    = gsn_csm_contour_map(wks,w_plane,opts_omg)

   ; Finally, draw the plot with everything overlaid (call 'draw' and 'frame')
     pres        = True 
     maximize_output(wks,pres)

  ; Tidy up before next iteration
    delete([/u_plane, v_plane, spd_plane, w_plane, vort_plane/])

  ; Tidy up (avoid dimension size errors)
   delete([/lon,lat,u0,v0,w0,vort/])

   ct = ct + 1

  end do      ; End time loop (do it = 0, numINPUT-1)

end