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

; ncl opt=\"x11\" km=\"8p8\" dat=\"04T00\" ens0=\"em07\" dist=2.0 ts=2 tf=30 rmw=1 lay=1 ar=1.5 
; nr=31 nt=73 ar0=1.0 plt=1 vlev0=2 vlev1=2 wlev0=5 wlev1=5 tc_hov_halong.ncl

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

begin

;==============================
; Prelash (reading in files)
;==============================

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

; 10/02/2020 --> modify 'st_rm' or similar to calculate storm centre for other cases

; Create arrays for data and output information
  llbox     = toint(dist*50)                   ; Calculate domain size (grid points)
  dsize     = (/llbox,llbox/)                  ; Domain size (grid points)

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 10000000000
  end setvalues

; Create array to hold (estimated) storm centre position at each time
  centre    = new( (/numINPUT,2/), "float")

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

;========================================
; Define cylindrical coordinate arrays
;========================================

; ar = distance (degrees) of outer radial circle (1.0)
; nr = total number of radial circles (21)

  pii        = 3.14159265
  radii      = fspan(0.0,ar,nr)      ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,nt)   ; Array of azimuth angles (0, 5, 10, etc)
  thetas_rad = (pii/180)*thetas      ; Degrees to radians
  rad_size   = dimsizes(radii)-1     ; Index for later calculation

; Define 1D arrays 
  vtan_azi   = new((/dimsizes(radii)/),float) ; Tangential wind
  vrad_azi   = new((/dimsizes(radii)/),float) ; Radial wind
  spd_azi    = new((/dimsizes(radii)/),float) ; Vector wind
  vvel_azi   = new((/dimsizes(radii)/),float) ; Vertical velocity
  vort_azi   = new((/dimsizes(radii)/),float) ; Relative vorticity

; Define 2D arrays
  vtan_plot0   = new((/dimsizes(radii),numINPUT/),float)
  vtan_plot0!0 = "rad"
  vtan_plot0!1 = "time"

  vrad_plot0   = new((/dimsizes(radii),numINPUT/),float)
  vrad_plot0!0 = "rad"
  vrad_plot0!1 = "time"

  vvel_plot0   = new((/dimsizes(radii),numINPUT/),float)
  vvel_plot0!0 = "rad"
  vvel_plot0!1 = "time"

  vort_plot0   = new((/dimsizes(radii),numINPUT/),float)
  vort_plot0!0 = "rad"
  vort_plot0!1 = "time"

  spd_plot0    = new((/dimsizes(radii),numINPUT/),float)
  spd_plot0!0  = "rad"
  spd_plot0!1  = "time"

; Add metadata for azimuthally averaged variables
  vtan_azi!0            = "rad"
  vtan_azi&rad          = radii
  vtan_azi@description  = "Azimuthally averaged tangential wind"
  vtan_azi@units        = "m s~S~-1~N~"

  vrad_azi!0            = "rad"
  vrad_azi&rad          = radii
  vrad_azi@description  = "Azimuthally averaged radial wind"
  vrad_azi@units        = "m s~S~-1~N~"

  spd_azi!0             = "rad"
  spd_azi&rad           = radii
  spd_azi@description   = "Azimuthally averaged vector wind"
  spd_azi@units         = "m s~S~-1~N~"

  vvel_azi!0            = "rad"
  vvel_azi&rad          = radii
  vvel_azi@description  = "Azimuthally averaged vertical velocity"
  vvel_azi@units        = "m s~S~-1~N~"

  vort_azi!0            = "rad"
  vort_azi&rad          = radii
  vort_azi@description  = "Azimuthally averaged relative vorticity"
  vort_azi@units        = "s~S~-1~N~"

; Create time array for later plotting 
  time_int   = new(numINPUT,"integer")
  time_int   = ispan(0,numINPUT-1,1)
  time_int!0 = "time"

;===============================
; Start loop over input files
;===============================

; Number of vertical (pressure) levels
  num_l      = 12

; Arrays to hold all 'xpos' and 'ypos' position arrays
  if (lay .eq. 1) then 
   xpos_all  = new((/dimsizes(radii),dimsizes(thetas),numINPUT/),float)
   ypos_all  = new((/dimsizes(radii),dimsizes(thetas),numINPUT/),float)
   xcen_all  = new((/dimsizes(radii),dimsizes(thetas),numINPUT/),float)
   ycen_all  = new((/dimsizes(radii),dimsizes(thetas),numINPUT/),float)
  else
   xpos_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),numINPUT/),float)
   ypos_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),numINPUT/),float)
   xcen_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),numINPUT/),float)
   ycen_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),numINPUT/),float)
  end if 

;==================================
; Get the variables we will need
;==================================

; Read in 'pb' stream data [u,v,w]
  b  = addfile(fili,"r")

; Counter variable (time)
  ct = 0                 

; Loop over times in file
  do it = ts, tf

    print("Working on time: "+time_str(it)+" (it = "+it+")")
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

 ; EDIT --> uncomment these lines when script is finished
 ; Calculate storm-relative winds
;   u = u0 - u_cyc(cn,it-1)
;   v = v0 - v_cyc(cn,it-1)

 ; Add metadata from original horizontal wind arrays
;   copy_VarCoords(u0,u)
;   copy_VarAtts(u0,u)
;   copy_VarCoords(v0,v)
;   copy_VarAtts(v0,v)

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

;========================================================================================
; Calculate smoothed wind field (proxy for sea level pressure or geopotential height)
;========================================================================================

    spd_smth   = smth9_Wrap(spd_plane, 0.5, 0.5, False)
    spd_smth   = smth9_Wrap(spd_smth, 0.5, 0.5, False)
    spd_smth   = smth9_Wrap(spd_smth, 0.5, 0.5, False)

;========================================================
; Manually find storm centre using relative vorticity
;========================================================

   ; Calculate maximum in smoothed 2D array
     vort_smth = smth9_Wrap(vort_plane, 0.5, 0.5, False)
     vort_smth = smth9_Wrap(vort_smth, 0.5, 0.5, False)
     vort_max0 = max(vort_smth)
     dims_v0   = dimsizes(vort_smth)

   ; Create 1D from 2D array
     vort_1d0  = ndtooned(vort_smth)

   ; Find index of maximum --> (0,0) = latitude, (0,1) = longitude
     inds_v0   = ind_resolve(maxind(vort_1d0),dims_v0)

   ; Retrieve lat/lon information from grid subset 
     ltN       = vort_plane&latitude
     lnN       = vort_plane&longitude

   ; Find lat/lon of vorticity maximum
     lat_max0  = ltN(0) + (dy * inds_v0(0,0) )
     lon_max0  = lnN(0) + (dx * inds_v0(0,1) )
     delete([/ltN, lnN, vort_1d0/])

   ; Calculate distance between current and previous storm centre
     dlat      = abs(lat_max0 - centre(it-1,0) )
     dlon      = abs(lon_max0 - centre(it-1,1) )

   ; If storm centre is further than a threshold distance from the last position, throw value out 
     if (ct .gt. 0) then 
      if (dlat .gt. 1.0 .or. dlon .gt. 1.0) then 

       print("Using previous storm position for current calculation...")
       print("Current storm position: "+lat_max0+","+lon_max0)
       print("Previous storm position: "+centre(it-1,0)+""+centre(it-1,1))

       lat_max0 = centre(it-1,0)
       lon_max0 = centre(it-1,1)

     ; Use larger search area for minimum windspeed if using previous storm position 
       lt0      = lat_max0 - 0.75
       lt1      = lat_max0 + 0.75
       ln0      = lon_max0 - 0.75
       ln1      = lon_max0 + 0.75

      else

     ; Find the minimum wind speed within a threshold distance of the vorticity maximum
       lt0      = lat_max0 - 0.25
       lt1      = lat_max0 + 0.25
       ln0      = lon_max0 - 0.25
       ln1      = lon_max0 + 0.25

      end if 

     else

     ; Find the minimum wind speed within a threshold distance of the vorticity maximum
       lt0      = lat_max0 - 0.25
       lt1      = lat_max0 + 0.25
       ln0      = lon_max0 - 0.25
       ln1      = lon_max0 + 0.25

     end if 

;================================================
; Refine this method using minimum wind speed
;================================================ 

   ; Calculate minimum in 2D array
     spd_pl_n  = spd_plane({lt0:lt1},{ln0:ln1})
     spd_min0  = min(spd_pl_n)
     dims_u0   = dimsizes(spd_pl_n)

   ; Create 1D from 2D array
     spd_1d0   = ndtooned(spd_pl_n)

   ; Find index of minimum
     inds_u0   = ind_resolve(minind(spd_1d0),dims_u0)

   ; Retrieve lat/lon information from grid subset 
     ltN       = spd_pl_n&latitude
     lnN       = spd_pl_n&longitude

     lat_max   = ltN(0) + (dy * inds_u0(0,0) )
     lon_max   = lnN(0) + (dx * inds_u0(0,1) )
     delete([/ltN, lnN, spd_pl_n, spd_1d0/])

   ; Input values into 'centre' array
     centre(it,0) = lat_max
     centre(it,1) = lon_max

;============================================
; Switch to cylindrical coordinate system
;============================================

   ; Create cylindrical grid centred on storm (lat_max, lon_max)
     do irad = 0, dimsizes(radii)-1
      r = radii(irad)
      do iang = 0, dimsizes(thetas)-1
       theta = thetas_rad(iang)
       xpos_all(irad,iang,ct) = lon_max + r * cos(theta)
       ypos_all(irad,iang,ct) = lat_max + r * sin(theta)
       xcen_all(irad,iang,ct) = lon_max
       ycen_all(irad,iang,ct) = lat_max
      end do
     end do

   ; Create new arrays for loop below 
     u_int    = new((/dimsizes(radii),dimsizes(thetas)/),float)
     v_int    = new((/dimsizes(radii),dimsizes(thetas)/),float)
     spd_int  = new((/dimsizes(radii),dimsizes(thetas)/),float)
     vort_int = new((/dimsizes(radii),dimsizes(thetas)/),float)
     w_int    = new((/dimsizes(radii),dimsizes(thetas)/),float)
     vr_int   = new((/dimsizes(radii),dimsizes(thetas)/),float)
     vt_int   = new((/dimsizes(radii),dimsizes(thetas)/),float)

   ; Linearly interpolate onto cylindrical grid
     do iang = 0, dimsizes(thetas_rad)-1
      theta = thetas_rad(iang)

    ; Horizontal wind components, relative vorticity and vertical velocity
      u_int(:,iang)    = linint2_points(lon,lat,u_plane,False,\
                                        xpos_all(:,iang,ct),ypos_all(:,iang,ct),0)
      v_int(:,iang)    = linint2_points(lon,lat,v_plane,False,\
                                        xpos_all(:,iang,ct),ypos_all(:,iang,ct),0)
      spd_int(:,iang)  = linint2_points(lon,lat,spd_plane,False,\
                                        xpos_all(:,iang,ct),ypos_all(:,iang,ct),0)
      vort_int(:,iang) = linint2_points(lon,lat,vort_plane,False,\
                                        xpos_all(:,iang,ct),ypos_all(:,iang,ct),0)      
      w_int(:,iang)    = linint2_points(lon,lat,w_plane,False,\
                                        xpos_all(:,iang,ct),ypos_all(:,iang,ct),0)

    ; Radial and tangential wind components         
      vr_int(:,iang)   = u_int(:,iang) * cos(theta) + v_int(:,iang) * sin(theta)
      vt_int(:,iang)   = -u_int(:,iang) * sin(theta) + v_int(:,iang) * cos(theta)

     end do

;===============================
; Mid-script plotting routine 
;===============================

     if (plt .eq. 1) then 

        out_loop = "$nep/halong/images/cxy_4panel_"+dat+"_"+ens0+"_"+time_arr(it)
        wks      = gsn_open_wks(opt, out_loop)
        gsn_define_colormap(wks,"prcp_new")

      ; Relative vorticity resources
        opts_rv                               = True
        opts_rv@cnFillOn                      = True
        opts_rv@sfXArray                      = xpos_all(:,:,ct)
        opts_rv@sfYArray                      = ypos_all(:,:,ct)
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
        opts_rv@gsnMaximize                   = False 

      ; No title or additional strings                                                     
        opts_rv@tiMainString                  = ""
        opts_rv@gsnLeftString                 = ""
        opts_rv@gsnRightString                = ""

      ; Turn individual label bar on/off                                                   
        opts_rv@lbLabelBarOn                  = True
        opts_rv@lbBoxEndCapStyle              = "TriangleBothEnds"
        opts_rv@lbOrientation                 = "Vertical"

        opts_rv@trYMaxF                       = lat_max + ar
        opts_rv@trYMinF                       = lat_max - ar
        opts_rv@trXMaxF                       = lon_max + ar
        opts_rv@trXMinF                       = lon_max - ar

      ; Do not draw the plot or advance the frame                                          
        opts_rv@gsnDraw                       = False
        opts_rv@gsnFrame                      = False
        opts_rv@gsnAddCyclic                  = False

      ; Background plot resources                                                          
        opts_xy                               = True

      ; No tickmarks on any axis                                                           
        opts_xy@tmYROn                        = False
        opts_xy@tmYLOn                        = True
        opts_xy@tmXTOn                        = False
        opts_xy@tmXBOn                        = True

        opts_xy@xyLineColor                   = -1
        opts_xy@gsnDraw                       = False
        opts_xy@gsnFrame                      = False

        opts_xy@trYMaxF                       = lat_max + ar
        opts_xy@trYMinF                       = lat_max - ar
        opts_xy@trXMaxF                       = lon_max + ar
        opts_xy@trXMinF                       = lon_max - ar
        opts_xy@gsnMaximize                   = False

      ; Resources for windspeed                                                            
        opts_wind                             = opts_rv
        delete([/opts_wind@cnLevels, opts_wind@cnFillColors/])
        opts_wind@cnFillOn                    = True
        opts_wind@cnLinesOn                   = False
        opts_wind@cnLineLabelsOn              = False
        opts_wind@cnInfoLabelOn               = False
        opts_wind@cnLevelSelectionMode        = "ExplicitLevels"
        opts_wind@cnFillPalette               = "radar_new"
        opts_wind@cnLevels                    = (/5.0, 10.0, 15.0, 20.0,\
                                                  25.0, 30.0, 35.0, 40.0, \
                                                  45.0, 50.0, 55.0, 60.0, 65.0/)
        opts_wind@cnFillColors                = (/-1,0,1,2,3,4,5,6,7,8,9,10,11,12/)
        opts_wind@gsnMaximize                 = False

      ; Resources for vertical velocity                                                    
        opts_vvel                             = opts_wind
        delete([/opts_vvel@cnLevels, opts_vvel@cnFillColors, opts_vvel@cnFillPalette/]) 
        opts_vvel@cnFillPalette               = "BlueDarkRed18"
        opts_vvel@cnLevels                    = (/-5.0, -3.0, -2.0, -1.5, -1.0, -0.5,\
                                                   0.5, 1.0, 1.5, 2.0, 3.0, 5.0/)
      ; Reverse order of colourmap (omega in pressure coords --> negative values = rising motion)
        opts_vvel@cnFillColors                = (/15,14,13,12,11,10,-1,7,6,5,4,3,2/)

      ; Resources for radial wind 
        opts_vrad                             = opts_wind
        delete([/opts_vrad@cnLevels, opts_vrad@cnFillColors, opts_vrad@cnFillPalette/])
        opts_vrad@cnFillPalette               = "BlueDarkRed18"
        opts_vrad@cnLevels                    = (/-20., -15., -10., -5., -3., -2., -1., -0.5, \
                                                   0, 0.5, 1., 2., 3., 5., 10., 15., 20./)
        opts_vrad@cnFillColors                = (/0,2,3,4,5,6,7,8,-1,-1,10,11,12,13,14,15,16,17/)

      ; Panel plot label resources                                                         
        txid                                   = new(4,graphic)
        amid                                   = new(4,graphic)

        panel_strings                          = (/ "a) Relative vorticity",\
                                                    "b) Vertical velocity",\
                                                    "c) Tangential wind",\
                                                    "d) Radial wind"/)

        txres                                  = True
        txres@txPerimOn                        = True
        txres@txFontHeightF                    = 0.0125
        txres@txBackgroundFillColor            = "White"

        amres                                  = True
        amres@amParallelPosF                   = -0.48
        amres@amOrthogonalPosF                 = -0.48
        amres@amJust                           = "TopLeft"

      ; Coordinates of panels 
        xf = (/0.05, 0.57, 0.05, 0.57/)
        yf = (/0.90, 0.90, 0.52, 0.52/)
        wf = 0.35
        hf = 0.32

      ; Panel 1 (relative vorticity)
        opts_rv@vpXF       = xf(0)
        opts_rv@vpYF       = yf(0)
        opts_rv@vpWidthF   = wf
        opts_rv@vpHeightF  = hf

        opts_xy0           = opts_xy
        opts_xy0@vpXF      = xf(0)
        opts_xy0@vpYF      = yf(0)
        opts_xy0@vpWidthF  = wf
        opts_xy0@vpHeightF = hf

        vort_int_plot      = vort_int * (10 ^ 4)

        plot_tl            = gsn_csm_xy(wks, xpos_all(:,:,ct), ypos_all(:,:,ct), opts_xy0)
        plot_tl0           = gsn_csm_contour(wks, vort_int_plot(:,:), opts_rv)
        overlay(plot_tl, plot_tl0)

      ; Overlay labels                                                                      
        txid(0)      = gsn_create_text(wks, panel_strings(0), txres)
        amid(0)      = gsn_add_annotation(plot_tl, txid(0), amres)

      ; Panel 2 (vertical velocity)                                                         
        opts_vvel@vpXF       = xf(1)
        opts_vvel@vpYF       = yf(1)
        opts_vvel@vpWidthF   = wf
        opts_vvel@vpHeightF  = hf

        opts_xy1             = opts_xy
        opts_xy1@vpXF        = xf(1)
        opts_xy1@vpYF        = yf(1)
        opts_xy1@vpWidthF    = wf
        opts_xy1@vpHeightF   = hf
 
        plot_tr              = gsn_csm_xy(wks, xpos_all(:,:,ct), \
                                          ypos_all(:,:,ct), opts_xy1)
        plot_tr0             = gsn_csm_contour(wks, w_int(:,:), opts_vvel)
        overlay(plot_tr, plot_tr0)

      ; Overlay labels                                                                      
        txid(1)      = gsn_create_text(wks, panel_strings(1), txres)
        amid(1)      = gsn_add_annotation(plot_tr, txid(1), amres)

      ; Panel 3 (tangential wind)                                                           
        opts_wind@vpXF       = xf(2)
        opts_wind@vpYF       = yf(2)
        opts_wind@vpWidthF   = wf
        opts_wind@vpHeightF  = hf

        opts_xy2             = opts_xy
        opts_xy2@vpXF        = xf(2)
        opts_xy2@vpYF        = yf(2)
        opts_xy2@vpWidthF    = wf
        opts_xy2@vpHeightF   = hf

        plot_bl              = gsn_csm_xy(wks, xpos_all(:,:,ct), \
                                          ypos_all(:,:,ct), opts_xy2)
        plot_bl0             = gsn_csm_contour(wks, vt_int(:,:), opts_wind)
        overlay(plot_bl, plot_bl0)

      ; Overlay labels                                                                      
        txid(2)      = gsn_create_text(wks, panel_strings(2), txres)
        amid(2)      = gsn_add_annotation(plot_bl, txid(2), amres)

      ; Panel 4 (radial wind)
        opts_strn0           = opts_vrad
        opts_strn0@vpXF      = xf(3)
        opts_strn0@vpYF      = yf(3)
        opts_strn0@vpWidthF  = wf
        opts_strn0@vpHeightF = hf
 
        opts_xy3             = opts_xy
        opts_xy3@vpXF        = xf(3)
        opts_xy3@vpYF        = yf(3)
        opts_xy3@vpWidthF    = wf
        opts_xy3@vpHeightF   = hf 

        plot_br              = gsn_csm_xy(wks, xpos_all(:,:,ct), \
                                          ypos_all(:,:,ct), opts_xy3)
        plot_br0             = gsn_csm_contour(wks, vr_int(:,:), opts_vrad)
        overlay(plot_br, plot_br0)

      ; Overlay labels                                                                      
        txid(3)      = gsn_create_text(wks, panel_strings(3), txres)
        amid(3)      = gsn_add_annotation(plot_br, txid(3), amres)

      ; Tidy up before final plots later in script                                          
        delete([/opts_vvel, opts_wind, opts_rv, opts_vrad, \
                 opts_xy, opts_xy0, opts_xy1, opts_xy2, opts_xy3, \
                 txid, amid, panel_strings/])

      ; Add radial circles to plot 
        dtr = 0.017453292519943

      ; Get customizations for circles                                                 
      ; Zoom into plot --> can't remember how to do this (29/10/2019)                  
        radius   = radii(rad_size)
        xcenter  = lon_max
        ycenter  = lat_max
        out_thck = 2
        in_thck  = 1
        spacing  = 0.05
        dpattrn  = 2
        out_clr  = 1
        in_clr   = 1

;=============================================                                       
; Construct and attach outer circle to plot                                          
;=============================================                                       

      ; Calculate arrays for outer circle                                              
        degrees  = ispan(0, 360, 5)
        xcos     = cos(dtr * degrees)
        xsin     = sin(dtr * degrees)

      ; Array of points (x,y) representing the outer circle                            
        xc       = xcenter + (radius * xcos)
        yc       = ycenter + (radius * xsin)

      ; Resources for outer circle                                                     
        lnres    = True
        lnres@gsLineThicknessF  = out_thck
        lnres@gsLineColor       = out_clr
        lnres@gsLineDashPattern = 0

      ; Attach outer circle to each panel plot                                         
        plot_outer0 = gsn_add_polyline(wks, plot_tl, xc, yc, lnres)
        plot_outer1 = gsn_add_polyline(wks, plot_tr, xc, yc, lnres)
        plot_outer2 = gsn_add_polyline(wks, plot_bl, xc, yc, lnres)
        plot_outer3 = gsn_add_polyline(wks, plot_br, xc, yc, lnres)

        if (spacing .gt. radius) then
         print("add_radial_circles: spacing is > radius, can't draw inner circles.")
         return
        end if

;==================================================                                    
; Now construct and attach inner circles to plot                                       
;==================================================                                    

    ; Draw inner circles if desired                                                    
      if (spacing .gt. 0) then

       count = 0
       size  = toint( (radius-spacing) / spacing )
       plot_inner = new( (/4,size/), graphic)

       do r1 = spacing, radius, spacing

       ; Calculate arrays for inner circles                                            
         xc = xcenter + (r1 * xcos)
         yc = ycenter + (r1 * xsin)

       ; Resources for inner circles                                                   
         delete([/lnres@gsLineThicknessF,lnres@gsLineColor,lnres@gsLineDashPattern/])
         lnres@gsLineColor       = in_clr
         lnres@gsLineThicknessF  = in_thck
         lnres@gsLineDashPattern = dpattrn

       ; Attach inner circles to plot                                                  
         plot_inner(0,count) = gsn_add_polyline(wks, plot_tl, xc, yc, lnres)
         plot_inner(1,count) = gsn_add_polyline(wks, plot_tr, xc, yc, lnres)
         plot_inner(2,count) = gsn_add_polyline(wks, plot_bl, xc, yc, lnres)
         plot_inner(3,count) = gsn_add_polyline(wks, plot_br, xc, yc, lnres)

         count = count + 1

       end do

      end if

    ; Tidy up                                                                          
      delete([/lnres, degrees, xcos, xsin, xc, yc/])

    ; Create a main title at the top of the page                                       
      title_res                   = True
      title_str                   = "Panel plot at "+p(vlev0)+" hPa: "+ens0+" ("+dat+")"
      title_res@txFontHeightF     = 0.02
      title_res@txFont            = "Helvetica"
      xpos_title                  = 0.50
      ypos_title                  = 0.94

      gsn_text_ndc(wks, title_str, xpos_title, ypos_title, title_res)

    ; Call 'draw' and 'frame'                                                          
      maximize_output(wks, True)

     end if 

;====================================
; End mid-script plotting routine 
;====================================

    ; Average variables azimuthally around the storm 
      do irad = 0, dimsizes(radii)-1
       r = radii(irad)
       vtan_azi(irad) = avg(vt_int(irad,:))
       vrad_azi(irad) = avg(vr_int(irad,:))
       vvel_azi(irad) = avg(w_int(irad,:))
       spd_azi(irad)  = avg(spd_int(irad,:))
       vort_azi(irad) = avg(vort_int(irad,:))
      end do

    ; Input azimuthal averages into larger arrays before plotting  
      vtan_plot0(:,it)  = vtan_azi(:)
      vrad_plot0(:,it)  = vrad_azi(:)
      vvel_plot0(:,it)  = vvel_azi(:)
      spd_plot0(:,it)   = spd_azi(:)
      vort_plot0(:,it)  = vort_azi(:)

    ; Tidy up
      delete([/u_int,v_int,w_int,vt_int,vr_int,spd_int,vort_int/])
      delete([/u_plane,v_plane,w_plane,vort_plane,spd_plane/])

  ; Advance counter variable (time)
    ct  = ct + 1

 ; Tidy up (avoid dimension size errors)
   delete([/lon,lat,u0,v0,w0,vort/])

  end do     ; End time loop (do it = ts, tf)

; Exit script if outputting 4-panel plots for each time (troubleshooting)
  if (plt .eq. 1) then
   exit()
  end if 

;========================
; Calculate tendencies
;========================

; Initialise new arrays
  vtan_tend0 	  = vtan_plot0
  vrad_tend0 	  = vrad_plot0

; Add metadata
  vtan_tend0!0    = "rad"
  vtan_tend0&rad  = radii
  vtan_tend0!1    = "time"
  vtan_tend0&time = time_int

  vrad_tend0!0    = "rad"
  vrad_tend0&rad  = radii
  vrad_tend0!1    = "time"
  vrad_tend0&time = time_int

  print("Calculating tendencies...")

; Loop over times to calculate tendencies 
  do it = ts, tf

   tp1  = min((/it+1,tf/))
   tm1  = max((/it-1,ts/))

   vtan_tend0(:,it) = (vtan_plot0(:,tp1) - vtan_plot0(:,tm1)) / 2 ; Units [m s-1 h-1]
   vrad_tend0(:,it) = (vrad_plot0(:,tp1) - vrad_plot0(:,tm1)) / 2 ; Units [m s-1 h-1]

  end do
; End time loop (do it = ts, tf)

;======================================
; Reorder dimensions before plotting 
;======================================

 ; Arrays now ordered [TIME * RADIUS] rather than [RADIUS * TIME]
   vtan_tend = vtan_tend0(time|:,rad|:)
   vrad_tend = vrad_tend0(time|:,rad|:)

   vtan_plot = vtan_plot0(time|:,rad|:)
   vrad_plot = vrad_plot0(time|:,rad|:)
   vvel_plot = vvel_plot0(time|:,rad|:)
   vort_plot = vort_plot0(time|:,rad|:)

;=====================================================
; Create a panel plot of Figs. 6a and 6b from N11 
;=====================================================

  ; Output file location and type
    output = "$nep/halong/images/hov_paper_"+dat+"_"+ens0+"_"+km
    wks    = gsn_open_wks(opt,output)

  ; Set up panel plot
    panel = new(3,graphic)

    gsn_define_colormap(wks,"n11")

;==========================
; Options for plotting 
;==========================

  ; Tangential wind
    opts_vtan                              = True
    opts_vtan@cnFillOn                     = True
    opts_vtan@cnLineLabelInterval          = 2.0
    opts_vtan@cnLineLabelFontHeightF       = 0.012
    opts_vtan@cnLineLabelBackgroundColor   = "transparent"
    opts_vtan@cnLineLabelPlacementMode     = "constant"
    opts_vtan@cnLinesOn                    = False
    opts_vtan@cnInfoLabelOn                = False
    opts_vtan@cnLevelSelectionMode         = "ExplicitLevels"

    opts_vtan@cnLevels                     = (/0.5, 1.0, 1.5, 2.0, 2.5, 3.0, \
                                               3.5, 4.0, 4.5, 5.0/)
    opts_vtan@cnFillColors                 = (/0,2,3,4,5,6,7,8,9,10,11/)

    opts_vtan@gsnPaperOrientation          = "landscape"
    opts_vtan@tiMainString                 = ""
    opts_vtan@tiMainFontHeightF            = 0.0125
    opts_vtan@gsnLeftString                = ""
    opts_vtan@gsnRightString               = ""
    opts_vtan@gsnMaximize                  = True
    opts_vtan@lbLabelBarOn                 = True                  ; Labelbar on/off
    opts_vtan@lbBoxEndCapStyle             = "TriangleBothEnds"    ; Labelbar end shape
    opts_vtan@pmLabelBarWidthF             = 0.32                  ; Labelbar width
    opts_vtan@pmLabelBarHeightF            = 0.08                  ; Labelbar height
    opts_vtan@lbLabelFontHeightF           = 0.0125                ; Labelbar font size
    opts_vtan@lbLabelFont                  = "Helvetica"           ; Labelbar font
    opts_vtan@lbPerimOn                    = False                 ; Perimeter on/off
    opts_vtan@gsnDraw                      = False                 ; Do not draw plot
    opts_vtan@gsnFrame                     = False                 ; Do not advance frame
    opts_vtan@vpWidthF                     = 0.30                  ; Match width in N11
    opts_vtan@vpHeightF                    = 0.50                  ; Match height in N11

  ; Additional plotting resources
    opts_vtan@tiYAxisString                = "Time (h)"
    opts_vtan@tiXAxisString                = "Radius (km)"
    opts_vtan@tiXAxisFontHeightF           = 0.018
    opts_vtan@trYLog                       = False

  ; Edit these values (10/02/2020)
    if (dat .eq. "02T12") then 
     opts_vtan@trYMinF                     = 20
     opts_vtan@trYMaxF                     = 32
    elseif (dat .eq. "03T00") then 
     opts_vtan@trYMinF                     = 22
     opts_vtan@trYMaxF                     = 34
    elseif (dat .eq. "04T00") then 
     opts_vtan@trYMinF                     = 14
     opts_vtan@trYMaxF                     = 30
    end if 

  ; Turn off top x-axis and right y-axis tickmarks
    opts_vtan@tmXTOn                       = "False"
    opts_vtan@tmYROn                       = "False"

  ; Maximise plot size
    opts_vtan@gsnMaximize                  = True

    opts_vtan@gsnAddCyclic                 = False
    opts_vtan@tmXBMode                     = "Explicit"                ; Choose tm mode
    opts_vtan@tmXBValues                   = (/0.0,0.5,1.0,1.5,2.0/)   ; Tm positions
    opts_vtan@tmXBLabels                   = (/"0","50","100",\
    					       "150","200"/)	       ; Tm values
    opts_vtan@tmXBMinorOn                  = False

    opts_vtan@pmLabelBarOrthogonalPosF     = 0.10    ; Move lb up/down (higher/lower)
    opts_vtan@tiXAxisOffsetYF              = 0.10    ; Move x-axis title up/down (h/l)

    opts_vtan@trXMinF                      = 0.0                       ; Set x-axis min
    opts_vtan@trXMaxF                      = ar0                       ; Set x-axis max
    opts_vtan@tmYLMode                     = "Manual"
    opts_vtan@tmYLTickStartF               = opts_vtan@trYMinF         ; See L730
    opts_vtan@tmYLTickEndF                 = opts_vtan@trYMaxF         ; See L731
    opts_vtan@tmYLTickSpacingF             = 6
    opts_vtan@tmYLMinorOn                  = False

  ; Resources to overlay contours of mean tangential wind
    opts_cont				   = True

    opts_cont@vpWidthF                     = 0.30                  ; Match width in N11
    opts_cont@vpHeightF                    = 0.50                  ; Match height in N11

    opts_cont@cnFillOn                     = False       ; Filled contours on/off
    opts_cont@cnLineColor                  = "blue3"     ; Contour line colour
    opts_cont@cnLevelSelectionMode         = "ExplicitLevels"
    opts_cont@cnLevels                     = (/10,20,30,40,50,60,70/)

    opts_cont@tiYAxisString                = "Time (h)"
    opts_cont@tiXAxisString                = "Radius (km)"
    opts_cont@tiXAxisFontHeightF           = 0.018
    opts_cont@trYLog                       = False

  ; Edit these values (10/02/2020) 
    opts_cont@trYMinF                      = ts
    opts_cont@trYMaxF                      = tf

  ; Turn off top x-axis and right y-axis tickmarks
    opts_cont@tmXTOn                       = "False"
    opts_cont@tmYROn                       = "False"

  ; Maximise plot size
    opts_cont@gsnMaximize                  = True

    opts_cont@gsnAddCyclic                 = False
    opts_cont@tmXBMode                     = "Explicit"                ; Choose tm mode
    opts_cont@tmXBValues                   = (/0.0,0.5,1.0,1.5/)       ; Tm positions
    opts_cont@tmXBLabels                   = (/"0","50","100","150"/)  ; Tm values
    opts_cont@tmXBMinorOn                  = False

    opts_cont@pmLabelBarOrthogonalPosF     = 0.10    ; Move lb up/down (higher/lower)
    opts_cont@tiXAxisOffsetYF              = 0.10    ; Move x-axis title up/down (h/l)

    opts_cont@trXMinF                      = 0.0                       ; Set x-axis min
    opts_cont@trXMaxF                      = 1.0                       ; Set x-axis max

    opts_cont@cnLineThicknessF             = 2.0         ; Contour line thickness
    opts_cont@cnInfoLabelOn                = False       ; Info label on/off
    opts_cont@gsnDraw                      = False       ; Do not draw the plot
    opts_cont@gsnFrame                     = False       ; Do not advance the frame
    opts_cont@cnLineLabelsOn               = True        ; Turn line labels on/off
    opts_cont@cnLineLabelFontColor         = "blue3"     ; Line label colour
    opts_cont@cnLineLabelFontHeightF       = 0.01        ; Default 0.006 for vpWidth=0.30
    opts_cont@cnLineLabelFontThicknessF    = 2.0         ; Line label thickness
    opts_cont@tiMainString                 = ""
    opts_cont@gsnLeftString                = ""
    opts_cont@gsnRightString               = ""

  ; Resources to overlay zero line
    opts_z                             = True
    opts_z@cnFillOn                    = False
    opts_z@cnLineColor                 = "gray82"
    opts_z@cnLevels                    = (/0.0/)
    opts_z@cnLineThicknessF            = 2.0
    opts_z@cnInfoLabelOn               = False
    opts_z@gsnDraw                     = False       ; Do not draw the plot
    opts_z@gsnFrame                    = False       ; Do no advance the frame
    opts_z@cnLineLabelsOn              = False       ; Turn off line labels
    opts_z@tiMainString                = ""
    opts_z@gsnLeftString               = ""
    opts_z@gsnRightString              = ""

  ; Vertical velocity resources
    opts_vvel                                = opts_vtan
    delete([/opts_vvel@cnLevels,opts_vvel@cnFillColors/])
    opts_vvel@cnFillPalette                  = "BlueDarkRed18"
    opts_vvel@cnLevels                       = (/-5., -3., -2., -1., \
                                                  0, 1., 2., 3., 5./)
    opts_vvel@cnFillColors                   = (/16,14,12,10,-1,-1,7,5,3,0/)

  ; Relative vorticity resources
    opts_vort                                = opts_vtan
    delete([/opts_vort@cnLevels, opts_vort@cnFillColors/])
    opts_vort@cnLevels                       = (/3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, \
                                                 15.0, 20.0, 25.0, 30.0, 35.0/)
    opts_vort@cnFillColors                   = (/0,2,3,4,5,6,7,8,9,10,11,12,13/)

  ; Radius of maximum wind (RMW) resources
    opts_rmw                               = True
    opts_rmw@xyLineColors                  = "black"
    opts_rmw@xyLineThicknesses             = 4.0
    opts_rmw@xyDashPatterns                = 0.0
    opts_rmw@gsnDraw                       = False       ; Do not draw the plot
    opts_rmw@gsnFrame                      = False       ; Do not advance the frame
    opts_rmw@tiMainString                  = ""
    opts_rmw@gsnLeftString                 = ""
    opts_rmw@gsnRightString                = ""

;========================
; Construct panel plot
;========================

  ; Read in RMW data from text file (first 3 values equal to zero)
;    if (rmw .eq. 1) then 
;     spd_plot = asciiread("$ar/text/mlev_"+dat+"_"+ens0+"_rmw_calc1.txt",(/numINPUT/),"float")
;    elseif (rmw .eq. 0) then 
;     spd_plot = asciiread("$ar/text/"+dat+"_"+ens0+"_rmw_calc1.txt",(/numINPUT/),"float")
;    end if 

    ; Tangential wind (LHS panel)
      vtan0 = gsn_csm_contour(wks,vtan_tend,opts_vtan)    ; Tangential wind tendency (filled)
      zero0 = gsn_csm_contour(wks,vtan_tend,opts_z)    	  ; Zero line (vtan tendency)
      zero  = ColorNegDashZeroPosContour(zero0,"transparent","gray69","transparent")
      cont0 = gsn_csm_contour(wks,vtan_plot,opts_cont)    ; Tangential wind
;      rmw0  = gsn_csm_xy(wks,spd_plot,time_int,opts_rmw)  ; RMW
      overlay(vtan0,zero)                                 ; Overlay vtan tendency zero line
;      overlay(vtan0,rmw0)				  ; Overlay RMW
      overlay(vtan0,cont0)                                ; Overlay tangential wind
      panel(0) = vtan0                                    ; Plot first panel
      delete([/zero0,zero/])                              ; Tidy up 

    ; Relative vorticity
      vort_plot = vort_plot * (10 ^ 4)
      vort_plot@units = "10~S~-4~N~ s~S~-1~N~"

      vort0 = gsn_csm_contour(wks,vort_plot,opts_vort)    ; Relative vorticity
      zero0 = gsn_csm_contour(wks,vtan_tend,opts_z)       ; Zero line (vtan tendency)          
      zero  = ColorNegDashZeroPosContour(zero0,"transparent","gray69","transparent")
      cont0 = gsn_csm_contour(wks,vtan_plot,opts_cont)    ; Tangential wind                    
;      rmw0  = gsn_csm_xy(wks,spd_plot,time_int,opts_rmw)  ; RMW                                
      overlay(vort0,zero)
      overlay(vort0,cont0)                                ; Overlay tangential wind            
;      overlay(vort0,rmw0)                                 ; Overlay RMW                        
      panel(1) = vort0                                    ; Plot second panel                  
      delete([/zero0,zero/])                              ; Tidy up 

    ; Vertical velocity (RHS panel)
      vvel0 = gsn_csm_contour(wks,vvel_plot,opts_vvel)    ; Vertical velocity
      zero0 = gsn_csm_contour(wks,vvel_plot,opts_z)       ; Zero line (vertical velocity)
      zero  = ColorNegDashZeroPosContour(zero0,"transparent","gray69","transparent")
;      rmw0  = gsn_csm_xy(wks,spd_plot,time_int,opts_rmw)  ; RMW
      overlay(vvel0,zero)                                 ; Overlay vvel zero line
;      overlay(vvel0,rmw0)                                 ; Overlay RMW
      panel(2) = vvel0                                    ; Plot third panel
      delete([/zero0,zero/])                              ; Tidy up

;=============================================
; Finally, put everything together and plot 
;=============================================

    optsP                       = True
    optsP@gsnFrame              = False    ; Do not advance the frame
    optsP@gsnPanelLabelBar      = False    ; Turn on panel labelbar
    optsP@txString              = "Simulation "+ens0+" ("+dat+")"

    optsP@gsnPanelFigureStrings = (/"a) Tangential wind","b) Relative vorticity",\
                                    "c) Vertical velocity"/)

    optsP@gsnMaximize           = True
    optsP@gsnPanelTop           = 0.90
    optsP@gsnPanelBottom        = 0.10
    optsP@amJust                = "TopLeft"
    optsP@gsnPanelFigureStringsFontHeightF = 0.0099 ; Reduce label size (default 0.01)

    gsn_panel(wks,panel,(/1,3/),optsP)             ; Draw as a single plot
    frame(wks)

  ; Print output file path to screen
    print(output)

end 