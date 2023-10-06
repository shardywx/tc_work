; Calculate composite plots of relative vorticity and strain in different VC phases

; Run using:

; ncl opt=\"pdf\" dist=3.0 nt=73 nr=41 ar=2.0 typ=\"sym\" 
; full=1 cn0=\"slp\" lay=0 mlev0=13 mlev1=25 tc_vstr_xyc_test.ncl

; 'opt'   = output file format ("pdf" or "x11")
; 'dist'  = size of box following storm (degrees)
; 'typ'	  = VC phase ("sym","asym","stoa","atos")
; 'full'  = top 18 (1); top 8 (-5); top 6 (-6)
; 'ar'    = distance in degrees of outer radius (1.0, 2.0, etc)
; 'nr'    = number of radial circles between r = 0 and r = 'ar' (21, 41, 61, etc)
; 'nt'    = number of azimuth angles in cylindrical grid
; 'cn0'   = offline storm centre: "slp", "vort", "geo", "geo_sm"

; 'pc' stream --> 'it'
; 'centre'    --> 'it-1'
; 'pd' stream --> 'it-2'

; Load main NCL functions and procedures
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"

; Also load user-defined functions and procedures
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"

begin

;==============================
; Prelash (reading in files)
;==============================

; For now, explicitly tell NCL which dates and simulations to analyse
  dat  = (/"02T12", "03T00", "03T12"/)

  ens0 = (/"em11", "em09", "em00"/)

  ts0  = (/50, 74, 33/)

  tf0  = (/53, 75, 36/)

;=====================================================
; Define cylindrical coordinate arrays for plotting 
;=====================================================

; Calculate total number of time intervals analysed
  num_int = sum( (tf0 - ts0) + 1)
  print("Total number of times in sample: "+num_int)

  thetas     = new(nt,float)
  pii        = 3.14159265
  radii      = fspan(0.0,ar,nr)      ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,nt)   ; Array of azimuth angles (0, 5, 10, etc)
  thetas_rad = (pii/180)*thetas      ; Degrees to radians
  rad_size   = dimsizes(radii)-1     ; Index for later calculation

; 3D arrays for all variables                                                                
  u_int0    = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)  ; Zonal wind         
  v_int0    = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)  ; Meridional wind    
  w_int0    = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)  ; Vertical velocity  
  str_int0  = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)  ; Strain rate        
  vort_int0 = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)  ; Relative vorticity

; Array to hold all 'xpos' and 'ypos' position arrays
  xpos_all  = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)
  ypos_all  = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)
  xcen_all  = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)
  ycen_all  = new((/dimsizes(radii),dimsizes(thetas),num_int/),float)

  u_int0@description   = "Zonal wind"
  u_int0@units         = "m s~S~-1~N~"
  u_int0!0             = "rad"
  u_int0!1             = "azi"
  u_int0!2             = "time"

  v_int0@description   = "Meridional wind"
  v_int0@units         = "m s~S~-1~N~"
  v_int0!0             = "rad"
  v_int0!1             = "azi"
  v_int0!2             = "time"

  w_int0@description   = "Vertical velocity"
  w_int0@units         = "m s~S~-1~N~"
  w_int0!0             = "rad"
  w_int0!1             = "azi"
  w_int0!2             = "time"

  str_int0@description = "Strain rate"
  str_int0@units       = "s~S~-2~N~"
  str_int0!0           = "rad"
  str_int0!1           = "azi"
  str_int0!2           = "time"

  vort_int0@description = "Relative vorticity"
  vort_int0@units       = "s~S~-1~N~"
  vort_int0!0           = "rad"
  vort_int0!1           = "azi"
  vort_int0!2           = "time"

; Initialise counter variable 
  ct                    = 0

;==========================================================
; Now start loop over these chosen input files and times
;==========================================================

 do st = 0, dimsizes(dat)-1 

  dt0 = tf0(st) - ts0(st)

; Find 'pc' and 'pd' files (model height level data)
  diri        = "/nfs/a319/earshar/"+dat(st)+"/"+ens0(st)
  fili_prefix = "201607"+dat(st)+"00Z_NPTK_4p4_L80_ra1t_"+ens0(st)
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
  time_d      = numINFO_d[1]		; Array of times from 'pd' stream file 

;============================================================
; Calculate storm motion using built-in function ('st_rm')
;============================================================

; 'dat'      = initialisation time (02T12, 03T00, ...)
; 'ens0'     = ensemble simulation (em00, em01, ...)
; 'diri'     = path to input files (see above)
; 'dist'     = size of box following storm (degrees)
; 'mins'     = analyse 1-h (0) or 5-min (1) data

  storm_rel  = st_rm(dat(st), ens0(st), "$ar/text/", dist, 0)
  u_cyc      = storm_rel[0]		   ; Zonal wind 
  v_cyc      = storm_rel[1]		   ; Meridional wind
  vel_cyc    = storm_rel[2]		   ; Vector wind 
  numTIMES   = storm_rel[3]                ; Number of times in file
  lat_arr    = storm_rel[4]                ; Storm-following domain subset (lat)
  lon_arr    = storm_rel[5]                ; Storm-following domain subset (lon)
  centre     = storm_rel[6]		   ; Storm centre at each hour (lat/lon)

;==========================================================================
; Create date/time string arrays for all times in file (use 'pc' stream)
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
  do ct0 = 0, numINPUT_c-1
   if (minute(ct0) .gt. 30) then
     hour(ct0) = hour(ct0)+1
   end if
  end do
  delete(ct0)

; Create arrays to hold strings for output plots
  date_str  = new(numINPUT_c,string)
  time_str  = new(numINPUT_c,string)
  time_arr  = new(numINPUT_c,string)
  title_arr = new(numINPUT_c,string)

  date_str  = sprinti("%0.2iUTC ", hour) + \
              sprinti("%0.2i ", day) \
            + month_abbr(month)
  time_str  = sprinti("%0.2i UTC ", hour) + \
              sprinti("%0.2i ", day) \
            + month_abbr(month)
  time_arr  = sprinti("%0.2i", day) + \
              month_abbr(month) + \
              "_" + sprinti("%0.2iZ", hour)

;==================================================================================
;;;;;;;;;;;;;;;;;;;;;;;;;;; Start of main part of script ;;;;;;;;;;;;;;;;;;;;;;;;;;
;==================================================================================

  print_clock("Starting main part of script...")

;=========================
; Start multiple loops
;=========================

   print_clock("Working on ensemble member: "+ens0(st))

   llbox       = toint(dist*50)	           ; Size of lat/lon arrays (see line below)
   dsize       = (/llbox,llbox/)           ; Size of lat/lon arrays (depends on 'dist')
   d0	       = dsize(0)
   d1	       = dsize(1)

   setvalues NhlGetWorkspaceObjectId
     "wsMaximumSize" : 1000000000
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

;===============================
; Start loop over input files 
;===============================

   c  = addfile(fili_c,"r")         	 ; Read in 'pc' stream [u,v,vort]
   d  = addfile(fili_d,"r")         	 ; Read in 'pd' stream [w,p,t,theta,kmh,tau]

;==================================
; Get the variables we will need
;==================================

   do it = ts0(st), tf0(st)		 ; Loop over times in file

    print("Working on time: "+time_str(it)+" (ct = "+ct+")" )
    itt = it
    title_arr(it) = "Valid at "+time_str(it)+" (T+"+itt+")"

;=====================================================
; Read in variables from 'pc' stream (model levels)
;=====================================================

    t0    = lat_arr(cn,it-1,0)
    t1    = lat_arr(cn,it-1,1)
    n0    = lon_arr(cn,it-1,0)
    n1    = lon_arr(cn,it-1,1)
    print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

  ; Read in basic variables from 'pc' stream
    lon1C  = c->longitude_1({n0:n1})              ; longitude ['d1' grid points]
    lat1C  = c->latitude_1({t0:t1})               ; latitude  ['d0' grid points]
    lonC   = c->longitude({n0:n1})                ; longitude ['d1' grid points]
    latC   = c->latitude({t0:t1})                 ; latitude  ['d0' grid points]
    hybC   = c->hybrid_ht(:) 		          ; 63 model ('theta') levels

  ; Calculate model level height information
    hy0    = sprintf("%0.0f",hybC(mlev0))
    hy1    = sprintf("%0.0f",hybC(mlev1))

    if (lay .eq. 1) then 
     print_clock("Working on layer between " +hy0+" and " +hy1+ "m AGL)" )
    elseif (lay .eq. 0) then 
     print_clock("Working on model level (" +hy0+ "m AGL)" )     
    end if 

  ; Read in horizontal wind components 
    u_varname = "x-wind"
    v_varname = "y-wind"

  ; Either read in multiple values over a layer 
    if (lay .eq. 1) then
     u0 = c->$u_varname$(it,mlev0:mlev1,{t0:t1},{n0:n1})
     v0 = c->$v_varname$(it,mlev0:mlev1,{t0:t1},{n0:n1})
  ; Or read in values on a single level
    else
     u0 = c->$u_varname$(it,mlev0,{t0:t1},{n0:n1})
     v0 = c->$v_varname$(it,mlev0,{t0:t1},{n0:n1})
    end if

;===================================
; Calculate storm-relative winds
;===================================

    u = u0 - u_cyc(cn,it-1)
    v = v0 - v_cyc(cn,it-1)

  ; Add metadata from original horizontal wind arrays
    copy_VarCoords(u0,u)
    copy_VarAtts(u0,u)

    copy_VarCoords(v0,v)
    copy_VarAtts(v0,v)

;==================================================================
; Calculate relative vorticity using centered finite differences
;==================================================================

  ; Option '2' --> boundary points estimated using one-sided difference scheme
    vort = uv2vr_cfd(u, v, latC, lonC, 2)
    copy_VarCoords(u, vort)
    vort@units = "s~S~-1~N~"
    vort@name = "Relative vorticity on model levels"

;===========================================================
; Calculate rate of strain (following Nguyen et al. 2011)
;===========================================================

    dim1   = "latitude"
    dim2   = "longitude"

  ; Calculate horizontal derivatives of 'u' (du_dx, du_dy)
    du     = grad_latlon_cfd(u,u&$dim1$,u&$dim2$,False,False)
    du_dy  = du[0]    ; Latitudinal derivative
    du_dx  = du[1]    ; Longitudinal derivative
    delete(du)

  ; Calculate horizontal derivatives of 'v' (dv_dx, dv_dy)
    dv     = grad_latlon_cfd(v,v&$dim1$,v&$dim2$,False,False)
    dv_dy  = dv[0]    ; Latitudinal derivative
    dv_dx  = dv[1]    ; Longitudinal derivative
    delete(dv)

  ; Calculate strain --> sqrt[ (du_dx - dv_dy)^2 + (dv_dx + du_dy)^2 ]
    str1a   = du_dx - dv_dy   ; Stretching deformation
    str2a   = dv_dx + du_dy   ; Shearing deformation
    str1    = (str1a)^2
    str2    = (str2a)^2

    if (lay .eq. 1) then 
     str0   = sqrt( str1(:,:,:) + str2(:,:,:) )
    else
     str0   = sqrt( str1(:,:) + str2(:,:) )
    end if

;======================================================
; Read in variables from 'pd' stream (model levels)
;======================================================

  ; Read in basic variables from 'pd' stream 
    lon1D  = d->longitude_1({n0:n1})                ; longitude [xxx grid points]
    lat1D  = d->latitude_1({t0:t1})                 ; latitude  [xxx grid points]
    lonD   = d->longitude({n0:n1})                  ; longitude [xxx grid points]
    latD   = d->latitude({t0:t1})                   ; latitude  [xxx grid points]
    hybD   = d->hybrid_ht(:)			    ; 63 model ('rho') levels

  ; Read in temperature, pressure, theta, vertical velocity and geopotential height
  ; All variables below -- [hybD | 63] * [latD | 150] * [lon1D | 151]
    if (lay .eq. 1) then 
     w     = d->dz_dt(it-2,mlev0:mlev1,{t0:t1},{n0:n1})
    elseif (lay .eq. 0) then 
     w     = d->dz_dt(it-2,mlev0,{t0:t1},{n0:n1})
    end if 

    dx0    = lonD(1) - lonD(0) ; Grid spacing (longitude)
    dy0    = latD(1) - latD(0) ; Grid spacing (latitude)

;==========================
; Loop over model levels 
;==========================

    if (lay .eq. 1) then
     vort_plane = dim_avg_n_Wrap(vort(:,:,:),0)
     w_plane    = dim_avg_n_Wrap(w(:,:,:),0)
     str_plane  = dim_avg_n_Wrap(str0(:,:,:),0)
     copy_VarCoords(u(0,:,:),str_plane)
    else
     vort_plane = vort
     w_plane    = w
     str_plane  = str0
     copy_VarCoords(u,str_plane)
    end if

  ; Read in precalculated storm centre position (change to 'it'?)
    lat_max    = centre(cn,it-1,0)
    lon_max    = centre(cn,it-1,1)

   ;======================================================================
   ; Create lat/lon arrays for later switch to cylindrical coordinates
   ;======================================================================

   ;===============
   ; 'pc' stream
   ;===============
  
   ; Indices of lat/lon grid for 'u' and 'v' in 'pc' stream
     latC_plane  = ind(t0 .le. latC .and. latC .le. t1)
     lonC_plane  = ind(n0 .le. lonC .and. lonC .le. n1)

   ; Array sizes
     latC_size   = dimsizes(latC_plane)
     lonC_size   = dimsizes(lonC_plane)
     
   ; Array of latitude/longitude points	corresponding to the indices above
     latC_sub    = latC(latC_plane(0):latC_plane(latC_size-1))
     lonC_sub    = lonC(lonC_plane(0):lonC_plane(lonC_size-1))

   ;===============
   ; 'pd' stream
   ;===============

   ; Indices of lat/lon grid for 'str_xz' in 'pd' stream
     latD_plane  = ind(t0 .le. latD .and. latD .le. t1)
     lonD_plane  = ind(n0 .le. lonD .and. lonD .le. n1)

   ; Array sizes
     latD_size   = dimsizes(latD_plane)
     lonD_size   = dimsizes(lonD_plane)

   ; Array of latitude/longitude points corresponding to the indices above
     latD_sub    = latD(latD_plane(0):latD_plane(latD_size-1))
     lonD_sub    = lonD(lonD_plane(0):lonD_plane(lonD_size-1))

   ; Indices of lat/lon grid for 'str_yz' in 'pd' stream
     lat1D_plane = ind(t0 .le. lat1D .and. lat1D .le. t1)
     lon1D_plane = ind(n0 .le. lon1D .and. lon1D .le. n1)

   ; Array sizes
     lat1D_size  = dimsizes(lat1D_plane)
     lon1D_size  = dimsizes(lon1D_plane)

   ; Array of latitude/longitude points corresponding to the indices above
     lat1D_sub   = lat1D(lat1D_plane(0):lat1D_plane(lat1D_size-1))
     lon1D_sub   = lon1D(lon1D_plane(0):lon1D_plane(lon1D_size-1))

;===========================================
; Switch to cylindrical coordinate system
;===========================================

   ; Create new arrays describing cylindrical coordinate system (2D for now)
     xpos     = new((/dimsizes(radii),dimsizes(thetas)/),float)
     ypos     = new((/dimsizes(radii),dimsizes(thetas)/),float)

     xcen     = new((/dimsizes(radii),dimsizes(thetas)/),float)
     ycen     = new((/dimsizes(radii),dimsizes(thetas)/),float)

     xpos!0   = "rad"
     xpos!1   = "azi"
     xpos&rad = radii
     xpos&azi = thetas

     ypos!0   = "rad"
     ypos!1   = "azi"
     ypos&rad = radii
     ypos&azi = thetas

     xcen!0   = "rad"
     xcen!1   = "azi"
     xcen&rad = radii
     xcen&azi = thetas

     ycen!0   = "rad"
     ycen!1   = "azi"
     ycen&rad = radii
     ycen&azi = thetas

   ; New arrays containing regular variables we want to plot
     w_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vertical velocity
     vort_int = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Relative vorticity
     str_int  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Strain rate

   ; Create cylindrincal coordinate array from original latitude/longitude grid
     do irad = 0, dimsizes(radii)-1
       r = radii(irad)
       do iang = 0, dimsizes(thetas)-1
          theta = thetas_rad(iang)   

           xpos(irad,iang) = lon_max + r*cos(theta)
 	   ypos(irad,iang) = lat_max + r*sin(theta)
	   xcen(irad,iang) = lon_max
	   ycen(irad,iang) = lat_max

       end do
     end do

   ; Fill 'xpos_all' and 'ypos_all' with latitude and longitude values
     xpos_all(:,:,ct) = xpos(:,:)
     ypos_all(:,:,ct) = ypos(:,:)
     xcen_all(:,:,ct) = xcen(:,:)
     ycen_all(:,:,ct) = ycen(:,:)     

   ; Fill new arrays with interpolated values on the cylindrical grid 
   ; Which latitude/longitude grids are each of these variables defined on? 
   ; vort_plane  -- latC,  lonC
   ; str_plane   -- latC,  lonC
   ; w_plane     -- latD,  lon1D

     do iang = 0, dimsizes(thetas_rad)-1

      vort_int(:,iang) = linint2_points(lonC_sub,latC_sub,vort_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      str_int(:,iang)  = linint2_points(lonC_sub,latC_sub,str_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      w_int(:,iang)    = linint2_points(lon1D_sub,latD_sub,w_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)

     end do

   ; Read interpolated grids into larger arrays 
     w_int0(:,:,ct)    = (/w_int(:,:)/)
     vort_int0(:,:,ct) = (/vort_int(:,:)/)
     str_int0(:,:,ct)  = (/str_int(:,:)/)

   ; Tidy up
     delete([/w_int,vort_int,str_int/])
     delete([/w_plane,vort_plane,str_plane/])
     delete([/xpos,ypos,xcen,ycen/])
     delete([/latC_plane,lonC_plane,latD_plane,lonD_plane,lon1D_plane,lat1D_plane/])
     delete([/latC_sub,lonC_sub,latD_sub,lonD_sub,lat1D_sub,lon1D_sub/])

     ct  = ct + 1		    ; Counter variable (time)

  ; Tidy up before next iteration (all variables from 'pc' and 'pd' streams)
    delete([/lonC,latC,lon1C,lat1C,u,v,vort/])
    delete([/lonD,latD,lon1D,lat1D,w/])

   end do     ; End time loop (do it = 0, times-1)

  end do      ; End loop over input file times (do st = 0, dimsizes(dat)-1)

;=================================================
; Produce time-averaged arrays of key variables
;=================================================

 ; Average along time dimension 
   w_arr    = dim_avg_n_Wrap(w_int0(:,:,:), 2)
   str_arr  = dim_avg_n_Wrap(str_int0(:,:,:), 2)
   vort_arr = dim_avg_n_Wrap(vort_int0(:,:,:), 2)

 ; EDIT - how do we give these final arrays fixed lat/lon grids that we can use to produce plots?
 ; At the moment there's no indication of radius or azimuth angle  

;==========================
; Options for plotting 
;==========================

 ; Relative vorticity
   opts_rv                               = True
;   opts_rv@sfXArray                      = xpos_all(:,:,0)
;   opts_rv@sfYArray                      = ypos_all(:,:,0)
   opts_rv@cnFillOn                      = True
   opts_rv@cnLineLabelInterval           = 2.0
   opts_rv@cnLineLabelFontHeightF        = 0.012
   opts_rv@cnLineLabelBackgroundColor    = "transparent"
   opts_rv@cnLineLabelPlacementMode      = "constant"
   opts_rv@cnLinesOn                     = False ; Contour lines off              
   opts_rv@cnInfoLabelOn                 = False ; Contour labels off             
   opts_rv@cnLevelSelectionMode          = "ExplicitLevels"
   opts_rv@cnFillPalette                 = "prcp_new"
   opts_rv@cnLevels                      = (/5., 10., 15., 20., 25., 30., \
                                             40., 50., 60., 70., 80./)
   opts_rv@cnFillColors                  = (/-1,2,3,4,5,6,7,8,9,10,11,12/)

   opts_rv@gsnPaperOrientation           = "landscape"
   opts_rv@tiMainString                  = ""
   opts_rv@tiMainFontHeightF             = 0.0125
   opts_rv@gsnLeftString                 = ""
   opts_rv@gsnRightString                = ""

 ; Labelbar resources
   opts_rv@lbLabelBarOn                  = True
;   opts_rv@pmLabelBarWidthF              = 0.52
;   opts_rv@pmLabelBarHeightF             = 0.13
;   opts_rv@lbLabelFontHeightF            = 0.0125
   opts_rv@lbPerimOn                     = False
   opts_rv@lbLabelFont                   = "Helvetica"
   opts_rv@lbBoxEndCapStyle              = "TriangleBothEnds"    ; Labelbar end shape
   opts_rv@lbTitleFontHeightF            = 0.02
   opts_rv@lbTitleDirection              = "Across"
   opts_rv@lbTitlePosition               = "Bottom"

 ; Plot size and position resources 
   opts_rv@vpXF                          = 0.07
   opts_rv@vpYF                          = 0.80
   opts_rv@vpWidthF                      = 0.35
   opts_rv@vpHeightF                     = 0.35

   if (lay .eq. 0) then

    opts_rv@lbTitleString                = hy0+\
                                              " m relative vorticity (10~S~-4~N~ s~S~-1~N~)"

   elseif (lay .eq. 1) then

    opts_rv@lbTitleString             = "Layer-averaged ("+hy0+"-"+hy1+\
                                        " hPa) ~C~     relative vorticity (10~S~-4~N~ s~S~-1~N~)"

   end if

 ; Additional plotting resources                                                        
   opts_rv@pmTickMarkDisplayMode         = "Always"   ; Turn on map tickmarks           
   opts_rv@tmXBMajorLengthF              = 0.005      ; Change tickmark length          
   opts_rv@tmXTOn                        = "False"    ; No tickmarks on top x-axis      
   opts_rv@tmYROn                        = "False"    ; No tickmarks on right y-axis    
   opts_rv@gsnAddCyclic                  = False
   opts_rv@gsnDraw                       = False      ; Do not draw the plot            
   opts_rv@gsnFrame                      = False      ; Do not advance the frame

;=========================================================                                  
; Now plot rate of strain using different colour table                                      
;=========================================================                                  

 ; Rate of strain                                                                         
   opts_str                              = opts_rv

 ; Contour line resources
   opts_str@cnFillOn                     = True
   opts_str@cnLineLabelInterval          = 2.0
   opts_str@cnLineLabelFontHeightF       = 0.012
   opts_str@cnLineLabelBackgroundColor   = "transparent"
   opts_str@cnLineLabelPlacementMode     = "constant"
   opts_str@cnLinesOn                    = False ; Contour lines off                      
   opts_str@cnInfoLabelOn                = False ; Contour labels off                     
   opts_str@cnFillPalette                = "circ_budget"
   opts_str@cnLevelSelectionMode         = "ExplicitLevels"
   delete([/opts_str@cnLevels, opts_str@cnFillColors/])
   opts_str@cnLevels                     = (/10,12,14,16,18,20,22,25,30/)
   opts_str@cnFillColors                 = (/-1,11,13,15,16,17,18,19,20,21/)

 ; Main plot resources
   opts_str@gsnPaperOrientation          = "landscape"
   opts_str@tiMainString                 = ""
   opts_str@tiMainFontHeightF            = 0.0125
   opts_str@gsnLeftString                = ""
   opts_str@gsnRightString               = ""

 ; Labelbar resources
   opts_str@lbLabelBarOn                 = True             ; Individual labelbar 
;   opts_str@pmLabelBarWidthF             = 0.52
;   opts_str@pmLabelBarHeightF            = 0.13
;   opts_str@lbLabelFontHeightF           = 0.0125
   opts_str@lbPerimOn                    = False
   opts_str@lbLabelFont                  = "Helvetica"

 ; Plot size and position resources
   opts_str@vpXF                         = 0.55
   opts_str@vpYF                         = 0.80
   opts_str@vpWidthF                     = 0.35
   opts_str@vpHeightF                    = 0.35

   if (lay .eq. 0) then
    opts_str@lbTitleString = hy0+" m strain rate (10~S~-4~N~ s~S~-2~N~)"
   elseif (lay .eq. 1) then
    opts_str@lbTitleString = "Layer-averaged ("+hy0+"-"+hy1+\
                             " m) ~C~         strain rate (10~S~-4~N~ s~S~-2~N~)"
   end if

   opts_str@lbBoxEndCapStyle             = "TriangleBothEnds"    ; Labelbar end shape     
   opts_str@lbTitleFontHeightF           = 0.02
   opts_str@lbTitleDirection             = "Across"
   opts_str@lbTitlePosition              = "Bottom"

;=============                                                                              
; Plot data                                                                                 
;=============                                                                              

 ; Output file location and type
   output = "$sam/nepartak/images/cyl_coords/rv_str_comp_"+typ+"_full"+full+"_"+cn0
   wks = gsn_open_wks(opt,output)

 ; Relative vorticity
   vort_arr    = vort_arr * (10 ^ 4)
   vort_arr@units = "10~S~-4~N~ s~S~-1~N~"
   plot_rv     = gsn_csm_contour(wks,vort_arr,opts_rv) ; Relative vorticity         

 ; Strain rate
   str_arr     = str_arr * (10 ^ 4)
   str_arr@units = "10~S~-4~N~ s~S~-2~N~"
   plot_str    = gsn_csm_contour(wks,str_arr,opts_str) ; Strain rate

 ; Troubleshooting 
   printMinMax(vort_arr, False)
   printMinMax(str_arr, False)

 ; Plot multiple panels without using 'gsn_panel' function
   pres = True
   maximize_output(wks,pres)

end 