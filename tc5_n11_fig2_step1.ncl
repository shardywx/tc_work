; Script to calculate the RMW and max 'v' at each time interval, and write values to text file 

; Run using:

; ncl dat=\"02T12\" opt=\"x11\" ens0=\"em11\" dist=3.0 ts0=947 tf0=995
; ar=1.5 nr=31 calc=1 cn0=\"slp\" tc5_n11_fig2_step1.ncl

; T+12    = 131 (144)
; T+24    = 275 (288)
; T+36    = 419 (432)
; T+48    = 563 (576)
; T+60    = 707 (720)
; T+72    = 851 (864)
; T+84    = 995 (1008)

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'opt'   = output file format ("pdf" or "x11")
; 'ens0'  = ensemble member (em00 -- em11)
; 'dist'  = size of box following storm (degrees)
; 'ts0'	  = starting time (0–120)
; 'tf0'   = ending time (0–120)
; 'ar'    = distance in degrees of outer radius (1.0, 2.0, etc)
; 'nr'    = number of radial circles between r = 0 and r = 'ar' (21, 41, 61, etc)
; 'calc'  = centre calculated offline (1) or using 980 m vort (2)

; 'pc' stream --> 'it'
; 'centre'    --> 'it-1'
; 'pd' stream --> 'it-2'

load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_rm.ncl"

begin

;==============================
; Prelash (reading in files)
;==============================

; Find 'pb' stream data to read in and plot
  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0+"/vc"
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0
  fili_c1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pc.nc")
  fili_d1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pd.nc")

  fili_c      = diri+"/"+fili_c1
  fili_d      = diri+"/"+fili_d1

  if ( .not. fileexists (fili_d) ) then
   print("Input file (pd stream) missing! Exiting script...")
   exit()
  elseif ( .not. fileexists (fili_c) ) then
   print("Input file (pc stream) missing! Exiting script...")
   exit()
  end if

; Calculate number of times in file using built-in function ('nc_times')
  numINFO_c   = nc_times(fili_c)
  numINFO_d   = nc_times(fili_d)
  numINPUT_c  = numINFO_c[0]		; Number of times in 'pc' stream
  time_c      = numINFO_c[1]		; Array of times from 'pc' stream file
  numINPUT_d  = numINFO_d[0]            ; Number of times in 'pd' stream
  time_d      = numINFO_d[1]            ; Array of times from 'pd' stream file

; Define number of vertical levels in input files 
  num_l	      = 63

; Calculate storm motion using built-in function ('st_rm')
; 'dat'     = initialisation time (02T12, 03T00, etc)
; 'ens0'    = ensemble simulation (em00-em11)
; 'diri'    = path to input files (see above)
; 'dist'    = size of box following storm (degrees)
; 'mins'    = analyse 1-h (0) or 5-min (1) data

  storm_rel = st_rm(dat, ens0, "$ar/text/", dist, 1)
  u_cyc     = storm_rel[0]             ; Zonal wind
  v_cyc     = storm_rel[1]             ; Meridional wind
  vel_cyc   = storm_rel[2]             ; Vector wind
  lat_arr   = storm_rel[4]             ; Latitude array subset (following storm)
  lon_arr   = storm_rel[5]             ; Longitude array subset (following storm)
  centre    = storm_rel[6]             ; Storm track (position) array

; Create arrays for data and output information
  llbox     = toint(dist*50)
  dsize     = (/llbox,llbox/)
  d0        = dsize(0)
  d1        = dsize(1)

  print("Working on ensemble member: "+ens0)

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 10000000000
  end setvalues

;==========================================================================
; Create date/time string arrays for all times in file (use 'pc' stream)
;==========================================================================

  month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
                    "Jul","Aug","Sep","Oct","Nov","Dec"/)

  utc_date   = cd_calendar(time_c, 0)

  year       = tointeger(utc_date(:,0))
  month      = tointeger(utc_date(:,1))
  day        = tointeger(utc_date(:,2))
  hour       = tointeger(utc_date(:,3))
  minute     = tointeger(utc_date(:,4))
  second     = utc_date(:,5)

; Create array to hold finished date strings
  hr_min = new(numINPUT_c,string)
  min0   = (/4,9,14,19,24,29,34,39,44,49,54,59/)      ; Minute array (before)
  min1   = (/5,10,15,20,25,30,35,40,45,50,55,0/)      ; Minute array (after)

; Before correcting, find the indices of the elements where 'minute' = 59
  hr_ind = ind(minute .eq. 59)

; Correct for errors in the code (round up values of 'minute' to multiples of 5)
  do i = 0, numINPUT_c-1
   do m = 0, dimsizes(min0)-1

    if (minute(i) .eq. min0(m) ) then
     minute(i) = min1(m)
    end if

   end do
  end do

; Edit all incorrect values of 'hour' (where 'minute' previously = 59)
  do j = 0, dimsizes(hr_ind)-1
   hour(hr_ind(j)) = hour(hr_ind(j)) + 1
  end do
  delete(j)

; Create string array of all values in 'minute'
  mins = tostring(minute)

; Create 'hhmm' strings for each time (for output and title strings)
  do i = 0, numINPUT_c-1

 ; Add zeros in front of single-digit hours
   if (hour(i) .lt. 10) then

    if (minute(i) .lt. 10) then
     hr_min(i) = "0"+hour(i)+"0"+mins(i)
    elseif (minute(i) .ge. 10) then
     hr_min(i) = "0"+hour(i)+mins(i)
    end if

   else

  ; Also add zeros in front of single-digit minutes
    if (minute(i) .lt. 10) then
     hr_min(i) = hour(i)+"0"+mins(i)
    elseif (minute(i) .ge. 10) then
     hr_min(i) = hour(i)+mins(i)
    end if

   end if

  end do
  delete(i)

  time_str  = new(numINPUT_c,string)
  time_arr  = new(numINPUT_c,string)
  title_arr = new(numINPUT_c,string)

  time_str  = hr_min + " UTC " + sprinti("%0.2i ", day) + month_abbr(month)
  time_arr  = sprinti("%0.2i", day) + month_abbr(month) + "_" + hr_min + "Z"
  title_arr = "Valid at "+time_str

;=======================================
; Choose which storm track to read in
;=======================================

  if (cn0 .eq. "slp") then
   cn = 0
  elseif (cn0 .eq. "slpf") then
   cn = 1
  end if

;==================================================================================
;;;;;;;;;;;;;;;;;;;;;;;;;;; Start of main part of script ;;;;;;;;;;;;;;;;;;;;;;;;;;
;==================================================================================

; Define cylindrical coordinate arrays
  pii        = 3.14159265
  radii      = fspan(0.0,ar,nr)      ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,73)   ; Array of azimuth angles (0, 5, 10, etc)
  thetas_rad = (pii/180)*thetas      ; Degrees to radians
  rad_size   = dimsizes(radii)-1     ; Index for later calculation

; Define arrays for later use in radius-time plots
  spd_azi               = new((/dimsizes(radii)/),float)        ; Vector wind [ {V} ]
  spd_azi!0             = "rad"
  spd_azi&rad           = radii
  spd_azi@description   = "Azimuthally averaged vector wind"
  spd_azi@units         = "m s~S~-1~N~"

  vtan_azi              = new((/dimsizes(radii)/),float)        ; Vector wind [ {V} ]
  vtan_azi!0            = "rad"
  vtan_azi&rad          = radii
  vtan_azi@description  = "Azimuthally averaged tangential wind"
  vtan_azi@units        = "m s~S~-1~N~"

; Define 3D arrays (hold values from above arrays at ALL times)

; Max. vector wind (2D)
  spd_max          = new((/num_l,numINPUT_c/),float)
  spd_max!0        = "lev"
  spd_max!1        = "time"

; Max. tangential wind (2D)
  vtan_max         = new((/num_l,numINPUT_c/),float)
  vtan_max!0       = "lev"
  vtan_max!1       = "time"

; Max. tangential wind (1D)
  vtan_max0        = new((/numINPUT_c/),float)
  vtan_max0!0      = "time"

; RMW array (3D)
  rmw_arr          = new((/2,num_l,numINPUT_c/),float)
  rmw_arr!0	   = "type"
  rmw_arr!1	   = "lev"
  rmw_arr!2	   = "time"

  rmw_ind	   = new((/2,numINPUT_c/),float)
  rmw_ind!1	   = "time"

;=========================
; Start multiple loops
;=========================

; Start loop over input files
  c  = addfile(fili_c,"r")         ; Read in 'pc' stream [u,v,vort]
  d  = addfile(fili_d,"r")	   ; Read in 'pd' stream [p]

  do it = ts0, tf0                 ; Loop over times in file

   print("Working on time: "+time_str(it))
   itt = it + 1
   title_arr(it) = "Valid at "+time_str(it)+" (T+"+itt+")"

 ; Get domain boundaries at each time from storm track (calculated off-line)
   t0     = lat_arr(cn,it,0)
   t1     = lat_arr(cn,it,1)
   n0     = lon_arr(cn,it,0)
   n1     = lon_arr(cn,it,1)
   print_clock("Grid (pc stream) extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

  ; Standard latitude/longitude grid
    lon    = c->longitude({n0:n1})                ; longitude ['d1' grid points]
    lat    = c->latitude({t0:t1})                 ; latitude  ['d0' grid points]

    dx0    = lon(1) - lon(0) ; Grid spacing (longitude)
    dy0    = lat(1) - lat(0) ; Grid spacing (latitude)

  ; Model levels (vertical)
    hybC   = c->hybrid_ht(:)                      ; 63 model ('theta') levels

    u_varname = "x-wind"
    v_varname = "y-wind"

  ; Read in horizontal wind components
  ; Both variables below -- [hybC | 63] * [latC | 150] * [lonC | 150]
    u0  = c->$u_varname$(it,:,{t0:t1},{n0:n1}) ; Zonal wind
    v0  = c->$v_varname$(it,:,{t0:t1},{n0:n1}) ; Meridional wind

;===================================                                                 
; Calculate storm-relative winds                                                     
;===================================                                                 

    u   = u0 - u_cyc(cn,it)
    v   = v0 - v_cyc(cn,it)

  ; Add metadata from original horizontal wind arrays                                
    copy_VarCoords(u0,u)
    copy_VarAtts(u0,u)

    copy_VarCoords(v0,v)
    copy_VarAtts(v0,v)

;==================================================================
; Calculate relative vorticity using centered finite differences
;==================================================================

  ; Option '2' --> boundary points estimated using one-sided difference scheme       
    vort = uv2vr_cfd(u, v, lat, lon, 2)
    copy_VarCoords(u, vort)
    vort@units = "s~S~-1~N~"
    vort@name = "Relative vorticity on model levels"

;======================================================                              
; Read in variables from 'pd' stream (model levels)                                  
;======================================================                              

  ; Read in basic variables from 'pd' stream                                         
    hybD   = d->hybrid_ht(:)                        ; 63 model ('rho') levels        

  ; Read in temperature, pressure, vertical velocity and geopotential height         
  ; All variables below -- [hybD | 63] * [latD | 150] * [lon1D | 151]                
    p      = d->p(it,:,{t0:t1},{n0:n1})             ; Pressure (Pa)                  
    z      = d->ht(it,:,{t0:t1},{n0:n1})            ; Geopotential height (m)

;==============================================================                      
; Make sure all arrays are same size (for later calculations)                        
;==============================================================                      

  ; 2-D array (4 * 3)                                                                 
    size  = (/dimsizes(u),dimsizes(v),dimsizes(p),dimsizes(vort)/)
    dsize = dimsizes(size)
    
  ; Reduce size of any dimension larger than the rest (both u and v)                  
    do sz  = 0, 2
     nind = min( size(:,sz) )
     if ( any (size(:,sz) .gt. nind) ) then
      size(:,sz) = nind
     end if
    end do

  ; Resize arrays based upon results of calculations above                            
    u2    = u(0:size(0,0)-1,0:size(0,1)-1,0:size(0,2)-1)
    v2    = v(0:size(1,0)-1,0:size(1,1)-1,0:size(1,2)-1)
    p2    = p(0:size(2,0)-1,0:size(2,1)-1,0:size(2,2)-1)
    vort2 = vort(0:size(3,0)-1,0:size(3,1)-1,0:size(3,2)-1)
    lat2  = lat(0:nind-1)
    lon2  = lon(0:nind-1)
    delete([/u,v,p,vort,lat,lon/])

  ; Reset names to avoid further code changes                                         
    u     = u2
    v     = v2
    vort  = vort2
    p     = p2
    lat   = lat2
    lon   = lon2
    delete([/u2,v2,p2,vort2,lat2,lon2/])

;=====================================
; Calculate storm centre at ~ 1 km
;=====================================

  ; Model level closest to 1 km AGL (980 m)
    ilev = 13

  ; Read in relative vorticity on chosen model level
    vort_plane0 = vort(ilev,:,:)

  ; Calculate location of vorticity maximum (storm centre) on chosen model level
    vort_max0               = max(vort_plane0)
    dims0                   = dimsizes(vort_plane0)
    vort1d0                 = ndtooned(vort_plane0)
    inds0                   = ind_resolve(maxind(vort1d0),dims0)
    delete([/vort_plane0, vort_max0, dims0, vort1d0/])

;==========================
; Loop over model levels 
;==========================

    do ilev = 0, num_l-1 ; Loop over all model levels

     hy = sprintf("%0.0f",hybD(ilev))
     print("Working on model level " + ilev+ " ("+hy+" m AGL)")

   ; Read in variables on chosen model level
     u_plane       = u(ilev,:,:)
     v_plane       = v(ilev,:,:)
     vort_plane	   = vort(ilev,:,:)
     prs_plane     = p(ilev,:,:)

   ; Calculate maximum relative vorticity on model level
     vort_max      = max(vort_plane)
     vort_smth     = smth9_Wrap(vort_plane, 0.5, 0.5, True)
     vort_max_smth = max(vort_smth)

   ; Find index of relative vorticity maximum 
     dims          = dimsizes(vort_plane)
     vort1d        = ndtooned(vort_plane)
     inds          = ind_resolve(maxind(vort1d),dims)
     vort1d_sm     = ndtooned(vort_smth)
     inds_sm       = ind_resolve(maxind(vort1d_sm),dims)

   ; Find storm centre
     if (calc .eq. 1) then			 ; Use offline values ('centre')

      lat_max       = centre(cn,it,0)
      lon_max       = centre(cn,it,1)

     elseif (calc .eq. 2) then			 ; Use vorticity maximum at 980 m (mlev=13)

      lat_max       = t0 + (dy0 * inds0(0,0))
      lon_max       = n0 + (dx0 * inds0(0,1))

     end if 

   ;==============================================================
   ; Create arrays for later switch to cylindrical coordinates
   ;==============================================================

   ;===============
   ; 'pc' stream
   ;===============

   ; Indices of lat/lon grid for 'u' and 'v' in 'pc' stream
     lat_plane  = ind(t0 .le. lat .and. lat .le. t1)
     lon_plane  = ind(n0 .le. lon .and. lon .le. n1)

   ; Array sizes
     lat_size   = dimsizes(lat_plane)
     lon_size   = dimsizes(lon_plane)

   ; Array of latitude/longitude points corresponding to the indices above
     lat_sub    = lat(lat_plane(0):d0-1)
     lon_sub    = lon(lon_plane(0):d1-1)

     delete([/vort1d,inds,vort1d_sm,vort_smth,inds_sm/])

;===========================================
; Switch to cylindrical coordinate system
;===========================================

   ; Create new arrays describing cylindrical coordinate system (2D for now)
     xpos     = new((/dimsizes(radii),dimsizes(thetas)/),float)
     ypos     = new((/dimsizes(radii),dimsizes(thetas)/),float)

   ; New arrays containing variables we want to plot
     u_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Zonal wind
     v_int    = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Meridional wind
     spd_int  = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Vector wind (RMW)
     vtan_int = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Tangential wind
     vrad_int = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Radial wind 

;=========================
; Add variable metadata
;=========================

     copy_VarAtts(u_plane,u_int)
     u_int@description = "Zonal wind"
     u_int@units       = "m s~S~-1~N~"

     copy_VarAtts(v_plane,v_int)
     v_int@description = "Meridional wind"
     v_int@units       = "m s~S~-1~N~"

     xpos!0    = "rad"
     xpos!1    = "azi"
     xpos&rad  = radii
     xpos&azi  = thetas

     ypos!0    = "rad"
     ypos!1    = "azi"
     ypos&rad  = radii
     ypos&azi  = thetas

     u_int!0   = "rad"
     u_int!1   = "azi"
     u_int&rad = radii
     u_int&azi = thetas

     v_int!0   = "rad"
     v_int!1   = "azi"
     v_int&rad = radii
     v_int&azi = thetas

     spd_int!0   = "rad"
     spd_int!1   = "azi"
     spd_int&rad = radii
     spd_int&azi = thetas

     vtan_int!0   = "rad"
     vtan_int!1   = "azi"
     vtan_int&rad = radii
     vtan_int&azi = thetas

   ; Latitude/longitude (small grid)
     do irad = 0, dimsizes(radii)-1
       r = radii(irad)
       do iang = 0, dimsizes(thetas)-1
         theta = thetas_rad(iang)
           xpos(irad,iang) = lon_max + r*cos(theta)
           ypos(irad,iang) = lat_max + r*sin(theta)
       end do
     end do

  ; Loop over azimuth angles to fill entire circular array of data points
    do iang = 0, dimsizes(thetas_rad)-1
      u_int(:,iang)    = linint2_points(lon_sub,lat_sub,u_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
      v_int(:,iang)    = linint2_points(lon_sub,lat_sub,v_plane,False,\
                                        xpos(:,iang),ypos(:,iang),0)
    end do

    spd_int = sqrt((u_int ^ 2) + (v_int ^ 2))
    copy_VarCoords(u_int,spd_int)
    copy_VarMeta(u_int,spd_int)
    spd_int@units = "m s~S~-1~N~"

  ; Calculate tangential wind at each point in array
    do iang = 0, dimsizes(thetas)-1
      theta = thetas_rad(iang)
      vtan_int(:,iang) = -u_int(:,iang)*sin(theta) + v_int(:,iang)*cos(theta) ; Tangential wind
    end do  

  ;===================================================
  ; Average variables azimuthally around the storm 
  ;===================================================

    do irad = 0, dimsizes(radii)-1
      spd_azi(irad)     = avg(spd_int(irad,:))	 ; Vector wind (for RMW)
      vtan_azi(irad)    = avg(vtan_int(irad,:))	 ; Tangential wind (for mass flux plots)
    end do

  ; Calculate maximum tangential and total windspeed on each model level
    vtan_max(ilev,it)	= max(vtan_azi)
    spd_max(ilev,it)	= max(spd_azi)

  ; Calculate RMW on each model level
    dims_s	        = dimsizes(spd_azi)
    spd_1d		= ndtooned(spd_azi)
    inds_s		= ind_resolve(maxind(spd_1d),dims_s)
    rmw		        = radii(inds_s(0,0))

  ; Input values of RMW on each level into larger array
    rmw_arr(0,ilev,it)  = rmw		    ; Radius (degrees)
    rmw_arr(1,ilev,it)  = (/inds_s(0,0)/) ; Index corresponding to radius

  ; Tidy up
    delete([/u_int,v_int,spd_int,vtan_int,xpos,ypos/])
    delete([/u_plane,v_plane,vort_plane,prs_plane/])
    delete([/lat_plane,lon_plane,lat_sub,lon_sub/])
    delete([/dims_s,spd_1d,inds_s,rmw/])

    end do		    ; End model level loop (do ilev = 0, num_l-1)

 ; Tidy up (avoid dimension size errors)
   delete([/lon,lat,u,v,u0,v0,vort,p/])

 end do     ; End time loop (do it = ts0, tf0)

;==========================================================
; At each time, calculate the maximum averaged windspeed
;==========================================================

 ; 'rmw_arr' contains the RMW on each pressure level, at each time
 ; Now calculate which level contains the largest maximum 

 ; rmw_arr --> [TYPE * LEV * TIME]

 ; (1) Calculate the largest maximum wind speed at each time (of all model levels)
 ; (2) Calculate the index of this maximum 
 ; (3) Calculate the corresponding radius of maximum wind 
	
   do it = ts0, tf0 

     print("Working on time: "+time_str(it)+" (it = "+it+")")

   ; Calculate maximum total and tangential wind (ignore upper troposphere)
     spd_max0        = max(spd_max(0:40,it))              ; Total wind
     vtan_max0(it)   = max(vtan_max(0:40,it))		  ; Tangential wind

   ; Find index of model level on which maximum tangential wind occurs
     dims_s          = dimsizes(spd_max(0:40,it))	  ; Size of model level array
     spd_1d          = ndtooned(spd_max(0:40,it))  	  ; Make 1D array
     inds_s  	     = ind_resolve(maxind(spd_1d),dims_s) ; Index of maximum total wind (mlev)
     h0		     = hybD(inds_s(0,0)) 		  ; Model level of maximum

     rmw_ind(0,it)   = (/rmw_arr(0,inds_s(0,0),it)/)	  ; Input values into large array
     rmw_ind(1,it)   = h0
     print("RMW ("+h0+" m) occurs at "+rmw_ind(0,it)+" degrees from centre")

     delete([/dims_s,spd_1d,inds_s/])    	 	  ; Tidy up after each iteration

   end do  	    ; End time loop (do it = ts0, tf0)

;=================================
; Write out values to text file 
;=================================

   print_clock("Writing out variables to text files...")
   diri_out = "$ar/text/ml5_"

 ; Write out RMW
   rmw_out0 = diri_out+dat+"_"+ens0+"_rmw_"+cn0+".txt"
   asciiwrite(rmw_out0,rmw_ind(0,:))

 ; Write out model level corresponding to RMW
   rmw_out1 = diri_out+dat+"_"+ens0+"_rmw_lev_"+cn0+".txt"
   asciiwrite(rmw_out1,rmw_ind(1,:))

 ; Write out maximum tangential wind 
   vtan_out = diri_out+dat+"_"+ens0+"_vtan_max_"+cn0+".txt"
   asciiwrite(vtan_out,vtan_max0(:)) 
   
end 