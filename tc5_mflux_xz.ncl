; Script to produce time-height plots of vertical mass flux 

; Run using:

; ncl dat=\"02T12\" opt=\"x11\" ens0=\"em11\" dist=3.0 ts0=539 tf0=767
; calc=1 cn0=\"slp\" ar=1.0 nr=21 nz=121 pos=1 tc5_mflux_xz.ncl

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'opt'   = output file format ("pdf" or "x11")
; 'ens0'  = ensemble member (em00 -- em11)
; 'dist'  = size of box following storm (degrees)
; 'ts0'   = starting time index (0-119)
; 'tf0'   = ending time index (0-119)
; 'calc'  = centre calculated on each level (0), calculated offline (1) or using 980 m vort (2)
; 'cn0'   = method of storm tracking ["slp", "geo", "geo_sm", "vort"]
; 'ar'    = radius of outermost radial circle on cylindrical grid (degrees)
; 'nr'    = number of radial circles on cylindrical grid
; 'nz'    = number of grid points in 360º (azimuthal spacing)
; 'pos'   = ignore (1) or include (0) negative values when calculating mass flux
 
; Load main NCL libraries ($NCARG_ROOT = /apps/applications/ncl/6.5.0/1/default)
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"

; Also load additional user-written functions
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/nc_levs.ncl"

begin

;=================
; Read in files 
;=================

; Find 'pc' and 'pd' files (model height level data)
  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0+"/vc"
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0
  fili_c1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pc.nc")
  fili_d1     =	systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pd.nc")
  fili_j1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pj.nc")
  fili_k1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pk.nc")

  fili_c      = diri+"/"+fili_c1
  fili_d      = diri+"/"+fili_d1
  fili_j      = diri+"/"+fili_j1
  fili_k      = diri+"/"+fili_k1

; Calculate number of times in each file using built-in function ('nc_times')         
  numINFO_c   = nc_times(fili_c)
  numINFO_d   = nc_times(fili_d)
  numINFO_j   = nc_times(fili_j)
  numINFO_k   = nc_times(fili_k)
  numINPUT_c  = numINFO_c[0]            ; Number of times in 'pc' stream              
  time_c      = numINFO_c[1]            ; Array of times from 'pc' stream file        
  numINPUT_d  = numINFO_d[0]            ; Number of times in 'pd' stream              
  time_d      = numINFO_d[1]            ; Array of times from 'pd' stream file        
  numINPUT_j  = numINFO_j[0]            ; Number of times in 'pj' stream              
  time_j      = numINFO_j[1]            ; Array of times from 'pj' stream file        
  numINPUT_k  = numINFO_k[0]            ; Number of times in 'pk' stream              
  time_k      = numINFO_k[1]            ; Array of times from 'pk' stream file 

; Calculate storm motion using built-in function ('st_rm')
; 'dat'     = initialisation time (02T12, 03T00, ...)
; 'ens0'    = ensemble simulation (em00, em01, ...)
; 'diri'    = path to input files (see above)
; 'dist'    = size of box following storm (degrees)
; 'mins'    = analyse 1-h (0) or 5-min (1) data

; Call external function 'st_rm'
  storm_rel = st_rm(dat, ens0, "$ar/text/", dist, 1)
  u_cyc     = storm_rel[0]		   ; Zonal wind 
  v_cyc     = storm_rel[1]		   ; Meridional wind
  vel_cyc   = storm_rel[2]		   ; Vector wind 
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

;==================================================================================
;;;;;;;;;;;;;;;;;;;;;;;;;;; Start of main part of script ;;;;;;;;;;;;;;;;;;;;;;;;;;
;==================================================================================

  print_clock("Starting main part of script...")

;========================================
; Define cylindrical coordinate arrays
;========================================

  pii        = 3.14159265
  radii      = fspan(0.0,ar,nr)      ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,nz)   ; Array of azimuth angles (0, 5, 10, etc)
  thetas_rad = (pii/180)*thetas      ; Degrees to radians
  rad_size   = dimsizes(radii)-1     ; Index for later calculation
  num_l      = 63      	   	     ; Number of model levels

; Define constants
  rd	= 287.0	  ; Gas constant for dry air (J/kg/K^2)
  cp0	= 1004.0  ; Specific heat of dry air at constant pressure (J/kg/K)
  kp	= 0.286   ; For Exner function calculations (0.286)
  pref	= 1000.0  ; Reference pressure p0 (1000 hPa)

;===========================
; Create necessary arrays 
;===========================

; Array to hold all 'xpos' and 'ypos' positions
  xpos_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),numINPUT_c/),float)
  ypos_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),numINPUT_c/),float)

; Arrays to hold lat/lon at centre of cylindrical grid
  xcen_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),numINPUT_c/),float)
  ycen_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),numINPUT_c/),float)

; Distances between points in cylindrical grid 
  rad_arr0  = new((/num_l,dimsizes(radii),dimsizes(thetas),numINPUT_c/),float)

;========================================
; Define 3D arrays for mean quantities 
;========================================

; Mean vertical velocity {w}
  mean_vvel0       = new((/num_l,dimsizes(radii),numINPUT_c/),float)
  mean_vvel0!0     = "lev"
  mean_vvel0!1     = "rad"
  mean_vvel0!2     = "time"
  mean_vvel0@description = "Azimuthally averaged vertical velocity"
  mean_vvel0@units = "m s~S~-1~N~"

; Mean mass flux {flux}
  mean_flux0       = new((/num_l,dimsizes(radii),numINPUT_c/),float)
  mean_flux0!0     = "lev"
  mean_flux0!1     = "rad"
  mean_flux0!2     = "time"
  mean_flux0@description = "Azimuthally averaged mass flux"
  mean_flux0@units = "kg m~S~-2~N~ s~S~-1~N~"

; Mean relative vorticity {zeta}
  mean_vort0       = new((/num_l,dimsizes(radii),numINPUT_c/),float)
  mean_vort0!0     = "lev"
  mean_vort0!1     = "rad"
  mean_vort0!2     = "time"
  mean_vort0@description = "Azimuthally averaged relative vorticity"
  mean_vort0@units = "s~S~-1~N~"

;=============================================
; Also define 4D arrays for mean quantities                                           
;=============================================                                        

; Mean vertical velocity {w}                                                     
  all_vvel0        = new((/num_l,dimsizes(radii),dimsizes(thetas),numINPUT_c/),float)
  all_vvel0!0      = "lev"
  all_vvel0!1      = "rad"
  all_vvel0!2      = "azi"
  all_vvel0!3      = "time"
  all_vvel0@description = "Azimuthally averaged vertical velocity"
  all_vvel0@units  = "m s~S~-1~N~"

; Mean mass flux {flux}                                                          
  all_flux0        = new((/num_l,dimsizes(radii),dimsizes(thetas),numINPUT_c/),float)
  all_flux0!0      = "lev"
  all_flux0!1      = "rad"
  all_flux0!2      = "azi"
  all_flux0!3      = "time"
  all_flux0@description = "Azimuthally averaged mass flux"
  all_flux0@units  = "kg m~S~-2~N~ s~S~-1~N~"

; Mean relative vorticity {zeta}                                                 
  all_vort0        = new((/num_l,dimsizes(radii),dimsizes(thetas),numINPUT_c/),float)
  all_vort0!0      = "lev"
  all_vort0!1      = "rad"
  all_vort0!2      = "azi"
  all_vort0!3      = "time"
  all_vort0@description = "Azimuthally averaged relative vorticity"
  all_vort0@units  = "s~S~-1~N~"

;=========================
; Start multiple loops
;=========================

  print_clock("Working on ensemble member: "+ens0)

  ct          = ts0                       ; Counter variable
  llbox       = toint(dist*50)  	   ; Size of lat/lon arrays (see line below)
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
  elseif (cn0 .eq. "slpf") then
   cn = 1
  end if 

;===============================
; Start loop over input files 
;===============================

  c  = addfile(fili_c,"r")         	 ; Read in 'pc' stream [u,v,vort]
  d  = addfile(fili_d,"r")         	 ; Read in 'pd' stream [w,p]
  j  = addfile(fili_j,"r")               ; Read in 'pj' stream [t,theta,mixing ratio]

  do it = ts0, tf0		; Loop over times in file

    print("Working on time: "+time_str(it))
    title_arr(it) = "Valid at "+time_str(it)+" (T+"+it+")"

;=====================================================
; Read in variables from 'pc' stream (model levels)
;=====================================================

    t0    = lat_arr(cn,it,0)
    t1    = lat_arr(cn,it,1)
    n0    = lon_arr(cn,it,0)
    n1    = lon_arr(cn,it,1)
    print_clock("Grid (pc stream) extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

  ; Read in basic variables from 'pc' stream
    lonC   = c->longitude({n0:n1})                ; longitude ['d1' grid points]
    latC   = c->latitude({t0:t1})                 ; latitude  ['d0' grid points]
    hybC   = c->hybrid_ht(:) 		          ; 63 model ('theta') levels

    u_varname = "x-wind"
    v_varname = "y-wind"

  ; Read in horizontal wind components 
  ; Both variables below -- [hybC | 63] * [latC | 150] * [lonC | 150]
    u0  = c->$u_varname$(it,:,{t0:t1},{n0:n1}) ; Zonal wind	
    v0  = c->$v_varname$(it,:,{t0:t1},{n0:n1}) ; Meridional wind 

  ; Calculate storm-relative winds
    u = u0 - u_cyc(cn,it)
    v = v0 - v_cyc(cn,it)

  ; Add metadata from original horizontal wind arrays
    copy_VarCoords(u0,u)
    copy_VarAtts(u0,u)

    copy_VarCoords(v0,v)
    copy_VarAtts(v0,v)

  ; Calculate relative vorticity using centered finite differences
  ; Option '2' --> boundary points estimated using one-sided difference scheme
    vort = uv2vr_cfd(u, v, latC, lonC, 2)
    copy_VarCoords(u, vort)
    vort@units = "s~S~-1~N~"
    vort@name = "Relative vorticity on model levels"

;======================================================
; Read in variables from 'pd' stream (model levels) 
;======================================================

  ; Read in basic variables from 'pd' stream 
    lon1D  = d->longitude_1({n0:n1})                ; longitude [xxx grid points]
    lat1D  = d->latitude_1({t0:t1})                 ; latitude  [xxx grid points]
    lonD   = d->longitude({n0:n1})                  ; longitude [xxx grid points]
    latD   = d->latitude({t0:t1})                   ; latitude  [xxx grid points]
    hybD   = d->hybrid_ht(:)			    ; 63 model ('rho') levels

  ; Read in temperature, pressure, vertical velocity and geopotential height
  ; All variables below -- [hybD | 63] * [latD | 150] * [lon1D | 151]
    p      = d->p(it,:,{t0:t1},{n0:n1})           ; Pressure (Pa)

;======================================                                            
; Read in variables from 'pj' stream                                               
;======================================                                            

  ; Read in basic variables from 'pj' stream                                       
    hybJ   = j->hybrid_ht(:)                        ; 63 model ('rho') levels      

  ; Read in temperature, potential temperature and mixing ratio                    
  ; All variables below  -- [hybJ | 63] * [latJ | 150]  * [lonJ | 150]             
    t      = j->temp(it,:,{t0:t1},{n0:n1})          ; Temperature                  

  ; Calculate dry air density from pressure and temperature
    rho	   = p / (rd * t)

  ; Add metadata
    copy_VarCoords(p, rho)
    rho@description = "Dry air density"
    rho@units = "kg m~S~-3~N~"

;======================================                                               
; Read in variables from 'pk' stream                                                  
;======================================                                               

  ; Read in basic variables from 'pk' stream                                          
    hybK   = k->hybrid_ht(:)                        ; 63 model ('rho') levels         

  ; Read in vertical velocity                                 
    w      = k->dz_dt(it,:,{t0:t1},{n0:n1})         ; Vertical velocity (m s-1)       

  ; Calculate vertical mass flux
    mflux = rho * w
    copy_VarCoords(w,mflux)
    copy_VarAtts(w,mflux)
    mflux@units         = "kg m~S~-2~N~ s~S~-1~N~"
    mflux@standard_name = "vertical_mass_flux"
    mflux@long_name     = "Vertical mass flux"
    mflux@title         = "Vertical mass flux"
    mflux@name          = "vert_mflux"
 
  ; Convert units of pressure before continuing
    p  	   = p / 100

  ; Add metadata
    p@units = "hPa"

    dx0     = lonD(1) - lonD(0) ; Grid spacing (longitude)
    dy0     = latD(1) - latD(0) ; Grid spacing (latitude)

;==================================================================
; Make sure all arrays are same size (make this into a function)
;==================================================================

  ; 2D array (8 * 3) -- in function, set dim-0 equal to number of variables we've read in
    size = (/dimsizes(u),dimsizes(v),dimsizes(vort),dimsizes(t),\
             dimsizes(p),dimsizes(w),dimsizes(z),dimsizes(mflux)/)

  ; Reduce size of any dimension larger than the rest 
    do sz  = 0, 1
     nind = min( size(:,sz) )
     if ( any (size(:,sz) .gt. nind) ) then
      size(:,sz) = nind
     end if
    end do    

  ; Resize arrays based upon results of calculations above
  ; 'size' = [8 * 3]
    u2    = u(:,0:size(0,1)-1,0:size(0,2)-1)
    v2    = v(:,0:size(1,1)-1,0:size(1,2)-1)
    vort2 = vort(:,0:size(2,1)-1,0:size(2,2)-1)
    t2    = t(:,0:size(3,1)-1,0:size(3,2)-1)
    p2    = p(:,0:size(4,1)-1,0:size(4,2)-1)
    w2    = w(:,0:size(6,1)-1,0:size(6,2)-1)
    z2    = z(:,0:size(5,1)-1,0:size(5,2)-1)
    flux2 = mflux(:,0:size(6,1)-1,0:size(6,2)-1)
    lat2C = latC(0:nind-1)
    lon2C = lonC(0:nind-1)
    lat2D = latD(0:nind-1)
    lon2D = lonD(0:nind-1)
    lat3D = lat1D(0:nind-1)
    lon3D = lon1D(0:nind-1)
    delete([/u,v,vort,t,p,z,w,mflux,latC,lonC,latD,lonD,lat1D,lon1D/])

  ; Reset names to avoid further code changes
    u     = u2
    v     = v2
    vort  = vort2
    t     = t2
    z     = z2
    p     = p2
    w     = w2
    mflux = flux2
    latC  = lat2C
    lonC  = lon2C
    latD  = lat2D
    lonD  = lon2D
    lat1D = lat3D
    lon1D = lon3D
    delete([/u2,v2,vort2,t2,z2,p2,w2,flux2,lat2C,lon2C,lat2D,lon2D,lat3D,lon3D/])

;==========================
; Loop over model levels 
;==========================

    do ilev = 0, num_l-1 ; Loop over chosen model levels

     hy = sprintf("%0.0f",hybD(ilev))
     print("Working on model level " + ilev+ " ("+hy+" m AGL)")

   ; Read in variables on model levels
     w_plane		    = w(ilev,:,:)	   ; Vertical velocity
     prs_plane		    = p(ilev,:,:)	   ; Pressure 
     vort_plane	   	    = vort(ilev,:,:)	   ; Relative vorticity
     flux_plane		    = mflux(ilev,:,:)	   ; Mass flux 
     
   ; Calculate location of vorticity maximum (storm centre) on model level
     vort_max      	    = max(vort_plane)
     vort_smth     	    = smth9_Wrap(vort_plane, 0.5, 0.5, True)
     vort_max_smth 	    = max(vort_smth)
   
   ; Find index of relative vorticity maximum
     dims          	    = dimsizes(vort_plane)                 ; Size of 'vort_plane'
     vort1d        	    = ndtooned(vort_plane)                 ; Make 1D array from 'vort_plane'
     inds          	    = ind_resolve(maxind(vort1d),dims)     ; Index of vort max
     vort1d_sm     	    = ndtooned(vort_smth)
     inds_sm       	    = ind_resolve(maxind(vort1d_sm),dims)  ; Index of smoothed vort max

   ;==================================================================
   ; Either calculate the centre, or read in pre-calculated values
   ;==================================================================

   if (calc .eq. 0) then

   ; Find latitude and longitude of storm centre using vorticity calculation above
     lat_max0 	            = t0 + (dy0 * inds(0,0))
     lon_max0 	   	    = n0 + (dx0 * inds(0,1))

     print("Centre (vort): "+lat_max0+" degrees N, "+lon_max0+" degrees E")

   ; Create smaller grid around location identified using vorticity maximum
     lt1 = lat_max0 - 0.15
     lt2 = lat_max0 + 0.15
     ln1 = lon_max0 - 0.15
     ln2 = lon_max0 + 0.15

   ; Calculate minimum pressure in region surrounding vorticity maximum 
     prs_plane0 = prs_plane({lt1:lt2},{ln1:ln2})       ; Smaller grid (1.0 x 1.0 deg)
     prs_min0   = min(prs_plane0)                      ; Find minimum on smaller grid
     dims_p0    = dimsizes(prs_plane0)                 ; Size of smaller grid
     prs_1d0    = ndtooned(prs_plane0)                 ; Create 1-D array
     inds_p0    = ind_resolve(minind(prs_1d0),dims_p0) ; Find index of minimum

     lat_max    = lt1 + (dy0 * inds_p0(0,0))
     lon_max    = ln1 + (dx0 * inds_p0(0,1))
     print("Centre (vort + pres): "+lat_max+" degrees N, "+lon_max+" degrees E")
     delete([/prs_plane0,prs_min0,dims_p0,prs_1d0,inds_p0/])

   ; Distance between centre calculated above using relative vorticity, and that calculated offline
     dlat       = abs(lat_max - centre(cn,it-1,0))
     dlon	= abs(lon_max - centre(cn,it-1,1))

   ; If our storm centre calculations deviate too much (> 0.75º) from the values
   ; we read in using the NCL function 'st_rm', calculate a new centre using
   ; the minimum pressure. Constrain the new calculation using the values from
   ; 'st_rm', to prevent spurious centres being identified and analysed. 

     if (dlat .gt. 1.00 .or. dlon .gt. 1.00) then 

      lt1 = centre(cn,it,0) - 0.25
      lt2 = centre(cn,it,0) + 0.25
      ln1 = centre(cn,it,1) - 0.25
      ln2 = centre(cn,it,1) + 0.25

      prs_plane0 = prs_plane({lt1:lt2},{ln1:ln2})       ; Smaller grid (1.0 x 1.0 deg)
      prs_min0   = min(prs_plane0)                      ; Find minimum on smaller grid
      dims_p0    = dimsizes(prs_plane0)                 ; Size of smaller grid
      prs_1d0    = ndtooned(prs_plane0)                 ; Create 1-D array
      inds_p0    = ind_resolve(minind(prs_1d0),dims_p0) ; Find index of minimum

      lat_max 	 = lt1 + (dy0 * inds_p0(0,0))
      lon_max 	 = ln1 + (dx0 * inds_p0(0,1))
      print("New centre: "+lat_max+" degrees N, "+lon_max+" degrees E")
      delete([/prs_plane0,prs_min0,dims_p0,prs_1d0,inds_p0/])

     end if 

   elseif (calc .eq. 1) then

      lat_max    = centre(cn,it,0)		; Does this relate to 'pc' or 'pd' stream???
      lon_max    = centre(cn,it,1)

   end if

   ;======================================================================
   ; Create lat/lon arrays for later switch to cylindrical coordinates
   ;======================================================================

   ;===============
   ; 'pc' stream
   ;===============
  
   ; Indices of lat/lon grid for 'u' and 'v' in 'pc' stream
     latC_plane  = ind(t0.le.latC.and.latC.le.t1)
     lonC_plane  = ind(n0.le.lonC.and.lonC.le.n1)

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
     latD_plane  = ind(t0.le.latD.and.latD.le.t1)
     lonD_plane  = ind(n0.le.lonD.and.lonD.le.n1)

   ; Array sizes
     latD_size   = dimsizes(latD_plane)
     lonD_size   = dimsizes(lonD_plane)

   ; Array of latitude/longitude points corresponding to the indices above
     latD_sub    = latD(latD_plane(0):latD_plane(latD_size-1))
     lonD_sub    = lonD(lonD_plane(0):lonD_plane(lonD_size-1))

   ; Indices of lat/lon grid for 'str_yz' in 'pd' stream
     lat1D_plane = ind(t0.le.lat1D.and.lat1D.le.t1)
     lon1D_plane = ind(n0.le.lon1D.and.lon1D.le.n1)

   ; Array sizes
     lat1D_size  = dimsizes(lat1D_plane)
     lon1D_size  = dimsizes(lon1D_plane)

   ; Array of latitude/longitude points corresponding to the indices above
     lat1D_sub   = lat1D(lat1D_plane(0):lat1D_plane(lat1D_size-1))
     lon1D_sub   = lon1D(lon1D_plane(0):lon1D_plane(lon1D_size-1))

     delete([/vort1d,inds,vort1d_sm,vort_smth,inds_sm/])

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
     flux_int = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Mass flux 

   ; Create cylindrincal coordinate array from original latitude/longitude grid
     do irad = 0, dimsizes(radii)-1
       r = radii(irad)
       do iang = 0, dimsizes(thetas)-1
          theta = thetas_rad(iang)   

           xpos(irad,iang) = lon_max + r * cos(theta)
 	   ypos(irad,iang) = lat_max + r * sin(theta)
           xcen(irad,iang) = lon_max
           ycen(irad,iang) = lat_max
	   
       end do
     end do

   ; Fill 'xpos_all' and 'ypos_all' with latitude & longitude values 
     xpos_all(ilev,:,:,it) = xpos(:,:)
     ypos_all(ilev,:,:,it) = ypos(:,:)
     xcen_all(ilev,:,:,it) = xcen(:,:)
     ycen_all(ilev,:,:,it) = ycen(:,:)

   ; Fill new arrays with interpolated values on the cylindrical grid 
   ; Which latitude/longitude grids are each of these variables defined on? 
   ; vort_plane	 -- latC,  lonC
   ; w_plane     -- latD,  lon1D
   ; flux_plane	 -- latD,  lon1D

     do iang = 0, dimsizes(thetas_rad)-1
      vort_int(:,iang) = linint2_points(lonC_sub,latC_sub,vort_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      w_int(:,iang)    = linint2_points(lon1D_sub,latD_sub,w_plane,False,\
                                          xpos(:,iang),ypos(:,iang),0)
      flux_int(:,iang) = linint2_points(lon1D_sub,latD_sub,flux_plane,False,\
					  xpos(:,iang),ypos(:,iang),0)
     end do

   ; Input into large arrays 
     all_vort0(ilev,:,:,ct)  = (/vort_int(:,:)/)
     all_vvel0(ilev,:,:,ct)  = (/w_int(:,:)/)

     if (pos .eq. 1) then 
      all_flux0(ilev,:,:,ct) = where( (/flux_int(:,:)/) .gt. 0, (/flux_int(:,:)/), 0)
     elseif (pos .eq. 0) then 
      all_flux0(ilev,:,:,ct) = (/flux_int(:,:)/)
     end if 

   ; Calculate mean (azimuthally-averaged) variables 
     do irad = 0, dimsizes(radii)-1

       mean_vvel0(ilev,irad,ct)    = avg(w_int(irad,:))      ; [ {w} ]
       mean_vort0(ilev,irad,ct)    = avg(vort_int(irad,:))   ; [ {zeta} ]
       mean_flux0(ilev,irad,ct)    = avg(flux_int(irad,:))   ; [ {mflux} ]

     end do

   ; Calculate radial distances between grid points 
     do irad = 1, dimsizes(radii)-1
      rp1 = min((/irad+1,dimsizes(radii)-1/))
      rm1 = max((/irad-1,0/))
      do th0 = 0, dimsizes(thetas)-1
       tp1 = min((/th0+1,dimsizes(thetas)-1/))
       tm1 = max((/th0-1,0/))
 
     ; Radial distance between each point and nearest neighbour. Arguments in 'gc_latlon':
     ; 'ypos(rm1,th0)' = latitude of first point
     ; 'xpos(rm1,th0)' = longitude of first point
     ; 'ypos(rp1,th0)' = latitude of second point 
     ; 'xpos(rp1,th0)' = longitude of second point
     ; '2'             = avoid extra interpolated points in between 
     ; '3'             = return the distance in metres
     ; See https://www.ncl.ucar.edu/Document/Functions/Built-in/gc_latlon.shtml for more info
       rad_arr0(ilev,irad,th0,ct) = gc_latlon(ypos(rm1,th0),xpos(rm1,th0),\
                                              ypos(rp1,th0),xpos(rp1,th0),2,3)
      end do
     end do

   ; Tidy up
     delete([/flux_int,w_int,vort_int/])
     delete([/w_plane,flux_plane,vort_plane/])
     delete([/xpos,ypos,xcen,ycen/])
     delete([/latC_plane,lonC_plane,latD_plane,lonD_plane,lon1D_plane,lat1D_plane/])
     delete([/latC_sub,lonC_sub,latD_sub,lonD_sub,lat1D_sub,lon1D_sub/])

    end do			    ; End model level loop (do ilev = 0, num_l-1)

   ct  = ct + 1  		    ; Counter variable (time)

 ; Tidy up before next iteration (all variables from 'pc' and 'pd' streams)
   delete([/lonC,latC,u,v,vort/])
   delete([/lonD,latD,lon1D,lat1D,t,p,w,z,rho,mflux/])

  end do     ; End time loop (do it = 0, times-1)

;========================================
; Calculate area integral of mass flux 
;========================================

  flux_sum  = all_flux0(:,:,:,:) * rad_arr0
  flux_int0 = new((/num_l,numTIMES/),float)

  do ilev = 0, num_l-1
   do it = ts0, tf0
    flux_int0(ilev,it) = sum(flux_sum(ilev,:,:,it) )
   end do
  end do 
  
;=================================================
; Tidy up and produce final arrays for plotting
;=================================================

; Add metadata 
  hyb_plot        = hybC / 1000
  hyb_plot@units  = "km"

  flux_int0!0     = "lev"
  flux_int0&lev   = hyb_plot
  flux_int0!1     = "time"
  flux_int0&time  = time_int

; Convert units for plotting (Kilroy et al. 2016)
  flux_int0       = flux_int0 / (10 ^ 7)

  printMinMax(flux_int0,False)

;==========================
; Options for plotting 
;==========================

  output = "$sam/nepartak/images/cyl_coords/mflux_xz_"+dat+"_"+ens0+"_"+ts0+"_"+tf0+"_"+cn0
  wks    = gsn_open_wks(opt,output)

  panel  = new(1,graphic)
  gsn_define_colormap(wks,"prcp_new")

; Mass flux resources 
  opts_flux                             = True
  opts_flux@cnFillOn                    = True
  opts_flux@cnLineLabelInterval         = 2.0
  opts_flux@cnLineLabelFontHeightF      = 0.012
  opts_flux@cnLineLabelBackgroundColor  = "transparent"
  opts_flux@cnLineLabelPlacementMode    = "constant"
  opts_flux@cnLinesOn                   = False
  opts_flux@cnInfoLabelOn               = False
  opts_flux@cnLevelSelectionMode        = "ExplicitLevels"

  opts_flux@cnLevels                    = (/0.0, 0.1, 0.2, 0.3, 0.4, 0.5, \
                                            0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8/)
  opts_flux@cnFillColors                = (/0,3,4,5,6,7,8,9,10,11,12,13,14,15/)
  opts_flux@gsnPaperOrientation         = "landscape"
  opts_flux@tiMainFontHeightF           = 0.0125
  opts_flux@gsnLeftString               = ""
  opts_flux@gsnRightString              = ""
  opts_flux@gsnMaximize                 = False
  opts_flux@lbLabelBarOn                = False
  opts_flux@gsnDraw                     = False
  opts_flux@gsnFrame                    = False
  opts_flux@gsnAddCyclic                = False

; Y-axis options
  opts_flux@tiYAxisString               = "Height (km)"
  opts_flux@tmYLMode                    = "Explicit"
  opts_flux@tmYROn                      = "False" ; Turn off tickmarks (right)
  opts_flux@gsnYAxisIrregular2Linear    = True 
  opts_flux@trYMinF                     = 0.0     ; Min
  opts_flux@trYMaxF                     = 16.0    ; Max
  dy_axes                               = toint(opts_flux@trYMaxF - opts_flux@trYMinF)

; Set y-axis labels (height in km)
  opts_flux@tmYLValues                  = fspan(opts_flux@trYMinF, opts_flux@trYMaxF, dy_axes)
  opts_flux@tmYLLabels                  = sprintf("%0.0f",opts_flux@tmYLValues)

; X-axis options
  opts_flux@tiXAxisString               = "Time (h)"
  opts_flux@tmXBMode                    = "Explicit"
  opts_flux@tmXTOn                      = "False" ; Turn off top tickmarks (top)
  opts_flux@trXMinF                     = ts0
  opts_flux@trXMaxF                     = tf0

; Create labels based on length of integration                                             
  tt0 = tf0 - ts0
  
  if (tt0 .eq. 6) then 
   opts_flux@tmXBValues                 = (/ts0, ts0+3, tf0/)
  elseif (tt0 .eq. 12) then 
   opts_flux@tmXBValues                 = (/ts0, ts0+6, tf0/)  
  elseif (tt0 .eq. 24) then
   opts_flux@tmXBValues                 = (/ts0, ts0+12, tf0/)
  elseif (tt0 .eq. 36) then
   opts_flux@tmXBValues                 = (/ts0, ts0+12, ts0+24, tf0/)
  elseif (tt0 .eq. 48)
   opts_flux@tmXBValues                 = (/ts0, ts0+12, ts0+24, ts0+36, tf0/)
  elseif (tt0 .eq. 60)
   opts_flux@tmXBValues                 = (/ts0, ts0+12, ts0+24, ts0+36, ts0+48, tf0/)
  elseif (tt0 .eq. 72)
   opts_flux@tmXBValues                 = (/ts0, ts0+12, ts0+24, ts0+36, \
                                             ts0+48, ts0+60, tf0/)
  elseif (tt0 .eq. 84)
   opts_flux@tmXBValues                 = (/ts0, ts0+12, ts0+24, ts0+36, ts0+48, \
                                             ts0+60, ts0+72, tf0/)
  elseif (tt0 .eq. 96)
   opts_flux@tmXBValues                 = (/ts0, ts0+12, ts0+24, ts0+36, ts0+48, \
                                             ts0+60, ts0+72, ts0+84, tf0/)
  end if

; Create x-axis labels from values in 'time_str'
  opts_flux@tmXBLabels = time_str(opts_flux@tmXBValues)

;================================
; Produce filled contour plots 
;================================

  panel(0) = gsn_csm_contour(wks, flux_int0, opts_flux)

; Panel plot resources 
  optsP                  = True 
  optsP@gsnFrame         = False
  optsP@gsnPanelLabelBar = False
  optsP@txString         = ""
  optsP@gsnMaximize      = True
  optsP@gsnPanelTop      = 0.95
  optsP@gsnPanelBottom   = 0.05

  gsn_panel(wks, panel, (/1,1/), optsP)
  frame(wks)

end