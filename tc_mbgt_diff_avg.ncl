; Script to read in the diffusivity and smooth out high-frequency fluctuations

; Run using:

; ncl dat=\"02T12\" opt=\"x11\" ens0=\"em11\" dist=3.0 ts0=599 tf0=743 clr=5
; calc=1 full=1 rad0=1 ave=0 wrt=0 new0=0 tc_mbgt_diff_avg.ncl

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'opt'   = output file format ("pdf" or "x11")
; 'ens0'  = ensemble member (em00 -- em11)
; 'dist'  = size of box following storm (degrees)
; 'ts0'   = starting time index (0-119)
; 'tf0'   = ending time index (0-119)
; 'clr'   = colour map for plots ("ncl_default" [1];  "bdr_extra" [2]; "amwg256" [3])
; 'calc'  = centre calculated on each level (0), calculated offline (1) or using 980 m vort (2)
; 'full'  = plot 13 panels inc. u,v,w (2), 10 panels w/mbgt terms only (1) or 4 panels (0)
; 'rad0'  = radius of final plots: 150 km (0) or 200 km (1)
; 'ave'   = calculate 3/6/12-h average (1) rather than analyse individual time interval (0)
; 'wrt'   = write out averages to text file (1) or leave alone (0 – default)
; 'new0'  = reduced (1) or original (0) scale for plotting u,v,w,AAM

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

; Find 'pc' and 'pd' files (model height level data)
  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0+"/vc"
  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0

  fili_k1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pk.nc")
  fili_k      = diri+"/"+fili_k1

; Calculate number of times in each file using built-in function ('nc_times')
  numINFO_k   = nc_times(fili_k)
  numINPUT_k  = numINFO_k[0]            ; Number of times in 'pk' stream
  time_k      = numINFO_k[1]            ; Array of times from 'pk' stream file

  print("There are "+numINPUT_k+" output time intervals in this file!")

; Call external function 'st_rm'
  storm_rel   = st_rm(dat, ens0, "$ar/text/", dist, 1)
  u_cyc       = storm_rel[0]                 ; Zonal wind
  v_cyc       = storm_rel[1]                 ; Meridional wind
  vel_cyc     = storm_rel[2]                 ; Vector wind
  numTIMES    = storm_rel[3]                 ; Number of times in file
  lat_arr     = storm_rel[4]                 ; Storm-following domain subset (lat)
  lon_arr     = storm_rel[5]                 ; Storm-following domain subset (lon)
  centre      = storm_rel[6]                 ; Storm centre at each hour (lat/lon)

; Create arrays for data and output information
  llbox       = toint(dist*50)               ; Calculate domain size (grid points)
  dsize       = (/llbox,llbox/)              ; Domain size (grid points)
  d0          = dsize(0)
  d1          = dsize(1)
  ct          = ts0                          ; Counter variable
  nh	      = (tf0 - ts0) / 12	     ; Number of hours analysed (for averaging later)
  num_l       = 63                    	     ; Number of model levels

  print_clock("Working on ensemble member: "+ens0)

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 1000000000
  end setvalues

;================================================================================
; Create new array to hold fixed values of storm centre position for each hour
;================================================================================

; Copy relevant arrays
  centre_fix  = centre
  ltarr_fix   = lat_arr
  lnarr_fix   = lon_arr 

; Fill 'centre' array with 'hourly' values rather than every 5 minutes
  do n = 1, 82
   n0 = (n * 12) - 1
   n1 = n0 + 11
   n2 = n0 + 5
   centre_fix(n0:n1,0) = centre(n2,0)
   centre_fix(n0:n1,1) = centre(n2,1)
  end do 
  delete([/n0,n1,n2/])

; Fill new 'lat/lon_arr' arrays with domain subset boundaries
  ltarr_fix(:,0) = centre_fix(:,0) - dist
  ltarr_fix(:,1) = centre_fix(:,0) + dist
  lnarr_fix(:,0) = centre_fix(:,1) - dist
  lnarr_fix(:,1) = centre_fix(:,1) + dist

;==========================================================================
; Create date/time string arrays for all times in file (use 'pk' stream)
;==========================================================================

  month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
                    "Jul","Aug","Sep","Oct","Nov","Dec"/)

  utc_date = cd_calendar(time_k, 0)

  year   = tointeger(utc_date(:,0))
  month  = tointeger(utc_date(:,1))
  day    = tointeger(utc_date(:,2))
  hour   = tointeger(utc_date(:,3))
  minute = tointeger(utc_date(:,4))
  second = utc_date(:,5)

; Create array to hold finished date strings
  hr_min = new(numINPUT_k,string)
  min0   = (/4,9,14,19,24,29,34,39,44,49,54,59/)      ; Minute array (before)
  min1   = (/5,10,15,20,25,30,35,40,45,50,55,0/)      ; Minute array (after)

; Before correcting, find the indices of the elements where 'minute' = 59
  hr_ind = ind(minute .eq. 59)

; Correct for errors in the code (round up values of 'minute' to multiples of 5)
  do i = 0, numINPUT_k-1
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
  do i = 0, numINPUT_k-1

   if (hour(i) .lt. 10) then                    ; Add zeros in front of single-digit hours

    if (minute(i) .lt. 10) then
     hr_min(i) = "0"+hour(i)+"0"+mins(i)
    elseif (minute(i) .ge. 10) then
     hr_min(i) = "0"+hour(i)+mins(i)
    end if

   else

    if (minute(i) .lt. 10) then
     hr_min(i) = hour(i)+"0"+mins(i)            ; Also add zeros in front of single-digit minutes
    elseif (minute(i) .ge. 10) then
     hr_min(i) = hour(i)+mins(i)
    end if

   end if

  end do
  delete(i)

  time_str  = new(numINPUT_k,string)
  time_arr  = new(numINPUT_k,string)
  title_arr = new(numINPUT_k,string)

  time_str  = hr_min + " UTC " + sprinti("%0.2i ", day) + month_abbr(month)
  time_arr  = sprinti("%0.2i", day) + month_abbr(month) + "_" + hr_min + "Z"
  title_arr = "Valid at "+time_str

;==================================================================================
;;;;;;;;;;;;;;;;;;;;;;;;;;; Start of main part of script ;;;;;;;;;;;;;;;;;;;;;;;;;;
;==================================================================================

; LHS (1) : d{v}/dt          - Local tendency of the mean tangential wind

; RHS (1) : -{u} * {f+vort}  - Mean radial influx of absolute vertical vorticity
; RHS (2) : -{w} * (d{v}/dz) - Mean vertical advection of mean tangential momentum
; RHS (3) : -{u' * vort'}    - Eddy radial vorticity flux
; RHS (4) : -{w' * (dv'/dz)} - Vertical advection of eddy tangential momentum
; RHS (5) : {D_v}            - Combined diffusive and PBL tendency

;===============================================
; Array to hold key terms over domain subsets
;===============================================

; Values every 5 minutes 
  taux_arr0  = new((/numINPUT_k,num_l,d0,d1/),float)
  tauy_arr0  = new((/numINPUT_k,num_l,d0,d1/),float) 
  kmh_arr0   = new((/numINPUT_k,num_l,d0,d1/),float)

; Hourly averages, output every 5 minutes
  taux_arr   = new((/numINPUT_k,num_l,d0,d1/),float)
  tauy_arr   = new((/numINPUT_k,num_l,d0,d1/),float)
  kmh_arr    = new((/numINPUT_k,num_l,d0,d1/),float)

;==================================
; Get the variables we will need
;==================================

  print_clock("Starting main part of script...")	

  k     = addfile(fili_k,"r")              ; Read in 'pk' stream [tau,PV,w]

  do it = ts0, tf0			   ; Loop over times in file

   print("Working on time: "+time_str(it))

    t0    = ltarr_fix(it,0)
    t1    = ltarr_fix(it,1)
    n0    = lnarr_fix(it,0)
    n1    = lnarr_fix(it,1)
    print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

;======================================
; Read in variables from 'pk' stream
;======================================

  ; Read in basic variables from 'pk' stream
    hybK   = k->hybrid_ht(:)                        ; 63 model ('rho') levels

  ; Alternative latitude/longitude grid
    lon1   = k->longitude_1({n0:n1})                ; longitude
    lat1   = k->latitude_1({t0:t1})                 ; latitude

  ; Read in diffusivity and wind stress
  ; 'str_xz' -- [hybD | 63] * [latD | 150]  * [lonD | 151]
  ; 'str_yz' -- [hybD | 63] * [lat1D | 150] * [lon1D | 151]
  ; 'kmh'    -- [hybD | 63] * [latD | 150]  * [lon1D | 151]
    str_xz = k->taux(it,:,{t0:t1},{n0:n1})          ; Wind stress, x-comp (N m/2)
    str_yz = k->tauy(it,:,{t0:t1},{n0:n1})          ; Wind stress, y-comp (N m/2)
    kmh    = k->unspecified(it,:,{t0:t1},{n0:n1})   ; Horizontal diffusivity (m^2 s-1)

  ; Fill pre-existing arrays with subset of each diagnostic
    taux_arr0(it,:,:,:)   = str_xz
    tauy_arr0(it,:,:,:)   = str_yz
    kmh_arr0(it,:,:,:)    = kmh

  end do

;===============================================================
; Now loop through these arrays and calculate hourly averages 
;===============================================================

  do it = 0, nh-1			; Loop over number of hours 

    ts = ts0 + (it * nh)
    tf = ts + 11

    taux_arr(ts:tf,:,:,:) = dim_avg_n_Wrap( taux_arr0(ts:tf,:,:,:), 0)
    tauy_arr(ts:tf,:,:,:) = dim_avg_n_Wrap( tauy_arr0(ts:tf,:,:,:), 0)
    kmh_arr(ts:tf,:,:,:)  = dim_avg_n_Wrap( kmh_arr0(ts:tf,:,:,:), 0)
  
  end do 

end 