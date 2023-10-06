; Script to compute the storm track for given TCs

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_rm.ncl"

; ncl dat=\"02T12\" int=1 opt=\"x11\" ens0=\"em11\" dist=3.0 tc_fc_times.ncl

begin

;==============================================
; Find 'pb' stream data to read in and plot 
;==============================================

  diri        = "/nfs/a319/earshar/"+dat+"/"+ens0+"/vc"

  fili_prefix = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0
  fili_b1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pb.nc")
  fili_b      = diri+"/"+fili_b1

; Get number of times in file using built-in function ('nc_times')
  numINFO     = nc_times(fili_b)
  numTIMES    = numINFO[0]		; Number of times
  time	      = numINFO[1]		; Array containing timestamp info at each time]
  print("Total number of time intervals = "+numTIMES)

;============================================================
; Calculate storm motion using built-in function ('st_rm')
;============================================================

; 'ens0'    = ensemble simulation
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
  llbox     = toint(dist*50)                   ; Calculate domain size (grid points)
  dsize     = (/llbox,llbox/)                  ; Domain size (grid points)
 
; Extract string for each ensemble member
  str1      = str_split(fili_b(0),"_")
  ens_str   = str1(5)

  print("Working on ensemble member "+ens_str)

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 10000000000
  end setvalues

;=====================================================
; Create correct date strings for each output time
;=====================================================

  month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
                    "Jul","Aug","Sep","Oct","Nov","Dec"/)

  utc_date = cd_calendar(time, 0)

  year   = tointeger(utc_date(:,0))
  month  = tointeger(utc_date(:,1))
  day    = tointeger(utc_date(:,2))
  hour   = tointeger(utc_date(:,3))
  minute = tointeger(utc_date(:,4))
  second = utc_date(:,5)

; Create array to hold finished date strings
  hr_min = new(numTIMES,string)

  min0   = (/4,9,14,19,24,29,34,39,44,49,54,59/)      ; Minute array (before)
  min1   = (/5,10,15,20,25,30,35,40,45,50,55,0/)      ; Minute array (after)

; Before correcting, find the indices of the elements where 'minute' = 59
  hr_ind = ind(minute .eq. 59)

; Correct for errors in the code (round up values of 'minute' to multiples of 5)
  ind_f = numTIMES-1

  do i = 0, ind_f
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

; Create string array of all values in 'minute'
  mins = tostring(minute)

; Create 'hhmm' strings for each time (for output and title strings)
  do i = 0, ind_f

   if (hour(i) .lt. 10) then                  ; Add zeros in front of single-digit hours

    if (minute(i) .lt. 10) then
     hr_min(i) = "0"+hour(i)+"0"+mins(i)
    elseif (minute(i) .ge. 10) then
     hr_min(i) = "0"+hour(i)+mins(i)
    end if

   else

    if (minute(i) .lt. 10) then
     hr_min(i) = hour(i)+"0"+mins(i)          ; Also add zeros in front of single-digit minutes
    elseif (minute(i) .ge. 10) then
     hr_min(i) = hour(i)+mins(i)
    end if

   end if

  end do

  time_str   = new(numTIMES,string)
  time_arr   = new(numTIMES,string)
  title_arr  = new(numTIMES,string)

  time_str   = hr_min + " UTC " + sprinti("%0.2i ", day) + month_abbr(month)
  time_arr   = sprinti("%0.2i", day) + month_abbr(month) + "_" + hr_min + "Z"
  title_arr  = "Valid at "+time_str

  print(time_str)
  print(time_arr)

end