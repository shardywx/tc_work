; Script to write out max. 10-m windspeed and min. MSLP following a TC 

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

; Find number of relevant input files in specific directory

  diri	      = "$sam/um/gl/ens/"
  fili_prefix = "20160704T0000Z_engl_"+ens+"_pb"
  fili_b1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"*.nc")
  fili_b      = diri+fili_b1
  numINPUT    = dimsizes(fili_b)
  numTIMES    = (numINPUT-1) * 6            ; Total number of times in all files [120]

  centre      = new((/numTIMES,2/),"float") ; Array to hold lat/lon at each time interval
  time_arr    = new(numTIMES,string)        ; Array to hold strings for second loop
  title_arr   = new(numTIMES,string)        ; Array to hold strings for plot titles
  ct	      = 0			    ; Counter variable

  spd_box     = new(numTIMES,float)	    ; Array to hold windspeed values 

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 10000000000
  end setvalues

;==========================
; Loop over input files
;==========================

  do nf = 0, numINPUT-1        ; 12 ensemble members in directory
  f = addfile(fili_b(nf),"r")

 ; Extract string for each ensemble member
   str1    = str_split(fili_b(nf),"_")
   ens_str = str1(3)

   print_clock("Working on file number "+nf)

;==================================
; Get the variables we will need
;==================================

   time  = f->t            ; Times in file (6 - every hour)
   times = dimsizes(time)

   lon   = f->longitude    ; longitude (800 points --> 0.225 to 359.775 degrees E)
   lat   = f->latitude     ; latitude  (600 points --> -89.85 to 89.85 degrees N)
   lon1  = f->longitude_1  ; longitude (800 points --> 0.0 to 359.55 degrees E)
   lat1  = f->latitude_1   ; latitude (601 points --> -90.0 to 90.0 degrees N)

   pres  = f->p            ; Pressure (levels)
   plevs = dimsizes(pres)  ; Size of pressure level array
   pres@units = "hPa"

   lat_0 = 6.05
   lat_1 = 30.0
   lon_0 = 110.0
   lon_1 = 150.0

   u     = f->u(:,:,{lat_0:lat_1},{lon_0:lon_1})    ; Zonal windspeed (m/s)
   v     = f->v(:,:,{lat_0:lat_1},{lon_0:lon_1})    ; Meridional windspeed (m/s)
   vort  = f->rvor(:,:,{lat_0:lat_1},{lon_0:lon_1}) ; Relative vorticity (/s)
   z     = f->ht(:,:,{lat_0:lat_1},{lon_0:lon_1})   ; Geopotential height (m)

   z     = z/10
   z@units = "dam"         ; Convert to decametres
   vort    = vort * (10 ^ -6) ; Convert to /s

   dx    = lon(1) - lon(0) 			    ; Grid spacing (longitude)
   dy    = lat(1) - lat(0) 			    ; Grid spacing (latitude)

 ; Pressure levels (pres)
 ; (0) 1000, (1) 950, (2) 900, (3) 850, (4) 800, (5) 750, (6) 700, (7) 650
 ; (8) 600, (9) 550, (10) 500, (11) 450, (12) 400, (13) 350
 ; (14) 300, (15) 250, (16) 200, (17) 150, (18) 100

;=====================================================
; Create correct date strings for each output time
;=====================================================

   month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
                     "Jul","Aug","Sep","Oct","Nov","Dec"/)

   times = dimsizes(time) ; Files are not all same size
   utc_date = cd_calendar(time, 0)

   year   = tointeger(utc_date(:,0))
   month  = tointeger(utc_date(:,1))
   day    = tointeger(utc_date(:,2))
   hour   = tointeger(utc_date(:,3))
   minute = tointeger(utc_date(:,4))
   second = utc_date(:,5)

 ; Correct for errors in the code (round up value of hour when minutes = 59)
   do it = 0, times-1
     if (minute(it).gt.30) then
       hour(it) = hour(it)+1
     end if
   end do

   date_str = new(times,string)
   time_str = new(times,string)
   out_str  = new(times,string)

;===========================
; Loop over times in file
;===========================

   do it = 0, times-1

     date_str(it) = sprinti("%0.2iUTC ", hour(it)) + \
                    sprinti("%0.2i ", day(it)) \
                    + month_abbr(month(it))
     time_str(it) = sprinti("%0.2i UTC ", hour(it)) + \
                    sprinti("%0.2i ", day(it)) \
                    + month_abbr(month(it))
     out_str(it)  = sprinti("%0.2i", day(it)) + \
                    month_abbr(month(it)) + \
                    "_" + sprinti("%0.2iZ", hour(it))

     time_arr(ct) = out_str(it)

     print_clock("Working on time: "+time_str(it))
     title_arr(ct) = "Valid at: "+time_str(it)

;================================================================================
; Find cyclone centre on subset of global grid (use both height and vorticity)
;================================================================================

     ilev	 = 2
     p		 = pres(ilev)
     print("Working on "+p+" hPa")

     u_plane     = u(it,ilev,:,:)
     v_plane     = v(it,ilev,:,:)
     vort_plane  = vort(it,ilev,:,:)
     geo_plane	 = z(it,ilev,:,:)
     dsize	 = dimsizes(vort_plane)

   ; RELATIVE VORTICITY
     vort_max      = max(vort_plane)
     vort_smth     = smth9_Wrap(vort_plane, 0.5, 0.5, True)
     vort_max_smth = max(vort_smth)

     dims          = dimsizes(vort_plane)
     vort1d        = ndtooned(vort_plane)
     inds          = ind_resolve(maxind(vort1d),dims)
     vort1d_sm     = ndtooned(vort_smth)
     inds_sm       = ind_resolve(maxind(vort1d_sm),dims)

     lat_max950 = lat_0 + (dy * inds(0,0)) ; Latitude of max. vorticity
     lon_max950 = lon_0 + (dx * inds(0,1)) ; Longitude of max. vorticity
     print("Cyclone centre (vort): "+lat_max950+" degrees N, "+lon_max950+" degrees E")

   ; GEOPOTENTIAL HEIGHT
     geo_min        = min(geo_plane)
     geo_smth       = smth9_Wrap(geo_plane, 0.5, 0.5, True)
     geo_min_smth   = min(geo_smth)

     dims_h         = dimsizes(geo_plane)
     geo1d          = ndtooned(geo_plane)
     inds_h         = ind_resolve(minind(geo1d),dims_h)
     geo1d_sm       = ndtooned(geo_smth)
     inds_h_sm      = ind_resolve(minind(geo1d_sm),dims_h)

     lat_min950 = lat_0 + (dy * inds_h(0,0))
     lon_min950 = lon_0 + (dx * inds_h(0,1))
     print("Cyclone centre (geo): "+lat_min950+" degrees N, "+lon_min950+" degrees E")

   ;========================================================
   ; EMPLOY A SAFETY NET IF WE IDENTIFY THE WRONG CYCLONE
   ;========================================================
     if (ct.eq.0) then
      centre(ct,0) = lat_min950
      centre(ct,1) = lon_min950
      print("First time in file!")

     else
      lt0  = centre(ct-1,0)             ; Previous TC latitude
      ln0  = centre(ct-1,1)             ; Previous TC longitude
      print("Previous TC latitude = "+lt0+" ; Previous TC longitude = "+ln0)
      d_lat = abs(lat_min950-lt0)       ; Change in latitude (t1 - t0)
      d_lon = abs(lon_min950-ln0)       ; Change in longitude
      print("Change in latitude = "+d_lat+" ; Change in longitude = "+d_lon)

      if (d_lat.gt.dist.or.d_lon.gt.dist) then     ; Recalculate TC centre if incorrect

     ; Create smaller grid [0.75 degrees]
       lt1 = lt0 - 0.75
       lt2 = lt0 + 0.75
       ln1 = ln0 - 0.75
       ln2 = ln0 + 0.75

       geo_plane0 = z(it,ilev,{lt1:lt2},{ln1:ln2})       ; Smaller grid (1.0 x 1.0 deg)
       geo_min0   = min(geo_plane0)                      ; Find minimum on smaller grid
       dims_h0    = dimsizes(geo_plane0)                 ; Size of smaller grid
       geo_1d0    = ndtooned(geo_plane0)                 ; Create 1-D array
       inds_h0    = ind_resolve(minind(geo_1d0),dims_h0) ; Find index of minimum

       lat_min950 = lt1 + (dy * inds_h0(0,0))
       lon_min950 = ln1 + (dy * inds_h0(0,1))
       print("New centre (hgt): "+lat_min950+" degrees N, "+lon_min950+" degrees E")
       centre(ct,0) = lat_min950
       centre(ct,1) = lon_min950

       delete([/geo_plane0, geo_min0, dims_h0, geo_1d0, inds_h0/])

      else
       print("Cyclone centre: "+lat_min950+" degrees N, "+lon_min950+" degrees E")
       centre(ct,0) = lat_min950
       centre(ct,1) = lon_min950
      end if

     end if

    ; Calculate 10-m wind speed
      speed = sqrt(u_plane^2.0 + v_plane^2.0)
      copy_VarCoords(u_plane, speed)
      copy_VarAtts(u_plane, speed)
      copy_VarMeta(u_plane, speed)
      speed@units = "m s~S~-1~N~"

    ; Tidy up before looping over pressure levels
      delete([/vort1d, inds, vort1d_sm, vort_smth, inds_sm/])
      delete([/geo1d, inds_h, geo1d_sm, geo_smth, inds_h_sm/])
      delete([/u_plane, v_plane, vort_plane, geo_plane/])

;================================
; Add min/max values to arrays 
;================================

      rad = 1
 
      spd_box(ct) = max(speed({lat_min950-rad:lat_min950+rad},\
                              {lon_min950-rad:lon_min950+rad}))

      print("900 hPa windspeed = "+spd_box(ct))

      ct = ct + 1 ; Counter variable (time)

      end do      ; End time loop (do it = 0, times-1)

    ; Tidy up before next iteration
      delete([/time,times,u,v,z,vort/])
      delete([/year,month,day,hour,minute,second,date_str,time_str,out_str,utc_date/])

    end do        ; End input file loop (do nf = 0, numINPUT-1)

;=================================
; Write out values to text file 
;=================================

    diri = "$sam/nepartak/images/cart_coords/gl_um/ens/gl_"+ens+"_"

  ; Windspeed
    spd_out  = diri+"uv"+p+".txt"
    asciiwrite(spd_out,spd_box)

end 