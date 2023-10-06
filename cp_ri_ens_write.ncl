; Script to compute the storm track for given TCs

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

; Find number of relevant input files in specific directory
  diri	      = "/nfs/a299/TCs/neptark/reg_netcdf"
  fili_prefix = "20160704T0000Z_ra1t_em"
  fili_a1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"*_pa.nc")
  fili_a      = diri+"/"+fili_a1
  numINPUT    = dimsizes(fili_a)
  centre      = new((/121,2/),"float")	; Array to hold lat/lon at each time interval

; Also add 'pb' stream files for vorticity information
  fili_b1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"*_pb.nc")
  fili_b      = diri+"/"+fili_b1
  numINPUT_b  = dimsizes(fili_b) 

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 10000000000
  end setvalues

;==========================
; Loop over input files
;==========================

  do nf = 0, numINPUT-1     ; 12 ensemble members in directory
  f = addfile(fili_a(nf),"r")

 ; Extract string for each ensemble member
   str1    = str_split(fili_a(nf),"_")
   ens_str = str1(3)

   print("Working on ensemble member "+ens_str)

;==================================
; Get the variables we will need
;==================================

 ;================================================
 ; Number of times in different variable arrays
 ;================================================

   time  = f->t            ; Times in file [119]
   times = dimsizes(time)

   time1 = f->t_1          ; Times in file [120]
   t1    = dimsizes(time1)

   time2 = f->t_2          ; Times in file [120]
   t2    = dimsizes(time2)

   time3 = f->t_3	   ; Times in file [120]
   t3	 = dimsizes(time3)

   time4 = f->t_4	   ; Times in file [121]
   t4	  = dimsizes(time4)

 ;=======================
 ; All other variables 
 ;=======================

   lon   = f->longitude    ; longitude (1098 points --> 109.04 to 152.92 degrees E)
   lat   = f->latitude     ; latitude  (810 points --> 1.8 to 34.16 degrees N)
   lon1  = f->longitude_1  ; longitude (1098 points --> 109.02 to 152.90 degrees E)
   lat1  = f->latitude_1   ; latitude (811 points --> 1.78 to 34.18 degrees N)

   u_varname = "x-wind"
   v_varname = "y-wind"

   u     = f->$u_varname$  ; 10-m zonal velocity (m/s) [time4 | z | lat_1 | lon_1]
   v     = f->$v_varname$  ; 10-m meridional velocity (m/s) [time4 | z | lat_1 | lon_1]
   t     = f->temp         ; Temperature (K) [time4 | z | lat | lon]
   rh    = f->rh           ; Relative humidity (%) [time4 | z | lat | lon]

   ll_ht = f->ll_ht        ; Cloud base height (kft) [time2 | z | lat | lon]
   omg   = f->dz_dt        ; Vertical velocity (m/s) [time2 | z | lat | lon]
   mslp  = f->p            ; Mean sea level pressure (Pa) [time2 | z | lat | lon]

   mslp  = mslp / 100      ; Convert to hPa
   mslp@units = "hPa"

   dx    = lon(1) - lon(0) ; Grid spacing (longitude)
   dy    = lat(1) - lat(0) ; Grid spacing (latitude)

  ; Pressure levels (pres)
  ; (0)   1000
  ; (1)   950
  ; (2)   900
  ; (3)   850
  ; (4)   800
  ; (5)   750
  ; (6)   700
  ; (7)   650
  ; (8)   600
  ; (9)   550
  ; (10)  500
  ; (11)  450
  ; (12)  400
  ; (13)  350
  ; (14)  300
  ; (15)  250
  ; (16)  200
  ; (17)  150
  ; (18)  100

  ; Degrees to radians (for calculation below)
    pii = 3.14159265
    lat_rad = (pii/180)*lat    ; u,v grid
    lon_rad = (pii/180)*lon

    lat1_rad = (pii/180)*lat1  ; height, temperature grid
    lon1_rad = (pii/180)*lon1

;=============================================================
; Create arrays to hold values of min. MSLP/ max. windspeed 
;=============================================================

 ; Arrays to hold values at each time interval
   slp_box  = new((t2/int)+1,float)
   spd_box  = new((t2/int)+1,float)
   
   ct       = 0 ; Counter variable

;===========================
; Loop over times in file
;===========================

     do	it = 0,	t2, int        ; Loop over times in file

    ; Output file string (time interval)
      intA  = sprinti("%0.1i",int)+"h"

      print("Working on forecast hour: "+it)
      title = it+" h forecast, initialised at 00 UTC 4 July"

;================================================================================
; Find cyclone centre on subset of global grid (use both height and vorticity)
;================================================================================

    ; Account for different sizes of arrays [u,v one element larger than MSLP]
      if (it.eq.0) then
        u_plane   = u(it,0,:,:)
        v_plane   = v(it,0,:,:)
        slp_plane = mslp(it,0,:,:)
      else
        u_plane    = u(it,0,:,:)        ; Zonal velocity            (m/s)
        v_plane    = v(it,0,:,:)        ; Meridional velocity       (m/s)
        slp_plane  = mslp(it-1,0,:,:)   ; MSLP                      (hPa)
      end if

    ; MSLP
      slp_min       = min(slp_plane)
      slp_smth      = smth9_Wrap(slp_plane, 0.5, 0.5, True)
      slp_max_smth  = max(slp_smth)

      dims          = dimsizes(slp_plane)
      slp1d         = ndtooned(slp_plane)
      inds          = ind_resolve(minind(slp1d),dims)
      slp1d_sm      = ndtooned(slp_smth)
      inds_sm       = ind_resolve(maxind(slp1d_sm),dims)

      lat_minslp = lat(0) + (dy * inds(0,0))  ; Distance between lat. points - 0.04
      lon_minslp = lon(0) + (dx * inds(0,1))  ; Distance between lon. points - 0.04
      print("Cyclone centre (mslp): "+lat_minslp+" degrees N, "+lon_minslp+" degrees E")

   ; EMPLOY A SAFETY NET IF WE IDENTIFY THE WRONG CYCLONE
     if (ct.eq.0) then
      centre(ct,0) = lat_minslp
      centre(ct,1) = lon_minslp

     else
      lt0  = centre(ct-1,0)             ; Previous TC latitude
      ln0  = centre(ct-1,1)             ; Previous TC longitude
      print("Previous TC latitude = "+lt0+" ; Previous TC longitude = "+ln0)
      d_lat = abs(lat_minslp-lt0)       ; Change in latitude (t1 - t0)
      d_lon = abs(lon_minslp-ln0)       ; Change in longitude
      print("Change in latitude = "+d_lat+" ; Change in longitude = "+d_lon)

      if (d_lat.gt.dist.or.d_lon.gt.dist) then     ; Recalculate TC centre if incorrect

     ; Create smaller grid [0.5 degrees]
       lt1 = lt0 - 0.5
       lt2 = lt0 + 0.5
       ln1 = ln0 - 0.5
       ln2 = ln0 + 0.5

       slp_plane0 = mslp(it-1,0,{lt1:lt2},{ln1:ln2})	 ; Smaller grid (1.0 x 1.0 deg)
       slp_min0   = min(slp_plane0)                      ; Find minimum on smaller grid
       dims_p0    = dimsizes(slp_plane0)                 ; Size of smaller grid
       slp_1d0    = ndtooned(slp_plane0)                 ; Create 1-D array
       inds_p0    = ind_resolve(minind(slp_1d0),dims_p0) ; Find index of minimum

       lat_minslp = lt1 + (dy * inds_p0(0,0))
       lon_minslp = ln1 + (dy * inds_p0(0,1))
       print("New centre: "+lat_minslp+" degrees N, "+lon_minslp+" degrees E")
       centre(ct,0) = lat_minslp
       centre(ct,1) = lon_minslp

       delete([/slp_plane0, slp_min0, dims_p0, slp_1d0, inds_p0/])

      else
       print("Cyclone centre: "+lat_minslp+" degrees N, "+lon_minslp+" degrees E")
       centre(ct,0) = lat_minslp
       centre(ct,1) = lon_minslp
      end if

     end if

    ; Tidy up
      delete([/slp1d,inds,slp1d_sm,slp_smth,inds_sm/])

    ; Calculate 10-m wind speed
      speed = sqrt(u_plane^2.0 + v_plane^2.0)
      copy_VarCoords(u_plane, speed)
      copy_VarAtts(u_plane, speed)
      copy_VarMeta(u_plane, speed)
      speed@units = "m s~S~-1~N~"

;================================
; Add min/max values to arrays 
;================================

    rad	= 1

    slp_box(ct) = min(slp_plane({lat_minslp-rad:lat_minslp+rad},\
				{lon_minslp-rad:lon_minslp+rad}))
    spd_box(ct) = max(speed({lat_minslp-rad:lat_minslp+rad},\
    		            {lon_minslp-rad:lon_minslp+rad}))

    ct = ct + 1 ; Counter variable (time)

    end do     ; End time loop

;====================================================================
; Calculate 24 h difference in windspeed - how many periods of RI? 
;====================================================================

  ; Find all periods of RI 
    do it = 0, 96, 6

     tf = it			     ; First time 
     tn = it+24	 		     ; Next time 
     ri = spd_box(tn) - spd_box(tf)  ; Compute 24 h change in windspeed 
     rw = -ri

     if (ri.ge.15) then
      print("RI occurs in simulation "+ens_str+" between t = "+tf+" and t = "+tn)
      print("The windspeed increased by "+ri+" m/s in this 24-h period.")
     end if 

    end do 

  ; Find all periods of weakening (ignore landfall)
    do it = 0, 72, 6

     tf = it                         ; First time
     tn = it+24                      ; Next time
     ri = spd_box(tn) - spd_box(tf)  ; Compute 24 h change in windspeed
     rw	= -ri

     if (ri.le.0) then
      print("Weakening occurs in simulation "+ens_str+" between t = "+tf+" and t = "+tn)
      print("The windspeed decreased by "+rw+" m/s in this 24-h period.")
     end if

    end do

;=================================
; Write out values to text file 
;=================================

    diri = "$sam/nepartak/images/cart_coords/cp_um/ens/cp_"+intA+"_"

  ; 'ens_str' is the number of the ensemble member (000 -- 011)

  ; MSLP
    mslp_out = diri+ens_str+"_slp.txt"
    asciiwrite(mslp_out,slp_box)

  ; Windspeed
    spd_out  = diri+ens_str+"_spd.txt"
    asciiwrite(spd_out,spd_box)

  end do	; End input file loop (do nf = 0, numINPUT-1)

end 