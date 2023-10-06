; Script to calculate circulation budget terms every hour for CP ensemble simulations
; To run script, use command line options ['ens', 'opt', rad, run, nave]
; e.g. [ncl p1=1 p2=2 int=1 'opt="pdf"' rad=3 run=1 tc_ens_xz_cbgt.ncl]

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$LIB/gsn_csm.ncl"

begin

;=========================================================
; Calculate storm motion [for system-relative velocity]
;=========================================================

  ; Prelash (define arrays for loops below)
    ens_arr  = (/ens0,ens1/)         ; Ensemble member string array
    centre   = new((/2,78,2/),float) ; Cyclone centre at each time interval
    x_cyc    = new((/2,78/),float)   ; Cyclone longitude
    y_cyc    = new((/2,78/),float)   ; Cyclone latitude
    u_cyc    = new((/2,78/),float)   ; Zonal cyclone velocity
    v_cyc    = new((/2,78/),float)   ; Meridional cyclone velocity
    vel_cyc  = new((/2,78/),float)   ; Cyclone velocity (vector)
    pt	     = 0		     ; Counter variable (ensemble member)

 do en = p1, p2, int

    print_clock("Working on ensemble member: "+ens_arr(en))

  ; List all files to read in and analyse
    diri        = "$sam/um/cp/ens/20160704T0000Z_ra1t_"
    fili_pb     = diri+ens_arr(en)+"_pb_"
    fili_pc     = diri+ens_arr(en)+"_pc_"

  ; 00Z 4 Jul [00], 06Z 4 Jul [01], 12Z 4 Jul [02], 18Z 4 Jul [03], 00Z 5 Jul [04],
  ; 06Z 5 Jul [05], 12Z 5 Jul [06], 18Z 5 Jul [07], 00Z 6 Jul [08], 06Z 6 Jul [09],
  ; 12Z 6 Jul [10], 18Z 6 Jul [11], 00Z 7 Jul [12], 06Z 7 Jul [13], 12Z 7 Jul [14],
  ; 18Z 7 Jul [15], 00Z 8 Jul [16], 06Z 8 Jul [17], 12Z 8 Jul [18], 18Z 8 Jul [19],
  ; 00Z 9 Jul [20]

    input_pb    = (/fili_pb+"04.nc",fili_pb+"05.nc",fili_pb+"06.nc",fili_pb+"07.nc",\
    		    fili_pb+"08.nc",fili_pb+"09.nc",fili_pb+"10.nc",\
    		    fili_pb+"11.nc",fili_pb+"12.nc",fili_pb+"13.nc",\
		    fili_pb+"14.nc",fili_pb+"15.nc",fili_pb+"16.nc"/)
    input_pc    = (/fili_pc+"04.nc",fili_pc+"05.nc",fili_pc+"06.nc",fili_pc+"07.nc",\
    		    fili_pc+"08.nc",fili_pc+"09.nc",fili_pc+"10.nc",\
                    fili_pc+"11.nc",fili_pc+"12.nc",fili_pc+"13.nc",\
                    fili_pc+"14.nc",fili_pc+"15.nc",fili_pc+"16.nc"/)

    numINPUT    = dimsizes(input_pb)        ; Number of input files
    numTIMES    = numINPUT * 6              ; Total number of times in all files
    ct          = 0                         ; Counter variable (time)

    setvalues NhlGetWorkspaceObjectId
      "wsMaximumSize" : 1000000000
    end setvalues

;==========================
; Loop over input files
;==========================

  do nf = 0, numINPUT-1     ; Loop over all times in file [00Z - 12Z 5 July]
  f = addfile(input_pb(nf),"r")

;==================================
; Get the variables we will need
;==================================

  time  = f->t            ; Times in file (6 per file - hourly)

  lon1  = f->longitude_1  ; longitude (1098 points --> 109.04 to 152.92 degrees E)
  lat1  = f->latitude_1   ; latitude  (810 points --> 1.8 to 34.16 degrees N)
  lon   = f->longitude    ; longitude (1098 points --> 109.02 to 152.90 degrees E)
  lat   = f->latitude     ; latitude (811 points --> 1.78 to 34.18 degrees N)

  pres  = f->p            ; Pressure (levels)
  plevs = dimsizes(pres)  ; Size of pressure level array

  u     = f->u            ; Zonal velocity (m/s) [time | pres | lat | lon]
  v     = f->v            ; Meridional velocity (m/s) [time | pres | lat | lon]
  vort  = f->rvor         ; Relative vorticity (/s) [time | pres | lat | lon]
  z     = f->ht           ; Geopotential height (m) [time | pres | lat1 | lon1]

  z     = z/10
  z@units = "dam"	  ; Convert to decametres
  vort  = vort * (10 ^ -6) ; Convert to /s

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

     print_clock("Working on time: "+time_str(it))

;================================================================================
; Find cyclone centre on subset of global grid (use both height and vorticity)
;================================================================================

   ; Create array subsets to speed up calculations [400 * 489]
     lat_0 = 6.95   ; 9.05   ; 6.92
     lat_1 = 30.02  ; 27.01  ; 25.01  ; 23.52  ; 26.02
     lon_0 = 114.00 ; 120.48
     lon_1 = 142.02 ; 147.02

   ; Choose a level in the lower troposphere (950 hPa)
     ilev = 1

   ; First, concentrate on a subset of the global grid
     u_plane    = u(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
     v_plane    = v(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
     geo_plane  = z(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
     vort_plane = vort(it,ilev,{lat_0:lat_1},{lon_0:lon_1})
     dsize	= dimsizes(u_plane)

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
     print("Cyclone centre (hgt): "+lat_min950+" degrees N, "+lon_min950+" degrees E")

   ;========================================================
   ; EMPLOY A SAFETY NET IF WE IDENTIFY THE WRONG CYCLONE
   ;========================================================   
     if (ct.eq.0) then
      centre(pt,ct,0) = lat_min950
      centre(pt,ct,1) = lon_min950

     else
      lt0   = centre(pt,ct-1,0)		; Previous TC latitude
      ln0   = centre(pt,ct-1,1)		; Previous TC longitude
      print("lt0 = "+lt0+" ; ln0 = "+ln0)
      d_lat = abs(lat_min950-lt0)	; Change in latitude (t1 - t0)
      d_lon = abs(lon_min950-ln0)	; Change in longitude (note extra minus sign)
      print("d_lat = "+d_lat+" ; d_lon = "+d_lon)

      if (d_lat.gt.1.or.d_lon.gt.1) then        ; Recalculate TC centre if incorrect

     ; Create smaller grid [0.5 degrees]
       inc = 0.5
       lt1 = lt0 - inc
       lt2 = lt0 + inc
       ln1 = ln0 - inc
       ln2 = ln0 + inc

       geo_plane0 = z(it,ilev,{lt1:lt2},{ln1:ln2})       ; Smaller grid (1.0 x 1.0 deg)
       geo_min0	  = min(geo_plane0)		         ; Find minimum on smaller grid
       dims_h0	  = dimsizes(geo_plane0)	         ; Size of smaller grid
       geo_1d0	  = ndtooned(geo_plane0)	         ; Create 1-D array
       inds_h0    = ind_resolve(minind(geo_1d0),dims_h0) ; Find index of minimum
       
       lat_min950 = lt1 + (dy * inds_h0(0,0))
       lon_min950 = ln1 + (dy * inds_h0(0,1))
       print("New centre (hgt): "+lat_min950+" degrees N, "+lon_min950+" degrees E") 
       centre(pt,ct,0) = lat_min950
       centre(pt,ct,1) = lon_min950

       delete([/geo_plane0, geo_min0, dims_h0, geo_1d0, inds_h0/])

      else
       print("Cyclone centre (hgt): "+lat_min950+" degrees N, "+lon_min950+" degrees E")
       centre(pt,ct,0) = lat_min950
       centre(pt,ct,1) = lon_min950
      end if

     end if

    ; Tidy up before looping over pressure levels
      delete([/vort1d, inds, vort1d_sm, vort_smth, inds_sm/])
      delete([/u_plane, v_plane, vort_plane, geo_plane/])

      ct = ct + 1

    end do    ; End of time loop (do it = 0, times-1)

  ; Tidy up to avoid array dimension errors [but leave alone on final iteration]
     if (pt.eq.0) then
      delete([/time,times,z,u,v,vort/])
      delete([/utc_date,year,month,day,hour,\
               minute,second,date_str,time_str,out_str/])
     else
      if (nf.lt.numINPUT-1)
        delete([/time,times,z,u,v,vort/])
        delete([/utc_date,year,month,day,hour,\
                minute,second,date_str,time_str,out_str/])
      end if
     end if

   end do      ; End of input file loop (do nf = 0, numINPUT-1) 

   pt = pt + 1 ; Counter variable (ensemble member)

  end do       ; End of ensemble loop (do en = p1, p2)

;============================================================================
; Calculate distance between two points --> then calculate system velocity
;============================================================================

   print_clock("Starting system velocity calculations...")
   pt = 0   ; Counter variable (ensemble member)

 do en = p1, p2, int

 ; Loop to calculate system velocity from hourly cyclone position [lat,lon]
   do iang = 1, numTIMES-1
    x_cyc(pt,iang)   = gc_latlon(centre(pt,iang-1,0),centre(pt,iang-1,1),\
    		       centre(pt,iang-1,0),centre(pt,iang,1),2,3)
    y_cyc(pt,iang)   = gc_latlon(centre(pt,iang-1,0),centre(pt,iang-1,1),\
    		       centre(pt,iang,0),centre(pt,iang-1,1),2,3)
    u_cyc(pt,iang)   = x_cyc(pt,iang)/(3600)
    v_cyc(pt,iang)   = y_cyc(pt,iang)/(3600)
    vel_cyc(pt,iang) = sqrt(u_cyc(pt,iang)^2.0 + v_cyc(pt,iang)^2.0)
   end do

 ; Set initial system velocity [0] equal to that at the second time [1]
   u_cyc(pt,0)   = u_cyc(pt,1)
   v_cyc(pt,0)   = v_cyc(pt,1)
   vel_cyc(pt,0) = vel_cyc(pt,1)

 ; Add metadata to new arrays
   copy_VarAtts(u,u_cyc(pt,:))
   u_cyc@units         = "m s~S~-1~N~"
   u_cyc@description   = "Cyclone zonal velocity"

   copy_VarAtts(v,v_cyc(pt,:))
   v_cyc@units         = "m s~S~-1~N~"
   v_cyc@description   = "Cyclone meridional velocity"

   vel_cyc@units       = "m s~S~-1~N~"
   vel_cyc@description = "Cyclone speed"

   u_cyc(pt,:) = -(u_cyc(pt,:)) ; Easterly motion is negative

   pt = pt + 1 	 		; Counter variable (ensemble member)

  end do      ; End of ensemble loop (do en = p1, p2)

 ; Tidy up
   delete([/lat,lon,lat1,lon1,pres,plevs,u,v,z,vort/])
   delete([/time,times,utc_date,year,month,day,hour,\
            minute,second,date_str,time_str,out_str/])

;==================================================================================
;;;;;;;;;;;;;;;;;;;;;;;;;;; Start of main part of script ;;;;;;;;;;;;;;;;;;;;;;;;;;
;==================================================================================

 print_clock("Starting main part of script...")

  ; Prelash [define arrays for use in loops below]

    time_arr    = new((/2,numTIMES/),string)      ; Array to hold strings for second loop
    title_arr   = new((/2,numTIMES/),string)      ; Array to hold strings for plot titles
    pos_arr     = new((/2,numTIMES,2/),float)     ; Array for cyclone position [lat/lon]
    r_earth     = 6.37e6                          ; Earth's radius (constant)
    pt		= 0				  ; Counter variable (ensemble member)

  ; Create array to hold all values of relative vorticity, to compute tendencies
    circ_size   = (/2,19,dsize(0),dsize(1)/) ; [2 * 19 * 400 * 489]
    circ_arr    = new((/numTIMES,circ_size(0),circ_size(1),\
    		        circ_size(2),circ_size(3)/),float)
    circ_arr!0  = "time"         ; Time [dimension 0]
    circ_arr!1	= "ens"       	 ; Ensemble member [dimension 1]
    circ_arr!2	= "pres"	 ; Pressure levels [dimension 2]

  ; Create array to hold values from frictional circulation tendency
    dtaux_dz  = new((/78,19,dsize(0),dsize(1)/),float) ; Friction (x-component)
    dtauy_dz  = new((/78,19,dsize(0),dsize(1)/),float) ; Friction (y-component)

    temp_arr	= new((/78,19,dsize(0),dsize(1)/),float)
    temp_arr!0  = "time"

  ; Create arrays to hold circulation tendency values along each box edge
  ; [ensemble member (2) * box ensemble (49) * time (numTIMES) * levs (19)]
    circ_north = new((/2,49,numTIMES,19/),float)
    circ_south = new((/2,49,numTIMES,19/),float)
    circ_west  = new((/2,49,numTIMES,19/),float)
    circ_east  = new((/2,49,numTIMES,19/),float)
    flux_north = new((/2,49,numTIMES,19/),float)
    flux_south = new((/2,49,numTIMES,19/),float)
    flux_west  = new((/2,49,numTIMES,19/),float)
    flux_east  = new((/2,49,numTIMES,19/),float)
    tilt_north = new((/2,49,numTIMES,19/),float)
    tilt_south = new((/2,49,numTIMES,19/),float)
    tilt_west  = new((/2,49,numTIMES,19/),float)
    tilt_east  = new((/2,49,numTIMES,19/),float)
    fric_north = new((/2,49,numTIMES,19/),float)
    fric_south = new((/2,49,numTIMES,19/),float)
    fric_west  = new((/2,49,numTIMES,19/),float)
    fric_east  = new((/2,49,numTIMES,19/),float)

;==================================================
; Start loop over ensemble members [em01, em05]
;==================================================

 do en = p1, p2, int

    print_clock("Working on ensemble member: "+ens_arr(en))

  ; List all files to read in and analyse
    diri        = "$sam/um/cp/ens/20160704T0000Z_ra1t_"
    fili_pb     = diri+ens_arr(en)+"_pb_"
    fili_pc     = diri+ens_arr(en)+"_pc_"
    fili_pf	= diri+ens_arr(en)+"_pf_"

  ; 00Z 4 Jul [00], 06Z 4 Jul [01], 12Z 4 Jul [02], 18Z 4 Jul [03], 00Z 5 Jul [04],
  ; 06Z 5 Jul [05], 12Z 5 Jul [06], 18Z 5 Jul [07], 00Z 6 Jul [08], 06Z 6 Jul [09],
  ; 12Z 6 Jul [10], 18Z 6 Jul [11], 00Z 7 Jul [12], 06Z 7 Jul [13], 12Z 7 Jul [14],
  ; 18Z 7 Jul [15], 00Z 8 Jul [16], 06Z 8 Jul [17], 12Z 8 Jul [18], 18Z 8 Jul [19],
  ; 00Z 9 Jul [20]

    input_pb    = (/fili_pb+"04.nc",fili_pb+"05.nc",fili_pb+"06.nc",fili_pb+"07.nc",\
                    fili_pb+"08.nc",fili_pb+"09.nc",fili_pb+"10.nc",\
                    fili_pb+"11.nc",fili_pb+"12.nc",fili_pb+"13.nc",\
                    fili_pb+"14.nc",fili_pb+"15.nc",fili_pb+"16.nc"/)
    input_pc    = (/fili_pc+"04.nc",fili_pc+"05.nc",fili_pc+"06.nc",fili_pc+"07.nc",\
                    fili_pc+"08.nc",fili_pc+"09.nc",fili_pc+"10.nc",\
                    fili_pc+"11.nc",fili_pc+"12.nc",fili_pc+"13.nc",\
                    fili_pc+"14.nc",fili_pc+"15.nc",fili_pc+"16.nc"/)
    input_pf    = (/fili_pf+"04.nc",fili_pf+"05.nc",fili_pf+"06.nc",fili_pf+"07.nc",\
                    fili_pf+"08.nc",fili_pf+"09.nc",fili_pf+"10.nc",\
                    fili_pf+"11.nc",fili_pf+"12.nc",fili_pf+"13.nc",\
                    fili_pf+"14.nc",fili_pf+"15.nc",fili_pf+"16.nc"/)

    numINPUT    = dimsizes(input_pb)        ; Number of input files
    numTIMES    = numINPUT * 6              ; Total number of times in all files
    ct          = 0                         ; Counter variable

    setvalues NhlGetWorkspaceObjectId
      "wsMaximumSize" : 1000000000
    end setvalues

;===============================
; Start loop over input files 
;===============================

  do nf = 0, numINPUT-1
  f = addfile(input_pb(nf),"r")		; Read in 'pb' stream [u,v,w,vort,div]
  g = addfile(input_pc(nf),"r")		; Read in 'pc' stream [friction]
  h = addfile(input_pf(nf),"r")		; Read in 'pf' stream [density on model levs]

;==============================================
; Get the variables we will need [pb stream]
;==============================================

  time  = f->t            ; Times in file (6 per file - hourly)
  times = dimsizes(time)

  lon1  = f->longitude_1  ; longitude (1098 points --> 109.04 to 152.92 degrees E)
  lat1  = f->latitude_1   ; latitude  (810 points --> 1.8 to 34.16 degrees N)
  lon   = f->longitude    ; longitude (1098 points --> 109.02 to 152.90 degrees E)
  lat   = f->latitude     ; latitude (811 points --> 1.78 to 34.18 degrees N)

  pres  = f->p		  ; Pressure (levels)
  plevs	= dimsizes(pres)  ; Size of pressure level array
  pres@units = "hPa"

  u     = f->u            ; Zonal velocity (m/s) [time | pres | lat | lon]
  v     = f->v            ; Meridional velocity (m/s) [time | pres | lat | lon]
  vort  = f->rvor         ; Relative vorticity (/s) [time | pres | lat | lon]
  div   = f->div          ; Divergence (/s) [time | pres | lat | lon]

  t     = f->temp	  ; Temperature (K) [time | pres | lat1 | lon1]
  z     = f->ht   	  ; Geopotential height (m) [time | pres | lat1 | lon1]
  w     = f->dz_dt	  ; Vertical velocity (m/s) [time | pres | lat1 | lon1]

; Convert to units we want
  z     = z/10
  z@units = "dam" ; Convert to decametres

  vort  = vort * (10 ^ -6) ; Convert to /s
  div   = div  * (10 ^ -6) ; Convert to /s

  dx    = lon(1) - lon(0) ; Grid spacing (longitude)
  dy    = lat(1) - lat(0) ; Grid spacing (latitude)

 ; Pressure levels (pres)
 ; (0)   1000	  	
 ; (1)   950		
 ; (2)   900		
 ; (3)   850
 ; (4)	 800
 ; (5)	 750		
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

 ; Calculate absolute vorticity
   d = dimsizes(vort)
   f0  = new((/d(0),d(1),d(2),d(3)/),float)  ; Earth's vorticity
   avo = new((/d(0),d(1),d(2),d(3)/),float)  ; Absolute vorticity

 ; Calculate the Coriolis parameter at all grid points
   do ilon = 0, d(3) - 1
      do ilat = 0, d(2) - 1
         f0(:,:,ilat,ilon) = 2 * (7.27 * 10 ^ -5) * sin(lat_rad(ilat))
      end do
   end do
   copy_VarMeta(vort,f0)

   avo           = vort + f0 ; Absolute vorticity
   copy_VarMeta(vort,avo)

   avo@long_name = "Vorticity (planetary)"
   avo@standard_name = "atmosphere_planetary_vorticity"
   avo@title     = "Absolute vorticity"
   avo@name      = "avo"

;==============================================
; Get the variables we will need [pc stream]
;==============================================

 ; Get additional variables from the model-level file
   lon1F = g->longitude_1  ; longitude (1098 points --> 109.04 to 152.92 degrees E)
   lat1F = g->latitude_1   ; latitude  (811 points --> 1.78 to 34.18 degrees N)
   lonF  = g->longitude    ; longitude (1098 points --> 109.02 to 152.9 degrees E)
   latF  = g->latitude     ; latitude (810 points --> 1.8 to 34.16 degrees N)
   hyb   = g->hybrid_ht    ; Array of 63 levels, from 5 m up to ~ 20 km
   strx  = g->taux         ; Wind stress, x-comp (N m/2) [t_1 | hyb_1 | latF  | lon1F]
   stry  = g->tauy         ; Wind stress, y-comp (N m/2) [t_1 | hyb_1 | lat1F | lonF ]

   dxF   = lonF(1) - lonF(0) ; Grid spacing (longitude)
   dyF   = latF(1) - latF(0) ; Grid spacing (latitude)   

;=========================================================================
; Get density from pf stream [same size as 'strx', smaller than 'stry']
;=========================================================================

   rho1  = h->field27	   ; Dry air density (kg/m-3) [t | hyb | lat | lon]

;  lon    = h->longitude    ; longitude (1098 points --> 109.04 to 152.92 degrees E)
;  lat    = f->latitude     ; latitude  (810 points --> 1.8 to 34.16 degrees N)

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

     time_arr(pt,ct) = out_str(it)

     print_clock("Working on time: "+time_str(it))
     title_arr(pt,ct) = "Valid at: "+time_str(it)

;==========================================================================
; Define new arrays for terms in vorticity equation (no friction for now)
;==========================================================================

  ; Create array subsets to speed up calculations [400 * 489]
    lat_0 = 6.95   ; 9.05   ; 6.92
    lat_1 = 30.02  ; 27.01  ; 25.01  ; 23.52  ; 26.02
    lon_0 = 114.00 ; 120.48
    lon_1 = 142.02 ; 147.02

    geo_sub    = z(it,:,{lat_0:lat_1},{lon_0:lon_1})
    tk_sub     = t(it,:,{lat_0:lat_1},{lon_0:lon_1})
    f0_sub     = f0(it,:,{lat_0:lat_1},{lon_0:lon_1})
    omega_sub  = w(it,:,{lat_0:lat_1},{lon_0:lon_1})
    vort_sub   = vort(it,:,{lat_0:lat_1},{lon_0:lon_1})
    avo_sub    = avo(it,:,{lat_0:lat_1},{lon_0:lon_1})
    div_sub    = div(it,:,{lat_0:lat_1},{lon_0:lon_1})
    u_sub      = u(it,:,{lat_0:lat_1},{lon_0:lon_1})
    v_sub      = v(it,:,{lat_0:lat_1},{lon_0:lon_1})

    temp_arr(ct,:,:,:) = tk_sub		; Add values to larger array

  ; Storm-relative velocity
    u_rel      = u_sub - u_cyc(pt,ct)
    v_rel      = v_sub - v_cyc(pt,ct)

  ; Add metadata
    avo_sub@units  = "s~S~-1~N~"
    f0_sub@units   = "s~S~-1~N~"

    d_sub     = dimsizes(vort_sub)    ; Size of arrays
    lat_sub   = lat({lat_0:lat_1})    ; Subset of lat array
    lon_sub   = lon({lon_0:lon_1})    ; Subset of lon array
    lat_subR  = lat_sub * (pii/180)   ; Convert to radians for calculation below
    lon_subR  = lon_sub * (pii/180)

   ; Troubleshooting (look for abnormally large values and exit the script if found)
     do p1a = 0, d_sub(0)-1    ; pressure
       do l1 = 0, d_sub(1)-1   ; latitude
         do l2 = 0, d_sub(2)-1 ; longitude
           if (u_rel(p1a,l1,l2) .gt. 200)
            print("Quitting because of storm-relative velocity error!!!")
             exit()
           end if
         end do
       end do
     end do

    tilt_x   = new(d_sub,float) ; Tilting (x-component)
    tilt_y   = new(d_sub,float) ; Tilting (y-component)

    adv_fl_x = new(d_sub,float)	; Advective flux (x-component)
    adv_fl_y = new(d_sub,float)	; Advective flux (y-component)

  ; Add basic metadata (modify details later)
    copy_VarMeta(vort_sub,tilt_x)
    copy_VarAtts(vort_sub,tilt_x)
    copy_VarCoords(vort_sub,tilt_x)

    copy_VarMeta(vort_sub,tilt_y)
    copy_VarAtts(vort_sub,tilt_y)
    copy_VarCoords(vort_sub,tilt_y)

    copy_VarMeta(vort_sub,adv_fl_x)
    copy_VarAtts(vort_sub,adv_fl_x)
    copy_VarCoords(vort_sub,adv_fl_x)

    copy_VarMeta(vort_sub,adv_fl_y)
    copy_VarAtts(vort_sub,adv_fl_y)
    copy_VarCoords(vort_sub,adv_fl_y)

    copy_VarMeta(u_sub,u_rel)
    copy_VarAtts(u_sub,u_rel)
    copy_VarCoords(u_sub,u_rel)

    copy_VarMeta(v_sub,v_rel)
    copy_VarAtts(v_sub,v_rel)
    copy_VarCoords(v_sub,v_rel)

;==========================================================
; Define latitude and longitude for finite differencing
;==========================================================

    y2d     = new((/d_sub(1),d_sub(2)/),float)
    y2d!0   = "lat"
    y2d!1   = "lon"
    y2d&lat =  lat_sub
    y2d&lon =  lon_sub

    y2d     = conform_dims(dimsizes(y2d),lat_sub,0) ; Array of latitude values
    coslat  = cos(y2d*pii/180)                      ; Same array in radians

  ; Determine horizontal resolution of data
    delta_lon = y2d
    delta_lat = y2d
    delta_lon = lon_sub(1) - lon_sub(0) ; Zonal grid spacing
    delta_lat = lat_sub(1) - lat_sub(0) ; Meridional grid spacing
    delta_lon = sqrt(delta_lon^2)
    delta_lat = sqrt(delta_lat^2)

    lons = lon_sub(1) - lon_sub(0)
    lats = lat_sub(1) - lat_sub(0)
    lons = sqrt(lons^2)
    lats = sqrt(lats^2)

  ; Calculate dx and dy (changes with longitude but not with latitude)
    ddx = y2d
    ddx = r_earth*(pii/180.0)*delta_lon*cos(y2d*pii/180)
    ddy = y2d
    ddy = r_earth*(pii/180)*delta_lat

  ; 3-D latitude array
    coslat3d = conform_dims(dimsizes(u_sub),coslat,(/1,2/))

;=================================
; Friction-related calculations 
;=================================

  ; Array of hybrid levels to approximately match pressure levels every 50 hPa
  ; Array [/1000, 950, 900, 850, 800, 750, etc, ..., 100/]
    lev_arr = (/1, 8, 13, 16, 19, 21, 23, 25, 27, 29, \
                31, 33, 35, 37, 39, 42, 47, 51, 56/)
    lev_arr@units = "m"

    l0    = 0  ; Subset of vertical levels (19 levels)
    l1    = 18 ; Level 56 of original array ~ 100 hPa (Int. Standard Atmosphere)

  ; Arrays to hold subset of 'strx' and 'stry', on levels in 'lev_arr'
    fric_x_hgt  = new((/dimsizes(lev_arr),dsize(0),dsize(1)/),float)
    fric_y_hgt  = new((/dimsizes(lev_arr),dsize(0),dsize(1)/),float)

  ; Array of pressure levels with the same dimensions as 'lev_arr'
    pres       = (/1000,950,900,850,800,750,700,650,600,550,\
                   500,450,400,350,300,250,200,150,100/)
    pres@units = "hPa"
    plevs      = dimsizes(pres)

    cp         = 0    ; Counter variable

  ; Create array of pressure values with same dimensions as 'fric_x_hgt'
    pres_fr    = conform_dims(dimsizes(fric_x_hgt), pres, 0)
    pres_fr = pres_fr * 100
    pres_fr@units = "Pa"

;===============================================================
; Tilting term (w * dv_dp, -w * du_dp)
;===============================================================

  ; Calculate vertical derivatives of horizontal [storm-relative] velocity
    du_dp = u_rel
    dv_dp = v_rel

  ; Convert pressure values to Pa for calculation below
    pres = pres * 100
    pres@units = "Pa"

  ; Fill arrays with vertical derivatives at each pressure level
  ; (0) 1000, (1) 950, ..., (17) 150, (18) 100
    do plevsN = 0, plevs-1

       kp1 = min((/plevsN+1,plevs-1/)) ; Account for top and bottom levels
       km1 = max((/plevsN-1,0/))

       du_dp(plevsN,:,:) = (u_rel(kp1,:,:) - \
                            u_rel(km1,:,:))\
                            /((pres(km1))-(pres(kp1)))
       dv_dp(plevsN,:,:) = (v_rel(kp1,:,:) - \
                            v_rel(km1,:,:))\
                            /((pres(km1))-(pres(kp1)))

    end do

  ; Convert pressure values back to hPa
    pres = pres / 100
    pres@units = "hPa"

  ; Calculate [dv_dp * w]
    tilt_y = dv_dp * omega_sub    ; X-component of relative vorticity

  ; Calculate [-du_dp * w]
    tilt_x = -(du_dp * omega_sub) ; Y-component of relative vorticity

  ; Metadata
    tilt_x@units = "m s~S~-2~N~"
    tilt_x@standard_name = "tilting_x"
    tilt_x@long_name = "Tilting (x-component)"
    tilt_x@title = "Tilting (x-component)"
    tilt_x@name = "tilting_x"

    tilt_y@units = "m s~S~-2~N~"
    tilt_y@standard_name = "tilting_y"
    tilt_y@long_name = "Tilting (y-component)"
    tilt_y@title = "Tilting (y-component)"
    tilt_y@name = "tilting_y"

;=======================================================
; Advective flux (horizontal advection + stretching) 
;=======================================================

  ; x-component
    adv_fl_x = u_rel * avo_sub 

  ; y-component
    adv_fl_y = v_rel * avo_sub

  ; Metadata
    adv_fl_x@standard_name = "adv_flux_x"
    adv_fl_x@long_name = "Advective flux (x-component)"
    adv_fl_x@title = "Advective flux (x-component)"
    adv_fl_x@name = "ad_fl_x"

    adv_fl_y@standard_name = "adv_flux_y"
    adv_fl_y@long_name = "Advective flux (y-component)"
    adv_fl_y@title = "Advective flux (y-component)"
    adv_fl_y@name = "ad_fl_y"

;================================================================================
; Find cyclone centre on subset of global grid (use both height and vorticity)
;================================================================================

   ; Choose a level in the lower troposphere (950 hPa)
     ilev = 1

   ; First, concentrate on a subset of the global grid 
     geo_plane     = geo_sub(ilev,{lat_0:lat_1},{lon_0:lon_1})
     vort_plane    = vort_sub(ilev,{lat_0:lat_1},{lon_0:lon_1})

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
     geo_min	    = min(geo_plane)
     geo_smth	    = smth9_Wrap(geo_plane, 0.5, 0.5, True)
     geo_min_smth   = min(geo_smth)
      
     dims_h	    = dimsizes(geo_plane)
     geo1d	    = ndtooned(geo_plane)
     inds_h	    = ind_resolve(minind(geo1d),dims_h)
     geo1d_sm	    = ndtooned(geo_smth)
     inds_h_sm	    = ind_resolve(minind(geo1d_sm),dims_h)            

     lat_min950 = lat_0 + (dy * inds_h(0,0))
     lon_min950 = lon_0 + (dx * inds_h(0,1))
     print("Cyclone centre (hgt): "+lat_min950+" degrees N, "+lon_min950+" degrees E")

   ;========================================================
   ; EMPLOY A SAFETY NET IF WE IDENTIFY THE WRONG CYCLONE
   ;========================================================
     if (ct.eq.0) then
      centre(pt,ct,0) = lat_min950
      centre(pt,ct,1) = lon_min950

     else
      lt0   = centre(pt,ct-1,0)         ; Previous TC latitude
      ln0   = centre(pt,ct-1,1)         ; Previous TC longitude
      print("lt0 = "+lt0+" ; ln0 = "+ln0)
      d_lat = abs(lat_min950-lt0)       ; Change in latitude (t1 - t0)
      d_lon = abs(lon_min950-ln0)       ; Change in longitude (note extra minus sign)
      print("d_lat = "+d_lat+" ; d_lon = "+d_lon)

      if (d_lat.gt.1.or.d_lon.gt.1) then        ; Recalculate TC centre if incorrect

     ; Create smaller grid [0.5 degrees]
       lt1 = lt0 - 0.5
       lt2 = lt0 + 0.5
       ln1 = ln0 - 0.5
       ln2 = ln0 + 0.5

       geo_plane0 = z(it,ilev,{lt1:lt2},{ln1:ln2})       ; Smaller grid (1.0 x 1.0 deg)
       geo_min0   = min(geo_plane0)                      ; Find minimum on smaller grid
       dims_h0    = dimsizes(geo_plane0)                 ; Size of smaller grid
       geo_1d0    = ndtooned(geo_plane0)                 ; Create 1-D array
       inds_h0    = ind_resolve(minind(geo_1d0),dims_h0) ; Find index of minimum

       lat_min950 = lt1 + (dy * inds_h0(0,0))
       lon_min950 = ln1 + (dy * inds_h0(0,1))
       print("New centre (hgt): "+lat_min950+" degrees N, "+lon_min950+" degrees E")
       centre(pt,ct,0) = lat_min950
       centre(pt,ct,1) = lon_min950

       delete([/geo_plane0, geo_min0, dims_h0, geo_1d0, inds_h0/])

      else
       print("Cyclone centre (hgt): "+lat_min950+" degrees N, "+lon_min950+" degrees E")
       centre(pt,ct,0) = lat_min950
       centre(pt,ct,1) = lon_min950
      end if

     end if

   ; Get information about box size for later ensemble calculations
     r = rad   ; Box radius (degrees)
     box = ddx({lat_min950-r:lat_min950+r},{lon_min950-r:lon_min950+r})
     box_size = dimsizes(box)
     size0    = box_size(0)
     size1    = box_size(1)

   ; Tidy up before looping over pressure levels
     delete([/vort1d, inds, vort1d_sm, vort_smth, inds_sm/])
     delete([/vort_plane, geo_plane/])

;============================================================
; Loop over pressure levels to calculate vertical profile
;============================================================

    cp = 0               ; Counter variable

    do ilev = 0, plevs-1 ; 1000,950,900,850,800,750,700,650,600,550,...,100

   ; Pressure levels ['pb' stream]
     p = pres(ilev)
     print_clock("Working on: " + p + " hPa")

   ;========================================
   ; Model levels ['pc' and 'pf' streams]
   ;========================================

     hy = sprintf("%0.0f",hyb(ilev))
     print_clock("Working on hybrid level " + ilev+ " ("+hy+" m AGL)")
     m0 = lev_arr(ilev)	      ; Selected model levels, chosen to match pressure levels

   ; Calculate components of friction term
     ip1 = min((/m0+1,lev_arr(l1)/)) ; Account for top and bottom levels
     im1 = max((/m0-1,lev_arr(l0)/))

   ; Current level [wind stress]
     taux   = strx(it,m0,{lat_0:lat_1},{lon_0:lon_1})
     tauy   = stry(it,m0,{lat_0:lat_1},{lon_0:lon_1})

   ; Density ['pf' stream]
     rho    = rho1(it,m0,{lat_0:lat_1},{lon_0:lon_1})
     rho_a = (1/rho) ; Specific volume (inverse density)

   ; Levels above/below for vertical derivatives [wind stress]
     tau_x1 = strx(it,ip1,{lat_0:lat_1},{lon_0:lon_1})
     tau_x0 = strx(it,im1,{lat_0:lat_1},{lon_0:lon_1})
     tau_y1 = stry(it,ip1,{lat_0:lat_1},{lon_0:lon_1})
     tau_y0 = stry(it,im1,{lat_0:lat_1},{lon_0:lon_1})

   ; Vertical derivative of 'strx'
     dtaux_dz(ct,ilev,:,:) = (tau_x1 - tau_x0) / (hyb(ip1) - hyb(im1))

   ; Vertical derivative of 'stry'
     dtauy_dz(ct,ilev,:,:) = (tau_y1 - tau_y0) / (hyb(ip1) - hyb(im1))

   ; X-component of friction
     fric_x_plane           = rho_a(:,:) * dtaux_dz(ct,ilev,:,:)
     copy_VarCoords(taux,fric_x_plane)
     copy_VarMeta(taux,fric_x_plane)
     fric_x_plane@units     = "m s~S~-2~N~"
     fric_x_plane@title     = "x-component of friction"
     fric_x_plane@long_name = "x-component of friction"
     fric_x_plane@name      = "fric_x"     

   ; Y-component of friction
     fric_y_plane           = rho_a(:,:) * dtauy_dz(ct,ilev,:,:)
     copy_VarCoords(tauy,fric_y_plane)
     copy_VarMeta(tauy,fric_y_plane)
     fric_y_plane@units     = "m s~S~-2~N~"
     fric_y_plane@title     = "y-component of friction"
     fric_y_plane@long_name = "y-component of friction"
     fric_y_plane@name      = "fric_y"

   ;=================================
   ; Pressure levels ['pb' stream]
   ;=================================

   ; Storm-relative velocity
     u_plane       = u_rel(ilev,{lat_0:lat_1},{lon_0:lon_1})
     v_plane       = v_rel(ilev,{lat_0:lat_1},{lon_0:lon_1})

   ; Other variables
     vort_plane    = vort_sub(ilev,{lat_0:lat_1},{lon_0:lon_1})
     avo_plane	   = avo_sub(ilev,{lat_0:lat_1},{lon_0:lon_1})
     div_plane	   = div_sub(ilev,{lat_0:lat_1},{lon_0:lon_1})
     geo_plane     = geo_sub(ilev,{lat_0:lat_1},{lon_0:lon_1})
     t_plane       = tk_sub(ilev,{lat_0:lat_1},{lon_0:lon_1})

     vort_max      = max(vort_plane)
     vort_smth     = smth9_Wrap(vort_plane, 0.5, 0.5, True)
     vort_max_smth = max(vort_smth)

     zmin          = toint(min(geo_plane)) ; Minimum geopotential on pressure level
     zmax          = toint(max(geo_plane)) ; Maximum geopotential on pressure level

     dims          = dimsizes(vort_plane)
     vort1d        = ndtooned(vort_plane)
     inds          = ind_resolve(maxind(vort1d),dims)
     vort1d_sm     = ndtooned(vort_smth)
     inds_sm       = ind_resolve(maxind(vort1d_sm),dims)

     lat_max = lat_0 + (dy * inds(0,0))
     lon_max = lon_0 + (dx * inds(0,1))

     lat_plane = ind(lat_0.le.lat.and.lat.le.lat_1) ; Indices of latitude array
     lon_plane = ind(lon_0.le.lon.and.lon.le.lon_1) ; Indices of longitude array

   ; Modify latitude/longitude arrays
     lat_size = dimsizes(lat_plane)		  	  	; Array sizes
     lon_size = dimsizes(lon_plane)

     lat_sub = lat(lat_plane(0):lat_plane(lat_size-1)) ; Create subset of points 
     lon_sub = lon(lon_plane(0):lon_plane(lon_size-1)) ; between values given above
     delete([/vort1d, inds, vort1d_sm, vort_smth, inds_sm/])

   ; Account for cases in which another maximum is chosen in error 
     lat_diff = abs(lat_max-lat_min950)
     lon_diff = abs(lon_max-lon_min950)

     if (lon_diff.gt.1.0.or.lat_diff.gt.1.0) then
        lat_max = lat_min950
        lon_max = lon_min950
     end if
     print("Cyclone centre: "+lat_max+" degrees N, "+lon_max+" degrees E")

   ; Add lat/lon values to array for later calculations
     pos_arr(pt,ct,0) = lat_max
     pos_arr(pt,ct,1) = lon_max

   ; Terms in the circulation budget
   ; 'tilt_x', 'tilt_y' contain values on all pressure levels for each time
   ; 'tilt_x_plane' and 'tilt_y_plane' contain values on single pressure levels
     tilt_x_plane               = tilt_x(ilev,:,:)       ; Tilting (x-component)
     tilt_y_plane               = tilt_y(ilev,:,:)       ; Tilting (y-component)

;====================================================
; Define arrays to hold ensemble of box positions
;====================================================

     ddx_box    = new((/49,size0,size1/),float)
     ddy_box    = new((/49,size0,size1/),float)
     tilt_x_box = new((/49,size0,size1/),float)
     tilt_y_box = new((/49,size0,size1/),float)
     vort_box   = new((/49,size0,size1/),float)
     avo_box    = new((/49,size0,size1/),float)
     div_box    = new((/49,size0,size1/),float)
     u_box      = new((/49,size0,size1/),float)
     v_box      = new((/49,size0,size1/),float)
     fric_x_box = new((/49,size0,size1/),float)
     fric_y_box = new((/49,size0,size1/),float)

     lat_d      = new(49,float)
     lon_d      = new(49,float)
     lat_u      = new(49,float)
     lon_u      = new(49,float)

     d_ave      = new(49,float)
     dz0        = new(49,float)
     dz1        = new(49,float)
     div_ave    = new(49,float)
     vort_ave   = new(49,float)
     vort_avg   = new(49,float)
     avo_avg    = new(49,float)
     u_avg      = new(49,float)
     v_avg      = new(49,float)

     vort_avg1  = new(49,float)
     avo_avg1   = new(49,float)
     u_avg1     = new(49,float)
     v_avg1     = new(49,float)

;=====================================================================================
; Create array of box positions either side of the original [Davis & Galarneau 2009]
;=====================================================================================

     lat_box = (/lat_max-(3*dy),lat_max-(2*dy),lat_max-(dy),lat_max,\
     	         lat_max+(dy),lat_max+(2*dy),lat_max+(3*dy)/)
     lon_box = (/lon_max-(3*dx),lon_max-(2*dx),lon_max-(dx),lon_max,\
                 lon_max+(dx),lon_max+(2*dx),lon_max+(3*dx)/)
    
;=============================================
; Calculate terms in the circulation budget
;=============================================

   ; Radius of box for calculation (degrees)
     r = rad
     r1 = sprintf("%0.1fdeg",2*r) ; Create string to add to output file

;===================================================================
; Loop over perturbed box positions, either side of the TC centre
;===================================================================

   do i = 0, dimsizes(lon_box)-1    ; Longitude
     do j = 0, dimsizes(lat_box)-1  ; Latitude

      a	  = (dimsizes(lon_box) * i) + j

      printVarSummary(ddx_box(a,:,:))
      printVarSummary(ddx({lat_box(j)-r:lat_box(j)+r},{lon_box(i)-r:lon_box(i)+r}))

      ddx_box(a,:,:)    = ddx({lat_box(j)-r:lat_box(j)+r},{lon_box(i)-r:lon_box(i)+r})
      ddy_box(a,:,:)    = ddy({lat_box(j)-r:lat_box(j)+r},{lon_box(i)-r:lon_box(i)+r})
      tilt_x_box(a,:,:) = tilt_x_plane({lat_box(j)-r:lat_box(j)+r},\
      		      	               {lon_box(i)-r:lon_box(i)+r})
      tilt_y_box(a,:,:) = tilt_y_plane({lat_box(j)-r:lat_box(j)+r},\
			               {lon_box(i)-r:lon_box(i)+r})

      vort_box(a,:,:)   = vort_plane({lat_box(j)-r:lat_box(j)+r},\
				     {lon_box(i)-r:lon_box(i)+r})
      avo_box(a,:,:)    = avo_plane({lat_box(j)-r:lat_box(j)+r},\
      		      	            {lon_box(i)-r:lon_box(i)+r})
      div_box(a,:,:)    = div_plane({lat_box(j)-r:lat_box(j)+r},\
				    {lon_box(i)-r:lon_box(i)+r})
      u_box(a,:,:)      = u_plane({lat_box(j)-r:lat_box(j)+r},\
      		      	          {lon_box(i)-r:lon_box(i)+r})
      v_box(a,:,:)      = v_plane({lat_box(j)-r:lat_box(j)+r},\
			 	  {lon_box(i)-r:lon_box(i)+r})
      fric_x_box(a,:,:) = fric_x_plane({lat_box(j)-r:lat_box(j)+r},\
                                  {lon_box(i)-r:lon_box(i)+r})
      fric_y_box(a,:,:)	= fric_y_plane({lat_box(j)-r:lat_box(j)+r},\
                                  {lon_box(i)-r:lon_box(i)+r})

    ; Calculate distance between four corner points of the 2 x 2 degree box
      lat_d(a) = lat_box(j)-r
      lon_d(a) = lon_box(i)-r
      lat_u(a) = lat_box(j)+r
      lon_u(a) = lon_box(i)+r

    ;=====================================================
    ; Calculate distances around each perimeter (m)
    ;=====================================================

    ; Meridional boundaries (same both sides)
      dm0 = r_earth * ((pii/180) * (2*r))
      dm0@units = "m"

    ; 'Bottom' boundary
      dz0(a) = r_earth * ((pii/180) * (2*r)) * cos((pii/180)*lat_d(a))
      dz0@units = "m"

    ; 'Top' boundary
      dz1(a) = r_earth * ((pii/180) * (2*r)) * cos((pii/180)*lat_u(a))
      dz1@units = "m"

    ; Use approximation to calculate the area of the box
      d_ave(a) = ((dz0(a) + dz1(a)) / 2) * dm0
      d_ave@units = "m~S~-2~N~"

    ; Calculate divergence averaged over the box
      div_ave(a) = wgt_areaave(div_box(a,:,:), 1.0, 1.0, 0) ; Divergence
      div_ave@units = "s~S~-1~N~"

    ; Calculate area-averaged vorticity 
      vort_ave(a) = wgt_areaave(vort_box(a,:,:), 1.0, 1.0, 0) ; Relative vorticity
      vort_ave@units = "s~S~-1~N~"

    ; Number of grid points along each side of box surrounding cyclone
      numPOINTS_x0 = dimsizes(ddx_box(a,:,:))
      numPOINTS_x  = numPOINTS_x0(1)   ; Longitude points
      numPOINTS_y  = numPOINTS_x0(0)   ; Latitude points

    ;==========================================================
    ; New arrays to hold integral sums for each side of box
    ;==========================================================

    ; Now calculate terms along the four edges of the box, for each ensemble member
    ; Hence four arrays required for each term in the budget

    ; All arrays have dimensions [49 * 13] -> ALL members on a SINGLE pressure level

    ; Relative vorticity (to compute perimeter average)
      vort_x_final0 = new((/49,numPOINTS_x/),float)
      vort_x_final1 = new((/49,numPOINTS_x/),float)
      vort_y_final0 = new((/49,numPOINTS_y/),float)
      vort_y_final1 = new((/49,numPOINTS_y/),float)

    ; Absolute vorticity (to compute perimeter average)
      avo_x_final0 = new((/49,numPOINTS_x/),float)
      avo_x_final1 = new((/49,numPOINTS_x/),float)
      avo_y_final0 = new((/49,numPOINTS_y/),float)
      avo_y_final1 = new((/49,numPOINTS_y/),float) 

    ; Zonal velocity (to compute perimeter average)
      ux_final0 = new((/49,numPOINTS_x/),float)
      ux_final1 = new((/49,numPOINTS_x/),float)
      uy_final0 = new((/49,numPOINTS_y/),float)
      uy_final1 = new((/49,numPOINTS_y/),float)

    ; Meridional velocity (to compute perimeter average)
      vx_final0 = new((/49,numPOINTS_x/),float)
      vx_final1 = new((/49,numPOINTS_x/),float)
      vy_final0 = new((/49,numPOINTS_y/),float)
      vy_final1 = new((/49,numPOINTS_y/),float)

    ;==========================================================================
    ; Integrate over each perimeter to calculate average relative vorticity
    ;==========================================================================

    ; Bottom
      do ilon = 0, numPOINTS_x-1
        vort_x_final0(a,ilon) = vort_box(a,0,ilon)
        avo_x_final0(a,ilon)  = avo_box(a,0,ilon)
        ux_final0(a,ilon)     = u_box(a,0,ilon)
        vx_final0(a,ilon)     = v_box(a,0,ilon)
      end do 

      rvx0 = avg(vort_x_final0(a,:))
      avx0 = avg(avo_x_final0(a,:))
      ux0  = avg(ux_final0(a,:))
      vx0  = avg(vx_final0(a,:))

    ; Right
      do ilat = 0, numPOINTS_y-1
        vort_y_final0(a,ilat) = vort_box(a,ilat,numPOINTS_x-1)
        avo_y_final0(a,ilat)  = avo_box(a,ilat,numPOINTS_x-1)
        uy_final0(a,ilat)     = u_box(a,ilat,numPOINTS_x-1)
        vy_final0(a,ilat)     = v_box(a,ilat,numPOINTS_x-1)
      end do

      rvy0 = avg(vort_y_final0(a,:))
      avy0 = avg(avo_y_final0(a,:))
      uy0  = avg(uy_final0(a,:))
      vy0  = avg(vy_final0(a,:))

    ; Top
      do ilon = 0, numPOINTS_x-1
        vort_x_final1(a,ilon) = vort_box(a,numPOINTS_y-1,ilon)
        avo_x_final1(a,ilon)  = avo_box(a,numPOINTS_y-1,ilon)
        ux_final1(a,ilon)     = u_box(a,numPOINTS_y-1,ilon)
        vx_final1(a,ilon)     = v_box(a,numPOINTS_y-1,ilon)
      end do

      rvx1 = avg(vort_x_final1(a,:))
      avx1 = avg(avo_x_final1(a,:))
      ux1 = avg(ux_final1(a,:))
      vx1 = avg(vx_final1(a,:))

    ; Left
      do ilat = 0, numPOINTS_y-1
        vort_y_final1(a,ilat) = vort_box(a,ilat,0)
        avo_y_final1(a,ilat)  = avo_box(a,ilat,0)
        uy_final1(a,ilat)     = u_box(a,ilat,0)
        vy_final1(a,ilat)     = v_box(a,ilat,0)
      end do

      rvy1 = avg(vort_y_final1(a,:))
      avy1 = avg(avo_y_final1(a,:))
      uy1 = avg(uy_final1(a,:))
      vy1 = avg(vy_final1(a,:))

    ; Averages using method above
      vort_avg(a) = (rvx0 + rvx1 + rvy0 + rvy1) / 4
      avo_avg(a)  = (avx0 + avx1 + avy0 + avy1) / 4
      u_avg(a)    = (ux0 + ux1 + uy0 + uy1) / 4
      v_avg(a)    = (vx0 + vx1 + vy0 + vy1) / 4

      if (u_avg(a).gt.100)
        print("Should be quitting because of box-average velocity error!!!")
        print("u [bottom] = "+ux0)
        print("u [top] = "+ux1)
        print("u [right] = "+uy0)
        print("u [left] = "+uy1)
        print("u_avg(a) = "+u_avg(a))
        exit()
      end if

;====================================================================================
; Second round of calculations (Horizontal advection, relative vorticity, tilting)
;====================================================================================

   ; Calculate terms on northern/southern boundaries
     tilt_x0     = dim_sum_n_Wrap((tilt_x_box(a,:,:) * ddx_box(a,:,:)),1)
     flux_x0     = dim_sum_n_Wrap((-avo_box(a,:,:) * v_box(a,:,:) * \
                                   ddx_box(a,:,:)),1)
     circ_x0     = dim_sum_n_Wrap((u_box(a,:,:) * ddx_box(a,:,:)),1)
     fric_x0     = dim_sum_n_Wrap((fric_x_box(a,:,:) * ddx_box(a,:,:)),1) ; PREVIOUS '0'
     tilt_b      = tilt_x0(0)              ; Tilting (south)
     tilt_t      = tilt_x0(numPOINTS_y-1)  ; Tilting (north)		  ; PREVIOUS x-1
     flux_b      = flux_x0(0)              ; Horizontal advection (south)
     flux_t      = -flux_x0(numPOINTS_y-1) ; Horizontal advection (north) ; PREVIOUS x-1
     circ_b      = circ_x0(0)              ; Circulation (south)
     circ_t      = -circ_x0(numPOINTS_y-1) ; Circulation (north)	  ; PREVIOUS x-1
     fric_b      = fric_x0(0)              ; Friction (south)
     fric_t      = -fric_x0(numPOINTS_y-1) ; Friction (north)		  ; PREVIOUS x-1

   ; Calculate terms on eastern/western boundaries
     tilt_y0    = dim_sum_n_Wrap((tilt_y_box(a,:,:) * ddy_box(a,:,:)),0)
     flux_y0    = dim_sum_n_Wrap((-avo_box(a,:,:) * u_box(a,:,:) * \
                                  ddy_box(a,:,:)),0)		          ; PREVIOUS '1'
     circ_y0    = dim_sum_n_Wrap((v_box(a,:,:) * ddy_box(a,:,:)),0)       ; PREVIOUS
     fric_y0    = dim_sum_n_Wrap((fric_y_box(a,:,:) * ddy_box(a,:,:)),0)
     tilt_l     = tilt_y0(0)              ; Tilting (west)
     tilt_r     = tilt_y0(numPOINTS_x-1)  ; Tilting (east)		  ; PREVIOUS y-1
     flux_l     = flux_y0(0)              ; Horizontal advection (west)
     flux_r     = -flux_y0(numPOINTS_x-1) ; Horizontal advection (east)	  ; PREVIOUS y-1
     circ_l     = -circ_y0(0)             ; Circulation (west)
     circ_r     = circ_y0(numPOINTS_x-1)  ; Circulation (east)            ; PREVIOUS y-1
     fric_l     = -fric_y0(0)             ; Friction (west)
     fric_r     = fric_y0(numPOINTS_x-1)  ; Friction (east)               ; PREVIOUS y-1

;===================================================
; Calculate final tendencies for all sides of box
;===================================================

   ; Advective flux
     flux_nT    = flux_t / d_ave(a)
     flux_sT	= flux_b / d_ave(a)
     flux_wT	= flux_l / d_ave(a)
     flux_eT	= flux_r / d_ave(a)

   ; Tilting 
     tilt_nT   = tilt_t / d_ave(a)
     tilt_sT   = tilt_b / d_ave(a)
     tilt_wT   = tilt_l / d_ave(a)
     tilt_eT   = tilt_r / d_ave(a)

   ; Circulation
     circ_nT   = circ_t
     circ_sT   = circ_b
     circ_wT   = circ_l
     circ_eT   = circ_r

   ; Friction
     fric_nT   = fric_t / d_ave(a)
     fric_sT   = fric_b / d_ave(a)
     fric_wT   = fric_l / d_ave(a)
     fric_eT   = fric_r / d_ave(a)

   ; Enter values into predefined arrays
     flux_north(pt,a,ct,cp)    = flux_nT		     ; Advective flux (N)	
     flux_south(pt,a,ct,cp)    = flux_sT		     ; Advective flux (S)
     flux_west(pt,a,ct,cp)     = flux_wT		     ; Advective flux (W)
     flux_east(pt,a,ct,cp)     = flux_eT     		     ; Advective flux (E)

     tilt_north(pt,a,ct,cp)    = tilt_nT                     ; Tilting (N)
     tilt_south(pt,a,ct,cp)    = tilt_sT                     ; Tilting (S)
     tilt_west(pt,a,ct,cp)     = tilt_wT                     ; Tilting (W)
     tilt_east(pt,a,ct,cp)     = tilt_eT                     ; Tilting (E)

     fric_north(pt,a,ct,cp)    = fric_nT                     ; Friction (N)
     fric_south(pt,a,ct,cp)    = fric_sT                     ; Friction (S)
     fric_west(pt,a,ct,cp)     = fric_wT                     ; Friction (W)
     fric_east(pt,a,ct,cp)     = fric_eT                     ; Friction (E)

     circ_north(pt,a,ct,cp)    = circ_nT                     ; Circulation (N)
     circ_south(pt,a,ct,cp)    = circ_sT                     ; Circulation (S)
     circ_west(pt,a,ct,cp)     = circ_wT                     ; Circulation (W)
     circ_east(pt,a,ct,cp)     = circ_eT                     ; Circulation (E)

   ; Tidy up (delete variables from within lat/lon loops before next iteration)
     delete([/flux_nT, tilt_nT, circ_nT, fric_nT/])
     delete([/flux_sT, tilt_sT, circ_sT, fric_sT/])
     delete([/flux_wT, tilt_wT, circ_wT, fric_wT/])
     delete([/flux_eT, tilt_eT, circ_eT, fric_eT/])
     delete([/vort_x_final0, vort_x_final1, vort_y_final0, vort_y_final1/])
     delete([/avo_x_final0,  avo_x_final1,  avo_y_final0,  avo_y_final1/])
     delete([/uy_final0,     uy_final1,     ux_final0,     ux_final1/])
     delete([/vy_final0,     vy_final1,     vx_final0,     vx_final1/])
     delete([/numPOINTS_x0, numPOINTS_x, numPOINTS_y/])

     end do     ; End latitude loop (do j = 0, dimsizes(lat_box)-1)
    end do	; End longitude loop (do i = 0, dimsizes(lon_box)-1)

    cp                   = cp + 1       ; Counter variable (pressure)

  ; Tidy up (delete variables from within pressure level loop)
    delete([/u_plane, v_plane, tilt_x_plane, tilt_y_plane/])
    delete([/fric_x_plane, fric_y_plane,rho,rho_a/])
    delete([/lat_plane, lon_plane, vort_plane, avo_plane, div_plane/])
    delete([/ddx_box, ddy_box, u_box, v_box, fric_x_box, fric_y_box/])
    delete([/tilt_x_box, tilt_y_box, vort_box, avo_box, div_box/])
    delete([/div_ave,vort_ave,vort_avg,avo_avg,u_avg,v_avg/])

   end do ; End pressure level loop (do ilev = 0, plevs-1)

   ct                  = ct + 1  ; Counter variable (time)

  end do                ; End time loop (do it = 0, times-1)

; Tidy up (delete variables from within input file loop before next iteration)
  delete([/lon,lat,lon1,lat1,u,v,vort,div,f0,avo,t,z,w/])
  delete([/time,times,utc_date,year,month,day,hour,minute,\
	   second,date_str,time_str,out_str/])
  delete([/lon1F,lat1F,lonF,latF,hyb,strx,stry,dxF,dyF,rho1/])

  end do      ; End input file loop (do nf = 0, numINPUT-1)

  pt = pt + 1 ; Counter variable (ensemble member)

 end do	      ; End ensemble member loop (do en = p1, p2, int)

;=====================================================
; Calculate term averages over all ensemble members
;=====================================================

   print_clock("Finished ensemble calculations! Now averaging over all members...")

 ; Advective flux
   flux_nf = new((/2,numTIMES,plevs/),"float")
   flux_sf = new((/2,numTIMES,plevs/),"float")
   flux_wf = new((/2,numTIMES,plevs/),"float")
   flux_ef = new((/2,numTIMES,plevs/),"float")

 ; Tilting
   tilt_nf = new((/2,numTIMES,plevs/),"float")
   tilt_sf = new((/2,numTIMES,plevs/),"float")
   tilt_wf = new((/2,numTIMES,plevs/),"float")
   tilt_ef = new((/2,numTIMES,plevs/),"float")

 ; Circulation
   circ_nf = new((/2,numTIMES,plevs/),"float")
   circ_sf = new((/2,numTIMES,plevs/),"float")
   circ_wf = new((/2,numTIMES,plevs/),"float")
   circ_ef = new((/2,numTIMES,plevs/),"float")

 ; Friction
   fric_nf = new((/2,numTIMES,plevs/),"float")
   fric_sf = new((/2,numTIMES,plevs/),"float")
   fric_wf = new((/2,numTIMES,plevs/),"float")
   fric_ef = new((/2,numTIMES,plevs/),"float")

 ; Divide by 1000000 [following 'tc_ens_xz_circ']
   circ_north  = circ_north / 1000000
   circ_south  = circ_south / 1000000
   circ_west   = circ_west  / 1000000	
   circ_east   = circ_east  / 1000000

   do en = 0, 1    ; Loop over ensemble members

  ; Advective flux 
    flux_nf(en,:,:) = dim_avg_n_Wrap(flux_north(en,:,:,:),0)
    flux_sf(en,:,:) = dim_avg_n_Wrap(flux_south(en,:,:,:),0)
    flux_wf(en,:,:) = dim_avg_n_Wrap(flux_west(en,:,:,:),0)
    flux_ef(en,:,:) = dim_avg_n_Wrap(flux_east(en,:,:,:),0)
 
  ; Tilting
    tilt_nf(en,:,:) = dim_avg_n_Wrap(tilt_north(en,:,:,:),0)
    tilt_sf(en,:,:) = dim_avg_n_Wrap(tilt_south(en,:,:,:),0)
    tilt_wf(en,:,:) = dim_avg_n_Wrap(tilt_west(en,:,:,:),0)
    tilt_ef(en,:,:) = dim_avg_n_Wrap(tilt_east(en,:,:,:),0)

  ; Circulation
    circ_nf(en,:,:) = dim_avg_n_Wrap(circ_north(en,:,:,:),0)
    circ_sf(en,:,:) = dim_avg_n_Wrap(circ_south(en,:,:,:),0)
    circ_wf(en,:,:) = dim_avg_n_Wrap(circ_west(en,:,:,:),0)
    circ_ef(en,:,:) = dim_avg_n_Wrap(circ_east(en,:,:,:),0)

  ; Friction
    fric_nf(en,:,:) = dim_avg_n_Wrap(fric_north(en,:,:,:),0)
    fric_sf(en,:,:) = dim_avg_n_Wrap(fric_south(en,:,:,:),0)
    fric_wf(en,:,:) = dim_avg_n_Wrap(fric_west(en,:,:,:),0)
    fric_ef(en,:,:) = dim_avg_n_Wrap(fric_east(en,:,:,:),0)

   end do	   ; End loop over ensemble members

;=========================================
; Integrate over Nepartak's life cycle
;=========================================

 ; Each calculation produces an array of 19 values [pressure levels]
 ; We do this calculation [numTIMES] times
 ; Create arrays with [numTIMES] * 19 elements

   flux_intNA = new((/2,numTIMES,plevs/),float)
   flux_intSA = new((/2,numTIMES,plevs/),float)
   flux_intWA = new((/2,numTIMES,plevs/),float)
   flux_intEA = new((/2,numTIMES,plevs/),float)

   tilt_intNA = new((/2,numTIMES,plevs/),float)
   tilt_intSA = new((/2,numTIMES,plevs/),float)
   tilt_intWA = new((/2,numTIMES,plevs/),float)
   tilt_intEA = new((/2,numTIMES,plevs/),float)

   tot_intNA  = new((/2,numTIMES,plevs/),float)
   tot_intSA  = new((/2,numTIMES,plevs/),float)
   tot_intWA  = new((/2,numTIMES,plevs/),float)
   tot_intEA  = new((/2,numTIMES,plevs/),float)

   fric_intNA = new((/2,numTIMES,plevs/),float)
   fric_intSA = new((/2,numTIMES,plevs/),float)
   fric_intWA = new((/2,numTIMES,plevs/),float)
   fric_intEA = new((/2,numTIMES,plevs/),float)

   csize     = dimsizes(flux_intNA)
   csize0    = csize(1)

   time_int = new(csize0,"integer")
   time_int = ispan(0,csize0-1,1)
   time_int!0 = "time"

   flux_intNA@units = "10~S~-5~N~ s~S~-1~N~"
   flux_intSA@units = "10~S~-5~N~ s~S~-1~N~"
   flux_intWA@units = "10~S~-5~N~ s~S~-1~N~"
   flux_intEA@units = "10~S~-5~N~ s~S~-1~N~"

   tilt_intNA@units = "10~S~-5~N~ s~S~-1~N~"
   tilt_intSA@units = "10~S~-5~N~ s~S~-1~N~"
   tilt_intWA@units = "10~S~-5~N~ s~S~-1~N~"
   tilt_intEA@units = "10~S~-5~N~ s~S~-1~N~"

   tot_intNA@units  = "10~S~-5~N~ s~S~-1~N~"
   tot_intSA@units  = "10~S~-5~N~ s~S~-1~N~"
   tot_intWA@units  = "10~S~-5~N~ s~S~-1~N~"
   tot_intEA@units  = "10~S~-5~N~ s~S~-1~N~"

   fric_intNA@units = "10~S~-5~N~ s~S~-1~N~"
   fric_intSA@units = "10~S~-5~N~ s~S~-1~N~"
   fric_intWA@units = "10~S~-5~N~ s~S~-1~N~"
   fric_intEA@units = "10~S~-5~N~ s~S~-1~N~"

   print_clock("Starting plots...")

   do en = 0, 1	         ; Loop over ensemble members [em01, em05]
    do ic = 0, csize0-1  ; 1 h integrations over the whole of Nepartak's life cycle

     if (ic.lt.csize0-1) then
      cs = ic            ; Start index
      cf = cs + 1        ; Finish index
     else
      cs = ic 
      cf = ic
     end if

      f_str = sprinti("%0.3i",cs)+"h"    ; Formatted string for output file

      flux_intNA(en,ic,:) = dim_sum_n_Wrap((flux_nf(en,cs:cf,:) * 3600), 0)
      flux_intSA(en,ic,:) = dim_sum_n_Wrap((flux_sf(en,cs:cf,:) * 3600), 0)
      flux_intWA(en,ic,:) = dim_sum_n_Wrap((flux_wf(en,cs:cf,:) * 3600), 0)
      flux_intEA(en,ic,:) = dim_sum_n_Wrap((flux_ef(en,cs:cf,:) * 3600), 0)

      tilt_intNA(en,ic,:) = dim_sum_n_Wrap((tilt_nf(en,cs:cf,:) * 3600), 0)
      tilt_intSA(en,ic,:) = dim_sum_n_Wrap((tilt_sf(en,cs:cf,:) * 3600), 0)
      tilt_intWA(en,ic,:) = dim_sum_n_Wrap((tilt_wf(en,cs:cf,:) * 3600), 0)
      tilt_intEA(en,ic,:) = dim_sum_n_Wrap((tilt_ef(en,cs:cf,:) * 3600), 0)

      fric_intNA(en,ic,:) = dim_sum_n_Wrap((fric_nf(en,cs:cf,:) * 3600), 0)
      fric_intSA(en,ic,:) = dim_sum_n_Wrap((fric_sf(en,cs:cf,:) * 3600), 0)
      fric_intWA(en,ic,:) = dim_sum_n_Wrap((fric_wf(en,cs:cf,:) * 3600), 0)
      fric_intEA(en,ic,:) = dim_sum_n_Wrap((fric_ef(en,cs:cf,:) * 3600), 0)

      flux_intNA(en,ic,:) = flux_intNA(en,ic,:) * (10 ^ 5)
      flux_intSA(en,ic,:) = flux_intSA(en,ic,:) * (10 ^ 5)
      flux_intWA(en,ic,:) = flux_intWA(en,ic,:) * (10 ^ 5)
      flux_intEA(en,ic,:) = flux_intEA(en,ic,:) * (10 ^ 5)

      tilt_intNA(en,ic,:) = tilt_intNA(en,ic,:) * (10 ^ 5)
      tilt_intSA(en,ic,:) = tilt_intSA(en,ic,:) * (10 ^ 5)
      tilt_intWA(en,ic,:) = tilt_intWA(en,ic,:) * (10 ^ 5)
      tilt_intEA(en,ic,:) = tilt_intEA(en,ic,:) * (10 ^ 5)

      fric_intNA(en,ic,:) = fric_intNA(en,ic,:) * (10 ^ 5)
      fric_intSA(en,ic,:) = fric_intSA(en,ic,:) * (10 ^ 5)
      fric_intWA(en,ic,:) = fric_intWA(en,ic,:) * (10 ^ 5)
      fric_intEA(en,ic,:) = fric_intEA(en,ic,:) * (10 ^ 5)

    end do               ; End time loop (do ic = 0, csize0-1)
   end do		 ; End ensemble member loop (do en = 0, 1)

   tot_intNA = flux_intNA + tilt_intNA + fric_intNA
   tot_intSA = flux_intSA + tilt_intSA + fric_intSA
   tot_intWA = flux_intWA + tilt_intWA + fric_intWA
   tot_intEA = flux_intEA + tilt_intEA + fric_intEA

 ; Add metadata before time-height plotting
   tot_intNA!0	   = "ens"
   tot_intNA&ens   = ens_arr(p1:p2:int)
   tot_intNA!1     = "time"
   tot_intNA&time  = time_int
   tot_intNA!2     = "lev"
   tot_intNA&lev   = pres

   tot_intSA!0     = "ens"
   tot_intSA&ens   = ens_arr(p1:p2:int)
   tot_intSA!1     = "time"
   tot_intSA&time  = time_int
   tot_intSA!2     = "lev"
   tot_intSA&lev   = pres

   tot_intWA!0     = "ens"
   tot_intWA&ens   = ens_arr(p1:p2:int)
   tot_intWA!1     = "time"
   tot_intWA&time  = time_int
   tot_intWA!2     = "lev"
   tot_intWA&lev   = pres

   tot_intEA!0     = "ens"
   tot_intEA&ens   = ens_arr(p1:p2:int)
   tot_intEA!1     = "time"
   tot_intEA&time  = time_int
   tot_intEA!2     = "lev"
   tot_intEA&lev   = pres

   flux_intNA!0	  = "ens"
   flux_intNA&ens  = ens_arr(p1:p2:int)
   flux_intNA!1    = "time"
   flux_intNA&time = time_int
   flux_intNA!2    = "lev"
   flux_intNA&lev  = pres

   flux_intSA!0   = "ens"
   flux_intSA&ens  = ens_arr(p1:p2:int)
   flux_intSA!1    = "time"
   flux_intSA&time = time_int
   flux_intSA!2    = "lev"
   flux_intSA&lev  = pres

   flux_intWA!0   = "ens"
   flux_intWA&ens  = ens_arr(p1:p2:int)
   flux_intWA!1    = "time"
   flux_intWA&time = time_int
   flux_intWA!2    = "lev"
   flux_intWA&lev  = pres

   flux_intEA!0   = "ens"
   flux_intEA&ens  = ens_arr(p1:p2:int)
   flux_intEA!1    = "time"
   flux_intEA&time = time_int
   flux_intEA!2    = "lev"
   flux_intEA&lev  = pres

   tilt_intNA!0	  = "ens"
   tilt_intNA&ens  = ens_arr(p1:p2:int)
   tilt_intNA!1    = "time"
   tilt_intNA&time = time_int
   tilt_intNA!2    = "lev"
   tilt_intNA&lev  = pres

   tilt_intSA!0   = "ens"
   tilt_intSA&ens  = ens_arr(p1:p2:int)
   tilt_intSA!1    = "time"
   tilt_intSA&time = time_int
   tilt_intSA!2    = "lev"
   tilt_intSA&lev  = pres

   tilt_intWA!0   = "ens"
   tilt_intWA&ens  = ens_arr(p1:p2:int)
   tilt_intWA!1    = "time"
   tilt_intWA&time = time_int
   tilt_intWA!2    = "lev"
   tilt_intWA&lev  = pres

   tilt_intEA!0   = "ens"
   tilt_intEA&ens  = ens_arr(p1:p2:int)
   tilt_intEA!1    = "time"
   tilt_intEA&time = time_int
   tilt_intEA!2    = "lev"
   tilt_intEA&lev  = pres

   fric_intNA!0    = "ens"
   fric_intNA&ens  = ens_arr(p1:p2:int)
   fric_intNA!1    = "time"
   fric_intNA&time = time_int
   fric_intNA!2    = "lev"
   fric_intNA&lev  = pres	;;; CHANGE TO MODEL LEVELS (?) ;;;

   fric_intSA!0    = "ens"
   fric_intSA&ens  = ens_arr(p1:p2:int)
   fric_intSA!1    = "time"
   fric_intSA&time = time_int
   fric_intSA!2    = "lev"
   fric_intSA&lev  = pres       ;;; CHANGE TO MODEL LEVELS (?) ;;;

   fric_intWA!0    = "ens"
   fric_intWA&ens  = ens_arr(p1:p2:int)
   fric_intWA!1    = "time"
   fric_intWA&time = time_int
   fric_intWA!2    = "lev"
   fric_intWA&lev  = pres       ;;; CHANGE TO MODEL LEVELS (?) ;;;

   fric_intEA!0    = "ens"
   fric_intEA&ens  = ens_arr(p1:p2:int)
   fric_intEA!1    = "time"
   fric_intEA&time = time_int
   fric_intEA!2    = "lev"
   fric_intEA&lev  = pres       ;;; CHANGE TO MODEL LEVELS (?) ;;;

;======================================
; Reorder dimensions before plotting
;======================================

 ; Define new arrays for final plotting 
   circ_intN = new((/plevs,2,numTIMES/),float)
   circ_intS = new((/plevs,2,numTIMES/),float)
   circ_intW = new((/plevs,2,numTIMES/),float)
   circ_intE = new((/plevs,2,numTIMES/),float)

   flux_intN = new((/plevs,2,numTIMES/),float)
   flux_intS = new((/plevs,2,numTIMES/),float)
   flux_intW = new((/plevs,2,numTIMES/),float)
   flux_intE = new((/plevs,2,numTIMES/),float)

   tilt_intN = new((/plevs,2,numTIMES/),float)
   tilt_intS = new((/plevs,2,numTIMES/),float)
   tilt_intW = new((/plevs,2,numTIMES/),float)
   tilt_intE = new((/plevs,2,numTIMES/),float)

   tot_intN  = new((/plevs,2,numTIMES/),float)
   tot_intS  = new((/plevs,2,numTIMES/),float)
   tot_intW  = new((/plevs,2,numTIMES/),float)
   tot_intE  = new((/plevs,2,numTIMES/),float)

   fric_intN = new((/plevs,2,numTIMES/),float)
   fric_intS = new((/plevs,2,numTIMES/),float)
   fric_intW = new((/plevs,2,numTIMES/),float)
   fric_intE = new((/plevs,2,numTIMES/),float)

 ; Add metadata to circulation array
   circ_nf!0     = "ens"
   circ_nf&ens   = ens_arr(p1:p2:int)
   circ_nf!1     = "time"
   circ_nf&time  = time_int
   circ_nf!2     = "lev"
   circ_nf&lev   = pres

   circ_sf!0     = "ens"
   circ_sf&ens   = ens_arr(p1:p2:int)
   circ_sf!1     = "time"
   circ_sf&time  = time_int
   circ_sf!2     = "lev"
   circ_sf&lev   = pres

   circ_wf!0     = "ens"
   circ_wf&ens   = ens_arr(p1:p2:int)
   circ_wf!1     = "time"
   circ_wf&time  = time_int
   circ_wf!2     = "lev"
   circ_wf&lev   = pres

   circ_ef!0     = "ens"
   circ_ef&ens   = ens_arr(p1:p2:int)
   circ_ef!1     = "time"
   circ_ef&time  = time_int
   circ_ef!2     = "lev"
   circ_ef&lev   = pres

 ; Reorder from [time * pressure] to [pressure * time]
   do en = 0, 1   ; Loop over ensemble members

    circ_intN(:,en,:) = circ_nf(lev|:,ens|en,time|:)
    circ_intS(:,en,:) = circ_sf(lev|:,ens|en,time|:)
    circ_intW(:,en,:) = circ_wf(lev|:,ens|en,time|:)
    circ_intE(:,en,:) = circ_ef(lev|:,ens|en,time|:)

    flux_intN(:,en,:) = flux_intNA(lev|:,ens|en,time|:)
    flux_intS(:,en,:) = flux_intSA(lev|:,ens|en,time|:)
    flux_intW(:,en,:) = flux_intWA(lev|:,ens|en,time|:)
    flux_intE(:,en,:) = flux_intEA(lev|:,ens|en,time|:)

    tilt_intN(:,en,:) = tilt_intNA(lev|:,ens|en,time|:)
    tilt_intS(:,en,:) = tilt_intSA(lev|:,ens|en,time|:)
    tilt_intW(:,en,:) = tilt_intWA(lev|:,ens|en,time|:)
    tilt_intE(:,en,:) = tilt_intEA(lev|:,ens|en,time|:)

    tot_intN(:,en,:)  = tot_intNA(lev|:,ens|en,time|:)
    tot_intS(:,en,:)  = tot_intSA(lev|:,ens|en,time|:)
    tot_intW(:,en,:)  = tot_intWA(lev|:,ens|en,time|:)
    tot_intE(:,en,:)  = tot_intEA(lev|:,ens|en,time|:)

    fric_intN(:,en,:) = fric_intNA(lev|:,ens|en,time|:)
    fric_intS(:,en,:) = fric_intSA(lev|:,ens|en,time|:)
    fric_intW(:,en,:) = fric_intWA(lev|:,ens|en,time|:)
    fric_intE(:,en,:) = fric_intEA(lev|:,ens|en,time|:)

   end do	  ; End ensemble member loop (do en = 0, 1)

;==============================================
; Calculate running average for budget terms
;==============================================

 ; Running average over time dimension (2)
 ; User input: only calculates average when "run" = 0 when calling script
 ; User input: "nave" [number of points used to calculate running average]

   if (run.eq.0) then
    tot_intN  = runave_n_Wrap(tot_intN, nave, 0, 2)
    tot_intS  = runave_n_Wrap(tot_intS, nave, 0, 2)
    tot_intW  = runave_n_Wrap(tot_intW, nave, 0, 2)
    tot_intE  = runave_n_Wrap(tot_intE, nave, 0, 2)

    flux_intN = runave_n_Wrap(hadv_intN, nave, 0, 2)
    flux_intS = runave_n_Wrap(hadv_intS, nave, 0, 2)
    flux_intW = runave_n_Wrap(hadv_intW, nave, 0, 2)
    flux_intE = runave_n_Wrap(hadv_intE, nave, 0, 2)

    tilt_intN = runave_n_Wrap(tilt_intN, nave, 0, 2)
    tilt_intS = runave_n_Wrap(tilt_intS, nave, 0, 2)
    tilt_intW = runave_n_Wrap(tilt_intW, nave, 0, 2)
    tilt_intE = runave_n_Wrap(tilt_intE, nave, 0, 2)

    fric_intN = runave_n_Wrap(fric_intN, nave, 0, 2)
    fric_intS = runave_n_Wrap(fric_intS, nave, 0, 2)
    fric_intW = runave_n_Wrap(fric_intW, nave, 0, 2)
    fric_intE = runave_n_Wrap(fric_intE, nave, 0, 2)

    output   = "$sam/nepartak/images/cart_coords/ens_"+r1+"_xz_cbgtF_sidesF_"+nave+"ave"
   else
    output   = "$sam/nepartak/images/cart_coords/ens_"+r1+"_xzF_cbgtF_sidesF_"+\
    	       ens_arr(p1)+"_"+ens_arr(p2)
   end if
 
;====================================================
; Open workstation for plotting (time-height plots)
;====================================================

     wks = gsn_open_wks(opt,output)

;======================================
; Create an x-y time series plot
;======================================

     panel = new(8,graphic)
     gsn_define_colormap(wks,"BlueDarkRed18")

   ; Circulation budget terms
     opts_bdgt                             = True
     opts_bdgt@cnFillOn                    = True
     opts_bdgt@cnLineLabelInterval         = 2.0
     opts_bdgt@cnLineLabelFontHeightF      = 0.012
     opts_bdgt@cnLineLabelBackgroundColor  = "transparent"
     opts_bdgt@cnLineLabelPlacementMode    = "constant"
     opts_bdgt@cnLinesOn                   = False
     opts_bdgt@cnInfoLabelOn               = False
     opts_bdgt@cnLevelSelectionMode        = "ExplicitLevels"
     opts_bdgt@cnLevels                    = (/-1.0, -0.5, -0.2, -0.1, -0.05, \
                                                0.05, 0.1, 0.2, 0.5, 1.0/)
     opts_bdgt@cnFillColors                = (/3,4,6,8,10,0,12,14,16,18,19/)
     opts_bdgt@gsnPaperOrientation         = "landscape"
     opts_bdgt@tiMainFontHeightF           = 0.0125
     opts_bdgt@gsnLeftString               = ""
     opts_bdgt@gsnRightString              = ""

     opts_bdgt@gsnMaximize                 = True
     opts_bdgt@tmXTOn                      = "False" ; Turn off top x-axis TM
     opts_bdgt@tmYROn                      = "False" ; Turn off right y-axis TM
     opts_bdgt@tiYAxisString               = "Pressure (hPa)"
     opts_bdgt@tiYAxisFontHeightF          = 0.030   ; Larger y-axis title
     opts_bdgt@tiXAxisString               = "Time"  ; X-axis title
     opts_bdgt@tiXAxisFontHeightF          = 0.030   ; Larger x-axis title
     opts_bdgt@vpWidthF                    = 0.6     ; Plot width
     opts_bdgt@vpHeightF                   = 0.5     ; Plot height
     opts_bdgt@tmXBMode                    = "Explicit"
     opts_bdgt@tmXBValues                  = (/0,23,47,71/)
     opts_bdgt@tmXBLabels                  = (/"00Z 5 Jul","00Z 6 Jul",\
     					       "00Z 7 Jul","00Z 8 Jul"/)
     opts_bdgt@tmXBLabelFontHeightF        = 0.025
     opts_bdgt@tmXBMinorValues             = ispan(2,78,3) ; Minor tickmarks (every 3h)
     opts_bdgt@trYMaxF                     = 1000          ; Ignore BL
     opts_bdgt@trYMinF                     = 150           ; Include upper troposphere

     opts_bdgt@tiMainString                = ""
     opts_bdgt@gsnDraw                     = False   ; Do not draw the plot
     opts_bdgt@gsnFrame                    = False   ; Do not advance the frame
     opts_bdgt@lbLabelBarOn                = False   ; Keep label bar on

   ; Additional resources for black contours of circulation (LHS)
     opts_cont                             = True           ; New resources
     opts_cont@cnFillOn                    = False          ; Turn off filled contours
     opts_cont@cnLineColor                 = "black"        ; Line colour
     opts_cont@gsnDraw                     = False          ; Do not draw the plot
     opts_cont@gsnFrame                    = False          ; Do not advance the frame
     opts_cont@cnLevelSelectionMode        = "ExplicitLevels"
     opts_cont@cnLevels                     = (/14.0, 16.0, 18.0, 20.0/)
     opts_cont@cnLineThicknessF             = 2             ; Contour line thickness
     opts_cont@cnLineLabelsOn               = False         ; Turn off contour labels
     opts_cont@gsnContourNegLineDashPattern = 2             ; Negative contours dashed
     opts_cont@gsnContourPosLineDashPattern = 0             ; Positive contours solid
     opts_cont@cnInfoLabelOn                = False         ; Turn off info label
     circ_intN@units                        = ""            ; Turn off string above plot

;=================
; Plot the data
;=================

     flux_north0   = gsn_csm_pres_hgt(wks,flux_intN(:,0,:),opts_bdgt)
     flux_south0   = gsn_csm_pres_hgt(wks,flux_intS(:,0,:),opts_bdgt)
     flux_west0    = gsn_csm_pres_hgt(wks,flux_intW(:,0,:),opts_bdgt)
     flux_east0    = gsn_csm_pres_hgt(wks,flux_intE(:,0,:),opts_bdgt)

     flux_north1   = gsn_csm_pres_hgt(wks,flux_intN(:,1,:),opts_bdgt)
     flux_south1   = gsn_csm_pres_hgt(wks,flux_intS(:,1,:),opts_bdgt)
     flux_west1    = gsn_csm_pres_hgt(wks,flux_intW(:,1,:),opts_bdgt)
     flux_east1    = gsn_csm_pres_hgt(wks,flux_intE(:,1,:),opts_bdgt)

     panel(0)  = flux_north0
     panel(1)  = flux_south0
     panel(2)  = flux_west0
     panel(3)  = flux_east0
     panel(4)  = flux_north1
     panel(5)  = flux_south1
     panel(6)  = flux_west1
     panel(7)  = flux_east1

     optsP                       = True
     optsP@gsnFrame              = False    ; Do not advance the frame
     optsP@gsnPanelLabelBar      = True     ; Turn on panel labelbar
     optsP@pmLabelBarWidthF      = 0.6
     optsP@pmLabelBarHeightF     = 0.15
     optsP@lbLabelFontHeightF    = 0.01
     optsP@lbPerimOn             = False
     optsP@lbLabelFont           = "Helvetica"
     optsP@lbTitleString = "Circulation tendency (10~S~-5~N~ s~S~-1~N~)"
     optsP@lbTitleFontHeightF    = 0.015
     optsP@lbTitleDirection      = "Across"
     optsP@lbTitlePosition       = "Bottom"
     optsP@txString              = "Nepartak: advective flux divergence"
     optsP@gsnPanelFigureStrings = (/"a) "+ens_arr(0)+" (North)", \
                                     "b) "+ens_arr(0)+" (Sorth)", \
                                     "c) "+ens_arr(0)+" (West)", \
				     "d) "+ens_arr(0)+" (East)", \
                                     "e) "+ens_arr(1)+" (North)", \
                                     "f) "+ens_arr(1)+" (South)", \
                                     "g) "+ens_arr(1)+" (West)", \
				     "h) "+ens_arr(1)+" (East)"/)
     optsP@gsnMaximize           = True
     optsP@gsnPanelTop           = 0.98
     optsP@gsnPanelBottom        = 0.02
     optsP@gsnPaperOrientation	 = "landscape"
     optsP@amJust                = "TopLeft"
     optsP@gsnPanelFigureStringsFontHeightF = 0.01 ; Label size (default 0.01)

     gsn_panel(wks,panel,(/2,4/),optsP)              ; Draw as a single plot
     frame(wks)
end