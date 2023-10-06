; Produce a box and whisker plot of theta-e for non-VC and VC composites 

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_rm.ncl"

; ncl opt=\"x11\" r0=0.5 sc=5 sub=1 vc_thetae_boxplot.ncl

; 'opt'   = output file format (pdf, png, x11, ...)
; 'r0'    = radius of box to calculate area average
; 'sc'    = storm centre calculation method (1-5); not currently perfect 
; 'sub'   = plot subset of VC forecasts (1) or comparison of VC and non-VC forecasts (0)

begin

;==============================
; Prelash (reading in files)
;==============================

; Read in data for the VC and no-VC groups                                             
  dat_novc  = asciiread("$ar/text/novc_dat.txt",16,"string")
  dat_vc    = asciiread("$ar/text/vc_dat.txt",29,"string")

  ens0_novc = asciiread("$ar/text/novc_sim.txt",16,"string")
  ens0_vc   = asciiread("$ar/text/vc_sim.txt",29,"string")

; Start and end times for the basic VC and non-VC datasets (add option for user interaction)
  ts0_novc  = new(dimsizes(dat_novc), integer)
  ts0_novc  = 30
  ts0_vc    = new(dimsizes(dat_vc), integer)
  ts0_vc    = 30

  tf0_novc  = ts0_novc + 60
  tf0_vc    = ts0_vc + 60

; Now create new arrays to hold all data for both VC phases                            
  var_novc  = new( (/16, 61/), "float")
  var_vc    = new( (/29, 61/), "float")

  vsize     = (/ dimsizes(var_novc), dimsizes(var_vc)/)

; Add metadata 
  var_novc!0    = "sim"
  var_novc!1    = "time"

  var_vc!0      = "sim"
  var_vc!1      = "time"

;================================================================
; Also read in data for the VC subset of the VC forecast data 
;================================================================

; 27/11/2019 --> add the subset of times on which VCs are occurring, during the VC forecasts
  dat_vcs = (/"02T12","02T12","02T12","02T12","02T12","02T12","02T12","02T12","02T12",\
              "03T00","03T00","03T00","03T00",\
              "03T12","03T12","03T12","03T12","03T12","03T12","03T12","03T12",\
              "04T00","04T00","04T00","04T00","04T00","04T00","04T00","04T00"/)

  ts0_vcs = (/46, 67, 53, 39, 66, 57, 34, 62, 44,\
              51, 58, 62, 41,\
              27, 64, 48, 28, 40, 55, 33, 36,\
              59, 59, 29, 37, 52, 40, 48, 35/)

  tf0_vcs = (/86, 94, 86, 66, 96, 96, 68, 96, 90,\
              85, 90, 95, 91,\
              96, 96, 77, 66, 91, 90, 72, 74,\
              91, 88, 82, 66, 76, 71, 85, 64/)

; Calculate how many times we need                                                             
  nsim    = (tf0_vcs - ts0_vcs) + 1
  ntot    = sum(nsim)

; Sum of each of the values in 'nsim' for use when filling array in loop below                 
  ns      = new(dimsizes(nsim),"integer")

; Calculate sum of each of the values in 'nsim' for use when filling array in loop below       
  do nt = 0, dimsizes(nsim)-1
   ns(nt) = sum(nsim(0:nt) ) - 1
  end do

; Add one to the values above                                                                  
  ns1   = ns + 1

; Create 1D array to hold all VC-only data points                                              
  vc_sub   = new(ntot, "float")

;====================================================================
; Finally, read in data for the non-VC subset of the VC forecasts
;====================================================================

; First calculate total number of times per forecast                                           
  npos0   = (96 - 24) + 1

; Then work out how many non-VC times we have, using the VC data                               
  nsim_no = nsim

; Loop through forecasts and subtract from 73                                                  
  do nn   = 0, dimsizes(nsim)-1
   nsim_no(nn) = npos0 - nsim(nn)
  end do

; Finally, calculate the total number of times                                                 
  ntot_no = sum(nsim_no)

; Sum of each of the values in 'nsim' for use when filling array in loop below                 
  ns_no   = new(dimsizes(nsim_no),"integer")

; Calculate sum of each of the values in 'nsim' for use when filling array in loop below       
  do nt = 0, dimsizes(nsim_no)-1
   ns_no(nt) = sum(nsim_no(0:nt) ) - 1
  end do

; Add one to the values above
  ns1_no   = ns_no + 1

; Create 1D array to hold all non-VC-only data points
  novc_sub = new(ntot_no, "float")

; Increase memory, before reading in forecast data 
  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 10000000000
  end setvalues

;===================================== 
;=====================================                                                 
; Read in data for NO-VC forecasts                                                     
;=====================================                                                 
;=====================================

  print("Reading in non-VC forecast data...")

; Loop over all 16 forecasts
  do st = 0, vsize(0,0)-1

 ; Print out useful information to the screen
   print_clock("Working on simulation "+ens0_novc(st)+" ("+dat_novc(st)+")")
   
 ; Define path for theta-e data
   diri       = "/nfs/a319/earshar/"+dat_novc(st)+"/"+ens0_novc(st)
   fili_p     = "201607"+dat_novc(st)+"00Z_NPTK_4p4_L80_ra1t_"+ens0_novc(st)
   fili_b1    = systemfunc("cd "+diri+" ; ls "+fili_p+"_pb.nc")
   fili_b     = diri+"/"+fili_b1

 ; Initialise counter variable (resets for each forecast)
   dt         = 0

 ; Calculate number of times in file 
   numINFO    = nc_times(fili_b)
   time       = numINFO[1]
   dsize      = dimsizes(time)

 ; Create date/time string arrays for all times in file 
   month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
                     "Jul","Aug","Sep","Oct","Nov","Dec"/)

   utc_date   = cd_calendar(time, 0)

   year       = tointeger(utc_date(:,0))
   month      = tointeger(utc_date(:,1))
   day        = tointeger(utc_date(:,2))
   hour       = tointeger(utc_date(:,3))
   minute     = tointeger(utc_date(:,4))
   second     = utc_date(:,5)

 ; Correct for errors in the code (round up value of hour when minutes = 59)        
   do ct = 0, dsize-1
    if (minute(ct).gt.30) then
      hour(ct) = hour(ct)+1
    end if
   end do
   delete(ct)

 ; Create arrays to hold strings for output plots (117 elements)                    
   date_str  = new(dsize,string)
   time_str  = new(dsize,string)
   time_arr  = new(dsize,string)
   title_arr = new(dsize,string)

   date_str  = sprinti("%0.2iUTC ", hour) + \
               sprinti("%0.2i ", day) \
             + month_abbr(month)
   time_str  = sprinti("%0.2i UTC ", hour) + \
               sprinti("%0.2i ", day) \
             + month_abbr(month)
   time_arr  = sprinti("%0.2i", day) + \
               month_abbr(month) + \
               "_" + sprinti("%0.2iZ", hour)

;=============================
; Get the variables we need 
;=============================

   b = addfile(fili_b,"r")

   do it = ts0_novc(st), tf0_novc(st)

  ; Read in storm position data for entire forecast 
    input      = "$ar/text/cp_new_"+dat_novc(st)
    lat_all    = asciiread(input+"_"+ens0_novc(st)+"_lat.txt",120,"float")
    lon_all    = asciiread(input+"_"+ens0_novc(st)+"_lon.txt",120,"float")

  ; Retrieve storm centre for desired time 
    lat0       = lat_all(it-1)
    lon0       = lon_all(it-1)

  ; Coordinates of box surrounding storm centre
    t0         = lat0-r0
    t1         = lat0+r0
    n0         = lon0-r0
    n1         = lon0+r0

  ; Pressure levels 
    pres       = b->p(3:5)
    plevs      = dimsizes(pres)
    pres@units = "hPa"

  ; Read in variables (layer-averaged between 800 and 900 hPa)

  ; Temperature (K), relative humidity (%) 
    t          = b->temp(it-4,3:5,{t0:t1},{n0:n1})
    rh         = b->rh(it-4,3:5,{t0:t1},{n0:n1})

;=========================================
; Calculate theta-e using NCL functions 
;========================================= 

  ; Create array with same dimensions as 't' for mixing ratio calculation                
    tdims   = dimsizes(t)
    pres0   = conform_dims(tdims,pres,0)
    delete(tdims)

  ; Now calculate the mixing ratio from pressure, temperature and RH                     
    mixr    = mixhum_ptrh(pres0,t,rh,1)

  ; Then calculate the lifting condensation level temperature using RH                   
  ; Option '(/1,1/)' specifies units of input and output temperatures (K)                
    tlcl    = tlcl_rh_bolton(t,rh,(/1,1/))

  ; Finally calculate equivalent potential temperature                                   
  ; Option '(/0,1,0,1/)' -- input p (hPa), temp (K), mixing ratio (kg/kg), theta-e (K)   
    th_e0   = pot_temp_equiv_tlcl(pres0,t,tlcl,mixr,(/0,1,0,1/))

  ; Calculate average over the layer 
    th_e    = dim_avg_n_Wrap(th_e0(:,:,:),0)

  ; Finally, calculate area-averaged theta-e across the box surrounding the storm centre 
    var_novc(st,dt) = wgt_areaave_Wrap(th_e, 1.0, 1.0, 0)

  ; Advance counter variable
    dt = dt + 1

  ; Tidy up 
    delete([/t,rh,pres0,mixr,tlcl,th_e0,th_e/])

 ; End loop over times (do it = ts0_novc(st), tf0_novc(st) ) 
   end do

 ; Tidy up 
   delete([/time,utc_date,year,month,day,hour,minute,second/])
   delete([/date_str,time_str,time_arr,title_arr/])

; End loop over all 16 forecasts (do st = 0, vsize(0,0)-1)
  end do

;=================================
;=================================                                                     
; Read in data for VC forecasts                                                        
;=================================                                                     
;=================================

  print("Now reading in VC forecast data...")

; Loop over all 29 forecasts                                                           
  do st = 0, vsize(1,0)-1

 ; Print out useful information to the screen
   print_clock("Working on simulation "+ens0_vc(st)+" ("+dat_vc(st)+")")

 ; Define path for theta-e data                                                         
   diri       = "/nfs/a319/earshar/"+dat_vc(st)+"/"+ens0_vc(st)
   fili_p     = "201607"+dat_vc(st)+"00Z_NPTK_4p4_L80_ra1t_"+ens0_vc(st)
   fili_b1    = systemfunc("cd "+diri+" ; ls "+fili_p+"_pb.nc")
   fili_b     = diri+"/"+fili_b1

 ; Initialise counter variable (resets for each forecast)                               
   dt         = 0

 ; Calculate number of times in file                                                    
   numINFO    = nc_times(fili_b)
   time       = numINFO[1]
   dsize      = dimsizes(time)

 ; Create date/time string arrays for all times in file                                 
   month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
                     "Jul","Aug","Sep","Oct","Nov","Dec"/)

   utc_date   = cd_calendar(time, 0)

   year       = tointeger(utc_date(:,0))
   month      = tointeger(utc_date(:,1))
   day        = tointeger(utc_date(:,2))
   hour       = tointeger(utc_date(:,3))
   minute     = tointeger(utc_date(:,4))
   second     = utc_date(:,5)

 ; Correct for errors in the code (round up value of hour when minutes = 59)            
   do ct = 0, dsize-1
    if (minute(ct).gt.30) then
      hour(ct) = hour(ct)+1
    end if
   end do
   delete(ct)

 ; Create arrays to hold strings for output plots (117 elements)                        
   date_str  = new(dsize,string)
   time_str  = new(dsize,string)
   time_arr  = new(dsize,string)
   title_arr = new(dsize,string)

   date_str  = sprinti("%0.2iUTC ", hour) + \
               sprinti("%0.2i ", day) \
             + month_abbr(month)
   time_str  = sprinti("%0.2i UTC ", hour) + \
               sprinti("%0.2i ", day) \
             + month_abbr(month)
   time_arr  = sprinti("%0.2i", day) + \
               month_abbr(month) + \
               "_" + sprinti("%0.2iZ", hour)

;=============================                                                          
; Get the variables we need                                                             
;=============================                                                          

   b = addfile(fili_b,"r")

   do it = ts0_vc(st), tf0_vc(st)

  ; Read in storm position data for entire forecast                                     
    input      = "$ar/text/cp_new_"+dat_vc(st)
    lat_all    = asciiread(input+"_"+ens0_vc(st)+"_lat.txt",120,"float")
    lon_all    = asciiread(input+"_"+ens0_vc(st)+"_lon.txt",120,"float")

  ; Retrieve storm centre for desired time                                              
    lat0       = lat_all(it-1)
    lon0       = lon_all(it-1)

  ; Coordinates of box surrounding storm centre                                         
    t0         = lat0-r0
    t1         = lat0+r0
    n0         = lon0-r0
    n1         = lon0+r0

  ; Pressure levels                                                                     
    pres       = b->p(3:5)
    plevs      = dimsizes(pres)
    pres@units = "hPa"

  ; Read in variables (layer-averaged between 800 and 900 hPa)                          

  ; Temperature (K), relative humidity (%)                                              
    t          = b->temp(it-4,3:5,{t0:t1},{n0:n1})
    rh         = b->rh(it-4,3:5,{t0:t1},{n0:n1})

;=========================================                                              
; Calculate theta-e using NCL functions                                                 
;========================================= 

  ; Create array with same dimensions as 't' for mixing ratio calculation                
    tdims   = dimsizes(t)
    pres0   = conform_dims(tdims,pres,0)
    delete(tdims)

  ; Now calculate the mixing ratio from pressure, temperature and RH                     
    mixr    = mixhum_ptrh(pres0,t,rh,1)

  ; Then calculate the lifting condensation level temperature using RH                   
  ; Option '(/1,1/)' specifies units of input and output temperatures (K)                
    tlcl    = tlcl_rh_bolton(t,rh,(/1,1/))

  ; Finally calculate equivalent potential temperature                                   
  ; Option '(/0,1,0,1/)' -- input p (hPa), temp (K), mixing ratio (kg/kg), theta-e (K)   
    th_e0   = pot_temp_equiv_tlcl(pres0,t,tlcl,mixr,(/0,1,0,1/))

  ; Calculate average over the layer                                                      
    th_e    = dim_avg_n_Wrap(th_e0(:,:,:),0)

  ; Finally, calculate area-averaged theta-e across the box surrounding the storm centre 
    var_vc(st,dt) = wgt_areaave_Wrap(th_e, 1.0, 1.0, 0)

  ; Advance counter variable                                                              
    dt = dt + 1

  ; Tidy up 
    delete([/t,rh,pres0,mixr,tlcl,th_e0,th_e/])

 ; End loop over times (do it = ts0_vc(st), tf0_vc(st) ) 
   end do

 ; Tidy up 
   delete([/time,utc_date,year,month,day,hour,minute,second/])
   delete([/date_str,time_str,time_arr,title_arr/])

; End loop over all 29 forecasts (do st = 0, vsize(0,0)-1)                                
  end do

;=================================
;=================================
; VC subset (from VC forecasts)
;=================================
;=================================

  print("Now reading in VC subset data...")

; Initialise counter variable
  dt    = 0

; Loop over all 29 forecasts                                                             
  do st = 0, vsize(1,0)-1

 ; Print out useful information to the screen                                            
   print_clock("Working on simulation "+ens0_vc(st)+" ("+dat_vc(st)+")")

 ; Define path for theta-e data                                                         
   diri       = "/nfs/a319/earshar/"+dat_vc(st)+"/"+ens0_vc(st)
   fili_p     = "201607"+dat_vc(st)+"00Z_NPTK_4p4_L80_ra1t_"+ens0_vc(st)
   fili_b1    = systemfunc("cd "+diri+" ; ls "+fili_p+"_pb.nc")
   fili_b     = diri+"/"+fili_b1

 ; Calculate number of times in file                                                    
   numINFO    = nc_times(fili_b)
   time       = numINFO[1]
   dsize      = dimsizes(time)

 ; Create date/time string arrays for all times in file                                 
   month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
                     "Jul","Aug","Sep","Oct","Nov","Dec"/)

   utc_date   = cd_calendar(time, 0)

   year       = tointeger(utc_date(:,0))
   month      = tointeger(utc_date(:,1))
   day        = tointeger(utc_date(:,2))
   hour       = tointeger(utc_date(:,3))
   minute     = tointeger(utc_date(:,4))
   second     = utc_date(:,5)

 ; Correct for errors in the code (round up value of hour when minutes = 59)            
   do ct = 0, dsize-1
    if (minute(ct).gt.30) then
      hour(ct) = hour(ct)+1
    end if
   end do
   delete(ct)

 ; Create arrays to hold strings for output plots (117 elements)                        
   date_str  = new(dsize,string)
   time_str  = new(dsize,string)
   time_arr  = new(dsize,string)
   title_arr = new(dsize,string)

   date_str  = sprinti("%0.2iUTC ", hour) + \
               sprinti("%0.2i ", day) \
             + month_abbr(month)
   time_str  = sprinti("%0.2i UTC ", hour) + \
               sprinti("%0.2i ", day) \
             + month_abbr(month)
   time_arr  = sprinti("%0.2i", day) + \
               month_abbr(month) + \
               "_" + sprinti("%0.2iZ", hour)

;=============================                                                          
; Get the variables we need                                                             
;=============================                                                          

   b = addfile(fili_b,"r")

   do it = ts0_vcs(st), tf0_vcs(st)

  ; Read in storm position data for entire forecast                                     
    input      = "$ar/text/cp_new_"+dat_vc(st)
    lat_all    = asciiread(input+"_"+ens0_vc(st)+"_lat.txt",120,"float")
    lon_all    = asciiread(input+"_"+ens0_vc(st)+"_lon.txt",120,"float")

  ; Retrieve storm centre for desired time                                              
    lat0       = lat_all(it-1)
    lon0       = lon_all(it-1)

  ; Coordinates of box surrounding storm centre                                         
    t0         = lat0-r0
    t1         = lat0+r0
    n0         = lon0-r0
    n1         = lon0+r0

  ; Pressure levels                                                                     
    pres       = b->p(3:5)
    plevs      = dimsizes(pres)
    pres@units = "hPa"

  ; Read in variables (layer-averaged between 800 and 900 hPa)                          

  ; Temperature (K), relative humidity (%)                                              
    t          = b->temp(it-4,3:5,{t0:t1},{n0:n1})
    rh         = b->rh(it-4,3:5,{t0:t1},{n0:n1})

;=========================================                                              
; Calculate theta-e using NCL functions                                                 
;========================================= 

  ; Create array with same dimensions as 't' for mixing ratio calculation                  
    tdims   = dimsizes(t)
    pres0   = conform_dims(tdims,pres,0)
    delete(tdims)

  ; Now calculate the mixing ratio from pressure, temperature and RH                       
    mixr    = mixhum_ptrh(pres0,t,rh,1)

  ; Then calculate the lifting condensation level temperature using RH                     
  ; Option '(/1,1/)' specifies units of input and output temperatures (K)                  
    tlcl    = tlcl_rh_bolton(t,rh,(/1,1/))

  ; Finally calculate equivalent potential temperature                                     
  ; Option '(/0,1,0,1/)' -- input p (hPa), temp (K), mixing ratio (kg/kg), theta-e (K)     
    th_e0   = pot_temp_equiv_tlcl(pres0,t,tlcl,mixr,(/0,1,0,1/))

  ; Calculate average over the layer                                                      
    th_e    = dim_avg_n_Wrap(th_e0(:,:,:),0)

  ; Finally, calculate area-averaged theta-e across the box surrounding the storm centre   
    vc_sub(dt) = wgt_areaave_Wrap(th_e, 1.0, 1.0, 0)

  ; Advance counter variable                                                              
    dt = dt + 1

  ; Tidy up                                                                                
    delete([/t,rh,pres0,mixr,tlcl,th_e0,th_e/])

 ; End loop over times (do it = ts0_vc(st), tf0_vc(st) )                                   
   end do

 ; Tidy up                                                                                 
   delete([/time,utc_date,year,month,day,hour,minute,second/])
   delete([/date_str,time_str,time_arr,title_arr/])

; End loop over all 29 forecasts (do st = 0, vsize(0,0)-1)                                
  end do

;======================================                                                        
;======================================                                                        
; Non-VC subset (from VC forecasts)                                                           
;======================================                                                        
;======================================  

; Initialise counter variables
  dt = 0 

; Useful information
  print_clock("Now reading in non-VC subset data...")

; Loop over all 29 forecasts                                                                   
  do st = 0, vsize(1,0)-1

 ; Print out useful information to the screen                                                   
   print_clock("Working on simulation "+ens0_vc(st)+" ("+dat_vc(st)+")")

 ; Define path for theta-e data                                                                 
   diri       = "/nfs/a319/earshar/"+dat_vc(st)+"/"+ens0_vc(st)
   fili_p     = "201607"+dat_vc(st)+"00Z_NPTK_4p4_L80_ra1t_"+ens0_vc(st)
   fili_b1    = systemfunc("cd "+diri+" ; ls "+fili_p+"_pb.nc")
   fili_b     = diri+"/"+fili_b1

 ; Calculate number of times in file                                                            
   numINFO    = nc_times(fili_b)
   time       = numINFO[1]
   dsize      = dimsizes(time)

 ; Create date/time string arrays for all times in file                                         
   month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
                     "Jul","Aug","Sep","Oct","Nov","Dec"/)

   utc_date   = cd_calendar(time, 0)

   year       = tointeger(utc_date(:,0))
   month      = tointeger(utc_date(:,1))
   day        = tointeger(utc_date(:,2))
   hour       = tointeger(utc_date(:,3))
   minute     = tointeger(utc_date(:,4))
   second     = utc_date(:,5)

 ; Correct for errors in the code (round up value of hour when minutes = 59)                    
   do ct = 0, dsize-1
    if (minute(ct).gt.30) then
      hour(ct) = hour(ct)+1
    end if
   end do
   delete(ct)

 ; Create arrays to hold strings for output plots (117 elements)                                
   date_str  = new(dsize,string)
   time_str  = new(dsize,string)
   time_arr  = new(dsize,string)
   title_arr = new(dsize,string)

   date_str  = sprinti("%0.2iUTC ", hour) + \
               sprinti("%0.2i ", day) \
             + month_abbr(month)
   time_str  = sprinti("%0.2i UTC ", hour) + \
               sprinti("%0.2i ", day) \
             + month_abbr(month)
   time_arr  = sprinti("%0.2i", day) + \
               month_abbr(month) + \
               "_" + sprinti("%0.2iZ", hour)

;=============================                                                                  
; Get the variables we need                                                                     
;=============================                                                                  

   b = addfile(fili_b,"r")

   do it = 24, 96

  ; Read in storm position data for entire forecast                                             
    input      = "$ar/text/cp_new_"+dat_vc(st)
    lat_all    = asciiread(input+"_"+ens0_vc(st)+"_lat.txt",120,"float")
    lon_all    = asciiread(input+"_"+ens0_vc(st)+"_lon.txt",120,"float")

  ; Retrieve storm centre for desired time                                                      
    lat0       = lat_all(it-1)
    lon0       = lon_all(it-1)

  ; Coordinates of box surrounding storm centre                                                 
    t0         = lat0-r0
    t1         = lat0+r0
    n0         = lon0-r0
    n1         = lon0+r0

  ; Pressure levels                                                                             
    pres       = b->p(3:5)
    plevs      = dimsizes(pres)
    pres@units = "hPa"

  ; Read in variables (layer-averaged between 800 and 900 hPa)                                  

  ; Temperature (K), relative humidity (%)                                                      
    t          = b->temp(it-4,3:5,{t0:t1},{n0:n1})
    rh         = b->rh(it-4,3:5,{t0:t1},{n0:n1})

;=========================================                                                      
; Calculate theta-e using NCL functions                                                         
;=========================================                                                      

  ; Create array with same dimensions as 't' for mixing ratio calculation                       
    tdims   = dimsizes(t)
    pres0   = conform_dims(tdims,pres,0)
    delete(tdims)

  ; Now calculate the mixing ratio from pressure, temperature and RH                            
    mixr    = mixhum_ptrh(pres0,t,rh,1)

  ; Then calculate the lifting condensation level temperature using RH                          
  ; Option '(/1,1/)' specifies units of input and output temperatures (K)                       
    tlcl    = tlcl_rh_bolton(t,rh,(/1,1/))

  ; Finally calculate equivalent potential temperature                                          
  ; Option '(/0,1,0,1/)' -- input p (hPa), temp (K), mixing ratio (kg/kg), theta-e (K)          
    th_e0   = pot_temp_equiv_tlcl(pres0,t,tlcl,mixr,(/0,1,0,1/))

  ; Calculate average over the layer                                                            
    th_e    = dim_avg_n_Wrap(th_e0(:,:,:),0)

  ; Finally, calculate area-averaged theta-e across the box surrounding the storm centre        
    if (it .lt. ts0_vcs(st) .or. it .gt. tf0_vcs(st) ) then 
     novc_sub(dt) = wgt_areaave_Wrap(th_e, 1.0, 1.0, 0)
   ; Only advance counter variable if time corresponds to non-VC
     dt = dt + 1
    end if 

  ; Tidy up                                                                                     
    delete([/t,rh,pres0,mixr,tlcl,th_e0,th_e/])

 ; End loop over times (do it = ts0_vc(st), tf0_vc(st) )                                        
   end do

 ; Tidy up                                                                                      
   delete([/time,utc_date,year,month,day,hour,minute,second/])
   delete([/date_str,time_str,time_arr,title_arr/])

; End loop over all 29 forecasts (do st = 0, vsize(0,0)-1)                                      
  end do

;=============================================                                                 
; Reshape to form 1D arrays before plotting                                                    
;=============================================                                                 

; Data from all forecasts in which VCs didn't occur                                           
  novc_the = ndtooned(var_novc(:,:) )
  nvc_size = dimsizes(novc_the)

; Data from all forecasts in which VCs occurred
  vc_the   = ndtooned(var_vc(:,:) )
  vc_size  = dimsizes(vc_the)

; Data from the VC subset (all VC times during VC forecasts)
  vcs_the  = vc_sub 
  vcs_size = dimsizes(vcs_the)

; Data from the non-VC subset (non-VC times during VC forecasts)
  nvcs_the = novc_sub
  nvcs_size= dimsizes(nvcs_the)

; Print out some useful information to screen 
  printMinMax(vc_the,False)
  printMinMax(vcs_the,False)
  printMinMax(nvcs_the,False)
  printMinMax(novc_the,False)

;======================================================                                        
; Produce box plot comparing VC and no-VC composites                                           
;======================================================                                        

; Create array to hold all box plot values                                                     
; 2 VC phases                                                                                  
; 5 values per phase [minimum, lower quartile, median, upper quartile, maximum]                
  box_y    = new( (/2,5/), "float")

; Options for 'stat_dispersion' function                                                       
  opts     = True
  opts@PrintStat = False

; Calculate stats for both VC phases                                                           
  if (sub .eq. 1) then
   ptb_novc = stat_dispersion(nvcs_the, opts)
   ptb_vc   = stat_dispersion(vcs_the, opts)
  else
   ptb_novc = stat_dispersion(novc_the, opts)
   ptb_vc   = stat_dispersion(vc_the, opts)
  end if

; No-VC composite                                                                              
  box_y(0,0) = ptb_novc(2)
  box_y(0,1) = ptb_novc(6)
  box_y(0,2) = ptb_novc(8)
  box_y(0,3) = ptb_novc(10)
  box_y(0,4) = ptb_novc(14)

; VC composite                                                                                 
  box_y(1,0) = ptb_vc(2)
  box_y(1,1) = ptb_vc(6)
  box_y(1,2) = ptb_vc(8)
  box_y(1,3) = ptb_vc(10)
  box_y(1,4) = ptb_vc(14)

;========================                                                                      
; Produce the box plot                                                                         
;========================                                                                      

; Define output file path                                                                      
  if (sub .eq. 1) then
   output = "$nep/nepartak/images/thetae_boxplot_sc"+sc+"_vc_novc_"+r0+"d_sub"
  else
   output = "$nep/nepartak/images/thetae_boxplot_sc"+sc+"_vc_novc_"+r0+"deg"
  end if

; Open workstation for plotting                                                                
  wks = gsn_open_wks(opt, output)

; Create set of integer values denoting box label positions                                    
  x0  = ispan(1,2,1)

; General box plot resources                                                                   
  p_res = True
  p_res@tmXBLabels      = (/"Non-VC","VC"/)
  p_res@gsnMaximize     = True
  p_res@tiYAxisFontHeightF = 0.0225
  p_res@tiYAxisString   = "800-900 hPa layer-averaged ~F33~q~B~e"

; Move y-axis title slightly left (away from plot)
  p_res@tiYAxisOffsetXF = -0.01

; Possibly set range (27/11/2019)
  p_res@trYMaxF  = 312.0
  p_res@trYMinF  = 304.0

; Box resources                                                                                
  bres = True
  bres@boxWidth  = 0.5
  bres@boxColors = (/"royalblue1", "purple3"/)

; Line resources (increase line thickness --> 05/12/2019)
  lres = True
  lres@gsLineThicknessF   = 3.0

; Produce the boxplot                                                                          
; Rightmost dimension needs to contain box plot reference data (min,median,max,...)            
  plot = boxplot(wks, x0, box_y(:,:), bres, p_res, lres)

  draw(plot)
  frame(wks)

; Print output file path to screen
  print(output)

end