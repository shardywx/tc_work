; Produce scatter and box & whisker plots for chosen variables (v,pmin,vorticity ratio)

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"

load "$sam/ncl_func/attach_vert_axis.ncl"
load "$sam/ncl_func/attach_horiz_axis.ncl"
load "$sam/ncl_func/cartesian_axis.ncl"

; ncl opt=\"x11\" zs=\"False\" z0=500 wind=\"vtan\" box=1 mth=\"p2\" trb=0 tend=2 w0=0.3 mlev0=17
; mlev1=46 sub=0 lgd0=\"tr\" t0=1 alt=1 sc=5 r0=15 r1=30 wt=0.5 wrt=1 var=\"v\" x=\"flux\" miss=1
; grp=\"8\" pm=-1 rthr=20 rthr2=5 zthr=10 sct=0 thr=0 boxv=\"v\" v0=2.0 vc_scatter_20dec.ncl

; 'opt'   = output file format ("pdf" or "x11")
; 'mth'   = method for calculating storm centre
; 'trb'   = turn troubleshooting on (1) or off (0)
; 'tend'  = method for calculating centred differences (2 or 1)
; 'alt'   = use original (0) or alternative (1) VC phase data
; 'plt'   = all phases ("all"); main four VC phases ("vc"); VC/No-VC ("no")
; 'sub'   = use a subset of the data points (e.g. only positive tendencies)
; 't0'    = wind speed / MSLP tendency threshold (e.g. 0,1,2), if 'sub=1'  
; 'zs'    = tendencies calculated on single (True) or any (False) model levels
; 'pm'    = discard data either side of RMW (1), inside the RMW (0), or between 0 and 'x'km (-1)

begin

;=========================================
; Now read in the data from text files 
;=========================================

; Create output radius constraint string (25/11/2019)
  if (rthr2 .gt. 0) then 
   rstr_out = rthr2+"km_"+rthr+"km"
  else
   rstr_out = rthr+"km"
  end if 

  fili_vc   = (/"$ar/text/vc_sym_dat_group"+grp+".txt",\
                "$ar/text/vc_s2a_dat_group"+grp+".txt",\
                "$ar/text/vc_asym_dat_group"+grp+".txt",\
                "$ar/text/vc_a2s_dat_group"+grp+".txt"/)

; Get file size information for each VC phase using 'systemfunc'                            
  rsize0 = new(4,"string")
  do il  = 0, 3
   rsize0(il) = systemfunc("wc -l < "+fili_vc(il) )
  end do
  rsize  = toint(rsize0)
  print(rsize)

  print("Reading in data from text files...")

; Read in correct data for all VC phases 
  dat_ring  = asciiread("$ar/text/vc_sym_dat_group"+grp+".txt",(/rsize(0)/),"string")
  dat_r2m   = asciiread("$ar/text/vc_s2a_dat_group"+grp+".txt",(/rsize(1)/),"string")
  dat_mono  = asciiread("$ar/text/vc_asym_dat_group"+grp+".txt",(/rsize(2)/),"string")
  dat_m2r   = asciiread("$ar/text/vc_a2s_dat_group"+grp+".txt",(/rsize(3)/),"string")

  ens0_ring = asciiread("$ar/text/vc_sym_sim_group"+grp+".txt",(/rsize(0)/),"string")
  ens0_r2m  = asciiread("$ar/text/vc_s2a_sim_group"+grp+".txt",(/rsize(1)/),"string")
  ens0_mono = asciiread("$ar/text/vc_asym_sim_group"+grp+".txt",(/rsize(2)/),"string")
  ens0_m2r  = asciiread("$ar/text/vc_a2s_sim_group"+grp+".txt",(/rsize(3)/),"string")

  ts0_ring  = asciiread("$ar/text/vc_sym_ts_group"+grp+".txt",(/rsize(0)/),"integer")
  ts0_r2m   = asciiread("$ar/text/vc_s2a_ts_group"+grp+".txt",(/rsize(1)/),"integer")
  ts0_mono  = asciiread("$ar/text/vc_asym_ts_group"+grp+".txt",(/rsize(2)/),"integer")
  ts0_m2r   = asciiread("$ar/text/vc_a2s_ts_group"+grp+".txt",(/rsize(3)/),"integer")

  tf0_ring  = asciiread("$ar/text/vc_sym_tf_group"+grp+".txt",(/rsize(0)/),"integer")
  tf0_r2m   = asciiread("$ar/text/vc_s2a_tf_group"+grp+".txt",(/rsize(1)/),"integer")
  tf0_mono  = asciiread("$ar/text/vc_asym_tf_group"+grp+".txt",(/rsize(2)/),"integer")
  tf0_m2r   = asciiread("$ar/text/vc_a2s_tf_group"+grp+".txt",(/rsize(3)/),"integer")

; Calculate maximum number of times to read in for any one simulation, for each VC phase
  td_ring   = (tf0_ring - ts0_ring) + 1
  nts_ring  = max(td_ring)

  td_r2m    = (tf0_r2m - ts0_r2m) + 1
  nts_r2m   = max(td_r2m)

  td_mono   = (tf0_mono - ts0_mono) + 1
  nts_mono  = max(td_mono)

  td_m2r    = (tf0_m2r - ts0_m2r) + 1
  nts_m2r   = max(td_m2r)

; 12/11/2019 --> better to use total number of actual times per VC phase?
  sum_ring  = sum(td_ring)
  sum_r2m   = sum(td_r2m)
  sum_mono  = sum(td_mono)
  sum_m2r   = sum(td_m2r)

; Now create new arrays to hold all data for each VC phase
; 29/10/2019 --> added an extra element in dimension 2 (for mass flux)
  var_ring  = new( (/rsize(0), nts_ring, 5/), "float")
  var_r2m   = new( (/rsize(1), nts_r2m, 5/), "float")
  var_mono  = new( (/rsize(2), nts_mono, 5/), "float")
  var_m2r   = new( (/rsize(3), nts_m2r, 5/), "float")

; Dimension size for each VC phase array 
  dsize     = (/ rsize(0), rsize(1), \
                 rsize(2), rsize(3) /)

; Add metadata to arrays above 
  var_ring!0    = "sim"
  var_ring!1    = "time"
  var_ring!2    = "diag"
  var_ring&diag = (/"Wind speed tendency", "MSLP tendency", \
                    "RVP","RVP tendency", "Mass flux"/)

  var_r2m!0     = "sim"
  var_r2m!1     = "time"
  var_r2m!2     = "diag"
  var_r2m&diag  = (/"Wind speed tendency", "MSLP tendency", \
                    "RVP","RVP tendency", "Mass flux"/)

  var_mono!0    = "sim"
  var_mono!1    = "time"
  var_mono!2    = "diag"
  var_mono&diag = (/"Wind speed tendency", "MSLP tendency", \
                    "RVP","RVP tendency", "Mass flux"/)

  var_m2r!0     = "sim"
  var_m2r!1     = "time"
  var_m2r!2     = "diag"
  var_m2r&diag  = (/"Wind speed tendency", "MSLP tendency", \
                    "RVP","RVP tendency", "Mass flux"/)

; Number of times in each text file 
  ntimes_v  = 121
  ntimes_p  = 120
  ntimes_vc = 119

; String indicating how we calculated tendencies 
  if (tend .eq. 2) then 
   tstr = "_tend2"
  else
   tstr = ""
  end if 

;====================================================== 
; Read in 'pd' stream data to grab height level info 
;====================================================== 

; Define file path  
  diri    = "/nfs/a319/earshar/02T12/em11"
  fili_p  = "20160702T1200Z_NPTK_4p4_L80_ra1t_em11"
  fili_d0 = systemfunc("cd "+diri+" ; ls "+fili_p+"_pd*.nc")
  fili_d  = diri+"/"+fili_d0

; Read in data and extract height-level information 
  d       = addfile(fili_d,"r")
  hybD    = d->hybrid_ht(:)
  hy0     = sprintf("%0.0f",hybD(mlev0))
  hy1     = sprintf("%0.0f",hybD(mlev1))

;==================================
; Start loop over ringlike phase
;==================================

  do st = 0, dsize(0)-1 

 ; Read in data for each of our desired parameters 
 ; Tangential wind speed tendency (121 values)
 ; 18/10/2019 --> read in single-level values, if desired 

 ; Calculate tendencies using tangential wind on single model level, only 
   if (zs .eq. "True") then 
  ; Tangential wind 
    if (wind .eq. "vtan") then 
     diri_v  = "$ar/text/ml_"+dat_ring(st)+"_"+ens0_ring(st)+"_v"+z0+"_tend_"+mth+tstr+".txt"
  ; Total wind
    else
     diri_v  = "$ar/text/ml_"+dat_ring(st)+"_"+ens0_ring(st)+"_vtot"+z0+"_tend_"+mth+tstr+".txt"
    end if 
   elseif (zs .eq. "max") then 
    diri_v   = "$ar/text/ml_"+dat_ring(st)+"_"+ens0_ring(st)+"_vtan_max_p2.txt"
 ; OR, calculate tendencies using tangential wind on any model level 
   else
    diri_v   = "$ar/text/ml_"+dat_ring(st)+"_"+ens0_ring(st)+"_vtan_tend_"+mth+tstr+".txt"
   end if 

 ; Mean sea level pressure (120 values; need to calculate tendency inline)
   diri_p    = "$ar/text/cp_new_"+dat_ring(st)+"_"+ens0_ring(st)+"_slp.txt"

 ; Radial vorticity profile (119 values)
   diri_rvp  = "$ar/text/"+dat_ring(st)+"_"+ens0_ring(st)+"_sc"+sc+"_inner_core.txt"

 ; Radial vorticity profile tendency (119 values)
   diri_ten  = "$ar/text/"+dat_ring(st)+"_"+ens0_ring(st)+"_sc"+sc+"_inner_core_grad.txt"

 ; Mass flux (119 values)
   diri_flux = "$ar/text/"+dat_ring(st)+"_"+ens0_ring(st)+"_sc"+sc+"_lay1_"+hy0+"_"+hy1+\
               "_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+".txt"

;==========================================================
; Read in both wind speed tendency and MSLP (11/10/2019)
;==========================================================

  ; Wind speed tendency
  ; Check whether file exists, and if so, read data from it 
    if (.not. fileexists(diri_v) ) then 
     tend_all0 = new( ntimes_v, "float")
     print("Data missing for simulation "+ens0_ring(st)+" ("+dat_ring(st)+")" )
    else
     print("Working on RI period: "+st+" ("+ens0_ring(st)+", "+dat_ring(st)+")")
     tend_all0 = asciiread(diri_v, ntimes_v, "float")

;     if (zs .eq."max") then
;      do it = 0, ntimes_v-1
;       tp1 = min((/it+1,ntimes_v-1/))
;       tm1 = max((/it-1,0/))
;       tend_all(it) = (tend_all0(tp1) - tend_all0(tm1) ) / 2
;      end do     
;     else
      tend_all  = tend_all0(2::)
;     end if 
     delete(tend_all0)
    end if 

  ; Mean sea level pressure 
  ; Check whether file containing MSLP data exists
    if (.not. fileexists(diri_p) ) then
     mslp_all1 = new( ntimes_p, "float")
     print("Data missing for simulation "+ens0_ring(st)+" ("+dat_ring(st)+")" )
    else
     mslp_all1 = asciiread(diri_p, ntimes_p, "float")
     mslp_all0 = mslp_all1(1::)
     delete(mslp_all1)
    end if

  ; Define new array to hold MSLP values (119 values)
    mslp_all   = new( ntimes_vc, "float")

  ; Calculate MSLP tendency inline 
    do it = 0, ntimes_vc-1
     tp1 = min((/it+1,ntimes_vc-1/))
     tm1 = max((/it-1,0/))
     mslp_all(it) = ( mslp_all0(tp1) - mslp_all0(tm1) ) / 2
    end do

  ; Radial vorticity profile (and tendency)
    rvp_all     = asciiread(diri_rvp, ntimes_vc, "float")
    ten_all     = asciiread(diri_ten, ntimes_vc, "float")

  ; Finally, all data into correct arrays [wind speed + MSLP tendency; radial vorticity profile]
    npts = ( tf0_ring(st) - ts0_ring(st) )
    var_ring(st,0:npts,0) = (/tend_all( ts0_ring(st):tf0_ring(st) )/)
    var_ring(st,0:npts,1) = (/mslp_all( ts0_ring(st):tf0_ring(st) )/)
    var_ring(st,0:npts,2) = (/rvp_all( ts0_ring(st):tf0_ring(st) )/)
    var_ring(st,0:npts,3) = (/ten_all( ts0_ring(st):tf0_ring(st) )/)

    if (x .eq. "flux") then
     flux_all   = asciiread(diri_flux, ntimes_vc, "float")
     var_ring(st,0:npts,4) = (/flux_all( ts0_ring(st):tf0_ring(st) )/)
     delete(diri_flux)
    end if

  end do  

  delete([/diri_v, diri_p, diri_rvp, diri_ten/])

;===================================================
; Start loop over ringlike to monopole transition
;=================================================== 

  do st = 0, dsize(1)-1

 ; Read in data for each of our desired parameters                                          
 ; Tangential wind speed tendency (121 values)                                              
   if (zs .eq. "True") then
  ; Tangential wind
    if (wind .eq. "vtan") then 
     diri_v = "$ar/text/ml_"+dat_r2m(st)+"_"+ens0_r2m(st)+"_v"+z0+"_tend_"+mth+tstr+".txt"
  ; Total wind
    else
     diri_v = "$ar/text/ml_"+dat_r2m(st)+"_"+ens0_r2m(st)+"_vtot"+z0+"_tend_"+mth+tstr+".txt"
    end if 
   else
    diri_v  = "$ar/text/ml_"+dat_r2m(st)+"_"+ens0_r2m(st)+"_vtan_tend_"+mth+tstr+".txt"
   end if

 ; Mean sea level pressure (120 values; need to calculate tendency inline)                  
   diri_p   = "$ar/text/cp_new_"+dat_r2m(st)+"_"+ens0_r2m(st)+"_slp.txt"

 ; Radial vorticity profile (119 values)                                                    
   diri_rvp = "$ar/text/"+dat_r2m(st)+"_"+ens0_r2m(st)+"_sc"+sc+"_inner_core.txt"

 ; Radial vorticity profile tendency (119 values)                                           
   diri_ten = "$ar/text/"+dat_r2m(st)+"_"+ens0_r2m(st)+"_sc"+sc+"_inner_core_grad.txt"

 ; Mass flux (119 values)                                                                      
   diri_flux = "$ar/text/"+dat_r2m(st)+"_"+ens0_r2m(st)+"_sc"+sc+"_lay1_"+hy0+"_"+hy1+\
               "_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+".txt"

;==========================================================                                 
; Read in both wind speed tendency and MSLP (11/10/2019)                                    
;==========================================================                                 

  ; Wind speed tendency                                                                     
  ; Check whether file exists, and if so, read data from it                                 
    if (.not. fileexists(diri_v) ) then
     tend_all0 = new( ntimes_v, "float")
     print("Data missing for simulation "+ens0_r2m(st)+" ("+dat_r2m(st)+")" )
    else
     print("Working on RI period: "+st+" ("+ens0_r2m(st)+", "+dat_r2m(st)+")")
     tend_all0 = asciiread(diri_v, ntimes_v, "float")
     tend_all  = tend_all0(2::)
     delete(tend_all0)
    end if

  ; Mean sea level pressure                                                                 
  ; Check whether file containing MSLP data exists                                          
    if (.not. fileexists(diri_p) ) then
     mslp_all1 = new( ntimes_p, "float")
     print("Data missing for simulation "+ens0_r2m(st)+" ("+dat_r2m(st)+")" )
    else
     mslp_all1 = asciiread(diri_p, ntimes_p, "float")
     mslp_all0 = mslp_all1(1::)
     delete(mslp_all1)
    end if

  ; Define new array to hold MSLP values (119 values)                                       
    mslp_all   = new( ntimes_vc, "float")

  ; Calculate MSLP tendency inline                                                          
    do it = 0, ntimes_vc-1
     tp1 = min((/it+1,ntimes_vc-1/))
     tm1 = max((/it-1,0/))
     mslp_all(it) = ( mslp_all0(tp1) - mslp_all0(tm1) ) / 2
    end do

  ; Radial vorticity profile (and tendency)                                                 
    rvp_all     = asciiread(diri_rvp, ntimes_vc, "float")
    ten_all     = asciiread(diri_ten, ntimes_vc, "float")

  ; Finally, all data into correct arrays [wind speed + MSLP tendency; radial vorticity profile]
    npts = ( tf0_r2m(st) - ts0_r2m(st) )
    var_r2m(st,0:npts,0) = (/tend_all( ts0_r2m(st):tf0_r2m(st) )/)
    var_r2m(st,0:npts,1) = (/mslp_all( ts0_r2m(st):tf0_r2m(st) )/)
    var_r2m(st,0:npts,2) = (/rvp_all( ts0_r2m(st):tf0_r2m(st) )/)
    var_r2m(st,0:npts,3) = (/ten_all( ts0_r2m(st):tf0_r2m(st) )/)

    if (x .eq. "flux") then
     flux_all   = asciiread(diri_flux, ntimes_vc, "float")
     var_r2m(st,0:npts,4) = (/flux_all( ts0_r2m(st):tf0_r2m(st) )/)
     delete(diri_flux)
    end if

  end do
  delete([/diri_v, diri_p, diri_rvp, diri_ten/])

;==================================
; Start loop over monopole phase 
;==================================

  do st = 0, dsize(2)-1

 ; Read in data for each of our desired parameters                                         
 ; Tangential wind speed tendency (121 values)                                             
   if (zs .eq. "True") then
  ; Tangential wind
    if (wind .eq. "vtan") then 
     diri_v = "$ar/text/ml_"+dat_mono(st)+"_"+ens0_mono(st)+"_v"+z0+"_tend_"+mth+tstr+".txt"
  ; Total wind
    else
     diri_v = "$ar/text/ml_"+dat_mono(st)+"_"+ens0_mono(st)+"_vtot"+z0+"_tend_"+mth+tstr+".txt"
    end if 
   else
    diri_v  = "$ar/text/ml_"+dat_mono(st)+"_"+ens0_mono(st)+"_vtan_tend_"+mth+tstr+".txt"
   end if

 ; Mean sea level pressure (120 values; need to calculate tendency inline)                 
   diri_p   = "$ar/text/cp_new_"+dat_mono(st)+"_"+ens0_mono(st)+"_slp.txt"

 ; Radial vorticity profile (119 values)                                                   
   diri_rvp = "$ar/text/"+dat_mono(st)+"_"+ens0_mono(st)+"_sc"+sc+"_inner_core.txt"

 ; Radial vorticity profile tendency (119 values)                                          
   diri_ten = "$ar/text/"+dat_mono(st)+"_"+ens0_mono(st)+"_sc"+sc+"_inner_core_grad.txt"

 ; Mass flux (119 values)                                                                      
   diri_flux = "$ar/text/"+dat_mono(st)+"_"+ens0_mono(st)+"_sc"+sc+"_lay1_"+hy0+"_"+hy1+\
               "_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+".txt"

;==========================================================                                
; Read in both wind speed tendency and MSLP (11/10/2019)                                   
;==========================================================                                

  ; Wind speed tendency                                                                    
  ; Check whether file exists, and if so, read data from it                                
    if (.not. fileexists(diri_v) ) then
     tend_all0 = new( ntimes_v, "float")
     print("Data missing for simulation "+ens0_mono(st)+" ("+dat_mono(st)+")" )
    else
     print("Working on RI period: "+st+" ("+ens0_mono(st)+", "+dat_mono(st)+")")
     tend_all0 = asciiread(diri_v, ntimes_v, "float")
     tend_all  = tend_all0(2::)
     delete(tend_all0)
    end if

  ; Mean sea level pressure                                                                
  ; Check whether file containing MSLP data exists                                         
    if (.not. fileexists(diri_p) ) then
     mslp_all1 = new( ntimes_p, "float")
     print("Data missing for simulation "+ens0_mono(st)+" ("+dat_mono(st)+")" )
    else
     mslp_all1 = asciiread(diri_p, ntimes_p, "float")
     mslp_all0 = mslp_all1(1::)
     delete(mslp_all1)
    end if

  ; Define new array to hold MSLP values (119 values)                                      
    mslp_all   = new( ntimes_vc, "float")

  ; Calculate MSLP tendency inline                                                         
    do it = 0, ntimes_vc-1
     tp1 = min((/it+1,ntimes_vc-1/))
     tm1 = max((/it-1,0/))
     mslp_all(it) = ( mslp_all0(tp1) - mslp_all0(tm1) ) / 2
    end do

  ; Radial vorticity profile (and tendency)                                                
    rvp_all     = asciiread(diri_rvp, ntimes_vc, "float")
    ten_all     = asciiread(diri_ten, ntimes_vc, "float")

  ; Finally, all data into correct arrays [wind speed + MSLP tendency; radial vorticity profile]
    npts = ( tf0_mono(st) - ts0_mono(st) )
    var_mono(st,0:npts,0) = (/tend_all( ts0_mono(st):tf0_mono(st) )/)
    var_mono(st,0:npts,1) = (/mslp_all( ts0_mono(st):tf0_mono(st) )/)
    var_mono(st,0:npts,2) = (/rvp_all( ts0_mono(st):tf0_mono(st) )/)
    var_mono(st,0:npts,3) = (/ten_all( ts0_mono(st):tf0_mono(st) )/)

    if (x .eq. "flux") then
     flux_all   = asciiread(diri_flux, ntimes_vc, "float")
     var_mono(st,0:npts,4) = (/flux_all( ts0_mono(st):tf0_mono(st) )/)
     delete(diri_flux)
    end if

  end do
  delete([/diri_v, diri_p, diri_rvp, diri_ten/])

;===================================================                                     
; Start loop over monopole to ringlike transition                                        
;===================================================   

  do st = 0, dsize(3)-1

 ; Read in data for each of our desired parameters                                          
 ; Tangential wind speed tendency (121 values)                                              
   if (zs .eq. "True") then
  ; Tangential wind
    if (wind .eq. "vtan") then 
     diri_v = "$ar/text/ml_"+dat_m2r(st)+"_"+ens0_m2r(st)+"_v"+z0+"_tend_"+mth+tstr+".txt"
  ; Total wind
    else
     diri_v = "$ar/text/ml_"+dat_m2r(st)+"_"+ens0_m2r(st)+"_vtot"+z0+"_tend_"+mth+tstr+".txt"
    end if 
   else
    diri_v  = "$ar/text/ml_"+dat_m2r(st)+"_"+ens0_m2r(st)+"_vtan_tend_"+mth+tstr+".txt"
   end if

 ; Mean sea level pressure (120 values; need to calculate tendency inline)                  
   diri_p   = "$ar/text/cp_new_"+dat_m2r(st)+"_"+ens0_m2r(st)+"_slp.txt"

 ; Radial vorticity profile (119 values)                                                    
   diri_rvp = "$ar/text/"+dat_m2r(st)+"_"+ens0_m2r(st)+"_sc"+sc+"_inner_core.txt"

 ; Radial vorticity profile tendency (119 values)                                           
   diri_ten = "$ar/text/"+dat_m2r(st)+"_"+ens0_m2r(st)+"_sc"+sc+"_inner_core_grad.txt"

 ; Mass flux (119 values)                                                                      
   diri_flux = "$ar/text/"+dat_m2r(st)+"_"+ens0_m2r(st)+"_sc"+sc+"_lay1_"+hy0+"_"+hy1+\
               "_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+".txt"

;==========================================================                                 
; Read in both wind speed tendency and MSLP (11/10/2019)                                    
;==========================================================                                 

  ; Wind speed tendency                                                                     
  ; Check whether file exists, and if so, read data from it                                 
    if (.not. fileexists(diri_v) ) then
     tend_all0 = new( ntimes_v, "float")
     print("Data missing for simulation "+ens0_m2r(st)+" ("+dat_m2r(st)+")" )
    else
     print("Working on RI period: "+st+" ("+ens0_m2r(st)+", "+dat_m2r(st)+")")
     tend_all0 = asciiread(diri_v, ntimes_v, "float")
     tend_all  = tend_all0(2::)
     delete(tend_all0)
    end if

  ; Mean sea level pressure                                                                 
  ; Check whether file containing MSLP data exists                                          
    if (.not. fileexists(diri_p) ) then
     mslp_all1 = new( ntimes_p, "float")
     print("Data missing for simulation "+ens0_m2r(st)+" ("+dat_m2r(st)+")" )
    else
     mslp_all1 = asciiread(diri_p, ntimes_p, "float")
     mslp_all0 = mslp_all1(1::)
     delete(mslp_all1)
    end if

  ; Define new array to hold MSLP values (119 values)                                       
    mslp_all   = new( ntimes_vc, "float")

  ; Calculate MSLP tendency inline                                                          
    do it = 0, ntimes_vc-1
     tp1 = min((/it+1,ntimes_vc-1/))
     tm1 = max((/it-1,0/))
     mslp_all(it) = ( mslp_all0(tp1) - mslp_all0(tm1) ) / 2
    end do

  ; Radial vorticity profile (and tendency)                                                 
    rvp_all     = asciiread(diri_rvp, ntimes_vc, "float")
    ten_all     = asciiread(diri_ten, ntimes_vc, "float")

  ; Finally, all data into correct arrays [wind speed + MSLP tendency; radial vorticity profile]
    npts = ( tf0_m2r(st) - ts0_m2r(st) )
    var_m2r(st,0:npts,0) = (/tend_all( ts0_m2r(st):tf0_m2r(st) )/)
    var_m2r(st,0:npts,1) = (/mslp_all( ts0_m2r(st):tf0_m2r(st) )/)
    var_m2r(st,0:npts,2) = (/rvp_all( ts0_m2r(st):tf0_m2r(st) )/)
    var_m2r(st,0:npts,3) = (/ten_all( ts0_m2r(st):tf0_m2r(st) )/)

    if (x .eq. "flux") then
     flux_all   = asciiread(diri_flux, ntimes_vc, "float")
     var_m2r(st,0:npts,4) = (/flux_all( ts0_m2r(st):tf0_m2r(st) )/)
     delete(diri_flux)
    end if

  end do
  delete([/diri_v, diri_p, diri_rvp, diri_ten/])

;====================================================================                       
; Find indices of data points where wind speed exceeds a threshold                          
;====================================================================                       

 if (wrt .eq. 1) then 

  print("Finding indices where wind speed/MSLP exceeds a threshold...")

; Reshape 2D to 1D array (wind speed tendency)                                              
  ring_v_1d    = ndtooned(var_ring(:,:,0) )
  r2m_v_1d     = ndtooned(var_r2m(:,:,0) )
  mono_v_1d    = ndtooned(var_mono(:,:,0) )
  m2r_v_1d     = ndtooned(var_m2r(:,:,0) )

; Repeat for MSLP tendency 
  ring_slp_1d  = ndtooned(var_ring(:,:,1) )
  r2m_slp_1d   = ndtooned(var_r2m(:,:,1) )
  mono_slp_1d  = ndtooned(var_mono(:,:,1) )
  m2r_slp_1d   = ndtooned(var_m2r(:,:,1) )

; 29/10/2019 --> do I need to do the same for mass flux?

; Find indices of tendencies either side of user-defined threshold value
; 17/10/2019 --> either look at upper or lower end of distributions 

; Wind speed tendency 
  if (var .eq. "v") then
 ; Upper end of distribution (largest increases in wind speed)
   if (v0 .gt. 0) then 
    inds_ring_1d = ind(ring_v_1d .gt. v0)
    inds_r2m_1d  = ind(r2m_v_1d  .gt. v0)
    inds_mono_1d = ind(mono_v_1d .gt. v0)
    inds_m2r_1d  = ind(m2r_v_1d  .gt. v0)
 ; Lower end of distribution (largest decreases in wind speed)
   else
    inds_ring_1d = ind(ring_v_1d .lt. v0)
    inds_r2m_1d  = ind(r2m_v_1d  .lt. v0)
    inds_mono_1d = ind(mono_v_1d .lt. v0)
    inds_m2r_1d  = ind(m2r_v_1d  .lt. v0)
   end if 
; MSLP tendency
  else
 ; Lower end of distribution (largest decreases in MSLP)
   if (v0 .lt. 0) then 
    inds_ring_1d = ind(ring_slp_1d .lt. v0)
    inds_r2m_1d  = ind(r2m_slp_1d  .lt. v0)
    inds_mono_1d = ind(mono_slp_1d .lt. v0)
    inds_m2r_1d  = ind(m2r_slp_1d  .lt. v0)
 ; Upper end of distribution (largest increases in MSLP)
   else
    inds_ring_1d = ind(ring_slp_1d .gt. v0)
    inds_r2m_1d  = ind(r2m_slp_1d  .gt. v0)
    inds_mono_1d = ind(mono_slp_1d .gt. v0)
    inds_m2r_1d  = ind(m2r_slp_1d  .gt. v0)
   end if 
  end if 

; Resolve list of indices to their multi-dimensional representation                         
  inds_ring    = ind_resolve(inds_ring_1d, (/dsize(0), nts_ring/) )
  inds_r2m     = ind_resolve(inds_r2m_1d, (/dsize(1), nts_r2m/) )
  inds_mono    = ind_resolve(inds_mono_1d, (/dsize(2), nts_mono/) )
  inds_m2r     = ind_resolve(inds_m2r_1d, (/dsize(3), nts_m2r/) )

; For each of these sets of points, find the corresponding initialisation time ('dat')
; 16/10/2019 --> error occurred when running with v0=-3.0 ("slp"). Missing values? Account for this.
  dat_th_ring  = dat_ring(inds_ring(:,0) )
  dat_th_r2m   = dat_r2m(inds_r2m(:,0) )
  dat_th_mono  = dat_mono(inds_mono(:,0) )
  dat_th_m2r   = dat_m2r(inds_m2r(:,0) )  

; Do the same for simulation ('ens0')
  ens_th_ring  = ens0_ring(inds_ring(:,0) )
  ens_th_r2m   = ens0_r2m(inds_r2m(:,0) )
  ens_th_mono  = ens0_mono(inds_mono(:,0) )
  ens_th_m2r   = ens0_m2r(inds_m2r(:,0) )

; Finally, calculate the correct time indices using 'inds_xxxx(:,1)' --> index
; EDIT FROM HERE (16/10/2019) --> check that I've calculated these time indices correctly  
  time_th_ring = ts0_ring(inds_ring(:,0) ) + inds_ring(:,1)
  time_th_r2m  = ts0_r2m(inds_r2m(:,0) ) + inds_r2m(:,1)
  time_th_mono = ts0_mono(inds_mono(:,0) ) + inds_mono(:,1)
  time_th_m2r  = ts0_m2r(inds_m2r(:,0) ) + inds_m2r(:,1)  

; Size of each of these arrays
  thr_size     = (/dimsizes(time_th_ring), dimsizes(time_th_r2m), \
                   dimsizes(time_th_mono), dimsizes(time_th_m2r)  /)
  tot_size     = sum(thr_size)

; Extract start/end indices for each set of data points 
  t_s0         = thr_size(0)
  t_s1         = t_s0 + thr_size(1)
  t_s2         = t_s1 + thr_size(2)
  t_s3         = t_s2 + thr_size(3)

; Combine all values into single array to write out to text files
  dat_th_all   = new( (/tot_size/), "string")
  ens_th_all   = new( (/tot_size/), "string")
  time_th_all  = new( (/tot_size/), "integer")

; Read values from each VC phase into correct array elements  
  dat_th_all(0:t_s0-1)     = dat_th_ring
  dat_th_all(t_s0:t_s1-1)  = dat_th_r2m
  dat_th_all(t_s1:t_s2-1)  = dat_th_mono
  dat_th_all(t_s2:t_s3-1)  = dat_th_m2r

  ens_th_all(0:t_s0-1)     = ens_th_ring
  ens_th_all(t_s0:t_s1-1)  = ens_th_r2m
  ens_th_all(t_s1:t_s2-1)  = ens_th_mono
  ens_th_all(t_s2:t_s3-1)  = ens_th_m2r

  time_th_all(0:t_s0-1)    = time_th_ring
  time_th_all(t_s0:t_s1-1) = time_th_r2m
  time_th_all(t_s1:t_s2-1) = time_th_mono
  time_th_all(t_s2:t_s3-1) = time_th_m2r

; Write these values out to text files 
  if (zs .eq. "True") then 
  ; Tangential wind
   if (wind .eq. "vtan") then 
    dir_out = "$ar/text/vc_scatter_sc"+sc+"_"+var+"_"+v0+"_"+z0+"_"
  ; Total wind
   else
    dir_out = "$ar/text/vc_scatter_sc"+sc+"_"+var+wind+"_"+v0+"_"+z0+"_"
   end if 
  else
   dir_out  = "$ar/text/vc_scatter_sc"+sc+"_"+var+"_"+v0+"_"
  end if 

  dat_out   = dir_out+"dat.txt"
  ens_out   = dir_out+"ens.txt"
  time_out  = dir_out+"time.txt"

; Only if files don't already exist
  if (.not. fileexists(dat_out) ) then 
   asciiwrite(dat_out, dat_th_all)
  end if 
  if (.not. fileexists(ens_out) ) then 
   asciiwrite(ens_out, ens_th_all)
  end if 
  if (.not. fileexists(time_out) ) then 
   asciiwrite(time_out, time_th_all)
  end if 

 end if 

;===================
; Troubleshooting
;===================

 if (trb .eq. 1) then

  arr_ring = var_ring(:,:,0)
  arr_r2m  = var_r2m(:,:,0)
  arr_mono = var_mono(:,:,0)
  arr_m2r  = var_m2r(:,:,0)

; Size of all 2D arrays (how many pairs of points reach the threshold)
  msize0   = (/dimsizes(inds_ring), dimsizes(inds_r2m), dimsizes(inds_mono), dimsizes(inds_m2r)/)

; Print out values above/threshold in ringlike phase
  do nn = 0, msize0(0,0)-1
   arr_ind_ring = inds_ring(nn,:)
   print(arr_ring(arr_ind_ring(0), arr_ind_ring(1) ) )
  end do

; Print out values in ringlike to monopole transition
  do nn = 0, msize0(1,0)-1
   arr_ind_r2m  = inds_r2m(nn,:)
   print(arr_r2m(arr_ind_r2m(0), arr_ind_r2m(1) ) )
  end do

; Print out values from monopole phase  
  do nn = 0, msize0(2,0)-1
   arr_ind_mono = inds_mono(nn,:)
   print(arr_mono(arr_ind_mono(0), arr_ind_mono(1) ) )
  end do

; Print out values from monopole to ringlike transition
  do nn = 0, msize0(3,0)-1
   arr_ind_m2r  = inds_m2r(nn,:)
   print(arr_m2r(arr_ind_m2r(0), arr_ind_m2r(1) ) )
  end do

  exit()

 end if

;==============================================================================
; Reshape to form 1D arrays before plotting (29/10/2019 --> added mass flux) 
;==============================================================================
 
   if (res .eq. "1h") then
    fili_max  = (/"$ar/text/vtan_max_scatter_ring_group"+grp+"_sc"+sc+"_w"+w0+"_ten.txt",\
                  "$ar/text/vtan_max_scatter_r2m_group"+grp+"_sc"+sc+"_w"+w0+"_ten.txt",\
                  "$ar/text/vtan_max_scatter_mono_group"+grp+"_sc"+sc+"_w"+w0+"_ten.txt",\
                  "$ar/text/vtan_max_scatter_m2r_group"+grp+"_sc"+sc+"_w"+w0+"_ten.txt"/)
   else
    fili_max  = (/"$ar/text/vtan5_max_ring_group"+grp+"_sc"+sc+"_w"+w0+"_ten.txt",\
                  "$ar/text/vtan5_max_r2m_group"+grp+"_sc"+sc+"_w"+w0+"_ten.txt",\
                  "$ar/text/vtan5_max_mono_group"+grp+"_sc"+sc+"_w"+w0+"_ten.txt",\
                  "$ar/text/vtan5_max_m2r_group"+grp+"_sc"+sc+"_w"+w0+"_ten.txt"/)
   end if 
 
 ; Get size information for text files containing tendency of maximum tangential wind information
   xsize0    = new(4,"string")
   do il  = 0, 3
    xsize0(il) = systemfunc("wc -l < "+fili_max(il) )
   end do
   xsize  = toint(xsize0)

;===============================
; Ringlike phase (14/12/2019)
;===============================
   if (zs .eq. "max") then 
    ring_uv1= asciiread(fili_max(0),xsize(0),"float")
    ring_uv0= ndtooned(ring_uv1)
    
   ; Calculate number of non-missing values 
     ct = 0 
     do n = 0, xsize(0)-1
      if (.not. ismissing(ring_uv0(n) ) ) then 
       ct = ct + 1       
      end if 
     end do 
     print("Ring-like phase array contains "+ct+" values")
     xsize_out = ct
     ring_uv   = new(xsize_out,"float")

   ; Create array to hold only non-missing values 
     ct = 0 
     do n = 0, xsize(0)-1
      if (.not. ismissing(ring_uv0(n) ) ) then
       ring_uv(ct) = ring_uv0(n)
       ct = ct + 1 
      end if 
     end do 

 ; Otherwise, do as before and just read in values we already have 
   else

    ring_uv = ndtooned(var_ring(:,:,0) )

   end if 

   ring_slp = ndtooned(var_ring(:,:,1) )
   ring_rvp = ndtooned(var_ring(:,:,2) )
   ring_ten = ndtooned(var_ring(:,:,3) )
   ring_fl  = ndtooned(var_ring(:,:,4) ) 
   r_size   = dimsizes(ring_uv)

;================================================
; Ringlike to monopole transition (14/12/2019)
;================================================
   if (zs .eq. "max") then
    r2m_uv1= asciiread(fili_max(1),xsize(1),"float")
    r2m_uv0= ndtooned(r2m_uv1)

   ; Calculate number of non-missing values                                    
     ct = 0
     do n = 0, xsize(1)-1
      if (.not. ismissing(r2m_uv0(n) ) ) then
       ct = ct + 1
      end if
     end do
     print("Ring-like to monopole array contains "+ct+" values")
     xsize_out = ct
     r2m_uv   = new(xsize_out,"float")

   ; Create array to hold only non-missing values                              
     ct = 0
     do n = 0, xsize(1)-1
      if (.not. ismissing(r2m_uv0(n) ) ) then
       r2m_uv(ct) = r2m_uv0(n)
       ct = ct + 1
      end if
     end do

 ; Otherwise, do as before and just read in values we already have             
   else

    r2m_uv  = ndtooned(var_r2m(:,:,0) )

   end if 

   r2m_slp  = ndtooned(var_r2m(:,:,1) )
   r2m_rvp  = ndtooned(var_r2m(:,:,2) )
   r2m_ten  = ndtooned(var_r2m(:,:,3) )
   r2m_fl   = ndtooned(var_r2m(:,:,4) )
   r2m_size = dimsizes(r2m_uv)

;================================
; Monopole phase (14/12/2019)
;================================
   if (zs .eq. "max") then
    mono_uv1= asciiread(fili_max(2),xsize(2),"float")
    mono_uv0= ndtooned(mono_uv1)

   ; Calculate number of non-missing values                                      
     ct = 0
     do n = 0, xsize(2)-1
      if (.not. ismissing(mono_uv0(n) ) ) then
       ct = ct + 1
      end if
     end do
     print("Monopole phase array contains "+ct+" values")
     xsize_out = ct
     mono_uv   = new(xsize_out,"float")

   ; Create array to hold only non-missing values                                
     ct = 0
     do n = 0, xsize(2)-1
      if (.not. ismissing(mono_uv0(n) ) ) then
       mono_uv(ct) = mono_uv0(n)
       ct = ct + 1
      end if
     end do

 ; Otherwise, do as before and just read in values we already have               
   else

    mono_uv = ndtooned(var_mono(:,:,0) )

   end if

   mono_slp = ndtooned(var_mono(:,:,1) )
   mono_rvp = ndtooned(var_mono(:,:,2) )
   mono_ten = ndtooned(var_mono(:,:,3) )
   mono_fl  = ndtooned(var_mono(:,:,4) )
   m_size   = dimsizes(mono_uv)

;=================================================                                                 
; Monopole to ring-like transition (14/12/2019)                                                     
;=================================================                                
   if (zs .eq. "max") then
    m2r_uv1= asciiread(fili_max(3),xsize(3),"float")
    m2r_uv0= ndtooned(m2r_uv1)

   ; Calculate number of non-missing values                                      
     ct = 0
     do n = 0, xsize(3)-1
      if (.not. ismissing(m2r_uv0(n) ) ) then
       ct = ct + 1
      end if
     end do
     print("Monopole to ring-like transition array contains "+ct+" values")
     xsize_out = ct
     m2r_uv   = new(xsize_out,"float")

   ; Create array to hold only non-missing values                                
     ct = 0
     do n = 0, xsize(3)-1
      if (.not. ismissing(m2r_uv0(n) ) ) then
       m2r_uv(ct) = m2r_uv0(n)
       ct = ct + 1
      end if
     end do

 ; Otherwise, do as before and just read in values we already have               
   else

    m2r_uv  = ndtooned(var_m2r(:,:,0) )

   end if

   m2r_slp  = ndtooned(var_m2r(:,:,1) )
   m2r_rvp  = ndtooned(var_m2r(:,:,2) )
   m2r_ten  = ndtooned(var_m2r(:,:,3) )
   m2r_fl   = ndtooned(var_m2r(:,:,4) )
   m2r_size = dimsizes(m2r_uv)

;==============================================
; Find and remove missing values (if needed)
;==============================================

 if (miss .eq. 1) then

 ; File paths for each set of options 
   if (pm .eq. 1) then
    diri_msg = (/"$ar/text/dvdt_rmw_±"+thr+"_sc"+sc+"_group"+grp+"_w"+w0+"_0.005km_ring.txt",\
                 "$ar/text/dvdt_rmw_±"+thr+"_sc"+sc+"_group"+grp+"_w"+w0+"_0.005km_r2m.txt" ,\
                 "$ar/text/dvdt_rmw_±"+thr+"_sc"+sc+"_group"+grp+"_w"+w0+"_0.005km_mono.txt",\
                 "$ar/text/dvdt_rmw_±"+thr+"_sc"+sc+"_group"+grp+"_w"+w0+"_0.005km_m2r.txt"/)
   elseif (pm .eq. 0) then
    diri_msg = (/"$ar/text/dvdt_rmw_"+thr+"_sc"+sc+"_group"+grp+"_w"+w0+"_0.005km_ring.txt",\
                 "$ar/text/dvdt_rmw_"+thr+"_sc"+sc+"_group"+grp+"_w"+w0+"_0.005km_r2m.txt" ,\
                 "$ar/text/dvdt_rmw_"+thr+"_sc"+sc+"_group"+grp+"_w"+w0+"_0.005km_mono.txt",\
                 "$ar/text/dvdt_rmw_"+thr+"_sc"+sc+"_group"+grp+"_w"+w0+"_0.005km_m2r.txt"/)
   elseif (pm .eq. -1) then
    diri_msg = (/"$ar/text/dvdt_r"+rstr_out+"_sc"+sc+"_group"+grp+"_w"+w0+"_0.005km_ring.txt",\
                 "$ar/text/dvdt_r"+rstr_out+"_sc"+sc+"_group"+grp+"_w"+w0+"_0.005km_r2m.txt" ,\
                 "$ar/text/dvdt_r"+rstr_out+"_sc"+sc+"_group"+grp+"_w"+w0+"_0.005km_mono.txt",\
                 "$ar/text/dvdt_r"+rstr_out+"_sc"+sc+"_group"+grp+"_w"+w0+"_0.005km_m2r.txt"/)
   elseif (pm .eq. -2) then
    diri_msg = (/"$ar/text/dvdt_r"+rstr_out+"_sc"+sc+"_group"+grp+"_w"+w0+"_0.005km_ring.txt",\
                 "$ar/text/dvdt_r"+rstr_out+"_sc"+sc+"_group"+grp+"_w"+w0+"_0.005km_r2m.txt" ,\
                 "$ar/text/dvdt_r"+rstr_out+"_sc"+sc+"_group"+grp+"_w"+w0+"_0.005km_mono.txt",\
                 "$ar/text/dvdt_r"+rstr_out+"_sc"+sc+"_group"+grp+"_w"+w0+"_0.005km_m2r.txt"/)
   elseif (pm .eq. -3) then
    diri_msg = (/"$ar/text/dvdt_z"+zthr+"km_sc"+sc+"_group"+grp+"_w"+w0+"_0.005km_ring.txt",\
                 "$ar/text/dvdt_z"+zthr+"km_sc"+sc+"_group"+grp+"_w"+w0+"_0.005km_r2m.txt" ,\
                 "$ar/text/dvdt_z"+zthr+"km_sc"+sc+"_group"+grp+"_w"+w0+"_0.005km_mono.txt",\
                 "$ar/text/dvdt_z"+zthr+"km_sc"+sc+"_group"+grp+"_w"+w0+"_0.005km_m2r.txt"/)
   elseif (pm .eq. -4) then
    diri_msg = (/"$ar/text/dvdt_z"+zthr+"km_r20km_sc"+sc+"_group"+grp+"_w"+w0+"_0.005km_ring.txt",\
                 "$ar/text/dvdt_z"+zthr+"km_r20km_sc"+sc+"_group"+grp+"_w"+w0+"_0.005km_r2m.txt" ,\
                 "$ar/text/dvdt_z"+zthr+"km_r20km_sc"+sc+"_group"+grp+"_w"+w0+"_0.005km_mono.txt",\
                 "$ar/text/dvdt_z"+zthr+"km_r20km_sc"+sc+"_group"+grp+"_w"+w0+"_0.005km_m2r.txt"/)
   elseif (pm .eq. -5) then
    diri_msg = (/"$ar/text/dvdt_z"+zthr+"km_r"+rstr_out+"_sc"+sc+"_group"+grp+\
                 "_w"+w0+"_0.005km_ring.txt",\
                 "$ar/text/dvdt_z"+zthr+"km_r"+rstr_out+"_sc"+sc+"_group"+grp+\
                 "_w"+w0+"_0.005km_r2m.txt" ,\
                 "$ar/text/dvdt_z"+zthr+"km_r"+rstr_out+"_sc"+sc+"_group"+grp+\
                 "_w"+w0+"_0.005km_mono.txt",\
                 "$ar/text/dvdt_z"+zthr+"km_r"+rstr_out+"_sc"+sc+"_group"+grp+\
                 "_w"+w0+"_0.005km_m2r.txt"/)
   end if

 ; Create array to hold indices of the data points we're removing in each instance 
   size_msg0 = new(4,"string")
   do il  = 0, 3
    size_msg0(il) = systemfunc("wc -l < "+diri_msg(il) )
   end do
   size_msg  = toint(size_msg0)

 ; Create new arrays with reduced number of data points, to hold no missing values 
   ring_nms = new( (/5, sum_ring/), "float")
   r2m_nms  = new( (/5, sum_r2m/), "float")
   mono_nms = new( (/5, sum_mono/), "float")
   m2r_nms  = new( (/5, sum_m2r/), "float")

;===================================================
; First, create new arrays with no missing values
;===================================================

   print("Removing missing values...")

 ; Ring-like phase 
   ct = 0 
   do n = 0, r_size-1
    if (.not. ismissing(ring_uv(n) ) ) then 
     ring_nms(0,ct) = ring_uv(n)
     ring_nms(1,ct) = ring_slp(n)
     ring_nms(2,ct) = ring_rvp(n)
     ring_nms(3,ct) = ring_ten(n)
     ring_nms(4,ct) = ring_fl(n)
     ct = ct + 1
    end if 
   end do 

 ; Ring-like to monopole transition 
   ct = 0
   do n = 0, r2m_size-1
    if (.not.ismissing(r2m_uv(n) ) ) then
     r2m_nms(0,ct)  = r2m_uv(n)
     r2m_nms(1,ct)  = r2m_slp(n)
     r2m_nms(2,ct)  = r2m_rvp(n)
     r2m_nms(3,ct)  = r2m_ten(n)
     r2m_nms(4,ct)  = r2m_fl(n)
     ct = ct + 1
    end if
   end do

 ; Monopole phase 
   ct = 0 
   do n = 0, m_size-1
    if (.not. ismissing(mono_uv(n) ) ) then
     mono_nms(0,ct) = mono_uv(n)
     mono_nms(1,ct) = mono_slp(n)
     mono_nms(2,ct) = mono_rvp(n)
     mono_nms(3,ct) = mono_ten(n)
     mono_nms(4,ct) = mono_fl(n)
     ct = ct + 1
    end if
   end do

 ; Monopole to ring-like transition 
   ct = 0
   do n = 0, m2r_size-1
    if (.not. ismissing(m2r_uv(n) ) ) then
     m2r_nms(0,ct)  = m2r_uv(n)
     m2r_nms(1,ct)  = m2r_slp(n)
     m2r_nms(2,ct)  = m2r_rvp(n)
     m2r_nms(3,ct)  = m2r_ten(n)
     m2r_nms(4,ct)  = m2r_fl(n)
     ct = ct + 1
    end if
   end do

;======================================================
; Now remove selected values based on r/z thresholds
;======================================================

   print("Now removing values based on r/z thresholds...")

 ; Read in indices of values to remove, for each VC phase                                         
   ring_msg = asciiread(diri_msg(0), size_msg(0), "integer")
   r2m_msg  = asciiread(diri_msg(1), size_msg(1), "integer")
   mono_msg = asciiread(diri_msg(2), size_msg(2), "integer")
   m2r_msg  = asciiread(diri_msg(3), size_msg(3), "integer")

 ; Create another set of arrays with the selected dv/dt values also removed
   ring_fin = new( (/5, sum_ring - size_msg(0)/), "float")
   r2m_fin  = new( (/5, sum_r2m - size_msg(1)/), "float")
   mono_fin = new( (/5, sum_mono - size_msg(2)/), "float")
   m2r_fin  = new( (/5, sum_m2r - size_msg(3)/), "float")

 ; Set elements where dv/dt occurred within the RMW to missing 
   do n = 0, dimsizes(ring_msg)-1
    ring_nms(:,ring_msg(n) ) = ring_nms@_FillValue
   end do 

   do n = 0, dimsizes(r2m_msg)-1
    r2m_nms(:,r2m_msg(n) ) = r2m_nms@_FillValue
   end do

   do n = 0, dimsizes(mono_msg)-1
    mono_nms(:,mono_msg(n) ) = mono_nms@_FillValue
   end do

   do n = 0, dimsizes(m2r_msg)-1
    m2r_nms(:,m2r_msg(n) ) = m2r_nms@_FillValue
   end do

 ; Now finally create arrays with values removed (12/11/2019 --> quicker way to do this?)

 ; Ring-like phase
   ct = 0
   do n = 0, sum_ring-1
    if (.not. ismissing(ring_nms(0,n) ) ) then 
     ring_fin(0,ct) = ring_nms(0,n)
     ring_fin(1,ct) = ring_nms(1,n)
     ring_fin(2,ct) = ring_nms(2,n)
     ring_fin(3,ct) = ring_nms(3,n)
     ring_fin(4,ct) = ring_nms(4,n)
     ct = ct + 1
    end if 
   end do 
   ringsize_f = dimsizes(ring_fin(0,:) ) 

 ; Ring-like to monopole transition
   ct = 0
   do n = 0, sum_r2m-1
    if (.not. ismissing(r2m_nms(0,n) ) ) then
     r2m_fin(0,ct) = r2m_nms(0,n)
     r2m_fin(1,ct) = r2m_nms(1,n)
     r2m_fin(2,ct) = r2m_nms(2,n)
     r2m_fin(3,ct) = r2m_nms(3,n)
     r2m_fin(4,ct) = r2m_nms(4,n)
     ct = ct + 1
    end if
   end do
   r2msize_f = dimsizes(r2m_fin(0,:) )

 ; Monopole phase 
   ct = 0
   do n = 0, sum_mono-1
    if (.not. ismissing(mono_nms(0,n) ) ) then
     mono_fin(0,ct) = mono_nms(0,n)
     mono_fin(1,ct) = mono_nms(1,n)
     mono_fin(2,ct) = mono_nms(2,n)
     mono_fin(3,ct) = mono_nms(3,n)
     mono_fin(4,ct) = mono_nms(4,n)
     ct = ct + 1
    end if
   end do
   monosize_f = dimsizes(mono_fin(0,:) )

 ; Monopole to ring-like transition
   ct = 0
   do n = 0, sum_m2r-1
    if (.not. ismissing(m2r_nms(0,n) ) ) then
     m2r_fin(0,ct) = m2r_nms(0,n)
     m2r_fin(1,ct) = m2r_nms(1,n)
     m2r_fin(2,ct) = m2r_nms(2,n)
     m2r_fin(3,ct) = m2r_nms(3,n)
     m2r_fin(4,ct) = m2r_nms(4,n)
     ct = ct + 1
    end if
   end do
   m2rsize_f = dimsizes(m2r_fin(0,:) )

 end if ;  End IF statement [if (miss .eq. 1) ] --> 

 ; Size information needed to create single array to hold all values 
 ; 12/11/2019 --> depending on whether we ignore intensification inside the RMW 
   if (miss .eq. 1) then 
    all_size = ringsize_f + r2msize_f + monosize_f + m2rsize_f
    s0       = ringsize_f
    s1       = s0 + r2msize_f
    s2       = s1 + monosize_f
    s3       = s2 + m2rsize_f
   else 
    all_size = r_size + r2m_size + m_size + m2r_size  
    s0       = r_size
    s1       = s0 + r2m_size
    s2       = s1 + m_size
    s3       = s2 + m2r_size
   end if 

 ; If reading in our new wind tendency data, remove all the missing values before plotting 


 ; Create new arrays to hold all values for each variable 
   all_uv    = new( (/all_size/), "float")
   all_slp   = new( (/all_size/), "float")
   all_rvp   = new( (/all_size/), "float")
   all_ten   = new( (/all_size/), "float")
   all_flux  = new( (/all_size/), "float")

 ; Fill each of these arrays with the values from each VC phase 
   if (miss .eq. 1) then 
    all_uv(0:s0-1)    = ring_fin(0,:)
    all_uv(s0:s1-1)   = r2m_fin(0,:)
    all_uv(s1:s2-1)   = mono_fin(0,:)
    all_uv(s2:s3-1)   = m2r_fin(0,:)

    all_slp(0:s0-1)   = ring_fin(1,:)
    all_slp(s0:s1-1)  = r2m_fin(1,:)
    all_slp(s1:s2-1)  = mono_fin(1,:)
    all_slp(s2:s3-1)  = m2r_fin(1,:)

    all_rvp(0:s0-1)   = ring_fin(2,:)
    all_rvp(s0:s1-1)  = r2m_fin(2,:)
    all_rvp(s1:s2-1)  = mono_fin(2,:)
    all_rvp(s2:s3-1)  = m2r_fin(2,:)

    all_ten(0:s0-1)   = ring_fin(3,:)
    all_ten(s0:s1-1)  = r2m_fin(3,:)
    all_ten(s1:s2-1)  = mono_fin(3,:)
    all_ten(s2:s3-1)  = m2r_fin(3,:)

    all_flux(0:s0-1)  = ring_fin(4,:)
    all_flux(s0:s1-1) = r2m_fin(4,:)
    all_flux(s1:s2-1) = mono_fin(4,:)
    all_flux(s2:s3-1) = m2r_fin(4,:)
   else
    all_uv(0:s0-1)    = ring_uv
    all_uv(s0:s1-1)   = r2m_uv
    all_uv(s1:s2-1)   = mono_uv
    all_uv(s2:s3-1)   = m2r_uv

    all_slp(0:s0-1)   = ring_slp
    all_slp(s0:s1-1)  = r2m_slp
    all_slp(s1:s2-1)  = mono_slp
    all_slp(s2:s3-1)  = m2r_slp

    all_rvp(0:s0-1)   = ring_rvp
    all_rvp(s0:s1-1)  = r2m_rvp
    all_rvp(s1:s2-1)  = mono_rvp
    all_rvp(s2:s3-1)  = m2r_rvp

    all_ten(0:s0-1)   = ring_ten
    all_ten(s0:s1-1)  = r2m_ten
    all_ten(s1:s2-1)  = mono_ten
    all_ten(s2:s3-1)  = m2r_ten

    all_flux(0:s0-1)  = ring_fl
    all_flux(s0:s1-1) = r2m_fl
    all_flux(s1:s2-1) = mono_fl
    all_flux(s2:s3-1) = m2r_fl
   end if 

 ; Add metadata so we can colour-code the markers by VC phase
   all_uv!0  = "vc_phase"
   all_slp!0 = "vc_phase"
   all_rvp!0 = "vc_phase"
   all_ten!0 = "vc_phase"
   all_flux!0= "vc_phase"

 ; Change units before plotting 
   all_rvp        = all_rvp * 10
   all_rvp@units  = "10~S~-2~N~ h~S~-1~N~"
   all_ten        = all_ten * 10
   all_ten@units  = "10~S~-2~N~ h~S~-1~N~"   ; EDIT FROM HERE --> check units and replot
   all_flux@units = "kg m~S~-2~N~ s~S~-1~N~"

;===========================================================
; If desired, remove selected data points before plotting 
;===========================================================

 ; I.e. ignore all instances of weakening tangential wind, or increasing MSLP 
   if (sub .eq. 1) then 
    if (var .eq. "v") then 
     all_uv = where(all_uv .lt. t0, all_uv@_FillValue, all_uv)
    elseif (var .eq. "slp") then 
     all_slp = where(all_slp .gt. t0, all_slp@_FillValue, all_slp)
    end if 
   end if    

;================================================
; Open workstation and define output file path
;================================================

; Define output string detailing which values we have ignored
  if (pm .eq. -1) then
   pstr = "r"+rthr
  elseif (pm .eq. -2) then 
   pstr = rstr_out
  elseif (pm .eq. -3) then
   pstr = "z"+zthr
  elseif (pm .eq. -4) then
   pstr = "r20_z"+zthr
  elseif (pm .eq. -5) then 
   pstr = rstr_out+"_"+zthr
  else
   pstr = ""
  end if

; Define additional output string to tidy up output file names 
  if (miss .eq. 1) then 
   mstr = "miss1"
  else
   mstr = ""
  end if 

  if (sub .eq. 1) then 
   sstr = "_sub"+t0+"_"
  else
   sstr = ""
  end if 

 ; Calculate maximum tendency using data on single model level
   if (zs .eq. "True") then 

  ; Tangential wind
    if (wind .eq. "vtan") then 
     output = "$nep/nepartak/images/vc_scatter_group"+grp+"_sc"+sc+"_"+var\
              +sstr+z0+"_"+mstr+"_"+thr+"_"+pstr

  ; Total wind
    else
     output = "$nep/nepartak/images/vc_scatter_group"+grp+"_sc"+sc+\
              "_"+var+wind+sstr+z0+"_"+mstr+"_"+thr+"_"+pstr
    end if 

   elseif (zs .eq. "max") then 
     output = "$nep/nepartak/images/vc_scatter_group"+grp+"_sc"+sc+\
              "_vtan_max"

 ; OR calculate maximum tendency using data on all levels
   else
    output  = "$nep/nepartak/images/vc_scatter_group"+grp+"_sc"+sc+\
              "_"+var+sstr+"_"+mstr+"_"+thr+"_"+pstr
   end if 
   wks      = gsn_open_wks(opt, output)
   print(output)

;==================================================
; Now produce scatter plot using all data points 
;==================================================

 ; Plotting resources 
   vc_opts                       = True
   vc_opts@gsnDraw               = False
   vc_opts@gsnFrame              = False
   vc_opts@gsnMaximize           = False

 ; Y-axis                                                                                 
   vc_opts@tiYAxisString         = "Wind speed tendency (m s~S~-1~N~ h~S~-1~N~)"
   vc_opts@tiYAxisFontHeightF    = 0.0150

   if (sub .ne. 1) then  
    vc_opts@trYMaxF              = 5.0
    vc_opts@trYMinF              = -5.0
   else
    vc_opts@trYMaxF              = 5.0
    vc_opts@trYMinF              = t0
   end if 

 ; X-axis (define our own labels; set limits)
   if (x .eq. "flux") then 
    vc_opts@tiXAxisString        = "Difference in vertical mass flux ("+all_flux@units+")"
    vc_opts@trXMinF              = -3.0
    vc_opts@trXMaxF              = 2.0
   elseif (x .eq. "rvp") then 
    vc_opts@tiXAxisString        = "Radial vorticity profile ("+all_rvp@units+")"
    vc_opts@trXMinF              = -10.0
    vc_opts@trXMaxF              = 10.0
   else
    vc_opts@tiXAxisString        = "Tendency of radial vorticity profile ("+all_ten@units+")"
    vc_opts@trXMinF              = -1.5
    vc_opts@trXMaxF              = 1.5
   end if 
   vc_opts@tiXAxisFontHeightF    = 0.0150

 ; Additional resources                                                              
   vc_opts@xyMarkLineMode        = "Markers"
   vc_opts@xyMonoMarkerColor     = True
   vc_opts@xyMarkerColor         = "RoyalBlue1"
   vc_opts@xyMonoMarkerSize      = True
   vc_opts@xyMonoMarkerThickness = True

 ; Marker size and thickness 
   vc_opts@xyMarkerSizeF         = 0.005
   vc_opts@xyMarkerThicknessF    = 0.5
   vc_opts@xyMonoMarker          = True 
   vc_opts@xyMarker              = 16

 ; Turn off selected axes and tickmarks
   vc_opts@tmXTOn                = False
   vc_opts@tmYROn                = False

 ; Panel position and size 
   if (var .eq. "both") then 
    vc_opts@vpXF                 = 0.09
    vc_opts@vpYF                 = 0.70
    vc_opts@vpWidthF             = 0.35
    vc_opts@vpHeightF            = 0.32
   else
    vc_opts@vpXF                 = 0.12
    vc_opts@vpYF                 = 0.88
    vc_opts@vpWidthF             = 0.75
    vc_opts@vpHeightF            = 0.75
   end if 

 ; Additional marker resources 
   vc_opts1                     = True 
   vc_opts1@gsMarkerIndex       = 16
   vc_opts1@gsMarkerColor       = "purple4"
   vc_opts1@gsMarkerSizeF       = 0.005
   vc_opts1@gsMarkerThicknessF  = 0.5

   vc_opts2                     = vc_opts1
   vc_opts2@gsMarkerColor       = "red3"
   
   vc_opts3                     = vc_opts1
   vc_opts3@gsMarkerColor       = "orange3"

 ; Plot first panel (tangential wind speed tendency vs radial vorticity profile)
   if (x .eq. "flux") then 
    vc_plot_uv0 = gsn_csm_xy(wks, all_flux(0:s0-1), all_uv(0:s0-1), vc_opts)
    uv1         = gsn_add_polymarker(wks, vc_plot_uv0, all_flux(s0:s1-1), all_uv(s0:s1-1), vc_opts1)
    uv2         = gsn_add_polymarker(wks, vc_plot_uv0, all_flux(s1:s2-1), all_uv(s1:s2-1), vc_opts2)
    uv3         = gsn_add_polymarker(wks, vc_plot_uv0, all_flux(s2:s3-1), all_uv(s2:s3-1), vc_opts3)
   elseif (x .eq. "rvp") then 
    vc_plot_uv0 = gsn_csm_xy(wks, all_rvp(0:s0-1), all_uv(0:s0-1), vc_opts)
    uv1         = gsn_add_polymarker(wks, vc_plot_uv0, all_rvp(s0:s1-1), all_uv(s0:s1-1), vc_opts1)
    uv2         = gsn_add_polymarker(wks, vc_plot_uv0, all_rvp(s1:s2-1), all_uv(s1:s2-1), vc_opts2)
    uv3         = gsn_add_polymarker(wks, vc_plot_uv0, all_rvp(s2:s3-1), all_uv(s2:s3-1), vc_opts3)
   else
    vc_plot_uv0 = gsn_csm_xy(wks, all_ten(0:s0-1), all_uv(0:s0-1), vc_opts)
    uv1         = gsn_add_polymarker(wks, vc_plot_uv0, all_ten(s0:s1-1), all_uv(s0:s1-1), vc_opts1)
    uv2         = gsn_add_polymarker(wks, vc_plot_uv0, all_ten(s1:s2-1), all_uv(s1:s2-1), vc_opts2)
    uv3         = gsn_add_polymarker(wks, vc_plot_uv0, all_ten(s2:s3-1), all_uv(s2:s3-1), vc_opts3)
   end if 

 ; Resources for second panel  
   vc_opts4                     = vc_opts
   vc_opts4@tiYAxisString       = "MSLP tendency (hPa h~S~-1~N~)"

   if (var .eq. "both") then 
    vc_opts4@vpXF               = 0.57
    vc_opts4@vpYF               = 0.70
    vc_opts4@vpWidthF           = 0.35
    vc_opts4@vpHeightF          = 0.32
   else
    vc_opts4@vpXF               = 0.12
    vc_opts4@vpYF               = 0.88
    vc_opts4@vpWidthF           = 0.75
    vc_opts4@vpHeightF          = 0.75
   end if 

   if (sub .ne. 1) then 
    vc_opts4@trYMaxF            = 6.0
    vc_opts4@trYMinF            = -6.0
   else
    vc_opts4@trYMaxF            = t0
    vc_opts4@trYMinF            = -6.0
   end if 

 ; Plot second panel (MSLP tendency vs radial vorticity profile)
   if (x .eq. "flux") then 
    vc_plot_slp0 = gsn_csm_xy(wks, all_flux(0:s0-1), all_slp(0:s0-1), vc_opts4)
    slp1         = gsn_add_polymarker(wks,vc_plot_slp0,all_flux(s0:s1-1),all_slp(s0:s1-1), vc_opts1)
    slp2         = gsn_add_polymarker(wks,vc_plot_slp0,all_flux(s1:s2-1),all_slp(s1:s2-1), vc_opts2)
    slp3         = gsn_add_polymarker(wks,vc_plot_slp0,all_flux(s2:s3-1),all_slp(s2:s3-1), vc_opts3)
   elseif (x .eq. "rvp") then 
    vc_plot_slp0 = gsn_csm_xy(wks, all_rvp(0:s0-1), all_slp(0:s0-1), vc_opts4)
    slp1         = gsn_add_polymarker(wks,vc_plot_slp0,all_rvp(s0:s1-1),all_slp(s0:s1-1), vc_opts1)
    slp2         = gsn_add_polymarker(wks,vc_plot_slp0,all_rvp(s1:s2-1),all_slp(s1:s2-1), vc_opts2)
    slp3         = gsn_add_polymarker(wks,vc_plot_slp0,all_rvp(s2:s3-1),all_slp(s2:s3-1), vc_opts3)
   else
    vc_plot_slp0 = gsn_csm_xy(wks, all_ten(0:s0-1), all_slp(0:s0-1), vc_opts4)
    slp1         = gsn_add_polymarker(wks,vc_plot_slp0,all_ten(s0:s1-1),all_slp(s0:s1-1), vc_opts1)
    slp2         = gsn_add_polymarker(wks,vc_plot_slp0,all_ten(s1:s2-1),all_slp(s1:s2-1), vc_opts2)
    slp3         = gsn_add_polymarker(wks,vc_plot_slp0,all_ten(s2:s3-1),all_slp(s2:s3-1), vc_opts3)
   end if

;===============================
; Before plotting, add legend
;===============================

 ; Display legend, control font and label height                                        
   lg_opts                            = True
   lg_opts@pmLegendDisplayMode        = "Always"
   lg_opts@lgAutoManage               = False
   lg_opts@lgLabelFont                = "Helvetica"
   lg_opts@lgLabelFontHeightF         = 0.050
    
 ; Choose item type (markers or lines)
   lg_opts@lgItemType                 = "Markers"

 ; Labels and marker colours
   marker_colours                     = (/vc_opts@xyMarkerColor, vc_opts1@gsMarkerColor, \ 
                                          vc_opts2@gsMarkerColor, vc_opts3@gsMarkerColor/)
   lg_opts@lgMarkerColors             = marker_colours
   lg_opts@lgMarkerIndexes            = (/16, 16, 16, 16/) 
   lg_opts@lgMarkerSizeF              = vc_opts@xyMarkerSizeF * 2
   lg_opts@lgMarkerThicknessF         = vc_opts@xyMarkerThicknessF * 2
   msize                              = dimsizes(marker_colours)

 ; Label position within box                                                            
   lg_opts@lgLabelPosition            = "Right"
   lg_opts@lgItemPlacement            = "ExplicitPlacement"
   lg_opts@lgItemPositions            = fspan(0.10, 0.90, msize)
   lg_opts@lgItemOrder                = (/3, 2, 1, 0/)

 ; Legend box size                                                                      
   lg_opts@vpWidthF                   = 0.100
   lg_opts@vpHeightF                  = 0.0650

   lg_opts@lgPerimColor               = "black"
   lg_opts@lgPerimThicknessF          = 4.0
   lg_opts@lgPerimFill                = "SolidFill"
   lg_opts@lgPerimFillColor           = "white"

 ; Customise legend labels depending on number and choice of input variables
   lab_arr                            = (/"ring", "r2m", "mono", "m2r"/)
   lgd                                = gsn_create_legend(wks, msize, lab_arr, lg_opts)

 ; Add legend to plot (top left or right hand corner)                                   
   am_opts                            = True

   if (lgd0 .eq. "tr") then
    am_opts@amJust                    = "TopRight"
    am_opts@amParallelPosF            = 0.5
    am_opts@amOrthogonalPosF          = -0.5
   elseif (lgd0 .eq. "tl") then
    am_opts@amJust                    = "TopLeft"
    am_opts@amParallelPosF            = -0.5
    am_opts@amOrthogonalPosF          = -0.5
   elseif (lgd0 .eq. "br") then
    am_opts@amJust                    = "BottomRight"
    am_opts@amParallelPosF            = 0.5
    am_opts@amOrthogonalPosF          = 0.5
   elseif (lgd0 .eq. "bl") then
    am_opts@amJust                    = "BottomLeft"
    am_opts@amParallelPosF            = -0.5
    am_opts@amOrthogonalPosF          = 0.5
   end if

 ; Add legend to selected panels 
 ; 13/10/2019 --> for multiple panels, do you have to define separate resources?
   if (var .eq. "v") then 
    annotate_left                     = gsn_add_annotation(vc_plot_uv0, lgd, am_opts)
   elseif (var .eq. "slp") then
    annotate_right                    = gsn_add_annotation(vc_plot_slp0, lgd, am_opts)
   else
    annotate_left                     = gsn_add_annotation(vc_plot_uv0, lgd, am_opts)
    annotate_right                    = gsn_add_annotation(vc_plot_slp0, lgd, am_opts)
   end if 

 ; Either draw the scatter plot, or skip this part of the script 
   if (sct .eq. 1) then 

;===============================================
; Now switch to Cartesian grid in both panels
;===============================================

    if (var .eq. "v") then 
     c_axis_uv  = cartesian_axis(wks, vc_plot_uv0)
    elseif (var .eq. "slp") then 
     c_axis_slp = cartesian_axis(wks, vc_plot_slp0)
    else
     c_axis_uv  = cartesian_axis(wks, vc_plot_uv0)
     c_axis_slp = cartesian_axis(wks, vc_plot_slp0)
    end if 

;===================================
; Finally, draw selected panel(s)
;===================================

    if (var .eq. "v") then 
     draw(vc_plot_uv0)
    elseif (var .eq. "slp") then 
     draw(vc_plot_slp0)
    else
     draw(vc_plot_uv0)
     draw(vc_plot_slp0)
    end if 

    frame(wks)

 ; End IF statement (if sct .eq. 1) 
   end if 

;===================================                                                      
; Also draw box plot, if desired                                                          
;=================================== 

; Create array to hold all box plot values                                                  
; 4 VC phases                                                                               
; 5 values per phase [minimum, lower quartile, median, upper quartile, maximum]             
  box_y    = new( (/4,5/), "float")

; Options for 'stat_dispersion' function                                                    
  opts     = True
  opts@PrintStat = False

; Change units of radial vorticity profile back before plotting here 
  all_rvp   = all_rvp / 100
  all_rvp@units  = "h~S~-1~N~"

; Calculate stats for each VC phase                                                          
  if (boxv .eq. "v") then 
   ptb_ring = stat_dispersion(all_uv(0:s0-1), opts)
   ptb_r2m  = stat_dispersion(all_uv(s0:s1-1), opts)
   ptb_mono = stat_dispersion(all_uv(s1:s2-1), opts)
   ptb_m2r  = stat_dispersion(all_uv(s2:s3-1), opts)
  elseif (boxv .eq. "slp") then 
   ptb_ring = stat_dispersion(all_slp(0:s0-1), opts)
   ptb_r2m  = stat_dispersion(all_slp(s0:s1-1), opts)
   ptb_mono = stat_dispersion(all_slp(s1:s2-1), opts)
   ptb_m2r  = stat_dispersion(all_slp(s2:s3-1), opts)
  elseif (boxv .eq. "rvp") then 
   ptb_ring = stat_dispersion(all_rvp(0:s0-1), opts)
   ptb_r2m  = stat_dispersion(all_rvp(s0:s1-1), opts)
   ptb_mono = stat_dispersion(all_rvp(s1:s2-1), opts)
   ptb_m2r  = stat_dispersion(all_rvp(s2:s3-1), opts)
  elseif (boxv .eq. "ten") then 
   ptb_ring = stat_dispersion(all_ten(0:s0-1), opts)
   ptb_r2m  = stat_dispersion(all_ten(s0:s1-1), opts)
   ptb_mono = stat_dispersion(all_ten(s1:s2-1), opts)
   ptb_m2r  = stat_dispersion(all_ten(s2:s3-1), opts)
  elseif (boxv .eq. "flux") then 
   ptb_ring = stat_dispersion(all_flux(0:s0-1), opts)
   ptb_r2m  = stat_dispersion(all_flux(s0:s1-1), opts)
   ptb_mono = stat_dispersion(all_flux(s1:s2-1), opts)
   ptb_m2r  = stat_dispersion(all_flux(s2:s3-1), opts)
  end if 

; Ringlike phase                                                                             
  box_y(0,0) = ptb_ring(2)
  box_y(0,1) = ptb_ring(6)
  box_y(0,2) = ptb_ring(8)
  box_y(0,3) = ptb_ring(10)
  box_y(0,4) = ptb_ring(14)

; Ringlike to monopole transition                                                            
  box_y(1,0) = ptb_r2m(2)
  box_y(1,1) = ptb_r2m(6)
  box_y(1,2) = ptb_r2m(8)
  box_y(1,3) = ptb_r2m(10)
  box_y(1,4) = ptb_r2m(14)

; Monopole phase                                                                             
  box_y(2,0) = ptb_mono(2)
  box_y(2,1) = ptb_mono(6)
  box_y(2,2) = ptb_mono(8)
  box_y(2,3) = ptb_mono(10)
  box_y(2,4) = ptb_mono(14)

; Monopole to ringlike transition                                                            
  box_y(3,0) = ptb_m2r(2)
  box_y(3,1) = ptb_m2r(6)
  box_y(3,2) = ptb_m2r(8)
  box_y(3,3) = ptb_m2r(10)
  box_y(3,4) = ptb_m2r(14)

  print(box_y(:,2) )

;===================================================
; Set output file path, and produce the box plot
;===================================================

  if (miss .eq. 1) then
   mstr = "miss1"
  else
   mstr = ""
  end if

  if (sub .eq. 1) then
   sstr = "_sub"+t0+"_"
  else
   sstr = ""
  end if

; Choose variable to plot 
  if (boxv .eq. "v") then 

 ; Maximum tendency of mean tangential wind (one level) 
   if (zs .eq. "True") then 
   ; Tangential wind
    if (wind .eq. "vtan") then 
     out_box = "$nep/nepartak/images/vc_boxplot_group"+grp+"_sc"+sc+"_"+var+\
               sstr+z0+"_"+mstr+"_"+thr+"_"+pstr
   ; Total wind
    else
     out_box = "$nep/nepartak/images/vc_boxplot_group"+grp+"_sc"+sc+"_"+var\
               +wind+sstr+z0+"_"+mstr+"_"+thr+"_"+pstr
    end if 

 ; 14/12/2019 --> tendency of maximum tangential wind (any level)
   elseif (zs .eq. "max") then
     out_box = "$nep/nepartak/images/vc_boxplot_group"+grp+"_sc"+sc+"_vtan_max"
 
 ; Maximum tendency of mean tangential wind (any level)
   else
    out_box  = "$nep/nepartak/images/vc_boxplot_group"+grp+"_sc"+sc+"_"+var+\
               sstr+"_"+mstr+"_"+thr+"_"+pstr
   end if 
; 31/10/2019 --> may need further editing to take into account 't0'
  else
    out_box  = "$nep/nepartak/images/vc_boxplot_group"+grp+"_sc"+sc+"_"+boxv+\
               "_"+mstr+"_"+thr+"_"+pstr
  end if 

; Open workstation for plotting 
  wks = gsn_open_wks(opt, out_box)

; Create set of integer values denoting box label positions                                  
  x0  = ispan(1,4,1)

; General box plot resources                                                                 
  pres = True
  pres@tmXBLabels         = (/"Ring","R to M","Mono","M to R"/)
  pres@gsnMaximize        = True
  pres@tiYAxisFontHeightF = 0.0225
;  pres@tiYAxisString      = "Tendency of maximum mean tangential wind (m s~S~-1~N~)"
  pres@tiYAxisString      = "dvtan~B~max~N~/dt (m s~S~-1~N~ h~S~-1~N~)"

; Move y-axis title slightly left (away from plot)
  pres@tiYAxisOffsetXF    = -0.01

  if (boxv .eq. "rvp") then 
   pres@trYMaxF  = 1.00
   pres@trYMinF  = 0.40
  elseif (boxv .eq. "ten") then 
   pres@trYMaxF  = 2.0
   pres@trYMinF  = -2.0  
  elseif (boxv .eq. "v") then 
   pres@trYMaxF  = 5.0
   pres@trYMinF  = -5.0
  end if 

; Box resources                                                                              
  bres = True
  bres@boxWidth  = 0.5
  bres@boxColors = (/"royalblue1", "purple3", "red3", "orange3"/)

; Line resources (increase line thickness)
  lres = True
  lres@gsLineThicknessF   = 3.0

; Produce the boxplot                                                                        
; Rightmost dimension needs to contain box plot reference data (min,median,max,...)          
  plot = boxplot(wks, x0, box_y(:,:), bres, pres, lres)

; Add a horizontal line through y = 0
  opts_l                   = True 
  opts_l@gsLineThicknessF  = 2.0
  opts_l@gsLineDashPattern = 1
  opts_l@gsLineColor       = "black"

  axis = gsn_add_polyline(wks, plot, (/0, 5/), (/0, 0/), opts_l)

  draw(plot)
  frame(wks)

  print(out_box)

end 
