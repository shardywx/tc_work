; Calculate 1-h wind speed tendencies during all 24-h RI periods in our [4 x 12] MetUM simulations

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

; ncl opt=\"x11\" sub=0 t0=1 mth=\"p2\" trb=0 tend=2 plt=\"all\" mlev0=17 mlev1=46 alt=1 sc=5
; box=1 lgd0=\"tr\" grp=\"7\" wt0=1 wt=0.5 w0=0.3 r0=15 r1=30 diff=0 vc_mflux_scatter.ncl

; 'opt'   = output file format ("pdf", "x11", etc)
; 'mth'   = method for calculating storm centre
; 'trb'   = turn troubleshooting on (1) or off (0)
; 'tend'  = method for calculating centred differences (2 or 1)
; 'alt'   = use original (0) or alternative (1) VC phase data
; 'plt'   = all phases ("all"); main four VC phases ("vc"); VC/No-VC ("no")
; 'sub'   = use a subset of the data points (e.g. only positive tendencies)
; 't0'    = wind speed / MSLP tendency threshold (e.g. 0,1,2), if 'sub=1'  
; 'wt0'   = use threshold method (1) to constrain radial integration for mass flux calculations

begin

;=========================================
; Now read in the data from text files 
;=========================================

; Follow notation of scripts including 'tc_mbgt_tan_cxz_paper.ncl'
; 11/10/2019. Group 5 data --> vc_arr_all(:,0). Group 6 data --> vc_arr_all(:,1)
  vc_arr_all = (/ (/56,38,51,30,22/), (/55,35,49,27,23/), \
                  (/69,39,64,37,27/), (/52,33,43,26,19/) /)

; Link the chosen group to the way we read in data below  
  if (grp .eq. "5") then 
   vc_ind = 0
  elseif (grp .eq. "6") then 
   vc_ind = 1
  elseif (grp .eq. "7") then 
   vc_ind = 2
  elseif (grp .eq. "8") then
   vc_ind = 3
  elseif (grp .eq. "9") then 
   vc_ind = 4
  end if 

; Read in correct data for all VC phases 
  dat_ring  = asciiread("$ar/text/vc_sym_dat_group"+grp+".txt",(/vc_arr_all(0,vc_ind)/),"string")
  dat_r2m   = asciiread("$ar/text/vc_s2a_dat_group"+grp+".txt",(/vc_arr_all(1,vc_ind)/),"string")
  dat_mono  = asciiread("$ar/text/vc_asym_dat_group"+grp+".txt",(/vc_arr_all(2,vc_ind)/),"string")
  dat_m2r   = asciiread("$ar/text/vc_a2s_dat_group"+grp+".txt",(/vc_arr_all(3,vc_ind)/),"string")

  ens0_ring = asciiread("$ar/text/vc_sym_sim_group"+grp+".txt",(/vc_arr_all(0,vc_ind)/),"string")
  ens0_r2m  = asciiread("$ar/text/vc_s2a_sim_group"+grp+".txt",(/vc_arr_all(1,vc_ind)/),"string")
  ens0_mono = asciiread("$ar/text/vc_asym_sim_group"+grp+".txt",(/vc_arr_all(2,vc_ind)/),"string")
  ens0_m2r  = asciiread("$ar/text/vc_a2s_sim_group"+grp+".txt",(/vc_arr_all(3,vc_ind)/),"string")

  ts0_ring  = asciiread("$ar/text/vc_sym_ts_group"+grp+".txt",(/vc_arr_all(0,vc_ind)/),"integer")
  ts0_r2m   = asciiread("$ar/text/vc_s2a_ts_group"+grp+".txt",(/vc_arr_all(1,vc_ind)/),"integer")
  ts0_mono  = asciiread("$ar/text/vc_asym_ts_group"+grp+".txt",(/vc_arr_all(2,vc_ind)/),"integer")
  ts0_m2r   = asciiread("$ar/text/vc_a2s_ts_group"+grp+".txt",(/vc_arr_all(3,vc_ind)/),"integer")

  tf0_ring  = asciiread("$ar/text/vc_sym_tf_group"+grp+".txt",(/vc_arr_all(0,vc_ind)/),"integer")
  tf0_r2m   = asciiread("$ar/text/vc_s2a_tf_group"+grp+".txt",(/vc_arr_all(1,vc_ind)/),"integer")
  tf0_mono  = asciiread("$ar/text/vc_asym_tf_group"+grp+".txt",(/vc_arr_all(2,vc_ind)/),"integer")
  tf0_m2r   = asciiread("$ar/text/vc_a2s_tf_group"+grp+".txt",(/vc_arr_all(3,vc_ind)/),"integer")

; Calculate maximum number of times to read in for any one simulation, for each VC phase
  td_ring   = (tf0_ring - ts0_ring) + 1
  nts_ring  = max(td_ring)

  td_r2m    = (tf0_r2m - ts0_r2m) + 1
  nts_r2m   = max(td_r2m)

  td_mono   = (tf0_mono - ts0_mono) + 1
  nts_mono  = max(td_mono)

  td_m2r    = (tf0_m2r - ts0_m2r) + 1
  nts_m2r   = max(td_m2r)

; Now create new arrays to hold all data for each VC phase
  var_ring  = new( (/vc_arr_all(0,vc_ind), nts_ring, 3/), "float")
  var_r2m   = new( (/vc_arr_all(1,vc_ind), nts_r2m, 3/), "float")
  var_mono  = new( (/vc_arr_all(2,vc_ind), nts_mono, 3/), "float")
  var_m2r   = new( (/vc_arr_all(3,vc_ind), nts_m2r, 3/), "float")

; Dimension size information for each VC phase
  dsize     = (/ vc_arr_all(0,vc_ind), vc_arr_all(1,vc_ind), \
                 vc_arr_all(2,vc_ind), vc_arr_all(3,vc_ind) /)

; Add metadata to arrays above 
  var_ring!0    = "sim"
  var_ring!1    = "time"
  var_ring!2    = "diag"
  var_ring&diag = (/"Mass flux", "RVP","RVP tendency"/)

  var_r2m!0     = "sim"
  var_r2m!1     = "time"
  var_r2m!2     = "diag"
  var_r2m&diag  = (/"Mass flux", "RVP","RVP tendency"/)

  var_mono!0    = "sim"
  var_mono!1    = "time"
  var_mono!2    = "diag"
  var_mono&diag = (/"Mass flux", "RVP","RVP tendency"/)

  var_m2r!0     = "sim"
  var_m2r!1     = "time"
  var_m2r!2     = "diag"
  var_m2r&diag  = (/"Mass flux", "RVP","RVP tendency"/)

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

;============================
; Loop over ringlike phase
;============================

  do st = 0, dsize(0)-1 

 ; Radial vorticity profile (119 values)
   diri_rvp = "$ar/text/"+dat_ring(st)+"_"+ens0_ring(st)+"_sc"+sc+"_inner_core.txt"

 ; Radial vorticity profile tendency (119 values)
   diri_ten = "$ar/text/"+dat_ring(st)+"_"+ens0_ring(st)+"_sc"+sc+"_inner_core_grad.txt"

 ; Radial vorticity profile (and tendency)
   rvp_all     = asciiread(diri_rvp, ntimes_vc, "float")
   ten_all     = asciiread(diri_ten, ntimes_vc, "float")

 ; Finally, all data into correct arrays [radial vorticity profile]
   npts = ( tf0_ring(st) - ts0_ring(st) )
   var_ring(st,0:npts,1) = (/rvp_all( ts0_ring(st):tf0_ring(st) )/)
   var_ring(st,0:npts,2) = (/ten_all( ts0_ring(st):tf0_ring(st) )/)

  end do  

  delete([/diri_rvp, diri_ten/])

;=============================================
; Loop over ringlike to monopole transition
;=============================================

  do st = 0, dsize(1)-1

 ; Radial vorticity profile (119 values)                                                    
   diri_rvp = "$ar/text/"+dat_r2m(st)+"_"+ens0_r2m(st)+"_sc"+sc+"_inner_core.txt"

 ; Radial vorticity profile tendency (119 values)                                           
   diri_ten = "$ar/text/"+dat_r2m(st)+"_"+ens0_r2m(st)+"_sc"+sc+"_inner_core_grad.txt"

 ; Radial vorticity profile (and tendency)                                                 
   rvp_all     = asciiread(diri_rvp, ntimes_vc, "float")
   ten_all     = asciiread(diri_ten, ntimes_vc, "float")

 ; Finally, all data into correct arrays [radial vorticity profile]
   npts = ( tf0_r2m(st) - ts0_r2m(st) )
   var_r2m(st,0:npts,1) = (/rvp_all( ts0_r2m(st):tf0_r2m(st) )/)
   var_r2m(st,0:npts,2) = (/ten_all( ts0_r2m(st):tf0_r2m(st) )/)

  end do
  delete([/diri_rvp, diri_ten/])

;============================
; Loop over monopole phase 
;============================

  do st = 0, dsize(2)-1

 ; Radial vorticity profile (119 values)                                                   
   diri_rvp = "$ar/text/"+dat_mono(st)+"_"+ens0_mono(st)+"_sc"+sc+"_inner_core.txt"

 ; Radial vorticity profile tendency (119 values)                                          
   diri_ten = "$ar/text/"+dat_mono(st)+"_"+ens0_mono(st)+"_sc"+sc+"_inner_core_grad.txt"

 ; Radial vorticity profile (and tendency)                                                
   rvp_all     = asciiread(diri_rvp, ntimes_vc, "float")
   ten_all     = asciiread(diri_ten, ntimes_vc, "float")

 ; Finally, all data into correct arrays [radial vorticity profile]
   npts = ( tf0_mono(st) - ts0_mono(st) )
   var_mono(st,0:npts,1) = (/rvp_all( ts0_mono(st):tf0_mono(st) )/)
   var_mono(st,0:npts,2) = (/ten_all( ts0_mono(st):tf0_mono(st) )/)

  end do
  delete([/diri_rvp, diri_ten/])

;=============================================
; Loop over monopole to ringlike transition                                        
;=============================================

  do st = 0, dsize(3)-1

 ; Radial vorticity profile (119 values)                                                    
   diri_rvp = "$ar/text/"+dat_m2r(st)+"_"+ens0_m2r(st)+"_sc"+sc+"_inner_core.txt"

 ; Radial vorticity profile tendency (119 values)                                           
   diri_ten = "$ar/text/"+dat_m2r(st)+"_"+ens0_m2r(st)+"_sc"+sc+"_inner_core_grad.txt"

 ; Radial vorticity profile (and tendency)                                                 
   rvp_all     = asciiread(diri_rvp, ntimes_vc, "float")
   ten_all     = asciiread(diri_ten, ntimes_vc, "float")

 ; Finally, all data into correct arrays [radial vorticity profile]
   npts = ( tf0_m2r(st) - ts0_m2r(st) )
   var_m2r(st,0:npts,1) = (/rvp_all( ts0_m2r(st):tf0_m2r(st) )/)
   var_m2r(st,0:npts,2) = (/ten_all( ts0_m2r(st):tf0_m2r(st) )/)

  end do
  delete([/diri_rvp, diri_ten/])

;=============================================
; Reshape to form 1D arrays before plotting 
;=============================================

 ; Ringlike phase 
   ring_rvp = ndtooned(var_ring(:,:,2) )
   r_size   = dimsizes(ring_rvp)

 ; Ringlike to monopole transition
   r2m_rvp  = ndtooned(var_r2m(:,:,2) )
   r2m_size = dimsizes(r2m_rvp)

 ; Monopole phase 
   mono_rvp = ndtooned(var_mono(:,:,2) )
   m_size   = dimsizes(mono_rvp)

 ; Monopole to ringlike transition
   m2r_rvp  = ndtooned(var_m2r(:,:,2) )
   m2r_size = dimsizes(m2r_rvp)

;==========================================
; Read in mass flux data from text files 
;==========================================

 ; Either read in vertical velocity data (troubleshooting)
   if (trb .eq. 1) then 
    if (r0 .eq. 10) then 
     fili     = (/"group"+grp+"_vvel_ring_sc"+sc+"_lay1_"+hy0+"_"\
                  +hy1+"_w"+w0+"_"+r0+"km_"+r1+"km.txt",\
                  "group"+grp+"_vvel_r2m_sc"+sc+"_lay1_"+hy0+"_"\
                  +hy1+"_w"+w0+"_"+r0+"km_"+r1+"km.txt",\
                  "group"+grp+"_vvel_mono_sc"+sc+"_lay1_"+hy0+"_"\
                  +hy1+"_w"+w0+"_"+r0+"km_"+r1+"km.txt",\
                  "group"+grp+"_vvel_m2r_sc"+sc+"_lay1_"+hy0+"_"\
                  +hy1+"_w"+w0+"_"+r0+"km_"+r1+"km.txt"/)
    elseif (r0 .eq. 15 .and. r1 .eq. 30) then 
     r2 = 10
     r3 = 25 
     fili     = (/"group"+grp+"_vvel_ring_sc"+sc+"_lay1_"+hy0+\
                  "_"+hy1+"_w"+w0+"_"+r0+"km_"+r1+"km.txt",\
                  "group"+grp+"_vvel_r2m_sc"+sc+"_lay1_"+hy0+\
                  "_"+hy1+"_w"+w0+"_"+r2+"km_"+r3+"km.txt",\
                  "group"+grp+"_vvel_mono_sc"+sc+"_lay1_"+hy0+\
                  "_"+hy1+"_w"+w0+"_"+r2+"km_"+r3+"km.txt",\
                 "group"+grp+"_vvel_m2r_sc"+sc+"_lay1_"+hy0+\
                 "_"+hy1+"_w"+w0+"_"+r0+"km_"+r1+"km.txt"/)  
    end if 
 ; Or, read in mass flux data (28/10/2019 --> analysis for paper)
   else

  ; No vertical velocity threshold (28/10/2019 --> unlikely to use this option)
    if (wt0 .ne. 1) then 

     if (r0 .eq. 10) then
      fili     = (/"group"+grp+"_ring_sc"+sc+"_lay1_"+hy0+"_"+hy1+"_w"+w0+"_"+r0+"km_"+r1+"km.txt",\
                   "group"+grp+"_r2m_sc"+sc+"_lay1_"+hy0+"_"+hy1+"_w"+w0+"_"+r0+"km_"+r1+"km.txt",\
                   "group"+grp+"_mono_sc"+sc+"_lay1_"+hy0+"_"+hy1+"_w"+w0+"_"+r0+"km_"+r1+"km.txt",\
                   "group"+grp+"_m2r_sc"+sc+"_lay1_"+hy0+"_"+hy1+"_w"+w0+"_"+r0+"km_"+r1+"km.txt"/)
     elseif (r0 .eq. 15 .and. r1 .eq. 30) then
      r2 = 10
      r3 = 25
      fili     = (/"group"+grp+"_ring_sc"+sc+"_lay1_"+hy0+"_"+hy1+"_w"+w0+"_"+r0+"km_"+r1+"km.txt",\
                   "group"+grp+"_r2m_sc"+sc+"_lay1_"+hy0+"_"+hy1+"_w"+w0+"_"+r2+"km_"+r3+"km.txt",\
                   "group"+grp+"_mono_sc"+sc+"_lay1_"+hy0+"_"+hy1+"_w"+w0+"_"+r2+"km_"+r3+"km.txt",\
                  "group"+grp+"_m2r_sc"+sc+"_lay1_"+hy0+"_"+hy1+"_w"+w0+"_"+r0+"km_"+r1+"km.txt"/)
     end if

  ; Vertical velocity threshold (28/10/2019 --> most recent analysis)
    else

     if (diff .eq. 1) then 

    ; Difference in mass flux between two layers
      fili     = (/"group"+grp+"_ring_sc"+sc+"_lay1_"+hy0+\
                   "_"+hy1+"_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+".txt",\
                   "group"+grp+"_r2m_sc"+sc+"_lay1_"+hy0+\
                   "_"+hy1+"_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+".txt",\
                   "group"+grp+"_mono_sc"+sc+"_lay1_"+hy0+\
                   "_"+hy1+"_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+".txt",\
                   "group"+grp+"_m2r_sc"+sc+"_lay1_"+hy0+\
                   "_"+hy1+"_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+".txt"/)

     else

    ; Read in upper and lower layers separately (28/10/2019 --> most recent analysis)
      fili_l   = (/"group"+grp+"_ring_sc"+sc+"_lay1_"+hy0+\
                   "_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+".txt",\
                   "group"+grp+"_r2m_sc"+sc+"_lay1_"+hy0+\
                   "_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+".txt",\
                   "group"+grp+"_mono_sc"+sc+"_lay1_"+hy0+\
                   "_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+".txt",\
                   "group"+grp+"_m2r_sc"+sc+"_lay1_"+hy0+\
                   "_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+".txt"/) 

      fili_u   = (/"group"+grp+"_ring_sc"+sc+"_lay1_"+hy1+\
                   "_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+".txt",\
                   "group"+grp+"_r2m_sc"+sc+"_lay1_"+hy1+\
                   "_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+".txt",\
                   "group"+grp+"_mono_sc"+sc+"_lay1_"+hy1+\
                   "_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+".txt",\
                   "group"+grp+"_m2r_sc"+sc+"_lay1_"+hy1+\
                   "_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+".txt"/)

     end if 

    end if 

   end if 

 ; Read in data using 'asciiread'
   if (diff .eq. 0) then 

   ; Ringlike phase
     ring_u_flux = asciiread("$ar/text/"+fili_u(0),(/r_size/),"float")
     ring_l_flux = asciiread("$ar/text/"+fili_l(0),(/r_size/),"float")
   
   ; Ringlike to monopole transition
     r2m_u_flux  = asciiread("$ar/text/"+fili_u(1),(/r2m_size/),"float")
     r2m_l_flux  = asciiread("$ar/text/"+fili_l(1),(/r2m_size/),"float")

   ; Monopole phase  
     mono_u_flux = asciiread("$ar/text/"+fili_u(2),(/m_size/),"float")
     mono_l_flux = asciiread("$ar/text/"+fili_l(2),(/m_size/),"float")

   ; Monopole to ringlike transition  
     m2r_u_flux  = asciiread("$ar/text/"+fili_u(3),(/m2r_size/),"float")
     m2r_l_flux  = asciiread("$ar/text/"+fili_l(3),(/m2r_size/),"float")

   else

   ; Ringlike phase
     ring_flux = asciiread("$ar/text/"+fili(0),(/r_size/),"float")

   ; Ringlike to monopole transition
     r2m_flux  = asciiread("$ar/text/"+fili(1),(/r2m_size/),"float")

   ; Monopole phase
     mono_flux = asciiread("$ar/text/"+fili(2),(/m_size/),"float")

   ; Monopole to ringlike transition
     m2r_flux  = asciiread("$ar/text/"+fili(3),(/m2r_size/),"float")

   end if 

;===========================================================
; Before reading into arrays below, apply running average
;===========================================================

 ; Check that 'ring_u_flux' array is 1D --> otherwise, potentially change final option from '0' 
   if (diff .eq. 0) then 

  ; Ringlike phase 
    ring_u_ave = runave_n_Wrap(ring_u_flux, 3, 0, 0)
    ring_l_ave = runave_n_Wrap(ring_l_flux, 3, 0, 0)

  ; Ringlike to monopole transition
    r2m_u_ave  = runave_n_Wrap(r2m_u_flux, 3, 0, 0)
    r2m_l_ave  = runave_n_Wrap(r2m_l_flux, 3, 0, 0)

  ; Monopole phase 
    mono_u_ave = runave_n_Wrap(mono_u_flux, 3, 0, 0)
    mono_l_ave = runave_n_Wrap(mono_l_flux, 3, 0, 0)

  ; Monopole to ringlike transition
    m2r_u_ave  = runave_n_Wrap(m2r_u_flux, 3, 0, 0)
    m2r_l_ave  = runave_n_Wrap(m2r_l_flux, 3, 0, 0)

  ; Calculate difference between upper and lower layers 
    ring_flux  = ring_u_ave - ring_l_ave
    r2m_flux   = r2m_u_ave - r2m_l_ave
    mono_flux  = mono_u_ave - mono_l_ave
    m2r_flux   = m2r_u_ave - m2r_l_ave

   end if 

   print(ring_u_flux)
   exit()

;===========================================================
; Now continue and read all data into arrays for plotting 
;===========================================================

 ; Size information needed to create single array to hold all values 
   all_size = r_size + r2m_size + m_size + m2r_size
   
   s0       = r_size
   s1       = s0 + r2m_size
   s2       = s1 + m_size
   s3       = s2 + m2r_size

 ; Create new arrays to hold all values for each variable 
   all_fl   = new( (/all_size/), "float")
   all_rvp  = new( (/all_size/), "float")

 ; Fill each of these arrays with the values from each VC phase 
   all_fl(0:s0-1)     = ring_flux
   all_fl(s0:s1-1)    = r2m_flux
   all_fl(s1:s2-1)    = mono_flux
   all_fl(s2:s3-1)    = m2r_flux

   all_rvp(0:s0-1)    = ring_rvp
   all_rvp(s0:s1-1)   = r2m_rvp
   all_rvp(s1:s2-1)   = mono_rvp
   all_rvp(s2:s3-1)   = m2r_rvp

 ; Add metadata so we can colour-code the markers by VC phase
   all_fl!0   = "vc_phase"
   all_rvp!0  = "vc_phase"

 ; Change units before plotting 
   all_rvp    = all_rvp * 100
   all_rvp@units = "10~S~-2~N~ h~S~-1~N~"

;===========================================================
; If desired, remove selected data points before plotting 
;===========================================================

   if (sub .eq. 1) then 
    all_fl = where(all_fl .lt. t0, all_fl@_FillValue, all_fl)
   end if    

;================================================
; Open workstation and define output file path
;================================================

   if (wt0 .eq. 1) then 
    output = "$sam/nepartak/images/vc_mflux_scatter_group"+grp+"_sc"+sc+\
             "_sub"+t0+"_"+r0+"km_"+r1+"km_wt"+wt
   else
    output = "$sam/nepartak/images/vc_mflux_scatter_group"+grp+"_sc"+sc+\
             "_sub"+t0+"_"+r0+"km_"+r1+"km"
   end if 
   wks     = gsn_open_wks(opt, output)

;==================================================
; Now produce scatter plot using all data points 
;==================================================

 ; Plotting resources 
   vc_opts                       = True
   vc_opts@gsnDraw               = False
   vc_opts@gsnFrame              = False
   vc_opts@gsnMaximize           = False

 ; Y-axis
   if (trb .ne. 1) then 
    vc_opts@tiYAxisString        = "Mass flux (kg m~S~-2~N~ s~S~-1~N~)"
   else
    vc_opts@tiYAxisString        = "Vertical velocity (m s~S~-1~N~)"   
   end if 
   vc_opts@tiYAxisFontHeightF    = 0.0150

   if (sub .ne. 1) then  
    if (trb .ne. 1) then 
     vc_opts@trYMaxF             = 5.0
     vc_opts@trYMinF             = -5.0
    else
     vc_opts@trYMaxF             = 10.0
     vc_opts@trYMinF             = -10.0
    end if 
   else
    vc_opts@trYMaxF              = 5.0
    vc_opts@trYMinF              = t0
   end if 

 ; X-axis (define our own labels; set limits)
   vc_opts@tiXAxisString         = "Tendency of radial vorticity profile ("+all_rvp@units+")"
   vc_opts@tiXAxisFontHeightF    = 0.0150
   vc_opts@trXMinF               = -10.0
   vc_opts@trXMaxF               = 10.0

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
   vc_opts@vpXF                 = 0.12
   vc_opts@vpYF                 = 0.88
   vc_opts@vpWidthF             = 0.75
   vc_opts@vpHeightF            = 0.75

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

 ; Plot mass flux vs tendency of radial vorticity profile
   vc_plot_fl0   = gsn_csm_xy(wks, all_rvp(0:s0-1), all_fl(0:s0-1), vc_opts)
   fl1           = gsn_add_polymarker(wks, vc_plot_fl0, all_rvp(s0:s1-1), all_fl(s0:s1-1), vc_opts1)
   fl2           = gsn_add_polymarker(wks, vc_plot_fl0, all_rvp(s1:s2-1), all_fl(s1:s2-1), vc_opts2)
   fl3           = gsn_add_polymarker(wks, vc_plot_fl0, all_rvp(s2:s3-1), all_fl(s2:s3-1), vc_opts3)

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
   annotate                           = gsn_add_annotation(vc_plot_fl0, lgd, am_opts)

 ; Now switch to Cartesian grid in both panels
   c_axis_uv  = cartesian_axis(wks, vc_plot_fl0)
 
 ; Finally, draw selected panel(s)
   draw(vc_plot_fl0)
   frame(wks)

;===================================
; Also draw box plot, if desired 
;===================================

 if (box .eq. 1) then 

; Create array to hold all box plot values                                                  
; 4 VC phases                                                                               
; 5 values per phase [minimum, lower quartile, median, upper quartile, maximum]             
  box_y    = new( (/4,5/), "float")

; Options for 'stat_dispersion' function                                                    
  opts     = True
  opts@PrintStat = False

; Choose variable to analyse (14/10/2019 --> possibly develop this part of script)
  plot_int0 = all_fl(:)
  
; Calculate stats for each VC phase 
  ptb_ring  = stat_dispersion(all_fl(0:s0-1), opts)
  ptb_r2m   = stat_dispersion(all_fl(s0:s1-1), opts)
  ptb_mono  = stat_dispersion(all_fl(s1:s2-1), opts)
  ptb_m2r   = stat_dispersion(all_fl(s2:s3-1), opts)

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

;========================                                              
; Produce the box plot                                                 
;========================     

; Set output file path
  out_box = "$sam/nepartak/images/vc_mflux_boxplot_group"+grp+"_sc"+sc+\
            "_sub"+t0+"_"+r0+"km_"+r1+"km"
  wks     = gsn_open_wks(opt, out_box)

; Create set of integer values denoting box label positions                         
  x0  = ispan(1,4,1)

; General box plot resources                                                        
  pres = True
  pres@tmXBLabels      = (/"Ring","R to M","Mono","M to R"/)
  pres@gsnMaximize     = True
  pres@tiYAxisFontHeightF = 0.020

; Box resources                                                                     
  bres = True
  bres@boxWidth  = 0.5

; Line resources                                                                    
  lres = False

; Produce the boxplot                                                               
; Rightmost dimension needs to contain box plot reference data (min,median,max,...)
  plot = boxplot(wks, x0, box_y(:,:), bres, pres, lres)

  draw(plot)
  frame(wks)

 end if 

end 
