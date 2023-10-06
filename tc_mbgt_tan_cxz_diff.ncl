; Produce plots of tangential momentum budget terms for VCs and ERCs 

; Run using:                                                                                  

; ncl opt=\"png\" dist=2.0 nt=73 nr=11 ar=0.5 sc=5 rad0=0.5 w_check=1 thr=2 grp=8 z0=0
; w0=0.3 zoom=0 typ=\"ring\" paper=1 res=\"5min\" ts=36 tf=48 
; col=3 row=3 p0=\"ring\" p1=\"mono\" run=1 tc_mbgt_tan_cxz_diff.ncl   

; 'opt'   = output file format ("pdf" or "x11")                                               
; 'dist'  = size of box following storm (degrees)
; 'typ'   = VC phase ("sym","asym","stoa","atos")
; 'thr'   = if plotting "v" or "slp" data, threshold above/below which we're plotting (e.g. 2/-2)
; 'grp'   = chosen set of data points to analyse: ("4", "4a", "4b", "4c", "5", "6")
; 'ar'    = distance in degrees of outer radius (1.0, 2.0, etc)                               
; 'nr'    = number of radial circles between r = 0 and r = 'ar' (21, 41, 61, etc)             
; 'nt'    = number of azimuth angles in cylindrical grid                                      
; 'cn0'   = offline storm centre: "slp", "vort", "geo", "geo_sm"                              
; 'rad0'  = radius (degrees) for final plots

; 'z0'    = 0 (0.002 / 0.005), 6 (6.128), 8 (8.011), 10 (10.624)

; Select data to plot
; "ring"      = ringlike phase
; "r2m"       = ringlike to monopole transition
; "mono"      = monopole phase
; "m2r"       = monopole to ringlike transition
; "mean"      = mean of all 4 phases

; "mean_r"    = ring - mean
; "mean_r2"   = r2m  - mean
; "mean_m"    = mono - mean
; "mean_m2"   = m2r  - mean

; "r2m_diff"  = r2m  - ring 
; "mono_diff" = mono - r2m 
; "m2r_diff"  = m2r  - mono
; "ring_diff" = ring - m2r

; "novc_diff" = mean - novc

; "paper0_diff" = mono - ring (individual terms)
; "paper1_diff" = ring - mono (individual terms)
; "paper2_diff" = mono - ring (combined terms)

; 'pc' stream --> 'it'                                                                        
; 'centre'    --> 'it-1'                                                                      
; 'pd' stream --> 'it-2'                                                                      

; Load main NCL functions and procedures                                                      
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"

; Also load user-defined functions and procedures                                             
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_centre.ncl"
load "$sam/ncl_func/setup_cyl.ncl"
load "$sam/ncl_func/draw_cyl.ncl"

begin

;===================================                                                          
; Reading in data from text files                                                             
;===================================                                                          

; Read in any 'pc' stream file to retrieve vertical level metadata 
  if (res .eq. "1h") then 
   diri      = "/nfs/a319/earshar/02T12/em11"
  else
   diri      = "/nfs/a319/earshar/02T12/em11/vc/"
  end if 

  fili_pr    = "20160702T1200Z_NPTK_4p4_L80_ra1t_em11"
  fili_c1    = systemfunc("cd "+diri+" ; ls "+fili_pr+"_pc*.nc")
  fili_c     = diri+"/"+fili_c1
  c          = addfile(fili_c,"r")

; Read in model level array
  hyb_levs   = c->hybrid_ht(:)
  num_levs   = dimsizes(hyb_levs)

  hyb_plot   = hyb_levs / 1000
  hyb_plot@units = "km"

; Define some cylindrical coordinate arrays 
  pii        = 3.14159265
  radii      = fspan(0.0, 1.0, 21)
  num_rad    = dimsizes(radii)

; Array of strings describing diagnostics                                                     
  diag_arr   = (/"lhs1","rhs1","rhs2","rhs3","rhs4","rhs5a","rhs5b","rhsT",\
                 "mean","eddy","vflux","vadv","vrad","vtan","vvel","aam"/)
  num_diag   = dimsizes(diag_arr)

; Modify string if we've used running average when calculating storm-relative winds 
  if (run .eq. 1) then 
    diag_arr = diag_arr+"_run12"
  end if 

; Array of strings describing VC phase
  if (typ .eq. "v") then 
   vc_arr    = (/"v"/)
  elseif (typ .eq. "slp") then 
   vc_arr    = (/"slp"/)
  else
   vc_arr    = (/"ring", "r2m", "mono", "m2r"/)
  end if 

; 'sc' array (sc always = 5 for 'mono')
  sc_arr     = (/sc,sc,5,sc/)

  num_vc     = dimsizes(vc_arr)

; Minimum height at which the threshold was applied (EDIT later if required)
; 09/10/2019 --> difference in 'zstr' results from using 'pc' vs 'pd' stream height data 
  if (z0 .eq. 0) then 
   zstr      = "0.005km"
  elseif (z0 .eq. 6) then 
   zstr      = "6.128km"
  elseif (z0 .eq. 8) then 
   zstr      = "8.011km"
  elseif (z0 .eq. 10) then 
   zstr      = "10.624km"
  end if 
  
; Time string information to help read in data  
  if (each .eq. 1) then
   if (typ .eq. "ring") then
    time_vc = ispan(647,670,1)
   elseif (typ .eq. "r2m") then
    time_vc = ispan(683,706,1)
   elseif (typ .eq. "mono") then
    time_vc = ispan(743,778,1)
   elseif (typ .eq. "m2r") then
    time_vc = ispan(839,862,1)
   end if 
 ; 'num_vc' here is the number of individual times (01/04/20)
 ; 'num_vc' in the original script is the number of VC phases
   num_vc = dimsizes(time_vc)
  end if

; Define 4D array to hold data from all VC phases, and for all diagnostics 
; [4 * 16 * 21 * 59] for regular VC phase analysis; [1 * 16 * 21 * 59] for 'v' or 'slp'
  mbgt_arr  = new( (/num_vc, num_diag, num_levs, num_rad/), float)

;===========================================================================
; Read in 2D composite arrays [59 * 21] of velocity vector and MBGT terms
;===========================================================================

; File path depends on whether we look at 1-h or 5-min data
  if (res .eq. "1h") then 
   diri_text = "$ar/text/"
  else
   diri_text = "$ar/text/5min_"
  end if 

; If plotting tail of distribution data...
  if (typ .eq. "v" .or. typ .eq. "slp") then 
   print("Working on : "+typ+" dataset")
   do diag   = 0, num_diag-1
  ; Read in 2D array [59 * 22]
  ; If using this option again, remember arrays are different shape to the newer ones (01/04/20)
    text_arr = asciiread(diri_text+diag_arr(diag)+"_sc"+sc+\
                         "_group"+grp+"_"+w0+"w_"+zstr+"_"+typ+thr+".txt", \
                         (/num_levs, num_rad+1/), "float")

  ; Ignore the leftmost column (index), and input everything else 
    mbgt_arr(0,diag,:,:) = text_arr(:,1:num_rad)
   end do

; If plotting regular VC phase data
  else

   do vc      = 0, num_vc-1
    print("Working on VC phase: "+vc_arr(vc) )
    do diag   = 0, num_diag-1

   ; Read in 2D array [59 * 22]
     if (res .eq. "1h" .and. grp .ne. "10" .and. grp .ne. "8") then 
      text_arr = asciiread(diri_text+diag_arr(diag)+"_sc"+sc_arr(vc)+\
                           "_group"+grp+"_"+w0+"w_"+zstr+"_"+vc_arr(vc)+".txt", \
                           (/num_levs, num_rad+1/), "float")

   ; 25/12/2019 --> have edited the text files to contain the correct number of columns
     else
      text_arr = asciiread(diri_text+diag_arr(diag)+"_sc"+sc_arr(vc)+\
                           "_group"+grp+"_"+w0+"w_"+zstr+"_"+vc_arr(vc)+".txt", \
                           (/num_levs, num_rad/), "float")
     end if 

   ; Ignore the leftmost column (index), and input everything else (only for 1-h data)
     if (res .eq. "1h" .and. grp .ne. "10" .and. grp .ne. "8") then 
      mbgt_arr(vc,diag,:,:) = text_arr(:,1:num_rad)
     else
      mbgt_arr(vc,diag,:,:) = text_arr(:,:)      
     end if 

    end do  ; End loop over diagnostics
   end do   ; End loop over VC phase

  end if 

; Define 3D array to hold data from 'non-VC' phase, for all diagnostics                     
; [16 * 21 * 59]
  novc_arr  = new( (/num_diag, num_levs, num_rad/), float)

; Also read in 'novc' data 
; [59 * 21] --> edited text files to contain correct number of columns (25/12/2019) 
; 12/12/2019 --> currently don't have 'novc' data for sc=4
  if (res .eq. "1h" .and. sc .eq. 5) then 
   do diag    = 0, num_diag-1
    novc_arr0 = asciiread(diri_text+diag_arr(diag)+"_sc"+sc+\
                          "_"+w0+"w_"+zstr+"_novc_"+ts+"_"+tf+".txt", \
                          (/num_levs, num_rad/), "float")
    novc_arr(diag,:,:) = novc_arr0(:,:)
   end do
  end if 

; Add metadata
  mbgt_arr!0    = "vc"
  mbgt_arr&vc   = vc_arr
  mbgt_arr!1    = "diag"
  mbgt_arr&diag = diag_arr
  mbgt_arr!2    = "lev"
  mbgt_arr&lev  = hyb_plot
  mbgt_arr!3    = "rad"
  mbgt_arr&rad  = radii

;==========================================
; Calculate mean fields etc, if required 
;==========================================

; Where are the data for all these terms located?                                             
; (0) LHS_1, (1) RHS_1, (2) RHS_2, (3) RHS_3, (4) RHS_4, (5) RHS_5A, (6) RHS_5B, (7) RHS_T    
; (8) MEAN, (9) EDDY, (10) VFLUX, (11) VADV, (12) VRAD, (13) VTAN, (14) VVEL, (15) AAM 

 if (typ .eq. "v") then 

  plot_arr = mbgt_arr(0,:,:,:)
  plot_tt  = "Wind speed tendency data"
 
 elseif (typ .eq. "slp") then 

  plot_arr = mbgt_arr(0,:,:,:)
  plot_tt  = "MSLP tendency data"

 else

; Split into individal VC phases again before plotting 
  ring_arr      = mbgt_arr(0,:,:,:)
  r2m_arr       = mbgt_arr(1,:,:,:)
  mono_arr      = mbgt_arr(2,:,:,:)
  m2r_arr       = mbgt_arr(3,:,:,:)

; Calculate mean + eddy terms for each of these phases
  if (typ .eq. "mean_eddy" .or. typ .eq. "paper2_diff" .or. typ .eq. "novc_mean" .or. \
      typ .eq. "novc_ring" .or. typ .eq. "novc_r2m" .or. typ .eq. "novc_mono" .or. \
      typ .eq. "novc_m2r") then 

  ; Ring-like phase
    if (fric .eq. 0) then 
     ring_sum    = ring_arr(8,:,:) + ring_arr(9,:,:)
    else 
     ring_sum    = ring_arr(8,:,:) + ring_arr(9,:,:) + ring_arr(6,:,:)
    end if 
    ring_sum!0   = "lev"
    ring_sum&lev = hyb_plot
    ring_sum!1   = "rad"
    ring_sum&rad = radii

  ; Ring-like to monopole transition
    if (fric .eq. 0) then 
     r2m_sum    = r2m_arr(8,:,:) + r2m_arr(9,:,:)
    else
     r2m_sum    = r2m_arr(8,:,:) + r2m_arr(9,:,:) + r2m_arr(6,:,:)
    end if 
    r2m_sum!0   = "lev"
    r2m_sum&lev = hyb_plot
    r2m_sum!1   = "rad"
    r2m_sum&rad = radii

  ; Monopole phase
    if (fric .eq. 0) then 
     mono_sum    = mono_arr(8,:,:) + mono_arr(9,:,:)
    else
     mono_sum    = mono_arr(8,:,:) + mono_arr(9,:,:) + mono_arr(6,:,:)
    end if 
    mono_sum!0   = "lev"
    mono_sum&lev = hyb_plot
    mono_sum!1   = "rad"
    mono_sum&rad = radii

  ; Monopole to ring-like transition
    if (fric .eq. 0) then 
     m2r_sum    = m2r_arr(8,:,:) + m2r_arr(9,:,:)
    else
     m2r_sum    = m2r_arr(8,:,:) + m2r_arr(9,:,:) + m2r_arr(6,:,:)
    end if 
    m2r_sum!0   = "lev"
    m2r_sum&lev = hyb_plot
    m2r_sum!1   = "rad"
    m2r_sum&rad = radii

  ; Non-VC composite 
    if (fric .eq. 0) then 
     novc_sum    = novc_arr(8,:,:) + novc_arr(9,:,:)
    else
     novc_sum    = novc_arr(8,:,:) + novc_arr(9,:,:) + novc_arr(6,:,:)
    end if  
    novc_sum!0   = "lev"
    novc_sum&lev = hyb_plot
    novc_sum!1   = "rad"
    novc_sum&rad = radii

  end if 

; Difference between ring-like + monopole phases for combined mean/eddy terms
  if (typ .eq. "paper2_diff") then 

  ; Mean + eddy
    if (p0 .eq. "ring") then 
     diff_sum     = mono_sum - ring_sum
    elseif (p0 .eq. "m2r") then 
     diff_sum     = mono_sum - m2r_sum
    elseif (p0 .eq. "r2m") then 
     diff_sum     = mono_sum - r2m_sum
    end if 
    diff_sum!0    = "lev"
    diff_sum&lev  = hyb_plot
    diff_sum!1    = "rad"
    diff_sum&rad  = radii

  ; Mean 
    if (p0 .eq. "ring") then 
     diff_mean    = mono_arr(8,:,:) - ring_arr(8,:,:)
    elseif (p0 .eq. "m2r") then 
     diff_mean    = mono_arr(8,:,:) - m2r_arr(8,:,:)
    elseif (p0 .eq. "r2m") then 
     diff_mean    = mono_arr(8,:,:) - r2m_arr(8,:,:)
    end if 
    diff_mean!0   = "lev"
    diff_mean&lev = hyb_plot
    diff_mean!1   = "rad"
    diff_mean&rad = radii

  ; Eddy
    if (p0 .eq."ring") then
     diff_eddy    = mono_arr(9,:,:) - ring_arr(9,:,:)
    elseif (p0 .eq. "m2r") then
     diff_eddy    = mono_arr(9,:,:) - r2m_arr(9,:,:)
    elseif (p0 .eq. "r2m") then
     diff_eddy    = mono_arr(9,:,:) - m2r_arr(9,:,:)
    end if 
    diff_eddy!0   = "lev"
    diff_eddy&lev = hyb_plot
    diff_eddy!1   = "rad"
    diff_eddy&rad = radii

  elseif (typ .eq. "novc_ring") then 

  ; Sum of mean + eddy 
    diff_sum      = ring_sum - novc_sum
    diff_sum!0    = "lev"
    diff_sum&lev  = hyb_plot
    diff_sum!1    = "rad"
    diff_sum&rad  = radii

  ; Mean
    diff_mean     = ring_arr(8,:,:) - novc_arr(8,:,:)
    diff_mean!0   = "lev"
    diff_mean&lev = hyb_plot
    diff_mean!1   = "rad"
    diff_mean&rad = radii

  ; Eddy
    diff_eddy     = ring_arr(9,:,:) - novc_arr(9,:,:)
    diff_eddy!0   = "lev"
    diff_eddy&lev = hyb_plot
    diff_eddy!1   = "rad"
    diff_eddy&rad = radii

  elseif (typ .eq. "novc_r2m") then 

  ; Sum of mean + eddy                                                              
    diff_sum      = r2m_sum - novc_sum
    diff_sum!0    = "lev"
    diff_sum&lev  = hyb_plot
    diff_sum!1    = "rad"
    diff_sum&rad  = radii

  ; Mean                                                                            
    diff_mean     = r2m_arr(8,:,:) - novc_arr(8,:,:)
    diff_mean!0   = "lev"
    diff_mean&lev = hyb_plot
    diff_mean!1   = "rad"
    diff_mean&rad = radii

  ; Eddy                                                                            
    diff_eddy     = r2m_arr(9,:,:) - novc_arr(9,:,:)
    diff_eddy!0   = "lev"
    diff_eddy&lev = hyb_plot
    diff_eddy!1   = "rad"
    diff_eddy&rad = radii

  elseif (typ .eq. "novc_mono") then 

  ; Sum of mean + eddy                                                              
    diff_sum      = mono_sum - novc_sum
    diff_sum!0    = "lev"
    diff_sum&lev  = hyb_plot
    diff_sum!1    = "rad"
    diff_sum&rad  = radii

  ; Mean                                                                            
    diff_mean     = mono_arr(8,:,:) - novc_arr(8,:,:)
    diff_mean!0   = "lev"
    diff_mean&lev = hyb_plot
    diff_mean!1   = "rad"
    diff_mean&rad = radii

  ; Eddy                                                                            
    diff_eddy     = mono_arr(9,:,:) - novc_arr(9,:,:)
    diff_eddy!0   = "lev"
    diff_eddy&lev = hyb_plot
    diff_eddy!1   = "rad"
    diff_eddy&rad = radii

  elseif (typ .eq. "novc_m2r") then 

  ; Sum of mean + eddy                                                              
    diff_sum      = m2r_sum - novc_sum
    diff_sum!0    = "lev"
    diff_sum&lev  = hyb_plot
    diff_sum!1    = "rad"
    diff_sum&rad  = radii

  ; Mean                                                                            
    diff_mean     = m2r_arr(8,:,:) - novc_arr(8,:,:)
    diff_mean!0   = "lev"
    diff_mean&lev = hyb_plot
    diff_mean!1   = "rad"
    diff_mean&rad = radii

  ; Eddy                                                                            
    diff_eddy     = m2r_arr(9,:,:) - novc_arr(9,:,:)
    diff_eddy!0   = "lev"
    diff_eddy&lev = hyb_plot
    diff_eddy!1   = "rad"
    diff_eddy&rad = radii

  end if 

; Calculate mean of all four VC phases 
  mean_arr      = dim_avg_n_Wrap(mbgt_arr(0:3,:,:,:), 0)

; Select data to plot 
; "ring", "r2m", "mono", "m2r", "mean", "mean_r", "mean_r2", "mean_m", "mean_m2"
; "ring_diff", "r2m_diff", "mono_diff", "m2r_diff" 
  if (typ .eq. "ring")

    plot_arr = ring_arr
    plot_tt  = "Ringlike phase"

  elseif (typ .eq. "r2m") then

    plot_arr = r2m_arr
    plot_tt  = "Ringlike to monopole transition"

  elseif (typ .eq. "mono") then

    plot_arr = mono_arr
    plot_tt  = "Monopole phase"

  elseif (typ .eq. "m2r") then 

    plot_arr = m2r_arr   
    plot_tt  = "Monopole to ringlike transition"

  elseif (typ .eq. "novc") then 

    plot_arr = novc_arr
    plot_tt  = "No VC"

  elseif (typ .eq. "mean") then 

    plot_arr = mean_arr
    plot_tt  = "Mean of all VC phases"

  elseif (typ .eq. "mean_r") then 

    plot_arr = ring_arr - mean_arr
    plot_tt  = "Ringlike phase - mean"

  elseif (typ .eq. "mean_r2") then 

    plot_arr = r2m_arr - mean_arr
    plot_tt  = "R2M transition - mean"

  elseif (typ .eq. "mean_m") then

    plot_arr = mono_arr - mean_arr
    plot_tt  = "Monopole phase - mean"

  elseif (typ .eq. "mean_m2") then

    plot_arr = m2r_arr - mean_arr
    plot_tt  = "M2R transition - mean"

  elseif (typ .eq. "ring_diff")

    plot_arr = r2m_arr - ring_arr
    plot_tt  = "R2M - Ring"

  elseif (typ .eq. "r2m_diff") then

    plot_arr = mono_arr - r2m_arr
    plot_tt  = "Mono - R2M"

  elseif (typ .eq. "mono_diff") then

    plot_arr = m2r_arr - mono_arr
    plot_tt  = "M2R - Mono"

  elseif (typ .eq. "m2r_diff") then

    plot_arr = ring_arr - m2r_arr
    plot_tt  = "Ring - M2R"

  elseif (typ .eq. "paper0_diff") then 

    plot_arr = mono_arr - ring_arr
    plot_tt  = "Mono - Ring"

  elseif (typ .eq. "paper1_diff") then 

    plot_arr = ring_arr - mono_arr
    plot_tt  = "Ring - Mono"

  elseif (typ .eq. "novc_ring") then

    plot_arr = ring_arr - novc_arr
    plot_tt  = ""

  elseif (typ .eq. "novc_r2m") then

    plot_arr = r2m_arr - novc_arr
    plot_tt  = ""

  elseif (typ .eq. "novc_mono") then

    plot_arr = mono_arr - novc_arr
    plot_tt  = ""

  elseif (typ .eq. "novc_m2r") then

    plot_arr = m2r_arr - novc_arr
    plot_tt  = ""

  elseif (typ .eq. "novc_mean") then 
  
    plot_arr = mean_arr - novc_arr
    plot_tt  = ""    

  end if   

 end if 

; Add metadata to new arrays before plotting 
  if (typ .ne. "comp" .and. typ .ne. "mean_eddy" .and. \
      typ .ne. "rad" .and. typ .ne. "lhs" .and. typ .ne. "paper2_diff" .and. \
      typ .ne. "novc_rad" .and. typ .ne. "novc_lhs") then 
   plot_arr!0    = "diag"
   plot_arr&diag = diag_arr
   plot_arr!1    = "lev"
   plot_arr&lev  = hyb_plot
   plot_arr!2    = "rad"
   plot_arr&rad  = radii
; 25/12/2019 --> plots look completely wrong without this step
  elseif (typ .eq. "novc_rad" .or. typ .eq. "novc_lhs") then 
   novc_arr!0    = "diag"
   novc_arr&diag = diag_arr
   novc_arr!1    = "lev"
   novc_arr&lev  = hyb_plot
   novc_arr!2    = "rad"
   novc_arr&rad  = radii
  end if 

;================================================
; Calculate radius of maximum tangential wind 
;================================================

; Ring-like phase 
  vtan_ring_1d  = ndtooned(ring_arr(13,:,:) )
  dsize_ring    = dimsizes(ring_arr(13,:,:) )
  ring_rmw      = ind_resolve(maxind(vtan_ring_1d),dsize_ring)
  ring_r_rmw    = radii(ring_rmw(0,1) )
  ring_z_rmw    = hyb_plot(ring_rmw(0,0) )

; Ring-like to monopole transition 
  vtan_r2m_1d   = ndtooned(r2m_arr(13,:,:) )
  dsize_r2m     = dimsizes(r2m_arr(13,:,:) )
  r2m_rmw       = ind_resolve(maxind(vtan_r2m_1d),dsize_r2m)
  r2m_r_rmw     = radii(r2m_rmw(0,1) )
  r2m_z_rmw     = hyb_plot(r2m_rmw(0,0) )

; Monopole phase     
  vtan_mono_1d  = ndtooned(mono_arr(13,:,:) )
  dsize_mono    = dimsizes(mono_arr(13,:,:) )
  mono_rmw      = ind_resolve(maxind(vtan_mono_1d),dsize_mono)
  mono_r_rmw    = radii(mono_rmw(0,1) )
  mono_z_rmw    = hyb_plot(mono_rmw(0,0) )

; Monopole to ring-like transition
  vtan_m2r_1d   = ndtooned(m2r_arr(13,:,:) )
  dsize_m2r     = dimsizes(m2r_arr(13,:,:) )
  m2r_rmw       = ind_resolve(maxind(vtan_m2r_1d),dsize_m2r)
  m2r_r_rmw     = radii(m2r_rmw(0,1) )
  m2r_z_rmw     = hyb_plot(m2r_rmw(0,0) )

; Non-VC composite 
  if (typ .eq. "novc_ring" .or. typ .eq. "novc_r2m" .or. typ .eq. "novc_mono" \
      .or. typ .eq. "novc_m2r") then 
   vtan_novc_1d  = ndtooned(novc_arr(13,:,:) )
   dsize_novc    = dimsizes(novc_arr(13,:,:) )
   novc_rmw      = ind_resolve(maxind(vtan_novc_1d),dsize_novc)
   novc_r_rmw    = radii(novc_rmw(0,1) )
   novc_z_rmw    = hyb_plot(novc_rmw(0,0) )
  end if 

;============================================
; Define output file path and plot details 
;============================================

; Create option for alternative resources if we're plotting difference fields  
  if (typ .ne. "mean" .and. typ .ne. "novc_ring" .and. typ .ne. "novc_r2m" .and. \
      typ .ne. "novc_mono" .and. typ .ne. "novc_m2r" .and. typ .ne. "novc_mean" .and. \
      typ .ne. "ring" .and. typ .ne. "r2m" .and. typ .ne. "mono" .and. \
      typ .ne. "m2r" .and. typ .ne. "v" .and. typ .ne. "slp" .and. typ .ne. "novc") then
   opts_alt     = True 
  else
   opts_alt     = False
  end if    

  ; Momentum budget terms 
    opts_mbgt                              = True
    opts_mbgt@cnFillOn                     = True
    opts_mbgt@cnLineLabelInterval          = 2.0
    opts_mbgt@cnLineLabelFontHeightF       = 0.012
    opts_mbgt@cnLineLabelBackgroundColor   = "transparent"
    opts_mbgt@cnLineLabelPlacementMode     = "constant"
    opts_mbgt@cnLinesOn                    = False
    opts_mbgt@cnInfoLabelOn                = False
    opts_mbgt@cnLevelSelectionMode         = "ExplicitLevels"
    opts_mbgt@cnFillPalette                = "BlueDarkRed18"
  ; Select contour levels based upon type of output plot 
    if (typ .eq. "rad" .or. typ .eq. "novc_rad") then 
     opts_mbgt@cnFillPalette               = "BlueDarkRed18"
     opts_mbgt@cnLevels                    = (/-10., -5., -3., -2., -1., -0.5, -0.2, \
                                                 0, 0.2, 0.5, 1., 2., 3., 5., 10./)
     opts_mbgt@cnFillColors                = (/0,2,3,4,5,6,7,-1,-1,10,11,12,13,14,15,17/)  
    elseif (typ .eq. "novc_lhs") then 
      opts_mbgt@cnLevels                   = (/-15.0, -10.0, -5.0, -3.0, -2.0, -1.0, -0.5, -0.2, \
                                                0, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 15.0/)
      opts_mbgt@cnFillColors               = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17/)
    else
     if (res .eq. "1h") then
      opts_mbgt@cnLevels                   = (/-20.0, -15.0, -10.0, -5.0, -3.0, -2.0, -1.0, -0.5, \
                                                0, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 15.0, 20.0/)
      opts_mbgt@cnFillColors               = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17/)
     else
      opts_mbgt@cnLevels                   = (/-25.0, -20.0, -15.0, -10.0, -5.0, -2.0, -1.0, -0.5, \
                                                0, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 25.0/)
      opts_mbgt@cnFillColors               = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17/)
     end if
    end if 

    opts_mbgt@gsnPaperOrientation          = "landscape"
    opts_mbgt@tiMainString                 = ""
    opts_mbgt@tiMainFontHeightF            = 0.0125
    opts_mbgt@gsnLeftString                = ""
    opts_mbgt@gsnRightString               = ""
    opts_mbgt@gsnMaximize                  = False
    opts_mbgt@lbLabelBarOn                 = False
    opts_mbgt@gsnDraw                      = False
    opts_mbgt@gsnFrame                     = False

  ; Additional plotting resources 
    opts_mbgt@tiYAxisString                = "Height (km)"
    opts_mbgt@tiXAxisString                = "Radial distance from cyclone centre (km)"

  ; Turn off tickmarks on top x-axis and right y-axis
    opts_mbgt@tmXTOn                       = "False"
    opts_mbgt@tmYROn                       = "False"

  ; Do not maximise the individual panels 
    opts_mbgt@gsnMaximize                  = False  
    opts_mbgt@gsnAddCyclic                 = False

  ; Set x-axis limits 
    opts_mbgt@tmXBMode                     = "Explicit"
    if (zoom .eq. 1 .or. typ .eq. "rad" .or. typ .eq. "lhs" .or. \
        typ .eq. "novc_rad" .or. typ .eq. "novc_lhs") then 
     opts_mbgt@tmXBValues                  = (/0.0,0.25,0.5/)
     opts_mbgt@tmXBLabels                  = (/"0","25","50"/)
     opts_mbgt@trXMaxF                     = 0.5
    else
     opts_mbgt@tmXBValues                  = (/0.0,0.25,0.5,0.75,1.0/)
     opts_mbgt@tmXBLabels                  = (/"0","25","50","75","100"/)
     opts_mbgt@trXMaxF                     = 1.0
    end if 

  ; Linearise the y-axis, and set min/max values
    opts_mbgt@gsnYAxisIrregular2Linear     = True
    opts_mbgt@trYMinF                      = 0.0 
    opts_mbgt@tmYLMode                     = "Explicit"

    if (zoom .eq. 1) then 
     opts_mbgt@trYMaxF                     = 8.0 
     opts_mbgt@tmYLValues                  = fspan(2.0, 8.0, 4)
     opts_mbgt@tmYLLabels                  = (/"2.0","4.0","6.0","8.0"/)
    elseif (typ .eq. "rad" .or. typ .eq. "lhs" .or. \
            typ .eq. "novc_rad" .or. typ .eq. "novc_lhs") then 
     opts_mbgt@trYMaxF                     = 5.0
     opts_mbgt@tmYLValues                  = fspan(1.0, 5.0, 5)
     opts_mbgt@tmYLLabels                  = (/"1.0","2.0","3.0","4.0", "5.0"/)
    else
     opts_mbgt@trYMaxF                     = 17.6
     opts_mbgt@tmYLValues                  = fspan(2.0, 18.0, 9)
     opts_mbgt@tmYLLabels                  = (/"2.0","4.0","6.0","8.0","10.0",\
                                               "12.0","14.0","16.0","18.0"/)
    end if 

  ; Resources to overlay radial wind (following M18)
    opts_vrad                              = True
    opts_vrad@cnFillOn                     = False
    opts_vrad@cnLineColor                  = "black"
    opts_vrad@cnLevelSelectionMode         = "ExplicitLevels"
    opts_vrad@cnLevels                     = (/-1.2, 1.2/)
    opts_vrad@gsnContourNegLineDashPattern = 16
    opts_vrad@gsnContourPosLineDashPattern = 0
    opts_vrad@cnLineThicknessF             = 2.5

    opts_vrad@cnInfoLabelOn                = False
    opts_vrad@gsnDraw                      = False
    opts_vrad@gsnFrame                     = False
    opts_vrad@cnLineLabelsOn               = False
    opts_vrad@tiMainString                 = ""
    opts_vrad@gsnLeftString                = ""
    opts_vrad@gsnRightString               = ""

  ; Resources to overlay absolute angular momentum 
    opts_aam                               = True 
    opts_aam@cnFillOn                      = False
    opts_aam@cnLineColor                   = "black"
    opts_aam@cnLevelSelectionMode          = "ExplicitLevels"
    opts_aam@cnLevels                      = (/1.0, 1.5, 2.0/)
    opts_aam@gsnContourNegLineDashPattern  = 16
    opts_aam@gsnContourPosLineDashPattern  = 0
    opts_aam@cnLineThicknessF              = 4.0

    opts_aam@cnInfoLabelOn                 = False
    opts_aam@gsnDraw                       = False
    opts_aam@gsnFrame                      = False
    opts_aam@cnLineLabelsOn                = False
    opts_aam@tiMainString                  = ""
    opts_aam@gsnLeftString                 = ""
    opts_aam@gsnRightString                = ""

  ; Resources to overlay vertical velocity (following M18) 
    opts_vvel                              = True
    opts_vvel@cnFillOn                     = False
    opts_vvel@cnLineColor                  = "yellow"
    opts_vvel@cnLevelSelectionMode         = "ExplicitLevels"
    opts_vvel@cnLevels                     = (/0.5/)
    opts_vvel@cnLineThicknessF             = 5.0

    opts_vvel@cnInfoLabelOn                = False
    opts_vvel@gsnDraw                      = False
    opts_vvel@gsnFrame                     = False
    opts_vvel@cnLineLabelsOn               = False
    opts_vvel@tiMainString                 = ""
    opts_vvel@gsnLeftString                = ""
    opts_vvel@gsnRightString               = ""

  ; Resources to overlay zero line
    opts_zero                             = True
    opts_zero@cnFillOn                    = False
    opts_zero@cnLineColor                 = "grey50"
    opts_zero@cnLineDashPattern           = 3
    opts_zero@cnLevelSelectionMode        = "ExplicitLevels"
    opts_zero@cnLevels                    = 0.0
    opts_zero@cnLineThicknessF            = 3.0
    opts_zero@cnInfoLabelOn               = False
    opts_zero@gsnDraw                     = False
    opts_zero@gsnFrame                    = False 
    opts_zero@gsnContourZeroLineThicknessF = 2.0
    opts_zero@cnLineLabelsOn              = False
    opts_zero@tiMainString                = ""
    opts_zero@gsnLeftString               = ""
    opts_zero@gsnRightString              = "" 

  ; Radial wind 
    opts_vrad2                             = opts_mbgt
    delete(opts_vrad2@cnLevels)
    delete(opts_vrad2@cnFillColors)
    opts_vrad2@cnFillPalette               = "BlueDarkRed18"
    opts_vrad2@cnLevels                    = (/-10., -5., -3., -2., -1., -0.5, -0.2, \
                                                0, 0.2, 0.5, 1., 2., 3., 5., 10./)
    opts_vrad2@cnFillColors                = (/0,2,3,4,5,6,7,-1,-1,10,11,12,13,14,15,17/)
    opts_vrad2@lbLabelBarOn                = True
    opts_vrad2@lbBoxEndCapStyle            = "TriangleBothEnds"
    opts_vrad2@lbOrientation               = "Vertical"
    if (paper .eq. 0) then 
     opts_vrad2@vpXF                       = 0.05
     opts_vrad2@vpYF                       = 0.95
     opts_vrad2@vpWidthF                   = 0.16
     opts_vrad2@vpHeightF                  = 0.16
    end if 

  ; Tangential wind 
    opts_vtan2                             = opts_mbgt
    if (.not. opts_alt) then 
     delete(opts_vtan2@cnLevels)
     delete(opts_vtan2@cnFillColors)
     opts_vtan2@cnFillPalette               = "radar_new"
     opts_vtan2@cnLevels                    = (/5.0, 10.0, 15.0, 20.0,\
                                                25.0, 30.0, 35.0, 40.0, \
                                                45.0, 50.0, 55.0, 60.0, \
                                                65.0/)
     opts_vtan2@cnFillColors                = (/-1,0,1,2,3,4,5,6,7,8,9,10,11,12/)
    end if 
    opts_vtan2@lbLabelBarOn                = True
    opts_vtan2@lbBoxEndCapStyle            = "TriangleBothEnds"
    opts_vtan2@lbOrientation               = "Vertical"
    if (paper .eq. 0) then 
     opts_vtan2@vpXF                       = 0.30
     opts_vtan2@vpYF                       = 0.95
     opts_vtan2@vpWidthF                   = 0.16
     opts_vtan2@vpHeightF                  = 0.16
    else
     opts_vtan2@vpXF                       = 0.37
     opts_vtan2@vpYF                       = 0.95
     opts_vtan2@vpWidthF                   = 0.20
     opts_vtan2@vpHeightF                  = 0.20
    end if 

  ; Vertical velocity 
    opts_vvel2                             = opts_mbgt
    if (.not. opts_alt) then 
     delete(opts_vvel2@cnLevels)
     delete(opts_vvel2@cnFillColors)
     opts_vvel2@cnFillPalette               = "n11"
     opts_vvel2@cnLevels                    = (/0.1, 0.2, 0.3, 0.4, 0.5, 0.6, \
                                                0.7, 0.8, 1.0, 1.2, 1.5, 2.0/)
     opts_vvel2@cnFillColors                = (/-1,0,1,2,3,4,5,6,7,8,9,10,11/)
    end if 
    opts_vvel2@lbLabelBarOn                = True
    opts_vvel2@lbBoxEndCapStyle            = "TriangleBothEnds"
    opts_vvel2@lbOrientation               = "Vertical"
    if (paper .eq. 0) then 
     opts_vvel2@vpXF                       = 0.55
     opts_vvel2@vpYF                       = 0.95
     opts_vvel2@vpWidthF                   = 0.16
     opts_vvel2@vpHeightF                  = 0.16
    else
     opts_vvel2@vpXF                       = 0.69
     opts_vvel2@vpYF                       = 0.95
     opts_vvel2@vpWidthF                   = 0.20
     opts_vvel2@vpHeightF                  = 0.20
    end if 

  ; Absolute angular momentum
    opts_aam2                               = opts_vtan2
    if (.not. opts_alt) then  
     delete(opts_aam2@cnLevels)
     delete(opts_aam2@cnFillColors)
     opts_aam2@cnLevels                     = (/0.6, 0.8, 1.0, 1.2, 1.5, \
                                               2.0, 2.5, 3.0, 3.5, \
                                               4.0, 4.5, 5.0/)
     opts_aam2@cnFillColors                 = (/-1,0,1,2,3,4,5,6,7,8,9,10,11/)
    end if 
    opts_aam2@lbLabelBarOn                  = True
    opts_aam2@lbBoxEndCapStyle              = "TriangleBothEnds"
    opts_aam2@lbOrientation                 = "Vertical"
    if (paper .eq. 0) then
     opts_aam2@vpXF                         = 0.80
     opts_aam2@vpYF                         = 0.95
     opts_aam2@vpWidthF                     = 0.16
     opts_aam2@vpHeightF                    = 0.16
    else
     opts_aam2@vpXF                         = 0.05
     opts_aam2@vpYF                         = 0.67
     opts_aam2@vpWidthF                     = 0.20
     opts_aam2@vpHeightF                    = 0.20
    end if 

;=============================================================== 
; Resources to overlay panel labels without using 'gsn_panel'
;=============================================================== 

  ; Create option for alternative resources if we're plotting difference fields
  ; 10/10/2019 --> add IF statement below into script 
    if (typ .eq. "mean" .or. typ .eq. "ring" .or. \
        typ .eq. "r2m" .or. typ .eq. "mono" .or. typ .eq. "m2r" .or. \
        typ .eq. "v" .or. typ .eq. "slp" .or. typ .eq. "novc") then

     if (paper .eq. 0) then 

      txid                        = new(12,graphic)
      amid                        = new(12,graphic)

    ; Label text resources                                                                   
      txres                       = True
      txres@txPerimOn             = True
      txres@txFontHeightF         = 0.010
      txres@txBackgroundFillColor = "White"

    ; Label position resources
    ; 10/10/2019 --> (amParallelPosF = -0.5/0.5 left/right edge; amOrthogonalPosF = -0.5 top edge) 
      amres                       = True
      amres@amParallelPosF        = 0.5 
      amres@amOrthogonalPosF      = -0.5
      amres@amJust                = "TopRight"

    ; Panel strings 
      panel_strings               = (/"a) Radial wind", "b) Tangential wind", \
                                      "c) Vertical velocity", "d) AAM", \
                                      "e) Mean h. advection", "f) Mean v. advection", \
                                      "g) Eddy h. advection", "h) Eddy v. advection", \
                                      "i) Total mean", "j) Total eddy", \
                                      "k) Sum of RHS terms", "l) Local tendency (LHS)"/)

   ; paper = 1 (8 rather than 12 panels)
     else

      txid                        = new(8,graphic)
      amid                        = new(8,graphic)

    ; Label text resources
      txres                       = True
      txres@txPerimOn             = True
      txres@txFontHeightF         = 0.0125
      txres@txBackgroundFillColor = "White"

    ; Label position resources
    ; 10/10/2019 --> (amParallelPosF = -0.5/0.5 left/right edge; amOrthogonalPosF = -0.5 top edge)
      amres                       = True
      amres@amParallelPosF        = 0.5
      amres@amOrthogonalPosF      = -0.5
      amres@amJust                = "TopRight"

    ; Panel strings
      panel_strings               = (/"a) Radial wind", "b) Tangential wind", \
                                      "c) Vertical velocity", "d) AAM", \
                                      "e) Mean tendency", "f) Eddy tendency", \
                                      "g) Total RHS", "h) Total LHS"/)

     end if 
     
  ; MEAN vs EDDY plots 
    elseif (typ .eq. "mean_eddy") then 

     txid                        = new(6,graphic)
     amid                        = new(6,graphic)

   ; Label text resources
     txres                       = True
     txres@txPerimOn             = True
     txres@txFontHeightF         = 0.0150
     txres@txBackgroundFillColor = "White"

   ; Label position resources 
   ; 10/10/2019 --> (amParallelPosF = -0.5/0.5 left/right edge; amOrthogonalPosF = -0.5 top edge) 
     amres                       = True
     amres@amParallelPosF        = 0.5
     amres@amOrthogonalPosF      = -0.5
     amres@amJust                = "TopRight"

   ; Panel strings 
     panel_strings               = (/"a) Mean", "b) Eddy", "c) Mean + eddy", \
                                     "d) Mean", "e) Eddy", "f) Mean + eddy"/)

  ; Comparison plots (mean vs VC phase) 
    elseif (typ .eq. "comp") then 

     txid                        = new(16,graphic)
     amid                        = new(16,graphic)

   ; Label text resources
     txres                       = True
     txres@txPerimOn             = True
     txres@txFontHeightF         = 0.010
     txres@txBackgroundFillColor = "White"

   ; Label position resources
   ; 10/10/2019 --> (amParallelPosF = -0.5/0.5 left/right edge; amOrthogonalPosF = -0.5 top edge)
     amres                       = True
     amres@amParallelPosF        = 0.5
     amres@amOrthogonalPosF      = -0.5
     amres@amJust                = "TopRight"

     panel_strings               = (/"a)", "b)", "c)", "d)", \
                                     "e)", "f)", "g)", "h)", \
                                     "i)", "j)", "k)", "l)", \
                                     "m)", "n)", "o)", "p)"/)

  ; Plot radial wind overlaid with AAM contours 
    elseif (typ .eq. "rad" .or. typ .eq. "lhs") then 

     txid                        = new(4,graphic)
     amid                        = new(4,graphic)

   ; Label text resources 
     txres                       = True
     txres@txPerimOn             = True
     txres@txFontHeightF         = 0.015
     txres@txBackgroundFillColor = "White"

   ; Label position resources 
   ; (amParallelPosF = -0.5/0.5 left/right edge; amOrthogonalPosF = -0.5 top edge)
     amres                       = True
     amres@amParallelPosF        = 0.5
     amres@amOrthogonalPosF      = -0.5
     amres@amJust                = "TopRight"
     panel_strings               = (/"a) Ring", "b) Ring to mono", "c) Mono", "d) Mono to ring"/)

  ; Radial wind + AAM or local tangential wind tendency (Non-VC)
    elseif (typ .eq. "novc_rad" .or. typ .eq. "novc_lhs") then 

     if (col .eq. 2) then 
      txid                       = new(2,graphic)
      amid                       = new(2,graphic)
     else
      txid                       = new(3,graphic)
      amid                       = new(3,graphic)
     end if 

   ; Label text resources 
     txres                       = True
     txres@txPerimOn             = True
     txres@txFontHeightF         = 0.015
     txres@txBackgroundFillColor = "White"

   ; Label position resources   
   ; (amParallelPosF = -0.5/0.5 left/right edge; amOrthogonalPosF = -0.5 top edge) 
     amres                       = True
     amres@amParallelPosF        = 0.5
     amres@amOrthogonalPosF      = -0.5
     amres@amJust                = "TopRight"

   ; Label for second panel
     if (p0 .eq. "ring") then 
      bstr = "b) Ring"
     elseif (p0 .eq. "r2m") then 
      bstr = "b) Ring to mono"
     elseif  (p0 .eq. "mono") then 
      bstr = "b) Mono"
     else  
      bstr = "b) Mono to ring"
     end if 

   ; Label for third panel (only if required)
     if (col .ne. 2) then 
      if (p1 .eq. "ring") then
       cstr = "c) Ring"
      elseif (p1 .eq. "r2m") then
       cstr = "c) Ring to mono"
      elseif (p1 .eq. "mono") then
       cstr = "c) Mono"
      elseif (p1 .eq. "mr2") then 
       cstr = "c) Mono to ring"
      end if
     end if 

     if (col .eq. 2) then 
      panel_strings              = (/"a) Non-VC", bstr/)
     else
      panel_strings              = (/"a) Non-VC", bstr, cstr/)
     end if 

  ; Difference plots (combined mean + eddy terms)
    elseif (typ .eq. "paper2_diff" .or. typ .eq. "novc_mean" .or. typ .eq. "novc_ring" .or. \
            typ .eq. "novc_r2m" .or. typ .eq. "novc_mono" .or. typ .eq. "novc_m2r") then 

     txid                        = new(3,graphic)
     amid                        = new(3,graphic)

   ; Label text resources
     txres                       = True
     txres@txPerimOn             = True
     txres@txFontHeightF         = 0.0150
     txres@txBackgroundFillColor = "White"

   ; Label position resources
     amres                       = True
     amres@amParallelPosF        = 0.5
     amres@amOrthogonalPosF      = -0.5
     amres@amJust                = "TopRight"

     panel_strings               = (/"a) Mean", "b) Eddy", "c) Mean + eddy"/)

  ; Difference plots (VC vs non-VC)
    else

     txid                        = new(4,graphic)
     amid                        = new(4,graphic)

   ; Label text resources
     txres                       = True
     txres@txPerimOn             = True
     if (typ .eq. "novc_mean" .or. typ .eq. "novc_ring" .or. typ .eq. "novc_r2m" .or. \
         typ .eq. "novc_mono" .or. typ .eq. "novc_m2r"  .or. \ 
         typ .eq. "paper1_diff" .or. typ .eq. "paper0_diff") then
      txres@txFontHeightF        = 0.0125
     else
      txres@txFontHeightF        = 0.010
     end if 
     txres@txBackgroundFillColor = "White"

   ; Label position resources
   ; 10/10/2019 --> (amParallelPosF = -0.5/0.5 left/right edge; amOrthogonalPosF = -0.5 top edge) 
     amres                       = True
     amres@amParallelPosF        = 0.5
     amres@amOrthogonalPosF      = -0.5
     amres@amJust                = "TopRight"

     panel_strings               = (/"a) V_mz", "b) V_mv", "c) V_ez", "d) V_ev"/)

    end if 

;=================================
; Plot the momentum budget data
;=================================

  ; Calculating the terms from Montgomery et al. (2018) - their Eq.(4) 
  ; LHS (1) : d{v}/dt          - Mean tangential momentum tendency
  ; RHS (1) : -{u} * {f+vort}  - Mean radial influx of absolute vertical vorticity
  ; RHS (2) : -{w} * (d{v}/dz) - Mean vertical advection of mean tangential momentum
  ; RHS (3) : -{u' * vort'}    - Eddy radial vorticity flux
  ; RHS (4) : -{w' * (dv'/dz)} - Vertical advection of eddy tangential momentum
  ; RHS (5) : {D_v}            - Combined diffusive and PBL tendency 

  ; Change output file string depending on whether we use 5-min or 1-h data 
    if (res .eq. "1h") then 
     out_str = "mbgt"
    else
     out_str = "mbgt5"
    end if 

  ; Output file location and type  
    if (zoom .eq. 1) then 
      if (typ .eq. "v" .or. typ .eq. "slp") then 
       output = "$nep/nepartak/images/cyl_coords/"+out_str+"_tan_cxz_diff_group"+grp+\
                "_sc"+sc+"_"+z0+"km_"+w0+"w_"+typ+thr+"_zoom"
      elseif (typ .eq. "novc") then 
       output = "$nep/nepartak/images/cyl_coords/"+out_str+"_tan_cxz_diff_group"+grp+\
                "_sc"+sc+"_"+z0+"km_"+w0+"w_"+typ+"_"+ts+"_"+tf+"_zoom"
      else
       output = "$nep/nepartak/images/cyl_coords/"+out_str+"_tan_cxz_diff_group"+grp+\
                "_sc"+sc+"_"+z0+"km_"+w0+"w_"+typ+"_zoom"
      end if 
    else
      if (typ .eq. "v" .or. typ .eq. "slp") then
       output = "$nep/nepartak/images/cyl_coords/"+out_str+"_tan_cxz_diff_group"+grp+\
                "_sc"+sc+"_"+z0+"km_"+w0+"w_"+typ+thr
      elseif (typ .eq. "novc" .or. typ .eq. "novc_mean" .or. typ .eq. "novc_ring" .or. \
              typ .eq. "novc_r2m" .or. typ .eq. "novc_mono" .or. typ .eq. "novc_m2r") then 
       output = "$nep/nepartak/images/cyl_coords/"+out_str+"_tan_cxz_diff_group"+grp+\
                "_sc"+sc+"_"+z0+"km_"+w0+"w_"+typ+"_"+ts+"_"+tf
      elseif (typ .eq. "novc_rad" .or. typ .eq. "novc_lhs") then 
       output = "$nep/nepartak/images/cyl_coords/"+out_str+"_tan_cxz_diff_group"+grp+\
                "_sc"+sc+"_"+z0+"km_"+w0+"w_"+typ+"_"+p0+"_"+p1
      else
       output = "$nep/nepartak/images/cyl_coords/"+out_str+"_tan_cxz_diff_group"+grp+\
                "_sc"+sc+"_"+z0+"km_"+w0+"w_"+typ
      end if
    end if 
    wks    = gsn_open_wks(opt,output)

  ; Where are the data for all these terms located? 
  ; (0) LHS_1, (1) RHS_1, (2) RHS_2, (3) RHS_3, (4) RHS_4, (5) RHS_5A, (6) RHS_5B, (7) RHS_T
  ; (8) MEAN, (9) EDDY, (10) VFLUX, (11) VADV, (12) VRAD, (13) VTAN, (14) VVEL, (15) AAM

;============================================================================
; Tidy this part up, but customise the number/size/position of panels here 
;============================================================================

  if (typ .eq. "rad" .or. typ .eq. "lhs" .or. \
      typ .eq. "paper1_diff" .or. typ .eq. "paper0_diff") then 
   xpos_plot = (/0.19, 0.57, 0.19, 0.57/)
   ypos_plot = (/0.90, 0.90, 0.55, 0.55/)
   xsize = 0.27
   ysize = 0.27
  elseif (typ .eq. "paper2_diff" .or. typ .eq. "novc_ring" .or. \
          typ .eq. "novc_r2m" .or. typ .eq. "novc_mono" .or. typ .eq. "novc_m2r") then 
   xpos_plot = (/0.07, 0.39, 0.71/)
   ypos_plot = (/0.97, 0.97, 0.97/)
   xsize     = 0.24
   ysize     = 0.24
  elseif (typ .eq. "novc_rad" .or. typ .eq. "novc_lhs") then 
   if (col .eq. 2) then 
    xpos_plot= (/0.22, 0.54/)
    ypos_plot= (/0.97, 0.97/)
   else
   xpos_plot = (/0.07, 0.39, 0.71/)
   ypos_plot = (/0.97, 0.97, 0.97/)
   end if 
   xsize     = 0.24
   ysize     = 0.24
  end if 

;=======================
; Start plotting data
;=======================

; Standard momentum budget plots for single VC phases or similar
  if (typ .eq. "mean" .or. typ .eq. "novc" .or. \
      typ .eq. "ring" .or. typ .eq. "r2m" .or. typ .eq. "mono" .or. \
      typ .eq. "m2r" .or. typ .eq. "v" .or. typ .eq. "slp") then

   if (paper .eq. 0) then 
     print("Producing 12-panel plot...")

  ; Before plotting, calculate total RHS tendency and add metadata
    rhs_tot     = plot_arr(1,:,:) + plot_arr(2,:,:) + plot_arr(3,:,:) + \
                  plot_arr(4,:,:) + plot_arr(6,:,:)
    rhs_tot!0   = "lev"
    rhs_tot&lev = hyb_plot
    rhs_tot!1   = "rad"
    rhs_tot&rad = radii

  ; Radial wind (filled contours) 
    panel_vrad  = gsn_csm_contour(wks,plot_arr(12,:,:),opts_vrad2) ; [u]                            
    panel_zero  = gsn_csm_contour(wks,plot_arr(12,:,:),opts_zero)  ; Zero line contour              
    zero0       = ColorNegDashZeroPosContour(panel_zero,"transparent","black","transparent")
    overlay(panel_vrad, zero0)
    txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
    amid(0)     = gsn_add_annotation(panel_vrad, txid(0), amres)

  ; Tangential wind (filled contours)                                                          
    panel_vtan  = gsn_csm_contour(wks,plot_arr(13,:,:),opts_vtan2) ; [v]                            
    txid(1)     = gsn_create_text(wks, panel_strings(1), txres)
    amid(1)     = gsn_add_annotation(panel_vtan, txid(1), amres)

  ; Vertical velocity (filled contours)                                                        
    panel_vvel  = gsn_csm_contour(wks,plot_arr(14,:,:),opts_vvel2) ; [w]                            
    txid(2)     = gsn_create_text(wks, panel_strings(2), txres)
    amid(2)     = gsn_add_annotation(panel_vvel, txid(2), amres)

  ; Absolute angular momentum (filled contours)                                                
    panel_aam   = gsn_csm_contour(wks,plot_arr(15,:,:),opts_aam2)  ; [AAM]
    txid(3)     = gsn_create_text(wks, panel_strings(3), txres)
    amid(3)     = gsn_add_annotation(panel_aam, txid(3), amres)

  ; Mean radial influx of absolute vertical vorticity [V_mzeta]
    opts_mbgt@vpXF       = 0.05
    opts_mbgt@vpYF       = 0.74
    opts_mbgt@vpWidthF   = 0.19
    opts_mbgt@vpHeightF  = 0.19
    panel0      = gsn_csm_contour(wks,plot_arr(1,:,:),opts_mbgt) ; [V_mzeta]                       
    cont0       = gsn_csm_contour(wks,plot_arr(1,:,:),opts_zero) ; Draw contours                   

    panel0_ovr  = gsn_csm_contour(wks,plot_arr(12,:,:),opts_vrad) ; [u]                             
    panel0_ovr2 = gsn_csm_contour(wks,plot_arr(14,:,:),opts_vvel) ; [w]                             

    overlay(panel0, cont0)                                   ; Overlay contours                
    overlay(panel0, panel0_ovr)
    overlay(panel0, panel0_ovr2)
    txid(4)     = gsn_create_text(wks, panel_strings(4), txres)
    amid(4)     = gsn_add_annotation(panel0, txid(4), amres)
    delete([/cont0/])

  ; Mean vertical advection of mean tangential momentum [V_mv]                                 
    opts_mbgt@vpXF       = 0.30
    opts_mbgt@vpYF       = 0.74
    opts_mbgt@vpWidthF   = 0.19
    opts_mbgt@vpHeightF  = 0.19
    panel1      = gsn_csm_contour(wks,plot_arr(2,:,:),opts_mbgt) ; [V_mv]                          
    cont0       = gsn_csm_contour(wks,plot_arr(2,:,:),opts_zero) ; Draw contours                   

    panel1_ovr  = gsn_csm_contour(wks,plot_arr(12,:,:),opts_vrad) ; [U]                             
    panel1_ovr2 = gsn_csm_contour(wks,plot_arr(14,:,:),opts_vvel) ; [w]                             

    overlay(panel1, cont0)                                   ; Overlay contours                
    overlay(panel1, panel1_ovr)
    overlay(panel1, panel1_ovr2)
    txid(5)     = gsn_create_text(wks, panel_strings(5), txres)
    amid(5)     = gsn_add_annotation(panel1, txid(5), amres)
    delete([/cont0/])

  ; Eddy radial vorticity flux [V_ezeta]                                                       
    opts_mbgt@vpXF       = 0.55
    opts_mbgt@vpYF       = 0.74
    opts_mbgt@vpWidthF   = 0.19
    opts_mbgt@vpHeightF  = 0.19
    panel2      = gsn_csm_contour(wks,plot_arr(3,:,:),opts_mbgt) ; [V_ezeta]                       
    cont0       = gsn_csm_contour(wks,plot_arr(3,:,:),opts_zero) ; Draw contours                   

    panel2_ovr  = gsn_csm_contour(wks,plot_arr(12,:,:),opts_vrad) ; [U]                             
    panel2_ovr2 = gsn_csm_contour(wks,plot_arr(14,:,:),opts_vvel) ; [w]                             

    overlay(panel2, cont0)                                   ; Overlay contours                
    overlay(panel2, panel2_ovr)
    overlay(panel2, panel2_ovr2)
    txid(6)     = gsn_create_text(wks, panel_strings(6), txres)
    amid(6)     = gsn_add_annotation(panel2, txid(6), amres)
    delete([/cont0/])

  ; Vertical advection of eddy tangential momentum [V_ev]                                      
    opts_mbgt@vpXF       = 0.80
    opts_mbgt@vpYF       = 0.74
    opts_mbgt@vpWidthF   = 0.19
    opts_mbgt@vpHeightF  = 0.19
    panel3      = gsn_csm_contour(wks,plot_arr(4,:,:),opts_mbgt) ; [V_ev]                          
    cont0       = gsn_csm_contour(wks,plot_arr(4,:,:),opts_zero) ; Draw contours                   

    panel3_ovr  = gsn_csm_contour(wks,plot_arr(12,:,:),opts_vrad) ; [U]                             
    panel3_ovr2 = gsn_csm_contour(wks,plot_arr(14,:,:),opts_vvel) ; [w]                             

    overlay(panel3, cont0)                                   ; Overlay contours                
    overlay(panel3, panel3_ovr)
    overlay(panel3, panel3_ovr2)
    txid(7)     = gsn_create_text(wks, panel_strings(7), txres)
    amid(7)     = gsn_add_annotation(panel3, txid(7), amres)
    delete([/cont0/])

  ; Vertical diffusive tendency of tangential momentum [V_dz]                                  
  ; 14/12/2019 --> combined mean terms 
    opts_mbgt@vpXF       = 0.05 ; 0.19
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.19
    opts_mbgt@vpHeightF  = 0.19

    panel5      = gsn_csm_contour(wks,plot_arr(8,:,:),opts_mbgt)   ; [V_dz]
    cont0       = gsn_csm_contour(wks,plot_arr(8,:,:),opts_zero)   ; Draw contours
    panel5_ovr  = gsn_csm_contour(wks,plot_arr(12,:,:),opts_vrad)  ; [U]
    panel5_ovr2 = gsn_csm_contour(wks,plot_arr(14,:,:),opts_vvel)  ; [w]

    overlay(panel5, cont0)                                         ; Overlay contours               
    overlay(panel5, panel5_ovr)
    overlay(panel5, panel5_ovr2)
    txid(8)     = gsn_create_text(wks, panel_strings(8), txres)
    amid(8)     = gsn_add_annotation(panel5, txid(8), amres)
    delete([/cont0/])

  ; Combined eddy terms 
    opts_mbgt@vpXF       = 0.30
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.19
    opts_mbgt@vpHeightF  = 0.19

    panel6      = gsn_csm_contour(wks,plot_arr(9,:,:),opts_mbgt)   ; [V_dz]                 
    cont0       = gsn_csm_contour(wks,plot_arr(9,:,:),opts_zero)   ; Draw contours          
    panel6_ovr  = gsn_csm_contour(wks,plot_arr(12,:,:),opts_vrad)  ; [U]                    
    panel6_ovr2 = gsn_csm_contour(wks,plot_arr(14,:,:),opts_vvel)  ; [w] 

    overlay(panel6, cont0)                                         ; Overlay contours   
    overlay(panel6, panel6_ovr)
    overlay(panel6, panel6_ovr2)
    txid(9)     = gsn_create_text(wks, panel_strings(9), txres)
    amid(9)     = gsn_add_annotation(panel6, txid(9), amres)
    delete([/cont0/])

  ; Total RHS [All V tend]                                                                    
    opts_mbgt@vpXF       = 0.55 ; 0.43
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.19
    opts_mbgt@vpHeightF  = 0.19

    panel7      = gsn_csm_contour(wks,rhs_tot(:,:),opts_mbgt)    ; [All V tend]
    cont0       = gsn_csm_contour(wks,rhs_tot(:,:),opts_zero)    ; Draw contours
    panel7_ovr  = gsn_csm_contour(wks,plot_arr(12,:,:),opts_vrad)   ; [U]
    panel7_ovr2 = gsn_csm_contour(wks,plot_arr(14,:,:),opts_vvel)   ; [w]

    overlay(panel7, cont0)                                   ; Overlay contours                
    overlay(panel7, panel7_ovr)
    overlay(panel7, panel7_ovr2)
    txid(10)    = gsn_create_text(wks, panel_strings(10), txres)
    amid(10)    = gsn_add_annotation(panel7, txid(10), amres)
    delete([/cont0/])

  ; Local tendency of the mean tangential wind [V_t]                                           
    opts_mbgt@vpXF       = 0.80 ; 0.67
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.19
    opts_mbgt@vpHeightF  = 0.19

  ; Multiply to get correct units 
;    if (res .eq. "5min" .and. sc .eq. 6) then 
;     plot_arr(0,:,:) = plot_arr(0,:,:) * 12
    if (res .eq. "5min" .and. sc .eq. 4 .and. typ .eq. "mono") then 
     plot_arr(0,:,:) = plot_arr(0,:,:) * 12
;    elseif (res .eq. "5min" .and. sc .eq. 5 .and. typ .ne. "ring") then 
;     plot_arr(0,:,:) = plot_arr(0,:,:) * 12
    else
     plot_arr(0,:,:) = plot_arr(0,:,:)
    end if 

    panel8      = gsn_csm_contour(wks,plot_arr(0,:,:),opts_mbgt) ; [V_t]                           
    cont0       = gsn_csm_contour(wks,plot_arr(0,:,:),opts_zero) ; Draw contours                   

    panel8_ovr  = gsn_csm_contour(wks,plot_arr(12,:,:),opts_vrad) ; [U]
    panel8_ovr2 = gsn_csm_contour(wks,plot_arr(14,:,:),opts_vvel) ; [w]

    overlay(panel8, cont0)                                   ; Overlay contours                
    overlay(panel8, panel8_ovr)
    overlay(panel8, panel8_ovr2)
    txid(11)    = gsn_create_text(wks, panel_strings(11), txres)
    amid(11)    = gsn_add_annotation(panel8, txid(11), amres)
    delete([/cont0/])

 ; Same as above, but with 'paper=1' (8 rather than 12 panels)
   elseif (paper .eq. 1) then
    print("Producing 8-panel plot...")

  ; Calculate total RHS tendency 
    rhs_tot     = plot_arr(1,:,:) + plot_arr(2,:,:) + plot_arr(3,:,:) + \
                  plot_arr(4,:,:) + plot_arr(6,:,:)
    rhs_tot!0   = "lev"
    rhs_tot&lev = hyb_plot
    rhs_tot!1   = "rad"
    rhs_tot&rad = radii

    if (row .eq. 3) then 

     xpos_plot = (/0.05, 0.37, 0.69, 0.05, 0.37, 0.69, 0.26, 0.55/)
     ypos_plot = (/0.95, 0.95, 0.95, 0.67, 0.67, 0.67, 0.39, 0.39/)
     xsize     = 0.20
     ysize     = 0.20

    elseif (row .eq. 2) then 

     xpos_plot = (/0.05, 0.30, 0.55, 0.80, 0.05, 0.30, 0.55, 0.80/)
     ypos_plot = (/0.95, 0.95, 0.95, 0.95, 0.72, 0.72, 0.72, 0.72/)
     xsize     = 0.16
     ysize     = 0.16

    end if 

  ; TOP ROW 

  ; Radial wind (filled contours)
    opts_vrad2@vpXF                       = xpos_plot(0)
    opts_vrad2@vpYF                       = ypos_plot(0)
    opts_vrad2@vpWidthF                   = xsize
    opts_vrad2@vpHeightF                  = ysize
    panel_vrad  = gsn_csm_contour(wks,plot_arr(12,:,:),opts_vrad2) ; [u]                            
    panel_zero  = gsn_csm_contour(wks,plot_arr(12,:,:),opts_zero)  ; Zero line contour              
    zero0       = ColorNegDashZeroPosContour(panel_zero,"transparent","black","transparent")
    overlay(panel_vrad, zero0)
    txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
    amid(0)     = gsn_add_annotation(panel_vrad, txid(0), amres)

  ; Tangential wind (filled contours)
    opts_vtan2@vpXF                       = xpos_plot(1)
    opts_vtan2@vpYF                       = ypos_plot(1)
    opts_vtan2@vpWidthF                   = xsize
    opts_vtan2@vpHeightF                  = ysize
    panel_vtan  = gsn_csm_contour(wks,plot_arr(13,:,:),opts_vtan2) ; [v]                            
    txid(1)     = gsn_create_text(wks, panel_strings(1), txres)
    amid(1)     = gsn_add_annotation(panel_vtan, txid(1), amres)

  ; Vertical velocity (filled contours)
    opts_vvel2@vpXF                       = xpos_plot(2)
    opts_vvel2@vpYF                       = ypos_plot(2)
    opts_vvel2@vpWidthF                   = xsize
    opts_vvel2@vpHeightF                  = ysize
    panel_vvel  = gsn_csm_contour(wks,plot_arr(14,:,:),opts_vvel2) ; [w]                            
    txid(2)     = gsn_create_text(wks, panel_strings(2), txres)
    amid(2)     = gsn_add_annotation(panel_vvel, txid(2), amres)

  ; MIDDLE ROW 

  ; Absolute angular momentum (filled contours)
    opts_aam2@vpXF                         = xpos_plot(3)
    opts_aam2@vpYF                         = ypos_plot(3)
    opts_aam2@vpWidthF                     = xsize
    opts_aam2@vpHeightF                    = ysize
    panel_aam   = gsn_csm_contour(wks,plot_arr(15,:,:),opts_aam2)  ; [AAM]
    txid(3)     = gsn_create_text(wks, panel_strings(3), txres)
    amid(3)     = gsn_add_annotation(panel_aam, txid(3), amres)

  ; Combined mean terms 
    opts_mbgt@vpXF       = xpos_plot(4)
    opts_mbgt@vpYF       = ypos_plot(4)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel0      = gsn_csm_contour(wks,plot_arr(8,:,:),opts_mbgt) ; [V_mzeta]
    cont0       = gsn_csm_contour(wks,plot_arr(8,:,:),opts_zero) ; Draw contours

    panel0_ovr  = gsn_csm_contour(wks,plot_arr(12,:,:),opts_vrad) ; [u]                             
    panel0_ovr2 = gsn_csm_contour(wks,plot_arr(14,:,:),opts_vvel) ; [w]                             

    overlay(panel0, cont0)                                   ; Overlay contours
    overlay(panel0, panel0_ovr)
    overlay(panel0, panel0_ovr2)
    txid(4)     = gsn_create_text(wks, panel_strings(4), txres)
    amid(4)     = gsn_add_annotation(panel0, txid(4), amres)
    delete([/cont0/])

  ; Combined eddy terms 
    opts_mbgt@vpXF       = xpos_plot(5)
    opts_mbgt@vpYF       = ypos_plot(5)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel1      = gsn_csm_contour(wks,plot_arr(9,:,:),opts_mbgt) ; [V_mv]
    cont0       = gsn_csm_contour(wks,plot_arr(9,:,:),opts_zero) ; Draw contours

    panel1_ovr  = gsn_csm_contour(wks,plot_arr(12,:,:),opts_vrad) ; [U]                             
    panel1_ovr2 = gsn_csm_contour(wks,plot_arr(14,:,:),opts_vvel) ; [w]                             

    overlay(panel1, cont0)                                   ; Overlay contours
    overlay(panel1, panel1_ovr)
    overlay(panel1, panel1_ovr2)
    txid(5)     = gsn_create_text(wks, panel_strings(5), txres)
    amid(5)     = gsn_add_annotation(panel1, txid(5), amres)
    delete([/cont0/])

  ; BOTTOM ROW 

  ; Total RHS [All V tend]
    opts_mbgt@vpXF       = xpos_plot(6)
    opts_mbgt@vpYF       = ypos_plot(6)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize

    panel6      = gsn_csm_contour(wks,rhs_tot(:,:),opts_mbgt)    ; [All V tend]
    cont0       = gsn_csm_contour(wks,rhs_tot(:,:),opts_zero)    ; Draw contours
    panel6_ovr  = gsn_csm_contour(wks,plot_arr(12,:,:),opts_vrad)   ; [U]
    panel6_ovr2 = gsn_csm_contour(wks,plot_arr(14,:,:),opts_vvel)   ; [w]

    overlay(panel6, cont0)                                   ; Overlay contours
    overlay(panel6, panel6_ovr)
    overlay(panel6, panel6_ovr2)
    txid(6)    = gsn_create_text(wks, panel_strings(6), txres)
    amid(6)    = gsn_add_annotation(panel6, txid(6), amres)
    delete([/cont0/])

  ; Local tendency of the mean tangential wind [V_t]
    opts_mbgt@vpXF       = xpos_plot(7)
    opts_mbgt@vpYF       = ypos_plot(7)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize

  ; Multiply to get correct units
;    if (res .eq. "5min" .and. sc .eq. 6) then   
;     plot_arr(0,:,:) = plot_arr(0,:,:) * 12
    if (res .eq. "5min" .and. sc .eq. 4 .and. typ .eq. "mono") then
     plot_arr(0,:,:) = plot_arr(0,:,:) * 12
;    elseif (res .eq. "5min" .and. sc .eq. 5 .and. typ .ne. "ring") then
;     plot_arr(0,:,:) = plot_arr(0,:,:) * 12
    else
     plot_arr(0,:,:) = plot_arr(0,:,:)
    end if

    panel7      = gsn_csm_contour(wks,plot_arr(0,:,:),opts_mbgt) ; [V_t]
    cont0       = gsn_csm_contour(wks,plot_arr(0,:,:),opts_zero) ; Draw contours

    panel7_ovr  = gsn_csm_contour(wks,plot_arr(12,:,:),opts_vrad) ; [U]
    panel7_ovr2 = gsn_csm_contour(wks,plot_arr(14,:,:),opts_vvel) ; [w]

    overlay(panel7, cont0)                                   ; Overlay contours
    overlay(panel7, panel7_ovr)
    overlay(panel7, panel7_ovr2)
    txid(7)    = gsn_create_text(wks, panel_strings(7), txres)
    amid(7)    = gsn_add_annotation(panel7, txid(7), amres)
    delete([/cont0/])

   end if 

  elseif (typ .eq. "comp") then 

;==============
;   TOP ROW
;==============

  ; 10/10/2019 --> Problem: hadn't renamed 'panel0', etc for other rows --> overwriting 
  ; 10/10/2019 --> Solution: create separate objects for each panel

  ; Mean radial influx of absolute vertical vorticity [V_mzeta]                           
    opts_mbgt@vpXF       = 0.13
    opts_mbgt@vpYF       = 0.99
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel0      = gsn_csm_contour(wks,mean_arr(1,:,:),opts_mbgt)  ; [V_mzeta]             
    cont0       = gsn_csm_contour(wks,mean_arr(1,:,:),opts_zero)  ; Draw contours         

    panel0_ovr  = gsn_csm_contour(wks,mean_arr(12,:,:),opts_vrad) ; [u]                   
    panel0_ovr2 = gsn_csm_contour(wks,mean_arr(14,:,:),opts_vvel) ; [w]                   

  ; Overlay contours                                                                      
    overlay(panel0, cont0)
    overlay(panel0, panel0_ovr)
    overlay(panel0, panel0_ovr2)
    txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
    amid(0)     = gsn_add_annotation(panel0, txid(0), amres)
    delete([/cont0/])

  ; Mean vertical advection of mean tangential momentum [V_mv]                            
    opts_mbgt@vpXF       = 0.34
    opts_mbgt@vpYF       = 0.99
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel1      = gsn_csm_contour(wks,mean_arr(2,:,:),opts_mbgt)  ; [V_mv]                
    cont0       = gsn_csm_contour(wks,mean_arr(2,:,:),opts_zero)  ; Draw contours         

    panel1_ovr  = gsn_csm_contour(wks,mean_arr(12,:,:),opts_vrad) ; [U]                   
    panel1_ovr2 = gsn_csm_contour(wks,mean_arr(14,:,:),opts_vvel) ; [w]                   

  ; Overlay contours                                                                      
    overlay(panel1, cont0)
    overlay(panel1, panel1_ovr)
    overlay(panel1, panel1_ovr2)
    txid(1)     = gsn_create_text(wks, panel_strings(1), txres)
    amid(1)     = gsn_add_annotation(panel1, txid(1), amres)
    delete([/cont0/])

  ; Eddy radial vorticity flux [V_ezeta]                                                  
    opts_mbgt@vpXF       = 0.55
    opts_mbgt@vpYF       = 0.99
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel2      = gsn_csm_contour(wks,mean_arr(3,:,:),opts_mbgt)  ; [V_ezeta]             
    cont0       = gsn_csm_contour(wks,mean_arr(3,:,:),opts_zero)  ; Draw contours         

    panel2_ovr  = gsn_csm_contour(wks,mean_arr(12,:,:),opts_vrad) ; [U]                   
    panel2_ovr2 = gsn_csm_contour(wks,mean_arr(14,:,:),opts_vvel) ; [w]                   

  ; Overlay contours                                                                      
    overlay(panel2, cont0)
    overlay(panel2, panel2_ovr)
    overlay(panel2, panel2_ovr2)
    txid(2)     = gsn_create_text(wks, panel_strings(2), txres)
    amid(2)     = gsn_add_annotation(panel2, txid(2), amres)
    delete([/cont0/])

  ; Vertical advection of eddy tangential momentum [V_ev]                                 
    opts_mbgt@vpXF       = 0.76
    opts_mbgt@vpYF       = 0.99
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel3      = gsn_csm_contour(wks,mean_arr(4,:,:),opts_mbgt)  ; [V_ev]                
    cont0       = gsn_csm_contour(wks,mean_arr(4,:,:),opts_zero)  ; Draw contours         

    panel3_ovr  = gsn_csm_contour(wks,mean_arr(12,:,:),opts_vrad) ; [U]                   
    panel3_ovr2 = gsn_csm_contour(wks,mean_arr(14,:,:),opts_vvel) ; [w]                   

  ; Overlay contours                                                                      
    overlay(panel3, cont0)
    overlay(panel3, panel3_ovr)
    overlay(panel3, panel3_ovr2)
    txid(3)     = gsn_create_text(wks, panel_strings(3), txres)
    amid(3)     = gsn_add_annotation(panel3, txid(3), amres)
    delete([/cont0/])

;=================
;   SECOND ROW 
;=================

  ; Mean radial influx of absolute vertical vorticity [V_mzeta]
    opts_mbgt@vpXF       = 0.13
    opts_mbgt@vpYF       = 0.80
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15

    plot_arr      = r2m_arr - mean_arr
    plot_arr!0    = "diag"
    plot_arr&diag = diag_arr
    plot_arr!1    = "lev"
    plot_arr&lev  = hyb_plot
    plot_arr!2    = "rad"
    plot_arr&rad  = radii

    panel4      = gsn_csm_contour(wks,plot_arr(1,:,:),opts_mbgt)  ; [V_mzeta]
    cont0       = gsn_csm_contour(wks,plot_arr(1,:,:),opts_zero)  ; Draw contours

    panel4_ovr  = gsn_csm_contour(wks,ring_arr(12,:,:),opts_vrad) ; [u]
    panel4_ovr2 = gsn_csm_contour(wks,ring_arr(14,:,:),opts_vvel) ; [w]

  ; Overlay contours
    overlay(panel4, cont0)
    overlay(panel4, panel4_ovr)
    overlay(panel4, panel4_ovr2)
    txid(4)     = gsn_create_text(wks, panel_strings(4), txres)
    amid(4)     = gsn_add_annotation(panel4, txid(4), amres)
    delete([/cont0/])

  ; Mean vertical advection of mean tangential momentum [V_mv]
    opts_mbgt@vpXF       = 0.34
    opts_mbgt@vpYF       = 0.80
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel5      = gsn_csm_contour(wks,plot_arr(2,:,:),opts_mbgt)  ; [V_mv]
    cont0       = gsn_csm_contour(wks,plot_arr(2,:,:),opts_zero)  ; Draw contours

    panel5_ovr  = gsn_csm_contour(wks,ring_arr(12,:,:),opts_vrad) ; [U]
    panel5_ovr2 = gsn_csm_contour(wks,ring_arr(14,:,:),opts_vvel) ; [w]

  ; Overlay contours
    overlay(panel5, cont0)
    overlay(panel5, panel5_ovr)
    overlay(panel5, panel5_ovr2)
    txid(5)     = gsn_create_text(wks, panel_strings(5), txres)
    amid(5)     = gsn_add_annotation(panel5, txid(5), amres)
    delete([/cont0/])

  ; Eddy radial vorticity flux [V_ezeta]
    opts_mbgt@vpXF       = 0.55
    opts_mbgt@vpYF       = 0.80
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel6      = gsn_csm_contour(wks,plot_arr(3,:,:),opts_mbgt)  ; [V_ezeta]               
    cont0       = gsn_csm_contour(wks,plot_arr(3,:,:),opts_zero)  ; Draw contours           

    panel6_ovr  = gsn_csm_contour(wks,ring_arr(12,:,:),opts_vrad) ; [U]                     
    panel6_ovr2 = gsn_csm_contour(wks,ring_arr(14,:,:),opts_vvel) ; [w]                     

  ; Overlay contours                                                                        
    overlay(panel6, cont0)
    overlay(panel6, panel6_ovr)
    overlay(panel6, panel6_ovr2)
    txid(6)     = gsn_create_text(wks, panel_strings(6), txres)
    amid(6)     = gsn_add_annotation(panel6, txid(6), amres)
    delete([/cont0/])

  ; Vertical advection of eddy tangential momentum [V_ev]                                 
    opts_mbgt@vpXF       = 0.76
    opts_mbgt@vpYF       = 0.80
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel7      = gsn_csm_contour(wks,plot_arr(4,:,:),opts_mbgt)  ; [V_ev]                
    cont0       = gsn_csm_contour(wks,plot_arr(4,:,:),opts_zero)  ; Draw contours         

    panel7_ovr  = gsn_csm_contour(wks,ring_arr(12,:,:),opts_vrad) ; [U]                   
    panel7_ovr2 = gsn_csm_contour(wks,ring_arr(14,:,:),opts_vvel) ; [w]                   

  ; Overlay contours                                                                      
    overlay(panel7, cont0)
    overlay(panel7, panel7_ovr)
    overlay(panel7, panel7_ovr2)
    txid(7)     = gsn_create_text(wks, panel_strings(7), txres)
    amid(7)     = gsn_add_annotation(panel7, txid(7), amres)
    delete([/cont0, plot_arr/])

;==============
;  THIRD ROW 
;==============

  ; Mean radial influx of absolute vertical vorticity [V_mzeta]                            
    opts_mbgt@vpXF       = 0.13
    opts_mbgt@vpYF       = 0.61
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15

    plot_arr      = mono_arr - mean_arr
    plot_arr!0    = "diag"
    plot_arr&diag = diag_arr
    plot_arr!1    = "lev"
    plot_arr&lev  = hyb_plot
    plot_arr!2    = "rad"
    plot_arr&rad  = radii

    panel8      = gsn_csm_contour(wks,plot_arr(1,:,:),opts_mbgt)  ; [V_mzeta]              
    cont0       = gsn_csm_contour(wks,plot_arr(1,:,:),opts_zero)  ; Draw contours          

    panel8_ovr  = gsn_csm_contour(wks,r2m_arr(12,:,:),opts_vrad) ; [u]                    
    panel8_ovr2 = gsn_csm_contour(wks,r2m_arr(14,:,:),opts_vvel) ; [w]                    

  ; Overlay contours                                                                       
    overlay(panel8, cont0)
    overlay(panel8, panel8_ovr)
    overlay(panel8, panel8_ovr2)
    txid(8)     = gsn_create_text(wks, panel_strings(8), txres)
    amid(8)     = gsn_add_annotation(panel8, txid(8), amres)
    delete([/cont0/])

  ; Mean vertical advection of mean tangential momentum [V_mv]                             
    opts_mbgt@vpXF       = 0.34
    opts_mbgt@vpYF       = 0.61
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel9      = gsn_csm_contour(wks,plot_arr(2,:,:),opts_mbgt)  ; [V_mv]                 
    cont0       = gsn_csm_contour(wks,plot_arr(2,:,:),opts_zero)  ; Draw contours          

    panel9_ovr  = gsn_csm_contour(wks,r2m_arr(12,:,:),opts_vrad) ; [U]                    
    panel9_ovr2 = gsn_csm_contour(wks,r2m_arr(14,:,:),opts_vvel) ; [w]                    

  ; Overlay contours                                                                       
    overlay(panel9, cont0)
    overlay(panel9, panel9_ovr)
    overlay(panel9, panel9_ovr2)
    txid(9)     = gsn_create_text(wks, panel_strings(9), txres)
    amid(9)     = gsn_add_annotation(panel9, txid(9), amres)
    delete([/cont0/])

  ; Eddy radial vorticity flux [V_ezeta]
    opts_mbgt@vpXF       = 0.55
    opts_mbgt@vpYF       = 0.61
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel10      = gsn_csm_contour(wks,plot_arr(3,:,:),opts_mbgt)  ; [V_ezeta]               
    cont0        = gsn_csm_contour(wks,plot_arr(3,:,:),opts_zero)  ; Draw contours           

    panel10_ovr  = gsn_csm_contour(wks,r2m_arr(12,:,:),opts_vrad) ; [U]                     
    panel10_ovr2 = gsn_csm_contour(wks,r2m_arr(14,:,:),opts_vvel) ; [w]                     

  ; Overlay contours                                                                        
    overlay(panel10, cont0)
    overlay(panel10, panel10_ovr)
    overlay(panel10, panel10_ovr2)
    txid(10)    = gsn_create_text(wks, panel_strings(10), txres)
    amid(10)    = gsn_add_annotation(panel10, txid(10), amres)
    delete([/cont0/])

  ; Vertical advection of eddy tangential momentum [V_ev]                                   
    opts_mbgt@vpXF       = 0.76
    opts_mbgt@vpYF       = 0.61
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel11      = gsn_csm_contour(wks,plot_arr(4,:,:),opts_mbgt)  ; [V_ev]                  
    cont0        = gsn_csm_contour(wks,plot_arr(4,:,:),opts_zero)  ; Draw contours           

    panel11_ovr  = gsn_csm_contour(wks,r2m_arr(12,:,:),opts_vrad) ; [U]                     
    panel11_ovr2 = gsn_csm_contour(wks,r2m_arr(14,:,:),opts_vvel) ; [w]                     

  ; Overlay contours                                                                        
    overlay(panel11, cont0)
    overlay(panel11, panel11_ovr)
    overlay(panel11, panel11_ovr2)
    txid(11)     = gsn_create_text(wks, panel_strings(11), txres)
    amid(11)     = gsn_add_annotation(panel11, txid(11), amres)
    delete([/cont0, plot_arr/])

;===============
;  BOTTOM ROW 
;===============

  ; Mean radial influx of absolute vertical vorticity [V_mzeta]                            
    opts_mbgt@vpXF       = 0.13
    opts_mbgt@vpYF       = 0.42
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15

    plot_arr      = m2r_arr - mean_arr
    plot_arr!0    = "diag"
    plot_arr&diag = diag_arr
    plot_arr!1    = "lev"
    plot_arr&lev  = hyb_plot
    plot_arr!2    = "rad"
    plot_arr&rad  = radii

    panel12     = gsn_csm_contour(wks,plot_arr(1,:,:),opts_mbgt)  ; [V_mzeta]              
    cont0       = gsn_csm_contour(wks,plot_arr(1,:,:),opts_zero)  ; Draw contours          

    panel12_ovr  = gsn_csm_contour(wks,mono_arr(12,:,:),opts_vrad) ; [u]                    
    panel12_ovr2 = gsn_csm_contour(wks,mono_arr(14,:,:),opts_vvel) ; [w]                    

  ; Overlay contours                                                                       
    overlay(panel12, cont0)
    overlay(panel12, panel12_ovr)
    overlay(panel12, panel12_ovr2)
    txid(12)    = gsn_create_text(wks, panel_strings(12), txres)
    amid(12)    = gsn_add_annotation(panel12, txid(12), amres)
    delete([/cont0/])

  ; Mean vertical advection of mean tangential momentum [V_mv]                             
    opts_mbgt@vpXF       = 0.34
    opts_mbgt@vpYF       = 0.42
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel13      = gsn_csm_contour(wks,plot_arr(2,:,:),opts_mbgt)  ; [V_mv]                 
    cont0       = gsn_csm_contour(wks,plot_arr(2,:,:),opts_zero)  ; Draw contours          

    panel13_ovr  = gsn_csm_contour(wks,mono_arr(12,:,:),opts_vrad) ; [U]                    
    panel13_ovr2 = gsn_csm_contour(wks,mono_arr(14,:,:),opts_vvel) ; [w]                    

  ; Overlay contours                                                                       
    overlay(panel13, cont0)
    overlay(panel13, panel13_ovr)
    overlay(panel13, panel13_ovr2)
    txid(13)    = gsn_create_text(wks, panel_strings(13), txres)
    amid(13)    = gsn_add_annotation(panel13, txid(13), amres)
    delete([/cont0/])

  ; Eddy radial vorticity flux [V_ezeta]
    opts_mbgt@vpXF       = 0.55
    opts_mbgt@vpYF       = 0.42
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel14      = gsn_csm_contour(wks,plot_arr(3,:,:),opts_mbgt)  ; [V_ezeta]               
    cont0       = gsn_csm_contour(wks,plot_arr(3,:,:),opts_zero)  ; Draw contours           

    panel14_ovr  = gsn_csm_contour(wks,mono_arr(12,:,:),opts_vrad) ; [U]                     
    panel14_ovr2 = gsn_csm_contour(wks,mono_arr(14,:,:),opts_vvel) ; [w]                     

  ; Overlay contours                                                                        
    overlay(panel14, cont0)
    overlay(panel14, panel14_ovr)
    overlay(panel14, panel14_ovr2)
    txid(14)    = gsn_create_text(wks, panel_strings(14), txres)
    amid(14)    = gsn_add_annotation(panel14, txid(14), amres)
    delete([/cont0/])

  ; Vertical advection of eddy tangential momentum [V_ev]                                   
    opts_mbgt@vpXF       = 0.76
    opts_mbgt@vpYF       = 0.42
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel15      = gsn_csm_contour(wks,plot_arr(4,:,:),opts_mbgt)  ; [V_ev]                  
    cont0       = gsn_csm_contour(wks,plot_arr(4,:,:),opts_zero)  ; Draw contours           

    panel15_ovr  = gsn_csm_contour(wks,mono_arr(12,:,:),opts_vrad) ; [U]                     
    panel15_ovr2 = gsn_csm_contour(wks,mono_arr(14,:,:),opts_vvel) ; [w]                     

  ; Overlay contours                                                                        
    overlay(panel15, cont0)
    overlay(panel15, panel15_ovr)
    overlay(panel15, panel15_ovr2)
    txid(15)    = gsn_create_text(wks, panel_strings(15), txres)
    amid(15)    = gsn_add_annotation(panel15, txid(15), amres)
    delete([/cont0, plot_arr/])

  elseif (typ .eq. "mean_eddy") then 

    if (col .eq. 3) then 

     xpos_plot = (/0.07, 0.39, 0.71, 0.07, 0.39, 0.71/)
     ypos_plot = (/0.97, 0.97, 0.97, 0.65, 0.65, 0.65/)
     xsize     = 0.24
     ysize     = 0.24

    elseif (col .eq. 4) then 

     xpos_plot = (/0.05, 0.30, 0.55, 0.80, 0.05, 0.30, 0.55, 0.80/)
     ypos_plot = (/0.95, 0.95, 0.95, 0.95, 0.72, 0.72, 0.72, 0.72/)
     xsize     = 0.16
     ysize     = 0.16

    end if 

  ;===================
  ; Ring-like phase 
  ;===================

  ; Mean
    opts_mbgt@vpXF       = xpos_plot(0)
    opts_mbgt@vpYF       = ypos_plot(0)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    if (p0 .eq. "ring") then 
     panel0      = gsn_csm_contour(wks,ring_arr(8,:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,ring_arr(8,:,:),opts_zero)
     panel0_ovr  = gsn_csm_contour(wks,ring_arr(12,:,:),opts_vrad)
     panel0_ovr2 = gsn_csm_contour(wks,ring_arr(14,:,:),opts_vvel)
    elseif (p0 .eq. "r2m") then 
     panel0      = gsn_csm_contour(wks,r2m_arr(8,:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,r2m_arr(8,:,:),opts_zero)
     panel0_ovr  = gsn_csm_contour(wks,r2m_arr(12,:,:),opts_vrad)
     panel0_ovr2 = gsn_csm_contour(wks,r2m_arr(14,:,:),opts_vvel)
    elseif (p0 .eq. "mono") then 
     panel0      = gsn_csm_contour(wks,mono_arr(8,:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,mono_arr(8,:,:),opts_zero)
     panel0_ovr  = gsn_csm_contour(wks,mono_arr(12,:,:),opts_vrad)
     panel0_ovr2 = gsn_csm_contour(wks,mono_arr(14,:,:),opts_vvel)
    elseif (p0 .eq. "m2r") then 
     panel0      = gsn_csm_contour(wks,m2r_arr(8,:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,m2r_arr(8,:,:),opts_zero)
     panel0_ovr  = gsn_csm_contour(wks,m2r_arr(12,:,:),opts_vrad)
     panel0_ovr2 = gsn_csm_contour(wks,m2r_arr(14,:,:),opts_vvel)
    end if 

  ; Overlay contours
    overlay(panel0, cont0)
    overlay(panel0, panel0_ovr)
    overlay(panel0, panel0_ovr2)
    txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
    amid(0)     = gsn_add_annotation(panel0, txid(0), amres)
    delete([/cont0/])

  ; Eddies
    opts_mbgt@vpXF       = xpos_plot(1)
    opts_mbgt@vpYF       = ypos_plot(1)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    if (p0 .eq. "ring") then
     panel1      = gsn_csm_contour(wks,ring_arr(9,:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,ring_arr(9,:,:),opts_zero)
     panel1_ovr  = gsn_csm_contour(wks,ring_arr(12,:,:),opts_vrad)
     panel1_ovr2 = gsn_csm_contour(wks,ring_arr(14,:,:),opts_vvel)
   elseif (p0 .eq. "r2m") then
     panel1      = gsn_csm_contour(wks,r2m_arr(9,:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,r2m_arr(9,:,:),opts_zero)
     panel1_ovr  = gsn_csm_contour(wks,r2m_arr(12,:,:),opts_vrad)
     panel1_ovr2 = gsn_csm_contour(wks,r2m_arr(14,:,:),opts_vvel)
    elseif (p0 .eq. "mono") then
     panel1      = gsn_csm_contour(wks,mono_arr(9,:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,mono_arr(9,:,:),opts_zero)
     panel1_ovr  = gsn_csm_contour(wks,mono_arr(12,:,:),opts_vrad)
     panel1_ovr2 = gsn_csm_contour(wks,mono_arr(14,:,:),opts_vvel)
   elseif (p0 .eq. "m2r") then
     panel1      = gsn_csm_contour(wks,m2r_arr(9,:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,m2r_arr(9,:,:),opts_zero)
     panel1_ovr  = gsn_csm_contour(wks,m2r_arr(12,:,:),opts_vrad)
     panel1_ovr2 = gsn_csm_contour(wks,m2r_arr(14,:,:),opts_vvel)
    end if

  ; Overlay contours
    overlay(panel1, cont0)
    overlay(panel1, panel1_ovr)
    overlay(panel1, panel1_ovr2)
    txid(1)     = gsn_create_text(wks, panel_strings(1), txres)
    amid(1)     = gsn_add_annotation(panel1, txid(1), amres)
    delete([/cont0/])

  ; Mean + eddies 
    opts_mbgt@vpXF       = xpos_plot(2)
    opts_mbgt@vpYF       = ypos_plot(2)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    if (p0 .eq. "ring") then
     panel2      = gsn_csm_contour(wks,ring_sum(:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,ring_sum(:,:),opts_zero)
     panel2_ovr  = gsn_csm_contour(wks,ring_arr(12,:,:),opts_vrad)
     panel2_ovr2 = gsn_csm_contour(wks,ring_arr(14,:,:),opts_vvel)
   elseif (p0 .eq. "r2m") then
     panel2      = gsn_csm_contour(wks,r2m_sum(:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,r2m_sum(:,:),opts_zero)
     panel2_ovr  = gsn_csm_contour(wks,r2m_arr(12,:,:),opts_vrad)
     panel2_ovr2 = gsn_csm_contour(wks,r2m_arr(14,:,:),opts_vvel)
    elseif (p0 .eq. "mono") then
     panel2      = gsn_csm_contour(wks,mono_sum(:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,mono_sum(:,:),opts_zero)
     panel2_ovr  = gsn_csm_contour(wks,mono_arr(12,:,:),opts_vrad)
     panel2_ovr2 = gsn_csm_contour(wks,mono_arr(14,:,:),opts_vvel)
   elseif (p0 .eq. "m2r") then
     panel2      = gsn_csm_contour(wks,m2r_sum(:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,m2r_sum(:,:),opts_zero)
     panel2_ovr  = gsn_csm_contour(wks,m2r_arr(12,:,:),opts_vrad)
     panel2_ovr2 = gsn_csm_contour(wks,m2r_arr(14,:,:),opts_vvel)
    end if

  ; Overlay contours 
    overlay(panel2, cont0)
    overlay(panel2, panel2_ovr)
    overlay(panel2, panel2_ovr2)
    txid(2)     = gsn_create_text(wks, panel_strings(2), txres)
    amid(2)     = gsn_add_annotation(panel2, txid(2), amres)
    delete([/cont0/])

  ;====================
  ; Monopole phase 
  ;====================

  ; Mean                                                                             
    opts_mbgt@vpXF       = xpos_plot(3)
    opts_mbgt@vpYF       = ypos_plot(3)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    if (p1 .eq. "ring") then
     panel3      = gsn_csm_contour(wks,ring_arr(8,:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,ring_arr(8,:,:),opts_zero)
     panel3_ovr  = gsn_csm_contour(wks,ring_arr(12,:,:),opts_vrad)
     panel3_ovr2 = gsn_csm_contour(wks,ring_arr(14,:,:),opts_vvel)
    elseif (p1 .eq. "r2m") then
     panel3      = gsn_csm_contour(wks,r2m_arr(8,:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,r2m_arr(8,:,:),opts_zero)
     panel3_ovr  = gsn_csm_contour(wks,r2m_arr(12,:,:),opts_vrad)
     panel3_ovr2 = gsn_csm_contour(wks,r2m_arr(14,:,:),opts_vvel)
    elseif (p1 .eq. "mono") then
     panel3      = gsn_csm_contour(wks,mono_arr(8,:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,mono_arr(8,:,:),opts_zero)
     panel3_ovr  = gsn_csm_contour(wks,mono_arr(12,:,:),opts_vrad)
     panel3_ovr2 = gsn_csm_contour(wks,mono_arr(14,:,:),opts_vvel)
    elseif (p1 .eq. "m2r") then
     panel3      = gsn_csm_contour(wks,m2r_arr(8,:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,m2r_arr(8,:,:),opts_zero)
     panel3_ovr  = gsn_csm_contour(wks,m2r_arr(12,:,:),opts_vrad)
     panel3_ovr2 = gsn_csm_contour(wks,m2r_arr(14,:,:),opts_vvel)
    end if

  ; Overlay contours                                                                 
    overlay(panel3, cont0)
    overlay(panel3, panel3_ovr)
    overlay(panel3, panel3_ovr2)
    txid(3)     = gsn_create_text(wks, panel_strings(3), txres)
    amid(3)     = gsn_add_annotation(panel3, txid(3), amres)
    delete([/cont0/])

  ; Eddies                                                                           
    opts_mbgt@vpXF       = xpos_plot(4)
    opts_mbgt@vpYF       = ypos_plot(4)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    if (p1 .eq. "ring") then
     panel4      = gsn_csm_contour(wks,ring_arr(9,:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,ring_arr(9,:,:),opts_zero)
     panel4_ovr  = gsn_csm_contour(wks,ring_arr(12,:,:),opts_vrad)
     panel4_ovr2 = gsn_csm_contour(wks,ring_arr(14,:,:),opts_vvel)
   elseif (p1 .eq. "r2m") then
     panel4      = gsn_csm_contour(wks,r2m_arr(9,:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,r2m_arr(9,:,:),opts_zero)
     panel4_ovr  = gsn_csm_contour(wks,r2m_arr(12,:,:),opts_vrad)
     panel4_ovr2 = gsn_csm_contour(wks,r2m_arr(14,:,:),opts_vvel)
    elseif (p1 .eq. "mono") then
     panel4      = gsn_csm_contour(wks,mono_arr(9,:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,mono_arr(9,:,:),opts_zero)
     panel4_ovr  = gsn_csm_contour(wks,mono_arr(12,:,:),opts_vrad)
     panel4_ovr2 = gsn_csm_contour(wks,mono_arr(14,:,:),opts_vvel)
   elseif (p1 .eq. "m2r") then
     panel4      = gsn_csm_contour(wks,m2r_arr(9,:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,m2r_arr(9,:,:),opts_zero)
     panel4_ovr  = gsn_csm_contour(wks,m2r_arr(12,:,:),opts_vrad)
     panel4_ovr2 = gsn_csm_contour(wks,m2r_arr(14,:,:),opts_vvel)
    end if

  ; Overlay contours                                                                 
    overlay(panel4, cont0)
    overlay(panel4, panel4_ovr)
    overlay(panel4, panel4_ovr2)
    txid(4)     = gsn_create_text(wks, panel_strings(4), txres)
    amid(4)     = gsn_add_annotation(panel4, txid(4), amres)
    delete([/cont0/])

  ; Mean + eddies                                                                            
    opts_mbgt@vpXF       = xpos_plot(5)
    opts_mbgt@vpYF       = ypos_plot(5)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    if (p1 .eq. "ring") then
     panel5      = gsn_csm_contour(wks,ring_sum(:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,ring_sum(:,:),opts_zero)
     panel5_ovr  = gsn_csm_contour(wks,ring_arr(12,:,:),opts_vrad)
     panel5_ovr2 = gsn_csm_contour(wks,ring_arr(14,:,:),opts_vvel)
   elseif (p1 .eq. "r2m") then
     panel5      = gsn_csm_contour(wks,r2m_sum(:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,r2m_sum(:,:),opts_zero)
     panel5_ovr  = gsn_csm_contour(wks,r2m_arr(12,:,:),opts_vrad)
     panel5_ovr2 = gsn_csm_contour(wks,r2m_arr(14,:,:),opts_vvel)
    elseif (p1 .eq. "mono") then
     panel5      = gsn_csm_contour(wks,mono_sum(:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,mono_sum(:,:),opts_zero)
     panel5_ovr  = gsn_csm_contour(wks,mono_arr(12,:,:),opts_vrad)
     panel5_ovr2 = gsn_csm_contour(wks,mono_arr(14,:,:),opts_vvel)
   elseif (p1 .eq. "m2r") then
     panel5      = gsn_csm_contour(wks,m2r_sum(:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,m2r_sum(:,:),opts_zero)
     panel5_ovr  = gsn_csm_contour(wks,m2r_arr(12,:,:),opts_vrad)
     panel5_ovr2 = gsn_csm_contour(wks,m2r_arr(14,:,:),opts_vvel)
    end if

  ; Overlay contours
    overlay(panel5, cont0)
    overlay(panel5, panel5_ovr)
    overlay(panel5, panel5_ovr2)
    txid(5)     = gsn_create_text(wks, panel_strings(5), txres)
    amid(5)     = gsn_add_annotation(panel5, txid(5), amres)
    delete([/cont0/])

; Radial wind or local tangential wind tendency for non-VC + single VC phase
  elseif (typ .eq. "novc_rad") then 

   ;=====================
   ; Non-VC composite 
   ;=====================
     opts_mbgt@vpXF       = xpos_plot(0)
     opts_mbgt@vpYF       = ypos_plot(0)
     opts_mbgt@vpWidthF   = xsize
     opts_mbgt@vpHeightF  = ysize

   ; Produce filled contours of radial wind
     panel0      = gsn_csm_contour(wks,novc_arr(12,:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,novc_arr(12,:,:),opts_zero)

   ; Overlay AAM and vertical velocity
     panel0_ovr  = gsn_csm_contour(wks,novc_arr(15,:,:),opts_aam)
     panel0_ovr2 = gsn_csm_contour(wks,novc_arr(14,:,:),opts_vvel)

   ; Overlay contours
     overlay(panel0, cont0)
     overlay(panel0, panel0_ovr)
     overlay(panel0, panel0_ovr2)
     txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
     amid(0)     = gsn_add_annotation(panel0, txid(0), amres)
     delete([/cont0/])

   ;=================================
   ; VC composite (selected phase)
   ;================================= 
     opts_mbgt@vpXF       = xpos_plot(1)
     opts_mbgt@vpYF       = ypos_plot(1)
     opts_mbgt@vpWidthF   = xsize
     opts_mbgt@vpHeightF  = ysize

   ; Produce filled contours of radial wind, and overlay AAM + vertical velocity
     if (p0 .eq. "ring") then 
      panel1      = gsn_csm_contour(wks,ring_arr(12,:,:),opts_mbgt)
      cont0       = gsn_csm_contour(wks,ring_arr(12,:,:),opts_zero)
      panel1_ovr  = gsn_csm_contour(wks,ring_arr(15,:,:),opts_aam)
      panel1_ovr2 = gsn_csm_contour(wks,ring_arr(14,:,:),opts_vvel)
     elseif (p0 .eq. "r2m") then 
      panel1      = gsn_csm_contour(wks,r2m_arr(12,:,:),opts_mbgt)
      cont0       = gsn_csm_contour(wks,r2m_arr(12,:,:),opts_zero)
      panel1_ovr  = gsn_csm_contour(wks,r2m_arr(15,:,:),opts_aam)
      panel1_ovr2 = gsn_csm_contour(wks,r2m_arr(14,:,:),opts_vvel)
     elseif (p0 .eq. "mono") then 
      panel1      = gsn_csm_contour(wks,mono_arr(12,:,:),opts_mbgt)
      cont0       = gsn_csm_contour(wks,mono_arr(12,:,:),opts_zero)
      panel1_ovr  = gsn_csm_contour(wks,mono_arr(15,:,:),opts_aam)
      panel1_ovr2 = gsn_csm_contour(wks,mono_arr(14,:,:),opts_vvel)
     elseif (p0 .eq. "m2r") then 
      panel1      = gsn_csm_contour(wks,m2r_arr(12,:,:),opts_mbgt)
      cont0       = gsn_csm_contour(wks,m2r_arr(12,:,:),opts_zero)
      panel1_ovr  = gsn_csm_contour(wks,m2r_arr(15,:,:),opts_aam)
      panel1_ovr2 = gsn_csm_contour(wks,m2r_arr(14,:,:),opts_vvel)
     end if 

   ; Overlay contours
     overlay(panel1, cont0)
     overlay(panel1, panel1_ovr)
     overlay(panel1, panel1_ovr2)
     txid(1)     = gsn_create_text(wks, panel_strings(1), txres)
     amid(1)     = gsn_add_annotation(panel1, txid(1), amres)
     delete([/cont0/])

   ;==================================
   ; 2nd VC composite (if selected)
   ;================================== 
   if (col .eq. 3) then 
     opts_mbgt@vpXF       = xpos_plot(2)
     opts_mbgt@vpYF       = ypos_plot(2)
     opts_mbgt@vpWidthF   = xsize
     opts_mbgt@vpHeightF  = ysize

   ; Produce filled contours of radial wind, and overlay AAM + vertical velocity
     if (p1 .eq. "ring") then
      panel2      = gsn_csm_contour(wks,ring_arr(12,:,:),opts_mbgt)
      cont0       = gsn_csm_contour(wks,ring_arr(12,:,:),opts_zero)
      panel2_ovr  = gsn_csm_contour(wks,ring_arr(15,:,:),opts_aam)
      panel2_ovr2 = gsn_csm_contour(wks,ring_arr(14,:,:),opts_vvel)
     elseif (p1 .eq. "r2m") then
      panel2      = gsn_csm_contour(wks,r2m_arr(12,:,:),opts_mbgt)
      cont0       = gsn_csm_contour(wks,r2m_arr(12,:,:),opts_zero)
      panel2_ovr  = gsn_csm_contour(wks,r2m_arr(15,:,:),opts_aam)
      panel2_ovr2 = gsn_csm_contour(wks,r2m_arr(14,:,:),opts_vvel)
     elseif (p1 .eq. "mono") then
      panel2      = gsn_csm_contour(wks,mono_arr(12,:,:),opts_mbgt)
      cont0       = gsn_csm_contour(wks,mono_arr(12,:,:),opts_zero)
      panel2_ovr  = gsn_csm_contour(wks,mono_arr(15,:,:),opts_aam)
      panel2_ovr2 = gsn_csm_contour(wks,mono_arr(14,:,:),opts_vvel)
     elseif (p1 .eq. "m2r") then
      panel2      = gsn_csm_contour(wks,m2r_arr(12,:,:),opts_mbgt)
      cont0       = gsn_csm_contour(wks,m2r_arr(12,:,:),opts_zero)
      panel2_ovr  = gsn_csm_contour(wks,m2r_arr(15,:,:),opts_aam)
      panel2_ovr2 = gsn_csm_contour(wks,m2r_arr(14,:,:),opts_vvel)
     end if

   ; Overlay contours
     overlay(panel2, cont0)
     overlay(panel2, panel2_ovr)
     overlay(panel2, panel2_ovr2)
     txid(2)     = gsn_create_text(wks, panel_strings(2), txres)
     amid(2)     = gsn_add_annotation(panel2, txid(2), amres)
     delete([/cont0/])

   end if 
 ; End IF statement (plot either 2 or 3 panels)

; Local tangential wind tendency (LHS) for non-VC composite 
  elseif (typ .eq. "novc_lhs") then 

   ;=====================                                                               
   ; Non-VC composite                                                                   
   ;=====================                                                               
     opts_mbgt@vpXF       = xpos_plot(0)
     opts_mbgt@vpYF       = ypos_plot(0)
     opts_mbgt@vpWidthF   = xsize
     opts_mbgt@vpHeightF  = ysize

   ; Produce filled contours of radial wind                                             
     panel0      = gsn_csm_contour(wks,novc_arr(0,:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,novc_arr(0,:,:),opts_zero)

   ; Overlay AAM and vertical velocity                                                  
     panel0_ovr  = gsn_csm_contour(wks,novc_arr(15,:,:),opts_aam)
     panel0_ovr2 = gsn_csm_contour(wks,novc_arr(14,:,:),opts_vvel)

   ; Overlay contours                                                                   
     overlay(panel0, cont0)
     overlay(panel0, panel0_ovr)
     overlay(panel0, panel0_ovr2)
     txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
     amid(0)     = gsn_add_annotation(panel0, txid(0), amres)
     delete([/cont0/])

   ;=================================                                                   
   ; VC composite (selected phase)                                                      
   ;=================================                                                   
     opts_mbgt@vpXF       = xpos_plot(1)
     opts_mbgt@vpYF       = ypos_plot(1)
     opts_mbgt@vpWidthF   = xsize
     opts_mbgt@vpHeightF  = ysize

   ; Produce filled contours of radial wind, and overlay AAM + vertical velocity        
     if (p0 .eq. "ring") then
      panel1      = gsn_csm_contour(wks,ring_arr(0,:,:),opts_mbgt)
      cont0       = gsn_csm_contour(wks,ring_arr(0,:,:),opts_zero)
      panel1_ovr  = gsn_csm_contour(wks,ring_arr(15,:,:),opts_aam)
      panel1_ovr2 = gsn_csm_contour(wks,ring_arr(14,:,:),opts_vvel)
     elseif (p0 .eq. "r2m") then
      panel1      = gsn_csm_contour(wks,r2m_arr(0,:,:),opts_mbgt)
      cont0       = gsn_csm_contour(wks,r2m_arr(0,:,:),opts_zero)
      panel1_ovr  = gsn_csm_contour(wks,r2m_arr(15,:,:),opts_aam)
      panel1_ovr2 = gsn_csm_contour(wks,r2m_arr(14,:,:),opts_vvel)
     elseif (p0 .eq. "mono") then
      panel1      = gsn_csm_contour(wks,mono_arr(0,:,:),opts_mbgt)
      cont0       = gsn_csm_contour(wks,mono_arr(0,:,:),opts_zero)
      panel1_ovr  = gsn_csm_contour(wks,mono_arr(15,:,:),opts_aam)
      panel1_ovr2 = gsn_csm_contour(wks,mono_arr(14,:,:),opts_vvel)
     elseif (p0 .eq. "m2r") then
      panel1      = gsn_csm_contour(wks,m2r_arr(0,:,:),opts_mbgt)
      cont0       = gsn_csm_contour(wks,m2r_arr(0,:,:),opts_zero)
      panel1_ovr  = gsn_csm_contour(wks,m2r_arr(15,:,:),opts_aam)
      panel1_ovr2 = gsn_csm_contour(wks,m2r_arr(14,:,:),opts_vvel)
     end if

   ; Overlay contours                                                                   
     overlay(panel1, cont0)
     overlay(panel1, panel1_ovr)
     overlay(panel1, panel1_ovr2)
     txid(1)     = gsn_create_text(wks, panel_strings(1), txres)
     amid(1)     = gsn_add_annotation(panel1, txid(1), amres)
     delete([/cont0/])

   ;==================================
   ; 2nd VC composite (if selected) 
   ;==================================
   if (col .eq. 3) then
     opts_mbgt@vpXF       = xpos_plot(2)
     opts_mbgt@vpYF       = ypos_plot(2)
     opts_mbgt@vpWidthF   = xsize
     opts_mbgt@vpHeightF  = ysize

   ; Produce filled contours of radial wind, and overlay AAM + vertical velocity
     if (p1 .eq. "ring") then
      panel2      = gsn_csm_contour(wks,ring_arr(0,:,:),opts_mbgt)
      cont0       = gsn_csm_contour(wks,ring_arr(0,:,:),opts_zero)
      panel2_ovr  = gsn_csm_contour(wks,ring_arr(15,:,:),opts_aam)
      panel2_ovr2 = gsn_csm_contour(wks,ring_arr(14,:,:),opts_vvel)
     elseif (p1 .eq. "r2m") then
      panel2      = gsn_csm_contour(wks,r2m_arr(0,:,:),opts_mbgt)
      cont0       = gsn_csm_contour(wks,r2m_arr(0,:,:),opts_zero)
      panel2_ovr  = gsn_csm_contour(wks,r2m_arr(15,:,:),opts_aam)
      panel2_ovr2 = gsn_csm_contour(wks,r2m_arr(14,:,:),opts_vvel)
     elseif (p1 .eq. "mono") then
      panel2      = gsn_csm_contour(wks,mono_arr(0,:,:),opts_mbgt)
      cont0       = gsn_csm_contour(wks,mono_arr(0,:,:),opts_zero)
      panel2_ovr  = gsn_csm_contour(wks,mono_arr(15,:,:),opts_aam)
      panel2_ovr2 = gsn_csm_contour(wks,mono_arr(14,:,:),opts_vvel)
     elseif (p1 .eq. "m2r") then
      panel2      = gsn_csm_contour(wks,m2r_arr(0,:,:),opts_mbgt)
      cont0       = gsn_csm_contour(wks,m2r_arr(0,:,:),opts_zero)
      panel2_ovr  = gsn_csm_contour(wks,m2r_arr(15,:,:),opts_aam)
      panel2_ovr2 = gsn_csm_contour(wks,m2r_arr(14,:,:),opts_vvel)
     end if

   ; Overlay contours
     overlay(panel2, cont0)
     overlay(panel2, panel2_ovr)
     overlay(panel2, panel2_ovr2)
     txid(2)     = gsn_create_text(wks, panel_strings(2), txres)
     amid(2)     = gsn_add_annotation(panel2, txid(2), amres)
     delete([/cont0/])

   end if
 ; End IF statement on plotting 2 or 3 panels

; Difference between ring-like and monopole phases, for the combined mean + eddy terms 
  elseif (typ .eq. "paper2_diff") then 

  ; Mean                                                                          
    opts_mbgt@vpXF       = xpos_plot(0)
    opts_mbgt@vpYF       = ypos_plot(0)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel0      = gsn_csm_contour(wks,diff_mean(:,:),opts_mbgt)
    cont0       = gsn_csm_contour(wks,diff_mean(:,:),opts_zero)
    panel0_ovr  = gsn_csm_contour(wks,mono_arr(12,:,:),opts_vrad)
    panel0_ovr2 = gsn_csm_contour(wks,mono_arr(14,:,:),opts_vvel)

  ; Overlay contours                                                              
    overlay(panel0, cont0)
    overlay(panel0, panel0_ovr)
    overlay(panel0, panel0_ovr2)
    txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
    amid(0)     = gsn_add_annotation(panel0, txid(0), amres)
    delete([/cont0/])

  ; Eddies                                                                        
    opts_mbgt@vpXF       = xpos_plot(1)
    opts_mbgt@vpYF       = ypos_plot(1)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel1      = gsn_csm_contour(wks,diff_eddy(:,:),opts_mbgt)
    cont0       = gsn_csm_contour(wks,diff_eddy(:,:),opts_zero)
    panel1_ovr  = gsn_csm_contour(wks,mono_arr(12,:,:),opts_vrad)
    panel1_ovr2 = gsn_csm_contour(wks,mono_arr(14,:,:),opts_vvel)

  ; Overlay contours                                                              
    overlay(panel1, cont0)
    overlay(panel1, panel1_ovr)
    overlay(panel1, panel1_ovr2)
    txid(1)     = gsn_create_text(wks, panel_strings(1), txres)
    amid(1)     = gsn_add_annotation(panel1, txid(1), amres)
    delete([/cont0/])

  ; Mean + eddies                                                                 
    opts_mbgt@vpXF       = xpos_plot(2)
    opts_mbgt@vpYF       = ypos_plot(2)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel2      = gsn_csm_contour(wks,diff_sum(:,:),opts_mbgt)
    cont0       = gsn_csm_contour(wks,diff_sum(:,:),opts_zero)
    panel2_ovr  = gsn_csm_contour(wks,mono_arr(12,:,:),opts_vrad)
    panel2_ovr2 = gsn_csm_contour(wks,mono_arr(14,:,:),opts_vvel)

  ; Overlay contours                                                              
    overlay(panel2, cont0)
    overlay(panel2, panel2_ovr)
    overlay(panel2, panel2_ovr2)
    txid(2)     = gsn_create_text(wks, panel_strings(2), txres)
    amid(2)     = gsn_add_annotation(panel2, txid(2), amres)
    delete([/cont0/])

  elseif (typ .eq. "novc_ring") then 

  ; Mean
    opts_mbgt@vpXF       = xpos_plot(0)
    opts_mbgt@vpYF       = ypos_plot(0)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel0      = gsn_csm_contour(wks,diff_mean(:,:),opts_mbgt)
    cont0       = gsn_csm_contour(wks,diff_mean(:,:),opts_zero)
    panel0_ovr  = gsn_csm_contour(wks,ring_arr(12,:,:),opts_vrad)
    panel0_ovr2 = gsn_csm_contour(wks,ring_arr(14,:,:),opts_vvel)

  ; Overlay contours
    overlay(panel0, cont0)
    overlay(panel0, panel0_ovr)
    overlay(panel0, panel0_ovr2)
    txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
    amid(0)     = gsn_add_annotation(panel0, txid(0), amres)
    delete([/cont0/])

  ; Eddies
    opts_mbgt@vpXF       = xpos_plot(1)
    opts_mbgt@vpYF       = ypos_plot(1)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel1      = gsn_csm_contour(wks,diff_eddy(:,:),opts_mbgt)
    cont0       = gsn_csm_contour(wks,diff_eddy(:,:),opts_zero)
    panel1_ovr  = gsn_csm_contour(wks,ring_arr(12,:,:),opts_vrad)
    panel1_ovr2 = gsn_csm_contour(wks,ring_arr(14,:,:),opts_vvel)

  ; Overlay contours
    overlay(panel1, cont0)
    overlay(panel1, panel1_ovr)
    overlay(panel1, panel1_ovr2)
    txid(1)     = gsn_create_text(wks, panel_strings(1), txres)
    amid(1)     = gsn_add_annotation(panel1, txid(1), amres)
    delete([/cont0/])

  ; Mean + eddies
    opts_mbgt@vpXF       = xpos_plot(2)
    opts_mbgt@vpYF       = ypos_plot(2)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel2      = gsn_csm_contour(wks,diff_sum(:,:),opts_mbgt)
    cont0       = gsn_csm_contour(wks,diff_sum(:,:),opts_zero)
    panel2_ovr  = gsn_csm_contour(wks,ring_arr(12,:,:),opts_vrad)
    panel2_ovr2 = gsn_csm_contour(wks,ring_arr(14,:,:),opts_vvel)

  ; Overlay contours
    overlay(panel2, cont0)
    overlay(panel2, panel2_ovr)
    overlay(panel2, panel2_ovr2)
    txid(2)     = gsn_create_text(wks, panel_strings(2), txres)
    amid(2)     = gsn_add_annotation(panel2, txid(2), amres)
    delete([/cont0/])

  elseif (typ .eq. "novc_r2m") then 

  ; Mean                                                                           
    opts_mbgt@vpXF       = xpos_plot(0)
    opts_mbgt@vpYF       = ypos_plot(0)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel0      = gsn_csm_contour(wks,diff_mean(:,:),opts_mbgt)
    cont0       = gsn_csm_contour(wks,diff_mean(:,:),opts_zero)
    panel0_ovr  = gsn_csm_contour(wks,r2m_arr(12,:,:),opts_vrad)
    panel0_ovr2 = gsn_csm_contour(wks,r2m_arr(14,:,:),opts_vvel)

  ; Overlay contours                                                               
    overlay(panel0, cont0)
    overlay(panel0, panel0_ovr)
    overlay(panel0, panel0_ovr2)
    txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
    amid(0)     = gsn_add_annotation(panel0, txid(0), amres)
    delete([/cont0/])

  ; Eddies                                                                         
    opts_mbgt@vpXF       = xpos_plot(1)
    opts_mbgt@vpYF       = ypos_plot(1)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel1      = gsn_csm_contour(wks,diff_eddy(:,:),opts_mbgt)
    cont0       = gsn_csm_contour(wks,diff_eddy(:,:),opts_zero)
    panel1_ovr  = gsn_csm_contour(wks,r2m_arr(12,:,:),opts_vrad)
    panel1_ovr2 = gsn_csm_contour(wks,r2m_arr(14,:,:),opts_vvel)

  ; Overlay contours                                                               
    overlay(panel1, cont0)
    overlay(panel1, panel1_ovr)
    overlay(panel1, panel1_ovr2)
    txid(1)     = gsn_create_text(wks, panel_strings(1), txres)
    amid(1)     = gsn_add_annotation(panel1, txid(1), amres)
    delete([/cont0/])

  ; Mean + eddies                                                                  
    opts_mbgt@vpXF       = xpos_plot(2)
    opts_mbgt@vpYF       = ypos_plot(2)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel2      = gsn_csm_contour(wks,diff_sum(:,:),opts_mbgt)
    cont0       = gsn_csm_contour(wks,diff_sum(:,:),opts_zero)
    panel2_ovr  = gsn_csm_contour(wks,r2m_arr(12,:,:),opts_vrad)
    panel2_ovr2 = gsn_csm_contour(wks,r2m_arr(14,:,:),opts_vvel)

  ; Overlay contours                                                               
    overlay(panel2, cont0)
    overlay(panel2, panel2_ovr)
    overlay(panel2, panel2_ovr2)
    txid(2)     = gsn_create_text(wks, panel_strings(2), txres)
    amid(2)     = gsn_add_annotation(panel2, txid(2), amres)
    delete([/cont0/])

  elseif (typ .eq. "novc_mono") then 

  ; Mean                                                                           
    opts_mbgt@vpXF       = xpos_plot(0)
    opts_mbgt@vpYF       = ypos_plot(0)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel0      = gsn_csm_contour(wks,diff_mean(:,:),opts_mbgt)
    cont0       = gsn_csm_contour(wks,diff_mean(:,:),opts_zero)
    panel0_ovr  = gsn_csm_contour(wks,mono_arr(12,:,:),opts_vrad)
    panel0_ovr2 = gsn_csm_contour(wks,mono_arr(14,:,:),opts_vvel)

  ; Overlay contours                                                               
    overlay(panel0, cont0)
    overlay(panel0, panel0_ovr)
    overlay(panel0, panel0_ovr2)
    txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
    amid(0)     = gsn_add_annotation(panel0, txid(0), amres)
    delete([/cont0/])

  ; Eddies                                                                         
    opts_mbgt@vpXF       = xpos_plot(1)
    opts_mbgt@vpYF       = ypos_plot(1)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel1      = gsn_csm_contour(wks,diff_eddy(:,:),opts_mbgt)
    cont0       = gsn_csm_contour(wks,diff_eddy(:,:),opts_zero)
    panel1_ovr  = gsn_csm_contour(wks,mono_arr(12,:,:),opts_vrad)
    panel1_ovr2 = gsn_csm_contour(wks,mono_arr(14,:,:),opts_vvel)

  ; Overlay contours                                                               
    overlay(panel1, cont0)
    overlay(panel1, panel1_ovr)
    overlay(panel1, panel1_ovr2)
    txid(1)     = gsn_create_text(wks, panel_strings(1), txres)
    amid(1)     = gsn_add_annotation(panel1, txid(1), amres)
    delete([/cont0/])

  ; Mean + eddies                                                                  
    opts_mbgt@vpXF       = xpos_plot(2)
    opts_mbgt@vpYF       = ypos_plot(2)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel2      = gsn_csm_contour(wks,diff_sum(:,:),opts_mbgt)
    cont0       = gsn_csm_contour(wks,diff_sum(:,:),opts_zero)
    panel2_ovr  = gsn_csm_contour(wks,mono_arr(12,:,:),opts_vrad)
    panel2_ovr2 = gsn_csm_contour(wks,mono_arr(14,:,:),opts_vvel)

  ; Overlay contours                                                               
    overlay(panel2, cont0)
    overlay(panel2, panel2_ovr)
    overlay(panel2, panel2_ovr2)
    txid(2)     = gsn_create_text(wks, panel_strings(2), txres)
    amid(2)     = gsn_add_annotation(panel2, txid(2), amres)
    delete([/cont0/])

  elseif (typ .eq. "novc_m2r") then 

  ; Mean                                                                           
    opts_mbgt@vpXF       = xpos_plot(0)
    opts_mbgt@vpYF       = ypos_plot(0)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel0      = gsn_csm_contour(wks,diff_mean(:,:),opts_mbgt)
    cont0       = gsn_csm_contour(wks,diff_mean(:,:),opts_zero)
    panel0_ovr  = gsn_csm_contour(wks,m2r_arr(12,:,:),opts_vrad)
    panel0_ovr2 = gsn_csm_contour(wks,m2r_arr(14,:,:),opts_vvel)

  ; Overlay contours                                                               
    overlay(panel0, cont0)
    overlay(panel0, panel0_ovr)
    overlay(panel0, panel0_ovr2)
    txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
    amid(0)     = gsn_add_annotation(panel0, txid(0), amres)
    delete([/cont0/])

  ; Eddies                                                                         
    opts_mbgt@vpXF       = xpos_plot(1)
    opts_mbgt@vpYF       = ypos_plot(1)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel1      = gsn_csm_contour(wks,diff_eddy(:,:),opts_mbgt)
    cont0       = gsn_csm_contour(wks,diff_eddy(:,:),opts_zero)
    panel1_ovr  = gsn_csm_contour(wks,m2r_arr(12,:,:),opts_vrad)
    panel1_ovr2 = gsn_csm_contour(wks,m2r_arr(14,:,:),opts_vvel)

  ; Overlay contours                                                               
    overlay(panel1, cont0)
    overlay(panel1, panel1_ovr)
    overlay(panel1, panel1_ovr2)
    txid(1)     = gsn_create_text(wks, panel_strings(1), txres)
    amid(1)     = gsn_add_annotation(panel1, txid(1), amres)
    delete([/cont0/])

  ; Mean + eddies                                                                  
    opts_mbgt@vpXF       = xpos_plot(2)
    opts_mbgt@vpYF       = ypos_plot(2)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel2      = gsn_csm_contour(wks,diff_sum(:,:),opts_mbgt)
    cont0       = gsn_csm_contour(wks,diff_sum(:,:),opts_zero)
    panel2_ovr  = gsn_csm_contour(wks,m2r_arr(12,:,:),opts_vrad)
    panel2_ovr2 = gsn_csm_contour(wks,m2r_arr(14,:,:),opts_vvel)

  ; Overlay contours                                                               
    overlay(panel2, cont0)
    overlay(panel2, panel2_ovr)
    overlay(panel2, panel2_ovr2)
    txid(2)     = gsn_create_text(wks, panel_strings(2), txres)
    amid(2)     = gsn_add_annotation(panel2, txid(2), amres)
    delete([/cont0/])

  elseif (typ .eq. "novc_mean") then 

  ; Mean radial influx of absolute vertical vorticity [V_mzeta]
    opts_mbgt@vpXF       = xpos_plot(0)
    opts_mbgt@vpYF       = ypos_plot(0)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel0      = gsn_csm_contour(wks,plot_arr(1,:,:),opts_mbgt)  ; [V_mzeta]             
    cont0       = gsn_csm_contour(wks,mean_arr(1,:,:),opts_zero)  ; Draw contours         

    panel0_ovr  = gsn_csm_contour(wks,mean_arr(12,:,:),opts_vrad) ; [u]                   
    panel0_ovr2 = gsn_csm_contour(wks,mean_arr(14,:,:),opts_vvel) ; [w]                   

  ; Overlay contours                                                                      
    overlay(panel0, cont0)
    overlay(panel0, panel0_ovr)
    overlay(panel0, panel0_ovr2)
    txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
    amid(0)     = gsn_add_annotation(panel0, txid(0), amres)
    delete([/cont0/])

  ; Mean vertical advection of mean tangential momentum [V_mv]                            
    opts_mbgt@vpXF       = xpos_plot(1)
    opts_mbgt@vpYF       = ypos_plot(1)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel1      = gsn_csm_contour(wks,plot_arr(2,:,:),opts_mbgt)  ; [V_mv]                
    cont0       = gsn_csm_contour(wks,mean_arr(2,:,:),opts_zero)  ; Draw contours         

    panel1_ovr  = gsn_csm_contour(wks,mean_arr(12,:,:),opts_vrad) ; [U]                   
    panel1_ovr2 = gsn_csm_contour(wks,mean_arr(14,:,:),opts_vvel) ; [w]                   

  ; Overlay contours                                                                      
    overlay(panel1, cont0)
    overlay(panel1, panel1_ovr)
    overlay(panel1, panel1_ovr2)
    txid(1)     = gsn_create_text(wks, panel_strings(1), txres)
    amid(1)     = gsn_add_annotation(panel1, txid(1), amres)
    delete([/cont0/])

  ; Eddy radial vorticity flux [V_ezeta]                                                  
    opts_mbgt@vpXF       = xpos_plot(2)
    opts_mbgt@vpYF       = ypos_plot(2)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel2      = gsn_csm_contour(wks,plot_arr(3,:,:),opts_mbgt)  ; [V_ezeta]             
    cont0       = gsn_csm_contour(wks,mean_arr(3,:,:),opts_zero)  ; Draw contours         

    panel2_ovr  = gsn_csm_contour(wks,mean_arr(12,:,:),opts_vrad) ; [U]                   
    panel2_ovr2 = gsn_csm_contour(wks,mean_arr(14,:,:),opts_vvel) ; [w]                   

  ; Overlay contours                                                                      
    overlay(panel2, cont0)
    overlay(panel2, panel2_ovr)
    overlay(panel2, panel2_ovr2)
    txid(2)     = gsn_create_text(wks, panel_strings(2), txres)
    amid(2)     = gsn_add_annotation(panel2, txid(2), amres)
    delete([/cont0/])

  ; Vertical advection of eddy tangential momentum [V_ev]                                 
    opts_mbgt@vpXF       = xpos_plot(3)
    opts_mbgt@vpYF       = ypos_plot(3)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel3      = gsn_csm_contour(wks,plot_arr(4,:,:),opts_mbgt)  ; [V_ev]                
    cont0       = gsn_csm_contour(wks,mean_arr(4,:,:),opts_zero)  ; Draw contours         

    panel3_ovr  = gsn_csm_contour(wks,mean_arr(12,:,:),opts_vrad) ; [U]                   
    panel3_ovr2 = gsn_csm_contour(wks,mean_arr(14,:,:),opts_vvel) ; [w]                   

  ; Overlay contours                                                                      
    overlay(panel3, cont0)
    overlay(panel3, panel3_ovr)
    overlay(panel3, panel3_ovr2)
    txid(3)     = gsn_create_text(wks, panel_strings(3), txres)
    amid(3)     = gsn_add_annotation(panel3, txid(3), amres)
    delete([/cont0/])

; Radial wind overlaid with absolute angular momentum
  elseif (typ .eq. "rad") then 

   if (row .eq. 2) then 

   ; Ring-like phase 
     opts_mbgt@vpXF       = xpos_plot(0)
     opts_mbgt@vpYF       = ypos_plot(0)
     opts_mbgt@vpWidthF   = xsize
     opts_mbgt@vpHeightF  = ysize

   ; Produce filled contours of radial wind 
     panel0      = gsn_csm_contour(wks,ring_arr(12,:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,ring_arr(12,:,:),opts_zero) 

   ; Overlay AAM and vertical velocity
     panel0_ovr  = gsn_csm_contour(wks,ring_arr(15,:,:),opts_aam)  
     panel0_ovr2 = gsn_csm_contour(wks,ring_arr(14,:,:),opts_vvel) 

   ; Overlay contours
     overlay(panel0, cont0)
     overlay(panel0, panel0_ovr)
     overlay(panel0, panel0_ovr2)
     txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
     amid(0)     = gsn_add_annotation(panel0, txid(0), amres)
     delete([/cont0/])

   ; Ring-like to monopole transition
     opts_mbgt@vpXF       = xpos_plot(1)
     opts_mbgt@vpYF       = ypos_plot(1)
     opts_mbgt@vpWidthF   = xsize
     opts_mbgt@vpHeightF  = ysize
 
   ; Produce filled contours of radial wind 
     panel1      = gsn_csm_contour(wks,r2m_arr(12,:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,r2m_arr(12,:,:),opts_zero)

   ; Overlay AAM and vertical velocity
     panel1_ovr  = gsn_csm_contour(wks,r2m_arr(15,:,:),opts_aam)
     panel1_ovr2 = gsn_csm_contour(wks,r2m_arr(14,:,:),opts_vvel)

   ; Overlay contours
     overlay(panel1, cont0)
     overlay(panel1, panel1_ovr)
     overlay(panel1, panel1_ovr2)
     txid(1)     = gsn_create_text(wks, panel_strings(1), txres)
     amid(1)     = gsn_add_annotation(panel1, txid(1), amres)
     delete([/cont0/])

   ; Monopole phase
     opts_mbgt@vpXF       = xpos_plot(2)
     opts_mbgt@vpYF       = ypos_plot(2)
     opts_mbgt@vpWidthF   = xsize
     opts_mbgt@vpHeightF  = ysize

   ; Produce filled contours of radial wind 
     panel2      = gsn_csm_contour(wks,mono_arr(12,:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,mono_arr(12,:,:),opts_zero)

   ; Overlay AAM and vertical velocity 
     panel2_ovr  = gsn_csm_contour(wks,mono_arr(15,:,:),opts_aam) 
     panel2_ovr2 = gsn_csm_contour(wks,mono_arr(14,:,:),opts_vvel)

   ; Overlay contours
     overlay(panel2, cont0)
     overlay(panel2, panel2_ovr)
     overlay(panel2, panel2_ovr2)
     txid(2)     = gsn_create_text(wks, panel_strings(2), txres)
     amid(2)     = gsn_add_annotation(panel2, txid(2), amres)
     delete([/cont0/])

   ; Monopole to ring-like transition
     opts_mbgt@vpXF       = xpos_plot(3)
     opts_mbgt@vpYF       = ypos_plot(3)
     opts_mbgt@vpWidthF   = xsize
     opts_mbgt@vpHeightF  = ysize

   ; Produce filled contours of radial wind 
     panel3      = gsn_csm_contour(wks,m2r_arr(12,:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,m2r_arr(12,:,:),opts_zero)

   ; Overlay AAM and vertical velocity 
     panel3_ovr  = gsn_csm_contour(wks,m2r_arr(15,:,:),opts_aam)
     panel3_ovr2 = gsn_csm_contour(wks,m2r_arr(14,:,:),opts_vvel)

   ; Overlay contours                                                                        
     overlay(panel3, cont0)
     overlay(panel3, panel3_ovr)
     overlay(panel3, panel3_ovr2)
     txid(3)     = gsn_create_text(wks, panel_strings(3), txres)
     amid(3)     = gsn_add_annotation(panel3, txid(3), amres)
     delete([/cont0/])

 ; Difference plot between two VC phases (single panel)
   else

   ; Ring-like phase                                                                      
     opts_mbgt@vpXF        = 0.05
     opts_mbgt@vpYF        = 0.95
     opts_mbgt@gsnMaximize = True 

   ; Calculate difference between monopole and ring-like phases

   ; Radial wind 
     vrad_diff   = mono_arr(12,:,:) - ring_arr(12,:,:)
     copy_VarCoords(ring_arr(12,:,:), vrad_diff)
     copy_VarAtts(ring_arr(12,:,:), vrad_diff)

   ; Absolute angular momentum 
     aam_diff    = mono_arr(15,:,:) - ring_arr(15,:,:)
     copy_VarCoords(ring_arr(15,:,:), aam_diff)
     copy_VarAtts(ring_arr(15,:,:), aam_diff)

   ; Vertical velocity
     vvel_diff   = mono_arr(14,:,:) - ring_arr(14,:,:)
     copy_VarCoords(ring_arr(14,:,:), vvel_diff)
     copy_VarAtts(ring_arr(14,:,:), vvel_diff)

   ; Produce filled contours of radial wind                                               
     panel0      = gsn_csm_contour(wks,vrad_diff,opts_mbgt)
     cont0       = gsn_csm_contour(wks,vrad_diff,opts_zero)

   ; Overlay AAM and vertical velocity                                                    
     panel0_ovr  = gsn_csm_contour(wks,aam_diff,opts_aam)
     panel0_ovr2 = gsn_csm_contour(wks,vvel_diff,opts_vvel)

   ; Overlay contours                                                                     
     overlay(panel0, cont0)
     overlay(panel0, panel0_ovr)
     overlay(panel0, panel0_ovr2)
     txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
     amid(0)     = gsn_add_annotation(panel0, txid(0), amres)
     delete([/cont0/])

   end if 

; Local tangential wind tendency (LHS)
  elseif (typ .eq. "lhs") then 

    r2m_arr(0,:,:)  = r2m_arr(0,:,:) * 12
    mono_arr(0,:,:) = mono_arr(0,:,:) * 12
    m2r_arr(0,:,:)  = m2r_arr(0,:,:) * 12

  ; Ring-like phase                                                                           
    opts_mbgt@vpXF       = xpos_plot(0)
    opts_mbgt@vpYF       = ypos_plot(0)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel0      = gsn_csm_contour(wks,ring_arr(0,:,:),opts_mbgt)  ; [V_mzeta]                
    cont0       = gsn_csm_contour(wks,ring_arr(0,:,:),opts_zero)  ; Draw contours           

  ; Overlay AAM and vertical velocity                                                         
    panel0_ovr  = gsn_csm_contour(wks,ring_arr(15,:,:),opts_aam)   ; [AAM]                    
    panel0_ovr2 = gsn_csm_contour(wks,ring_arr(14,:,:),opts_vvel)  ; [w]                      

  ; Overlay contours                                                                          
    overlay(panel0, cont0)
    overlay(panel0, panel0_ovr)
    overlay(panel0, panel0_ovr2)
    txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
    amid(0)     = gsn_add_annotation(panel0, txid(0), amres)
    delete([/cont0/])

  ; Ring-like to monopole transition                                                          
    opts_mbgt@vpXF       = xpos_plot(1)
    opts_mbgt@vpYF       = ypos_plot(1)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel1      = gsn_csm_contour(wks,r2m_arr(0,:,:),opts_mbgt)  ; [V_mv]                    
    cont0       = gsn_csm_contour(wks,r2m_arr(0,:,:),opts_zero)  ; Draw contours             

  ; Overlay AAM and vertical velocity                                                         
    panel1_ovr  = gsn_csm_contour(wks,r2m_arr(15,:,:),opts_aam)   ; [AAM]                     
    panel1_ovr2 = gsn_csm_contour(wks,r2m_arr(14,:,:),opts_vvel)  ; [w]                       

  ; Overlay contours                                                                          
    overlay(panel1, cont0)
    overlay(panel1, panel1_ovr)
    overlay(panel1, panel1_ovr2)
    txid(1)     = gsn_create_text(wks, panel_strings(1), txres)
    amid(1)     = gsn_add_annotation(panel1, txid(1), amres)
    delete([/cont0/])

  ; Monopole phase                                                                            
    opts_mbgt@vpXF       = xpos_plot(2)
    opts_mbgt@vpYF       = ypos_plot(2)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel2      = gsn_csm_contour(wks,mono_arr(0,:,:),opts_mbgt)  ; [V_ezeta]               
    cont0       = gsn_csm_contour(wks,mono_arr(0,:,:),opts_zero)  ; Draw contours           

  ; Overlay AAM and vertical velocity                                                         
    panel2_ovr  = gsn_csm_contour(wks,mono_arr(15,:,:),opts_aam)   ; [AAM]                    
    panel2_ovr2 = gsn_csm_contour(wks,mono_arr(14,:,:),opts_vvel)  ; [w]                     

  ; Overlay contours                                                                          
    overlay(panel2, cont0)
    overlay(panel2, panel2_ovr)
    overlay(panel2, panel2_ovr2)
    txid(2)     = gsn_create_text(wks, panel_strings(2), txres)
    amid(2)     = gsn_add_annotation(panel2, txid(2), amres)
    delete([/cont0/])

  ; Monopole to ring-like transition                                                          
    opts_mbgt@vpXF       = xpos_plot(3)
    opts_mbgt@vpYF       = ypos_plot(3)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel3      = gsn_csm_contour(wks,m2r_arr(0,:,:),opts_mbgt)  ; [V_ev]                    
    cont0       = gsn_csm_contour(wks,m2r_arr(0,:,:),opts_zero)  ; Draw contours             

  ; Overlay AAM and vertical velocity                                                         
    panel3_ovr  = gsn_csm_contour(wks,m2r_arr(15,:,:),opts_aam)   ; [AAM]                     
    panel3_ovr2 = gsn_csm_contour(wks,m2r_arr(14,:,:),opts_vvel)  ; [w]                       

  ; Overlay contours                                                                          
    overlay(panel3, cont0)
    overlay(panel3, panel3_ovr)
    overlay(panel3, panel3_ovr2)
    txid(3)     = gsn_create_text(wks, panel_strings(3), txres)
    amid(3)     = gsn_add_annotation(panel3, txid(3), amres)
    delete([/cont0/])


; Which options for 'typ' does this part cover?  
  else

  ; Mean radial influx of absolute vertical vorticity [V_mzeta]
    opts_mbgt@vpXF       = xpos_plot(0)
    opts_mbgt@vpYF       = ypos_plot(0)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel0      = gsn_csm_contour(wks,plot_arr(1,:,:),opts_mbgt)  ; [V_mzeta]
    cont0       = gsn_csm_contour(wks,plot_arr(1,:,:),opts_zero)  ; Draw contours

    panel0_ovr  = gsn_csm_contour(wks,plot_arr(12,:,:),opts_vrad) ; [u]
    panel0_ovr2 = gsn_csm_contour(wks,plot_arr(14,:,:),opts_vvel) ; [w]

  ; Overlay contours
    overlay(panel0, cont0)
    overlay(panel0, panel0_ovr)
    overlay(panel0, panel0_ovr2)
    txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
    amid(0)     = gsn_add_annotation(panel0, txid(0), amres)
    delete([/cont0/])

  ; Mean vertical advection of mean tangential momentum [V_mv]
    opts_mbgt@vpXF       = xpos_plot(1)
    opts_mbgt@vpYF       = ypos_plot(1)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel1      = gsn_csm_contour(wks,plot_arr(2,:,:),opts_mbgt)  ; [V_mv]
    cont0       = gsn_csm_contour(wks,plot_arr(2,:,:),opts_zero)  ; Draw contours

    panel1_ovr  = gsn_csm_contour(wks,plot_arr(12,:,:),opts_vrad) ; [U]
    panel1_ovr2 = gsn_csm_contour(wks,plot_arr(14,:,:),opts_vvel) ; [w]

  ; Overlay contours
    overlay(panel1, cont0)
    overlay(panel1, panel1_ovr)
    overlay(panel1, panel1_ovr2)
    txid(1)     = gsn_create_text(wks, panel_strings(1), txres)
    amid(1)     = gsn_add_annotation(panel1, txid(1), amres)
    delete([/cont0/])

  ; Eddy radial vorticity flux [V_ezeta]
    opts_mbgt@vpXF       = xpos_plot(2)
    opts_mbgt@vpYF       = ypos_plot(2)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel2      = gsn_csm_contour(wks,plot_arr(3,:,:),opts_mbgt)  ; [V_ezeta]
    cont0       = gsn_csm_contour(wks,plot_arr(3,:,:),opts_zero)  ; Draw contours

    panel2_ovr  = gsn_csm_contour(wks,plot_arr(12,:,:),opts_vrad) ; [U]
    panel2_ovr2 = gsn_csm_contour(wks,plot_arr(14,:,:),opts_vvel) ; [w]

  ; Overlay contours
    overlay(panel2, cont0)
    overlay(panel2, panel2_ovr)
    overlay(panel2, panel2_ovr2)
    txid(2)     = gsn_create_text(wks, panel_strings(2), txres)
    amid(2)     = gsn_add_annotation(panel2, txid(2), amres)
    delete([/cont0/])

  ; Vertical advection of eddy tangential momentum [V_ev]
    opts_mbgt@vpXF       = xpos_plot(3)
    opts_mbgt@vpYF       = ypos_plot(3)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel3      = gsn_csm_contour(wks,plot_arr(4,:,:),opts_mbgt)  ; [V_ev]
    cont0       = gsn_csm_contour(wks,plot_arr(4,:,:),opts_zero)  ; Draw contours

    panel3_ovr  = gsn_csm_contour(wks,plot_arr(12,:,:),opts_vrad) ; [U]
    panel3_ovr2 = gsn_csm_contour(wks,plot_arr(14,:,:),opts_vvel) ; [w]

  ; Overlay contours
    overlay(panel3, cont0)
    overlay(panel3, panel3_ovr)
    overlay(panel3, panel3_ovr2)
    txid(3)     = gsn_create_text(wks, panel_strings(3), txres)
    amid(3)     = gsn_add_annotation(panel3, txid(3), amres)
    delete([/cont0/])

  end if 

;==============================================
; Overlay mean position of RMW (if required)
;==============================================

  ; Plotting resources 
    nres0                    = True
    nres0@gsMarkerIndex      = 12
    nres0@gsMarkerSizeF      = 10.0
    nres0@gsMarkerColor      = "black"
    nres0@gsMarkerThicknessF = 3.0

; If overlaying on certain type of panel plot (8 or 12 panels)
  if (typ .eq. "mean" .or. typ .eq. "ring" .or. \
      typ .eq. "r2m" .or. typ .eq. "mono" .or. typ .eq. "m2r" .or. \
      typ .eq. "v" .or. typ .eq. "slp" .or. typ .eq. "novc") then

 ; Customise to change number of panels 
   mean_marker  = new(8, graphic)

   if (paper .eq. 1 .or. paper .eq. 0)

  ; Overlay marker onto existing tangential wind plot 
    if (typ .eq. "ring") then 
     mean_marker(0)  = gsn_add_polymarker(wks, panel_vrad, ring_r_rmw, ring_z_rmw, nres0)
     mean_marker(1)  = gsn_add_polymarker(wks, panel_vtan, ring_r_rmw, ring_z_rmw, nres0)
     mean_marker(2)  = gsn_add_polymarker(wks, panel_vvel, ring_r_rmw, ring_z_rmw, nres0)
     mean_marker(3)  = gsn_add_polymarker(wks, panel_aam,  ring_r_rmw, ring_z_rmw, nres0)
     mean_marker(4)  = gsn_add_polymarker(wks, panel0,     ring_r_rmw, ring_z_rmw, nres0)
     mean_marker(5)  = gsn_add_polymarker(wks, panel1,     ring_r_rmw, ring_z_rmw, nres0)
     mean_marker(6)  = gsn_add_polymarker(wks, panel6,     ring_r_rmw, ring_z_rmw, nres0)
     mean_marker(7)  = gsn_add_polymarker(wks, panel7,     ring_r_rmw, ring_z_rmw, nres0)
    elseif (typ .eq. "r2m") then 
     mean_marker(0)  = gsn_add_polymarker(wks, panel_vrad, r2m_r_rmw, r2m_z_rmw, nres0)
     mean_marker(1)  = gsn_add_polymarker(wks, panel_vtan, r2m_r_rmw, r2m_z_rmw, nres0)
     mean_marker(2)  = gsn_add_polymarker(wks, panel_vvel, r2m_r_rmw, r2m_z_rmw, nres0)
     mean_marker(3)  = gsn_add_polymarker(wks, panel_aam,  r2m_r_rmw, r2m_z_rmw, nres0)
     mean_marker(4)  = gsn_add_polymarker(wks, panel0,     r2m_r_rmw, r2m_z_rmw, nres0)
     mean_marker(5)  = gsn_add_polymarker(wks, panel1,     r2m_r_rmw, r2m_z_rmw, nres0)
     mean_marker(6)  = gsn_add_polymarker(wks, panel6,     r2m_r_rmw, r2m_z_rmw, nres0)
     mean_marker(7)  = gsn_add_polymarker(wks, panel7,     r2m_r_rmw, r2m_z_rmw, nres0)
    elseif (typ .eq. "mono") then 
     mean_marker(0)  = gsn_add_polymarker(wks, panel_vrad, mono_r_rmw, mono_z_rmw, nres0)
     mean_marker(1)  = gsn_add_polymarker(wks, panel_vtan, mono_r_rmw, mono_z_rmw, nres0)
     mean_marker(2)  = gsn_add_polymarker(wks, panel_vvel, mono_r_rmw, mono_z_rmw, nres0)
     mean_marker(3)  = gsn_add_polymarker(wks, panel_aam,  mono_r_rmw, mono_z_rmw, nres0)
     mean_marker(4)  = gsn_add_polymarker(wks, panel0,     mono_r_rmw, mono_z_rmw, nres0)
     mean_marker(5)  = gsn_add_polymarker(wks, panel1,     mono_r_rmw, mono_z_rmw, nres0)
     mean_marker(6)  = gsn_add_polymarker(wks, panel6,     mono_r_rmw, mono_z_rmw, nres0)
     mean_marker(7)  = gsn_add_polymarker(wks, panel7,     mono_r_rmw, mono_z_rmw, nres0)
    elseif (typ .eq. "m2r") then 
     mean_marker(0)  = gsn_add_polymarker(wks, panel_vrad, m2r_r_rmw, m2r_z_rmw, nres0)
     mean_marker(1)  = gsn_add_polymarker(wks, panel_vtan, m2r_r_rmw, m2r_z_rmw, nres0)
     mean_marker(2)  = gsn_add_polymarker(wks, panel_vvel, m2r_r_rmw, m2r_z_rmw, nres0)
     mean_marker(3)  = gsn_add_polymarker(wks, panel_aam,  m2r_r_rmw, m2r_z_rmw, nres0)
     mean_marker(4)  = gsn_add_polymarker(wks, panel0,     m2r_r_rmw, m2r_z_rmw, nres0)
     mean_marker(5)  = gsn_add_polymarker(wks, panel1,     m2r_r_rmw, m2r_z_rmw, nres0)
     mean_marker(6)  = gsn_add_polymarker(wks, panel6,     m2r_r_rmw, m2r_z_rmw, nres0)
     mean_marker(7)  = gsn_add_polymarker(wks, panel7,     m2r_r_rmw, m2r_z_rmw, nres0)
    end if 

   end if 

; Otherwise overlay on our 4-panel radial wind plot 
  elseif (typ .eq. "rad" .or. typ .eq. "lhs") then 

   if (row .ne. 1) then 
    mean_marker    = new(4, graphic)
    mean_marker(0) = gsn_add_polymarker(wks, panel0, ring_r_rmw, ring_z_rmw, nres0)
    mean_marker(1) = gsn_add_polymarker(wks, panel1, r2m_r_rmw,  r2m_z_rmw, nres0)
    mean_marker(2) = gsn_add_polymarker(wks, panel2, mono_r_rmw, mono_z_rmw, nres0)
    mean_marker(3) = gsn_add_polymarker(wks, panel3,  m2r_r_rmw,  m2r_z_rmw, nres0)
   else
    mean_marker    = new(1, graphic)
    mean_marker(0) = gsn_add_polymarker(wks, panel0, ring_r_rmw, ring_z_rmw, nres0)
   end if 

; Difference plots of combined mean + eddy terms 
  elseif (typ .eq. "paper2_diff") then 

   mean_marker  = new(3, graphic)

   mean_marker(0)  = gsn_add_polymarker(wks, panel0, mono_r_rmw, mono_z_rmw, nres0)
   mean_marker(1)  = gsn_add_polymarker(wks, panel1, mono_r_rmw, mono_z_rmw, nres0)
   mean_marker(2)  = gsn_add_polymarker(wks, panel2, mono_r_rmw, mono_z_rmw, nres0)

; Comparison of non-VC composite and single VC phase (radial wind + AAM, or local vtan tendency)
  elseif (typ .eq. "novc_rad" .or. typ .eq. "novc_lhs") then 

   if (col .eq. 2) then 
    mean_marker    = new(2, graphic)
   else
    mean_marker    = new(3, graphic)
   end if 

   mean_marker(0)  = gsn_add_polymarker(wks, panel0, novc_r_rmw, novc_z_rmw, nres0)
   if (p0 .eq. "ring") then 
    mean_marker(1) = gsn_add_polymarker(wks, panel1, ring_r_rmw, ring_z_rmw, nres0)
   elseif (p0 .eq. "r2m") then 
    mean_marker(1) = gsn_add_polymarker(wks, panel1, r2m_r_rmw, r2m_z_rmw, nres0)   
   elseif (p0 .eq. "mono") then 
    mean_marker(1) = gsn_add_polymarker(wks, panel1, mono_r_rmw, mono_z_rmw, nres0)
   elseif (p0 .eq. "m2r") then 
    mean_marker(1) = gsn_add_polymarker(wks, panel1, m2r_r_rmw, m2r_z_rmw, nres0)
   end if 

 ; Add label for third panel (if required)
   if (col .ne. 2) then 
    if (p1 .eq. "ring") then
     mean_marker(2) = gsn_add_polymarker(wks, panel2, ring_r_rmw, ring_z_rmw, nres0)
    elseif (p1 .eq. "r2m") then
     mean_marker(2) = gsn_add_polymarker(wks, panel2, r2m_r_rmw, r2m_z_rmw, nres0)
    elseif (p1 .eq. "mono") then
     mean_marker(2) = gsn_add_polymarker(wks, panel2, mono_r_rmw, mono_z_rmw, nres0)
    elseif (p1 .eq. "m2r") then 
     mean_marker(2) = gsn_add_polymarker(wks, panel2, m2r_r_rmw, m2r_z_rmw, nres0)
    end if
   end if 

; Overlay on mean-eddy panel plot
  elseif (typ .eq. "mean_eddy") then 

   mean_marker  = new(6, graphic)

 ; Top row 
   if (p0 .eq. "ring") then 
    mean_marker(0) = gsn_add_polymarker(wks, panel0, ring_r_rmw, ring_z_rmw, nres0)
    mean_marker(1) = gsn_add_polymarker(wks, panel1, ring_r_rmw, ring_z_rmw, nres0)
    mean_marker(2) = gsn_add_polymarker(wks, panel2, ring_r_rmw, ring_z_rmw, nres0)
   elseif (p0 .eq. "r2m") then 
    mean_marker(0) = gsn_add_polymarker(wks, panel0, r2m_r_rmw, r2m_z_rmw, nres0)
    mean_marker(1) = gsn_add_polymarker(wks, panel1, r2m_r_rmw, r2m_z_rmw, nres0)
    mean_marker(2) = gsn_add_polymarker(wks, panel2, r2m_r_rmw, r2m_z_rmw, nres0)
   elseif (p0 .eq. "mono") then 
    mean_marker(0) = gsn_add_polymarker(wks, panel0, mono_r_rmw, mono_z_rmw, nres0)
    mean_marker(1) = gsn_add_polymarker(wks, panel1, mono_r_rmw, mono_z_rmw, nres0)
    mean_marker(2) = gsn_add_polymarker(wks, panel2, mono_r_rmw, mono_z_rmw, nres0)
   else 
    mean_marker(0) = gsn_add_polymarker(wks, panel0, m2r_r_rmw, m2r_z_rmw, nres0)
    mean_marker(1) = gsn_add_polymarker(wks, panel1, m2r_r_rmw, m2r_z_rmw, nres0)
    mean_marker(2) = gsn_add_polymarker(wks, panel2, m2r_r_rmw, m2r_z_rmw, nres0)
   end if 

 ; Bottom row 
   if (p1 .eq. "ring") then
    mean_marker(3) = gsn_add_polymarker(wks, panel3, ring_r_rmw, ring_z_rmw, nres0)
    mean_marker(4) = gsn_add_polymarker(wks, panel4, ring_r_rmw, ring_z_rmw, nres0)
    mean_marker(5) = gsn_add_polymarker(wks, panel5, ring_r_rmw, ring_z_rmw, nres0)
   elseif (p1 .eq. "r2m") then
    mean_marker(3) = gsn_add_polymarker(wks, panel3, r2m_r_rmw, r2m_z_rmw, nres0)
    mean_marker(4) = gsn_add_polymarker(wks, panel4, r2m_r_rmw, r2m_z_rmw, nres0)
    mean_marker(5) = gsn_add_polymarker(wks, panel5, r2m_r_rmw, r2m_z_rmw, nres0)
   elseif (p1 .eq. "mono") then
    mean_marker(3) = gsn_add_polymarker(wks, panel3, mono_r_rmw, mono_z_rmw, nres0)
    mean_marker(4) = gsn_add_polymarker(wks, panel4, mono_r_rmw, mono_z_rmw, nres0)
    mean_marker(5) = gsn_add_polymarker(wks, panel5, mono_r_rmw, mono_z_rmw, nres0)
   else
    mean_marker(3) = gsn_add_polymarker(wks, panel3, m2r_r_rmw, m2r_z_rmw, nres0)
    mean_marker(4) = gsn_add_polymarker(wks, panel4, m2r_r_rmw, m2r_z_rmw, nres0)
    mean_marker(5) = gsn_add_polymarker(wks, panel5, m2r_r_rmw, m2r_z_rmw, nres0)
   end if

  end if 

;==================================
; Add labelbar at bottom of page
;==================================

    lopt                        = True
    lopt@lbAutoManage           = False
    lopt@lbOrientation          = "horizontal"
    lopt@lbLabelAlignment       = "InteriorEdges"
    lopt@lbFillPattern          = "SolidFill"
    lopt@lbMonoFillPattern      = True
    lopt@lbBoxEndCapStyle       = "TriangleBothEnds"

  ; Define colourmap and set labelbar colours 
    gsn_define_colormap(wks,"BlueDarkRed18")
    if (typ .ne. "rad" .and. typ .ne. "novc_rad") then
     lopt@lbFillColors           = (/2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19/)
    else
     lopt@lbFillColors           = (/2,4,5,6,7,8,9,0,0,12,13,14,15,16,17,19/)
    end if 

  ; Customise labelbar (number of boxes, strings, size, etc) 
    if (typ .ne. "rad" .and. typ .ne. "novc_rad" .and. \
        typ .ne. "lhs" .and. typ .ne. "novc_lhs") then 

     if (res .eq. "5min") then 
      label_str                 = (/"-25.0","-20.0","-15.0","-10.0","-5.0","-2.0","-1.0",\
                                    "-0.5","0.0", "0.5", "1.0", "2.0", \
                                    "5.0", "10.0", "15.0", "20.0", "25.0", "30.0"/)
     else
      label_str                 = (/"-20.0","-15.0","-10.0","-5.0","-3.0","-2.0","-1.0",\
                                    "-0.5","0.0","0.5", "1.0", "2.0", "3.0", \
                                    "5.0", "10.0", "15.0", "20.0", "25.0"/)
     end if 

    elseif (typ .eq. "rad" .or. typ .eq. "novc_rad") then 

     if (res .eq. "5min") then 
      label_str                 = (/"-10.0","-5.0","-3.0","-2.0","-1.0", "-0.5", \
                                    "-0.2","0.0", "0.2", "0.5", "1.0", "2.0", "3.0", \
                                    "5.0", "10.0", "15.0"/)
     else
      label_str                 = (/"-5.0","-3.0","-2.0","-1.0","-0.5","-0.2","-0.1",\
                                    "0.0","0.1","0.2","0.5","1.0","2.0","3.0","5.0","10.0"/)
     end if 

    elseif (typ .eq. "lhs") then

      label_str                 = (/"-20.0","-15.0","-10.0","-5.0","-3.0","-2.0","-1.0",\
                                    "-0.5","0.0","0.5", "1.0", "2.0", "3.0", \
                                    "5.0", "10.0", "15.0", "20.0", "25.0"/)

    elseif (typ .eq. "novc_lhs") then 

      label_str                 = (/"-15.0","-10.0","-5.0","-3.0","-2.0","-1.0",\
                                    "-0.5","-0.2","0.0","0.2","0.5","1.0",\
                                    "2.0","3.0","5.0","10.0","15.0","20.0"/)

    end if 

    nlabels                     = dimsizes(label_str)
    lopt@lbLabelFontHeightF     = 0.010
    lopt@txFont                 = "Helvetica"
    lopt@lbTitleOn              = True
    if (typ .ne. "rad" .and. typ .ne. "lhs" .and. \
        typ .ne. "novc_rad" .and. typ .ne. "novc_lhs") then 
     lopt@lbTitleString         = "Momentum budget (m s~S~-1~N~ h~S~-1~N~)"
    elseif (typ .eq. "rad" .or. typ .eq. "novc_rad") then 
     lopt@lbTitleString         = "Radial wind (m s~S~-1~N~)"
    elseif (typ .eq. "lhs" .or. typ .eq. "novc_lhs") then 
     lopt@lbTitleString         = "Tangential wind tendency (m s~S~-1~N~ h~S~-1~N~)"
    end if 
    lopt@lbTitleFontHeightF     = 0.0125
    lopt@lbTitlePosition        = "Bottom"
    
  ; If 'paper' option is active, then change labelbar size + position accordingly 
    if (typ .eq. "mean" .or. typ .eq. "novc" .or. \
        typ .eq. "ring" .or. typ .eq. "r2m" .or. typ .eq. "mono" .or. \
        typ .eq. "m2r" .or. typ .eq. "v" .or. typ .eq. "slp") then
      if (paper .eq. 1) then 
        if (row .eq. 3) then 
         lopt@vpWidthF             = 0.65
         lopt@vpHeightF            = 0.08
         xpos_lbar                 = 0.16
         ypos_lbar                 = 0.10
        else
         lopt@vpWidthF             = 0.70
         lopt@vpHeightF            = 0.10
         xpos_lbar                 = 0.16
         ypos_lbar                 = 0.48
        end if 
      elseif (paper .eq. 0) then
        lopt@vpWidthF              = 0.75
        lopt@vpHeightF             = 0.12
        xpos_lbar                  = 0.12
        ypos_lbar                  = 0.25
      end if 
    end if 

  ; Change labelbar size + position depending on type of plot 
    if (typ .eq. "comp") then
      lopt@vpWidthF              = 0.7
      lopt@vpHeightF             = 0.1
      xpos_lbar                  = 0.16
      ypos_lbar                  = 0.22
    elseif (typ .eq. "mean_eddy") then 
      lopt@vpWidthF              = 0.75
      lopt@vpHeightF             = 0.12
      xpos_lbar                  = 0.13
      ypos_lbar                  = 0.34
    elseif (typ .eq. "paper2_diff" .or. typ .eq. "novc_mean" .or. typ .eq. "novc_ring" .or. \
            typ .eq. "novc_r2m" .or. typ .eq. "novc_mono" .or. typ .eq. "novc_m2r") then 
      lopt@vpWidthF              = 0.75
      lopt@vpHeightF             = 0.12
      xpos_lbar                  = 0.13
      ypos_lbar                  = 0.67
    elseif (typ .eq. "novc_lhs" .or. typ .eq. "novc_rad") then 
      lopt@vpWidthF              = 0.75
      lopt@vpHeightF             = 0.12
      xpos_lbar                  = 0.13
      ypos_lbar                  = 0.67
    elseif (typ .eq. "paper1_diff" .or. typ .eq. "paper0_diff") then 
      lopt@vpWidthF              = 0.7
      lopt@vpHeightF             = 0.1
      xpos_lbar                  = 0.16
      ypos_lbar                  = 0.20
    elseif (typ .eq. "rad" .or. typ .eq. "lhs") then 
     if (row .ne. 1) then 
      lopt@vpWidthF              = 0.7
      lopt@vpHeightF             = 0.1
      xpos_lbar                  = 0.16
      ypos_lbar                  = 0.20
     else
      lopt@vpWidthF              = 0.7
      lopt@vpHeightF             = 0.1
      xpos_lbar                  = 0.16
      ypos_lbar                  = 0.12
     end if 
    end if 

  ; Use routine to overlay labelbar
    gsn_labelbar_ndc(wks, nlabels, label_str, xpos_lbar, ypos_lbar, lopt)

;=================================
; Also add title at top of page 
;=================================

    if (typ .ne. "comp" .and. typ .ne. "mean_eddy" .and. paper .ne. 1 .and. paper .ne. 0) then 

     title_res                  = True
     title_str                  = plot_tt
     title_res@txFontHeightF    = 0.02
     title_res@txFont           = "Helvetica"
     xpos_title                 = 0.50
     ypos_title                 = 0.94

   ; Use routine to overlay labelbar 
     gsn_text_ndc(wks, title_str, xpos_title, ypos_title, title_res)

    end if 

  ; Finally, draw the plot with everything overlaid (call 'draw' and 'frame')
    pres        = True 
    maximize_output(wks,pres)

  ; Print output file path to screen 
    print(output)

end 