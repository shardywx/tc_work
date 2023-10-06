; Read in vertical mass flux data, and produce scatter/box plots

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

; ncl opt=\"x11\" sub=0 t0=1 mth=\"p2\" trb=0 tend=2 plt=\"all\" mlev0=17 mlev1=45 alt=1 sc=5
; box=1 lgd0=\"tr\" grp=\"8\" wt0=1 wt=0 w0=0.3 r0=0 r1=25 run0=3 vc_mflux_scatter_paper.ncl

; 'opt'   = output file format ("pdf", "x11", etc)
; 'mth'   = method for calculating storm centre
; 'trb'   = turn troubleshooting on (1) or off (0)
; 'tend'  = method for calculating centred differences (2 or 1)
; 'alt'   = use original (0) or alternative (1) VC phase data
; 'plt'   = all phases ("all"); main four VC phases ("vc"); VC/No-VC ("no")
; 'sub'   = use a subset of the data points (e.g. only positive tendencies)
; 't0'    = wind speed / MSLP tendency threshold (e.g. 0,1,2), if 'sub=1'  
; 'wt0'   = use threshold method (1) to constrain radial integration for mass flux calculations
; 'r0'    = inner radius for calculation
; 'r1'    = outer radius for calculation 
; 'run0'  = length of time (h) over which to apply running average (3,6,12)

begin

;=========================================
; Now read in the data from text files 
;=========================================

; File paths 
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

; Now create new arrays to hold all data for each VC phase
  var_ring  = new( (/rsize(0), nts_ring, 3/), "float")
  var_r2m   = new( (/rsize(1), nts_r2m, 3/), "float")
  var_mono  = new( (/rsize(2), nts_mono, 3/), "float")
  var_m2r   = new( (/rsize(3), nts_m2r, 3/), "float")

; Dimension size information for each VC phase
  dsize     = (/ rsize(0), rsize(1), \
                 rsize(2), rsize(3) /)

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
   diri_rvp  = "$ar/text/"+dat_ring(st)+"_"+ens0_ring(st)+"_sc"+sc+"_inner_core.txt"

 ; Radial vorticity profile tendency (119 values)
   diri_ten  = "$ar/text/"+dat_ring(st)+"_"+ens0_ring(st)+"_sc"+sc+"_inner_core_grad.txt"

 ; Mass flux (119 values)

 ; Upper layer 
   diri_up   = "$ar/text/flux_upper_"+dat_ring(st)+"_"+ens0_ring(st)+"_sc"+sc+\
               "_lay1_"+hy0+"_"+hy1+"_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+"_run"+run0+".txt"

 ; Lower layer 
   diri_low  = "$ar/text/flux_lower_"+dat_ring(st)+"_"+ens0_ring(st)+"_sc"+sc+\
               "_lay1_"+hy0+"_"+hy1+"_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+"_run"+run0+".txt"

 ; Radial vorticity profile (and tendency)
   rvp_all     = asciiread(diri_rvp,  ntimes_vc, "float")
   ten_all     = asciiread(diri_ten,  ntimes_vc, "float")
   flux_up     = asciiread(diri_up,   ntimes_vc, "float")
   flux_low    = asciiread(diri_low,  ntimes_vc, "float")

 ; Calculate difference in vertical mass flux between layers 
   flux_all0   = flux_up - flux_low

 ; Express difference as a percentage of the mass flux in the lower layer
   flux_all    = (flux_all0 / flux_low) * 100
   flux_all@units = "%"

 ; Finally, all data into correct arrays [radial vorticity profile + mass flux]
   npts = ( tf0_ring(st) - ts0_ring(st) )
   var_ring(st,0:npts,0) = (/flux_all( ts0_ring(st):tf0_ring(st) )/)
   var_ring(st,0:npts,1) = (/rvp_all( ts0_ring(st):tf0_ring(st) )/)
   var_ring(st,0:npts,2) = (/ten_all( ts0_ring(st):tf0_ring(st) )/)

  end do  

  delete([/diri_rvp, diri_ten, diri_low, diri_up/])

;=============================================
; Loop over ringlike to monopole transition
;=============================================

  do st = 0, dsize(1)-1

 ; Radial vorticity profile (119 values)                                                    
   diri_rvp  = "$ar/text/"+dat_r2m(st)+"_"+ens0_r2m(st)+"_sc"+sc+"_inner_core.txt"

 ; Radial vorticity profile tendency (119 values)                                           
   diri_ten  = "$ar/text/"+dat_r2m(st)+"_"+ens0_r2m(st)+"_sc"+sc+"_inner_core_grad.txt"

 ;==============
 ; Mass flux 
 ;==============

 ; Upper layer                                                                               
   diri_up   = "$ar/text/flux_upper_"+dat_r2m(st)+"_"+ens0_r2m(st)+"_sc"+sc+\
               "_lay1_"+hy0+"_"+hy1+"_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+"_run"+run0+".txt"

 ; Lower layer                                                                               
   diri_low  = "$ar/text/flux_lower_"+dat_r2m(st)+"_"+ens0_r2m(st)+"_sc"+sc+\
               "_lay1_"+hy0+"_"+hy1+"_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+"_run"+run0+".txt"

 ; Radial vorticity profile (and tendency)                                                 
   rvp_all     = asciiread(diri_rvp,  ntimes_vc, "float")
   ten_all     = asciiread(diri_ten,  ntimes_vc, "float")
   flux_up     = asciiread(diri_up,   ntimes_vc, "float")
   flux_low    = asciiread(diri_low,  ntimes_vc, "float")

 ; Calculate difference in vertical mass flux between layers                                
   flux_all0   = flux_up - flux_low

 ; Express difference as a percentage of the mass flux in the lower layer                   
   flux_all    = (flux_all0 / flux_low) * 100
   flux_all@units = "%"

 ; Finally, all data into correct arrays [radial vorticity profile]
   npts = ( tf0_r2m(st) - ts0_r2m(st) )
   var_r2m(st,0:npts,0) = (/flux_all( ts0_r2m(st):tf0_r2m(st) )/)
   var_r2m(st,0:npts,1) = (/rvp_all( ts0_r2m(st):tf0_r2m(st) )/)
   var_r2m(st,0:npts,2) = (/ten_all( ts0_r2m(st):tf0_r2m(st) )/)

  end do
  delete([/diri_rvp, diri_ten, diri_low, diri_up/])

;============================
; Loop over monopole phase 
;============================

  do st = 0, dsize(2)-1

 ; Radial vorticity profile (119 values)                                                   
   diri_rvp  = "$ar/text/"+dat_mono(st)+"_"+ens0_mono(st)+"_sc"+sc+"_inner_core.txt"

 ; Radial vorticity profile tendency (119 values)                                          
   diri_ten  = "$ar/text/"+dat_mono(st)+"_"+ens0_mono(st)+"_sc"+sc+"_inner_core_grad.txt"

 ; Mass flux 

 ; Upper layer
   diri_up   = "$ar/text/flux_upper_"+dat_mono(st)+"_"+ens0_mono(st)+"_sc"+sc+\
               "_lay1_"+hy0+"_"+hy1+"_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+"_run"+run0+".txt"

 ; Lower layer
   diri_low  = "$ar/text/flux_lower_"+dat_mono(st)+"_"+ens0_mono(st)+"_sc"+sc+\
               "_lay1_"+hy0+"_"+hy1+"_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+"_run"+run0+".txt"

 ; Radial vorticity profile (and tendency)                                                
   rvp_all     = asciiread(diri_rvp,  ntimes_vc, "float")
   ten_all     = asciiread(diri_ten,  ntimes_vc, "float")
   flux_up     = asciiread(diri_up,   ntimes_vc, "float")
   flux_low    = asciiread(diri_low,  ntimes_vc, "float")

 ; Calculate difference in vertical mass flux between layers                                
   flux_all0   = flux_up - flux_low

 ; Express difference as a percentage of the mass flux in the lower layer                   
   flux_all    = (flux_all0 / flux_low) * 100
   flux_all@units = "%"

 ; Finally, all data into correct arrays [radial vorticity profile]
   npts = ( tf0_mono(st) - ts0_mono(st) )
   var_mono(st,0:npts,0) = (/flux_all( ts0_mono(st):tf0_mono(st) )/)
   var_mono(st,0:npts,1) = (/rvp_all( ts0_mono(st):tf0_mono(st) )/)
   var_mono(st,0:npts,2) = (/ten_all( ts0_mono(st):tf0_mono(st) )/)

  end do
  delete([/diri_rvp, diri_ten, diri_low, diri_up/])

;=============================================
; Loop over monopole to ringlike transition                                        
;=============================================

  do st = 0, dsize(3)-1

 ; Radial vorticity profile (119 values)                                                    
   diri_rvp  = "$ar/text/"+dat_m2r(st)+"_"+ens0_m2r(st)+"_sc"+sc+"_inner_core.txt"

 ; Radial vorticity profile tendency (119 values)                                           
   diri_ten  = "$ar/text/"+dat_m2r(st)+"_"+ens0_m2r(st)+"_sc"+sc+"_inner_core_grad.txt"

 ; Mass flux 

 ; Upper layer
   diri_up   = "$ar/text/flux_upper_"+dat_m2r(st)+"_"+ens0_m2r(st)+"_sc"+sc+\
               "_lay1_"+hy0+"_"+hy1+"_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+"_run"+run0+".txt"

 ; Lower layer
   diri_low  = "$ar/text/flux_lower_"+dat_m2r(st)+"_"+ens0_m2r(st)+"_sc"+sc+\
               "_lay1_"+hy0+"_"+hy1+"_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+"_run"+run0+".txt"

 ; Radial vorticity profile (and tendency)                                                 
   rvp_all     = asciiread(diri_rvp,  ntimes_vc, "float")
   ten_all     = asciiread(diri_ten,  ntimes_vc, "float")
   flux_up     = asciiread(diri_up,   ntimes_vc, "float")
   flux_low    = asciiread(diri_low,  ntimes_vc, "float")

 ; Calculate difference in vertical mass flux between layers                                
   flux_all0   = flux_up - flux_low

 ; Express difference as a percentage of the mass flux in the lower layer                   
   flux_all    = (flux_all0 / flux_low) * 100
   flux_all@units = "%"

 ; Finally, all data into correct arrays [radial vorticity profile]
   npts = ( tf0_m2r(st) - ts0_m2r(st) )
   var_m2r(st,0:npts,0) = (/flux_all( ts0_m2r(st):tf0_m2r(st) )/)
   var_m2r(st,0:npts,1) = (/rvp_all( ts0_m2r(st):tf0_m2r(st) )/)
   var_m2r(st,0:npts,2) = (/ten_all( ts0_m2r(st):tf0_m2r(st) )/)

  end do
  delete([/diri_rvp, diri_ten, diri_low, diri_up/])

;================================================================
; Reshape 2D arrays (sim * time) to 1D (time) before plotting 
;================================================================

 ; Ringlike phase 
   ring_flux= ndtooned(var_ring(:,:,0) )
   ring_rvp = ndtooned(var_ring(:,:,2) )
   r_size   = dimsizes(ring_rvp)

 ; Ringlike to monopole transition
   r2m_flux = ndtooned(var_r2m(:,:,0) )
   r2m_rvp  = ndtooned(var_r2m(:,:,2) )
   r2m_size = dimsizes(r2m_rvp)

 ; Monopole phase 
   mono_flux=ndtooned(var_mono(:,:,0) )
   mono_rvp = ndtooned(var_mono(:,:,2) )
   m_size   = dimsizes(mono_rvp)

 ; Monopole to ringlike transition
   m2r_flux = ndtooned(var_m2r(:,:,0) )
   m2r_rvp  = ndtooned(var_m2r(:,:,2) )
   m2r_size = dimsizes(m2r_rvp)

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
   all_flux = new( (/all_size/), "float")
   all_rvp  = new( (/all_size/), "float")

 ; Fill each of these arrays with the values from each VC phase 
   all_flux(0:s0-1)   = ring_flux
   all_flux(s0:s1-1)  = r2m_flux
   all_flux(s1:s2-1)  = mono_flux
   all_flux(s2:s3-1)  = m2r_flux

   all_rvp(0:s0-1)    = ring_rvp
   all_rvp(s0:s1-1)   = r2m_rvp
   all_rvp(s1:s2-1)   = mono_rvp
   all_rvp(s2:s3-1)   = m2r_rvp

 ; Add metadata so we can colour-code the markers by VC phase
   all_flux!0         = "vc_phase"
   all_rvp!0          = "vc_phase"

 ; Change units before plotting 
   all_rvp            = all_rvp * 100
   all_rvp@units      = "10~S~-2~N~ h~S~-1~N~"

;===========================================================
; If desired, remove selected data points before plotting 
;===========================================================

   if (sub .eq. 1) then 
    all_flux = where(all_flux .lt. t0, all_flux@_FillValue, all_flux)
   end if    

;================================================
; Open workstation and define output file path
;================================================

   output = "$nep/nepartak/images/vc_mflux_scatter_paper_group"+grp+"_sc"+sc+\
             "_sub"+t0+"_"+r0+"km_"+r1+"km_wt"+wt
   wks    = gsn_open_wks(opt, output)

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
   if (scat .eq. 1) then  
    vc_plot_fl0   = gsn_csm_xy(wks, all_rvp(0:s0-1), all_flux(0:s0-1), vc_opts)
    fl1           = gsn_add_polymarker(wks,vc_plot_fl0,all_rvp(s0:s1-1),all_flux(s0:s1-1),vc_opts1)
    fl2           = gsn_add_polymarker(wks,vc_plot_fl0,all_rvp(s1:s2-1),all_flux(s1:s2-1),vc_opts2)
    fl3           = gsn_add_polymarker(wks,vc_plot_fl0,all_rvp(s2:s3-1),all_flux(s2:s3-1),vc_opts3)
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
   if (scat .eq. 1) then 
    annotate                          = gsn_add_annotation(vc_plot_fl0, lgd, am_opts)

  ; Now switch to Cartesian grid in both panels
    c_axis_uv  = cartesian_axis(wks, vc_plot_fl0)
 
  ; Finally, draw selected panel(s)
    draw(vc_plot_fl0)
    frame(wks)
   end if 

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
  plot_int0 = all_flux(:)
  
; Calculate stats for each VC phase 
  ptb_ring  = stat_dispersion(all_flux(0:s0-1), opts)
  ptb_r2m   = stat_dispersion(all_flux(s0:s1-1), opts)
  ptb_mono  = stat_dispersion(all_flux(s1:s2-1), opts)
  ptb_m2r   = stat_dispersion(all_flux(s2:s3-1), opts)

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

; Print some useful information to screen 
  print("Printing mean values from each phase: "+box_y(:,2) )

;========================                                              
; Produce the box plot                                                 
;========================     

; Set output file path
;  out_box = "$nep/nepartak/images/vc_mflux_boxplot_paper_group"\
;            +grp+"_sc"+sc+"_sub"+t0+"_"+r0+"km_"+r1+"km"
  out_box = "./mflux_boxplot"
  wks     = gsn_open_wks(opt, out_box)

; Create set of integer values denoting box label positions                         
  x0  = ispan(1,4,1)

; General box plot resources                                                        
  pres = True
  pres@tmXBLabels      = (/"Ring","R to M","Mono","M to R"/)
  pres@gsnMaximize     = True
  pres@tiYAxisFontHeightF = 0.0225

; Y-axis title (+ position)
  pres@tiYAxisString   = "Mass flux difference (%)"
  pres@tiYAxisOffsetXF = -0.01

; Y-axis limits 
  pres@trYMaxF         = 60.0
  pres@trYMinF         = -90.0

; Box resources                                                                     
  bres = True
  bres@boxWidth  = 0.5
  bres@boxColors = (/"royalblue1", "purple3", "red3", "orange3"/)

; Line resources (06/12/2019 --> make boxes thicker)
  lres = True
  lres@gsLineThicknessF = 3.0

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

; Print output file path to screen 
  print(out_box)

 end if 

end 
