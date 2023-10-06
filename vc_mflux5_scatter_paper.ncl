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
; box=1 lgd0=\"tr\" grp=\"8\" wt0=1 wt=0 w0=0.3 r0=0 r1=25 run0=3 vc_mflux5_scatter_paper.ncl

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

;========================================= 
; Now read in the data from text files 
;========================================= 

; File paths 
  diri_vc   = "$ar/text/flux5_"

; Upper layer
  up_ring   = diri_vc+"upper_ring_sc"+sc+"_"+hy0+"_"+hy1+"_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+".txt"
  up_r2m    = diri_vc+"upper_r2m_sc"+sc+"_"+hy0+"_"+hy1+"_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+".txt"
  up_mono   = diri_vc+"upper_mono_sc"+sc+"_"+hy0+"_"+hy1+"_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+".txt"
  up_m2r    = diri_vc+"upper_m2r_sc"+sc+"_"+hy0+"_"+hy1+"_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+".txt"

; Lower layer
  lo_ring   = diri_vc+"lower_ring_sc"+sc+"_"+hy0+"_"+hy1+"_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+".txt"
  lo_r2m    = diri_vc+"lower_r2m_sc"+sc+"_"+hy0+"_"+hy1+"_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+".txt"
  lo_mono   = diri_vc+"lower_mono_sc"+sc+"_"+hy0+"_"+hy1+"_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+".txt"
  lo_m2r    = diri_vc+"lower_m2r_sc"+sc+"_"+hy0+"_"+hy1+"_w"+w0+"_"+r0+"km_"+r1+"km_wt"+wt+".txt"

; Group all data in single [2 * 4] array 
  fili_vc   = (/ (/up_ring, up_r2m, up_mono, up_m2r/), \
                 (/lo_ring, lo_r2m, lo_mono, lo_m2r/) /)

; Get file size information for each VC phase using 'systemfunc'
  rsize0 = new(4,"string")
  do il  = 0, 3
   rsize0(il) = systemfunc("wc -l < "+fili_vc(0,il) )
  end do
  rsize  = toint(rsize0)

; Read in data (upper layer)
  upflux_ring = asciiread(fili_vc(0,0),rsize(0),"float")
  upflux_r2m  = asciiread(fili_vc(0,1),rsize(1),"float")
  upflux_mono = asciiread(fili_vc(0,2),rsize(2),"float")
  upflux_m2r  = asciiread(fili_vc(0,3),rsize(3),"float")

; Read in data (lower layer)
  loflux_ring = asciiread(fili_vc(1,0),rsize(0),"float")
  loflux_r2m  = asciiread(fili_vc(1,1),rsize(1),"float")
  loflux_mono = asciiread(fili_vc(1,2),rsize(2),"float")
  loflux_m2r  = asciiread(fili_vc(1,3),rsize(3),"float")

; Now create new arrays to hold all data for each VC phase
  var_ring    = new( (/rsize(0), 3/), "float")
  var_r2m     = new( (/rsize(1), 3/), "float")
  var_mono    = new( (/rsize(2), 3/), "float")
  var_m2r     = new( (/rsize(3), 3/), "float")

; Add metadata to arrays above 
  var_ring!0    = "time"
  var_ring!1    = "diag"
  var_ring&diag = (/"Mass flux", "RVP","RVP tendency"/)

  var_r2m!0     = "time"
  var_r2m!1     = "diag"
  var_r2m&diag  = (/"Mass flux", "RVP","RVP tendency"/)

  var_mono!0    = "time"
  var_mono!1    = "diag"
  var_mono&diag = (/"Mass flux", "RVP","RVP tendency"/)

  var_m2r!0     = "time"
  var_m2r!1     = "diag"
  var_m2r&diag  = (/"Mass flux", "RVP","RVP tendency"/)

;=======================================================
; Calculate mass flux difference for ring-like phase
;=======================================================

;; Radial vorticity profile (119 values)
;  diri_rvp  = "$ar/text/"+dat_ring(st)+"_"+ens0_ring(st)+"_sc"+sc+"_inner_core.txt"

;; Radial vorticity profile tendency (119 values)
;  diri_ten  = "$ar/text/"+dat_ring(st)+"_"+ens0_ring(st)+"_sc"+sc+"_inner_core_grad.txt"

;; Radial vorticity profile (and tendency) 
;  rvp_all     = asciiread(diri_rvp,  ntimes_vc, "float")
;  ten_all     = asciiread(diri_ten,  ntimes_vc, "float")

; Calculate difference in vertical mass flux between layers 
  diff_ring   = upflux_ring - loflux_ring

; Express difference as a percentage of the mass flux in the lower layer
  flux_plot   = (diff_ring / loflux_ring) * 100
  flux_plot@units = "%"

; Finally, all data into correct arrays [radial vorticity profile + mass flux]
  var_ring(:,0) = (/flux_plot/)
  delete(flux_plot)
;  var_ring(:,1) = (/rvp_all( ts0_ring(st):tf0_ring(st) )/)
;  var_ring(:,2) = (/ten_all( ts0_ring(st):tf0_ring(st) )/)

;====================================
; Ring-like to monopole transition
;====================================

; Calculate difference in vertical mass flux between layers
  diff_r2m   = upflux_r2m - loflux_r2m

; Express difference as a percentage of the mass flux in the lower layer
  flux_plot   = (diff_r2m / loflux_r2m) * 100
  flux_plot@units = "%"

; Finally, all data into correct arrays [radial vorticity profile + mass flux] 
  var_r2m(:,0) = (/flux_plot/)
  delete(flux_plot)

;===================
; Monopole phase
;===================

; Calculate difference in vertical mass flux between layers                             
  diff_mono   = upflux_mono - loflux_mono

; Express difference as a percentage of the mass flux in the lower layer                
  flux_plot   = (diff_mono / loflux_mono) * 100
  flux_plot@units = "%"

; Finally, all data into correct arrays [radial vorticity profile + mass flux]          
  var_mono(:,0) = (/flux_plot/)
  delete(flux_plot)

;=====================================
; Monopole to ring-like transition
;=====================================

; Calculate difference in vertical mass flux between layers                             
  diff_m2r   = upflux_m2r - loflux_m2r

; Express difference as a percentage of the mass flux in the lower layer                
  flux_plot   = (diff_m2r / loflux_m2r) * 100
  flux_plot@units = "%"

; Finally, all data into correct arrays [radial vorticity profile + mass flux]          
  var_m2r(:,0) = (/flux_plot/)
  delete(flux_plot)

;================================================
; Open workstation and define output file path
;================================================

   output = "$nep/nepartak/images/vc_mflux5_scatter_paper_group"+grp+"_sc"+sc+\
             "_sub"+t0+"_"+r0+"km_"+r1+"km_wt"+wt
   wks    = gsn_open_wks(opt, output)

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

; Calculate stats for each VC phase 
  ptb_ring  = stat_dispersion(var_ring(:,0), opts)
  ptb_r2m   = stat_dispersion(var_r2m(:,0), opts)
  ptb_mono  = stat_dispersion(var_mono(:,0), opts)
  ptb_m2r   = stat_dispersion(var_m2r(:,0), opts)

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
  out_box = "./vc_mflux5_boxplot"
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
  pres@trYMaxF         = 90.0
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
