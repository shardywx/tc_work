; Script to plot key circulation budget terms for different VC phases

; Run using:

; ncl opt=\"x11\" sc=5 full=9 flux=1 tc_cbgt_cxz_plot.ncl 

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"

load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/nc_levs.ncl"
load "$sam/ncl_func/st_centre.ncl"
load "$sam/ncl_func/setup_cyl.ncl"
load "$sam/ncl_func/ring_mono.ncl"

begin

;==============================
; Prelash (reading in files) 
;==============================

; Also retrieve pressure-level array from 'pb' stream
  diri       = "/nfs/a319/earshar/02T12/em11/"
  fili_pr    = "20160702T1200Z_NPTK_4p4_L80_ra1t_em11"
  fili_b1    = systemfunc("cd "+diri+" ; ls "+fili_pr+"_pb.nc")
  fili_b     = diri+"/"+fili_b1

; Retrieve number of levels [1]; array of levels [0] 
  nLEVS_b    = nc_levs(fili_b)
  num_l      = nLEVS_b[1]
  levs_b     = nLEVS_b[0]
  levs_b@units = "hPa"

; Symmetric phase
  str_sym    = asciiread("$ar/text/cbgt_cxz_sym_full"+full+"_sc"+sc+"_str.txt",(/num_l/),"float")
  hadv_sym   = asciiread("$ar/text/cbgt_cxz_sym_full"+full+"_sc"+sc+"_hadv.txt",(/num_l/),"float")
  tilt_sym   = asciiread("$ar/text/cbgt_cxz_sym_full"+full+"_sc"+sc+"_tilt.txt",(/num_l/),"float")

; Symmetric to asymmetric transition
  str_s2a    = asciiread("$ar/text/cbgt_cxz_s2a_full"+full+"_sc"+sc+"_str.txt",(/num_l/),"float")
  hadv_s2a   = asciiread("$ar/text/cbgt_cxz_s2a_full"+full+"_sc"+sc+"_hadv.txt",(/num_l/),"float")
  tilt_s2a   = asciiread("$ar/text/cbgt_cxz_s2a_full"+full+"_sc"+sc+"_tilt.txt",(/num_l/),"float")

; Asymmetric phase
  str_asym   = asciiread("$ar/text/cbgt_cxz_asym_full"+full+"_sc"+sc+"_str.txt",(/num_l/),"float")
  hadv_asym  = asciiread("$ar/text/cbgt_cxz_asym_full"+full+"_sc"+sc+"_hadv.txt",(/num_l/),"float")
  tilt_asym  = asciiread("$ar/text/cbgt_cxz_asym_full"+full+"_sc"+sc+"_tilt.txt",(/num_l/),"float")

; Asymmetric to symmetric transition
  str_a2s    = asciiread("$ar/text/cbgt_cxz_a2s_full"+full+"_sc"+sc+"_str.txt",(/num_l/),"float")
  hadv_a2s   = asciiread("$ar/text/cbgt_cxz_a2s_full"+full+"_sc"+sc+"_hadv.txt",(/num_l/),"float")
  tilt_a2s   = asciiread("$ar/text/cbgt_cxz_a2s_full"+full+"_sc"+sc+"_tilt.txt",(/num_l/),"float")

; Also add stretching and horizontal advection (advective flux)
  flux_sym   = str_sym + hadv_sym
  flux_s2a   = str_s2a + hadv_s2a
  flux_asym  = str_asym + hadv_asym
  flux_a2s   = str_a2s + hadv_a2s

;================================== 
; Create array for each VC phase   
;==================================

; Symmetric
  if (flux .eq. 0) then 
   arr_sym     = (/str_sym, hadv_sym/) ; , tilt_sym/)
   arr_sym!1   = "lev"
  else
   arr_sym     = (/flux_sym/)
   arr_sym!0   = "lev"
  end if 
  arr_sym&lev  = levs_b

; Symmetric to asymmetric transition
  if (flux .eq. 0) then 
   arr_s2a     = (/str_s2a, hadv_s2a/) ; , tilt_s2a/)
   arr_s2a!1   = "lev"
  else
   arr_s2a     = (/flux_s2a/)
   arr_s2a!0   = "lev"
  end if 
  arr_s2a&lev  = levs_b

; Asymmetric phase
  if (flux .eq. 0) then 
   arr_asym    = (/str_asym, hadv_asym/) ; , tilt_asym/)
   arr_asym!1  = "lev"
  else
   arr_asym    = (/flux_asym/)
   arr_asym!0  = "lev"
  end if 
  arr_asym&lev = levs_b

; Asymmetric to symmetric transition
  if (flux .eq. 0) then 
   arr_a2s     = (/str_a2s, hadv_a2s/) ; , tilt_a2s/)
   arr_a2s!1   = "lev"
  else
   arr_a2s     = (/flux_a2s/)
   arr_a2s!0   = "lev"
  end if 
  arr_a2s&lev  = levs_b

;=============================
; Output file type and path
;=============================

  if (flux .eq. 0) then 
   output      = "$sam/nepartak/images/cbgt_cxz_plot_full"+full+"_sc"+sc
  else
   output      = "$sam/nepartak/images/cbgt_cxz_plot_flux_full"+full+"_sc"+sc
  end if 
  wks          = gsn_open_wks(opt, output)

;====================== 
; Plotting resources 
;======================

; Y-axis resources 
  opts_cbgt                              = True
  opts_cbgt@tiYAxisString                = "Pressure (hPa)"
  opts_cbgt@tmYROn                       = False
  opts_cbgt@tmYRMode                     = "Explicit"

  opts_cbgt@trYMinF                      = 150
  opts_cbgt@trYMaxF                      = 900
  opts_cbgt@trYReverse                   = "True"

  if (flux .eq. 0) then 
   opts_cbgt@trXMinF                     = -3
   opts_cbgt@trXMaxF                     = 3
  else
   opts_cbgt@trXMinF                     = -5
   opts_cbgt@trXMaxF                     = 5
  end if 
  opts_cbgt@tiXAxisOn                    = False
  opts_cbgt@tiXAxisString                = "Circulation budget (m s~S~-1~N~ h~S~-1~N~)"
  opts_cbgt@tmXTOn                       = False
  opts_cbgt@tmXBMode                     = "Explicit"  

  opts_cbgt@gsnDraw                      = False
  opts_cbgt@gsnFrame                     = False
  opts_cbgt@tmXMajorGridLineDashPattern  = 2
  opts_cbgt@tmYMajorGridLineDashPattern  = 2

; Line plot resources for circulation budget terms 
  if (flux .eq. 0) then 
   opts_cbgt@xyLineColors                = (/"royalblue1","orange2"/)
   opts_cbgt@xyDashPatterns              = (/0, 0/)
   opts_cbgt@xyLineThicknesses           = (/2.0, 2.0/)
  else
   opts_cbgt@xyLineColors                = (/"royalblue1"/)
   opts_cbgt@xyDashPatterns              = (/0/)
   opts_cbgt@xyLineThicknesses           = (/2.0/)
  end if 

; Panel plot label resources 
  txid                                   = new(4,graphic)
  amid                                   = new(4,graphic)
  panel_strings                          = (/ "a) SYM", "b) S2A", "c) ASYM", "d) A2S"/)

  txres                                  = True
  txres@txPerimOn                        = True
  txres@txFontHeightF                    = 0.0125
  txres@txBackgroundFillColor            = "White"

  amres                                  = True
  amres@amParallelPosF                   = -0.48
  amres@amOrthogonalPosF                 = -0.48
  amres@amJust                           = "TopLeft"

; Coordinates of panels 
  xf = (/0.12, 0.58, 0.12, 0.58/)
  yf = (/0.90, 0.90, 0.49, 0.49/)
  wf = 0.32
  hf = 0.28

; Plot each panel 
  panel    = new(4,graphic)

; Panel 1 (SYM)
  opts_cbgt1           = opts_cbgt
  opts_cbgt1@vpXF      = xf(0)
  opts_cbgt1@vpYF      = yf(0)
  opts_cbgt1@vpWidthF  = wf
  opts_cbgt1@vpHeightF = hf

  panel(0) = gsn_csm_xy(wks, arr_sym, levs_b, opts_cbgt1)
  txid(0)  = gsn_create_text(wks, panel_strings(0), txres)
  amid(0)  = gsn_add_annotation(panel(0), txid(0), amres)

; Panel 2 (S2A)
  opts_cbgt2           = opts_cbgt
  opts_cbgt2@vpXF      = xf(1)
  opts_cbgt2@vpYF      = yf(1)
  opts_cbgt2@vpWidthF  = wf
  opts_cbgt2@vpHeightF = hf

  panel(1) = gsn_csm_xy(wks, arr_s2a, levs_b, opts_cbgt2)
  txid(1)  = gsn_create_text(wks, panel_strings(1), txres)
  amid(1)  = gsn_add_annotation(panel(1), txid(1), amres)

; Panel 3 (ASYM)
  opts_cbgt3           = opts_cbgt
  opts_cbgt3@vpXF      = xf(2)
  opts_cbgt3@vpYF      = yf(2)
  opts_cbgt3@vpWidthF  = wf
  opts_cbgt3@vpHeightF = hf

  panel(2) = gsn_csm_xy(wks, arr_asym, levs_b, opts_cbgt3)
  txid(2)  = gsn_create_text(wks, panel_strings(2), txres)
  amid(2)  = gsn_add_annotation(panel(2), txid(2), amres)

; Panel 4 (A2S)
  opts_cbgt4           = opts_cbgt
  opts_cbgt4@vpXF      = xf(3)
  opts_cbgt4@vpYF      = yf(3)
  opts_cbgt4@vpWidthF  = wf
  opts_cbgt4@vpHeightF = hf

  panel(3) = gsn_csm_xy(wks, arr_a2s, levs_b, opts_cbgt4)
  txid(3)  = gsn_create_text(wks, panel_strings(3), txres)
  amid(3)  = gsn_add_annotation(panel(3), txid(3), amres)

;==============================================
; Overlay vertical line at x=0 on each panel 
;============================================== 

; New graphical object to hold lines for each panel 
  lines = new( dimsizes(panel), graphic)

; Line resources 
; See https://www.ncl.ucar.edu/Document/Graphics/Images/dashpatterns.png for dash patterns
  opts_l = True
  opts_l@gsLineThicknessF  = 2.0
  opts_l@gsLineDashPattern = 1 
  opts_l@gsLineColor       = "black" 

;; Overlay vertical line onto each panel
;  do n = 0, dimsizes(panel)-1
;    lines(n) = gsn_add_polyline(wks,panel(n),(/0.0,opts_cbgt@trYMinF/),\
;                                (/0.0,opts_cbgt@trYMaxF/),opts_l)
;  end do 

  overlay_0  = gsn_add_polyline(wks,panel(0),(/0.0,opts_cbgt@trYMinF/),\
                                (/0.0,opts_cbgt@trYMaxF/),opts_l)
  overlay_1  = gsn_add_polyline(wks,panel(1),(/0.0,opts_cbgt@trYMinF/),\
                                (/0.0,opts_cbgt@trYMaxF/),opts_l)
  overlay_2  = gsn_add_polyline(wks,panel(2),(/0.0,opts_cbgt@trYMinF/),\
                                (/0.0,opts_cbgt@trYMaxF/),opts_l)
  overlay_3  = gsn_add_polyline(wks,panel(3),(/0.0,opts_cbgt@trYMinF/),\
                                (/0.0,opts_cbgt@trYMaxF/),opts_l)

  draw(panel)
  frame(wks)

end 