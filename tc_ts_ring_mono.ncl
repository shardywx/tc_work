; Script to investigate vacillation cycles in developing TCs
; Plot the 'symmetry ratio' of the inner-core vorticity vs pressure/windspeed tendency

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"

; ncl opt=\"x11\" dat=\"02T12\" ens0=\"em11\" mlev0=13 ts=30 tf=90 lay=1
; tend0=1 lgd0=\"br\" mlev1=25 cn0=\"slp\" wrt=0 mth=\"p2\" 
; sc=5 rvp=0 run=1 run0=3 vmax=1.0 vmin=0.5 std=1 paper=1 tc_ts_ring_mono.ncl

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'ens0'  = ensemble simulation (em00-em11)
; 'opt'   = output file format ("pdf" or "x11")
; 'run'	  = apply running average (1) or leave alone (0)
; 'ts'	  = start time for plot
; 'tf'	  = end time for plot 
; 'mlev0' = chosen model level for analysis (13, 16, 22, etc)
; 'tend0' = calculate windspeed/MSLP tendency using t2-t0 (2) or t1-t0 (1)
; 'vmax'  = maximum standard deviation for final plot 
; 'vmin'  = minimum standard deviation for final plot
; 'lgd0'  = place legend in top right ("tr") or top left ("tl") of final plot
; 'lay'	  = read in layer-averaged data (1) or single-layer data (0)
; 'mlev1' = upper model level if 'lay=1' (25, 28, 34, etc)
; 'wrt'   = only write out MSLP and vtan tendency (1) or plot as normal (0)
; 'rvp'   = write out the radial vorticity profile, and then quit (1)
; 'std'   = use standard deviation as metric (1) rather than relative vorticity ratio

; 'pc' stream --> 'it'   [121 times]
; 'centre'    --> 'it-1' [120 times]
; 'pd' stream --> 'it-2' [119 times]

begin

;======================================
; Read in the data from netCDF files
;======================================

; Get model level information from 'pd' stream, for output file string
  diri_f  = "/nfs/a319/earshar/"+dat+"/"+ens0
  fili_pr = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0
  fili_d1 = systemfunc("cd "+diri_f+" ; ls "+fili_pr+"_pd*.nc")
  fili_d  = diri_f+"/"+fili_d1

  if ( .not. fileexists (fili_d) ) then
   print("Input	file missing! Exiting script...")
   exit()
  end if

; Retrieve model level array ('hybD') from 'pd' stream
  d	  = addfile(fili_d,"r")	
  hybD	  = d->hybrid_ht(:)

; EDIT -- why is there one less time (118) in em08, 03T00?
; time    = d->t
  size 	  = 119 ; dimsizes(time)

; Create height string (e.g. 980m) for output file based upon chosen model level (e.g. 13) 
  hy0	  = sprintf("%0.0f", hybD(mlev0))

  if (lay .eq. 1) then 
   hy1    = sprintf("%0.0f", hybD(mlev1))
   print("Working on layer average: "+hy0+" m to "+hy1+" m")
  else
   print("Working on height: "+hy0+" m")
  end if 

; Create new arrays 
  plotR   = new( (/1,119/), float)
  if (paper .eq. 1) then 
   plotL  = new( (/3,119/), float)
  else
   plotL  = new( (/2,119/), float)
  end if 

; Print out relevant information to screen about storm tracking
  sc_arr      = (/"relative vorticity",\                          ; sc = 0 
                  "pressure",\                                    ; sc = 1 
                  "smoothed pressure",\                           ; sc = 2   
                  "extra smoothed pressure",\                     ; sc = 3 
                  "vort_min during symmetric phase",\             ; sc = 4 
                  "windspeed minimum",\                           ; sc = 5 
                  "geopotential"/)                                ; sc = 6 

  if (sc .lt. 3) then 
   print("Wrong option chosen - exiting!")
   exit()
  end if                   

  print("Using "+sc_arr(sc)+" to calculate storm centre position...")

;====================================
; Read in the data from	text files
;====================================

; Directory path
  diri = "$ar/text/"

; Define part of output file string (relative vorticity vs PV)
; EDIT -- remove references to this line 
  vort_str = "vort"

; 10/10/2019 --> Renamed 'structure parameter' as 'inner-core radial vorticity profile'
; 119 values
; 02/10/2019 --> EDIT FROM HERE: calculate gradients at every time interval
  sym_arr         = asciiread(diri+"paper/cp_sym_sc"+sc+"a_"+dat+"_"+ens0+".txt", (/size/), "float")

; Also read in relative vorticity standard deviation
  std_arr         = asciiread(diri+"cp_std_sc"+sc+"_"+dat+"_"+ens0+".txt", (/size/), "float")
  std_mean        = avg(std_arr)
  std_max         = max(std_arr)

; Create output file string ('std' vs 'ring')
  if (std .eq. 1) then 
   out_str = "std"
  else
   out_str = "ring"
  end if 

; Read radial vorticity profile into array
; < 1 --> ringlike (symmetric) 
; = 1 --> monopole (asymmetric) 
  if (std .eq. 1) then 
   plotR(0,:)     = std_arr(:) / std_mean
  else
   plotR(0,:)     = sym_arr(:)
  end if 

; Calculate standard deviation of radial vorticity profile
  std_dev         = dim_stddev_n(plotR(0,ts:tf), 0)
  sd_round        = sprintf("%0.2f",std_dev)
  sd_string       = "SD = "+sd_round
  print(sd_string+" ("+dat+", "+ens0+")" )

; 15/10/2019 --> write out these values to text files for additional analysis 
;  if (rvp .eq. 1) then 
;   diri_sd = diri+"rvp_sd_"+dat+"_"+ens0+".txt"
;   asciiwrite(diri_sd, sd_round)
;   exit()
;  end if 

;=====================================================================
; Also calculate standard deviation using data from all simulations 
;=====================================================================

; Size of array to calculate SD using all data points 
; [61 times * 46 sims] -- 2 simulations missing 
  num_t = 61
  num_s = 48
  num_l = num_t * num_s

; Create array of size 'num_l' [T+30 to T+90 for all sims]
  sym_arr_all     = new( (/num_s, num_t/), float) 

; Initialisation date array
  dat_arr = (/"02T12","03T00","03T12","04T00"/)

; Simulation array
  ens_arr = (/"em00","em01","em02","em03","em04","em05","em06","em07","em08","em09","em10","em11"/)

; Troubleshooting to avoid later errors
;  if (ts .ne. 30 .or. tf .ne. 90) then 
;   print("Wrong input file size - change 'ts' or 'tf'! Exiting...")
;   exit()
;  end if 

; Read in all data
; Fill elements with missing values if input file doesn't exist
  ii = 0 
  do id  = 0, dimsizes(dat_arr)-1
   do is = 0, dimsizes(ens_arr)-1

     fili               = diri+"paper/cp_sym_sc"+sc+"a_"+dat_arr(id)+"_"+ens_arr(is)+".txt"

     if (.not. fileexists (fili) ) then 
      sym_all           = new( (/size/), "float")
      sym_arr_all(ii,:) = sym_all(ts:tf)
     else
      sym_all           = asciiread(fili, (/size/), "float")
      sym_arr_all(ii,:) = sym_all(ts:tf)
     end if 

     ii = ii + 1

   end do 
  end do 

; Calculate standard deviation
  std_dev_all     = stddev(sym_arr_all(:,:) )
  sd_round_all    = sprintf("%0.2f",std_dev_all)
  sd_string_all   = "SD = "+sd_round_all
  print(sd_string_all+" (all simulations)" )

; Write out standard deviation to text file (15/10/2019)
  if (rvp .eq. 1) then
   diri_sd = diri+"rvp_sd_all_"+dat+"_"+ens0+".txt"
   asciiwrite(diri_sd, sd_round_all)
   exit()
  end if

;==================================
; Read in additional diagnostics 
;==================================

; TANGENTIAL WIND AND TENDENCY
; Read in maximum tangential wind and tendency from text files 
; 121 values (T+0 to T+120)
; Only one model height level (calculated from maximum mean vtan, which exists only on one level)
; 18/10/2019 --> rerun this script using single-level wind maxima (e.g. v500, v600, ...)
  vtan_arr0       = asciiread(diri+"ml_"+dat+"_"+ens0+"_vtan_max_"+mth+".txt",(/121/),"float")
  vtan_tend0      = asciiread(diri+"ml_"+dat+"_"+ens0+"_vtan_tend_"+mth+".txt",(/121/),"float")

; MEAN SEA LEVEL PRESSURE 
; Also read in MSLP data from text files (not calculated using same tracking method, but fine)
; 120 values (T+1 to T+120)
  mslp_arr0       = asciiread(diri+"cp_new_"+dat+"_"+ens0+"_slp.txt",(/120/),"float")

; Make all input arrays the same size (119 values)
; 10/10/2019 --> all arrays start from T+2, which should be taken into account when plotting
  vtan_arr        = vtan_arr0(2::)
  vtan_tend       = vtan_tend0(2::)
  mslp_arr        = mslp_arr0(1::)

; WINDSPEED -- calculate y-axis limits using range of data in text file
  xvt             = max(vtan_arr(ts-1:tf-1))+5
  nvt             = min(vtan_arr(ts-1:tf-1))-5

; MSLP -- calculate y-axis limits using range of data in text file 
  xslp 	          = max(mslp_arr(ts-1:tf-1))+10
  nslp 		  = min(mslp_arr(ts-1:tf-1))-10

; Calculate MSLP tendency from 'mslp_arr'
  mslp_tend       = mslp_arr
  do it = 0, size-1

   tp1 = min((/it+1,size-1/))	
   tm1 = max((/it-1,0/))

   if (tend0 .eq. 2) then 
    mslp_tend(it) = (mslp_arr(tp1) - mslp_arr(tm1)) / 2
   elseif (tend0 .eq. 1) then 
    mslp_tend(it) = (mslp_arr(it) - mslp_arr(tm1)) / 1
   end if 

  end do 

; Read MSLP and windspeed tendencies into 'plotL' (either 2 or 3 variables)
  if (paper .eq. 1) then 
   plotL(0,:)      = mslp_tend(:)
   plotL(1,:)      = vtan_tend(:)
   plotL(2,:)      = vtan_arr(:)
  else
   plotL(0,:)      = mslp_tend(:)
   plotL(1,:)      = vtan_tend(:)
  end if 

  plotR@_FillValue = -999
  plotL@_FillValue = -999

  psize0    	   = dimsizes(plotL)

;=========================================
; Sort out time array for x-axis labels
;=========================================

  time		   = ispan(0,size-1,1)
  times	    	   = dimsizes(time)

;======================================
; Create an x-y time series plot
;======================================

; Options for plotting radial vorticity profile (ringlike vs monopole)
  p0_opts		       	       = True
  p0_opts@gsnMaximize          	       = False
  p0_opts@vpXF                         = 0.20
  p0_opts@vpYF                         = 0.36
  p0_opts@vpWidthF                     = 0.55
  p0_opts@vpHeightF                    = 0.22
  p0_opts@tmYROn	       	       = False	       ; Turn off tick marks
  p0_opts@tmXTOn	       	       = False	       ; Turn off tick marks

; Set y-axis limits
  if (std .ne. 1) then 
   p0_opts@trYMaxF	       	       = vmax
   p0_opts@trYMinF	       	       = vmin
  else
   p0_opts@trYMaxF                     = 1.5
   p0_opts@trYMinF                     = 0.5   
  end if 

; Account for input arrays starting at T+2, not T+0
  p0_opts@trXMinF	       	       = ts             ; Set x-axis minimum
  p0_opts@trXMaxF	       	       = tf             ; Set x-axis maximum
 
  p0_opts@tmXBMode		       = "Explicit"
  p0_opts@tmXBValues	       	       = ispan(p0_opts@trXMinF, p0_opts@trXMaxF, 12)
  p0_opts@tmXBLabels                   = "T+"+tostring(p0_opts@tmXBValues)

; Main title string and font size
  if (paper .eq. 1) then 
   p0_opts@tiMainString      	       = ""
  else
   p0_opts@tiMainString                = "Simulation "+ens0+" ("+dat+")"
  end if 
  p0_opts@tiMainFontHeightF            = 0.0225

; X-axis (title string; title font size; label font height)
  if (paper .eq. 1) then
   p0_opts@tiXAxisString       	       = ""
  else
   p0_opts@tiXAxisString               = "Time"
  end if 
  p0_opts@tiXAxisFontHeightF   	       = 0.0150
  p0_opts@tmXBLabelFontHeightF 	       = 0.0125

; Y-axis (title string; title font size)
  if (std .eq. 1) then 
   p0_opts@tiYAxisString               = "<-- Symmetric     Asymmetric -->"
   p0_opts@tiYAxisFontHeightF          = 0.0125
  else
   p0_opts@tiYAxisString               = "<-- Ringlike        Monopole -->"
   p0_opts@tiYAxisFontHeightF          = 0.0150
  end if 

; Dashed gridlines 
  p0_opts@tmXMajorGridLineDashPattern  = 2
  p0_opts@tmYMajorGridLineDashPattern  = 2
  p0_opts@gsnDraw		       = False
  p0_opts@gsnFrame		       = False
  p0_opts@gsnPaperOrientation	       = "Portrait"

  p0_opts@xyMarkLineMode	       = "MarkLines"
  p0_opts@xyLineColors		       = "red3"
  p0_opts@xyLineThicknesses	       = 3.5
  p0_opts@xyMarkers		       = 16
  p0_opts@xyMarkerColors      	       = "red3"
  p0_opts@xyMarkerSizes		       = 5.0
  p0_opts@xyMarkerThicknesses          = 3.0

; Resources for second panel 
  p1_opts                              = p0_opts 
  p1_opts@vpXF                         = 0.20
  p1_opts@vpYF                         = 0.66
  p1_opts@vpWidthF                     = 0.55
  p1_opts@vpHeightF                    = 0.22
  p1_opts@tiMainString                 = ""

; Resources for third panel 
  p2_opts                              = p0_opts
  p2_opts@vpXF                         = 0.20
  p2_opts@vpYF                         = 0.96
  p2_opts@vpWidthF                     = 0.55
  p2_opts@vpHeightF                    = 0.22
  p2_opts@tiMainString                 = ""

; MSLP tendency resources 
  var1_opts                            = True
  var1_opts@tmXBLabelFontHeightF       = 0.01
  var1_opts@trYMaxF                    = 6.0
  var1_opts@trYMinF                    = -6.0
  var1_opts@tiMainString               = ""
  var1_opts@tiYAxisString              = "Pressure tendency (hPa h~S~-1~N~)"
  var1_opts@tiYAxisFontHeightF         = 0.0125
  var1_opts@xyLineColor                = "royalblue1"
  var1_opts@xyLineThicknessF           = 3.5
  var1_opts@xyDashPatterns             = 0.0
  var1_opts@xyMarkLineMode             = "MarkLines"
  var1_opts@xyMarker                   = 9
  var1_opts@xyMarkerColor              = "royalblue1"
  var1_opts@xyMarkerSizeF              = 0.006
  var1_opts@xyMarkerThicknessF         = 1.0

; Windspeed tendency resources 
  var0_opts                            = var1_opts
  var0_opts@tiMainString               = ""
  var0_opts@tiYAxisString              = "Wind speed tendency (m s~S~-1~N~ h~S~-1~N~)"
  var0_opts@xyLineColor                = "royalblue1"
  var0_opts@xyMarkerColor              = "royalblue1"
  var0_opts@trYMaxF                    = 4.0
  var0_opts@trYMinF                    = -4.0

; Wind speed resources 
  var2_opts                            = var1_opts
  var2_opts@tiMainString               = ""
  var2_opts@tiYAxisString              = "Wind speed (m s~S~-1~N~)"
  var2_opts@xyLineColor                = "black"
  var2_opts@xyMarkerColor              = "black"
  var2_opts@trYMaxF                    = 70.0
  var2_opts@trYMinF                    = 30.0

; Options for plotting MSLP/tangential wind tendency
  tend_optsL   	       		       = True
  tend_optsL@gsLineColor               = "royalblue1"
  tend_optsL@gsLineThicknessF          = 1.5
 
  tend_optsM			       = True
  tend_optsM@gsMarkerIndex             = 7
  tend_optsM@gsMarkerColor             = "royalblue1"
  tend_optsM@gsMarkerSizeF             = 10.0
  tend_optsM@gsMarkerThicknessF        = 3.0

; Apply running average 
  if (run .eq. 1) then 
  
   plotR(0,:) = runave_n_Wrap(plotR(0,:), run0, 0, 0)
   plotL(0,:) = runave_n_Wrap(plotL(0,:), 6, 0, 0)
   plotL(1,:) = runave_n_Wrap(plotL(1,:), 6, 0, 0)
   plotL(2,:) = runave_n_Wrap(plotL(2,:), 6, 0, 0)

  end if 

; Either write out MSLP and tangential wind tendency to text files           
; 02/10/2019 --> output file format = "inner_core.txt"
  if (wrt .eq. 1) then

   diri_out       = "$ar/text/"+dat+"_"+ens0+"_sc"+sc+"_"

   slp_tend_out   = diri_out+"mslp_tend.txt"
   vtan_tend_out  = diri_out+"vtan_tend.txt"
   inner_core_out = diri_out+"inner_core.txt"

 ; 'Structure parameter' and MSLP and windspeed tendencies
   asciiwrite(slp_tend_out, mslp_tend)
   asciiwrite(vtan_tend_out, vtan_tend)
   asciiwrite(inner_core_out, plotR(0,:) )

; Or plot time series of asymmetry and pressure/tangential wind tendency
  else 

;===============================================
; Output file path and final plotting options
;===============================================

; Output file path and type 
  if (run .eq. 1) then 

   if (paper .eq. 1) then 
    output = "$nep/nepartak/images/vc/time_series_paper_"+out_str+"_sc"+sc+"_lay_"+hy0+\
     	      "m_"+hy1+"m_"+dat+"_"+ens0+"_ave"+run0
   else
    output = "$nep/nepartak/images/vc/time_series_"+out_str+"_sc"+sc+"_lay_"+hy0+\
              "m_"+hy1+"m_"+dat+"_"+ens0+"_ave"+run0
   end if 

  else

   if (paper .eq. 1) then 
    output = "$nep/nepartak/images/vc/time_series_paper_"+out_str+"_sc"+sc+"_lay_"+hy0+\
    	     "m_"+hy1+"m_"+dat+"_"+ens0
   else
    output = "$nep/nepartak/images/vc/time_series_"+out_str+"_sc"+sc+"_lay_"+hy0+\
             "m_"+hy1+"m_"+dat+"_"+ens0
   end if 

  end if

; Open workstation for plotting 
  wks = gsn_open_wks(opt,output)

; Wind speed + inner-core structure
  res_spd = gsn_csm_xy2(wks, time, plotR(0,:), plotL(2,:), p2_opts, var2_opts)

; MSLP tendency + inner-core structure
  res_slp = gsn_csm_xy2(wks, time, plotR(0,:), plotL(0,:), p0_opts, var1_opts)

; Windspeed tendency + inner-core structure
  res_ten = gsn_csm_xy2(wks, time, plotR(0,:), plotL(1,:), p1_opts, var0_opts)

; Create graphical object to overlay horizontal line onto plot 
  line    = new(1,graphic)

; Line resources 
  opts_l  = True
  opts_l@gsLineThicknessF  = 2.0       ; Line thickness
  opts_l@gsLineDashPattern = 1         ; Dash pattern 
  opts_l@gsLineColor       = "black"   ; Line colour 

; Overlay horizontal line through y=0
; EDIT -- second to last option ( avg(vmax + vmin) )
  line    = gsn_add_polyline(wks, res_slp, (/p0_opts@trXMinF, p0_opts@trXMaxF/), \
                             (/0.75, 0.75/), opts_l)

  line1   = gsn_add_polyline(wks, res_ten, (/p1_opts@trXMinF, p1_opts@trXMaxF/), \
                             (/0.75, 0.75/), opts_l)
  
  line2   = gsn_add_polyline(wks, res_spd, (/p2_opts@trXMinF, p2_opts@trXMaxF/), \
                             (/0.75, 0.75/), opts_l)

;===============================
; Add and customise a legend
;===============================

  ; Display legend, control font and label height 
    lg_opts                            = True
    lg_opts@pmLegendDisplayMode        = "Always" 
    lg_opts@lgAutoManage               = False    
    lg_opts@lgLabelFont                = "Helvetica"
    lg_opts@lgLabelFontHeightF         = 0.025
  
  ; Labels and line colours ('SP' and 'tendency')
    lg_opts@lgLineColors               = (/p0_opts@xyLineColors, var1_opts@xyLineColor/)
    lg_opts@lgDashIndexes              = (/0, 0/)
    lsize                              = dimsizes(lg_opts@lgDashIndexes)

  ; Label position within box
    lg_opts@lgLabelPosition            = "Right"
    lg_opts@lgItemPlacement            = "ExplicitPlacement"
    lg_opts@lgItemPositions            = fspan(0.20, 0.80, lsize)

    lg_opts@lgLineThicknessF           = 4.0
    lg_opts@lgLineDashSegLenF          = 0.15

  ; Legend box size
    lg_opts@vpHeightF                  = 0.045
    lg_opts@vpWidthF                   = 0.25

    lg_opts@lgPerimColor               = "black"
    lg_opts@lgPerimThicknessF          = 3.0
    lg_opts@lgPerimFill                = "SolidFill"
    lg_opts@lgPerimFillColor           = "white"

  ; Customise legend labels depending on number and choice of input variables               
    if (std .eq. 1) then 
     str                               = "standard deviation"
     lab_arr                           = (/"~F33~z~F21~ "+str,"Sea level pressure tendency"/)
    else
     lab_arr                           = (/"Vorticity ratio (R)",\
                                           "    Min. sea level ~C~ pressure tendency"/)
    end if 
    lgd                                = gsn_create_legend(wks, lsize, lab_arr, lg_opts)

  ; Add legend to plot (top left or right hand corner)
    am_opts                            = True

    if (lgd0 .eq. "tr") then 
     am_opts@amJust                    = "TopRight"
     am_opts@amParallelPosF            = 0.5
     am_opts@amOrthogonalPosF          = -0.5
    elseif (lgd0 .eq. "tl") then 
     am_opts@amJust                    = "TopLeft"
     am_opts@amParallelPosF	       = -0.5
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

  ; Additional legend resources --> tangential wind tendency
    lg_opts2                           = lg_opts
    if (std .eq. 1) then 
     lab_arr2                          = (/"~F33~z~F21~ "+str,"Tangential wind ~C~    tendency"/)
    else
     lab_arr2                          = (/"Vorticity ratio (R)",\
                                           "Tangential wind ~C~     tendency"/)
    end if 

  ; Change line colours to match panel plot 
    lg_opts2@lgLineColors              = (/p0_opts@xyLineColors, var0_opts@xyLineColor/)

    lgd2                               = gsn_create_legend(wks, lsize, lab_arr2, lg_opts2)

  ; Additional legend resources --> tangential wind 
    lg_opts3                           = lg_opts
    if (std .eq. 1) then 
     lab_arr3                          = (/"~F33~z~F21~ "+str,"Tangential wind"/)
    else
     lab_arr3                          = (/"Vorticity ratio (R)","Tangential wind"/)
    end if 

  ; Change line colours to match panel plot 
    lg_opts3@lgLineColors              = (/p0_opts@xyLineColors, var2_opts@xyLineColor/)
    lgd3                               = gsn_create_legend(wks, lsize, lab_arr3, lg_opts3)

  ; Add legend to first panel (tangential wind tendency)
    annotate_1                         = gsn_add_annotation(res_ten, lgd2, am_opts)

  ; Also add legend to second panel (MSLP)
    annotate_2                         = gsn_add_annotation(res_slp, lgd, am_opts)

  ; Add legend to third panel (tangential wind)
    annotate_3                         = gsn_add_annotation(res_spd, lgd3, am_opts)

;=============================================
; Also add standard deviation label to plot 
;=============================================

  ; Resources 
    sd_res                   = True
    sd_str                   = sd_string
    sd_res@txFontHeightF     = 0.015
    sd_res@txFont            = "Helvetica"
    xpos_sd                  = 0.80
    ypos_sd                  = 0.96

;    gsn_text_ndc(wks, sd_str, xpos_sd, ypos_sd, sd_res)

  ; Also draw box around the label (FIX LATER)
    poly_x   = (/0.75, 0.83, 0.83, 0.75/)
    poly_y   = (/0.94, 0.94, 0.98, 0.98/)

    poly_res = True
    poly_res@gsFillColor = -1
    
;    gsn_polygon_ndc(wks, poly_x, poly_y, poly_res)

  ; Draw the plot and overlay all additional objects 
    draw(res_spd)
    draw(res_ten)
    draw(res_slp)
    frame(wks)

 end if            ; End IF statement on plotting or writing out 

; Print output file path to screen
  print(output)

end
