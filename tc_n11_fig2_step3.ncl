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

; ncl opt=\"x11\" dat=\"02T12\" ens0=\"em11\" mlev0=13 ts=24 tf=96 run=1 lay=1
; tend0=1 vmax=1.0 vmin=-1.0 var=\"slp\" lgd0=\"tr\" mlev1=25 
; cn0=\"slp\" two=1 wrt=1 mth=\"p2\" tc_n11_fig2_step3.ncl

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
; 'var'	  = plot either pressure tendency ("slp") or maximum windspeed ("vtan")
; 'lgd0'  = place legend in top right ("tr") or top left ("tl") of final plot
; 'lay'	  = read in layer-averaged data (1) or single-layer data (0)
; 'mlev1' = upper model level if 'lay=1' (25, 28, 34, etc)
; 'wrt'   = only write out MSLP and vtan tendency (1) or plot as normal (0)

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
  time	  = d->t
  size 	  = dimsizes(time)

; Create height string (e.g. 980m) for output file based upon chosen model level (e.g. 13) 
  hy0	  = sprintf("%0.0f", hybD(mlev0))

  if (lay .eq. 1) then 
   hy1    = sprintf("%0.0f", hybD(mlev1))
   print("Working on layer average: "+hy0+" m to "+hy1+" m")
  else
   print("Working on height: "+hy0+" m")
  end if 

; Create new arrays 
  plotR   = new( (/2,119/), float)
  plotL   = new( (/2,119/), float)

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

; VORTICITY STANDARD DEVIATION
; 120 values 
  vort_arr        = asciiread(diri+"cp_sd_sc"+sc+"_"+dat+"_"+ens0+"_vort_lay_980m_3380m.txt", \
                              (/size/),"float")

;; Ignore the first element of the array 
;  vort_arr        = vort_arr0(1:119)

; Estimate and remove the linear trend of the time series 
  vort_dtrend     = dtrend_msg_n(time, vort_arr, True, False, 0)
  plotR(0,:)      = vort_dtrend

; SYMMETRY RATIO 
; 119 values
  sym_arr         = asciiread(diri+"cp_sym_sc"+sc+"_"+dat+"_"+ens0+".txt", (/size/), "float")

; Read symmetry ratio into array, and subtract 1 --> 0 is the new origin 
; < 0 --> ringlike (symmetric) 
; > 0 --> monopole (asymmetric) 
  plotR(1,:)      = sym_arr(:) - 1

; TANGENTIAL WIND AND TENDENCY
; Read in maximum tangential wind and tendency from text files (121 values)
; Only one model height level (calculated from maximum mean vtan, which exists only on one level)
  vtan_arr0       = asciiread(diri+"ml_"+dat+"_"+ens0+"_vtan_max_"+mth+".txt",(/121/),"float")
  vtan_tend0      = asciiread(diri+"ml_"+dat+"_"+ens0+"_vtan_tend_"+mth+".txt",(/121/),"float")

; MEAN SEA LEVEL PRESSURE 
; Also read in MSLP data from text files (not calculated using same tracking method, but fine)
; 120 values (T+1 to T+120)
  mslp_arr0       = asciiread(diri+"cp_new_"+dat+"_"+ens0+"_slp.txt",(/120/),"float")

; Make all input arrays the same size (119 values)
  vtan_arr        = vtan_arr0(2::)
  vtan_tend       = vtan_tend0(2::)
  mslp_arr        = mslp_arr0(1::)

; Calculate y-axis limits using range of data in text file
  xvt             = max(vtan_arr(ts-1:tf-1))+5
  nvt             = min(vtan_arr(ts-1:tf-1))-5

; Calculate y-axis limits using range of data in text file 
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

; Read tangential wind (and tendency) into new array
  if (var .eq. "slp") then 
   plotL(0,:)  	   = mslp_arr(:)			; MSLP
   plotL(1,:)	   = mslp_tend(:)			; MSLP tendency
  elseif (var .eq. "vtan") then 
   plotL(0,:)      = vtan_arr(:)                        ; Tangential wind
   plotL(1,:)      = vtan_tend(:)                       ; Tangential wind tendency
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

; Options for plotting inner-core symmetry
  p0_opts		       	      = True
  p0_opts@gsnMaximize          	      = False
  p0_opts@vpXF                        = 0.15
  p0_opts@vpYF                        = 0.93
  p0_opts@vpWidthF             	      = 0.70
  p0_opts@vpHeightF		      = 0.35
  p0_opts@tmYROn	       	      = False	       ; Turn off tick marks
  p0_opts@tmXTOn	       	      = False	       ; Turn off tick marks
  p0_opts@trYMaxF	       	      = vmax	       ; Set y-axis maximum
  p0_opts@trYMinF	       	      = vmin  	       ; Set y-axis minimum

; Account for input arrays starting at T+2, not T+0
  p0_opts@trXMinF	       	      = ts             ; Set x-axis minimum
  p0_opts@trXMaxF	       	      = tf             ; Set x-axis maximum
 
  p0_opts@tmXBMode		      = "Explicit"
  p0_opts@tmXBValues	       	      = ispan(p0_opts@trXMinF, p0_opts@trXMaxF, 12)
  p0_opts@tmXBLabels                  = "T+"+tostring(p0_opts@tmXBValues)
  p0_opts@tiMainString         	      = "Simulation "+ens0+" ("+dat+")"

  if (var .eq. "vtan") then 

   p0_opts@tiYAxisString              = "Symmetry ratio"

  elseif (var .eq. "slp") then 

   p0_opts@tiYAxisString              = "Symmetry ratio"

  end if 

  p0_opts@tiXAxisString        	      = "Time"
  p0_opts@tiXAxisFontHeightF   	      = 0.0200    ; Reduce title font size
  p0_opts@tmXBLabelFontHeightF 	      = 0.0150

  p0_opts@tiYAxisFontHeightF   	      = 0.0200    ; y-axis font size
  p0_opts@tiMainFontHeightF    	      = 0.0225    ; title font size

  p0_opts@tmXMajorGridLineDashPattern = 2    ; Dashed gridlines
  p0_opts@tmYMajorGridLineDashPattern = 2
  p0_opts@gsnDraw		      = False
  p0_opts@gsnFrame		      = False
  p0_opts@gsnPaperOrientation	      = "Portrait"

  p0_opts@xyMarkLineMode	      = "MarkLines"
  p0_opts@xyLineColors		      = "red3"
  p0_opts@xyLineThicknesses	      = 3.0
  p0_opts@xyMarkers		      = 16
  p0_opts@xyMarkerColors      	      = "red3"
  p0_opts@xyMarkerSizes		      = 5.0
  p0_opts@xyMarkerThicknesses         = 3.0

; Resources for vorticity standard deviation 
  p1_opts                             = p0_opts 
  delete([/p1_opts@trYMinF,p1_opts@trYMaxF/])
  p1_opts@trYMaxF                     = 6.0
  p1_opts@trYMinF                     = -6.0
  p1_opts@tiYAxisString               = "10~S~-4~N~ s~S~-1~N~ / hPa h~S~-1~N~"
  p1_opts@vpXF                        = 0.15
  p1_opts@vpYF                        = 0.46
  p1_opts@vpWidthF                    = 0.70
  p1_opts@vpHeightF                   = 0.35
  p1_opts@tiMainString                = ""

; Options for plotting tangential wind or MSLP

  if (var .eq. "vtan") then 
   var1_opts			       = True
   var1_opts@tmXBLabelFontHeightF      = 0.01
   var1_opts@trYMaxF                   = xvt
   var1_opts@trYMinF                   = nvt
   var1_opts@tiMainString              = ""
   var1_opts@tiYAxisString             = "Windspeed tendency (m s~S~-1~N~ h~S~-1~N~)"
   var1_opts@tiYAxisFontHeightF        = 0.0125
   var1_opts@xyLineColors              = "black"
   var1_opts@xyLineThicknesses         = 4.0
   var1_opts@xyDashPatterns            = 0.0 
  elseif (var .eq. "slp")
   var1_opts                           = True
   var1_opts@tmXBLabelFontHeightF      = 0.01
   var1_opts@trYMaxF                   = 5.0
   var1_opts@trYMinF                   = -5.0
   var1_opts@tiMainString              = ""
   var1_opts@tiYAxisString             = "MSLP tendency (hPa h~S~-1~N~)"
   var1_opts@tiYAxisFontHeightF        = 0.0125
   var1_opts@xyLineColor               = "royalblue1"
   var1_opts@xyLineThicknessF          = 4.0
   var1_opts@xyDashPatterns            = 0.0
   var1_opts@xyMarkLineMode            = "MarkLines"
   var1_opts@xyMarker                  = 12
   var1_opts@xyMarkerColor             = "royalblue1"
   var1_opts@xyMarkerSizeF             = 0.01
   var1_opts@xyMarkerThicknessF        = 2.0
  end if 

; Options for plotting MSLP/tangential wind tendency
  tend_optsL   	       		       = True
  tend_optsL@gsLineColor               = "royalblue1"
  tend_optsL@gsLineThicknessF          = 2.0
 
  tend_optsM			       = True
  tend_optsM@gsMarkerIndex             = 7
  tend_optsM@gsMarkerColor             = "royalblue1"
  tend_optsM@gsMarkerSizeF             = 10.0
  tend_optsM@gsMarkerThicknessF        = 3.0

; Apply running average 
  if (run .eq. 1) then 
  
   plotR(0,:) = runave_n_Wrap(plotR(0,:), 6, 0, 0)
   plotR(1,:) = runave_n_Wrap(plotR(1,:), 6, 0, 0)
   plotL(1,:) = runave_n_Wrap(plotL(1,:), 6, 0, 0)

  end if 

; Either write out MSLP and tangential wind tendency to text files           
  if (wrt .eq. 1) then

   diri_out      = "$ar/text/"+dat+"_"+ens0+"_"+cn0+"_"

   slp_tend_out  = diri_out+"mslp_tend.txt"
   vtan_tend_out = diri_out+"vtan_tend.txt"
   vort_asym_out = diri_out+"vort_asym.txt"

 ; EDIT -- 'plotR' has 2 dimensions --> is this always the case?
 ; plotR(0,:)      = vort_dtrend
 ; plotR(1,:)      = sym_arr(:) - 1
   if (var .eq. "slp") then 
    asciiwrite(slp_tend_out, mslp_tend)
    asciiwrite(vort_asym_out, plotR(:) )
   elseif (var .eq. "vtan") then 
    asciiwrite(vtan_tend_out, vtan_tend)
    asciiwrite(vort_asym_out, plotR(:) )
   end if

; Or plot time series of asymmetry and pressure/tangential wind tendency
  else 

;===============================================
; Output file path and final plotting options
;===============================================

; Output file path and type 
  if (run .eq. 1) then 

   if (var .eq. "vtan") then 

    output = "$sam/nepartak/images/vc/n11_fig1a_sc"+sc+"_"+vort_str+"_lay_"+hy0+\
     	      "m_"+hy1+"m_"+dat+"_"+ens0+"_ave"

   elseif (var .eq. "slp")

    output = "$sam/nepartak/images/vc/n11_fig2_sc"+sc+"_"+vort_str+"_lay_"+hy0+\
             "m_"+hy1+"m_"+dat+"_"+ens0+"_ave"

   end if 

  elseif (run .eq. 0) then

   if (var .eq. "vtan") then 

    output = "$sam/nepartak/images/vc/n11_fig1a_sc"+sc+"_"+vort_str+"_lay_"+hy0+\
    	     "m_"+hy1+"m_"+dat+"_"+ens0

   elseif (var .eq. "slp") 

    output = "$sam/nepartak/images/vc/n11_fig2_sc"+sc+"_"+vort_str+"_lay_"+hy0+\
             "m_"+hy1+"m_"+dat+"_"+ens0

   end if 

  end if

; Open workstation for plotting 
  wks = gsn_open_wks(opt,output)

; SYMMETRY RATIO (NEW)  
; First, compare the new symmetry ratio [ plotR(1,:) ] and MSLP tendency [ plotL(1,:) ]
  if (two .eq. 1) then
   result = gsn_csm_xy2(wks,time,plotR(1,:),plotL(1,:),p0_opts,var1_opts) ; Only plot tendency
  elseif (two .eq. 0) then
   result = gsn_csm_xy2(wks,time,plotR,plotL(0:1,:),p0_opts,var1_opts) ; Plot both MSLP & tendency
  end if 

; VORTICITY STANDARD DEVIATION (EXISTING)
  res_sd  = gsn_csm_xy2(wks,time,plotR(0,:),plotL(1,:),p1_opts,var1_opts) ; 

; Create graphical object to overlay horizontal line onto plot 
  line    = new(1,graphic)

; Line resources 
  opts_l  = True
  opts_l@gsLineThicknessF  = 2.0       ; Line thickness
  opts_l@gsLineDashPattern = 1         ; Dash pattern 
  opts_l@gsLineColor       = "black"   ; Line colour 

; Overlay horizontal line through y=0
  line    = gsn_add_polyline(wks, result, (/p0_opts@trXMinF, p0_opts@trXMaxF/), \
                             (/0.0, 0.0/), opts_l)

  line1   = gsn_add_polyline(wks, res_sd, (/p1_opts@trXMinF, p1_opts@trXMaxF/), \
                             (/0.0, 0.0/), opts_l)

;===============================
; Add and customise a legend
;===============================

    lg_opts                            = True
    lg_opts@pmLegendDisplayMode        = "Always"          ; Display a legend
    lg_opts@lgAutoManage               = False             ; Gives us control
    lg_opts@lgLabelFont                = "Helvetica"       ; Label font
    lg_opts@lgLabelFontHeightF         = 0.030             ; Height of labels

  ; Customise labels and line colours depending on number of labels
    if (two .eq. 1) then
     lg_opts@lgLineColors               = (/p0_opts@xyLineColors, \      ; Asymmetry       
                                            tend_optsL@gsLineColor/)     ; Tendency        
     lg_opts@lgDashIndexes              = (/0, 0/)
    elseif (two .eq. 0) then
     lg_opts@lgLineColors               = (/p0_opts@xyLineColors, \      ; Asymmetry       
                                            var1_opts@xyLineColors, \    ; MSLP/windspeed  
                                            tend_optsL@gsLineColor/)     ; Tendency        
     lg_opts@lgDashIndexes              = (/0, 0, 0/)
    end if
    
    lsize                              = dimsizes(lg_opts@lgDashIndexes)

    lg_opts@lgLabelPosition            = "Right"           ; Label position within box
    lg_opts@lgItemPlacement            = "ExplicitPlacement"
    lg_opts@lgItemPositions            = fspan(0.20, 0.80, lsize)

    lg_opts@lgLineThicknessF           = 4.0
    lg_opts@lgLineDashSegLenF          = 0.20
    lg_opts@vpWidthF                   = 0.275

  ; Customise legend box size depending on number of labels                                 
    if (two .eq. 1) then
     lg_opts@vpHeightF                 = 0.0625
    elseif (two .eq. 0) then
     lg_opts@vpHeightF                 = 0.125
    end if

    lg_opts@lgPerimColor               = "black"
    lg_opts@lgPerimThicknessF          = 3.0
    lg_opts@lgPerimFill                = "SolidFill"
    lg_opts@lgPerimFillColor           = "white"

  ; Customise legend labels depending on number and choice of input variables               
     if (var .eq. "vtan") then

      if (two .eq. 1) then
      lab_arr                          = (/"Relative vorticity standard deviation",\
                                           "Total wind tendency"/)
      elseif (two .eq. 0) then
      lab_arr                          = (/"Relative vorticity standard deviation",\
                                           "Max. total wind", "Total wind tendency"/)
      end if

     elseif (var .eq. "slp") then

      if (two .eq. 1) then
      lab_arr                          = (/"Inner-core symmetry parameter",\
                                           "MSLP tendency"/)
      elseif (two .eq. 0) then
      lab_arr                          = (/"Inner-core symmetry parameter",\
                                           "Minimum MSLP", "MSLP tendency"/)
      end if

     end if

    lgd                                = gsn_create_legend(wks,lsize,lab_arr,lg_opts)

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

  ; Additional legend resources 
    am_opts2                           = am_opts
    lg_opts2                           = lg_opts

    lab_arr2                           = (/"Relative vorticity standard deviation",\
                                           "MSLP tendency"/)

    lgd2                               = gsn_create_legend(wks,lsize,lab_arr2,lg_opts2)

  ; Add legend to first panel
    annotate_1                         = gsn_add_annotation(result, lgd, am_opts)

  ; Also add legend to second panel 
    annotate_2                         = gsn_add_annotation(res_sd, lgd2, am_opts2)

; Draw the plot and overlay all additional objects 
  draw(result)
  draw(res_sd)
  frame(wks)

 end if            ; End IF statement on plotting or writing out 

end
