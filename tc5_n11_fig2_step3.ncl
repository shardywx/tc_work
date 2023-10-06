; Plot vorticity asymmetry and pressure tendency (Fig. 2 from N11)
; using input from 'tc5_n11_fig2_step1.ncl' and 'tc5_n11_fig2_step2.ncl'

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/times_arr.ncl"

; ncl opt=\"x11\" dat=\"02T12\" ens0=\"em11\" mlev0=13 ts=491 tf=983 run=1 vort0=\"rv\" lay=1
; trend=1 tend0=1 vmax=10.0 vmin=-10.0 var=\"slp\" lgd0=\"tr\" mean0=True mlev1=25 
; cn0=\"slp\" tc5_n11_fig2_step3.ncl

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'ens0'  = ensemble simulation (em00-em11)
; 'opt'   = output file format ("pdf" or "x11")
; 'run'	  = apply running average (1) or leave alone (0)
; 'ts'	  = start time for plot
; 'tf'	  = end time for plot 
; 'mlev0' = chosen model level for analysis (13, 16, 22, etc)
; 'vort0' = choose relative ("rv") or potential ("pv") vorticity
; 'trend' = use NCL function to remove linear trend (1) or leave alone (0)
; 'tend0' = calculate windspeed/MSLP tendency using t2-t0 (2) or t1-t0 (1)
; 'mean0' = if trend=1, set mean=True to subtract the mean as well as the linear trend 
; 'vmax'  = maximum standard deviation for final plot 
; 'vmin'  = minimum standard deviation for final plot
; 'var'	  = plot either pressure tendency ("slp") or maximum windspeed ("vtan")
; 'lgd0'  = place legend in top right ("tr") or top left ("tl") of final plot
; 'lay'	  = read in layer-averaged data (1) or single-layer data (0)
; 'mlev1' = upper model level if 'lay=1' (25, 28, 34, etc)
; 'cn0'   = method of storm tracking ["slp" or "slpf"]

begin

;======================================
; Read in the data from netCDF files
;======================================

; Get model level information from 'pd' stream, for output file string
  diri_f  = "/nfs/a319/earshar/"+dat+"/"+ens0+"/vc"
  fili_pr = "201607"+dat+"00Z_NPTK_4p4_L80_ra1t_"+ens0
  fili_d1 = systemfunc("cd "+diri_f+" ; ls "+fili_pr+"_pd.nc")
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
  print("Working on height: "+hy0+" m")

  if (lay .eq. 1) then 
   hy1    = sprintf("%0.0f", hybD(mlev1))
   print("Working on layer average: "+hy0+" m to "+hy1+" m")
  end if 

;====================================
; Read in the data from	text files
;====================================

; Directory path
  diri = "$ar/text/"

; Define part of output file string (relative vorticity vs PV)
  if (vort0 .eq. "rv") then
   vort_str = "vort"
  elseif (vort0 .eq. "pv") then
   vort_str = "pv"
  end if

; Read in the standard deviation of either PV or relative vorticity 
  if (vort0 .eq. "pv") then 

 ; Potential vorticity
   if (lay .eq. 0) then    
    vort_arr      = asciiread(diri+"ml5_"+dat+"_"+ens0+"_"+cn0+"_pv_"+hy0+"m.txt",(/size/),"float")
   elseif (lay .eq. 1) then 
    vort_arr      = asciiread(diri+"ml5_"+dat+"_"+ens0+"_"+cn0+"_pv_lay_"+hy0+"m_"+hy1+"m.txt",\
    		    (/size/),"float")
   end if 

; Relative vorticity
  elseif (vort0 .eq. "rv") then 

   if (lay .eq. 0) then 
    vort_arr      = asciiread(diri+"ml5_"+dat+"_"+ens0+"_"+cn0+"_vort_"+hy0+"m.txt",\
                             (/size/),"float")
   elseif (lay .eq. 1) then 
    vort_arr      = asciiread(diri+"ml5_"+dat+"_"+ens0+"_"+cn0+"_vort_lay_"+hy0+"m_"+hy1+"m.txt",\
    		             (/size/),"float")
   end if 

  end if 

; Estimate and remove the (linear) trend of the time series 
  if (trend .eq. 1) then 
   vort_dtrend    = dtrend_msg_n(time, vort_arr, mean0, False, 0)
  end if 

; Read in maximum tangential wind and tendency from text files
; Need to calculate 'vtan_arr' using storm track methods other than 'geo_sm'
  vtan_arr        = asciiread(diri+"ml5_"+dat+"_"+ens0+"_vtan_max_"+cn0+".txt",\
                             (/size/),"float")
;  tend_arr        = asciiread(diri+"ml5_"+dat+"_"+ens0+"_"+cn0+"_"+hy0+"m_tend.txt",\
;                             (/size/),"float")

; Calculate y-axis limits using range of data in text file
  xvt              = max( vtan_arr(ts-1:tf-1) ) + 5
  nvt              = min( vtan_arr(ts-1:tf-1) ) - 5

; Read relative vorticity/PV standard deviation into new array ('plotR')
  if (trend .eq. 1) then 
   plotR           = vort_dtrend(:)
  else
   plotR           = vort_arr(:)
  end if

; Also read in MSLP data from text files 
; Investigate how changing the method of calculating the MSLP minimum changes our results
  mslp_arr         = asciiread(diri+"cp5_"+dat+"_"+ens0+"_slp.txt",(/size/),"float")

; Calculate y-axis limits using range of data in text file 
  xslp	    	   = max(mslp_arr(ts-1:tf-1))+10
  nslp		   = min(mslp_arr(ts-1:tf-1))-10

; Calculate MSLP and windspeed tendencies from 'mslp_arr' and 'vtan_arr'
  mslp_tend        = mslp_arr
  vtan_tend        = vtan_arr
  do it = 0, size-1

   tp1 = min((/it+1,size-1/))	
   tm1 = max((/it-1,0/))

   if (tend0 .eq. 2) then 

  ; Divide by 2, multiply by 12 
    mslp_tend(it) = (mslp_arr(tp1) - mslp_arr(tm1)) * 6
    vtan_tend(it) = (vtan_arr(tp1) - vtan_arr(tm1)) * 6

   elseif (tend0 .eq. 1) then 

  ; Divide by 1, multiply by 12 
    mslp_tend(it) = (mslp_arr(it) - mslp_arr(tm1)) * 12
    vtan_tend(it) = (vtan_arr(it) - vtan_arr(tm1)) * 12 
   end if 

  end do 

; Read tangential wind (and tendency) into new array
  plotL		   = new((/2,size/),"float")
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

; Options for plotting vorticity standard deviation 
  pv_opts		       	      = True
  pv_opts@gsnMaximize          	      = True
  pv_opts@vpWidthF             	      = 0.8
  pv_opts@vpHeightF		      = 0.5
  pv_opts@tmYROn	       	      = False	       ; Turn off tick marks
  pv_opts@tmXTOn	       	      = False	       ; Turn off tick marks
  pv_opts@trYMaxF	       	      = vmax	       ; Set maximum standard deviation
  pv_opts@trYMinF	       	      = vmin  	       ; Set minimum standard deviation
  pv_opts@trXMinF	       	      = ts             ; Set x-axis minimum
  pv_opts@trXMaxF	       	      = tf             ; Set x-axis maximum
 
  pv_opts@tmXBMode		      = "Explicit"
  pv_opts@tmXBValues	       	      = ispan(ts, tf, 144)
  pv_opts@tmXBLabels           	      = (/"T+42","T+54","T+66","T+78"/)

  pv_opts@tiMainString         	      = "Simulation "+ens0+" ("+dat+")"

  if (var .eq. "vtan") then 

   if (vort0 .eq. "pv") then 
    pv_opts@tiYAxisString              = "PVU / m s~S~-1~N~ h~S~-1~N~"
   elseif (vort0 .eq. "rv") then 
    pv_opts@tiYAxisString              = "10~S~-4~N~ s~S~-1~N~ / m s~S~-1~N~ h~S~-1~N~"
   end if 

  elseif (var .eq. "slp") then 

   if (vort0 .eq. "pv") then
    pv_opts@tiYAxisString              = "PVU / hPa h~S~-1~N~"
   elseif (vort0 .eq. "rv") then
    pv_opts@tiYAxisString              = "10~S~-4~N~ s~S~-1~N~ / hPa h~S~-1~N~"
   end if

  end if 

  pv_opts@tiXAxisString        	      = "Time"
  pv_opts@tiXAxisFontHeightF   	      = 0.0200    ; Reduce title font size
  pv_opts@tmXBLabelFontHeightF 	      = 0.0150

  pv_opts@tiYAxisFontHeightF   	      = 0.0200    ; y-axis font size
  pv_opts@tiMainFontHeightF    	      = 0.0225    ; title font size

  pv_opts@tmXMajorGridLineDashPattern = 2    ; Dashed gridlines
  pv_opts@tmYMajorGridLineDashPattern = 2
  pv_opts@gsnDraw		      = False
  pv_opts@gsnFrame		      = False
  pv_opts@gsnPaperOrientation	      = "Portrait"

  pv_opts@xyMarkLineMode	      = "MarkLines"
  pv_opts@xyLineColors		      = "red3"
  pv_opts@xyLineThicknesses	      = 3.0
  pv_opts@xyMarkers		      = 16
  pv_opts@xyMarkerColors      	      = "red3"
  pv_opts@xyMarkerSizes		      = 5.0
  pv_opts@xyMarkerThicknesses         = 3.0

; Options for plotting tangential wind

  if (var .eq. "vtan") then 
   var1_opts			       = True
   var1_opts@tmXBLabelFontHeightF      = 0.01
   var1_opts@trYMaxF                   = xvt
   var1_opts@trYMinF                   = nvt
   var1_opts@tiMainString              = ""
   var1_opts@tiYAxisString             = "Maximum windspeed (m s~S~-1~N~)"
   var1_opts@tiYAxisFontHeightF        = 0.0125
   var1_opts@xyLineColors              = "black"
   var1_opts@xyLineThicknesses         = 4.0
   var1_opts@xyDashPatterns            = 0.0 
  elseif (var .eq. "slp")
   var1_opts                           = True
   var1_opts@tmXBLabelFontHeightF      = 0.01
   var1_opts@trYMaxF                   = xslp
   var1_opts@trYMinF                   = nslp
   var1_opts@tiMainString              = ""
   var1_opts@tiYAxisString             = "Minimum MSLP (hPa)"
   var1_opts@tiYAxisFontHeightF        = 0.0125
   var1_opts@xyLineColors              = "black"
   var1_opts@xyLineThicknesses         = 4.0
   var1_opts@xyDashPatterns            = 0.0
  end if 

; Options for plotting MSLP/tangential wind tendency
  tend_optsL   	       		       = True
  tend_optsL@gsLineColor               = "royalblue1"
  tend_optsL@gsLineThicknessF          = 2.0
 
  tend_optsM			       = True
  tend_optsM@gsMarkerIndex             = 7
  tend_optsM@gsMarkerColor             = "royalblue1"
  tend_optsM@gsMarkerSizeF             = 5.0
  tend_optsM@gsMarkerThicknessF        = 3.0

; Apply running average 
  if (run .eq. 1) then 
  
   plotR(:)   = runave_n_Wrap(plotR(:), 24, 0, 0)

   if (var .eq. "vtan") then 
    plotL(1,:) = runave_n_Wrap(plotL(1,:), 36, 0, 0)
   elseif (var .eq. "slp") then 
    plotL(1,:) = runave_n_Wrap(plotL(1,:), 24, 0, 0)
   end if 

  end if 

;=============================================
; Print out how asymmetry changes with time
;=============================================

 ; Call external function 'times_arr'
 ; Returns 'time_str' (12 UTC 4 July), 'time_arr' (), and 'title_arr' ()
   all_times  = times_arr(dat, ens0)

 ; Return strings e.g. 0600 UTC 04 July 
   time_str   = all_times[0]

 ; Return strings e.g. 04Jul_0600Z
   time_arr   = all_times[1]

   do it = ts, tf

    print("Vort-SD = "+plotR(it)+" at "+time_str(it)+" (index = "+it+")" )

   end do 

;===================================================
; Find local maxima/minima in vort-SD and overlay 
;===================================================

 ; 'False'    = Is data cyclic (True) or not (False)?
 ; '0.15'     = Tolerance level when finding local maxima (higher value = stricter condition)
 ; '1'        = Return indices of maxima, rather than the values themselves (0)

 ; Missing values cause the function to crash, so ignore these parts of the 'plotR' array
   sd_max0    = local_max_1d(plotR(ts-1:tf-1), False, 0.50, 1)
   sd_min0    = local_min_1d(plotR(ts-1:tf-1), False, -0.50, 1)

 ; Now calculate the correct indices after taking into account the starting time
   sd_max     = sd_max0 + ts
   sd_min     = sd_min0 + ts

   print("Indices of local maxima in asymmetry: " +sd_max)
   print("Indices of local minima in asymmetry: " +sd_min)

;===============================================
; Also find local minima in pressure tendency
;===============================================

 ; 'False'    = Is data cyclic (True) or not (False)?
 ; '0.15'     = Tolerance level when finding local maxima (higher value = stricter condition)
 ; '1'        = Return indices of maxima, rather than the values themselves (0)

 if (var .eq. "slp") then 

 ; Calculate the indices corresponding to the pressure tendency minima
 ; Missing values cause	the function to	crash, so ignore these parts of	the 'plotL' array
   pt_min0    = local_min_1d(plotL(1,ts-1:tf-1), False, -0.08, 1)

 ; Also calculate the absolute values of the minima
   pt_min1    = local_min_1d(plotL(1,ts-1:tf-1), False, -0.08, 0)

 ; Now calculate the correct indices after taking into account the starting time
   pt_min     = pt_min0 + ts 

 elseif (var .eq. "vtan")

 ; Calculate the indices corresponding to the tangential wind tendency minima
 ; Missing values cause the function to crash, so ignore these parts of the 'plotL' array
   pt_min0    = local_min_1d(plotL(1,ts-1:tf-1), False, -0.10, 1)
   pt_max0    = local_max_1d(plotL(1,ts-1:tf-1), False, 0.10, 1)

 ; Now calculate the correct indices after taking into account the starting time
   pt_min     = pt_min0 + ts
   pt_max     = pt_max0 + ts

 end if 

;===============================================
; Output file path and final plotting options
;===============================================

; Output file path and type 
  if (run .eq. 1) then 

   if (var .eq. "vtan") then 

    if (lay .eq. 1) then 
     output = "$sam/nepartak/images/vc/n11_fig1a_5min_"+vort_str+"_lay_"+hy0+\
     	      "m_"+hy1+"m_"+dat+"_"+ens0+"_ave"
    elseif (lay .eq. 0) then 
     output = "$sam/nepartak/images/vc/n11_fig1a_5min_"+vort_str+"_"+hy0+"m_"+dat+"_"+ens0+"_ave"
    end if 

   elseif (var .eq. "slp")

    if (lay .eq. 1) then 
     output = "$sam/nepartak/images/vc/n11_fig2_5min_"+vort_str+"_lay_"+hy0+\
              "m_"+hy1+"m_"+dat+"_"+ens0+"_ave"
    elseif (lay .eq. 0) then 
     output = "$sam/nepartak/images/vc/n11_fig2_5min_"+vort_str+"_"+hy0+"m_"+dat+"_"+ens0+"_ave"
    end if 

   end if 

  elseif (run .eq. 0) then

   if (var .eq. "vtan") then 

    if (lay .eq. 1) then 
     output = "$sam/nepartak/images/vc/n11_fig1a_5min_"+vort_str+"_lay_"+hy0+\
     	      "m_"+hy1+"m_"+dat+"_"+ens0
    elseif (lay .eq. 0) then 
     output = "$sam/nepartak/images/vc/n11_fig1a_5min_"+vort_str+"_"+hy0+"m_"+dat+"_"+ens0
    end if 

   elseif (var .eq. "slp") 

    if (lay .eq. 1) then 
     output = "$sam/nepartak/images/vc/n11_fig2_5min_"+vort_str+"_lay_"+hy0+\
              "m_"+hy1+"m_"+dat+"_"+ens0
    elseif (lay .eq. 0) then 
     output = "$sam/nepartak/images/vc/n11_fig2_5min_"+vort_str+"_"+hy0+"m_"+dat+"_"+ens0
    end if 

   end if 

  end if

  wks = gsn_open_wks(opt,output)

; First, plot 2 graphs
  if (two .eq. 1) then 
   result = gsn_csm_xy2(wks,time,plotR,plotL(1,:),pv_opts,var1_opts)    ; Only plot tendency
  elseif (two .eq. 0) then
   result = gsn_csm_xy2(wks,time,plotR,plotL(0:1,:),pv_opts,var1_opts)  ; Plot both MSLP & tendency
  end if 

; Now overlay 3rd graph (line & markers separately)
  tendL  = gsn_add_polyline(wks,result,time(:),plotL(1,:),tend_optsL)
  tendM  = gsn_add_polymarker(wks,result,time(:),plotL(1,:),tend_optsM)

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
     lg_opts@lgLineColors               = (/pv_opts@xyLineColors, \      ; Asymmetry
                                            tend_optsL@gsLineColor/)     ; Tendency
     lg_opts@lgDashIndexes              = (/0, 0/)
    elseif (two .eq. 0) then 
     lg_opts@lgLineColors               = (/pv_opts@xyLineColors, \      ; Asymmetry
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
    if (vort0 .eq. "rv") then 

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
      lab_arr                          = (/"Relative vorticity standard deviation",\
                                           "MSLP tendency"/)
      elseif (two .eq. 0) then 
      lab_arr                          = (/"Relative vorticity standard deviation",\
                                           "Minimum MSLP", "MSLP tendency"/)
      end if 
     end if 					   

    elseif (vort0 .eq. "pv") then 

     if (var .eq. "vtan") then 
      lab_arr                          = (/"PV standard deviation", "Max. total wind", \
                                            "Total wind tendency"/)
     elseif (var .eq. "slp") then 
      lab_arr                          = (/"PV standard deviation", "Minimum MSLP", \
                                            "MSLP tendency"/)
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

    annotate_1                         = gsn_add_annotation(result,lgd,am_opts)

;========================================================================
; Also overlay lines denoting the asymmetric (and symmetric) VC phases
;========================================================================

;====================
; Asymmetric phase
;====================

  a_size		 = dimsizes(sd_max)
  vc_asym  		 = new(a_size, graphic)

; Plotting resources
  ares                   = True
  ares@gsLineColor       = "orange2"
  ares@gsLineThicknessF  = 5.0
  ares@gsLineDashPattern = 16
  ares@tfPolyDrawOrder   = "Draw"

; Define the asymmetric phases from the standard deviation and the pressure/wind tendency

;  if (var .eq. "slp") then    	    		  ; Plot if dp/dt is more negative than threshold

;   do i = 0, a_size-1
;    if ( plotL(1,sd_max(i) ) .lt. -1.5) then
;     print("x-axis index (asym) = "+sd_max(i))
;     vc_asym(i) = gsn_add_polyline(wks, result, (/sd_max(i), sd_max(i)/), (/vmin, vmax/), \
;                                   ares)
;    end if
;   end do

;  elseif (var .eq. "vtan") then		  ; Plot if dv/dt is closer to zero than threshold

;   print("We use pressure, not windspeed tendency to define VCs!")

;  end if 

;===================
; Symmetric phase
;===================

  s_size		 = dimsizes(sd_min)
  vc_sym   		 = new(s_size, graphic)

; Plotting resources
  sres                   = ares
  sres@gsLineDashPattern = 0

;  if (var .eq. "slp") then                        ; Plot if dp/dt is smaller than threshold

;   do i = 0, s_size-1
;    if ( abs( plotL(1,sd_min(i) ) ) .lt. 1.5) then 
;     print("x-axis index (sym) = "+sd_min(i))
;     vc_sym(i) = gsn_add_polyline(wks, result, (/sd_min(i), sd_min(i)/), (/vmin, vmax/), \
;     	      	 		  sres)
;    end if 
;   end do

;  elseif (var .eq. "vtan") then                   ; Plot if dv/dt is larger than threshold

;   print("We use pressure, not windspeed tendency to define VCs!")

;  end if 

; Draw the plot and overlay all additional objects 

  draw(result)
  frame(wks)

end
