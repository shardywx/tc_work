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

; ncl opt=\"x11\" lgd0=\"br\" z0=1 z1=6 sc=5 r0=0 r1=50 ts0=563 tf0=995 
; diff=1 p2=\"vr\" trun=24 vrun=12 frun=12 vc_ts_mflux.ncl

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'ens0'  = ensemble simulation (em00-em11)
; 'opt'   = output file format ("pdf" or "x11")
; 'z0'    = height of lower level in km
; 'z1'    = height of upper level in km 
; 'lgd0'  = place legend in top right ("tr") or top left ("tl") of final plot
; 'diff'  = calculate difference between layers (1) or plot both layers 
; 'r0'    = inner radius in km (0,5,10)
; 'r1'    = outer radius in km (30,40,50)
; 'p2'    = variable to overlay onto vertical mass flux plot ("vr", "ten")
; 'trun'  = length of time for running average on tangential wind tendency 
; 'vrun'  = length of time for running average on relative vorticity ratio
; 'frun'  = length of time for running average on vertical mass flux (from 'mflux5_write')

begin

;======================================
; Read in the data from netCDF files
;======================================

; Get some useful information from the netCDF file (5-min data)
  diri_f  = "/nfs/a319/earshar/02T12/em11/vc"
  fili_pr = "20160702T1200Z_NPTK_4p4_L80_ra1t_em11"

  fili_d1 = systemfunc("cd "+diri_f+" ; ls "+fili_pr+"_pd*.nc")
  fili_d  = diri_f+"/"+fili_d1

  if ( .not. fileexists (fili_d) ) then
   print("Input	file missing! Exiting script...")
   exit()
  end if

; Calculate number of times in each file 
  numINFO_d   = nc_times(fili_d)
  numINPUT_d  = numINFO_d[0]
  time_d      = numINFO_d[1]

;=======================================================================  
; Need to read in time information from 5-min simulation (for x-axis) 
;=======================================================================  

  month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
                    "Jul","Aug","Sep","Oct","Nov","Dec"/)

  utc_date = cd_calendar(time_d, 0)

  year   = tointeger(utc_date(:,0))
  month  = tointeger(utc_date(:,1))
  day    = tointeger(utc_date(:,2))
  hour   = tointeger(utc_date(:,3))
  minute = tointeger(utc_date(:,4))
  second = utc_date(:,5)

; Create arrays to hold finished date strings                                        
  hr_min = new(numINPUT_d,string)
  min0   = (/4,9,14,19,24,29,34,39,44,49,54,59/)      ; Minute array (before)        
  min1   = (/5,10,15,20,25,30,35,40,45,50,55,0/)      ; Minute array (after)         

; Before correcting, find the indices of the elements where 'minute' = 59            
  hr_ind = ind(minute .eq. 59)

; Correct for errors in the code (round up values of 'minute' to multiples of 5)     
  do i = 0, numINPUT_d-1
   do m = 0, dimsizes(min0)-1

    if (minute(i) .eq. min0(m) ) then
     minute(i) = min1(m)
    end if

   end do
  end do

; Edit all incorrect values of 'hour' (where 'minute' previously = 59)               
  do j = 0, dimsizes(hr_ind)-1
   hour(hr_ind(j)) = hour(hr_ind(j)) + 1
  end do
  delete(j)

; Create string array of all values in 'minute'                                      
  mins = tostring(minute)

; Create 'hhmm' strings for each time (for output and title strings)                 
  do i = 0, numINPUT_d-1

 ; Add zeros in front of single-digit hours                                          
   if (hour(i) .lt. 10) then

    if (minute(i) .lt. 10) then
     hr_min(i) = "0"+hour(i)+"0"+mins(i)
    elseif (minute(i) .ge. 10) then
     hr_min(i) = "0"+hour(i)+mins(i)
    end if

   else

  ; Also add zeros in front of single-digit minutes                                  
    if (minute(i) .lt. 10) then
     hr_min(i) = hour(i)+"0"+mins(i)
    elseif (minute(i) .ge. 10) then
     hr_min(i) = hour(i)+mins(i)
    end if

   end if

  end do
  delete(i)

; Create arrays to hold time information for ALL times 
  time_str0  = new(numINPUT_d,string)
  time_arr0  = new(numINPUT_d,string)
  title_arr0 = new(numINPUT_d,string)

; Construct these strings for x-axis and output file paths
  time_str0  = hr_min + " UTC " + sprinti("%0.2i ", day) + month_abbr(month)
  time_arr0  = hr_min+ "Z " + sprinti("%0.2i", day) + " " + month_abbr(month)
  title_arr0 = "Valid at "+time_str0

; Now get the subset of dates we're interested in 
; 18 UTC 3 July (it = 347) --> 00 UTC 6 July (it = 995)
  time_str   = time_str0(ts0:tf0)
  time_arr   = time_arr0(ts0:tf0)
  title_arr  = title_arr0(ts0:tf0)

;=====================================
; Now grad model height information
;=====================================

; Retrieve model level array from 'pc' stream data
  d	  = addfile(fili_d,"r")	
  hybD	  = toint(d->hybrid_ht(:) )

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

;==========================
; Select your two levels
;==========================

; First level
  if (z0 .eq. 1.0) then
   zstr0   = "1125"
  elseif (z0 .eq. 1.5) then
   zstr0   = "1620"
  elseif (z0 .eq. 2.0) then
   zstr0   = "2205"
  elseif (z0 .eq. 3.0) then
   zstr0   = "3125"
  elseif (z0 .eq. 4.0) then
   zstr0   = "4205"
  elseif (z0 .eq. 6.0) then
   zstr0   = "6128"
  elseif (z0 .eq. 7.0) then
   zstr0   = "7226"
  elseif (z0 .eq. 8.0) then
   zstr0   = "8011"
  elseif (z0 .eq. 9.0) then
   zstr0   = "9268"
  elseif (z0 .eq. 10.0) then
   zstr0   = "10624"
  end if

; Second level
  if (z1 .eq. 1.0) then
   zstr1   = "1125"
  elseif (z1 .eq. 1.5) then
   zstr1   = "1620"
  elseif (z1 .eq. 2.0) then
   zstr1   = "2205"
  elseif (z1 .eq. 3.0) then
   zstr1   = "3125"
  elseif (z1 .eq. 4.0) then
   zstr1   = "4205"
  elseif (z1 .eq. 6.0) then
   zstr1   = "6128"
  elseif (z1 .eq. 7.0) then
   zstr1   = "7226"
  elseif (z1 .eq. 8.0) then
   zstr1   = "8011"
  elseif (z1 .eq. 9.0) then
   zstr1   = "9268"
  elseif (z1 .eq. 10.0) then
   zstr1   = "10624"
  end if

; Print out information to screen about model level                                          
  print("First layer centred on: "+zstr0+" m")
  print("Second layer centred on: "+zstr1+" m")

;====================================                                             
; Read in the data from text files                                                
;====================================                                             

; Directory path                                                                  
  diri   = "$ar/text/flux5_all_sc5_lay1_"

; Create radius strings to differentiate between text files 
  ri     = r0+"km"
  ro     = r1+"km"

; Select the two files we want 
  fili_0 = systemfunc("cd $ar/text ; ls "+diri+zstr0+"_w0.3_"+ri+"_"+ro+"_*run"+frun+".txt")
  fili_1 = systemfunc("cd $ar/text ; ls "+diri+zstr1+"_w0.3_"+ri+"_"+ro+"_*run"+frun+".txt")

; Merge into a single variable 
  fili   = (/fili_0, fili_1/)

; Number of files                                                                 
  nfil   = dimsizes(fili)

; Number of times in each file                                                    
  ntime0 = systemfunc("wc -l < "+fili(0) )
  ntimes = toint(ntime0)

; Create array to hold all data                                                   
  mflux_arr = new( (/nfil, ntimes/), "float")

; Read in data from all text files                                                
  do nf = 0, nfil-1
   mflux_arr(nf,:) = asciiread(fili(nf),ntimes,"float")
  end do

; Add useful metadata                                                             
  mflux_arr!0 = "level"
  mflux_arr&level = (/zstr0+"m",zstr1+"m"/)
  mflux_arr!1 = "time"

;==========================================================================
; Also read in maximum tangential wind tendency and vorticity ratio data
;==========================================================================

; Maximum tangential wind
  ntot = (tf0 - ts0) + 1
  vtan_arr    = asciiread("$ar/text/vtan5_max_all_sc5_w0.3_max.txt",ntot,"float")
  vtan_arr!0  = "time"

; Calculate tangential wind tendency (5-min data --> divide by 2, multiply by 12)
  vten_plt    = vtan_arr
  do it = 1, ntot-2
   vten_plt(it) = (vtan_arr(it+1) - vtan_arr(it-1) ) * 6  
  end do 

; Add metadata 
  vten_plt@units = "m s~S~-1~N~ h~S~-1~N~"

; Ignore first and last values 
  vten_plt(0)      = vten_plt@_FillValue
  vten_plt(ntot-1) = vten_plt@_FillValue

; Calculate 1-h running average to smooth some of the higher-frequency fluctuations 
  vten_arr         = runave_n_Wrap(vten_plt, trun, 0, 0)

; Relative vorticity ratio 
  vr_arr   = asciiread("$ar/text/cp5_vr_sc5.txt",ntot,"float")
  vr_arr!0 = "time"

; Calculate 1-h running average on the vorticity ratio too 
  vr_arr   = runave_n_Wrap(vr_arr, vrun, 0, 0)

;================================================================
; Calculate the difference between mass flux within two layers
;================================================================

; Calculate the difference between layers 
  mflux_diff  = mflux_arr(1,:) - mflux_arr(0,:)

; Express difference as percentage of the mass flux in the lower layer 
  mflux_diff  = (mflux_diff / mflux_arr(0,:) ) * 100
  mflux_diff@units = "%"

;======================================
; Create an x-y time series plot
;======================================

; Sort out time array for x-axis labels
  time        = ispan(0,ntimes-1,1)
  times       = dimsizes(time)

; Time series plotting options
  p0_opts		       	       = True
  p0_opts@gsnMaximize          	       = True

; Explicitly set plot dimensions (comment out if not required)
  p0_opts@vpWidthF                     = 0.55
  p0_opts@vpHeightF                    = 0.30

; Turn off tick marks on RHS y-axis and top x-axis
  p0_opts@tmYROn	       	       = False
  p0_opts@tmXTOn	       	       = False

; Set y-axis limits (EDIT)
  if (diff .eq. 1) then 
   p0_opts@trYMaxF	       	       = 100.0
   p0_opts@trYMinF	       	       = -100.0
  end if 

; Set x-axis limits 
  p0_opts@trXMinF	       	       = time(216)
  p0_opts@trXMaxF	       	       = time(576)
 
  p0_opts@tmXBMode		       = "Explicit"
  p0_opts@tmXBValues	       	       = ispan(p0_opts@trXMinF, p0_opts@trXMaxF, 72)
  p0_opts@tmXBLabels                   = (/"T+48", "T+54", "T+60", "T+66", "T+72", "T+78"/)
  p0_opts@tmXBLabelFontHeightF         = 0.0100

; Main title string and font size
  p0_opts@tiMainString      	       = ""
  p0_opts@tiMainFontHeightF            = 0.0225

; X-axis (title string; title font size; label font height)
  p0_opts@tiXAxisString                = "Time"
  p0_opts@tiXAxisFontHeightF   	       = 0.0150

; Y-axis (title string; title font size)
  if (diff .eq. 1) then 
   p0_opts@tiYAxisString               = "Vertical mass flux difference (%)"
  else
   p0_opts@tiYAxisString               = "Vertical mass flux (kg m~S~-2~N~ s~S~-1~N~)"
  end if
  p0_opts@tiYAxisFontHeightF           = 0.0125

; Dashed gridlines 
  p0_opts@tmXMajorGridLineDashPattern  = 2
  p0_opts@tmYMajorGridLineDashPattern  = 2
  p0_opts@gsnDraw		       = False
  p0_opts@gsnFrame		       = False
  p0_opts@gsnPaperOrientation	       = "Portrait"

  p0_opts@xyMarkLineMode	       = "MarkLines"
  p0_opts@xyLineColors		       = "red3"
  p0_opts@xyLineThicknesses	       = 3.0
  p0_opts@xyMarkers		       = 16
  p0_opts@xyMarkerColors      	       = "red3"
  p0_opts@xyMarkerSizes		       = 4.0
  p0_opts@xyMarkerThicknesses          = 2.0

; If we're plotting two curves, use this set of options 
  p1_opts                              = True
  p1_opts@gsLineThicknessF             = 3.0
  p1_opts@gsLineDashPattern            = 1
  p1_opts@gsLineColor                  = "black"

; Maximum tangential wind tendency 
  p2_opts                              = True
  p2_opts@xyLineThicknessF             = 4.0
  p2_opts@xyDashPatterns               = 0.0
  p2_opts@xyLineColor                  = "black"
  p2_opts@trYMaxF                      = 8.0
  p2_opts@trYMinF                      = -8.0

  p2_opts@tiYAxisString                = "Tangential wind tendency (m s~S~-1~N~ h~S~-1~N~)"
  p2_opts@tiYAxisFontHeightF           = 0.0125

; Relative vorticity ratio (19/12/2019)
  p3_opts                              = p2_opts
  p3_opts@xyLineColor                  = "royalblue1"
  p3_opts@trYMaxF                      = 1.0
  p3_opts@trYMinF                      = 0.5
  p3_opts@xyLineThicknessF             = 5.0
  p3_opts@tiYAxisString                = "Vorticity ratio"

;===============================================
; Output file path and final plotting options
;===============================================

; Output file path and type 
  output     = "$nep/nepartak/images/vc/mflux5_time_series_sc"+sc+\
 	       "_"+zstr0+"m_"+zstr1+"m_frun"+frun+"_trun"+trun+"_vrun"+vrun+"_"+p2

; Open workstation for plotting 
  wks        = gsn_open_wks(opt,output)

; Vertical mass flux 
  if (diff .eq. 1) then 
 ; Overlay the tendency of the maximum tangential wind 
   if (p2 .eq. "ten") then  
    flux_plot0 = gsn_csm_xy2(wks, time(216::), mflux_diff(216::), vten_arr, p0_opts, p2_opts)
 ; Or the relative vorticity ratio 
   elseif (p2 .eq. "vr") then 
    flux_plot0 = gsn_csm_xy2(wks, time(216::), mflux_diff(216::), vr_arr, p0_opts, p3_opts)
 ; Or overlay both (?)
   else
    print("Work in progress...")
   end if 
  else
   flux_plot0  = gsn_csm_xy(wks, time, mflux_arr(0,:), p0_opts)
   flux_plot1  = gsn_add_polyline(wks, flux_plot0, time, mflux_arr(1,:), p1_opts)
   overlay(flux_plot0, flux_plot1)
  end if 

;; Maximum tangential wind tendency 
;  ten_plot0   = gsn_add_polyline(wks, flux_plot0, time(216::), vten_arr, p2_opts)

; Create graphical object to overlay horizontal line onto plot 
  line       = new(1,graphic)

; Line resources (thickness; pattern; colour)
  opts_l     = True
  opts_l@gsLineThicknessF  = 2.0
  opts_l@gsLineDashPattern = 1  
  opts_l@gsLineColor       = "black"

; Overlay horizontal line through y=0
  line2   = gsn_add_polyline(wks, flux_plot0, (/p0_opts@trXMinF, p0_opts@trXMaxF/), \
                             (/0.0, 0.0/), opts_l)

;===============================
; Add and customise a legend
;===============================

; Display legend, control font and label height 
  lg_opts                            = True
  lg_opts@pmLegendDisplayMode        = "Always" 
  lg_opts@lgAutoManage               = False    
  lg_opts@lgLabelFont                = "Helvetica"
  lg_opts@lgLabelFontHeightF         = 0.025
  
; Labels and line colours 
  if (diff .ne. 1) then 
   lg_opts@lgLineColors              = (/p0_opts@xyLineColors, p1_opts@gsLineColor/)
   lg_opts@lgDashIndexes             = (/0, 0/)
   lsize                             = dimsizes(lg_opts@lgDashIndexes)
  end if 

; Label position within box
  lg_opts@lgLabelPosition            = "Right"
  lg_opts@lgItemPlacement            = "ExplicitPlacement"
  if (diff .ne. 1) then 
   lg_opts@lgItemPositions           = fspan(0.20, 0.80, lsize)
  end if 

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
  if (diff .ne. 1) then 
   lab_arr                           = (/zstr0+" m vertical mass flux", \
                                         zstr1+" m vertical mass flux"/)
   lgd                               = gsn_create_legend(wks, lsize, lab_arr, lg_opts)
  end if 

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

; Add legend to plot
  if (diff .ne. 1) then 
   annotate                          = gsn_add_annotation(flux_plot0, lgd, am_opts)
  end if 

; Draw the plot and overlay all additional objects 
  draw(flux_plot0)
  frame(wks)

; Print output file path to screen
  print(output) 

end
