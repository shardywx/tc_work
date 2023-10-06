; Script to reproduce Figure 6 from Nguyen et al. (2011)

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_rm.ncl"

; ncl opt=\"x11\" sc=5 thr=18 ts=30 tf=36 ti=6 vc_size_boxplot_fig.ncl

; 'opt'   = output file format (pdf, png, x11, ...)
; 'sc'    = storm centre calculation method (1-5); not currently perfect 
; 'sub'   = plot subset of VC forecasts (1) or comparison of VC and non-VC forecasts (0)
; 'ts'    = start time index for non-VC composite (30)                                      
; 'tf'    = end time index for non-VC composite (36)                                        
; 'ti'    = number of hours to include in VC composite, before the start of each VC (6)

begin

;=================================
; Read in data from text files
;=================================

; Directory path 
  diri     = "$ar/text/"

; Get file size
  fili_vc  = (/diri+"vc_"+thr+"_tm"+ti+".txt", \
               diri+"nvc_"+thr+"_"+ts+"_"+tf+".txt"/)
  rsize0   = new(2,"string")
  do il    = 0, 1
   rsize0(il) = systemfunc("wc -l < "+fili_vc(il) )
  end do
  rsize    = toint(rsize0)

; Read in data for forecasts in which VCs occur   
  vc_arr   = asciiread(fili_vc(0), rsize(0), "float")

; Read in data for forecasts in which VCs do not occur 
  nvc_arr  = asciiread(fili_vc(1), rsize(1), "float")

;======================================================                                        
; Produce box plot comparing VC and no-VC composites                                           
;======================================================                                        

; Create array to hold all box plot values                                                     
; 2 VC phases                                                                                  
; 5 values per phase [minimum, lower quartile, median, upper quartile, maximum]                
  box_y   = new( (/2,5/), "float")

; Options for 'stat_dispersion' function                                                       
  opts       = True
  opts@PrintStat = False

; Calculate stats for both VC phases                                                           
  ptb_novc   = stat_dispersion(nvc_arr, opts)
  ptb_vc     = stat_dispersion(vc_arr, opts)

; No-VC composite                                                                              
  box_y(0,0) = ptb_novc(2)
  box_y(0,1) = ptb_novc(6)
  box_y(0,2) = ptb_novc(8)
  box_y(0,3) = ptb_novc(10)
  box_y(0,4) = ptb_novc(14)

; VC composite                                                                                 
  box_y(1,0) = ptb_vc(2)
  box_y(1,1) = ptb_vc(6)
  box_y(1,2) = ptb_vc(8)
  box_y(1,3) = ptb_vc(10)
  box_y(1,4) = ptb_vc(14)

; Set variable y-axis limits
  ymax0 = (/ptb_novc(14), ptb_vc(14)/)
  ymin0 = (/ptb_novc(2), ptb_vc(2)/)

  ymax  = max(ymax0) + 1.0
  ymin  = max(ymin0) - 1.0

;========================                                                                      
; Produce the box plot                                                                         
;========================                                                                      

; Define output file path
  output = "$nep/nepartak/images/size_boxplot_vc_novc_"+thr+"_"+ts+"_"+tf+"_tm"+ti

; Open workstation for plotting                                                                
  wks = gsn_open_wks(opt, output)

; Create set of integer values denoting box label positions                                    
  x0  = ispan(1,2,1)

; General box plot resources                                                                   
  p_res = True
  p_res@tmXBLabels      = (/"Non-VC","VC"/)
  p_res@gsnMaximize     = True
  p_res@tiYAxisFontHeightF = 0.0225
  p_res@tiYAxisString   = "Outer radius of TS-strength winds (~S~o~N~)"

; Move y-axis title slightly left (away from plot)
  p_res@tiYAxisOffsetXF  = -0.01

; Set range (27/11/2019)
  p_res@trYMaxF  = 3.0
  p_res@trYMinF  = 0

; Box resources                                                                                
  bres = True
  bres@boxWidth  = 0.5
  bres@boxColors = (/"royalblue1", "purple3"/)

; Line resources                                                                               
  lres = True
  lres@gsLineThicknessF = 3.0

; Produce the boxplot                                                                          
; Rightmost dimension needs to contain box plot reference data (min,median,max,...)            
  plot = boxplot(wks, x0, box_y(:,:), bres, p_res, lres)

; Print output file path to screen
  print(output)

  draw(plot)
  frame(wks)

end