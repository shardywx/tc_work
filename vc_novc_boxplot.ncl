; Why do VCs develop in some forecasts and not others? Create box plots of key variables 

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

load "$sam/ncl_func/attach_vert_axis.ncl"
load "$sam/ncl_func/attach_horiz_axis.ncl"
load "$sam/ncl_func/cartesian_axis.ncl"

; ncl opt=\"x11\" sc=5 r0=0.25 sub=1 grp=\"8\" vc_novc_boxplot.ncl

; 'opt'   = output file format ("pdf" or "x11")
; 'sc'    = method for calculating centred differences (2 or 1)
; 'sub'   = subset the VC forecasts (1) or use only all times from a given forecast (0)
; 'r0'    = radius of box for area-averaged SST calculations (0.25, 0.5, etc)

begin

;=========================================
; Now read in the data from text files 
;=========================================

; Follow notation of scripts including 'tc_mbgt_tan_cxz_paper.ncl'
; 11/10/2019. Group 5 data --> vc_arr_all(:,0). Group 6 data --> vc_arr_all(:,1)
  vc_arr_all = (/ (/56,38,51,30,22/), (/55,35,49,27,23/), \
                  (/69,39,64,37,27/), (/52,33,43,26,19/) /)

; Read in data for the VC and no-VC groups
  dat_novc  = asciiread("$nep/best_track/sst/novc_dat.txt",16,"string")
  dat_vc    = asciiread("$nep/best_track/sst/vc_dat.txt",29,"string")
  
  ens0_novc = asciiread("$nep/best_track/sst/novc_sim.txt",16,"string")
  ens0_vc   = asciiread("$nep/best_track/sst/vc_sim.txt",29,"string")

  ts0_novc  = new(dimsizes(dat_novc), integer)
  ts0_novc  = 30
  ts0_vc    = new(dimsizes(dat_vc), integer)
  ts0_vc    = 30

  tf0_novc  = ts0_novc + 60
  tf0_vc    = ts0_vc + 60

; Now create new arrays to hold all data for both VC phases
  var_novc  = new( (/16, 61/), "float")
  var_vc    = new( (/29, 61/), "float")

  dsize     = (/dimsizes(var_novc), dimsizes(var_vc)/)

; Add metadata to arrays above 
  var_novc!0    = "sim"
  var_novc!1    = "time"

  var_vc!0      = "sim"
  var_vc!1      = "time"

; Number of times in each text file 
  ntimes        = 120

; Radius of box used for area-averaged SST calculations 
  if (r0 .eq. 0.25) then 
   rstr         = "0.25deg_"
  else
   rstr         = ""  
  end if 

;=====================================
; Read in data for NO-VC forecasts
;=====================================

; Loop over all 16 forecasts
  do st = 0, dsize(0,0)-1 

 ; Define path for SST data 
   diri_sst = "$nep/best_track/sst/"+dat_novc(st)+"_"+ens0_novc(st)+"_"+rstr+"sst.txt"

 ; Read in SST data using 'asciiread'
   sst_data       = asciiread(diri_sst, ntimes, "float")
   var_novc(st,:) = (/sst_data(ts0_novc(st):tf0_novc(st) )/)

  end do  

;=================================
; Read in data for VC forecasts
;=================================

; Loop over all 29 forecasts
  do st = 0, dsize(1,0)-1

 ; Define path for SST data
   diri_sst = "$nep/best_track/sst/"+dat_vc(st)+"_"+ens0_vc(st)+"_"+rstr+"sst.txt"

 ; Read in SST data using 'asciiread'
   sst_data       = asciiread(diri_sst, ntimes, "float")
   var_vc(st,:) = (/sst_data(ts0_vc(st):tf0_vc(st) )/)

  end do

;==================================================
; Times for which VCs occurred for each forecast
;==================================================

; First delete existing arrays before creating these (no need to delete 'ens0_vc'/'ens0_novc')
  delete([/dat_vc, dat_novc, ts0_vc, ts0_novc, tf0_vc, tf0_novc/])  

; From 'vc_ri_for_ncl.xslx'
  dat_vc  = (/"02T12","02T12","02T12","02T12","02T12","02T12","02T12","02T12","02T12",\
              "03T00","03T00","03T00","03T00",\
              "03T12","03T12","03T12","03T12","03T12","03T12","03T12","03T12",\
              "04T00","04T00","04T00","04T00","04T00","04T00","04T00","04T00"/)

  ts0_vc  = (/46, 67, 53, 39, 66, 57, 34, 62, 44,\
              51, 58, 62, 41,\
              27, 64, 48, 28, 40, 55, 33, 36,\
              59, 59, 29, 37, 52, 40, 48, 35/)

  tf0_vc  = (/86, 94, 86, 66, 96, 96, 68, 96, 90,\
              85, 90, 95, 91,\
              96, 96, 77, 66, 91, 90, 72, 74,\
              91, 88, 82, 66, 76, 71, 85, 64/)

; Calculate how many times we need 
  nsim    = (tf0_vc - ts0_vc) + 1
  ntot    = sum(nsim)

; Sum of each of the values in 'nsim' for use when filling array in loop below 
  ns      = new(dimsizes(nsim),"integer")

; Calculate sum of each of the values in 'nsim' for use when filling array in loop below 
  do nt = 0, dimsizes(nsim)-1
   ns(nt) = sum(nsim(0:nt) ) - 1
  end do 

; Add one to the values above 
  ns1   = ns + 1

; Create 1D array to hold all VC-only data points 
  vc_only   = new(ntot, "float")

;==========================================================
; Read in data for the specific times VCs were occurring                                    
;==========================================================

; Loop over all 29 forecasts                                                                
  do st = 1, dsize(1,0)-1

 ; Define path for SST data                                                                 
   diri_sst  = "$nep/best_track/sst/"+dat_vc(st)+"_"+ens0_vc(st)+"_"+rstr+"sst.txt"

 ; Calculate indices for data to go in this array                                           

 ; Read in SST data using 'asciiread'                                                       
   sst_data  = asciiread(diri_sst, ntimes, "float")
   vc_only(ns1(st-1):ns(st) ) = (/sst_data(ts0_vc(st):tf0_vc(st) )/)

  end do

; Fill in gaps for st = 0                                                                   
  diri_sst         = "$nep/best_track/sst/"+dat_vc(0)+"_"+ens0_vc(0)+"_"+rstr+"sst.txt"
  sst_data         = asciiread(diri_sst, ntimes, "float")
  vc_only(0:ns(0)) = (/sst_data(ts0_vc(0):tf0_vc(0) )/)

;================================================================
; Times for which VCs were NOT occurring in the same forecasts
;================================================================

; First calculate total number of times per forecast
  npos0   = (96 - 24) + 1 

; Then work out how many non-VC times we have, using the VC data 
  nsim_no = nsim

; Loop through forecasts and subtract from 73
  do nn   = 0, dimsizes(nsim)-1
   nsim_no(nn) = npos0 - nsim(nn)
  end do 

; Finally, calculate the total number of times 
  ntot_no = sum(nsim_no)

; Sum of each of the values in 'nsim' for use when filling array in loop below
  ns_no   = new(dimsizes(nsim_no),"integer")

; Calculate sum of each of the values in 'nsim' for use when filling array in loop below
  do nt = 0, dimsizes(nsim_no)-1
   ns_no(nt) = sum(nsim_no(0:nt) ) - 1
  end do

; Add one to the values above
  ns1_no  = ns_no + 1

; Create 1D array to hold all non-VC-only data points
  novc_only = new(ntot_no, "float")

;===============================================================
; Read in data for the specific times VCs were NOT occurring
;===============================================================

; Initialise counter variable 
  ct = 0

; Loop over all 29 forecasts                                                         
  do st = 0, dsize(1,0)-1

 ; Define path for SST data                                                          
   diri_sst  = "$nep/best_track/sst/"+dat_vc(st)+"_"+ens0_vc(st)+"_"+rstr+"sst.txt"

 ; Read in SST data using 'asciiread'
   sst_data  = asciiread(diri_sst, ntimes, "float")

 ; Loop over each individual time, unlike for our VC cases 
 ; Read in data for all times VCs weren't occurring 
   do tt = 24, 96
    if (tt .lt. ts0_vc(st) .or. tt .gt. tf0_vc(st) ) then 
     novc_only(ct) = (/sst_data(tt)/)

   ; Only advance counter variable if time corresponds to non-VC
     ct = ct + 1 
    end if 

 ; End loop over individual times 
   end do

; End loop over all 29 forecasts
  end do

;=============================================
; Reshape to form 1D arrays before plotting
;=============================================

 ; Data from all forecasts in which VCs occurred 
   vc_sst   = ndtooned(var_vc(:,:) )
   vc_size  = dimsizes(vc_sst)

 ; Data from all forecasts in which VCs didn't occur 
   novc_sst = ndtooned(var_novc(:,:) )
   nvc_size = dimsizes(novc_sst)

 ; Data from all times in which VCs occured (subset of 'vc_sst')
   vc_all   = vc_only

 ; Non-VC data from all forecasts in which VCs occurred (inverse of 'vc_all')
   novc_all = novc_only

;======================================================                                      
; Produce box plot comparing VC and no-VC composites
;======================================================

; Create array to hold all box plot values                                                  
; 2 VC phases                                                                               
; 5 values per phase [minimum, lower quartile, median, upper quartile, maximum]             
  box_y    = new( (/2,5/), "float")

; Options for 'stat_dispersion' function                                                    
  opts     = True
  opts@PrintStat = False

; Calculate stats for both VC phases
  if (sub .eq. 1) then
   ptb_novc = stat_dispersion(novc_all, opts)
   ptb_vc   = stat_dispersion(vc_all, opts)
  else
   ptb_novc = stat_dispersion(novc_sst, opts)
   ptb_vc   = stat_dispersion(vc_sst, opts)
  end if 

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

  print(box_y(:,2) )
  
;========================                                                                    
; Produce the box plot                                                                       
;========================                                                                    

; Define output file path 
  if (sub .eq. 1) then 
   output = "$nep/best_track/boxplot_group_sc"+sc+"_vc_novc_sub"
  else
   output = "$nep/best_track/boxplot_group_sc"+sc+"_vc_novc"
  end if 

; Open workstation for plotting 
  wks = gsn_open_wks(opt, output)

; Create set of integer values denoting box label positions                                  
  x0  = ispan(1,2,1)

; General box plot resources                                                                 
  pres = True
  pres@tmXBLabels      = (/"No VC","VC"/)
  pres@gsnMaximize     = True
  pres@tiYAxisFontHeightF = 0.020

  pres@trYMaxF  = 32.0
  pres@trYMinF  = 28.0

; Box resources
; EDIT --> make the lines thicker
  bres = True
  bres@boxWidth  = 0.5
  bres@boxColors = (/"royalblue1", "purple3"/)

; Line resources                                                                             
  lres = False

; Produce the boxplot                                                                        
; Rightmost dimension needs to contain box plot reference data (min,median,max,...)          
  plot = boxplot(wks, x0, box_y(:,:), bres, pres, lres)

  draw(plot)
  frame(wks)

end 
