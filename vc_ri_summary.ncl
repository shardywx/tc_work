; Calculate 1-h wind speed tendencies during all 24-h RI periods in our [4 x 12] MetUM simulations

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"

; ncl opt=\"x11\" mth=\"p2\" trb=0 tend=2 plt=\"all\" sub=1 
; t0=0 alt=1 var=\"wind\" vc_ri_summary.ncl

; 'opt'   = output file format ("pdf" or "x11")
; 'mth'   = method for calculating storm centre
; 'trb'   = turn troubleshooting on (1) or off (0)
; 'tend'  = method for calculating centred differences (2 or 1)
; 'alt'   = use original (0) or alternative (1) VC phase data
; 'plt'   = all phases ("all"); main four VC phases ("vc"); VC/No-VC ("no")
; 'sub'   = use a subset of the data points (e.g. only positive tendencies)
; 't0'    = wind speed / MSLP tendency threshold (e.g. 0,1,2), if 'sub=1'  
; 'var'   = analyse wind speed ("wind") or MSLP ("slp") tendency

begin

;======================================
; Read in the data from text files 
;======================================

; File path 
  diri  = "$ar/text/vc_ri_"

; Number of values in text files (no. of 24-h RI periods in our dataset)
  dsize = 38  

  dat  = asciiread(diri+"dat.txt",(/dsize/),"string")
  ens0 = asciiread(diri+"sim.txt",(/dsize/),"string")
  ts0  = asciiread(diri+"ts.txt",(/dsize/),"integer")
  tf0  = ts0 + 24

; Arrays to hold all tendencies (+ VC phase information)
  tsize    = 25
  tend_arr = new( (/dsize, tsize/), "float") 
  vc_arr   = new( (/dsize, tsize/), "string")

; Add metadata 
  vc_arr@description   = "VC phase"
  vc_arr@name          = "VC phase"     
  tend_arr@description = "Wind speed tendency" 
  tend_arr@name        = "Wind speed tendency"

; Number of times in each text file 
  ntimes_v   = 121
  ntimes_p   = 120

; String indicating how we calculated tendencies 
  if (tend .eq. 2) then 
   tstr = "_tend2"
  else
   tstr = ""
  end if 

; Tidy up before main loops 
  delete(diri)

; 26/09/2019 --> edited 'ml_03T00_em02_vtan_max_p2.txt' (L59; 3.426 --> 53.426)
;            --> this change reduced the maximum tendencies to ~5, rather than ~25 m/s/h

;==============
; Start loop
;==============

  do st = 0, dimsizes(dat)-1 

 ; Select the correct text file
   diri_v = "$ar/text/ml_"+dat(st)+"_"+ens0(st)+"_vtan_tend_"+mth+tstr+".txt"
   diri_p = "$ar/text/cp_new_"+dat(st)+"_"+ens0(st)+"_slp.txt"

 ; Read in either wind speed or MSLP data 
 ; 10/10/2019 edits
   if (var .eq. "wind") then 

  ; Check whether file exists, and if so, read data from it 
    if (.not. fileexists(diri_v) ) then 
     tend_all = new( ntimes_v, "float")
     print("Data missing for simulation "+ens0(st)+" ("+dat(st)+")" )
    else
     print("Working on RI period: "+st+" ("+ens0(st)+", "+dat(st)+")")
     tend_all = asciiread(diri_v, ntimes_v, "float")
    end if 

  ; Input wind speed values into 'tend_arr'
    tend_arr(st,:) = tend_all(ts0(st):tf0(st) )

   else

  ; Check whether file exists, and if so, read data from it
    if (.not. fileexists(diri_p) ) then
     mslp_all = new( ntimes_p, "float")
     tend_all = new( ntimes_p, "float")
     print("Data missing for simulation "+ens0(st)+" ("+dat(st)+")" )
    else
     print("Working on RI period: "+st+" ("+ens0(st)+", "+dat(st)+")")
     mslp_all = asciiread(diri_p, ntimes_p, "float")
     tend_all = mslp_all
    end if

  ; Calculate MSLP tendency inline 
    do it = 0, ntimes_p-1

     tp1 = min((/it+1,ntimes_p-1/))
     tm1 = max((/it-1,0/))      

     tend_all(it) = ( mslp_all(tp1) - mslp_all(tm1) ) / 2 
 
    end do 

  ; Input MSLP tendency values into 'tend_arr'
    tend_arr(st,:) = tend_all(ts0(st)-1:tf0(st)-1 )
   
   end if 

 ; Also read in corresponding VC phase information 
   diri_vc     = "$ar/text/vc_ri/"+ens0(st)+"_"+dat(st)+".txt"
   diri_vc_alt = "$ar/text/vc_ri/"+ens0(st)+"_"+dat(st)+"_alt.txt"

 ; Read in alternative VC phase data by setting 'alt=1'
 ; 04/10/2019 --> edited text files to match our Excel spreadsheet data (vc_phase_summary_export)
   if (alt .eq. 1) then 

    if (.not. fileexists(diri_vc_alt) ) then 
     vc_all = asciiread(diri_vc, tsize, "string")
    else
     vc_all = asciiread(diri_vc_alt, tsize, "string")
    end if 

   else

    vc_all = asciiread(diri_vc, tsize, "string")

   end if 

 ; Also copy corresponding VC phase information 
   vc_arr(st,:) = vc_all

 ; Troubleshooting 
   if (trb .eq. 1) then 
    printMinMax(tend_arr(st,:),False)
   end if 

; End loop over initialisation times [do st = 0, dimsizes(dat)-1]
  end do 

;=======================================
; Categorise our data before plotting
;======================================= 

; Create 1D arrays before finding indices 
  tend_1d = ndtooned(tend_arr)


; Do the same for our VC phase information 
  vc_1d   = ndtooned(vc_arr)

; Array of VC phase names corresponding to data in our text files 
  vc_name = (/"no_vc","r","r2m","m","m2r","no"/)

; Find indices corresponding to each VC phase 
  novc    = ind_resolve(ind(vc_1d .eq. vc_name(0) ), dimsizes(vc_arr) )
  ring    = ind_resolve(ind(vc_1d .eq. vc_name(1) ), dimsizes(vc_arr) )
  r2m     = ind_resolve(ind(vc_1d .eq. vc_name(2) ), dimsizes(vc_arr) )
  mono    = ind_resolve(ind(vc_1d .eq. vc_name(3) ), dimsizes(vc_arr) )
  m2r     = ind_resolve(ind(vc_1d .eq. vc_name(4) ), dimsizes(vc_arr) )
  misc    = ind_resolve(ind(vc_1d .eq. vc_name(5) ), dimsizes(vc_arr) )

; Troubleshooting (if required)
  if (trb .eq. 1) then 
   print(dimsizes((novc) )+" no-vc data points." )
   print(dimsizes((ring) )+" ringlike phase data points." )
   print(dimsizes((r2m) )+" r2m transition data points." )
   print(dimsizes((mono) )+" monopole phase data points." )
   print(dimsizes((m2r) )+" m2r transition data points." )
   print(dimsizes((misc) )+" additional data points." )
  end if 

;======================================================================== 
; For each VC phase, calculate the corresponding wind speed tendencies
;========================================================================

; Get size information from above arrays
  vc_size0    = (/ dimsizes(novc), dimsizes(ring), dimsizes(r2m), \
                   dimsizes(mono), dimsizes(m2r), dimsizes(misc) /)
  vc_size     = vc_size0(:,0)

; Define arrays to hold data points for plotting
  novc_plot   = new( (/vc_size(0)/) ,"float")
  ring_plot   = new( (/vc_size(1)/) ,"float")
  r2m_plot    = new( (/vc_size(2)/) ,"float")
  mono_plot   = new( (/vc_size(3)/) ,"float")
  m2r_plot    = new( (/vc_size(4)/) ,"float")
  misc_plot   = new( (/vc_size(5)/) ,"float")

; Follow the above steps for an array containing all VC phases
  vc_size_all = sum(vc_size(1:5) )
  vc_all_plot = new( (/vc_size_all/) ,"float")

; No vacillation cycle 
  do it = 0, vc_size(0)-1
   ind_plot = novc(it,:)
   novc_plot(it) = tend_arr(ind_plot(0), ind_plot(1) )
  end do 

; Ringlike phase 
  do it = 0, vc_size(1)-1
   ind_plot = ring(it,:)
   ring_plot(it) = tend_arr(ind_plot(0), ind_plot(1) )
  end do

; Ringlike to monopole transition
  do it = 0, vc_size(2)-1
   ind_plot = r2m(it,:)
   r2m_plot(it)  = tend_arr(ind_plot(0), ind_plot(1) )
  end do

; Monopole phase 
  do it = 0, vc_size(3)-1
   ind_plot = mono(it,:)
   mono_plot(it) = tend_arr(ind_plot(0), ind_plot(1) )
  end do

; Monopole to ringlike transition 
  do it = 0, vc_size(4)-1
   ind_plot = m2r(it,:) 
   m2r_plot(it)  = tend_arr(ind_plot(0), ind_plot(1) )
  end do

; Additional points (no clear VC phase)
  do it = 0, vc_size(5)-1
   ind_plot = misc(it,:)
   misc_plot(it) = tend_arr(ind_plot(0), ind_plot(1) )
  end do

; Work out array indices for following calculation 
  v1 = vc_size(1) + vc_size(2)
  v2 = v1 + vc_size(3)
  v3 = v2 + vc_size(4)
  v4 = v3 + vc_size(5)

; Finally, add all these VC-phase data points into a single array 
  vc_all_plot(0:vc_size(1)-1)  = ring_plot(:)
  vc_all_plot(vc_size(1):v1-1) = r2m_plot(:)
  vc_all_plot(v1:v2-1)         = mono_plot(:)
  vc_all_plot(v2:v3-1)         = m2r_plot(:)
  vc_all_plot(v3:v4-1)         = misc_plot(:)

;===============================================================
; If desired, remove negative wind speed tendency data points
;===============================================================

  if (sub .eq. 1) then 
   if (var .eq. "wind") then 
    novc_plot   = where(novc_plot .lt. t0, novc_plot@_FillValue, novc_plot)
    ring_plot   = where(ring_plot .lt. t0, ring_plot@_FillValue, ring_plot)
    r2m_plot    = where(r2m_plot .lt. t0, r2m_plot@_FillValue, r2m_plot)
    mono_plot   = where(mono_plot .lt. t0, mono_plot@_FillValue, mono_plot)
    m2r_plot    = where(m2r_plot .lt. t0, m2r_plot@_FillValue, m2r_plot)
    misc_plot   = where(misc_plot .lt. t0, misc_plot@_FillValue, misc_plot)
    vc_all_plot = where(vc_all_plot .lt. t0, vc_all_plot@_FillValue, vc_all_plot) 
   else
    novc_plot   = where(novc_plot .gt. t0, novc_plot@_FillValue, novc_plot)
    ring_plot   = where(ring_plot .gt. t0, ring_plot@_FillValue, ring_plot)
    r2m_plot    = where(r2m_plot .gt. t0, r2m_plot@_FillValue, r2m_plot)
    mono_plot   = where(mono_plot .gt. t0, mono_plot@_FillValue, mono_plot)
    m2r_plot    = where(m2r_plot .gt. t0, m2r_plot@_FillValue, m2r_plot)
    misc_plot   = where(misc_plot .gt. t0, misc_plot@_FillValue, misc_plot)
    vc_all_plot = where(vc_all_plot .gt. t0, vc_all_plot@_FillValue, vc_all_plot)
   end if 
  end if 

;==================================
; Customise and produce box plot 
;==================================

; Create array to hold all box plot values
; 6 VC 'phases' ['novc', 'ring', 'r2m', 'mono', 'm2r', 'none', 'all_vc']
; 5 values per phase [minimum, lower quartile, median, upper quartile, maximum] 
  box_y      = new( (/7,5/), "float")  

; Options for 'stat_dispersion' function
  opts       = True
  opts@PrintStat = False

; 01/10/2019 --> ADD CODE HERE TO CHANGE TYPE OF PLOT, ETC 

; Calculate statistics for all phases using inbuilt function 
  stat_arr   = new( (/7,30/), "float")
  bsize = dimsizes(box_y)

; 01/10/2019 --> can't create single, larger array to hold novc_plot b/c sizes are all different  
  stat_arr(0,:) = stat_dispersion(novc_plot, opts)
  stat_arr(1,:) = stat_dispersion(ring_plot, opts)
  stat_arr(2,:) = stat_dispersion(r2m_plot, opts)
  stat_arr(3,:) = stat_dispersion(mono_plot, opts)
  stat_arr(4,:) = stat_dispersion(m2r_plot, opts)
  stat_arr(5,:) = stat_dispersion(misc_plot, opts)
  stat_arr(6,:) = stat_dispersion(vc_all_plot, opts)

; Fill 'box_y' with selected values [min, lower quartile, median, upper quartile, max]
  box_y(:,0)    = stat_arr(:,2)
  box_y(:,1)    = stat_arr(:,6)
  box_y(:,2)    = stat_arr(:,8)
  box_y(:,3)    = stat_arr(:,10)
  box_y(:,4)    = stat_arr(:,14)

;========================
; Produce the box plot 
;========================

; Set output file path 
; 01/10/2019 --> FUTURE SCRIPT EDIT: check if directory exists, and create one if not 
  if (sub .eq. 0) then 
   output = "$sam/nepartak/images/vc_ri/boxplot_"+var+"_tend"+tend+"_"+plt+"_alt"+alt
  else 
   output = "$sam/nepartak/images/vc_ri/boxplot_"+var+"_tend"+tend+"_"+plt+\
            "_alt"+alt+"_sub"+sub+"_"+t0+"_t0"
  end if 
  wks = gsn_open_wks(opt, output)

; Create set of integer values denoting box label positions 
  if (plt .eq. "all") then 
   x0  = ispan(1, bsize(0)-1, 1)
  elseif (plt .eq. "vc") then 
   x0  = ispan(1, 4, 1)
  else
   x0  = (/1, 2/)
  end if 

; General box plot resources  
  pres = True 
  if (plt .eq. "all") then 
   pres@tmXBLabels = (/"No VC","Ring","R to M","Mono","M to R", "Misc"/)  
  elseif (plt .eq. "vc") then 
   pres@tmXBLabels = (/"Ring","R to M","Mono","M to R"/)
  else
   pres@tmXBLabels = (/"No VC","All VC"/)
  end if 
  pres@gsnMaximize        = True
  pres@tiYAxisFontHeightF = 0.020

  if (var .eq. "wind") then 
   pres@tmYLLabels  = (/"Tangential wind tendency (m s~S~-1~N~ h~S~-1~N~)"/)
  else
   pres@tmYLLabels  = (/"MSLP tendency (hPa h~S~-1~N~)"/)
  end if 

; Box resources 
  bres = True
  bres@boxWidth  = 0.5

; Line resources 
  lres = False

; Produce the boxplot    
; Rightmost dimension needs to contain box plot reference data (min,median,max,...) 
  if (plt .eq. "all") then 
   plot = boxplot(wks, x0, box_y(0:5,:), bres, pres, lres)
  elseif (plt .eq. "vc") then
   plot = boxplot(wks, x0, box_y(1:4,:), bres, pres, lres) 
  else
   plot = boxplot(wks, x0, box_y(0::6,:), bres, pres, lres)
  end if 

  draw(plot)
  frame(wks)

end
