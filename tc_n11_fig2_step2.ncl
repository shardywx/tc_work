; Script calculates and writes max. standard deviation of vorticity
; and tangential wind tendency, to text files

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"

; ncl dat=\"02T12\" ens0=\"em11\" tend=2 mth=\"p2\" var=1 lev=\"same\" z0=500 tc_n11_fig2_step2.ncl

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'ens0'  = ensemble simulation (em00-em11)
; 'tend'  = calculate wind tendency using [t_1 - t_-1 / 2] (tend=2) or [t_1 - t_0 / 1] (tend=1)
; 'mth'   = method of in-script storm tracking ["v", "v1", "p", "p1", "p2"]
; 'var'   = tangential wind (0) or total wind (1)
; 'lev'   = maximum wind speed calculated on single model level ("same") or any level ("diff")
; 'z0'    = if lev = "same", z0 = approx. height of chosen model level (500, 600, 720, 845 m)

; 'pc' stream --> 'it'
; 'centre'    --> 'it-1'
; 'pd' stream --> 'it-2'

begin

;==================================
; Write out values to text file 
;==================================

; Output file directory path (layer-averaged or single layer)
  diri_out = "$ar/text/ml_"+dat+"_"+ens0

  numTIMES = 121

; Read in the maximum tangential wind from text file (calculated on model levels)
; This value is the maximum tangential wind on ANY model level (level can change with time)
; 18/10/2019 --> added option to calculate total (in additional to tangential) wind tendency
  if (lev .eq. "same") then 
 ; Tangential wind 
   if (var .eq. 0) then 
    vtan_max  = asciiread(diri_out+"_v"+z0+"_max_"+mth+".txt",\  
 	 	          (/numTIMES/),"float")	
 ; Total wind 
   else
    vtan_max  = asciiread(diri_out+"_vtot"+z0+"_max_"+mth+".txt",\
                          (/numTIMES/),"float")
   end if 
  else
   vtan_max   = asciiread("$ar/text/ml_"+dat+"_"+ens0+"_vtan_max_"+mth+".txt",\
                           (/numTIMES/),"float")
  end if 

; 1D array to hold wind tendency
  vtan_tend   = new(numTIMES,float)
  vtan_tend!0 = "time"

; Calculate the tendency of the maximum tangential wind
  do it = 0, numTIMES-1

   tp1 = min((/it+1,numTIMES-1/))
   tm1 = max((/it-1,0/))  

   if (tend .eq. 2) then 
    vtan_tend(it) = (vtan_max(tp1) - vtan_max(tm1)) / 2
   elseif (tend .eq. 1) then 
    vtan_tend(it) = (vtan_max(it) - vtan_max(tm1)) / 1
   end if 
  
  end do 

; Write out tangential wind tendency (calculated from text file)
  if (lev .eq. "same") then 
 ; Tangential wind 
   if (var .eq. 0) then 
    ten_out = diri_out+"_v"+z0+"_tend_"+mth+"_tend"+tend+".txt"
 ; Total wind
   else
    ten_out = diri_out+"_vtot"+z0+"_tend_"+mth+"_tend"+tend+".txt"
   end if 
  else
   ten_out = diri_out+"_vtan_tend_"+mth+"_tend"+tend+".txt"
  end if 
  asciiwrite(ten_out,vtan_tend)

end 