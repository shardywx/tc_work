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

; ncl dat=\"02T12\" ens0=\"em11\" wrt=0 sc=5 tc_grad_ring_mono.ncl

; 'dat'   = initialisation time string ("02T12", "03T00", ...)
; 'ens0'  = ensemble simulation (em00-em11)
; 'sc'    = method of calculating storm centre (relates to previous scripts)

; 'pc' stream --> 'it'
; 'centre'    --> 'it-1'
; 'pd' stream --> 'it-2'

begin

;==================================
; Write out values to text file 
;==================================

; Output file directory path (layer-averaged or single layer)
  diri_out = "$ar/text/"+dat+"_"+ens0+"_sc"+sc+"_inner_core_grad.txt"

  numTIMES = 119

; Read in the inner-core radial profile of relative vorticity
  inner_core = asciiread("$ar/text/"+dat+"_"+ens0+"_sc"+sc+"_inner_core.txt",\
			  (/numTIMES/),"float")	

; 1D array to hold gradient of inner-core radial profile 
  core_grad   = new(numTIMES,float)
  core_grad!0 = "time"

  print("Working on simulation: "+ens0+" ("+dat+")" )

; Calculate the gradient/tendency of the inner-core radial profile 
  do it = 0, numTIMES-1

   tp1 = min((/it+1,numTIMES-1/))
   tm1 = max((/it-1,0/))  

;  print("it = "+it+" ; inner_core(tp1) = "+inner_core(tp1)+" ; inner_core(tm1) = "+inner_core(tm1))

   core_grad(it) = (inner_core(tp1) - inner_core(tm1)) / 2
  
  end do 

; Write out gradient/tendency
  if (wrt .eq. 1) then 
   asciiwrite(diri_out, core_grad)
  end if 

end 