load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"

load "$sam/ncl_func/attach_vert_axis.ncl"
load "$sam/ncl_func/attach_horiz_axis.ncl"
load "$sam/ncl_func/cartesian_axis.ncl"

; Run using 

; ncl sc=5 grp=8 w0=0.3 vc_5min_scatter.ncl

begin

;==========================
; Define all file paths 
;==========================

; Directory path 
  diri = "$ar/text/vtan5_max_"

; Maximum mean tangential wind 
  fili_max  = (/diri+"ring_group"+grp+"_sc"+sc+"_w"+w0+"_max.txt",\
                diri+"r2m_group"+grp+"_sc"+sc+"_w"+w0+"_max.txt",\
                diri+"mono_group"+grp+"_sc"+sc+"_w"+w0+"_max.txt",\
                diri+"m2r_group"+grp+"_sc"+sc+"_w"+w0+"_max.txt"/)

; Tendency of maximum mean tangential wind 
  fili_ten  = (/diri+"ring_group"+grp+"_sc"+sc+"_w"+w0+"_ten.txt",\
                diri+"r2m_group"+grp+"_sc"+sc+"_w"+w0+"_ten.txt",\
                diri+"mono_group"+grp+"_sc"+sc+"_w"+w0+"_ten.txt",\
                diri+"m2r_group"+grp+"_sc"+sc+"_w"+w0+"_ten.txt"/)

; Height of vmax
  fili_hgt  = (/diri+"ring_group"+grp+"_sc"+sc+"_w"+w0+"_hgt.txt",\
                diri+"r2m_group"+grp+"_sc"+sc+"_w"+w0+"_hgt.txt",\
                diri+"mono_group"+grp+"_sc"+sc+"_w"+w0+"_hgt.txt",\
                diri+"m2r_group"+grp+"_sc"+sc+"_w"+w0+"_hgt.txt"/)

; Radius of vmax
  fili_rad  = (/diri+"ring_group"+grp+"_sc"+sc+"_w"+w0+"_rad.txt",\
                diri+"r2m_group"+grp+"_sc"+sc+"_w"+w0+"_rad.txt",\
                diri+"mono_group"+grp+"_sc"+sc+"_w"+w0+"_rad.txt",\
                diri+"m2r_group"+grp+"_sc"+sc+"_w"+w0+"_rad.txt"/)

; Radial wind at vmax 
  fili_vrad = (/diri+"ring_group"+grp+"_sc"+sc+"_w"+w0+"_vrad.txt",\
                diri+"r2m_group"+grp+"_sc"+sc+"_w"+w0+"_vrad.txt",\
                diri+"mono_group"+grp+"_sc"+sc+"_w"+w0+"_vrad.txt",\
                diri+"m2r_group"+grp+"_sc"+sc+"_w"+w0+"_vrad.txt"/)

; Coriolis parameter at vmax (for AAM calculation)
  fili_f0   = (/diri+"ring_group"+grp+"_sc"+sc+"_w"+w0+"_f0.txt",\
                diri+"r2m_group"+grp+"_sc"+sc+"_w"+w0+"_f0.txt",\
                diri+"mono_group"+grp+"_sc"+sc+"_w"+w0+"_f0.txt",\
                diri+"m2r_group"+grp+"_sc"+sc+"_w"+w0+"_f0.txt"/)

; File size information 
  rsize0      = new(4,"string")
  do il       = 0, 3
   rsize0(il) = systemfunc("wc -l <"+fili_ten(il) )
  end do 
  rsize       = toint(rsize0)

;========================
; Now read in all data 
;========================

; Define arrays to hold all variables for each VC phase
  var_ring = new( (/6,rsize(0)/), "float")
  var_r2m  = new( (/6,rsize(1)/), "float")
  var_mono = new( (/6,rsize(2)/), "float")
  var_m2r  = new( (/6,rsize(3)/), "float")

; Read all variables into same array, for each VC phase 

; Ring-like phase 
  var_ring(0,:) = asciiread(fili_max(0),(/rsize(0)/),"float")
  var_ring(1,:) = asciiread(fili_ten(0),(/rsize(0)/),"float")
  var_ring(2,:) = asciiread(fili_hgt(0),(/rsize(0)/),"float")
  var_ring(3,:) = asciiread(fili_rad(0),(/rsize(0)/),"float")
  var_ring(4,:) = asciiread(fili_vrad(0),(/rsize(0)/),"float")
  var_ring(5,:) = asciiread(fili_f0(0),(/rsize(0)/),"float")

; Ring-like to monopole transition
  var_r2m(0,:) = asciiread(fili_max(1),(/rsize(1)/),"float")
  var_r2m(1,:) = asciiread(fili_ten(1),(/rsize(1)/),"float")
  var_r2m(2,:) = asciiread(fili_hgt(1),(/rsize(1)/),"float")
  var_r2m(3,:) = asciiread(fili_rad(1),(/rsize(1)/),"float")
  var_r2m(4,:) = asciiread(fili_vrad(1),(/rsize(1)/),"float")
  var_r2m(5,:) = asciiread(fili_f0(1),(/rsize(1)/),"float")

; Monopole phase 
  var_mono(0,:) = asciiread(fili_max(2),(/rsize(2)/),"float")
  var_mono(1,:) = asciiread(fili_ten(2),(/rsize(2)/),"float")
  var_mono(2,:) = asciiread(fili_hgt(2),(/rsize(2)/),"float")
  var_mono(3,:) = asciiread(fili_rad(2),(/rsize(2)/),"float")
  var_mono(4,:) = asciiread(fili_vrad(2),(/rsize(2)/),"float")
  var_mono(5,:) = asciiread(fili_f0(2),(/rsize(2)/),"float")

; Monopole to ring-like transition
  var_m2r(0,:) = asciiread(fili_max(3),(/rsize(3)/),"float")
  var_m2r(1,:) = asciiread(fili_ten(3),(/rsize(3)/),"float")
  var_m2r(2,:) = asciiread(fili_hgt(3),(/rsize(3)/),"float")
  var_m2r(3,:) = asciiread(fili_rad(3),(/rsize(3)/),"float")
  var_m2r(4,:) = asciiread(fili_vrad(3),(/rsize(3)/),"float")
  var_m2r(5,:) = asciiread(fili_f0(3),(/rsize(3)/),"float")

;======================================= 
; Calculate mean of desired variables 
;=======================================

; Create new arrays 
  ring_out      = new(6,"float")
  r2m_out       = new(6,"float")
  mono_out      = new(6,"float")
  m2r_out       = new(6,"float")

; Calculate mean values 
  do vr = 0, 5
   ring_out(vr) = avg(var_ring(vr,:) )
   r2m_out(vr)  = avg(var_r2m(vr,:) )
   mono_out(vr) = avg(var_mono(vr,:) )
   m2r_out(vr)  = avg(var_m2r(vr,:) )
  end do 

; Organise by variable rather than VC phase 
  max_out  = (/ring_out(0),r2m_out(0),mono_out(0),m2r_out(0)/)
  ten_out  = (/ring_out(1),r2m_out(1),mono_out(1),m2r_out(1)/)
  hgt_out  = (/ring_out(2),r2m_out(2),mono_out(2),m2r_out(2)/)
  rad_out  = (/ring_out(3),r2m_out(3),mono_out(3),m2r_out(3)/)
  vrad_out = (/ring_out(4),r2m_out(4),mono_out(4),m2r_out(4)/)
  f0_out   = (/ring_out(5),r2m_out(5),mono_out(5),m2r_out(5)/)

; Print out mean values to screen 
  print("Maximum mean tangential wind: ")
  print("Ring = "+max_out(0)+" , R2M = "+max_out(1)+" ; Mono = "+max_out(2)+" M2R = "+max_out(3))

  print("Tendency of maximum mean tangential wind: ")
  print("Ring = "+ten_out(0)+" , R2M = "+ten_out(1)+" ; Mono = "+ten_out(2)+" M2R = "+ten_out(3))

  print("Height of maximum mean tangential wind: ")
  print("Ring = "+hgt_out(0)+" , R2M = "+hgt_out(1)+" ; Mono = "+hgt_out(2)+" M2R = "+hgt_out(3))

  print("Radius of maximum mean tangential wind: ")
  print("Ring = "+rad_out(0)+" , R2M = "+rad_out(1)+" ; Mono = "+rad_out(2)+" M2R = "+rad_out(3))

  print("Radial wind at vmax: ")
  print("Ring = "+vrad_out(0)+" , R2M = "+vrad_out(1)+\
        " ; Mono = "+vrad_out(2)+" M2R = "+vrad_out(3))
  
end 