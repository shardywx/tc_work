; Calculate how many times the maximum wind speed increase occurred inside the RMW

; Run using:                                                                                  

; ncl sc=5 grp=8 z0=0 w0=0.3 rmw=0 var0=\"rz\" thr=0.0
; rthr=0.20 rthr2=0.05 zthr=10.0 vc_dvdt_rmw.ncl

; 'thr'   = distance between the RMW and dv/dt max, above/below which to ignore data points 
; 'rthr'  = ignore data points outside of this radius (e.g. 0.25)
; 'zthr'  = ignore data points above this height (e.g. 10)
; 'grp'   = chosen set of data points to analyse: ("4", "4a", "4b", "4c", "5", "6")
; 'ar'    = distance in degrees of outer radius (1.0, 2.0, etc)                               
; 'nt'    = number of azimuth angles in cylindrical grid                                      
; 'cn0'   = offline storm centre: "slp", "vort", "geo", "geo_sm"                              
; 'rmw'   = calculate subset based purely on radial distance (0) or use RMW (1)

; 'z0'    = 0 (0.002 / 0.005), 6 (6.128), 8 (8.011), 10 (10.624)

; 'pc' stream --> 'it'                                                                        
; 'centre'    --> 'it-1'                                                                      
; 'pd' stream --> 'it-2'                                                                      

; Load main NCL functions and procedures                                                      
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"

; Also load user-defined functions and procedures                                             
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_centre.ncl"
load "$sam/ncl_func/setup_cyl.ncl"
load "$sam/ncl_func/draw_cyl.ncl"

begin

;===================================                                                          
; Reading in data from text files                                                             
;===================================                                                          

; Read in any 'pc' stream file to retrieve vertical level metadata 
  diri       = "/nfs/a319/earshar/02T12/em11"
  fili_pr    = "20160702T1200Z_NPTK_4p4_L80_ra1t_em11"
  fili_c1    = systemfunc("cd "+diri+" ; ls "+fili_pr+"_pc*.nc")
  fili_c     = diri+"/"+fili_c1
  c          = addfile(fili_c,"r")
  hyb_levs   = c->hybrid_ht(:)
  num_levs   = dimsizes(hyb_levs)

  hyb_plot   = hyb_levs / 1000
  hyb_plot@units = "km"

; Define some cylindrical coordinate arrays 
  pii        = 3.14159265
  radii      = fspan(0.0, 1.0, 21)
  num_rad    = dimsizes(radii)

; Array of strings describing VC phase
  vc_arr     = (/"ring", "r2m", "mono", "m2r"/)
  num_vc     = dimsizes(vc_arr)

; Minimum height at which the threshold was applied (EDIT later if required)
; 09/10/2019 --> difference in 'zstr' results from using 'pc' vs 'pd' stream height data 
  if (z0 .eq. 0) then 
   zstr      = "0.005km"
  elseif (z0 .eq. 6) then 
   zstr      = "6.128km"
  elseif (z0 .eq. 8) then 
   zstr      = "8.011km"
  elseif (z0 .eq. 10) then 
   zstr      = "10.624km"
  end if 

;=====================================================================
; Read in height/radius of maximum wind speed tendency (overlaying)
;=====================================================================

; Location of text files 
  diri_scat_z  = (/"$ar/text/vtan_scatter_ring_group"+grp+"_sc"+sc+"_w"+w0+"_hgt.txt",\
                   "$ar/text/vtan_scatter_r2m_group"+grp+"_sc"+sc+"_w"+w0+"_hgt.txt",\
                   "$ar/text/vtan_scatter_mono_group"+grp+"_sc"+sc+"_w"+w0+"_hgt.txt",\
                   "$ar/text/vtan_scatter_m2r_group"+grp+"_sc"+sc+"_w"+w0+"_hgt.txt"/)

  diri_scat_r  = (/"$ar/text/vtan_scatter_ring_group"+grp+"_sc"+sc+"_w"+w0+"_rad.txt",\
                   "$ar/text/vtan_scatter_r2m_group"+grp+"_sc"+sc+"_w"+w0+"_rad.txt",\
                   "$ar/text/vtan_scatter_mono_group"+grp+"_sc"+sc+"_w"+w0+"_rad.txt",\
                   "$ar/text/vtan_scatter_m2r_group"+grp+"_sc"+sc+"_w"+w0+"_rad.txt"/)

; Get file size information for each VC phase using 'systemfunc'                            
  rsize0 = new(4,"string")
  do il  = 0, 3
   rsize0(il) = systemfunc("wc -l < "+diri_scat_z(il) )
  end do
  
  rsize  = toint(rsize0)

; Number of elements in each text file
  ring_max_arr = new( (/2, rsize(0)/), "float")
  r2m_max_arr  = new( (/2, rsize(1)/), "float")
  mono_max_arr = new( (/2, rsize(2)/), "float")
  m2r_max_arr  = new( (/2, rsize(3)/), "float")

; Size of above arrays 
  ring_size    = dimsizes(ring_max_arr)
  r2m_size     = dimsizes(r2m_max_arr)
  mono_size    = dimsizes(mono_max_arr)
  m2r_size     = dimsizes(m2r_max_arr)

; Read in height and radius index data for each VC phase 
  ring_max_z_arr      = asciiread( diri_scat_z(0), (/ring_size(1)/), "integer")
  ring_max_r_arr      = asciiread( diri_scat_r(0), (/ring_size(1)/), "integer")

  r2m_max_z_arr       = asciiread( diri_scat_z(1), (/r2m_size(1)/), "integer")
  r2m_max_r_arr       = asciiread( diri_scat_r(1), (/r2m_size(1)/), "integer")

  mono_max_z_arr      = asciiread( diri_scat_z(2), (/mono_size(1)/), "integer")
  mono_max_r_arr      = asciiread( diri_scat_r(2), (/mono_size(1)/), "integer") 

  m2r_max_z_arr       = asciiread( diri_scat_z(3), (/m2r_size(1)/), "integer")
  m2r_max_r_arr       = asciiread( diri_scat_r(3), (/m2r_size(1)/), "integer")

; Now retrieve actual heights and radii 

; Ringlike phase
  do tt = 0, ring_size(1)-1
   ring_max_arr(0,tt) = hyb_plot(ring_max_z_arr(tt) )
   ring_max_arr(1,tt) = radii(ring_max_r_arr(tt) )
  end do 

; Ringlike to monopole transition
  do tt = 0, r2m_size(1)-1
   r2m_max_arr(0,tt)  = hyb_plot(r2m_max_z_arr(tt) )
   r2m_max_arr(1,tt)  = radii(r2m_max_r_arr(tt) )
  end do 

; Monopole phase
  do tt = 0, mono_size(1)-1
   mono_max_arr(0,tt) = hyb_plot(mono_max_z_arr(tt) )
   mono_max_arr(1,tt) = radii(mono_max_r_arr(tt) )
  end do 

; Monopole to ringlike transition 
  do tt = 0, m2r_size(1)-1
   m2r_max_arr(0,tt)  = hyb_plot(m2r_max_z_arr(tt) )
   m2r_max_arr(1,tt)  = radii(m2r_max_r_arr(tt) )
  end do 

;==============================================
; Also read in corresponding RMW information
;============================================== 

; EDIT to read in either tangential or total wind 

  diri_rmw_z   = (/"$ar/text/rmw_z_vtan_ring_group"+grp+"_sc"+sc+".txt",\
                   "$ar/text/rmw_z_vtan_r2m_group"+grp+"_sc"+sc+".txt",\
                   "$ar/text/rmw_z_vtan_mono_group"+grp+"_sc"+sc+".txt",\
                   "$ar/text/rmw_z_vtan_m2r_group"+grp+"_sc"+sc+".txt"/)

  diri_rmw_r   = (/"$ar/text/rmw_vtan_ring_group"+grp+"_sc"+sc+".txt",\
                   "$ar/text/rmw_vtan_r2m_group"+grp+"_sc"+sc+".txt",\
                   "$ar/text/rmw_vtan_mono_group"+grp+"_sc"+sc+".txt",\
                   "$ar/text/rmw_vtan_m2r_group"+grp+"_sc"+sc+".txt"/)

; Read in RMW height and radius index data for each VC phase
  ring_rmw_z_arr      = asciiread( diri_rmw_z(0), (/rsize(0)/), "float")
  ring_rmw_r_arr      = asciiread( diri_rmw_r(0), (/rsize(0)/), "float")

  r2m_rmw_z_arr       = asciiread( diri_rmw_z(1), (/rsize(1)/), "float")
  r2m_rmw_r_arr       = asciiread( diri_rmw_r(1), (/rsize(1)/), "float")

  mono_rmw_z_arr      = asciiread( diri_rmw_z(2), (/rsize(2)/), "float")
  mono_rmw_r_arr      = asciiread( diri_rmw_r(2), (/rsize(2)/), "float")

  m2r_rmw_z_arr       = asciiread( diri_rmw_z(3), (/rsize(3)/), "float")
  m2r_rmw_r_arr       = asciiread( diri_rmw_r(3), (/rsize(3)/), "float")

; Number of elements in each text file
  ring_rmw_arr        = new( (/2, rsize(0)/), "float")
  r2m_rmw_arr         = new( (/2, rsize(1)/), "float")
  mono_rmw_arr        = new( (/2, rsize(2)/), "float")
  m2r_rmw_arr         = new( (/2, rsize(3)/), "float")

; Now retrieve actual heights and radii

; Ringlike phase 
  do tt = 0, rsize(0)-1
   ring_rmw_arr(0,tt) = tofloat( ring_rmw_z_arr(tt) ) / 1000
   ring_rmw_arr(1,tt) = ring_rmw_r_arr(tt)
  end do

; Ringlike to monopole transition
  do tt = 0, rsize(1)-1
   r2m_rmw_arr(0,tt)  = tofloat( r2m_rmw_z_arr(tt) ) / 1000
   r2m_rmw_arr(1,tt)  = r2m_rmw_r_arr(tt)
  end do

; Monopole phase
  do tt = 0, rsize(2)-1
   mono_rmw_arr(0,tt) = tofloat( mono_rmw_z_arr(tt) ) / 1000
   mono_rmw_arr(1,tt) = mono_rmw_r_arr(tt)
  end do

; Monopole to ringlike transition
  do tt = 0, rsize(3)-1
   m2r_rmw_arr(0,tt)  = tofloat( m2r_rmw_z_arr(tt) ) / 1000
   m2r_rmw_arr(1,tt)  = m2r_rmw_r_arr(tt)
  end do

;================================================================
; Calculate difference between RMW and radius of maximum dv/dt
;================================================================

; Ringlike phase
  ring_z_diff = ring_rmw_arr(0,:) - ring_max_arr(0,:)
  ring_r_diff = ring_rmw_arr(1,:) - ring_max_arr(1,:)

; Ringlike to monopole transition
  r2m_z_diff  = r2m_rmw_arr(0,:) - r2m_max_arr(0,:)
  r2m_r_diff  = r2m_rmw_arr(1,:) - r2m_max_arr(1,:)

; Monopole phase
  mono_z_diff = mono_rmw_arr(0,:) - mono_max_arr(0,:)
  mono_r_diff = mono_rmw_arr(1,:) - mono_max_arr(1,:)

; Monopole to ringlike transition
  m2r_z_diff  = m2r_rmw_arr(0,:) - m2r_max_arr(0,:)
  m2r_r_diff  = m2r_rmw_arr(1,:) - m2r_max_arr(1,:)

; Average differences between height of RMW and maximum dv/dt
  ring_z_ave  = avg(ring_z_diff)
  r2m_z_ave   = avg(r2m_z_diff)
  mono_z_ave  = avg(mono_z_diff)
  m2r_z_ave   = avg(m2r_z_diff)

; Average difference between RMW and location of maximum dv/dt for each VC phase 
  ring_r_ave  = avg(ring_r_diff)
  r2m_r_ave   = avg(r2m_r_diff)
  mono_r_ave  = avg(mono_r_diff)
  m2r_r_ave   = avg(m2r_r_diff)

; How many times does the maximum dv/dt occur some distance away from the RMW?
  ring_r_ind  = ind(ring_r_diff .lt. -thr .or. ring_r_diff .gt. thr)
  r2m_r_ind   = ind(r2m_r_diff  .lt. -thr .or. r2m_r_diff  .gt. thr)
  mono_r_ind  = ind(mono_r_diff .lt. -thr .or. mono_r_diff .gt. thr)
  m2r_r_ind   = ind(m2r_r_diff  .lt. -thr .or. m2r_r_diff  .gt. thr)

; Constrain by radius (ignore values outside a threshold radius)
  ring_r_thr  = ind(ring_max_arr(1,:) .ge. rthr .or. ring_max_arr(1,:) .le. rthr2)
  r2m_r_thr   = ind(r2m_max_arr(1,:) .ge. rthr .or. r2m_max_arr(1,:) .le. rthr2)
  mono_r_thr  = ind(mono_max_arr(1,:) .ge. rthr .or. mono_max_arr(1,:) .le. rthr2)
  m2r_r_thr   = ind(m2r_max_arr(1,:) .ge. rthr .or. m2r_max_arr(1,:) .le. rthr2)

; Now constrain by height (ignore values above a threshold height)
  ring_z_thr  = ind(ring_max_arr(0,:) .gt. zthr)
  r2m_z_thr   = ind(r2m_max_arr(0,:)  .gt. zthr)
  mono_z_thr  = ind(mono_max_arr(0,:) .gt. zthr)
  m2r_z_thr   = ind(m2r_max_arr(0,:)  .gt. zthr)

; Now constrain by both height and radius (combining both steps above)
  ring_rz     = ind(ring_max_arr(0,:) .gt. zthr .and. ring_max_arr(1,:) .ge. rthr)
  r2m_rz      = ind(r2m_max_arr(0,:) .gt. zthr .and. r2m_max_arr(1,:) .ge. rthr)
  mono_rz     = ind(mono_max_arr(0,:) .gt. zthr .and. mono_max_arr(1,:) .ge. rthr)
  m2r_rz      = ind(m2r_max_arr(0,:) .gt. zthr .and. m2r_max_arr(1,:) .ge. rthr)

; Finally, combine all these steps and add three constraints 
  ring_rrz    = ind(ring_max_arr(0,:) .gt. zthr .and. (ring_max_arr(1,:) .ge. rthr \
                    .or. ring_max_arr(1,:) .le. rthr2) )
  r2m_rrz     = ind(r2m_max_arr(0,:) .gt. zthr .and. (r2m_max_arr(1,:) .ge. rthr \
                    .or. r2m_max_arr(1,:) .le. rthr2) )
  mono_rrz    = ind(mono_max_arr(0,:) .gt. zthr .and. (mono_max_arr(1,:) .ge. rthr \
                    .or. mono_max_arr(1,:) .le. rthr2) )
  m2r_rrz     = ind(m2r_max_arr(0,:) .gt. zthr .and. (m2r_max_arr(1,:) .ge. rthr \
                    .or. m2r_max_arr(1,:) .le. rthr2) )

;===========================================================================
; Write out these data points to text files (for use in composite scripts)
;===========================================================================

; Create output string (km) from radial distance in degrees 
  rthr  = rthr * 100
  rthrO = sprintf("%0.0fkm", rthr)
  
  rthr2 = rthr2 * 100
  rthrI = sprintf("%0.0fkm", rthr2)

; Make output string a bit tidier 
  if (rthr2 .gt. 0)
   rstr_out = rthrI+"_"+rthrO 
  else
   rstr_out = rthrO
  end if 

; Write out indices of instances where max. dv/dt occurs inside the RMW 
  if (rmw .eq. 1) then 
   diri_out = "$ar/text/dvdt_rmw_±"+thr+"_sc"+sc+"_group"+grp+"_w"+w0+"_"+zstr+"_"
   asciiwrite(diri_out+"ring.txt",ring_r_ind)
   asciiwrite(diri_out+"r2m.txt",r2m_r_ind)
   asciiwrite(diri_out+"mono.txt",mono_r_ind)
   asciiwrite(diri_out+"m2r.txt",m2r_r_ind)
  else
 ; Indices where max dv/dt occurs inside some defined radius ("r")
   if (var0 .eq. "r") then 
    diri_out = "$ar/text/dvdt_r"+rstr_out+"_sc"+sc+"_group"+grp+"_w"+w0+"_"+zstr+"_"
    asciiwrite(diri_out+"ring.txt",ring_r_thr)
    asciiwrite(diri_out+"r2m.txt",r2m_r_thr)
    asciiwrite(diri_out+"mono.txt",mono_r_thr)
    asciiwrite(diri_out+"m2r.txt",m2r_r_thr)

 ; Indices where max dv/dt occurs below a defined height ("z") 
   elseif (var0 .eq. "z") then 
    diri_out = "$ar/text/dvdt_z"+zthr+"km_sc"+sc+"_group"+grp+"_w"+w0+"_"+zstr+"_"
    asciiwrite(diri_out+"ring.txt",ring_z_thr)
    asciiwrite(diri_out+"r2m.txt",r2m_z_thr)
    asciiwrite(diri_out+"mono.txt",mono_z_thr)
    asciiwrite(diri_out+"m2r.txt",m2r_z_thr)
 ; Indices where max dv/dt occurs below a defined height ("z") and inside a given radius ("r")
 ; "rrz" is as "rz", but both inner and outer radii are defined 
   elseif (var0 .eq. "rz" .or. var0 .eq. "rrz") then
    diri_out = "$ar/text/dvdt_z"+zthr+"km_r"+rstr_out+"_sc"+sc+"_group"+grp+"_w"+w0+"_"+zstr+"_"
    asciiwrite(diri_out+"ring.txt",ring_rz)
    asciiwrite(diri_out+"r2m.txt",r2m_rz)
    asciiwrite(diri_out+"mono.txt",mono_rz)
    asciiwrite(diri_out+"m2r.txt",m2r_rz)
   end if 
  end if 

end 