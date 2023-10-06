; Calculate tangential wind tendency and overlay locations of maximum dv/dt

; Run using:                                                                                  

; ncl opt=\"x11\" sc=5 w_check=1 grp=8 z0=0 w0=0.3 var=\"vt\" 
; mark=0 res=\"1h\" mean=0 rmw=1 tc_mbgt_cxz_lhs_scatter.ncl

; 'opt'   = output file format ("pdf" or "x11")                                               
; 'grp'   = chosen set of data points to analyse: ("4", "4a", "4b", "4c", "5", "6")
; 'cn0'   = offline storm centre: "slp", "vort", "geo", "geo_sm"
; 'mean'  = overlay position of mean maximum dv/dt
; 'rmw'   = overlay position of mean RMW 

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
  if (res .eq. "1h") then
   diri      = "/nfs/a319/earshar/02T12/em11"
  else
   diri      = "/nfs/a319/earshar/02T12/em11/vc/"
  end if

  fili_pr    = "20160702T1200Z_NPTK_4p4_L80_ra1t_em11"
  fili_c1    = systemfunc("cd "+diri+" ; ls "+fili_pr+"_pc*.nc")
  fili_c     = diri+"/"+fili_c1
  c          = addfile(fili_c,"r")

; Read in model level array
  hyb_levs   = c->hybrid_ht(:)
  num_levs   = dimsizes(hyb_levs)

  hyb_plot   = hyb_levs / 1000
  hyb_plot@units = "km"

; Define some cylindrical coordinate arrays 
  pii        = 3.14159265
  radii      = fspan(0.0, 1.0, 21)
  num_rad    = dimsizes(radii)

; Array of strings describing diagnostics                                                     
  diag_arr   = (/"lhs1","rhs1","rhs2","rhs3","rhs4","rhs5a","rhs5b","rhsT",\
                 "mean","eddy","vflux","vadv","vrad","vtan","vvel","aam"/)
  num_diag   = dimsizes(diag_arr)

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

; Define 4D array to hold data from all VC phases, and for all diagnostics 
; [4 * 16 * 21 * 59] for regular VC phase analysis; [1 * 16 * 21 * 59] for 'v' or 'slp'
  mbgt_arr  = new( (/num_vc, num_diag, num_levs, num_rad/), float)

;=========================================================================== 
; Read in 2D composite arrays [59 * 21] of velocity vector and MBGT terms
;=========================================================================== 

; Define path based on output frequency 
  if (res .eq. "1h") then
   diri_text = "$ar/text/"
  else
   diri_text = "$ar/text/5min_"
  end if

  do vc      = 0, num_vc-1
    print("Working on VC phase: "+vc_arr(vc) )
    do diag   = 0, num_diag-1
   ; Read in 2D array [59 * 21]
     if (res .eq. "1h") then 
      text_arr = asciiread(diri_text+diag_arr(diag)+"_sc"+sc+\
                           "_group"+grp+"_"+w0+"w_"+zstr+"_"+vc_arr(vc)+".txt", \
                           (/num_levs, num_rad+1/), "float")
    ; Ignore the leftmost column (index), and input everything else
      mbgt_arr(vc,diag,:,:) = text_arr(:,1:num_rad)
     else
    ; Read in 2D array [59 * 21]
      text_arr = asciiread(diri_text+diag_arr(diag)+"_sc"+sc+\
                           "_group"+grp+"_"+w0+"w_"+zstr+"_"+vc_arr(vc)+".txt", \
                           (/num_levs, num_rad/), "float")
    ; Ignore the leftmost column (index), and input everything else                             
      mbgt_arr(vc,diag,:,:) = text_arr(:,:)
     end if 
    end do
  end do

; Add metadata
  mbgt_arr!0    = "vc"
  mbgt_arr&vc   = vc_arr
  mbgt_arr!1    = "diag"
  mbgt_arr&diag = diag_arr
  mbgt_arr!2    = "lev"
  mbgt_arr&lev  = hyb_plot
  mbgt_arr!3    = "rad"
  mbgt_arr&rad  = radii

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
  do il = 0, 3
   rsize0(il) = systemfunc("wc -l < "+diri_scat_z(il) )
  end do
  rsize = toint(rsize0)

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

; Finally, calculate mean position of maximum dv/dt for each phase 
  ring_z_ave = avg(ring_max_arr(0,:) )
  ring_r_ave = avg(ring_max_arr(1,:) )

  r2m_z_ave  = avg(r2m_max_arr(0,:) )
  r2m_r_ave  = avg(r2m_max_arr(1,:) )

  mono_z_ave = avg(mono_max_arr(0,:) )
  mono_r_ave = avg(mono_max_arr(1,:) )

  m2r_z_ave  = avg(m2r_max_arr(0,:) )
  m2r_r_ave  = avg(m2r_max_arr(1,:) )

;==========================================
; Calculate mean fields etc, if required 
;==========================================

; Split into individal VC phases again before plotting 
  ring_arr      = mbgt_arr(0,:,:,:)
  r2m_arr       = mbgt_arr(1,:,:,:)
  mono_arr      = mbgt_arr(2,:,:,:)
  m2r_arr       = mbgt_arr(3,:,:,:)

; Also calculate mean radius of maximum wind for each phase

; Ring-like phase
  vtan_ring_1d  = ndtooned(ring_arr(13,:,:) )
  dsize_ring    = dimsizes(ring_arr(13,:,:) )
  ring_rmw      = ind_resolve(maxind(vtan_ring_1d),dsize_ring)
  ring_r_rmw    = radii(ring_rmw(0,1) )
  ring_z_rmw    = hyb_plot(ring_rmw(0,0) )

; Ring-like to monopole transition 
  vtan_r2m_1d   = ndtooned(r2m_arr(13,:,:) )
  dsize_r2m     = dimsizes(r2m_arr(13,:,:) )
  r2m_rmw       = ind_resolve(maxind(vtan_r2m_1d),dsize_r2m)
  r2m_r_rmw     = radii(r2m_rmw(0,1) )
  r2m_z_rmw     = hyb_plot(r2m_rmw(0,0) )

; Monopole phase 
  vtan_mono_1d  = ndtooned(mono_arr(13,:,:) )
  dsize_mono    = dimsizes(mono_arr(13,:,:) )
  mono_rmw      = ind_resolve(maxind(vtan_mono_1d),dsize_mono)
  mono_r_rmw    = radii(mono_rmw(0,1) )
  mono_z_rmw    = hyb_plot(mono_rmw(0,0) )

; Monopole to ring-like transition
  vtan_m2r_1d   = ndtooned(m2r_arr(13,:,:) )
  dsize_m2r     = dimsizes(m2r_arr(13,:,:) )
  m2r_rmw       = ind_resolve(maxind(vtan_m2r_1d),dsize_m2r)
  m2r_r_rmw     = radii(m2r_rmw(0,1) )
  m2r_z_rmw     = hyb_plot(m2r_rmw(0,0) )

;============================================
; Define output file path and plot details 
;============================================

  ; Leftover from alternative plotting options in 'tc_mbgt_tan_cxz_diff.ncl'
    opts_alt                               = False 

  ; Momentum budget terms 
    opts_mbgt                              = True
    opts_mbgt@cnFillOn                     = True
    opts_mbgt@cnLineLabelInterval          = 2.0
    opts_mbgt@cnLineLabelFontHeightF       = 0.012
    opts_mbgt@cnLineLabelBackgroundColor   = "transparent"
    opts_mbgt@cnLineLabelPlacementMode     = "constant"
    opts_mbgt@cnLinesOn                    = False
    opts_mbgt@cnInfoLabelOn                = False
    opts_mbgt@cnLevelSelectionMode         = "ExplicitLevels"
    opts_mbgt@cnFillPalette                = "BlueDarkRed18"
    opts_mbgt@cnLevels                     = (/-10.0, -5.0, -3.0, -2.0, -1.0, -0.5, -0.2, \
                                                0, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0/)
    opts_mbgt@cnFillColors                 = (/0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17/)

    opts_mbgt@gsnPaperOrientation          = "landscape"
    opts_mbgt@tiMainString                 = ""
    opts_mbgt@tiMainFontHeightF            = 0.0125
    opts_mbgt@gsnLeftString                = ""
    opts_mbgt@gsnRightString               = ""
    opts_mbgt@gsnMaximize                  = False
    opts_mbgt@lbLabelBarOn                 = False
    opts_mbgt@gsnDraw                      = False
    opts_mbgt@gsnFrame                     = False

  ; Additional plotting resources 
    opts_mbgt@tiYAxisString                = "Height (km)"
    opts_mbgt@tiXAxisString                = "Radial distance from cyclone centre (km)"
    opts_mbgt@tmXTOn                       = "False"
    opts_mbgt@tmYROn                       = "False"
    opts_mbgt@gsnMaximize                  = False  
    opts_mbgt@gsnAddCyclic                 = False
    opts_mbgt@tmXBMode                     = "Explicit"
    opts_mbgt@tmXBValues                   = (/0.0,0.25,0.5,0.75,1.0/)
    opts_mbgt@tmXBLabels                   = (/"0","25","50","75","100"/)
    opts_mbgt@tmYLMode                     = "Explicit"
    if (paper .eq. 0) then 
     opts_mbgt@tmYLValues                  = fspan(2.0, 18.0, 9)
     opts_mbgt@tmYLLabels                  = (/"2.0","4.0","6.0","8.0","10.0",\
                                               "12.0","14.0","16.0","18.0"/)
     opts_mbgt@trYMaxF                     = 17.6
     opts_mbgt@trXMaxF                     = 1.0
    else
     opts_mbgt@tmYLValues                  = fspan(1.0, 5.0, 5)
     opts_mbgt@tmYLLabels                  = (/"1.0","2.0","3.0","4.0","5.0"/)
     opts_mbgt@trYMaxF                     = 5.0
     opts_mbgt@trXMaxF                     = 0.5
    end if 
    opts_mbgt@gsnYAxisIrregular2Linear     = True       ; Linearise the y-axis
    opts_mbgt@trYMinF                      = 0.0        ; y-axis minimum value  

  ; Resources to overlay radial wind (following M18)
    opts_vrad                              = True
    opts_vrad@cnFillOn                     = False
    opts_vrad@cnLineColor                  = "black"
    opts_vrad@cnLevelSelectionMode         = "ExplicitLevels"
    opts_vrad@cnLevels                     = (/-1.2, 1.2/)
    opts_vrad@gsnContourNegLineDashPattern = 16
    opts_vrad@gsnContourPosLineDashPattern = 0
    opts_vrad@cnLineThicknessF             = 2.5

    opts_vrad@cnInfoLabelOn                = False
    opts_vrad@gsnDraw                      = False
    opts_vrad@gsnFrame                     = False
    opts_vrad@cnLineLabelsOn               = False
    opts_vrad@tiMainString                 = ""
    opts_vrad@gsnLeftString                = ""
    opts_vrad@gsnRightString               = ""

  ; Resources to overlay vertical velocity (following M18) 
    opts_vvel                              = True
    opts_vvel@cnFillOn                     = False
    opts_vvel@cnLineColor                  = "yellow"
    opts_vvel@cnLevelSelectionMode         = "ExplicitLevels"
    opts_vvel@cnLevels                     = (/0.5/)
    opts_vvel@cnLineThicknessF             = 5.0

    opts_vvel@cnInfoLabelOn                = False
    opts_vvel@gsnDraw                      = False
    opts_vvel@gsnFrame                     = False
    opts_vvel@cnLineLabelsOn               = False
    opts_vvel@tiMainString                 = ""
    opts_vvel@gsnLeftString                = ""
    opts_vvel@gsnRightString               = ""

  ; Resources to overlay zero line
    opts_zero                             = True
    opts_zero@cnFillOn                    = False
    opts_zero@cnLineColor                 = "grey50"
    opts_zero@cnLineDashPattern           = 3
    opts_zero@cnLevelSelectionMode        = "ExplicitLevels"
    opts_zero@cnLevels                    = 0.0
    opts_zero@cnLineThicknessF            = 3.0
    opts_zero@cnInfoLabelOn               = False
    opts_zero@gsnDraw                     = False
    opts_zero@gsnFrame                    = False 
    opts_zero@gsnContourZeroLineThicknessF = 2.0
    opts_zero@cnLineLabelsOn              = False
    opts_zero@tiMainString                = ""
    opts_zero@gsnLeftString               = ""
    opts_zero@gsnRightString              = "" 

  ; Radial wind 
    opts_vrad2                             = opts_mbgt
    delete(opts_vrad2@cnLevels)
    delete(opts_vrad2@cnFillColors)
    opts_vrad2@cnFillPalette               = "BlueDarkRed18"
    opts_vrad2@cnLevels                    = (/-10., -5., -3., -2., -1., -0.5, -0.2, \
                                                0, 0.2, 0.5, 1., 2., 3., 5., 10./)
    opts_vrad2@cnFillColors                = (/0,2,3,4,5,6,7,-1,-1,10,11,12,13,14,15,17/)
    opts_vrad2@lbLabelBarOn                = True
    opts_vrad2@lbBoxEndCapStyle            = "TriangleBothEnds"
    opts_vrad2@lbOrientation               = "Vertical"

  ; Tangential wind 
    opts_vtan2                             = opts_mbgt
    if (.not. opts_alt) then 
     delete(opts_vtan2@cnLevels)
     delete(opts_vtan2@cnFillColors)
     opts_vtan2@cnFillPalette               = "radar_new"
     opts_vtan2@cnLevels                    = (/5.0, 10.0, 15.0, 20.0,\
                                                25.0, 30.0, 35.0, 40.0, \
                                                45.0, 50.0, 55.0, 60.0, \
                                                65.0/)
     opts_vtan2@cnFillColors                = (/-1,0,1,2,3,4,5,6,7,8,9,10,11,12/)
    end if 
    opts_vtan2@lbLabelBarOn                = False
    opts_vtan2@lbBoxEndCapStyle            = "TriangleBothEnds"
    opts_vtan2@lbOrientation               = "Vertical"

  ; Vertical velocity 
    opts_vvel2                             = opts_mbgt
    if (.not. opts_alt) then 
     delete(opts_vvel2@cnLevels)
     delete(opts_vvel2@cnFillColors)
     opts_vvel2@cnFillPalette               = "n11"
     opts_vvel2@cnLevels                    = (/0.1, 0.2, 0.3, 0.4, 0.5, 0.6, \
                                                0.7, 0.8, 1.0, 1.2, 1.5, 2.0/)
     opts_vvel2@cnFillColors                = (/-1,0,1,2,3,4,5,6,7,8,9,10,11/)
    else
     opts_vvel2@cnLevels                    = (/-1.0, -0.5, -0.3, -0.2, -0.15, -0.1, -0.05, \
                                                0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 1.0/)
    end if 
    opts_vvel2@lbLabelBarOn                = True
    opts_vvel2@lbBoxEndCapStyle            = "TriangleBothEnds"
    opts_vvel2@lbOrientation               = "Vertical"

;=============================================================== 
; Resources to overlay panel labels without using 'gsn_panel'
;=============================================================== 

  ; Create option for alternative resources if we're plotting difference fields
     txid                        = new(4,graphic)
     amid                        = new(4,graphic)

   ; Label text resources
     txres                       = True
     txres@txPerimOn             = True
     txres@txFontHeightF         = 0.0175
     txres@txBackgroundFillColor = "White"

   ; Label position resources
   ; 10/10/2019 --> (amParallelPosF = -0.5/0.5 left/right edge; amOrthogonalPosF = -0.5 top edge) 
     amres                       = True
     amres@amParallelPosF        = 0.5
     amres@amOrthogonalPosF      = -0.5
     amres@amJust                = "TopRight"

     if (paper .eq. 0) then 
      panel_strings              = (/"a) Ring", "b) Ring to mono", "c) Mono", "d) Mono to ring"/)
     else
      panel_strings              = (/"a) Mono to ring", "b) Ring", "c) Ring to mono", "d) Mono"/)
     end if 

;====================================
; Plot tangential wind or tendency
;====================================

  ; Change output file string depending on whether we use 5-min or 1-h data 
    if (res .eq. "1h") then
     out_str = "mbgt"
    else
     out_str = "mbgt5"
    end if

  ; Also change units before plotting (12/12/2019)
    if (res .eq. "5min") then 
     r2m_arr(0,:,:)  = r2m_arr(0,:,:) * 12 
     mono_arr(0,:,:) = mono_arr(0,:,:) * 12
     m2r_arr(0,:,:)  = m2r_arr(0,:,:) * 12
    end if 

  ; Position each of the panels 
    xpos_plot = (/0.09, 0.48, 0.09, 0.48/)
    ypos_plot = (/0.90, 0.90, 0.50, 0.50/)
    size      = 0.30

  ; Output file location and type  
    output = "$nep/nepartak/images/cyl_coords/"+out_str+"_lhs_scatter_"+var+"_group"+grp+\
             "_sc"+sc+"_"+z0+"km_"+w0+"w"
    wks    = gsn_open_wks(opt,output)

  ; Where are the data for all these terms located? 
  ; (0) LHS_1, (12) VRAD, (13) VTAN, (14) VVEL

; Overlay on plots of tangential wind 
  if (var .eq. "v") then 

  ; Ringlike phase 
    if (paper .eq. 0) then 
     opts_vtan2@vpXF      = xpos_plot(0)
     opts_vtan2@vpYF      = ypos_plot(0)
    else
     opts_vtan2@vpXF      = xpos_plot(1)
     opts_vtan2@vpYF      = ypos_plot(1)
    end if 
    opts_vtan2@vpWidthF   = size
    opts_vtan2@vpHeightF  = size
    panel0      = gsn_csm_contour(wks,ring_arr(13,:,:),opts_vtan2)
    txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
    amid(0)     = gsn_add_annotation(panel0, txid(0), amres)

  ; Ringlike to monopole transition
    if (paper .eq. 0) then
     opts_vtan2@vpXF      = xpos_plot(1)
     opts_vtan2@vpYF      = ypos_plot(1)
    else
     opts_vtan2@vpXF      = xpos_plot(2)
     opts_vtan2@vpYF      = ypos_plot(2)      
    end if 
    opts_vtan2@vpWidthF   = size
    opts_vtan2@vpHeightF  = size
    panel1      = gsn_csm_contour(wks,r2m_arr(13,:,:),opts_vtan2)
    txid(1)     = gsn_create_text(wks, panel_strings(1), txres)
    amid(1)     = gsn_add_annotation(panel1, txid(1), amres)

  ; Monopole phase
    if (paper .eq. 0) then 
     opts_vtan2@vpXF      = xpos_plot(2)
     opts_vtan2@vpYF      = ypos_plot(2)
    else
     opts_vtan2@vpXF      = xpos_plot(3)
     opts_vtan2@vpYF      = ypos_plot(3)
    end if 
    opts_vtan2@vpWidthF   = size
    opts_vtan2@vpHeightF  = size
    panel2      = gsn_csm_contour(wks,mono_arr(13,:,:),opts_vtan2)
    txid(2)     = gsn_create_text(wks, panel_strings(2), txres)
    amid(2)     = gsn_add_annotation(panel2, txid(2), amres)

  ; Monopole to ringlike transition
    if (paper .eq. 0) then 
     opts_vtan2@vpXF      = xpos_plot(3)
     opts_vtan2@vpYF      = ypos_plot(3)
    else
     opts_vtan2@vpXF      = xpos_plot(0)
     opts_vtan2@vpYF      = ypos_plot(0)
    end if 
    opts_vtan2@vpWidthF   = size
    opts_vtan2@vpHeightF  = size
    panel3      = gsn_csm_contour(wks,m2r_arr(13,:,:),opts_vtan2)
    txid(3)     = gsn_create_text(wks, panel_strings(3), txres)
    amid(3)     = gsn_add_annotation(panel3, txid(3), amres)

; Overlay on plots of tangential wind tendency
  else 

  ; Ringlike phase 
    opts_mbgt@vpXF       = xpos_plot(0)
    opts_mbgt@vpYF       = ypos_plot(0)
    opts_mbgt@vpWidthF   = size
    opts_mbgt@vpHeightF  = size

    if (paper .eq. 0) then 
     panel0      = gsn_csm_contour(wks,ring_arr(0,:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,ring_arr(0,:,:),opts_zero)
     panel0_ovr  = gsn_csm_contour(wks,ring_arr(12,:,:),opts_vrad)
     panel0_ovr2 = gsn_csm_contour(wks,ring_arr(14,:,:),opts_vvel)
    else
     panel0      = gsn_csm_contour(wks,m2r_arr(0,:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,m2r_arr(0,:,:),opts_zero)
     panel0_ovr  = gsn_csm_contour(wks,m2r_arr(12,:,:),opts_vrad)
     panel0_ovr2 = gsn_csm_contour(wks,m2r_arr(14,:,:),opts_vvel)
    end if 

  ; Overlay contours
    overlay(panel0, cont0)
    overlay(panel0, panel0_ovr)
    overlay(panel0, panel0_ovr2)
    txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
    amid(0)     = gsn_add_annotation(panel0, txid(0), amres)
    delete([/cont0/])

  ; Ringlike to monopole transition
    opts_mbgt@vpXF       = xpos_plot(1)
    opts_mbgt@vpYF       = ypos_plot(1)
    opts_mbgt@vpWidthF   = size
    opts_mbgt@vpHeightF  = size

    if (paper .eq. 0) then 
     panel1      = gsn_csm_contour(wks,r2m_arr(0,:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,r2m_arr(0,:,:),opts_zero)
     panel1_ovr  = gsn_csm_contour(wks,r2m_arr(12,:,:),opts_vrad)
     panel1_ovr2 = gsn_csm_contour(wks,r2m_arr(14,:,:),opts_vvel)
    else
     panel1      = gsn_csm_contour(wks,ring_arr(0,:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,ring_arr(0,:,:),opts_zero)
     panel1_ovr  = gsn_csm_contour(wks,ring_arr(12,:,:),opts_vrad)
     panel1_ovr2 = gsn_csm_contour(wks,ring_arr(14,:,:),opts_vvel)
    end if 

  ; Overlay contours                                                              
    overlay(panel1, cont0)
    overlay(panel1, panel1_ovr)
    overlay(panel1, panel1_ovr2)
    txid(1)     = gsn_create_text(wks, panel_strings(1), txres)
    amid(1)     = gsn_add_annotation(panel1, txid(1), amres)
    delete([/cont0/])


  ; Monopole phase
    opts_mbgt@vpXF       = xpos_plot(2)
    opts_mbgt@vpYF       = ypos_plot(2)
    opts_mbgt@vpWidthF   = size
    opts_mbgt@vpHeightF  = size

    if (paper .eq. 0) then 
     panel2      = gsn_csm_contour(wks,mono_arr(0,:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,mono_arr(0,:,:),opts_zero)
     panel2_ovr  = gsn_csm_contour(wks,mono_arr(12,:,:),opts_vrad)
     panel2_ovr2 = gsn_csm_contour(wks,mono_arr(14,:,:),opts_vvel)
    else
     panel2      = gsn_csm_contour(wks,r2m_arr(0,:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,r2m_arr(0,:,:),opts_zero)
     panel2_ovr  = gsn_csm_contour(wks,r2m_arr(12,:,:),opts_vrad)
     panel2_ovr2 = gsn_csm_contour(wks,r2m_arr(14,:,:),opts_vvel)
    end if 

  ; Overlay contours                                                              
    overlay(panel2, cont0)
    overlay(panel2, panel2_ovr)
    overlay(panel2, panel2_ovr2)
    txid(2)     = gsn_create_text(wks, panel_strings(2), txres)
    amid(2)     = gsn_add_annotation(panel2, txid(2), amres)
    delete([/cont0/])


  ; Monopole to ringlike transition
    opts_mbgt@vpXF       = xpos_plot(3)
    opts_mbgt@vpYF       = ypos_plot(3)
    opts_mbgt@vpWidthF   = size
    opts_mbgt@vpHeightF  = size

    if (paper .eq. 0) then 
     panel3      = gsn_csm_contour(wks,m2r_arr(0,:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,m2r_arr(0,:,:),opts_zero)
     panel3_ovr  = gsn_csm_contour(wks,m2r_arr(12,:,:),opts_vrad)
     panel3_ovr2 = gsn_csm_contour(wks,m2r_arr(14,:,:),opts_vvel)
    else
     panel3      = gsn_csm_contour(wks,mono_arr(0,:,:),opts_mbgt)
     cont0       = gsn_csm_contour(wks,mono_arr(0,:,:),opts_zero)
     panel3_ovr  = gsn_csm_contour(wks,mono_arr(12,:,:),opts_vrad)
     panel3_ovr2 = gsn_csm_contour(wks,mono_arr(14,:,:),opts_vvel)
    end if 

  ; Overlay contours                                                              
    overlay(panel3, cont0)
    overlay(panel3, panel3_ovr)
    overlay(panel3, panel3_ovr2)
    txid(3)     = gsn_create_text(wks, panel_strings(3), txres)
    amid(3)     = gsn_add_annotation(panel3, txid(3), amres)
    delete([/cont0/])

  end if 

;==================================
; Add labelbar at bottom of page
;==================================

    lopt                        = True

  ; Define colourmap and set labelbar colours 
    if (var .eq. "v") then 

   ; Define colourmap 
     gsn_define_colormap(wks,"radar_new")
     lopt@lbFillColors          = (/0,2,3,4,5,6,7,8,9,10,11,12,13,14/)

   ; Customise labelbar (number of boxes, strings, size, etc)
     label_str                 = (/"5","10","15","20","25","30","35","40","45","50","55","60","65"/)

    else

   ; Define colourmap
     gsn_define_colormap(wks,"BlueDarkRed18")
     lopt@lbFillColors          = (/2,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19/)

   ; Customise labelbar (number of boxes, strings, size, etc) 
     label_str                  = (/"-20.0","-15.0","-10.0","-5.0","-2.0","-1.0","-0.5","0.0",\
                                     "0.5", "1.0", "2.0", "5.0", "10.0", "15.0", "20.0", "25.0"/)

    end if 

    nlabels                     = dimsizes(label_str)
    lopt@lbLabelFontHeightF     = 0.010
    lopt@txFont                 = "Helvetica"
    lopt@lbTitleOn              = True
    if (var .eq. "v") then 
     lopt@lbTitleString         = "Tangential wind ~C~ speed (m s~S~-1~N~)"
     lopt@lbTitleFontHeightF    = 0.014
    else
     lopt@lbTitleString         = "   Tangential wind ~C~ tendency (m s~S~-1~N~ h~S~-1~N~)"
     lopt@lbTitleFontHeightF    = 0.012
    end if 
    lopt@lbTitleOffsetF         = -0.02
    lopt@lbTitlePosition        = "Top"

  ; Labelbar size, orientation, etc 
    lopt@lbAutoManage           = False
    lopt@lbLabelAlignment       = "InteriorEdges"
    lopt@lbFillPattern          = "SolidFill"
    lopt@lbMonoFillPattern      = True
    lopt@lbBoxEndCapStyle       = "TriangleBothEnds"

    lopt@lbOrientation          = "vertical"
    lopt@vpWidthF               = 0.08
    lopt@vpHeightF              = 0.72
    xpos_lbar                   = 0.81
    ypos_lbar                   = 0.95

  ; Use routine to overlay labelbar
    gsn_labelbar_ndc(wks, nlabels, label_str, xpos_lbar, ypos_lbar, lopt)

;=================================
; Also add title at top of page 
;=================================

    if (var .eq. "v") then 
     plot_title = ""
    else
     plot_title = ""
    end if 

    title_res                  = True
    title_str                  = plot_title
    title_res@txFontHeightF    = 0.02
    title_res@txFont           = "Helvetica"
    xpos_title                 = 0.50
    ypos_title                 = 0.94

  ; Use routine to overlay labelbar 
    gsn_text_ndc(wks, title_str, xpos_title, ypos_title, title_res)

;====================================================================================
; Overlay markers for each phase, denoting position of maximum wind speed increase
;==================================================================================== 

    m_colours  = "black"
    m_indices  = 9
    m_sizes    = 5.0
    m_thick    = 3.0

  ; Properties would change, if resources above were arrays rather than single values  
    markers    = new( (/4, dimsizes(m_colours)/), graphic)

    if (paper .eq. 0) then 
     do i = 0, dimsizes(m_colours)-1
      mres0                    = True
       mres0@gsMarkerIndex      = m_indices(i)
       mres0@gsMarkerSizeF      = m_sizes(i)
       mres0@gsMarkerColor      = m_colours(i)
       mres0@gsMarkerThicknessF = m_thick(i)
       markers(0,i)             = gsn_add_polymarker(wks, panel0, ring_max_arr(1,:),\
                                                     ring_max_arr(0,:), mres0)
       markers(1,i)             = gsn_add_polymarker(wks, panel1, r2m_max_arr(1,:),\
                                                     r2m_max_arr(0,:), mres0)
       markers(2,i)             = gsn_add_polymarker(wks, panel2, mono_max_arr(1,:),\
                                                     mono_max_arr(0,:), mres0)
       markers(3,i)             = gsn_add_polymarker(wks, panel3, m2r_max_arr(1,:),\
                                                     m2r_max_arr(0,:), mres0)
     end do
    end if 

;===========================================================================
; Also overlay markers denoting the mean radius and height for each phase 
;===========================================================================

    if (mean .eq. 1) then 

      n_colours    = "grey45"
      n_indices    = 12
      n_sizes      = 8.0
      n_thick      = 5.0

      mean_markers = new(4, graphic)

      nres0                    = True 
      nres0@gsMarkerIndex      = n_indices
      nres0@gsMarkerSizeF      = n_sizes
      nres0@gsMarkerColor      = n_colours
      nres0@gsMarkerThicknessF = n_thick

    ; Overlay markers onto existing plot 
      mean_markers(0)          = gsn_add_polymarker(wks, panel0, ring_r_ave, ring_z_ave, nres0)
      mean_markers(1)          = gsn_add_polymarker(wks, panel1, r2m_r_ave, r2m_z_ave, nres0)
      mean_markers(2)          = gsn_add_polymarker(wks, panel2, mono_r_ave, mono_z_ave, nres0)
      mean_markers(3)          = gsn_add_polymarker(wks, panel3, m2r_r_ave, m2r_z_ave, nres0)

    end if 

;==================================================================
; If required, also overlay the RMW (from mean tangential wind)
;==================================================================

    if (rmw .eq. 1) then 

      n_colours    = "black"
      n_indices    = 12
      n_sizes      = 12.0
      n_thick      = 8.0

      mean_markers = new(4, graphic)
 
      nres0                    = True
      nres0@gsMarkerIndex      = n_indices
      nres0@gsMarkerSizeF      = n_sizes
      nres0@gsMarkerColor      = n_colours
      nres0@gsMarkerThicknessF = n_thick

    ; Overlay markers onto existing plot
      if (paper .eq. 0) then 
       mean_markers(0)         = gsn_add_polymarker(wks,panel0, ring_r_rmw, ring_z_rmw, nres0)
       mean_markers(1)         = gsn_add_polymarker(wks,panel1, r2m_r_rmw, r2m_z_rmw, nres0)
       mean_markers(2)         = gsn_add_polymarker(wks,panel2, mono_r_rmw, mono_z_rmw, nres0)
       mean_markers(3)         = gsn_add_polymarker(wks,panel3, m2r_r_rmw, m2r_z_rmw, nres0) 
      else
       mean_markers(0)         = gsn_add_polymarker(wks,panel0, m2r_r_rmw, m2r_z_rmw, nres0)
       mean_markers(1)         = gsn_add_polymarker(wks,panel1, ring_r_rmw, ring_z_rmw, nres0)
       mean_markers(2)         = gsn_add_polymarker(wks,panel2, r2m_r_rmw, r2m_z_rmw, nres0)
       mean_markers(3)         = gsn_add_polymarker(wks,panel3, mono_r_rmw, mono_z_rmw, nres0)
      end if 

    end if 

  ; Finally, draw the plot with everything overlaid (call 'draw' and 'frame')
    pres        = True 
    maximize_output(wks,pres)

  ; Print output file path to screen 
    print(output)

end 