; Produce plots of radial momentum budget terms for ERCs 

; Run using:                                                                                  

; ncl opt=\"x11\" dist=2.0 nt=73 nr=11 ar=0.5 sc=5 rad0=0.5 z0=0 w0=0.3 zoom=0 paper=0
; col=3 p0=\"ring\" p1=\"mono\" typ=\"em04\" t0=79 t1=81
; ph=\"erc1\" trb=0 erc_mbgt_rad_cxz_diff.ncl   

; 'opt'   = output file format ("pdf" or "x11")                                               
; 'dist'  = size of box following storm (degrees)
; 'typ'   = ERC phase ("erc1","erc2","erc3")
; 'ar'    = distance in degrees of outer radius (1.0, 2.0, etc)                               
; 'nr'    = number of radial circles between r = 0 and r = 'ar' (21, 41, 61, etc)             
; 'nt'    = number of azimuth angles in cylindrical grid                                      
; 'cn0'   = offline storm centre: "slp", "vort", "geo", "geo_sm"                              
; 'rad0'  = radius (degrees) for final plots

; 'z0'    = 0 (0.002 / 0.005), 6 (6.128), 8 (8.011), 10 (10.624)

; Select data to plot (selects individual storms, for now)
; "em04"      = em04, 03T00
; "em11"      = em11, 03T00

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
  diag_arr   = (/"lhs1", "lhs2", "lhs3", "agrd", \
                 "rhs1", "rhs2", "rhs3", "rhs4", "rhs5a", "rhs5b", "rhsT", \
                 "mean", "eddy", "vflux", "vadv", "vrad", "vtan", "vvel", "aam"/)
  num_diag   = dimsizes(diag_arr)

; Minimum height at which the vertical velocity threshold was applied
; 15/01/2020 --> all files will probably have the option below, so no need to edit in future
  zstr      = "0.005km"

; File path         
  diri_text = "$ar/text/rad/"

; Count the number of times to analyse
  fili_vc   = systemfunc("cd "+diri_text+" ; ls vvel_*"+typ+"_T*.txt")
  num_vc    = dimsizes(fili_vc)

; Retrieve time information depending on simulation we're analysing 
  if (typ .eq. "em04") then 
   time_vc  = ispan(78,108,1)
  elseif (typ .eq. "em11") then
   time_vc  = ispan(72,102,1)
  end if  

; Define 4D array to hold data from all ERC phases, and for all diagnostics 
; [ntimes * 20 * 21 * 59] for regular ERC phase analysis
  mbgt_arr  = new( (/num_vc, num_diag, num_levs, num_rad/), float)

;===========================================================================
; Read in 2D composite arrays [59 * 21] of velocity vector and MBGT terms
;===========================================================================

; If plotting regular ERC phase data
  do erc      = 0, num_vc-1
    print("Working on time: T+"+time_vc(erc) )
    do diag   = 0, num_diag-1
   ; Read in data for each momentum budget term
     text_arr = asciiread(diri_text+diag_arr(diag)+"_sc5_"\
                          +w0+"w_"+zstr+"_"+typ+"_T"+time_vc(erc)+".txt", \
                          (/num_levs, num_rad/), "float")
   ; Input data into larger array
     mbgt_arr(erc,diag,:,:) = text_arr(:,:)      
    end do
  ; End loop over diagnostics
  end do
; End loop over ERC phase 

; Add metadata
  mbgt_arr!0    = "time"
  mbgt_arr&time = time_vc
  mbgt_arr!1    = "diag"
  mbgt_arr&diag = diag_arr
  mbgt_arr!2    = "lev"
  mbgt_arr&lev  = hyb_plot
  mbgt_arr!3    = "rad"
  mbgt_arr&rad  = radii

;==========================================
; Calculate mean fields etc, if required 
;==========================================

; Where are the data for all these terms located?

; (0) LHS_1, (1) LHS_2, (2) LHS_3, (3) AGRD, 
; (4) RHS_1, (5) RHS_2, (6) RHS_3, (7) RHS_4, (8) RHS_5A, (9) RHS_5B, (10) RHS_T    
; (11) MEAN, (12) EDDY, (13) VFLUX, (14) VADV, (15) VRAD, (16) VTAN, (17) VVEL, (18) AAM 

; Identify individual ERC phases, and calculate time averages for each 

; Find index of 'mbgt_arr&time' that matches each input time (T+79, T+84, etc)
  erc_ind = (/ind(mbgt_arr&time .eq. t0), ind(mbgt_arr&time .eq. t1)/)

; Calculate time average
  erc_arr = dim_avg_n_Wrap(mbgt_arr(erc_ind(0):erc_ind(1),:,:,:), 0)

; Select ERC phase to produce plots for 
  if (ph .eq. "erc1")
    plot_tt  = "Intensification phase"
  elseif (ph .eq. "erc2") then
    plot_tt  = "Weakening phase"
  elseif (ph .eq. "erc3") then
    plot_tt  = "Reintensification phase"
  end if   

; Add metadata to new arrays before plotting 
  erc_arr!0    = "diag"
  erc_arr&diag = diag_arr
  erc_arr!1    = "lev"
  erc_arr&lev  = hyb_plot
  erc_arr!2    = "rad"
  erc_arr&rad  = radii

;=================================================================
; Calculate RHS without agradient force terms (troubleshooting)
;================================================================= 

  rhs_test      = erc_arr(11,:,:) - erc_arr(4,:,:)
  rhs_test!0    = "lev"
  rhs_test&lev  = hyb_plot
  rhs_test!1    = "rad"
  rhs_test&rad  = radii

;======================================================
; Calculate total LHS tendency (material derivative)
;======================================================

  lhs_tot       = erc_arr(0,:,:) + erc_arr(1,:,:) + erc_arr(2,:,:) 
  lhs_tot!0     = "lev"
  lhs_tot&lev   = hyb_plot
  lhs_tot!1     = "rad"
  lhs_tot&rad   = radii

; Calculate radius of maximum tangential wind 
  vtan_erc_1d  = ndtooned(erc_arr(16,:,:) )
  dsize_erc    = dimsizes(erc_arr(16,:,:) )
  erc_rmw      = ind_resolve(maxind(vtan_erc_1d),dsize_erc)
  erc_r_rmw    = radii(erc_rmw(0,1) )
  erc_z_rmw    = hyb_plot(erc_rmw(0,0) )

;============================================
; Define output file path and plot details 
;============================================

; Create option for alternative resources if we're plotting difference fields  
  if (ph .ne. "erc1" .and. ph .ne. "erc2" .and. ph .ne. "erc3")
   opts_alt     = True 
  else
   opts_alt     = False
  end if    

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
  ; Select contour levels based upon type of output plot 
    if (typ .eq. "rad" .or. typ .eq. "novc_rad") then 
     opts_mbgt@cnFillPalette               = "BlueDarkRed18"
     opts_mbgt@cnLevels                    = (/-10., -5., -3., -2., -1., -0.5, -0.2, \
                                                 0, 0.2, 0.5, 1., 2., 3., 5., 10./)
     opts_mbgt@cnFillColors                = (/0,2,3,4,5,6,7,-1,-1,10,11,12,13,14,15,17/)  
    elseif (typ .eq. "novc_lhs") then 
      opts_mbgt@cnLevels                   = (/-15.0, -10.0, -5.0, -3.0, -2.0, -1.0, -0.5, -0.2, \
                                                0, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 15.0/)
      opts_mbgt@cnFillColors               = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17/)
    else
     opts_mbgt@cnLevels                    = (/-20.0, -15.0, -10.0, -5.0, -3.0, -2.0, -1.0, -0.5, \
                                                0, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 15.0, 20.0/)
     opts_mbgt@cnFillColors                = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17/)
    end if 

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

  ; Turn off tickmarks on top x-axis and right y-axis
    opts_mbgt@tmXTOn                       = "False"
    opts_mbgt@tmYROn                       = "False"

  ; Do not maximise the individual panels 
    opts_mbgt@gsnMaximize                  = False  
    opts_mbgt@gsnAddCyclic                 = False

  ; Set x-axis limits 
    opts_mbgt@tmXBMode                     = "Explicit"
    if (zoom .eq. 1 .or. typ .eq. "rad" .or. typ .eq. "lhs" .or. \
        typ .eq. "novc_rad" .or. typ .eq. "novc_lhs") then 
     opts_mbgt@tmXBValues                  = (/0.0,0.25,0.5/)
     opts_mbgt@tmXBLabels                  = (/"0","25","50"/)
     opts_mbgt@trXMaxF                     = 0.5
    else
     opts_mbgt@tmXBValues                  = (/0.0,0.25,0.5,0.75,1.0/)
     opts_mbgt@tmXBLabels                  = (/"0","25","50","75","100"/)
     opts_mbgt@trXMaxF                     = 1.0
    end if 

  ; Linearise the y-axis, and set min/max values
    opts_mbgt@gsnYAxisIrregular2Linear     = True
    opts_mbgt@trYMinF                      = 0.0 
    opts_mbgt@tmYLMode                     = "Explicit"

    if (zoom .eq. 1) then 
     opts_mbgt@trYMaxF                     = 8.0 
     opts_mbgt@tmYLValues                  = fspan(2.0, 8.0, 4)
     opts_mbgt@tmYLLabels                  = (/"2.0","4.0","6.0","8.0"/)
    elseif (typ .eq. "rad" .or. typ .eq. "lhs" .or. \
            typ .eq. "novc_rad" .or. typ .eq. "novc_lhs") then 
     opts_mbgt@trYMaxF                     = 5.0
     opts_mbgt@tmYLValues                  = fspan(1.0, 5.0, 5)
     opts_mbgt@tmYLLabels                  = (/"1.0","2.0","3.0","4.0", "5.0"/)
    else
     opts_mbgt@trYMaxF                     = 17.6
     opts_mbgt@tmYLValues                  = fspan(2.0, 18.0, 9)
     opts_mbgt@tmYLLabels                  = (/"2.0","4.0","6.0","8.0","10.0",\
                                               "12.0","14.0","16.0","18.0"/)
    end if 

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

  ; Resources to overlay absolute angular momentum 
    opts_aam                               = True 
    opts_aam@cnFillOn                      = False
    opts_aam@cnLineColor                   = "black"
    opts_aam@cnLevelSelectionMode          = "ExplicitLevels"
    opts_aam@cnLevels                      = (/1.0, 1.5, 2.0/)
    opts_aam@gsnContourNegLineDashPattern  = 16
    opts_aam@gsnContourPosLineDashPattern  = 0
    opts_aam@cnLineThicknessF              = 4.0

    opts_aam@cnInfoLabelOn                 = False
    opts_aam@gsnDraw                       = False
    opts_aam@gsnFrame                      = False
    opts_aam@cnLineLabelsOn                = False
    opts_aam@tiMainString                  = ""
    opts_aam@gsnLeftString                 = ""
    opts_aam@gsnRightString                = ""

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
    delete([/opts_vrad2@cnLevels,opts_vrad2@cnFillColors/])
    opts_vrad2@cnFillPalette               = "BlueDarkRed18"
    opts_vrad2@cnLevels                    = (/-10., -5., -3., -2., -1., -0.5, -0.2, \
                                                0, 0.2, 0.5, 1., 2., 3., 5., 10./)
    opts_vrad2@cnFillColors                = (/0,2,3,4,5,6,7,-1,-1,10,11,12,13,14,15,17/)
    opts_vrad2@lbLabelBarOn                = True
    opts_vrad2@lbBoxEndCapStyle            = "TriangleBothEnds"
    opts_vrad2@lbOrientation               = "Vertical"
    if (paper .eq. 0) then 
     opts_vrad2@vpXF                       = 0.05
     opts_vrad2@vpYF                       = 0.95
     opts_vrad2@vpWidthF                   = 0.16
     opts_vrad2@vpHeightF                  = 0.16
    end if 

  ; Tangential wind 
    opts_vtan2                             = opts_mbgt
    if (.not. opts_alt) then 
     delete([/opts_vtan2@cnLevels,opts_vtan2@cnFillColors/])
     opts_vtan2@cnFillPalette               = "radar_new"
     opts_vtan2@cnLevels                    = (/5.0, 10.0, 15.0, 20.0,\
                                                25.0, 30.0, 35.0, 40.0, \
                                                45.0, 50.0, 55.0, 60.0, \
                                                65.0/)
     opts_vtan2@cnFillColors                = (/-1,0,1,2,3,4,5,6,7,8,9,10,11,12/)
    end if 
    opts_vtan2@lbLabelBarOn                = True
    opts_vtan2@lbBoxEndCapStyle            = "TriangleBothEnds"
    opts_vtan2@lbOrientation               = "Vertical"
    if (paper .eq. 0) then 
     opts_vtan2@vpXF                       = 0.30
     opts_vtan2@vpYF                       = 0.95
     opts_vtan2@vpWidthF                   = 0.16
     opts_vtan2@vpHeightF                  = 0.16
    else
     opts_vtan2@vpXF                       = 0.37
     opts_vtan2@vpYF                       = 0.95
     opts_vtan2@vpWidthF                   = 0.20
     opts_vtan2@vpHeightF                  = 0.20
    end if 

  ; Vertical velocity 
    opts_vvel2                             = opts_mbgt
    if (.not. opts_alt) then 
     delete(opts_vvel2@cnLevels)
     delete(opts_vvel2@cnFillColors)
     opts_vvel2@cnFillPalette               = "n11"
     opts_vvel2@cnLevels                    = (/0.1, 0.2, 0.3, 0.4, 0.5, 0.6, \
                                                0.7, 0.8, 1.0, 1.2, 1.5, 2.0/)
     opts_vvel2@cnFillColors                = (/-1,0,1,2,3,4,5,6,7,8,9,10,11/)
    end if 
    opts_vvel2@lbLabelBarOn                = True
    opts_vvel2@lbBoxEndCapStyle            = "TriangleBothEnds"
    opts_vvel2@lbOrientation               = "Vertical"
    if (paper .eq. 0) then 
     opts_vvel2@vpXF                       = 0.55
     opts_vvel2@vpYF                       = 0.95
     opts_vvel2@vpWidthF                   = 0.16
     opts_vvel2@vpHeightF                  = 0.16
    else
     opts_vvel2@vpXF                       = 0.69
     opts_vvel2@vpYF                       = 0.95
     opts_vvel2@vpWidthF                   = 0.20
     opts_vvel2@vpHeightF                  = 0.20
    end if 

  ; Absolute angular momentum
    opts_aam2                               = opts_vtan2
    if (.not. opts_alt) then  
     delete(opts_aam2@cnLevels)
     delete(opts_aam2@cnFillColors)
     opts_aam2@cnLevels                     = (/0.6, 0.8, 1.0, 1.2, 1.5, \
                                               2.0, 2.5, 3.0, 3.5, \
                                               4.0, 4.5, 5.0/)
     opts_aam2@cnFillColors                 = (/-1,0,1,2,3,4,5,6,7,8,9,10,11/)
    end if 
    opts_aam2@lbLabelBarOn                  = True
    opts_aam2@lbBoxEndCapStyle              = "TriangleBothEnds"
    opts_aam2@lbOrientation                 = "Vertical"
    if (paper .eq. 0) then
     opts_aam2@vpXF                         = 0.80
     opts_aam2@vpYF                         = 0.95
     opts_aam2@vpWidthF                     = 0.16
     opts_aam2@vpHeightF                    = 0.16
    else
     opts_aam2@vpXF                         = 0.05
     opts_aam2@vpYF                         = 0.67
     opts_aam2@vpWidthF                     = 0.20
     opts_aam2@vpHeightF                    = 0.20
    end if 

;=============================================================== 
; Resources to overlay panel labels without using 'gsn_panel'
;=============================================================== 

  ; Create option for alternative resources if we're plotting difference fields
  ; 10/10/2019 --> add IF statement below into script 
    if (paper .eq. 0) then 

      txid                        = new(10,graphic)
      amid                        = new(10,graphic)

    ; Label text resources                                                                   
      txres                       = True
      txres@txPerimOn             = True
      txres@txFontHeightF         = 0.010
      txres@txBackgroundFillColor = "White"

    ; Label position resources
    ; 10/10/2019 --> (amParallelPosF = -0.5/0.5 left/right edge; amOrthogonalPosF = -0.5 top edge) 
      amres                       = True
      amres@amParallelPosF        = 0.5 
      amres@amOrthogonalPosF      = -0.5
      amres@amJust                = "TopRight"

    ; Panel strings 
      panel_strings               = (/"a) Radial wind", "b) Tangential wind", \
                                      "c) Vertical velocity", "d) AAM", \
                                      "e) Agradient force", "f) Mean v. advection", \
                                      "g) Eddy v. advection", "h) Diffusion", \
                                      "i) Sum of LHS terms", "j) Sum of RHS terms"/)

    end if 

;=================================
; Plot the momentum budget data
;=================================

  ; Calculating the terms from Montgomery et al. (2018) - their Eq.(8)                  
  ; LHS (1) : d{u}/dt          - Local tendency of the mean radial wind                 
  ; LHS (2) : {u} * d{u}/dr    - Mean radial advection of rad mom                       
  ; LHS (3) : {(u' * du'/dr) + (v'/r * du'/dy)}  - Mean horiz. adv. of eddy rad mom     

  ; RHS (1) : -{w} * d{u}/dz   - Mean vertical advection of mean rad mom                
  ; RHS (2) : -{w' * d{u'}/dz} - Eddy vertical advection of eddy rad mom                
  ; RHS (3) : {v}^2/r + f{v} - 1/{rho} * d{p}/dr - Mean agradient force                 
  ; RHS (4) : {v'^2/r - 1/rho * dp'/dr}          - Eddy agradient force                 
  ; RHS (5) : {F_r}            - Mean radial/vertical diffusive tendency of rad mom 

  ; Output file location and type  
    if (zoom .eq. 1) then 
      output = "$nep/nepartak/images/cyl_coords/erc_rad_cxz_"+typ+"_"\
               +t0+"_"+t1+"_zoom"
    else
      output = "$nep/nepartak/images/cyl_coords/erc_rad_cxz_"+typ+"_"\
               +t0+"_"+t1
    end if 
    wks    = gsn_open_wks(opt,output)

  ; Where are the data for all these terms located? 

  ; (0) LHS_1, (1) LHS_2, (2) LHS_3, (3) LHS_T, (4) AGRD,                                     
  ; (5) RHS_1, (6) RHS_2, (7) RHS_3, (8) RHS_4, (9) RHS_5A, (10) RHS_5B, (11) RHS_T           
  ; (12) MEAN, (13) EDDY, (14) VFLUX, (15) VADV, (16) VRAD, (17) VTAN, (18) VVEL, (19) AAM 

;============================================================================
; Tidy this part up, but customise the number/size/position of panels here 
;============================================================================

; 22/01/2020 --> this part shouldn't affect the new ERC plots 
  if (typ .eq. "rad" .or. typ .eq. "lhs" .or. \
      typ .eq. "paper1_diff" .or. typ .eq. "paper0_diff") then 
   xpos_plot = (/0.19, 0.57, 0.19, 0.57/)
   ypos_plot = (/0.90, 0.90, 0.55, 0.55/)
   xsize = 0.27
   ysize = 0.27
  elseif (typ .eq. "paper2_diff" .or. typ .eq. "novc_ring" .or. \
          typ .eq. "novc_r2m" .or. typ .eq. "novc_mono" .or. typ .eq. "novc_m2r") then 
   xpos_plot = (/0.07, 0.39, 0.71/)
   ypos_plot = (/0.97, 0.97, 0.97/)
   xsize     = 0.24
   ysize     = 0.24
  elseif (typ .eq. "novc_rad" .or. typ .eq. "novc_lhs") then 
   if (col .eq. 2) then 
    xpos_plot= (/0.22, 0.54/)
    ypos_plot= (/0.97, 0.97/)
   else
   xpos_plot = (/0.07, 0.39, 0.71/)
   ypos_plot = (/0.97, 0.97, 0.97/)
   end if 
   xsize     = 0.24
   ysize     = 0.24
  end if 

;=======================
; Start plotting data
;=======================

; Standard momentum budget plots for single VC phases or similar
  if (paper .eq. 0) then 
    print("Producing 10-panel plot...")

  ; Radial wind (filled contours) 
    panel_vrad  = gsn_csm_contour(wks,erc_arr(15,:,:),opts_vrad2) ; [u]                            
    panel_zero  = gsn_csm_contour(wks,erc_arr(15,:,:),opts_zero)  ; Zero line contour              
    zero0       = ColorNegDashZeroPosContour(panel_zero,"transparent","black","transparent")
    overlay(panel_vrad, zero0)
    txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
    amid(0)     = gsn_add_annotation(panel_vrad, txid(0), amres)

  ; Tangential wind (filled contours)                                                          
    panel_vtan  = gsn_csm_contour(wks,erc_arr(16,:,:),opts_vtan2) ; [v]                            
    txid(1)     = gsn_create_text(wks, panel_strings(1), txres)
    amid(1)     = gsn_add_annotation(panel_vtan, txid(1), amres)

  ; Vertical velocity (filled contours)                                                        
    panel_vvel  = gsn_csm_contour(wks,erc_arr(17,:,:),opts_vvel2) ; [w]                            
    txid(2)     = gsn_create_text(wks, panel_strings(2), txres)
    amid(2)     = gsn_add_annotation(panel_vvel, txid(2), amres)

  ; Absolute angular momentum (filled contours)                                                
    panel_aam   = gsn_csm_contour(wks,erc_arr(18,:,:),opts_aam2)  ; [AAM]
    txid(3)     = gsn_create_text(wks, panel_strings(3), txres)
    amid(3)     = gsn_add_annotation(panel_aam, txid(3), amres)

  ; Total agradient force (mean + eddy)
    opts_mbgt@vpXF       = 0.16
    opts_mbgt@vpYF       = 0.74
    opts_mbgt@vpWidthF   = 0.19
    opts_mbgt@vpHeightF  = 0.19
    panel0      = gsn_csm_contour(wks,erc_arr(3,:,:),opts_mbgt) ; [V_mzeta]                       
    cont0       = gsn_csm_contour(wks,erc_arr(3,:,:),opts_zero) ; Draw contours                   

    panel0_ovr  = gsn_csm_contour(wks,erc_arr(15,:,:),opts_vrad) ; [u]                             
    panel0_ovr2 = gsn_csm_contour(wks,erc_arr(17,:,:),opts_vvel) ; [w]                             

    overlay(panel0, cont0)                                   ; Overlay contours                
    overlay(panel0, panel0_ovr)
    overlay(panel0, panel0_ovr2)
    txid(4)     = gsn_create_text(wks, panel_strings(4), txres)
    amid(4)     = gsn_add_annotation(panel0, txid(4), amres)
    delete([/cont0/])

  ; Mean vertical advection
    opts_mbgt@vpXF       = 0.42 ; 0.55
    opts_mbgt@vpYF       = 0.74
    opts_mbgt@vpWidthF   = 0.19
    opts_mbgt@vpHeightF  = 0.19
    panel2      = gsn_csm_contour(wks,erc_arr(4,:,:),opts_mbgt) ; [V_ezeta]                       
    cont0       = gsn_csm_contour(wks,erc_arr(4,:,:),opts_zero) ; Draw contours                   

    panel2_ovr  = gsn_csm_contour(wks,erc_arr(15,:,:),opts_vrad) ; [U]                             
    panel2_ovr2 = gsn_csm_contour(wks,erc_arr(17,:,:),opts_vvel) ; [w]

    overlay(panel2, cont0)                                   ; Overlay contours                
    overlay(panel2, panel2_ovr)
    overlay(panel2, panel2_ovr2)
    txid(5)     = gsn_create_text(wks, panel_strings(5), txres)
    amid(5)     = gsn_add_annotation(panel2, txid(5), amres)
    delete([/cont0/])

  ; Eddy vertical advection
    opts_mbgt@vpXF       = 0.68 ; 0.80
    opts_mbgt@vpYF       = 0.74
    opts_mbgt@vpWidthF   = 0.19
    opts_mbgt@vpHeightF  = 0.19
    panel3      = gsn_csm_contour(wks,erc_arr(5,:,:),opts_mbgt) ; [V_ev]                          
    cont0       = gsn_csm_contour(wks,erc_arr(5,:,:),opts_zero) ; Draw contours                   

    panel3_ovr  = gsn_csm_contour(wks,erc_arr(15,:,:),opts_vrad) ; [U]                             
    panel3_ovr2 = gsn_csm_contour(wks,erc_arr(17,:,:),opts_vvel) ; [w]                             

    overlay(panel3, cont0)                                   ; Overlay contours                
    overlay(panel3, panel3_ovr)
    overlay(panel3, panel3_ovr2)
    txid(6)     = gsn_create_text(wks, panel_strings(6), txres)
    amid(6)     = gsn_add_annotation(panel3, txid(6), amres)
    delete([/cont0/])

  ; Diffusion                                                                        
  ; 22/01/2020 -> vertical for now; add horizontal back in later                     
    opts_mbgt@vpXF       = 0.16
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.19
    opts_mbgt@vpHeightF  = 0.19
    panel1      = gsn_csm_contour(wks,erc_arr(9,:,:),opts_mbgt) ; [V_mv]           
    cont0       = gsn_csm_contour(wks,erc_arr(9,:,:),opts_zero) ; Draw contours    

    panel1_ovr  = gsn_csm_contour(wks,erc_arr(15,:,:),opts_vrad) ; [U]              
    panel1_ovr2 = gsn_csm_contour(wks,erc_arr(17,:,:),opts_vvel) ; [w]              

    overlay(panel1, cont0)                                   ; Overlay contours      
    overlay(panel1, panel1_ovr)
    overlay(panel1, panel1_ovr2)
    txid(7)     = gsn_create_text(wks, panel_strings(7), txres)
    amid(7)     = gsn_add_annotation(panel1, txid(7), amres)
    delete([/cont0/])

  ; Total LHS 
    opts_mbgt@vpXF       = 0.42
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.19
    opts_mbgt@vpHeightF  = 0.19

    panel5      = gsn_csm_contour(wks,lhs_tot(:,:),opts_mbgt)      ; [V_dz]
    cont0       = gsn_csm_contour(wks,lhs_tot(:,:),opts_zero)      ; Draw contours
    panel5_ovr  = gsn_csm_contour(wks,erc_arr(15,:,:),opts_vrad)  ; [U]
    panel5_ovr2 = gsn_csm_contour(wks,erc_arr(17,:,:),opts_vvel)  ; [w]

    overlay(panel5, cont0)                                         ; Overlay contours               
    overlay(panel5, panel5_ovr)
    overlay(panel5, panel5_ovr2)
    txid(8)     = gsn_create_text(wks, panel_strings(8), txres)
    amid(8)     = gsn_add_annotation(panel5, txid(8), amres)
    delete([/cont0/])

  ; Total RHS
    opts_mbgt@vpXF       = 0.68
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.19
    opts_mbgt@vpHeightF  = 0.19

    if (trb .eq. 1) then 
     panel7     = gsn_csm_contour(wks,rhs_test(:,:),opts_mbgt)   ; [All V tend]
     cont0      = gsn_csm_contour(wks,rhs_test(:,:),opts_zero)   ; Draw contours
    else
     panel7     = gsn_csm_contour(wks,erc_arr(10,:,:),opts_mbgt)   ; [All V tend]
     cont0      = gsn_csm_contour(wks,erc_arr(10,:,:),opts_zero)   ; Draw contours 
    end if 
    panel7_ovr  = gsn_csm_contour(wks,erc_arr(15,:,:),opts_vrad)   ; [U]
    panel7_ovr2 = gsn_csm_contour(wks,erc_arr(17,:,:),opts_vvel)   ; [w]

    overlay(panel7, cont0)                                   ; Overlay contours                
    overlay(panel7, panel7_ovr)
    overlay(panel7, panel7_ovr2)
    txid(9)     = gsn_create_text(wks, panel_strings(9), txres)
    amid(9)     = gsn_add_annotation(panel7, txid(9), amres)
    delete([/cont0/])

 ; Same as above, but with 'paper=1' (8 rather than 12 panels)
 ; 22/01/2020 --> EDIT FROM HERE 
   elseif (paper .eq. 1) then
    print("Producing 8-panel plot...")

  ; Calculate total RHS tendency 
  ; 29/01/2020 --> need to update lines below (still relevant to old 'mbgt_tan' script
    rhs_tot     = erc_arr(1,:,:) + erc_arr(2,:,:) + erc_arr(3,:,:) + \
                  erc_arr(4,:,:) + erc_arr(6,:,:)
    rhs_tot!0   = "lev"
    rhs_tot&lev = hyb_plot
    rhs_tot!1   = "rad"
    rhs_tot&rad = radii

    if (row .eq. 3) then 

     xpos_plot = (/0.05, 0.37, 0.69, 0.05, 0.37, 0.69, 0.26, 0.55/)
     ypos_plot = (/0.95, 0.95, 0.95, 0.67, 0.67, 0.67, 0.39, 0.39/)
     xsize     = 0.20
     ysize     = 0.20

    elseif (row .eq. 2) then 

     xpos_plot = (/0.05, 0.30, 0.55, 0.80, 0.05, 0.30, 0.55, 0.80/)
     ypos_plot = (/0.95, 0.95, 0.95, 0.95, 0.72, 0.72, 0.72, 0.72/)
     xsize     = 0.16
     ysize     = 0.16

    end if 

  ; TOP ROW 

  ; Radial wind (filled contours)
    opts_vrad2@vpXF                       = xpos_plot(0)
    opts_vrad2@vpYF                       = ypos_plot(0)
    opts_vrad2@vpWidthF                   = xsize
    opts_vrad2@vpHeightF                  = ysize
    panel_vrad  = gsn_csm_contour(wks,erc_arr(12,:,:),opts_vrad2) ; [u]                            
    panel_zero  = gsn_csm_contour(wks,erc_arr(12,:,:),opts_zero)  ; Zero line contour              
    zero0       = ColorNegDashZeroPosContour(panel_zero,"transparent","black","transparent")
    overlay(panel_vrad, zero0)
    txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
    amid(0)     = gsn_add_annotation(panel_vrad, txid(0), amres)

  ; Tangential wind (filled contours)
    opts_vtan2@vpXF                       = xpos_plot(1)
    opts_vtan2@vpYF                       = ypos_plot(1)
    opts_vtan2@vpWidthF                   = xsize
    opts_vtan2@vpHeightF                  = ysize
    panel_vtan  = gsn_csm_contour(wks,erc_arr(13,:,:),opts_vtan2) ; [v]                            
    txid(1)     = gsn_create_text(wks, panel_strings(1), txres)
    amid(1)     = gsn_add_annotation(panel_vtan, txid(1), amres)

  ; Vertical velocity (filled contours)
    opts_vvel2@vpXF                       = xpos_plot(2)
    opts_vvel2@vpYF                       = ypos_plot(2)
    opts_vvel2@vpWidthF                   = xsize
    opts_vvel2@vpHeightF                  = ysize
    panel_vvel  = gsn_csm_contour(wks,erc_arr(14,:,:),opts_vvel2) ; [w]                            
    txid(2)     = gsn_create_text(wks, panel_strings(2), txres)
    amid(2)     = gsn_add_annotation(panel_vvel, txid(2), amres)

  ; MIDDLE ROW 

  ; Absolute angular momentum (filled contours)
    opts_aam2@vpXF                         = xpos_plot(3)
    opts_aam2@vpYF                         = ypos_plot(3)
    opts_aam2@vpWidthF                     = xsize
    opts_aam2@vpHeightF                    = ysize
    panel_aam   = gsn_csm_contour(wks,erc_arr(15,:,:),opts_aam2)  ; [AAM]
    txid(3)     = gsn_create_text(wks, panel_strings(3), txres)
    amid(3)     = gsn_add_annotation(panel_aam, txid(3), amres)

  ; Combined mean terms 
    opts_mbgt@vpXF       = xpos_plot(4)
    opts_mbgt@vpYF       = ypos_plot(4)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel0      = gsn_csm_contour(wks,erc_arr(8,:,:),opts_mbgt) ; [V_mzeta]
    cont0       = gsn_csm_contour(wks,erc_arr(8,:,:),opts_zero) ; Draw contours

    panel0_ovr  = gsn_csm_contour(wks,erc_arr(12,:,:),opts_vrad) ; [u]                             
    panel0_ovr2 = gsn_csm_contour(wks,erc_arr(14,:,:),opts_vvel) ; [w]                             

    overlay(panel0, cont0)                                   ; Overlay contours
    overlay(panel0, panel0_ovr)
    overlay(panel0, panel0_ovr2)
    txid(4)     = gsn_create_text(wks, panel_strings(4), txres)
    amid(4)     = gsn_add_annotation(panel0, txid(4), amres)
    delete([/cont0/])

  ; Combined eddy terms 
    opts_mbgt@vpXF       = xpos_plot(5)
    opts_mbgt@vpYF       = ypos_plot(5)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize
    panel1      = gsn_csm_contour(wks,erc_arr(9,:,:),opts_mbgt) ; [V_mv]
    cont0       = gsn_csm_contour(wks,erc_arr(9,:,:),opts_zero) ; Draw contours

    panel1_ovr  = gsn_csm_contour(wks,erc_arr(12,:,:),opts_vrad) ; [U]                             
    panel1_ovr2 = gsn_csm_contour(wks,erc_arr(14,:,:),opts_vvel) ; [w]                             

    overlay(panel1, cont0)                                   ; Overlay contours
    overlay(panel1, panel1_ovr)
    overlay(panel1, panel1_ovr2)
    txid(5)     = gsn_create_text(wks, panel_strings(5), txres)
    amid(5)     = gsn_add_annotation(panel1, txid(5), amres)
    delete([/cont0/])

  ; BOTTOM ROW 

  ; Total RHS [All V tend]
    opts_mbgt@vpXF       = xpos_plot(6)
    opts_mbgt@vpYF       = ypos_plot(6)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize

    panel6      = gsn_csm_contour(wks,rhs_tot(:,:),opts_mbgt)    ; [All V tend]
    cont0       = gsn_csm_contour(wks,rhs_tot(:,:),opts_zero)    ; Draw contours
    panel6_ovr  = gsn_csm_contour(wks,erc_arr(12,:,:),opts_vrad)   ; [U]
    panel6_ovr2 = gsn_csm_contour(wks,erc_arr(14,:,:),opts_vvel)   ; [w]

    overlay(panel6, cont0)                                   ; Overlay contours
    overlay(panel6, panel6_ovr)
    overlay(panel6, panel6_ovr2)
    txid(6)    = gsn_create_text(wks, panel_strings(6), txres)
    amid(6)    = gsn_add_annotation(panel6, txid(6), amres)
    delete([/cont0/])

  ; Local tendency of the mean tangential wind [V_t]
    opts_mbgt@vpXF       = xpos_plot(7)
    opts_mbgt@vpYF       = ypos_plot(7)
    opts_mbgt@vpWidthF   = xsize
    opts_mbgt@vpHeightF  = ysize

    panel7      = gsn_csm_contour(wks,erc_arr(0,:,:),opts_mbgt) ; [V_t]
    cont0       = gsn_csm_contour(wks,erc_arr(0,:,:),opts_zero) ; Draw contours

    panel7_ovr  = gsn_csm_contour(wks,erc_arr(12,:,:),opts_vrad) ; [U]
    panel7_ovr2 = gsn_csm_contour(wks,erc_arr(14,:,:),opts_vvel) ; [w]

    overlay(panel7, cont0)                                   ; Overlay contours
    overlay(panel7, panel7_ovr)
    overlay(panel7, panel7_ovr2)
    txid(7)    = gsn_create_text(wks, panel_strings(7), txres)
    amid(7)    = gsn_add_annotation(panel7, txid(7), amres)
    delete([/cont0/])

  end if 

;==============================================
; Overlay mean position of RMW (if required)
;==============================================

  ; Plotting resources 
    nres0                    = True
    nres0@gsMarkerIndex      = 12
    nres0@gsMarkerSizeF      = 10.0
    nres0@gsMarkerColor      = "black"
    nres0@gsMarkerThicknessF = 3.0

; If overlaying on certain type of panel plot (8 or 12 panels)
  if (typ .eq. "erc1" .or. typ .eq. "erc2" .or. typ .eq. "erc3") then

 ; Customise to change number of panels 
   mean_marker  = new(8, graphic)

   if (paper .eq. 1 .or. paper .eq. 0)

  ; Overlay marker onto existing tangential wind plot 
    if (typ .eq. "erc1") then 
     mean_marker(0)  = gsn_add_polymarker(wks, panel_vrad, erc1_r_rmw, erc1_z_rmw, nres0)
     mean_marker(1)  = gsn_add_polymarker(wks, panel_vtan, erc1_r_rmw, erc1_z_rmw, nres0)
     mean_marker(2)  = gsn_add_polymarker(wks, panel_vvel, erc1_r_rmw, erc1_z_rmw, nres0)
     mean_marker(3)  = gsn_add_polymarker(wks, panel_aam,  erc1_r_rmw, erc1_z_rmw, nres0)
     mean_marker(4)  = gsn_add_polymarker(wks, panel0,     erc1_r_rmw, erc1_z_rmw, nres0)
     mean_marker(5)  = gsn_add_polymarker(wks, panel1,     erc1_r_rmw, erc1_z_rmw, nres0)
     mean_marker(6)  = gsn_add_polymarker(wks, panel6,     erc1_r_rmw, erc1_z_rmw, nres0)
     mean_marker(7)  = gsn_add_polymarker(wks, panel7,     erc1_r_rmw, erc1_z_rmw, nres0)
    elseif (typ .eq. "erc2") then 
     mean_marker(0)  = gsn_add_polymarker(wks, panel_vrad, erc2_r_rmw, erc2_z_rmw, nres0)
     mean_marker(1)  = gsn_add_polymarker(wks, panel_vtan, erc2_r_rmw, erc2_z_rmw, nres0)
     mean_marker(2)  = gsn_add_polymarker(wks, panel_vvel, erc2_r_rmw, erc2_z_rmw, nres0)
     mean_marker(3)  = gsn_add_polymarker(wks, panel_aam,  erc2_r_rmw, erc2_z_rmw, nres0)
     mean_marker(4)  = gsn_add_polymarker(wks, panel0,     erc2_r_rmw, erc2_z_rmw, nres0)
     mean_marker(5)  = gsn_add_polymarker(wks, panel1,     erc2_r_rmw, erc2_z_rmw, nres0)
     mean_marker(6)  = gsn_add_polymarker(wks, panel6,     erc2_r_rmw, erc2_z_rmw, nres0)
     mean_marker(7)  = gsn_add_polymarker(wks, panel7,     erc2_r_rmw, erc2_z_rmw, nres0)
    elseif (typ .eq. "erc3") then 
     mean_marker(0)  = gsn_add_polymarker(wks, panel_vrad, erc3_r_rmw, erc3_z_rmw, nres0)
     mean_marker(1)  = gsn_add_polymarker(wks, panel_vtan, erc3_r_rmw, erc3_z_rmw, nres0)
     mean_marker(2)  = gsn_add_polymarker(wks, panel_vvel, erc3_r_rmw, erc3_z_rmw, nres0)
     mean_marker(3)  = gsn_add_polymarker(wks, panel_aam,  erc3_r_rmw, erc3_z_rmw, nres0)
     mean_marker(4)  = gsn_add_polymarker(wks, panel0,     erc3_r_rmw, erc3_z_rmw, nres0)
     mean_marker(5)  = gsn_add_polymarker(wks, panel1,     erc3_r_rmw, erc3_z_rmw, nres0)
     mean_marker(6)  = gsn_add_polymarker(wks, panel6,     erc3_r_rmw, erc3_z_rmw, nres0)
     mean_marker(7)  = gsn_add_polymarker(wks, panel7,     erc3_r_rmw, erc3_z_rmw, nres0)
    end if 

   end if 

; Otherwise overlay on our 3-panel radial wind plot 
  elseif (typ .eq. "rad" .or. typ .eq. "lhs") then 

   if (row .ne. 1) then 
    mean_marker    = new(3, graphic)
    mean_marker(0) = gsn_add_polymarker(wks, panel0, erc1_r_rmw, erc1_z_rmw, nres0)
    mean_marker(1) = gsn_add_polymarker(wks, panel1, erc2_r_rmw,  erc2_z_rmw, nres0)
    mean_marker(2) = gsn_add_polymarker(wks, panel2, erc3_r_rmw, erc3_z_rmw, nres0)
   else
    mean_marker    = new(1, graphic)
    mean_marker(0) = gsn_add_polymarker(wks, panel0, erc1_r_rmw, erc1_z_rmw, nres0)
   end if 

; Overlay on mean-eddy panel plot
  elseif (typ .eq. "mean_eddy") then 

   mean_marker  = new(6, graphic)

 ; Top row 
   if (p0 .eq. "erc1") then 
    mean_marker(0) = gsn_add_polymarker(wks, panel0, erc1_r_rmw, erc1_z_rmw, nres0)
    mean_marker(1) = gsn_add_polymarker(wks, panel1, erc1_r_rmw, erc1_z_rmw, nres0)
    mean_marker(2) = gsn_add_polymarker(wks, panel2, erc1_r_rmw, erc1_z_rmw, nres0)
   elseif (p0 .eq. "erc2") then 
    mean_marker(0) = gsn_add_polymarker(wks, panel0, erc2_r_rmw, erc2_z_rmw, nres0)
    mean_marker(1) = gsn_add_polymarker(wks, panel1, erc2_r_rmw, erc2_z_rmw, nres0)
    mean_marker(2) = gsn_add_polymarker(wks, panel2, erc2_r_rmw, erc2_z_rmw, nres0)
   elseif (p0 .eq. "erc3") then 
    mean_marker(0) = gsn_add_polymarker(wks, panel0, erc3_r_rmw, erc3_z_rmw, nres0)
    mean_marker(1) = gsn_add_polymarker(wks, panel1, erc3_r_rmw, erc3_z_rmw, nres0)
    mean_marker(2) = gsn_add_polymarker(wks, panel2, erc3_r_rmw, erc3_z_rmw, nres0)
   end if 

 ; Bottom row 
   if (p1 .eq. "erc1") then
    mean_marker(3) = gsn_add_polymarker(wks, panel3, erc1_r_rmw, erc1_z_rmw, nres0)
    mean_marker(4) = gsn_add_polymarker(wks, panel4, erc1_r_rmw, erc1_z_rmw, nres0)
    mean_marker(5) = gsn_add_polymarker(wks, panel5, erc1_r_rmw, erc1_z_rmw, nres0)
   elseif (p1 .eq. "erc2") then
    mean_marker(3) = gsn_add_polymarker(wks, panel3, erc2_r_rmw, erc2_z_rmw, nres0)
    mean_marker(4) = gsn_add_polymarker(wks, panel4, erc2_r_rmw, erc2_z_rmw, nres0)
    mean_marker(5) = gsn_add_polymarker(wks, panel5, erc2_r_rmw, erc2_z_rmw, nres0)
   elseif (p1 .eq. "erc3") then
    mean_marker(3) = gsn_add_polymarker(wks, panel3, erc3_r_rmw, erc3_z_rmw, nres0)
    mean_marker(4) = gsn_add_polymarker(wks, panel4, erc3_r_rmw, erc3_z_rmw, nres0)
    mean_marker(5) = gsn_add_polymarker(wks, panel5, erc3_r_rmw, erc3_z_rmw, nres0)
   end if

  end if 

;==================================
; Add labelbar at bottom of page
;==================================

    lopt                        = True
    lopt@lbAutoManage           = False
    lopt@lbOrientation          = "horizontal"
    lopt@lbLabelAlignment       = "InteriorEdges"
    lopt@lbFillPattern          = "SolidFill"
    lopt@lbMonoFillPattern      = True
    lopt@lbBoxEndCapStyle       = "TriangleBothEnds"

  ; Define colourmap and set labelbar colours 
    gsn_define_colormap(wks,"BlueDarkRed18")
    if (typ .ne. "rad" .and. typ .ne. "novc_rad") then
     lopt@lbFillColors           = (/2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19/)
    else
     lopt@lbFillColors           = (/2,4,5,6,7,8,9,0,0,12,13,14,15,16,17,19/)
    end if 

  ; Customise labelbar (number of boxes, strings, size, etc) 
    label_str                   = (/"-20.0","-15.0","-10.0","-5.0","-3.0","-2.0","-1.0",\
                                    "-0.5","0.0","0.5", "1.0", "2.0", "3.0", \
                                    "5.0", "10.0", "15.0", "20.0", "25.0"/)

    nlabels                     = dimsizes(label_str)
    lopt@lbLabelFontHeightF     = 0.010
    lopt@txFont                 = "Helvetica"
    lopt@lbTitleOn              = True
    lopt@lbTitleString          = "Momentum budget (m s~S~-1~N~ h~S~-1~N~)"
    lopt@lbTitleFontHeightF     = 0.0125
    lopt@lbTitlePosition        = "Bottom"
    
  ; If 'paper' option is active, then change labelbar size + position accordingly 
    if (ph .eq. "erc1" .or. ph .eq. "erc2" .or. ph .eq. "erc3") then
      if (paper .eq. 1) then 
        if (row .eq. 3) then 
         lopt@vpWidthF             = 0.65
         lopt@vpHeightF            = 0.08
         xpos_lbar                 = 0.16
         ypos_lbar                 = 0.10
        else
         lopt@vpWidthF             = 0.70
         lopt@vpHeightF            = 0.10
         xpos_lbar                 = 0.16
         ypos_lbar                 = 0.48
        end if 
      elseif (paper .eq. 0) then
        lopt@vpWidthF              = 0.75
        lopt@vpHeightF             = 0.12
        xpos_lbar                  = 0.12
        ypos_lbar                  = 0.25
      end if 
    end if 

  ; Use routine to overlay labelbar
    gsn_labelbar_ndc(wks, nlabels, label_str, xpos_lbar, ypos_lbar, lopt)

  ; Finally, draw the plot with everything overlaid (call 'draw' and 'frame')
    pres        = True 
    maximize_output(wks,pres)

  ; Print output file path to screen 
    print(output)

end 