; Script to calculate the radial momentum budget as in Montgomery et al. (2018)

; Run using:

; ncl 'opt="pdf"' 'typ="stoa"' clr=4 full=1 rad0=1 grid=0 tc_mbgt_rad_comp.ncl

; 'typ'   = input directory ("sym", "asym", "ri", "sst", "stoa", "atos")
; 'opt'   = output file format ("pdf" or "x11")
; 'clr'	  = colour map for plots (option 4 is used in Roger Smith's papers):
; 	    "ncl_default" [1],  "bdr_extra" [2], "amwg256" [3], "GMT_polar" [4]
; 'full'  = plot 14 panels inc. u,v,w,mom (2), 10 panels w/mbgt terms only (1) or 4 panels (0) 
; 'rad0'  = radius of final plots: 150 km (0) or 200 km (1)
; 'grid'  = overlay NDC grid (1) or leave alone (0)

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"

begin

;==============================
; Prelash (reading in files)
;==============================

; Important directory information 
  diri        = "$ar/text/comp/"+typ+"/"
  fili_prefix = "mbgt_rad_"

; Read in chosen terms

; Momentum budget terms (radial advection of mean radial momentum; mean agradient force)
  fili_l2     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"lhs2*.txt")
  fili_r3     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"rhs3*.txt")

; Radial wind 
  fili_vr     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"vrad*.txt")

; Information from original scripts (size of arrays)
  nlev        = 59                      ; Number of vertical (model) levels
  ntimes      = dimsizes(fili_r3)	; Number of time averages to read in 
  ndiag	      = 3			; Number of diagnostics (source terms)
  nrad	      = 41			; Number of radial circles out from storm centre

; Create 4D array to hold all incoming data
  plot_arr    = new( (/nlev, ntimes, ndiag, nrad/), "float")

; Read in radial wind data
  vr_test     = asciiread( diri+fili_vr(2), (/41,59/), "float")  

;===============================
; Add metadata before ploting 
;===============================

 ; Radial direction metadata
   arad		      = 2.0
   nrad		      = 41			
   radii	      = fspan(0.0, arad, nrad)

 ; Vertical direction (model level) metadata
   fili_pc    	      = "/nfs/a319/earshar/02T12/em00/20160702T1200Z_NPTK_4p4_L80_ra1t_em00_pc.nc"
   c		      = addfile(fili_pc,"r")
   hyb_c	      = c->hybrid_ht(:)
   hyb_levs	      = hyb_c / 1000
   hyb_levs@units     = "km"   

   vr_test!1	      = "lev"
   vr_test&lev	      = hyb_levs
   
   vr_test!0	      = "rad"
   vr_test&rad	      = radii
   vr_test@units      = "m s~S~-1~N~"
   vr_test@name	      = "Radial wind"

;==========================
; Options for plotting 
;==========================

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

    if (clr.eq.1) then					       
     opts_mbgt@cnFillPalette		   = "ncl_default"
     opts_mbgt@cnFillColors                = (/2,8,16,28,40,52,64,76,88,100,118,-1,-1,\
                                               134,150,166,178,186,194,\
                                               202,214,226,238,250/)
     opts_mbgt@cnLevels                    = (/-150., -100., -50., -20., -10., -5., \
                                                -3., -2., -1., -0.5, -0.2, \
                                                0, 0.2, 0.5, 1., 2., 3., \
                                                5., 10., 20., 50., 100., 150./)

    elseif (clr.eq.2) then 
     opts_mbgt@cnFillPalette               = "bdr_extra"
     opts_mbgt@cnFillColors                = (/0,1,2,3,4,5,6,7,8,9,10,-1,-1,\
                                               13,14,15,16,17,18,19,20,21,22,23/)
     opts_mbgt@cnLevels                    = (/-150., -100., -50., -20., -10., -5., \
                                                -3., -2., -1., -0.5, -0.2, \
                                                0, 0.2, 0.5, 1., 2., 3., \
                                                5., 10., 20., 50., 100., 150./)

    elseif (clr.eq.3) then
     opts_mbgt@cnFillPalette               = "amwg256"
     opts_mbgt@cnFillColors                = (/0,6,12,18,27,36,45,54,66,78,90,-1,-1,\
                                               148,156,166,178,186,194,\
                                               202,214,226,238,250/)
     opts_mbgt@cnLevels                    = (/-150., -100., -50., -20., -10., -5., \
                                                -3., -2., -1., -0.5, -0.2, \
                                                0, 0.2, 0.5, 1., 2., 3., \
                                                5., 10., 20., 50., 100., 150./)
    elseif (clr.eq.4) then
     opts_mbgt@cnFillPalette               = "GMT_polar"
     opts_mbgt@cnFillColors                = (/3,5,7,-1,11,14,18/)
     opts_mbgt@cnLevels                    = (/-20,-5,0,0.5,5,20/)
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

    opts_mbgt@tmXTOn                       = "False" ; Turn off top x-axis TM
    opts_mbgt@tmYROn                       = "False" ; Turn off right y-axis TM
    opts_mbgt@gsnMaximize                  = False   ; Maximise plot size
    opts_mbgt@gsnAddCyclic                 = False
    opts_mbgt@tmXBMode                     = "Explicit"  ; Set tick marks explicitly

    if (rad0 .eq. 0) then
     opts_mbgt@tmXBValues                   = (/0.0,0.25,0.5,0.75,1.0,1.25,1.5/)
     opts_mbgt@tmXBLabels                   = (/"0","25","50","75","100","125","150"/)
    else if (rad0 .eq. 1) then
     opts_mbgt@tmXBValues                   = (/0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0/)
     opts_mbgt@tmXBLabels                   = (/"0","25","50","75","100","125","150","175","200"/)
    end if
    end if

    opts_mbgt@tmYLMode			   = "Explicit"
    opts_mbgt@tmYLValues		   = fspan(2.0, 18.0, 9)
    opts_mbgt@tmYLLabels		   = (/"2.0","4.0","6.0","8.0","10.0",\
    					       "12.0","14.0","16.0","18.0"/)
    opts_mbgt@gsnYAxisIrregular2Linear	   = True	; Linearise the y-axis 

    opts_mbgt@trYMinF			   = 0.0	; y-axis minimum value
    opts_mbgt@trYMaxF			   = 17.6	; y-axis maximum value 

;==================================================================
; Also add resources for mean quantities (u,v,w) when 'full' = 2
;==================================================================

  ; Reading in colour table using 'cnFillPalette' removes the first 2 colours (black/white)
  ; The indices are thus all 2 values less than they'd normally be

  ; Radial wind 
    opts_vrad2                             = opts_mbgt
    delete(opts_vrad2@cnLevels)
    delete(opts_vrad2@cnFillColors)
    opts_vrad2@cnFillPalette	           = "BlueDarkRed18"
    opts_vrad2@cnLevels			   = (/-10., -5., -3., -2., -1., -0.5, -0.2, \
                                                0, 0.2, 0.5, 1., 2., 3., 5., 10./)
    opts_vrad2@cnFillColors		   = (/0,2,3,4,5,6,7,-1,-1,10,11,12,13,14,15,17/)
    opts_vrad2@lbLabelBarOn                = True
    opts_vrad2@lbBoxEndCapStyle            = "TriangleBothEnds"
    opts_vrad2@lbOrientation		   = "Vertical"
    opts_vrad2@vpXF                        = 0.25
    opts_vrad2@vpYF                        = 0.80
    opts_vrad2@vpWidthF                    = 0.55
    opts_vrad2@vpHeightF                   = 0.55

;=================
; Plot the data
;=================

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
    output = "$sam/nepartak/images/cyl_coords/mbgt_rad_xz_comp_"+typ
    wks = gsn_open_wks(opt,output)

    vr_test0 = vr_test(lev|:,rad|:)

  ; Radial wind (filled contours)
    panel_vrad  = gsn_csm_contour(wks,vr_test0(:,:),opts_vrad2) ; [u]

  ; Finally, produce the plot
    pres        = True

  ; Call 'draw' and 'frame'
    maximize_output(wks,pres)

end