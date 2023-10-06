; Script to quickly plot labelbar without plotting contours 

; Run using ncl 'opt="x11"' clr=1 tc_colourbar.ncl

  ; Open workstation
    output = "lb_test_"+clr
    wks    = gsn_open_wks(opt,output)

  ; Plotting options 
    opts_mbgt                              = True
    opts_mbgt@cnFillOn                     = True
    opts_mbgt@cnLineLabelInterval          = 2.0
    opts_mbgt@cnLineLabelFontHeightF       = 0.012
    opts_mbgt@cnLineLabelBackgroundColor   = "transparent"
    opts_mbgt@cnLineLabelPlacementMode     = "constant"
    opts_mbgt@cnLinesOn                    = False
    opts_mbgt@cnInfoLabelOn                = False
    opts_mbgt@cnLevelSelectionMode         = "ExplicitLevels"

;    if (clr .eq. 1) then
;     opts_mbgt@cnFillPalette               = "ncl_default"
;     opts_mbgt@cnFillColors                = (/2,8,16,28,40,52,64,76,88,100,118,-1,-1,\
;                                               134,150,166,178,186,194,\
;                                               202,214,226,238,250/)
;     opts_mbgt@cnLevels                    = (/-150., -100., -50., -20., -10., -5., \
;                                                -3., -2., -1., -0.5, -0.2, \
;                                                0, 0.2, 0.5, 1., 2., 3., \
;                                                5., 10., 20., 50., 100., 150./)

;    elseif (clr .eq. 2) then
     opts_mbgt@cnFillPalette               = "bdr_extra"
     opts_mbgt@cnFillColors                = (/0,1,2,3,4,5,6,7,8,9,10,-1,-1,\
                                               13,14,15,16,17,18,19,20,21,22,23/)
     opts_mbgt@cnLevels                    = (/-150., -100., -50., -20., -10., -5., \
                                                -3., -2., -1., -0.5, -0.2, \
                                                0, 0.2, 0.5, 1., 2., 3., \
                                                5., 10., 20., 50., 100., 150./)

;    elseif (clr .eq. 3) then
;     opts_mbgt@cnFillPalette               = "amwg256"
;     opts_mbgt@cnFillColors                = (/0,6,12,18,27,36,45,54,66,78,90,-1,-1,\
;                                               148,156,166,178,186,194,\
;                                               202,214,226,238,250/)
;     opts_mbgt@cnLevels                    = (/-150., -100., -50., -20., -10., -5., \
;                                                -3., -2., -1., -0.5, -0.2, \
;                                                0, 0.2, 0.5, 1., 2., 3., \
;                                                5., 10., 20., 50., 100., 150./)

;    elseif (clr .eq. 4) then
;     opts_mbgt@cnFillPalette               = "GMT_polar"
;     opts_mbgt@cnFillColors                = (/3,5,7,-1,11,14,18/)
;     opts_mbgt@cnLevels                    = (/-20,-5,0,0.5,5,20/)
;    end if

    opts_mbgt@lbLabelBarOn		   = True
    opts_mbgt@lbAutoManage		   = False
    opts_mbgt@lbFillColors		   = opts_mbgt@cnFillColors
    opts_mbgt@lbPerimOn		           = False
    opts_mbgt@lbMonoFillPattern		   = True
    opts_mbgt@lbBoxEndCapStyle		   = "TriangleBothEnds"
    opts_mbgt@lbOrientation		   = "Horizontal"
    opts_mbgt@vpWidthF			   = 0.6
    opts_mbgt@vpHeightF			   = 0.15	
    opts_mbgt@lbTitleString		   = "Radial momentum budget (m s~S~-1~N~ h~S~-1~N~)"
    opts_mbgt@lbTitleFontHeightF	   = 0.020
    opts_mbgt@lbTitlePosition		   = "Bottom"
    opts_mbgt@lbTitleOffsetF		   = 0.1
    opts_mbgt@gsnDraw			   = False
    opts_mbgt@gsnFrame			   = False

    nboxes				   = dimsizes(opts_mbgt@cnLevels)
;    labels				   = new(dimsizes(opts_mbgt@cnLevels),"string")

    print("Colour table = "+opts_mbgt@cnFillPalette)

;    if (clr .eq. 4) then 
;     labels = (/"-20","-5","0","0.5","5","20"/)    
;    else
     labels = (/"-150.0","-100.0","-50.0","-20.0","-10.0","-5.0",\
     	        "-3.0","-2.0","-1.0","-0.5","-0.2",\
     	        "0.0","0.2","0.5","1.0","2.0","3.0","5.0","10.0",\
		"20.0","50.0","100.0","150.0"/)
;    end if 
    
    gsn_labelbar_ndc(wks,nboxes,labels,0.2,0.8,opts_mbgt)
    frame(wks)