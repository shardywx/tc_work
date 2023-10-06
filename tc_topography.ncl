; Plot terrain height (topography) for RA1T model domain
; Ocean fill code from --> https://www.ncl.ucar.edu/Applications/Scripts/topo_3.ncl

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/rad_circ.ncl"

; Run using:

; ncl opt=\"x11\" tc_topography.ncl

; 'opt'   = output file format ("pdf" or "x11")

begin

;==============================
; Prelash (reading in files)
;==============================

; Define file path
  fili  = "/nfs/a319/earshar/ra1t_topography.nc"

; Read in data
  f     = addfile(fili,"r")

; Only interested in terrain height (topography) --> create 2D array 
  lat   = f->latitude
  lon   = f->longitude
  topo  = f->ht(0,0,:,:)

; Set all values below -100 to missing (remove all ocean elevation values)
  topo  = where( topo .lt. -100., topo@_FillValue, topo)

; Domain boundaries (for final plot)
  tsize = dimsizes(lat)
  nsize = dimsizes(lon)

  t0    = lat(0)
  t1    = lat(tsize-1)
  n0    = lon(0)
  n1    = lon(nsize-1)

;============================================ 
; Define output file path and plot details
;============================================

; Output file location and type
  output = "$sam/nepartak/images/basic/ra1t_topography"
  wks = gsn_open_wks(opt,output)

; Choose colourmap
  gsn_define_colormap(wks,"MPL_terrain")

;==========================
; Options for plotting
;==========================

; Relative vorticity
    opts_topo                             = True
    opts_topo@cnFillOn                    = True
    opts_topo@cnLineLabelInterval         = 2.0
    opts_topo@cnLineLabelFontHeightF      = 0.012
    opts_topo@cnLineLabelBackgroundColor  = "transparent"
    opts_topo@cnLineLabelPlacementMode    = "constant"
    opts_topo@cnLinesOn                   = False
    opts_topo@cnInfoLabelOn               = False
    opts_topo@cnLevelSelectionMode        = "ExplicitLevels"

    opts_topo@cnLevels                     = (/0., 500., 1000., 1500., 2000., 2500., 3000/)
    opts_topo@cnFillColors                 = (/0, 46, 58, 70, 82, 94, 106, 118/)

    opts_topo@gsnPaperOrientation          = "landscape"
    opts_topo@tiMainString                 = ""
    opts_topo@tiMainFontHeightF            = 0.0125
    opts_topo@gsnLeftString                = ""
    opts_topo@gsnRightString               = ""
    opts_topo@gsnMaximize                  = True
    opts_topo@gsnDraw                      = False
    opts_topo@gsnFrame                     = False

  ; Additional plotting resources 
    opts_topo@mpDataBaseVersion            = "Ncarg4_1" ; More recent database
    opts_topo@mpDataSetName                = "Earth..4" ; High resolution
    opts_topo@mpGeophysicalLineColor       = "black"    ; Colour borders black
    opts_topo@mpGeophysicalLineThicknessF  = 1.0        ; Border line thickness
    opts_topo@pmTickMarkDisplayMode        = "Always"   ; Turn on map tickmarks
    opts_topo@tmXBMajorLengthF             = 0.005      ; Change tickmark length
    opts_topo@gsnMaximize                  = True       ; Maximise plot size
    opts_topo@gsnAddCyclic                 = False

  ; Dashed gridlines, every 4º
    opts_topo@mpGridAndLimbOn              = True
    opts_topo@mpGridLatSpacingF            = 5.0
    opts_topo@mpGridLonSpacingF            = 5.0
    opts_topo@mpGridLineDashPattern        = 2

  ; EDIT -- change tickmark spacing to 4º, rather than 5º

  ; Turn off tickmarks on selected axes
    opts_topo@tmXTOn                       = "False"
    opts_topo@tmYROn                       = "False"

  ; Colour fill the ocean, and do this after we have plotted the topography 
    opts_topo@mpFillOn                     = True
    opts_topo@mpOceanFillColor             = "PaleTurquoise3"
    opts_topo@mpLandFillColor              = "transparent"
    opts_topo@mpFillDrawOrder              = "PostDraw"

    opts_topo@mpLimitMode                  = "Corners"
    opts_topo@mpLeftCornerLatF             = t0
    opts_topo@mpLeftCornerLonF             = n0
    opts_topo@mpRightCornerLatF            = t1
    opts_topo@mpRightCornerLonF            = n1

  ; Do not draw the plot or advance the frame 
    opts_topo@gsnDraw                      = False
    opts_topo@gsnFrame                     = False

  ; Add labelbar
    opts_topo@lbLabelBarOn                 = True 
    opts_topo@lbBoxEndCapStyle             = "TriangleHighEnd"
    opts_topo@lbLabelFontHeightF           = 0.0115
    opts_topo@lbLabelFont                  = "Helvetica"
    opts_topo@lbPerimOn                    = False
    opts_topo@lbTitleOffsetF               = -2.00
    opts_topo@pmLabelBarOrthogonalPosF     = -0.20
    opts_topo@pmLabelBarWidthF             = 0.6
    opts_topo@pmLabelBarHeightF            = 0.15
    opts_topo@lbLabelFontHeightF           = 0.01
    opts_topo@lbTitleOn                    = True
    opts_topo@lbTitleString                = "Orography (m)"
    opts_topo@lbTitleFontHeightF           = 0.015
    opts_topo@lbTitleDirection             = "Across"

  ; Plot terrain height 
    plot_topo   = gsn_csm_contour_map(wks, topo, opts_topo)
    draw(plot_topo)
    frame(wks)

end