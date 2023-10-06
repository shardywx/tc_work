; Script to calculate circulation following a TC, and produce a time-height plot 

; To run script, use command below, choosing appropriate values for each option
; ncl 'opt="pdf"' 'ens="em07"' rad=0.48 dist=1.8 nbox=25 tc_cp_ens_xz_cbgt_core.ncl

; 'opt' = output file format (x11 or PDF); 'ens0/1' = ensemble simulation (00-11)
; 'rad' = radius of box around TC for budget calculations (integer multiple of dx,dy)
; 'dist' = extent of gridded data either side of cyclone (subset - speeds up script)
; 'nbox' = number of ensemble box positions for circ. budget calculations (49,25,9,1)

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"

begin

;==================================================================================
;;;;;;;;;;;;;;;;;;;;;;;;;;; Start of main part of script ;;;;;;;;;;;;;;;;;;;;;;;;;;
;==================================================================================
   
    print_clock("Starting main part of script...")

  ; List all files to read in and analyse
  ; Simulation runs from 00Z 4 to 9 July [120 h]
    diri        = "$ens/reg_netcdf/"
    fili_prefix = "20160704T0000Z_ra1t_em07_pb"
    fili_a      = systemfunc("cd "+diri+" ; ls "+fili_prefix+"*.nc")
    fili        = diri+"/"+fili_a
    numINPUT    = dimsizes(fili)            ; Number of input files [21]
    numTIMES    = 120                       ; Total number of times in all files [120]
    ct	        = 0		            ; Counter variable

    setvalues NhlGetWorkspaceObjectId
      "wsMaximumSize" : 1000000000
    end setvalues

;===============================
; Start loop over input files 
;===============================

  f      = addfile(fili,"r")

;==================================================================
; Get the variables we will need (at each time interval in turn)
;==================================================================

   it    = 0

   time  = f->t(it)        ; Time in file
   times = dimsizes(time)

   lon1  = f->longitude_1  ; longitude  [160 grid points]
   lat1  = f->latitude_1   ; latitude	[160 grid points]
   lon   = f->longitude    ; longitude	[160 grid points]
   lat   = f->latitude     ; latitude	[160 grid points]

   ltsize = dimsizes(lat)
   lnsize = dimsizes(lon)

   pres  = f->p		      ; Pressure (levels)
   plevs = dimsizes(pres)     ; Size of pressure level array
   pres@units = "hPa"

   u     = f->u(it,:,:,:)     ; Zonal velocity (m/s)
   v     = f->v(it,:,:,:)     ; Meridional velocity (m/s)
   z     = f->ht(it,:,:,:)    ; Geopotential height (m)

 ; Convert to units we want
   z     = z/10
   z@units = "dam" ; Convert to decametres

   dx    = lon(1) - lon(0) ; Grid spacing (longitude)
   dy    = lat(1) - lat(0) ; Grid spacing (latitude)

;=====================================================
; Create correct date strings for each output time
;=====================================================

   month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
                     "Jul","Aug","Sep","Oct","Nov","Dec"/)

   times = dimsizes(time) ; Files are not all same size
   utc_date = cd_calendar(time, 0)

   year   = tointeger(utc_date(:,0))
   month  = tointeger(utc_date(:,1))
   day    = tointeger(utc_date(:,2))
   hour   = tointeger(utc_date(:,3))
   minute = tointeger(utc_date(:,4))
   second = utc_date(:,5)

;====================================================
; Open workstation for plotting (time-height plots)
;====================================================

     output = "$sam/nepartak/images/cart_coords/ra1t_domain"
     wks = gsn_open_wks(opt,output)

;========================
; Options for plotting 
;========================

   ; Circulation
     opts_terr                             = True
     opts_terr@mpOutlineOn		   = True
     opts_terr@mpFillOn			   = True
     opts_terr@mpFillColors		   = (/"background","deepskyblue","green3",\
     					       "blue1"/)
     opts_terr@mpGridAndLimbOn		   = True
     opts_terr@mpLimitMode		   = "Corners"
     opts_terr@mpLeftCornerLatF		   = lat(0)
     opts_terr@mpLeftCornerLonF		   = lon(0)
     opts_terr@mpRightCornerLatF	   = lat(ltsize-1)
     opts_terr@mpRightCornerLonF	   = lon(lnsize-1)

     print("Grid extends from "+lat(0)+ " to "+lat(ltsize-1)+" deg N, and from "+lon(0)+\
     	   " to "+lon(lnsize-1)+" deg E")

     opts_terr@mpDataBaseVersion           = "Ncarg4_1"     ; More recent database
     opts_terr@mpDataSetName               = "Earth..4"     ; High resolution
     opts_terr@mpOutlineBoundarySets       = "National"     ; National borders

     opts_terr@mpGridSpacingF              = 5.0            ; Spacing between grid lines
     opts_terr@mpGridLineDashPattern	   = 2		    ; Dash pattern
     opts_terr@pmTickMarkDisplayMode       = "Always"       ; Turn on map tickmarks
     opts_terr@tmXBMajorLengthF            = 0.005          ; Change tickmark length
     opts_terr@tmXBMinorOn                 = False          ; Minor tm (x-axis)
     opts_terr@tmYLMinorOn                 = False          ; Minor tm (y-axis)
     opts_terr@tmXTOn                      = "False"        ; Tm (top x-axis)
     opts_terr@tmYROn                      = "False"        ; Tm (right y-axis)

;=================
; Plot the data
;=================

     plot_terr = gsn_csm_map(wks,opts_terr)
     draw(plot_terr)
     frame(wks)

end