; Script to calculate the tangential momentum budget as in Huang et al. (2018)

; Run using:

; ncl opt=\"x11\" w0=0.3 z0=0 dist=1.2 clr=5 calc=1 rad0=0 nt=73
; nr=22 ar=1.05 lgd0=\"tr\" new0=0 cn0=\"geo_sm\" sm=0 sc=5 w_check=1
; nts=36 ntf=54 mean0=1 grp=\"7\" run=1 plt=0 typ=\"test\" tc_mbgt_tan_cxz_rhs.ncl

; 'opt'   = output file format ("pdf" or "x11")
; 'dist'  = size of box following storm (degrees)
; 'clr'	  = colour map for plots (option 4 is used in Roger Smith's papers):
; 	    "ncl_default" [1],  "bdr_extra" [2], "amwg256" [3], "GMT_polar" [4]
; 'calc'  = centre calculated on each level (0), calculated offline (1) or using 980 m vort (2)
; 'rad0'  = radius of final plots: 150 km (0) or 200 km (1)
; 'typ'	  = VC phase ("sym","asym","stoa","atos")
; 'grp'   = which set of times to analyse: "6", "5", "4c", "4b", "4a", "4"
; 'ar'    = distance in degrees of outer radius (1.0, 2.0, etc)
; 'nr'    = number of radial circles between r = 0 and r = 'ar' (21, 41, 61, etc)
; 'nt'    = number of azimuth angles in cylindrical grid
; 'new0'  = reduced scale for plotting tangential wind and AAM (1)
; 'mean0' = plot eddy and mean terms (1) instead of diffusive tendencies (0)
; 'cn0'   = offline storm centre: "slp", "vort", "geo", "geo_sm"
; 'nts'   = start time for no-VC composite 
; 'ntf'   = finish time for no-VC composite

; 500 m [10], 1 km [14], 1.5 km [17], 2 km [20], 3 km [24], 4 km [28], 5 km [31], 6 km [34]  
; 7 km [37], 8 km [39], 9 km [42], 10 km [45], 11 km [46], 12 km [48], 13 km [50], 14 km [52]
; 15 km [54], 16 km [56]

; 'pc' stream --> 'it'
; 'centre'    --> 'it-1'
; 'pd' stream --> 'it-2'

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$sam/ncl_scripts/nepartak/archer_march2018/gsn_csm.ncl"
load "$sam/ncl_func/st_rm.ncl"
load "$sam/ncl_func/nc_times.ncl"
load "$sam/ncl_func/st_centre.ncl"
load "$sam/ncl_func/setup_cyl.ncl"
load "$sam/ncl_func/ring_mono.ncl"

begin

;==============================
; Prelash (reading in files)
;==============================

; Size of arrays to read in, if we're looking at tails of v/slp distributions
  if (typ .eq. "v") then
   if (thr .eq. 2) then
    arr_size = 121
   elseif (thr .eq. 2.5) then
    arr_size = 73
   elseif (thr .eq. 3) then
    arr_size = 41
   elseif (thr .eq. 3.5) then
    arr_size = 27
   elseif (thr .eq. 4) then
    arr_size = 17
   elseif (thr .eq. -1) then 
    arr_size = 92
   elseif (thr .eq. -1.5) then 
    arr_size = 52
   elseif (thr .eq. -2) then 
    arr_size = 26
   elseif (thr .eq. -2.5) then 
    arr_size = 12
   end if
  elseif (typ .eq. "slp") then
   if (thr .eq. -1) then
    arr_size = 92
   elseif (thr .eq. -1.5) then
    arr_size = 52
   elseif (thr .eq. -2) then
    arr_size = 26
   elseif (thr .eq. -2.5) then
    arr_size = 12
   elseif (thr .eq. 2) then 
    arr_size = 121
   elseif (thr .eq. 2.5) then 
    arr_size = 73
   elseif (thr .eq. 3) then 
    arr_size = 41
   elseif (thr .eq. 3.5) then 
    arr_size = 27
   elseif (thr .eq. 4) then 
    arr_size = 17
   end if
  end if

; Retrieve file size information using 'systemfunc'                                    
  if (typ .ne. "erc1" .and. typ .ne. "erc2" .and. typ .ne. "erc3" .and. typ .ne. "test") then 
   fili  = "$ar/text/vc_"+typ+"_dat_group"+grp+".txt"
   size0 = systemfunc("wc -l < "+fili)
   size  = toint(size0)
  end if 

; Read in data from text files according to VC phase 
  if (typ .eq. "sym") then 
   typ0   = "ring"
  elseif (typ .eq. "s2a") then 
   typ0   = "r2m"
  elseif (typ .eq. "asym") then 
   typ0   = "mono"
  elseif (typ .eq. "a2s") then 
   typ0   = "m2r"
  elseif (typ .eq. "novc") then 
   typ0 = "novc"
  elseif (typ .eq. "erc1") then 
   typ0 = "erc1"
  elseif (typ .eq. "erc2") then 
   typ0 = "erc2"
  elseif (typ .eq. "erc3") then 
   typ0 = "erc3"
  elseif (typ .eq. "test") then 
   typ0 = "test"
  end if 

; Which set of simulations do we want to analyse? 

; No VCs
  if (typ .eq. "novc") then 
   novc_size = 16
   dat  = asciiread("$ar/text/novc_dat.txt",(/novc_size/),"string")
   ens0 = asciiread("$ar/text/novc_sim.txt",(/novc_size/),"string")
   ts0  = new(novc_size,integer)
   ts0  = nts0
   tf0  = new(novc_size,integer)
   tf0  = ntf0
; Wind speed or MSLP tendency
  elseif (typ .eq. "v" .or. typ .eq. "slp") then 
   dat  = asciiread("$ar/text/vc_scatter_sc"+sc+"_"+typ+"_"+thr+"_dat.txt",(/arr_size/),"string")
   ens0 = asciiread("$ar/text/vc_scatter_sc"+sc+"_"+typ+"_"+thr+"_ens.txt",(/arr_size/),"string")
   ts0  = asciiread("$ar/text/vc_scatter_sc"+sc+"_"+typ+"_"+thr+"_time.txt",(/arr_size/),"integer")
   tf0  = ts0 
; ERC intensifying phase 
  elseif (typ .eq. "erc1") then 
   dat  = (/"03T00"/)
   ens0 = (/"em04"/)
   ts0  = (/78/)
   tf0  = (/108/)
; ERC weakening phase 
  elseif (typ .eq. "erc2") then 
   dat  = (/"03T00"/)
   ens0 = (/"em11"/)
   ts0  = (/72/)
   tf0  = (/102/)
; Test mode
  elseif (typ .eq. "test") then 
   dat  = (/"03T00"/)
   ens0 = (/"em04"/)
   ts0  = (/92/)
   tf0  = (/93/)
; Any other VC phase 
  else
   dat  = asciiread("$ar/text/vc_"+typ+"_dat_group"+grp+".txt",size,"string")
   ens0 = asciiread("$ar/text/vc_"+typ+"_sim_group"+grp+".txt",size,"string")
   ts0  = asciiread("$ar/text/vc_"+typ+"_ts_group"+grp+".txt",size,"integer")
   tf0  = asciiread("$ar/text/vc_"+typ+"_tf_group"+grp+".txt",size,"integer")
  end if 

;========================================================================
; Before main loops, calculate max number of times in any one VC phase 
;======================================================================== 

; Ending time indices of all VC phases 
  tx   = tf0 + 1

; Starting time indices of all VC phases 
  tn   = ts0 - 1

; Starting minus ending time indices of all VC phases
  td   = (tx - tn) + 1

; Maximum number of times in any one VC phase
  nts  = max(td)

; Number of individual simulations to analyse for composite 
  nsim = dimsizes(dat)

; Total number of elements in time-dimension
  ntot = nts * nsim 

;======================================================================================
; Create array to hold all values for composite (previous used with 'wrt=1' option)
;======================================================================================

; Dimensions of 5D array (levs; times; radial circles; diagnostics; sims)
  nlev      = 59
  ndiag     = 10

  plot_arr0 = new( (/nlev, nts, nr, ndiag, nsim/), "float")

;==========================================================
; Now start loop over these chosen input files and times
;==========================================================

 do st = 0, nsim-1 

  dt0 = tf0(st) - ts0(st)

; Find 'pc' and 'pd' files (model height level data)
  diri        = "/nfs/a319/earshar/"+dat(st)+"/"+ens0(st)
  fili_p      = "201607"+dat(st)+"00Z_NPTK_4p4_L80_ra1t_"+ens0(st)
  fili_c1     = systemfunc("cd "+diri+" ; ls "+fili_p+"_pc*.nc")
  fili_c      = diri+"/"+fili_c1
  fili_d1     =	systemfunc("cd "+diri+" ; ls "+fili_p+"_pd*.nc")
  fili_d      = diri+"/"+fili_d1

; Calculate number of times in each file using built-in function ('nc_times')
  numINFO_c   = nc_times(fili_c)
  numINFO_d   = nc_times(fili_d)

; Number/array of times from each file
  numINPUT_c  = numINFO_c[0]
  time_c      = numINFO_c[1]
  numINPUT_d  = numINFO_d[0]
  time_d      = numINFO_d[1]

;============================================================
; Calculate storm motion using built-in function ('st_rm')
;============================================================

; 'dat'       = initialisation time (02T12, 03T00, ...)
; 'ens0'      = ensemble simulation (em00, em01, ...)
; 'diri'      = path to input files (see above)
; 'dist'      = size of box following storm (degrees)
; 'mins'      = analyse 1-h (0) or 5-min (1) data

; Zonal; meridional; vector wind; times; lat; lon; centre
  storm_rel = st_rm(dat(st), ens0(st), "$ar/text/", dist, 0)
; Choose whether to calculate a running average on the storm motion data, or leave alone
  if (run .eq. 1) then 
   u_cyc    = runave_Wrap(storm_rel[0],3,0)
   v_cyc    = runave_Wrap(storm_rel[1],3,0)
  else
   u_cyc    = storm_rel[0]
   v_cyc    = storm_rel[1]
  end if 
  lat_arr   = storm_rel[4]
  lon_arr   = storm_rel[5]
  centre    = storm_rel[6]

;==========================================================================
; Create date/time string arrays for all times in file (use 'pc' stream)
;==========================================================================

  month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun",\
                    "Jul","Aug","Sep","Oct","Nov","Dec"/)

  utc_date = cd_calendar(time_c, 0)

  year   = tointeger(utc_date(:,0))
  month  = tointeger(utc_date(:,1))
  day    = tointeger(utc_date(:,2))
  hour   = tointeger(utc_date(:,3))
  minute = tointeger(utc_date(:,4))
  second = utc_date(:,5)

; Correct for errors in the code (round up value of hour when minutes = 59)
  do ct0 = 0, numINPUT_c-1
   if (minute(ct0).gt.30) then
     hour(ct0) = hour(ct0)+1
   end if
  end do
  delete(ct0)

; Create arrays to hold strings for output plots
  date_str  = new(numINPUT_c,string)
  time_str  = new(numINPUT_c,string)
  time_arr  = new(numINPUT_c,string)
  title_arr = new(numINPUT_c,string)

  date_str = sprinti("%0.2iUTC ", hour) + \
             sprinti("%0.2i ", day) \
           + month_abbr(month)
  time_str = sprinti("%0.2i UTC ", hour) + \
             sprinti("%0.2i ", day) \
           + month_abbr(month)
  time_arr = sprinti("%0.2i", day) + \
             month_abbr(month) + \
             "_" + sprinti("%0.2iZ", hour)

;==================================================================================
;;;;;;;;;;;;;;;;;;;;;;;;;;; Start of main part of script ;;;;;;;;;;;;;;;;;;;;;;;;;;
;==================================================================================

  print_clock("Starting main part of script...")

;========================================
; Define cylindrical coordinate arrays
;========================================

  pii        = 3.14159265
  radii      = fspan(0.0,ar,nr)      ; Array of radial distances ('Start','End','Number')
  thetas     = fspan(0.0,360.0,nt)   ; Array of azimuth angles (0, 5, 10, etc)
  thetas_rad = (pii/180)*thetas      ; Degrees to radians
  rad_size   = dimsizes(radii)-1     ; Index for later calculation
  num_l      = 59      	   	     ; Number of model levels (edit to evolve smoothly)

; Define constants
  rd	     = 287.0	             ; Gas constant for dry air (J/kg/K^2)
  cp0	     = 1004.0                ; Specific heat of dry air at constant pressure (J/kg/K)
  kp	     = 0.286                 ; For Exner function calculations (0.286)
  pref	     = 1000.0                ; Reference pressure p0 (1000 hPa)

;========================================================================================
;========================================================================================
; Calculating the terms from Montgomery et al. (2018) - their Eq.(4)
; Ignore the perturbation PGF term because it is several orders of magnitude smaller

; LHS (1) : d{v}/dt          - Local tendency of the mean tangential wind

; RHS (1) : -{u} * {f+vort}  - Mean radial influx of absolute vertical vorticity
; RHS (2) : -{w} * (d{v}/dz) - Mean vertical advection of mean tangential momentum
; RHS (3) : -{u' * vort'}    - Eddy radial vorticity flux
; RHS (4) : -{w' * (dv'/dz)} - Vertical advection of eddy tangential momentum
; RHS (5) : {D_v}            - Combined diffusive and PBL tendency

;========================================================================================
;========================================================================================

; Arrays to hold all 'xpos' and 'ypos' position arrays
  xpos_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)
  ypos_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)

; Arrays to hold lat/lon at centre of cylindrical grid 
  xcen_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)
  ycen_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)

; Arrays for perturbation (eddy) variables (4D)
  ptb_vrad0 = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; [u']
  ptb_vtan0 = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; [v']
  ptb_vvel0 = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; [w']
  ptb_vort0 = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; [vort']
  ptb_dvdz0 = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; [dv'/dz]

;===========================================================
; Arrays to hold expanded 'v_int'/'vort_int' arrays (4D)
;===========================================================

  u_int0    = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)  ; Zonal wind
  v_int0    = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)  ; Meridional wind
  w_int0    = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)  ; Vertical velocity
  f_int0    = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)  ; Coriolis parameter
  kmh_int0  = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)  ; Horiz. diffusivity 
  vort_int0 = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)  ; Relative vorticity
  avo_int0  = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)  ; Absolute vorticity
  rho_int0  = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)  ; Density
  v_rad0    = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)  ; Radial wind
  v_tan0    = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)  ; Tangential wind

  difv_x0   = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)  ; Vert. diffusion (x)
  difv_y0   = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)  ; Vert. diffusion (y)
  difv_tot0 = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)  ; Vert. diff. (tot)

  difh_x0   = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)  ; Horiz. diffusion (x)
  difh_y0   = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)  ; Horiz. diffusion (y)

  u_int0@description = "Zonal velocity"
  u_int0@units       = "m s~S~-1~N~"
  u_int0!0           = "lev"
  u_int0!1           = "rad"
  u_int0!2           = "azi"
  u_int0!3           = "time"

  v_int0@description = "Meridional velocity"
  v_int0@units       = "m s~S~-1~N~"
  v_int0!0           = "lev"
  v_int0!1           = "rad"
  v_int0!2           = "azi"
  v_int0!3           = "time"

  w_int0@description = "Vertical velocity"
  w_int0@units       = "m s~S~-1~N~"
  w_int0!0           = "lev"
  w_int0!1           = "rad"
  w_int0!2           = "azi"
  w_int0!3           = "time"

  vort_int0@description = "Relative vorticity"
  vort_int0@units    = "s~S~-1~N~"
  vort_int0!0           = "lev"
  vort_int0!1           = "rad"
  vort_int0!2           = "azi"
  vort_int0!3           = "time"

  avo_int0@description = "Relative vorticity"
  avo_int0@units       = "s~S~-1~N~"
  avo_int0!0           = "lev"
  avo_int0!1           = "rad"
  avo_int0!2           = "azi"
  avo_int0!3           = "time"

  f_int0@description = "Coriolis parameter"
  f_int0@units       = "s~S~-1~N~"
  f_int0!0           = "lev"
  f_int0!1           = "rad"
  f_int0!2           = "azi"
  f_int0!3           = "time"

  kmh_int0@description = "Horizontal diffusivity"
  kmh_int0@units       = "m~S~2~N~ s~S~-1~N~"
  kmh_int0!0           = "lev"
  kmh_int0!1           = "rad"
  kmh_int0!2           = "azi"
  kmh_int0!3           = "time"

  difv_y0@description  = "Vertical diffusion (y-component)"
  difv_y0@units        = "m s~S~-2~N~"
  difv_y0!0            = "lev"
  difv_y0!1            = "rad"
  difv_y0!2            = "azi"
  difv_y0!3            = "time"

  difv_x0@description  = "Vertical diffusion (x-component)"
  difv_x0@units        = "m s~S~-2~N~"
  difv_x0!0            = "lev"
  difv_x0!1            = "rad"
  difv_x0!2            = "azi"
  difv_x0!3            = "time"

  difh_x0@description = "Horizontal diffusion (x-component)"
  difh_x0@units       = "m s~S~-2~N~"
  difh_x0!0           = "lev"
  difh_x0!1           = "rad"
  difh_x0!2           = "azi"
  difh_x0!3           = "time"

  difh_y0@description = "Horizontal diffusion (y-component)"
  difh_y0@units       = "m s~S~-2~N~"
  difh_y0!0           = "lev"
  difh_y0!1           = "rad"
  difh_y0!2           = "azi"
  difh_y0!3           = "time"

  v_rad0@description  = "Radial velocity"
  v_rad0@units        = "m s~S~-1~N~"
  v_rad0!0            = "lev"
  v_rad0!1            = "rad"
  v_rad0!2            = "azi"
  v_rad0!3            = "time"

  v_tan0@description  = "Tangential velocity"
  v_tan0@units        = "m s~S~-1~N~"
  v_tan0!0            = "lev"
  v_tan0!1            = "rad"
  v_tan0!2            = "azi"
  v_tan0!3            = "time"

;==============================================================
; Arrays for RHS source terms containing eddy variables (4D)
;==============================================================

; Eddy radial vorticity flux: RHS (3)
; {(u' * vort')}
; ptb_vez1  = (ptb_vrad0 * ptb_vort0)
  ptb_vez1  = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; [R3]

; Vertical advection of eddy tangential momentum: RHS (4)
; [-w' * dv'/dz]
; ptb_vev1  = - (ptb_vvel0 * ptb_dvdz0)
  ptb_vev1  = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; [R4]

; Arrays for RHS source terms after azimuthally averaging (3D)
  ptb_vez0  = new((/num_l,dimsizes(radii),nts/),float)                  ; [R3]
  ptb_vev0  = new((/num_l,dimsizes(radii),nts/),float)                  ; [R4]

;==================================================
; Arrays for horizontal diffusive tendency terms
;==================================================

; Horizontal diffusivity (K_mh)
  k_mh0     = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)

; Arrays for fractional expressions in diffusive tendency source terms below (4D)
; (1) du/dy, (2) d(v/r)/dr, (3) u/r, (4) r, (5) d/dr, (6) d/dy
  dudy0     = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; (1)
  dvdr0     = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; (2)
  v_r0      = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; (3)
  r_0       = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; (4)
  ddr0      = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; (5)
  ddy0      = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; (6)

; Additional arrays for terms outside main loop
  du	    = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ;
  dy	    = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ;
  dr	    = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ;
  r1	    = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ;

; Stress tensor (tau_ry)
; { K_mh * ( (1/r * du/dy) + d(v/r)/dr ) }
; tau_ry1   = k_mh0 * ( ( (1 / r_0) * dudy0) + ( v_r0 * ddr0) )
  tau_ry1   = new((/num_l,dimsizes(radii),nts/),float)                  ; [3D array]

; Total horizontal diffusive tendency of tangential momentum (using 'tau_ry')
; ( (1/(r^2) * {rho}) * d( r ^ 2 * {rho} * {tau_ry} ) /dr)
; v_dh0     = ( ( 1/( r_0 ^ 2 * mean_rho0) ) *
;             (ddr0 * ( r_0 ^ 2 * mean_rho0 * tau_ry1) ) )
  v_dh0     = new((/num_l,dimsizes(radii),nts/),float)                  ; [3D array]

; Total vertical diffusive tendency of tangential momentum
; v_dz0     = (1 /{rho0}) * ( (d/dz) * ({rho0} * {tau_yz}) )
  v_dz0     = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; [4D array]

;========================================
; Define 4D arrays for mean quantities
;========================================

; Mean tangential wind {v}
  mean_vtan0       = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)
  mean_vtan0!0     = "lev"
  mean_vtan0!1     = "rad"
  mean_vtan0!2	   = "azi"
  mean_vtan0!3     = "time"
  mean_vtan0@description = "Azimuthally averaged tangential wind"
  mean_vtan0@units = "m s~S~-1~N~"

; Mean radial wind {u}
  mean_vrad0       = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)
  mean_vrad0!0     = "lev"
  mean_vrad0!1     = "rad"	
  mean_vrad0!2     = "azi"
  mean_vrad0!3     = "time"
  mean_vrad0@description = "Azimuthally averaged radial wind"
  mean_vrad0@units = "m s~S~-1~N~"

; Mean vertical velocity {w}
  mean_vvel0       = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)
  mean_vvel0!0     = "lev"
  mean_vvel0!1     = "rad"
  mean_vvel0!2     = "azi"
  mean_vvel0!3     = "time"
  mean_vvel0@description = "Azimuthally averaged vertical velocity"
  mean_vvel0@units = "m s~S~-1~N~"

; Mean absolute angular momentum {aam}
  mean_aam0        = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)
  mean_aam0!0      = "lev"
  mean_aam0!1      = "rad"
  mean_aam0!2      = "azi"
  mean_aam0!3      = "time"
  mean_aam0@description = "Azimuthally averaged absolute angular momentum"
  mean_aam0@units  = "m s~S~-1~N~"

; Mean horizontal diffusivity {k_mh}
  mean_kmh0        = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)
  mean_kmh0!0      = "lev"
  mean_kmh0!1      = "rad"
  mean_kmh0!2      = "azi"
  mean_kmh0!3      = "time"
  mean_kmh0@description = "Azimuthally averaged horizontal diffusivity"
  mean_kmh0@units  = "m~S~2~N~ s~S~-1~N~"

; Mean relative vorticity {zeta}
  mean_vort0       = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)
  mean_vort0!0     = "lev"
  mean_vort0!1     = "rad"
  mean_vort0!2     = "azi"
  mean_vort0!3     = "time"
  mean_vort0@description = "Azimuthally averaged relative vorticity"
  mean_vort0@units = "s~S~-1~N~"

; Mean absolute vorticity {zeta + f0}
  mean_avo0        = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)
  mean_avo0!0      = "lev"
  mean_avo0!1      = "rad"
  mean_avo0!2      = "azi"
  mean_avo0!3      = "time"
  mean_avo0@description = "Azimuthally averaged absolute vorticity"
  mean_avo0@units  = "s~S~-1~N~"

; Mean Coriolis force {f0}
  mean_cor0        = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)
  mean_cor0!0      = "lev"
  mean_cor0!1      = "rad"
  mean_cor0!2      = "azi"
  mean_cor0!3      = "time"
  mean_cor0@description = "Azimuthally averaged Coriolis force"
  mean_cor0@units  = "s~S~-1~N~"

; Mean density {rho}
  mean_rho0        = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)
  mean_rho0!0      = "lev"
  mean_rho0!1      = "rad"
  mean_rho0!2      = "azi"
  mean_rho0!3      = "time"
  mean_rho0@description = "Azimuthally averaged density"
  mean_rho0@units  = "kg m~S~-3~N~"

; Mean vertical diffusion (total)
  mean_difv0          = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)
  mean_difv0!0        = "lev"
  mean_difv0!1        = "rad"
  mean_difv0!2        = "azi"
  mean_difv0!3        = "time"
  mean_difv0@description = "Azimuthally averaged vertical diffusion (tot)"
  mean_difv0@units    = "m s~S~-1~N~"

; Vertical derivative of mean tangential wind [ d{v} / dz ]
  mean_dvdz0       = new((/num_l,dimsizes(radii),nts/),float)
  mean_dvdz0!0     = "lev"
  mean_dvdz0!1     = "rad"
  mean_dvdz0!2     = "time"
  mean_dvdz0@description = "Azimuthally averaged 'dv_dz'"
  mean_dvdz0@units = "s~S~-1~N~"

;=========================
; Start multiple loops
;=========================

  print_clock("Working on ensemble member: "+ens0(st))

; Set size of arrays 
  llbox	      = toint(dist*50)
  dsize	      = (/llbox,llbox/)
  d0          = dsize(0)
  d1	      = dsize(1)

; Also set time-dependent counter variable (08/10/2019 EDIT)
; 'ct' refreshes to zero (0) at the end of each simulation loop (do st = ...)
  ct          = 0

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 1000000000
  end setvalues

;=======================================
; Choose which storm track to read in
;=======================================

  if (cn0 .eq. "slp") then
   cn = 0
  elseif (cn0 .eq. "vort") then
   cn = 1
  elseif (cn0 .eq. "geo") then
   cn = 2
  elseif (cn0 .eq. "geo_sm") then
   cn = 3
  end if

;==========================================================  
; Read in required variables from 'pc' and 'pd' streams
;========================================================== 

; Create array to hold storm track positions using different methods  
  centre_new  = new((/7, num_l, nts, 2/),"float")

  c  = addfile(fili_c,"r")         	 ; Read in 'pc' stream [u,v,vort]
  d  = addfile(fili_d,"r")         	 ; Read in 'pd' stream [w,p,t,theta,kmh,tau]

;==================================
; Get the variables we will need
;==================================

  do it = ts0(st), tf0(st)		 ; Loop over times in file

    print_clock("Looping between T+"+ts0(st)+" and T+"+tf0(st)+" ("+ens0(st)+", "+dat(st)+")")

    print("Working on time: "+time_str(it)+" (T+"+it+")" )
    title_arr(it) = "Valid at "+time_str(it)+" (T+"+it+")"

;=====================================================
; Read in variables from 'pc' stream (model levels)
;=====================================================

    t0    = lat_arr(cn,it-1,0)
    t1    = lat_arr(cn,it-1,1)
    n0    = lon_arr(cn,it-1,0)
    n1    = lon_arr(cn,it-1,1)
    print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

  ; Combine values above into array for input into external function
    ll_arr = (/t0,t1,n0,n1/) 

  ; Read in basic variables from 'pc' stream
    lon1C  = c->longitude_1({n0:n1})              ; longitude ['d1' grid points]
    lat1C  = c->latitude_1({t0:t1})               ; latitude  ['d0' grid points]
    lonC   = c->longitude({n0:n1})                ; longitude ['d1' grid points]
    latC   = c->latitude({t0:t1})                 ; latitude  ['d0' grid points]
    hybC   = c->hybrid_ht(:) 		          ; 63 model ('theta') levels

    u_varname = "x-wind"
    v_varname = "y-wind"

  ; Read in horizontal wind components 
  ; Both variables below -- [hybC | 63] * [latC | 150] * [lonC | 150]
    u0  = c->$u_varname$(it,:,{t0:t1},{n0:n1}) ; Zonal wind	
    v0  = c->$v_varname$(it,:,{t0:t1},{n0:n1}) ; Meridional wind 

;===================================
; Calculate storm-relative winds
;===================================

    u = u0 - u_cyc(cn,it-1)
    v = v0 - v_cyc(cn,it-1)

  ; Add metadata from original horizontal wind arrays
    copy_VarCoords(u0,u)
    copy_VarAtts(u0,u)

    copy_VarCoords(v0,v)
    copy_VarAtts(v0,v)

;==================================================================
; Calculate relative vorticity using centered finite differences
;==================================================================

  ; Option '2' --> boundary points estimated using one-sided difference scheme
    vort = uv2vr_cfd(u, v, latC, lonC, 2)
    copy_VarCoords(u, vort)
    vort@units = "s~S~-1~N~"
    vort@name = "Relative vorticity on model levels"

  ; Coriolis parameter and absolute vorticity
    f1   = coriolis_param(latC)
    f0   = conform_dims(dimsizes(vort),f1,1)
    copy_VarMeta(vort,f0)
    copy_VarCoords(vort,f0)
    f0@description = "Coriolis parameter"
    f0@name        = "Coriolis parameter"

    avo  = vort + f0
    copy_VarMeta(vort,avo)
    copy_VarAtts(vort,avo)
    copy_VarCoords(vort,avo)
    avo@description = "Absolute vorticity"
    avo@name        = "Absolute vorticity on model levels"
    avo@long_name   = "Absolute vorticity"

;===========================================================  
; Calculate rate of strain (following Nguyen et al. 2011) 
;===========================================================  

    dim1   = "latitude"
    dim2   = "longitude"

  ; Calculate horizontal derivatives of 'u' (du_dx, du_dy)
    du0    = grad_latlon_cfd(u,u&$dim1$,u&$dim2$,False,False)
    du_dy  = du0[0]
    du_dx  = du0[1]
    delete(du0)

  ; Calculate horizontal derivatives of 'v' (dv_dx, dv_dy) 
    dv0    = grad_latlon_cfd(v,v&$dim1$,v&$dim2$,False,False)
    dv_dy  = dv0[0]
    dv_dx  = dv0[1]
    delete(dv0)

  ; Calculate strain rate --> sqrt[ (du_dx - dv_dy)^2 + (dv_dx + du_dy)^2 ]  

  ; Stretching deformation 
    str1a   = du_dx - dv_dy

  ; Shearing deformation
    str2a   = dv_dx + du_dy
    str1    = (str1a)^2
    str2    = (str2a)^2

  ; Strain rate (3D array)
    str0    = sqrt( str1(:,:,:) + str2(:,:,:) )

;==========================================================
; Calculate divergence using one-sided difference scheme
;==========================================================

    div = uv2dv_cfd(u, v, latC, lonC, 2)
    copy_VarCoords(u, div)
    div@units = "s~S~-1~N~"
    div@name = "Divergence on pressure levels"

;======================================================
; Read in variables from 'pd' stream (model levels)
;======================================================

  ; Read in basic variables from 'pd' stream 
    lon1D  = d->longitude_1({n0:n1})                ; longitude [xxx grid points]
    lat1D  = d->latitude_1({t0:t1})                 ; latitude  [xxx grid points]
    lonD   = d->longitude({n0:n1})                  ; longitude [xxx grid points]
    latD   = d->latitude({t0:t1})                   ; latitude  [xxx grid points]
    hybD   = d->hybrid_ht(:)			    ; 63 model ('rho') levels

  ; Read in diffusivity and wind stress
  ; 'str_xz' -- [hybD | 63] * [latD | 150]  * [lonD | 151]
  ; 'str_yz' -- [hybD | 63] * [lat1D | 150] * [lon1D | 151]
  ; 'kmh'    -- [hybD | 63] * [latD | 150]  * [lon1D | 151]
    str_xz = d->taux(it-2,:,{t0:t1},{n0:n1})        ; Wind stress, x-comp (N m/2)
    str_yz = d->tauy(it-2,:,{t0:t1},{n0:n1})        ; Wind stress, y-comp (N m/2)
    kmh	   = d->unspecified(it-2,:,{t0:t1},{n0:n1}) ; Horizontal diffusivity (m^2 s-1)

  ; Read in temperature, pressure, theta, vertical velocity and geopotential height
  ; All variables below -- [hybD | 63] * [latD | 150] * [lon1D | 151]
    t      = d->temp(it-2,:,{t0:t1},{n0:n1})        ; Temperature (K)
    p      = d->p(it-2,:,{t0:t1},{n0:n1})           ; Pressure (Pa)
    th     = d->theta(it-2,:,{t0:t1},{n0:n1})	    ; Potential temperature (K)
    w      = d->dz_dt(it-2,:,{t0:t1},{n0:n1})	    ; Vertical velocity (m s-1)
    z      = d->ht(it-2,:,{t0:t1},{n0:n1})    	    ; Geopotential height (m)

  ; Calculate dry air density from pressure and temperature
    rho	   = p / (rd * t)

  ; Add metadata
    copy_VarCoords(p, rho)
    rho@description = "Dry air density"
    rho@units = "kg m~S~-3~N~"

  ; Also calculate specific volume (for	use in diffusive tendency source terms)
    rho_a  = 1 / rho
    copy_VarCoords(rho, rho_a)
    rho_a@description = "Specific volume"
    rho_a@units = "m~S~3~N~ kg~S~-1~N~"

  ; Convert units of pressure and geopotential height before continuing
    p  = p / 100
    z  = z / 10

  ; Add metadata
    p@units = "hPa"
    z@units = "dam"

    dx0     = lonD(1) - lonD(0) ; Grid spacing (longitude)
    dy0     = latD(1) - latD(0) ; Grid spacing (latitude)

;==========================
; Loop over model levels 
;==========================

  ; Loop over chosen model levels
    do ilev = 0, num_l-1

   ; Create output strings (model level info)
     hy  = sprintf("%0.0f",hybD(ilev) )
     hy1 = toint(hy)
     hy0 = sprintf("%05g",hy1) 

     hy2 = sprintf("%0.0f",hybD(z0) )
     hy3 = tofloat(hy2) / 1000
     print_clock("Working on model level " + ilev+ " ("+hy+" m AGL)")

   ; Calculate components of vertical diffusion term
   ; accounting for top and bottom levels
     ip1 = min((/ilev+1,num_l-1/))
     im1 = max((/ilev-1,0/))

   ; Wind stress components on current model level
     tau_xz = str_xz(ilev,0:d0-1,0:d1-1)
     tau_yz = str_yz(ilev,0:d0-1,0:d1-1)

   ; Wind stress components on model levels above/below current level
     tau_x1 = str_xz(ip1,0:d0-1,0:d1-1)
     tau_x0 = str_xz(im1,0:d0-1,0:d1-1)
     tau_y1 = str_yz(ip1,0:d0-1,0:d1-1)
     tau_y0 = str_yz(im1,0:d0-1,0:d1-1)
     dz	    = hybD(ip1) - hybD(im1)

   ; Vertical derivative of x-component of wind stress
     dtaux_dz = (tau_x1 - tau_x0) / dz

   ; Vertical derivative of y-component of wind stress
     dtauy_dz = (tau_y1 - tau_y0) / dz

   ; X-component of vertical diffusive tendency
     difv_xplane             = rho_a(ilev,:,:) * dtaux_dz
     copy_VarCoords(tau_xz,difv_xplane)
     copy_VarMeta(tau_xz,difv_xplane)
     difv_xplane@units       = "m s~S~-2~N~"
     difv_xplane@title       = "vertical diffusive tendency (x-comp)"
     difv_xplane@long_name   = "vertical diffusive tendency (x-comp)"
     difv_xplane@name        = "difv_x"

   ; Y-component of vertical diffusive tendency
     difv_yplane             = rho_a(ilev,:,:) * dtauy_dz
     copy_VarCoords(tau_yz,difv_yplane)
     copy_VarMeta(tau_yz,difv_yplane)
     difv_yplane@units       = "m s~S~-2~N~"
     difv_yplane@title       = "vertical diffusive tendency (y-comp)"
     difv_yplane@long_name   = "vertical diffusive tendency (y-comp)"
     difv_yplane@name        = "difv_y"

   ; Read in other variables on model levels
     u_plane 	            = u(ilev,:,:)	   ; Zonal wind 
     v_plane		    = v(ilev,:,:)	   ; Meridional wind 
     w_plane		    = w(ilev,:,:)	   ; Vertical velocity
     f_plane		    = f0(ilev,:,:)	   ; Coriolis parameter
     rho_plane              = rho(ilev,:,:)	   ; Density
     prs_plane              = p(ilev,:,:)	   ; Pressure
     str_plane              = str0(ilev,:,:)       ; Strain rate
     vort_plane	   	    = vort(ilev,:,:)	   ; Relative vorticity
     avo_plane		    = avo(ilev,:,:)	   ; Absolute vorticity
     kmh_plane		    = kmh(ilev,:,:)	   ; Horizontal diffusivity (K_mh)
     div_plane              = div(ilev,:,:)        ; Divergence

   ; Calculate vector windspeed and smooth 
     spd_plane              = sqrt( (u_plane ^ 2) + (v_plane ^ 2) )
     copy_VarCoords(u_plane, spd_plane)
     vort_pl_smth           = smth9_Wrap(vort_plane, 0.5, 0.5, True)

   ; Combine all arrays above into single, larger array for input into 'setup_cyl'

   ; 'pc' stream  
   ; full = 0 [u, v, spd, vort, avo, str, slp, div] 
   ; full = 1 [u, v, spd, vort, avo, f, str, div]
     pc_plane               = (/u_plane, v_plane, spd_plane, vort_plane, \
                                avo_plane, f_plane, str_plane, div_plane/)

   ; 'pd' stream 
   ; full = 0 [w_plane, prs_plane, pv_plane, flux_plane]
   ; full = 1 [w_plane, prs_plane, ..., ..., km_phlane, difv_xplane, difv_yplane, rho_plane]
     pd_plane               = (/w_plane, prs_plane, prs_plane, prs_plane, kmh_plane, \
                                difv_xplane, difv_yplane, rho_plane/)
     
;=======================================================
; Calculate storm centre position (external function)
;=======================================================

   ; Set radial distance (º) when looking for speed/vort min (in ringlike phase) 
     r0            = 0.12

   ; 08/10/2019 --> 'ct' = 53, 54, etc (counter variable doesn't start at 0...)
     print_clock("Calculating storm centre position...")
     centre_arr    = st_centre(ilev, it, cn, r0, centre, vort_plane, prs_plane, spd_plane, \
                               centre_new, dy0, dx0, "comp", ct, "pd", "off", "1h")

   ; Updated storm track information contained in 'centre_new'
     centre_new    = centre_arr[0] 

;====================================================================== 
; Make sure that the storm 'centre' is not within the eyewall updraft
;====================================================================== 

   if (w_check .eq. 1) then 

   ; If vertical velocity is above a threshold value at the storm centre, try again

   ; Retrieve coordinates of storm centre
     cen_loc0 = centre_new(sc,ilev,ct,0)
     cen_loc1 = centre_new(sc,ilev,ct,1)

   ; Coordinates of grid points either side
     cen_0n = cen_loc0 - dy0
     cen_0x = cen_loc0 + dy0
     cen_1n = cen_loc1 - dy0
     cen_1x = cen_loc1 + dy0

   ; Calculate vertical velocity and divergence at these grid points
     cen_w    = w_plane( {cen_0n:cen_0x}, {cen_1n:cen_1x} )
     cen_d    = div_plane( {cen_0n:cen_0x}, {cen_1n:cen_1x} )

   ; If average over several grid points is > threshold value, search for centre again
     ave_w    = avg(cen_w)
     ave_d    = avg(cen_d)
     ave_out  = sprintf("%0.1f",ave_w)

   ; Coordinates of grid points either side (bigger grid than above)
     delete([/cen_0n, cen_0x, cen_1n, cen_1x/])
     cen_0n = cen_loc0 - (2 * dy0)
     cen_0x = cen_loc0 + (2 * dy0)
     cen_1n = cen_loc1 - (2 * dy0)
     cen_1x = cen_loc1 + (2 * dy0)

     if (ilev .gt. z0 .and. ave_w .gt. w0) then

      print("ave. vertical velocity = "+ave_w+" ; ave. divergence = "+ave_d)

    ; Create smaller grid around previous centre  
      w_pl_new  = w_plane({cen_0n:cen_0x},{cen_1n:cen_1x})

    ; Find minimum vertical velocity on this smaller grid  
      w_min_new = min(w_pl_new)

    ; Reshape to 1D array, and find index of minimum vertical velocity 
      dims_w    = dimsizes(w_pl_new)
      w_1d      = ndtooned(w_pl_new)
      inds_w    = ind_resolve(minind(w_1d), dims_w)

    ; Retrieve lat/lon information from subset of grid
      ltN       = w_pl_new&latitude
      lnN       = w_pl_new&longitude_1

      lat_min_w = ltN(0) + (dy0 * inds_w(0,0) )
      lon_min_w = lnN(0) + (dx0 * inds_w(0,1) )
      delete([/ltN, lnN, w_1d, w_pl_new/])

    ; Print new storm centre to screen 
      print("Centre (w_min): "+lat_min_w+" degrees N, "+lon_min_w+" degrees E")

    ; Replace values in 'centre_new' before calling 'setup_cyl' below
      centre_new(sc,ilev,ct,0) = lat_min_w
      centre_new(sc,ilev,ct,1) = lon_min_w

     end if

   ; Tidy up 
     delete([/cen_w, cen_d, cen_0n, cen_0x, cen_1n, cen_1x/])
     
   end if 

;========================================================
; Switch to cylindrical coordinates (external function)
;========================================================

   ; Last option = 0 --> basic variables output 
   ; Last option = 1 --> additional MBGT variables output  
   ; 08/10/2019 --> 'nts' option changed from 'numTIMES'. Currently does not impact output
     print_clock("Translating to cylindrical grid...")
     cyl_arr  = setup_cyl(it, dist, ll_arr, lonC, latC, lonD, latD, lon1D, lat1D, \
                          pc_plane, pd_plane, nts, centre_new, \
                          sc, sm, ar, nr, nt, ilev, num_l, 1, "comp", ct, \
                          0, xpos_all, ypos_all, xcen_all, ycen_all)

     xpos_all = cyl_arr[0]
     ypos_all = cyl_arr[1]
     xcen_all = cyl_arr[2]
     ycen_all = cyl_arr[3]
     lat_max  = cyl_arr[4]
     lon_max  = cyl_arr[5]
     radii    = cyl_arr[6]
     rad_size = cyl_arr[7]
     u_int    = cyl_arr[8]
     v_int    = cyl_arr[9]
     spd_int  = cyl_arr[10]
     vort_int = cyl_arr[11]
     avo_int  = cyl_arr[12]

     f_int    = cyl_arr[13]
     vr_int   = cyl_arr[14]
     vt_int   = cyl_arr[15]
     w_int    = cyl_arr[16]
     prs_int  = cyl_arr[17]
     kmh_int  = cyl_arr[18]
     dvx_int  = cyl_arr[19]
     dvy_int  = cyl_arr[20]
     rho_int  = cyl_arr[21]
     dvt_int  = cyl_arr[22]
     str_int  = cyl_arr[27]
     div_int  = cyl_arr[29]

     xpos     = cyl_arr[23]
     ypos     = cyl_arr[24]
     xcen     = cyl_arr[25]
     ycen     = cyl_arr[26]

     print("Storm centre = "+lat_max+" ; "+lon_max) 

   if (plt .eq. 0) then 

;===================================== 
; START MID-SCRIPT PLOTTING ROUTINE 
;===================================== 

    if (ilev .gt. 35 .and. ilev .lt. 48) then

   ; Plot the data as we go  
     if (typ .eq. "v" .or. typ .eq. "slp") then 
      output_loop = "$nep/nepartak/images/cyl_coords/cxy_4panel_"+w0+"w_"+dat(st)+\
                    "_"+ens0(st)+"_"+typ+"_"+hy0+"m_sc"+sc+"_T"+it
     else
      output_loop = "$nep/nepartak/images/cyl_coords/cxy_4panel_"+w0+"w_"+dat(st)+\
                    "_"+ens0(st)+"_"+typ0+"_"+hy0+"m_sc"+sc+"_T"+it
     end if 
     wks = gsn_open_wks(opt, output_loop)

     gsn_define_colormap(wks,"prcp_new")

   ; Relative vorticity resources  
     opts_rv                               = True
     opts_rv@cnFillOn                      = True
     opts_rv@sfXArray                      = xpos_all(ilev,:,:,ct)
     opts_rv@sfYArray                      = ypos_all(ilev,:,:,ct)
     opts_rv@cnLineLabelInterval           = 2.0
     opts_rv@cnLineLabelFontHeightF        = 0.012
     opts_rv@cnLineLabelBackgroundColor    = "transparent"
     opts_rv@cnLineLabelPlacementMode      = "constant"
     opts_rv@cnLinesOn                     = False
     opts_rv@cnInfoLabelOn                 = False
     opts_rv@cnLevelSelectionMode          = "ExplicitLevels"
     opts_rv@cnLevels                      = (/2., 3., 5., 10., 15., 20., 25., 30., \
                                               40., 50., 60., 70., 80./)
     opts_rv@cnFillColors                  = (/0,3,4,5,6,7,8,9,10,11,12,13,14,15/)
     opts_rv@gsnMaximize                   = False

   ; No title or additional strings 
     opts_rv@tiMainString                  = ""
     opts_rv@gsnLeftString                 = ""
     opts_rv@gsnRightString                = ""

   ; Turn individual label bar on/off
     opts_rv@lbLabelBarOn                  = True
     opts_rv@lbBoxEndCapStyle              = "TriangleBothEnds"
     opts_rv@lbOrientation                 = "Vertical"

      opts_rv@trYMaxF                       = lat_max + ar
      opts_rv@trYMinF                       = lat_max - ar 
      opts_rv@trXMaxF                       = lon_max + ar 
      opts_rv@trXMinF                       = lon_max - ar

   ; Do not draw the plot or advance the frame 
      opts_rv@gsnDraw                       = False
      opts_rv@gsnFrame                      = False
      opts_rv@gsnAddCyclic                  = False

   ; Background plot resources 
     opts_xy                                = True

   ; No tickmarks on any axis
      opts_xy@tmYROn                        = False
      opts_xy@tmYLOn                        = True
      opts_xy@tmXTOn                        = False
      opts_xy@tmXBOn                        = True

      opts_xy@xyLineColor                   = -1
      opts_xy@gsnDraw                       = False
      opts_xy@gsnFrame                      = False

      opts_xy@trYMaxF                       = lat_max + ar
      opts_xy@trYMinF                       = lat_max - ar
      opts_xy@trXMaxF                       = lon_max + ar
      opts_xy@trXMinF                       = lon_max - ar
      opts_xy@gsnMaximize                   = False

   ; Resources for windspeed 
      opts_wind                             = opts_rv
      delete([/opts_wind@cnLevels, opts_wind@cnFillColors/])
      opts_wind@cnFillOn                    = True
      opts_wind@cnLinesOn                   = False
      opts_wind@cnLineLabelsOn              = False
      opts_wind@cnInfoLabelOn               = False
      opts_wind@cnLevelSelectionMode        = "ExplicitLevels"
      opts_wind@cnFillPalette               = "radar_new"
      opts_wind@cnLevels                    = (/5.0, 10.0, 15.0, 20.0,\
                                                25.0, 30.0, 35.0, 40.0, \
                                                45.0, 50.0, 55.0, 60.0, 65.0/)
      opts_wind@cnFillColors                = (/-1,0,1,2,3,4,5,6,7,8,9,10,11,12/)
      opts_wind@gsnMaximize                 = False

   ; Resources for vertical velocity       
      opts_vvel                             = opts_wind
      delete([/opts_vvel@cnLevels, opts_vvel@cnFillColors, opts_vvel@cnFillPalette/])

     opts_vvel@cnFillPalette               = "BlueDarkRed18"
     opts_vvel@cnLevels                    = (/-3.0, -2.0, -1.5, -1.0, -0.5, -0.2, -0.1, -0.05,\
                                                0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0/)
     opts_vvel@cnFillColors                = (/1,2,3,4,5,6,7,8,-1,9,10,11,12,13,14,15,16/)

   ; Resources for strain rate                                                                    
     opts_strn                             = opts_vvel
     delete([/opts_strn@cnLevels, opts_strn@cnFillColors, opts_strn@cnFillPalette/])
     opts_strn@cnFillPalette               = "circ_budget"
     opts_strn@cnLevels                    = (/10,12,14,16,18,20,22,25,30/)
     opts_strn@cnFillColors                = (/-1,11,13,15,16,17,18,19,20,21/)

   ; Resources for divergence                                                                     
     opts_divg                             = opts_vvel
     delete([/opts_divg@cnLevels, opts_divg@cnFillColors, opts_divg@cnFillPalette/])
     opts_divg@cnFillPalette               = "BlueDarkRed18"
     opts_divg@cnLevels                    = (/-100, -50, -30, -20, -10, -5, \
                                               5, 10, 20, 30, 50, 100/)
     opts_divg@cnFillColors                = (/1,2,3,4,6,8,-1,9,11,13,14,15,16/)

   ; Resources for radial wind                                                                    
     opts_vrad                             = opts_vvel
     delete([/opts_vrad@cnLevels, opts_vrad@cnFillColors, opts_vrad@cnFillPalette/])
     opts_vrad@cnFillPalette               = "BlueDarkRed18"
     opts_vrad@cnLevels                    = (/-20., -15., -10., -5., -3., -2., -1., -0.5, \
                                                0, 0.5, 1., 2., 3., 5., 10., 15., 20./)
     opts_vrad@cnFillColors                = (/0,2,3,4,5,6,7,8,-1,-1,10,11,12,13,14,15,16,17/)

   ; Panel plot label resources 
      txid                                   = new(4,graphic)
      amid                                   = new(4,graphic)

      panel_strings                          = (/ "a) Relative vorticity",\
                                                  "b) Vertical velocity",\
                                                  "c) Tangential wind",\
                                                  "d) Radial wind"/)

      txres                                  = True
      txres@txPerimOn                        = True
      txres@txFontHeightF                    = 0.0125
      txres@txBackgroundFillColor            = "White"

      amres                                  = True
      amres@amParallelPosF                   = -0.48
      amres@amOrthogonalPosF                 = -0.48
      amres@amJust                           = "TopLeft"

   ; Coordinates of panels 
      xf = (/0.05, 0.57, 0.05, 0.57/)
      yf = (/0.90, 0.90, 0.52, 0.52/)
      wf = 0.35
      hf = 0.32

   ; Panel 1 (relative vorticity)
      opts_rv@vpXF       = xf(0)
      opts_rv@vpYF       = yf(0)
      opts_rv@vpWidthF   = wf
      opts_rv@vpHeightF  = hf

      opts_xy0           = opts_xy
      opts_xy0@vpXF      = xf(0)
      opts_xy0@vpYF      = yf(0)
      opts_xy0@vpWidthF  = wf
      opts_xy0@vpHeightF = hf

      vort_int_plot      = vort_int * (10 ^ 4)

      plot_tl            = gsn_csm_xy(wks, xpos_all(ilev,:,:,ct), \
                                      ypos_all(ilev,:,:,ct), opts_xy0)
      plot_tl0           = gsn_csm_contour(wks, vort_int_plot(:,:), opts_rv)
      overlay(plot_tl, plot_tl0)

   ; Overlay labels 
      txid(0)      = gsn_create_text(wks, panel_strings(0), txres)
      amid(0)      = gsn_add_annotation(plot_tl, txid(0), amres)

   ; Panel 2 (vertical velocity)
      opts_vvel@vpXF       = xf(1)
      opts_vvel@vpYF       = yf(1)
      opts_vvel@vpWidthF   = wf
      opts_vvel@vpHeightF  = hf

      opts_xy1             = opts_xy
      opts_xy1@vpXF        = xf(1)
      opts_xy1@vpYF        = yf(1)
      opts_xy1@vpWidthF    = wf
      opts_xy1@vpHeightF   = hf

      plot_tr              = gsn_csm_xy(wks, xpos_all(ilev,:,:,ct), \
                                        ypos_all(ilev,:,:,ct), opts_xy1)
      plot_tr0             = gsn_csm_contour(wks, w_int(:,:), opts_vvel)
      overlay(plot_tr, plot_tr0)

   ; Overlay labels
      txid(1)      = gsn_create_text(wks, panel_strings(1), txres)
      amid(1)      = gsn_add_annotation(plot_tr, txid(1), amres)

   ; Panel 3 (tangential wind) 
      opts_wind@vpXF       = xf(2)
      opts_wind@vpYF       = yf(2)
      opts_wind@vpWidthF   = wf
      opts_wind@vpHeightF  = hf

      opts_xy2             = opts_xy
      opts_xy2@vpXF        = xf(2)
      opts_xy2@vpYF        = yf(2)
      opts_xy2@vpWidthF    = wf
      opts_xy2@vpHeightF   = hf

      plot_bl              = gsn_csm_xy(wks, xpos_all(ilev,:,:,ct), \
                                        ypos_all(ilev,:,:,ct), opts_xy2)
      plot_bl0             = gsn_csm_contour(wks, vt_int(:,:), opts_wind)
      overlay(plot_bl, plot_bl0)

   ; Overlay labels  
      txid(2)      = gsn_create_text(wks, panel_strings(2), txres)
      amid(2)      = gsn_add_annotation(plot_bl, txid(2), amres)

   ; Panel 4 (radial wind)
     opts_vrad0           = opts_vrad
     opts_vrad0@vpXF      = xf(3)
     opts_vrad0@vpYF      = yf(3)
     opts_vrad0@vpWidthF  = wf
     opts_vrad0@vpHeightF = hf

     opts_xy3             = opts_xy
     opts_xy3@vpXF        = xf(3)
     opts_xy3@vpYF        = yf(3)
     opts_xy3@vpWidthF    = wf
     opts_xy3@vpHeightF   = hf

     plot_br              = gsn_csm_xy(wks, xpos_all(ilev,:,:,ct), \
                                       ypos_all(ilev,:,:,ct), opts_xy3)
     plot_br0             = gsn_csm_contour(wks, vr_int(:,:), opts_vrad0)
     overlay(plot_br, plot_br0)

   ; Overlay labels
     txid(3)      = gsn_create_text(wks, panel_strings(3), txres)
     amid(3)      = gsn_add_annotation(plot_br, txid(3), amres)

   ; Tidy up before final plots later in script
     delete([/opts_vvel, opts_wind, opts_rv, opts_strn, opts_divg, opts_vrad, \
              opts_xy, opts_xy0, opts_xy1, opts_xy2, opts_xy3, \
              txid, amid, panel_strings/])

;==================================================
; Overlay markers for each storm centre estimate 
;==================================================

; (0) vort, (1) pressure, (2) pressure_smth, (3), pressure_smth2, (4) vort_min,
; (5) speed_min, (6) offline 
  m_colours  = (/"royalblue1", "orange4", "orange4", "orange4", "royalblue1", \
                 "black", "red4"/)
  m_indices  = (/12, 12, 16, 8, 8, 8, 12/)
  m_sizes    = (/10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0/)
  m_thick    = (/5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0/)

  markers    = new( (/4, dimsizes(m_colours)/), graphic)

  do i = 0, dimsizes(m_colours)-1
    mres0                    = True
    mres0@gsMarkerIndex      = m_indices(i)
    mres0@gsMarkerSizeF      = m_sizes(i)
    mres0@gsMarkerColor      = m_colours(i)
    mres0@gsMarkerThicknessF = m_thick(i)
    markers(0,i)             = gsn_add_polymarker(wks, plot_tl, centre_new(i,ilev,ct,1),\
                                                  centre_new(i,ilev,ct,0), mres0)
    markers(1,i)             = gsn_add_polymarker(wks, plot_tr, centre_new(i,ilev,ct,1),\
                                                  centre_new(i,ilev,ct,0), mres0)
    markers(2,i)             = gsn_add_polymarker(wks, plot_bl, centre_new(i,ilev,ct,1),\
                                                  centre_new(i,ilev,ct,0), mres0)
    markers(3,i)             = gsn_add_polymarker(wks, plot_br, centre_new(i,ilev,ct,1),\
                                                  centre_new(i,ilev,ct,0), mres0)
  end do

;===================================
; Add a legend (only on one panel) 
;===================================

; Display legend; give us control; label font; label height
  lg_opts                            = True
  lg_opts@pmLegendDisplayMode        = "Always"
  lg_opts@lgAutoManage               = False
  lg_opts@lgLabelFont                = "Helvetica"
  lg_opts@lgLabelFontHeightF         = 0.060

; Customise labels and line colours 
  lg_opts@lgItemType                 = "Markers"
  lg_opts@lgMarkerColors             = m_colours
  lg_opts@lgMarkerIndexes            = m_indices
  lg_opts@lgMarkerSizeF              = mres0@gsMarkerSizeF
  lg_opts@lgMarkerThicknessF         = mres0@gsMarkerThicknessF
  lsize                              = dimsizes(lg_opts@lgMarkerColors)

  lg_opts@lgLabelPosition            = "Right"
  lg_opts@lgItemPlacement            = "ExplicitPlacement"
  lg_opts@lgItemPositions            = fspan(0.05, 0.95, lsize)
  lg_opts@lgLeftMarginF              = 0.01
  lg_opts@lgRightMarginF             = 0.15

  lg_opts@vpWidthF                   = 0.150
  lg_opts@vpHeightF                  = 0.090
  lg_opts@lgPerimColor               = "black"
  lg_opts@lgPerimThicknessF          = 3.0
  lg_opts@lgPerimFill                = "SolidFill"
  lg_opts@lgPerimFillColor           = "white"

  lab_arr                            = (/"vort", "pres", \
                                         "pres_sm", "pres_sm2", \
                                         "vort_min", "windspeed_min", \
                                         "offline ("+cn0+")"/)
  lgd                                = gsn_create_legend(wks,lsize,lab_arr,lg_opts)

; Add legend to plot (top left or right hand corner)  
  am_opts                            = True

  if (lgd0 .eq. "tr") then
   am_opts@amJust                    = "TopRight"
   am_opts@amParallelPosF            = 0.5
   am_opts@amOrthogonalPosF          = -0.5
  elseif (lgd0 .eq. "tl") then
   am_opts@amJust                    = "TopLeft"
   am_opts@amParallelPosF            = -0.5
   am_opts@amOrthogonalPosF          = -0.5
  elseif (lgd0 .eq. "br") then
   am_opts@amJust                    = "BottomRight"
   am_opts@amParallelPosF            = 0.5
   am_opts@amOrthogonalPosF          = 0.5
  elseif (lgd0 .eq. "bl") then
   am_opts@amJust                    = "BottomLeft"
   am_opts@amParallelPosF            = -0.5
   am_opts@amOrthogonalPosF          = 0.5
  end if

  annotate_0                         = gsn_add_annotation(plot_bl, lgd, am_opts)

;============================== 
; Add radial circles to plot
;============================== 

      dtr = 0.017453292519943

    ; Get customizations for circles
    ; Zoom into plot --> can't remember how to do this (29/10/2019)
      radius   = radii(rad_size)
      xcenter  = lon_max
      ycenter  = lat_max
      out_thck = 2
      in_thck  = 1
      spacing  = 0.05
      dpattrn  = 2
      out_clr  = 1
      in_clr   = 1

;=============================================
; Construct and attach outer circle to plot  
;=============================================

    ; Calculate arrays for outer circle 
      degrees  = ispan(0, 360, 5)
      xcos     = cos(dtr * degrees)
      xsin     = sin(dtr * degrees)

    ; Array of points (x,y) representing the outer circle
      xc       = xcenter + (radius * xcos)
      yc       = ycenter + (radius * xsin)

    ; Resources for outer circle
      lnres    = True
      lnres@gsLineThicknessF  = out_thck
      lnres@gsLineColor       = out_clr
      lnres@gsLineDashPattern = 0

    ; Attach outer circle to each panel plot
      plot_outer0 = gsn_add_polyline(wks, plot_tl, xc, yc, lnres)
      plot_outer1 = gsn_add_polyline(wks, plot_tr, xc, yc, lnres)
      plot_outer2 = gsn_add_polyline(wks, plot_bl, xc, yc, lnres)
      plot_outer3 = gsn_add_polyline(wks, plot_br, xc, yc, lnres)

      if (spacing .gt. radius) then
       print("add_radial_circles: spacing is > radius, can't draw inner circles.")
       return
      end if

;==================================================
; Now construct and attach inner circles to plot 
;==================================================

    ; Draw inner circles if desired
      if (spacing .gt. 0) then

       count = 0
       size  = toint( (radius-spacing) / spacing )
       plot_inner = new( (/4,size/), graphic)

       do r = spacing, (radius-spacing), spacing

       ; Calculate arrays for inner circles
         xc = xcenter + (r * xcos)
         yc = ycenter + (r * xsin)

       ; Resources for inner circles
         delete([/lnres@gsLineThicknessF,lnres@gsLineColor,lnres@gsLineDashPattern/])
         lnres@gsLineColor       = in_clr
         lnres@gsLineThicknessF  = in_thck
         lnres@gsLineDashPattern = dpattrn

       ; Attach inner circles to plot
         plot_inner(0,count) = gsn_add_polyline(wks, plot_tl, xc, yc, lnres)
         plot_inner(1,count) = gsn_add_polyline(wks, plot_tr, xc, yc, lnres)
         plot_inner(2,count) = gsn_add_polyline(wks, plot_bl, xc, yc, lnres)
         plot_inner(3,count) = gsn_add_polyline(wks, plot_br, xc, yc, lnres)

         count = count + 1

       end do

      end if

    ; Tidy up
      delete([/lnres, degrees, xcos, xsin, xc, yc/])

    ; Create a main title at the top of the page
      title_res                   = True
      title_str                   = "Panel plot at "+hy+" m AGL: "+ens0(st)+\
                                    ", "+dat(st)+" (T+"+it+")"
      title_res@txFontHeightF     = 0.02
      title_res@txFont            = "Helvetica"
      xpos_title                  = 0.50
      ypos_title                  = 0.94

      gsn_text_ndc(wks, title_str, xpos_title, ypos_title, title_res)

    ; Call 'draw' and 'frame'
      maximize_output(wks, True)

;=====================================
; END MID-SCRIPT PLOTTING ROUTINE 
;=====================================

    end if ; Model levels

   end if  ; Plotting routine 

   ; Fill 'u_dz0' array with values from 'difv_tot' (vertical diffusive tendency)
     v_dz0(ilev,:,:,ct)   = dvt_int(:,:)

;===================================================
; Calculate mean (azimuthally-averaged) variables 
;===================================================

     do irad = 0, dimsizes(radii)-1
      mean_vtan0(ilev,irad,:,ct)    = avg(vt_int(irad,:))     ; [ {v} ]
      mean_vrad0(ilev,irad,:,ct)    = avg(vr_int(irad,:))     ; [ {u} ]
      mean_vvel0(ilev,irad,:,ct)    = avg(w_int(irad,:))      ; [ {w} ]
      mean_vort0(ilev,irad,:,ct)    = avg(vort_int(irad,:))   ; [ {zeta} ]
      mean_avo0(ilev,irad,:,ct)     = avg(avo_int(irad,:))    ; [ {zeta + f} ]
      mean_difv0(ilev,irad,:,ct)    = avg(dvt_int(irad,:))    ; [ {diff_xz + diff_yz} ]
      mean_rho0(ilev,irad,:,ct)     = avg(rho_int(irad,:))    ; [ {rho} ]
      mean_kmh0(ilev,irad,:,ct)     = avg(kmh_int(irad,:))    ; [ {k_mh} ]
      mean_cor0(ilev,irad,:,ct)     = avg(f_int(irad,:))      ; [ {f} ]
    end do

    fcor = mean_cor0 * mean_vtan0

  ; Read 'v_int' etc into larger arrays 
    u_int0(ilev,:,:,ct)    = (/u_int(:,:)/)
    v_int0(ilev,:,:,ct)    = (/v_int(:,:)/)
    w_int0(ilev,:,:,ct)    = (/w_int(:,:)/)
    vort_int0(ilev,:,:,ct) = (/vort_int(:,:)/)
    avo_int0(ilev,:,:,ct)  = (/avo_int(:,:)/)
    f_int0(ilev,:,:,ct)    = (/f_int(:,:)/)
    kmh_int0(ilev,:,:,ct)  = (/kmh_int(:,:)/)
    rho_int0(ilev,:,:,ct)  = (/rho_int(:,:)/)
    difv_x0(ilev,:,:,ct)   = (/dvx_int(:,:)/)
    difv_y0(ilev,:,:,ct)   = (/dvy_int(:,:)/)
    difv_tot0(ilev,:,:,ct) = (/dvt_int(:,:)/)
    v_rad0(ilev,:,:,ct)    = (/vr_int(:,:)/)
    v_tan0(ilev,:,:,ct)    = (/vt_int(:,:)/)

  ; Tidy up
    delete([/u_int,v_int,w_int,vt_int,vr_int/])
    delete([/rho_int,dvx_int,dvy_int,f_int,kmh_int,vort_int,avo_int/])
    delete([/u_plane,v_plane,w_plane,f_plane,vort_plane,avo_plane/])
    delete([/difv_xplane,difv_yplane,kmh_plane/])

   end do			    ; End model level loop (do ilev = 0, num_l-1)

  ct  = ct + 1		    ; Counter variable (time)

; Tidy up before next iteration (all variables from 'pc' and 'pd' streams)
  delete([/lonC,latC,lon1C,lat1C,u,v,vort,avo,f0,f1/])
  delete([/lonD,latD,lon1D,lat1D,str_xz,str_yz,kmh,t,p,w,z,rho,rho_a/])

 end do     ; End time loop (do it = 0, times-1)

;====================================
; Calculate perturbation variables 
;====================================

  print_clock("Finished looping over time intervals! Now calculating perturbation terms...")

; 08/10/2019 --> STILL WITHIN SIMULATION LOOP (do st = 0, nsim-1)
  ptb_vtan0(:,:,:,:) = v_tan0(:,:,:,:)    - mean_vtan0(:,:,:,:)    ; [v']
  ptb_vvel0(:,:,:,:) = w_int0(:,:,:,:)    - mean_vvel0(:,:,:,:)    ; [w']
  ptb_vrad0(:,:,:,:) = v_rad0(:,:,:,:)    - mean_vrad0(:,:,:,:)    ; [u']
  ptb_vort0(:,:,:,:) = vort_int0(:,:,:,:) - mean_vort0(:,:,:,:)    ; [zeta']

  k_mh0(:,:,:,:)     = kmh_int0(:,:,:,:)                           ; [kmh]

;===========================================================
; Calculate term involving azimuthal derivatives [du/dy]
;===========================================================

; Calculate the azimuthal distance ('dy') between all points in array
; Arguments 1 and 2 represent the lat/lon of the first point, respectively
; Arguments 3 and 4 represent the lat/lon of the second point, respectively
; Argument 5 ('2') tells the script to avoid extra interpolated points
; Argument 6 ('3') returns the distance in metres

  dth  = dimsizes(thetas)

  dy(:,:,0:dth-2,:)  = gc_latlon( ypos_all(:,:,0:dth-2,:),xpos_all(:,:,0:dth-2,:),\
                   	  	  ypos_all(:,:,1:dth-1,:),xpos_all(:,:,1:dth-1,:), 2, 3)

; Replace zeros with positive values 
  dy  = where( dy .ne. 0, dy, 2000)

; Calculate radial wind increment in azimuthal direction
  du(:,:,0:dth-2,:)  = v_rad0(:,:,1:dth-1,:) - v_rad0(:,:,0:dth-2,:)

; Calculate [du/dy -- for horizontal diffusive tendency]
  dudy0(:,:,:,:) = du / dy

;===========================================================
; Calculate term involving radial derivatives [d(v/r)/dr]
;===========================================================

; Calculate the radial distance ('dr') between all points in array
; Arguments 1 and 2 represent the lat/lon of the first point, respectively
; Arguments 3 and 4 represent the lat/lon of the second point, respectively
; Argument 5 ('2') tells the script to avoid extra interpolated points
; Argument 6 ('3') returns the distance in metres

  drad = dimsizes(radii)

; Radial distance between all points in array
  dr(:,0:drad-2,:,:)  = gc_latlon( ypos_all(:,0:drad-2,:,:),xpos_all(:,0:drad-2,:,:),\
                        	   ypos_all(:,1:drad-1,:,:),xpos_all(:,1:drad-1,:,:), 2, 3)

; Replace zeros with positive values
  dr = where( dr .ne. 0, dr, 2000)

  r1(:,0:drad-2,:,:)  = gc_latlon( ycen_all(:,1:drad-1,:,:),xcen_all(:,1:drad-1,:,:),\
      	               	   	   ypos_all(:,1:drad-1,:,:),xpos_all(:,1:drad-1,:,:), 2, 3)

; Replace zeros with positive values
  r1 = where( r1 .ne. 0, r1, 2000)

; Calculate components of horizontal diffusive tendency source term

; Calculate [v/r]
  v_r0(:,:,:,:)		  = v_tan0(:,:,:,:) / r1

; Calculate [d(v/r)/dr]
  dvdr0(:,:,:,:)	  = v_r0(:,:,:,:) / dr

; Calculate [r]
  r_0(:,:,:,:)		  = r1

; Calculate [d/dr]
  ddr0(:,:,:,:)		  = 1 / dr

;=============================================
; Calculate absolute angular momentum (AAM)
;=============================================

; Calculate increments
  vm   = mean_vtan0(:,:,:,:)
  fm   = mean_cor0(:,:,:,:) * (r1 ^ 2)

; Calculate absolute angular momentum components (AAM)
  aam0 = vm * r1
  aam1 = 0.5 * fm

; Calculate absolute angular momentum [ rv + 0.5 * (fr ^ 2) ]
  mean_aam0(:,:,:,:)  = aam0 + aam1

;=======================================================
; Remove unncessary 3rd dimension from mean variables
;=======================================================

  mean_vtan1 = mean_vtan0(:,:,0,:)
  mean_vrad1 = mean_vrad0(:,:,0,:)
  mean_vvel1 = mean_vvel0(:,:,0,:)
  mean_vort1 = mean_vort0(:,:,0,:)
  mean_avo1  = mean_avo0(:,:,0,:)
  mean_rho1  = mean_rho0(:,:,0,:)
  mean_aam1  = mean_aam0(:,:,0,:)

;================================================
; Calculate vertical derivatives of mean terms 
;================================================

; Redefine model level array for vertical derivative calculations
  hybC   = c->hybrid_ht(:)         ; 63 model ('theta') levels
  mlevs  = dimsizes(hybC)	   ; Size of model level array
  hybC@units = "m"		   ; Model level units 

  print_clock("Calculating vertical derivatives of mean terms...")

; Loop over times       
  ft = 0 
  do it = ts0(st), tf0(st)

 ; Loop over model levels
   do ilev = 0, num_l-1	
	
    mp1 = min((/ilev+1,num_l-1/))
    mm1 = max((/ilev-1,0/))

   ; Loop over radial circles 
     do irad = 0, dimsizes(radii)-1

   ;=========================================================
   ; Calculate vertical derivative of mean tangential wind
   ;=========================================================

   ; [ d{v}/dz ]
     mean_dvdz0(ilev,irad,ft)      = (mean_vtan1(mp1,irad,ft) - \
                                      mean_vtan1(mm1,irad,ft)) / \
                                     (hybC(mp1) - hybC(mm1))

   ; End loop over radial circles (do irad = 0, ...)
     end do

 ; End loop over model levels (do ilev = 0, num_l-1)
   end do

 ; Advance counter variable
   ft = ft + 1 

; End loop over times (do it = ts0(st), tf0(st) )
  end do

;================================================
; Calculate vertical derivatives of eddy terms 
;================================================

  print_clock("Calculating vertical derivatives of eddy terms...")

; Initialise counter variable 
  ft = 0

  do it = ts0(st), tf0(st)

   do ilev = 0, num_l-1

    mp1 = min((/ilev+1,num_l-1/))
    mm1 = max((/ilev-1,0/))

     do irad = 0, dimsizes(radii)-1

      do th0 = 0, dimsizes(thetas)-1

    ;==========================================================
    ; Calculate vertical derivative of eddy tangential wind
    ;==========================================================

    ; [ dv'/dz ]
      ptb_dvdz0(ilev,irad,th0,ft) = (ptb_vtan0(mp1,irad,th0,ft) - \
                                     ptb_vtan0(mm1,irad,th0,ft)) / \
                                    (hybC(mp1) - hybC(mm1))

      end do   	            ; End loop over azimuth angles (do th0 = 0, ...)

     end do                 ; End loop over radial circles (do irad = 0, ...)

   end do                   ; End loop over model levels (do ilev = 0, num_l-1)

 ; Advance counter variable 
   ft = ft + 1

  end do                    ; End loop over times (do it = ts0, tf0)

;==========================================================
; Finally, calculate source terms and azimuthal averages 
;==========================================================

  print_clock("Calculating source terms and azimuthal averages...")

;=============================================================
; Mean radial influx of absolute vertical vorticity: RHS (1)
;=============================================================

  rhs1_st = -( mean_vrad1 * mean_avo1 )

;===============================================================
; Mean vertical advection of mean tangential momentum: RHS (2)
;===============================================================

  rhs2_st = -( mean_vvel1 * mean_dvdz0 )

;=======================================
; Eddy radial vorticity flux: RHS (3)
;=======================================

  ptb_vez1 = -(ptb_vrad0 * ptb_vort0)

; Calculate the azimuthal average of 'ptb_vez1' and input into 'rhs3_st'
  rhs3_st  = dim_avg_n(ptb_vez1(:,:,:,:), 2)

;===========================================================
; Vertical advection of eddy tangential momentum: RHS (4)
;===========================================================

  ptb_vev1 = -(ptb_vvel0 * ptb_dvdz0)

; Calculate the azimuthal average of 'ptb_vev1' and input into 'rhs4_st'
  rhs4_st  = dim_avg_n(ptb_vev1(:,:,:,:), 2)

;===================================================================================
; Combined mean tangential/vert diffusive tendency of tangential momentum: RHS (5)
;===================================================================================

; Stress tensor
  tau_ry1   = dim_avg_n( ( k_mh0 * ( ( (1 / r_0) * dudy0) + dvdr0 ) ), 2)

;=============================================================================
; Use stress tensor to calculate horizontal diffusive tendency source term
;=============================================================================

; Azimuthally average the 'r' and 'd/dr' expressions
  r_1   = dim_avg_n(r_0(:,:,:,:),2)
  ddr1  = dim_avg_n(ddr0(:,:,:,:),2)

  v_dh0 = ( ( 1/( r_1 ^ 2 * mean_rho1) ) * \
          (ddr1 * ( r_1 ^ 2 * mean_rho1 * tau_ry1) ) )

  rhs5a_st = v_dh0

; Vertical diffusive tendency source term has already been calculated (L827)
; Now calculate the azimuthal average using 'dim_avg_n'
  rhs5b_st = dim_avg_n( v_dz0(:,:,:,:),2 )

;=================================================
; Tidy up and produce final arrays for plotting 
;=================================================

; Calculating the terms from Montgomery et al. (2018) - their Eq.(4)
; Ignore the perturbation PGF term because it is several orders of magnitude smaller

; LHS (1) : d{v}/dt          - Mean tangential wind tendency

; RHS (1) : -{u} * {f+vort}  - Mean radial advection
; RHS (2) : -{w} * (d{v}/dz) - Mean vertical advection
; RHS (3) : -{u' * vort'}    - Eddy radial advection
; RHS (4) : -{w' * (dv'/dz)} - Eddy vertical advection
; RHS (5) : {D_v}            - Combined diffusive and PBL tendency

; Define new arrays for terms in momentum budget
  rhs1A  = new((/num_l,dimsizes(radii),nts/),float)   ; -{u} * {f+vort}
  rhs2A  = new((/num_l,dimsizes(radii),nts/),float)   ; -{w} * (d{v}/dz)
  rhs3A  = new((/num_l,dimsizes(radii),nts/),float)   ; -{u' * vort'}
  rhs4A  = new((/num_l,dimsizes(radii),nts/),float)   ; -{w' * (dv'/dz)}
  rhs5aA = new((/num_l,dimsizes(radii),nts/),float)   ; {D_r}
  rhs5bA = new((/num_l,dimsizes(radii),nts/),float)   ; {D_z}
  rhsTA  = new((/num_l,dimsizes(radii),nts/),float)   ; Sum of RHS terms

;========================================
; Make units consistent for all terms 
;========================================

  rhs1A = (/rhs1_st/)
  rhs1A = rhs1A * 3600
  rhs1A@units = "m s~S~-1~N~ h~S~-1~N~"

  rhs2A = (/rhs2_st/)
  rhs2A = rhs2A * 3600
  rhs2A@units = "m s~S~-1~N~ h~S~-1~N~"
       
  rhs3A = (/rhs3_st/)
  rhs3A = rhs3A * 3600
  rhs3A@units = "m s~S~-1~N~ h~S~-1~N~"

  rhs4A = (/rhs4_st/)
  rhs4A = rhs4A * 3600
  rhs4A@units = "m s~S~-1~N~ h~S~-1~N~"

  rhs5aA = (/rhs5a_st/)
  rhs5aA = rhs5aA * 3600 
  rhs5aA@units = "m s~S~-1~N~ h~S~-1~N~"

  rhs5bA = (/rhs5b_st/)
  rhs5bA = rhs5bA * 3600
  rhs5bA@units = "m s~S~-1~N~ h~S~-1~N~"

  hyb_plot   = hybC / 1000
  hyb_plot@units = "km"

; Add metadata for all the other important arrays (instantaneous)
  rhs1A!0    = "lev"
  rhs1A&lev  = hyb_plot
  rhs1A!1    = "rad"
  rhs1A&rad  = radii
  rhs1A!2    = "time"

  rhs2A!0    = "lev"
  rhs2A&lev  = hyb_plot
  rhs2A!1    = "rad"
  rhs2A&rad  = radii
  rhs2A!2    = "time"

  rhs3A!0    = "lev"
  rhs3A&lev  = hyb_plot
  rhs3A!1    = "rad"
  rhs3A&rad  = radii
  rhs3A!2    = "time"

  rhs4A!0    = "lev"
  rhs4A&lev  = hyb_plot
  rhs4A!1    = "rad"
  rhs4A&rad  = radii
  rhs4A!2    = "time"

  rhs5aA!0   = "lev"
  rhs5aA&lev = hyb_plot
  rhs5aA!1   = "rad"
  rhs5aA&rad = radii
  rhs5aA!2   = "time"

  rhs5bA!0   = "lev"
  rhs5bA&lev = hyb_plot
  rhs5bA!1   = "rad"
  rhs5bA&rad = radii
  rhs5bA!2   = "time"

  rhsTA!0    = "lev"
  rhsTA&lev  = hyb_plot
  rhsTA!1    = "rad"
  rhsTA&rad  = radii
  rhsTA!2    = "time"

;==============================================
; Switch order of dimensions before plotting
;==============================================

 ; 'gsn_csm_contour' plots the leftmost (rightmost) dimension on the y-axis (x-axis)
 ; Switch 'rad' with 'time' so that 'rad' is the rightmost dimension

 ; Mean radial wind -- also add metadata
   mean_vrad       = (/mean_vrad1(lev|:,time|:,rad|:)/)
   mean_vrad!0     = "lev"
   mean_vrad&lev   = hyb_plot
   mean_vrad!1     = "time"
   mean_vrad!2     = "rad"
   mean_vrad&rad   = radii
   mean_vrad@units = "m s~S~-1~N~"

 ; Mean tangential wind -- also add metadata
   mean_vtan       = (/mean_vtan1(lev|:,time|:,rad|:)/)
   mean_vtan!0     = "lev"
   mean_vtan&lev   = hyb_plot
   mean_vtan!1     = "time"
   mean_vtan!2     = "rad"
   mean_vtan&rad   = radii
   mean_vtan@units = "m s~S~-1~N~"

 ; Mean vertical velocity -- also add metadata
   mean_vvel       = (/mean_vvel1(lev|:,time|:,rad|:)/)
   mean_vvel!0     = "lev"
   mean_vvel&lev   = hyb_plot
   mean_vvel!1     = "time"
   mean_vvel!2     = "rad"
   mean_vvel&rad   = radii
   mean_vvel@units = "m s~S~-1~N~"

 ; Mean absolute angular momentum
   mean_aam        = (/mean_aam1(lev|:,time|:,rad|:)/)
   mean_aam!0      = "lev"
   mean_aam&lev    = hyb_plot
   mean_aam!1      = "time"
   mean_aam!2      = "rad"
   mean_aam&rad    = radii
   mean_aam        = mean_aam / (10 ^ 6)
   mean_aam@units  = "10~S~6~N~ m~S~2~N~ s~S~-1~N~"

;==================================
; Reorder arrays before plotting
;==================================
       
   rhs1  	  = rhs1A(lev|:,time|:,rad|:)
   rhs2  	  = rhs2A(lev|:,time|:,rad|:)
   rhs3  	  = rhs3A(lev|:,time|:,rad|:)
   rhs4  	  = rhs4A(lev|:,time|:,rad|:)
   rhs5a 	  = rhs5aA(lev|:,time|:,rad|:)
   rhs5b 	  = rhs5bA(lev|:,time|:,rad|:)

;=========================================
; Finally, smooth data before plotting
;=========================================

 ; Mean hadv, mean vadv, eddy hadv, eddy vadv, h-diff, v-diff
   rhs1  	  = smth9_Wrap(rhs1, 0.5, 0.5, True)
   rhs2  	  = smth9_Wrap(rhs2, 0.5, 0.5, True)
   rhs3  	  = smth9_Wrap(rhs3, 0.5, 0.5, True)
   rhs4  	  = smth9_Wrap(rhs4, 0.5, 0.5, True)
   rhs5a 	  = smth9_Wrap(rhs5a, 0.5, 0.5, True)
   rhs5b 	  = smth9_Wrap(rhs5b, 0.5, 0.5, True)

 ; Calculate total RHS tendency
   rhsT       = rhs1 + rhs2 + rhs3 + rhs4 + rhs5a + rhs5b
   rhsT!0     = "lev"
   rhsT&lev   = hyb_plot
   rhsT!1     = "time"
   rhsT!2     = "rad"
   rhsT&rad   = radii
   rhsT@units = "m s~S~-1~N~ h~S~-1~N~"

 ; Calculate sum of all mean terms (vadv + hadv + diff)
 ; EDIT -- start from here, and output these terms (write out, and plot)
   mean       = rhs1 + rhs2 ; + rhs5a + rhs5b
   mean!0     = "lev"
   mean&lev   = hyb_plot
   mean!1     = "time"
   mean!2     = "rad"
   mean&rad   = radii
   mean@units = "m s~S~-1~N~ h~S~-1~N~"

 ; Calculate sum of all eddy terms (vadv + hadv)
   eddy       = rhs3 + rhs4
   eddy!0     = "lev"
   eddy&lev   = hyb_plot
   eddy!1     = "time"
   eddy!2     = "rad"
   eddy&rad   = radii
   eddy@units = "m s~S~-1~N~ h~S~-1~N~"

 ; Calculate total horizontal advection
   hadv       = rhs1 + rhs3
   hadv!0     = "lev"
   hadv&lev   = hyb_plot
   hadv!1     = "time"
   hadv!2     = "rad"
   hadv&rad   = radii
   hadv@units = "m s~S~-1~N~ h~S~-1~N~"

 ; Calculate total vertical advection
   vadv       = rhs2 + rhs4
   vadv!0     = "lev"
   vadv&lev   = hyb_plot
   vadv!1     = "time"
   vadv!2     = "rad"
   vadv&rad   = radii
   vadv@units = "m s~S~-1~N~ h~S~-1~N~"

 ; Calculate total diffusive tendency
   diff       = rhs5a + rhs5b
   diff!0     = "lev"
   diff&lev   = hyb_plot
   diff!1     = "time"
   diff!2     = "rad"
   diff&rad   = radii
   diff@units = "m s~S~-1~N~ h~S~-1~N~"

 ; Fill 'plot_arr' with values from time averages
 ; 08/10/2019 --> edited the arrays below from 'mean' to 'instantaneous'
   mbgt_arr0  = (/rhs1, rhs2, rhs3, rhs4, \
                  rhs5a, rhs5b, mean_vrad, \
                  mean_vtan, mean_vvel, mean_aam/)

 ; Fill composite array with values from terms in 'mbgt_arr0'
 ; mbgt_arr  = (/ndiag, nlev, ntimes, nr/)
 ; plot_arr0 = (/nlev, nr, ndiag, nsim, nts/) where 'ntimes' = 'nsim' * 'nts'
   do mb = 0, ndiag-1 
    plot_arr0(:,:,:,mb,st) = mbgt_arr0(mb,:,:,:)
   end do 

 ; Tidy up before next iteration
   delete([/utc_date,year,month,day,hour,minute,second,date_str,time_str,time_arr,title_arr/])
   delete([/numINFO_c,numINFO_d,numINPUT_c,time_c,numINPUT_d,time_d/])
   delete([/storm_rel,u_cyc,v_cyc,lat_arr,lon_arr,centre/])

  end do		; End loop over sets of dates (do st = 0, nsim-1)

  print_clock("Finished looping over input data! Now calculating averages...")

  if (typ .eq. "test") then
   exit()
  end if 

;=================================================================================
; Before plotting, average over all chosen time periods (unless analysing ERCs)
;=================================================================================

; If we're analysing ERCs rather than VCs, write out data at each time individually 
; 15/01/2020 --> add this step for ERCs specifically 
  if (typ .eq. "erc1" .or. typ .eq. "erc2" .or. typ .eq. "erc3" .or. typ .eq. "test") then 
   
 ; Extract individual forcing terms from 'plot_arr'
   rhs1_c0  = plot_arr0(:,:,:,0,:)
   rhs2_c0  = plot_arr0(:,:,:,1,:)
   rhs3_c0  = plot_arr0(:,:,:,2,:)
   rhs4_c0  = plot_arr0(:,:,:,3,:)
   rhs5a_c0 = plot_arr0(:,:,:,4,:)
   rhs5b_c0 = plot_arr0(:,:,:,5,:)
   vrad_c0  = plot_arr0(:,:,:,6,:)
   vtan_c0  = plot_arr0(:,:,:,7,:)
   vvel_c0  = plot_arr0(:,:,:,8,:)
   aam_c0   = plot_arr0(:,:,:,9,:)

 ; Calculate additional terms before averaging                                       

 ; Total RHS                                                                         
   rhsT_c0    = rhs1_c0 + rhs2_c0 + rhs3_c0 + rhs4_c0 + rhs5b_c0 + rhs5a_c0
   copy_VarCoords(rhs1_c0,rhsT_c0)
   copy_VarMeta(rhs1_c0,rhsT_c0)
   rhsT_c0@name = "Total RHS"

 ; All mean terms                                                                    
   mean_c0    = rhs1_c0 + rhs2_c0
   copy_VarCoords(rhs1_c0,mean_c0)
   copy_VarMeta(rhs1_c0,mean_c0)
   mean_c0@name = "Sum of mean terms"

 ; All eddy terms                                                                    
   eddy_c0    = rhs3_c0 + rhs4_c0
   copy_VarCoords(rhs2_c0,eddy_c0)
   copy_VarMeta(rhs2_c0,eddy_c0)
   eddy_c0@name = "Sum of eddy terms"

 ; Radial vorticity flux (eddy + mean)                                       
   vflux_c0   = rhs1_c0 + rhs3_c0
   copy_VarCoords(rhs1_c0,vflux_c0)
   copy_VarMeta(rhs1_c0,vflux_c0)
   vflux_c0@name = "Total radial vorticity flux"

 ; Vertical advection (eddy + mean)                                          
   vadv_c0   = rhs2_c0 + rhs4_c0
   copy_VarCoords(rhs1_c0,vadv_c0)
   copy_VarMeta(rhs1_c0,vadv_c0)
   vflux_c0@name = "Total vertical advection"

;=======================================================
; Add metadata before writing out (does this matter?)
;=======================================================

   rhs1_c0!0      = "lev"
   rhs1_c0&lev    = hyb_plot
   rhs1_c0!2      = "rad"
   rhs1_c0&rad    = radii

   rhs2_c0!0      = "lev"
   rhs2_c0&lev    = hyb_plot
   rhs2_c0!2      = "rad"
   rhs2_c0&rad    = radii

   rhs3_c0!0      = "lev"
   rhs3_c0&lev    = hyb_plot
   rhs3_c0!2      = "rad"
   rhs3_c0&rad    = radii

   rhs4_c0!0      = "lev"
   rhs4_c0&lev    = hyb_plot
   rhs4_c0!2      = "rad"
   rhs4_c0&rad    = radii

   rhs5a_c0!0     = "lev"
   rhs5a_c0&lev   = hyb_plot
   rhs5a_c0!2     = "rad"
   rhs5a_c0&rad   = radii

   rhs5b_c0!0     = "lev"
   rhs5b_c0&lev   = hyb_plot
   rhs5b_c0!2     = "rad"
   rhs5b_c0&rad   = radii

   rhsT_c0!0      = "lev"
   rhsT_c0&lev    = hyb_plot
   rhsT_c0!2      = "rad"
   rhsT_c0&rad    = radii

   mean_c0!0      = "lev"
   mean_c0&lev    = hyb_plot
   mean_c0!2      = "rad"
   mean_c0&rad    = radii

   eddy_c0!0      = "lev"
   eddy_c0&lev    = hyb_plot
   eddy_c0!2      = "rad"
   eddy_c0&rad    = radii

   vflux_c0!0     = "lev"
   vflux_c0&lev   = hyb_plot
   vflux_c0!2     = "rad"
   vflux_c0&rad   = radii

   vadv_c0!0      = "lev"
   vadv_c0&lev    = hyb_plot
   vadv_c0!2      = "rad"
   vadv_c0&rad    = radii

   vrad_c0!0      = "lev"
   vrad_c0&lev    = hyb_plot
   vrad_c0!2      = "rad"
   vrad_c0&rad    = radii

   vtan_c0!0      = "lev"
   vtan_c0&lev    = hyb_plot
   vtan_c0!2      = "rad"
   vtan_c0&rad    = radii

   vvel_c0!0      = "lev"
   vvel_c0&lev    = hyb_plot
   vvel_c0!2      = "rad"
   vvel_c0&rad    = radii

   aam_c0!0       = "lev"
   aam_c0&lev     = hyb_plot
   aam_c0!2       = "rad"
   aam_c0&rad     = radii

 ; Troubleshooting 
   if (plt .eq. 0) then 
    exit()
   end if 

;=====================================================
; Loop over times and write out values at each time 
;=====================================================

 ; Set counter variable (time)
   ct = 0

   do it = ts0(0), tf0(0)

    print("Working on time: T+"+it)

  ; Additional metadata for 'write_matrix' procedure [21 * 59]
    ncol           = nr-1
    nrow           = nlev
    fmt            = nrow + "f10.3"

  ; Resources for writing out
    opts_mat       = True
    opts_mat@row   = True
    opts_mat@title = ""

  ; Name all variables we want to write
    opt_arr = (/"rhs1", "rhs2", "rhs3", "rhs4", "rhs5a", "rhs5b", "rhsT", \
                "mean", "eddy", "vflux", "vadv", "vrad", "vtan", "vvel", "aam"/)

  ; Array of all composite terms
    mbgt_arr = (/rhs1_c0, rhs2_c0, rhs3_c0, rhs4_c0, rhs5a_c0, rhs5b_c0, rhsT_c0, \
                 mean_c0, eddy_c0, vflux_c0, vadv_c0, vrad_c0, vtan_c0, vvel_c0, aam_c0/)

  ; Call procedure and write time-averaged
    print_clock("Writing 2D arrays to matrices for further analysis...")
    do mb = 0, dimsizes(opt_arr)-1

     print("Working on "+opt_arr(mb) )

   ; Set file path, etc
     t_dir  = "./"
     t_out  = opt_arr(mb)+"_sc"+sc+"_group"+grp+"_"+w0+"w_"+hy3+"km_"+typ0+"_T"+it+".txt"

     t_path = t_dir + t_out
     print_clock("Writing to file: "+t_path)
     opts_mat@fout = t_path

   ; Check for existing files and then write
     if (typ .ne. "test") then 
      system("/bin/rm -f " + t_path)
      write_matrix( mbgt_arr(mb,:,ct,0:nr-2,0), fmt, opts_mat)
     end if 

    end do
  ; End loop over diagnostics 

    ct = ct + 1
  ; Advance counter variable (time)

   end do
 ; End loop over times 

 ; Exit the script 
   if (plt .ne. 1) then 

   ; Tidy up text files by deleting indices in first column ('0,1,2,3,...')
     system("sed -i -r 's/.{10}//' ./*"+typ+"_T*.txt")

   ; Move text files to correct directory 
     system("mv ./*"+typ+"_T*.txt $ar/text/tan/")

   ; Now, exit script  
     exit()

 ; Or, continue to plotting (delete existing arrays before continuing)
   else
    delete([/rhs1_c0,rhs2_c0,rhs3_c0,rhs4_c0,rhs5a_c0,rhs5b_c0,vrad_c0,vtan_c0,vvel_c0,aam_c0/])
    delete([/mbgt_arr,rhsT_c0,mean_c0,eddy_c0,vflux_c0,vadv_c0/])
   end if 

  end if 

; Extract individual forcing terms from 'plot_arr' and average over times to create 3D arrays
  rhs1_c0  = dim_avg_n_Wrap(plot_arr0(:,:,:,0,:), 1)
  rhs2_c0  = dim_avg_n_Wrap(plot_arr0(:,:,:,1,:), 1)
  rhs3_c0  = dim_avg_n_Wrap(plot_arr0(:,:,:,2,:), 1)
  rhs4_c0  = dim_avg_n_Wrap(plot_arr0(:,:,:,3,:), 1)
  rhs5a_c0 = dim_avg_n_Wrap(plot_arr0(:,:,:,4,:), 1)
  rhs5b_c0 = dim_avg_n_Wrap(plot_arr0(:,:,:,5,:), 1)
  vrad_c0  = dim_avg_n_Wrap(plot_arr0(:,:,:,6,:), 1)
  vtan_c0  = dim_avg_n_Wrap(plot_arr0(:,:,:,7,:), 1)
  vvel_c0  = dim_avg_n_Wrap(plot_arr0(:,:,:,8,:), 1)
  aam_c0   = dim_avg_n_Wrap(plot_arr0(:,:,:,9,:), 1)

; Now average over all times to create 2D arrays 
  rhs1_c   = dim_avg_n_Wrap(rhs1_c0, 2)
  rhs2_c   = dim_avg_n_Wrap(rhs2_c0, 2)
  rhs3_c   = dim_avg_n_Wrap(rhs3_c0, 2)
  rhs4_c   = dim_avg_n_Wrap(rhs4_c0, 2)
  rhs5a_c  = dim_avg_n_Wrap(rhs5a_c0, 2)
  rhs5b_c  = dim_avg_n_Wrap(rhs5b_c0, 2)
  vrad_c   = dim_avg_n_Wrap(vrad_c0, 2)
  vtan_c   = dim_avg_n_Wrap(vtan_c0, 2)
  vvel_c   = dim_avg_n_Wrap(vvel_c0, 2)
  aam_c    = dim_avg_n_Wrap(aam_c0, 2)

  rhs1_c@name  = "Mean radial vorticity flux"
  rhs2_c@name  = "Mean vertical advection of mean tangential momentum"
  rhs3_c@name  = "Eddy radial vorticity flux"
  rhs4_c@name  = "Eddy vertical advection of eddy tangential momentum"
  rhs5a_c@name = "Horizontal diffusive tendency of tangential momentum"
  rhs5b_c@name = "Vertical diffusive tendency of tangential momentum"

  vrad_c@name  = "Radial wind"
  vtan_c@name  = "Tangential wind"
  vvel_c@name  = "Vertical velocity"
  aam_c@name   = "Absolute angular momentum"

;===============================================
; Calculate additional terms before averaging
;===============================================

; Total RHS
; 09/10/2019 --> EDIT TO REMOVE 'rhs5a'
  rhsT_c    = rhs1_c + rhs2_c + rhs3_c + rhs4_c + rhs5b_c ; + rhs5a_c
  copy_VarCoords(rhs1_c,rhsT_c)
  copy_VarMeta(rhs1_c,rhsT_c)
  rhsT_c@name = "Total RHS"

; All mean terms
  mean_c    = rhs1_c + rhs2_c
  copy_VarCoords(rhs1_c,mean_c)
  copy_VarMeta(rhs1_c,mean_c)
  mean_c@name = "Sum of mean terms"

; All eddy terms
  eddy_c    = rhs3_c + rhs4_c
  copy_VarCoords(rhs2_c,eddy_c)
  copy_VarMeta(rhs2_c,eddy_c)
  eddy_c@name = "Sum of eddy terms"

; Radial vorticity flux (eddy + mean)
  vflux_c   = rhs1_c + rhs3_c
  copy_VarCoords(rhs1_c,vflux_c)
  copy_VarMeta(rhs1_c,vflux_c)
  vflux_c@name = "Total radial vorticity flux"

; Vertical advection (eddy + mean)
  vadv_c   = rhs2_c + rhs4_c
  copy_VarCoords(rhs1_c,vadv_c)
  copy_VarMeta(rhs1_c,vadv_c)
  vflux_c@name = "Total vertical advection"
 
;===============================
; Add metadata before ploting 
;===============================

  rhs1_c!0      = "lev"
  rhs1_c&lev    = hyb_plot
  rhs1_c!1	= "rad"
  rhs1_c&rad	= radii

  rhs2_c!0      = "lev"
  rhs2_c&lev    = hyb_plot
  rhs2_c!1	= "rad"
  rhs2_c&rad	= radii

  rhs3_c!0      = "lev"
  rhs3_c&lev    = hyb_plot
  rhs3_c!1	= "rad"
  rhs3_c&rad	= radii

  rhs4_c!0      = "lev"
  rhs4_c&lev    = hyb_plot
  rhs4_c!1	= "rad"
  rhs4_c&rad	= radii

  rhs5a_c!0     = "lev"
  rhs5a_c&lev   = hyb_plot
  rhs5a_c!1	= "rad"
  rhs5a_c&rad	= radii

  rhs5b_c!0     = "lev"
  rhs5b_c&lev   = hyb_plot
  rhs5b_c!1	= "rad"
  rhs5b_c&rad	= radii

  rhsT_c!0      = "lev"
  rhsT_c&lev    = hyb_plot
  rhsT_c!1	= "rad"
  rhsT_c&rad	= radii

  mean_c!0      = "lev"
  mean_c&lev    = hyb_plot
  mean_c!1	= "rad"
  mean_c&rad	= radii

  eddy_c!0      = "lev"
  eddy_c&lev    = hyb_plot
  eddy_c!1	= "rad"
  eddy_c&rad	= radii

  vflux_c!0     = "lev"
  vflux_c&lev   = hyb_plot
  vflux_c!1     = "rad"
  vflux_c&rad   = radii

  vadv_c!0      = "lev"
  vadv_c&lev    = hyb_plot
  vadv_c!1      = "rad"
  vadv_c&rad    = radii

  vrad_c!0      = "lev"
  vrad_c&lev    = hyb_plot
  vrad_c!1	= "rad"
  vrad_c&rad	= radii

  vtan_c!0      = "lev"
  vtan_c&lev    = hyb_plot
  vtan_c!1	= "rad"
  vtan_c&rad	= radii

  vvel_c!0      = "lev"
  vvel_c&lev    = hyb_plot
  vvel_c!1	= "rad"
  vvel_c&rad	= radii

  aam_c!0       = "lev"
  aam_c&lev     = hyb_plot
  aam_c!1	= "rad"
  aam_c&rad	= radii

;========================================================================
; Write values to output matrices (text files) using in-built function
;========================================================================

; Additional metadata for 'write_matrix' procedure [21 * 59]
  ncol           = nr-1
  nrow           = nlev
  fmt            = nrow + "f10.3"

; Resources for writing out 
  opts_mat       = True
  opts_mat@row   = True
  opts_mat@title = ""

; Name all variables we want to write 
  opt_arr = (/"rhs1", "rhs2", "rhs3", "rhs4", "rhs5a", "rhs5b", "rhsT", \
              "mean", "eddy", "vflux", "vadv", "vrad", "vtan", "vvel", "aam"/)

; Array of all composite terms
  mbgt_arr = (/rhs1_c, rhs2_c, rhs3_c, rhs4_c, rhs5a_c, rhs5b_c, rhsT_c, \
               mean_c, eddy_c, vflux_c, vadv_c, vrad_c, vtan_c, vvel_c, aam_c/)

; Call procedure and write time-averaged
  print_clock("Writing 2D arrays to matrices for further analysis...")
  do mb = 0, dimsizes(opt_arr)-1

    print("Working on "+opt_arr(mb) )

  ; Set file path, etc
    t_dir  = "./"

    if (w_check .eq. 1) then 
     if (typ .eq. "v" .or. typ .eq. "slp") then 
      t_out = opt_arr(mb)+"_sc"+sc+"_group"+grp+"_"+w0+"w_"+hy3+"km_"+typ+thr+".txt"
     elseif (typ .eq. "novc") then 
      t_out = opt_arr(mb)+"_sc"+sc+"_group"+grp+"_"+w0+"w_"+hy3+"km_"+typ0+"_"+nts0+"_"+ntf0+".txt"
     else
      t_out = opt_arr(mb)+"_sc"+sc+"_group"+grp+"_"+w0+"w_"+hy3+"km_"+typ0+".txt" 
      end if 
    else
     if (typ .eq. "v" .or. typ .eq. "slp") then
      t_out = opt_arr(mb)+"_sc"+sc+"_group"+grp+"_orig_"+typ+thr+".txt"
     else
      t_out = opt_arr(mb)+"_sc"+sc+"_group"+grp+"_orig_"+typ0+".txt"
     end if 
    end if 

    t_path = t_dir + t_out
    print_clock("Writing to file: "+t_path)
    opts_mat@fout = t_path

  ; Check for existing files and then write 
    if (typ .ne. "test") then 
     system("/bin/rm -f " + t_path)
     write_matrix( mbgt_arr(mb,:,:), fmt, opts_mat)
    end if 

  end do

;==========================
; Options for plotting 
;==========================

  if (plt .eq. 1) then 

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

    if (clr .eq. 1) then					       
     opts_mbgt@cnFillPalette		   = "ncl_default"
     opts_mbgt@cnFillColors                = (/2,8,16,28,40,52,64,76,88,100,118,-1,-1,\
                                               134,150,166,178,186,194,\
                                               202,214,226,238,250/)
     opts_mbgt@cnLevels                    = (/-150., -100., -50., -20., -10., -5., \
                                                -3., -2., -1., -0.5, -0.2, \
                                                0, 0.2, 0.5, 1., 2., 3., \
                                                5., 10., 20., 50., 100., 150./)

    elseif (clr .eq. 2) then 
     opts_mbgt@cnFillPalette               = "bdr_extra"
     opts_mbgt@cnFillColors                = (/0,1,2,3,4,5,6,7,8,9,10,-1,-1,\
                                               13,14,15,16,17,18,19,20,21,22,23/)
     opts_mbgt@cnLevels                    = (/-150., -100., -50., -20., -10., -5., \
                                                -3., -2., -1., -0.5, -0.2, \
                                                0, 0.2, 0.5, 1., 2., 3., \
                                                5., 10., 20., 50., 100., 150./)

    elseif (clr .eq. 3) then
     opts_mbgt@cnFillPalette               = "amwg256"
     opts_mbgt@cnFillColors                = (/0,6,12,18,27,36,45,54,66,78,90,-1,-1,\
                                               148,156,166,178,186,194,\
                                               202,214,226,238,250/)
     opts_mbgt@cnLevels                    = (/-150., -100., -50., -20., -10., -5., \
                                                -3., -2., -1., -0.5, -0.2, \
                                                0, 0.2, 0.5, 1., 2., 3., \
                                                5., 10., 20., 50., 100., 150./)
    elseif (clr .eq. 4) then
     opts_mbgt@cnFillPalette               = "GMT_polar"
     opts_mbgt@cnFillColors                = (/3,5,7,-1,11,14,18/)
     opts_mbgt@cnLevels                    = (/-20,-5,0,0.5,5,20/)
    elseif (clr .eq. 5) then 
     opts_mbgt@cnFillPalette               = "BlueDarkRed18"
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
;    opts_mbgt@tiXAxisFontHeightF           = 0.025   ; X-axis title font height 

  ; Turn off tickmarks on top x-axis and right y-axis
    opts_mbgt@tmXTOn                       = "False"
    opts_mbgt@tmYROn                       = "False"

  ; Do not maximise the individual panels 
    opts_mbgt@gsnMaximize                  = False
    opts_mbgt@gsnAddCyclic                 = False


    opts_mbgt@tmXBMode                     = "Explicit"

    if (rad0 .eq. 0) then
     opts_mbgt@tmXBValues                   = (/0.0,0.25,0.5,0.75,1.0,1.25,1.5/)
     opts_mbgt@tmXBLabels                   = (/"0","25","50","75","100","125","150"/)
    elseif (rad0 .eq. 1) then
     opts_mbgt@tmXBValues                   = (/0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0/)
     opts_mbgt@tmXBLabels                   = (/"0","25","50","75","100","125","150","175","200"/)
    end if

    opts_mbgt@tmYLMode			   = "Explicit"
    opts_mbgt@tmYLValues		   = fspan(2.0, 18.0, 9)
    opts_mbgt@tmYLLabels		   = (/"2.0","4.0","6.0","8.0","10.0",\
    					       "12.0","14.0","16.0","18.0"/)
;    opts_mbgt@tmYLLabelFontHeightF	   = 0.02
    opts_mbgt@gsnYAxisIrregular2Linear	   = True	; Linearise the y-axis 

  ; y-axis limits
    opts_mbgt@trYMinF			   = 0.0
    opts_mbgt@trYMaxF			   = 17.6

  ; x-axis limits
    opts_mbgt@trXMinF                      = 0.0
    opts_mbgt@trXMaxF                      = 1.0

  ; Resources to overlay radial wind (following M18)
    opts_vrad				   = True
    opts_vrad@cnFillOn			   = False
    opts_vrad@cnLineColor		   = "black"
    opts_vrad@cnLevelSelectionMode	   = "ExplicitLevels"
    opts_vrad@cnLevels                     = (/-1.2, 1.2/)
    opts_vrad@gsnContourNegLineDashPattern = 16
    opts_vrad@gsnContourPosLineDashPattern = 0
    opts_vrad@cnLineThicknessF             = 2.5

    opts_vrad@cnInfoLabelOn                = False
    opts_vrad@gsnDraw			   = False
    opts_vrad@gsnFrame			   = False
    opts_vrad@cnLineLabelsOn		   = False
    opts_vrad@tiMainString                 = ""
    opts_vrad@gsnLeftString                = ""
    opts_vrad@gsnRightString               = ""

  ; Resources to overlay vertical velocity (following M18)
    opts_vvel 	 	 	  	   = True 
    opts_vvel@cnFillOn                     = False

  ; Set contour line colour to match M18 if we choose their colour scheme
    if (clr .eq. 4 .or. clr .eq. 5) then
     opts_vvel@cnLineColor		   = "yellow"
    else
     opts_vvel@cnLineColor                 = "grey45"
    end if 

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
    opts_zero@cnLineLabelsOn               = False
    opts_zero@tiMainString                = ""
    opts_zero@gsnLeftString               = ""
    opts_zero@gsnRightString              = ""

;==================================================================
; Also add resources for mean quantities (u,v,w) when 'full' = 2
;==================================================================

  ; Reading in colour table using 'cnFillPalette' removes the first 2 colours (black/white)
  ; This means that the indices are all 2 values lower than they'd normally be

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
    opts_vrad2@vpXF			   = 0.05
    opts_vrad2@vpYF			   = 0.95
    opts_vrad2@vpWidthF                    = 0.16
    opts_vrad2@vpHeightF		   = 0.16

  ; Tangential wind 
    opts_vtan2				   = opts_mbgt
    delete(opts_vtan2@cnLevels)
    delete(opts_vtan2@cnFillColors)
    opts_vtan2@cnFillPalette		   = "radar_new"

    opts_vtan2@cnLevels                    = (/5.0, 10.0, 15.0, 20.0,\
                                               25.0, 30.0, 35.0, 40.0, \
                                               45.0, 50.0, 55.0, 60.0, \
                                               65.0/)
    opts_vtan2@cnFillColors                = (/-1,0,1,2,3,4,5,6,7,8,9,10,11,12/)

    opts_vtan2@lbLabelBarOn                = True
    opts_vtan2@lbBoxEndCapStyle            = "TriangleBothEnds"
    opts_vtan2@lbOrientation               = "Vertical"
    opts_vtan2@vpXF                        = 0.30
    opts_vtan2@vpYF                        = 0.95
    opts_vtan2@vpWidthF                    = 0.16
    opts_vtan2@vpHeightF                   = 0.16

  ; Vertical velocity
    opts_vvel2                             = opts_mbgt
    delete(opts_vvel2@cnLevels)
    delete(opts_vvel2@cnFillColors)
    opts_vvel2@cnFillPalette               = "n11"
    opts_vvel2@cnLevels                    = (/0.1, 0.2, 0.3, 0.4, 0.5, 0.6, \
                                               0.7, 0.8, 1.0, 1.2, 1.5, 2.0/)
    opts_vvel2@cnFillColors                = (/-1,0,1,2,3,4,5,6,7,8,9,10,11/)
    opts_vvel2@lbLabelBarOn                = True
    opts_vvel2@lbBoxEndCapStyle            = "TriangleBothEnds"
    opts_vvel2@lbOrientation               = "Vertical"
    opts_vvel2@vpXF                        = 0.55
    opts_vvel2@vpYF                        = 0.95
    opts_vvel2@vpWidthF                    = 0.16
    opts_vvel2@vpHeightF                   = 0.16

  ; Absolute angular momentum
    opts_aam				   = opts_vtan2
    delete(opts_aam@cnLevels)
    delete(opts_aam@cnFillColors)

    opts_aam@cnLevels                      = (/0.6, 0.8, 1.0, 1.2, 1.5, \
                                               2.0, 2.5, 3.0, 3.5, \
                                               4.0, 4.5, 5.0/)
    opts_aam@cnFillColors                  = (/-1,0,1,2,3,4,5,6,7,8,9,10,11/)

    opts_aam@lbLabelBarOn                  = True
    opts_aam@lbBoxEndCapStyle              = "TriangleBothEnds"
    opts_aam@lbOrientation                 = "Vertical"
    opts_aam@vpXF                          = 0.80
    opts_aam@vpYF                          = 0.95
    opts_aam@vpWidthF                      = 0.16
    opts_aam@vpHeightF                     = 0.16

;===============================================================
; Resources to overlay panel labels without using 'gsn_panel'
;===============================================================

    txid			= new(12,graphic)
    amid	    	        = new(12,graphic)

  ; Label text resources 
    txres                       = True
    txres@txPerimOn             = True
    txres@txFontHeightF         = 0.010
    txres@txBackgroundFillColor = "White"

  ; Label position resources 
  ; amParallelPosF = -0.5/0.5 left/right edge; amOrthogonalPosF = -0.5 top edge
    amres 	   	        = True
    amres@amParallelPosF	= 0.5
    amres@amOrthogonalPosF	= -0.5
    amres@amJust		= "TopRight"

  ; Panel strings
    panel_strings               = (/"a) Radial wind", "b) Tangential wind", \
                                    "c) Vertical velocity", "d) AAM", \
                                    "e) Mean h. advection", "f) Mean v. advection", \
                                    "g) Eddy h. advection", "h) Eddy v. advection", \
                                    "i) Total mean", "j) Total eddy", \
                                    "k) Sum of RHS terms", "l) Local tendency (LHS)"/)

;=================
; Plot the data 
;=================

  ; Calculating the terms from Montgomery et al. (2018) - their Eq.(4)
  ; Ignore the perturbation PGF term because it is several orders of magnitude smaller

  ; LHS (1) : d{v}/dt          - Mean tangential momentum tendency

  ; RHS (1) : -{u} * {f+vort}  - Mean radial influx of absolute vertical vorticity
  ; RHS (2) : -{w} * (d{v}/dz) - Mean vertical advection of mean tangential momentum
  ; RHS (3) : -{u' * vort'}    - Eddy radial vorticity flux
  ; RHS (4) : -{w' * (dv'/dz)} - Vertical advection of eddy tangential momentum
  ; RHS (5) : {D_v}            - Combined diffusive and PBL tendency

  ; Output file location and type (including w and z threshold information)
    if (typ .eq. "v" .or. typ .eq. "slp") then 
     output = "$nep/nepartak/images/cyl_coords/mbgt_tan_cxz_"+typ+"_thr"+thr+"_group"+grp+\
              "_sc"+sc+"_"+hy3+"km_"+w0+"w_mean"+mean0
    else
     if (typ .eq. "novc") then 
     output = "$nep/nepartak/images/cyl_coords/mbgt_tan_cxz_"+typ0+"_"+nts+"_"+ntf+"_group"+grp+\
              "_sc"+sc+"_"+hy3+"km_"+w0+"w_mean"+mean0
     elseif (typ .eq. "erc1" .or. typ .eq. "erc2" .or. typ .eq. "test") then
      output = "$nep/nepartak/images/cyl_coords/erc_tan_cxz_"+typ+"_"\
               +ts0+"_"+tf0
     else
     output = "$nep/nepartak/images/cyl_coords/mbgt_tan_cxz_"+typ0+"_group"+grp+\
              "_sc"+sc+"_"+hy3+"km_"+w0+"w_mean"+mean0   
     end if 
    end if 
    wks    = gsn_open_wks(opt,output)

  ; Radial wind (filled contours)
    panel_vrad  = gsn_csm_contour(wks,vrad_c(:,:),opts_vrad2) ; [u]
    panel_zero  = gsn_csm_contour(wks,vrad_c(:,:),opts_zero)  ; Zero line contour
    zero0       = ColorNegDashZeroPosContour(panel_zero,"transparent","black","transparent")
    overlay(panel_vrad, zero0)
    txid(0)     = gsn_create_text(wks, panel_strings(0), txres)
    amid(0)     = gsn_add_annotation(panel_vrad, txid(0), amres)

  ; Tangential wind (filled contours)
    panel_vtan  = gsn_csm_contour(wks,vtan_c(:,:),opts_vtan2) ; [v]
    txid(1)     = gsn_create_text(wks, panel_strings(1), txres)
    amid(1)     = gsn_add_annotation(panel_vtan, txid(1), amres)

  ; Vertical velocity (filled contours)
    panel_vvel  = gsn_csm_contour(wks,vvel_c(:,:),opts_vvel2) ; [w]
    txid(2)     = gsn_create_text(wks, panel_strings(2), txres)
    amid(2)     = gsn_add_annotation(panel_vvel, txid(2), amres)

  ; Absolute angular momentum (filled contours)
    panel_aam   = gsn_csm_contour(wks,aam_c(:,:),opts_aam)   ; [AAM]
    txid(3)     = gsn_create_text(wks, panel_strings(3), txres)
    amid(3)     = gsn_add_annotation(panel_aam, txid(3), amres)

  ; Mean radial influx of absolute vertical vorticity [V_mzeta]
    opts_mbgt@vpXF       = 0.05
    opts_mbgt@vpYF       = 0.74
    opts_mbgt@vpWidthF   = 0.19
    opts_mbgt@vpHeightF	 = 0.19
    panel0      = gsn_csm_contour(wks,rhs1_c(:,:),opts_mbgt) ; [V_mzeta]
    cont0       = gsn_csm_contour(wks,rhs1_c(:,:),opts_zero) ; Draw contours

    panel0_ovr  = gsn_csm_contour(wks,vrad_c(:,:),opts_vrad) ; [u]
    panel0_ovr2 = gsn_csm_contour(wks,vvel_c(:,:),opts_vvel) ; [w]

    overlay(panel0, cont0)				     ; Overlay contours
    overlay(panel0, panel0_ovr)
    overlay(panel0, panel0_ovr2)
    txid(4)     = gsn_create_text(wks, panel_strings(4), txres)
    amid(4)     = gsn_add_annotation(panel0, txid(4), amres)
    delete([/cont0/])

  ; Mean vertical advection of mean tangential momentum [V_mv]
    opts_mbgt@vpXF       = 0.30
    opts_mbgt@vpYF       = 0.74
    opts_mbgt@vpWidthF   = 0.19
    opts_mbgt@vpHeightF  = 0.19
    panel1      = gsn_csm_contour(wks,rhs2_c(:,:),opts_mbgt) ; [V_mv]
    cont0       = gsn_csm_contour(wks,rhs2_c(:,:),opts_zero) ; Draw contours

    panel1_ovr  = gsn_csm_contour(wks,vrad_c(:,:),opts_vrad) ; [U]
    panel1_ovr2 = gsn_csm_contour(wks,vvel_c(:,:),opts_vvel) ; [w]

    overlay(panel1, cont0)                                   ; Overlay contours    
    overlay(panel1, panel1_ovr)
    overlay(panel1, panel1_ovr2)
    txid(5)     = gsn_create_text(wks, panel_strings(5), txres)
    amid(5)     = gsn_add_annotation(panel1, txid(5), amres)
    delete([/cont0/])

  ; Eddy radial vorticity flux [V_ezeta]
    opts_mbgt@vpXF       = 0.55
    opts_mbgt@vpYF       = 0.74
    opts_mbgt@vpWidthF   = 0.19
    opts_mbgt@vpHeightF  = 0.19
    panel2      = gsn_csm_contour(wks,rhs3_c(:,:),opts_mbgt) ; [V_ezeta]
    cont0       = gsn_csm_contour(wks,rhs3_c(:,:),opts_zero) ; Draw contours

    panel2_ovr  = gsn_csm_contour(wks,vrad_c(:,:),opts_vrad) ; [U]
    panel2_ovr2 = gsn_csm_contour(wks,vvel_c(:,:),opts_vvel) ; [w]

    overlay(panel2, cont0)                                   ; Overlay contours
    overlay(panel2, panel2_ovr)
    overlay(panel2, panel2_ovr2)
    txid(6)     = gsn_create_text(wks, panel_strings(6), txres)
    amid(6)     = gsn_add_annotation(panel2, txid(6), amres)
    delete([/cont0/])

  ; Vertical advection of eddy tangential momentum [V_ev]
    opts_mbgt@vpXF       = 0.80
    opts_mbgt@vpYF       = 0.74
    opts_mbgt@vpWidthF   = 0.19
    opts_mbgt@vpHeightF  = 0.19
    panel3      = gsn_csm_contour(wks,rhs4_c(:,:),opts_mbgt) ; [V_ev]
    cont0       = gsn_csm_contour(wks,rhs4_c(:,:),opts_zero) ; Draw contours

    panel3_ovr  = gsn_csm_contour(wks,vrad_c(:,:),opts_vrad) ; [U]
    panel3_ovr2 = gsn_csm_contour(wks,vvel_c(:,:),opts_vvel) ; [w]

    overlay(panel3, cont0)                                   ; Overlay contours
    overlay(panel3, panel3_ovr)
    overlay(panel3, panel3_ovr2)
    txid(7)     = gsn_create_text(wks, panel_strings(7), txres)
    amid(7)     = gsn_add_annotation(panel3, txid(7), amres)
    delete([/cont0/])

  ; Horizontal diffusive tendency of tangential momentum [V_dh]
    opts_mbgt@vpXF       = 0.05
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.19
    opts_mbgt@vpHeightF  = 0.19

    if (mean0 .eq. 1) then 
     panel4     = gsn_csm_contour(wks,mean_c(:,:),opts_mbgt) ; [V_dh]
     cont0      = gsn_csm_contour(wks,mean_c(:,:),opts_zero) ; Draw contours
    elseif (mean0 .eq. 0) then  
     panel4     = gsn_csm_contour(wks,rhs5a_c(:,:),opts_mbgt) ; [V_dh]
     cont0      = gsn_csm_contour(wks,rhs5a_c(:,:),opts_zero) ; Draw contours
    end if 

    panel4_ovr  = gsn_csm_contour(wks,vrad_c(:,:),opts_vrad) ; [U]
    panel4_ovr2 = gsn_csm_contour(wks,vvel_c(:,:),opts_vvel) ; [w]

    overlay(panel4, cont0)                                   ; Overlay contours
    overlay(panel4, panel4_ovr)
    overlay(panel4, panel4_ovr2)
    txid(8)     = gsn_create_text(wks, panel_strings(8), txres)
    amid(8)     = gsn_add_annotation(panel4, txid(8), amres)
    delete([/cont0/])

  ; Vertical diffusive tendency of tangential momentum [V_dz]
    opts_mbgt@vpXF       = 0.30
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.19
    opts_mbgt@vpHeightF  = 0.19

    if (mean0 .eq. 1) then 
     panel5     = gsn_csm_contour(wks,eddy_c(:,:),opts_mbgt)  ; [V_dz]
     cont0      = gsn_csm_contour(wks,eddy_c(:,:),opts_zero)  ; Draw contours
    elseif (mean0 .eq. 0) then 
     panel5     = gsn_csm_contour(wks,rhs5b_c(:,:),opts_mbgt) ; [V_dz]
     cont0      = gsn_csm_contour(wks,rhs5b_c(:,:),opts_zero) ; Draw contours
    end if 

    panel5_ovr  = gsn_csm_contour(wks,vrad_c(:,:),opts_vrad)  ; [U]
    panel5_ovr2 = gsn_csm_contour(wks,vvel_c(:,:),opts_vvel)  ; [w]

    overlay(panel5, cont0)                                    ; Overlay contours
    overlay(panel5, panel5_ovr)
    overlay(panel5, panel5_ovr2)
    txid(9)     = gsn_create_text(wks, panel_strings(9), txres)
    amid(9)     = gsn_add_annotation(panel5, txid(9), amres)
    delete([/cont0/])

  ; Total RHS [All V tend]
    opts_mbgt@vpXF       = 0.55
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.19
    opts_mbgt@vpHeightF  = 0.19
    panel6      = gsn_csm_contour(wks,rhsT_c(:,:),opts_mbgt) ; [All V tend]
    cont0       = gsn_csm_contour(wks,rhsT_c(:,:),opts_zero) ; Draw contours

    panel6_ovr  = gsn_csm_contour(wks,vrad_c(:,:),opts_vrad) ; [U]
    panel6_ovr2 = gsn_csm_contour(wks,vvel_c(:,:),opts_vvel) ; [w]

    overlay(panel6, cont0)                                   ; Overlay contours
    overlay(panel6, panel6_ovr)
    overlay(panel6, panel6_ovr2)
    txid(10)    = gsn_create_text(wks, panel_strings(10), txres)
    amid(10)    = gsn_add_annotation(panel6, txid(10), amres)
    delete([/cont0/])

  ; Local tendency of the mean tangential wind [V_t]
    opts_mbgt@vpXF       = 0.80
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.19
    opts_mbgt@vpHeightF  = 0.19

  ; No LHS tendency calculated here, so fill with RHS tendency as placeholder
    panel7      = gsn_csm_contour(wks,rhsT_c(:,:),opts_mbgt) ; [V_t]
    cont0       = gsn_csm_contour(wks,rhsT_c(:,:),opts_zero) ; Draw contours

    panel7_ovr  = gsn_csm_contour(wks,vrad_c(:,:),opts_vrad) ; [U]
    panel7_ovr2 = gsn_csm_contour(wks,vvel_c(:,:),opts_vvel) ; [w]

    overlay(panel7, cont0)                                   ; Overlay contours
    overlay(panel7, panel7_ovr)
    overlay(panel7, panel7_ovr2)
    txid(11)    = gsn_create_text(wks, panel_strings(11), txres)
    amid(11)    = gsn_add_annotation(panel7, txid(11), amres)
    delete([/cont0/])

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    pres			= True 
    maximize_output(wks,pres)				   ; Calls 'draw' and 'frame'

; End IF statement (choose whether to plot in addition to writing out to text files) 
  end if 

end 