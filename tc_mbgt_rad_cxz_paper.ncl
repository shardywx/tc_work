; Script to calculate the radial momentum budget as in Huang et al. (2018)

; Run using:

; ncl opt=\"png\" w0=0.3 z0=0 dist=1.2 clr=5 calc=1 rad0=0 nt=73  
; nr=22 ar=1.05 lgd0=\"tr\" new0=0 mean0=0 cn0=\"geo_sm\" sm=0 sc=5 w_check=1
; nts=36 ntf=54 grp=\"7\" plt=1 ueh=0 typ=\"test\" tc_mbgt_rad_cxz_paper.ncl 

; 'opt'   = output file format ("pdf" or "x11")
; 'dist'  = size of box following storm (degrees)
; 'clr'	  = colour map for plots (option 4 is used in Roger Smith's papers):
; 	    "ncl_default" [1],  "bdr_extra" [2], "amwg256" [3], "GMT_polar" [4], "BlueDarkRed18" [5]
; 'calc'  = centre calculated on each level (0) or equal to 950 hPa centre (1)
; 'rad0'  = radius of final plots: 150 km (0) or 200 km (1)
; 'ave'	  = calculate 3/6/12-h average (1) rather than analyse individual time interval (0)
; 'nt'	  = number of azimuth angles in cylindrical grid 
; 'wrt'	  = write out averages to text file (1) or leave alone (0 – default)
; 'grid'  = overlay NDC grid (1) or leave alone (0)
; 'typ'	  = VC phase ("sym","asym","stoa","atos")
; 'full'  = test (3), all time periods (2), shortened time (1) or shortest time (0)
; 'ueh'   = plot eddy horizontal advection of radial momentum (1) or leave out, as in M18 (0)
; 'ar'     = distance in degrees of outer radius (1.0, 2.0, etc)
; 'nr'     = number of radial circles between r = 0 and r = 'ar' (21, 41, 61, etc)

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
   ts0  = (/88/)
   tf0  = (/89/)
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
  ndiag     = 13

  plot_arr0 = new( (/nlev, nts, nr, ndiag, nsim/), "float")

;==========================================================
; Now start loop over these chosen input files and times
;==========================================================

 do st = 0, dimsizes(dat)-1 

  dt0 = tf0(st) - ts0(st)

; Find 'pc' and 'pd' files (model height level data)
  diri        = "/nfs/a319/earshar/"+dat(st)+"/"+ens0(st)
  fili_prefix = "201607"+dat(st)+"00Z_NPTK_4p4_L80_ra1t_"+ens0(st)
  fili_c1     = systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pc*.nc")
  fili_c      = diri+"/"+fili_c1
  fili_d1     =	systemfunc("cd "+diri+" ; ls "+fili_prefix+"_pd*.nc")
  fili_d      = diri+"/"+fili_d1

; Calculate number of times in each file using built-in function ('nc_times')
  numINFO_c   = nc_times(fili_c)
  numINFO_d   = nc_times(fili_d)
  numINPUT_c  = numINFO_c[0]		; Number of times in 'pc' stream
  time_c      = numINFO_c[1]		; Array of times from 'pc' stream file
  numINPUT_d  = numINFO_d[0]		; Number of times in 'pd' stream
  time_d      = numINFO_d[1]		; Array of times from 'pd' stream file 

;============================================================  
; Calculate storm motion using built-in function ('st_rm')
;============================================================  

; 'dat'     = initialisation time (02T12, 03T00, ...)
; 'ens0'    = ensemble simulation (em00, em01, ...)
; 'diri'    = path to input files (see above)
; 'dist'    = size of box following storm (degrees)
; 'mins'    = analyse 1-h (0) or 5-min (1) data

; Zonal; meridional; vector wind; times; lat; lon; centre 
  storm_rel = st_rm(dat(st), ens0(st), "$ar/text/", dist, 0)
  u_cyc     = storm_rel[0]
  v_cyc     = storm_rel[1]
  vel_cyc   = storm_rel[2]
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
  do ct = 0, numINPUT_c-1
   if (minute(ct).gt.30) then
     hour(ct) = hour(ct)+1
   end if
  end do
  delete(ct)

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
  rd	= 287.0	  ; Gas constant for dry air (J/kg/K^2)
  cp0	= 1004.0  ; Specific heat of dry air at constant pressure (J/kg/K)
  kp	= 0.286   ; For Exner function calculations (0.286)
  pref	= 1000.0  ; Reference pressure p0 (1000 hPa)

;========================================================================================
;========================================================================================
; Calculating the terms from Montgomery et al. (2018) - their Eq.(8)
; LHS (1) : d{u}/dt	     - Local tendency of the mean radial velocity
; LHS (2) : {u} * d{u}/dr    - Mean radial advection of	radial momentum
; LHS (3) : {(u' * du'/dr) + (v'/r * du'/dy)}  - Mean horiz. adv. of eddy radial momentum

; RHS (1) : -{w} * d{u}/dz   - Mean vertical advection of mean radial momentum
; RHS (2) : -{w' * d{u'}/dz} - Eddy vertical advection of eddy radial momentum
; RHS (3) : {v}^2/r + f{v} - 1/{rho} * d{p}/dr - Mean agradient force 
; RHS (4) : {v'^2/r - 1/rho * dp'/dr}	       - Eddy agradient force 
; RHS (5) : {F_r}	     - Mean radial/vertical diffusive tendency of radial momentum

; All terms except LHS (1) and RHS (5) are per unit mass (?)
;========================================================================================
;========================================================================================

; Array to hold all 'xpos' and 'ypos' position arrays
  xpos_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)
  ypos_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)

; Arrays to hold lat/lon at centre of cylindrical grid
  xcen_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)
  ycen_all  = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)

; Arrays for perturbation (eddy) variables (4D)
  ptb_vrad0 = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; [u']
  ptb_vtan0 = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; [v']
  ptb_vvel0 = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; [w']
  ptb_pres0 = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; [p']

; Additional variables (same dimensions but not eddy quantities)
  rho0      = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; [ro]

; Arrays for eddy derivative terms (4D)
; (1) du'/dr, (2) du'/dy, (3) du'/dz (4) dp'/dr
  ptb_dudr0 = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; (1)
  ptb_dudy0 = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; (2)
  ptb_dudz0 = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; (3)
  ptb_dpdr0 = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; (4)

; Arrays for eddy 1/r terms (4D)
; (1) v'/r, (2) v'^2 / r
  ptb_vt_r0  = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; (1)
  ptb_vttr0  = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; (2)

; Array for mean 1/r term (3D)
; (1) {v}^2 / r
  mean_vttr1 = new((/num_l,dimsizes(radii),nts/),float) ; (1) 

; Vertical derivative of mean radial wind [ d{u} / dz ]
  mean_dudz1          = new((/num_l,dimsizes(radii),nts/),float)
  mean_dudz1!0        = "lev"
  mean_dudz1!1        = "rad"
  mean_dudz1!2        = "time"
  mean_dudz1@description = "Azimuthally averaged 'du_dz'"
  mean_dudz1@units = "s~S~-1~N~"

;===========================================================
; Arrays to hold expanded 'v_int'/'vort_int' arrays (4D)
;===========================================================

  u_int0    = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)  ; Zonal wind
  v_int0    = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)  ; Meridional wind
  w_int0    = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)  ; Vertical velocity
  f_int0    = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)  ; Coriolis parameter
  kmh_int0  = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)  ; Horiz. diffusivity
  prs_int0  = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)  ; Pressure 
  rho_int0  = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)  ; Density
  v_rad0    = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)  ; Radial wind
  v_tan0    = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)  ; Tangential wind

  difv_x0   = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)  ; Vert. diffusion (x)
  difv_y0   = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)  ; Vert. diffusion (y)
  difv_tot0 = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)  ; Vert. diff. (tot)

  difh_x0   = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)  ; Horiz. diffusion (x)
  difh_y0   = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)  ; Horiz. diffusion (y)

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

  prs_int0@description = "Pressure"
  prs_int0@units       = "hPa"
  prs_int0!0           = "lev"
  prs_int0!1           = "rad"
  prs_int0!2           = "azi"
  prs_int0!3           = "time"

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

  difh_x0@description  = "Horizontal diffusion (x-component)"
  difh_x0@units        = "m s~S~-2~N~"
  difh_x0!0            = "lev"
  difh_x0!1            = "rad"
  difh_x0!2            = "azi"
  difh_x0!3            = "time"

  difh_y0@description  = "Horizontal diffusion (y-component)"
  difh_y0@units        = "m s~S~-2~N~"
  difh_y0!0            = "lev"
  difh_y0!1            = "rad"
  difh_y0!2            = "azi"
  difh_y0!3            = "time"

  v_rad0@description   = "Radial velocity"
  v_rad0@units         = "m s~S~-1~N~"
  v_rad0!0             = "lev"
  v_rad0!1             = "rad"
  v_rad0!2             = "azi"
  v_rad0!3             = "time"

  v_tan0@description   = "Tangential velocity"
  v_tan0@units         = "m s~S~-1~N~"
  v_tan0!0             = "lev"
  v_tan0!1             = "rad"
  v_tan0!2             = "azi"
  v_tan0!3             = "time"

;=======================================================================
; Arrays for LHS and RHS source terms containing eddy variables (4D)
;=======================================================================

; Mean horizontal advection of eddy radial momentum: LHS (3)
; {(u' * du'/dr) + (v'/r * du'/dy)}
; ptb_ueh1  = (ptb_vrad0 * ptb_dudr0) + ( (ptb_vtan0 / r) * (ptb_dudy0) )
  ptb_ueh1  = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; [L1]

; Eddy vertical advection of eddy radial momentum: RHS (2)
; [-w' * du'/dz]
; ptb_uev1  = - (ptb_vvel0 * ptb_dudz0) 
  ptb_uev1  = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; [R2]

; Eddy agradient force: RHS (3)
; [ {v'^2/r - 1/rho * dp'/dr} ]
; ptb_uegf1 = ( ( (ptb_vtan0) ^ 2) / r) - ( (1/rho) * ptb_dpdr0)
  ptb_uegf1 = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; [R3]

; Arrays for LHS and RHS source terms after azimuthally averaging (4D)
  ptb_ueh0  = new((/num_l,dimsizes(radii),nts/),float)		  ; [L1]
  ptb_uev0  = new((/num_l,dimsizes(radii),nts/),float)  		  ; [R2]
  ptb_uegf0 = new((/num_l,dimsizes(radii),nts/),float)	          ; [R3]

;==================================================
; Arrays for horizontal diffusive tendency terms 
;==================================================

; Horizontal diffusivity (K_mh)
  k_mh0	    = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)

; Arrays for fractional expressions in diffusive tendency source terms below (4D)
; (1) du/dr, (2) dv/dy, (3) u/r, (4) r, (5) d/dr
  dudr0     = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; (1)
  dvdy0	    = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; (2)
  u_r0      = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; (3)
  r_0       = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; (4)
  ddr0	    = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; (5)

; Additional arrays for terms outside main loop
  du        = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ;
  dv	    = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ;
  dl        = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ;
  dr        = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ;
  r1        = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ;

  du0	    = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; Radial 'u' increment
  du1	    = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; Radial 'u'' increment
  dp1	    = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; Radial 'p'' increment

  dmp	    = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; Radial 'p' increment
  dmu	    = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; Radial 'u' increment 

; Radial stress tensor: RHS (5a)
; { (2 * K_mh) * (du/dr) }
; tau_rr1   = (2 * k_mh0) * (dudr0)
  tau_rr1   = new((/num_l,dimsizes(radii),nts/),float)		     ; [3D array]

; Azimuthal stress tensor: RHS (5b)
; { (2 * K_mh) * ( (1/r * dv/dy) + u/r ) }
; tau_yy1   = (2 * k_mh0) * ( ( (1 / r_0) * dvdy0) + u_r0 ) 
  tau_yy1   = new((/num_l,dimsizes(radii),nts/),float)		     ; [3D array]

; Total horizontal diffusive tendency of radial momentum 
; ( (1/r*{rho}) * d( r * {rho} * {tau_rr} )/dr) - tau_yy/r
; u_dh0	    = ((1/(r_0*mean_rho0)) * (ddr0*(r_0*mean_rho0*tau_rr1)))-((1/r_0)*(tau_yy1))
  u_dh0     = new((/num_l,dimsizes(radii),nts/),float)		     ; [3D array]

; Total vertical diffusive tendency of radial momentum 
; u_dz0	    = (1 /{rho0}) * ( (d/dz) * ({rho0} * {tau_rz}) )  
  u_dz0	    = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float) ; [4D array]

;========================================
; Define 3D arrays for mean quantities 
;========================================

; Mean tangential wind {v}
  mean_vtan0       = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)
  mean_vtan0!0     = "lev"
  mean_vtan0!1     = "rad"
  mean_vtan0!2     = "azi"
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

; Mean pressure {p}
  mean_pres0       = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)
  mean_pres0!0     = "lev"
  mean_pres0!1     = "rad"
  mean_pres0!2	   = "azi"
  mean_pres0!3     = "time"
  mean_pres0@description = "Azimuthally averaged pressure"
  mean_pres0@units = "hPa"

; Mean density {rho}
  mean_rho0        = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)
  mean_rho0!0      = "lev"
  mean_rho0!1      = "rad"
  mean_rho0!2	   = "azi"
  mean_rho0!3      = "time"
  mean_rho0@description = "Azimuthally averaged density"
  mean_rho0@units  = "kg m~S~-3~N~"

; Mean horizontal diffusivity {k_mh}
  mean_kmh0        = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)
  mean_kmh0!0      = "lev"
  mean_kmh0!1      = "rad"
  mean_kmh0!2      = "azi"
  mean_kmh0!3      = "time"
  mean_kmh0@description = "Azimuthally averaged horizontal diffusivity"
  mean_kmh0@units  = "m~S~2~N~ s~S~-1~N~"

; Mean Coriolis force {f0}
  mean_cor0        = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)
  mean_cor0!0      = "lev"
  mean_cor0!1      = "rad"
  mean_cor0!2      = "azi"
  mean_cor0!3      = "time"
  mean_cor0@description = "Azimuthally averaged Coriolis force"
  mean_cor0@units  = "s~S~-1~N~"

; Mean vertical diffusion (total)
  mean_difv0          = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)
  mean_difv0!0        = "lev"
  mean_difv0!1        = "rad"
  mean_difv0!2        = "azi"
  mean_difv0!3        = "time"
  mean_difv0@description = "Azimuthally averaged vertical diffusion (tot)"
  mean_difv0@units    = "m s~S~-1~N~"

; Radial derivative of mean radial wind [ d{u} / dr ]
  mean_dudr0          = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)
  mean_dudr0!0        = "lev"
  mean_dudr0!1        = "rad"
  mean_dudr0!2	      = "azi"
  mean_dudr0!3        = "time"
  mean_dudr0@description = "Azimuthally averaged 'du_dr'"

; Radial derivative of pressure [ d{p} / dr ] 
  mean_dpdr0          = new((/num_l,dimsizes(radii),dimsizes(thetas),nts/),float)
  mean_dpdr0!0        = "lev"
  mean_dpdr0!1        = "rad"
  mean_dpdr0!2	      = "azi"
  mean_dpdr0!3        = "time"
  mean_dpdr0@description = "Azimuthally averaged 'dp_dr'"

;=========================
; Start multiple loops
;=========================

    print_clock("Working on ensemble member: "+ens0(st))

  ; Set counter variable and size of arrays 
    llbox	= toint(dist*50)	   
    dsize	= (/llbox,llbox/)          
    d0		= dsize(0)
    d1		= dsize(1)

  ; Set time-dependent counter variable (21/01/2020)
  ; previously --> ct = ts0(st)
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

;===============================
; Start loop over input files 
;===============================

  ; Read in 'pc' stream [u,v,vort]
    c  = addfile(fili_c,"r")

  ; Read in 'pd' stream [w,p,t,theta,kmh,tau]
    d  = addfile(fili_d,"r")

  ; Create array to hold storm track positions using different methods
    centre_new = new((/7, num_l, nts, 2/),float)

  ; EDIT -- CHECK THIS IS FINE
    icp_arr    = new((/nts/),float)    

;==================================
; Get the variables we will need
;==================================

 ; Loop over times in file
   do it = ts0(st), tf0(st)

    print_clock("Looping between T+"+ts0(st)+" and T+"+tf0(st)+" ("+ens0(st)+", "+dat(st)+")")

    print("Working on time: "+time_str(it))
    title_arr(it) = "Valid at "+time_str(it)+" (T+"+it+")"

;=====================================================
; Read in variables from 'pc' stream (model levels)
;=====================================================

    t0    = lat_arr(cn,it-1,0)
    t1    = lat_arr(cn,it-1,1)
    n0    = lon_arr(cn,it-1,0)
    n1    = lon_arr(cn,it-1,1)
    print_clock("Grid extends from "+t0+" to "+t1+" deg N, and "+n0+" to "+n1+" deg E")

  ; Array for input into external function
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
    copy_VarCoords(u,vort)
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
    du_0   = grad_latlon_cfd(u,u&$dim1$,u&$dim2$,False,False)
    du_dy  = du_0[0]
    du_dx  = du_0[1]
    delete(du_0)

  ; Calculate horizontal derivatives of 'v' (dv_dx, dv_dy)                         
    dv_0   = grad_latlon_cfd(v,v&$dim1$,v&$dim2$,False,False)
    dv_dy  = dv_0[0]
    dv_dx  = dv_0[1]
    delete(dv_0)

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

    dx0      = lonD(1) - lonD(0) ; Grid spacing (longitude)
    dy0      = latD(1) - latD(0) ; Grid spacing (latitude)

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
     avo_plane              = avo(ilev,:,:)        ; Absolute vorticity
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
     r0 = 0.12

     print_clock("Calculating storm centre position...")
     centre_arr    = st_centre(ilev, it, cn, r0, centre, vort_plane, prs_plane, spd_plane, \
                               centre_new, dy0, dx0, "comp", ct, "pd", "off", "1h")

   ; Updated storm centre info
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

;=========================================================
; Switch to cylindrical coordinates (external function)
;=========================================================

   ; Call 'setup_cyl' 
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
     dvr_int  = cyl_arr[28]

     xpos     = cyl_arr[23]
     ypos     = cyl_arr[24]
     xcen     = cyl_arr[25]
     ycen     = cyl_arr[26]

     print("Storm centre = "+lat_max+" ; "+lon_max)

   ; Fill 'u_dz0' array with values from 'difv_tot' (vertical diffusive tendency)
     u_dz0(ilev,:,:,ct)   = dvr_int(:,:)    

;===================================================
; Calculate mean (azimuthally-averaged) variables 
;===================================================

     do irad = 0, dimsizes(radii)-1

        mean_vtan0(ilev,irad,:,ct)    = avg(vt_int(irad,:))      ; [ {v} ]
	mean_vrad0(ilev,irad,:,ct)    = avg(vr_int(irad,:))      ; [ {u} ]
	mean_vvel0(ilev,irad,:,ct)    = avg(w_int(irad,:))      ; [ {w} ]
	mean_pres0(ilev,irad,:,ct)    = avg(prs_int(irad,:))    ; [ {p} ]
	mean_difv0(ilev,irad,:,ct)    = avg(dvr_int(irad,:))   ; [ {diff_xz + diff_yz} ]
	mean_rho0(ilev,irad,:,ct)     = avg(rho_int(irad,:))    ; [ {rho} ]
	mean_kmh0(ilev,irad,:,ct)     = avg(kmh_int(irad,:))    ; [ {k_mh} ]
	mean_cor0(ilev,irad,:,ct)     = avg(f_int(irad,:))      ; [ {f} ]

     end do

   ; Test out the centered finite difference function 

   ; Calculate radial distance between all points in our array                        
;     drad       = dimsizes(radii)
;     xy_coords  = gc_latlon( ypos(0:drad-2,:),xpos(0:drad-2,:),\
;                             ypos(1:drad-1,:),xpos(1:drad-1,:), 2, 3)

   ; 29/01/2020. Function works with scalar argument, but produces 'inf' values with vector argument
;     du_dr  = center_finite_diff_n(vr_int(0:drad-2,:),5500.0,True,0,0)

     fcor = mean_cor0 * mean_vtan0

   ; Read 'v_int' etc into larger arrays
     u_int0(ilev,:,:,ct)    = (/u_int(:,:)/)
     v_int0(ilev,:,:,ct)    = (/v_int(:,:)/)
     w_int0(ilev,:,:,ct)    = (/w_int(:,:)/)
     prs_int0(ilev,:,:,ct)  = (/prs_int(:,:)/)
     f_int0(ilev,:,:,ct)    = (/f_int(:,:)/)
     kmh_int0(ilev,:,:,ct)  = (/kmh_int(:,:)/)
     rho_int0(ilev,:,:,ct)  = (/rho_int(:,:)/)
     difv_x0(ilev,:,:,ct)   = (/dvx_int(:,:)/)
     difv_y0(ilev,:,:,ct)   = (/dvy_int(:,:)/)
     difv_tot0(ilev,:,:,ct) = (/dvr_int(:,:)/)
     v_rad0(ilev,:,:,ct)    = (/vr_int(:,:)/)
     v_tan0(ilev,:,:,ct)    = (/vt_int(:,:)/)

   ; Tidy up
     delete([/u_int,v_int,w_int,vt_int,vr_int/])
     delete([/prs_int,rho_int,dvx_int,dvy_int,dvr_int,f_int,kmh_int/])
     delete([/u_plane,v_plane,w_plane,f_plane,vort_plane,avo_plane,prs_plane/])
     delete([/difv_xplane,difv_yplane,kmh_plane/])

     end do
   ; End model level loop (do ilev = 0, num_l-1)

    ct  = ct + 1
  ; Advance counter variable (time)

  ; Tidy up before next iteration (all variables from 'pc' and 'pd' streams)
    delete([/lonC,latC,lon1C,lat1C,u,v,vort,avo,f0,f1/])
    delete([/lonD,latD,lon1D,lat1D,str_xz,str_yz,kmh,t,p,w,z,rho,rho_a/])

   end do
 ; End time loop (do it = ts0(st), tf0(st) ) 

;====================================
; Calculate perturbation variables
;====================================

  print_clock("Finished looping over time intervals! Now calculating perturbation terms...")

  ptb_vtan0(:,:,:,:) = v_tan0(:,:,:,:)    - mean_vtan0(:,:,:,:)    ; [v']
  ptb_vvel0(:,:,:,:) = w_int0(:,:,:,:)    - mean_vvel0(:,:,:,:)    ; [w']
  ptb_vrad0(:,:,:,:) = v_rad0(:,:,:,:)    - mean_vrad0(:,:,:,:)    ; [u']
  ptb_pres0(:,:,:,:) = prs_int0(:,:,:,:)  - mean_pres0(:,:,:,:)    ; [p']

  k_mh0(:,:,:,:)     = kmh_int0(:,:,:,:)                           ; [kmh]
  rho0(:,:,:,:)	     = rho_int0(:,:,:,:)			   ; [rho]

;==================================================================
; Calculate all horizontal derivatives using built-in functions
;==================================================================

; Define arrays to hold each of these terms 
  mean_du_dr = new( (/num_l,dimsizes(radii),dimsizes(thetas),nts/), "float")
  ptb_du_dr  = new( (/num_l,dimsizes(radii),dimsizes(thetas),nts/), "float")
  mean_dp_dr = new( (/num_l,dimsizes(radii),dimsizes(thetas),nts/), "float")
  ptb_dp_dr  = new( (/num_l,dimsizes(radii),dimsizes(thetas),nts/), "float")
  du_dr      = new( (/num_l,dimsizes(radii),dimsizes(thetas),nts/), "float")
  dw_dr      = new( (/num_l,dimsizes(radii),dimsizes(thetas),nts/), "float")

; Loop over all times, vertical levels and grid points to calculate derivatives 
  do it = 0, nts-1
   print("Working on time: "+it)
   do il = 0, num_l-1
    do th0 = 0, dimsizes(thetas)-1

   ; Calculate radial distance between all points in our array
     do rad0 = 0, dimsizes(radii)-2

      xy_coords = gc_latlon( ypos_all(il,rad0,th0,it),\
                             xpos_all(il,rad0,th0,it),\
                             ypos_all(il,rad0+1,th0,it),\
                             xpos_all(il,rad0+1,th0,it), 2, 3)

    ; Arrays for which to calculate radial derivatives 
      mean_u_arr  = mean_vrad0(il,:,th0,it)
      ptb_u_arr   = ptb_vrad0(il,:,th0,it)
      u_arr       = v_rad0(il,:,th0,it)
      mean_p_arr  = mean_pres0(il,:,th0,it)
      ptb_p_arr   = ptb_pres0(il,:,th0,it)
      w_arr       = w_int0(il,:,th0,it)

    ; Calculate finite differences in radial direction 
      mean_du_dr(il,:,th0,it) = (center_finite_diff_n( mean_u_arr,xy_coords,True,0,0) )
      ptb_du_dr(il,:,th0,it)  = (center_finite_diff_n( ptb_u_arr,xy_coords,True,0,0) )
      du_dr(il,:,th0,it)      = (center_finite_diff_n( u_arr,xy_coords,False,0,0) )
      mean_dp_dr(il,:,th0,it) = (center_finite_diff_n( mean_p_arr,xy_coords,False,0,0) )
      ptb_dp_dr(il,:,th0,it)  = (center_finite_diff_n( ptb_p_arr,xy_coords,False,0,0) )
      dw_dr(il,:,th0,it)      = (center_finite_diff_n( w_arr,xy_coords,False,0,0) )

;    ; 28/01/2020 --> why does setting opt1=5500 get rid of the values of infinity?
;      if ( any (mean_du_dr(il,:,th0,it) .gt. mean_du_dr@_FillValue) ) then 
;       print("Something wrong with azimuthal region: "+th0+"...")
;       print(xy_coords)
;       print(mean_vrad0(il,:,th0,it))
;       print(mean_du_dr(il,:,th0,it))
;      end if 

     end do ; End loop over radial bands 

    end do 
   end do 
  end do 
  delete([/mean_u_arr,ptb_u_arr,u_arr,mean_p_arr,ptb_p_arr,w_arr,xy_coords/])

; Also define arrays to hold azimuthal derivatives 
  dv_dy0     = new( (/num_l,dimsizes(radii),dimsizes(thetas),nts/), "float")
  ptb_du_dy  = new( (/num_l,dimsizes(radii),dimsizes(thetas),nts/), "float")

  do it = 0, nts-1
   print("Working on time: "+it)
   do il = 0, num_l-1
    do rad0 = 0, dimsizes(radii)-1

    ; Calculate azimuthal distance between all points in our array 
      do dth     = 0, dimsizes(thetas)-2

       xy_coords  = gc_latlon( ypos_all(il,rad0,dth,it),\
                               xpos_all(il,rad0,dth,it),\
                               ypos_all(il,rad0,dth+1,it),\
                               xpos_all(il,rad0,dth+1,it), 2, 3)

     ; Arrays to calculate derivatives for 
       ptb_u_arr   = ptb_vrad0(il,rad0,:,it)
       v_arr       = v_tan0(il,rad0,:,it)
    
     ; Avoid dividing by zero at the storm centre 
       xy_coords   = where( xy_coords .ne. 0, xy_coords, 400)

     ; Calculate finite differences in radial direction
       ptb_du_dy(il,rad0,:,it) = (center_finite_diff( ptb_u_arr,xy_coords,True,0) )
       dv_dy0(il,rad0,:,it)    = (center_finite_diff_n( v_arr,xy_coords,True,0,0) )

       if ( any (ptb_du_dy(il,rad0,:,it) .gt. ptb_du_dy@_FillValue) ) then
        print("Something wrong with azimuthal angle: "+thetas(dth)+"...")
        print(xy_coords)
        print(ptb_du_dy(il,rad0,:,it))
       end if

     end do
   ; End loop over azimuthal angles 

    end do
   end do
  end do
  delete([/ptb_u_arr,v_arr,xy_coords/])

; Get correct units for pressure tendency terms 
  mean_dp_dr = mean_dp_dr * 100
  ptb_dp_dr  = ptb_dp_dr * 100

;====================================================================
; Calculate distance between all grid points (and centre of grid)
;====================================================================

  drad = dimsizes(radii)

; Calculate distance between each grid point and the centre of the grid 
  r1(:,0:drad-2,:,:)  = gc_latlon( ycen_all(:,1:drad-1,:,:),xcen_all(:,1:drad-1,:,:),\
                                   ypos_all(:,1:drad-1,:,:),xpos_all(:,1:drad-1,:,:), 2, 3)

; Replace zeros with positive values  
  r1 = where( r1 .ne. 0, r1, 1000)

; Calculate distance between each neighbouring grid point 
  dr(:,0:drad-2,:,:)  = gc_latlon( ypos_all(:,0:drad-2,:,:),xpos_all(:,0:drad-2,:,:),\
                                   ypos_all(:,1:drad-1,:,:),xpos_all(:,1:drad-1,:,:), 2, 3)

; Replace zeros with positive values
  dr = where( dr .ne. 0, dr, 1000)

;======================================================================
; Calculate eddy terms involving radial derivatives [du'/dr, dp'/dr]
;======================================================================

; Calculate [du'/dr]
  ptb_dudr0(:,:,:,:)    = ptb_du_dr

; [du'dy]
  ptb_dudy0(:,:,:,:)    = ptb_du_dy

; Calculate [du/dr -- for horizontal diffusive tendency]
  dudr0(:,:,:,:)        = du_dr

; Calculate [dp'/dr]
  ptb_dpdr0(:,:,:,:)    = ptb_dp_dr

; Calculate [u/r -- for horizontal diffusive tendency]
  u_r0(:,:,:,:)         = v_rad0(:,:,:,:) / r1

; Calculate [r -- for horizontal diffusive tendency]
  r_0(:,:,:,:)          = r1

; Calculate [d/dr -- for horizontal diffusive tendency]
; EDIT --> not 100% sure of the meaning of this expression
  ddr0(:,:,:,:)         = 1 / dr

; Calculate [v'/r]
  ptb_vt_r0(:,:,:,:)    = ptb_vtan0(:,:,:,:) / r1

; Calculate [v'^2/r]
  ptb_vttr0(:,:,:,:)    = (ptb_vtan0(:,:,:,:) ^ 2) / r1

;=========================================================================
; Calculate mean terms involving radial derivatives [ d{u}/dr, d{p}/dr ]
;=========================================================================

    ; Calculate the radial distance ('dr') between all points in array
    ; Arguments 1 and 2 represent the lat/lon of the first point, respectively
    ; Arguments 3 and 4 represent the lat/lon of the second point, respectively
    ; Argument 5 ('2') tells the script to avoid extra interpolated points
    ; Argument 6 ('3') returns the distance in metres

    ; 'r1' and 'dr' have already been defined

;=============================================
; Calculate absolute angular momentum (AAM)
;=============================================

    ; Calculate increments 
      vm   = mean_vtan0(:,:,:,:)
      vmm  = vm ^ 2
      fm   = mean_cor0(:,:,:,:) * (r1 ^ 2)

    ; Calculate absolute angular momentum components (AAM)
      aam0 = vm * r1
      aam1 = 0.5 * fm

    ; Finally, calculate absolute angular momentum [ rv + 0.5 * (fr ^ 2) ]
      mean_aam0(:,:,:,:)  = aam0 + aam1

;====================================================
; Calculate additional components including radius
;====================================================

    ; Calculate [ d{u}/dr ]
      mean_dudr0          = mean_du_dr

    ; Calculate [ d{p}/dr ]
      mean_dpdr0          = mean_dp_dr

    ; Calculate [ {v}^2/r ]
      mean_vttr0          = ( mean_vtan0(:,:,:,:) ^ 2) / r1

;=======================================================
; Remove unncessary 3rd dimension from mean variables
;=======================================================

      mean_vtan1 = mean_vtan0(:,:,0,:)
      mean_vrad1 = mean_vrad0(:,:,0,:)
      mean_vvel1 = mean_vvel0(:,:,0,:)
      mean_rho1  = mean_rho0(:,:,0,:)
      mean_aam1  = mean_aam0(:,:,0,:)
      mean_cor1	 = mean_cor0(:,:,0,:)

      mean_dudr1 = mean_dudr0(:,:,0,:)
      mean_dpdr1 = mean_dpdr0(:,:,0,:)
      mean_vttr1 = mean_vttr0(:,:,0,:)

;================================================
; Calculate vertical derivatives of mean terms 
;================================================

    ; Redefine model level array for vertical derivative calculations
      hybC   = c->hybrid_ht(:)         ; 63 model ('theta') levels
      mlevs  = dimsizes(hybC)	       ; Size of model level array
      hybC@units = "m"		       ; Model level units 


    ; Loop over times 
      ft    = 0
      do it = ts0(st), tf0(st)

      ; Loop over model levels
	do ilev = 0, num_l-1
	
        mp1 = min((/ilev+1,num_l-1/))
        mm1 = max((/ilev-1,0/))
 
        ; Loop over radial circles 
          do irad = 0, dimsizes(radii)-1

          ;=====================================================
	  ; Calculate vertical derivative of mean radial wind
	  ;=====================================================

          ; [ d{u}/dz ]
	    mean_dudz1(ilev,irad,ft)    = (mean_vrad1(mp1,irad,ft) - \
                                           mean_vrad1(mm1,irad,ft)) / \
                                          (hybC(mp1) - hybC(mm1))

          end do
        ; End loop over radial circles (do irad = 0, ...)

        end do
      ; End loop over model levels (do ilev = 0, num_l-1) 

      ; Advance counter variable
        ft = ft + 1 

       end do
     ; End loop over times (do it = ts0(st), tf0(st) )

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

           ;=====================================================
           ; Calculate vertical derivative of eddy radial wind
           ;=====================================================

           ; [ du'/dz ]
             ptb_dudz0(ilev,irad,th0,ft) = (ptb_vrad0(mp1,irad,th0,ft) - \
                                            ptb_vrad0(mm1,irad,th0,ft)) / \
                                           (hybC(mp1) - hybC(mm1))

           end do
         ; End loop over azimuth angles (do th0 = 0, ...)
 
          end do
        ; End loop over radial circles (do irad = 0, ...)

        end do
      ; End loop over model levels (do ilev = 0, num_l-1) 

      ; Advance counter variable  
        ft = ft + 1

       end do
     ; End loop over times (do it = ts0, tf0)

;==========================================================
; Finally, calculate source terms and azimuthal averages 
;==========================================================

  print_clock("Calculating source terms and azimuthal averages...")

;===================================================================
; Mean radial advection of radial momentum per unit mass: LHS (2)
;===================================================================
  lhs2_st = mean_vrad1 * mean_dudr1

;=============================================================
; Mean horizontal advection of eddy radial momentum: LHS (3)
;=============================================================
  ptb_ueh1 = (ptb_vrad0 * ptb_dudr0) + (ptb_vt_r0 * ptb_dudy0)

; Calculate the azimuthal average of 'ptb_ueh1' and input into 'lhs3_st
  lhs3_st = dim_avg_n(ptb_ueh1(:,:,:,:), 2)

;=========================================================================
; Mean vertical advection of mean radial momentum per unit mass: RHS (1)
;=========================================================================
  rhs1_st = -(mean_vvel1 * mean_dudz1)

;=========================================================================
; Eddy vertical advection of eddy radial momentum per unit mass: RHS (2)
;=========================================================================
  ptb_uev1 = -(ptb_vvel0 * ptb_dudz0)

; Calculate the azimuthal average of 'ptb_uev1' and input into 'rhs2_st'
  rhs2_st  = dim_avg_n(ptb_uev1(:,:,:,:),2)

;==============================================
; Mean agradient force per unit mass: RHS (3)
;==============================================

; Account for elements where 'mean_rho1' is equal to zero (should be none)
  mean_rho1 = where(mean_rho1 .ne. 0, mean_rho1, mean_rho1@_FillValue)

  rhs3_st   = mean_vttr1 + (mean_cor1 * mean_vtan1) - ( (1 / mean_rho1) * mean_dpdr1 )

; Split agradient force term into components
  mean_cor  = mean_cor1 * mean_vtan1
  mean_cen  = mean_vttr1
  mean_pgf  = - ( (1 / mean_rho1) * mean_dpdr1 )

;==============================================
; Eddy agradient force per unit mass: RHS (4)
;==============================================

; Account for elements where 'rho0' is equal to zero
  rho0  = where(rho0 .ne. 0, rho0, rho0@_FillValue)
     
; Now calculate eddy agradient force 
  ptb_uegf1 = ptb_vttr0 - ( (1 / rho0) * ptb_dpdr0 )

; Calculate the azimuthal average of 'ptb_uegf1' and input into 'rhs4_st'
  rhs4_st   = dim_avg_n(ptb_uegf1(:,:,:,:),2)

;===================================================================================
; Combined mean radial and vertical diffusive tendency of radial momentum: RHS (5)
;===================================================================================

; Stress tensor (radial component) 
  tau_rr1   = dim_avg_n( ((2 * k_mh0) * (du_dr)), 2)
     
; Stress tensor (azimuthal component)
  tau_yy1	 = dim_avg_n( ((2 * k_mh0) * ( ( (1 / r_0) * dv_dy0) + u_r0 )), 2)

;=============================================================================
; Use stress tensors to calculate horizontal diffusive tendency source term
;=============================================================================

; Azimuthally average the 'r' and 'd/dr' expressions
  r_1   = dim_avg_n(r_0(:,:,:,:),2)
  ddr1  = dim_avg_n(ddr0(:,:,:,:),2)

  u_dh0 = ( ( 1/(r_1 * mean_rho1) ) * (ddr1 * (r_1 * mean_rho1 * tau_rr1)) ) - \
          ( (1/r_1) * (tau_yy1) )

  rhs5a_st = u_dh0	       

; Vertical diffusive tendency source term has already been calculated (L827)
; Now calculate the azimuthal average using 'dim_avg_n'
  rhs5b_st = dim_avg_n( u_dz0(:,:,:,:),2 )

;=================================================
; Tidy up and produce final arrays for plotting 
;=================================================

    ; Calculating the terms from Montgomery et al. (2018) - their Eq.(8)
    ; LHS (1) : d{u}/dt          - Local tendency of the mean radial velocity
    ; LHS (2) : {u} * d{u}/dr    - Mean radial advection of radial momentum
    ; LHS (3) : {(u' * du'/dr) + (v'/r * du'/dy)}  - Mean horiz. adv. of eddy radial momentum

    ; RHS (1) : -{w} * d{u}/dz   - Mean vertical advection of mean radial momentum
    ; RHS (2) : -{w' * d{u'}/dz} - Eddy vertical advection of eddy radial momentum
    ; RHS (3) : {v}^2/r + f{v} - 1/{rho} * d{p}/dr - Mean agradient force
    ; RHS (4) : {v'^2/r - 1/rho * dp'/dr}          - Eddy agradient force
    ; RHS (5) : {F_r}            - Mean radial/vertical diffusive tendency of radial momentum

    ; Define new arrays for terms in momentum budget 
      lhs1A    = new((/num_l,dimsizes(radii),nts/),float)   ; Local 'u' tendency
      lhs2A    = new((/num_l,dimsizes(radii),nts/),float)   ; Mean r-adv of rad mom
      lhs3A    = new((/num_l,dimsizes(radii),nts/),float)   ; Mean h-adv of eddy mom
      lhsTA    = new((/num_l,dimsizes(radii),nts/),float)   ; Sum of LHS terms

      rhs1A    = new((/num_l,dimsizes(radii),nts/),float)   ; Mean v-adv of mean mom
      rhs2A    = new((/num_l,dimsizes(radii),nts/),float)   ; Eddy v-adv of eddy mom
      rhs3A    = new((/num_l,dimsizes(radii),nts/),float)   ; Mean agradient force 
      rhs4A    = new((/num_l,dimsizes(radii),nts/),float)   ; Eddy agradient force
      rhs5aA   = new((/num_l,dimsizes(radii),nts/),float)   ; Mean horiz diffusion
      rhs5bA   = new((/num_l,dimsizes(radii),nts/),float)   ; Mean vert diffusion
      rhsTA    = new((/num_l,dimsizes(radii),nts/),float)   ; Sum of RHS terms

    ; Components of mean agradient force terms
      mean_cor = new((/num_l,dimsizes(radii),nts/),float)   ; Coriolis force 
      mean_cen = new((/num_l,dimsizes(radii),nts/),float)   ; Centrifugal force
      mean_pgf = new((/num_l,dimsizes(radii),nts/),float)   ; Pressure gradient force

    ;=================================================
    ; Calculate the mean radial wind tendency (LHS)
    ;=================================================

    ; Loop over total number of times in arrays
      do it = 0, nts-1

        cp1 = min((/it+1,nts-1/))
        cm1 = max((/it-1,0/))
 
        lhs1A(:,:,it) = (mean_vrad1(:,:,cp1) - mean_vrad1(:,:,cm1)) / 2

      end do
    ; End loop over times (do it = 0, nts-1)

    ; Add metadata and make units consistent for all terms
      lhs1A@units = "m s~S~-1~N~ h~S~-1~N~"

      lhs2A = (/lhs2_st/)
      lhs2A = lhs2A * 3600 
      lhs2A@units = "m s~S~-1~N~ h~S~-1~N~"

      lhs3A = (/lhs3_st/)
      lhs3A = lhs3A * 3600
      lhs3A@units = "m s~S~-1~N~ h~S~-1~N~"

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

    ; Also add metadata for the components of the agradient force term
      mean_cor       = mean_cor * 3600
      mean_cen       = mean_cen * 3600
      mean_pgf       = mean_pgf * 3600
      mean_cor@units = "m s~S~-1~N~ h~S~-1~N~"
      mean_cen@units = "m s~S~-1~N~ h~S~-1~N~"
      mean_pgf@units = "m s~S~-1~N~ h~S~-1~N~"

      mean_cen!0     = "lev"
      mean_cen&lev   = hyb_plot
      mean_cen!1     = "rad"
      mean_cen&rad   = radii
      mean_cen!2     = "time"

      mean_cor!0     = "lev"
      mean_cor&lev   = hyb_plot
      mean_cor!1     = "rad"
      mean_cor&rad   = radii
      mean_cor!2     = "time"

      mean_pgf!0     = "lev"
      mean_pgf&lev   = hyb_plot
      mean_pgf!1     = "rad"
      mean_pgf&rad   = radii
      mean_pgf!2     = "time"

    ; Add metadata for all the other important arrays (instantaneous)
      lhs1A!0    = "lev"
      lhs1A&lev  = hyb_plot
      lhs1A!1    = "rad"
      lhs1A&rad  = radii
      lhs1A!2    = "time"

      lhs2A!0    = "lev"
      lhs2A&lev  = hyb_plot
      lhs2A!1    = "rad"
      lhs2A&rad  = radii
      lhs2A!2    = "time"

      lhs3A!0    = "lev"
      lhs3A&lev  = hyb_plot
      lhs3A!1    = "rad"
      lhs3A&rad  = radii
      lhs3A!2    = "time"

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

      lhsTA!0    = "lev"
      lhsTA&lev  = hyb_plot
      lhsTA!1    = "rad"
      lhsTA&rad  = radii
      lhsTA!2    = "time"

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
      mean_aam	      = mean_aam / (10 ^ 6)
      mean_aam@units  = "10~S~6~N~ m~S~2~N~ s~S~-1~N~"

    ;==================================
    ; Reorder arrays before plotting
    ;==================================
       
      lhs1  	  = lhs1A(lev|:,time|:,rad|:)
      lhs2  	  = lhs2A(lev|:,time|:,rad|:)
      lhs3  	  = lhs3A(lev|:,time|:,rad|:)

      rhs1  	  = rhs1A(lev|:,time|:,rad|:)
      rhs2  	  = rhs2A(lev|:,time|:,rad|:)
      rhs3  	  = rhs3A(lev|:,time|:,rad|:)
      rhs4  	  = rhs4A(lev|:,time|:,rad|:)
      rhs5a 	  = rhs5aA(lev|:,time|:,rad|:)
      rhs5b 	  = rhs5bA(lev|:,time|:,rad|:)

    ; Components of the mean agradient force term
      cor_plot    = mean_cor(lev|:,time|:,rad|:)
      cen_plot    = mean_cen(lev|:,time|:,rad|:)
      pgf_plot    = mean_pgf(lev|:,time|:,rad|:)

;=========================================
; Finally, smooth data before plotting
;=========================================

    ; 21/01/2020 --> add description for LHS terms (following H18)
      lhs1       = smth9_Wrap(lhs1, 0.5, 0.5, True)
      lhs2       = smth9_Wrap(lhs2, 0.5, 0.5, True)
      lhs3  	 = smth9_Wrap(lhs3, 0.5, 0.5, True)

    ; Mean hadv, mean vadv, eddy hadv, eddy vadv, h-diff, v-diff 
      rhs1  	 = smth9_Wrap(rhs1, 0.5, 0.5, True)
      rhs2  	 = smth9_Wrap(rhs2, 0.5, 0.5, True)
      rhs3  	 = smth9_Wrap(rhs3, 0.5, 0.5, True)
      rhs4  	 = smth9_Wrap(rhs4, 0.5, 0.5, True)
      rhs5a 	 = smth9_Wrap(rhs5a, 0.5, 0.5, True)
      rhs5b 	 = smth9_Wrap(rhs5b, 0.5, 0.5, True)

    ; Calculate total RHS tendency
      rhsT       = rhs1 + rhs2 + rhs3 + rhs4 + rhs5a + rhs5b
      rhsT!0     = "lev"
      rhsT&lev   = hyb_plot
      rhsT!1     = "time"
      rhsT!2     = "rad"
      rhsT&rad   = radii
      rhsT@units = "m s~S~-1~N~ h~S~-1~N~"

    ; Calculate total LHS tendency
      lhsT       = lhs1 + lhs2 + lhs3
      lhsT!0     = "lev"
      lhsT&lev   = hyb_plot
      lhsT!1     = "time"
      lhsT!2     = "rad"
      lhsT&rad   = radii
      lhsT@units = "m s~S~-1~N~ h~S~-1~N~"

    ; Calculate sum of all mean terms (vadv + hadv + diff)
      mean       = lhs2 + rhs1 + rhs5a + rhs5b
      mean!0     = "lev"
      mean&lev   = hyb_plot
      mean!1     = "time"
      mean!2     = "rad"
      mean&rad   = radii
      mean@units = "m s~S~-1~N~ h~S~-1~N~"

    ; Calculate sum of all eddy terms (vadv + hadv)
      eddy       = rhs2 + lhs3
      eddy!0     = "lev"
      eddy&lev   = hyb_plot
      eddy!1     = "time"
      eddy!2     = "rad"
      eddy&rad   = radii
      eddy@units = "m s~S~-1~N~ h~S~-1~N~"

    ; Fill 'plot_arr' with values from time averages
      if (agrd .eq. 1) then 
       mbgt_arr0 = (/mean_cor, mean_cen, mean_pgf, rhs3/)
      else
       mbgt_arr0 = (/lhs1, lhs2, lhs3, rhs1, rhs2, \
                     rhs3, rhs4, rhs5a, rhs5b, mean_vrad, \
                     mean_vtan, mean_vvel, mean_aam/)
      end if 

    ; Fill composite array with values from terms in 'mbgt_arr'
      do mb = 0, ndiag-1 
       plot_arr0(:,:,:,mb,st) = mbgt_arr0(mb,:,:,:)
      end do 
      delete(mbgt_arr0)

   ; Tidy up before next iteration
     delete([/utc_date,year,month,day,hour,minute,second,date_str,time_str,time_arr,title_arr/])
     delete([/numINFO_c,numINFO_d,numINPUT_c,time_c,numINPUT_d,time_d/])
     delete([/storm_rel,u_cyc,v_cyc,vel_cyc,lat_arr,lon_arr,centre/])

   end do
 ; End loop over sets of dates (do st = 0, dimsizes(dat)-1)

   print_clock("Finished looping over input data! Now calculating averages...")

;=================================================================================         
; Before plotting, average over all chosen time periods (unless analysing ERCs)            
;=================================================================================         

; If we're analysing ERCs rather than VCs, write out data at each time individually        
; 15/01/2020 --> add this step for ERCs specifically                                       
  if (typ .eq. "erc1" .or. typ .eq. "erc2" .or. typ .eq. "erc3" .or. typ .eq. "test") then

 ; Extract individual forcing terms from 'plot_arr' 
   if (agrd .eq. 0) then 
   lhs1_c0  = plot_arr0(:,:,:,0,:)
   lhs2_c0  = plot_arr0(:,:,:,1,:)
   lhs3_c0  = plot_arr0(:,:,:,2,:)

   rhs1_c0  = plot_arr0(:,:,:,3,:)
   rhs2_c0  = plot_arr0(:,:,:,4,:)
   rhs3_c0  = plot_arr0(:,:,:,5,:)
   rhs4_c0  = plot_arr0(:,:,:,6,:)
   rhs5a_c0 = plot_arr0(:,:,:,7,:)
   rhs5b_c0 = plot_arr0(:,:,:,8,:)

   vrad_c0  = plot_arr0(:,:,:,9,:)
   vtan_c0  = plot_arr0(:,:,:,10,:)
   vvel_c0  = plot_arr0(:,:,:,11,:)
   aam_c0   = plot_arr0(:,:,:,12,:)
 
   lhs1_c0@name  = "Local tendency of mean radial wind"
   lhs2_c0@name  = "Mean radial advection of radial momentum"
   lhs3_c0@name  = "Eddy horizontal advection of eddy radial momentum"

   rhs1_c0@name  = "Mean vertical advection of mean radial momentum"
   rhs2_c0@name  = "Eddy vertical advection of eddy radial momentum"
   rhs3_c0@name  = "Mean agradient force"
   rhs4_c0@name  = "Eddy agradient force"
   rhs5a_c0@name = "Radial diffusive tendency of radial momentum"
   rhs5b_c0@name = "Vertical diffusive tendency of radial momentum"

   vrad_c0@name  = "Radial wind"
   vtan_c0@name  = "Tangential wind"
   vvel_c0@name  = "Vertical velocity"
   aam_c0@name   = "Absolute angular momentum"

   else
    cor_plot0  = plot_arr0(:,:,:,0,:)
    cen_plot0  = plot_arr0(:,:,:,1,:)
    pgf_plot0  = plot_arr0(:,:,:,2,:)
    agr_plot0  = plot_arr0(:,:,:,3,:)
   end if 

;===============================================
; Calculate additional terms before averaging
;===============================================

    if (agrd .eq. 0) then 

  ; Total LHS
    lhsT_c0       = lhs1_c0 + lhs2_c0 + lhs3_c0
    copy_VarCoords(lhs1_c0,lhsT_c0)
    copy_VarMeta(lhs1_c0,lhsT_c0)
    lhsT_c0@name  = "Total LHS"

  ; Total RHS (21/01/2020 --> ignore horizontal component, which is tiny)
    rhsT_c0       = rhs1_c0 + rhs2_c0 + rhs3_c0 + rhs4_c0 + rhs5b_c0 ; + rhs5a_c0
    copy_VarCoords(rhs1_c0,rhsT_c0)
    copy_VarMeta(rhs1_c0,rhsT_c0)
    rhsT_c0@name  = "Total RHS"

  ; All mean terms
    mean_c0       = lhs2_c0 + rhs1_c0 + rhs5a_c0 + rhs5b_c0
    copy_VarCoords(rhs1_c0,mean_c0)
    copy_VarMeta(rhs1_c0,mean_c0)
    mean_c0@name  = "Sum of mean terms"

  ; All eddy terms
    eddy_c0       = rhs2_c0 + lhs3_c0
    copy_VarCoords(rhs2_c0,eddy_c0)
    copy_VarMeta(rhs2_c0,eddy_c0)
    eddy_c0@name  = "Sum of eddy terms"

  ; Radial vorticity flux (eddy + mean)
    vflux_c0      = lhs2_c0 + lhs3_c0
    copy_VarCoords(rhs1_c0,vflux_c0)
    copy_VarMeta(rhs1_c0,vflux_c0)
    vflux_c0@name = "Total radial vorticity flux"

  ; Vertical advection (eddy + mean)
    vadv_c0       = rhs1_c0 + rhs2_c0
    copy_VarCoords(rhs1_c0,vadv_c0)
    copy_VarMeta(rhs1_c0,vadv_c0)
    vflux_c0@name = "Total vertical advection"

  ; Agradient force (eddy + mean)
    agrd_c0       = rhs3_c0 + rhs4_c0
    copy_VarCoords(rhs1_c0,agrd_c0)
    copy_VarMeta(rhs1_c0,agrd_c0)
    agrd_c0@name  = "Total agradient force"

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

    if (agrd .eq. 0) then 

    ; Name all variables we want to write
      opt_arr = (/"lhs1", "lhs2", "lhs3", "lhsT", "agrd", \
                  "rhs1", "rhs2", "rhs3", "rhs4", "rhs5a", "rhs5b", "rhsT", \
                  "mean", "eddy", "vflux", "vadv", "vrad", "vtan", "vvel", "aam"/)

    ; Array of all composite terms                                                            
      mbgt_arr = (/lhs1_c0, lhs2_c0, lhs3_c0, lhsT_c0, agrd_c0, \
                   rhs1_c0, rhs2_c0, rhs3_c0, rhs4_c0, rhs5a_c0, rhs5b_c0, rhsT_c0, \
                   mean_c0, eddy_c0, vflux_c0, vadv_c0, vrad_c0, vtan_c0, vvel_c0, aam_c0/)

    else
    ; Name all variables we want to write
      opt_arr = (/"cor", "cen", "pgf", "agrd"/)
    ; Array of all composite terms
      mbgt_arr = (/cor_plot0, cen_plot0, pgf_plot0, agr_plot/)
    end if 

  ; Call procedure and write time-averaged                                                  
    print_clock("Writing 2D arrays to matrices for further analysis...")
    do mb = 0, dimsizes(opt_arr)-1

     print("Working on "+opt_arr(mb) )

   ; Set file path, etc                                                                     
     t_dir  = "./rad_mbgt_"
     t_out  = opt_arr(mb)+"_sc"+sc+"_group"+grp+"_"+w0+"w_"+hy3+"km_"+typ0+"_T"+it+".txt"

     t_path = t_dir + t_out
     print_clock("Writing to file: "+t_path)
     opts_mat@fout = t_path

   ; Check for existing files and then write
   ; 21/01/2020 --> run script for one extra radial circle to avoid missing AAM values  
     system("/bin/rm -f " + t_path)
     write_matrix( mbgt_arr(mb,:,ct,0:nr-2,0), fmt, opts_mat)

    end do
  ; End loop over diagnostics                                                               

    ct = ct + 1
  ; Advance counter variable (time)                                                         

   end do
 ; End loop over times                                                                      

 ; Tidy up if necessary
   delete(mbgt_arr)

  end if
; End IF statement on writing out data at each time interval (for ERC analysis)

;==================================================================
; Average over time dimension to calculate composite diagnostics
;================================================================== 

; START EDITING FROM HERE (03/03/2020)

  lhs1_c1  = dim_avg_n_Wrap(lhs1_c0(:,:,:,:),1)
  lhs2_c1  = dim_avg_n_Wrap(lhs2_c0(:,:,:,:),1)
  lhs3_c1  = dim_avg_n_Wrap(lhs3_c0(:,:,:,:),1)

  rhs1_c1  = dim_avg_n_Wrap(rhs1_c0(:,:,:,:),1)
  rhs2_c1  = dim_avg_n_Wrap(rhs2_c0(:,:,:,:),1)
  rhs3_c1  = dim_avg_n_Wrap(rhs3_c0(:,:,:,:),1)
  rhs4_c1  = dim_avg_n_Wrap(rhs4_c0(:,:,:,:),1)
  rhs5a_c1 = dim_avg_n_Wrap(rhs5a_c0(:,:,:,:),1)
  rhs5b_c1 = dim_avg_n_Wrap(rhs5b_c0(:,:,:,:),1)

  lhsT_c1  = dim_avg_n_Wrap(lhsT_c0(:,:,:,:),1)
  rhsT_c1  = dim_avg_n_Wrap(rhsT_c0(:,:,:,:),1)
  mean_c1  = dim_avg_n_Wrap(mean_c0(:,:,:,:),1)
  eddy_c1  = dim_avg_n_Wrap(eddy_c0(:,:,:,:),1)
  
  vflux_c1 = dim_avg_n_Wrap(vflux_c0(:,:,:,:),1)
  vadv_c1  = dim_avg_n_Wrap(vadv_c0(:,:,:,:),1)
  agrd_c1  = dim_avg_n_Wrap(agrd_c0(:,:,:,:),1)

  vrad_c1  = dim_avg_n_Wrap(vrad_c0(:,:,:,:),1)
  vtan_c1  = dim_avg_n_Wrap(vtan_c0(:,:,:,:),1)
  vvel_c1  = dim_avg_n_Wrap(vvel_c0(:,:,:,:),1)
  aam_c1   = dim_avg_n_Wrap(aam_c0(:,:,:,:),1)

; Now also average over all simulations to create 2D arrays 
  lhs1_c   = dim_avg_n_Wrap(lhs1_c1, 2)
  lhs2_c   = dim_avg_n_Wrap(lhs2_c1, 2)
  lhs3_c   = dim_avg_n_Wrap(lhs3_c1, 2)

  rhs1_c   = dim_avg_n_Wrap(rhs1_c1, 2)
  rhs2_c   = dim_avg_n_Wrap(rhs2_c1, 2)
  rhs3_c   = dim_avg_n_Wrap(rhs3_c1, 2)
  rhs4_c   = dim_avg_n_Wrap(rhs4_c1, 2)
  rhs5a_c  = dim_avg_n_Wrap(rhs5a_c1, 2)
  rhs5b_c  = dim_avg_n_Wrap(rhs5b_c1, 2)

  lhsT_c   = dim_avg_n_Wrap(lhsT_c1, 2)
  rhsT_c   = dim_avg_n_Wrap(rhsT_c1, 2)
  mean_c   = dim_avg_n_Wrap(mean_c1, 2)
  eddy_c   = dim_avg_n_Wrap(eddy_c1, 2)

  vflux_c  = dim_avg_n_Wrap(vflux_c1, 2)
  vadv_c   = dim_avg_n_Wrap(vadv_c1, 2)
  agrd_c   = dim_avg_n_Wrap(agrd_c1, 2)

  vrad_c   = dim_avg_n_Wrap(vrad_c1, 2)
  vtan_c   = dim_avg_n_Wrap(vtan_c1, 2)
  vvel_c   = dim_avg_n_Wrap(vvel_c1, 2)
  aam_c    = dim_avg_n_Wrap(aam_c1, 2)

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
  opt_arr  = (/"lhs1", "lhs2", "lhs3", "lhsT", "agrd", \
               "rhs1", "rhs2", "rhs3", "rhs4", "rhs5a", "rhs5b", "rhsT", \
               "mean", "eddy", "vflux", "vadv", "vrad", "vtan", "vvel", "aam"/)

; Array of all composite terms
  mbgt_arr = (/lhs1_c, lhs2_c, lhs3_c, lhsT_c, agrd_c, \
               rhs1_c, rhs2_c, rhs3_c, rhs4_c, rhs5a_c, rhs5b_c, rhsT_c, \
               mean_c, eddy_c, vflux_c, vadv_c, vrad_c, vtan_c, vvel_c, aam_c/)

; Call procedure and write time-averaged terms to text files 
  if (typ .ne. "erc1" .and. typ .ne. "erc2" .and. typ .ne. "erc3" .and. typ .ne. "test") then 
  print_clock("Writing 2D arrays to matrices for further analysis...")
  do mb = 0, dimsizes(opt_arr)-1

    print("Working on "+opt_arr(mb) ) 

  ; Set file path, etc                                                                       
    t_dir  = "./rad_mbgt_"
    t_out  = opt_arr(mb)+"_sc"+sc+"_group"+grp+"_"+w0+"w_"+hy3+"km_"+typ0+".txt"

    t_path = t_dir + t_out
    print_clock("Writing to file: "+t_path)
    opts_mat@fout = t_path

  ; Check for existing files and then write                                                  
    system("/bin/rm -f " + t_path)
    write_matrix( mbgt_arr(mb,:,0:nr-2), fmt, opts_mat)

   end do
 ; End loop over diagnostics                                                                 
   end if 

;===============================
; Add metadata before ploting 
;===============================

  lhs1_c!0       = "lev"
  lhs1_c&lev     = hyb_plot
  lhs1_c!1	 = "rad"
  lhs1_c&rad	 = radii

  lhs2_c!0       = "lev"
  lhs2_c&lev     = hyb_plot
  lhs2_c!1	 = "rad"
  lhs2_c&rad	 = radii

  lhs3_c!0       = "lev"
  lhs3_c&lev     = hyb_plot
  lhs3_c!1	 = "rad"
  lhs3_c&rad	 = radii

  rhs1_c!0       = "lev"
  rhs1_c&lev     = hyb_plot
  rhs1_c!1	 = "rad"
  rhs1_c&rad	 = radii

  rhs2_c!0       = "lev"
  rhs2_c&lev     = hyb_plot
  rhs2_c!1	 = "rad"
  rhs2_c&rad	 = radii

  rhs3_c!0       = "lev"
  rhs3_c&lev     = hyb_plot
  rhs3_c!1	 = "rad"
  rhs3_c&rad	 = radii

  rhs4_c!0       = "lev"
  rhs4_c&lev     = hyb_plot
  rhs4_c!1	 = "rad"
  rhs4_c&rad	 = radii

  rhs5a_c!0      = "lev"
  rhs5a_c&lev    = hyb_plot
  rhs5a_c!1	 = "rad"
  rhs5a_c&rad	 = radii

  rhs5b_c!0      = "lev"
  rhs5b_c&lev    = hyb_plot
  rhs5b_c!1	 = "rad"
  rhs5b_c&rad	 = radii

  lhsT_c!0       = "lev"
  lhsT_c&lev     = hyb_plot
  lhsT_c!1	 = "rad"
  lhsT_c&rad	 = radii

  rhsT_c!0       = "lev"
  rhsT_c&lev     = hyb_plot
  rhsT_c!1	 = "rad"
  rhsT_c&rad	 = radii

  mean_c!0       = "lev"
  mean_c&lev     = hyb_plot
  mean_c!1	 = "rad"
  mean_c&rad	 = radii

  eddy_c!0       = "lev"
  eddy_c&lev     = hyb_plot
  eddy_c!1	 = "rad"
  eddy_c&rad	 = radii

  vflux_c!0      = "lev"
  vflux_c&lev    = hyb_plot
  vflux_c!1      = "rad"
  vflux_c&rad    = radii

  vadv_c!0       = "lev"
  vadv_c&lev     = hyb_plot
  vadv_c!1       = "rad"
  vadv_c&rad     = radii

  agrd_c!0       = "lev"
  agrd_c&lev     = hyb_plot
  agrd_c!1       = "rad"
  agrd_c&rad     = radii

  vrad_c!0       = "lev"
  vrad_c&lev     = hyb_plot
  vrad_c!1	 = "rad"
  vrad_c&rad	 = radii

  vtan_c!0       = "lev"
  vtan_c&lev     = hyb_plot
  vtan_c!1	 = "rad"
  vtan_c&rad	 = radii

  vvel_c!0       = "lev"
  vvel_c&lev     = hyb_plot
  vvel_c!1	 = "rad"
  vvel_c&rad	 = radii

  aam_c!0        = "lev"
  aam_c&lev      = hyb_plot
  aam_c!1	 = "rad"
  aam_c&rad	 = radii

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
     opts_mbgt@cnFillPalette               = "ncl_default"
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
     opts_mbgt@cnLevels                    = (/-10.0, -5.0, -3.0, -2.0, -1.0, -0.5, -0.2, \
                                                0, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0/)
     opts_mbgt@cnFillColors                = (/0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17/)
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

  ; Turn off top x-axis; right y-xxis; maximise plot size
    opts_mbgt@tmXTOn                       = "False"
    opts_mbgt@tmYROn                       = "False"
    opts_mbgt@gsnMaximize                  = False  
    opts_mbgt@gsnAddCyclic                 = False

  ; Set tick marks explicitly
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

  ; Linearise the y-axis and set limits
    opts_mbgt@gsnYAxisIrregular2Linear	   = True
    opts_mbgt@trYMinF			   = 0.0
    opts_mbgt@trYMaxF			   = 17.6

  ; Set x-axis limits 
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
    opts_zero@gsnDraw                     = False       ; Do not draw the plot
    opts_zero@gsnFrame                    = False       ; Do no advance the frame
    opts_zero@gsnContourZeroLineThicknessF = 2.0
    opts_zero@cnLineLabelsOn               = False      ; Turn off line labels
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
    delete([/opts_vrad2@cnLevels,opts_vrad2@cnFillColors/])
    opts_vrad2@cnFillPalette	           = "BlueDarkRed18"
    opts_vrad2@cnLevels			   = (/-10., -5., -3., -2., -1., -0.5, -0.2, \
                                                0, 0.2, 0.5, 1., 2., 3., 5., 10./)
    opts_vrad2@cnFillColors		   = (/0,2,3,4,5,6,7,-1,-1,10,11,12,13,14,15,17/)
    opts_vrad2@lbLabelBarOn                = True
    opts_vrad2@lbBoxEndCapStyle            = "TriangleBothEnds"
    opts_vrad2@lbOrientation		   = "Vertical"
    opts_vrad2@vpXF			   = 0.07
    opts_vrad2@vpYF			   = 0.90
    opts_vrad2@vpWidthF                    = 0.15
    opts_vrad2@vpHeightF		   = 0.15

  ; Tangential wind 
    opts_vtan2				   = opts_mbgt
    delete([/opts_vtan2@cnLevels,opts_vtan2@cnFillColors/])
    opts_vtan2@cnFillPalette		   = "radar_new"

    if (new0 .eq. 1) then
      opts_vtan2@cnLevels                    = (/5.0, 10.0, 15.0, 17.5, 20.0, \
                                                22.5, 25.0, 27.5, 30.0, \
                                                32.5, 35.0, 37.5, 40.0/)
      opts_vtan2@cnFillColors                = (/-1,0,1,2,3,4,5,6,7,8,9,10,11,12/)
    else
      opts_vtan2@cnLevels                    = (/5.0, 10.0, 15.0, 20.0,\
                                                25.0, 30.0, 35.0, 40.0, \
                                                45.0, 50.0, 55.0, 60.0, \
                                                65.0/)
      opts_vtan2@cnFillColors                = (/-1,0,1,2,3,4,5,6,7,8,9,10,11,12/)
    end if

    opts_vtan2@lbLabelBarOn                = True
    opts_vtan2@lbBoxEndCapStyle            = "TriangleBothEnds"
    opts_vtan2@lbOrientation               = "Vertical"
    opts_vtan2@vpXF                        = 0.31
    opts_vtan2@vpYF                        = 0.90
    opts_vtan2@vpWidthF                    = 0.15
    opts_vtan2@vpHeightF                   = 0.15

  ; Vertical velocity
    opts_vvel2                             = opts_mbgt
    delete([/opts_vvel2@cnLevels,opts_vvel2@cnFillColors/])
    opts_vvel2@cnFillPalette               = "n11"
    opts_vvel2@cnLevels                    = (/0.1, 0.2, 0.3, 0.4, 0.5, 0.6, \
                                               0.7, 0.8, 1.0, 1.2, 1.5, 2.0/)
    opts_vvel2@cnFillColors                = (/-1,0,1,2,3,4,5,6,7,8,9,10,11/)
    opts_vvel2@lbLabelBarOn                = True
    opts_vvel2@lbBoxEndCapStyle            = "TriangleBothEnds"
    opts_vvel2@lbOrientation               = "Vertical"
    opts_vvel2@vpXF                        = 0.55
    opts_vvel2@vpYF                        = 0.90
    opts_vvel2@vpWidthF                    = 0.15
    opts_vvel2@vpHeightF                   = 0.15

  ; Absolute angular momentum
    opts_aam				   = opts_vtan2
    delete([/opts_aam@cnLevels,opts_aam@cnFillColors/])

    if (new0 .eq. 2) then
      opts_aam@cnLevels                      = (/0.5, 1.0, 1.25, 1.5, 1.75, 2.0, \
                                                 2.25, 2.5, 2.75, 3.0, 3.25, 3.5/)
      opts_aam@cnFillColors                  = (/-1,0,1,2,3,4,5,6,7,8,9,10,11/)
    elseif (new0 .eq. 1) then
      opts_aam@cnLevels                      = (/1.0, 1.5, 1.75, 2.0, \
                                                 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 4.0, 4.5/)
      opts_aam@cnFillColors                  = (/-1,0,1,2,3,4,5,6,7,8,9,10,11/)
    else
      opts_aam@cnLevels                      = (/1.0, 1.5, 2.0, 2.5, 3.0, \
                                                 3.5, 4.0, 4.5, 5.0, \
                                                 5.5, 6.0, 6.5/)
      opts_aam@cnFillColors                  = (/-1,0,1,2,3,4,5,6,7,8,9,10,11/)
    end if

    opts_aam@lbLabelBarOn                  = True
    opts_aam@lbBoxEndCapStyle              = "TriangleBothEnds"
    opts_aam@lbOrientation                 = "Vertical"
    opts_aam@vpXF                          = 0.79
    opts_aam@vpYF                          = 0.90
    opts_aam@vpWidthF                      = 0.15
    opts_aam@vpHeightF                     = 0.15

;===============================================================
; Resources to overlay panel labels without using 'gsn_panel'
;===============================================================

    txid			= new(14,graphic)
    amid	    	        = new(14,graphic)

  ; Label text resources 
    txres                       = True
    txres@txPerimOn             = True
    txres@txFontHeightF         = 0.010
    txres@txBackgroundFillColor = "White"

  ; Label position resources 
    amres 	   	        = True
    amres@amParallelPosF	= 0.5			; Right edge
    amres@amOrthogonalPosF	= -0.5			; Top edge
    amres@amJust		= "TopRight"

    if (ueh .eq. 0) then 
     panel_strings		= (/"a) vrad", "b) vtan", "c) vvel", "d) aam", \
                                     "e) U_magf", "f) U_mr", "g) U_mv", \
                                     "h) U_eagf", "i) U_dh", \
                                     "j) U_dz", "k) U_ev", "l) U_t", \
                                     "m) U_t + U_eh + U_mr", "n) All U tend"/)
    elseif (ueh .eq. 1) then 
     panel_strings              = (/"a) vrad", "b) vtan", "c) vvel", "d) aam", \
                                     "e) U_magf", "f) U_mr", "g) U_mv", \
                                     "h) U_eagf", "i) U_eh", \
                                     "j) U_dz + U_dh", "k) U_ev", "l) U_t", \
                                     "m) U_t + U_eh + U_mr", "n) All U tend"/)
    end if 

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
    output = "$nep/nepartak/images/cyl_coords/mbgt_rad_cxz_"+typ0+"_group"+grp+\
             "_sc"+sc+"_"+hy3+"km_"+w0+"w_mean"+mean0
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

  ; Mean agradient force [U_magf]
    opts_mbgt@vpXF       = 0.05
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel0      = gsn_csm_contour(wks,rhs3_c(:,:),opts_mbgt) ; [U_magf]
    cont0       = gsn_csm_contour(wks,rhs3_c(:,:),opts_zero) ; Draw contours

    panel0_ovr  = gsn_csm_contour(wks,vrad_c(:,:),opts_vrad) ; [u]
    panel0_ovr2 = gsn_csm_contour(wks,vvel_c(:,:),opts_vvel) ; [w]

    overlay(panel0, cont0)                                   ; Overlay contours
    overlay(panel0, panel0_ovr)
    overlay(panel0, panel0_ovr2)
    txid(4)     = gsn_create_text(wks, panel_strings(4), txres)
    amid(4)     = gsn_add_annotation(panel0, txid(4), amres)
    delete([/cont0/])

  ; Mean radial advection of mean radial momentum [U_mr]
    opts_mbgt@vpXF       = 0.24
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel1      = gsn_csm_contour(wks,lhs2_c(:,:),opts_mbgt)      ; [U_mr]
    cont0       = gsn_csm_contour(wks,lhs2_c(:,:),opts_zero) 	  ; Draw contours

    panel1_ovr  = gsn_csm_contour(wks,vrad_c(:,:),opts_vrad)   	  ; [U]
    panel1_ovr2 = gsn_csm_contour(wks,vvel_c(:,:),opts_vvel)      ; [w]

    overlay(panel1, cont0)                                        ; Overlay contours
    overlay(panel1, panel1_ovr)
    overlay(panel1, panel1_ovr2)
    txid(5)     = gsn_create_text(wks, panel_strings(5), txres)
    amid(5)     = gsn_add_annotation(panel1, txid(5), amres)
    delete([/cont0/])

  ; Minus the mean vertical advection of mean radial momentum [U_mv]
    opts_mbgt@vpXF       = 0.44
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel2      = gsn_csm_contour(wks,rhs1_c(:,:),opts_mbgt)      ; [U_mv]
    cont0       = gsn_csm_contour(wks,rhs1_c(:,:),opts_zero) 	  ; Draw contours

    panel2_ovr  = gsn_csm_contour(wks,vrad_c(:,:),opts_vrad)	  ; [U]
    panel2_ovr2 = gsn_csm_contour(wks,vvel_c(:,:),opts_vvel)	  ; [w]

    overlay(panel2, cont0)                                        ; Overlay contours
    overlay(panel2, panel2_ovr)
    overlay(panel2, panel2_ovr2)
    txid(6)     = gsn_create_text(wks, panel_strings(6), txres)
    amid(6)     = gsn_add_annotation(panel2, txid(6), amres)
    delete([/cont0/])

  ; Eddy agradient force [U_eagf]
    opts_mbgt@vpXF       = 0.64
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel3      = gsn_csm_contour(wks,rhs4_c(:,:),opts_mbgt)       ; [U_eagf]
    cont0       = gsn_csm_contour(wks,rhs4_c(:,:),opts_zero) 	   ; Draw contours
    
    panel3_ovr = gsn_csm_contour(wks,vrad_c(:,:),opts_vrad)	   ; [U]
    panel3_ovr2 = gsn_csm_contour(wks,vvel_c(:,:),opts_vvel) 	   ; [w]

    overlay(panel3, cont0)                                         ; Overlay contours
    overlay(panel3, panel3_ovr)
    overlay(panel3, panel3_ovr2)
    txid(7)     = gsn_create_text(wks, panel_strings(7), txres)
    amid(7)     = gsn_add_annotation(panel3, txid(7), amres)
    delete([/cont0/])

  ; Radial diffusive tendency of radial momentum [U_dh]
    opts_mbgt@vpXF       = 0.83
    opts_mbgt@vpYF       = 0.70
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15

    if (ueh .eq. 0) then     			  
	
     panel4      = gsn_csm_contour(wks,rhs5a_c(:,:),opts_mbgt)     ; [U_dh]
     cont0       = gsn_csm_contour(wks,rhs5a_c(:,:),opts_zero) 	   ; Draw contours

    elseif (ueh .eq. 1) then

     panel4      = gsn_csm_contour(wks,lhs3_c(:,:),opts_mbgt)      ; [U_eh]
     cont0       = gsn_csm_contour(wks,lhs3_c(:,:),opts_zero)      ; Draw contours

    end if 

    panel4_ovr  = gsn_csm_contour(wks,vrad_c(:,:),opts_vrad)	  ; [U]
    panel4_ovr2 = gsn_csm_contour(wks,vvel_c(:,:),opts_vvel) 	  ; [w]

    overlay(panel4, cont0)                                        ; Overlay contours
    overlay(panel4, panel4_ovr)
    overlay(panel4, panel4_ovr2)
    txid(8)     = gsn_create_text(wks, panel_strings(8), txres)
    amid(8)     = gsn_add_annotation(panel4, txid(8), amres)
    delete([/cont0/])

  ; Vertical diffusive tendency of radial momentum [U_dz]
    opts_mbgt@vpXF       = 0.05
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel5      = gsn_csm_contour(wks,rhs5b_c(:,:),opts_mbgt)     ; [U_dz]
    cont0       = gsn_csm_contour(wks,rhs5b_c(:,:),opts_zero) 	  ; Draw contours

    panel5_ovr  = gsn_csm_contour(wks,vrad_c(:,:),opts_vrad)	  ; [U]
    panel5_ovr2 = gsn_csm_contour(wks,vvel_c(:,:),opts_vvel) 	  ; [w]

    overlay(panel5, cont0)                                        ; Overlay contours
    overlay(panel5, panel5_ovr)
    overlay(panel5, panel5_ovr2)
    txid(9)     = gsn_create_text(wks, panel_strings(9), txres)
    amid(9)     = gsn_add_annotation(panel5, txid(9), amres)
    delete([/cont0/])

  ; Minus the eddy vertical advection of eddy radial momentum [U_ev]
    opts_mbgt@vpXF       = 0.24
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel6      = gsn_csm_contour(wks,rhs2_c(:,:),opts_mbgt)      ; [U_ev]
    cont0       = gsn_csm_contour(wks,rhs2_c(:,:),opts_zero) 	  ; Draw contours

    panel6_ovr  = gsn_csm_contour(wks,vrad_c(:,:),opts_vrad)	  ; [U]
    panel6_ovr2 = gsn_csm_contour(wks,vvel_c(:,:),opts_vvel) 	  ; [w]

    overlay(panel6, cont0)                                        ; Overlay contours
    overlay(panel6, panel6_ovr)
    overlay(panel6, panel6_ovr2)
    txid(10)    = gsn_create_text(wks, panel_strings(10), txres)
    amid(10)    = gsn_add_annotation(panel6, txid(10), amres)
    delete([/cont0/])

  ; Local tendency of the mean radial wind [U_t]
    opts_mbgt@vpXF       = 0.44
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel7      = gsn_csm_contour(wks,lhs1_c(:,:),opts_mbgt)      ; [U_t]
    cont0       = gsn_csm_contour(wks,lhs1_c(:,:),opts_zero) 	  ; Draw contours

    panel7_ovr  = gsn_csm_contour(wks,vrad_c(:,:),opts_vrad)	  ; [U]
    panel7_ovr2 = gsn_csm_contour(wks,vvel_c(:,:),opts_vvel) 	  ; [w]

    overlay(panel7, cont0)                                        ; Overlay contours
    overlay(panel7, panel7_ovr)
    overlay(panel7, panel7_ovr2)
    txid(11)    = gsn_create_text(wks, panel_strings(11), txres)
    amid(11)    = gsn_add_annotation(panel7, txid(11), amres)
    delete([/cont0/])

  ; Total LHS [U_t + U_eh + U_mr]
    opts_mbgt@vpXF       = 0.64
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel8      = gsn_csm_contour(wks,lhsT_c(:,:),opts_mbgt)      ; [U_t + U_eh + U_mr]
    cont0       = gsn_csm_contour(wks,lhsT_c(:,:),opts_zero) 	  ; Draw contours

    panel8_ovr  = gsn_csm_contour(wks,vrad_c(:,:),opts_vrad)	  ; [U]
    panel8_ovr2 = gsn_csm_contour(wks,vvel_c(:,:),opts_vvel) 	  ; [w]

    overlay(panel8, cont0)                                        ; Overlay contours
    overlay(panel8, panel8_ovr)
    overlay(panel8, panel8_ovr2)
    txid(12)    = gsn_create_text(wks, panel_strings(12), txres)
    amid(12)    = gsn_add_annotation(panel8, txid(12), amres)
    delete([/cont0/])

  ; Total RHS [All U tend]
    opts_mbgt@vpXF       = 0.83
    opts_mbgt@vpYF       = 0.50
    opts_mbgt@vpWidthF   = 0.15
    opts_mbgt@vpHeightF  = 0.15
    panel9      = gsn_csm_contour(wks,rhsT_c(:,:),opts_mbgt)      ; [All U tend]
    cont0       = gsn_csm_contour(wks,rhsT_c(:,:),opts_zero) 	  ; Draw contours

    panel9_ovr  = gsn_csm_contour(wks,vrad_c(:,:),opts_vrad)	  ; [U]
    panel9_ovr2 = gsn_csm_contour(wks,vvel_c(:,:),opts_vvel) 	  ; [w]

    overlay(panel9, cont0)                                        ; Overlay contours
    overlay(panel9, panel9_ovr)
    overlay(panel9, panel9_ovr2)
    txid(13)    = gsn_create_text(wks, panel_strings(13), txres)
    amid(13)    = gsn_add_annotation(panel9, txid(13), amres)
    delete([/cont0/])

;====================================================
; Finally, draw the plot with everything overlaid
;====================================================

    pres			= True 
    maximize_output(wks,pres)				   ; Calls 'draw' and 'frame'

  end if 
; End IF statement on whether to produce plots 

end 