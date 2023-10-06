;============================================================
; Loop over pressure levels to calculate vertical profile
;============================================================

   cp = 0                            ; Counter variable

   do ilev = 0, plevs-1              ; Loop over chosen pressure levels

    p = pres(ilev)
    print("Working on: " + p + " hPa")

   ;========================================
   ; Model levels ['pc' and 'pf' streams]
   ;========================================

     hy = sprintf("%0.0f",hyb(ilev))
     print("Working on hybrid level " + ilev+ " ("+hy+" m AGL)")
     m0 = lev_arr(ilev)       ; Selected model levels, chosen to match pressure levels

   ; Calculate components of friction term
     ip1 = min((/m0+1,lev_arr(l1)/)) ; Account for top and bottom levels
     im1 = max((/m0-1,lev_arr(l0)/))

   ; Current level [wind stress]
     taux   = strx(m0,0:d1,0:d2)
     tauy   = stry(m0,0:d1,0:d2)

   ; Density ['pf' stream]
     rho    = rho1(m0,0:d1,0:d2)
     rho_a = (1/rho) ; Specific volume (inverse density)

   ; Levels above/below for vertical derivatives [wind stress]
     tau_x1 = strx(ip1,0:d1,0:d2)
     tau_x0 = strx(im1,0:d1,0:d2)
     tau_y1 = stry(ip1,0:d1,0:d2)
     tau_y0 = stry(im1,0:d1,0:d2)

   ; Vertical derivative of 'strx'
     dtaux_dz = (tau_x1 - tau_x0) / (hyb(ip1) - hyb(im1))

   ; Vertical derivative of 'stry'
     dtauy_dz = (tau_y1 - tau_y0) / (hyb(ip1) - hyb(im1))

   ; X-component of friction
     fric_x_plane           = rho_a(:,:) * dtaux_dz
     copy_VarCoords(taux,fric_x_plane)
     copy_VarMeta(taux,fric_x_plane)
     fric_x_plane@units     = "m s~S~-2~N~"
     fric_x_plane@title     = "x-component of friction"
     fric_x_plane@long_name = "x-component of friction"
     fric_x_plane@name      = "fric_x"

   ; Y-component of friction
     fric_y_plane           = rho_a(:,:) * dtauy_dz
     copy_VarCoords(tauy,fric_y_plane)
     copy_VarMeta(tauy,fric_y_plane)
     fric_y_plane@units     = "m s~S~-2~N~"
     fric_y_plane@title     = "y-component of friction"
     fric_y_plane@long_name = "y-component of friction"
     fric_y_plane@name      = "fric_y"

;====================================================
; Define arrays to hold ensemble of box positions
;====================================================

     fric_x_box = new((/lsize,size0,size1/),float)
     fric_y_box = new((/lsize,size0,size1/),float)

;===================================================================
; Loop over perturbed box positions, either side of the TC centre
;===================================================================

   do i = 0, dimsizes(lon_box)-1    ; Longitude
    do j = 0, dimsizes(lat_box)-1  ; Latitude

     a    = (dimsizes(lon_box) * i) + j
     print_clock("Working on box number "+a)

;===================================================
; Troubleshooting (size of arrays on LHS and RHS)
;===================================================

     dim1 = dimsizes(fric_x_plane({lat_box(j)-r:lat_box(j)+r},\
                                  {lon_box(i)-r:lon_box(i)+r}))

     dim2 = dimsizes(fric_y_plane({lat_box(j)-r:lat_box(j)+r},\
                                  {lon_box(i)-r:lon_box(i)+r}))

   ; 'size0' = size of dimension 1 of fric_x_box(a,:,:) - LATITUDE
   ; 'size1' = size of dimension 2 of fric_x_box(a,:,:) - LONGITUDE
   ; 'dim1(0)' = size of dimension 1 of fric_x_plane({lat_box(j)-r...}) - LATITUDE
   ; 'dim1(1)' = size of dimension 2 of fric_x_plane({lat_box(j)-r...}) - LONGITUDE
   ; 'dim2(0)' = size of dimension 1 of fric_y_plane({lat_box(j)-r...}) - LATITUDE
   ; 'dim2(1)' = size of dimension 2 of fric_y_plane({lat_box(j)-r...}) - LONGITUDE

     print("LHS 'fric_x/y' dimensions: "+size0+" lat. points; "+size1+" lon. points")
     print("RHS 'fric_x' dimensions: "+dim1(0)+" lat. points; "+dim1(1)+ " lon. points")
     print("RHS 'fric_y' dimensions: "+dim2(0)+" lat. points; "+dim2(1)+ " lon. points")

     dxy = 0.02

    ;=========================================================
    ; Make sure all arrays are the same size (fric_x_plane)
    ;=========================================================

      if (dim1(0) .lt. size0) then
        print("RHS = "+dim1(0)+" ; LHS = "+size0)
        print("RHS latitude array is smaller than LHS!")
        fric_x_box(a,:,:) = fric_x_plane({lat_box(j)-r:lat_box(j)+r+dxy},\
                                         {lon_box(i)-r:lon_box(i)+r})
      elseif (dim1(0) .gt. size0) then
        print("RHS = "+dim1(0)+" ; LHS = "+size0)
        print("RHS latitude array is bigger than LHS!")
        fric_x_box(a,:,:) = fric_x_plane({lat_box(j)-r:lat_box(j)+r-dxy},\
                                         {lon_box(i)-r:lon_box(i)+r})
      elseif (dim1(1) .lt. size1) then
        print("RHS = "+dim1(1)+" ; LHS = "+size1)
        print("RHS longitude array is smaller than LHS!")
        fric_x_box(a,:,:) = fric_x_plane({lat_box(j)-r:lat_box(j)+r},\
                                         {lon_box(i)-r:lon_box(i)+r+dxy})
      elseif (dim1(1) .gt. size1) then
        print("RHS = "+dim1(1)+" ; LHS = "+size1)
        print("RHS longitude array is bigger than LHS!")
        fric_x_box(a,:,:) = fric_x_plane({lat_box(j)-r:lat_box(j)+r},\
                                         {lon_box(i)-r:lon_box(i)+r-dxy})
      else
        fric_x_box(a,:,:) = fric_x_plane({lat_box(j)-r:lat_box(j)+r},\
                                         {lon_box(i)-r:lon_box(i)+r})
      end if

    ;=========================================================
    ; Make sure all arrays are the same size (fric_y_plane)
    ;=========================================================
      if (dim2(0) .lt. size0) then
        print("RHS = "+dim2(0)+" ; LHS = "+size0)
        print("RHS latitude array is smaller than LHS!")
        fric_y_box(a,:,:) = fric_y_plane({lat_box(j)-r:lat_box(j)+r+dxy},\
                                         {lon_box(i)-r:lon_box(i)+r})
      elseif (dim2(0) .gt. size0) then
        print("RHS = "+dim2(0)+" ; LHS = "+size0)
        print("RHS latitude array is bigger than LHS!")
        fric_y_box(a,:,:) = fric_y_plane({lat_box(j)-r:lat_box(j)+r-dxy},\
                                         {lon_box(i)-r:lon_box(i)+r})
      elseif (dim2(1) .lt. size1) then
        print("RHS = "+dim2(1)+" ; LHS = "+size1)
        print("RHS longitude array is smaller than LHS!")
        fric_y_box(a,:,:) = fric_y_plane({lat_box(j)-r:lat_box(j)+r},\
                                         {lon_box(i)-r:lon_box(i)+r+dxy})
      elseif (dim2(1) .gt. size1) then
        print("RHS = "+dim2(1)+" ; LHS = "+size1)
        print("RHS longitude array is bigger than LHS!")
        fric_y_box(a,:,:) = fric_y_plane({lat_box(j)-r:lat_box(j)+r},\
                                         {lon_box(i)-r:lon_box(i)+r-dxy})
      else
        fric_y_box(a,:,:) = fric_y_plane({lat_box(j)-r:lat_box(j)+r},\
                                         {lon_box(i)-r:lon_box(i)+r})
      end if

;====================================================================================
; Second round of calculations (Horizontal advection, relative vorticity, tilting)
;====================================================================================

   ; Calculate terms on northern/southern boundaries
     fric_x0     = dim_sum_n_Wrap((fric_x_box(a,:,:) * ddx_box(a,:,:)),1)
     fric_s      = fric_x0(0)              ; Friction (south)
     fric_n      = -fric_x0(numPOINTS_y-1) ; Friction (north)

   ; Calculate terms on eastern/western boundaries
     fric_y0     = dim_sum_n_Wrap((fric_y_box(a,:,:) * ddy_box(a,:,:)),0)
     fric_w      = -fric_y0(0)             ; Friction (west)
     fric_e      = fric_y0(numPOINTS_x-1)  ; Friction (east)

;========================
; Further calculations
;========================

   ; Calculate horizontal advection from line integrals (4 sides of box)
     fric_tenT   = fric_s + fric_n + fric_w + fric_e
     fric_tenT   = fric_tenT / d_ave(a) ; Divide by box area [Rios-Berrios et al. 2016]

     fric_ten(a,it,cp)      = fric_tenT              ; Friction

     delete([/fric_tenT,fric_x0,fric_y0/])
     delete([/fric_x_plane, fric_y_plane,fric_x_box,fric_y_box/])
     delete([/dtaux_dz, dtauy_dz/])

   ; WHOLE SET OF VARIABLES TO DELETE AT END OF LOOP (BEFORE EDITING)
   ; Tidy up (delete variables from within time loop before next iteration)
     delete([/lon,lat,lon1,lat1,u,v,vort,div,f0,avo,t,z,w,box/])
     delete([/time,times,utc_date,year,month,day,hour,minute,second/])
     delete([/lon1F,lat1F,lonF,latF,hyb,strx,stry,dxF,dyF,rho1/])

;=====================================================
; Calculate term averages over all ensemble members
;=====================================================

     fric_final = new((/numINPUT_a,plevs/),"float")

   ; Average over ensemble members (do we still need this step?)
     fric_final(:,:) = dim_avg_n_Wrap(fric_ten(:,:,:),0)

;=========================================
; Integrate over Nepartak's life cycle
;=========================================

 ; Each calculation produces an array of [plevs] values [pressure levels]
 ; We do this calculation [numINPUT_a] times
 ; Create arrays with [numINPUT_a] * [plevs] elements

   fric_intA = new((/numINPUT_a,plevs/),float)

   do ic = 0, csize0-1  ; 1 h integrations over the whole of Nepartak's life cycle

    cp1 = min((/ic+1,csize0-1/))
    cm1 = max((/ic-1,0/))

     f_str = sprinti("%0.3i",ic)+"h"    ; Formatted string for output file

     fric_intA(ic,:) = dim_sum_n_Wrap((fric_final(ic:cp1,:) * 3600), 0)

     if (fct .eq. 4) then
      fric_intA(ic,:) = fric_intA(ic,:) * (10 ^ 4)
     elseif (fct .eq. 5) then
      fric_intA(ic,:) = fric_intA(ic,:) * (10 ^ 5)
     end if

    end do               ; End time loop (do ic = 0, csize0-1)

 ; Add metadata after calculating tendencies
   if (fct .eq. 4) then
    fric_intA@units = "10~S~-4~N~ s~S~-2~N~"
   elseif (fct .eq. 5) then
    fric_intA@units = "10~S~-5~N~ s~S~-2~N~"
   end if

   tot_intA  = str_intA + hadv_intA + tilt_intA + fric_intA

   fric_intA!0    = "ens"
   fric_intA&ens  = ens_arr
   fric_intA!1    = "time"
   fric_intA&time = time_int
   fric_intA!2    = "lev"
   fric_intA&lev  = pres        ;;; CHANGE TO MODEL LEVELS (?) ;;;

;======================================
; Reorder dimensions before plotting
;======================================

 ; Define new arrays for final plotting
   fric_int = new((/plevs,numINPUT_a/),float)

 ; Reorder from [time * pressure] to [pressure * time]
   fric_int(:,:) = fric_intA(lev|:,ens|:,time|:)

;=================
; Plot the data
;=================

     plot_fric0   = gsn_csm_pres_hgt(wks,fric_int(:,0,:),opts_bdgt)

   ; First row
     panel(5)  = plot_fric0