; Script to calculate all 24-h periods of RI for a 12-member RA1T MetUM ensemble

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; ncl dat=\"03T00\"  cat=\"slow\" spd=20 ts=48 tf=72 nep_cp_ens_ri.ncl

begin

; Read in 10-m windspeed data from text files 
  diri	      = "$ar/text"
  fili_a1     = systemfunc("cd "+diri+" ; ls cp_"+dat+"*10m.txt")
  fili_a      = diri+"/"+fili_a1	                ; All file paths 
  numINPUT    = dimsizes(fili_a)			; Number of simulations
  numTIMES    = 120					; Number of times 

  setvalues NhlGetWorkspaceObjectId
    "wsMaximumSize" : 10000000000
  end setvalues

  print_clock("Working on initialisation time: "+dat)

; Troubleshooting 
  if (tf .gt. 95) then 
   print("Too close to end of simulation: can't calculate 24-h change...")
   exit()  
  end if 

;==========================
; Loop over input files
;==========================

  do nf = 0, numINPUT-1     ; 12 ensemble members in directory

 ; Extract string for each ensemble member
   str1    = str_split(fili_a(nf),"_")
   ens_str = str1(2)

   print("Working on ensemble member: "+ens_str)

 ; Read in 10-m windspeed from text files 
   vec_10m = asciiread(fili_a(nf),120,"float")

;====================================================================
; Calculate 24 h difference in windspeed - how many periods of RI? 
;====================================================================

  ; Loop over all times in file (stop 24 h before end)
    do it = ts, tf

     t1  = it			     ; Index of current time (T+0)
     t2  = it+24 		     ; Index of time 24 h ahead (T+24)
     dv  = vec_10m(t2) - vec_10m(t1) ; Compute 24 h change in windspeed 
     dv0 = abs(dv)

     if (cat .eq. "fast") then 

      if (dv .ge. spd) then

       dvv = new(25,"float")
       ct0 = 0 
       
       do ct = t1, t2
        dvv(ct0) = vec_10m(ct+1) - vec_10m(ct)
	ct0 = ct0 + 1
       end do 

       xt = num(dvv .gt. 0.0)		; Calculate the number of 1-h periods with +ve tendency
       xu = num(dvv .lt. -1.0)		; Calculate number of 1-h periods with -ve tendency 

       if ( xt .ge. thr0 .and. xu .lt. thr1) then 

        print("RI occurs in simulation "+ens_str+" between t = "+t1+" and t = "+t2)
	print("1-h periods in which windspeed increases: "+xt)
	print("1-h periods in which windspeed decreases by less than 0.5 m/s/h: "+xu)
        print("The windspeed increased by "+dv+" m/s in this 24-h period.")

       end if 

      end if 

     elseif (cat .eq. "slow") then 

      if (dv0 .le. spd) then
       print("Steady-state period occurs in simulation "+ens_str+" between t = "+t1+" and t = "+t2)
       print("The change in windspeed was only "+dv0+" m/s in this 24-h period.")
      end if

     end if 

     delete([/dv,dv0/])

    end do			     ; End loop over times (do it = 0, numTIMES-1) 
    	
  end do			     ; End loop over input files (do nf = 0, numINPUT-1)

end 