; Only label single contour line                                                                     
;==================================                                                                  

;    do pp = 0, 1            ; Loop over both panels                                                 

;    ; Retrieve the contour levels (and put into 'levels')                                           
;      getvalues slp_plot_gl                                                                         
;        "cnLevels" : levels                                                                         
;      end getvalues                                                                                 

;    ; Create an array to set a flag for each contour level                                          
;      nlevs = dimsizes(levels)                                                                      
;      flags = new( (/nlevs/), string)                                                               

;    ; Loop through each level and set the corresponding                                             
;      do ii = 0, nlevs-1                                                                            
;       if (levels(ii) .ne. 25) then                                                                 
;        flags(ii) = "LineOnly"                                                                      
;       else                                                                                         
;        flags(ii) = "LineAndLabel"                                                                  
;       end if                                                                                       
;      end do                                                                                        

;    ; Set new resource in the plot we've already created, and redraw                                
;      setvalues slp_plot_gl                                                                         
;       "cnLevelFlags" : flags                                                                       
;      end setvalues                                                                                 

;    end do