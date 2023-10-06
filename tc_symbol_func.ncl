  ; Add NH tropical cyclone symbol to index
    tc0 = NhlNewMarker(wks, "p", 137, 0.0, 0.0, 1.0, 1.0, 0.0)

    mres0               = True
    mres0@gsMarkerIndex = tc0         ; TC symbol
    mres0@gsMarkerSizeF = 20.0        ; Marker size
    mres0@gsMarkerColor = "red3"      ; Marker colour
    mres0@gsMarkerThicknessF = 5.0    ; Marker thickness (increase)

    marker    = gsn_add_polymarker(wks,spd_plot,centre(ct,1),centre(ct,0),mres0)
