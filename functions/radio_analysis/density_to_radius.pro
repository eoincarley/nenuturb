function density_to_radius, dens, $
          model = model, $
          fold = fold
  
    ; dens in cm^-3        
          
    if keyword_set(fold) then fold=fold else fold=1       

    radius = (findgen(1000)*(10. - 1.)/999.) + 1. 				; Rsun in heliocentric distance

    if model eq 'saito' then begin
        ;-----------------;
        ;		Saito QS
        ;
        c1 = 1.36e6			  ; Params for Saito quiet sun background
        c2 = 1.68e8
        d1 = 2.14
        d2 = 6.13

        ;--------------------;
        ;   Saito Eq Hole
        ;
        ;c1 = 5.27e6
        ;c2 = 3.54e6
        ;d1 = 3.3
        ;d2 = 5.8

        ;--------------------;
        ;   Saito Polar Hole
        ;
        ;c1 = 3.15e6
        ;c2 = 1.60e6
        ;d1 = 4.71
        ;d2 = 3.01

        den_model = c1*radius^(-1.0*d1) + c2*radius^(-1.0*d2)  
    endif  

    if model eq 'baum' then begin 
        den_model =  1.0e8*(  2.99*(radius)^(-16.0) + 1.55*(radius)^(-6.0) + 0.036*(radius)^(-1.5)  ) 
    endif 
    
    if model eq 'newkirk' then begin   
        N=4.2e4
        den_model =  N*10.0^(4.32/radius) 
    endif  
    
    if model eq 'mann' then begin     
        RESOLVE_ROUTINE, 'sw_mann'
        sw_mann, radius, n = den_model
    endif 

    if model eq 'leblanc' then begin     
        c1 = 3.3e5
        c2 = 4.1e6
        c3 = 8.0e7
        den_model = c1*radius^(-2.0) + c2*radius^(-4.0) + c3*radius^(-6.0)
    endif 

    if model eq 'st_hilaire_0' then begin     ; For fundamental emission
         den_model = 5.0e6*(radius-1.0)^(-2.0)
    endif 

    if model eq 'st_hilaire_1' then begin     ; For 1st harmonic emission
         den_model = 1.2e6*(radius-1.0)^(-2.0)
    endif 
  
    if model eq 'tcd' then begin     
        restore,'~/Data/2011_sep_22/density_mag/tcd_model_20110922.sav'
        indices = where(tcd_rads gt 1.01)
        radius = tcd_rads[indices]
        den_model = tcd_ne[indices]
    endif  

    ; Plane parallel hydrostatic equilibrium model.
    if model eq 'hydro_stat' then begin 
        n0 = 2.0e9*1e6   ;per cubic meter
        k = 1.38e-23    ;J/K
        T = 1.1e6 ;1.0e6         ;K
        m = 1.66e-27    ;kg (amu)
        gearth=9.88     ;m/s/s
        gsun = 27.94*gearth     ;m/s/s
        H = (k*T)/(m*gsun)      ;m
        H = H/6.95e8            ;Rsun
        den_model = n0*exp(-1.0*((radius-1.0)/H))   ;m^-3
        den_model = den_model/1e6           ;cm^-3
    endif

    den_model = den_model*fold  
    rads = interpol(radius, den_model, dens)

    return, rads 	;Rsun
	
	
END
