function rel_energy, c_fraction, proton=proton

	; c_fraction is light speed fraction.

	; Default electron
	e_mass = 9.10938188D-31  ;kg
	e_charge = 1.602D-19	 ;C
	c = 2.99792458D8		 ;m/s
	rest_E = e_mass*(c^2.0)  ;J

	; Make proton if needed
	if keyword_set(proton) then begin
		p_mass = 1.67e-27		 ;kg
		e_charge = 1.602D-19	 ;C
		c = 2.99792458D8		 ;m/s
		rest_E = p_mass*(c^2.0)  ;J
	endif	


	kin_e = rest_E/sqrt(1.0 - c_fraction^2.0) - rest_E
	rel_kine = (kin_e/(e_charge)) / 1.0D3  ;relativistic kinetic energy (keV)
	print, 'Relativistic kinetic energy: '+string(rel_kine)+' keV'

	totE = rest_E/sqrt(1.0 - c_fraction^2.0) 
	rel_tote =  (totE/(e_charge)) / 1.0D3  ;total relativistic eneregy (keV), including electron rest mass
	print, 'Relativistic kinetic energy including electron rest mass: '+string(rel_tote)+' keV'


	rel_e = [rel_kine, rel_tote]
	return, rel_e

END
