function freq_to_dens, freq

	; frequency of plasma oscillation in Hz
	
	e_mass = 9.10938291D-31*1000.0d ; grams (cgs)
	e_q = 4.8032e-10			;statcouloumb (cgs)
	omega = freq*2.0d*!PI
	
	
	n_e = ( omega^2.0 * e_mass) / (4.0d*!PI*e_q^2.0)
	
	
	return, n_e ;returned in cm^-3

END
