pro setup_ps, name
  
   set_plot,'ps'
   !p.font=0
   !p.charsize=1.8
   device, filename = name, $
          /color, $
          /helvetica, $
          /inches, $
          xsize=19, $
          ysize=7, $
          /encapsulate, $
          yoffset=5, $
          bits_per_pixel = 16

end

pro read_nfar_data, file, t0, t1, data=data, utimes=utimes, freq=freq
   	READ_NU_SPEC, file, data,time,freq,beam,ndata,nt,dt,nf,df,ns, $
                tmin=t0*60.0, tmax=t1*60.0, fmin=20.0, fmax=33.0
        utimes=anytim(file2time(file), /utim) + time
        data = reverse(data, 2)
        freq = reverse(freq)
end

pro finterp_typeII, save=save, postscript=postscript

	path = '/databf2/nenufar-tf/ES11/2019/03/20190320_104900_20190320_125000_SUN_TRACKING_BHR/'
        file = 'SUN_TRACKING_20190320_104936_0.spectra'

	read_nfar_data, path+file, 39, 49, data=data, utimes=utimes, freq=freq
	   

	if keyword_set(postscript) then begin
		setup_ps, './nfar_'+time2file(utimes[0])+'_typeII.eps'
	endif else begin
		!p.charsize=1.8
		window, xs=1200, ys=700
	endelse	

	data = alog10(data)
	;data = constbacksub(data, /auto)
	;data = smooth(data,3)
	;------------------------------------------;
	;	Empty template to get black ticks
	;
	loadct, 0
        utplot, utimes, freq, yr=[57,20], /xs, /ys, xtitle='Time (UT)', ytitle='Frequency (MHz)', $
		title='NenuFAR-ES11 '+time2file(utimes[0], /date), pos=[0.1, 0.12, 0.9, 0.9], /normal, color=150
	
        ;------------------------------------------;
	;	     Plot spectrogram
	;
	loadct, 74
        reverse_ct
	spectro_plot, sigrange(data), utimes, freq, /xs, /ys, $
		ytitle=' ', xtitle=' ', yr=[57, 20], /noerase, XTICKFORMAT="(A1)", YTICKFORMAT="(A1)", $
		position=[0.1, 0.12, 0.9, 0.9], /normal

	;-----------------------------------------;
	;	         Flatten	
	;
	data = congrid(data, 10000, 1984)
	;shape = size(data)
	;prof = reform(reform(transpose(data), shape[4], 1))
	;delvar, data

	;-----------------------------------------;
	;	Each profile is evenly sampled
	;	in f, but unevenly in space.
	; 	This gets an even sample in space
	;	by interpolation.
	
	npoints=((freq*1e6/1.)/8980.0)^2.0 	
	rads = density_to_radius(npoints, model='newkirk')
	even_rads = interpol([rads[0], rads[-1]], n_elements(freq))
	nt=n_elements(data[*,0])-1
	for i=0, nt do begin
		prof = data[i, *]
		even_prof = interpol(prof, rads, even_rads)

		if i eq 0 then begin
			profiles = even_prof
		endif else begin
			profiles = [profiles, even_prof]
		endelse
		progress_percent, i, 0, nt
	endfor	
	def = even_rads[2]-even_rads[1]
	loadct, 0
	window, 1, xs=700, ys=700
	;plot, prof

	power = FFT_PowerSpectrum(profiles, def, FREQ=pfreq, /TUKEY, WIDTH=0.001, SIGNIFICANCE=signif)
	
	window, 2, xs=700, ys=700
	
	pfreq = alog10(pfreq)
	power = alog10(power)
	ind0 = closest(pfreq, 1.0)
	ind1 = closest(pfreq, 2.3)
	pfreq = pfreq[ind0:ind1]
	power = power[ind0:ind1]

	plot, pfreq, power, /xs, /ys, ytitle='log!L10!N(Power Rs!U-1!N)', xtitle='log!L10!N(k Rs!U-1!N)';, yr=[1e8, 1e12]
	
	result = linfit(pfreq, power)
	print, result[1]
	pfsim = interpol([pfreq[0], pfreq[-1]], 100)
	powsim = result[0] + result[1]*pfsim
	set_line_color
	plot, pfreq, power, /xs, /ys, ytitle='log!L10!N(Power Rs!U-1!N)', xtitle='log!L10!N(k Rs!U-1!N)', $
		title='Spectral index:'+ string(result[1]);, yr=[1e8, 1e12]
	oplot, pfsim, powsim, color=10

	if keyword_set(postscript) then begin
		device, /close
		set_plot, 'x'
	endif	
	  

stop
END
