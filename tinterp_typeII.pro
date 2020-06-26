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

pro interp_typeII, save=save, postscript=postscript

	path = '/databf2/nenufar-tf/ES11/2019/03/20190320_104900_20190320_125000_SUN_TRACKING_BHR/'
        file = 'SUN_TRACKING_20190320_104936_0.spectra'

	read_nfar_data, path+file, 37.5, 45, data=data, utimes=utimes, freq=freq
	   

	if keyword_set(postscript) then begin
		setup_ps, './nfar_'+time2file(utimes[0])+'_typeII.eps'
	endif else begin
		!p.charsize=1.8
		window, xs=1200, ys=700
	endelse	

	data = smooth(data,3)
	;------------------------------------------;
	;	Empty template to get black ticks
	;
	loadct, 0
        utplot, utimes, freq, yr=[57,20], /xs, /ys, xtitle='Time (UT)', ytitle='Frequency (MHz)', $
		title='NenuFAR-ES11 '+time2file(utimes[0], /date), pos=[0.1, 0.12, 0.9, 0.9], /normal, color=150
	
        ;------------------------------------------;
	;		Plot spectrogram
	;
	loadct, 74
        reverse_ct
	spectro_plot, sigrange(data), utimes, freq, /xs, /ys, $
		ytitle=' ', xtitle=' ', yr=[57, 20], /noerase, XTICKFORMAT="(A1)", YTICKFORMAT="(A1)", $
		position=[0.1, 0.12, 0.9, 0.9], /normal

	stop
	point, x, y
	x0 = closest(utimes, x[0])
	x1 = closest(utimes, x[1])
	y0 = closest(freq, y[0])
        y1 = closest(freq, y[1])
	
	xpoints = interpol([x0,x1],1000)
	ypoints = interpol([y0,y1],1000)
	prof = interpolate(data, xpoints, ypoints)
	
	loadct, 0
	window, 1, xs=700, ys=700
	plot, prof

	fpoints = interpol([y[0], y[1]], 1000)

	f = FFT_PowerSpectrum(prof, 1.27e-4, FREQ=freq, $
  		/TUKEY, WIDTH=0.01, SIGNIFICANCE=signif)
	window, 2, xs=700, ys=700
	plot, freq, f, /xs, /ys, /ylog, /xlog, xr=[1e1, 1e4], yr=[1e8, 1e12]
	
	xr = 10^interpol([1,4], 100)
	yr = 1e15*xr^(-5/3.)
	set_line_color
	oplot, xr, yr, color=5, linestyle=1


	if keyword_set(postscript) then begin
		device, /close
		set_plot, 'x'
	endif	
	  

stop
END
