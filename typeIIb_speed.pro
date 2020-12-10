pro setup_ps, name, xsize=xsize, ysize=ysize
  
   set_plot,'ps'
   !p.font=0
   !p.charsize=1.5
   device, filename = name, $
          /color, $
          /helvetica, $
          /inches, $
          xsize=xsize, $
          ysize=ysize, $
          /encapsulate, $
          yoffset=5, $
          bits_per_pixel = 16

end

pro read_nfar_data, file, t0, t1, f0, f1, data=data, utimes=utimes, freq=freq
   	READ_NU_SPEC, file, data, time, freq, beam, ndata, nt, dt, nf, df, ns, $
                tmin=t0*60.0, tmax=t1*60.0, fmin=f0, fmax=f1, fflat=3, ntimes=16, nchannels=512
        
	utimes=anytim(file2time(file), /utim) + time
        data = reverse(data, 2)
        freq = reverse(freq)
end

pro plot_spectro, data, time, freq, f0, f1, posit

	loadct, 0
        utplot, time, freq, yr=[f1,f0], /xs, /ys, xtitle='Time (UT)', ytitle='Frequency (MHz)', $
                title=' ', pos=posit, /normal, color=150, $
                xr=[time[0], time[-1]]
	;------------------------------------------;
        ;            Plot spectrogram
        ;
        loadct, 74
        reverse_ct
        spectro_plot, sigrange(data), time, freq, /xs, /ys, $
                ytitle=' ', xtitle=' ', yr=[f1, f0], /noerase, XTICKFORMAT="(A1)", YTICKFORMAT="(A1)", $
                position=posit, /normal, xr=[time[0], time[-1]]


END

pro typeIIb_speed, save=save, plot_ipsd=plot_ipsd, postscript=postscript, rebin=rebin

	; PSD of first type II. Code working.

	; This version of the code excludes certain powerlaw fits based on their p-value.

	; The version with drift extracts an intensity profile along a drifting structure within the radio burst. 
	
	; More appropriate when structures actually drift in frequency time.

	path = '/databf2/nenufar-tf/ES11/2019/03/20190320_104900_20190320_125000_SUN_TRACKING_BHR/'
        file = 'SUN_TRACKING_20190320_104936_0.spectra'

	
	t0 = 33.0
        t1 = 35.0
        f0 = 20.0
        f1 = 25.0
	read_nfar_data, path+file, t0, t1, f0, f1, data=data, utimes=utimes, freq=freq
	   

	!p.charsize=1.8
	window, 0, xs=1600, ys=600

	posit = [0.05, 0.15, 0.7, 0.9]
      
	all_data = data ;smooth(data, 3)
        all_times = utimes
        all_freqs = freq
	hbspeeds = 0
	hbalphas = 0	
	
	;------------------------------------------;
	;	     Plot spectrogram
	;
        data = all_data
	utimes = all_times
	freq = all_freqs
	plot_spectro, data, utimes, freq, f0, f1, posit

	;------------------------------------------;
	;
	;	Extract profile along drift
	;
	loadct, 0
	point, tpoints, fpoints, /data
	
	npoints=((fpoints*1e6/1.)/8980.0)^2.0
        rads = density_to_radius(npoints, model='newkirk')
	

	;----------------------------------;
	;	Fit burst
	;	
	window, 1, xs=500, ys=500
	tsec = tpoints - tpoints[0]
	fdrift = (linfit(tsec, fpoints))[1]
	plot, tsec, rads, xtitle='Time (s)', ytitle='Heliocentric distance', /xs, /ys
	result = linfit(tsec, rads)
	ymodel = result[0] + result[1]*tsec
	set_line_color
	oplot, tsec, ymodel, color=5
	speed = result[1]*695e3
	print, 'Speed (km/s)'+string(speed)
	print, 'Drift rate (MHz/s)'+string(fdrift)

	stop	
END
