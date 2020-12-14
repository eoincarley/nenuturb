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
                tmin=t0*60.0, tmax=t1*60.0, fmin=f0, fmax=f1, fflat=3, ntimes=16, $
		ex_chan=[0], /fill, /exactfreq
       
       	stop	
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

pro typeIIa_hbkins, save=save, plot_ipsd=plot_ipsd, postscript=postscript, rebin=rebin

	; PSD of first type II. Code working.

	; This version of the code excludes certain powerlaw fits based on their p-value.

	; The version with drift extracts an intensity profile along a drifting structure within the radio burst. 
	
	; More appropriate when structures actually drift in frequency time.

	path = '/databf2/nenufar-tf/ES11/2019/03/20190320_104900_20190320_125000_SUN_TRACKING_BHR/'
        file = 'SUN_TRACKING_20190320_104936_0.spectra'

	
	t0 = 33.0
        t1 = 35.0
        f0 = 32.0
        f1 = 55.0
	read_nfar_data, path+file, t0, t1, f0, f1, data=data, utimes=utimes, freq=freq
	   

	!p.charsize=1.8
	window, xs=1600, ys=600

	posit = [0.05, 0.15, 0.7, 0.9]
      
	all_data = data ;smooth(data, 3)
        all_times = utimes
        all_freqs = freq
	hbspeeds = 0
	hbalphas = 0	
	for i=0, 8 do begin
	wset, 0 
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
	cursor, tpoint, fpoint, /data
	;save, tpoint, fpoint, filename='herringbone_tfpoint_H.sav'
	hdata = data
	hfreq = freq

	findex = (where(freq le fpoint))[0]
	tindex = (where(utimes ge tpoint))[0]

	tpix = 5.0 ; Search for peak in the next f channel at +/- 7 pix of peak in previous channel.
	it0 = tindex-tpix
	it1 = tindex+tpix
	iprof = data[it0:it1, findex]
	utprof = utimes[it0:it1]

	itpeak = where(iprof eq max(iprof)) + it0
	
	set_line_color
	plots, utimes[itpeak], freq[findex], /data, psym=1, symsize=2, color=0	

	iburst = data[itpeak, findex]
	iburst_max = iburst
	isample = iburst_max
	fburst = freq[findex]
	tburst = utimes[itpeak]
	print, iburst_max
	while isample gt iburst_max*0.5 do begin	
		it0 = itpeak-tpix
        	it1 = itpeak+tpix
		iprof = data[it0:it1, findex]
        	utprof = utimes[it0:it1]
        	itpeak = where(iprof eq max(iprof)) + it0
		plots, utimes[itpeak], freq[findex], /data, psym=1, symsize=0.1, color=5

		isample = data[itpeak, findex]
		iburst = [iburst, isample]
		fburst = [fburst, freq[findex]]
		tburst = [tburst, utimes[itpeak]]
		findex = findex-1
		;print, isample
	endwhile

	freq = fburst[1:-1]
	prof = iburst[1:-1]
	tburst = tburst[1:-1]
	
	;---------------------------------------;
	;	Now perform PSD on iburst
	;	
	;-----------------------------------------;
        ;       Each profile is evenly sampled
        ;       in f, but unevenly in space.
        ;       This gets an even sample in space
        ;       by interpolation.
        npoints = ((freq*1e6/2.)/8980.0)^2.0
        rads = density_to_radius(npoints, model='newkirk')
        even_rads = interpol([rads[0], rads[-1]], n_elements(freq))
        nt = n_elements(data[*,0])-1
        def = abs(even_rads[2]-even_rads[1])
        pspecerr = 0.05
	wavenum0 = 1.0+alog10(2.0*!pi)
        wavenum1 = 3.5+alog10(2.0*!pi)
	rsunMm = 696.34 ; Mm
	
	stop
        ;----------------------------------------------;
        ;
	;  Get evenly sampled in space and perform PSD
        ;
        even_prof = interpol(prof, rads, even_rads)
        even_prof = even_prof/max(even_prof)
	
	;loadct, 0
        ;wset, 1
	;plot, even_rads, even_prof, /xs, /ys, /normal, $
        ;        xtitle=' ', ytitle='Intensity', XTICKFORMAT="(A1)", xticklen=-1e-8

        ;axis, xaxis=0, xr = [even_rads[0], even_rads[-1]], /xs, xtitle='Heliocentric distance (R!Ls!N)'
        ;axis, xaxis=1, xr = [even_rads[0], even_rads[-1]]*rsunMm, /xs, xtitle='(Mm)'
	
	;stop
        power = FFT_PowerSpectrum(even_prof, def, FREQ=pfreq,$
                /tukey, width=0.0025, sig_level=0.01, SIGNIFICANCE=signif)

        pfreq = alog10(pfreq*2.0*!pi) ; x 2pi to get wavenumber from 1/lambda
	power = alog10(power)
        ind0 = closest(pfreq, wavenum0)
        ind1 = closest(pfreq, wavenum1)
        pfreq = pfreq[ind0:ind1]
        power = power[ind0:ind1]
        sigcutoff = alog10(signif[0])

	;----------------------------------;
    	;       Fit PSD and plot
    	;
	set_line_color
	plot, 10^pfreq, 10^power, /xlog, /ylog, /xs, /ys, ytitle='PSD', $
            xtitle=' ', thick=2, xr = [10^wavenum0, 10^wavenum1], yr=10.0^[-7.0, -2.0], $
                /noerase, position=[0.76, 0.15, 0.96, 0.9], psym=10, $
            XTICKFORMAT="(A1)", xticklen=1e-10, /normal
	
    	result = fit_psd(pfreq, power, pspecerr=pspecerr)
    	pvalue = result[4]
    	pfsim = interpol([pfreq[0], pfreq[-1]], 100)
    	powsim = result[0] + result[1]*pfsim
    	rsusMm = 696.34 ; Mm
	turbindex = result[1]

	;----------------------------;
        ;    Plor 5/3 and 7/3 PSDs.
        ;
        powturb = result[0]-0.0 + (-5/3.)*pfsim
        oplot, 10^pfsim, 10^powturb, linestyle=5, color=7, thick=8
	oplot, 10^pfsim, 10^powsim, color=5, thick=4

    	;----------------------------;
    	;
    	;    Plot 99% confidence
    	;
    	meansig = 10^sigcutoff
    	print, '99% confidence thresh: '+string(meansig)
    	oplot, [10^wavenum0, 10^wavenum1], [meansig, meansig], linestyle=5, color=0


    	sindfit = string(round(result[1]*100.0)/100., format='(f6.2)')
    	alpha = cgsymbol('alpha')
    	legend,[alpha+':'+sindfit, alpha+'!L5/3!N'], linestyle=[0, 5], color=[5, 7], $
            box=0, /top, /right, thick=[4, 4]

    	axis, xaxis=0, xr = [10.0^wavenum0, 10.0^wavenum1], /xlog, /xs, xtitle='Wavenumber (R!U-1!N)'
    	axis, xaxis=1, xr = [10.0^wavenum0/rsunMm, 10^wavenum1/rsunMm], /xlog, /xs, xtitle='(Mm!U-1!N)'

	;-------------------------------;
	;	Fit burst
	;	
	window, 1, xs=500, ys=500
	tburst = [tburst[20], tburst[-20]]
	rads = [rads[20], rads[-20]]
	tsec = tburst - tburst[0]
	plot, tsec, rads, xtitle='Time (s)', ytitle='Heliocentric distance', /xs, /ys
	result = linfit(tsec, rads)
	ymodel = result[0] + result[1]*tsec
	set_line_color
	oplot, tsec, ymodel, color=5
	speed = result[1]*695.0e6/2.98e8
	print, 'Speed (c)'+string(speed)

	hbspeeds = [hbspeeds, speed]
	hbalphas = [hbalphas, turbindex]
        
	endfor	

	stop	
END
