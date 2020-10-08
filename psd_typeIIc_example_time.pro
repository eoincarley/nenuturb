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

; /fflat or fflat=1 = division of the output dynamic spectrum by the Stokes I average spectrum
;           fflat=2 = division of the output dynamic spectrum by the Stokes I median spectrum
;           fflat=3 = division of the output dynamic spectrum by the Stokes I empirical bandpass correction     => normalisation spectrum stored in corrf
; /tflat or tflat=1 = division of the output dynamic spectrum by the Stokes I average time profile
;           tflat=2 = division of the output dynamic spectrum by the Stokes I median time profile
;           tflat=3 = division of the output dynamic spectrum by the Stokes I 6-min profile computed at output timescale => gain variation stored in corrt
;           tflat=4 = same as tflat=3 + blanking of 3 spectra around 6-min gain jumps
;           tflat=5 = on-the-fly division of the dynamic spectrum by the Stokes I 6-min profile computed in a first reading of the data and stored in corrt (with abscissa in t


   	READ_NU_SPEC, file, data,time,freq,beam,ndata,nt,dt,nf,df,ns, $
                tmin=t0*60.0, tmax=t1*60.0, fmin=f0, fmax=f1, fflat=1, ntimes=4;, fclean=6
        utimes=anytim(file2time(file), /utim) + time
        data = reverse(data, 2)
        freq = reverse(freq)
end

function plot_alpha_time, utimes, sindices
	
	set_line_color
        sturb0 = -5/3.
        sturb1 = -7/3.
        alpha = cgsymbol('alpha')
        utplot, utimes, sindices, pos=[0.12, 0.54, 0.95, 0.74], $
                /noerase, /xs, /ys, yr=[-3, -1.0], $
                psym=1, symsize=0.8, color=5, xr=[utimes[0], utimes[-1]], $
                xtitle='Time (UT)', ytitle='PSD spectral index ('+alpha+')'

        outplot, [utimes[0], utimes[-1]], [sturb0, sturb0], linestyle=5, thick=4, color=7
        outplot, [utimes[0], utimes[-1]], [sturb1, sturb1], linestyle=0, thick=4, color=6
        outplot, utimes, sindices, psym=1, symsize=0.5, color=5
        legend,[alpha+'!L5/3!N', alpha+'!L7/3!N'], linestyle=[5, 0], color=[7, 6], $
                box=0, /top, /right, charsize=1.4, thick=[4, 6]


end

function fit_psd, frequency, power, pspecerr=pspecerr
	
	start = [-1, -1]
        fit = 'p[0] + p[1]*x'
        err = power
        err[*] = pspecerr*abs(power)

	; Note here that if I supply the weights keyword with poisson weighting the fit
	; gives a p-value of 1 for eveything, e.g. perfect fit within error. Every fit
	; is then accepted below, even a few power spectra that are not power law.
	; If I give an errors (without the weights) then not everything is accepted. Some
	; fits are rejected. However, the choice of err here is slightly arbitrary.

        p = mpfitexpr(fit, frequency, power, err, weights=1/power^2, start, perror = perror, $
		yfit=yfit, bestnorm=bestnorm, dof=dof, /quiet)
        perror = perror * SQRT(BESTNORM / DOF)
        aerr = perror[1]*2.0 ; 2-sigma uncertainty on the slope
        ierr = perror[0]*2.0 ; 2-sigma uncertainty on the intercept

	;----------------------------------------;
	; This is the probability the chi^2 score 
	; is worse than calculated value.
	; If <5% reject the fit
	pvalue = (1.0-CHISQR_PDF(bestnorm, DOF))*100.0	

	return, [p[0], p[1], aerr, ierr, pvalue]

end


function plot_alpha_hist, sindices

	alpha = cgsymbol('alpha')
	set_line_color
        plothist, sindices, bin=0.025, $
                xtitle='PSD spectral index '+alpha, ytitle='Count', $
                pos = [0.59, 0.18, 0.95, 0.42], /noerase, color=0, yr=[0, 250], thick=4
        meanalpha = string(round(median(sindices)*100.)/100.0, format='(f5.2)')
        oplot, [meanalpha, meanalpha], [0, 600.0], color=5, thick=5
        oplot, [-1.66, -1.66], [0, 600.0], color=7, thick=4, linestyle=5
	oplot, [-2.33, -2.33], [0, 600.0], color=6, thick=4, linestyle=5

        legend,[cgsymbol('mu')+'!L'+alpha+'!N: '+meanalpha, alpha+'!L5/3!N', alpha+'!L7/3!N'], linestyle=[0, 5, 5], color=[5, 7, 6], $
                box=0, /top, /left, charsize=1.4, thick=[5,4,4]
 	
end


function apply_response, data, freq

	restore, 'nfar_response.sav'
	ind0 = where(rfreq eq freq[0])
	ind1 = where(rfreq eq freq[-1])
	response = response[ind0:ind1]
	
	for i=0, n_elements(data[*, 0])-1 do begin
		data[i, *] = data[i, *]/response
	endfor	
	
	return, data

END


pro psd_typeIIc_example_time, save=save, plot_ipsd=plot_ipsd, postscript=postscript, rebin=rebin

	; PSD of first type II. Code working.

	; This code take a single sample in type IIc as an example of how the PSD is performed for one time.

	; This one step through time to find a good example.

	; psd_typeIIc_example.pro plots the good example.

	path = '/databf2/nenufar-tf/ES11/2019/03/20190320_104900_20190320_125000_SUN_TRACKING_BHR/'
        file = 'SUN_TRACKING_20190320_104936_0.spectra'

	t0 = 42.0
	t1 = 43.0	
	f0 = 21.0
	f1 = 33.0
	read_nfar_data, path+file, t0, t1, f0, f1, data=data, utimes=utimes, freq=freq
	   

	if keyword_set(postscript) then begin
		setup_ps, './eps/nfar_PSD_lin_typeIIc.eps', xsize=10, ysize=14
	endif else begin
		!p.charsize=1.8
		window, xs=1200, ys=500
	endelse	
	

	;------------------------------------------;
	;	Empty template to get black ticks
	;
	posit=[0.08, 0.15, 0.28, 0.9]
	loadct, 0
        utplot, utimes, freq, yr=[f1,f0], /xs, /ys, xtitle='Time (UT)', ytitle='Frequency (MHz)', $
		title='NenuFAR-ES11 '+time2file(utimes[0], /date), pos=posit, /normal, color=150, $
		xr=[utimes[0], utimes[-1]], tick_unit=20.0
	
        ;------------------------------------------;
	;
	;	     Plot spectrogram
	;
	loadct, 74
        reverse_ct
	spectro_plot, sigrange(data), utimes, freq, /xs, /ys, $
		ytitle=' ', xtitle=' ', yr=[f1, f0], /noerase, XTICKFORMAT="(A1)", YTICKFORMAT="(A1)", $
		position=posit, /normal, xr=[utimes[0], utimes[-1]]

	tsample = anytim('2019-03-20T11:31:52', /utim)
	outplot, [tsample, tsample], [f0, f1], linestyle=0, thick=3


	;-----------------------------------------;
	;	Each profile is evenly sampled
	;	in f, but unevenly in space.
	; 	This gets an even sample in space
	;	by interpolation.	
	npoints=((freq*1e6/2.)/8980.0)^2.0 	
	rads = density_to_radius(npoints, model='newkirk')
	even_rads = interpol([rads[0], rads[-1]], n_elements(freq))
	nt=n_elements(data[*,0])-1
	def = even_rads[2]-even_rads[1]
	loadct, 0
	pspecerr = 0.05
	wavenum0 = 1.0+alog10(2.0*!pi)
        wavenum1 = 2.5+alog10(2.0*!pi)
        rsunMm = 696.34 ; Mm	
	
	;----------------------------------------;
	;	Get profile and plot.
	;	
	for tindex=0, n_elements(utimes)-1 do begin
		wset, 0
		;tindex = (where(utimes ge tsample))[0]
		prof = data[tindex, *]
	
		plot, freq, prof, /xs, /ys, pos=[0.35, 0.15, 0.68, 0.9], /normal, $
			xtitle='Frequency (MHz)', ytitle='Intensity', title=anytim(utimes[tindex], /cc)
	
		;----------------------------------------------;
		;  Get evenly sampled in space and perform PSD
		;	
		even_prof = interpol(prof, rads, even_rads)
		even_prof = even_prof/max(even_prof)

		power = FFT_PowerSpectrum(even_prof, def, FREQ=pfreq,$ 
			/tukey, width=0.001, sig_level=0.01, SIGNIFICANCE=signif)

		pfreq = alog10(pfreq*2.0*!pi)
        	power = alog10(power)
        	ind0 = closest(pfreq, wavenum0)
        	ind1 = closest(pfreq, wavenum1)
		pfreq = pfreq[ind0:ind1]
		power = power[ind0:ind1]
		sigcutoff = alog10(signif[0])
		
		;----------------------------------;
		;	Fit PSD and plot
		; 
		result = fit_psd(pfreq, power, pspecerr=pspecerr)
		pvalue = result[4]
		
	        ;print, 'Reduced chi square value: ' + string(chisq)
	        ;print, 'Prob random variables has better chi: '+ string(pvalue)+'%'
	
		pfsim = interpol([pfreq[0], pfreq[-1]], 100)
		powsim = result[0] + result[1]*pfsim
			        
		plot, 10^pfreq, 10^power, /xlog, /ylog, /xs, /ys, ytitle='PSD', $
                	xtitle=' ', thick=2, $
                	yr=10^[-6, -2], /noerase, position=[0.75, 0.15, 0.99, 0.9], psym=10, $
                	XTICKFORMAT="(A1)", xticklen=1e-10

	        axis, xaxis=0, xr = [10.0^wavenum0, 10.0^wavenum1], /xlog, /xs, xtitle='Wavenumber (R!U-1!N)'
	        axis, xaxis=1, xr = [10.0^wavenum0/rsunMm, 10^wavenum1/rsunMm], /xlog, /xs, xtitle='(Mm!U-1!N)'

		;----------------------------;
	        ;    Plot 99% confidence 
	        ;
	        set_line_color
	        meansig = 10^sigcutoff
	        ;print, '99% confidence thresh: '+string(meansig)
	        oplot, 10^[pfsim[0], pfsim[-1]], [meansig, meansig], linestyle=5, color=1

	        xerr = dblarr(n_elements(pfreq))
	        yerr = pspecerr*abs(power) ;replicate(0.1, n_elements(pfreq))

	        ;oploterror, pfreq, power, xerr, yerr
	        set_line_color
	        oplot, 10^pfsim, 10^powsim, color=5, thick=4
	        oplot, [10^wavenum0, 10^wavenum1], [sigcutoff, sigcutoff], color=1, linestyle=1

	        sindfit = string(round(result[1]*100.0)/100., format='(f6.2)')
	        alpha = cgsymbol('alpha')
	        legend,[alpha+':'+sindfit], linestyle=[0], color=[5], $
	                box=0, /top, /right, thick=[4]
	
		wait, 0.1
	;	if utimes[tindex] ge anytim('2019-03-20T11:32:17', /utim) then stop
	endfor
stop
END
