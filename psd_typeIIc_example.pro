pro setup_ps, name, xsize=xsize, ysize=ysize
  
   set_plot,'ps'
   !p.font=0
   !p.charsize=1.2
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
                tmin=t0*60.0, tmax=t1*60.0, fmin=f0, fmax=f1, fflat=1;, fclean=6
        utimes=anytim(file2time(file), /utim) + time
        data = reverse(data, 2)
        freq = reverse(freq)
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


pro psd_typeIIc_example, save=save, plot_ipsd=plot_ipsd, postscript=postscript, rebin=rebin

	; PSD of first type II. Code working.

	; This code take a single sample in type IIc as an example of how the PSD is performed for one time.


	path = '/databf2/nenufar-tf/ES11/2019/03/20190320_104900_20190320_125000_SUN_TRACKING_BHR/'
        file = 'SUN_TRACKING_20190320_104936_0.spectra'

	t0 = 42.0
	t1 = 43.0	
	f0 = 21.0
	f1 = 33.0
	read_nfar_data, path+file, t0, t1, f0, f1, data=data, utimes=utimes, freq=freq
	   

	if keyword_set(postscript) then begin
		setup_ps, './eps/nfar_PSD_typeIIc_example.eps', xsize=12, ysize=5
	endif else begin
		!p.charsize=1.8
		window, xs=1200, ys=500
	endelse	
	
	if keyword_set(rebin) then begin	
		nfbin = (size(data))[2]
		data = data[0:9999, *]
		utimes = utimes[0:9999]
		tbin = 2000
		data = rebin(data, tbin, nfbin)
		utimes = congrid(utimes, tbin)
		ntsteps=1
	endif else begin
		ntsteps=10
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

	tsample = anytim('2019-03-20T11:32:17.100', /utim)
	outplot, [tsample, tsample], [f0, f1], linestyle=0, thick=5
	;xyouts, 0.1, 0.85, 'a', charthick=4, /normal

	;-----------------------------------------;
	;	Each profile is evenly sampled
	;	in f, but unevenly in space.
	; 	This gets an even sample in space
	;	by interpolation.	
	npoints=((freq*1e6/1.)/8980.0)^2.0 	
	rads = density_to_radius(npoints, model='newkirk')
	even_rads = interpol([rads[0], rads[-1]], n_elements(freq))
	nt=n_elements(data[*,0])-1
	def = even_rads[2]-even_rads[1]
	loadct, 0
	pspecerr = 0.05
	
	;----------------------------------------;
	;	Get profile and plot.
	;	
	tindex = (where(utimes ge tsample))[0]
	prof = data[tindex, *]
	
	
	;----------------------------------------------;
	;  Get evenly sampled in space and perform PSD
	;	
	even_prof = interpol(prof, rads, even_rads)
	even_prof = even_prof/max(even_prof)	

        plot, even_rads, even_prof/1e7, /xs, /ys, pos=[0.37, 0.15, 0.68, 0.9], /normal, /noerase, $
                xtitle=' ', ytitle='Intensity', XTICKFORMAT="(A1)", xticklen=1e-10

	axis, xaxis=0, xr = [even_rads[0], even_rads[-1]], /xs, xtitle='Heliocentric distance (R!Ls!N)'
	axis, xaxis=1, xr = [even_rads[0], even_rads[-1]]*696.34, /xs, xtitle='(Mm)'

	power = FFT_PowerSpectrum(even_prof, def, FREQ=pfreq,$ 
		/tukey, width=0.001, sig_level=0.01, SIGNIFICANCE=signif)
	
	pfreq = alog10(pfreq)
	power = alog10(power)
	ind0 = closest(pfreq, 1.0)
	ind1 = closest(pfreq, 2.5)
	pfreq = pfreq[ind0:ind1]
	power = power[ind0:ind1]
	sigcutoff = alog10(signif[0])
	
	;xyouts, 0.4, 0.85, 'b', charthick=4, /normal

	;----------------------------------;
	;	Fit PSD and plot
	; 
	result = fit_psd(pfreq, power, pspecerr=pspecerr)
	pvalue = result[4]
		
        ;print, 'Reduced chi square value: ' + string(chisq)
        ;print, 'Prob random variables has better chi: '+ string(pvalue)+'%'
	
	pfsim = interpol([pfreq[0], pfreq[-1]], 100)
	powsim = result[0] + result[1]*pfsim

			
	plot, 10^pfreq*2.0*!pi, 10^power, /xlog, /ylog, /xs, /ys, ytitle='PSD', $
		xtitle=' ', thick=2, $
		yr=10^[-6, -2], /noerase, position=[0.75, 0.15, 0.99, 0.9], psym=10, $
		XTICKFORMAT="(A1)", xticklen=1e-10
                	
	axis, xaxis=0, xr = [10.0, 10.0^2.5]*2.0*!pi, /xlog, /xs, xtitle='Wavenumber (R!U-1!N)'
        axis, xaxis=1, xr = [10.0/696.34, 10^2.5/696.34]*2.0*!pi, /xlog, /xs, xtitle='(Mm!U-1!N)'
	
	;----------------------------;
        ;    Plot 99% confidence 
        ;
	set_line_color
        meansig = 10^sigcutoff
        print, '99% confidence thresh: '+string(meansig)
        oplot, 10^[pfsim[0], pfsim[-1]]*2.0*!pi, [meansig, meansig], linestyle=5, color=1


	xerr = dblarr(n_elements(pfreq))
	yerr = pspecerr*abs(power) ;replicate(0.1, n_elements(pfreq))
	
	;oploterror, pfreq, power, xerr, yerr
	set_line_color
	oplot, 10^pfsim*2.0*!pi, 10^powsim, color=5, thick=4
	oplot, 10^[1.0, 2.5], [sigcutoff, sigcutoff], color=1, linestyle=1
	
	sindfit = string(round(result[1]*100.0)/100., format='(f6.2)')
	alpha = cgsymbol('alpha')
	legend,[alpha+':'+sindfit], linestyle=[0], color=[5], $
                box=0, /top, /right, thick=[4]
	
	;xyouts, 0.79, 0.85, 'c', charthick=4, /normal

	if keyword_set(postscript) then device, /close
	set_plot, 'x'
	
stop
END
