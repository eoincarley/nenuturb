pro setup_ps, name
  
   set_plot,'ps'
   !p.font=0
   !p.charsize=1.5
   device, filename = name, $
          /color, $
          /helvetica, $
          /inches, $
          xsize=10, $
          ysize=14, $
          /encapsulate, $
          yoffset=5, $
          bits_per_pixel = 16

end

pro read_nfar_data, file, t0, t1, f0, f1, data=data, utimes=utimes, freq=freq
   	READ_NU_SPEC, file, data,time,freq,beam,ndata,nt,dt,nf,df,ns, $
                tmin=t0*60.0, tmax=t1*60.0, fmin=f0, fmax=f1
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
                psym=1, symsize=0.5, color=5, xr=[utimes[0], utimes[-1]], $
                xtitle='Time (UT)', ytitle='PSD spectral index ('+alpha+')'
        outplot, [utimes[0], utimes[-1]], [sturb0, sturb0], linestyle=5, thick=4, color=7
        outplot, [utimes[0], utimes[-1]], [sturb1, sturb1], linestyle=5, thick=4, color=6
        outplot, utimes, sindices, psym=1, symsize=0.5, color=5
        legend,[alpha+'!L5/3!N', alpha+'!L7/3!N'], linestyle=[5, 6], color=[7, 6], $
                box=0, /top, /right, charsize=1.4, thick=[4, 6]


end

function fit_psd, frequency, power, signif
	
	start = [-1, -1]
        fit = 'p[0] + p[1]*x'
        err = power
        err[*] = 0.3
        p = mpfitexpr(fit, frequency, power, err, start, perror = perror, $
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

function plot_mean_psd, powers, pfreqs

	mp = mean(powers, dim=2)
        mf = mean(pfreqs, dim=2)
	p = fit_psd(mf, mp)
	aerr = p[2]
	ierr = p[3]

        pfsim = interpol([mf[0], mf[-1]], 100)
        powsim = p[0] + p[1]*pfsim
        set_line_color
        plot, mf, mp, /xs, /ys, ytitle='log!L10!N(PSD)', xtitle='log!L10!N(k) R!U-1!N', $
              pos = [0.12, 0.23, 0.48, 0.47], /noerase, thick=5
        oplot, pfsim, powsim, color=5, thick=6

        powturb = p[0]-0.35 + (-5/3.)*pfsim
        oplot, pfsim, powturb, linestyle=5, color=7, thick=4

	powturb = p[0] + (-7/3.)*pfsim
        oplot, pfsim, powturb, linestyle=5, color=6, thick=4


	alpha = cgsymbol('alpha')
        aerrstr = string(round(aerr*100.0)/100., format='(f4.2)')
        sindfit = string(round(p[1]*100.0)/100., format='(f6.2)')
        legend,[alpha+':'+sindfit+'+/-'+aerrstr, alpha+'!L5/3!N', alpha+'!L7/3!N'], linestyle=[0,5,5], color=[5, 7, 6], $
                box=0, /top, /right, charsize=1.4, thick=[4,4,4]

        loadct, 0
        powturb = p[0]+ierr + (p[1]+aerr)*(pfsim)
        oplot, pfsim, powturb, linestyle=5, color=50, thick=2
        powturb = p[0]-ierr + (p[1]-aerr)*(pfsim)
        oplot, pfsim, powturb, linestyle=5, color=50, thick=2
end

function plot_alpha_hist, sindices

	alpha = cgsymbol('alpha')
	set_line_color
        plothist, sindices, bin=0.025, $
                xtitle='PSD spectral index '+alpha, ytitle='Count', $
                pos = [0.59, 0.23, 0.95, 0.47], /noerase, color=0, yr=[0, 250], thick=4
        meanalpha = string(round(median(sindices)*100.)/100.0, format='(f5.2)')
        oplot, [meanalpha, meanalpha], [0, 600.0], color=5, thick=5
        oplot, [-1.66, -1.66], [0, 600.0], color=7, thick=4, linestyle=5
	oplot, [-2.33, -2.33], [0, 600.0], color=6, thick=4, linestyle=5

        legend,[cgsymbol('mu')+'!L'+alpha+'!N: '+meanalpha, alpha+'!L5/3!N', alpha+'!L7/3!N'], linestyle=[0, 5, 5], color=[5, 7, 6], $
                box=0, /top, /left, charsize=1.4, thick=[5,4,4]

end


pro psd_typeIIa_lin, save=save, postscript=postscript, rebin=rebin

	; PSD of first type II. Code working.

	path = '/databf2/nenufar-tf/ES11/2019/03/20190320_104900_20190320_125000_SUN_TRACKING_BHR/'
        file = 'SUN_TRACKING_20190320_104936_0.spectra'

	t0 = 32.8
	t1 = 34.2	
	f0 = 21.0
	f1 = 24.5
	read_nfar_data, path+file, t0, t1, f0, f1, data=data, utimes=utimes, freq=freq
	   

	if keyword_set(postscript) then begin
		setup_ps, './eps/nfar_'+time2file(utimes[0])+'_PSD_lin.eps'
	endif else begin
		!p.charsize=1.8
		window, xs=800, ys=1200
	endelse	

	;data = 10.0*alog10(data)
	;data = constbacksub(data, /auto)
	;data = smooth(data,3)

	if keyword_set(rebin) then begin	
		nfbin = (size(data))[2]
		data = data[0:15999, *]
		data = rebin(data, tbin, nfbin)
		utimes = congrid(utimes, tbin)
		ntsteps=1
	endif else begin
		ntsteps=10
	endelse	
		;stop
	;------------------------------------------;
	;	Empty template to get black ticks
	;
	posit=[0.12, 0.75, 0.95, 0.95]
	loadct, 0
        utplot, utimes, freq, yr=[f1,f0], /xs, /ys, xtitle=' ', ytitle='Frequency (MHz)', $
		title='NenuFAR-ES11 '+time2file(utimes[0], /date), pos=posit, /normal, color=150, $
		xr=[utimes[0], utimes[-1]], XTICKFORMAT="(A1)"
	
        ;------------------------------------------;
	;	     Plot spectrogram
	;
	loadct, 74
        reverse_ct
	spectro_plot, sigrange(data), utimes, freq, /xs, /ys, $
		ytitle=' ', xtitle=' ', yr=[f1, f0], /noerase, XTICKFORMAT="(A1)", YTICKFORMAT="(A1)", $
		position=posit, /normal, xr=[utimes[0], utimes[-1]]
	;-----------------------------------------;
	;	         Flatten	
	;
	;data = congrid(data, 1e4, 1984)
	;utimes = congrid(utimes, 1e4)

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
	sindices = fltarr(nt+1)
	stimes = dblarr(nt+1)
	vsave=0
	loadct, 0
	window, 1, xs=600, ys=600	
	set_line_color
	for i=0, nt, ntsteps do begin
		prof = data[i, *]
		even_prof = interpol(prof, rads, even_rads)
		even_prof = even_prof/max(even_prof)

		power = FFT_PowerSpectrum(even_prof, def, FREQ=pfreq, /tukey, width=0.002, SIGNIFICANCE=signif)

		pfreq = alog10(pfreq)
		power = alog10(power)
		ind0 = closest(pfreq, 1.0)
		ind1 = closest(pfreq, 2.5)
		pfreq = pfreq[ind0:ind1]
		power = power[ind0:ind1]
		signif = signif[ind0:ind1]
		;wset, 1
		;plot, pfreq, power, /xs, /ys, ytitle='log!L10!N(Power Rs!U-1!N)', $
		;xtitle='log!L10!N(k Rs!U-1!N)';, yr=[1e8, 1e12]

		result = fit_psd(pfreq, power, signif)
		pvalue = result[4]
		
		;print, ' ' 
                ;print, 'Reduced chi square value: ' + string(chisq)
                ;print, 'Prob random variables has better chi: '+ string(pvalue)+'%'
                ;print, '---'    
	
		pfsim = interpol([pfreq[0], pfreq[-1]], 100)
		powsim = result[0] + result[1]*pfsim

		if pvalue gt 1.0 then begin	
			
			plot, pfreq, power, /xs, /ys, ytitle='log!L10!N(PSD Rs!U-1!N)', xtitle='log!L10!N(k Rs!U-1!N)', $
                        	title=anytim(utimes[i], /cc)+'  S:'+string(result[1], format='(f5.2)');, yr=[1e8, 1e12]
                	oplot, pfsim, powsim, color=10
                	xyouts, 0.6, 0.8, pvalue, /normal
		
			sindices[i] = result[1]
			stimes[i] = utimes[i]
			if vsave eq 0 then begin
				powers = [power]
				pfreqs = [pfreq]
				vsave = 1
			endif else begin
				powers = [ [powers], [[power]] ]
				pfreqs = [ [pfreqs], [[pfreq]] ]
			endelse	
		endif 	
	endfor

	sindices = sindices[where(sindices ne 0)]	
	stimes = stimes[where(stimes ne 0)]
	wset, 0
	;-----------------------------------;
	;
        ;       Plot alpha time series
        ;
	result = plot_alpha_time(stimes, sindices)


	;-----------------------------------;
	;
	;	Plot mean PSD
	;
	result = plot_mean_psd(powers, pfreqs)


	;-----------------------------------;
	;   Plot hist of spectral indice
	;
	result = plot_alpha_hist(sindices)

	if keyword_set(postscript) then begin
		device, /close
		set_plot, 'x'
	endif	
	  

stop
END
