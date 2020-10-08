pro compute_all_psds, data, utimes, freq, $
	sindices=sindices, stimes=stimes, pfreqs=pfreqs, powers=powers, $
	pspecerr=pspecerr, sigcuts=sigcuts, psdsmooth=psdsmooth, pvalcutoff=pvalcutoff, component=component

	
	;-----------------------------------------;
        ;       Each profile is evenly sampled
        ;       in f, but unevenly in space.
        ;       This gets an even sample in space
        ;       by interpolation.
        npoints=((freq*1e6/component)/8980.0)^2.0
        rads = density_to_radius(npoints, model='newkirk')
        even_rads = interpol([rads[0], rads[-1]], n_elements(freq))
        nt=n_elements(data[*,0])-1
        def = even_rads[2]-even_rads[1]
        sindices = fltarr(nt+1)
        stimes = dblarr(nt+1)
        vsave=0
        loadct, 0

        set_line_color
        pspecerr = 0.05

        wavenum0 = 1.0+alog10(2.0*!pi)
        wavenum1 = 3.0	;2.5+alog10(2.0*!pi)
        for i=0, nt do begin
                prof = data[i, *]
                even_prof = interpol(prof, rads, even_rads)
                even_prof = even_prof/max(even_prof)

                power = FFT_PowerSpectrum(even_prof, def, FREQ=pfreq,$
                        /tukey, width=psdsmooth, sig_level=0.01, SIGNIFICANCE=signif)


                pfreq = alog10(pfreq*2.0*!pi) ; x 2pi to get wavenumber from 1/lambda
                power = alog10(power)
                ind0 = closest(pfreq, wavenum0 )
                ind1 = closest(pfreq, wavenum1 )
                pfreq = pfreq[ind0:ind1]
                power = power[ind0:ind1]
                sigcutoff = alog10(signif[0])

                result = fit_psd(pfreq, power, pspecerr=pspecerr)
                pvalue = result[4]

                ;print, ' '
                ;print, 'Reduced chi square value: ' + string(chisq)
                ;print, 'Prob random variables has better chi: '+ string(pvalue)+'%'
                ;print, '---'

                pfsim = interpol([pfreq[0], pfreq[-1]], 100)
                powsim = result[0] + result[1]*pfsim

		if pvalue gt pvalcutoff and result[1] lt -1.0 then begin

                        if keyword_set(plot_ipsd) then begin
                        plot, pfreq, power, /xs, /ys, ytitle='log!L10!N(PSD Rs!U-1!N)', $
                                xtitle='log!L10!N(k Rs!U-1!N)', $
                                title=anytim(utimes[i], /cc)+'  S:'+string(result[1], format='(f5.2)'), $
                                yr=[-6, -2];, /noerase, color=colors[i], psym=1

                        xerr = dblarr(n_elements(pfreq))
                        yerr = pspecerr*abs(power) ;replicate(0.1, n_elements(pfreq))
                        oploterror, pfreq, power, xerr, yerr
                        oplot, pfsim, powsim, color=10
                        oplot, [wavenum0, wavenum1], [sigcutoff, sigcutoff], color=200
                        ;wait, 0.001
                        stop
                        xyouts, 0.6, 0.8, pvalue, /normal
                        endif

                        sindices[i] = result[1]
                        stimes[i] = utimes[i]
                        if vsave eq 0 then begin
                                powers = [power]
                                pfreqs = [pfreq]
                                sigcuts = sigcutoff
                                vsave = 1
                        endif else begin
                                powers = [ [powers], [[power]] ]
                                pfreqs = [ [pfreqs], [[pfreq]] ]
                                sigcuts = [sigcuts, sigcutoff]
                        endelse
                endif
        endfor
			
END
