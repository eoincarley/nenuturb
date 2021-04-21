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

pro plot_spectro, data, time, freq, f0, f1, posit

        loadct, 0
        utplot, time, freq, yr=[f1,f0], /xs, /ys, xtitle='Time (UT)', ytitle='Frequency (MHz)', $
                title=' ', pos=posit, /normal, color=150, $
                xr=[time[0], time[-1]], charsize=2
        ;------------------------------------------;
        ;            Plot spectrogram
        ;
        loadct, 74
        reverse_ct
	
	;data = alog10(data)
        spectro_plot, data>1e5<15e6, time, freq, /xs, /ys, $
                ytitle=' ', xtitle=' ', yr=[f1, f0], /noerase, XTICKFORMAT="(A1)", YTICKFORMAT="(A1)", $
                position=posit, /normal, xr=[time[0], time[-1]]


END

pro plot_powerspec, freq, prof, pos = pos, ytickfmt = ytickfmt, ytitle=ytitle, xtitle=xtitle, xtiupper=xtiupper

	;---------------------------------------;
        ;       Now perform PSD on iburst
        ;       
        ;-----------------------------------------;
        ;       Each profile is evenly sampled
        ;       in f, but unevenly in space.
        ;       This gets an even sample in space
        ;       by interpolation.
        npoints=((freq*1e6/1.)/8980.0)^2.0
        rads = density_to_radius(npoints, model='saito')
        even_rads = interpol([rads[0], rads[-1]], n_elements(freq))
        ;nt = n_elements(data[*,0])-1
        def = even_rads[2]-even_rads[1]
        pspecerr = 0.05
        wavenum0 = 1.0+alog10(2.0*!pi)
        wavenum1 = 3.5+alog10(2.0*!pi)
        rsunMm = 696.34 ; Mm

        loadct, 0
        ;----------------------------------------------;
        ;
        ;  Get evenly sampled in space and perform PSD
        ;
        even_prof = interpol(prof, rads, even_rads)
        even_prof = even_prof/max(even_prof)

        ;plot, even_rads, even_prof, /xs, /ys, pos=[0.48, 0.15, 0.7, 0.9], /normal, /noerase, $
        ;       xtitle=' ', ytitle='Intensity', XTICKFORMAT="(A1)", xticklen=-1e-8


        ;axis, xaxis=0, xr = [even_rads[0], even_rads[-1]], /xs, xtitle='Heliocentric distance (R!Ls!N)'
        ;axis, xaxis=1, xr = [even_rads[0], even_rads[-1]]*rsunMm, /xs, xtitle='(Mm)'

        power = FFT_PowerSpectrum(even_prof, def, FREQ=pfreq,$
                /tukey, width=0.002, sig_level=0.01, SIGNIFICANCE=signif)

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
        result = fit_psd(pfreq, power, pspecerr=pspecerr)
        pvalue = result[4]
        pfsim = interpol([pfreq[0], pfreq[-1]], 100)
        powsim = result[0] + result[1]*pfsim
        rsusMm = 696.34 ; Mm
	print, result[1]
        set_line_color
        plot, 10^pfreq, 10^power, /xlog, /ylog, /xs, /ys, ytitle=ytitle, $
            xtitle=' ', thick=2, xr = [10^wavenum0, 10^wavenum1], yr=10.0^[-7.0, -2.0], $
            /noerase, psym=10, ytickformat=ytickfmt, $
            XTICKFORMAT="(A1)", xticklen=1e-10, /normal, $
	    pos = pos 

        ;----------------------------;
        ;    Plor 5/3 and 7/3 PSDs.
        ;
        powturb = result[0]-0.5 + (-5/3.)*pfsim
        oplot, 10^pfsim, 10^powturb, linestyle=5, color=7, thick=8
        oplot, 10^pfsim, 10^powsim, color=5, thick=4

        ;----------------------------;
        ;
        ;    Plot 99% confidence
        ;
        meansig = 10^sigcutoff
        print, '99% confidence thresh: '+string(meansig)
        oplot, [10^wavenum0, 10^wavenum1], [meansig, meansig], linestyle=5, color=0

        ;xerr = dblarr(n_elements(pfreq))
        ;yerr = pspecerr*abs(power) ;replicate(0.1, n_elements(pfreq))
        ;oploterror, pfreq, power, xerr, yerr

        sindfit = string(round(result[1]*100.0)/100., format='(f6.2)')
        alpha = cgsymbol('alpha')
        ;legend,[alpha+':'+sindfit, alpha+'!L5/3!N'], linestyle=[0, 5], color=[5, 7], $
        ;   box=0, /top, /right, thick=[4, 4]

        axis, xaxis=0, xr = [10.0^wavenum0, 10.0^wavenum1], /xlog, /xs, xtitle=xtitle, xtickname=['10!U1!N', '10!U2!N', '10!U3!N']
        axis, xaxis=1, xr = [10.0^wavenum0/rsunMm, 10^wavenum1/rsunMm], /xlog, /xs, xtitle=xtiupper, xtickname=[' ', '10!U0!N', '10!U1!N']



END

pro typeIIc_drifter_plots, postscript=postscript

	if keyword_set(postscript) then begin
                setup_ps, './eps/psd_typeIc_drifters.eps', xsize=6, ysize=9
		!p.charsize=1.8
        endif else begin
                !p.charsize=1.8
                window, xs=600, ys=900
        endelse

	posit = [0.15, 0.6, 0.95, 0.94]
	restore,'./typeIIc_drifters/typeIIc_drifter_spectrogram.sav'
	plot_spectro, data, utimes, freq, f0, f1, posit

	;---------------------------------;
	; 	
	;        First burst
	;
	restore,'./typeIIc_drifters/typeIIc_drifter_20190320_1132_1.sav', /verb
	tburst0 = tburst
	fburst0 = fburst
	iburst0 = iburst
	set_line_color
        plots, tburst, fburst, /data, psym=1, symsize=0.2, color=5
	
	;---------------------------------;
        ;
        ;        Second burst
        ;
        restore,'./typeIIc_drifters/typeIIc_drifter_20190320_1132_2.sav', /verb
        tburst1 = tburst
        fburst1 = fburst
	iburst1 = iburst
	set_line_color
        plots, tburst, fburst, /data, psym=1, symsize=0.1, color=5


	;---------------------------------;
        ;
        ;        Third burst
        ;
        restore,'./typeIIc_drifters/typeIIc_drifter_20190320_1133.sav', /verb
        tburst2 = tburst
        fburst2 = fburst
	iburst2 = iburst
	set_line_color
        plots, tburst, fburst, /data, psym=1, symsize=0.1, color=5

	plot_powerspec, fburst0, iburst0, pos = [0.15, 0.15, 0.4, 0.46], ytitle='PSD'
	plot_powerspec, fburst1, iburst1, pos = [0.42, 0.15, 0.67, 0.46], ytickfmt="(A1)", xtitle='Wavenumber (R!U-1!N)', xtiupper='(Mm!U-1!N)'
	plot_powerspec, fburst2, iburst2, pos = [0.69, 0.15, 0.94, 0.46], ytickfmt="(A1)"
		
	if keyword_set(postscript) then device, /close
        set_plot, 'x'	
	
	
END
