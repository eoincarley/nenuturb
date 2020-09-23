function plot_mean_psd, powers, pfreqs, pspecerr, sigcuts


        ;setup_ps, './eps/nfar_mean_PSD_lin_typeIIc.eps', xsize=7, ysize=7

        mp = mean(powers, dim=2)
        mf = mean(pfreqs, dim=2)
        p = fit_psd(mf, mp, pspecerr=pspecerr)
        aerr = p[2]
        ierr = p[3]

        pfsim = interpol([mf[0], mf[-1]], 100)
        powsim = p[0] + p[1]*pfsim
        set_line_color

        ;--------------------------;
        ;      Plot mean PSD
        ;
	mf = 10^mf
        mp = 10^mp
        powsim = 10^powsim
        pfsim = 10^pfsim

        plot, mf, mp, /xs, /ys, ytitle='PSD', $
              pos = [0.15, 0.15, 0.9, 0.9], /noerase, thick=5, XTICKINTERVAL=0.5, /ylog, /xlog, $
              xr=[10.0, 10.0^2.5], XTICKFORMAT="(A1)", xticklen=1e-10

        axis, xaxis=0, xr = [10.0, 10.0^2.5], /xlog, /xs, xtitle='Wavenumber (R!U-1!N)'
        axis, xaxis=1, xr = [10.0/696.34, 10^2.5/696.34], /xlog, /xs, xtitle='(Mm!U-1!N)'

        oplot, pfsim, powsim, color=5, thick=8
	
	;----------------------------;
        ;    Plor 5/3 and 7/3 PSDs.
        ;
        powturb = p[0]-0.45 + (-5/3.)*alog10(pfsim)
        oplot, pfsim, 10^powturb, linestyle=5, color=7, thick=8

        powturb = p[0]+0.2 + (-7/3.)*alog10(pfsim)
        oplot, pfsim, 10^powturb, linestyle=5, color=6, thick=8
	
	;----------------------------;
	;    Plot 99% confidence 
	;
        meansig = 10^mean(sigcuts)
        print, '99% confidence thresh: '+string(meansig)
	oplot, [pfsim[0], pfsim[-1]], [meansig, meansig], linestyle=5, color=1
	
	;----------------------------;
        ;    Plot errors 
        ;
        alpha = cgsymbol('alpha')
        aerrstr = string(round(aerr*100.0)/100., format='(f4.2)')
        sindfit = string(round(p[1]*100.0)/100., format='(f6.2)')
        legend,[alpha+':'+sindfit+'+/-'+aerrstr, alpha+'!L5/3!N', alpha+'!L7/3!N'], linestyle=[0,5,5], color=[5, 7, 6], $
                box=0, /top, /right, charsize=1.6, thick=[4,4,4]

        loadct, 0
        powturb = p[0]+ierr + (p[1]+aerr)*alog10(pfsim)
        oplot, pfsim, 10^powturb, linestyle=1, color=50, thick=4
        powturb = p[0]-ierr + (p[1]-aerr)*alog10(pfsim)
        oplot, pfsim, 10^powturb, linestyle=1, color=50, thick=4



        ;axis, xaxis=1, xtickv = [10.0/696.34, 100.0/696.34], /xs
        ;device, /close
        ;set_plot, 'x'  

end
