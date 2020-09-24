; Used in each of the psd_typeIIx_lin_v2 codes.

function plot_all_psd, pfreqs, powers, times, powrange=powrange

        ; Plot all spectra over time
        ntimes = n_elements(powers[0,*])
        colors = interpol([0,255], ntimes)
       
       	wavenum0=1.0+alog10(2.0*!pi)
	wavenum1=2.5+alog10(2.0*!pi)	
	loadct, 0
        plot, 10^[wavenum0, wavenum1], 10^powrange, xr=10^[wavenum0, wavenum1], /nodata, yr=10^powrange, $
		/xs, /ys, ytitle='PSD', /ylog, /xlog, $
        	xtitle='Wavenumber Rs!U-1!N', pos = [0.12, 0.18, 0.48, 0.42], /noerase

        loadct, 72
        reverse_ct
        pfq = 10^pfreqs
        pow = 10^powers
        for i=ntimes-1, 0, -1 do begin
                oplot, pfq[*, i], pow[*, i], psym=1, color=colors[i], symsize=0.3
        endfor

        trange = (times - times[0])/60.0
        cgCOLORBAR, range=[trange[0], trange[-1]],  POSITION=[0.12, 0.43, 0.48, 0.45], $
                title='Mins after '+anytim(times[0], /cc, /trun), /top, charsize=1.0

END
