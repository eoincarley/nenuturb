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
   	
	
	restore, '~/nenuturb/savs/calibration_factor.sav'
	freq=freq0 & corrf=corrf0
	
	;READ_NU_SPEC, file, data, time, freq, beam, ndata, nt, dt, nf, df, ns, tmin=t0*60.0, tmax=t1*60.0, ntimes=8, ex_chan=[0], /fill
	;data=data/rebin(reform(cor,1,nf),nt,nf)
	;w=where(freq ge f0 and freq le f1)
	;data=data[*,w]
	;freq = freq[w]

	;READ_NU_SPEC, file, data, time, freq, beam, ndata, nt, dt, nf, df, ns, $
	;	jd0, h0, corrf, tmin=t0*60.0, tmax=t1*60.0, ntimes=8, nchannels=512, $
	;	fmin=f0, fmax=f1, /exactfreq, fflat=4, /fill
	READ_NU_SPEC, file, data, time, freq, beam, ndata, nt, dt, nf, df, ns, $
                tmin=t0*60.0, tmax=t1*60.0, fmin=f0, fmax=f1, fflat=3, ntimes=8, ex_chan=[0], /fill, /exactfreq
        
	utimes=anytim(file2time(file), /utim) + time
        data = reverse(data, 2)
        freq = reverse(freq)

end

function remove_spikes, prof

	; Harcode removal of spikes, no fancy business.
	endprof = prof[800:-1]
	ind = where(endprof gt 3e7)
	endprof[ind] = median(endprof)
	prof[800:-1] = endprof
	return, prof

END

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

pro psd_typeIIa_hb_avg, save=save, plot_ipsd=plot_ipsd, postscript=postscript, rebin=rebin

	; PSD of herringbone in type IIa	
	; More appropriate when structures actually drift in frequency time.

	path = '/volumes/plasma/nenufar/nenufar-tf/ES11/2019/03/20190320_104900_20190320_125000_SUN_TRACKING_BHR/'
  file = 'SUN_TRACKING_20190320_104936_0.spectra'

	
	t0 = 33.61
  t1 = 33.64
  f0 = 33.0
  f1 = 41.0
	read_nfar_data, path+file, t0, t1, f0, f1, data=data, utimes=utimes, freq=freq
	   

	if keyword_set(postscript) then begin
		setup_ps, './eps/psd_hbone_drift_H.eps', xsize=18, ysize=5.5
	endif else begin
		!p.charsize=1.8
		window, xs=1400, ys=600 ;xs=1600, ys=600
	endelse	

	posit=[0.05, 0.15, 0.42, 0.9]
  ;------------------------------------------;
	;	     Plot spectrogram
	;
  plot_spectro, data, utimes, freq, f0, f1, posit

  for i=0, n_elements(utimes)-1 do begin

  prof = transpose(data[i, *])

	prof = remove_spikes(prof)
	
	;---------------------------------------;
	;      	Now perform PSD on iburst
	;	
	;-----------------------------------------;
  ;       Each profile is evenly sampled
  ;       in f, but unevenly in space.
  ;       This gets an even sample in space
  ;       by interpolation.
  npoints=((freq*1e6/2.)/8980.0)^2.0
  rads = density_to_radius(npoints, model='newkirk')
  even_rads = interpol([rads[0], rads[-1]], n_elements(freq))
  nt = n_elements(data[*,0])-1
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
  plot, even_rads, even_prof, /xs, /ys, pos=[0.48, 0.15, 0.7, 0.9], /normal, /noerase, $
      xtitle=' ', ytitle='Intensity', XTICKFORMAT="(A1)", xticklen=-1e-8;, xtickv=[2.1, 2.15, 2.2]


  axis, xaxis=0, xr = [even_rads[0], even_rads[-1]], xticks=4, xtickv=[2.1, 2.15, 2.2], xminor=2, /xs, xtitle='Heliocentric distance (R!Ls!N)'
  axis, xaxis=1, xr = [even_rads[0], even_rads[-1]]*rsunMm, /xs, xtitle='(Mm)'

  power = FFT_PowerSpectrum(even_prof, def, FREQ=pfreq, $
        /tukey, width=0.001, sig_level=0.01, SIGNIFICANCE=signif)

  if i eq 0 then powers = [power] else powers = [[powers], [power]]
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

	set_line_color
  plot, 10^pfreq, 10^power, /xlog, /ylog, /xs, /ys, ytitle='PSD', $
            xtitle=' ', thick=2, xr = 10^[wavenum0,wavenum1], yr=10.0^[-7.0, -2.0], $
	        /noerase, position=[0.76, 0.15, 0.96, 0.9], psym=10, $
            XTICKFORMAT="(A1)", xticklen=1e-10, /normal

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
  oplot, 10^[wavenum0, wavenum1], [meansig, meansig], linestyle=5, color=0

    	;xerr = dblarr(n_elements(pfreq))
    	;yerr = pspecerr*abs(power) ;replicate(0.1, n_elements(pfreq))
    	;oploterror, pfreq, power, xerr, yerr

  sindfit = string(round(result[1]*100.0)/100., format='(f6.2)')
  alpha = cgsymbol('alpha')
  legend,[alpha+':'+sindfit, alpha+'!L5/3!N'], linestyle=[0, 5], color=[5, 7], $
            box=0, /top, /right, thick=[4, 4]

  axis, xaxis=0, xr = 10^[wavenum0, wavenum1], /xlog, /xs, xtitle='Wavenumber (R!U-1!N)'
  axis, xaxis=1, xr = [10^wavenum0/rsunMm, 10^wavenum1/rsunMm], /xlog, /xs, xtitle='(Mm!U-1!N)'

  if keyword_set(postscript) then device, /close
  ;set_plot, 'x'   	
	endfor

  ;--------------------------------------
  ;
  ;     Plot average spec.
  ;
  ;
  window, 1, xs=500, ys=500
  set_line_color

  power = mean(powers, dim=2)
  power = alog10(power)
  ind0 = closest(pfreq, wavenum0)
  ind1 = closest(pfreq, wavenum1)
  power = power[ind0:ind1]
  
  set_line_color
  plot, 10^pfreq, 10^power, /xlog, /ylog, /xs, /ys, ytitle='PSD', $
          thick=2, xr = 10^[wavenum0, wavenum1], $
          /noerase, psym=10, color=1, yr=[1e-7, 1e-2], xtitle='Wavenumber (1/R)'


  result = linfit(pfreq, power)
  pfsim = interpol([pfreq[0], pfreq[-1]], 100)
  powsim = result[0] + result[1]*pfsim

  ;----------------------------;
  ;    Plor 5/3 and 7/3 PSDs.
  ;
  powturb = result[0]+0.8 + (-5/3.)*pfsim
  oplot, 10^pfsim, 10^powturb, linestyle=5, color=7, thick=8
  ;oplot, 10^pfsim, 10^powsim, color=5, thick=4


stop	
END
