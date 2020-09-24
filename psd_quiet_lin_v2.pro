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
                tmin=t0*60.0, tmax=t1*60.0, fmin=f0, fmax=f1, fflat=3
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

        p = mpfitexpr(fit, frequency, power, err, weights=1/err^2, start, perror = perror, $
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


pro psd_quiet_lin_v2, save=save, plot_ipsd=plot_ipsd, postscript=postscript, rebin=rebin

	; PSD of first type II. Code working.

	; This version of the code excludes certain powerlaw fits based on their p-value.


	path = '/databf2/nenufar-tf/ES11/2019/03/20190320_104900_20190320_125000_SUN_TRACKING_BHR/'
        file = 'SUN_TRACKING_20190320_104936_0.spectra'

	t0 = 114.0
	t1 = 116.0	
	f0 = 21.0
	f1 = 33.0
	read_nfar_data, path+file, t0, t1, f0, f1, data=data, utimes=utimes, freq=freq
	   

	if keyword_set(postscript) then begin
		setup_ps, './eps/nfar_PSD_lin_typeIIc.eps', xsize=10, ysize=14
	endif else begin
		!p.charsize=1.8
		window, xs=800, ys=1200
	endelse	

	;data = 10.0*alog10(data)
	;data = constbacksub(data, /auto)
	;data = smooth(data,3)
	
	if keyword_set(rebin) then begin	
		nfbin = (size(data))[2]
		data = data[0:19999, *]
		utimes = utimes[0:19999]
		tbin = 5000
		data = rebin(data, tbin, nfbin)
		utimes = congrid(utimes, tbin)
		ntsteps=1
	endif else begin
		ntsteps=10
	endelse	
	
	;data = apply_response(data, freq)

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

	if ~keyword_set(postscript) then $
                window, 1, xs=600, ys=600
        ;-----------------------------------------;
        ;         Get all PSDs and fits
        ;
        compute_all_psds, data, utimes, freq, $
                sindices=sindices, stimes=stimes, pfreqs=pfreqs, powers=powers, $
                pspecerr=pspecerr, sigcuts=sigcuts, ntsteps=ntsteps, psdsmooth=0.001, pvalcutoff=0.0


	if ~keyword_set(postscript) then wset, 0
	sindices = sindices[where(sindices ne 0)]	
	stimes = stimes[where(stimes ne 0)]
	;-----------------------------------;
	;
        ;       Plot alpha time series
        ;
	result = plot_alpha_time(stimes, sindices)

	;-----------------------------------;
	;
	;   	Plot hist of spectral indices
	;
	result = plot_alpha_hist(sindices)

	;----------------------------------;
	;
	;	Plot all psd
	;
	result = plot_all_psd(pfreqs, powers, stimes, pow=[-8.0, -1.0])

	if keyword_set(postscript) then begin
		device, /close
		set_plot, 'x'
	endif	

	loadct, 0
	if ~keyword_set(postscript) then window, 1, xs=400, ys=400  
 	
	if keyword_set(postscript) then $
                 setup_ps, './eps/nfar_mean_PSD_lin_quiet_time.eps', xsize=5, ysize=5
	;-----------------------------------;
        ;
        ;       Plot mean PSD
        ;
        result = plot_mean_psd(powers, pfreqs, pspecerr, sigcuts)
 	
	if keyword_set(postscript) then begin
                device, /close
                set_plot, 'x'
        endif
	
	
	stop
END
