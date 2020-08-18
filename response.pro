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
   	READ_NU_SPEC, file, data,time,freq,beam,ndata,nt,dt,nf,df,ns, $
                tmin=t0*60.0, tmax=t1*60.0, fmin=f0, fmax=f1
        utimes=anytim(file2time(file), /utim) + time
        data = reverse(data, 2)
        freq = reverse(freq)
end


pro response, points=points

	; This gets a spectrum from a quite time in the event and saves it to be used as a response 
	; curve to normalise the dynamic spectra before a PSD is taken.

	path = '/databf2/nenufar-tf/ES11/2019/03/20190320_104900_20190320_125000_SUN_TRACKING_BHR/'
        file = 'SUN_TRACKING_20190320_104936_0.spectra'

	t0 = 118.15
	t1 = 119.0	
	f0 = 21.0
	f1 = 55.0
	read_nfar_data, path+file, t0, t1, f0, f1, data=data, utimes=utimes, freq=freq
	   

	if keyword_set(postscript) then begin
		setup_ps, './eps/nfar_PSD_lin_typeIIc.eps', xsize=10, ysize=14
	endif else begin
		!p.charsize=1.8
		window, xs=400, ys=400
	endelse	

	;data = alog10(data)
	;data = constbacksub(data, /auto)
	;data = smooth(data,3)

	if keyword_set(rebin) then begin	
		nfbin = (size(data))[2]
		data = data[0:39999, *]
		tbin = 10000
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
	posit=[0.12, 0.12, 0.95, 0.95]
	loadct, 0
        utplot, utimes, freq, yr=[f1,f0], /xs, /ys, xtitle=' ', ytitle='Frequency (MHz)', $
		title='NenuFAR-ES11 '+time2file(utimes[0], /date), pos=posit, /normal, color=150, $
		xr=[utimes[0], utimes[-1]]
	
        ;------------------------------------------;
	;	     Plot spectrogram
	;
	loadct, 74
        reverse_ct
	spectro_plot, sigrange(data), utimes, freq, /xs, /ys, $
		ytitle=' ', xtitle=' ', yr=[f1, f0], /noerase, XTICKFORMAT="(A1)", YTICKFORMAT="(A1)", $
		position=posit, /normal, xr=[utimes[0], utimes[-1]]

	loadct, 0
	window, 2, xs=400, ys=400
	window, 1, xs=400, ys=400
	spectra = min(data[0:100, *], dim=1)
	
 	for i=0, n_elements(spectra)-1 do begin
		if spectra[i] gt 1.5e6 then spectra[i]=spectra[i-30]
	endfor	


	;-----------------------------------------;
	;
	; 	Firstly check the quiet time PSD
	;
	;-----------------------------------------;
        ;       Each profile is evenly sampled
        ;       in f, but unevenly in space.
        ;       This gets an even sample in space
        ;       by interpolation.       
        ;
	npoints=((freq*1e6/1.)/8980.0)^2.0
        rads = density_to_radius(npoints, model='newkirk')
        even_rads = interpol([rads[0], rads[-1]], n_elements(freq))
        nt=n_elements(data[*,0])-1
        def = even_rads[2]-even_rads[1]
        sindices = fltarr(nt+1)
        stimes = dblarr(nt+1)

	prof = mean(data[0:100, *], dim=1)
        even_prof = interpol(prof, rads, even_rads)
       	even_prof = even_prof/max(even_prof)
	power = FFT_PowerSpectrum(even_prof, def, FREQ=pfreq,$
                        /tukey, width=0.001, sig_level=0.01, SIGNIFICANCE=signif)
	
        pfreq = alog10(pfreq)
        power = alog10(power)
        ind0 = closest(pfreq, 1.0)
        ind1 = closest(pfreq, 2.5)
        pfreq = pfreq[ind0:ind1]
        power = power[ind0:ind1]
       	wset, 2	
	plot, pfreq, power, /xs, /ys, ytitle='log!L10!N(Power Rs!U-1!N)', $
                xtitle='log!L10!N(k Rs!U-1!N)';, yr=[1e8, 1e12]	

	wset, 1
	;spectra[where(spectra gt 8e5)]=4e5
	spectra = smooth(spectra, 100, /edge_mirror)
	response = spectra/max(spectra)
	plot, freq, response, /xs, /ys
	rfreq=freq
	save, response, rfreq, filename='nfar_response.sav'


stop
END
