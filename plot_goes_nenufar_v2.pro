pro setup_ps, name
  
   set_plot,'ps'
   !p.font=0
   !p.charsize=1.0
   device, filename = name, $
          /color, $
          /helvetica, $
          /inches, $
          xsize=12, $
          ysize=8, $
          /encapsulate, $
          yoffset=5, $
          bits_per_pixel = 16

end

pro plot_nfar, data, utimes, freq, scl=scl, pos=pos

        data = reverse(data, 2)
        freq = reverse(freq)
	data = smooth(data, 3)
        ;data=alog10(data)

	loadct, 74
        reverse_ct
        spectro_plot, bytscl(data, scl[0], scl[1]), utimes, freq, /xs, /ys, $
                ytitle=' ', xtitle=' ', yr=[57, 20], /noerase, XTICKFORMAT="(A1)", YTICKFORMAT="(A1)", $
                position=pos, /normal, title=' ', xr=[utimes[0], utimes[-1]]

        ;------------------------------------------;
        ;     Empty template to get black ticks
        ;
        loadct, 0
        utplot, utimes, freq, yr=[57,20], /xs, /ys, /nodata, /noerase, $
                xtitle='Time (UT)', ytitle='Frequency (MHz)', $
                title='  ', pos=pos, /normal, color=250, xr=[utimes[0], utimes[-1]]

end



pro plot_goes_nenufar_v2, save=save, postscript=postscript

	!p.charsize=1.5
	path = '/databf2/nenufar-tf/ES11/2019/03/20190320_104900_20190320_125000_SUN_TRACKING_BHR/'
        file = 'SUN_TRACKING_20190320_104936_0.spectra'


	if keyword_set(postscript) then begin
                setup_ps, './nfar_goes_20190320.eps'
        endif else begin
                window, xs=1200, ys=700
        endelse 

        ;------------------------------------;
	; 	First spectrogram
	;
	READ_NU_SPEC, path+file, data,time,freq,beam,ndata,nt,dt,nf,df,ns, ntimes=10, tmin=29*60.0, tmax=45*60.0
        utimes0=anytim(file2time(file), /utim) + time
	
	if keyword_set(rebin) then begin
                nfbin = (size(data))[2]
                data = data[0:11999, *]
                tbin = 3000
                data = rebin(data, tbin, nfbin)
                utimes = congrid(utimes, tbin)
                ntsteps=1
        endif
	plot_nfar, data, utimes0, freq, scl=[4.5e5, 5e7], pos=[0.2, 0.45, 0.85, 0.7]
	
	;------------------------------------;
	;	Second spectrogram
	;
	;READ_NU_SPEC, path+file, data,time,freq,beam,ndata,nt,dt,nf,df,ns, tmin=39*60.0, tmax=51*60.0
        ;utimes1=anytim(file2time(file), /utim) + time
	if keyword_set(rebin) then begin
                nfbin = (size(data))[2]
                data = data[0:11999, *]
                tbin = 3000
                data = rebin(data, tbin, nfbin)
                utimes = congrid(utimes, tbin)
                ntsteps=1
        endif
	;plot_nfar, data, utimes1, freq, scl=[4.5e5, 4e7], pos=[0.2, 0.12, 0.85, 0.37]
		
	;-------------------------;
        ;       Plot GOES
        ;
        pos = [0.2, 0.75, 0.85, 0.95]
        plot_goes_20190320, '2019-03-20T10:00', '2019-03-20T13:00', position=pos, tlines=[utimes0[0], utimes0[-1]]


	if keyword_set(postscript) then begin
		device, /close
		set_plot, 'x'
	endif	
	  
	if keyword_set(save) then $
                save, data, utimes, freq, filename='nfar_'+time2file(utimes[0])+'_sun_tracking.sav'

stop
END
