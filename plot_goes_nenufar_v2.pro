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

;**********************************************;
;                               Plot GOES

pro plot_goes_20190320, t1, t2, position=position, download=download, tlines=tlines, imgtims=imgtims

        if keyword_set(download) then begin
                use_network
                goes_obj = ogoes()
                goes_obj->set, tstart=anytim(t1, /atim), tend=anytim(t2, /atim)
                goes = goes_obj->getdata()
                times = goes_obj->getdata(/times)
                utbase = goes_obj->get(/utbase)
                save, goes, times, utbase, filename='goes_20190320.sav'
        endif else begin
                restore, 'goes_20190320.sav', /verb
        endelse

        yrange = [1e-9, 1e-5]
        ;--------------------------------;
        ;                        Xray
        set_line_color
        utplot, times, goes[*,0], utbase, $
                        thick = 3, $
                        ytit = 'Watts m!U-2!N', $
                        xtit = ' ', $
                        color = 3, $
                        ;xrange = [t1, t2], $
                        /xs, $
                        /ys, $
                        yrange = yrange, $
                        /ylog, $
                        /normal, $
                        /noerase, $
                        position=position

        outplot, times, goes[*, 1], utbase, color=5, thick=3, linestyle=0

        axis, yaxis=1, yrange=yrange, /ys, ytickname=[' ','A','B','C','M']
        axis, yaxis=0, yrange=yrange, /ys

        plots, times, 1e-8, linestyle=2
        plots, times, 1e-7, linestyle=2
        plots, times, 1e-6, linestyle=2

        outplot, anytim([tlines[0], tlines[0]], /cc), yrange, linestyle=0, thick=2
        outplot, anytim([tlines[1], tlines[1]], /cc), yrange, linestyle=0, thick=2
        ;outplot, anytim([tlines[2], tlines[2]], /cc), yrange, linestyle=0, thick=3

	
	; Plot vertical lines at times of images
	outplot, anytim([imgtims[0], imgtims[0]], /cc), yrange, linestyle=1, thick=0.5
        outplot, anytim([imgtims[1], imgtims[1]], /cc), yrange, linestyle=1, thick=0.5
        outplot, anytim([imgtims[2], imgtims[2]], /cc), yrange, linestyle=1, thick=0.5



        legend, ['GOES15 0.1-0.8nm','GOES15 0.05-0.4nm'], $
                        linestyle=[0,0], $
                        color=[3,5], $
                        box=0, $
                        ;pos = [0.12, 0.98], $
                        /normal, $
                        charsize=0.8, $
                        thick=3, $
                        /right


END

pro plot_nfar, data, utimes, freq, scl=scl, pos=pos, imgtims=imgtims

        data = reverse(data, 2)
        freq = reverse(freq)
	data = smooth(data, 3)
        ;data = alog10(data)

	loadct, 74
        reverse_ct
        ;spectro_plot, bytscl(data, scl[0], scl[1]), utimes, freq, /xs, /ys, $
        spectro_plot, sigrange(alog10(data)), utimes, freq, /xs, /ys, $
		ytitle=' ', xtitle=' ', yr=[60, 20], /noerase, XTICKFORMAT="(A1)", YTICKFORMAT="(A1)", $
                position=pos, /ylog , /normal, title=' ', ytickname=['60', '50', '40', '30', '20'], xr=[utimes[0], utimes[-1]], yticklen=-5e-10

        ;------------------------------------------;
        ;     Empty template to get black ticks
        ;
        loadct, 0
        utplot, utimes, freq, yr=[60,20], /xs, /ys, /ylog, /nodata, /noerase, $
                xtitle='Time (UT)', yticklen=-5e-3, $
                title='  ', pos=pos, /normal, color=250, ytickname=['60', '50', '40', '30', '20'], ytickv=[60, 50, 40, 30, 20], xr=[utimes[0], utimes[-1]]

	yrange = [60.0, 20.0]
	; Plot vertical lines at times of images
        outplot, anytim([imgtims[0], imgtims[0]], /cc), yrange, linestyle=1, thick=1
        outplot, anytim([imgtims[1], imgtims[1]], /cc), yrange, linestyle=1, thick=1
        outplot, anytim([imgtims[2], imgtims[2]], /cc), yrange, linestyle=1, thick=1


	data = cgScaleVector(Findgen(101), 100, 1000)
   	ticks = LogLevels([20,60], /fine)
   	nticks = N_Elements(ticks)
	set_line_color
   	cgPlot, data, /YLOG, YRANGE=[60,20], YTICKS=nticks-1, $
      		YTICKV=Reverse(ticks), ticklen=-5e-3,  pos=pos, color=0, AXESCOLOR=0, /noerase, XTICKFORMAT="(A1)", ytitle='Frequency (MHz)' 

end



pro plot_goes_nenufar_v2, save=save, postscript=postscript

	!p.charsize=1.5
	path = '/databf2/nenufar-tf/ES11/2019/03/20190320_104900_20190320_125000_SUN_TRACKING_BHR/'
        file = 'SUN_TRACKING_20190320_104936_0.spectra'
	imgtims = anytim('2019-03-20T'+['11:17:57', '11:23:09', '11:29:57'], /cc)

	if keyword_set(postscript) then begin
                setup_ps, './eps/nfar_goes_20190320.eps'
        endif else begin
		loadct, 0
		!p.color=0
		!p.background=200
                window, xs=1200, ys=700
        endelse 

        ;------------------------------------;
	; 	First spectrogram
	;
	t0 = 28.0 & t1 = 44.0 
	restore, 'calibration_factor.sav'
        freq=freq0 & corrf=corrf0

        ;READ_NU_SPEC, path+file, data, time, freq, beam, ndata, nt, dt, nf, df, ns, $
        ;       jd0, h0, corrf, tmin=t0*60.0, tmax=t1*60.0, ntimes=8, nchannels=512, /exactfreq, fflat=4, /fill

	READ_NU_SPEC, path+file, data, time, freq, beam, ndata, nt, dt, nf, df, ns, $
		ntimes=10, tmin=28.0*60.0, tmax=44.0*60.0, fflat=3
        utimes0=anytim(file2time(file), /utim) + time
	
	plot_nfar, data, utimes0,  freq, scl=[4.5e5, 5e7], pos=[0.2, 0.45, 0.85, 0.7], imgtims=imgtims
	
		
	;-------------------------;
        ;       Plot GOES
        ;
        pos = [0.2, 0.75, 0.85, 0.95]
        plot_goes_20190320, '2019-03-20T10:00', '2019-03-20T13:00', position=pos, tlines=[utimes0[0], utimes0[-1]], imgtims=imgtims


	if keyword_set(postscript) then begin
		device, /close
		set_plot, 'x'
	endif	
	  
	if keyword_set(save) then $
                save, data, utimes, freq, filename='nfar_'+time2file(utimes[0])+'_sun_tracking.sav'

stop
END
