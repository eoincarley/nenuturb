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
   	
	READ_NU_SPEC, file, data, time, freq, beam, ndata, nt, dt, nf, df, ns, $
      tmin=t0*60.0, tmax=t1*60.0, fmin=f0, fmax=f1, fflat=3, ntimes=8, ex_chan=[0], /fill, /exactfreq
       

	print, 'Time resolution: '+string(dt)
	utimes=anytim(file2time(file), /utim) + time
        data = reverse(data, 2)
        freq = reverse(freq)
end

pro plot_spectro, data, time, freq, f0, f1, posit

	loadct, 0
  utplot, time, freq, yr=[f1,f0], /xs, /ys, xtitle='Time (UT)', ytitle='Frequency (MHz)', $
                title=' ', pos=posit, /normal, color=150, $
                xr=[time[0], time[-1]], /noerase
	;------------------------------------------;
  ;            Plot spectrogram
  ;
  loadct, 74
  reverse_ct
  spectro_plot, sigrange(data), time, freq, /xs, /ys, $
          ytitle=' ', xtitle=' ', yr=[f1, f0], /noerase, XTICKFORMAT="(A1)", YTICKFORMAT="(A1)", $
          position=posit, /normal, xr=[time[0], time[-1]]


END

pro typeIIb_stria, postscript=postscript


  ; Taking a look at the stria drift in the type IIb

  path = '/volumes/plasma/nenufar/nenufar-tf/ES11/2019/03/20190320_104900_20190320_125000_SUN_TRACKING_BHR/'
  file = 'SUN_TRACKING_20190320_104936_0.spectra'

	t0 = 33.61
  t1 = 33.66
  f0 = 21.83
  f1 = 21.96
	read_nfar_data, path+file, t0, t1, f0, f1, data=data, utimes=utimes, freq=freq
	   

	if keyword_set(postscript) then begin
		setup_ps, './eps/psd_hbone_drift_F.eps', xsize=18, ysize=5.5
	endif else begin
		!p.charsize=1.8
		window, 0, xs=700, ys=1400
	endelse	

	posit=[0.15, 0.55, 0.9, 0.9]

  ;------------------------------------------;
	;	     Plot spectrogram
	;
  plot_spectro, data, utimes, freq, f0, f1, posit

  ;------------------------------------------;
  ;
  ;   Interpolate across the map
  ;
  map = make_map(data)

  inter = 20.

  ti = utimes - utimes[0]
  nt = n_elements(ti)
  dti = ti[1]-ti[0]
  ti2 = findgen(nt*inter)*dti/inter+ti[0]

  fi = freq ;t3_freq_r[t3_freq_range]
  nf = n_elements(fi)
  dfi = fi[1]-fi[0]
  fi2 = findgen(nf*inter)*dfi/inter+fi[0]

  map_arr2 = interp2d(data,ti,fi,ti2,fi2,/grid)

  timnew = utimes[0] + ti2
  freqnew = fi2
  data = reverse(map_arr2, 2)
  posit=[0.15, 0.15, 0.9, 0.45]
  plot_spectro, smooth(data, 15, /edge_mirror), timnew, freqnew, f0, f1, posit

  if keyword_set(postscript) then device, /close

stop	
END
