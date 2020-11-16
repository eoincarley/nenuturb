pro setup_ps, name
  
   set_plot,'ps'
   !p.font=0
   !p.charsize=1.8
   device, filename = name, $
          /color, $
          /helvetica, $
          /inches, $
          xsize=19, $
          ysize=7, $
          /encapsulate, $
          yoffset=5, $
          bits_per_pixel = 16

end

pro plot_nenufar_typeII, save=save, postscript=postscript

	!p.charsize=1.5
	path = '/databf2/nenufar-tf/ES11/2019/03/20190320_104900_20190320_125000_SUN_TRACKING_BHR/'
        file = 'SUN_TRACKING_20190320_104936_0.spectra'
	;read_nu_spec, path+file, /info

	;for i=0, 2 do begin	
	  READ_NU_SPEC, path+file, data,time,freq,beam,ndata,nt,dt,nf,df,ns, $
		tmin=36*60.0, tmax=46*60.0
		  ;tmin=(i*15.0+11.0)*60.0, tmax=((i+1.0)*15.0+11.0)*60.0 
		 

	  utimes=anytim(file2time(file), /utim) + time

	  ;trebin = (size(data))[1]/6
	  ;frebin = (size(data))[2]/1
	  data = reverse(data, 2)
	  freq = reverse(freq)
	  ;-------------------------------------------;
	  ;	      Rebin if necessary
	  ;
	  ;data = constbacksub(data, /auto)
	  ;data = rebin(data, trebin, frebin)
	  ;freq = rebin(freq, frebin)
	  ;utimes = rebin(utimes, trebin)
	   
	  data = smooth(data,3)
	  ;for i=0, nf-1 do $
	  ;	data[*, i] = data[*, i] - smooth(data[*,i],1000)

	  if keyword_set(postscript) then begin
		setup_ps, './nfar_'+time2file(utimes[0])+'_typeII.eps'
	  endif else begin
		window, xs=1200, ys=700
	  endelse	

	  ;------------------------------------------;
	  ;	Empty template to get black ticks
	  ;
	  loadct, 0
          utplot, utimes, freq, yr=[57,20], /xs, /ys, xtitle='Time (UT)', ytitle='Frequency (MHz)', $
		title='NenuFAR-ES11 '+time2file(utimes[0], /date), pos=[0.1, 0.12, 0.9, 0.9], /normal, color=150
	
          ;------------------------------------------;
	  ;		Plot spectrogram
	  ;
	  loadct, 74
          reverse_ct
	  spectro_plot, sigrange(data), utimes, freq, /xs, /ys, $
		ytitle=' ', xtitle=' ', yr=[57, 20], /noerase, XTICKFORMAT="(A1)", YTICKFORMAT="(A1)", $
		position=[0.1, 0.12, 0.9, 0.9], /normal

	
	  if keyword_set(postscript) then begin
		device, /close
		set_plot, 'x'
	  endif	
	  
	  if keyword_set(save) then $
                save, data, utimes, freq, filename='nfar_'+time2file(utimes[0])+'_sun_tracking.sav'
;	endfor

stop
END
