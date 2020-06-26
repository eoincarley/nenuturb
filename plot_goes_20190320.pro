;**********************************************;
;				Plot GOES

pro plot_goes_20190320, t1, t2, position=position, download=download, tlines=tlines
	
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
	;			 Xray
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

	outplot, anytim([tlines[0], tlines[0]], /cc), yrange, linestyle=0, thick=3
        outplot, anytim([tlines[1], tlines[1]], /cc), yrange, linestyle=5, thick=5
	outplot, anytim([tlines[2], tlines[2]], /cc), yrange, linestyle=0, thick=3		

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
