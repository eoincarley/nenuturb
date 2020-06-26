;**********************************************;
;				Plot GOES

pro plot_goes_obj, t1, t2

	use_network
	goes_obj = ogoes()
	goes_obj->set, tstart=anytim(t1, /atim), tend=anytim(t2, /atim)
	goes = goes_obj->getdata() 
	times = goes_obj->getdata(/times) 
	utbase = goes_obj->get(/utbase)

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
			yrange = [1e-9,1e-2], $
			/ylog, $
			/normal, $
			/noerase
		
	outplot, times, goes[*, 1], utbase, color=5, thick=3
	
	axis, yaxis=1, yrange=[1e-9, 1e-2], /ys, ytickname=[' ','A','B','C','M','X',' ',' ']
	axis, yaxis=0, yrange=[1e-9, 1e-2], /ys

	plots, times, 1e-8, linestyle=2
	plots, times, 1e-7, linestyle=2
	plots, times, 1e-6, linestyle=2
	plots, times, 1e-5, linestyle=2
	plots, times, 1e-4, linestyle=2
	plots, times, 1e-3, linestyle=2
			
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
