;-----------------------------------------------------------------------
 pro NU_PLOTS_ES11, data, freq, time, titre, ndata, L=L, FZ=FZ, TZ=TZ
;-----------------------------------------------------------------------
; /L to plot linear polarization (default = I & V)
; /FZ to plot overall spectrum, FZ=[fmin,fmax] to plot overall spectrum + a zoom in frequency
; /TZ to plot overall time series, TZ=[tmin,tmax] to plot overall time series + a zoom in time
; NB: do not use fmin=0 or tmin=0

  !p.multi=[0,1,4]

  if keyword_set(FZ) then begin
    nf=n_elements(freq)
    s=rebin(10.*alog10(data(*,*,0)>1),1,nf)
    plot,freq,s,xtit='Frequency (MHz)',ytit='Stokes I (dB)',tit=titre,/xsty,/ynoz
    v=rebin(data(*,*,1)/data(*,*,0),1,nf)
    plot,freq,v,xtit='Frequency (MHz)',ytit='Stokes V/I',yra=[min(v)>(-1.),max(v)<1.],/xsty,/ysty
    oplot,[0,100],[0,0],line=1
    if keyword_set(L) then begin
      l=rebin(data(*,*,2)/data(*,*,0),1,nf)
      plot,freq,l,xtit='Frequency (MHz)',ytit='Stokes L/I',/xsty,/ynoz
    endif
    if n_elements(ndata) gt 0 then begin
      n=rebin(ndata(*,*)*100,1,nf)
      plot,freq,n,xtit='Frequency (MHz)',ytit='Pix.Ok = 1-Flag (%)',yra=[-5,105],/xsty,/ysty
      oplot,[0,100],[0,0],line=1
      oplot,[0,100],[100,100],line=1
    endif
  endif
  if (keyword_set(L) and n_elements(ndata) eq 0) or (not(keyword_set(L)) and n_elements(ndata) gt 0) then !p.multi=[0,1,4]

  if n_elements(FZ) gt 1 then begin
    fmin=FZ(0) & fmax=FZ(1)
    w=where(freq ge fmin and freq le fmax)
    plot,freq(w),s(w),xtit='Frequency (MHz)',ytit='Stokes I (dB)',tit=titre,/ynoz,/xsty
    plot,freq(w),v(w),xtit='Frequency (MHz)',ytit='Stokes V/I',yra=[min(v(w))>(-1.),max(v(w))<1.],/xsty,/ysty
    oplot,[0,100],[0,0],line=1
    if keyword_set(L) then plot,freq(w),l(w),xtit='Frequency (MHz)',ytit='Stokes L/I',/xsty,/ynoz
    if n_elements(ndata) gt 0 then begin
      plot,freq(w),n(w),xtit='Frequency (MHz)',ytit='Pix.Ok = 1-Flag (%)',yra=[-5,105],/xsty,/ysty
      oplot,[0,100],[0,0],line=1
      oplot,[0,100],[100,100],line=1
    endif
  endif
  if (keyword_set(L) and n_elements(ndata) eq 0) or (not(keyword_set(L)) and n_elements(ndata) gt 0) then !p.multi=[0,1,4]

  if keyword_set(TZ) then begin
    nt=n_elements(time)
    s=rebin(10.*alog10(data(*,*,0)>1),nt,1)
    plot,time/3600.,s,xtit='Time since file start (h)',ytit='Stokes I (dB)',tit=titre,/xsty,/ynoz
    v=rebin(data(*,*,1)/data(*,*,0),nt,1)
    plot,time/3600.,v,xtit='Time since file start (h)',ytit='Stokes V/I',yra=[min(v)>(-1.),max(v)<1.],/xsty,/ysty
    oplot,[0,100],[0,0],line=1
    if keyword_set(L) then begin
      l=rebin(data(*,*,2)/data(*,*,0),nt,1)
      plot,time/3600.,l,xtit='Time since file start (h)',ytit='Stokes L/I',/xsty,/ynoz
    endif
    if n_elements(ndata) gt 0 then begin
      n=rebin(ndata(*,*)*100,nt,1)
      plot,time/3600.,n,xtit='Time since file start (h)',ytit='Pix.Ok = 1-Flag (%)',yra=[-5,105],/xsty,/ysty
      oplot,[0,100],[0,0],line=1
      oplot,[0,100],[100,100],line=1
    endif
  endif
  if (keyword_set(L) and n_elements(ndata) eq 0) or (not(keyword_set(L)) and n_elements(ndata) gt 0) then !p.multi=[0,1,4]

  if n_elements(TZ) gt 1 then begin
    tmin=TZ(0) & tmax=TZ(1)
    w=where(time ge tmin and time le tmax)
    plot,time(w)/60.,s(w),xtit='Time since file start (min)',ytit='Stokes I (dB)',tit=titre,/xsty,/ynoz
    plot,time(w)/60.,v(w),xtit='Time since file start (min)',ytit='Stokes V/I',yra=[min(v(w))>(-1.),max(v(w))<1.],/xsty,/ysty
    oplot,[0,100],[0,0],line=1
    if keyword_set(L) then plot,time(w)/60.,l(w),xtit='Time since file start (min)',ytit='Stokes L/I',/xsty,/ynoz
    if n_elements(ndata) gt 0 then begin
      plot,time(w)/60.,n(w),xtit='Time since file start (min)',ytit='Pix.Ok = 1-Flag (%)',yra=[-5,105],/xsty,/ysty
      oplot,[0,100],[0,0],line=1
      oplot,[0,100],[100,100],line=1
    endif
  endif

return
end

