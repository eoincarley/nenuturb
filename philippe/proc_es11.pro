;---------------------------------------------------------
  pro PROC_ES11, pathobs, files
;---------------------------------------------------------

if strmid(pathobs,strlen(pathobs)-1,1) ne '/' then pathobs=pathobs+'/'
if n_elements(files) eq 0 then begin
  spawn,'ls -1 '+pathobs+'*.spectra > list'
  files=['rien'] & buf=''
  on_ioerror,suite
  openr,v,'list',/get_lun
encore:
  readf,v,buf
  files=[files,DELPATH(buf)]
  goto,encore
suite:
  close,v & free_lun,v
  files=files[1:*]
endif

onoff=['B0','B1','B2','B3']

for ifile=0,n_elements(files)-1 do begin

  name=strmid(files(ifile),0,strlen(files(ifile))-8)
  set_ps,name+'.ps',/portrait
  device,/color
  !p.charsize=1.3
  loadct,0
  ionoff=fix(strmid(name,strlen(name)-1,1))

; first 5 min, raw, high spectral resolution

  READ_NU_SPEC, pathobs+files(ifile), data,time,freq,beam,ndata, nt,dt,nf,df,ns, jd0,h0, corrf,corrt,fref, tmax=300,ntimes=100,nstokes=2,writetxt=1
  if n_elements(where(data(*,*,0) eq 0)) ge n_elements(data(*,*,0))/2 then begin
    plot,[0,1],[0,1],/nodata & xyouts,0.5,0.5,'Corrupted file',charsize=3.,align=0.5
    goto,makepdf
  endif
  df0=df
  ntimes=round(0.250*100/dt)
  print,'ntimes=',ntimes
  NU_PLOTS_ES11, data, freq, time, name+', '+onoff(ionoff)+', first 5 min, raw', FZ=[mean(freq)-1,mean(freq)+1]

; raw

  READ_NU_SPEC, pathobs+files(ifile), data,time,freq,beam,ndata, nt,dt,nf,df,ns, jd0,h0, corrf,corrt,fref,nchannels=8,ntimes=ntimes,nstokes=3,block_inc=1000
  nplots=ceil(max(time)/3600.)		; number of plots (1 / hour)
  NU_PLOTS_ES11, data, freq, time, name+', '+onoff(ionoff)+', entire file, raw', /L, /FZ, TZ=[1,600]

  for iplot=0,nplots-1 do begin
    w=where(time ge iplot*3600. and time le (iplot+1)*3600., nw)
    if nw gt 10 then begin
      nw=long(nw/2)*2 & w=w(0:nw-1)
      !p.multi=[0,1,3]
      SPDYNPS, rebin(10.*alog10(data(w,*,0)),nw/2,nf),min(time(w))/60,max(time(w))/60,min(freq),max(freq),'Time since file start (min)','Frequency (MHz)',name+', '+onoff(ionoff)+', raw, Stokes I',0,1,-0.5,0.,0.98,0,'dB'
      SPDYNPS, rebin(data(w,*,1)/data(w,*,0),nw/2,nf),min(time(w))/60,max(time(w))/60,min(freq),max(freq),'Time since file start (min)','Frequency (MHz)',name+', '+onoff(ionoff)+', raw, Stokes V/I',0,0,0,0.02,0.98,0,'V/I'
      SPDYNPS, rebin(data(w,*,2)/data(w,*,0),nw/2,nf),min(time(w))/60,max(time(w))/60,min(freq),max(freq),'Time since file start (min)','Frequency (MHz)',name+', '+onoff(ionoff)+', raw, Stokes L/I',0,0,0,0.02,0.98,0,'L/I'
    endif
  endfor

goto,makepdf

  x=reform(rebin(data(*,*,0),1,nf/4,1))
  xf=rebin(freq,nf/4)
  restore,'q5ref.sav'
  FIND_FILTER, pathobs+files(ifile), jd0, filterN, filterT
  filt5=min(filterN)
  x5=INTERPOL(q5(*,filt5),fq5,xf)
  LE_AUTO_,x/(10.^(x5/10.)),31,6.,0,xnet,xpar,EDGE=2	; ********** ADJUST **********
  xpar=1b-dilate(erode(1-xpar,[1,1]),[1,1])		; omit 1-channels gaps
  w=where(xpar eq 0)
  if w(0) ne -1 then ex_beamlets=w else ex_beamlets=[]
  print,'ex_beamlets = ',ex_beamlets

; proc

  FIND_BITMODE,pathobs+files(ifile),nbits
  if nbits eq 8 then ex_chan=[0,round(195.3125/df0)/2] else ex_chan=[0]
  print,'ex_chan = ',ex_chan

  READ_NU_SPEC, pathobs+files(ifile), data,time,freq,beam,ndata, nt,dt,nf,df,ns, jd0,h0, corrf,corrt,fref, $
	nchannels=4,ntimes=ntimes,nstokes=3,ex_chan=ex_chan,ex_beamlets=ex_beamlets,fclean=[4.,101],/bclean,/tclean,tflat=[4,36.5,70.],/fill,writefits=2,writetxt=2,block_inc=1000

  w=where(rebin(ndata,1,nf) eq 0, nw)
  if w(0) ne -1 then begin
    m=reform(median(data,dimension=2),nt,1,ns)
    for i=0,nw-1 do data(*,w(i),*)=m
  endif

  nplots=ceil(max(time)/3600.)		; number of plots (1 / hour)
  titre=name+', '+onoff(ionoff)+', proc'
  NU_PLOTS_ES11, data, freq, time, titre, ndata, /L, /FZ, TZ=[1,600]

  loadct,13
  !p.multi=[0,1,4]
  s=size(corrt)
  plot,time/3600.,corrt(*,0),xtit='Time since file start (h)',ytit='Gain(t) correction : corrt(0)',tit=titre,/xsty,/ynoz
  if s(0) gt 1 then begin
    plot,time/3600.,corrt(*,1),xtit='Time since file start (h)',ytit='data freq.-integ. time profile & fit',tit=titre,/xsty,/ynoz
    oplot,time/3600.,corrt(*,2),color=250
    plot,time/3600.,corrt(*,3),xtit='Time since file start (h)',ytit='[a].log(t) + b',/xsty,/ynoz
    plot,time/3600.,corrt(*,4),xtit='Time since file start (h)',ytit='a.log(t) + [b]',/xsty,/ynoz
  endif

  w=where(time ge 1 and time le 600)
  plot,time(w)/60.,corrt(w,0),xtit='Time since file start (min)',ytit='Gain(t) correction : corrt(0)',tit=titre,/xsty,/ynoz
  if s(0) gt 1 then begin
    plot,time(w)/60.,corrt(w,1),xtit='Time since file start (min)',ytit='data freq.-integ. time profile & fit',tit=titre,/xsty,/ynoz
    oplot,time(w)/60.,corrt(w,2),color=250
  endif
  loadct,0

  !p.multi=[0,1,4]
  for iplot=0,nplots-1 do begin
    w=where(time ge iplot*3600. and time le (iplot+1)*3600., nw)
    if nw gt 10 then begin
      nw=long(nw/2)*2 & w=w(0:nw-1)
      SPDYNPS, rebin(10.*alog10(data(w,*,0)),nw/2,nf),min(time(w))/60,max(time(w))/60,min(freq),max(freq),'Time since file start (min)','Frequency (MHz)',titre+', Stokes I',0,1,-0.5,0.,0.98,0,'dB'
      SPDYNPS, rebin(data(w,*,1)/data(w,*,0),nw/2,nf),min(time(w))/60,max(time(w))/60,min(freq),max(freq),'Time since file start (min)','Frequency (MHz)',titre+', Stokes V/I',0,0,0,0.02,0.98,0,'V/I'
      SPDYNPS, rebin(data(w,*,2)/data(w,*,0),nw/2,nf),min(time(w))/60,max(time(w))/60,min(freq),max(freq),'Time since file start (min)','Frequency (MHz)',titre+', Stokes L/I',0,0,0,0.02,0.98,0,'L/I'
      SPDYNPS, rebin(ndata(w,*)*100,nw/2,nf),min(time(w))/60,max(time(w))/60,min(freq),max(freq),'Time since file start (min)','Frequency (MHz)',titre+', pixOK=1-flag',0,0,0,0.,1.,1,'%'
    endif
  endfor

makepdf:

  device,/close
  ps_pdf,name+'.ps',/rem

endfor

exit_ps
end

