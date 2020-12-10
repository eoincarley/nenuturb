set_ps,'test_spec_es11.ps',/landscape
device,/color
!p.charsize=1.6
loadct,13
!p.multi=[0,2,2]

filename = '/databf2/nenufar-tf/ES11/2019/03/20190320_104900_20190320_125000_SUN_TRACKING_BHR/SUN_TRACKING_20190320_104936_0.spectra'
t0 = 33.0 & t1 = 34.5 & f0 = 21.0 & f1 = 24.0

READ_NU_SPEC, filename, data, time, freq, beam, ndata, nt, dt, nf, df, ns, tmin=t0*60.0, tmax=t1*60.0, ntimes=2, fmin=f0, fmax=f1, /exactfreq
x=10.*alog10(data)
plot,freq,rebin(x,1,nf),/ynoz,xtit='Frequency (MHz)',ytit='dB',tit='raw'
for i=0,511 do oplot,([i,i]+0.5)*0.1953125,[0,100],line=1
plot,freq,rebin(x(0:1./dt,*),1,nf),yra=[57,79],xtit='Frequency (MHz)',ytit='dB',tit='raw',/ysty
for i=0,511 do oplot,([i,i]+0.5)*0.1953125,[0,100],line=1
for i=0,nt-15/dt,15./dt do oplot,freq,rebin(x(i:i+1./dt,*),1,nf)
x0=x

READ_NU_SPEC, filename, data, time, freq, beam, ndata, nt, dt, nf, df, ns, tmin=t0*60.0, tmax=t1*60.0, ntimes=2, fmin=f0, fmax=f1, /exactfreq, fflat=3, ex_chan=[0], /fill
x=10.*alog10(data)
plot,freq,rebin(x,1,nf),/ynoz,xtit='Frequency (MHz)',ytit='dB',tit='fflat=3, ex_chan=[0], /fill'
for i=0,511 do oplot,([i,i]+0.5)*0.1953125,[0,100],line=1
plot,freq,rebin(x(0:1./dt,*),1,nf),yra=[57,79],xtit='Frequency (MHz)',ytit='dB',tit='fflat=3, ex_chan=[0], /fill',/ysty
for i=0,511 do oplot,([i,i]+0.5)*0.1953125,[0,100],line=1
for i=0,nt-15/dt,15./dt do oplot,freq,rebin(x(i:i+1./dt,*),1,nf)
x1=x

plot,freq,rebin(x0,1,nf),/ynoz,xtit='Frequency (MHz)',ytit='dB',tit='comparison'
oplot,freq,rebin(x1,1,nf),color=250
for i=0,511 do oplot,([i,i]+0.5)*0.1953125,[0,100],line=1
plot,freq,rebin(x0(0:1./dt,*),1,nf),yra=[57,79],xtit='Frequency (MHz)',ytit='dB',tit='comparison',/ysty
for i=0,511 do oplot,([i,i]+0.5)*0.1953125,[0,100],line=1
for i=0,nt-15/dt,15./dt do oplot,freq,rebin(x(i:i+1./dt,*),1,nf)
oplot,freq,rebin(x1(0:1./dt,*),1,nf),color=250
for i=0,nt-15/dt,15./dt do oplot,freq,rebin(x1(i:i+1./dt,*),1,nf),color=250

device,/close
exit_ps
ps_pdf,/rem
end

