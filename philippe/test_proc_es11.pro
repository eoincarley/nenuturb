;----------------------------------------------------
 pro PLOTSPDYN, data,t0,t1,f0,f1,title,col=col
;----------------------------------------------------
!p.multi=[0,2,2]
if keyword_set(col) then begin
  loadct,23
  SPDYNPS, data,t0,t1,f0,f1,'Time (min)','Frequency (MHz)',title,0,0,0,0.01,0.95,1,'lin'
  SPDYNPS, data,t0,t1,f0,f1,'Time (min)','Frequency (MHz)',title+' -bck',0,1,-1.,0.,0.95,1,'lin'
  SPDYNPS, 10.*alog10(data),t0,t1,f0,f1,'Time (min)','Frequency (MHz)',title,0,0,0,0.02,0.98,1,'dB'
  SPDYNPS, 10.*alog10(data),t0,t1,f0,f1,'Time (min)','Frequency (MHz)',title+' -bck',0,1,-1.,0.,0.95,1,'dB'
endif else begin
  loadct,0
  SPDYNPS, data,t0,t1,f0,f1,'Time (min)','Frequency (MHz)',title,0,0,0,0.,0.9,0,'lin'
  SPDYNPS, data,t0,t1,f0,f1,'Time (min)','Frequency (MHz)',title+' -bck',0,1,-1.,0.,0.9,0,'lin'
  SPDYNPS, 10.*alog10(data),t0,t1,f0,f1,'Time (min)','Frequency (MHz)',title,0,0,0,0.03,0.92,0,'dB'
  SPDYNPS, 10.*alog10(data),t0,t1,f0,f1,'Time (min)','Frequency (MHz)',title+' -bck',0,1,-1.,0.,0.9,0,'dB'
endelse
return
end

;----------------------------------------
path = '/databf2/nenufar-tf/ES11/2019/03/20190320_104900_20190320_125000_SUN_TRACKING_BHR/'
file = 'SUN_TRACKING_20190320_104936_0.spectra'        
filename = path+file

;set_ps,'test_proc_es11.ps',/landscape
;device,/color

window, 0, xs=700, ys=700
!p.charsize=1.3
!p.multi=[0,2,2]
loadct,13

; computation of a robust background from the 5% intensity quantile of a quiet interval
READ_NU_SPEC, filename, data,time,freq,beam,ndata, nt,dt,nf,df,ns, jd0,h0, tmax=4*60.
m5=fltarr(nf)
for i=0,nf-1 do m5(i)=dyn_n(data(*,i),0.05)
m=m5/(1-gauss_cvf(0.05)/sqrt(df*1000.*dt/2.))		; /2=apod, m <=> SEFD
plot_io,freq,m,/xsty,xtit='Frequency (MHz)',ytit='Background',tit='0-14 min > start'

; ad-hoc correction curve from a method similar to fflat=3 but applied to the robust background
nchannels=h0.fftlen & nbeamlets=h0. nbeamlets
cor=reform(m,nchannels,nbeamlets)
cor=cor/rebin(reform(median(cor,dimension=1),1,nbeamlets),nchannels,nbeamlets)
cor=reform(cor,nf)
plot,freq,cor,/xsty,xtit='Frequency (MHz)',ytit='Correction',/ynoz
plot_io,freq,m/cor,/xsty,xtit='Frequency (MHz)',ytit='Background corr.'
stop

; Determination of the effective area at time of observation
; Coordinates => Sun elev ~42° on 2019/03/20, azim ~180° => use of nenufar_gain.pro (on /cep/lofar/nenufar/pro/general)
a=[-1] & for f=20,60,5 do begin & nenufar_gain,42,180,f,ae & a=[a,ae] & endfor & a=a[1:*]
f=findgen(9)*5.+20
plot,f,a,psym=4,xra=[19,61],/xsty,/ynoz,xtit='Frequency (MHz)',ytit='A!Deff!N (m!U2!N)'
x=poly_fit(f,a,3)
ae=x(0)+x(1)*freq+x(2)*(freq^2)+x(3)*(freq^3)
oplot,freq,ae

; calibration curve corrf0 for frequency ramp freq0
T=60.*(300./freq)^2.55
k=1.38e3
SEFD=2.*k*T/ae
corrf0 = m/SEFD & freq0=freq

t0 = 33.0 & t1 = 34.5 & f0 = 21.0 & f1 = 24.0

; raw data
;READ_NU_SPEC, filename, data, time, freq, beam, ndata, nt, dt, nf, df, ns, tmin=t0*60.0, tmax=t1*60.0, ntimes=2, fmin=f0, fmax=f1, /exactfreq
;PLOTSPDYN, data,t0,t1,f0,f1,'raw'

; fflat=3
;READ_NU_SPEC, filename, data, time, freq, beam, ndata, nt, dt, nf, df, ns, tmin=t0*60.0, tmax=t1*60.0, ntimes=2, fmin=f0, fmax=f1, /exactfreq, fflat=3
;PLOTSPDYN, data,t0,t1,f0,f1,'fflat=3'

; ad-hoc correction
;READ_NU_SPEC, filename, data, time, freq, beam, ndata, nt, dt, nf, df, ns, tmin=t0*60.0, tmax=t1*60.0, ntimes=2
;data=data/rebin(reform(cor,1,nf),nt,nf)
;w=where(freq ge f0 and freq le f1)
;data=data(*,w)
;PLOTSPDYN, data,t0,t1,f0,f1,'cor'

; fflat=4
freq=freq0 & corrf=corrf0
READ_NU_SPEC, filename, data, time, freq, beam, ndata, nt, dt, nf, df, ns, jd0,h0, corrf, tmin=t0*60.0, tmax=t1*60.0, ntimes=2, fmin=f0, fmax=f1, /exactfreq, fflat=4
PLOTSPDYN, data,t0,t1,f0,f1,'fflat=4 (Jy)'
stop
;--- same but in addition excluding channel 0 ---

; computation of a robust background from the 5% intensity quantile of a quiet interval
READ_NU_SPEC, filename, data,time,freq,beam,ndata, nt,dt,nf,df,ns, jd0,h0, tmax=14*60., ex_chan=[0], /fill
m5=fltarr(nf)
for i=0,nf-1 do m5(i)=dyn_n(data(*,i),0.05)
m=m5/(1-gauss_cvf(0.05)/sqrt(df*1000.*dt/2.))		; /2=apod, m <=> SEFD
plot_io,freq,m,/xsty,xtit='Frequency (MHz)',ytit='Background',tit='0-14 min > start, ex_chan=[0]'

; ad-hoc correction curve from a method similar to fflat=3 but applied to the robust background
nchannels=h0.fftlen & nbeamlets=h0. nbeamlets
cor=reform(m,nchannels,nbeamlets)
cor=cor/rebin(reform(median(cor,dimension=1),1,nbeamlets),nchannels,nbeamlets)
cor=reform(cor,nf)
plot,freq,cor,/xsty,xtit='Frequency (MHz)',ytit='Correction',/ynoz
plot_io,freq,m/cor,/xsty,xtit='Frequency (MHz)',ytit='Background corr.'

; Determination of the effective area at time of observation
; Coordinates => Sun elev ~42° on 2019/03/20, azim ~180° => use of nenufar_gain.pro (on /cep/lofar/nenufar/pro/general)
a=[-1] & for f=20,60,5 do begin & nenufar_gain,42,180,f,ae & a=[a,ae] & endfor & a=a[1:*]
f=findgen(9)*5.+20
plot,f,a,psym=4,xra=[18,62],/xsty,/ynoz,xtit='Frequency (MHz)',ytit='A!Deff!N (m!U2!N)'
x=poly_fit(f,a,3)
ae=x(0)+x(1)*freq+x(2)*(freq^2)+x(3)*(freq^3)
oplot,freq,ae

; calibration curve corrf0 for frequency ramp freq0
T=60.*(300./freq)^2.55
k=1.38e3
SEFD=2.*k*T/ae
corrf0 = m/SEFD & freq0=freq

t0 = 33.0 & t1 = 34.5 & f0 = 21.0 & f1 = 24.0

; raw data
;READ_NU_SPEC, filename, data, time, freq, beam, ndata, nt, dt, nf, df, ns, tmin=t0*60.0, tmax=t1*60.0, ntimes=2, fmin=f0, fmax=f1, /exactfreq, ex_chan=[0], /fill
;PLOTSPDYN, data,t0,t1,f0,f1,'raw, ex_chan=[0]'

; fflat=3
;READ_NU_SPEC, filename, data, time, freq, beam, ndata, nt, dt, nf, df, ns, tmin=t0*60.0, tmax=t1*60.0, ntimes=2, fmin=f0, fmax=f1, /exactfreq, fflat=3, ex_chan=[0], /fill
;PLOTSPDYN, data,t0,t1,f0,f1,'fflat=3, ex_chan=[0]'

; ad-hoc correction
;READ_NU_SPEC, filename, data, time, freq, beam, ndata, nt, dt, nf, df, ns, tmin=t0*60.0, tmax=t1*60.0, ntimes=2, ex_chan=[0], /fill
;data=data/rebin(reform(cor,1,nf),nt,nf)
;w=where(freq ge f0 and freq le f1)
;data=data(*,w)
;PLOTSPDYN, data,t0,t1,f0,f1,'cor, ex_chan=[0]'

; fflat=4
freq=freq0 & corrf=corrf0
READ_NU_SPEC, filename, data, time, freq, beam, ndata, nt, dt, nf, df, ns, jd0,h0, corrf, tmin=t0*60.0, tmax=t1*60.0, ntimes=2, fmin=f0, fmax=f1, /exactfreq, fflat=4, ex_chan=[0], /fill
PLOTSPDYN, data,t0,t1,f0,f1,'fflat=4 (Jy), ex_chan=[0]'

;device,/close
;ps_pdf,'test_proc_es11.ps',/rem
;exit_ps
end
