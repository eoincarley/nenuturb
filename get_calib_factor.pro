pro get_calib_factor

path = '/databf2/nenufar-tf/ES11/2019/03/20190320_104900_20190320_125000_SUN_TRACKING_BHR/'
file = 'SUN_TRACKING_20190320_104936_0.spectra'
filename = path+file


READ_NU_SPEC, filename, data,time,freq,beam,ndata, nt,dt,nf,df,ns, jd0,h0, tmax=14*60., ex_chan=[0], /fill
m5 = fltarr(nf)
for i=0, nf-1 do m5[i] = dyn_n(data[*,i], 0.05)
m=m5/(1-gauss_cvf(0.05)/sqrt(df*1000.*dt/2.))           ; /2=apod, m <=> SEFD

; Determination of the effective area at time of observation
; Coordinates => Sun elev ~42° on 2019/03/20, azim ~180° => use of nenufar_gain.pro (on /cep/lofar/nenufar/pro/general)
a = [-1]
for f=20, 60, 5 do begin
	nenufar_gain, 42, 180, f, ae
	a=[a,ae]
endfor
a = a[1:*]

f = findgen(9)*5.+20
x = poly_fit(f,a,3)
ae = x(0)+x(1)*freq+x(2)*(freq^2)+x(3)*(freq^3)

; calibration curve corrf0 for frequency ramp freq0
T = 60.*(300./freq)^2.55
k = 1.38e3
SEFD = 2.*k*T/ae
corrf0 = m/SEFD
freq0 = freq
save, freq0, corrf0, filename = 'calibration_factor.sav'

END
