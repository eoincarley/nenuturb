

; A code to get del/n from delI/I

function calcspeed, fburst, tburst, fh=fh
  
  npoints = ((fburst*1e6/fh)/8980.0)^2.0
  rads = density_to_radius(npoints, model='newkirk')
  
   
  rads = smooth(rads, 10)
  ; First get the drifter speed
  tsec = tburst - tburst[0]
  tsec = smooth(tsec, 10)

  tsec = tsec[20:-20]
  rads = rads[20:-20]
  ;window, 1, xs=500, ys=500
  ;plot, tsec, rads, xtitle='Time (s)', ytitle='Heliocentric distance', /xs, /ys
  
  ;------------------------;
  ;   Fit distance time
  result = linfit(tsec, rads)
  ymodel = result[0] + result[1]*tsec
  set_line_color
  ;oplot, tsec, ymodel, color=5
  speed = result[1]*695.0e6/2.998e8
  print, 'Speed (c)'+string(speed)

  return, speed ; c

END

function calc_di_i, iburst, fburst, npoints=npoints
  
  plot, fburst, iburst, xtitle='Frequency (MHz)', ytitle='Intensity (Arb. Units)', $
      /xs, /ys, charsize=2.0, pos=[0.15, 0.15, 0.9, 0.9]

  ;--------------------------;
  ;
  ;   Smooth and calc di_i
  ;
  sg = savgol(npoints, npoints, 0, 3)
  i_sm2 = convol(iburst, sg, /edge_truncate)

  set_line_color
  oplot, fburst, i_sm2, color=5
  di = mean((iburst-i_sm2)^2)
  _i = mean(iburst)^2

  di_i = sqrt(di/_i)

  print, 'DelI/I (%): '+string(100*di_i)
  return, di_i

END


function calc_dn_n, di_i, vb

  T = 2e6
  kb = 1.38e-23
  me = 9.1e-31
  vth = sqrt((kb*T)/me)
  c = 2.998e8
  vb = vb*c
  print, 'Electron beam velocity (Mm/s): '+string(vb/1e6)
  print, 'Electron thermal velocity (Mm/s): '+string(vth/1e6)

  dn_n = di_i*(vth^2/vb^2)
  print, 'Deln/n (%): '+string(100.0*dn_n)

END



pro calcdeln

  ;---------------------------------------;
  ;
  ; Read HB intensity profile. 
  ; Made during psd_typeIIa_drift.pro
  ;
  print, '----------------------'
  print, ' '
  print, 'Herringbone: '
  print, ' '
  restore, '~/nenuturb/savs/hbfluxprofile.sav'
  iburst = prof 
      
  ;--------------------------;
  ;
  ;     Calc speed.
  ;
  vb = calcspeed(fburst, tburst, fh=2.0)  

  ;--------------------------;
  ;
  ;   Smooth and calc di_i
  ;
  !p.background=255
  !p.color=0
  window, 0, xs=1000, ys=500
  di_i = calc_di_i(iburst, fburst, npoints=200)

  ;----------------------;
  ;
  ;     Get dn_n
  null = calc_dn_n(di_i, vb) 
 

  ;------------------------;
  ;------------------------;
  ;------------------------;
  ;
  ;   Now do the same for
  ;   drifter in type IIa
  print, '----------------------'
  print, ' '
  print, 'Type IIa drifter: '
  print, ' '
  restore, '~/nenuturb/savs/a_drifter_fluxprof.sav'

  ;------------------------;
  ;
  ;     Calc speed.
  ;
  vb = calcspeed(fburst, tburst, fh=1.0)

  ;--------------------------;
  ;
  ;   Smooth and calc di_i
  ;
  window, 1, xs=1000, ys=500
  di_i = calc_di_i(iburst, fburst, npoints=100)

  ;----------------------;
  ;
  ;     Get dn_n
  null = calc_dn_n(di_i, vb)

stop
END
