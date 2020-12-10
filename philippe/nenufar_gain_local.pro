; ------------------------------------------------------------------------------------------
  pro NenuFAR_gain_local, elev,azim, freq, ae, g_Array,g_NE_SW,g_NW_SE, $
			accurate=accurate, plot=plot, verbose=verbose, help=help, hlp=hlp
; ------------------------------------------------------------------------------------------
; Theoretical gain of NenuFAR for a pointing at elev,azim
;
; INPUTS
; elev & azim in degrees, freq in MHz
;
; OUTPUTS
; Effective area  ae in m^2
; Gains g_Array (array term), g_NE_SW,g_NW_SE (array & antenna terms)  in linear values

if keyword_set(hlp) then begin
  print
  print,'NenuFAR_gain, [IN]   elev,azim, freq, $'
  print,'              [OUT]  ae, g_Array,g_NE_SW,g_NW_SE, $'
  print,'              [KEY]  /accurate, /plot, /verbose'
  print
  return
endif

if keyword_set(help) then begin
  print
  print,' Computes the theoretical gain of present NenuFAR configuration for a pointing at elev,azim'
  print
  print,' INPUTS'
  print,' elev, azim = elevation, azimuth (deg.)'
  print,' freq = frequency (MHz)'
  print
  print,' KEYWORDS'
  print,' /accurate => accurate (<0.1%) and longer calculation of Aeff'
  print,' /plot     => display Mini-Arrays & antennas layout'
  print,' /verbose  => display results'
  print
  print,' OUTPUTS'
  print,' g_Array = array gain                 [linear value]'
  print,' g_NE_SW = array & NE_SW antenna gain [linear value]'
  print,' g_NW_SE = array & NW_SE antenna gain [linear value]'
  print
  return
endif

; Antennas in a Mini-Array
restore,'/cep/lofar/nenufar/pro/general/antennas_arrays/amrlss.sav'

; Coordinates of Mini-Array centers
r=read_csv('/cep/lofar/nenufar/pro/general/antennas_arrays/Coordonnees_MR00-MR55_20191105_OK_LD.csv')
nMR=n_elements(r.field01)
MRc=transpose([[r.field08],[r.field09],[r.field10]])
MRc=Mrc-rebin(reform(mean(MRc,dim=2),3,1),3,nMR)

; Mini-Array Rotations
MRrot=r.field05

antennas=fltarr(3,19*nMR)

if KEYWORD_SET(PLOT) then begin
  window,0,xs=1200,ys=800
  plot,MRc(0,*),MRc(1,*),psym=4,/iso,xtit='m',ytit='m',tit='NenuFAR MA 00-'+strtrim(fix(nMR-1),2),xra=[-250,250],/xsty,charsize=1.5
  xyouts,60.,205,'Projection perp. to (El,Az) = '+strtrim(fix(elev),2)+'!Uo!N, '+strtrim(fix(azim),2)+'!Uo!N',color=255,charsize=1.5
  for i=0,nMR-1 do begin
    MRrot_i=MRrot(i)*!dtor
    Mrot=[[cos(MRrot_i),-sin(MRrot_i)],[sin(MRrot_i),cos(MRrot_i)]]
    a_i=fltarr(3,19)
    a_i(0:1,*)=Mrot#a(0:1,*)
    a_i=a_i+rebin(reform(MRc(*,i),3,1),3,19)
    antennas(*,i*19:(i+1)*19-1)=a_i
    oplot,a_i(0,*),a_i(1,*),psym=1
    xyouts,MRc(0,i),MRc(1,i)-15,strtrim(i,2)
  endfor
  wait, 0.5
endif else begin
  for i=0,nMR-1 do begin
    MRrot_i=MRrot(i)*!dtor
    Mrot=[[cos(MRrot_i),-sin(MRrot_i)],[sin(MRrot_i),cos(MRrot_i)]]
    a_i=fltarr(3,19)
    a_i(0:1,*)=Mrot#a(0:1,*)
    a_i=a_i+rebin(reform(MRc(*,i),3,1),3,19)
    antennas(*,i*19:(i+1)*19-1)=a_i
  endfor
endelse

ph=azim*!dtor
th=elev*!dtor
u=[sin(ph)*cos(th),cos(ph)*cos(th),sin(th)]
for i=0,19*nMR-1 do antennas(*,i)=antennas(*,i)-total(antennas(*,i)*u)*u

if KEYWORD_SET(PLOT) then begin
  oplot,antennas(0,*),antennas(1,*),psym=1,color=255
  oplot,[0,0],[-300,300],line=1
  oplot,[-300,300],[0,0],line=1
  wait, 0.5
endif

lambda=300./freq
aa=antennas/lambda
lg=0.2
if KEYWORD_SET(ACCURATE) then lg=0.05
if KEYWORD_SET(PLOT) then window,1,xs=1200,ys=800
ARRAY_AREA, 19*nMR, aa, 3., ae, lgrid=lg, PLOT=PLOT
g_Array=ae*4*!pi
READ_CONV_GAIN,'/cep/lofar/nenufar/pro/general/antennas_arrays/','NE_SW',freq, g_NE_SW
g_NE_SW=g_NE_SW(round(azim),round(elev))
READ_CONV_GAIN,'/cep/lofar/nenufar/pro/general/antennas_arrays/','NW_SE',freq, g_NW_SE
g_NW_SE=g_NW_SE(round(azim),round(elev))
ae=ae*lambda^2

if KEYWORD_SET(VERBOSE) then begin
  print
  print,'Elevation = ',strtrim(elev,2),'°,   Azimuth = ',strtrim(azim,2),'°,   Frequency = ',strtrim(freq,2),' MHz'
  print,'Ae Array           = ',ae,' m^2'
  print,'Gain Array         = ',g_Array,' = ',10.*alog10(g_Array),' dB'
  print,'Gain Array * NE_SW = ',g_Array*g_NE_SW,' = ',10.*alog10(g_Array*g_NE_SW),' dB'
  print,'Gain Array * NW_SE = ',g_Array*g_NW_SE,' = ',10.*alog10(g_Array*g_NW_SE),' dB'
  print
endif

return
end
