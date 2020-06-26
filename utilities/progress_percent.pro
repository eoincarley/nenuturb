pro progress_percent,i,start_i,end_i

; copyright  ©  2015 Hamish 'BFTD' Reid

; Pro to read out the progress in a loop assuming integer start and end indexes.
; The assumption is also increments of += 1 to i

; INPUTS

; i - the indexing number of the loop
; start_i - the start index in the loop
; end_i - the end index in the loop

; SAMPLE CALL

; progress,i,start_i,end_i



st = floor(start_i)
et = floor(end_i)

if i eq st then print,'0%',format='(A,$)'
if ((i-st)*1. mod ((et-st)/10.)) lt (((i-st)*1.-1) mod ((et-st)/10.)) then begin
	if floor((i-st)*100./(et-st)) lt 100 then print,floor((i-st)*100./(et-st)),'%',format='(I2,A,$)'
	if floor((i-st)*100./(et-st)) ge 100 and i eq et then begin
		print,floor((i-st)*100./(et-st)),'%',format='(I3,A)'
		return
	endif
endif else if i ge et then begin
	print,'100%',format='(A)'
	return
endif else if (et-st) gt 100 then if ((i-st)*1. mod ((et-st)/50.)) lt (((i-st)*1.-1) mod ((et-st)/50.)) then begin
	print,'.',format='(A,$)'
endif

if (et-st) lt 100 then	print,'.',format='(A,$)'

end
