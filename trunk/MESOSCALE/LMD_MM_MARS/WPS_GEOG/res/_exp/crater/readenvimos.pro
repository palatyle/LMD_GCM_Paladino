pro readenvimos,hdrfile,mos,lines,samples,easternmostlongitude,westernmostlongitude,$
                minimumlatitude,maximumlatitude



datfile  = strmid(hdrfile,0,strlen(hdrfile) - 3) + 'dat'

openr,1,hdrfile
a = ''
readf,1,a
readf,1,a
nband = 1L
vue = 0
while(not(eof(1))) do begin
  readf,1,a
  b = strmid(a,0,3)
  print,b
  if b eq 'sam' then begin
    c = strsplit(a,/extract)
    samples = long(c[2])
  endif
  if b eq 'lin' then begin
    c = strsplit(a,/extract)
    lines = long(c[2])
  endif
  if b eq 'map' then begin
    c = strsplit(a,/extract)
    westernmostlongitude = float(c[7])
    minimumlatitude = float(c[8])
    deltalat = float(c[9])
    deltalon = float(c[10])
    easternmostlongitude = westernmostlongitude + deltalon*(samples - 1)
    maximumlatitude = minimumlatitude + deltalat*(lines - 1)
  endif
  if (b eq 'ban' and vue eq 0) then begin
    c = strsplit(a,/extract)
    nband = long(c[2])
    vue = 1
  endif
  if b eq 'dat' then begin
    c = strsplit(a,/extract)
    print,c,c[3]
    data_type = fix(c[3])
  endif
endwhile
close,1


openr,1,datfile
if nband eq 1L then begin
  if data_type eq 1 then mos = bytarr(samples,lines)
  if data_type eq 2 then mos = intarr(samples,lines)
  if data_type eq 4 then mos = fltarr(samples,lines)
  help,mos
endif

if nband gt 1L then begin
  if data_type eq 1 then mos = bytarr(samples,lines,nband)
  if data_type eq 2 then mos = intarr(samples,lines,nband)
  if data_type eq 4 then mos = fltarr(samples,lines,nband)



endif

readu,1,mos



close,1

end