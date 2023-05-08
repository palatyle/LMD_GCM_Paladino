function getget, name, charvar, count=count, offset=offset, stride=stride, coordmean=coordmean, anomaly=anomaly
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; name: char
; charvar: char
; count: 4 elements array   ; The dimensions of the sub-image
; offset: 4 elements array  ; A variable that contains the offset for the sub-image
; stride: 4 elements array  ; Create a variable to be used as a value for the STRIDE keyword
; coordmean: 1 to 4 elements array (dimension 1 is 1)
; anomaly: a non-null element that would be filled with mean
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


nosave = 1  ;; comment to get an IDL data file
notime = 1  ;; comment to get mention of time taken
if (n_elements(notime) eq 0) then timetime = SYSTIME(1)

;
; MINIMAL ARGUMENTS
;
if (n_elements(name) eq 0) then stop
if (n_elements(charvar) eq 0) then stop

;
; IDL DATA DIRECTORY
;
name_idl = charvar+name+'.idl'
  if (n_elements(coordmean) ne 0) then begin
    nc = n_elements(coordmean)
    for i=0,nc-1 do begin
     name_idl = string(coordmean(i),'(I0)')+name_idl   
    endfor
  endif
name_idl = './idl_dat/'+name_idl

;
; TEST IF FILE EXISTS
;
openr,unit,name_idl,/get_lun,error=err

;
; IF IDL FILE EXISTS, THAT IS RATHER SIMPLE
;
IF (err EQ 0) THEN BEGIN

   print, 'read from file'
   free_lun,unit
   restore, filename=name_idl

;
; IF IDL FILE DOES NOT EXIST, WELL, READ NETCDF
;
ENDIF ELSE BEGIN

  ;
  ; check save directory
  ;
  if (n_elements(nosave) eq 0) then SPAWN, 'mkdir -p idl_dat'  

  ;
  ; open file
  ;
  id=ncdf_open(name) 

  ;
  ; get dimensions
  ;
  var = intarr(4)
  NCDF_DIMINQ, id, NCDF_DIMID(id, 'west_east'    ), toto, titi & var(0) = titi
  NCDF_DIMINQ, id, NCDF_DIMID(id, 'south_north'  ), toto, titi & var(1) = titi 
  NCDF_DIMINQ, id, NCDF_DIMID(id, 'bottom_top'   ), toto, titi & var(2) = titi
  NCDF_DIMINQ, id, NCDF_DIMID(id, 'Time'         ), toto, titi & var(3) = titi 
  ;print, var

  ;
  ; define default data windows
  ;
  if (n_elements(offset) le 3) then offset=[0,0,0,0]
  if (n_elements(stride) le 3) then stride=[1,1,1,1]
  if (n_elements(count) le 3) then count=intarr(4)
  w = where(count eq 0) & if (w(0) ne -1) then count[w] = var[w] / stride[w]

  ;
  ; read variable
  ;
  varid=ncdf_varid(id,charvar) ;& print, 'get '+charvar 
  if (n_elements(notime) eq 0) then timetime = SYSTIME(1) ;& help, /memory
  ncdf_varget, id, varid, invar, count=count, offset=offset, stride=stride & if (n_elements(notime) eq 0) then print, SYSTIME(1) - timetime, ' s' ;& help, /memory
  ncdf_close, id

  ;
  ; OPTIONAL: ANOMALY (value is returned in the keyword)
  ; 
  if (n_elements(anomaly) ne 0) then begin
    anomaly=TOTAL(TOTAL(invar,1),1)/float(var[0])/float(var[1])
    for i=0,var[0]-1 do for j=0,var[1]-1 do invar(i,j,*,*) = TEMPORARY(invar(i,j,*,*)) - anomaly  ; coute cher
  endif

  ;
  ; OPTIONAL: MEAN 
  ;
  if (n_elements(coordmean) ne 0) then begin
  for i=0,nc-1 do begin
   print, 'mean over dimension '+string(coordmean(i),'(I0)')+' - '+string(var(coordmean(i)),'(I0)')+' elements'
   invar=total(TEMPORARY(invar),coordmean(i))/var(coordmean(i))
  endfor
  endif

  ;
  ; SAVE
  ;
  if (n_elements(nosave) eq 0) then save, invar, filename=name_idl

ENDELSE

if (n_elements(notime) eq 0) then print, SYSTIME(1) - timetime, ' s'

;
; RETURN FIELD
;
return, invar

end
