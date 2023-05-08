pro getcdf, $
	file=file, $
	charvar=charvar, $
	invar=invar	

T = SYSTIME(1)

cdfid = ncdf_open(file)
varid=ncdf_varid(cdfid,charvar)
ncdf_varget, cdfid, varid, invar

print, 'got '+charvar+' in '+file+' within ', SYSTIME(1) - T, ' seconds'


end
