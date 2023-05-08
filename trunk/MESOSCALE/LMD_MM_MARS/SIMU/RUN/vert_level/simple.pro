pro simple, z_scale


nlev=101 
;
; z-scale doit etre suffisamment bas pour que l avant-dernier niveau soit assez haut
; mais pas trop sinon instable (il faut avoir des niveaux espaces en haut semble-t-il)
;
z_scale = .25
;z_scale = .15 marche pas (niveaux pas assez espaces en haut ?)
ptop = 0.2
npbl = 5
	;z_scale = .15
	;ptop = 0.05 ;marche pas

z_scale = .20 & ptop = 0.1 & npbl = 6
z_scale = .25 & ptop = 0.5 & npbl = 6 & nlev=61  ;; HOLDEN NESTS

psurf=610.

pbl = [1.0000,0.9995,0.9980,0.9950,0.9850,0.9700,0.9400,0.9000] 
x=[findgen(nlev)]
nlev = nlev - npbl + 1 

etas=[findgen(nlev)]
FOR k=0, nlev-1 DO BEGIN
     etas(k) = (exp(-float(k)/float(nlev-1)/z_scale) - exp(-1./z_scale)) / (1.-exp(-1./z_scale)) 
ENDFOR
etas(nlev-1)=0.

etas = [pbl(0:npbl-1),etas(1:nlev-1)]
nlev = nlev + npbl - 1 
print, nlev

;pressions = ptop + etas * (psurf-ptop)
;zzz = 10. * alog ( psurf / pressions )
;print, zzz

;
; HARD PARAM
;
!p.multi=[0,2,2]
;set_plot, 'ps'
;device, file='plot.ps', /landscape



;
; PRINT
;
print, 'eta levels'
print, etas
;
print, 'pressure'
pi=etas*(psurf-ptop)+ptop
print, pi
print, 'ptop', ptop
;
print, 'pseudo-altitude'
pseudo=10.*alog(610./pi)
print, 1000.*10.*alog(610./pi)


;
; PLOT
;
plot, x, etas, title='ETA', xtickinterval=tinv 
!psym=7 & oplot, x, etas
!psym=0
plot, x, pseudo, title='ALTITUDE (km)', xtickinterval=tinv
!psym=7 & oplot, x, pseudo 
!psym=0
plot, x, pi, title='PRESSURE (Pa)', ylog=1, xtickinterval=tinv
!psym=7 & oplot, x, pi 
!psym=0
;
diff=pseudo - shift(pseudo,-1) & diff=-diff(0:nlev-2)
;
plot, x, diff, title='ECART (km, L=theor x=calc)', xtickinterval=tinv
!psym=7 & oplot, x, diff
!psym=0

;stop

;
; PRINT
;
openw, 1, 'levels'
for i=0, nlev-1 do begin
        printf, 1, etas(i);,',' 
endfor
close, 1

end
