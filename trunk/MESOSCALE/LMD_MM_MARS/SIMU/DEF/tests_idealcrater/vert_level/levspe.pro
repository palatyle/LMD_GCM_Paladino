pro levspe


;
; TWEAK PARAM
;
@param.idl


;;
;; HARD PARAM
;;
;!p.multi=[0,2,2]
set_plot, 'ps'
device, file='plot.ps';, /landscape
tinv=5      ;; intervalle plot
hache=10.   ;; scale height
psurf=610.  ;; pression surface -- ne change pas les eta levels

;
; PRELIM
;
param = 'e_vert='+string(nlev,'(I0)')+' ztop='+string(altmax,'(I0)')+' e='+string(epsilon,'(I0)')+' c='+string(elong_cos,'(F4.2)')
;
epsilon = epsilon / (nlev-1) / 100.
;
exposant = !pi/nlev/elong_cos
;
; alpha (plus grand ecart en km) est determine pour que max(altitudes)=altmax
;
alpha =  altmax / ( (nlev-1)/2. - sin(2.*exposant*(nlev-1))/4./exposant + epsilon*(nlev-1)^2/2. )
print, 'alpha (km)', alpha
;
x=[findgen(nlev)]


;
; CALC
;
ecart = alpha - alpha*cos(exposant*x)^2 + epsilon*alpha*x 
altitudes = alpha*x/2. - alpha*sin(2.*exposant*x)/4./exposant + epsilon*alpha*x^2/2. 
logpressions = alog(psurf) - altitudes/hache
pressions=exp(logpressions)
ptop=psurf*exp(-altmax/hache) & print, ptop  
etas=(pressions-ptop)/(psurf-ptop)
etas(nlev-1)=0.


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
print, 'pseudo-altitude (meters)'
pseudo=10.*alog(610./pi)
print, 1000.*10.*alog(610./pi)


;
; PLOT
;
plot, x, etas, title='ETA '+param, xtickinterval=tinv 
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
plot, x, ecart, title='ECART (km, L=theor x=calc)', xtickinterval=tinv
!psym=7 & oplot, x, diff
!psym=0


;
; PRINT
;
;openw, 1, 'levels'
for i=0, nlev-1 do begin
;        printf, 1, etas(i);,',' 
print, etas(i);,','
endfor
close, 1

end
