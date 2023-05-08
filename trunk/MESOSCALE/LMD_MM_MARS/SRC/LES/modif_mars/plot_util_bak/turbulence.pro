pro turbulence

;;----------------------------------------------------
;;----------------------------------------------------
;;----------------------------------------------------


retrieve='true'
;retrieve='false'



;; CALCUL DE LA MOYENNE
;; 18 pas assez 74 trop
meansmooth=37
;meansmooth=18
;meansmooth=74

;; CALCUL DE LA MOYENNE DES FLUX
;smoothampl=74
;smoothampl=2
;smoothampl=10
;smoothampl=18
smoothampl=meansmooth

saveps='false'
;saveps='true'

yu=20

;;----------------------------------------------------
;;----------------------------------------------------
;;----------------------------------------------------


p0=610. & t0=220. & r_cp=1/4.4 & grav=3.72 & R=192.

;
; graphics definition
;
if (saveps eq 'false') then begin
   !p.multi=[0,3,2]
   !P.CHARSIZE=2.
endif else begin
   !p.charthick = 2.0
   !p.thick = 3.0
   !x.thick = 2.0
   !y.thick = 2.0
endelse

;
; input files
;
OPENR, 22, 'input_coord' & READF, 22, lonu & READF, 22, latu & READF, 22, lsu & READF, 22, lctu & CLOSE, 22
OPENR, 23, 'input_more' & READF, 23, hgtu, tsurfu & CLOSE, 23
print, lonu, latu, lsu, lctu, hgtu, tsurfu


;
; retrieve
;
if (retrieve eq 'true') then begin


;
; get fields
;
getcdf, file='u'+string(yu,'(I0)')+'.nc', charvar='U', invar=u
getcdf, file='v'+string(yu,'(I0)')+'.nc', charvar='V', invar=v
getcdf, file='t'+string(yu,'(I0)')+'.nc', charvar='T', invar=t
getcdf, file='w'+string(yu,'(I0)')+'.nc', charvar='W', invar=w
getcdf, file='p'+string(yu,'(I0)')+'.nc', charvar='PTOT', invar=p
	;
	; no vertical staggering for w
	;
	wplus=shift(w,0,0,-1,0)
	w=(w+wplus)/2.
	w=w(*,*,0:n_elements(w(0,0,*,0))-2,*)


ndim=n_elements(t(*,0,0,0))
lat=fltarr(ndim,ndim)
lon=fltarr(ndim,ndim)
  for i=0,ndim-1 do begin
  for j=0,ndim-1 do begin
    lat(i,j)=latu+float(j)*100./59000.
    lon(i,j)=lonu+float(i)*100./59000.
  endfor
  endfor
what_I_plot=0. & what_I_plot2=0. & what_I_plot3=0. & what_I_plot4=0. & what_I_plot5=0. & what_I_plot6=0.
maxwhat_I_plot=0. & maxwhat_I_plot2=0. & maxwhat_I_plot3=0.
minwhat_I_plot=0.


;
; loop on west-east coordinate
;
nloop=n_elements(t(*,0,0,0))
for i=0,nloop-1 do begin

	indp=[i,0]	;;indp=[75,0]

	;
	; 2D field for plot
	;
	t2=reform(t(indp[0],indp[1],*,*))
t2_save=t0+reform(t(indp[0],indp[1],*,*))
treal2_save=t2_save*(reform(p(indp[0],indp[1],*,*))/p0)^r_cp
	u2=reform(u(indp[0],indp[1],*,*))
	v2=reform(v(indp[0],indp[1],*,*))
	w2=reform(w(indp[0],indp[1],*,*))
	nz=n_elements(t2(*,0))
	nt=n_elements(t2(0,*))

nt=nt-2
t2=t2[*,0:nt-1]
u2=u2[*,0:nt-1]
v2=v2[*,0:nt-1]
w2=w2[*,0:nt-1]

        	;;
        	;; local time
        	;;
        	;localtime=21.+100.*findgen(nt)/3700.+lon(indp[0],indp[1])/15.-24. 
                localtime=lctu+100.*findgen(nt)/3700.
		;;print, localtime
		
	;
	; 1. perturbations
	;
	for k=0,nz-1 do begin
		t2(k,*)=t2(k,*)-SMOOTH(reform(t2(k,*)),meansmooth,/EDGE_TRUNCATE)
		u2(k,*)=u2(k,*)-SMOOTH(reform(u2(k,*)),meansmooth,/EDGE_TRUNCATE)
		v2(k,*)=v2(k,*)-SMOOTH(reform(v2(k,*)),meansmooth,/EDGE_TRUNCATE)
		w2(k,*)=w2(k,*);-SMOOTH(reform(w2(k,*)),meansmooth,/EDGE_TRUNCATE)
		;treal2_save(k,*)=treal2_save(k,*)-SMOOTH(reform(treal2_save(k,*)),meansmooth,/EDGE_TRUNCATE)
	endfor
	;
	; 2. nonlinear products
	;
	vertical_eddy_heat_flux=w2*t2
	uprime2=u2^2
	vprime2=v2^2
	wprime2=w2^2
	;
	; 3. Reynolds averaging
	;
	for k=0,nz-1 do begin
		vertical_eddy_heat_flux(k,*)=SMOOTH(reform(vertical_eddy_heat_flux(k,*)),smoothampl,/EDGE_TRUNCATE)
	        uprime2(k,*)=SMOOTH(reform(uprime2(k,*)),smoothampl,/EDGE_TRUNCATE)
	        vprime2(k,*)=SMOOTH(reform(vprime2(k,*)),smoothampl,/EDGE_TRUNCATE)
	        wprime2(k,*)=SMOOTH(reform(wprime2(k,*)),smoothampl,/EDGE_TRUNCATE)
	endfor
	tke=0.5*(uprime2+vprime2+wprime2)


what_I_plot=what_I_plot+vertical_eddy_heat_flux/nloop
what_I_plot2=what_I_plot2+tke/nloop
what_I_plot3=what_I_plot3+wprime2/nloop
what_I_plot4=what_I_plot4+t2_save/nloop
what_I_plot5=what_I_plot5+treal2_save/nloop
if (i eq nloop-1) then  what_I_plot6=t2

;;help, t2
;;help, treal2_save
;treal2_save=treal2_save(0:n_elements(t2(*,0))-1,0:n_elements(t2(0,*))-1)
;if (i eq nloop-1) then  what_I_plot6=treal2_save

yeah=max(vertical_eddy_heat_flux)
yeah2=max(tke)
yeah3=max(wprime2)
yeah4=min(vertical_eddy_heat_flux)

if (yeah gt maxwhat_I_plot) then maxwhat_I_plot=yeah
if (yeah2 gt maxwhat_I_plot2) then maxwhat_I_plot2=yeah2
if (yeah3 gt maxwhat_I_plot3) then maxwhat_I_plot3=yeah3
if (yeah4 lt minwhat_I_plot) then minwhat_I_plot=yeah4

endfor





;
; get model height
;
getcdf, file='ph'+string(yu,'(I0)')+'.nc', charvar='PHTOT', invar=ph
height=reform(ph(*,0,*,*))   
height=total(height,1)/n_elements(height(*,0,0))
height=reform(height)
height=total(height(*,1:nt-1),2)/(nt-1)
height=height/1000./3.72
heightp=height(0:n_elements(height(*))-2) 
	;
	; altitude above ground
	;	
        heightp=heightp-hgtu/1000.


;
; save
;
save, what_I_plot, $
	what_I_plot2, $
	what_I_plot3, $
	what_I_plot4, $
	what_I_plot5, $
	what_I_plot6, $
	localtime, $
	heightp, $
	lat, $
	lon, $
	hgtu, $
	filename='turbulence.dat' 

endif else begin
	restore, filename='turbulence.dat'	
endelse	

help, what_I_plot
help, localtime
help, heightp


latwrite=string(lat(0,yu),'(F5.1)')




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (saveps eq 'true') then begin
  set_plot, 'ps' & device, filename='plot/TKE.ps'
endif

lev=findgen(15)+1.
contour, $
        transpose(what_I_plot2), $
        localtime, $
        heightp, $
        xtitle='Local Time (h)', $
        xrange=[8.,18.], $
        xtickinterval=1., $
        ytitle='Altitude above surface (km)', $
        yrange=[0.,10.], $
        ytickinterval=1., $
        title="Turbulent Kinetic Energy 0.5[<u'!U2!N>+<v'!U2!N>+<w'!U2!N>] (m!U2!N.s!U-2!N)", $
        levels=lev, $
        c_labels=findgen(n_elements(lev))*0.+1., $
	subtitle='zonal average at lat. '+latwrite+' (max value is '+string(maxwhat_I_plot2,'(F4.1)')+' m!U2!N.s!U-2!N)'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (saveps eq 'true') then begin
  device, /close
  set_plot, 'ps' & device, filename='plot/zTKE.ps'
endif

lev=findgen(15)+1.
contour, $
        0.5*transpose(what_I_plot3), $
        localtime, $
        heightp, $
        xtitle='Local Time (h)', $
        xrange=[8.,18.], $
        xtickinterval=1., $
        ytitle='Altitude above surface (km)', $
        yrange=[0.,10.], $
        ytickinterval=1., $
        title="Vertical TKE 0.5<w'!U2!N> (m!U2!N.s!U-2!N)", $
        levels=lev, $
        c_labels=findgen(n_elements(lev))*0.+1., $
        subtitle='zonal average at lat. '+latwrite+' (max value is '+string(0.5*maxwhat_I_plot3,'(F4.1)')+' m!U2!N.s!U-2!N)'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;




if (saveps eq 'true') then begin
  device, /close
  set_plot, 'ps' & device, filename='plot/HF.ps'
endif

;lev=[-0.15,-0.05,-0.02,0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95]
lev = -1. + findgen(20)/10.
contour, $
        transpose(what_I_plot), $
        localtime, $
        heightp, $ 
	xtitle='Local Time (h)', $
	xrange=[8.,18.], $
	xtickinterval=1., $
	ytitle='Altitude above surface (km)', $
	yrange=[0.,10.], $
	ytickinterval=1., $
	title="Vertical Eddy Heat Flux <w'!7h!3'> (K.m.s!U-1!N)", $
	levels=lev, $	
	C_LINESTYLE = (lev LT 0.0), $
	c_labels=findgen(n_elements(lev))*0.+1., $
        subtitle='zonal average at lat. '+latwrite+' (min/max values are '+string(minwhat_I_plot,'(F5.1)')+'/'+string(maxwhat_I_plot,'(F3.1)')+' K.m.s!U-1!N)'





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
goto, no_pert_temp
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (saveps eq 'true') then begin
  device, /close
  set_plot, 'ps' & device, filename='plot/pert_temp.ps'
endif

lev=[-4.,-2.,-1.,-0.5,0.5,1.,2.,4.]
contour, $
	transpose(what_I_plot6), $
        localtime, $
        heightp, $
        xtitle='Local Time (h)', $
        xrange=[13.5,14.], $
        xtickinterval=0.1, $
        ytitle='Altitude above surface (km)', $
        yrange=[0.,0.5], $
        ytickinterval=0.1, $
        title="Potential Temperature Perturbation !7h!3' (K)", $
        levels=lev, $
;nlevels=20, $
        C_LINESTYLE = (lev LT 0.0), $
        c_labels=findgen(n_elements(lev))*0.+1., $
;/cell_fill, $
        subtitle='latitude '+latwrite+' longitude '+string(lon(nloop-1,0),'(F6.1)')
;xyouts, 14.007, heightp(0), '1'
;xyouts, 14.016, heightp(1), '2'
;xyouts, 14.007, heightp(2), '3'
;xyouts, 14.007, heightp(3), '4'
;xyouts, 14.007, heightp(4), '5'
;xyouts, 14.007, heightp(5), '6'
;xyouts, 14.007, heightp(6), '7'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
no_pert_temp:
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;





if (saveps eq 'true') then begin
  device, /close
  set_plot, 'ps' & device, filename='plot/temp_pot.ps'
endif

user_lt=10. & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt))))
plot, $
        what_I_plot4(*,yeah(0)), $
        heightp,$
	xtitle='Potential Temperature (K)',$
;	xrange=[215,240], $
        xtickinterval=5., $
        ytitle='Altitude above surface (km)', $
	yrange=[0.,10.], $
        ytickinterval=1., $
	title="Potential Temperature !7h!3 profile (K)", $
        subtitle='zonal average at lat. '+latwrite
user_lt=12. & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt))))
        oplot, what_I_plot4(*,yeah(0)), heightp, linestyle=1
user_lt=14. & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt))))
        oplot, what_I_plot4(*,yeah(0)), heightp
user_lt=16. & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt))))
        oplot, what_I_plot4(*,yeah(0)), heightp, linestyle=1
user_lt=18. & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt))))
        oplot, what_I_plot4(*,yeah(0)), heightp

;xyouts, 217.3, 0.5, '10:00'
;xyouts, 225.8, 1.5, '12:00'
;xyouts, 231.8, 3.0, '14:00'
;xyouts, 234.8, 0.7, '16:00'
;xyouts, 232.3, 2.0, '18:00'


if (saveps eq 'true') then begin
  device, /close
  set_plot, 'ps' & device, filename='plot/temp.ps'
endif

user_lt=14.+55./60. & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt))))
plot, $
        what_I_plot5(*,yeah(0)), $
        heightp,$
        xtitle='Temperature (K)',$
        xrange=[210,250], $
        xtickinterval=5., $
        ytitle='Altitude above surface (km)', $
        yrange=[0.,2.], $
        ytickinterval=0.2, $
        title="LMD LES Near-surface T profile (K)", $
        subtitle='zonal average at lat. '+latwrite, $
	linestyle=2	
user_lt=10.+24./60. & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt))))
	oplot, what_I_plot5(*,yeah(0)), heightp
user_lt=10.+59./60. & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt))))
	oplot, what_I_plot5(*,yeah(0)), heightp, linestyle=2
user_lt=11.+59./60. & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt))))
	oplot, what_I_plot5(*,yeah(0)), heightp
user_lt=15.+22./60. & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt))))
	oplot, what_I_plot5(*,yeah(0)), heightp
user_lt=17.+13./60. & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt))))
        oplot, what_I_plot5(*,yeah(0)), heightp, linestyle=2

;xyouts, 221,   0.2,  '10:24'
;xyouts, 226.5,   0.4,  '10:59'
;xyouts, 231.8, 0.2,  '11:59'
;xyouts, 233,   0.8,  '14:55'
;xyouts, 241,   0.1,  '15:22'
;xyouts, 236.5, 0.05, '17:13'
 

;;
;xyouts, 250.05, heightp(0), '1'
;xyouts, 250.45, heightp(1), '2'
;xyouts, 251.00, heightp(2), '3'
;xyouts, 250.45, heightp(3), '4'
;xyouts, 250.45, heightp(4), '5'
;xyouts, 250.45, heightp(5), '6'
;xyouts, 250.45, heightp(6), '7'
;xyouts, 250.45, heightp(7), '8'
;xyouts, 250.45, heightp(8), '9'
;xyouts, 250.25, heightp(9), '10'
;xyouts, 250.25, heightp(10), '11'
;xyouts, 250.25, heightp(11), '12'
;xyouts, 250.25, heightp(12), '13'
;xyouts, 250.25, heightp(13), '14'
;xyouts, 250.25, heightp(14), '15'
;xyouts, 250.25, heightp(15), '16'



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;goto, no_staticstab
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (saveps eq 'true') then begin
  device, /close
  set_plot, 'ps' & device, filename='plot/staticstab.ps'
endif

user_lt=17. & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt))))

   ;;; recompute height @ given local time 
   caca=heightp+hgtu/1000.
   height=reform(ph(*,0,*,*))
   height=total(height,1)/n_elements(height(*,0,0))
   height=reform(height)
   height=reform(height(*,yeah(0)))
   height=height/1000./3.72
   heightp=height(0:n_elements(height(*))-2) 
   print, 'new minus old', heightp - caca
   ;;; recompute height @ given local time 

        press=reform(p(*,0,*,*))
        press=total(press,1)/n_elements(press(*,0,0))
        press=reform(press)
        press=reform(press(*,yeah(0)))
        press=press(0:n_elements(press(*))-2)

staticstab = DERIV( heightp , reform(what_I_plot5(*,yeah(0))) ) + 1000.*3.72/844.6  ;; 4.9 dans Hinson
   staticstab = staticstab(0:n_elements(staticstab)-5)
   heightp = heightp(0:n_elements(heightp)-5)

plot, $
        staticstab, $
        heightp,$
        xtitle='Static stability (K/km)',$
        xrange=[-1.,5.], $
        xtickinterval=0.5, $
        ytitle='Altitude above surface (km)', $
        yrange=[min(heightp),max(heightp)], $
        ytickinterval=1., $
        title="LMD LES Static stability @ LT 17h (K/km)", $
        subtitle='zonal average at lat. '+latwrite
oplot, $
        staticstab, $
        heightp,$
        psym=5
oplot, $
       findgen(n_elements(heightp))*0. + 1.,$
       heightp,$
       linestyle=2
oplot, $
       findgen(n_elements(heightp))*0. + 2.,$
       heightp,$
       linestyle=2

;;;;;;;;;;
w = where((staticstab gt 1.5) and ((heightp- hgtu/1000.) gt 1.))
t_top_plus = what_I_plot5(w(0),yeah(0))
z_top_plus = heightp(w(0))
p_top_plus = press(w(0))
t_top_moins = what_I_plot5(w(0)-1,yeah(0))
z_top_moins = heightp(w(0)-1)
p_top_moins = press(w(0)-1)
pbl_depth = (z_top_plus*alog(p_top_plus) + z_top_moins*alog(p_top_moins))/(alog(p_top_plus) + alog(p_top_moins)) - hgtu/1000.
xyouts, 3., 1.5 + (max(heightp) + min(heightp)) / 3., 'Ls = '+string(lsu,'(F5.1)')+'!Uo!N', CHARSIZE=1
xyouts, 3., 1. + (max(heightp) + min(heightp)) / 3., 'Lat = '+string(latu,'(F5.1)')+'!Uo!N', CHARSIZE=1
xyouts, 3., 0.5 + (max(heightp) + min(heightp)) / 3., 'LonE = '+string(lonu,'(F6.1)')+'!Uo!N', CHARSIZE=1
xyouts, 3., (max(heightp) + min(heightp)) / 3., 'T!Dt!N = '+string(t_top_plus,'(I0)')+'/'+string(t_top_moins,'(I0)')+' K ', CHARSIZE=1
xyouts, 3., -0.5 + (max(heightp) + min(heightp)) / 3., 'p!Dt!N = '+string(p_top_plus,'(I0)')+'/'+string(p_top_moins,'(I0)')+' Pa', CHARSIZE=1
xyouts, 3., -1. + (max(heightp) + min(heightp)) / 3., 'z!Dt!N = '+string(z_top_plus,'(F4.1)')+'/'+string(z_top_moins,'(F4.1)')+' km', CHARSIZE=1
xyouts, 3., -1.5 + (max(heightp) + min(heightp)) / 3., 'z!Ds!N = '+string(hgtu/1000.,'(F4.1)')+' km', CHARSIZE=1
xyouts, 3., -2. + (max(heightp) + min(heightp)) / 3., 'D = '+string(pbl_depth,'(F4.1)')+' km', CHARSIZE=1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
no_staticstab:
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (saveps eq 'true') then begin
  device, /close
endif
stop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


user_lt=13.
yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt))))

user_h=0.1
walt=where(abs(heightp-user_h) eq (min(abs(heightp-user_h))))

;mapfield=reform(t[*,*,walt(0),w(0)])
;help, mapfield
;contour, mapfield

section=reform(w[*,0,*,yeah(0)])

section=reform(w[*,0,*,160])

lev=[-12.,-8.,-4.,4.,8.,12.]
lev=[-5.,-4.,-3.,-2.,-1.,1.,2.,3.,4.,5.]
contour, $
	section, $
	(lon(*,0)-lon(0,0))*59., $
	heightp, $
	levels=lev , $
	C_LINESTYLE = (lev LT 0.0)

device, /close

end
