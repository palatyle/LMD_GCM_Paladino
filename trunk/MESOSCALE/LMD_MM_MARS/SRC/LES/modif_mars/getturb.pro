pro getturb, saveps=saveps



;----------------------------------
; USE: getturb
;      getturb, saveps='false'     
;----------------------------------


history_interval_s = 100.
smoothampl=3700/history_interval_s
smoothampl=0.

;
; constantes
;
p0=610. & t0=220. & r_cp=1/4.4 & grav=3.72 & R=192.

;
; graphics definition
;
if (n_elements(saveps) eq 0) then saveps='true'
if (saveps eq 'false') then begin
   ;!p.multi=[0,3,2] 
   !P.CHARSIZE=2.
   WINDOW, /PIXMAP & WDELETE & DEVICE,BYPASS_TRANSLATION=0,DECOMPOSED=0,RETAIN=2
endif else begin
   PREF_SET, 'IDL_PATH', '/home/spiga/Save/SOURCES/IDL/fsc_psconfig:<IDL_DEFAULT>', /COMMIT
endelse

;
; retrieve fields
;
openr,unit,'getturb.dat',/get_lun,error=err
IF (err ne 0) THEN BEGIN

;
; input files
;
OPENR, 22, 'input_coord' & READF, 22, lonu & READF, 22, latu & READF, 22, lsu & READF, 22, lctu & CLOSE, 22
OPENR, 23, 'input_more' & READF, 23, hgtu, tsurfu & CLOSE, 23

;
; get fields
;
domain='d01' & filesWRF = FindFile('wrfout_'+domain+'_????-??-??_??:??:??') & nf=n_elements(filesWRF)

;
; get dimensions
;
id=ncdf_open(filesWRF(0))
NCDF_DIMINQ, id, NCDF_DIMID(id, 'west_east'    ), toto, nx & NCDF_DIMINQ, id, NCDF_DIMID(id, 'south_north'  ), toto, ny
NCDF_DIMINQ, id, NCDF_DIMID(id, 'bottom_top'   ), toto, nz & NCDF_DIMINQ, id, NCDF_DIMID(id, 'Time'         ), toto, nt
NCDF_CLOSE, id 

;
; prepare loop
;
nloop1 = nf-1 & nloop2 = nt & yeye = 0 
localtime = lctu + history_interval_s*findgen(nloop1*nloop2)/3700.
wt  = fltarr(nz,nloop1*nloop2) & tke = fltarr(nz,nloop1*nloop2) & ztke = fltarr(nz,nloop1*nloop2) & t = fltarr(nz,nloop1*nloop2)
p = fltarr(nz) & ph = fltarr(nz) & pht = fltarr(nz,nloop1*nloop2) & pt = fltarr(nz,nloop1*nloop2) & stst = fltarr(nz,nloop1*nloop2)

;
; loop loop
;
for loop  = 0, nloop1-1 do begin
                                          timetime = SYSTIME(1)
for loop2 = 0, nloop2-1 do begin

   ; t' = t - <t>
   ; ------------
 yeyeye = 1.
 tprime = getget(filesWRF(loop), 'T', anomaly=yeyeye, count=[0,0,0,1], offset=[0,0,0,loop2])  ;; t' = t - <t>
 t(*,yeye) = t0 + TEMPORARY(yeyeye)
   ; w' = w   
   ; ------
 wprime  = getget(filesWRF(loop), 'W', count=[0,0,0,1], offset=[0,0,0,loop2])     
   ; tke = 0.5 ( <u'^2> + <v'^2> + <w'^2> ) ; u' = u ; v' = v  
   ; --------------------------------------------------------
 ztke(*,yeye) = 0.5 * TOTAL(TOTAL(wprime^2,1),1) / float(nx) / float(ny)
 tke(*,yeye) = ztke(*,yeye) + $
               0.5 * ( $
        TOTAL(TOTAL(getget(filesWRF(loop), 'U', count=[0,0,0,1], offset=[0,0,0,loop2])^2,1),1) + $  
        TOTAL(TOTAL(getget(filesWRF(loop), 'V', count=[0,0,0,1], offset=[0,0,0,loop2])^2,1),1)   $  
        ) / float(nx) / float(ny)  
   ; <w't'>
   ; ------
 wt(*,yeye)  = TOTAL(TOTAL(TEMPORARY(tprime)  * TEMPORARY(wprime),1),1) / float(nx) / float(ny)  
   ; p & ph
   ; ------
 if (loop + loop2 eq 0) then nloopbeware = nloop2-1 else nloopbeware = nloop2  ; 1ere valeur vaut 0
 pht(*,yeye) = TOTAL(TOTAL(getget(filesWRF(loop), 'PHTOT',  count=[0,0,0,1], offset=[0,0,0,loop2]),1),1) / float(nx) / float(ny) / 1000. / 3.72
 pt(*,yeye)  = TOTAL(TOTAL(getget(filesWRF(loop), 'PTOT' ,  count=[0,0,0,1], offset=[0,0,0,loop2]),1),1) / float(nx) / float(ny)
 ph = TEMPORARY(ph) + pht(*,yeye) / nloopbeware / nloop1
 p  = TEMPORARY(p ) + pt(*,yeye) / nloop2 / nloop1
   ; static stability
   ; ----------------
 stst(*,yeye) = DERIV( reform(pht(*,yeye)) - hgtu/1000. , reform(t(*,yeye)) * ( reform(pt(*,yeye)) /p0 )^r_cp ) + 1000.*grav / (R / r_cp)  ;; 4.9 dans Hinson
   ; loop count
   ; ---------
 yeye = TEMPORARY(yeye) + 1 

endfor
                                          print, 'file '+string(loop+1,'(I0)'), SYSTIME(1) - timetime, ' s'
endfor
h  = TEMPORARY(ph)  - hgtu/1000.  ;; altitude above ground
ht = TEMPORARY(pht) - hgtu/1000.  
;
; save
;
save, wt, tke, ztke, h, ht, t, p, pt, stst, localtime, filename='getturb.dat'

ENDIF ELSE BEGIN

print, 'OK, file is here'
restore, filename='getturb.dat'

ENDELSE

;
; smooth smooth
;
wt  = SMOOTH(TEMPORARY(wt),  [0,smoothampl], /EDGE_TRUNCATE)
tke = SMOOTH(TEMPORARY(tke), [0,smoothampl], /EDGE_TRUNCATE)
ztke = SMOOTH(TEMPORARY(ztke), [0,smoothampl], /EDGE_TRUNCATE)


;*******************;
;*******************;
; PLOTS PLOTS PLOTS ;
;*******************;
;*******************;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
zefield       =  transpose(tke) 
zex           =  localtime 
zey           =  h
set_name      =  'TKE.ps'
set_title     =  "Turbulent Kinetic Energy 0.5[<u'!U2!N>+<v'!U2!N>+<w'!U2!N>] (m!U2!N.s!U-2!N)"
set_titlex    =  'Local Time (h)'
set_titley    =  'Altitude above surface (km)'
set_subtitle  =  'Mean over the simulation domain'
set_xrange    =  [8.,18.]
set_yrange    =  [0.,10.]
set_tickx     =  1.
set_ticky     =  1.
minval        =  0.
maxval        =  20. ;15.
nlev          =  maxval-minval
pal           =  22
rrr           =  'no'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (saveps eq 'true') then PS_Start, FILENAME=set_name
!P.Charsize = 1.2
;; 0. levels
lev = minval + (maxval-minval)*findgen(nlev+1)/float(nlev) & if (minval ne 0.) then lev = lev[where(lev ne 0.)]
;;; 1. background
loadct, 0 & contour, /NODATA, zefield, zex, zey, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0 
;; 2. color field
loadct, pal & if (rrr eq 'yes') then TVLCT, r, g, b, /Get & if (rrr eq 'yes') then TVLCT, Reverse(r), Reverse(g), Reverse(b) 
            contour, zefield, zex, zey, levels=lev, c_labels=findgen(n_elements(lev))*0.+1., /overplot, /cell_fill
;; 3. contour field
loadct, 0 & contour, zefield, zex, zey, levels=lev, c_labels=findgen(n_elements(lev))*0.+1., /noerase, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, C_LINESTYLE = (lev LT 0.0)
;; 4. choose output
if (saveps eq 'true') then PS_End, /PNG
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
zefield       =  transpose(ztke)
;zefield       =  100.*transpose(ztke)/transpose(tke) & zefield[where(transpose(tke) lt 1.)]=0.
zex           =  localtime
zey           =  h
set_name      =  'zTKE.ps'
set_title     =  "Vertical Turbulent Kinetic Energy 0.5[<w'!U2!N>] (m!U2!N.s!U-2!N)"
set_titlex    =  'Local Time (h)'
set_titley    =  'Altitude above surface (km)'
set_subtitle  =  'Mean over the simulation domain'
set_xrange    =  [8.,18.]
set_yrange    =  [0.,10.]
set_tickx     =  1.
set_ticky     =  1.
minval        =  0.
maxval        =  14. ;10.
nlev          =  maxval-minval
pal           =  22
rrr           =  'no'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (saveps eq 'true') then PS_Start, FILENAME=set_name
!P.Charsize = 1.2
;; 0. levels
lev = minval + (maxval-minval)*findgen(nlev+1)/float(nlev) & if (minval ne 0.) then lev = lev[where(lev ne 0.)]
;;; 1. background
loadct, 0 & contour, /NODATA, zefield, zex, zey, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0
;; 2. color field
loadct, pal & if (rrr eq 'yes') then TVLCT, r, g, b, /Get & if (rrr eq 'yes') then TVLCT, Reverse(r), Reverse(g), Reverse(b)
            contour, zefield, zex, zey, levels=lev, c_labels=findgen(n_elements(lev))*0.+1., /overplot, /cell_fill
;; 3. contour field
loadct, 0 & contour, zefield, zex, zey, levels=lev, c_labels=findgen(n_elements(lev))*0.+1., /noerase, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, C_LINESTYLE = (lev LT 0.0)
;; 4. choose output
if (saveps eq 'true') then PS_End, /PNG
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
zefield       =  transpose(wt)
zex           =  localtime
zey           =  h
set_name      =  'HF.ps'
set_title     =  "Vertical Eddy Heat Flux <w'!7h!3'> (K.m.s!U-1!N)"
set_titlex    =  'Local Time (h)'
set_titley    =  'Altitude above surface (km)'
set_subtitle  =  'Mean over the simulation domain'
set_xrange    =  [8.,18.]
set_yrange    =  [0.,10.]
set_tickx     =  1.
set_ticky     =  1.
minval        =  -2.
maxval        =  2.
nlev          =  floor(maxval-minval)*10
pal           =  33 ;4 
rrr           =  'no'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (saveps eq 'true') then PS_Start, FILENAME=set_name
!P.Charsize = 1.2
;; 0. levels
lev = minval + (maxval-minval)*findgen(nlev+1)/float(nlev) & if (minval ne 0.) then lev = lev[where(lev ne 0.)]
;;; 1. background
loadct, 0 & contour, /NODATA, zefield, zex, zey, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0
;; 2. color field
loadct, pal & if (rrr eq 'yes') then TVLCT, r, g, b, /Get & if (rrr eq 'yes') then TVLCT, Reverse(r), Reverse(g), Reverse(b)
            ;;;--------------------------------------------------------------------------------------------------------------------------------
            ;;; WHITE ZONE - 1. get location of interval in the CT - 2. change the CT to have a white zone
            ulim=0.09 & dlim=-0.09 & w=where(lev le dlim) & n1=w[n_elements(w)-1] & w=where(lev ge ulim) & n2=w[0] & yy=BYTSCL(lev) & nd=yy[n1] & nu=yy[n2]-5
            nu = nd + (nu-nd)/2  ;; otherwise the interval is too large (because we removed 0)
            TVLCT, r, g, b, /Get & r[nd:nu]=255 & g[nd:nu]=255 & b[nd:nu]=255 & TVLCT, r, g, b
            ;;;--------------------------------------------------------------------------------------------------------------------------------
            contour, zefield, zex, zey, levels=lev, c_labels=findgen(n_elements(lev))*0.+1., /overplot, /cell_fill
;; 3. contour field
loadct, 0 & contour, zefield, zex, zey, levels=lev, c_labels=findgen(n_elements(lev))*0.+1., /noerase, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, C_LINESTYLE = (lev LT 0.0)
;; 4. choose output
if (saveps eq 'true') then PS_End, /PNG
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
zefield       =  t
zey           =  ht 
set_name      =  'T.ps'
set_title     =  "Potential Temperature (K)"
set_titlex    =  'Potential Temperature (K)'
set_titley    =  'Altitude above surface (km)'
set_subtitle  =  'Mean over the simulation domain'
set_xrange    =  [min(t),max(t)]
set_yrange    =  [0.,10.]
set_tickx     =  5.
set_ticky     =  1.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
localtimes = [9,10,11,12,13,14,15,16,17,18]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (saveps eq 'true') then PS_Start, FILENAME=set_name
!P.Charsize = 1.2
altlin=0 & loadct, 0
user_lt=localtimes[0] & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt)))) & nntt = yeah(0)
plot, zefield(*,nntt), zey(*,nntt), xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, linestyle=altlin
for ll = 1, n_elements(localtimes)-1 do begin
  CASE altlin OF
  0: altlin=1
  1: altlin=0
  ENDCASE  
  user_lt=localtimes[ll] & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt)))) & nntt = yeah(0)
  if (nntt ne -1) then oplot, zefield(*,nntt), zey(*,nntt), linestyle=altlin
endfor
if (saveps eq 'true') then PS_End, /PNG
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
zefield       =  stst
zey           =  ht 
set_name      =  'STST.ps'
set_title     =  'Static stability (K.m!U-1!N)'
set_titlex    =  'Static stability (K.m!U-1!N)'
set_titley    =  'Altitude above surface (km)'
set_subtitle  =  'Mean over the simulation domain'
set_xrange    =  [-1.,5.]
set_yrange    =  [0.,10.]
set_tickx     =  1.
set_ticky     =  1.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
localtimes = [9,11,13,15,17]
localtimes = [17]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (saveps eq 'true') then PS_Start, FILENAME=set_name
!P.Charsize = 1.2
altlin=0 & loadct, 0
user_lt=localtimes[0] & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt)))) & nntt = yeah(0)
plot, zefield(*,nntt), zey(*,nntt), xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, linestyle=altlin
oplot, zefield(*,nntt), zey(*,nntt), psym=5
for ll = 1, n_elements(localtimes)-1 do begin
  CASE altlin OF
  0: altlin=1 
  1: altlin=0
  ENDCASE
  user_lt=localtimes[ll] & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt)))) & nntt = yeah(0)
  if (nntt ne -1) then oplot, zefield(*,nntt), zey(*,nntt), linestyle=altlin
  ;if (nntt ne -1) then oplot, zefield(*,nntt), zey(*,nntt), psym=5
endfor
oplot, 0.*zefield(*,nntt) + 1.5, zey(*,nntt), linestyle=2  
if (saveps eq 'true') then PS_End, /PNG
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



stop




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
goto, no_staticstab
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

