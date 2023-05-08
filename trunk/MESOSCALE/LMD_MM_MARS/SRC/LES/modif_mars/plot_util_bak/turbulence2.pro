pro turbulence2


retrieve='true'

fast='true'
;fast='false'


p0=610. & t0=220. & r_cp=1/4.4 & grav=3.72 & R=192. & g_cp=1000.*grav/844.6
;;;;;;;;;;;;;;;;;;;;;;;
;                     ;
; HINSON - HINSON !!! ;
;                     ;
r_cp = 0.25           ; 
;g_cp = 4.9           ;
;;;;;;;;;;;;;;;;;;;;;;;



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



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
; get fields and dimensions
;
getcdf, file='t.nc', charvar='T', invar=t
getcdf, file='p.nc', charvar='PTOT', invar=p
getcdf, file='ph.nc', charvar='PHTOT', invar=ph
nx=n_elements(t(*,0,0,0))
ny=n_elements(t(0,*,0,0))
nz=n_elements(t(0,0,*,0))
nt=n_elements(t(0,0,0,*))
   ;;;; pour eviter de detecter des fluctuations plus basses que PBL top ;;;
     smooth_gridpoint = 2 & smooth_time = 2 & smooth_vert = 2  
     ;smooth_gridpoint = 0 & smooth_time = 0 & smooth_vert = 0  
     t = smooth(t,[smooth_gridpoint,smooth_gridpoint,smooth_vert,smooth_time],/EDGE_TRUNCATE)
     p = smooth(p,[smooth_gridpoint,smooth_gridpoint,smooth_vert,smooth_time],/EDGE_TRUNCATE)
     ;ph = smooth(ph,[smooth_gridpoint,smooth_gridpoint,smooth_vert,smooth_time],/EDGE_TRUNCATE)
   ;;;; pour eviter de detecter des fluctuations plus basses que PBL top ;;;

;
; localtime loop
;
localtime = lctu + 100.*findgen(nt)/3700.
  ;
  ; radio-occultations
  ;
  user_lt=17.0 & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt))))
  user_lt=17.5 & yeah2=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt))))
  ;;
  ;; PBL growth
  ;;
  ;user_lt=14.5 & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt))))
  ;user_lt=18.0 & yeah2=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt))))
  ;

  user_lt=14.0 & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt))))
  user_lt=14.5 & yeah2=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt))))



deb=yeah(0)
fin=yeah2(0)

;
; calculate over domain
;
pbl_depth_tab=findgen(nx,ny,fin-deb+1)
nloop = nx*ny

;
; LOOOOOOP
;
for tt=deb,fin do begin

 tprof_mean = 0.
 heightp_mean = 0.
 staticstab_mean = 0.
 nbad = 0

for i=0,nx-1 do begin
for j=0,ny-1 do begin

 height=reform( ph( i,j,*,tt ) ) / 1000. / grav
 heightp=height(0:nz-1) ;no stagger

 thprof = t0 + reform( t( i,j,*,tt ) )
 pprof = reform( p( i,j,*,tt ) )
 tprof = thprof * ( pprof/p0 )^r_cp

 staticstab = DERIV( heightp , tprof ) + g_cp  

if (fast eq 'true') then begin

; MOINS BON mais PLUS RAPIDE
tprof_mean = tprof_mean + tprof/nloop
heightp_mean = heightp_mean + heightp/nloop
staticstab_mean = staticstab_mean + staticstab/nloop

endif else begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  z_in = heightp
  profile_in = staticstab
  resol = 500.  ;200 un peu juste
  ;********************************************************************
  ; ---- numerical recipes in C - section 3.3 ----
  ;
  ; Calculate interpolating cubic spline 
        yspline = SPL_INIT(z_in,profile_in)
  ; Define the X values P at which we desire interpolated Y values
        xspline=min(z_in)+findgen(floor(max(z_in)-min(z_in))*resol)/resol
  ; Calculate the interpolated Y values
        result=spl_interp(z_in,profile_in,yspline,xspline)
  ;********************************************************************
  staticstab = result
  heightp = xspline

  heightp_mean = heightp_mean + heightp/nloop
  staticstab_mean = staticstab_mean + staticstab/nloop
    yspline = SPL_INIT(z_in,tprof)
    tprof = spl_interp(z_in,tprof,yspline,xspline)
  tprof_mean = tprof_mean + tprof/nloop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

endelse

lim=3.
lim=1.5
w = where( (staticstab ge lim) and (heightp-heightp(0) ge lim) )
pbl_depth_tab(i,j,tt-deb) = heightp(w(0)) - heightp(0)

   ; garde-fou supplementaire
   if (heightp(w(0)) - heightp(0) lt 2.) then begin
     pbl_depth_tab(i,j,tt-deb) = !VALUES.F_NAN
     ;plot, profile_in, z_in, xrange=[-1.,5.], psym=5, title='bad!' & oplot, staticstab, heightp & oplot, staticstab, findgen(nz)*0. + heightp(w(0)), linestyle=2
     print, heightp(w) - heightp(0)
     nbad = nbad+1
   endif

;plot, profile_in, z_in, xrange=[-1.,5.], psym=5 & oplot, staticstab, heightp & oplot, staticstab, findgen(nz)*0. + heightp(w(0)), linestyle=2 
;pause

endfor
;print, i
endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  print, '-------------------------------------------------------'
  print, 'Local Time : ', localtime(tt)
  print, 'nbad ', nbad
  print, 'PBL height - min ', min(pbl_depth_tab(*,*,tt-deb))
  print, 'PBL height - max ', max(pbl_depth_tab(*,*,tt-deb))
  print, 'PBL height - mean ', mean(pbl_depth_tab(*,*,tt-deb))
    ;; static stability of mean profile
    staticstab = DERIV( heightp_mean , tprof_mean ) + g_cp

	  ;z_in = heightp_mean
	  ;profile_in = staticstab
	  ;resol = 500  ;200 un peu juste
	  ;;********************************************************************
	  ;; ---- numerical recipes in C - section 3.3 ----
	  ;;
	  ;; Calculate interpolating cubic spline 
	  ;      yspline = SPL_INIT(z_in,profile_in)
	  ;; Define the X values P at which we desire interpolated Y values
	  ;      xspline=min(z_in)+findgen(floor(max(z_in))*resol)/resol
	  ;; Calculate the interpolated Y values
	  ;      result=spl_interp(z_in,profile_in,yspline,xspline)
	  ;;********************************************************************
	  ;staticstab = result
	  ;w = where((staticstab ge 1.5) and (xspline-xspline(0) gt 2.))
	  ;print, 'PBL height - mean profile', xspline(w(0)) - xspline(0)

    w = where((staticstab ge 1.5) and (heightp_mean-heightp_mean(0) gt 1.5))
    print, 'PBL height - mean profile', heightp_mean(w(0))-heightp_mean(0)
    ;; mean static stability
    w = where((staticstab_mean ge 1.5) and (heightp_mean-heightp_mean(0) gt 1.5))
    print, 'PBL height - mean static stability', heightp_mean(w(0)) - heightp_mean(0)

   ;plot, staticstab, heightp_mean, xrange=[-1.,5.]
   ;oplot, staticstab, findgen(nz)*0. + heightp_mean(w(0)), linestyle=2
   ;oplot, staticstab_mean, heightp_mean
   ;oplot, staticstab_mean, findgen(nz)*0. + heightp_mean(w(0)), linestyle=2
   ;pause
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
endfor
endif

  print, '-------------------------------------------------------'
  print, 'PBL height - min ', min(pbl_depth_tab(*,*,*))
  print, 'PBL height - max ', max(pbl_depth_tab(*,*,*))
  print, 'PBL height - mean ', mean(pbl_depth_tab(*,*,*))
  print, '-------------------------------------------------------'

stop

!p.charthick = 2.0
!p.thick = 3.0
!x.thick = 2.0
!y.thick = 2.0
set_plot, 'ps' & device, filename='plot/map.ps', /color
;map_latlon, reform(pbl_depth_tab(*,*,0))


what_I_plot = reform(pbl_depth_tab(*,*,0))
pal=4
colors=128
title_user='PBL depth'
map_latlon, $
        what_I_plot, $                          ; 2D field
        lon, $                                  ; 1D latitude
        lat, $                                  ; 1D longitude
;        minfield=minfield_init, $               ; minimum value of plotted field (=0: calculate)
;        maxfield=maxfield_init, $               ; maximum value of plotted field (=0: calculate)
;        overcontour=overcontour, $              ; another 2D field to overplot with contour lines (=0: no)
;        overvector_x=overvector_x, $            ; wind vector - x component (=0: no)
;        overvector_y=overvector_y, $            ; wind vector - y component (=0: no)
        ct=pal, $                               ; color table (33-rainbow is default)
        colors=colors, $                        ; number of colors/levels (32 is default)
        title=title_user;, $                     ; title of the plot ('' is default)
;        format=format                           ; format of colorbar annotations ('(F6.2)' is default)



device, /close
end
