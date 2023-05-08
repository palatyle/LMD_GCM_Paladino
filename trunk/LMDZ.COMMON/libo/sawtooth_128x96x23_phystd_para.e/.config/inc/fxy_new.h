!
! $Header$
!
c--------------------------------------------------------------
         REAL ripx
         REAL fx,fxprim,fy,fyprim,ri,rj,bigy
c
c....stretching in x...
c
        ripx(  ri )= (ri-1.0) *2.*pi/REAL(iim) 
        fx  (  ri )= ripx(ri) + transx  +
     *         alphax * SIN( ripx(ri)+transx-pxo ) - pi
        fxprim(ri) = 2.*pi/REAL(iim)  *
     *        ( 1.+ alphax * COS( ripx(ri)+transx-pxo ) )

c....stretching in y...
c
        bigy(rj)   = 2.* (REAL(jjp1)-rj ) *pi/jjm
        fy(rj)     =  ( bigy(rj) + transy  +
     *        alphay * SIN( bigy(rj)+transy-pyo ) ) /2.  - pi/2.
        fyprim(rj) = ( pi/jjm ) * ( 1.+
     *           alphay * COS( bigy(rj)+transy-pyo ) )

c       fy(rj)= pyo-pisjjm*(rj-transy)+coefalpha*SIN(depisjm*(rj-
c     *  transy ))
c       fyprim(rj)= pisjjm-pisjjm*coefy2* COS(depisjm*(rj-transy)) 
c--------------------------------------------------------------
