!
! $Id: inigeom.F 1403 2010-07-01 09:02:53Z fairhead $
!
c
c
      SUBROUTINE inigeom
c
c     Auteur :  P. Le Van
c
c   ............      Version  du 01/04/2001     ........................
c
c  Calcul des elongations cuij1,.cuij4 , cvij1,..cvij4  aux memes en-
c     endroits que les aires aireij1,..aireij4 .

c  Choix entre f(y) a derivee sinusoid. ou a derivee tangente hyperbol.
c
c
      USE comconst_mod, ONLY: rad,g,omeg,pi
      USE logic_mod, ONLY: fxyhypb,ysinus
      USE serre_mod, ONLY: clon,clat,grossismx,grossismy,dzoomx,dzoomy,
     .		alphax,alphay,taux,tauy,transx,transy,pxo,pyo
      use fxhyp_m, only: fxhyp
      use fyhyp_m, only: fyhyp

      IMPLICIT NONE
c
#include "dimensions.h"
#include "paramet.h"
#include "comgeom2.h"
#include "comdissnew.h"

c-----------------------------------------------------------------------
c   ....  Variables  locales   ....
c
      INTEGER  i,j,itmax,itmay,iter
      REAL cvu(iip1,jjp1),cuv(iip1,jjm)
      REAL ai14,ai23,airez,rlatp,rlatm,xprm,xprp,un4rad2,yprp,yprm
      REAL eps,x1,xo1,f,df,xdm,y1,yo1,ydm
      REAL coslatm,coslatp,radclatm,radclatp
      REAL cuij1(iip1,jjp1),cuij2(iip1,jjp1),cuij3(iip1,jjp1),
     *     cuij4(iip1,jjp1)
      REAL cvij1(iip1,jjp1),cvij2(iip1,jjp1),cvij3(iip1,jjp1),
     *     cvij4(iip1,jjp1)
      REAL rlonvv(iip1),rlatuu(jjp1)
      REAL rlatu1(jjm),yprimu1(jjm),rlatu2(jjm),yprimu2(jjm) ,
     *     yprimv(jjm),yprimu(jjp1)
      REAL gamdi_gdiv, gamdi_grot, gamdi_h
 
      REAL rlonm025(iip1),xprimm025(iip1), rlonp025(iip1),
     ,  xprimp025(iip1)
      SAVE rlatu1,yprimu1,rlatu2,yprimu2,yprimv,yprimu
      SAVE rlonm025,xprimm025,rlonp025,xprimp025

      REAL      SSUM
c
c
c   ------------------------------------------------------------------
c   -                                                                -
c   -    calcul des coeff. ( cu, cv , 1./cu**2,  1./cv**2  )         -
c   -                                                                -
c   ------------------------------------------------------------------
c
c      les coef. ( cu, cv ) permettent de passer des vitesses naturelles
c      aux vitesses covariantes et contravariantes , ou vice-versa ...
c
c
c     on a :  u (covariant) = cu * u (naturel)   , u(contrav)= u(nat)/cu
c             v (covariant) = cv * v (naturel)   , v(contrav)= v(nat)/cv
c
c       on en tire :  u(covariant) = cu * cu * u(contravariant)
c                     v(covariant) = cv * cv * v(contravariant)
c
c
c     on a l'application (  x(X) , y(Y) )   avec - im/2 +1 <  X  < im/2
c                                                          =     =
c                                           et   - jm/2    <  Y  < jm/2
c                                                          =     =
c
c      ...................................................
c      ...................................................
c      .  x  est la longitude du point  en radians       .
c      .  y  est la  latitude du point  en radians       .
c      .                                                 .
c      .  on a :  cu(i,j) = rad * COS(y) * dx/dX         .
c      .          cv( j ) = rad          * dy/dY         .
c      .        aire(i,j) =  cu(i,j) * cv(j)             .
c      .                                                 .
c      . y, dx/dX, dy/dY calcules aux points concernes   .
c      .                                                 .
c      ...................................................
c      ...................................................
c
c
c
c                                                           ,
c    cv , bien que dependant de j uniquement,sera ici indice aussi en i
c          pour un adressage plus facile en  ij  .
c
c
c
c  **************  aux points  u  et  v ,           *****************
c      xprimu et xprimv sont respectivement les valeurs de  dx/dX
c      yprimu et yprimv    .  .  .  .  .  .  .  .  .  .  .  dy/dY
c      rlatu  et  rlatv    .  .  .  .  .  .  .  .  .  .  .la latitude
c      cvu    et   cv      .  .  .  .  .  .  .  .  .  .  .    cv
c
c  **************  aux points u, v, scalaires, et z  ****************
c      cu, cuv, cuscal, cuz sont respectiv. les valeurs de    cu
c
c
c
c         Exemple de distribution de variables sur la grille dans le
c             domaine de travail ( X,Y ) .
c     ................................................................
c                  DX=DY= 1
c
c   
c        +     represente  un  point scalaire ( p.exp  la pression )
c        >     represente  la composante zonale du  vent
c        V     represente  la composante meridienne du vent
c        o     represente  la  vorticite
c
c       ----  , car aux poles , les comp.zonales covariantes sont nulles
c
c
c
c         i ->
c
c         1      2      3      4      5      6      7      8
c  j
c  v  1   + ---- + ---- + ---- + ---- + ---- + ---- + ---- + --
c
c         V   o  V   o  V   o  V   o  V   o  V   o  V   o  V  o
c
c     2   +   >  +   >  +   >  +   >  +   >  +   >  +   >  +  >
c
c         V   o  V   o  V   o  V   o  V   o  V   o  V   o  V  o
c
c     3   +   >  +   >  +   >  +   >  +   >  +   >  +   >  +  >
c
c         V   o  V   o  V   o  V   o  V   o  V   o  V   o  V  o
c
c     4   +   >  +   >  +   >  +   >  +   >  +   >  +   >  +  >
c
c         V   o  V   o  V   o  V   o  V   o  V   o  V   o  V  o
c
c     5   + ---- + ---- + ---- + ---- + ---- + ---- + ---- + --
c
c
c      Ci-dessus,  on voit que le nombre de pts.en longitude est egal
c                 a   IM = 8
c      De meme ,   le nombre d'intervalles entre les 2 poles est egal
c                 a   JM = 4
c
c      Les points scalaires ( + ) correspondent donc a des valeurs
c       entieres  de  i ( 1 a IM )   et  de  j ( 1 a  JM +1 )   .
c
c      Les vents    U       ( > ) correspondent a des valeurs  semi-
c       entieres  de i ( 1+ 0.5 a IM+ 0.5) et entieres de j ( 1 a JM+1)
c
c      Les vents    V       ( V ) correspondent a des valeurs entieres
c       de     i ( 1 a  IM ) et semi-entieres de  j ( 1 +0.5  a JM +0.5)
c
c
c
      WRITE(6,3) 
 3    FORMAT( // 10x,' ....  INIGEOM  date du 01/06/98   ..... ',
     * //5x,'   Calcul des elongations cu et cv  comme sommes des 4 ' /
     *  5x,' elong. cuij1, .. 4  , cvij1,.. 4  qui les entourent , aux 
     * '/ 5x,' memes endroits que les aires aireij1,...j4   . ' / )
c
c
      IF( nitergdiv.NE.2 ) THEN
        gamdi_gdiv = coefdis/ ( REAL(nitergdiv) -2. )
      ELSE
        gamdi_gdiv = 0.
      ENDIF
      IF( nitergrot.NE.2 ) THEN
        gamdi_grot = coefdis/ ( REAL(nitergrot) -2. )
      ELSE
        gamdi_grot = 0.
      ENDIF
      IF( niterh.NE.2 ) THEN
        gamdi_h = coefdis/ ( REAL(niterh) -2. )
      ELSE
        gamdi_h = 0.
      ENDIF

      WRITE(6,*) ' gamdi_gd ',gamdi_gdiv,gamdi_grot,gamdi_h,coefdis,
     *  nitergdiv,nitergrot,niterh
c
      pi    = 2.* ASIN(1.)
c
      WRITE(6,990) 

c     ----------------------------------------------------------------
c
      IF( .NOT.fxyhypb )   THEN
c
c
       IF( ysinus )  THEN
c
        WRITE(6,*) ' ***  Inigeom ,  Y = Sinus ( Latitude ) *** '
c
c   .... utilisation de f(x,y )  avec  y  =  sinus de la latitude  .....

        CALL  fxysinus (rlatu,yprimu,rlatv,yprimv,rlatu1,yprimu1,
     ,                    rlatu2,yprimu2,
     ,  rlonu,xprimu,rlonv,xprimv,rlonm025,xprimm025,rlonp025,xprimp025)

       ELSE
c
        WRITE(6,*) '*** Inigeom ,  Y = Latitude  , der. sinusoid . ***'

c  .... utilisation  de f(x,y) a tangente sinusoidale , y etant la latit. ...
c
 
        pxo   = clon *pi /180.
        pyo   = 2.* clat* pi /180.
c
c  ....  determination de  transx ( pour le zoom ) par Newton-Raphson ...
c
        itmax = 10
        eps   = .1e-7
c
        xo1 = 0.
        DO 10 iter = 1, itmax
        x1  = xo1
        f   = x1+ alphax *SIN(x1-pxo)
        df  = 1.+ alphax *COS(x1-pxo)
        x1  = x1 - f/df
        xdm = ABS( x1- xo1 )
        IF( xdm.LE.eps )GO TO 11
        xo1 = x1
 10     CONTINUE
 11     CONTINUE
c
        transx = xo1

        itmay = 10
        eps   = .1e-7
C
        yo1  = 0.
        DO 15 iter = 1,itmay
        y1   = yo1
        f    = y1 + alphay* SIN(y1-pyo)
        df   = 1. + alphay* COS(y1-pyo)
        y1   = y1 -f/df
        ydm  = ABS(y1-yo1)
        IF(ydm.LE.eps) GO TO 17
        yo1  = y1
 15     CONTINUE
c
 17     CONTINUE
        transy = yo1

        CALL fxy ( rlatu,yprimu,rlatv,yprimv,rlatu1,yprimu1,
     ,              rlatu2,yprimu2,
     ,  rlonu,xprimu,rlonv,xprimv,rlonm025,xprimm025,rlonp025,xprimp025)

       ENDIF
c
      ELSE
c
c   ....  Utilisation  de fxyhyper , f(x,y) a derivee tangente hyperbol.
c   .....................................................................

      WRITE(6,*)'*** Inigeom , Y = Latitude  , der.tg. hyperbolique ***'
 
      CALL fyhyp(rlatu, yprimu, rlatv, rlatu2, yprimu2, rlatu1, yprimu1)
      CALL fxhyp(xprimm025, rlonv, xprimv, rlonu, xprimu, xprimp025)
  
      ENDIF
c
c  -------------------------------------------------------------------

c
      rlatu(1)    =     ASIN(1.)
      rlatu(jjp1) =  - rlatu(1)
c
c
c   ....  calcul  aux  poles  ....
c
      yprimu(1)      = 0.
      yprimu(jjp1)   = 0.
c
c
      un4rad2 = 0.25 * rad * rad
c
c   --------------------------------------------------------------------
c   --------------------------------------------------------------------
c   -                                                                  -
c   -  calcul  des aires ( aire,aireu,airev, 1./aire, 1./airez  )      -
c   -      et de   fext ,  force de coriolis  extensive  .             -
c   -                                                                  -
c   --------------------------------------------------------------------
c   --------------------------------------------------------------------
c
c
c
c   A 1 point scalaire P (i,j) de la grille, reguliere en (X,Y) , sont
c   affectees 4 aires entourant P , calculees respectivement aux points
c            ( i + 1/4, j - 1/4 )    :    aireij1 (i,j)
c            ( i + 1/4, j + 1/4 )    :    aireij2 (i,j)
c            ( i - 1/4, j + 1/4 )    :    aireij3 (i,j)
c            ( i - 1/4, j - 1/4 )    :    aireij4 (i,j)
c
c           ,
c   Les cotes de chacun de ces 4 carres etant egaux a 1/2 suivant (X,Y).
c   Chaque aire centree en 1 point scalaire P(i,j) est egale a la somme
c   des 4 aires  aireij1,aireij2,aireij3,aireij4 qui sont affectees au
c   point (i,j) .
c   On definit en outre les coefficients  alpha comme etant egaux a
c    (aireij / aire), c.a.d par exp.  alpha1(i,j)=aireij1(i,j)/aire(i,j)
c
c   De meme, toute aire centree en 1 point U est egale a la somme des
c   4 aires aireij1,aireij2,aireij3,aireij4 entourant le point U .
c    Idem pour  airev, airez .
c
c       On a ,pour chaque maille :    dX = dY = 1
c
c
c                             . V
c
c                 aireij4 .        . aireij1
c
c                   U .       . P      . U
c
c                 aireij3 .        . aireij2
c
c                             . V
c
c
c
c
c
c   ....................................................................
c
c    Calcul des 4 aires elementaires aireij1,aireij2,aireij3,aireij4
c   qui entourent chaque aire(i,j) , ainsi que les 4 elongations elemen
c   taires cuij et les 4 elongat. cvij qui sont calculees aux memes 
c     endroits  que les aireij   .    
c
c   ....................................................................
c
c     .......  do 35  :   boucle sur les  jjm + 1  latitudes   .....
c
c
      DO 35 j = 1, jjp1
c
      IF ( j. eq. 1 )  THEN
c
      yprm           = yprimu1(j)
      rlatm          = rlatu1(j)
c
      coslatm        = COS( rlatm )
      radclatm       = 0.5* rad * coslatm
c
      DO 30 i = 1, iim
      xprp           = xprimp025( i )
      xprm           = xprimm025( i )
      aireij2( i,1 ) = un4rad2 * coslatm  * xprp * yprm
      aireij3( i,1 ) = un4rad2 * coslatm  * xprm * yprm
      cuij2  ( i,1 ) = radclatm * xprp
      cuij3  ( i,1 ) = radclatm * xprm
      cvij2  ( i,1 ) = 0.5* rad * yprm
      cvij3  ( i,1 ) = cvij2(i,1)
  30  CONTINUE
c
      DO  i = 1, iim
      aireij1( i,1 ) = 0.
      aireij4( i,1 ) = 0.
      cuij1  ( i,1 ) = 0.
      cuij4  ( i,1 ) = 0.
      cvij1  ( i,1 ) = 0.
      cvij4  ( i,1 ) = 0.
      ENDDO
c
      END IF
c
      IF ( j. eq. jjp1 )  THEN
       yprp               = yprimu2(j-1)
       rlatp              = rlatu2 (j-1)
ccc       yprp             = fyprim( REAL(j) - 0.25 )
ccc       rlatp            = fy    ( REAL(j) - 0.25 )
c
      coslatp             = COS( rlatp )
      radclatp            = 0.5* rad * coslatp
c
      DO 31 i = 1,iim
        xprp              = xprimp025( i )
        xprm              = xprimm025( i )
        aireij1( i,jjp1 ) = un4rad2 * coslatp  * xprp * yprp
        aireij4( i,jjp1 ) = un4rad2 * coslatp  * xprm * yprp
        cuij1(i,jjp1)     = radclatp * xprp
        cuij4(i,jjp1)     = radclatp * xprm
        cvij1(i,jjp1)     = 0.5 * rad* yprp
        cvij4(i,jjp1)     = cvij1(i,jjp1)
 31   CONTINUE
c
       DO   i    = 1, iim
        aireij2( i,jjp1 ) = 0.
        aireij3( i,jjp1 ) = 0.
        cvij2  ( i,jjp1 ) = 0.
        cvij3  ( i,jjp1 ) = 0.
        cuij2  ( i,jjp1 ) = 0.
        cuij3  ( i,jjp1 ) = 0.
       ENDDO
c
      END IF
c

      IF ( j .gt. 1 .AND. j .lt. jjp1 )  THEN
c
        rlatp    = rlatu2 ( j-1 )
        yprp     = yprimu2( j-1 )
        rlatm    = rlatu1 (  j  )
        yprm     = yprimu1(  j  )
cc         rlatp    = fy    ( REAL(j) - 0.25 )
cc         yprp     = fyprim( REAL(j) - 0.25 )
cc         rlatm    = fy    ( REAL(j) + 0.25 )
cc         yprm     = fyprim( REAL(j) + 0.25 )

         coslatm  = COS( rlatm )
         coslatp  = COS( rlatp )
         radclatp = 0.5* rad * coslatp
         radclatm = 0.5* rad * coslatm
c
         ai14            = un4rad2 * coslatp * yprp
         ai23            = un4rad2 * coslatm * yprm
         DO 32 i = 1,iim
         xprp            = xprimp025( i )
         xprm            = xprimm025( i )
      
         aireij1 ( i,j ) = ai14 * xprp
         aireij2 ( i,j ) = ai23 * xprp
         aireij3 ( i,j ) = ai23 * xprm
         aireij4 ( i,j ) = ai14 * xprm
         cuij1   ( i,j ) = radclatp * xprp
         cuij2   ( i,j ) = radclatm * xprp
         cuij3   ( i,j ) = radclatm * xprm
         cuij4   ( i,j ) = radclatp * xprm
         cvij1   ( i,j ) = 0.5* rad * yprp
         cvij2   ( i,j ) = 0.5* rad * yprm
         cvij3   ( i,j ) = cvij2(i,j)
         cvij4   ( i,j ) = cvij1(i,j)
  32     CONTINUE
c
      END IF
c
c    ........       periodicite   ............
c
         cvij1   (iip1,j) = cvij1   (1,j)
         cvij2   (iip1,j) = cvij2   (1,j)
         cvij3   (iip1,j) = cvij3   (1,j)
         cvij4   (iip1,j) = cvij4   (1,j)
         cuij1   (iip1,j) = cuij1   (1,j)
         cuij2   (iip1,j) = cuij2   (1,j)
         cuij3   (iip1,j) = cuij3   (1,j)
         cuij4   (iip1,j) = cuij4   (1,j)
         aireij1 (iip1,j) = aireij1 (1,j )
         aireij2 (iip1,j) = aireij2 (1,j )
         aireij3 (iip1,j) = aireij3 (1,j )
         aireij4 (iip1,j) = aireij4 (1,j )
        
  35  CONTINUE
c
c    ..............................................................
c
      DO 37 j = 1, jjp1
      DO 36 i = 1, iim
      aire    ( i,j )  = aireij1(i,j) + aireij2(i,j) + aireij3(i,j) +
     *                          aireij4(i,j)
      alpha1  ( i,j )  = aireij1(i,j) / aire(i,j)
      alpha2  ( i,j )  = aireij2(i,j) / aire(i,j)
      alpha3  ( i,j )  = aireij3(i,j) / aire(i,j)
      alpha4  ( i,j )  = aireij4(i,j) / aire(i,j)
      alpha1p2( i,j )  = alpha1 (i,j) + alpha2 (i,j)
      alpha1p4( i,j )  = alpha1 (i,j) + alpha4 (i,j)
      alpha2p3( i,j )  = alpha2 (i,j) + alpha3 (i,j)
      alpha3p4( i,j )  = alpha3 (i,j) + alpha4 (i,j)
  36  CONTINUE
c
c
      aire    (iip1,j) = aire    (1,j)
      alpha1  (iip1,j) = alpha1  (1,j)
      alpha2  (iip1,j) = alpha2  (1,j)
      alpha3  (iip1,j) = alpha3  (1,j)
      alpha4  (iip1,j) = alpha4  (1,j)
      alpha1p2(iip1,j) = alpha1p2(1,j)
      alpha1p4(iip1,j) = alpha1p4(1,j)
      alpha2p3(iip1,j) = alpha2p3(1,j)
      alpha3p4(iip1,j) = alpha3p4(1,j)
  37  CONTINUE
c

      DO 42 j = 1,jjp1
      DO 41 i = 1,iim
      aireu       (i,j)= aireij1(i,j) + aireij2(i,j) + aireij4(i+1,j) +
     *                                aireij3(i+1,j)
      unsaire    ( i,j)= 1./ aire(i,j)
      unsair_gam1( i,j)= unsaire(i,j)** ( - gamdi_gdiv )
      unsair_gam2( i,j)= unsaire(i,j)** ( - gamdi_h    )
      airesurg   ( i,j)= aire(i,j)/ g
  41  CONTINUE
      aireu     (iip1,j)  = aireu  (1,j)
      unsaire   (iip1,j)  = unsaire(1,j)
      unsair_gam1(iip1,j) = unsair_gam1(1,j)
      unsair_gam2(iip1,j) = unsair_gam2(1,j)
      airesurg   (iip1,j) = airesurg(1,j)
  42  CONTINUE
c
c
      DO 48 j = 1,jjm
c
        DO i=1,iim
         airev     (i,j) = aireij2(i,j)+ aireij3(i,j)+ aireij1(i,j+1) +
     *                           aireij4(i,j+1)
        ENDDO
         DO i=1,iim
          airez         = aireij2(i,j)+aireij1(i,j+1)+aireij3(i+1,j) +
     *                           aireij4(i+1,j+1)
          unsairez(i,j) = 1./ airez
          unsairz_gam(i,j)= unsairez(i,j)** ( - gamdi_grot )
          fext    (i,j)   = airez * SIN(rlatv(j))* 2.* omeg
         ENDDO
        airev     (iip1,j)  = airev(1,j)
        unsairez  (iip1,j)  = unsairez(1,j)
        fext      (iip1,j)  = fext(1,j)
        unsairz_gam(iip1,j) = unsairz_gam(1,j)
c
  48  CONTINUE
c
c
c    .....      Calcul  des elongations cu,cv, cvu     .........
c
      DO    j   = 1, jjm
       DO   i  = 1, iim
       cv(i,j) = 0.5 *( cvij2(i,j)+cvij3(i,j)+cvij1(i,j+1)+cvij4(i,j+1))
       cvu(i,j)= 0.5 *( cvij1(i,j)+cvij4(i,j)+cvij2(i,j)  +cvij3(i,j) )
       cuv(i,j)= 0.5 *( cuij2(i,j)+cuij3(i,j)+cuij1(i,j+1)+cuij4(i,j+1))
       unscv2(i,j) = 1./ ( cv(i,j)*cv(i,j) )
       ENDDO
       DO   i  = 1, iim
       cuvsurcv (i,j)    = airev(i,j)  * unscv2(i,j)
       cvsurcuv (i,j)    = 1./cuvsurcv(i,j)
       cuvscvgam1(i,j)   = cuvsurcv (i,j) ** ( - gamdi_gdiv )
       cuvscvgam2(i,j)   = cuvsurcv (i,j) ** ( - gamdi_h )
       cvscuvgam(i,j)    = cvsurcuv (i,j) ** ( - gamdi_grot )
       ENDDO
       cv       (iip1,j)  = cv       (1,j)
       cvu      (iip1,j)  = cvu      (1,j)
       unscv2   (iip1,j)  = unscv2   (1,j)
       cuv      (iip1,j)  = cuv      (1,j)
       cuvsurcv (iip1,j)  = cuvsurcv (1,j)
       cvsurcuv (iip1,j)  = cvsurcuv (1,j)
       cuvscvgam1(iip1,j) = cuvscvgam1(1,j)
       cuvscvgam2(iip1,j) = cuvscvgam2(1,j)
       cvscuvgam(iip1,j)  = cvscuvgam(1,j)
      ENDDO

      DO  j     = 2, jjm
        DO   i  = 1, iim
        cu(i,j) = 0.5*(cuij1(i,j)+cuij4(i+1,j)+cuij2(i,j)+cuij3(i+1,j))
        unscu2    (i,j)  = 1./ ( cu(i,j) * cu(i,j) )
        cvusurcu  (i,j)  =  aireu(i,j) * unscu2(i,j)
        cusurcvu  (i,j)  = 1./ cvusurcu(i,j)
        cvuscugam1 (i,j) = cvusurcu(i,j) ** ( - gamdi_gdiv ) 
        cvuscugam2 (i,j) = cvusurcu(i,j) ** ( - gamdi_h    ) 
        cuscvugam (i,j)  = cusurcvu(i,j) ** ( - gamdi_grot )
        ENDDO
        cu       (iip1,j)  = cu(1,j)
        unscu2   (iip1,j)  = unscu2(1,j)
        cvusurcu (iip1,j)  = cvusurcu(1,j)
        cusurcvu (iip1,j)  = cusurcvu(1,j)
        cvuscugam1(iip1,j) = cvuscugam1(1,j)
        cvuscugam2(iip1,j) = cvuscugam2(1,j)
        cuscvugam (iip1,j) = cuscvugam(1,j)
      ENDDO

c
c   ....  calcul aux  poles  ....
c
      DO    i      =  1, iip1
        cu    ( i, 1 )  =   0.
        unscu2( i, 1 )  =   0.
        cvu   ( i, 1 )  =   0.
c
        cu    (i, jjp1) =   0.
        unscu2(i, jjp1) =   0.
        cvu   (i, jjp1) =   0.
      ENDDO
c
c    ..............................................................
c
      DO j = 1, jjm
        DO i= 1, iim
         airvscu2  (i,j) = airev(i,j)/ ( cuv(i,j) * cuv(i,j) )
         aivscu2gam(i,j) = airvscu2(i,j)** ( - gamdi_grot )
        ENDDO
         airvscu2  (iip1,j)  = airvscu2(1,j)
         aivscu2gam(iip1,j)  = aivscu2gam(1,j)
      ENDDO

      DO j=2,jjm
        DO i=1,iim
         airuscv2   (i,j)    = aireu(i,j)/ ( cvu(i,j) * cvu(i,j) )
         aiuscv2gam (i,j)    = airuscv2(i,j)** ( - gamdi_grot ) 
        ENDDO
         airuscv2  (iip1,j)  = airuscv2  (1,j)
         aiuscv2gam(iip1,j)  = aiuscv2gam(1,j)
      ENDDO

c
c   calcul des aires aux  poles :
c   -----------------------------
c
      apoln       = SSUM(iim,aire(1,1),1)
      apols       = SSUM(iim,aire(1,jjp1),1)
      unsapolnga1 = 1./ ( apoln ** ( - gamdi_gdiv ) )
      unsapolsga1 = 1./ ( apols ** ( - gamdi_gdiv ) )
      unsapolnga2 = 1./ ( apoln ** ( - gamdi_h    ) )
      unsapolsga2 = 1./ ( apols ** ( - gamdi_h    ) )
c
c-----------------------------------------------------------------------
c     gtitre='Coriolis version ancienne'
c     gfichier='fext1'
c     CALL writestd(fext,iip1*jjm)
c
c   changement F. Hourdin calcul conservatif pour fext
c   constang contient le produit a * cos ( latitude ) * omega
c
      DO i=1,iim
         constang(i,1) = 0.
      ENDDO
      DO j=1,jjm-1
        DO i=1,iim
         constang(i,j+1) = rad*omeg*cu(i,j+1)*COS(rlatu(j+1))
        ENDDO
      ENDDO
      DO i=1,iim
         constang(i,jjp1) = 0.
      ENDDO
c
c   periodicite en longitude
c
      DO j=1,jjm
        fext(iip1,j)     = fext(1,j)
      ENDDO
      DO j=1,jjp1
        constang(iip1,j) = constang(1,j)
      ENDDO

c fin du changement

c
c-----------------------------------------------------------------------
c
       WRITE(6,*) '   ***  Coordonnees de la grille  *** '
       WRITE(6,995)
c
       WRITE(6,*) '   LONGITUDES  aux pts.   V  ( degres )  '
       WRITE(6,995)
        DO i=1,iip1
         rlonvv(i) = rlonv(i)*180./pi
        ENDDO
       WRITE(6,400) rlonvv
c
       WRITE(6,995)
       WRITE(6,*) '   LATITUDES   aux pts.   V  ( degres )  '
       WRITE(6,995)
        DO i=1,jjm
         rlatuu(i)=rlatv(i)*180./pi
        ENDDO
       WRITE(6,400) (rlatuu(i),i=1,jjm)
c
        DO i=1,iip1
          rlonvv(i)=rlonu(i)*180./pi
        ENDDO
       WRITE(6,995)
       WRITE(6,*) '   LONGITUDES  aux pts.   U  ( degres )  '
       WRITE(6,995)
       WRITE(6,400) rlonvv
       WRITE(6,995)

       WRITE(6,*) '   LATITUDES   aux pts.   U  ( degres )  '
       WRITE(6,995)
        DO i=1,jjp1
         rlatuu(i)=rlatu(i)*180./pi
        ENDDO
       WRITE(6,400) (rlatuu(i),i=1,jjp1)
       WRITE(6,995)
c
444    format(f10.3,f6.0)
400    FORMAT(1x,8f8.2)
990    FORMAT(//)
995    FORMAT(/)
c
      RETURN
      END
