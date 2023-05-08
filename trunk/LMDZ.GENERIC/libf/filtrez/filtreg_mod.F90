MODULE filtreg_mod

CONTAINS

      SUBROUTINE inifilr
!
!    ... H. Upadhyaya, O.Sharma   ...
!
      USE logic_mod, ONLY: fxyhypb,ysinus
      USE serre_mod, ONLY: alphax

      IMPLICIT NONE
!
!     version 3 .....

!     Correction  le 28/10/97    P. Le Van .
!  -------------------------------------------------------------------
#include "dimensions.h"
#include "paramet.h"
#include "parafilt.h"
!  -------------------------------------------------------------------
#include "comgeom.h"
#include "coefils.h"

      REAL  dlonu(iim),dlatu(jjm)
      REAL  rlamda( iim ),  eignvl( iim )
!

      REAL    lamdamax,pi,cof
      INTEGER i,j,modemax,imx,k,kf,ii
      REAL dymin,dxmin,colat0
      REAL eignft(iim,iim), coff
      REAL matriceun,matriceus,matricevn,matricevs,matrinvn,matrinvs 
      COMMON/matrfil/matriceun(iim,iim,nfilun),matriceus(iim,iim,nfilus) &
                   , matricevn(iim,iim,nfilvn),matricevs(iim,iim,nfilvs) &
                   ,  matrinvn(iim,iim,nfilun),matrinvs (iim,iim,nfilus)
#ifdef CRAY
      INTEGER   ISMIN
      EXTERNAL  ISMIN
      INTEGER iymin 
      INTEGER ixmineq
#endif
      EXTERNAL  inifgn
!
! ------------------------------------------------------------
!   This routine computes the eigenfunctions of the laplacien
!   on the stretched grid, and the filtering coefficients
!      
!  We designate:
!   eignfn   eigenfunctions of the discrete laplacien
!   eigenvl  eigenvalues
!   jfiltn   indexof the last scalar line filtered in NH
!   jfilts   index of the first line filtered in SH
!   modfrst  index of the mode from where modes are filtered
!   modemax  maximum number of modes ( im )
!   coefil   filtering coefficients ( lamda_max*cos(rlat)/lamda )
!   sdd      SQRT( dx )
!      
!     the modes are filtered from modfrst to modemax
      
!-----------------------------------------------------------
!

       pi       = 2. * ASIN( 1. )

       DO i = 1,iim
        dlonu(i) = xprimu( i )
       ENDDO
!
       CALL inifgn(eignvl)
!
        print *,' EIGNVL '
        PRINT 250,eignvl
250     FORMAT( 1x,5e13.6)
!
! compute eigenvalues and eigenfunctions
!
!
!.................................................................
!
!  compute the filtering coefficients for scalar lines and 
!  meridional wind v-lines
!
!  we filter all those latitude lines where coefil < 1
!  NO FILTERING AT POLES
!
!  colat0 is to be used  when alpha (stretching coefficient)
!  is set equal to zero for the regular grid case 
!
!    .......   Calcul  de  colat0   .........
!     .....  colat0 = minimum de ( 0.5, min dy/ min dx )   ...
!
!
      DO 45 j = 1,jjm
         dlatu( j ) = rlatu( j ) - rlatu( j+1 )
 45   CONTINUE
!
#ifdef CRAY
      iymin   = ISMIN( jjm, dlatu, 1 )
      ixmineq = ISMIN( iim, dlonu, 1 )
      dymin   = dlatu( iymin )
      dxmin   = dlonu( ixmineq )
#else
      dxmin   =  dlonu(1)
       DO  i  = 2, iim
        dxmin = MIN( dxmin,dlonu(i) )
       ENDDO
      dymin  = dlatu(1)
       DO j  = 2, jjm
        dymin = MIN( dymin,dlatu(j) )
       ENDDO
#endif
!
!
      colat0  =  MIN( 0.5, dymin/dxmin )
!
      IF( .NOT.fxyhypb.AND.ysinus )  THEN
           colat0 = 0.6
!         ...... a revoir  pour  ysinus !   .......
           alphax = 0.
      ENDIF
!
      PRINT 50, colat0,alphax
  50  FORMAT(/15x,' Inifilr colat0 alphax ',2e16.7)
!
      IF(alphax.EQ.1. )  THEN
        PRINT *,' Inifilr  alphax doit etre  <  a 1.  Corriger '
         STOP
      ENDIF
!
      lamdamax = iim / ( pi * colat0 * ( 1. - alphax ) )

!c                        ... Correction  le 28/10/97  ( P.Le Van ) ..
!
      DO 71 i = 2,iim
       rlamda( i ) = lamdamax/ SQRT( ABS( eignvl(i) ) )
 71   CONTINUE
!

      DO 72 j = 1,jjm
	    DO 73 i = 1,iim
	    coefilu( i,j )  = 0.0
	    coefilv( i,j )  = 0.0
	    coefilu2( i,j ) = 0.0
	    coefilv2( i,j ) = 0.0
 73     CONTINUE
 72   CONTINUE

!
!    ... Determination de jfiltnu,jfiltnv,jfiltsu,jfiltsv ....
!    .........................................................
!
       modemax = iim

!ccc    imx = modemax - 4 * (modemax/iim)

       imx  = iim
!
       PRINT *,' TRUNCATION AT ',imx
!
      DO 75 j = 2, jjm/2+1
       cof = COS( rlatu(j) )/ colat0
	    IF ( cof .LT. 1. ) THEN
          IF( rlamda(imx) * COS(rlatu(j) ).LT.1. ) jfiltnu= j
        ENDIF

       cof = COS( rlatu(jjp1-j+1) )/ colat0
	    IF ( cof .LT. 1. ) THEN
          IF( rlamda(imx) * COS(rlatu(jjp1-j+1) ).LT.1. ) &
            jfiltsu= jjp1-j+1
        ENDIF
 75   CONTINUE
!
      DO 76 j = 1, jjm/2
       cof = COS( rlatv(j) )/ colat0
	    IF ( cof .LT. 1. ) THEN
          IF( rlamda(imx) * COS(rlatv(j) ).LT.1. ) jfiltnv= j
        ENDIF

       cof = COS( rlatv(jjm-j+1) )/ colat0
	    IF ( cof .LT. 1. ) THEN
          IF( rlamda(imx) * COS(rlatv(jjm-j+1) ).LT.1. ) &
             jfiltsv= jjm-j+1
        ENDIF
 76   CONTINUE
!                                 

      IF( jfiltnu.LE.0 .OR. jfiltnu.GT. jjm/2 +1 )  THEN
        PRINT *,' jfiltnu en dehors des valeurs acceptables ' ,jfiltnu
        STOP
      ENDIF

      IF( jfiltsu.LE.0 .OR. jfiltsu.GT.  jjm  +1 )  THEN
        PRINT *,' jfiltsu en dehors des valeurs acceptables ' ,jfiltsu
        STOP
      ENDIF

      IF( jfiltnv.LE.0 .OR. jfiltnv.GT. jjm/2    )  THEN
        PRINT *,' jfiltnv en dehors des valeurs acceptables ' ,jfiltnv
        STOP
      ENDIF

      IF( jfiltsv.LE.0 .OR. jfiltsv.GT.     jjm  )  THEN
        PRINT *,' jfiltsv en dehors des valeurs acceptables ' ,jfiltsv
        STOP
      ENDIF

       PRINT *,' jfiltnv jfiltsv jfiltnu jfiltsu ' , &
                 jfiltnv,jfiltsv,jfiltnu,jfiltsu

!                                 
!   ... Determination de coefilu,coefilv,n=modfrstu,modfrstv ....
!................................................................
!
!
      DO 77 j = 1,jjm
	  modfrstu( j ) = iim
	  modfrstv( j ) = iim
 77   CONTINUE
!
      DO 84 j = 2,jfiltnu
       DO 81 k = 2,modemax
	     cof = rlamda(k) * COS( rlatu(j) )
         IF ( cof .LT. 1. ) GOTO 82
 81    CONTINUE
      GOTO 84
 82   modfrstu( j ) = k
!
	  kf = modfrstu( j )
	   DO 83 k = kf , modemax
        cof = rlamda(k) * COS( rlatu(j) )
	    coefilu(k,j) = cof - 1.
	    coefilu2(k,j) = cof*cof - 1.
 83    CONTINUE
 84   CONTINUE
!                                 
!
      DO 89 j = 1,jfiltnv
!
       DO 86 k = 2,modemax
	    cof = rlamda(k) * COS( rlatv(j) )
         IF ( cof .LT. 1. ) GOTO 87
 86    CONTINUE
      GOTO 89
 87   modfrstv( j ) = k
!
	   kf = modfrstv( j )
	   DO 88 k = kf , modemax
        cof = rlamda(k) * COS( rlatv(j) )
	    coefilv(k,j) = cof - 1.
	    coefilv2(k,j) = cof*cof - 1.
 88    CONTINUE
!
 89    CONTINUE
!
      DO 94 j = jfiltsu,jjm
       DO 91 k = 2,modemax
	    cof = rlamda(k) * COS( rlatu(j) )
         IF ( cof .LT. 1. ) GOTO 92
 91    CONTINUE
      GOTO 94
 92   modfrstu( j ) = k
!
        kf = modfrstu( j )
	 DO 93 k = kf , modemax
          cof = rlamda(k) * COS( rlatu(j) )
	  coefilu(k,j) = cof - 1.
	  coefilu2(k,j) = cof*cof - 1.
 93      CONTINUE
 94    CONTINUE
!                                 
      DO 99 j = jfiltsv,jjm
       DO 96 k = 2,modemax
	     cof = rlamda(k) * COS( rlatv(j) )
         IF ( cof .LT. 1. ) GOTO 97
 96    CONTINUE
      GOTO 99
 97   modfrstv( j ) = k
!
       kf = modfrstv( j )
	   DO 98 k = kf , modemax
        cof = rlamda(k) * COS( rlatv(j) )
	    coefilv(k,j) = cof - 1.
	    coefilv2(k,j) = cof*cof - 1.
 98    CONTINUE
 99   CONTINUE
!

       IF(jfiltnv.GE.jjm/2 .OR. jfiltnu.GE.jjm/2)THEN

         IF(jfiltnv.EQ.jfiltsv)jfiltsv=1+jfiltnv
         IF(jfiltnu.EQ.jfiltsu)jfiltsu=1+jfiltnu

          PRINT *,'jfiltnv jfiltsv jfiltnu jfiltsu' , &
              jfiltnv,jfiltsv,jfiltnu,jfiltsu
       ENDIF

       PRINT *,'   Modes premiers  v  '
       PRINT 334,modfrstv
       PRINT *,'   Modes premiers  u  '
       PRINT 334,modfrstu

     
      IF( nfilun.LT. jfiltnu )  THEN
       PRINT *,' le parametre nfilun utilise pour la matrice ', &
         ' matriceun  est trop petit ! ' 
       PRINT *,'Le changer dans parafilt.h et le mettre a  ',jfiltnu
        PRINT *,' Pour information, nfilun,nfilus,nfilvn,nfilvs ' &
       ,'doivent etre egaux successivement a  ',jfiltnu,jjm-jfiltsu+1 &
        ,jfiltnv,jjm-jfiltsv+1
               STOP
      ENDIF
      IF( nfilun.GT. jfiltnu+ 2 )  THEN
           PRINT *,' le parametre nfilun utilise pour la matrice ', &
      ' matriceun est trop grand ! Gachis de memoire ! ' 
       PRINT *,'Le changer dans parafilt.h et le mettre a  ',jfiltnu
        PRINT *,' Pour information, nfilun,nfilus,nfilvn,nfilvs ' &
       ,'doivent etre egaux successivement a  ',jfiltnu,jjm-jfiltsu+1 &
        ,jfiltnv,jjm-jfiltsv+1
!              STOP
      ENDIF
      IF( nfilus.LT. jjm - jfiltsu +1 )  THEN
            PRINT *,' le parametre nfilus utilise pour la matrice ', &
         ' matriceus  est trop petit !  '
       PRINT *,' Le changer dans parafilt.h et le mettre a  ', &
       jjm - jfiltsu + 1
        PRINT *,' Pour information , nfilun,nfilus,nfilvn,nfilvs ' &
       ,'doivent etre egaux successivement a  ',jfiltnu,jjm-jfiltsu+1 &
        ,jfiltnv,jjm-jfiltsv+1
               STOP
      ENDIF
      IF( nfilus.GT. jjm - jfiltsu + 3 )  THEN
           PRINT *,' le parametre nfilus utilise pour la matrice ', &
       ' matriceus  est trop grand ! ' 
       PRINT *,' Le changer dans parafilt.h et le mettre a  ' , &
       jjm - jfiltsu + 1
        PRINT *,' Pour information , nfilun,nfilus,nfilvn,nfilvs ' &
       ,'doivent etre egaux successivement a  ',jfiltnu,jjm-jfiltsu+1 &
        ,jfiltnv,jjm-jfiltsv+1
!              STOP
      ENDIF
      IF( nfilvn.LT. jfiltnv )  THEN
            PRINT *,' le parametre nfilvn utilise pour la matrice ', &
         ' matricevn  est trop petit ! '  
       PRINT *,'Le changer dans parafilt.h et le mettre a  ',jfiltnv
        PRINT *,' Pour information , nfilun,nfilus,nfilvn,nfilvs ' &
       ,'doivent etre egaux successivement a  ',jfiltnu,jjm-jfiltsu+1 &
        ,jfiltnv,jjm-jfiltsv+1
               STOP
      ENDIF
      IF( nfilvn.GT. jfiltnv+ 2 )  THEN
           PRINT *,' le parametre nfilvn utilise pour la matrice ', &
      ' matricevn est trop grand !  Gachis de memoire ! ' 
       PRINT *,'Le changer dans parafilt.h et le mettre a  ',jfiltnv
        PRINT *,' Pour information , nfilun,nfilus,nfilvn,nfilvs ' &
       ,'doivent etre egaux successivement a  ',jfiltnu,jjm-jfiltsu+1 &
        ,jfiltnv,jjm-jfiltsv+1
!              STOP
      ENDIF
      IF( nfilvs.LT. jjm - jfiltsv +1 )  THEN
            PRINT *,' le parametre nfilvs utilise pour la matrice ', &
         ' matricevs  est trop petit !  Le changer dans parafilt.h '
       PRINT *,' Le changer dans parafilt.h et le mettre a  ' &
       , jjm - jfiltsv + 1
        PRINT *,' Pour information , nfilun,nfilus,nfilvn,nfilvs ' &
       ,'doivent etre egaux successivement a  ',jfiltnu,jjm-jfiltsu+1 &
        ,jfiltnv,jjm-jfiltsv+1
               STOP
      ENDIF
      IF( nfilvs.GT. jjm - jfiltsv + 3 )  THEN
           PRINT *,' le parametre nfilvs utilise pour la matrice ', &
       ' matricevs  est trop grand ! Gachis de memoire ! '
       PRINT *,' Le changer dans parafilt.h et le mettre a  ' &
         ,  jjm - jfiltsv + 1
        PRINT *,' Pour information , nfilun,nfilus,nfilvn,nfilvs ' &
       ,'doivent etre egaux successivement a  ',jfiltnu,jjm-jfiltsu+1 &
        ,jfiltnv,jjm-jfiltsv+1
!              STOP
      ENDIF

!  
!   ...................................................................
!
!   ... Calcul de la matrice filtre 'matriceu'  pour les champs situes
!                       sur la grille scalaire                 ........
!   ...................................................................
!
        DO j = 2, jfiltnu

         DO i=1,iim
          coff = coefilu(i,j)
          IF( i.LT.modfrstu(j) ) coff = 0.
           DO k=1,iim
            eignft(i,k) = eignfnv(k,i) * coff
           ENDDO
         ENDDO
#ifdef CRAY
         CALL MXM( eignfnv,iim,eignft,iim,matriceun(1,1,j),iim )
#else
#ifdef BLAS
         CALL SGEMM ('N', 'N', iim, iim, iim, 1.0, &
         eignfnv, iim, eignft, iim, 0.0, matriceun(1,1,j), iim)
#else
         DO k = 1, iim
         DO i = 1, iim
            matriceun(i,k,j) = 0.0
            DO ii = 1, iim
               matriceun(i,k,j) = matriceun(i,k,j) &
                                + eignfnv(i,ii)*eignft(ii,k)
            ENDDO
         ENDDO
         ENDDO
#endif
#endif

        ENDDO

        DO j = jfiltsu, jjm

         DO i=1,iim
          coff = coefilu(i,j)
          IF( i.LT.modfrstu(j) ) coff = 0.
            DO k=1,iim
             eignft(i,k) = eignfnv(k,i) * coff
            ENDDO
         ENDDO
#ifdef CRAY
         CALL MXM(eignfnv,iim,eignft,iim,matriceus(1,1,j-jfiltsu+1),iim)
#else
#ifdef BLAS
      CALL SGEMM ('N', 'N', iim, iim, iim, 1.0, &
                 eignfnv, iim, eignft, iim, 0.0, &
                 matriceus(1,1,j-jfiltsu+1), iim)
#else
         DO k = 1, iim
         DO i = 1, iim
            matriceus(i,k,j-jfiltsu+1) = 0.0
            DO ii = 1, iim
               matriceus(i,k,j-jfiltsu+1) = matriceus(i,k,j-jfiltsu+1) &
                                          + eignfnv(i,ii)*eignft(ii,k)
            ENDDO
         ENDDO
         ENDDO
#endif
#endif

        ENDDO

!   ...................................................................
!
!   ... Calcul de la matrice filtre 'matricev'  pour les champs situes
!                       sur la grille   de V ou de Z           ........
!   ...................................................................
!
        DO j = 1, jfiltnv

         DO i = 1, iim
          coff = coefilv(i,j)
          IF( i.LT.modfrstv(j) ) coff = 0.
           DO k = 1, iim
            eignft(i,k) = eignfnu(k,i) * coff
           ENDDO
         ENDDO
#ifdef CRAY
          CALL MXM( eignfnu,iim,eignft,iim,matricevn(1,1,j),iim )
#else
#ifdef BLAS
      CALL SGEMM ('N', 'N', iim, iim, iim, 1.0, &
           eignfnu, iim, eignft, iim, 0.0, matricevn(1,1,j), iim)
#else
         DO k = 1, iim
         DO i = 1, iim
            matricevn(i,k,j) = 0.0
            DO ii = 1, iim
               matricevn(i,k,j) = matricevn(i,k,j) &
                                + eignfnu(i,ii)*eignft(ii,k)
            ENDDO
         ENDDO
         ENDDO
#endif
#endif

        ENDDO

        DO j = jfiltsv, jjm

         DO i = 1, iim
          coff = coefilv(i,j)
          IF( i.LT.modfrstv(j) ) coff = 0.
            DO k = 1, iim
             eignft(i,k) = eignfnu(k,i) * coff
            ENDDO
         ENDDO
#ifdef CRAY
         CALL MXM(eignfnu,iim,eignft,iim,matricevs(1,1,j-jfiltsv+1),iim)
#else
#ifdef BLAS
      CALL SGEMM ('N', 'N', iim, iim, iim, 1.0, &
                 eignfnu, iim, eignft, iim, 0.0, &
                 matricevs(1,1,j-jfiltsv+1), iim)
#else
         DO k = 1, iim
         DO i = 1, iim
            matricevs(i,k,j-jfiltsv+1) = 0.0
            DO ii = 1, iim
               matricevs(i,k,j-jfiltsv+1) = matricevs(i,k,j-jfiltsv+1) &
                                + eignfnu(i,ii)*eignft(ii,k)
            ENDDO
         ENDDO
         ENDDO
#endif
#endif

        ENDDO

!   ...................................................................
!
!   ... Calcul de la matrice filtre 'matrinv'  pour les champs situes
!              sur la grille scalaire , pour le filtre inverse ........
!   ...................................................................
!
        DO j = 2, jfiltnu

         DO i = 1,iim
          coff = coefilu(i,j)/ ( 1. + coefilu(i,j) )
          IF( i.LT.modfrstu(j) ) coff = 0.
           DO k=1,iim
            eignft(i,k) = eignfnv(k,i) * coff
           ENDDO
         ENDDO
#ifdef CRAY
          CALL MXM( eignfnv,iim,eignft,iim,matrinvn(1,1,j),iim )
#else
#ifdef BLAS
      CALL SGEMM ('N', 'N', iim, iim, iim, 1.0, &
           eignfnv, iim, eignft, iim, 0.0, matrinvn(1,1,j), iim)
#else
         DO k = 1, iim
         DO i = 1, iim
            matrinvn(i,k,j) = 0.0
            DO ii = 1, iim
               matrinvn(i,k,j) = matrinvn(i,k,j) &
                                + eignfnv(i,ii)*eignft(ii,k)
            ENDDO
         ENDDO
         ENDDO
#endif
#endif

        ENDDO

        DO j = jfiltsu, jjm

         DO i = 1,iim
          coff = coefilu(i,j) / ( 1. + coefilu(i,j) )
          IF( i.LT.modfrstu(j) ) coff = 0.
            DO k=1,iim
             eignft(i,k) = eignfnv(k,i) * coff
            ENDDO
         ENDDO
#ifdef CRAY
         CALL MXM(eignfnv,iim,eignft,iim,matrinvs(1,1,j-jfiltsu+1),iim)
#else
#ifdef BLAS
      CALL SGEMM ('N', 'N', iim, iim, iim, 1.0, &
       eignfnv, iim, eignft, iim, 0.0, matrinvs(1,1,j-jfiltsu+1), iim)
#else
         DO k = 1, iim
         DO i = 1, iim
            matrinvs(i,k,j-jfiltsu+1) = 0.0
            DO ii = 1, iim
               matrinvs(i,k,j-jfiltsu+1) = matrinvs(i,k,j-jfiltsu+1) &
                                + eignfnv(i,ii)*eignft(ii,k)
            ENDDO
         ENDDO
         ENDDO
#endif
#endif

        ENDDO

!   ...................................................................

!
334    FORMAT(1x,24i3)
755    FORMAT(1x,6f10.3,i3)

       RETURN
       END SUBROUTINE inifilr

END MODULE  filtreg_mod
