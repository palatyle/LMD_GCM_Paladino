!
! $Id $
!
MODULE filtreg_mod

  REAL, DIMENSION(:,:,:), ALLOCATABLE :: matriceun,matriceus,matricevn
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: matricevs,matrinvn,matrinvs

CONTAINS

  SUBROUTINE inifilr
#ifdef CPP_PARA
  USE mod_filtre_fft, ONLY : use_filtre_fft,Init_filtre_fft
  USE mod_filtre_fft_loc, ONLY : Init_filtre_fft_loc=>Init_filtre_fft    !
#endif
  USE logic_mod, ONLY: fxyhypb,ysinus
  USE serre_mod, ONLY: alphax
    !    ... H. Upadhyaya, O.Sharma   ...
    !
    IMPLICIT NONE
    !
    !     version 3 .....

    !     Correction  le 28/10/97    P. Le Van .
    !  -------------------------------------------------------------------
#include "dimensions.h"
#include "paramet.h"
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

    LOGICAL, SAVE :: first_call_inifilr = .TRUE.

#ifdef CRAY
    INTEGER   ISMIN
    EXTERNAL  ISMIN
    INTEGER iymin 
    INTEGER ixmineq
#endif
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
    !   modfrst  index of the mode from WHERE modes are filtered
    !   modemax  maximum number of modes ( im )
    !   coefil   filtering coefficients ( lamda_max*COS(rlat)/lamda )
    !   sdd      SQRT( dx )
    !      
    !     the modes are filtered from modfrst to modemax
    !      
    !-----------------------------------------------------------
    !

    pi       = 2. * ASIN( 1. )

    DO i = 1,iim
       dlonu(i) = xprimu( i )
    ENDDO
    !
    CALL inifgn(eignvl)
    !
    PRINT *,'inifilr: EIGNVL '
    PRINT 250,eignvl
250 FORMAT( 1x,5e14.6)
    !
    ! compute eigenvalues and eigenfunctions
    !
    !
    !.................................................................
    !
    !  compute the filtering coefficients for scalar lines and 
    !  meridional wind v-lines
    !
    !  we filter all those latitude lines WHERE coefil < 1
    !  NO FILTERING AT POLES
    !
    !  colat0 is to be used  when alpha (stretching coefficient)
    !  is set equal to zero for the regular grid CASE 
    !
    !    .......   Calcul  de  colat0   .........
    !     .....  colat0 = minimum de ( 0.5, min dy/ min dx )   ...
    !
    !
    DO j = 1,jjm
       dlatu( j ) = rlatu( j ) - rlatu( j+1 )
    ENDDO
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
    ! For a regular grid, we want the filter to start at latitudes
    ! corresponding to lengths dx of the same size as dy (in terms
    ! of angles: dx=2*dy) => at colat0=0.5 (i.e. colatitude=30 degrees
    !  <=> latitude=60 degrees).
    ! Same idea for the zoomed grid: start filtering polewards as soon
    ! as length dx becomes of the same size as dy 
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

    !                        ... Correction  le 28/10/97  ( P.Le Van ) ..
    !
    DO i = 2,iim
       rlamda( i ) = lamdamax/ SQRT( ABS( eignvl(i) ) )
    ENDDO
    !

    DO j = 1,jjm
       DO i = 1,iim
          coefilu( i,j )  = 0.0
          coefilv( i,j )  = 0.0
          coefilu2( i,j ) = 0.0
          coefilv2( i,j ) = 0.0
       ENDDO
    ENDDO

    !
    !    ... Determination de jfiltnu,jfiltnv,jfiltsu,jfiltsv ....
    !    .........................................................
    !
    modemax = iim

!!!!    imx = modemax - 4 * (modemax/iim)

    imx  = iim
    !
    PRINT *,'inifilr: TRUNCATION AT ',imx
    !
! Ehouarn: set up some defaults
    jfiltnu=2 ! avoid north pole
    jfiltsu=jjm ! avoid south pole (which is at jjm+1)
    jfiltnv=1 ! NB: no poles on the V grid
    jfiltsv=jjm

    DO j = 2, jjm/2+1
       cof = COS( rlatu(j) )/ colat0
       IF ( cof .LT. 1. ) THEN
          IF( rlamda(imx) * COS(rlatu(j) ).LT.1. ) THEN
            jfiltnu= j
          ENDIF
       ENDIF

       cof = COS( rlatu(jjp1-j+1) )/ colat0
       IF ( cof .LT. 1. ) THEN
          IF( rlamda(imx) * COS(rlatu(jjp1-j+1) ).LT.1. ) THEN
               jfiltsu= jjp1-j+1
          ENDIF
       ENDIF
    ENDDO
    !
    DO j = 1, jjm/2
       cof = COS( rlatv(j) )/ colat0
       IF ( cof .LT. 1. ) THEN
          IF( rlamda(imx) * COS(rlatv(j) ).LT.1. ) THEN
            jfiltnv= j
          ENDIF
       ENDIF

       cof = COS( rlatv(jjm-j+1) )/ colat0
       IF ( cof .LT. 1. ) THEN
          IF( rlamda(imx) * COS(rlatv(jjm-j+1) ).LT.1. ) THEN
               jfiltsv= jjm-j+1
          ENDIF
       ENDIF
    ENDDO
    !                                 

    IF( jfiltnu.GT. jjm/2 +1 )  THEN
       PRINT *,' jfiltnu en dehors des valeurs acceptables ' ,jfiltnu
       STOP
    ENDIF

    IF( jfiltsu.GT.  jjm  +1 )  THEN
       PRINT *,' jfiltsu en dehors des valeurs acceptables ' ,jfiltsu
       STOP
    ENDIF

    IF( jfiltnv.GT. jjm/2    )  THEN
       PRINT *,' jfiltnv en dehors des valeurs acceptables ' ,jfiltnv
       STOP
    ENDIF

    IF( jfiltsv.GT.     jjm  )  THEN
       PRINT *,' jfiltsv en dehors des valeurs acceptables ' ,jfiltsv
       STOP
    ENDIF

    PRINT *,'inifilr: jfiltnv jfiltsv jfiltnu jfiltsu ' , &
         jfiltnv,jfiltsv,jfiltnu,jfiltsu

    IF(first_call_inifilr) THEN
       ALLOCATE(matriceun(iim,iim,jfiltnu))
       ALLOCATE(matriceus(iim,iim,jjm-jfiltsu+1))
       ALLOCATE(matricevn(iim,iim,jfiltnv))
       ALLOCATE(matricevs(iim,iim,jjm-jfiltsv+1))
       ALLOCATE( matrinvn(iim,iim,jfiltnu))
       ALLOCATE( matrinvs(iim,iim,jjm-jfiltsu+1))
       first_call_inifilr = .FALSE.
    ENDIF

    !                                 
    !   ... Determination de coefilu,coefilv,n=modfrstu,modfrstv ....
    !................................................................
    !
    !
    DO j = 1,jjm
    !default initialization: all modes are retained (i.e. no filtering)
       modfrstu( j ) = iim
       modfrstv( j ) = iim
    ENDDO
    !
    DO j = 2,jfiltnu
       DO k = 2,modemax
          cof = rlamda(k) * COS( rlatu(j) )
          IF ( cof .LT. 1. ) GOTO 82
       ENDDO
       GOTO 84
82     modfrstu( j ) = k
       !
       kf = modfrstu( j )
       DO k = kf , modemax
          cof = rlamda(k) * COS( rlatu(j) )
          coefilu(k,j) = cof - 1.
          coefilu2(k,j) = cof*cof - 1.
       ENDDO
84     CONTINUE
    ENDDO
    !                                 
    !
    DO j = 1,jfiltnv
       !
       DO k = 2,modemax
          cof = rlamda(k) * COS( rlatv(j) )
          IF ( cof .LT. 1. ) GOTO 87
       ENDDO
       GOTO 89
87     modfrstv( j ) = k
       !
       kf = modfrstv( j )
       DO k = kf , modemax
          cof = rlamda(k) * COS( rlatv(j) )
          coefilv(k,j) = cof - 1.
          coefilv2(k,j) = cof*cof - 1.
       ENDDO
89     CONTINUE
    ENDDO
    !
    DO j = jfiltsu,jjm
       DO k = 2,modemax
          cof = rlamda(k) * COS( rlatu(j) )
          IF ( cof .LT. 1. ) GOTO 92
       ENDDO
       GOTO 94
92     modfrstu( j ) = k
       !
       kf = modfrstu( j )
       DO k = kf , modemax
          cof = rlamda(k) * COS( rlatu(j) )
          coefilu(k,j) = cof - 1.
          coefilu2(k,j) = cof*cof - 1.
       ENDDO
94     CONTINUE
    ENDDO
    !                                 
    DO j = jfiltsv,jjm
       DO k = 2,modemax
          cof = rlamda(k) * COS( rlatv(j) )
          IF ( cof .LT. 1. ) GOTO 97
       ENDDO
       GOTO 99
97     modfrstv( j ) = k
       !
       kf = modfrstv( j )
       DO k = kf , modemax
          cof = rlamda(k) * COS( rlatv(j) )
          coefilv(k,j) = cof - 1.
          coefilv2(k,j) = cof*cof - 1.
       ENDDO
99     CONTINUE
    ENDDO
    !

    IF(jfiltnv.GE.jjm/2 .OR. jfiltnu.GE.jjm/2)THEN
! Ehouarn: and what are these for??? Trying to handle a limit case
!          where filters extend to and meet at the equator?
       IF(jfiltnv.EQ.jfiltsv)jfiltsv=1+jfiltnv
       IF(jfiltnu.EQ.jfiltsu)jfiltsu=1+jfiltnu

       PRINT *,'jfiltnv jfiltsv jfiltnu jfiltsu' , &
            jfiltnv,jfiltsv,jfiltnu,jfiltsu
    ENDIF

    PRINT *,'   Modes premiers  v  '
    PRINT 334,modfrstv
    PRINT *,'   Modes premiers  u  '
    PRINT 334,modfrstu

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
       ENDDO ! of DO i=1,iim
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
       ENDDO ! of DO k = 1, iim
#endif
#endif

    ENDDO ! of DO j = 2, jfiltnu

    DO j = jfiltsu, jjm

       DO i=1,iim
          coff = coefilu(i,j)
          IF( i.LT.modfrstu(j) ) coff = 0.
          DO k=1,iim
             eignft(i,k) = eignfnv(k,i) * coff
          ENDDO
       ENDDO ! of DO i=1,iim
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
       ENDDO ! of DO k = 1, iim
#endif
#endif

    ENDDO ! of DO j = jfiltsu, jjm

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

    ENDDO ! of DO j = 1, jfiltnv

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

    ENDDO ! of DO j = jfiltsv, jjm

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

    ENDDO ! of DO j = 2, jfiltnu

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

    ENDDO ! of DO j = jfiltsu, jjm

#ifdef CPP_PARA
    IF (use_filtre_fft) THEN
       CALL Init_filtre_fft(coefilu,modfrstu,jfiltnu,jfiltsu,  &
                           coefilv,modfrstv,jfiltnv,jfiltsv)
       CALL Init_filtre_fft_loc(coefilu,modfrstu,jfiltnu,jfiltsu,  &
                           coefilv,modfrstv,jfiltnv,jfiltsv)
    ENDIF
#endif
    !   ...................................................................

    !
334 FORMAT(1x,24i3)
755 FORMAT(1x,6f10.3,i3)

    RETURN
  END SUBROUTINE inifilr

END MODULE filtreg_mod
