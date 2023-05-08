










!
! $Id: fxysinus.F 1403 2010-07-01 09:02:53Z fairhead $
!
      SUBROUTINE fxysinus (rlatu,yprimu,rlatv,yprimv,rlatu1,yprimu1,
     ,                    rlatu2,yprimu2,
     ,  rlonu,xprimu,rlonv,xprimv,rlonm025,xprimm025,rlonp025,xprimp025)

      USE comconst_mod, ONLY: pi

      IMPLICIT NONE
c
c     Calcul  des longitudes et des latitudes  pour une fonction f(x,y)
c            avec y = Asin( j )  .
c
c     Auteur  :  P. Le Van
c
c
!-----------------------------------------------------------------------
!   INCLUDE 'dimensions.h'
!
!   dimensions.h contient les dimensions du modele
!   ndm est tel que iim=2**ndm
!-----------------------------------------------------------------------

      INTEGER iim,jjm,llm,ndm

      PARAMETER (iim= 128,jjm=96,llm=23,ndm=1)

!-----------------------------------------------------------------------
!
! $Header$
!
!
!  ATTENTION!!!!: ce fichier include est compatible format fixe/format libre
!                 veillez  n'utiliser que des ! pour les commentaires
!                 et  bien positionner les & des lignes de continuation
!                 (les placer en colonne 6 et en colonne 73)
!
!
!-----------------------------------------------------------------------
!   INCLUDE 'paramet.h'

      INTEGER  iip1,iip2,iip3,jjp1,llmp1,llmp2,llmm1
      INTEGER  kftd,ip1jm,ip1jmp1,ip1jmi1,ijp1llm
      INTEGER  ijmllm,mvar
      INTEGER jcfil,jcfllm

      PARAMETER( iip1= iim+1,iip2=iim+2,iip3=iim+3                       &
     &    ,jjp1=jjm+1-1/jjm)
      PARAMETER( llmp1 = llm+1,  llmp2 = llm+2, llmm1 = llm-1 )
      PARAMETER( kftd  = iim/2 -ndm )
      PARAMETER( ip1jm  = iip1*jjm,  ip1jmp1= iip1*jjp1 )
      PARAMETER( ip1jmi1= ip1jm - iip1 )
      PARAMETER( ijp1llm= ip1jmp1 * llm, ijmllm= ip1jm * llm )
      PARAMETER( mvar= ip1jmp1*( 2*llm+1) + ijmllm )
      PARAMETER( jcfil=jjm/2+5, jcfllm=jcfil*llm )

!-----------------------------------------------------------------------

       INTEGER i,j

       REAL rlatu(jjp1), yprimu(jjp1),rlatv(jjm), yprimv(jjm),
     , rlatu1(jjm), yprimu1(jjm), rlatu2(jjm), yprimu2(jjm)
       REAL rlonu(iip1),xprimu(iip1),rlonv(iip1),xprimv(iip1),
     , rlonm025(iip1),xprimm025(iip1), rlonp025(iip1),xprimp025(iip1)

!
! $Header$
!
c-----------------------------------------------------------------------
c INCLUDE 'fxyprim.h'
c
c    ................................................................
c    ................  Fonctions in line  ...........................
c    ................................................................
c
      REAL  fy, fx, fxprim, fyprim
      REAL  ri, rj
c
c
      fy(rj)=ASIN(1.+2.*((1.-rj)/REAL(jjm)))
      fyprim(rj)=1./SQRT((rj-1.)*(jjm+1.-rj))

      fx    ( ri ) = 2.*pi/REAL(iim) * ( ri - 0.5*  REAL(iim) - 1. )
c     fx    ( ri ) = 2.*pi/REAL(iim) * ( ri - 0.5* ( REAL(iim) + 1.) )
      fxprim( ri ) = 2.*pi/REAL(iim)
c
c
c    La valeur de pi est passee par le common/const/ou /const2/ .
c    Sinon, il faut la calculer avant d'appeler ces fonctions .
c
c   ----------------------------------------------------------------
c     Fonctions a changer eventuellement, selon x(x) et y(y) choisis .
c   -----------------------------------------------------------------
c
c    .....  ici, on a l'application particuliere suivante   ........
c
c                **************************************
c                **     x = 2. * pi/iim *  X         **
c                **     y =      pi/jjm *  Y         **
c                **************************************
c
c   ..................................................................
c   ..................................................................
c
c
c
c-----------------------------------------------------------------------


c    ......  calcul  des  latitudes  et de y'   .....
c
       DO j = 1, jjm + 1 
          rlatu(j) = fy    ( REAL( j )        )
         yprimu(j) = fyprim( REAL( j )        )
       ENDDO


       DO j = 1, jjm

         rlatv(j)  = fy    ( REAL( j ) + 0.5  )
         rlatu1(j) = fy    ( REAL( j ) + 0.25 ) 
         rlatu2(j) = fy    ( REAL( j ) + 0.75 ) 

        yprimv(j)  = fyprim( REAL( j ) + 0.5  ) 
        yprimu1(j) = fyprim( REAL( j ) + 0.25 )
        yprimu2(j) = fyprim( REAL( j ) + 0.75 )

       ENDDO

c
c     .....  calcul   des  longitudes et de  x'   .....
c
       DO i = 1, iim + 1
           rlonv(i)     = fx    (   REAL( i )          )
           rlonu(i)     = fx    (   REAL( i ) + 0.5    )
        rlonm025(i)     = fx    (   REAL( i ) - 0.25  )
        rlonp025(i)     = fx    (   REAL( i ) + 0.25  )

         xprimv  (i)    = fxprim (  REAL( i )          )
         xprimu  (i)    = fxprim (  REAL( i ) + 0.5    )
        xprimm025(i)    = fxprim (  REAL( i ) - 0.25   )
        xprimp025(i)    = fxprim (  REAL( i ) + 0.25   )
       ENDDO

c
       RETURN
       END

