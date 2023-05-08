










!
! $Header$
!
      SUBROUTINE rotatst (klevel,x, y, rot )
c
c  P. Le Van
c
c    *****************************************************************
c     .. calcule le rotationnel a tous les niveaux d'1 vecteur de comp. x et y ..
c         x  et  y etant des composantes  covariantes  .....
c    *****************************************************************
c        x  et y     sont des arguments d'entree pour le s-prog
c        rot          est  un argument  de sortie pour le s-prog
c
      IMPLICIT NONE
c
      INTEGER klevel
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

      REAL rot( ip1jm,klevel )
      REAL x( ip1jmp1,klevel ), y( ip1jm,klevel )
      INTEGER  l, ij
c
c
      DO 5 l = 1,klevel
c
      DO 1 ij = 1, ip1jm - 1
      rot( ij,l )  =  (  y( ij+1 , l )  -  y( ij,l )   +
     *                 x(ij +iip1, l )  -  x( ij,l )  )
   1  CONTINUE
c
c    .... correction pour rot( iip1,j,l)  ....
c
c    ....   rot(iip1,j,l)= rot(1,j,l) ...
CDIR$ IVDEP
      DO 2 ij = iip1, ip1jm, iip1
      rot( ij,l ) = rot( ij -iim,l )
   2  CONTINUE
c
   5  CONTINUE
      RETURN
      END
