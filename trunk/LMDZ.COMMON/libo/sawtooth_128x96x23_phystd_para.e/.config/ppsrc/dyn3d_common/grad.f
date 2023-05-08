










!
! $Header$
!
      SUBROUTINE  grad(klevel, pg,pgx,pgy )
c
c      P. Le Van
c
c    ******************************************************************
c     .. calcul des composantes covariantes en x et y du gradient de g
c
c    ******************************************************************
c             pg        est un   argument  d'entree pour le s-prog
c       pgx  et  pgy    sont des arguments de sortie pour le s-prog
c
      IMPLICIT NONE
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
      INTEGER klevel
      REAL  pg( ip1jmp1,klevel )
      REAL pgx( ip1jmp1,klevel ) , pgy( ip1jm,klevel )
      INTEGER  l,ij
c
c
      DO 6 l = 1,klevel
c
      DO 2  ij = 1, ip1jmp1 - 1
      pgx( ij,l ) = pg( ij +1,l ) - pg( ij,l )
   2  CONTINUE
c
c    .... correction pour  pgx(ip1,j,l)  ....
c    ...    pgx(iip1,j,l)= pgx(1,j,l)  ....
CDIR$ IVDEP
      DO 3  ij = iip1, ip1jmp1, iip1
      pgx( ij,l ) = pgx( ij -iim,l )
   3  CONTINUE
c
      DO 4 ij = 1,ip1jm
      pgy( ij,l ) = pg( ij,l ) - pg( ij +iip1,l )
   4  CONTINUE
c
   6  CONTINUE
      RETURN
      END
