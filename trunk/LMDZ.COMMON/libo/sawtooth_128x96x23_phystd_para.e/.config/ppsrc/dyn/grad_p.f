










      SUBROUTINE  grad_p(klevel, pg,pgx,pgy )
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
      USE parallel_lmdz
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
      INTEGER :: ijb,ije,jjb,jje
c
c
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO 6 l = 1,klevel
c
      ijb=ij_begin
      ije=ij_end
      DO 2  ij = ijb, ije - 1
        pgx( ij,l ) = pg( ij +1,l ) - pg( ij,l )
   2  CONTINUE
c
c    .... correction pour  pgx(ip1,j,l)  ....
c    ...    pgx(iip1,j,l)= pgx(1,j,l)  ....
CDIR$ IVDEP
      DO 3  ij = ijb+iip1-1, ije, iip1
        pgx( ij,l ) = pgx( ij -iim,l )
   3  CONTINUE
c
      ijb=ij_begin-iip1
      ije=ij_end
      if (pole_nord) ijb=ij_begin
      if (pole_sud)  ije=ij_end-iip1
      
      DO 4 ij = ijb,ije
        pgy( ij,l ) = pg( ij,l ) - pg( ij +iip1,l )
   4  CONTINUE
c
   6  CONTINUE
c$OMP END DO NOWAIT

      RETURN
      END
