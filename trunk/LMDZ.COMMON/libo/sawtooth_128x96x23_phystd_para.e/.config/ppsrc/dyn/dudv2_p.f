










      SUBROUTINE dudv2_p ( teta, pkf, bern, du, dv  )
      USE parallel_lmdz
      IMPLICIT NONE
c
c=======================================================================
c
c   Auteur:  P. Le Van
c   -------
c
c   Objet:
c   ------
c
c   *****************************************************************
c   ..... calcul du terme de pression (gradient de p/densite )   et
c          du terme de ( -gradient de la fonction de Bernouilli ) ...
c   *****************************************************************
c          Ces termes sont ajoutes a  d(ucov)/dt et a d(vcov)/dt  ..
c
c
c    teta , pkf, bern  sont des arguments d'entree  pour le s-pg  ....
c    du et dv          sont des arguments de sortie pour le s-pg  ....
c
c=======================================================================
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

      REAL teta( ip1jmp1,llm ),pkf( ip1jmp1,llm ) ,bern( ip1jmp1,llm ),
     *         du( ip1jmp1,llm ),  dv( ip1jm,llm )
      INTEGER  l,ij,ijb,ije
c
c
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO 5 l = 1,llm
c
      ijb=ij_begin
      ije=ij_end
      if (pole_nord) ijb=ijb+iip1
      if (pole_sud)  ije=ije-iip1

      DO 2  ij  = ijb, ije - 1
       du(ij,l) = du(ij,l) + 0.5* ( teta( ij,l ) + teta( ij+1,l ) ) *
     * ( pkf( ij,l ) - pkf(ij+1,l) )  + bern(ij,l) - bern(ij+1,l)
   2  CONTINUE
c
c
c    .....  correction  pour du(iip1,j,l),  j=2,jjm   ......
c    ...          du(iip1,j,l) = du(1,j,l)                 ...
c
CDIR$ IVDEP
      DO 3 ij = ijb+iip1-1, ije, iip1
      du( ij,l ) = du( ij - iim,l )
   3  CONTINUE
c
c
      if (pole_nord) ijb=ijb-iip1

      DO 4 ij  = ijb,ije
      dv( ij,l) = dv(ij,l) + 0.5 * ( teta(ij,l) + teta( ij+iip1,l ) ) *
     *                             ( pkf(ij+iip1,l) - pkf(  ij,l  ) )
     *                           +   bern( ij+iip1,l ) - bern( ij  ,l )
   4  CONTINUE
c
   5  CONTINUE
c$OMP END DO NOWAIT 
c
      RETURN
      END
