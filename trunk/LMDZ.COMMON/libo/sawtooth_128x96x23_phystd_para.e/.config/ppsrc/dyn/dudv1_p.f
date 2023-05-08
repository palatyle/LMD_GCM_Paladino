










      SUBROUTINE dudv1_p ( vorpot, pbaru, pbarv, du, dv )
      USE parallel_lmdz
      IMPLICIT NONE
c
c-----------------------------------------------------------------------
c
c   Auteur:   P. Le Van
c   -------
c
c   Objet:
c   ------
c   calcul du terme de  rotation
c   ce terme est ajoute a  d(ucov)/dt et a d(vcov)/dt  ..
c   vorpot, pbaru et pbarv sont des arguments d'entree  pour le s-pg ..
c   du  et dv              sont des arguments de sortie pour le s-pg ..
c
c-----------------------------------------------------------------------

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

      REAL vorpot( ip1jm,llm ) ,pbaru( ip1jmp1,llm ) ,
     *     pbarv( ip1jm,llm ) ,du( ip1jmp1,llm ) ,dv( ip1jm,llm )
      INTEGER  l,ij,ijb,ije
c
c
      
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)       
      DO 10 l = 1,llm
c
      ijb=ij_begin
      ije=ij_end
      
      if (pole_nord) ijb=ij_begin+iip1
      if (pole_sud)  ije=ij_end-iip1
      
      DO 2  ij = ijb, ije-1 
      du( ij,l ) = 0.125 *(  vorpot(ij-iip1, l) + vorpot( ij, l)  ) *
     *                    (   pbarv(ij-iip1, l) + pbarv(ij-iim,  l) +
     *                        pbarv(   ij  , l) + pbarv(ij+ 1 ,  l)   )
   2  CONTINUE
   
 
c
      if (pole_nord) ijb=ij_begin
      
      DO 3 ij = ijb, ije-1 
      dv( ij+1,l ) = - 0.125 *(  vorpot(ij, l)  + vorpot(ij+1, l)  ) *
     *                        (   pbaru(ij, l)  +  pbaru(ij+1   , l) +
     *                       pbaru(ij+iip1, l)  +  pbaru(ij+iip2, l)  )
   3  CONTINUE
c
c    .... correction  pour  dv( 1,j,l )  .....
c    ....   dv(1,j,l)= dv(iip1,j,l) ....
c
CDIR$ IVDEP
      DO 4 ij = ijb, ije, iip1
      dv( ij,l ) = dv( ij + iim, l )
   4  CONTINUE
c
  10  CONTINUE
c$OMP END DO NOWAIT
      RETURN
      END
