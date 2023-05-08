










      SUBROUTINE vitvert_p ( convm , w )
c
      USE parallel_lmdz
      USE comvert_mod, ONLY: bp
      IMPLICIT NONE

c=======================================================================
c
c   Auteurs:  P. Le Van , F. Hourdin .
c   -------
c
c   Objet:
c   ------
c
c    *******************************************************************
c  .... calcul de la vitesse verticale aux niveaux sigma  ....
c    *******************************************************************
c     convm   est un argument  d'entree pour le s-pg  ......
c       w     est un argument de sortie pour le s-pg  ......
c
c    la vitesse verticale est orientee de  haut en bas .
c    au sol, au niveau sigma(1),   w(i,j,1) = 0.
c    au sommet, au niveau sigma(llm+1) , la vit.verticale est aussi
c    egale a 0. et n'est pas stockee dans le tableau w  .
c
c
c=======================================================================

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

      REAL w(ip1jmp1,llm),convm(ip1jmp1,llm)
      INTEGER   l, ij,ijb,ije


      ijb=ij_begin
      ije=ij_end+iip1
      
      if (pole_sud) ije=ij_end
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      DO 2  l = 1,llmm1

      DO 1 ij = ijb,ije
      w( ij, l+1 ) = convm( ij, l+1 ) - bp(l+1) * convm( ij, 1 )
   1  CONTINUE

   2  CONTINUE
c$OMP END DO
c$OMP MASTER
      DO 5 ij  = ijb,ije
      w(ij,1)  = 0.
5     CONTINUE
c$OMP END MASTER
c$OMP BARRIER
      RETURN
      END
