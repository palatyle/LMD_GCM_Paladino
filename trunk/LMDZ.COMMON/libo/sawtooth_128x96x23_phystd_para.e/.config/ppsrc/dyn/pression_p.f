










      SUBROUTINE pression_p( ngrid, ap, bp, ps, p )
      USE parallel_lmdz
c

c      Auteurs : P. Le Van , Fr.Hourdin  .

c  ************************************************************************
c     Calcule la pression p(l) aux differents niveaux l = 1 ( niveau du
c     sol) a l = llm +1 ,ces niveaux correspondant aux interfaces des (llm) 
c     couches , avec  p(ij,llm +1) = 0.  et p(ij,1) = ps(ij)  .      
c  ************************************************************************
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
c
      INTEGER ngrid
      INTEGER l,ij
 
      REAL ap( llmp1 ), bp( llmp1 ), ps( ngrid ), p( ngrid,llmp1 ) 
      
      INTEGER ijb,ije

      
      ijb=ij_begin-iip1
      ije=ij_end+2*iip1
      
      if (pole_nord) ijb=ij_begin
      if (pole_sud)  ije=ij_end

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
      DO    l    = 1, llmp1
        DO  ij   = ijb, ije
         p(ij,l) = ap(l) + bp(l) * ps(ij)
        ENDDO
      ENDDO
c$OMP END DO NOWAIT   
      RETURN
      END
