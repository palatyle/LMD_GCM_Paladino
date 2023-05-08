










      SUBROUTINE bernoui_p (ngrid,nlay,pphi,pecin,pbern)
      USE parallel_lmdz
      IMPLICIT NONE

c=======================================================================
c
c   Auteur:   P. Le Van
c   -------
c
c   Objet:
c   ------
c     calcul de la fonction de Bernouilli aux niveaux s  .....
c     phi  et  ecin  sont des arguments d'entree pour le s-pg .......
c          bern       est un  argument de sortie pour le s-pg  ......
c
c    fonction de Bernouilli = bern = filtre de( geopotentiel + 
c                              energ.cinet.)
c
c=======================================================================
c
c-----------------------------------------------------------------------
c   Decalrations:
c   -------------
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
c   Arguments:
c   ----------
c
      INTEGER nlay,ngrid
      REAL pphi(ngrid,nlay),pecin(ngrid,nlay),pbern(ngrid,nlay)
c
c   Local:
c   ------
c
      INTEGER   ij,l,ijb,ije,jjb,jje
c
c-----------------------------------------------------------------------
c   calcul de Bernouilli:
c   ---------------------
c
      ijb=ij_begin
      ije=ij_end+iip1
      if (pole_sud) ije=ij_end

      jjb=jj_begin
      jje=jj_end+1
      if (pole_sud) jje=jj_end

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)                
      DO l=1,llm
    
        DO 4 ij = ijb,ije
          pbern( ij,l ) =  pphi( ij,l ) + pecin( ij,l )
   4    CONTINUE
       
       ENDDO
c$OMP END DO NOWAIT
c
c-----------------------------------------------------------------------
c   filtre:
c   -------
c

        
        CALL filtreg_p( pbern,jjb,jje, jjp1, llm, 2,1, .true., 1 )
c
c-----------------------------------------------------------------------
      
      
      RETURN
      END
