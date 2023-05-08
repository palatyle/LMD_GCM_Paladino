










!
! $Header$
!
      SUBROUTINE geopot (ngrid, teta, pk, pks, phis, phi )
      IMPLICIT NONE

c=======================================================================
c
c   Auteur:  P. Le Van
c   -------
c
c   Objet:
c   ------
c
c    *******************************************************************
c    ....   calcul du geopotentiel aux milieux des couches    .....
c    *******************************************************************
c
c     ....   l'integration se fait de bas en haut  ....
c
c     .. ngrid,teta,pk,pks,phis sont des argum. d'entree pour le s-pg ..
c              phi               est un  argum. de sortie pour le s-pg .
c
c  This computation (with teta = cp T / pk !) is identical to 
c     delta phi = R/RMD T/p delta p         (r=R/RMD=cpp*kappa)
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

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

c   Arguments:
c   ----------

      INTEGER ngrid
      REAL teta(ngrid,llm),pks(ngrid),phis(ngrid),pk(ngrid,llm) ,
     *       phi(ngrid,llm)


c   Local:
c   ------

      INTEGER  l, ij


c-----------------------------------------------------------------------
c     calcul de phi au niveau 1 pres du sol  .....

      DO   1  ij  = 1, ngrid
      phi( ij,1 ) = phis( ij ) + teta(ij,1) * ( pks(ij) - pk(ij,1) )
   1  CONTINUE

c     calcul de phi aux niveaux superieurs  .......

      DO  l = 2,llm
        DO  ij    = 1,ngrid
        phi(ij,l) = phi(ij,l-1) + 0.5 * ( teta(ij,l)  + teta(ij,l-1) ) 
     *                              *   (  pk(ij,l-1) -  pk(ij,l)    )
        ENDDO
      ENDDO

      RETURN
      END
