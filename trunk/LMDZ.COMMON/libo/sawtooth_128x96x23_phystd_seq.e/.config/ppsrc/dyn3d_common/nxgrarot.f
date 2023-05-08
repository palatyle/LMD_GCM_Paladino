










!
! $Header$
!
      SUBROUTINE nxgrarot (klevel,xcov, ycov, lr, grx, gry )
c   ***********************************************************
c
c    Auteur :  P.Le Van  
c
c                                 lr
c      calcul de  ( nXgrad (rot) )   du vect. v  ....
c
c       xcov et ycov  etant les compos. covariantes de  v
c   ***********************************************************
c     xcov , ycov et lr  sont des arguments  d'entree pour le s-prog
c      grx   et  gry     sont des arguments de sortie pour le s-prog
c
c
      IMPLICIT NONE
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
!
! $Header$
!
!  Attention : ce fichier include est compatible format fixe/format libre
!                 veillez à n'utiliser que des ! pour les commentaires
!                 et à bien positionner les & des lignes de continuation 
!                 (les placer en colonne 6 et en colonne 73)
!-----------------------------------------------------------------------
! INCLUDE comdissipn.h

      REAL  tetaudiv, tetaurot, tetah, cdivu, crot, cdivh
!
      COMMON/comdissipn/ tetaudiv(llm),tetaurot(llm),tetah(llm)   ,     &
     &                        cdivu,      crot,         cdivh

!
!    Les parametres de ce common proviennent des calculs effectues dans 
!             Inidissip  .
!
!-----------------------------------------------------------------------
c
      INTEGER klevel
      REAL xcov( ip1jmp1,klevel ), ycov( ip1jm,klevel )
      REAL  grx( ip1jmp1,klevel ),  gry( ip1jm,klevel )
c
      REAL rot(ip1jm,llm)

      INTEGER l,ij,iter,lr
c
c
c
      CALL SCOPY ( ip1jmp1*klevel, xcov, 1, grx, 1 )
      CALL SCOPY (  ip1jm*klevel, ycov, 1, gry, 1 )
c
      DO 10 iter = 1,lr
      CALL  rotat (klevel,grx, gry, rot )
      CALL filtreg( rot, jjm, klevel, 2,1, .false.,2)
      CALL nxgrad (klevel,rot, grx, gry )
c
      DO 5  l = 1, klevel
      DO 2 ij = 1, ip1jm
      gry( ij,l ) = - gry( ij,l ) * crot
   2  CONTINUE
      DO 3 ij = 1, ip1jmp1
      grx( ij,l ) = - grx( ij,l ) * crot
   3  CONTINUE
   5  CONTINUE
c
  10  CONTINUE
      RETURN
      END
