










!
! $Header$
!
      SUBROUTINE test_period ( ucov, vcov, teta, q, p, phis )
c
c     Auteur : P. Le Van  
c    ---------
c  ....  Cette routine teste la periodicite en longitude des champs   ucov,
c                           teta, q , p et phis                 .......... 
c
      USE infotrac
c     IMPLICIT NONE
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
c    ......  Arguments   ......
c
      REAL ucov(ip1jmp1,llm), vcov(ip1jm,llm), teta(ip1jmp1,llm) ,
     ,      q(ip1jmp1,llm,nqtot), p(ip1jmp1,llmp1), phis(ip1jmp1)
c
c   .....  Variables  locales  .....
c
      INTEGER ij,l,nq
c
      DO l = 1, llm
         DO ij = 1, ip1jmp1, iip1
          IF( ucov(ij,l).NE.ucov(ij+iim,l) )  THEN
          PRINT *,'STOP dans test_period car ---  UCOV  ---  n est pas',  
     ,  ' periodique en longitude ! '
          PRINT *,' l,  ij = ', l, ij, ij+iim
          STOP
          ENDIF
          IF( teta(ij,l).NE.teta(ij+iim,l) )  THEN
          PRINT *,'STOP dans test_period car ---  TETA  ---  n est pas',  
     ,   ' periodique en longitude ! '
          PRINT *,' l,  ij = ', l, ij, ij+iim
     ,      , teta(ij,l),   teta(ij+iim,l)
          STOP
          ENDIF
         ENDDO

         do ij=1,iim
          if (teta(ij,l).ne.teta(1,l)
     s     .or.teta(ip1jm+ij,l).ne.teta(ip1jm+1,l) ) then
          PRINT *,'STOP dans test_period car ---  TETA  ---  n est pas',  
     ,  ' constant aux poles ! '
          print*,'teta(',1 ,',',l,')=',teta(1 ,l)
          print*,'teta(',ij,',',l,')=',teta(ij,l)
          print*,'teta(',ip1jm+1 ,',',l,')=',teta(ip1jm+1 ,l)
          print*,'teta(',ip1jm+ij,',',l,')=',teta(ip1jm+ij,l)
          stop
          endif
         enddo
      ENDDO

c
      DO l = 1, llm
         DO ij = 1, ip1jm, iip1
          IF( vcov(ij,l).NE.vcov(ij+iim,l) )  THEN
          PRINT *,'STOP dans test_period car ---  VCOV  ---  n est pas',  
     ,   ' periodique en longitude !'
          PRINT *,' l,  ij = ', l, ij, ij+iim,vcov(ij+iim,l),vcov(ij,l)
          vcov(ij+iim,l)=vcov(ij,l)
c         STOP
          ENDIF
         ENDDO
      ENDDO
      
c
      DO nq =1, nqtot
        DO l =1, llm
          DO ij = 1, ip1jmp1, iip1
          IF( q(ij,l,nq).NE.q(ij+iim,l,nq) )  THEN
          PRINT *,'STOP dans test_period car ---  Q  ---  n est pas ',  
     ,   'periodique en longitude !'
          PRINT *,' nq , l,  ij = ', nq, l, ij, ij+iim
          STOP
          ENDIF
          ENDDO
        ENDDO
      ENDDO
c
       DO l = 1, llm
         DO ij = 1, ip1jmp1, iip1
          IF( p(ij,l).NE.p(ij+iim,l) )  THEN
          PRINT *,'STOP dans test_period car ---  P  ---  n est pas',  
     ,    ' periodique en longitude !'
          PRINT *,' l ij = ',l, ij, ij+iim
          STOP
          ENDIF
          IF( phis(ij).NE.phis(ij+iim) )  THEN
          PRINT *,'STOP dans test_period car ---  PHIS  ---  n est pas',  
     ,   ' periodique en longitude !  l, IJ = ', l, ij,ij+iim
          PRINT *,' ij = ', ij, ij+iim
          STOP
          ENDIF
         ENDDO
         do ij=1,iim
          if (p(ij,l).ne.p(1,l)
     s     .or.p(ip1jm+ij,l).ne.p(ip1jm+1,l) ) then
          PRINT *,'STOP dans test_period car ---  P     ---  n est pas',  
     ,  ' constant aux poles ! '
          print*,'p(',1 ,',',l,')=',p(1 ,l)
          print*,'p(',ij,',',l,')=',p(ij,l)
          print*,'p(',ip1jm+1 ,',',l,')=',p(ip1jm+1 ,l)
          print*,'p(',ip1jm+ij,',',l,')=',p(ip1jm+ij,l)
          stop
          endif
         enddo
       ENDDO
c
c
         RETURN
         END
