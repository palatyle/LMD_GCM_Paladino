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
#include "dimensions.h"
#include "paramet.h"
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
