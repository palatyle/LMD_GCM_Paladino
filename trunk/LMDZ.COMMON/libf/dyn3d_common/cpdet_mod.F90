module cpdet_mod

implicit none

! ADAPTATION OF GCM TO CP(T)
!======================================================================
! S. Lebonnois, 10/2010
!
! Cp must be computed using cpdet(t) to be valid
!
! The Exner function is still pk = RCPD*(play/pref)**RKAPPA
! (RCPD=cpp, RKAPPA=kappa)
!
! One goes from T to teta (potential temperature) using t2tpot(t,teta,pk)
! One goes from teta to T using tpot2t(teta,t,pk)
!
!======================================================================

contains

      SUBROUTINE ini_cpdet
      
      USE control_mod, ONLY: cpofT 
      USE comconst_mod, ONLY: nu_venus,t0_venus
      IMPLICIT none
!======================================================================
! Initialization of nu_venus and t0_venus
!======================================================================

      if (cpofT) then
          nu_venus=0.35
          t0_venus=460.
      else
          nu_venus=0.
          t0_venus=0.
      endif

      return
      end subroutine ini_cpdet

!======================================================================
!======================================================================

      FUNCTION cpdet(t)

      USE control_mod, ONLY: cpofT 
      USE comconst_mod, ONLY: cpp,t0_venus,nu_venus
      IMPLICIT none

! for cpp, nu_venus and t0_venus:

      real,intent(in) :: t
      real cpdet

      if (cpofT) then
          cpdet = cpp*(t/t0_venus)**nu_venus
      else
          cpdet = cpp
      endif

      return
      end function cpdet
      
!======================================================================
!======================================================================

      SUBROUTINE t2tpot(npoints, yt, yteta, ypk)
!======================================================================
! Arguments:
!
! yt   --------input-R- Temperature
! yteta-------output-R- Temperature potentielle
! ypk  --------input-R- Fonction d'Exner: RCPD*(pplay/pref)**RKAPPA
!
!======================================================================

      USE control_mod, ONLY: cpofT 
      USE comconst_mod, ONLY: cpp,t0_venus,nu_venus

      IMPLICIT NONE

      integer,intent(in) :: npoints
      REAL,intent(in) :: yt(npoints), ypk(npoints)
      REAL,intent(out) :: yteta(npoints)
      
!----------------------
! ATMOSPHERE PROFONDE
      real ypklim,ypklim2,mmm0
      real ratio_mod(npoints)
      integer k
! theta_mod = theta * ratio_mod => ratio_mod = mmm0/mmm(p)

      ypklim  = cpp*(  6e6/9.2e6)**0.1914
      ypklim2 = cpp*(8.9e6/9.2e6)**0.1914
      mmm0 = 43.44
      
      DO k = 1, npoints
        ratio_mod(k) = 1.
! ATM PROFONDE DESACTIVEE !
        if (   (1 .EQ. 0)  .AND.(ypk(k).gt.ypklim)) then
	   ratio_mod(k) = mmm0 /                                        &
     &  (mmm0+0.56*(log(ypklim/cpp)-log(ypk(k)/cpp))                    &
     &            /(log(ypklim/cpp)-log(ypklim2/cpp)))
           ratio_mod(k) = max(ratio_mod(k),mmm0/(mmm0+0.56))
	endif
      ENDDO
!----------------------

      if (cpofT) then
          yteta = yt**nu_venus                                          &
     &            - nu_venus * t0_venus**nu_venus * log(ypk/cpp)
          yteta = yteta**(1./nu_venus)
!----------------------
! ATMOSPHERE PROFONDE
          yteta = yteta*ratio_mod
!----------------------

      else
          yteta = yt * cpp/ypk
      endif

      return
      end subroutine t2tpot

!======================================================================
!======================================================================

      SUBROUTINE tpot2t(npoints,yteta, yt, ypk)
!======================================================================
! Arguments:
!
! yteta--------input-R- Temperature potentielle
! yt   -------output-R- Temperature
! ypk  --------input-R- Fonction d'Exner: RCPD*(pplay/pref)**RKAPPA
!
!======================================================================

      USE control_mod, ONLY: cpofT 
      USE comconst_mod, ONLY: cpp,nu_venus,t0_venus

      IMPLICIT NONE

      integer,intent(in) :: npoints
      REAL,intent(in) :: yteta(npoints), ypk(npoints)
      REAL,intent(out) :: yt(npoints)
      
!----------------------
! ATMOSPHERE PROFONDE
      real ypklim,ypklim2,mmm0
      real ratio_mod(npoints)
      integer k
! theta_mod = theta * ratio_mod => ratio_mod = mmm0/mmm(p)

      ypklim  = cpp*(  6e6/9.2e6)**0.1914
      ypklim2 = cpp*(8.9e6/9.2e6)**0.1914
      mmm0 = 43.44

      DO k = 1, npoints
        ratio_mod(k) = 1.
! ATM PROFONDE DESACTIVEE !
        if (   (1 .EQ. 0)  .AND.(ypk(k).gt.ypklim)) then
	   ratio_mod(k) = mmm0 /                                        &
     &  (mmm0+0.56*(log(ypklim/cpp)-log(ypk(k)/cpp))                    &
     &            /(log(ypklim/cpp)-log(ypklim2/cpp)))
           ratio_mod(k) = max(ratio_mod(k),mmm0/(mmm0+0.56))
	endif
      ENDDO
!----------------------

      if (cpofT) then

!----------------------
! ATMOSPHERE PROFONDE
!----------------------
         yt = (yteta/ratio_mod)**nu_venus                               &
!----------------------
     &       + nu_venus * t0_venus**nu_venus * log(ypk/cpp)
         yt = yt**(1./nu_venus)
      else
          yt = yteta * ypk/cpp
      endif
  
      return
      end subroutine tpot2t

!======================================================================
!======================================================================
! Routines pour les calculs paralleles
!======================================================================
#ifdef CPP_PARA
!======================================================================

      SUBROUTINE t2tpot_p(nlon,nlev, yt, yteta, ypk)
! Parallel version of t2tpot, for an arbitrary number of columns
      USE control_mod, only : cpofT
      USE parallel_lmdz, only : OMP_CHUNK
      USE comconst_mod, ONLY: cpp,nu_venus,t0_venus

      IMPLICIT none

      integer,intent(in) :: nlon,nlev
      real,intent(in) :: yt(nlon,nlev)
      real,intent(out) :: yteta(nlon,nlev)
      real,intent(in) :: ypk(nlon,nlev)
! local variable:
      integer :: l
      
!----------------------
! ATMOSPHERE PROFONDE
      real ypklim,ypklim2,mmm0
      real ratio_mod(nlon,nlev)
      integer k
! theta_mod = theta * ratio_mod => ratio_mod = mmm0/mmm(p)

      ypklim  = cpp*(  6e6/9.2e6)**0.1914
      ypklim2 = cpp*(8.9e6/9.2e6)**0.1914
      mmm0 = 43.44

      DO k = 1, nlon
      DO l = 1, nlev
        ratio_mod(k,l) = 1.
! ATM PROFONDE DESACTIVEE !
        if (   (1 .EQ. 0)  .AND.(ypk(k,l).gt.ypklim)) then
	   ratio_mod(k,l) = mmm0 /                                      &
     &  (mmm0+0.56*(log(ypklim/cpp)-log(ypk(k,l)/cpp))                  &
     &            /(log(ypklim/cpp)-log(ypklim2/cpp)))
           ratio_mod(k,l) = max(ratio_mod(k,l),mmm0/(mmm0+0.56))
	endif
      ENDDO
      ENDDO
!----------------------

      if (cpofT) then
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        do l=1,nlev
          yteta(:,l)=yt(:,l)**nu_venus                                  &
     &                     -nu_venus*t0_venus**nu_venus*                &
     &                          log(ypk(:,l)/cpp)
          yteta(:,l)=yteta(:,l)**(1./nu_venus)
!----------------------
! ATMOSPHERE PROFONDE
          yteta(:,l)=yteta(:,l)*ratio_mod(:,l)
!----------------------
        enddo
!$OMP END DO
      else
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        do l=1,nlev
          yteta(:,l)=yt(:,l)*cpp/ypk(:,l)
        enddo
!$OMP END DO
      endif ! of if (cpofT)

      end subroutine t2tpot_p

!======================================================================
!======================================================================

      SUBROUTINE t2tpot_glo_p(yt, yteta, ypk)
! Parallel version of t2tpot, over the full dynamics (scalar) grid
! (more efficient than multiple calls to t2tpot_p() with slices of data)
      USE parallel_lmdz, only : jj_begin,jj_end,OMP_CHUNK
      USE control_mod, only : cpofT
      USE comconst_mod, ONLY: cpp,nu_venus,t0_venus

      IMPLICIT none
! for iip1, jjp1 and llm
#include "dimensions.h"
#include "paramet.h"

      real,intent(in) :: yt(iip1,jjp1,llm)
      real,intent(out) :: yteta(iip1,jjp1,llm)
      real,intent(in) :: ypk(iip1,jjp1,llm)
! local variable:
      integer :: j,l
      integer :: jjb,jje
      
!----------------------
! ATMOSPHERE PROFONDE
      real ypklim,ypklim2,mmm0
      real ratio_mod(iip1,jjp1,llm)
      integer i
! theta_mod = theta * ratio_mod => ratio_mod = mmm0/mmm(p)

      ypklim  = cpp*(  6e6/9.2e6)**0.1914
      ypklim2 = cpp*(8.9e6/9.2e6)**0.1914
      mmm0 = 43.44

      DO i = 1, iip1
      DO j = 1, jjp1
      DO l = 1, llm
        ratio_mod(i,j,l) = 1.
! ATM PROFONDE DESACTIVEE !
        if (   (1 .EQ. 0)  .AND.(ypk(i,j,l).gt.ypklim)) then
	   ratio_mod(i,j,l) = mmm0 /                                    &
     &  (mmm0+0.56*(log(ypklim/cpp)-log(ypk(i,j,l)/cpp))                &
     &            /(log(ypklim/cpp)-log(ypklim2/cpp)))
           ratio_mod(i,j,l) = max(ratio_mod(i,j,l),mmm0/(mmm0+0.56))
	endif
      ENDDO
      ENDDO
      ENDDO
!----------------------

      jjb=jj_begin
      jje=jj_end

      if (cpofT) then
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        do l=1,llm
          yteta(:,jjb:jje,l)=yt(:,jjb:jje,l)**nu_venus                  &
     &                     -nu_venus*t0_venus**nu_venus*                &
     &                          log(ypk(:,jjb:jje,l)/cpp)
          yteta(:,jjb:jje,l)=yteta(:,jjb:jje,l)**(1./nu_venus)
!----------------------
! ATMOSPHERE PROFONDE
          yteta(:,jjb:jje,l)=yteta(:,jjb:jje,l)*ratio_mod(:,jjb:jje,l)
!----------------------
        enddo
!$OMP END DO
      else
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        do l=1,llm
          yteta(:,jjb:jje,l)=yt(:,jjb:jje,l)*cpp/ypk(:,jjb:jje,l)
        enddo
!$OMP END DO
      endif ! of if (cpofT)

      end subroutine t2tpot_glo_p

!======================================================================
!======================================================================

      SUBROUTINE tpot2t_p(nlon,nlev,yteta,yt,ypk)
! Parallel version of tpot2t, for an arbitrary number of columns
      USE control_mod, only : cpofT
      USE parallel_lmdz, only : OMP_CHUNK
      USE comconst_mod, ONLY: cpp,nu_venus,t0_venus

      IMPLICIT none

      integer,intent(in) :: nlon,nlev
      real,intent(out) :: yt(nlon,nlev)
      real,intent(in) :: yteta(nlon,nlev)
      real,intent(in) :: ypk(nlon,nlev)

! local variable:
      integer :: l

!----------------------
! ATMOSPHERE PROFONDE
      real ypklim,ypklim2,mmm0
      real ratio_mod(nlon,nlev)
      integer k
! theta_mod = theta * ratio_mod => ratio_mod = mmm0/mmm(p)

      ypklim  = cpp*(  6e6/9.2e6)**0.1914
      ypklim2 = cpp*(8.9e6/9.2e6)**0.1914
      mmm0 = 43.44

      DO k = 1, nlon
      DO l = 1, nlev
        ratio_mod(k,l) = 1.
! ATM PROFONDE DESACTIVEE !
        if (   (1 .EQ. 0)  .AND.(ypk(k,l).gt.ypklim)) then
	   ratio_mod(k,l) = mmm0 /                                      &
     &  (mmm0+0.56*(log(ypklim/cpp)-log(ypk(k,l)/cpp))                  &
     &            /(log(ypklim/cpp)-log(ypklim2/cpp)))
           ratio_mod(k,l) = max(ratio_mod(k,l),mmm0/(mmm0+0.56))
	endif
      ENDDO
      ENDDO
!----------------------

      if (cpofT) then
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        do l=1,nlev
!----------------------
! ATMOSPHERE PROFONDE
          yt(:,l)=(yteta(:,l)/ratio_mod(:,l))**nu_venus                 &
!----------------------
     &                  +nu_venus*t0_venus**nu_venus*                   &
     &                       log(ypk(:,l)/cpp)
          yt(:,l)=yt(:,l)**(1./nu_venus)
        enddo
!$OMP END DO
      else
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        do l=1,nlev
          yt(:,l)=yteta(:,l)*ypk(:,l)/cpp
        enddo
!$OMP END DO
      endif ! of if (cpofT)
      end subroutine tpot2t_p

!======================================================================
!======================================================================

      SUBROUTINE tpot2t_glo_p(yteta,yt,ypk)
! Parallel version of tpot2t, over the full dynamics (scalar) grid
! (more efficient than multiple calls to tpot2t_p() with slices of data)
      USE parallel_lmdz, only : jj_begin,jj_end,OMP_CHUNK
      USE control_mod, only : cpofT
      USE comconst_mod, ONLY: cpp,nu_venus,t0_venus

      IMPLICIT none
! for iip1, jjp1 and llm
#include "dimensions.h"
#include "paramet.h"

      real,intent(out) :: yt(iip1,jjp1,llm)
      real,intent(in) :: yteta(iip1,jjp1,llm)
      real,intent(in) :: ypk(iip1,jjp1,llm)
! local variable:
      integer :: j,l
      integer :: jjb,jje
      
!----------------------
! ATMOSPHERE PROFONDE
      real ypklim,ypklim2,mmm0
      real ratio_mod(iip1,jjp1,llm)
      integer i
! theta_mod = theta * ratio_mod => ratio_mod = mmm0/mmm(p)

      ypklim  = cpp*(  6e6/9.2e6)**0.1914
      ypklim2 = cpp*(8.9e6/9.2e6)**0.1914
      mmm0 = 43.44

      DO i = 1, iip1
      DO j = 1, jjp1
      DO l = 1, llm
        ratio_mod(i,j,l) = 1.
! ATM PROFONDE DESACTIVEE !
        if (   (1 .EQ. 0)  .AND.(ypk(i,j,l).gt.ypklim)) then
	   ratio_mod(i,j,l) = mmm0 /                                    &
     &  (mmm0+0.56*(log(ypklim/cpp)-log(ypk(i,j,l)/cpp))                &
     &            /(log(ypklim/cpp)-log(ypklim2/cpp)))
           ratio_mod(i,j,l) = max(ratio_mod(i,j,l),mmm0/(mmm0+0.56))
	endif
      ENDDO
      ENDDO
      ENDDO
!----------------------

      jjb=jj_begin
      jje=jj_end

      if (cpofT) then
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        do l=1,llm
!----------------------
! ATMOSPHERE PROFONDE
          yt(:,jjb:jje,l)=                                              &
     &          (yteta(:,jjb:jje,l)/ratio_mod(:,jjb:jje,l))**nu_venus   &
!----------------------
     &                  +nu_venus*t0_venus**nu_venus*                   &
     &                       log(ypk(:,jjb:jje,l)/cpp)
          yt(:,jjb:jje,l)=yt(:,jjb:jje,l)**(1./nu_venus)
        enddo
!$OMP END DO
      else
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        do l=1,llm
          yt(:,jjb:jje,l)=yteta(:,jjb:jje,l)*ypk(:,jjb:jje,l)/cpp
        enddo
!$OMP END DO
      endif ! of if (cpofT)
      end subroutine tpot2t_glo_p

!======================================================================
#endif
!======================================================================
! Fin routines specifiques parallele
!======================================================================
!======================================================================
!
! ATTENTION
!
! Si un jour on a besoin, il faudra coder les routines 
!    dt2dtpot / dtpto2dt 
!
!======================================================================
!======================================================================
end module cpdet_mod
