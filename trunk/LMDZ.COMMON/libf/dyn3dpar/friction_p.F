!
! $Id: friction_p.F 1437 2010-09-30 08:29:10Z emillour $
!
c=======================================================================
      SUBROUTINE friction_p(ucov,vcov,pdt)
      USE parallel_lmdz
      USE control_mod
#ifdef CPP_IOIPSL
      USE IOIPSL
#else
! if not using IOIPSL, we still need to use (a local version of) getin
      USE ioipsl_getincom
#endif
      USE comconst_mod, ONLY: pi
      IMPLICIT NONE

!=======================================================================
!
!   Friction for the Newtonian case:
!   --------------------------------
!    2 possibilities (depending on flag 'friction_type'
!     friction_type=0 : A friction that is only applied to the lowermost
!                       atmospheric layer
!     friction_type=1 : Friction applied on all atmospheric layer (but
!       (default)       with stronger magnitude near the surface; see
!                       iniacademic.F)
!=======================================================================

#include "dimensions.h"
#include "paramet.h"
#include "comgeom2.h"
#include "iniprint.h"
#include "academic.h"

! arguments:
      REAL,INTENT(inout) :: ucov( iip1,jjp1,llm )
      REAL,INTENT(inout) :: vcov( iip1,jjm,llm )
      REAL,INTENT(in) :: pdt ! time step

! local variables:
      REAL modv(iip1,jjp1),zco,zsi
      REAL vpn,vps,upoln,upols,vpols,vpoln
      REAL u2(iip1,jjp1),v2(iip1,jjm)
      INTEGER  i,j,l
      REAL,PARAMETER :: cfric=1.e-5
      LOGICAL,SAVE :: firstcall=.true.
      INTEGER,SAVE :: friction_type=1
      CHARACTER(len=20) :: modname="friction_p"
      CHARACTER(len=80) :: abort_message
!$OMP THREADPRIVATE(firstcall,friction_type)
      integer :: jjb,jje

!$OMP SINGLE
      IF (firstcall) THEN
        ! set friction type
        call getin("friction_type",friction_type)
        if ((friction_type.lt.0).or.(friction_type.gt.1)) then
          abort_message="wrong friction type"
          write(lunout,*)'Friction: wrong friction type',friction_type
          call abort_gcm(modname,abort_message,42)
        endif
        firstcall=.false.
      ENDIF
!$OMP END SINGLE COPYPRIVATE(friction_type,firstcall)

      if (friction_type.eq.0) then ! friction on first layer only
!$OMP SINGLE
c   calcul des composantes au carre du vent naturel
      jjb=jj_begin
      jje=jj_end+1
      if (pole_sud) jje=jj_end
      
      do j=jjb,jje
         do i=1,iip1
            u2(i,j)=ucov(i,j,1)*ucov(i,j,1)*unscu2(i,j)
         enddo
      enddo
      
      jjb=jj_begin-1
      jje=jj_end+1
      if (pole_nord) jjb=jj_begin
      if (pole_sud) jje=jj_end-1
      
      do j=jjb,jje
         do i=1,iip1
            v2(i,j)=vcov(i,j,1)*vcov(i,j,1)*unscv2(i,j)
         enddo
      enddo

c   calcul du module de V en dehors des poles
      jjb=jj_begin
      jje=jj_end+1
      if (pole_nord) jjb=jj_begin+1
      if (pole_sud) jje=jj_end-1
      
      do j=jjb,jje
         do i=2,iip1
            modv(i,j)=sqrt(0.5*(u2(i-1,j)+u2(i,j)+v2(i,j-1)+v2(i,j)))
         enddo
         modv(1,j)=modv(iip1,j)
      enddo

c   les deux composantes du vent au pole sont obtenues comme
c   premiers modes de fourier de v pres du pole
      if (pole_nord) then
      
        upoln=0.
        vpoln=0.
     
        do i=2,iip1
           zco=cos(rlonv(i))*(rlonu(i)-rlonu(i-1))
           zsi=sin(rlonv(i))*(rlonu(i)-rlonu(i-1))
           vpn=vcov(i,1,1)/cv(i,1)
           upoln=upoln+zco*vpn
           vpoln=vpoln+zsi*vpn
        enddo
        vpn=sqrt(upoln*upoln+vpoln*vpoln)/pi
        do i=1,iip1
c          modv(i,1)=vpn
           modv(i,1)=modv(i,2)
        enddo

      endif
      
      if (pole_sud) then
      
        upols=0.
        vpols=0.
        do i=2,iip1
           zco=cos(rlonv(i))*(rlonu(i)-rlonu(i-1))
           zsi=sin(rlonv(i))*(rlonu(i)-rlonu(i-1))
           vps=vcov(i,jjm,1)/cv(i,jjm)
           upols=upols+zco*vps
           vpols=vpols+zsi*vps
        enddo
        vps=sqrt(upols*upols+vpols*vpols)/pi
        do i=1,iip1
c        modv(i,jjp1)=vps
         modv(i,jjp1)=modv(i,jjm)
        enddo
      
      endif
      
c   calcul du frottement au sol.

      jjb=jj_begin
      jje=jj_end
      if (pole_nord) jjb=jj_begin+1
      if (pole_sud) jje=jj_end-1

      do j=jjb,jje
         do i=1,iim
            ucov(i,j,1)=ucov(i,j,1)
     s      -cfric*pdt*0.5*(modv(i+1,j)+modv(i,j))*ucov(i,j,1)
         enddo
         ucov(iip1,j,1)=ucov(1,j,1)
      enddo
      
      jjb=jj_begin
      jje=jj_end
      if (pole_sud) jje=jj_end-1
      
      do j=jjb,jje
         do i=1,iip1
            vcov(i,j,1)=vcov(i,j,1)
     s      -cfric*pdt*0.5*(modv(i,j+1)+modv(i,j))*vcov(i,j,1)
         enddo
         vcov(iip1,j,1)=vcov(1,j,1)
      enddo
!$OMP END SINGLE
      endif ! of if (friction_type.eq.0)

      if (friction_type.eq.1) then
       ! for ucov() 
        jjb=jj_begin
        jje=jj_end
        if (pole_nord) jjb=jj_begin+1
        if (pole_sud) jje=jj_end-1

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        do l=1,llm
          ucov(1:iip1,jjb:jje,l)=ucov(1:iip1,jjb:jje,l)*
     &                                  (1.-pdt*kfrict(l))
        enddo
!$OMP END DO NOWAIT
        
       ! for vcoc()
        jjb=jj_begin
        jje=jj_end
        if (pole_sud) jje=jj_end-1
        
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        do l=1,llm
          vcov(1:iip1,jjb:jje,l)=vcov(1:iip1,jjb:jje,l)*
     &                                  (1.-pdt*kfrict(l))
        enddo
!$OMP END DO
      endif ! of if (friction_type.eq.1)

      RETURN
      END

