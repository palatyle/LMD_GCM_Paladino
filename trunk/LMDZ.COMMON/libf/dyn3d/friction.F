!
! $Id: friction.F 1437 2010-09-30 08:29:10Z emillour $
!
c=======================================================================
      SUBROUTINE friction(ucov,vcov,pdt)

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
      REAL,INTENT(out) :: ucov( iip1,jjp1,llm )
      REAL,INTENT(out) :: vcov( iip1,jjm,llm )
      REAL,INTENT(in) :: pdt ! time step

! local variables:

      REAL modv(iip1,jjp1),zco,zsi
      REAL vpn,vps,upoln,upols,vpols,vpoln
      REAL u2(iip1,jjp1),v2(iip1,jjm)
      INTEGER  i,j,l
      REAL,PARAMETER :: cfric=1.e-5
      LOGICAL,SAVE :: firstcall=.true.
      INTEGER,SAVE :: friction_type=1
      CHARACTER(len=20) :: modname="friction"
      CHARACTER(len=80) :: abort_message
      
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

      if (friction_type.eq.0) then
c   calcul des composantes au carre du vent naturel
      do j=1,jjp1
         do i=1,iip1
            u2(i,j)=ucov(i,j,1)*ucov(i,j,1)*unscu2(i,j)
         enddo
      enddo
      do j=1,jjm
         do i=1,iip1
            v2(i,j)=vcov(i,j,1)*vcov(i,j,1)*unscv2(i,j)
         enddo
      enddo

c   calcul du module de V en dehors des poles
      do j=2,jjm
         do i=2,iip1
            modv(i,j)=sqrt(0.5*(u2(i-1,j)+u2(i,j)+v2(i,j-1)+v2(i,j)))
         enddo
         modv(1,j)=modv(iip1,j)
      enddo

c   les deux composantes du vent au pole sont obtenues comme
c   premiers modes de fourier de v pres du pole
      upoln=0.
      vpoln=0.
      upols=0.
      vpols=0.
      do i=2,iip1
         zco=cos(rlonv(i))*(rlonu(i)-rlonu(i-1))
         zsi=sin(rlonv(i))*(rlonu(i)-rlonu(i-1))
         vpn=vcov(i,1,1)/cv(i,1)
         vps=vcov(i,jjm,1)/cv(i,jjm)
         upoln=upoln+zco*vpn
         vpoln=vpoln+zsi*vpn
         upols=upols+zco*vps
         vpols=vpols+zsi*vps
      enddo
      vpn=sqrt(upoln*upoln+vpoln*vpoln)/pi
      vps=sqrt(upols*upols+vpols*vpols)/pi
      do i=1,iip1
c        modv(i,1)=vpn
c        modv(i,jjp1)=vps
         modv(i,1)=modv(i,2)
         modv(i,jjp1)=modv(i,jjm)
      enddo

c   calcul du frottement au sol.
      do j=2,jjm
         do i=1,iim
            ucov(i,j,1)=ucov(i,j,1)
     s      -cfric*pdt*0.5*(modv(i+1,j)+modv(i,j))*ucov(i,j,1)
         enddo
         ucov(iip1,j,1)=ucov(1,j,1)
      enddo
      do j=1,jjm
         do i=1,iip1
            vcov(i,j,1)=vcov(i,j,1)
     s      -cfric*pdt*0.5*(modv(i,j+1)+modv(i,j))*vcov(i,j,1)
         enddo
         vcov(iip1,j,1)=vcov(1,j,1)
      enddo
      endif ! of if (friction_type.eq.0)

      if (friction_type.eq.1) then
        do l=1,llm
          ucov(:,:,l)=ucov(:,:,l)*(1.-pdt*kfrict(l))
          vcov(:,:,l)=vcov(:,:,l)*(1.-pdt*kfrict(l))
        enddo
      endif
      
      RETURN
      END

