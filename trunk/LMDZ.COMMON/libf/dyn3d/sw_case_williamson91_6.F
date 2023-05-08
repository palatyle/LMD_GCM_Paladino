!
! $Id $
!
      SUBROUTINE sw_case_williamson91_6(vcov,ucov,teta,masse,ps)

c=======================================================================
c
c   Author:    Thomas Dubos      original: 26/01/2010
c   -------
c
c   Subject:
c   ------
c   Realise le cas-test 6 de Williamson et al. (1991) : onde de Rossby-Haurwitz
c
c   Method:
c   --------
c
c   Interface:
c   ----------
c
c      Input:
c      ------
c
c      Output:
c      -------
c
c=======================================================================
      USE comvert_mod, ONLY: ap,bp,preff
      USE comconst_mod, ONLY: cpp,omeg,rad

      IMPLICIT NONE
c-----------------------------------------------------------------------
c   Declararations:
c   ---------------

#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
#include "iniprint.h"

c   Arguments:
c   ----------

c   variables dynamiques
      REAL vcov(ip1jm,llm),ucov(ip1jmp1,llm) ! vents covariants
      REAL teta(ip1jmp1,llm)                 ! temperature potentielle
      REAL ps(ip1jmp1)                       ! pression  au sol
      REAL masse(ip1jmp1,llm)                ! masse d'air
      REAL phis(ip1jmp1)                     ! geopotentiel au sol

c   Local:
c   ------

      REAL p (ip1jmp1,llmp1  )               ! pression aux interfac.des couches
      REAL pks(ip1jmp1)                      ! exner au  sol
      REAL pk(ip1jmp1,llm)                   ! exner au milieu des couches
      REAL pkf(ip1jmp1,llm)                  ! exner filt.au milieu des couches
      REAL alpha(ip1jmp1,llm),beta(ip1jmp1,llm)

      REAL :: sinth,costh,costh2, Ath,Bth,Cth, lon,dps
      INTEGER i,j,ij

      REAL, PARAMETER    :: rho=1 ! masse volumique de l'air (arbitraire)
      REAL, PARAMETER    :: K    = 7.848e-6  ! K = \omega
      REAL, PARAMETER    :: gh0  = 9.80616 * 8e3 
      INTEGER, PARAMETER :: R0=4, R1=R0+1, R2=R0+2         ! mode 4
c NB : rad = 6371220 dans W91 (6371229 dans LMDZ)
c      omeg = 7.292e-5 dans W91 (7.2722e-5 dans LMDZ)
 
      IF(0==0) THEN
c Williamson et al. (1991) : onde de Rossby-Haurwitz
         teta = preff/rho/cpp
c geopotentiel (pression de surface)
         do j=1,jjp1
            costh2 = cos(rlatu(j))**2
            Ath = (R0+1)*(costh2**2) + (2*R0*R0-R0-2)*costh2 - 2*R0*R0
            Ath = .25*(K**2)*(costh2**(R0-1))*Ath
            Ath = .5*K*(2*omeg+K)*costh2 + Ath 
            Bth = (R1*R1+1)-R1*R1*costh2
            Bth = 2*(omeg+K)*K/(R1*R2) * (costh2**(R0/2))*Bth
            Cth = R1*costh2 - R2
            Cth = .25*K*K*(costh2**R0)*Cth
            do i=1,iip1
               ij=(j-1)*iip1+i
               lon = rlonv(i)
               dps = Ath + Bth*cos(R0*lon) + Cth*cos(2*R0*lon)
               ps(ij) = rho*(gh0 + (rad**2)*dps)
            enddo
         enddo
         write(lunout,*) 'W91 ps', MAXVAL(ps), MINVAL(ps)
c vitesse zonale ucov
         do j=1,jjp1
            costh  = cos(rlatu(j))
            costh2 = costh**2
            Ath = rad*K*costh
            Bth = R0*(1-costh2)-costh2
            Bth = rad*K*Bth*(costh**(R0-1))
            do i=1,iip1
               ij=(j-1)*iip1+i
               lon = rlonu(i)
               ucov(ij,1) = (Ath + Bth*cos(R0*lon))
            enddo
         enddo
         write(lunout,*) 'W91 u', MAXVAL(ucov(:,1)), MINVAL(ucov(:,1))
         ucov(:,1)=ucov(:,1)*cu
c vitesse meridienne vcov
         do j=1,jjm
            sinth  = sin(rlatv(j))
            costh  = cos(rlatv(j))
            Ath = -rad*K*R0*sinth*(costh**(R0-1))
            do i=1,iip1
               ij=(j-1)*iip1+i
               lon = rlonv(i)
               vcov(ij,1) = Ath*sin(R0*lon)
            enddo
         enddo
         write(lunout,*) 'W91 v', MAXVAL(vcov(:,1)), MINVAL(vcov(:,1))
         vcov(:,1)=vcov(:,1)*cv
         
c         ucov=0
c         vcov=0
      ELSE
c test non-tournant, onde se propageant en latitude
         do j=1,jjp1
            do i=1,iip1
               ij=(j-1)*iip1+i
               ps(ij) = 1e5*(1 + .1*exp(-100*(1+sin(rlatu(j)))**2) )
            enddo
         enddo
         
c     rho = preff/(cpp*teta)
         teta = .01*preff/cpp   ! rho = 100 ; phi = ps/rho = 1e3 ; c=30 m/s = 2600 km/j = 23 degres / j
         ucov=0.
         vcov=0.
      END IF      
      
      CALL pression ( ip1jmp1, ap, bp, ps, p       )
      CALL massdair(p,masse)

      END
c-----------------------------------------------------------------------
