      subroutine sponge(ucov,vcov,h,pext,dt,mode)
      implicit none
#include "dimensions.h"
#include "paramet.h"
#include "comdissip.h"
#include "comvert.h"
#include "comgeom2.h"
#include "sponge.h"
      real ucov(iip1,jjp1,llm),vcov(iip1,jjm,llm)
      real h(iip1,jjp1,llm),pext(iip1,jjp1)
      real dt
      real sig_s(llm) !sigma au milieu des couches
      save sig_s

c   Local:
c   ------

      REAL vm,um,hm,ptot(jjp1)
      real cst(llm)
      integer mode

      INTEGER l,i,j,l0

      real ssum

      real echelle,zkm
      logical firstcall
      save firstcall,cst,l0
      data firstcall/.true./

      if (firstcall) then
       do l=1,llm
          sig_s(l)=((ap(l)+ap(l+1))/preff+bp(l)+bp(l+1))/2.
       enddo

       IF (mode.eq.3) then
         l0=llm-2
         echelle=10.
 
         PRINT*
         print*,'sponge mode',mode
         print*,'hsponge',hsponge
         print*,'tetasponge n intervient pas'
         print*,'Constantes de dissipations fixees comme les anglais'
         print*,'Coeffs pour la couche en eponge'
         print*,'Z (km)  tau'
         cst(llm)=dt/(0.5*88775)
         cst(llm-1)=dt/(88775)
         cst(l0)=dt/(2*88775)
         do l=l0,llm
            zkm=-echelle*log(sig_s(l))
            print*,zkm,dt/cst(l),cst(l)
         enddo
         firstcall=.false.
         PRINT*
       ELSE
         l0=1
         echelle=10.
 
         PRINT*
         print*,'sponge mode',mode
         print*,'hsponge tetasponge ',hsponge,tetasponge
         print*,'Coeffs pour la couche en eponge'
         print*,'Z (km)  tau'
         do l=1,llm
            zkm=-echelle*log(sig_s(l))
            cst(l)=.5*(1.+tanh((zkm-hsponge)/echelle))
            cst(l)= max(tetasponge*1.e-15,cst(l))
            print*,zkm,1./cst(l)*tetasponge,cst(l)*dt/tetasponge
            cst(l)=cst(l)*dt/tetasponge
         enddo
         firstcall=.false.
         PRINT*
       ENDIF
      endif

c-----------------------------------------------------------------------
c   calcul de la dissipation:
c   -------------------------

      do j=1,jjp1
         ptot(j)=ssum(iim,pext(1,j),1)
      enddo
 
c   temperature potentielle
      do l=l0,llm
         do j=1,jjp1
            hm=0.
            do i=1,iim
               hm=hm+h(i,j,l)*pext(i,j)
            enddo
            hm=hm/ptot(j)
            do i=1,iim
               h(i,j,l)=h(i,j,l)-cst(l)*(h(i,j,l)-hm)
            enddo
            h(iip1,j,l)=h(1,j,l)
         enddo
      enddo

c   vent zonal
      do l=l0,llm
         do j=2,jjm
            um=0.
            if(mode.ge.1) then
               do i=1,iim
                  um=um+0.5*ucov(i,j,l)*(pext(i,j)+pext(i+1,j))
     s               /cu(i,j)
               enddo
               um=um/ptot(j)
            endif
            do i=1,iim
               ucov(i,j,l)=ucov(i,j,l)-cst(l)*(ucov(i,j,l)-um*cu(i,j))
            enddo
            ucov(iip1,j,l)=ucov(1,j,l)
         enddo
      enddo

c  vent meridien
      do l=l0,llm
         do j=1,jjm
            vm=0.
            if(mode.ge.2) then
               do i=1,iim
                  vm=vm+vcov(i,j,l)*(pext(i,j)+pext(i,j+1))
     s               /cv(i,j)
               enddo
               vm=vm/(ptot(j)+ptot(j+1))
            endif
            do i=1,iim
               vcov(i,j,l)=vcov(i,j,l)-cst(l)*(vcov(i,j,l)-vm*cv(i,j))
            enddo
            vcov(iip1,j,l)=vcov(1,j,l)
         enddo
      enddo

      RETURN
      end
