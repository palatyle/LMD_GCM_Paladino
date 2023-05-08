!Mars importe directement de la version martienne 
      subroutine groupeun(jjmax,llmax,q)
      implicit none

c *********************************************************
c     lissage fort des champs aux poles pour permettre 
c     de tourner avec plafond haut.
c     F. Hourdin, 1997.
c     (commentaire et test sur ngroup par Francois Forget, 05/2000)
c
c     Technique :
c
c     1) Les points sont regroupes par paquet de 2**ngroup 
c     aux poles (e.g. ngroup=3 -> paquet de 8), puis 2**(ngroup-1)
c     2**(ngroup-2), etc... aux latitude plus basse en s'eloignant du pole
c
c     2) Plus ngroup est eleve, plus le lissage est efficace
c
c     3) MAIS, il faut iim divisible par 2**ngroup !!!
c
c
c *********************************************************

#include "dimensions.h"
#include "paramet.h"
#include "comconst.h"
#include "comgeom2.h"

      integer jjmax,llmax
      real q(iip1,jjmax,llmax)

      integer ngroup
      parameter (ngroup=2) 

      logical intensive

      real airen,airecn,qn
      real aires,airecs,qs

      integer i,j,l,ig,j1,j2,i0,jd
      logical firstcall
      DATA firstcall/.true./

c -------------------------------------------------------
c   test sur ngroup :

      if (firstcall) then
        IF(MOD(iim,2**ngroup).NE.0) then
           write(*,*) 'Problem in groupeun.F'
           write(*,*) 'iim= ', iim, ' ngroup=', ngroup
           write(*,*) 'iim must be divisible by par 2**ngroup= ',
     &                 2**ngroup
           write(*,*) ' YOU must change ngroup in groupeun.F'
           write(*,*) ' Have fun... '
           stop
        END IF
        write(*,*) 'in groupeun.F, ngroup=', ngroup
        firstcall = .false.
      end if
      

c Champs 3D
      intensive=.true.
      if (intensive) then
      jd=jjp1-jjmax
      do l=1,llmax
       j1=1+jd
       j2=2
       do ig=1,ngroup
         do j=j1-jd,j2-jd
c           print*,'groupe ',ig,'  j= ',j,2**(ngroup-ig+1),'pts groupes'
            do i0=1,iim,2**(ngroup-ig+1)
               airen=0.
               airecn=0.
               qn=0.
               aires=0.
               airecs=0.
               qs=0.
               do i=i0,i0+2**(ngroup-ig+1)-1
                  airen=airen+aire(i,j)
                  aires=aires+aire(i,jjp1-j+1)
                  qn=qn+q(i,j,l)*aire(i,j)
                  qs=qs+q(i,jjp1-j+1-jd,l)*aire(i,jjp1-j+1)
               enddo
               airecn=0.
               airecs=0.
               do i=i0,i0+2**(ngroup-ig+1)-1
                  q(i,j,l)=qn/airen
                  q(i,jjp1-j+1-jd,l)=qs/aires
               enddo
            enddo
            q(iip1,j,l)=q(1,j,l)
            q(iip1,jjp1-j+1-jd,l)=q(1,jjp1-j+1-jd,l)
         enddo
         j1=j2+1
         j2=j2+2**ig
       enddo
      enddo

c  ---------------
      else
c   Cas sans multiplication par les aires.
c  ---------------

      jd=jjp1-jjmax
      do l=1,llmax
      j1=1+jd
      j2=2
      do ig=1,ngroup
         do j=j1-jd,j2-jd
c           print*,'groupe ',ig,'  j= ',j,2**(ngroup-ig+1),'pts groupes'
            do i0=1,iim,2**(ngroup-ig+1)
               airen=0.
               airecn=0.
               qn=0.
               aires=0.
               airecs=0.
               qs=0.
               do i=i0,i0+2**(ngroup-ig+1)-1
                  airen=airen+aire(i,j)
                  aires=aires+aire(i,jjp1-j+1)
                  qn=qn+q(i,j,l)
                  qs=qs+q(i,jjp1-j+1-jd,l)
               enddo
               airecn=0.
               airecs=0.
               do i=i0,i0+2**(ngroup-ig+1)-1
                  q(i,j,l)=qn*aire(i,j)/airen
                  q(i,jjp1-j+1-jd,l)=qs*aire(i,jjp1-j+1)/aires
               enddo
            enddo
            q(iip1,j,l)=q(1,j,l)
            q(iip1,jjp1-j+1-jd,l)=q(1,jjp1-j+1-jd,l)
         enddo
         j1=j2+1
         j2=j2+2**ig
      enddo
      enddo

      endif

c     print*,'On groupe jusque j=',j1-1,rlatu(j1-1)*180./pi

      return
      end
