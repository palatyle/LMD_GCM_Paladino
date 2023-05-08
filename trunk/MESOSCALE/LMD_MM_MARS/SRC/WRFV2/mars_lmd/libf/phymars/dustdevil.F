      SUBROUTINE dustdevil(ngrid,nlay,nq, pplev,pu,pv,pt, ptsurf,pq2,
     &                pdqdev,pdqs_dev)
      IMPLICIT NONE

c=======================================================================
c
c
c  VERSION SPECIAL TRACEURS :
c  Parameterization of dust devil activities
c  to estimate dust lifting
c  The dust devil activity is estimated based on 
c  Renno et al. 1998 (JAS 55, 3244-3252)  
c
c  It is proportional to (1-b)*Fs
c
c  With b= [ps**(rcp+1) - ptop**(rcp+1)] / [(ps-ptop)*(rcp+1)* ps**rcp]
c  with ptop pressure of the top of the boundary layer
c       rcp = R/cp
c
c  and Fs the surface sensible heat flux = Cd*|U|*(T(1) -Tsurf)
c       
c=======================================================================

c-----------------------------------------------------------------------
c   declarations:
c   -------------

#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"      
c#include "comconst.h"        ! TEMPORAIRE AVEC ANLDEVIL !!!!
#include "planete.h"
#include "comgeomfi.h"
#include "tracer.h"
c   arguments:
c   ----------

      INTEGER ngrid,nlay
      REAL pplev(ngrid,nlay+1)
      REAL pt(ngrid,nlay)
      REAL pu(ngrid,nlay)
      REAL pv(ngrid,nlay)
      REAL pq2(ngrid,nlay+1)
      REAL ptsurf(ngrid)

c    Traceurs :
      integer nq 
      real pdqdev(ngrid,nlay,nq) 
      real pdqs_dev(ngrid,nq) 
      
c   local:
c   ------

      INTEGER ig,l,iq
      real Cd, z1
      save Cd

      LOGICAL firstcall
      SAVE firstcall


      REAL devila(ngridmx)
      integer ltop(ngridmx)
      real b,rho,Fs,wind



      REAL  q2top , seuil
      SAVE  q2top , seuil
      DATA q2top/.5/ ! value of q2 at the top of PBL
      DATA seuil/.3/ ! value of minimum dust devil activity for dust lifting


      DATA firstcall/.true./

c   TEMPORAIRE AVEC ANLDEVIL : *************
c        real b_diag(ngridmx)
c       real localtime(ngridmx)
c       common/temporaire/localtime
c      real ztop(ngridmx),magwind(ngridmx),t1(ngridmx)
c      real rcp ,cpp
c      rcp = kappa
c      cpp = r/rcp
c   **********************************
       

c-----------------------------------------------------------------------
c    initialisation
c    --------------

      IF (firstcall) THEN

        write(*,*) 'In dustdevil :'
        write(*,*) '    q2top= ',q2top,'     seuil= ', seuil 
c       un petit test de coherence:
         IF(ngrid.NE.ngridmx) THEN
            PRINT*,'STOP dans coefdifv'
            PRINT*,'probleme de dimensions :'
            PRINT*,'ngrid  =',ngrid
            PRINT*,'ngridmx  =',ngridmx
            STOP
         ENDIF

c A rough estimation of the horizontal drag coefficient Cd
c (the scale heigh is taken to be 13 km since we are typically
c dealing with daytime temperature around 250K.
c 
         z1= -0.5*13.e3*log(pplev(1,2)/pplev(1,1))
         Cd = (0.4/log(z1/z0))**2

         firstcall=.false.

c        Temporaire
c        open(77,file='devil')
      
      ENDIF

c-----------------------------------------------------------------------
c Initialisation
      do iq=1,nq
       do l=1,nlay
           do ig=1,ngrid
             pdqdev(ig,l,iq)= 0
           end do
       end do
      end do


c-----------------------------------------------------------------------
c      Determining the top of the boundary layer
c      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do ig=1,ngrid
         do  l=2,nlay-1
            if (pq2(ig,l).lt.q2top)then
              ltop(ig)=l
              goto 99
            end if
         enddo
 99      continue

c        ***************************************
cc        if (ptsurf(ig).gt.255)then
c         write(*,*) 'tsurf, ztop (km): ', ig,
c     &   ptsurf(ig), -12.*log(pplev(ig,ltop(ig))/pplev(ig,1))
c         if ((ptsurf(ig).gt.50.).and.(
c     &      (-12.*log(pplev(ig,ltop(ig))/pplev(ig,1))).gt.60.))then
c            do l=1,nlay
c             write(*,*) l,
c     &       -12.*log(pplev(ig,l)/pplev(ig,1)),pq2(ig,l)
c            end do
c         end if
cc        end if
c        ***************************************
     
      enddo

c        ***************************************
c        do ig=100,148
c           write(*,*)'time,tsurf,ztop', localtime(ig),ptsurf(ig),
c    &      -12.*log(pplev(ig,ltop(ig))/pplev(ig,1))
c        end do
c        ***************************************


c   Calculation : dust devil intensity
c   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do ig=1,ngrid

c --------------------------------------------------
c 1) Version 1 : sensible heat flux using actual wind :
c        Wind magnitude:
c        wind = sqrt(pu(ig,1)**2+pv(ig,1)**2)
c
c --------------------------------------------------
c 2) Version 2 : sensible heat flux using  wind = 15 m/s
         wind = 15.
c ----------------------------------------------------
c        Density :
         rho=pplev(ig,1)/(R*pt(ig,1))

c        Sensible heat flux (W.m-2) (>0 if up)
         Fs= rho*cpp*Cd * wind
     &       * (ptsurf(ig) -pt(ig,1))
         b= (pplev(ig,1)**(rcp+1) - pplev(ig,ltop(ig))**(rcp+1)) /
     &    ( (pplev(ig,1)-pplev(ig,ltop(ig)))*(rcp+1)*pplev(ig,1)**rcp)

c        b_diag(ig) = b     ! TEMPORAIRE (pour diagnostique)

c   Energy flux available to drive dust devil (W.m-2) : (1-b)*Fs
c   Dust devil activity : 
         devila(ig)= max( 0. , (1-b)*Fs - seuil ) 
      enddo
c   
c     lifted dust (kg m-2 s-1)  (<0 when lifting)
c     ~~~~~~~~~~  
      do iq=1,nq
         do ig=1,ngrid
           pdqs_dev(ig,iq)= - alpha_devil(iq) * devila(ig)
         end do
      end do 

c     Injection of dust in the atmosphere (up to the top of pbl)
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do iq=1,nq
       do ig=1,ngrid
         if (devila(ig).ne.0.) then
           do l=1,ltop(ig)
             pdqdev(ig,l,iq)=-pdqs_dev(ig,iq)*g/
     &        (pplev(ig,1)-pplev(ig,ltop(ig)))
           end do
         end if
       end do
      end do

c *********************************************************     
c     TEMPORAIRE AVEC ANLDEVIL:
c     IF (ngrid.gt.1) THEN
c      do ig=2,ngridmx-1
c       write(77,88) lati(ig)*180./pi,localtime(ig),
c    &        -12.*log(pplev(ig,ltop(ig))/pplev(ig,1)),
c    &   devil(ig),min(sqrt(pu(ig,1)**2+pv(ig,1)**2),40.),
c    &   ptsurf(ig)-pt(ig,1),ptsurf(ig),pplev(ig,1)
c      end do    
c88    format (f7.3,1x,f7.3,1x,f6.3,1x,f6.4,1x,f7.4,1x,
c    &        f7.3,1x,f7.3,1x,f9.3)
c      do ig=1,ngridmx
c       ztop(ig) = -12.*log(pplev(ig,ltop(ig))/pplev(ig,1))
c       magwind(ig) = sqrt(pu(ig,1)**2+pv(ig,1)**2)
c       t1(ig) =max(ptsurf(ig)- pt(ig,1),0.)
c      end do

c       call WRITEDIAGFI(ngridmx,'dqs_dev','dqs devil',
c    &               'kg.m-2.s-1',2,pdqs_dev)
c       call WRITEDIAGFI(ngridmx,'wind','wind',
c    &               'm.s-1',2,magwind)
c       call WRITEDIAGFI(ngridmx,'ztop','top pbl',
c    &               'km',2,ztop)
c       call WRITEDIAGFI(ngridmx,'tsurf','tsurf',
c    &               'K',2,ptsurf)
c       call WRITEDIAGFI(ngridmx,'T1','T(1)',
c    &               'K',2,t1)
c       call WRITEDIAGFI(ngridmx,'b','b',
c    &               ' ',2,b_diag)
c     END If
c *********************************************************     
         
      RETURN
      END


