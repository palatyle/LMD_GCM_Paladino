      SUBROUTINE vlz_fi(ngrid,q,pente_max,masse,w,wq)
c
c     Auteurs:   P.Le Van, F.Hourdin, F.Forget 
c
c    ********************************************************************
c     Shema  d'advection " pseudo amont " dans la verticale
c    pour appel dans la physique (sedimentation)
c    ********************************************************************
c    q rapport de melange (kg/kg)...
c    masse : masse de la couche Dp/g
c    w : masse d'atm ``transferee'' a chaque pas de temps (kg.m-2)
c    pente_max = 2 conseillee
c
c
c   --------------------------------------------------------------------
      IMPLICIT NONE
c
#include "dimensions.h"
#include "dimphys.h"

c
c
c   Arguments:
c   ----------
      integer ngrid
      real masse(ngrid,llm),pente_max
      REAL q(ngrid,llm)
      REAL w(ngrid,llm)
      REAL wq(ngrid,llm+1)
c
c      Local 
c   ---------
c
      INTEGER i,ij,l,j,ii
c

      real dzq(ngridmx,llm),dzqw(ngridmx,llm),adzqw(ngridmx,llm),dzqmax
      real newmasse
      real sigw, Mtot, MQtot
      integer m

      REAL      SSUM,CVMGP,CVMGT
      integer ismax,ismin


c    On oriente tout dans le sens de la pression c'est a dire dans le
c    sens de W

      do l=2,llm
         do ij=1,ngrid
            dzqw(ij,l)=q(ij,l-1)-q(ij,l)
            adzqw(ij,l)=abs(dzqw(ij,l))
         enddo
      enddo

      do l=2,llm-1
         do ij=1,ngrid
#ifdef CRAY
            dzq(ij,l)=0.5*
     ,      cvmgp(dzqw(ij,l)+dzqw(ij,l+1),0.,dzqw(ij,l)*dzqw(ij,l+1))
#else
            if(dzqw(ij,l)*dzqw(ij,l+1).gt.0.) then
                dzq(ij,l)=0.5*(dzqw(ij,l)+dzqw(ij,l+1))
            else
                dzq(ij,l)=0.
            endif
#endif
            dzqmax=pente_max*min(adzqw(ij,l),adzqw(ij,l+1))
            dzq(ij,l)=sign(min(abs(dzq(ij,l)),dzqmax),dzq(ij,l))
         enddo
      enddo

      do ij=1,ngrid
         dzq(ij,1)=0.
         dzq(ij,llm)=0.
      enddo
c ---------------------------------------------------------------
c   .... calcul des termes d'advection verticale  .......
c ---------------------------------------------------------------

c calcul de  - d( q   * w )/ d(sigma)    qu'on ajoute a  dq pour calculer dq
c
c      No flux at the model top:
       do ij=1,ngrid
          wq(ij,llm+1)=0.
       enddo

c      1) Compute wq where w > 0 (down) (ALWAYS FOR SEDIMENTATION)     
c      ===============================

       do l = 1,llm          ! loop different than when w<0
        do  ij = 1,ngrid

         if(w(ij,l).gt.0.)then

c         Regular scheme (transfered mass < 1 layer)
c         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if(w(ij,l).le.masse(ij,l))then
            sigw=w(ij,l)/masse(ij,l)
            wq(ij,l)=w(ij,l)*(q(ij,l)+0.5*(1.-sigw)*dzq(ij,l))
            

c         Extended scheme (transfered mass > 1 layer)
c         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          else 
            m=l
            Mtot = masse(ij,m)
            MQtot = masse(ij,m)*q(ij,m)
            if(m.ge.llm)goto 88
            do while(w(ij,l).gt.(Mtot+masse(ij,m+1)))
                m=m+1
                Mtot = Mtot + masse(ij,m)
                MQtot = MQtot + masse(ij,m)*q(ij,m)
                if(m.ge.llm)goto 88
            end do
 88         continue
            if (m.lt.llm) then
                sigw=(w(ij,l)-Mtot)/masse(ij,m+1)
                wq(ij,l)=(MQtot + (w(ij,l)-Mtot)*
     &          (q(ij,m+1)+0.5*(1.-sigw)*dzq(ij,m+1)) )
            else
                w(ij,l) = Mtot
                wq(ij,l) = Mqtot 
            end if
          end if
         end if
        enddo
       enddo

c      2) Compute wq where w < 0 (up) (NOT USEFUL FOR SEDIMENTATION)     
c      ===============================
       goto 99 ! SKIPPING THIS PART FOR SEDIMENTATION 

c      Surface flux up:
       do  ij = 1,ngrid
         if(w(ij,1).lt.0.) wq(ij,1)=0. ! warning : not always valid
       end do

       do l = 1,llm-1  ! loop different than when w>0
        do  ij = 1,ngrid
         if(w(ij,l+1).le.0)then

c         Regular scheme (transfered mass < 1 layer)
c         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if(-w(ij,l+1).le.masse(ij,l))then
            sigw=w(ij,l+1)/masse(ij,l)
            wq(ij,l+1)=w(ij,l+1)*(q(ij,l)-0.5*(1.+sigw)*dzq(ij,l))
c         Extended scheme (transfered mass > 1 layer)
c         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          else 
             m = l-1
             Mtot = masse(ij,m+1)
             MQtot = masse(ij,m+1)*q(ij,m+1)
             if (m.le.0)goto 77
             do while(-w(ij,l+1).gt.(Mtot+masse(ij,m)))
                m=m-1
                Mtot = Mtot + masse(ij,m+1)
                MQtot = MQtot + masse(ij,m+1)*q(ij,m+1)
                if (m.le.0)goto 77
             end do
 77          continue

             if (m.gt.0) then
                sigw=(w(ij,l+1)+Mtot)/masse(ij,m)
                wq(ij,l+1)= (MQtot + (-w(ij,l+1)-Mtot)*
     &          (q(ij,m)-0.5*(1.+sigw)*dzq(ij,m))  )
             else
c               wq(ij,l+1)= (MQtot + (-w(ij,l+1)-Mtot)*qm(ij,1))
                write(*,*) 'a rather weird situation in vlz_fi !'
                stop
             end if
          endif
         endif
        enddo
       enddo
 99    continue

      do l=1,llm
         do ij=1,ngrid

cccccccc lines below not used for sedimentation (No real flux)
ccccc       newmasse=masse(ij,l)+w(ij,l+1)-w(ij,l) 
ccccc       q(ij,l)=(q(ij,l)*masse(ij,l)+wq(ij,l+1)-wq(ij,l))
ccccc&         /newmasse
ccccc       masse(ij,l)=newmasse

            q(ij,l)=q(ij,l) +  (wq(ij,l+1)-wq(ij,l))/masse(ij,l)

         enddo
      enddo



      return
      end
