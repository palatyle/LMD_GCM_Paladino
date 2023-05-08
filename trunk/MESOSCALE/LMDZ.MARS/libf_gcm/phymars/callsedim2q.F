      SUBROUTINE callsedim2q(ngrid,nlay, ptimestep,
     $                pplev,zlev, pt,
     &                pq, pdqfi, pdqsed,pdqs_sed,nq)
      IMPLICIT NONE

c=======================================================================
c
c      Sedimentation of the  Martian dust 
c      using method with mass (iq=1) and number(iq=2) mixing ratio
c
c      F.Forget 1999
c
c
c=======================================================================

c-----------------------------------------------------------------------
c   declarations:
c   -------------

#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
#include "tracer.h"
c
c   arguments:
c   ----------

      INTEGER ngrid,nlay
      REAL ptimestep             ! pas de temps physique (s)
      REAL pplev(ngrid,nlay+1)   ! pression aux inter-couches (Pa)
      REAL pt(ngrid,nlay)        ! temperature au centre des couches (K)
      REAL zlev(ngrid,nlay+1)    ! altitude at layer boundaries


c    Traceurs :
      integer nq                 ! nombre de traceurs
      real pq(ngrid,nlay,nq)     ! traceur (kg/kg)
      real pdqfi(ngrid,nlay,nq)  ! tendance avant sedimentation (kg/kg.s-1)
      real pdqsed(ngrid,nlay,nq) ! tendande due a la sedimentation (kg/kg.s-1)
      real pdqs_sed(ngrid,nq)    ! flux en surface (kg.m-2.s-1)
      
c   local:
c   ------

      INTEGER l,ig

      LOGICAL firstcall
      SAVE firstcall

c    Traceurs :
c    ~~~~~~~~ 
      INTEGER iq
      real zq(ngridmx,nlayermx,2)
      real masse(ngridmx,nlayermx)      ! Layer mass (kg.m-2)
      real epaisseur(ngridmx,nlayermx)  ! Layer thickness (m)
      real wq(ngridmx,nlayermx+1)       ! moved tracer mass (kg.m-2)
      real r0(ngridmx,nlayermx)         ! geometric mean radius (m)


c     Discretisation en taille:
c     ~~~~~~~~~~~~~~
c       1) Discretisation pour representer la variation de Vitesse de chute
      integer nr,ir
      parameter (nr=7)   ! nombre de rayon pour le calcul du transport
      real rd(nr),qr(ngridmx,nlayermx,nr)
      real rdi(nr+1)    ! extreme and intermediate radii
      real Sq(ngridmx,nlayermx)
      real rdmin,rdmax,rdimin,rdimax
      data rdmin/1.e-7/
      data rdmax/30.e-6/
      data rdimin/1.e-10/
      data rdimax/1e-4/
      save rd, rdi

c       2) Sous discretisation pour bien integrer la loi lognormale
c          (calcul de q pour chaque gamme de rayon rd)

      integer nint, iint   
      parameter (nint=4)  ! nombre de point entre chaque rayon rdi
      real rr(nint,nr)
      save rr
      
      real reff(ngridmx,nlayermx,2)  ! for diagnostic only
      



      DATA firstcall/.true./

c    ** un petit test de coherence
c       --------------------------

      IF (firstcall) THEN
         IF(ngrid.NE.ngridmx) THEN
            Print*,'STOP dans coefdifv'
            Print*,'probleme de dimensions :'
            Print*,'ngrid  =',ngrid
            Print*,'ngridmx  =',ngridmx
            STOP
         ENDIF
         firstcall=.false.
 
         do ir=1,nr
             rd(ir)= rdmin*(rdmax/rdmin)**(float(ir-1)/float(nr-1))
c            rd(ir) = rdmin + (float(ir-1)/float(nr-1))*(rdmax-rdmin)
         end do
         rdi(1)=rdimin
         do ir=2,nr
           rdi(ir)= sqrt(rd(ir-1)*rd(ir))
         end do
         rdi(nr+1)=rdimax

         do ir=1,nr
           do iint=1,nint
             rr(iint,ir)=
     &       rdi(ir)*(rdi(ir+1)/rdi(ir))**(float(iint-1)/float(nint-1))
c             write(*,*) rr(iint,ir)
           end do
         end do
      ENDIF





c-----------------------------------------------------------------------
c    1. initialisation
c    -----------------

c    On "update" la valeur de q apres les autres parametrisations
c    ~~~~~~~~~~~~~~~~~~~~~~~~~~
      do iq=1,2
        do l=1,nlay
          do ig=1,ngrid
            zq(ig,l,iq)=pq(ig,l,iq)+pdqfi(ig,l,iq)*ptimestep
          enddo
        enddo
      enddo


c    Calcul preliminaires  de caracteristiques des couches
c    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c    Masse (kg.m-2), epaisseur(m), temps de traversee (s)  etc...

      do  l=1,nlay
        do ig=1, ngrid
          masse(ig,l)=(pplev(ig,l) - pplev(ig,l+1)) /g
          epaisseur(ig,l)= zlev(ig,l+1) - zlev(ig,l)
        end do
      end do


c     Calcul de la distribution de taille
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      do  l=1,nlay
         do ig=1, ngrid
              r0(ig,l)=
     &        (r3n_q*zq(ig,l,1)/max(zq(ig,l,min(2,nq)),0.01))**(1./3.)
              r0(ig,l)=min(max(r0(ig,l),1.e-10),500.e-6)
         end do
      end do

      DO iq=1,2  ! LOOP on Mass (iq=1) then Number(iq=2) mixing ratio

c        Calcul du rapport de melange pour les nir tailles considerees
c        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call zerophys(ngridmx*nlayermx, Sq)
          do ir=1,nr
            do l=1,nlay
              do ig=1,ngrid
c                **************** nouvelle methode :
c                Calcul de q pour chaque gamme de rayon rd
c                Sous discretisation pour bien integrer la loi lognormale
c                (Methode des trapezes)
                 qr(ig,l,ir)=0.5*(rr(2,ir)-rr(1,ir))*
     &             (rr(1,ir)**(5-3*iq))*
     &             exp(-(log(rr(1,ir)/r0(ig,l)))**2/(2*varian**2))
                 do iint=2,nint-1
                   qr(ig,l,ir)=qr(ig,l,ir) +
     &             0.5*(rr(iint+1,ir)-rr(iint-1,ir))*
     &             (rr(iint,ir)**(5-3*iq))*
     &             exp(-(log(rr(iint,ir)/r0(ig,l)))**2/(2*varian**2))
                 end do
                   qr(ig,l,ir)=qr(ig,l,ir) +
     &             0.5*(rr(nint,ir)-rr(nint-1,ir))*
     &             (rr(nint,ir)**(5-3*iq))*
     &             exp(-(log(rr(nint,ir)/r0(ig,l)))**2/(2*varian**2))

c                **************** Methode simple BUGGEE
c                qr(ig,l,ir)=(rd(ir)**(5-3*iq))*
c    &           exp( -(log(rd(ir)/r0(ig,l)))**2 / (2*varian**2) )
c                ******************************
                 Sq(ig,l)=Sq(ig,l)+qr(ig,l,ir) 
              enddo
            enddo
          enddo
          do ir=1,nr
            do l=1,nlay
              do ig=1,ngrid
                 qr(ig,l,ir)= zq(ig,l,iq)* qr(ig,l,ir)/Sq(ig,l)
              enddo
            enddo
c          write(123,*) rd(ir)*1.e6, qr(1,10,ir)
          enddo

c         Calcul du transport par sedimentation pour chaque traceur
c         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          call zerophys(ngridmx*nlayermx,zq(1,1,iq))
          call zerophys(ngridmx*nlayermx,pdqs_sed(1,iq))

c         write(*,*)
c         if(iq.eq.1) write(*,*) '5 fois  wq(6), wq(7), q(6)'
          
          do ir=1,nr
             call newsedim(ngrid,nlay,1,ptimestep,
     &       pplev,masse,epaisseur,pt,rd(ir),rho_dust,qr(1,1,ir),wq)

c            Calcul des tendances
c            ~~~~~~~~~~~~~~~~~~~~
             do ig=1,ngrid
               pdqs_sed(ig,iq) = pdqs_sed(ig,iq) + wq(ig,1)/ptimestep
             end do
             DO l = 1, nlay
               DO ig=1,ngrid
                 zq(ig,l,iq)=zq(ig,l,iq)+qr(ig,l,ir)
               ENDDO
             ENDDO
c         if(iq.eq.1) write(*,*) 'callsed: wq(6), wq(7), q(6)',
c    &               wq(1,6),wq(1,7),qr(1,6,ir)

c               if(iq.eq.1)write(*,*)'R=', rd(ir)*1.e6
c               write(*,*)
          enddo

          do l = 1, nlay
            do ig=1,ngrid
              pdqsed(ig,l,iq)=(zq(ig,l,iq)-
     $        (pq(ig,l,iq) + pdqfi(ig,l,iq)*ptimestep))/ptimestep
            enddo
          enddo

c     ***** Diagnostique *****************
c         call zerophys(ngridmx*nlayermx,reff(1,1,iq))
c         do ir=1,nr
c           do l = 1, nlay
c             do ig=1,ngrid
c              if(iq.eq.1) then 
c               reff(ig,l,iq)=reff(ig,l,iq) + qr(ig,l,ir)/rd(ir)
c              else
c               reff(ig,l,iq)=reff(ig,l,iq) + qr(ig,l,ir)*log(rd(ir))
c              endif
c             enddo
c           enddo
c         enddo
c         do l = 1, nlay
c           do ig=1,ngrid
c            if(iq.eq.1) then 
c             reff(ig,l,iq)=zq(ig,l,iq)/reff(ig,l,iq)
c            else
c             reff(ig,l,iq)=exp(reff(ig,l,iq)/zq(ig,l,iq))
c            endif
c           enddo
c         enddo
c     **** FIN Diagnostique *****************

      enddo  ! end loop on iq=1,2



c ********************************************
c     Diagnostique
c     if(ngrid.eq.1) then
c      do  l=1,nlay
c       do ig=1, ngrid
c         r0(ig,l)=
c    &    (r3n_q*zq(ig,l,1)/max(zq(ig,l,min(2,nq)),0.01))**(1./3.)
c         r0(ig,l)=min(max(r0(ig,l),1.e-10),500.e-6)
c
c         write(*,*)'l,r0, r02 ', l, r0(ig,l),reff(ig,l,1)*0.3626
c    &    ,reff(ig,l,2)
c    &    ,r0(ig,l)/(reff(ig,l,1)*0.3626) 
c
c       end do
c      end do
c      CALL writeg1d(ngrid,nlay,r0,'r0','m')
c     end if
c ********************************************


      RETURN
      END

