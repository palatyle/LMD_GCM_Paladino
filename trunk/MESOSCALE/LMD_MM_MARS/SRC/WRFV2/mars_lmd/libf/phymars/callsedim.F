      SUBROUTINE callsedim(ngrid,nlay, ptimestep,
     $                pplev,zlev, pt,
     &                pq, pdqfi, pdqsed,pdqs_sed,nq)
      IMPLICIT NONE

c=======================================================================
c      Sedimentation of the  Martian aerosols
c      depending on their density and radius
c
c      F.Forget 1999
c
c=======================================================================

c-----------------------------------------------------------------------
c   declarations:
c   -------------

#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
#include "tracer.h"
#include "callkeys.h"

#include "fisice.h"
c
c   arguments:
c   ----------

      INTEGER ngrid,nlay
      REAL ptimestep            ! pas de temps physique (s)
      REAL pplev(ngrid,nlay+1)   ! pression aux inter-couches (Pa)
      REAL pt(ngrid,nlay)    ! temperature au centre des couches (K)
      REAL zlev(ngrid,nlay+1) ! altitude at layer boundaries


c    Traceurs :
      integer nq             ! nombre de traceurs
      real pq(ngrid,nlay,nq)  ! traceur (kg/kg)
      real pdqfi(ngrid,nlay,nq)  ! tendance avant sedimentation (kg/kg.s-1)
      real pdqsed(ngrid,nlay,nq) ! tendande due a la sedimentation (kg/kg.s-1)
      real pdqs_sed(ngrid,nq)    ! flux en surface (kg.m-2.s-1)
      
c   local:
c   ------

      INTEGER l,ig, iq

      real zqi(ngridmx,nlayermx)
      real masse (ngridmx,nlayermx)
      real epaisseur (ngridmx,nlayermx)
      real wq(ngridmx,nlayermx+1) 
c      real dens(ngridmx,nlayermx) ! Mean density of the ice part. accounting for dust core
      real rfall(ngridmx,nlayermx)


      LOGICAL firstcall
      SAVE firstcall
      DATA firstcall/.true./

c    ** un petit test de coherence
c       --------------------------

      IF (firstcall) THEN
         IF(ngrid.NE.ngridmx) THEN
            PRINT*,'STOP dans coefdifv'
            PRINT*,'probleme de dimensions :'
            PRINT*,'ngrid  =',ngrid
            PRINT*,'ngridmx  =',ngridmx
            STOP
         ENDIF
         firstcall=.false.
      ENDIF

c-----------------------------------------------------------------------
c    1. initialisation
c    -----------------
c
c     Calcul preliminaires  de caracteristiques des couches
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c    Masse (kg.m-2), epaisseur(m), temps de traversee (s)  etc...


      do  l=1,nlay
        do ig=1, ngrid
          masse(ig,l)=(pplev(ig,l) - pplev(ig,l+1)) /g 
          epaisseur(ig,l)= zlev(ig,l+1) - zlev(ig,l)
        end do
      end do

      do iq=1, nq

        if(radius(iq).gt.1.e-9) then   ! no sedim for gaz (defined by radius=0)

c         On "update" la valeur de q apres les autres parametrisations
c         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          do l=1,nlay
            do ig=1,ngrid
              zqi(ig,l)=pq(ig,l,iq)+pdqfi(ig,l,iq)*ptimestep
              if (iceparty.and.iq.eq.nq-1) then
c               On affecte un rayon moyen aux poussieres a chaque altitude du type :
c               r(z)=r0*exp(-z/H) avec r0=0.8 micron et H=18 km.
c               ''''''''''''''''''''''''''''''''''''''''''''''''
                rfall(ig,l)=max( rice(ig,l)*1.5,rdust(ig,l) )
c modif FranckMM pour ameliorer cycle H2O: rfall= 20 microns
                rfall(ig,l)=min(rfall(ig,l),1.e-4)
!mars commente pour l'instant               rfall(ig,l)=20.e-6
                
              endif
            enddo
          enddo

c         Calcul du transport par sedimentation pour chaque traceur
c         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if(iceparty.and.iq.eq.nq-1) then
            call newsedim(ngrid,nlay,ngrid*nlay,ptimestep,
     &      pplev,masse,epaisseur,pt,rfall,rho_ice,zqi,wq)
          else
            call newsedim(ngrid,nlay,1,ptimestep,
     &      pplev,masse,epaisseur,pt,radius(iq),rho_dust,zqi,wq)
          endif

c         Calcul des tendances
c         ~~~~~~~~~~~~~~~~~~~~
          do ig=1,ngrid 
            if(iceparty.and.iq.eq.nq-1) then
                pdqs_sed(ig,nq) = wq(ig,1)/ptimestep 
            else
                pdqs_sed(ig,iq) = wq(ig,1)/ptimestep 
            endif
          end do
          DO l = 1, nlay
            DO ig=1,ngrid
c	      zqi(ig,l)=max(zqi(ig,l), 1.e-10)
              pdqsed(ig,l,iq)=(zqi(ig,l)-
     $        (pq(ig,l,iq) + pdqfi(ig,l,iq)*ptimestep))/ptimestep
            ENDDO
          ENDDO
        endif
      enddo
 
      RETURN
      END

