      SUBROUTINE lectflux(irec,massem,pbarun,pbarvn,wn,tetan,phin,
     s     nrec,avant,airefi,
     s     zmfu, zmfd, zen_u, zde_u,zen_d, zde_d, coefkz,
     s     yu1,yv1,ftsol,pctsrf,
     s     frac_impa,frac_nucl,phis)

      IMPLICIT NONE

#include "dimensions.h"
#include "paramet.h"

#include "comvert.h"
#include "comconst.h"
#include "comgeom2.h"

#include "tracstoke.h"

      integer irec,nrec,i,j

      integer ngridmx,ig,l
      parameter (ngridmx=iim*(jjm-1)+2)
      INTEGER nbsrf
      PARAMETER (nbsrf=4) ! nombre de sous-fractions pour une maille

      real zmfd(ngridmx,llm),zde_d(ngridmx,llm),zen_d(ngridmx,llm)
      real zmfu(ngridmx,llm),zde_u(ngridmx,llm),zen_u(ngridmx,llm)
      real coefkz(ngridmx,llm)
      real frac_impa(ngridmx,llm),frac_nucl(ngridmx,llm)
      real yu1(ngridmx), yv1(ngridmx)
      real ftsol(ngridmx,nbsrf),pctsrf(ngridmx,nbsrf)
      integer imfu,imfd,ien_u,ide_u,
     s      ien_d,ide_d,
     s      icoefkz,izu1,izv1,
     s      itsol,ipsf,
     s      ilei, ilec
      parameter(imfu=1,imfd=llm+1,ien_u=2*llm+1,ide_u=3*llm+1,
     s      ien_d=4*llm+1,ide_d=5*llm+1,
     s      icoefkz=6*llm+1,
     s      ilei=7*llm+1,ilec=8*llm+1,
     s      izu1=9*llm+1,izv1=9*llm+2,
     s      itsol=9*llm+3,ipsf=9*llm+3+nbsrf)
      logical avant

      real massefi(ngridmx,llm)

      real massem(ip1jmp1,llm),tetan(ip1jmp1,llm)
      real pbarun(iip1,jjp1,llm),pbarvn(iip1,jjm,llm)
      real pbarvst(iip1,jjp1,llm)
      real wn(iip1,jjp1,llm),phin(iip1,jjp1,llm)
      real phis(iip1,jjp1)

      real airefi(ngridmx)

      real xlecn(ngridmx,9*llm+2+2*nbsrf)

      real zcontrole(ngridmx),zmass,tmpdyn(iip1,jjp1),zflux

      real ziadvtrac,zrec,ziadvtrac2,zrec2
      real zim,zjm,zlm,zklon,zklev

      real zpi

      zpi=2.*asin(1.)


c==================================================================
c   Si le numero du record est 0 alors: INITIALISATION
c==================================================================
c
      print*,'ENTREE DANS LECTFLUX'
        print*,'IREC=',IREC
      if(irec.eq.0) then

        print*,'IREC==',0

C test         call inigeom

c==================================================================
c   ouverture des fichiers
c==================================================================

c   Fichier fluxmass
#ifdef CRAY
         CALL ASSIGN("assign -N ieee -F null f:fluxmass")
#endif
      open(47,file='fluxmass',form='unformatted',
     s     access='direct'
     s     ,recl=4*(6*ijp1llm))
      read(47,rec=1) zrec,dtvr,ziadvtrac,zim,zjm,zlm,
     s   rlonu,rlonv,rlatu,rlatv,aire
     s    ,phis
      print*,'zrec,dtvr,ziadvtrac,zim,zjm,zlm'
      print*,zrec,dtvr,ziadvtrac,zim,zjm,zlm
      print*,rlonv



c  Fichier physique
c  Fichier lessivage (supprime les donnees utiles sont dans "physique")
#ifdef CRAY
         CALL ASSIGN("assign -N ieee -F null f:physique")
#endif
      open(49,file='physique',form='unformatted',
     s     access='direct'
     s     ,recl=4*ngridmx*(9*llm+2+2*nbsrf))
      read(49,rec=1) zrec2,ziadvtrac2,zklon,zklev
      print*,'Entete du fichier physique'
      print*,zrec2,ziadvtrac2,zklon,zklev

      nrec=zrec
      print*,'nrec=',nrec

      istdyn=ziadvtrac
      istphy=ziadvtrac2

c==================================================================
c   Fin des initialisations
      else ! irec=0
c==================================================================


c-----------------------------------------------------------------------
c   Lecture des fichiers fluxmass et  physique:
c   -----------------------------------------------------

c  Variables dynamiques
         read(47,rec=irec) massem,pbarun,pbarvst,wn,tetan,phin
        do l=1,llm
           do j=1,jjm
              do i=1,iip1
                 pbarvn(i,j,l)=pbarvst(i,j,l)
              enddo
           enddo
        enddo

c  Variables physiques
         read(49,rec=irec) ((xlecn(ig,l),ig=1,ngridmx),
     s                                    l=1,9*llm+2+2*nbsrf)

       do l=1,llm
          do ig=1,ngridmx
             coefkz(ig,l)=xlecn(ig,icoefkz+l-1)
             frac_impa(ig,l)=xlecn(ig,ilei+l-1)
             frac_nucl(ig,l)=xlecn(ig,ilec+l-1)
          enddo
       enddo
       do l=1,nbsrf
          do ig=1,ngridmx
             ftsol(ig,l)=xlecn(ig,itsol+l-1)
             pctsrf(ig,l)=xlecn(ig,ipsf+l-1)
          enddo
       enddo
       do ig=1,ngridmx
          yv1(ig)=xlecn(ig,izv1)
          yu1(ig)=xlecn(ig,izu1)
       enddo
C
      if(avant) then
c   Simu directe
       do l=1,llm
          do ig=1,ngridmx
             zmfu(ig,l)=xlecn(ig,imfu+l-1)
             zmfd(ig,l)=xlecn(ig,imfd+l-1)
             zde_u(ig,l)=xlecn(ig,ide_u+l-1)
             zen_u(ig,l)=xlecn(ig,ien_u+l-1)
             zde_d(ig,l)=xlecn(ig,ide_d+l-1)
             zen_d(ig,l)=xlecn(ig,ien_d+l-1)
          enddo
       enddo
      else
c   Simu retro
       do l=1,llm
          do ig=1,ngridmx
             zmfd(ig,l)=-xlecn(ig,imfu+l-1)
             zmfu(ig,l)=-xlecn(ig,imfd+l-1)
             zen_d(ig,l)=xlecn(ig,ide_u+l-1)
             zde_d(ig,l)=xlecn(ig,ien_u+l-1)
             zen_u(ig,l)=xlecn(ig,ide_d+l-1)
             zde_u(ig,l)=xlecn(ig,ien_d+l-1)
          enddo
       enddo
      endif

c-----------------------------------------------------------------------
c   PETIT CONTROLE SUR LES FLUX CONVECTIFS...
c-----------------------------------------------------------------------

      print*,'Ap redec irec'

         call gr_dyn_fi(llm,iip1,jjp1,ngridmx,massem,massefi)

         do ig=1,ngridmx
            zcontrole(ig)=1.
         enddo
c   zmass=(max(massem(ig,l),massem(ig,l-1))/airefi(ig)
         do l=2,llm
            do ig=1,ngridmx
               zmass=max(massefi(ig,l),massefi(ig,l-1))/airefi(ig)
               zflux=max(abs(zmfu(ig,l)),abs(zmfd(ig,l)))*dtphys
               if(zflux.gt.0.9*zmass) then
                 zcontrole(ig)=min(zcontrole(ig),0.9*zmass/zflux)
               endif
            enddo
         enddo

         do ig=1,ngridmx
            if(zcontrole(ig).lt.0.99999) then
               print*,'ATTENTION !!! on reduit les flux de masse '
               print*,'convectifs au point ig=',ig
            endif
         enddo

         call gr_fi_dyn(1,ngridmx,iip1,jjp1,zcontrole,tmpdyn)

         do l=1,llm
            do ig=1,ngridmx
               zmfu(ig,l)=zmfu(ig,l)*zcontrole(ig)
               zmfd(ig,l)=zmfd(ig,l)*zcontrole(ig)
               zen_u(ig,l)=zen_u(ig,l)*zcontrole(ig)
               zde_u(ig,l)=zde_u(ig,l)*zcontrole(ig)
               zen_d(ig,l)=zen_d(ig,l)*zcontrole(ig)
               zde_d(ig,l)=zde_d(ig,l)*zcontrole(ig)
            enddo
         enddo


      endif ! irec=0


      RETURN
      END


