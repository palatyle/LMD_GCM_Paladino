!
       SUBROUTINE coefkzmin(knon,ypaprs,ypplay,yu,yv,yt,yq,ycdragm
     .   ,km,kn)

      USE dimphy
      IMPLICIT NONE

      include "YOMCST.h"

c.......................................................................
c  Entrees modifies en attendant une version ou les zlev, et zlay soient
c  disponibles.

      REAL  ycdragm(klon)

      REAL yu(klon,klev), yv(klon,klev)
      REAL yt(klon,klev), yq(klon,klev)
      REAL ypaprs(klon,klev+1), ypplay(klon,klev)
      REAL yustar(klon)
      real yzlay(klon,klev),yzlev(klon,klev+1),yteta(klon,klev)

      integer i

c.......................................................................
c
c  En entree :
c  -----------
c
c zlev : altitude a chaque niveau (interface inferieure de la couche
c        de meme indice)
c ustar : u*
c
c teta : temperature potentielle au centre de chaque couche
c        (en entree : la valeur au debut du pas de temps)
c
c  en sortier :
c  ------------
c
c km : diffusivite turbulente de quantite de mouvement (au bas de chaque
c      couche)
c      (en sortie : la valeur a la fin du pas de temps)
c kn : diffusivite turbulente des scalaires (au bas de chaque couche)
c      (en sortie : la valeur a la fin du pas de temps)
c
c.......................................................................

      real ustar(klon)
      real kmin,qmin,pblhmin(klon),coriol(klon)
      REAL zlev(klon,klev+1)
      REAL teta(klon,klev)

      REAL km(klon,klev+1)
      REAL kn(klon,klev+1)
      integer knon


      integer nlay,nlev
      integer ig,k

      real,parameter :: kap=0.4

      nlay=klev
      nlev=klev+1
c.......................................................................
c  en attendant une version ou les zlev, et zlay soient
c  disponibles.
c  Debut de la partie qui doit etre unclue a terme dans clmain.
c
         do i=1,knon
            yzlay(i,1)=RD*yt(i,1)/(0.5*(ypaprs(i,1)+ypplay(i,1)))
     .                *(ypaprs(i,1)-ypplay(i,1))/RG
         enddo
         do k=2,klev
            do i=1,knon
               yzlay(i,k)=yzlay(i,k-1)+RD*0.5*(yt(i,k-1)+yt(i,k))
     s                /ypaprs(i,k)*(ypplay(i,k-1)-ypplay(i,k))/RG
            enddo
         enddo
         do k=1,klev
            do i=1,knon
cATTENTION:on passe la temperature potentielle virt. pour le calcul de K
             yteta(i,k)=yt(i,k)*(ypaprs(i,1)/ypplay(i,k))**rkappa
     s          *(1.+0.61*yq(i,k))
            enddo
         enddo
         do i=1,knon
            yzlev(i,1)=0.
            yzlev(i,klev+1)=2.*yzlay(i,klev)-yzlay(i,klev-1)
         enddo
         do k=2,klev
            do i=1,knon
               yzlev(i,k)=0.5*(yzlay(i,k)+yzlay(i,k-1))
            enddo
         enddo

      yustar(1:knon) =SQRT(ycdragm(1:knon)*
     $       (yu(1:knon,1)*yu(1:knon,1)+yv(1:knon,1)*yv(1:knon,1)))

c  Fin de la partie qui doit etre unclue a terme dans clmain.

Cette routine est ecrite pour avoir en entree ustar, teta et zlev
c  Ici, on a inclut le calcul de ces trois variables dans la routine
c  coefkzmin en attendant une nouvelle version de la couche limite
c  ou ces variables seront disponibles.

c Debut de la routine coefkzmin proprement dite.

      ustar=yustar
      teta=yteta
      zlev=yzlev

      do ig=1,knon
         coriol(ig)=1.e-4
         pblhmin(ig)=0.07*ustar(ig)/max(abs(coriol(ig)),2.546e-5)
      enddo
         
      do k=2,klev
         do ig=1,knon
            if (teta(ig,2).gt.teta(ig,1)) then
               qmin=ustar(ig)*(max(1.-zlev(ig,k)/pblhmin(ig),0.))**2
               kmin=kap*zlev(ig,k)*qmin
            else
               kmin=0. ! kmin n'est utilise que pour les SL stables.
            endif 
            kn(ig,k)=kmin
            km(ig,k)=kmin
         enddo
      enddo


      return
      end
