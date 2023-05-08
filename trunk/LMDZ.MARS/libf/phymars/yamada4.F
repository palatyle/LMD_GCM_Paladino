!************************************************************
!************************************************************
!
!      YAMADA4 EARTH =>>> MARS VERSION
!      Modifications by: A.C. 02-03-2012 (marked by 'MARS')
!      Original version given by F.H. 01-03-2012
!
!************************************************************
!************************************************************
      SUBROUTINE yamada4(ngrid,nlay,nq,dt,g,rconst,plev,temp
     s   ,zlev,zlay,u,v,phc,pq,cd,q2,km,kn,kq,ustar
     s   ,iflag_pbl)
      use tracer_mod, only: noms
      use turb_mod, only: l0
      IMPLICIT NONE
!.......................................................................
! MARS
#include "callkeys.h"
!.......................................................................
!
! dt : pas de temps
! g  : g
! zlev : altitude a chaque niveau (interface inferieure de la couche
!        de meme indice)
! zlay : altitude au centre de chaque couche
! u,v : vitesse au centre de chaque couche
!       (en entree : la valeur au debut du pas de temps)
! phc : temperature potentielle au centre de chaque couche
!        (en entree : la valeur au debut du pas de temps)
! cd : cdrag
!      (en entree : la valeur au debut du pas de temps)
! q2 : $q^2$ au bas de chaque couche
!      (en entree : la valeur au debut du pas de temps)
!      (en sortie : la valeur a la fin du pas de temps)
! km : diffusivite turbulente de quantite de mouvement (au bas de chaque
!      couche)
!      (en sortie : la valeur a la fin du pas de temps)
! kn : diffusivite turbulente des scalaires (au bas de chaque couche)
!      (en sortie : la valeur a la fin du pas de temps)
!
!  iflag_pbl doit valoir entre 6 et 9
!      l=6, on prend  systematiquement une longueur d equilibre
!    iflag_pbl=6 : MY 2.0
!    iflag_pbl=7 : MY 2.0.Fournier
!    iflag_pbl=8 : MY 2.5
!    iflag_pbl>=9 : MY 2.5 avec diffusion verticale
!
!.......................................................................

      REAL,    INTENT(IN)    :: dt,g,rconst
      REAL,    INTENT(IN)    :: u(ngrid,nlay)
      REAL,    INTENT(IN)    :: v(ngrid,nlay)
      REAL,    INTENT(IN)    :: phc(ngrid,nlay)
      REAL,    INTENT(IN)    :: cd(ngrid)
      REAL,    INTENT(IN)    :: temp(ngrid,nlay)
      REAL,    INTENT(IN)    :: plev(ngrid,nlay+1)
      REAL,    INTENT(IN)    :: ustar(ngrid)
      REAL,    INTENT(IN)    :: zlev(ngrid,nlay+1)
      REAL,    INTENT(IN)    :: zlay(ngrid,nlay)
      INTEGER, INTENT(IN)    :: iflag_pbl,ngrid
      INTEGER, INTENT(IN)    :: nlay
      INTEGER, INTENT(IN)    :: nq
      REAL,    INTENT(INOUT) :: q2(ngrid,nlay+1)
      REAL,    INTENT(OUT)   :: km(ngrid,nlay+1)
      REAL,    INTENT(OUT)   :: kn(ngrid,nlay+1)
      REAL,    INTENT(OUT)   :: kq(ngrid,nlay+1)

      REAL kmin,qmin,pblhmin(ngrid),coriol(ngrid)
      REAL unsdz(ngrid,nlay)
      REAL unsdzdec(ngrid,nlay+1)
      REAL kmpre(ngrid,nlay+1),tmp2
      REAL mpre(ngrid,nlay+1)
      REAL ff(ngrid,nlay+1),delta(ngrid,nlay+1)
      REAL aa(ngrid,nlay+1),aa0,aa1,qpre

      LOGICAL first
      INTEGER ipas,nlev
      SAVE first,ipas
!FH/IM     DATA first,ipas/.true.,0/
      DATA first,ipas/.false.,0/
      INTEGER ig,k


      REAL ri,zrif,zalpha,zsm,zsn
      REAL rif(ngrid,nlay+1),sm(ngrid,nlay+1),alpha(ngrid,nlay)

      REAL m2(ngrid,nlay+1),dz(ngrid,nlay+1),zq,n2(ngrid,nlay+1)
      REAL dtetadz(ngrid,nlay+1)
      REAL m2cstat,mcstat,kmcstat
      REAL l(ngrid,nlay+1)
      REAL sq(ngrid),sqz(ngrid),zz(ngrid,nlay+1)
      INTEGER iter

      REAL ric,rifc,b1,kap
      SAVE ric,rifc,b1,kap
      DATA ric,rifc,b1,kap/0.195,0.191,16.6,0.4/
      REAL frif,falpha,fsm
      REAL fl,zzz,zl0,zq2,zn2

      REAL rino(ngrid,nlay+1),smyam(ngrid,nlay),styam(ngrid,nlay)
     s  ,lyam(ngrid,nlay),knyam(ngrid,nlay)
     s  ,w2yam(ngrid,nlay),t2yam(ngrid,nlay)
      LOGICAL,SAVE :: firstcall=.true.
      frif(ri)=0.6588*(ri+0.1776-sqrt(ri*ri-0.3221*ri+0.03156))
      falpha(ri)=1.318*(0.2231-ri)/(0.2341-ri)
      fsm(ri)=1.96*(0.1912-ri)*(0.2341-ri)/((1.-ri)*(0.2231-ri))
      fl(zzz,zl0,zq2,zn2)=
     s     max(min(l0(ig)*kap*zlev(ig,k)/(kap*zlev(ig,k)+l0(ig))
     s     ,0.5*sqrt(q2(ig,k))/sqrt(max(n2(ig,k),1.e-10))) ,1.)


! MARS
      REAL,SAVE :: q2min,q2max,knmin,kmmin
      DATA q2min,q2max,knmin,kmmin/1.E-10,1.E+2,1.E-5,1.E-5/
      INTEGER ico2,iq
      SAVE ico2
      REAL m_co2, m_noco2, A , B
      SAVE A, B
      REAL teta(ngrid,nlay)
      REAL pq(ngrid,nlay,nq)
      REAL kminfact
      INTEGER i
      REAL ztimestep
      INTEGER :: ndt

      nlev=nlay+1

c.......................................................................
c  Initialization
c.......................................................................

      !! firstcall: OK absolute
      if(firstcall) then
        ico2=0
        if (tracer) then
!     Prepare Special treatment if one of the tracers is CO2 gas
           do iq=1,nq
             if (noms(iq).eq."co2") then
                ico2=iq
                m_co2 = 44.01E-3  ! CO2 molecular mass (kg/mol)
                m_noco2 = 33.37E-3  ! Non condensible mol mass (kg/mol)
!               Compute A and B coefficient use to compute
!               mean molecular mass Mair defined by
!               1/Mair = q(ico2)/m_co2 + (1-q(ico2))/m_noco2
!               1/Mair = A*q(ico2) + B
                A =(1/m_co2 - 1/m_noco2)
                B=1/m_noco2
             end if
           enddo
        endif
        firstcall=.false.
      endif !of if firstcall

      !! AS: moved out of firstcall to allow nesting+evoluting timestep
      ndt=ceiling(3840./(3699.*24./dt))

c.......................................................................
c  Special treatment for co2
c.......................................................................

!      if (ico2.ne.0) then
!!     Special case if one of the tracers is CO2 gas
!         DO k=1,nlay
!           DO ig=1,ngrid
!            teta(ig,k) = phc(ig,k)*(A*pq(ig,k,ico2)+B)
!           ENDDO
!         ENDDO
!       else
          teta(:,:)=phc(:,:)
!       end if
      
      if (.not.(iflag_pbl.ge.6.and.iflag_pbl.le.10)) then
           stop'probleme de coherence dans appel a MY'
      endif

      ipas=ipas+1
! MARS
!      if (0.eq.1.and.first) then
!      do ig=1,1000
!         ri=(ig-800.)/500.
!         if (ri.lt.ric) then
!            zrif=frif(ri)
!         else
!            zrif=rifc
!         endif
!         if(zrif.lt.0.16) then
!            zalpha=falpha(zrif)
!            zsm=fsm(zrif)
!         else
!            zalpha=1.12
!            zsm=0.085
!         endif
!     print*,ri,rif,zalpha,zsm
!      enddo
!      endif

!.......................................................................
!  les increments verticaux
!.......................................................................
!
!!!!!! allerte !!!!!
!!!!!! zlev n'est pas declare a nlev !!!!!
!!!!!! ---->
! MARS
!
!                                                      DO ig=1,ngrid
!            zlev(ig,nlev)=zlay(ig,nlay)
!     &             +( zlay(ig,nlay) - zlev(ig,nlev-1) )
!                                                      ENDDO
!!!!! <----
!!!!! allerte !!!!!

      DO k=1,nlay
                                                      DO ig=1,ngrid
        unsdz(ig,k)=1.E+0/(zlev(ig,k+1)-zlev(ig,k))
                                                      ENDDO
      ENDDO
                                                      DO ig=1,ngrid
      unsdzdec(ig,1)=1.E+0/(zlay(ig,1)-zlev(ig,1))
                                                      ENDDO
      DO k=2,nlay
                                                      DO ig=1,ngrid
        unsdzdec(ig,k)=1.E+0/(zlay(ig,k)-zlay(ig,k-1))
                                                     ENDDO
      ENDDO
                                                      DO ig=1,ngrid
      unsdzdec(ig,nlay+1)=1.E+0/(zlev(ig,nlay+1)-zlay(ig,nlay))
                                                     ENDDO
!
!.......................................................................

      do k=2,nlay
                                                          do ig=1,ngrid
         dz(ig,k)=zlay(ig,k)-zlay(ig,k-1)
         m2(ig,k)=((u(ig,k)-u(ig,k-1))**2+(v(ig,k)-v(ig,k-1))**2)
     s             /(dz(ig,k)*dz(ig,k))
         dtetadz(ig,k)=(teta(ig,k)-teta(ig,k-1))/dz(ig,k)
         n2(ig,k)=g*2.*dtetadz(ig,k)/(teta(ig,k-1)+teta(ig,k))
!        n2(ig,k)=0.
         ri=n2(ig,k)/max(m2(ig,k),1.e-10)
         if (ri.lt.ric) then
            rif(ig,k)=frif(ri)
         else
            rif(ig,k)=rifc
         endif
         if(rif(ig,k).lt.0.16) then
            alpha(ig,k)=falpha(rif(ig,k))
            sm(ig,k)=fsm(rif(ig,k))
         else
            alpha(ig,k)=1.12
            sm(ig,k)=0.085
         endif
         zz(ig,k)=b1*m2(ig,k)*(1.-rif(ig,k))*sm(ig,k)
!     print*,'RIF L=',k,rif(ig,k),ri*alpha(ig,k)


                                                          enddo
      enddo


!====================================================================
!   Au premier appel, on determine l et q2 de facon iterative.
! iterration pour determiner la longueur de melange


      if (first.or.iflag_pbl.eq.6) then
                                                          do ig=1,ngrid
! MARS 
!      l0(ig)=10.
      l0(ig)=160.
                                                          enddo
      do k=2,nlay-1
                                                          do ig=1,ngrid
        l(ig,k)=l0(ig)*kap*zlev(ig,k)/(kap*zlev(ig,k)+l0(ig))
                                                          enddo
      enddo

      do iter=1,10
                                                          do ig=1,ngrid
         sq(ig)=1.e-10
         sqz(ig)=1.e-10
                                                          enddo
         do k=2,nlay-1
                                                          do ig=1,ngrid
           q2(ig,k)=l(ig,k)**2*zz(ig,k)
           l(ig,k)=fl(zlev(ig,k),l0(ig),q2(ig,k),n2(ig,k))
           zq=sqrt(q2(ig,k))
           sqz(ig)=sqz(ig)+zq*zlev(ig,k)*(zlay(ig,k)-zlay(ig,k-1))
           sq(ig)=sq(ig)+zq*(zlay(ig,k)-zlay(ig,k-1))
                                                          enddo
         enddo
                                                          do ig=1,ngrid
         l0(ig)=0.2*sqz(ig)/sq(ig)
!        l0(ig)=30.
                                                          enddo
!      print*,'ITER=',iter,'  L0=',l0

      enddo

!     print*,'Fin de l initialisation de q2 et l0'

      endif ! first

!====================================================================
!  Calcul de la longueur de melange.
!====================================================================

!   Mise a jour de l0
                                                          do ig=1,ngrid
      sq(ig)=1.e-10
      sqz(ig)=1.e-10
                                                          enddo
      do k=2,nlay-1
                                                          do ig=1,ngrid
        zq=sqrt(q2(ig,k))
        sqz(ig)=sqz(ig)+zq*zlev(ig,k)*(zlay(ig,k)-zlay(ig,k-1))
        sq(ig)=sq(ig)+zq*(zlay(ig,k)-zlay(ig,k-1))
                                                          enddo
      enddo
                                                          do ig=1,ngrid
      l0(ig)=0.2*sqz(ig)/sq(ig)
!        l0(ig)=30.
                                                          enddo
!      print*,'ITER=',iter,'  L0=',l0
!   calcul de l(z)
      do k=2,nlay
                                                          do ig=1,ngrid
         l(ig,k)=fl(zlev(ig,k),l0(ig),q2(ig,k),n2(ig,k))
         if(first) then
           q2(ig,k)=l(ig,k)**2*zz(ig,k)
         endif
                                                          enddo
      enddo

!====================================================================
!   Yamada 2.0
!====================================================================
      if (iflag_pbl.eq.6) then

      do k=2,nlay
                                                          do ig=1,ngrid
         q2(ig,k)=l(ig,k)**2*zz(ig,k)
                                                          enddo
      enddo


      else if (iflag_pbl.eq.7) then
!====================================================================
!   Yamada 2.Fournier
!====================================================================

!  Calcul de l,  km, au pas precedent
      do k=2,nlay
                                                          do ig=1,ngrid
c        print*,'SMML=',sm(ig,k),l(ig,k)
         delta(ig,k)=q2(ig,k)/(l(ig,k)**2*sm(ig,k))
         kmpre(ig,k)=l(ig,k)*sqrt(q2(ig,k))*sm(ig,k)
         mpre(ig,k)=sqrt(m2(ig,k))
c        print*,'0L=',k,l(ig,k),delta(ig,k),km(ig,k)
                                                          enddo
      enddo

      do k=2,nlay-1
                                                          do ig=1,ngrid
        m2cstat=max(alpha(ig,k)*n2(ig,k)+delta(ig,k)/b1,1.e-12)
        mcstat=sqrt(m2cstat)

!        print*,'M2 L=',k,mpre(ig,k),mcstat
!
!  -----{puis on ecrit la valeur de q qui annule l'equation de m
!        supposee en q3}
!
        IF (k.eq.2) THEN
          kmcstat=1.E+0 / mcstat
     &    *( unsdz(ig,k)*kmpre(ig,k+1)
     &                        *mpre(ig,k+1)
     &      +unsdz(ig,k-1)
     &              *cd(ig)
     &              *( sqrt(u(ig,3)**2+v(ig,3)**2)
     &                -mcstat/unsdzdec(ig,k)
     &                -mpre(ig,k+1)/unsdzdec(ig,k+1) )**2)
     &      /( unsdz(ig,k)+unsdz(ig,k-1) )
        ELSE
          kmcstat=1.E+0 / mcstat
     &    *( unsdz(ig,k)*kmpre(ig,k+1)
     &                        *mpre(ig,k+1)
     &      +unsdz(ig,k-1)*kmpre(ig,k-1)
     &                          *mpre(ig,k-1) )
     &      /( unsdz(ig,k)+unsdz(ig,k-1) )
        ENDIF
!       print*,'T2 L=',k,tmp2
        tmp2=kmcstat
     &      /( sm(ig,k)/q2(ig,k) )
     &      /l(ig,k)

! MARS
!        q2(ig,k)=max(tmp2,1.e-12)**(2./3.)
        q2(ig,k)=max(q2min,max(tmp2,1.e-12)**(2./3.))

!       print*,'Q2 L=',k,q2(ig,k)
!
                                                          enddo
      enddo

      else if (iflag_pbl.ge.8) then
!====================================================================
!   Yamada 2.5 a la Didi
!====================================================================

      ztimestep=dt/real(ndt)
      do i=1,ndt

!  Calcul de l,  km, au pas precedent
      do k=2,nlay
       do ig=1,ngrid
!        print*,'SMML=',sm(ig,k),l(ig,k)
         delta(ig,k)=q2(ig,k)/(l(ig,k)**2*sm(ig,k))
         if (delta(ig,k).lt.1.e-20) then
!     print*,'ATTENTION   L=',k,'   Delta=',delta(ig,k)
            delta(ig,k)=1.e-20
         endif
         km(ig,k)=l(ig,k)*sqrt(q2(ig,k))*sm(ig,k)
         aa0=
     s   (m2(ig,k)-alpha(ig,k)*n2(ig,k)-delta(ig,k)/b1)
         aa1=
     s   (m2(ig,k)*(1.-rif(ig,k))-delta(ig,k)/b1)
! abder      print*,'AA L=',k,aa0,aa1,aa1/max(m2(ig,k),1.e-20)
         aa(ig,k)=aa1*ztimestep/(delta(ig,k)*l(ig,k))
!     print*,'0L=',k,l(ig,k),delta(ig,k),km(ig,k)
         qpre=sqrt(q2(ig,k))
         if (iflag_pbl.eq.8 ) then
            if (aa(ig,k).gt.0.) then
               q2(ig,k)=(qpre+aa(ig,k)*qpre*qpre)**2
            else
               q2(ig,k)=(qpre/(1.-aa(ig,k)*qpre))**2
            endif
         else ! iflag_pbl=9
            if (aa(ig,k)*qpre.gt.0.9) then
               q2(ig,k)=(qpre*10.)**2
            else
               q2(ig,k)=(qpre/(1.-aa(ig,k)*qpre))**2
            endif
         endif

! MARS
         q2(ig,k)=min(max(q2(ig,k),q2min),q2max)
!         q2(ig,k)=min(max(q2(ig,k),1.e-10),1.e4)

!     print*,'Q2 L=',k,q2(ig,k),qpre*qpre
       enddo
      enddo

! MARS
      q2(:,nlay+1)=q2(:,nlay)

      if (iflag_pbl .eq. 9) then
      do k=2,nlay
      do ig=1,ngrid
        zq=sqrt(q2(ig,k))
        km(ig,k)=l(ig,k)*zq*sm(ig,k)
        kn(ig,k)=km(ig,k)*alpha(ig,k)
        kq(ig,k)=l(ig,k)*zq*0.2
      enddo
      enddo
      ! boundary conditions for km
      km(:,nlay+1)=0
      km(:,1)=km(:,2) ! km(:,1)=0
      ! boundary conditions for kn
      kn(:,nlay+1)=0
      kn(:,1)=kn(:,2) ! kn(:,1)=0
      ! boundary conditions for kq
      kq(:,nlay+1)=0  ! zero at top of atmosphere
      kq(:,1)=kq(:,2) ! no gradient at surface

      q2(:,1)=q2(:,2)
       call vdif_q2(ztimestep,g,rconst,ngrid,nlay,plev,temp,kq,q2)

      endif ! of if iflag_pbl eq 9

      enddo !of i=1,ndt

      endif ! Fin du cas 8

!     print*,'OK8'

!====================================================================
!   Calcul des coefficients de melange
!====================================================================
      if (iflag_pbl .ne. 9) then
      do k=2,nlay
!     print*,'k=',k
                                                          do ig=1,ngrid
!abde      print*,'KML=',l(ig,k),q2(ig,k),sm(ig,k)
         zq=sqrt(q2(ig,k))
         km(ig,k)=l(ig,k)*zq*sm(ig,k)
         kn(ig,k)=km(ig,k)*alpha(ig,k)
         kq(ig,k)=l(ig,k)*zq*0.2
!     print*,'KML=',km(ig,k),kn(ig,k)
                                                          enddo
      enddo

! MARS
      km(:,nlay+1)=km(:,nlay)
      kn(:,nlay+1)=kn(:,nlay)
      kq(:,nlay+1)=kq(:,nlay)

! Transport diffusif vertical de la TKE.
!      if (iflag_pbl.ge.9) then
!!       print*,'YAMADA VDIF'
!        q2(:,1)=q2(:,2)
!        call vdif_q2(dt,g,rconst,ngrid,nlay,plev,temp,kq,q2)
!      endif

      endif

!   Traitement des cas noctrunes avec l'introduction d'une longueur
!   minilale.
!
!====================================================================
!   Traitement particulier pour les cas tres stables.
!   D'apres Holtslag Boville.

! MARS
!       callkmin=.true.
!       call getin("callkmin",callkmin)
!       IF (callkmin) THEN
                                                          do ig=1,ngrid
!      coriol(ig)=1.e-4
!      pblhmin(ig)=0.07*ustar(ig)/max(abs(coriol(ig)),2.546e-5)

       if (ngrid .eq. 1) then
       kminfact=0.3
       else
       kminfact=0.45
       endif

       pblhmin(ig)=kminfact*0.07*MAX(ustar(ig),1.e-3)/1.e-4
                                                   enddo
!      print*,'pblhmin ',pblhmin
!CTest a remettre 21 11 02
! test abd 13 05 02      if(0.eq.1) then
!      if(0.eq.1) then
      do k=2,nlay
         do ig=1,ngrid
            if (teta(ig,2).gt.teta(ig,1)) then
               qmin=ustar(ig)*(max(1.-zlev(ig,k)/pblhmin(ig),0.))**2
!               kmin=kap*zlev(ig,k)*qmin
               kmin=fl(zlev(ig,k),l0(ig),qmin**2,n2(ig,k))*qmin
            else
               kmin=-1. ! kmin n'est utilise que pour les SL stables.
            endif
            if (kn(ig,k).lt.kmin.or.km(ig,k).lt.kmin) then
!               print*,'Seuil min Km K=',k,kmin,km(ig,k),kn(ig,k)
!     s           ,sqrt(q2(ig,k)),pblhmin(ig),qmin/sm(ig,k)
!               kn(ig,k)=kmin
!               km(ig,k)=kmin
!               kq(ig,k)=kmin

               kn(ig,k)=kmin*alpha(ig,k)
               km(ig,k)=kmin
               kq(ig,k)=kmin*0.2
!   la longueur de melange est suposee etre l= kap z
!   K=l q Sm d'ou q2=(K/l Sm)**2
!               q2(ig,k)=(qmin/sm(ig,k))**2
               q2(ig,k)=(kmin/
     &     (fl(zlev(ig,k),l0(ig),qmin**2,n2(ig,k))*sm(ig,k)))**2
            endif
         enddo
      enddo
!      endif

!      ENDIF

!   Diagnostique pour stokage

      if(1.eq.0)then
      rino=rif
      smyam(1:ngrid,1)=0.
      styam(1:ngrid,1)=0.
      lyam(1:ngrid,1)=0.
      knyam(1:ngrid,1)=0.
      w2yam(1:ngrid,1)=0.
      t2yam(1:ngrid,1)=0.

      smyam(1:ngrid,2:nlay)=sm(1:ngrid,2:nlay)
      styam(1:ngrid,2:nlay)=sm(1:ngrid,2:nlay)*alpha(1:ngrid,2:nlay)
      lyam(1:ngrid,2:nlay)=l(1:ngrid,2:nlay)
      knyam(1:ngrid,2:nlay)=kn(1:ngrid,2:nlay)

!   Estimations de w'2 et T'2 d'apres Abdela et McFarlane

      w2yam(1:ngrid,2:nlay)=q2(1:ngrid,2:nlay)*0.24
     s    +lyam(1:ngrid,2:nlay)*5.17*kn(1:ngrid,2:nlay)
     s    *n2(1:ngrid,2:nlay)/sqrt(q2(1:ngrid,2:nlay))

      t2yam(1:ngrid,2:nlay)=9.1*kn(1:ngrid,2:nlay)
     s    *dtetadz(1:ngrid,2:nlay)**2
     s    /sqrt(q2(1:ngrid,2:nlay))*lyam(1:ngrid,2:nlay)
      endif

!     print*,'OKFIN'
      first=.false.
      return
      end
      SUBROUTINE vdif_q2(timestep,gravity,rconst,ngrid,nlay
     & ,plev,temp,kmy,q2)
      IMPLICIT NONE
!.......................................................................
! MARS
#include "callkeys.h"
!.......................................................................
!
! dt : pas de temps
!
      REAL plev(ngrid,nlay+1)
      REAL temp(ngrid,nlay)
      REAL timestep
      REAL gravity,rconst
      REAL kstar(ngrid,nlay+1),zz
      REAL kmy(ngrid,nlay+1)
      REAL q2(ngrid,nlay+1)
      REAL deltap(ngrid,nlay+1)
      REAL denom(ngrid,nlay+1),alpha(ngrid,nlay+1),beta(ngrid,nlay+1)
      INTEGER ngrid,nlay

      INTEGER i,k

! 	print*,'RD=',rconst
      do k=1,nlay
         do i=1,ngrid
! test
!       print*,'i,k',i,k
! 	print*,'temp(i,k)=',temp(i,k)
! 	print*,'(plev(i,k)-plev(i,k+1))=',plev(i,k),plev(i,k+1)
            zz=(plev(i,k)+plev(i,k+1))*gravity/(rconst*temp(i,k))
            kstar(i,k)=0.125*(kmy(i,k+1)+kmy(i,k))*zz*zz
     s      /(plev(i,k)-plev(i,k+1))*timestep
         enddo
      enddo

      do k=2,nlay
         do i=1,ngrid
            deltap(i,k)=0.5*(plev(i,k-1)-plev(i,k+1))
         enddo
      enddo
      do i=1,ngrid
         deltap(i,1)=0.5*(plev(i,1)-plev(i,2))
         deltap(i,nlay+1)=0.5*(plev(i,nlay)-plev(i,nlay+1))
         denom(i,nlay+1)=deltap(i,nlay+1)+kstar(i,nlay)
         alpha(i,nlay+1)=deltap(i,nlay+1)*q2(i,nlay+1)/denom(i,nlay+1)
         beta(i,nlay+1)=kstar(i,nlay)/denom(i,nlay+1)
      enddo

      do k=nlay,2,-1
         do i=1,ngrid
            denom(i,k)=deltap(i,k)+(1.-beta(i,k+1))*
     s      kstar(i,k)+kstar(i,k-1)
!   correction d'un bug 10 01 2001
            alpha(i,k)=(q2(i,k)*deltap(i,k)
     s      +kstar(i,k)*alpha(i,k+1))/denom(i,k)
            beta(i,k)=kstar(i,k-1)/denom(i,k)
         enddo
      enddo

!  Si on recalcule q2(1)
      if(1.eq.0) then
      do i=1,ngrid
         denom(i,1)=deltap(i,1)+(1-beta(i,2))*kstar(i,1)
         q2(i,1)=(q2(i,1)*deltap(i,1)
     s      +kstar(i,1)*alpha(i,2))/denom(i,1)
      enddo
      endif
!   sinon, on peut sauter cette boucle...

      do k=2,nlay+1
         do i=1,ngrid
            q2(i,k)=alpha(i,k)+beta(i,k)*q2(i,k-1)
         enddo
      enddo

      return
      end
      SUBROUTINE vdif_q2e(timestep,gravity,rconst,ngrid,nlay,
     &   plev,temp,kmy,q2)
      IMPLICIT NONE
!.......................................................................
! MARS
#include "callkeys.h"
!.......................................................................
!
! dt : pas de temps

      REAL plev(ngrid,nlay+1)
      REAL temp(ngrid,nlay)
      REAL timestep
      REAL gravity,rconst
      REAL kstar(ngrid,nlay+1),zz
      REAL kmy(ngrid,nlay+1)
      REAL q2(ngrid,nlay+1)
      REAL deltap(ngrid,nlay+1)
      REAL denom(ngrid,nlay+1),alpha(ngrid,nlay+1),beta(ngrid,nlay+1)
      INTEGER ngrid,nlay

      INTEGER i,k

      do k=1,nlay
         do i=1,ngrid
            zz=(plev(i,k)+plev(i,k+1))*gravity/(rconst*temp(i,k))
            kstar(i,k)=0.125*(kmy(i,k+1)+kmy(i,k))*zz*zz
     s      /(plev(i,k)-plev(i,k+1))*timestep
         enddo
      enddo

      do k=2,nlay
         do i=1,ngrid
            deltap(i,k)=0.5*(plev(i,k-1)-plev(i,k+1))
         enddo
      enddo
      do i=1,ngrid
         deltap(i,1)=0.5*(plev(i,1)-plev(i,2))
         deltap(i,nlay+1)=0.5*(plev(i,nlay)-plev(i,nlay+1))
      enddo

      do k=nlay,2,-1
         do i=1,ngrid
            q2(i,k)=q2(i,k)+
     s      ( kstar(i,k)*(q2(i,k+1)-q2(i,k))
     s       -kstar(i,k-1)*(q2(i,k)-q2(i,k-1)) )
     s      /deltap(i,k)
         enddo
      enddo

      do i=1,ngrid
         q2(i,1)=q2(i,1)+
     s   ( kstar(i,1)*(q2(i,2)-q2(i,1))
     s                                      )
     s   /deltap(i,1)
         q2(i,nlay+1)=q2(i,nlay+1)+
     s   ( 
     s    -kstar(i,nlay)*(q2(i,nlay+1)-q2(i,nlay)) )
     s   /deltap(i,nlay+1)
      enddo

      return
      end
