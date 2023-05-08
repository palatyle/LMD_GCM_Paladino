!
! $Id: thermcell_plume.F90 1403 2010-07-01 09:02:53Z fairhead $
!
      SUBROUTINE thermcell_plume(itap,ngrid,klev,ptimestep,ztv,zthl,po,zl,rhobarz,  &
     &           zlev,pplev,pphi,zpspsk,alim_star,alim_star_tot,  &
     &           lalim,f0,detr_star,entr_star,f_star,csc,ztva,  &
     &           ztla,zqla,zqta,zha,zw2,w_est,ztva_est,zqsatth,lmix,lmix_bis,linter &
     &           ,lev_out,lunout1,igout)

!--------------------------------------------------------------------------
!thermcell_plume: calcule les valeurs de qt, thetal et w dans l ascendance
!--------------------------------------------------------------------------

      IMPLICIT NONE

#include "YOMCST.h"
#include "YOETHF.h"
#include "FCTTRE.h"
#include "iniprint.h"
#include "thermcell.h"

      INTEGER itap
      INTEGER lunout1,igout
      INTEGER ngrid,klev
      REAL ptimestep
      REAL ztv(ngrid,klev)
      REAL zthl(ngrid,klev)
      REAL po(ngrid,klev)
      REAL zl(ngrid,klev)
      REAL rhobarz(ngrid,klev)
      REAL zlev(ngrid,klev+1)
      REAL pplev(ngrid,klev+1)
      REAL pphi(ngrid,klev)
      REAL zpspsk(ngrid,klev)
      REAL alim_star(ngrid,klev)
      REAL f0(ngrid)
      INTEGER lalim(ngrid)
      integer lev_out                           ! niveau pour les print
      real zcon2(ngrid)
    
      real alim_star_tot(ngrid)

      REAL ztva(ngrid,klev)
      REAL ztla(ngrid,klev)
      REAL zqla(ngrid,klev)
      REAL zqta(ngrid,klev)
      REAL zha(ngrid,klev)

      REAL detr_star(ngrid,klev)
      REAL coefc
      REAL entr_star(ngrid,klev)
      REAL detr(ngrid,klev)
      REAL entr(ngrid,klev)

      REAL csc(ngrid,klev)

      REAL zw2(ngrid,klev+1)
      REAL w_est(ngrid,klev+1)
      REAL f_star(ngrid,klev+1)
      REAL wa_moy(ngrid,klev+1)

      REAL ztva_est(ngrid,klev)
      REAL zqla_est(ngrid,klev)
      REAL zqsatth(ngrid,klev)
      REAL zta_est(ngrid,klev)
      REAL zdw2
      REAL zw2modif
      REAL zeps

      REAL linter(ngrid)
      INTEGER lmix(ngrid)
      INTEGER lmix_bis(ngrid)
      REAL    wmaxa(ngrid)

      INTEGER ig,l,k

      real zdz,zfact,zbuoy,zalpha,zdrag
      real zcor,zdelta,zcvm5,qlbef
      real Tbef,qsatbef
      real dqsat_dT,DT,num,denom
      REAL REPS,RLvCp,DDT0
      PARAMETER (DDT0=.01)
      logical Zsat
      LOGICAL active(ngrid),activetmp(ngrid)
      REAL fact_gamma,fact_epsilon,fact_gamma2
      REAL c2(ngrid,klev)
      REAL a1,m

      REAL zw2fact,expa
      Zsat=.false.
! Initialisation
      RLvCp = RLVTT/RCPD
     
      
         fact_epsilon=0.002
         a1=2./3.
         fact_gamma=0.9
         zfact=fact_gamma/(1+fact_gamma)
         fact_gamma2=zfact
         expa=0.
      

! Initialisations des variables reeles
if (1==1) then
      ztva(:,:)=ztv(:,:)
      ztva_est(:,:)=ztva(:,:)
      ztla(:,:)=zthl(:,:)
      zqta(:,:)=po(:,:)
      zha(:,:) = ztva(:,:)
else
      ztva(:,:)=0.
      ztva_est(:,:)=0.
      ztla(:,:)=0.
      zqta(:,:)=0.
      zha(:,:) =0.
endif

      zqla_est(:,:)=0.
      zqsatth(:,:)=0.
      zqla(:,:)=0.
      detr_star(:,:)=0.
      entr_star(:,:)=0.
      alim_star(:,:)=0.
      alim_star_tot(:)=0.
      csc(:,:)=0.
      detr(:,:)=0.
      entr(:,:)=0.
      zw2(:,:)=0.
      w_est(:,:)=0.
      f_star(:,:)=0.
      wa_moy(:,:)=0.
      linter(:)=1.
      linter(:)=1.

! Initialisation des variables entieres
      lmix(:)=1
      lmix_bis(:)=2
      wmaxa(:)=0.
      lalim(:)=1

!-------------------------------------------------------------------------
! On ne considere comme actif que les colonnes dont les deux premieres
! couches sont instables.
!-------------------------------------------------------------------------
      active(:)=ztv(:,1)>ztv(:,2)

!-------------------------------------------------------------------------
! Definition de l'alimentation a l'origine dans thermcell_init
!-------------------------------------------------------------------------
      do l=1,klev-1
         do ig=1,ngrid
            if (ztv(ig,l)> ztv(ig,l+1) .and. ztv(ig,1)>=ztv(ig,l) ) then
               alim_star(ig,l)=MAX((ztv(ig,l)-ztv(ig,l+1)),0.)  &
     &                       *sqrt(zlev(ig,l+1)) 
               lalim(:)=l+1
               alim_star_tot(ig)=alim_star_tot(ig)+alim_star(ig,l)
            endif
         enddo
      enddo
      do l=1,klev
         do ig=1,ngrid 
            if (alim_star_tot(ig) > 1.e-10 ) then
               alim_star(ig,l)=alim_star(ig,l)/alim_star_tot(ig)
            endif
         enddo
      enddo
      alim_star_tot(:)=1.


!------------------------------------------------------------------------------
! Calcul dans la premiere couche
! On decide dans cette version que le thermique n'est actif que si la premiere
! couche est instable.
! Pourrait etre change si on veut que le thermiques puisse se déclencher
! dans une couche l>1
!------------------------------------------------------------------------------
do ig=1,ngrid
! Le panache va prendre au debut les caracteristiques de l'air contenu
! dans cette couche.
    if (active(ig)) then
    ztla(ig,1)=zthl(ig,1) 
    zqta(ig,1)=po(ig,1)
    zqla(ig,1)=zl(ig,1)
!cr: attention, prise en compte de f*(1)=1
    f_star(ig,2)=alim_star(ig,1)
    zw2(ig,2)=2.*RG*(ztv(ig,1)-ztv(ig,2))/ztv(ig,2)  &
&                     *(zlev(ig,2)-zlev(ig,1))  &
&                     *0.4*pphi(ig,1)/(pphi(ig,2)-pphi(ig,1))
    w_est(ig,2)=zw2(ig,2)
    endif
enddo
!

!==============================================================================
!boucle de calcul de la vitesse verticale dans le thermique
!==============================================================================
do l=2,klev-1
!==============================================================================


! On decide si le thermique est encore actif ou non
! AFaire : Il faut sans doute ajouter entr_star a alim_star dans ce test
    do ig=1,ngrid
       active(ig)=active(ig) &
&                 .and. zw2(ig,l)>1.e-10 &
&                 .and. f_star(ig,l)+alim_star(ig,l)>1.e-10
    enddo



! Premier calcul de la vitesse verticale a partir de la temperature
! potentielle virtuelle
!     if (1.eq.1) then
!         w_est(ig,3)=zw2(ig,2)* &
!    &      ((f_star(ig,2))**2) &
!    &      /(f_star(ig,2)+alim_star(ig,2))**2+ &
!    &      2.*RG*(ztva(ig,2)-ztv(ig,2))/ztv(ig,2) &
!    &      *(zlev(ig,3)-zlev(ig,2))
!      endif


!---------------------------------------------------------------------------
! calcul des proprietes thermodynamiques et de la vitesse de la couche l
! sans tenir compte du detrainement et de l'entrainement dans cette
! couche
! Ici encore, on doit pouvoir ajouter entr_star (qui peut etre calculer
! avant) a l'alimentation pour avoir un calcul plus propre
!---------------------------------------------------------------------------

     call thermcell_condens(ngrid,active, &
&          zpspsk(:,l),pplev(:,l),ztla(:,l-1),zqta(:,l-1),zqla_est(:,l))

    do ig=1,ngrid
        if(active(ig)) then
        ztva_est(ig,l) = ztla(ig,l-1)*zpspsk(ig,l)+RLvCp*zqla_est(ig,l)
        zta_est(ig,l)=ztva_est(ig,l)
        ztva_est(ig,l) = ztva_est(ig,l)/zpspsk(ig,l)
        ztva_est(ig,l) = ztva_est(ig,l)*(1.+RETV*(zqta(ig,l-1)  &
     &      -zqla_est(ig,l))-zqla_est(ig,l))

         if (1.eq.0) then 
!calcul de w_est sans prendre en compte le drag 
             w_est(ig,l+1)=zw2(ig,l)*  &
     &                   ((f_star(ig,l))**2)  &
     &                   /(f_star(ig,l)+alim_star(ig,l))**2+  &
     &                   2.*RG*(ztva_est(ig,l)-ztv(ig,l))/ztv(ig,l)  &
     &                   *(zlev(ig,l+1)-zlev(ig,l))
         else

           zdz=zlev(ig,l+1)-zlev(ig,l)
           zalpha=f0(ig)*f_star(ig,l)/sqrt(w_est(ig,l))/rhobarz(ig,l)
           zbuoy=RG*(ztva_est(ig,l)-ztv(ig,l))/ztv(ig,l)
           zdrag=fact_epsilon/(zalpha**expa)
           zw2fact=zbuoy/zdrag*a1
           w_est(ig,l+1)=(w_est(ig,l)-zw2fact)*exp(-2.*zdrag/(1+fact_gamma)*zdz) &
      &    +zw2fact

         endif

             if (w_est(ig,l+1).lt.0.) then
                w_est(ig,l+1)=zw2(ig,l)
             endif
       endif
    enddo

!-------------------------------------------------
!calcul des taux d'entrainement et de detrainement
!-------------------------------------------------

     do ig=1,ngrid
        if (active(ig)) then
          zdz=zlev(ig,l+1)-zlev(ig,l)
          zbuoy=RG*(ztva_est(ig,l)-ztv(ig,l))/ztv(ig,l)
 
! estimation de la fraction couverte par les thermiques
          zalpha=f0(ig)*f_star(ig,l)/sqrt(w_est(ig,l))/rhobarz(ig,l)

!calcul de la soumission papier 
! Calcul  du taux d'entrainement entr_star (epsilon)
           entr_star(ig,l)=f_star(ig,l)*zdz * (  zfact * MAX(0.,  &     
     &     a1*zbuoy/w_est(ig,l+1) &
     &     - fact_epsilon/zalpha**expa  ) &
     &     +0. )
           
!calcul du taux de detrainment (delta)
!           detr_star(ig,l)=f_star(ig,l)*zdz * (                           &
!     &      MAX(1.e-3, &
!     &      -fact_gamma2*a1*zbuoy/w_est(ig,l+1)        &
!     &      +0.01*(max(zqta(ig,l-1)-po(ig,l),0.)/(po(ig,l))/(w_est(ig,l+1)))**0.5    &    
!     &     +0. ))

          m=0.5

          detr_star(ig,l)=1.*f_star(ig,l)*zdz *                    &
    &     MAX(5.e-4,-fact_gamma2*a1*(1./w_est(ig,l+1))*((1.-(1.-m)/(1.+70*zqta(ig,l-1)))*zbuoy        &
    &     -40*(1.-m)*(max(zqta(ig,l-1)-po(ig,l),0.))/(1.+70*zqta(ig,l-1)) )   )

!           detr_star(ig,l)=f_star(ig,l)*zdz * (                           &
!     &      MAX(0.0, &
!     &      -fact_gamma2*a1*zbuoy/w_est(ig,l+1)        &
!     &      +20*(max(zqta(ig,l-1)-po(ig,l),0.))**1*(zalpha/w_est(ig,l+1))**0.5    &    
!     &     +0. ))


! En dessous de lalim, on prend le max de alim_star et entr_star pour
! alim_star et 0 sinon
        if (l.lt.lalim(ig)) then
          alim_star(ig,l)=max(alim_star(ig,l),entr_star(ig,l))
          entr_star(ig,l)=0.
        endif

!attention test
!        if (detr_star(ig,l).gt.(f_star(ig,l)+alim_star(ig,l)+entr_star(ig,l))) then       
!            detr_star(ig,l)=f_star(ig,l)+alim_star(ig,l)+entr_star(ig,l)
!        endif
! Calcul du flux montant normalise
      f_star(ig,l+1)=f_star(ig,l)+alim_star(ig,l)+entr_star(ig,l)  &
     &              -detr_star(ig,l)

      endif
   enddo

!----------------------------------------------------------------------------
!calcul de la vitesse verticale en melangeant Tl et qt du thermique
!---------------------------------------------------------------------------
   activetmp(:)=active(:) .and. f_star(:,l+1)>1.e-10
   do ig=1,ngrid
       if (activetmp(ig)) then 
           Zsat=.false.
           ztla(ig,l)=(f_star(ig,l)*ztla(ig,l-1)+  &
     &            (alim_star(ig,l)+entr_star(ig,l))*zthl(ig,l))  &
     &            /(f_star(ig,l+1)+detr_star(ig,l))
           zqta(ig,l)=(f_star(ig,l)*zqta(ig,l-1)+  &
     &            (alim_star(ig,l)+entr_star(ig,l))*po(ig,l))  &
     &            /(f_star(ig,l+1)+detr_star(ig,l))

        endif
    enddo

   call thermcell_condens(ngrid,activetmp,zpspsk(:,l),pplev(:,l),ztla(:,l),zqta(:,l),zqla(:,l))


   do ig=1,ngrid
      if (activetmp(ig)) then
        if (prt_level.ge.20) print*,'coucou calcul detr 4512: ig, l', ig, l
! on ecrit de maniere conservative (sat ou non)
!          T = Tl +Lv/Cp ql
           ztva(ig,l) = ztla(ig,l)*zpspsk(ig,l)+RLvCp*zqla(ig,l)
           ztva(ig,l) = ztva(ig,l)/zpspsk(ig,l)
!on rajoute le calcul de zha pour diagnostiques (temp potentielle)
           zha(ig,l) = ztva(ig,l)
           ztva(ig,l) = ztva(ig,l)*(1.+RETV*(zqta(ig,l)  &
     &              -zqla(ig,l))-zqla(ig,l))

!on ecrit zqsat 
           zqsatth(ig,l)=qsatbef  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          zw2(ig,l+1)=&
!     &                 zw2(ig,l)*(1-fact_epsilon/(1.+fact_gamma)*2.*(zlev(ig,l+1)-zlev(ig,l))) &
!     &                 +2.*RG*(ztva(ig,l)-ztv(ig,l))/ztv(ig,l)  &
!     &                 *1./(1.+fact_gamma) &
!     &                 *(zlev(ig,l+1)-zlev(ig,l))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! La meme en plus modulaire :
           zbuoy=RG*(ztva(ig,l)-ztv(ig,l))/ztv(ig,l)
           zdz=zlev(ig,l+1)-zlev(ig,l)


           zeps=(entr_star(ig,l)+alim_star(ig,l))/(f_star(ig,l)*zdz)

if (1==0) then
           zw2modif=zw2(ig,l)*(1-fact_epsilon/(1.+fact_gamma)*2.*zdz)
           zdw2=2.*zbuoy/(1.+fact_gamma)*zdz
           zw2(ig,l+1)=zw2modif+zdw2
else
           zdrag=fact_epsilon/(zalpha**expa)
           zw2fact=zbuoy/zdrag*a1
           zw2(ig,l+1)=(zw2(ig,l)-zw2fact)*exp(-2.*zdrag/(1+fact_gamma)*zdz) &
      &    +zw2fact


endif

      endif
   enddo

        if (prt_level.ge.20) print*,'coucou calcul detr 460: ig, l',ig, l
!
!---------------------------------------------------------------------------
!initialisations pour le calcul de la hauteur du thermique, de l'inversion et de la vitesse verticale max 
!---------------------------------------------------------------------------

   do ig=1,ngrid
            if (zw2(ig,l+1)>0. .and. zw2(ig,l+1).lt.1.e-10) then
!               stop'On tombe sur le cas particulier de thermcell_dry'
                print*,'On tombe sur le cas particulier de thermcell_plume'
                zw2(ig,l+1)=0.
                linter(ig)=l+1
            endif

        if (zw2(ig,l+1).lt.0.) then 
           linter(ig)=(l*(zw2(ig,l+1)-zw2(ig,l))  &
     &               -zw2(ig,l))/(zw2(ig,l+1)-zw2(ig,l))
           zw2(ig,l+1)=0.
        endif

           wa_moy(ig,l+1)=sqrt(zw2(ig,l+1)) 

        if (wa_moy(ig,l+1).gt.wmaxa(ig)) then
!   lmix est le niveau de la couche ou w (wa_moy) est maximum
!on rajoute le calcul de lmix_bis
            if (zqla(ig,l).lt.1.e-10) then
               lmix_bis(ig)=l+1
            endif
            lmix(ig)=l+1
            wmaxa(ig)=wa_moy(ig,l+1)
        endif
   enddo

!=========================================================================
! FIN DE LA BOUCLE VERTICALE
      enddo
!=========================================================================

!on recalcule alim_star_tot
       do ig=1,ngrid
          alim_star_tot(ig)=0.
       enddo
       do ig=1,ngrid
          do l=1,lalim(ig)-1
          alim_star_tot(ig)=alim_star_tot(ig)+alim_star(ig,l)
          enddo
       enddo
       

        if (prt_level.ge.20) print*,'coucou calcul detr 470: ig, l', ig, l


      return 
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE thermcellV1_plume(itap,ngrid,klev,ptimestep,ztv,zthl,po,zl,rhobarz,  &
&           zlev,pplev,pphi,zpspsk,alim_star,alim_star_tot,  &
&           lalim,f0,detr_star,entr_star,f_star,csc,ztva,  &
&           ztla,zqla,zqta,zha,zw2,w_est,ztva_est,zqsatth,lmix,lmix_bis,linter &
&           ,lev_out,lunout1,igout)

!--------------------------------------------------------------------------
!thermcell_plume: calcule les valeurs de qt, thetal et w dans l ascendance
! Version conforme a l'article de Rio et al. 2010.
! Code ecrit par Catherine Rio, Arnaud Jam et Frederic Hourdin
!--------------------------------------------------------------------------

      IMPLICIT NONE

#include "YOMCST.h"
#include "YOETHF.h"
#include "FCTTRE.h"
#include "iniprint.h"
#include "thermcell.h"

      INTEGER itap
      INTEGER lunout1,igout
      INTEGER ngrid,klev
      REAL ptimestep
      REAL ztv(ngrid,klev)
      REAL zthl(ngrid,klev)
      REAL po(ngrid,klev)
      REAL zl(ngrid,klev)
      REAL rhobarz(ngrid,klev)
      REAL zlev(ngrid,klev+1)
      REAL pplev(ngrid,klev+1)
      REAL pphi(ngrid,klev)
      REAL zpspsk(ngrid,klev)
      REAL alim_star(ngrid,klev)
      REAL f0(ngrid)
      INTEGER lalim(ngrid)
      integer lev_out                           ! niveau pour les print
    
      real alim_star_tot(ngrid)

      REAL ztva(ngrid,klev)
      REAL ztla(ngrid,klev)
      REAL zqla(ngrid,klev)
      REAL zqta(ngrid,klev)
      REAL zha(ngrid,klev)

      REAL detr_star(ngrid,klev)
      REAL coefc
      REAL entr_star(ngrid,klev)
      REAL detr(ngrid,klev)
      REAL entr(ngrid,klev)

      REAL csc(ngrid,klev)

      REAL zw2(ngrid,klev+1)
      REAL w_est(ngrid,klev+1)
      REAL f_star(ngrid,klev+1)
      REAL wa_moy(ngrid,klev+1)

      REAL ztva_est(ngrid,klev)
      REAL zqla_est(ngrid,klev)
      REAL zqsatth(ngrid,klev)
      REAL zta_est(ngrid,klev)
      REAL ztemp(ngrid),zqsat(ngrid)
      REAL zdw2
      REAL zw2modif
      REAL zw2fact
      REAL zeps(ngrid,klev)

      REAL linter(ngrid)
      INTEGER lmix(ngrid)
      INTEGER lmix_bis(ngrid)
      REAL    wmaxa(ngrid)

      INTEGER ig,l,k

      real zdz,zbuoy(ngrid,klev),zalpha,gamma(ngrid,klev),zdqt(ngrid,klev),zw2m
      real zbuoybis
      real zcor,zdelta,zcvm5,qlbef,zdz2
      real betalpha,zbetalpha
      real eps, afact
      REAL REPS,RLvCp,DDT0
      PARAMETER (DDT0=.01)
      logical Zsat
      LOGICAL active(ngrid),activetmp(ngrid)
      REAL fact_gamma,fact_epsilon,fact_gamma2,fact_epsilon2
      REAL c2(ngrid,klev)
      Zsat=.false.
! Initialisation

      RLvCp = RLVTT/RCPD
      fact_epsilon=0.002
      betalpha=0.9 
      afact=2./3.            

      zbetalpha=betalpha/(1.+betalpha)


! Initialisations des variables reeles
if (1==0) then
      ztva(:,:)=ztv(:,:)
      ztva_est(:,:)=ztva(:,:)
      ztla(:,:)=zthl(:,:)
      zqta(:,:)=po(:,:)
      zha(:,:) = ztva(:,:)
else
      ztva(:,:)=0.
      ztva_est(:,:)=0.
      ztla(:,:)=0.
      zqta(:,:)=0.
      zha(:,:) =0.
endif

      zqla_est(:,:)=0.
      zqsatth(:,:)=0.
      zqla(:,:)=0.
      detr_star(:,:)=0.
      entr_star(:,:)=0.
      alim_star(:,:)=0.
      alim_star_tot(:)=0.
      csc(:,:)=0.
      detr(:,:)=0.
      entr(:,:)=0.
      zw2(:,:)=0.
      zbuoy(:,:)=0.
      gamma(:,:)=0.
      zeps(:,:)=0.
      w_est(:,:)=0.
      f_star(:,:)=0.
      wa_moy(:,:)=0.
      linter(:)=1.
!     linter(:)=1.
! Initialisation des variables entieres
      lmix(:)=1
      lmix_bis(:)=2
      wmaxa(:)=0.
      lalim(:)=1


!-------------------------------------------------------------------------
! On ne considere comme actif que les colonnes dont les deux premieres
! couches sont instables.
!-------------------------------------------------------------------------
      active(:)=ztv(:,1)>ztv(:,2)

!-------------------------------------------------------------------------
! Definition de l'alimentation a l'origine dans thermcell_init
!-------------------------------------------------------------------------
      do l=1,klev-1
         do ig=1,ngrid
            if (ztv(ig,l)> ztv(ig,l+1) .and. ztv(ig,1)>=ztv(ig,l) ) then
               alim_star(ig,l)=MAX((ztv(ig,l)-ztv(ig,l+1)),0.)  &
     &                       *sqrt(zlev(ig,l+1)) 
               lalim(ig)=l+1
               alim_star_tot(ig)=alim_star_tot(ig)+alim_star(ig,l)
            endif
         enddo
      enddo
      do l=1,klev
         do ig=1,ngrid 
            if (alim_star_tot(ig) > 1.e-10 ) then
               alim_star(ig,l)=alim_star(ig,l)/alim_star_tot(ig)
            endif
         enddo
      enddo
      alim_star_tot(:)=1.



!------------------------------------------------------------------------------
! Calcul dans la premiere couche
! On decide dans cette version que le thermique n'est actif que si la premiere
! couche est instable.
! Pourrait etre change si on veut que le thermiques puisse se déclencher
! dans une couche l>1
!------------------------------------------------------------------------------
do ig=1,ngrid
! Le panache va prendre au debut les caracteristiques de l'air contenu
! dans cette couche.
    if (active(ig)) then
    ztla(ig,1)=zthl(ig,1) 
    zqta(ig,1)=po(ig,1)
    zqla(ig,1)=zl(ig,1)
!cr: attention, prise en compte de f*(1)=1
    f_star(ig,2)=alim_star(ig,1)
    zw2(ig,2)=2.*RG*(ztv(ig,1)-ztv(ig,2))/ztv(ig,2)  &
&                     *(zlev(ig,2)-zlev(ig,1))  &
&                     *0.4*pphi(ig,1)/(pphi(ig,2)-pphi(ig,1))
    w_est(ig,2)=zw2(ig,2)
    endif
enddo
!

!==============================================================================
!boucle de calcul de la vitesse verticale dans le thermique
!==============================================================================
do l=2,klev-1
!==============================================================================


! On decide si le thermique est encore actif ou non
! AFaire : Il faut sans doute ajouter entr_star a alim_star dans ce test
    do ig=1,ngrid
       active(ig)=active(ig) &
&                 .and. zw2(ig,l)>1.e-10 &
&                 .and. f_star(ig,l)+alim_star(ig,l)>1.e-10
    enddo



!---------------------------------------------------------------------------
! calcul des proprietes thermodynamiques et de la vitesse de la couche l
! sans tenir compte du detrainement et de l'entrainement dans cette
! couche
! C'est a dire qu'on suppose 
! ztla(l)=ztla(l-1) et zqta(l)=zqta(l-1)
! Ici encore, on doit pouvoir ajouter entr_star (qui peut etre calculer
! avant) a l'alimentation pour avoir un calcul plus propre
!---------------------------------------------------------------------------

   ztemp(:)=zpspsk(:,l)*ztla(:,l-1)
   call thermcell_qsat(ngrid,active,pplev(:,l),ztemp,zqta(:,l-1),zqsat(:))

    do ig=1,ngrid 
!       print*,'active',active(ig),ig,l
        if(active(ig)) then 
        zqla_est(ig,l)=max(0.,zqta(ig,l-1)-zqsat(ig))
        ztva_est(ig,l) = ztla(ig,l-1)*zpspsk(ig,l)+RLvCp*zqla_est(ig,l)
        zta_est(ig,l)=ztva_est(ig,l)
        ztva_est(ig,l) = ztva_est(ig,l)/zpspsk(ig,l)
        ztva_est(ig,l) = ztva_est(ig,l)*(1.+RETV*(zqta(ig,l-1)  &
     &      -zqla_est(ig,l))-zqla_est(ig,l))

!------------------------------------------------
!AJAM:nouveau calcul de w�  
!------------------------------------------------
              zdz=zlev(ig,l+1)-zlev(ig,l)
              zbuoy(ig,l)=RG*(ztva_est(ig,l)-ztv(ig,l))/ztv(ig,l)

              zw2fact=fact_epsilon*2.*zdz/(1.+betalpha)
              zdw2=(afact)*zbuoy(ig,l)/(fact_epsilon)
              w_est(ig,l+1)=Max(0.0001,exp(-zw2fact)*(w_est(ig,l)-zdw2)+zdw2)
 

             if (w_est(ig,l+1).lt.0.) then
                w_est(ig,l+1)=zw2(ig,l)
             endif
       endif
    enddo


!-------------------------------------------------
!calcul des taux d'entrainement et de detrainement
!-------------------------------------------------

     do ig=1,ngrid
        if (active(ig)) then

          zw2m=max(0.5*(w_est(ig,l)+w_est(ig,l+1)),0.1)
          zw2m=w_est(ig,l+1)
          zdz=zlev(ig,l+1)-zlev(ig,l)
          zbuoy(ig,l)=RG*(ztva_est(ig,l)-ztv(ig,l))/ztv(ig,l)
!          zbuoybis=zbuoy(ig,l)+RG*0.1/300.
          zbuoybis=zbuoy(ig,l)
          zalpha=f0(ig)*f_star(ig,l)/sqrt(w_est(ig,l+1))/rhobarz(ig,l)
          zdqt(ig,l)=max(zqta(ig,l-1)-po(ig,l),0.)/po(ig,l)

          
          entr_star(ig,l)=f_star(ig,l)*zdz*  zbetalpha*MAX(0.,  &
    &     afact*zbuoybis/zw2m - fact_epsilon )


          detr_star(ig,l)=f_star(ig,l)*zdz                        &
    &     *MAX(1.e-3, -afact*zbetalpha*zbuoy(ig,l)/zw2m          &
    &     + 0.012*(zdqt(ig,l)/zw2m)**0.5 )
          
! En dessous de lalim, on prend le max de alim_star et entr_star pour
! alim_star et 0 sinon
        if (l.lt.lalim(ig)) then
          alim_star(ig,l)=max(alim_star(ig,l),entr_star(ig,l))
          entr_star(ig,l)=0.
        endif

! Calcul du flux montant normalise
      f_star(ig,l+1)=f_star(ig,l)+alim_star(ig,l)+entr_star(ig,l)  &
     &              -detr_star(ig,l)

      endif
   enddo


!----------------------------------------------------------------------------
!calcul de la vitesse verticale en melangeant Tl et qt du thermique
!---------------------------------------------------------------------------
   activetmp(:)=active(:) .and. f_star(:,l+1)>1.e-10
   do ig=1,ngrid
       if (activetmp(ig)) then 
           Zsat=.false.
           ztla(ig,l)=(f_star(ig,l)*ztla(ig,l-1)+  &
     &            (alim_star(ig,l)+entr_star(ig,l))*zthl(ig,l))  &
     &            /(f_star(ig,l+1)+detr_star(ig,l))
           zqta(ig,l)=(f_star(ig,l)*zqta(ig,l-1)+  &
     &            (alim_star(ig,l)+entr_star(ig,l))*po(ig,l))  &
     &            /(f_star(ig,l+1)+detr_star(ig,l))

        endif
    enddo

   ztemp(:)=zpspsk(:,l)*ztla(:,l)
   call thermcell_qsat(ngrid,activetmp,pplev(:,l),ztemp,zqta(:,l),zqsatth(:,l))

   do ig=1,ngrid
      if (activetmp(ig)) then
        if (prt_level.ge.20) print*,'coucou calcul detr 4512: ig, l', ig, l
! on ecrit de maniere conservative (sat ou non)
!          T = Tl +Lv/Cp ql
           zqla(ig,l)=max(0.,zqta(ig,l)-zqsatth(ig,l))
           ztva(ig,l) = ztla(ig,l)*zpspsk(ig,l)+RLvCp*zqla(ig,l)
           ztva(ig,l) = ztva(ig,l)/zpspsk(ig,l)
!on rajoute le calcul de zha pour diagnostiques (temp potentielle)
           zha(ig,l) = ztva(ig,l)
           ztva(ig,l) = ztva(ig,l)*(1.+RETV*(zqta(ig,l)  &
     &              -zqla(ig,l))-zqla(ig,l))
           zbuoy(ig,l)=RG*(ztva(ig,l)-ztv(ig,l))/ztv(ig,l)
           zdz=zlev(ig,l+1)-zlev(ig,l)
           zeps(ig,l)=(entr_star(ig,l)+alim_star(ig,l))/(f_star(ig,l)*zdz)

            zw2fact=fact_epsilon*2.*zdz/(1.+betalpha)
            zdw2=afact*zbuoy(ig,l)/(fact_epsilon)
            zw2(ig,l+1)=Max(0.0001,exp(-zw2fact)*(zw2(ig,l)-zdw2)+zdw2) 
      endif
   enddo

        if (prt_level.ge.20) print*,'coucou calcul detr 460: ig, l',ig, l
!
!---------------------------------------------------------------------------
!initialisations pour le calcul de la hauteur du thermique, de l'inversion et de la vitesse verticale max 
!---------------------------------------------------------------------------

   do ig=1,ngrid
            if (zw2(ig,l+1)>0. .and. zw2(ig,l+1).lt.1.e-10) then
!               stop'On tombe sur le cas particulier de thermcell_dry'
                print*,'On tombe sur le cas particulier de thermcell_plume'
                zw2(ig,l+1)=0.
                linter(ig)=l+1
            endif

        if (zw2(ig,l+1).lt.0.) then 
           linter(ig)=(l*(zw2(ig,l+1)-zw2(ig,l))  &
     &               -zw2(ig,l))/(zw2(ig,l+1)-zw2(ig,l))
           zw2(ig,l+1)=0.
        endif

           wa_moy(ig,l+1)=sqrt(zw2(ig,l+1)) 

        if (wa_moy(ig,l+1).gt.wmaxa(ig)) then
!   lmix est le niveau de la couche ou w (wa_moy) est maximum
!on rajoute le calcul de lmix_bis
            if (zqla(ig,l).lt.1.e-10) then
               lmix_bis(ig)=l+1
            endif
            lmix(ig)=l+1
            wmaxa(ig)=wa_moy(ig,l+1)
        endif
   enddo

!=========================================================================
! FIN DE LA BOUCLE VERTICALE
      enddo
!=========================================================================

!on recalcule alim_star_tot
       do ig=1,ngrid
          alim_star_tot(ig)=0.
       enddo
       do ig=1,ngrid
          do l=1,lalim(ig)-1
          alim_star_tot(ig)=alim_star_tot(ig)+alim_star(ig,l)
          enddo
       enddo
       

        if (prt_level.ge.20) print*,'coucou calcul detr 470: ig, l', ig, l

#undef wrgrads_thermcell
#ifdef wrgrads_thermcell
         call wrgradsfi(1,klev,entr_star(igout,1:klev),'esta      ','esta      ')
         call wrgradsfi(1,klev,detr_star(igout,1:klev),'dsta      ','dsta      ')
         call wrgradsfi(1,klev,zbuoy(igout,1:klev),'buoy      ','buoy      ')
         call wrgradsfi(1,klev,zdqt(igout,1:klev),'dqt      ','dqt      ')
         call wrgradsfi(1,klev,w_est(igout,1:klev),'w_est     ','w_est     ')
         call wrgradsfi(1,klev,w_est(igout,2:klev+1),'w_es2     ','w_es2     ')
         call wrgradsfi(1,klev,zw2(igout,1:klev),'zw2A      ','zw2A      ')
#endif


     return 
     end

