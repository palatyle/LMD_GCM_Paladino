         subroutine snuages(ngrid,tab1,tab2,tab3,tab4,tab5,
     &   x1,xnz,xnrad,ihor,fluxi,taused,precip)

!*
!
!1   c quantite de noyaux  de la grille de rayon r a l'altitude z 
!2   c quantite de glace 1 de la grille de rayon r a l'altitude z 
!3   c quantite de glace 2 de la grille de rayon r a l'altitude z 
!4   c quantite de glace 3 de la grille de rayon r a l'altitude z 
!5   c quantite d'aerosols de la grille de rayon r a l'altitude z 
c
c
!x1  timestep.

c    Notes : taused n'est actuellement pas calculé !
c    taused = dz/vitesse(z)
c
c
*--------------------------------------------------------------*
*                                                               *
*             ENTRE 0 ET 1000 KILOMETRES                        *     
*                                                               * 
*     la dimension fractale est en tableau, attention au        *
*     raccordement entre le regime moleculaire et le regime     *
*     fluide                                                    *
*                                                               *
*     Modele microphysique:    Cabane et al.,1992 /             *
*     Modele version fractale: Cabane et al.,1993 /             *
*                                                               *
*--------------------------------------------------------------*
*  VERSION DU 2 JUIN 1993  --- AUT 1994 --- 11/04/96
*
*  changer: altitude de production z0=/taux de production ctot= 
*         : la charge/micron, ne
*         : df(h),rf... 
*  raccordement aknc 
*
*  declaration des blocs communs
*------------------------------

         use dimphy          
         IMPLICIT NONE
#include  "dimensions.h"
#include  "microtab.h"
#include  "varmuphy.h"

         common/x2ctps/li,lf,dt
         common/x2con/c,c1,c2,c3,caer
         common/cldpart/rad,mmu
         common/x2frac/rfg,dfg
         common/x2effets/ xsaison

  
*  declaration des variables communes
*  ----------------------------------

         integer xnz,xnrad,ngrid
         integer i,j,k,ihor,ibidon
         integer ihormx,imx
         integer li,lf
         real xsaison
         real dt
         real c(nz,nrad,2)
         real c1(nz,nrad,2),c2(nz,nrad,2)
         real c3(nz,nrad,2)
         real caer(nz,nrad,2)
         real rad(nz,nrad), mmu(nz,nrad)
         real rtemp(nz,nrad),mmutemp(nz,nrad)
         real knu,knu2,w
         real rfg(nz),dfg(nz,nrad)
         real ddt,dtau,vitesse2
         real vmax,vmin,rmax,rmin
         real fluxi(ngrid,nz,3)
         real taused(ngrid,nz)
         real precip(ngrid,5)
         real rmx


*   variables internes
*   ------------------

         integer h,ti,itime,icoal
         real tab1(nz,nrad),tab2(nz,nrad),
     &        tab3(nz,nrad),tab4(nz,nrad)
         real tab5(nz,nrad)
         real x1,xvolume,xmasse,xnoyaux 
         real knd1,knd2,knd3,knd4,knd5,knd6

         real dice1,dice2,dice3,dice4

          data itime/0/
          data icoal/0/
          save itime,icoal
          save vmax,vmin,rmax,rmin


*  controles
*  ---------

         if (nrad.ne.xnrad)  stop 'nrad.ne.xnrad'
         if (nz.ne.xnz)      stop 'nz.ne.xnz'

         do i=1,nz
           do j=1,nrad 
             c(i,j,1)=tab1(i,j)
             c(i,j,2)=0.0 
             c1(i,j,1)=tab2(i,j)
             c1(i,j,2)=0.0
             c2(i,j,1)=tab3(i,j)
             c2(i,j,2)=0.0 
             c3(i,j,1)=tab4(i,j)
             c3(i,j,2)=0.0 
             caer(i,j,1)=tab5(i,j)
             caer(i,j,2)=0.0 
           enddo 
         enddo 


*  initialisation
*  --------------

*  Parametres physiques de Titan 
  
         if (itime.eq.0) then 
           call init2
           print*,'init2 dans snuages.F'
           itime=1
         endif

         aknc=2.92         !<--------Df=3

*  Propriete des gouttes nuageuses/ calculee ici pour la sedimentation
*  r(i,j)=r  (cst pour la colonne). Cela evite les accumulations dans
*  certains niveaux du a des differences de vitesse de sedimentation...

*  ATTENTION, CES DEFINITIONS SONT INTERNES/ VOIR MUPHYS.F POUR DEFINITIONS
*  EXTERNES (UTILISEES DANS OPTCV ET OPTCI, PAR EX.)

         do i=1,nz
           xnoyaux=0.
           xvolume=0.
           xmasse=0.
           do j=1,nrad
             dfg(i,j)=3.00
             rfg(j)=6.6e-8     !<--- arbitraire pour df=3
             xnoyaux=xnoyaux+c(i,j,1)
             xvolume=xvolume+c1(i,j,1)+c2(i,j,1)+c3(i,j,1)+
     &       v_e(j)*c(i,j,1)
             xmasse =xmasse+
     &       c1(i,j,1)*rhoi_ch4+
     &       c2(i,j,1)*rhoi_c2h6+
     &       c3(i,j,1)*rhoi_c2h2+ 
     &       v_e(j)*c(i,j,1)*rhol
           enddo

           do j=1,nrad

             if (xnoyaux .le. 1.e-25 .or. xvolume.eq.0.) then   !  utile ?
               dfg(i,j)=3.00
               rtemp(i,j) =  1.e-9
               mmutemp(i,j)= rhol     !on prend la masse vol des aerosols
             else
c         Mais pourquoi se compliquer la vie alors que de toute facon on 
c         reforce la masse volumique a son veritable calcul ^^
c
               if ((c1(i,j,1)+c2(i,j,1)+c3(i,j,1))/xvolume.le.0.1) then !  glace / total  
                 dfg(i,j)=3.00
                 rtemp(i,j) = ( (xvolume/xnoyaux)*0.75/pi)**(1./3.)
c                 mmutemp(i,j) = 1000.
                 mmutemp(i,j)= xmasse/xvolume 
                 if(rtemp(i,j).gt.3.e-4) rtemp(i,j)=3.e-4
               else
                 dfg(i,j)=3.00
                 rtemp(i,j) = ( (xvolume/xnoyaux)*0.75/pi)**(1./3.)
                 mmutemp(i,j)= xmasse/xvolume 
                 if(rtemp(i,j).gt.3.e-4) rtemp(i,j)=3.e-4
               endif
             endif
           enddo
         enddo

         do i=1,nz
           do j=1,nrad
             rad(i,j)=rtemp(i,1)   ! mm valeur qqsoit j
             mmu(i,j)=mmutemp(i,1)
           enddo
         enddo


*  Coefficients de coagulation

         if (icoal.eq.1) call calcoag2             ! 2 <-- 1
         vmin=1.e+6
         vmax=0.
         rmin=1.e+6
         rmax=0.
         ihormx=1
         imx=1


*  choix interne du temps d iteration
*  ----------------------------------

         dt=x1

*  iteration du modele sur le temps
*  ---------------------------------

         li=1  
         lf=2

*  li=1  lf=2

         call sedifn_fast(ihor,dice1,dice2,dice3,dice4) ! 1 <-- 1
!        call sedifn                                 ! 2 <-- 1

         li=3-li         ! li devient 2
         lf=3-lf         ! lf devient 1

         if (icoal.eq.1) then 
*  li=2  lf=1
           call coagul2(ihor)                          ! 1 <-- 2

           li=3-li         ! li devient 1
           lf=3-lf         ! lf devient 2

         endif

         do i=1,nz
           vmin=vitesse2(i,1,1)
           do j=1,nrad
             fluxi(ihor,nz+1-i,1)=fluxi(ihor,nz+1-i,1)
     &             -vmin*c1(i,j,li)*rhoi_ch4 
             fluxi(ihor,nz+1-i,2)=fluxi(ihor,nz+1-i,2)
     &            -vmin*c2(i,j,li)*rhoi_c2h6
             fluxi(ihor,nz+1-i,3)=fluxi(ihor,nz+1-i,3)
     &            -vmin*c3(i,j,li)*rhoi_c2h2
           enddo
         enddo
  
c        En theorie, les diceX sont NEGATIF (en sedimentant on ne fait que perdre de la glace)
c        Les precipitations sont comptees positivement. (ET ON NE PREND QUE DES VALEURS POSITIVES !)

         precip(ihor,1)=AMAX1(-dice1/rhoi_ch4,0.)     ! m3/m2 = m :)
         precip(ihor,2)=AMAX1(-dice2/rhoi_c2h6,0.)
         precip(ihor,3)=AMAX1(-dice3/rhoi_c2h2,0.)
         precip(ihor,4)=AMAX1(-dice4/rhol,0.)

         do i=1,nz
           do j=1,nrad
             tab1(i,j)=c(i,j,li)      ! li=1
             tab2(i,j)=c1(i,j,li)     ! li=1
             tab3(i,j)=c2(i,j,li)     ! li=1
             tab4(i,j)=c3(i,j,li)     ! li=1
           enddo
         enddo 

         return  

         end 



*______________________________________________________________________

         real function lambda2(j,indic)
*
*------------------------------------------------------------------*
*   fonction calculant le libre parcours moyen des molecules        *
*   atmospheriques( rayon =ra) se trouvant dans la couche no j.     *
*   pour indic=0  ...... la particule se trouve a la frontiere entre*
*                         les couches j et j-1                      *
*   pour indic=1  ...... la particule se trouve au milieu de la     *
*                          la couche j                              *
*------------------------------------------------------------------*
*
*  declaration des blocs communs
*------------------------------
         use dimphy
         IMPLICIT NONE
#include  "dimensions.h"
#include  "microtab.h"
#include  "varmuphy.h"

         common/x2frac/rfg,dfg

*  declaration des variables communes
*  ----------------------------------

         real rfg(nz),dfg(nz,nrad)

*   declaration des variables internes
*   ----------------------------------

         integer indic,j
         real pp,ra

         ra=1.75e-10

*  traitement
*  ----------

         if (indic.eq.0) then
           pp=pb(j)
         else
           if (indic.ne.1) then
             print*,'erreur argument fonction lambda'
             stop 
           endif
           pp=p(j)
         endif

         lambda2=kbz*t(j)/(4*sqrt(2.)*pi*(ra**2)*pp)
         end

*******************************************************************************

         real function knu2(j,k,indic)
*
*--------------------------------------------------------------*
*   fonction calculant le nombre de knudsen d'une particule    *
*   d'aerosol de rayon rad(k) se trouvant dans la couche no j    *
*   indic ......  idem function lambda                         *
*--------------------------------------------------------------*
*
*  declaration des blocs communs
*------------------------------
         use dimphy
         IMPLICIT NONE
#include  "dimensions.h"
#include  "microtab.h"
#include  "varmuphy.h"

         common/cldpart/rad,mmu
         common/x2frac/rfg,dfg


*  declaration des variables communes
*  ----------------------------------

         real rad(nz,nrad),mmu(nz,nrad)
         real rfg(nz),dfg(nz,nrad)


*   declaration des variables internes
*   ----------------------------------

         integer indic,j,k 
         real  lambda2,rfk

*  traitement
*  ----------

         if (indic.ne.0 .and.indic.ne.1) then
           print*,'erreur argument fonction knu'
           stop 
         endif
    

         rfk=(rad(j,k)**(3./dfg(j,k)))*((rfg(k))**(1.-3./dfg(j,k)))
         knu2=lambda2(j,indic)/rfk
 
         end

*****************************************************************************

         real function nud2(j,indic)
*
*--------------------------------------------------------------*
*   fonction calculant la viscosite dynamique (en USI) de l air *
*   d apres la formule de Sutherlant a l altitude j             *
*   indic  ......... idem fonction lambda                       *
*--------------------------------------------------------------*
*
         use dimphy 
         IMPLICIT NONE
#include  "dimensions.h"
#include  "microtab.h"
#include  "varmuphy.h"          
         integer indic,j
         real nud0,c,tt
         real rfg(nz),dfg(nz,nrad)

         common/x2frac/rfg,dfg

*
         nud0=1.74e-5
         c=109.

         if(indic.ne.0.and.indic.ne.1) then
           print*,'erreur argument fonction nud'
           stop 
         endif

         if (indic.eq.0) tt=tb(j)
         if (indic.eq.1) tt=t(j)
         nud2=nud0*sqrt(tt/293)*(1+c/293)/(1+c/tt)
         end

****************************************************************************

         real function vitesse2(j,k,indic)
*
*-----------------------------------------------------------------*
*   fonction calculant la vitesse de chute d une particule de rayon*
*   k se trouvant a l altitude j  suivant la valeur du nombre de   *
*    Knudsen                                                       *
*   indic ....... idem function lambda                             *
*-----------------------------------------------------------------*
*

*  declaration des blocs communs
*------------------------------
         use dimphy
         IMPLICIT NONE
#include  "dimensions.h"
#include  "microtab.h"
#include  "varmuphy.h"          
         common/cldpart/rad,mmu
         common/x2frac/rfg,dfg
         common/x2ctps/li,lf,dt

*  declaration des variables communes
*  ----------------------------------

         real rad(nz,nrad),mmu(nz,nrad)
         real rfg(nz),dfg(nz,nrad)


*   declaration des variables internes
*   ----------------------------------

         integer indic,j,k,li,lf
         real w,g,m,a0,zz,knu2,nud2,knud,tt,rhoh
         real vlimite,akncx,rbis,rfk
         real dt

*   traitement
*   ----------

         if (indic.ne.0.and.indic.ne.1) then
           print*,'erreur argument fonction vitesse'
           stop 
         endif

         if(indic.eq.0) then
           zz=z(j)+dz(j)/2.
           tt=tb(j)
           rhoh=rhob(j)
         endif
         if(indic.eq.1) then
            zz=z(j)
            tt=t(j)
            rhoh=rho(j)
         endif

         g=g0*(rtit/(rtit+zz))**2
         a0=0.74
         m=(ach4(j)*mch4+aar(j)*mar+an2(j)*mn2)/nav
         knud=knu2(j,k,indic)

c        akncx=aknc
c        if(df(k).gt.2.5) akncx=2.7

c        if(knud.ge.akncx) then
c        rbis=(rad(j,k)**(3.-6./dfg(j,k)))*((rfg(k))**(-2.+6./dfg(j,k)))
c        w=a0*g*rbis*mmu(j,k)/(rhoh*sqrt(8*kbz*tt/(pi*m)))
c        endif

          rfk=(rad(j,k)**(3./dfg(j,k)))*((rfg(k))**(1.-3./dfg(j,k)))
          w=2./9.*rfk**(dfg(j,k)-1.)*rfg(k)**(3.-dfg(j,k))*g*mmu(j,k)
     &       /nud2(j,indic)

          w=w*(1+1.2517*knud+0.4*knud*exp(-1.1/knud))

          w=w!*3.   ! on tient compte de la largeur de distribution... a affiner
          vitesse2=w 

         
         end
***********************************************************************
         real function kd2(h)
*
*--------------------------------------------------------------------*
*   cette fonction calcule le coefficient du terme de  eddy diffusion *
*   a l altitude j                                                    *
*--------------------------------------------------------------------*
*
         use dimphy
         IMPLICIT NONE
#include  "dimensions.h"
#include  "microtab.h"
#include  "varmuphy.h"

         common/x2frac/rfg,dfg

         real zbx
         real rfg(nz),dfg(nz,nrad)

         integer h

         zbx=z(h)+dz(h)/2.
         if(zbx.le.42000.) then
c           kd2=1.64e+12*(pb(h)/(kbz*tb(h)))**(-1./2.)
           kd2=4.
         else
           kd2=1.64e+12*(pb(h)/(kbz*tb(h)))**(-1./2.)
           kd2=0.0*kd2
         endif
         kd2=00.


         return
         end


*____________________________________________________________________________

         subroutine init2
*
*--------------------------------------------------------------------*
*   cette routine effectue  :                                         *
*		1)  interpolation a partir des donnees initiales des  *
*                    valeurs de p,t,rho,ach4,aar,an2  sur la grille   *
*                2) initialisation des constantes (common/x2phys/)      *
*	         3) initialisation des variables temporelles (common  *
*                     /temps/)                                        *
*                4) definition des grilles en rayon et verticale      *
*                5)  initialisation de c(z,r,t) avec les donnees du   *
*                      fichier unit=1                                 *
*                                                                     *
*   les donnees sont des valeurs caracterisques de l atmosphere de    *
*     TITAN  ( voir Lelouch and co )                                  *
*--------------------------------------------------------------------*

*  declaration des blocs communs
*------------------------------
         use dimphy
         IMPLICIT NONE
#include  "dimensions.h"
#include  "microtab.h"
#include  "varmuphy.h"          
         common/x2ctps/li,lf,dt
         common/cldpart/rad,mmu

*  declaration des variables communes
*  ----------------------------------

         integer li,lf
         real dt
         real rad(nz,nrad), mmu(nz,nrad)


*  declaration des variables internes
*  ----------------------------------
         integer nzd,i,ii
         parameter (nzd=254)
         integer limsup,liminf,j1,j2
         real zd(nzd),ach4d(nzd),rap
         real m


*  initialisation des variables temporelles
*  ----------------------------------------

         li=1
         lf=2


*  interpolation de xch4,xar et xn2 sur la grille
*  ----------------------------------------------

*  donnees initiales (Lellouch et al,87) 
*  ------------------------------------- 

c        print*,'****** init'
         do 1 i=1,168
           zd(i)=(1000.-5*(i-1))*1000.
1        continue
         do 2 i=1,78
           zd(168+i)=(160.-2*(i-1))*1000.
2        continue
         do 3 i=1,4
           zd(246+i)=(5.-(i-1))*1000.
3        continue
         do 4 i=1,4
           zd(250+i)=(1.5-(i-1)*0.5)*1000.
4        continue

         data (ach4d(i),i=1,168)/168*1.5e-2/
         data (ach4d(i),i=169,254)/63*1.5e-2,1.6e-2,1.8e-2,1.8e-2,
     &   1.9e-2,2.e-2,2.1e-2,2.3e-2,2.5e-2,2.8e-2,3.1e-2,3.6e-2,
     &   4.1e-2,4.7e-2,5.7e-2,6.7e-2,7.5e-2,7*8.e-2/

         liminf=0
         limsup=0

*  interpolation des taux de melange de ch4,ar,n2  
*-----------------------------------------------   

         do 20 j1=1,nz
           do 21 j2=1,nzd
             if( zd(j2).le.z(j1)) goto 22
21           continue
22           liminf=j2
             if (zd(liminf).eq.z(j1) )then
               ach4(j1)=ach4d(liminf)
               goto 20
             endif
             if (j2.ne.1) then
               limsup=j2-1
             else
               limsup=j2
             endif
             if (limsup.eq.liminf) then
               ach4(j1)=ach4(limsup)
             else
               ach4(j1)=ach4d(liminf)-(ach4d(limsup)-ach4d(liminf))/
     s         (zd(limsup)-zd(liminf))*(zd(liminf)-z(j1))
             endif
20       continue

*   rap= aar/an2  cst sur l altitude

         rap=0.191
         do 23 i=1,nz
           an2(i)=(1.-ach4(i))/(1.+rap)
           aar(i)=rap*an2(i)
23       continue

         do 24 i=1,nz
           m=ach4(i)*mch4+an2(i)*mn2+aar(i)*mar
           rho(i)=p(i)*m/(rgp*t(i))
24       continue

         do 34 i=1,nz
           m=ach4(i)*mch4+an2(i)*mn2+aar(i)*mar
           rhob(i)=pb(i)*m/(rgp*tb(i))
c          print*,pb(i),m,rgp,tb(i),rhob(i),rho(i)
34       continue

*  fin d interpolation des taux de melange
*----------------------------------------  

c        print*,'**** fin init'
540      continue
         return

500       print*,'erreur lecture initialisation de c...erreur=',ii
          stop

         end

*__________________________________________________________________________

         subroutine sedifn
*
*------------------------------------------------------------------*
*   cette routine calcule l evolution de la fonction de distribution*
*   c(z,r,t) pour les phenomenes de sedimentation et de diffusion   *
*------------------------------------------------------------------*
*
*
*  declaration des blocs communs
*------------------------------
         use dimphy
         IMPLICIT NONE
#include  "dimensions.h"
#include  "microtab.h"
#include  "varmuphy.h"          

         common/x2ctps/li,lf,dt
         common/x2con/ctmp,c1,c2,c3,caer
         common/cldpart/rad,mmu
         common/x2frac/rfg,dfg


*  declaration des variables communes
*  ----------------------------------

         integer li,lf
         real  dt
         real  c(nz,nrad,2)
         real  c1(nz,nrad,2)
         real  c2(nz,nrad,2)
         real  c3(nz,nrad,2)
         real caer(nz,nrad,2)
         real  ctmp(nz,nrad,2)
         real  rad(nz,nrad),mmu(nz,nrad)
         real rfg(nz),dfg(nz,nrad)


*  declaration des variables internes
*  ----------------------------------

         real w,w1,dzbX,dc
         integer i,j,k,ii,nb
         double precision sigma,theta,hc,l,rap,cmp,wp 
         double precision fs(nz+1),ft(nz+1)
         real as(nz),bs(nz),cs(nz),ds(nz)
         double precision asi(nztop:nz),bsi(nztop:nz),
     &                    csi(nztop:nz)
         double precision dsi(nztop:nz),xsol(nztop:nz)
         real vitesse2,kd2

         external dtridgl

*  resolution
*------------
  
         do 5 ii=1,4
           do  k=1,nrad  
             do  j=nztop,nz
               if( ii.eq.1 ) c(j,k,li)=ctmp(j,k,li)
               if( ii.eq.2 ) c(j,k,li)=c1(j,k,li)
               if( ii.eq.3 ) c(j,k,li)=c2(j,k,li)
               if( ii.eq.4 ) c(j,k,li)=c3(j,k,li)
             enddo
           enddo

           do 10 k=1,nrad  
             do 20 j=nztop,nz
               if (j.eq.1) goto 20
*  calcul de la vitesse corrigee

               dzbX=(dz(j)+dz(j-1))/2.
               w= -1*vitesse2(j,k,0)

               if (kd2(j).ne.0.) then
                 theta=0.5*(w*dzbX/kd2(j)+log(rho(j-1)/rho(j)))
                 if (theta.ne.0) then
                   sigma=1./dtanh(theta)-1./theta
                 else
                   sigma=1.
                 endif
               else
                 sigma=1.
               endif

               if(c(j,k,li).eq.0.) then
                 rap=10.
               else
                 rap=c(j-1,k,li)/c(j,k,li)
                 if( rap.gt.10.) rap=10.
                 if( rap.lt.0.1) rap=0.1
               endif

               if (rap.gt.0.9 .and. rap.lt.1.1) then
                 w1=w
               else
                 if(w.ne.0) then
                   hc=dzbX/dlog(rap)
                   l=dzbX/(w*dt)*(dexp(-w*dt/hc)-1.)/(1.-rap)
                   wp=w*1.d0
                   cmp=dlog(-wp)+abs(sigma)*dlog(l)
                   if (cmp.gt.38) then
                     goto 20
                   endif
                   w1=-dexp(cmp)
                 else
                   w1=0.
                 endif
               endif

*   calcul des flux aux interfaces


               if (kd2(j).ne.0.) then
                 if (theta.ne.0.) then
                   ft(j)=(w1+log(rho(j-1)/rho(j))*kd2(j)/dzbX)/
     &             (dexp(2.*theta)-1.)
                   fs(j)=ft(j)*dexp(2.*theta)
                 else
                   ft(j)=kd2(j)/dzbX
                   fs(j)=kd2(j)/dzbX
                 endif
               else
                 if (w1.lt.0.)then
                   ft(j)=-w1
                   fs(j)=0.
                 else
                   ft(j)=0.
                   fs(j)=w1
                 endif
               endif

20           continue

*  conditions aux limites pour les flux aux interfaces

             fs(1)=0.
             ft(1)=0.
             fs(nz+1)=0.
             ft(nz+1)=-w1

*  calcul des coefficients de l equation discrete

             do 21 j=nztop,nz
               as(j)=-dz(j)/dt
               bs(j)=-ft(j)
               cs(j)=ft(j+1)+fs(j)-dz(j)/dt
               ds(j)=-fs(j+1)
               if ( cs(j).gt.0) goto 100
21           continue

*  cas explicite (mu=0) : calcul de la fonction c(z,r,t+1)

             do 22 j=nztop,nz-1

               if (j.eq.nztop) then
                 dc=(cs(nztop)*c(nztop,k,li)+ds(nztop)
     &                        *c(nztop+1,k,li))/as(nztop)
                 c(nztop,k,lf)=dc
                 goto 22
               endif

               dc=(bs(j)*c(j-1,k,li)+cs(j)*c(j,k,li)+ds(j)*c(j+1,k,li))
     s            /as(j)
               c(j,k,lf)=dc

22           continue

             dc=(bs(nz)*c(nz-1,k,li)+cs(nz)*c(nz,k,li))/as(nz)
             c(nz,k,lf)=dc

             if (nztop.ne.1) then
               do 32 j=1,nztop-1
                 c(j,k,lf)=c(j,k,li)
32             continue
             endif

             goto 10

100          continue

*  cas implicite (mu=1) : calcul de la fonction c(z,r,t+1)

             do 101 j=nztop,nz
               asi(j)=ft(j)
               bsi(j)=-(ft(j+1)+fs(j)+dz(j)/dt)
               csi(j)=fs(j+1)
               dsi(j)=-dz(j)/dt*c(j,k,li)
101          continue

*  inversion de la matrice tridiagonale 

             nb=nz-nztop+1

             call dtridgl(nb,asi,bsi,csi,dsi,xsol) 

             do 102 j=nztop,nz
               c(j,k,lf)=xsol(j)
102          continue

             if (nztop.ne.1) then
               do 110 j=1,nztop-1
                 c(j,k,lf)=c(j,k,li)
110            continue
             endif



10         continue

           do  k=1,nrad  
             do  j=nztop,nz
               if( ii.eq.1 ) ctmp(j,k,lf)=c(j,k,lf)
               if( ii.eq.2 ) c1(j,k,lf)  =c(j,k,lf)
               if( ii.eq.3 ) c2(j,k,lf)  =c(j,k,lf)
               if( ii.eq.4 ) c3(j,k,lf)  =c(j,k,lf)
             enddo
           enddo

5        continue

         return

         end

*__________________________________________________________________________

         subroutine sedifn_fast(ihor,dice1,dice2,dice3,dice4)
*
*------------------------------------------------------------------    *
*   cette routine calcule l evolution de la fonction de distribution   *
*   c(z,r,t) pour les phenomenes de sedimentation  {pas de diffusion}  *
*   dice1 = delta glace CH4                                            *
*   dice2 = delta glace C2H6                                           *
*   dice3 = delta glace C2H2                                           *
*   dice4 = delta Volume noyaux                                        *
*------------------------------------------------------------------    *
*
*  declaration des blocs communs
*------------------------------
         use dimphy
         IMPLICIT NONE
#include  "dimensions.h"
#include  "microtab.h"
#include  "varmuphy.h"          

         common/x2ctps/li,lf,dt
         common/x2con/c,c1,c2,c3,caer
         common/cldpart/rad,mmu
         common/x2frac/rfg,dfg


*  declaration des variables communes
*  ----------------------------------

         integer li,lf,i,j,k,ihor
         integer jinf,jsup,jj,iiter
         real dt
         real c(nz,nrad,2)
         real c1(nz,nrad,2)
         real c2(nz,nrad,2)
         real c3(nz,nrad,2)
         real caer(nz,nrad,2)
         real ci(nz,nrad,2)
         real ci1(nz,nrad,2)
         real ci2(nz,nrad,2)
         real ci3(nz,nrad,2)
         real rad(nz,nrad),mmu(nz,nrad)
         real rfg(nz),dfg(nz,nrad)
         real puit(nz)
c ------ echange est cree sur la taille maxi mais n'est utilisee
c        que sur la dim geree par le proc (klon ou jjm+1)
         integer ngrid
         parameter (ngrid=(jjm-1)*iim+2)  ! = taille maximum
         real echange(nz,nz,ngrid)
c pas genial mais vu que c est tres local, pas de soucis a priori en parallele.
         real bilan1,bilan2,bilan3,bilan4,bilan5
         real bilan11,bilan12,bilan13,bilan14,bilan15
         real compte,compte2,xepl
         real dice1,dice2,dice3,dice4
        


*  declaration des variables internes
*  ----------------------------------

         real vitesse2,kd2
         real w,ws,wi,zs,zi,alpha,v0,deltazs,deltazi
         real zni,znip1,xcnt,ichx,arg1,arg2,xft,xf
         real argexp,explim

         save echange


*  resolution
*------------
  

         bilan1=0.
         bilan2=0.
         bilan3=0.
         bilan4=0.
         bilan5=0.
         do  k=1,nrad  
           do  j=nztop,nz
             ci(j,k,li)=  c(j,k,li)*dzb(j)   ! li
             ci1(j,k,li)=c1(j,k,li)*dzb(j)
             ci2(j,k,li)=c2(j,k,li)*dzb(j)
             ci3(j,k,li)=c3(j,k,li)*dzb(j)
             bilan5=bilan5+ci(j,k,li)
             bilan1=bilan1+ci1(j,k,li)
             bilan2=bilan2+ci2(j,k,li)
             bilan3=bilan3+ci3(j,k,li)
             bilan4=bilan4+
     &       ci(j,k,li)*4./3.*pi*rf(k)**3.*vrat_e**(k-imono)

             ci(j,k,lf)= 0.                 ! lf
             ci1(j,k,lf)=0.
             ci2(j,k,lf)=0.
             ci3(j,k,lf)=0.
           enddo
         enddo

*  calcul de la matrice d echange 
*----------------------------------------------------------------

         do j=nztop,nz
           do i=nztop,nz
             echange(i,j,ihor)=0.
           enddo
         enddo

         do 20 i=nztop,nz
           puit(i)=0.
           do 30 k=1,1 
             xcnt=0.
             ICHX=1         ! extrapolation 0: lineaire 1: exponentielle
             IF(ICHX.eq.0  .or. ICHX.eq.2) THEN 
               ws=vitesse2(i,k,0)
               if(i.lt.nz) wi=vitesse2(i+1,k,0)
               if(i.eq.nz) wi=vitesse2(i  ,k,1)
               w=(ws+wi)/2.
               zni  =zb(i)-w*dt
               znip1=zb(i)-dzb(i)-w*dt
             ENDIF

             explim=30.

             IF(ICHX.eq.1 .or.  ICHX.eq.2) THEN 
               ws=vitesse2(i,k,0)
               zs=zb(i)
               wi=vitesse2(i,k,1)
               zi=z(i)

c              if(wi.eq.ws)  wi=ws/1.001   ! sinon ca plante !
               if(abs(wi-ws)/wi .le. 1.e-3)  wi=ws/1.001   ! sinon ca plante !

               if(wi.ne.0.) alpha= alog(ws/wi) /(zs-zi)
               argexp=alpha*zs
               if(argexp.lt.-explim) argexp=-explim
               if(argexp.gt. explim) argexp=+explim
               v0   = ws/exp(argexp)

               argexp=alpha*zb(i)
               if(argexp.lt.-explim) argexp=-explim
               if(argexp.gt. explim) argexp=+explim
               arg1=1.+v0*alpha*exp(argexp)*dt

               argexp=alpha*(zb(i)-dzb(i))
               if(argexp.lt.-explim) argexp=-explim
               if(argexp.gt. explim) argexp=+explim
               arg2=1.+v0*alpha*exp(argexp)*dt

               iiter=0
101            continue 

               if (iiter.le.25) then
                 if(arg1.le.0.  .or. arg2.le.0.) then 
                   print*,ihor,i,iiter, 'ajustement vitesse',arg1,arg2
                   print*,ws,wi,     ' w1 w2 anc valeurs'
                   print*,alpha,     ' alpha anc valeurs'
                   print*,rad(i,k),  'r(j,k)'
                   print*,mmu(i,k),  ' mmu(j,k)'
                   print*,t(i),      ' t(j,k)'
                   print*,arg1,arg2, ' arg1 arg2 anc valeurs'
                   arg2=1+(arg2-1.)/2.
                   arg1=1+(arg1-1.)/2.
                   iiter=iiter+1
                   print*,arg1,arg2, ' arg1 arg2 nle valeurs'
                   goto 101
                 endif
               else 
                 stop 
               endif

               if(arg1.lt.0.) print*,'arg1:',arg1
               if(arg2.lt.0.) print*,'arg2:',arg2

               deltazs=-alog(arg1)/alpha
               deltazi=-alog(arg2)/alpha

               zni  =zb(i)+deltazs
               znip1=zb(i)-dzb(i)+deltazi

             ENDIF

             if(zni.ne.znip1) xft=zni/(zni-znip1)
             if(zni.eq.znip1 .and. zni.le.0.) xft=0.
             if(zni.eq.znip1 .and. zni.gt.0.) then
               xft=0.
               print*,'zni..znip1', zni,znip1
             endif

*  Si des aerosols touchent le sol (zni < 0 )  alors on fixe
*  le niveau a 0, et on elimine les aerosols correspondants
*-----------------------------------------------------------

             if(znip1 .lt. 0.) znip1=0.
             if(zni   .lt. 0.) zni  =0.
            
             if(zni.le.0.  .and.  znip1 .le. 0.) then 
c            print*,'voie 1 / disparition complete'
               xft=0.
               xf=1.
               xcnt=xcnt+xf
               puit(i)=puit(i)+xf
             endif

             if(zni.gt.0.  .and.  znip1 .le. 0.) then 
c            print*,'voie 2 / disparitipon partielle'
               xf=(1-xft)
               xcnt=xcnt+xf
               puit(i)=puit(i)+xf
             endif

             if(zni.gt.0.  .and.  znip1 .gt. 0.) then 
c            print*,'voie 3 / pas de disparition'
               xft=1.
               xf=(1-xft)
               xcnt=xcnt+xf
               puit(i)=puit(i)+xf
             endif

             jsup=nz+1
             jinf=nz+1
             do j=nztop,nz
               if(zni.le.zb(j)  .and. zni.ge.zb(j)-dzb(j))   jsup=j
               if(znip1.le.zb(j).and. znip1.ge.zb(j)-dzb(j)) jinf=j
             enddo
             if(zni   .ge. 0. .and. zni   .lt. 1.e-3)   jsup=nz
             if(znip1 .ge. 0. .and. znip1 .lt. 1.e-3)   jinf=nz
           

*  Volume inclu dans un seul niveau
*----------------------------------
  
             if (jsup .eq. jinf. and. jsup.ge. nz+1) then 
               xcnt=xcnt+1.
               print*,'cas impossible'
               print*,'alpha= ',alpha
               print*,'ws wi ',ws,wi
               print*,'deltazs deltazi ',deltazs,deltazi 
               print*,' r(i,k) mmu(i,k)', rad(i,k),mmu(i,k)
               print*,' t(j,k)',t(i)
               print*,zni,znip1,jsup,jinf
               print*,'STOP'
               STOP
             endif

             if (jsup .eq. jinf. and. jsup.le. nz) then 
               xf=1.
               xcnt=xcnt+xft*xf
               if(jinf.le.nz) then 
                 echange(jinf,i,ihor)=echange(jinf,i,ihor)+xft*xf
               endif
             endif

*  Volume a cheval sur 2  niveaux
*--------------------------------

             if (jinf .eq. jsup+1) then
               xf=(zni-zb(jsup)+dzb(jsup))/(zni-znip1)
               xcnt=xcnt+xf*xft
               if(jsup.le.nz) then 
                 echange(jsup,i,ihor)=echange(jsup,i,ihor)+xft*xf
               endif
               xf=(zb(jinf)-znip1)/(zni-znip1)
               xcnt=xcnt+xf*xft
               if(jinf.le.nz) then 
                 echange(jinf,i,ihor)=echange(jinf,i,ihor)+xft*xf
               endif
             endif

*  Volume a cheval sur 3 ou plus de niveaux
*------------------------------------------

             if (jinf .gt. jsup+1) then

c             print*,' voie C / dans N  cases'
               xf=(zni-zb(jsup)+dzb(jsup))/(zni-znip1)
               xcnt=xcnt+xf*xft
               if(jsup.le.nz) then 
                 echange(jsup,i,ihor)=echange(jsup,i,ihor)+xft*xf
               endif

               xf=(zb(jinf)-znip1)/(zni-znip1)
               xcnt=xcnt+xf*xft
               if(jinf.le.nz) then 
                 echange(jinf,i,ihor)=echange(jinf,i,ihor)+xft*xf
               endif

               do jj=jsup+1,jinf-1
                 xf=(zb(jj)-zb(jj+1))/(zni-znip1)
                 xcnt=xcnt+xf*xft
                 if(jj.le.nz) then 
                   echange(jj,i,ihor)=echange(jj,i,ihor)+xft*xf
                 endif
               enddo

             endif 


*  et sur les rayons...
*---------------------
30         continue

*  fin de la grande boucle sur les altitudes...
*----------------------------------------------
20       continue

*  Calcul etat final     Cfinal = Echange*initial 
*----------------------------------------------

         compte=0.
         compte2=0.
         do  j=1,nz
           xepl=0.
           do  jj=1,nz
             xepl=xepl+echange(jj,j,ihor)
             compte=compte+echange(jj,j,ihor)
             compte2=compte2+echange(jj,j,ihor)
           enddo
           compte2=compte2+puit(j)
         enddo
  
         if(abs(compte2-nz) .gt. 1.e-4) 
     &   print*,'Matrice calculee#',ihor,'tx d expl (/55):',
     &   compte,compte2 


*  Fin du calcul de la matrice d*echange
*----------------------------------------------

         do  j=nztop,nz
           do k=1,nrad

             do  jj=nztop,nz
               ci(j,k,lf)=ci(j,k,lf)+
     &         echange(j,jj,ihor)*ci(jj,k,li)
               ci1(j,k,lf)=ci1(j,k,lf)+
     &         echange(j,jj,ihor)*ci1(jj,k,li)
               ci2(j,k,lf)=ci2(j,k,lf)+
     &         echange(j,jj,ihor)*ci2(jj,k,li)
               ci3(j,k,lf)=ci3(j,k,lf)+
     &         echange(j,jj,ihor)*ci3(jj,k,li)
             enddo
           enddo
         enddo


*  Controles et affichage Bilan
*----------------------------------------------
  
         bilan11=0.
         bilan12=0.
         bilan13=0.
         bilan14=0.
         bilan15=0.
         do  k=1,nrad  
           do  j=nztop,nz
             c(j,k,lf) =ci(j,k,lf)/dzb(j)
             c1(j,k,lf)=ci1(j,k,lf)/dzb(j)
             c2(j,k,lf)=ci2(j,k,lf)/dzb(j)
             c3(j,k,lf)=ci3(j,k,lf)/dzb(j)
             bilan15=bilan15+ ci(j,k,lf)
             bilan11=bilan11+ci1(j,k,lf)
             bilan12=bilan12+ci2(j,k,lf)
             bilan13=bilan13+ci3(j,k,lf)
             bilan14=bilan14+
     &       ci(j,k,lf)*4./3.*pi*rf(k)**3.*vrat_e**(k-imono)
           enddo
         enddo

c           print*,'sedifn_fast Bilans:'        
c     &     ,bilan11,bilan12,bilan13        
c           print*,'Bilan1:',bilan1,bilan11
c           print*,'Bilan2:',bilan2,bilan12
c           print*,'Bilan3:',bilan3,bilan13
c           print*,'Bilan5:',bilan5,bilan15
  
         dice1=0.
         dice2=0.
         dice3=0.
         dice4=0.
         dice1=(bilan11-bilan1)*rhoi_ch4     !glace 1    m^3.m^-2 * kg.m^-3   pourquoi ????
         dice2=(bilan12-bilan2)*rhoi_c2h6    !glace 2
         dice3=(bilan13-bilan3)*rhoi_c2h2    !glace 3
         dice4=(bilan14-bilan4)*rhol         !noyaux

         return
         end


         subroutine coagul2(ihor)

*********************************************************
*  ce programme calcule la nouvelle concentration dans   *
*  le a ieme intervalle de rayon, a l'altitude h, a      *
*  l'instant t+dt                                        *
*********************************************************

         use dimphy
         IMPLICIT NONE
#include  "dimensions.h"
#include  "microtab.h"
#include  "varmuphy.h"          

*  declaration des blocs communs
*------------------------------

         common/x2ctps/li,lf,dt
         common/x2con/c,c1,c2,c3,caer
         common/cldpart/rad,mmu

*  declaration des variables
*  --------------------------

         integer li,lf
         real dt
         real c(nz,nrad,2), c1(nz,nrad,2), c2(nz,nrad,2),
     &        c3(nz,nrad,2)
         real caer(nz,nrad,2)
         real rad(nz,nrad),mmu(nz,nrad)

*  declaration des variables propres au ss-programme
*  -------------------------------------------------

         integer h,a,ihor,i
         real pr,pe,eta1,eta2
         real  sum11,sum12,sum13,sum21,sum22,sum23
         real rfb,rfx,rpr

*   traitement
*   ----------


         sum11=0.
         sum12=0.
         sum13=0.
         sum21=0.
         sum22=0.
         sum23=0.
         
         do 10 h=nztop,nz
           do 11 a=1,nrad  ! boucle aerosol secs
             call pertpro2(ihor,h,a,pe,pr)

             if(1+dt*pe .gt. 1.1) print*,a,1+dt*pe,' scav'
             caer(h,a,lf)=(caer(h,a,li)+pr*dt)/(1+dt*pe)
             c(h,a,lf)=c(h,a,li)
             c1(h,a,lf)=c1(h,a,li)
             c2(h,a,lf)=c2(h,a,li)
c
c           eta1=0.
c           eta2=0.
c           if(c(h,a,li) .ne. 0.)  eta1=c1(h,a,li)/c(h,a,li)
c           if(c(h,a,li) .ne. 0.)  eta2=c2(h,a,li)/c(h,a,li)
c            c(h,a,lf) =( c(h,a,li)+pr*dt)/(1+dt*pe)
c            c1(h,a,lf)=(c1(h,a,li)+eta1*pr*dt)/(1+dt*pe)
c            c2(h,a,lf)=(c2(h,a,li)+eta2*pr*dt)/(1+dt*pe)
c
c
c           sum11=sum11+c(h,a,li)
c           sum12=sum12+c1(h,a,li)
c           sum13=sum13+c2(h,a,li)
c
c           sum21=sum21+ c(h,a,lf)
c           sum22=sum22+c1(h,a,lf)
c           sum23=sum23+c2(h,a,lf)
c
c
11         continue
10       continue


         if (nztop.ne.1) then
           do 12 h=1,nztop-1
             do 12 a=1,nrad
               c(h,a,lf)=c(h,a,li)
               c1(h,a,lf)=c1(h,a,li)
               c2(h,a,lf)=c2(h,a,li)
               caer(h,a,lf)=caer(h,a,li)
12         continue
         endif

         return
         end
       
       
*__________________________________________________________________________

         subroutine  calcoag2

***************************************************************
*                                                              *
*   Ce programme calcule les coefficients de collection  d'une *
*  particule de rayon x avec une particule de rayon b a une    *
*  altitude donnee h                                           *
***************************************************************  

*  declaration des blocs communs
*------------------------------
         use dimphy
         IMPLICIT NONE
#include  "dimensions.h"
#include  "microtab.h"
#include  "varmuphy.h"          
         
         common/x2ctps/li,lf,dt
         common/x2con/c,c1,c2,c3,caer
         common/cldpart/rad,mmu
         common/x2coag/k
         common/x2frac/rfg,dfg

*  declaration des variables
*  --------------------------

         integer li,lf
         real dt
         real knu2,nud2,k(nz,nrad,nrad)
         real c(nz,nrad,2), c1(nz,nrad,2), c2(nz,nrad,2),
     &        c3(nz,nrad,2)
         real caer(nz,nrad,2)
         real rad(nz,nrad),mmu(nz,nrad)
         real rfg(nz),dfg(nz,nrad)


*  declaration des variables propres au ss-programme
*  -------------------------------------------------

         integer h,b,x,ihor,i
         real nua,lambb,lambx,knb,knx,alphab,alphax,d,e,f,kcg
         real db,dx,rm,dm,deltab,deltax,del,g,beta,gx,gb
         real rfb,rfx,rpr
         real*8 ne,qe,epso
         real*8 corelec,yy

         real kco,vx,vb,vitesse2,sto,ee,a,dd,bb,p0,t0,l0,ccol
         real st(37),ef(37)
         real vitesse,knu
         external vitesse,knu 


*  initialisation
*  --------------



*   -nombres de STOCKES

         data(st(i),i=1,37)/1.35,1.5,1.65,1.85,2.05,2.25,2.5,2.8,3.1,
     s    3.35,3.6,3.95,4.3,4.7,5.05,5.45,5.9,6.4,7.,7.6,8.3,9.05,9.9,
     s       10.9,11.1,13.5,15.3,17.25,20.5,24.5,30.4,39.3,48,57,86.,
     s       187.,600./

*   -coef. d'efficacite de collection

         ef(1)=3.75
         ef(2)=8.75
         do 11 i=3,37
           ef(i)=ef(i-1)+2.5
11       continue

         do 2 i=1,37
           ef(i)=ef(i)*1e-2
2        continue

         qe=1.6e-19
         ne=-30e+6
         epso=1e-9/(36*pi)

         d=1.257
         e=0.4
         f=-1.1

*   iteration sur z
  
         do 1 h=1,nz

           nua=nud2(h,1)      

*   iteration sur les rayons

           do 1 b=1,nrad  ! boucle aerosols secs : indice '' 

             knb=knu(h,b,1)      ! knu et vitesse se trouvent dans brume.F
             vb=vitesse(h,b,1)   ! ils concernent les aerosols 

             do 1 x=1,nrad  !boucles gouttes : indice '2'

               knx=knu2(h,x,1)
               vx=vitesse2(h,x,1)

**   COAGULATION  ****************************************************
**  --------------****************************************************
*  calcul du terme correcteur 'slip-flow'

               alphab=d+e*exp(f/knb)
               alphax=d+e*exp(f/knx)

*  calcul du coefficient de diffusion


               rfb=(r_e(b)**(3./df(b))) *((rf(b)) **(1.-3./df(b)))
               rfx=(rad(h,x)**
     &         (3./dfg(h,x)))*((rfg(x))**(1.-3./dfg(h,x)))

               db=kbz*t(h)*(1+alphab*knb)/(6*pi*nua*rfb)
               dx=kbz*t(h)*(1+alphax*knx)/(6*pi*nua*rfx)

*  calcul du coefficient de coagulation

               rpr=rfb+rfx
               kcg=4*pi*rpr*(db+dx)

*  calcul de la vitesse thermique

               gb=sqrt(6*kbz*t(h)/(rhol*pi**2*r_e(b)**3))
               gx=sqrt(6*kbz*t(h)/(rhol*pi**2*rad(h,x)**3))

*  calcul du libre parcours apparent des aerosols

               lambb=8*db/(pi*gb)
               lambx=8*dx/(pi*gx)

*calcul  du terme correcteur beta

               rm=rpr/2.
               dm=(dx+db)/2.
               g=sqrt(gx**2+gb**2)
               deltab=(((2*rfb+lambb)**3-(4*rfb**2+lambb**2)**1.5)
     s         /(6*rfb*lambb)-2*rfb)*sqrt(2.)
               deltax=(((2*rfx+lambx)**3-(4*rfx**2+lambx**2)**1.5)
     s         /(6*rfx*lambx)-2*rfx)*sqrt(2.)
               del=sqrt(deltab**2+deltax**2)
               beta=1/((rm/(rm+del/2))+(4*dm/(g*rm)))

*  calcul du coefficient de coagulation corrige

               kcg=kcg*beta

**   COALESCENCE  **************************************************
**   -------------**************************************************

               kco=0.

               if ( b.eq. x) continue ! goto 9


*  calcul du nombre de Stockes de la petite particule

               sto=2*rhol*rfx**2*abs(vx-vb)/(9*nua*rfb)

*  calcul du coef. de Cunningham-Millikan

               a=1.246
               bb=0.42
               dd=0.87
               l0=0.653e-7
               p0=101325.
               t0=288.

               ee=1+
     &         (l0*t(h)*p0*(a+bb*exp(-dd*rfx*t0*p(h)/(l0*t(h)*p0))))
     s         /(rfx*t0*p(h))

*  calcul du nombre de Stockes corrige

               sto=sto*ee

               if (sto .le. 1.2) goto 9

               if (sto .ge. 600.) then 
                 ccol=1.
                 goto 8
               endif

*   recherche du coefficient de collection

               do 3 i=1,37
                 if (sto .gt. st(i)) then
                   goto 3
                 endif
                 if (sto .eq. st(i)) then
                   ccol=ef(i+1)
                 else
                   ccol=ef(i)
                 endif
                 goto 8
3              continue

*   calcul du coefficient de coalescence

8              kco=pi*(rfb+rfx)**2*ccol*abs(vb-vx)

9              continue

**   CORRECTION ELECTRICITE *******************************
**   ------------------------******************************

c        yy=1.d0*ne**2*r(x)*r(b)*qe**2
c     &  /(1.d0*kbz*t(h)*(r(b)+r(x))*4*pi*epso)
c        corelec=0.
c        if (yy.lt.50.) corelec=yy/(exp(yy)-1.)

               corelec=1.

c        b=aerosol
c        x=gouttes

               k(h,b,x)=(kcg+kco)*corelec

c
c        ATTENTION, IL N'Y A PLUS DE SYMETRIE...
c        k(ihor,h,x,b)=k(ihor,h,b,x)
c
c


1        continue
         return
         end

*______________________________________________________________________

         subroutine pertpro2(ihor,h,a,l_,pr_)

*****************************************************************************
*                                                                           *
*  ce programme permet le calcul du terme de production (pr) et de perte (l)*
*  pour le phenomene de coagulation                                         *
*  dans le a ieme intervalle de rayon a une altitude h                      *
****************************************************************************

*  declaration des blocs communs
*------------------------------
         use dimphy          
         IMPLICIT NONE
#include  "dimensions.h"
#include  "microtab.h"
#include  "varmuphy.h"

         common/x2ctps/li,lf,dt
         common/x2con/c,c1,c2,c3,caer 
         common/cldpart/rad,mmu
         common/x2coag/k


*  declaration des variables
*  --------------------------

         integer li,lf
         real dt
         real dr(nrad),dv(nrad)
         real k(nz,nrad,nrad)
         real c(nz,nrad,2), c1(nz,nrad,2), c2(nz,nrad,2),
     &        c3(nz,nrad,2)
         real caer(nz,nrad,2)
         real rad(nz,nrad),mmu(nz,nrad)


*  declaration des variables propres au ss-programme
*  -------------------------------------------------

         integer h,b,a,x,ihor,i
         real*8 pr,ss,s,l
         real pr_,l_,vol,del

*  traitement
*  -----------



*   production
*+++++++++++++
         s=0.d0
         ss=0.d0
         pr=0.d0

c   Pas de production d'aerosols par scavenging !!!

         pr=0.d0

*   perte
*-  - - - -

         l=0.d0


         do 10 x=1,nrad   ! boucle sur les gouttes
           l=l+k(h,a,x)*c(h,x,li)
           if(l.ne.0.d0) then 
             print*,a,x,k(h,a,x),c(h,x,li),l,' : detail coal'
           endif
  10     continue 


#ifdef  CRAY
         l_=l
         pr_=pr
#else
         l_=sngl(l)
         pr_=sngl(pr)
#endif

         return

         end

