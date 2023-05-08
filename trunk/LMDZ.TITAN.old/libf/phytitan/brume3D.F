         subroutine brume(ngrid,tab1,x1,
     &                    xnz,xnrad,taused,ihor,
     &                    pmu0,pfract,
     &                    precip)


*1  c nombre de particules de la grille de rayon r a l'altitude z 
*2  dt, pas de temps, en heure.  

*--------------------------------------------------------------*
*                                                              *
*            ENTRE 0 ET 1000 KILOMETRES                        *     
*                                                              * 
*    la dimension fractale est en tableau, attention au        *
*    raccordement entre le regime moleculaire et le regime     *
*    fluide                                                    *
*                                                              *
*    Modele microphysique:    Cabane et al.,1992 /             *
*    Modele version fractale: Cabane et al.,1993 /             *
*                                                              *
*--------------------------------------------------------------*
* VERSION DU 2 JUIN 1993  --- AUT 1994 --- 11/04/96
*
* changer: altitude de production z0=/taux de production ctot= 
*        : la charge/micron, ne
*        : df(h),rf... 
* raccordement aknc 
*
* declaration des blocs communs
*------------------------------

         use dimphy
         IMPLICIT NONE
#include "dimensions.h"
#include "microtab.h"
#include "varmuphy.h"


         real pmu0,pfract

         common/ctps/li,lf,dt
         common/con/c
         common/coag/k
         common/effets/xsaison

 
* declaration des variables communes
* ----------------------------------

         integer xnz,xnrad,ngrid
         integer li,lf,ihor
         real dt,g
         real c0(nz,nrad),c(nz,nrad,2)
         real k(nz,nrad,nrad),knu
         real xsaison
         real taused(nz,nrad)
         real precip(ngrid,5)

        

*  variables internes
*  ------------------

         integer h,ti,itime,i,j
         real tab1(nz,nrad)
         real x1 
         real somme,v1,dice3
         
         real vitesse

         save itime
         data itime/0/

* controles
* ---------

         if (nrad.ne.xnrad)  stop 'nrad.ne.xnrad'
         if (nz.ne.xnz)    stop 'nz.ne.xnz'

         do i=1,nz
           do j=1,nrad 
             c(i,j,1)=tab1(i,j)
             c(i,j,2)=0.0 
           enddo 
         enddo 

         dt=x1 

* initialisation unique
* --------------

c         if (itime.eq.0) then 
c           ITIME=1
c         endif

* initialisation
* --------------

           call init
           call calcoag

* effet saisonnier
* ----------------


         xsaison=0.
         xsaison=pmu0*4.*pfract
                              !=Pi si fract=1/2 (equinoxe) et 
                              !    si mu0(ihor)=1 sous le soleil
                              !    exactement.

c        xsaison=0.
c        if (ihor.le.9.or.ihor.ge.41) xsaison=8.  ! rapport des surfaces
c        xsaison=1.

         do i=1,nz,1
           do j=1,nrad
             v1=vitesse(i,j,0)
             g=g0*(rtit/(rtit+z(i)))**2
             taused(i,j)=rgp*t(i)/(mn2*g)/v1
           enddo
         enddo

         call coagul

         call production(ihor)
 
         li=3-li
         lf=3-lf

         call sedif(dice3)
c        En theorie, dice3 est NEGATIF (en sedimentant on ne fait que perdre des aerosols)
c        Les precipitations sont comptees positivement. (ET ON NE PREND QUE DES VALEURS POSITIVES)
         precip(ihor,5)=AMAX1(-dice3/rhol,0.)    ! m3/m2=m 

         li=3-li
         lf=3-lf

         do i=1,nz
           do j=1,nrad
             tab1(i,j)=c(i,j,li)     ! li=1
           enddo
         enddo 
         
         return  

         end 
*________________________________________________________________________ 

         subroutine coagul

*********************************************************
* ce programme calcule la nouvelle concentration dans   *
* le a ieme intervalle de rayon, a l'altitude h, a      *
* l'instant t+dt                                        *
*********************************************************
         use dimphy
         IMPLICIT NONE
#include "dimensions.h"
#include "microtab.h"
#include "varmuphy.h"


* declaration des blocs communs
*------------------------------

         common/ctps/li,lf,dt
         common/con/c

* declaration des variables
* --------------------------

         integer li,lf
         real dt
         real c(nz,nrad,2)

* declaration des variables propres au ss-programme
* -------------------------------------------------

         integer h,a
         real pr,pe

*  traitement
*  ----------

         do h=nztop,nz
           do a=1,nrad
             call pertpro(h,a,pe,pr)
c            if((1+dt*pe).lt.0.) stop 'denom.eq.0'
             c(h,a,lf)=(c(h,a,li)+pr*dt)/(1+dt*pe)
           enddo
         enddo 

         if (nztop.ne.1) then
           do h=1,nztop-1
             do a=1,nrad
               c(h,a,lf)=c(h,a,li)
             enddo
           enddo 
         endif

         return
         end
      
      
*__________________________________________________________________________

         subroutine  calcoag

***************************************************************
*                                                             *
*  Ce programme calcule les coefficients de collection  d'une *
* particule de rayon x avec une particule de rayon b a une    *
* altitude donnee h                                           *
*************************************************************** 

* declaration des blocs communs
*------------------------------
         use dimphy
         IMPLICIT NONE
#include "dimensions.h"
#include "microtab.h"
#include "varmuphy.h"

         common/ctps/li,lf,dt
         common/con/c
         common/coag/k

* declaration des variables
* --------------------------

         integer li,lf,i
         real dt
         real c(nz,nrad,2)
         real knu,nud,k(nz,nrad,nrad)

* declaration des variables propres au ss-programme
* -------------------------------------------------

         integer h,b,x
         real nua,lambb,lambx,knb,knx,alphab,alphax,d,e,f,kcg
         real db,dx,rm,dm,deltab,deltax,del,g,beta,gx,gb
         real rfx,rfb,rpr
         real*8 ne,qe,epso
         real*8 corelec,yy

         real kco,vx,vb,vitesse,sto,ee,a,dd,bb,p0,t0,l0,ccol
         real st(37),ef(37)

* initialisation
* --------------
c        print*,'**** calcoag'

*  -nombres de STOCKES

         data(st(i),i=1,37)/1.35,1.5,1.65,1.85,2.05,2.25,2.5,2.8,3.1,
     s    3.35,3.6,3.95,4.3,4.7,5.05,5.45,5.9,6.4,7.,7.6,8.3,9.05,9.9,
     s       10.9,11.1,13.5,15.3,17.25,20.5,24.5,30.4,39.3,48,57,86.,
     s       187.,600./

*  -coef. d'efficacite de collection

         ef(1)=3.75
         ef(2)=8.75
         do i=3,37
           ef(i)=ef(i-1)+2.5
         enddo

         do i=1,37
           ef(i)=ef(i)*1e-2
         enddo

         qe=1.6e-19
         ne=-30.e+6         ! "vieille valeur!"
         ne=-15.e+6          ! pour fitter DISR !

         epso=1e-9/(36*pi)

         d=1.257
         e=0.4
         f=-1.1

*  iteration sur z
 
         do 1 h=1,nz
           nua=nud(h,1)      

*  iteration sur les rayons

           do 1 b=1,nrad

             knb=knu(h,b,1)
             vb=vitesse(h,b,1)

             do 1 x=1,b

               knx=knu(h,x,1)
               vx=vitesse(h,x,1)

**  COAGULATION  ****************************************************
** --------------****************************************************
* calcul du terme correcteur 'slip-flow'

               alphab=d+e*exp(f/knb)
               alphax=d+e*exp(f/knx)

* calcul du coefficient de diffusion

               rfb=(r_e(b)**(3./df(b)))*((rf(b))**(1.-3./df(b)))
               rfx=(r_e(x)**(3./df(x)))*((rf(x))**(1.-3./df(x)))
               db=kbz*t(h)*(1+alphab*knb)/(6*pi*nua*rfb)
               dx=kbz*t(h)*(1+alphax*knx)/(6*pi*nua*rfx)

* calcul du coefficient de coagulation

               rpr=rfb+rfx
               kcg=4*pi*rpr*(db+dx)

* calcul de la vitesse thermique

               gx=sqrt(6*kbz*t(h)/(rhol*pi**2*r_e(x)**3))
               gb=sqrt(6*kbz*t(h)/(rhol*pi**2*r_e(b)**3))

* calcul du libre parcours apparent des aerosols

               lambb=8*db/(pi*gb)
               lambx=8*dx/(pi*gx)

*calcul du terme correcteur beta

               rm=rpr/2.
               dm=(dx+db)/2.
               g=sqrt(gx**2+gb**2)
               deltab=(((2*rfb+lambb)**3-(4*rfb**2+lambb**2)**1.5)
     s         /(6*rfb*lambb)-2*rfb)*sqrt(2.)
               deltax=(((2*rfx+lambx)**3-(4*rfx**2+lambx**2)**1.5)
     s         /(6*rfx*lambx)-2*rfx)*sqrt(2.)
               del=sqrt(deltab**2+deltax**2)
               beta=1/((rm/(rm+del/2))+(4*dm/(g*rm)))

* calcul du coefficient de coagulation corrige

               kcg=kcg*beta
c       print*,' kcg:', rfb,rfx,knb,knx,db,dx
c       print*,' beta:', gx,gb,lambb,lambx,rm,dm
c       print*,' beta:', g,deltab,deltax,del,beta
c       print*,' '

**  COALESCENCE  **************************************************
**  -------------**************************************************

               kco=0.

               if (b.eq.x) goto 9

* calcul du nombre de Stockes de la petite particule

               sto=2*rhol*rfx**2*abs(vx-vb)/(9*nua*rfb)

* calcul du coef. de Cunningham-Millikan

               a=1.246
               bb=0.42
               dd=0.87
               l0=0.653e-7
               p0=101325.
               t0=288.

               ee=1+(l0*t(h)*p0*
     &         (a+bb*exp(-dd*rfx*t0*p(h)/(l0*t(h)*p0))))
     s         /(rfx*t0*p(h))

* calcul du nombre de Stockes corrige

               sto=sto*ee

               if (sto .le. 1.2) goto 9

               if (sto .ge. 600.) then 
                 ccol=1.
                 goto 8
               endif

*  recherche du coefficient de collection

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

*  calcul du coefficient de coalescence

8              kco=pi*(rfb+rfx)**2*ccol*abs(vb-vx)

9              continue

**  CORRECTION ELECTRICITE *******************************
**  ------------------------******************************

               yy=1.d0*ne**2*r_e(x)*r_e(b)*qe**2
     &         /(1.d0*kbz*t(h)*(r_e(b)+r_e(x))*4*pi*epso)


               corelec=0.
               if (yy.lt.50.) corelec=yy/(exp(yy)-1.)
               if (yy.le.1.e-3)  corelec=1.

               k(h,b,x)=(kcg+kco)*corelec
               k(h,x,b)=k(h,b,x)
 

1        continue
         return
         end

*______________________________________________________________________

         real function lambda(j,indic)
*
*------------------------------------------------------------------*
*  fonction calculant le libre parcours moyen des molecules        *
*  atmospheriques( rayon =ra) se trouvant dans la couche no j.     *
*  pour indic=0  ...... la particule se trouve a la frontiere entre*
*                        les couches j et j-1                      *
*  pour indic=1  ...... la particule se trouve au milieu de la     *
*                         la couche j                              *
*------------------------------------------------------------------*
*
* declaration des blocs communs
*------------------------------
         IMPLICIT NONE
#include "dimensions.h"
#include "microtab.h"
#include "varmuphy.h"

* declaration des variables communes
* ----------------------------------

         integer i,j

*  declaration des variables internes
*  ----------------------------------

         integer indic
         real pp,ra

         ra=1.75e-10

* traitement
* ----------

         if (indic.eq.0) then
           pp=pb(j)
         else
           if (indic.ne.1) then
             print*,'erreur argument fonction lambda'
             return
           endif
           pp=p(j)
         endif

         lambda=kbz*t(j)/(4*sqrt(2.)*pi*(ra**2)*pp)
         end

*******************************************************************************

         real function knu(j,k,indic)
*
*--------------------------------------------------------------*
*  fonction calculant le nombre de knudsen d'une particule     *
*  d'aerosol de rayon r_e(k) se trouvant dans la couche no j   *
*  indic ......  idem function lambda                          *
*--------------------------------------------------------------*
*
* declaration des blocs communs
*------------------------------
         IMPLICIT NONE
#include "dimensions.h"
#include "microtab.h"
#include "varmuphy.h"

*  declaration des variables internes
*  ----------------------------------

         integer indic,j,k
         real lambda,rfk

* traitement
* ----------

         if (indic.ne.0 .and.indic.ne.1) then
           print*,'erreur argument fonction knu'
           return
         endif

         rfk=(r_e(k)**(3./df(k)))*((rf(k))**(1.-3./df(k)))
         knu=lambda(j,indic)/rfk
         end

*****************************************************************************

         real function nud(j,indic)
*
*--------------------------------------------------------------*
*  fonction calculant la viscosite dynamique (en USI) de l'air *
*  d'apres la formule de Sutherlant a l'altitude j             *
*  indic  ......... idem fonction lambda                       *
*--------------------------------------------------------------*
*
         IMPLICIT NONE
#include "dimensions.h"
#include "microtab.h"
#include "varmuphy.h"


         integer indic,j
         real nud0,c,tt
*
         nud0=1.74e-5
         c=109.
 
         if(indic.ne.0.and.indic.ne.1) then
           print*,'erreur argument fonction nud'
           return
         endif
 
         if(indic.eq.0) tt=tb(j)
         if (indic.eq.1) tt=t(j)
 
         nud=nud0*sqrt(tt/293)*(1+c/293)/(1+c/tt)
         end

****************************************************************************

         real function vitesse(j,k,indic)
*
*-----------------------------------------------------------------*
*   fonction calculant la vitesse de chute d'une particule de rayon*
*   k se trouvant a l'altitude j  suivant la valeur du nombre de   *
*    Knudsen                                                       *
*   indic ....... idem function lambda                             *
*-----------------------------------------------------------------*
*
*  declaration des blocs communs
*------------------------------
         IMPLICIT NONE
#include  "dimensions.h"
#include  "microtab.h"
#include  "varmuphy.h"

*   declaration des variables internes
*   ----------------------------------

         integer indic,j,k
         real w,g,m,a0,zz,nud,knud,tt,rhoh
         real rbis, rfk,vb,zbx
         real akncx,knu

*   traitement
*   ----------

         if (indic.ne.0.and.indic.ne.1) then
           print*,'erreur argument fonction vitesse'
           return
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
         knud=knu(j,k,indic)


         rfk=(r_e(k)**(3./df(k)))*((rf(k))**(1.-3./df(k)))
c        rfk=r_e(7)


          w=2./9.*rfk**(df(k)-1.)*rf(k)**(3.-df(k))*g*rhol/nud(j,indic)

         w=w*(1+1.2517*knud+0.4*knud*exp(-1.1/knud))


c        if (p(j).lt.500..and.k.eq.nrad) then
c           w=0.
c        endif

         vitesse=w

         end
***********************************************************************

         real function kd(h)
*
*--------------------------------------------------------------------*
*   cette fonction calcule le coefficient du terme de 'eddy diffusion'*
*   a l altitude j						     *
*--------------------------------------------------------------------*
*
         IMPLICIT NONE
#include  "dimensions.h"
#include  "microtab.h"
#include  "varmuphy.h"

         real zbx

         integer h

         zbx=z(h)+dz(h)/2.

c ATTENTION !! 
c toutes ces definitions sont contradictoires, 
c pour mettre 0 au bout du compte...
c A NETTOYER !!

         if(zbx.le.42000.) then
           kd=1.64e+12*(pb(h)/(kbz*tb(h)))**(-1./2.)
           kd=4.
         else
           kd=1.64e+12*(pb(h)/(kbz*tb(h)))**(-1./2.)
         endif

         if(zbx.le.50000.) then
           kd=1.64e+12*(pb(h)/(kbz*tb(h)))**(-1./2.)
         endif

         kd=0.0*kd

         return
         end


*____________________________________________________________________________

         subroutine init
*
*--------------------------------------------------------------------*
*   cette routine effectue  :                                         *
*                1) interpolation a partir des donnees initiales des  *
*                    valeurs de p,t,rho,ach4,aar,an2  sur la grille   *
*                2) initialisation des constantes (common/phys/)      *
*                3) initialisation des variables temporelles (common  *
*                     /temps/)                                        *
*                4) definition des grilles en rayon et verticale      *
*                5)  initialisation de c(z,r,t) avec les donnees du   *
*                      fichier unit=1                                 *
*                                                                     *
*   les donnees sont des valeurs caracterisques de l'atmosphere de    *
*     TITAN  ( voir Lelouch and co )                                  *
*--------------------------------------------------------------------*

*  declaration des blocs communs
*------------------------------
         use dimphy
         IMPLICIT NONE
#include  "dimensions.h"
#include  "microtab.h"
#include  "varmuphy.h"


         common/ctps/li,lf,dt
         common/con/c

*  declaration des variables communes
*  ----------------------------------

         real c(nz,nrad,2)
         integer li,lf
         integer i,ii
         real dt

*  declaration des variables internes
*  ----------------------------------

         integer nzd
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

!        do 20 j1=1,nz
!           do 21 j2=1,nzd
!               if( zd(j2).le.z(j1)) goto 22
!21	    continue
!22         liminf=j2

          do 20 j1=1,nz
             do 21 j2=1,nzd
                 if( zd(j2).le.z(j1)) goto 22
  21         continue
  22         if (j2.ge.254) j2=254
             liminf = j2

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
     s       (zd(limsup)-zd(liminf))*(zd(liminf)-z(j1))
          endif
20      continue

*   rap= aar/an2  cst sur l'altitude

         rap=0.02
c         rap=0.191
         do 23 i=1,nz
           an2(i)=(1.-ach4(i))/(1.+rap)
           aar(i)=rap*an2(i)
23       continue

         do 24 i=1,nz
           m=ach4(i)*mch4+an2(i)*mn2+aar(i)*mar
           rho(i)=p(i)*m/(rgp*t(i))
24       continue

         do i=1,nz
           m=ach4(i)*mch4+an2(i)*mn2+aar(i)*mar
           rhob(i)=pb(i)*m/(rgp*tb(i))
c          print*,pb(i),m,rgp,tb(i),rhob(i),rho(i)
         enddo

*  fin d'interpolation des taux de melange
*----------------------------------------  

c        print*,'**** fin init'
540      continue
         return

500       print*,'erreur lecture initialisation de c...erreur=',ii
          stop

        end

*____________________________________________________________________________

         subroutine pertpro(h,a,l_,pr_)

*****************************************************************************
*                                                                           *
*  ce programme permet le calcul du terme de production (pr) et de perte (l)*
*  pour le phenomene de coagulation                                         *
*  dans le a ieme intervalle de rayon a une altitude h                      *
*****************************************************************************

*  declaration des blocs communs
*------------------------------
         use dimphy
         IMPLICIT NONE
#include  "dimensions.h"
#include  "microtab.h"
#include  "varmuphy.h"


         common/ctps/li,lf,dt
         common/con/c
         common/coag/k

*  declaration des variables
*  --------------------------

         integer li,lf
         real dt
         real c(nz,nrad,2),k(nz,nrad,nrad)

*  declaration des variables propres au ss-programme
*  -------------------------------------------------

         integer h,b,a,x,i
         real*8 pr,ss,s,l
         real pr_,l_,vol,del

*  traitement
*  -----------

*   production
*+++++++++++++
         s=0.d0
         ss=0.d0
         pr=0.

         if (a .eq. 1) goto 2  
         b=a-1

         if (c(h,b,lf) .eq. 0 .and. c(h,b,li) .eq. 0) goto 2
  
         do 1 i=1,b

         if(c(h,i,li) .eq. 0 .and. c(h,i,lf) .eq. 0) goto 1
     
         if (i .ne. b)del=1.
         if (i .eq. b) del=.5

         s=(v_e(i)*1.d0)*del*(k(h,b,i)*1.d0)*(c(h,i,li)*1.d0)+s
         ss=(v_e(i)*1.d0)*del*(k(h,b,i)*1.d0)*(c(h,i,lf)*1.d0)+ss

c        if (a.eq.2) print*,'SS>',v_e(i),k(h,b,i),c(h,b,lf)
c        if (a.eq.2) print*,'SS>',del*v_e(i)*k(h,b,i)*c(h,b,lf)


1        continue

*   calcul du terme de production

          pr=(c(h,b,lf)*s/(vrat_e-1.)+c(h,b,li)*ss)/v_e(a)
c         if (a.eq.2) print*,'PR>',s,ss,c(h,b,lf),v_e(a)

2         continue


*   perte
*-  - - - -

          l=0


*     condition limite : pas de perte dans le dernier intervalle

         if (a .eq. nrad) goto 9

         do 10 x=1,nrad

         if (c(h,x,li) .eq. 0) goto 10

         if (a .lt. x) vol=1.
         if (a .eq. x) vol=.5*vrat_e/(vrat_e-1)
         if (a .gt. x) vol=v_e(x)/(v_e(a)*(vrat_e-1))

         l=l+k(h,a,x)*c(h,x,li)*vol*1.d0

      
10       continue
9        continue

#ifdef  CRAY
         l_=l
         pr_=pr
#else
         l_=sngl(l)
         pr_=sngl(pr)
#endif
c        l_=sngl(l)
c        pr_=sngl(pr)
c          if (a.eq.2) print*,'pr_,l_',h,a,pr_,l_
c          if (a.eq.2) print*,'-----------------------'
c          if (a.eq.2) STOP



         return

         end

*_____________________________________________________________________________

         subroutine production(ihor)
*
*--------------------------------------------------------------------*
*   routine calculant le terme de production des molecules organiques *
*   composant les aerosols . rini= rayon des aerosols initiaux        *
*--------------------------------------------------------------------*
*
         use dimphy
         IMPLICIT NONE
#include  "dimensions.h"
#include  "microtab.h"
#include  "varmuphy.h"
#include  "clesphys.h"

c#include  "aerprod.h"


         integer ndz,i,ihor,k,k1
         real c(nz,nrad,2)
         integer li,lf
         real dt
         real zprod,zy,c0,ctot,prod,rini,rfron
         real xsaison,p0
         common/ctps/li,lf,dt
         common /con/c
         common/effets/xsaison

! Pressure level of aerosol production
          p0=p_prodaer

          do i=1,nz-1
           if (pb(i).lt.p0.and.pb(i+1).gt.p0) zprod=(z(i)+z(i+1))/2.
          enddo

         ctot=3.5e-13*tx ! ATTENTION, ??COHERENT AVEC INITPAR??
         ctot=ctot*xsaison  ! 

c        z0=385.e+3
         zy=20.e+3
         rini=1.3e-9
         rini=r_e(1)   
         ndz=50
*
         do 10 i=1,nrad 
           if(rini.lt.r_e(i)) goto 100
10       continue
100      continue
         if (i.eq.1) then
           rini=r_e(1)
         else
           rfron=(r_e(i)+r_e(i-1))/2
           if (rini .lt.rfron) then 
             rini=r_e(i-1)
             i=i-1
           else
             rini=r_e(i)
           endif
         endif
*
         c0=ctot/(sqrt(2.*pi)*zy)
         c0=c0*3./(4.*pi*rhol*rini**3)
*
         do 20 k=nztop,nz
         prod=0.
         do 201 k1=1,ndz
           prod=prod+c0*exp(-0.5*(((z(k)+dz(k)/2.-k1*dz(k)/(2.*ndz)
     s     -zprod)/zy)**2))*dt/ndz
201      continue

          if (prod .le. 1) prod=0.
            c(k,i,lf)=c(k,i,lf)+prod
20       continue

!        do 30 k=nztop,nz
!        c(k,i,lf)=c(k,i,lf)+prodaer(ihor,nz-k+1,1)
!     & *3./(4.*pi*rhol*rini**3)
!30      continue

         return
         end

*-------------------------------------------------------------------*



        subroutine sedif(dice3)
*
*------------------------------------------------------------------*
*   cette routine calcule l'evolution de la fonction de distribution*
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

         common/ctps/li,lf,dt
         common/con/c

*  declaration des variables communes
*  ----------------------------------

         integer li,lf
         integer i,j,k,nb
         real dt
         real c(nz,nrad,2),dice3,bilan4,bilan14

*  declaration des variables internes
*  ----------------------------------

         real w,w1,dzbX,dc
         double precision sigma,theta,hc,l,rap,cmp,wp 
         double precision fs(nz+1),ft(nz+1)
         real as(nz),bs(nz),cs(nz),ds(nz)
         double precision asi(nztop:nz),bsi(nztop:nz),
     &                    csi(nztop:nz)
         double precision dsi(nztop:nz),xsol(nztop:nz)
         real vitesse,kd

         external dtridgl
*  resolution
*------------

            
            bilan4=0.
            do  k=1,nrad
            do  j=nztop,nz
              bilan4=bilan4+c(j,k,li)*dzb(j)*
     &        4./3.*pi*rf(k)**3.*vrat_e**(k-imono)
            enddo
            enddo


            do 10 k=1,nrad  
            do 20 j=nztop,nz

            if (j.eq.1) goto 20

*  calcul de la vitesse corrigee

              dzbX=(dz(j)+dz(j-1))/2.
              w=-1*vitesse(j,k,0)
              if (kd(j).ne.0.) then
                theta=0.5*(w*dzbX/kd(j)+log(rho(j-1)/rho(j)))
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

c               if(w1.ge.0. .or. w1.le.0.) then 
c                continue
c               else
c                print*,j,k,w1,hc,dzbx,rap,dlog(rap),' w1'
c               endif

*   calcul des flux aux interfaces


             if (kd(j).ne.0.) then
               if (theta.ne.0.) then
                 ft(j)=(w1+log(rho(j-1)/rho(j))*kd(j)/dzbX)/(dexp(2.*
     s           theta)-1.)
                 fs(j)=ft(j)*dexp(2.*theta)
               else
                 ft(j)=kd(j)/dzbX
                 fs(j)=kd(j)/dzbX
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

20          continue

*  conditions aux limites pour les flux aux interfaces

            fs(1)=0.
            ft(1)=0.
            fs(nz+1)=0.
            ft(nz+1)=-w1

*  calcul des coefficients de l'equation discrete

            do 21 j=nztop,nz
              as(j)=-dz(j)/dt
              bs(j)=-ft(j)
              cs(j)=ft(j+1)+fs(j)-dz(j)/dt
              ds(j)=-fs(j+1)

              if ( cs(j).gt.0) goto 100
21         continue

*  cas explicite (mu=0) : calcul de la fonction c(z,r,t+1)

           do 22 j=nztop,nz-1

           if (j.eq.nztop) then
             dc=(cs(nztop)*c(nztop,k,li)+ds(nztop)
     &                    *c(nztop+1,k,li))/as(nztop)
             c(nztop,k,lf)=dc
             goto 22
           endif

           dc=(bs(j)*c(j-1,k,li)+cs(j)*c(j,k,li)+ds(j)*c(j+1,k,li))
     s         /as(j)
           c(j,k,lf)=dc


22          continue

            dc=(bs(nz)*c(nz-1,k,li)+cs(nz)*c(nz,k,li))/as(nz)
            c(nz,k,lf)=dc

           if (nztop.ne.1) then
           do 32 j=1,nztop-1
             c(j,k,lf)=c(j,k,li)
32         continue
           endif

           goto 10

100        continue

*  cas implicite (mu=1) : calcul de la fonction c(z,r,t+1)

          do 101 j=nztop,nz
            asi(j)=ft(j)
            bsi(j)=-(ft(j+1)+fs(j)+dz(j)/dt)
            csi(j)=fs(j+1)
            dsi(j)=-dz(j)/dt*c(j,k,li)
            xsol(j)=0.
101       continue

*  inversion de la matrice tridiagonale 

         nb=nz-nztop+1

         call dtridgl(nb,asi,bsi,csi,dsi,xsol) 


         do 102 j=nztop,nz
          c(j,k,lf)=xsol(j)
102      continue

         if (nztop.ne.1) then
           do 110 j=1,nztop-1
           c(j,k,lf)=c(j,k,li)
110        continue
         endif



10       continue

         bilan14=0.
         do  k=1,nrad
         do  j=nztop,nz
           bilan14=bilan14+c(j,k,lf)*dzb(j)*
     &     4./3.*pi*rf(k)**3.*vrat_e**(k-imono)
         enddo
         enddo


         dice3=(bilan14-bilan4)*rhol


         return

          end

