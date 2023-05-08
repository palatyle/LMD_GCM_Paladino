         subroutine cnuages(
     &   tq,tqc1,tqc2,tqc3,tqcn,gaz1,gaz2,gaz3,   ! aerosol/glace/gas
     &   ddt)

c
c
c  °
c              °
c  SERT A APPELE LA ROUTINE MICROPHYSIQUE DES NUAGES 
c
c     °
c  ICI ON NE FAIT QUE LA NUCLEATION/CONDENSATION 
c  ET GESTION°DES NOYAUX. LA SEDIMENTATION EST DANS
c  SNUAGES.F
c            °   °
c
c        °
c                °
c              °   °
c                      
c                  °   °
c                   °  
c                    ° °
c                   ° \|/  °
c                    (@ @)°
c-----------------oOo--O--oOo--------------------------
c                        
c
c
c Interface entre physiq.F et les routines n_<nom_compose>.F 
c
c Date: 3 Nov 2003
c
c   EN ENTREE/SORTIE DE LA ROUTINE
c ------------------------------------
c
c     Les aerosols, noyaux (tq,tqcn) sont en nbre/m^2 dans la colonne
c     Les  condensats (tqc1,tqc2)   sont en volume/m^2 dans la colonne
c     Le gaz (gaz1, gaz2) est en fraction molaire
c
c   EN APPEL DES ROUTINES NUAGES
c ------------------------------------
c
c     Les aerosols et noyaux doivent etre en nombre /kg d'air
c     Les condensats doivent etre en volume / kg d'air
c     Le gaz en kg/kg d'air
c
c   LES TENDANCES ET DIFFERENCES SONT HOMOGENES AUX QUANTITES
c ------------------------------------------------------------
c
c
c------------------------------------------------------

         use dimphy
         IMPLICIT NONE
#include "dimensions.h"
#include "microtab.h"
#include "varmuphy.h"

         integer NG1,NG,NL 
         parameter (NG1=1,NG=NG1,NL=llm)
  
c*************************************
c  declaration des variables internes *
c*************************************

c        INTERNE!        *
c-----------------------*

         real  tqc1(NG,NL,nrad)
         real  tqc2(NG,NL,nrad)
         real  tqc3(NG,NL,nrad)
         real  tqcn(NG,NL,nrad)
         real  tq(NG,NL,nrad)
*
         real  tdqc1(NG,NL,nrad)
         real  tdqc2(NG,NL,nrad)
         real  tdqc3(NG,NL,nrad)
         real  tdq(NG,NL,nrad,ntype-2+1)
         real  tdqcn(NG,NL,nrad,ntype-2+1)
*
         real  gaz1(NG,NL)
         real  gaz2(NG,NL)
         real  gaz3(NG,NL)
         real  dgaz1(NG,NL)
         real  dgaz2(NG,NL)
         real  dgaz3(NG,NL)
*
         real  ppch4(NG,NL)
         real  ppc2h6(NG,NL)
         real  ppn2(NG,NL)
*
         real  pmixch4(NL)
         real  pmixc2h6(NL)
         real  pmixn2(NL)
c   composition initiale estimée (interne)
         real  xprime1(NG,NL)
         real  xprime2(NG,NL)
         real  xprime3(NG,NL)
c   composition calculée (output)
         real  x1(NL)
         real  x2(NL)
         real  x3(NL)
c   moyenne "glissante" pondéréee (output + mémoire)
         real  x1o(NL)
         real  x2o(NL)
         real  x3o(NL)
         real  icefrac(NL)
         real  dmn2(NL+1)

         real  ppch4t,ppc2h6t,ppn2t
         real  psatch4,psatc2h6,psatn2
         real  xprime(3),x(3),frac
         real  melange
         real  sum,sum0

*                    RAPPEL: NG=1

         real ddt
         real masspaer

         common/mixing/x1,x2,x3,icefrac,
     &         pmixch4,pmixc2h6,pmixn2,
     &         x1o,x2o,x3o



c  FORMAT MICRO DES NUAGES
c------------------------*

         real  especes(NG,NL,3*nrad+1)          
         real  condens(NG,NL,nrad)          
         real  gg(NL),xmair
         real  effg     ! effg est une fonction(z), z en m.

         integer jsup,jinf,h,i,j,k,ndim
         integer ival1,ival2,ival3

         integer iprem

         save iprem,xprime
         data iprem/0/

         ndim=3*nrad+1

         do j=1,NL
           gg(j)=effg(z(j))
         enddo

*********************************************
*********************************************
*  Appel de la condensation du methane
*********************************************
*********************************************



*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
*                   Bilan avant sur le methane                         *
*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

         do i=1,ng1  ! ng1=1 !!
           do j=1,NL

*  RAZ des delta d'especes communes.
*----------------------------------
             do k=1,nrad
               tdqcn(i,j,k,1)= 0.
               tdq( i,j,k,1) = 0.
             enddo

             xmair=(pb(j+1)-pb(j))/gg(j)/dzb(j)

             do k=1,nrad
               especes(i,j,k)=tq(i,j,k)         /xmair       ! aerosols, noyaux,
               especes(i,j,k+nrad)=tqc1(i,j,k)  /xmair       ! methane condense,
               condens(i,j,k)=(tqc2(i,j,k)+tqc3(i,j,k))/xmair! autre(s) condensat(s)
               especes(i,j,k+2*nrad)=tqcn(i,j,k)/xmair       ! nombre/kg &  volume/kg
             enddo
             especes(i,j,3*nrad+1)=gaz1(i,j)*mch4/mair       ! methane gazeux kg/kg
           enddo
         enddo


1001   format(7(1x,e12.6),' avN2CH4C2H6')
1003   format(7(1x,e12.6),' miN2CH4C2H6')
1002   format(7(1x,e12.6),' apN2CH4C2H6')
*
         call n_methane(ng1,ndim,nrad,ddt,
     &                  p,t,r_e,especes,condens)


*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
*                   Bilan apres sur le methane                         *
*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

         do i=1,ng1
           do j=1,NL

             xmair=(pb(j+1)-pb(j))/gg(j)/dzb(j)

*  ici ce sont les tendances a sortir de nuages.F pour le methane....
*-------------------------------------------------------------------
             sum=0.   
             do k=1,nrad
               tdqc1(i,j,k)=(especes(i,j,k+nrad)*xmair-tqc1(i,j,k) )
               sum =sum+tdqc1(i,j,k)/xmair*rhoi_ch4
               tqc1(i,j,k) = especes(i,j,k+nrad)*xmair
             enddo

             dgaz1(i,j)= especes(i,j,3*nrad+1)*mair/mch4-gaz1(i,j)
             gaz1(i,j)= especes(i,j,3*nrad+1)*mair/mch4

c         dgaz1(i,j)=-sum*xmuair/16.
c          gaz1(i,j)=gaz1(i,j)+dgaz1(i,j) 

*  Premiere tendance sur les variables communes (aerosols et noyaux)
*------------------------------------------------------------------
             do k=1,nrad
               tdqcn(i,j,k,1)=(especes(i,j,k+2*nrad)*xmair-tqcn(i,j,k))
               tdq( i,j,k,1) =(especes(i,j,k)*xmair       -tq(i,j,k))  
             enddo
           enddo
         enddo



*     attention, si il y a  de l'ethane sur les noyaux... il est impossible
*     de les restituer - en revanche on peut en creer de nouveaux ! !!
*     Le corrolaire de la condition ci dessus est que si il est impossible de 
*     restituer des noyaux, le nombre d'aeorsols ne peux pas augmenter, il
*     peut en revanche diminuer


*********************************************
*********************************************
*  Appel de la condensation de l'ethane
*********************************************
*********************************************

         do i=1,ng1
           do j=1,NL
             do k=1,nrad
               tdqcn(i,j,k,2)= 0.
               tdq( i,j,k,2) = 0.
             enddo

             xmair=(pb(j+1)-pb(j))/gg(j)/dzb(j)

             do k=1,nrad
               especes(i,j,k)=tq(i,j,k)/xmair
               especes(i,j,k+nrad)=tqc2(i,j,k)/xmair           ! ethane condense
               condens(i,j,k)=(tqc1(i,j,k)+tqc3(i,j,k))/xmair  ! autres condensats
               especes(i,j,k+2*nrad)=tqcn(i,j,k)/xmair
             enddo

             especes(i,j,3*nrad+1)=gaz2(i,j)*mc2h6/mair       ! ethane gazeux

           enddo
         enddo



         call n_ethane(ng1,ndim,nrad,ddt,
     &                 p,t,r_e,especes,condens)

         do i=1,ng1
           do j=1,NL
             xmair=(pb(j+1)-pb(j))/gg(j)/dzb(j)

*  ici ce sont les tendances a sortir de nuages.F pour l'ethane....
*-----------------------------------------------------------------

             sum=0.
             do k=1,nrad
               tdqc2(i,j,k)=(especes(i,j,k+nrad)*xmair-tqc2(i,j,k) )
               sum =sum+tdqc2(i,j,k)/xmair*rhoi_c2h6
               tqc2(i,j,k) = especes(i,j,k+nrad)*xmair
             enddo

             dgaz2(i,j)=especes(i,j,3*nrad+1)*mair/mc2h6 - gaz2(i,j)
             gaz2(i,j)=especes(i,j,3*nrad+1)*mair/mc2h6 

c         dgaz2(i,j)=-sum*xmuair/30.
c          gaz2(i,j)=gaz2(i,j)+dgaz2(i,j) 


*  Deuxieme tendance sur les variables communes (aerosols et noyaux)
*------------------------------------------------------------------
   
             do k=1,nrad
               tdqcn(i,j,k,2)=(especes(i,j,k+2*nrad)*xmair-tqcn(i,j,k))
               tdq(i,j,k,2)  =(especes(i,j,k)*xmair       -tq(i,j,k))
             enddo

           enddo
         enddo


*********************************************
*********************************************
*  Appel de la condensation de l'acethylene
*********************************************
*********************************************

         do i=1,ng1
           do j=1,NL
             do k=1,nrad
               tdqcn(i,j,k,3)= 0.
               tdq( i,j,k,3) = 0.
             enddo
             xmair=(pb(j+1)-pb(j))/gg(j)/dzb(j)

             do k=1,nrad
               especes(i,j,k)=tq(i,j,k)/xmair
               especes(i,j,k+nrad)=tqc3(i,j,k)/xmair           ! acethylene condense
               condens(i,j,k)=(tqc1(i,j,k)+tqc2(i,j,k))/xmair  ! autres condensats
               especes(i,j,k+2*nrad)=tqcn(i,j,k)/xmair
             enddo

             especes(i,j,3*nrad+1)=gaz3(i,j)*mc2h2/mair      ! acethylene gazeux

           enddo
         enddo

         call n_acethylene(ng1,ndim,nrad,ddt,
     &                     p,t,r_e,especes,condens)


         do i=1,ng1
           do j=1,NL

             xmair=(pb(j+1)-pb(j))/gg(j)/dzb(j)

*  ici ce sont les tendances a sortir de nuages.F pour l'ethane....
*-----------------------------------------------------------------

             sum=0.
             do k=1,nrad
               tdqc3(i,j,k)=(especes(i,j,k+nrad)*xmair-tqc3(i,j,k) )
               sum =sum+tdqc3(i,j,k)/xmair*rhoi_c2h2
               tqc3(i,j,k) = especes(i,j,k+nrad)*xmair
             enddo

             dgaz3(i,j)=especes(i,j,3*nrad+1)*mair/mc2h2 - gaz3(i,j)
             gaz3(i,j)=especes(i,j,3*nrad+1)*mair/mc2h2 

c         dgaz3(i,j)=-sum*xmuair/26.
c         gaz3(i,j)=gaz3(i,j)+dgaz3(i,j) 

*  Troisieme tendance sur les variables communes (aerosols et noyaux)
*------------------------------------------------------------------


             do k=1,nrad
               tdqcn(i,j,k,3)=(especes(i,j,k+2*nrad)*xmair-tqcn(i,j,k))
               tdq(i,j,k,3)  =(especes(i,j,k)*xmair       -tq(i,j,k))
             enddo

           enddo
         enddo


*  FIN DES APPELS DE NUAGES ET BILAN DES TENDANCES...
*------------------------------------------------------------------


         do i=1,ng1
           do j=1,NL
             do k=1,nrad

*          Ici on test l'activité nuageuse : si on a l'association
*          tdqcX(i,j,k) = 0 et tqcX(i,j,k) = 0 alors ivalX reste à 0 (pas d'ativité)
*                                              sinon ivalX passe à 1 (activité)
*------------------------------------------------------------------------------------

               ival1=0
               ival2=0
               ival3=0

               if(tdqc1(i,j,k).ne.0. .or. tqc1(i,j,k).gt.0.)  ival1=1
               if(tdqc2(i,j,k).ne.0. .or. tqc2(i,j,k).gt.0.)  ival2=1
               if(tdqc3(i,j,k).ne.0. .or. tqc3(i,j,k).gt.0.)  ival3=1

*          Ici on definit la tendances des noyaux en faisant deux choses: 
*          -1 On ecarte les cas  tdqcn(i,j,k,X)=0 si ils sont associés à une 
*          absence d'activité nuageuse de l'espèce (tdqcX(i,j,k)=0.)
*          -2 Sélectionne la tendance la plus élevée. Si aucune activité nuageuse
*           n'exits dans cette case (ivalX=0 pour les 3 especes), alors on
*           retrouve la valeur -1.e40 que l'on mets alors à 0. 
*----------------------------------------------------------------------
c23456789012345678901234567890123456789012345678901234567890123456789012

               tdqcn(i,j,k,ntype-1)=-1.e40                      ! plus petite valeur possible 

               if(ival1.eq.1.and.
     &         tdqcn(i,j,k,1).ge.tdqcn(i,j,k,ntype-1)) ! Si activité de l'espece 1
     &         tdqcn(i,j,k,ntype-1)=tdqcn(i,j,k,1)
               if(ival2.eq.1.and.
     &         tdqcn(i,j,k,2).ge.tdqcn(i,j,k,ntype-1)) ! Si activité de l'espece 2
     &         tdqcn(i,j,k,ntype-1)=tdqcn(i,j,k,2)
               if(ival3.eq.1.and.
     &         tdqcn(i,j,k,3).ge.tdqcn(i,j,k,ntype-1)) ! Si activité de l'espece 3
     &         tdqcn(i,j,k,ntype-1)=tdqcn(i,j,k,3)

               if(tdqcn(i,j,k,ntype-1).le.-0.99e39) 
     &         tdqcn(i,j,k,ntype-1)=0.

               tdq(i,j,k,ntype-1)=1.e40                    ! plus grande valeur possible 

               if(ival1.eq.1 .and. tdq(i,j,k,1).le.tdq(i,j,k,ntype-1)) ! Si activité de l'espece 1
     &           tdq(i,j,k,ntype-1)=tdq(i,j,k,1)
               if(ival2.eq.1 .and. tdq(i,j,k,2).le.tdq(i,j,k,ntype-1)) ! Si activité de l'espece 2
     &           tdq(i,j,k,ntype-1)=tdq(i,j,k,2)
               if(ival3.eq.1 .and. tdq(i,j,k,3).le.tdq(i,j,k,ntype-1)) ! Si activité de l'espece 3
     &           tdq(i,j,k,ntype-1)=tdq(i,j,k,3)

               if(tdq(i,j,k,ntype-1).ge.0.99e39) tdq(i,j,k,ntype-1)=0.

                tqcn(i,j,k)=tqcn(i,j,k)+tdqcn(i,j,k,ntype-1)           ! Alors on ajoute les tendances (positive pour qcn ?)
                tq(i,j,k)  =tq(i,j,k)  +tdq(i,j,k,ntype-1)             !

                if(tqcn(i,j,k).le.0.) tqcn(i,j,k)=0.                   ! et on régularise les tableaux noyaux et aerosols.
                if(tq(i,j,k)  .le.0.) tq(i,j,k)=0.                     !

             enddo
           enddo
         enddo

        continue

1202    format(i2,1x,i2,6(1x,e12.4) )
       return 
       end
