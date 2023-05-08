










!
! $Id: interpre.F 1403 2010-07-01 09:02:53Z fairhead $
!
       subroutine interpre(q,qppm,w,fluxwppm,masse,
     s            apppm,bpppm,massebx,masseby,pbaru,pbarv,
     s            unatppm,vnatppm,psppm)

      USE control_mod
      USE comvert_mod, ONLY: ap,bp
      USE comconst_mod, ONLY: g

       implicit none

!-----------------------------------------------------------------------
!   INCLUDE 'dimensions.h'
!
!   dimensions.h contient les dimensions du modele
!   ndm est tel que iim=2**ndm
!-----------------------------------------------------------------------

      INTEGER iim,jjm,llm,ndm

      PARAMETER (iim= 128,jjm=96,llm=23,ndm=1)

!-----------------------------------------------------------------------
c#include "paramr2.h"
!
! $Header$
!
!
!  ATTENTION!!!!: ce fichier include est compatible format fixe/format libre
!                 veillez  n'utiliser que des ! pour les commentaires
!                 et  bien positionner les & des lignes de continuation
!                 (les placer en colonne 6 et en colonne 73)
!
!
!-----------------------------------------------------------------------
!   INCLUDE 'paramet.h'

      INTEGER  iip1,iip2,iip3,jjp1,llmp1,llmp2,llmm1
      INTEGER  kftd,ip1jm,ip1jmp1,ip1jmi1,ijp1llm
      INTEGER  ijmllm,mvar
      INTEGER jcfil,jcfllm

      PARAMETER( iip1= iim+1,iip2=iim+2,iip3=iim+3                       &
     &    ,jjp1=jjm+1-1/jjm)
      PARAMETER( llmp1 = llm+1,  llmp2 = llm+2, llmm1 = llm-1 )
      PARAMETER( kftd  = iim/2 -ndm )
      PARAMETER( ip1jm  = iip1*jjm,  ip1jmp1= iip1*jjp1 )
      PARAMETER( ip1jmi1= ip1jm - iip1 )
      PARAMETER( ijp1llm= ip1jmp1 * llm, ijmllm= ip1jm * llm )
      PARAMETER( mvar= ip1jmp1*( 2*llm+1) + ijmllm )
      PARAMETER( jcfil=jjm/2+5, jcfllm=jcfil*llm )

!-----------------------------------------------------------------------
!
! $Header$
!
!-----------------------------------------------------------------------
! INCLUDE comdissip.h

      COMMON/comdissip/                                                 &
     &    niterdis,coefdis,tetavel,tetatemp,gamdissip


      INTEGER niterdis

      REAL tetavel,tetatemp,coefdis,gamdissip

!-----------------------------------------------------------------------
!
! $Header$
!
!CDK comgeom2
      COMMON/comgeom/                                                   &
     & cu(iip1,jjp1),cv(iip1,jjm),unscu2(iip1,jjp1),unscv2(iip1,jjm)  , &
     & aire(iip1,jjp1),airesurg(iip1,jjp1),aireu(iip1,jjp1)           , &
     & airev(iip1,jjm),unsaire(iip1,jjp1),apoln,apols                 , &
     & unsairez(iip1,jjm),airuscv2(iip1,jjm),airvscu2(iip1,jjm)       , &
     & aireij1(iip1,jjp1),aireij2(iip1,jjp1),aireij3(iip1,jjp1)       , &
     & aireij4(iip1,jjp1),alpha1(iip1,jjp1),alpha2(iip1,jjp1)         , &
     & alpha3(iip1,jjp1),alpha4(iip1,jjp1),alpha1p2(iip1,jjp1)        , &
     & alpha1p4(iip1,jjp1),alpha2p3(iip1,jjp1),alpha3p4(iip1,jjp1)    , &
     & fext(iip1,jjm),constang(iip1,jjp1), rlatu(jjp1),rlatv(jjm),      &
     & rlonu(iip1),rlonv(iip1),cuvsurcv(iip1,jjm),cvsurcuv(iip1,jjm)  , &
     & cvusurcu(iip1,jjp1),cusurcvu(iip1,jjp1)                        , &
     & cuvscvgam1(iip1,jjm),cuvscvgam2(iip1,jjm),cvuscugam1(iip1,jjp1), &
     & cvuscugam2(iip1,jjp1),cvscuvgam(iip1,jjm),cuscvugam(iip1,jjp1) , &
     & unsapolnga1,unsapolnga2,unsapolsga1,unsapolsga2                , &
     & unsair_gam1(iip1,jjp1),unsair_gam2(iip1,jjp1)                  , &
     & unsairz_gam(iip1,jjm),aivscu2gam(iip1,jjm),aiuscv2gam(iip1,jjm)  &
     & , xprimu(iip1),xprimv(iip1)


      REAL                                                               &
     & cu,cv,unscu2,unscv2,aire,airesurg,aireu,airev,apoln,apols,unsaire &
     & ,unsairez,airuscv2,airvscu2,aireij1,aireij2,aireij3,aireij4     , &
     & alpha1,alpha2,alpha3,alpha4,alpha1p2,alpha1p4,alpha2p3,alpha3p4 , &
     & fext,constang,rlatu,rlatv,rlonu,rlonv,cuvscvgam1,cuvscvgam2     , &
     & cvuscugam1,cvuscugam2,cvscuvgam,cuscvugam,unsapolnga1           , &
     & unsapolnga2,unsapolsga1,unsapolsga2,unsair_gam1,unsair_gam2     , &
     & unsairz_gam,aivscu2gam,aiuscv2gam,cuvsurcv,cvsurcuv,cvusurcu    , &
     & cusurcvu,xprimu,xprimv

c---------------------------------------------------
c Arguments     
      real   apppm(llm+1),bpppm(llm+1)
      real   q(iip1,jjp1,llm),qppm(iim,jjp1,llm)
c---------------------------------------------------
      real   masse(iip1,jjp1,llm) 
      real   massebx(iip1,jjp1,llm),masseby(iip1,jjm,llm)      
      real   w(iip1,jjp1,llm)
      real   fluxwppm(iim,jjp1,llm)
      real   pbaru(iip1,jjp1,llm )
      real   pbarv(iip1,jjm,llm)
      real   unatppm(iim,jjp1,llm)
      real   vnatppm(iim,jjp1,llm)
      real   psppm(iim,jjp1)
c---------------------------------------------------
c Local
      real   vnat(iip1,jjp1,llm)
      real   unat(iip1,jjp1,llm)
      real   fluxw(iip1,jjp1,llm)
      real   smass(iip1,jjp1)
c----------------------------------------------------
      integer l,ij,i,j

c       CALCUL DE LA PRESSION DE SURFACE
c       Les coefficients ap et bp sont passés en common 
c       Calcul de la pression au sol en mb optimisée pour 
c       la vectorialisation
                   
         do j=1,jjp1
             do i=1,iip1
                smass(i,j)=0.
             enddo
         enddo

         do l=1,llm
             do j=1,jjp1
                 do i=1,iip1
                    smass(i,j)=smass(i,j)+masse(i,j,l)
                 enddo
             enddo
         enddo
      
         do j=1,jjp1
             do i=1,iim
                 psppm(i,j)=smass(i,j)/aire(i,j)*g*0.01
             end do
         end do                        
       
c RECONSTRUCTION DES CHAMPS CONTRAVARIANTS
c Le programme ppm3d travaille avec les composantes
c de vitesse et pas les flux, on doit donc passer de l'un à l'autre
c Dans le même temps, on fait le changement d'orientation du vent en v
      do l=1,llm
          do j=1,jjm
              do i=1,iip1
                  vnat(i,j,l)=-pbarv(i,j,l)/masseby(i,j,l)*cv(i,j)             
              enddo
          enddo
          do  i=1,iim
          vnat(i,jjp1,l)=0.
          enddo
          do j=1,jjp1
              do i=1,iip1
                  unat(i,j,l)=pbaru(i,j,l)/massebx(i,j,l)*cu(i,j)
              enddo
          enddo
      enddo
              
c CALCUL DU FLUX MASSIQUE VERTICAL
c Flux en l=1 (sol) nul
      fluxw=0.        
      do l=1,llm
           do j=1,jjp1
              do i=1,iip1              
               fluxw(i,j,l)=w(i,j,l)*g*0.01/aire(i,j)
C               print*,i,j,l,'fluxw(i,j,l)=',fluxw(i,j,l),
C     c                      'w(i,j,l)=',w(i,j,l)
              enddo
           enddo
      enddo
      
c INVERSION DES NIVEAUX
c le programme ppm3d travaille avec une 3ème coordonnée inversée par rapport
c de celle du LMDZ: z=1<=>niveau max, z=llm+1<=>surface
c On passe donc des niveaux du LMDZ à ceux de Lin
     
      do l=1,llm+1
          apppm(l)=ap(llm+2-l)
          bpppm(l)=bp(llm+2-l)         
      enddo 
     
      do l=1,llm
          do j=1,jjp1
             do i=1,iim     
                 unatppm(i,j,l)=unat(i,j,llm-l+1)
                 vnatppm(i,j,l)=vnat(i,j,llm-l+1)
                 fluxwppm(i,j,l)=fluxw(i,j,llm-l+1)
                 qppm(i,j,l)=q(i,j,llm-l+1)                              
             enddo
          enddo                                
      enddo
   
      return
      end






