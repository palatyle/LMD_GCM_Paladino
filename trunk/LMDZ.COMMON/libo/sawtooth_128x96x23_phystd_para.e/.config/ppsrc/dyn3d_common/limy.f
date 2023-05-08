










!
! $Header$
!
      SUBROUTINE limy(s0,sy,sm,pente_max)
c
c     Auteurs:   P.Le Van, F.Hourdin, F.Forget 
c
c    ********************************************************************
c     Shema  d'advection " pseudo amont " .
c    ********************************************************************
c     q,w sont des arguments d'entree  pour le s-pg ....
c     dq 	       sont des arguments de sortie pour le s-pg ....
c
c
c   --------------------------------------------------------------------
      USE comconst_mod, ONLY: pi

      IMPLICIT NONE
c
!-----------------------------------------------------------------------
!   INCLUDE 'dimensions.h'
!
!   dimensions.h contient les dimensions du modele
!   ndm est tel que iim=2**ndm
!-----------------------------------------------------------------------

      INTEGER iim,jjm,llm,ndm

      PARAMETER (iim= 128,jjm=96,llm=23,ndm=1)

!-----------------------------------------------------------------------
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
!CDK comgeom
      COMMON/comgeom/                                                   &
     & cu(ip1jmp1),cv(ip1jm),unscu2(ip1jmp1),unscv2(ip1jm),             &
     & aire(ip1jmp1),airesurg(ip1jmp1),aireu(ip1jmp1),                  &
     & airev(ip1jm),unsaire(ip1jmp1),apoln,apols,                       &
     & unsairez(ip1jm),airuscv2(ip1jm),airvscu2(ip1jm),                 &
     & aireij1(ip1jmp1),aireij2(ip1jmp1),aireij3(ip1jmp1),              &
     & aireij4(ip1jmp1),alpha1(ip1jmp1),alpha2(ip1jmp1),                &
     & alpha3(ip1jmp1),alpha4(ip1jmp1),alpha1p2(ip1jmp1),               &
     & alpha1p4(ip1jmp1),alpha2p3(ip1jmp1),alpha3p4(ip1jmp1),           &
     & fext(ip1jm),constang(ip1jmp1),rlatu(jjp1),rlatv(jjm),            &
     & rlonu(iip1),rlonv(iip1),cuvsurcv(ip1jm),cvsurcuv(ip1jm),         &
     & cvusurcu(ip1jmp1),cusurcvu(ip1jmp1),cuvscvgam1(ip1jm),           &
     & cuvscvgam2(ip1jm),cvuscugam1(ip1jmp1),                           &
     & cvuscugam2(ip1jmp1),cvscuvgam(ip1jm),cuscvugam(ip1jmp1),         &
     & unsapolnga1,unsapolnga2,unsapolsga1,unsapolsga2,                 &
     & unsair_gam1(ip1jmp1),unsair_gam2(ip1jmp1),unsairz_gam(ip1jm),    &
     & aivscu2gam(ip1jm),aiuscv2gam(ip1jm),xprimu(iip1),xprimv(iip1)

!
        REAL                                                            &
     & cu,cv,unscu2,unscv2,aire,airesurg,aireu,airev,unsaire,apoln     ,&
     & apols,unsairez,airuscv2,airvscu2,aireij1,aireij2,aireij3,aireij4,&
     & alpha1,alpha2,alpha3,alpha4,alpha1p2,alpha1p4,alpha2p3,alpha3p4 ,&
     & fext,constang,rlatu,rlatv,rlonu,rlonv,cuvscvgam1,cuvscvgam2     ,&
     & cvuscugam1,cvuscugam2,cvscuvgam,cuscvugam,unsapolnga1,unsapolnga2&
     & ,unsapolsga1,unsapolsga2,unsair_gam1,unsair_gam2,unsairz_gam    ,&
     & aivscu2gam ,aiuscv2gam,cuvsurcv,cvsurcuv,cvusurcu,cusurcvu,xprimu&
     & , xprimv
!
c
c
c   Arguments:
c   ----------
      real pente_max
      real s0(ip1jmp1,llm),sy(ip1jmp1,llm),sm(ip1jmp1,llm)
c
c      Local 
c   ---------
c
      INTEGER i,ij,l
c
      REAL q(ip1jmp1,llm)
      REAL airej2,airejjm,airescb(iim),airesch(iim)
      real sigv,dyq(ip1jmp1),dyqv(ip1jm)
      real adyqv(ip1jm),dyqmax(ip1jmp1)
      REAL qbyv(ip1jm,llm)

      REAL qpns,qpsn,appn,apps,dyn1,dys1,dyn2,dys2
      Logical extremum,first
      save first

      real convpn,convps,convmpn,convmps
      real sinlon(iip1),sinlondlon(iip1)
      real coslon(iip1),coslondlon(iip1)
      save sinlon,coslon,sinlondlon,coslondlon
c
c
      REAL      SSUM
      integer ismax,ismin
      EXTERNAL  SSUM, convflu,ismin,ismax
      EXTERNAL filtreg

      data first/.true./

      if(first) then
         print*,'SCHEMA AMONT NOUVEAU'
         first=.false.
         do i=2,iip1
            coslon(i)=cos(rlonv(i))
            sinlon(i)=sin(rlonv(i))
            coslondlon(i)=coslon(i)*(rlonu(i)-rlonu(i-1))/pi
            sinlondlon(i)=sinlon(i)*(rlonu(i)-rlonu(i-1))/pi
         enddo
         coslon(1)=coslon(iip1)
         coslondlon(1)=coslondlon(iip1)
         sinlon(1)=sinlon(iip1)
         sinlondlon(1)=sinlondlon(iip1)
      endif

c

      do l = 1, llm
c
         DO ij=1,ip1jmp1
               q(ij,l) = s0(ij,l) / sm ( ij,l )
               dyq(ij) = sy(ij,l) / sm ( ij,l )
         ENDDO
c
c   --------------------------------
c      CALCUL EN LATITUDE
c   --------------------------------

c   On commence par calculer la valeur du traceur moyenne sur le premier cercle
c   de latitude autour du pole (qpns pour le pole nord et qpsn pour
c    le pole nord) qui sera utilisee pour evaluer les pentes au pole.

      airej2 = SSUM( iim, aire(iip2), 1 )
      airejjm= SSUM( iim, aire(ip1jm -iim), 1 ) 
      DO i = 1, iim
      airescb(i) = aire(i+ iip1) * q(i+ iip1,l)
      airesch(i) = aire(i+ ip1jm- iip1) * q(i+ ip1jm- iip1,l)
      ENDDO
      qpns   = SSUM( iim,  airescb ,1 ) / airej2
      qpsn   = SSUM( iim,  airesch ,1 ) / airejjm

c   calcul des pentes aux points v

      do ij=1,ip1jm
         dyqv(ij)=q(ij,l)-q(ij+iip1,l)
         adyqv(ij)=abs(dyqv(ij))
      ENDDO

c   calcul des pentes aux points scalaires

      do ij=iip2,ip1jm
         dyqmax(ij)=min(adyqv(ij-iip1),adyqv(ij))
         dyqmax(ij)=pente_max*dyqmax(ij)
      enddo

c   calcul des pentes aux poles

c   calcul des pentes limites aux poles

c     print*,dyqv(iip1+1)
c     appn=abs(dyq(1)/dyqv(iip1+1))
c     print*,dyq(ip1jm+1)
c     print*,dyqv(ip1jm-iip1+1)
c     apps=abs(dyq(ip1jm+1)/dyqv(ip1jm-iip1+1))
c     do ij=2,iim
c        appn=amax1(abs(dyq(ij)/dyqv(ij)),appn)
c        apps=amax1(abs(dyq(ip1jm+ij)/dyqv(ip1jm-iip1+ij)),apps)
c     enddo
c     appn=min(pente_max/appn,1.)
c     apps=min(pente_max/apps,1.)


c   cas ou on a un extremum au pole

c     if(dyqv(ismin(iim,dyqv,1))*dyqv(ismax(iim,dyqv,1)).le.0.)
c    &   appn=0.
c     if(dyqv(ismax(iim,dyqv(ip1jm-iip1+1),1)+ip1jm-iip1+1)*
c    &   dyqv(ismin(iim,dyqv(ip1jm-iip1+1),1)+ip1jm-iip1+1).le.0.)
c    &   apps=0.

c   limitation des pentes aux poles
c     do ij=1,iip1
c        dyq(ij)=appn*dyq(ij)
c        dyq(ip1jm+ij)=apps*dyq(ip1jm+ij)
c     enddo

c   test
c      do ij=1,iip1
c         dyq(iip1+ij)=0.
c         dyq(ip1jm+ij-iip1)=0.
c      enddo
c      do ij=1,ip1jmp1
c         dyq(ij)=dyq(ij)*cos(rlatu((ij-1)/iip1+1))
c      enddo

      if(dyqv(ismin(iim,dyqv,1))*dyqv(ismax(iim,dyqv,1)).le.0.)
     &   then
         do ij=1,iip1
            dyqmax(ij)=0.
         enddo
      else
         do ij=1,iip1
            dyqmax(ij)=pente_max*abs(dyqv(ij))
         enddo
      endif

      if(dyqv(ismax(iim,dyqv(ip1jm-iip1+1),1)+ip1jm-iip1+1)*
     & dyqv(ismin(iim,dyqv(ip1jm-iip1+1),1)+ip1jm-iip1+1).le.0.)
     &then
         do ij=ip1jm+1,ip1jmp1
            dyqmax(ij)=0.
         enddo
      else
         do ij=ip1jm+1,ip1jmp1
            dyqmax(ij)=pente_max*abs(dyqv(ij-iip1))
         enddo
      endif

c   calcul des pentes limitees

      do ij=1,ip1jmp1
         if(dyqv(ij)*dyqv(ij-iip1).gt.0.) then
            dyq(ij)=sign(min(abs(dyq(ij)),dyqmax(ij)),dyq(ij))
         else
            dyq(ij)=0.
         endif
      enddo

         DO ij=1,ip1jmp1
               sy(ij,l) = dyq(ij) * sm ( ij,l )
        ENDDO

      enddo ! fin de la boucle sur les couches verticales

      RETURN
      END
