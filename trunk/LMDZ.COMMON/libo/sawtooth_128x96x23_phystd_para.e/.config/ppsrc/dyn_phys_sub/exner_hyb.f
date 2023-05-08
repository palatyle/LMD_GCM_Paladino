










      SUBROUTINE  exner_hyb ( ngrid, ps, p,beta, pks, pk, pkf )
c
c     Auteurs :  F. Forget , Y. Wanherdrick
c P.Le Van  , Fr. Hourdin  .
c    ..........
c
c    ....  ngrid, ps,p             sont des argum.d'entree  au sous-prog ...
c    .... beta, pks,pk,pkf   sont des argum.de sortie au sous-prog ...
c
c   ************************************************************************
c    Calcule la fonction d'Exner pk = Cp * (p/preff) ** kappa , aux milieux des 
c    couches .   Pk(l) sera calcule aux milieux  des couches l ,entre les
c    pressions p(l) et p(l+1) ,definis aux interfaces des llm couches .
c   ************************************************************************
c    .. N.B : Au sommet de l'atmosphere,  p(llm+1) = 0. , et ps et pks sont
c    la pression et la fonction d'Exner  au  sol  .
c
c     WARNING : CECI est une version speciale de exner_hyb originale
c               Utilis‰ dans la version martienne pour pouvoir 
c               tourner avec des coordonn‰es verticales complexe
c              => Il ne verifie PAS la condition la proportionalit‰ en 
c              ‰nergie totale/ interne / potentielle (F.Forget 2001)
c    ( voir note de Fr.Hourdin )  ,
c
      USE comvert_mod, ONLY: preff
      USE comconst_mod, ONLY: jmp1,kappa,cpp

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

      INTEGER  ngrid
      REAL p(ngrid,llmp1),pk(ngrid,llm),pkf(ngrid,llm)
      REAL ps(ngrid),pks(ngrid), beta(ngrid,llm)

c    .... variables locales   ...

      INTEGER l, ij
      REAL dum1

      REAL ppn(iim),pps(iim)
      REAL xpn, xps
      REAL SSUM
      EXTERNAL filtreg, SSUM
      
c     -------------
c     Calcul de pks
c     -------------
   
      DO   ij  = 1, ngrid
        pks(ij) = cpp * ( ps(ij)/preff ) ** kappa
      ENDDO

      DO  ij   = 1, iim
        ppn(ij) = aire(   ij   ) * pks(  ij     )
        pps(ij) = aire(ij+ip1jm) * pks(ij+ip1jm )
      ENDDO
      xpn      = SSUM(iim,ppn,1) /apoln
      xps      = SSUM(iim,pps,1) /apols

      DO ij   = 1, iip1
        pks(   ij     )  =  xpn
        pks( ij+ip1jm )  =  xps
      ENDDO
c
c
c    .... Calcul de pk  pour la couche l 
c    --------------------------------------------
c
      dum1 = cpp * (2*preff)**(-kappa) 
      DO l = 1, llm-1
        DO   ij   = 1, ngrid
         pk(ij,l) = dum1 * (p(ij,l) + p(ij,l+1))**kappa
        ENDDO
      ENDDO

c    .... Calcul de pk  pour la couche l = llm ..
c    (on met la meme distance (en log pression)  entre Pk(llm)
c    et Pk(llm -1) qu'entre Pk(llm-1) et Pk(llm-2)

      DO   ij   = 1, ngrid
         pk(ij,llm) = pk(ij,llm-1)**2 / pk(ij,llm-2)
      ENDDO


c    calcul de pkf
c    -------------
      CALL SCOPY   ( ngrid * llm, pk, 1, pkf, 1 )
      CALL filtreg ( pkf, jmp1, llm, 2, 1, .TRUE., 1 )
      
c    EST-CE UTILE ?? : calcul de beta
c    --------------------------------
      DO l = 2, llm
        DO   ij   = 1, ngrid
          beta(ij,l) = pk(ij,l) / pk(ij,l-1)   
        ENDDO
      ENDDO

      RETURN
      END
