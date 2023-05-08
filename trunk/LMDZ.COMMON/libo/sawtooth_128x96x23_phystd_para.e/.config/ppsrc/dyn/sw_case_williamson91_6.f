










!
! $Id $
!
      SUBROUTINE sw_case_williamson91_6(vcov,ucov,teta,masse,ps)

c=======================================================================
c
c   Author:    Thomas Dubos      original: 26/01/2010
c   -------
c
c   Subject:
c   ------
c   Realise le cas-test 6 de Williamson et al. (1991) : onde de Rossby-Haurwitz
c
c   Method:
c   --------
c
c   Interface:
c   ----------
c
c      Input:
c      ------
c
c      Output:
c      -------
c
c=======================================================================
      USE comvert_mod, ONLY: ap,bp,preff
      USE comconst_mod, ONLY: cpp,omeg,rad

      IMPLICIT NONE
c-----------------------------------------------------------------------
c   Declararations:
c   ---------------

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
!
! $Header$
!
!
! gestion des impressions de sorties et de débogage
! lunout:    unité du fichier dans lequel se font les sorties 
!                           (par defaut 6, la sortie standard)
! prt_level: niveau d'impression souhaité (0 = minimum)
!
      INTEGER lunout, prt_level
      COMMON /comprint/ lunout, prt_level

c   Arguments:
c   ----------

c   variables dynamiques
      REAL vcov(ip1jm,llm),ucov(ip1jmp1,llm) ! vents covariants
      REAL teta(ip1jmp1,llm)                 ! temperature potentielle
      REAL ps(ip1jmp1)                       ! pression  au sol
      REAL masse(ip1jmp1,llm)                ! masse d'air
      REAL phis(ip1jmp1)                     ! geopotentiel au sol

c   Local:
c   ------

      REAL p (ip1jmp1,llmp1  )               ! pression aux interfac.des couches
      REAL pks(ip1jmp1)                      ! exner au  sol
      REAL pk(ip1jmp1,llm)                   ! exner au milieu des couches
      REAL pkf(ip1jmp1,llm)                  ! exner filt.au milieu des couches
      REAL alpha(ip1jmp1,llm),beta(ip1jmp1,llm)

      REAL :: sinth,costh,costh2, Ath,Bth,Cth, lon,dps
      INTEGER i,j,ij

      REAL, PARAMETER    :: rho=1 ! masse volumique de l'air (arbitraire)
      REAL, PARAMETER    :: K    = 7.848e-6  ! K = \omega
      REAL, PARAMETER    :: gh0  = 9.80616 * 8e3 
      INTEGER, PARAMETER :: R0=4, R1=R0+1, R2=R0+2         ! mode 4
c NB : rad = 6371220 dans W91 (6371229 dans LMDZ)
c      omeg = 7.292e-5 dans W91 (7.2722e-5 dans LMDZ)
 
      IF(0==0) THEN
c Williamson et al. (1991) : onde de Rossby-Haurwitz
         teta = preff/rho/cpp
c geopotentiel (pression de surface)
         do j=1,jjp1
            costh2 = cos(rlatu(j))**2
            Ath = (R0+1)*(costh2**2) + (2*R0*R0-R0-2)*costh2 - 2*R0*R0
            Ath = .25*(K**2)*(costh2**(R0-1))*Ath
            Ath = .5*K*(2*omeg+K)*costh2 + Ath 
            Bth = (R1*R1+1)-R1*R1*costh2
            Bth = 2*(omeg+K)*K/(R1*R2) * (costh2**(R0/2))*Bth
            Cth = R1*costh2 - R2
            Cth = .25*K*K*(costh2**R0)*Cth
            do i=1,iip1
               ij=(j-1)*iip1+i
               lon = rlonv(i)
               dps = Ath + Bth*cos(R0*lon) + Cth*cos(2*R0*lon)
               ps(ij) = rho*(gh0 + (rad**2)*dps)
            enddo
         enddo
         write(lunout,*) 'W91 ps', MAXVAL(ps), MINVAL(ps)
c vitesse zonale ucov
         do j=1,jjp1
            costh  = cos(rlatu(j))
            costh2 = costh**2
            Ath = rad*K*costh
            Bth = R0*(1-costh2)-costh2
            Bth = rad*K*Bth*(costh**(R0-1))
            do i=1,iip1
               ij=(j-1)*iip1+i
               lon = rlonu(i)
               ucov(ij,1) = (Ath + Bth*cos(R0*lon))
            enddo
         enddo
         write(lunout,*) 'W91 u', MAXVAL(ucov(:,1)), MINVAL(ucov(:,1))
         ucov(:,1)=ucov(:,1)*cu
c vitesse meridienne vcov
         do j=1,jjm
            sinth  = sin(rlatv(j))
            costh  = cos(rlatv(j))
            Ath = -rad*K*R0*sinth*(costh**(R0-1))
            do i=1,iip1
               ij=(j-1)*iip1+i
               lon = rlonv(i)
               vcov(ij,1) = Ath*sin(R0*lon)
            enddo
         enddo
         write(lunout,*) 'W91 v', MAXVAL(vcov(:,1)), MINVAL(vcov(:,1))
         vcov(:,1)=vcov(:,1)*cv
         
c         ucov=0
c         vcov=0
      ELSE
c test non-tournant, onde se propageant en latitude
         do j=1,jjp1
            do i=1,iip1
               ij=(j-1)*iip1+i
               ps(ij) = 1e5*(1 + .1*exp(-100*(1+sin(rlatu(j)))**2) )
            enddo
         enddo
         
c     rho = preff/(cpp*teta)
         teta = .01*preff/cpp   ! rho = 100 ; phi = ps/rho = 1e3 ; c=30 m/s = 2600 km/j = 23 degres / j
         ucov=0.
         vcov=0.
      END IF      
      
      CALL pression ( ip1jmp1, ap, bp, ps, p       )
      CALL massdair(p,masse)

      END
c-----------------------------------------------------------------------
