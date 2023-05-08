










!
! $Id: friction.F 1437 2010-09-30 08:29:10Z emillour $
!
c=======================================================================
      SUBROUTINE friction(ucov,vcov,pdt)

      USE control_mod
! if not using IOIPSL, we still need to use (a local version of) getin
      USE ioipsl_getincom
      USE comconst_mod, ONLY: pi
      
      IMPLICIT NONE

!=======================================================================
!
!   Friction for the Newtonian case:
!   --------------------------------
!    2 possibilities (depending on flag 'friction_type'
!     friction_type=0 : A friction that is only applied to the lowermost
!                       atmospheric layer
!     friction_type=1 : Friction applied on all atmospheric layer (but
!       (default)       with stronger magnitude near the surface; see
!                       iniacademic.F)
!=======================================================================

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
!
! $Id: academic.h 1437 2010-09-30 08:29:10Z emillour $
!
      common/academic/tetarappel,knewt_t,kfrict,knewt_g,clat4
      real :: tetarappel(ip1jmp1,llm)
      real :: knewt_t(llm)
      real :: kfrict(llm)
      real :: knewt_g
      real :: clat4(ip1jmp1)

! arguments:
      REAL,INTENT(out) :: ucov( iip1,jjp1,llm )
      REAL,INTENT(out) :: vcov( iip1,jjm,llm )
      REAL,INTENT(in) :: pdt ! time step

! local variables:

      REAL modv(iip1,jjp1),zco,zsi
      REAL vpn,vps,upoln,upols,vpols,vpoln
      REAL u2(iip1,jjp1),v2(iip1,jjm)
      INTEGER  i,j,l
      REAL,PARAMETER :: cfric=1.e-5
      LOGICAL,SAVE :: firstcall=.true.
      INTEGER,SAVE :: friction_type=1
      CHARACTER(len=20) :: modname="friction"
      CHARACTER(len=80) :: abort_message
      
      IF (firstcall) THEN
        ! set friction type
        call getin("friction_type",friction_type)
        if ((friction_type.lt.0).or.(friction_type.gt.1)) then
          abort_message="wrong friction type"
          write(lunout,*)'Friction: wrong friction type',friction_type
          call abort_gcm(modname,abort_message,42)
        endif
        firstcall=.false.
      ENDIF

      if (friction_type.eq.0) then
c   calcul des composantes au carre du vent naturel
      do j=1,jjp1
         do i=1,iip1
            u2(i,j)=ucov(i,j,1)*ucov(i,j,1)*unscu2(i,j)
         enddo
      enddo
      do j=1,jjm
         do i=1,iip1
            v2(i,j)=vcov(i,j,1)*vcov(i,j,1)*unscv2(i,j)
         enddo
      enddo

c   calcul du module de V en dehors des poles
      do j=2,jjm
         do i=2,iip1
            modv(i,j)=sqrt(0.5*(u2(i-1,j)+u2(i,j)+v2(i,j-1)+v2(i,j)))
         enddo
         modv(1,j)=modv(iip1,j)
      enddo

c   les deux composantes du vent au pole sont obtenues comme
c   premiers modes de fourier de v pres du pole
      upoln=0.
      vpoln=0.
      upols=0.
      vpols=0.
      do i=2,iip1
         zco=cos(rlonv(i))*(rlonu(i)-rlonu(i-1))
         zsi=sin(rlonv(i))*(rlonu(i)-rlonu(i-1))
         vpn=vcov(i,1,1)/cv(i,1)
         vps=vcov(i,jjm,1)/cv(i,jjm)
         upoln=upoln+zco*vpn
         vpoln=vpoln+zsi*vpn
         upols=upols+zco*vps
         vpols=vpols+zsi*vps
      enddo
      vpn=sqrt(upoln*upoln+vpoln*vpoln)/pi
      vps=sqrt(upols*upols+vpols*vpols)/pi
      do i=1,iip1
c        modv(i,1)=vpn
c        modv(i,jjp1)=vps
         modv(i,1)=modv(i,2)
         modv(i,jjp1)=modv(i,jjm)
      enddo

c   calcul du frottement au sol.
      do j=2,jjm
         do i=1,iim
            ucov(i,j,1)=ucov(i,j,1)
     s      -cfric*pdt*0.5*(modv(i+1,j)+modv(i,j))*ucov(i,j,1)
         enddo
         ucov(iip1,j,1)=ucov(1,j,1)
      enddo
      do j=1,jjm
         do i=1,iip1
            vcov(i,j,1)=vcov(i,j,1)
     s      -cfric*pdt*0.5*(modv(i,j+1)+modv(i,j))*vcov(i,j,1)
         enddo
         vcov(iip1,j,1)=vcov(1,j,1)
      enddo
      endif ! of if (friction_type.eq.0)

      if (friction_type.eq.1) then
        do l=1,llm
          ucov(:,:,l)=ucov(:,:,l)*(1.-pdt*kfrict(l))
          vcov(:,:,l)=vcov(:,:,l)*(1.-pdt*kfrict(l))
        enddo
      endif
      
      RETURN
      END

