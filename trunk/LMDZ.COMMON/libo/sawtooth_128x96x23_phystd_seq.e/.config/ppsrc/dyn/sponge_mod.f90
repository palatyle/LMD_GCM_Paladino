










module sponge_mod

implicit none

! sponge parameters (set/read via conf_gcm.F)
logical,save :: callsponge  ! do we use a sponge on upper layers
integer,save :: mode_sponge ! sponge mode
integer,save :: nsponge ! number of sponge layers 
real,save :: tetasponge  ! sponge time scale (s) at topmost layer


contains

      subroutine sponge(ucov,vcov,h,ps,dt,mode)

! Sponge routine: Quench ucov, vcov and potential temperature near the
!                 top of the model
! Depending on 'mode' relaxation of variables is towards:
! mode = 0 : h -> h_mean , ucov -> 0 , vcov -> 0
! mode = 1 : h -> h_mean , ucov -> ucov_mean , vcov -> 0
! mode >= 2 : h -> h_mean , ucov -> ucov_mean , vcov -> vcov_mean
! Number of layer over which sponge is applied is 'nsponge' (read from def file)
! Time scale for quenching at top level is given by 'tetasponge' (read from
! def file) and doubles as level indexes decrease.
! Quenching is modeled as: A(t)=Am+A0exp(-lambda*t)
! where Am is the zonal average of the field (or zero), and lambda the inverse
! of the characteristic quenching/relaxation time scale
! Thus, assuming Am to be time-independent, field at time t+dt is given by:
! A(t+dt)=A(t)-(A(t)-Am)*(1-exp(-lambda*dt))

      USE comvert_mod, ONLY: ap,bp,preff,scaleheight

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

! Arguments:
!------------
      real,intent(inout) :: ucov(iip1,jjp1,llm) ! covariant zonal wind
      real,intent(inout) :: vcov(iip1,jjm,llm) ! covariant meridional wind
      real,intent(inout) :: h(iip1,jjp1,llm) ! potential temperature
!      real,intent(in) :: pext(iip1,jjp1) ! extensive pressure
      real,intent(in) :: ps(iip1,jjp1) ! surface pressure (Pa)
      real,intent(in) :: dt   ! time step
      integer,intent(in) :: mode  ! sponge mode

!   Local:
!   ------

      real,save :: sig_s(llm) !sigma au milieu des couches
      REAL vm,um,hm,ptot(jjp1)
      real,save :: cst(llm)
      real :: pext(iip1,jjp1) ! extensive pressure

      INTEGER l,i,j
      integer,save :: l0 ! layer down to which sponge is applied

      real ssum

      real zkm
      logical,save :: firstcall=.true.



      if (firstcall) then

       ! build approximative sigma levels at midlayer
        do l=1,llm
          sig_s(l)=((ap(l)+ap(l+1))/preff+bp(l)+bp(l+1))/2.
        enddo

        l0=llm-nsponge+1
 
        write(lunout,*)
        write(lunout,*)'sponge mode',mode
        write(lunout,*)'nsponge tetasponge ',nsponge,tetasponge
        write(lunout,*)'Coeffs for the sponge layer'
        write(lunout,*)'Z (km)         tau             cst'
        do l=llm,l0,-1
          ! double time scale with every level, starting from the top
          cst(l)=1.-exp(-dt/(tetasponge*2**(llm-l)))
        enddo

        do l=l0,llm
           zkm=-scaleheight*log(sig_s(l))
           print*,zkm,tetasponge*2**(llm-l),cst(l)
        enddo
        PRINT*

        firstcall=.false.
      endif ! of if (firstcall)

!-----------------------------------------------------------------------
!   compute sponge relaxation:
!   -------------------------

      pext(1:iip1,1:jjp1)=ps(1:iip1,1:jjp1)*aire(1:iip1,1:jjp1)

      do j=1,jjp1
         ptot(j)=sum(pext(1:iim,j))
      enddo
 
!   potential temperature
      do l=l0,llm
         do j=1,jjp1
            ! compute zonal average
            hm=0.
            do i=1,iim
               hm=hm+h(i,j,l)*pext(i,j)
            enddo
            hm=hm/ptot(j)
            ! update h()
            do i=1,iim
               h(i,j,l)=h(i,j,l)-cst(l)*(h(i,j,l)-hm)
            enddo
            h(iip1,j,l)=h(1,j,l)
         enddo
      enddo

!   zonal wind
      do l=l0,llm
         do j=2,jjm
            um=0.
            if(mode.ge.1) then ! compute zonal average
               do i=1,iim
                  um=um+0.5*ucov(i,j,l)*(pext(i,j)+pext(i+1,j))/cu(i,j)
               enddo
               um=um/ptot(j)
            endif
            ! update ucov()
            do i=1,iim
               ucov(i,j,l)=ucov(i,j,l)-cst(l)*(ucov(i,j,l)-um*cu(i,j))
            enddo
            ucov(iip1,j,l)=ucov(1,j,l)
         enddo
      enddo

!  meridional wind
      do l=l0,llm
         do j=1,jjm
            vm=0.
            if(mode.ge.2) then ! compute zonal average
               do i=1,iim
                  vm=vm+vcov(i,j,l)*(pext(i,j)+pext(i,j+1))/cv(i,j)
               enddo
               vm=vm/(ptot(j)+ptot(j+1))
            endif
            ! update vcov
            do i=1,iim
               vcov(i,j,l)=vcov(i,j,l)-cst(l)*(vcov(i,j,l)-vm*cv(i,j))
            enddo
            vcov(iip1,j,l)=vcov(1,j,l)
         enddo
      enddo

      end subroutine sponge
      
end module sponge_mod

