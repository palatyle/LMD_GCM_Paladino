










      subroutine groupe_p(pext,pbaru,pbarv,pbarum,pbarvm,wm)
      USE parallel_lmdz
      USE comconst_mod, ONLY: ngroup
      implicit none

c   sous-programme servant a fitlrer les champs de flux de masse aux
c   poles en "regroupant" les mailles 2 par 2 puis 4 par 4 etc. au fur
c   et a mesure qu'on se rapproche du pole.
c
c   en entree: pext, pbaru et pbarv
c
c   en sortie:  pbarum,pbarvm et wm.
c
c   remarque, le wm est recalcule a partir des pbaru pbarv et on n'a donc
c   pas besoin de w en entree.

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

!      integer ngroup
!      parameter (ngroup=3)


      real pbaru(iip1,jjp1,llm),pbarv(iip1,jjm,llm)
      real pext(iip1,jjp1,llm)

      real pbarum(iip1,jjp1,llm),pbarvm(iip1,jjm,llm)
      real wm(iip1,jjp1,llm)

      real,save :: zconvm(iip1,jjp1,llm)
      real,save :: zconvmm(iip1,jjp1,llm)

      real uu

      integer i,j,l

      logical firstcall,groupe_ok
      save firstcall,groupe_ok
c$OMP THREADPRIVATE(firstcall,groupe_ok)

      data firstcall/.true./
      data groupe_ok/.true./

      integer ijb,ije,jjb,jje
      
      if (iim==1) then
         groupe_ok=.false.
      endif

      if (firstcall) then
         if (groupe_ok) then
           if(mod(iim,2**ngroup).ne.0) stop'probleme du nombre de point'
         endif
         firstcall=.false.
      endif

c   Champs 1D

      call convflu_p(pbaru,pbarv,llm,zconvm)

c
c      call scopy(ijp1llm,zconvm,1,zconvmm,1)
c      call scopy(ijmllm,pbarv,1,pbarvm,1)
      
      jjb=jj_begin
      jje=jj_end

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      do l=1,llm
        zconvmm(:,jjb:jje,l)=zconvm(:,jjb:jje,l)
      enddo
c$OMP END DO NOWAIT

      if (groupe_ok) then
         call groupeun_p(jjp1,llm,jjb,jje,zconvmm)
      endif
      
      jjb=jj_begin-1
      jje=jj_end
      if (pole_nord) jjb=jj_begin
      if (pole_sud)  jje=jj_end-1
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      do l=1,llm
        pbarvm(:,jjb:jje,l)=pbarv(:,jjb:jje,l)
      enddo
c$OMP END DO NOWAIT

      if (groupe_ok) then
         call groupeun_p(jjm,llm,jjb,jje,pbarvm)
      endif

c   Champs 3D
   
      jjb=jj_begin
      jje=jj_end
      if (pole_nord) jjb=jj_begin+1
      if (pole_sud)  jje=jj_end-1
      
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      do l=1,llm
         do j=jjb,jje
            uu=pbaru(iim,j,l)
            do i=1,iim
               uu=uu+pbarvm(i,j,l)-pbarvm(i,j-1,l)-zconvmm(i,j,l)
               pbarum(i,j,l)=uu
c     zconvm(i,j,l ) =  xflu(i-1,j,l)-xflu(i,j,l)+
c    *                      yflu(i,j,l)-yflu(i,j-1,l)
            enddo
            pbarum(iip1,j,l)=pbarum(1,j,l)
         enddo
      enddo
c$OMP END DO NOWAIT

c    integration de la convergence de masse de haut  en bas ......
   
      jjb=jj_begin
      jje=jj_end

c$OMP BARRIER
c$OMP MASTER      
      do  l = llm-1,1,-1
          do j=jjb,jje
             do i=1,iip1
                zconvmm(i,j,l)=zconvmm(i,j,l)+zconvmm(i,j,l+1)
             enddo
          enddo
      enddo

      if (.not. pole_sud) then
        zconvmm(:,jj_end+1,:)=0
cym	wm(:,jj_end+1,:)=0
      endif
      
c$OMP END MASTER
c$OMP BARRIER      

      CALL vitvert_p(zconvmm(1,1,1),wm(1,1,1))

      return
      end

