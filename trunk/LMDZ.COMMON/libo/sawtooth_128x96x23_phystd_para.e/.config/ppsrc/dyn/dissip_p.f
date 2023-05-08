










!
! $Id: $
!
      SUBROUTINE dissip_p( vcov,ucov,teta,p, dv,du,dh )
c
      USE parallel_lmdz
      USE write_field_p
      USE comconst_mod, ONLY: dtdiss
      IMPLICIT NONE


c ..  Avec nouveaux operateurs star :  gradiv2 , divgrad2, nxgraro2  ...
c                                 (  10/01/98  )

c=======================================================================
c
c   Auteur:  P. Le Van
c   -------
c
c   Objet:
c   ------
c
c   Dissipation horizontale
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

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
! $Id: comdissnew.h 1319 2010-02-23 21:29:54Z fairhead $
!
!
!  ATTENTION!!!!: ce fichier include est compatible format fixe/format libre
!                 veillez � n'utiliser que des ! pour les commentaires
!                 et � bien positionner les & des lignes de continuation 
!                 (les placer en colonne 6 et en colonne 73)
!
!-----------------------------------------------------------------------
! INCLUDE 'comdissnew.h'

      COMMON/comdissnew/ lstardis,nitergdiv,nitergrot,niterh,tetagdiv,  &
     &                   tetagrot,tetatemp,coefdis, vert_prof_dissip

      LOGICAL lstardis
      INTEGER nitergdiv, nitergrot, niterh

! For the Earth model:
      integer vert_prof_dissip ! vertical profile of horizontal dissipation
!     Allowed values:
!     0: rational fraction, function of pressure
!     1: tanh of altitude

      REAL     tetagdiv, tetagrot,  tetatemp, coefdis

!
! ... Les parametres de ce common comdissnew sont  lues par defrun_new 
!              sur le fichier  run.def    ....
!
!-----------------------------------------------------------------------
!
! $Header$
!
!  Attention : ce fichier include est compatible format fixe/format libre
!                 veillez à n'utiliser que des ! pour les commentaires
!                 et à bien positionner les & des lignes de continuation 
!                 (les placer en colonne 6 et en colonne 73)
!-----------------------------------------------------------------------
! INCLUDE comdissipn.h

      REAL  tetaudiv, tetaurot, tetah, cdivu, crot, cdivh
!
      COMMON/comdissipn/ tetaudiv(llm),tetaurot(llm),tetah(llm)   ,     &
     &                        cdivu,      crot,         cdivh

!
!    Les parametres de ce common proviennent des calculs effectues dans 
!             Inidissip  .
!
!-----------------------------------------------------------------------

c   Arguments:
c   ----------

      REAL,INTENT(IN) :: vcov(ip1jm,llm) ! covariant meridional wind
      REAL,INTENT(IN) :: ucov(ip1jmp1,llm) ! covariant zonal wind
      REAL,INTENT(IN) :: teta(ip1jmp1,llm) ! potentail temperature
      REAL,INTENT(IN) :: p(ip1jmp1,llmp1) ! pressure
      ! tendencies (.../s) on covariant winds and potential temperature
      REAL,INTENT(OUT) :: dv(ip1jm,llm)
      REAL,INTENT(OUT) :: du(ip1jmp1,llm)
      REAL,INTENT(OUT) :: dh(ip1jmp1,llm)

c   Local:
c   ------

      REAL gdx(ip1jmp1,llm),gdy(ip1jm,llm)
      REAL grx(ip1jmp1,llm),gry(ip1jm,llm)
      REAL te1dt(llm),te2dt(llm),te3dt(llm)
      REAL deltapres(ip1jmp1,llm)

      INTEGER l,ij

      REAL  SSUM
      integer :: ijb,ije
c-----------------------------------------------------------------------
c   initialisations:
c   ----------------

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
      DO l=1,llm
         te1dt(l) = tetaudiv(l) * dtdiss
         te2dt(l) = tetaurot(l) * dtdiss
         te3dt(l) = tetah(l)    * dtdiss
      ENDDO
c$OMP END DO NOWAIT
c      CALL initial0( ijp1llm, du )
c      CALL initial0( ijmllm , dv )
c      CALL initial0( ijp1llm, dh )
     
      ijb=ij_begin
      ije=ij_end

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
      DO l=1,llm
        du(ijb:ije,l)=0
        dh(ijb:ije,l)=0
      ENDDO
c$OMP END DO NOWAIT
      
      if (pole_sud) ije=ij_end-iip1

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
      DO l=1,llm
        dv(ijb:ije,l)=0
      ENDDO
c$OMP END DO NOWAIT
     
c-----------------------------------------------------------------------
c   Calcul de la dissipation:
c   -------------------------

c   Calcul de la partie   grad  ( div ) :
c   -------------------------------------
      
     
      
      IF(lstardis) THEN
c      IF (.FALSE.) THEN
         CALL gradiv2_p( llm,ucov,vcov,nitergdiv,gdx,gdy )
      ELSE
         CALL gradiv_p ( llm,ucov,vcov,nitergdiv,gdx,gdy )
      ENDIF


      ijb=ij_begin
      ije=ij_end
      if (pole_sud) ije=ij_end-iip1

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
      DO l=1,llm
         if (pole_nord) then
           DO ij = 1, iip1
              gdx(     ij ,l) = 0.
           ENDDO
         endif
         
         if (pole_sud) then
           DO ij = 1, iip1
              gdx(ij+ip1jm,l) = 0.
           ENDDO
         endif
         
         if (pole_nord) ijb=ij_begin+iip1
         DO ij = ijb,ije
            du(ij,l) = du(ij,l) - te1dt(l) *gdx(ij,l)
         ENDDO

         if (pole_nord) ijb=ij_begin
         DO ij = ijb,ije
            dv(ij,l) = dv(ij,l) - te1dt(l) *gdy(ij,l)
         ENDDO

       ENDDO
c$OMP END DO NOWAIT
c   calcul de la partie   n X grad ( rot ):
c   ---------------------------------------

      IF(lstardis) THEN
c      IF (.FALSE.) THEN
         CALL nxgraro2_p( llm,ucov, vcov, nitergrot,grx,gry )
      ELSE
         CALL nxgrarot_p( llm,ucov, vcov, nitergrot,grx,gry )
      ENDIF



      ijb=ij_begin
      ije=ij_end
      if (pole_sud) ije=ij_end-iip1

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)       
      DO l=1,llm
         
         if (pole_nord) then
           DO ij = 1, iip1
              grx(ij,l) = 0.
           ENDDO
         endif
         
         if (pole_nord) ijb=ij_begin+iip1
         DO ij = ijb,ije
            du(ij,l) = du(ij,l) - te2dt(l) * grx(ij,l)
         ENDDO
         
         if (pole_nord) ijb=ij_begin
         DO ij =  ijb, ije
            dv(ij,l) = dv(ij,l) - te2dt(l) * gry(ij,l)
         ENDDO
      
      ENDDO
c$OMP END DO NOWAIT

c   calcul de la partie   div ( grad ):
c   -----------------------------------

        
      IF(lstardis) THEN
c      IF (.FALSE.) THEN
    
      ijb=ij_begin
      ije=ij_end

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
       DO l = 1, llm
          DO ij = ijb, ije
            deltapres(ij,l) = AMAX1( 0.,  p(ij,l) - p(ij,l+1) )
          ENDDO
       ENDDO
c$OMP END DO NOWAIT
         CALL divgrad2_p( llm,teta, deltapres  ,niterh, gdx )
      ELSE
         CALL divgrad_p ( llm,teta, niterh, gdx        )
      ENDIF

c      call write_field3d_p('gdx2',reshape(gdx,(/iip1,jmp1,llm/)))
c      stop

      ijb=ij_begin
      ije=ij_end
      
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      DO l = 1,llm
         DO ij = ijb,ije
            dh( ij,l ) = dh( ij,l ) - te3dt(l) * gdx( ij,l )
         ENDDO
      ENDDO
c$OMP END DO NOWAIT

      RETURN
      END
