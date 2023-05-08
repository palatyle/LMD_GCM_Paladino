










!
! $Header$
!
      SUBROUTINE groupeun(jjmax,llmax,q)
      
      USE comconst_mod, ONLY: ngroup
      IMPLICIT NONE

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

      INTEGER jjmax,llmax
      REAL q(iip1,jjmax,llmax)

!      INTEGER ngroup
!      PARAMETER (ngroup=3)

      REAL airecn,qn
      REAL airecs,qs

      INTEGER i,j,l,ig,ig2,j1,j2,i0,jd

c--------------------------------------------------------------------c 
c Strategie d'optimisation                                           c
c stocker les valeurs systematiquement recalculees                   c
c et identiques d'un pas de temps sur l'autre. Il s'agit des         c
c aires des cellules qui sont sommees. S'il n'y a pas de changement  c
c de grille au cours de la simulation tout devrait bien se passer.   c
c Autre optimisation : determination des bornes entre lesquelles "j" c
c varie, au lieu de faire un test à chaque fois...
c--------------------------------------------------------------------c 

      INTEGER j_start, j_finish

      REAL, SAVE :: airen_tab(iip1,jjp1,0:1)
      REAL, SAVE :: aires_tab(iip1,jjp1,0:1)

      LOGICAL, SAVE :: first = .TRUE.
!      INTEGER,SAVE :: i_index(iim,ngroup)
      INTEGER      :: offset
!      REAL         :: qsum(iim/ngroup)

      IF (first) THEN
         CALL INIT_GROUPEUN(airen_tab, aires_tab)
         first = .FALSE.
      ENDIF


c Champs 3D
      jd=jjp1-jjmax
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
         j1=1+jd
         j2=2
         DO ig=1,ngroup

c     Concerne le pole nord
            j_start  = j1-jd
            j_finish = j2-jd
            DO ig2=1,ngroup-ig+1
              offset=2**(ig2-1)
              DO j=j_start, j_finish
!CDIR NODEP
!CDIR ON_ADB(q)
                 DO i0=1,iim,2**ig2
                   q(i0,j,l)=q(i0,j,l)+q(i0+offset,j,l) 
                 ENDDO
              ENDDO
            ENDDO
            
            DO j=j_start, j_finish
!CDIR NODEP
!CDIR ON_ADB(q)
               DO i=1,iim
                 q(i,j,l)=q(i-MOD(i-1,2**(ngroup-ig+1)),j,l)
               ENDDO
            ENDDO

            DO j=j_start, j_finish
!CDIR ON_ADB(airen_tab)
!CDIR ON_ADB(q)
               DO i=1,iim
                 q(i,j,l)=q(i,j,l)*airen_tab(i,j,jd)
               ENDDO
               q(iip1,j,l)=q(1,j,l)
            ENDDO
       
!c     Concerne le pole sud
            j_start  = j1-jd
            j_finish = j2-jd
            DO ig2=1,ngroup-ig+1
              offset=2**(ig2-1)
              DO j=j_start, j_finish
!CDIR NODEP
!CDIR ON_ADB(q)
                 DO i0=1,iim,2**ig2
                   q(i0,jjp1-j+1-jd,l)= q(i0,jjp1-j+1-jd,l)
     &                                 +q(i0+offset,jjp1-j+1-jd,l) 
                 ENDDO
              ENDDO
            ENDDO


            DO j=j_start, j_finish
!CDIR NODEP
!CDIR ON_ADB(q)
               DO i=1,iim
                 q(i,jjp1-j+1-jd,l)=q(i-MOD(i-1,2**(ngroup-ig+1)),
     &                                jjp1-j+1-jd,l)
               ENDDO
            ENDDO

            DO j=j_start, j_finish
!CDIR ON_ADB(aires_tab)
!CDIR ON_ADB(q)
               DO i=1,iim
                 q(i,jjp1-j+1-jd,l)=q(i,jjp1-j+1-jd,l)*  
     &                              aires_tab(i,jjp1-j+1,jd)
               ENDDO
               q(iip1,jjp1-j+1-jd,l)=q(1,jjp1-j+1-jd,l)
            ENDDO

        
            j1=j2+1
            j2=j2+2**ig
         ENDDO
      ENDDO
!$OMP END DO NOWAIT

      RETURN
      END
      
      
      
      
      SUBROUTINE INIT_GROUPEUN(airen_tab, aires_tab)

      USE comconst_mod, ONLY: ngroup
      IMPLICIT NONE

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

!      INTEGER ngroup
!      PARAMETER (ngroup=3)

      REAL airen,airecn
      REAL aires,airecs

      INTEGER i,j,l,ig,j1,j2,i0,jd

      INTEGER j_start, j_finish

      REAL :: airen_tab(iip1,jjp1,0:1)
      REAL :: aires_tab(iip1,jjp1,0:1)

      DO jd=0, 1
         j1=1+jd
         j2=2
         DO ig=1,ngroup
            
!     c     Concerne le pole nord
            j_start = j1-jd
            j_finish = j2-jd
            DO j=j_start, j_finish
               DO i0=1,iim,2**(ngroup-ig+1)
                  airen=0.
                  DO i=i0,i0+2**(ngroup-ig+1)-1
                     airen = airen+aire(i,j)
                  ENDDO
                  DO i=i0,i0+2**(ngroup-ig+1)-1
                     airen_tab(i,j,jd) = 
     &                    aire(i,j) / airen
                  ENDDO
               ENDDO
            ENDDO
            
!     c     Concerne le pole sud
            j_start = j1-jd
            j_finish = j2-jd
            DO j=j_start, j_finish
               DO i0=1,iim,2**(ngroup-ig+1)
                  aires=0.
                  DO i=i0,i0+2**(ngroup-ig+1)-1
                     aires=aires+aire(i,jjp1-j+1)
                  ENDDO
                  DO i=i0,i0+2**(ngroup-ig+1)-1
                     aires_tab(i,jjp1-j+1,jd) = 
     &                    aire(i,jjp1-j+1) / aires
                  ENDDO
               ENDDO
            ENDDO
            
            j1=j2+1
            j2=j2+2**ig
         ENDDO
      ENDDO
      
      RETURN
      END
