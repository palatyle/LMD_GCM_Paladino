










      SUBROUTINE divgrad_p (klevel,h, lh, divgra_out )
      USE parallel_lmdz
      USE times
      IMPLICIT NONE
c
c=======================================================================
c
c  Auteur :   P. Le Van
c  ----------
c
c                              lh
c      calcul de  (div( grad ))   de h  .....
c      h  et lh  sont des arguments  d'entree pour le s-prog
c      divgra     est  un argument  de sortie pour le s-prog
c
c=======================================================================
c
c   declarations:
c   -------------
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
c
      INTEGER klevel
      REAL h( ip1jmp1,klevel ), divgra_out( ip1jmp1,klevel )
      REAL,SAVE :: divgra( ip1jmp1,llm )

c
      REAL ghy(ip1jm,llm), ghx(ip1jmp1,llm)

      INTEGER  l,ij,iter,lh
c
      INTEGER ijb,ije,jjb,jje
c
c
c      CALL SCOPY ( ip1jmp1*klevel,h,1,divgra,1 )
      
      ijb=ij_begin
      ije=ij_end
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l = 1, klevel      
      divgra(ijb:ije,l)=h(ijb:ije,l)
      ENDDO
c$OMP END DO NOWAIT
c

c
      DO 10 iter = 1,lh

      jjb=jj_begin
      jje=jj_end
      CALL filtreg_p ( divgra,jjb,jje,jjp1,klevel,2,1,.true.,1  )

c      call exchange_Hallo(divgra,ip1jmp1,llm,0,1)
c$OMP BARRIER
c$OMP MASTER      
      call suspend_timer(timer_dissip)
      call exchange_Hallo(divgra,ip1jmp1,llm,1,1)
      call resume_timer(timer_dissip)
c$OMP END MASTER
c$OMP BARRIER       
      CALL    grad_p (klevel,divgra, ghx  , ghy          )

c$OMP BARRIER
c$OMP MASTER   
      call suspend_timer(timer_dissip)
      call exchange_Hallo(ghy,ip1jm,llm,1,0)
      call resume_timer(timer_dissip)
c$OMP END MASTER
c$OMP BARRIER            

      CALL  diverg_p (klevel,  ghx , ghy  , divgra       )

      jjb=jj_begin
      jje=jj_end
      CALL filtreg_p( divgra,jjb,jje,jjp1,klevel,2,1,.true.,1)

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO 5 l = 1,klevel
      DO 4  ij = ijb, ije
      divgra_out( ij,l ) = - cdivh * divgra( ij,l )
   4  CONTINUE
   5  CONTINUE
c$OMP END DO NOWAIT
c
  10  CONTINUE
      RETURN
      END
