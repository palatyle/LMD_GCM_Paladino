










      SUBROUTINE tourpot_p ( vcov, ucov, massebxy, vorpot )
      USE parallel_lmdz
      IMPLICIT NONE

c=======================================================================
c
c   Auteur:  P. Le Van
c   -------
c
c   Objet:
c   ------
c
c    *******************************************************************
c    .........      calcul du tourbillon potentiel             .........
c    *******************************************************************
c
c     vcov,ucov,fext et pbarxyfl sont des argum. d'entree pour le s-pg .
c             vorpot            est  un argum.de sortie pour le s-pg .
c
c=======================================================================

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

      REAL  rot( ip1jm,llm )
      REAL vcov( ip1jm,llm ),ucov( ip1jmp1,llm )
      REAL massebxy( ip1jm,llm ),vorpot( ip1jm,llm )

      INTEGER l, ij ,ije,ijb,jje,jjb


      ijb=ij_begin-iip1
      ije=ij_end
      
      if (pole_nord) ijb=ij_begin
      
      
c  ... vorpot = ( Filtre( d(vcov)/dx - d(ucov)/dy ) + fext ) /psbarxy ..



c    ........  Calcul du rotationnel du vent V  puis filtrage  ........
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO 5 l = 1,llm

      if (pole_sud)  ije=ij_end-iip1-1
      DO 2 ij = ijb, ije 
      rot( ij,l ) = vcov(ij+1,l)-vcov(ij,l)+ucov(ij+iip1,l)-ucov(ij,l)
   2  CONTINUE

c    ....  correction pour  rot( iip1,j,l )  .....
c    ....     rot(iip1,j,l) = rot(1,j,l)    .....

CDIR$ IVDEP

      if (pole_sud)  ije=ij_end-iip1
     
      DO 3 ij = ijb+iip1-1, ije, iip1
      rot( ij,l ) = rot( ij -iim, l )
   3  CONTINUE

   5  CONTINUE
c$OMP END DO NOWAIT
      jjb=jj_begin-1
      jje=jj_end
      
      if (pole_nord) jjb=jjb+1
      if (pole_sud)  jje=jje-1
      CALL  filtreg_p( rot, jjb,jje,jjm, llm, 2, 1, .FALSE., 1 )

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)     
      DO 10 l = 1, llm
      
      if (pole_sud)  ije=ij_end-iip1-1  
      
      DO 6 ij = ijb, ije
      vorpot( ij,l ) = ( rot(ij,l) + fext(ij) ) / massebxy(ij,l)
   6  CONTINUE

c    ..... correction pour  vorpot( iip1,j,l)  .....
c    ....   vorpot(iip1,j,l)= vorpot(1,j,l) ....
CDIR$ IVDEP
      if (pole_sud)  ije=ij_end-iip1
      DO 8 ij = ijb+iip1-1, ije, iip1
      vorpot( ij,l ) = vorpot( ij -iim,l )
   8  CONTINUE

  10  CONTINUE
c$OMP END DO NOWAIT
      RETURN
      END
