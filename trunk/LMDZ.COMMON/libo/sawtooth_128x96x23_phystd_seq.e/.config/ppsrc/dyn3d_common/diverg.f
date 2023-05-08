










!
! $Header$
!
      SUBROUTINE diverg(klevel,x,y,div)
c
c     P. Le Van
c
c  *********************************************************************
c  ... calcule la divergence a tous les niveaux d'1 vecteur de compos. 
c     x et y...
c              x et y  etant des composantes covariantes   ...
c  *********************************************************************
      IMPLICIT NONE
c
c      x  et  y  sont des arguments  d'entree pour le s-prog
c        div      est  un argument  de sortie pour le s-prog
c
c
c   ---------------------------------------------------------------------
c
c    ATTENTION : pendant ce s-pg , ne pas toucher au COMMON/scratch/  .
c
c   ---------------------------------------------------------------------
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
c    ..........          variables en arguments    ...................
c
      INTEGER klevel
      REAL x( ip1jmp1,klevel ),y( ip1jm,klevel ),div( ip1jmp1,klevel )
      INTEGER   l,ij
c
c    ...............     variables  locales   .........................

      REAL aiy1( iip1 ) , aiy2( iip1 )
      REAL sumypn,sumyps
c    ...................................................................
c
      REAL      SSUM
c
c
      DO 10 l = 1,klevel
c
        DO  ij = iip2, ip1jm - 1
         div( ij + 1, l )     =  
     *   cvusurcu( ij+1 ) * x( ij+1,l ) - cvusurcu( ij ) * x( ij , l) +
     *   cuvsurcv(ij-iim) * y(ij-iim,l) - cuvsurcv(ij+1) * y(ij+1,l) 
        ENDDO
c
c     ....  correction pour  div( 1,j,l)  ......
c     ....   div(1,j,l)= div(iip1,j,l) ....
c
CDIR$ IVDEP
        DO  ij = iip2,ip1jm,iip1
         div( ij,l ) = div( ij + iim,l )
        ENDDO
c
c     ....  calcul  aux poles  .....
c
        DO  ij  = 1,iim
         aiy1(ij) =    cuvsurcv(    ij       ) * y(     ij     , l )
         aiy2(ij) =    cuvsurcv( ij+ ip1jmi1 ) * y( ij+ ip1jmi1, l )
        ENDDO
        sumypn = SSUM ( iim,aiy1,1 ) / apoln
        sumyps = SSUM ( iim,aiy2,1 ) / apols
c
        DO  ij = 1,iip1
         div(     ij    , l ) = - sumypn
         div( ij + ip1jm, l ) =   sumyps
        ENDDO
  10  CONTINUE
c

ccc        CALL filtreg( div, jjp1, klevel, 2, 2, .TRUE., 1 )
      
c
        DO l = 1, klevel
           DO ij = iip2,ip1jm
            div(ij,l) = div(ij,l) * unsaire(ij) 
          ENDDO
        ENDDO
c
       RETURN
       END
