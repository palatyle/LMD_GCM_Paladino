










!
! $Header$
!
      SUBROUTINE pbar ( pext, pbarx, pbary, pbarxy )
      IMPLICIT NONE

c=======================================================================
c
c   Auteur:  P. Le Van
c   -------
c
c   Objet:
c   ------
c
c **********************************************************************
c calcul des moyennes en x et en y de (pression au sol*aire variable) ..
c *********************************************************************
c
c          pext               est  un argum. d'entree  pour le s-pg ..
c     pbarx,pbary et pbarxy  sont des argum. de sortie pour le s-pg ..
c
c   Methode:
c   --------
c
c    A chaque point scalaire P (i,j) est affecte 4 coefficients d'aires
c       alpha1(i,j)  calcule  au point ( i+1/4,j-1/4 )
c       alpha2(i,j)  calcule  au point ( i+1/4,j+1/4 )
c       alpha3(i,j)  calcule  au point ( i-1/4,j+1/4 )
c       alpha4(i,j)  calcule  au point ( i-1/4,j-1/4 )
c
c    Avec  alpha1(i,j) = aire(i+1/4,j-1/4)/ aire(i,j)        
c
c    N.B .  Pour plus de details, voir s-pg  ...  iniconst ...
c
c
c
c   alpha4 .         . alpha1    . alpha4
c    (i,j)             (i,j)       (i+1,j)
c
c             P .        U .          . P
c           (i,j)       (i,j)         (i+1,j)
c
c   alpha3 .         . alpha2    .alpha3 
c    (i,j)              (i,j)     (i+1,j)
c
c             V .        Z .          . V
c           (i,j)
c
c   alpha4 .         . alpha1    .alpha4
c   (i,j+1)            (i,j+1)   (i+1,j+1) 
c
c             P .        U .          . P
c          (i,j+1)                    (i+1,j+1)
c
c
c
c
c                       On  a :
c
c    pbarx(i,j) = Pext(i  ,j) * ( alpha1(i  ,j) + alpha2(i,j))      +
c                 Pext(i+1,j) * ( alpha3(i+1,j) + alpha4(i+1,j) )
c     localise  au point  ... U (i,j) ...
c
c    pbary(i,j) = Pext(i,j  ) * ( alpha2(i,j  ) + alpha3(i,j  )     +
c                 Pext(i,j+1) * ( alpha1(i,j+1) + alpha4(i,j+1)  
c     localise  au point  ... V (i,j) ...
c
c  pbarxy(i,j)= Pext(i,j) *alpha2(i,j) + Pext(i+1,j) *alpha3(i+1,j) +
c               Pext(i,j+1)*alpha1(i,j+1)+ Pext(i+1,j+1)*alpha4(i+1,j+1)
c     localise  au point  ... Z (i,j) ...
c
c
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

      REAL pext( ip1jmp1 ),  pbarx ( ip1jmp1 )
      REAL pbary(  ip1jm  ),  pbarxy(  ip1jm  )

      INTEGER   ij



      DO 1 ij = 1, ip1jmp1 - 1
      pbarx( ij ) = pext(ij) * alpha1p2(ij) + pext(ij+1)*alpha3p4(ij+1)
   1  CONTINUE

c    .... correction pour pbarx( iip1,j) .....

c    ...    pbarx(iip1,j)= pbarx(1,j) ...
CDIR$ IVDEP
      DO 2 ij = iip1, ip1jmp1, iip1
      pbarx( ij ) = pbarx( ij - iim )
   2  CONTINUE


      DO 3 ij = 1,ip1jm
      pbary( ij ) = pext(   ij  )   * alpha2p3(   ij   )     +
     *              pext( ij+iip1 ) * alpha1p4( ij+iip1 )
   3  CONTINUE


      DO 5 ij = 1, ip1jm - 1
      pbarxy( ij ) = pext(ij)*alpha2(ij) + pext(ij+1)*alpha3(ij+1) +
     *   pext(ij+iip1)*alpha1(ij+iip1) + pext(ij+iip2)*alpha4(ij+iip2)
   5  CONTINUE


c    ....  correction pour     pbarxy( iip1,j )  ........

CDIR$ IVDEP

      DO 7 ij = iip1, ip1jm, iip1
      pbarxy( ij ) = pbarxy( ij - iim )
   7  CONTINUE


      RETURN
      END
