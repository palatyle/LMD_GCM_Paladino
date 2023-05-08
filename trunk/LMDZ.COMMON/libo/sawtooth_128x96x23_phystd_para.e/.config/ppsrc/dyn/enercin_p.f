










      SUBROUTINE enercin_p ( vcov, ucov, vcont, ucont, ecin )
      USE parallel_lmdz
      IMPLICIT NONE

c=======================================================================
c
c   Auteur: P. Le Van
c   -------
c
c   Objet:
c   ------
c
c *********************************************************************
c .. calcul de l'energie cinetique aux niveaux s  ......
c *********************************************************************
c  vcov, vcont, ucov et ucont sont des arguments d'entree pour le s-pg .
c  ecin         est  un  argument de sortie pour le s-pg
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

      REAL vcov( ip1jm,llm ),vcont( ip1jm,llm ),
     * ucov( ip1jmp1,llm ),ucont( ip1jmp1,llm ),ecin( ip1jmp1,llm )

      REAL ecinni( iip1 ),ecinsi( iip1 )

      REAL ecinpn, ecinps
      INTEGER     l,ij,i,ijb,ije

      EXTERNAL    SSUM
      REAL        SSUM



c                 . V
c                i,j-1

c      alpha4 .       . alpha1


c        U .      . P     . U
c       i-1,j    i,j      i,j

c      alpha3 .       . alpha2


c                 . V
c                i,j

c    
c  L'energie cinetique au point scalaire P(i,j) ,autre que les poles, est :
c       Ecin = 0.5 * U(i-1,j)**2 *( alpha3 + alpha4 )  +
c              0.5 * U(i  ,j)**2 *( alpha1 + alpha2 )  +
c              0.5 * V(i,j-1)**2 *( alpha1 + alpha4 )  +
c              0.5 * V(i,  j)**2 *( alpha2 + alpha3 )

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO 5 l = 1,llm
      
      ijb=ij_begin
      ije=ij_end+iip1
      
      IF (pole_nord) ijb=ij_begin+iip1
      IF (pole_sud)  ije=ij_end-iip1
      
      DO 1  ij = ijb, ije -1
      ecin( ij+1, l )  =    0.5  *
     * (   ucov( ij   ,l ) * ucont( ij   ,l ) * alpha3p4( ij +1 )   +
     *     ucov( ij+1 ,l ) * ucont( ij+1 ,l ) * alpha1p2( ij +1 )   +
     *     vcov(ij-iim,l ) * vcont(ij-iim,l ) * alpha1p4( ij +1 )   +
     *     vcov( ij+ 1,l ) * vcont( ij+ 1,l ) * alpha2p3( ij +1 )   )
   1  CONTINUE

c    ... correction pour  ecin(1,j,l)  ....
c    ...   ecin(1,j,l)= ecin(iip1,j,l) ...

CDIR$ IVDEP
      DO 2 ij = ijb, ije, iip1
      ecin( ij,l ) = ecin( ij + iim, l )
   2  CONTINUE

c     calcul aux poles  .......

      IF (pole_nord) THEN
    
        DO  i = 1, iim
         ecinni(i) = vcov(    i  ,  l) * 
     *               vcont(    i    ,l) * aire(   i   )
        ENDDO

        ecinpn = 0.5 * SSUM( iim,ecinni,1 ) / apoln

        DO ij = 1,iip1
          ecin(   ij     , l ) = ecinpn
        ENDDO
   
      ENDIF

      IF (pole_sud) THEN
    
        DO  i = 1, iim
         ecinsi(i) = vcov(i+ip1jmi1,l)* 
     *               vcont(i+ip1jmi1,l) * aire(i+ip1jm)
        ENDDO

        ecinps = 0.5 * SSUM( iim,ecinsi,1 ) / apols

        DO ij = 1,iip1
          ecin( ij+ ip1jm, l ) = ecinps
        ENDDO
   
      ENDIF

      
   5  CONTINUE
c$OMP END DO NOWAIT
      RETURN
      END
