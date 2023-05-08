










      SUBROUTINE massbar_p(  masse, massebx, masseby )
     
c
c **********************************************************************
c
c  Calcule les moyennes en x et  y de la masse d'air dans chaque maille.
c **********************************************************************
c    Auteurs : P. Le Van , Fr. Hourdin  .
c   ..........
c
c  ..  masse                 est  un argum. d'entree  pour le s-pg ...
c  ..  massebx,masseby      sont des argum. de sortie pour le s-pg ...
c     
c
      USE parallel_lmdz
      IMPLICIT NONE
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
c
      REAL    masse( ip1jmp1,llm ), massebx( ip1jmp1,llm )  ,
     *      masseby(   ip1jm,llm )
      INTEGER ij,l,ijb,ije
c
c
c   Methode pour calculer massebx et masseby .
c   ----------------------------------------
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
c                       On  a :
c
c    massebx(i,j) = masse(i  ,j) * ( alpha1(i  ,j) + alpha2(i,j))   +
c                   masse(i+1,j) * ( alpha3(i+1,j) + alpha4(i+1,j) )
c     localise  au point  ... U (i,j) ...
c
c    masseby(i,j) = masse(i,j  ) * ( alpha2(i,j  ) + alpha3(i,j  )  +
c                   masse(i,j+1) * ( alpha1(i,j+1) + alpha4(i,j+1)  
c     localise  au point  ... V (i,j) ...
c
c
c=======================================================================
      
      
      
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)  
      DO   100    l = 1 , llm
c
        ijb=ij_begin
        ije=ij_end+iip1
        if (pole_sud) ije=ije-iip1
        
        DO  ij = ijb, ije - 1
         massebx(ij,l) =  masse( ij, l) * alpha1p2( ij  )     + 
     *                   masse(ij+1, l) * alpha3p4(ij+1 )
        ENDDO

c    .... correction pour massebx( iip1,j) .....
c    ...    massebx(iip1,j)= massebx(1,j) ...
c
CDIR$ IVDEP

        

        DO  ij = ijb+iim, ije+iim, iip1
         massebx( ij,l ) = massebx( ij - iim,l )
        ENDDO


      
        ijb=ij_begin-iip1
        ije=ij_end+iip1
        if (pole_nord) ijb=ij_begin
        if (pole_sud) ije=ij_end-iip1

         DO  ij = ijb,ije
         masseby( ij,l ) = masse(  ij   , l ) * alpha2p3(   ij    )  +
     *                     masse(ij+iip1, l ) * alpha1p4( ij+iip1 )
         ENDDO

100   CONTINUE
c$OMP END DO NOWAIT
c
      RETURN
      END
