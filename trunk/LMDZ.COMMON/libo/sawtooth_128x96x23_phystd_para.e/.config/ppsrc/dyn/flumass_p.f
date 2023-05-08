










      SUBROUTINE flumass_p(massebx,masseby, vcont, ucont, pbaru, pbarv)
      USE parallel_lmdz
      IMPLICIT NONE

c=======================================================================
c
c   Auteurs:  P. Le Van, F. Hourdin  .
c   -------
c
c   Objet:
c   ------
c
c *********************************************************************
c     .... calcul du flux de masse  aux niveaux s ......
c *********************************************************************
c   massebx,masseby,vcont et ucont sont des argum. d'entree pour le s-pg .
c       pbaru  et pbarv            sont des argum.de sortie pour le s-pg .
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

      REAL massebx( ip1jmp1,llm ),masseby( ip1jm,llm ) ,
     * vcont( ip1jm,llm ),ucont( ip1jmp1,llm ),pbaru( ip1jmp1,llm ),
     * pbarv( ip1jm,llm )

      REAL apbarun( iip1 ),apbarus( iip1 )

      REAL sairen,saireun,saires,saireus,ctn,cts,ctn0,cts0
      INTEGER  l,ij,i
      INTEGER ijb,ije
      
      EXTERNAL   SSUM
      REAL       SSUM
      
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)    
      DO  5 l = 1,llm

        ijb=ij_begin
        ije=ij_end+iip1
      
        if (pole_nord) ijb=ij_begin+iip1
        if (pole_sud)  ije=ij_end-iip1
        
        DO  1 ij = ijb,ije
          pbaru( ij,l ) = massebx( ij,l ) * ucont( ij,l )
   1    CONTINUE

        ijb=ij_begin-iip1
        ije=ij_end+iip1
      
        if (pole_nord) ijb=ij_begin
        if (pole_sud)  ije=ij_end-iip1
        
        DO 3 ij = ijb,ije
          pbarv( ij,l ) = masseby( ij,l ) * vcont( ij,l )
   3    CONTINUE

   5  CONTINUE
c$OMP END DO NOWAIT
c    ................................................................
c     calcul de la composante du flux de masse en x aux poles .......
c    ................................................................
c     par la resolution d'1 systeme de 2 equations .

c     la premiere equat.decrivant le calcul de la divergence en 1 point i
c     du pole,ce calcul etant itere de i=1 a i=im .
c                 c.a.d   ,
c     ( ( 0.5*pbaru(i)-0.5*pbaru(i-1) - pbarv(i))/aire(i)   =
c                                           - somme de ( pbarv(n) )/aire pole

c     l'autre equat.specifiant que la moyenne du flux de masse au pole est =0.
c     c.a.d    somme de pbaru(n)*aire locale(n) = 0.

c     on en revient ainsi a determiner la constante additive commune aux pbaru
c     qui representait pbaru(0,j,l) dans l'equat.du calcul de la diverg.au pt
c     i=1 .
c     i variant de 1 a im
c     n variant de 1 a im

      IF (pole_nord) THEN
     
        sairen = SSUM( iim,  aire(   1     ), 1 )
        saireun= SSUM( iim, aireu(   1     ), 1 )

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)    
        DO l = 1,llm
 
          ctn =  SSUM( iim, pbarv(    1     ,l),  1 )/ sairen
      
          pbaru(1,l)=pbarv(1,l) - ctn * aire(1)
        
          DO i = 2,iim
            pbaru(i,l) = pbaru(i- 1,l )    +
     *                   pbarv(i,l) - ctn * aire(i )
          ENDDO
        
          DO i = 1,iim
            apbarun(i) = aireu(    i   ) * pbaru(   i    , l)
          ENDDO
      
          ctn0 = -SSUM( iim,apbarun,1 )/saireun
        
          DO i = 1,iim
            pbaru(   i    , l) = 2. * ( pbaru(   i    , l) + ctn0 )
          ENDDO
       
          pbaru(   iip1 ,l ) = pbaru(    1    ,l )
        
        ENDDO
c$OMP END DO NOWAIT              

      ENDIF

      
      IF (pole_sud) THEN
  
        saires = SSUM( iim,  aire( ip1jm+1 ), 1 )
        saireus= SSUM( iim, aireu( ip1jm+1 ), 1 )

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)    
        DO  l = 1,llm
 
          cts =  SSUM( iim, pbarv(ip1jmi1+ 1,l),  1 )/ saires
          pbaru(ip1jm+1,l)= - pbarv(ip1jmi1+1,l) + cts * aire(ip1jm+1)
   
          DO i = 2,iim
            pbaru(i+ ip1jm,l) = pbaru(i+ip1jm-1,l)    -
     *                          pbarv(i+ip1jmi1,l)+cts*aire(i+ip1jm)
          ENDDO
        
          DO i = 1,iim
            apbarus(i) = aireu(i +ip1jm) * pbaru(i +ip1jm, l)
          ENDDO

          cts0 = -SSUM( iim,apbarus,1 )/saireus

          DO i = 1,iim
            pbaru(i+ ip1jm, l) = 2. * ( pbaru(i +ip1jm, l) + cts0 )
          ENDDO

          pbaru( ip1jmp1,l ) = pbaru( ip1jm +1,l )
       
        ENDDO
c$OMP END DO NOWAIT         
      ENDIF
      
      RETURN
      END
