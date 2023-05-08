










!
! $Header$
!
        subroutine interpost(q,qppm)

       implicit none


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

c Arguments   
      real   q(iip1,jjp1,llm)
      real   qppm(iim,jjp1,llm)
c Local
      integer l,i,j
  
c RE-INVERSION DES NIVEAUX
c le programme ppm3d travaille avec une 3ème coordonnée inversée par rapport
c de celle du LMDZ: z=1<=>niveau max, z=llm+1<=>surface
c On passe donc des niveaux de Lin à ceux du LMDZ
           
        do l=1,llm
          do j=1,jjp1
             do i=1,iim
                 q(i,j,l)=qppm(i,j,llm-l+1)
             enddo
          enddo
         enddo
            
c BOUCLAGE EN LONGITUDE PAS EFFECTUE DANS PPM3D

         do l=1,llm
           do j=1,jjp1
            q(iip1,j,l)=q(1,j,l)
           enddo
         enddo
  
      
       return

       end
