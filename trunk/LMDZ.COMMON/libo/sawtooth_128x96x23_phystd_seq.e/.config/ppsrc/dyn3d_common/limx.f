










!
! $Header$
!
      SUBROUTINE limx(s0,sx,sm,pente_max)
c
c     Auteurs:   P.Le Van, F.Hourdin, F.Forget 
c
c    ********************************************************************
c     Shema  d'advection " pseudo amont " .
c    ********************************************************************
c     nq,iq,q,pbaru,pbarv,w sont des arguments d'entree  pour le s-pg ....
c
c
c   --------------------------------------------------------------------
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
c
c   Arguments:
c   ----------
      real pente_max
      REAL s0(ip1jmp1,llm),sm(ip1jmp1,llm)
      real sx(ip1jmp1,llm)
c
c      Local 
c   ---------
c
      INTEGER ij,l,j,i,iju,ijq,indu(ip1jmp1),niju
      integer n0,iadvplus(ip1jmp1,llm),nl(llm)
c
      REAL q(ip1jmp1,llm)
      real dxq(ip1jmp1,llm)


      REAL new_m,zm
      real dxqu(ip1jmp1)
      real adxqu(ip1jmp1),dxqmax(ip1jmp1)

      Logical extremum,first
      save first

      REAL      SSUM,CVMGP,CVMGT
      integer ismax,ismin
      EXTERNAL  SSUM, convflu,ismin,ismax
      EXTERNAL filtreg

      data first/.true./


       DO  l = 1,llm
         DO  ij=1,ip1jmp1
               q(ij,l) = s0(ij,l) / sm ( ij,l )
               dxq(ij,l) = sx(ij,l) /sm(ij,l)
         ENDDO
       ENDDO

c   calcul de la pente a droite et a gauche de la maille

      do l = 1, llm
         do ij=iip2,ip1jm-1
            dxqu(ij)=q(ij+1,l)-q(ij,l)
         enddo
         do ij=iip1+iip1,ip1jm,iip1
            dxqu(ij)=dxqu(ij-iim)
         enddo

         do ij=iip2,ip1jm
            adxqu(ij)=abs(dxqu(ij))
         enddo

c   calcul de la pente maximum dans la maille en valeur absolue

         do ij=iip2+1,ip1jm
            dxqmax(ij)=pente_max*min(adxqu(ij-1),adxqu(ij))
         enddo

         do ij=iip1+iip1,ip1jm,iip1
            dxqmax(ij-iim)=dxqmax(ij)
         enddo

c   calcul de la pente avec limitation

         do ij=iip2+1,ip1jm
            if(     dxqu(ij-1)*dxqu(ij).gt.0.
     &         .and. dxq(ij,l)*dxqu(ij).gt.0.) then
              dxq(ij,l)=
     &         sign(min(abs(dxq(ij,l)),dxqmax(ij)),dxq(ij,l))
            else
c   extremum local
               dxq(ij,l)=0.
            endif
         enddo
         do ij=iip1+iip1,ip1jm,iip1
            dxq(ij-iim,l)=dxq(ij,l)
         enddo

         DO  ij=1,ip1jmp1
               sx(ij,l) = dxq(ij,l)*sm(ij,l)
         ENDDO

       ENDDO

      RETURN
      END
