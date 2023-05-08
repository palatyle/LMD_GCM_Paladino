










!
! $Header: /home/cvsroot/LMDZ4/libf/filtrez/inifgn.F,v 1.1.1.1 2004-05-19 12:53:09 lmdzadmin Exp $
!
      SUBROUTINE inifgn(dv)
c  
c    ...  H.Upadyaya , O.Sharma  ... 
c
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
      REAL vec(iim,iim),vec1(iim,iim)
      REAL dlonu(iim),dlonv(iim)
      REAL du(iim),dv(iim),d(iim)
      REAL pi
      INTEGER i,j,k,imm1,nrot
C
!
! $Header$
!
      COMMON/coefils/jfiltnu,jfiltsu,jfiltnv,jfiltsv,sddu(iim),sddv(iim)&
     & ,unsddu(iim),unsddv(iim),coefilu(iim,jjm),coefilv(iim,jjm),      &
     & modfrstu(jjm),modfrstv(jjm),eignfnu(iim,iim),eignfnv(iim,iim)    &
     & ,coefilu2(iim,jjm),coefilv2(iim,jjm)
!c
      INTEGER jfiltnu,jfiltsu,jfiltnv,jfiltsv,modfrstu,modfrstv
      REAL    sddu,sddv,unsddu,unsddv,coefilu,coefilv,eignfnu,eignfnv
      REAL    coefilu2,coefilv2
c
      EXTERNAL SSUM, acc,eigen,jacobi
      REAL SSUM
c

      imm1  = iim -1
      pi = 2.* ASIN(1.)
C
      DO 5 i=1,iim
       dlonu(i)=  xprimu( i )
       dlonv(i)=  xprimv( i )
   5  CONTINUE

      DO 12 i=1,iim
      sddv(i)   = SQRT(dlonv(i))
      sddu(i)   = SQRT(dlonu(i))
      unsddu(i) = 1./sddu(i)
      unsddv(i) = 1./sddv(i)
  12  CONTINUE
C
      DO 17 j=1,iim
      DO 17 i=1,iim
      vec(i,j)     = 0.
      vec1(i,j)    = 0.
      eignfnv(i,j) = 0.
      eignfnu(i,j) = 0.
  17  CONTINUE
c
c
      eignfnv(1,1)    = -1.
      eignfnv(iim,1)  =  1.
      DO 20 i=1,imm1
      eignfnv(i+1,i+1)= -1.
      eignfnv(i,i+1)  =  1.
  20  CONTINUE
      DO 25 j=1,iim
      DO 25 i=1,iim
      eignfnv(i,j) = eignfnv(i,j)/(sddu(i)*sddv(j))
  25  CONTINUE
      DO 30 j=1,iim
      DO 30 i=1,iim
      eignfnu(i,j) = -eignfnv(j,i)
  30  CONTINUE
c
      DO j = 1, iim
      DO i = 1, iim
        vec (i,j) = 0.0
        vec1(i,j) = 0.0
       DO k = 1, iim
        vec (i,j) = vec(i,j)  + eignfnu(i,k) * eignfnv(k,j)
        vec1(i,j) = vec1(i,j) + eignfnv(i,k) * eignfnu(k,j)
       ENDDO
      ENDDO
      ENDDO

c
      CALL jacobi(vec,iim,iim,dv,eignfnv,nrot)
      CALL acc(eignfnv,d,iim)
      CALL eigen_sort(dv,eignfnv,iim,iim)
c
      CALL jacobi(vec1,iim,iim,du,eignfnu,nrot)
      CALL acc(eignfnu,d,iim)
      CALL eigen_sort(du,eignfnu,iim,iim)

cc   ancienne version avec appels IMSL
c
c     CALL MXM(eignfnu,iim,eignfnv,iim,vec,iim)
c     CALL MXM(eignfnv,iim,eignfnu,iim,vec1,iim)
c     CALL EVCSF(iim,vec,iim,dv,eignfnv,iim)
c     CALL acc(eignfnv,d,iim)
c     CALL eigen(eignfnv,dv)
c
c     CALL EVCSF(iim,vec1,iim,du,eignfnu,iim)
c     CALL acc(eignfnu,d,iim)
c     CALL eigen(eignfnu,du)

      RETURN
      END

