










      SUBROUTINE tidal_forces (t, du, dv)

      IMPLICIT NONE
c
c=======================================================================
c
c   Auteur:  B. Charnay  (10/2010)
c   -------
c
c   Objet:
c   ------
c
c   *****************************************************************
c   ..... calcul du gradient horizontal du potentiel gravitationnel du aux forces de marees causees par Saturne
c   ..... Formule tiree de Tokano 2002
c   *****************************************************************
c          Ces termes sont ajoutes a  d(ucov)/dt et a d(vcov)/dt  ..
c
c
c    du et dv          sont des arguments de sortie pour le s-pg  ....
c
c=======================================================================
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
!#include "comorbit.h"
      REAL t        ! jour de l'annee
      REAL du( ip1jmp1,llm ),  dv( ip1jm,llm )

c     variables locales
      REAL Vo
      PARAMETER (Vo=-4.691e-6)
      INTEGER  l,ij,i,k
      REAL n                ! 2pi/periode de rotation siderale (en jours)
      REAL a0               ! angle à l'instant initial entre Titan et le perihelie
      PARAMETER (a0=0.)      

c     cos et sin de la latitude et longitude, calcules au premiers appel
      REAL coslonv(ip1jm),sinlonv(ip1jm)
      REAL sinlatv(ip1jm),coslatv(ip1jm)
      REAL coslonu(ip1jmp1),sinlonu(ip1jmp1)
      REAL sinlatu(ip1jmp1),coslatu(ip1jmp1)      

      LOGICAl first      

      SAVE coslonv,coslonu,sinlonu,sinlonv
      SAVE coslatv,coslatu,sinlatu,sinlatv
      SAVE first, n

      DATA first /.true./

! Calcul des sin et cos aux points consideres

      IF(first) THEN
         first=.false.
         n=2*3.145!*(1+1/673.)
         do i=1,iip1
          do k=1,jjm
            coslonv(i+(k-1)*iip1)=cos(rlonv(i))
            sinlonv(i+(k-1)*iip1)=sin(rlonv(i))
            coslatv(i+(k-1)*iip1)=cos(rlatv(k))
            sinlatv(i+(k-1)*iip1)=sin(rlatv(k))
          ENDDO
         ENDDO



         do i=1,iip1
          do k=1,jjp1
            coslonu(i+(k-1)*iip1)=cos(rlonu(i))
            sinlonu(i+(k-1)*iip1)=sin(rlonu(i))
            coslatu(i+(k-1)*iip1)=cos(rlatu(k))
            sinlatu(i+(k-1)*iip1)=sin(rlatu(k))
          ENDDO
         ENDDO



      ENDIF


! Tendance du aux forces de maree

      DO l = 1,llm

      DO ij  = 1, ip1jmp1 

       du(ij,l) = cu(ij)*Vo
     $    *(3*sinlonu(ij)*coslonu(ij)*coslatu(ij)*cos(n*t+a0) 
     $    -2*coslatu(ij)*(2*coslonu(ij)**2-1)*sin(n*t+a0))        
      ENDDO

      DO ij  = 1, ip1jm 
       dv(ij,l) = cv(ij)*Vo
     $    *(3*sinlatv(ij)*coslatv(ij)*coslonv(ij)**2*cos(n*t+a0) 
     $    + 4*coslatv(ij)*sinlatv(ij)*sinlonv(ij)*coslonv(ij)
     $    *sin(n*t+a0))      
      ENDDO

      ENDDO


c
      RETURN
      END
