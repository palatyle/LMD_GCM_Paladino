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
#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
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
