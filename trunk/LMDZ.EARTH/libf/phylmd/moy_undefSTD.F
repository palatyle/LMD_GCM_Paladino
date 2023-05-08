!
! $Id: moy_undefSTD.F 1403 2010-07-01 09:02:53Z fairhead $
!
      SUBROUTINE moy_undefSTD(itap,freq_outNMC,freq_moyNMC)
      USE netcdf
      USE dimphy
      USE phys_state_var_mod ! Variables sauvegardees de la physique
      IMPLICIT none
c
c====================================================================
c
c I. Musat : 09.2004
c
c Moyenne - a des frequences differentes - des valeurs bien definies
c         (.NE.missing_val) des variables interpolees a un niveau de
c         pression.
c 1) les variables de type "day" (nout=1) ou "mth" (nout=2) sont sommees
c    tous les pas de temps de la physique
c
c 2) les variables de type "NMC" (nout=3) sont calculees a partir
c    des valeurs instantannees toutes les 6 heures
c
c
c NB: mettre "inst(X)" dans le write_hist*NMC.h !
c====================================================================
cym#include "dimensions.h"
cym      integer jjmp1
cym      parameter (jjmp1=jjm+1-1/jjm)
cym#include "dimphy.h"
c
c
c variables Input
c     INTEGER nlevSTD, klevSTD, itap
c     PARAMETER(klevSTD=17)
      INTEGER itap
c
c variables locales
c     INTEGER i, k, nout, n
c     PARAMETER(nout=3) !nout=1 day/nout=2 mth/nout=3 NMC
      INTEGER i, k, n
c     REAL dtime, freq_outNMC(nout), freq_moyNMC(nout)
      REAL freq_outNMC(nout), freq_moyNMC(nout)
c
c variables Output
c     REAL tnondef(klon,klevSTD,nout)
c     REAL tsumSTD(klon,klevSTD,nout)
c
      REAL missing_val
c
      missing_val=nf90_fill_real
c
      DO n=1, nout
c 
c calcul 1 fois par jour
c
       IF(MOD(itap,NINT(freq_outNMC(n)/dtime)).EQ.0) THEN
c
        DO k=1, nlevSTD
         DO i=1, klon
          IF(tnondef(i,k,n).NE.(freq_moyNMC(n))) THEN
           tsumSTD(i,k,n)=tsumSTD(i,k,n)/
     $     (freq_moyNMC(n)-tnondef(i,k,n))
cIM BEG
          usumSTD(i,k,n)=usumSTD(i,k,n)/
     $     (freq_moyNMC(n)-tnondef(i,k,n))
          vsumSTD(i,k,n)=vsumSTD(i,k,n)/
     $     (freq_moyNMC(n)-tnondef(i,k,n))
          wsumSTD(i,k,n)=wsumSTD(i,k,n)/
     $     (freq_moyNMC(n)-tnondef(i,k,n))
          phisumSTD(i,k,n)=phisumSTD(i,k,n)/
     $     (freq_moyNMC(n)-tnondef(i,k,n))
          qsumSTD(i,k,n)=qsumSTD(i,k,n)/
     $     (freq_moyNMC(n)-tnondef(i,k,n))
          rhsumSTD(i,k,n)=rhsumSTD(i,k,n)/
     $     (freq_moyNMC(n)-tnondef(i,k,n))
          uvsumSTD(i,k,n)=uvsumSTD(i,k,n)/
     $     (freq_moyNMC(n)-tnondef(i,k,n))
          vqsumSTD(i,k,n)=vqsumSTD(i,k,n)/
     $     (freq_moyNMC(n)-tnondef(i,k,n))
          vTsumSTD(i,k,n)=vTsumSTD(i,k,n)/
     $     (freq_moyNMC(n)-tnondef(i,k,n))
          wqsumSTD(i,k,n)=wqsumSTD(i,k,n)/
     $     (freq_moyNMC(n)-tnondef(i,k,n))
          vphisumSTD(i,k,n)=vphisumSTD(i,k,n)/
     $     (freq_moyNMC(n)-tnondef(i,k,n))
          wTsumSTD(i,k,n)=wTsumSTD(i,k,n)/
     $     (freq_moyNMC(n)-tnondef(i,k,n))
          u2sumSTD(i,k,n)=u2sumSTD(i,k,n)/
     $     (freq_moyNMC(n)-tnondef(i,k,n))
          v2sumSTD(i,k,n)=v2sumSTD(i,k,n)/
     $     (freq_moyNMC(n)-tnondef(i,k,n))
          T2sumSTD(i,k,n)=T2sumSTD(i,k,n)/
     $     (freq_moyNMC(n)-tnondef(i,k,n))
          O3sumSTD(i,k,n)=O3sumSTD(i,k,n)/
     $     (freq_moyNMC(n)-tnondef(i,k,n))
          O3daysumSTD(i,k,n)=O3daysumSTD(i,k,n)/
     $     (freq_moyNMC(n)-tnondef(i,k,n))
cIM END
          ELSE
           tsumSTD(i,k,n)=missing_val
           usumSTD(i,k,n)=missing_val
           vsumSTD(i,k,n)=missing_val
           wsumSTD(i,k,n)=missing_val
           phisumSTD(i,k,n)=missing_val
           qsumSTD(i,k,n)=missing_val
           rhsumSTD(i,k,n)=missing_val
           uvsumSTD(i,k,n)=missing_val
           vqsumSTD(i,k,n)=missing_val
           vTsumSTD(i,k,n)=missing_val
           wqsumSTD(i,k,n)=missing_val
           vphisumSTD(i,k,n)=missing_val
           wTsumSTD(i,k,n)=missing_val
           u2sumSTD(i,k,n)=missing_val
           v2sumSTD(i,k,n)=missing_val
           T2sumSTD(i,k,n)=missing_val
           O3sumSTD(i,k,n)=missing_val
           O3daysumSTD(i,k,n)=missing_val
          ENDIF !tnondef(i,k,n).NE.(freq_moyNMC(n))
         ENDDO !i
        ENDDO !k
       ENDIF !MOD(itap,NINT(freq_outNMC(n)/dtime)).EQ.0
c
      ENDDO !n
c
      RETURN
      END  
