!
! $Id: ini_undefSTD.F 1403 2010-07-01 09:02:53Z fairhead $
!

      SUBROUTINE ini_undefSTD(itap,
     $           freq_outNMC)
      USE dimphy
      USE phys_state_var_mod ! Variables sauvegardees de la physique
      IMPLICIT none
c
c====================================================================
c
c I. Musat : 09.2004
c
c Initialisation - a des frequences differentes : 
c
c 1) des variables moyennees sur la journee "day" ou sur le mois "mth"
c    calculees a partir des valeurs "instantannees" de la physique
c
c 2) des variables moyennes mensuelles "NMC" calculees a partir des val.
c    toutes les 6 heures
c
c nout=1 !var. journaliere "day" moyenne sur tous les pas de temps
c              ! de la physique
c nout=2 !var. mensuelle "mth" moyennee sur tous les pas de temps
c              ! de la physique
c nout=3 !var. mensuelle "NMC" moyennee toutes les 6heures
c
c NB: mettre "inst(X)" dans le write_hist*NMC.h !
c====================================================================
c
cym #include "dimensions.h"
cym      integer jjmp1
cym      parameter (jjmp1=jjm+1-1/jjm)
cym #include "dimphy.h"
c variables Input/Output
c     INTEGER nlevSTD, klevSTD, itap
      INTEGER itap
c     PARAMETER(klevSTD=17)
c     REAL dtime
c
c variables locales
c     INTEGER i, k, nout, n
      INTEGER i, k, n
c     PARAMETER(nout=3) !nout=1 day/nout=2 mth/nout=3 NMC
      REAL freq_outNMC(nout)
c
c variables Output
c     REAL tnondef(klon,klevSTD,nout)
c     REAL tsumSTD(klon,klevSTD,nout)
c
      DO n=1, nout
c
c initialisation variables en debut de la journee ou du mois
c
       IF(MOD(itap,NINT(freq_outNMC(n)/dtime)).EQ.1.) THEN
        DO k=1, nlevSTD
         DO i=1, klon
          tnondef(i,k,n)=0.
          tsumSTD(i,k,n)=0.
          usumSTD(i,k,n)=0.
          vsumSTD(i,k,n)=0.
          wsumSTD(i,k,n)=0.
          phisumSTD(i,k,n)=0.
          qsumSTD(i,k,n)=0.
          rhsumSTD(i,k,n)=0.
          uvsumSTD(i,k,n)=0.
          vqsumSTD(i,k,n)=0.
          vTsumSTD(i,k,n)=0.
          wqsumSTD(i,k,n)=0.
          vphisumSTD(i,k,n)=0.
          wTsumSTD(i,k,n)=0.
          u2sumSTD(i,k,n)=0.
          v2sumSTD(i,k,n)=0.
          T2sumSTD(i,k,n)=0.
          O3sumSTD(i,k,n)=0.
          O3daysumSTD(i,k,n)=0.
         ENDDO !i
        ENDDO !k
c
       ENDIF !MOD(itap,NINT(freq_outNMC(n)/dtime)).EQ.1.
c
      ENDDO !n
c
      RETURN
      END  
