!
! $Id: undefSTD.F 1418 2010-07-19 15:11:24Z jghattas $
!
      SUBROUTINE undefSTD(itap,freq_calNMC, read_climoz)
      USE netcdf
      USE dimphy
      USE phys_state_var_mod ! Variables sauvegardees de la physique
      IMPLICIT none
c
c====================================================================
c
c I. Musat : 09.2004
c
c Calcul * du nombre de pas de temps (FLOAT(ecrit_XXX)-tnondef)) 
c          ou la variable tlevSTD est bien definie (.NE.missing_val), 
c et 
c        * de la somme de tlevSTD => tsumSTD
c
c nout=1 !var. journaliere "day" moyenne sur tous les pas de temps
c        ! de la physique
c nout=2 !var. mensuelle "mth" moyennee sur tous les pas de temps
c        ! de la physique
c nout=3 !var. mensuelle "NMC" moyennee toutes les ecrit_hf
c
c
c NB: mettre "inst(X)" dans le write_hist*NMC.h !
c====================================================================
c
cym#include "dimensions.h"
cym      integer jjmp1
cym      parameter (jjmp1=jjm+1-1/jjm)
cym#include "dimphy.h"
c variables Input
c
c     INTEGER nlevSTD, klevSTD, itap
c     PARAMETER(klevSTD=17)
      INTEGER itap
c     REAL dtime
c
c variables locales
c     INTEGER i, k, nout, n
c     PARAMETER(nout=3) !nout=1 : day; =2 : mth; =3 : NMC
      INTEGER i, k, n
      REAL freq_calNMC(nout)
      INTEGER read_climoz
c
c variables Output
c     REAL tlevSTD(klon,klevSTD), tsumSTD(klon,klevSTD,nout)
c     LOGICAL oknondef(klon,klevSTD,nout)
c     REAL tnondef(klon,klevSTD,nout)
c
      REAL missing_val
c
      missing_val=nf90_fill_real
c
      DO n=1, nout
c
c
c calcul variables tous les freq_calNMC(n)/dtime pas de temps 
c de la physique
c
       IF(MOD(itap,NINT(freq_calNMC(n)/dtime)).EQ.0) THEN
        DO k=1, nlevSTD
         DO i=1, klon
          IF(tlevSTD(i,k).EQ.missing_val) THEN
c          IF(oknondef(i,k,n)) THEN          
            tnondef(i,k,n)=tnondef(i,k,n)+1.
c          ENDIF !oknondef(i,k)
c
          ELSE IF(tlevSTD(i,k).NE.missing_val) THEN
           tsumSTD(i,k,n)=tsumSTD(i,k,n)+tlevSTD(i,k)
           usumSTD(i,k,n)=usumSTD(i,k,n)+ulevSTD(i,k)
           vsumSTD(i,k,n)=vsumSTD(i,k,n)+vlevSTD(i,k)
           wsumSTD(i,k,n)=wsumSTD(i,k,n)+wlevSTD(i,k)
           phisumSTD(i,k,n)=phisumSTD(i,k,n)+philevSTD(i,k)
           qsumSTD(i,k,n)=qsumSTD(i,k,n)+qlevSTD(i,k)
           rhsumSTD(i,k,n)=rhsumSTD(i,k,n)+rhlevSTD(i,k)
           uvsumSTD(i,k,n)=uvsumSTD(i,k,n)+uvSTD(i,k)
           vqsumSTD(i,k,n)=vqsumSTD(i,k,n)+vqSTD(i,k)
           vTsumSTD(i,k,n)=vTsumSTD(i,k,n)+vTSTD(i,k)
           wqsumSTD(i,k,n)=wqsumSTD(i,k,n)+wqSTD(i,k)
           vphisumSTD(i,k,n)=vphisumSTD(i,k,n)+vphiSTD(i,k)
           wTsumSTD(i,k,n)=wTsumSTD(i,k,n)+wTSTD(i,k)
           u2sumSTD(i,k,n)=u2sumSTD(i,k,n)+u2STD(i,k)
           v2sumSTD(i,k,n)=v2sumSTD(i,k,n)+v2STD(i,k)
           T2sumSTD(i,k,n)=T2sumSTD(i,k,n)+T2STD(i,k)
           O3sumSTD(i,k,n)=O3sumSTD(i,k,n)+O3STD(i,k)
           IF (read_climoz==2) 
     &          O3daysumSTD(i,k,n)=O3daysumSTD(i,k,n)+O3daySTD(i,k)

          ENDIF 
         ENDDO !i
        ENDDO !k
c
       ENDIF !MOD(itap,NINT(freq_calNMC(n)/dtime)).EQ.0
c
      ENDDO !n
c
      RETURN
      END  
