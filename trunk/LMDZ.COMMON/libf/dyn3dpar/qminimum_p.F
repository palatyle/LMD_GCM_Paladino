      SUBROUTINE qminimum_p( q,nq,deltap )
      USE parallel_lmdz
      USE comvert_mod, ONLY: presnivs
      IMPLICIT none
c
c  -- Objet : Traiter les valeurs trop petites (meme negatives)
c             pour l'eau vapeur et l'eau liquide
c
#include "dimensions.h"
#include "paramet.h"
c
      INTEGER nq
      REAL q(ip1jmp1,llm,nq), deltap(ip1jmp1,llm)
c
      INTEGER iq_vap, iq_liq
      PARAMETER ( iq_vap = 1 ) ! indice pour l'eau vapeur
      PARAMETER ( iq_liq = 2 ) ! indice pour l'eau liquide
      REAL seuil_vap, seuil_liq
      PARAMETER ( seuil_vap = 1.0e-10 ) ! seuil pour l'eau vapeur
      PARAMETER ( seuil_liq = 1.0e-11 ) ! seuil pour l'eau liquide
c
c  NB. ....( Il est souhaitable mais non obligatoire que les valeurs des
c            parametres seuil_vap, seuil_liq soient pareilles a celles 
c            qui  sont utilisees dans la routine    ADDFI       )
c     .................................................................
c
      INTEGER i, k, iq
      REAL zx_defau, zx_abc, zx_pump(ip1jmp1), pompe
c
      REAL SSUM
      EXTERNAL SSUM
c
      INTEGER imprim
      SAVE imprim
      DATA imprim /0/
c$OMP THREADPRIVATE(imprim)
      INTEGER ijb,ije
      INTEGER Index_pump(ip1jmp1)
      INTEGER nb_pump
c
c Quand l'eau liquide est trop petite (ou negative), on prend
c l'eau vapeur de la meme couche et la convertit en eau liquide
c (sans changer la temperature !)
c

      ijb=ij_begin
      ije=ij_end

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)       
      DO 1000 k = 1, llm
      DO 1040 i = ijb, ije
            if (seuil_liq - q(i,k,iq_liq) .gt. 0.d0 ) then
               q(i,k,iq_vap) = q(i,k,iq_vap) + q(i,k,iq_liq) - seuil_liq
               q(i,k,iq_liq) = seuil_liq
            endif
 1040 CONTINUE
 1000 CONTINUE
c$OMP END DO NOWAIT
c$OMP BARRIER
c --->  SYNCHRO OPENMP ICI

c
c Quand l'eau vapeur est trop faible (ou negative), on complete
c le defaut en prennant de l'eau vapeur de la couche au-dessous.
c
      iq = iq_vap
c
      DO k = llm, 2, -1
ccc      zx_abc = dpres(k) / dpres(k-1)
c$OMP DO SCHEDULE(STATIC)
      DO i = ijb, ije
         if ( seuil_vap - q(i,k,iq) .gt. 0.d0 ) then
            q(i,k-1,iq) =  q(i,k-1,iq) - ( seuil_vap - q(i,k,iq) ) *
     &           deltap(i,k) / deltap(i,k-1)
            q(i,k,iq)   =  seuil_vap  
         endif
      ENDDO
c$OMP END DO NOWAIT
      ENDDO
c$OMP BARRIER
c
c Quand il s'agit de la premiere couche au-dessus du sol, on
c doit imprimer un message d'avertissement (saturation possible).
c
      nb_pump=0
c$OMP DO SCHEDULE(STATIC)
      DO i = ijb, ije
         zx_pump(i) = AMAX1( 0.0, seuil_vap - q(i,1,iq) )
         q(i,1,iq)  = AMAX1( q(i,1,iq), seuil_vap )
         IF (zx_pump(i) > 0.0) THEN
            nb_pump = nb_pump+1
            Index_pump(nb_pump)=i
         ENDIF
      ENDDO
c$OMP END DO  
!      pompe = SSUM(ije-ijb+1,zx_pump(ijb),1)

      IF (imprim.LE.100 .AND. nb_pump .GT. 0 ) THEN
         PRINT *, 'ATT!:on pompe de l eau au sol'
         DO i = 1, nb_pump
               imprim = imprim + 1
               PRINT*,'  en ',index_pump(i),zx_pump(index_pump(i))
         ENDDO
      ENDIF
c
      RETURN
      END
