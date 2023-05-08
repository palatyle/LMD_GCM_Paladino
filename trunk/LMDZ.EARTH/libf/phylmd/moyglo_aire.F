!
! $Header$
!
      SUBROUTINE moyglo_pondaire(nhori, champ, aire, 
     .            ok_msk, msk, moyglo)
c
       USE dimphy
       IMPLICIT none
c
c ==================================================================
c I. Musat, 07.2004
c
c Calcul moyenne globale ponderee par l'aire totale, avec ou sans masque
c
c moyenne = Somme_(champ* aire)/Somme_aire
c
c ==================================================================
c
#include "dimensions.h"
cym#include "dimphy.h"
       INTEGER i, nhori
       REAL champ(klon), aire(klon), msk(klon)
       LOGICAL ok_msk 
       REAL moyglo
c
c var locale
       REAL airetot
c
c      PRINT*,'moyglo_pondaire nhori',nhori
c
       airetot=0.
       moyglo=0.
c
       IF(ok_msk) THEN
        DO i=1, nhori
c        IF(msk(i).EQ.1.) THEN 
         IF(msk(i).GT.0.) THEN 
c
c aire totale
          airetot=airetot+aire(i)*msk(i)
c
c ponderation par la masse
          moyglo=moyglo+champ(i)* aire(i)*msk(i) 
         ENDIF
        ENDDO
c
       ELSE !ok_msk
        DO i=1, nhori
c
c aire totale
          airetot=airetot+aire(i)
c
c ponderation par la masse
          moyglo=moyglo+champ(i)* aire(i) 
        ENDDO 
c
       ENDIF
c 
c moyenne ponderee par l'aire
       moyglo=moyglo/airetot
c
       RETURN 
       END
c
       SUBROUTINE moyglo_pondaima(nhori, nvert, champ,
     . aire, pbord, moyglo)
       USE dimphy
       IMPLICIT none
c ==================================================================
c I. Musat, 07.2004
c
c Calcul moyenne globale ponderee par la masse d'air, divisee par l'aire
c totale avec ou sans masque
c
c moyenne = Somme_(champ* masse_dair)/Somme_aire
c
c ==================================================================
#include "dimensions.h"
cym#include "dimphy.h"
#include "YOMCST.h"
       INTEGER i, k, nhori, nvert
       REAL champ(klon,klev), aire(klon)
       REAL pbord(klon,klev+1)
       REAL moyglo
c
c var locale
       REAL airetot
c
c      PRINT*,'moyglo_pondaima RG, nhori, nvert',RG,nhori,nvert
c
c ponderation par la masse
       moyglo=0.
       DO k=1, nvert
       DO i=1, nhori
        moyglo=moyglo+
     .  champ(i,k)*(pbord(i,k)-pbord(i,k+1))/RG*aire(i)
       ENDDO
       ENDDO
c
c aire totale
       airetot=0.
       DO i=1, nhori
        airetot=airetot+aire(i)
       ENDDO
c
c moyenne par mettre carre avec ponderation par la masse
       moyglo=moyglo/airetot
c
      RETURN
      END
c
       SUBROUTINE moyglo_pondmass(nhori, nvert, champ,
     . aire, pbord, moyglo)
       USE dimphy
       IMPLICIT none
c ==================================================================
c I. Musat, 07.2004
c
c Calcul moyenne globale ponderee par la masse d'air, divisee par la
c masse totale d'air, avec ou sans masque
c
c moyenne = Somme_(champ* masse_dair)/Somme_(masse_dair)
c
c ==================================================================
#include "dimensions.h"
cym#include "dimphy.h"
#include "YOMCST.h"
       INTEGER i, k, nhori, nvert
       REAL champ(klon,klev), aire(klon)
       REAL pbord(klon,klev+1)
       REAL moyglo
c
c var locale
       REAL massetot
c
c      PRINT*,'moyglo_pondmass RG, nhori, nvert',RG,nhori,nvert
c
c ponderation par la masse
       moyglo=0.
       DO k=1, nvert
       DO i=1, nhori
        moyglo=moyglo+
     .  champ(i,k)*(pbord(i,k)-pbord(i,k+1))/RG*aire(i)
       ENDDO
       ENDDO
c
c masse totale
       massetot=0.
       DO k=1, nvert
       DO i=1, nhori
        massetot=massetot+
     .  (pbord(i,k)-pbord(i,k+1))/RG*aire(i)
       ENDDO
       ENDDO
c
c moyenne par mettre carre avec ponderation par la masse
       moyglo=moyglo/massetot
c
      RETURN
      END
c
