!
! $Header$
!
      SUBROUTINE homogene(paprs, q, dq, u,v, du, dv)
      USE dimphy
      IMPLICIT NONE
c==============================================================
c Schema ad hoc du melange vertical pour les vitesses u et v,
c a appliquer apres le schema de convection (fiajc et fiajh).
c
c paprs:input, pression demi-couche (inter-couche)
c q:    input, vapeur d'eau (kg/kg)
c dq:   input, incrementation de vapeur d'eau (de la convection)
c u:    input, vitesse u
c v:    input, vitesse v
c
c du:   output, incrementation pour u
c dv:   output, incrementation pour v
c==============================================================
cym#include "dimensions.h"
cym#include "dimphy.h"
c
      REAL paprs(klon,klev+1)
      REAL q(klon,klev), dq(klon,klev)
      REAL u(klon,klev), du(klon,klev)
      REAL v(klon,klev), dv(klon,klev)
c
      REAL zm_dq(klon) ! quantite totale de l'eau deplacee
      REAL zm_q(klon)  ! quantite totale de la vapeur d'eau
      REAL zm_u(klon)  ! moyenne de u (brassage parfait et total)
      REAL zm_v(klon)  ! moyenne de v (brassage parfait et total)
      REAL z_frac(klon) ! fraction du brassage parfait et total
      REAL zm_dp(klon)
c
      REAL zx
      INTEGER i, k
      REAL frac_max
      PARAMETER (frac_max=0.1)
      REAL seuil
      PARAMETER (seuil=1.0e-10)
      LOGICAL faisrien
      PARAMETER (faisrien=.true.)
c
      DO k = 1, klev
      DO i = 1, klon
         du(i,k) = 0.0
         dv(i,k) = 0.0
      ENDDO
      ENDDO
c
      IF (faisrien) RETURN
c
      DO i = 1, klon
         zm_dq(i)=0.
         zm_q(i) =0.
         zm_u(i)=0.
         zm_v(i)=0.
         zm_dp(i)=0.
      ENDDO
      DO k = 1, klev
      DO i = 1, klon
      IF (ABS(dq(i,k)).GT.seuil) THEN
         zx = paprs(i,k) - paprs(i,k+1)
         zm_dq(i) = zm_dq(i) + ABS(dq(i,k))*zx
         zm_q(i) = zm_q(i) + q(i,k)*zx
         zm_dp(i) = zm_dp(i) + zx
         zm_u(i) = zm_u(i) + u(i,k)*zx
         zm_v(i) = zm_v(i) + v(i,k)*zx
      ENDIF
      ENDDO
      ENDDO
c
c Hypothese principale: apres la convection, la vitesse de chaque
c couche est composee de deux parties: celle (1-z_frac) de la vitesse 
c original et celle (z_frac) de la vitesse moyenne qui serait la
c vitesse de chaque couche si le brassage etait parfait et total.
c La fraction du brassage est calculee par le rapport entre la quantite
c totale de la vapeur d'eau deplacee (ou condensee) et la quantite
c totale de la vapeur d'eau. Et cette fraction est limitee a frac_max 
c (Est-ce vraiment raisonnable ? Z.X. Li, le 07-09-1995).
c
      DO i = 1, klon
      IF (zm_dp(i).GE.1.0E-15 .AND. zm_q(i).GE.1.0E-15) THEN
         z_frac(i)=MIN(frac_max,zm_dq(i)/zm_q(i))
         zm_u(i)=zm_u(i)/zm_dp(i)
         zm_v(i)=zm_v(i)/zm_dp(i)
      ENDIF
      ENDDO
      DO k = 1, klev
      DO i = 1, klon
      IF (zm_dp(i).GE.1.e-15 .AND. zm_q(i).GE.1.e-15
     .                         .AND. ABS(dq(i,k)).GT.seuil) THEN
         du(i,k) = u(i,k)*(1.-z_frac(i)) + zm_u(i)*z_frac(i) - u(i,k)
         dv(i,k) = v(i,k)*(1.-z_frac(i)) + zm_v(i)*z_frac(i) - v(i,k)
      ENDIF
      ENDDO
      ENDDO
c
      RETURN
      END
