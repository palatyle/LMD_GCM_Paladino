!
! $Id: o3cm.F 1403 2010-07-01 09:02:53Z fairhead $
!
      SUBROUTINE o3cm (amb, bmb, sortie, ntab)
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS) date: 19930818
c Objet: Ce programme calcule le contenu en ozone "sortie"
c        (unite: cm.atm) entre deux niveaux "amb" et "bmb" (unite: mb)
c        "ntab" est le nombre d'intervalles pour l'integration, sa
c        valeur depend bien sur de l'epaisseur de la couche et de
c        la precision qu'on souhaite a obtenir
c======================================================================
      REAL amb, bmb, sortie
      INTEGER ntab
c======================================================================
      INTEGER n
      REAL xtab(500), xa, xb, ya, yb, xincr
c======================================================================
      external mbtozm
      CHARACTER (LEN=20) :: modname=''
      CHARACTER (LEN=80) :: abort_message
c======================================================================
c la fonction en ligne w(x) donne le profil de l'ozone en fonction
c de l'altitude (unite: cm.atm / km)
c (Green 1964, Appl. Opt. 3: 203-208)
      REAL wp, xp, h, x, w, con
      PARAMETER (wp=0.218, xp=23.25, h=4.63, con=1.0)
      w(x) = wp/h * EXP((x-xp)/h)/ (con+EXP((x-xp)/h))**2
c======================================================================
      IF (ntab .GT. 499) THEN
        abort_message = 'BIG ntab'
        CALL abort_gcm (modname,abort_message,1)
      ENDIF
      xincr = (bmb-amb) / REAL(ntab)
      xtab(1) = amb
      DO n = 2, ntab
         xtab(n) = xtab(n-1) + xincr
      ENDDO
      xtab(ntab+1) = bmb
      sortie = 0.0
      DO n = 1, ntab
         CALL mbtozm(xtab(n), xa)
         CALL mbtozm(xtab(n+1), xb)
         xa = xa / 1000.
         xb = xb / 1000.
         ya = w(xa)
         yb = w(xb)
         sortie = sortie + (ya+yb)/2.0 * ABS(xb-xa)
      ENDDO
      RETURN
      END
      SUBROUTINE mbtozm(rmb,zm)
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS)
c Objet: transformer une hauteur de mb (rmb) en metre (zm)
c======================================================================
      REAL rmb, zm
c======================================================================
      REAL gama, tzero, pzero, g, r
      PARAMETER (gama=6.5e-3, tzero=288., pzero=1013.25)
      PARAMETER (g=9.81, r=287.0)
      zm = tzero/gama * ( 1.-(rmb/pzero)**(r*gama/g) )
      RETURN
      END
