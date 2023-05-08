      FUNCTION AZENER(V,J)
      IMPLICIT REAL (A-H,O-Z)
      REAL J,JP
      VP=V+0.5
      JP=J+1.
      E=2359.61*VP-14.456*VP**2+7.51E-3*VP**3+
     $(1.9987-1.87E-2*VP)*J*JP-(5.8E-6-1.E-9*VP)*J**2*JP**2
      AZENER=E
      RETURN
      END
