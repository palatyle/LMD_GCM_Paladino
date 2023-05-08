      SUBROUTINE grille_m(imdep, jmdep, xdata, ydata, entree,
     .                    imar, jmar, x, y, sortie)
c=======================================================================
c z.x.li (le 1 avril 1994) (voir aussi A. Harzallah et L. Fairhead)
c
c Methode naive pour transformer un champ d'une grille fine a une
c grille grossiere. Je considere que les nouveaux points occupent
c une zone adjacente qui comprend un ou plusieurs anciens points
c
c Aucune ponderation est consideree (voir grille_p)
c
c           (c)
c        ----d-----
c        | . . . .|
c        |        |
c     (b)a . * . .b(a)
c        |        |
c        | . . . .|
c        ----c-----
c           (d)
C=======================================================================
c INPUT:
c        imdep, jmdep: dimensions X et Y pour depart
c        xdata, ydata: coordonnees X et Y pour depart
c        entree: champ d'entree a transformer
c OUTPUT:
c        imar, jmar: dimensions X et Y d'arrivee
c        x, y: coordonnees X et Y d'arrivee
c        sortie: champ de sortie deja transforme
C=======================================================================
      IMPLICIT none

      INTEGER imdep, jmdep
      REAL xdata(imdep),ydata(jmdep) 
      REAL entree(imdep,jmdep)
c
      INTEGER imar, jmar
      REAL x(imar),y(jmar)
      REAL sortie(imar,jmar)
c
      INTEGER i, j, ii, jj
      REAL a(2200),b(2200),c(1100),d(1100)
      REAL number(2200,1100)
      REAL distans(2200*1100)
      INTEGER i_proche, j_proche, ij_proche
#ifdef CRAY
      INTEGER ISMIN
#else
      REAL zzmin
#endif
c
      IF (imar.GT.2200 .OR. jmar.GT.1100) THEN
         PRINT*, 'imar ou jmar trop grand', imar, jmar
         CALL ABORT
      ENDIF
c
c Calculer les limites des zones des nouveaux points
c

      a(1) = x(1) - (x(2)-x(1))/2.0
      b(1) = (x(1)+x(2))/2.0
      DO i = 2, imar-1
         a(i) = b(i-1)
         b(i) = (x(i)+x(i+1))/2.0
      ENDDO
      a(imar) = b(imar-1)
      b(imar) = x(imar) + (x(imar)-x(imar-1))/2.0

      c(1) = y(1) - (y(2)-y(1))/2.0
      d(1) = (y(1)+y(2))/2.0
      DO j = 2, jmar-1
         c(j) = d(j-1)
         d(j) = (y(j)+y(j+1))/2.0
      ENDDO
      c(jmar) = d(jmar-1)
      d(jmar) = y(jmar) + (y(jmar)-y(jmar-1))/2.0

      DO i = 1, imar
      DO j = 1, jmar
         number(i,j) = 0.0
         sortie(i,j) = 0.0
      ENDDO
      ENDDO
c
c Determiner la zone sur laquelle chaque ancien point se trouve
c
c
c  .....  Modif  P. Le Van ( 23/08/95 )  ....

      DO ii = 1, imar
      DO jj = 1, jmar
        DO i = 1, imdep
         IF( ( xdata(i)-a(ii).GE.1.e-5.AND.xdata(i)-b(ii).LE.1.e-5 ).OR.
     .     (   xdata(i)-a(ii).LE.1.e-5.AND.xdata(i)-b(ii).GE.1.e-5 )   )
     .           THEN
          DO j = 1, jmdep
          IF( (ydata(j)-c(jj).GE.1.e-5.AND.ydata(j)-d(jj).LE.1.e-5 ).OR.
     .      (  ydata(j)-c(jj).LE.1.e-5.AND.ydata(j)-d(jj).GE.1.e-5 )   )
     .           THEN
               number(ii,jj) = number(ii,jj) + 1.0
               sortie(ii,jj) = sortie(ii,jj) + entree(i,j)
          ENDIF
          ENDDO
         ENDIF
        ENDDO
      ENDDO
      ENDDO
c
c Si aucun ancien point tombe sur une zone, c'est un probleme
c

      DO i = 1, imar
      DO j = 1, jmar
         IF (number(i,j) .GT. 0.001) THEN
         sortie(i,j) = sortie(i,j) / number(i,j)
         ELSE
         PRINT*, 'probleme,i,j=', i,j
ccc         CALL ABORT
         CALL dist_sphe(x(i),y(j),xdata,ydata,imdep,jmdep,distans)
#ifdef CRAY
         ij_proche = ISMIN(imdep*jmdep,distans,1)
#else
         ij_proche = 1
         zzmin = distans(ij_proche)
         DO ii = 2, imdep*jmdep
            IF (distans(ii).LT.zzmin) THEN
               zzmin = distans(ii)
               ij_proche = ii
            ENDIF
         ENDDO
#endif
         j_proche = (ij_proche-1)/imdep + 1
         i_proche = ij_proche - (j_proche-1)*imdep
         PRINT*, "solution:", ij_proche, i_proche, j_proche
         sortie(i,j) = entree(i_proche,j_proche)
         ENDIF
      ENDDO
      ENDDO

      RETURN
      END


      SUBROUTINE grille_p(imdep, jmdep, xdata, ydata, entree,
     .                    imar, jmar, x, y, sortie)
c=======================================================================
c z.x.li (le 1 avril 1994) (voir aussi A. Harzallah et L. Fairhead)
c
c Methode naive pour transformer un champ d'une grille fine a une
c grille grossiere. Je considere que les nouveaux points occupent
c une zone adjacente qui comprend un ou plusieurs anciens points
c
c Consideration de la distance des points (voir grille_m)
c
c           (c)
c        ----d-----
c        | . . . .|
c        |        |
c     (b)a . * . .b(a)
c        |        |
c        | . . . .|
c        ----c-----
c           (d)
C=======================================================================
c INPUT:
c        imdep, jmdep: dimensions X et Y pour depart
c        xdata, ydata: coordonnees X et Y pour depart
c        entree: champ d'entree a transformer
c OUTPUT:
c        imar, jmar: dimensions X et Y d'arrivee
c        x, y: coordonnees X et Y d'arrivee
c        sortie: champ de sortie deja transforme
C=======================================================================
      IMPLICIT none

      INTEGER imdep, jmdep
      REAL xdata(imdep),ydata(jmdep) 
      REAL entree(imdep,jmdep)
c
      INTEGER imar, jmar
      REAL x(imar),y(jmar)
      REAL sortie(imar,jmar)
c
      INTEGER i, j, ii, jj
      REAL a(400),b(400),c(200),d(200)
      REAL number(400,200)
      INTEGER indx(400,200), indy(400,200)
      REAL dist(400,200), distsom(400,200)
c
      IF (imar.GT.400 .OR. jmar.GT.200) THEN
         PRINT*, 'imar ou jmar trop grand', imar, jmar
         CALL ABORT
      ENDIF
c
      IF (imdep.GT.400 .OR. jmdep.GT.200) THEN
         PRINT*, 'imdep ou jmdep trop grand', imdep, jmdep
         CALL ABORT
      ENDIF
c
c calculer les bords a et b de la nouvelle grille
c
      a(1) = x(1) - (x(2)-x(1))/2.0
      b(1) = (x(1)+x(2))/2.0
      DO i = 2, imar-1
         a(i) = b(i-1)
         b(i) = (x(i)+x(i+1))/2.0
      ENDDO
      a(imar) = b(imar-1)
      b(imar) = x(imar) + (x(imar)-x(imar-1))/2.0

c
c calculer les bords c et d de la nouvelle grille
c
      c(1) = y(1) - (y(2)-y(1))/2.0
      d(1) = (y(1)+y(2))/2.0
      DO j = 2, jmar-1
         c(j) = d(j-1)
         d(j) = (y(j)+y(j+1))/2.0
      ENDDO
      c(jmar) = d(jmar-1)
      d(jmar) = y(jmar) + (y(jmar)-y(jmar-1))/2.0

c
c trouver les indices (indx,indy) de la nouvelle grille sur laquelle
c un point de l'ancienne grille est tombe.
c
c
c  .....  Modif  P. Le Van ( 23/08/95 )  ....

      DO ii = 1, imar
      DO jj = 1, jmar
        DO i = 1, imdep
         IF( ( xdata(i)-a(ii).GE.1.e-5.AND.xdata(i)-b(ii).LE.1.e-5 ).OR.
     .     (   xdata(i)-a(ii).LE.1.e-5.AND.xdata(i)-b(ii).GE.1.e-5 )   )
     .           THEN
          DO j = 1, jmdep
          IF( (ydata(j)-c(jj).GE.1.e-5.AND.ydata(j)-d(jj).LE.1.e-5 ).OR.
     .      (  ydata(j)-c(jj).LE.1.e-5.AND.ydata(j)-d(jj).GE.1.e-5 )   )
     .           THEN
               indx(i,j) = ii
               indy(i,j) = jj
          ENDIF
          ENDDO
         ENDIF
        ENDDO
      ENDDO
      ENDDO
c
c faire une verification
c

      DO i = 1, imdep
      DO j = 1, jmdep
         IF (indx(i,j).GT.imar .OR. indy(i,j).GT.jmar) THEN
            PRINT*, 'Probleme grave,i,j,indx,indy=',
     .              i,j,indx(i,j),indy(i,j)
            CALL abort
         ENDIF
      ENDDO
      ENDDO

c
c calculer la distance des anciens points avec le nouveau point,
c on prend ensuite une sorte d'inverse pour ponderation.
c

      DO i = 1, imar
      DO j = 1, jmar
         number(i,j) = 0.0
         distsom(i,j) = 0.0
      ENDDO
      ENDDO
      DO i = 1, imdep
      DO j = 1, jmdep
         dist(i,j) = SQRT ( (xdata(i)-x(indx(i,j)))**2
     .                     +(ydata(j)-y(indy(i,j)))**2 )
         distsom(indx(i,j),indy(i,j)) = distsom(indx(i,j),indy(i,j))
     .                                  + dist(i,j)
         number(indx(i,j),indy(i,j)) = number(indx(i,j),indy(i,j)) +1.
      ENDDO
      ENDDO
      DO i = 1, imdep
      DO j = 1, jmdep
         dist(i,j) = 1.0 - dist(i,j)/distsom(indx(i,j),indy(i,j))
      ENDDO
      ENDDO

      DO i = 1, imar
      DO j = 1, jmar
         number(i,j) = 0.0
         sortie(i,j) = 0.0
      ENDDO
      ENDDO
      DO i = 1, imdep
      DO j = 1, jmdep
         sortie(indx(i,j),indy(i,j)) = sortie(indx(i,j),indy(i,j))
     .                                 + entree(i,j) * dist(i,j)
         number(indx(i,j),indy(i,j)) = number(indx(i,j),indy(i,j))
     .                                 + dist(i,j)
      ENDDO
      ENDDO
      DO i = 1, imar
      DO j = 1, jmar
         IF (number(i,j) .GT. 0.001) THEN
         sortie(i,j) = sortie(i,j) / number(i,j)
         ELSE
         PRINT*, 'probleme,i,j=', i,j
         CALL ABORT
         ENDIF
      ENDDO
      ENDDO

      RETURN
      END




      SUBROUTINE mask_c_o(imdep, jmdep, xdata, ydata, relief,
     .                    imar, jmar, x, y, mask)
c=======================================================================
c z.x.li (le 1 avril 1994): A partir du champ de relief, on fabrique
c                           un champ indicateur (masque) terre/ocean
c                           terre:1; ocean:0
c
c Methode naive (voir grille_m)
C=======================================================================
      IMPLICIT none

      INTEGER imdep, jmdep
      REAL xdata(imdep),ydata(jmdep) 
      REAL relief(imdep,jmdep)
c
      INTEGER imar, jmar
      REAL x(imar),y(jmar)
      REAL mask(imar,jmar)
c
      INTEGER i, j, ii, jj
      REAL a(2200),b(2200),c(1100),d(1100)
      REAL num_tot(2200,1100), num_oce(2200,1100)
c
      IF (imar.GT.2200 .OR. jmar.GT.1100) THEN
         PRINT*, 'imar ou jmar trop grand', imar, jmar
         CALL ABORT
      ENDIF
c

      a(1) = x(1) - (x(2)-x(1))/2.0
      b(1) = (x(1)+x(2))/2.0
      DO i = 2, imar-1
         a(i) = b(i-1)
         b(i) = (x(i)+x(i+1))/2.0
      ENDDO
      a(imar) = b(imar-1)
      b(imar) = x(imar) + (x(imar)-x(imar-1))/2.0

      c(1) = y(1) - (y(2)-y(1))/2.0
      d(1) = (y(1)+y(2))/2.0
      DO j = 2, jmar-1
         c(j) = d(j-1)
         d(j) = (y(j)+y(j+1))/2.0
      ENDDO
      c(jmar) = d(jmar-1)
      d(jmar) = y(jmar) + (y(jmar)-y(jmar-1))/2.0

      DO i = 1, imar
      DO j = 1, jmar
         num_oce(i,j) = 0.0
         num_tot(i,j) = 0.0
      ENDDO
      ENDDO

c
c  .....  Modif  P. Le Van ( 23/08/95 )  ....

      DO ii = 1, imar
      DO jj = 1, jmar
        DO i = 1, imdep
         IF( ( xdata(i)-a(ii).GE.1.e-5.AND.xdata(i)-b(ii).LE.1.e-5 ).OR.
     .     (   xdata(i)-a(ii).LE.1.e-5.AND.xdata(i)-b(ii).GE.1.e-5 )   )
     .           THEN
          DO j = 1, jmdep
          IF( (ydata(j)-c(jj).GE.1.e-5.AND.ydata(j)-d(jj).LE.1.e-5 ).OR.
     .      (  ydata(j)-c(jj).LE.1.e-5.AND.ydata(j)-d(jj).GE.1.e-5 )   )
     .           THEN
               num_tot(ii,jj) = num_tot(ii,jj) + 1.0
               IF (.NOT. ( relief(i,j) - 0.9. GE. 1.e-5 ) )
     .             num_oce(ii,jj) = num_oce(ii,jj) + 1.0
          ENDIF
          ENDDO
         ENDIF
        ENDDO
      ENDDO
      ENDDO
c
c
c
      DO i = 1, imar
      DO j = 1, jmar
         IF (num_tot(i,j) .GT. 0.001) THEN
           IF ( num_oce(i,j)/num_tot(i,j) - 0.5 .GE. 1.e-5 ) THEN
              mask(i,j) = 0.
           ELSE
              mask(i,j) = 1.
           ENDIF
         ELSE
         PRINT*, 'probleme,i,j=', i,j
         CALL ABORT
         ENDIF
      ENDDO
      ENDDO

      RETURN
      END
c
c


      SUBROUTINE rugosite(imdep, jmdep, xdata, ydata, entree,
     .                    imar, jmar, x, y, sortie, mask)
c=======================================================================
c z.x.li (le 1 avril 1994): Transformer la longueur de rugosite d'une
c grille fine a une grille grossiere. Sur l'ocean, on impose une valeur
c fixe (0.001m).
c
c Methode naive (voir grille_m)
C=======================================================================
      IMPLICIT none

      INTEGER imdep, jmdep
      REAL xdata(imdep),ydata(jmdep) 
      REAL entree(imdep,jmdep)
c
      INTEGER imar, jmar
      REAL x(imar),y(jmar)
      REAL sortie(imar,jmar), mask(imar,jmar)
c
      INTEGER i, j, ii, jj
      REAL a(400),b(400),c(400),d(400)
      REAL num_tot(400,400)
      REAL distans(400*400)
      INTEGER i_proche, j_proche, ij_proche
#ifdef CRAY
      INTEGER ISMIN
#else
      REAL zzmin
#endif
c
      IF (imar.GT.400 .OR. jmar.GT.400) THEN
         PRINT*, 'imar ou jmar trop grand', imar, jmar
         CALL ABORT
      ENDIF
c

      a(1) = x(1) - (x(2)-x(1))/2.0
      b(1) = (x(1)+x(2))/2.0
      DO i = 2, imar-1
         a(i) = b(i-1)
         b(i) = (x(i)+x(i+1))/2.0
      ENDDO
      a(imar) = b(imar-1)
      b(imar) = x(imar) + (x(imar)-x(imar-1))/2.0

      c(1) = y(1) - (y(2)-y(1))/2.0
      d(1) = (y(1)+y(2))/2.0
      DO j = 2, jmar-1
         c(j) = d(j-1)
         d(j) = (y(j)+y(j+1))/2.0
      ENDDO
      c(jmar) = d(jmar-1)
      d(jmar) = y(jmar) + (y(jmar)-y(jmar-1))/2.0

      DO i = 1, imar
      DO j = 1, jmar
         num_tot(i,j) = 0.0
         sortie(i,j) = 0.0
      ENDDO
      ENDDO

c
c
c  .....  Modif  P. Le Van ( 23/08/95 )  ....

      DO ii = 1, imar
      DO jj = 1, jmar
        DO i = 1, imdep
         IF( ( xdata(i)-a(ii).GE.1.e-5.AND.xdata(i)-b(ii).LE.1.e-5 ).OR.
     .     (   xdata(i)-a(ii).LE.1.e-5.AND.xdata(i)-b(ii).GE.1.e-5 )   )
     .           THEN
          DO j = 1, jmdep
          IF( (ydata(j)-c(jj).GE.1.e-5.AND.ydata(j)-d(jj).LE.1.e-5 ).OR.
     .      (  ydata(j)-c(jj).LE.1.e-5.AND.ydata(j)-d(jj).GE.1.e-5 )   )
     .           THEN
              sortie(ii,jj)  = sortie(ii,jj) + LOG(entree(i,j))
              num_tot(ii,jj) = num_tot(ii,jj) + 1.0
          ENDIF
          ENDDO
         ENDIF
        ENDDO
      ENDDO
      ENDDO
c

      DO i = 1, imar
      DO j = 1, jmar
       IF (NINT(mask(i,j)).EQ.1) THEN
         IF (num_tot(i,j) .GT. 0.0) THEN
            sortie(i,j) = sortie(i,j) / num_tot(i,j)
            sortie(i,j) = EXP(sortie(i,j))
         ELSE
            PRINT*, 'probleme,i,j=', i,j
ccc            CALL ABORT
         CALL dist_sphe(x(i),y(j),xdata,ydata,imdep,jmdep,distans)
#ifdef CRAY
         ij_proche = ISMIN(imdep*jmdep,distans,1)
#else
         ij_proche = 1
         zzmin = distans(ij_proche)
         DO ii = 2, imdep*jmdep
            IF (distans(ii).LT.zzmin) THEN
               zzmin = distans(ii)
               ij_proche = ii
            ENDIF
         ENDDO
#endif
         j_proche = (ij_proche-1)/imdep + 1
         i_proche = ij_proche - (j_proche-1)*imdep
         PRINT*, "solution:", ij_proche, i_proche, j_proche
         sortie(i,j) = entree(i_proche,j_proche)
         ENDIF
       ELSE
         sortie(i,j) = 0.001
       ENDIF
      ENDDO
      ENDDO

      RETURN
      END





      SUBROUTINE sea_ice(imdep, jmdep, xdata, ydata, glace01,
     .                    imar, jmar, x, y, frac_ice)
c=======================================================================
c z.x.li (le 1 avril 1994): Transformer un champ d'indicateur de la
c glace (1, sinon 0) d'une grille fine a un champ de fraction de glace
c (entre 0 et 1) dans une grille plus grossiere.
c
c Methode naive (voir grille_m)
C=======================================================================
      IMPLICIT none

      INTEGER imdep, jmdep
      REAL xdata(imdep),ydata(jmdep) 
      REAL glace01(imdep,jmdep)
c
      INTEGER imar, jmar
      REAL x(imar),y(jmar)
      REAL frac_ice(imar,jmar)
c
      INTEGER i, j, ii, jj
      REAL a(400),b(400),c(400),d(400)
      REAL num_tot(400,400), num_ice(400,400)
      REAL distans(400*400)
      INTEGER i_proche, j_proche, ij_proche
#ifdef CRAY
      INTEGER ISMIN
#else
      REAL zzmin
#endif
c
      IF (imar.GT.400 .OR. jmar.GT.400) THEN
         PRINT*, 'imar ou jmar trop grand', imar, jmar
         CALL ABORT
      ENDIF
c

      a(1) = x(1) - (x(2)-x(1))/2.0
      b(1) = (x(1)+x(2))/2.0
      DO i = 2, imar-1
         a(i) = b(i-1)
         b(i) = (x(i)+x(i+1))/2.0
      ENDDO
      a(imar) = b(imar-1)
      b(imar) = x(imar) + (x(imar)-x(imar-1))/2.0

      c(1) = y(1) - (y(2)-y(1))/2.0
      d(1) = (y(1)+y(2))/2.0
      DO j = 2, jmar-1
         c(j) = d(j-1)
         d(j) = (y(j)+y(j+1))/2.0
      ENDDO
      c(jmar) = d(jmar-1)
      d(jmar) = y(jmar) + (y(jmar)-y(jmar-1))/2.0

      DO i = 1, imar
      DO j = 1, jmar
         num_ice(i,j) = 0.0
         num_tot(i,j) = 0.0
      ENDDO
      ENDDO

c
c
c  .....  Modif  P. Le Van ( 23/08/95 )  ....

      DO ii = 1, imar
      DO jj = 1, jmar
        DO i = 1, imdep
         IF( ( xdata(i)-a(ii).GE.1.e-5.AND.xdata(i)-b(ii).LE.1.e-5 ).OR.
     .     (   xdata(i)-a(ii).LE.1.e-5.AND.xdata(i)-b(ii).GE.1.e-5 )   )
     .           THEN
          DO j = 1, jmdep
          IF( (ydata(j)-c(jj).GE.1.e-5.AND.ydata(j)-d(jj).LE.1.e-5 ).OR.
     .      (  ydata(j)-c(jj).LE.1.e-5.AND.ydata(j)-d(jj).GE.1.e-5 )   )
     .           THEN
             num_tot(ii,jj) = num_tot(ii,jj) + 1.0
              IF (NINT(glace01(i,j)).EQ.1 ) 
     .       num_ice(ii,jj) = num_ice(ii,jj) + 1.0
          ENDIF
          ENDDO
         ENDIF
        ENDDO
      ENDDO
      ENDDO
c
c

      DO i = 1, imar
      DO j = 1, jmar
         IF (num_tot(i,j) .GT. 0.001) THEN
           IF (num_ice(i,j).GT.0.001) THEN
            frac_ice(i,j) = num_ice(i,j) / num_tot(i,j)
           ELSE
              frac_ice(i,j) = 0.0
           ENDIF
         ELSE
           PRINT*, 'probleme,i,j=', i,j
ccc           CALL ABORT
         CALL dist_sphe(x(i),y(j),xdata,ydata,imdep,jmdep,distans)
#ifdef CRAY
         ij_proche = ISMIN(imdep*jmdep,distans,1)
#else
         ij_proche = 1
         zzmin = distans(ij_proche)
         DO ii = 2, imdep*jmdep
            IF (distans(ii).LT.zzmin) THEN
               zzmin = distans(ii)
               ij_proche = ii
            ENDIF
         ENDDO
#endif
         j_proche = (ij_proche-1)/imdep + 1
         i_proche = ij_proche - (j_proche-1)*imdep
         PRINT*, "solution:", ij_proche, i_proche, j_proche
         IF (NINT(glace01(i_proche,j_proche)).EQ.1 ) THEN
            frac_ice(i,j) = 1.0
         ELSE
            frac_ice(i,j) = 0.0
         ENDIF
         ENDIF
      ENDDO
      ENDDO

      RETURN
      END



      SUBROUTINE rugsoro(imrel, jmrel, xrel, yrel, relief,
     .                    immod, jmmod, xmod, ymod, rugs)
c=======================================================================
c Calculer la longueur de rugosite liee au relief en utilisant
c l'ecart-type dans une maille de 1x1
C=======================================================================
      IMPLICIT none
c
#ifdef CRAY
      INTEGER ISMIN
#else
      REAL zzmin
#endif
c
      REAL amin, AMAX
c
      INTEGER imrel, jmrel
      REAL xrel(imrel),yrel(jmrel)
      REAL relief(imrel,jmrel)
c
      INTEGER immod, jmmod
      REAL xmod(immod),ymod(jmmod)
      REAL rugs(immod,jmmod)
c
      INTEGER imtmp, jmtmp
      PARAMETER (imtmp=360,jmtmp=180)
      REAL xtmp(imtmp), ytmp(jmtmp)
      REAL*8 cham1tmp(imtmp,jmtmp), cham2tmp(imtmp,jmtmp)
      REAL zzzz
c
      INTEGER i, j, ii, jj
      REAL a(2200),b(2200),c(1100),d(1100)
      REAL number(2200,1100)
c
      REAL distans(400*400)
      INTEGER i_proche, j_proche, ij_proche
c
      IF (immod.GT.2200 .OR. jmmod.GT.1100) THEN
         PRINT*, 'immod ou jmmod trop grand', immod, jmmod
         CALL ABORT
      ENDIF
c
c Calculs intermediares:
c
      xtmp(1) = -180.0 + 360.0/FLOAT(imtmp) / 2.0
      DO i = 2, imtmp
         xtmp(i) = xtmp(i-1) + 360.0/FLOAT(imtmp)
      ENDDO
      DO i = 1, imtmp
         xtmp(i) = xtmp(i) /180.0 * 4.0*ATAN(1.0)
      ENDDO
      ytmp(1) = -90.0 + 180.0/FLOAT(jmtmp) / 2.0
      DO j = 2, jmtmp
         ytmp(j) = ytmp(j-1) + 180.0/FLOAT(jmtmp)
      ENDDO
      DO j = 1, jmtmp
         ytmp(j) = ytmp(j) /180.0 * 4.0*ATAN(1.0)
      ENDDO
c
      a(1) = xtmp(1) - (xtmp(2)-xtmp(1))/2.0
      b(1) = (xtmp(1)+xtmp(2))/2.0
      DO i = 2, imtmp-1
         a(i) = b(i-1)
         b(i) = (xtmp(i)+xtmp(i+1))/2.0
      ENDDO
      a(imtmp) = b(imtmp-1)
      b(imtmp) = xtmp(imtmp) + (xtmp(imtmp)-xtmp(imtmp-1))/2.0

      c(1) = ytmp(1) - (ytmp(2)-ytmp(1))/2.0
      d(1) = (ytmp(1)+ytmp(2))/2.0
      DO j = 2, jmtmp-1
         c(j) = d(j-1)
         d(j) = (ytmp(j)+ytmp(j+1))/2.0
      ENDDO
      c(jmtmp) = d(jmtmp-1)
      d(jmtmp) = ytmp(jmtmp) + (ytmp(jmtmp)-ytmp(jmtmp-1))/2.0

      DO i = 1, imtmp
      DO j = 1, jmtmp
         number(i,j) = 0.0
         cham1tmp(i,j) = 0.0
         cham2tmp(i,j) = 0.0
      ENDDO
      ENDDO
c
c
c
c  .....  Modif  P. Le Van ( 23/08/95 )  ....

      DO ii = 1, imtmp
      DO jj = 1, jmtmp
        DO i = 1, imrel
         IF( ( xrel(i)-a(ii).GE.1.e-5.AND.xrel(i)-b(ii).LE.1.e-5 ).OR.
     .     (   xrel(i)-a(ii).LE.1.e-5.AND.xrel(i)-b(ii).GE.1.e-5 )   )
     .           THEN
          DO j = 1, jmrel
          IF( (yrel(j)-c(jj).GE.1.e-5.AND.yrel(j)-d(jj).LE.1.e-5 ).OR.
     .      (  yrel(j)-c(jj).LE.1.e-5.AND.yrel(j)-d(jj).GE.1.e-5 )   )
     .           THEN
              number(ii,jj) = number(ii,jj) + 1.0
              cham1tmp(ii,jj) = cham1tmp(ii,jj) + relief(i,j)
              cham2tmp(ii,jj) = cham2tmp(ii,jj) 
     .                              + relief(i,j)*relief(i,j)
          ENDIF
          ENDDO
         ENDIF
        ENDDO
      ENDDO
      ENDDO
c
c
      DO i = 1, imtmp
      DO j = 1, jmtmp
         IF (number(i,j) .GT. 0.001) THEN
         cham1tmp(i,j) = cham1tmp(i,j) / number(i,j)
         cham2tmp(i,j) = cham2tmp(i,j) / number(i,j)
         zzzz=cham2tmp(i,j)-cham1tmp(i,j)**2
         if (zzzz .lt. 0.0) then
           if (zzzz .gt. -7.5) then
             zzzz = 0.0
             print*,'Pb rugsoro, -7.5 < zzzz < 0, => zzz = 0.0'
           else
              stop 'Pb rugsoro, zzzz <-7.5'
           endif
         endif
         cham2tmp(i,j) = SQRT(zzzz)
         ELSE
         PRINT*, 'probleme,i,j=', i,j
         CALL ABORT
         ENDIF
      ENDDO
      ENDDO
c
      amin = cham2tmp(1,1)
      AMAX = cham2tmp(1,1)
      DO j = 1, jmtmp
      DO i = 1, imtmp
         IF (cham2tmp(i,j).GT.AMAX) AMAX = cham2tmp(i,j)
         IF (cham2tmp(i,j).LT.amin) amin = cham2tmp(i,j)
      ENDDO
      ENDDO
      PRINT*, 'Ecart-type 1x1:', amin, AMAX
c
c
c
      a(1) = xmod(1) - (xmod(2)-xmod(1))/2.0
      b(1) = (xmod(1)+xmod(2))/2.0
      DO i = 2, immod-1
         a(i) = b(i-1)
         b(i) = (xmod(i)+xmod(i+1))/2.0
      ENDDO
      a(immod) = b(immod-1)
      b(immod) = xmod(immod) + (xmod(immod)-xmod(immod-1))/2.0

      c(1) = ymod(1) - (ymod(2)-ymod(1))/2.0
      d(1) = (ymod(1)+ymod(2))/2.0
      DO j = 2, jmmod-1
         c(j) = d(j-1)
         d(j) = (ymod(j)+ymod(j+1))/2.0
      ENDDO
      c(jmmod) = d(jmmod-1)
      d(jmmod) = ymod(jmmod) + (ymod(jmmod)-ymod(jmmod-1))/2.0
c
      DO i = 1, immod
      DO j = 1, jmmod
         number(i,j) = 0.0
         rugs(i,j) = 0.0
      ENDDO
      ENDDO
c
c
c
c  .....  Modif  P. Le Van ( 23/08/95 )  ....

      DO ii = 1, immod
      DO jj = 1, jmmod
        DO i = 1, imtmp
         IF( ( xtmp(i)-a(ii).GE.1.e-5.AND.xtmp(i)-b(ii).LE.1.e-5 ).OR.
     .     (   xtmp(i)-a(ii).LE.1.e-5.AND.xtmp(i)-b(ii).GE.1.e-5 )   )
     .           THEN
          DO j = 1, jmtmp
          IF( (ytmp(j)-c(jj).GE.1.e-5.AND.ytmp(j)-d(jj).LE.1.e-5 ).OR.
     .      (  ytmp(j)-c(jj).LE.1.e-5.AND.ytmp(j)-d(jj).GE.1.e-5 )   )
     .           THEN
              number(ii,jj) = number(ii,jj) + 1.0
              rugs(ii,jj) = rugs(ii,jj)
     .                       + LOG(MAX(0.001,real(cham2tmp(i,j))))
          ENDIF
          ENDDO
         ENDIF
        ENDDO
      ENDDO
      ENDDO
c
c
      DO i = 1, immod
      DO j = 1, jmmod
         IF (number(i,j) .GT. 0.001) THEN
         rugs(i,j) = rugs(i,j) / number(i,j)
         rugs(i,j) = EXP(rugs(i,j))
         ELSE
         PRINT*, 'probleme,i,j=', i,j
ccc         CALL ABORT
         CALL dist_sphe(xmod(i),ymod(j),xtmp,ytmp,imtmp,jmtmp,distans)
#ifdef CRAY
         ij_proche = ISMIN(imtmp*jmtmp,distans,1)
#else
         ij_proche = 1
         zzmin = distans(ij_proche)
         DO ii = 2, imtmp*jmtmp
            IF (distans(ii).LT.zzmin) THEN
               zzmin = distans(ii)
               ij_proche = ii
            ENDIF
         ENDDO
#endif
         j_proche = (ij_proche-1)/imtmp + 1
         i_proche = ij_proche - (j_proche-1)*imtmp
         PRINT*, "solution:", ij_proche, i_proche, j_proche
         rugs(i,j) = LOG(MAX(0.001,real(cham2tmp(i_proche,j_proche))))
         ENDIF
      ENDDO
      ENDDO
c
      amin = rugs(1,1)
      AMAX = rugs(1,1)
      DO j = 1, jmmod
      DO i = 1, immod
         IF (rugs(i,j).GT.AMAX) AMAX = rugs(i,j)
         IF (rugs(i,j).LT.amin) amin = rugs(i,j)
      ENDDO
      ENDDO
      PRINT*, 'Ecart-type du modele:', amin, AMAX
c
      DO j = 1, jmmod
      DO i = 1, immod
         rugs(i,j) = rugs(i,j) / AMAX * 20.0
      ENDDO
      ENDDO
c
      amin = rugs(1,1)
      AMAX = rugs(1,1)
      DO j = 1, jmmod
      DO i = 1, immod
         IF (rugs(i,j).GT.AMAX) AMAX = rugs(i,j)
         IF (rugs(i,j).LT.amin) amin = rugs(i,j)
      ENDDO
      ENDDO
      PRINT*, 'Longueur de rugosite du modele:', amin, AMAX
c
      RETURN
      END
c
      SUBROUTINE dist_sphe(rf_lon,rf_lat,rlon,rlat,im,jm,distance)
c
c Auteur: Laurent Li (le 30 decembre 1996)
c
c Ce programme calcule la distance minimale (selon le grand cercle)
c entre deux points sur la terre
c
c Input:
      INTEGER im, jm ! dimensions
      REAL rf_lon ! longitude du point de reference (degres)
      REAL rf_lat ! latitude du point de reference (degres)
      REAL rlon(im), rlat(jm) ! longitude et latitude des points
c
c Output:
      REAL distance(im,jm) ! distances en metre
c
      REAL rlon1, rlat1
      REAL rlon2, rlat2
      REAL dist
      REAL pa, pb, p, pi
c
      REAL radius
      PARAMETER (radius=6371229.)
c
      pi = 4.0 * ATAN(1.0)
c
      DO 9999 j = 1, jm
      DO 9999 i = 1, im
c
      rlon1=rf_lon
      rlat1=rf_lat
      rlon2=rlon(i)
      rlat2=rlat(j)
      pa = pi/2.0 - rlat1*pi/180.0 ! dist. entre pole n et point a
      pb = pi/2.0 - rlat2*pi/180.0 ! dist. entre pole n et point b
      p = (rlon1-rlon2)*pi/180.0 ! angle entre a et b (leurs meridiens)
c
      dist = ACOS( COS(pa)*COS(pb) + SIN(pa)*SIN(pb)*COS(p))
      dist = radius * dist
      distance(i,j) = dist
c
 9999 CONTINUE
c
      END
