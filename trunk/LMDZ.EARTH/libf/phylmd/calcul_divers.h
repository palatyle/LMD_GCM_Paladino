c
c $Header$
c

c     Initialisations diverses au "debut" du mois
      IF(debut) THEN
         nday_rain(:)=0.

c        surface terre
         paire_ter(:)=0.
         DO i=1, klon
            IF(pctsrf(i,is_ter).GT.0.) THEN
               paire_ter(i)=airephy(i)*pctsrf(i,is_ter)
            ENDIF 
         ENDDO
      ENDIF

cIM   Calcul une fois par jour : total_rain, nday_rain
      IF(MOD(itap,INT(un_jour/dtime)).EQ.0) THEN
         DO i = 1, klon
            total_rain(i)=rain_fall(i)+snow_fall(i)  
            IF(total_rain(i).GT.0.) nday_rain(i)=nday_rain(i)+1.
         ENDDO
      ENDIF
