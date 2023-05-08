c
c $Header$
c
       nsrf=3
       DO nreg=1, nbregdyn
       DO i=1, klon

c       IF (debut) THEN
         IF(rlon(i).LT.0.) THEN
           rlonPOS(i)=rlon(i)+360.
         ELSE
           rlonPOS(i)=rlon(i)  
         ENDIF
c       ENDIF

        pct_ocean(i,nreg)=0

c test si c'est 1 point d'ocean
        IF(pctsrf(i,nsrf).EQ.1.) THEN

         IF(nreg.EQ.1) THEN

c TROP
          IF(rlat(i).GE.-30.AND.rlat(i).LE.30.) THEN
           pct_ocean(i,nreg)=1
          ENDIF

c PACIFIQUE NORD
          ELSEIF(nreg.EQ.2) THEN
           IF(rlat(i).GE.40.AND.rlat(i).LE.60.) THEN
            IF(rlonPOS(i).GE.160..AND.rlonPOS(i).LE.235.) THEN 
             pct_ocean(i,nreg)=1
            ENDIF
           ENDIF
c CALIFORNIE ST-CU
         ELSEIF(nreg.EQ.3) THEN
          IF(rlonPOS(i).GE.220..AND.rlonPOS(i).LE.250.) THEN
           IF(rlat(i).GE.15.AND.rlat(i).LE.35.) THEN
            pct_ocean(i,nreg)=1
           ENDIF
          ENDIF
c HAWAI
        ELSEIF(nreg.EQ.4) THEN 
         IF(rlonPOS(i).GE.180..AND.rlonPOS(i).LE.220.) THEN
          IF(rlat(i).GE.15.AND.rlat(i).LE.35.) THEN
           pct_ocean(i,nreg)=1
          ENDIF
         ENDIF
c WARM POOL
        ELSEIF(nreg.EQ.5) THEN 
         IF(rlonPOS(i).GE.70..AND.rlonPOS(i).LE.150.) THEN
          IF(rlat(i).GE.-5.AND.rlat(i).LE.20.) THEN
           pct_ocean(i,nreg)=1
          ENDIF
         ENDIF
        ENDIF !nbregdyn
c TROP
c        IF(rlat(i).GE.-30.AND.rlat(i).LE.30.) THEN
c         pct_ocean(i)=.TRUE.
c         WRITE(*,*) 'pct_ocean =',i, rlon(i), rlat(i)
c          ENDIF !lon
c         ENDIF !lat

        ENDIF !pctsrf
       ENDDO !klon
       ENDDO !nbregdyn
cIM 190504      ENDIF !ok_regdyn
 
cIM somme de toutes les nhistoW BEG
      IF (debut) THEN
      DO nreg = 1, nbregdyn
       DO k = 1, kmaxm1
        DO l = 1, lmaxm1
         DO iw = 1, iwmax
          nhistoWt(k,l,iw,nreg)=0.
         ENDDO !iw
        ENDDO !l
       ENDDO !k
      ENDDO !nreg
      ENDIF !(debut) THEN
cIM 190504 BEG
