      SUBROUTINE aerave ( ndata,
     & longdata,epdata,omegdata,gdata,          
     &            longref,epref,temp,nir,longir
     &            ,epir,omegir,gir,qref,omegaref        )
      IMPLICIT NONE
c.......................................................................
c
c R.Fournier 02/1996 
c (modif F.Forget 02/1996)
c le spectre est decoupe en "nir" bandes et cette routine calcule
c les donnees radiatives moyenne sur chaque bande : l'optimisation
c est faite pour une temperature au sol "temp" et une epaisseur
c optique de l'atmosphere "epref" a la longueur d'onde "longref"
c
c dans la version actuelle, les ponderations sont independantes de
c l'epaisseur optique : c'est a dire que "omegir", "gir"
c et "epir/epre" sont independants de "epref".
c en effet les ponderations sont choisies pour une solution exacte
c en couche mince et milieu isotherme. 
c
c entree
c
c    ndata : taille des champs data
c    longdata,epdata,omegdata,gdata : proprietes radiative de l'aerosol
c                  (longdata longueur d'onde en METRES)
c  * longref : longueur d'onde a laquelle l'epaisseur optique
c              est connue
c  * epref : epaisseur optique a longref
c  * temp : temperature choisie pour la ponderation (Planck)
c  * nir : nombre d'intervals dans la discretisation spectrale
c           du GCM
c  * longir : longueurs d'onde definissant ces intervals
c
c sortie
c
c  * epir : epaisseur optique moyenne pour chaque interval
c  * omegir : "scattering albedo" moyen pour chaque interval
c  * gir : "assymetry factor" moyen pour chaque interval
c  * qref : extinction coefficient at reference wavelength
c  * omegaref : single scat. albedo at reference wavelength
c
c.......................................................................
c
      REAL longref
      REAL epref
      REAL temp
      INTEGER nir
      REAL longir(nir+1)
      REAL epir(nir)
      REAL omegir(nir)
      REAL gir(nir)
c
c.......................................................................
c
      INTEGER iir,nirmx
      PARAMETER (nirmx=100)
      INTEGER idata,ndata
c
c.......................................................................
c
      REAL emit
      REAL totalemit(nirmx)
      REAL longdata(ndata),epdata(ndata)
     &    ,omegdata(ndata),gdata(ndata)
      REAL qextcorrdata(ndata)
      INTEGER ibande,nbande
      PARAMETER (nbande=1000)
      REAL long,deltalong
      INTEGER ilong
      INTEGER i1,i2
      REAL c1,c2
      REAL factep,qextcorr,omeg,g
      REAL qref,omegaref
c
c.......................................................................
c
      DOUBLE PRECISION tmp1
      REAL tmp2,tmp3
c
c
      long=longref
c
c********************************************************
c interpolation
      ilong=1
      DO idata=2,ndata
        IF (long.gt.longdata(idata)) ilong=idata
      ENDDO
      i1=ilong
      i2=ilong+1
      IF (i2.gt.ndata) i2=ndata
      IF (long.lt.longdata(1)) i2=1
      IF (i1.eq.i2) THEN
        c1=1.E+0
        c2=0.E+0
      ELSE
        c1=(longdata(i2)-long) / (longdata(i2)-longdata(i1))
        c2=(longdata(i1)-long) / (longdata(i1)-longdata(i2))
      ENDIF
c********************************************************
c
      qref=c1*epdata(i1)+c2*epdata(i2)
      omegaref=c1*omegdata(i1)+c2*omegdata(i2)
      factep=qref/epref
      DO idata=1,ndata
        qextcorrdata(idata)=epdata(idata)/factep
      ENDDO
c
c.......................................................................
c
      DO iir=1,nir
c
c.......................................................................
c
        deltalong=(longir(iir+1)-longir(iir)) / nbande
        totalemit(iir)=0.E+0
        epir(iir)=0.E+0
        omegir(iir)=0.E+0
        gir(iir)=0.E+0
c
c.......................................................................
c
        DO ibande=1,nbande
c
c.......................................................................
c
          long=longir(iir) + (ibande-0.5E+0) * deltalong
          CALL blackl(DBLE(long),DBLE(temp),tmp1)
          emit=REAL(tmp1)
c
c.......................................................................
c
c********************************************************
c interpolation
      ilong=1
      DO idata=2,ndata
        IF (long.gt.longdata(idata)) ilong=idata
      ENDDO
      i1=ilong
      i2=ilong+1
      IF (i2.gt.ndata) i2=ndata
      IF (long.lt.longdata(1)) i2=1
      IF (i1.eq.i2) THEN
        c1=1.E+0
        c2=0.E+0
      ELSE
        c1=(longdata(i2)-long) / (longdata(i2)-longdata(i1))
        c2=(longdata(i1)-long) / (longdata(i1)-longdata(i2))
      ENDIF
c********************************************************
c
          qextcorr=c1*qextcorrdata(i1)+c2*qextcorrdata(i2)
          omeg=c1*omegdata(i1)+c2*omegdata(i2)
          g=c1*gdata(i1)+c2*gdata(i2)
c
c.......................................................................
c
          totalemit(iir)=totalemit(iir)+deltalong*emit
          epir(iir)=epir(iir)+deltalong*emit*qextcorr
          omegir(iir)=omegir(iir)+deltalong*emit*omeg*qextcorr
          gir(iir)=gir(iir)+deltalong*emit*omeg*qextcorr*g
c
c.......................................................................
c
        ENDDO
c
c.......................................................................
c
        gir(iir)=gir(iir)/omegir(iir)
        omegir(iir)=omegir(iir)/epir(iir)
        epir(iir)=epir(iir)/totalemit(iir)
c
c.......................................................................
c
      ENDDO
c
c......................................................................
c
c     Diagnostic de controle si on moyenne sur tout le spectre vis ou IR :
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     tmp2=0.E+0
c     DO iir=1,nir
c       tmp2=tmp2+totalemit(iir)
c     ENDDO
c     tmp3=5.67E-8 * temp**4
c     IF (abs((tmp2-tmp3)/tmp3).gt.0.05E+0) THEN
c       PRINT *,'!!!! <---> il manque du Planck (voir moyenne.F)'
c       PRINT *,'somme des bandes :',tmp2,'--- Planck:',tmp3
c     ENDIF
c
c......................................................................
c
      RETURN
      END
