










      subroutine aerave( ndata,
     &     longdata,epdata,omegdata,gdata,          
     &     longref,epref,temp,nir,longir
     &     ,epir,omegir,gir,qref )


      implicit none

!==================================================================
!     
!     Purpose
!     -------
!     Calculate mean values of aerosol quantities in each band.
!     
!     Authors
!     ------- 
!     R. Fournier (1996) 
!     F. Forget (1996)
!     
!     Called by
!     ---------
!     suaer_corrk.F90
!     
!     Calls
!     -----
!     blackl.F
!     
!==================================================================

!     
!     R.Fournier 02/1996 
!     (modif F.Forget 02/1996)
!     le spectre est decoupe en "nir" bandes et cette routine calcule
!     les donnees radiatives moyenne sur chaque bande : l'optimisation
!     est faite pour une temperature au sol "temp" et une epaisseur
!     optique de l'atmosphere "epref" a la longueur d'onde "longref"
!     
!     dans la version actuelle, les ponderations sont independantes de
!     l'epaisseur optique : c'est a dire que "omegir", "gir"
!     et "epir/epre" sont independants de "epref".
!     en effet les ponderations sont choisies pour une solution exacte
!     en couche mince et milieu isotherme. 
!     
!     entree
!     
!     ndata : taille des champs data
!     longdata,epdata,omegdata,gdata : proprietes radiative de l'aerosol
!     (longdata longueur d'onde en METRES)
!     * longref : longueur d'onde a laquelle l'epaisseur optique
!     est connue
!     * epref : epaisseur optique a longref
!     * temp : temperature choisie pour la ponderation (Planck)
!     * nir : nombre d'intervals dans la discretisation spectrale
!     du GCM
!     * longir : longueurs d'onde definissant ces intervals
!     
!     sortie
!     
!     * epir : epaisseur optique moyenne pour chaque interval
!     * omegir : "scattering albedo" moyen pour chaque interval
!     * gir : "assymetry factor" moyen pour chaque interval
!     * qref : extinction coefficient at reference wavelength

!=======================================================================
! output

      REAL longref
      REAL epref
      REAL temp
      INTEGER nir
      REAL*8 longir(nir+1)
      REAL epir(nir)
      REAL omegir(nir)
      REAL gir(nir)

!=======================================================================

      INTEGER iir,nirmx
      PARAMETER (nirmx=100)
      INTEGER idata,ndata

      DOUBLE PRECISION tmp1
      REAL tmp2,tmp3

!=======================================================================
! input

      REAL emit
      REAL totalemit(nirmx)
      REAL longdata(ndata),epdata(ndata)
     &     ,omegdata(ndata),gdata(ndata)
      REAL qextcorrdata(ndata)
      INTEGER ibande,nbande
      PARAMETER (nbande=1000)

      REAL long,deltalong
      INTEGER ilong
      INTEGER i1,i2
      REAL c1,c2
      REAL factep,qextcorr,omeg,g,qref

      long=longref


!=======================================================================
!     pre-interpolation
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

      qref=c1*epdata(i1)+c2*epdata(i2)
      factep=qref/epref
      DO idata=1,ndata
         qextcorrdata(idata)=epdata(idata)/factep
      ENDDO
!=======================================================================

      DO iir=1,nir

         deltalong=(longir(iir+1)-longir(iir)) / nbande
         totalemit(iir)=0.E+0
         epir(iir)=0.E+0
         omegir(iir)=0.E+0
         gir(iir)=0.E+0

         DO ibande=1,nbande

            long=longir(iir) + (ibande-0.5E+0) * deltalong
            CALL blackl(DBLE(long),DBLE(temp),tmp1)
            emit=REAL(tmp1)

!=======================================================================
!     interpolation
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
!=======================================================================

            qextcorr=c1*qextcorrdata(i1)+c2*qextcorrdata(i2)
            omeg=c1*omegdata(i1)+c2*omegdata(i2)
            g=c1*gdata(i1)+c2*gdata(i2)

            totalemit(iir)=totalemit(iir)+deltalong*emit
            epir(iir)=epir(iir)+deltalong*emit*qextcorr
            omegir(iir)=omegir(iir)+deltalong*emit*omeg*qextcorr
            gir(iir)=gir(iir)+deltalong*emit*omeg*qextcorr*g

         ENDDO

         gir(iir)=gir(iir)/omegir(iir)
         omegir(iir)=omegir(iir)/epir(iir)
         epir(iir)=epir(iir)/totalemit(iir)

      ENDDO

      return
      end
