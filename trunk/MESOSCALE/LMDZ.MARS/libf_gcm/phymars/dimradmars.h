c-----------------------------------------------------------------------
c------------------------------------------------------------------------
c   INCLUDE 'dimradmars.h'

c   Declaration and initialisation or radiative transfer calculations
c------------------------------------------------------------------------
c------------------------------------------------------------------------

c Splitting of horizontal grid (to reduce program size on workstation)
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c nflev: number of vertical layer
c ndlon: nombre de points de la grille horizontale
c NDLO2 et ndomainsz pour le decoupage de l'appel a la physique
c ATTENTION:  Il faut  1 < ndomainsz =< ngridmx

      INTEGER  NFLEV,NDLON, ndomainsz, NDLO2

       parameter (ndomainsz=ngridmx)
c     parameter (ndomainsz=(ngridmx-1)/20 + 1)
c      parameter (ndomainsz=(ngridmx-1)/5 + 1)

      parameter (NFLEV=nlayermx,NDLON=ndomainsz)	! avec decoupage 
      parameter (NDLO2=NDLON)


c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      REAL longrefir,longrefvis
      REAL long1vis,long2vis,long3vis, long1ir,long2ir
      REAL long1co2,long2co2
      REAL sunfr(2)
      integer nir, nuco2
      INTEGER npademx,nabsmx,nt_pademx, NSUN


c Definition of spectral intervals at thermal infrared wavelengths (LW)
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      parameter (nir=4) ! Total number of thermal IR bands
      parameter (nuco2=2) ! number of bands in CO2 bands
      PARAMETER (long1ir=5.E-6 , long2ir=200.E-6)
      PARAMETER (long1co2=1.E+0 / 865.E+2 , long2co2=1.E+0 / 500.E+2)

c  Warning : the "nir" thermal IR bands are not ordered by wavelength:
c      iir=1 : central 15um CO2 bands \    
c      iir=2 : CO2 band wings    [long1co2-long2co2] MINUS central band
c      iir=3 : 9 um band [long1ir - long1co2]
c      iir=4 : Far IR    [long2co2 - long2ir]
    
c  Definition of spectral interval at solar wavelengths (SW)
c  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      PARAMETER (NSUN=2)   ! do not change that !
c  Boundaries of spectral intervals (m) : 
      PARAMETER (long1vis=0.1E-6 , long2vis=0.5E-6 , long3vis=5.E-6)
c  Fraction of solar energy in solar band #1 [long1vis-long2vis]
      DATA sunfr(1) / 0.274490 /  
c  Fraction of solar energy in solar band #2 [long2vis-long3vis]
      DATA sunfr(2) / 0.725509 /

c Reference wavelengths used to compute reference dust optical depth (m)
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      PARAMETER (longrefir=9.E-6,longrefvis=0.67E-6)

c Number of kind of tracer radiative properties 
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c (ex: naerkind=1 if you use one dust mode without ice ...)
c (ex: naerkind=2 if you use one dust mode and active ice ...)
      integer naerkind        
      parameter (naerkind=1)

c Various initialisation for LW radiative code
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c npademx : nombre de coef de pade 
c nabsmx : ?
c nt_pademx : nombre d'intervalles de temperature pour pade

      PARAMETER (npademx=4,nabsmx=2,nt_pademx=19)
