c=======================================================================
      SUBROUTINE datareadnc(relief,phisinit,alb,ith,
     .                    zmea,zstd,zsig,zgam,zthe)
c=======================================================================
c
c
c   Author: F. Hourdin      01/1997
c   -------
c
c   Object: To read data from Martian surface to use in a GCM
c   ------                from NetCDF file "surface.nc"
c
c
c   Arguments:
c   ----------
c
c     Inputs:
c     ------
c
c     Outputs:
c     --------
c
c=======================================================================
c   donnees ALBEDO, INERTIE THERMIQUE, RELIEF:
c
c       Ces donnees sont au format NetCDF dans le fichier "surface.nc"
c
c   360 valeurs en longitude (de -179.5 a 179.5)
c   180 valeurs en latitudes (de 89.5 a -89.5)
c
c   Pour les passer au format de la grille, on utilise "interp_horiz.F"
c
c   Il faut donc que ces donnees soient au format grille scalaire
c               (imold+1 jmold+1)
c       avec passage des coordonnees de la "boite" (rlonu, rlatv)
c
c   On prend imd (d pour donnees!) 
c           imd = 360 avec copie de la 1ere valeur sur la imd+1 
c                   (rlonud de -179 a -181)
c           jmd = 179 
c                   (rlatvd de 89 a -89)
c=======================================================================

      implicit none

#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
#include "comconst.h"
#include "netcdf.inc"
#include "datafile.h"

c=======================================================================
c   Declarations:
C=======================================================================

      INTEGER    imd,jmd,imdp1,jmdp1
      parameter    (imd=360,jmd=179,imdp1=361,jmdp1=180)

      INTEGER    iimp1
      parameter    (iimp1=iim+1-1/iim)

      CHARACTER relief*3

      REAL        zdata(imd*jmdp1)
      REAL        zdataS(imdp1*jmdp1)
      REAL        pfield(iimp1*jjp1)

      REAL        alb(iimp1*jjp1)
      REAL        ith(iimp1*jjp1)
      REAL        phisinit(iimp1*jjp1)

      REAL        zmea(imdp1*jmdp1)
      REAL        zstd(imdp1*jmdp1)
      REAL        zsig(imdp1*jmdp1)
      REAL        zgam(imdp1*jmdp1)
      REAL        zthe(imdp1*jmdp1)

      INTEGER   lnblnk, ierr
      EXTERNAL    lnblnk

      INTEGER   unit,nvarid

      INTEGER    i,j,k

      INTEGER klatdat,ngridmxgdat
      PARAMETER (klatdat=180,ngridmxgdat=360)

c    on passe une grille en rlonu rlatv et im+1 jm a interp_horiz)

      REAL longitude(imd),latitude(jmdp1) ! Pour lecture des donnees
      REAL rlonud(imdp1),rlatvd(jmd)

      CHARACTER*20 string
      DIMENSION string(4)

      CHARACTER*80    file

#include "lmdstd.h"
#include "fxyprim.h"

      pi=2.*ASIN(1.)

c=======================================================================
c    rlonud, rlatvd
c=======================================================================

c-----------------------------------------------------------------------
c    Lecture NetCDF des donnees latitude et longitude
c-----------------------------------------------------------------------
      write(*,*) 'ouverture du fichier surface.nc'

      ierr = NF_OPEN (datafile(1:lnblnk(datafile))//'/surface.nc',
     &  NF_NOWRITE,unit)
      IF (ierr.NE.NF_NOERR) THEN
        write(*,*)'Error : cannot open file surface.nc '
        write(*,*)'(in phymars/datareadnc.F)'
        write(*,*)'It should be in :',datafile(1:lnblnk(datafile)),'/'
        write(*,*)'1) You can change this directory address in '
        write(*,*)'   file phymars/datafile.h'
        write(*,*)'2) If necessary surface.nc (and other datafiles)'
        write(*,*)'   can be obtained online on:'
        write(*,*)' http://www.lmd.jussieu.fr/~forget/datagcm/datafile'
        CALL ABORT
      ENDIF

c
c Lecture des latitudes (coordonnees):
c
      ierr = NF_INQ_VARID (unit, "latitude", nvarid)
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(unit, nvarid, latitude)
#else
      ierr = NF_GET_VAR_REAL(unit, nvarid, latitude)
#endif
c
c Lecture des longitudes (coordonnees):
c
      ierr = NF_INQ_VARID (unit, "longitude", nvarid)
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(unit, nvarid, longitude)
#else
      ierr = NF_GET_VAR_REAL(unit, nvarid, longitude)
#endif

c-----------------------------------------------------------------------
c    Passage au format boites scalaires
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c       longitude(imd)        -->      rlonud(imdp1) 
c-----------------------------------------------------------------------

c Passage en coordonnees boites scalaires et en radian
      do i=1,imd 
          rlonud(i)=(longitude(i)+.5)*pi/180.
      enddo

c Repetition de la valeur im+1
      rlonud(imdp1)=rlonud(1) + 2*pi

c-----------------------------------------------------------------------
c        latitude(jmdp1)         -->        rlonvd(jmd)
c-----------------------------------------------------------------------

c Passage en coordonnees boites scalaires et en radian
      do j=1,jmd 
          rlatvd(j)=(latitude(j)-.5)*pi/180.
      enddo

c=======================================================================
c   lecture NetCDF de albedo, thermal, relief, zdtm (pour francois Lott)
c=======================================================================

      string(1) = 'albedo'
      string(2) = 'thermal'
      if (relief.ne.'pla') then
        write(*,*) ' La topographie est celle de MOLA'
        relief = 'MOL'
          string(3) = 'z'//relief
      else
          string(3) = 'zMOL'  ! pour qu''il lise qqchose sur le fichier
                            ! remise a 0 derriere
      endif
      string(4) = 'zMOL'    ! lecture pour calcul topog. sous-maille
 

      DO k=1,4
          write(*,*) 'string',k,string(k)
          
c-----------------------------------------------------------------------
c    initialisation
c-----------------------------------------------------------------------
      call initial0(iimp1*jjp1,pfield)
      call initial0(imd*jmdp1,zdata)
      call initial0(imdp1*jmdp1,zdataS)

c-----------------------------------------------------------------------
c    Lecture NetCDF  
c-----------------------------------------------------------------------

      ierr = NF_INQ_VARID (unit, string(k), nvarid)
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(unit, nvarid, zdata)
#else
      ierr = NF_GET_VAR_REAL(unit, nvarid, zdata)
#endif

c-----------------------------------------------------------------------
c        Cas particulier "Francois Lott" ( k=4 ) (relief sous-maille)
c-----------------------------------------------------------------------
      if (k.eq.4) then

          call multscal(imd*jmdp1,zdata,1000.,zdata)
          call multscal(imd,longitude,pi/180.,longitude)
          call multscal(jmdp1,latitude,pi/180.,latitude)

          call grid_noro1(360, 180, longitude, latitude, zdata,
     .         iim, jjp1, rlonv, rlatu, zmea,zstd,zsig,zgam,zthe)

          CALL dump2d(iip1,jjp1,zmea,'zmea')
          CALL dump2d(iip1,jjp1,zstd,'zstd')
          CALL dump2d(iip1,jjp1,zsig,'zsig')
          CALL dump2d(iip1,jjp1,zgam,'zgam')
          CALL dump2d(iip1,jjp1,zthe,'zthe')

      endif

c-----------------------------------------------------------------------
c   Passage de zdata en grille (imdp1 jmdp1)
c-----------------------------------------------------------------------
      do j=1,jmdp1
          do i=1,imd
              zdataS(i+imdp1*(j-1)) = zdata(i+ngridmxgdat*(j-1))
          enddo
          zdataS(imdp1+imdp1*(j-1)) = zdata(1+ngridmxgdat*(j-1))
      enddo

c-----------------------------------------------------------------------
c    Interpolation
c-----------------------------------------------------------------------
      call interp_horiz(zdataS,pfield,imd,jmd,
     .    iim, jjm,1,rlonud,rlatvd,rlonu,rlatv) 

c-----------------------------------------------------------------------
c    Periodicite    
c-----------------------------------------------------------------------

      do j=1,jjp1
         pfield(iimp1*j) =  pfield(1+iimp1*(j-1))
      enddo 
 
c-----------------------------------------------------------------------
c    Sauvegarde des champs    
c-----------------------------------------------------------------------

      if (k.eq.1) then                    ! albedo
         do i=1,iimp1*jjp1
              alb(i) = pfield(i)
          enddo
      elseif (k.eq.2) then                ! thermal
         do i=1,iimp1*jjp1
              ith(i) = pfield(i)
          enddo
      elseif (k.eq.3) then                ! relief
        if (relief.eq.'pla') then
              call initial0(iimp1*jjp1,phisinit)
        else
             do i=1,iimp1*jjp1
                  phisinit(i) = pfield(i)
              enddo
        endif
      endif

      ENDDO

c-----------------------------------------------------------------------
c    Traitement Phisinit
c-----------------------------------------------------------------------

      DO i=1,iimp1*jjp1
            phisinit(i)=1000.*phisinit(i)
      ENDDO
      CALL dump2d(iimp1,jjp1,phisinit,'Altitude en m')
      CALL multscal(iimp1*jjp1,phisinit,g,phisinit)

c-----------------------------------------------------------------------
c    FIN
c-----------------------------------------------------------------------

      END
