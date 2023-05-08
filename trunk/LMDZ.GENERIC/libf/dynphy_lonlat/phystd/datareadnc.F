c=======================================================================
      SUBROUTINE datareadnc(relief,filename,phisinit,alb,ith,
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

      use datafile_mod, only: datadir, surfdir
! to use  'getin'
      USE ioipsl_getincom 
      USE comconst_mod, ONLY: g,pi
      implicit none

#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
#include "netcdf.inc"

c=======================================================================
c   Declarations:
C=======================================================================

      INTEGER    imd,jmd,imdp1,jmdp1
      parameter    (imd=360,jmd=179,imdp1=361,jmdp1=180)

      INTEGER    iimp1
      parameter    (iimp1=iim+1-1/iim)

      character(len=3),intent(inout) :: relief*3
      character(len=*),intent(in) :: filename ! surface.nc file
      real,intent(out) :: phisinit(iimp1*jjp1) ! surface geopotential
      real,intent(out) :: alb(iimp1*jjp1) ! albedo
      real,intent(out) :: ith(iimp1*jjp1) ! thermal inertia
      real,intent(out) :: zmea(imdp1*jmdp1)
      real,intent(out) :: zstd(imdp1*jmdp1)
      real,intent(out) :: zsig(imdp1*jmdp1)
      real,intent(out) :: zgam(imdp1*jmdp1)
      real,intent(out) :: zthe(imdp1*jmdp1)
      
      REAL        zdata(imd*jmdp1)
      REAL        zdataS(imdp1*jmdp1)
      REAL        pfield(iimp1*jjp1)

      INTEGER   ierr

      INTEGER   unit,nvarid

      INTEGER    i,j,k

      INTEGER klatdat,ngridmixgdat
      PARAMETER (klatdat=180,ngridmixgdat=360)

c    on passe une grille en rlonu rlatv et im+1 jm a interp_horiz)

      REAL longitude(imd),latitude(jmdp1) ! Pour lecture des donnees
      REAL rlonud(imdp1),rlatvd(jmd)

      CHARACTER*20 string
      DIMENSION string(4)
!#include "fxyprim.h"

      pi=2.*ASIN(1.)

c=======================================================================
c    rlonud, rlatvd
c=======================================================================

c-----------------------------------------------------------------------
c    Lecture NetCDF des donnees latitude et longitude
c-----------------------------------------------------------------------
      ierr = NF_OPEN (trim(datadir)//'/'//trim(surfdir)//'/'//
     &                trim(adjustl(filename)),
     &                NF_NOWRITE,unit)
      IF (ierr.NE.NF_NOERR) THEN
        ! In ye old days this file was stored in datadir;
        ! let's be retro-compatible
        ierr = NF_OPEN (trim(datadir)//'/'//
     &                trim(adjustl(filename)),
     &                NF_NOWRITE,unit)
        
      ENDIF
      IF (ierr.NE.NF_NOERR) THEN
        write(*,*)'Error : cannot open file '//trim(filename)
        write(*,*)'(in phystd/datareadnc.F)'
        write(*,*)'It should be in :',trim(datadir),'/',trim(surfdir)
        write(*,*)'Check that your path to datagcm:',trim(datadir)
        write(*,*)' is correct. You can change it in callphys.def with:'
        write(*,*)' datadir = /absolute/path/to/datagcm'
        write(*,*)'If necessary surface.nc (and other datafiles)'
        write(*,*)' can be obtained online on:'
        write(*,*)' http://www.lmd.jussieu.fr/~lmdz/planets/'//
     &             'LMDZ.GENERIC/datagcm/'
        STOP
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
      pfield(1:iimp1*jjp1)=0
      zdata(1:imd*jmdp1)=0
      zdataS(1:iimp1*jjp1)=0

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

          zdata(:)=1000.*zdata(:)
          longitude(:)=(pi/180.)*longitude(:)
          latitude(:)=(pi/180.)*latitude(:)

          call grid_noro1(360, 180, longitude, latitude, zdata,
     .         iim, jjp1, rlonv, rlatu, zmea,zstd,zsig,zgam,zthe)

          !CALL dump2d(iip1,jjp1,zmea,'zmea')
          !CALL dump2d(iip1,jjp1,zstd,'zstd')
          !CALL dump2d(iip1,jjp1,zsig,'zsig')
          !CALL dump2d(iip1,jjp1,zgam,'zgam')
          !CALL dump2d(iip1,jjp1,zthe,'zthe')

      endif

c-----------------------------------------------------------------------
c   Passage de zdata en grille (imdp1 jmdp1)
c-----------------------------------------------------------------------
      do j=1,jmdp1
          do i=1,imd
              zdataS(i+imdp1*(j-1)) = zdata(i+ngridmixgdat*(j-1))
          enddo
          zdataS(imdp1+imdp1*(j-1)) = zdata(1+ngridmixgdat*(j-1))
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
             phisinit(1:iimp1*jjp1)=0
        else
             phisinit(1:iimp1*jjp1)=pfield(1:iimp1*jjp1)
        endif
      endif

      ENDDO

c-----------------------------------------------------------------------
c    Traitement Phisinit
c-----------------------------------------------------------------------

      phisinit(1:iimp1*jjp1)=1000.*phisinit(1:iimp1*jjp1)
      !CALL dump2d(iimp1,jjp1,phisinit,'Altitude en m')
      phisinit(:)=g*phisinit(:)

c-----------------------------------------------------------------------
c    FIN
c-----------------------------------------------------------------------

      END
