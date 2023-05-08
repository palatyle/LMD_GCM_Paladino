c=======================================================================
      SUBROUTINE datareadnc(relief,phisinit,alb,ith,z0,
     &                    zmea,zstd,zsig,zgam,zthe,
     &                    hmons,summit,base,zavg)
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

! to use  'getin'
      use ioipsl_getincom, only: getin 
      use comconst_mod, only: g,pi
      use datafile_mod, only: datadir
      use avg_horiz_mod, only: avg_horiz
      use mvc_horiz_mod, only: mvc_horiz

      implicit none

      include "dimensions.h"
      include "paramet.h"
      include "comgeom.h"
      include "netcdf.inc"

c=======================================================================
c   Declarations:
C=======================================================================

      INTEGER    imd,jmd,imdp1,jmdp1
      parameter    (imd=360,jmd=179,imdp1=361,jmdp1=180)

      INTEGER    iimp1
      parameter    (iimp1=iim+1-1/iim)

! Arguments:
      CHARACTER(len=3),intent(inout) :: relief
      REAL,intent(out) :: phisinit(iimp1*jjp1)
      REAL,intent(out) :: alb(iimp1*jjp1)
      REAL,intent(out) :: ith(iimp1*jjp1)
      REAL,intent(out) :: z0(iimp1*jjp1)
      REAL,intent(out) :: zmea(imdp1*jmdp1)
      REAL,intent(out) :: zstd(imdp1*jmdp1)
      REAL,intent(out) :: zsig(imdp1*jmdp1)
      REAL,intent(out) :: zgam(imdp1*jmdp1)
      REAL,intent(out) :: zthe(imdp1*jmdp1)
      REAL,intent(out) :: hmons(imdp1*jmdp1) !CW17,361*180 hmons
      REAL,intent(out) :: summit(imdp1*jmdp1)
      REAL,intent(out) :: base(imdp1*jmdp1)
      REAL,intent(out) :: zavg(imdp1*jmdp1)
      
! Local variables:
      REAL        zdata(imd*jmdp1)
      REAL        zdataS(imdp1*jmdp1)
      REAL        pfield(iimp1*jjp1)

      INTEGER     ierr

      INTEGER   unit,nvarid

      INTEGER    i,j,k

      INTEGER klatdat,ngridmxgdat
      PARAMETER (klatdat=180,ngridmxgdat=360)

c    on passe une grille en rlonu rlatv et im+1 jm a interp_horiz)

      REAL longitude(imd),latitude(jmdp1) ! Pour lecture des donnees
      REAL rlonud(imdp1),rlatvd(jmd)

      CHARACTER*20 string
      DIMENSION string(0:7)


!#include "lmdstd.h"
!#include "fxyprim.h"

      pi=2.*ASIN(1.)

c=======================================================================
c    rlonud, rlatvd
c=======================================================================

c-----------------------------------------------------------------------
c    Lecture NetCDF des donnees latitude et longitude
c-----------------------------------------------------------------------
      write(*,*) 'datareadnc: opening file surface.nc'

      datadir="/u/lmdz/WWW/planets/mars/datadir" ! default path to surface.nc
      call getin("datadir",datadir) ! but users may specify another path
      
      ierr = NF_OPEN (trim(datadir)//'/surface.nc',
     &  NF_NOWRITE,unit)
      IF (ierr.NE.NF_NOERR) THEN
        write(*,*)'Error : cannot open file surface.nc '
        write(*,*)'(in phymars/datareadnc.F)'
        write(*,*)'It should be in :',trim(datadir),'/'
        write(*,*)'1) You can set this path in the 
     & callphys.def file:'
        write(*,*)'   datadir=/path/to/the/datafiles'
        write(*,*)'2) If necessary, surface.nc (and other datafiles)'
        write(*,*)'   can be obtained online on:'
        write(*,*)'http://www.lmd.jussieu.fr/~lmdz/planets/mars/datadir'
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

      string(0) = 'z0'
      string(1) = 'albedo'
      string(2) = 'thermal'
      if (relief.ne.'pla') then
        write(*,*) ' MOLA topography'
        relief = 'MOL'
          string(3) = 'z'//relief
      else
          string(3) = 'zMOL'  ! pour qu''il lise qqchose sur le fichier
                            ! remise a 0 derriere
      endif
      string(4) = 'zMOL'    ! lecture pour calcul topog. sous-maille
 

      DO k=0,4
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
      if (ierr.ne.nf_noerr) then
        write(*,*) 'datareadnc error, cannot find ',trim(string(k))
        write(*,*) ' in file ',trim(datadir),'/surface.nc'
        stop
      endif
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(unit, nvarid, zdata)
#else
      ierr = NF_GET_VAR_REAL(unit, nvarid, zdata)
#endif
      if (ierr.ne.nf_noerr) then
        write(*,*) 'datareadnc: error failed loading ',trim(string(k))
        stop
      endif

c-----------------------------------------------------------------------
c        Cas particulier "Francois Lott" ( k=4 ) (relief sous-maille)
c-----------------------------------------------------------------------
      if (k.eq.4) then

          zdata(:)=1000.*zdata(:)
          longitude(:)=(pi/180.)*longitude(:)
          latitude(:)=(pi/180.)*latitude(:)

          call grid_noro1(360, 180, longitude, latitude, zdata,
     .         iim, jjp1, rlonv, rlatu, zmea,zstd,zsig,zgam,zthe)
     
          !CW17
          call avg_horiz(imd,jmdp1,iim,jjm,longitude,
     .                latitude,rlonv,rlatu,zdata,pfield)

          do j=1,jjp1 !CW17 49, iimp1=65, the last column = first column
             pfield(iimp1*j) =  pfield(1+iimp1*(j-1))
          enddo 
          do i=1,iimp1*jjp1
c             if (pfield(i) .ne. -999999.) then
                zavg(i) = pfield(i)
c             else
c               zavg(i)=-999999.
c             endif
          enddo

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

      if (k.eq.0) then                    ! z0
         z0(1:iimp1*jjp1)=pfield(1:iimp1*jjp1)*.01
         ! multiplied by 0.01 to have z0 in m
      elseif (k.eq.1) then                    ! albedo
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

      phisinit(1:iimp1*jjp1)=1000.*phisinit(1:iimp1*jjp1)
      phisinit(:)=g*phisinit(:)

c-----------------------------------------------------------------------
c    FIN
c-----------------------------------------------------------------------
 
c=======================================================================
c<<<<<<<<<<<<<<<<<<<<<<<add new dataset>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! name of dataset
      string(5) = 'hmons' !subgrid hmons
c     the following data could be useful in future, but currently,
c     we don't need to put them into restartfi.nc
      string(6) = 'summit' !subgrid summit
      string(7) = 'base'   !base of summit (not subgrid)
c     string(8) = 'bottom' !subgrid base
      do k=5,7
          write(*,*) 'string',k,string(k)  
c-----------------------------------------------------------------------
c    Lecture NetCDF  
c-----------------------------------------------------------------------

          ierr = NF_INQ_VARID (unit, string(k), nvarid)
          if (ierr.ne.nf_noerr) then
            write(*,*) 'datareadnc error, cannot find ',trim(string(k))
            write(*,*) ' in file ',trim(datadir),'/surface.nc'
            stop
          endif

c-----------------------------------------------------------------------
c    initialisation
c-----------------------------------------------------------------------
          call initial0(iimp1*jjp1,pfield)
          call initial0(imd*jmdp1,zdata)
          call initial0(imdp1*jmdp1,zdataS)

#ifdef NC_DOUBLE
          ierr = NF_GET_VAR_DOUBLE(unit, nvarid, zdata)
#else
          ierr = NF_GET_VAR_REAL(unit, nvarid, zdata)
#endif
          if (ierr.ne.nf_noerr) then
            write(*,*) 'datareadnc: error failed loading ',
     .                  trim(string(k))
            stop
          endif

c-----------------------------------------------------------------------
c   Passage de zdata en grille (imdp1 jmdp1)
c-----------------------------------------------------------------------
          do j=1,jmdp1     !  180
              do i=1,imd   !  360
                  !copy zdata to zdataS, line by line
                  !      i+ 361 *(j-1)          i+  360*(j-1)
                  zdataS(i+imdp1*(j-1)) = zdata(i+ngridmxgdat*(j-1))
              enddo
              ! the last column = the first column of zdata
              zdataS(imdp1+imdp1*(j-1)) = zdata(1+ngridmxgdat*(j-1))
          enddo

          if (k .eq. 5 .or. k .eq. 6 .or. k .eq. 7) then 
              ! hmons, summit, keep the maximum of subgrids 
              call mvc_horiz(imd,jmdp1,iim,jjm,longitude,
     .                latitude,rlonv,rlatu,zdata,pfield)
          endif


c-----------------------------------------------------------------------
c    Periodicite    
c-----------------------------------------------------------------------

          do j=1,jjp1 ! 49, iimp1=65, the last column = first column
             pfield(iimp1*j) =  pfield(1+iimp1*(j-1))
          enddo 
 
c-----------------------------------------------------------------------
c    Sauvegarde des champs    
c-----------------------------------------------------------------------

          if (k.eq.5) then ! hmons
             do i=1,iimp1*jjp1
                  if (pfield(i) .ne. -999999.) then
                    hmons(i) = pfield(i)
                  else
                    hmons(i)=0.
                  endif
             enddo
          endif
          
          if (k.eq.6) then  ! summit
             do i=1,iimp1*jjp1
                  if (pfield(i) .ne. -999999.) then             
                    summit(i) = pfield(i)
                  else
                    summit(i)=0.
                  endif
             enddo
          endif

          if (k.eq.7) then  ! base
             do i=1,iimp1*jjp1
                  if (pfield(i) .ne. -999999.) then             
                    base(i) = pfield(i)
                  else
                    base(i)=0.
                  endif
             enddo
          endif

      enddo
c<<<<<<<<<<<<<<<<<<<<<<<<<done add new dataset>>>>>>>>>>>>>>>>>>>>>>>>>>
c=======================================================================
      END
