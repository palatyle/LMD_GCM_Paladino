      SUBROUTINE SUAER
      implicit none
C     
C     Purpose.
C     --------
C     initialize yomaer, the common that contains the
C     radiative characteristics of the aerosols
c     
C     AUTHOR.
C     -------
c     Richard Fournier (1996) Francois Forget (1996)
c     Frederic Hourdin
C     Jean-jacques morcrette *ECMWF*
c     MODIF Francois Forget (2000)
c     MODIF Franck Montmessin (add water ice)
C     ------------------------------------------------------------------
C     
C-----------------------------------------------------------------------
C     
#include "dimensions.h"
#include "dimphys.h"
#include "dimradmars.h"
#include "yomaer.h"
     
c     Aerosol Spectral properties :
#include "aerdust.h"

#ifdef ICE
#include "aerice.h"
#endif
C-----------------------------------------------------------------------
c     
      INTEGER iaer,isun,iir,n

c I/O  of "aerave" (subroutine averaging spectrally sing.scat.parameters) 
      REAL tsun             ! Sun brightness temperature (for SW)
      REAL tsol             ! Surface reference brightness temp (for LW)
      REAL longref          ! reference wavelengths
      REAL longsun(nsun+1)  ! solar band boundaries
      REAL longir(nir+1)    ! IR band boundaries
      REAL epref            ! reference extinction ep at wavelength "longref"
      REAL epav(nir)        ! average ep  (= <Qext>/Qext(longref) if epref=1)
      REAL omegav(nir)      ! Average sing.scat.albedo
      REAL gav(nir)         ! Average assymetry parameter
      REAL QIRTM9sQrefIR    ![Qext averaged over IRTM 9um band]/Qext(longrefir)
                            ! with longrefir defined in dimradmars.h 


C-----------------------------------------------------------------------
c     quelques initialisations a 0
      call zerophys(naerkind*nsun,gvis)
      call zerophys(naerkind*nsun,omegavis)
      call zerophys(naerkind*nsun,QVISsQREF)

      call zerophys(naerkind*nir,gIR)
      call zerophys(naerkind*nir,omegaIR)
      call zerophys(naerkind*nir,QIRsQREF)
c     
C-----------------------------------------------------------------------
C     
C     ----------------------------------------------------------------
C     *       1.    SHORTWAVE COEFFICIENTS
C     ----------------------------------------------------------------
C     
c     Computing average optical properties on both solar bands
c     For pure dust (naerkind=1)
c     
      do iaer=1,naerkind

      tsun=6000.E+0
      longsun(1)=long1vis
      longsun(2)=long2vis
      longsun(3)=long3vis
      longref=longrefvis
      epref=1.E+0         


c     Here, epav is <Qext>/Qext(longrefvis) since epref=1 :

      if (iaer.eq.1) then
      CALL aerave ( ndustvis,
     &     longdustvis,epdustvis,omegdustvis,gdustvis,
     &     longref,epref,tsun
     &     ,nsun,longsun, epav,omegav,gav )
      elseif (iaer.eq.2) then
#ifdef ICE
      CALL aerave ( nicevis,
     &     longicevis,epicevis,omegicevis,gicevis,
     &     longref,epref,tsun
     &     ,nsun,longsun, epav,omegav,gav )
#endif
      endif

       do isun=1,nsun
            QVISsQREF(isun,iaer)=epav(isun)
            gvis(isun,iaer)=gav(isun)
            omegavis(isun,iaer)=omegav(isun)
c TEST
c	    if (iaer.eq.2)  omegavis(isun,iaer)=.86   
c	    if (iaer.eq.2)  gvis(isun,iaer)=-1.  
c	    if (iaer.eq.2) then   !TEST
c              QVISsQREF(isun,iaer)=QVISsQREF(isun,1)
c              gvis(isun,iaer)=gvis(isun,1)
c              omegavis(isun,iaer)=omegavis(isun,1)
c            endif
c END TEST

       enddo

c     

C     ----------------------------------------------------------------
C     *       2.    LONGWAVE COEFFICIENTS
C     ----------------------------------------------------------------

      if (iaer.eq.1) then

c     Computing average optical properties on both solar bands
c     For dust (iaer=1)
c
c     Calcul preliminaire
c     ~~~~~~~~~~~~~~~~~~~
c     Ratio betwen Qext averaged over IRTM 9um band
c     and Qext(longrefir) (longrefir is defined in dimradmars.h)
c     -> useful because the ratio of extinction "solsir"
c      is defined between 0.67um and the IRTM 9um band 
c      (for which it has been estimated)

      tsol=215.D+0
c     IRTM band (T9):
      longir(1)=8.3E-6
      longir(2)=9.7E-6
      longref=longrefir
      epref=1.E+0

c     Here, epav is QavIRTM9/Qext(longrefir) since epref=1 :

       write(*,*) 'Call test 9 micron'
      CALL aerave ( ndustir,
     &     longdustir,epdustir,omegdustir,gdustir,
     &     longref,epref,tsol
     &     ,1,longir,epav,omegav,gav )
       write(*,*) 'OK test 9 micron'
      QIRTM9sQrefIR=epav(1)
      
c     Average scaterring properties of 3 IR bands defined as :
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     [long1ir - long1co2] , [long1co2-long2co2], [long2co2 - long2ir]
c

      tsol=215.D+0
      longir(1)=long1ir
      longir(2)=long1co2
      longir(3)=long2co2
      longir(4)=long2ir
      longref=longrefir
      epref=1.E+0

c     Here, epav is <QIR>/Qext(longrefir) since epref=1
      CALL aerave ( ndustir,
     &     longdustir,epdustir,omegdustir,gdustir,
     &     longref,epref,tsol
     &     ,nir-1,longir,epav,omegav,gav )

c     Computing  <QIR>/Qext(longrefvis)
      DO iir=1,4
         epav(iir)=  epav(iir) / (QIRTM9sQrefIR * solsir)
      ENDDO

      elseif (iaer.eq.2) then

c     Average scaterring properties of 3 IR bands defined as :
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     [long1ir - long1co2] , [long1co2-long2co2], [long2co2 - long2ir]
c

      tsol=215.D+0
      longir(1)=long1ir
      longir(2)=long1co2
      longir(3)=long2co2
      longir(4)=long2ir
      longref=longrefir
      epref=1.E+0

c     Here, epav is <QIR>/Qext(longrefir) since epref=1
#ifdef ICE
      CALL aerave ( niceir,
     &     longiceir,epiceir,omegiceir,giceir,
     &     longref,epref,tsol
     &     ,nir-1,longir,epav,omegav,gav )

c     Computing  <QIR>/Qext(longrefvis)
      DO iir=1,4
         epav(iir)=  epav(iir) / solsirice
      ENDDO
#endif

      endif
c  Single scattering properties in each of the "nir" bands (cf. dimramars.h)
c  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c  iir=1 : central 15um CO2 bands   
      QIRsQREF(1,iaer)=epav(2)
      omegaIR(1,iaer)=omegav(2)
      gIR(1,iaer)=gav(2)

c  iir=2 : CO2 band wings (same properties than for central part)
      QIRsQREF(2,iaer)=epav(2)
      omegaIR(2,iaer)=omegav(2)
      gIR(2,iaer)=gav(2)

c  iir=3 : 9 um band [long1ir - long1co2]
      QIRsQREF(3,iaer)=epav(1)
      omegaIR(3,iaer)=omegav(1)
      gIR(3,iaer)=gav(1)

c  iir=4 : Far IR    [long2co2 - long2ir]
      QIRsQREF(4,iaer)=epav(3)
      omegaIR(4,iaer)=omegav(3)
      gIR(4,iaer)=gav(3)

c      if (iaer.eq.2) then  !TEST
c	do iir=1,4
c	  QIRsQREF(iir,iaer)=QIRsQREF(iir,1)
c          omegaIR(iir,iaer)= omegaIR(iir,1)
c          gIR(iir,iaer)=gIR(iir,1)
c        enddo
c      endif

C     ----------------------------------------------------------------
C     Output on screen 
C     ----------------------------------------------------------------
      if (iaer.eq.1) then
        PRINT*,'PURE DUST PROPERTIES :'
        PRINT*
        PRINT*,'Rapport Solaire/IR :',solsir
        PRINT*
      elseif (iaer.eq.2) then
#ifdef ICE
        PRINT*,'ICE PROPERTIES :'
        PRINT*
        PRINT*,'Rapport Solaire/IR :',solsirice
        PRINT*
#endif
      endif

      PRINT *,'Les donnees spectrales :'
      PRINT *,'Solaire (SW) ---->'
      PRINT *,'<Qext>/Qext(0.67um) ; omega ; g'
      DO isun=1,nsun
         PRINT *,QVISsQREF(isun,iaer),omegavis(isun,iaer)
     &    ,gvis(isun,iaer)
      ENDDO
      PRINT *,'Thermal IR (LW) ---->'
      PRINT *,'<Qext>/Qext(0.67um) ; omega ; g'
      DO iir=1,nir
         PRINT *,QIRsQREF(iir,iaer),omegaIR(iir,iaer),gIR(iir,iaer)
      ENDDO
c     

      print *,'Dans le co2 on prend <Qabs>/Qext(0.67um) =',
     &  QIRsQREF(1,iaer)*(1-omegaIR(1,iaer))  
      write(*,*)

      enddo  ! Loop on iaer
      


      RETURN
      END
