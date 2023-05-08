      SUBROUTINE load_ksi(ksive)
      
      use dimphy
      IMPLICIT none

#include "YOMCST.h"
#include "comcstVE.h"
C
C     ------------------------------------------------------------------
C
C     PURPOSE.
C     --------
C
c     This routine loads the longwave matrix of factors Ksi
c     
c     The Ksi matrixes have been computed by Vincent Eymet
C
C     AUTHOR.
C     -------
C        Sebastien Lebonnois
C
C     MODIFICATIONS.
C     --------------
C        
c   New ksi matrix: possibility of different cloud model fct of lat   05/2014
C     ------------------------------------------------------------------
C
C* ARGUMENTS:
C
c inputs
      real   psurf(klon)           ! Surface pressure
c outputs
      real   ksive(0:klev+1,0:klev+1,nnuve,nbmat)  ! ksi matrixes in Vincent's file

c local variables
      integer i,j,isza,ips,band,pve,sve,nlve
      integer lat,Nb,m,mat
      character*9 tmp1
      character*100 file
      CHARACTER*2 str2
      CHARACTER*3 str3
      CHARACTER*10 format_lect 
      real   lambda(nnuve)            ! wavelenght in table (mu->m, middle of interval)
      real   lambdamin(nnuve),lambdamax(nnuve) ! in microns
      real   dlambda                  ! in microns 

      nlve = klev

cc      GG modif below 
c----------------------------------
c   Initialisation of values to 0 
c     (for all vertical levels)
c----------------------------------

      ksive(0:klev+1,0:klev+1,1:nnuve,1:nbmat) = 0.0

c ------------------------
c Loading the ksi file
c ------------------------

      file = "ksi_global.txt"
      open(10,file=file)
      
      read(10,*)
      read(10,*) nlatve
      read(10,*)

      write(*,*) 'This is subroutine load_ksi'
      write(*,*) 'Nb of lat bands:',nlatve
      
      do lat=1,nlatve
        read(10,*) !line for lat range
        read(10,*) indexve(lat)
        read(10,*) nbpsve(lat)
        read(10,*)
        read(10,*) nbszave(lat)
        read(10,*)
	
	do isza=1,nbszave(lat)
	  do ips=1,nbpsve(lat)
	  
        read(10,*) (tmp1,j=1,3),mat    !line for matrix number
        read(10,*) (tmp1,j=1,2),pve
        psurfve(ips,lat) = pve*1.e5    ! pve in bar, psurfve in Pa
        read(10,*) (tmp1,j=1,3),sve 
	szave(isza,lat) = cos(sve*RPI/180.) ! conversion in mu0
        read(10,*)
        read(10,*) m,Nb
cc      GG modif below  
        if (nlve.le.78.and.m.ne.nlve) then
         write(*,*) 'This is subroutine load_ksi'
         print*,'Dimension problem between ksi.txt and nlve'
         print*,'N levels = ',m,nlve
         stop
        endif
        if (Nb.ne.nnuve) then
         write(*,*) 'This is subroutine load_ksi'
         print*,'Dimension problem between ksi.txt and nnuve'
         print*,'N freq = ',Nb,nnuve
         stop
        endif
c     Now reading ksi matrix index "mat"
        if ((m+2).ge.100) then
          write(str3,'(i3.3)') m+2
          format_lect='('//str3//'e17.9)'
        else
          write(str2,'(i2.2)') m+2
          format_lect='('//str2//'e17.9)'
        endif
        do band=1,Nb
         read(10,*) lambdamin(band),lambdamax(band)
         do i=0,m+1
            read(10,format_lect) (ksive(i,j,band,mat),j=0,m+1) ! no unit
         enddo                  ! i
        enddo                     ! band
c       print*,"Matrice ",mat," lue"
c       print*,"   psurf=",psurfve(ips,lat)," bars"
        if (mat+1.gt.nbmat) then
         write(*,*) 'This is subroutine load_ksi'
         print*,'Dimension problem between ksi.txt and nbmat'
         print*,'(max number of matrixes)'
         print*,'nbmat (in comcstVE.h) should be raised'
         stop
        endif

	  enddo    ! ips
	enddo      ! isza
      enddo        ! lat
      
      write(*,*) 'Total nb of matrixes:',mat
      
      close(10)

c central wavelength and wavelength width
      do band=1,nnuve
         lambda(band)=(lambdamin(band)+lambdamax(band))/2.*1.e-6   ! in m
         dlambda     =(lambdamax(band)-lambdamin(band))            ! in microns
c        print*,band,lambdamin(band),dlambda,lambdamax(band)

c sign convention for ksi,
c and taking into account the wavelength width (in microns):
         do mat=1,nbmat
         do i=0,nlve+1
           do j=0,nlve+1
              ksive(i,j,band,mat) = +ksive(i,j,band,mat)*dlambda    ! in µm
           enddo
         enddo
         enddo
c computing coeff al and bl for Planck luminance 
         al(band) = 2.*RHPLA*RCLUM*RCLUM/(lambda(band))**5.
c in W/m²/m
c We need W/m²/µm :
     .                * 1.e-6
         bl(band) = RHPLA*RCLUM/(RKBOL*lambda(band))
      enddo
      
      print*,"LOAD_KSI OK"

      return
      end

