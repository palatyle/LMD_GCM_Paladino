c**********************************************************************

      subroutine hrtherm(ig,nlayer,
     .      euvmod,rm,nespeuv,tx,iz,zenit,zday,jtot)


c     feb 2002        fgg           first version
c     nov 2002        fgg           second version

c**********************************************************************

      use param_v4_h, only: ninter,nabs,jfotsout,fluxtop,freccen

      implicit none

c     common variables and constants
      include "callkeys.h"


c    local parameters and variables

      real       xabsi(nabs,nlayer) 			!densities
      real       jergs(ninter,nabs,nlayer)
      
      integer    i,j,k,indexint          !indexes
      character  dn


c     input and output variables

      integer    ig  ,euvmod,nlayer 
      integer    nespeuv
      real       rm(nlayer,nespeuv)              !density matrix (cm^-3)
      real       jtot(nlayer)                    !output: heating rate(erg/s)
      real       tx(nlayer)                      !temperature
      real       zenit
      real       iz(nlayer)
      real       zday

      ! tracer indexes for the EUV heating:
!!! ATTENTION. These values have to be identical to those in chemthermos.F90
!!! If the values are changed there, the same has to be done here  !!!
      integer,parameter :: i_co2  =  1
      integer,parameter :: i_co   =  2
      integer,parameter :: i_o    =  3
      integer,parameter :: i_o1d  =  4
      integer,parameter :: i_o2   =  5
      integer,parameter :: i_o3   =  6
      integer,parameter :: i_h    =  7
      integer,parameter :: i_h2   =  8
      integer,parameter :: i_oh   =  9
      integer,parameter :: i_ho2  = 10
      integer,parameter :: i_h2o2 = 11
      integer,parameter :: i_h2o  = 12
      integer,parameter :: i_n    = 13
      integer,parameter :: i_n2d  = 14
      integer,parameter :: i_no   = 15
      integer,parameter :: i_no2  = 16
      integer,parameter :: i_n2   = 17
!      integer,parameter :: i_co2=1
!      integer,parameter :: i_o2=2
!      integer,parameter :: i_o=3
!      integer,parameter :: i_co=4
!      integer,parameter :: i_h=5
!      integer,parameter :: i_h2=8
!      integer,parameter :: i_h2o=9
!      integer,parameter :: i_h2o2=10
!      integer,parameter :: i_o3=12
!      integer,parameter :: i_n2=13
!      integer,parameter :: i_n=14
!      integer,parameter :: i_no=15
!      integer,parameter :: i_no2=17

c*************************PROGRAM STARTS*******************************

      !If nighttime, photoabsorption coefficient set to 0
      if(zenit.gt.140.) then
         dn='n'
         else
         dn='d'
      end if
      if(dn.eq.'n') then
        do i=1,nlayer                                    
	      jtot(i)=0.
        enddo       
        return
      endif 

      !initializations
      jergs(:,:,:)=0.
      xabsi(:,:)=0.
      jtot(:)=0.
      !All number densities to a single array, xabsi(species,layer)
      do i=1,nlayer
         xabsi(1,i)  = rm(i,i_co2)
         xabsi(2,i)  = rm(i,i_o2)
         xabsi(3,i)  = rm(i,i_o)
         xabsi(4,i)  = rm(i,i_h2o)
         xabsi(5,i)  = rm(i,i_h2)
         xabsi(6,i)  = rm(i,i_h2o2)
         !Only if O3, N or ion chemistry requested
         if(euvmod.ge.1) then
            xabsi(7,i)  = rm(i,i_o3)
         endif
         !Only if N or ion chemistry requested
         if(euvmod.ge.2) then
            xabsi(8,i)  = rm(i,i_n2)
            xabsi(9,i)  = rm(i,i_n)
            xabsi(10,i) = rm(i,i_no)
            xabsi(13,i) = rm(i,i_no2)
         endif
         xabsi(11,i) = rm(i,i_co)
         xabsi(12,i) = rm(i,i_h)
      end do

      !Calculation of photoabsortion coefficient
      call jthermcalc_e107(ig,nlayer,euvmod,
     .           rm,nespeuv,tx,iz,zenit,zday)

      !Total photoabsorption coefficient
      do i=1,nlayer
         jtot(i)=0.
        do j=1,nabs
          do indexint=1,ninter
            jergs(indexint,j,i) = jfotsout(indexint,j,i) 
     $              * xabsi (j,i) * fluxtop(indexint)  
     $              / (0.5e9 * freccen(indexint))
            jtot(i)=jtot(i)+jergs(indexint,j,i)
          end do
        end do
      end do

      end

