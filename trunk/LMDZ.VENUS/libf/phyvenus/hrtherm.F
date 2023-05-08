c**********************************************************************

      subroutine hrtherm(ig,euvmod,rm,nespeuv,tx,iz,zenit,zday,jtot)


c     feb 2002        fgg           first version
c     nov 2002        fgg           second version

c**********************************************************************
      use dimphy
      use conc
      implicit none

c     common variables and constants


#include "param.h"
#include "param_v4.h"
#include "clesphys.h"


c    local parameters and variables

      real       xabsi(nabs,klev) 			!densities (cm^-3)
      real       jergs(ninter,nabs,klev)
      
      integer    i,j,k,indexint          !indexes
      character  dn


c     input and output variables

      integer    ig  ,euvmod 
      integer    nespeuv 
      real       rm(klev,nespeuv)              !density matrix (cm^-3)
      real       jtot(klev)                    !output: heating rate(erg/s cm3)
      real       tx(klev)                      !temperature
      real       zenit
      real       iz(klev)
      real       zday

      ! tracer indexes for the EUV heating:
!!! ATTENTION. These values have to be identical to those in euvheat.F90
!!! If the values are changed there, the same has to be done here  !!!

      integer,parameter :: i_co2=1
      integer,parameter :: i_n2=13
      integer,parameter :: i_n=14
      integer,parameter :: i_o=3
      integer,parameter :: i_co=4


c*************************PROGRAM STARTS*******************************

      !If nighttime, photoabsorption coefficient set to 0
      if(zenit.gt.90.) then  !140 in the martian routine
         dn='n'
         else
         dn='d'
      end if
      if(dn.eq.'n') then
        do i=1,klev                                    
	      jtot(i)=0.
        enddo       
        return
      endif 

      !initializations
      jergs(:,:,:)=0.
      xabsi(:,:)=0.
      jtot(:)=0.
      !All number densities to a single array, xabsi(species,layer)
      ! WARNING xabs(nabs,nlev), j=1,nabs --> the values of j should 
      !         be the same for xabs than for jfotsout(indexint,j,i) 
      !
      do i=1,klev
         xabsi(1,i)  = rm(i,i_co2)
         xabsi(3,i)  = rm(i,i_o)
         xabsi(8,i)  = rm(i,i_n2)
         xabsi(11,i)  = rm(i,i_co)

c         xabsi(6,i)  = rm(i,i_h2o2)
         !Only if O3, N or ion chemistry requested
c         if(euvmod.ge.1) then
c            xabsi(7,i)  = rm(i,i_o)
c         endif
         !Only if N or ion chemistry requested
c         if(euvmod.ge.2) then
c            xabsi(8,i)  = rm(i,i_n2)
c            xabsi(9,i)  = rm(i,i_n)
c            xabsi(10,i) = rm(i,i_no)
c            xabsi(13,i) = rm(i,i_no2)
c         endif
      end do

      !Calculation of photoabsortion coefficient
      if(solvarmod.eq.0) then
         call jthermcalc(ig,euvmod,rm,nespeuv,tx,iz,zenit)
      else if (solvarmod.eq.1) then
         call jthermcalc_e107(ig,euvmod,rm,nespeuv,tx,iz,zenit,zday)
         do indexint=1,ninter
            fluxtop(indexint)=1.
         enddo
      endif

      !Total photoabsorption coefficient    !  erg/(s*cm3) 
      do i=1,klev
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

      return

      end

