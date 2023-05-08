c**********************************************************************

      subroutine hrtherm
     $ (co2x,o2x,o3px,h2x,h2ox,h2o2x,aux1,aux2,tx,nz,iz,date,zenit,jtot)


c     feb 2002        fgg           first version
c     nov 2002        fgg           second version

c**********************************************************************

      implicit none

c     common variables and constants

      include 'param.h'
      include 'param_v3.h'
      include 'callkeys.h'


c    local parameters and variables

      real       xabsi(nabs,nzmax) 			!densities
      real       nada
      real       jergs(ninter,nabs,nzmax)
      
      integer    i,j,k,indexint          !indexes
      character  dn


c     input and output variables

      integer    nz                          !number of layers
      real       co2x(nz)                    !density of CO2(cm^-3)
      real       o2x(nz)                     !density of O2(cm^-3)
      real       o3px(nz)                    !density of O(3P)(cm^-3)
      real       h2x(nz)                     !density of H2(cm^-3)
      real       h2ox(nz)                    !density of H2O(cm^-3)
      real       h2o2x(nz)                   !density of H2O2(cm^-3)
      real       aux1(nz)	  	             !auxiliar variable
      real       aux2(nz)                    !auxiliar variable
      real       jtot(nz)                    !output: heating rate(erg/s)
      real       tx(nz)                      !temperature
      real       date
      real       zenit
      real       iz(nz)

      logical firstcall
      save firstcall
      data firstcall /.true./

c*************************PROGRAM STARTS*******************************

c      if (firstcall) then
c        if(.not. thermochem) call param_read
c        firstcall= .false.
c      endif

      if(zenit.gt.100.) then
         dn='n'
         else
         dn='d'
      end if
      if(dn.eq.'n') then
        do i=1,nz                                    
	      jtot(i)=0.
        enddo       
        return
      endif 
     
      do i=1,nz
         xabsi(1,i) = co2x(i)
         xabsi(2,i) = o2x(i)
         xabsi(3,i) = o3px(i)
         xabsi(4,i) = h2ox(i)
         xabsi(5,i) = h2x(i)
         xabsi(6,i) = h2o2x(i)
         jtot(i) = 0.
      end do

      if(.not. thermochem) then
        call jthermcalc
     $ (co2x,o2x,o3px,h2x,h2ox,h2o2x,aux1,aux2,tx,nz,iz,date,zenit)
      endif

      do i=1,nz
        do j=1,nabs
          do indexint=1,33
            jergs(indexint,j,i) = jfotsout(indexint,j,i) 
     $              * xabsi (j,i) * fluxtop(indexint)  
     $              / (0.5e9 * freccen(indexint))
            jtot(i)=jtot(i)+jergs(indexint,j,i)
          end do
        end do
      end do

      return

      end

