c**********************************************************************

      subroutine column
     $ (co2x,o2x,o3px,h2x,h2ox,h2o2x,tx,nz,iz,zenit,co2colx,o2colx,
     $ o3pcolx,h2o2colx)


c     nov 2002        fgg           first version

c**********************************************************************

      implicit none


c     common variables and constants

      include 'param.h'
      include 'param_v3.h'



c    local parameters and variables



c     input and output variables

      integer    nz
      real       co2colx(nz)                 !column density of CO2 (cm^-2)
      real       o2colx(nz)                  !column density of O2(cm^-2)
      real       o3pcolx(nz)                 !column density of O(3P)(cm^-2)
      real       h2o2colx(nz)                !column density of H2O2(cm^-2)
      real       co2x(nz)                    !density of CO2(cm^-3)
      real       o2x(nz)                     !density of O2(cm^-3)
      real       o3px(nz)                    !density of O(3P)(cm^-3)
      real       h2x(nz)                     !density of H2(cm^-3)
      real       h2ox(nz)                    !density of H2O(cm^-3)
      real       h2o2x(nz)                   !density of H2O2(cm^-3)
      real       jtot(nz)                    !output: heating rate(erg/s)
      real       tx(nz)                      !temperature profile
      real       iz(nz+1)                    !height profile
      real       zenit

c     local variables

      real       xx
      real       grav(nzmax)
      real       ntx(nzmax)
      real       pmmixx(nzmax)
      real       Hmix(nzmax)
      real       nint(nzmax)
      real       co2vmrx(nzmax)
      real       o2vmrx(nzmax)
      real       o3pvmrx(nzmax)
      real       h2vmrx(nzmax)
      real       h2ovmrx(nzmax)
      real       h2o2vmrx(nzmax)
      real       Hco2
      real       Ho3p
      real       Ho2
      real       Hh2o2
      real       chapm(nzmax)



c     local variables

      integer    i,j,k,icol,indexint          !indexes
      real*8     znt

c External functions

       external    ch
       real        ch


c*************************PROGRAM STARTS*******************************

      do i=1,nz
         xx = ( radio + iz(i) ) * 1.e5
         ntx(i)=co2x(i)+o2x(i)+o3px(i)+h2x(i)+h2ox(i)+h2o2x(i)
         co2vmrx(i)=co2x(i)/ntx(i)
         o2vmrx(i)=o2x(i)/ntx(i)
         o3pvmrx(i)=o3px(i)/ntx(i)
         h2vmrx(i)=h2x(i)/ntx(i)
         h2ovmrx(i)=h2ox(i)/ntx(i)
         h2o2vmrx(i)=h2o2x(i)/ntx(i)
         pmmixx(i) = (pmco2*co2x(i)+pmo2*o2x(i) 
     $                + pmo3p*o3px(i))/ntx(i)
         grav(i) = gg * masa /(xx**2)
         Hmix(i) = kboltzman*tx(i)*n_avog / (pmmixx(i)*grav(i))
         nint(i) = -Hmix(i)*ntx(i)*
     $   (exp((iz(i)-iz(i+1))*1e5/Hmix(i))-1)
      end do

      xx = kboltzman * tx(nz) * n_avog / grav(nz)
      Ho3p = xx / pmo3p 
      Hco2 = xx / pmco2
      Ho2  = xx / pmo2
      Hh2o2 = xx / pmh2o2
      o3pcolx(nz) = o3px(nz) * Ho3p
      co2colx(nz) = co2x(nz) * Hco2
      h2o2colx(nz) = h2o2x(nz) * Hh2o2
      o2colx(nz) = o2x(nz) * Ho2
      do i=1,nz-1
         co2colx(i) = co2colx(nz)
         o2colx(i) = o2colx(nz)
         o3pcolx(i) = o3pcolx(nz)
         h2o2colx(i) = h2o2colx(nz)
         do j=i,nz-1
            co2colx(i) = co2colx(i) + nint(j) 
     $                   * ((co2vmrx(j+1) + co2vmrx(j)) / 2.)
            o2colx(i) = o2colx(i) + nint(j)
     $                  * ((o2vmrx(j+1) + o2vmrx(j)) / 2.)
            o3pcolx(i) = o3pcolx(i) + nint(j)
     $                   * ((o3pvmrx(j+1) + o3pvmrx(j)) /2.)
            h2o2colx(i) = h2o2colx(i) +nint(j)
     $                     * ((h2o2vmrx(j+1) + h2o2vmrx(j)) / 2.)
         end do
      end do

      znt=zenit*3.1416/180.
      if(zenit.gt.80.0) then
         do i=1,nz
            chapm(i) =  
     $           ch( (iz(i)+radio)/(Hmix(i)*1.e-5),znt)
         enddo
      end if
      znt=zenit*3.1416/180.
      do i=1,nz
         if(zenit.le.80.0) then
            co2colx(i)=co2colx(i)/cos(znt)
            o2colx(i)=o2colx(i)/cos(znt)
            o3pcolx(i)=o3pcolx(i)/cos(znt)
            h2o2colx(i)=h2o2colx(i)/cos(znt)
         else
            co2colx(i)=co2colx(i)*chapm(i)
            o2colx(i)=o2colx(i)*chapm(i)
            o3pcolx(i)=o3pcolx(i)*chapm(i)
            h2o2colx(i)=h2o2colx(i)*chapm(i)
         end if
      end do


      return


      end





