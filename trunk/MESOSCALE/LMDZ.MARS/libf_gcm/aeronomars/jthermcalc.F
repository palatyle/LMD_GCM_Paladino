c**********************************************************************

      subroutine jthermcalc
     $ (co2x,o2x,o3px,h2x,h2ox,h2o2x,aux1,aux2,tx,nz,iz,date,zenit)


c     feb 2002        fgg           first version
c     nov 2002        fgg           second version
c
c modified from paramhr.F
c MAC July 2003
c**********************************************************************

      implicit none

c     common variables and constants

      include 'param.h'
      include 'param_v3.h'

c    local parameters and variables

      real       xabsi(nabs,nzmax) !densities
      real       nada
      real       co2colx2(nzmax,ninter)
      real       o2colx2(nzmax,ninter)
      real       o3pcolx2(nzmax,ninter)
      real       co2colx(nzmax)              !column density of CO2 (cm^-2)
      real       o2colx(nzmax)               !column density of O2(cm^-2)
      real       o3pcolx(nzmax)              !column density of O(3P)(cm^-2)
      real       h2o2colx(nzmax)             !column density of H2O2(cm^-2)
      real       t2(nzmax)
      real       coltemp(nzmax)
      real       dt(nzmax)
      real       sigma(ninter,nzmax)
      real       alfa(ninter,nzmax)
      real       xx
      real       grav(nzmax)
      
      integer    i,j,k,icol,indexint          !indexes
      character  dn

c     input and output variables

      integer    nz                          !number of layers
      real       co2x(nz)                    !density of CO2(cm^-3)
      real       o2x(nz)                     !density of O2(cm^-3)
      real       o3px(nz)                    !density of O(3P)(cm^-3)
      real       h2x(nz)                     !density of H2(cm^-3)
      real       h2ox(nz)                    !density of H2O(cm^-3)
      real       h2o2x(nz)                   !density of H2O2(cm^-3)
      real       aux1(nz)	             !auxiliar variable
      real       aux2(nz)                    !auxiliar variable
      real       tx(nz)                      !temperature
      real       date
      real       zenit
      real       iz(nz)


c     variables used in interpolation

      real aux3(nz2), aux4(nz2)
                                              !auxiliar variables
      real       limdown                      !limits for interpolation
      real       limup                        !        ""

c*************************PROGRAM STARTS*******************************

      if(zenit.gt.100.) then
         dn='n'
         else
         dn='d'
      end if
      if(dn.eq.'n') then
        return
      endif 
     
      do i=1,nz
         xabsi(1,i) = co2x(i)
         xabsi(2,i) = o2x(i)
         xabsi(3,i) = o3px(i)
         xabsi(4,i) = h2ox(i)
         xabsi(5,i) = h2x(i)
         xabsi(6,i) = h2o2x(i)
      end do

      do i=1,nz
         t2(i)=tx(i)
         if(t0(i).lt.195.0) t0(i)=195.0
         if(t0(i).gt.295.0) t0(i)=295.0
         if(t2(i).lt.195.0) t2(i)=195.0
         if(t2(i).gt.295.0) t2(i)=295.0
      end do


c    COLUMN CALCULATION 

      call column(co2x,o2x,o3px,h2x,h2ox,h2o2x,tx,nz,iz,zenit,co2colx
     $,o2colx,o3pcolx,h2o2colx)

      coltemp(nz)=co2colx(nz)*abs(t2(nz)-t0(nz))
      dt(nz)=coltemp(nz)/co2colx(nz)
      do i=nz-1,1,-1
        coltemp(i)=coltemp(i+1)+ 
     $         ( co2x(i) + co2x(i+1) ) * 0.5 
     $         * 1e5 * (iz(i+1)-iz(i)) * abs(t2(i)-t0(i))
        dt(i)=coltemp(i)/co2colx(i)
      end do


c     Calculo la seccion eficaz de CO2 a la temperatura t0 a cada altura
c CALCULATION OF CO2 CROSS-SECTION AT TEMPERATURE T0 AT EACH HEIGHT

      do i=1,nz
         do indexint=24,32
           sigma(indexint,i)=co2crsc195(indexint-23)*
     $         (1+((co2crsc295(indexint-23)/co2crsc195(indexint-23))-1.)
     $         *(t0(i)-195.)/100.)
           if(t0(i).ne.295.) then
             alfa(indexint,i)=((co2crsc295(indexint-23)
     $                        /sigma(indexint,i))-1.)/(295.-t0(i))
           else 
             alfa(indexint,i)=0.
           end if
         end do
      end do


c     Comienza la interpolacion
c INTERPOLATION BEGINS
c************************************************
c     o2,0.1,5.0
c************************************************
      
      indexint=1
      limdown=1e-20
      limup=1e26
       do i=1,nz
         aux2(nz-i+1) = o2colx(i) + o3pcolx(i)
       end do
       do i=1,nz2
         aux3(i) = jabsifotsint(indexint,2,nz2-i+1)
         aux4(i) = c23(nz2-i+1)
       enddo
       call interpfast 
     $        ( aux1, aux2, nz,  aux3, aux4, nz2 , limdown , limup)
       do i=1,nz
         jfotsout(indexint,2,i) = aux1(nz-i+1)
       enddo

c************************************************
c     o3p,0.1,5.0
c************************************************

      indexint=1
      limdown=1e-20
      limup=1e26
         do i=1,nz
            aux2(nz-i+1) = o2colx(i) + o3pcolx(i)
         end do
         do i=1,nz2
            aux3(i) = jabsifotsint(indexint,3,nz2-i+1)
            aux4(i) = c23(nz2-i+1)
         enddo
         call interpfast 
     $        ( aux1, aux2, nz,  aux3, aux4, nz2 , limdown , limup)
         do i=1,nz
            jfotsout(indexint,3,i) = aux1(nz-i+1)
         enddo


c************************************************
c     h2,0.1,5.0
c************************************************

      indexint=1
      limdown=1e-20
      limup=1e26
         do i=1,nz
            aux2(nz-i+1) = o2colx(i) + o3pcolx(i)
         end do
         do i=1,nz2
            aux3(i) = jabsifotsint(indexint,5,nz2-i+1)
            aux4(i) = c23(nz2-i+1)
         enddo
         call interpfast 
     $        ( aux1, aux2, nz,  aux3, aux4, nz2 , limdown , limup )
         do i=1,nz
            jfotsout(indexint,5,i) = aux1(nz-i+1)
         enddo


c************************************************
c     interpola para el o2 entre 5 y 90.8nm
c************************************************

      limdown=1e-20
      limup=1e12
      do indexint=2,16
         co2colx2(nz,indexint)= co2colx(nz) 
     $              * crscabsi2(1,indexint)
         o2colx2(nz,indexint)= o2colx(nz) 
     $              * crscabsi2(2,indexint)
         o3pcolx2(nz,indexint)= o3pcolx(nz) 
     $              * crscabsi2(3,indexint)
         do i=nz-1,1,-1
            co2colx2(i,indexint) = co2colx2(i+1,indexint) + 
     $         ( xabsi(1,i) + xabsi(1,i+1) ) * 0.5 
     $           * 1e5 * (iz(i+1)-iz(i)) * crscabsi2(1,indexint)
            o2colx2(i,indexint) = o2colx2(i+1,indexint) + 
     $         ( xabsi(2,i) + xabsi(2,i+1) ) * 0.5 
     $         * 1e5 * (iz(i+1)-iz(i)) * crscabsi2(2,indexint)
            o3pcolx2(i,indexint) = o3pcolx2(i+1,indexint) + 
     $         ( xabsi(3,i) + xabsi(3,i+1) ) * 0.5 
     $         * 1e5 * (iz(i+1)-iz(i)) * crscabsi2(3,indexint)
         end do
         do i=1,nz
            aux2(nz-i+1) = co2colx2(i,indexint)+
     $      o2colx2(i,indexint)+
     $      o3pcolx2(i,indexint)
         end do
         do i=1,nz2
            aux3(i) = jabsifotsint(indexint,2,nz2-i+1)
            aux4(i) = c123(indexint,nz2-i+1)
         enddo
          call interpfast 
     $        ( aux1, aux2, nz,  aux3, aux4, nz2 , limdown , limup )
          do i=1,nz
            jfotsout(indexint,2,i) = aux1(nz-i+1) 
          enddo
      end do
      

c************************************************
c     idem para o3p entre 5 y 90.8nm
c************************************************

      do indexint=2,16
         do i=1,nz
            aux2(nz-i+1) = co2colx2(i,indexint) +
     $      o2colx2(i,indexint) +
     $      o3pcolx2(i,indexint)
         end do
         do i=1,nz2
            aux3(i) = jabsifotsint(indexint,3,nz2-i+1)
            aux4(i) = c123(indexint,nz2-i+1)
         enddo
         call interpfast 
     $        ( aux1, aux2, nz,  aux3, aux4, nz2 , limdown , limup )
         do i=1,nz
            jfotsout(indexint,3,i) = aux1(nz-i+1)
         enddo
      end do


c************************************************
c     idem para co2 entre 5 y 90.8nm
c************************************************

      do indexint=2,16
         do i=1,nz
            aux2(nz-i+1) = co2colx2(i,indexint) +
     $      o2colx2(i,indexint) +
     $      o3pcolx2(i,indexint)
         end do
         do i=1,nz2
            aux3(i) = jabsifotsint(indexint,1,nz2-i+1)
            aux4(i) = c123(indexint,nz2-i+1)
         enddo
         call interpfast 
     $        ( aux1, aux2, nz,  aux3, aux4, nz2 , limdown , limup )
         do i=1,nz
            jfotsout(indexint,1,i) = aux1(nz-i+1)
         enddo
      end do


c************************************************
c     idem para h2 entre 5 y 80.6nm
c************************************************

      do indexint=2,15
         do i=1,nz
            aux2(nz-i+1) = co2colx2(i,indexint) +
     $      o2colx2(i,indexint) +
     $      o3pcolx2(i,indexint)
         end do
         do i=1,nabs
         end do
         do i=1,nz2
            aux3(i) = jabsifotsint(indexint,5,nz2-i+1)
            aux4(i) = c123(indexint,nz2-i+1)
         enddo
         call interpfast 
     $        ( aux1, aux2, nz,  aux3, aux4, nz2 , limdown , limup )
         do i=1,nz
            jfotsout(indexint,5,i) = aux1(nz-i+1)
         enddo
      end do


c************************************************
c      co2,90.9,119.5
c************************************************

      limdown=1e-20
      limup=1e26
      do indexint=17,24
         do i=1,nz
            aux2(nz-i+1) = co2colx(i) + o2colx(i)
         end do
         do i=1,nz2
            aux3(i) = jabsifotsint(indexint,1,nz2-i+1)
            aux4(i) = c12(nz2-i+1)
         enddo
         call interpfast 
     $        ( aux1, aux2, nz,  aux3, aux4, nz2 , limdown , limup )
         do i=1,nz
c           print*,indexint,i,aux2(i),aux1(i),aux3(i),aux4(i)
           jfotsout(indexint,1,i) = aux1(nz-i+1)
           if(indexint.eq.24) then
           if(sigma(indexint,i)*alfa(indexint,i)*coltemp(i).lt.60.) then
             jfotsout(indexint,1,i) = aux1(nz-i+1)
     $            *exp(-sigma(indexint,i)*alfa(indexint,i)*coltemp(i))*
     $           (1+alfa(indexint,i)*(t2(i)-t0(i)))
           else 
             jfotsout(indexint,1,i)=0.
           end if
           end if
c           print*,indexint,i,(jfotsout(indexint,j,i),j=1,nabs)
         enddo
      end do


c************************************************
c     o2,90.9,119.5
c************************************************
      
      do indexint=17,24
         do i=1,nz
            aux2(nz-i+1) = o2colx(i) + co2colx(i)
         end do
         do i=1,nz2
            aux3(i) = jabsifotsint(indexint,2,nz2-i+1)
            aux4(i) = c12(nz2-i+1)
         enddo
         call interpfast 
     $        ( aux1, aux2, nz,  aux3, aux4, nz2 , limdown , limup )
         do i=1,nz
            jfotsout(indexint,2,i) = aux1(nz-i+1)
            if(indexint.eq.24) then
           if(sigma(indexint,i)*alfa(indexint,i)*coltemp(i).lt.60.) then
               jfotsout(indexint,2,i) = aux1(nz-i+1)
     $            *exp(-sigma(indexint,i)*alfa(indexint,i)*coltemp(i))
               else
                  jfotsout(indexint,2,i)=0.
               end if
            end if
         enddo
      end do


c************************************************
c      co2,119.6,202.5
c************************************************

      do indexint=25,31
         do i=1,nz
            aux2(nz-i+1) = co2colx(i) + o2colx(i)
         end do
         do i=1,nz2
            aux3(i) = jabsifotsint(indexint,1,nz2-i+1)
            aux4(i) = c12(nz2-i+1)
         enddo
         call interpfast 
     $        ( aux1, aux2, nz,  aux3, aux4, nz2 , limdown , limup )
         do i=1,nz
           if(sigma(indexint,i)*alfa(indexint,i)*coltemp(i).lt.60.) then
               jfotsout(indexint,1,i) = aux1(nz-i+1)
     $            *exp(-sigma(indexint,i)*alfa(indexint,i)*coltemp(i))*
     $           (1+alfa(indexint,i)*(t2(i)-t0(i)))
            else
               jfotsout(indexint,1,i) = 0.
            end if
         enddo
      end do


c************************************************
c      o2,119.6,202.5
c************************************************

      do indexint=25,31
         do i=1,nz
            aux2(nz-i+1) = co2colx(i) + o2colx(i)
         end do
         do i=1,nz2
            aux3(i) = jabsifotsint(indexint,2,nz2-i+1)
            aux4(i) = c12(nz2-i+1)
         enddo
         call interpfast 
     $        ( aux1, aux2, nz,  aux3, aux4, nz2 , limdown , limup )
         do i=1,nz
           if(sigma(indexint,i)*alfa(indexint,i)*coltemp(i).lt.60.) then
            jfotsout(indexint,2,i) = aux1(nz-i+1)
     $            *exp(-sigma(indexint,i)*alfa(indexint,i)*coltemp(i))
            else
               jfotsout(indexint,2,i)=0.
            end if
         enddo
      end do


c************************************************
c      h2o,119.6,202.5
c************************************************

      do indexint=25,31
         do i=1,nz
            aux2(nz-i+1) = co2colx(i) + o2colx(i)
         end do
         do i=1,nz2
            aux3(i) = jabsifotsint(indexint,4,nz2-i+1)
            aux4(i) = c12(nz2-i+1)
         enddo
         call interpfast 
     $        ( aux1, aux2, nz,  aux3, aux4, nz2 , limdown , limup )
         do i=1,nz
           if(sigma(indexint,i)*alfa(indexint,i)*coltemp(i).lt.60.) then
            jfotsout(indexint,4,i) = aux1(nz-i+1) 
     $            *exp(-sigma(indexint,i)*alfa(indexint,i)*coltemp(i))
            else
               jfotsout(indexint,4,i) = 0.
            end if
         enddo
      end do


c************************************************
c      h2o2,119.6,202.5
c************************************************

      do indexint=25,31
         do i=1,nz
            aux2(nz-i+1) = co2colx(i) + o2colx(i)
         end do
         do i=1,nz2
            aux3(i) = jabsifotsint(indexint,6,nz2-i+1)
            aux4(i) = c12(nz2-i+1)
         enddo
         call interpfast 
     $        ( aux1, aux2, nz,  aux3, aux4, nz2 , limdown , limup )
         do i=1,nz
           if(sigma(indexint,i)*alfa(indexint,i)*coltemp(i).lt.60.) then
            jfotsout(indexint,6,i) = aux1(nz-i+1) 
     $            *exp(-sigma(indexint,i)*alfa(indexint,i)*coltemp(i))
            else
               jfotsout(indexint,6,i) = 0.
            end if
         enddo
      end do


c************************************************
c     co2,202.6,210.0
c************************************************

      indexint=32
         do i=1,nz
            aux2(nz-i+1) = co2colx(i)
         end do
         do i=1,nz2
            aux3(i) = jabsifotsint(indexint,1,nz2-i+1)
            aux4(i) = c1(nz2-i+1)
         enddo
         call interpfast 
     $        ( aux1, aux2, nz,  aux3, aux4, nz2 , limdown , limup )
         do i=1,nz
           if(sigma(indexint,i)*alfa(indexint,i)*coltemp(i).lt.60.) then
               jfotsout(indexint,1,i) = aux1(nz-i+1)
     $            *exp(-sigma(indexint,i)*alfa(indexint,i)*coltemp(i))*
     $           (1+alfa(indexint,i)*(t2(i)-t0(i)))
            else
               jfotsout(indexint,1,i)=0.
            end if
         enddo


c************************************************
c     h2o2,202.6,210.0
c************************************************

      indexint=32
         do i=1,nz
            aux2(nz-i+1) = co2colx(i)
         end do
         do i=1,nz2
            aux3(i) = jabsifotsint(indexint,6,nz2-i+1)
            aux4(i) = c1(nz2-i+1)
         enddo
         call interpfast 
     $        ( aux1, aux2, nz,  aux3, aux4, nz2 , limdown , limup )
         do i=1,nz
           if(sigma(indexint,i)*alfa(indexint,i)*coltemp(i).lt.60.) then
            jfotsout(indexint,6,i) = aux1(nz-i+1) 
     $            *exp(-sigma(indexint,i)*alfa(indexint,i)*coltemp(i))
            else
               jfotsout(indexint,6,i)=0.
            end if
         enddo


c************************************************
c     h2o2,210.1,337.5
c************************************************

      indexint=33
      limdown=1.e-20
      limup=1e25
         do i=1,nz
            aux2(nz-i+1) = h2o2colx(i)
         end do
         do i=1,nz2
            aux3(i) = jabsifotsint(indexint,6,nz2-i+1)
            aux4(i) = ch2o2(nz2-i+1)
         enddo
         call interpfast 
     $        ( aux1, aux2, nz,  aux3, aux4, nz2 , limdown , limup )
         do i=1,nz
            jfotsout(indexint,6,i) = aux1(nz-i+1)
         enddo

      return

      end






















