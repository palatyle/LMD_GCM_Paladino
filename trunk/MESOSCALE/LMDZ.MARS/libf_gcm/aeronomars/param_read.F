      subroutine param_read

      implicit none


c     common variables and constants

#include "param.h"
#include "param_v3.h"
#include "datafile.h"


c     local variables

      integer    i,j,k,inter                          !indexes
      integer ierr,lnblnk 
      real       nada
      external lnblnk

     
c*************************+PROGRAM STARTS**************************


c     data

      data anchint/4.9,5.3,6.1,12.8,2.0,6.9,11.8,0.4,1.9,4.9,0.9,4.0,
     $6.9,4.9,4.9,10.2,5.6,1.9,0.9,0.9,0.9,0.9,3.7
     $,13.1,3.9,6.0,4.8,23.9,8.4,8.5,26.8,7.5,127.6/

      data (crscabsi2(1,j),j=1,ninter2) /0.0,4.42e-18,7.51e-18
     $,1.498e-17,2.344e-17,2.388e-17,2.927E-17,3.161e-17,3.161e-17,
     $3.4e-17,3.4e-17,2.588e-17,2.596e-17,2.248e-17,3.183e-17,
     $1.284e-17/

      data (crscabsi2(2,j),j=1,ninter2) /2.736e-19,1.055e-18,4.789e-18
     $,1.063e-17,1.559e-17,1.698e-17,2.164e-17,2.241e-17,2.427e-17,
     $2.444e-17,2.502e-17,2.516e-17,2.907e-17,3.649e-17,3.066e-17,
     $2.093e-17/

      data (crscabsi2(3,j),j=1,ninter2) /2.603e-19,6.985e-19,3.97e-18
     $,6.621e-18,8.433e-18,9.011e-18,9.862e-18,1.033e-17,1.012e-17,
     $1.033e-17,1.033e-17,1.034e-17,9.0e-18,3.564e-18,3.572e-18
     $,3.47E-18/

      data freccen /3.4,7.5,14.5,23.0,30.3,34.1,49.6,50.5,52.5,56.0,
     $59.0,61.5,68.7,73.1,78.4,83.1,92.4,97.5,99.3,100.1,100.7,102.1,
     $104.5,116.8,121.3,127.0,130.6,153.7,162.8,171.4
     $,195.6,206.3,273.5/

      data co2crsc195/3.691e-19,4.44216e-20,3.86945e-19,5.94208e-19,
     $2.93217e-19,7.58769e-20,8.60192e-21,4.20007e-24,2.29996e-26/

      data co2crsc295/3.691e-19,5.21572e-20,4.23488e-19,6.54728e-19,
     $3.30227e-19,1.03183e-19,1.55722e-20,1.72317e-23,7.0e-25/


c     Reads tabulated functions

      open(210, status = 'old',
     $file=datafile(1:lnblnk(datafile))//'/EUVDAT/coln.dat',iostat=ierr)

      IF (ierr.NE.0) THEN
       write(*,*)'cant find directory EUVDAT and content coln.dat'
       write(*,*)'(in aeronomars/param_read.F)'
       write(*,*)'It should be in :',datafile(1:lnblnk(datafile)),'/'
       write(*,*)'1) You can change this directory address in '
       write(*,*)'   file phymars/datafile.h'
       write(*,*)'2) If necessary, EUVDAT (and other datafiles)'
       write(*,*)'   can be obtained online on:'
       write(*,*)'   http://www.lmd.jussieu.fr/~forget/datagcm/datafile'
       STOP
      ENDIF

      open(220,file=datafile(1:lnblnk(datafile))//'/EUVDAT/j2_an.dat')
      open(230,file=datafile(1:lnblnk(datafile))//'/EUVDAT/j3_an.dat')
      open(240,file=datafile(1:lnblnk(datafile))//'/EUVDAT/j1_an.dat')
      open(250,file=datafile(1:lnblnk(datafile))//'/EUVDAT/j2_bn.dat')
      open(260,file=datafile(1:lnblnk(datafile))//'/EUVDAT/j2_cn.dat')
      open(270,file=datafile(1:lnblnk(datafile))//'/EUVDAT/j1_bn.dat')
      open(280,file=datafile(1:lnblnk(datafile))//'/EUVDAT/j1_cn.dat')
      open(290,file=datafile(1:lnblnk(datafile))//'/EUVDAT/j1_dn.dat')
      open(150,file=datafile(1:lnblnk(datafile))//'/EUVDAT/j4n.dat')
      open(160,file=datafile(1:lnblnk(datafile))//'/EUVDAT/j5n.dat')
      open(170,file=datafile(1:lnblnk(datafile))//'/EUVDAT/j6n.dat')

      do i=210,290,10
         read(i,*)
         read(i,*)
      end do

      do i=150,170,10
         read(i,*)
         read(i,*)
      end do

      do i=nz2,1,-1
        read(210,*) c23(i),(c123(j,i),j=2,ninter2),c12(i),c1(i),ch2o2(i)
      end do

      do i=nz2,1,-1
         read(220,*) (jabsifotsint(j,2,i),j=1,ninter2)
      end do
      
      do i=nz2,1,-1
         read(230,*) (jabsifotsint(j,3,i),j=1,ninter2)
      end do

      do i=nz2,1,-1
         read(240,*) (jabsifotsint(j,1,i),j=2,ninter2)
      end do

      do i=nz2,1,-1
         read(250,*) (jabsifotsint(j,2,i),j=17,24)
      end do

      do i=nz2,1,-1
         read(260,*) (jabsifotsint(j,2,i),j=25,31)
      end do

      do i=nz2,1,-1
         read(270,*) (jabsifotsint(j,1,i),j=17,24)
      end do

      do i=nz2,1,-1
         read(280,*) (jabsifotsint(j,1,i),j=25,31)
      end do

      do i=nz2,1,-1
         read(290,*) jabsifotsint(32,1,i)
      end do

      do i=nz2,1,-1
         read(160,*) (jabsifotsint(j,5,i),j=1,15)
      end do

      do i=nz2,1,-1
         read(150,*) (jabsifotsint(j,4,i),j=25,31)
      end do

      do i=nz2,1,-1
         read(170,*) (jabsifotsint(j,6,i),j=25,33)
      end do

      do i=210,290,10
         close(i)
      end do

      do i=150,170,10
         close(i)
      end do


c     reads t0

      open(120,file=datafile(1:lnblnk(datafile))//'/EUVDAT/t0.dat')
      do i=1,201
         read(120,*)t0(i)
      end do
      close(120)

      open(100,file=datafile(1:lnblnk(datafile))//'/EUVDAT/flujo.dat')
      do i=1,ninter
         read(100,*) inter,fluxtophr(i)
      end do
      close(100)


      open(99,file=datafile(1:lnblnk(datafile))//'/EUVDAT/varflujo.dat')
      read(99,*)
      do i=1,24
         read(99,*) inter,ct1(i),p1(i),ct2(i),p2(i),nada
      end do
      close(99)

c     eficiencias de disociacion (de Torr et al, 1979)

      do inter=1,11
         efdisco2(inter) = 0.
         efdiso2(inter) = 0.
      end do

      do inter=1,15
         efdish2(inter) = 1.
      end do

      efdisco2(12) = 0.183
      efdiso2(12) = 0.003

      efdisco2(13) = 0.163
      efdiso2(13) = 0.170

      efdisco2(14) = 0.243
      efdiso2(14) =0.180

      efdisco2(15) = 0.323
      efdiso2(15) =0.653

      efdisco2(16) = 0.235
      efdiso2(16) =0.616

      do inter=17,32
         efdisco2(inter) = 1.0
      end do

      efdiso2(17) = 0.399
      do inter=18,20
         efdiso2(inter) = 0.261
      end do

      do inter=21,23
         efdiso2(inter) = 0.755
      end do


      do inter=24,31
         efdiso2(inter) = 1.
      end do

      do inter=25,31
         efdish2o(inter) = 1.
      end do

      do inter=25,33
         efdish2o2(inter) = 1.
      end do


      return


      end

