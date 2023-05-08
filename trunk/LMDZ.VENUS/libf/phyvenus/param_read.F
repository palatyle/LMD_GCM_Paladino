      subroutine param_read

      use dimphy
      use conc
      implicit none


c     common variables and constants
      include 'param.h'
      include 'param_v4.h'
c      include 'datafile.h'
 
 
c     local variables

      integer    i,j,k,inter                          !indexes
      integer ierr
      real       nada
      character (len=100) :: datafile="HIGHATM"
      
c*************************+PROGRAM STARTS**************************


c     data for the UV heating tabulation

      data (crscabsi2(1,j),j=1,16) /5.61031E-19,1.59677E-18,4.7072E-18,
     $     1.48254e-17,2.07445e-17,2.573e-17,2.901e-17,3.083e-17,
     $     3.217e-17,3.539e-17,3.658e-17,3.63e-17,3.41239e-17,
     $     2.71019e-17,4.93677e-17,1.64e-17/

      data (crscabsi2(2,j),j=1,16) /0.27250E-18,0.11650E-17,0.39250E-17,
     $     0.10630E-16,0.15590E-16,0.17180E-16,0.19270E-16,0.22860E-16,
     $     0.24270E-16,0.24440E-16,0.25020E-16,0.26600E-16,0.25400E-16,
     $     0.35800E-16,0.25590E-16,0.16740E-16/

      data (crscabsi2(3,j),j=1,16) /0.2776E-18,0.9792E-18,0.3313E-17,
     $     0.6621E-17,0.8481E-17,0.9146E-17,0.9414E-17,0.1039E-16,
     $     0.1012E-16,0.1033E-16,0.1033E-16,0.1033E-16,0.8268E-17,
     $     0.6563E-17,0.3506E-17,0.3470E-17/

      data (crscabsi2(5,j),j=1,16) /.5E-20,.1077607E-19,.5670491E-19,
     $     .3322716E-18,.1054509E-17,.1700005E-17,.3171188E-17,
     $     .4734241E-17,.5108741E-17,.6022236E-17,.6741537E-17,
     $     .7277079E-17,.9070787E-17,.9708916E-17,.4026281E-17,0.0/

      data (crscabsi2(8,j),j=1,16) /0.0, 7.44175e-19, 2.23167e-18,
     $    8.46200e-18,1.18275e-17,1.54900e-17,2.32475e-17,2.41373e-17,
     $     2.55482e-17,2.38431e-17,2.28600e-17,2.35067e-17,2.56000e-17,
     $     2.64636e-17,2.86260e-17,3.26561e-17/

      data(crscabsi2(9,j),j=1,16) /3.48182e-20,3.37038e-19,1.03077e-18,
     $     4.01364e-18,6.45e-18,7.8e-18,1.0e-17,1.13500e-17,1.15500e-17,
     $     1.18000e-17,1.17500e-17,1.16000e-17,1.28667e-17,1.18500e-17,
     $     1.11000e-17,9.50000e-18/

      data(crscabsi2(10,j),j=1,16) /0.0,9.39833e-19,2.87714e-18,
     $     9.66900e-18,1.37063e-17,1.61820e-17,2.30450e-17,2.63373e-17,
     $     2.63773e-17,2.67677e-17,2.64100e-17,2.53000e-17,2.18100e-17,
     $     2.04941e-17,2.28160e-17,2.93550e-17/

      data(crscabsi2(11,j),j=1,16) /0.0,9.58555e-19,2.52767e-18,
     $     8.29700e-18,1.21850e-17,1.40500e-17,1.97025e-17,2.12018e-17,
     $     2.14673e-17,2.20331e-17,2.21500e-17,2.21600e-17,2.33200e-17,
     $     2.67800e-17,2.56400e-17,3.58561e-17/

      data(crscabsi2(12,j),j=1,16) /0.0,1.0e-20,2.5e-20,1.30400e-19,
     $     2.93800e-19,4.36000e-19,8.49400e-19,1.29400e-18,1.40500e-18,
     $     1.67600e-18,1.93400e-18,2.12200e-18,2.75800e-18,3.48400e-18,
     $     4.17200e-18,5.26000e-18/

      data(crscabsi2(13,j),j=1,16) /0.0,1.60e-18,4.99111e-18,1.48496e-17
     $     ,2.17395e-17,2.55857e-17,2.87754e-17,3.65571e-17,3.85691e-17,
     $     4.16286e-17,4.15117e-17,4.05901e-17,3.64000e-17,2.99670e-17,
     $     2.46796e-17,2.51789e-17/

      data freccen /3.4,7.5,14.5,23.0,30.3,34.1,49.6,50.5,52.5,56.0,
     $59.0,61.5,68.7,73.1,78.4,83.1,92.4,97.5,99.3,100.1,100.7,102.1,
     $104.5,116.8,121.3,127.0,130.6,153.7,162.8,171.4
     $,195.6,206.3,222.0,236.0,289.0,600./

      data co2crsc195/2.05864e-17,5.90557e-20,3.1027e-19,6.70653e-19,
     $4.55132e-19,8.87122e-20,1.32138e-20,7.22244e-23,2.88002e-26/

      data co2crsc295/2.05897e-17,6.71104e-20,3.45509e-19,7.45711e-19,
     $4.82752e-19,1.11594e-19,1.98308e-20,1.3853e-22,2.1414e-25/

c     Reads tabulated functions

      !Tabulated column amount
      open(210, status = 'old',
c    $file=trim(datafile)//'/EUVDAT/coln.dat',iostat=ierr)
     $file=trim(datafile)//'/EUVDAT/param_v5/coln.dat',iostat=ierr)

      IF (ierr.NE.0) THEN 
       write(*,*)'cant find directory EUVDAT containing param_v5 subdir'
       write(*,*)'(in aeronomars/param_read.F)'
       write(*,*)'It should be in :', trim(datafile),'/'
       write(*,*)'1) You can change this directory address in '
       write(*,*)'   callphys.def with datadir=/path/to/dir'
       write(*,*)'2) If necessary, EUVDAT (and other datafiles)'
       write(*,*)'   can be obtained online on:'
       write(*,*)'   http://www.lmd.jussieu.fr/~forget/datagcm/datafile'
       STOP
      ENDIF
 
      !Tabulated photoabsorption coefficients
      open(220,file=trim(datafile)//'/EUVDAT/param_v5/j2_an.dat')
      open(230,file=trim(datafile)//'/EUVDAT/param_v5/j3_an.dat')
      open(240,file=trim(datafile)//'/EUVDAT/param_v5/j1_an.dat')
      open(250,file=trim(datafile)//'/EUVDAT/param_v5/j2_bn.dat')
      open(260,file=trim(datafile)//'/EUVDAT/param_v5/j2_cn.dat')
      open(300,file=trim(datafile)//'/EUVDAT//param_v5/j2_dn.dat')
      open(270,file=trim(datafile)//'/EUVDAT//param_v5/j1_bn.dat')
      open(280,file=trim(datafile)//'/EUVDAT//param_v5/j1_cn.dat')
      open(290,file=trim(datafile)//'/EUVDAT//param_v5/j1_dn.dat')
      open(150,file=trim(datafile)//'/EUVDAT//param_v5/j4n.dat')
      open(160,file=trim(datafile)//'/EUVDAT//param_v5/j5n.dat')
      open(170,file=trim(datafile)//'/EUVDAT//param_v5/j6n.dat')
      open(180,file=trim(datafile)//'/EUVDAT//param_v5/j7n.dat')
      open(390,file=trim(datafile)//'/EUVDAT//param_v5/j8_an.dat')
      open(400,file=trim(datafile)//'/EUVDAT//param_v5/j8_bn.dat')
      open(410,file=trim(datafile)//'/EUVDAT//param_v5/j9n.dat')
      open(420,file=trim(datafile)//'/EUVDAT//param_v5/j10_an.dat')
      open(430,file=trim(datafile)//'/EUVDAT//param_v5/j10_bn.dat')
      open(440,file=trim(datafile)//'/EUVDAT//param_v5/j10_cn.dat')
      open(450,file=trim(datafile)//'/EUVDAT//param_v5/j11_an.dat')
      open(460,file=trim(datafile)//'/EUVDAT//param_v5/j11_bn.dat')
      open(470,file=trim(datafile)//'/EUVDAT//param_v5/j11_cn.dat')
      open(480,file=trim(datafile)//'/EUVDAT//param_v5/j12n.dat')
      open(490,file=trim(datafile)//'/EUVDAT//param_v5/j13_an.dat')
      open(500,file=trim(datafile)//'/EUVDAT//param_v5/j13_bn.dat')
      open(510,file=trim(datafile)//'/EUVDAT//param_v5/j13_cn.dat')

      
      do i=210,300,10
         read(i,*)
         read(i,*)
      end do

      do i=150,180,10
         read(i,*)
         read(i,*)
      end do

      do i=390,510,10
         read(i,*)
         read(i,*)
      enddo

      do i=nz2,1,-1
         read(210,*) (c1_16(i,j),j=1,16),c17_24(i),c25_29(i),c30_31(i),
     $        c32(i),c33(i),c34(i),c35(i),c36(i)
      end do

      do i=nz2,1,-1
         read(220,*) (jabsifotsintpar(i,2,j),j=1,16)
      end do
      
      do i=nz2,1,-1
         read(230,*) (jabsifotsintpar(i,3,j),j=1,16)
      end do

      do i=nz2,1,-1
         read(240,*) (jabsifotsintpar(i,1,j),j=1,16)
      end do

      do i=nz2,1,-1
         read(250,*) (jabsifotsintpar(i,2,j),j=17,24)
      end do


      do i=nz2,1,-1
         read(260,*) (jabsifotsintpar(i,2,j),j=25,31)
      end do

      do i=nz2,1,-1
         read(270,*) (jabsifotsintpar(i,1,j),j=17,24)
      end do

      do i=nz2,1,-1
         read(280,*) (jabsifotsintpar(i,1,j),j=25,31)
      end do

      do i=nz2,1,-1
         read(290,*) jabsifotsintpar(i,1,32)
      end do

      do i=nz2,1,-1
         read(300,*) (jabsifotsintpar(i,2,j),j=32,34)
      end do

      do i=nz2,1,-1
         read(160,*) (jabsifotsintpar(i,5,j),j=1,15)
      end do

      do i=nz2,1,-1
         read(150,*) (jabsifotsintpar(i,4,j),j=25,31)
      end do

      do i=nz2,1,-1
         read(170,*) (jabsifotsintpar(i,6,j),j=25,35)
      end do

      do i=nz2,1,-1
         read(180,*) (jabsifotsintpar(i,7,j),j=34,36)
      end do

      do i=nz2,1,-1
         read(390,*) (jabsifotsintpar(i,8,j),j=2,16)
      enddo

      do i=nz2,1,-1
         read(400,*) (jabsifotsintpar(i,8,j),j=17,24)
      enddo

      do i=nz2,1,-1
         read(410,*) (jabsifotsintpar(i,9,j),j=1,16)
      enddo

      do i=nz2,1,-1
         read(420,*) (jabsifotsintpar(i,10,j),j=2,16)
      enddo

      do i=nz2,1,-1
         read(430,*) (jabsifotsintpar(i,10,j),j=17,24)
      enddo

      do i=nz2,1,-1
         read(440,*) (jabsifotsintpar(i,10,j),j=25,32)
      enddo

      do i=nz2,1,-1
         read(450,*) (jabsifotsintpar(i,11,j),j=2,16)
      enddo

      do i=nz2,1,-1
         read(460,*) (jabsifotsintpar(i,11,j),j=17,24)
      enddo

      do i=nz2,1,-1
         read(470,*) (jabsifotsintpar(i,11,j),j=25,29)
      enddo
      
      do i=nz2,1,-1
         read(480,*) (jabsifotsintpar(i,12,j),j=2,16)
      enddo

      do i=nz2,1,-1
         read(490,*) (jabsifotsintpar(i,13,j),j=2,16)
      enddo
      
      do i=nz2,1,-1
         read(500,*) (jabsifotsintpar(i,13,j),j=17,24)
      enddo
      
      do i=nz2,1,-1
         read(510,*) (jabsifotsintpar(i,13,j),j=25,36)
      enddo

      do i=210,300,10
         close(i)
      end do

      do i=150,180,10
         close(i)
      end do

      do i=390,510,10
         close(i)
      enddo


c     set t0

      do i=1,nz2
         t0(i)=195.
      end do


      do i=1,ninter
         fluxtop(i)=1.
      end do

      !Parameters for the variation of the solar flux with 11 years cycle
      open(100,file=trim(datafile)//'/EUVDAT/param_v5/varflujo.dat')
      read(100,*)
      do i=1,24
         read(100,*) inter,ct1(i),p1(i),ct2(i),p2(i),nada
      end do
      close(100)

c     dissociation and ionization efficiencies

!      do inter=1,ninter
!         efdisco2(inter)=0.
!         efdiso2(inter)=0.
!         efdish2(inter)=0.
!         efdish2o(inter)=0.
!         efdish2o2(inter)=0.
!         efdiso3(inter)=0.
!         efdisco(inter)=0.
!         efdisn2(inter)=0.
!         efdisno(inter)=0.
!         efdisno2(inter)=0.
!         efionco2(inter,1)=0.
!         efionco2(inter,2)=0.
!         efionco2(inter,3)=0.
!         efionco2(inter,4)=0.
!         efiono3p(inter)=0.
!         efionn2(inter,1)=0.
!         efionn2(inter,2)=0.
!         efionco(inter,1)=0.
!         efionco(inter,2)=0.
!         efionco(inter,3)=0.
!         efionn(inter)=0.
!         efionh(inter)=0.
!         efionno(inter)=0.
!      enddo


c     CO2, O2, NO

!      open(120,file=trim(datafile)//'/EUVDAT/param_v5/efdis_inter.dat')
!      read(120,*)
!      do i=1,21
!         read(120,*)inter,efdisco2(inter),efdiso2(inter),efdisno(inter)
!      do inter=8,28
!         read(120,*)i,efdisco2(inter),efdiso2(inter),efdisno(inter)
!      enddo
!      do inter=29,ninter
!         efdisco2(inter)=1.
!         efdiso2(inter)=1.
!         efdisno(inter)=1.
!      enddo


c     N2

!      efdisn2(15)=0.1
!      do inter=16,ninter
!         efdisn2(inter)=1.
!      enddo


c     CO

!      efdisco(16)=0.5
!      do inter=17,ninter
!         efdisco(inter)=1.
!      enddo

      
c     O, N, H

!      do inter=1,ninter
!         efdiso(inter)=0.
!         efdisn(inter)=0.
!         efdish(inter)=0.
!      enddo


c     H2O, H2O2, O3, NO2

!      do inter=25,31
!         efdish2o(inter)=1.
!      enddo
!      do inter=25,35
!         efdish2o2(inter)=1.
!      enddo
!      do inter=34,36
!         efdiso3(inter)=1.
!      enddo
!      do inter=27,36
!         efdisno2(inter)=1.
!      enddo
!      do inter=1,15
!         efdish2(inter)=1.
!      enddo
         
      !4 possible channels for CO2 ionization
!      do inter=14,16
!         efionco2(inter,1)=1.-efdisco2(inter)
!      enddo
!      efionco2(13,1)=0.805*(1.-efdisco2(13))
!      efionco2(13,2)=0.195*(1.-efdisco2(13))
!      do inter=11,12
!         efionco2(inter,3)=1.-efdisco2(inter)
!      enddo
!      efionco2(10,3)=0.9*(1.-efdisco2(10))
!      efionco2(10,4)=0.1*(1.-efdisco2(10))
!      do inter=2,9
!         efionco2(inter,4)=1.-efdisco2(inter)
!      enddo

      !For O(3p) total ionization under 91.1 nm
!      do inter=1,16
!         efiono3p(inter)=1.d0
!      enddo

      !2 channels for N2 ionization
!      do inter=9,15
!         efionn2(inter,1)=1.-efdisn2(inter)
!      enddo
!      do inter=2,8
!         efionn2(inter,2)=1.-efdisn2(inter)
!      enddo
      
      !3 channels for CO ionization
!      do inter=11,16
!         efionco(inter,1)=1.-efdisco(inter)
!      enddo
!      efionco(10,1)=0.87*(1.-efdisco(10))
!      efionco(10,2)=0.13*(1.-efdisco(10))
!      do inter=8,9
!         efionco(inter,2)=1.-efdisco(inter)
!      enddo
!      efionco(7,2)=0.1*(1.-efdisco(7))
!      efionco(7,3)=0.9*(1.-efdisco(7))
!      do inter=2,6
!         efionco(inter,3)=1.-efdisco(inter)
!      enddo

      !Total ionization under 85 nm for N
!      do inter=1,16
!         efionn(inter)=1.
!      enddo

      !NO
!      do inter=2,28
!         efionno(inter)=1.-efdisno(inter)
!      enddo

      !Total ionization under 90 nm for H
!      do inter=3,16
!         efionh(inter)=1.
!      enddo


      return


      end

