      subroutine param_read_e107

      use param_v4_h, only: jfotsout,crscabsi2,
     .    c1_16,c17_24,c25_29,c30_31,c32,c33,c34,c35,c36,
     .    co2crsc195,co2crsc295,t0,
     .    jabsifotsintpar,ninter,nz2,
     .    nabs,e107,date_e107,e107_tab,
     .    coefit0,coefit1,coefit2,coefit3,coefit4,
     .    efdisco2,efdiso2,efdish2o,
     .    efdish2o2,efdish2,efdiso3,
     .    efdiso,efdisn,efdish,
     .    efdisno,efdisn2,efdisno2,
     .    efdisco,efionco2,efiono2,efionn2,
     .    efionco,efiono3p,efionn,
     .    efionno,efionh,
     .    fluxtop,ct1,ct2,p1,p2

      use datafile_mod, only: datadir

      implicit none


c     common variables and constants
      include "callkeys.h"
 
 
c     local variables

      integer    i,j,k,inter                          !indexes
      integer ierr,lnblnk 
      real       nada
      character*13  filename
  
      
c*************************+PROGRAM STARTS**************************


c     Reads tabulated functions

      !Tabulated column amount
      open(210, status = 'old',
c    $file=trim(datadir)//'/EUVDAT/coln.dat',iostat=ierr)
     $file=trim(datadir)//'/EUVDAT/param_v6/coln.dat',iostat=ierr)

      IF (ierr.NE.0) THEN 
       write(*,*)'cant find directory EUVDAT containing param_v6 subdir'
       write(*,*)'(in aeronomars/param_read.F)'
       write(*,*)'It should be in :', trim(datadir),'/'
       write(*,*)'1) You can change this directory address in '
       write(*,*)'   callphys.def with datadir=/path/to/dir'
       write(*,*)'2) If necessary, EUVDAT (and other datafiles)'
       write(*,*)'   can be obtained online on:'
       write(*,*)' http://www.lmd.jussieu.fr/~lmdz/planets/mars/datadir'
       STOP
      ENDIF
 
      !Tabulated photoabsorption coefficients
      open(220,file=trim(datadir)//'/EUVDAT/param_v6/trans2_an.dat')
      open(230,file=trim(datadir)//'/EUVDAT/param_v6/trans3_an.dat')
      open(240,file=trim(datadir)//'/EUVDAT/param_v6/trans1_an.dat')
      open(250,file=trim(datadir)//'/EUVDAT/param_v6/trans2_bn.dat')
      open(260,file=trim(datadir)//'/EUVDAT/param_v6/trans2_cn.dat')
      open(300,file=trim(datadir)//'/EUVDAT/param_v6/trans2_dn.dat')
      open(270,file=trim(datadir)//'/EUVDAT/param_v6/trans1_bn.dat')
      open(280,file=trim(datadir)//'/EUVDAT/param_v6/trans1_cn.dat')
      open(290,file=trim(datadir)//'/EUVDAT/param_v6/trans1_dn.dat')
      open(150,file=trim(datadir)//'/EUVDAT/param_v6/trans4n.dat')
      open(160,file=trim(datadir)//'/EUVDAT/param_v6/trans5n.dat')
      open(170,file=trim(datadir)//'/EUVDAT/param_v6/trans6n.dat')
      open(180,file=trim(datadir)//'/EUVDAT/param_v6/trans7n.dat')
      open(390,file=trim(datadir)//'/EUVDAT/param_v6/trans8_an.dat')
      open(400,file=trim(datadir)//'/EUVDAT/param_v6/trans8_bn.dat')
      open(410,file=trim(datadir)//'/EUVDAT/param_v6/trans9n.dat')
      open(420,file=trim(datadir)//'/EUVDAT/param_v6/trans10_an.dat')
      open(430,file=trim(datadir)//'/EUVDAT/param_v6/trans10_bn.dat')
      open(440,file=trim(datadir)//'/EUVDAT/param_v6/trans10_cn.dat')
      open(450,file=trim(datadir)//'/EUVDAT/param_v6/trans11_an.dat')
      open(460,file=trim(datadir)//'/EUVDAT/param_v6/trans11_bn.dat')
      open(470,file=trim(datadir)//'/EUVDAT/param_v6/trans11_cn.dat')
      open(480,file=trim(datadir)//'/EUVDAT/param_v6/trans12n.dat')
      open(490,file=trim(datadir)//'/EUVDAT/param_v6/trans13_an.dat')
      open(500,file=trim(datadir)//'/EUVDAT/param_v6/trans13_bn.dat')
      open(510,file=trim(datadir)//'/EUVDAT/param_v6/trans13_cn.dat')

      
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
      open(620,file=trim(datadir)//'/EUVDAT/param_v6/fit_js_e107.dat')
      do i=1,ninter
         read(620,*)
         do j=1,nabs
            read(620,*)coefit0(i,j),coefit1(i,j),coefit2(i,j),
     $           coefit3(i,j),coefit4(i,j)
         enddo
      enddo
      close(620)


      !Tabulated values for E10.7
      if(solvaryear.eq.23) then
         filename="e107_MY23.dat"
      else if(solvaryear.eq.24) then
         filename="e107_MY24.dat"
      else if(solvaryear.eq.25) then
         filename="e107_MY25.dat"
      else if (solvaryear.eq.26) then
         filename="e107_MY26.dat"
      else if (solvaryear.eq.27) then
         filename="e107_MY27.dat"
      else if (solvaryear.eq.28) then
         filename="e107_MY28.dat"
      else if (solvaryear.eq.29) then
         filename="e107_MY29.dat"
      else if(solvaryear.eq.30) then
         filename="e107_MY30.dat"
      else if(solvaryear.eq.31) then
         filename="e107_MY31.dat"
      else if(solvaryear.eq.32) then
         filename="e107_MY32.dat"
      else if(solvaryear.eq.33) then
         filename="e107_MY33.dat"
      else if(solvaryear.eq.34) then
         filename="e107_MY34.dat"
      else
         write(*,*)"param_read_e107: "
         write(*,*)"bad value for solvaryear in callphys.def"
         write(*,*)"solvaryear must be between 24 and 33"
         stop
      endif
      
      open(640,file=trim(datadir)//'/EUVDAT/param_v6/'//filename)
      read(640,*)
      do i=1,669
         read(640,*)date_e107(i),e107_tab(i)
         write(*,*)'param_read_e107/292',date_e107(i),e107_tab(i),i
      enddo
      close(640)
         

c     dissociation and ionization efficiencies

      do inter=1,ninter
         efdisco2(inter)=0.
         efdiso2(inter)=0.
         efdish2(inter)=0.
         efdish2o(inter)=0.
         efdish2o2(inter)=0.
         efdiso3(inter)=0.
         efdisco(inter)=0.
         efdisn2(inter)=0.
         efdisno(inter)=0.
         efdisno2(inter)=0.
         efionco2(inter,1)=0.
         efionco2(inter,2)=0.
         efionco2(inter,3)=0.
         efionco2(inter,4)=0.
         efiono2(inter,1)=0.
         efiono2(inter,2)=0.
         efiono3p(inter)=0.
         efionn2(inter,1)=0.
         efionn2(inter,2)=0.
         efionco(inter,1)=0.
         efionco(inter,2)=0.
         efionco(inter,3)=0.
         efionn(inter)=0.
         efionh(inter)=0.
         efionno(inter)=0.
      enddo


c     CO2, O2, NO

!      open(120,file=trim(datadir)//'/EUVDAT/param_v5/efdis_inter.dat')
!      read(120,*)
!!      do i=1,21
!!         read(120,*)inter,efdisco2(inter),efdiso2(inter),efdisno(inter)
!      do inter=8,28
!         read(120,*)i,efdisco2(inter),efdiso2(inter),efdisno(inter)
!      enddo
!      do inter=29,ninter
!         efdisco2(inter)=1.
!         efdiso2(inter)=1.
!         efdisno(inter)=1.
!      enddo
!      close(120)


c     N2

      efdisn2(15)=0.1
      do inter=16,ninter
         efdisn2(inter)=1.
      enddo


c     CO

      efdisco(16)=0.5
      do inter=17,ninter
         efdisco(inter)=1.
      enddo

      
c     O, N, H

      do inter=1,ninter
         efdiso(inter)=0.
         efdisn(inter)=0.
         efdish(inter)=0.
      enddo


c     H2O, H2O2, O3, NO2

      do inter=25,31
         efdish2o(inter)=1.
      enddo
      do inter=25,35
         efdish2o2(inter)=1.
      enddo
      do inter=34,36
         efdiso3(inter)=1.
      enddo
      do inter=27,36
         efdisno2(inter)=1.
      enddo
      do inter=1,15
         efdish2(inter)=1.
      enddo
         
      !4 possible channels for CO2 ionization
      open(130,file=trim(datadir)//'/EUVDAT'//
     $     '/co2ion_branchingratio_schunkandnagy2000_param.dat')
      do inter=1,16
         read(130,*)i,nada,efionco2(inter,1),efionco2(inter,2),
     $        efionco2(inter,3),efionco2(inter,4)
         !Multiply the relative efficiency of each channel by the total ionization efficiency (second column)
         efdisco2(inter)=1.-nada
         efionco2(inter,1)=(1.-efdisco2(inter))*efionco2(inter,1)
         efionco2(inter,2)=(1.-efdisco2(inter))*efionco2(inter,2)
         efionco2(inter,3)=(1.-efdisco2(inter))*efionco2(inter,3)
         efionco2(inter,4)=(1.-efdisco2(inter))*efionco2(inter,4)
      enddo
      close(130)
      do inter=17,36
         efdisco2(inter)=1.
      enddo

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


      !2 possible channels for O2 ionization
      open(131,file=trim(datadir)//'/EUVDAT'//
     $     '/o2ion_branchingratio_schunkandnagy2000_param.dat')
      do inter=1,23
         read(131,*)i,nada,efiono2(inter,1),efiono2(inter,2)
         !Multiply the relative efficiency of each channel by the total ionization efficiency (second column)
         efdiso2(inter)=1.-nada
         efiono2(inter,1)=(1.-efdiso2(inter))*efiono2(inter,1)
         efiono2(inter,2)=(1.-efdiso2(inter))*efiono2(inter,2)
      enddo
      close(131)
      do inter=24,36
         efdiso2(inter)=1.
      enddo


      !For O(3p) total ionization under 91.1 nm
      do inter=1,16
         efiono3p(inter)=1.d0
      enddo


      !2 channels for N2 ionization
      open(132,file=trim(datadir)//'/EUVDAT'//
     $     '/n2ion_branchingratio_schunkandnagy2000_param.dat')
      do inter=1,15
         read(132,*)i,nada,efionn2(inter,1),efionn2(inter,2)
         !Multiply the relative efficiency of each channel by the total ionization efficiency (second column)
         efdisn2(inter)=1.-nada
         efionn2(inter,1)=(1.-efdisn2(inter))*efionn2(inter,1)
         efionn2(inter,2)=(1.-efdisn2(inter))*efionn2(inter,2)
      enddo
      close(132)
      do inter=16,36
         efdisn2(inter)=1.
      enddo

!      do inter=9,15
!         efionn2(inter,1)=1.-efdisn2(inter)
!      enddo
!      do inter=2,8
!         efionn2(inter,2)=1.-efdisn2(inter)
!      enddo
 
     
      !3 channels for CO ionization
       open(133,file=trim(datadir)//'/EUVDAT'//
     $     '/coion_branchingratio_schunkandnagy2000_param.dat')
      do inter=1,16
         read(133,*)i,nada,efionco(inter,1),efionco(inter,2),
     $        efionco(inter,3)
         !Multiply the relative efficiency of each channel by the total ionization efficiency (second column)
         efdisco(inter)=1.-nada
         efionco(inter,1)=(1.-efdisco(inter))*efionco(inter,1)
         efionco(inter,2)=(1.-efdisco(inter))*efionco(inter,2)
         efionco(inter,3)=(1.-efdisco(inter))*efionco(inter,3)
      enddo
      close(133)
      do inter=17,36
         efdisco(inter)=1.
      enddo

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
      do inter=1,16
         efionn(inter)=1.
      enddo

      !NO
      do inter=2,28
         efionno(inter)=1.-efdisno(inter)
      enddo

      !Total ionization under 90 nm for H
      do inter=3,16
         efionh(inter)=1.
      enddo


      return


      end

