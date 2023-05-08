c**********************************************************************
      
      subroutine paramfoto
     $(lswitch,zdens,tx,timestep,zenit,nz,co2x,o2x,o3px,cox,hx,ohx,
     & ho2x,h2x,h2ox,h2o2x,o1dx)

c     oct 2002      fgg
c**********************************************************************

      implicit none

      include 'param.h'
      include 'param_v3.h'

c     arguments

      integer nz,lswitch
      real zdens(nz)
      real tx(nz)
      real co2x(nz),o2x(nz),o3px(nz),cox(nz),hx(nz),ohx(nz)
      real    h2x(nz),h2ox(nz),h2o2x(nz),o1dx(nz),ho2x(nz)
      real    zenit
      real    timestep

c     local variables

      real*8    dt,deltat
      integer nt, ntmax
      real*8    tauco2(nzmax,nreact)
      real*8    tauo2(nzmax,nreact)
      real*8    tauo3p(nzmax,nreact)
      real*8    tauco(nzmax,nreact)
      real*8    tauh(nzmax,nreact)
      real*8    tauoh(nzmax,nreact)
      real*8    tauho2(nzmax,nreact)
      real*8    tauh2(nzmax,nreact)
      real*8    tauh2o(nzmax,nreact)
      real*8    tauo1d(nzmax,nreact)
      real*8    tauh2o2(nzmax,nreact)
      real*8    co2x8,o2x8,o3px8,cox8,hx8,ohx8,ho2x8,h2x8,h2ox8
      real*8    h2o2x8,o1dx8
      real*8    o1dx8aux,ohx8aux,ho2x8aux,dif_o1d,dif_oh,dif_ho2
      real*8    timestep8
      real*8    jdistot8(nabs,nzmax)
      real*8    jdistot8_b(2,nzmax)
      real*8    tx8

      real*8    tmin(nzmax)
      integer i,j,k,nada,n

      real*8    auxp,auxl
      
      character*1 o1d_eq(nzmax)
      character*1 oh_eq(nzmax)
      character*1 ho2_eq(nzmax)
      character*1 dn

      logical firstcall
      save firstcall
      data firstcall /.true./

      if (firstcall) then
        firstcall= .false.
        print*,'FOTO'
      endif

c 123  format(1x,f10.7,11(tr3,e12.7))

c**********************************************************************

C     Start: altitude loop
      timestep8=dble(timestep)

      do i=nz,lswitch,-1
         

c     Concentrations to real*8

         co2x8=dble(co2x(i))
         o2x8=dble(o2x(i))
         o3px8=dble(o3px(i))
         cox8=dble(cox(i))
         hx8=dble(hx(i))
         ohx8=dble(ohx(i))
         ho2x8=dble(ho2x(i))
         h2x8=dble(h2x(i))
         h2ox8=dble(h2ox(i))
         h2o2x8=dble(h2o2x(i))
         o1dx8=dble(o1dx(i))
         tx8=dble(tx(i))


c     reaction rates         

         call getch(tx8)

c         print*,i,co2x(i),o2x(i),o3px(i),cox(i),hx(i),ohx(i),ho2x(i),
c     &            h2x(i),h2ox(i),h2o2x(i),o1dx(i),tx(i),zdens(i)


c     photodissotiation rates

         jdistot(1,i) = 0.
         jdistot(2,i) = 0.
         jdistot(3,i) = 0.
         jdistot(4,i) = 0.
         jdistot(5,i) = 0.
         jdistot(6,i) = 0.
         jdistot_b(1,i) = 0.
         jdistot_b(2,i) = 0.

         call phdisrate(zenit,i)

         do j=1,6
            jdistot8(j,i) = dble(jdistot(j,i))
         end do
         jdistot8_b(1,i) = dble(jdistot_b(1,i))
         jdistot8_b(2,i) = dble(jdistot_b(2,i))
         
         
c     Lifetimes calculation

         do j=1,nreact
            tauco2(i,j) = 1.d30
            tauo2(i,j) = 1.d30
            tauo3p(i,j) = 1.d30
            tauco(i,j) = 1.d30
            tauh(i,j) = 1.d30
            tauoh(i,j) = 1.d30
            tauho2(i,j) = 1.d30
            tauh2(i,j) = 1.d30
            tauh2o(i,j) = 1.d30
            tauo1d(i,j) = 1.d30
            tauh2o2(i,j) = 1.d30
         end do

         if(jdistot8(1,i).gt.1d-30) tauco2(i,1) = 1.d0 / jdistot8(1,i)

         if(ch2*o2x8*co2x8.gt.1d-30) 
     $      tauh(i,2) = 1.d0 / (ch2 * o2x8 * co2x8)  
         if(ch2*hx8*co2x8.gt.1d-30)
     $      tauo2(i,2) = 1.d0 / (ch2 * hx8 * co2x8)
 
         if(ch3*o3px8.gt.1d-30) tauho2(i,3) = 1.d0 / (ch3 * o3px8)
         if(ch3*ho2x8.gt.1d-30) tauo3p(i,3) = 1.d0 / (ch3 * ho2x8)
 
         if(ch4*cox8.gt.1d-30) tauoh(i,4) = 1.d0 / (ch4 * cox8)
         if(ch4*ohx8.gt.1d-30) tauco(i,4) = 1.d0 / (ch4 * ohx8)
 
         if(ch5*ho2x8.gt.1d-30)tauho2(i,5) = 1.d0 / (2.d0 * ch5 * ho2x8)

         if(jdistot8(6,i).gt.1d-30) tauh2o2(i,6) = 1.d0 / jdistot8(6,i)

         if(ch7*ohx8.gt.1d-30) tauho2(i,7) = 1.d0 / (ch7 * ohx8)
         if(ch7*ho2x8.gt.1d-30) tauoh(i,7) = 1.d0 / (ch7 * ho2x8)
 
         if(jdistot8(4,i).gt.1d-30) tauh2o(i,8) = 1.d0 / jdistot8(4,i)

         if(ch9*o1dx8.gt.1d-30) tauh2o(i,9) = 1.d0 / (ch9 * o1dx8)
         if(ch9*h2ox8.gt.1d-30) tauo1d(i,9) = 1.d0 / (ch9 * h2ox8)

         if(ch10*o3px8*co2x8.gt.1d-30) 
     $     tauo3p(i,10) = 1.d0 / (2.d0 * ch10 * o3px8 * co2x8)

         if(ch11*o3px8.gt.1d-30) tauoh(i,11) = 1.d0 / (ch11 * o3px8)
         if(ch11*ohx8.gt.1d-30) tauo3p(i,11) = 1.d0 / (ch11 * ohx8)

         if(jdistot8(2,i).gt.1d-30) tauo2(i,12) = 1.d0 / jdistot8(2,i)

         if(ch13*hx8.gt.1d-30) tauho2(i,13) = 1.d0 / (ch13 * hx8)
         if(ch13*ho2x8.gt.1d-30) tauh(i,13) = 1.d0 / (ch13 * ho2x8)

         if(ch14*o1dx8.gt.1d-30) tauh2(i,14) = 1.d0 / (ch14 * o1dx8)
         if(ch14*h2x8.gt.1d-30) tauo1d(i,14) = 1.d0 / (ch14 * h2x8)

         if(ch15*ohx8.gt.1d-30) tauh2(i,15) = 1.d0 / (ch15 * ohx8)
         if(ch15*h2x8.gt.1d-30) tauoh(i,15) = 1.d0 / (ch15 * h2x8)

         if(jdistot8_b(1,i).gt.1d-30) tauco2(i,16) =1.d0/jdistot8_b(1,i)

         if(jdistot8_b(2,i).gt.1d-30) tauo2(i,17) =1.d0/jdistot8_b(2,i)

         if(ch18*ohx8.gt.1d-30) tauh2o2(i,18) = 1.d0 / (ch18 * ohx8)
         if(ch18*h2o2x8.gt.1d-30) tauoh(i,18) = 1.d0 / (ch18 * h2o2x8)

         if(ch19*co2x8.gt.1d-30) tauo1d(i,19) = 1.d0 / (ch19 * co2x8)

         if(ch20*o2x8.gt.1d-30) tauo1d(i,20) = 1.d0 / (ch20 * o2x8)

         if(ch21*o2x8*co2x8.gt.1d-30) tauo3p(i,21) = 1.d0 / 
     $        (ch21 * o2x8 * co2x8)
         if(ch21*o3px8*co2x8.gt.1d-30) tauo2(i,21) = 1.d0 / 
     $        (ch21 * o3px8 * co2x8)

         if(jdistot8(5,i).gt.1d-30) tauh2(i,22) = 1.d0 / jdistot8(5,i)

         !minimum lifetime
         tminco2(i) = 1.d30
         tmino2(i) = 1.d30
         tmino3p(i) = 1.d30
         tminco(i) = 1.d30
         tminh(i) = 1.d30
         tminoh(i) = 1.d30
         tminho2(i) = 1.d30
         tminh2(i) = 1.d30
         tminh2o(i) = 1.d30
         tmino1d(i) = 1.d30
         tminh2o2(i) = 1.d30
         
         do j=1,nreact
            tminco2(i) = min(tminco2(i),tauco2(i,j))
            tmino2(i) = min(tmino2(i),tauo2(i,j))
            tmino3p(i) = min(tmino3p(i),tauo3p(i,j))
            tminco(i) = min(tminco(i),tauco(i,j))
            tminh(i) = min(tminh(i),tauh(i,j))
            tminoh(i) = min(tminoh(i),tauoh(i,j))
            tminho2(i) = min(tminho2(i),tauho2(i,j))
            tminh2(i) = min(tminh2(i),tauh2(i,j))
            tminh2o(i) = min(tminh2o(i),tauh2o(i,j))
            tmino1d(i) = min(tmino1d(i),tauo1d(i,j))
            tminh2o2(i) = min(tminh2o2(i),tauh2o2(i,j))
         end do

         tmin(i) = 1.d30
         tmin(i) = min(tmin(i),tminco2(i))
         tmin(i) = min(tmin(i),tmino2(i))
         tmin(i) = min(tmin(i),tmino3p(i))
         tmin(i) = min(tmin(i),tminco(i))
         tmin(i) = min(tmin(i),tminh(i))
         tmin(i) = min(tmin(i),tminh2(i))
         tmin(i) = min(tmin(i),tminh2o(i))
         tmin(i) = min(tmin(i),tminh2o2(i))


         !Internal timestep
         if(tmin(i)/10.d0.gt.timestep8*3600.d0) then
            deltat = timestep8 * 3600.d0
         else
            deltat = tmin(i)/10.d0
         end if
         if(deltat.ne.timestep8*3600.d0) then
            n=int(timestep8*3600.d0/deltat)+1
            deltat=timestep8*3600.d0/n
         endif


         !Photochemical equilibrium for O1D,OH and HO2?
         if(tmino1d(i).lt.deltat/5.d0) then
            o1d_eq(i)='Y'
         else 
            o1d_eq(i)='N'
c            o1d_eq(i)='Y'
         endif
         if(tminoh(i).lt.deltat/5.d0) then
           oh_eq(i)='Y'
         else
           oh_eq(i)='N'
c           oh_eq(i)='Y'
         end if
         if(tminho2(i).lt.deltat/5.d0) then
            ho2_eq(i)='Y'
         else
            ho2_eq(i)='N'
c            ho2_eq(i)='Y'
         endif


C     Start: temporal loop
c        do dt=deltat,timestep8*3600.d0,deltat
         do nt=1,int(timestep8*3600.d0/deltat)
            dt=nt*deltat

            !Productions and losses
            do j=1,nreact
               Lco2(i,j) = 0.d0
               Pco2(i,j) = 0.d0
               Lo2(i,j) = 0.d0
               Po2(i,j) = 0.d0
               Lo3p(i,j) = 0.d0
               Po3p(i,j) = 0.d0
               Lco(i,j) = 0.d0
               Pco(i,j) = 0.d0
               Ph(i,j) = 0.d0
               Lh(i,j) = 0.d0
               Poh(i,j) = 0.d0
               Loh(i,j) = 0.d0
               Pho2(i,j) = 0.d0
               Lho2(i,j) = 0.d0
               Ph2(i,j) = 0.d0
               Lh2(i,j) = 0.d0
               Ph2o(i,j) = 0.d0
               Lh2o(i,j) = 0.d0
               Po1d(i,j) = 0.d0
               Lo1d(i,j) = 0.d0
               Ph2o2(i,j) = 0.d0
               Lh2o2(i,j) = 0.d0
            end do
            Pco2tot(i) = 0.d0
            Lco2tot(i) = 0.d0
            Po2tot(i) = 0.d0
            Lo2tot(i) = 0.d0
            Po3ptot(i) = 0.d0
            Lo3ptot(i) = 0.d0
            Pcotot(i) = 0.d0
            Lcotot(i) = 0.d0
            Phtot(i) = 0.d0
            Lhtot(i) = 0.d0
            Pohtot(i) = 0.d0
            Lohtot(i) = 0.d0
            Pho2tot(i) = 0.d0
            Lho2tot(i) = 0.d0
            Ph2tot(i) = 0.d0
            Lh2tot(i) = 0.d0
            Ph2otot(i) = 0.d0
            Lh2otot(i) = 0.d0
            Po1dtot(i) = 0.d0
            Lo1dtot(i) = 0.d0
            Ph2o2tot(i) = 0.d0
            Lh2o2tot(i) = 0.d0


c           R1: CO2 + hv -> CO + O

            Lco2(i,1) = jdistot8(1,i)
            Pco(i,1) = co2x8 * jdistot8(1,i)
            Po3p(i,1) = co2x8 * jdistot8(1,i)


c           R16(1b): CO2 + hv -> CO + O1D

            Lco2(i,16) = jdistot8_b(1,i)
            Pco(i,16) = co2x8 * jdistot8_b(1,i)
            Po1d(i,16) = co2x8 * jdistot8_b(1,i)


c           R2: H + O2 + CO2 -> HO2 + CO2

            Lh(i,2) = ch2 * o2x8 * co2x8
            Lo2(i,2) = ch2 * hx8 * co2x8
            Pho2(i,2) = ch2 * hx8 * o2x8 * co2x8


c           R3: O + HO2 -> OH + O2

            Lo3p(i,3) = ch3 * ho2x8 
            Lho2(i,3) = ch3 * o3px8
            Poh(i,3) = ch3 * o3px8 * ho2x8
            Po2(i,3) = ch3 * o3px8 * ho2x8


c           R4: CO + OH -> CO2 + H

            Lco(i,4) = ch4 * ohx8
            Loh(i,4) = ch4 * cox8
            Pco2(i,4) = ch4 * cox8 * ohx8
            Ph(i,4) = ch4 * cox8 * ohx8


c           R5: 2HO2 -> H2O2 + O2

            Lho2(i,5) = 2.d0 * ch5 * ho2x8 
            Po2(i,5) = ch5 * ho2x8 * ho2x8 
            Ph2o2(i,5) = ch5 * ho2x8 * ho2x8


c           R6: H2O2 + hv -> 2OH

            Lh2o2(i,6) = jdistot8(6,i) 
            Poh(i,6) = jdistot8(6,i) * h2o2x8


c           R7: OH + HO2 -> H2O + O2

            Loh(i,7) = ch7 * ho2x8 
            Lho2(i,7) = ch7 * ohx8
            Po2(i,7) = ch7 * ohx8 * ho2x8
            Ph2o(i,7) = ch7 * ohx8 * ho2x8


c           R8: H20 + hv -> H + OH

            Lh2o(i,8) = jdistot8(4,i)
            Ph(i,8) = jdistot8(4,i) * h2ox8
            Poh(i,8) = jdistot8(4,i) * h2ox8


c           R9: O(1D) + H2O -> 2OH

            Lo1d(i,9) = ch9 * h2ox8
            Lh2o(i,9) = ch9 * o1dx8
            Poh(i,9) = 2.d0 * ch9 * o1dx8 * h2ox8 
         

c           R10: 2O + CO2 -> O2 + CO2

            Lo3p(i,10) = 2.d0 * ch10 * o3px8 * co2x8  
            Po2(i,10) = ch10 * o3px8 * o3px8 * co2x8


c           R11: O + OH -> O2 + H

            Lo3p(i,11) = ch11 * ohx8
            Loh(i,11) = ch11 * o3px8
            Po2(i,11) = ch11 * o3px8 * ohx8
            Ph(i,11) = ch11 * o3px8 * ohx8


c           R12: O2 + hv -> 2O

            Lo2(i,12) = jdistot8(2,i)
            Po3p(i,12) = 2.d0 * jdistot8(2,i) * o2x8


c           R17(12b): O2 + hv -> O + O1D

            Lo2(i,17) = jdistot8_b(2,i)
            Po3p(i,17) = jdistot8_b(2,i) * o2x8
            Po1d(i,17) = jdistot8_b(2,i) * o2x8


c           R13: H + HO2 -> H2 + O2

            Lh(i,13) = ch13 * ho2x8 
            Lho2(i,13) = ch13 * hx8
            Ph2(i,13) = ch13 * hx8 * ho2x8
            Po2(i,13) = ch13 * hx8 * ho2x8


c           R14: O(1D) + H2 -> H + OH

            Lo1d(i,14) = ch14 * h2x8
            Lh2(i,14) = ch14 * o1dx8
            Ph(i,14) = ch14 * o1dx8 * h2x8
            Poh(i,14) = ch14 * o1dx8 * h2x8


c           R15: OH + H2 -> H + H20

            Loh(i,15) = ch15 * h2x8
            Lh2(i,15) = ch15 * ohx8
            Ph(i,15) = ch15 * ohx8 * h2x8
            Ph2o(i,15) = ch15 * ohx8 * h2x8


c           R18: OH + H2O2 -> H2O + HO2

            Loh(i,18) = ch18 * h2o2x8
            Lh2o2(i,18) = ch18 * ohx8
            Ph2o(i,18) = ch18 * ohx8 * h2o2x8
            Pho2(i,18) = ch18 * ohx8 * h2o2x8


c           R19: O(1D) + CO2 -> O + CO2

            Lo1d(i,19) = ch19 * co2x8
            Po3p(i,19) = ch19 * o1dx8 * co2x8


c           R20: O(1D) + O2 -> O + O2

            Lo1d(i,20) = ch20 * o2x8
            Po3p(i,20) = ch20 * o1dx8 * o2x8


c           R21: O + O2 + CO2 -> O3 + CO2

            Lo3p(i,21) = ch21 * o2x8 * co2x8
            Lo2(i,21) = ch21 * o3px8 * co2x8


c           R22: H2 + hv -> 2H

            Lh2(i,22) = jdistot8(5,i)
            Ph(i,22) = 2.d0 * jdistot8(5,i) * h2x8 


            do j=1,nreact
               Pco2tot(i) = Pco2tot(i) + Pco2(i,j)
               Lco2tot(i) = Lco2tot(i) + Lco2(i,j)
               Po2tot(i) = Po2tot(i) + Po2(i,j)
               Lo2tot(i) = Lo2tot(i) + Lo2(i,j)
               Po3ptot(i) = Po3ptot(i) + Po3p(i,j)
               Lo3ptot(i) = Lo3ptot(i) + Lo3p(i,j)
               Pcotot(i) = Pcotot(i) + Pco(i,j)
               Lcotot(i) = Lcotot(i) + Lco(i,j)
               Phtot(i) = Phtot(i) + Ph(i,j)
               Lhtot(i) = Lhtot(i) + Lh(i,j)
               Pohtot(i) = Pohtot(i) + Poh(i,j)
               Lohtot(i) = Lohtot(i) + Loh(i,j)
               Pho2tot(i) = Pho2tot(i) + Pho2(i,j)
               Lho2tot(i) = Lho2tot(i) + Lho2(i,j)
               Ph2tot(i) = Ph2tot(i) + Ph2(i,j)
               Lh2tot(i) = Lh2tot(i) + Lh2(i,j)
               Ph2otot(i) = Ph2otot(i) + Ph2o(i,j)
               Lh2otot(i) = Lh2otot(i) + Lh2o(i,j)
               Po1dtot(i) = Po1dtot(i) + Po1d(i,j)
               Lo1dtot(i) = Lo1dtot(i) + Lo1d(i,j)
               Ph2o2tot(i) = Ph2o2tot(i) + Ph2o2(i,j)
               Lh2o2tot(i) = Lh2o2tot(i) + Lh2o2(i,j)
            end do


            !New concentrations, implicit scheme
            co2x8 = (co2x8 + Pco2tot(i) * deltat) /
     $                (1.d0 + Lco2tot(i) * deltat)
            o2x8 = (o2x8 + Po2tot(i) * deltat) /
     $               (1.d0 + Lo2tot(i) * deltat)
            o3px8 = (o3px8 + Po3ptot(i) * deltat) /
     $                (1.d0 + Lo3ptot(i) * deltat)
            cox8 = (cox8 + Pcotot(i) * deltat) /
     $               (1.d0 + Lcotot(i) * deltat)
            hx8 = (hx8 + Phtot(i) * deltat) /
     $              (1.d0 + Lhtot(i) * deltat)
             h2x8 = (h2x8 + Ph2tot(i) * deltat) /
     $               (1.d0 + Lh2tot(i) * deltat)
            h2ox8 = (h2ox8 + Ph2otot(i) * deltat) / 
     $                (1.d0 + Lh2otot(i) * deltat)
            h2o2x8 = (h2o2x8 + Ph2o2tot(i) * deltat) /
     $                (1.d0 + Lh2o2tot(i) * deltat)

            if(o1d_eq(i).eq.'Y') then
               auxp=jdistot8_b(1,i)*co2x8+jdistot8_b(2,i)*o2x8
               auxl=ch9*h2ox8+ch14*h2x8+ch19*co2x8+ch20*o2x8
               if(auxl.ne.0.) then
                  o1dx8=auxp/auxl
               else
                  o1dx8=o1dx8+Po1dtot(i)*deltat
               endif
            else
               o1dx8 = (o1dx8 + Po1dtot(i) * deltat) /
     $                   (1.d0 + Lo1dtot(i) * deltat)
            end if

            if(oh_eq(i).eq.'Y') then
               auxp=ch3*o3px8*ho2x8+jdistot(4,i)*h2ox8+
     $              2.d0*ch9*o1dx8*h2ox8+ch14*o1dx8*h2x8+
     $              2.d0*jdistot(6,i)*h2o2x8
               auxl=ch4*cox8+ch7*ho2x8+ch11*o3px8+ch15*h2x8+
     $              ch18*h2o2x8
               if(auxl.ne.0.) then
                  ohx8=auxp/auxl
               else
                  ohx8=ohx8+Pohtot(i)*deltat
               end if
            else
               ohx8 = (ohx8 + Pohtot(i) * deltat) /
     $                  (1.d0 + Lohtot(i) * deltat)
            end if

            if(ho2_eq(i).eq.'Y') then
               auxp=ch2*hx8*o2x8*co2x8+ch18*ohx8*h2o2x8
               auxl=ch3*o3px8+2.d0*ch5*ho2x8+ch7*ohx8+ch13*hx8
               if(auxl.ne.0.) then
                  ho2x8=auxp/auxl
               else
                  ho2x8=ho2x8+Pho2tot(i)*deltat
               end if
            else
               ho2x8 = (ho2x8 + Pho2tot(i) * deltat) /
     $                  (1.d0 + Lho2tot(i) * deltat)
            end if


C        End: temporal loop
         end do
         
         
         co2x(i)  = real(co2x8)
         o2x(i)   = real(o2x8)
         o3px(i)  = real(o3px8)
         cox(i)   = real(cox8)
         hx(i)    = real(hx8)
         ohx(i)   = real(ohx8)
         ho2x(i)  = real(ho2x8)
         h2x(i)   = real(h2x8)
         h2ox(i)  = real(h2ox8)
         h2o2x(i) = real(h2o2x8)
         o1dx(i)  = real(o1dx8)

c         co2x(i)=(co2x8)
c         o2x(i)=(o2x8)
c         o3px(i)=(o3px8)
c         cox(i)=(cox8)
c         hx(i)=(hx8)
c         ohx(i)=(ohx8)
c         ho2x(i)=(ho2x8)
c         h2x(i)=(h2x8)
c         h2ox(i)=(h2ox8)
c         h2o2x(i)=(h2o2x8)
c         o1dx(i)=(o1dx8)

C     End: altitude loop
      end do

      return

      end


