c**********************************************************************

      subroutine jthermcalc(ig,chemthermod,rm,nesptherm,tx,iz,zenit)


c     feb 2002        fgg           first version
c     nov 2002        fgg           second version
c     
c     mar 2014        gg            update for Venus GCM 
c
c modified from paramhr.F
c MAC July 2003
c**********************************************************************
      use dimphy
      use conc
c      use chemparam_mod
      implicit none

c     common variables and constants
      include "param.h"
      include "param_v4.h"

c     input and output variables

      integer    ig
      integer    chemthermod 
      integer    nesptherm                  !Number of species considered
      real       rm(klev,nesptherm)         !Densities (cm-3)
      real       tx(klev)                   !temperature
      real       zenit                      !SZA
      real       iz(klev)                   !Local altitude


c    local parameters and variables

      real       co2colx(klev)              !column density of CO2 (cm^-2)
      real       o3pcolx(klev)              !column density of O(3P)(cm^-2)
      real       n2colx(klev)               !N2 column density (cm-2)
      real       cocolx(klev)               !CO column density (cm-2)
c      real       o2colx(klev)               !column density of O2(cm^-2)
c      real       h2colx(klev)               !H2 column density (cm-2)
c      real       h2ocolx(klev)              !H2O column density (cm-2)
c      real       h2o2colx(klev)             !column density of H2O2(cm^-2)
c      real       o3colx(klev)               !O3 column density (cm-2)
c     real       hcolx(klev)                !H column density (cm-2)
c     real       no2colx(klev)              !NO2 column density (cm-2)
c     real       nocolx(klev)               !NO column density (cm-2)
      real       t2(klev)
      real       coltemp(klev)
      real       sigma(ninter,klev)
      real       alfa(ninter,klev)
    
      integer    i,j,k,indexint                 !indexes
      character  dn



c     variables used in interpolation

      real*8      auxcoltab(nz2)
      real*8      auxjco2(nz2)
c     real*8      auxjo2(nz2)
      real*8      auxjo3p(nz2)
c     real*8      auxjh2o(nz2)
c     real*8      auxjh2(nz2)
c     real*8      auxjh2o2(nz2)
c     real*8      auxjo3(nz2)
      real*8      auxjn2(nz2)
c      real*8      auxjn(nz2)
c     real*8      auxjno(nz2)
      real*8      auxjco(nz2)
c     real*8      auxjh(nz2)
c     real*8      auxjno2(nz2)
      real*8      wp(klev),wm(klev)
      real*8      auxcolinp(klev)
      integer     auxind(klev)
      integer     auxi
      integer     ind
      real*8      cortemp(klev)

      real*8      limdown                      !limits for interpolation
      real*8      limup                        !        ""

 
! Tracer indexes in the thermospheric chemistry:
!!! ATTENTION. These values have to be identical to those in euvheat.F90
!!! If the values are changed there, the same has to be done here  !!!

      integer,parameter :: ix_co2=1
      integer,parameter :: ix_n2=13
c      integer,parameter :: i_n=14
      integer,parameter :: ix_o=3
      integer,parameter :: ix_co=4


c*************************PROGRAM STARTS*******************************
      
      if(zenit.gt.140.) then
         dn='n'
         else
         dn='d'
      end if
      if(dn.eq.'n') then
        return
      endif
      
      !Initializing the photoabsorption coefficients
      jfotsout(:,:,:)=0.

      !Auxiliar temperature to take into account the temperature dependence
      !of CO2 cross section
      do i=1,klev
         t2(i)=tx(i)
         if(t2(i).lt.195.0) t2(i)=195.0
         if(t2(i).gt.295.0) t2(i)=295.0
      end do

      !Calculation of column amounts 
      call column(ig,chemthermod,rm,nesptherm,tx,iz,zenit,
     $     co2colx,o3pcolx, n2colx,cocolx)

      !Auxiliar column to include the temperature dependence 
      !of CO2 cross section
      coltemp(klev)=co2colx(klev)*abs(t2(klev)-t0(klev))
      do i=klev-1,1,-1
        coltemp(i)=!coltemp(i+1)+     PQ SE ELIMINA? REVISAR 
     $         ( rm(i,ix_co2) + rm(i+1,ix_co2) ) * 0.5 
     $         * 1e5 * (iz(i+1)-iz(i)) * abs(t2(i)-t0(i))
      end do
      
      !Calculation of CO2 cross section at temperature t0(i)
      do i=1,klev
         do indexint=24,32
           sigma(indexint,i)=co2crsc195(indexint-23)
           alfa(indexint,i)=((co2crsc295(indexint-23)
     $          /sigma(indexint,i))-1.)/(295.-t0(i))
        end do
      end do

! Interpolation to the tabulated photoabsorption coefficients for each species
! in each spectral interval


c     auxcolinp-> Actual atmospheric column
c     auxj*-> Tabulated photoabsorption coefficients
c     auxcoltab-> Tabulated atmospheric columns

ccccccccccccccccccccccccccccccc
c     0.1,5.0 (int 1)
c
c     Absorption by: 
c     CO2, O2, O, H2, N
ccccccccccccccccccccccccccccccc

c     Input atmospheric column
      indexint=1
      do i=1,klev
         auxcolinp(klev-i+1) = co2colx(i)*crscabsi2(1,indexint) +
c     $        o2colx(i)*crscabsi2(2,indexint) + 
     $        o3pcolx(i)*crscabsi2(3,indexint) 
c     $        h2colx(i)*crscabsi2(5,indexint) + 
      end do
      limdown=1.e-20
      limup=1.e26


c     Interpolations

      do i=1,nz2
         auxi = nz2-i+1
         !CO2 tabulated coefficient
         auxjco2(i) = jabsifotsintpar(auxi,1,indexint)
         !O2 tabulated coefficient
c         auxjo2(i) = jabsifotsintpar(auxi,2,indexint)
         !O3p tabulated coefficient
         auxjo3p(i) = jabsifotsintpar(auxi,3,indexint)
         !H2 tabulated coefficient
c         auxjh2(i) = jabsifotsintpar(auxi,5,indexint)
         !N tabulated coefficient
c         auxjn(i) = jabsifotsintpar(auxi,9,indexint)
         !Tabulated column
         auxcoltab(i) = c1_16(auxi,indexint)
      enddo

      !Only if chemthermod.ge.2
      !N tabulated coefficient
c      if(chemthermod.ge.2) then
c         do i=1,nz2
c            auxjn(i) = jabsifotsintpar(nz2-i+1,9,indexint)
c         enddo
c      endif

      call interfast 
     $     (wm,wp,auxind,auxcolinp,klev,auxcoltab,nz2,limdown,limup)
      do i=1,klev
         ind=auxind(i)
         auxi=klev-i+1
         !CO2 interpolated coefficient
         jfotsout(indexint,1,auxi) = wm(i)*auxjco2(ind+1) +
     $        wp(i)*auxjco2(ind)
         !O2 interpolated coefficient
c         jfotsout(indexint,2,auxi) = wm(i)*auxjo2(ind+1) +
c     $        wp(i)*auxjo2(ind)
         !O3p interpolated coefficient
          jfotsout(indexint,3,auxi) = wm(i)*auxjo3p(ind+1) +
     $        wp(i)*auxjo3p(ind)
          !H2 interpolated coefficient
c          jfotsout(indexint,5,auxi) = wm(i)*auxjh2(ind+1) +
c     $         wp(i)*auxjh2(ind)
c      !N interpolated coefficient
c            jfotsout(indexint,9,auxi) =  wm(i)*auxjn(ind+1) +
c     $         wp(i)*auxjn(ind)

C      print*, '--- L214 jthermcal.F ---'
C      print*, jfotsout(indexint,1,auxi)
c      STOP         

      enddo

      !Only if chemthermod.ge.2
      !N interpolated coefficient
c      if(chemthermod.ge.2) then
c         do i=1,klev
c            ind=auxind(i)
c            jfotsout(indexint,9,klev-i+1) =  wm(i)*auxjn(ind+1) +
c     $         wp(i)*auxjn(ind)
c         enddo
c      endif

	
c     End interval 1


ccccccccccccccccccccccccccccccc
c     5-80.5nm (int 2-15)
c
c     Absorption by:
c     CO2, O2, O, H2, N2, N, 
c     NO, CO, H, NO2
ccccccccccccccccccccccccccccccc

c     Input atmospheric column
      do indexint=2,15
         do i=1,klev
            auxcolinp(klev-i+1) = co2colx(i)*crscabsi2(1,indexint)+
     $           o3pcolx(i)*crscabsi2(3,indexint)+
     $           n2colx(i)*crscabsi2(8,indexint)+
     $           cocolx(i)*crscabsi2(11,indexint)

c     $           o2colx(i)*crscabsi2(2,indexint)+
c     $           h2colx(i)*crscabsi2(5,indexint)+
c     $           nocolx(i)*crscabsi2(10,indexint)+
c     $           hcolx(i)*crscabsi2(12,indexint)+
c     $           no2colx(i)*crscabsi2(13,indexint)
         end do

c     Interpolations

         do i=1,nz2
            auxi = nz2-i+1
            !O2 tabulated coefficient
c            auxjo2(i) = jabsifotsintpar(auxi,2,indexint)
            !O3p tabulated coefficient
            auxjo3p(i) = jabsifotsintpar(auxi,3,indexint)
            !CO2 tabulated coefficient
            auxjco2(i) = jabsifotsintpar(auxi,1,indexint)
            !H2 tabulated coefficient
c            auxjh2(i) = jabsifotsintpar(auxi,5,indexint)
            !N2 tabulated coefficient
            auxjn2(i) = jabsifotsintpar(auxi,8,indexint)
            !N tabulated coefficient
c            auxjn(i) = jabsifotsintpar(auxi,9,indexint)
            !CO tabulated coefficient
            auxjco(i) = jabsifotsintpar(auxi,11,indexint)
            !H tabulated coefficient
c            auxjh(i) = jabsifotsintpar(auxi,12,indexint)
            !tabulated column
            auxcoltab(i) = c1_16(auxi,indexint)
         enddo
c         !Only if chemthermod.ge.2
c         if(chemthermod.ge.2) then
c            do i=1,nz2
c               auxi = nz2-i+1
c               !N tabulated coefficient
c               auxjn(i) = jabsifotsintpar(auxi,9,indexint)
c               !NO tabulated coefficient
c               auxjno(i) = jabsifotsintpar(auxi,10,indexint)
c               !NO2 tabulated coefficient
c               auxjno2(i) = jabsifotsintpar(auxi,13,indexint)
c            enddo
c         endif

          call interfast(wm,wp,auxind,auxcolinp,klev,
     $        auxcoltab,nz2,limdown,limup)
          do i=1,klev
             ind=auxind(i)
             auxi = klev-i+1
             !O2 interpolated coefficient
c             jfotsout(indexint,2,auxi) = wm(i)*auxjo2(ind+1) +
c     $            wp(i)*auxjo2(ind)
             !O3p interpolated coefficient
             jfotsout(indexint,3,auxi) = wm(i)*auxjo3p(ind+1) +
     $            wp(i)*auxjo3p(ind)
             !CO2 interpolated coefficient
             jfotsout(indexint,1,auxi) = wm(i)*auxjco2(ind+1) +
     $            wp(i)*auxjco2(ind)
             !H2 interpolated coefficient
c             jfotsout(indexint,5,auxi) = wm(i)*auxjh2(ind+1) +
c     $            wp(i)*auxjh2(ind)
             !N2 interpolated coefficient
             jfotsout(indexint,8,auxi) = wm(i)*auxjn2(ind+1) +
     $            wp(i)*auxjn2(ind)             
             !N interpolated coefficient
c             jfotsout(indexint,9,auxi) = wm(i)*auxjn(ind+1) +
c     $            wp(i)*auxjn(ind)
             !CO interpolated coefficient
             jfotsout(indexint,11,auxi) = wm(i)*auxjco(ind+1) +
     $            wp(i)*auxjco(ind)
             !H interpolated coefficient
c             jfotsout(indexint,12,auxi) = wm(i)*auxjh(ind+1) +
c     $            wp(i)*auxjh(ind)             
          enddo
          !Only if chemthermod.ge.2
c          if(chemthermod.ge.2) then
c             do i=1,klev
c                ind=auxind(i)
c                auxi = klev-i+1
c                !N interpolated coefficient
c                jfotsout(indexint,9,auxi) = wm(i)*auxjn(ind+1) +
c     $               wp(i)*auxjn(ind)
c                !NO interpolated coefficient
c                jfotsout(indexint,10,auxi)=wm(i)*auxjno(ind+1) +
c     $               wp(i)*auxjno(ind)
c                !NO2 interpolated coefficient
c                jfotsout(indexint,13,auxi)=wm(i)*auxjno2(ind+1)+
c     $               wp(i)*auxjno2(ind)
c             enddo
c          endif   
      end do

c     End intervals 2-15

	
ccccccccccccccccccccccccccccccc
c     80.6-90.8nm (int16)
c
c     Absorption by:
c     CO2, O2, O, N2, N, NO,
c     CO, H, NO2
ccccccccccccccccccccccccccccccc

c     Input atmospheric column
      indexint=16
      do i=1,klev
         auxcolinp(klev-i+1) = co2colx(i)*crscabsi2(1,indexint)+
c     $        o2colx(i)*crscabsi2(2,indexint)+
     $        o3pcolx(i)*crscabsi2(3,indexint)+
     $        n2colx(i)*crscabsi2(8,indexint)+
     $        cocolx(i)*crscabsi2(11,indexint)
c     $        hcolx(i)*crscabsi2(12,indexint)+
c     $        no2colx(i)*crscabsi2(13,indexint)
      end do

c     Interpolations

      do i=1,nz2
         auxi = nz2-i+1
         !O2 tabulated coefficient
c         auxjo2(i) = jabsifotsintpar(auxi,2,indexint)
         !CO2 tabulated coefficient
         auxjco2(i) = jabsifotsintpar(auxi,1,indexint)
         !O3p tabulated coefficient
         auxjo3p(i) = jabsifotsintpar(auxi,3,indexint)
         !N2 tabulated coefficient
         auxjn2(i) = jabsifotsintpar(auxi,8,indexint)
         !CO tabulated coefficient
         auxjco(i) = jabsifotsintpar(auxi,11,indexint)
c        !N tabulated coefficient
c         auxjn(i) = jabsifotsintpar(auxi,9,indexint)

         !NO tabulated coefficient
c         auxjno(i) = jabsifotsintpar(auxi,10,indexint)
         !H tabulated coefficient
c         auxjh(i) = jabsifotsintpar(auxi,12,indexint)
         !NO2 tabulated coefficient
c         auxjno2(i) = jabsifotsintpar(auxi,13,indexint)

         !Tabulated column
         auxcoltab(i) = c1_16(auxi,indexint)
      enddo
      !Only if chemthermod.ge.2
c      if(chemthermod.ge.2) then
c         do i=1,nz2
c            auxi = nz2-i+1
c            !N tabulated coefficient
c            auxjn(i) = jabsifotsintpar(auxi,9,indexint)
c            !NO tabulated coefficient
c            auxjno(i) = jabsifotsintpar(auxi,10,indexint)
c            !NO2 tabulated coefficient
c            auxjno2(i) = jabsifotsintpar(auxi,13,indexint)
c         enddo
c      endif

      call interfast
     $     (wm,wp,auxind,auxcolinp,klev,auxcoltab,nz2,limdown,limup)
      do i=1,klev
         ind=auxind(i)
         auxi = klev-i+1
         !O2 interpolated coefficient
c         jfotsout(indexint,2,auxi) = wm(i)*auxjo2(ind+1) +
c     $            wp(i)*auxjo2(ind)
         !CO2 interpolated coefficient
         jfotsout(indexint,1,auxi) = wm(i)*auxjco2(ind+1) +
     $            wp(i)*auxjco2(ind) 
         !O3p interpolated coefficient
         jfotsout(indexint,3,auxi) = wm(i)*auxjo3p(ind+1) +
     $            wp(i)*auxjo3p(ind) 
         !N2 interpolated coefficient
         jfotsout(indexint,8,auxi) = wm(i)*auxjn2(ind+1) +
     $            wp(i)*auxjn2(ind) 
         !CO interpolated coefficient
         jfotsout(indexint,11,auxi) = wm(i)*auxjco(ind+1) +
     $            wp(i)*auxjco(ind)  
         !N interpolated coefficient
c          jfotsout(indexint,9,auxi) = wm(i)*auxjn(ind+1) +
c     $           wp(i)*auxjn(ind) 

c         !H interpolated coefficient
c         jfotsout(indexint,12,auxi) = wm(i)*auxjh(ind+1) +
c     $            wp(i)*auxjh(ind)         
      enddo
      !Only if chemthermod.ge.2
c      if(chemthermod.ge.2) then
c         do i=1,klev
c            ind=auxind(i)
c            auxi = klev-i+1
c            !N interpolated coefficient
c            jfotsout(indexint,9,auxi) = wm(i)*auxjn(ind+1) +
c     $           wp(i)*auxjn(ind) 
c            !NO interpolated coefficient
c            jfotsout(indexint,10,auxi) = wm(i)*auxjno(ind+1) +
c     $           wp(i)*auxjno(ind)
c            !NO2 interpolated coefficient
c            jfotsout(indexint,13,auxi) = wm(i)*auxjno2(ind+1) +
c     $           wp(i)*auxjno2(ind)
c         enddo
c      endif
c     End interval 16

	
ccccccccccccccccccccccccccccccc
c     90.9-119.5nm (int 17-24)
c
c     Absorption by:
c     CO2, O2, N2, NO, CO, NO2
ccccccccccccccccccccccccccccccc

c     Input column

      do i=1,klev
         auxcolinp(klev-i+1) = co2colx(i) +  n2colx(i) +
     $         + cocolx(i) 
      end do

      do indexint=17,24

c     Interpolations

         do i=1,nz2
            auxi = nz2-i+1
            !CO2 tabulated coefficient
            auxjco2(i) = jabsifotsintpar(auxi,1,indexint)
            !O2 tabulated coefficient
c            auxjo2(i) = jabsifotsintpar(auxi,2,indexint)
            !N2 tabulated coefficient
            auxjn2(i) = jabsifotsintpar(auxi,8,indexint)
            !CO tabulated coefficient
            auxjco(i) = jabsifotsintpar(auxi,11,indexint)            
            !Tabulated column
            auxcoltab(i) = c17_24(auxi)
         enddo
         !Only if chemthermod.ge.2
c         if(chemthermod.ge.2) then
c            do i=1,nz2
c               auxi = nz2-i+1
c               !NO tabulated coefficient
c               auxjno(i) = jabsifotsintpar(auxi,10,indexint)
c               !NO2 tabulated coefficient
c               auxjno2(i) = jabsifotsintpar(auxi,13,indexint)
c            enddo
c         endif

         call interfast
     $     (wm,wp,auxind,auxcolinp,klev,auxcoltab,nz2,limdown,limup)
         !Correction to include T variation of CO2 cross section
         if(indexint.eq.24) then
            do i=1,klev
               auxi = klev-i+1
               if(sigma(indexint,auxi)*
     $              alfa(indexint,auxi)*coltemp(auxi)
     $              .lt.60.) then
                  cortemp(i)=exp(-sigma(indexint,auxi)*
     $                alfa(indexint,auxi)*coltemp(auxi))
               else 
                  cortemp(i)=0.
               end if
            enddo
         else
            do i=1,klev
               cortemp(i)=1.
            enddo
         end if
         do i=1,klev           
            ind=auxind(i)
            auxi = klev-i+1
            !O2 interpolated coefficient
c            jfotsout(indexint,2,auxi) = (wm(i)*auxjo2(ind+1) +
c     $           wp(i)*auxjo2(ind)) * cortemp(i)
            !CO2 interpolated coefficient
            jfotsout(indexint,1,auxi) = (wm(i)*auxjco2(ind+1) +
     $           wp(i)*auxjco2(ind)) * cortemp(i)
            if(indexint.eq.24) jfotsout(indexint,1,auxi)=
     $           jfotsout(indexint,1,auxi)*
     $           (1+alfa(indexint,auxi)*
     $           (t2(auxi)-t0(auxi)))
            !N2 interpolated coefficient
            jfotsout(indexint,8,auxi) = (wm(i)*auxjn2(ind+1) +
     $            wp(i)*auxjn2(ind)) * cortemp(i)            
            !CO interpolated coefficient
            jfotsout(indexint,11,auxi) = (wm(i)*auxjco(ind+1) +
     $            wp(i)*auxjco(ind)) * cortemp(i)            
         enddo
         !Only if chemthermod.ge.2
c         if(chemthermod.ge.2) then
c            do i=1,klev
c               ind=auxind(i)
c               auxi = klev-i+1
c               !NO interpolated coefficient
c               jfotsout(indexint,10,auxi)=(wm(i)*auxjno(ind+1) +
c     $              wp(i)*auxjno(ind)) * cortemp(i)
c               !NO2 interpolated coefficient
c               jfotsout(indexint,13,auxi)=(wm(i)*auxjno2(ind+1)+
c     $              wp(i)*auxjno2(ind)) * cortemp(i)
c            enddo
c         endif               
      end do
c     End intervals 17-24

	
ccccccccccccccccccccccccccccccc
c     119.6-167.0nm (int 25-29)
c
c     Absorption by:
c     CO2, O2, H2O, H2O2, NO,
c     CO, NO2
ccccccccccccccccccccccccccccccc

c     Input atmospheric column

      do i=1,klev
c         auxcolinp(klev-i+1) = co2colx(i) + o2colx(i) + h2ocolx(i) + 
c     $        h2o2colx(i) + nocolx(i) + cocolx(i) + no2colx(i)
         auxcolinp(klev-i+1) = co2colx(i)  + cocolx(i) 

      end do

      do indexint=25,29

c     Interpolations

         do i=1,nz2
            auxi = nz2-i+1
            !CO2 tabulated coefficient
            auxjco2(i) = jabsifotsintpar(auxi,1,indexint)
            !O2 tabulated coefficient
c            auxjo2(i) = jabsifotsintpar(auxi,2,indexint)
            !H2O tabulated coefficient
c            auxjh2o(i) = jabsifotsintpar(auxi,4,indexint)
            !H2O2 tabulated coefficient
c            auxjh2o2(i) = jabsifotsintpar(auxi,6,indexint)            
            !CO tabulated coefficient
            auxjco(i) = jabsifotsintpar(auxi,11,indexint)            
            !Tabulated column
            auxcoltab(i) = c25_29(auxi)
         enddo
         !Only if chemthermod.ge.2
c         if(chemthermod.ge.2) then
c            do i=1,nz2
c               auxi = nz2-i+1
c               !NO tabulated coefficient
c               auxjno(i) = jabsifotsintpar(auxi,10,indexint)
c               !NO2 tabulated coefficient
c               auxjno2(i) = jabsifotsintpar(auxi,13,indexint)
c            enddo
c         endif
         call interfast
     $     (wm,wp,auxind,auxcolinp,klev,auxcoltab,nz2,limdown,limup)
         do i=1,klev
            ind=auxind(i)
            auxi = klev-i+1
            !Correction to include T variation of CO2 cross section
            if(sigma(indexint,auxi)*alfa(indexint,auxi)*
     $           coltemp(auxi).lt.60.) then
               cortemp(i)=exp(-sigma(indexint,auxi)*
     $              alfa(indexint,auxi)*coltemp(auxi))
            else 
               cortemp(i)=0.
            end if
            !CO2 interpolated coefficient
            jfotsout(indexint,1,auxi) = (wm(i)*auxjco2(ind+1) +
     $           wp(i)*auxjco2(ind)) * cortemp(i) *
     $           (1+alfa(indexint,auxi)*
     $           (t2(auxi)-t0(auxi)))
            !O2 interpolated coefficient
c            jfotsout(indexint,2,auxi) = (wm(i)*auxjo2(ind+1) +
c     $           wp(i)*auxjo2(ind)) * cortemp(i)
            !H2O interpolated coefficient
c            jfotsout(indexint,4,auxi) = (wm(i)*auxjh2o(ind+1) +
c     $           wp(i)*auxjh2o(ind)) * cortemp(i)
            !H2O2 interpolated coefficient
c            jfotsout(indexint,6,auxi) = (wm(i)*auxjh2o2(ind+1) +
c     $           wp(i)*auxjh2o2(ind)) * cortemp(i)            
            !CO interpolated coefficient
            jfotsout(indexint,11,auxi) = (wm(i)*auxjco(ind+1) +
     $           wp(i)*auxjco(ind)) * cortemp(i)
         enddo
         !Only if chemthermod.ge.2
c         if(chemthermod.ge.2) then
c            do i=1,klev
c               ind=auxind(i)
c               auxi = klev-i+1
c               !NO interpolated coefficient
c               jfotsout(indexint,10,auxi)=(wm(i)*auxjno(ind+1) +
c     $              wp(i)*auxjno(ind)) * cortemp(i)
c               !NO2 interpolated coefficient
c               jfotsout(indexint,13,auxi)=(wm(i)*auxjno2(ind+1)+
c     $              wp(i)*auxjno2(ind)) * cortemp(i)
c            enddo
c         endif

      end do

c     End intervals 25-29


cccccccccccccccccccccccccccccccc
c     167.1-202.5nm (int 30-31)
c    
c     Absorption by:
c     CO2, O2, H2O, H2O2, NO,
c     NO2
cccccccccccccccccccccccccccccccc

c     Input atmospheric column

      do i=1,klev
         auxcolinp(klev-i+1) = co2colx(i)
      end do

c     Interpolation

      do indexint=30,31

         do i=1,nz2
            auxi = nz2-i+1
            !CO2 tabulated coefficient
            auxjco2(i) = jabsifotsintpar(auxi,1,indexint)
            !O2 tabulated coefficient
c            auxjo2(i) = jabsifotsintpar(auxi,2,indexint)
            !H2O tabulated coefficient
c            auxjh2o(i) = jabsifotsintpar(auxi,4,indexint)
            !H2O2 tabulated coefficient
c            auxjh2o2(i) = jabsifotsintpar(auxi,6,indexint)            
            !Tabulated column
c            auxcoltab(i) = c30_31(auxi)
         enddo
         !Only if chemthermod.ge.2
c         if(chemthermod.ge.2) then
c            do i=1,nz2
c               auxi = nz2-i+1
c               !NO tabulated coefficient
c               auxjno(i) = jabsifotsintpar(auxi,10,indexint)
c               !NO2 tabulated coefficient
c               auxjno2(i) = jabsifotsintpar(auxi,13,indexint)
c            enddo
c         endif

         call interfast
     $     (wm,wp,auxind,auxcolinp,klev,auxcoltab,nz2,limdown,limup)
         do i=1,klev
            ind=auxind(i)
            auxi = klev-i+1
            !Correction to include T variation of CO2 cross section
            if(sigma(indexint,auxi)*alfa(indexint,auxi)*
     $           coltemp(auxi).lt.60.) then
               cortemp(i)=exp(-sigma(indexint,auxi)*
     $              alfa(indexint,auxi)*coltemp(auxi))
            else 
               cortemp(i)=0.
            end if
            !CO2 interpolated coefficient
            jfotsout(indexint,1,auxi) = (wm(i)*auxjco2(ind+1) +
     $           wp(i)*auxjco2(ind)) * cortemp(i) *
     $           (1+alfa(indexint,auxi)*
     $           (t2(auxi)-t0(auxi)))
            !O2 interpolated coefficient
c            jfotsout(indexint,2,auxi) = (wm(i)*auxjo2(ind+1) +
c     $            wp(i)*auxjo2(ind)) * cortemp(i)
            !H2O interpolated coefficient
c            jfotsout(indexint,4,auxi) = (wm(i)*auxjh2o(ind+1) +
c     $            wp(i)*auxjh2o(ind)) * cortemp(i)
            !H2O2 interpolated coefficient
c            jfotsout(indexint,6,auxi) = (wm(i)*auxjh2o2(ind+1) +
c     $            wp(i)*auxjh2o2(ind)) * cortemp(i)            
         enddo
c         !Only if chemthermod.ge.2
c         if(chemthermod.ge.2) then
c            do i=1,klev 
c               ind=auxind(i)
c               auxi = klev-i+1
c               !NO interpolated coefficient
c               jfotsout(indexint,10,auxi)=(wm(i)*auxjno(ind+1) +
c     $              wp(i)*auxjno(ind)) * cortemp(i)
c               !NO2 interpolated coefficient
c               jfotsout(indexint,13,auxi)=(wm(i)*auxjno2(ind+1)+
c     $              wp(i)*auxjno2(ind)) * cortemp(i)
c            enddo
c         endif

      end do

c     End intervals 30-31


ccccccccccccccccccccccccccccccc
c     202.6-210.0nm (int 32)
c
c     Absorption by:
c     CO2, O2, H2O2, NO, NO2
ccccccccccccccccccccccccccccccc

c     Input atmospheric column

      indexint=32
      do i=1,klev
         auxcolinp(klev-i+1) =co2colx(i)
      end do

c     Interpolation

      do i=1,nz2
         auxi = nz2-i+1
         !CO2 tabulated coefficient
         auxjco2(i) = jabsifotsintpar(auxi,1,indexint)
         !O2 tabulated coefficient
c         auxjo2(i) = jabsifotsintpar(auxi,2,indexint)
         !H2O2 tabulated coefficient
c         auxjh2o2(i) = jabsifotsintpar(auxi,6,indexint)         
         !Tabulated column
         auxcoltab(i) = c32(auxi)
      enddo
      !Only if chemthermod.ge.2
c      if(chemthermod.ge.2) then
c         do i=1,nz2
c            auxi = nz2-i+1
c            !NO tabulated coefficient
c            auxjno(i) = jabsifotsintpar(auxi,10,indexint)
c            !NO2 tabulated coefficient
c            auxjno2(i) = jabsifotsintpar(auxi,13,indexint)
c         enddo
c      endif
      call interfast
     $     (wm,wp,auxind,auxcolinp,klev,auxcoltab,nz2,limdown,limup)
      do i=1,klev
         ind=auxind(i)
         auxi = klev-i+1
         !Correction to include T variation of CO2 cross section
         if(sigma(indexint,klev-i+1)*alfa(indexint,auxi)*
     $        coltemp(auxi).lt.60.) then
            cortemp(i)=exp(-sigma(indexint,auxi)*
     $           alfa(indexint,auxi)*coltemp(auxi))
         else 
            cortemp(i)=0.
         end if
         !CO2 interpolated coefficient
         jfotsout(indexint,1,auxi) = (wm(i)*auxjco2(ind+1) +
     $        wp(i)*auxjco2(ind)) * cortemp(i) *
     $        (1+alfa(indexint,auxi)*
     $        (t2(auxi)-t0(auxi)))
         !O2 interpolated coefficient
c         jfotsout(indexint,2,auxi) = (wm(i)*auxjo2(ind+1) +
c     $        wp(i)*auxjo2(ind)) * cortemp(i)
         !H2O2 interpolated coefficient
c         jfotsout(indexint,6,auxi) = (wm(i)*auxjh2o2(ind+1) +
c     $        wp(i)*auxjh2o2(ind)) * cortemp(i)         
      enddo
      !Only if chemthermod.ge.2
c      if(chemthermod.ge.2) then
c         do i=1,klev
c            auxi = klev-i+1
c            ind=auxind(i)
c            !NO interpolated coefficient
c            jfotsout(indexint,10,auxi) = (wm(i)*auxjno(ind+1) +
c     $           wp(i)*auxjno(ind)) * cortemp(i)
           !NO2 interpolated coefficient
c            jfotsout(indexint,13,auxi) = (wm(i)*auxjno2(ind+1) +
c     $           wp(i)*auxjno2(ind)) * cortemp(i)
c         enddo
c      endif

c     End of interval 32


ccccccccccccccccccccccccccccccc
c     210.1-231.0nm (int 33)
c     
c     Absorption by:
c     O2, H2O2, NO2
ccccccccccccccccccccccccccccccc

c     Input atmospheric column

c      indexint=33
c      do i=1,klev
c         auxcolinp(klev-i+1) = o2colx(i) + h2o2colx(i) + no2colx(i)
c      end do

c     Interpolation

c      do i=1,nz2
c         auxi = nz2-i+1
         !O2 tabulated coefficient
c         auxjo2(i) = jabsifotsintpar(auxi,2,indexint)
         !H2O2 tabulated coefficient
c         auxjh2o2(i) = jabsifotsintpar(auxi,6,indexint)
         !Tabulated column
c        auxcoltab(i) = c33(auxi)
c     enddo
      !Only if chemthermod.ge.2
c     if(chemthermod.ge.2) then
c        do i=1,nz2
            !NO2 tabulated coefficient
c           auxjno2(i) = jabsifotsintpar(nz2-i+1,13,indexint)
c        enddo
c     endif
c     call interfast
c    $     (wm,wp,auxind,auxcolinp,klev,auxcoltab,nz2,limdown,limup)
c     do i=1,klev
c        ind=auxind(i)
c        auxi = klev-i+1
         !O2 interpolated coefficient
c        jfotsout(indexint,2,auxi) = wm(i)*auxjo2(ind+1) +
c    $        wp(i)*auxjo2(ind)
c        !H2O2 interpolated coefficient
c        jfotsout(indexint,6,auxi) = wm(i)*auxjh2o2(ind+1) +
c    $        wp(i)*auxjh2o2(ind)         
c     enddo
      !Only if chemthermod.ge.2
c     if(chemthermod.ge.2) then
c        do i=1,klev
c           ind=auxind(i)
c           !NO2 interpolated coefficient
c           jfotsout(indexint,13,klev-i+1) = wm(i)*auxjno2(ind+1) +
c    $           wp(i)*auxjno2(ind)
c        enddo
c     endif

c     End of interval 33


ccccccccccccccccccccccccccccccc
c     231.1-240.0nm (int 34)
c
c     Absorption by:
c     O2, H2O2, O3, NO2
ccccccccccccccccccccccccccccccc

c     Input atmospheric column

c     indexint=34
c     do i=1,klev
c        auxcolinp(klev-i+1) = h2o2colx(i) + o2colx(i) + o3colx(i) + 
c    $        no2colx(i)
c     end do

c     Interpolation

c     do i=1,nz2
c        auxi = nz2-i+1
         !O2 tabulated coefficient
c        auxjo2(i) = jabsifotsintpar(auxi,2,indexint)
         !H2O2 tabulated coefficient
c        auxjh2o2(i) = jabsifotsintpar(auxi,6,indexint)
         !O3 tabulated coefficient
c        auxjo3(i) = jabsifotsintpar(auxi,7,indexint)         
         !Tabulated column
c        auxcoltab(i) = c34(nz2-i+1)
c     enddo
      !Only if chemthermod.ge.2
c     if(chemthermod.ge.2) then
c        do i=1,nz2
            !NO2 tabulated coefficient
c           auxjno2(i) = jabsifotsintpar(nz2-i+1,13,indexint)
c        enddo
c     endif
c     call interfast
c    $     (wm,wp,auxind,auxcolinp,klev,auxcoltab,nz2,limdown,limup)
c     do i=1,klev
c        ind=auxind(i)
c        auxi = klev-i+1
         !O2 interpolated coefficient
c        jfotsout(indexint,2,auxi) = wm(i)*auxjo2(ind+1) +
c    $        wp(i)*auxjo2(ind)
         !H2O2 interpolated coefficient
c        jfotsout(indexint,6,auxi) = wm(i)*auxjh2o2(ind+1) +
c    $        wp(i)*auxjh2o2(ind)
         !O3 interpolated coefficient
c        jfotsout(indexint,7,auxi) = wm(i)*auxjo3(ind+1) +
c    $        wp(i)*auxjo3(ind)         
c     enddo
      !Only if chemthermod.ge.2
c     if(chemthermod.ge.2) then
c        do i=1,klev
c           ind=auxind(i)
            !NO2 interpolated coefficient
c           jfotsout(indexint,13,klev-i+1) = wm(i)*auxjno2(ind+1) +
c    $           wp(i)*auxjno2(ind)
c        enddo
c     endif

c     End of interval 34      


ccccccccccccccccccccccccccccccc
c     240.1-337.7nm (int 35)
c
c     Absorption by:
c     H2O2, O3, NO2
ccccccccccccccccccccccccccccccc

c     Input atmospheric column

      indexint=35
c      do i=1,klev
c         auxcolinp(klev-i+1) = o3colx(i) 
c      end do
c
c     Interpolation

c      do i=1,nz2
c         auxi = nz2-i+1
         !H2O2 tabulated coefficient
c         auxjh2o2(i) = jabsifotsintpar(auxi,6,indexint)
         !O3 tabulated coefficient
c         auxjo3(i) = jabsifotsintpar(auxi,7,indexint)
         !Tabulated column
c         auxcoltab(i) = c35(auxi)
c      enddo
      !Only if chemthermod.ge.2
c      if(chemthermod.ge.2) then
c         do i=1,nz2
c            !NO2 tabulated coefficient
c            auxjno2(i) = jabsifotsintpar(nz2-i+1,13,indexint)
c         enddo
c      endif
c      call interfast
c     $     (wm,wp,auxind,auxcolinp,klev,auxcoltab,nz2,limdown,limup)
c      do i=1,klev
c         ind=auxind(i)
c         auxi = klev-i+1
c         !H2O2 interpolated coefficient
c         jfotsout(indexint,6,auxi) = wm(i)*auxjh2o2(ind+1) +
c     $        wp(i)*auxjh2o2(ind)
c         !O3 interpolated coefficient
c         jfotsout(indexint,7,auxi) = wm(i)*auxjo3(ind+1) +
c     $        wp(i)*auxjo3(ind)         
c      enddo
c      if(chemthermod.ge.2) then
c         do i=1,klev
c            ind=auxind(i)
c            !NO2 interpolated coefficient
c            jfotsout(indexint,13,klev-i+1) = wm(i)*auxjno2(ind+1) +
c     $           wp(i)*auxjno2(ind)
c         enddo
c      endif

c     End of interval 35

ccccccccccccccccccccccccccccccc
c     337.8-800.0 nm (int 36)
c     
c     Absorption by:
c     O3, NO2
ccccccccccccccccccccccccccccccc

c     Input atmospheric column

      indexint=36
c      do i=1,klev
c         auxcolinp(klev-i+1) = o3colx(i) 
c      end do

c     Interpolation

c      do i=1,nz2
c         auxi = nz2-i+1
         !O3 tabulated coefficient
c         auxjo3(i) = jabsifotsintpar(auxi,7,indexint)         
         !Tabulated column
c         auxcoltab(i) = c36(auxi)
c      enddo
      !Only if chemthermod.ge.2
c      if(chemthermod.ge.2) then
c         do i=1,nz2
c            !NO2 tabulated coefficient
c            auxjno2(i) = jabsifotsintpar(nz2-i+1,13,indexint)
c         enddo
c      endif
c      call interfast
c     $     (wm,wp,auxind,auxcolinp,klev,auxcoltab,nz2,limdown,limup)
c      do i=1,klev
c         ind=auxind(i)
c         !O3 interpolated coefficient
c         jfotsout(indexint,7,klev-i+1) = wm(i)*auxjo3(ind+1) +
c     $        wp(i)*auxjo3(ind)         
c      enddo
c      !Only if chemthermod.ge.2
c      if(chemthermod.ge.2) then
c         do i=1,klev
c            ind=auxind(i)
c           !NO2 interpolated coefficient
c            jfotsout(indexint,13,klev-i+1) = wm(i)*auxjno2(ind+1) +
c     $           wp(i)*auxjno2(ind)
c         enddo
c      endif

c     End of interval 36

c     End of interpolation to obtain photoabsorption rates

	
      return
      end



c**********************************************************************
c**********************************************************************

      subroutine column(ig,chemthermod,rm,nesptherm,tx,iz,zenit,
     $     co2colx,o3pcolx, n2colx, cocolx)

c     mar 2014        gg            adapted to Venus GCM
c     nov 2002        fgg           first version

c**********************************************************************
      use dimphy
      use conc
      implicit none


c     common variables and constants
c#include "tracer.h"
#include "param.h"
#include "param_v4.h"
#include "clesphys.h"
#include "mmol.h"


c    local parameters and variables

c     input and output variables

      integer    ig
      integer    chemthermod
      integer    nesptherm                   !# of species undergoing chemistry, input
      real       rm(klev,nesptherm)         !densities (cm-3), input
      real       tx(klev)                   !temperature profile, input
      real       iz(klev+1)                 !height profile, input
      real       zenit                      !SZA, input
      real       co2colx(klev)              !column density of CO2 (cm^-2), output
      real       o3pcolx(klev)              !column density of O(3P)(cm^-2), output
      real       n2colx(klev)               !N2 column density (cm-2), output
      real       cocolx(klev)               !CO column density (cm-2), output

c      real       o2colx(klev)               !column density of O2(cm^-2), output
c      real       h2colx(klev)               !H2 column density (cm-2), output
c      real       h2ocolx(klev)              !H2O column density (cm-2), output
c      real       h2o2colx(klev)             !column density of H2O2(cm^-2), output
c      real       o3colx(klev)               !O3 column density (cm-2), output
c      real       nocolx(klev)               !NO column density (cm-2), output
c      real       hcolx(klev)                !H column density (cm-2), output
c      real       no2colx(klev)              !NO2 column density (cm-2), output
      

c     local variables

      real       xx
      real       grav(klev)
      real       Hco2,Ho3p,Ho2,Hh2,Hh2o,Hh2o2
      real       Ho3,Hn2,Hn,Hno,Hco,Hh,Hno2

      real       co2x(klev)
      real       o3px(klev)
      real       cox(klev)
      real       n2x(klev)
      real       nx(klev)

c      real       o2x(klev)
c      real       o3x(klev)
c      real       hx(klev)
c      real       h2x(klev)
c      real       h2ox(klev)
c      real       h2o2x(klev)
c      real       nox(klev)
c      real       no2x(klev)

      integer    i,j,k,icol,indexint          !indexes

c     variables for optical path calculation

      integer    nz3
!      parameter  (nz3=nz*2)

      integer    jj
      real*8      esp(klev*2)
      real*8      ilayesp(klev*2)
      real*8      szalayesp(klev*2)
      integer     nlayesp
      real*8      zmini
      real*8      depth
      real*8      espco2, espo2, espo3p, esph2, esph2o, esph2o2,espo3
      real*8      espn2,espn,espno,espco,esph,espno2
      real*8      rcmnz, rcmmini
      real*8      szadeg

! Tracer indexes in the thermospheric chemistry:
!!! ATTENTION. These values have to be identical to those in euvheat.F90
!!! If the values are changed there, the same has to be done here  !!!

      integer,parameter :: ix_co2=1
      integer,parameter :: ix_n2=13
      integer,parameter :: ix_o=3
      integer,parameter :: ix_co=4

c*************************PROGRAM STARTS*******************************

      nz3 = klev*2
      do i=1,klev
        xx = ( radio + iz(i) ) * 1.e5    ! conversion [km] ---> [cm]
         grav(i) = gg * masa /(xx**2)    ! [cm/s2]
      end do

      !Scale heights  H = kT /Mg  -->  [cm]
      xx = kboltzman * tx(klev) * n_avog / grav(klev)   ! g cm mol-1

      Ho3p  = xx / mmolo
      Hco2  = xx / mmolco2
      Hco   = xx / mmolco
      Hn2   = xx / mmoln2
      Hn    = xx / mmoln

      !Only if O3 chem. required
c      if(chemthermod.ge.1) 
!     $     Ho3   = xx / mmol(igcm_o3)
c     $     Ho3   = xx / mmolo3
c      !Only if N or ion chem.
c      if(chemthermod.ge.2) then
c         Hn2   = xx / mmoln2
c         Hn    = xx / mmoln
c         Hno   = xx / mmolno
c         Hno2  = xx / mmolno2
c      endif
      ! first loop in altitude : initialisation
      do i=klev,1,-1
         !Column initialisation
         co2colx(i)  = 0.
         o3pcolx(i)  = 0.
         n2colx(i)   = 0.
         cocolx(i)   = 0.

         !--Densities [cm-3]
         co2x(i)  = rm(i,ix_co2)
         o3px(i)  = rm(i,ix_o)
         cox(i)   = rm(i,ix_co)
         n2x(i)   = rm(i,ix_n2)

c         write(*,*), '--jthermcalc--', co2x(i)

         !Only if O3 chem. required
c         if(chemthermod.ge.1) 
c     $        o3x(i)   = rm(i,i_o3)
c         !Only if Nitrogen of ion chemistry requested
c         if(chemthermod.ge.2) then
c            n2x(i)   = rm(i,i_n2)
c            nx(i)    = rm(i,i_n)
c            nox(i)   = rm(i,i_no)
c            no2x(i)  = rm(i,i_no2)
c         endif
      enddo     ! end first loop
      ! second loop in altitude : column calculations
      do i=klev,1,-1
         !Routine to calculate the geometrical length of each layer
         call espesor_optico_A(ig,i,zenit,iz(i),nz3,iz,esp,ilayesp,
     $         szalayesp,nlayesp, zmini)
         if(ilayesp(nlayesp).eq.-1) then
            co2colx(i)=1.e25
            o3pcolx(i)=1.e25
            n2colx(i)=1.e25
            cocolx(i)=1.e25

c            o2colx(i)=1.e25
c            o3pcolx(i)=1.e25
c            h2colx(i)=1.e25
c            h2ocolx(i)=1.e25
c            h2o2colx(i)=1.e25
c            o3colx(i)=1.e25
c            ncolx(i)=1.e25
c            nocolx(i)=1.e25
c            cocolx(i)=1.e25
c            hcolx(i)=1.e25
c            no2colx(i)=1.e25
         else
            rcmnz = ( radio + iz(klev) ) * 1.e5    ! km --> cm
            rcmmini = ( radio + zmini ) * 1.e5
            !Column calculation taking into account the geometrical depth 
            !calculated before
            do j=1,nlayesp
               jj=ilayesp(j)
               !Top layer
               if(jj.eq.klev) then
                  if(zenit.le.60.) then 
                     o3pcolx(i)=o3pcolx(i)+o3px(klev)*Ho3p*esp(j)
     $                    *1.e-5
                     co2colx(i)=co2colx(i)+co2x(klev)*Hco2*esp(j)
     $                    *1.e-5
                     cocolx(i)=cocolx(i)+cox(klev)*Hco*esp(j)
     $                    *1.e-5
                     n2colx(i)=n2colx(i)+n2x(klev)*Hn2*esp(j)
     $                    *1.e-5

c                     h2o2colx(i)=h2o2colx(i)+
c     $                    h2o2x(klev)*Hh2o2*esp(j)*1.e-5
c                     o2colx(i)=o2colx(i)+o2x(klev)*Ho2*esp(j)
c     $                    *1.e-5
c                     h2colx(i)=h2colx(i)+h2x(klev)*Hh2*esp(j)
c     $                    *1.e-5
c                     h2ocolx(i)=h2ocolx(i)+h2ox(klev)*Hh2o*esp(j)
c     $                    *1.e-5                     
c                     cocolx(i)=cocolx(i)+cox(klev)*Hco*esp(j)
c     $                    *1.e-5
c                     hcolx(i)=hcolx(i)+hx(klev)*Hh*esp(j)
c     $                    *1.e-5
                     !Only if O3 chemistry required
c                     if(chemthermod.ge.1) o3colx(i)=
c     $                    o3colx(i)+o3x(klev)*Ho3*esp(j)
c     $                    *1.e-5
                     !Only if N or ion chemistry requested
c                     if(chemthermod.ge.2) then
c                        n2colx(i)=n2colx(i)+n2x(klev)*Hn2*esp(j)
c     $                    *1.e-5

c                     endif
                  else if(zenit.gt.60.) then
                     espco2 =sqrt((rcmnz+Hco2)**2 -rcmmini**2) - esp(j)
                     espo3p = sqrt((rcmnz+Ho3p)**2 -rcmmini**2)- esp(j)
                     espco  = sqrt((rcmnz+Hco)**2 -rcmmini**2) - esp(j)
                     espn2 =sqrt((rcmnz+Hn2)**2-rcmmini**2)-esp(j)
                     espn  =sqrt((rcmnz+Hn)**2-rcmmini**2)  - esp(j)

c                     espo2  = sqrt((rcmnz+Ho2)**2 -rcmmini**2) - esp(j)
c                     esph2  = sqrt((rcmnz+Hh2)**2 -rcmmini**2) - esp(j)
c                     esph2o = sqrt((rcmnz+Hh2o)**2 -rcmmini**2)- esp(j)
c                     esph2o2= sqrt((rcmnz+Hh2o2)**2-rcmmini**2)- esp(j)
c                     esph   = sqrt((rcmnz+Hh)**2 -rcmmini**2)  - esp(j)
                     !Only if O3 chemistry required
c                     if(chemthermod.ge.1)                     
c     $                   espo3=sqrt((rcmnz+Ho3)**2-rcmmini**2)-esp(j)
c                     !Only if N or ion chemistry requested
c                     if(chemthermod.ge.2) then
c                        espn2 =sqrt((rcmnz+Hn2)**2-rcmmini**2)-esp(j)
c                        espn  =sqrt((rcmnz+Hn)**2-rcmmini**2)  - esp(j)
c                        espno =sqrt((rcmnz+Hno)**2-rcmmini**2) - esp(j)
c                        espno2=sqrt((rcmnz+Hno2)**2-rcmmini**2)- esp(j)
c                     endif

                     co2colx(i) = co2colx(i) + espco2*co2x(klev)
                     o3pcolx(i) = o3pcolx(i) + espo3p*o3px(klev)
                     cocolx(i)  = cocolx(i)  + espco*cox(klev)
                     n2colx(i)  = n2colx(i)  + espn2*n2x(klev)

c                     o2colx(i)  = o2colx(i)  + espo2*o2x(klev)
c                     h2colx(i)  = h2colx(i)  + esph2*h2x(klev)
c                     h2ocolx(i) = h2ocolx(i) + esph2o*h2ox(klev)
c                     h2o2colx(i)= h2o2colx(i)+ esph2o2*h2o2x(klev)
c                     cocolx(i)  = cocolx(i)  + espco*cox(klev)
c                     hcolx(i)   = hcolx(i)   + esph*hx(klev)
                     !Only if O3 chemistry required
c                     if(chemthermod.ge.1)                      
c     $                  o3colx(i) = o3colx(i)  + espo3*o3x(klev)
c                     !Only if N or ion chemistry requested
c                     if(chemthermod.ge.2) then
c                        n2colx(i)  = n2colx(i)  + espn2*n2x(klev)
c                        ncolx(i)   = ncolx(i)   + espn*nx(klev)
c                        nocolx(i)  = nocolx(i)  + espno*nox(klev)
c                        no2colx(i) = no2colx(i) + espno2*no2x(klev)
c                     endif
                  endif !Of if zenit.lt.60
               !Other layers
               else 
                  co2colx(i)  = co2colx(i) + 
     $                 esp(j) * (co2x(jj)+co2x(jj+1)) / 2.
                  o3pcolx(i)  = o3pcolx(i) + 
     $                 esp(j) * (o3px(jj)+o3px(jj+1)) / 2.
                  cocolx(i)   = cocolx(i) + 
     $                 esp(j) * (cox(jj)+cox(jj+1)) / 2.
                  n2colx(i)   = n2colx(i) + 
     $                 esp(j) * (n2x(jj)+n2x(jj+1)) / 2.

c
c                  o2colx(i)   = o2colx(i) + 
c     $                 esp(j) * (o2x(jj)+o2x(jj+1)) / 2.
c                  h2colx(i)   = h2colx(i) + 
c     $                 esp(j) * (h2x(jj)+h2x(jj+1)) / 2.
c                  h2ocolx(i)  = h2ocolx(i) + 
c     $                 esp(j) * (h2ox(jj)+h2ox(jj+1)) / 2.
c                  h2o2colx(i) = h2o2colx(i) + 
c     $                 esp(j) * (h2o2x(jj)+h2o2x(jj+1)) / 2.
c                  hcolx(i)    = hcolx(i) + 
c     $                 esp(j) * (hx(jj)+hx(jj+1)) / 2.
                  !Only if O3 chemistry required
c                  if(chemthermod.ge.1) 
c     $                 o3colx(i) = o3colx(i) + 
c     $                 esp(j) * (o3x(jj)+o3x(jj+1)) / 2.
c                  !Only if N or ion chemistry requested
c                  if(chemthermod.ge.2) then
c                     n2colx(i)   = n2colx(i) + 
c     $                 esp(j) * (n2x(jj)+n2x(jj+1)) / 2.
c                     ncolx(i)    = ncolx(i) + 
c     $                    esp(j) * (nx(jj)+nx(jj+1)) / 2.
c                     nocolx(i)   = nocolx(i) + 
c     $                    esp(j) * (nox(jj)+nox(jj+1)) / 2.
c                     no2colx(i)  = no2colx(i) + 
c     $                    esp(j) * (no2x(jj)+no2x(jj+1)) / 2.
c                  endif

               endif  !Of if jj.eq.klev
            end do    !Of do j=1,nlayesp
         endif        !Of ilayesp(nlayesp).eq.-1
      enddo           !Of do i=klev,1,-1
      return


      end


c**********************************************************************
c**********************************************************************

      subroutine interfast(wm,wp,nm,p,nlayer,pin,nl,limdown,limup)
C
C subroutine to perform linear interpolation in pressure from 1D profile 
C escin(nl) sampled on pressure grid pin(nl) to profile
C escout(nlayer) on pressure grid p(nlayer).
C
      real*8 wm(nlayer),wp(nlayer),p(nlayer)
      integer nm(nlayer)
      real*8 pin(nl)
      real*8 limup,limdown
      integer nl,nlayer,n1,n,np,nini
      nini=1
      do n1=1,nlayer
         if(p(n1) .gt. limup .or. p(n1) .lt. limdown) then
            wm(n1) = 0.d0
            wp(n1) = 0.d0
         else
            do n = nini,nl-1
               if (p(n1).ge.pin(n).and.p(n1).le.pin(n+1)) then
                  nm(n1)=n
                  np=n+1
                  wm(n1)=abs(pin(n)-p(n1))/(pin(np)-pin(n))
                  wp(n1)=1.d0 - wm(n1)
                  nini = n
                  exit
               endif
            enddo
         endif
      enddo
      return
      end


c**********************************************************************
c**********************************************************************

      subroutine espesor_optico_A (ig,capa, szadeg,z,
     @                   nz3,iz,esp,ilayesp,szalayesp,nlayesp, zmini)

c     fgg              nov 03      Adaptation to Martian model
c     malv             jul 03      Corrected z grid. Split in alt & frec codes
c     fgg              feb 03      first version
*************************************************************************
      use dimphy
      implicit none


c     common variables and constants
#include "param.h"
#include "param_v4.h"

c     arguments

      real        szadeg                ! I. SZA [rad]
      real        z                     ! I. altitude of interest [km]
      integer     nz3,ig                   ! I. dimension of esp, ylayesp, etc...
                                        !  (=2*klev= max# of layers in ray path)
      real     iz(klev+1)              ! I. Altitude of each layer
      real*8        esp(nz3)            ! O. layer widths after geometrically 
                                        !    amplified; in [cm] except at TOA
                                        !    where an auxiliary value is used
      real*8        ilayesp(nz3)        ! O. Indexes of layers along ray path
      real*8        szalayesp(nz3)      ! O. SZA [deg]    "     "       "
      integer       nlayesp
!      real*8        nlayesp             ! O. # layers along ray path at this z
      real*8        zmini               ! O. Minimum altitud of ray path [km]


c     local variables and constants

        integer     j,i,capa
        integer     jmin                  ! index of min.altitude along ray path
        real*8      szarad                ! SZA [deg]
        real*8      zz
        real*8      diz(klev+1)
        real*8      rkmnz                 ! distance TOA to center of Planet [km]
        real*8      rkmmini               ! distance zmini to center of P [km] 
        real*8      rkmj                  ! intermediate distance to C of P [km]

c external function
        external  grid_R8          ! Returns index of layer containing the altitude
                                ! of interest, z; for example, if 
                                ! zkm(i)=z or zkm(i)<z<zkm(i+1) => grid(z)=i 
        integer   grid_R8

*************************************************************************     
        szarad = dble(szadeg)*3.141592d0/180.d0
        zz=dble(z)
        do i=1,klev
           diz(i)=dble(iz(i))
        enddo
        do j=1,nz3 
          esp(j) = 0.d0
          szalayesp(j) = 777.d0
          ilayesp(j) = 0
        enddo
        nlayesp = 0

        ! First case: szadeg<60
        ! The optical thickness will be given by  1/cos(sza)
        ! We deal with 2 different regions:
        !   1: First, all layers between z and ztop ("upper part of ray")
        !   2: Second, the layer at ztop
        if(szadeg.lt.60.d0) then

           zmini = zz
           if(abs(zz-diz(klev)).lt.1.d-3) goto 1357
           ! 1st Zone: Upper part of ray
           !
           do j=grid_R8(zz,diz,klev),klev-1
             nlayesp = nlayesp + 1 
             ilayesp(nlayesp) = j
             esp(nlayesp) = (diz(j+1)-diz(j)) / cos(szarad)        ! [km]
             esp(nlayesp) = esp(nlayesp) * 1.d5                    ! [cm]
             szalayesp(nlayesp) = szadeg
           end do

           ! 
           ! 2nd Zone: Top layer
 1357      continue
           nlayesp = nlayesp + 1 
           ilayesp(nlayesp) = klev
           esp(nlayesp) = 1.d0 / cos(szarad)         ! aux. non-dimens. factor
           szalayesp(nlayesp) = szadeg


        ! Second case:  60 < szadeg < 90
        ! The optical thickness is evaluated.
        !    (the magnitude of the effect of not using cos(sza) is 3.e-5 
        !     for z=60km & sza=30, and 5e-4 for z=60km & sza=60, approximately)
        ! We deal with 2 different regions:
        !   1: First, all layers between z and ztop ("upper part of ray")
        !   2: Second, the layer at ztop ("uppermost layer")
        else if(szadeg.le.90.d0.and.szadeg.ge.60.d0) then

           zmini=(radio+zz)*sin(szarad)-radio
           rkmmini = radio + zmini

           if(abs(zz-diz(klev)).lt.1.d-4) goto 1470

           ! 1st Zone: Upper part of ray
           !
           do j=grid_R8(zz,diz,klev),klev-1
              nlayesp = nlayesp + 1 
              ilayesp(nlayesp) = j
              esp(nlayesp) = 
     #             sqrt( (radio+diz(j+1))**2 - rkmmini**2 ) -
     #             sqrt( (radio+diz(j))**2 - rkmmini**2 )           ! [km]
              esp(nlayesp) = esp(nlayesp) * 1.d5                    ! [cm]
              rkmj = radio+diz(j)
              szalayesp(nlayesp) = asin( rkmmini/rkmj )             ! [rad]
              szalayesp(nlayesp) = szalayesp(nlayesp) * 180.d0/3.141592 ! [deg]
           end do
 1470      continue
           ! 2nd Zone:  Uppermost layer of ray.
           !
           nlayesp = nlayesp + 1 
           ilayesp(nlayesp) = klev
           rkmnz = radio+diz(klev)
           esp(nlayesp)  =  sqrt( rkmnz**2 - rkmmini**2 )       ! aux.factor[km]
           esp(nlayesp)  =  esp(nlayesp) * 1.d5                 ! aux.f. [cm]
           szalayesp(nlayesp) = asin( rkmmini/rkmnz )           ! [rad]
           szalayesp(nlayesp) = szalayesp(nlayesp) * 180.d0/3.141592! [deg]


        ! Third case:  szadeg > 90
        ! The optical thickness is evaluated.
        ! We deal with 5 different regions:
        !   1: all layers between z and ztop ("upper part of ray")
        !   2: the layer at ztop ("uppermost layer")
        !   3: the lowest layer, at zmini
        !   4: the layers increasing from zmini to z (here SZA<90)
        !   5: the layers decreasing from z to zmini (here SZA>90)
        else if(szadeg.gt.90.d0) then

           zmini=(radio+zz)*sin(szarad)-radio
           rkmmini = radio + zmini

           if(zmini.lt.diz(1)) then         ! Can see the sun?  No => esp(j)=inft
             nlayesp = nlayesp + 1 
             ilayesp(nlayesp) = - 1     ! Value to mark "no sun on view"
!             esp(nlayesp) = 1.e30

           else
              jmin=grid_R8(zmini,diz,klev)+1
              

              if(abs(zz-diz(klev)).lt.1.d-4) goto 9876

              ! 1st Zone: Upper part of ray
              !
              do j=grid_R8(zz,diz,klev),klev-1
                nlayesp = nlayesp + 1 
                ilayesp(nlayesp) = j
                esp(nlayesp) = 
     $                sqrt( (radio+diz(j+1))**2 - rkmmini**2 ) -
     $                sqrt( (radio+diz(j))**2 - rkmmini**2 )          ! [km]
                esp(nlayesp) = esp(nlayesp) * 1.d5                    ! [cm]
                rkmj = radio+diz(j)
                szalayesp(nlayesp) = asin( rkmmini/rkmj )              ! [rad]
                szalayesp(nlayesp) = szalayesp(nlayesp) *180.d0/3.141592      ! [deg]
              end do

 9876         continue
              ! 2nd Zone:  Uppermost layer of ray.
              !
              nlayesp = nlayesp + 1 
              ilayesp(nlayesp) = klev
              rkmnz = radio+diz(klev)
              esp(nlayesp) =  sqrt( rkmnz**2 - rkmmini**2 )      !aux.factor[km]
              esp(nlayesp) = esp(nlayesp) * 1.d5                 !aux.f.[cm]
              szalayesp(nlayesp) = asin( rkmmini/rkmnz )           ! [rad]
              szalayesp(nlayesp) = szalayesp(nlayesp) *180.d0/3.141592 ! [deg]

              ! 3er Zone: Lowestmost layer of ray
              !
              if ( jmin .ge. 2 ) then      ! If above the planet's surface
                j=jmin-1
                nlayesp = nlayesp + 1 
                ilayesp(nlayesp) = j
                esp(nlayesp) = 2. * 
     $                 sqrt( (radio+diz(j+1))**2 -rkmmini**2 )       ! [km]
                esp(nlayesp) = esp(nlayesp) * 1.d5                   ! [cm]
                rkmj = radio+diz(j+1)
                szalayesp(nlayesp) = asin( rkmmini/rkmj ) ! [rad]
                szalayesp(nlayesp) = szalayesp(nlayesp) *180.d0/3.141592 ! [deg]
              endif

              ! 4th zone: Lower part of ray, increasing from zmin to z
              !    ( layers with SZA < 90 deg )
              do j=jmin,grid_R8(zz,diz,klev)-1
                nlayesp = nlayesp + 1 
                ilayesp(nlayesp) = j
                esp(nlayesp) = 
     $                    sqrt( (radio+diz(j+1))**2 - rkmmini**2 )
     $                  - sqrt( (radio+diz(j))**2 - rkmmini**2 )       !  [km]
                esp(nlayesp) = esp(nlayesp) * 1.d5                     ! [cm]
                rkmj = radio+diz(j)
                szalayesp(nlayesp) = asin( rkmmini/rkmj )              ! [rad]
                szalayesp(nlayesp) = szalayesp(nlayesp) *180.d0/3.141592 ! [deg]
              end do

              ! 5th zone: Lower part of ray, decreasing from z to zmin
              !    ( layers with SZA > 90 deg )
              do j=grid_R8(zz,diz,klev)-1, jmin, -1
                nlayesp = nlayesp + 1 
                ilayesp(nlayesp) = j
                esp(nlayesp) = 
     $                    sqrt( (radio+diz(j+1))**2 - rkmmini**2 )
     $                  - sqrt( (radio+diz(j))**2 - rkmmini**2 )        ! [km]
                esp(nlayesp) = esp(nlayesp) * 1.d5                      ! [cm]
                rkmj = radio+diz(j)
                szalayesp(nlayesp) = 3.141592 - asin( rkmmini/rkmj )          ! [rad]
                szalayesp(nlayesp) = szalayesp(nlayesp)*180.d0/3.141592 ! [deg]
              end do

           end if

        end if

        return

        end


c**********************************************************************
c***********************************************************************

        function grid_R8 (z, zgrid, nz)

c Returns the index where z is located within vector zgrid
c The vector zgrid must be monotonously increasing, otherwise program stops.
c If z is outside zgrid limits, or zgrid dimension is nz<2, the program stops. 
c
c FGG     Aug-2004     Correct z.lt.zgrid(i) to .le. 
c MALV    Jul-2003
c***********************************************************************

        implicit none

c Arguments 
        integer   nz
        real*8      z
        real*8      zgrid(nz)
        integer   grid_R8

c Local  
        integer   i, nz1, nznew

c*** CODE START 

        if ( z .lt. zgrid(1) .or. z.gt.zgrid(nz) ) then 
           write (*,*) ' GRID/ z outside bounds of zgrid '
           write (*,*) ' z,zgrid(1),zgrid(nz) =', z,zgrid(1),zgrid(nz)
           stop ' Serious error in GRID.F '
        endif
        if ( nz .lt. 2 ) then 
           write (*,*) ' GRID/ zgrid needs 2 points at least ! '
           stop ' Serious error in GRID.F '
        endif
        if ( zgrid(1) .ge. zgrid(nz) ) then 
           write (*,*) ' GRID/ zgrid must increase with index'
           stop ' Serious error in GRID.F '
        endif

        nz1 = 1
        nznew = nz/2
        if ( z .gt. zgrid(nznew) ) then
           nz1 = nznew
           nznew = nz
        endif
        do i=nz1+1,nznew
           if ( z. eq. zgrid(i) ) then
              grid_R8=i
              return
              elseif ( z .le. zgrid(i) ) then
              grid_R8 = i-1
              return
           endif
        enddo
        grid_R8 = nz
        return

        end



!c***************************************************
!c***************************************************

      subroutine flujo(date)


!c     fgg           nov 2002     first version
!c***************************************************
      use dimphy
      use conc
      implicit none


!     common variables and constants
      include "param.h"
      include 'param_v4.h'
      include "clesphys.h"

!     Arguments

       real date
c       integer, parameter :: dateyr = 2006

!     Local variable and constants
!    dist_sol : distance venus - soleil

      real, parameter :: dist_sol=0.72333
      integer i
      integer inter
      real    nada

!c*************************************************

      if(date.lt.1985.) date=1985.
      if(date.gt.2001.) date=2001.
      
      do i=1,ninter
         fluxtop(i)=1.
         !Variation of solar flux with 11 years solar cycle
         !For more details, see Gonzalez-Galindo et al. 2005
         !To be improved in next versions
        if(i.le.24 .and.solvarmod.eq.0) then
          fluxtop(i)=(((ct1(i)+p1(i)*date)/2.)                  
     $        *sin(2.*3.1416/11.*(date-1985.-3.1416))          
     $        +(ct2(i)+p2(i)*date)+1.)*fluxtop(i)
          
        end if
 !     The solar flux calculated  
 !      is  corrected for
 !     the actual Venus-Sun dist
        fluxtop(i)=fluxtop(i)*(1./dist_sol)**2


       end do
     
      return
      end
