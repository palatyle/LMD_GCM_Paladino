c**********************************************************************

      subroutine jthermcalc_e107
     $     (ig,nlayer,chemthermod,rm,nesptherm,tx,iz,zenit,zday)


c     feb 2002        fgg           first version
c     nov 2002        fgg           second version
c
c modified from paramhr.F
c MAC July 2003
c**********************************************************************

      use param_v4_h, only: jfotsout,crscabsi2,
     .    c1_16,c17_24,c25_29,c30_31,c32,c33,c34,c35,c36,
     .    co2crsc195,co2crsc295,t0,
     .    jabsifotsintpar,ninter,nz2,
     .    nabs,e107,date_e107,e107_tab,
     .    coefit0,coefit1,coefit2,coefit3,coefit4

      implicit none

      include "callkeys.h"

c     input and output variables

      integer    ig,nlayer
      integer    chemthermod
      integer    nesptherm                      !Number of species considered
      real       rm(nlayer,nesptherm)         !Densities (cm-3)
      real       tx(nlayer)                   !temperature
      real       zenit                          !SZA
      real       iz(nlayer)                   !Local altitude
      real       zday                           !Martian day after Ls=0


c    local parameters and variables

      real       co2colx(nlayer)              !column density of CO2 (cm^-2)
      real       o2colx(nlayer)               !column density of O2(cm^-2)
      real       o3pcolx(nlayer)              !column density of O(3P)(cm^-2)
      real       h2colx(nlayer)               !H2 column density (cm-2)
      real       h2ocolx(nlayer)              !H2O column density (cm-2)
      real       h2o2colx(nlayer)             !column density of H2O2(cm^-2)
      real       o3colx(nlayer)               !O3 column density (cm-2)
      real       n2colx(nlayer)               !N2 column density (cm-2)
      real       ncolx(nlayer)                !N column density (cm-2)
      real       nocolx(nlayer)               !NO column density (cm-2)
      real       cocolx(nlayer)               !CO column density (cm-2)
      real       hcolx(nlayer)                !H column density (cm-2)
      real       no2colx(nlayer)              !NO2 column density (cm-2)
      real       t2(nlayer)
      real       coltemp(nlayer)
      real       sigma(ninter,nlayer)
      real       alfa(ninter,nlayer)
      real       realday
      
      integer    i,j,k,indexint                 !indexes
      character  dn
      integer    tinf,tsup



c     variables used in interpolation

      real*8      auxcoltab(nz2)
      real*8      auxjco2(nz2)
      real*8      auxjo2(nz2)
      real*8      auxjo3p(nz2)
      real*8      auxjh2o(nz2)
      real*8      auxjh2(nz2)
      real*8      auxjh2o2(nz2)
      real*8      auxjo3(nz2)
      real*8      auxjn2(nz2)
      real*8      auxjn(nz2)
      real*8      auxjno(nz2)
      real*8      auxjco(nz2)
      real*8      auxjh(nz2)
      real*8      auxjno2(nz2)
      real*8      wp(nlayer),wm(nlayer)
      real*8      auxcolinp(nlayer)
      integer     auxind(nlayer)
      integer     auxi
      integer     ind
      real*8      cortemp(nlayer)

      real*8      limdown                      !limits for interpolation
      real*8      limup                        !        ""

      !!!ATTENTION. Here i_co2 has to have the same value than in chemthermos.F90
      !!!If the value is changed there, if has to be changed also here !!!!
      integer,parameter :: i_co2=1


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
      do i=1,nlayer
         t2(i)=tx(i)
         if(t2(i).lt.195.0) t2(i)=195.0
         if(t2(i).gt.295.0) t2(i)=295.0
      end do

      !Calculation of column amounts 
      call column(ig,nlayer,chemthermod,rm,nesptherm,tx,iz,zenit,
     $     co2colx,o2colx,o3pcolx,h2colx,h2ocolx,
     $     h2o2colx,o3colx,n2colx,ncolx,nocolx,cocolx,hcolx,no2colx)

      !Auxiliar column to include the temperature dependence 
      !of CO2 cross section
      coltemp(nlayer)=co2colx(nlayer)*abs(t2(nlayer)-t0(nlayer))
      do i=nlayer-1,1,-1
        coltemp(i)=!coltemp(i+1)+     PQ SE ELIMINA? REVISAR 
     $         ( rm(i,i_co2) + rm(i+1,i_co2) ) * 0.5 
     $         * 1e5 * (iz(i+1)-iz(i)) * abs(t2(i)-t0(i))
      end do
      
      !Calculation of CO2 cross section at temperature t0(i)
      do i=1,nlayer
         do indexint=24,32
           sigma(indexint,i)=co2crsc195(indexint-23)
           alfa(indexint,i)=((co2crsc295(indexint-23)
     $          /sigma(indexint,i))-1.)/(295.-t0(i))
        end do
      end do

      if (solvarmod==0) then
        e107=fixed_euv_value
      else
        !E10.7 for the day: linear interpolation to tabulated values
        realday=mod(zday,669.)
        if(realday.lt.date_e107(1)) then
         e107=e107_tab(1)
        else if(realday.ge.date_e107(669)) then
         e107=e107_tab(669)   
        else if(realday.ge.date_e107(1).and.
     $        realday.lt.date_e107(669)) then
         do i=1,668
            if(realday.ge.date_e107(i).and.
     $           realday.lt.date_e107(i+1)) then
               tinf=i
               tsup=i+1
               e107=e107_tab(tinf)+(realday-date_e107(tinf))*
     $              (e107_tab(tsup)-e107_tab(tinf))
            endif
         enddo
        endif
      endif ! of if (solvarmod==0)

      !Photoabsorption coefficients at TOA as a function of E10.7
      do j=1,nabs
         do indexint=1,ninter
            jfotsout(indexint,j,nlayer)=coefit0(indexint,j)+
     $           coefit1(indexint,j)*e107+coefit2(indexint,j)*e107**2+
     $           coefit3(indexint,j)*e107**3+coefit4(indexint,j)*e107**4
         enddo
      enddo
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
      do i=1,nlayer
         auxcolinp(nlayer-i+1) = co2colx(i)*crscabsi2(1,indexint) +
     $        o2colx(i)*crscabsi2(2,indexint) + 
     $        o3pcolx(i)*crscabsi2(3,indexint) + 
     $        h2colx(i)*crscabsi2(5,indexint) + 
     $        ncolx(i)*crscabsi2(9,indexint)
      end do
      limdown=1.e-20
      limup=1.e26


c     Interpolations

      do i=1,nz2
         auxi = nz2-i+1
         !CO2 tabulated coefficient
         auxjco2(i) = jabsifotsintpar(auxi,1,indexint)
         !O2 tabulated coefficient
         auxjo2(i) = jabsifotsintpar(auxi,2,indexint)
         !O3p tabulated coefficient
         auxjo3p(i) = jabsifotsintpar(auxi,3,indexint)
         !H2 tabulated coefficient
         auxjh2(i) = jabsifotsintpar(auxi,5,indexint)
         !Tabulated column
         auxcoltab(i) = c1_16(auxi,indexint)
      enddo
      !Only if chemthermod.ge.2
      !N tabulated coefficient
      if(chemthermod.ge.2) then
         do i=1,nz2
            auxjn(i) = jabsifotsintpar(nz2-i+1,9,indexint)
         enddo
      endif

      call interfast 
     $     (wm,wp,auxind,auxcolinp,nlayer,auxcoltab,nz2,limdown,limup)
      do i=1,nlayer
         ind=auxind(i)
         auxi=nlayer-i+1
         !CO2 interpolated coefficient
         jfotsout(indexint,1,auxi) = jfotsout(indexint,1,nlayer) *
     $        (wm(i)*auxjco2(ind+1) + wp(i)*auxjco2(ind))
         !O2 interpolated coefficient
         jfotsout(indexint,2,auxi) = jfotsout(indexint,2,nlayer) *
     $        (wm(i)*auxjo2(ind+1) + wp(i)*auxjo2(ind))
         !O3p interpolated coefficient
         jfotsout(indexint,3,auxi) = jfotsout(indexint,3,nlayer) *
     $        (wm(i)*auxjo3p(ind+1) + wp(i)*auxjo3p(ind))
         !H2 interpolated coefficient
         jfotsout(indexint,5,auxi) = jfotsout(indexint,5,auxi) *
     $        (wm(i)*auxjh2(ind+1) + wp(i)*auxjh2(ind))
      enddo
      !Only if chemthermod.ge.2
      !N interpolated coefficient
      if(chemthermod.ge.2) then
         do i=1,nlayer
            ind=auxind(i)
            jfotsout(indexint,9,nlayer-i+1) =  
     $           jfotsout(indexint,9,nlayer) * 
     $           (wm(i)*auxjn(ind+1) + wp(i)*auxjn(ind))
         enddo
      endif
         

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
         do i=1,nlayer
            auxcolinp(nlayer-i+1) = co2colx(i)*crscabsi2(1,indexint)+
     $           o2colx(i)*crscabsi2(2,indexint)+
     $           o3pcolx(i)*crscabsi2(3,indexint)+
     $           h2colx(i)*crscabsi2(5,indexint)+
     $           n2colx(i)*crscabsi2(8,indexint)+
     $           ncolx(i)*crscabsi2(9,indexint)+
     $           nocolx(i)*crscabsi2(10,indexint)+
     $           cocolx(i)*crscabsi2(11,indexint)+
     $           hcolx(i)*crscabsi2(12,indexint)+
     $           no2colx(i)*crscabsi2(13,indexint)
         end do

c     Interpolations

         do i=1,nz2
            auxi = nz2-i+1
            !O2 tabulated coefficient
            auxjo2(i) = jabsifotsintpar(auxi,2,indexint)
            !O3p tabulated coefficient
            auxjo3p(i) = jabsifotsintpar(auxi,3,indexint)
            !CO2 tabulated coefficient
            auxjco2(i) = jabsifotsintpar(auxi,1,indexint)
            !H2 tabulated coefficient
            auxjh2(i) = jabsifotsintpar(auxi,5,indexint)
            !N2 tabulated coefficient
            auxjn2(i) = jabsifotsintpar(auxi,8,indexint)
            !CO tabulated coefficient
            auxjco(i) = jabsifotsintpar(auxi,11,indexint)
            !H tabulated coefficient
            auxjh(i) = jabsifotsintpar(auxi,12,indexint)
            !tabulated column
            auxcoltab(i) = c1_16(auxi,indexint)
         enddo
         !Only if chemthermod.ge.2
         if(chemthermod.ge.2) then
            do i=1,nz2
               auxi = nz2-i+1
               !N tabulated coefficient
               auxjn(i) = jabsifotsintpar(auxi,9,indexint)
               !NO tabulated coefficient
               auxjno(i) = jabsifotsintpar(auxi,10,indexint)
               !NO2 tabulated coefficient
               auxjno2(i) = jabsifotsintpar(auxi,13,indexint)
            enddo
         endif

          call interfast(wm,wp,auxind,auxcolinp,nlayer,
     $        auxcoltab,nz2,limdown,limup)
          do i=1,nlayer
             ind=auxind(i)
             auxi = nlayer-i+1
             !O2 interpolated coefficient
             jfotsout(indexint,2,auxi) = 
     $            jfotsout(indexint,2,nlayer) *
     $            (wm(i)*auxjo2(ind+1) + wp(i)*auxjo2(ind))
             !O3p interpolated coefficient
             jfotsout(indexint,3,auxi) = 
     $            jfotsout(indexint,3,nlayer) *
     $            (wm(i)*auxjo3p(ind+1) + wp(i)*auxjo3p(ind))
             !CO2 interpolated coefficient
             jfotsout(indexint,1,auxi) = 
     $            jfotsout(indexint,1,nlayer) *
     $            (wm(i)*auxjco2(ind+1) + wp(i)*auxjco2(ind))
             !H2 interpolated coefficient
             jfotsout(indexint,5,auxi) = 
     $            jfotsout(indexint,5,nlayer) *
     $            (wm(i)*auxjh2(ind+1) + wp(i)*auxjh2(ind))
             !N2 interpolated coefficient
             jfotsout(indexint,8,auxi) = 
     $            jfotsout(indexint,8,nlayer) *
     $            (wm(i)*auxjn2(ind+1) + wp(i)*auxjn2(ind))
             !CO interpolated coefficient
             jfotsout(indexint,11,auxi) = 
     $            jfotsout(indexint,11,nlayer) *
     $            (wm(i)*auxjco(ind+1) + wp(i)*auxjco(ind))
             !H interpolated coefficient
             jfotsout(indexint,12,auxi) = 
     $            jfotsout(indexint,12,nlayer) *
     $            (wm(i)*auxjh(ind+1) + wp(i)*auxjh(ind))
          enddo
          !Only if chemthermod.ge.2
          if(chemthermod.ge.2) then
             do i=1,nlayer
                ind=auxind(i)
                auxi = nlayer-i+1
                !N interpolated coefficient
                jfotsout(indexint,9,auxi) = 
     $               jfotsout(indexint,9,nlayer) *
     $               (wm(i)*auxjn(ind+1) + wp(i)*auxjn(ind))
                !NO interpolated coefficient
                jfotsout(indexint,10,auxi)=
     $               jfotsout(indexint,10,nlayer) *
     $               (wm(i)*auxjno(ind+1) + wp(i)*auxjno(ind))
                !NO2 interpolated coefficient
                jfotsout(indexint,13,auxi)=
     $               jfotsout(indexint,13,nlayer) *
     $               (wm(i)*auxjno2(ind+1) + wp(i)*auxjno2(ind))
             enddo
          endif   
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
      do i=1,nlayer
         auxcolinp(nlayer-i+1) = co2colx(i)*crscabsi2(1,indexint)+
     $        o2colx(i)*crscabsi2(2,indexint)+
     $        o3pcolx(i)*crscabsi2(3,indexint)+
     $        n2colx(i)*crscabsi2(8,indexint)+
     $        ncolx(i)*crscabsi2(9,indexint)+
     $        nocolx(i)*crscabsi2(10,indexint)+
     $        cocolx(i)*crscabsi2(11,indexint)+
     $        hcolx(i)*crscabsi2(12,indexint)+
     $        no2colx(i)*crscabsi2(13,indexint)
      end do

c     Interpolations

      do i=1,nz2
         auxi = nz2-i+1
         !O2 tabulated coefficient
         auxjo2(i) = jabsifotsintpar(auxi,2,indexint)
         !CO2 tabulated coefficient
         auxjco2(i) = jabsifotsintpar(auxi,1,indexint)
         !O3p tabulated coefficient
         auxjo3p(i) = jabsifotsintpar(auxi,3,indexint)
         !N2 tabulated coefficient
         auxjn2(i) = jabsifotsintpar(auxi,8,indexint)
         !CO tabulated coefficient
         auxjco(i) = jabsifotsintpar(auxi,11,indexint)
         !H tabulated coefficient
         auxjh(i) = jabsifotsintpar(auxi,12,indexint)
         !NO2 tabulated coefficient
         auxjno2(i) = jabsifotsintpar(auxi,13,indexint)
         !Tabulated column
         auxcoltab(i) = c1_16(auxi,indexint)
      enddo
      !Only if chemthermod.ge.2
      if(chemthermod.ge.2) then
         do i=1,nz2
            auxi = nz2-i+1
            !N tabulated coefficient
            auxjn(i) = jabsifotsintpar(auxi,9,indexint)
            !NO tabulated coefficient
            auxjno(i) = jabsifotsintpar(auxi,10,indexint)
            !NO2 tabulated coefficient
            auxjno2(i) = jabsifotsintpar(auxi,13,indexint)
         enddo
      endif

      call interfast
     $     (wm,wp,auxind,auxcolinp,nlayer,auxcoltab,nz2,limdown,limup)
      do i=1,nlayer
         ind=auxind(i)
         auxi = nlayer-i+1
         !O2 interpolated coefficient
         jfotsout(indexint,2,auxi) = jfotsout(indexint,2,nlayer) *
     $        (wm(i)*auxjo2(ind+1) + wp(i)*auxjo2(ind))
         !CO2 interpolated coefficient
         jfotsout(indexint,1,auxi) = jfotsout(indexint,1,nlayer) *
     $        (wm(i)*auxjco2(ind+1) + wp(i)*auxjco2(ind))
         !O3p interpolated coefficient
         jfotsout(indexint,3,auxi) = jfotsout(indexint,3,nlayer) *
     $        (wm(i)*auxjo3p(ind+1) + wp(i)*auxjo3p(ind))
         !N2 interpolated coefficient
         jfotsout(indexint,8,auxi) = jfotsout(indexint,8,nlayer) *
     $        (wm(i)*auxjn2(ind+1) + wp(i)*auxjn2(ind))
         !CO interpolated coefficient
         jfotsout(indexint,11,auxi) = 
     $        jfotsout(indexint,11,nlayer) *
     $        (wm(i)*auxjco(ind+1) + wp(i)*auxjco(ind))
         !H interpolated coefficient
         jfotsout(indexint,12,auxi) = 
     $        jfotsout(indexint,12,nlayer) *
     $        (wm(i)*auxjh(ind+1) + wp(i)*auxjh(ind))
      enddo
      !Only if chemthermod.ge.2
      if(chemthermod.ge.2) then
         do i=1,nlayer
            ind=auxind(i)
            auxi = nlayer-i+1
            !N interpolated coefficient
            jfotsout(indexint,9,auxi) = 
     $           jfotsout(indexint,9,nlayer) *
     $           (wm(i)*auxjn(ind+1) + wp(i)*auxjn(ind))
            !NO interpolated coefficient
            jfotsout(indexint,10,auxi) = 
     $           jfotsout(indexint,10,nlayer) * 
     $           (wm(i)*auxjno(ind+1) + wp(i)*auxjno(ind))
            !NO2 interpolated coefficient
            jfotsout(indexint,13,auxi) = 
     $           jfotsout(indexint,13,nlayer) *
     $           (wm(i)*auxjno2(ind+1) + wp(i)*auxjno2(ind))
         enddo
      endif
c     End interval 16


ccccccccccccccccccccccccccccccc
c     90.9-119.5nm (int 17-24)
c
c     Absorption by:
c     CO2, O2, N2, NO, CO, NO2
ccccccccccccccccccccccccccccccc

c     Input column

      do i=1,nlayer
         auxcolinp(nlayer-i+1) = co2colx(i) + o2colx(i) + n2colx(i) +
     $        nocolx(i) + cocolx(i) + no2colx(i)
      end do

      do indexint=17,24

c     Interpolations

         do i=1,nz2
            auxi = nz2-i+1
            !CO2 tabulated coefficient
            auxjco2(i) = jabsifotsintpar(auxi,1,indexint)
            !O2 tabulated coefficient
            auxjo2(i) = jabsifotsintpar(auxi,2,indexint)
            !N2 tabulated coefficient
            auxjn2(i) = jabsifotsintpar(auxi,8,indexint)
            !CO tabulated coefficient
            auxjco(i) = jabsifotsintpar(auxi,11,indexint)            
            !Tabulated column
            auxcoltab(i) = c17_24(auxi)
         enddo
         !Only if chemthermod.ge.2
         if(chemthermod.ge.2) then
            do i=1,nz2
               auxi = nz2-i+1
               !NO tabulated coefficient
               auxjno(i) = jabsifotsintpar(auxi,10,indexint)
               !NO2 tabulated coefficient
               auxjno2(i) = jabsifotsintpar(auxi,13,indexint)
            enddo
         endif

         call interfast
     $     (wm,wp,auxind,auxcolinp,nlayer,auxcoltab,nz2,limdown,limup)
         !Correction to include T variation of CO2 cross section
         if(indexint.eq.24) then
            do i=1,nlayer
               auxi = nlayer-i+1
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
            do i=1,nlayer
               cortemp(i)=1.
            enddo
         end if
         do i=1,nlayer           
            ind=auxind(i)
            auxi = nlayer-i+1
            !O2 interpolated coefficient
            jfotsout(indexint,2,auxi) = 
     $           jfotsout(indexint,2,nlayer) *
     $           (wm(i)*auxjo2(ind+1) + wp(i)*auxjo2(ind)) * 
     $           cortemp(i)
            !CO2 interpolated coefficient
            jfotsout(indexint,1,auxi) = 
     $           jfotsout(indexint,1,nlayer) *
     $           (wm(i)*auxjco2(ind+1) + wp(i)*auxjco2(ind)) 
     $           * cortemp(i)
            if(indexint.eq.24) jfotsout(indexint,1,auxi)=
     $           jfotsout(indexint,1,auxi)*
     $           (1+alfa(indexint,auxi)*
     $           (t2(auxi)-t0(auxi)))
            !N2 interpolated coefficient
            jfotsout(indexint,8,auxi) = 
     $           jfotsout(indexint,8,nlayer) *
     $           (wm(i)*auxjn2(ind+1) + wp(i)*auxjn2(ind)) * 
     $           cortemp(i)            
            !CO interpolated coefficient
            jfotsout(indexint,11,auxi) = 
     $           jfotsout(indexint,11,nlayer) *
     $           (wm(i)*auxjco(ind+1) + wp(i)*auxjco(ind)) * 
     $           cortemp(i)            
         enddo
         !Only if chemthermod.ge.2
         if(chemthermod.ge.2) then
            do i=1,nlayer
               ind=auxind(i)
               auxi = nlayer-i+1
               !NO interpolated coefficient
               jfotsout(indexint,10,auxi)=
     $              jfotsout(indexint,10,nlayer) *
     $              (wm(i)*auxjno(ind+1) + wp(i)*auxjno(ind)) * 
     $              cortemp(i)
               !NO2 interpolated coefficient
               jfotsout(indexint,13,auxi)=
     $              jfotsout(indexint,13,nlayer) *
     $              (wm(i)*auxjno2(ind+1)+ wp(i)*auxjno2(ind)) * 
     $              cortemp(i)
            enddo
         endif               
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

      do i=1,nlayer
         auxcolinp(nlayer-i+1) = co2colx(i) + o2colx(i) + h2ocolx(i) + 
     $        h2o2colx(i) + nocolx(i) + cocolx(i) + no2colx(i)
      end do

      do indexint=25,29

c     Interpolations

         do i=1,nz2
            auxi = nz2-i+1
            !CO2 tabulated coefficient
            auxjco2(i) = jabsifotsintpar(auxi,1,indexint)
            !O2 tabulated coefficient
            auxjo2(i) = jabsifotsintpar(auxi,2,indexint)
            !H2O tabulated coefficient
            auxjh2o(i) = jabsifotsintpar(auxi,4,indexint)
            !H2O2 tabulated coefficient
            auxjh2o2(i) = jabsifotsintpar(auxi,6,indexint)            
            !CO tabulated coefficient
            auxjco(i) = jabsifotsintpar(auxi,11,indexint)            
            !Tabulated column
            auxcoltab(i) = c25_29(auxi)
         enddo
         !Only if chemthermod.ge.2
         if(chemthermod.ge.2) then
            do i=1,nz2
               auxi = nz2-i+1
               !NO tabulated coefficient
               auxjno(i) = jabsifotsintpar(auxi,10,indexint)
               !NO2 tabulated coefficient
               auxjno2(i) = jabsifotsintpar(auxi,13,indexint)
            enddo
         endif
         call interfast
     $     (wm,wp,auxind,auxcolinp,nlayer,auxcoltab,nz2,limdown,limup)
         do i=1,nlayer
            ind=auxind(i)
            auxi = nlayer-i+1
            !Correction to include T variation of CO2 cross section
            if(sigma(indexint,auxi)*alfa(indexint,auxi)*
     $           coltemp(auxi).lt.60.) then
               cortemp(i)=exp(-sigma(indexint,auxi)*
     $              alfa(indexint,auxi)*coltemp(auxi))
            else 
               cortemp(i)=0.
            end if
            !CO2 interpolated coefficient
            jfotsout(indexint,1,auxi) = 
     $           jfotsout(indexint,1,nlayer) *
     $           (wm(i)*auxjco2(ind+1) + wp(i)*auxjco2(ind)) * 
     $           cortemp(i) *
     $           (1+alfa(indexint,auxi)*
     $           (t2(auxi)-t0(auxi)))
            !O2 interpolated coefficient
            jfotsout(indexint,2,auxi) = 
     $           jfotsout(indexint,2,nlayer) *
     $           (wm(i)*auxjo2(ind+1) + wp(i)*auxjo2(ind)) * 
     $           cortemp(i)
            !H2O interpolated coefficient
            jfotsout(indexint,4,auxi) = 
     $           jfotsout(indexint,4,nlayer) *
     $           (wm(i)*auxjh2o(ind+1) + wp(i)*auxjh2o(ind)) * 
     $           cortemp(i)
            !H2O2 interpolated coefficient
            jfotsout(indexint,6,auxi) = 
     $           jfotsout(indexint,6,nlayer) *
     $           (wm(i)*auxjh2o2(ind+1) + wp(i)*auxjh2o2(ind)) * 
     $           cortemp(i)            
            !CO interpolated coefficient
            jfotsout(indexint,11,auxi) = 
     $           jfotsout(indexint,11,nlayer) *
     $           (wm(i)*auxjco(ind+1) + wp(i)*auxjco(ind)) * 
     $           cortemp(i)
         enddo
         !Only if chemthermod.ge.2
         if(chemthermod.ge.2) then
            do i=1,nlayer
               ind=auxind(i)
               auxi = nlayer-i+1
               !NO interpolated coefficient
               jfotsout(indexint,10,auxi)=
     $              jfotsout(indexint,10,nlayer) *
     $              (wm(i)*auxjno(ind+1) + wp(i)*auxjno(ind)) * 
     $              cortemp(i)
               !NO2 interpolated coefficient
               jfotsout(indexint,13,auxi)=
     $              jfotsout(indexint,13,nlayer) *
     $              (wm(i)*auxjno2(ind+1) + wp(i)*auxjno2(ind)) * 
     $              cortemp(i)
            enddo
         endif

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

      do i=1,nlayer
         auxcolinp(nlayer-i+1) = co2colx(i) + o2colx(i) + h2ocolx(i) + 
     $        h2o2colx(i) + nocolx(i) + no2colx(i)
      end do

c     Interpolation

      do indexint=30,31

         do i=1,nz2
            auxi = nz2-i+1
            !CO2 tabulated coefficient
            auxjco2(i) = jabsifotsintpar(auxi,1,indexint)
            !O2 tabulated coefficient
            auxjo2(i) = jabsifotsintpar(auxi,2,indexint)
            !H2O tabulated coefficient
            auxjh2o(i) = jabsifotsintpar(auxi,4,indexint)
            !H2O2 tabulated coefficient
            auxjh2o2(i) = jabsifotsintpar(auxi,6,indexint)            
            !Tabulated column
            auxcoltab(i) = c30_31(auxi)
         enddo
         !Only if chemthermod.ge.2
         if(chemthermod.ge.2) then
            do i=1,nz2
               auxi = nz2-i+1
               !NO tabulated coefficient
               auxjno(i) = jabsifotsintpar(auxi,10,indexint)
               !NO2 tabulated coefficient
               auxjno2(i) = jabsifotsintpar(auxi,13,indexint)
            enddo
         endif

         call interfast
     $     (wm,wp,auxind,auxcolinp,nlayer,auxcoltab,nz2,limdown,limup)
         do i=1,nlayer
            ind=auxind(i)
            auxi = nlayer-i+1
            !Correction to include T variation of CO2 cross section
            if(sigma(indexint,auxi)*alfa(indexint,auxi)*
     $           coltemp(auxi).lt.60.) then
               cortemp(i)=exp(-sigma(indexint,auxi)*
     $              alfa(indexint,auxi)*coltemp(auxi))
            else 
               cortemp(i)=0.
            end if
            !CO2 interpolated coefficient
            jfotsout(indexint,1,auxi) = 
     $           jfotsout(indexint,1,nlayer) *
     $           (wm(i)*auxjco2(ind+1) + wp(i)*auxjco2(ind)) * 
     $           cortemp(i) *
     $           (1+alfa(indexint,auxi)*
     $           (t2(auxi)-t0(auxi)))
            !O2 interpolated coefficient
            jfotsout(indexint,2,auxi) = 
     $           jfotsout(indexint,2,nlayer) *
     $           (wm(i)*auxjo2(ind+1) + wp(i)*auxjo2(ind)) * 
     $           cortemp(i)
            !H2O interpolated coefficient
            jfotsout(indexint,4,auxi) = 
     $           jfotsout(indexint,4,nlayer) *
     $           (wm(i)*auxjh2o(ind+1) + wp(i)*auxjh2o(ind)) * 
     $           cortemp(i)
            !H2O2 interpolated coefficient
            jfotsout(indexint,6,auxi) = 
     $           jfotsout(indexint,6,nlayer) *
     $           (wm(i)*auxjh2o2(ind+1) + wp(i)*auxjh2o2(ind)) * 
     $           cortemp(i)            
         enddo
         !Only if chemthermod.ge.2
         if(chemthermod.ge.2) then
            do i=1,nlayer 
               ind=auxind(i)
               auxi = nlayer-i+1
               !NO interpolated coefficient
               jfotsout(indexint,10,auxi)=
     $              jfotsout(indexint,10,nlayer) *
     $              (wm(i)*auxjno(ind+1) +wp(i)*auxjno(ind)) * 
     $              cortemp(i)
               !NO2 interpolated coefficient
               jfotsout(indexint,13,auxi)=
     $              jfotsout(indexint,13,auxi) *
     $              (wm(i)*auxjno2(ind+1)+wp(i)*auxjno2(ind)) * 
     $              cortemp(i)
            enddo
         endif

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
      do i=1,nlayer
         auxcolinp(nlayer-i+1) =co2colx(i) + o2colx(i) + h2o2colx(i) + 
     $        nocolx(i) + no2colx(i)
      end do

c     Interpolation

      do i=1,nz2
         auxi = nz2-i+1
         !CO2 tabulated coefficient
         auxjco2(i) = jabsifotsintpar(auxi,1,indexint)
         !O2 tabulated coefficient
         auxjo2(i) = jabsifotsintpar(auxi,2,indexint)
         !H2O2 tabulated coefficient
         auxjh2o2(i) = jabsifotsintpar(auxi,6,indexint)         
         !Tabulated column
         auxcoltab(i) = c32(auxi)
      enddo
      !Only if chemthermod.ge.2
      if(chemthermod.ge.2) then
         do i=1,nz2
            auxi = nz2-i+1
            !NO tabulated coefficient
            auxjno(i) = jabsifotsintpar(auxi,10,indexint)
            !NO2 tabulated coefficient
            auxjno2(i) = jabsifotsintpar(auxi,13,indexint)
         enddo
      endif
      call interfast
     $     (wm,wp,auxind,auxcolinp,nlayer,auxcoltab,nz2,limdown,limup)
      do i=1,nlayer
         ind=auxind(i)
         auxi = nlayer-i+1
         !Correction to include T variation of CO2 cross section
         if(sigma(indexint,nlayer-i+1)*alfa(indexint,auxi)*
     $        coltemp(auxi).lt.60.) then
            cortemp(i)=exp(-sigma(indexint,auxi)*
     $           alfa(indexint,auxi)*coltemp(auxi))
         else 
            cortemp(i)=0.
         end if
         !CO2 interpolated coefficient
         jfotsout(indexint,1,auxi) = 
     $        jfotsout(indexint,1,nlayer) *
     $        (wm(i)*auxjco2(ind+1)+wp(i)*auxjco2(ind)) * 
     $        cortemp(i) *
     $        (1+alfa(indexint,auxi)*
     $        (t2(auxi)-t0(auxi)))
         !O2 interpolated coefficient
         jfotsout(indexint,2,auxi) = 
     $        jfotsout(indexint,2,nlayer) *
     $        (wm(i)*auxjo2(ind+1) + wp(i)*auxjo2(ind)) * 
     $        cortemp(i)
         !H2O2 interpolated coefficient
         jfotsout(indexint,6,auxi) = 
     $        jfotsout(indexint,6,nlayer) *
     $        (wm(i)*auxjh2o2(ind+1) + wp(i)*auxjh2o2(ind)) * 
     $        cortemp(i)         
      enddo
      !Only if chemthermod.ge.2
      if(chemthermod.ge.2) then
         do i=1,nlayer
            auxi = nlayer-i+1
            ind=auxind(i)
            !NO interpolated coefficient
            jfotsout(indexint,10,auxi) = 
     $           jfotsout(indexint,10,nlayer) * 
     $           (wm(i)*auxjno(ind+1) + wp(i)*auxjno(ind)) * 
     $           cortemp(i)
           !NO2 interpolated coefficient
            jfotsout(indexint,13,auxi) = 
     $           jfotsout(indexint,13,nlayer) *
     $           (wm(i)*auxjno2(ind+1) + wp(i)*auxjno2(ind)) * 
     $           cortemp(i)
         enddo
      endif

c     End of interval 32


ccccccccccccccccccccccccccccccc
c     210.1-231.0nm (int 33)
c     
c     Absorption by:
c     O2, H2O2, NO2
ccccccccccccccccccccccccccccccc

c     Input atmospheric column

      indexint=33
      do i=1,nlayer
         auxcolinp(nlayer-i+1) = o2colx(i) + h2o2colx(i) + no2colx(i)
      end do

c     Interpolation

      do i=1,nz2
         auxi = nz2-i+1
         !O2 tabulated coefficient
         auxjo2(i) = jabsifotsintpar(auxi,2,indexint)
         !H2O2 tabulated coefficient
         auxjh2o2(i) = jabsifotsintpar(auxi,6,indexint)
         !Tabulated column
         auxcoltab(i) = c33(auxi)
      enddo
      !Only if chemthermod.ge.2
      if(chemthermod.ge.2) then
         do i=1,nz2
            !NO2 tabulated coefficient
            auxjno2(i) = jabsifotsintpar(nz2-i+1,13,indexint)
         enddo
      endif
      call interfast
     $     (wm,wp,auxind,auxcolinp,nlayer,auxcoltab,nz2,limdown,limup)
      do i=1,nlayer
         ind=auxind(i)
         auxi = nlayer-i+1
         !O2 interpolated coefficient
         jfotsout(indexint,2,auxi) = jfotsout(indexint,2,nlayer) *
     $        (wm(i)*auxjo2(ind+1) + wp(i)*auxjo2(ind))
         !H2O2 interpolated coefficient
         jfotsout(indexint,6,auxi) = jfotsout(indexint,6,nlayer) *
     $        (wm(i)*auxjh2o2(ind+1) + wp(i)*auxjh2o2(ind))
      enddo
      !Only if chemthermod.ge.2
      if(chemthermod.ge.2) then
         do i=1,nlayer
            ind=auxind(i)
            !NO2 interpolated coefficient
            jfotsout(indexint,13,nlayer-i+1) = 
     $           jfotsout(indexint,13,nlayer) *
     $           (wm(i)*auxjno2(ind+1) + wp(i)*auxjno2(ind))
         enddo
      endif

c     End of interval 33


ccccccccccccccccccccccccccccccc
c     231.1-240.0nm (int 34)
c
c     Absorption by:
c     O2, H2O2, O3, NO2
ccccccccccccccccccccccccccccccc

c     Input atmospheric column

      indexint=34
      do i=1,nlayer
         auxcolinp(nlayer-i+1) = h2o2colx(i) + o2colx(i) + o3colx(i) + 
     $        no2colx(i)
      end do

c     Interpolation

      do i=1,nz2
         auxi = nz2-i+1
         !O2 tabulated coefficient
         auxjo2(i) = jabsifotsintpar(auxi,2,indexint)
         !H2O2 tabulated coefficient
         auxjh2o2(i) = jabsifotsintpar(auxi,6,indexint)
         !O3 tabulated coefficient
         auxjo3(i) = jabsifotsintpar(auxi,7,indexint)         
         !Tabulated column
         auxcoltab(i) = c34(nz2-i+1)
      enddo
      !Only if chemthermod.ge.2
      if(chemthermod.ge.2) then
         do i=1,nz2
            !NO2 tabulated coefficient
            auxjno2(i) = jabsifotsintpar(nz2-i+1,13,indexint)
         enddo
      endif
      call interfast
     $     (wm,wp,auxind,auxcolinp,nlayer,auxcoltab,nz2,limdown,limup)
      do i=1,nlayer
         ind=auxind(i)
         auxi = nlayer-i+1
         !O2 interpolated coefficient
         jfotsout(indexint,2,auxi) = jfotsout(indexint,2,nlayer) * 
     $        (wm(i)*auxjo2(ind+1) + wp(i)*auxjo2(ind))
         !H2O2 interpolated coefficient
         jfotsout(indexint,6,auxi) = jfotsout(indexint,6,nlayer) *
     $        (wm(i)*auxjh2o2(ind+1) + wp(i)*auxjh2o2(ind))
         !O3 interpolated coefficient
         jfotsout(indexint,7,auxi) = jfotsout(indexint,7,nlayer) *
     $        (wm(i)*auxjo3(ind+1) + wp(i)*auxjo3(ind))
      enddo
      !Only if chemthermod.ge.2
      if(chemthermod.ge.2) then
         do i=1,nlayer
            ind=auxind(i)
            !NO2 interpolated coefficient
            jfotsout(indexint,13,nlayer-i+1) = 
     $           jfotsout(indexint,13,nlayer) *
     $           (wm(i)*auxjno2(ind+1) + wp(i)*auxjno2(ind))
         enddo
      endif

c     End of interval 34      


ccccccccccccccccccccccccccccccc
c     240.1-337.7nm (int 35)
c
c     Absorption by:
c     H2O2, O3, NO2
ccccccccccccccccccccccccccccccc

c     Input atmospheric column

      indexint=35
      do i=1,nlayer
         auxcolinp(nlayer-i+1) = h2o2colx(i) + o3colx(i) + no2colx(i)
      end do

c     Interpolation

      do i=1,nz2
         auxi = nz2-i+1
         !H2O2 tabulated coefficient
         auxjh2o2(i) = jabsifotsintpar(auxi,6,indexint)
         !O3 tabulated coefficient
         auxjo3(i) = jabsifotsintpar(auxi,7,indexint)
         !Tabulated column
         auxcoltab(i) = c35(auxi)
      enddo
      !Only if chemthermod.ge.2
      if(chemthermod.ge.2) then
         do i=1,nz2
            !NO2 tabulated coefficient
            auxjno2(i) = jabsifotsintpar(nz2-i+1,13,indexint)
         enddo
      endif
      call interfast
     $     (wm,wp,auxind,auxcolinp,nlayer,auxcoltab,nz2,limdown,limup)
      do i=1,nlayer
         ind=auxind(i)
         auxi = nlayer-i+1
         !H2O2 interpolated coefficient
         jfotsout(indexint,6,auxi) = jfotsout(indexint,6,nlayer) *
     $        (wm(i)*auxjh2o2(ind+1) + wp(i)*auxjh2o2(ind))
         !O3 interpolated coefficient
         jfotsout(indexint,7,auxi) = jfotsout(indexint,7,nlayer) *
     $        (wm(i)*auxjo3(ind+1) + wp(i)*auxjo3(ind))
      enddo
      if(chemthermod.ge.2) then
         do i=1,nlayer
            ind=auxind(i)
            !NO2 interpolated coefficient
            jfotsout(indexint,13,nlayer-i+1) = 
     $           jfotsout(indexint,13,nlayer) *
     $           (wm(i)*auxjno2(ind+1) + wp(i)*auxjno2(ind))
         enddo
      endif

c     End of interval 35

ccccccccccccccccccccccccccccccc
c     337.8-800.0 nm (int 36)
c     
c     Absorption by:
c     O3, NO2
ccccccccccccccccccccccccccccccc

c     Input atmospheric column

      indexint=36
      do i=1,nlayer
         auxcolinp(nlayer-i+1) = o3colx(i) + no2colx(i)
      end do

c     Interpolation

      do i=1,nz2
         auxi = nz2-i+1
         !O3 tabulated coefficient
         auxjo3(i) = jabsifotsintpar(auxi,7,indexint)         
         !Tabulated column
         auxcoltab(i) = c36(auxi)
      enddo
      !Only if chemthermod.ge.2
      if(chemthermod.ge.2) then
         do i=1,nz2
            !NO2 tabulated coefficient
            auxjno2(i) = jabsifotsintpar(nz2-i+1,13,indexint)
         enddo
      endif
      call interfast
     $     (wm,wp,auxind,auxcolinp,nlayer,auxcoltab,nz2,limdown,limup)
      do i=1,nlayer
         ind=auxind(i)
         !O3 interpolated coefficient
         jfotsout(indexint,7,nlayer-i+1) = 
     $        jfotsout(indexint,7,nlayer) *
     $        (wm(i)*auxjo3(ind+1) + wp(i)*auxjo3(ind))
      enddo
      !Only if chemthermod.ge.2
      if(chemthermod.ge.2) then
         do i=1,nlayer
            ind=auxind(i)
            !NO2 interpolated coefficient
            jfotsout(indexint,13,nlayer-i+1) = 
     $           jfotsout(indexint,13,nlayer) *
     $           (wm(i)*auxjno2(ind+1) + wp(i)*auxjno2(ind))
         enddo
      endif

c     End of interval 36

c     End of interpolation to obtain photoabsorption rates


      return

      end






