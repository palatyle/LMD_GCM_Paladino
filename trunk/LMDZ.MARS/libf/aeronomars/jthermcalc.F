c**********************************************************************

      subroutine jthermcalc(ig,nlayer,chemthermod,
     .      rm,nesptherm,tx,iz,zenit)


c     feb 2002        fgg           first version
c     nov 2002        fgg           second version
c
c modified from paramhr.F
c MAC July 2003
c**********************************************************************

      use param_v4_h, only: jfotsout,crscabsi2,
     .    c1_16,c17_24,c25_29,c30_31,c32,c33,c34,c35,c36,
     .    co2crsc195,co2crsc295,t0,
     .    jabsifotsintpar,ninter,nz2

      implicit none

c     input and output variables

      integer    ig,nlayer
      integer    chemthermod
      integer    nesptherm                      !Number of species considered
      real       rm(nlayer,nesptherm)         !Densities (cm-3)
      real       tx(nlayer)                   !temperature
      real       zenit                          !SZA
      real       iz(nlayer)                   !Local altitude


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
      
      integer    i,j,k,indexint                 !indexes
      character  dn



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
         ! Ehouarn: test
          if ((ind+1.gt.nz2).or.(ind.le.0)) then
            write(*,*) "jthercalc error: ind=",ind,ig,zenit
            write(*,*) " auxind(1:nlayer)=",auxind
            write(*,*) " auxcolinp(:nlayer)=",auxcolinp
            write(*,*) " co2colx(:nlayer)=",co2colx
            write(*,*) " o2colx(:nlayer)=",o2colx
            write(*,*) " o3pcolx(:nlayer)=",o3pcolx
            write(*,*) " h2colx(:nlayer)=",h2colx
            write(*,*) " ncolx(:nlayer)=",ncolx
            write(*,*) " auxcoltab(1:nz2)=",auxcoltab
            write(*,*) " limdown=",limdown
            write(*,*) " limup=",limup
#ifndef MESOSCALE
            call abort_physic('jthermcalc','error',1)
#endif
          endif
         !CO2 interpolated coefficient
         jfotsout(indexint,1,auxi) = wm(i)*auxjco2(ind+1) +
     $        wp(i)*auxjco2(ind)
         !O2 interpolated coefficient
         jfotsout(indexint,2,auxi) = wm(i)*auxjo2(ind+1) +
     $        wp(i)*auxjo2(ind)
         !O3p interpolated coefficient
          jfotsout(indexint,3,auxi) = wm(i)*auxjo3p(ind+1) +
     $        wp(i)*auxjo3p(ind)
          !H2 interpolated coefficient
          jfotsout(indexint,5,auxi) = wm(i)*auxjh2(ind+1) +
     $         wp(i)*auxjh2(ind)
      enddo
      !Only if chemthermod.ge.2
      !N interpolated coefficient
      if(chemthermod.ge.2) then
         do i=1,nlayer
            ind=auxind(i)
            jfotsout(indexint,9,nlayer-i+1) =  wm(i)*auxjn(ind+1) +
     $         wp(i)*auxjn(ind)
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
             jfotsout(indexint,2,auxi) = wm(i)*auxjo2(ind+1) +
     $            wp(i)*auxjo2(ind)
             !O3p interpolated coefficient
             jfotsout(indexint,3,auxi) = wm(i)*auxjo3p(ind+1) +
     $            wp(i)*auxjo3p(ind)
             !CO2 interpolated coefficient
             jfotsout(indexint,1,auxi) = wm(i)*auxjco2(ind+1) +
     $            wp(i)*auxjco2(ind)
             !H2 interpolated coefficient
             jfotsout(indexint,5,auxi) = wm(i)*auxjh2(ind+1) +
     $            wp(i)*auxjh2(ind)
             !N2 interpolated coefficient
             jfotsout(indexint,8,auxi) = wm(i)*auxjn2(ind+1) +
     $            wp(i)*auxjn2(ind)             
             !CO interpolated coefficient
             jfotsout(indexint,11,auxi) = wm(i)*auxjco(ind+1) +
     $            wp(i)*auxjco(ind)
             !H interpolated coefficient
             jfotsout(indexint,12,auxi) = wm(i)*auxjh(ind+1) +
     $            wp(i)*auxjh(ind) 
          enddo
          !Only if chemthermod.ge.2
          if(chemthermod.ge.2) then
             do i=1,nlayer
                ind=auxind(i)
                auxi = nlayer-i+1
                !N interpolated coefficient
                jfotsout(indexint,9,auxi) = wm(i)*auxjn(ind+1) +
     $               wp(i)*auxjn(ind)
                !NO interpolated coefficient
                jfotsout(indexint,10,auxi)=wm(i)*auxjno(ind+1) +
     $               wp(i)*auxjno(ind)
                !NO2 interpolated coefficient
                jfotsout(indexint,13,auxi)=wm(i)*auxjno2(ind+1)+
     $               wp(i)*auxjno2(ind)
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
         jfotsout(indexint,2,auxi) = wm(i)*auxjo2(ind+1) +
     $            wp(i)*auxjo2(ind)
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
         !H interpolated coefficient
         jfotsout(indexint,12,auxi) = wm(i)*auxjh(ind+1) +
     $            wp(i)*auxjh(ind)         
      enddo
      !Only if chemthermod.ge.2
      if(chemthermod.ge.2) then
         do i=1,nlayer
            ind=auxind(i)
            auxi = nlayer-i+1
            !N interpolated coefficient
            jfotsout(indexint,9,auxi) = wm(i)*auxjn(ind+1) +
     $           wp(i)*auxjn(ind) 
            !NO interpolated coefficient
            jfotsout(indexint,10,auxi) = wm(i)*auxjno(ind+1) +
     $           wp(i)*auxjno(ind)
            !NO2 interpolated coefficient
            jfotsout(indexint,13,auxi) = wm(i)*auxjno2(ind+1) +
     $           wp(i)*auxjno2(ind)
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
            jfotsout(indexint,2,auxi) = (wm(i)*auxjo2(ind+1) +
     $           wp(i)*auxjo2(ind)) * cortemp(i)
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
         if(chemthermod.ge.2) then
            do i=1,nlayer
               ind=auxind(i)
               auxi = nlayer-i+1
               !NO interpolated coefficient
               jfotsout(indexint,10,auxi)=(wm(i)*auxjno(ind+1) +
     $              wp(i)*auxjno(ind)) * cortemp(i)
               !NO2 interpolated coefficient
               jfotsout(indexint,13,auxi)=(wm(i)*auxjno2(ind+1)+
     $              wp(i)*auxjno2(ind)) * cortemp(i)
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
            jfotsout(indexint,1,auxi) = (wm(i)*auxjco2(ind+1) +
     $           wp(i)*auxjco2(ind)) * cortemp(i) *
     $           (1+alfa(indexint,auxi)*
     $           (t2(auxi)-t0(auxi)))
            !O2 interpolated coefficient
            jfotsout(indexint,2,auxi) = (wm(i)*auxjo2(ind+1) +
     $           wp(i)*auxjo2(ind)) * cortemp(i)
            !H2O interpolated coefficient
            jfotsout(indexint,4,auxi) = (wm(i)*auxjh2o(ind+1) +
     $           wp(i)*auxjh2o(ind)) * cortemp(i)
            !H2O2 interpolated coefficient
            jfotsout(indexint,6,auxi) = (wm(i)*auxjh2o2(ind+1) +
     $           wp(i)*auxjh2o2(ind)) * cortemp(i)            
            !CO interpolated coefficient
            jfotsout(indexint,11,auxi) = (wm(i)*auxjco(ind+1) +
     $           wp(i)*auxjco(ind)) * cortemp(i)
         enddo
         !Only if chemthermod.ge.2
         if(chemthermod.ge.2) then
            do i=1,nlayer
               ind=auxind(i)
               auxi = nlayer-i+1
               !NO interpolated coefficient
               jfotsout(indexint,10,auxi)=(wm(i)*auxjno(ind+1) +
     $              wp(i)*auxjno(ind)) * cortemp(i)
               !NO2 interpolated coefficient
               jfotsout(indexint,13,auxi)=(wm(i)*auxjno2(ind+1)+
     $              wp(i)*auxjno2(ind)) * cortemp(i)
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
            jfotsout(indexint,1,auxi) = (wm(i)*auxjco2(ind+1) +
     $           wp(i)*auxjco2(ind)) * cortemp(i) *
     $           (1+alfa(indexint,auxi)*
     $           (t2(auxi)-t0(auxi)))
            !O2 interpolated coefficient
            jfotsout(indexint,2,auxi) = (wm(i)*auxjo2(ind+1) +
     $            wp(i)*auxjo2(ind)) * cortemp(i)
            !H2O interpolated coefficient
            jfotsout(indexint,4,auxi) = (wm(i)*auxjh2o(ind+1) +
     $            wp(i)*auxjh2o(ind)) * cortemp(i)
            !H2O2 interpolated coefficient
            jfotsout(indexint,6,auxi) = (wm(i)*auxjh2o2(ind+1) +
     $            wp(i)*auxjh2o2(ind)) * cortemp(i)            
         enddo
         !Only if chemthermod.ge.2
         if(chemthermod.ge.2) then
            do i=1,nlayer 
               ind=auxind(i)
               auxi = nlayer-i+1
               !NO interpolated coefficient
               jfotsout(indexint,10,auxi)=(wm(i)*auxjno(ind+1) +
     $              wp(i)*auxjno(ind)) * cortemp(i)
               !NO2 interpolated coefficient
               jfotsout(indexint,13,auxi)=(wm(i)*auxjno2(ind+1)+
     $              wp(i)*auxjno2(ind)) * cortemp(i)
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
         jfotsout(indexint,1,auxi) = (wm(i)*auxjco2(ind+1) +
     $        wp(i)*auxjco2(ind)) * cortemp(i) *
     $        (1+alfa(indexint,auxi)*
     $        (t2(auxi)-t0(auxi)))
         !O2 interpolated coefficient
         jfotsout(indexint,2,auxi) = (wm(i)*auxjo2(ind+1) +
     $        wp(i)*auxjo2(ind)) * cortemp(i)
         !H2O2 interpolated coefficient
         jfotsout(indexint,6,auxi) = (wm(i)*auxjh2o2(ind+1) +
     $        wp(i)*auxjh2o2(ind)) * cortemp(i)         
      enddo
      !Only if chemthermod.ge.2
      if(chemthermod.ge.2) then
         do i=1,nlayer
            auxi = nlayer-i+1
            ind=auxind(i)
            !NO interpolated coefficient
            jfotsout(indexint,10,auxi) = (wm(i)*auxjno(ind+1) +
     $           wp(i)*auxjno(ind)) * cortemp(i)
           !NO2 interpolated coefficient
            jfotsout(indexint,13,auxi) = (wm(i)*auxjno2(ind+1) +
     $           wp(i)*auxjno2(ind)) * cortemp(i)
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
         jfotsout(indexint,2,auxi) = wm(i)*auxjo2(ind+1) +
     $        wp(i)*auxjo2(ind)
         !H2O2 interpolated coefficient
         jfotsout(indexint,6,auxi) = wm(i)*auxjh2o2(ind+1) +
     $        wp(i)*auxjh2o2(ind)         
      enddo
      !Only if chemthermod.ge.2
      if(chemthermod.ge.2) then
         do i=1,nlayer
            ind=auxind(i)
            !NO2 interpolated coefficient
            jfotsout(indexint,13,nlayer-i+1) = wm(i)*auxjno2(ind+1) +
     $           wp(i)*auxjno2(ind)
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
         jfotsout(indexint,2,auxi) = wm(i)*auxjo2(ind+1) +
     $        wp(i)*auxjo2(ind)
         !H2O2 interpolated coefficient
         jfotsout(indexint,6,auxi) = wm(i)*auxjh2o2(ind+1) +
     $        wp(i)*auxjh2o2(ind)
         !O3 interpolated coefficient
         jfotsout(indexint,7,auxi) = wm(i)*auxjo3(ind+1) +
     $        wp(i)*auxjo3(ind)         
      enddo
      !Only if chemthermod.ge.2
      if(chemthermod.ge.2) then
         do i=1,nlayer
            ind=auxind(i)
            !NO2 interpolated coefficient
            jfotsout(indexint,13,nlayer-i+1) = wm(i)*auxjno2(ind+1) +
     $           wp(i)*auxjno2(ind)
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
         jfotsout(indexint,6,auxi) = wm(i)*auxjh2o2(ind+1) +
     $        wp(i)*auxjh2o2(ind)
         !O3 interpolated coefficient
         jfotsout(indexint,7,auxi) = wm(i)*auxjo3(ind+1) +
     $        wp(i)*auxjo3(ind)         
      enddo
      if(chemthermod.ge.2) then
         do i=1,nlayer
            ind=auxind(i)
            !NO2 interpolated coefficient
            jfotsout(indexint,13,nlayer-i+1) = wm(i)*auxjno2(ind+1) +
     $           wp(i)*auxjno2(ind)
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
         jfotsout(indexint,7,nlayer-i+1) = wm(i)*auxjo3(ind+1) +
     $        wp(i)*auxjo3(ind)         
      enddo
      !Only if chemthermod.ge.2
      if(chemthermod.ge.2) then
         do i=1,nlayer
            ind=auxind(i)
            !NO2 interpolated coefficient
            jfotsout(indexint,13,nlayer-i+1) = wm(i)*auxjno2(ind+1) +
     $           wp(i)*auxjno2(ind)
         enddo
      endif

c     End of interval 36

c     End of interpolation to obtain photoabsorption rates


      return

      end



c**********************************************************************
c**********************************************************************

      subroutine column(ig,nlayer,chemthermod,rm,nesptherm,tx,iz,zenit,
     $     co2colx,o2colx,o3pcolx,h2colx,h2ocolx,h2o2colx,o3colx,
     $     n2colx,ncolx,nocolx,cocolx,hcolx,no2colx)

c     nov 2002        fgg           first version

c**********************************************************************

      use tracer_mod, only: igcm_o, igcm_co2, igcm_o2, igcm_h2,
     &                      igcm_h2o_vap, igcm_h2o2, igcm_co, igcm_h,
     &                      igcm_o3, igcm_n2, igcm_n, igcm_no, igcm_no2,
     &                      mmol
      use param_v4_h, only: radio,gg,masa,kboltzman,n_avog

      implicit none


c     common variables and constants
      include 'callkeys.h'



c    local parameters and variables



c     input and output variables

      integer    ig,nlayer
      integer    chemthermod
      integer    nesptherm                      !# of species undergoing chemistry, input
      real       rm(nlayer,nesptherm)         !densities (cm-3), input
      real       tx(nlayer)                   !temperature profile, input
      real       iz(nlayer+1)                 !height profile, input
      real       zenit                          !SZA, input
      real       co2colx(nlayer)              !column density of CO2 (cm^-2), output
      real       o2colx(nlayer)               !column density of O2(cm^-2), output
      real       o3pcolx(nlayer)              !column density of O(3P)(cm^-2), output
      real       h2colx(nlayer)               !H2 column density (cm-2), output
      real       h2ocolx(nlayer)              !H2O column density (cm-2), output
      real       h2o2colx(nlayer)             !column density of H2O2(cm^-2), output
      real       o3colx(nlayer)               !O3 column density (cm-2), output
      real       n2colx(nlayer)               !N2 column density (cm-2), output
      real       ncolx(nlayer)                !N column density (cm-2), output
      real       nocolx(nlayer)               !NO column density (cm-2), output
      real       cocolx(nlayer)               !CO column density (cm-2), output
      real       hcolx(nlayer)                !H column density (cm-2), output
      real       no2colx(nlayer)              !NO2 column density (cm-2), output


c     local variables

      real       xx
      real       grav(nlayer)
      real       Hco2,Ho3p,Ho2,Hh2,Hh2o,Hh2o2
      real       Ho3,Hn2,Hn,Hno,Hco,Hh,Hno2

      real       co2x(nlayer)
      real       o2x(nlayer)
      real       o3px(nlayer)
      real       cox(nlayer)
      real       hx(nlayer)
      real       h2x(nlayer)
      real       h2ox(nlayer)
      real       h2o2x(nlayer)
      real       o3x(nlayer)
      real       n2x(nlayer)
      real       nx(nlayer)
      real       nox(nlayer)
      real       no2x(nlayer)

      integer    i,j,k,icol,indexint          !indexes

c     variables for optical path calculation

      integer    nz3
!      parameter  (nz3=nz*2)

      integer    jj
      real*8      esp(nlayer*2)
      real*8      ilayesp(nlayer*2)
      real*8      szalayesp(nlayer*2)
      integer     nlayesp
      real*8      zmini
      real*8      depth
      real*8      espco2, espo2, espo3p, esph2, esph2o, esph2o2,espo3
      real*8      espn2,espn,espno,espco,esph,espno2
      real*8      rcmnz, rcmmini
      real*8      szadeg


      ! Tracer indexes in the thermospheric chemistry:
      !!! ATTENTION. These values have to be identical to those in chemthermos.F90
      !!! If the values are changed there, the same has to be done here  !!!
      integer,parameter :: i_co2  =  1
      integer,parameter :: i_co   =  2
      integer,parameter :: i_o    =  3
      integer,parameter :: i_o1d  =  4
      integer,parameter :: i_o2   =  5
      integer,parameter :: i_o3   =  6
      integer,parameter :: i_h    =  7
      integer,parameter :: i_h2   =  8
      integer,parameter :: i_oh   =  9
      integer,parameter :: i_ho2  = 10
      integer,parameter :: i_h2o2 = 11
      integer,parameter :: i_h2o  = 12
      integer,parameter :: i_n    = 13
      integer,parameter :: i_n2d  = 14
      integer,parameter :: i_no   = 15
      integer,parameter :: i_no2  = 16
      integer,parameter :: i_n2   = 17
!      integer,parameter :: i_co2=1
!      integer,parameter :: i_o2=2
!      integer,parameter :: i_o=3
!      integer,parameter :: i_co=4
!      integer,parameter :: i_h=5
!      integer,parameter :: i_h2=8
!      integer,parameter :: i_h2o=9
!      integer,parameter :: i_h2o2=10
!      integer,parameter :: i_o3=12
!      integer,parameter :: i_n2=13
!      integer,parameter :: i_n=14
!      integer,parameter :: i_no=15
!      integer,parameter :: i_no2=17


c*************************PROGRAM STARTS*******************************

      nz3 = nlayer*2
      do i=1,nlayer
         xx = ( radio + iz(i) ) * 1.e5
         grav(i) = gg * masa /(xx**2)
      end do

      !Scale heights
      xx = kboltzman * tx(nlayer) * n_avog / grav(nlayer)
      Ho3p  = xx / mmol(igcm_o) 
      Hco2  = xx / mmol(igcm_co2)
      Ho2   = xx / mmol(igcm_o2)
      Hh2   = xx / mmol(igcm_h2)
      Hh2o  = xx / mmol(igcm_h2o_vap)
      Hh2o2 = xx / mmol(igcm_h2o2)
      Hco   = xx / mmol(igcm_co)
      Hh    = xx / mmol(igcm_h)
      !Only if O3 chem. required
      if(chemthermod.ge.1) 
     $     Ho3   = xx / mmol(igcm_o3)
      !Only if N or ion chem.
      if(chemthermod.ge.2) then
         Hn2   = xx / mmol(igcm_n2)
         Hn    = xx / mmol(igcm_n)
         Hno   = xx / mmol(igcm_no)
         Hno2  = xx / mmol(igcm_no2)
      endif
      ! first loop in altitude : initialisation
      do i=nlayer,1,-1
         !Column initialisation
         co2colx(i)  = 0.
         o2colx(i)   = 0.
         o3pcolx(i)  = 0.
         h2colx(i)   = 0.
         h2ocolx(i)  = 0.
         h2o2colx(i) = 0.
         o3colx(i)   = 0.
         n2colx(i)   = 0.
         ncolx(i)    = 0.
         nocolx(i)   = 0.
         cocolx(i)   = 0.
         hcolx(i)    = 0.
         no2colx(i)  = 0.
         !Densities
         co2x(i)  = rm(i,i_co2)
         o2x(i)   = rm(i,i_o2)
         o3px(i)  = rm(i,i_o)
         h2x(i)   = rm(i,i_h2)
         h2ox(i)  = rm(i,i_h2o)
         h2o2x(i) = rm(i,i_h2o2)
         cox(i)   = rm(i,i_co)
         hx(i)    = rm(i,i_h)
         !Only if O3 chem. required
         if(chemthermod.ge.1) 
     $        o3x(i)   = rm(i,i_o3)
         !Only if Nitrogen of ion chemistry requested
         if(chemthermod.ge.2) then
            n2x(i)   = rm(i,i_n2)
            nx(i)    = rm(i,i_n)
            nox(i)   = rm(i,i_no)
            no2x(i)  = rm(i,i_no2)
         endif
      enddo
      ! second loop in altitude : column calculations
      do i=nlayer,1,-1
         !Routine to calculate the geometrical length of each layer
         call espesor_optico_A(ig,i,nlayer,zenit,iz(i),nz3,iz,esp,
     $         ilayesp,szalayesp,nlayesp, zmini)
         if(ilayesp(nlayesp).eq.-1) then
            co2colx(i)=1.e25
            o2colx(i)=1.e25
            o3pcolx(i)=1.e25
            h2colx(i)=1.e25
            h2ocolx(i)=1.e25
            h2o2colx(i)=1.e25
            o3colx(i)=1.e25
            n2colx(i)=1.e25
            ncolx(i)=1.e25
            nocolx(i)=1.e25
            cocolx(i)=1.e25
            hcolx(i)=1.e25
            no2colx(i)=1.e25
         else
            rcmnz = ( radio + iz(nlayer) ) * 1.e5
            rcmmini = ( radio + zmini ) * 1.e5
            !Column calculation taking into account the geometrical depth 
            !calculated before
            do j=1,nlayesp
               jj=ilayesp(j)
               !Top layer
               if(jj.eq.nlayer) then
                  if(zenit.le.60.) then 
                     o3pcolx(i)=o3pcolx(i)+o3px(nlayer)*Ho3p*esp(j)
     $                    *1.e-5
                     co2colx(i)=co2colx(i)+co2x(nlayer)*Hco2*esp(j)
     $                    *1.e-5
                     h2o2colx(i)=h2o2colx(i)+
     $                    h2o2x(nlayer)*Hh2o2*esp(j)*1.e-5
                     o2colx(i)=o2colx(i)+o2x(nlayer)*Ho2*esp(j)
     $                    *1.e-5
                     h2colx(i)=h2colx(i)+h2x(nlayer)*Hh2*esp(j)
     $                    *1.e-5
                     h2ocolx(i)=h2ocolx(i)+h2ox(nlayer)*Hh2o*esp(j)
     $                    *1.e-5                     
                     cocolx(i)=cocolx(i)+cox(nlayer)*Hco*esp(j)
     $                    *1.e-5
                     hcolx(i)=hcolx(i)+hx(nlayer)*Hh*esp(j)
     $                    *1.e-5
                     !Only if O3 chemistry required
                     if(chemthermod.ge.1) o3colx(i)=
     $                    o3colx(i)+o3x(nlayer)*Ho3*esp(j)
     $                    *1.e-5
                     !Only if N or ion chemistry requested
                     if(chemthermod.ge.2) then
                        n2colx(i)=n2colx(i)+n2x(nlayer)*Hn2*esp(j)
     $                    *1.e-5
                        ncolx(i)=ncolx(i)+nx(nlayer)*Hn*esp(j)
     $                       *1.e-5
                        nocolx(i)=nocolx(i)+nox(nlayer)*Hno*esp(j)
     $                       *1.e-5
                        no2colx(i)=no2colx(i)+no2x(nlayer)*Hno2*esp(j)
     $                       *1.e-5
                     endif
                  else if(zenit.gt.60.) then
                     espco2 =sqrt((rcmnz+Hco2)**2 -rcmmini**2) - esp(j)
                     espo2  = sqrt((rcmnz+Ho2)**2 -rcmmini**2) - esp(j)
                     espo3p = sqrt((rcmnz+Ho3p)**2 -rcmmini**2)- esp(j)
                     esph2  = sqrt((rcmnz+Hh2)**2 -rcmmini**2) - esp(j)
                     esph2o = sqrt((rcmnz+Hh2o)**2 -rcmmini**2)- esp(j)
                     esph2o2= sqrt((rcmnz+Hh2o2)**2-rcmmini**2)- esp(j)
                     espco  = sqrt((rcmnz+Hco)**2 -rcmmini**2) - esp(j)
                     esph   = sqrt((rcmnz+Hh)**2 -rcmmini**2)  - esp(j)
                     !Only if O3 chemistry required
                     if(chemthermod.ge.1)                     
     $                   espo3=sqrt((rcmnz+Ho3)**2-rcmmini**2)-esp(j)
                     !Only if N or ion chemistry requested
                     if(chemthermod.ge.2) then
                        espn2 =sqrt((rcmnz+Hn2)**2-rcmmini**2)-esp(j)
                        espn  =sqrt((rcmnz+Hn)**2-rcmmini**2)  - esp(j)
                        espno =sqrt((rcmnz+Hno)**2-rcmmini**2) - esp(j)
                        espno2=sqrt((rcmnz+Hno2)**2-rcmmini**2)- esp(j)
                     endif
                     co2colx(i) = co2colx(i) + espco2*co2x(nlayer)
                     o2colx(i)  = o2colx(i)  + espo2*o2x(nlayer)
                     o3pcolx(i) = o3pcolx(i) + espo3p*o3px(nlayer)
                     h2colx(i)  = h2colx(i)  + esph2*h2x(nlayer)
                     h2ocolx(i) = h2ocolx(i) + esph2o*h2ox(nlayer)
                     h2o2colx(i)= h2o2colx(i)+ esph2o2*h2o2x(nlayer)
                     cocolx(i)  = cocolx(i)  + espco*cox(nlayer)
                     hcolx(i)   = hcolx(i)   + esph*hx(nlayer)
                     !Only if O3 chemistry required
                     if(chemthermod.ge.1)                      
     $                  o3colx(i) = o3colx(i)  + espo3*o3x(nlayer)
                     !Only if N or ion chemistry requested
                     if(chemthermod.ge.2) then
                        n2colx(i)  = n2colx(i)  + espn2*n2x(nlayer)
                        ncolx(i)   = ncolx(i)   + espn*nx(nlayer)
                        nocolx(i)  = nocolx(i)  + espno*nox(nlayer)
                        no2colx(i) = no2colx(i) + espno2*no2x(nlayer)
                     endif
                  endif !Of if zenit.lt.60
               !Other layers
               else 
                  co2colx(i)  = co2colx(i) + 
     $                 esp(j) * (co2x(jj)+co2x(jj+1)) / 2.
                  o2colx(i)   = o2colx(i) + 
     $                 esp(j) * (o2x(jj)+o2x(jj+1)) / 2.
                  o3pcolx(i)  = o3pcolx(i) + 
     $                 esp(j) * (o3px(jj)+o3px(jj+1)) / 2.
                  h2colx(i)   = h2colx(i) + 
     $                 esp(j) * (h2x(jj)+h2x(jj+1)) / 2.
                  h2ocolx(i)  = h2ocolx(i) + 
     $                 esp(j) * (h2ox(jj)+h2ox(jj+1)) / 2.
                  h2o2colx(i) = h2o2colx(i) + 
     $                 esp(j) * (h2o2x(jj)+h2o2x(jj+1)) / 2.
                  cocolx(i)   = cocolx(i) + 
     $                 esp(j) * (cox(jj)+cox(jj+1)) / 2.
                  hcolx(i)    = hcolx(i) + 
     $                 esp(j) * (hx(jj)+hx(jj+1)) / 2.
                  !Only if O3 chemistry required
                  if(chemthermod.ge.1) 
     $                 o3colx(i) = o3colx(i) + 
     $                 esp(j) * (o3x(jj)+o3x(jj+1)) / 2.
                  !Only if N or ion chemistry requested
                  if(chemthermod.ge.2) then
                     n2colx(i)   = n2colx(i) + 
     $                 esp(j) * (n2x(jj)+n2x(jj+1)) / 2.
                     ncolx(i)    = ncolx(i) + 
     $                    esp(j) * (nx(jj)+nx(jj+1)) / 2.
                     nocolx(i)   = nocolx(i) + 
     $                    esp(j) * (nox(jj)+nox(jj+1)) / 2.
                     no2colx(i)  = no2colx(i) + 
     $                    esp(j) * (no2x(jj)+no2x(jj+1)) / 2.
                  endif
               endif  !Of if jj.eq.nlayer
            end do    !Of do j=1,nlayesp
         endif        !Of ilayesp(nlayesp).eq.-1
      enddo           !Of do i=nlayer,1,-1

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
      real*8,intent(out) :: wm(nlayer),wp(nlayer) ! interpolation weights
      integer,intent(out) :: nm(nlayer) ! index of nearest point
      real*8,intent(in) :: pin(nl),p(nlayer)
      real*8,intent(in) :: limup,limdown
      integer,intent(in) :: nl,nlayer
      integer :: n1,n,np,nini
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

      end


c**********************************************************************
c**********************************************************************

      subroutine espesor_optico_A (ig,capa,nlayer, szadeg,z,
     @                   nz3,iz,esp,ilayesp,szalayesp,nlayesp, zmini)

c     fgg              nov 03      Adaptation to Martian model
c     malv             jul 03      Corrected z grid. Split in alt & frec codes
c     fgg              feb 03      first version
*************************************************************************

      use param_v4_h, only: radio
      implicit none

c     arguments

      real        szadeg                ! I. SZA [rad]
      real        z                     ! I. altitude of interest [km]
      integer     nz3,ig,nlayer         ! I. dimension of esp, ylayesp, etc...
                                        !  (=2*nlayer= max# of layers in ray path)
      real     iz(nlayer+1)              ! I. Altitude of each layer
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
        real*8      diz(nlayer+1)
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
        do i=1,nlayer
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
           if(abs(zz-diz(nlayer)).lt.1.d-3) goto 1357
           ! 1st Zone: Upper part of ray
           !
           do j=grid_R8(zz,diz,nlayer),nlayer-1
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
           ilayesp(nlayesp) = nlayer
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

           if(abs(zz-diz(nlayer)).lt.1.d-4) goto 1470

           ! 1st Zone: Upper part of ray
           !
           do j=grid_R8(zz,diz,nlayer),nlayer-1
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
           ilayesp(nlayesp) = nlayer
           rkmnz = radio+diz(nlayer)
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
           !zmini should be lower than zz, as SZA<90. However, in situations
           !where SZA is very close to 90, rounding errors can make zmini
           !slightly higher than zz, causing problems in the determination
           !of the jmin index. A correction is implemented in the determination
           !of jmin, some lines below
           rkmmini = radio + zmini

           if(zmini.lt.diz(1)) then         ! Can see the sun?  No => esp(j)=inft
             nlayesp = nlayesp + 1 
             ilayesp(nlayesp) = - 1     ! Value to mark "no sun on view"
!             esp(nlayesp) = 1.e30

           else
              jmin=grid_R8(zmini,diz,nlayer)+1
              !Correction for possible rounding errors when SZA very close 
              !to 90 degrees
              if(jmin.gt.grid_R8(zz,diz,nlayer)) then
                 write(*,*)'jthermcalc warning: possible rounding error'
                 write(*,*)'point,sza,layer:',ig,szadeg,capa
                 jmin=grid_R8(zz,diz,nlayer)
              endif

              if(abs(zz-diz(nlayer)).lt.1.d-4) goto 9876

              ! 1st Zone: Upper part of ray
              !
              do j=grid_R8(zz,diz,nlayer),nlayer-1
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
              ilayesp(nlayesp) = nlayer
              rkmnz = radio+diz(nlayer)
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
              do j=jmin,grid_R8(zz,diz,nlayer)-1
                nlayesp = nlayesp + 1 
                ilayesp(nlayesp) = j
                esp(nlayesp) = 
     $                    sqrt( (radio+diz(j+1))**2 - rkmmini**2 )
     $                  - sqrt( (radio+diz(j))**2 - rkmmini**2 )       ! [km]
                esp(nlayesp) = esp(nlayesp) * 1.d5                     ! [cm]
                rkmj = radio+diz(j)
                szalayesp(nlayesp) = asin( rkmmini/rkmj )              ! [rad]
                szalayesp(nlayesp) = szalayesp(nlayesp) *180.d0/3.141592 ! [deg]
              end do

              ! 5th zone: Lower part of ray, decreasing from z to zmin
              !    ( layers with SZA > 90 deg )
              do j=grid_R8(zz,diz,nlayer)-1, jmin, -1
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

        if ( z .lt. zgrid(1)  ) then 
           write (*,*) ' GRID/ z outside bounds of zgrid '
           write (*,*) ' z,zgrid(1),zgrid(nz) =', z,zgrid(1),zgrid(nz)
           z = zgrid(1) 
           write(*,*) 'WARNING: error in grid_r8 (jthermcalc.F)'
           write(*,*) 'Please check values of z and zgrid above'
        endif
        if (z .gt. zgrid(nz) ) then
           write (*,*) ' GRID/ z outside bounds of zgrid '
           write (*,*) ' z,zgrid(1),zgrid(nz) =', z,zgrid(1),zgrid(nz)
           z = zgrid(nz) 
           write(*,*) 'WARNING: error in grid_r8 (jthermcalc.F)'
           write(*,*) 'Please check values of z and zgrid above'
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

      use comsaison_h, only: dist_sol
      use param_v4_h, only: ninter, 
     .                      fluxtop, ct1, ct2, p1, p2
      implicit none


!     common variables and constants
      include "callkeys.h"


!     Arguments

      real date


!     Local variable and constants

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
        if(i.le.24.and.solvarmod.eq.0) then
            fluxtop(i)=(((ct1(i)+p1(i)*date)/2.)                  
     $           *sin(2.*3.1416/11.*(date-1985.-3.1416))          
     $           +(ct2(i)+p2(i)*date)+1.)*fluxtop(i)
         end if
         fluxtop(i)=fluxtop(i)*(1.52/dist_sol)**2
      end do
     
      return
      end
