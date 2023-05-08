      SUBROUTINE chemthermos(ig,nlayer,lswitch,chemthermod,zycol,ztemp, &
           zdens,zpress,zlocal,zenit,ptimestep,zday,em_no,em_o2)

      use tracer_mod, only: nqmx, igcm_co2, igcm_co, igcm_o, igcm_o1d,  &
                            igcm_o2, igcm_h, igcm_h2, igcm_oh, igcm_ho2,&
                            igcm_h2o2, igcm_h2o_vap, igcm_o3, igcm_n2,  &
                            igcm_n, igcm_no, igcm_no2, igcm_n2d,        &
                            igcm_co2plus, igcm_o2plus, igcm_coplus,     &
                            igcm_cplus, igcm_nplus, igcm_noplus,        &
                            igcm_n2plus, igcm_hplus, igcm_hco2plus,     &
                            igcm_elec, igcm_oplus

      use param_v4_h, only: Pno, Po2
      USE comcstfi_h
      IMPLICIT NONE
!=======================================================================
!   subject: 
!   --------
!   Computing chemical variations in the thermosphere
!
!   author:  MAC July 2003
!   ------
!   modifications:
!   -------------
!     Ehouarn Sept 2008: added handling of tracers by their names
!     Francisco Sept 2009: added ionosphere and N chemistry
!     Francisco Mar 2012: Different chemical schemes possible according to 
!                         traceur.def
!=======================================================================
!
!    0.  Declarations :
!    ------------------
!
#include "callkeys.h"
!-----------------------------------------------------------------------
!    Input/Output
!    ------------
      integer :: lswitch,ig,chemthermod,nlayer
      real :: zycol(nlayer,nqmx)
      real :: ztemp(nlayer)
      real :: zdens(nlayer)
      real :: zpress(nlayer)			! in mbar
      real :: zlocal(nlayer)
      real :: zenit
      real :: ptimestep
      real :: zday 
      real :: em_no(nlayer)
      real :: em_o2(nlayer)
!
!    Local variables :
!    -----------------
      integer :: l
      integer,save :: nesptherm
      real,allocatable :: rm(:,:) 		!number density (cm-3)
      logical,save :: firstcall=.true.

! Tracer indexes in the GCM:
      integer,save :: g_co2=0
      integer,save :: g_co=0
      integer,save :: g_o=0
      integer,save :: g_o1d=0
      integer,save :: g_o2=0
      integer,save :: g_o3=0
      integer,save :: g_h=0
      integer,save :: g_h2=0
      integer,save :: g_oh=0
      integer,save :: g_ho2=0
      integer,save :: g_h2o2=0
      integer,save :: g_n2=0
      integer,save :: g_h2o_vap=0
      integer,save :: g_n=0
      integer,save :: g_no=0
      integer,save :: g_no2=0
      integer,save :: g_n2d=0
      integer,save :: g_co2plus=0
      integer,save :: g_oplus=0
      integer,save :: g_o2plus=0
      integer,save :: g_coplus=0
      integer,save :: g_cplus=0
      integer,save :: g_nplus=0
      integer,save :: g_noplus=0
      integer,save :: g_n2plus=0
      integer,save :: g_hplus=0
      integer,save :: g_hco2plus=0
      integer,save :: g_elec=0
      
! Tracer indexes in the thermospheric chemistry:
!!! IMPORTANT. These values have to be the same than in
!!! jthermcal.F, column.F, paramfoto_compact.F, euvheat.F
!!! and hrtherm.F
!!! If some value is changed here, the same change has to 
!!! be made into the other subroutines
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
!      integer,parameter :: i_oh=6
!      integer,parameter :: i_ho2=7
!      integer,parameter :: i_h2=8
!      integer,parameter :: i_h2o=9
!      integer,parameter :: i_h2o2=10
!      integer,parameter :: i_o1d=11
!      integer,parameter :: i_o3=12
!      integer,parameter :: i_n2=13
!      integer,parameter :: i_n=14
!      integer,parameter :: i_no=15
!      integer,parameter :: i_n2d=16
!      integer,parameter :: i_no2=17
      integer,parameter :: i_co2plus=18
      integer,parameter :: i_oplus=19
      integer,parameter :: i_o2plus=20
      integer,parameter :: i_coplus=21
      integer,parameter :: i_cplus=22
      integer,parameter :: i_nplus=23
      integer,parameter :: i_noplus=24
      integer,parameter :: i_n2plus=25
      integer,parameter :: i_hplus=26
      integer,parameter :: i_hco2plus=27
      integer,parameter :: i_elec=28

      
! Initializations at first call
      if (firstcall) then
         !Number of species to consider
         nesptherm=0
         ! get the indexes of the tracers we'll need
         ! It is not really necessary to check if the species are present,
         ! as it is already done in calchim.F90. But we include it to double-check

         ! Species always included if photochemistry
         g_co2=igcm_co2
         if (g_co2.eq.0) then
            write(*,*) "chemthermos: Error; no CO2 tracer !!!"
            write(*,*) "CO2 must always be present if thermochem=T"
            stop
         else
            nesptherm=nesptherm+1  
         endif
         g_co=igcm_co
         if (g_co.eq.0) then
            write(*,*) "chemthermos: Error; no CO tracer !!!"
            write(*,*) "CO must always be present if thermochem=T"
            stop
         else
            nesptherm=nesptherm+1  
         endif
         g_o=igcm_o
         if (g_o.eq.0) then
            write(*,*) "chemthermos: Error; no O tracer !!!"
            write(*,*) "O must always be present if thermochem=T"
            stop
         else
            nesptherm=nesptherm+1  
         endif
         g_o1d=igcm_o1d
         if (g_o1d.eq.0) then
            write(*,*) "chemthermos: Error; no O1D tracer !!!"
            write(*,*) "O1D must always be present if thermochem=T"
            stop
         else
            nesptherm=nesptherm+1  
         endif
         g_o2=igcm_o2
         if (g_o2.eq.0) then
            write(*,*) "chemthermos: Error; no O2 tracer !!!"
            write(*,*) "O2 must always be present if thermochem=T"
            stop
         else
            nesptherm=nesptherm+1  
         endif
         g_h=igcm_h
         if (g_h.eq.0) then
            write(*,*) "chemthermos: Error; no H tracer !!!"
            write(*,*) "H must always be present if thermochem=T"
            stop
         else
            nesptherm=nesptherm+1  
         endif
         g_h2=igcm_h2
         if (g_h2.eq.0) then
            write(*,*) "chemthermos: Error; no H2 tracer !!!"
            write(*,*) "H2 must always be present if thermochem=T"
            stop
         else
            nesptherm=nesptherm+1  
         endif
         g_oh=igcm_oh
         if (g_oh.eq.0) then
            write(*,*) "chemthermos: Error; no OH tracer !!!"
            write(*,*) "OH must always be present if thermochem=T"
            stop
         else
            nesptherm=nesptherm+1  
         endif
         g_ho2=igcm_ho2
         if (g_ho2.eq.0) then
            write(*,*) "chemthermos: Error; no HO2 tracer !!!"
            write(*,*) "HO2 must always be present if thermochem=T"
            stop
         else
            nesptherm=nesptherm+1  
         endif
         g_h2o2=igcm_h2o2
         if (g_h2o2.eq.0) then
            write(*,*) "chemthermos: Error; no H2O2 tracer !!!"
            write(*,*) "H2O2 must always be present if thermochem=T"
            stop
         else
            nesptherm=nesptherm+1  
         endif
         g_h2o_vap=igcm_h2o_vap
         if (g_h2o_vap.eq.0) then
            write(*,*) "chemthermos: Error; no water vapor tracer !!!"
            write(*,*) "H2O must always be present if thermochem=T"
            stop
         else
            nesptherm=nesptherm+1  
         endif
         !Verify if O3 chemistry wanted
         g_o3=igcm_o3
         if(chemthermod.ge.1) then
            if (g_o3.eq.0) then
               write(*,*) "chemthermos: Error; no O3 tracer !!!"
               write(*,*) "O3 must be present if NO is in traceur.def"
               stop
            else
               nesptherm=nesptherm+1
            endif
         else
            if(g_o3 == 0) then
               write(*,*) "chemthermos: No O3 chemistry"
            else               
               write(*,*) "chemthermos: O3 chemistry included"
               nesptherm=nesptherm+1  
               chemthermod=1
            endif
         endif    ! Of if(chemthermod.ge.1)
            
         ! N chemistry
         if(chemthermod.ge.2) then
            g_n2=igcm_n2
            if (g_n2.eq.0) then
               write(*,*) "chemthermos: Error; no N2 tracer !!!"
               write(*,*) "N2 must be present if NO is in traceur.def"
               stop
            else
               nesptherm=nesptherm+1  
            endif
            g_n=igcm_n
            if (g_n.eq.0) then
               write(*,*) "chemthermos: Error; no N tracer !!!"
               write(*,*) "N must be present if NO is in traceur.def"
               stop
            else
               nesptherm=nesptherm+1  
            endif
            g_no=igcm_no
            if (g_no.eq.0) then
               write(*,*) "chemthermos: Error; no NO tracer !!!"
               write(*,*) "NO must be present if N chemistry wanted"
               stop
            else
               nesptherm=nesptherm+1  
            endif
            g_no2=igcm_no2
            if (g_no2.eq.0) then
               write(*,*) "chemthermos: Error; no NO2 tracer !!!"
               write(*,*) "NO must be present if NO is in traceur.def"
               stop
            else
               nesptherm=nesptherm+1  
            endif
            g_n2d=igcm_n2d
            if (g_n2d.eq.0) then
               write(*,*) "chemthermos: Error; no N2D tracer !!!"
               write(*,*) "NO must be present if NO is in traceur.def"
               stop
            else
               nesptherm=nesptherm+1  
            endif
         endif   ! Of if(chemthermod.ge.2)
         
         !Ion chemistry
         if(chemthermod.eq.3) then
            g_co2plus=igcm_co2plus
            if (g_co2plus.eq.0 .and. chemthermod.ge.3) then
               write(*,*) "chemthermos: Error; no CO2+ tracer !!!"
               write(*,*) "CO2+ must be present if chemthermod=3"
               stop
            else
               nesptherm=nesptherm+1  
            endif
            g_oplus=igcm_oplus
            if (g_oplus.eq.0 .and. chemthermod.ge.3) then
               write(*,*) "chemthermos: Error; no O+ tracer !!!"
               write(*,*) "O+ must be present if chemthermod=3"
               stop
            else
               nesptherm=nesptherm+1  
            endif
            g_o2plus=igcm_o2plus
            if (g_o2plus.eq.0 .and. chemthermod.ge.3) then
               write(*,*) "chemthermos: Error; no O2+ tracer !!!"
               write(*,*) "O2+ must be present if chemthermod=3"
               stop
            else
               nesptherm=nesptherm+1  
            endif
            g_coplus=igcm_coplus
            if (g_coplus.eq.0) then
               write(*,*) "chemthermos: Error; no CO+ tracer !!!"
               write(*,*) "CO+ must be present if chemthermod=3"
               stop
            else
               nesptherm=nesptherm+1  
            endif
            g_cplus=igcm_cplus
            if (g_cplus.eq.0) then
               write(*,*) "chemthermos: Error; no C+ tracer !!!"
               write(*,*) "C+ must be present if chemthermod=3"
               stop
            else
               nesptherm=nesptherm+1  
            endif
            g_nplus=igcm_nplus
            if (g_nplus.eq.0) then
               write(*,*) "chemthermos: Error; no N+ tracer !!!"
               write(*,*) "N+ must be present if chemthermod=3"
               stop
            else
               nesptherm=nesptherm+1  
            endif
            g_noplus=igcm_noplus
            if (g_noplus.eq.0) then
               write(*,*) "chemthermos: Error; no NO+ tracer !!!"
               write(*,*) "NO+ must be present if chemthermod=3"
               stop
            else
               nesptherm=nesptherm+1  
            endif
            g_n2plus=igcm_n2plus
            if (g_n2plus.eq.0) then
               write(*,*) "chemthermos: Error; no N2+ tracer !!!"
               write(*,*) "N2+ must be present if chemthermod=3"
               stop
            else
               nesptherm=nesptherm+1  
            endif
            g_hplus=igcm_hplus
            if (g_hplus.eq.0) then
               write(*,*) "chemthermos: Error; no H+ tracer !!!"
               write(*,*) "H+ must be present if chemthermod=3"
               stop
            else
               nesptherm=nesptherm+1  
            endif
            g_hco2plus=igcm_hco2plus
            if (g_hco2plus.eq.0 .and. chemthermod.ge.3) then
               write(*,*) "chemthermos: Error; no HCO2+ tracer !!!"
               write(*,*) "HCO2+ must be present if chemthermod=3"
               stop
            else
               nesptherm=nesptherm+1  
            endif
            g_elec=igcm_elec
            if (g_elec.eq.0) then
               write(*,*) "chemthermos: Error; no e- tracer !!!"
               write(*,*) "e- must be present if chemthermod=3"
               stop
            else
               nesptherm=nesptherm+1  
            endif
         endif

         !Check if number of species is as expected
         select case(chemthermod)
         case(0)
            if(nesptherm.ne.11) then
               write(*,*)"chemthermos: Error:" 
               write(*,*)"Number of tracers not correct"
               write(*,*)"There are",nesptherm," tracers"
               write(*,*)"There should be 11 tracers"
               stop
            else
               write(*,*)"chemthermos:",nesptherm," species"
            endif
         case(1)            
            if(nesptherm.ne.12) then
               write(*,*)"chemthermos: Error:" 
               write(*,*)"Number of tracers not correct"
               write(*,*)"There are",nesptherm," tracers"
               write(*,*)"There should be 12 tracers"
               stop
            else
               write(*,*)"chemthermos:",nesptherm," species"
            endif
         case(2)
            if(nesptherm.ne.17) then
               write(*,*)"chemthermos: Error:" 
               write(*,*)"Number of tracers not correct"
               write(*,*)"There are",nesptherm," tracers"
               write(*,*)"There should be 17 tracers"
               stop
            else
               write(*,*)"chemthermos:",nesptherm," species"
            endif
         case(3)
            if(nesptherm.ne.28) then
               write(*,*)"chemthermos: Error:" 
               write(*,*)"Number of tracers not correct"
               write(*,*)"There are",nesptherm," tracers"
               write(*,*)"There should be 28 tracers"
               stop
            else
               write(*,*)"chemthermos:",nesptherm," species"
            endif
         end select
         firstcall= .false.
         write(*,*)'chemthermos: chemthermod=',chemthermod
      endif                     ! of if (firstcall)
      
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      !Allocate density vector
      allocate(rm(nlayer,nesptherm))

      do l=1,nlayer
        rm(l,i_co2)     = max(zycol(l,g_co2)*zdens(l),1.e-30)
        rm(l,i_co)      = max(zycol(l,g_co)*zdens(l),1.e-30) 
        rm(l,i_o)       = max(zycol(l,g_o)*zdens(l),1.e-30)
        rm(l,i_o1d)     = max(zycol(l,g_o1d)*zdens(l),1.e-30) 
        rm(l,i_o2)      = max(zycol(l,g_o2)*zdens(l),1.e-30)
        rm(l,i_h)       = max(zycol(l,g_h)*zdens(l),1.e-30) 
        rm(l,i_h2)      = max(zycol(l,g_h2)*zdens(l),1.e-30)
        rm(l,i_oh)      = max(zycol(l,g_oh)*zdens(l),1.e-30) 
        rm(l,i_ho2)     = max(zycol(l,g_ho2)*zdens(l),1.e-30)
        rm(l,i_h2o2)    = max(zycol(l,g_h2o2)*zdens(l),1.e-30) 
        rm(l,i_h2o)     = max(zycol(l,g_h2o_vap)*zdens(l),1.e-30)
        if(chemthermod.ge.1)  &
             rm(l,i_o3)      = max(zycol(l,g_o3)*zdens(l),1.e-30)
        if(chemthermod.ge.2) then
           rm(l,i_n2)      = max(zycol(l,g_n2)*zdens(l),1.e-30)
           rm(l,i_n)       = max(zycol(l,g_n)*zdens(l),1.e-30)
           rm(l,i_no)      = max(zycol(l,g_no)*zdens(l),1.e-30)
           rm(l,i_co)      = max(zycol(l,g_co)*zdens(l),1.e-30)
           rm(l,i_n2d)     = max(zycol(l,g_n2d)*zdens(l),1.e-30)
           rm(l,i_no2)     = max(zycol(l,g_no2)*zdens(l),1.e-30)
        endif
        if(chemthermod.eq.3) then
           rm(l,i_co2plus) = max(zycol(l,g_co2plus)*zdens(l),1.e-30)
           rm(l,i_oplus)   = max(zycol(l,g_oplus)*zdens(l),1.e-30)
           rm(l,i_o2plus)  = max(zycol(l,g_o2plus)*zdens(l),1.e-30)
           rm(l,i_coplus)  = max(zycol(l,g_coplus)*zdens(l),1.e-30)
           rm(l,i_cplus)   = max(zycol(l,g_cplus)*zdens(l),1.e-30)
           rm(l,i_nplus)   = max(zycol(l,g_nplus)*zdens(l),1.e-30)
           rm(l,i_noplus)  = max(zycol(l,g_noplus)*zdens(l),1.e-30)
           rm(l,i_n2plus)  = max(zycol(l,g_n2plus)*zdens(l),1.e-30)
           rm(l,i_hplus)   = max(zycol(l,g_hplus)*zdens(l),1.e-30)
           rm(l,i_hco2plus)= max(zycol(l,g_hco2plus)*zdens(l),1.e-30)
           rm(l,i_elec)    = max(zycol(l,g_elec)*zdens(l),1.e-30)
        endif
        if(chemthermod.gt.3.or.chemthermod.lt.0) then
           write(*,*)'chemthermos: bad value for chemthermod. Stop'
           stop
        endif
      enddo

      !Solar flux calculation

      !Photoabsorption coefficients
      call jthermcalc_e107(ig,nlayer,chemthermod,rm,nesptherm,ztemp,zlocal,zenit,zday)

      !Chemistry
      call paramfoto_compact  &
           (ig,nlayer,chemthermod,lswitch,ztemp,ptimestep,zenit,zlocal,rm,nesptherm)

      !Concentrations back to GCM
      do l=lswitch,nlayer
        zycol(l,g_co2)     = max(rm(l,i_co2)     / zdens(l) , 1.e-30)
        zycol(l,g_co)      = max(rm(l,i_co)      / zdens(l) , 1.e-30)
        zycol(l,g_o2)      = max(rm(l,i_o2)      / zdens(l) , 1.e-30)
        zycol(l,g_h2)      = max(rm(l,i_h2)      / zdens(l) , 1.e-30)
        zycol(l,g_h)       = max(rm(l,i_h)       / zdens(l) , 1.e-30)
        zycol(l,g_oh)      = max(rm(l,i_oh)      / zdens(l) , 1.e-30)
        zycol(l,g_ho2)     = max(rm(l,i_ho2)     / zdens(l) , 1.e-30)
        zycol(l,g_h2o_vap) = max(rm(l,i_h2o)     / zdens(l) , 1.e-30)
        zycol(l,g_h2o2)    = max(rm(l,i_h2o2)    / zdens(l) , 1.e-30)
        zycol(l,g_o1d)     = max(rm(l,i_o1d)     / zdens(l) , 1.e-30)
        zycol(l,g_o)       = max(rm(l,i_o)       / zdens(l) , 1.e-30)
        if(chemthermod.ge.1) &
             zycol(l,g_o3) = max(rm(l,i_o3)      / zdens(l) , 1.e-30)
        if(chemthermod.ge.2) then
           zycol(l,g_n2)      = max(rm(l,i_n2)      / zdens(l) , 1.e-30)
           zycol(l,g_n)       = max(rm(l,i_n)       / zdens(l) , 1.e-30)
           zycol(l,g_no)      = max(rm(l,i_no)      / zdens(l) , 1.e-30)
           zycol(l,g_no2)     = max(rm(l,i_no2)     / zdens(l) , 1.e-30)
           zycol(l,g_n2d)     = max(rm(l,i_n2d)     / zdens(l) , 1.e-30)
        endif
        if(chemthermod.ge.3) then
           zycol(l,g_co2plus) = max(rm(l,i_co2plus) / zdens(l) , 1.e-30)
           zycol(l,g_oplus)   = max(rm(l,i_oplus)   / zdens(l) , 1.e-30)
           zycol(l,g_o2plus)  = max(rm(l,i_o2plus)  / zdens(l) , 1.e-30)
           zycol(l,g_coplus)  = max(rm(l,i_coplus)  / zdens(l) , 1.e-30)
           zycol(l,g_cplus)   = max(rm(l,i_cplus)   / zdens(l) , 1.e-30)
           zycol(l,g_nplus)   = max(rm(l,i_nplus)   / zdens(l) , 1.e-30)
           zycol(l,g_noplus)  = max(rm(l,i_noplus)  / zdens(l) , 1.e-30)
           zycol(l,g_n2plus)  = max(rm(l,i_n2plus)  / zdens(l) , 1.e-30)
           zycol(l,g_hplus)   = max(rm(l,i_hplus)   / zdens(l) , 1.e-30)
           zycol(l,g_hco2plus)= max(rm(l,i_hco2plus)/ zdens(l) , 1.e-30)
           zycol(l,g_elec)    = max(rm(l,i_elec)    / zdens(l) , 1.e-30)
        endif
        em_no(l)=Pno(l,45)
        em_o2(l)=Po2(l,10)*0.75
      enddo                     !nlayer

      !Deallocations
      deallocate(rm)

      return
      end 




