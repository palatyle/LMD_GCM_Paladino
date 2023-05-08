      SUBROUTINE euvheat(nlon, nlev,nqmx, pt,pplev,pplay,zzlay, &
        	   mu0,ptimestep,ptime,zday,       &
                   pq, pdq, pdteuv)

        use chemparam_mod
	use dimphy
        use conc, only:  rnew, cpnew

      IMPLICIT NONE
!=======================================================================
!   subject:
!   --------
!   Computing heating rate due to EUV absorption
!
!   author:  MAC 2002
!   ------
!
!   input:
!   ----- 
!   mu0(klon)           
!   pplay(ngrid,nlayer)   pressure at middle of layers (Pa)
!   zzlay                 ! altitude at the middle of the layers (m)
!
!   output:
!   -------
!
!   pdteuv(ngrid,nlayer)      Heating rate (K/s)
!
!=======================================================================
!
!    0.  Declarations :
!    ------------------
!
#include "YOMCST.h"
#include "clesphys.h"
#include "param.h"
#include "param_v4.h"
!#include "chimiedata.h"
#include "mmol.h"
!-----------------------------------------------------------------------
!    Input/Output
!    ------------


      integer :: nlon
      integer :: nlev 
      integer :: nqmx

      real :: pt(nlon,nlev)
      real :: pplev(nlon,nlev+1)
      real :: pplay(nlon,nlev)
      real :: zzlay(nlon,nlev)
      real :: mu0(nlon)

      real :: ptimestep,ptime
      real :: zday
      real :: pq(nlon,nlev,nqmx)
      real :: pdq(nlon,nlev,nqmx)
      real :: pdteuv(nlon,nlev)
!
!    Local variables :
!    -----------------

      integer,save :: nespeuv=17    ! Number of species considered (11, 12 or 17 (with nitrogen))
      integer,save :: nspeuv_vgcm    ! Number of species considered currently considered into VGCM


      INTEGER :: l,ig,n
      integer,save :: euvmod = 0     !0: 4 (main) species  1: O3 chemistry 2: N chemistry, 3: C/O/H
      real, allocatable, save :: rm(:,:)   !  number density (cm-3)
      real :: zq(nlon,nlev,nqmx) ! local updated tracer quantity
      real :: zt(nlon,nlev)      ! local updated atmospheric temperature
      real :: zlocal(nlev)
      real :: zenit
      real :: jtot(nlev)
      real :: dens		! Total number density (cm-3)
      real :: tx(nlev)
       
! tracer indexes for the EUV heating:
!!! ATTENTION. These values have to be identical to those in hrterm.F90
!!! If the values are changed there, the same has to be done here  !!!
      
      integer,parameter :: ix_co2=1
      integer,parameter :: ix_o=3
      integer,parameter :: ix_co=4
      integer,parameter :: ix_n2=13


! Tracer indexes in the GCM:
      integer,save :: g_co2=0
      integer,save :: g_o=0
      integer,save :: g_co=0
      integer,save :: g_n2=0
      
      logical,save :: firstcall=.true.

! Initializations and sanity checks:


      if (firstcall) then
         nspeuv_vgcm=0
!        ! identify the indexes of the tracers we'll need
         g_co2=i_co2
         if (g_co2.eq.0) then
            write(*,*) "euvheat: Error; no CO2 tracer !!!"
            write(*,*) "CO2 is always needed if calleuv=.true."
            stop
         else
            nspeuv_vgcm=nspeuv_vgcm+1
         endif
         g_o=i_o
         if (g_o.eq.0) then
            write(*,*) "euvheat: Error; no O tracer !!!"
!            write(*,*) "O is always needed if calleuv=.true."
            stop
         else
           nspeuv_vgcm=nspeuv_vgcm+1 
         endif
         g_co=i_co
         if (g_co.eq.0) then
            write(*,*) "euvheat: Error; no CO tracer !!!"
!            write(*,*) "CO is always needed if calleuv=.true."
            stop
         else
            nspeuv_vgcm=nspeuv_vgcm+1
         endif
         ! n2
         g_n2=i_n2
            if (g_n2.eq.0) then
               write(*,*) "euvheat: Error; no N2 tracer !!!"
!               write(*,*) "N2 needed if NO is in traceur.def"
               stop
            else
               nspeuv_vgcm=nspeuv_vgcm+1
            endif

!         g_o2=igcm_o2
!         if (g_o2.eq.0) then
!            write(*,*) "euvheat: Error; no O2 tracer !!!"
!            write(*,*) "O2 is always needed if calleuv=.true."
!            stop
!         else
!            nespeuv=nespeuv+1
!         endif
!         g_h2=igcm_h2
!         if (g_h2.eq.0) then
!            write(*,*) "euvheat: Error; no H2 tracer !!!"
!            write(*,*) "H2 is always needed if calleuv=.true."
!            stop
!         else
!            nespeuv=nespeuv+1
!         endif
!         g_oh=igcm_oh
!         if (g_oh.eq.0) then
!            write(*,*) "euvheat: Error; no OH tracer !!!"
!            write(*,*) "OH must always be present if thermochem=T"
!            stop
!         else
!            nespeuv=nespeuv+1  
!         endif
!         g_ho2=igcm_ho2
!         if (g_ho2.eq.0) then
!            write(*,*) "euvheat: Error; no HO2 tracer !!!"
!            write(*,*) "HO2 must always be present if thermochem=T"
!            stop
!         else
!            nespeuv=nespeuv+1  
!         endif
!         g_h2o2=igcm_h2o2
!         if (g_h2o2.eq.0) then
!            write(*,*) "euvheat: Error; no H2O2 tracer !!!"
!            write(*,*) "H2O2 is always needed if calleuv=.true."
!            stop
!         else
!            nespeuv=nespeuv+1
!         endif
!         g_h2o=igcm_h2o_vap
!         if (g_h2o.eq.0) then
!            write(*,*) "euvheat: Error; no water vapor tracer !!!"
!            write(*,*) "H2O is always needed if calleuv=.true."
!            stop
!         else
!            nespeuv=nespeuv+1
!         endif
!         g_o1d=igcm_o1d
!         if (g_o1d.eq.0) then
!            write(*,*) "euvheat: Error; no O1D tracer !!!"
!            write(*,*) "O1D must always be present if thermochem=T"
!            stop
!         else
!            nespeuv=nespeuv+1  
!         endif
!         g_h=igcm_h
!         if (g_h.eq.0) then
!            write(*,*) "euvheat: Error; no H tracer !!!"
!            write(*,*) "H is always needed if calleuv=.true."
!            stop
!         else
!            nespeuv=nespeuv+1
!         endif
         
!         euvmod = 3            !Default: C/O/H chemistry 
!         !Check if O3 is present
!         g_o3=igcm_o3
!         if (g_o3.eq.0) then
!            write(*,*) "euvheat: Error; no O3 tracer !!!"
!            write(*,*) "O3 must be present if calleuv=.true."
!            stop
!         else 
!            nespeuv=nespeuv+1
!            euvmod=1
!         endif

         !Nitrogen species
         !NO is used to determine if N chemistry is wanted
         !euvmod=2 -> N chemistry
!         g_no=igcm_no
!         if (g_no.eq.0) then
!            write(*,*) "euvheat: no NO tracer"
!            write(*,*) "No N species in UV heating"
!         else if(g_no.ne.0) then
!            nespeuv=nespeuv+1
!            euvmod=2
!         endif
         ! N
!         g_n=igcm_n
!         if(euvmod == 2) then
!            if (g_n.eq.0) then
!               write(*,*) "euvheat: Error; no N tracer !!!"
!               write(*,*) "N needed if NO is in traceur.def"
!               stop
!            else if(g_n.ne.0) then
!               nespeuv=nespeuv+1
!            endif
!         else
!            if(g_n /= 0) then
!               write(*,*) "euvheat: Error: N present, but NO not!!!"
!               write(*,*) "Both must be in traceur.def"
!               stop
!            endif
!         endif   !Of if(euvmod==2)
!         !NO2
!         g_no2=igcm_no2
!         if(euvmod == 2) then
!            if (g_no2.eq.0) then
!               write(*,*) "euvheat: Error; no NO2 tracer !!!"
!               write(*,*) "NO2 needed if NO is in traceur.def"
!               stop
!            else if(g_no2.ne.0) then
!               nespeuv=nespeuv+1
!            endif
!         else
!            if(g_no2 /= 0) then
!               write(*,*) "euvheat: Error: NO2 present, but NO not!!!"
!               write(*,*) "Both must be in traceur.def"
!               stop
!            endif
!         endif   !Of if(euvmod==2)
         !N2D
!         g_n2d=igcm_n2d
!         if(euvmod == 2) then
!            if (g_n2d.eq.0) then
!               write(*,*) "euvheat: Error; no N2D tracer !!!"
!               write(*,*) "N2D needed if NO is in traceur.def"
!               stop
!            else
!               nespeuv=nespeuv+1  
!            endif
!         else
!            if(g_n2d /= 0) then
!               write(*,*) "euvheat: Error: N2D present, but NO not!!!"
!               write(*,*) "Both must be in traceur.def"
!               stop
!            endif
!         endif   !Of if(euvmod==2)

         !Check if nespeuv is appropriate for the value of euvmod
!         select case(euvmod)
!         case(0)
!            if(nespeuv.ne.11) then
!               write(*,*)'euvheat: Wrong number of tracers!'
!               stop
!            else
!               write(*,*)'euvheat: Computing absorption by',nespeuv, &
!                    ' species'
!            endif
!         case(1)
!            if(nespeuv.ne.12) then
!               write(*,*)'euvheat: Wrong number of tracers!',nespeuv
!               stop
!            else
!               write(*,*)'euvheat: Computing absorption by',nespeuv,  &
!                    ' species'
!            endif
!         case(2)
!            if(nespeuv.ne.17) then
!               write(*,*)'euvheat: Wrong number of tracers!'
!               stop
!            else
!               write(*,*)'euvheat: Computing absorption by',nespeuv,  &
!                    ' species'
!            endif
!         end select


      !Allocate density vector
      allocate(rm(nlev,nespeuv))

         firstcall= .false.
      endif                     ! of if (firstcall)

!      write(*,*),  "CHECK n species currently used into VGCM",  nspeuv_vgcm


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ! build local updated values of tracers (if any) and temperature

      do l=1,nlev
        do ig=1,nlon

! chemical species
          zq(ig,l,g_co2)=pq(ig,l,g_co2)  
          zq(ig,l,g_co)=pq(ig,l,g_co)   
          zq(ig,l,g_o)=pq(ig,l,g_o)     
          zq(ig,l,g_n2)=pq(ig,l,g_n2)    

! atmospheric temperature
         zt(ig,l)=pt(ig,l)

!      write(*,*),  "CHECK update densities L332 euv",   zq(ig,l,g_co2)


        enddo
      enddo
      
      !Solar flux calculation
      if(solvarmod.eq.0) call flujo(solarcondate) 

                                         ! Not recommended for long runs 
                                         ! (e.g. to build the MCD), as the
                                         ! solar conditions at the end will 
                                         ! be different to the ones initially
                                         ! set
      
      do ig=1,nlon
         zenit=acos(mu0(ig))*180./acos(-1.)  !convers from rad to deg
         
         do l=1,nlev
         
!   Conversion to number density
	   	    
!!!  use R specific = R/MolarMass

            dens=pplay(ig,l)/(rnew(ig,l)*zt(ig,l)) / 1.66e-21   ! [g mol-1] [cm-3] 	   

            rm(l,ix_co2)  = zq(ig,l,g_co2) * dens / M_tr(g_co2)   ! [cm-3] 
            rm(l,ix_o)    = zq(ig,l,g_o)   * dens / M_tr(g_o)
            rm(l,ix_co)   = zq(ig,l,g_co) * dens  / M_tr(g_co)
            rm(l,ix_n2)   = zq(ig,l,g_n2) * dens  / M_tr(g_n2)

!      	    write(*,*),  "CHECK n density", l, rm(l,ix_co2)


         enddo

!        zlocal(1)=-log(pplay(ig,1)/pplev(ig,1))
!     &            *Rnew(ig,1)*zt(ig,1)/g
         zlocal(1)=zzlay(ig,1)
         zlocal(1)=zlocal(1)/1000.    ! conversion m ---> km
         tx(1)=zt(ig,1)
	 
         do l=2,nlev
            tx(l)=zt(ig,l)
            zlocal(l)=zzlay(ig,l)/1000.
         enddo

        !Routine to calculate the UV heating
         call hrtherm (ig,euvmod,rm,nespeuv,tx,zlocal,zenit,zday,jtot)
	 

        !Calculates the UV heating from the total photoabsorption coefficient
        do l=1,nlev
!       jtot conversion from erg/(s*cm3) ---> J/(s*m3)
          pdteuv(ig,l)=euveff*jtot(l)/10.                  &
               /(cpnew(ig,l)*pplay(ig,l)/(rnew(ig,l)*zt(ig,l)))

        enddo	
      enddo  ! of do ig=1,nlon

      !Deallocations
      !deallocate(rm)

      return
      end 
