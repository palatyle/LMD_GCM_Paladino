      subroutine calchim(ptimestep,pplay,pplev,pt,pdt,dist_sol,mu0,
     $                   zzlay,zday,pq,pdq,dqchim,dqschim,dqcloud,
     $                   dqscloud)
c
      implicit none
c
c=======================================================================
c
c   subject:
c   --------
c
c  Prepare the call for the photochemical module, and send back the
c  tendencies from photochemistry in the chemical species mass mixing ratios
c
c   Author:   Sebastien Lebonnois (08/11/2002)
c   -------
c    update 12/06/2003 for water ice clouds and compatibility with dust
c    update 07/2003 for coupling with thermosphere (Monica Angelats-i-Coll)
c    update 03/05/2005 cosmetic changes (Franck Lefevre)
c
c   Arguments:
c   ----------
c
c  Input:
c
c    ptimestep                  timestep (s)
c    pplay(ngridmx,nlayermx)    Pressure at the middle of the layers (Pa)
c    pplev(ngridmx,nlayermx+1)  Intermediate pressure levels (Pa)
c    pt(ngridmx,nlayermx)       Temperature (K)
c    pdt(ngridmx,nlayermx)      Temperature tendency (K)
c    dist_sol                   distance of the sun (AU)
c    mu0(ngridmx)               cos of solar zenith angle (=1 when sun at zenith)
c    pq(ngridmx,nlayermx,nqmx)  Advected fields, ie chemical species here
c    pdq(ngridmx,nlayermx,nqmx) Previous tendencies on pq
c
c  Output:
c
c    dqchim(ngridmx,nlayermx,nqmx) ! tendencies on pq due to chemistry
c    dqschim(ngridmx,nqmx)         ! tendencies on qsurf 
c
c=======================================================================

c    Declarations :
c    --------------

#include "dimensions.h"
#include "dimphys.h"
#include "chimiedata.h"
#include "tracer.h"
#include "comcstfi.h"
#include "callkeys.h"
#include "fisice.h"
#include "conc.h"

c Arguments :
c -----------

c   inputs:
c   -------

      real    ptimestep
      real    pplay(ngridmx,nlayermx)    ! pressure at the middle of the layers
      real    zzlay(ngridmx,nlayermx)    ! pressure at the middle of the layers
      real    pplev(ngridmx,nlayermx+1)  ! intermediate pressure levels
      real    pt(ngridmx,nlayermx)       ! temperature
      real    pdt(ngridmx,nlayermx)      ! temperature tendency
      real    dist_sol                   ! distance of the sun (AU)
      real    mu0(ngridmx)       ! cos of solar zenith angle (=1 when sun at zenith)
      real    pq(ngridmx,nlayermx,nqmx)  ! tracers mass mixing ratio
      real    pdq(ngridmx,nlayermx,nqmx) ! previous tendencies
      real    zday                       ! date (time since Ls=0, in martian days)

c   outputs:
c   --------

      real dqchim(ngridmx,nlayermx,nqmx) ! tendencies on pq due to chemistry
      real dqschim(ngridmx,nqmx)         ! tendencies on qsurf 
      real dqcloud(ngridmx,nlayermx,nqmx)! tendencies on pq due to condensation
      real dqscloud(ngridmx,nqmx)        ! tendencies on qsurf 

c Local variables :
c -----------------

      character*5 str5
      integer  ig,l,iq, i_co2, i_o
      integer  foundswitch, lswitch
      real zq(ngridmx,nlayermx,nqmx) ! pq+pdq*ptimestep before chemistry
                                     ! new mole fraction after
      real colden(ngridmx,nqmx)      ! Column densities (cm-2)
      real zt(ngridmx,nlayermx)      ! temperature
c
c  for each column of atmosphere:
c
      real zpress(nlayermx)       !  Pressure (mbar)
      real zdens(nlayermx)        !  Density  (cm-3)
      real ztemp(nlayermx)        !  Temperature (K)
      real zlocal(nlayermx)       !  Altitude (km)
      real zycol(nlayermx,nqmx)   !  Composition (mole fractions)
      real szacol                 !  Solar zenith angle
      real jo3(nlayermx)          !  Photodissociation rate O3->O1D (s-1)
c
c for output:
c
      real zdens3d(ngridmx,nlayermx) !  Density  (cm-3)
      real jo3_3d(ngridmx,nlayermx)  !  Photodissociation rate O3->O1D (s-1)
      real surfice(ngridmx,nlayermx) !  Surface of ice particules (um2/cm3)
      logical output

      logical firstcall
      data    firstcall/.true./
      save    firstcall
c
c  scheme A: 1 ; scheme B: 2
c
      integer scheme
      data  scheme/2/
c
c=======================================================================
c     initialization of the chemistry (first call only)
c=======================================================================
c
      if (firstcall) then
c
         if (photochem) then
            print*,'INIT CHEMISTRY'
            if (scheme  .eq.  1) then
               print*,'Scheme A : A METTRE A JOUR !!'
               stop
c              call init_chimie_A
            else
               print*,'Scheme B'
               call init_chimie_B
            end if
         end if
         firstcall = .false.
      end if
c
c=======================================================================
c     loop over grid
c=======================================================================
c
      do ig = 1,ngridmx
c
c     local updates
c 
         foundswitch = 0
         do l = 1,nlayermx
            do iq = nqchem_min,nqmx
               zq(ig,l,iq) = pq(ig,l,iq) + pdq(ig,l,iq)*ptimestep
               zt(ig,l)    = pt(ig,l)    + pdt(ig,l)   *ptimestep
               zycol(l,iq) = zq(ig,l,iq) * mmean(ig,l)/mmol(iq)
            enddo
            zpress(l) = pplay(ig,l)/100.
            ztemp(l)  = zt(ig,l)
            zdens(l)  = zpress(l)/(kb*1.e4*ztemp(l))
            zlocal(l) = zzlay(ig,l)/1000.
c
c     search for switch index between regions
c
            if (photochem .and. thermochem) then
               if (foundswitch .eq. 0 .and. pplay(ig,l).lt.1.e-3) then
                  lswitch = l
                  foundswitch=1
               end if
            end if
            if ( .not. photochem) then
               lswitch = 22
            end if
            if (.not.  thermochem) then
               lswitch = min(33,nlayermx+1)
            end if
c
c     ice surface area  in microns^2/cm^3
c
c     = 4 pi r^2 * [ zq * mugaz/NA / (rhoice*4/3 pi r^3) ] *zdens
c     = 3/r * [ zq * mugaz/NA / rhoice ] *zdens
c     with r in microns, rhoice = 0.92e-12 g microns^-3 and zdens in cm^-3
c
            if (iceparty) then
               zycol(l,nqmx-1) = (3.e-6/rice(ig,l))*zq(ig,l,nqmx-1)
     $                           *(mugaz/6.022e23)*zdens(l)/0.92e-12
c              write(*,*) "rice=",rice(ig,l)," m / zdens=",zdens(l),
c    $        " cm-3 / icesurf=",zycol(l,nqmx-1)," microns^2/cm^3"
               surfice(ig,l) = zycol(l,nqmx-1)
            end if
c
         end do
c
         szacol = acos(mu0(ig))*180./pi
c
c=======================================================================
c     call chemical subroutine
c=======================================================================
c
        if (photochem) then
           if (scheme .eq. 1) then
              print*,'Scheme A : A METTRE A JOUR !!'
c             call photochemist_A(zycol,szacol,ptimestep,
c    $                            zpress,ztemp,zdens,dist_sol)
           else
              call photochemist_B(lswitch,zycol,szacol,ptimestep,
     $                            zpress,ztemp,zdens,dist_sol,jo3)
           end if
        end if
        if (thermochem) then
           call chemthermos(ig,lswitch,zycol,ztemp,zdens,zpress,
     $                      zlocal,szacol,ptimestep,zday)
        end if
c
c=======================================================================
c     tendencies
c=======================================================================
c
c     must be 0. for water ice:
c
         if (iceparty) then
            do l = 1,nlayermx
               dqchim(ig,l,nqmx-1) = 0.
            end do
         end if
c
c     tendency for CO2 = - sum of others for lower atmosphere
c     tendency for O   = - sum of others for upper atmosphere
c
         do l = 1,nlayermx
            if (l .lt. lswitch) then
               do iq = nqchem_min,nqmx
                  if ((noms(iq) .ne. "co2") .and. 
     $                (noms(iq) .ne. "ice")) then
                     dqchim(ig,l,iq) = (zycol(l,iq)*mmol(iq)/mmean(ig,l)
     $                                - zq(ig,l,iq))/ptimestep
                  else if (noms(iq) .eq. "co2") then
                     i_co2 = iq
                     dqchim(ig,l,iq) = 0.
                  end if
                  dqschim(ig,iq) = 0.
               end do
               do iq = nqchem_min,nqmx
                  if (noms(iq) .ne. "co2") then
                     dqchim(ig,l,i_co2) = dqchim(ig,l,i_co2) 
     $                                  - dqchim(ig,l,iq)
                  end if
                end do
             else if (l .ge. lswitch) then
                do iq = nqchem_min,nqmx
                   if ((noms(iq).ne."o") .and.(noms(iq) .ne."ice")) then
                      dqchim(ig,l,iq) = (zycol(l,iq)*mmol(iq)
     $                                  /mmean(ig,l)
     $                                 - zq(ig,l,iq))/ptimestep
                   else if (noms(iq) .eq. "o") then
                      i_o = iq
                      dqchim(ig,l,iq) = 0.
                   end if
                enddo
                do iq = nqchem_min,nqmx
                   if (noms(iq) .ne. "o") then
                      dqchim(ig,l,i_o) = dqchim(ig,l,i_o) 
     $                                 - dqchim(ig,l,iq)
                   end if
                end do
             end if
          end do
c
c     dust
c
          if (nqchem_min .gt. 1) then
             do iq = 1,nqchem_min-1
                do l = 1,nlayermx
                   dqchim(ig,l,iq) = 0.
                end do
                dqschim(ig,iq) = 0.
             end do
          end if
c
c     condensation of h2o2
c
          call perosat(ig,ptimestep,pplev,pplay,
     $                 ztemp,zycol,dqcloud,dqscloud)
c
c     for outputs
c
          do iq = nqchem_min,nqmx
             colden(ig,iq) = 0.
             do l = 1,nlayermx
c
c     column density converted in cm-2
c     pplev en pa, mugaz en g.mol-1 et g en m.s-2
c     not for ice
c
                if (noms(iq) .ne. "h2o2") then
                  colden(ig,iq) = colden(ig,iq) + zycol(l,iq)
     $                         *6.022e22*(pplev(ig,l)-pplev(ig,l+1))
     $                         /(mmean(ig,l)*g)
                else   ! for H2O2, remove condensation from zycol
                  colden(ig,iq) = colden(ig,iq) + (zycol(l,iq) +
     $               dqcloud(ig,l,iq)*ptimestep*mmean(ig,l)/mmol(iq))
     $                         *6.022e22*(pplev(ig,l)-pplev(ig,l+1))
     $                         /(mmean(ig,l)*g)
                end if
c
c     local densities, for outputs (put in zq)
c     not for ice
c
                zq(ig,l,iq) = zycol(l,iq)*zdens(l)
c                        for H2O2, remove condensation from zycol
                if (noms(iq) .eq. "h2o2") then
                   zq(ig,l,iq) = zdens(l)*(zycol(l,iq) +
     $               dqcloud(ig,l,iq)*ptimestep*mmean(ig,l)/mmol(iq))
                end if
             end do
          end do
c
c     density and j(o3->o1d), for outputs
c
          zdens3d(ig,1) = zdens(1)
          jo3_3d(ig,1) = jo3(1)
          do l = 2,nlayermx
             zdens3d(ig,l) = zdens(l)
             jo3_3d(ig,l) = jo3(l)
          end do
c
c=======================================================================
c     end of loop over grid
c=======================================================================
c
      end do
c
c=======================================================================
c     write outputs
c=======================================================================
c
      output = .true.

      if (output) then
         if (ngridmx .gt. 1) then
c           call writediagfi(ngridmx,'dens','atm dens.','cm-3',3,zdens3d(1,1))
c           call writediagfi(ngridmx,'jo3','j o3->o1d','s-1',3,jo3_3d(1,1))
c           call writediagfi(ngridmx,'sice','ice surf.','um2/cm3',3,surfice(1,1))
            do iq = nqchem_min,nqmx
               if (noms(iq) .ne. "ice") then
                  write(str5(1:5),'(a5)') noms(iq)
c                  call writediagfi(ngridmx,'n_'//str5,'density',
c     $                             'cm-3',3,zq(1,1,iq))
cc                 call writediagfi(ngridmx,'dqch_'//str5,'density','cm-3',3,dqchim(1,1,iq))
cc                 if (noms(iq) .eq. "h2o2" .or. noms(iq) .eq. "h2o") then
cc                    call writediagfi(ngridmx,'cl_'//str5,'density','cm-3',3,dqcloud(1,1,iq))
cc                 end if
c                  call writediagfi(ngridmx,'c_'//str5,'col. dens.',
c     $                             'cm-2',2,colden(1,iq))
               end if
            end do
c
            if (callstats) then
c
c              convert to mole.cm-2 for the column densities
c
               do iq = nqchem_min,nqmx
                  do ig = 1,ngridmx
                     colden(ig,iq) = colden(ig,iq)/6.022e23
                  end do   
               end do 
c
c              call wstats(ngridmx,"jo3","jo3->o1d","s-1",3,jo3_3d)
c
c               do iq = nqchem_min,nqmx
c                  if (noms(iq) .ne. "ice") then
c                     write(str5(1:5),'(a5)') noms(iq)
c                     call wstats(ngridmx,"n_"//str5,"density",
c     &                           "cm-3",3,zq(1,1,iq))
c                     call wstats(ngridmx,"c_"//str5,"col. dens.",
c     &                           "mol cm-2",2,colden(1,iq))
c                  end if
c               end do
            end if
         end if
c
      endif
c
      end
