c*****************************************************************
c
c     Photochemical routines for Scheme B
c
c     Author: Franck Lefevre
c     ------
c
c     update 12/06/2003: compatibility with dust and iceparty
c        (S. Lebonnois)
c
c*****************************************************************
c
      subroutine photochemist_B(lswitch,zycol, sza, ptimestep, press,
c     $                          temp, dens, dist_sol)
     $                          temp, dens, dist_sol, jo3)
c
c
      implicit none
c
#include "dimensions.h"
#include "dimphys.h"
#include "chimiedata.h" 
#include "callkeys.h"
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     input/output:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      real zycol(nlayermx,nqmx)  ! chemical species volume mixing ratio
c if iceparty, zycol(l,nqmx-1) is ice surface area (in microns^2/cm^3)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     inputs:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      real sza                   ! solar zenith angle (deg)
      real ptimestep             ! physics timestep (s)
      real press(nlayermx)       ! pressure (hpa)
      real temp(nlayermx)        ! temperature (k)
      real dens(nlayermx)        ! density (cm-3)
      real dist_sol              ! sun distance (au) 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     output:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      real jo3(nlayermx)  ! photodissociation rate O3 -> O1D
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     local:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      integer phychemrat         ! ratio physics timestep/chemistry timestep
      integer istep
      integer i_o3,lev
      integer nesp
      integer lswitch
c
      parameter (nesp = 15)      ! number of species in the chemistry
c
      real ctimestep             ! chemistry timestep (s) 
      real rm(nlayermx,nesp)     ! species volume mixing ratio 
      real j(nlayermx,nd)        ! interpolated photolysis rates (s-1)
      real icesurf(nlayermx)     ! ice surface area (cm^2/cm^3)
      real rmo3(nlayermx)        ! ozone mixing ratio
c
c     reaction rates
c
      real a001(nlayermx), a002(nlayermx), a003(nlayermx)
      real b001(nlayermx), b002(nlayermx), b003(nlayermx), 
     $     b004(nlayermx), b005(nlayermx), b006(nlayermx)
      real c001(nlayermx), c002(nlayermx), c003(nlayermx), 
     $     c004(nlayermx), c005(nlayermx), c006(nlayermx), 
     $     c007(nlayermx), c008(nlayermx), c009(nlayermx), 
     $     c010(nlayermx), c011(nlayermx), c012(nlayermx),
     $     c013(nlayermx), c014(nlayermx), c015(nlayermx), 
     $     c016(nlayermx), c017(nlayermx), c018(nlayermx)
      real d001(nlayermx), d002(nlayermx), d003(nlayermx)
      real e001(nlayermx), e002(nlayermx)
      real h001(nlayermx), h002(nlayermx)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ctimestep : chemistry timestep (s)                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      ctimestep = 10.*60.
c
      phychemrat = int(ptimestep/ctimestep)
c     print*,"PHYCHEMRAT = ",phychemrat
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     initialisation of chemical species and families            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call gcmtochim(zycol, lswitch, nesp, rm) 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ice surface area for heterogenous chemistry
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (iceparty) then
       do lev=1,lswitch-1
        icesurf(lev)= zycol(lev,nqmx-1)*1.e-8 ! microns^2 -> cm^2
       enddo
      else
       do lev=1,lswitch-1
        icesurf(lev)=0.0e0
       enddo
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     compute reaction rates                                     c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call chemrates(lswitch, dens, press, temp, icesurf,
     $               a001, a002, a003,
     $               b001, b002, b003, b004, b005, b006,
     $               c001, c002, c003, c004, c005, c006,
     $               c007, c008, c009, c010, c011, c012,
     $               c013, c014, c015, c016, c017, c018,
     $               d001, d002, d003,
     $               e001, e002,
     $               h001, h002)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     interpolation of photolysis rates in the lookup table      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      i_o3 = 5 
c
      do lev=1,lswitch-1
        rmo3(lev)=rm(lev,i_o3)
      enddo

      call phot(lswitch, press, temp, sza, dist_sol, rmo3, j)
c
      do lev=1,lswitch-1
        jo3(lev) = j(lev,5)
      enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     loop over time within the physical timestep                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do istep = 1,phychemrat
         call chimie_B(lswitch,nesp, rm, j, dens, ctimestep, 
     $                 press, temp, sza, dist_sol,
     $                 a001, a002, a003,
     $                 b001, b002, b003, b004, b005, b006,
     $                 c001, c002, c003, c004, c005, c006,
     $                 c007, c008, c009, c010, c011, c012,
     $                 c013, c014, c015, c016, c017, c018,
     $                 d001, d002, d003,
     $                 e001, e002,
     $                 h001, h002)
c
c        print*,'appel chimie ok iteration ', istep
c
      end do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     save chemical species for the gcm                          c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call chimtogcm(zycol, lswitch, nesp, rm) 
c
      return
      end
c
c*****************************************************************
c
      subroutine chimie_B(lswitch, nesp, rm, j, dens, dt, 
     $                    press, t, sza, dist_sol,
     $                    a001, a002, a003,
     $                    b001, b002, b003, b004, b005, b006,
     $                    c001, c002, c003, c004, c005, c006,
     $                    c007, c008, c009, c010, c011, c012,
     $                    c013, c014, c015, c016, c017, c018,
     $                    d001, d002, d003,
     $                    e001, e002,
     $                    h001, h002)
c
c*****************************************************************
c
      implicit none 
c
#include "dimensions.h"
#include "dimphys.h"
#include "chimiedata.h" 
#include "callkeys.h" 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     input/output:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      integer lswitch         ! interface level between chemistries
      integer nesp            ! number of species
c
      real rm(nlayermx,nesp)  ! volume mixing ratios 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     inputs:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real dens(nlayermx)     ! density (cm-3)
      real dt                 ! chemistry timestep (s) 
      real j(nlayermx,nd)     ! interpolated photolysis rates (s-1) 
      real press(nlayermx)    ! pressure (hpa)
      real t(nlayermx)        ! temperature (k)
      real sza                ! solar zenith angle (deg)
      real dist_sol           ! sun distance (au) 
c
c     reaction rates
c
      real a001(nlayermx), a002(nlayermx), a003(nlayermx)
      real b001(nlayermx), b002(nlayermx), b003(nlayermx), 
     $     b004(nlayermx), b005(nlayermx), b006(nlayermx)
      real c001(nlayermx), c002(nlayermx), c003(nlayermx), 
     $     c004(nlayermx), c005(nlayermx), c006(nlayermx), 
     $     c007(nlayermx), c008(nlayermx), c009(nlayermx), 
     $     c010(nlayermx), c011(nlayermx), c012(nlayermx),
     $     c013(nlayermx), c014(nlayermx), c015(nlayermx), 
     $     c016(nlayermx), c017(nlayermx), c018(nlayermx)
      real d001(nlayermx), d002(nlayermx), d003(nlayermx)
      real e001(nlayermx), e002(nlayermx)
      real h001(nlayermx), h002(nlayermx)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     local:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      integer HETERO
      parameter (HETERO=0)

      integer iesp, iter, l, niter
      integer i_co2, i_co, i_o2, i_h2, i_h2o, i_h2o2, i_hox, i_ox,
     $        i_o1d, i_o, i_o3, i_h, i_oh, i_ho2, i_n2
      integer j_o2_o, j_o2_o1d, j_co2_o, j_co2_o1d,
     $        j_o3_o1d, j_o3_o, j_h2o, j_hdo, j_h2o2,
     $        j_ho2, j_no2
c
      parameter (niter = 5)          ! iterations in the chemical scheme
c
      real*8 cc0(nlayermx,nesp)        ! initial number density (cm-3)
      real*8 cc(nlayermx,nesp)         ! final number density (cm-3)
      real*8 co2(nlayermx)             ! co2 number density (cm-3)
      real*8 nox(nlayermx)             ! nox number density (cm-3)
      real*8 no(nlayermx)              ! no  number density (cm-3)
      real*8 no2(nlayermx)             ! no2 number density (cm-3)
      real*8 production(nlayermx,nesp) ! production rate
      real*8 loss(nlayermx,nesp)       ! loss rate
c
      real c007h(nlayermx)           ! c007 rate modified
c                                    ! by heterogenous reactions h001 and h002
c
      real ro_o3, rh_ho2, roh_ho2, rno2_no
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     tracer numbering in the chemistry
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      i_co2    =  1
      i_co     =  2
      i_o      =  3
      i_o1d    =  4
      i_o2     =  5
      i_o3     =  6
      i_h      =  7
      i_h2     =  8
      i_oh     =  9
      i_ho2    = 10 
      i_h2o2   = 11
      i_h2o    = 12
      i_n2     = 13
      i_hox    = 14
      i_ox     = 15
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     numbering of photolysis rates
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      j_o2_o     =  1      ! o2 + hv   -> o + o
      j_o2_o1d   =  2      ! o2 + hv   -> o + o(1d)
      j_co2_o    =  3      ! co2 + hv  -> co + o
      j_co2_o1d  =  4      ! co2 + hv  -> co + o(1d)
      j_o3_o1d   =  5      ! o3 + hv   -> o2 + o(1d)
      j_o3_o     =  6      ! o3 + hv   -> o2 + o
      j_h2o      =  7      ! h2o + hv  -> h + oh
      j_hdo      =  8      ! hdo + hv  -> d + oh
      j_h2o2     =  9      ! h2o2 + hv -> oh + oh
      j_ho2      =  10     ! ho2 + hv  -> oh + o
      j_no2      =  11     ! no2 + hv  -> no + o
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     volume mixing ratio -> density conversion
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do iesp = 1,nesp
         do l = 1,lswitch-1
            cc0(l,iesp) = rm(l,iesp)*dens(l)
            cc(l,iesp)  = cc0(l,iesp)
         end do
      end do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     co2 and nox number densities (cm-3)    
c      
c     co2 mixing ratio: 0.953
c     nox mixing ratio: 6.e-10 (yung and demore, 1999)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do l = 1,lswitch-1
         co2(l) = cc(l,i_co2)
         nox(l) = 6.e-10*dens(l)
      end do 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     loop over iterations in the chemical scheme
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do iter = 1,niter
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     modification of c007
c        by heterogenous reactions h001 and h002
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (HETERO.eq.1) then
       do l = 1,lswitch-1
         c007h(l) = c007(l)
     $            + h001(l)/max(cc(l,i_oh), 1.e-30*dens(l))
     $            + h002(l)/max(cc(l,i_ho2),1.e-30*dens(l))
c         write(*,*) "C007 = ",c007(l)," / C007H = ",c007h(l)
       end do
      else
       do l = 1,lswitch-1
         c007h(l) = c007(l)
       end do
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        nox species
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        no2/no
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do l = 1,lswitch-1
c
         rno2_no = (d002(l)*cc(l,i_o3) + d003(l)*cc(l,i_ho2))
     $            /(j(l,j_no2) +
     $              d001(l)*max(cc(l,i_o),1.e-30*dens(l)))
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        no
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         no(l) = nox(l)/(1. + rno2_no)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        no2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         no2(l) = no(l)*rno2_no
c
      end do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        hox species
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        photochemical equilibrium for oh and ho2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        h/ho2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do l = 1,lswitch-1
c
         rh_ho2 = (c001(l)*cc(l,i_o)
     $           + c004(l)*cc(l,i_h)
     $           + c005(l)*cc(l,i_h)
     $           + c006(l)*cc(l,i_h)
     $           + c007h(l)*cc(l,i_oh)
     $           + 2.*c008(l)*cc(l,i_ho2)
     $           + c015(l)*cc(l,i_o3)
     $           + 2.*c016(l)*cc(l,i_ho2)
     $           + j(l,j_ho2))
     $          /(c011(l)*cc(l,i_o2))
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        oh/ho2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         roh_ho2 = (c001(l)*cc(l,i_o)
     $            + c003(l)*cc(l,i_o3)*rh_ho2
     $            + 2.*c004(l)*cc(l,i_h)
     $            + 2.*c008(l)*cc(l,i_ho2)
     $            + c015(l)*cc(l,i_o3)
     $            + d003(l)*no(l)
     $            + j(l,j_ho2)
     $            + j(l,j_h2o)*cc(l,i_h2o)
     $             /max(cc(l,i_ho2),dens(l)*1.e-30))
     $            /(c002(l)*cc(l,i_o)
     $            + c007h(l)*cc(l,i_ho2)
     $            + e001(l)*cc(l,i_co))
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        h
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         cc(l,i_h) = cc(l,i_hox)
     $                 /(1. + (1. + roh_ho2)/rh_ho2)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        ho2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         cc(l,i_ho2) = cc(l,i_h)/rh_ho2
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        oh
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         cc(l,i_oh) = cc(l,i_ho2)*roh_ho2
c
      end do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     modification of c007
c        by heterogenous reactions h001 and h002
c           after update of cc(l,i_oh) and cc(l,i_ho2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (HETERO.eq.1) then
       do l = 1,lswitch-1
         c007h(l) = c007(l)
     $            + h001(l)/max(cc(l,i_oh), 1.e-30*dens(l))
     $            + h002(l)/max(cc(l,i_ho2),1.e-30*dens(l))
c         write(*,*) "C007 = ",c007(l)," / C007H = ",c007h(l)
       end do
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        ox species
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        day:
c        - o1d at photochemical equilibrium 
c        - o3 at photochemical equilibrium
c        - continuity equation for ox
c        night:
c        - o1d = 0
c        - continuity equation for o3
c        - continuity equation for o
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        
      if (sza .le. 95.) then
c
         do l = 1,lswitch-1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        o(1d)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
            cc(l,i_o1d) = (j(l,j_co2_o1d)*co2(l)
     $                   + j(l,j_o2_o1d)*cc(l,i_o2)
     $                   + j(l,j_o3_o1d)*cc(l,i_o3))
     $                   /( b001(l)*co2(l)
     $                   + b002(l)*cc(l,i_h2o)
     $                   + b003(l)*cc(l,i_h2)
     $                   + b004(l)*cc(l,i_o2)
     $                   + b005(l)*cc(l,i_o3)
     $                   + b006(l)*cc(l,i_o3))
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        o/o3
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
            ro_o3 = (j(l,j_o3_o1d) + j(l,j_o3_o)
     $              + a003(l)*cc(l,i_o)
     $              + c003(l)*cc(l,i_h)
     $              + c014(l)*cc(l,i_oh)
     $              + c015(l)*cc(l,i_ho2))
     $              /(a001(l)*cc(l,i_o2))
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        o3
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
            cc(l,i_o3) = cc(l,i_ox)/(1. + ro_o3)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        o
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
            cc(l,i_o) = cc(l,i_o3)*ro_o3
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        ox = o + o3
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
            production(l,i_ox) =
     $                   + j(l,j_co2_o)*co2(l)
     $                   + j(l,j_co2_o1d)*co2(l)
     $                   + j(l,j_ho2)*cc(l,i_ho2)
     $                   + 2.*j(l,j_o2_o)*cc(l,i_o2)
     $                   + 2.*j(l,j_o2_o1d)*cc(l,i_o2)
     $                   + c006(l)*cc(l,i_h)*cc(l,i_ho2)
     $                   + c013(l)*cc(l,i_oh)*cc(l,i_oh)
     $                   + d003(l)*cc(l,i_ho2)*no(l)
c
            loss(l,i_ox) = 2.*a002(l)*cc(l,i_o)*cc(l,i_o)
     $                   + 2.*a003(l)*cc(l,i_o)*cc(l,i_o3)
     $                   + c001(l)*cc(l,i_ho2)*cc(l,i_o)
     $                   + c002(l)*cc(l,i_oh)*cc(l,i_o)
     $                   + c003(l)*cc(l,i_h)*cc(l,i_o3)
     $                   + c012(l)*cc(l,i_o)*cc(l,i_h2o2)
     $                   + c014(l)*cc(l,i_o3)*cc(l,i_oh)
     $                   + c015(l)*cc(l,i_o3)*cc(l,i_ho2)
     $                   + d001(l)*cc(l,i_o)*no2(l)
     $                   + e002(l)*cc(l,i_o)*cc(l,i_co)
c
            loss(l,i_ox) = loss(l,i_ox)/cc(l,i_ox)
c
         end do
c
      else
c
         do l = 1,lswitch-1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        o(1d)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
            cc(l,i_o1d) = 0.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        o3
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
            production(l,i_o3) = a001(l)*cc(l,i_o2)*cc(l,i_o)
c
            loss(l,i_o3) = a003(l)*cc(l,i_o)
     $                   + c003(l)*cc(l,i_h)
     $                   + c014(l)*cc(l,i_oh)
     $                   + c015(l)*cc(l,i_ho2)
c

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        o
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
            production(l,i_o) = c006(l)*cc(l,i_h)*cc(l,i_ho2)
     $                        + c013(l)*cc(l,i_oh)*cc(l,i_oh)
c
            loss(l,i_o)  =  a001(l)*cc(l,i_o2)
     $                   + 2.*a002(l)*cc(l,i_o)
     $                   + a003(l)*cc(l,i_o3)
     $                   + c001(l)*cc(l,i_ho2)
     $                   + c002(l)*cc(l,i_oh)
     $                   + c012(l)*cc(l,i_h2o2)
     $                   + e002(l)*cc(l,i_co)
c
         end do
      end if
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     other species
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do l = 1,lswitch-1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        co
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         production(l,i_co) = j(l,j_co2_o)*co2(l)
     $                      + j(l,j_co2_o1d)*co2(l)
c
         loss(l,i_co) = e001(l)*cc(l,i_oh)
     $                + e002(l)*cc(l,i_o)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        o2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         production(l,i_o2) = 
     $                  j(l,j_o3_o)*cc(l,i_o3)
     $                + j(l,j_o3_o1d)*cc(l,i_o3)
     $                + a002(l)*cc(l,i_o)*cc(l,i_o)
     $                + 2.*a003(l)*cc(l,i_o)*cc(l,i_o3)
     $                + 2.*b005(l)*cc(l,i_o1d)*cc(l,i_o3)
     $                + b006(l)*cc(l,i_o1d)*cc(l,i_o3)
     $                + c001(l)*cc(l,i_o)*cc(l,i_ho2)
     $                + c002(l)*cc(l,i_o)*cc(l,i_oh)
     $                + c003(l)*cc(l,i_h)*cc(l,i_o3)
     $                + c005(l)*cc(l,i_h)*cc(l,i_ho2)
     $                + c007h(l)*cc(l,i_oh)*cc(l,i_ho2)
     $                + c008(l)*cc(l,i_ho2)*cc(l,i_ho2)
     $                + c014(l)*cc(l,i_o3)*cc(l,i_oh)
     $                + 2.*c015(l)*cc(l,i_o3)*cc(l,i_ho2)
     $                + c016(l)*cc(l,i_ho2)*cc(l,i_ho2)
     $                + d001(l)*cc(l,i_o)*no2(l)
c
         loss(l,i_o2) = j(l,j_o2_o)
     $                + j(l,j_o2_o1d)
     $                + a001(l)*cc(l,i_o)
     $                + c011(l)*cc(l,i_h)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        h2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         production(l,i_h2) = c005(l)*cc(l,i_h)*cc(l,i_ho2)
     $                      + c018(l)*cc(l,i_h)*cc(l,i_h)
c
         loss(l,i_h2) = b003(l)*cc(l,i_o1d)
     $                + c010(l)*cc(l,i_oh)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        h2o
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         production(l,i_h2o) = 
     $                   c006(l)*cc(l,i_h)*cc(l,i_ho2)
     $                 + c007h(l)*cc(l,i_oh)*cc(l,i_ho2)
     $                 + c009(l)*cc(l,i_oh)*cc(l,i_h2o2)
     $                 + c010(l)*cc(l,i_oh)*cc(l,i_h2)
     $                 + c013(l)*cc(l,i_oh)*cc(l,i_oh)
c
         loss(l,i_h2o) = j(l,j_h2o)
     $                 + b002(l)*cc(l,i_o1d)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        h2o2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         production(l,i_h2o2) =
     $                    c008(l)*cc(l,i_ho2)*cc(l,i_ho2)
     $                  + c016(l)*cc(l,i_ho2)*cc(l,i_ho2)
     $                  + c017(l)*cc(l,i_oh)*cc(l,i_oh)
c
         loss(l,i_h2o2) = j(l,j_h2o2)
     $                  + c009(l)*cc(l,i_oh)
     $                  + c012(l)*cc(l,i_o)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        hox = h + oh + ho2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         production(l,i_hox) = 
     $                   2.*j(l,j_h2o)*cc(l,i_h2o)
     $                 + 2.*j(l,j_h2o2)*cc(l,i_h2o2)
     $                 + 2.*b002(l)*cc(l,i_o1d)*cc(l,i_h2o)
     $                 + 2.*b003(l)*cc(l,i_o1d)*cc(l,i_h2)
     $                 + 2.*c012(l)*cc(l,i_o)*cc(l,i_h2o2)
c
         loss(l,i_hox) = 2.*c005(l)*cc(l,i_h)*cc(l,i_ho2)
     $                 + 2.*c006(l)*cc(l,i_h)*cc(l,i_ho2)
     $                 + 2.*c007h(l)*cc(l,i_oh)*cc(l,i_ho2)
     $                 + 2.*c008(l)*cc(l,i_ho2)*cc(l,i_ho2)
     $                 + 2.*c013(l)*cc(l,i_oh)*cc(l,i_oh)
     $                 + 2.*c016(l)*cc(l,i_ho2)*cc(l,i_ho2)
     $                 + 2.*c017(l)*cc(l,i_oh)*cc(l,i_oh)
     $                 + 2.*c018(l)*cc(l,i_h)*cc(l,i_h)
c
         loss(l,i_hox) = loss(l,i_hox)/cc(l,i_hox)
c
      end do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     update number densities
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     long-lived species
c
      do l = 1,lswitch-1
         cc(l,i_co)  = (cc0(l,i_co) + production(l,i_co)*dt)
     $              /(1. + loss(l,i_co)*dt)
         cc(l,i_o2)  = (cc0(l,i_o2) + production(l,i_o2)*dt)
     $              /(1. + loss(l,i_o2)*dt)
         cc(l,i_h2)  = (cc0(l,i_h2) + production(l,i_h2)*dt)
     $              /(1. + loss(l,i_h2)*dt)
         cc(l,i_h2o2)= (cc0(l,i_h2o2) + production(l,i_h2o2)*dt)
     $              /(1. + loss(l,i_h2o2)*dt)
         cc(l,i_h2o) = (cc0(l,i_h2o) + production(l,i_h2o)*dt)
     $              /(1. + loss(l,i_h2o)*dt)
         cc(l,i_hox) = (cc0(l,i_hox) + production(l,i_hox)*dt)
     $              /(1. + loss(l,i_hox)*dt)
      end do
c
c     ox species
c
      if (sza .le. 95.) then
         do l = 1,lswitch-1
            cc(l,i_ox) = (cc0(l,i_ox) + production(l,i_ox)*dt)
     $                   /(1. + loss(l,i_ox)*dt)
         end do
      else
         do l = 1,lswitch-1
            cc(l,i_o)  = (cc0(l,i_o) + production(l,i_o)*dt)
     $                   /(1. + loss(l,i_o)*dt)
            cc(l,i_o3) = (cc0(l,i_o3) + production(l,i_o3)*dt)
     $                   /(1. + loss(l,i_o3)*dt)
            cc(l,i_ox) = cc(l,i_o) + cc(l,i_o3)
         end do
      end if 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     end of loop over iterations
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      end do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     density -> volume mixing ratio conversion
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do iesp = 1,nesp
         do l = 1,lswitch-1
            rm(l,iesp) = max(cc(l,iesp)/dens(l), 1.e-30)
         end do
      end do
c
      return
      end
c
c*****************************************************************
c
      subroutine phot(lswitch, press, temp, sza, dist_sol, rmo3, j)
c
c*****************************************************************
c
      implicit none
c
#include "dimensions.h"
#include "dimphys.h"
#include "chimiedata.h" 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     inputs:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       
      integer lswitch         ! interface level between chemistries
      real press(nlayermx)             ! pressure (hPa)
      real temp(nlayermx)              ! temperature (K)
      real sza                         ! solar zenith angle (deg)
      real dist_sol                    ! sun distance (AU) 
      real rmo3(nlayermx)              ! ozone mixing ratio
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     output:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real j(nlayermx,nd)              ! interpolated photolysis rates (s-1) 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     local:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      integer icol, ij, indsza, indcol, indozo, indtemp, 
     $        iozo, isza, it, iztop, l
c
      real col(nlayermx)               ! overhead air column   (molecule cm-2)
      real colo3(nlayermx)             ! overhead ozone column (molecule cm-2)
      real poids(2,2,2,2)              ! 4D interpolation weights 
      real tref                        ! temperature  at 1.9 hPa in the gcm (K)
      real table_temp(ntemp)           ! temperatures at 1.9 hPa in jmars   (K)
      real cinf, csup, cicol, 
     $     ciozo, cisza, citemp
      real colo3min, dp, scaleh
      real ratio_o3(nlayermx)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     day/night criterion
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (sza .le. 95.) then 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     temperatures at 1.9 hPa in jmars
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      table_temp(1) = 226.2
      table_temp(2) = 206.2
      table_temp(3) = 186.2
      table_temp(4) = 169.8
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     interpolation in solar zenith angle
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      indsza = nsza - 1
      do isza = 1,nsza
         if (szatab(isza) .ge. sza) then
            indsza = min(indsza,isza - 1)
            indsza = max(indsza, 1)
         end if
      end do
      cisza = (sza - szatab(indsza))
     $       /(szatab(indsza + 1) - szatab(indsza))
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     air and ozone columns
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     air column at model top: assign standard value
c
      scaleh = 10.
      iztop = nint(scaleh*log(press(1)/press(lswitch-1)))
      iztop = min(iztop,200)
      col(lswitch-1) = colairtab(iztop)
c
c     ozone column at model top
c
      colo3(lswitch-1) = 0.
c
c     air and ozone columns for other levels
c
      do l = lswitch-2,1,-1
         dp = (press(l) - press(l+1))*100.
         col(l) = col(l+1) + factor*dp
         col(l) = min(col(l), colairtab(0))
         colo3(l) = colo3(l+1)
     $            + factor*(rmo3(l+1) + rmo3(l))*0.5*dp
      end do
c
c     ratio ozone column/minimal theoretical column (ro3 = 7.171e-10)
c
      do l = 1,lswitch-1
         colo3min    = col(l)*7.171e-10
         ratio_o3(l) = colo3(l)/colo3min
         ratio_o3(l) = min(ratio_o3(l), table_ozo(nozo)*10.)
         ratio_o3(l) = max(ratio_o3(l), table_ozo(1)*10.)
      end do
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     temperature dependence
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     1) search for temperature at 1.9 hPa (tref): vertical interpolation
c
c      print*,'****** recherche de T a 1.9 hPa ******'
      tref = temp(1)
      do l = (lswitch-1)-1,1,-1
         if (press(l) .gt. 1.9) then
            cinf = (press(l) - 1.9)
     $            /(press(l) - press(l+1))
            csup = 1. - cinf
            tref = cinf*temp(l+1) + csup*temp(l)
c
c            print*, press(l+1), temp(l+1)
c            print*, press(l  ), temp(l  )
c            print*, 'tref = ', tref
            go to 10
         end if
      end do 
 10   continue
c
c     2) interpolation in temperature
c
      tref = min(tref,table_temp(1))
      tref = max(tref,table_temp(ntemp))
c
      do it = 2, ntemp
         if (table_temp(it) .le. tref) then
            citemp = (log(tref) - log(table_temp(it)))
     $              /(log(table_temp(it-1)) - log(table_temp(it)))
            indtemp = it - 1
c            print*,'interpolation entre ', table_temp(it), ' et ',
c     $              table_temp(it-1) 
c            print*,'indtemp = ', indtemp, ' citemp = ', citemp
            goto 20
         end if
      end do
  20  continue
c      print*,'****** fin interpolation en T ********'
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     loop over vertical levels
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do l = 1,lswitch-1
c
c     interpolation in air column
c
         do icol = 0,200
            if (colairtab(icol) .lt. col(l)) then
               cicol = (log(col(l)) - log(colairtab(icol)))
     $                /(log(colairtab(icol-1)) - log(colairtab(icol)))
               indcol = icol - 1
               goto 30
            end if
         end do
 30      continue
c
cc    interpolation in ozone column
c
         indozo = nozo - 1
         do iozo = 1,nozo
            if (table_ozo(iozo)*10. .ge. ratio_o3(l)) then
               indozo = min(indozo, iozo - 1)
               indozo = max(indozo, 1)
            end if
         end do
         ciozo = (ratio_o3(l) - table_ozo(indozo)*10.)
     $          /(table_ozo(indozo + 1)*10. - table_ozo(indozo)*10.)
c
c         write(6,'(i2,4e13.4,5f10.4)') 
c    $         l, press(l), col(l), 
c    $         colairtab(indcol), colairtab(indcol+1), 
c    $         colo3(l)/2.687e15, 
c    $         ratio_o3(l), table_ozo(indozo), table_ozo(indozo+1), 
c    $         ciozo
c 
cc    4-dimensional interpolation weights
c
         poids(1,1,1,1) = citemp*(1.-cisza)*cicol*(1.-ciozo)
         poids(1,1,1,2) = citemp*(1.-cisza)*cicol*ciozo
         poids(1,1,2,1) = citemp*(1.-cisza)*(1.-cicol)*(1.-ciozo)
         poids(1,1,2,2) = citemp*(1.-cisza)*(1.-cicol)*ciozo
         poids(1,2,1,1) = citemp*cisza*cicol*(1.-ciozo)
         poids(1,2,1,2) = citemp*cisza*cicol*ciozo
         poids(1,2,2,1) = citemp*cisza*(1.-cicol)*(1.-ciozo)
         poids(1,2,2,2) = citemp*cisza*(1.-cicol)*ciozo
         poids(2,1,1,1) = (1.-citemp)*(1.-cisza)*cicol*(1.-ciozo)
         poids(2,1,1,2) = (1.-citemp)*(1.-cisza)*cicol*ciozo
         poids(2,1,2,1) = (1.-citemp)*(1.-cisza)*(1.-cicol)*(1.-ciozo)
         poids(2,1,2,2) = (1.-citemp)*(1.-cisza)*(1.-cicol)*ciozo
         poids(2,2,1,1) = (1.-citemp)*cisza*cicol*(1.-ciozo)
         poids(2,2,1,2) = (1.-citemp)*cisza*cicol*ciozo
         poids(2,2,2,1) = (1.-citemp)*cisza*(1.-cicol)*(1.-ciozo)
         poids(2,2,2,2) = (1.-citemp)*cisza*(1.-cicol)*ciozo
c
cc    4-dimensional interpolation in the lookup table
c
         do ij = 1,nd
           j(l,ij) =
     $     poids(1,1,1,1)*jphot(indtemp,indsza,indcol,indozo,ij)
     $   + poids(1,1,1,2)*jphot(indtemp,indsza,indcol,indozo+1,ij)
     $   + poids(1,1,2,1)*jphot(indtemp,indsza,indcol+1,indozo,ij)
     $   + poids(1,1,2,2)*jphot(indtemp,indsza,indcol+1,indozo+1,ij)
     $   + poids(1,2,1,1)*jphot(indtemp,indsza+1,indcol,indozo,ij)
     $   + poids(1,2,1,2)*jphot(indtemp,indsza+1,indcol,indozo+1,ij)
     $   + poids(1,2,2,1)*jphot(indtemp,indsza+1,indcol+1,indozo,ij)
     $   + poids(1,2,2,2)*jphot(indtemp,indsza+1,indcol+1,indozo+1,ij)
     $   + poids(2,1,1,1)*jphot(indtemp+1,indsza,indcol,indozo,ij)
     $   + poids(2,1,1,2)*jphot(indtemp+1,indsza,indcol,indozo+1,ij)
     $   + poids(2,1,2,1)*jphot(indtemp+1,indsza,indcol+1,indozo,ij)
     $   + poids(2,1,2,2)*jphot(indtemp+1,indsza,indcol+1,indozo+1,ij)
     $   + poids(2,2,1,1)*jphot(indtemp+1,indsza+1,indcol,indozo,ij)
     $   + poids(2,2,1,2)*jphot(indtemp+1,indsza+1,indcol,indozo+1,ij)
     $   + poids(2,2,2,1)*jphot(indtemp+1,indsza+1,indcol+1,indozo,ij)
     $   + poids(2,2,2,2)*jphot(indtemp+1,indsza+1,indcol+1,indozo+1,ij)
         end do
c
cc    correction for sun distance
c
         do ij = 1,nd
            j(l,ij) = j(l,ij)*(1.52/dist_sol)**2.
         end do 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     end of loop over vertical levels
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      end do
c
      else
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     night
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do ij = 1,nd
         do l = 1,lswitch-1
            j(l,ij) = 0.
         end do
      end do
c
      end if
c
      return
      end
c
c*****************************************************************
c
      subroutine gcmtochim(zycol, lswitch, nesp, rm)
c
c*****************************************************************
c
      implicit none
c
#include "dimensions.h"
#include "dimphys.h"
#include "callkeys.h"
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     inputs:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      real zycol(nlayermx,nqmx)! species volume mixing ratio in the gcm
c
      integer nesp            ! number of species in the chemistry
      integer lswitch         ! interface level between chemistries
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     outputs:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real rm(nlayermx,nesp)   ! species volume mixing ratio 
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     local:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      integer l
      integer i_co2, i_co, i_o2, i_h2, i_h2o, i_h2o2, i_hox, i_ox,
     $        i_o1d, i_o, i_o3, i_h, i_oh, i_ho2, i_n2
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     tracers numbering in the gcm
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     co2      =  nqchem_min           
c     co       =  nqchem_min + 1          
c     o        =  nqchem_min + 2
c     o(1d)    =  nqchem_min + 3
c     o2       =  nqchem_min + 4
c     o3       =  nqchem_min + 5
c     h        =  nqchem_min + 6
c     h2       =  nqchem_min + 7
c     oh       =  nqchem_min + 8
c     ho2      =  nqchem_min + 9
c     h2o2     =  nqchem_min + 10
c     n2       =  nqchem_min + 11
c     ar       =  nqchem_min + 12
c     h2o      =  nqmx
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     tracers numbering in the chemistry
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      i_co2    =  1
      i_co     =  2
      i_o      =  3
      i_o1d    =  4
      i_o2     =  5
      i_o3     =  6
      i_h      =  7
      i_h2     =  8
      i_oh     =  9
      i_ho2    = 10
      i_h2o2   = 11
      i_h2o    = 12
      i_n2     = 13
      i_hox    = 14
      i_ox     = 15
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     initialise chemical species
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do l = 1,lswitch-1
       rm(l,i_co2)  = max(zycol(l, nqchem_min),      1.e-30)
       rm(l,i_co)   = max(zycol(l, nqchem_min + 1),  1.e-30)
       rm(l,i_o)    = max(zycol(l, nqchem_min + 2),  1.e-30)
       rm(l,i_o1d)  = max(zycol(l, nqchem_min + 3),  1.e-30)
       rm(l,i_o2)   = max(zycol(l, nqchem_min + 4),  1.e-30)
       rm(l,i_o3)   = max(zycol(l, nqchem_min + 5),  1.e-30)
       rm(l,i_h)    = max(zycol(l, nqchem_min + 6),  1.e-30)
       rm(l,i_h2)   = max(zycol(l, nqchem_min + 7),  1.e-30)
       rm(l,i_oh)   = max(zycol(l, nqchem_min + 8),  1.e-30)
       rm(l,i_ho2)  = max(zycol(l, nqchem_min + 9),  1.e-30)
       rm(l,i_h2o2) = max(zycol(l, nqchem_min + 10), 1.e-30)
       rm(l,i_n2)   = max(zycol(l, nqchem_min + 11), 1.e-30)
       rm(l,i_h2o)  = max(zycol(l, nqmx),            1.e-30)
c 
      end do 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     initialise chemical families                     c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do l = 1,lswitch-1
         rm(l,i_hox) = rm(l,i_h) 
     $               + rm(l,i_oh)
     $               + rm(l,i_ho2)
         rm(l,i_ox)  = rm(l,i_o) 
     $               + rm(l,i_o3)
      end do
c
      return
      end
c
c*****************************************************************
c
      subroutine chimtogcm(zycol, lswitch, nesp, rm)
c
c*****************************************************************
c
      implicit none
c
#include "dimensions.h"
#include "dimphys.h"
#include "callkeys.h"
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     inputs:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      integer nesp               ! number of species in the chemistry
      integer lswitch            ! interface level between chemistries
c      
      real rm(nlayermx,nesp)     ! species volume mixing ratio 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     output:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      real zycol(nlayermx,nqmx)  ! species volume mixing ratio in the gcm
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     local:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      integer l
      integer i_co2, i_co, i_o2, i_h2, i_h2o, i_h2o2, i_hox, i_ox,
     $        i_o1d, i_o, i_o3, i_h, i_oh, i_ho2, i_n2
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     tracers numbering in the gcm
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     co2      =  nqchem_min
c     co       =  nqchem_min + 1
c     o        =  nqchem_min + 2
c     o(1d)    =  nqchem_min + 3
c     o2       =  nqchem_min + 4
c     o3       =  nqchem_min + 5
c     h        =  nqchem_min + 6
c     h2       =  nqchem_min + 7
c     oh       =  nqchem_min + 8
c     ho2      =  nqchem_min + 9
c     h2o2     =  nqchem_min + 10
c     n2       =  nqchem_min + 11
c     ar       =  nqchem_min + 12
c     h2o      =  nqmx
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     tracers numbering in the chemistry
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      i_co2    =  1
      i_co     =  2
      i_o      =  3
      i_o1d    =  4
      i_o2     =  5
      i_o3     =  6
      i_h      =  7
      i_h2     =  8
      i_oh     =  9
      i_ho2    = 10
      i_h2o2   = 11
      i_h2o    = 12
      i_n2     = 13
      i_hox    = 14
      i_ox     = 15
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     save mixing ratios for the gcm
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do l = 1,lswitch-1
         zycol(l, nqchem_min)      = rm(l,i_co2) 
         zycol(l, nqchem_min + 1)  = rm(l,i_co) 
         zycol(l, nqchem_min + 2)  = rm(l,i_o) 
         zycol(l, nqchem_min + 3)  = rm(l,i_o1d)
         zycol(l, nqchem_min + 4)  = rm(l,i_o2) 
         zycol(l, nqchem_min + 5)  = rm(l,i_o3) 
         zycol(l, nqchem_min + 6)  = rm(l,i_h)  
         zycol(l, nqchem_min + 7)  = rm(l,i_h2) 
         zycol(l, nqchem_min + 8)  = rm(l,i_oh) 
         zycol(l, nqchem_min + 9)  = rm(l,i_ho2) 
         zycol(l, nqchem_min + 10) = rm(l,i_h2o2)
         zycol(l, nqchem_min + 11) = rm(l,i_n2)
         zycol(l, nqmx)            = rm(l,i_h2o)
      end do 
c
c  water ice: zycol(nqmx-1) has not been changed if iceparty on
c  
      return
      end
c
c*****************************************************************
c
      subroutine chemrates(lswitch, dens, press, t, icesurf,
     $                     a001, a002, a003,
     $                     b001, b002, b003, b004, b005, b006,
     $                     c001, c002, c003, c004, c005, c006,
     $                     c007, c008, c009, c010, c011, c012,
     $                     c013, c014, c015, c016, c017, c018,
     $                     d001, d002, d003,
     $                     e001, e002,
     $                     h001, h002)
c
c*****************************************************************
c
      implicit none
c
#include "dimensions.h"
#include "dimphys.h"
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     inputs:                                                    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      integer lswitch                  ! interface level between chemistries

      real dens(nlayermx)              ! density (cm-3)
      real press(nlayermx)             ! pressure (hpa)
      real t(nlayermx)                 ! temperature (k)
      real icesurf(nlayermx)           ! ice surface area (cm^2/cm^3)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     outputs:                                                   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real a001(nlayermx), a002(nlayermx), a003(nlayermx)
      real b001(nlayermx), b002(nlayermx), b003(nlayermx), 
     $     b004(nlayermx), b005(nlayermx), b006(nlayermx)
      real c001(nlayermx), c002(nlayermx), c003(nlayermx), 
     $     c004(nlayermx), c005(nlayermx), c006(nlayermx), 
     $     c007(nlayermx), c008(nlayermx), c009(nlayermx), 
     $     c010(nlayermx), c011(nlayermx), c012(nlayermx),
     $     c013(nlayermx), c014(nlayermx), c015(nlayermx), 
     $     c016(nlayermx), c017(nlayermx), c018(nlayermx)
      real d001(nlayermx), d002(nlayermx), d003(nlayermx)
      real e001(nlayermx), e002(nlayermx)
      real h001(nlayermx), h002(nlayermx)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     local:                                                     c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real ak0, ak1, rate, xpo
c
      integer l
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     compute reaction rates from JPL 1997, otherwise mentioned
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do l = 1,lswitch-1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        oxygen compounds
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccc      a001: o + o2 + co2 -> o3 + co2
c
c        jpl 2003
c
         a001(l) = 2.5
     $             *6.0e-34*(t(l)/300.)**(-2.4)*dens(l)
c
c        nair et al., 1994
c
c        a001(l) = 1.3e-34*exp(724./t(l))*dens(l)
c
ccc      a002: o + o + co2 -> o2 + co2
c
c        Tsang and Hampson, J. Chem. Phys. Ref. Data, 15, 1087, 1986
c
         a002(l) = 2.5*5.2e-35*exp(900./t(l))*dens(l)
c
c        Campbell and Gray, Chem. Phys. Lett., 18, 607, 1973
c
c        a002(l) = 1.2e-32*(300./t(l))**(2.0)*dens(l)
c
ccc      a003: o + o3 -> o2 + o2
c
c        jpl 2003
c
         a003(l) = 8.0e-12*exp(-2060./t(l))
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        reactions with o(1d)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccc      b001: o(1d) + co2  -> o + co2
c
c        jpl 2003
c
         b001(l) = 7.4e-11*exp(120./t(l))
c
ccc      b002: o(1d) + h2o  -> oh + oh
c
c        jpl 2003
c
         b002(l) = 2.2e-10
c    
ccc      b003: o(1d) + h2  -> oh + h
c
c        jpl 2003
c
         b003(l) = 1.1e-10
c
c        nair et al., 1994
c
c        b003(l) = 1.0e-10
c    
ccc      b004: o(1d) + o2  -> o + o2
c
c        jpl 2003
c
         b004(l) = 3.2e-11*exp(70./t(l))
c    
ccc      b005: o(1d) + o3  -> o2 + o2
c
c        jpl 2003
c
         b005(l) = 1.2e-10
c    
ccc      b006: o(1d) + o3  -> o2 + o + o
c
c        jpl 2003
c
         b006(l) = 1.2e-10
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        hydrogen compounds    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccc      c001: o + ho2 -> oh + o2
c
c        jpl 2003
c 
         c001(l) = 3.0e-11*exp(200./t(l))
c        c001(l) = 0.75*3.0e-11*exp(200./t(l))
c
ccc      c002: o + oh -> o2 + h
c
c        jpl 2003
c 
         c002(l) = 2.2e-11*exp(120./t(l))
c
ccc      c003: h + o3 -> oh + o2
c
c        jpl 2003
c 
         c003(l) = 1.4e-10*exp(-470./t(l))
c
ccc      c004: h + ho2 -> oh + oh
c
c        jpl 2003
c 
         c004(l) = 8.1e-11*0.90
c
ccc      c005: h + ho2 -> h2 + o2
c
c        jpl 2003
c 
         c005(l) = 8.1e-11*0.08
c
ccc      c006: h + ho2 -> h2o + o
c
c        jpl 2003
c 
         c006(l) = 8.1e-11*0.02
c
ccc      c007: oh + ho2 -> h2o + o2
c
c        jpl 2003
c 
         c007(l) = 4.8e-11*exp(250./t(l))
c        c007(l) = 0.75*4.8e-11*exp(250./t(l))
c
ccc      c008: ho2 + ho2 -> h2o2 + o2
c
c        jpl 2003
c 
c        c008(l) = 2.3e-13*exp(600./t(l))
c
c        christensen et al., grl, 13, 2002
c
         c008(l) = 1.5e-12*exp(19./t(l))
c
ccc      c009: oh + h2o2 -> h2o + ho2
c
c        jpl 2003
c 
         c009(l) = 2.9e-12*exp(-160./t(l))
c
ccc      c010: oh + h2 -> h2o + h
c
c        jpl 2003
c 
         c010(l) = 5.5e-12*exp(-2000./t(l))
c
ccc      c011: h + o2 + co2 -> ho2 + co2
c
c        jpl 2003
c 
         ak0 = 2.5*5.7e-32*(t(l)/300.)**(-1.6)
         ak1 = 7.5e-11*(t(l)/300.)**(0.0)
c
         rate = (ak0*dens(l))/(1. + ak0*dens(l)/ak1)
         xpo = 1./(1. + alog10((ak0*dens(l))/ak1)**2)
         c011(l) = rate*0.6**xpo
c
ccc      c012: o + h2o2 -> oh + ho2
c
c        jpl 2003
c 
         c012(l) = 1.4e-12*exp(-2000./t(l))
c
ccc      c013: oh + oh -> h2o + o
c
c        jpl 2003
c 
         c013(l) = 4.2e-12*exp(-240./t(l))
c
ccc      c014: oh + o3 -> ho2 + o2
c
c        jpl 2003
c
         c014(l) = 1.7e-12*exp(-940./t(l))
c
c        jpl 2000
c
c        c014(l) = 1.5e-12*exp(-880./t(l))
c
c        nair et al., 1994 (jpl 1997)
c
c        c014(l) = 1.6e-12*exp(-940./t(l))
c
ccc      c015: ho2 + o3 -> oh + o2 + o2
c
c        jpl 2003
c
         c015(l) = 1.0e-14*exp(-490./t(l))
c
c        jpl 2000
c
c        c015(l) = 2.0e-14*exp(-680./t(l))
c
c        nair et al., 1994 (jpl 1997)
c
c        c015(l) = 1.1e-14*exp(-500./t(l))
c
ccc      c016: ho2 + ho2 + co2 -> h2o2 + o2 + co2
c
c        jpl 2003
c
         c016(l) = 2.5*1.7e-33
     $              *exp(1000./t(l))*dens(l)
c
ccc      c017: oh + oh + co2 -> h2o2 + co2
c
c        jpl 2003
c
         ak0 = 2.5*6.9e-31*(t(l)/300.)**(-1.0)
         ak1 = 2.6e-11*(t(l)/300.)**(0.0)
c
c        jpl 1997
c
c        ak0 = 2.5*6.2e-31*(t(l)/300.)**(-1.0)
c        ak1 = 2.6e-11*(t(l)/300.)**(0.0)
c
c        nair et al., 1994
c
c        ak0 = 2.5*7.1e-31*(t(l)/300.)**(-0.8)
c        ak1 = 1.5e-11*(t(l)/300.)**(0.0)
c
         rate = (ak0*dens(l))/(1. + ak0*dens(l)/ak1)
         xpo = 1./(1. + alog10((ak0*dens(l))/ak1)**2)
         c017(l) = rate*0.6**xpo
c
ccc      c018: h + h + co2 -> h2 + co2
c
c        baulch et al., 1992
c
         c018(l) = 2.5*8.85e-33*(t(l)/298.)**(-0.6)*dens(l)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        nitrogen compounds
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccc      d001: no2 + o -> no + o2
c
c        jpl 2003
c
         d001(l) = 5.6e-12*exp(180./t(l))
c
ccc      d002: no + o3 -> no2 + o2
c
c        jpl 2003
c
         d002(l) = 3.0e-12*exp(-1500./t(l))
c
ccc      d003: no + ho2 -> no2 + oh
c
c        jpl 2003
c
         d003(l) = 3.5e-12*exp(250./t(l))
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        carbon compounds
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccc      e001: oh + co -> co2 + h
c
c        jpl 2003
c
c        e001(l) = 1.5e-13*(1 + 0.6*press(l)/1013.)
c
c        mccabe et al., grl, 28, 3135, 2001
c
         e001(l) = 1.57e-13 + 3.54e-33*dens(l)
c
ccc      e002: o + co + m -> co2 + m
c
c        tsang and hampson, 1986.
c
         e002(l) = 2.5*6.5e-33*exp(-2184./t(l))*dens(l)
c
c        baulch et al., butterworths, 1976.
c
c        e002(l) = 1.6e-32*exp(-2184./t(l))*dens(l)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        heterogenous chemistry
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccc      h001: ho2 + ice surface -> (ho2)ads
ccc           (ho2)ads + oh -> h2o + o2
ccc      total : oh + ho2 -> h2o + o2, first order
ccc         K = h001*[ho2],
ccc      h001 = (S*v*gamma1)/4. (s-1), v=100*sqrt(3kT/m_ho2) (cm s-1)
c
c        cooper and abbatt, 1996
c
         h001(l) = icesurf(l) * 25.*sqrt(3*8.31*t(l)/33.) * 0.025
c         h001(l) = 0.0e0
c
ccc      h002: oh + ice surface -> (oh)ads
ccc           (oh)ads + ho2 -> h2o + o2
ccc      total : oh + ho2 -> h2o + o2, first order
ccc         K = h002*[oh],
ccc      h002 = (S*v*gamma2)/4. (s-1), v=100*sqrt(3kT/m_oh) (cm s-1)
c
c        cooper and abbatt, 1996
c
         h002(l) = icesurf(l) * 25.*sqrt(3*8.31*t(l)/17.) * 0.03
c         h002(l) = 0.0e0
c
c        write(*,*) "level ",l," / icesurf=",icesurf(l),
c     $       " / 0.25vg1=",(25.*sqrt(3*8.31*t(l)/33.) * 0.025),
c     $       " / 0.25vg2=",(25.*sqrt(3*8.31*t(l)/17.) * 0.03)
c
      end do
c
      return
      end
