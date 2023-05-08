      SUBROUTINE concentrations(ngrid,nlayer,nq,
     &                          pplay,pt,pdt,pq,pdq,ptimestep)
                                             
      use tracer_h, only: igcm_co2, igcm_co, igcm_o, igcm_o1d,
     &                      igcm_o2, igcm_o3, igcm_h, igcm_h2,
     &                      igcm_oh, igcm_ho2, igcm_n2, igcm_ar,
     &                      igcm_h2o_vap, igcm_n, igcm_no, igcm_no2,
     &                      igcm_n2d, igcm_ch4, 
     &                      igcm_ch3, igcm_ch, igcm_3ch2, igcm_1ch2,      
     &                      igcm_cho, igcm_ch2o, igcm_ch3o,               
     &                      igcm_c, igcm_c2, igcm_c2h, igcm_c2h2,         
     &                      igcm_c2h3, igcm_c2h4, igcm_c2h6, igcm_ch2co,  
     &                      igcm_ch3co, igcm_hcaer,
     &                      igcm_h2o2, mmol

      use conc_mod, only: mmean, Akknew, rnew, cpnew
      USE comcstfi_mod                   
      use callkeys_mod
      implicit none

!=======================================================================
! CALCULATION OF MEAN MOLECULAR MASS, Cp, Akk and R
!
! mmean(ngrid,nlayer)	amu
! cpnew(ngrid,nlayer)	J/kg/K
! rnew(ngrid,nlayer)	J/kg/K
! akknew(ngrid,nlayer)	coefficient of thermal concduction
!
! version: April 2012 - Franck Lefevre
!=======================================================================

!     declarations
 
#include "chimiedata.h"
!     input/output

      integer,intent(in) :: ngrid ! number of atmospheric columns
      integer,intent(in) :: nlayer ! number of atmospheric layers
      integer,intent(in) :: nq ! number of tracers
      real,intent(in) :: pplay(ngrid,nlayer)
      real,intent(in) :: pt(ngrid,nlayer)
      real,intent(in) :: pdt(ngrid,nlayer)
      real,intent(in) :: pq(ngrid,nlayer,nq)
      real,intent(in) :: pdq(ngrid,nlayer,nq)
      real,intent(in) :: ptimestep

!     local variables

      integer       :: i, l, ig, iq
      integer, save :: nbq
      integer,allocatable,save :: niq(:)
      real          :: ni(nq), ntot
      real          :: zq(ngrid, nlayer, nq)
      real          :: zt(ngrid, nlayer)
      real,allocatable,save    :: aki(:)
      real,allocatable,save    :: cpi(:)

      logical, save :: firstcall = .true.


      if (firstcall) then

         ! allocate local saved arrays:
         allocate(aki(nq))
         allocate(cpi(nq))
         allocate(niq(nq))
!        find index of chemical tracers to use
!        initialize thermal conductivity and specific heat coefficients
!        !? values are estimated

         nbq = 0 ! to count number of tracers used in this subroutine

         if (igcm_co2 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_co2
            aki(nbq) = 3.072e-4
            cpi(nbq) = 0.834e3
         end if
         if (igcm_co /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_co
            aki(nbq) = 4.87e-4
            cpi(nbq) = 1.034e3
         end if
         if (igcm_o /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_o
            aki(nbq) = 7.59e-4
            cpi(nbq) = 1.3e3
         end if
         if (igcm_o1d /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_o1d
            aki(nbq) = 7.59e-4  !?
            cpi(nbq) = 1.3e3    !?
         end if
         if (igcm_o2 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_o2
            aki(nbq) = 5.68e-4
            cpi(nbq) = 0.9194e3
         end if
         if (igcm_o3 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_o3
            aki(nbq) = 3.00e-4  !?
            cpi(nbq) = 0.800e3  !?
         end if
         if (igcm_h /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_h
            aki(nbq) = 0.0
            cpi(nbq) = 20.780e3
         end if
         if (igcm_h2 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_h2
            aki(nbq) = 36.314e-4
            cpi(nbq) = 14.266e3
         end if
         if (igcm_oh /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_oh
            aki(nbq)  = 7.00e-4 !?
            cpi(nbq)  = 1.045e3
         end if
         if (igcm_ho2 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_ho2
            aki(nbq) = 0.0
            cpi(nbq) = 1.065e3  !?
         end if
         if (igcm_h2o2 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_h2o2
            aki(nbq) = 0.0
            cpi(nbq) = 1.065e3  !?
         end if
         if (igcm_h2o_vap /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_h2o_vap
            aki(nbq) = 0.0
            cpi(nbq) = 1.870e3
         end if
         if (igcm_n /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_n
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_n2d /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_n2d
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_no /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_no
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_no2 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_no2
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if (igcm_n2 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_n2
            aki(nbq) = 5.6e-4
            cpi(nbq) = 1.034e3
         end if
         if(igcm_ch4 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_ch4
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_ch3 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_ch3
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_ch /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_ch
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_1ch2 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_1ch2
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_3ch2 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_3ch2
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_cho /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_cho
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_ch2o /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_ch2o
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_ch3o /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_ch3o
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_c /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_c
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_c2 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_c2
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_c2h /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_c2h
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_c2h2 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_c2h2
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_c2h3 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_c2h3
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_c2h4 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_c2h4
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_c2h6 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_c2h6
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_ch2co /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_ch2co
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_ch3co /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_ch3co
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_hcaer /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_hcaer
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if (igcm_ar /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_ar
            aki(nbq) = 0.0      !?
            cpi(nbq) = 1.000e3  !?
         end if


         ! tell the world about it:
         write(*,*) "concentrations: firstcall, nbq=",nbq
         write(*,*) "  niq(1:nbq)=",niq(1:nbq)
         write(*,*) "  aki(1:nbq)=",aki(1:nbq)
         write(*,*) "  cpi(1:nbq)=",cpi(1:nbq)

         firstcall = .false.

      end if ! if (firstcall)

!     update temperature

      do l = 1,nlayer
         do ig = 1,ngrid
            zt(ig,l) = pt(ig,l) + pdt(ig,l)*ptimestep
         end do
      end do

!     update tracers

      do l = 1,nlayer
         do ig = 1,ngrid
            do i = 1,nbq
               iq = niq(i) 
               zq(ig,l,iq) = max(1.e-30, pq(ig,l,iq)
     $                                 + pdq(ig,l,iq)*ptimestep)
            end do
         end do
      end do

!     mmean : mean molecular mass
!     rnew  : specific gas constant

      mmean(:,:)  = 0.
      do l = 1,nlayer
         do ig = 1,ngrid
            do i = 1,nbq
               iq = niq(i)
               mmean(ig,l) = mmean(ig,l) + zq(ig,l,iq)/mmol(iq)
            end do
            mmean(ig,l) = 1./mmean(ig,l)
            rnew(ig,l) = 8.314/mmean(ig,l)*1.e3     ! J/kg/K		
         end do
      end do


!     cpnew  : specicic heat
!     akknew : thermal conductivity cofficient
      cpnew(:,:)  = 0.
      akknew(:,:) = 0.

      do l = 1,nlayer
         do ig = 1,ngrid
            ntot = pplay(ig,l)/(1.381e-23*zt(ig,l))*1.e-6  ! in #/cm3
            do i = 1,nbq
               iq = niq(i) 
               ni(iq) = ntot*zq(ig,l,iq)*mmean(ig,l)/mmol(iq)
               cpnew(ig,l) = cpnew(ig,l) + ni(iq)*cpi(i)
               akknew(ig,l) = akknew(ig,l) + ni(iq)*aki(i)
            end do 
            cpnew(ig,l) = cpnew(ig,l)/ntot
            akknew(ig,l) = akknew(ig,l)/ntot
         end do
      end do

      return
      end 
