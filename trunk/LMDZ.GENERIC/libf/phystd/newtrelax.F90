subroutine newtrelax(ngrid,nlayer,mu0,sinlat,popsk,temp,pplay,pplev,dtrad,firstcall) 
        
  use comcstfi_mod, only: rcp, pi
  use callkeys_mod, only: tau_relax
  implicit none

#include "netcdf.inc"

!==================================================================
!     
!     Purpose
!     -------
!     Alternative Newtonian radiative transfer scheme.
!     
!     Authors
!     -------
!     R. Wordsworth (2010)
!     
!==================================================================
 
 
  ! Input
  integer,intent(in) :: ngrid, nlayer
  logical,intent(in) :: firstcall
  real,intent(in) :: mu0(ngrid)            ! cosine of sun incident angle
  real,intent(in) :: sinlat(ngrid)         ! sine of latitude
  real,intent(in) :: temp(ngrid,nlayer)    ! temperature at each layer (K)
  real,intent(in) :: pplay(ngrid,nlayer)   ! pressure at each layer (Pa)
  real,intent(in) :: pplev(ngrid,nlayer+1) ! pressure at each level (Pa)
  real,intent(in) :: popsk(ngrid,nlayer)   ! pot. T to T converter

  ! Output
  real,intent(out) :: dtrad(ngrid,nlayer) 

  ! Internal
  real Trelax_V, Trelax_H
  real,allocatable,dimension(:,:),save :: Trelax
!$OMP THREADPRIVATE(Trelax)

  real T_trop ! relaxation temperature at tropopause (K)
  real T_surf ! relaxation temperature at surface (K)
  real dT_EP  ! Equator-Pole relaxation temperature difference (K)

  real sig, f_sig, sig_trop
  integer l,ig


  logical tidallocked
  parameter (tidallocked = .true.)

  ! Setup relaxation temperature  
  if(firstcall) then

     ALLOCATE(Trelax(ngrid,nlayer))

     print*,'-----------------------------------------------------'
     print*,'| ATTENTION: You are using a Newtonian cooling scheme'
     print*,'| for the radiative transfer. This means that ALL'
     print*,'| other physics subroutines must be switched off.'
     print*,'-----------------------------------------------------'

     if(tidallocked)then
        do ig=1,ngrid

           T_surf = 126. + 239.*mu0(ig)
           T_trop = 140. + 52.*mu0(ig)
           do l=1,nlayer

              if(mu0(ig).le.0.0)then ! night side
                 Trelax(ig,l)=0.0
              else                   ! day side
                 Trelax(ig,l) = T_surf*popsk(ig,l)
                 if (Trelax(ig,l).lt.T_trop) Trelax(ig,l) = T_trop
              endif

           enddo
        enddo

     else

        T_trop = 200.
        T_surf = 288.
        dT_EP  = 70.

        sig_trop=(T_trop/T_surf)**(1./rcp)

        do l=1,nlayer
           do ig=1,ngrid

              ! vertically varying component
              Trelax_V = T_surf*popsk(ig,l)
              if (Trelax_V.lt.T_trop) Trelax_V = T_trop
              
              ! horizontally varying component
              sig = pplay(ig,l)/pplev(ig,1)
              if(sig.ge.sig_trop)then
                 f_sig=sin((pi/2)*((sig-sig_trop)/(1-sig_trop)))
              else
                 f_sig=0.0
              endif
              Trelax_H = -f_sig*dT_EP*(sinlat(ig)**2 - 1./3.)
              
              Trelax(ig,l) = Trelax_V + Trelax_H            
           
           enddo
        enddo

     endif

  endif

  ! Calculate radiative forcing
  do l=1,nlayer
     do ig=1,ngrid
        dtrad(ig,l) = -(temp(ig,l) - Trelax(ig,l)) / tau_relax
        if(temp(ig,l).gt.500.)then ! Trelax(ig,l))then
           print*,'ig=',ig
           print*,'l=',l
           print*,'temp=',temp(ig,l)
           print*,'Trelax=',Trelax(ig,l)
        endif
     enddo
  enddo

  call writediagfi(ngrid,'Tref','rad forc temp','K',3,Trelax)
  !call writediagfi(ngrid,'ThetaZ','stellar zenith angle','deg',2,mu0)

end subroutine newtrelax
