










      subroutine totalcloudfrac(ngrid,nlayer,nq,rneb,totalrneb,pplev,pq,tau)

      use watercommon_h
      use comdiurn_h
      USE tracer_h, only: igcm_h2o_ice
      USE callkeys_mod, ONLY: CLFfixval
      implicit none

!==================================================================
!     
!     Purpose
!     -------
!     Calculates the total cloud fraction 
!     
!     Authors
!     -------
!     Adapted from the LMDTERRE code by B Charnay (2010)
!     
!==================================================================

      integer,intent(in) :: ngrid        ! number of atmospheric columns
      integer,intent(in) :: nlayer       ! number of atmospheric layers
      integer,intent(in) :: nq           ! number of tracers
      real,intent(in) :: rneb(ngrid,nlayer)    ! cloud fraction     
      real,intent(out) :: totalrneb(ngrid)       ! total cloud fraction 
      real,intent(in) :: pplev(ngrid,nlayer+1) ! inter-layer pressure (Pa)
      real,intent(in) :: pq(ngrid,nlayer,nq)   ! tracers (.../kg_of_air)
      real,intent(in) :: tau(ngrid,nlayer)

      real, dimension(nlayer+1) :: masse
      integer, parameter          :: recovery=7
      integer ltau_max
      real massetot

! hypothesis behind recovery. value: 
! 1 = random recovery
! 2 = maximal recovery
! 3 = minimal recovery
! 4 = fixed recovery
! 5 = recovery on the thicker layer
!     Local variables
      integer ig, l
      real clear,tau_min
      real, parameter ::  tau_c=0.1 !threshold of optical depth for the calculation of total cloud fraction 
      real rneb2(nlayer)


      do ig=1,ngrid
         totalrneb(ig) = 0.

         if (recovery.eq.1) then
            clear = (1.-rneb(ig,1))
            do l=2,nlayer      
               clear = clear*(1.-rneb(ig,l))
            enddo
            totalrneb(ig) = 1.-clear

         elseif (recovery.eq.2) then
            totalrneb(ig) = rneb(ig,1)
            do l=2,14 !nlayer    
               totalrneb(ig) = max(rneb(ig,l),totalrneb(ig))
            enddo
            
         elseif (recovery.eq.3) then
            totalrneb(ig) = rneb(ig,1)
            do l=2,nlayer    
               totalrneb(ig) = min(rneb(ig,l),totalrneb(ig))
            enddo
         
         elseif (recovery.eq.4) then
            totalrneb(ig) = CLFfixval

         elseif (recovery.eq.5) then
            totalrneb(ig) = rneb(ig,1)            
            do l=1,nlayer
               masse(l)=pq(ig,l,igcm_h2o_ice)*(pplev(ig,l)-pplev(ig,l+1))
            enddo
            ltau_max=maxloc(masse,dim=1)
            totalrneb(ig) = rneb(ig,ltau_max)

         elseif (recovery.eq.6) then
            totalrneb(ig) = 0.            
            do l=1,nlayer
               masse(l)=pq(ig,l,igcm_h2o_ice)*(pplev(ig,l)-pplev(ig,l+1))
               masse(l)=max(masse(l),0.)
            enddo
            massetot=sum(masse,dim=1)
            do l=1,nlayer
               totalrneb(ig) = totalrneb(ig)+rneb(ig,l)*masse(l)/massetot
            enddo

         elseif (recovery.eq.7) then

            rneb2(:)=rneb(ig,1:nlayer)
	    tau_min=MIN(tau_c,MAXVAL(tau(ig,1:nlayer))/2.)
            do l=1,nlayer
               if(tau(ig,l)<tau_min) rneb2(l)=0.      	
            enddo
            totalrneb(ig)=maxval(rneb2(1:nlayer))

         endif                  ! (recovery=)   

         totalrneb(ig) = min(1.,totalrneb(ig))
         totalrneb(ig) = max(0.,totalrneb(ig))
         
      enddo                     ! (ig=)
      
      
    end subroutine totalcloudfrac
