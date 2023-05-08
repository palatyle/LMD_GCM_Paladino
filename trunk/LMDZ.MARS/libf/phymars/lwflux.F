       subroutine lwflux (ig0,kdlon,kflev,dp
     .                   ,bsurf,btop,blev,blay,dbsublay
     .                   ,tlay, tlev, dt0      ! pour sortie dans g2d uniquement
     .                   ,emis
     .                   , tautotal,omegtotal,gtotal
     .                   ,coolrate,fluxground,fluxtop
     .                   ,netrad)

c----------------------------------------------------------------------
c     LWFLUX     computes the fluxes
c----------------------------------------------------------------------

      use dimradmars_mod, only: ndlo2, nir, ndlon, nuco2, nflev
      use yomlw_h, only: nlaylte, xi, xi_ground, gcp
      implicit none
 
      include "callkeys.h"
      include "comg1d.h"

c----------------------------------------------------------------------
c         0.1   arguments
c               ---------
c                                                            inputs:
c                                                            -------
      integer,intent(in) :: ig0
      integer,intent(in) :: kdlon           ! part of ngrid
      integer,intent(in) :: kflev           ! part of nlayer

      real,intent(in) :: dp (ndlo2,kflev)   ! layer pressure thickness (Pa)

      real,intent(in) :: bsurf (ndlo2,nir) ! surface spectral planck function
      real,intent(in) :: blev (ndlo2,nir,kflev+1) ! level   spectral planck function
      real,intent(in) :: blay (ndlo2,nir,kflev) ! layer   spectral planck function
      real,intent(in) :: btop (ndlo2,nir) ! top spectral planck function
      real,intent(in) :: dbsublay (ndlo2,nir,2*kflev)  ! layer gradient spectral planck
                                                       ! function in sub layers

      real,intent(in) :: dt0 (ndlo2) ! surface temperature discontinuity
      real,intent(in) :: tlay (ndlo2,kflev) ! layer temperature
      real,intent(in) :: tlev (ndlo2,kflev+1) ! level temperature

      real,intent(in) :: emis (ndlo2) ! surface emissivity

      real,intent(in) :: tautotal(ndlo2,kflev,nir)  ! \   Total single scattering
      real,intent(in) :: omegtotal(ndlo2,kflev,nir) !  >  properties (Addition of the
      real,intent(in) :: gtotal(ndlo2,kflev,nir)    ! /   NAERKIND aerosols prop.)


c                                                            outputs:
c                                                            --------
      real,intent(out) :: coolrate(ndlo2,kflev) ! radiative cooling rate (K/s) 
      real,intent(out) :: netrad (ndlo2,kflev) ! radiative budget (W/m2) 
      real,intent(out) :: fluxground(ndlo2) ! downward flux on the ground
                                            ! for surface radiative budget
      real,intent(out) :: fluxtop(ndlo2) ! upward flux on the top of atm ("OLR")


c----------------------------------------------------------------------
c         0.2   local arrays
c               ------------

      integer ja,jl,j,i,ig1d,ig,l
      real  ksidb (ndlon,nuco2+1,0:nflev+1,0:nflev+1) ! net exchange rate (W/m2)

      real dpsgcp (0:nflev+1,0:nflev+1)    !  dp/(g.cp)
      real temp (0:nflev+1,0:nflev+1)       

      real fluxdiff(ndlon,2,nflev+1)  ! diffusion flux: upward(1) downward(2) 

      real*4 reel4

c     To compute IR flux in the atmosphere  (For diagnostic only !!)
      logical computeflux
      real coefd(kdlon,nuco2,nflev+1,nflev+1) 
      real coefu(kdlon,nuco2,0:nflev,nflev+1)
      real flw_up(kdlon,nflev+1), flw_dn(kdlon,nflev+1) ! fluxes (W/m2)


      ksidb(:,:,:,:)=0

c----------------------------------------------------------------------
c         1.1   exchanges (layer i <--> all layers up to i) 
c               -------------------------------------------

      do i = 1,nlaylte
        do j = i+1,nlaylte
          do ja = 1,nuco2
            do jl = 1,kdlon

      ksidb(jl,ja,i,j) = xi(ig0+jl,ja,i,j)
     .                 * (blay(jl,ja,j)-blay(jl,ja,i))
c                                                        ksidb reciprocity
c                                                        -----------------
      ksidb(jl,ja,j,i) = -ksidb(jl,ja,i,j)
      
            enddo
          enddo
        enddo
      enddo

c----------------------------------------------------------------------
c         1.2   exchanges (ground <--> all layers) 
c               ----------------------------------

      do i = 1,nlaylte
        do ja = 1,nuco2
          do jl = 1,kdlon

      ksidb(jl,ja,i,0) = xi(ig0+jl,ja,0,i)
     .                 * (bsurf(jl,ja)-blay(jl,ja,i))
c                                                        ksidb reciprocity
c                                                        -----------------
      ksidb(jl,ja,0,i) = -ksidb(jl,ja,i,0)

          enddo
        enddo
      enddo

c--------------------------------------------------------
c  Here we add the neighbour contribution 
c           for exchanges between ground and first layer
c--------------------------------------------------------

        do ja = 1,nuco2
          do jl = 1,kdlon

      ksidb(jl,ja,1,0) = ksidb(jl,ja,1,0) 
     .                 - xi_ground(ig0+jl,ja)
     .                 * (blev(jl,ja,1)-blay(jl,ja,1))

cc                                                       ksidb reciprocity
cc                                                       -----------------
      ksidb(jl,ja,0,1) = - ksidb(jl,ja,1,0)

          enddo
        enddo

c----------------------------------------------------------------------
c         1.3   exchanges (layer i <--> space) 
c               ------------------------------

      do i = 1,nlaylte
        do ja = 1,nuco2
          do jl = 1,kdlon

      ksidb(jl,ja,i,nlaylte+1) = xi(ig0+jl,ja,i,nlaylte+1) 
     .                       * (-blay(jl,ja,i))
c                                                        ksidb reciprocity
c                                                        -----------------
      ksidb(jl,ja,nlaylte+1,i) = - ksidb(jl,ja,i,nlaylte+1)

          enddo
        enddo
      enddo

c----------------------------------------------------------------------
c         1.4   exchanges (ground <--> space) 
c               -----------------------------

      do ja = 1,nuco2
        do jl = 1,kdlon

      ksidb(jl,ja,0,nlaylte+1) = xi(ig0+jl,ja,0,nlaylte+1) 
     .                       * (-bsurf(jl,ja))

c                                                        ksidb reciprocity
c                                                        -----------------
      ksidb(jl,ja,nlaylte+1,0) = - ksidb(jl,ja,0,nlaylte+1)

        enddo
      enddo

c----------------------------------------------------------------------
c         2.0   sum of band 1 and 2 of co2 contribution 
c               ---------------------------------------

      do i = 0,nlaylte+1
        do j = 0,nlaylte+1
          do jl = 1,kdlon

            ksidb(jl,3,i,j)= ksidb(jl,1,i,j) + ksidb(jl,2,i,j)

          enddo
        enddo
      enddo

c----------------------------------------------------------------------
c         3.0   Diffusion
c               ---------

      i = nlaylte+1
      do jl = 1,kdlon
        fluxdiff(jl,1,i)   = 0.
        fluxdiff(jl,2,i)   = 0.
      enddo

      call lwdiff (kdlon,kflev
     .         ,bsurf,btop,dbsublay
     .         ,tautotal,omegtotal,gtotal
     .         ,emis,fluxdiff)

c----------------------------------------------------------------------
c         4.0   Radiative Budget for each layer i
c               ---------------------------------

      do i = 1,nlaylte
        do jl = 1,kdlon
            netrad(jl,i) = 0.
        enddo
      enddo

      do i = 1,nlaylte
        do j = 0,nlaylte+1
          do jl = 1,kdlon

            netrad(jl,i) = netrad(jl,i) + ksidb(jl,3,i,j)

          enddo
        enddo
      enddo
c                                                 diffusion contribution
c                                                 ----------------------
      do i = 1,nlaylte
        do jl = 1,kdlon

          netrad(jl,i) = netrad(jl,i) 
     .                 - fluxdiff(jl,1,i+1) - fluxdiff(jl,2,i+1)   
     .                 + fluxdiff(jl,1,i)   + fluxdiff(jl,2,i)

        enddo
      enddo

c----------------------------------------------------------------------
c         4.0   cooling rate for each layer i
c               -----------------------------

      do i = 1,nlaylte
        do jl = 1,kdlon

          coolrate(jl,i) = gcp * netrad(jl,i) / dp(jl,i)

        enddo
      enddo

c----------------------------------------------------------------------
c         5.0   downward flux (all layers --> ground): "fluxground" 
c               ---------------------------------------------------

      do jl = 1,kdlon
          fluxground(jl) = 0.
      enddo

      do i = 1,nlaylte
        do ja = 1,nuco2
          do jl = 1,kdlon

      fluxground(jl) = fluxground(jl)
     .               + xi(ig0+jl,ja,0,i) * (blay(jl,ja,i))

          enddo
        enddo
      enddo
   
      do jl = 1,kdlon
        fluxground(jl) = fluxground(jl) - fluxdiff(jl,2,1)  
      enddo

c----------------------------------------------------------------------
c         6.0   outgoing flux (all layers --> space): "fluxtop" 
c               ---------------------------------------------------

      do jl = 1,kdlon
          fluxtop(jl) = 0.
      enddo

      do i = 0,nlaylte
        do jl = 1,kdlon
           fluxtop(jl) = fluxtop(jl)- ksidb(jl,3,i,nlaylte+1)
        enddo
      enddo
   
      do jl = 1,kdlon
        fluxtop(jl) = fluxtop(jl) + fluxdiff(jl,1,nlaylte+1)  
      enddo

c----------------------------------------------------------------------
c         6.5 ONLY FOR DIAGNOSTIC : Compute IR flux in the atmosphere 
c             -------------------
c     The broadband fluxes (W.m-2) at every level  from surface level (l=1) 
c     up the top of the  upper layer (here: l=nlaylte+1) are:
c     upward : flw_up(ig1d,l)   ;   downward : flw_dn(ig1d,j) 
c
      computeflux = .false.

      IF (computeflux) THEN  ! not used by the GCM only for diagnostic !
c      upward flux 
c      ~~~~~~~~~~~
       do i = 0,nlaylte
        do j = 1,nlaylte+1
         do ja = 1,nuco2
          do jl = 1,kdlon
            coefu(jl,ja,i,j) =0.
            do l=j,nlaylte+1 
              coefu(jl,ja,i,j)=coefu(jl,ja,i,j)+xi(ig0+jl,ja,l,i)
            end do

          enddo
         enddo
        enddo
       enddo
       do j = 1,nlaylte+1
        do jl = 1,kdlon
           flw_up(jl,j) = 0.
           do ja = 1,nuco2
             flw_up(jl,j)=flw_up(jl,j)+bsurf(jl,ja)*coefu(jl,ja,0,j)
             do i=1,j-1
              flw_up(jl,j)=flw_up(jl,j)+blay(jl,ja,i)*coefu(jl,ja,i,j) 
             end do 
           end do
           flw_up(jl,j)=flw_up(jl,j) + fluxdiff(jl,1,j)
        end do
       end do

c      downward flux
c      ~~~~~~~~~~~~~
       do i = 1,nlaylte+1
        do j = 1,nlaylte+1
         do ja = 1,nuco2
          do jl = 1,kdlon
            coefd(jl,ja,i,j) =0.
            do l=0,j-1
              coefd(jl,ja,i,j)=coefd(jl,ja,i,j)+xi(ig0+jl,ja,l,i)
            end do
          enddo
         enddo
        enddo
       enddo
       do j = 1,nlaylte+1
        do jl = 1,kdlon
           flw_dn(jl,j) = 0.
           do ja = 1,nuco2
             do i=j,nlaylte
              flw_dn(jl,j)=flw_dn(jl,j)+blay(jl,ja,i)*coefd(jl,ja,i,j) 
             end do 
           end do
           flw_dn(jl,j)=flw_dn(jl,j) - fluxdiff(jl,2,j)
        end do
       end do
      END IF

c----------------------------------------------------------------------
c         7.0   outputs Grads 2D
c               ----------------

c ig1d: point de la grille physique ou on veut faire la sortie
c ig0+1:  point du decoupage de la grille physique

      if (callg2d) then

      ig1d = kdlon/2 + 1
c     ig1d = kdlon

      if ((ig0+1).LE.ig1d .and. ig1d.LE.(ig0+kdlon)
     .    .OR.  kdlon.EQ.1   ) then

          ig = ig1d-ig0
        print*, 'Sortie g2d: ig1d, ig, ig0', ig1d, ig, ig0

c--------------------------------------------
c   Ouverture de g2d.dat
c--------------------------------------------
      if (g2d_premier) then
        open (47,file='g2d.dat'
clmd     &       ,form='unformatted',access='direct',recl=4)
     &        ,form='unformatted',access='direct',recl=1
     &        ,status='unknown')
        g2d_irec=0
        g2d_appel=0
        g2d_premier=.false.
      endif
        g2d_appel = g2d_appel+1

c--------------------------------------------
c   Sortie g2d des xi proches + distants 
c--------------------------------------------
cl                               if (nflev .NE. 500) then
      do ja = 1,nuco2
        do j = 0,nlaylte+1
          do i = 0,nlaylte+1
            g2d_irec=g2d_irec+1
            reel4 = xi(ig1d,ja,i,j)
            write(47,rec=g2d_irec) reel4
          enddo
        enddo
      enddo
cl                               endif

c------------------------------------------------------
c   Writeg2d des ksidb 
c------------------------------------------------------
      do ja = 1,nuco2
c       ja=1
        do j = 0,nlaylte+1
          do i = 0,nlaylte+1
            g2d_irec=g2d_irec+1
            reel4 = ksidb(ig,ja,i,j)
            write(47,rec=g2d_irec) reel4
          enddo
        enddo
      enddo

      do j = 0,nlaylte+1
        do i = 0,nlaylte+1
          g2d_irec=g2d_irec+1
          reel4 = ksidb(ig,3,i,j)
          write(47,rec=g2d_irec) reel4
        enddo
      enddo

c------------------------------------------------------
c  Writeg2d dpsgcp
c------------------------------------------------------

        do j = 1 , nlaylte
          do i = 0 , nlaylte+1
            dpsgcp(i,j) = dp(ig,j) / gcp
          enddo
        enddo

        do i = 0 , nlaylte+1
c         dpsgcp(i,0) = 0.0002  ! (rapport ~ entre 1000 et 10000 pour le sol)
          dpsgcp(i,0) = 1.      ! (pour regler l'echelle des sorties)
          dpsgcp(i,nlaylte+1) = 0.
        enddo

c     print*
c     print*,'gcp: ',gcp
c     print*
c       do i = 0 , nlaylte+1
c     print*,i,'dp: ',dp(ig,i)
c       enddo
c     print*
c       do i = 0 , nlaylte+1
c     print*,i,'dpsgcp: ',dpsgcp(i,1)
c       enddo
  
      do j = 0,nlaylte+1
        do i = 0,nlaylte+1
          g2d_irec=g2d_irec+1
          reel4 = dpsgcp(i,j)
          write(47,rec=g2d_irec) reel4
        enddo
      enddo

c------------------------------------------------------
c  Writeg2d temperature
c------------------------------------------------------

        do j = 1 , nlaylte
          do i = 0 , nlaylte+1
            temp(i,j) = tlay(ig,j)
          enddo
        enddo

        do i = 0 , nlaylte+1
          temp(i,0) = tlev(ig,1)+dt0(ig)     ! temperature surface
          temp(i,nlaylte+1) = 0.               ! temperature espace  (=0)
        enddo

      do j = 0,nlaylte+1
        do i = 0,nlaylte+1
          g2d_irec=g2d_irec+1
          reel4 = temp(i,j)
          write(47,rec=g2d_irec) reel4
        enddo
      enddo

        write(76,*) 'ig1d, ig, ig0', ig1d, ig, ig0
        write(76,*) 'nlaylte', nlaylte
        write(76,*) 'nflev', nflev
        write(76,*) 'kdlon', kdlon
        write(76,*) 'ndlo2', ndlo2
        write(76,*) 'ndlon', ndlon
      do ja=1,4
        write(76,*) 'bsurf', ja, bsurf(ig,ja)
        write(76,*) 'btop', ja, btop(ig,ja)

        do j=1,nlaylte+1
          write(76,*) 'blev', ja, j, blev(ig,ja,j)
        enddo

        do j=1,nlaylte
          write(76,*) 'blay', ja, j, blay(ig,ja,j)
        enddo

        do j=1,2*nlaylte
          write(76,*) 'dbsublay', ja, j, dbsublay(ig,ja,j)
        enddo
      enddo

      endif
c************************************************************************
      endif  !   callg2d

      end
