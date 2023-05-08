c**************************************************************************
c
      subroutine nltecool(ngrid,nlayer,nq,pplay,pt,pq,dtnlte)
c
c This code was designed as a delivery for the "Martian Environment Models" 
c project ( ESA contract 11369/95/nl/jg CCN2 )
c Computes non-LTE heating rates from CO2 emission at 15 um
c in the Martian upper atmosphere.
c Uses a simplified model consisting of two excited levels with two
c emission bands, one of them stronger than the other, which correspond
c to the behaviours of the 626 fundamental band and the isotopic fund.bands.
c It uses a cool-to-space approximation with tabulated escape functions.
c These escape functions have been precomputed for the strong and weak bands,
c and are given as a function of pressure in separate files.
c The output values are the heating rates (actually, cooling, since they
c are always negative) for the two bands, i.e., the total cooling is the
c sum of them.
c Miguel A. Lopez Valverde
c Instituto de Astrofisica de Andalucia (CSIC), Granada, Spain
c
c Version 1b.  See description above.  22-March-2000.
c Adapted as a subroutine for use in GCM -- PLR/SRL 6/2000
c Version 1c.  Inclusion of VMR in the tabulation of escape functions. 
c              Table contains now only 1 input file -- Miguel 11/Jul/2000
c Version 1d  data contained in original input file "nlte_escape.dat"
c now stored in include file "nltedata.h" Y.Wanherdrick 09/2000

c       jul 2011 fgg   Modified to allow variable O
c     
c***************************************************************************

      use tracer_mod, only: igcm_co2, igcm_co, igcm_o, igcm_n2, mmol
      use conc_mod, only: mmean
      implicit none

#include "nltedata.h" ! (Equivalent to the reading of the "nlte_escape.dat" file)
#include "chimiedata.h"
#include "callkeys.h"

c Input and output variables
c
      integer     ngrid                           ! no. of horiz. gridpoints
      integer     nlayer                          ! no. of atmospheric layers
      integer     nq                              ! no. of tracers
      real        pplay(ngrid,nlayer)             ! input pressure grid
      real        pt(ngrid,nlayer)                ! input temperatures
      real        pq(ngrid,nlayer,nq)                ! input mmrs
      real        dtnlte(ngrid,nlayer)            ! output temp. tendencies

c
c Standard atmosphere variables
c
      real        nt                              ! number density [cm-3]
      real        co2(nlayer)                     !  "   of CO2
      real        o3p(nlayer)                     !  "   of atomic oxygen
      real        n2co(nlayer)                    !  "  of N2 + CO
      real        pyy(nlayer)                     ! auxiliary pressure grid

c
c Vectors and indexes for the tabulation of escape functions and VMR
c
c                 np                              ! # data points in tabulation
c                 pnb(np)                         ! Pressure in tabulation
c                 ef1(np)                         ! Esc.funct.#1, tabulated
c                 ef2(np)                         ! Esc.funct.#2, tabulated
c                 co2vmr(np)                      ! CO2 VMR tabulated
c                 o3pvmr(np)                      ! CO2 VMR tabulated
c                 n2covmr(np)                     ! N2+CO VMR tabulated
      real        escf1(nlayer)                   ! Esc.funct.#1, interpolated
      real        escf2(nlayer)                   ! Esc.funct.#2, interpolated


c
c Local Constants
c
      real       nu1, nu2                         ! freq. of energy levels
      real       imr1, imr2                       ! isotopic abundances
      real       hplanck, gamma, vlight           ! physical constants
      real       ee
      real       rfvt                             ! collisional rate
      real       rfvto3p                          !     "
      real       rfvv                             !     "

c
c Local variables for the main loop
c
      real       n1, n2, co2t                     ! ground populations
      real       l1, p1, p12                      ! prod & losses
      real       l2, p2, p21
      real       tt                               ! dummy variable
      real       c1, c2                           ! molecular constants
      real       ae1, ae2                         ! einstein spontaneous emission
      real       a1, a2, a12, a21
      real       pl1, pl2
      real       el1, el2
      real       hr1, hr2                         ! heating rate due to each band
      real       hr(nlayer)                       ! total heating rate

c
c Indexes
c
      integer    i
      integer    j,ii

c
c Rate coefficients
c
      real       k19xca, k19xcb
      real       k19cap1, k19cap2
      real       k19cbp1, k19cbp2
      real       d19c, d19cp1, d19cp2
      real       k20xc, k20cp1, k20cp2
      real       k21xc, k21cp2

      logical    firstcall
      data       firstcall/.true./
      save       firstcall,ef1,ef2,co2vmr,n2covmr,o3pvmr,pnb

c
c Data
c
      data       nu1, nu2, hplanck, gamma, vlight, ee/
     1     667.38, 662.3734, 6.6261e-27, 1.191e-5, 3.e10, 1.438769/

c*************************************************************************
c       PROGRAM  STARTS
c*************************************************************************

      imr1 = 0.987
      imr2 = 0.00408 + 0.0112
      rfvt = 0.1
      rfvto3p = 1.0
      rfvv = 0.1

      if(firstcall) then

         do i=1,np
            pnb(i)=1.0e-4*exp(pnb(i)) ! p into Pa
         end do

         firstcall = .false.

      endif

c
c MAIN LOOP, for each gridpoint and altitude:
c
      do j=1,ngrid  ! loop over grid points
c
c set up local pressure grid
c
         do ii=1,nlayer
            pyy(ii)=pplay(j,ii)
         enddo
!
! Interpolate escape functions and VMR to the desired grid
!
         call interp1(escf2,pyy,nlayer,ef2,pnb,np)
         call interp1(escf1,pyy,nlayer,ef1,pnb,np)
         if(nltemodel.eq.0) then
            call interp3(co2,o3p,n2co,pyy,nlayer,
     &           co2vmr,o3pvmr,n2covmr,pnb,np)
         endif
        
         do i=1,nlayer  ! loop over layers
C
C test if p lies outside range (p > 3.5 Pa)
C changed to 1 Pa since transition will always be higher than this
C
            if(pyy(i) .gt. 1.0 .or. pyy(i) .lt. 4.0e-6) then 
               hr(i)=0.0
               dtnlte(j,i)=0.0
            else
c
c           if(pt(j,i).lt.1.0)print*,pt(j,i)
               nt = pyy(i)/(1.381e-17*pt(j,i)) ! nt in cm-3
               !Dynamic composition
               if(nltemodel.eq.1) then
                  co2(i)=pq(j,i,igcm_co2)*mmean(j,i)/mmol(igcm_co2)
                  o3p(i)=pq(j,i,igcm_o)*mmean(j,i)/mmol(igcm_o)
                  n2co(i)=pq(j,i,igcm_co)*mmean(j,i)/mmol(igcm_co) +
     $                 pq(j,i,igcm_n2)*mmean(j,i)/mmol(igcm_n2)
               endif

               !Mixing ratio to density
               co2(i)=co2(i)*nt                ! CO2 density in cm-3
               o3p(i)=o3p(i)*nt                ! O3p density in cm-3
               n2co(i)=n2co(i)*nt              ! N2+CO in cm-3
c molecular populations
               n1 = co2(i) * imr1
               n2 = co2(i) * imr2
               co2t = n1 + n2

c intermediate collisional rates
               tt = pt(j,i)*pt(j,i)

               if (pt(j,i).le.175.0) then
                  k19xca = 3.3e-15
                  k19xcb = 7.6e-16
               else
                  k19xca = 4.2e-12 * exp( -2988.0/pt(j,i) + 303930.0/tt)
                  k19xcb = 2.1e-12 * exp( -2659.0/pt(j,i) + 223052.0/tt)
               endif
               k19xca = k19xca * rfvt
               k19xcb = k19xcb * rfvt
               k19cap1 = k19xca * 2.0 * exp( -ee*nu1/pt(j,i) )
               k19cap2 = k19xca * 2.0 * exp( -ee*nu2/pt(j,i) )
               k19cbp1 = k19xcb * 2.0 * exp( -ee*nu1/pt(j,i) )
               k19cbp2 = k19xcb * 2.0 * exp( -ee*nu2/pt(j,i) )
               d19c = k19xca*co2t + k19xcb*n2co(i)
               d19cp1 = k19cap1*co2t + k19cbp1*n2co(i)
               d19cp2 = k19cap2*co2t + k19cbp2*n2co(i)
                                !
               k20xc = 3.e-12 * rfvto3p
               k20cp1 = k20xc * 2.0 * exp( -ee/pt(j,i) * nu1 )
               k20cp2 = k20xc * 2.0 * exp( -ee/pt(j,i) * nu2 )
                                !
               k21xc = 2.49e-11 * 0.5 * rfvv
               k21cp2 = k21xc * exp( - ee/pt(j,i) * (nu2-nu1) )
                                !
               l1 = d19c + k20xc*o3p(i) + k21cp2*n2
               p1 = ( d19cp1 + k20cp1*o3p(i) ) * n1
               p12 = k21xc*n1
                                !
               l2 = d19c + k20xc*o3p(i) + k21xc*n1
               p2 = ( d19cp2 + k20cp2*o3p(i) ) * n2
               p21 = k21cp2*n2

c radiative rates
               ae1 = 1.3546 * 1.66 / 4.0 * escf1(i)
               ae2 = ( 1.3452 + 1.1878 ) * 1.66 / 4.0 * escf2(i)
               l1 = l1 + ae1
               l2 = l2 + ae2

c solving the system
               c1 = gamma*nu1**3. * 0.5
               c2 = gamma*nu2**3. * 0.5
               a1 = c1 * p1 / (n1*l1)
               a2 = c2 * p2 / (n2*l2)
               a12 = (nu1/nu2)**3. * n2/n1 * p12/l1
               a21 = (nu2/nu1)**3. * n1/n2 * p21/l2
               el2 = (a2 + a21 * a1 ) / ( 1.0 - a21 * a12 )
               el1 = a1 + a12 * el2
               pl1 = el1 * n1 / c1
               pl2 = el2 * n2 / c2

c  heating rate
               hr1 = - hplanck*vlight * nu1 * ae1 * pl1
               hr2 = - hplanck*vlight * nu2 * ae2 * pl2
               hr(i) = hr1 + hr2
               dtnlte(j,i)=0.1*hr(i)*pt(j,i)/(4.4*pyy(i)) ! dtnlte in K s-1
c              write(7,25)pxx(i),hr1,hr2,hr(i),qt
c  25         format(' ',1p5e12.4)

            endif

         enddo  ! end loop over layers
      enddo     ! end loop over grid points
c     close(7)
c
        return
        end

c***********************************************************************

      subroutine interp1(escout,p,nlayer,escin,pin,nl)
C
C subroutine to perform linear interpolation in pressure from 1D profile 
C escin(nl) sampled on pressure grid pin(nl) to profile
C escout(nlayer) on pressure grid p(nlayer).
C
      real escout(nlayer),p(nlayer)
      real escin(nl),pin(nl),wm,wp
      integer nl,nlayer,n1,n,nm,np
      do n1=1,nlayer
         if(p(n1) .gt. 3.5 .or. p(n1) .lt. 4.0e-6) then
            escout(n1) = 0.0
         else
            do n = 1,nl-1
               if (p(n1).le.pin(n).and.p(n1).ge.pin(n+1)) then
                  nm=n
                  np=n+1
                  wm=abs(pin(np)-p(n1))/(pin(nm)-pin(np))
                  wp=1.0 - wm
               endif
            enddo
            escout(n1) = escin(nm)*wm + escin(np)*wp
         endif
      enddo
      return
      end

c***********************************************************************

      subroutine interp3(esco1,esco2,esco3,p,nlayer,
     1     esci1,esci2,esci3,pin,nl)
C
C subroutine to perform 3 simultaneous linear interpolations in pressure from 
C 1D profiles esci1-3(nl) sampled on pressure grid pin(nl) to 1D profiles 
C esco1-3(nlayer) on pressure grid p(ngrid,nlayer).
C
      real esco1(nlayer),esco2(nlayer),esco3(nlayer),p(nlayer)
      real esci1(nl),    esci2(nl),    esci3(nl),    pin(nl),wm,wp
      integer nl,nlayer,n1,n,nm,np
      do n1=1,nlayer
         if (p(n1).gt. 3.5 .or. p(n1) .lt. 4.0e-6) then
            esco1(n1)=0.0
            esco2(n1)=0.0
            esco3(n1)=0.0
         else 
            do n = 1,nl-1
               if (p(n1).le.pin(n).and.p(n1).ge.pin(n+1)) then
                  nm=n
                  np=n+1
                  wm=abs(pin(np)-p(n1))/(pin(nm)-pin(np))
                  wp=1.0 - wm
               endif
            enddo
            esco1(n1) = esci1(nm)*wm + esci1(np)*wp
            esco2(n1) = esci2(nm)*wm + esci2(np)*wp
            esco3(n1) = esci3(nm)*wm + esci3(np)*wp
         endif
      enddo
      return
      end
