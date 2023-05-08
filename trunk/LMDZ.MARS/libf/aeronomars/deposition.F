      subroutine deposition(ngrid, nlayer, nq,
     &                      ig, ig_vl1, pplay, pplev, zzlay, zzlev,
     $                      zu, zv, zt, zycol, ptimestep, co2ice)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     dry deposition of chemical species
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      use surfdat_h, only: z0 ! surface roughness
      use conc_mod, only: rnew ! specific gas constant
      implicit none
c
c
c     input
c 
      integer,intent(in) :: ngrid ! number of atmospheric columns
      integer,intent(in) :: nlayer ! number of atmospheric layers
      integer,intent(in) :: nq ! number of tracers
      integer ig                     ! grid point index
      integer ig_vl1                 ! viking 1 grid point
      real    pplay(ngrid,nlayer)    ! pressure at the middle of the layers (pa)
      real    pplev(ngrid,nlayer+1)  ! pressure at layer boundaries (pa)
      real    zzlay(ngrid,nlayer)    ! altitude at the middle of the layers (m)
      real    zzlev(ngrid,nlayer+1)  ! altitude at layer boundaries (m)
      real    zu(ngrid,nlayer)       ! u component of the wind (m.s-1)
      real    zv(ngrid,nlayer)       ! v component of the wind (m.s-1)
      real    zt(ngrid,nlayer)       ! temperature (k)
      real    zycol(nlayer,nq)       ! composition (volume mixing ratio)
      real    ptimestep              ! physical timestep (s)
      real    co2ice(ngrid)          ! co2 ice surface layer (kg.m-2)
c
c     local
c
      real    ubar                       ! wind module (m.s-1)
      real    ustar                      ! friction velocity (m.s-1)
      real    karman                     ! von karman constant
      real    ra                         ! aerodynamic resistance (s.m-1)
      real    rb                         ! quasi-laminar layer resistance (s.m-1)
      real    vd                         ! dry deposition velocity (cm.s-1)
      real    deltaz                     ! thickness of first layer (m)
      real    mu0                        ! dynamic viscosity of co2 at 293.15 k (pa.s)
      real    mu                         ! dynamic viscosity of co2 (pa.s)
      real    nu                         ! kinematic viscosity of co2 (cm2.s-1)
      real    rho                        ! density
      real    d                          ! molecular diffusivity of methane (cm2.s-1)
      real    schmidt                    ! schmidt number
      real    prandtl                    ! prandtl number
      real    p0                         ! reference pressure (pa)
      real    loss                       ! loss rate (s-1)                       
      integer iq
c
      real    nuch4
      real    tau
      real    gravity
      real    gam
      real    dp
      real    cd
c
      data    karman / 0.4 /
      data    prandtl / 0.72 /
      data    mu0 / 14.8e-6 /
      data    p0 /1.e5/
c
c     deposition is only active on surface uncovered by ice
c
c     if ((.not. watercaptag(ig)) .and. (co2ice(ig) .eq. 0.)) then
c
c        wind module (m.s-1)
c
         ubar = (zu(ig,1)**2. + zv(ig,1)**2.)**0.5
c
c        friction velocity (m.s-1)
c
         ustar = ubar*karman/log(zzlay(ig,1)/z0(ig))
c
c        aerodynamic resistance (s.m-1)
c
         ra = 1./(karman*ustar)*log(zzlay(ig,1)/z0(ig))
c
c        molecular diffusivity of methane in *air*  (cm2.s-1)
c        massman, atmospheric environment, 32, 1111-1127, 1998
c
         d = 0.1952*(p0/pplay(ig,1))*(zt(ig,1)/273.15)**1.81
c
c        dynamic viscosity: temperature dependance (pa.s)
c        sutherland's formula
c
         mu = mu0*(293.15 + 240.)/(zt(ig,1) + 240.)
     $           *(zt(ig,1)/293.15)**(3./2.)
c
c        density (kg.m-3)
c
         rho = pplay(ig,1)/(rnew(ig,1)*zt(ig,1))
c
c        kinematic viscosity (cm2.s-1) 
c
         nu = mu/rho*1.e4
c
c        schmidt number
c
         schmidt = nu/d
c
c        quasi-laminar layer resistance (s.m-1)
c
         rb = (2./(karman*ustar))*(schmidt/prandtl)**2./3.
c
c        dry deposition velocity (m.s-1)
c
         vd = 1./(ra + rb)
c
c        thickness of first layer (m)
c
         deltaz = zzlev(ig,2) - zzlev(ig,1)
c
c        loss rate (s-1)
c
         loss = vd/deltaz
c
c        test
c
c        loss = 1./(3600.*6.)
c        vd = loss*deltaz
c
         nuch4 = sqrt(8.*8.31*zt(ig,1)
     $                /(3.1416*16.e-3))
         tau = 6.*3600.
         gravity = 3.7
         dp = pplev(ig,1) - pplev(ig,2)
         cd = (karman/log(zzlay(ig,1)/z0(ig)))**2.
c
         gam = (4./nuch4)*dp
     $        /(tau*gravity*rho - 1./(cd*ubar))
c
c        methane index
c
         iq = 12
c
c        update methane in first layer
c
c        zycol(1,iq) = zycol(1,iq)*exp(-loss*ptimestep)
c
         if (ig .eq. ig_vl1) then
            print*,'**** deposition ****'
            print*,'z0     = ', z0(ig), 'm'
            print*,'deltaz = ', deltaz, ' m'
            print*,'deltap = ', dp, 'pa'
            print*,'zzlay  = ', zzlay(ig, 1), ' m'
            print*,'pplay  = ', pplay(ig, 1), ' pa'
            print*,'t      = ', zt(ig,1), ' k'
            print*,'u      = ', zu(ig,1), ' m.s-1'
            print*,'v      = ', zv(ig,1), ' m.s-1'
            print*,'rho    = ', rho, ' kg.m-3'
            print*,'ubar   = ', ubar, ' m.s-1'
            print*,'d      = ', d, ' cm2.s-1'
            print*,'mu     = ', mu, ' pa.s'
            print*,'nu     = ', nu, ' cm2.s-1'
            print*,'schmi  = ', schmidt
            print*,'Ra     = ', ra, ' s.m-1'
            print*,'Rb     = ', rb, ' s.m-1'
            print*,'vd     = ', vd*100., 'cm.s-1'
            print*,'tau dep= ', 1./loss, 's'
            print*,'R      = ', rnew(ig,1)
            print*,'nuch4  = ', nuch4, 'm.s-1'
            print*,'tau    = ', tau, ' s'
            print*,'gamma  = ', gam
            print*,'taugrho= ',tau*gravity*rho
            print*,'1surcdu= ', 1./(cd*ubar)

            print*,'********************'
         end if 
c
c     end if
c
      return
      end
