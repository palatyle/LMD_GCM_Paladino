      subroutine inichim_newstart(ngrid, nq, pq, qsurf, ps, &
                                  flagh2o, flagthermo)

      use tracer_mod
      USE vertical_layers_mod, ONLY: aps,bps
      USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat, nbp_lev
      USE datafile_mod, ONLY: datadir
      implicit none

!=======================================================================
!
!  Purpose:
!  --------
!
!  Initialization of the chemistry for use in the building of a new start file
!  used by program newstart.F
!  also used by program testphys1d.F
!
!  Authors: 
!  -------
!  Initial version 11/2002 by Sebastien Lebonnois
!  Revised 07/2003 by Monica Angelats-i-Coll to use input files
!  Modified 10/2008 Identify tracers by their names Ehouarn Millour
!  Modified 11/2011 Addition of methane Franck Lefevre
!  Rewritten 04/2012 Franck Lefevre
!
!  Arguments:
!  ----------
!
!  pq(nbp_lon+1,nbp_lat,nbp_lev,nq)  Advected fields, ie chemical species here
!  qsurf(ngrid,nq)     Amount of tracer on the surface (kg/m2)
!  ps(nbp_lon+1,nbp_lat)           Surface pressure (Pa)
!  flagh2o                 flag for initialisation of h2o (1: yes / 0: no)
!  flagthermo              flag for initialisation of thermosphere only (1: yes / 0: no)
!
!=======================================================================

      include "callkeys.h"

! inputs :

      integer,intent(in) :: ngrid         ! number of atmospheric columns in the physics
      integer,intent(in) :: nq                    ! number of tracers
      real,intent(in) :: ps(nbp_lon+1,nbp_lat)            ! surface pressure in the gcm (Pa)   
      integer,intent(in) :: flagh2o               ! flag for h2o initialisation
      integer,intent(in) :: flagthermo            ! flag for thermosphere initialisation only

! outputs :

      real,intent(out) :: pq(nbp_lon+1,nbp_lat,nbp_lev,nq)  ! advected fields, ie chemical species
      real,intent(out) :: qsurf(ngrid,nq)     ! surface values (kg/m2) of tracers

! local :

      integer :: iq, i, j, l, n
      integer :: count, ierr, dummy
      real    :: mmean(nbp_lon+1,nbp_lat,nbp_lev)             ! mean molecular mass (g)
      real    :: pgcm                             ! pressure at each layer in the gcm (Pa)

      integer, parameter         :: nalt = 252    ! number of altitudes in the initialization files
      integer                    :: nspe          ! number of species in the initialization files
      integer, allocatable       :: niq(:)        ! local index of species in initialization files
      real, dimension(nalt)      :: tinit, zzfile ! temperature in initialization files
      real, dimension(nalt)      :: pinit         ! pressure in initialization files
      real, dimension(nalt)      :: densinit      ! total number density in initialization files
      real, allocatable          :: vmrinit(:,:)  ! mixing ratios in initialization files
      real, allocatable          :: vmrint(:)     ! mixing ratio interpolated onto the gcm vertical grid
      real                       :: vmr

      character(len=20)          :: txt           ! to store some text
      logical                    :: flagnitro     ! checks if N species present

! 1. identify tracers by their names: (and set corresponding values of mmol)

! 1.1 initialize tracer indexes to zero:
      nqmx=nq ! initialize value of nqmx
      
      igcm_dustbin(1:nq)=0
      igcm_co2_ice=0
      igcm_ccnco2_mass=0
      igcm_ccnco2_number=0
      igcm_dust_mass=0
      igcm_dust_number=0
      igcm_ccn_mass=0
      igcm_ccn_number=0
      igcm_dust_submicron=0
      igcm_h2o_vap=0
      igcm_h2o_ice=0
      igcm_co2=0
      igcm_co=0
      igcm_o=0
      igcm_o1d=0
      igcm_o2=0
      igcm_o3=0
      igcm_h=0
      igcm_h2=0
      igcm_oh=0
      igcm_ho2=0
      igcm_h2o2=0
      igcm_ch4=0
      igcm_n2=0
      igcm_ar=0
      igcm_ar_n2=0
      igcm_n=0
      igcm_no=0
      igcm_no2=0
      igcm_n2d=0
      igcm_he=0
      igcm_co2plus=0
      igcm_oplus=0
      igcm_o2plus=0
      igcm_coplus=0
      igcm_cplus=0
      igcm_nplus=0
      igcm_noplus=0
      igcm_n2plus=0
      igcm_hplus=0
      igcm_hco2plus=0
      igcm_hcoplus=0
      igcm_elec=0

! 1.2 find dust tracers

      count = 0

      if (dustbin > 0) then
         do iq = 1,nqmx
            txt = " "
            write(txt,'(a4,i2.2)') 'dust', count + 1
            if (noms(iq) == txt) then
               count = count + 1
               igcm_dustbin(count) = iq
               mmol(iq) = 100.
            end if
         end do !do iq=1,nqmx
      end if ! of if (dustbin.gt.0)

      if (doubleq) then
         do iq = 1,nqmx
            if (noms(iq) == "dust_mass") then
               igcm_dust_mass = iq
               count = count + 1
            end if
            if (noms(iq) == "dust_number") then
               igcm_dust_number = iq
               count = count + 1
            end if
         end do
      end if ! of if (doubleq)

      if (scavenging) then
         do iq = 1,nqmx
            if (noms(iq) == "ccn_mass") then
               igcm_ccn_mass = iq
               count = count + 1
            end if
            if (noms(iq) == "ccn_number") then
               igcm_ccn_number = iq
               count = count + 1
            end if
         end do
      end if ! of if (scavenging)

      if (submicron) then
         do iq=1,nqmx
            if (noms(iq) == "dust_submicron") then
               igcm_dust_submicron = iq
               mmol(iq) = 100.
               count = count + 1
            end if
         end do
      end if ! of if (submicron)

! 1.3 find chemistry and water tracers

      do iq = 1,nqmx
         if (noms(iq) == "co2") then
            igcm_co2 = iq
            mmol(igcm_co2) = 44.
            count = count + 1
        end if
        if (noms(iq) == "co") then
           igcm_co = iq
           mmol(igcm_co) = 28.
           count = count + 1
        end if
        if (noms(iq) == "o") then
           igcm_o = iq
           mmol(igcm_o) = 16.
           count = count + 1
        end if
        if (noms(iq) == "o1d") then
           igcm_o1d = iq
           mmol(igcm_o1d) = 16.
           count = count + 1
        end if
        if (noms(iq) == "o2") then
           igcm_o2 = iq
           mmol(igcm_o2) = 32.
           count = count + 1
        end if
        if (noms(iq) == "o3") then
           igcm_o3 = iq
           mmol(igcm_o3) = 48.
           count = count + 1
        end if
        if (noms(iq) == "h") then
           igcm_h = iq
           mmol(igcm_h) = 1.
           count = count + 1
        end if
        if (noms(iq) == "h2") then
           igcm_h2 = iq
           mmol(igcm_h2) = 2.
           count = count + 1
        end if
        if (noms(iq) == "oh") then
           igcm_oh = iq
           mmol(igcm_oh) = 17.
           count = count + 1
        end if
        if (noms(iq) == "ho2") then
           igcm_ho2 = iq
           mmol(igcm_ho2) = 33.
           count = count + 1
        end if
        if (noms(iq) == "h2o2") then
           igcm_h2o2 = iq
           mmol(igcm_h2o2) = 34.
           count = count + 1
        end if
        if (noms(iq) == "ch4") then
           igcm_ch4 = iq
           mmol(igcm_ch4) = 16.
           count = count + 1
        end if
        if (noms(iq) == "n2") then
           igcm_n2 = iq
           mmol(igcm_n2) = 28.
           count = count + 1
        end if
        if (noms(iq) == "n") then
           igcm_n = iq
           mmol(igcm_n) = 14.
           count = count + 1
        end if
        if (noms(iq) == "n2d") then
           igcm_n2d = iq
           mmol(igcm_n2d) = 14.
           count = count + 1
        end if
        if (noms(iq) == "no") then
           igcm_no = iq
           mmol(igcm_no) = 30.
           count = count + 1
        end if
        if (noms(iq) == "no2") then
           igcm_no2 = iq
           mmol(igcm_no2) = 46.
           count = count + 1
        end if
        if (noms(iq) == "ar") then
           igcm_ar = iq
           mmol(igcm_ar) = 40.
           count = count + 1
        end if
        if (noms(iq) == "h2o_vap") then
           igcm_h2o_vap = iq
           mmol(igcm_h2o_vap) = 18.
           count = count + 1
        end if
        if (noms(iq) == "h2o_ice") then
           igcm_h2o_ice = iq
           mmol(igcm_h2o_ice) = 18.
           count = count + 1
        end if
        if (noms(iq).eq."he") then
          igcm_he=iq
          mmol(igcm_he)=4.
          count=count+1
        endif

! 1.4 find ions

        if (noms(iq) == "co2plus") then
           igcm_co2plus = iq
           mmol(igcm_co2plus) = 44.
           count = count + 1
        end if
        if (noms(iq) == "oplus") then
           igcm_oplus = iq
           mmol(igcm_oplus) = 16.
           count = count + 1
        end if
        if (noms(iq) == "o2plus") then
           igcm_o2plus = iq
           mmol(igcm_o2plus) = 32.
           count = count + 1
        end if
        if (noms(iq) == "coplus") then
           igcm_coplus = iq
           mmol(igcm_coplus) = 28.
           count = count + 1
        end if
        if (noms(iq) == "cplus") then
           igcm_cplus = iq
           mmol(igcm_cplus) = 12.
           count = count + 1
        end if
        if (noms(iq) == "nplus") then
           igcm_nplus = iq
           mmol(igcm_nplus) = 14.
           count = count + 1
        end if
        if (noms(iq) == "noplus") then
           igcm_noplus = iq
           mmol(igcm_noplus) = 30.
           count = count + 1
        end if
        if (noms(iq) == "n2plus") then
           igcm_n2plus = iq
           mmol(igcm_n2plus) = 28.
           count = count + 1
        end if
        if (noms(iq) == "hplus") then
           igcm_hplus = iq
           mmol(igcm_hplus) = 1.
           count = count + 1
        end if
        if (noms(iq) == "hco2plus") then
           igcm_hco2plus = iq
           mmol(igcm_hco2plus) = 45.
           count = count + 1
        end if
        if (noms(iq) == "hcoplus") then
           igcm_hcoplus = iq
           mmol(igcm_hcoplus) = 29.
           count = count + 1
        end if
        if (noms(iq) == "elec") then
           igcm_elec = iq
           mmol(igcm_elec) = 1./1822.89
           count = count + 1
        end if

! 1.5 find idealized non-condensible tracer

        if (noms(iq) == "Ar_N2") then
           igcm_ar_n2 = iq
           mmol(igcm_ar_n2) = 30.
           count = count + 1
        end if

      end do ! of do iq=1,nqmx
      
! 1.6 check that we identified all tracers:

      if (count /= nqmx) then
         write(*,*) "inichim_newstart: found only ",count," tracers"
         write(*,*) "                  expected ",nqmx
         do iq = 1,count
            write(*,*) '      ', iq, ' ', trim(noms(iq))
         end do
         stop
      else
         write(*,*) "inichim_newstart: found all expected tracers"
         do iq = 1,nqmx
            write(*,*) '      ', iq, ' ', trim(noms(iq))
         end do
      end if

! 1.7 check if nitrogen species are present:

      if(igcm_no == 0) then
         !check that no N species is in traceur.def
         if(igcm_n /= 0 .or. igcm_no2 /= 0 .or. igcm_n2d /= 0) then
            write(*,*)'inichim_newstart error:'
            write(*,*)'N, NO2 and/or N2D are in traceur.def, but not NO'
            write(*,*)'stop'
            stop
         endif
         flagnitro = .false.
         nspe = 14
      else
         !check that all N species are in traceur.def
         if(igcm_n == 0 .or. igcm_no2 == 0 .or. igcm_n2d == 0) then
            write(*,*)'inichim_newstart error:'
            write(*,*)'if NO is in traceur.def, N, NO2 and N2D must also be'
            write(*,*)'stop'
            stop
         endif
         flagnitro = .true.
         nspe = 18
      endif

! 1.8 allocate arrays

      allocate(niq(nspe))
      allocate(vmrinit(nalt,nspe))
      allocate(vmrint(nspe))

! 2. load in chemistry data for initialization:

! order of major species in initialization file:
!
!    1: co2 
!    2: ar
!    3: n2  
!    4: o2  
!    5: co  
!    6: o   
!    7: h2
!
! order of minor species in initialization file:
!
!    1: h  
!    2: oh 
!    3: ho2 
!    4: h2o 
!    5: h2o2 
!    6: o1d
!    7: o3
!
! order of nitrogen species in initialization file:
!
!    1: n
!    2: no
!    3: no2
!    4: n2d

! major species:

      niq(1) = igcm_co2
      niq(2) = igcm_ar
      niq(3) = igcm_n2
      niq(4) = igcm_o2
      niq(5) = igcm_co
      niq(6) = igcm_o
      niq(7) = igcm_h2

! minor species:

      niq(8)  = igcm_h
      niq(9)  = igcm_oh
      niq(10) = igcm_ho2
      niq(11) = igcm_h2o_vap
      niq(12) = igcm_h2o2
      niq(13) = igcm_o1d 
      niq(14) = igcm_o3

! nitrogen species:

      if (flagnitro) then
         niq(15) = igcm_n
         niq(16) = igcm_no
         niq(17) = igcm_no2
         niq(18) = igcm_n2d         
      end if

! 2.1 open initialization files

      open(210, iostat=ierr,file=trim(datadir)//'/atmosfera_LMD_may.dat')
      if (ierr /= 0) then
         write(*,*)'Error : cannot open file atmosfera_LMD_may.dat '
         write(*,*)'(in aeronomars/inichim_newstart.F)'
         write(*,*)'It should be in :', trim(datadir),'/'
         write(*,*)'1) You can change this path in callphys.def with'
         write(*,*)'   datadir=/path/to/datafiles/'
         write(*,*)'2) If necessary atmosfera_LMD_may.dat (and others)'
         write(*,*)'   can be obtained online on:'
         write(*,*)' http://www.lmd.jussieu.fr/~lmdz/planets/mars/datadir'
         stop
      end if
      open(220, iostat=ierr,file=trim(datadir)//'/atmosfera_LMD_min.dat')
      if (ierr /= 0) then
         write(*,*)'Error : cannot open file atmosfera_LMD_min.dat '
         write(*,*)'(in aeronomars/inichim_newstart.F)'
         write(*,*)'It should be in :', trim(datadir),'/'
         write(*,*)'1) You can change this path in callphys.def with'
         write(*,*)'   datadir=/path/to/datafiles/'
         write(*,*)'2) If necessary atmosfera_LMD_min.dat (and others)'
         write(*,*)'   can be obtained online on:'
         write(*,*)' http://www.lmd.jussieu.fr/~lmdz/planets/mars/datadir'
         stop
      end if
      if(flagnitro) then
         open(230, iostat=ierr,file=trim(datadir)//'/atmosfera_LMD_nitr.dat')
         if (ierr.ne.0) then
            write(*,*)'Error : cannot open file atmosfera_LMD_nitr.dat '
            write(*,*)'(in aeronomars/inichim_newstart.F)'
            write(*,*)'It should be in :', trim(datadir),'/'
            write(*,*)'1) You can change this path in callphys.def with'
            write(*,*)'   datadir=/path/to/datafiles/'
            write(*,*)'2) If necessary atmosfera_LMD_nitr.dat (and others)'
            write(*,*)'   can be obtained online on:'
            write(*,*)' http://www.lmd.jussieu.fr/~lmdz/planets/mars/datadir'
            STOP
         endif
      endif   ! Of if(flagnitro)

! 2.2 read initialization files

! major species

      read(210,*)
      do l = 1,nalt
         read(210,*) dummy, tinit(l), pinit(l), densinit(l), &
                     (vmrinit(l,n), n = 1,7)
         pinit(l) = pinit(l)*100.              ! conversion in Pa
         pinit(l) = log(pinit(l))              ! for the vertical interpolation
      end do
      close(210)

! minor species

      read(220,*)
      do l = 1,nalt
         read(220,*) dummy, (vmrinit(l,n), n = 8,14)
      end do 
      close(220)

! nitrogen species

      if (flagnitro) then
         read(230,*)
         do l = 1,nalt
            read(230,*) dummy, (vmrinit(l,n), n = 15,18)
         end do
         close(230)
      end if
      
! 3. initialization of tracers

      do i = 1,nbp_lon+1
         do j = 1,nbp_lat
            do l = 1,nbp_lev

               pgcm = aps(l) + bps(l)*ps(i,j)  ! gcm pressure
               pgcm = log(pgcm)                ! for the vertical interpolation
               mmean(i,j,l) = 0.

! 3.1 vertical interpolation

               do n = 1,nspe
                  call intrplf(pgcm,vmr,pinit,vmrinit(:,n),nalt)
                  vmrint(n) = vmr
                  iq = niq(n)
                  mmean(i,j,l) = mmean(i,j,l) + vmrint(n)*mmol(iq)
               end do

! 3.2 attribute mixing ratio: - all layers or only thermosphere
!                             - with our without h2o 

               if (flagthermo == 0 .or. (flagthermo == 1 .and. exp(pgcm) < 0.1)) then
                  do n = 1,nspe
                     iq = niq(n)
                     if (iq /= igcm_h2o_vap .or. flagh2o == 1) then
                        pq(i,j,l,iq) = vmrint(n)*mmol(iq)/mmean(i,j,l)
                     end if
                  end do
               end if

            end do
         end do
      end do

! set surface values of chemistry tracers to zero

      if (flagthermo == 0) then
         ! NB: no problem for "surface water vapour" tracer which is always 0
         do n = 1,nspe
            iq = niq(n)
            qsurf(1:ngrid,iq) = 0.
         end do
      end if

! 3.3 initialization of tracers not contained in the initialization files

! methane : 10 ppbv

      if (igcm_ch4 /= 0) then
         vmr = 10.e-9       
         do i = 1,nbp_lon+1
            do j = 1,nbp_lat
               do l = 1,nbp_lev
                  pq(i,j,l,igcm_ch4) = vmr*mmol(igcm_ch4)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngrid,igcm_ch4) = 0.
      end if

! ions: 0

      if (igcm_co2plus /= 0) then
         !check that all required ions are in traceur.def
         if (igcm_o2plus == 0 .or. igcm_oplus == 0 .or. igcm_coplus == 0          &
              .or. igcm_cplus == 0 .or. igcm_nplus == 0 .or. igcm_noplus == 0    & 
              .or. igcm_n2plus == 0 .or. igcm_hplus == 0 .or. igcm_hco2plus == 0 &
              .or. igcm_hcoplus == 0 .or. igcm_elec == 0) then
            write(*,*)'inichim_newstart error:'
            write(*,*)'if co2plus is in traceur.def, all other ions must also be'
            write(*,*)'o2plus, oplus, coplus, cplus, nplus, noplus, n2plus'
            write(*,*)'hplus, hco2plus, hcoplus, and elec'
            write(*,*)'stop'
            stop
         end if

         do i = 1,nbp_lon+1
            do j = 1,nbp_lat
               do l = 1,nbp_lev
                  ! all ions to 0     
                  pq(i,j,l,igcm_co2plus)  = 0.
                  pq(i,j,l,igcm_o2plus)   = 0.
                  pq(i,j,l,igcm_oplus)    = 0.
                  pq(i,j,l,igcm_coplus)   = 0.
                  pq(i,j,l,igcm_cplus)    = 0.
                  pq(i,j,l,igcm_nplus)    = 0.
                  pq(i,j,l,igcm_noplus)   = 0.
                  pq(i,j,l,igcm_n2plus)   = 0.
                  pq(i,j,l,igcm_hplus)    = 0.
                  pq(i,j,l,igcm_hco2plus) = 0.
                  pq(i,j,l,igcm_hcoplus)  = 0.
                  pq(i,j,l,igcm_elec)     = 0.
               end do
            end do
         end do

         ! surface value to 0

         qsurf(1:ngrid,igcm_co2plus)  = 0.
         qsurf(1:ngrid,igcm_o2plus)   = 0.
         qsurf(1:ngrid,igcm_oplus)    = 0.
         qsurf(1:ngrid,igcm_coplus)   = 0.
         qsurf(1:ngrid,igcm_cplus)    = 0.
         qsurf(1:ngrid,igcm_nplus)    = 0.
         qsurf(1:ngrid,igcm_noplus)   = 0.
         qsurf(1:ngrid,igcm_n2plus)   = 0.
         qsurf(1:ngrid,igcm_hplus)    = 0.
         qsurf(1:ngrid,igcm_hco2plus) = 0.
         qsurf(1:ngrid,igcm_hcoplus)  = 0.
         qsurf(1:ngrid,igcm_elec)     = 0.

      else

         if (igcm_o2plus /= 0 .or. igcm_oplus /= 0 .or. igcm_coplus /= 0          &
              .or. igcm_cplus /= 0 .or. igcm_nplus /= 0 .or. igcm_noplus /= 0    & 
              .or. igcm_n2plus /= 0 .or. igcm_hplus /= 0 .or. igcm_hco2plus /= 0 &
              .or. igcm_hcoplus /= 0 .or. igcm_elec /= 0) then
            write(*,*)'inichim_newstart error:'
            write(*,*)'some ions are in traceur.def, but not co2plus'
            write(*,*)'stop'
            stop
         end if
      end if    ! of if(igcm_co2 /= 0)
      
      ! deallocations

      deallocate(niq)
      deallocate(vmrinit)
      deallocate(vmrint)

      end
