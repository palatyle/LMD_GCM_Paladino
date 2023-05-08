! *****************************************************
! control parameters of the Martian thermal plume model
! *****************************************************
! Reference paper:
! A. Colaïtis, A. Spiga, F. Hourdin, C. Rio, F. Forget, and E. Millour. 
! A thermal plume model for the Martian convective boundary layer. 
! Journal of Geophysical Research (Planets), 118:1468-1487, July 2013.
! http://dx.doi.org/10.1002/jgre.20104
! http://arxiv.org/abs/1306.6215
! -----------------------------------------------------------------------
! Author : A. Colaitis 2011-01-05 (with updates 2011-2013)
! Institution : Laboratoire de Meteorologie Dynamique (LMD) Paris, France
! -----------------------------------------------------------------------
! Corresponding author : A. Spiga aymeric.spiga_AT_upmc.fr
! -----------------------------------------------------------------------

       module comtherm_h
       implicit none

       !------------------------------------------------------------
       !------------------------------------------------------------
       ! nsplit_thermals       ! Sub-timestep for the thermal plume model. 
                               ! Dependent on the timestep chosen for radiative transfer.
                               ! It is recommended to run with 96 timestep per day and 
                               ! call radiative transfer at each timestep, 
                               ! configuration in which thermals can run
                               ! very well with a sub-timestep of 10.
#ifdef MESOSCALE
       ! ---- MESOSCALE (timesteps < 200s)
       ! -- only a few subtimestep are needed because physical timestep is low
       INTEGER,PARAMETER :: nsplit_thermals = 4
#else
       ! ---- GCM
       ! -- we recommend here a value for 96 physical timesteps per day in the GCM
       !    (a larger value will not change results...
       !    but would cost more for negligible benefit in stability and accuracy)
       INTEGER,PARAMETER :: nsplit_thermals = 10
       ! -- if there is less than 96 physical timesteps per day we recommend
       !INTEGER,PARAMETER :: nsplit_thermals=35
#endif

       !------------------------------------------------------------
       !------------------------------------------------------------ 
       ! r_aspect_thermals     ! Mainly control the shape of the temperature profile
                               ! at the bottom of the mixed layer. Decreasing it goes toward
                               ! a convective-adjustment like profile.
                               ! (see paragraph 45 of paper and appendix S4)
       REAL,PARAMETER :: r_aspect_thermals = 1.

       !------------------------------------------------------------
       !------------------------------------------------------------ 
       ! qtransport_thermals  ! logical to activate tracer transport in thermals
       !
       LOGICAL,PARAMETER :: qtransport_thermals = .true. 

       !------------------------------------------------------------
       !------------------------------------------------------------ 
       ! dtke_thermals  ! logical to activate TKE transport in thermals
                        ! -- still experimental, for testing purposes only.
                        ! -- not used in current thermal plume models both on Earth and Mars.
       LOGICAL,PARAMETER :: dtke_thermals = .false.

       !------------------------------------------------------------
       !------------------------------------------------------------ 
       ! thermverbose  ! make thermal plume model more verbose
       LOGICAL,PARAMETER :: thermverbose = .false.


       ! ------------------------------------------------------------------------------------
       ! -------------- TUNING PARAMETERS FOR MARTIAN THERMALS MODEL ------------------------
       ! ------------------------------------------------------------------------------------
       ! Detrainment
       REAL,PARAMETER :: ad = 0.0004    ! D_2 in paper, see paragraph 44
       REAL,PARAMETER :: bd = -0.6697   ! D_1 in paper, see paragraph 44
       ! Entrainment
       REAL,PARAMETER :: ae = 0.03683   ! E_1 in paper, see paragraph 43
       REAL,PARAMETER :: be = 0.631631  ! E_2 in paper, see paragraph 43
       ! Downdraft
       REAL,PARAMETER :: fdfu=-0.8      ! downdraft to updraft mass flux ratio
                                        ! see paper paragraph 48
       REAL,PARAMETER :: omega=-0.03    ! omega. see paper paragraph 48
       ! Vertical velocity equation
       REAL,PARAMETER :: a1=1.          ! a in paper, see paragraph 41
       REAL,PARAMETER :: b1=0.0001      ! b in paper, see paragraph 41
       ! Inversion layer (same as a1 b1 actually)
       REAL :: a1inv=1.         ! a1 coeff in inversion layer
       REAL :: b1inv=0.0001     ! b1 coeff in inversion layer
       ! ------------------------------------------------------------------------------------
       ! ------------------------------------------------------------------------------------
       ! ------------------------------------------------------------------------------------

       end module comtherm_h

