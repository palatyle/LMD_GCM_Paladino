! $Id$
!
MODULE aero_mod
  ! Declaration des indices pour les aerosols 

  ! Total number of aerosols
  INTEGER, PARAMETER :: naero_tot = 10 

  ! Identification number used in aeropt_2bands and aeropt_5wv
  ! corresponding to naero_tot
  INTEGER, PARAMETER :: id_ASBCM    = 1
  INTEGER, PARAMETER :: id_ASPOMM   = 2
  INTEGER, PARAMETER :: id_ASSO4M   = 3
  INTEGER, PARAMETER :: id_CSSO4M   = 4
  INTEGER, PARAMETER :: id_SSSSM    = 5
  INTEGER, PARAMETER :: id_CSSSM    = 6
  INTEGER, PARAMETER :: id_ASSSM    = 7
  INTEGER, PARAMETER :: id_CIDUSTM  = 8
  INTEGER, PARAMETER :: id_AIBCM    = 9
  INTEGER, PARAMETER :: id_AIPOMM   = 10

  ! Total number of aerosols actually used in LMDZ 
  ! 1 =  ASBCM
  ! 2 =  ASPOMM
  ! 3 =  ASSO4M ( = SO4) 
  ! 4 =  CSSO4M 
  ! 5 =  SSSSM 
  ! 6 =  CSSSM
  ! 7 =  ASSSM
  ! 8 =  CIDUSTM
  ! 9 =  AIBCM
  !10 =  AIPOMM
  INTEGER, PARAMETER :: naero_spc = 10

  ! Corresponding names for the aerosols
  CHARACTER(len=7),DIMENSION(naero_spc) :: name_aero=(/&
       "ASBCM  ", &
       "ASPOMM ", &
       "SO4    ", &
       "CSSO4M ", &
       "SSSSM  ", &
       "CSSSM  ", &
       "ASSSM  ", &
       "CIDUSTM", &
       "AIBCM  ", &
       "AIPOMM " /)


  ! Number of aerosol groups
  ! 1 = ZERO    
  ! 2 = AER total    
  ! 3 = NAT    
  ! 4 = BC    
  ! 5 = SO4    
  ! 6 = POM    
  ! 7 = DUST    
  ! 8 = SS    
  ! 9 = NO3    
  INTEGER, PARAMETER :: naero_grp = 9 

  ! Number of  wavelengths
  INTEGER, PARAMETER :: nwave = 5

  ! Number of modes spectral bands
  INTEGER, parameter :: nbands = 2


END MODULE aero_mod
