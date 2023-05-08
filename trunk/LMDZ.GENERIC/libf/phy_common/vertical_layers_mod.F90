! $Id: $

MODULE vertical_layers_mod

   REAL,SAVE             :: preff  ! reference surface pressure (Pa)
   REAL,SAVE             :: scaleheight ! atmospheric reference scale height (km)
   REAL,SAVE,ALLOCATABLE :: ap(:) ! hybrid (pressure contribution) coordinate 
                                  ! at layer interfaces (Pa)
   REAL,SAVE,ALLOCATABLE :: bp(:) ! hybrid (sigma contribution) coordinate 
                                  ! at layer interfaces
   REAL,SAVE,ALLOCATABLE :: aps(:) ! hybrid (pressure contribution) coordinate 
                                   ! at mid-layer (Pa)
   REAL,SAVE,ALLOCATABLE :: bps(:) ! hybrid (sigma contribution) coordinate 
                                   ! at mid-layer
   REAL,SAVE,ALLOCATABLE :: presnivs(:) ! reference pressure at mid-layer (Pa),
                                        ! based on preff, ap and bp
   REAL,SAVE,ALLOCATABLE :: pseudoalt(:) ! pseudo-altitude of model layers (km),
                                         ! based on preff and scaleheight
   
!$OMP THREADPRIVATE(preff,scaleheight,ap,bp,aps,bps,presnivs,pseudoalt)


CONTAINS

  SUBROUTINE init_vertical_layers(nlayer,preff_,scaleheight_,ap_,bp_,&
                                 aps_,bps_,presnivs_, pseudoalt_)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nlayer ! number of atmospheric layers
    REAL,INTENT(IN)    :: preff_ ! reference surface pressure (Pa)
    REAL,INTENT(IN)    :: scaleheight_ ! atmospheric scale height (km)
    REAL,INTENT(IN)    :: ap_(nlayer+1) ! hybrid coordinate at interfaces
    REAL,INTENT(IN)    :: bp_(nlayer+1) ! hybrid coordinate at interfaces
    REAL,INTENT(IN)    :: aps_(nlayer) ! hybrid coordinate at mid-layer
    REAL,INTENT(IN)    :: bps_(nlayer) ! hybrid coordinate at mid-layer
    REAL,INTENT(IN)    :: presnivs_(nlayer) ! Appproximative pressure of atm. layers (Pa)
    REAL,INTENT(IN)    :: pseudoalt_(nlayer) ! pseudo-altitude of atm. layers (km)
  
    ALLOCATE(ap(nlayer+1))
    ALLOCATE(bp(nlayer+1))
    ALLOCATE(aps(nlayer))
    ALLOCATE(bps(nlayer))
    ALLOCATE(presnivs(nlayer))
    ALLOCATE(pseudoalt(nlayer))
  
    preff = preff_
    scaleheight=scaleheight_
    ap(:) = ap_(:)
    bp(:) = bp_(:)
    aps(:) = aps_(:)
    bps(:) = bps_(:)
    presnivs(:) = presnivs_(:)
    pseudoalt(:) = pseudoalt_(:)

  END SUBROUTINE init_vertical_layers

END MODULE vertical_layers_mod
