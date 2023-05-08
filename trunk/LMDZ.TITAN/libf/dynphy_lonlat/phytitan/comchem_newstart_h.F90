MODULE comchem_newstart_h

! Stores data and routines for chemistry management specific to newstart

! Author : Jan Vatant d'Ollone (2018)

IMPLICIT NONE  
  
  ! Variable and allocatables for regriding chemistry pressure in newstart
  ! ----------------------------------------------------------------------
    
  INTEGER :: nlaykimold ! Number of upper atm. layers for chemistry in the start_archive file
  REAL, ALLOCATABLE, DIMENSION(:) :: preskimold ! Pressure grid of upper chemistry in the start_archive file
    
  ! Allocatable arrays for newstart fields
  ! --------------------------------------
  
  ! Nouvelle grille physique, ancienne grille verticale
  REAL,ALLOCATABLE :: ykim_up_oldv(:,:,:) ! (mol/mol)

  ! Nouvelle grille scalaire, ancienne grille verticale  
  REAL,ALLOCATABLE :: ykim_upS(:,:,:,:)  ! (mol/mol)

  ! Ancienne grille scalaire, ancienne grille verticale  
  REAL,ALLOCATABLE :: ykim_upoldS(:,:,:,:) ! (mol/mol)
  
END MODULE comchem_newstart_h
