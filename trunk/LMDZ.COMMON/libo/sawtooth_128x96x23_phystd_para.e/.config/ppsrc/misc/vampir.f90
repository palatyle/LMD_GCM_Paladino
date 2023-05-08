










module Vampir

  INTEGER,parameter :: VTcaldyn=1
  INTEGER,parameter :: VTintegre=2
  INTEGER,parameter :: VTadvection=3
  INTEGER,parameter :: VTdissipation=4
  INTEGER,parameter :: VThallo=5
  INTEGER,parameter :: VTphysiq=6
  INTEGER,parameter :: VTinca=7
  
  INTEGER,parameter :: nb_inst=7
  INTEGER :: MPE_begin(nb_inst)
  INTEGER :: MPE_end(nb_inst)
  
contains

  subroutine InitVampir
    implicit none


  end subroutine InitVampir

  subroutine VTb(number)
    implicit none
    INTEGER :: number

  end subroutine VTb

  subroutine VTe(number)
    implicit none
    INTEGER :: Number


  end subroutine VTe
  
end module Vampir
