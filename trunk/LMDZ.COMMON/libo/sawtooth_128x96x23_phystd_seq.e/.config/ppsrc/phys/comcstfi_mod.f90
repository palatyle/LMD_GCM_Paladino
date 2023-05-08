










MODULE comcstfi_mod
IMPLICIT NONE
      
      REAL,SAVE :: pi ! something like 3.14159
      REAL,SAVE :: rad ! radius of the planet (m)
      REAL,SAVE :: g ! gravity (m/s2)
      REAL,SAVE :: r ! reduced gas constant (r=8.314511/(mugaz/1000.0))
      REAL,SAVE :: cpp ! Cp of the atmosphere
      REAL,SAVE :: rcp ! r/cpp
      REAL,SAVE :: mugaz ! molar mass of the atmosphere (g/mol)
      REAL,SAVE :: omeg ! planet rotation rate (rad/s)
      REAL,SAVE :: avocado ! something like 6.022e23
!$OMP THREADPRIVATE(pi,rad,g,r,cpp,rcp,mugaz,omeg,avocado)

END MODULE comcstfi_mod
