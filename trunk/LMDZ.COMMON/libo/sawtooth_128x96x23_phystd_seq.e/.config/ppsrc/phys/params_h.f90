










      module params_h

!  for use with kcm subroutines
!  RDW 22/09/11

      implicit none

      double precision rc,cp_n,m_n,m_v,rmn
      save cp_n,m_n,m_v,rmn

      logical ideal_v
      save ideal_v

      parameter (rc  = 8.314462d0) ! ideal gas constant [J mol^-1 K^-1]

      end module params_h
