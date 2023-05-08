c-----------------------------------------------------------------------
c   INCLUDE 'nlteparams.h'
c
c   Parameters which govern the transition from LTE to NLTE radiation
c   tendencies.
c-----------------------------------------------------------------------

      real ptrans           ! central pressure for transition (Pa)
      parameter (ptrans = 0.1)
      real zw               ! half-width for transition (scale heights)
      parameter (zw = 0.5)
      real pminte           ! pressure up to which LTE is calculated (Pa)
      parameter (pminte = 0.4*ptrans)
c     almost one scale height above transition in worst case is very safe
      real zwi
      parameter (zwi = 2./zw)

c-----------------------------------------------------------------------
