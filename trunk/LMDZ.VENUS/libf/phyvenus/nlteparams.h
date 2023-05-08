c-----------------------------------------------------------------------
c   INCLUDE 'nlteparams.h'
c
c   Parameters which govern the transition from LTE to NLTE radiation
c   tendencies.
c-----------------------------------------------------------------------

      real ptrans           ! central pressure for transition (Pa)
c     "Nominal" values in Gilli+2017
c      parameter (ptrans = 0.5)
c     "New" values by D.Quirino  (better agreement with SOIR/SPICAV temperature)
      parameter (ptrans = 0.2)
      real zw               ! half-width for transition (scale heights)
      parameter (zw = 0.5)
      real pminte           ! pressure up to which LTE is calculated (Pa)
      parameter (pminte = 0.3*ptrans)
c     almost one scale height above transition in worst case is very safe
      real zwi
      parameter (zwi = 2./zw)

c-----------------------------------------------------------------------

c     Number of layers on which LTE calculations (in lw and sw) are performed 
c     (Computed in nlthermeq) :
      INTEGER nlaylte    
      common/nltepar/nlaylte
