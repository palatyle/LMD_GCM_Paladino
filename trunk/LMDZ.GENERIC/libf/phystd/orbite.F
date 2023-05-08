      subroutine orbite(pls,pdist_star,pdecli,pright_ascenc)

      use planete_mod, only: p_elips, e_elips, timeperi, obliquit
      use comcstfi_mod, only: pi
      implicit none
!==================================================================
!     
!     Purpose
!     -------
!     Distance from star and declination as a function of the stellar
!     longitude Ls
!     
!     Inputs
!     ------
!     pls          Ls
!
!     Outputs
!     -------
!     pdist_star    Distance Star-Planet in UA
!     pdecli        declinaison ( in radians )
!     pright_ascenc right ascension ( in radians )
!
!=======================================================================

c   Declarations:
c   -------------

c arguments:
c ----------

      REAL pday,pdist_star,pdecli,pright_ascenc,pls,i

c-----------------------------------------------------------------------

c Star-Planet Distance

      pdist_star = p_elips/(1.+e_elips*cos(pls+timeperi))

c Stellar declination

c ********************* version before 01/01/2000 *******

      pdecli = asin (sin(pls)*sin(obliquit*pi/180.))

c********************* version after 01/01/2000 *******
c     i=obliquit*pi/180.
c     pdecli=asin(sin(pls)*sin(i)/sqrt(sin(pls)**2+
c    & cos(pls)**2*cos(i)**2))
c ******************************************************

c right ascencion
      If((pls.lt.pi/2.d0)) then
         pright_ascenc= atan(tan(pls)*cos(obliquit*pi/180.))
      else if((pls.gt.pi/2.d0).and.(pls.lt.3.d0*pi/2.d0)) then
         pright_ascenc= pi+atan(tan(pls)*cos(obliquit*pi/180.))
      else if((pls.gt.3.d0*pi/2.d0)) then
         pright_ascenc= 2.d0*pi+atan(tan(pls)*cos(obliquit*pi/180.))
      else if (Abs(pls-pi/2.d0).le.1.d-10) then
         pright_ascenc= pi/2.d0 
      else 
         pright_ascenc=-pi/2.d0 
      end if
      	 
      RETURN
      END
