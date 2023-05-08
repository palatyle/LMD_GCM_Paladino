      SUBROUTINE orbite(pls,pdist_sol,pdecli)
      IMPLICIT NONE

c=======================================================================
c
c   Objet:
c   ------
c
c   Distance from sun and declimation as a function of the solar
c   longitude Ls
c
c   Interface:
c   ----------
c
c
c
c   Arguments:
c   ----------
c
c   Input:
c   ------
c   pls          Ls
c
c   Output:
c   -------
c   pdist_sol     Distance Sun-Planet in UA
c   pdecli        declinaison ( en radians )
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

#include "planete.h"
#include "comcstfi.h"

c arguments:
c ----------

      REAL pday,pdist_sol,pdecli,pls,i

c-----------------------------------------------------------------------

c Distance Sun-Planet

      pdist_sol=p_elips/(1.+e_elips*cos(pls+timeperi))

c Solar declination

c ********************* version before 01/01/2000 *******

      pdecli= asin (sin(pls)*sin(obliquit*pi/180.))

c ********************* version after 01/01/2000 *******
c     i=obliquit*pi/180.
c     pdecli=asin(sin(pls)*sin(i)/sqrt(sin(pls)**2+
c    & cos(pls)**2*cos(i)**2))
c ******************************************************



      RETURN
      END
