      SUBROUTINE orbite(pls,pdist_sol,pdecli)
      IMPLICIT NONE

c=======================================================================
c
c   Objet:
c   ------
c
c   Distance from sun and declination as a function of the solar
c   longitude Ls
c
c   Interface:
c   ----------
c
c    common 'comorbit.h' initialized by comorbit.f
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

#include "comorbit.h"

c arguments:
c ----------

      REAL pday,pdist_sol,pdecli,pls

c-----------------------------------------------------------------------

c Distance Sun-Planet

      pdist_sol=p_elips/(1.+e_elips*cos(pls+timeperi))

c Solar declination

      pdecli= asin (sin(pls)*sin(obliquit*pi/180.))

      RETURN
      END
