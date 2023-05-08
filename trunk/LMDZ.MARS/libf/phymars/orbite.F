      SUBROUTINE orbite(pls,pdist_sol,pdecli)
      USE planete_h, ONLY: e_elips, p_elips, obliquit, timeperi
      USE comcstfi_h, ONLY: pi
      IMPLICIT NONE

c=======================================================================
c
c   Compute distance to Sun and declination as a function of the solar
c   longitude Ls
c
c   Arguments:
c   ----------
c
c   Input:
c   ------
c   pls          Ls (radians)
c
c   Output:
c   -------
c   pdist_sol     Distance Sun-Planet in UA
c   pdecli        declinaison ( in radians )
c
c=======================================================================

c arguments:
c ----------

      REAL,INTENT(IN) :: pls
      REAL,INTENT(OUT) :: pdist_sol,pdecli

c-----------------------------------------------------------------------

c Distance Sun-Planet

      pdist_sol=p_elips/(1.+e_elips*cos(pls+timeperi))

c Solar declination

      pdecli= asin (sin(pls)*sin(obliquit*pi/180.))

      END
