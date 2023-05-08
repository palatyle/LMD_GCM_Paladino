      subroutine gr_int_dyn(champin,champdyn,iim,jp1)
      implicit none
c=======================================================================
c   passage d'un champ interpole a un champ sur grille scalaire
c=======================================================================
c-----------------------------------------------------------------------
c   declarations:
c   -------------

      INTEGER iim
      integer ip1, jp1
      REAL champin(iim, jp1)
      REAL champdyn(iim+1, jp1)

      INTEGER i, j
      real polenord, polesud

c-----------------------------------------------------------------------
c   calcul:
c   -------

      ip1 = iim + 1
      polenord = 0.
      polesud = 0.
      do i = 1, iim
        polenord = polenord + champin (i, 1)
        polesud = polesud + champin (i, jp1)
      enddo
      polenord = polenord / iim
      polesud = polesud / iim
      do j = 1, jp1
        do i = 1, iim
          if (j .eq. 1) then
            champdyn(i, j) = polenord
          else if (j .eq. jp1) then
            champdyn(i, j) = polesud
          else
            champdyn(i, j) = champin (i, j)
          endif
        enddo
        champdyn(ip1, j) = champdyn(1, j)
      enddo

      RETURN
      END

