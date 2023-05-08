*********************************************************
      double precision function sig(t)
      implicit none
*    this function computes the surface tension (N.m)   *
*   between water ice and air as a function of temp.    *
*********************************************************

      real t

      sig = (141. - 0.15 * dble(t)) * 1.e-3

      return
      end

