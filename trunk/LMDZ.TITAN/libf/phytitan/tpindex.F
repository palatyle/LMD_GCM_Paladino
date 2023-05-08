      subroutine tpindex(pw,tw,pref,tref,LCOEF,MT,MP)

!==================================================================
!     
!     Purpose
!     -------
!     Interpolate K-coefficients to the given P,T.
!
!     Notes
!     -----
!     The interpolation is the usual one in two dimensions given
!     in "Numerical Recipes", where the "X" are P, the "Y" are
!     T, and the F(X,Y) are the CO2 K-coefficients.
!
!     The interpolating box is:
!
!           (PL,TU)                        (PR,TU)
!
!                          (TW,PW)
!
!           
!           (PL,TL)                        (PR,TL)
!
!      PL  - Pressure left
!      PR  - Pressure right
!      TL  - Temperature lower
!      TU  - Temperature upper
!      PW  - Pressure wanted
!      TW  - Temperature wanted
!
!     Inputs
!     ------
!     PW                 - The pressure to interpolate to
!     TW                 - The temperature to interpolate to
!     Pref(NP)           - The pressure grid array
!     Tref(NT)           - The temperature grid array
!   
!     Outputs
!     -------
!     TI                 - Interpolation term (pressure)
!     UI                 - Interpolation term (temperature)
!     MT                 - Temperature index (bottom left temperature)
!                          of bounding box
!     MP                 - Pressure index (bottom left pressure)
!                          of bounding box
!
!     Authors
!     -------
!     Adapted from the NASA Ames code by R. Wordsworth (2009)
!     
!==================================================================

      use radinc_h

      implicit none

      real*8 Tref(L_NTREF)
      real*8 pref(L_PINT)

      integer MT, MP, N, M, NP
      real*8  PW, TW
      real*8  PWL, LCOEF(4), T, U

C======================================================================C
 
!     Get the upper and lower temperature grid indicies that bound the
!     requested temperature. If the requested temperature is outside
!     the T-grid, set up to extrapolate from the appropriate end.
!     TW : temperature to be interpolated
!     TREF : grid array
!     MT : index of TREF for bounding new temperature
!     U : new index (real) for temperature interpolated

      IF(TW.LE.TREF(1)) THEN
        MT = 1
        IF (TW.LT.TREF(1)) THEN
         write(*,*) 'tpindex: Caution! Temperature of upper levels lower 
     $ than ref temperature for k-coef: k-coeff fixed for upper levels'
         write(*,*) "         TW=",TW
         write(*,*) "         TREF(1)=",TREF(1)
        ENDIF
      ELSE
        do n=1,L_NTREF-1
          if(tw.gt.Tref(n) .and. TW.LE.TREF(N+1)) then
            MT = n
            goto 10
          end if
        end do

        MT = L_NTREF-1
      
   10   continue
      END IF

      !TB15 : case low temp : MT=1: fixed TW right above tref(1)
      IF (MT.eq.1) THEN
         TW=tref(1)*1.00
!         write(*,*) 'tpindex: Caution! Temperature of upper levels lower 
!     $than ref temperature for k-coef: k-coeff fixed for upper levels'
!         write(*,*) "         TW=",TW
!         write(*,*) "         TREF(1)=",TREF(1)
      ENDIF

      U = (TW-TREF(MT))/(TREF(MT+1)-TREF(MT))

!     Get the upper and lower pressure grid indicies that bound the
!     requested pressure. If the requested pressure is outside
!     the P-grid, set up to extrapolate from the appropriate end.

      pwl = log10(pw)

      do n=2,L_PINT-1
        if(pwl.le.Pref(n)) then
          MP = n-1
          goto 20
        end if
      end do

      MP = L_PINT-1

   20 continue
      
      !TB15 : case low pressure : n=2 : fixed pwl, right above pref(1)
      IF (MP.eq.1) THEN
        IF (PWL.LT.PREF(1)) THEN
         write(*,*) 'tpindex: Caution! Pressure of upper levels lower 
     $than ref pressure for k-coef: k-coeff fixed for upper levels'
         write(*,*) "         PWL=",PWL
         write(*,*) "         PREF(1)=",PREF(1)
        ENDIF
        PWL=Pref(1)*1.00
      ENDIF

!     interpolated pressure
      T = (PWL-PREF(MP))/(PREF(MP+1)-PREF(MP))

!  Fill in the interpolation coefficients
      LCOEF(1) = (1.0-T)*(1.0-U)
      LCOEF(2) = T*(1.0-U)
      LCOEF(3) = T*U
      LCOEF(4) = (1.0-T)*U

   30 CONTINUE

      return
      end
