










      subroutine calc_rayleigh

!==================================================================
!     
!     Purpose
!     -------
!     Average the Rayleigh scattering in each band, weighting the 
!     average by the blackbody function at temperature tstellar.
!     Works for an arbitrary mix of gases (currently CO2, N2 and 
!     H2 are possible).
!     
!     Authors
!     ------- 
!     Robin Wordsworth (2010)
!     Jeremy Leconte (2012): Added option for variable gas. Improved water rayleigh (Bucholtz 1995).
!
!     Called by
!     ---------
!     setspv.F
!     
!     Calls
!     -----
!     none
!     
!==================================================================

      use radinc_h, only: L_NSPECTV
      use radcommon_h, only: WAVEV, BWNV, DWNV, tstellar, tauray, taurayvar, scalep
      use gases_h
      use comcstfi_mod, only: g, mugaz

      implicit none

      real*8 wl
      integer N,Nfine,ifine,igas
      parameter(Nfine=500.0)
      real*8 :: P0 = 9.423D+6   ! Rayleigh scattering reference pressure in pascals.    

      logical typeknown
      real*8 tauvar,tauvarvar,tausum,tausumvar,tauwei,tauweivar,bwidth,bstart
      double precision df

      real*8 tauconsti(ngasmx)
      real*8 tauvari(ngasmx)

      integer icantbewrong

      ! tau0/p0=tau/p (Hansen 1974)
      ! Then calculate tau0 = (8*pi^3*p_1atm)/(3*N0^2) * 4*delta^2/(g*mugaz*lambda^4)
      ! where delta=n-1 and N0 is an amagat

      typeknown=.false.
      do igas=1,ngasmx
         if(igas.eq.vgas)then
            print*,'variable gas is ',trim(gnom(igas)),' in Rayleigh scattering '
	 endif
	 if((igas/=vgas).and.(gfrac(igas).lt.5.e-2))then
            print*,'Ignoring ',trim(gnom(igas)),' in Rayleigh scattering '// &
            'as its mixing ratio is less than 0.05.' 
            ! ignore variable gas in Rayleigh calculation
            ! ignore gases of mixing ratio < 0.05 in Rayleigh calculation
            tauconsti(igas) = 0.0
         else
            if(igas.eq.igas_CO2) then
               tauconsti(igas) = (8.7/g)*1.527*scalep/P0
            elseif(igas.eq.igas_N2)then
               tauconsti(igas) = (9.81/g)*8.569E-3*scalep/(P0/93.0)
            elseif(igas.eq.igas_H2O)then
!               tauconsti(igas) = (10.0/g)*9.22E-3*scalep/101325.0 
               tauconsti(igas) = 1.98017E-10/(g*mugaz*100.) 
            elseif(igas.eq.igas_H2)then
               tauconsti(igas) = (10.0/g)*0.0241*scalep/101325.0
               !tauconsti(igas) = (10.0/g)*0.0487*scalep/(101325.0)
               ! uses n=1.000132 from Optics, Fourth Edition.
               ! was out by a factor 2!
            elseif(igas.eq.igas_He)then
               tauconsti(igas) = (10.0/g)*0.00086*scalep/101325.0
            else
               print*,'Error in calc_rayleigh: Gas species not recognised!'
               call abort
            endif

            if((gfrac(igas).eq.1.0).and.(vgas.eq.0))then
               print*,'Rayleigh scattering is for a pure ',trim(gnom(igas)),' atmosphere.'
               typeknown=.true.
            endif

         endif
      enddo

      if(.not.typeknown)then
         print*,'Rayleigh scattering is for a mixed gas atmosphere.'
         typeknown=.true.
      endif

 
      do N=1,L_NSPECTV
         
         tausum = 0.0
         tauwei = 0.0
         tausumvar = 0.0
         tauweivar = 0.0
         bstart = 10000.0/BWNV(N+1)
         bwidth = (10000.0/BWNV(N)) - (10000.0/BWNV(N+1))
         do ifine=1,Nfine
            wl=bstart+dble(ifine)*bwidth/Nfine

            tauvar=0.0
            tauvarvar=0.0
            do igas=1,ngasmx
               if((igas/=vgas).and.(gfrac(igas).lt.5.e-2))then
                  ! ignore variable gas in Rayleigh calculation
                  ! ignore gases of mixing ratio < 0.05 in Rayleigh calculation
                  tauvari(igas)   = 0.0
               else
                  if(igas.eq.igas_CO2)then
                     tauvari(igas) = (1.0+0.013/wl**2)/wl**4
                  elseif(igas.eq.igas_N2)then
                     tauvari(igas) = (1.0+0.0113/wl**2+0.00013/wl**4)/wl**4
                  elseif(igas.eq.igas_H2O)then
!                     tauvari(igas) = 1.0/wl**4 ! to be improved...
                     tauvari(igas) = (5.7918E6/(2.38E2-1/wl**2)+1.679E5/(57.36E0-1/wl**2))**2/wl**4 
                  elseif(igas.eq.igas_H2)then
                     tauvari(igas) = 1.0/wl**4 
                  elseif(igas.eq.igas_He)then
                     tauvari(igas) = 1.0/wl**4 
                  else
                     print*,'Error in calc_rayleigh: Gas species not recognised!'
                     call abort
                  endif
               endif

               if(igas.eq.vgas) then
                  tauvarvar=tauvarvar+tauconsti(igas)*tauvari(igas)
                  tauvar=tauvar+0.0*0.0*gfrac(igas)
               else
                  tauvar=tauvar+tauconsti(igas)*tauvari(igas)*gfrac(igas)
               endif

            enddo

            call blackl(dble(wl*1e-6),dble(tstellar),df)
            df=df*bwidth/Nfine
            tauwei=tauwei+df
            tausum=tausum+tauvar*df
            tauweivar=tauweivar+df
            tausumvar=tausumvar+tauvarvar*df
         
         enddo
         TAURAY(N)=tausum/tauwei
         TAURAYVAR(N)=tausumvar/tauweivar
         ! we multiply by scalep here (100) because plev, which is used in optcv,
         ! is in units of mBar, so we need to convert.

      end do

      IF (L_NSPECTV > 6) THEN
        icantbewrong = L_NSPECTV-6
        print*,'At 1 atm and lambda = ',WAVEV(icantbewrong),' um'
        print*,'tau_R = ',TAURAY(icantbewrong)*1013.25
        print*,'sig_R = ',TAURAY(icantbewrong)*g*mugaz*1.67e-27*100, &
               'cm^2 molecule^-1'
      ENDIF

    end subroutine calc_rayleigh
