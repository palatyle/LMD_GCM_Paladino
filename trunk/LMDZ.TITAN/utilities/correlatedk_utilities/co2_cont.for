c     This program computes RT CIA spectra of CO2 at low density limit 
C     at temperatures between 200 and 800K.
C     
C     ============================================================================
C     Copyright (C) 1997   Marcin Gruszka 
C     ============================================================================
C     
C     Copyright Notice:
C     You may use this program for your scientific applications,
C     but please do not	distribute it yourself.
C     Direct all your inquires to the  author: e-mail aborysow@stella.nbi.dk
C     Final request: If you publish  your work which benefited from this
C     program, please acknowledge using this program and quote 
C     the original paper describing the procedure:
C     M. Gruszka and A. Borysow,
C     Roto-translational collision-induced absorption of CO2
C     for the Atmosphere of Venus at Frequencies from 0 to 250 cm$^{-1}$
C     and at temperature from 200K to 800K,
C     Icarus, vol. 129, pp. 172-177, (1997). 

C     
C     The paper (Icarus) describing this computer model, is based on  original paper by
C     M. Gruszka and A. Borysow,
C     "Computer Simulation of the Far Infrared Collision-induced absorption spectra 
C     of gaseous CO2", Mol. Phys., vol. 93, pp. 1007-1016, 1998. 
C     ============================================================================
C     
C     NOTE: The program gives data at spacing corresponding to
c     Freq(max)-Freq(min)/(npoint-1),
c     Since reasons are unknown (author's secret), it is recommended to
c     request one more point that neccesary and get the right freq. step
c     ==================================================================

c     
c     The following input parametres are required (*):
c     ------------------------------------------------
c     (*) some input parameters are Fortran REAL numbers
c     and should be entered with a decimal point,
c     not like INTEGER!
c     
c     1. Temperature: (in K) (REAL)
c     -the model is restricted to the temperature range from
c     200 K to 800 K.
c     
c     2. Frequency limits: (in cm^-1) (REAL,REAL)
c     -the model has been designed to reproduce the MD and the
c     experimental data up to 250 cm^-1. There is no restriction
c     on requesting a higher value, but the accuracy may decrease 
c     significantly.
c     
c     3. Number of points in the output file: (less than Npmx) (INTEGER)
c     -the frequency grid is controled by the number of points
c     in the output file.
c     (i.e. if the upper limit = 250.0 and the number of points = 50.0,
c     the spectrum will be computed every 5 cm^-1 )
c     
c     4. The output in automaticaly directed to output.dat file.
c     
c     STRUCTURE OF THE PROGRAM:
c     -------------------------------------------------------------- 
c     The MAIN program consists of three modules:
c     
c     1. getdat - provides the user interface, 
c     in this subroutine the Temperature, Frequency Limits and
c     the number of points are read from the keyboard.
c     
c     2. getspc - the main subroutine, where the spectrum is computed.
c     The three input parameters (temp, range and n; all REAL numbers)
c     are transmited by value and the absorption coefficient abcoef(..) 
c     is returned. The default lenght of abcoef(..) is Npmx and if 
c     n < Npmx, the empty places are set to 0. If a denser grid is
c     required please change the n appropriaetly in the data declaration
c     of the MAIN program and inside the getspc.
c     
c     3. putdat - generates the output file.
c     The frequency (in cm^-1) and 
c     the absorption coefficient (in cm^-1 amagat^-2) are written
c     to the 'output.dat' file in the followinf format:
c     format(f10.3,e20.7)
c     
c===========================================================================
c     
c     
c     ---------------------------------------------------------------------------------------------------------------
c     2008/10/15
c     
c     Modifications of the code by V. Eymet as follows:
c     + only the "getspc" routine remains
c     + inputs: temp (temperature in K), wn (wavenumber in cm^-1); output: abcoef (absorption coefficient, in cm^-1 amagat^-2)
c     ---------------------------------------------------------------------------------------------------------------
c
c
c
c     ---------------------------------------------------------------------------------------------------------------
c     2009/2/28
c     
c     Modifications of the code by R. Wordsworth as follows:
c     + "baranov" routine added to calculate CIA dimer spectrum between 1100 and 1600 cm^-1
c     + this uses recorded data rather than an analytical function, so bilinear interpolation is needed
c     + inputs and outputs are identical to those of getspc
c     ---------------------------------------------------------------------------------------------------------------



      subroutine getspc(temp,wn,abcoef) ! computes the spectrum
      implicit none
      real*8 temp               !temperature
      real*8 wn        !wavenumber
      real*8 abcoef       !absorption coefficient
c     -------------------------------------------------------
c     The A, B and C parameters as given in the paper:
      integer ndeg
      parameter (ndeg=3)
      real*8 ah(ndeg),bh(ndeg),at(ndeg),bt(ndeg),gam(ndeg)
c     tau^{L}_{1}:
      data ah /0.1586914382D+13,-0.9344296879D+01,0.6943966881D+00/
c     tau^{L}_{2}:
      data bh /0.1285676961D-12,0.9420973263D+01,-0.7855988401D+00/
c     tau^{H}_{1}:
      data at /0.3312598766D-09,0.7285659464D+01,-0.6732642658D+00/
c     tau^{H}_{2}:
      data bt /0.1960966173D+09,-0.6834613750D+01,0.5516825232D+00/
c     gamma_{1}:
      data gam /0.1059151675D+17,-0.1048630307D+02, 0.7321430968D+00/
c     **************************************************************
c     The S - shape function parameters:
      real*8 w1,w2
      data w1,w2 /50.0d0, 100.0d0/
c     --------------------------------------------------------------
c     local variables:
c     ----------------------------------
      real*8 incrmt
      integer p1,p2,n
      real*8 a1,b1,a2,b2,gamma
      real*8 bc1,bc2,bcbc
      real*8 scon,mtot,spunit,gm0con
      real*8 icm, frq, lntemp
      real*8 xk1,x1,x2
      external xk1
      incrmt=250.0D+0
      n=1
      p1=idint(dnint(w1/incrmt))
      p2=idint(dnint(w2/incrmt))
      lntemp=dlog(temp)
      a1=1.0d0/(ah(1)*dexp(ah(2)*lntemp+ah(3)*lntemp*lntemp))
      b1=bh(1)*dexp(bh(2)*lntemp+bh(3)*lntemp*lntemp)
      a2=1.0d0/(at(1)*dexp(at(2)*lntemp+at(3)*lntemp*lntemp))
      b2=bt(1)*dexp(bt(2)*lntemp+bt(3)*lntemp*lntemp)
      gamma=gam(1)*dexp(gam(2)*lntemp+gam(3)*lntemp*lntemp)
      icm=0.1885d0
      frq=wn
      x1=b1*dsqrt(a1*a1+icm*icm*frq*frq)
      bc1=dexp(a1*b1)*a1*XK1(x1)/(a1*a1+icm*icm*frq*frq)
      x2=b2*dsqrt(a2*a2+icm*icm*frq*frq)
      bc2=dexp(a2*b2)*a2*XK1(x2)/(a2*a2+icm*icm*frq*frq)
      if (p1.ge.1) then
         bcbc=bc1
      else if (p2.lt.1) then
         bcbc=bc2
      else
         bcbc=dexp( (1.0d0-dble(1-p1)/dble(p2-p1))*dlog(bc1)+
     &        dble(1-p1)/dble(p2-p1)*dlog(bc2) )
      endif
      spunit=1.296917d55
      gm0con=1.259009d-6
      scon=spunit/temp
      mtot=gm0con*temp*gamma*1.0d-56
      frq=wn
      abcoef=scon*mtot*bcbc*frq*frq
      return
      end
c     
c     
      FUNCTION XK1(X)
C     MODIFIED BESSEL FUNCTION K1(X) TIMES X
C     PRECISION IS BETTER THAN 2.2e-7 EVERYWHERE.
C     ABRAMOWITZ AND S,TEGUN, P.379; TABLES P.417.
      implicit double precision (a-h,o-z)
!      IF(X-2.) 10,10,20
      IF(X-2..le.0.0)then
! 10   T=(X/3.75)**2
      T=(X/3.75)**2
      FI1=X*((((((.00032411*T+.00301532)*T+.02658733)*T+.15084934)
     1     *T+.51498869)*T+.87890594)*T+.5)
      T=(X/2.)**2
      P=(((((-.00004686*T-.00110404)*T-.01919402)*T-.18156897)*T-
     1     .67278579)*T+.15443144)*T+1.
      XK1=X*dLOG(X/2)*FI1+P
!      RETURN
      else
! 20   T=2./X
      T=2./X
      P=(((((-.00068245*T+.00325614)*T-.00780353)*T+.01504268)*T-
     1     .03655620)*T+.23498619)*T+1.25331414
      X=dMIN1(X,330.d0)
      XK1=dSQRT(X)*dEXP(-X)*P
      endif
      RETURN
      END




c     ----------------------------------------------------------------------------------------------
c     Added by RDW for inclusion of Baranov (2004) data
      subroutine baranov(temp,wn,abcoef)
      ! computes the dimer spectrum

      implicit none
      real*8 temp               !temperature
      real*8 wn                 !wavenumber
      real*8 abcoef             !absorption coefficient

      integer nS,nT
!      parameter(nS=1713)
      parameter(nS=3000)
      parameter(nT=9)

      real*8 wn_arr(nS)
      real*8 temp_arr(nT)
      real*8 dim_arr(nS,nT)

      integer pop
      logical firstcall

      save wn_arr, temp_arr, dim_arr

      character*100 dt_file
      integer strlen,ios
      character*100 label
      label='subroutine baranov'

      dt_file='/san/home/rdword/kspectrum/line_data/CO2_dimer_data_NEW'
!      dt_file='./data/CO2_dimer_data'
      open(33,file=dt_file(1:strlen(dt_file)),
     &     form='unformatted',
     &     status='old',iostat=ios)
      if (ios.ne.0) then        ! file not found
         write(*,*) 'Error from ',label(1:strlen(label)),' :'
         write(*,*) 'Data file could not be found:'
         write(*,*) dt_file(1:strlen(dt_file))
         stop
      else
         read(33) wn_arr
         read(33) temp_arr
         read(33) dim_arr
      endif
      close(33)


!	print*,wn_arr(1)!
!	print*,wn_arr(3000)
!	print*,temp_arr
!	print*,dim_arr(1500,1)
!	stop


      call bilinear(wn_arr,temp_arr,nS,nT,dim_arr,wn,temp,abcoef)

 111  continue
      return
      end
c     ----------------------------------------------------------------------------------------------

c     ----------------------------------------------------------------------------------------------
c     Added by RDW for inclusion of Baranov (2004) data
      subroutine bilinear(x_arr,y_arr,nX,nY,f2d_arr,x,y,f)

      implicit none

      integer nX,nY,i,j,a,b

      real*8 x,y,x1,x2,y1,y2
      real*8 f,f11,f12,f21,f22,fA,fB
      real*8 x_arr(nX)
      real*8 y_arr(nY)
      real*8 f2d_arr(nX,nY)
      
      integer strlen
      character*100 label
      label='subroutine bilinear'

c     1st check we're within the wavenumber range
      if ((x.lt.x_arr(2)).or.(x.gt.x_arr(nX-2))) then
         f=0.0D+0
c     print*,'x=',x,x_Arr(2),x_arr(nX-2)
c     stop('Outside CIA dimer wavenumber range!\n')
         goto 112
      else
         
c     in the x (wavenumber) direction 1st
         i=1
 10      if (i.lt.(nX+1)) then ! what passes for a 'while' loop in this sorry language...
            if (x_arr(i).gt.x) then
               x1=x_arr(i-1)
               x2=x_arr(i)
               a=i-1
               i=9999
            endif
            i=i+1
            goto 10
         endif
      endif
      
      if ((y.lt.y_arr(1)).or.(y.gt.y_arr(nY))) then
         write(*,*) 'Error from ',label(1:strlen(label)),' :'
         write(*,*) 'Outside CIA dimer temperature range!'
         if(y.lt.y_arr(1))then
            y=y_arr(1)+0.01
         endif
         if(y.gt.y_arr(nY))then
            y=y_arr(nY)-0.01
         endif
!         stop
      else
c     in the y (temperature) direction 2nd
         j=1
 20      if (j.lt.(nY+1)) then
            if (y_arr(j).gt.y) then
               y1=y_arr(j-1)
               y2=y_arr(j)
               b=j-1
               j=9999
            endif
            j=j+1
            goto 20
         endif
      endif
      
      f11=f2d_arr(a,b)
      f21=f2d_arr(a+1,b)
      f12=f2d_arr(a,b+1)
      f22=f2d_arr(a+1,b+1)
      
c     1st in x-direction
      fA=f11*(x2-x)/(x2-x1)+f21*(x-x1)/(x2-x1)
      fB=f12*(x2-x)/(x2-x1)+f22*(x-x1)/(x2-x1)
      
c     then in y-direction
      f=fA*(y2-y)/(y2-y1)+fB*(y-y1)/(y2-y1)
      
 112  continue
      return
      end
c     ----------------------------------------------------------------------------------------------
      
