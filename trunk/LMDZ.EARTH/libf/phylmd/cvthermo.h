!
! $Header$
!
c Thermodynamical constants for convectL:

      real cpd, cpv, cl, rrv, rrd, lv0, g, rowl, t0
      real clmcpv, clmcpd, cpdmcp, cpvmcpd, cpvmcl
      real eps, epsi, epsim1
      real ginv, hrd
      real grav

      COMMON /cvthermo/ cpd, cpv, cl, rrv, rrd, lv0, g, rowl, t0
     :                 ,clmcpv, clmcpd, cpdmcp, cpvmcpd, cpvmcl 
     :                 ,eps, epsi, epsim1, ginv, hrd, grav

c$OMP THREADPRIVATE(/cvthermo/)
