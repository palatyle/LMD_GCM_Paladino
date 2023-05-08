c-----------------------------------------------------------------------
c INCLUDE fisice.h

      COMMON/fisice/dqsf,rice,rdust,zcondicea,pclc

      REAL dqsf(nqmx)       ! tendance of tracer on surface (e.g. kg.m-2)
      REAL rice(ngridmx,nlayermx) ! Estimated ice crystal radius (m)
      REAL rdust(ngridmx,nlayermx) ! Prescribed dust radius in each layer (m) 
c Variables used to define water ice scavenging by atmos. CO2 condensing
      REAL zcondicea(ngridmx,nlayermx)
c Some diagnostic variables 

      REAL pclc(ngridmx)  ! cloud cover ratio (used by radiative transfer, may depend of nlay in the future)
c-----------------------------------------------------------------------
