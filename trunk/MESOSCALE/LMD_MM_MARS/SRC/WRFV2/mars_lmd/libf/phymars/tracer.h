c-----------------------------------------------------------------------
c INCLUDE 'tracer.h'

      character*10  noms(nqmx)  ! name of the tracer
      real mmol(nqmx)           ! mole mass of tracer (g/mol-1) 
      real radius(nqmx)   ! dust and ice particle radius (m)
      real qext(nqmx)     ! Single Scat. Extinction coeff at 0.67 um
      real alpha_lift(nqmx) ! saltation vertical flux/horiz flux ratio (m-1)
      real alpha_devil(nqmx) ! lifting coeeficient by dust devil

      real varian      ! Characteristic variance of log-normal distribution
      real r3n_q     ! used to compute r0 from number and mass mixing ratio
      real qextrhor(nqmx) ! Intermediate for computing opt. depth from q
      real rho_dust     ! Mars dust density (kg.m-3)
      real rho_ice     ! Water ice density (kg.m-3)
      real ref_r0        ! for computing reff=ref_r0*r0 (in log.n. distribution)

      real dryness(ngridmx)!"Dryness coefficient" for grnd water ice sublimation
      
      COMMON/tracer/radius,qext,alpha_lift,alpha_devil,noms,mmol,
     & varian,r3n_q,qextrhor,rho_dust,rho_ice,ref_r0,dryness
c-----------------------------------------------------------------------
