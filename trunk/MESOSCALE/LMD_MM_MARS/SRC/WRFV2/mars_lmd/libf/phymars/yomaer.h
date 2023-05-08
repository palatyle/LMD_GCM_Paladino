C*    Radiative characteristics of the aerosols
C
C   Shortwave
c   ~~~~~~~~~
c
c tauvis: dust optical depth at reference wavelength  ("longrefvis" set
c in dimradmars.h : typically longrefvis = 0.67E-6 m, as measured by Viking )

c For the "naerkind" kind of aerosol radiative properties : 
C QVISsQREF  :  Qext / Qext("longrefvis")   <--- For both solar bands
C omegavis   :  sinle scattering albedo     <--- For both solar bands
C gvis       :  assymetry factor            <--- For both solar bands
c
C   Longwave
c   ~~~~~~~~
c
c For the "naerkind" kind of aerosol radiative properties : 
c QIRsQREF :  Qext / Qext("longrefvis")     <--- For the nir bandes IR
c omegaIR  :  mean single scattering albedo <--- For the nir bandes IR
c gIR      :  mean assymetry factor         <--- For the nir bandes IR
c
      real tauvis
      real QVISsQREF(nsun,naerkind)
      real omegavis(nsun,naerkind),gvis(nsun,naerkind)
      real QIRsQREF(nir,naerkind)
      real omegaIR(nir,naerkind),gIR(nir,naerkind)

      COMMON/YOMAER/tauvis,QVISsQREF,omegavis,gvis,
     &              QIRsQREF,omegaIR,gIR
