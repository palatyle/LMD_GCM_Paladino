c***********************************************

c	param.par

c	Parameters for paramhr.f
c***********************************************

	integer    ninter
	parameter  (ninter=33)

	integer    nabs
	parameter  (nabs=6)

	integer    nz2
	parameter  (nz2=203)

	integer    ninter2
        parameter  (ninter2=16)

	real       kboltzman                  !cte Boltzman
	parameter  (kboltzman = 1.381e-16)
	
	real       n_avog                    !# de Avogadro
	parameter  (n_avog = 6.023e23)

	real       gg                        !cte gravitacion
	parameter  (gg = 6.67259e-8)

	real       masa                      !masa de Marte(g)
	parameter  (masa = 6.4163e26)

	real       radio
	parameter  (radio = 3390.)           !radio de Marte(km)

	real       pmco2
	parameter  (pmco2 = 44.)

	real       pmo3p
	parameter  (pmo3p = 16.)

	real       pmo2
	parameter  (pmo2 = 32.)

	real       pmh2o2
	parameter  (pmh2o2 = 36.)

	real       cpco2
	parameter  (cpco2=20.0e6)

	real       cpo2
	parameter  (cpo2=0.9194e7)

	real       cpo3p
	parameter  (cpo3p=1.30e7)

	integer	   nreact
	parameter  (nreact=22)

	integer    nzmax
	parameter  (nzmax=210)


