
!***********************************************

!	param.par

!	Parameters for paramhr.f
!***********************************************

	integer    ninter
	parameter  (ninter=36)

	integer    nabs
	parameter  (nabs=13)
!	parameter  (nabs=5)

	integer    nz2
	parameter  (nz2=253)

	integer    ninter2
        parameter  (ninter2=16)

	real       kboltzman                  !cte Boltzman
	parameter  (kboltzman = 1.381e-16)
	
	real       n_avog                    !# de Avogadro
	parameter  (n_avog = 6.023e23)

	real       gg                        !cte gravitacion
	parameter  (gg = 6.67259e-8)

	real       masa                      !masa de Venus(g)
	parameter  (masa = 4.8676e27)

	real       radio
	parameter  (radio = 6051.8)           !radio de Venus(km)

	integer	   nreact
	parameter  (nreact=93)




