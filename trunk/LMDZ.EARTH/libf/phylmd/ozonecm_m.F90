! $Header$
module ozonecm_m

  IMPLICIT NONE

contains

  function ozonecm(rlat, paprs, rjour)

    ! The ozone climatology is based on an analytic formula which fits the
    ! Krueger and Mintzner (1976) profile, as well as the variations with
    ! altitude and latitude of the maximum ozone concentrations and the total
    ! column ozone concentration of Keating and Young (1986). The analytic
    ! formula have been established by J.-F. Royer (CRNM, Meteo France), who
    ! also provided us the code.

    ! A. J. Krueger and R. A. Minzner, A Mid-Latitude Ozone Model for the
    ! 1976 U.S. Standard Atmosphere, J. Geophys. Res., 81, 4477, (1976).

    ! Keating, G. M. and D. F. Young, 1985: Interim reference models for the
    ! middle atmosphere, Handbook for MAP, vol. 16, 205-229.

    USE dimphy, only: klon, klev
    use assert_m, only: assert

    REAL, INTENT (IN) :: rlat(:) ! (klon)
    REAL, INTENT (IN) :: paprs(:, :) ! (klon,klev+1)
    REAL, INTENT (IN) :: rjour

    REAL ozonecm(klon,klev)
    ! "ozonecm(j, k)" is the column-density of ozone in cell "(j, k)", that is
    ! between interface "k" and interface "k + 1", in kDU.

    ! Variables local to the procedure:

    REAL tozon ! equivalent pressure of ozone above interface "k", in Pa
    real pi, pl
    INTEGER i, k

    REAL field(klon,klev+1)
    ! "field(:, k)" is the column-density of ozone between interface
    ! "k" and the top of the atmosphere (interface "llm + 1"), in kDU.

    real, PARAMETER:: ps=101325.
    REAL, parameter:: an = 360., zo3q3 = 4E-8
    REAL, parameter:: dobson_unit = 2.1415E-5 ! in kg m-2
    REAL gms, zslat, zsint, zcost, z, ppm, qpm, a
    REAL asec, bsec, aprim, zo3a3

    !----------------------------------------------------------

    call assert((/size(rlat), size(paprs, 1)/) == klon, "ozonecm klon")
    call assert(size(paprs, 2) == klev + 1, "ozonecm klev")

    pi = 4. * atan(1.)
    DO k = 1, klev
       DO i = 1, klon
          zslat = sin(pi / 180. * rlat(i))
          zsint = sin(2 * pi * (rjour + 15.) / an)
          zcost = cos(2 * pi * (rjour + 15.) / an)
          z = 0.0531 + zsint * (-0.001595+0.009443*zslat) &
               + zcost * (-0.001344-0.00346*zslat) &
               + zslat**2 * (.056222 + zslat**2 &
               * (-.037609+.012248*zsint+.00521*zcost+.008890*zslat))
          zo3a3 = zo3q3/ps/2.
          z = z - zo3q3*ps
          gms = z
          ppm = 800. - (500.*zslat+150.*zcost)*zslat
          qpm = 1.74E-5 - (7.5E-6*zslat+1.7E-6*zcost)*zslat
          bsec = 2650. + 5000.*zslat**2
          a = 4.0*(bsec)**(3./2.)*(ppm)**(3./2.)*(1.0+(bsec/ps)**(3./2.))
          a = a/(bsec**(3./2.)+ppm**(3./2.))**2
          aprim = (2.666666*qpm*ppm-a*gms)/(1.0-a)
          aprim = amax1(0., aprim)
          asec = (gms-aprim)*(1.0+(bsec/ps)**(3./2.))
          asec = amax1(0.0, asec)
          aprim = gms - asec/(1.+(bsec/ps)**(3./2.))
          pl = paprs(i, k)
          tozon = aprim / (1. + 3. * (ppm / pl)**2) &
               + asec / (1. + (bsec / pl)**(3./2.)) + zo3a3 * pl * pl
          ! Convert from Pa to kDU:
          field(i, k) = tozon / 9.81 / dobson_unit / 1e3
       END DO
    END DO

    field(:,klev+1) = 0.
    forall (k = 1: klev) ozonecm(:,k) = field(:,k) - field(:,k+1)
    ozonecm = max(ozonecm, 1e-12)

  END function ozonecm

end module ozonecm_m
