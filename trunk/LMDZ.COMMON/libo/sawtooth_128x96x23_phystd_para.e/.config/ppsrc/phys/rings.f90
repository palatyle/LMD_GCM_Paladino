










SUBROUTINE rings(ngrid, declin, ptime, rad, flat, eclipse)
! Calculates Saturn's rings shadowing
! Includes rings opacities measured by Cassini/UVIS
! Authors: M. Sylvestre, M. Capderou, S. Guerlet, A. Spiga

    use comdiurn_h, only: sinlat, sinlon, coslat, coslon
    use geometry_mod, only: latitude ! (rad)
 
    implicit none   

    INTEGER, INTENT(IN) :: ngrid  ! horizontal grid dimension
    REAL, INTENT(IN) :: declin    ! latitude of the subsolar point
    REAL, INTENT(IN) :: ptime     ! UTC time in sol fraction : ptime=0.5 at noon
    REAL, INTENT(IN) :: rad       ! equatorial radius of the planet
    REAL, INTENT(IN) :: flat      ! flattening of the planet 
    REAL, DIMENSION(ngrid), INTENT(OUT) :: eclipse ! absorption of the light by the rings    
    
    REAL :: rpol   ! polar radius of the planet
    REAL :: e      ! shape excentricity of the planet : (1-e*e) = (1-f)*(1-f)    
    INTEGER, PARAMETER :: nb_a = 4 ! number of subdivisions of the A ring
    INTEGER, PARAMETER :: nb_b = 3 ! number of subdivisions of the B ring
    INTEGER, PARAMETER :: nb_c = 3 ! number of subdivisions of the C ring
    INTEGER, PARAMETER :: nb_ca = 2 ! number of subdivisions in the Cassini division
    INTEGER :: i

    ! arrays for the rings. TBD: dynamical?
    REAL, DIMENSION(nb_a) :: A_Rint ! internal radii of the subdivisions of the A ring 
    REAL, DIMENSION(nb_a) :: A_Rext ! external radii of the subdivisions of the A ring
    REAL, DIMENSION(nb_b) :: B_Rint ! internal radii of the subdivisions of the B ring
    REAL, DIMENSION(nb_b) :: B_Rext ! external radii of the subdivisions of the B ring 
    REAL, DIMENSION(nb_c) :: C_Rint ! internal radii of the subdivisions of the C ring
    REAL, DIMENSION(nb_c) :: C_Rext ! external radii of the subdivisions of the C ring 
    REAL, DIMENSION(nb_ca) :: Ca_Rint ! internal radii of the subdivisions of the Cassini Division
    REAL, DIMENSION(nb_ca) :: Ca_Rext ! external radii of the subdivisions of the Cassini Division

    ! Opacities of the rings : for each one we can give different opacities for each part
    REAL, DIMENSION(nb_a) :: tau_A ! opacity of the A ring
    REAL, DIMENSION(nb_b) :: tau_B ! opacity of the B ring
    REAL, DIMENSION(nb_c) :: tau_C ! opacity of the C ring
    REAL, DIMENSION(nb_ca) :: tau_Ca ! opacity of the Cassini Division 

    ! Parameters used to calculate if a point is under a ring subdivision's shadow
    REAL :: phi_S                             ! subsolar point longitude
    REAL, PARAMETER :: pi=acos(-1.0)    
    REAL, DIMENSION(:), ALLOCATABLE:: x, y, z ! cartesian coordinates of the points on the planet
    REAL :: xs, ys, zs                        ! cartesian coordinates of the points of the subsolar point
    REAL, DIMENSION(:), ALLOCATABLE :: k
    REAL, DIMENSION(:), ALLOCATABLE :: N      ! parameter to compute cartesian coordinates on a ellipsoidal planet
    REAL, DIMENSION(:), ALLOCATABLE :: r      ! distance at which the incident ray of sun crosses the equatorial plane
                                              ! measured from the center of the planet   
    REAL :: Ns                                ! (same for the subsolar point)
   
    ! equinox --> no shadow (AS: why is this needed?)
    if(declin .eq. 0.) then
        eclipse(:) = 0.
        return 
    endif 

! 1) INITIALIZATION

    ! Generic
    rpol = (1.- flat)*rad
    e = sqrt(2*flat - flat**2)
    ALLOCATE(x(ngrid))
    ALLOCATE(y(ngrid))
    ALLOCATE(z(ngrid))
    ALLOCATE(k(ngrid))
    ALLOCATE(N(ngrid))
    ALLOCATE(r(ngrid))
    eclipse(:) = 2000.

! Model of the rings with Cassini/UVIS opacities

    ! Size of the rings
    A_Rint(1) = 2.03*rad
    A_Rext(1) = 2.06*rad
    A_Rint(2) = 2.06*rad
    A_Rext(2) = 2.09*rad
    A_Rint(3) = 2.09*rad
    A_Rext(3) = 2.12*rad
    A_Rint(4) = 2.12*rad
    A_Rext(4) = 2.27*rad

    B_Rint(1) = 1.53*rad
    B_Rext(1) = 1.64*rad
    B_Rint(2) = 1.64*rad
    B_Rext(2) = 1.83*rad
    B_Rint(3) = 1.83*rad
    B_Rext(3) = 1.95*rad
    
    C_Rint(1) = 1.24*rad
    C_Rext(1) = 1.29*rad
    C_Rint(2) = 1.29*rad
    C_Rext(2) = 1.43*rad
    C_Rint(3) = 1.43*rad
    C_Rext(3) = 1.53*rad

    Ca_Rint(1) = 1.95*rad
    Ca_Rext(1) = 1.99*rad
    Ca_Rint(2) = 1.99*rad
    Ca_Rext(2) = 2.03*rad


    ! Opacities of the rings
    tau_A(1) = 1.24
    tau_A(2) = 0.81
    tau_A(3) = 0.67
    tau_A(4) = 0.58
                
    tau_B(1) = 1.29
    tau_B(2) = 5.13 
    tau_B(3) = 2.84 
    
    tau_C(1) = 0.06
    tau_C(2) = 0.10
    tau_C(3) = 0.14

    tau_Ca(1) = 0.06
    tau_Ca(2) = 0.24

    ! Convert to cartesian coordinates
    N(:) = rad / sqrt(1-(e**2)*sinlat(:)**2)
    x(:) = N(:)*coslat(:)*coslon(:)
    y(:) = N(:)*coslat(:)*sinlon(:)
    z(:) = N(:)*(1-e**2)*sinlat(:)

! 2) LOCATION OF THE SUBSOLAR POINT 
 
    ! subsolar longitude is deduced from time fraction ptime
    ! SG: the minus sign is important! ... otherwise subsolar point adopts a reverse rotation
    phi_S = -(ptime - 0.5)*2.*pi 
!    write(*,*) 'subsol point coords : ', declin*180./pi, phi_S*180./pi

    ! subsolar latitude is declin (declination of the sun)
    ! now convert in cartesian coordinates : 
    Ns = rad/sqrt(1-(e**2)*sin(declin)**2)
    xs = Ns*cos(declin)*cos(phi_S)
    ys = Ns*cos(declin)*sin(phi_S)
    zs = Ns*(1-e**2)*sin(declin)

! 3) WHERE DOES THE INCIDENT RAY OF SUN CROSS THE EQUATORIAL PLAN ?

    k(:) = -z(:)/zs
    r(:) = (k(:)*xs + x(:))**2 + (k(:)*ys + y(:))**2 
    r(:) = sqrt(r(:))

! 4) SO WHERE ARE THE SHADOW OF THESE RINGS ?

    ! Summer hemisphere is not under the shadow of the rings
    where(latitude(:)*declin .gt. 0.)
       eclipse(:) = 1000.
    end where

    ! No shadow of the rings by night
    where(x(:)*xs + y(:)*ys + z(:)*zs .lt. 0.)
       eclipse(:) = 1000.
    end where

    ! if the incident rays of sun cross a ring, there is a shadow
    do i=1, nb_A 
        where(r(:) .ge. A_Rint(i) .and. r(:) .le. A_Rext(i) .and. eclipse(:) .ne. 1000.)
            eclipse(:) = 1. - exp(-tau_A(i)/abs(sin(declin)))
        end where
    end do 

    do i=1, nb_B 
        where(r(:) .ge. B_Rint(i) .and. r(:) .le. B_Rext(i) .and. eclipse(:) .ne. 1000.)
            eclipse(:) = 1. - exp(-tau_B(i)/abs(sin(declin)))
        end where
    enddo
    
    do i=1, nb_C 
        where(r(:) .ge. C_Rint(i) .and. r(:) .le. C_Rext(i) .and. eclipse(:) .ne. 1000.)
            eclipse(:) = 1. - exp(-tau_C(i)/abs(sin(declin)))
        end where
    enddo

    do i=1, nb_ca
        where(r(:) .ge. Ca_Rint(i) .and. r(:) .le. Ca_Rext(i) .and. eclipse(:) .ne. 1000.)
            eclipse(:) = 1. - exp(-tau_Ca(i)/abs(sin(declin)))
        end where
    enddo

    ! At the other places and the excluded ones, eclipse is 0. 
    where(eclipse(:) .eq. 2000. .or. eclipse(:) .eq. 1000.)
        eclipse(:) = 0. 
    end where 

! 5) CLEAN THE PLACE
    DEALLOCATE(x)
    DEALLOCATE(y)
    DEALLOCATE(z)
    DEALLOCATE(k)
    DEALLOCATE(N)
    DEALLOCATE(r)

END SUBROUTINE rings
