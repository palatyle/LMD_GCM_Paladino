  module compo_hedin83_mod2
  !!05/2013 Laura Salmi
  !!03/2014 revision Gabriella Gilli
  !!Calcul des vmr pour CO2, CO, O, N et N2 en s'appuyant sur les tables donnees dans l'article Hedin (1983)
  implicit none      
  contains

        subroutine compo_hedin83_init2

        implicit none
#include "YOMCST.h"
#include "hedin.h"

	REAL :: alpha
        REAL :: mu(musize),z(zsize)
        REAL :: T(musize,zsize)
	REAL :: ntot(musize, zsize)
	REAL :: nco2(musize, zsize),nn2(musize, zsize)
	REAL :: nco(musize, zsize), no(musize, zsize)
	REAL :: nn(musize, zsize),nhe(musize, zsize)
	REAL :: sig_co2(zsize), sig_o(zsize), sig_n2(zsize),sig_co(zsize)
	integer :: i,j
	REAL factor_ox3p


! Initialisation des sza et de l'atlitude
        mu(1)=1.
        do i=2,musize-1
              mu(i)=mu(i-1)-1./9. 
        enddo
        mu(musize)=0.

        z(1)=100.
        do i=2,zsize
              z(i)=z(i-1)+5.
        enddo

!Lecture des variables (tables de Hedin (1983))

       open( 11, file='HIGHATM/noon.txt')
      read (11, *)
      do i=1,zsize
         read (11,*) T(1,i), nco2(1,i),no(1,i), nco(1,i), nn2(1,i),nn(1,i), nhe(1,i)
      enddo
      close (11)
       open( 12, file='HIGHATM/midnight.txt')
      read (12, *)
      do i=1,zsize
         read (12,*) T(musize,i), nco2(musize,i),no(musize,i), nco(musize,i), nn2(musize,i), nn(musize,i), nhe(musize,i)
      enddo
      close (12)

!Dependance en SZA

	do i=2,musize-1
          alpha=0.5*tanh(4.*(mu(i)-0.5))+0.5
	  do j=1,zsize 
		T(i,j)=alpha*T(1,j)+(1.-alpha)*T(musize,j)
		nco2(i,j)=alpha*nco2(1,j)+(1.-alpha)*nco2(musize,j)
	  	no(i,j)=alpha*no(1,j)+(1.-alpha)*no(musize,j)
	  	nco(i,j)=alpha*nco(1,j)+(1.-alpha)*nco(musize,j)
	  	nn2(i,j)=alpha*nn2(1,j)+(1.-alpha)*nn2(musize,j)
	  	nn(i,j)=alpha*nn(1,j)+(1.-alpha)*nn(musize,j)
	  	nhe(i,j)=alpha*nhe(1,j)+(1.-alpha)*nhe(musize,j)

	  enddo
         enddo

!! Test: effect of varying atomic oxygen abundances

         factor_ox3p = 2.
           DO i=1,musize
             DO j=1,zsize
               no(i,j)    =  factor_ox3p* no(i,j)
             ENDDO
           ENDDO






!! Conversion en volume mixture ratio

	   do i=1,musize
	   	  do j=1,zsize 
		  ntot(i,j)=nco2(i,j)+nn2(i,j)+nn(i,j)+nco(i,j)+no(i,j)+nhe(i,j)+2e-3*nco2(i,j) !for o2 and no si on change la proportion changer aussi dans euvheat
		  pres_hedin(i,j)=ntot(i,j)*1e6*RKBOL*T(i,j)
	 	 co2_hedin(i,j)=nco2(i,j)/ntot(i,j)
	 	 co_hedin(i,j)=nco(i,j)/ntot(i,j)
	 	 n2_hedin(i,j)=nn2(i,j)/ntot(i,j)
	 	 n_hedin(i,j)=nn(i,j)/ntot(i,j)
	 	 o_hedin(i,j)=no(i,j)/ntot(i,j)
	 	 mu_hedin(i)=mu(i)

	  enddo
         enddo

!            print*, pres_hedin(1,:)
!	    print*, ' '
!           print*, T(musize,:)  
!	    print*, ' '
!            print*, T(1,:)
!		print*, nn2(10,:)
!		print*, ' '
!		print*,nco2(10,:)
!		print*, ' '
!		print*,nco(10,:)
!		print*, ' '
!		print*,no(10,:)

!		print*, " "
!		print*, ntot(1,:)
!		print*, co2_hedin(:,10)
!		print*, mu_hedin(:)
!		print*, " "
!		print*, n2_hedin(1,:)
!		print*, " "
!		print*, o_hedin(1,:)
!		print*, " "
!		print*, co_hedin(1,:)
!		print*, ' '
!	 stop
		
       end subroutine compo_hedin83_init2



SUBROUTINE compo_hedin83_mod(pression,mu0,co2vmr_gcm,covmr_gcm,ovmr_gcm,n2vmr_gcm,nvmr_gcm)


!!Interpolation des profils de Hedin (1983) sur la grille du modele

	use dimphy
	implicit none
	include 'hedin.h'
	REAL, intent(in) :: pression(klon,klev), mu0(klon)
	integer :: i,k,iz,jmu,jz,z_ok, mu_ok
	REAL :: factp(klev)
	REAL, intent(out) :: co2vmr_gcm(klon,klev),covmr_gcm(klon,klev)
	REAL, intent(out) ::ovmr_gcm(klon,klev), n2vmr_gcm(klon,klev), nvmr_gcm(klon,klev)
	REAL :: ang0, ang_hedin
!        print*, pres_hedin(1,:)
!	print*, ' '
!	print*, pres_hedin(10,:)
!	stop
       do i=1,klon
	 ang0=acos(mu0(i))*180./3.1415
	 do jmu=1,musize
	    ang_hedin=acos(mu_hedin(jmu))*180./3.1415
         if (ang_hedin.le.ang0 ) then   
!          if (mu_hedin(jmu) .le. mu0(i) ) then  
              mu_ok = jmu	   
  	 endif
	 enddo 
!	print*, ang_hedin, ang0, mu_ok
!	STOP

	do k=1,klev
	 z_ok=2		
         do jz=2,zsize-1
           if (pres_hedin(mu_ok,jz).ge.pression(i,k)) then
                 z_ok = jz+1
           endif     
          enddo

        factp(k) = (log10(pression(i,k))-log10(pres_hedin(mu_ok,z_ok-1)))/(log10(pres_hedin(mu_ok,z_ok))-log10(pres_hedin(mu_ok,z_ok-1)))

      	
        ovmr_gcm(i,k) = (10.**(log10(o_hedin(mu_ok,z_ok))*factp(k)+log10(o_hedin(mu_ok,z_ok-1))*(1.-factp(k))))
        n2vmr_gcm(i,k) = 10.**(log10(n2_hedin(mu_ok,z_ok))*factp(k)+log10(n2_hedin(mu_ok,z_ok-1))*(1.-factp(k)))
        nvmr_gcm(i,k) = 10.**(log10(n_hedin(mu_ok,z_ok))*factp(k)+log10(n_hedin(mu_ok,z_ok-1))*(1.-factp(k)))
	covmr_gcm(i,k) = 10.**(log10(co_hedin(mu_ok,z_ok))*factp(k)+log10(co_hedin(mu_ok,z_ok-1))*(1.-factp(k)))
        co2vmr_gcm(i,k) =1.- (ovmr_gcm(i,k)+n2vmr_gcm(i,k)+covmr_gcm(i,k))
!        co2vmr_gcm(i,k) =1.- (ovmr_gcm(i,k)+n2vmr_gcm(i,k)+covmr_gcm(i,k)+nvmr_gcm(i,k))

	enddo	
       enddo
!      stop
	

END SUBROUTINE compo_hedin83_mod


  end module compo_hedin83_mod2

