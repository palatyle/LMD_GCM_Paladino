        subroutine rdf() 



#include "dimensions.h"
#include "microtab.h"
#include "clesphys.h"

        common/part/v,rayon,vrat,dr,dv
 
* declaration des variables communes
* ----------------------------------

        integer xnz,xnrad,xnztop
        integer li,lf,h
* nrad dans microtab.h
        real v(nrad),rayon(nrad),vrat,dr(nrad),dv(nrad)
        
* controles


       
* definition de la grille de rayon
* --------------------------------
        print*,'NOUVELLE GRILLE RAYON BASEE SUR 40 BINS'
        print*,'ATTENTION; TRAVAIL AVEC UN NOUVEAU RAYON'
        print*,'             DE MONOMER :               '
        cbase=nint(40./nrad)*1. 
        pi=3.1415926
        rayon(1)=13e-10*2**((cbase/2.-.5)/3.)   
        rayon(1)=13e-10*2.**.3333333333333    !<*****
	vrat=2.**cbase
	v(1)=4./3.*pi*rayon(1)**3

	do 9 i=1,nrad-1
	   rayon(i+1)=rayon(1)*vrat**(i/3.)
	   v(i+1)=v(1)*vrat**i
9	continue

	do 10 i=1,nrad
	  dv(i)=((vrat-1.)/(vrat+1.))*2.*v(i)
	  dr(i)=(2./(vrat+1))**(1./3.)*(vrat**(1./3.)-1.)*rayon(i)
10      continue

 



* parametres fractals : rf & df(h)
*------------------------------------

        rf0=0.066e-6      !Rayon monomere...OBLIGATOIRE meme en df=3!
                          !meme si dans ce cas, sa valeur n'a aucune
			  !importance
        do h=1,nrad
        rf(h)=rf0                             !<*********
        enddo 

        print*, rf(5),' METRES '             

        do h=1,nrad
          df(h)=3.        !Df pour petites particules
          if(rayon(h).ge.rf(h)) df(h)=df_GP      
        enddo 



        aknc=2.92         !<--------Df=3 
        aknc=6.86         !<--------Df=2 


        print*,'tcorrect=',tcorrect,' tx=',tx
        print*,'Df aerosols /1 a nqtot/'
        write(*,*) (rf(h),h=1,nrad)
        write(*,*) (df(h),h=1,nrad)
        write(*,*) (rayon(h),h=1,nrad)



        return
	end

