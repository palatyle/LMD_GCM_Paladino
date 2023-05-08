	 subroutine cmie(lambda,xn,xk,rad,qext,qsca,qabs,gg)

*  COPIE SUR B&H 
*  lambda: longueur d'onde
*  xn, xk: indice de refraction reel et imaginaire 
*  rad: rayone de la particule 
*  qsca,qext: coefficient de diffusion et d'absorption
*  gg: parametre d'assymetrie
 

       common/angle/theta  
       parameter (nang=451)
       complex refrel,s1(20000),s2(20000)
       real rad,lambda,s11(20000),theta(10000)
       real si1(20000),si2(20000)
      
c      rad=.525
c      lambda=.6328
       refrel=cmplx(xn,xk)
       x=2.*3.14159265*rad/lambda
       dang=1.570796327/dfloat(nang-1)

          call intmie2(x,refrel,nang,s1,s2,qext,qsca)

          surf=3.1415926*rad**2.*1e4
          qext=qext*surf
          qsca=qsca*surf
          qabs=qext-qsca
 

       do 355 j=1,2*nang-1
         s11(j)=cabs(s2(j))*cabs(s2(j))+cabs(s1(j))*cabs(s1(j))
         scarre=s11(j)
         s11(j)=s11(j)/(2*x**2.)
c       print*,scarre,theta(j)*180./3.1415926
 355   continue 
 

       tot1=0.
         g1=0.
         g2=0.
 
*  s11(j): fonction de phase 

       do 365 j=1,2*nang-2,2
        tot1=2*dang/6*(s11(j)*sin(theta(j))+s11(j+2)*sin(theta(j+2))
     &  +4*s11(j+1)*sin(theta(j+1)))*2*3.141592+tot1
 
        g1=2*dang/6*(s11(j)*sin(theta(j))+s11(j+2)*sin(theta(j+2))
     &  +4*s11(j+1)*sin(theta(j+1)))*2*3.141592+g1
 
        g2=2*dang/6*(s11(j)*sin(theta(j))*cos(theta(j))
     &  +s11(j+2)*sin(theta(j+2))*cos(theta(j+2))
     &  +4*s11(j+1)*sin(theta(j+1))*cos(theta(j+1)))*2*3.141592+g2
 
 365    continue
  
        gg=g2/g1
         do j=1,2*nang-1
c      print*,s11(j)/s11(1),s11(j),qsca
c      print*,si1(j),si2(j),theta(j) 

         enddo
c       print*,'******************************'
c      verification que integrale de la fonction de phase= qsca. 
c      write(*,*) tot,tot1,'int'
c      write(*,*) 'asymetrie:',g,s,gg
c       print*,'******************************'
        return
       end

c------------------------------------------------------------
        subroutine intmie2(x,refrel,nang,s1,s2,qext,qsca)

       common/angle/theta  
      
        real amu(10000),theta(10000),pi(10000)
        real tau(10000),pi0(10000),pi1(10000)
        complex d(300000),y,refrel,xi,xi0,xi1,an,bn,s1(20000),s2(20000)
        double precision psi0,psi1,psi,dn,dx

       dx=x            !dx en double precision ....            
       y=x*refrel


       xstop=x+4*x**.3333+2.
       nstop=xstop
       ymod=cabs(y)
       nmx=amax1(xstop,ymod)+15
c      print*,nmx,xstop,nstop,ymod
       dang=1.570796327/dfloat(nang-1)

         do 555 j=1,2*nang-1           
         theta(j)=(dfloat(j)-1.)*dang
 555     amu(j)=cos(theta(j))

       d(nmx)=cmplx(0.,0.)  
       nn=nmx-1
          
         do 120 n=1,nn
         rn=nmx-n+1
         d(nmx-n)=(rn/y)-(1./(d(nmx-n+1)+rn/y))
 120     continue
  
         do 666 j=1,nang
         pi0(j)=0.       ! fonction de legendre 
 666     pi1(j)=1.
   
       nn=2*nang-1
        
         do 777 j=1,nn
         s1(j)=cmplx(0.,0.)
 777     s2(j)=cmplx(0.,0.)


       psi0=dcos(dx)      !initialisation des fonctions de Bessel
       psi1=dsin(dx)
       chi0=-sin(x)
       chi1=cos(x)

       apsi0=psi0        !psi en double prec. et apsi en simple
       apsi1=psi1

       xi0=cmplx(apsi0,-chi0)
       xi1=cmplx(apsi1,-chi1)

       qsca=0.

       n=1

c      *************debut de l'iteration sur n *************
 200   dn=n
       rn=n
       fn=(2.*rn+1.)/(rn*(rn+1.))
      
       psi=(2.*dn-1.)*psi1/dx-psi0     ! calcul des fct de Bessel  
       chi=(2.*rn-1.)*chi1/x-chi0      ! relation recurrente a 2 niveaux
       apsi=psi
       xi=cmplx(apsi,-chi)

       an=(d(n)/refrel+rn/x)*apsi-apsi1
       an=an/((d(n)/refrel+rn/x)*xi-xi1)
       bn=(refrel*d(n)+rn/x)*apsi-apsi1
       bn=bn/((refrel*d(n)+rn/x)*xi-xi1)
      
       qsca=qsca+(2*rn+1.)*(cabs(an)*cabs(an)+cabs(bn)*cabs(bn))
  
c      print*,rn,cabs(an),cabs(bn) 
c     ***************debut de la boucle sur les angles******* 

       do 789 j=1,nang
       jj=2*nang-j
          pi(j)=pi1(j)                           !
          tau(j)=rn*amu(j)*pi(j)-(rn+1.)*pi0(j)  ! fonction de legendre
          s1(j)=s1(j)+fn*(an*pi(j)+bn*tau(j))
          s2(j)=s2(j)+fn*(an*tau(j)+bn*pi(j))
          p=(-1)**(n-1)
          t=(-1)**n

       if (j.eq.jj) goto 789
            s1(jj)=s1(jj)+fn*(an*pi(j)*p+bn*tau(j)*t)
            s2(jj)=s2(jj)+fn*(an*tau(j)*t+bn*pi(j)*p)
 789   continue

       psi0=psi1
       psi1=psi
       apsi1=psi1           ! double prec=simple

       chi0=chi1
       chi1=chi
       xi1=cmplx(apsi1,-chi1)

       n=n+1
       rn=n

       do 999 j=1,nang
        pi1(j)=((2.*rn-1.)/(rn-1.))*amu(j)*pi(j)
        pi1(j)=pi1(j)-rn*pi0(j)/(rn-1.)
 999    pi0(j)=pi(j)

       if (n-1-nstop) 200,300,300
 300   qsca=(2./(x*x))*qsca
       qext=(4./(x*x))*real(s1(1))
       
       return
       end
 
    
      
