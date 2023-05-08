         subroutine cfffv11(lambda,xn,xk,rad,DFS,NB,fsca,fext,fabs,fg)
*
*   NEW VERSION MARCH, 31, 1999.
*   WITHOUT MAXWELL-GARNETT APPROACH FOR IM{M} > 0.1
*   
 
 
       parameter (nang=91)

       complex mtest,refrel,s1(2000),s2(2000)
       real rad,lambda,s11(2000),theta(10000)
       real s110u(2000),s111u(2000)
       real s110d(2000),s111d(2000)
       real pol(2000),pp(2,2000),strS(2000)
       real NB
       real xxn

* COMMON WITH INTMIE

       common/angle11/theta  

* COMMON WITH CORRINT11

       common/pack1/cjup1,cju,cjup1_,cju_,xrap,xechj,xechjp1
       common/pack2/cjdp1,cjd,cjdp1_,cjd_,xku,xkd,jindex
       common/pack3/xexp5,xexp6

* COMMON WITH FFFV11

       DATA IND /0/


       if(nang*2.gt.2000) 
     & stop 'INCREASE THE SIZE OF THE ARRAYS IN CFFFV11' 

c      refrel=(xn,xk)
       refrel=cmplx(xn,xk)
       nindex=int(xn*100.)                
       x=2.*3.14159265*rad/lambda
       dang=1.570796327/dfloat(nang-1)
       pi=3.14159265
       NC=20 
          DF=1.99999
          alpha=1.1
          alphaS=(1.4-1.1)/(2.5-2.)*DFS-.1  !APPROXIMATED!
                                            !EXACT FOR D=2 AND 2.5!

******************************************************************
* STEP # 1
******************************************************************
* CALL THE ROUTINE THAT GIVE THE SCALING FACTOR'S PSI_s AND PSI_e
* FOR EACH BOUND OF THE INTERVAL  nX(j) < n*X < nX(j+1)  AND   
* xku < xk < xkd. 
* THESE FACTORS ARE GIVEN ASSUMING : N=20, DF=2 AND KNOWING n
*
* NB: HERE PSI ARE NOTED  cju,cjd,cjup1,cjdp1.....
******************************************************************

          call corrint11(x,xn,xk,NC)


******************************************************************
* STEP # 2
******************************************************************
* FOR EACH OF THE FOUR "POINTS" (nX,k), COMPUTE THE MIE CROSS-SECTIONS
* AND PHASE FUNCTION FOR A SPHERE AS ONE MONOMER. 
* DERIVE, THANKS TO THE PSI, THE MONOMER CROSS-SECTIONS
* COMPUTE ALSO THE ACTUAL MIE CROSS SECTION (I.E FOR THE ACTUAL 
* VALUES OF n,X and k)
******************************************************************

* STEP #2.1
******************************************
* if Mie & monomer in Rayleigh scattering: 
******************************************
* the phase function shape does not depends
* on x neither the PSI's: for each values of
* k, only one computation is to be done 
*****************************************

       IF (x*xn.le.0.085) THEN 


* pour k=xku
*            # x=xech(j)
*            # x=xech(j+1)

c      refrel=(xn,xku)
       refrel=cmplx(xn,xku)

       call intmie(x,refrel,nang,s1,s2,qext_,qsca_)
        qsca1u=alog10(qsca_*400.)
        qext1u=alog10(qext_*20.)
        qsca2u=alog10(qsca_*400.)
        qext2u=alog10(qext_*20.)

* pour k=xkd
*            # x=xech(j)
*            # x=xech(j+1)

c      refrel=(xn,xkd)
       refrel=cmplx(xn,xkd)

       call intmie(x,refrel,nang,s1,s2,qext_,qsca_)
        qsca1d=alog10(qsca_*400.)
        qext1d=alog10(qext_*20.)
        qsca2d=alog10(qsca_*400.)
        qext2d=alog10(qext_*20.)

        xrap=0.

* Mie exact
*

c      refrel=(xn,xk)
       refrel=cmplx(xn,xk)

       call intmie(x,refrel,nang,s1,s2,qext_,qsca_)
        do j=1,2*nang-1
        s11(j)=cabs(s2(j))*cabs(s2(j))+cabs(s1(j))*cabs(s1(j))
        s11(j)=s11(j)/(2.*pi*x**2.*qsca_)
        enddo


       ELSE

* STEP #2.2
************************************************
* if Mie & monomer in the range X about 1 to 5 *
************************************************

* pour k=xkd
*            # x=xech(j)
*            # x=xech(j+1)


        qsca1d=alog10(cjd)
        qext1d=alog10(cjd_)

        qsca2d=alog10(cjdp1)
        qext2d=alog10(cjdp1_)

* pour k=xku
*            # x=xech(j)
*            # x=xech(j+1)


        qsca1u=alog10(cju)
        qext1u=alog10(cju_)

        qsca2u=alog10(cjup1)
        qext2u=alog10(cjup1_)
 

* Mie exact:
*  

c      refrel=(xn,xk)
       refrel=cmplx(xn,xk)
       call intmie(x,refrel,nang,s1,s2,qext_,qsca_)
        do j=1,2*nang-1
        s11(j)=cabs(s2(j))*cabs(s2(j))+cabs(s1(j))*cabs(s1(j))
        enddo

       ENDIF

* THIS IS MIE:

          qabs_ =qext_ -qsca_ 

* THIS IS FRACTAL 20 MONOMERS
* d
        qabs2d=alog10(10.**qext2d-10.**qsca2d)
        qabs1d=alog10(10.**qext1d-10.**qsca1d)
* u
        qabs2u=alog10(10.**qext2u-10.**qsca2u)
        qabs1u=alog10(10.**qext1u-10.**qsca1u)

 
***********************************************************************
*       QUICK FIX FOR NEGATIVE ABSORPTION WHEN APPEARS [OUT OF n RANGE]
**********************************************************************
        GOTO 505
* This avoid a "crash" if one among four of the qsca is larger than qext
* Averaging the w=qsca/qext, and applied to the "wrong" qsca
* NOTE this occurs only if real refractive index > 2

        w1=10.**qsca2d/10.**qext2d 
        w2=10.**qsca1d/10.**qext1d 
        w3=10.**qsca2u/10.**qext2u 
        w4=10.**qsca1u/10.**qext1u 


        IF(10.**qext2d.lt.10.**qsca2d) 
     &   qsca2d=alog10((10.**qext2d)*(w2+w3+w4)/3.)
        IF(10.**qext1d.lt.10.**qsca1d)
     &   qsca1d=alog10((10.**qext1d)*(w1+w3+w4)/3.)
        IF(10.**qext2u.lt.10.**qsca2u)
     &   qsca2u=alog10((10.**qext2u)*(w2+w1+w4)/3.)
        IF(10.**qext1u.lt.10.**qsca1u)
     &   qsca1u=alog10((10.**qext1u)*(w2+w3+w1)/3.)

        qabs2d=alog10(10.**qext2d-10.**qsca2d)
        qabs1d=alog10(10.**qext1d-10.**qsca1d)
        qabs2u=alog10(10.**qext2u-10.**qsca2u)
        qabs1u=alog10(10.**qext1u-10.**qsca1u)
 505    continue

************************************************
* STEP # 5
************************************************
* INTERPOLATION OVER THE TWO VARIABLES n*X AND K 
************************************************

* X*n VARIABLE/ interpolation log-log. 

          sextu=((qext2u-qext1u)*xrap+qext1u)
          sscau=((qsca2u-qsca1u)*xrap+qsca1u)
          sabsu=((qabs2u-qabs1u)*xrap+qabs1u)
          sextd=((qext2d-qext1d)*xrap+qext1d)
          sscad=((qsca2d-qsca1d)*xrap+qsca1d)
          sabsd=((qabs2d-qabs1d)*xrap+qabs1d)

* K VARIABLE/ interpolation log-log
          xki=xk 
          if (xk.lt.1.e-2) xki=1.e-2
          if (xk.gt.1.) xki=1.
          rki=-alog10(xki)
          rku=-alog10(xku)
          rkd=-alog10(xkd)
          ssca=(sscau-sscad)*(rki-rkd)+sscad
          sext=(sextu-sextd)*(rki-rkd)+sextd
          sabs=(sabsu-sabsd)*(rki-rkd)+sabsd

          sscap=(sscau-sscad)*(rki-rkd)+sscad
          sextp=(sextu-sextd)*(rki-rkd)+sextd
          sabsp=(sabsu-sabsd)*(rki-rkd)+sabsd

* storage Mie down and up
c         refrel=(xn,1.)
          refrel=cmplx(xn,1.)
          call intmie(x,refrel,nang,s1,s2,qextmd_,qscamd_)
          qabsmd_=qextmd_-qscamd_

c         refrel=(xn,.1)
          refrel=cmplx(xn,.1)
          call intmie(x,refrel,nang,s1,s2,qextmu_,qscamu_)
          qabsmu_=qextmu_-qscamu_

******************************************************
* STEP # 5.1 :   
******************************************************
* SPECIAL CASE FOR LARGE k (>.1) : MAXWELL GARNETT
******************************************************


          if (xk.gt.0.1) then         ! THIS CONCERNS THE LARGE IMAGINARY REFRACTIVE INDEXES

*Prepare Maxwell Garnett approx. (see B&H) 
*-----------------------------------------

          xbulk=x*(NC*1.)**.3333333
          xsphe=x*(NC*1.)**.5
          frac=(xbulk/xsphe)**3.
              xkff=xk*frac
              xnff=1.+(xn-1.)*frac

* For small x, aggregate scattering= N**2 monmer scatt.:
* -----------------------------------------------------

          mtest=cmplx(xn,xk)
          call intmie(x,mtest,nang,s1,s2,qe,qs)
          sscaR=alog10(qs*NC*(xsphe/x)**2.)

 
*      For small x, aggregate abs. is proportionnal to N
* -----------------------------------------------------
          mtest=cmplx(xn,xk)
          call intmie(x,mtest,nang,s1,s2,qe,qs)
          sabsR=alog10((qe-qs)*20.)


          endif 


 
*****************************************************************
*     FINALLY, WE GET THE VALUES FOR AN AGGREGATE OF 20 MONOMERS
*          FOR THE ACTUAL X,k, AND n
*****************************************************************


         ssca=10.**ssca
         sext=10.**sext
         sabs=10.**sabs
 
* Sharp cross over between Rayleigh 
           xup=0.3             !            |<----------
           xdo=0.3             !  --------->|

          if (xk.gt.0.1) then       ! large imaginary refractive index

              if (x*xn.ge.xup) then  ! large size parameter
                 sext=sabs+ssca    
                 scal1=(10.**sabsd)/qabsmd_
                 scal2=(10.**sabsu)/qabsmu_
                 scalx=(alog10(xk)+1.)*(scal1-scal2)+scal2
                 sabs=qabs_*scalx
                 ssca=sext-sabs
               endif

              if (x*xn.le.xdo) then  ! PURE RAYLEIGH
                 sabs=10.**sabsR
                 ssca=10.**sscaR 
                 sext=sabs+ssca    
              endif

         
          else

                 sext=sabs+ssca   !self-consistent ext sca abs set

          endif


*****************************************************************
* STEP # 6
*****************************************************************
* NOW, USE USE THE N-LAWS TO COMPUTE THE Q(N) FROM Q(N=20) AND
* THE Q(MIE) [EQUATIONS (9) AND (10)].
*****************************************************************

         g0=0.
         g1=0.
         g2=0.

       surf=rad**2.*3.1415926*1.e18


* STEP # 6.1
*****************************************************************
* FIRST; COMPUTE THE STRUCTURE TERM THAT IS ONLY USED TO COMPUTE
* THE ASYMETRY PARAMETER THAT DEPENDS ON THE SHAPE OF THE PHASE
* FUNCTION.
*****************************************************************


       do 367 j=1,2*nang-1
         z=sin(theta(j)/2.)
         xx_=alphaS*x*(NB*1.)**(1./DFS)
         call structure(NB,xx_,z,str,DFS)
         if (z.eq.0.) str=(NB*1.)**2.
         strS(j)=str
  367   continue


       do 365 j=1,2*nang-2,2

        g1=2*dang/6*(s11(j)*sin(theta(j))*strS(j)
     &              +s11(j+2)*sin(theta(j+2))*strS(j+2)
     &            +4*s11(j+1)*sin(theta(j+1))*strS(j+1))
     &                          *2*3.141592+g1

        g2=2*dang/6*(s11(j)*sin(theta(j))*cos(theta(j))*strS(j)
     &        +s11(j+2)*sin(theta(j+2))*cos(theta(j+2))*strS(j+2)
     &      +4*s11(j+1)*sin(theta(j+1))*cos(theta(j+1))*strS(j+1))
     &           *2.*3.141592+g2

 365    continue


* STEP # 6.2
*****************************************************************
* USE THE EQUATIONS (10) AND (9)  FOR MEDIUM RANGE OR THE RAYLEIGH
* RANGE. 
*****************************************************************


* 1: FOR MEDIUM X RANGE: 1< X < 10 to 50  [EQ 10]
* 2: FOR VERY LARGE X RANGE:  X >> 50     [EQ 11]
*    N^a LOG(N^b) IS REPLACED BY A PURE N^a LAW. THIS OCCURS
*    WHEN N>50 FOR ABSORPTION AND WHEN N>25/x FOR SCATTERING,
*    AND CORRESPOND TO THE GEOMETRIC RANGE. THE EXPONENT a IS
*    SET BY EMPIRIC CONSIDERATIONS (SEE PAPER).

* XEXP5 AND XEXP6 ARE THE EXPONENT FACTOR USED IN EQ 10
* XEXP1 AND XEXP2 ARE THE EXPONENT FACTOR USED VERY LARGE BEHAVIOUR

         xkf=1.
         if (xk.le.1.e-2) xkf=xk/1.e-2


        rho=(NC*1.)**xexp5
        rho_=(NB*1.)**xexp5

* EQ 10 RESULT FOR ABSORPTION:
*-----------------------------

        sabs1=(sabs*xkf-qabs_)/alog10(20.)*alog10(NB*1.)
     &                        *rho_/rho+qabs_


* NB: IF NB>NTA OR NB>NTS , EFFICIENCY FACTOR FOLLOW A POWER LAW
*     BEHAVIOUR. ELSE, EQUATION 10 IS VALIDE.

        NTA=50

       if(NB.gt.NTA) then
       rhox=(NTA*1.)**xexp5

* COMPUTE sabs1 FOR NTA TO START THE POWER LAW:
        sabs1=(sabs*xkf-qabs_)/alog10(20.)*alog10(NTA*1.)
     &  *rhox/rho+qabs_

* COMPUTE xexp1 EXPONENT FACTOR: 
                        xexp1=.85
        if(x*xn.lt.2.2) xexp1=(.85-.96)*1.42*(x*xn-1.5)+.96
        if(x*xn.lt.1.5) xexp1=.96

* COMPUTE sabs1 VALUE FOR LARGE X:
        sabs1=sabs1*(NB*1./(NTA*1.))**xexp1
       endif

* HERE SABS1 MEDIUM RANGE IS KNOWN

        rho=(NC*1.)**xexp6
        rho_=(NB*1.)**xexp6


* EQ 10 RESULT FOR SCATTERING:
*-----------------------------

        ssca1=(ssca-qsca_)/alog10(20.)*alog10(NB*1.)
     &        *rho_/rho+qsca_

        NTS=int((20./x)**(1./.65))

* NB: IF NB>NTA OR NB>NTS , EFFICIENCY FACTOR FOLLOW A POWER LAW
*     BEHAVIOUR. ELSE, EQUATION 10 IS VALIDE.

       if(NB.gt.NTS) then
        rhox=(NTS*1.)**xexp6

* COMPUTE sabs1 FOR NTS TO START THE POWER LAW:
        ssca1=(ssca-qsca_)/alog10(20.)*alog10(NTS*1.)
     &  *rhox/rho+qsca_

* COMPUTE xexp2 EXPONENT FACTOR: 
                     xexp2=.95
        if(x.lt.1.0) xexp2=(.95-1.1)*2*(x-.5)+1.1

* COMPUTE sabs1 VALUE FOR LARGE X:
        ssca1=ssca1*(NB*1./(NTS*1.))**xexp2 
       endif

* HERE SSCA1 MEDIUM RANGE IS KNOWN


        sext1=ssca1+sabs1





* 3: FOR THE RAYLEIGH RANGE [EQ 9]

        ssca2=(alog10(ssca)-alog10(qsca_))/alog10(20.)
     &            *alog10(NB*1.)+alog10(qsca_)
        sext2=(alog10(sext)-alog10(qext_))/alog10(20.)
     &            *alog10(NB*1.)+alog10(qext_)
        sabs2=(alog10(sabs*xkf)-alog10(qabs_))/alog10(20.)
     &            *alog10(NB*1.)+alog10(qabs_)



************************************************
* STEP # 7
************************************************
* THE CROSS OVER BETWEEN RAYLEYGH AND MEDIUM RANGE
************************************************

        rho_=(NB*1.)**.66666
* 1< X < 50.
        sext1=alog10(sext1/rho_)
        ssca1=alog10(ssca1/rho_)
        sabs1=alog10(sabs1/rho_)
* RAYLEIGH
        sext2=alog10((10.**sext2)/rho_)
        ssca2=alog10((10.**ssca2)/rho_)
        sabs2=alog10((10.**sabs2)/rho_)



        


* 7.1: SCATTERING AND EXTINCTION  + CROSS OVER
************************************************


* RAYLEIGH RANGE:
            fsca=10.**(ssca2)
            fext=10.**(sext2)
            fg=g2/g1

         xx=.5*alpha*x*(NB*1.)**.33333
         xup=10.
         xdo=.1 

         f=alog10(xx/xdo)/alog10(xup/xdo)
         g=1.-f


* MEDIUM X
        if (xx.le.xup.and.xx.ge.xdo) then
           fsca=10.**(f*ssca1+g*ssca2)
           fext=10.**(f*sext1+g*sext2)
           fg=g2/g1
        endif

* LARGE X
        if (xx.gt.xup) then
           fsca=10.**(ssca1)
           fext=10.**(sext1)
           fg=g2/g1
        endif

* 7.2: ABSORPTION + CROSS OVER
************************************************

* RAYLEIGH RANGE

         fabs=10.**(sabs2)

         xx=.5*alpha*x*(NB*1.)**.5
         f=alog10(xx/xdo)/alog10(xup/xdo)
         g=1.-f

* MEDIUM X
        if (xx.le.xup.and.xx.ge.xdo) then
           fabs=10.**(f*sabs1+g*sabs2)
        endif

* LARGE X
        if (xx.gt.xup) then
           fabs=10.**(sabs1)
        endif

         fext=fsca+fabs
           

***********************************************
* STEP # 8
***********************************************
*   DISPLAY  THE RESULTS : PHASE FUNCTION OR 
*  EFFICIENCY COEFFICIENTS 
***********************************************

         IF (IND.eq.1) THEN

* 8.1:    CROSS SECTIONS  &  PHASE FUNCTION
***********************************************

         tot =0.
         scoef=NB**(2./3.)*rad**2.*3.1415926

        do 366 j=1,2*nang-1
 
         pp(1,j)=theta(j)*180./3.1415926
         pp(2,j)=s11(j)*2*3.1415936*strS(j)*1./g1
 366    continue
* 
        ELSE

* 8.2:       EFFICIENCY COEFFICIENTS
***********************************************

       scoef=NB**(2./3.)*rad**2.*3.1415926  ! METERS^2
        fsca=fsca*scoef
        fext=fext*scoef
        fabs=fabs*scoef

c new fashion...

        fabs=(qext_-qsca_)*NB/(NB*1.)**(2./3.)*scoef
        fext=fsca+fabs
 

        ENDIF

        return
       end

c------------------------------------------------------------
        subroutine intmie(x,refrel,nang,s1,s2,
     &                     qext,qsca)
**********************************************************
* THIS ROUTINE COMES FROM BORHEN AND HUFFMAN BOOK:
* "ABSORPTION AND SCATTERING OF LIGHT BY SMALL PARTICLES"
*  WILEY INTERSCIENCE PUBLICATION, 1983
**********************************************************
       common/angle11/theta  
      
        real amu(10000),theta(10000),pi(10000)
        real tau(10000),pi0(10000),pi1(10000)
        complex d(300000),y,refrel,xi,xi0,xi1,an,bn,s1(20000),s2(20000)
        complex s_(20000)
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
c      d(nmx)=(0.,0.)  
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
c        s_(j)=(0.,0.)
c        s1(j)=(0.,0.)
c777     s2(j)=(0.,0.)
         s_(j)=cmplx(0.,0.)
         s1(j)=cmplx(0.,0.)
 777     s2(j)=cmplx(0.,0.)
       psi0=dcos(dx)      !initialisation des fonctions de Bessel
       psi1=dsin(dx)
       chi0=-sin(x)
       chi1=cos(x)
       apsi0=psi0        !psi en double prec. et apsi en simple
       apsi1=psi1
c      xi0=(apsi0,-chi0)
c      xi1=(apsi1,-chi1)
       xi0=cmplx(apsi0,-chi0)
       xi1=cmplx(apsi1,-chi1)
       qsca=0.
       qsca_=0.
       n=1
c      *************debut de l'iteration sur n *************
 200   dn=n
       rn=n
       fn=(2.*rn+1.)/(rn*(rn+1.))
      
       psi=(2.*dn-1.)*psi1/dx-psi0     ! calcul des fct de Bessel  
       chi=(2.*rn-1.)*chi1/x-chi0      ! relation recurrente a 2 niveaux
       apsi=psi
c      xi=(apsi,-chi)
       xi=cmplx(apsi,-chi)
       an=(d(n)/refrel+rn/x)*apsi-apsi1
       an=an/((d(n)/refrel+rn/x)*xi-xi1)
       bn=(refrel*d(n)+rn/x)*apsi-apsi1
       bn=bn/((refrel*d(n)+rn/x)*xi-xi1)
       qsca=qsca+(2*rn+1.)*(cabs(an)*cabs(an)+cabs(bn)*cabs(bn))
  
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
c      xi1=(apsi1,-chi1)
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
       qabs=qext-qsca
       
       
       return
       end
 
    
       subroutine corrint11(x,xn,xk,NB)
      parameter (nech=28)
      parameter (nex=6)
      real xech(nech)
      real A0(2*nech), B0(2*nech), C0(2*nech), D0(2*nech)
      real A1(2*nech), B1(2*nech), C1(2*nech), D1(2*nech)
      real A2(2*nech), B2(2*nech), C2(2*nech), D2(2*nech)
      real A3(2*nech), B3(2*nech), C3(2*nech), D3(2*nech)
      real x,xn,xk,correct,correct_
      integer NB,ifirst
      real ww1(nex),ww2(nex)
      real xx1(nex),xx2(nex)
      real yy1(nex),yy2(nex)
      real zz1(nex),zz2(nex)
      real xldo(nex)
      real xexp5,xexp6
       common/pack1/cjup1_,cju_,cjup1,cju,xrapl,xechj,xechjp1
       common/pack2/cjdp1_,cjd_,cjdp1,cjd,xku,xkd,jindex
       common/pack3/xexp5,xexp6
*      * ATTENTION: ORDRE INVERSE DANS cfffpf, cfffcs,...  *
      DATA xech/0.05,0.1,0.25,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,
     & 1.4,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0,4.25,
     & 4.5,4.75,5.0/
      DATA A0/ .1117E+00, .2284E+00, .4473E+00,-.3351E+00,-.7323E+00,
     &        -.8070E+00,-.5968E+00,-.2956E+00, .1171E-01, .2515E+00,
     &         .3607E+00, .3653E+00, .2658E+00, .1145E+00,-.1505E+00,
     &         .1982E-01, .8496E-01,-.5869E-01,-.7623E-01,-.5410E-02,
     &        -.2850E-01,-.5501E-01,-.5267E-01,-.2317E-01,-.3344E-01,
     &        -.3936E-01,-.3094E-01,-.3099E-01,
     &         .6136E-03, .7863E-02, .9823E-01, .6704E-01,-.8727E-01,
     &        -.1741E+00,-.1402E+00,-.3270E-01, .1074E+00, .2397E+00,
     &         .3254E+00, .3654E+00, .3495E+00, .2869E+00, .7161E-01,
     &         .6914E-01, .1319E+00, .7578E-01, .2340E-01, .4031E-01,
     &         .4202E-01, .2560E-01, .1679E-01, .2965E-01, .3108E-01,
     &         .2983E-01, .3023E-01, .2988E-01/
      DATA B0/-.5034E+00,-.1033E+01,-.2186E+01, .3670E+00, .2014E+01,
     &         .2646E+01, .2279E+01, .1496E+01, .5713E+00,-.2137E+00,
     &        -.6203E+00,-.7014E+00,-.4274E+00, .3641E-01, .8890E+00,
     &         .3026E+00, .4121E-01, .4774E+00, .5016E+00, .2491E+00,
     &         .3097E+00, .3754E+00, .3536E+00, .2500E+00, .2734E+00,
     &         .2838E+00, .2453E+00, .2429E+00,
     &        -.2493E-02,-.3224E-01,-.4251E+00,-.5115E+00, .7309E-01,
     &         .5198E+00, .5940E+00, .3982E+00, .3812E-01,-.3514E+00,
     &        -.6344E+00,-.7949E+00,-.7802E+00,-.6081E+00, .5837E-01,
     &         .4640E-01,-.1963E+00,-.4932E-01, .9157E-01, .1454E-01,
     &        -.6767E-02, .3216E-01, .4963E-01,-.9453E-04,-.1184E-01,
     &        -.1405E-01,-.2189E-01,-.2395E-01/
      DATA C0/ .6106E+00, .1257E+01, .2904E+01, .1547E+01, .1531E+00,
     &        -.5216E+00,-.4001E+00, .1061E+00, .7939E+00, .1421E+01,
     &         .1776E+01, .1883E+01, .1693E+01, .1334E+01, .6370E+00,
     &         .1107E+01, .1329E+01, .9812E+00, .9691E+00, .1171E+01,
     &         .1119E+01, .1070E+01, .1086E+01, .1164E+01, .1145E+01,
     &         .1135E+01, .1167E+01, .1164E+01,
     &         .2614E-02, .3429E-01, .4919E+00, .1051E+01, .6861E+00,
     &         .3474E+00, .2567E+00, .3720E+00, .6355E+00, .9464E+00,
     &         .1195E+01, .1361E+01, .1388E+01, .1283E+01, .7888E+00,
     &         .8176E+00, .1038E+01, .9433E+00, .8490E+00, .9234E+00,
     &         .9487E+00, .9256E+00, .9171E+00, .9601E+00, .9731E+00,
     &         .9774E+00, .9868E+00, .9892E+00/
      DATA A1/ .8917E-02, .1683E-01, .1060E-01,-.2180E+00,-.4334E+00,
     &        -.6977E+00,-.9510E+00,-.1114E+01,-.1154E+01,-.1186E+01,
     &        -.1170E+01,-.1184E+01,-.1214E+01,-.1243E+01,-.1324E+01,
     &        -.1570E+01,-.1304E+01,-.7215E+00,-.9031E+00,-.1164E+01,
     &        -.5731E+00,-.4235E+00,-.8438E+00,-.3892E+00, .3259E+00,
     &         .2680E+00, .6864E+00, .1542E+01,
     &        -.6903E-04,-.1081E-02,-.2849E-01,-.2557E+00,-.4515E+00,
     &        -.6935E+00,-.9313E+00,-.1098E+01,-.1165E+01,-.1222E+01,
     &        -.1239E+01,-.1285E+01,-.1353E+01,-.1432E+01,-.1606E+01,
     &        -.1809E+01,-.1613E+01,-.1139E+01,-.1174E+01,-.1452E+01,
     &        -.1075E+01,-.7739E+00,-.1122E+01,-.8658E+00, .8014E-02,
     &         .2460E+00, .7393E+00, .1802E+01/
      DATA B1/-.4265E-01,-.8221E-01,-.9239E-01, .7328E+00, .1571E+01,
     &         .2661E+01, .3805E+01, .4705E+01, .5192E+01, .5586E+01,
     &         .5713E+01, .5844E+01, .5940E+01, .5939E+01, .5617E+01,
     &         .6022E+01, .5165E+01, .2970E+01, .3196E+01, .4005E+01,
     &         .1932E+01, .1164E+01, .2469E+01, .8862E+00,-.1572E+01,
     &        -.1359E+01,-.2718E+01,-.5456E+01,
     &         .2345E-03, .3756E-02, .1058E+00, .9941E+00, .1787E+01,
     &         .2813E+01, .3909E+01, .4823E+01, .5400E+01, .5881E+01,
     &         .6131E+01, .6373E+01, .6598E+01, .6769E+01, .6816E+01,
     &         .7097E+01, .6380E+01, .4542E+01, .4294E+01, .5046E+01,
     &         .3642E+01, .2361E+01, .3303E+01, .2337E+01,-.6819E+00,
     &        -.1528E+01,-.3110E+01,-.6542E+01/
      DATA C1/ .5537E-01, .1094E+00, .1933E+00,-.3295E+00,-.9664E+00,
     &        -.1836E+01,-.2791E+01,-.3579E+01,-.4033E+01,-.4384E+01,
     &        -.4467E+01,-.4506E+01,-.4480E+01,-.4338E+01,-.3594E+01,
     &        -.3592E+01,-.2857E+01,-.9067E+00,-.8569E+00,-.1483E+01,
     &         .2426E+00, .1010E+01,-.2250E-01, .1272E+01, .3286E+01,
     &         .3064E+01, .4112E+01, .6241E+01,
     &        -.1499E-03,-.2538E-02,-.7950E-01,-.7859E+00,-.1437E+01,
     &        -.2302E+01,-.3255E+01,-.4084E+01,-.4636E+01,-.5081E+01,
     &        -.5287E+01,-.5434E+01,-.5522E+01,-.5524E+01,-.5105E+01,
     &        -.5006E+01,-.4323E+01,-.2649E+01,-.2219E+01,-.2726E+01,
     &        -.1514E+01,-.3276E+00,-.9927E+00,-.1745E+00, .2306E+01,
     &         .2984E+01, .4190E+01, .6871E+01/
      DATA A2/ .8210E-03, .6438E-03,-.2567E-01,-.2701E+00,-.4978E+00,
     &        -.8074E+00,-.1155E+01,-.1458E+01,-.1649E+01,-.1731E+01,
     &        -.1683E+01,-.1675E+01,-.1646E+01,-.1544E+01,-.1234E+01,
     &        -.1855E+01,-.1307E+01,-.2446E+00,-.5417E+00,-.1624E+01,
     &        -.2091E+00, .2351E+00,-.1552E+01,-.5738E+00, .1313E+01,
     &         .4078E-01, .1116E+01, .5153E+01,
     &        -.7454E-04,-.1152E-02,-.2981E-01,-.2758E+00,-.5027E+00,
     &        -.8111E+00,-.1158E+01,-.1463E+01,-.1658E+01,-.1747E+01,
     &        -.1709E+01,-.1708E+01,-.1687E+01,-.1596E+01,-.1308E+01,
     &        -.1910E+01,-.1371E+01,-.3399E+00,-.6009E+00,-.1651E+01,
     &        -.3374E+00, .1413E+00,-.1566E+01,-.7407E+00, .1151E+01,
     &         .9063E-01, .1134E+01, .5181E+01/
      DATA B2/-.4020E-02,-.4571E-02, .9034E-01, .1031E+01, .1929E+01,
     &         .3185E+01, .4662E+01, .6060E+01, .7090E+01, .7687E+01,
     &         .7703E+01, .7758E+01, .7644E+01, .7176E+01, .5232E+01,
     &         .6690E+01, .4972E+01, .1198E+01, .1616E+01, .5180E+01,
     &         .4740E+00,-.1502E+01, .4302E+01, .9918E+00,-.5496E+01,
     &        -.1222E+01,-.4838E+01,-.1792E+02,
     &         .2569E-03, .4046E-02, .1111E+00, .1065E+01, .1964E+01,
     &         .3219E+01, .4696E+01, .6101E+01, .7148E+01, .7768E+01,
     &         .7824E+01, .7907E+01, .7820E+01, .7389E+01, .5527E+01,
     &         .6926E+01, .5219E+01, .1545E+01, .1852E+01, .5278E+01,
     &         .8977E+00,-.1184E+01, .4303E+01, .1494E+01,-.5022E+01,
     &        -.1536E+01,-.4991E+01,-.1814E+02/
      DATA C2/ .5370E-02, .8378E-02,-.5663E-01,-.7911E+00,-.1516E+01,
     &        -.2549E+01,-.3786E+01,-.4980E+01,-.5871E+01,-.6372E+01,
     &        -.6307E+01,-.6228E+01,-.5961E+01,-.5361E+01,-.3045E+01,
     &        -.3777E+01,-.2316E+01, .8933E+00, .8813E+00,-.1974E+01,
     &         .1810E+01, .3680E+01,-.9658E+00, .1690E+01, .7053E+01,
     &         .3513E+01, .6392E+01, .1679E+02,
     &        -.1729E-03,-.2839E-02,-.8489E-01,-.8456E+00,-.1578E+01,
     &        -.2617E+01,-.3861E+01,-.5067E+01,-.5977E+01,-.6502E+01,
     &        -.6476E+01,-.6425E+01,-.6184E+01,-.5615E+01,-.3369E+01,
     &        -.4053E+01,-.2592E+01, .5416E+00, .6159E+00,-.2112E+01,
     &         .1417E+01, .3365E+01,-.9937E+00, .1261E+01, .6647E+01,
     &         .3797E+01, .6509E+01, .1696E+02/
      DATA itime /0/
         data xldo/0.1,0.5,1.0,2.0,4.0,8.0/
         data  xx1/0.66,0.70,0.60,0.50,0.50,0.59/
         data  xx2/0.66,0.70,0.81,0.60,0.59,0.62/
         data  yy1/0.66,0.70,0.52,0.44,0.59,0.59/
         data  yy2/0.66,0.70,0.74,0.55,0.62,0.62/
         data  zz1/0.66,0.70,0.42,0.48,0.55,0.59/
         data  zz2/0.66,0.70,0.65,0.55,0.63,0.62/

      save ifirst

         if (ifirst.eq.0) then
         print*,' IFIRST', ifirst
          do i=1,nech
             xech(i)=xech(i)*1.7
          enddo
          ifirst=1
         endif

**    1: compute the exponent of the law N^a LOG N
**    with the index n and the size parameter x
****************************************************
           xexp5=0.66666
           xexp6=0.66666
******     INTERPOLATION WITH THE VALUE OF REAL REFR. INDEX 
           do i=1,nex
           if(xn.lt.1.7) then
           ww1(i)=(yy1(i)-xx1(i))/.3*(xn-1.4)+xx1(i)
           ww2(i)=(yy2(i)-xx2(i))/.3*(xn-1.4)+xx2(i)
           endif
           if(xn.ge.1.7) then
           ww1(i)=(zz1(i)-yy1(i))/.3*(xn-1.7)+yy1(i)
           ww2(i)=(zz2(i)-yy2(i))/.3*(xn-1.7)+yy2(i)
           endif
           enddo
******     INTERPOLATION  WITH THE SIZE PARAMETER
********  XEXP5 AND XEXP6 ARE THE EXPONENT TO USE IN N^a LOGN LAW
           do i=2,nex
           if(x.lt.xldo(i).and.x.ge.xldo(i-1)) then
            rap=xldo(i)-xldo(i-1)
            xexp5=(ww1(i)-ww1(i-1))/rap*(x-xldo(i-1))+ww1(i-1)
            xexp6=(ww2(i)-ww2(i-1))/rap*(x-xldo(i-1))+ww2(i-1)
           endif
           enddo
           if(x.ge.xldo(nex)) then
              xexp5=ww1(nex)
              xexp6=ww2(nex)
           endif
 
         xxn=xn*x

**   location in the defined range of n
********************************************
**   Out of the range
            xni=xn
            dxn=0.
            if (xn.lt.1.4) then 
                 dxn=(xn-1.4)
                 xn=1.4
            endif
            if (xn.gt.2.0) then 
                 dxn=(xn-2.0)
                 xn=2.0
            endif
           
**   location in the defined ranges of n*X
********************************************
**   Out of the range
            if (xxn.ge.xech(nech)) xxn=xech(nech)*.99 
**   Rayleigh range
            if (xxn.lt.xech(1)) xxn=xech(1)
              do i=1,nech
             if(xech(i).le.xxn) j=i
              enddo
**  Location inside the slab j to j+1
            xechj=xech(j)
            xechjp1=xech(j+1)
            jindex=j
            xrap=(xxn-xech(j))/(xech(j+1)-xech(j))
            xrapl=(alog10(xxn)-alog10(xech(j)))
     &        /(alog10(xech(j+1))-alog10(xech(j)))
             if(xxn.gt.1.7) xrapl=xrap
***************************************************
** Calculation of the (parabolic) coefficients   **
***************************************************
          xki=xk
          if (xk.lt.1.e-2) xki=1.e-2
          if (xk.gt.1.) xki=1.
          rki=-alog10(xki)                                                      
****   Computation of the parabolic coefficients
**     For index j
*       f(x0+dx)=f(x0)+[df/dx](x0)*dx
* with f=ax^2+bx+c  ---> f'=[f-c+ax^2]/x=2ax+b
*
       if(rki.gt.2.) then 
        cjd =A2(j)     *xn**2.+B2(j)     *xn+C2(j)
        cjd_=A2(j+nech)*xn**2.+B2(j+nech)*xn+C2(j+nech) 
        cju =A2(j)     *xn**2.+B2(j)     *xn+C2(j)
        cju_=A2(j+nech)*xn**2.+B2(j+nech)*xn+C2(j+nech) 
* neutre si dxn=0  :  
        cjd =cjd +(cjd -C2(j)     +A2(j)     *xn**2.)/xn*dxn
        cjd_=cjd_+(cjd_-C2(j+nech)+A2(j+nech)*xn**2.)/xn*dxn
        cju =cju +(cju -C2(j)     +A2(j)     *xn**2.)/xn*dxn
        cju_=cju_+(cju_-C2(j+nech)+A2(j+nech)*xn**2.)/xn*dxn
         xku=1.e-2
         xkd=1.e-2
       endif
       if(rki.le.2.) then
        cjd =A1(j)     *xn**2.+B1(j)     *xn+C1(j)
        cjd_=A1(j+nech)*xn**2.+B1(j+nech)*xn+C1(j+nech) 
        cju =A2(j)     *xn**2.+B2(j)     *xn+C2(j)
        cju_=A2(j+nech)*xn**2.+B2(j+nech)*xn+C2(j+nech) 
* neutre si dxn=0  :
        cjd =cjd +(cjd -C1(j)     +A1(j)     *xn**2.)/xn*dxn
        cjd_=cjd_+(cjd_-C1(j+nech)+A1(j+nech)*xn**2.)/xn*dxn
        cju =cju +(cju -C2(j)     +A2(j)     *xn**2.)/xn*dxn
        cju_=cju_+(cju_-C2(j+nech)+A2(j+nech)*xn**2.)/xn*dxn
         xku=1.e-2
         xkd=1.e-1
       endif
       if(rki.le.1.) then 
        cju =A1(j)     *xn**2.+B1(j)     *xn+C1(j)
        cju_=A1(j+nech)*xn**2.+B1(j+nech)*xn+C1(j+nech) 
        cjd =A0(j)     *xn**2.+B0(j)     *xn+C0(j)
        cjd_=A0(j+nech)*xn**2.+B0(j+nech)*xn+C0(j+nech) 
* neutre si dxn=0  :
        cjd =cjd +(cjd -C0(j)     +A0(j)     *xn**2.)/xn*dxn
        cjd_=cjd_+(cjd_-C0(j+nech)+A0(j+nech)*xn**2.)/xn*dxn
        cju =cju +(cju -C1(j)     +A1(j)     *xn**2.)/xn*dxn
        cju_=cju_+(cju_-C1(j+nech)+A1(j+nech)*xn**2.)/xn*dxn
         xku=1.e-1
         xkd=1.e-0
       endif
 
 
**     For index j+1
       if(rki.gt.2.) then 
        cjdp1 =A2(j+1)     *xn**2.+B2(j+1)     *xn+C2(j+1)
        cjdp1_=A2(j+1+nech)*xn**2.+B2(j+1+nech)*xn+C2(j+1+nech) 
        cjup1 =A2(j+1)     *xn**2.+B2(j+1)     *xn+C2(j+1)
        cjup1_=A2(j+1+nech)*xn**2.+B2(j+1+nech)*xn+C2(j+1+nech) 
* neutre si dxn=0  :
        cjdp1 =cjdp1 +(cjdp1 -C2(j+1)     +A2(j+1)     *xn**2.)/xn*dxn
        cjdp1_=cjdp1_+(cjdp1_-C2(j+1+nech)+A2(j+1+nech)*xn**2.)/xn*dxn
        cjup1 =cjup1 +(cjup1 -C2(j+1)     +A2(j+1)     *xn**2.)/xn*dxn
        cjup1_=cjup1_+(cjup1_-C2(j+1+nech)+A2(j+1+nech)*xn**2.)/xn*dxn
         xku=1.e-2
         xkd=1.e-2
       endif
       if(rki.le.2.) then
        cjdp1 =A1(j+1)     *xn**2.+B1(j+1)     *xn+C1(j+1)
        cjdp1_=A1(j+1+nech)*xn**2.+B1(j+1+nech)*xn+C1(j+1+nech) 
        cjup1 =A2(j+1)     *xn**2.+B2(j+1)     *xn+C2(j+1)
        cjup1_=A2(j+1+nech)*xn**2.+B2(j+1+nech)*xn+C2(j+1+nech) 
* neutre si dxn=0  :
        cjdp1 =cjdp1 +(cjdp1 -C1(j+1)     +A1(j+1)     *xn**2.)/xn*dxn
        cjdp1_=cjdp1_+(cjdp1_-C1(j+1+nech)+A1(j+1+nech)*xn**2.)/xn*dxn
        cjup1 =cjup1 +(cjup1 -C2(j+1)     +A2(j+1)     *xn**2.)/xn*dxn
        cjup1_=cjup1_+(cjup1_-C2(j+1+nech)+A2(j+1+nech)*xn**2.)/xn*dxn
         xku=1.e-2
         xkd=1.e-1
       endif
       if(rki.le.1.) then 
        cjup1 =A1(j+1)     *xn**2.+B1(j+1)     *xn+C1(j+1)
        cjup1_=A1(j+1+nech)*xn**2.+B1(j+1+nech)*xn+C1(j+1+nech) 
        cjdp1 =A0(j+1)     *xn**2.+B0(j+1)     *xn+C0(j+1)
        cjdp1_=A0(j+1+nech)*xn**2.+B0(j+1+nech)*xn+C0(j+1+nech) 
* neutre si dxn=0  :
        cjdp1 =cjdp1 +(cjdp1 -C0(j+1)     +A0(j+1)     *xn**2.)/xn*dxn
        cjdp1_=cjdp1_+(cjdp1_-C0(j+1+nech)+A0(j+1+nech)*xn**2.)/xn*dxn
        cjup1 =cjup1 +(cjup1 -C1(j+1)     +A1(j+1)     *xn**2.)/xn*dxn
        cjup1_=cjup1_+(cjup1_-C1(j+1+nech)+A1(j+1+nech)*xn**2.)/xn*dxn
         xku=1.e-1
         xkd=1.e-0
       endif
 
 
**     Computation of the monomer-factor
 
         cjup1 =cjup1*20.   
         cjup1_=cjup1_*20.
         cju   =cju*20.
         cju_  =cju_*20.
         cjdp1 =cjdp1*20.   
         cjdp1_=cjdp1_*20.
         cjd   =cjd*20.
         cjd_  =cjd_*20.
         xn=xni              ! restitution de xn

         return
         end
        
        subroutine structure(NB,X,Z,STRUCT,DF)
        implicit real (a-h,o-z)
c       integer NB
        real NB
        real X,Z,D
 
        D=DF
        if (DF.eq.2) D=2.0001
        if (z.eq.0.) z=1.e-4
 
         STRUCT=1.
         NOSTRUCT=0       ! if asymetry parameters are not needed
                          ! just skip the computation: save your time!
        if (NOSTRUCT.eq.1)  goto 102 
        u0=5.
        pi=3.1415926
* If convergence test in on (end of the loop):
        xacc=1.e-3
* Else, computation is done once: accuracy is generally about 1%
        
* The structure factor is computed in order to evaluate the asymetry
* parameter (not for cross sections). We need to compute the integral 
* of the following function:
*
*        sin(2XZu)exp(-1/2u**2) for u between 0 and 5.
*
*   where  X,Z are provided through the subroutine calling
*   A=4*pi (normalization factor for D=2 --> Botet et al., 1995)
*
* And STRUCT is given as:
*
* STRUCT=N*[1+(N-1)*2*pi/(A*X*Z) INTEGRAL[sin(2XZu)exp(-1/2u**2)du]
*
* The number of oscillations for sin(2XZu) between 0 and U is:
*  n=UXZ/pi.... let integer (simpson integration).
         
         A=-5.026*D+22.618
         A=A/(2.*pi)                    !<---- here is A/(2*PI)
         nplt=6*int(5.*Z*X/3.1415926+1.)
*        nplt=int(5.*Z*X/3.1415926+1.)
*        This is  the number of periods for the sinus in the range 
*        u E [0,5]. Integration is done with 6 points per period.
*        Accuracy is about 1% on the final result STRUCT. 
*
*        The minimum value for nplt is set to 17 to do computation 
*        in the Rayleigh range ( when Z*X  reaches 0)
         STRUCT=1.e-10
 101     if ((nplt/2)*1..eq.nplt*.5) nplt=nplt+1  !---> odd
         if (nplt.lt.17) nplt=17   
         dint=u0/(nplt-1)
         STRUCT_OLD=STRUCT
         STRUCT=0.
       
***   EXTENDED  SIMPSON INTEGRATION
         iint=0
         ucr=iint*dint
         um1=sin(2.*x*z*ucr)*exp(-.5*ucr**D)*ucr**(D-2.)
          STRUCT=STRUCT+um1
         iint=nplt-1
         ucr=iint*dint
         um1=sin(2.*x*z*ucr)*exp(-.5*ucr**D)*ucr**(D-2.)
          STRUCT=STRUCT+um1
         
         do iint=1,nplt-2,2
         ucr=iint*dint
         um1=4.*sin(2.*x*z*ucr)*exp(-.5*ucr**D)*ucr**(D-2.)
          STRUCT=STRUCT+um1
         enddo 
   
         do iint=2,nplt-3,2
         ucr=iint*dint
         um1=2.*sin(2.*x*z*ucr)*exp(-.5*ucr**D)*ucr**(D-2.)
          STRUCT=STRUCT+um1
         enddo 
         STRUCT=dint/3.*STRUCT
         ERR=abs(STRUCT_OLD-STRUCT)/STRUCT  
         nplt=int(nplt*2)
C        UNCOMMENT THE IF STATEMENT TO
c        SET ON THE CONVERGENCE TEST: 
c        if (ERR.gt.xacc) GOTO 101 
         STRUCT=(STRUCT/(x*z*a)*(NB-1)+1.)*NB

 102     continue        
         return 
         end

      
