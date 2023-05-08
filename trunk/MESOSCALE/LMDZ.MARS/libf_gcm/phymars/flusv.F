      SUBROUTINE flusv(KDLON,nsf,n,omega,g,tau,emis,bh,bsol,fah,fdh)
      IMPLICIT NONE
c.......................................................................
c
c  calcul des flux ascendant et descendant aux interfaces entre n couches
c  * dans l'infrarouge 
c  * B est une fonction lineaire de $\tau$ a l'interieur de chaque couche
c  * le B du sol peut etre different de celui qui correspond au profil
c    de la n-ieme couche
c  * l'hypothese est une hypothese a deux flux isotropes sur chaque 
c    hemisphere ("hemispheric constant") + "source function technique"
c    (voir Toon et al. 1988)
c  * le flux descendant en haut de l'atmosphere est nul
c  * les couches sont numerotees du haut de l'atmosphere vers le sol
c
c  in :   * KDLON      ---> dimension de vectorisation
c         * nsf        ---> nsf=0 ==> "hemispheric constant"
c                           nsf>0 ==> "hemispheric constant" + "source function"
c         * n          ---> nombre de couches
c         * omega(i)   ---> single scattering albedo pour la i-eme couche
c         * g(i)       ---> asymmetry parameter pour la i-eme couche
c         * tau(i)     ---> epaisseur optique de la i-eme couche
c         * emis       ---> emissivite du sol
c         * bh(i)      ---> luminance du corps noir en haut de la i-eme
c                           couche, bh(n+1) pour la valeur au sol qui
c                           correspond au profil de la n-ieme couche
c         * bsol       ---> luminance du corps noir au sol
c
c  out :  * fah(i)     ---> flux ascendant en haut de la i-eme couche,
c                           fah(n+1) pour le sol 
c         * fdh(i)     ---> flux descendant en haut de la i-eme couche,
c                           fdh(n+1) pour le sol 
c
c.......................................................................
c  include des dimensions locales
c 
#include "dimensions.h"
#include "dimphys.h"
#include "dimradmars.h"
c.......................................................................
c  declaration des arguments 
c
      INTEGER KDLON,nsf,n
      REAL omega(NDLO2,n),g(NDLO2,n),tau(NDLO2,n),emis(NDLO2)
     &,bh(NDLO2,n+1),bsol(NDLO2),fah(NDLO2,n+1),fdh(NDLO2,n+1)
c.......................................................................
c  declaration des variables locales
c 
      REAL pi
      PARAMETER (pi=3.141592653589793E+0)
      INTEGER iv,i,j
      REAL beta,gama1,gama2,amu1,grgama,b0,b1
      REAL a(NDLON,4*nlayermx),b(NDLON,4*nlayermx)
     &    ,d(NDLON,4*nlayermx),e(NDLON,4*nlayermx)
     &    ,y(NDLON,4*nlayermx)
     &    ,alambda(NDLON,2*nlayermx)
     &    ,e1(NDLON,2*nlayermx),e2(NDLON,2*nlayermx)
     &    ,e3(NDLON,2*nlayermx),e4(NDLON,2*nlayermx)
     &    ,cah(NDLON,2*nlayermx),cab(NDLON,2*nlayermx)
     &    ,cdh(NDLON,2*nlayermx),cdb(NDLON,2*nlayermx)
      REAL grg(NDLON,2*nlayermx),grh(NDLON,2*nlayermx)
     &    ,grj(NDLON,2*nlayermx),grk(NDLON,2*nlayermx)
     &    ,alpha1(NDLON,2*nlayermx),alpha2(NDLON,2*nlayermx)
     &    ,sigma1(NDLON,2*nlayermx),sigma2(NDLON,2*nlayermx)
      INTEGER nq
      PARAMETER (nq=8)
      REAL x(nq),w(nq),gri(NDLON,nq)
      DATA x/1.9855071751231860E-2 , 0.1016667612931866E+0
     &     , 0.2372337950418355E+0 , 0.4082826787521751E+0
     &     , 0.5917173212478250E+0 , 0.7627662049581645E+0
     &     , 0.8983332387068134E+0 , 0.9801449282487682E+0/
      DATA w/5.0614268145185310E-2 , 0.1111905172266872E+0
     &     , 0.1568533229389437E+0 , 0.1813418916891810E+0
     &     , 0.1813418916891810E+0 , 0.1568533229389437E+0
     &     , 0.1111905172266872E+0 , 5.0614268145185310E-2/
c....................................................................... 
c
c....................................................................... 
      do 10001 i=1,n
                                                   do 99999 iv=1,KDLON
      beta=(1.E+0-g(iv,i))/2.E+0
      gama1=2.E+0*(1.E+0-omega(iv,i)*(1.E+0-beta))
      gama2=2.E+0*omega(iv,i)*beta
      amu1=5.E-1
      alambda(iv,i)=sqrt(gama1**2-gama2**2)
      grgama=(gama1-alambda(iv,i))/gama2
c
c ici une petite bidouille : si l'epaisseur optique d'une couche
c est trop faible, $dB \over d\tau$ devient tres grand et le schema
c s'ecroule. dans ces cas, on fait l'hypothese de couche isotherme.
c
      if (tau(iv,i).gt.1.E-3) then
      b0=bh(iv,i)
      b1=(bh(iv,i+1)-b0)/tau(iv,i)
      else
      b0=(bh(iv,i)+bh(iv,i+1))/2.E+0
      b1=0.E+0
      endif
c
      e1(iv,i)=1.E+0+grgama*exp(-alambda(iv,i)*tau(iv,i))
      e2(iv,i)=1.E+0-grgama*exp(-alambda(iv,i)*tau(iv,i))
      e3(iv,i)=grgama+exp(-alambda(iv,i)*tau(iv,i))
      e4(iv,i)=grgama-exp(-alambda(iv,i)*tau(iv,i))
      cah(iv,i)=2.E+0*pi*amu1*(b0+b1/(gama1+gama2))
      cab(iv,i)=2.E+0*pi*amu1*(b0+b1*(tau(iv,i)+1.E+0/(gama1+gama2)))
      cdh(iv,i)=2.E+0*pi*amu1*(b0-b1/(gama1+gama2))
      cdb(iv,i)=2.E+0*pi*amu1*(b0+b1*(tau(iv,i)-1.E+0/(gama1+gama2)))  
c
      grg(iv,i)=(1.E+0/amu1-alambda(iv,i))
      grh(iv,i)=grgama*(alambda(iv,i)+1.E+0/amu1)
      grj(iv,i)=grh(iv,i)
      grk(iv,i)=grg(iv,i)
      alpha1(iv,i)=2.E+0*pi*(b0+b1*(1.E+0/(gama1+gama2)-amu1))
      alpha2(iv,i)=2.E+0*pi*b1
      sigma1(iv,i)=2.E+0*pi*(b0-b1*(1.E+0/(gama1+gama2)-amu1))
      sigma2(iv,i)=alpha2(iv,i)
c
99999                                              continue
10001 continue
c.......................................................................
                                                   do 99998 iv=1,KDLON
      a(iv,1)=0.E+0
      b(iv,1)=e1(iv,1)
      d(iv,1)=-e2(iv,1)
      e(iv,1)=-cdh(iv,1)
99998                                              continue
c
      do 10002 i=1,n-1 
      j=2*i+1
                                                   do 99997 iv=1,KDLON
      a(iv,j)=e2(iv,i)*e3(iv,i)-e4(iv,i)*e1(iv,i)
      b(iv,j)=e1(iv,i)*e1(iv,i+1)-e3(iv,i)*e3(iv,i+1)
      d(iv,j)=e3(iv,i)*e4(iv,i+1)-e1(iv,i)*e2(iv,i+1)
      e(iv,j)=e3(iv,i)*(cah(iv,i+1)-cab(iv,i))
     &       +e1(iv,i)*(cdb(iv,i)-cdh(iv,i+1))
99997                                              continue
10002 continue
c
      do 10003 i=1,n-1
      j=2*i
                                                   do 99996 iv=1,KDLON
      a(iv,j)=e2(iv,i+1)*e1(iv,i)-e3(iv,i)*e4(iv,i+1)
      b(iv,j)=e2(iv,i)*e2(iv,i+1)-e4(iv,i)*e4(iv,i+1)
      d(iv,j)=e1(iv,i+1)*e4(iv,i+1)-e2(iv,i+1)*e3(iv,i+1)
      e(iv,j)=e2(iv,i+1)*(cah(iv,i+1)-cab(iv,i))
     &       +e4(iv,i+1)*(cdb(iv,i)-cdh(iv,i+1))
99996                                              continue
10003 continue
c     
      j=2*n
                                                   do 99995 iv=1,KDLON
      a(iv,j)=e1(iv,n)-(1.E+0-emis(iv))*e3(iv,n)
      b(iv,j)=e2(iv,n)-(1.E+0-emis(iv))*e4(iv,n)
      d(iv,j)=0.E+0
      e(iv,j)=emis(iv)*pi*bsol(iv)-cab(iv,n)+(1.E+0-emis(iv))*cdb(iv,n)
99995                                              continue
c.......................................................................
      call sys3v(KDLON,2*n,a,b,d,e,y)
c....................................................................... 
      do 10004 i=1,n
                                                   do 99994 iv=1,KDLON
      grg(iv,i)=grg(iv,i)*(y(iv,2*i-1)+y(iv,2*i))
      grh(iv,i)=grh(iv,i)*(y(iv,2*i-1)-y(iv,2*i))
      grj(iv,i)=grj(iv,i)*(y(iv,2*i-1)+y(iv,2*i))
      grk(iv,i)=grk(iv,i)*(y(iv,2*i-1)-y(iv,2*i))
99994                                              continue
10004 continue
c.......................................................................
c les valeurs de flux "hemispheric constant"
c 
      IF (nsf.eq.0) THEN
      do 10005 i=1,n
                                                   do 99993 iv=1,KDLON
      fah(iv,i)=e3(iv,i)*y(iv,2*i-1)-e4(iv,i)*y(iv,2*i)+cah(iv,i)
      fdh(iv,i)=e1(iv,i)*y(iv,2*i-1)-e2(iv,i)*y(iv,2*i)+cdh(iv,i)
99993                                              continue
10005 continue
                                                   do 99992 iv=1,KDLON
      fah(iv,n+1)=e1(iv,n)*y(iv,2*n-1)+e2(iv,n)*y(iv,2*n)+cab(iv,n)
      fdh(iv,n+1)=e3(iv,n)*y(iv,2*n-1)+e4(iv,n)*y(iv,2*n)+cdb(iv,n)
99992                                              continue
      GOTO 10100
      ENDIF 
c....................................................................... 
c passage a "source function" 
c 
c on applique une quadrature de dimension nq
c x est le vecteur des \mu de la quadrature
c w est le vecteur des poids correspondants
c nq est fixe en parameter
c x et w sont fixes en data
c 
c....................................................................... 
c on part d'en haut et on descent selon les nq angles pour calculer
c tous les flux descendants 
c
      do 10006 j=1,nq 
                                                   do 99991 iv=1,KDLON
      gri(iv,j)=0.E+0
99991                                              continue
10006 continue
                                                   do 99990 iv=1,KDLON
      fdh(iv,1)=0.E+0
99990                                              continue
      do 10007 i=1,n
      do 10008 j=1,nq
                                                   do 99989 iv=1,KDLON
      gri(iv,j)=gri(iv,j)*exp(-tau(iv,i)/x(j))
     &+grj(iv,i)/(alambda(iv,i)*x(j)+1.E+0)
     &*(1.E+0-exp(-tau(iv,i)*(alambda(iv,i)+1.E+0/x(j))))
     &+grk(iv,i)/(alambda(iv,i)*x(j)-1.E+0)
     &*(exp(-tau(iv,i)/x(j))-exp(-tau(iv,i)*alambda(iv,i)))
     &+sigma1(iv,i)*(1.E+0-exp(-tau(iv,i)/x(j)))
     &+sigma2(iv,i)*(x(j)*exp(-tau(iv,i)/x(j))+tau(iv,i)-x(j))
99989                                              continue
10008 continue
                                                   do 99988 iv=1,KDLON
      fdh(iv,i+1)=0.E+0
99988                                              continue
      do 10009 j=1,nq 
                                                   do 99987 iv=1,KDLON
      fdh(iv,i+1)=fdh(iv,i+1)+w(j)*x(j)*gri(iv,j)
99987                                              continue
10009 continue
10007 continue
c.......................................................................
c on applique la condition de reflexion a sol 
c
                                                   do 99986 iv=1,KDLON
      fah(iv,n+1)=(1.E+0-emis(iv))*fdh(iv,n+1)+pi*emis(iv)*bsol(iv)
99986                                              continue
      do 10010 j=1,nq
                                                   do 99985 iv=1,KDLON
      gri(iv,j)=2.E+0*fah(iv,n+1)
99985                                              continue
10010 continue 
c.......................................................................
c on remonte pour calculer tous les flux ascendants 
c  
      do 10011 i=n,1,-1 
      do 10012 j=1,nq 
                                                   do 99984 iv=1,KDLON
      gri(iv,j)=gri(iv,j)*exp(-tau(iv,i)/x(j))
     &+grg(iv,i)/(alambda(iv,i)*x(j)-1.E+0)
     &*(exp(-tau(iv,i)/x(j))-exp(-tau(iv,i)*alambda(iv,i)))
     &+grh(iv,i)/(alambda(iv,i)*x(j)+1.E+0)
     &*(1.E+0-exp(-tau(iv,i)*(alambda(iv,i)+1.E+0/x(j))))
     &+alpha1(iv,i)*(1.E+0-exp(-tau(iv,i)/x(j)))
     &+alpha2(iv,i)*(x(j)-(tau(iv,i)+x(j))*exp(-tau(iv,i)/x(j)))
99984                                              continue
10012 continue
                                                   do 99983 iv=1,KDLON
      fah(iv,i)=0.E+0
99983                                              continue
      do 10013 j=1,nq 
                                                   do 99982 iv=1,KDLON
      fah(iv,i)=fah(iv,i)+w(j)*x(j)*gri(iv,j)
99982                                              continue
10013 continue
10011 continue
c.......................................................................
10100 continue      
c.......................................................................
c
c.......................................................................
      RETURN
      END

c ***************************************************************


      SUBROUTINE sys3v(KDLON,n,a,b,d,e,y)
      IMPLICIT NONE
c.......................................................................
c
c  inversion d'un systeme lineaire tridiagonal
c
c  |  b1  d1               |   | y1 |   | e1 |
c  |  a2  b2  d2           |   | y2 |   | e2 |
c  |      a3  b3  d3       | * | y3 | = | e3 |
c  |             ....      |   |    |   |    |
c  |              an  bn   |   | yn |   | en |
c
c  in :   * KDLON        --> dimension de vectorisation
c         * n            --> dimension du systeme
c         * a,b,d,e      --> voir dessin
c
c  out :  * y            --> voir dessin
c
c.......................................................................
c  include des dimensions locales
c 
#include "dimensions.h"
#include "dimphys.h"
#include "dimradmars.h"
c.......................................................................
c  declaration des arguments 
c
      INTEGER KDLON,n
      REAL a(NDLO2,n),b(NDLO2,n),d(NDLO2,n),e(NDLO2,n),y(NDLO2,n)
c.......................................................................
c  declaration des variables locales
c 
      INTEGER iv,i
      REAL as(NDLON,4*nlayermx),ds(NDLON,4*nlayermx)
     &    ,x(NDLON,4*nlayermx)
c.......................................................................
c
c.......................................................................
                                                   do 99999 iv=1,KDLON
      as(iv,n)=a(iv,n)/b(iv,n)
      ds(iv,n)=e(iv,n)/b(iv,n)
99999                                              continue
      do 10001 i=n-1,1,-1
                                                   do 99998 iv=1,KDLON
      x(iv,i)=1.E+0/(b(iv,i)-d(iv,i)*as(iv,i+1))
      as(iv,i)=a(iv,i)*x(iv,i)
      ds(iv,i)=(e(iv,i)-d(iv,i)*ds(iv,i+1))*x(iv,i)
99998                                              continue
10001 continue 
                                                   do 99997 iv=1,KDLON
      y(iv,1)=ds(iv,1)
99997                                              continue
      do 10002 i=2,n 
                                                   do 99996 iv=1,KDLON
      y(iv,i)=ds(iv,i)-as(iv,i)*y(iv,i-1)
99996                                              continue
10002 continue
c.......................................................................
c
c.......................................................................
      RETURN
      END 
