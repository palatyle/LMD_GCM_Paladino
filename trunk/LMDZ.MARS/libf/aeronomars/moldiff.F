      subroutine moldiff(ngrid,nlayer,nq,
     &                   pplay,pplev,pt,pdt,pq,pdq,ptimestep,
     &                   zzlay,pdteuv,pdtconduc,pdqdiff)

      use tracer_mod, only: igcm_co2, igcm_co, igcm_o, igcm_o1d,
     &                      igcm_o2, igcm_o3, igcm_h, igcm_h2, igcm_oh,
     &                      igcm_ho2, igcm_h2o2, igcm_n2, igcm_ar,
     &                      igcm_h2o_vap, mmol
      use conc_mod, only: rnew, mmean
      USE comcstfi_h
      implicit none

c
c Input/Output
c
      integer,intent(in) :: ngrid ! number of atmospheric columns
      integer,intent(in) :: nlayer ! number of atmospheric layers
      integer,intent(in) :: nq ! number of advected tracers
      real ptimestep
      real pplay(ngrid,nlayer)
      real zzlay(ngrid,nlayer)
      real pplev(ngrid,nlayer+1)
      real pq(ngrid,nlayer,nq)
      real pdq(ngrid,nlayer,nq)
      real pt(ngrid,nlayer)
      real pdt(ngrid,nlayer)
      real pdteuv(ngrid,nlayer)
      real pdtconduc(ngrid,nlayer)
      real pdqdiff(ngrid,nlayer,nq)
c
c Local
c

      integer,parameter :: ncompmoldiff = 14
       real hco2(ncompmoldiff),ho

      integer ig,nz,l,n,nn,iq
      real del1,del2, tmean ,dalfinvdz, d
      real hh,dcoef,dcoef1,ptfac, ntot, dens, dens2, dens3
      real hp(nlayer)
      real tt(nlayer)
      real qq(nlayer,ncompmoldiff)
      real dmmeandz(nlayer)
      real qnew(nlayer,ncompmoldiff)
      real zlocal(nlayer)
      real alf(ncompmoldiff-1,ncompmoldiff-1)
      real alfinv(nlayer,ncompmoldiff-1,ncompmoldiff-1)
      real indx(ncompmoldiff-1)
      real b(nlayer,ncompmoldiff-1)
      real y(ncompmoldiff-1,ncompmoldiff-1)
      real aa(nlayer,ncompmoldiff-1,ncompmoldiff-1)
      real bb(nlayer,ncompmoldiff-1,ncompmoldiff-1)
      real cc(nlayer,ncompmoldiff-1,ncompmoldiff-1)
      real atri(nlayer-2)
      real btri(nlayer-2)
      real ctri(nlayer-2)
      real rtri(nlayer-2)
      real qtri(nlayer-2)
      real alfdiag(ncompmoldiff-1)
      real wi(ncompmoldiff), flux(ncompmoldiff), pote

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     tracer numbering in the molecular diffusion
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Atomic oxygen must always be the LAST species of the list as
c it is the dominant species at high altitudes. 
      integer,parameter :: i_co   = 1
      integer,parameter :: i_n2   = 2
      integer,parameter :: i_o2   = 3
      integer,parameter :: i_co2  = 4
      integer,parameter :: i_h2   = 5
      integer,parameter :: i_h    = 6
      integer,parameter :: i_oh   = 7
      integer,parameter :: i_ho2  = 8
      integer,parameter :: i_h2o  = 9
      integer,parameter :: i_h2o2 = 10
      integer,parameter :: i_o1d  = 11
      integer,parameter :: i_o3   = 12
      integer,parameter :: i_ar   = 13
      integer,parameter :: i_o    = 14

! Tracer indexes in the GCM:
      integer,save :: g_co2=0
      integer,save :: g_co=0
      integer,save :: g_o=0
      integer,save :: g_o1d=0
      integer,save :: g_o2=0
      integer,save :: g_o3=0
      integer,save :: g_h=0
      integer,save :: g_h2=0
      integer,save :: g_oh=0
      integer,save :: g_ho2=0
      integer,save :: g_h2o2=0
      integer,save :: g_n2=0
      integer,save :: g_ar=0
      integer,save :: g_h2o=0

      integer,save :: gcmind(ncompmoldiff) ! array of GCM indexes
      integer ierr

      logical,save :: firstcall=.true.
      real abfac(ncompmoldiff)
      real,save :: dij(ncompmoldiff,ncompmoldiff)

! Initializations at first call
      if (firstcall) then
        call moldiffcoeff(dij)
        print*,'MOLDIFF  EXO'
        
        ! identify the indexes of the tracers we'll need
        g_co2=igcm_co2
        if (g_co2.eq.0) then
          write(*,*) "moldiff: Error; no CO2 tracer !!!"
          stop
        endif
        g_co=igcm_co
        if (g_co.eq.0) then
          write(*,*) "moldiff: Error; no CO tracer !!!"
          stop
        endif
        g_o=igcm_o
        if (g_o.eq.0) then
          write(*,*) "moldiff: Error; no O tracer !!!"
          stop
        endif
        g_o1d=igcm_o1d
        if (g_o1d.eq.0) then
          write(*,*) "moldiff: Error; no O1D tracer !!!"
          stop
        endif
        g_o2=igcm_o2
        if (g_o2.eq.0) then
          write(*,*) "moldiff: Error; no O2 tracer !!!"
          stop
        endif
        g_o3=igcm_o3
        if (g_o3.eq.0) then
          write(*,*) "moldiff: Error; no O3 tracer !!!"
          stop
        endif
        g_h=igcm_h
        if (g_h.eq.0) then
          write(*,*) "moldiff: Error; no H tracer !!!"
          stop
        endif
        g_h2=igcm_h2
        if (g_h2.eq.0) then
          write(*,*) "moldiff: Error; no H2 tracer !!!"
          stop
        endif
        g_oh=igcm_oh
        if (g_oh.eq.0) then
          write(*,*) "moldiff: Error; no OH tracer !!!"
          stop
        endif
        g_ho2=igcm_ho2
        if (g_ho2.eq.0) then
          write(*,*) "moldiff: Error; no HO2 tracer !!!"
          stop
        endif
        g_h2o2=igcm_h2o2
        if (g_h2o2.eq.0) then
          write(*,*) "moldiff: Error; no H2O2 tracer !!!"
          stop
        endif
        g_n2=igcm_n2
        if (g_n2.eq.0) then
          write(*,*) "moldiff: Error; no N2 tracer !!!"
          stop
        endif
        g_ar=igcm_ar
        if (g_ar.eq.0) then
          write(*,*) "moldiff: Error; no AR tracer !!!"
          stop
        endif
        g_h2o=igcm_h2o_vap
        if (g_h2o.eq.0) then
          write(*,*) "moldiff: Error; no water vapor tracer !!!"
          stop
        endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   fill array to relate local indexes to gcm indexes
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        gcmind(i_co)	= g_co
        gcmind(i_n2)	= g_n2
        gcmind(i_o2)	= g_o2
        gcmind(i_co2)	= g_co2
        gcmind(i_h2)	= g_h2
        gcmind(i_h)	= g_h
        gcmind(i_oh)	= g_oh
        gcmind(i_ho2)	= g_ho2
        gcmind(i_h2o)	= g_h2o
        gcmind(i_h2o2)  = g_h2o2
        gcmind(i_o1d)   = g_o1d
        gcmind(i_o3)    = g_o3
        gcmind(i_o)     = g_o
        gcmind(i_ar)    = g_ar

        firstcall= .false.
      endif ! of if (firstcall)



c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      nz=nlayer

      do ig=1,ngrid

        do l=2,nz-1
          tt(l)=pt(ig,l)+pdt(ig,l)*ptimestep
     &                  +pdteuv(ig,l)*ptimestep
     &                  +pdtconduc(ig,l)*ptimestep

          do nn=1,ncompmoldiff
            qq(l,nn)=pq(ig,l,gcmind(nn))+pdq(ig,l,gcmind(nn))*ptimestep
            qq(l,nn)=max(qq(l,nn),1.e-30)
          enddo
          hp(l)=-log(pplay(ig,l+1)/pplay(ig,l-1))
          dmmeandz(l)=(mmean(ig,l+1)-mmean(ig,l-1))/hp(l)
        enddo

        tt(1)=pt(ig,1)  +pdt(ig,1)*ptimestep
     &                  +pdteuv(ig,1)*ptimestep
     &                  +pdtconduc(ig,1)*ptimestep
        tt(nz)=pt(ig,nz)+pdt(ig,nz)*ptimestep
     &                  +pdteuv(ig,nz)*ptimestep
     &                  +pdtconduc(ig,nz)*ptimestep

        do nn=1,ncompmoldiff
          qq(1,nn)=pq(ig,1,gcmind(nn))+pdq(ig,1,gcmind(nn))*ptimestep
          qq(nz,nn)=pq(ig,nz,gcmind(nn))+pdq(ig,nz,gcmind(nn))*ptimestep
          qq(1,nn)=max(qq(1,nn),1.e-30)
          qq(nz,nn)=max(qq(nz,nn),1.e-30)
        enddo
        hp(1)=-log(pplay(ig,2)/pplay(ig,1))
        dmmeandz(1)=(-3.*mmean(ig,1)+4.*mmean(ig,2)
     &               -mmean(ig,3))/(2.*hp(1))
        hp(nz)=-log(pplay(ig,nz)/pplay(ig,nz-1))
        dmmeandz(nz)=(3.*mmean(ig,nz)-4.*mmean(ig,nz-1)
     &                +mmean(ig,nz-2))/(2.*hp(nz))
c
c Setting-up matrix of alfa coefficients from Dickinson and Ridley 1972
c
      do l=1,nz
       if(abs(dmmeandz(l)) .lt. 1.e-5) dmmeandz(l)=0.0
        hh=rnew(ig,l)*tt(l)/g
        ptfac=(1.e5/pplay(ig,l))*(tt(l)/273)**1.75
        ntot=pplay(ig,l)/(1.381e-23*tt(l))      		! in #/m3

        do nn=1,ncompmoldiff-1            ! rows
          alfdiag(nn)=0.
          dcoef1=dij(nn,i_o)*ptfac
          do n=1,ncompmoldiff-1           ! columns
            y(nn,n)=0.
            dcoef=dij(nn,n)*ptfac
            alf(nn,n)=qq(l,nn)/ntot/1.66e-27
     &         *(1./(mmol(gcmind(n))*dcoef)-1./(mmol(g_o)*dcoef1))
            alfdiag(nn)=alfdiag(nn)
     &       +(1./(mmol(gcmind(n))*dcoef)-1./(mmol(g_o)*dcoef1))
     &        *qq(l,n)
          enddo
          dcoef=dij(nn,nn)*ptfac
          alfdiag(nn)=alfdiag(nn)
     &       -(1./(mmol(gcmind(nn))*dcoef)-1./(mmol(g_o)*dcoef1))
     &          *qq(l,nn)
          alf(nn,nn)=-(alfdiag(nn)
     &                 +1./(mmol(g_o)*dcoef1))/ntot/1.66e-27
          y(nn,nn)=1.
          b(l,nn)=-(dmmeandz(l)/mmean(ig,l)
     &              +mmol(gcmind(nn))/mmean(ig,l)-1.)
        enddo

c
c Inverting the alfa matrix
c
        call ludcmp_sp(alf,ncompmoldiff-1,ncompmoldiff-1,indx,d,ierr)

c       TEMPORAIRE *****************************
        if (ierr.ne.0) then
            write(*,*)'In moldiff: Problem in LUDCMP_SP with matrix alf'
            write(*,*) 'Singular matrix ?'
c           write(*,*) 'Matrix alf = ', alf
            write(*,*) 'ig, l=',ig, l
            write(*,*) 'No molecular diffusion this time !'
            pdqdiff(1:ngrid,1:nlayer,1:nq)=0
            return
c           stop
        end if
c       *******************************************
        do n=1,ncompmoldiff-1
       call lubksb_sp(alf,ncompmoldiff-1,ncompmoldiff-1,indx,y(1,n))
          do nn=1,ncompmoldiff-1
            alfinv(l,nn,n)=y(nn,n)/hh
          enddo
        enddo
      enddo

c
c Calculating coefficients of the system
c

c      zlocal(1)=-log(pplay(ig,1)/pplev(ig,1))* Rnew(ig,1)*tt(1)/g
      zlocal(1)=zzlay(ig,1)

      do l=2,nz-1
        del1=hp(l)*pplay(ig,l)/(g*ptimestep)
        del2=(hp(l)/2)**2*pplay(ig,l)/(g*ptimestep)
        do nn=1,ncompmoldiff-1
          do n=1,ncompmoldiff-1
            dalfinvdz=(alfinv(l+1,nn,n)-alfinv(l-1,nn,n))/hp(l)
            aa(l,nn,n)=-dalfinvdz/del1+alfinv(l,nn,n)/del2
     &                +alfinv(l-1,nn,n)*b(l-1,n)/del1    
            bb(l,nn,n)=-2.*alfinv(l,nn,n)/del2
            cc(l,nn,n)=dalfinvdz/del1+alfinv(l,nn,n)/del2
     &                -alfinv(l+1,nn,n)*b(l+1,n)/del1    
          enddo
        enddo

c        tmean=tt(l)
c        if(tt(l).ne.tt(l-1))
c     &        tmean=(tt(l)-tt(l-1))/log(tt(l)/tt(l-1))
c        zlocal(l)= zlocal(l-1)
c     &         -log(pplay(ig,l)/pplay(ig,l-1))*rnew(ig,l)*tmean/g
      zlocal(l)=zzlay(ig,l)
      enddo 

c      zlocal(nz)= zlocal(nz-1)
c     &         -log(pplay(ig,nz)/pplay(ig,nz-1))*rnew(ig,nz)*tmean/g
      zlocal(nz)=zzlay(ig,nz)
        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Escape velocity from Jeans equation for the flux of H and H2
c (Hunten 1973, eq. 5)
      
      do n=1,ncompmoldiff
        wi(n)=1.
        flux(n)=0.
        abfac(n)=1.
      enddo

       dens=pplay(ig,nz)/(rnew(ig,nz)*tt(nz))
c
c For H:
c
       pote=(3398000.+zlocal(nz))/
     &         (1.381e-23*tt(nz)/(1.6605e-27*mmol(g_h)*g))
       wi(i_h)=sqrt(2.*1.381e-23*tt(nz)/(1.6605e-27*mmol(g_h)))
     &             /(2.*sqrt(3.1415))*(1.+pote)*exp(-pote)
       flux(i_h)=qq(nz,i_h)*dens/(1.6605e-27*mmol(g_h))*wi(i_h)
       flux(i_h)=flux(i_h)*1.6606e-27
       abfac(i_h)=0.
c
c For H2:
c
       pote=(3398000.+zlocal(nz))/
     &         (1.381e-23*tt(nz)/(1.6605e-27*mmol(g_h2)*g))
       wi(i_h2)=sqrt(2.*1.381e-23*tt(nz)/(1.6605e-27*mmol(g_h2)))
     &              /(2.*sqrt(3.1415))*(1.+pote)*exp(-pote)
       flux(i_h2)=qq(nz,i_h2)*dens/(1.6605e-27*mmol(g_h2))*wi(i_h2)
       flux(i_h2)=flux(i_h2)*1.6606e-27
       abfac(i_h2)=0.

c ********* TEMPORAIRE : no escape for h and h2
c     do n=1,ncomptot
c       wi(n)=1.
c       flux(n)=0.
c       abfac(n)=1.
c     enddo
c ********************************************


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c Setting coefficients for tridiagonal matrix and solving the system
c

      do nn=1,ncompmoldiff-1
        do l=2,nz-1
          atri(l-1)=aa(l,nn,nn)
          btri(l-1)=bb(l,nn,nn)+1.
          ctri(l-1)=cc(l,nn,nn)
          rtri(l-1)=qq(l,nn)
          do n=1,ncompmoldiff-1
            rtri(l-1)=rtri(l-1)-(aa(l,nn,n)*qq(l-1,n)
     &                           +bb(l,nn,n)*qq(l,n)
     &                           +cc(l,nn,n)*qq(l+1,n))
          enddo
          rtri(l-1)=rtri(l-1)+(aa(l,nn,nn)*qq(l-1,nn)
     &                         +bb(l,nn,nn)*qq(l,nn)
     &                         +cc(l,nn,nn)*qq(l+1,nn))
        enddo

c
c Boundary conditions:
c	Escape flux for H and H2 at top
c	Diffusive equilibrium for the other species at top
c	Perfect mixing for all at bottom
c

        rtri(nz-2)=rtri(nz-2)
     &             -ctri(nz-2)*flux(nn)*mmol(gcmind(nn))/(dens*wi(nn))

        atri(nz-2)=atri(nz-2)
     &             -abfac(nn)*ctri(nz-2)/(3.-2.*hp(nz)*b(nz,nn))
        btri(nz-2)=btri(nz-2)
     &             +abfac(nn)*4.*ctri(nz-2)/(3.-2.*hp(nz)*b(nz,nn))

c        rtri(1)=rtri(1)-atri(1)*qq(1,nn)
        btri(1)=btri(1)+atri(1)

        call tridag_sp(atri,btri,ctri,rtri,qtri,nz-2)

        do l=2,nz-1
c          qnew(l,nn)=qtri(l-1)
          qnew(l,nn)=max(qtri(l-1),1.e-30)
        enddo

        qnew(nz,nn)=flux(nn)*mmol(gcmind(nn))/(dens*wi(nn))
     &               +abfac(nn)*(4.*qnew(nz-1,nn)-qnew(nz-2,nn))
     &                /(3.-2.*hp(nz)*b(nz,nn))
c        qnew(1,nn)=qq(1,nn)
        qnew(1,nn)=qnew(2,nn)
         
        qnew(nz,nn)=max(qnew(nz,nn),1.e-30)
        qnew(1,nn)=max(qnew(1,nn),1.e-30)

      enddo 	! loop on species

      DO l=1,nz
        if(zlocal(l).gt.65000.)then
        pdqdiff(ig,l,g_o)=0.
        do n=1,ncompmoldiff-1
          pdqdiff(ig,l,gcmind(n))=(qnew(l,n)-qq(l,n))
          pdqdiff(ig,l,g_o)=pdqdiff(ig,l,g_o)-(qnew(l,n)-qq(l,n))
          pdqdiff(ig,l,gcmind(n))=pdqdiff(ig,l,gcmind(n))/ptimestep
        enddo
		pdqdiff(ig,l,g_o)=pdqdiff(ig,l,g_o)/ptimestep
        endif
      ENDDO

c      do l=2,nz
c        do n=1,ncomptot-1
c          hco2(n)=qnew(l,n)*pplay(ig,l)/(rnew(ig,l)*tt(l)) /
c     &      (qnew(l-1,n)*pplay(ig,l-1)/(rnew(ig,l-1)*tt(l-1))) 
c          hco2(n)=-(zlocal(l)-zlocal(l-1))/log(hco2(n))/1000.
c        enddo
c        write(225,*),l,pt(1,l),(hco2(n),n=1,6)
c        write(226,*),l,pt(1,l),(hco2(n),n=7,12)
c      enddo

      enddo		! ig loop

      return
      end

c    ********************************************************************
c    ********************************************************************
c    ********************************************************************
 
      subroutine tridag_sp(a,b,c,r,u,n)
c      parameter (nmax=100)
c      dimension gam(nmax),a(n),b(n),c(n),r(n),u(n)
      real gam(n),a(n),b(n),c(n),r(n),u(n)
      if(b(1).eq.0.)then
        stop 'tridag_sp: error: b(1)=0 !!! '
      endif
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.) then
          stop 'tridag_sp: error: bet=0 !!! '
        endif
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      return
      end

c    ********************************************************************
c    ********************************************************************
c    ********************************************************************

      SUBROUTINE LUBKSB_SP(A,N,NP,INDX,B)

      implicit none

      integer i,j,n,np,ii,ll
      real sum
      real a(np,np),indx(np),b(np)

c      DIMENSION A(NP,NP),INDX(N),B(N)
      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14    CONTINUE
      RETURN
      END

c    ********************************************************************
c    ********************************************************************
c    ********************************************************************

      SUBROUTINE LUDCMP_SP(A,N,NP,INDX,D,ierr)

      implicit none

      integer n,np,nmax,i,j,k,imax
      real d,tiny,aamax
      real a(np,np),indx(np)
      integer ierr  ! error =0 if OK, =1 if problem

      PARAMETER (NMAX=100,TINY=1.0E-20)
c      DIMENSION A(NP,NP),INDX(N),VV(NMAX)
      real sum,vv(nmax),dum

      D=1.
      DO 12 I=1,N
        AAMAX=0.
        DO 11 J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
11      CONTINUE
        IF (AAMAX.EQ.0.) then
            write(*,*) 'In moldiff: Problem in LUDCMP_SP with matrix A'
            write(*,*) 'Singular matrix ?'
c           write(*,*) 'Matrix A = ', A
c           TO DEBUG :
            ierr =1
            return
c           stop
        END IF 

        VV(I)=1./AAMAX
12    CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUM=1./A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(A(N,N).EQ.0.)A(N,N)=TINY
      ierr =0
      RETURN
      END

