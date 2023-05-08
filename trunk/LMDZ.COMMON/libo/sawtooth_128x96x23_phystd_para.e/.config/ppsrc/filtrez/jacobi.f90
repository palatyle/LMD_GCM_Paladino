










!
! $Id: jacobi.F90 1289 2009-12-18 14:51:22Z emillour $
!
      SUBROUTINE JACOBI(A,N,NP,D,V,NROT)
      implicit none
! Arguments:
      integer,intent(in) :: N
      integer,intent(in) :: NP
      integer,intent(out) :: NROT
      real,intent(inout) :: A(NP,NP)
      real,intent(out) :: D(NP)
      real,intent(out) :: V(NP,NP)

! local variables:
      integer :: IP,IQ,I,J
      real :: SM,TRESH,G,H,T,THETA,C,S,TAU
      real :: B(N)
      real :: Z(N)
      
      DO IP=1,N
        DO IQ=1,N
          V(IP,IQ)=0.
        ENDDO
        V(IP,IP)=1.
      ENDDO
      DO IP=1,N
        B(IP)=A(IP,IP)
        D(IP)=B(IP)
        Z(IP)=0.
      ENDDO
      NROT=0
      DO I=1,50 ! 50? I suspect this should be NP
                !     but convergence is fast enough anyway
        SM=0.
        DO IP=1,N-1
          DO IQ=IP+1,N
            SM=SM+ABS(A(IP,IQ))
          ENDDO
        ENDDO
        IF(SM.EQ.0.)RETURN
        IF(I.LT.4)THEN
          TRESH=0.2*SM/N**2
        ELSE
          TRESH=0.
        ENDIF
        DO IP=1,N-1
          DO IQ=IP+1,N
            G=100.*ABS(A(IP,IQ))
            IF((I.GT.4).AND.(ABS(D(IP))+G.EQ.ABS(D(IP))) &
               .AND.(ABS(D(IQ))+G.EQ.ABS(D(IQ))))THEN
              A(IP,IQ)=0.
            ELSE IF(ABS(A(IP,IQ)).GT.TRESH)THEN
              H=D(IQ)-D(IP)
              IF(ABS(H)+G.EQ.ABS(H))THEN
                T=A(IP,IQ)/H
              ELSE
                THETA=0.5*H/A(IP,IQ)
                T=1./(ABS(THETA)+SQRT(1.+THETA**2))
                IF(THETA.LT.0.)T=-T
              ENDIF
              C=1./SQRT(1+T**2)
              S=T*C
              TAU=S/(1.+C)
              H=T*A(IP,IQ)
              Z(IP)=Z(IP)-H
              Z(IQ)=Z(IQ)+H
              D(IP)=D(IP)-H
              D(IQ)=D(IQ)+H
              A(IP,IQ)=0.
              DO J=1,IP-1
                G=A(J,IP)
                H=A(J,IQ)
                A(J,IP)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
             ENDDO
              DO J=IP+1,IQ-1
                G=A(IP,J)
                H=A(J,IQ)
                A(IP,J)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
              ENDDO
              DO J=IQ+1,N
                G=A(IP,J)
                H=A(IQ,J)
                A(IP,J)=G-S*(H+G*TAU)
                A(IQ,J)=H+S*(G-H*TAU)
              ENDDO
              DO J=1,N
                G=V(J,IP)
                H=V(J,IQ)
                V(J,IP)=G-S*(H+G*TAU)
                V(J,IQ)=H+S*(G-H*TAU)
              ENDDO
              NROT=NROT+1
            ENDIF
          ENDDO
        ENDDO
        DO IP=1,N
          B(IP)=B(IP)+Z(IP)
          D(IP)=B(IP)
          Z(IP)=0.
        ENDDO
      ENDDO ! of DO I=1,50
      STOP 'Jacobi: 50 iterations should never happen'
      RETURN
      END SUBROUTINE JACOBI
