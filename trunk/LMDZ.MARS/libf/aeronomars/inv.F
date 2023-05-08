      SUBROUTINE inv(A,N)

      implicit none

c=======================================================================
c   Scheme A
c
c   subject:
c   --------
c
c   Inversion of a matrix:
c   Combinaison des routines ludcmp.f et lubksb.f du Numerical Recipes
c=======================================================================

      real       TINY
      PARAMETER (TINY=1.0E-20)
      real       AAMAX,DUM,SUM
      integer    IMAX,I,II,J,K,L,LL,N
      real       A(N,N),INDX(N),VV(N,N)

      IMAX = 0
      do I=1,N
          AAMAX=0.
          do J=1,N
            IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
          enddo
          if (AAMAX.EQ.0.) then
              write(*,*) 'Singular matrix.'
              stop
          endif
          VV(I,1)=1./AAMAX
      enddo
      do J=1,N
        IF (J.GT.1) THEN
          do I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              do K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
              enddo
              A(I,J)=SUM
            ENDIF
          enddo
        ENDIF
        AAMAX=0.
        do I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            do K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
            enddo
            A(I,J)=SUM
          ENDIF
          DUM=VV(I,1)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
        enddo
        IF (J.NE.IMAX)THEN
          do K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
          enddo
          VV(IMAX,1)=VV(J,1)
        ENDIF
        INDX(J)=IMAX
        if(abs(A(J,J)).LT.TINY) then
              write(*,*) 'Pivot too small.'
              stop
        endif
        IF(J.NE.N)THEN
          DUM=1./A(J,J)
          do I=J+1,N
            A(I,J)=A(I,J)*DUM
          enddo
        ENDIF
      enddo

      do I=1,N
         do J=1,N
            VV(I,J) = 0.
         enddo
         VV(I,I) = 1.
      enddo

      do L=1,N
        II=0
        do I=1,N
          LL=INDX(I)
          SUM=VV(LL,L)
          VV(LL,L)=VV(I,L)
          IF (II.NE.0)THEN
            do J=II,I-1
              SUM=SUM-A(I,J)*VV(J,L)
            enddo
          ELSE IF (SUM.NE.0.) THEN
            II=I
          ENDIF
          VV(I,L)=SUM
        enddo
        do I=N,1,-1
          SUM=VV(I,L)
          IF(I.LT.N)THEN
            do J=I+1,N
              SUM=SUM-A(I,J)*VV(J,L)
            enddo
          ENDIF
          VV(I,L)=SUM/A(I,I)
        enddo
      enddo

      do I=1,N
        do L=1,N
          A(I,L)=VV(I,L)
        enddo
      enddo

      return
      END
