          SUBROUTINE eigen_sort(d,v,n,np)
          INTEGER n,np
          REAL d(np),v(np,np)
          INTEGER i,j,k
          REAL p

       DO i=1,n-1
          k=i
          p=d(i)
        DO j=i+1,n
           IF(d(j).ge.p) THEN
            k=j
            p=d(j)
           ENDIF
        ENDDO
          
        IF(k.ne.i) THEN
          d(k)=d(i)
          d(i)=p
         DO j=1,n
          p=v(j,i)
          v(j,i)=v(j,k)
          v(j,k)=p
         ENDDO
        ENDIF
       ENDDO

        RETURN
        END
