       subroutine scatter(n,a,index,b)
C
       implicit none
       integer N,INDEX(n),I
       real A(n),B(n)
c
c
       DO 100 I=1,N
          A(INDEX(I))=B(I)
100    CONTINUE
c
       return
       end 
c
c
c      peut etre vec_scatter(v,rindices,count,r) en f77hp
c      equivalent code:
c      do 99 i=1,count
c99    r(rindices(i))=v(i)
c      voir p. 11-13
c
