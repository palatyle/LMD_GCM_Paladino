C 
C $Header$
C
      subroutine iniinterp_horiz (imo,jmo,imn,jmn ,kllm,
     &       rlonuo,rlatvo,rlonun,rlatvn,
     &       ktotal,iik,jjk,jk,ik,intersec,airen)
   
      implicit none



c ---------------------------------------------------------
c Prepare l' interpolation des variables d'une grille LMDZ
c  dans une autre grille LMDZ en conservant la quantite
c  totale pour les variables intensives (/m2) : ex : Pression au sol
c
c   (Pour chaque case autour d'un point scalaire de la nouvelle
c    grille, on calcule la surface (en m2)en intersection avec chaque
c    case de l'ancienne grille , pour la future interpolation)
c
c on calcule aussi l' aire dans la nouvelle grille 
c
c
c   Auteur:  F.Forget 01/1995
c   -------
c
c ---------------------------------------------------------
c   Declarations:
c ==============
c
c  ARGUMENTS
c  """""""""
c INPUT
       integer imo, jmo ! dimensions ancienne grille
       integer imn,jmn  ! dimensions nouvelle grille
       integer kllm ! taille du tableau des intersections
       real rlonuo(imo+1)     !  Latitude et
       real rlatvo(jmo)       !  longitude des
       real rlonun(imn+1)     !  bord des
       real rlatvn(jmn)     !  cases "scalaires" (input)

c OUTPUT
       integer ktotal ! nombre totale d'intersections reperees
       integer iik(kllm), jjk(kllm),jk(kllm),ik(kllm)
       real intersec(kllm)  ! surface des intersections (m2)
       real airen (imn+1,jmn+1) ! aire dans la nouvelle grille


       
 
c Autres variables
c """"""""""""""""
       integer i,j,ii,jj,k
       real a(imo+1),b(imo+1),c(jmo+1),d(jmo+1)
       real an(imn+1),bn(imn+1),cn(jmn+1),dn(jmn+1)
       real aa, bb,cc,dd
       real pi

       pi      = 2.*ASIN( 1. )



c On repere les frontieres des cases :
c =================================== 
c
c Attention, on ruse avec des latitudes = 90 deg au pole.


c  ANcienne grile
c  """"""""""""""
      a(1) =   - rlonuo(imo+1)
      b(1) = rlonuo(1)
      do i=2,imo+1
         a(i) = rlonuo(i-1)
         b(i) =  rlonuo(i)
      end do

      d(1) = pi/2 
      do j=1,jmo
         c(j) = rlatvo(j) 
         d(j+1) = rlatvo(j)
      end do
      c(jmo+1) = -pi/2 
      

c  Nouvelle grille
c  """""""""""""""
      an(1) =  - rlonun(imn+1)
      bn(1) = rlonun(1)
      do i=2,imn+1
         an(i) = rlonun(i-1)
         bn(i) =  rlonun(i)
      end do

      dn(1) = pi/2 
      do j=1,jmn
         cn(j) = rlatvn(j)
         dn(j+1) = rlatvn(j)
      end do
      cn(jmn+1) = -pi/2 

c Calcul de la surface des cases scalaires de la nouvelle grille
c ==============================================================
      do ii=1,imn + 1
        do jj = 1,jmn+1
               airen(ii,jj) = (bn(ii)-an(ii))*(sin(dn(jj))-sin(cn(jj)))
        end do
      end do

c Calcul de la surface des intersections
c ======================================

c     boucle sur la nouvelle grille
c     """"""""""""""""""""""""""""
      ktotal = 0
      do jj = 1,jmn+1
       do j=1, jmo+1
          if((cn(jj).lt.d(j)).and.(dn(jj).gt.c(j)))then
              do ii=1,imn + 1
                do i=1, imo +1
                    if (  ((an(ii).lt.b(i)).and.(bn(ii).gt.a(i)))
     &        .or. ((an(ii).lt.b(i)-2*pi).and.(bn(ii).gt.a(i)-2*pi)
     &             .and.(b(i)-2*pi.lt.-pi) )
     &        .or. ((an(ii).lt.b(i)+2*pi).and.(bn(ii).gt.a(i)+2*pi)
     &             .and.(a(i)+2*pi.gt.pi) )
     &                     )then
                      ktotal = ktotal +1
                      iik(ktotal) =ii
                      jjk(ktotal) =jj
                      ik(ktotal) =i
                      jk(ktotal) =j

                      dd = min(d(j), dn(jj))
                      cc = cn(jj)
                      if (cc.lt. c(j))cc=c(j)
                      if((an(ii).lt.b(i)-2*pi).and.
     &                  (bn(ii).gt.a(i)-2*pi)) then 
                          bb = min(b(i)-2*pi,bn(ii))
                          aa = an(ii)
                          if (aa.lt.a(i)-2*pi) aa=a(i)-2*pi
                      else if((an(ii).lt.b(i)+2*pi).and.
     &                       (bn(ii).gt.a(i)+2*pi)) then
                          bb = min(b(i)+2*pi,bn(ii))
                          aa = an(ii)
                          if (aa.lt.a(i)+2*pi) aa=a(i)+2*pi
                      else 
                          bb = min(b(i),bn(ii))
                          aa = an(ii)
                          if (aa.lt.a(i)) aa=a(i)
                      end if
                      intersec(ktotal)=(bb-aa)*(sin(dd)-sin(cc))
                     end if
                end do
               end do
             end if
         end do
       end do       



c     TEST  INFO
c     DO k=1,ktotal 
c      ii = iik(k) 
c      jj = jjk(k)
c      i = ik(k)
c      j = jk(k)
c      if ((ii.eq.10).and.(jj.eq.10).and.(i.eq.10).and.(j.eq.10))then
c      if (jj.eq.2.and.(ii.eq.1))then
c          write(*,*) '**************** jj=',jj,'ii=',ii
c          write(*,*) 'i,j =',i,j
c          write(*,*) 'an,bn,cn,dn', an(ii), bn(ii), cn(jj),dn(jj)
c          write(*,*) 'a,b,c,d', a(i), b(i), c(j),d(j)
c          write(*,*) 'intersec(k)',intersec(k)
c          write(*,*) 'airen(ii,jj)=',airen(ii,jj)
c      end if
c     END DO 

      return
      end
