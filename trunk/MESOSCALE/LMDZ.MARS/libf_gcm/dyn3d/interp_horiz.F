      subroutine interp_horiz (varo,varn,imo,jmo,imn,jmn,lm,
     &  rlonuo,rlatvo,rlonun,rlatvn)  

c===========================================================
c  Interpolation Horizontales des variables d'une grille LMDZ
c (des points SCALAIRES au point SCALAIRES)
c  dans une autre grille LMDZ en conservant la quantite
c  totale pour les variables intensives (/m2) : ex : Pression au sol
c
c Francois Forget (01/1995)
c===========================================================

      IMPLICIT NONE 

c   Declarations:
c ==============
c
c  ARGUMENTS
c  """""""""
        
       INTEGER imo, jmo ! dimensions ancienne grille (input)
       INTEGER imn,jmn  ! dimensions nouvelle grille (input)

       REAL rlonuo(imo+1)     !  Latitude et
       REAL rlatvo(jmo)       !  longitude des
       REAL rlonun(imn+1)     !  bord des 
       REAL rlatvn(jmn)     !  cases "scalaires" (input)

       INTEGER lm ! dimension verticale (input)
       REAL varo (imo+1, jmo+1,lm) ! var dans l'ancienne grille (input)
       REAL varn (imn+1,jmn+1,lm) ! var dans la nouvelle grille (output)

c Autres variables
c """"""""""""""""
       INTEGER imnmx2,jmnmx2
       parameter (imnmx2=190,jmnmx2=100)
       REAL airetest(imnmx2+1,jmnmx2+1)
       INTEGER ii,jj,l

       REAL airen ((imnmx2+1)*(jmnmx2+1)) ! aire dans la nouvelle grille
       REAL airentotn	! aire totale pole nord dans la nouvelle grille
       REAL airentots	! aire totale pole sud dans la nouvelle grille
c    Info sur les ktotal intersection entre les cases new/old grille

c kmax: le nombre  max  d'intersections entre les 2 grilles horizontales
c On fixe kmax a la taille de la grille des donnees martiennes (360x179) 
c + des pouiemes (cas ou une maille est a cheval sur 2 ou 4 mailles)
c  Il y a un test dans iniinterp_h pour s'assurer que ktotal < kmax
       INTEGER kmax, k, ktotal
       parameter (kmax = 360*179 + 200000)
c      parameter (kmax = 360*179 + 40000)

       INTEGER iik(kmax), jjk(kmax),jk(kmax),ik(kmax)
       REAL intersec(kmax)
       REAL R
       REAL totn, tots
       REAL prev_sumdim
       save prev_sumdim
       

       logical firsttest, aire_ok
       save firsttest
       data firsttest /.true./
       data aire_ok /.true./

       integer imoS,jmoS,imnS,jmnS
       save imoS,jmoS,imnS,jmnS
       save ktotal,iik,jjk,jk,ik,intersec,airen
       REAL pi

c Test dimensions imnmx2 jmnmx2
c""""""""""""""""""""""""""""""
c test dimensionnement tableau airetest
      if (imn.GT.imnmx2.OR.jmn.GT.jmnmx2) then
         write(*,*) 'STOP pb dimensionnement tableau airetest'
         write(*,*) 'il faut imn < imnmx2 et jmn < jmnmx2'
         write(*,*) 'imn imnmx2', imn,imnmx2
         write(*,*) 'jmn jmnmx2', jmn,jmnmx2
         call exit(1)
      endif

c initialisation
c --------------
c Si c'est le premier appel,  on prepare l'interpolation
c en calculant pour chaque case autour d'un point scalaire de la
c nouvelle grille, la surface  de intersection avec chaque
c    case de l'ancienne grille.

c  This must also be done if we change the dimension
      if (imo+jmo+imn+jmn.ne.prev_sumdim) then
          firsttest=.true.
          prev_sumdim=imo+jmo+imn+jmn
      end if       

      if (firsttest) then 
        call iniinterp_h(imo,jmo,imn,jmn ,kmax,
     &       rlonuo,rlatvo,rlonun,rlatvn,
     &          ktotal,iik,jjk,jk,ik,intersec,airen)
       imoS=imo
       jmoS=jmo
       imnS=imn
       jmnS=jmn
      else
       if(imo.NE.imoS.OR.jmo.NE.jmoS.OR.imn.NE.imnS.OR.jmn.NE.jmnS) then
        call iniinterp_h(imo,jmo,imn,jmn ,kmax,
     &       rlonuo,rlatvo,rlonun,rlatvn,
     &          ktotal,iik,jjk,jk,ik,intersec,airen)
       imoS=imo
       jmoS=jmo
       imnS=imn
       jmnS=jmn
       end if
      end if

      do l=1,lm
       do jj =1 , jmn+1
        do ii=1, imn+1
          varn(ii,jj,l) =0.
        end do
       end do
      end do 
       
c Interpolation horizontale
c -------------------------
c boucle sur toute les ktotal intersections entre les cases
c de l'ancienne et la  nouvelle grille
c
     
      do k=1,ktotal
        do l=1,lm
         varn(iik(k),jjk(k),l) = varn(iik(k),jjk(k),l) 
     &   + varo(ik(k), jk(k),l)*intersec(k)/airen(iik(k)
     &   +(jjk(k)-1)*(imn+1))
        end do
      end do

c Une seule valeur au pole pour les variables ! :
c -----------------------------------------------
      DO l=1, lm
         totn =0.
         tots =0.


c moyenne du champ au poles (ponderee par les aires)
c"""""""""""""""""""""""""""""""
         airentotn=0.
         airentots=0.

         do ii =1, imn+1
            totn = totn + varn(ii,1,l)*airen(ii)
            tots = tots + varn (ii,jmn+1,l)*airen(jmn*(imn+1)+ii)
            airentotn=airentotn + airen(ii)
            airentots=airentots + airen(jmn*(imn+1)+ii)
         end do 

         do ii =1, imn+1
            varn(ii,1,l) = totn/airentotn
            varn(ii,jmn+1,l) = tots/airentots
         end do 

      ENDDO
           

c---------------------------------------------------------------
c  TEST  TEST  TEST  TEST  TEST  TEST  TEST  TEST  TEST  TEST 
      if (firsttest) then
      pi=2.*asin(1.)
      firsttest = .false.
      write (*,*) 'INTERP. HORIZ. : TEST SUR LES AIRES:'

      do jj =1 , jmn+1
        do ii=1, imn+1
          airetest(ii,jj) =0.
        end do
      end do 
      do k=1,ktotal
         airetest(iik(k),jjk(k))= airetest(iik(k),jjk(k)) +intersec(k) 
      end do
      do jj =1 , jmn+1
       do ii=1, imn+1
         r = airen(ii+(jj-1)*(imn+1))/airetest(ii,jj)
         if ((r.gt.1.001).or.(r.lt.0.999)) then
             write (*,*) '********** PROBLEME D'' AIRES !!!',
     &                   ' DANS L''INTERPOLATION HORIZONTALE'
             write(*,*)'ii,jj,airen,airetest',
     &          ii,jj,airen(ii+(jj-1)*(imn+1)),airetest(ii,jj)
             aire_ok = .false.
         end if
       end do
      end do
      if (aire_ok) write(*,*) 'INTERP. HORIZ. : AIRES OK'
      endif

c FIN TEST  FIN TEST  FIN TEST  FIN TEST  FIN TEST  FIN TEST  FIN TEST
c --------------------------------------------------------------------


        return
        end
