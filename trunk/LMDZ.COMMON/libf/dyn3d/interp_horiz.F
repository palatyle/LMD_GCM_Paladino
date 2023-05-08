c
c $Id: interp_horiz.F 1403 2010-07-01 09:02:53Z fairhead $
c
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
        
       integer imo, jmo ! dimensions ancienne grille (input)
       integer imn,jmn  ! dimensions nouvelle grille (input)

       real rlonuo(imo+1)     !  Latitude et
       real rlatvo(jmo)       !  longitude des
       real rlonun(imn+1)     !  bord des 
       real rlatvn(jmn)     !  cases "scalaires" (input)

       integer lm ! dimension verticale (input)
       real varo (imo+1, jmo+1,lm) ! var dans l'ancienne grille (input)
       real varn (imn+1,jmn+1,lm) ! var dans la nouvelle grille (output)

c Autres variables
c """"""""""""""""
       real airetest(imn+1,jmn+1)
       integer ii,jj,l

       real airen (imn+1,jmn+1) ! aire dans la nouvelle grille
c    Info sur les ktotal intersection entre les cases new/old grille
       integer kllm, k, ktotal
       parameter (kllm = 400*200*10)
       integer iik(kllm), jjk(kllm),jk(kllm),ik(kllm)
       real intersec(kllm)
       real R
       real totn, tots

       logical firstcall, firsttest, aire_ok
       save firsttest
       data firsttest /.true./
       data aire_ok /.true./

       



c initialisation
c --------------
c Si c'est le premier appel, on prepare l'interpolation
c en calculant pour chaque case autour d'un point scalaire de la
c nouvelle grille, la surface  de intersection avec chaque
c    case de l'ancienne grille.


        call iniinterp_horiz (imo,jmo,imn,jmn ,kllm,
     &       rlonuo,rlatvo,rlonun,rlatvn,
     &          ktotal,iik,jjk,jk,ik,intersec,airen)

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
      PRINT *, 'ktotal 1 = ', ktotal
     
      do k=1,ktotal
        do l=1,lm
         varn(iik(k),jjk(k),l) = varn(iik(k),jjk(k),l) 
     &        + varo(ik(k), jk(k),l)*intersec(k)/airen(iik(k),jjk(k))
        end do
      end do

c Une seule valeur au pole pour les variables ! :
c -----------------------------------------------
       do l=1, lm
         totn =0.
         tots =0.
           do ii =1, imn+1
             totn = totn + varn(ii,1,l)
             tots = tots + varn (ii,jmn+1,l)
           end do 
           do ii =1, imn+1
             varn(ii,1,l) = totn/REAL(imn+1)
             varn(ii,jmn+1,l) = tots/REAL(imn+1)
           end do 
       end do
           

c---------------------------------------------------------------
c  TEST  TEST  TEST  TEST  TEST  TEST  TEST  TEST  TEST  TEST 
!!       if (.not.(firsttest)) goto 99
!!       firsttest = .false.
!! !     write (*,*) 'INTERP. HORIZ. : TEST SUR LES AIRES:'
!!       do jj =1 , jmn+1
!!         do ii=1, imn+1
!!           airetest(ii,jj) =0.
!!         end do
!!       end do 
!!       PRINT *, 'ktotal = ', ktotal
!!       PRINT *, 'jmn+1 =', jmn+1, 'imn+1', imn+1
!! 
!!       do k=1,ktotal
!!          airetest(iik(k),jjk(k))= airetest(iik(k),jjk(k)) +intersec(k) 
!!       end DO
!! 
!! 
!!       PRINT *, 'fin boucle'
!!       do jj =1 , jmn+1
!!        do ii=1, imn+1
!!          r = airen(ii,jj)/airetest(ii,jj)
!!          if ((r.gt.1.001).or.(r.lt.0.999)) then
!! !             write (*,*) '********** PROBLEME D'' AIRES !!!',
!! !     &                   ' DANS L''INTERPOLATION HORIZONTALE'
!! !             write(*,*)'ii,jj,airen,airetest',
!! !     &          ii,jj,airen(ii,jj),airetest(ii,jj)
!!              aire_ok = .false.
!!          end if
!!        end do
!!       end do
!! !      if (aire_ok) write(*,*) 'INTERP. HORIZ. : AIRES OK'
!!  99   continue

c FIN TEST  FIN TEST  FIN TEST  FIN TEST  FIN TEST  FIN TEST  FIN TEST
c---------------------------------------------------------------








        return
        end
