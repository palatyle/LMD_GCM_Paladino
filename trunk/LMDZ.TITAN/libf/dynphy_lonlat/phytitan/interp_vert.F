c******************************************************
      SUBROUTINE   interp_vert(varo,varn,lmo,lmn,apso,bpso,
     &             aps,bps,ps,Nhoriz)
c
c interpolation lineaire pour passer
c a une nouvelle discretisation verticale pour
c les variables de GCM
c Francois Forget (01/1995)
c Modif pour coordonnees hybrides FF (03/2003)
c**********************************************************

      IMPLICIT NONE

c   Declarations:
c ==============
c
c  ARGUMENTS
c  """""""""

       integer lmo ! dimensions ancienne couches (input)
       integer lmn ! dimensions nouvelle couches (input)

       real apso(lmo),bpso(lmo)! anciennes coord hybrides midlayer (input)
       real aps(lmn), bps(lmn)! nouvelles coord hybrides (midlayer) (input)

       integer Nhoriz ! nombre de point horizontale (input)
       real ps(nhoriz) !pression de surface (input)

       real varo(Nhoriz,lmo) ! var dans l''ancienne grille (input)
       real varn(Nhoriz,lmn) ! var dans la nouvelle grille (output)

c Autres variables
c """"""""""""""""
       integer n, ln ,lo 
       real coef
       REAL sigmo(lmo) ! niveau sigma des variables dans les anciennes coord
       REAL sigmn(lmn) ! niveau sigma des variables dans les nouvelles coord

c run
c ====

      do n=1,Nhoriz
        do ln=1,lmn
            sigmn(ln)=aps(ln)/ps(n)+bps(ln)
        end do
        do lo=1,lmo
            sigmo(lo)=apso(lo)/ps(n)+bpso(lo)
        end do

        do ln=1,lmn
           if (sigmn(ln).ge.sigmo(1))then
             varn(n,ln) =  varo(n,1)  
           else if (sigmn(ln).le.sigmo(lmo)) then
             varn(n,ln) =  varo(n,lmo)
           else
              do lo =1,lmo-1 
                if ( (sigmn(ln).le.sigmo(lo)).and.
     &             (sigmn(ln).gt.sigmo(lo+1)) )then
                  coef = (sigmn(ln)-sigmo(lo))/(sigmo(lo+1)-sigmo(lo))
                   varn(n,ln)=varo(n,lo) +coef*(varo(n,lo+1)-varo(n,lo))
                end if
              end do           
           end if
         end do

      end do


      return
      end
