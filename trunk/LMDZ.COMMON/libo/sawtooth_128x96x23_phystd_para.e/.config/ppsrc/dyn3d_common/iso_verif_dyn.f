










        function iso_verif_noNaN_nostop(x,err_msg)
        implicit none
        ! si x est NaN, on affiche message
        ! d'erreur et return 1 si erreur

        ! input:
        real x
        character*(*) err_msg ! message d''erreur à afficher

        ! output 
        real borne
        parameter (borne=1e19)
        integer iso_verif_noNaN_nostop

        if ((x.gt.-borne).and.(x.lt.borne)) then
                iso_verif_noNAN_nostop=0
        else
            write(*,*) 'erreur detectee par iso_verif_nonNaN:'
            write(*,*) err_msg
            write(*,*) 'x=',x
            iso_verif_noNaN_nostop=1
        endif      

        return
        end

        function iso_verif_egalite_nostop
     :           (a,b,err_msg)
        implicit none
        ! compare a et b. Si pas egal, on affiche message
        ! d'erreur et stoppe
        ! pour egalite, on verifie erreur absolue et arreur relative

        ! input:
        real a, b
        character*(*) err_msg ! message d''erreur à afficher

        ! locals
        real errmax ! erreur maximale en absolu.
        real errmaxrel ! erreur maximale en relatif autorisée
        parameter (errmax=1e-8)
        parameter (errmaxrel=1e-3)

        ! output
        integer iso_verif_egalite_nostop

        iso_verif_egalite_nostop=0

        if (abs(a-b).gt.errmax) then
          if (abs((a-b)/max(max(abs(b),abs(a)),1e-18))
     :            .gt.errmaxrel) then
            write(*,*) 'erreur detectee par iso_verif_egalite:'
            write(*,*) err_msg
            write(*,*) 'a=',a
            write(*,*) 'b=',b
            iso_verif_egalite_nostop=1
          endif
        endif     
        
        return
        end        


        function iso_verif_aberrant_nostop
     :           (x,iso,q,err_msg)
        USE infotrac
        implicit none
        
        ! input:
        real x,q
        integer iso ! 2=HDO, 1=O18
        character*(*) err_msg ! message d''erreur à afficher

        ! locals
        real qmin,deltaD
        real deltaDmax,deltaDmin
        parameter (qmin=1e-11)
        parameter (deltaDmax=200.0,deltaDmin=-999.9)

        ! output
        integer iso_verif_aberrant_nostop

        iso_verif_aberrant_nostop=0

        ! verifier que HDO est raisonable
         if (q.gt.qmin) then
             deltaD=(x/q/tnat(iso)-1)*1000
             if ((deltaD.gt.deltaDmax).or.(deltaD.lt.deltaDmin)) then
                  write(*,*) 'erreur detectee par iso_verif_aberrant:'
                  write(*,*) err_msg
                  write(*,*) 'q=',q
                  write(*,*) 'deltaD=',deltaD
                  write(*,*) 'iso=',iso
                  iso_verif_aberrant_nostop=1
             endif !if ((deltaD.gt.deltaDmax).or.(deltaD.lt.deltaDmin)) then
          endif !if (q(i,k,iq).gt.qmin) then

       
        return
        end        

