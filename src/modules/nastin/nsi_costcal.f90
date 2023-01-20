!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!> nsi_costcal.f90
!> @file nsi_costcal.f90 
!> @fn nsi_costcal 
!> This subroutine calculates the costf by summing between subdomain and demostrating
!>

subroutine nsi_costcal

  use def_parame
  use def_elmtyp
  use def_master
  use def_kermod, only : costf,kfl_cost_type
!  use def_kermod, only : kfl_adj_prob
  use def_domain
  use mod_communications, only : PAR_SUM
  use def_nastin
  implicit none
  
  real(rp)    :: cd,cl,cdv,clv

  
!   if( kfl_stead_nsi == 1 ) then
  
    if (kfl_cost_type == 2) then
      call PAR_SUM(costf,  'IN MY CODE' )
      print*, "costf nsi", costf, (costf-5.94007169879792_rp)/0.000001_rp !!!!!!
    endif                            


    if (kfl_cost_type == 5) then 
    
      
      if (INOTSLAVE) then
      
        cl = (0.96_rp*(vbset(7,2)+vbset(7,3)+vbset(7,4)) - 0.28_rp*(vbset(6,2)+vbset(6,3)+vbset(6,4)))
        cd = -(0.96_rp*(vbset(6,2)+vbset(6,3)+vbset(6,4)) + 0.28_rp*(vbset(7,2)+vbset(7,3)+vbset(7,4)))

        clv = (0.96_rp*(vbset(10,2)+vbset(10,3)+vbset(10,4)) - 0.28_rp*(vbset(9,2)+vbset(9,3)+vbset(9,4)))
        cdv = -(0.96_rp*(vbset(9,2)+vbset(9,3)+vbset(9,4)) + 0.28_rp*(vbset(10,2)+vbset(10,3)+vbset(10,4)))
        
        open(10,file='functional.dat')
        open(11,file='draglift.dat')
      
        write (11,*) cl, cd
        write (10,*) cd
        
        close(10)
        close(11)
      endif
    
      if (INOTSLAVE) then
!          print *, "drag f. d. -10",vbset(6,3), (vbset(6,3) - 118.885959488902_rp)/0.00000001_rp
!          print *, "lift f. d. -10",vbset(7,3), (vbset(7,3) + 513.625124481609_rp)/0.00000001_rp
!          print *, "drag f. d.",vbset(6,3), (vbset(6,3) - 118.885997071852_rp)/0.000001_rp
!          print *, "lift f. d.",vbset(7,3), (vbset(7,3) + 513.625245178594_rp)/0.000001_rp
!          print *, "drag f. d. walli -8",vbset(6,3), -(vbset(6,3) - 116.822514733434_rp)/0.000001_rp
!          print *, "lift f. d. walli -8",vbset(7,3), -(vbset(7,3) + 506.489277314636_rp)/0.000001_rp
!          print *, "drag f. d. walli -12",vbset(6,3), (vbset(6,3) - 116.823229885800_rp)/0.00000001_rp
!          print *, "lift f. d. walli -12",vbset(7,3), (vbset(7,3) + 506.491404530935_rp)/0.00000001_rp
!          print *, "drag f. d. walli -10 ramp",vbset(6,4), (vbset(6,4) + 1.47034473542032_rp)/0.00000001_rp
!          print *, "lift f. d. walli -10 ramp",vbset(7,4), (vbset(7,4) - 98.0229820727916_rp)/0.00000001_rp
!          print *, "drag f. d. walli -10 ramp all",vbset(6,4), (vbset(6,4) + 1.48040209782823_rp)/0.00000001_rp
!          print *, "lift f. d. walli -10 ramp all",vbset(7,4), (vbset(7,4) - 98.6934728985589_rp)/0.00000001_rp
!          print *, "cylinder2d drag f. d. -10", vbset(6,5) + vbset(6,6), -(vbset(6,5) + vbset(6,6) + 5001570.34108898_rp)/0.001_rp
!          print *, "cylinder2d lift f. d. -10", vbset(7,5) + vbset(7,6), -(vbset(7,5) + vbset(7,6) + 14342.3261879901_rp)/0.001_rp
!          print *, "xforce NACA3D",vbset(6,4)!, (vbset(6,4) + 5324692.91402950_rp)/0.01_rp
!          print *, "zforce NACA3D",vbset(8,4)!, (vbset(8,4) - 1038696.26150939_rp)/0.01_rp
!         print *, "local x-force 30P",cd
!         print *, "local y-force 30P",cl
         
      endif

    endif
  
!   endif
 
  
!   
!   if (kfl_cost_type == 5 .and. kfl_adj_prob == 1) then
!     print *, "vbset(6,1:2)", -vbset(6,1) - vbset(6,2)
!   endif
  

  
end subroutine nsi_costcal
