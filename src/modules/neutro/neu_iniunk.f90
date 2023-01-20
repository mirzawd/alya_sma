!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Neutro
!> @{
!> @file    neu_iniunk.f90
!> @author  Guillaume Houzeaux
!> @brief   Initial value 
!> @details Set up the initial condition
!> @} 
!------------------------------------------------------------------------
subroutine neu_iniunk()
  use def_parame 
  use def_master
  use def_elmtyp
  use def_kermod
  use def_domain
  use def_neutro 
  use mod_ker_space_time_function
  implicit none
  integer(ip) :: idirection,ienergy,ipoin
  external :: neu_updunk

  if( kfl_rstar == 0 ) then  
     !
     ! Cheak if there are weak no-penetration conditions (Y esto????)
     !
     do idirection = 1,num_directions_neu
        do ienergy = 1,num_energies_neu
           do ipoin = 1,npoin
              neutr(ienergy,idirection,ipoin,nprev_neu) = bvess_neu(ienergy,idirection) % a(1,ipoin,1)
           end do
        end do
     end do
     !
     ! NEUTRONS(:,:,:,1) <= NEUTRONS(:,:,:,nprev_neu)
     !
     call neu_updunk(11_ip) 

  end if 


end subroutine neu_iniunk
