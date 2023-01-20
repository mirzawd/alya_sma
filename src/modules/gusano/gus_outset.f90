!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Gusano
!> @{
!> @file    gus_outset.f90
!> @author  houzeaux
!> @date    2020-10-26
!> @brief   Output on sets
!> @details Output on sets
!> @} 
!-----------------------------------------------------------------------

subroutine gus_outset()

  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_gusano
  use mod_iofile
  use mod_communications, only : PAR_BROADCAST
  use mod_output_postprocess, only : output_postprocess_boundary_sets_parall

  implicit none
  
  integer(ip) :: kbset,iboun,ipoin

  !----------------------------------------------------------------------
  !
  ! Boundary sets
  !
  !----------------------------------------------------------------------

  if( maxval(postp(1) % npp_setsb) > 0 ) then

     do kbset = 1,nbset 
        do iboun = 1,nboun
           if( lbset(iboun) == lbsec(kbset) ) then
              ipoin = lnodb(1,iboun)
              if( postp(1) % npp_setsb(1) /= 0 ) then
                 vbset(1,kbset) = press(ipoin,1)
              end if
              if( postp(1) % npp_setsb(2) /= 0 ) then
                 vbset(2,kbset) = flowr(ipoin,1) * exn1d_gus(ipoin)
              end if
           end if
        end do 
     end do
     call output_postprocess_boundary_sets_parall() 
 
  end if

end subroutine gus_outset
