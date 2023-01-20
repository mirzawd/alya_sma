!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Gusano
!> @{
!> @file    gus_inibcs.f90
!> @author  Guillaume Houzeaux
!> @brief   Impose boundary conditions
!> @details Impose boundary conditions\n
!> @} 
!-----------------------------------------------------------------------
subroutine gus_inibcs()
  use def_parame
  use def_inpout
  use def_master
  use def_kermod
  use def_domain
  use def_gusano
  use mod_communications,                only : PAR_MAX, PAR_INTERFACE_NODE_EXCHANGE
  use mod_gusano,                        only : gusano_memory_allocate
  use mod_communications_point_to_point, only : PAR_INTERFACE_NODE_EXCHANGE
  implicit none
  integer(ip) :: ielem,ipoin

  call gusano_memory_allocate('BOUNDARY CONDITIONS') 
  call gusano_memory_allocate('EXTERIOR NORMAL') 
  !
  ! Compute 1d normal (only required for end points)
  !
  do ielem = 1,nelem
     if( ltype(ielem) > 0 ) then
        ipoin            =  lnods(1,ielem)
        exn1d_gus(ipoin) = -1.0_rp
        ipoin            =  lnods(2,ielem)
        exn1d_gus(ipoin) =  1.0_rp
     end if
  end do
  call PAR_INTERFACE_NODE_EXCHANGE(exn1d_gus,'SUM')

  if( INOTMASTER ) then
     !
     ! Node codes
     ! 
     if( kfl_icodn > 0 ) then
        iffun     =  1
        kfl_fixno => kfl_fixno_gus
        kfl_funno => kfl_funno_gus
        kfl_funtn => kfl_funtn_gus
        bvess     => bvess_gus(:,:,1)
        tncod     => tncod_gus(1:)
        call reacod(IMPOSE_NODE_CODES)
        do ipoin = 1,npoin
           bvess_gus(1,ipoin,1) = exn1d_gus(ipoin) * bvess_gus(1,ipoin,1)
        end do
     end if
     !
     ! Boundary codes
     !
     if( kfl_icodb > 0 ) then
        iffun     =  1
        kfl_fixbo => kfl_fixbo_gus
        kfl_funbo => kfl_funbo_gus
        kfl_funtb => kfl_funtb_gus
        bvnat     => bvnat_gus(:,:,1)
        tbcod     => tbcod_gus(1:)
        call reacod(IMPOSE_BOUNDARY_CODES)
     end if
  end if

end subroutine gus_inibcs
