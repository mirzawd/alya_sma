!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup GusanoInput
!> @{
!> @file    gus_updbcs.f90
!> @date    29/01/2018
!> @author  Guillaume Houzeaux
!> @brief   Update boundary conditions
!> @details This routine updates the velocity boundary conditions
!> @} 
!-----------------------------------------------------------------------
subroutine gus_updbcs(itask)

  use def_elmtyp
  use def_master
  use def_kermod
  use def_domain
  use def_gusano
  use mod_elmgeo,                        only : elmgeo_number_nodes
  use mod_communications_point_to_point, only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_gus_transmssion_conditions,    only : gus_transmssion_conditions_initialization
  implicit none
  
  integer(ip),         intent(in) :: itask
  integer(ip)                     :: ielem,ipoin,inode,pelty

  select case( itask )

  case( ITASK_TURNON )

     !-------------------------------------------------------------------
     !
     ! Initialization
     !
     !-------------------------------------------------------------------
     !
     ! Declare first node of a coupling element as the Dirichlet one
     !
     call memgen(1_ip,npoin,0_ip)
     do ielem = 1,nelem
        pelty = ltype(ielem)
        if( pelty == DDDNE ) then
           ipoin = lnods(1,ielem)
           !
           ! Dirichlet node
           !
           gisca(ipoin) = -1
           !
           ! Neumann nodes
           !
           do inode = 2,elmgeo_number_nodes(pelty,lnods(:,ielem))
              ipoin = lnods(inode,ielem)
              gisca(ipoin) = -2
           end do
        end if
     end do
     call PAR_INTERFACE_NODE_EXCHANGE(gisca,'SUM')
     do ipoin = 1,npoin
        if( gisca(ipoin) /= 0 ) kfl_fixno_gus(1,ipoin) = gisca(ipoin)
     end do
     call memgen(3_ip,npoin,0_ip)
     
  case( ITASK_BEGSTE , ITASK_INIUNK )

     !-------------------------------------------------------------------
     !  
     ! Before a time step
     !     
     !-------------------------------------------------------------------

     do ipoin = 1,npoin
        if( kfl_fixno_gus(1,ipoin) > 0 ) flowr(ipoin,1) = bvess_gus(1,ipoin,1)
        if( kfl_fixno_gus(2,ipoin) > 0 ) press(ipoin,1) = bvess_gus(2,ipoin,1)
     end do
     !
     ! Compute transmission arrays
     !
     if( itask == ITASK_INIUNK ) call gus_transmssion_conditions_initialization()
     !
     ! Schur complement fixity
     !
     if( kfl_algor_gus == GUS_SCHUR_COMPLEMENT ) then
        do ipoin = 1,npoin
           kfl_fixsc_gus(1,ipoin) = 0
           if( lpoty(ipoin) /= 0 ) then 
              if( kfl_fixno(1,ipoin) == 0 ) then
                 kfl_fixsc_gus(1,ipoin) = 1
              end if
           end if
        end do
     end if
     
  case( ITASK_BEGITE )

     !-------------------------------------------------------------------
     !
     ! Before a global iteration
     !  
     !-------------------------------------------------------------------

  case( ITASK_BEGINN )

     !-------------------------------------------------------------------
     !
     ! Before an inner iteration
     !  
     !-------------------------------------------------------------------     

  end select

end subroutine gus_updbcs



