!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Coupling
!> @{
!> @name    Coupling arrays
!> @file    cou_memory.f90
!> @author  Guillaume Houzeaux
!> @date    03/03/2014
!> @brief   Allocate memory
!> @details Allocate memory and nullify pointers
!> @}
!------------------------------------------------------------------------

subroutine cou_memory(itask)
  use def_kintyp,          only : ip,rp
  use def_domain,          only : npoin
  use def_coupli,          only : mcoup
  use def_coupli,          only : mask_cou
  use def_coupli,          only : coupling_type
  use def_coupli,          only : memor_cou
  use def_coupli,          only : ncoup_implicit_n
  use def_coupli,          only : ncoup_implicit_d
  use def_coupli,          only : lcoup_implicit_n
  use def_coupli,          only : lcoup_implicit_d
  use def_coupli,          only : scala_cou
  use def_coupli,          only : nscal_cou
  use mod_coupling_memory, only : cou_initialization  
  use mod_memory,          only : memory_alloca
  use mod_memory,          only : memory_deallo
  implicit none 
  integer(ip), intent(in) :: itask

  select case ( itask )

  case ( 0_ip ) 

     call memory_deallo(memor_cou,'MASK_COU'        ,'cou_deallocate',mask_cou)
     call memory_deallo(memor_cou,'LCOUP_IMPLICIT_N','cou_deallocate',lcoup_implicit_n)
     call memory_deallo(memor_cou,'LCOUP_IMPLICIT_D','cou_deallocate',lcoup_implicit_d)     
     if( associated(coupling_type) ) deallocate( coupling_type )
 
  case ( 1_ip ) 

     allocate( coupling_type(mcoup) )
     call cou_initialization(coupling_type)

  case ( 2_ip ) 

     call memory_alloca(memor_cou,'MASK_COU','cou_memory',mask_cou,npoin)

  case ( 3_ip ) 

     call memory_alloca(memor_cou,'LCOUP_IMPLICIT_N','cou_memory',lcoup_implicit_n,ncoup_implicit_n)
     call memory_alloca(memor_cou,'LCOUP_IMPLICIT_D','cou_memory',lcoup_implicit_d,ncoup_implicit_d)

  case ( 4_ip ) 

     if( nscal_cou > 0 ) & 
          call memory_alloca(memor_cou,'SCALA_COU','cou_memory',scala_cou,mcoup,2_ip)
     
  end select

end subroutine cou_memory
