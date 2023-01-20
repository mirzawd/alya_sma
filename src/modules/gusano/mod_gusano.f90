!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Gusano
!> @{
!> @file    mod_gusano.f90
!> @author  houzeaux
!> @date    2020-04-22
!> @brief   Module for Gusano
!> @details Main gusano module with some useful tools
!-----------------------------------------------------------------------

module mod_gusano

  use def_master
  use def_domain
  use def_gusano
  use mod_memory
  
  implicit none
  private

  public :: gusano_main
  public :: gusano_memory_allocate
  public :: gusano_memory_deallocate
  public :: gusano_memory_reallocate

  character(10), parameter :: vacal = 'mod_gusano'
  
contains

  subroutine gusano_main(order)
    implicit none
    integer(ip), intent(in) :: order
    call Gusano(order)
  end subroutine gusano_main

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-09-22
  !> @brief   Reallocate array
  !> @details Reallocate a gusano array
  !> 
  !-----------------------------------------------------------------------

  subroutine gusano_memory_reallocate(wvari,NUMBER1,NUMBER2)
    character(len=*), intent(in)           :: wvari
    integer(ip),      intent(in), optional :: NUMBER1
    integer(ip),      intent(in), optional :: NUMBER2
    call gusano_memory_deallocate(wvari,NUMBER1,NUMBER2)
    call gusano_memory_allocate  (wvari,NUMBER1,NUMBER2)
  end subroutine gusano_memory_reallocate
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-09-22
  !> @brief   Allocate array
  !> @details Allocate a gusano array
  !> 
  !-----------------------------------------------------------------------

  subroutine gusano_memory_allocate(wvari,NUMBER1,NUMBER2)
    character(len=*), intent(in)           :: wvari
    integer(ip),      intent(in), optional :: NUMBER1
    integer(ip),      intent(in), optional :: NUMBER2
    call gusano_memory_allocate_deallocate(0_ip,wvari,NUMBER1,NUMBER2)
  end subroutine gusano_memory_allocate
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-09-22
  !> @brief   Deallocate array
  !> @details Deallocate a gusano array
  !> 
  !-----------------------------------------------------------------------

  subroutine gusano_memory_deallocate(wvari,NUMBER1,NUMBER2)
    character(len=*), intent(in)           :: wvari
    integer(ip),      intent(in), optional :: NUMBER1
    integer(ip),      intent(in), optional :: NUMBER2
    call gusano_memory_allocate_deallocate(1_ip,wvari,NUMBER1,NUMBER2)
  end subroutine gusano_memory_deallocate


  subroutine gusano_memory_allocate_deallocate(itask,wvari,NUMBER1,NUMBER2)

    integer(ip),      intent(in)           :: itask
    character(len=*), intent(in)           :: wvari
    integer(ip),      intent(in), optional :: NUMBER1
    integer(ip),      intent(in), optional :: NUMBER2

    select case(trim(wvari))

    case ( 'BOUNDARY CONDITIONS' )
       !
       ! Arrays read in reageo 
       !
       if( itask == 0 ) then
          call memory_alloca( mem_modul(1:2,modul),'KFL_FIXNO_GUS' ,vacal,kfl_fixno_gus , 2_ip,npoin       )   
          call memory_alloca( mem_modul(1:2,modul),'KFL_FIXBO_GUS' ,vacal,kfl_fixbo_gus , nboun            )            
          call memory_alloca( mem_modul(1:2,modul),'BVESS_GUS'     ,vacal,bvess_gus     , 2_ip, npoin,2_ip )
          call memory_alloca( mem_modul(1:2,modul),'BVNAT_GUS'     ,vacal,bvnat_gus     , 1_ip, nboun,2_ip )
          call memory_alloca( mem_modul(1:2,modul),'KFL_FUNNO_GUS' ,vacal,kfl_funno_gus , npoin            )
          call memory_alloca( mem_modul(1:2,modul),'KFL_FUNBO_GUS' ,vacal,kfl_funbo_gus , nboun            )
          call memory_alloca( mem_modul(1:2,modul),'KFL_FUNTN_GUS' ,vacal,kfl_funtn_gus , npoin            )
          call memory_alloca( mem_modul(1:2,modul),'KFL_FUNTB_GUS' ,vacal,kfl_funtb_gus , nboun            )
       else
          call memory_deallo( mem_modul(1:2,modul),'KFL_FIXNO_GUS' ,vacal,kfl_fixno_gus                    )   
          call memory_deallo( mem_modul(1:2,modul),'KFL_FIXBO_GUS' ,vacal,kfl_fixbo_gus                    )            
          call memory_deallo( mem_modul(1:2,modul),'BVESS_GUS'     ,vacal,bvess_gus                        )
          call memory_deallo( mem_modul(1:2,modul),'BVNAT_GUS'     ,vacal,bvnat_gus                        )
          call memory_deallo( mem_modul(1:2,modul),'KFL_FUNNO_GUS' ,vacal,kfl_funno_gus                    )
          call memory_deallo( mem_modul(1:2,modul),'KFL_FUNBO_GUS' ,vacal,kfl_funbo_gus                    )
          call memory_deallo( mem_modul(1:2,modul),'KFL_FUNTN_GUS' ,vacal,kfl_funtn_gus                    )
          call memory_deallo( mem_modul(1:2,modul),'KFL_FUNTB_GUS' ,vacal,kfl_funtb_gus                    )
       end if
       
   case ( 'VELOC' )
       !
       ! Arrays read in reageo 
       !
       if( itask == 0 ) then
          call memory_alloca( mem_modul(1:2,modul),'VELOC' ,vacal,veloc,ndime,npoin,1_ip )   
       else
          call memory_deallo( mem_modul(1:2,modul),'VELOC' ,vacal,veloc                  )     
       end if

    case ( 'EXTERIOR NORMAL','EXN1D_GUS' )
       !
       ! Exterior normal 
       !
       if( itask == 0 ) then
          call memory_alloca(mem_modul(1:2,modul),'EXN1D_GUS' ,'gus_memall',exn1d_gus ,npoin ) 
       else
          call memory_deallo(mem_modul(1:2,modul),'EXN1D_GUS' ,'gus_memall',exn1d_gus        )
       end if
       
    end select

  end subroutine gusano_memory_allocate_deallocate
  
end module mod_gusano
!> @}
