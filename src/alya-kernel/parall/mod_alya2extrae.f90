!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> Brige to DLB
!> @{
!> @file    mod_alya2extrae.f90
!> @author  houzeaux
!> @date    2019-01-07
!> @brief   Bridge to EXTRAE
!> @details Interfaces with EXTRAE
!-----------------------------------------------------------------------

module mod_alya2extrae

  use def_master
  use, intrinsic :: ISO_C_BINDING, only: C_CHAR, C_NULL_CHAR, C_PTR, C_LOC
  implicit none
  private

  public :: alya2extrae_initialization

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  bsc21943
  !> @date    2022-05-13
  !> @brief   Define descriptor for EXTRAE values
  !> @details Define descriptor for EXTRAE values
  !> 
  !-----------------------------------------------------------------------

  subroutine alya2extrae_initialization()

#ifdef ALYA_EXTRAE

    integer(8),                    dimension(mmodu*34+2)         :: values 
    character(KIND=C_CHAR,LEN=50), dimension(mmodu*34+2), target :: description_values
    type(C_PTR),                   dimension(mmodu*34+2)         :: description_values_ptrs
    character(KIND=C_CHAR,LEN=20)                                :: evt_desc = "Module execution" // C_NULL_CHAR
    integer(ip)                                                  :: itask,imodu,ii

    ii                          = 0
    !
    ! End of an event: 0
    !
    ii                          = ii + 1
    values(ii)                  = 0
    description_values(ii)      = "END " // C_NULL_CHAR
    description_values_ptrs(ii) = C_LOC(description_values(ii))
    !
    ! Coupling: 99
    !
    ii                          = ii + 1
    values(ii)                  = 99 
    description_values(ii)      = "COUPLING " // C_NULL_CHAR
    description_values_ptrs(ii) = C_LOC(description_values(ii))
    !
    ! Coupling init: 98
    !
    ii                          = ii + 1
    values(ii)                  = 98
    description_values(ii)      = "COUPLING INITIALIZATION" // C_NULL_CHAR
    description_values_ptrs(ii) = C_LOC(description_values(ii))

    do imodu = 1,mmodu
       if( kfl_modul(imodu) == 1 ) then
          !
          ! Module tasks: imodu*100+itask
          !
          do itask = 1,size(TASK_LONG_NAME)
             ii                          = ii + 1
             values(ii)                  = imodu*100+itask
             description_values(ii)      = namod(imodu) // ' ' // TASK_LONG_NAME(itask) // C_NULL_CHAR
             description_values_ptrs(ii) = C_LOC(description_values(ii))          
          end do
       end if
    end do
    call extrae_define_event_type (900, evt_desc, ii , values, description_values_ptrs)

#endif

  end subroutine alya2extrae_initialization

end module mod_alya2extrae
