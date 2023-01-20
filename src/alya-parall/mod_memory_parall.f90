!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!>
!> @defgroup Memory_Toolbox
!> @{
!> @name    ToolBox for memory management
!> @file    mod_memory_parall.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for memory management
!> @details Tools for memory of variables defined in def_kintyp_parall
!>
!------------------------------------------------------------------------

module mod_memory_parall

  use def_kintyp,        only : ip,rp,lg
  use def_kintyp_comm,   only : tAdj_par
   
  use mod_memory_tools
  use mod_memory_basic

  implicit none

  private

  interface memory_alloca
     module procedure &
          &           memory_alloca_tAdj_par
  end interface memory_alloca

  interface memory_deallo
     module procedure &
          &           memory_deallo_tAdj_par
  end interface memory_deallo

  public :: memory_alloca
  public :: memory_deallo

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-10
  !> @brief   Allotate edge structure
  !> @details Allotate edge structure
  !> 
  !-----------------------------------------------------------------------

  subroutine memory_alloca_tAdj_par(memor,vanam,vacal,varia,ndim1,wzero)
    !
    ! tadj_par
    !  
    character(*),              intent(in)    :: vanam         !< Variable name
    character(*),              intent(in)    :: vacal         !< Calling subroutine name
    integer(8),                intent(inout) :: memor(2)      !< Memory counter
    integer(ip),               intent(in)    :: ndim1
    type(tAdj_par), pointer                  :: varia(:)
    character(*),   optional,  intent(in)    :: wzero
    integer(4)                               :: istat
    integer(ip)                              :: idim1
    logical(lg)                              :: lzero

    if( ndim1 > 0 ) then

       if( kfl_alloc == 1 )    call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       allocate( varia(ndim1) , stat = istat )
       
       lzero = .true.
       if( present(wzero) ) then
          if(      trim(wzero) == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*ip
          if( lzero ) then
             do idim1 = 1,ndim1
                varia(idim1) % node1 = 0_ip
                varia(idim1) % node2 = 0_ip
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'tadj_par')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_tAdj_par

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-10
  !> @brief   Deallotate edge structure
  !> @details Deallotate edge structure
  !> 
  !-----------------------------------------------------------------------

  subroutine memory_deallo_tAdj_par(memor,vanam,vacal,varia)
    !
    ! tadj_par
    !
    character(*),   intent(in)            :: vanam         !< Variable name
    character(*),   intent(in)            :: vacal         !< Calling subroutine name
    integer(8),     intent(inout)         :: memor(2)      !< Memory counter
    type(tAdj_par), pointer               :: varia(:)
    integer(4)                            :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*ip

       deallocate( varia , stat = istat )

       nullify(varia)

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)

       end if
       call memory_info(memor,vanam,vacal,'tadj_par')

    else

       lbytm = 0
       
    end if

  end subroutine memory_deallo_tAdj_par
 
end module mod_memory_parall
!> @}
