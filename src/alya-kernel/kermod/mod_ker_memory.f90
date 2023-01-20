!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kermod
!> @{
!> @file    mod_ker_memory.f90
!> @author  houzeaux
!> @date    2020-11-13
!> @brief   Memory
!> @details Memory allocate and deallocate
!-----------------------------------------------------------------------

module mod_ker_memory

  use def_kintyp,      only : ip,rp,lg
  use def_kermod,      only : exch_loc_elemv
  use mod_memory,      only : memory_alloca
  use mod_memory,      only : memory_deallo
  use mod_memory_tools
  
  private

  integer(4) :: istat
  
  interface memory_alloca
     module procedure &
          &           memory_alloca_rexch_loc_elemv
  end interface memory_alloca

  interface memory_deallo
     module procedure &
          &           memory_deallo_rexch_loc_elemv
  end interface memory_deallo

  public :: memory_alloca
  public :: memory_deallo 

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-10
  !> @brief   Slip wall
  !> @details Allotate slip wall
  !> 
  !-----------------------------------------------------------------------
  
  subroutine memory_alloca_rexch_loc_elemv(memor,vanam,vacal,varia,ndim1,REALLOCATE)

    character(*),                    intent(in)     :: vanam         !< Variable name
    character(*),                    intent(in)     :: vacal         !< Calling subroutine name
    integer(8),                      intent(inout)  :: memor(2)      !< Memory counter
    type(exch_loc_elemv), pointer,   intent(inout)  :: varia(:)
    integer(ip),                     intent(in)     :: ndim1
    logical(lg),          optional,  intent(in)     :: REALLOCATE
    integer(ip)                                     :: idim1

    if( ndim1 > 0 ) then

       if( present(REALLOCATE) ) then
          if( REALLOCATE ) then
             call memory_deallo_rexch_loc_elemv(memor,vanam,vacal,varia)
          end if
       end if

       if( kfl_alloc == 1 ) call memory_deallo_rexch_loc_elemv(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       allocate( varia(ndim1) , stat = istat )

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(ip,8)
          do idim1 = lbound(varia,1),ubound(varia,1)
             varia(idim1) % nbogp   = 0
             varia(idim1) % fact    = 0.0_rp
             varia(idim1) % vel_aux = 0.0_rp
             varia(idim1) % velav   = 0.0_rp                
          end do
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'exch_loc_elemv')

    else

       nullify(varia)

    end if
    
  end subroutine memory_alloca_rexch_loc_elemv

  subroutine memory_deallo_rexch_loc_elemv(memor,vanam,vacal,varia)

    character(*),                    intent(in)     :: vanam         !< Variable name
    character(*),                    intent(in)     :: vacal         !< Calling subroutine name
    integer(8),                      intent(inout)  :: memor(2)      !< Memory counter
    type(exch_loc_elemv), pointer,   intent(inout)  :: varia(:)

       if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(ip,8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'exch_loc_elemv')

    else

       lbytm = 0
      
    end if

  end subroutine memory_deallo_rexch_loc_elemv
  
end module mod_ker_memory
