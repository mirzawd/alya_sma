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
!> @file    mod_memory.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for memory management
!> @details ToolBox for memory management: allocate, deallocate
!>          for basic types defined in def_kintyp_basic
!>
!------------------------------------------------------------------------

module mod_memory_spmat

  use def_kintyp_basic, only : ip,rp,lg,i1p,i2p,i3p,r1p,r2p,r3p,r4p,i1pp
  use def_spmat,        only : spmat
  use mod_memory_basic, only : memory_alloca
  use mod_memory_basic, only : memory_deallo
  use mod_memory_tools, only : memory_error
  use mod_memory_tools, only : memory_info
  use mod_memory_tools, only : memory_already_associated
  use mod_memory_tools, only : lbytm
  use mod_memory_tools, only : kfl_alloc
  implicit none

  private
  !
  ! MEMORY_ALLOC: allocate
  !
  interface memory_alloca
     !                      (:)                    (:,:)                 (:,:,:)               (:,:,:,:)
     !                --------------------------------------------------------------------------------------------
     module procedure memory_alloca_spmat  , memory_alloca_spmat_1                                                  ! SPMAT
  end interface memory_alloca

  interface memory_deallo
     !                      (:)                    (:,:)                 (:,:,:)               (:,:,:,:)
     !                -------------------------------------------------------------------------------------------
     module procedure memory_deallo_spmat  , memory_deallo_spmat_1                                                  ! SPMAT
  end interface memory_deallo

  public :: memory_alloca
  public :: memory_deallo

contains

  subroutine memory_alloca_spmat_1(memor,vanam,vacal,varia,ndim1,wzero,lboun)
    !
    ! Type(spmat)
    !

    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    type(spmat),  intent(inout), pointer  :: varia(:)
    integer(ip),  intent(in)              :: ndim1
    character(*), intent(in),    optional :: wzero
    integer(ip),  intent(in),    optional :: lboun
    integer(ip)                           :: idim1
    integer(4)                            :: istat
    logical(lg)                           :: lzero

    if( ndim1 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       istat=0
       if( present(lboun) ) then
          allocate( varia(lboun:lboun+ndim1-1) , stat = istat )
       else
          allocate( varia(ndim1) , stat = istat )
       end if

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*ip
          if( lzero ) then
             do idim1 = lbound(varia,1_ip),ubound(varia,1_ip)
                call varia(idim1) % init()
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type(spmat)')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_spmat_1

  subroutine memory_alloca_spmat(memor,vanam,vacal,varia,ndof1,ndim1,ndof2,wzero)
    !
    ! Type(spmat)
    !
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    type(spmat), intent(out)            :: varia
    integer(ip),  intent(in)            :: ndof1
    integer(ip),  intent(in)            :: ndim1
    integer(ip),  intent(in),  optional :: ndof2
    character(*), intent(in),  optional :: wzero
    integer(4)                          :: istat
    integer(4)                          :: aux_istat
    logical(lg)                         :: lzero
    integer(ip)                         :: ndof

    if( ndim1 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia % iA) .or. associated(varia % jA) .or. associated(varia % vA) ) then
          call memory_already_associated(vanam,vacal)
       end if

       if( present(ndof2) ) then
          ndof = ndof2
       else
          ndof = ndof1
       end if

       !call varia % alloca(ndof1,ndof,ndim1,MEMORY_COUNTER=memor)
       istat=0
       allocate( varia % iA(ndim1) , stat = aux_istat )
       istat = istat + aux_istat
       allocate( varia % jA(ndim1) , stat = aux_istat )
       istat = istat + aux_istat
       allocate( varia % vA(ndof1, ndof, ndim1) , stat = aux_istat )
       istat = istat + aux_istat
       lzero = .true.
       
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia % iA,kind=8) * ip + size(varia % jA,kind=8) * ip + size(varia % vA,kind=8) * rp 
          if( lzero ) then
             varia % nz    = ndim1
             varia % ndof1 = ndof1
             varia % ndof2 = ndof
             varia % nrows = 0
             varia % ncols = 0
             varia % iA(1:ndim1) = 0_ip
             varia % jA(1:ndim1) = 0_ip
             varia % vA(1:ndof1,1:ndof,1:ndim1) = 0_ip
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type(spmat)')

    else

       call varia % init()

    end if

  end subroutine memory_alloca_spmat

  subroutine memory_deallo_spmat_1(memor,vanam,vacal,varia)
    !
    ! type(spmat)
    !

    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(spmat), pointer                :: varia(:)
    integer(4)                          :: istat
    integer(ip)                         :: idim1
    if( associated(varia) ) then

       do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
          call memory_deallo_spmat( memor, vanam, vacal, varia(idim1) )
       end do
       lbytm = -size(varia,kind=8)*ip
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'real')
    end if

  end subroutine memory_deallo_spmat_1

  subroutine memory_deallo_spmat(memor,vanam,vacal,varia)
    !
    ! type(spmat)
    !

    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(spmat)                         :: varia
    integer(4)                          :: istat
    integer(4)                          :: istat_aux

    istat = 0
    varia % ndof1 = 0
    varia % ndof2 = 0
    varia % nrows = 0
    varia % ncols = 0

    if( associated(varia % iA) .or.  associated(varia % jA) .or. associated(varia % vA) ) then

       lbytm = 0
       if(associated(varia % iA) ) then
          lbytm = lbytm - size(varia % iA,kind=8) * ip 
          deallocate(varia % iA, stat=istat_aux)
          istat = istat + istat_aux
       endif
       if(associated(varia % jA) ) then
          lbytm = lbytm - size(varia % jA,kind=8) * ip 
          deallocate(varia % jA, stat=istat_aux)
          istat = istat + istat_aux
       endif
       if(associated(varia % vA) ) then
          lbytm = lbytm - size(varia % vA,kind=8) * rp 
          deallocate(varia % vA, stat=istat_aux)
          istat = istat + istat_aux
       endif

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia % iA)
          nullify (varia % jA)
          nullify (varia % vA)
       end if
       call memory_info(memor,vanam,vacal,'real')
    end if

  end subroutine memory_deallo_spmat

end module mod_memory_spmat
!> @}
