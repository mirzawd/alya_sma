!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!>
!> @addtogroup Domain
!> ToolBox for elemental data base
!> @{
!> @file    mod_element_data_base.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for elements data base
!> @details Different functions to manage element data base
!> @{
!>
!------------------------------------------------------------------------

module mod_element_data_base

  use def_parame
  use def_master
  use def_domain
  use def_kermod
  use mod_memchk
  use mod_memory
  use mod_messages, only : messages_live
  implicit none
  private

  integer(8) :: memor_db(2)

  interface element_data_base_memory
     module procedure &
          & element_data_base_memory_4, &
          & element_data_base_memory_8
  end interface element_data_base_memory

  public :: element_data_base_memory
  public :: element_data_base_save
  public :: element_data_base_initialization

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-05-28
  !> @brief   Element data base deallocate
  !> @details Element data base deallocate
  !> 
  !-----------------------------------------------------------------------

  subroutine element_data_base_deallocate()

    integer(ip) :: ielem
    
    if( associated(elmda) ) then
       do ielem = 1,nelem
          if( associated(elmda(ielem) % gpcar) ) deallocate(elmda(ielem) % gpcar)
          if( associated(elmda(ielem) % gpvol) ) deallocate(elmda(ielem) % gpvol)
          if( associated(elmda(ielem) % hleng) ) deallocate(elmda(ielem) % hleng)
          if( associated(elmda(ielem) % tragl) ) deallocate(elmda(ielem) % tragl)
          if( associated(elmda(ielem) % gphes) ) deallocate(elmda(ielem) % gphes)
       end do
       deallocate(elmda)
    end if

  end subroutine element_data_base_deallocate

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-05-28
  !> @brief   Element data base initialization
  !> @details Element data base initialization
  !> 
  !-----------------------------------------------------------------------

  subroutine element_data_base_initialization()

    nullify(elmda)
    nullify(elmda_gpvol)
    nullify(elmda_gpcar)

  end subroutine element_data_base_initialization

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-05-28
  !> @brief   Element data base memory
  !> @details Element data base memory
  !> 
  !-----------------------------------------------------------------------

  subroutine element_data_base_memory_8(memory)
    integer(8), intent(out) :: memory
    if( kfl_savda /= 0 ) then
       memory = memor_db(1)
    else
       memory = 0_8
    end if
  end subroutine element_data_base_memory_8

  subroutine element_data_base_memory_4(memory)
    integer(4), intent(out) :: memory
    if( kfl_savda /= 0 ) then
       memory = int(memor_db(1),4)
    else
       memory = 0_4
    end if
  end subroutine element_data_base_memory_4

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-05-28
  !> @brief   Save element data base
  !> @details This routine save the element data base:
  !>          GPCAR: Cartesian derivatives
  !>          GPVOL: Unit volume
  !>          GPHES: Hessian matrix
  !> 
  !-----------------------------------------------------------------------

  subroutine element_data_base_save()

    integer(ip) :: pelty,pnode,pgaus,plapl
    real(rp)    :: elcod(ndime,mnode)
    integer(ip) :: inode,idime,ipoin,ielem
    real(rp)    :: dummr(1,1,1)

    memor_db = 0_8

    if( INOTMASTER .and. kfl_savda == 1 ) then

       !-----------------------------------------------------------------
       !
       ! Element data base per element
       !
       !-----------------------------------------------------------------

       call  messages_live('SAVE ELEMENT DATA BASE')

       if( associated(elmda) ) then
          call element_data_base_deallocate()
       end if

       allocate(elmda(nelem))
       do ielem = 1,nelem
          nullify(elmda(ielem) % gpcar)
          nullify(elmda(ielem) % gpvol)
          nullify(elmda(ielem) % hleng)
          nullify(elmda(ielem) % tragl)
          nullify(elmda(ielem) % gphes)
       end do

       elements: do ielem = 1,nelem
          ! 
          ! Element properties and dimensions
          !
          pelty = ltype(ielem)
          pnode = nnode(pelty)
          pgaus = ngaus(pelty)
          plapl = llapl(pelty)
          !
          ! Gather
          !
          do inode = 1,pnode
             ipoin = lnods(inode,ielem)
             do idime = 1,ndime
                elcod(idime,inode) = coord(idime,ipoin)
             end do
          end do
          !
          ! Allocate memory
          !
          call memory_alloca(memor_db,'GPCAR','memgeo',elmda(ielem) % gpcar,ndime,mnode,pgaus)
          call memory_alloca(memor_db,'GPVOL','memgeo',elmda(ielem) % gpvol,pgaus)
          call memory_alloca(memor_db,'HLENG','memgeo',elmda(ielem) % hleng,ndime)
          call memory_alloca(memor_db,'TRAGL','memgeo',elmda(ielem) % tragl,ndime,ndime)
          !
          ! HLENG and TRAGL at center of gravity
          !
          call elmlen(&
               ndime,pnode,elmar(pelty)%dercg,elmda(ielem)%tragl,elcod,&
               hnatu(pelty),elmda(ielem)%hleng)

          if( plapl == 1 ) then
             !
             ! With Hessian: GPVOL, GPCAR, GPHES
             !
             call memory_alloca(memor_db,'TRAGL','memgeo',elmda(ielem) % gphes,ntens,mnode,pgaus)

             call elmcar(&
                  pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
                  elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,elmda(ielem) % gpvol,&
                  elmda(ielem) % gpcar,elmda(ielem) % gphes,ielem)  

          else
             !
             ! Without Hessian: GPVOL, GPCAR
             !
             call elmcar(&
                  pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
                  elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,elmda(ielem) % gpvol,&
                  elmda(ielem) % gpcar,dummr,ielem)

          end if

       end do elements

       if( kfl_data_base_array == 1 ) then
          call memory_alloca(memor_db,'ELMDA_GPVOL','savdom',elmda_gpvol,mgaus,nelem)
          call memory_alloca(memor_db,'ELMDA_GPCAR','savdom',elmda_gpcar,ndime,mnode,mgaus,nelem)
          do ielem = 1,nelem
             pelty = ltype(ielem)
             pnode = nnode(pelty)
             pgaus = ngaus(pelty)
             elmda_gpvol(1:pgaus,ielem)                 = elmda(ielem) % gpvol(1:pgaus)
             elmda_gpcar(1:ndime,1:mnode,1:pgaus,ielem) = elmda(ielem) % gpcar(1:ndime,1:mnode,1:pgaus)
          end do
       end if

       kfl_savda    = 2

    end if

    if( kfl_savda /= 0 ) then

       memor_dom(1) = memor_dom(1) + memor_db(1)
       memor_dom(2) = max(memor_dom(2),memor_dom(1))

    end if

  end subroutine element_data_base_save

end module mod_element_data_base
!> @}
