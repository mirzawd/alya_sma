!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!> @defgroup Check_Type_Toolbox
!! @{
!> @file    mod_chktyp.f90
!> @author  Guillaume Houzeaux
!> @date    18/09/2012
!> @brief   Check the existence and size of fortran types
!> @details Check the existence and size of fortran types
!
!-----------------------------------------------------------------------
module mod_chktyp

  use mod_memchk
  use mod_memory
  use def_parame
  use def_kintyp
  use def_master
  use def_domain
  implicit none

  private

   interface check_type
      module procedure chktyp_r1ptyp,chktyp_r2ptyp,chktyp_r3ptyp
   end interface check_type

   public :: check_type

contains 

  !-----------------------------------------------------------------------
  !
  !> @brief   Check r1p type
  !> @details Check r1p type
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine chktyp_r1ptyp(xarra,nposi,ndim1,vanam)
    
    type(r1p),    intent(in), pointer  :: xarra(:)   !< Type
    integer(ip),  intent(in), optional :: nposi      !< XARRA(NPOSI) % A(:)
    integer(ip),  intent(in), optional :: ndim1      !< XARRA(NPOSI) % A(NDIM1)
    character(*), intent(in), optional :: vanam

    if( .not. associated(xarra) ) then
       call runend('CHKTYP_R1PTYP: TYPE WAS NOT CREATED')
    else 
       if( present(nposi) ) then 
          if( nposi <= 0 ) then
             call runend('CHKTYP_R1PTYP: TYPE IS WRONG') 
          end if
          if( size(xarra,KIND=ip) < nposi ) then
             call runend('CHKTYP_R1PTYP: TYPE WAS NOT CREATED') 
          end if
          if( .not. associated(xarra(nposi) % a) ) then
             call runend('CHKTYP_R1PTYP: TYPE WAS NOT CREATED')                 
          end if
          if( present(ndim1) ) then 
             if( size( xarra(nposi) % a ,KIND=ip) /= ndim1 ) then
                call runend('CHKTYP_R1PTYP: WRONG TYPE 1ST DIMENSION')             
             end if
          end if
       end if
    end if
    
  end subroutine chktyp_r1ptyp

  !-----------------------------------------------------------------------
  !
  !> @brief   Check r2p type
  !> @details Check r2p type
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine chktyp_r2ptyp(xarra,nposi,ndim1,ndim2)
    
    type(r2p),   intent(in), pointer  :: xarra(:)   !< Type
    integer(ip), intent(in), optional :: nposi      !< XARRA(NPOSI) % A(:,:)
    integer(ip), intent(in), optional :: ndim1      !< XARRA(NPOSI) % A(NDIM1,:)
    integer(ip), intent(in), optional :: ndim2      !< XARRA(NPOSI) % A(NDIM1,NDIM2)

    if( .not. associated(xarra) ) then
       call runend('CHKTYP_R2PTYP: TYPE WAS NOT CREATED')
    else 
       if( present(nposi) ) then 
          if( nposi <= 0 ) then
             call runend('CHKTYP_R2PTYP: TYPE IS WRONG') 
          end if
          if( size(xarra,KIND=ip) < nposi ) then
             call runend('CHKTYP_R2PTYP: TYPE WAS NOT CREATED') 
          end if
          if( .not. associated(xarra(nposi) % a) ) then
             call runend('CHKTYP_R2PTYP: TYPE WAS NOT CREATED')                 
          end if
          if( present(ndim1) ) then 
             if( size( xarra(nposi) % a,1,KIND=ip) /= ndim1 ) then
                call runend('CHKTYP_R2PTYP: WRONG TYPE 1ST DIMENSION')             
             end if
          end if
          if( present(ndim2) ) then 
             if( size( xarra(nposi) % a,2,KIND=ip) /= ndim2 ) then
                call runend('CHKTYP_R2PTYP: WRONG TYPE 2ND DIMENSION')             
             end if
          end if
       end if
    end if
    
  end subroutine chktyp_r2ptyp

  !-----------------------------------------------------------------------
  !
  !> @brief   Check r3p type
  !> @details Check r3p type
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine chktyp_r3ptyp(xarra,nposi,ndim1,ndim2,ndim3)

    type(r3p),   intent(in)           :: xarra(:)   !< Type
    integer(ip), intent(in), optional :: nposi      !< XARRA(NPOSI) % A(:,:)
    integer(ip), intent(in), optional :: ndim1      !< XARRA(NPOSI) % A(NDIM1,:)
    integer(ip), intent(in), optional :: ndim2      !< XARRA(NPOSI) % A(NDIM1,NDIM2)
    integer(ip), intent(in), optional :: ndim3      !< XARRA(NPOSI) % A(NDIM1,NDIM2)
    integer(ip)                       :: nsize

    nsize = 1
    if( present(ndim1) ) nsize = nsize * ndim1 
    if( present(ndim2) ) nsize = nsize * ndim2 
    if( present(ndim3) ) nsize = nsize * ndim3
    
    if( present(nposi) .and. nsize /= 0 ) then 
       if( nposi <= 0 ) then
          call runend('CHKTYP_R3PTYP: TYPE IS WRONG') 
       end if
       if( .not. associated(xarra(nposi) % a) ) then
          call runend('CHKTYP_R3PTYP: TYPE WAS NOT CREATED')                 
       end if
       if( present(ndim1) ) then 
          if( size( xarra(nposi) % a,1,KIND=ip) /= ndim1 ) then
             call runend('CHKTYP_R3PTYP: WRONG TYPE 1ST DIMENSION')             
          end if
       end if
       if( present(ndim2) ) then 
          if( size( xarra(nposi) % a,2,KIND=ip) /= ndim2 ) then
             call runend('CHKTYP_R3PTYP: WRONG TYPE 2ND DIMENSION')             
          end if
       end if
       if( present(ndim3) ) then 
          if( size( xarra(nposi) % a,3,KIND=ip) /= ndim3 ) then
             call runend('CHKTYP_R3PTYP: WRONG TYPE 2ND DIMENSION')             
          end if
          if( present(ndim2) ) then 
             if( size( xarra(nposi) % a,2,KIND=ip) /= ndim2 ) then
                call runend('CHKTYP_R3PTYP: WRONG TYPE 2ND DIMENSION')             
             end if
          end if
          if( present(ndim3) ) then 
             if( size( xarra(nposi) % a,3,KIND=ip) /= ndim3 ) then
                call runend('CHKTYP_R3PTYP: WRONG TYPE 3RD DIMENSION')             
             end if
          end if
       end if
    end if

  end subroutine chktyp_r3ptyp

end module mod_chktyp
!> @}
