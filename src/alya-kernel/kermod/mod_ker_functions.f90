!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kermod
!> @{
!> @file    mod_ker_functions.f90
!> @author  houzeaux
!> @date    2020-03-30
!> @brief   Apply functions
!> @details Apply different types of functions
!-----------------------------------------------------------------------

module mod_ker_functions

  use def_kintyp,                  only : ip,rp,lg
  use def_master,                  only : FUNCTION_SPACE_TIME
  use def_master,                  only : FUNCTION_FIELD
  use def_master,                  only : FUNCTION_TIME
  use def_master,                  only : FUNCTION_DISCRETE
  use def_master,                  only : cutim
  use def_domain,                  only : coord
  use def_domain,                  only : xfiel
  use def_domain,                  only : ndime
  use def_domain,                  only : kfl_field
  use def_domain,                  only : k_tran_fiel
  use def_domain,                  only : x_tran_fiel
  use def_domain,                  only : lnnob
  use def_domain,                  only : lnodb
  use def_kermod,                  only : time_function
  use mod_ker_space_time_function, only : ker_space_time_function
  use mod_ker_discrete_function,   only : ker_discrete_function
  use mod_optional_argument,       only : optional_argument
  
  implicit none

  real(rp),  external   :: funcre
  
  interface ker_functions
     module procedure &
          &           ker_functions_scalar,&
          &           ker_functions_vector
  end interface ker_functions
          
contains

  subroutine ker_functions_vector(ienti,ifunc,itype,bvess,xx,TIME,COORDINATES,WHEREIN)
     
    integer(ip),                intent(in)  :: ienti
    integer(ip),                intent(in)  :: ifunc
    integer(ip),                intent(in)  :: itype 
    real(rp),                   intent(in)  :: bvess(:)
    real(rp),                   intent(out) :: xx(:)
    real(rp),         optional, intent(in)  :: TIME
    real(rp),         optional, intent(in)  :: COORDINATES(:)
    character(LEN=*), optional, intent(in)  :: WHEREIN
    integer(ip)                             :: idime,kk,ipoin,inodb
    real(rp)                                :: tt,rtime
    real(rp)                                :: xcoor(3)
    logical(lg)                             :: on_nodes

    on_nodes = .true.
    if( present(WHEREIN) ) then
       if( trim(WHEREIN) == 'ON BOUNDARIES' ) on_nodes = .false.
    end if
    
    if( present(COORDINATES) ) then
       !
       ! Coordinates gien by user 
       !
       xcoor(1:ndime) = COORDINATES(1:ndime)
       
    else if( on_nodes ) then
       !
       ! This is a node
       !
       xcoor(1:ndime) = coord(1:ndime,ienti)
       
    else
       !
       ! This is a boundary
       !
       xcoor = 0.0_rp
       do inodb = 1,lnnob(ienti)
          ipoin = lnodb(inodb,ienti)
          xcoor(1:ndime) = xcoor(1:ndime) + coord(1:ndime,ipoin)
       end do
       xcoor(1:ndime) = xcoor(1:ndime) / real(lnnob(ienti),rp)
    end if
    
    rtime = optional_argument(cutim,TIME)
    
    if( itype == 0 .or. ifunc == 0 ) then
       !
       ! Nothing to do
       !
       xx(1:ndime) = bvess(1:ndime)
       
    else if( itype == FUNCTION_SPACE_TIME ) then
       !
       ! Space time function
       !
       call ker_space_time_function(&
            ifunc,xcoor(1),xcoor(2),xcoor(ndime),rtime,xx(1:ndime))
       do idime = 1,ndime
          xx(idime) = xx(idime) * bvess(idime)
       end do

    else if( itype == FUNCTION_FIELD ) then

       if( kfl_field(4,ifunc) == 1 ) then
          !
          ! Constant field
          !
          do idime = 1,ndime
             xx(idime) = xfiel(ifunc) % a(idime,ienti,1) 
          end do
       else
          !
          ! Transient fields
          !
          kk = k_tran_fiel(ifunc)
          tt = x_tran_fiel(ifunc)
          do idime = 1,ndime
             xx(idime) = xfiel(ifunc) % a(idime,ienti,kk)   * tt + &
                  &         xfiel(ifunc) % a(idime,ienti,kk+1) * (1.0_rp-tt)
          end do
       end if

    else if( itype == FUNCTION_TIME ) then
       !
       ! Function: U(t)=f(t)*U(0)
       !
       do idime = 1,ndime
          xx(idime) = bvess(idime)  * funcre( &
               time_function(ifunc) % parameters,    &
               time_function(ifunc) % npara,         &
               time_function(ifunc) % kfl_type,      &
               rtime)
       end do

    else if( itype == FUNCTION_DISCRETE ) then
       !
       ! Discrete function
       !
       call ker_discrete_function(ifunc,cutim,xx(1:ndime))
       do idime = 1,ndime
          xx(idime) = xx(idime) * bvess(idime)
       end do
       
    end if

  end subroutine ker_functions_vector
  
  subroutine ker_functions_scalar(ienti,ifunc,itype,bvess,xx,TIME,COORDINATES,WHEREIN)

    integer(ip),                intent(in)  :: ienti
    integer(ip),                intent(in)  :: ifunc
    integer(ip),                intent(in)  :: itype 
    real(rp),                   intent(in)  :: bvess
    real(rp),                   intent(out) :: xx
    real(rp),         optional, intent(in)  :: TIME
    real(rp),         optional, intent(in)  :: COORDINATES(:)
    character(LEN=*), optional, intent(in)  :: WHEREIN
    integer(ip)                             :: kk,ipoin,inodb
    real(rp)                                :: tt,rtime
    real(rp)                                :: xcoor(3)
    logical(lg)                             :: on_nodes

    on_nodes = .true.
    if( present(WHEREIN) ) then
       if( trim(WHEREIN) == 'ON BOUNDARIES' ) on_nodes = .false.
    end if

    if( present(COORDINATES) ) then
       !
       ! Coordinates gien by user 
       !
       xcoor(1:ndime) = COORDINATES(1:ndime)
       
    else if( on_nodes ) then
       !
       ! This is a node
       !
       xcoor(1:ndime) = coord(1:ndime,ienti)
    else
       !
       ! This is a boundary
       !
       xcoor = 0.0_rp
       do inodb = 1,lnnob(ienti)
          ipoin = lnodb(inodb,ienti)
          xcoor(1:ndime) = xcoor(1:ndime) + coord(1:ndime,ipoin)
       end do
       xcoor(1:ndime) = xcoor(1:ndime) / real(lnnob(ienti),rp)
       
    end if

    rtime = optional_argument(cutim,TIME)
    
    if( itype == 0 .or. ifunc == 0 ) then
       !
       ! Nothing to do
       !
       return
       
    else if( itype == FUNCTION_SPACE_TIME ) then
       !
       ! Space time function
       !
       call ker_space_time_function(&
            ifunc,xcoor(1),xcoor(2),xcoor(ndime),rtime,xx)
       xx = xx * bvess

    else if( itype == FUNCTION_FIELD ) then

       if( kfl_field(4,ifunc) == 1 ) then
          !
          ! Constant field
          !
          xx = xfiel(ifunc) % a(1,ienti,1) 

       else
          !
          ! Transient fields
          !
          kk = k_tran_fiel(ifunc)
          tt = x_tran_fiel(ifunc)
          xx = xfiel(ifunc) % a(1,ienti,kk)   * tt + &
               &         xfiel(ifunc) % a(1,ienti,kk+1) * (1.0_rp-tt)

       end if

    else if( itype == FUNCTION_TIME ) then
       !
       ! Function: U(t)=f(t)*U(0)
       !
       xx = bvess * funcre( &
            time_function(ifunc) % parameters,    &
            time_function(ifunc) % npara,         &
            time_function(ifunc) % kfl_type,      &
            rtime)

    else if( itype == FUNCTION_DISCRETE ) then
       !
       ! Discrete function
       !
       call ker_discrete_function(ifunc,cutim,xx)
       xx = xx * bvess

    end if

  end subroutine ker_functions_scalar
  
end module mod_ker_functions
!> @}
