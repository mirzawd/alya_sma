!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kermod
!> @{
!> @file    ker_arefun.f90
!> @author  Guillaume Houzeaux
!> @date    29/10/2014
!> @brief   Define advection velocity
!> @details Define the advection ADVEC:
!>
!>          ADVEC(1:NDIME,1:NPOIN,1) ... current advection
!>          ADVEC(1:NDIME,1:NPOIN,2) ... last coupling advection
!>          ADVEC(1:NDIME,1:NPOIN,3) ... last time step advection
!>
!>          According to KFL_VEFUN, it is computed as:
!>
!>          kfl_vefun = 0 ... ADVEC => VELOC. Nothing to do here.
!>                    < 0 ... ADVEC => FIELD. Constant in time and defined in ker_memall
!>                    > 0 ... ADVEC comes from a user defined function
!>
!>          Therefore, if kfl_vefun, there is nothing to do
!> @} 
!-----------------------------------------------------------------------
subroutine ker_arefun(itask)
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use mod_memory
  use mod_chktyp,                  only : check_type
  use mod_ker_space_time_function, only : ker_space_time_function

  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,ifunc,ifiel
  integer(ip)             :: kk,icomp
  real(rp)                :: xx

  !----------------------------------------------------------------
  !
  ! AREAS Computed only for user defined functions
  !
  !----------------------------------------------------------------

  if( INOTMASTER .and. kfl_arfun /= 0 ) then

     if( itask == ITASK_BEGSTE .or. itask == ITASK_INIUNK .or. itask == ITASK_BEGRUN ) then
        
        if( kfl_arfun < 0 ) then
           !
           ! From Field
           !
           ifiel = -kfl_arfun
           if( kfl_field(4,ifiel) > 1 ) then
              kk = k_tran_fiel(ifiel)
              xx = x_tran_fiel(ifiel)

              do ipoin = 1,npoin
                 areas(ipoin,1) = xfiel(ifiel) % a(1,ipoin,kk) * xx + &
                      xfiel(ifiel) % a(1,ipoin,kk+1) * (1.0_rp-xx)
              end do
           else
              do ipoin = 1,npoin
                 areas(ipoin,1) = xfiel(ifiel) % a(1,ipoin,1) 
              end do              
           end if
           
        else if( kfl_arfun > 1000 ) then
           !
           ! Space time function
           !
           do ipoin = 1,npoin
              ifunc = kfl_arfun  - 1000     
              call ker_space_time_function(&
                   ifunc,coord(1,ipoin),coord(2,ipoin),coord(ndime,ipoin),cutim,areas(ipoin,1))
           end do 

        else if( kfl_arfun > 0 ) then
           !
           ! Programmer functions
           !
           select case ( kfl_arfun )

           case ( 99_ip ) 
              !
              ! Do not do anything
              !
              continue

           end select

        end if
        !
        ! Assume constant initial areastion
        !
        if( itask == ITASK_INIUNK .or. itask == ITASK_BEGRUN ) then
           do icomp = 2,memory_size(areas,2_ip)
              do ipoin = 1,npoin 
                 areas(ipoin,icomp) = areas(ipoin,1)
              end do
           end do
        end if

     else if( itask == ITASK_ENDSTE ) then

        !----------------------------------------------------------------
        !
        ! Save previous areastion
        ! KFL_ARFUN = 0, AREAS point to VELOC which should not be modified
        !
        !----------------------------------------------------------------

        do ipoin = 1,npoin 
           areas(ipoin,3) = areas(ipoin,1)
        end do

     end if

  end if

end subroutine ker_arefun
