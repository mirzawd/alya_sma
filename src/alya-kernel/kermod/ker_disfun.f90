!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kermod
!> @{
!> @file    ker_disfun.f90
!> @author  Guillaume Houzeaux
!> @date    29/10/2014
!> @brief   Define advection velocity
!> @details Define the advection ADVEC:
!>
!>          ADVEC(1:NDIME,1:NPOIN,1) ... current advection
!>          ADVEC(1:NDIME,1:NPOIN,2) ... last coupling advection
!>          ADVEC(1:NDIME,1:NPOIN,3) ... last time step advection
!>
!>          According to KFL_DIFUN, it is computed as:
!>
!>          kfl_difun = 0 ... ADVEC => VELOC. Nothing to do here.
!>                    < 0 ... ADVEC => FIELD. Constant in time and defined in ker_memall
!>                    > 0 ... ADVEC comes from a user defined function
!>
!>          Therefore, if kfl_difun, there is nothing to do
!> @} 
!-----------------------------------------------------------------------
subroutine ker_disfun(itask)
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use mod_memory
  use mod_chktyp,                  only : check_type
  use mod_ker_space_time_function, only : ker_space_time_function
  use mod_messages,                only : messages_live
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,idime,ifunc,ifiel,kk,icomp
  real(rp)                :: xx

  !----------------------------------------------------------------
  !
  ! DISPM Computed only for user defined functions
  !
  !----------------------------------------------------------------

  if( kfl_difun /= 0 ) then
     !
     ! Geometrical arrays should be updated!!!
     !
     kfl_domar = 1
     
     if( itask == ITASK_BEGSTE .or. itask == ITASK_INIUNK ) then        

        if( kfl_difun < 0 ) then
           !
           ! From Field
           !          
           ifiel = -kfl_difun
           if( kfl_field(4,ifiel) > 1 ) then
              kk = k_tran_fiel(ifiel)
              xx = x_tran_fiel(ifiel)
              do ipoin = 1,npoin
                 do idime = 1,ndime
                    dispm(idime,ipoin,1) = xfiel(ifiel) % a(idime,ipoin,kk) * xx * difun_facto + &
                         xfiel(ifiel) % a(idime,ipoin,kk+1) * (1.0_rp-xx) * difun_facto
                 end do
              end do
           else
              do ipoin = 1,npoin
                 do idime = 1,ndime
                    dispm(idime,ipoin,1) = xfiel(ifiel) % a(idime,ipoin,1) * difun_facto
                 end do
              end do
           end if

        else if( kfl_difun > 1000 ) then

           do ipoin = 1,npoin
              ifunc = kfl_difun  - 1000     
              call ker_space_time_function(&
                   ifunc,coord(1,ipoin),coord(2,ipoin),coord(ndime,ipoin),cutim,dispm(1:ndime,ipoin,1))
              do idime=1,ndime
                dispm(idime,ipoin,1) = dispm(idime,ipoin,1) * difun_facto
              enddo
           end do

        else

           call runend('KER_DISFUN: NOT CODED')

        end if
        !
        ! Assume constant initial dispmtion
        !
        if( itask == ITASK_INIUNK ) then
           do icomp = 2,memory_size(dispm,3_ip)
              do ipoin = 1,npoin 
                 do idime = 1,ndime
                    dispm(idime,ipoin,icomp) = dispm(idime,ipoin,1)
                 end do
              end do
           end do
        end if
        !
        ! VELOM: Compute mesh velocity
        !
        do ipoin = 1,npoin 
           do idime = 1,ndime
              velom(idime,ipoin) = (dispm(idime,ipoin,1)-dispm(idime,ipoin,3))*dtinv
           end do
        end do
        !
        ! COORD: move coordinates
        ! x^n+1 = x^0 + d^n+1
        ! x^n   = x^0 + d^n   =>
        ! x^n+1 = x^n + d^n+1 - d^n
        !
        do ipoin = 1,npoin 
           do idime = 1,ndime
              coord(idime,ipoin) = coord(idime,ipoin) + (dispm(idime,ipoin,1)-dispm(idime,ipoin,3))
           end do
        end do
        !
        ! Actualize geometry
        !
        call messages_live('MOVING MESH: RECOMPUTE MESH DEPENDENT ARRAYS','START SECTION')
        call domarr(2_ip)
        call messages_live('MOVING MESH','END SECTION')

     else if( itask == ITASK_ENDSTE ) then

        !----------------------------------------------------------------
        !
        ! Save previous dispmtion
        ! KFL_DIFUN = 0, DISPM point to VELOC which should not be modified
        !
        !----------------------------------------------------------------

        do ipoin = 1,npoin 
           do idime = 1,ndime
              dispm(idime,ipoin,3) = dispm(idime,ipoin,1)
           end do
        end do

     end if

  end if

end subroutine ker_disfun
