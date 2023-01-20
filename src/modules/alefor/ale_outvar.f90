!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Alefor
!> @{
!> @file    ale_outvar.f90
!> @author  Guillaume Houzeaux
!> @date    16/11/1966
!> @brief   Output a postprocess variable
!> @details Output a postprocess variable
!> @} 
!-----------------------------------------------------------------------
subroutine ale_outvar(ivari,imesh)
  use def_parame
  use def_master
  use def_domain
  use def_alefor
  use def_kermod,         only : kfl_adj_prob
  use mod_outvar,         only : outvar
  use mod_arrays,         only : arrays_name
  implicit none
  integer(ip), intent(in) :: ivari
  integer(ip), intent(in) :: imesh
  integer(ip)             :: idime,ipoin
  real(rp)                :: rutim
  !
  ! Define postprocess variable
  !
  rutim = cutim

  select case ( arrays_name(ivari) )  

  case( 'DISPM' )
     !
     ! DISPM: Mesh displacement
     !
     if( INOTMASTER ) then
        if( coupling('ALEFOR','SOLIDZ') >= 1 ) then
           call memgen(zero,ndime,npoin)
           do ipoin = 1,npoin
              do idime = 1,ndime
                 gevec(idime,ipoin) = displ(idime,ipoin,1) + coord(idime,ipoin) - coord_ale(idime,ipoin,2)
              end do
           end do
        else
!            gevec => dispm(:,:,1)
           !
           ! Alefor postprocess the total displacement with respect to the original configuration, like solidz
           !
           call memgen(zero,ndime,npoin)

           do ipoin = 1,npoin

              do idime = 1,ndime

                 if (kfl_adj_prob /= 1) then 
                   gevec(idime,ipoin) = coord(idime,ipoin) - coord_ori(idime,ipoin)
                 else
                   gevec(idime,ipoin) = dispm(idime,ipoin,1)
                 endif

              end do

           end do

        end if

     end if

  case( 'VELOM' )
     !
     ! Mesh velocity
     !
     if( INOTMASTER ) gevec => velom(:,:) 

  case( 'COALE' )
     !
     ! New coordinates
     !
     if( INOTMASTER ) gevec => coord_ale(:,:,1) 

  case( 'GROUP' )
     !
     ! Groups
     !
     if( INOTMASTER ) then
       do ipoin = 1,npoin
           rhsid(ipoin) = real(solve(1) % lgrou(ipoin),rp)
        end do
        gesca => rhsid
     end if

  case( 'BVESS' )
     !
     ! New coordinates
     !
     if( INOTMASTER ) gevec => bvess_ale(:,:,1)

  end select

  call outvar(&
       ivari,&
       ittim,rutim,postp(1) % wopos(:,ivari),MESH_ID=imesh)

end subroutine ale_outvar
