!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Algebraic_Solver
!> @{
!> @file    outdef.f90
!> @author  Guillaume Houzeaux
!> @brief   Groups information
!> @details Output groups information
!> @} 
!-----------------------------------------------------------------------
subroutine outdef(itask)
  use def_kintyp
  use def_master
  use def_domain
  use def_solver
  use mod_communications, only : PAR_SUM
  use mod_iofile,         only : iofile_flush_unit
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ngrou,igrou,ipoin
  integer(ip), pointer    :: lgrou(:)

  if( itask == 1 ) then
     !
     ! Write number of nodes per group
     !
     ngrou = solve_sol(1) % ngrou
     allocate(lgrou(0:ngrou))
     do igrou = 0,ngrou
        lgrou(igrou) = 0
     end do
 
     if( INOTMASTER ) then
        do ipoin = 1,npoi1
           igrou = solve_sol(1) % lgrou(ipoin)
           lgrou(igrou) = lgrou(igrou) + 1
        end do
        do ipoin = npoi2,npoi3
           igrou = solve_sol(1) % lgrou(ipoin)
           lgrou(igrou) = lgrou(igrou) + 1
        end do
     end if

     call PAR_SUM(ngrou+1_ip,lgrou(0:ngrou))

     if( INOTSLAVE ) then
        do igrou = 0,ngrou
           write(solve_sol(1) % lun_solve,1) igrou,lgrou(igrou)
        end do
        call iofile_flush_unit(solve_sol(1) % lun_solve)
     end if

     deallocate(lgrou)

  end if

1 format('# GROUPS ',i7,i9)

end subroutine outdef
