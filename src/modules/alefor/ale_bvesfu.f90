!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_bvesfu()
  use def_kintyp, only : ip,rp
  use mod_maths
  use def_domain
  use def_master
  use def_alefor
  use mod_ker_space_time_function,    only : ker_space_time_function

  implicit none

  integer(ip)                  :: ipoin,ibopo,iroty,itotv,ifunc
  real(rp)                     :: vefun(3),worma(3),worve(3)


  do ipoin = 1,npoin
     ibopo = lpoty(ipoin)
     if( kfl_funno_ale(ipoin) < 0 ) then 
        ifunc = -kfl_funno_ale(ipoin)      
        !print*,'ifunc',ifunc      
        call ker_space_time_function(&
             ifunc,coord(1,ipoin),coord(2,ipoin),coord(ndime,ipoin),cutim,vefun(1))
        call ker_space_time_function(&
             ifunc,coord(1,ipoin),coord(2,ipoin),coord(ndime,ipoin),cutim-dtime,vefun(2))
        !print*,'vefun1',vefun(1)
        !print*,'vefun2',vefun(2)
        bvess_ale(1:ndime,ipoin,1) = bvess_ale(1:ndime,ipoin,1)*(vefun(1) - vefun(2))
        !print*,'bvess_ale-1',bvess_ale(1:ndime,ipoin)

        if( ibopo > 0 ) then
           iroty = kfl_fixrs_ale(ibopo)
           !print*,'iroty',iroty
           if( iroty == -1 ) then
              !
              ! Boundary conditions in the tangent skew system
              !
              itotv=(ipoin-1)*ndime
              worve(1:ndime)=bvess_ale(1:ndime,ipoin,1)
              call mbvab0(worma,exnor(1,1,ibopo),worve,ndime,ndime)
              bvess_ale(1:ndime,ipoin,1)=worma(1:ndime)
           end if
        end if
     end if
  end do

end subroutine ale_bvesfu
