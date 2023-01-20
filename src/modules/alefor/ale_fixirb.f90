!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_fixirb
  !-----------------------------------------------------------------------
  !****f* alefor/ale_fixirb
  ! NAME
  !    ale_fixirb
  ! DESCRIPTION
  !    Obtain bvess_ale for rigid body
  ! USED BY
  !    domain
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_alefor

  implicit none
  integer(ip)             :: ipoin,idime,kpoin,iimbo


  if( INOTMASTER ) then
     if( associated(kfl_fixno_ale) ) then
        if ( kfl_mvext_ale == 1 ) then ! set x & y components of displacement = disp COG
           do ipoin = 1,npoin
              do idime = 1,2
                 if (kfl_fixno_ale(idime,ipoin) == 1) bvess_ale(idime,ipoin,1) = rbbou(1) % posil(idime,1) & 
                      - rbbou(1) % posil(idime,3)    !delta x_cog                
              end do
           end do
        end if
        do iimbo = 1,nrbod
           do kpoin = 1,rbbou(iimbo) % npoib
              ipoin = rbbou(iimbo) % lninv(kpoin)
              do idime = 1,ndime    
                 bvess_ale(idime,ipoin,1)     = rbbou(iimbo) % cooib(idime,kpoin) - coord(idime,ipoin)
                 kfl_fixno_ale(idime,ipoin) = 1
              end do
           end do
        end do
     end if 
  end if
  
end subroutine ale_fixirb
