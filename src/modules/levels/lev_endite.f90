!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_endite(itask)
  !-----------------------------------------------------------------------
  !****f* Levels/lev_endite
  ! NAME 
  !    lev_endite
  ! DESCRIPTION
  !    This routine checks convergence and performs updates of the
  !    temperature  at:
  !    - itask=1 The end of an internal iteration
  !    - itask=2 The end of the internal loop iteration
  ! USES
  !    lev_cvgunk
  !    lev_updunk
  ! USED BY
  !    lev_doiter
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_levels
  use def_kermod,     only : kfl_cutel
  use mod_ker_elmcut, only : ker_elmcut_free_surface
  use mod_messages, only : livinf
  implicit none
  integer(ip) :: itask

  select case(itask)

  case(1)
     !
     !  Compute convergence residual of the internal iteration (that is,
     !  || u(n,i,j) - u(n,i,j-1)|| / ||u(n,i,j)||) and update unknowns:
     !  u(n,i,j-1) <-- u(n,i,j) 
     !
     call lev_cvgunk(1_ip)
     call lev_updunk(3_ip)
     !   call lev_output(2_ip)

  case(2)
     !
     !  Compute convergence residual of the external iteration (that is,
     !  || u(n,i,*) - u(n,i-1,*)|| / ||u(n,i,*)||) and update unknowns:
     !  u(n,i-1,*) <-- u(n,i,*) 
     !
     call livinf(16_ip,' ',itinn(modul))
     call lev_cvgunk(2_ip)
     call lev_updunk(4_ip)

     if( kfl_cutel == 1 ) call ker_elmcut_free_surface()

  end select

end subroutine lev_endite
