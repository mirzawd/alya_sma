!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine exm_endite(itask)
  !-----------------------------------------------------------------------
  !****f* Exmedi/exm_endite
  ! NAME 
  !    exm_endite
  ! DESCRIPTION
  !    This routine checks convergence and updates unknowns at:
  !    - itask=1 The end of an internal iteration
  !    - itask=2 The end of the internal loop iteration
  ! USES
  !    exm_cvgunk
  ! USED BY
  !    exm_doiter
  !***
  !-----------------------------------------------------------------------
  use      def_master
  use      def_parame
  use      def_exmedi
  use      mod_messages,            only : livinf


  implicit none
  integer(ip) :: itask
  character(300)   :: messa
  integer(ip)      :: maxiter

  select case(itask)

  case(ITASK_ENDINN)
     !
     !  Compute convergence residual of the internal iteration (that is,
     !  || u(n,i,j) - u(n,i,j-1)|| / ||u(n,i,j)||) and update unknowns:
     !  u(n,i,j-1) <-- u(n,i,j)
     !
     call exm_cvgunk(ITASK_ENDINN)
     call exm_updunk(ITASK_ENDINN)     !  u(,ITER_K) <-- unkno  

     maxiter= solve_sol(1) % miter
     messa = &
          ' (SUBIT: '//trim(intost(itinn(modul)))//'/'//trim(intost(miinn_exm))//' IT: '//trim(intost(last_iters_exm))//'/'//trim(intost(maxiter))//')'
     call livinf(-3_ip,messa,one)
     call livinf(56_ip,' ',modul)


  case(ITASK_ENDITE)
     !
     !  Compute convergence residual of the external iteration (that is,
     !  || u(n,i,*) - u(n,i-1,*)|| / ||u(n,i,*)||) and update unknowns:
     !  u(n,i-1,*) <-- u(n,i,*)
     !
     call livinf(16_ip,' ',itinn(modul))
     call exm_cvgunk(ITASK_ENDITE)
     call exm_updunk(ITASK_ENDITE)     !  u(,ITER_AUX) <-- u(,ITER_K)

  end select

end subroutine exm_endite
