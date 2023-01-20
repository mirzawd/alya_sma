!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_endite()
  !-----------------------------------------------------------------------
  !****f* Alefor/ale_endite
  ! NAME 
  !    ale_endite
  ! DESCRIPTION
  !    End of iterations   
  ! USES
  !    ale_coupli
  ! USED BY
  !    ale_doiter
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_alefor
  implicit none
  !
  ! Coupling
  !
  call ale_coupli(ITASK_ENDITE)
  !
  ! Update (:,2) <= (:,1)
  !
  call ale_updunk(ITASK_ENDITE)

end subroutine ale_endite

