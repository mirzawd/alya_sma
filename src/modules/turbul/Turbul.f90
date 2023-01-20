!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine Turbul(order)
  !------------------------------------------------------------------------
  !****f* Turbul/turbul
  ! NAME
  !   Turbul
  ! DESCRIPTION
  !   This routine deals with the turbulence equation. The task done
  !   corresponds to the order given by the master.
  ! USES
  ! USES
  !    tur_turnon
  !    tur_timste
  !    tur_begste
  !    tur_doiter
  !    tur_concon
  !    tur_conblk
  !    tur_newmsh
  !    tur_endste
  !    tur_turnof
  ! USED BY
  !    Reapro
  !    Turnon
  !    Timste
  !    Begste
  !    Doiter
  !    Concon
  !    Conblk
  !    Newmsh
  !    Endste
  !    Turnof
  !***
  !------------------------------------------------------------------------
  use      def_master
  implicit none
  integer(ip), intent(in) :: order

  select case (order)

  case( ITASK_TURNON )        ; call tur_turnon()
  case( ITASK_SOLMEM )        ; call tur_solmem()
  case( ITASK_TIMSTE )        ; call tur_timste()
  case( ITASK_INIUNK )        ; call tur_iniunk()
  case( ITASK_BEGSTE )        ; call tur_begste()
  case( ITASK_DOITER )        ; call tur_doiter()
  case( ITASK_CONCOU )        ; call tur_concou()
  case( ITASK_CONBLK )        ; call tur_conblk()
  case( ITASK_ENDSTE )        ; call tur_endste()
  case( ITASK_OUTPUT )        ; call tur_output()
  case( ITASK_TURNOF )        ; call tur_turnof()
  case( ITASK_READ_RESTART  ) ; call tur_restar(ITASK_READ_RESTART)
  case( ITASK_WRITE_RESTART ) ; call tur_restar(ITASK_WRITE_RESTART)
     
  end select

end subroutine Turbul
!------------------------------------------------------------------------
!
!  tur_reapro -> Definition of the model (SA, k-e, etc.)
!  tur_turnon -> tur_reaphy -> Definition of constants (Cmu, Ce1, etc.)
!                ----------
!             -> tur_addarr -> Distance to the wall, nearest wall node
!                ----------
!  tur_begste -> tur_iniunk -> Initial values for turbulence variables
!  tur_doiter -> tur_begite -> tur_frivel -> Calculate U*=ustar_tur
!                              ----------
!                tur_solite -> tur_updibc -> Update b.c. 
!                              ----------
!                              tur_elmope -> Ensemble matrix A and b
!                                         -> tur_elmcoe -> Model coef.
!                                            ----------   
!                              Solve A*unkno=b
!                tur_endite -> tur_updunk -> untur=relax*unkno+(1-relax)*untur
!                              ker_proper -> Update mut=turmu
!                              ----------
!  tur_concou
!  tur_conblk
!  tur_newmsh
!  tur_endste
!  tur_turnof
!
!------------------------------------------------------------------------

