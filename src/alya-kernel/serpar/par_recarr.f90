!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_recarr(itask,ndimi,ndimr,ivari,rvari)
  !------------------------------------------------------------------------
  !****f* Parall/par_recarr
  ! NAME
  !    par_recarr
  ! DESCRIPTION
  ! OUTPUT
  ! USED BY
  !    Parall
  !***
  !------------------------------------------------------------------------
  use def_master
  implicit none
  integer(ip), intent(in) :: itask,ndimi,ndimr
  integer(ip), target     :: ivari(ndimi)
  real(rp),    target     :: rvari(ndimr)

  npari = 0
  nparr = 0
  nparc = 0

  select case ( itask )

  case ( 1_ip ) 

     npari =  ndimi
     parin => ivari 
     call par_receiv()

  case ( 2_ip )

     nparr =  ndimr
     parre => rvari
     call par_receiv()

  end select
  
end subroutine par_recarr
 
