!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_srinte(itask,ndimi,rvari)
  !------------------------------------------------------------------------
  !****f* Parall/par_srinte
  ! NAME
  !    par_srinte
  ! DESCRIPTION
  ! OUTPUT
  ! USED BY
  !    Parall
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only    :  ip,rp
  use def_master, only    :  npari,nparr,nparc,parin
  implicit none
  integer(ip), intent(in) :: itask,ndimi
  integer(ip), target     :: rvari(max(1_ip,ndimi))

  nparr = 0
  nparc = 0

  select case ( itask )

  case ( 1_ip ) 

     npari =  ndimi
     parin => rvari
     call par_sendin()
     nullify(parin)

  case ( 2_ip )

     npari =  ndimi
     parin => rvari
     call par_receiv()
     nullify(parin)

  end select
  
end subroutine par_srinte
 
