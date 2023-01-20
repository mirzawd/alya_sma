!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_parasr(wtask,ndims,rvars,ndimr,rvarr)
  !------------------------------------------------------------------------
  !****f* Parall/par_parasr
  ! NAME
  !    par_parasr
  ! DESCRIPTION
  !    Works with arrays to deal with parallelization
  ! OUTPUT
  ! USED BY
  !    Parall
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_master, only     :  npasr,nparr,parrs,parre
  use def_master, only     :  ISLAVE
  implicit none
  character(3), intent(in) :: wtask
  integer(ip),  intent(in) :: ndims,ndimr
  real(rp),     target     :: rvars(ndims)
  real(rp),     target     :: rvarr(ndimr)

  if( ISLAVE ) then

     npasr = 0
     nparr = 0 

     select case ( wtask )

     case ( 'SLX' ) 
        !
        ! par_slaves
        !
        npasr =  ndims
        nparr =  ndimr
        parrs => rvars
        parre => rvarr
        call par_slaves(1_ip)
        
     case default
        
        call runend('PARASI: WRONG CASE')
        
     end select

     npasr = 0
     nparr = 0

  end if
  
end subroutine par_parasr
 
