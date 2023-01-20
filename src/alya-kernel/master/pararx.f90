!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine pararx(wtask,ntype,ndimx,rvarx)
  !------------------------------------------------------------------------
  !****f* Parall/pararx
  ! NAME
  !    pararr
  ! DESCRIPTION
  !    Works with arrays to deal with parallelization
  ! OUTPUT
  ! USED BY
  !    Parall
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_master, only     :  npari,nparr,nparc
  use def_master, only     :  IPARALL
!  use def_master, only     :  IPARALL,pard1
  use def_master, only     :  NPOIN_TYPE,NBOPO_TYPE
  use mod_communications, only  :  PAR_INTERFACE_NODE_EXCHANGE
  implicit none
  character(3), intent(in) :: wtask
  integer(ip),  intent(in) :: ntype,ndimx
  complex(rp),  target     :: rvarx(ndimx)

  if( IPARALL ) then

     npari = 0
     nparc = 0
     nparr = 0

     select case ( wtask )

     case ( 'SLX' )
        !
        ! par_slexch
        !
        if( ntype == NPOIN_TYPE .or. ntype == NBOPO_TYPE ) then
           call runend('EXHCNAGE FOR COMPLEX SHOULD BE CODED')
           !call PAR_INTERFACE_NODE_EXCHANGE(pard1,rvarx)
        else
           call runend('PARARX: NOT CODED')
        end if

     case default

        call runend('PARARX: WRONG CASE')

     end select

     !nparr = 0

  end if

end subroutine pararx
 
