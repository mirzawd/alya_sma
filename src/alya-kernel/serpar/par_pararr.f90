!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_pararr(wtask,ntype,ndimr,rvarr)
  !------------------------------------------------------------------------
  !****f* Parall/par_pararr
  ! NAME
  !    par_pararr
  ! DESCRIPTION
  !    Works with arrays to deal with parallelization
  ! OUTPUT
  ! USED BY
  !    Parall
  !***
  !------------------------------------------------------------------------
  use def_kintyp,         only : ip,rp
  use def_master,         only : npari,nparr,nparc,parre,parr1
  use def_master,         only : party,parki,pardi,IPARALL,pard1,pari1
  use def_master,         only : NPOIN_TYPE,NBOPO_TYPE,NBOUN_TYPE,NELEM_TYPE,lntra
  use def_domain,         only : npoin,nbopo,nboun,nelem
  use mod_communications, only : PAR_MIN,PAR_MAX,PAR_SUM
  implicit none
  character(3), intent(in) :: wtask
  integer(ip),  intent(in) :: ntype,ndimr
  real(rp),     target     :: rvarr(ndimr)

  if( IPARALL ) then

     npari = 0
     nparc = 0
     nparr = 0

     select case ( wtask )

     case ( 'SND' ) 

        nparr =  ndimr
        parre => rvarr
        call par_sendin() 

     case ( 'RCV' )

        nparr =  ndimr
        parre => rvarr
        call par_receiv() 

     case ( 'GAT' ) 

        if( ntype == NPOIN_TYPE .or. ntype == NBOPO_TYPE .or. ntype == NBOUN_TYPE .or. ntype == NELEM_TYPE ) then
           if( ntype == NPOIN_TYPE ) pard1 = ndimr/npoin
           if( ntype == NELEM_TYPE ) pard1 = ndimr/nelem
           if( ntype == NBOPO_TYPE ) pard1 = ndimr/max(1_ip,nbopo)
           if( ntype == NBOUN_TYPE ) pard1 = ndimr/max(1_ip,nboun)
           party = ntype
           if( pard1 == 1 ) then
              parki =  2
              pardi =  1
           else
              parki =  6
              pardi =  1
           end if
        else
           party =  ntype
           parki =  2
           pardi =  1
           npari =  ndimr
        end if
        parr1 => rvarr
        call par_mygather()

     case ( 'IBI' ) 

        pari1 => lntra
        pard1 =  ndimr/npoin
        parr1 => rvarr(1:ndimr)
        call par_slexib(2_ip)

     case default

        call runend('PARARR: WRONG CASE')

     end select

     nparr = 0
     nullify(parre)

  end if

end subroutine par_pararr
 
