!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_parari(wtask,ntype,ndimi,rvari)
  !------------------------------------------------------------------------
  !****f* Parall/par_parari
  ! NAME
  !    par_parari
  ! DESCRIPTION
  !    Works with arrays to deal with parallelization
  ! OUTPUT
  ! USED BY
  !    Parall
  !***
  !------------------------------------------------------------------------
  use def_kintyp,         only : ip,rp
  use def_master,         only : npari,npari,nparc,nparr,parin,parin,pari1
  use def_master,         only : party,parki,pardi,IPARALL,pard1
  use def_master,         only : NPOIN_TYPE,NBOPO_TYPE
  use def_domain,         only : npoin,nbopo
  use mod_communications, only : PAR_MIN,PAR_MAX,PAR_SUM
  implicit none
  character(3), intent(in) :: wtask
  integer(ip),  intent(in) :: ntype,ndimi
  integer(ip),  target     :: rvari(ndimi)
 
  if( IPARALL ) then

     npari = 0
     nparr = 0
     nparc = 0

     select case ( wtask )

     case ( 'SND' ) 
        !
        ! par_sendin
        !
        npari =  ndimi
        parin => rvari 
        call par_sendin() 

     case ( 'RCV' )
        !
        ! par_receiv
        !
        npari =  ndimi
        parin => rvari
        call par_receiv() 

     case ( 'GAT' )
        !
        ! par_gather
        !
        if( ntype == NPOIN_TYPE .or. ntype == NBOPO_TYPE ) then
           if( ntype == NPOIN_TYPE ) pard1 = ndimi/npoin
           if( ntype == NBOPO_TYPE ) pard1 = ndimi/nbopo
           party = ntype
           if( pard1 == 1 ) then
              parki =  1
              pardi =  1
           else
              parki =  5
              pardi =  1
           end if
        else
           party =  ntype
           parki =  1
           pardi =  1
           npari =  ndimi
        end if
        pari1 => rvari
        call par_mygather() 

     case default
        
        call runend('PAR_PARARI: WRONG CASE')
        
     end select
     
     npari = 0
     nullify(parin)
     nullify(pari1)

  end if
  
end subroutine par_parari
 
