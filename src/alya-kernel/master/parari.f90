!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine parari(wtask,ntype,ndimi,rvari)
  !------------------------------------------------------------------------
  !****f* Parall/parari
  ! NAME
  !    parari
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
  use def_master,         only : NPOIN_TYPE,NFACE_TYPE,nfacg
  use def_domain,         only : npoin
  use mod_communications, only : PAR_MIN,PAR_MAX,PAR_SUM,PAR_INTERFACE_NODE_EXCHANGE
  implicit none
  character(3), intent(in) :: wtask
  integer(ip),  intent(in) :: ntype,ndimi
  integer(ip),  target     :: rvari(ndimi)
 
  if( IPARALL ) then

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

     case ( 'SLX' ) 
        !
        ! par_slexch for vectors(ndimi)
        !
        if( ntype == NPOIN_TYPE ) then
           
           pard1 = ndimi/npoin
           call PAR_INTERFACE_NODE_EXCHANGE(pard1,rvari,'SUM')

        else if( ntype == NFACE_TYPE ) then

           if( nfacg == 0 ) return
           pard1 = ndimi/nfacg
           party = ntype
           if( pard1 == 1 ) then
              parki =  1
              pardi =  1
           else
              parki =  6
              pardi =  2
           end if
           pari1 => rvari
           call par_lgface(2_ip)
        else
           call runend('PARARR: NOT CODED')
        end if

     case default
        
        call runend('PARARI: WRONG CASE') 
        
     end select
     
     npari = 0
     nullify(parin)
     nullify(pari1)

  end if
  
end subroutine parari
 
