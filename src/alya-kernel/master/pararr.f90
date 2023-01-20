!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine pararr(wtask,ntype,ndimr,rvarr)
  !------------------------------------------------------------------------
  !****f* Parall/pararr
  ! NAME
  !    pararr
  ! DESCRIPTION
  !    Works with arrays to deal with parallelization
  ! OUTPUT
  ! USED BY
  !    Parall
  !***
  !------------------------------------------------------------------------
  use def_kintyp,         only  :  ip,rp,lg
  use def_master,         only  :  npari,nparr,nparc,parre,parr1,kfl_paral
  use def_master,         only  :  party,parki,pardi,IPARALL,pard1
  use def_master,         only  :  ISEQUEN,INOTEMPTY
  use def_master,         only  :  NPOIN_TYPE,NFACE_TYPE,nfacg
  use def_domain,         only  :  npoin,nperi
  use def_master,         only  :  NEDGE_TYPE
  use def_domain,         only  :  nedge
  use mod_couplings,      only  :  COU_INTERPOLATE_NODAL_VALUES_go
  use mod_communications, only  :  PAR_INTERFACE_NODE_EXCHANGE
  use def_coupli,         only  :  ncoup_implicit
  use def_coupli,         only  :  ncoup_implicit_d
  use def_coupli,         only  :  lcoup_implicit_d
  use def_coupli,         only  :  ncoup_implicit_n
  use def_coupli,         only  :  lcoup_implicit_n
  use def_coupli,         only  :  RESIDUAL,UNKNOWN
  use def_coupli,         only  :  BETWEEN_SUBDOMAINS
  use def_master,         only  :  ISLAVE,mmodu
  use def_solver,         only  :  solve_sol
  use mod_communications, only  :  PAR_MIN,PAR_MAX,PAR_SUM
  use mod_periodicity,    only  :  periodicity_sequential
  use def_mpi
#include "def_mpi.inc" 
  implicit none 
  character(3), intent(in) :: wtask
  integer(ip),  intent(in) :: ntype,ndimr
  real(rp),     target     :: rvarr(ndimr)

  integer(ip)              :: kk,icoup
  integer(ip)              :: kcoup
  real(rp),     pointer    :: rvarr_tmp(:)
  real(rp),     pointer    :: rvarr_int(:)


  if( IPARALL .or. ( ISEQUEN .and. (wtask == 'SLX' .or. wtask == 'SOL' .or. wtask == 'SLA' ) ) ) then

     npari = 0
     nparc = 0

     select case ( wtask )

     case ( 'SND' )

        nparr =  ndimr
        parre => rvarr
        call par_sendin()

     case ( 'RCV' ) 

        nparr =  ndimr
        parre => rvarr
        call par_receiv()

     case ( 'SLX' , 'SOL' )

        if( ISEQUEN .and. nperi /= 0 .and. ntype == NPOIN_TYPE ) then
           pard1 = ndimr/npoin
           call periodicity_sequential(pard1,rvarr)
        end if
        !
        ! par_slexch
        !
        if( INOTEMPTY ) then

           if( ntype == NPOIN_TYPE ) then

              pard1 = ndimr/npoin
              if( ISLAVE ) call PAR_INTERFACE_NODE_EXCHANGE(pard1,rvarr,'SUM')
              
           else if( ntype == NFACE_TYPE ) then

              pard1 = ndimr/nfacg
              party = ntype
              if( pard1 == 1 ) then
                 parki =  2
                 pardi =  1
              else
                 parki =  5
                 pardi =  1
              end if

              if( ISLAVE ) then
                 parr1 => rvarr
                 call par_lgface(2_ip)
              end if

           else
              call runend('PARARR: NOT CODED')
           end if
           !
           ! Subdomain coupling 
           !
           if( ntype == NPOIN_TYPE .and. ncoup_implicit > 0 .and. INOTEMPTY ) then

              allocate( rvarr_tmp(ndimr),rvarr_int(ndimr) )
              rvarr_tmp(1:ndimr) = rvarr(1:ndimr)
              kk                 = ndimr/npoin
              !
              ! Neumann transmission condition
              !
              do kcoup = 1,ncoup_implicit_n
                 icoup = lcoup_implicit_n(kcoup)
                 rvarr_int(1:ndimr) = 0.0_rp
                 if( wtask == 'SOL' ) then
                    call COU_INTERPOLATE_NODAL_VALUES_go(icoup,kk,rvarr_int,rvarr_tmp,solve_sol(1) &
                         &% kfl_fixno)
                 else
                    call COU_INTERPOLATE_NODAL_VALUES_go(icoup,kk,rvarr_int,rvarr_tmp)
                 end if
                 rvarr(1:ndimr) = rvarr(1:ndimr) + rvarr_int(1:ndimr)
              end do
              deallocate(rvarr_int)
              rvarr_tmp(1:ndimr) = rvarr(1:ndimr)
              !
              ! Dirichlet transmission condition
              !              
              do kcoup = 1,ncoup_implicit_d
                 icoup = lcoup_implicit_d(kcoup)
                 if( wtask == 'SOL' ) then
                    call COU_INTERPOLATE_NODAL_VALUES_go(icoup,kk,rvarr,rvarr_tmp,solve_sol(1) %&
                         & kfl_fixno)
                 else
                    call COU_INTERPOLATE_NODAL_VALUES_go(icoup,kk,rvarr,rvarr_tmp)
                 end if
              end do
              deallocate(rvarr_tmp)
           end if
        end if

     case ( 'SLA' ) 
        
        if( ISEQUEN .and. nperi /= 0 .and. ntype == NPOIN_TYPE ) then
           pard1 = ndimr/npoin
           call periodicity_sequential(pard1,rvarr)
        end if
        
        if( ntype == NPOIN_TYPE ) then

           if( INOTEMPTY .and. ISLAVE ) then
              pard1 = ndimr/npoin
              call PAR_INTERFACE_NODE_EXCHANGE(pard1,rvarr,'SUM',wsynch='ASYNCHRONOUS')
           end if
           
        else
           call runend('PARARR: NOT CODED')
        end if
        
     case default

        call runend('PARARR: WRONG CASE')

     end select

     nparr = 0
     nullify(parre)
     nullify(parr1)

  end if

end subroutine pararr
