!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_number_mesh(itask)

  use def_kintyp
  use def_parame
  use def_domain
  use def_master
  use mod_memory
  use mod_communications, only : PAR_ALLGATHER
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications, only : PAR_SUM
  use mod_communications, only : PAR_MAX
  use mod_parall,         only : par_memor
  implicit none
  integer(ip), intent(in)  :: itask
  integer(ip)              :: ipoin,ii,kpoin
  integer(ip), pointer     :: lninv_tmp(:) 
  integer(ip), pointer     :: own_nodes(:) 
  integer(ip)              :: npoin_total_old

  nullify(lninv_tmp)
  nullify(own_nodes)

  select case ( itask )

  case ( 1_ip )
     !
     ! Copy old LNINV_LOC
     !
     call memory_alloca(par_memor,'LNINV_TMP','par_submsh',lninv_tmp,npoin_old) 
     do ipoin = 1,npoin_old
        lninv_tmp(ipoin) = lninv_loc(ipoin)
     end do
     call memory_deallo(par_memor,'LNINV_LOC','par_submsh',lninv_loc)
     call memory_alloca(par_memor,'LNINV_LOC','par_submsh',lninv_loc,npoin)
     do ipoin = 1,npoin_old
        lninv_loc(ipoin) = lninv_tmp(ipoin) 
     end do
     call memory_deallo(par_memor,'LNINV_TMP','par_submsh',lninv_tmp)

  case ( 2_ip )
     !
     ! Number missing nodes
     !
     if( INOTMASTER ) npoin_total_old = maxval( lninv_loc )
     call PAR_MAX(npoin_total_old)

     if( INOTMASTER ) then

        do ipoin = npoi1+1,npoi2-1
           lninv_loc(ipoin) = 0
        end do
        do ipoin = npoi3+1,npoin 
           lninv_loc(ipoin) = 0
        end do
        kpoin = 0
        do ipoin = 1,npoin
           if( lninv_loc(ipoin) == 0 ) then
              if( ipoin <= npoi1 .or. ( ipoin >= npoi2 .and. ipoin <= npoi3 ) ) then
                 kpoin = kpoin + 1
              end if
           end if
        end do
     else
        kpoin = 0
     end if

     call memory_alloca(par_memor,'ONW_NODES','par_submsh',own_nodes,npart+1_ip,'INITIALIZE',0_ip)
     call PAR_ALLGATHER(kpoin,own_nodes,1_4,'IN MY CODE')

     if( INOTMASTER ) then
        do ii = 1,npart
           own_nodes(ii) = own_nodes(ii-1) + own_nodes(ii)
        end do
        kpoin = own_nodes(kfl_paral-1) + npoin_total_old 

        do ipoin = 1,npoin
           if( lninv_loc(ipoin) == 0 ) then
              if( ipoin <= npoi1 .or. ( ipoin >= npoi2 .and. ipoin <= npoi3 ) ) then
                 kpoin = kpoin + 1
                 lninv_loc(ipoin) = kpoin
              end if
           end if
        end do

        call PAR_INTERFACE_NODE_EXCHANGE(lninv_loc,'SUM')

     end if

     !if( INOTMASTER ) npoin_total_old = maxval( lninv_loc )
     !call PAR_MAX(npoin_total_old)
     !if( IMASTER ) print*,'MAX=',npoin_total_old

     call memory_deallo(par_memor,'ONW_NODES','par_submsh',own_nodes)

  end select

end subroutine par_number_mesh
