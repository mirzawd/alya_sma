!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-------------------------------------------------------------------------||---!
!                                                                              !
!------------------------------------------------------------------------------!
subroutine allocate_block_system( mod_module ) !, module_k, react, bvnat) 
  !
  ! CREATE: 2015JAN29 
  ! 
  use def_kintyp, only: ip, rp 
  use def_kintyp_solvers, only: soltyp
  use def_kintyp, only: tymod 
  use def_domain, only: npoin, r_dom, npoin
  use mod_memory, only: memory_alloca  
  implicit none
  !integer(ip),  intent(in)    :: module_k
  !logical(ip),  intent(in)    :: react, bvnat 
  type(tymod),  intent(inout) :: mod_module 
  type(soltyp), pointer  :: solve_solAux
  integer(ip)            :: ndofn, ipoin
  integer(ip)            :: num_blocks, i_block, jblok, jzdom
  integer(ip)            :: ndofn_iblok, ndofn_jblok

  nullify(solve_solAux)
  solve_solAux => mod_module % solve(1)
  if( solve_solAux % block_num == 1 .and. solve_solAux % kfl_algso /= -999 ) then
     ndofn      = solve_solAux % ndofn
     num_blocks = solve_solAux % num_blocks ! ??
     !
     !                    if( ireaction > 0 ) then
     !call PAR_INTERFACE_NODE_EXCHANGE(solve_solAux % LPOIN_REACTION,'OR','IN CURRENT ZONE')
     !
     if( num_blocks == 1 .and. solve_solAux % kfl_react > 0 ) then
        !
        ! Monolithic system
        !
        allocate( solve_solAux % lpoin_block(npoin) )
        ndofn = solve_solAux % ndofn
        do ipoin = 1,npoin
           if( solve_solAux % LPOIN_REACTION(ipoin) ) then
              jzdom = r_dom(ipoin+1) - r_dom(ipoin)
              allocate( solve_solAux % lpoin_block(ipoin) % block2_num(1,1)                             )
              allocate( solve_solAux % lpoin_block(ipoin) % block1_num(1)                               )
              allocate( solve_solAux % lpoin_block(ipoin) % block1_num(1)   % rhs(ndofn)                )
              allocate( solve_solAux % lpoin_block(ipoin) % block2_num(1,1) % matrix(ndofn,ndofn,jzdom) )
              solve_solAux % lpoin_block(ipoin) % block1_num(1)   % rhs    = 0.0_rp
              solve_solAux % lpoin_block(ipoin) % block2_num(1,1) % matrix = 0.0_rp
           end if
        end do

     else if( solve_solAux % kfl_react > 0 ) then
        !
        ! ( n x n ) block system
        !
        allocate( solve_solAux % lpoin_block(npoin) )
        do ipoin = 1,npoin
           if( solve_solAux % lpoin_reaction(ipoin) ) then
              allocate( solve_solAux % lpoin_block(ipoin) % block2_num(num_blocks,num_blocks) )
              allocate( solve_solAux % lpoin_block(ipoin) % block1_num(num_blocks) )
              jzdom = r_dom(ipoin+1) - r_dom(ipoin)
              do i_block = 1,num_blocks
                 ndofn_iblok = solve_solAux % block_dimensions(i_block)
                 allocate( solve_solAux % lpoin_block(ipoin) % block1_num(i_block) % rhs(ndofn_iblok) )
                 solve_solAux % lpoin_block(ipoin) % block1_num(i_block) % rhs = 0.0_rp
                 do jblok = 1,num_blocks
                    ndofn_jblok = solve_solAux % block_dimensions(jblok)
                    allocate( solve_solAux % lpoin_block(ipoin) % block2_num(i_block,jblok) % matrix(ndofn_jblok,ndofn_iblok,jzdom) )
                    solve_solAux % lpoin_block(ipoin) % block2_num(i_block,jblok) % matrix = 0.0_rp
                 end do
              end do
           end if
        end do
     end if
     !
     !                    endif! ireaction
  endif ! block_num
  !
end subroutine allocate_block_system
!-------------------------------------------------------------------------||---!
!                                                                              !
!------------------------------------------------------------------------------!
