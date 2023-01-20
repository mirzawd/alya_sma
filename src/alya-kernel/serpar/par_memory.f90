!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_memory(itask)
  !-------------------------------------------------------------------------------
  !****f* parall/par_memory
  ! NAME
  !    par_memory
  ! DESCRIPTION
  !    Allocate memory for partition dimensions and arrays
  ! INPUT
  ! OUTPUT
  ! USED BY
  !    par_create_graph_arrays
  !***
  !-------------------------------------------------------------------------------
  use def_parame
  use def_domain
  use def_master
  use def_parall 
  use mod_parall, only : par_memor
  use mod_memory
  implicit none
  integer(ip), intent(in) :: itask

  if( IPARALL ) then

     select case(itask)

     case(1_ip)
        !
        ! Allocate memory for graph partition arrays
        ! 
        call memory_alloca(par_memor,'LEPAR_PAR','par_memory' , lepar_par ,  nelem )
        call memory_alloca(par_memor,'LEPER_PAR','par_memory' , leper_par ,  nelem )
        call memory_alloca(par_memor,'LEINV_PAR','par_memory' , leinv_par ,  nelem )
        call memory_alloca(par_memor,'LBPAR_PAR','par_memory' , lbpar_par ,  nboun )
        call memory_alloca(par_memor,'LBPER_PAR','par_memory' , lbper_par ,  nboun )
        call memory_alloca(par_memor,'LBINV_PAR','par_memory' , lbinv_par ,  nboun )
        call memory_alloca(par_memor,'LNPAR_PAR','par_memory' , lnpar_par ,  npoin )
        call memory_alloca(par_memor,'LNPER_PAR','par_memory' , lnper_par ,  npoin )
        call memory_alloca(par_memor,'LNINV_PAR','par_memory' , lninv_par ,  npoin )

     case(2_ip)

        call memory_alloca(par_memor,'LNINV_LOC','par_memory' , lninv_loc , npoin_total )
        call memory_alloca(par_memor,'XLNIN_LOC','par_memory' , xlnin_loc , npart_par+1_ip )

     case(3_ip)

        call memory_alloca(par_memor,'BADJ' ,'par_memory' , badj  , gnb+1_ip )
        call memory_alloca(par_memor,'BDOM' ,'par_memory' , bdom  , npoin_total-gni )
        call memory_alloca(par_memor,'BPOIN','par_memory' , bpoin , npoin_total-gni )

     case(4_ip)

        if( .not. READ_AND_RUN() ) then

           call memory_deallo(par_memor,'XLNIN_LOC','par_memory' , xlnin_loc )
           call memory_deallo(par_memor,'LNINV_LOC','par_memory' , lninv_loc )
           call memory_deallo(par_memor,'LNINV_PAR','par_memory' , lninv_par )
           call memory_deallo(par_memor,'LNPER_PAR','par_memory' , lnper_par )
           call memory_deallo(par_memor,'LNPAR_PAR','par_memory' , lnpar_par )
           call memory_deallo(par_memor,'LBINV_PAR','par_memory' , lbinv_par )
           call memory_deallo(par_memor,'LBPER_PAR','par_memory' , lbper_par )
           call memory_deallo(par_memor,'LBPAR_PAR','par_memory' , lbpar_par )
           call memory_deallo(par_memor,'LEINV_PAR','par_memory' , leinv_par )
           call memory_deallo(par_memor,'LEPER_PAR','par_memory' , leper_par )
           call memory_deallo(par_memor,'LEPAR_PAR','par_memory' , lepar_par )

        end if
        !
        ! XLNIN_LOC is used for restart: do not deallocate
        !
        call memory_deallo(par_memor,'GINDE_PAR','par_memory' , ginde_par )
        call memory_deallo(par_memor,'LNEIG_PAR','par_memory' , lneig_par )

     case(5_ip) 

        call memory_deallo(par_memor,'BADJ'      ,'par_memory' , badj      )
        call memory_deallo(par_memor,'BDOM'      ,'par_memory' , bdom      )
        call memory_deallo(par_memor,'BPOIN'     ,'par_memory' , bpoin     )
        call memory_deallo(par_memor,'NEIGHDOM'  ,'par_memory' , neighDom  )
        call memory_deallo(par_memor,'LCOMM_PAR' ,'par_memory' , lcomm_par )

     case(7_ip)

        call memory_deallo(par_memor ,'IADUAL'    ,'par_memory' , iaDual     )
        call memory_deallo(par_memor ,'JADUAL'    ,'par_memory' , jaDual     )
        call memory_deallo(par_memor ,'TRANSLDUAL','par_memory' , translDual )
        call memory_deallo(par_memor ,'COLOURS'   ,'par_memory' , colours    )

     case(8_ip)

        call memory_alloca(par_memor,'GINDE_PAR','par_memory' , ginde_par , 4_ip , npart_par+1_ip )
        call memory_alloca(par_memor,'LNEIG_PAR','par_memory' , lneig_par , npart_par      )
        call memory_alloca(par_memor,'LEIND_PAR','par_memory' , leind_par , npart_par+1_ip )
        call memory_alloca(par_memor,'LBIND_PAR','par_memory' , lbind_par , npart_par+1_ip )
        !
        ! Allocate memory for dimensions
        !
        call memory_alloca(par_memor,'NPOIN_PAR','par_memory' , npoin_par , npart_par )
        call memory_alloca(par_memor,'NELEM_PAR','par_memory' , nelem_par , npart_par )
        call memory_alloca(par_memor,'NBOUN_PAR','par_memory' , nboun_par , npart_par )
        call memory_alloca(par_memor,'SLFBO_PAR','par_memory' , slfbo_par , npart_par )

     case(11_ip)

        call memory_alloca(par_memor,'SOLVE % DISPL','par_memory' , solve_sol(1) % displ , npart_par+1_ip )
        call memory_alloca(par_memor,'SOLVE % LCOUN','par_memory' , solve_sol(1) % lcoun , npart_par+1_ip )

     case(-11_ip)

        call memory_deallo(par_memor,'SOLVE % DISPL','par_memory' , solve_sol(1) % displ )
        call memory_deallo(par_memor,'SOLVE % LCOUN','par_memory' , solve_sol(1) % lcoun )

     case(12_ip)

        call memory_alloca(par_memor,'SOLVE % XBIG'  ,'par_memory' , solve_sol(1) % xbig   , &
             solve_sol(1) % nbig             )
        call memory_alloca(par_memor,'SOLVE % LBIG'  ,'par_memory' , solve_sol(1) % lbig   , &
             solve_sol(1) % nbig             )
        call memory_alloca(par_memor,'SOLVE % XSMALL','par_memory' , solve_sol(1) % xsmall , &
             max(1_ip,solve_sol(1) % nsmall) )

     case(13_ip)

        call memory_alloca(par_memor,'SOLVE % DISP4'  ,'par_memory' , solve_sol(1) % disp4 , npart_par+1_ip )
        call memory_alloca(par_memor,'SOLVE % LCOU4'  ,'par_memory' , solve_sol(1) % lcou4 , npart_par+1_ip )

     case(-15_ip)

        call memory_deallo(par_memor,'SOLVE % LBIG','par_memory' , solve_sol(1) % lbig )

     case( 16_ip)

        call memory_alloca(par_memor,'XADJDOM','par_memory' , xadjDom , npart_par+1_ip )

     case(-16_ip)

        call memory_deallo(par_memor,'XADJDOM','par_memory' , xadjDom )
        call memory_deallo(par_memor,'ADJDOM' ,'par_memory' , adjDom  )

     case( 17_ip)

        call memory_alloca(par_memor,'ADJDOM','par_memory' , adjDom , xadjDom(npart_par+1_ip)-1_ip )

     case(18_ip)

        call memory_alloca(par_memor,'NELEM_TOT','par_memory' , nelem_tot , npart_par+1_ip )
        call memory_alloca(par_memor,'NPOIN_TOT','par_memory' , npoin_tot , npart_par+1_ip )
        call memory_alloca(par_memor,'NPOIA_TOT','par_memory' , npoia_tot , npart_par+1_ip )
        call memory_alloca(par_memor,'NBOUN_TOT','par_memory' , nboun_tot , npart_par+1_ip )

     case(20_ip)

        !
        ! Required for REREAD option to store/read lbper_par array from partition files (GROWSMARTER)
        !
        call memory_alloca(par_memor,'LBPER_PAR','par_memory' , lbper_par , nboun_total )

     end select

  end if

end subroutine par_memory
