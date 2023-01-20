!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine skygro()
  !-----------------------------------------------------------------------
  !****f* domain/skygro
  ! NAME
  !    skygro
  ! DESCRIPTION
  !    Set up the skyline or sparse format for deflated solvers
  ! USED BY
  !    domain
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_solver
  use mod_memory
  use mod_postpr
  use mod_communications_global, only : PAR_MIN
  implicit none
  integer(ip)          :: igrou,kskyl,kgrou,idofn,jdofn,ndofn,ipoin
  integer(ip)          :: izdom,jgrou,jpoin,ngrou,mgrou
  integer(ip), pointer :: iskyl(:)
  integer(ip), pointer :: idiag(:)
  integer(ip), pointer :: lgrou(:) 

  nullify(iskyl)
  nullify(idiag)
  nullify(lgrou)

  ndofn =  solve_sol(1) % ndofn
  ngrou =  solve_sol(1) % ngrou
  lgrou => solve_sol(1) % lgrou

  !----------------------------------------------------------------------
  !
  ! Skyline format
  !
  !----------------------------------------------------------------------

  if( ndofn == 1 ) then
     !
     ! NDOFN = 1
     !
     call memory_alloca(memit,'SOLVE % ISKYL','skygro',solve_sol(1) % iskyl,ngrou+1_ip)
     iskyl => solve_sol(1) % iskyl

     do igrou = 1,ngrou+1
        iskyl(igrou) = ngrou
     end do

     if( INOTMASTER ) then
        do ipoin = 1,npoin
           igrou = solve_sol(1) % lgrou(ipoin)
           if( igrou > 0 ) then
              do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                 jpoin = c_dom(izdom)           
                 jgrou = solve_sol(1) % lgrou(jpoin)
                 if( jgrou > 0 ) then
                    if( igrou >= jgrou .and. jgrou < iskyl(igrou+1) ) iskyl(igrou+1) = jgrou
                 end if
              end do
           end if
        end do
     end if

  else
     !
     ! NDOFN > 1
     !
     call memory_alloca(memit,'SOLVE % ISKYL','skygro',solve_sol(1) % iskyl,ndofn*ngrou+1)
     iskyl => solve_sol(1) % iskyl

     if( INOTMASTER ) then
        do igrou = 1,(ngrou*ndofn)+1
           iskyl(igrou) = ngrou*ndofn
        end do
        do ipoin = 1,npoin
           igrou = solve_sol(1) % lgrou(ipoin)
           if( igrou > 0 ) then
              do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                 jpoin = c_dom(izdom)           
                 jgrou = solve_sol(1) % lgrou(jpoin)
                 if( jgrou > 0 .and. igrou >= jgrou ) then
                    do idofn = 1,ndofn 
                       kgrou = (igrou-1)*ndofn+idofn+1  
                       do jdofn = 1,ndofn
                          mgrou = (jgrou-1)*ndofn+jdofn
                          if( mgrou < iskyl(kgrou) ) iskyl(kgrou) = mgrou
                       end do
                    end do
                 end if
              end do
           end if
        end do
     end if

  end if

  call PAR_MIN(ndofn*ngrou+1_ip,iskyl)

  solve_sol(1) % nskyl = 1
  iskyl(1)             = 1

  if(  solve_sol(1) % kfl_symme == 1                                .or. &
       solve_sol(1) % kfl_algso == SOL_SOLVER_CG                    .or. &
       solve_sol(1) % kfl_algso == SOL_SOLVER_DEFLATED_CG           .or. &
       solve_sol(1) % kfl_algso == SOL_SOLVER_PIPELINED_CG          .or. &
       solve_sol(1) % kfl_algso == SOL_SOLVER_A_DEF2                .or. &
       solve_sol(1) % kfl_algso == SOL_SOLVER_PIPELINED_DEFLATED_CG ) then
     ! 
     ! For the symmetric case, do not need idiag
     !
     if( ndofn == 1 ) then 

        do igrou = 1,ngrou
           kskyl                = igrou - iskyl(igrou+1) + 1
           solve_sol(1) % nskyl = solve_sol(1) % nskyl + kskyl
           iskyl(igrou+1)       = solve_sol(1) % nskyl
        end do

     else

        kgrou = 0_ip 
        do igrou = 1,ngrou
           do idofn = 1,ndofn
              kgrou                = kgrou + 1
              kskyl                = kgrou - iskyl(kgrou+1) + 1
              solve_sol(1) % nskyl = solve_sol(1) % nskyl + kskyl
              iskyl(kgrou+1)       = solve_sol(1) % nskyl
           end do
        end do

     end if

  else
     !
     ! For the nonsymmetric case, set idiag 
     !
     if( ndofn == 1 ) then

        call memory_alloca(memit,'SOLVE % IDIAG','skygro',solve_sol(1) % idiag,ngrou)
        idiag => solve_sol(1) % idiag
        do igrou = 1,ngrou
           kskyl                = igrou - iskyl(igrou+1)
           idiag(igrou)         = solve_sol(1) % nskyl + kskyl
           kskyl                = 2 * kskyl + 1  
           solve_sol(1) % nskyl = solve_sol(1) % nskyl + kskyl
           iskyl(igrou+1)       = solve_sol(1) % nskyl
        end do

     else

        call memory_alloca(memit,'SOLVE % IDIAG','skygro',solve_sol(1) % idiag,ngrou*ndofn)
        idiag  => solve_sol(1) % idiag

        kgrou = 0_ip
        do igrou = 1,ngrou
           do idofn = 1,ndofn 
              kgrou                = kgrou + 1
              kskyl                = kgrou - iskyl(kgrou+1)
              idiag(kgrou)         = solve_sol(1) % nskyl + kskyl
              kskyl                = 2 * kskyl + 1  
              solve_sol(1) % nskyl = solve_sol(1) % nskyl + kskyl
              iskyl(kgrou+1)       = solve_sol(1) % nskyl
           enddo
        end do

     end if

  end if

  solve_sol(1) % nskyl = solve_sol(1) % nskyl - 1 

end subroutine skygro

