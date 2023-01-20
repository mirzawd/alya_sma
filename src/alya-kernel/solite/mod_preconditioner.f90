!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @defgroup Preconditioner
!> @ingroup  Algebraic_Solver
!> @{
!> @name    preconditioner
!> @file    mod_preconditioner.f90
!> @author  Guillaume Houzeaux
!> @date    16/09/2015
!> @brief   ToolBox for preconditioners
!> @details ToolBox for preconditioners
!------------------------------------------------------------------------

module mod_preconditioner
  use def_domain,    only : npoin
  use def_kintyp,    only : ip,rp,lg
  use def_master,    only : INOTMASTER,npoi1,npoi2,npoi3  
  use def_solver,    only : solve_sol 
  use def_solver,    only : SYMMETRIC_GAUSS_SEIDEL 
  use def_solver,    only : STREAMWISE_GAUSS_SEIDEL
  use def_solver,    only : STREAMWISE_BIDIAGONAL
  use def_solver,    only : BIDIAGONAL             
  use def_solver,    only : memit
  use mod_memory,    only : memory_alloca
  use mod_memory,    only : memory_deallo
  use def_solver,    only : direct_solver_typ
  use mod_graphs,    only : graphs_number_along_vector

  implicit none   
  private
  real(rp), parameter :: zeror = epsilon(1.0_rp)

  public :: preconditioner_gauss_seidel
  public :: preconditioner_bidiagonal_initialize
  public :: preconditioner_bidiagonal_arrays
  public :: preconditioner_RAS_initialize

contains

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    07/10/2014
  !> @brief   Gauss-Seidel preconditioner
  !> @details According to ITASK:
  !>          ITASK = 1 ... Solve   M x = b
  !>                = 2 ... Compute   b = M x
  !>
  !>          \verbatim
  !>
  !>          Gauss-Seidel: (L+D) x = b => 
  !>          -------------
  !>
  !>                                    i-1
  !>                                    ---
  !>          x_i^{k+1} = 1/a_ii ( bi - \   a_ij * x_j^{k+1} )
  !>                                    /__
  !>                                    j=1
  !>
  !>
  !>          Streamwise Gauss-Seidel: 
  !>          ------------------------
  !>
  !>          Same as Gauss-Seidel but with renumbering
  !>          
  !>
  !>          Symmetric Gauss-Seidel: (L+D) D^-1 (U+D) x = b =>
  !>          -----------------------
  !>
  !>          1. (L+D) y = b
  !>          2. z = D y
  !>          3. (U+D) x = z
  !>      
  !>                                    i-1
  !>                                    ---
  !>          y_i^{k+1} = 1/a_ii ( bi - \   a_ij * y_j^{k+1} )
  !>                                    /__
  !>                                    j=1
  !>
  !>          z_i^{k+1} = D y_i^{k+1}
  !>                                            n
  !>                                           ---
  !>          x_i^{k+1} = 1/a_ii ( z_i^{k+1} - \   a_ij * x_j^{k+1} )
  !>                                           /__
  !>                                          j=i+1
  !>
  !>          \endverbatim
  !>
  !>          In parallel, we consider the following block preconditioning, where
  !>          i refers to internal nodes aad b to interface boundary nodes.
  !>
  !>          \verbatim
  !>          
  !>              +-        -+ +-  -+   +-  -+
  !>              | Mii  Mib | | xi |   | bi |
  !>              |          | |    | = |    | where Mbb = Diag(A)
  !>              |  0   Mbb | | xb |   | bb |
  !>              +-        -+ +-  -+   +-  -+
  !>
  !>          \endverbatim
  !> 
  !>           ITASK = 1 .... Solve M x = b:
  !>                          1. Mbb xb = bb 
  !>                          2. Mii xi = bi - Mib xb
  !>                            
  !>           ITASK = 2 .... Compute x = M b
  !>                          1. xb = Mbb b
  !>                          2. xi = Mii xi +  Mib xb
  !>
  !----------------------------------------------------------------------
  subroutine preconditioner_gauss_seidel(&
       itask,nbnodes,nbvar,aa,ja,ia,invdiag,xx,bb) 

    implicit none
    integer(ip), intent(in)          :: itask                !< What to do
    integer(ip), intent(in)          :: nbnodes              !< Number of nodes
    integer(ip), intent(in)          :: nbvar                !< Number of dof per node
    real(rp),    intent(in)          :: aa(nbvar,nbvar,*)    !< Matrix
    integer(ip), intent(in)          :: ja(*)                !< Matrix graph
    integer(ip), intent(in)          :: ia(*)                !< Matrix graph
    real(rp),    intent(in)          :: invdiag(nbvar,*)     !< Inverse diagonal
    real(rp),    intent(in)          :: bb(nbvar,*)          !< RHS
    real(rp),    intent(inout)       :: xx(nbvar,*)          !< Solution
    integer(ip)                      :: ii,jj,kk,ll
    integer(ip)                      :: ii_new,col
    real(rp),    pointer             :: yy(:,:)

    if( INOTMASTER ) then

       !-----------------------------------------------------------------
       !
       ! Subdomain interface treatment aad memory allocation: 
       !
       ! Solve xb = Mbb bb
       !
       !-----------------------------------------------------------------

       nullify(yy)        

       if( solve_sol(1) % kfl_renumbered_gs == STREAMWISE_GAUSS_SEIDEL ) then
          if( itask == 1 ) then
             do ii = npoi1+1,nbnodes
                xx(1:nbvar,ii) = bb(1:nbvar,ii) * invdiag(1:nbvar,ii)
             end do
          else          
             do ii = npoi1+1,nbnodes
                xx(1:nbvar,ii) = bb(1:nbvar,ii) / invdiag(1:nbvar,ii)
             end do
          end if
       else
          allocate(yy(nbvar,nbnodes))
          if( itask == 1 ) then
             do ii = npoi1+1,nbnodes
                xx(1:nbvar,ii) = bb(1:nbvar,ii) * invdiag(1:nbvar,ii)
                yy(1:nbvar,ii) = xx(1:nbvar,ii)
             end do
          else          
             do ii = npoi1+1,nbnodes
                xx(1:nbvar,ii) = bb(1:nbvar,ii) / invdiag(1:nbvar,ii)
                yy(1:nbvar,ii) = xx(1:nbvar,ii)
             end do
          end if
       end if

       !-----------------------------------------------------------------
       !
       ! Subdomain interface treatment aad memory allocation: 
       !
       ! Solve Mii xi = bi - Mib xb
       !
       !-----------------------------------------------------------------

       if( nbvar == 1 ) then

          !-------------------------------------------------------------
          !
          ! NBVAR = 1 
          !
          !-------------------------------------------------------------

          if( solve_sol(1) % kfl_renumbered_gs == STREAMWISE_GAUSS_SEIDEL ) then


             if( itask == 1 ) then

                !-------------------------------------------------------------
                !
                ! Streamwise Gauss-Seidel, NBVAR = 1, ITASK = 1
                !
                !-------------------------------------------------------------
                !
                ! 1. Solve (L+D) x = b
                !
                do ii = 1,npoi1
                   if( solve_sol(1) % lgrou_gs(ii) < 0 ) then 
                      xx(1,ii) = bb(1,ii) * invdiag(1,ii)
                   end if
                end do

                do ii_new = 1,npoi1
                   ii = solve_sol(1) % invpr_gs(ii_new)
                   if( solve_sol(1) % lgrou_gs(ii) > 0 ) then 
                      xx(1,ii) = bb(1,ii)
                      do jj  = ia(ii),ia(ii+1)-1
                         col = ja(jj) 
                         if( solve_sol(1) % permr_gs(col) < ii_new ) then
                            xx(1,ii) = xx(1,ii) - aa(1,1,jj) * xx(1,col) 
                         end if
                      end do
                      xx(1,ii) = xx(1,ii) * invdiag(1,ii)
                   end if
                end do

             else if( itask == 2 ) then

                !-------------------------------------------------------------
                !
                ! Streamwise Gauss-Seidel, NBVAR = 1, ITASK = 2
                !
                !-------------------------------------------------------------
                !
                ! 1. Compute x = (L+D) b
                !
                do ii = 1,npoi1
                   if( solve_sol(1) % lgrou_gs(ii) < 0 ) then 
                      xx(1,ii) = bb(1,ii) / invdiag(1,ii)
                   end if
                end do

                do ii_new = 1,npoi1
                   ! x = D b
                   ii = solve_sol(1) % invpr_gs(ii_new)
                   if( solve_sol(1) % lgrou_gs(ii) > 0 ) then                       
                      xx(1,ii) = bb(1,ii) / invdiag(1,ii) 
                      ! x = x + L b
                      do jj  = ia(ii),ia(ii+1)-1
                         col = ja(jj) 
                         if( solve_sol(1) % permr_gs(col) < ii_new ) then
                            xx(1,ii) = xx(1,ii) + aa(1,1,jj) * bb(1,col) 
                         end if
                      end do
                   end if
                end do

             end if

          else 

             !---------------------------------------------------------------
             !
             ! Gauss-Seidel aad symmetric Gauss-Seidel (SSOR), NBVAR = 1, ITASK = 1
             !
             !---------------------------------------------------------------
             !
             ! Gauss-Seidel:
             !
             ! 1. (L+D) y = b
             !
             ! Symmetric Gauss-Seidel:
             !
             ! (L+D) D^-1 (U+D) x = b => x = (U+D)^-1 D (L+D)^-1 b
             !                                          <-------->
             !                                               y
             !                                        <---------->
             !                                             z
             !                               <------------------->
             !                                           x
             ! 1. (L+D) y = b
             ! 2. z = D y
             ! 3. (U+D) x = z
             !
             if( itask == 1 ) then
                !
                ! (L+D) y = b
                !
                do ii = 1,npoi1
                   xx(1,ii) = bb(1,ii)
                   do jj  = ia(ii),ia(ii+1)-1
                      col = ja(jj) 
                      if( col < ii ) then
                         xx(1,ii) = xx(1,ii) - aa(1,1,jj) * xx(1,col) 
                      end if
                   end do
                   xx(1,ii) = xx(1,ii) * invdiag(1,ii)
                end do
                
                if( solve_sol(1) % kfl_renumbered_gs == SYMMETRIC_GAUSS_SEIDEL ) then
                   !
                   ! 2. z = D y
                   !              
                   do ii = 1,npoi1
                      yy(1,ii) = xx(1,ii) / invdiag(1,ii)
                   end do
                   !
                   ! 3. (U+D) x = z
                   !
                   do ii = npoi1,1,-1
                      xx(1,ii) = yy(1,ii) 
                      do jj  = ia(ii),ia(ii+1)-1
                         col = ja(jj) 
                         if( col > ii ) then
                            xx(1,ii) = xx(1,ii) - aa(1,1,jj) * xx(1,col) 
                         end if
                      end do
                      xx(1,ii) = xx(1,ii) * invdiag(1,ii)
                   end do
                end if

             else if( itask == 2 ) then

                !---------------------------------------------------------------
                !
                ! Gauss-Seidel aad symmetric Gauss-Seidel (SSOR), NBVAR = 1, ITASK = 2
                !
                !---------------------------------------------------------------
                !
                ! Gauss-Seidel:
                !
                ! 1. x = (L+D) b
                !
                ! Symmetric Gauss-Seidel:
                !
                ! x = (L+D) D^-1 (U+D) b
                ! 1. y = (L+D) b
                ! 2. y = D^-1 y
                ! 3. x = (U+D) y 
                !
                do ii = 1,npoi1
                   ! y = D b
                   xx(1,ii) = bb(1,ii) / invdiag(1,ii)
                   ! y = y + L b
                   do jj  = ia(ii),ia(ii+1)-1
                      col = ja(jj) 
                      if( col < ii ) then
                         xx(1,ii) = xx(1,ii) + aa(1,1,jj) * bb(1,col) 
                      end if
                   end do
                end do

                if( solve_sol(1) % kfl_renumbered_gs == SYMMETRIC_GAUSS_SEIDEL ) then
                   !
                   ! 2. y = D^-1 y
                   !              
                   do ii = 1,npoi1
                      yy(1,ii) = xx(1,ii) * invdiag(1,ii)
                   end do
                   !
                   ! 3. x = (U+D) y
                   !
                   do ii = npoi1,1,-1
                      xx(1,ii) = yy(1,ii) / invdiag(1,ii)
                      do jj  = ia(ii),ia(ii+1)-1
                         col = ja(jj) 
                         if( col > ii ) then
                            xx(1,ii) = xx(1,ii) + aa(1,1,jj) * yy(1,col)
                         end if
                      end do
                   end do

                end if

             end if

          end if

       else if( nbvar > 1 ) then

          !------------------------------------------------------------------
          !
          ! NBVAR > 1
          !
          !------------------------------------------------------------------

          if( solve_sol(1) % kfl_renumbered_gs == STREAMWISE_GAUSS_SEIDEL ) then

             if( itask == 1 ) then

                !---------------------------------------------------------------
                !
                ! Streamwise Gauss-Seidel, , NBVAR > 1, ITASK = 1
                !
                !---------------------------------------------------------------
                do ii = 1, npoi1
                   if(solve_sol(1)%lgrou_gs(ii)<0) then
                      xx(1:nbvar,ii) = bb(1:nbvar,ii) * invdiag(1:nbvar,ii)
                   end if
                end do

                do ii_new = 1,npoi1
                   ii = solve_sol(1) % invpr_gs(ii_new)
                   if( solve_sol(1)%lgrou_gs(ii)>0 ) then
                      xx(1:nbvar,ii) = bb(1:nbvar,ii)
                      do jj  = ia(ii),ia(ii+1)-1
                         col = ja(jj) 
                         if( solve_sol(1) % permr_gs(col) < ii_new ) then
                            do ll = 1,nbvar
                               do kk = 1,nbvar
                                  xx(kk,ii) = xx(kk,ii) - aa(ll,kk,jj) * xx(ll,col) 
                               end do
                            end do
                         end if
                      end do
                      xx(1:nbvar,ii) = xx(1:nbvar,ii) * invdiag(1:nbvar,ii)
                   end if
                end do

             else if( itask == 2 ) then

                !-------------------------------------------------------------
                !
                ! Streamwise Gauss-Seidel, NBVAR > 1, ITASK = 2
                !
                !-------------------------------------------------------------
                do ii = 1, npoi1
                   if( solve_sol(1)%lgrou_gs(ii)<0) then
                      xx(1:nbvar,ii) = bb(1:nbvar,ii)/invdiag(1:nbvar,ii)
                   end if
                end do

                do ii_new = 1,npoi1
                   ! x = D b
                   ii = solve_sol(1) % invpr_gs(ii_new)

                   if ( solve_sol(1)%lgrou_gs(ii)>0 ) then 
                      xx(1:nbvar,ii) = bb(1:nbvar,ii) / invdiag(1:nbvar,ii) 
                      !
                      ! x = x + L b
                      do jj  = ia(ii),ia(ii+1)-1
                         col = ja(jj) 
                         if( solve_sol(1) % permr_gs(col) < ii_new ) then
                            do ll = 1,nbvar
                               do kk = 1,nbvar                              
                                  xx(kk,ii) = xx(kk,ii) + aa(ll,kk,jj) * bb(ll,col) 
                               end do
                            end do
                         end if
                      end do
                   end if
                end do

             end if

          else 

             !---------------------------------------------------------------
             !
             ! Gauss-Seidel aad symmetric Gauss-Seidel (SSOR), NBVAR > 1, ITASK = 1
             !
             !---------------------------------------------------------------

             if( itask == 1 ) then
                !
                ! (L+D) y = b
                !
                do ii = 1,npoi1
                   xx(1:nbvar,ii) = bb(1:nbvar,ii)
                   do jj  = ia(ii),ia(ii+1)-1
                      col = ja(jj) 
                      if( col < ii ) then
                         do ll = 1,nbvar
                            do kk = 1,nbvar
                               xx(kk,ii) = xx(kk,ii) - aa(ll,kk,jj) * xx(ll,col) 
                            end do
                         end do
                      end if
                   end do
                   xx(1:nbvar,ii) = xx(1:nbvar,ii) * invdiag(1:nbvar,ii)
                end do

                if( solve_sol(1) % kfl_renumbered_gs == SYMMETRIC_GAUSS_SEIDEL ) then
                   !
                   ! 2. z = D y
                   !              
                   do ii = 1,npoi1
                      yy(1:nbvar,ii) = xx(1:nbvar,ii) / invdiag(1:nbvar,ii)
                   end do
                   !
                   ! 3. (U+D) x = z
                   !
                   do ii = npoi1,1,-1
                      xx(1:nbvar,ii) = yy(1:nbvar,ii) 
                      do jj  = ia(ii),ia(ii+1)-1
                         col = ja(jj) 
                         if( col > ii ) then
                            do ll = 1,nbvar
                               do kk = 1,nbvar
                                  xx(kk,ii) = xx(kk,ii) - aa(ll,kk,jj) * xx(ll,col) 
                               end do
                            end do
                         end if
                      end do
                      xx(1:nbvar,ii) = xx(1:nbvar,ii) * invdiag(1:nbvar,ii)
                   end do
                end if

             else if( itask == 2 ) then

                !---------------------------------------------------------------
                !
                ! Gauss-Seidel aad symmetric Gauss-Seidel (SSOR), NBVAR > 1, ITASK = 2
                !
                !---------------------------------------------------------------

                do ii = 1,npoi1
                   ! y = D b
                   xx(1:nbvar,ii) = xx(1:nbvar,ii) / invdiag(1:nbvar,ii)
                   ! y = y + L b
                   do jj  = ia(ii),ia(ii+1)-1
                      col = ja(jj) 
                      if( col < ii ) then
                         do ll = 1,nbvar
                            do kk = 1,nbvar
                               xx(kk,ii) = xx(kk,ii) + aa(ll,kk,jj) * bb(ll,col) 
                            end do
                         end do
                      end if
                   end do
                end do

                if( solve_sol(1) % kfl_renumbered_gs == SYMMETRIC_GAUSS_SEIDEL ) then
                   !
                   ! 2. y = D^-1 y
                   !              
                   do ii = 1,npoi1
                      yy(1:nbvar,ii) = xx(1:nbvar,ii) * invdiag(1:nbvar,ii)
                   end do
                   !
                   ! 3. x = (U+D) y
                   !
                   do ii = npoi1,1,-1
                      xx(1:nbvar,ii) = yy(1:nbvar,ii) / invdiag(1:nbvar,ii)
                      do jj  = ia(ii),ia(ii+1)-1
                         col = ja(jj) 
                         if( col > ii ) then
                            do ll = 1,nbvar
                               do kk = 1,nbvar
                                  xx(kk,ii) = xx(kk,ii) + aa(ll,kk,jj) * yy(ll,col) 
                               end do
                            end do
                         end if
                      end do
                   end do
                end if

             end if

          end if

       end if

       if( associated(yy) ) deallocate(yy)

    end if

  end subroutine preconditioner_gauss_seidel

  subroutine preconditioner_bidiagonal_arrays(nbnodes,nbvar,aa,ja,ia)

    integer(ip), intent(in)             :: nbnodes
    integer(ip), intent(in)             :: nbvar
    real(rp),    intent(in)             :: aa(nbvar,nbvar,*)
    integer(ip), intent(in)             :: ja(*)
    integer(ip), intent(in)             :: ia(*)
    integer(ip),                pointer :: idiag1(:,:)
    real(rp),                   pointer :: adiag1(:,:,:)
    integer(ip)                         :: izdom
    integer(ip)                         :: kpoin,mpoin,ipoin,igrou,kgrou,jgrou,jpoin,mgrou,ngrou
    integer(ip)                         :: kk, ll

    nullify(adiag1)
    nullify(idiag1)

    if( INOTMASTER ) then

       ngrou = solve_sol(1)%ngrou_gs
       idiag1 => solve_sol(1)%idiag1
       adiag1 => solve_sol(1)%adiag1

       if( nbvar == 1 ) then

          jgrou = 0
          do kpoin = 1,npoi1    

             ipoin = solve_sol(1) % invpr_gs(kpoin) ! ipoin(old) --> kpoin(new)
             igrou = solve_sol(1) % lgrou_gs(ipoin)    
    
             if( igrou == jgrou .and. igrou > 0 ) then

                izdom = ia(ipoin)
                mpoin = ja(izdom)
                do while( mpoin /= jpoin .and. izdom < ia(ipoin+1)-1 )
                   izdom = izdom + 1
                   mpoin = ja(izdom)
                end do
                if( mpoin /= jpoin ) then
                   print*,'GROUPS=',igrou,jgrou
                   print*,'NODES =',ipoin,jpoin,kpoin
                   call runend('CANNOT FIND POINT')
                end if

                adiag1(1,1,kpoin) = aa(1,1,izdom) 
                idiag1(1,  kpoin) = jpoin
                
             else
                !
                ! Seed
                !
                idiag1(1,  kpoin) = kpoin  
                adiag1(1,1,kpoin) = 0.0_rp
                jgrou             = igrou
             end if

             jpoin = ipoin 

          end do

          do ipoin = 1,npoi1
             if( idiag1(1,ipoin) == 0 ) then
                print*,'bordel de merde'
                stop
             end if
          end do


          return

          do igrou = 1, ngrou           
             kgrou = 0            
             do kpoin = 1, npoi1        ! renumbered npoi1 

                if( idiag1(1_ip,kpoin) == 0 ) then     
    
                   ipoin = solve_sol(1) % invpr_gs(kpoin) !ipoin(old)-->kpoin(new)
                   mgrou = solve_sol(1) % lgrou_gs(ipoin)        

                   if( solve_sol(1) % lgrou_gs(ipoin) <= 0 ) then                      

                      idiag1(1_ip,kpoin)  = kpoin                               
                      adiag1(1,1,kpoin)   = 0.0_rp
                      
                   else if( solve_sol(1) % lgrou_gs(ipoin) > 0 ) then                      

                      if( mgrou == igrou ) then

                         if( kgrou == 0 ) then

                            idiag1(1_ip,kpoin)  =  kpoin  
                            adiag1(1,1,kpoin)   =  0.0_rp
                            izdom = ia(ipoin)
                            mpoin = ja(izdom)

                         else  

                            izdom = ia(ipoin)
                            mpoin = ja(izdom)

                            do while (mpoin /= jpoin)
                               izdom = izdom + 1
                               mpoin = ja(izdom)
                            end do


                            adiag1(1,1,kpoin)   = aa(1,1,izdom) 
                            idiag1(1_ip,kpoin)  = mpoin                                   


                         end if
                         kgrou = kgrou + 1
                         jpoin = ipoin
                      end if
                   end if
                end if
             end do
          end do
          do ipoin = 1,npoi1
             if( idiag1(1_ip,ipoin) == 0 ) stop
             idiag1(1_ip,ipoin) = abs(idiag1(1_ip,ipoin))
          end do
!!$          do igrou = 1, ngrou           
!!$             kgrou = 0            
!!$             do kpoin = 1, npoi1        ! renumbered npoi1 
!!$
!!$                if(idiag1(1_ip,kpoin) == 0) then         
!!$                   ipoin = solve_sol(1) % invpr_gs(kpoin) !ipoin(old)-->kpoin(new)
!!$                   mgrou = solve_sol(1) % lgrou_gs(ipoin)        
!!$                   if( solve_sol(1) % lgrou_gs(ipoin) > 0 ) then                      
!!$
!!$                      if (mgrou == igrou) then
!!$                         if( kgrou == 0 ) then
!!$
!!$                            idiag1(1_ip,kpoin)       =-ipoin  
!!$                            adiag1(1,1,kpoin)      =  -1_rp
!!$                            izdom = ia(ipoin)
!!$                            mpoin = ja(izdom)
!!$
!!$                         else  
!!$
!!$                            izdom = ia(ipoin)
!!$                            mpoin = ja(izdom)
!!$
!!$                            do while (mpoin /= jpoin)
!!$                               izdom = izdom + 1
!!$                               mpoin = ja(izdom)
!!$                            end do
!!$
!!$
!!$                            adiag1(1,1,kpoin)   = aa(1,1,izdom) 
!!$                            idiag1(1_ip,kpoin)  = mpoin                                   
!!$
!!$
!!$                         end if
!!$                         kgrou = kgrou + 1
!!$                         jpoin = ipoin
!!$                      end if
!!$                   end if
!!$                end if
!!$             end do
!!$          end do

       else if(nbvar > 1) then



          do ipoin = 1, npoi1
             idiag1(1,ipoin) = 0
             adiag1(1:nbvar,1:nbvar,ipoin) = 0
          end do

          do igrou = 1, ngrou           
             kgrou = 0            
             do kpoin = 1, npoi1        ! renumbered npoi1 

                if(idiag1(1_ip,kpoin) == 0) then         
                   ipoin = solve_sol(1) % invpr_gs(kpoin) !ipoin(old)-->kpoin(new)
                   mgrou = solve_sol(1) % lgrou_gs(ipoin)                       
                     if( mgrou > 0 ) then
                      if (mgrou == igrou) then
                         if( kgrou == 0 ) then

                            idiag1(1_ip,kpoin)       =-ipoin  
                            adiag1(1:nbvar,1:nbvar,kpoin)      =-1_rp
                            izdom = ia(ipoin)
                            mpoin = ja(izdom)

                         else  

                            izdom = ia(ipoin)
                            mpoin = ja(izdom)

                            do while (mpoin /= jpoin)
                               izdom = izdom + 1
                               mpoin = ja(izdom)
                            end do

                            do ll = 1, nbvar
                               do kk = 1, nbvar
                                  adiag1(kk,ll,kpoin) = aa(ll,kk,izdom) 
                                  idiag1(1_ip,kpoin)  = mpoin                                   
                               end do
                            end do

                         end if
                         kgrou = kgrou + 1
                         jpoin = ipoin
                      end if
                  end if 
                end if
             end do
          end do
       end if !nbvar
    end if !INOTMASTER

  end subroutine preconditioner_bidiagonal_arrays

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux, Juan Carlos Cajas, Paula Córdoba
  !> @date    04/02/2015
  !> @brief   Bidiagonal solver initialization
  !> @details Allocate aad compute some arrays for the bidiagonal 
  !>          preconditioner in two ways:
  !>          1. 'Streamwise' ----> Mesh renumbering along 'advection'
  !>          2. 'Non-Streamwise'-> Non-renumbered mesh
  !>          According to ITASK:
  !>          ITASK = 1 ... Solve   M x = b
  !>                = 2 ... Compute   b = M x
  !>
  !>          \verbatim
  !>           
  !>          In both cases an iteration Gauss-Seidel is uzed where the lower matrix L corresponds to the subdiagonal in each case          !>          either the streamwise or the non-streamwise.
  !>
  !>          Gauss-Seidel: (L+D) x = b => 
  !>          -------------
  !>
  !>                                    i-1
  !>                                    ---
  !>          x_i^{k+1} = 1/a_ii ( bi - \   a_ij * x_j^{k+1} )
  !>                                    /__
  !>                                    j=1
  !>         \end verbatim
  !>
  !----------------------------------------------------------------------

  subroutine preconditioner_bidiagonal_initialize(itask,nbnodes,nbvar,invdiag,xx,bb)
    integer(ip), intent(in)             :: itask
    integer(ip), intent(in)             :: nbnodes
    integer(ip), intent(in)             :: nbvar
    real(rp),    intent(in)             :: bb(nbvar,*)
    real(rp),    intent(in)             :: invdiag(nbvar,*)
    integer(ip),                pointer :: idiag1(:,:)
    real(rp),                   pointer :: adiag1(:,:,:)
    real(rp),    intent(inout)          :: xx(nbvar,*)
    integer(ip)                         :: ipoin, kpoin
    integer(ip)                         :: kk, ll



    if( INOTMASTER ) then

       idiag1 => solve_sol(1)%idiag1
       adiag1 => solve_sol(1)%adiag1

       if( solve_sol(1) % kfl_renumbered_gs == STREAMWISE_BIDIAGONAL ) then
          if( itask == 1 ) then
             do ipoin = npoi1+1,nbnodes 
                xx(1:nbvar,ipoin) = bb(1:nbvar,ipoin) * invdiag(1:nbvar,ipoin)
             end do
          else          
             do ipoin = npoi1+1,nbnodes
                xx(1:nbvar,ipoin) = bb(1:nbvar,ipoin) / invdiag(1:nbvar,ipoin)
             end do
          end if
       else
          if( itask == 1 ) then
             do ipoin = npoi1+1,nbnodes
                xx(1:nbvar,ipoin) = bb(1:nbvar,ipoin) * invdiag(1:nbvar,ipoin)              
             end do
          else          
             do ipoin = npoi1+1,nbnodes
                xx(1:nbvar,ipoin) = bb(1:nbvar,ipoin) / invdiag(1:nbvar,ipoin)             
             end do
          end if
       end if

       if ( nbvar == 1 ) then

          if(solve_sol(1)%kfl_renumbered_gs == STREAMWISE_BIDIAGONAL) then

             if(itask == 1) then
                !
                ! Diagonal Preconditioner for groups with a single node
                !
                !do ipoin = 1,npoi1
                !   if( solve_sol(1) % lgrou_gs(ipoin) < 0 ) then
                !      xx(1,ipoin) = bb(1,ipoin) * invdiag(1,ipoin)
                !   end if
                !end do

                !do kpoin =  1, npoi1
                !   ipoin = solve_sol(1)%invpr_gs(kpoin)  
                !   if( solve_sol(1) % lgrou_gs(ipoin) > 0 ) then        
                !      xx(1,ipoin) = bb(1,ipoin)
                !      if( idiag1(1,kpoin) > 0 )then
                !         xx(1,ipoin) = xx(1,ipoin) - adiag1(1,1,kpoin)*xx(1,idiag1(1,kpoin))
                !      end if
                !      xx(1,ipoin) = xx(1,ipoin) * invdiag(1,ipoin) 
                !   end if
                !end do

                do ipoin = 1,npoi1
                   if( solve_sol(1) % lgrou_gs(ipoin) <= 0 ) then
                      xx(1,ipoin) = bb(1,ipoin) * invdiag(1,ipoin)
                   end if
                end do

                do kpoin = 1,npoi1
                   ipoin = solve_sol(1) % invpr_gs(kpoin) 
                   if( solve_sol(1) % lgrou_gs(ipoin) > 0 ) then
                      xx(1,ipoin) = invdiag(1,ipoin) * ( bb(1,ipoin) - adiag1(1,1,kpoin) * xx(1,idiag1(1,kpoin)) )
                   end if
                end do

             else if( itask == 2 ) then 

                do ipoin = 1,npoi1
                   if( solve_sol(1) % lgrou_gs(ipoin) <= 0 ) then
                      xx(1,ipoin) = bb(1,ipoin) / invdiag(1,ipoin)
                   end if
                end do

                do kpoin = 1,npoi1                   
                   ipoin = solve_sol(1)%invpr_gs(kpoin)   
                   if( solve_sol(1) % lgrou_gs(ipoin) > 0 ) then
                      xx(1,ipoin) = ( bb(1,ipoin) + adiag1(1,1,kpoin) * bb(1,idiag1(1,kpoin)) ) / invdiag(1,ipoin)
                   end if
                end do

             end if
          end if

       else if ( nbvar > 1 ) then

          if(solve_sol(1)%kfl_renumbered_gs == STREAMWISE_BIDIAGONAL) then

             if (itask == 1) then
                !
                ! Diagonal Preconditioner for groups with a single node
                !
                do ipoin = 1,npoi1
                   if( solve_sol(1) % lgrou_gs(ipoin) < 0 ) then
                      xx(1:nbvar,ipoin) = bb(1:nbvar,ipoin) * invdiag(1:nbvar,ipoin)
                   end if
                end do


                do kpoin =  1, npoi1 


                   ipoin = solve_sol(1)%invpr_gs(kpoin)   
                   if( solve_sol(1) % lgrou_gs(ipoin) > 0 ) then                           
                      xx(1:nbvar,ipoin) = bb(1:nbvar,ipoin)                 

                      if(idiag1(1_ip,kpoin) > 1) then 
                         do ll = 1, nbvar
                            do kk = 1, nbvar
                               xx(kk,ipoin) = xx(kk,ipoin) - adiag1(kk,ll,kpoin)*xx(ll,idiag1(1_ip,kpoin))
                            end do
                         end do
                      end if

                      xx(1:nbvar,ipoin) = xx(1:nbvar,ipoin) * invdiag(1:nbvar,ipoin) 
                   end if

                end do

             else if ( itask == 2 ) then               
                !
                ! Diagonal Preconditioner for groups with a single node
                !
                do ipoin = 1,npoi1
                   if( solve_sol(1) % lgrou_gs(ipoin) < 0 ) then
                      xx(1:nbvar,ipoin) = bb(1:nbvar,ipoin) / invdiag(1:nbvar,ipoin)
                   end if
                end do


                do kpoin =  1, npoi1                   
                   ipoin = solve_sol(1)%invpr_gs(kpoin)   
                   if( solve_sol(1) % lgrou_gs(ipoin) > 0 ) then                           
                      xx(1:nbvar,ipoin) = bb(1:nbvar,ipoin)                 

                      if(idiag1(1_ip,kpoin) > 1) then 
                         do ll = 1, nbvar
                            do kk = 1, nbvar
                               xx(kk,ipoin) = xx(kk,ipoin) + adiag1(kk,ll,kpoin)*bb(ll,idiag1(1_ip,kpoin))
                            end do
                         end do
                      end if
                      xx(1:nbvar,ipoin) = xx(1:nbvar,ipoin) / invdiag(1:nbvar,ipoin) 
                   end if

                end do

             end if

          end if
       end if
    end if



  end subroutine preconditioner_bidiagonal_initialize


  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    07/10/2014
  !> @brief   Initialize RAS preconditioner
  !> @details Initialize RAS preconditioner
  !>
  !----------------------------------------------------------------------

  subroutine preconditioner_RAS_initialize(ndof,ia,ja,direct_solver)
    integer(ip),             intent(in)          :: ndof
    integer(ip),             intent(in), pointer :: ia(:) 
    integer(ip),             intent(in), pointer :: ja(:) 
    type(direct_solver_typ), intent(inout)       :: direct_solver 

    direct_solver % ia    => ia
    direct_solver % ja    => ja
    direct_solver % nn    =  size(ia,KIND=ip)-1
    direct_solver % ndof  =  ndof
    direct_solver % nz    =  ia(direct_solver % nn+1)-1
    if( size(ja,KIND=ip) /= direct_solver % nz ) call runend('preconditioner_RAS_initialize: WRONG GRAPH')

    !call solver_direct_solver_initialize(direct_solver)

  end subroutine preconditioner_RAS_initialize

end module mod_preconditioner
!> @}

