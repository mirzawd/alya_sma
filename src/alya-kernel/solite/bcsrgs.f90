!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Algebraic_Solver
!> @addtogroup Krylov_Solver
!> @{                                                                   
!> @file    all_precon.f90
!> @author  Guillaume Houzeaux
!> @date    15/07/2015
!> @brief   Gauss-Seidel type preconditioning
!> @details Solve M x = b where:
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
!> @}                                                                   
!-----------------------------------------------------------------------

subroutine bcsrgs(nbnodes,nbvar,an,ja,ia,invdiag,xx,bb) 
  !----------------------------------------------------------------------
  !****f* mathru/bcsrax
  ! NAME 
  !     bcsrax
  ! DESCRIPTION
  !     Solve (L+D) D^-1 (U+D) x = b
  ! INPUT
  !    NBNODES .... Number of equations
  !    NBVAR ...... Number of variables
  !    AN ......... Matrix
  !    JA ......... List of elements
  !    IA ......... Pointer to list of elements
  !    BB ......... Vector
  ! OUTPUT
  !    XX ......... result vector
  ! USES  
  ! USED BY
  !***
  !----------------------------------------------------------------------
  use def_kintyp, only             :  ip,rp,lg
  use def_master, only             :  INOTMASTER,npoi1
  use def_solver, only             :  solve_sol
  implicit none
  integer(ip), intent(in)          :: nbnodes,nbvar
  real(rp),    intent(in)          :: an(nbvar,nbvar,*)
  integer(ip), intent(in)          :: ja(*),ia(*)
  real(rp),    intent(in)          :: invdiag(nbvar,*)
  real(rp),    intent(in)          :: bb(nbvar,*)
  real(rp),    intent(inout)       :: xx(nbvar,*)
  integer(ip)                      :: ii,jj,kk,ll,col,ii_new
!  integer(ip)                      :: igrou,jgrou
  real(rp),    pointer             :: yy(:,:)

  if( INOTMASTER ) then

     allocate(yy(nbvar,nbnodes))

     do ii = npoi1+1,nbnodes
        xx(1:nbvar,ii) = bb(1:nbvar,ii) * invdiag(1:nbvar,ii)
        yy(1:nbvar,ii) = xx(1:nbvar,ii)
     end do

     if( nbvar == 1 ) then

        if( solve_sol(1) % kfl_renumbered_gs == 1 ) then

           !-------------------------------------------------------------
           !
           ! Streamwise Gauss-Seidel
           !
           !-------------------------------------------------------------
           !
           ! 1. (L+D) y = b
           !
           !*OMP   PARALLEL DO SCHEDULE (GUIDED)               & 
           !*OMP   DEFAULT (NONE)                              &
           !*OMP   PRIVATE ( ii, jj, kk, ll, col, raux)        &
           !*OMP   SHARED ( nbnodes, nbvar, bb, xx, ia, ja, an)
           do ii_new = 1,npoi1
              ii       = solve_sol(1) % invpr_gs(ii_new)
              !igrou    = solve_sol(1) % lgrou_gs(ii)
              yy(1,ii) = bb(1,ii)
              do jj  = ia(ii),ia(ii+1)-1
                 col = ja(jj) 
                 !jgrou = solve_sol(1) % lgrou_gs(col)
                 !if( solve_sol(1) % permr_gs(col) < ii_new .and. igrou == jgrou ) then
                 if( solve_sol(1) % permr_gs(col) < ii_new ) then
                    yy(1,ii) = yy(1,ii) - an(1,1,jj) * yy(1,col) 
                 end if
              end do
              yy(1,ii) = yy(1,ii) * invdiag(1,ii)
           end do
           do ii = 1,npoi1
              xx(1,ii) = yy(1,ii) 
           end do
           goto 10

        else 

           !-------------------------------------------------------------
           !
           ! Gauss-Seidel and symmetric Gauss-Seidel (SSOR)
           !
           !-------------------------------------------------------------
           !
           ! Gauss-Seidel:
           !
           ! 1. (L+D) y = b
           !
           ! Symmetric Gauss-Seidel:
           !
           ! (L+D) D^-1 (U+D) x = b
           ! 1. (L+D) y = b
           ! 2. z = D y
           ! 3. (U+D) x = z
           !
           do ii = 1,npoi1
              yy(1,ii) = bb(1,ii)
              do jj  = ia(ii),ia(ii+1)-1
                 col = ja(jj) 
                 if( col < ii ) then
                    yy(1,ii) = yy(1,ii) - an(1,1,jj) * yy(1,col) 
                 end if
              end do
              yy(1,ii) = yy(1,ii) * invdiag(1,ii)
           end do
           if( solve_sol(1) % kfl_renumbered_gs == 0 ) then
              do ii = 1,npoi1
                 xx(1,ii) = yy(1,ii) 
              end do
              goto 10 
           else           
              !
              ! 2. z = D y
              !              
              do ii = 1,npoi1
                 yy(1,ii) = yy(1,ii) / invdiag(1,ii)
              end do
              !
              ! 3. (U+D) x = z
              !
              do ii = npoi1,1,-1
                 xx(1,ii) = yy(1,ii) * invdiag(1,ii)
                 
                 do jj  = ia(ii),ia(ii+1)-1
                    col = ja(jj) 
                    if( col > ii ) then
                       xx(1,ii) = xx(1,ii) - invdiag(1,ii) * an(1,1,jj) * xx(1,col) 
                    end if
                 end do
              end do
           end if
        end if

     else

        if( solve_sol(1) % kfl_renumbered_gs == 1 ) then

           !-------------------------------------------------------------
           !
           ! Streamwise Gauss-Seidel
           !
           !-------------------------------------------------------------
           !
           ! 1. (L+D) y = b
           !
           do ii_new = 1,npoi1
              ii = solve_sol(1) % invpr_gs(ii_new)
              yy(1:nbvar,ii) = bb(1:nbvar,ii)
              do jj  = ia(ii),ia(ii+1)-1
                 col = ja(jj) 
                 if( solve_sol(1) % permr_gs(col) < ii_new ) then
                    do kk = 1,nbvar
                       do ll = 1,nbvar
                          yy(kk,ii) = yy(kk,ii) - an(ll,kk,jj) * yy(ll,col) 
                       end do
                    end do
                 end if
              end do
              yy(1:nbvar,ii) = yy(1:nbvar,ii) * invdiag(1:nbvar,ii)
           end do
           do ii = 1,npoi1
              xx(1:nbvar,ii) = yy(1:nbvar,ii) 
           end do
           goto 10

        else 

           !-------------------------------------------------------------
           !
           ! Gauss-Seidel and symmetric Gauss-Seidel (SSOR)
           !
           !-------------------------------------------------------------
           !
           ! Gauss-Seidel:
           !
           ! 1. (L+D) y = b
           !
           ! Symmetric Gauss-Seidel:
           !
           ! (L+D) D^-1 (U+D) x = b
           ! 1. (L+D) y = b
           ! 2. z = D y
           ! 3. (U+D) x = z
           !
           do ii = 1,npoi1
              yy(1:nbvar,ii) = bb(1:nbvar,ii)
              do jj  = ia(ii),ia(ii+1)-1
                 col = ja(jj) 
                 if( col < ii ) then
                    do kk = 1,nbvar
                       do ll = 1,nbvar
                          yy(kk,ii) = yy(kk,ii) - an(ll,kk,jj) * yy(ll,col) 
                       end do
                    end do
                 end if
              end do
              yy(1:nbvar,ii) = yy(1:nbvar,ii) * invdiag(1:nbvar,ii)
           end do
           if( solve_sol(1) % kfl_renumbered_gs == 0 ) then
              do ii = 1,npoi1
                 xx(1:nbvar,ii) = yy(1:nbvar,ii) 
              end do
              goto 10 
           else           
              !
              ! 2. z = D y
              !              
              do ii = 1,npoi1
                 yy(1:nbvar,ii) = yy(1:nbvar,ii) / invdiag(1:nbvar,ii)
              end do
              !
              ! 3. (U+D) x = z
              !
              do ii = npoi1,1,-1
                 xx(1:nbvar,ii) = yy(1:nbvar,ii) * invdiag(1:nbvar,ii)
                 
                 do jj  = ia(ii),ia(ii+1)-1
                    col = ja(jj) 
                    if( col > ii ) then
                       do kk = 1,nbvar
                          do ll = 1,nbvar
                             xx(kk,ii) = xx(kk,ii) - an(ll,kk,jj) * xx(ll,col) 
                          end do
                       end do
                    end if
                 end do
                 xx(1:nbvar,ii) = xx(1:nbvar,ii) * invdiag(1:nbvar,ii)
              end do
           end if
        end if
        
     end if
     
10   continue
     deallocate(yy)
     
  end if

end subroutine bcsrgs

subroutine bcsrg2(nbnodes,nbvar,an,ja,ia,invdiag,xx,bb) 
  !----------------------------------------------------------------------
  !****f* mathru/bcsrax
  ! NAME 
  !     bcsrax
  ! DESCRIPTION
  !     Multiply a non symmetric matrix stored in BCSR by a vector
  !     WITH ALGEBRAIC LIMITER
  !     XX = A BB 
  ! INPUT
  !    NBNODES .... Number of equations
  !    NBVAR ...... Number of variables
  !    AN ......... Matrix
  !    JA ......... List of elements
  !    IA ......... Pointer to list of elements
  !    BB ......... Vector
  ! OUTPUT
  !    XX ......... result vector
  ! USES
  ! USED BY
  !***
  !----------------------------------------------------------------------
  use def_kintyp, only             :  ip,rp,lg
  use def_master, only             :  INOTMASTER,npoi1
  implicit none
  integer(ip), intent(in)          :: nbnodes,nbvar
  real(rp),    intent(in)          :: an(nbvar,nbvar,*)
  integer(ip), intent(in)          :: ja(*),ia(*)
  real(rp),    intent(in)          :: invdiag(nbvar,*)
  real(rp),    intent(in)          :: bb(nbvar,*)
  real(rp),    intent(inout)       :: xx(nbvar,*)
  integer(ip)                      :: ii,jj,kk,ll,col,nn
  real(rp),    pointer             :: yy(:,:),d2(:,:)

  if( INOTMASTER ) then

     allocate(yy(nbvar,nbnodes))
     allocate(d2(nbvar,nbnodes))

     do ii = 1,nbnodes
        do kk = 1,nbvar
           !d2(kk,ii) = invdiag(kk,ii)
        end do
     end do

     if( nbvar == 1 ) then
        !
        ! NBVAR = 1
        !
        do ii= 1, npoi1 
           jj = ia(ii)
           d2(1,ii) = 0.0_rp
           do while ( jj < ia(ii+1) )
              if( ja(jj) == ii .or. an(1,1,jj) > 0.0_rp ) then
                 d2(1,ii) = d2(1,ii) + an(1,1,jj)
              end if
              jj = jj + 1
           end do
           d2(1,ii) = 1.0_rp / d2(1,ii)
        end do

     else
        !
        ! NBVAR > 1
        !
        do ii= 1, npoi1 
           
           jj = ia(ii)
           ll = -1
           do while ( jj < ia(ii+1) .and. ll == -1 )
              if(ja(jj)==ii) then
                 ll = jj
              end if
              jj = jj + 1
           end do

           jj = ia(ii)
           do while ( jj < ia(ii+1) )
              if( jj == ll ) then
                 do kk = 1,nbvar
                    do nn = 1,nbvar
                       if( an(kk,nn,ll) > 0.0_rp .and. kk /= nn ) then
                          d2(nn,ii) = d2(nn,ii) + an(kk,nn,ll)
                       end if
                    end do
                 end do
              else
                 do kk = 1,nbvar
                    do nn = 1,nbvar
                       if( an(kk,nn,jj) > 0.0_rp ) then
                          d2(nn,ii) = d2(nn,ii) + an(kk,nn,jj)
                       end if
                    end do
                 end do
              end if
              jj = jj + 1
           end do

        end do

     end if
     !
     ! NBVAR = whatever
     !
     !*OMP   PARALLEL DO SCHEDULE (GUIDED)               & 
     !*OMP   DEFAULT (NONE)                              &
     !*OMP   PRIVATE ( ii, jj, kk, ll, col, raux)        &
     !*OMP   SHARED ( nbnodes, nbvar, bb, xx, ia, ja, an)

     !
     ! (L+D) D^-1 (U+D) x = b
     ! 1. (L+D) y = b
     ! 2. z = D y
     ! 3. (U+D) x = z
     !
     do ii = npoi1+1,nbnodes
        do kk = 1,nbvar
           xx(kk,ii) = bb(kk,ii) * d2(kk,ii)
           yy(kk,ii) = xx(kk,ii)
        end do
     end do

     if( nbvar == 1 ) then
        !
        ! 1. (L+D) y = b
        !           
        do ii = 1,npoi1
           !
           ! x_i^{k+1} = 1/a_ii * b_i
           !
           yy(1,ii) = bb(1,ii) * d2(1,ii)

           do jj  = ia(ii),ia(ii+1)-1
              col = ja(jj) 
              if( col < ii .and. an(1,1,jj) < 0.0_rp ) then
                 !
                 ! j<i: x_i^{k+1} = x_i^{k+1} - 1/a_ii * a_ij * x_j^{k+1}
                 !
                 yy(1,ii) = yy(1,ii) - d2(1,ii) * an(1,1,jj) * yy(1,col)
              end if
           end do
        end do
        !
        ! 2. z = D y
        !              
        do ii = 1,npoi1
           yy(1,ii) = yy(1,ii) / d2(1,ii)
        end do
        !
        ! 3. (U+D) x = z
        !
        do ii = npoi1,1,-1
           !
           ! x_i^{k+1} = 1/a_ii * b_i
           !
           xx(1,ii) = yy(1,ii) * d2(1,ii)

           do jj  = ia(ii),ia(ii+1)-1
              col = ja(jj) 
              if( col > ii .and. an(1,1,jj) < 0.0_rp ) then
                 !
                 ! j>i: x_i^{k+1} = x_i^{k+1} - 1/a_ii * a_ij * x_j^k
                 !
                 xx(1,ii) = xx(1,ii) - d2(1,ii) * an(1,1,jj) * xx(1,col) 
              end if
           end do
        end do

     else
        !
        ! 1. (L+D) y = b
        !           
        do ii = 1,npoi1
           !
           ! x_i^{k+1} = 1/a_ii * b_i
           !
           do kk = 1,nbvar
              yy(kk,ii) = bb(kk,ii) * d2(kk,ii)
           end do

           do jj  = ia(ii),ia(ii+1)-1
              col = ja(jj) 
              if( col < ii ) then
                 !
                 ! j<i: x_i^{k+1} = x_i^{k+1} - 1/a_ii * a_ij * x_j^{k+1}
                 !
                 do kk = 1,nbvar
                    do ll = 1,nbvar
                       yy(kk,ii) = yy(kk,ii) - d2(kk,ii) * an(ll,kk,jj) * yy(ll,col) 
                    end do
                 end do
              end if
           end do
        end do
        !
        ! 2. z = D y
        !              
        do ii = 1,npoi1
           do kk = 1,nbvar
              yy(kk,ii) = yy(kk,ii) / d2(kk,ii)
           end do
        end do
        !
        ! 3. (U+D) x = z
        !
        do ii = npoi1,1,-1
           !
           ! x_i^{k+1} = 1/a_ii * b_i
           !           
           do kk = 1,nbvar
              xx(kk,ii) = yy(kk,ii) * d2(kk,ii)
           end do

           do jj  = ia(ii),ia(ii+1)-1
              col = ja(jj) 
              if( col > ii ) then
                 !
                 ! j>i: x_i^{k+1} = x_i^{k+1} - 1/a_ii * a_ij * x_j^k
                 !
                 do kk = 1,nbvar
                    do ll = 1,nbvar
                       xx(kk,ii) = xx(kk,ii) - d2(kk,ii) * an(ll,kk,jj) * xx(ll,col) 
                    end do
                 end do
              end if
           end do
        end do

     end if

     deallocate(d2)
     deallocate(yy)

  end if

end subroutine bcsrg2


