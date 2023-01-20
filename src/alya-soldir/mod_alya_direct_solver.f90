!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Algebraic_Solver
!> @{
!> @name    Alya sparse direct solver
!> @file    mod_alya_direct_solver.f90
!> @author  Guillaume Houzeaux
!> @date    13/06/2016
!> @brief   Alya sparse direct solver module 
!> @details Alya sparse direct solver module 
!>          Input: iA_in, jA_in, An_in
!>          1. Initialization
!>          2. Ordering
!>             Allocate Permr and Invpr if necessary
!>             Compute Invpr, Permr, iA, jA
!>          3. Symbolic factorization:  alya_Symbolical_CSR_LU
!>             Compute iL,jL,iU,jU
!>          +---> 
!>          |   
!>          4. Numerical factorization: alya_Numerical_CSR_LU
!>             Permute An_in into An
!>             Compute Ln,Un
!>          5. Solve An x = b
!>          6. Partial clean
!>             Deallocate LN,UN 
!>              |
!>          <---+
!>          7. Clean
!>             Deallocate iL,jL,iU,jU
!>             Deallocate iA,jA,An,Invpr,Permr
!------------------------------------------------------------------------

module mod_alya_direct_solver

  use def_kintyp_basic, only :  ip,rp,i1p,lg
  use mod_memory_basic, only :  memory_alloca
  use mod_memory_basic, only :  memory_deallo
  use def_mat_csr,      only :  mat_csr
  use def_mat_sky,      only :  mat_sky
  use mod_skyline,      only :  chosol
  use mod_skyline,      only :  chofac
  implicit none
  private
  
  integer(8) :: alya_CSR_memor(2) = 0_8

  type :: alya_direct_solver  
     integer(ip), pointer      :: iL(:)  ! Lower part graph
     integer(ip), pointer      :: jL(:)  ! Lower part graph
     real(rp),    pointer      :: Ln(:)  ! Lower part values
     integer(ip), pointer      :: iU(:)  ! Upper part graph
     integer(ip), pointer      :: jU(:)  ! Upper part graph
     real(rp),    pointer      :: Un(:)  ! Upper part graph
     real(rp)                  :: fillin ! Matrix fillin
     integer(ip)               :: n
     integer(ip)               :: ndof
     integer(ip)               :: sing   ! singularity marker
     integer(8)                :: memor(2)
   contains
     procedure,           pass :: init
     procedure,           pass :: solve
     procedure,           pass :: deallo_symbolical
     procedure,           pass :: deallo_numerical
     procedure,           pass :: factor_symbolical
     procedure,           pass :: factor_numerical
  end type alya_direct_solver
  
  type :: alya_cholesky_solver
     type(mat_sky)             :: a
   contains
     procedure,           pass :: init             => init_cholesky
     procedure,           pass :: deallo           => deallo_cholesky
     procedure,           pass :: solve            => solve_cholesky    
     procedure,           pass :: factor_numerical => factor_numerical_cholesky
  end type alya_cholesky_solver
  
  public :: alya_Symbolical_CSR_LU            ! Symbolical factorization: iL,jL,iU,jU
  public :: alya_Numerical_CSR_LU             ! Numerical factorization; Ln,Un
  public :: alya_CSR_LUSol                    ! Solution: x
  public :: alya_Numerical_CSR_LU_Deallocate  ! Deallocate Un,Ln
  public :: alya_Symbolical_CSR_LU_Deallocate ! Deallocate 
  public :: alya_CSR_LUfin                    ! Deallocate everything

  public :: alya_cholesky_factorization       ! Cholesky factorization
  public :: alya_cholesky_solution            ! Cholesky solution

  public :: alya_direct_solver_initialization ! Initializaiton of the module

  public :: alya_CSR_memor                    ! Memory counter

  public :: alya_direct_solver  
  public :: alya_cholesky_solver  

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-12-11
  !> @brief   Initialization
  !> @details Initialization
  !> 
  !-----------------------------------------------------------------------

  subroutine init(self)
    class(alya_direct_solver), intent(inout) :: self

    self % memor  = 0_8
    self % fillin = 0
    self % n      = 0
    self % ndof   = 1
    self % sing   = 0
    nullify(self % iL)
    nullify(self % jL)
    nullify(self % Ln)
    nullify(self % iU)
    nullify(self % jU)
    nullify(self % Un)
    
  end subroutine init
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-12-11
  !> @brief   Deallocate
  !> @details Deallocate
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo_symbolical(self)
    class(alya_direct_solver), intent(inout) :: self

    call alya_Symbolical_CSR_LU_Deallocate(self % iL,self % jL,self % iU,self % jU,MEMORY_COUNTER=self % memor)
    
  end subroutine deallo_symbolical
  
  subroutine deallo_numerical(self)
    class(alya_direct_solver), intent(inout) :: self

    call alya_Numerical_CSR_LU_Deallocate(self % Ln,self % Un,MEMORY_COUNTER=self % memor)
    
  end subroutine deallo_numerical
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-12-11
  !> @brief   Direct solver initialization
  !> @details Initialize some variables of this module
  !> 
  !-----------------------------------------------------------------------

  subroutine alya_direct_solver_initialization()

    alya_CSR_memor = 0_8

  end subroutine alya_direct_solver_initialization
  
  !---------------------------------------------------------------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-12-11
  !> @brief   Direct solver initialization
  !> @details This routine performs Symbolical CSR LU factorization of a matrix A stored in the CSR format.
  !>          INPUT ARGUMENTS
  !>            NBROWS .... Number of rows of matrix A (dimension of matrix A)
  !>            IA ........ CSR format: index vector for beginning of a row block for matrix A
  !>            JA ........ CSR format: index vector for column numbers for matrix A
  !>          OUTPUT ARGUMENTS 
  !>            IL ........ CSR format: index vector for beginning of a row block for lower triangular matrix L
  !>            JL ........ CSR format: index vector for column numbers for lower triangular matrix L
  !>            IU ........ CSR format: index vector for beginning of a row block for upper triangular matrix U
  !>            JU ........ CSR format: index vector for column numbers for upper triangular matrix U
  !>
  !---------------------------------------------------------------------------------------------------------------------------

  subroutine factor_symbolical(self,a,ia,ja)
    class(alya_direct_solver),          intent(inout) :: self
    class(*),                           intent(in)    :: a
    integer(ip), optional, pointer,     intent(in)    :: ia(:)
    integer(ip), optional, pointer,     intent(in)    :: ja(:)

    select type ( a )
    class is ( mat_csr )
       self % n    = a % nrows
       self % ndof = a % ndof1
       if( present(ia) .and. present(ja) ) then
          call alya_Symbolical_CSR_LU(&
               a % nrows,iA,jA,self % iL,self % jL,self % iU,&
               self % jU,self % fillin,MEMORY_COUNTER=self % memor)
       else
          if( associated(a % ia) .and. associated(a % ja) ) then
             call alya_Symbolical_CSR_LU(&
                  a % nrows,a % iA,a % jA,self % iL,self % jL,self % iU,&
                  self % jU,self % fillin,MEMORY_COUNTER=self % memor)
          end if
       end if
    end select
    
  end subroutine factor_symbolical

  subroutine alya_Symbolical_CSR_LU(nbrows,iA,jA,iL,jL,iU,jU,rfillin,MEMORY_COUNTER)

    implicit none
    integer(ip), intent(in)              :: nbrows
    integer(ip), intent(in)              :: iA(*),jA(*)
    integer(ip), pointer                 :: iL(:),jL(:)
    integer(ip), pointer                 :: iU(:),jU(:)
    real(rp),    intent(out),   optional :: rfillin
    integer(8),  intent(inout), optional :: MEMORY_COUNTER(2)
    integer(ip)                          :: i,j,k,col,nxt,totalL,totalU
    integer(ip)                          :: rowIA_ini,rowIA_diag,rowIA_fin  
    integer(ip)                          :: minColL,maxColL,minColU,maxColU
    integer(ip), pointer                 :: firstL(:),firstU(:),seen(:),nextU(:),iwa(:)   
    integer(ip), pointer                 :: nzL(:),nzU(:)       !Local arrays for L/U matrices 
    type(i1p),   pointer                 :: ptL(:),ptU(:)
    integer(8)                           :: memor(2)

    if( present(MEMORY_COUNTER) ) then
       memor = MEMORY_COUNTER
    else
       memor = alya_CSR_memor
    end if

    if( associated(iL) .or. associated(jL) .or. associated(iU) .or. associated(jU) ) then
       write(*,*) 'Error!! Pointers already associated.'
    else
       !
       ! Alloc local work arrays
       !
       nullify(firstL)
       nullify(firstU)
       nullify(seen)
       nullify(nextU)
       nullify(iwa)
       nullify(nzL)
       nullify(nzU)
       nullify(ptL)
       nullify(ptU)

       call memory_alloca(memor,'FIRSTL','Symbolical_CSR_LU',firstL,nbrows)
       call memory_alloca(memor,'FIRSTU','Symbolical_CSR_LU',firstU,nbrows)
       call memory_alloca(memor,'SEEN'  ,'Symbolical_CSR_LU',seen,nbrows)
       call memory_alloca(memor,'NEXTU' ,'Symbolical_CSR_LU',nextU,nbrows)
       call memory_alloca(memor,'IWA'   ,'Symbolical_CSR_LU',iwa,nbrows)
       call memory_alloca(memor,'NZL'   ,'Symbolical_CSR_LU',nzL,nbrows)
       call memory_alloca(memor,'NZU'   ,'Symbolical_CSR_LU',nzU,nbrows)
       call memory_alloca(memor,'PTL'   ,'Symbolical_CSR_LU',ptL,nbrows)
       call memory_alloca(memor,'PTU'   ,'Symbolical_CSR_LU',ptU,nbrows)
       !
       ! Initialize the fill-in links and local L/U arrays
       !
       do i= 1, nbrows
          firstL(i) = -1
          firstU(i) = -1
          seen(i)   = -1
          nzL(i)    =  0
          nzU(i)    =  0
          nullify( ptL(i) % l , ptU(i) % l )   
       end do
       !
       ! Initialize total number of non-zero elements in L/U
       !
       totalL = 0
       totalU = 0
       !
       ! Main loop in rows
       !
       do i= 1, nbrows
          minColL = nbrows
          maxColL = -1
          minColU = nbrows
          maxColU = -1
          !
          !For all the elements in the i-th row of the matrix A
          !
          rowIA_ini  = iA(i)
          rowIA_fin  = iA(i+1)-1
          rowIA_diag = iA(i+1)
          !
          ! Find a diagonal element in the i-th row of the matrix A
          !
          loop1: do k = rowIA_ini, rowIA_fin
             if( jA(k) == i ) then
                rowIA_diag = k
                exit loop1
             end if
          end do loop1
          !
          ! For all the elements in the i-th row of the matrix A before diagonal
          !
          do j = rowIA_ini,rowIA_diag-1
             col = jA(j)
             if( seen(col) /= i ) then
                seen(col) = i
                minColL   = min(minColL,col)
                maxColL   = max(maxColL,col)
                !
                ! Compute Reachable Set from L
                !
                col = firstL(col)
                loop2:do while (col /= -1)
                   if( seen(col) /= i ) then
                      minColL   = min(minColL,col)
                      maxColL   = max(maxColL,col)
                      seen(col) = i
                      col = firstL(col)
                   else
                      exit loop2
                   end if
                end do loop2
             end if
          end do
          !
          ! For all the elements in the i-th row of the matrix A after diagonal
          !
          do j = rowIA_diag+1, rowIA_fin
             col       = jA(j)
             seen(col) = i
             minColU   = min(minColU,col)
             maxColU   = max(maxColU,col)
          end do
          !
          ! Compute Reachable Set from U
          ! 
          k = firstU(i)
          do while (k /= -1)
             !
             ! For all the elements in the k-th row of the matrix U without the diagonal element
             !
             do j= 2, nzU(k)
                col = ptU(k) % l(j)
                if( col > i ) then
                   minColU = min(minColU,col)
                   maxColU = max(maxColU,col)
                   seen(col) = i
                end if
             end do
             k = nextU(k)
          end do
          !
          ! For all the non-zero elements of the matrix L. L matrix is stored without the diagonal element
          !
          nxt = 0
          do j= minColL, maxColL
             if( seen(j) == i ) then
                iwa(nxt+1) = j
                nxt = nxt + 1
                if( firstL(j) == -1) firstL(j) = i
             end if
          end do
          !
          ! Allocate space for L values and copy iwa to ptL 
          !
          totalL = totalL + nxt
          nzL(i) = nxt
          call memory_alloca(memor,'PTL % L','Symbolical_CSR_LU',ptL(i) % l,nxt,'DO_NOT_INITIALIZE')
          !if( istat /= 0 ) call memerr(0_ip,'ptL(i)%l','Symbolical_CSR_LU',0_ip)
          do k= 1, nxt
             ptL(i) % l(k) = iwa(k)
          end do
          !
          ! For all the non-zero elements of the matrix U. The diagonal element must be included
          !
          iwa(1) = i       !The diagonal element 
          nxt    = 1
          do j= minColU, maxColU
             if( seen(j) == i ) then
                iwa(nxt+1) = j
                nxt = nxt + 1
             end if
          end do
          !
          ! Put the proper link for future fill-in generation 
          !
          if( nxt > 2 ) then
             col         = iwa(2)
             nextU(i)    = firstU(col)
             firstU(col) = i
          end if
          !
          ! Allocate space for U values and copy iwa to ptU 
          !
          totalU = totalU + nxt
          nzU(i) = nxt 
          call memory_alloca(memor,'PTU % L','Symbolical_CSR_LU',ptU(i) % l,nxt,'DO_NOT_INITIALIZE')
          !if( istat /= 0 ) call memerr(0_ip,'ptU(i)%l','Symbolical_CSR_LU',0_ip)
          do k= 1, nxt
             ptU(i) % l(k) = iwa(k)
          end do
       end do
       !
       ! Put L and U in CSR format 
       !
       call memory_alloca(memor,'IL','Symbolical_CSR_LU',iL,nbrows+1_ip)
       call memory_alloca(memor,'IU','Symbolical_CSR_LU',iU,nbrows+1_ip)
       call memory_alloca(memor,'JL','Symbolical_CSR_LU',jL,totalL)
       call memory_alloca(memor,'JU','Symbolical_CSR_LU',jU,totalU)
       !
       ! CSR format for the L matrix
       !
       iL(1) = 1
       do i= 1, nbrows
          iL(i+1) = iL(i) + nzL(i)
          nxt = iL(i) - 1
          do k= 1, nzL(i)
             jL(nxt+k) = ptL(i) % l(k)
          end do
       end do
       !
       ! CSR format for the U matrix
       !
       iU(1) = 1
       do i= 1, nbrows
          iU(i+1) = iU(i) + nzU(i)

          nxt = iU(i) - 1
          do k= 1, nzU(i)
             jU(nxt+k) =  ptU(i) % l(k)
          end do
       end do
       !
       ! Free local arrays
       !
       !do i= 1, nbrows
       !   deallocate(ptL(i) % l) 
       !   if( istat /= 0 ) call memerr(2_ip,'ptL(i)%l','Symbolical_CSR_LU',0_ip)
       !   deallocate(ptU(i) % l)
       !   if( istat /= 0 ) call memerr(2_ip,'ptU(i)%l','Symbolical_CSR_LU',0_ip)
       !end do

       call memory_deallo(memor,'PTU'   ,'Symbolical_CSR_LU',ptU)
       call memory_deallo(memor,'PTL'   ,'Symbolical_CSR_LU',ptL)
       call memory_deallo(memor,'NZU'   ,'Symbolical_CSR_LU',nzU)
       call memory_deallo(memor,'NZL'   ,'Symbolical_CSR_LU',nzL)
       call memory_deallo(memor,'IWA'   ,'Symbolical_CSR_LU',iwa)
       call memory_deallo(memor,'NEXTU' ,'Symbolical_CSR_LU',nextU)
       call memory_deallo(memor,'SEEN'  ,'Symbolical_CSR_LU',seen)
       call memory_deallo(memor,'FIRSTU','Symbolical_CSR_LU',firstU)
       call memory_deallo(memor,'FIRSTL','Symbolical_CSR_LU',firstL)

       if( present(rfillin) ) then
          rfillin = (real(iL(nbrows+1)-1,rp)+real(iU(nbrows+1)-1,rp))/real(iA(nbrows+1)-1,rp)
       end if

       if( present(MEMORY_COUNTER) ) then
          MEMORY_COUNTER = memor
       else
          alya_CSR_memor = memor 
       end if

    end if

  end subroutine Alya_Symbolical_CSR_LU

  subroutine alya_Numerical_CSR_LU1(nbrows,iA,jA,An,iL,jL,Ln,iU,jU,Un,sing,invpr,permr,iAold,jAold)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    Numerical_CSR_LU1
    ! DESCRIPTION
    !    This routine performs Numerical CSR LU factorization of a matrix A stored in the CSR Row (CSR) format.
    !    This routine deals with matrices whose number of degrees of freedom is 1.
    ! INPUT ARGUMENTS
    !    NBROWS .... Number of rows of matrix A (dimension of matrix A)
    !    IA ........ CSR format: index vector for beginning of a row block for matrix A
    !    JA ........ CSR format: index vector for column numbers for matrix A
    !    AN ........ Sparse complex matrix A in BCSR (Blocked CSR Row) format
    !    IL ........ CSR format: index vector for beginning of a row block for lower triangular matrix L
    !    JL ........ CSR format: index vector for column numbers for lower triangular matrix L
    !    IU ........ CSR format: index vector for beginning of a row block for upper triangular matrix U
    !    JU ........ CSR format: index vector for column numbers for upper triangular matrix U
    ! OUTPUT ARGUMENTS 
    !    LN ........ Sparse complex lower triangular matrix L in BCSR (Blocked CSR Row) format
    !    UN ........ Sparse complex upper triangular matrix U in BCSR (Blocked CSR Row) format
    !    SING ...... Singularity marker
    !---------------------------------------------------------------------------------------------------------------------------

    implicit none
    !
    ! Dummy arguments
    !
    integer(ip), intent(in)                     :: nbrows
    integer(ip), intent(in)                     :: iA(*), jA(*)
    real(rp),    intent(in)                     :: An(*)
    integer(ip), intent(in)                     :: iL(*), jL(*)
    real(rp),    intent(inout), pointer         :: Ln(:)
    integer(ip), intent(in)                     :: iU(*), jU(*)
    real(rp),    intent(inout), pointer         :: Un(:)
    integer(ip), intent(out)                    :: sing
    integer(ip), intent(in),  pointer, optional :: invpr(:)
    integer(ip), intent(in),  pointer, optional :: permr(:)
    integer(ip), intent(in),  pointer, optional :: iAold(:)
    integer(ip), intent(in),  pointer, optional :: jAold(:)
    !
    ! Local variables
    !
    integer(ip)                                 :: i,j,k,col,coli,dia
    integer(ip)                                 :: iold,jold
    real(rp)                                    :: pivot
    real(rp),    pointer                        :: wa(:) 
    !
    ! Allocate working array, wa
    !
    nullify(wa)
    call memory_alloca(alya_CSR_memor,'WA','Numerical_CSR_LU1',wa,nbrows,'DO_NOT_INITIALIZE')
    sing = 0
    !
    ! Main loop in rows
    !
    do i = 1, nbrows
       !
       ! Initialize wa with zeros and values of the matrix A 
       !
       do j = iL(i), iL(i+1)-1
          col     = jL(j)
          wa(col) = 0.0_rp
       end do
       do j = iU(i), iU(i+1)-1
          col     = jU(j)
          wa(col) = 0.0_rp
       end do
       if( present(invpr) ) then
          iold = invpr(i) 
          do jold = iAold(iold), iAold(iold+1)-1
             col     = permr(jAold(jold))
             wa(col) = An(jold) 
          end do
       else
          do j = iA(i), iA(i+1)-1
             col     = jA(j)
             wa(col) = An(j)
          end do
       end if
       !
       ! Factorize the i-th row of the matrix A
       ! For all the elements in the i-th row before diagonal
       !
       do j = iL(i), iL(i+1)-1
          col     = jL(j) 
          dia     = iU(col)
          pivot   = wa(col) * Un(dia)       !U diagonal is stored inverted 
          wa(col) = pivot
          !For all the elements in the col-th row after diagonal
          do k = iU(col)+1, iU(col+1)-1
             coli     = jU(k)
             wa(coli) = wa(coli) - pivot * Un(k)
          end do
       end do
       !
       ! Store factorized i-th row in L/U 
       !
       pivot = wa(i)
       if( pivot /= 0.0_rp ) then
          !Store the elements before diagonal in L
          do j     = iL(i), iL(i+1)-1
             col   = jL(j)
             Ln(j) = wa(col)
          end do
          !Store the diagonal element in U, at the first place in the i-th row of the matrix U
          dia     = iU(i) 
          Un(dia) = 1.0_rp / pivot       !Diagonal element is stored inverted
          !Store the elements after diagonal in U
          do j = iU(i)+1, iU(i+1)-1
             col   = jU(j)
             Un(j) = wa(col)
          end do
       else
          write(*,*) 'Singular pivot in row ', i
          sing = -1
       end if
    end do

    !Deallocate working array
    call memory_deallo(alya_CSR_memor,'WA','Numerical_CSR_LU1',wa)

  end subroutine Alya_Numerical_CSR_LU1

  subroutine alya_Numerical_CSR_LUM(nbrows,ndof,iA,jA,An,iL,jL,Ln,iU,jU,Un,sing,invpr,permr,iAold,jAold)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    Numerical_CSR_LU
    ! DESCRIPTION
    !    This routine performs Numerical CSR LU factorization of a matrix A stored in the CSR Row (CSR) format.
    ! INPUT ARGUMENTS
    !    NBROWS .... Number of rows of matrix A (dimension of matrix A)
    !    NDOF ...... Number of degrees of freedom in each node
    !    IA ........ CSR format: index vector for beginning of a row block for matrix A
    !    JA ........ CSR format: index vector for column numbers for matrix A
    !    AN ........ Sparse complex matrix A in BCSR (Blocked CSR Row) format
    !    IL ........ CSR format: index vector for beginning of a row block for lower triangular matrix L
    !    JL ........ CSR format: index vector for column numbers for lower triangular matrix L
    !    IU ........ CSR format: index vector for beginning of a row block for upper triangular matrix U
    !    JU ........ CSR format: index vector for column numbers for upper triangular matrix U
    ! OUTPUT ARGUMENTS 
    !    LN ........ Sparse complex lower triangular matrix L in BCSR (Blocked CSR Row) format
    !    UN ........ Sparse complex upper triangular matrix U in BCSR (Blocked CSR Row) format
    !    SING ...... Singularity marker
    !---------------------------------------------------------------------------------------------------------------------------
    implicit none
    !
    ! Dummy arguments
    !
    integer(ip), intent(in)                     :: nbrows, ndof
    integer(ip), intent(in)                     :: iA(*)
    integer(ip), intent(in)                     :: jA(*)
    real(rp),    intent(in)                     :: An(ndof,ndof,*)
    integer(ip), intent(in)                     :: iL(*)
    integer(ip), intent(in)                     :: jL(*)
    real(rp),    intent(out)                    :: Ln(ndof,ndof,*)
    integer(ip), intent(in)                     :: iU(*)
    integer(ip), intent(in)                     :: jU(*)
    real(rp),    intent(out)                    :: Un(ndof,ndof,*)
    integer(ip), intent(out)                    :: sing
    integer(ip), intent(in),  pointer, optional :: invpr(:)
    integer(ip), intent(in),  pointer, optional :: permr(:)
    integer(ip), intent(in),  pointer, optional :: iAold(:)
    integer(ip), intent(in),  pointer, optional :: jAold(:)
    !
    ! Local variables
    !
    integer(ip)                                 :: i,j,k,l,s,t,col,dia,coli
    integer(ip)                                 :: iold,jold
    real(rp)                                    :: pivot
    real(rp), pointer                           :: wa(:,:,:)     
    !
    ! Allocate working array, wa
    !
    nullify(wa)
    call memory_alloca(alya_CSR_memor,'WA','Numerical_CSR_LUM',wa,ndof,ndof,nbrows,'DO_NOT_INITIALIZE')

    sing = 0
    !
    ! Main loop in rows
    !
    do i= 1, nbrows
       !
       ! Initialize wa with zeros and values of the matrix A 
       !
       do j= iL(i), iL(i+1)-1
          col = jL(j)
          do l= 1, ndof
             do k= 1, ndof
                wa(k,l,col) = 0.0_rp
             end do
          end do
       end do
       do j= iU(i), iU(i+1)-1
          col = jU(j)
          do l= 1, ndof
             do k= 1, ndof
                wa(k,l,col) = 0.0_rp
             end do
          end do
       end do
       if( present(invpr) ) then
          iold = invpr(i) 
          do jold = iAold(iold), iAold(iold+1)-1
             col = permr(jAold(jold))
             do l= 1, ndof
                do k= 1, ndof
                   wa(k,l,col) = An(l,k,jold)
                end do
             end do
          end do
       else
          do j= iA(i), iA(i+1)-1
             col = jA(j)
             do l= 1, ndof
                do k= 1, ndof
                   wa(k,l,col) = An(l,k,j)
                end do
             end do
          end do
       end if
       !
       ! Factorize all the blocks in the i-th row of the matrix A
       ! For all the elements (blocks) in the i-th row before diagonal
       !
       do j= iL(i), iL(i+1)-1
          col = jL(j)
          dia = iU(col)
          !Factorize the block in column col
          do l= 1, ndof
             do k= 1, ndof      
                pivot = wa(k,l,col) * Un(l,l,dia)       !U diagonal is stored inverted
                wa(k,l,col) = pivot
                !For all the elements in the k-th row inside the current block
                do s= l+1, ndof
                   wa(k,s,col) = wa(k,s,col) - pivot * Un(l,s,dia)
                end do
             end do
          end do
          !For all the blocks after diagonal
          do s= iU(col)+1, iU(col+1)-1
             coli = jU(s)
             do l= 1, ndof
                do k= 1, ndof
                   do t= 1, ndof 
                      wa(k,t,coli) = wa(k,t,coli) - Un(l,t,s) * wa(k,l,col)
                   end do
                end do
             end do
          end do
       end do
       !
       ! Factorize the lower triangle of the diagonal block
       !
       do k= 1, ndof
          do l= 1, k-1
             pivot = wa(k,l,i) * wa(l,l,i)
             wa(k,l,i) = pivot
             !For all the elements in the k-th row inside the diagonal block
             do s= l+1, ndof
                wa(k,s,i) = wa(k,s,i) - pivot * wa(l,s,i)
             end do
          end do
          !Check if some of the diagonal elements of the diagonal block is zero
          if( wa(k,k,i) /= 0.0_rp ) then
             wa(k,k,i) = 1.0_rp / wa(k,k,i)
          else
             sing = -1
             print*,'singular matrix'
             exit
          end if
       end do
       !
       ! For all the blocks after diagonal
       !
       do s= iU(i)+1, iU(i+1)-1
          coli = jU(s)
          do k= 1, ndof
             do l= 1, k-1
                do t= 1, ndof 
                   wa(k,t,coli) = wa(k,t,coli) - wa(l,t,coli) * wa(k,l,i)
                end do
             end do
          end do
       end do
       !
       ! Store factorized i-th row of blocks in L/U matrices
       !
       if( sing == 0 ) then
          !Store the elements of the blocks before diagonal in L
          do j= iL(i), iL(i+1)-1
             col= jL(j)
             do l= 1, ndof
                do k= 1, ndof   
                   Ln(k,l,j) = wa(k,l,col)
                end do
             end do
          end do
          !Store the elements of the diagonal block in U
          !Note that the WHOLE diagonal block is stored in the matrix U
          col = iU(i)
          do l= 1, ndof
             do k= 1, ndof   
                Un(k,l,col) = wa(k,l,i)
             end do
          end do
          !Store the elements of the blocks after diagonal in U
          do j= iU(i)+1, iU(i+1)-1
             col= jU(j)
             do l= 1, ndof
                do k= 1, ndof
                   Un(k,l,j) = wa(k,l,col)    
                end do
             end do
          end do
       end if
    end do
    !
    ! Deallocate working array   
    !
    call memory_deallo(alya_CSR_memor,'WA','Numerical_CSR_LUM',wa)

  end subroutine Alya_Numerical_CSR_LUM
  
  !---------------------------------------------------------------------------------------------------------------------------
  ! NAME
  !    CSR_LU_Factorization
  ! DESCRIPTION
  !    This routine performs CSR LU factorization of a matrix A stored in the CSR Row (CSR) format.
  ! INPUT ARGUMENTS
  !    NBROWS .... Number of rows of matrix A (dimension of matrix A)
  !    NDOF ...... Number of degrees of freedom in each node
  !    IA ........ CSR format: index vector for beginning of a row block for matrix A
  !    JA ........ CSR format: index vector for column numbers for matrix A
  !    AN ........ Sparse complex matrix A in BCSR (Blocked CSR Row) format
  ! OUTPUT ARGUMENTS 
  !    IL ........ CSR format: index vector for beginning of a row block for lower triangular matrix L
  !    JL ........ CSR format: index vector for column numbers for lower triangular matrix L
  !    LN ........ Sparse complex lower triangular matrix L in BCSR (Blocked CSR Row) format
  !    IU ........ CSR format: index vector for beginning of a row block for upper triangular matrix U
  !    JU ........ CSR format: index vector for column numbers for upper triangular matrix U
  !    UN ........ Sparse complex upper triangular matrix U in BCSR (Blocked CSR Row) format
  !    SING ...... Singularity marker
  !---------------------------------------------------------------------------------------------------------------------------
  
  subroutine factor_numerical(self,a,invp,perm,ia,ja)
    class(alya_direct_solver),                   intent(inout) :: self 
    class(*),                                    intent(in)    :: a
    integer(ip),              optional, pointer, intent(in)    :: invp(:)
    integer(ip),              optional, pointer, intent(in)    :: perm(:)
    integer(ip),              optional, pointer, intent(in)    :: ia(:)
    integer(ip),              optional, pointer, intent(in)    :: ja(:)

    select type ( a )
    class is ( mat_csr )
       if( associated(a % ia) .and. associated(a % ja) ) then
          if( present(invp) .and. present(perm) .and. present(ia) .and. present(ja) ) then
             call alya_Numerical_CSR_LU(a % nrows,a % ndof1,iA,jA,a % vA,&
                  self % iL,self % jL,self % Ln,self % iU,self % jU,self % Un,self % sing,&
                  invp,perm,a % ia,a % ja)
          else
             call alya_Numerical_CSR_LU(a % nrows,a % ndof1,a % iA,a % jA,a % vA,&
                  self % iL,self % jL,self % Ln,self % iU,self % jU,self % Un,self % sing)
          end if
       end if
    end select
    
  end subroutine factor_numerical
  
  subroutine alya_Numerical_CSR_LU(nbrows,ndof,iA,jA,An,iL,jL,Ln,iU,jU,Un,sing,invpr,permr,iAold,jAold)

    implicit none
    !
    ! Dummy arguments
    !
    integer(ip), intent(in)                    :: nbrows  !< Number of rows of matrix A (dimension of matrix A)
    integer(ip), intent(in)                    :: ndof    !< Number of degrees of freedom in each node
    integer(ip), intent(in)                    :: iA(*)   !< CSR format: index vector for beginning of a row block for matrix A
    integer(ip), intent(in)                    :: jA(*)   !< CSR format: index vector for column numbers for matrix A
    real(rp),    intent(in)                    :: An(*)   !< Sparse complex matrix A in BCSR (Blocked CSR Row) format
    integer(ip), pointer                       :: iL(:)   !< CSR format: index vector for beginning of a row block for lower triangular matrix L
    integer(ip), pointer                       :: jL(:)   !< CSR format: index vector for column numbers for lower triangular matrix L
    real(rp),    pointer                       :: Ln(:)   !< Sparse complex lower triangular matrix L in BCSR (Blocked CSR Row) format
    integer(ip), pointer                       :: iU(:)   !< CSR format: index vector for beginning of a row block for upper triangular matrix U
    integer(ip), pointer                       :: jU(:)   !< CSR format: index vector for column numbers for upper triangular matrix U
    real(rp),    pointer                       :: Un(:)   !< Sparse complex upper triangular matrix U in BCSR (Blocked CSR Row) format
    integer(ip), intent(out)                   :: sing    !< Singularity marker
    integer(ip), intent(in), optional, pointer :: invpr(:)
    integer(ip), intent(in), optional, pointer :: permr(:)
    integer(ip), intent(in), optional, pointer :: iAold(:)
    integer(ip), intent(in), optional, pointer :: jAold(:)
    logical(lg)                                :: lpermute
    !
    ! Allocate Un and Ln if necessary
    !
    call alya_Numerical_CSR_LU_Allocate(nbrows,ndof,iL,iU,Ln,Un)
    !
    ! Permute or not
    !
    lpermute = present(invpr) .and. present(permr) .and. present(iAold) .and. present(jAold)
    !
    ! Perform Numerical factorization dependig on number of degrees of freedom
    !
    if( ndof == 1 ) then
       if( lpermute ) then
          call alya_Numerical_CSR_LU1(nbrows,iA,jA,An,iL,jL,Ln,iU,jU,Un,sing,invpr,permr,iAold,jAold)
       else
          call alya_Numerical_CSR_LU1(nbrows,iA,jA,An,iL,jL,Ln,iU,jU,Un,sing)
       end if
    else
       if( lpermute ) then 
          call alya_Numerical_CSR_LUM(nbrows,ndof,iA,jA,An,iL,jL,Ln,iU,jU,Un,sing,invpr,permr,iAold,jAold)  
       else
          call alya_Numerical_CSR_LUM(nbrows,ndof,iA,jA,An,iL,jL,Ln,iU,jU,Un,sing)  
       end if
    end if

  end subroutine Alya_Numerical_CSR_LU

  subroutine alya_CSR_Lsol(nbrows,ndof,invpR,iL,jL,Ln,iU,jU,Un,b,x)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    CSR_Lsol
    ! DESCRIPTION
    !    This routine solves a triangular linear system of equations where the system matrix is L - lower triangular matrix
    ! INPUT ARGUMENTS
    !    NBROWS .... Dimension of the sysytem
    !    NDOF ...... Number of degrees of freedom in each node
    !    INVPR ..... Permutation array: OLD = INVPR(NEW)
    !    IL ........ CSR format: index vector for beginning of a row block for lower triangular matrix L
    !    JL ........ CSR format: index vector for column numbers for lower triangular matrix L
    !    LN ........ Sparse complex lower triangular matrix L in BCSR (Blocked CSR Row) format
    !    IU ........ CSR format: index vector for beginning of a row block for upper triangular matrix U
    !    JU ........ CSR format: index vector for column numbers for upper triangular matrix U
    !    UN ........ Sparse complex upper triangular matrix U in BCSR (Blocked CSR Row) format
    !     B ........ RHS of the system
    ! OUTPUT ARGUMENTS 
    !     X ........ Solution vector of the system
    !---------------------------------------------------------------------------------------------------------------------------

    !
    ! Dummy arguments
    !
    integer(ip), intent(in)   :: nbrows,ndof
    integer(ip), pointer      :: invpR(:)
    integer(ip), intent(in)   :: iL(*),jL(*)
    real(rp),    intent(in)   :: Ln(ndof,ndof,*)
    integer(ip), intent(in)   :: iU(*),jU(*)
    real(rp),    intent(in)   :: Un(ndof,ndof,*)
    real(rp),    intent(in)   :: b(ndof,*)
    real(rp),    intent(out)  :: x(ndof,*)
    !
    ! Local variables
    !
    integer(ip)               :: i,j,k,l,col

    if( .not. associated(invpR) ) then
       !
       ! x <= b
       !
       do j = 1,nbrows
          x(1:ndof,j) = b(1:ndof,j)
       end do
    else 
       !
       ! x (new) <= b (old)
       !
       if( ndof == 1 ) then
          do j = 1,nbrows
             k = invpR(j)
             x(1,j) = b(1,k)
          end do
       else
          do j = 1,nbrows
             k = invpR(j)
             x(1:ndof,j) = b(1:ndof,k)
          end do
       end if
    end if

    if( ndof == 1 ) then
       do i = 1, nbrows
          do j = iL(i),iL(i+1)-1
             col = jL(j)
             x(1,i) = x(1,i) - Ln(1,1,j) * x(1,col)
          end do
       end do
    else
       do i = 1, nbrows
          do j = iL(i),iL(i+1)-1
             col = jL(j)
             do l = 1,ndof  
                do k = 1,ndof
                   x(k,i) = x(k,i) - Ln(k,l,j) * x(l,col)
                end do
             end do
          end do
          !Update using the lower triangle of the diagonal block stored in the U matrix
          j = iU(i)
          do l = 1,ndof
             do k = l+1,ndof
                x(k,i) = x(k,i) - Un(k,l,j) * x(l,i)
             end do
          end do
       end do
    end if

  end subroutine Alya_CSR_Lsol

  subroutine Alya_CSR_Usol(nbrows,ndof,invpC,iU,jU,Un,b,x)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    CSR_Usol
    ! DESCRIPTION
    !    This routine solves a triangular linear system of equations where the system matrix is U - upper triangular matrix
    ! INPUT ARGUMENTS
    !    NBROWS .... Dimension of the sysytem
    !    NDOF ...... Number of degrees of freedom in each node
    !    INVPC ..... 
    !    IU ........ CSR format: index vector for beginning of a row block for lower triangular matrix L
    !    JU ........ CSR format: index vector for column numbers for lower triangular matrix L
    !    UN ........ Sparse complex lower triangular matrix L in BCSR (Blocked CSR Row) format
    !     B ........ RHS of the system
    ! OUTPUT ARGUMENTS 
    !     X ........ Solution vector of the system
    !---------------------------------------------------------------------------------------------------------------------------

    !
    ! Dummy arguments
    !
    integer(ip), intent(in)   :: nbrows, ndof
    integer(ip), pointer      :: invpC(:)
    integer(ip), intent(in)   :: iU(*), jU(*)
    real(rp),    intent(in)   :: Un(ndof,ndof,*)
    real(rp),    intent(inout):: b(ndof,*)
    real(rp),    intent(out)  :: x(ndof,*)
    !
    ! Local variables
    !
    integer(ip)               :: i,j,k,l,col    

    if( ndof == 1 ) then
       DO i= nbrows, 1, -1
          !Update rhs with the block after diagonal block
          do j= iU(i)+1, iU(i+1)-1
             col = jU(j)
             b(1,i) = b(1,i) - Un(1,1,j) * b(1,col)
          end do
          !Solve the upper triangle in the diagonal of the U matrix
          j = iU(i)
          b(1,i) = b(1,i) * Un(1,1,j)
       END DO
    else
       DO i= nbrows, 1, -1
          !Update rhs with the block after diagonal block
          do j= iU(i)+1, iU(i+1)-1
             col = jU(j)
             do l= 1, ndof
                do k= 1, ndof
                   b(k,i) = b(k,i) - Un(k,l,j) * b(l,col)
                end do
             end do
          end do
          !Solve the upper triangle in the diagonal of the U matrix
          j = iU(i)
          do l= ndof, 1, -1
             b(l,i) = b(l,i) * Un(l,l,j)
             do k= 1, l-1
                b(k,i) = b(k,i) - Un(k,l,j) * b(l,i)
             end do
          end do
       END DO
    end if

    if( .not. associated(invpC) ) then
       !
       ! x <= b
       !
       do j= 1, nbrows
          do i= 1, ndof
             x(i,j) = b(i,j)
          end do
       end do
    else
       !
       ! x (old) <= b (new)
       !
       if( ndof == 1 ) then
          do j= 1, nbrows
             col = invpC(j)
             x(1,col) = b(1,j)
          end do
       else
          do j= 1, nbrows
             col = invpC(j)
             do i= 1, ndof
                x(i,col) = b(i,j)
             end do
          end do
       end if
    end if

  end subroutine Alya_CSR_Usol

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    CSR_LUsol
    ! DESCRIPTION
    !    This routine solves LU triangular linear systems of equations.
    !
    !          +---------+
    !    b --> | LUx = b | --> x
    !          +---------+
    !   (old)    (new)        (old)
    !
    ! INPUT ARGUMENTS
    !    NBROWS .... Dimension of the sysytem
    !    NDOF ...... Number of degrees of freedom in each node
    !    INVPR ..... 
    !    INVPC ..... 
    !    IU ........ CSR format: index vector for beginning of a row block for lower triangular matrix L
    !    JU ........ CSR format: index vector for column numbers for lower triangular matrix L
    !    UN ........ Sparse complex lower triangular matrix L in BCSR (Blocked CSR Row) format
    !     B ........ RHS of the system
    ! OUTPUT ARGUMENTS 
    !     X ........ Solution vector of the system
    !---------------------------------------------------------------------------------------------------------------------------

  subroutine solve(self,b,x,a,invp)
    class(alya_direct_solver),          intent(inout) :: self 
    real(rp),                  pointer, intent(in)    :: b(:)              !< RHS b
    real(rp),                  pointer, intent(inout) :: x(:)              !< Solve x
    class(*),                           intent(in)    :: a
    integer(ip),   optional,   pointer, intent(in)    :: invp(:)
    integer(ip),               pointer                :: inul(:)
    
    if( associated(x) .and. associated(b) .and. self % n > 0 ) then
       if( present(invp) ) then
          call alya_CSR_LUsol(&
               self % n,self % ndof,invp,invp,self % iL,&
               self % jL,self % Ln,self % iU,self % jU,self % Un,b,x)
       else
          nullify(inul)
          call alya_CSR_LUsol(&
               self % n,self % ndof,inul,inul,self % iL,&
               self % jL,self % Ln,self % iU,self % jU,self % Un,b,x)
       end if
    end if
    
  end subroutine solve
  
  subroutine alya_CSR_LUsol(nbrows,ndof,invpR,invpC,iL,jL,Ln,iU,jU,Un,b,x)

    !
    ! Dummy arguments
    !
    integer(ip), intent(in)  :: nbrows,ndof
    integer(ip), pointer     :: invpR(:),invpC(:)
    integer(ip), intent(in)  :: iL(*),jL(*)
    real(rp),    intent(in)  :: Ln(ndof,ndof,*)
    integer(ip), intent(in)  :: iU(*),jU(*)
    real(rp),    intent(in)  :: Un(ndof,ndof,*)
    real(rp),    intent(in)  :: b(ndof,*)
    real(rp),    intent(out) :: x(ndof,*)
    real(rp),    pointer     :: b_tmp(:,:)
    !
    ! Local variables
    !
    integer(ip)               :: j

    nullify(b_tmp)
    call memory_alloca(alya_CSR_memor,'B_TMP','CSR_LUsol',b_tmp,ndof,nbrows,'DO_NOT_INITIALIZE')

    b_tmp(1:ndof,1:nbrows) = b(1:ndof,1:nbrows)

    call alya_CSR_Lsol(nbrows,ndof,invpR,iL,jL,Ln,iU,jU,Un,b_tmp,x)
    call alya_CSR_Usol(nbrows,ndof,invpC,iU,jU,Un,x,b_tmp)

    do j = 1,nbrows
       x(1:ndof,j) = b_tmp(1:ndof,j)
    end do

    call memory_deallo(alya_CSR_memor,'B_TMP','CSR_LUsol',b_tmp)

  end subroutine Alya_CSR_LUsol

  subroutine alya_CSR_LUfin(iL,jL,Ln,iU,jU,Un,invpR,invpC)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    CSR_LUfin
    ! DESCRIPTION
    !    This routine solves LU triangular linear systems of equations.eallocates memoryd
    ! INPUT ARGUMENTS
    ! OUTPUT ARGUMENTS 
    !     X ........ Solution vector of the system
    !---------------------------------------------------------------------------------------------------------------------------

    !Dummy arguments
    integer(ip), pointer           :: iL(:),jL(:)
    real(rp),    pointer           :: Ln(:)
    integer(ip), pointer           :: iU(:),jU(:)
    real(rp),    pointer           :: Un(:)
    integer(ip), pointer, optional :: invpR(:)
    integer(ip), pointer, optional :: invpC(:)

    call memory_deallo(alya_CSR_memor,'IL','CSR_LUfin',iL)
    call memory_deallo(alya_CSR_memor,'JL','CSR_LUfin',jL)
    call memory_deallo(alya_CSR_memor,'LN','CSR_LUfin',Ln)
    call memory_deallo(alya_CSR_memor,'IU','CSR_LUfin',iU)
    call memory_deallo(alya_CSR_memor,'JU','CSR_LUfin',jU)
    call memory_deallo(alya_CSR_memor,'UN','CSR_LUfin',Un)
    if( present(invpR) ) call memory_deallo(alya_CSR_memor,'invpR','CSR_LUfin',invpR)
    if( present(invpC) ) call memory_deallo(alya_CSR_memor,'invpC','CSR_LUfin',invpC)

  end subroutine Alya_CSR_LUfin

  subroutine alya_Symbolical_CSR_LU_Deallocate(iL,jL,iU,jU,MEMORY_COUNTER)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    CSR_LUfin
    ! DESCRIPTION
    !    This routine solves LU triangular linear systems of equations.eallocates memoryd
    ! INPUT ARGUMENTS
    ! OUTPUT ARGUMENTS 
    !     X ........ Solution vector of the system
    !---------------------------------------------------------------------------------------------------------------------------

    !Dummy arguments
    integer(ip),           pointer       :: iL(:),jL(:)
    integer(ip),           pointer       :: iU(:),jU(:)
    integer(8),  optional, intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                           :: memor(2)

    if( present(MEMORY_COUNTER) ) then
       memor = MEMORY_COUNTER
    else
       memor = alya_CSR_memor
    end if

    call memory_deallo(memor,'iL','CSR_LUfin',iL)
    call memory_deallo(memor,'jL','CSR_LUfin',jL)
    call memory_deallo(memor,'iU','CSR_LUfin',iU)
    call memory_deallo(memor,'jU','CSR_LUfin',jU)

    if( present(MEMORY_COUNTER) ) then
       MEMORY_COUNTER = memor
    else
       alya_CSR_memor = memor 
    end if

  end subroutine Alya_Symbolical_CSR_LU_Deallocate

  subroutine alya_Numerical_CSR_LU_Deallocate(Ln,Un,MEMORY_COUNTER)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    CSR_LUfin
    ! DESCRIPTION
    !    This routine solves LU triangular linear systems of equations.eallocates memoryd
    ! INPUT ARGUMENTS
    ! OUTPUT ARGUMENTS 
    !     X ........ Solution vector of the system
    !---------------------------------------------------------------------------------------------------------------------------

    !Dummy arguments
    real(rp), pointer :: Ln(:),Un(:)
    integer(8),  optional, intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                           :: memor(2)

    if( present(MEMORY_COUNTER) ) then
       memor = MEMORY_COUNTER
    else
       memor = alya_CSR_memor
    end if

    call memory_deallo(memor,'LN','CSR_LUfin',Ln)
    call memory_deallo(memor,'UN','CSR_LUfin',Un)

    if( present(MEMORY_COUNTER) ) then
       MEMORY_COUNTER = memor
    else
       alya_CSR_memor = memor 
    end if

  end subroutine Alya_Numerical_CSR_LU_Deallocate

  subroutine alya_Numerical_CSR_LU_Allocate(nbrows,ndof,iL,iU,Ln,Un)

    integer(ip), intent(in)             :: nbrows
    integer(ip), intent(in)             :: ndof
    integer(ip), intent(in),    pointer :: iL(:),iU(:)
    real(rp),    intent(inout), pointer :: Ln(:),Un(:)
    integer(ip)                         :: nnzL,nnzU,ndof2
    !
    ! Allocate memory for L/U
    !
    ndof2 = ndof*ndof
    if( .not. associated(Ln) ) then
       nnzL = iL(nbrows+1)-1
       call memory_alloca(alya_CSR_memor,'LN','alya_Numerical_CSR_LU_Allocate',Ln,nnzL*ndof2)
    end if
    if( .not. associated(Un) ) then
       nnzU = iU(nbrows+1)-1
       call memory_alloca(alya_CSR_memor,'UN','alya_Numerical_CSR_LU_Allocate',Un,nnzU*ndof2) 
    end if

  end subroutine alya_Numerical_CSR_LU_Allocate

  subroutine alya_CSR_Permute_and_Copy_matrix(nbrows,ndof,iAin,jAin,iAout,jAout,invpR,Ain,Aout)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    CSR_SMVP
    ! DESCRIPTION
    !    This routine copy a matrix to another
    ! INPUT ARGUMENTS
    !    NBROWS .... Dimension of the sysytem
    !    NDOF ...... Number of degrees of freedom in each node
    !    INVPR ..... Permutation array: OLD = INVPR(NEW)
    !     B ........ RHS of the system
    ! OUTPUT ARGUMENTS 
    !     X ........ Solution vector of the system
    !---------------------------------------------------------------------------------------------------------------------------

    !
    ! Dummy arguments
    !
    integer(ip), intent(in)   :: nbrows            !< Number of rows
    integer(ip), intent(in)   :: ndof              !< Number of dof per row
    integer(ip), intent(in)   :: iAin(*)           !< Matrix graph
    integer(ip), intent(in)   :: jAin(*)           !< Matrix graph
    integer(ip), intent(in)   :: iAout(*)          !< Matrix graph
    integer(ip), intent(in)   :: jAout(*)          !< Matrix graph
    integer(ip), pointer      :: invpR(:)          !< Inverse permutation
    real(rp),    intent(in)   :: Ain(ndof,ndof,*)  !< Input matrix
    real(rp),    intent(out)  :: Aout(ndof,ndof,*) !< Output matrix
    integer(ip)               :: ii,jj,iz
    integer(ip)               :: jjold,iiold,kz

    if( .not. associated(invpR) ) then
       do ii = 1,nbrows
          do iz = iAin(ii),iAin(ii+1)-1
             Aout(1:ndof,1:ndof,iz) = Ain(1:ndof,1:ndof,iz)
          end do
       end do
    else
       do ii = 1,nbrows
          iiold = invpR(ii)
          do iz = iAout(ii),iAout(ii+1)-1
             jj    = jAout(iz)
             jjold = invpR(jj)
             kz    = iAin(iiold)
             do while( jAin(kz) /= jjold )
                kz = kz + 1
             end do
             Aout(1:ndof,1:ndof,iz) = Ain(1:ndof,1:ndof,kz)
          end do
       end do
    end if

  end subroutine Alya_CSR_Permute_and_Copy_matrix

  !-----------------------------------------------------------------------
  !****f* domain/skygro
  ! NAME
  !    skygro
  ! DESCRIPTION
  !    Set up the skyline structure
  ! USED BY
  !    domain
  !***
  !-----------------------------------------------------------------------

  subroutine alya_cholesky_initialization(nbrows,ndof,kfl_symmetric,ia,ja,nskyl,iskyl,idiag)

    integer(ip), intent(in)             :: nbrows
    integer(ip), intent(in)             :: ndof
    integer(ip), intent(in)             :: kfl_symmetric
    integer(ip), intent(out)            :: nskyl
    integer(ip), intent(in),    pointer :: ia(:)
    integer(ip), intent(in),    pointer :: ja(:)
    integer(ip), intent(inout), pointer :: iskyl(:)
    integer(ip), intent(inout), pointer :: idiag(:)
    integer(ip)                         :: kskyl,idof,jdof
    integer(ip)                         :: izdom,ii,jj,kk,ll
    !
    ! Allocate memory
    !
    if( .not. associated(iskyl) ) then
       call memory_alloca(alya_CSR_memor,'ISKYL','alya_cholesky_initialization',iskyl,ndof*nbrows+1)
    end if
    !
    ! Skyline format
    !
    do ii = 1,nbrows*ndof+1
       iskyl(ii) = nbrows*ndof
    end do
    do ii = 1,nbrows
       do izdom = ia(ii),ia(ii+1)-1
          jj = ja(izdom)  
          if( jj > 0 .and. ii >= jj ) then
             do idof = 1,ndof 
                kk = (ii-1)*ndof+idof+1  
                do jdof = 1,ndof
                   ll = (jj-1)*ndof+jdof
                   if( ll < iskyl(kk) ) iskyl(kk) = ll
                end do
             end do
          end if
       end do
    end do

    !call PAR_MIN(ndof*nbrows+1_ip,iskyl)

    nskyl    = 1
    iskyl(1) = 1

    if( kfl_symmetric == 1 ) then
       ! 
       ! For the symmetric case, do not need idiag
       !
       do kk = 1,nbrows*ndof
          kskyl       = kk - iskyl(kk+1) + 1
          nskyl       = nskyl + kskyl
          iskyl(kk+1) = nskyl
       end do

    else
       !
       ! For the nonsymmetric case, set idiag 
       !
       if( .not. associated(idiag) ) then
          call memory_alloca(alya_CSR_memor,'IDIAG','alya_cholesky_initialization',idiag,ndof*nbrows)
       end if
       do kk = 1,nbrows*ndof
          kskyl       = kk - iskyl(kk+1)
          idiag(kk)   = nskyl + kskyl
          kskyl       = 2 * kskyl + 1  
          nskyl       = nskyl + kskyl
          iskyl(kk+1) = nskyl
       end do

    end if
    nskyl = nskyl - 1 

  end subroutine alya_cholesky_initialization
  
  subroutine init_cholesky(self)
    class(alya_cholesky_solver), intent(inout) :: self 

    call self % a % init()
    
  end subroutine init_cholesky
  
  subroutine deallo_cholesky(self)
    class(alya_cholesky_solver), intent(inout) :: self 

    call self % a % deallo()
    
  end subroutine deallo_cholesky
   
  subroutine factor_numerical_cholesky(self,a)
    class(alya_cholesky_solver),                 intent(inout) :: self 
    class(*),                                    intent(in)    :: a
    integer(ip)                                                :: nn,info,iz
    
    select type ( a )
    class is ( mat_sky )
       nn = a % nrows * a % ndof1
       self % a % nskyl = a % nskyl
       self % a % nrows = a % nrows
       self % a % ndof1 = a % ndof1
       self % a % ndof2 = a % ndof2
       call self % a % alloca_matrix()
       do iz = 1,a % nskyl
          self % a % va(iz) = a % va(iz)
       end do
       call chofac( nn, a % nskyl, a % iskyl, self % a % va, info )
    end select
    
  end subroutine factor_numerical_cholesky

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-10-17
  !> @brief   Cholesky factorization
  !> @details Cholesky factorization
  !> 
  !-----------------------------------------------------------------------

  subroutine alya_cholesky_factorization( n, nnz, nz, a, info )

    integer(ip), intent(in)    :: n
    integer(ip), intent(in)    :: nnz
    integer(ip), intent(in)    :: nz( n+1 )
    real(rp),    intent(inout) :: a( nnz )
    integer(ip), intent(out)   :: info
    integer(ip)                :: i, j,  k
    integer(ip)                :: i0, j0
    integer(ip)                :: ipos,idif,irow,jcol
    real(rp)                   :: temp, rdiag

    call chofac( n, nnz, nz, a, info )
  
  end subroutine alya_cholesky_factorization

  subroutine solve_cholesky(self,b,x,a)
    class(alya_cholesky_solver),        intent(inout) :: self 
    real(rp),                  pointer, intent(in)    :: b(:)              !< RHS b
    real(rp),                  pointer, intent(inout) :: x(:)              !< Solve x
    class(*),                           intent(in)    :: a
    integer(ip)                                       :: nn,info,i
    
    select type ( a )
    class is ( mat_sky )
       nn = self % a % nrows * self % a % ndof1
       do i = 1,nn
          x(i) = b(i)
       end do
       call chosol( nn, a % nskyl, a % iskyl, 1_ip, self % a % va, x,nn, info )
    end select
    
  end subroutine solve_cholesky
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-10-17
  !> @brief   Cholesky solution
  !> @details Cholesky solution
  !> 
  !-----------------------------------------------------------------------

  subroutine alya_cholesky_solution( n, nnz, nz, nrhs, a, b, ldb, info )

    integer(ip)          :: info, ldb, n, nnz, nrhs
    integer(ip)          :: nz( n+1 )
    real(rp)             :: a( nnz ), b( ldb, * )
    integer(ip)          :: j, k
    integer(ip)          :: k0
    integer(ip)          :: jcnt, i
    real(rp)             :: temp

    
    call chosol( n, nnz, nz, nrhs, a, b, ldb, info )
    
  end subroutine alya_cholesky_solution

end module mod_alya_direct_solver
!> @} 
