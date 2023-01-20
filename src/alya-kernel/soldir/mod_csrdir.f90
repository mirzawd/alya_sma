!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_csrdir

  !------------------------------------------------------------------------
  !****f* soldir/mod_csrdir
  ! NAME 
  !    mod_csrdir
  ! DESCRIPTION
  !    Sparse direct solver module
  ! USES
  ! USED BY 
  !***
  !------------------------------------------------------------------------

  use def_kintyp, only : ip,rp,i1p
  use mod_memory, only : memory_alloca
  use mod_memory, only : memory_deallo
  implicit none 
  integer(8), save :: memso(2) = 0_8

contains

  subroutine Symbolical_CSR_LU(nbrows,iA,jA,iL,jL,iU,jU,memor_opt)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    Symbolical_CSR_LU
    ! DESCRIPTION
    !    This routine performs Symbolical CSR LU factorization of a matrix A stored in the CSR format.
    ! INPUT ARGUMENTS
    !    NBROWS .... Number of rows of matrix A (dimension of matrix A)
    !    IA ........ CSR format: index vector for beginning of a row block for matrix A
    !    JA ........ CSR format: index vector for column numbers for matrix A
    ! OUTPUT ARGUMENTS 
    !    IL ........ CSR format: index vector for beginning of a row block for lower triangular matrix L
    !    JL ........ CSR format: index vector for column numbers for lower triangular matrix L
    !    IU ........ CSR format: index vector for beginning of a row block for upper triangular matrix U
    !    JU ........ CSR format: index vector for column numbers for upper triangular matrix U
    ! USES
    ! USED BY 
    !***
    !---------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)              :: nbrows
    integer(ip), intent(in)              :: iA(*),jA(*)
    integer(ip), pointer                 :: iL(:),jL(:)
    integer(ip), pointer                 :: iU(:),jU(:)

    integer(8),  intent(inout), optional :: memor_opt(2)
    integer(8)                           :: memor(2)    
    integer(ip)                          :: i,j,k,col,nxt,totalL,totalU
    integer(ip)                          :: rowIA_ini,rowIA_diag,rowIA_fin  
    integer(ip)                          :: minColL,maxColL,minColU,maxColU
    integer(ip), pointer                 :: firstL(:),firstU(:),seen(:),nextU(:),iwa(:)   
    integer(ip), pointer                 :: nzL(:),nzU(:)       !Local arrays for L/U matrices 
    type(i1p),   pointer                 :: ptL(:),ptU(:)

    if( present(memor_opt) ) then
       memor = memor_opt
    else
       memor = memso
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
          !allocate( ptL(i) % l(nxt) )
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
          !allocate( ptU(i) % l(nxt) )
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

    end if
    
    if( present(memor_opt) ) then
       memor_opt = memor 
    else
       memso = memor
    end if
    
  end subroutine Symbolical_CSR_LU

  subroutine Numerical_CSR_LU1(nbrows,iA,jA,An,iL,jL,Ln,iU,jU,Un,sing)

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
    integer(ip), intent(in)   :: nbrows
    integer(ip), intent(in)   :: iA(*), jA(*)
    real(rp),    intent(in)   :: An(*)
    integer(ip), intent(in)   :: iL(*), jL(*)
    real(rp),    intent(out)  :: Ln(*)
    integer(ip), intent(in)   :: iU(*), jU(*)
    real(rp),    intent(out)  :: Un(*)
    integer(ip), intent(out)  :: sing
    !
    ! Local variables
    !
    integer(ip)               :: i,j,k,col,coli,dia  
    real(rp)                  :: pivot
    real(rp), pointer         :: wa(:)      
    !
    ! Allocate working array, wa
    !
    nullify(wa)
    call memory_alloca(memso,'WA','Numerical_CSR_LU1',wa,nbrows,'DO_NOT_INITIALIZE')

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
       do j = iA(i), iA(i+1)-1
          col     = jA(j)
          wa(col) = An(j)
       end do
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
    call memory_deallo(memso,'WA','Numerical_CSR_LU1',wa)

  end subroutine Numerical_CSR_LU1

  subroutine Numerical_CSR_LUM(nbrows,ndof,iA,jA,An,iL,jL,Ln,iU,jU,Un,sing)

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
    integer(ip), intent(in)   :: nbrows, ndof
    integer(ip), intent(in)   :: iA(*)
    integer(ip), intent(in)   :: jA(*)
    real(rp),    intent(in)   :: An(ndof,ndof,*)
    integer(ip), intent(in)   :: iL(*)
    integer(ip), intent(in)   :: jL(*)
    real(rp),    intent(out)  :: Ln(ndof,ndof,*)
    integer(ip), intent(in)   :: iU(*)
    integer(ip), intent(in)   :: jU(*)
    real(rp),    intent(out)  :: Un(ndof,ndof,*)
    integer(ip), intent(out)  :: sing
    !
    ! Local variables
    !
    integer(ip)               :: i,j,k,l,s,t,col,dia,coli
    real(rp)                  :: pivot
    real(rp), pointer         :: wa(:,:,:)      
    !
    ! Allocate working array, wa
    !
    nullify(wa)
    call memory_alloca(memso,'WA','Numerical_CSR_LUM',wa,ndof,ndof,nbrows,'DO_NOT_INITIALIZE')

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
       do j= iA(i), iA(i+1)-1
          col = jA(j)
          do l= 1, ndof
             do k= 1, ndof
                wa(k,l,col) = An(l,k,j)
             end do
          end do
       end do
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
    call memory_deallo(memso,'WA','Numerical_CSR_LUM',wa)

  end subroutine Numerical_CSR_LUM

  subroutine CSR_LU_Factorization(nbrows,ndof,iA,jA,An,iL,jL,Ln,iU,jU,Un,sing)

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

    implicit none
    !
    ! Dummy arguments
    !
    integer(ip), intent(in)   :: nbrows !< Number of rows of matrix A (dimension of matrix A)
    integer(ip), intent(in)   :: ndof   !< Number of degrees of freedom in each node
    integer(ip), intent(in)   :: iA(*)  !< CSR format: index vector for beginning of a row block for matrix A
    integer(ip), intent(in)   :: jA(*)  !< CSR format: index vector for column numbers for matrix A
    real(rp),    intent(in)   :: An(*)  !< Sparse complex matrix A in BCSR (Blocked CSR Row) format
    integer(ip), pointer      :: iL(:)  !< CSR format: index vector for beginning of a row block for lower triangular matrix L
    integer(ip), pointer      :: jL(:)  !< CSR format: index vector for column numbers for lower triangular matrix L
    real(rp),    pointer      :: Ln(:)  !< Sparse complex lower triangular matrix L in BCSR (Blocked CSR Row) format
    integer(ip), pointer      :: iU(:)  !< CSR format: index vector for beginning of a row block for upper triangular matrix U
    integer(ip), pointer      :: jU(:)  !< CSR format: index vector for column numbers for upper triangular matrix U
    real(rp),    pointer      :: Un(:)  !< Sparse complex upper triangular matrix U in BCSR (Blocked CSR Row) format
    integer(ip), intent(out)  :: sing   !< Singularity marker
    !
    ! Local variables
    !
    integer(ip)               :: ndof2,nnzL,nnzU

    if( associated(Ln) .or. associated(Un) ) then
       write(*,*) 'Error!! Pointers already associated.'
    else
       !
       ! Perform Symbolical factorization
       !
       call Symbolical_CSR_LU(nbrows,iA,jA,iL,jL,iU,jU)
       !
       ! Number of non-zero elements in L/U
       !
       nnzL = iL(nbrows+1)-1
       nnzU = iU(nbrows+1)-1
       !
       ! Allocate memory for L/U
       !
       ndof2 = ndof*ndof
       call memory_alloca(memso,'LN','CSR_LU_Factorization',Ln,nnzL*ndof2)
       call memory_alloca(memso,'UN','CSR_LU_Factorization',Un,nnzU*ndof2) 

       !do i= 1, nnzL*ndof2
       !   Ln(i) = 0.0_rp
       !end do
       !do i= 1, nnzU*ndof2
       !   Un(i) = 0.0_rp
       !end do
       !
       ! Perform Numerical factorization dependig on number of degrees of freedom
       !
       if( ndof == 1 ) then
          call Numerical_CSR_LU1(nbrows,iA,jA,An,iL,jL,Ln,iU,jU,Un,sing)
       else
          call Numerical_CSR_LUM(nbrows,ndof,iA,jA,An,iL,jL,Ln,iU,jU,Un,sing)  
       end if

    end if

  end subroutine CSR_LU_Factorization

  subroutine CSR_Lsol(nbrows,ndof,invpR,iL,jL,Ln,iU,jU,Un,b,x)

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

    implicit none
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
       do j = 1,nbrows
          do i = 1,ndof
             x(i,j) = b(i,j)
          end do
       end do
    else 
       do j = 1,nbrows
          k = invpR(j)
          do i = 1,ndof
             x(i,j) = b(i,k)
          end do
       end do
    end if

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

  end subroutine CSR_Lsol

  subroutine CSR_Usol(nbrows,ndof,invpC,iU,jU,Un,b,x)

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

    implicit none
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

    if( .not. associated(invpC) ) then
       do j= 1, nbrows
          do i= 1, ndof
             x(i,j) = b(i,j)
          end do
       end do
    else
       do j= 1, nbrows
          col = invpC(j)
          do i= 1, ndof
             x(i,col) = b(i,j)
          end do
       end do
    end if

  end subroutine CSR_Usol

  subroutine CSR_LUsol(nbrows,ndof,invpR,invpC,iL,jL,Ln,iU,jU,Un,b,x)

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

    implicit none
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
    integer(ip)               :: i,j

    nullify(b_tmp)
    call memory_alloca(memso,'B_TMP','CSR_LUsol',b_tmp,ndof,nbrows,'DO_NOT_INITIALIZE')

    do j = 1,nbrows
       do i = 1,ndof
          b_tmp(i,j) = b(i,j)
       end do
    end do

    call CSR_Lsol(nbrows,ndof,invpR,iL,jL,Ln,iU,jU,Un,b_tmp,x)
    call CSR_Usol(nbrows,ndof,invpC,iU,jU,Un,x,b_tmp)

    do j = 1,nbrows
       do i = 1,ndof
          x(i,j) = b_tmp(i,j)
       end do
    end do

    call memory_deallo(memso,'B_TMP','CSR_LUsol',b_tmp)

  end subroutine CSR_LUsol

  subroutine CSR_LUfin(iL,jL,Ln,iU,jU,Un,invpR,invpC)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    CSR_LUfin
    ! DESCRIPTION
    !    This routine solves LU triangular linear systems of equations.eallocates memoryd
    ! INPUT ARGUMENTS
    ! OUTPUT ARGUMENTS 
    !     X ........ Solution vector of the system
    !---------------------------------------------------------------------------------------------------------------------------

    implicit none

    !Dummy arguments
    integer(ip), pointer           :: iL(:),jL(:)
    real(rp),    pointer           :: Ln(:)
    integer(ip), pointer           :: iU(:),jU(:)
    real(rp),    pointer           :: Un(:)
    integer(ip), pointer, optional :: invpR(:)
    integer(ip), pointer, optional :: invpC(:)

    call memory_deallo(memso,'iL','CSR_LUfin',iL)
    call memory_deallo(memso,'jL','CSR_LUfin',jL)
    call memory_deallo(memso,'Ln','CSR_LUfin',Ln)
    call memory_deallo(memso,'iU','CSR_LUfin',iU)
    call memory_deallo(memso,'jU','CSR_LUfin',jU)
    call memory_deallo(memso,'Un','CSR_LUfin',Un)
    if( present(invpR) ) call memory_deallo(memso,'invpR','CSR_LUfin',invpR)
    if( present(invpC) ) call memory_deallo(memso,'invpC','CSR_LUfin',invpC)

  end subroutine CSR_LUfin

  subroutine Symbolical_CSR_LU_Deallocate(iL,jL,iU,jU,memor_opt)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    CSR_LUfin
    ! DESCRIPTION
    !    This routine solves LU triangular linear systems of equations.eallocates memoryd
    ! INPUT ARGUMENTS
    ! OUTPUT ARGUMENTS 
    !     X ........ Solution vector of the system
    !---------------------------------------------------------------------------------------------------------------------------

    implicit none

    !Dummy arguments
    integer(ip),                pointer  :: iL(:),jL(:)
    integer(ip),                pointer  :: iU(:),jU(:)
    integer(8),  intent(inout), optional :: memor_opt(2)
    integer(8)                           :: memor(2)
    
    if( present(memor_opt) ) then
       memor = memor_opt
    else
       memor = memso
    end if
    
    call memory_deallo(memor,'iL','CSR_LUfin',iL)
    call memory_deallo(memor,'jL','CSR_LUfin',jL)
    call memory_deallo(memor,'iU','CSR_LUfin',iU)
    call memory_deallo(memor,'jU','CSR_LUfin',jU)

    if( present(memor_opt) ) then
       memor_opt = memor 
    else
       memso = memor
    end if
    
  end subroutine Symbolical_CSR_LU_Deallocate

  subroutine CSR_SpMV(nbrows,ndof,iA,jA,invpR,An,x,y)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    CSR_SpMV
    ! DESCRIPTION
    !    This routine solves a triangular linear system of equations where the system matrix is L - lower triangular matrix
    ! INPUT ARGUMENTS
    !    NBROWS .... Dimension of the sysytem
    !    NDOF ...... Number of degrees of freedom in each node
    !    INVPR ..... Permutation array: NEW = INVPR(OLD)
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

    implicit none
    !
    ! Dummy arguments
    !
    integer(ip), intent(in)   :: nbrows,ndof
    integer(ip), intent(in)   :: iA(*),jA(*)
    integer(ip), pointer      :: invpR(:)
    real(rp),    intent(in)   :: An(ndof,ndof,*)
    real(rp),    intent(in)   :: x(ndof,*)
    real(rp),    intent(out)  :: y(ndof,*)
    integer(ip)               :: ii,jj,iz,ki,kj,jjold,iiold
    real(rp),    allocatable  :: y_tmp(:,:)

    if( .not. associated(invpR) ) then
       do ii = 1,nbrows
          y(1:ndof,ii) = 0.0_rp
          do iz = iA(ii),iA(ii+1)-1
             jj = jA(iz)
             do ki = 1,ndof
                do kj = 1,ndof
                   y(ki,ii) = y(ki,ii) + An(kj,ki,iz) * x(kj,jj)
                end do
             end do
          end do
       end do
    else
       !
       ! x is in old numbering old = invpr(new)
       ! 
       allocate(y_tmp(ndof,nbrows))
       do ii = 1,nbrows
          y_tmp(1:ndof,ii) = 0.0_rp
          do iz = iA(ii),iA(ii+1)-1
             jj = jA(iz)
             jjold = invpR(jj)
             do ki = 1,ndof
                do kj = 1,ndof
                   y_tmp(ki,ii) = y_tmp(ki,ii) + An(kj,ki,iz) * x(kj,jjold)
                end do
             end do
          end do
       end do
       do ii = 1,nbrows
          iiold = invpR(ii)
          y(1:ndof,iiold) = y_tmp(1:ndof,ii)
       end do
       deallocate(y_tmp)
    end if

  end subroutine CSR_SpMV

  subroutine CSR_Permute_and_Copy_matrix(nbrows,ndof,iAin,jAin,iAout,jAout,invpR,Ain,Aout)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    CSR_SpMV
    ! DESCRIPTION
    !    This routine copy a matrix to another
    ! INPUT ARGUMENTS
    !    NBROWS .... Dimension of the sysytem
    !    NDOF ...... Number of degrees of freedom in each node
    !    INVPR ..... Permutation array: NEW = INVPR(OLD)
    !     B ........ RHS of the system
    ! OUTPUT ARGUMENTS 
    !     X ........ Solution vector of the system
    !---------------------------------------------------------------------------------------------------------------------------

    implicit none
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

  end subroutine CSR_Permute_and_Copy_matrix

end module mod_csrdir
