!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_csrdirx

  !------------------------------------------------------------------------
  !****f* soldir/mod_csrdirx
  ! NAME 
  !    mod_csrdirx
  ! DESCRIPTION
  !    Complex sparse direct solver module
  ! USES
  ! USED BY 
  !***
  !------------------------------------------------------------------------

  use def_kintyp, only    : ip,rp
  use def_master, only    : kfl_paral 
  implicit none 

contains

  subroutine Symbolicalx_CSR_LU(nbrows,iA,jA,iL,jL,iU,jU)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    Symbolicalx_CSR_LU
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

    use def_kintyp,  only     :  ip,rp,i1p
    implicit none
    integer(ip), intent(in)   :: nbrows
    integer(ip), intent(in)   :: iA(*),jA(*)
    integer(ip), pointer      :: iL(:),jL(:)
    integer(ip), pointer      :: iU(:),jU(:)
    integer(4)                :: istat
    integer(ip)               :: i,j,k,col,nxt,totalL,totalU
    integer(ip)               :: rowIA_ini,rowIA_diag,rowIA_fin  
    integer(ip)               :: minColL,maxColL,minColU,maxColU
    integer(ip), pointer      :: firstL(:),firstU(:),seen(:),nextU(:),iwa(:)   
    integer(ip), pointer      :: nzL(:),nzU(:)       !Local arrays for L/U matrices 
    type(i1p),   allocatable  :: ptL(:),ptU(:)

    if( associated(iL) .or. associated(jL) .or. associated(iU) .or. associated(jU) ) then
       write(*,*) 'Error!! Pointers already associated.'
    else
       !
       ! Alloc local work arrays
       !
       allocate(firstL(nbrows),stat=istat)
       if( istat /= 0 ) call memerr(0_ip,'firstL','Symbolical_CSR_LU',0_ip)
       
       allocate(firstU(nbrows),stat=istat)
       if( istat /= 0 ) call memerr(0_ip,'firstU','Symbolical_CSR_LU',0_ip)

       allocate(seen(nbrows),  stat=istat)
       if( istat /= 0 ) call memerr(0_ip,'seen','Symbolical_CSR_LU',0_ip)

       allocate(nextU(nbrows), stat=istat)
       if( istat /= 0 ) call memerr(0_ip,'nextU','Symbolical_CSR_LU',0_ip)

       allocate(iwa(nbrows),   stat=istat)
       if( istat /= 0 ) call memerr(0_ip,'iwa','Symbolical_CSR_LU',0_ip)

       allocate(nzL(nbrows),   stat=istat)
       if( istat /= 0 ) call memerr(0_ip,'nzL','Symbolical_CSR_LU',0_ip)

       allocate(nzU(nbrows),   stat=istat)
       if( istat /= 0 ) call memerr(0_ip,'nzU','Symbolical_CSR_LU',0_ip)

       allocate( ptL(nbrows), stat=istat) 
       if( istat /= 0 ) call memerr(0_ip,'ptL','Symbolical_CSR_LU',0_ip)

       allocate( ptU(nbrows), stat=istat) 
       if( istat /= 0 ) call memerr(0_ip,'ptU','Symbolical_CSR_LU',0_ip)
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
          allocate( ptL(i) % l(nxt) , stat = istat )
          if( istat /= 0 ) call memerr(0_ip,'ptL(i)%l','Symbolical_CSR_LU',0_ip)
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
          allocate( ptU(i) % l(nxt) , stat = istat )
          if( istat /= 0 ) call memerr(0_ip,'ptU(i)%l','Symbolical_CSR_LU',0_ip)
          do k= 1, nxt
             ptU(i) % l(k) = iwa(k)
          end do
       end do
       !
       ! Put L and U in CSR format 
       !
       allocate(iL(nbrows+1),stat=istat)
       if( istat /= 0 ) call memerr(0_ip,'iL','Symbolical_CSR_LU',0_ip)

       allocate(iU(nbrows+1),stat=istat)
       if( istat /= 0 ) call memerr(0_ip,'iU','Symbolical_CSR_LU',0_ip)

       allocate(jL(totalL),stat=istat)
       if( istat /= 0 ) call memerr(0_ip,'jL','Symbolical_CSR_LU',0_ip)

       allocate(jU(totalU),stat=istat)
       if( istat /= 0 ) call memerr(0_ip,'jU','Symbolical_CSR_LU',0_ip)
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
       do i= 1, nbrows
          deallocate(ptL(i) % l) 
          if( istat /= 0 ) call memerr(2_ip,'ptL(i)%l','Symbolical_CSR_LU',0_ip)
          deallocate(ptU(i) % l)
          if( istat /= 0 ) call memerr(2_ip,'ptU(i)%l','Symbolical_CSR_LU',0_ip)
       end do
       deallocate(ptU,stat=istat)
       if( istat /= 0 ) call memerr(2_ip,'ptU','Symbolical_CSR_LU',0_ip)

       deallocate(ptL,stat=istat)
       if( istat /= 0 ) call memerr(2_ip,'ptL','Symbolical_CSR_LU',0_ip)

       deallocate(nzU,stat=istat)
       if( istat /= 0 ) call memerr(2_ip,'nzU','Symbolical_CSR_LU',0_ip)

       deallocate(nzL,stat=istat)
       if( istat /= 0 ) call memerr(2_ip,'nzL','Symbolical_CSR_LU',0_ip)

       deallocate(iwa,stat=istat)
       if( istat /= 0 ) call memerr(2_ip,'iwa','Symbolical_CSR_LU',0_ip)

       deallocate(nextU,stat=istat)
       if( istat /= 0 ) call memerr(2_ip,'nextU','Symbolical_CSR_LU',0_ip)

       deallocate(seen,stat=istat)
       if( istat /= 0 ) call memerr(2_ip,'seen','Symbolical_CSR_LU',0_ip)

       deallocate(firstU,stat=istat)
       if( istat /= 0 ) call memerr(2_ip,'firstU','Symbolical_CSR_LU',0_ip)

       deallocate(firstL,stat=istat)
       if( istat /= 0 ) call memerr(2_ip,'firstL','Symbolical_CSR_LU',0_ip)

    end if

  end subroutine Symbolicalx_CSR_LU

  subroutine Numericalx_CSR_LU1(nbrows,iA,jA,An,iL,jL,Ln,iU,jU,Un,sing)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    Numericalx_CSR_LU1
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

    !Declaration statements
    use def_kintyp,  only   : ip,rp
    implicit none

    !Dummy arguments
    integer(ip), intent(in)   :: nbrows
    integer(ip), intent(in)   :: iA(*), jA(*)
    complex(rp), intent(in)   :: An(*)
    integer(ip), intent(in)   :: iL(*), jL(*)
    complex(rp), intent(out)  :: Ln(*)
    integer(ip), intent(in)   :: iU(*), jU(*)
    complex(rp), intent(out)  :: Un(*)
    integer(ip), intent(out)  :: sing

    !Local variables
    integer(4)                :: istat
    integer(ip)               :: i,j,k,col,coli,dia  
    complex(rp)               :: pivot
    complex(rp), pointer      :: wa(:)      

    !Allocate working array, wa
    allocate(wa(nbrows), stat=istat)
    if( istat == 0 ) then
       sing = 0

       !Main loop in rows
       DO i= 1, nbrows
          !Initialize wa with zeros and values of the matrix A 	
          do j= iL(i), iL(i+1)-1
             col = jL(j)
             wa(col) = (0.0_rp,0.0_rp)
          end do
          do j= iU(i), iU(i+1)-1
             col = jU(j)
             wa(col) = (0.0_rp,0.0_rp)
          end do
          do j= iA(i), iA(i+1)-1
             col = jA(j)
             wa(col) = An(j)
          end do

          !Factorize the i-th row of the matrix A
          !For all the elements in the i-th row before diagonal
          do j= iL(i), iL(i+1)-1
             col = jL(j) 
             dia = iU(col)
             pivot = wa(col) * Un(dia)       !U diagonal is stored inverted 
             wa(col) = pivot
             !For all the elements in the col-th row after diagonal
             do k= iU(col)+1, iU(col+1)-1
                coli = jU(k)
                wa(coli) = wa(coli) - pivot * Un(k)
             end do
          end do

          !Store factorized i-th row in L/U 
          pivot = wa(i)
          if( pivot /= (0.0_rp,0.0_rp) ) then
             !Store the elements before diagonal in L
             do j= iL(i), iL(i+1)-1
                col = jL(j)
                Ln(j) = wa(col)
             end do
             !Store the diagonal element in U, at the first place in the i-th row of the matrix U
             dia = iU(i) 
             Un(dia) = (1.0_rp,0.0_rp) / pivot       !Diagonal element is stored inverted
             !Store the elements after diagonal in U
             do j= iU(i)+1, iU(i+1)-1
                col = jU(j)
                Un(j) = wa(col)
             end do
          else
             write(*,*) 'Singular pivot in row ', i
             sing = -1
          end if
       END DO

       !Deallocate working array
       deallocate(wa,stat=istat) 
    else
       write(*,*) 'Error!!! Out of memory!'
    end if

  end subroutine Numericalx_CSR_LU1

  subroutine Numericalx_CSR_LUM(nbrows,ndof,iA,jA,An,iL,jL,Ln,iU,jU,Un,sing)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    Numericalx_CSR_LUM
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

    !Declaration statements
    use def_kintyp,  only   : ip,rp
    implicit none

    !Dummy arguments
    integer(ip), intent(in)   :: nbrows, ndof
    integer(ip), intent(in)   :: iA(*), jA(*)
    complex(rp), intent(in)   :: An(ndof,ndof,*)
    integer(ip), intent(in)   :: iL(*), jL(*)
    complex(rp), intent(out)  :: Ln(ndof,ndof,*)
    integer(ip), intent(in)   :: iU(*), jU(*)
    complex(rp), intent(out)  :: Un(ndof,ndof,*)
    integer(ip), intent(out)  :: sing

    !Local variables
    integer(4)                :: istat
    integer(ip)               :: i,j,k,l,s,t,col,dia,coli
    complex(rp)               :: pivot
    complex(rp), pointer      :: wa(:,:,:)      

    !Allocate working array, wa
    allocate(wa(ndof,ndof,nbrows), stat=istat)

    if( istat == 0 ) then
       sing = 0

       !Main loop in rows
       DO i= 1, nbrows
          !Initialize wa with zeros and values of the matrix A 	
          do j= iL(i), iL(i+1)-1
             col = jL(j)
             do l= 1, ndof
                do k= 1, ndof
                   wa(k,l,col) = (0.0_rp,0.0_rp)
                end do
             end do
          end do
          do j= iU(i), iU(i+1)-1
             col = jU(j)
             do l= 1, ndof
                do k= 1, ndof
                   wa(k,l,col) = (0.0_rp,0.0_rp)
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

          !Factorize all the blocks in the i-th row of the matrix A
          !For all the elements (blocks) in the i-th row before diagonal
          DO j= iL(i), iL(i+1)-1
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
          END DO
          !Factorize the lower triangle of the diagonal block
          DO k= 1, ndof
             do l= 1, k-1
                pivot = wa(k,l,i) * wa(l,l,i)
                wa(k,l,i) = pivot
                !For all the elements in the k-th row inside the diagonal block
                do s= l+1, ndof
                   wa(k,s,i) = wa(k,s,i) - pivot * wa(l,s,i)
                end do
             end do
             !Check if some of the diagonal elements of the diagonal block is zero
             if( wa(k,k,i) /= (0.0_rp,0.0_rp) ) then
                wa(k,k,i) = (1.0_rp,0.0_rp) / wa(k,k,i)
             else
                sing = -1                
                exit
             end if
          END DO
          !For all the blocks after diagonal
          DO s= iU(i)+1, iU(i+1)-1
             coli = jU(s)
             do k= 1, ndof
                do l= 1, k-1
                   do t= 1, ndof 
                      wa(k,t,coli) = wa(k,t,coli) - wa(l,t,coli) * wa(k,l,i)
                   end do
                end do
             end do
          END DO

          !Store factorized i-th row of blocks in L/U matrices
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
       END DO

       !Deallocate working array   
       deallocate(wa,stat=istat) 
    else
       write(*,*) 'Error Out of memory!'
    end if

  end subroutine Numericalx_CSR_LUM


  subroutine CSR_LU_Factorizationx(nbrows,ndof,iA,jA,An,iL,jL,Ln,iU,jU,Un,sing)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    CSR_LU_Factorizationx
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

    !Declaration statements
    use def_kintyp,  only   : ip,rp
!    use def_master, only    : kfl_paral
    implicit none

    !Dummy arguments
    integer(ip), intent(in)   :: nbrows,ndof
    integer(ip), intent(in)   :: iA(*),jA(*)
    complex(rp), intent(in)   :: An(*)
    integer(ip), pointer      :: iL(:),jL(:)
    complex(rp), pointer      :: Ln(:)
    integer(ip), pointer      :: iU(:),jU(:)
    complex(rp), pointer      :: Un(:)
    integer(ip), intent(out)  :: sing

    !Local variables
    integer(4)                :: istat
    integer(ip)               :: i,ndof2,nnzL,nnzU!,ierror
!    real(rp)                  :: cpu_syma1,cpu_syma2,cpu_syma

    if( associated(Ln) .or. associated(Un) ) then
       write(*,*) 'Error!! Pointers already associated.'
    else
  !!call cputim(cpu_syma1)
       !Perform Symbolical factorization
       call Symbolicalx_CSR_LU(nbrows,iA,jA,iL,jL,iU,jU)
  !!call cputim(cpu_syma2)
  !!cpu_syma = cpu_syma2 - cpu_syma1
  !!write(*,*)'Time to symbolical:',cpu_syma, kfl_paral
       !Number of non-zero elements in L/U
       nnzL = iL(nbrows+1)-1
       nnzU = iU(nbrows+1)-1
       !Allocate memory for L/U
       ndof2 = ndof*ndof
       allocate(Ln(nnzL*ndof2),stat=istat) 
       allocate(Un(nnzU*ndof2),stat=istat) 
       do i= 1, nnzL*ndof2
          Ln(i) = (0.0_rp,0.0_rp)
       end do
       do i= 1, nnzU*ndof2
          Un(i) = (0.0_rp,0.0_rp)
       end do
  !!call cputim(cpu_syma1)
       !Perform Numerical factorization dependig on number of degrees of freedom
       if( ndof==1 ) then
          call Numericalx_CSR_LU1(nbrows,iA,jA,An,iL,jL,Ln,iU,jU,Un,sing)
       else
          call Numericalx_CSR_LUM(nbrows,ndof,iA,jA,An,iL,jL,Ln,iU,jU,Un,sing)  
       end if
  !!call cputim(cpu_syma2)
  !!cpu_syma = cpu_syma2 - cpu_syma1
  !!write(*,*)'Time to numerical:',cpu_syma, kfl_paral
       !Print L/U
       !call print_matrix(nbrows,ndof,nnzL,iL,jL,Ln,'L')
       !call print_matrix(nbrows,ndof,nnzU,iU,jU,Un,'U')
       !open(unit=2,status='replace',file='outL.out',iostat=ierror)
       !if( ierror == 0 ) then
       !   write(2,*) nbrows, ndof, nnzL 
       !   write(2,*) (iL(i), i = 1,nbrows+1)
       !   write(2,*) (jL(i), i = 1, nnzL)
       !   write(2,*) (Ln(i), i = 1, nnzL*ndof2)
       !end if
       !open(unit=3,status='replace',file='outU.out',iostat=ierror)
       !if( ierror == 0 ) then
       !   write(3,*) nbrows, ndof, nnzU 
       !   write(3,*) (iU(i), i = 1,nbrows+1)
       !   write(3,*) (jU(i), i = 1, nnzU)
       !   write(3,*) (Un(i), i = 1, nnzU*ndof2)
       !end if
    end if

  end subroutine CSR_LU_Factorizationx

  subroutine CSR_Lsolx(nbrows,ndof,invpR,iL,jL,Ln,iU,jU,Un,b,x)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    CSR_Lsolx
    ! DESCRIPTION
    !    This routine solves a triangular linear system of equations where the system matrix is L - lower triangular matrix
    ! INPUT ARGUMENTS
    !    NBROWS .... Dimension of the sysytem
    !    NDOF ...... Number of degrees of freedom in each node
    !    INVPR ..... 
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

    !Declaration statements
    use def_kintyp,  only   : ip,rp
    implicit none

    !Dummy arguments
    integer(ip), intent(in)   :: nbrows,ndof
    integer(ip), pointer      :: invpR(:)
    integer(ip), intent(in)   :: iL(*),jL(*)
    complex(rp), intent(in)   :: Ln(ndof,ndof,*)
    integer(ip), intent(in)   :: iU(*),jU(*)
    complex(rp), intent(in)   :: Un(ndof,ndof,*)
    complex(rp), intent(in)   :: b(ndof,*)
    complex(rp), intent(out)  :: x(ndof,*)

    !Local variables
    integer(ip)               :: i,j,k,l,col

    if( .not. associated(invpR) ) then
       do j= 1, nbrows
          do i= 1, ndof
             x(i,j) = b(i,j)
          end do
       end do
    else 
       do j= 1, nbrows
          k = invpR(j)
          do i= 1, ndof
             x(i,j) = b(i,k)
          end do
       end do
    end if

    DO i= 1, nbrows
       do j= iL(i), iL(i+1)-1
          col = jL(j)
          do l= 1, ndof  
             do k= 1, ndof
                x(k,i) = x(k,i) - Ln(k,l,j) * x(l,col)
             end do
          end do
       end do
       !Update using the lower triangle of the diagonal block stored in the U matrix
       j = iU(i)
       do l= 1, ndof
          do k= l+1, ndof
             x(k,i) = x(k,i) - Un(k,l,j) * x(l,i)
          end do
       end do
    END DO

  end subroutine CSR_Lsolx

  subroutine CSR_Usolx(nbrows,ndof,invpC,iU,jU,Un,b,x)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    CSR_Usolx
    ! DESCRIPTION
    !    This routine solves a triangular linear system of equations where the system matrix is U - upper triangular matrix
    ! INPUT ARGUMENTS
    !    NBROWS .... Dimension of the sysytem
    !    NDOF ...... Number of degrees of freedom in each node
    !    INVPC ..... 
    !    IU ........ CSR format: index vector for beginning of a row block for upper triangular matrix U
    !    JU ........ CSR format: index vector for column numbers for upper triangular matrix U
    !    UN ........ Sparse complex upper triangular matrix U in BCSR (Blocked CSR Row) format
    !     B ........ RHS of the system
    ! OUTPUT ARGUMENTS 
    !     X ........ Solution vector of the system
    !---------------------------------------------------------------------------------------------------------------------------

    !Declaration statements
    use def_kintyp,  only   : ip,rp
    implicit none

    !Dummy arguments
    integer(ip), intent(in)   :: nbrows, ndof
    integer(ip), pointer      :: invpC(:)
    integer(ip), intent(in)   :: iU(*), jU(*)
    complex(rp), intent(in)   :: Un(ndof,ndof,*)
    complex(rp), intent(inout):: b(ndof,*)
    complex(rp), intent(out)  :: x(ndof,*)

    !Local variables
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

  end subroutine CSR_Usolx

  subroutine CSR_LUsolx(nbrows,ndof,invpR,invpC,iL,jL,Ln,iU,jU,Un,b,x)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    CSR_LUsolx
    ! DESCRIPTION
    !    This routine solves LU triangular linear systems of equations.
    ! INPUT ARGUMENTS
    !    NBROWS .... Dimension of the sysytem
    !    NDOF ...... Number of degrees of freedom in each node
    !    INVPR ..... 
    !    INVPC ..... 
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

    !Declaration statements
    use def_kintyp,  only   : ip,rp
!    use def_master, only    : kfl_paral
    implicit none

    !Dummy arguments
    integer(ip), intent(in)  :: nbrows,ndof
    integer(ip), pointer     :: invpR(:),invpC(:)
    integer(ip), intent(in)  :: iL(*),jL(*)
    complex(rp), intent(in)  :: Ln(ndof,ndof,*)
    integer(ip), intent(in)  :: iU(*),jU(*)
    complex(rp), intent(in)  :: Un(ndof,ndof,*)
    complex(rp), intent(in)  :: b(ndof,*)
    complex(rp), intent(out) :: x(ndof,*)
    complex(rp), pointer     :: b_tmp(:,:)

    !Local variables
    integer(ip)               :: i,j
    integer(4)                :: istat
!    real(rp)                  :: cpu_syma1,cpu_syma2,cpu_syma

    allocate( b_tmp(ndof,nbrows) , stat = istat )
    do j = 1,nbrows
       do i = 1,ndof
          b_tmp(i,j) = b(i,j)
       end do
    end do
  !!call cputim(cpu_syma1)
    call CSR_Lsolx(nbrows,ndof,invpR,iL,jL,Ln,iU,jU,Un,b_tmp,x)
  !!call cputim(cpu_syma2)
  !!cpu_syma = cpu_syma2 - cpu_syma1
  !!write(*,*)'Time to Lsolve:',cpu_syma, kfl_paral
  !!call cputim(cpu_syma1)
    call CSR_Usolx(nbrows,ndof,invpC,iU,jU,Un,x,b_tmp)
  !!call cputim(cpu_syma2)
  !!cpu_syma = cpu_syma2 - cpu_syma1
  !!write(*,*)'Time to Usolve:',cpu_syma, kfl_paral

    do j = 1,nbrows
       do i = 1,ndof
          x(i,j) = b_tmp(i,j)
       end do
    end do

    deallocate( b_tmp , stat = istat )

  end subroutine CSR_LUsolx

  subroutine CSR_LUfinx(iL,jL,Ln,iU,jU,Un)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    CSR_LUfinx
    ! DESCRIPTION
    !    This routine deallocates memory for LU solver.
    ! INPUT ARGUMENTS
    ! OUTPUT ARGUMENTS 
    !---------------------------------------------------------------------------------------------------------------------------

    !Declaration statements
    use def_kintyp,  only   : ip,rp
    implicit none

    !Dummy arguments
    integer(ip), pointer :: iL(:),jL(:)
    complex(rp), pointer :: Ln(:)
    integer(ip), pointer :: iU(:),jU(:)
    complex(rp), pointer :: Un(:)

    deallocate(iL)
    deallocate(jL)
    deallocate(Ln)
    deallocate(iU)
    deallocate(jU)
    deallocate(Un)

  end subroutine CSR_LUfinx

end module mod_csrdirx
