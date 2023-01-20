!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Output
!> @{
!> @file    mod_matrix_market.f90
!> @author  houzeaux
!> @date    2018-09-22
!> @brief   Output
!> @details Output matrix market format
!>
!-----------------------------------------------------------------------

module mod_matrix_market

  use def_kintyp_basic,      only : ip,rp,lg
  use mod_optional_argument, only : optional_argument
  use mod_iofile_basic,      only : iofile_available_unit
  use mod_memory,            only : memory_size
  use def_mat_fmt,           only : mat
  use def_mat_csr,           only : mat_csr
  use def_mat_coo,           only : mat_coo
  use mod_strings,           only : add_extension
  implicit none
  private

  interface matrix_market_matrix
     module procedure matrix_market_matrix_mat,&
          &           matrix_market_matrix_csr_rp,&
          &           matrix_market_matrix_csr_ip
  end interface matrix_market_matrix

  interface matrix_market_vector
     module procedure matrix_market_vector_rp_0,&
          &           matrix_market_vector_rp_1,&
          &           matrix_market_vector_ip_1,&
          &           matrix_market_vector_ip_2
  end interface matrix_market_vector
  
  public :: matrix_market_matrix
  public :: matrix_market_vector

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-05-06
  !> @brief   Output in matrix market format
  !> @details Output a graph and optionally, matrix and arrays
  !> 
  !-----------------------------------------------------------------------

  subroutine matrix_market_matrix_mat(a,FILENAME,PERM)

    class(mat),                       intent(in) :: a
    character(*),  optional,          intent(in) :: FILENAME
    integer(ip),   optional, pointer, intent(in) :: PERM(:)
    
    select type ( a )

    class is ( mat_csr ) ; call a % output('MATRIX MARKET',FILENAME)
    class is ( mat_coo ) ; call a % output('MATRIX MARKET',FILENAME)
       
    end select


  end subroutine matrix_market_matrix_mat

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-05-06
  !> @brief   Output in matrix market format
  !> @details Output a matrix in matrix market format
  !> 
  !-----------------------------------------------------------------------

  subroutine matrix_market_matrix_csr_rp(ndofr,ndofc,n,ia,ja,a,FILENAME,PERM)

    integer(ip),                      intent(in) :: ndofr
    integer(ip),                      intent(in) :: ndofc
    integer(ip),                      intent(in) :: n
    integer(ip),             pointer, intent(in) :: ia(:)
    integer(ip),             pointer, intent(in) :: ja(:)
    real(rp),                pointer, intent(in) :: a(:)
    character(*),  optional,          intent(in) :: FILENAME
    integer(ip),   optional, pointer, intent(in) :: PERM(:)
    integer(ip)                                  :: iz,ii,jj,nz
    integer(ip)                                  :: idof1,idof2,iztot
    integer(ip)                                  :: nrows,ncols
    integer(4)                                   :: unit4
    character(150)                               :: filename_loc

    if( associated(ia) ) then
       !
       ! Dimensions
       !
       nrows = ndofr*n
       ncols = maxval(ja) * ndofc
       nz    = (ia(n+1)-1)*ndofr*ndofc
       !
       ! Unit and files
       !
       filename_loc = optional_argument('MATRIXMARKET_MATRIX.mtx',FILENAME)
       unit4        = iofile_available_unit(90_ip)
       open(unit=unit4,file=trim(filename_loc),status='unknown')
       !
       ! Matrix 
       !
       write(unit4,'(a)')         '%%MatrixMarket matrix coordinate real general'
       write(unit4,'(3(1x,i12))') nrows,ncols,nz

       do ii = 1,n
          do idof2 = 1,ndofr
             do iz = ia(ii),ia(ii+1)-1
                jj = ja(iz)
                do idof1 = 1,ndofc
                   iztot = (iz-1)*ndofc*ndofr+(idof2-1)*ndofr+idof1
                   if( present(PERM) ) then
                      write(unit4,'(2(1x,i9),1x,e13.6)') PERM((ii-1)*ndofr+idof2),PERM((jj-1)*ndofc+idof1),a(iztot)
                   else
                      write(unit4,'(2(1x,i9),1x,e13.6)') (ii-1)*ndofr+idof2,(jj-1)*ndofc+idof1,a(iztot)
                   end if
                end do
             end do
          end do
       end do

       close(unit=unit4)

    end if

  end subroutine matrix_market_matrix_csr_rp

  subroutine matrix_market_matrix_csr_ip(ndofr,ndofc,n,ia,ja,a,FILENAME,PERM,PERMA)

    integer(ip),                      intent(in) :: ndofr
    integer(ip),                      intent(in) :: ndofc
    integer(ip),                      intent(in) :: n
    integer(ip),             pointer, intent(in) :: ia(:)
    integer(ip),             pointer, intent(in) :: ja(:)
    integer(ip),   optional, pointer, intent(in) :: a(:)
    character(*),  optional,          intent(in) :: FILENAME
    integer(ip),   optional, pointer, intent(in) :: PERM(:)
    integer(ip),   optional, pointer, intent(in) :: PERMA(:)
    integer(ip)                                  :: iz,ii,jj,nz
    integer(ip)                                  :: idof1,idof2,iztot
    integer(ip)                                  :: nrows,ncols
    integer(4)                                   :: unit4
    character(150)                               :: filename_loc

    if( associated(ia) ) then
       !
       ! Dimensions
       !
       nrows = n * ndofr
       ncols = maxval(ja) * ndofc
       nz    = (ia(n+1)-1)*ndofr*ndofc
       !
       ! Unit and files
       !
       filename_loc = optional_argument('MATRIXMARKET_MATRIX.mtx',FILENAME)
       unit4        = iofile_available_unit(90_ip)
       open(unit=unit4,file=trim(filename_loc),status='unknown')
       !
       ! Matrix 
       !
       write(unit4,'(a)')         '%%MatrixMarket matrix coordinate real general'
       write(unit4,'(3(1x,i12))') nrows,ncols,nz

       do ii = 1,n
          do idof2 = 1,ndofr
             do iz = ia(ii),ia(ii+1)-1
                jj = ja(iz)
                do idof1 = 1,ndofc
                   iztot = (iz-1)*ndofc*ndofr+(idof2-1)*ndofr+idof1
                   if( present(a) ) then
                      if( present(PERM) ) then
                         if( present(PERMA) ) then
                            write(unit4,'(2(1x,i9),1x,i9)') PERM((ii-1)*ndofr+idof2),PERM((jj-1)*ndofc+idof1),PERMA(a(iztot))
                         else
                            write(unit4,'(2(1x,i9),1x,i9)') PERM((ii-1)*ndofr+idof2),PERM((jj-1)*ndofc+idof1),a(iztot)
                         end if
                      else
                         if( present(PERMA) ) then
                            write(unit4,'(2(1x,i9),1x,i9)') (ii-1)*ndofr+idof2,(jj-1)*ndofc+idof1,PERMA(a(iztot))
                         else
                            write(unit4,'(2(1x,i9),1x,i9)') (ii-1)*ndofr+idof2,(jj-1)*ndofc+idof1,a(iztot)
                         end if
                      end if
                   else
                      if( present(PERM) ) then
                         write(unit4,'(2(1x,i9),1x,i9)') PERM((ii-1)*ndofr+idof2),PERM((jj-1)*ndofc+idof1)
                      else
                         write(unit4,'(2(1x,i9),1x,i9)') (ii-1)*ndofr+idof2,(jj-1)*ndofc+idof1
                      end if
                   end if
                end do
             end do
          end do
       end do

       close(unit=unit4)

    end if

  end subroutine matrix_market_matrix_csr_ip
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-05-06
  !> @brief   Output in matrix market format
  !> @details Output a vector in matrix market format
  !> 
  !-----------------------------------------------------------------------

  subroutine matrix_market_vector_ip_1(x,FILENAME,PERM)

    integer(ip),             pointer, intent(in) :: x(:)
    character(*),  optional,          intent(in) :: FILENAME
    integer(ip),   optional, pointer, intent(in) :: PERM(:)
    integer(ip)                                  :: n,ndofn,ii
    character(150)                               :: filename_loc
    integer(4)                                   :: unit4
    !
    ! Dimensions
    !
    n     = memory_size(x)
    ndofn = 1
    !
    ! Unit and files
    !
    filename_loc = optional_argument('MATRIXMARKET_ARRAY.mtx',FILENAME)
    unit4        = iofile_available_unit(90_ip)
    open(unit=unit4,file=trim(filename_loc),status='unknown')
    !
    ! Vector 
    !
    write(unit4,'(a)') '%%MatrixMarket matrix array real general'
    write(unit4,'(3(1x,i12))') n * ndofn,1_ip
    if( present(PERM) ) then
       do ii = 1,n * ndofn
          write(unit4,'(1x,i9,1x,i9)') PERM(ii),x(ii)
       end do
    else
       do ii = 1,n * ndofn
          write(unit4,'(1x,i9,1x,i9)') ii,x(ii)
       end do
    end if
    
    close(unit=unit4)

  end subroutine matrix_market_vector_ip_1
  
  subroutine matrix_market_vector_ip_2(x,FILENAME,PERM,PERMX)

    integer(ip),             pointer, intent(in) :: x(:,:)
    character(*),  optional,          intent(in) :: FILENAME
    integer(ip),   optional, pointer, intent(in) :: PERM(:)
    integer(ip),   optional, pointer, intent(in) :: PERMX(:)
    integer(ip)                                  :: n,ndofn,ii
    character(150)                               :: filename_loc
    integer(4)                                   :: unit4
    !
    ! Dimensions
    !
    ndofn = memory_size(x,1_ip)
    n     = memory_size(x,2_ip)
    !
    ! Unit and files
    !
    filename_loc = optional_argument('MATRIXMARKET_ARRAY.mtx',FILENAME)
    unit4        = iofile_available_unit(90_ip)
    open(unit=unit4,file=trim(filename_loc),status='unknown')
    !
    ! Vector 
    !
    write(unit4,'(a)') '%%MatrixMarket matrix array real general'
    write(unit4,'(3(1x,i12))') n,ndofn
    if(      present(PERM) .and. (.not. present(PERMX)) ) then
       do ii = 1,n 
          write(unit4,'(1x,i9,10(1x,i9))') PERM(ii),x(:,ii)
       end do
    else if( present(PERM) .and. present(PERMX) ) then
       do ii = 1,n 
          write(unit4,'(1x,i9,10(1x,i9))') PERM(ii),PERMX(x(:,ii))
       end do
    else if( (.not. present(PERM)) .and. present(PERMX) ) then
       do ii = 1,n 
          write(unit4,'(1x,i9,10(1x,i9))') ii,PERMX(x(:,ii))
       end do
    else
       do ii = 1,n 
          write(unit4,'(1x,i9,10(1x,i9))') ii,x(:,ii)
       end do
    end if
    
    close(unit=unit4)

  end subroutine matrix_market_vector_ip_2
  
  subroutine matrix_market_vector_rp_1(x,FILENAME,PERM)

    real(rp),                pointer, intent(in) :: x(:)
    character(*),  optional,          intent(in) :: FILENAME
    integer(ip),   optional, pointer, intent(in) :: PERM(:)
    integer(ip)                                  :: n,ndofn,ii
    character(150)                               :: filename_loc
    integer(4)                                   :: unit4
    !
    ! Dimensions
    !
    n     = memory_size(x)
    ndofn = 1
    !
    ! Unit and files
    !
    filename_loc = optional_argument('MATRIXMARKET_ARRAY.mtx',FILENAME)
    unit4        = iofile_available_unit(90_ip)
    open(unit=unit4,file=trim(filename_loc),status='unknown')
    !
    ! Vector 
    !
    write(unit4,'(a)') '%%MatrixMarket matrix array real general'
    write(unit4,'(3(1x,i12))') n * ndofn,1_ip
    if( present(PERM) ) then
       do ii = 1,n * ndofn
          write(unit4,'(1(1x,i9),1x,e13.6)') PERM(ii),x(ii)
       end do
    else
       do ii = 1,n * ndofn
          write(unit4,'(1(1x,i9),1x,e13.6)') ii,x(ii)
       end do
    end if
    
    close(unit=unit4)

  end subroutine matrix_market_vector_rp_1
  
  subroutine matrix_market_vector_rp_0(ndofn,n,x,FILENAME,PERM)

    integer(ip),                      intent(in) :: ndofn
    integer(ip),                      intent(in) :: n
    real(rp),                         intent(in) :: x(*)
    character(*),  optional,          intent(in) :: FILENAME
    integer(ip),   optional, pointer, intent(in) :: PERM(:)
    integer(ip)                                  :: ii
    character(150)                               :: filename_loc
    integer(4)                                   :: unit4
    !
    ! Unit and files
    !
    filename_loc = optional_argument('MATRIXMARKET_ARRAY.mtx',FILENAME)
    unit4        = iofile_available_unit(90_ip)
    open(unit=unit4,file=trim(filename_loc),status='unknown')
    !
    ! Vector
    !
    write(unit4,'(a)') '%%MatrixMarket matrix array real general'
    write(unit4,'(3(1x,i12))') n * ndofn,1_ip
    if( present(PERM) ) then
       do ii = 1,n * ndofn
          write(unit4,'(1(1x,i9),1x,e13.6)') PERM(ii),x(ii)
       end do
    else
       do ii = 1,n * ndofn
          write(unit4,'(1(1x,i9),1x,e13.6)') ii,x(ii)
       end do
    end if
    
    close(unit=unit4)

  end subroutine matrix_market_vector_rp_0
  
end module mod_matrix_market
!> @}

