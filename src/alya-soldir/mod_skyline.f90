!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solver
!> @{
!> @file    mod_skyline.f90
!> @author  houzeaux
!> @date    2022-10-13
!> @brief   Module for skyline format based solvers
!> @details Subroutine to solve skyline factorization and solution
!-----------------------------------------------------------------------

module mod_skyline

  use def_kintyp, only :  ip,rp
  implicit none
  private

  public :: chofac
  public :: chosol
  public :: lufact
  public :: lusolv

contains

  subroutine chofac( n, nnz, nz, a, err )

    !  chofac computes the Cholesky factorization of a real symmetric
    !  positive definite matrix A stored in skyline format.

    integer(ip)          :: err, n, nnz
    integer(ip)          :: nz( n+1 )
    real(rp)             :: a( nnz )
    integer(ip)          :: i, j,  k
    integer(ip)          :: i0, j0
    integer(ip)          :: ipos,idif,irow,jcol
    real(rp)             :: temp, rdiag
    !
    ! test the input parameters.
    !
    err = 0
    if( n<0 ) then
       write(*,*)'n <0' 
       err = -1
    else if( nnz<0 ) then
       write(*,*)'nnz <0' 
       err = -2
    end if
    if( err.ne.0 ) then
       return
    end if
    !
    ! quick return if possible.
    !
    if( n==0.or.nnz==0 )  return
    !
    ! compute the cholesky factorization a = ll^t.
    !
    !
    ! initialize pointer
    !
    ipos = 2
    !
    ! compute first pivot
    !
    if(a(1)<0.0_rp)then
       write(*,*) 'first pivot negative'
       err=1
       return
    end if
    ! 
    a(1) = sqrt ( a(1))
    !    
    do  i = 2, n
       !
       ! set to 0 diagonal term
       !
       rdiag = 0.0_rp
       !
       ! column number of the first non zero in row i
       !
       jcol=i-(nz(i+1)-nz(i))+1
       !     
       ! compute elements jcol:(i-1) of line i.
       !     
       do  j = jcol , i-1
          !
          ! row number of first non zero in column j
          !
          irow=j-(nz(j+1)-nz(j))+1
          !
          ! check for non zero bounds to setup the pointers for the scalar product
          !
          if(jcol>irow)then
             idif = j-jcol
             i0 = nz(i)
             j0 = nz(j+1)-idif-1
          else
             idif = j-irow
             i0 = ipos-idif
             j0 = nz(j)
          end if
          !
          ! compute scalar product
          !
          temp = 0.0_rp
          do k= 1,idif
             temp = temp + a( i0 )*a( j0 )
             i0 = i0 + 1
             j0 = j0 + 1
          end do
          !
          a(ipos)=(a(ipos)-temp)/a(nz(j+1)-1)
          !
          ! accumulate in rdiag
          !
          rdiag=rdiag+a(ipos)*a(ipos)
          !
          ! move pointer
          !
          ipos=ipos+1
          !
       end do
       !
       ! compute the diagonal
       !
       a(ipos) = a(ipos)-rdiag
       if(a(ipos)<0.0_rp)then
          write(*,*)'pivot negative in line ',i
          err=i
          return
       end if

       a(ipos)=sqrt(a(ipos)) 
       ipos=ipos+1

    end do
    !
  end subroutine chofac

  subroutine chosol( n, nnz, nz, nrhs, a, b, ldb, err )
    !
    !  cho_solve solves a system of linear equations A*X = B with a symmetric
    !  positive definite matrix A stored in skyline format using the Cholesky
    !  factorization A = LL^T computed by chol_fact.
    !
    integer(ip)          :: err, ldb, n, nnz, nrhs
    integer(ip)          :: nz( n+1 )
    real(rp)             :: a( nnz ), b( ldb, * )
    integer(ip)          :: j, k
    integer(ip)          :: k0
    integer(ip)          :: jcnt, i
    real(rp)             :: temp
    !
    ! test the input parameters.
    !
    err = 0
    if( n<0 ) then
       write(*,*)'n <0' 
       err = -1
    else if( nnz<0 ) then
       write(*,*)'nnz <0' 
       err = -2
    else if( nrhs<0 ) then
       write(*,*)'nrhs <0' 
       err = -4
    else if( ldb<max( 1_ip, n ) ) then
       write(*,*)'ldb < max(1,n)' 
       err = -7
    end if
    if( err.ne.0 ) then
       return
    end if
    !
    ! quick return if possible.
    !
    if( n==0 .or. nnz==0 .or. nrhs==0 )  return
    !
    ! loop over nrhs
    !
    do  i = 1, nrhs
       !     
       ! forward substitution.
       !    
       b(1,i)=b(1,i)/a(1)
       jcnt = 2 
       
       do  j = 2, n
               
          k0 = j - (nz(j+1)-jcnt) + 1
               
          temp = b(j,i)
          do  k = k0, j-1
                  
             temp = temp - a(jcnt)*b(k,i)
             jcnt = jcnt + 1
                  
          end do
              
          b(j,i) = temp /a(jcnt)
          jcnt=jcnt+1 
               
       end do
       !     
       ! backward substitution.
       !     
       do  j = n, 1, -1
               
          jcnt=jcnt-1
          temp = b(j,i) / a(jcnt)
          b(j,i) = temp
               
          k0 = j - (nz(j+1)-nz(j)) + 1
          
          do  k =  j-1, k0, -1
                  
             jcnt = jcnt - 1
             b(k,i) = b(k,i) - a(jcnt)*temp
                  
          end do
               
       end do
            
    end do
    
  end subroutine chosol

  subroutine lufact( N, NNZ, NZ, A, IDIAG, ERR )
    !
    !  LU_fact computes the LU Doolittle factorization of a real
    !  matrix A stored in nonsymmetri! skyline format.
    !
    integer(ip)          :: err, n, nnz
    integer(ip)          :: nz( n+1 ), idiag( n )
    real(rp)             :: a( nnz )
    integer(ip)          :: i, j, k
    integer(ip)          :: i0, j0
    integer(ip)          :: idif, irow, iidiag, jcol, ipos
    real(rp)             :: temp,tol
    !
    ! test the input parameters.
    !
    tol  = 1.0e-14_rp
    err = 0
    if( n<0 ) then
       err = -1
    else if( nnz<0 ) then
       err = -2
    end if
    if( err.ne.0 ) then
       return
    end if
    !
    ! quick return if possible.
    !
    if( n==0.or.nnz==0 )  return
    !
    ! compute the lu factorization a = lu.
    ! initially l(1,1)=1 and u(1,1)=a(1,1) 
    ! 
    !
    ! compute first pivot
    !
    if(abs(a(1))<tol)then
       write(*,*) 'first pivot null'
       err=1
       return
    end if
    !
    ! set ipos
    !
    ipos=2
    !     
    ! outer loop on the order of the matrix
    !
    do i = 2, n
       !     
       ! compute the lower part l(i,j) 
       !
       ! column number of the first non zero in row i
       !
       jcol=i-(idiag(i)-nz(i))
       !
       ! fill row i
       !
       do j=jcol,i-1
          !
          ! row number of first non zero in column j
          !
          irow=j-(nz(j+1)-1-idiag(j))  
          !
          ! check for non zero bounds to setup the pointers for the scalar product
          !
          if(jcol>irow)then
             idif = j-jcol
             i0 = nz(i)
             j0 = nz(j+1)-idif
          else
             idif = j-irow
             i0 = ipos-idif
             j0 = idiag(j)+1 
          end if
          !
          ! compute scalar product
          !
          temp = 0.0_rp  
          do k= 1,idif
             temp = temp + a( i0 )*a( j0 )
             i0 = i0 + 1
             j0 = j0 + 1
          end do
          !
          a(ipos)=(a(ipos)-temp)/a(idiag(j))
          ipos=ipos+1 
          !
       end do
       !
       ! jump the diagonal
       !
       ipos=ipos+1
       !
       ! compute the upper part of the factorization
       !
       !
       ! row number of the first non zero in column i
       !
       irow=i-(nz(i+1)-1-idiag(i))
       !
       ! fill column i
       !
       do  j = irow,i-1
          !
          ! column number of first non zero in row j
          !
          jcol=j-(idiag(j)-nz(j))    
          !
          ! check for non zero bounds to setup the pointers
          !
          if(jcol>irow)then
             idif = j - jcol
             i0 = nz(j)
             j0 = ipos-idif
          else
             idif = j-irow
             i0 = idiag(j)-idif
             j0 = idiag(i)+1
          end if
          !
          ! compute scalar product
          !
          temp = 0.0_rp  
          do  k= 1,idif
             temp = temp + a( i0 )*a( j0 )
             i0 = i0 + 1
             j0 = j0 + 1
          end do
          
          a(ipos) = a(ipos)-temp
          ipos = ipos + 1
          
       end do
       !
       !
       ! diagonal term
       !
       !
       ! column number of first non zero in row i
       !
       jcol=i-(idiag(i)-nz(i))    
       !
       ! check for non zero bounds to setup the pointers
       !
       if(jcol>irow)then
          idif = i - jcol
          i0 = nz(i)
          j0 = ipos-idif
       else
          idif = i-irow
          i0 = idiag(i)-idif
          j0 = idiag(i)+1
       end if
       !
       ! compute scalar product
       !
       temp = 0.0_rp 
       do  k= 1,idif
          temp = temp + a( i0 )*a( j0 )
          i0 = i0 + 1
          j0 = j0 + 1
       end do
       !
       iidiag=idiag(i)
       a(iidiag) = a(iidiag)-temp
       !
       !     compute first pivot
       !
       if(abs(a(iidiag))<tol)then
          write(*,*) 'pivot null line',i
          err=i
          return
       end if
       
    end do
    
  end subroutine lufact

  subroutine lusolv( n, nnz, nz, nrhs, a, b, ldb, idiag, err )
    !
    !  LU_solve solves a system of linear equations A*X = B with a real non  symmetric
    !  matrix A stored in skyline format using the LU
    !  factorization A = LU computed by LU_fact.
    !  
    integer(ip)          :: err, ldb, n, nnz, nrhs
    integer(ip)          :: nz( n+1 ), idiag ( n )
    real(rp)             :: a( nnz ), b( ldb, nrhs )
    integer(ip)          :: j, k
    integer(ip)          :: k0
    integer(ip)          :: jcnt, iidiag, i
    real(rp)             :: temp
    !
    ! test the input parameters.
    !
    err = 0
    if( n < 0 ) then
       err = -1
    else if( nnz < 0 ) then
       err = -2
    else if( nrhs < 0 ) then
       err = -4
    else if( ldb < max( 1_ip, n ) ) then
       err = -7
    end if
    if( err /= 0 ) then
       return
    end if
    !
    ! quick return if possible.
    !
    if( n==0 .or. nnz==0 .or. nrhs==0 )   return
    !
    ! loop over nrhs
    !
    do  i = 1, nrhs
       !     
       ! forward substitution. l has 1 on its diagonal
       !     
       do  j = 2, n
          
          jcnt = nz(j)
          k0   = j - (idiag(j)-jcnt)  
               
          temp = b(j,i)
          do  k = k0, j-1
                  
             temp = temp - a(jcnt)*b(k,i) 
             jcnt = jcnt + 1
                  
          end do
          
          b(j,i) = temp
          
       end do
       !     
       ! backward substitution.
       !     
       do  j = n, 1, -1
             
          iidiag = idiag (j)
          temp   = b(j,i) / a(iidiag)
          b(j,i) = temp
          
          k0   = j - (nz(j+1)-1-iidiag) 
          jcnt = nz(j+1)
          
          do  k =  j-1, k0, -1
                  
             jcnt   = jcnt - 1
             b(k,i) = b(k,i) - a(jcnt)*temp
                  
          end do
               
       end do
            
    end do
    
  end subroutine lusolv

end module mod_skyline
!> @}
