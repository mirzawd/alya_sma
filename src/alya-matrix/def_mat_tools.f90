!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Maths
!> @{
!> @file    def_mat_tools.f90
!> @author  guillaume
!> @date    2021-01-26
!> @brief   Tools for matrices
!> @details Tools for matrices
!-----------------------------------------------------------------------

module def_mat_tools

  use def_kintyp_basic,      only : ip,rp,lg
  use def_mat,               only : mat
  use def_mat_dia,           only : mat_dia
  use def_mat_csr,           only : mat_csr
  use def_mat_sky,           only : mat_sky
  use def_mat_den,           only : mat_den
  use mod_optional_argument, only : optional_argument
  implicit none

  private
  real(rp), parameter :: epsil = epsilon(1.0_rp)

  public :: convert

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-26
  !> @brief   Matrix conversion
  !> @details Matrix conversion
  !> 
  !-----------------------------------------------------------------------
  
  subroutine convert(mat_in,mat_out,SYMMETRIC,TOLERANCE)

    class(*),                          intent(in)    :: mat_in
    class(*),                          intent(inout) :: mat_out
    logical(lg),    optional,          intent(in)    :: SYMMETRIC
    real(rp),       optional,          intent(in)    :: TOLERANCE
    integer(ip)                                      :: n,ndof,i,iz,j,k,l
    integer(ip)                                      :: idof,jdof,kskyl

    select type ( mat_in )
       
    class is ( mat_csr )
       !
       ! CSR => ...
       !
       select type ( mat_out )
       class is ( mat_sky ) ; call mat_out % csr2sky(mat_in,SYMMETRIC=SYMMETRIC)
       class default        ; call runend('CONVERT: CONVERSION NOT CODED')
       end select
       
    class is ( mat_sky )
       !
       ! SKY => ...
       !
       select type ( mat_out )
       class is ( mat_den ) ; call sky2den(mat_in,mat_out)
       class is ( mat_csr ) ; call sky2csr(mat_in,mat_out,TOLERANCE)
       class default        ; call runend('CONVERT: CONVERSION NOT CODED')
       end select
    end select
    
  end subroutine convert

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-26
  !> @brief   Matrix conversion
  !> @details SKY => DEN
  !> 
  !-----------------------------------------------------------------------
  
  subroutine sky2den(sky,den)
    
    class(mat_sky),           intent(in)    :: sky
    class(mat_den),           intent(inout) :: den
    integer(ip)                             :: ndof1,ndof2
    integer(ip)                             :: ndofn,nrows,ncols,jt
    integer(ip)                             :: ii,idof2,jj,kk,ll
    integer(ip)                             :: col,kskyl,idof1,irow
    integer(4)                              :: unit4
    character(150)                          :: filename_loc
    !
    ! Dimensions
    !
    ndof1 = sky % ndof1
    ndof2 = sky % ndof2
    ndofn = sky % ndof1 * sky % ndof2
    nrows = sky % nrows
    ncols = sky % nrows

    call den % init()
    den % nrows = nrows
    den % ndof1 = ndof1
    den % ndof2 = ndof2
    den % ncols = ncols
    
    call den % alloca()

    if( sky % symme ) then

       if( ndof1 /= 1 .or. ndof2 /= 1 ) stop 1
       do ii = 1,nrows
          ll  = sky % iskyl(ii+1)-sky % iskyl(ii)               
          ! column number of the first non zero in row i
          col = ii-ll
          do kk = 1,ll
             col                  = col + 1
             den % vA(1,1,ii,col) = sky % vA(sky % iskyl(ii)+kk-1)
             den % vA(1,1,col,ii) = den % vA(1,1,ii,col)
          end do
       end do

    else

       do ii = 1,nrows
          !
          ! Lower part
          !
          ! First column of row ii
          jj = ii-(sky % idiag(ii)-sky % iskyl(ii))
          do kk = 1,sky % idiag(ii)-sky % iskyl(ii)+1
             kskyl  = sky % idiag(ii)  + (jj-ii)
             den % vA(1,1,ii,jj) = sky % va(kskyl)
             jj = jj + 1
          end do
          jj = ii-(sky % iskyl(ii+1)-1-sky % idiag(ii))          
          do kskyl = sky % idiag(ii)+1,sky % iskyl(ii+1)-1
             den % vA(1,1,jj,ii) = sky % va(kskyl)
             jj = jj + 1
          end do 
       end do       

    end if

  end subroutine sky2den

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-26
  !> @brief   Matrix conversion
  !> @details SKY => CSR
  !> 
  !-----------------------------------------------------------------------
  
  subroutine sky2csr(sky,csr,TOLERANCE)

    class(mat_sky),           intent(in)    :: sky
    class(mat_csr),           intent(inout) :: csr
    real(rp),       optional, intent(in)    :: TOLERANCE
    integer(ip)                             :: ndof1,ndof2
    integer(ip)                             :: ndofn,nrows,ncols,jt
    integer(ip)                             :: ii,idof2,jj,kk,ll
    integer(ip)                             :: col,kskyl,idof1,irow
    integer(4)                              :: unit4
    character(150)                          :: filename_loc
    integer(ip),              allocatable   :: ia(:)
    real(rp)                                :: toler
    !
    ! Dimensions
    !
    nrows = sky % nrows
    if( nrows > 0 ) then
       ndof1 = sky % ndof1
       ndof2 = sky % ndof2
       ndofn = sky % ndof1 * sky % ndof2
       ncols = sky % nrows
       toler = optional_argument(epsil,TOLERANCE)

       call csr % init()
       csr % nrows = nrows
       csr % ndof1 = ndof1
       csr % ndof2 = ndof2
       csr % ncols = ncols

       if( .not. present(TOLERANCE) ) then
          if( sky % symme ) then
             csr % nz = 2*sky % nskyl-nrows
          else
             csr % nz = sky % nskyl
          end if
          call csr % alloca()
       end if
       
       allocate(ia(nrows+1))
       do ii = 1,nrows+1
          ia(ii) = 0
       end do
       
       if( sky % symme ) then

          if( ndof1 /= 1 .or. ndof2 /= 1 ) stop 1
          do ii = 1,nrows
             ll  = sky % iskyl(ii+1)-sky % iskyl(ii)               
             ! column number of the first non zero in row i
             col = ii-ll
             do kk = 1,ll
                col    = col + 1
                if( present(TOLERANCE) ) then
                   if( abs(sky % vA(sky % iskyl(ii)+kk-1))>toler ) then
                      ia(ii) = ia(ii)  + 1
                      if( ii /= col ) & 
                           ia(col) = ia(col) + 1
                   end if
                else
                   ia(ii) = ia(ii)  + 1
                   if( ii /= col ) & 
                        ia(col) = ia(col) + 1
                end if
             end do
          end do
          if( present(TOLERANCE) ) then
             csr % nz = sum(ia)
             call csr % alloca()
          end if
          
          do ii = 1,nrows+1
             csr % ia(ii) = ia(ii)
             ia(ii) = csr % ia(ii) - 1
          end do
          call csr % size_to_linked_list()
          do ii = 1,nrows+1
             ia(ii) = csr % ia(ii) -1
          end do
          do ii = 1,nrows
             ll  = sky % iskyl(ii+1)-sky % iskyl(ii)               
             ! column number of the first non zero in row i
             col = ii-ll
             do kk = 1,ll
                col                  = col + 1
                if( present(TOLERANCE) ) then
                   if( abs(sky % vA(sky % iskyl(ii)+kk-1))>toler ) then
                      ia(ii)               = ia(ii) + 1
                      csr % vA(1,1,ia(ii)) = sky % vA(sky % iskyl(ii)+kk-1)
                      csr % ja(ia(ii))     = col
                      if( ii /= col ) then
                         ia(col)               = ia(col) + 1
                         csr % ja(ia(col))     = ii
                         csr % vA(1,1,ia(col)) = sky % vA(sky % iskyl(ii)+kk-1)
                      end if
                   end if
                else                   
                   ia(ii)               = ia(ii) + 1
                   csr % vA(1,1,ia(ii)) = sky % vA(sky % iskyl(ii)+kk-1)
                   csr % ja(ia(ii))     = col
                   if( ii /= col ) then
                      ia(col)               = ia(col) + 1
                      csr % ja(ia(col))     = ii
                      csr % vA(1,1,ia(col)) = sky % vA(sky % iskyl(ii)+kk-1)
                   end if
                end if
             end do
          end do

       else

          do ii = 1,nrows
             !
             ! Lower part
             !
             ! First column of row ii
             jj = ii-(sky % idiag(ii)-sky % iskyl(ii))
             do kk = 1,sky % idiag(ii)-sky % iskyl(ii)+1
                kskyl  = sky % idiag(ii)  + (jj-ii)
                if( present(TOLERANCE) ) then
                   if( abs(sky % va(kskyl))>toler ) then
                      ia(ii) = ia(ii) + 1
                   end if
                else
                   ia(ii) = ia(ii) + 1
                end if
                jj = jj + 1
             end do
             jj = ii-(sky % iskyl(ii+1)-1-sky % idiag(ii))          
             do kskyl = sky % idiag(ii)+1,sky % iskyl(ii+1)-1
                if( present(TOLERANCE) ) then
                   if( abs(sky % va(kskyl))>toler ) then
                      ia(jj) = ia(jj) + 1
                   end if
                else
                   ia(jj) = ia(jj) + 1
                end if
                jj = jj + 1
             end do
          end do
          if( present(TOLERANCE) ) then
             csr % nz = sum(ia)
             call csr % alloca()
          end if
          do ii = 1,nrows+1
             csr % ia(ii) = ia(ii)
          end do
          call csr % size_to_linked_list()
          do ii = 1,nrows+1
             ia(ii) = csr % ia(ii) -1
          end do
          do ii = 1,nrows
             jj = ii-(sky % idiag(ii)-sky % iskyl(ii))
             do kk = 1,sky % idiag(ii)-sky % iskyl(ii)+1
                kskyl  = sky % idiag(ii)  + (jj-ii)
                if( present(TOLERANCE) ) then
                   if( abs(sky % va(kskyl))>toler ) then
                      ia(ii) = ia(ii) + 1
                      csr % ja(ia(ii)) = jj
                      csr % vA(1,1,ia(ii)) = sky % va(kskyl)
                   end if
                else
                   ia(ii) = ia(ii) + 1
                   csr % ja(ia(ii)) = jj
                   csr % vA(1,1,ia(ii)) = sky % va(kskyl)
                end if
                jj = jj + 1
             end do
             jj = ii-(sky % iskyl(ii+1)-1-sky % idiag(ii))          
             do kskyl = sky % idiag(ii)+1,sky % iskyl(ii+1)-1
                if( present(TOLERANCE) ) then
                   if( abs(sky % va(kskyl))>toler ) then
                      ia(jj) = ia(jj) + 1
                      csr % ja(ia(jj)) = ii
                      csr % vA(1,1,ia(jj)) = sky % va(kskyl)
                   end if
                else
                   ia(jj) = ia(jj) + 1
                   csr % ja(ia(jj)) = ii
                   csr % vA(1,1,ia(jj)) = sky % va(kskyl)
                end if
                jj = jj + 1
             end do

          end do

       end if
       deallocate(ia)
    end if

  end subroutine sky2csr

end module def_mat_tools
!> @}
