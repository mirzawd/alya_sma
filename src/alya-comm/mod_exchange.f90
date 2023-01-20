!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup IO
!> @{
!> @file    mod_exchange.f90
!> @author  houzeaux and Eduardo Perez
!> @date    2020-05-20
!> @brief   Exchange
!> @details Some tools for exchange
!-----------------------------------------------------------------------

module mod_exchange

  use def_kintyp_basic,   only : ip,rp,lg
  use mod_communications, only : PAR_BROADCAST
  implicit none
  private

  type ptr_inte_pack
     integer(ip), pointer :: p
  end type ptr_inte_pack
  
  type ptr_real_pack
     real(rp),    pointer :: p
  end type ptr_real_pack

  type ptr_logi_pack
     logical(lg), pointer :: p
  end type ptr_logi_pack

  type ptr_char_pack
     character(len=:), pointer :: p
  end type ptr_char_pack

  interface exchange_add
     module procedure &
          exchange_add_ip_s,&
          exchange_add_ip_1,&
          exchange_add_ip_2,&
          exchange_add_ip_3,&
          exchange_add_rp_s,&
          exchange_add_rp_1,&
          exchange_add_rp_2,&
          exchange_add_rp_3,&
          exchange_add_rp_4,&
          exchange_add_lg_s,&
          exchange_add_lg_1,&
          exchange_add_ch_s,&
          exchange_add_ch_1
  end interface exchange_add

  type(ptr_inte_pack), allocatable :: inte_pack(:)
  type(ptr_real_pack), allocatable :: real_pack(:)
  type(ptr_logi_pack), allocatable :: logi_pack(:)
  type(ptr_char_pack), allocatable :: char_pack(:)
  
  integer(ip)                      :: inte_size
  integer(ip)                      :: real_size
  integer(ip)                      :: logi_size
  integer(ip)                      :: char_size

  public :: exchange_init
  public :: exchange_add
  public :: exchange_end
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-20
  !> @brief   Initialize exchange
  !> @details Initialize exchange
  !> 
  !-----------------------------------------------------------------------

  subroutine exchange_init(INTE_PACK_SIZE,REAL_PACK_SIZE,LOGI_PACK_SIZE)

    integer(ip), optional, intent(in) :: INTE_PACK_SIZE
    integer(ip), optional, intent(in) :: REAL_PACK_SIZE
    integer(ip), optional, intent(in) :: LOGI_PACK_SIZE

    inte_size = 0
    real_size = 0
    logi_size = 0
    char_size = 0

    if( allocated(inte_pack) ) deallocate(inte_pack)
    if( allocated(real_pack) ) deallocate(real_pack)
    if( allocated(logi_pack) ) deallocate(logi_pack)
    if( allocated(char_pack) ) deallocate(char_pack)

    if( present(INTE_PACK_SIZE) ) allocate(inte_pack(INTE_PACK_SIZE))
    if( present(REAL_PACK_SIZE) ) allocate(inte_pack(REAL_PACK_SIZE))
    if( present(LOGI_PACK_SIZE) ) allocate(logi_pack(LOGI_PACK_SIZE))

  end subroutine exchange_init
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-20
  !> @brief   Finalize the task
  !> @details Finalize the task by writing or broadcasting
  !> 
  !-----------------------------------------------------------------------

  subroutine exchange_end()

    integer(ip)                    :: ii,nn
    integer(ip),      allocatable  :: buffi(:)
    real(rp),         allocatable  :: buffr(:)
    logical(lg),      allocatable  :: buffl(:)
    character(len=:), allocatable  :: buffc

    if( allocated(inte_pack) ) then
       allocate(buffi(inte_size))
       do ii = 1,inte_size
          buffi(ii) = inte_pack(ii) % p
       end do
       call PAR_BROADCAST(inte_size,buffi)
       do ii = 1,inte_size
          inte_pack(ii) % p = buffi(ii) 
       end do
       deallocate(buffi)           
    end if

    if( allocated(real_pack) ) then
       allocate(buffr(real_size))
       do ii = 1,real_size
          buffr(ii) = real_pack(ii) % p
       end do
       call PAR_BROADCAST(real_size,buffr)
       do ii = 1,real_size
          real_pack(ii) % p = buffr(ii) 
       end do
       deallocate(buffr)           
    end if

    if( allocated(logi_pack) ) then
       allocate(buffl(logi_size))
       do ii = 1,logi_size
          buffl(ii) = logi_pack(ii) % p
       end do
       call PAR_BROADCAST(logi_size,buffl)
       do ii = 1,logi_size
          logi_pack(ii) % p = buffl(ii) 
       end do
       deallocate(buffl)           
    end if

    if( allocated(char_pack) ) then
       do ii = 1,char_size
          nn    = len(char_pack(ii) % p)
          buffc = char_pack(ii) % p
          call PAR_BROADCAST(nn,buffc,'IN MY CODE')
          char_pack(ii) % p = buffc
          deallocate(buffc)           
       end do
    end if
    !
    ! Reinitialize and deallocate
    !
    call exchange_init()

  end subroutine exchange_end
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-20
  !> @brief   Pack integer variables
  !> @details Pack integer variables
  !> 
  !-----------------------------------------------------------------------
  
  subroutine exchange_add_ip_s(val)
    type (ptr_inte_pack), allocatable :: ptr_tmp(:)
    integer(ip),          target      :: val
    
    if( .not. allocated(inte_pack) )  then
       allocate(inte_pack(0))
       inte_size = 0
    end if
    inte_size = inte_size + 1

    if( inte_size > size(inte_pack) ) then
       allocate(ptr_tmp(2_ip*inte_size))
       if( inte_size >= 1 ) ptr_tmp(1:inte_size-1) = inte_pack(1:inte_size-1)
       call move_alloc(from=ptr_tmp, to=inte_pack)
    end if
    
    inte_pack(inte_size) % p => val

  end subroutine exchange_add_ip_s

  subroutine exchange_add_ip_1(val)
    type (ptr_inte_pack), allocatable :: ptr_tmp(:)
    integer(ip),          target      :: val(:)
    integer(ip)                       :: kk,ii
    
    kk = size(val)
    if( .not. allocated(inte_pack) )  then
       allocate(inte_pack(0))
       inte_size = 0
    end if
    inte_size = inte_size + kk

    if( inte_size > size(inte_pack) ) then   
       allocate(ptr_tmp(2_ip*inte_size))
       if( inte_size >= 1 ) ptr_tmp(1:inte_size-kk) = inte_pack(1:inte_size-kk)
       call move_alloc(from=ptr_tmp, to=inte_pack)
    end if
    
    do ii = 1,kk
       inte_pack(inte_size-kk+ii) % p => val(ii)
    end do
    
  end subroutine exchange_add_ip_1
  
  subroutine exchange_add_ip_2(val)
    type (ptr_inte_pack), allocatable :: ptr_tmp(:)
    integer(ip),          target      :: val(:,:)
    integer(ip)                       :: ii,jj,kk,ll

    kk = size(val)
    if( .not. allocated(inte_pack) )  then
       allocate(inte_pack(0))
       inte_size = 0
    end if
    inte_size = inte_size + kk

    if( inte_size > size(inte_pack) ) then   
       allocate(ptr_tmp(2_ip*inte_size))
       if( inte_size >= 1 ) ptr_tmp(1:inte_size-kk) = inte_pack(1:inte_size-kk)
       call move_alloc(from=ptr_tmp, to=inte_pack)
    end if

    ll = 0
    do ii = 1,size(val,2)
       do jj = 1,size(val,1)
          ll = ll + 1
          inte_pack(inte_size-kk+ll) % p => val(jj,ii)
       end do
    end do
    
  end subroutine exchange_add_ip_2
  
  subroutine exchange_add_ip_3(val)
    type (ptr_inte_pack), allocatable :: ptr_tmp(:)
    integer(ip),          target      :: val(:,:,:)
    integer(ip)                       :: ii,jj,kk,ll,mm

    kk = size(val)
    if( .not. allocated(inte_pack) )  then
       allocate(inte_pack(0))
       inte_size = 0
    end if
    inte_size = inte_size + kk

    if( inte_size > size(inte_pack) ) then   
       allocate(ptr_tmp(2_ip*inte_size))
       if( inte_size >= 1 ) ptr_tmp(1:inte_size-kk) = inte_pack(1:inte_size-kk)
       call move_alloc(from=ptr_tmp, to=inte_pack)
    end if

    ll = 0
    do mm = 1,size(val,3)
       do ii = 1,size(val,2)
          do jj = 1,size(val,1)
             ll = ll + 1
             inte_pack(inte_size-kk+ll) % p => val(jj,ii,mm)
          end do
       end do
    end do
    
  end subroutine exchange_add_ip_3
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-20
  !> @brief   Pack real variables
  !> @details Pack real variables
  !> 
  !-----------------------------------------------------------------------
 
  subroutine exchange_add_rp_s(val)
    type (ptr_real_pack), allocatable :: ptr_tmp(:)
    real(rp),             target      :: val
    
    if( .not. allocated(real_pack) )  then
       allocate(real_pack(0))
       real_size = 0
    end if
    real_size = real_size + 1

    if( real_size > size(real_pack) ) then 
       allocate(ptr_tmp(2_ip*real_size))
       if( real_size >= 1 ) ptr_tmp(1:real_size-1) = real_pack(1:real_size-1)
       call move_alloc(from=ptr_tmp, to=real_pack)
    end if
    
    real_pack(real_size) % p => val

  end subroutine exchange_add_rp_s

  subroutine exchange_add_rp_1(val)
    type (ptr_real_pack), allocatable :: ptr_tmp(:)
    real(rp),             target      :: val(:)
    integer(ip)                       :: kk,ii

    kk = size(val)
    if( .not. allocated(real_pack) )  then
       allocate(real_pack(0))
       real_size = 0
    end if
    real_size = real_size + kk
    
    if( real_size > size(real_pack) ) then 
       allocate(ptr_tmp(2_ip*real_size))
       if( real_size >= 1 ) ptr_tmp(1:real_size-kk) = real_pack(1:real_size-kk)
       call move_alloc(from=ptr_tmp, to=real_pack)
    end if

    do ii = 1,kk
      real_pack(real_size-kk+ii) % p => val(ii)
    end do
    
  end subroutine exchange_add_rp_1

  subroutine exchange_add_rp_2(val)
    type (ptr_real_pack), allocatable :: ptr_tmp(:)
    real(rp),             target      :: val(:,:)
    integer(ip)                       :: ii,jj,kk,ll

    kk = size(val)
    if( .not. allocated(real_pack) )  then
       allocate(real_pack(0))
       real_size = 0
    end if
    real_size = real_size + kk

    if( real_size > size(real_pack) ) then   
       allocate(ptr_tmp(2_ip*real_size))
       if( real_size >= 1 ) ptr_tmp(1:real_size-kk) = real_pack(1:real_size-kk)
       call move_alloc(from=ptr_tmp, to=real_pack)
    end if

    ll = 0
    do ii = 1,size(val,2)
       do jj = 1,size(val,1)
          ll = ll + 1
          real_pack(real_size-kk+ll) % p => val(jj,ii)
       end do
    end do
    
  end subroutine exchange_add_rp_2
    
  subroutine exchange_add_rp_3(val)
    type (ptr_real_pack), allocatable :: ptr_tmp(:)
    real(rp),             target      :: val(:,:,:)
    integer(ip)                       :: ii,jj,kk,ll,mm

    kk = size(val)
    if( .not. allocated(real_pack) )  then
       allocate(real_pack(0))
       real_size = 0
    end if
    real_size = real_size + kk

    if( real_size > size(real_pack) ) then   
       allocate(ptr_tmp(2_ip*real_size))
       if( real_size >= 1 ) ptr_tmp(1:real_size-kk) = real_pack(1:real_size-kk)
       call move_alloc(from=ptr_tmp, to=real_pack)
    end if

    ll = 0
    do mm = 1,size(val,3)
       do ii = 1,size(val,2)
          do jj = 1,size(val,1)
             ll = ll + 1
             real_pack(real_size-kk+ll) % p => val(jj,ii,mm)
          end do
       end do
    end do
    
  end subroutine exchange_add_rp_3
    
  subroutine exchange_add_rp_4(val)
    type (ptr_real_pack), allocatable :: ptr_tmp(:)
    real(rp),             target      :: val(:,:,:,:)
    integer(ip)                       :: ii,jj,kk,ll,mm,nn

    kk = size(val)
    if( .not. allocated(real_pack) )  then
       allocate(real_pack(0))
       real_size = 0
    end if
    real_size = real_size + kk

    if( real_size > size(real_pack) ) then   
       allocate(ptr_tmp(2_ip*real_size))
       if( real_size >= 1 ) ptr_tmp(1:real_size-kk) = real_pack(1:real_size-kk)
       call move_alloc(from=ptr_tmp, to=real_pack)
    end if

    ll = 0
    do nn = 1,size(val,4)
       do mm = 1,size(val,3)
          do ii = 1,size(val,2)
             do jj = 1,size(val,1)
                ll = ll + 1
                real_pack(real_size-kk+ll) % p => val(jj,ii,mm,nn)
             end do
          end do
       end do
    end do
    
  end subroutine exchange_add_rp_4
    
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-20
  !> @brief   Pack logical
  !> @details Pack logical
  !> 
  !-----------------------------------------------------------------------

  subroutine exchange_add_lg_s(val)
    type (ptr_logi_pack), allocatable :: ptr_tmp(:)
    logical(lg),          target      :: val
    
    if( .not. allocated(logi_pack) )  then
       allocate(logi_pack(0))
       logi_size = 0
    end if
    logi_size = logi_size + 1

    if( logi_size > size(logi_pack) ) then
       allocate(ptr_tmp(2_ip*logi_size))
       if( logi_size >= 1 ) ptr_tmp(1:logi_size-1) = logi_pack(1:logi_size-1)
       call move_alloc(from=ptr_tmp, to=logi_pack)
    end if
    
    logi_pack(logi_size) % p => val

  end subroutine exchange_add_lg_s

  subroutine exchange_add_lg_1(val)
    type (ptr_logi_pack), allocatable :: ptr_tmp(:)
    logical(lg),          target      :: val(:)
    integer(ip)                       :: kk,ii

    kk = size(val)
    if( .not. allocated(logi_pack) )  then
       allocate(logi_pack(0))
       logi_size = 0
    end if
    logi_size = logi_size + kk

    if( logi_size > size(logi_pack) ) then   
       allocate(ptr_tmp(2_ip*logi_size))
       if( logi_size >= 1 ) ptr_tmp(1:logi_size-kk) = logi_pack(1:logi_size-kk)
       call move_alloc(from=ptr_tmp, to=logi_pack)
    end if
    
    do ii = 1,kk
       logi_pack(logi_size-kk+ii) % p => val(ii)
    end do
    
  end subroutine exchange_add_lg_1

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-20
  !> @brief   Pack logical
  !> @details Pack logical
  !> 
  !-----------------------------------------------------------------------

  subroutine exchange_add_ch_s(val)
    type (ptr_char_pack), allocatable :: ptr_tmp(:)
    character(LEN=*),     target      :: val

    if( .not. allocated(char_pack) )  then
       allocate(char_pack(0))
       char_size = 0
    end if
    char_size = char_size + 1

    if( char_size > size(char_pack) ) then
       allocate(ptr_tmp(2_ip*char_size))
       if( char_size >= 1 ) ptr_tmp(1:char_size-1) = char_pack(1:char_size-1)
       call move_alloc(from=ptr_tmp, to=char_pack)
    end if
    
    char_pack(char_size) % p => val

  end subroutine exchange_add_ch_s

  subroutine exchange_add_ch_1(val)
    type (ptr_char_pack), allocatable :: ptr_tmp(:)
    character(LEN=*),     target      :: val(:)
    integer(ip)                       :: kk,ii

    kk = size(val)
    if( .not. allocated(char_pack) )  then
       allocate(char_pack(0))
       char_size = 0
    end if
    char_size = char_size + kk
    
    if( char_size > size(char_pack) ) then 
       allocate(ptr_tmp(2_ip*char_size))
       if( char_size >= 1 ) ptr_tmp(1:char_size-kk) = char_pack(1:char_size-kk)
       call move_alloc(from=ptr_tmp, to=char_pack)
    end if

    do ii = 1,kk
      char_pack(char_size-kk+ii) % p => val(ii)
    end do
    
  end subroutine exchange_add_ch_1
  
end module mod_exchange
!> @}
