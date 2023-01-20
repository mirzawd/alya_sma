!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_par_check_continuity

  use def_kintyp,         only : ip,rp  
  use def_master,         only : INOTMASTER,IMASTER,lun_outpu
  use def_domain,         only : npoin
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications, only : PAR_MAX
  implicit none 

  private 

  interface par_check_continuity 
     module procedure par_check_continuity_rp_i0,  &
          &           par_check_continuity_rp_i00, &
          &           par_check_continuity_rp_i1,  &
          &           par_check_continuity_rp_i2,  &
          &           par_check_continuity_rp_i3,  &
          &           par_check_continuity_rp_i31, &
          &           par_check_continuity_ip_i1, &
          &           par_check_continuity_ip_i2
  end interface par_check_continuity

  public :: par_check_continuity

contains

  subroutine par_check_continuity_rp_i00(ndofn,unkno)

    integer(ip), intent(in)    :: ndofn
    real(rp),    intent(inout) :: unkno(*)
    integer(ip)                :: ierro
    real(rp),    pointer       :: unkno_cpy(:)

    nullify(unkno_cpy)

    ierro = 0
    if( INOTMASTER ) then
       call PAR_INTERFACE_NODE_EXCHANGE(ndofn,unkno,'DIF')
       if( maxval(unkno(1:ndofn*npoin)) >= huge(1.0_rp) ) then
          ierro = 1
       end if
    end if
    call PAR_MAX(ierro)
    if( ierro > 0 ) then
       if( IMASTER ) write(6,'(a)') 'VARIABLE IS NOT CONTINUOUS'
       call runend('O.K.!')
    end if

  end subroutine par_check_continuity_rp_i00

  subroutine par_check_continuity_rp_i0(ndofn,unkno)

    integer(ip), intent(in)    :: ndofn
    real(rp),    intent(inout) :: unkno(ndofn,*)
    integer(ip)                :: ierro

    ierro = 0
    if( INOTMASTER ) then
       call PAR_INTERFACE_NODE_EXCHANGE(ndofn,unkno,'DIF')
       if( maxval(unkno(1:ndofn,1:npoin)) >= huge(1.0_rp) ) then
          ierro = 1
       end if
    end if
    call PAR_MAX(ierro)
    if( ierro > 0 ) then
       if( IMASTER ) write(6,'(a)') 'VARIABLE IS NOT CONTINUOUS'
       call runend('O.K.!')
    end if

  end subroutine par_check_continuity_rp_i0

  subroutine par_check_continuity_rp_i2(unkno)

    real(rp),    intent(inout), pointer :: unkno(:,:)
    integer(ip)                         :: ierro,ndofn

    ierro = 0
    if( INOTMASTER ) then
       ndofn = size(unkno,1)
       call PAR_INTERFACE_NODE_EXCHANGE(unkno,'DIF')
       if( maxval(unkno(1:ndofn,1:npoin)) >= huge(1.0_rp) ) then
          ierro = 1
       end if
    end if
    call PAR_MAX(ierro)
    if( ierro > 0 ) then
       if( IMASTER ) write(6,'(a)') 'VARIABLE IS NOT CONTINUOUS'
       call runend('O.K.!')
    end if

  end subroutine par_check_continuity_rp_i2

  subroutine par_check_continuity_rp_i3(unkno)

    real(rp),    intent(inout), pointer :: unkno(:,:,:)
    integer(ip)                         :: ierro,ndofn

    ierro = 0
    if( INOTMASTER ) then
       ndofn = size(unkno,1)
       call PAR_INTERFACE_NODE_EXCHANGE(unkno,'DIF')
       if( maxval(unkno(1:ndofn,1:npoin,1)) >= huge(1.0_rp) ) then
          ierro = 1
       end if
    end if
    call PAR_MAX(ierro)
    if( ierro > 0 ) then
       if( IMASTER ) write(6,'(a)') 'VARIABLE IS NOT CONTINUOUS'
       call runend('O.K.!')
    end if

  end subroutine par_check_continuity_rp_i3

  subroutine par_check_continuity_rp_i31(ndofn,unkno)

    integer(ip)                         :: ndofn
    real(rp),    intent(inout), pointer :: unkno(:,:,:)
    integer(ip)                         :: ierro
    real(rp),                   pointer :: unkno_tmp(:,:)
    
    nullify(unkno_tmp)
    ierro = 0
    if( INOTMASTER ) then
       allocate(unkno_tmp(ndofn,npoin))
       unkno_tmp(1:ndofn,1:npoin) = unkno(1:ndofn,1:npoin,1)
       call PAR_INTERFACE_NODE_EXCHANGE(ndofn,unkno_tmp,'DIF')
       unkno(:,:,1) = unkno_tmp(:,:)
       if( maxval(unkno(1:ndofn,1:npoin,1)) >= huge(1.0_rp) ) then
          ierro = 1
       end if
       deallocate(unkno_tmp)
    end if
    call PAR_MAX(ierro)
    if( ierro > 0 ) then
       if( IMASTER ) write(6,'(a)') 'VARIABLE IS NOT CONTINUOUS'
       call runend('O.K.!')
    end if

  end subroutine par_check_continuity_rp_i31

  subroutine par_check_continuity_rp_i1(unkno)

    real(rp),    intent(inout), pointer :: unkno(:)
    integer(ip)                         :: ierro,ndofn

    ierro = 0
    ndofn = size(unkno,1)
    if( INOTMASTER ) then
       call PAR_INTERFACE_NODE_EXCHANGE(ndofn,unkno,'DIF')
       if( maxval(unkno(1:npoin)) >= huge(1.0_rp) ) then
          ierro = 1
       end if
    end if
    call PAR_MAX(ierro)
    if( ierro > 0 ) then
       if( IMASTER ) write(6,'(a)') 'VARIABLE IS NOT CONTINUOUS'
       call runend('O.K.!')
    end if

  end subroutine par_check_continuity_rp_i1

  subroutine par_check_continuity_ip_i1(unkno)

    integer(ip), intent(inout), pointer :: unkno(:)
    integer(ip)                         :: ierro

    ierro = 0
    if( INOTMASTER ) then
       call PAR_INTERFACE_NODE_EXCHANGE(unkno,'DIF')
       if( maxval(unkno(1:npoin)) >= huge(1_ip) ) then
          ierro = 1
       end if
    end if
    call PAR_MAX(ierro)
    if( ierro > 0 ) then
       if( IMASTER ) write(6,'(a)') 'VARIABLE IS NOT CONTINUOUS'
       call runend('O.K.!')
    end if

  end subroutine par_check_continuity_ip_i1

  subroutine par_check_continuity_ip_i2(unkno)

    integer(ip), intent(inout), pointer :: unkno(:,:)
    integer(ip)                         :: ierro,ndofn

    ierro = 0
    if( INOTMASTER ) then
       ndofn = size(unkno,1)
       call PAR_INTERFACE_NODE_EXCHANGE(unkno,'DIF')
       if( maxval(unkno(1:ndofn,1:npoin)) >= huge(1_ip) ) then
          ierro = 1
       end if
    end if
    call PAR_MAX(ierro)
    if( ierro > 0 ) then
       if( IMASTER ) write(6,'(a)') 'VARIABLE IS NOT CONTINUOUS'
       call runend('O.K.!')
    end if

  end subroutine par_check_continuity_ip_i2

end module mod_par_check_continuity
