!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    mod_par_restart.f90
!> @author  houzeaux
!> @date    2019-10-21
!> @brief   Restart
!> @details Parallel restart
!-----------------------------------------------------------------------

module mod_par_parallel_restart 

  use def_kintyp,         only : ip,rp,lg
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_memory,         only : memory_size
  use def_master,         only : kfl_ptask
  use def_master,         only : READ_AND_RUN
  use def_master,         only : PART_AND_WRITE
  use def_master,         only : npari,nparr
  use def_master,         only : nparc,nparl
  use def_master,         only : parin
  use def_master,         only : npoi1,npoi2,npoi3
  use def_master,         only : INOTMASTER
  use def_master,         only : leinv_loc
  use def_master,         only : lninv_loc
  use def_master,         only : lbinv_loc
  use mod_communications, only : PAR_EXCHANGE
  use def_parall,         only : kfl_ascii_par
  use def_parall,         only : lun_aonlp_par
  use mod_parall,         only : par_memor
  use mod_parall,         only : commd
  use mod_parall,         only : PAR_COMM_MY_CODE_ARRAY
  use mod_domain,         only : domain_memory_allocate
  use mod_messages,       only : messages_live
  use def_domain
  implicit none

  private

  interface par_restart_read_write
     module procedure &
          par_restart_read_write_ip_1,par_restart_read_write_ip_2,par_restart_read_write_ip_3,&
          par_restart_read_write_rp_1,par_restart_read_write_rp_2,par_restart_read_write_rp_3
  end interface par_restart_read_write
  
  public :: par_parallel_restart 
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-10-21
  !> @brief   Write and read
  !> @details Write and read restart files
  !> 
  !-----------------------------------------------------------------------

  subroutine par_parallel_restart()

    integer(ip) :: ipass,ifiel,n1,n2,kfl_lgrou_dom
    logical(lg) :: if_read_data
    logical(lg) :: if_write_data

    !--------------------------------------------------------------------
    !
    ! Flags and memory
    !
    !--------------------------------------------------------------------

    if_read_data  = READ_AND_RUN()
    if_write_data = PART_AND_WRITE()
    n1            = size(kfl_field,DIM=1_ip,KIND=ip)
    n2            = size(kfl_field,DIM=2_ip,KIND=ip)
    kfl_lgrou_dom = memory_size(lgrou_dom)

    if( if_read_data ) then
       call messages_live('READ PARALLELIZATION RESTART FILES')
    else if( if_write_data ) then
       call messages_live('WRITE PARALLELIZATION RESTART FILES')
    end if
    
    !--------------------------------------------------------------------
    !
    ! Dimensions
    !
    !--------------------------------------------------------------------
    
    do ipass = 1,2 
       npari = 0 ; nparr = 0 ; nparc = 0 ; nparl = 0
       !
       ! Dimensions
       !
#ifdef NDIMEPAR
#else
       call PAR_EXCHANGE(ndime,             parin,npari,ipass,READ_DATA=if_read_data,WRITE_DATA=if_write_data)
#endif       
       call PAR_EXCHANGE(kfl_autbo,         parin,npari,ipass,READ_DATA=if_read_data,WRITE_DATA=if_write_data)
       call PAR_EXCHANGE(npoin,             parin,npari,ipass,READ_DATA=if_read_data,WRITE_DATA=if_write_data)
       call PAR_EXCHANGE(nelem,             parin,npari,ipass,READ_DATA=if_read_data,WRITE_DATA=if_write_data)
       call PAR_EXCHANGE(nboun,             parin,npari,ipass,READ_DATA=if_read_data,WRITE_DATA=if_write_data)
       call PAR_EXCHANGE(nperi,             parin,npari,ipass,READ_DATA=if_read_data,WRITE_DATA=if_write_data)
       call PAR_EXCHANGE(nhang,             parin,npari,ipass,READ_DATA=if_read_data,WRITE_DATA=if_write_data)
       call PAR_EXCHANGE(nfiel,             parin,npari,ipass,READ_DATA=if_read_data,WRITE_DATA=if_write_data)
       call PAR_EXCHANGE(n1,n2,kfl_field,   parin,npari,ipass,READ_DATA=if_read_data,WRITE_DATA=if_write_data)
       call PAR_EXCHANGE(mcodb,             parin,npari,ipass,READ_DATA=if_read_data,WRITE_DATA=if_write_data)
#ifndef PNODE_VALUE
       call PAR_EXCHANGE(mnode,             parin,npari,ipass,READ_DATA=if_read_data,WRITE_DATA=if_write_data)
#endif
       !
       ! Read in reageo
       !       
       continue
       ! 
       ! Read in reaset
       !
       call PAR_EXCHANGE(neset,             parin,npari,ipass,READ_DATA=if_read_data,WRITE_DATA=if_write_data)
       call PAR_EXCHANGE(nbset,             parin,npari,ipass,READ_DATA=if_read_data,WRITE_DATA=if_write_data)
       call PAR_EXCHANGE(nnset,             parin,npari,ipass,READ_DATA=if_read_data,WRITE_DATA=if_write_data)
       call PAR_EXCHANGE(neset_origi,       parin,npari,ipass,READ_DATA=if_read_data,WRITE_DATA=if_write_data)
       call PAR_EXCHANGE(nbset_origi,       parin,npari,ipass,READ_DATA=if_read_data,WRITE_DATA=if_write_data)
       call PAR_EXCHANGE(nnset_origi,       parin,npari,ipass,READ_DATA=if_read_data,WRITE_DATA=if_write_data)
       call PAR_EXCHANGE(kfl_lgrou_dom,     parin,npari,ipass,READ_DATA=if_read_data,WRITE_DATA=if_write_data)
       !
       ! Read in reabcs
       !
       call PAR_EXCHANGE(kfl_icodn,         parin,npari,ipass,READ_DATA=if_read_data,WRITE_DATA=if_write_data)
       call PAR_EXCHANGE(kfl_icodb,         parin,npari,ipass,READ_DATA=if_read_data,WRITE_DATA=if_write_data)
       !
       ! Read in reafie
       !
       continue
       !
       ! Communication dimensions
       !
       call PAR_EXCHANGE(npoi1,             parin,npari,ipass,READ_DATA=if_read_data,WRITE_DATA=if_write_data)    
       call PAR_EXCHANGE(npoi2,             parin,npari,ipass,READ_DATA=if_read_data,WRITE_DATA=if_write_data)    
       call PAR_EXCHANGE(npoi3,             parin,npari,ipass,READ_DATA=if_read_data,WRITE_DATA=if_write_data)    
       call PAR_EXCHANGE(commd % nneig,     parin,npari,ipass,READ_DATA=if_read_data,WRITE_DATA=if_write_data)    
       call PAR_EXCHANGE(commd % bound_dim, parin,npari,ipass,READ_DATA=if_read_data,WRITE_DATA=if_write_data)

       if( ipass == 1 ) then
          call memory_alloca(par_memor,'PARIN','par_sendat',parin,npari)
          if( READ_AND_RUN() ) then
             call par_restart_read_write(parin)
          end if
       end if

    end do

    if( PART_AND_WRITE() ) then
       call par_restart_read_write(parin)
    end if

    call memory_deallo(par_memor,'PARIN','par_sendat',parin)

    !--------------------------------------------------------------------
    !
    ! Communication arrays
    !
    !--------------------------------------------------------------------

    if( READ_AND_RUN() ) then
       call PAR_COMM_MY_CODE_ARRAY(1) % alloca(par_memor)
       !call PAR_ALLOCATE_COMMUNICATION_ARRAY(PAR_COMM_MY_CODE_ARRAY(1),par_memor)!,COMM_NAME='COMMD')
    end if

    do ipass = 1,2 
       npari = 0 ; nparr = 0 ; nparc = 0 ; nparl = 0
       
       call PAR_EXCHANGE(commd % neights,   parin,npari,ipass,READ_DATA=if_read_data,WRITE_DATA=if_write_data)
       call PAR_EXCHANGE(commd % bound_size,parin,npari,ipass,READ_DATA=if_read_data,WRITE_DATA=if_write_data)
       call PAR_EXCHANGE(commd % bound_perm,parin,npari,ipass,READ_DATA=if_read_data,WRITE_DATA=if_write_data)
       call PAR_EXCHANGE(commd % bound_invp,parin,npari,ipass,READ_DATA=if_read_data,WRITE_DATA=if_write_data)
 
       if( ipass == 1 ) then 
          call memory_alloca(par_memor,'PARIN','par_sendat',parin,npari)
          if( READ_AND_RUN() ) then
             call par_restart_read_write(parin)
          end if
       end if
    end do
    if( PART_AND_WRITE() ) then
       call par_restart_read_write(parin)
    end if

    call memory_deallo(par_memor,'PARIN','par_sendat',parin)
    
    !--------------------------------------------------------------------
    !
    ! Geometrical arrays
    !
    !--------------------------------------------------------------------

    if( INOTMASTER ) then
       if( READ_AND_RUN() ) then
          call domain_memory_allocate('GEOMETRY' )   
          if( kfl_lgrou_dom > 0 ) call domain_memory_allocate('LGROU_DOM')  
          call domain_memory_allocate('LEINV_LOC')  
          call domain_memory_allocate('LNINV_LOC')  
          call domain_memory_allocate('LBINV_LOC')  
          call domain_memory_allocate('LESET'    )      
          call domain_memory_allocate('LBSET'    )      
          call domain_memory_allocate('LNSET'    )      
          call domain_memory_allocate('KFL_CODNO')  
          call domain_memory_allocate('KFL_CODBO')
          call domain_memory_allocate('XFIEL % A')
       end if
       !
       ! Geometry
       !
       call par_restart_read_write(ltype     )    
       call par_restart_read_write(lelch     )         
       call par_restart_read_write(lnods     )         
       call par_restart_read_write(lesub     )         
       call par_restart_read_write(lmate     )         
       call par_restart_read_write(leinv_loc )

       call par_restart_read_write(coord     )     
       call par_restart_read_write(lnoch     )     
       call par_restart_read_write(lmast     )     
       call par_restart_read_write(lninv_loc ) 

       call par_restart_read_write(lnodb     )     
       call par_restart_read_write(ltypb     )     
       call par_restart_read_write(lboch     )     
       call par_restart_read_write(lelbo     )     
       call par_restart_read_write(lbinv_loc )
       !
       ! Others
       !
       if( kfl_lgrou_dom > 0 ) call par_restart_read_write(lgrou_dom )     
       call par_restart_read_write(leset     )     
       call par_restart_read_write(lbset     )     
       call par_restart_read_write(lnset     )     
       call par_restart_read_write(kfl_codno )     
       call par_restart_read_write(kfl_codbo )
       do ifiel = 1,nfiel
          call par_restart_read_write(xfiel(ifiel) % a )     
       end do
    end if
    
    !--------------------------------------------------------------------
    !
    ! We are now in normal run... just in case we do repartitioning
    !
    !--------------------------------------------------------------------

    kfl_ptask = 1
    
  end subroutine par_parallel_restart

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-10-18
  !> @brief   Write and read
  !> @details Write and read arrays
  !> 
  !-----------------------------------------------------------------------

  subroutine par_restart_read_write_ip_1(xx)
  
    integer(ip), pointer, intent(inout) :: xx(:)
    integer(ip)                      :: nn
    integer(ip)                      :: xx_tmp(2)
    
    nn = memory_size(xx)
    if( nn > 0 ) then
       call par_restart_read_write_go_ip(nn,xx)
    else
       call par_restart_read_write_go_ip(nn,xx_tmp)
    end if
    
  end subroutine par_restart_read_write_ip_1
    
  subroutine par_restart_read_write_ip_2(xx)
  
    integer(ip), pointer, intent(inout) :: xx(:,:)
    integer(ip)                      :: nn
    integer(ip)                      :: xx_tmp(2)
    
    nn = memory_size(xx)
    if( nn > 0 ) then
       call par_restart_read_write_go_ip(nn,xx)
    else
       call par_restart_read_write_go_ip(nn,xx_tmp)
    end if
    
  end subroutine par_restart_read_write_ip_2
    
  subroutine par_restart_read_write_ip_3(xx)
  
    integer(ip), pointer, intent(inout) :: xx(:,:,:)
    integer(ip)                      :: nn
    integer(ip)                      :: xx_tmp(2)
    
    nn = memory_size(xx)
    if( nn > 0 ) then
       call par_restart_read_write_go_ip(nn,xx)
    else
       call par_restart_read_write_go_ip(nn,xx_tmp)
    end if
    
  end subroutine par_restart_read_write_ip_3
    
  subroutine par_restart_read_write_rp_1(xx)
  
    real(rp), pointer, intent(inout) :: xx(:)
    integer(ip)                   :: nn
    real(rp)                      :: xx_tmp(2)
    
    nn = memory_size(xx)
    if( nn > 0 ) then
       call par_restart_read_write_go_rp(nn,xx)
    else
       call par_restart_read_write_go_rp(nn,xx_tmp)
    end if
    
  end subroutine par_restart_read_write_rp_1
    
  subroutine par_restart_read_write_rp_2(xx)
  
    real(rp), pointer, intent(inout) :: xx(:,:)
    integer(ip)                   :: nn
    real(rp)                      :: xx_tmp(2)
    
    nn = memory_size(xx)
    if( nn > 0 ) then
       call par_restart_read_write_go_rp(nn,xx)
    else
       call par_restart_read_write_go_rp(nn,xx_tmp)
    end if
    
  end subroutine par_restart_read_write_rp_2
    
  subroutine par_restart_read_write_rp_3(xx)
  
    real(rp), pointer, intent(inout) :: xx(:,:,:)
    integer(ip)                   :: nn
    real(rp)                      :: xx_tmp(2)
    
    nn = memory_size(xx)
    if( nn > 0 ) then
       call par_restart_read_write_go_rp(nn,xx)
    else
       call par_restart_read_write_go_rp(nn,xx_tmp)
    end if
    
  end subroutine par_restart_read_write_rp_3
    
  subroutine par_restart_read_write_go_ip(nn,xx)

    integer(ip), intent(inout) :: nn
    integer(ip), intent(inout) :: xx(*)
    integer(ip)                :: ii
    integer(4)                 :: iunit4

    iunit4 = int(lun_aonlp_par,KIND=4)
 
    if( kfl_ascii_par == 0 ) then
       if(    PART_AND_WRITE() ) then
          write(iunit4) nn
          if( nn > 0 ) write(iunit4) (xx(ii),ii=1,nn)
       else if( READ_AND_RUN() ) then
          read(iunit4) nn
          if( nn > 0 ) read(iunit4) (xx(ii),ii=1,nn)       
       end if
    else
       if(    PART_AND_WRITE() ) then
          write(iunit4,*) nn
          if( nn > 0 ) write(iunit4,*) (xx(ii),ii=1,nn)
       else if( READ_AND_RUN() ) then
          read(iunit4,*) nn
          if( nn > 0 ) read(iunit4,*) (xx(ii),ii=1,nn)       
       end if
    end if
    
  end subroutine par_restart_read_write_go_ip
  
  subroutine par_restart_read_write_go_rp(nn,xx)

    integer(ip), intent(inout) :: nn
    real(rp),    intent(inout) :: xx(*)
    integer(ip)                :: ii
    integer(4)                 :: iunit4

    iunit4 = int(lun_aonlp_par,KIND=4)
    
    if( kfl_ascii_par == 0 ) then
       if(    PART_AND_WRITE() ) then
          write(iunit4) nn
          if( nn > 0 ) write(iunit4) (xx(ii),ii=1,nn)
       else if( READ_AND_RUN() ) then
          read(iunit4) nn
          if( nn > 0 ) read(iunit4) (xx(ii),ii=1,nn)
       end if
    else
       if(    PART_AND_WRITE() ) then
          write(iunit4,*) nn
          if( nn > 0 ) write(iunit4,*) (xx(ii),ii=1,nn)
       else if( READ_AND_RUN() ) then
          read(iunit4,*) nn
          if( nn > 0 ) read(iunit4,*) (xx(ii),ii=1,nn)
       end if
    end if
    
  end subroutine par_restart_read_write_go_rp
  
end module mod_par_parallel_restart
!> @}
