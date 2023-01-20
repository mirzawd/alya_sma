!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    mod_color_communicators.f90
!> @author  houzeaux
!> @date    2019-06-13
!> @brief   Module for communicators
!> @details Module for color communicators
!-----------------------------------------------------------------------

module mod_par_color_communicators

  use def_kintyp,         only : ip,rp,lg
  use def_master,         only : ISEQUEN,kfl_paral,npart
  use mod_parall,         only : PAR_WORLD_SIZE
  use mod_parall,         only : PAR_COMM_WORLD_TO_CODE_PERM
  use mod_parall,         only : I_AM_IN_COLOR
  use mod_parall,         only : PAR_COMM_COLOR
  use mod_parall,         only : PAR_CPU_TO_COLOR
  use mod_parall,         only : PAR_COLOR_TO_CPU
  use mod_parall,         only : PAR_COMM_COLOR_PERM
  use mod_parall,         only : PAR_COMM_COLOR_ARRAY
  use mod_parall,         only : PAR_COMM_MY_CODE_ARRAY
  use mod_parall,         only : mcode
  use mod_parall,         only : msubd
  use mod_parall,         only : mzone
  use mod_parall,         only : mcolo
  use mod_parall,         only : par_memor
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_COMM_NULL
  use mod_communications, only : PAR_SEND_RECEIVE
  use def_domain,         only : nzone,nsubd
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_parall,         only : par_one_sided
  use mod_parall,         only : commd
  use mod_messages,       only : messages_live
  use def_mpi
  
  implicit none

  private

  public :: par_color_communicators_allocate
  public :: par_color_communicators_deallocate
  public :: par_color_communicators

contains

  !-----------------------------------------------------------------------
  !
  !> @author  houzeaux
  !> @date    2019-06-13
  !> @brief   Allocate communicators
  !> @details Allocate color communicators
  !
  !-----------------------------------------------------------------------

  subroutine par_color_communicators_allocate()

    integer(ip) :: icolo,jcolo,m1,cz

    m1    = mcolo+1_ip
    cz    = int(PAR_WORLD_SIZE,ip)

    allocate(PAR_COMM_COLOR(0:m1-1,0:m1-1))
    PAR_COMM_COLOR = PAR_COMM_NULL
    
    call memory_alloca(par_memor,'I_AM_IN_COLOR',              'par_color_communicators',I_AM_IN_COLOR,              m1,      'INITIALIZE',0_ip)
    call memory_alloca(par_memor,'PAR_CPU_TO_COLOR',           'par_color_communicators',PAR_CPU_TO_COLOR,           cz,      'INITIALIZE',0_ip)
    call memory_alloca(par_memor,'PAR_COLOR_TO_CPU',           'par_color_communicators',PAR_COLOR_TO_CPU,           m1,      'INITIALIZE',0_ip)
    call memory_alloca(par_memor,'PAR_COMM_COLOR_PERM',        'par_color_communicators',PAR_COMM_COLOR_PERM,        m1,m1,cz,'INITIALIZE',0_ip,0_ip,0_ip)
    call memory_alloca(par_memor,'PAR_COMM_WORLD_TO_CODE_PERM','par_color_communicators',PAR_COMM_WORLD_TO_CODE_PERM,2_ip,cz, 'INITIALIZE',1_ip,0_ip)
    !
    ! In sequential
    !
    if( ISEQUEN ) then
       allocate( PAR_COMM_MY_CODE_ARRAY(1) )
       call PAR_COMM_MY_CODE_ARRAY(1) % init(COMM_NAME='COMMD')
    end if
    !
    ! Color arrays
    !
    allocate( PAR_COMM_COLOR_ARRAY(0:mcolo) )
    do icolo = 0,mcolo        
       call PAR_COMM_COLOR_ARRAY(icolo) % init(COMM_NAME='PAR_COMM_COLOR_ARRAY')
    end do
    !
    ! Color communicator
    !
    do icolo = 0,mcolo
       do jcolo = 0,mcolo
          PAR_COMM_COLOR(icolo,jcolo) = PAR_COMM_NULL
       end do
    end do

  end subroutine par_color_communicators_allocate

  !-----------------------------------------------------------------------
  !
  !> @author  houzeaux 
  !> @date    2019-06-13
  !> @brief   Allocate communicators
  !> @details Allocate color communicators
  !
  !-----------------------------------------------------------------------

  subroutine par_color_communicators_deallocate()

    integer(ip) :: icolo
    
    call memory_deallo(par_memor,'I_AM_IN_COLOR',              'par_color_communicators',I_AM_IN_COLOR)
    deallocate(PAR_COMM_COLOR)
    !call memory_deallo(par_memor,'PAR_COMM_COLOR',             'par_color_communicators',PAR_COMM_COLOR)
    call memory_deallo(par_memor,'PAR_CPU_TO_COLOR',           'par_color_communicators',PAR_CPU_TO_COLOR)
    call memory_deallo(par_memor,'PAR_COLOR_TO_CPU',           'par_color_communicators',PAR_COLOR_TO_CPU)
    call memory_deallo(par_memor,'PAR_COMM_COLOR_PERM',        'par_color_communicators',PAR_COMM_COLOR_PERM)
    call memory_deallo(par_memor,'PAR_COMM_WORLD_TO_CODE_PERM','par_color_communicators',PAR_COMM_WORLD_TO_CODE_PERM)

    do icolo = 1,mcolo                  
       call PAR_COMM_COLOR_ARRAY(icolo) % deallo(COMM_NAME='PAR_COMM_COLOR_ARRAY')
    end do
    
  end subroutine par_color_communicators_deallocate

  !----------------------------------------------------------------------
  !> @author  Guillaume Houzeaux
  !> @date    01/02/2014
  !> @brief   Max number of zones and subdomains
  !> @details Compute the max number of zones and subdomains in order
  !>          to define the mapping (code,zone,subd) => color
  !>
  !>          \verbatim
  !>
  !>          I_AM_IN_COLOR(ICOLO) .................................... if I have color ICOLO (TRUE/FALSE)
  !>          PAR_COMM_COLOR(0:MCOLO,0:MCOLO) ......................... Intercolor communicator
  !>          PAR_COMM_COLOR_ARRAY(0:MCOLO) ........................... Color communication arrays
  !>          PAR_CPU_TO_COLOR(0:PAR_WORLD_SIZE) % L(:) ............... List of colors for each world partition
  !>          PAR_COLOR_TO_CPU(0:MCOLO) % L(:) ........................ List of world partitions for each color
  !>          PAR_COMM_COLOR_PERM(0:MCOLO,0:MCOLO,0:PAR_WORLD_SIZE) ... Ranks for each communicator
  !>          PAR_COMM_WORLD_TO_CODE_PERM(2,0:PAR_WORLD_SIZE) ......... Rank permutation from world to code
  !>
  !>          \endverbatim
  !> @}
  !----------------------------------------------------------------------

  subroutine par_color_communicators()
    use def_kintyp,         only    :  ip,rp,lg
    use mod_parall,         only    :  PAR_WORLD_SIZE
    use mod_parall,         only    :  PAR_MY_CODE_RANK, PAR_COMM_CURRENT
    use mod_parall,         only    :  PAR_COMM_WORLD_TO_CODE_PERM
    use mod_parall,         only    :  I_AM_IN_COLOR
    use mod_parall,         only    :  PAR_COMM_COLOR
    use mod_parall,         only    :  PAR_CPU_TO_COLOR
    use mod_parall,         only    :  PAR_COLOR_TO_CPU
    use mod_parall,         only    :  PAR_COMM_COLOR_PERM
    use mod_parall,         only    :  PAR_COMM_MY_CODE
    use mod_parall,         only    :  PAR_COMM_WORLD
    use mod_parall,         only    :  PAR_MY_WORLD_RANK
    use mod_parall,         only    :  PAR_COMM_MY_CODE_ARRAY
    use mod_parall,         only    :  mcode
    use mod_parall,         only    :  msubd
    use mod_parall,         only    :  mzone
    use mod_parall,         only    :  mcolo
    use mod_parall,         only    :  ncolo
    use mod_parall,         only    :  par_code_zone_subd_to_color
    use mod_parall,         only    :  par_color_to_subd
    use mod_parall,         only    :  par_color_to_zone
    use mod_parall,         only    :  par_color_to_code
    use mod_parall,         only    :  par_memor
    use mod_communications, only    :  PAR_MAX
    use mod_communications, only    :  PAR_SUM
    use mod_communications, only    :  PAR_ALLGATHER
    use mod_communications, only    :  PAR_ALLGATHERV
    use mod_communications, only    :  PAR_COMM_SPLIT
    use mod_communications, only    :  PAR_SUM_ALL,PAR_BARRIER
    use mod_communications, only    :  PAR_MAX_ALL
    use mod_communications, only    :  PAR_COMM_FREE
    use mod_communications, only    :  PAR_COMM_NULL
    use def_master,         only    :  current_code
    use def_master,         only    :  ioutp,lun_outpu
    use def_master,         only    :  I_AM_IN_ZONE
    use def_master,         only    :  I_AM_IN_SUBD
    use def_domain,         only    :  nzone,nsubd
    use mod_memory,         only    :  memory_alloca
    use mod_outfor,         only    :  outfor
    use mod_messages,       only    :  livinf
    use mod_iofile,         only    :  iofile_flush_unit

    implicit none

    integer(ip)                     :: icode,izone,isubd,icolo,jcolo,ierro
    integer(ip)                     :: my_new_rank,ipart,kcolo,icolor
    integer(ip)                     :: jsubd,jzone,ii,jcode,jpart
    integer(ip)                     :: color_current_code
    integer(ip)                     :: color_world
    logical(lg)                     :: isplit

    integer(ip),         pointer    :: COMM_MATRIX(:)
    integer(ip),         pointer    :: COMM_MATRIX_TMP(:)

    integer(ip),         pointer    :: number_colors(:)
    integer(ip),         pointer    :: number_parts(:)
    integer(ip),         pointer    :: my_list_of_colors(:)
    integer(ip),         pointer    :: list_of_colors(:)
    logical(lg),         pointer    :: compute_communicator(:,:)

    integer(4)                      :: recvcount4
    integer(ip),         pointer    :: sendbuf(:)
    integer(ip),         pointer    :: recvbuf(:)

    call livinf(0_ip,'PARALL: COMPUTE INTER-COLOR COMMUNICATORS',0_ip)
    !
    ! Nullify local pointers  and initialization
    !
    nullify(COMM_MATRIX)
    nullify(COMM_MATRIX_TMP)
    nullify(number_colors)
    nullify(number_parts)
    nullify(my_list_of_colors)
    nullify(list_of_colors)
    nullify(compute_communicator)

    nullify(sendbuf)
    nullify(recvbuf)

    PAR_COMM_COLOR_PERM = 0
    !
    ! Allocate local memory
    !
    allocate( COMM_MATRIX(0:mcolo) )
    COMM_MATRIX = 0
    !
    ! Colors that shoukd not be modified
    !
    color_current_code = par_code_zone_subd_to_color(current_code,0_ip,0_ip)
    color_world        = par_code_zone_subd_to_color(0_ip,0_ip,0_ip)

    ierro = 0

    if( current_code > mcode ) ierro = 1
    !   call PAR_MAX(ierro,'IN THE UNIVERSE')
    call PAR_MAX_ALL(ierro)
    call PAR_BARRIER()  !!DMM a GH
    if( ierro == 1 ) call runend('PAR_COLOR_COMMUNICATORS: CODE NUMBER EXCEED NUMBER OF CODES')

    !----------------------------------------------------------------------
    !
    !   Using zone and subdomain information, fill in the color arrays:
    !
    !   I_AM_IN_COLOR(icolo) = .true./.false.
    !   COMM_MATRIX(icolo)   = 1/0
    !
    !   I have zone 3,5 and subd 2; 0 is all
    !   A priori we will communicate only between zones or between subdomains
    !   Structure is prepared to do combinations
    !
    !   +---+---+---+---+---+---+
    ! 0 |   |   |   | x |   | x |
    !   +---+---+---+---+---+---+
    ! 1 |   |   |   |   |   |   |
    !   +---+---+---+---+---+---+
    ! 2 | x |   |   |   |   |   |
    !   +---+---+---+---+---+---+
    ! 3 |   |   |   |   |   |   |
    !   +---+---+---+---+---+---+
    !     0   1   2   3   4   5   => zone
    !
    !----------------------------------------------------------------------

    if( PAR_MY_CODE_RANK /= 0 ) then
       !
       ! Zones and subdomains inside my code
       !
       isubd = 0
       do izone = 1,nzone
          if( I_AM_IN_ZONE(izone) ) then
             icode                = current_code
             icolo                = par_code_zone_subd_to_color(icode,izone,isubd)
             COMM_MATRIX(icolo)   = 1
          end if
       end do
       izone = 0
       do isubd = 1,nsubd
          if( I_AM_IN_SUBD(isubd) ) then
             icode                = current_code
             icolo                = par_code_zone_subd_to_color(icode,izone,isubd)
             COMM_MATRIX(icolo)   = 1
          end if
       end do
       !
       ! I am in code, of course!
       !
       icolo = par_code_zone_subd_to_color(current_code,0_ip,0_ip)
       COMM_MATRIX(icolo) = 1
       !
       ! We are obviously in COMM_WORLD tool
       !
       COMM_MATRIX(0) = 1

    end if
    !
    ! Involve the master in all communications
    !
    allocate( COMM_MATRIX_TMP(0:mcolo) )
    COMM_MATRIX_TMP = COMM_MATRIX
    call PAR_MAX(COMM_MATRIX_TMP,'IN MY CODE')
    if( PAR_MY_CODE_RANK == 0 ) COMM_MATRIX = COMM_MATRIX_TMP
    deallocate( COMM_MATRIX_TMP )

    do icolo = 0,mcolo
       if( COMM_MATRIX(icolo) == 1 ) I_AM_IN_COLOR(icolo) = .true.
    end do

    !----------------------------------------------------------------------
    !
    ! PAR_CPU_TO_COLOR(IPART) % L(:) = List of colors of partition IPART
    ! PAR_COLOR_TO_CPU(ICOLO) % L(:) = List of partitions in color ICOLO
    !
    !----------------------------------------------------------------------
    !
    ! MY_LIST_OF_COLORS(1:NCOLO) = My list of colors
    !
    ncolo = 0
    do icolo = 0,mcolo
       if( COMM_MATRIX(icolo) == 1 ) ncolo = ncolo + 1
    end do
    allocate( my_list_of_colors(ncolo) )
    ncolo = 0
    do icolo = 0,mcolo
       if( COMM_MATRIX(icolo) == 1 ) then
          ncolo = ncolo + 1
          my_list_of_colors(ncolo) = icolo
       end if
    end do
    !
    ! NUMBER_COLORS(1:PAR_WORLD_SIZE) = Number of colors of all partitions
    !
    allocate( number_colors(0:PAR_WORLD_SIZE-1) )
    call PAR_ALLGATHER(ncolo,number_colors,1_4,'IN THE WORLD')
    !
    ! LIST_OF_COLORS(1:PAR_WORLD_SIZE) = All gather list of colors of all partitions
    !
    kcolo = 0
    do ipart = 0,PAR_WORLD_SIZE-1
       kcolo = kcolo + number_colors(ipart)
    end do
    allocate( list_of_colors(kcolo) )
    call PAR_ALLGATHERV(my_list_of_colors,list_of_colors,number_colors,'IN THE WORLD')
    !
    ! PAR_CPU_TO_COLOR(IPART) % L(:) = List of colors
    !
    kcolo = 0
    do ipart = 0,PAR_WORLD_SIZE-1
       if( number_colors(ipart) > 0 ) then
          call memory_alloca(par_memor,'PAR_CPU_TO_COLOR % L','par_color_communicators',PAR_CPU_TO_COLOR(ipart) % l,number_colors(ipart))
          do icolo = 1,number_colors(ipart)
             kcolo = kcolo + 1
             PAR_CPU_TO_COLOR(ipart) % l(icolo) = list_of_colors(kcolo)
          end do
          call heapsorti1(2_ip,number_colors(ipart),PAR_CPU_TO_COLOR(ipart) % l)
       else
          nullify(  PAR_CPU_TO_COLOR(ipart) % l )
       end if
    end do

    if( associated(number_colors)     ) deallocate( number_colors )
    if( associated(list_of_colors)    ) deallocate( list_of_colors )
    if( associated(my_list_of_colors) ) deallocate( my_list_of_colors )
    !
    ! PAR_COLOR_TO_CPU(ICOLO) % L(:) = List of partitions
    !
    allocate( number_parts(0:mcolo) )
    do icolo = 0,mcolo
       number_parts(icolo) = 0
    end do
    do ipart = 0,PAR_WORLD_SIZE-1
       if( associated(PAR_CPU_TO_COLOR(ipart) % l) ) then
          do kcolo = 1,size(PAR_CPU_TO_COLOR(ipart) % l,kind=ip)
             icolo = PAR_CPU_TO_COLOR(ipart) % l(kcolo)
             number_parts(icolo) = number_parts(icolo) + 1
          end do
       end if
    end do
    do icolo = 0,mcolo
       if( number_parts(icolo) > 0 ) then
          call memory_alloca(par_memor,'PAR_COLOR_TO_CPU % L','par_color_communicators',PAR_COLOR_TO_CPU(icolo) % l,number_parts(icolo))
          number_parts(icolo) = 0
       else
          nullify( PAR_COLOR_TO_CPU(icolo) % l )
       end if
    end do
    do ipart = 0,PAR_WORLD_SIZE-1
       if( associated(PAR_CPU_TO_COLOR(ipart) % l )) then
          do kcolo = 1,size(PAR_CPU_TO_COLOR(ipart) % l,kind=ip)
             icolo = PAR_CPU_TO_COLOR(ipart) % l(kcolo)
             number_parts(icolo) = number_parts(icolo) + 1
             PAR_COLOR_TO_CPU(icolo) % l(number_parts(icolo)) = ipart
          end do
       end if
    end do
    if( associated(number_parts) ) deallocate( number_parts )
    !
    ! Up to now COMM_MATRIX is the inter-code matrix
    ! Take max over all codes to get the full one
    !
    call PAR_MAX(COMM_MATRIX,'IN THE WORLD')
    allocate( compute_communicator(0:mcolo,0:mcolo) )
    do icolo = 0,mcolo
       do jcolo = icolo,mcolo
          compute_communicator(icolo,jcolo) = .false.
          if( COMM_MATRIX(icolo) == 1 .and. COMM_MATRIX(jcolo) == 1 ) compute_communicator(icolo,jcolo) = .true.
          compute_communicator(jcolo,icolo) = compute_communicator(icolo,jcolo)
       end do
    end do
    !
    ! These two communicators should not be split, they are original ones!
    !
    compute_communicator(color_world,color_world)               = .false.
    compute_communicator(color_current_code,color_current_code) = .false.

    !do jcolo = 0,mcolo
    !   do icolo = 0,mcolo
    !      if( compute_communicator(icolo,jcolo) ) then
    !         call PAR_COMM_NULL(PAR_COMM_COLOR(icolo,jcolo))
    !      end if
    !   end do
    !end do

    !----------------------------------------------------------------------
    !
    ! Split world communicator to compute PAR_COMM_COLOR(ICOLO,JCOLO)
    ! Save permutation array PAR_COMM_COLOR_PERM(icolo,jcolo,PAR_MY_WORLD_RANK)
    !
    !----------------------------------------------------------------------
    !
    ! Inter-zone and inter-subd communicators only
    !
    do icolo = 0,mcolo
       do jcolo = icolo,mcolo
          if( compute_communicator(icolo,jcolo) ) then
             !
             ! Existing combinations: Should I be involved?
             !
             if( I_AM_IN_COLOR(icolo) .or. I_AM_IN_COLOR(jcolo) ) then
                icolor = 1
             else
                icolor = 0
             end if
             isubd  = par_color_to_subd(icolo)
             izone  = par_color_to_zone(icolo)
             icode  = par_color_to_code(icolo)
             jsubd  = par_color_to_subd(jcolo)
             jzone  = par_color_to_zone(jcolo)
             jcode  = par_color_to_code(jcolo)
             isplit = .false.
             if( izone /= 0 .and. jzone /= 0 .and. isubd == 0 .and. jsubd == 0 ) isplit = .true.
             if( isubd /= 0 .and. jsubd /= 0 .and. izone == 0 .and. jzone == 0 ) isplit = .true.
             if( icode /= 0 .and. jcode /= 0 .and. icode /= jcode .and. izone == 0 .and. jzone == 0 .and. isubd == 0 .and. jsubd == 0 ) isplit = .true.
             if( isplit ) then
                call PAR_COMM_FREE(PAR_COMM_COLOR(icolo,jcolo))
                PAR_COMM_COLOR(icolo,jcolo) =  PAR_COMM_NULL
                call PAR_COMM_SPLIT(icolor,PAR_COMM_COLOR(icolo,jcolo),my_new_rank,'IN THE WORLD')
                if( my_new_rank /= -1 ) then
                   PAR_COMM_COLOR_PERM(icolo,jcolo,PAR_MY_WORLD_RANK) = my_new_rank
                else
                   PAR_COMM_COLOR_PERM(icolo,jcolo,PAR_MY_WORLD_RANK) = -1
                end if
             else
                PAR_COMM_COLOR_PERM(icolo,jcolo,PAR_MY_WORLD_RANK) = -1
             end if
          else
             PAR_COMM_COLOR_PERM(icolo,jcolo,PAR_MY_WORLD_RANK) = -1
          end if
          PAR_COMM_COLOR(jcolo,icolo) = PAR_COMM_COLOR(icolo,jcolo)
          PAR_COMM_COLOR_PERM(jcolo,icolo,PAR_MY_WORLD_RANK) = PAR_COMM_COLOR_PERM(icolo,jcolo,PAR_MY_WORLD_RANK)
       end do
    end do
    deallocate( COMM_MATRIX )
    deallocate( compute_communicator )
    !
    ! COMM WORLD communicator
    !
    PAR_COMM_COLOR(0,0)                                = PAR_COMM_WORLD
    PAR_COMM_COLOR_PERM(0,0,PAR_MY_WORLD_RANK)         = PAR_MY_WORLD_RANK
    !
    ! This code communicator
    ! At this point, PAR_COMM_COLOR_ARRAY is not fully usable. In fact, we miss
    ! the communicator for fringe geometry
    !
    icolo                                              = color_current_code
    PAR_COMM_COLOR(icolo,icolo)                        = PAR_COMM_MY_CODE
    PAR_COMM_COLOR_PERM(icolo,icolo,PAR_MY_WORLD_RANK) = PAR_MY_CODE_RANK
    PAR_COMM_MY_CODE_ARRAY(1) % PAR_COMM_WORLD         = PAR_COMM_MY_CODE
    PAR_COMM_MY_CODE_ARRAY(1) % RANK4                  = int(PAR_MY_CODE_RANK,4)
    PAR_COMM_MY_CODE_ARRAY(1) % SIZE4                  = int(npart+1,4)
    ! 
    ! Initialization of PAR_COMM_CURRENT
    !
    PAR_COMM_CURRENT = PAR_COMM_MY_CODE
    !
    ! Get the permutation of everybody
    !
    !  call PAR_SUM(PAR_COMM_COLOR_PERM,'IN THE UNIVERSE')
    call PAR_SUM_ALL(PAR_COMM_COLOR_PERM)
    !
    ! World to code rank
    !
    allocate( recvbuf(0:PAR_WORLD_SIZE*2-1) )
    allocate( sendbuf(2) )
    sendbuf(1) = current_code
    sendbuf(2) = int(PAR_MY_CODE_RANK,ip)
    recvcount4 = 2_4
    call PAR_ALLGATHER(sendbuf,recvbuf,recvcount4,'IN THE WORLD')
    ii = -1
    do ipart = 0,PAR_WORLD_SIZE-1
       ii = ii + 1
       PAR_COMM_WORLD_TO_CODE_PERM(1,ipart) = recvbuf(ii)
       ii = ii + 1
       PAR_COMM_WORLD_TO_CODE_PERM(2,ipart) = recvbuf(ii)
    end do
    deallocate( sendbuf )
    deallocate( recvbuf )

    !----------------------------------------------------------------------
    !
    ! Check code number is not repeated over the different codes of the world
    !
    !----------------------------------------------------------------------

    ierro = 0
    if( PAR_MY_CODE_RANK == 0 ) then
       do ipart = 0,PAR_WORLD_SIZE-1
          if( ipart /= PAR_MY_WORLD_RANK ) then
             jcode = PAR_COMM_WORLD_TO_CODE_PERM(1,ipart)
             jpart = PAR_COMM_WORLD_TO_CODE_PERM(2,ipart)
             if( jpart == 0 .and. jcode == current_code ) then
                ierro = 1
             end if
          end if
       end do
    end if
    !  call PAR_MAX(ierro,'IN THE UNIVERSE')
    call PAR_MAX_ALL(ierro)
    if( ierro == 1 ) call runend('PAR_COLOR_COMMUNICATORS: WRONG CODE NUMBER')

    !----------------------------------------------------------------------
    !
    ! Output
    !
    !----------------------------------------------------------------------

    if( PAR_MY_CODE_RANK == 0 ) then
       ioutp(1) = mcode; ioutp(2) = mzone; ioutp(3) = msubd
       ioutp(4) = nzone; ioutp(5) = nsubd
       call outfor(66_ip,lun_outpu,'')
       do icolo = 0,mcolo
          do jcolo = icolo,mcolo
             if( PAR_COMM_COLOR(icolo,jcolo) /= PAR_COMM_NULL ) then
                isubd  = par_color_to_subd(icolo)
                izone  = par_color_to_zone(icolo)
                icode  = par_color_to_code(icolo)
                jsubd  = par_color_to_subd(jcolo)
                jzone  = par_color_to_zone(jcolo)
                jcode  = par_color_to_code(jcolo)
                ioutp(1) = icode; ioutp(2) = izone; ioutp(3) = isubd
                ioutp(4) = jcode; ioutp(5) = jzone; ioutp(6) = jsubd
                call outfor(67_ip,lun_outpu,'')
             end if
          end do
       end do
       call iofile_flush_unit(lun_outpu)
    end if

  end subroutine par_color_communicators

end module mod_par_color_communicators
!> @}
