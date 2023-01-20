!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_opebcs
  !-----------------------------------------------------------------------
  !****f* outrut/mod_opebcs
  ! NAME
  !   mod_opebcs
  ! DESCRIPTION
  !   This routine manages the opebcsocess
  ! USES
  ! USED BY
  !   output_***
  !   outvar_***
  !***
  !-----------------------------------------------------------------------

  use def_master,            only : rp,ip,lg
  use def_master,            only : IMASTER
  use def_master,            only : IPARALL
  use def_master,            only : bc_nodes
  use def_master,            only : bc_bound
  use def_master,            only : coutp
  use def_master,            only : mem_modul
  use def_master,            only : modul
  use def_domain,            only : mcodb
  use def_domain,            only : mcono
  use def_domain,            only : kfl_icodb
  use def_domain,            only : kfl_icodn
  use def_domain,            only : tncod
  use def_domain,            only : tbcod
  use mod_memory,            only : memory_size
  use mod_outfor,            only : outfor
  use mod_exchange,          only : exchange_init
  use mod_exchange,          only : exchange_add
  use mod_exchange,          only : exchange_end
  use mod_optional_argument, only : optional_argument
  use mod_codes,             only : codes_read_on_boundaries
  use mod_codes,             only : codes_read_on_nodes
  implicit none

  private

  integer(ip) :: mcodc

  interface opebcs_initialization_structure
     module procedure     opebcs_initialization_structure_boundary,&
          &               opebcs_initialization_structure_node
  end interface opebcs_initialization_structure

  interface opebcs_initialization_variable
     module procedure     opebcs_initialization_variable_boundary,&
          &               opebcs_initialization_variable_boundary_s,&
          &               opebcs_initialization_variable_node_s,&
          &               opebcs_initialization_variable_node_1
  end interface opebcs_initialization_variable

  interface boundary_conditions_impose_edge_codes
     module procedure     boundary_conditions_impose_edge_codes_1,&
          &               boundary_conditions_impose_edge_codes_2
  end interface boundary_conditions_impose_edge_codes
  
  interface boundary_conditions_exchange
     module procedure     spbbcs,&
          &               spnbcs
  end interface boundary_conditions_exchange
  
  public :: boundary_conditions_read_node_codes  
  public :: boundary_conditions_read_boundary_codes   
  public :: boundary_conditions_impose_node_codes     ! Impose conditions on nodes 
  public :: boundary_conditions_impose_edge_codes     ! Impose conditions on edge 
  public :: boundary_conditions_number_codes          ! Return number of declared codes
  public :: boundary_conditions_exchange              ! Exchange boundary conditions
  public :: opebcs_initialization_structure           ! Allocate bc structure
  public :: opebcs_initialization_variable            ! Initialize bc structure
  public :: opbbcs
  public :: opnbcs
  public :: spbbcs
  public :: spnbcs
  public :: cpybcs
  public :: cpybcs_boundaries
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-11-09
  !> @brief   Number of boundary codes
  !> @details Number of boundary codes
  !> 
  !-----------------------------------------------------------------------
  
  integer(ip) function boundary_conditions_number_codes(txcod)

    class(*), intent(inout) :: txcod

    select type ( txcod )
    class is ( bc_bound )
       if( associated(txcod % l) ) then 
          boundary_conditions_number_codes = size(txcod % l)
       else
          boundary_conditions_number_codes = 0
       end if
     class is ( bc_nodes )
       if( associated(txcod % l) ) then 
          boundary_conditions_number_codes = size(txcod % l)
       else
          boundary_conditions_number_codes = 0
       end if
    end select
    
  end function boundary_conditions_number_codes
  
  subroutine opebcs_initialization()

    mcodc = 3_ip * mcodb

  end subroutine
 
  subroutine opnbcs(itask,nvari,ndofn,kfl_ibopo,tncod_xxx)

    !--------------------------------------------------------------------
    !
    ! NODE
    !
    !--------------------------------------------------------------------

    integer(ip),    intent(in)             :: itask,nvari,ndofn,kfl_ibopo
    integer(ip)                            :: icode,icono,idofn,ivari
    type(bc_nodes), intent(inout), pointer :: tncod_xxx(:)

    call opebcs_initialization()

    select case ( itask )

    case ( 0_ip )
       !
       ! Allocate type for NVARI variables
       !
       allocate( tncod_xxx(nvari) )
       do ivari = 1,nvari
          allocate( tncod_xxx(ivari) % l(mcodc) )
          tncod_xxx(ivari) % kfl_ibopo = kfl_ibopo
          tncod_xxx(ivari) % ndofn = ndofn  
          tncod_xxx(ivari) % ncode = 0
          do icode = 1,mcodc 
             allocate( tncod_xxx(ivari) % l(icode) % lcode(mcono) )
             allocate( tncod_xxx(ivari) % l(icode) % bvess(ndofn) )
             tncod_xxx(ivari) % l(icode) % kfl_fixno  = 0
             tncod_xxx(ivari) % l(icode) % kfl_value  = 0
             tncod_xxx(ivari) % l(icode) % kfl_funno  = 0
             tncod_xxx(ivari) % l(icode) % kfl_funtyp = 0
             tncod_xxx(ivari) % l(icode) % kfl_fixrs  = 0
             tncod_xxx(ivari) % l(icode) % tag        = ''
             tncod_xxx(ivari) % l(icode) % fname      = ''
             do icono = 1,mcono
                tncod_xxx(ivari) % l(icode) % lcode(icono) =  mcodb+1
             end do
             do idofn = 1,ndofn
                tncod_xxx(ivari) % l(icode) % bvess(idofn) =  0.0_rp
             end do
          end do
       end do

    case ( 1_ip )
       !
       ! Allocate type for NVARI variables
       !
       allocate( tncod_xxx(nvari) )

    case ( 2_ip )
       !
       ! Allocate memory for each code ICODE
       !
       allocate( tncod_xxx(nvari) % l(mcodc) )
       tncod_xxx(nvari) % kfl_ibopo = kfl_ibopo
       tncod_xxx(nvari) % ndofn = ndofn     
       tncod_xxx(nvari) % ncode = 0     
       do icode = 1,mcodc
          allocate( tncod_xxx(nvari) % l(icode) % lcode(mcono) )
          allocate( tncod_xxx(nvari) % l(icode) % bvess(ndofn) )
          tncod_xxx(nvari) % l(icode) % kfl_fixno  = 0
          tncod_xxx(nvari) % l(icode) % kfl_value  = 0
          tncod_xxx(nvari) % l(icode) % kfl_funno  = 0
          tncod_xxx(nvari) % l(icode) % kfl_funtyp = 0
          tncod_xxx(nvari) % l(icode) % kfl_fixrs  = 0
          tncod_xxx(nvari) % l(icode) % tag        = ''
          tncod_xxx(nvari) % l(icode) % fname      = ''
          do icono = 1,mcono
             tncod_xxx(nvari) % l(icode) % lcode(icono) =  mcodb+1
          end do
          do idofn = 1,ndofn
             tncod_xxx(nvari) % l(icode) % bvess(idofn) =  0.0_rp
          end do
       end do

    case ( 3_ip )
       !
       ! Deallocate whole type
       !
       do ivari = 1,size(tncod_xxx,KIND=ip)
          do icode = 1,mcodc
             deallocate( tncod_xxx(ivari) % l(icode) % lcode )
             deallocate( tncod_xxx(ivari) % l(icode) % bvess )
          end do
          deallocate( tncod_xxx(ivari) % l )
       end do
       deallocate( tncod_xxx )

    end select

  end subroutine opnbcs

  !--------------------------------------------------------------------
  !
  ! Initialization
  !
  !--------------------------------------------------------------------

  subroutine opebcs_initialization_structure_boundary(nvari,tbcod_xxx)
    implicit none
    integer(ip),             intent(in)    :: nvari
    integer(4)                             :: istat
    type(bc_bound), pointer, intent(inout) :: tbcod_xxx(:)

    if( associated(tbcod_xxx) ) call runend('STRUCTURE TBCOD ALREADY ASSOCIATED')
    allocate( tbcod_xxx(nvari) , stat = istat )
    if( istat /= 0 ) call runend('COULD NOT ALLOCATE STRUCTURE')

  end subroutine opebcs_initialization_structure_boundary

  subroutine opebcs_initialization_structure_node(nvari,tncod_xxx)
    implicit none
    integer(ip),             intent(in)    :: nvari
    integer(4)                             :: istat
    type(bc_nodes), pointer, intent(inout) :: tncod_xxx(:)

    if( associated(tncod_xxx) ) call runend('STRUCTURE TNCOD ALREADY ASSOCIATED')
    allocate( tncod_xxx(nvari) , stat = istat )
    if( istat /= 0 ) call runend('COULD NOT ALLOCATE STRUCTURE')

  end subroutine opebcs_initialization_structure_node

  subroutine opebcs_initialization_variable_boundary(ndofn,tbcod_xxx)
    implicit none
    integer(ip),             intent(in)    :: ndofn
    type(bc_bound), pointer, intent(inout) :: tbcod_xxx(:)
    integer(ip)                            :: icode,idofn,ivari


    call opebcs_initialization()

    if( .not. associated(tbcod_xxx) ) call runend('opebcs_initialization_variable_boundary: STRUCTURE TBCOD NOT ASSOCIATED')

    do ivari = 1,size(tbcod_xxx,KIND=ip)
       allocate( tbcod_xxx(ivari) % l(mcodc) )
       tbcod_xxx(ivari) % ndofn = ndofn     
       tbcod_xxx(ivari) % ncode = 0
       do icode = 1,mcodc
          allocate( tbcod_xxx(ivari) % l(icode) % bvnat(ndofn) )
          tbcod_xxx(ivari) % l(icode) % kfl_fixbo  =  0
          tbcod_xxx(ivari) % l(icode) % kfl_value  =  0
          tbcod_xxx(ivari) % l(icode) % lcode      =  mcodb+1
          tbcod_xxx(ivari) % l(icode) % kfl_funtyp =  0
          tbcod_xxx(ivari) % l(icode) % kfl_funbo  =  0
          do idofn = 1,ndofn
             tbcod_xxx(ivari) % l(icode) % bvnat(idofn) =  0.0_rp
          end do
       end do
    end do

  end subroutine opebcs_initialization_variable_boundary

  subroutine opebcs_initialization_variable_boundary_s(ndofn,tbcod_xxx)
    implicit none
    integer(ip),    intent(in)    :: ndofn
    type(bc_bound), intent(inout) :: tbcod_xxx
    integer(ip)                   :: icode,idofn


    call opebcs_initialization()

    allocate( tbcod_xxx % l(mcodc) )
    tbcod_xxx % ndofn = ndofn     
    tbcod_xxx % ncode = 0
    do icode = 1,mcodc
       allocate( tbcod_xxx % l(icode) % bvnat(ndofn) )
       tbcod_xxx % l(icode) % kfl_fixbo =  0
       tbcod_xxx % l(icode) % kfl_value =  0
       tbcod_xxx % l(icode) % lcode     =  mcodb+1
       tbcod_xxx % l(icode) % kfl_funtyp =  0
       tbcod_xxx % l(icode) % kfl_funbo =  0
       do idofn = 1,ndofn
          tbcod_xxx % l(icode) % bvnat(idofn) =  0.0_rp
       end do
    end do

  end subroutine opebcs_initialization_variable_boundary_s

  subroutine opebcs_initialization_variable_node_1(ndofn,tncod_xxx,kfl_ibopo)
    implicit none
    integer(ip),              intent(in)    :: ndofn
    type(bc_nodes), pointer,  intent(inout) :: tncod_xxx(:)
    integer(ip),    optional, intent(in)    :: kfl_ibopo
    integer(ip)                             :: kbopo
    integer(ip)                             :: icode,icono,idofn,ivari

    call opebcs_initialization()

    if( present(kfl_ibopo) ) then
       kbopo = kfl_ibopo
    else
       kbopo = 0
    end if

    if( .not. associated(tncod_xxx) ) call runend('opebcs_initialization_variable_node: STRUCTURE TNCOD NOT ASSOCIATED')
    do ivari = 1,size(tncod_xxx,KIND=ip)
       allocate( tncod_xxx(ivari) % l(mcodc) )
       tncod_xxx(ivari) % kfl_ibopo = kbopo
       tncod_xxx(ivari) % ndofn = ndofn     
       tncod_xxx(ivari) % ncode = 0     
       do icode = 1,mcodc
          allocate( tncod_xxx(ivari) % l(icode) % lcode(mcono) )
          allocate( tncod_xxx(ivari) % l(icode) % bvess(ndofn) )
          tncod_xxx(ivari) % l(icode) % kfl_fixno  = 0
          tncod_xxx(ivari) % l(icode) % kfl_value  = 0
          tncod_xxx(ivari) % l(icode) % kfl_funno  = 0
          tncod_xxx(ivari) % l(icode) % kfl_funtyp = 0
          tncod_xxx(ivari) % l(icode) % kfl_fixrs  = 0
          tncod_xxx(ivari) % l(icode) % tag        = ''
          tncod_xxx(ivari) % l(icode) % fname      = ''
          do icono = 1,mcono
             tncod_xxx(ivari) % l(icode) % lcode(icono) =  mcodb+1
          end do
          do idofn = 1,ndofn
             tncod_xxx(ivari) % l(icode) % bvess(idofn) =  0.0_rp
          end do
       end do
    end do

  end subroutine opebcs_initialization_variable_node_1

  subroutine opebcs_initialization_variable_node_s(ndofn,tncod_xxx,kfl_ibopo)
    implicit none
    integer(ip),              intent(in)    :: ndofn
    type(bc_nodes),           intent(inout) :: tncod_xxx
    integer(ip),    optional, intent(in)    :: kfl_ibopo
    integer(ip)                             :: kbopo,icode,idofn,icono



    call opebcs_initialization()

    if( present(kfl_ibopo) ) then
       kbopo = kfl_ibopo
    else
       kbopo = 0
    end if

    allocate( tncod_xxx % l(mcodc) )
    tncod_xxx % kfl_ibopo = kbopo
    tncod_xxx % ndofn = ndofn     
    tncod_xxx % ncode = 0     
    do icode = 1,mcodc
       allocate( tncod_xxx % l(icode) % lcode(mcono) )
       allocate( tncod_xxx % l(icode) % bvess(ndofn) )
       tncod_xxx % l(icode) % kfl_fixno  = 0
       tncod_xxx % l(icode) % kfl_value  = 0
       tncod_xxx % l(icode) % kfl_funno  = 0
       tncod_xxx % l(icode) % kfl_funtyp = 0
       tncod_xxx % l(icode) % kfl_fixrs  = 0
       tncod_xxx % l(icode) % tag        = ''
       tncod_xxx % l(icode) % fname      = ''
       do icono = 1,mcono
          tncod_xxx % l(icode) % lcode(icono) =  mcodb+1
       end do
       do idofn = 1,ndofn
          tncod_xxx % l(icode) % bvess(idofn) =  0.0_rp
       end do
    end do

  end subroutine opebcs_initialization_variable_node_s

  subroutine opbbcs(itask,nvari,ndofn,tbcod_xxx)

    !--------------------------------------------------------------------
    !
    ! BOUNDARY
    !
    !--------------------------------------------------------------------

    implicit none
    integer(ip),             intent(in)    :: itask,nvari,ndofn
    type(bc_bound), pointer, intent(inout) :: tbcod_xxx(:)
    integer(ip)                            :: icode,ivari,idofn



    call opebcs_initialization()

    select case ( itask )

    case ( 0_ip )
       !
       ! Allocate type for NVARI variables
       !
       allocate( tbcod_xxx(nvari) )
       do ivari = 1,nvari
          allocate( tbcod_xxx(ivari) % l(mcodc) )
          tbcod_xxx(ivari) % ndofn = ndofn  
          tbcod_xxx(ivari) % ncode = 0
          do icode = 1,mcodc
             allocate( tbcod_xxx(ivari) % l(icode) % bvnat(ndofn) )
             tbcod_xxx(ivari) % l(icode) % kfl_fixbo  = 0
             tbcod_xxx(ivari) % l(icode) % kfl_value  = 0
             tbcod_xxx(ivari) % l(icode) % kfl_funtyp = 0
             tbcod_xxx(ivari) % l(icode) % kfl_funbo  = 0
             do idofn = 1,ndofn
                tbcod_xxx(ivari) % l(icode) % bvnat(idofn) =  0.0_rp
             end do
          end do
       end do

    case ( 1_ip )
       !
       ! Allocate type for NVARI variables
       !
       allocate( tbcod_xxx(nvari) )

    case ( 2_ip )
       !
       ! Allocate memory for each code ICODE
       !
       allocate( tbcod_xxx(nvari) % l(mcodc) )
       tbcod_xxx(nvari) % ndofn = ndofn     
       tbcod_xxx(nvari) % ncode = 0
       do icode = 1,mcodc
          allocate( tbcod_xxx(nvari) % l(icode) % bvnat(ndofn) )
          tbcod_xxx(nvari) % l(icode) % kfl_fixbo  =  0
          tbcod_xxx(nvari) % l(icode) % kfl_value  =  0
          tbcod_xxx(nvari) % l(icode) % lcode      =  mcodb+1
          tbcod_xxx(nvari) % l(icode) % kfl_funtyp =  0
          tbcod_xxx(nvari) % l(icode) % kfl_funbo  =  0
          do idofn = 1,ndofn
             tbcod_xxx(nvari) % l(icode) % bvnat(idofn) =  0.0_rp
          end do
       end do

    case ( 3_ip )
       !
       ! Deallocate whole type
       !
       do ivari = 1,size(tbcod_xxx,KIND=ip)
          do icode = 1,mcodc
             deallocate( tbcod_xxx(ivari) % l(icode) % bvnat )
          end do
          deallocate( tbcod_xxx(ivari) % l )
       end do
       deallocate( tbcod_xxx )

    end select

  end subroutine opbbcs

  subroutine spnbcs(tncod_xxx,INCLUDE_CHARACTER)

    !--------------------------------------------------------------------
    !
    ! PARALL for TNCOD_XXX
    !
    !--------------------------------------------------------------------

    type(bc_nodes), pointer, intent(inout) :: tncod_xxx(:)
    logical(lg),    optional, intent(in)   :: INCLUDE_CHARACTER
    integer(ip)                            :: dummi,nvari
    integer(ip)                            :: ivari,ndofn,icode
    logical(lg)                            :: if_include_character

    if_include_character = optional_argument(.true.,INCLUDE_CHARACTER)
    
    call opebcs_initialization()

    if( IPARALL .and. associated(tncod_xxx) ) then

       nvari = size(tncod_xxx,KIND=ip)
       call exchange_init()

       do ivari = 1,nvari
          ndofn = tncod_xxx(ivari) % ndofn
          call exchange_add(tncod_xxx(ivari) % kfl_ibopo)
          call exchange_add(tncod_xxx(ivari) % ndofn    )
          call exchange_add(tncod_xxx(ivari) % ncode    )
          do icode = 1,mcodc
             call exchange_add(tncod_xxx(ivari) % l(icode) % kfl_fixno)
             call exchange_add(tncod_xxx(ivari) % l(icode) % kfl_value)
             call exchange_add(tncod_xxx(ivari) % l(icode) % kfl_funno)
             call exchange_add(tncod_xxx(ivari) % l(icode) % kfl_funtyp)
             call exchange_add(tncod_xxx(ivari) % l(icode) % kfl_fixrs)   
             call exchange_add(tncod_xxx(ivari) % l(icode) % lcode)
             call exchange_add(tncod_xxx(ivari) % l(icode) % bvess) 
             if( if_include_character ) then
                call exchange_add(tncod_xxx(ivari) % l(icode) % tag)
                call exchange_add(tncod_xxx(ivari) % l(icode) % fname)
             end if
          end do

       end do
       call exchange_end()
       !
       ! Master deallocates memory
       !
       if( IMASTER ) then
          call opnbcs(3_ip,dummi,dummi,dummi,tncod_xxx)
       end if

    end if

  end subroutine spnbcs

  subroutine spbbcs(tbcod_xxx)

    !--------------------------------------------------------------------
    !
    ! PARALL for TBCOD_XXX
    !
    !--------------------------------------------------------------------

    type(bc_bound), pointer, intent(inout) :: tbcod_xxx(:)
    integer(ip)                            :: dummi,nvari
    integer(ip)                            :: ivari,ndofn,icode

    call opebcs_initialization()

    if( kfl_icodb /= 0 .and. IPARALL .and. associated(tbcod_xxx) ) then

       nvari = size(tbcod_xxx,KIND=ip)

       call exchange_init()
       do ivari = 1,nvari
          ndofn = tbcod_xxx(ivari) % ndofn
          call exchange_add(tbcod_xxx(ivari) % ndofn)
          call exchange_add(tbcod_xxx(ivari) % ncode)
          do icode = 1,mcodc
             call exchange_add(tbcod_xxx(ivari) % l(icode) % kfl_fixbo)
             call exchange_add(tbcod_xxx(ivari) % l(icode) % kfl_value)
             call exchange_add(tbcod_xxx(ivari) % l(icode) % kfl_funtyp)
             call exchange_add(tbcod_xxx(ivari) % l(icode) % kfl_funbo)
             call exchange_add(tbcod_xxx(ivari) % l(icode) % lcode)
             call exchange_add(tbcod_xxx(ivari) % l(icode) % bvnat) 
          end do
       end do
       !
       call exchange_end()
       !
       ! Master deallocates memory
       !
       if( IMASTER ) then
          call opbbcs(3_ip,dummi,dummi,tbcod_xxx)
       end if

    end if

  end subroutine spbbcs

  subroutine cpybcs(ivari,ivaro,tncod_xxx)

    !--------------------------------------------------------------------
    !
    ! Copy TNCOD_XXX to TNCOD_YYY
    !
    !--------------------------------------------------------------------

    implicit none
    integer(ip),    intent(in)    :: ivari,ivaro
    type(bc_nodes), intent(inout) :: tncod_xxx(:)
    integer(ip)                   :: icono,ndofn,icode,idofn


    call opebcs_initialization()

    if( kfl_icodn /= 0) then

       ndofn = tncod_xxx(ivari) % ndofn
       tncod_xxx(ivaro) % kfl_ibopo = tncod_xxx(ivari) % kfl_ibopo 
       tncod_xxx(ivaro) % ndofn     = tncod_xxx(ivari) % ndofn      
       tncod_xxx(ivaro) % ncode     = tncod_xxx(ivari) % ncode      
       do icode = 1,mcodc                  
          tncod_xxx(ivaro) % l(icode) % kfl_fixno  = tncod_xxx(ivari) % l(icode) % kfl_fixno
          tncod_xxx(ivaro) % l(icode) % kfl_value  = tncod_xxx(ivari) % l(icode) % kfl_value
          tncod_xxx(ivaro) % l(icode) % kfl_funno  = tncod_xxx(ivari) % l(icode) % kfl_funno
          tncod_xxx(ivaro) % l(icode) % kfl_funtyp = tncod_xxx(ivari) % l(icode) % kfl_funtyp
          tncod_xxx(ivaro) % l(icode) % kfl_fixrs  = tncod_xxx(ivari) % l(icode) % kfl_fixrs
          tncod_xxx(ivaro) % l(icode) % tag        = tncod_xxx(ivari) % l(icode) % tag
          tncod_xxx(ivaro) % l(icode) % fname      = tncod_xxx(ivari) % l(icode) % fname
          do icono = 1,mcono
             tncod_xxx(ivaro) % l(icode) % lcode(icono) = tncod_xxx(ivari) % l(icode) % lcode(icono)
          end do
          do idofn = 1,ndofn
             tncod_xxx(ivaro) % l(icode) % bvess(idofn) = tncod_xxx(ivari) % l(icode) % bvess(idofn)
          end do
       end do

    end if

  end subroutine cpybcs

  subroutine cpybcs_boundaries(ivari,ivaro,tbcod_xxx)

    !--------------------------------------------------------------------
    !
    ! Copy TNCOD_XXX to TNCOD_YYY
    !
    !--------------------------------------------------------------------

    implicit none
    integer(ip),    intent(in)    :: ivari,ivaro
    type(bc_bound), intent(inout) :: tbcod_xxx(:)
    integer(ip)                   :: idofn,ndofn,icode


    call opebcs_initialization()

    if( kfl_icodb /= 0) then
       ndofn =  tbcod_xxx(ivari) % ndofn 
       tbcod_xxx(ivaro) % ndofn = tbcod_xxx(ivari) % ndofn
       tbcod_xxx(ivaro) % ncode = tbcod_xxx(ivari) % ncode  
       do icode = 1,mcodc
          tbcod_xxx(ivaro) % l(icode) % kfl_fixbo  = tbcod_xxx(ivari) % l(icode) % kfl_fixbo
          tbcod_xxx(ivaro) % l(icode) % kfl_value  = tbcod_xxx(ivari) % l(icode) % kfl_value
          tbcod_xxx(ivaro) % l(icode) % kfl_funtyp = tbcod_xxx(ivari) % l(icode) % kfl_funtyp
          tbcod_xxx(ivaro) % l(icode) % kfl_funbo  = tbcod_xxx(ivari) % l(icode) % kfl_funbo
          tbcod_xxx(ivaro) % l(icode) % lcode      = tbcod_xxx(ivari) % l(icode) % lcode
          do idofn = 1,ndofn
             tbcod_xxx(ivaro) % l(icode) % bvnat(idofn) = tbcod_xxx(ivari) % l(icode) % bvnat(idofn) 
          end do
       end do

    end if

  end subroutine cpybcs_boundaries

  subroutine boundary_conditions_read_node_codes(message,tncod_xxx)

    character(*),   optional,          intent(in)    :: message
    type(bc_nodes), optional, pointer, intent(inout) :: tncod_xxx(:)

    if( present(tncod_xxx) ) tncod => tncod_xxx
    if( present(message) ) then
       coutp(1) = trim(message)
    else
       coutp(1) ='UNDEFINED VARIABLE'
    end if
    call codes_read_on_nodes(100_ip)
    
  end subroutine boundary_conditions_read_node_codes

  subroutine boundary_conditions_read_boundary_codes(message,tbcod_xxx)

    character(*),   optional,          intent(in)    :: message
    type(bc_bound), optional, pointer, intent(inout) :: tbcod_xxx(:)

    if( present(tbcod_xxx) ) tbcod => tbcod_xxx
    if( present(message) ) then
       coutp(1) = trim(message)
    else
       coutp(1) ='UNDEFINED VARIABLE'
    end if
    call codes_read_on_boundaries(200_ip)
    
  end subroutine boundary_conditions_read_boundary_codes

  !-------------------------------------------------------------------
  ! 
  ! Impose node or edge codes (done only for untagged nodes)
  !
  !-------------------------------------------------------------------

  ! alf alf
  ! debug debug
  ! kfl_funtyp should also be included here?
  subroutine boundary_conditions_impose_node_codes(tncod_xxx,kfl_fixno,bvess,kfl_funno,kfl_fixrs,wherein)

    use def_master, only : lun_outpu
    use def_master, only : intost
    use def_master, only : IMPOSE_NODE_CODES
    use def_master, only : IMPOSE_BOUNDARY_CODES  
    use def_master, only : IMPOSE_EDGE_CODES     
    use def_domain, only : kfl_codno
    use def_domain, only : kfl_coded
    use def_domain, only : mcono,npoin
    use def_domain, only : lpoty
    use def_domain, only : meshe,xfiel
    use def_kermod, only : ndivi
    use mod_chktyp, only : check_type
  
    type(bc_nodes),   intent(in)                       :: tncod_xxx
    integer(ip),      intent(inout), pointer           :: kfl_fixno(:,:)
    real(rp),         intent(inout), pointer, optional :: bvess(:,:)
    integer(ip),      intent(inout), pointer, optional :: kfl_funno(:)
    integer(ip),      intent(inout), pointer, optional :: kfl_fixrs(:)
    character(len=*), intent(in),             optional :: wherein

    integer(ip)                                        :: ipoin,icode,idofn,ibopo
    integer(ip)                                        :: ndofn,ncode,ibves,kpoin,ierro,itask
    integer(ip)                                        :: jcode,icono,nmcod,kcode(mcono),ivcod
    integer(ip)                                        :: ntotn,mcono_tmp
    character(20)                                      :: messa
    integer(ip), pointer                               :: kfl_codes_tmp(:,:)
    logical(lg)                                        :: if_on_nodes

    ierro = 0_ip

    if_on_nodes = .true.
    if( present(wherein) ) then
       if( trim(wherein) == 'ON EDGES') if_on_nodes = .false.
    end if

    if( if_on_nodes ) then
       itask = IMPOSE_NODE_CODES
    else
       itask = IMPOSE_EDGE_CODES
    end if

    if( itask == IMPOSE_NODE_CODES ) then
       ntotn         =  npoin
       kfl_codes_tmp => kfl_codno
       mcono_tmp     =  mcono
    else if( itask == IMPOSE_EDGE_CODES ) then
       ntotn         =  meshe(ndivi) % nedge
       kfl_codes_tmp => kfl_coded
       mcono_tmp     =  2
    end if
    
    ndofn         = memory_size(kfl_fixno,1_ip)
    if( present(bvess) ) ibves = memory_size(bvess,1_ip)
    ibves = ntotn + 1

    do ncode = 1,tncod_xxx % ncode  
       !
       ! Do it only when the code is untagged
       !   
       if( itask == IMPOSE_NODE_CODES ) then          
          if(      tncod_xxx % l(ncode) % lcode(1) == mcodb+1 ) then
             nmcod = 0
          else if( tncod_xxx % l(ncode) % lcode(2) == mcodb+1 ) then
             nmcod = 1
          else if( tncod_xxx % l(ncode) % lcode(3) == mcodb+1 ) then
             nmcod = 2
          else 
             nmcod = 3
          end if
       else
          if(      tncod_xxx % l(ncode) % lcode(1) == mcodb+1 ) then
             nmcod = 0
          else if( tncod_xxx % l(ncode) % lcode(2) == mcodb+1 ) then
             nmcod = 1
          else 
             nmcod = 2
          end if
       end if
       !
       ! Check if value function exist
       !
       if( itask == IMPOSE_NODE_CODES ) then
          if( tncod_xxx % l(ncode) % kfl_value > 0 ) then
             ivcod = tncod_xxx % l(ncode) % kfl_value
             call check_type(xfiel,ivcod,ndofn,npoin)
          end if
       end if
       !
       ! Order codes
       !
       do jcode = 1,mcono_tmp
          kcode(jcode) = tncod_xxx % l(ncode) % lcode(jcode)
       end do
       call heapsorti1(2_ip,mcono_tmp,kcode)
       icode = tncod_xxx % l(ncode) % lcode(1)

       do ipoin = 1,ntotn
          icono = 0
          do jcode = 1,mcono_tmp
             if( kfl_codes_tmp(jcode,ipoin) == abs(kcode(jcode)) ) icono = icono + 1
          end do

          if( icono == mcono_tmp ) then

             kpoin = ipoin

             if( kpoin == 0 ) then

                messa = intost(ipoin)
                ierro = ierro + 1
                call outfor(2_ip,lun_outpu,&
                     'BOUNDARY CONDITION CANNOT BE IMPOSED ON INTERIOR NODE: '//trim(messa))
             else

                kfl_fixno(1,kpoin) = tncod_xxx % l(ncode) % kfl_fixno

                call codfix(ndofn,kfl_fixno(1,kpoin))

                if( present(bvess) ) then
                   if( ibves == ntotn ) then
                      if( tncod_xxx % l(ncode) % kfl_value == 0 ) then
                         bvess(kpoin,1) = tncod_xxx % l(ncode) % bvess(1)
                      else                          
                         ivcod = tncod_xxx % l(ncode) % kfl_value
                         bvess(kpoin,1) = xfiel(ivcod) % a(1,ipoin,1)
                      end if
                   else
                      if( tncod_xxx % l(ncode) % kfl_value == 0 ) then 
                         do idofn = 1,ndofn
                            bvess(idofn,kpoin) = tncod_xxx % l(ncode) % bvess(idofn)
                         end do
                      else
                         ivcod = tncod_xxx % l(ncode) % kfl_value
                         do idofn = 1,ndofn
                            bvess(idofn,kpoin) = xfiel(ivcod) % a(idofn,ipoin,1)
                         end do
                      end if
                   end if
                end if

                if( present(kfl_funno) ) then
                   kfl_funno(kpoin) = tncod_xxx % l(ncode) % kfl_funno                       
                   if( present(kfl_fixrs) ) then 
                      ibopo = lpoty(kpoin)
                      if( ibopo /= 0 ) then                         
                         kfl_fixrs(ibopo) = tncod_xxx % l(ncode) % kfl_fixrs
                         !if( kfl_fixrs(ibopo) == -2 .and. nskew > 0 ) call geofix(kpoin,ibopo)
                      end if
                   end if
                else if( present(kfl_fixrs) ) then
                   ibopo = lpoty(kpoin)
                   if( ibopo /= 0 ) then
                      kfl_fixrs(ibopo) = tncod_xxx % l(ncode) % kfl_fixrs
                      !if(kfl_fixrs(ibopo)==-2.and.nskew>0) call geofix(kpoin,ibopo)
                   end if
                end if
             end if

          end if

       end do

       !end if

    end do

  end subroutine boundary_conditions_impose_node_codes

  !-------------------------------------------------------------------
  ! 
  ! Impose node or edge codes (done only for untagged nodes)
  !
  !-------------------------------------------------------------------
  
  subroutine boundary_conditions_impose_edge_codes_1(tncod_xxx,kfl_fixno,bvess,kfl_funno,kfl_fixrs)

    type(bc_nodes), intent(in)                       :: tncod_xxx
    integer(ip),    intent(inout), pointer           :: kfl_fixno(:,:)
    real(rp),       intent(inout), pointer           :: bvess(:,:)
    integer(ip),    intent(inout), pointer, optional :: kfl_funno(:)
    integer(ip),    intent(inout), pointer, optional :: kfl_fixrs(:)

    call boundary_conditions_impose_node_codes(tncod_xxx,kfl_fixno,bvess,kfl_funno,kfl_fixrs,'ON EDGES')

  end subroutine boundary_conditions_impose_edge_codes_1
  
  subroutine boundary_conditions_impose_edge_codes_2(tncod_xxx,kfl_fixno,bvess,kfl_funno,kfl_fixrs)

    type(bc_nodes), intent(in)                       :: tncod_xxx
    integer(ip),    intent(inout), pointer           :: kfl_fixno(:,:)
    real(rp),       intent(inout), pointer           :: bvess(:,:,:)
    integer(ip),    intent(inout), pointer, optional :: kfl_funno(:)
    integer(ip),    intent(inout), pointer, optional :: kfl_fixrs(:)
    real(rp),                      pointer           :: bvess_loc(:,:)

    if( associated(bvess) ) then
       bvess_loc => bvess(:,:,1)
    else
       nullify(bvess_loc)   
    end if
    call boundary_conditions_impose_node_codes(tncod_xxx,kfl_fixno,bvess_loc,kfl_funno,kfl_fixrs,'ON EDGES')

  end subroutine boundary_conditions_impose_edge_codes_2
  
end module mod_opebcs

