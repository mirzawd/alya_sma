!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    mod_boundary_conditions.f90
!> @author  houzeaux
!> @date    2020-04-13
!> @brief   Boundary conditions
!> @details Module to read and impose boundary conditions
!-----------------------------------------------------------------------

module mod_boundary_conditions

  use def_kintyp
  use def_master
  use def_domain  
  use def_elmtyp,         only : BOEXT
  use def_kermod,         only : kfl_edge_elements
  use def_kermod,         only : ndivi
  use mod_elmgeo,         only : element_type 
  use mod_outfor,         only : outfor
  use mod_communications, only : PAR_GHOST_BOUNDARY_EXCHANGE
  use mod_communications, only : PAR_INTERFACE_EDGE_EXCHANGE
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_SUM
  use mod_communications, only : PAR_GATHER
  use mod_communications, only : PAR_GATHERV
  use mod_memory,         only : memory_size
  use mod_outfor,         only : outfor
  use mod_messages,       only : messages_live
  use mod_domain,         only : domain_memory_allocate  
  use mod_parall,         only : PAR_WORLD_SIZE
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_memory,         only : memory_resize
  use mod_maths,          only : maths_heap_sort
  
  use mod_opebcs,         only : boundary_conditions_initialization_structure => opebcs_initialization_structure 
  use mod_opebcs,         only : boundary_conditions_initialization_variable  => opebcs_initialization_variable
  use mod_opebcs,         only : boundary_conditions_read_node_codes
  use mod_opebcs,         only : boundary_conditions_read_boundary_codes   
  use mod_opebcs,         only : boundary_conditions_impose_node_codes
  use mod_opebcs,         only : boundary_conditions_impose_edge_codes

  implicit none

  private

  interface boundary_conditions_read_edge_codes 
     module procedure boundary_conditions_read_node_codes 
  end interface boundary_conditions_read_edge_codes

  public :: boundary_conditions_read_edge_codes  
  public :: boundary_conditions_read_node_codes  
  public :: boundary_conditions_read_boundary_codes   
  public :: boundary_conditions_impose_node_codes         ! Impose conditions on nodes
  public :: boundary_conditions_impose_edge_codes         ! Impose conditions on nodes
  public :: boundary_conditions_initialization_structure  ! Allocate bc structure
  public :: boundary_conditions_initialization_variable   ! Initialize bc structure
  
  public :: boundary_conditions_extrapolate               ! Extrapolate boundary conditions
  public :: boundary_conditions_boundaries_to_nodes       ! Extrapolate from boundaries to nodes
  public :: boundary_conditions_check_codes

  public :: boundary_conditions_number_boundary_codes     ! Number of boundary codes
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-13
  !> @brief   Extrapolate boundary conditions
  !> @details Extrapolate boundary conditions from boundaries
  !>          to node and edges
  !> 
  !-----------------------------------------------------------------------

  subroutine boundary_conditions_extrapolate()

    integer(ip) :: ipoin,iedge

    if( kfl_extra == 1 ) then
       !
       ! On nodes
       !
       call messages_live('EXTRAPOLATE BOUNDARY CODES TO NODE CODES')
       if( kfl_icodn == 0 ) then
          kfl_icodn = 1
          call domain_memory_allocate('KFL_CODNO')
          do ipoin = 1,npoin
             kfl_codno(1,ipoin) = mcodb+1
          end do
       end if
       call boundary_conditions_boundaries_to_nodes()
       !
       ! On edges
       !
       if( kfl_edge_elements == 1 ) then
          if( kfl_icode == 0 ) then
             kfl_icode = 1
             call domain_memory_allocate('KFL_CODED')
             do iedge = 1,meshe(ndivi) % nedge
                kfl_coded(:,iedge) = mcodb+1
             end do
          end if
          call boundary_conditions_boundaries_to_edges(meshe(ndivi))          
       end if

    end if
   
  end subroutine boundary_conditions_extrapolate

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-13
  !> @brief   Extrapolate
  !> @details Extrapolate from boundaries to nodes
  !> 
  !-----------------------------------------------------------------------

  subroutine boundary_conditions_boundaries_to_nodes()

    integer(ip)   :: ipoin,iboun,inodb,kodbo
    integer(ip)   :: kcode,icono,mcod1,isize
    character(20) :: mess1,mess2
    !
    ! Check if we have to extrapolate something
    !
    isize = memory_size(kfl_codbo)
    call PAR_MAX(isize)
    if( isize == 0 ) return

    mess1 = intost(mcono)
    !
    ! Parall: Exchange KFL_CODBO on fringe nodes
    !  
    if( ISLAVE ) then
       call memgen(1_ip,nboun_2,0_ip)
       do iboun = 1,nboun
          gisca(iboun) = kfl_codbo(iboun)
       end do
       call PAR_GHOST_BOUNDARY_EXCHANGE(gisca,'SUBSTITUTE','IN MY CODE')
    else 
       gisca => kfl_codbo
    end if
    !
    ! Extrapolate to a dummi array GIVEC
    !
    call memgen(1_ip,mcono,npoin_2)
    mcod1 = mcodb+1
    do ipoin = 1,npoin_2
       do icono = 1,mcono
          givec(icono,ipoin) = mcod1
       end do
    end do

    do iboun = 1,nboun_2
       if( lboch(iboun) /= BOEXT ) then
          kodbo = gisca(iboun)
          if( abs(kodbo) <= mcodb ) then
             do inodb = 1,nnode(ltypb(iboun))
                ipoin = lnodb(inodb,iboun)
                !
                ! Look for a place
                !
                kcode = 1
                mcono_loop2: do while( givec(kcode,ipoin) /= mcod1 )
                   if( kodbo == givec(kcode,ipoin) ) kcode = mcono + 1           
                   kcode = kcode + 1
                   if( kcode > mcono ) exit mcono_loop2
                end do mcono_loop2

                if( kcode <= mcono ) then
                   if( givec(kcode,ipoin) /= mcod1 ) then
                      mess2=intost(ipoin)
                      call outfor(2_ip,lun_outpu,'NODE '//trim(mess1)&
                           //' HAVE MORE THAN '//trim(mess2))
                   else
                      givec(kcode,ipoin) = kodbo
                   end if
                end if

             end do

          end if
       end if

    end do
    !
    ! Deallocate memory
    !
    if( ISLAVE ) then
       call memgen(3_ip,nboun_2,0_ip)
    else
       nullify(gisca)
    end if
    !
    ! Sort codes
    !
    do ipoin = 1,npoin
       kcode = 1
       kcode_loop2: do while( givec(kcode,ipoin) /= mcod1 )        
          kcode = kcode + 1
          if( kcode > mcono ) exit kcode_loop2
       end do kcode_loop2
       kcode = kcode - 1
       call heapsorti1(2_ip,kcode,givec(1,ipoin))
    end do
    !
    ! Write over KFL_CODNO
    ! 
    do ipoin = 1,npoin
       if( kfl_codno(1,ipoin) == mcod1 ) then
          do icono = 1,mcono
             kfl_codno(icono,ipoin) = givec(icono,ipoin)
          end do
       end if
    end do
    call memgen(3_ip,mcono,npoin_2)

  end subroutine boundary_conditions_boundaries_to_nodes

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-13
  !> @brief   Extrapolate
  !> @details Extrapolate from boundaries to nodes
  !> 
  !-----------------------------------------------------------------------

  subroutine boundary_conditions_boundaries_to_edges(meshe_in)

    type(mesh_type), intent(in) :: meshe_in
    integer(ip)                 :: iboun
    integer(ip)                 :: mcod1,isize
    integer(ip)                 :: iedge
    integer(ip)                 :: iedgg,icodb
    integer(ip)                 :: icodb_min,icodb_max
    integer(ip), pointer        :: codeg(:,:)

    if( kfl_edge_elements /= 0 ) then
       !
       ! Check if we have to extrapolate something
       !
       nullify(codeg)
       isize = memory_size(kfl_codbo)
       call PAR_MAX(isize)
       if( isize == 0 ) return
       mcod1 = mcodb+1
       !
       ! Parall: Exchange KFL_CODBO on fringe nodes
       !
       call memgen(1_ip,nboun_2,0_ip)
       do iboun = 1,nboun
          gisca(iboun) = kfl_codbo(iboun)
       end do
       if( IPARALL ) call PAR_GHOST_BOUNDARY_EXCHANGE(gisca,'SUBSTITUTE','IN MY CODE')
       !
       ! Inherit code of boundaries
       !
       do iboun = 1,meshe_in % nboun
          icodb = gisca(iboun)
          if( icodb /= mcod1 ) then
             do iedge = 1,meshe_in % lnneb(iboun)
                iedgg = meshe_in % ledgb(iedge,iboun)
                if(      kfl_coded(1,iedgg) == mcod1 ) then
                   kfl_coded(1,iedgg) = icodb
                else if( kfl_coded(2,iedgg) == mcod1 .and. icodb /= kfl_coded(1,iedgg) ) then
                   kfl_coded(2,iedgg) = icodb
                end if
             end do
          end if
       end do
       !
       ! Merge with my neighbors... the MERGE option requires 0 values as free positions...
       !
       do iedgg = 1,meshe_in % nedge
          if( kfl_coded(1,iedgg) == mcod1 ) kfl_coded(1,iedgg) = 0
          if( kfl_coded(2,iedgg) == mcod1 ) kfl_coded(2,iedgg) = 0
       end do
       call PAR_INTERFACE_EDGE_EXCHANGE(kfl_coded,'MERGE')
       do iedgg = 1,meshe_in % nedge
          if( kfl_coded(1,iedgg) == 0     ) kfl_coded(1,iedgg) = mcod1
          if( kfl_coded(2,iedgg) == 0     ) kfl_coded(2,iedgg) = mcod1
       end do
       !
       ! Deallocate memory
       !
       call memgen(3_ip,nboun_2,0_ip)
       !
       ! Sort codes
       !
       do iedgg = 1,meshe_in % nedge
          if( kfl_coded(1,iedgg) >= 0 .and. kfl_coded(2,iedgg) >= 0 ) then
             icodb_min          = min( kfl_coded(1,iedgg),kfl_coded(2,iedgg) )
             icodb_max          = max( kfl_coded(1,iedgg),kfl_coded(2,iedgg) )
             kfl_coded(1,iedgg) = icodb_min
             kfl_coded(2,iedgg) = icodb_max
          end if
       end do

    end if

  end subroutine boundary_conditions_boundaries_to_edges

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-14
  !> @brief   Check codes
  !> @details Check node and edge codes
  !> 
  !-----------------------------------------------------------------------

  subroutine boundary_conditions_check_codes()
    
    if( kfl_icodn /= 0 ) then
       call boundary_conditions_check_code(kfl_codno,'NODES')
    end if
    if( kfl_edge_elements == 1 ) then
       call boundary_conditions_check_code(kfl_coded,'EDGES')
    end if

  end subroutine boundary_conditions_check_codes

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-14
  !> @brief   Check a specific code
  !> @details Check a specific code, on node or edges
  !> 
  !-----------------------------------------------------------------------
  
  subroutine boundary_conditions_check_code(kfl_codes,on_what)

    use mod_std

    integer(ip), pointer, intent(in) :: kfl_codes(:,:)
    character(len=*),     intent(in) :: on_what
    integer(ip)                      :: ienti,icono,icode,kcomb,icomb,mcod1
    integer(ip)                      :: icomb_par,mcodm,nenti
    logical(lg)                      :: lcodt(-mcodb-1:mcodb+1)
    integer(ip)                      :: ncomb_par
    integer(ip)                      :: ncomb
    integer(4)                       :: ncomb4
    integer(4),  pointer             :: ncomb_gat4(:) 
    integer(ip), pointer             :: lcomb(:,:) 
    integer(ip), pointer             :: lcomb_gat(:,:) 

    mcodm = memory_size(kfl_codes,1_ip)
    nenti = memory_size(kfl_codes,2_ip)
    call PAR_MAX(mcodm)
    
    if( mcodm > 0 ) then

       nullify(lcomb)
       nullify(lcomb_gat)
       nullify(ncomb_gat4)

       !----------------------------------------------------------------------
       !
       ! Compute possible combinations LCOMB
       !
       !----------------------------------------------------------------------

       if( INOTMASTER ) then
          !
          ! NCOMB: Number of different imposed codes
          !
          mcod1 = mcodb + 1
          lcodt = .false.
          do ienti = 1,nenti
             do icono = 1,mcodm
                icode = kfl_codes(icono,ienti)
                lcodt(icode) = .true.           
             end do
          end do

          ncomb = count(lcodt(-mcodb:mcodb),KIND=ip)
          !
          ! NCOMB: Number of possible combinations
          !
          if( mcodm == 3 ) then
             ncomb =    ncomb &
                  &  +  ncomb * ( ncomb - 1 ) / 2 &
                  &  +  ncomb * ( ncomb - 1 ) * ( ncomb - 2 ) / 3
          else if( mcodm == 2 ) then
             ncomb =    ncomb &
                  &  +  ncomb * ( ncomb - 1 ) / 2
          else
             call runend('MOD_BOUNDARY_CONDITIONS: NOT POSSIBLE')
          end if

          call memory_alloca(memor_dom,'LCOMB','chkcod',lcomb,mcodm,ncomb)
          !
          ! Fill in combination table
          !
          if( ncomb > 0 ) then

             lcomb = mcod1
             kcomb = 0

             do ienti = 1,nenti

                if( kfl_codes(1,ienti) /= mcod1 ) then
                   icomb = 1
                   icode = lcomb(1,1) 
                   do while( icode /= mcod1 )
                      if(    kfl_codes(    1,ienti) == lcomb(     1,icomb) .and. &
                           & kfl_codes(    2,ienti) == lcomb(     2,icomb) .and. &
                           & kfl_codes(mcodm,ienti) == lcomb(mcodm,icomb) ) then
                         icomb = -1
                         icode = mcod1
                      else
                         icomb = icomb + 1
                         icode = lcomb(1,icomb)
                      end if
                   end do
                   !
                   ! One new combination has been found
                   !
                   if( icomb /= -1 ) then
                      kcomb                = kcomb + 1
                      lcomb(1:mcodm,kcomb) = kfl_codes(1:mcodm,ienti)
                   end if

                end if

             end do
             !
             ! Sort combinations
             !
             ncomb = kcomb
             do icomb = 1,ncomb
                if( lcomb(    2,icomb) == mcod1 ) lcomb(2,icomb)     = -mcod1 
                if( lcomb(mcodm,icomb) == mcod1 ) lcomb(mcodm,icomb) = -mcod1 
             end do
             call memory_resize(memor_dom,'LCOMB','chkcod',lcomb,mcodm,ncomb)
             if( ncomb > 0 ) call maths_heap_sort(2_ip,lcomb,NROWS=ncomb)
          end if

       end if

       !----------------------------------------------------------------------
       !
       ! Gather number of combinations NCOMB_GAT4 and list of combinations LCOMB_GAT
       !
       !----------------------------------------------------------------------

       if( IPARALL ) then
          if( IMASTER ) then
             ncomb = 0
             call memory_alloca(memor_dom,'NCOMB_GAT4','chkcod',ncomb_gat4,int(PAR_WORLD_SIZE,4),lboun=0_4)
          end if
          ncomb4 = int(ncomb,4)
          call PAR_GATHER(ncomb4,ncomb_gat4)
          if( IMASTER ) then
             ncomb_par  = int(sum(ncomb_gat4),ip)
             ncomb_gat4 = ncomb_gat4*int(mcodm,4)
             call memory_alloca(memor_dom,'LCOMB_GAT','chkcod',lcomb_gat,mcodm,ncomb_par)
          end if
          call PAR_GATHERV(lcomb,lcomb_gat,ncomb_gat4)
          if( IMASTER ) then
             call memory_alloca(memor_dom,'LCOMB','chkcod',lcomb,mcodm,ncomb_par)
             ncomb = 0    
             do icomb_par = 1,ncomb_par
                icomb = 1
                icode = lcomb(1,icomb)
                do while( icode > 0 ) 
                   if(          lcomb(    1,icomb) == lcomb_gat(    1,icomb_par) &
                        & .and. lcomb(    2,icomb) == lcomb_gat(    2,icomb_par) &
                        & .and. lcomb(mcodm,icomb) == lcomb_gat(mcodm,icomb_par) ) then
                      icode = -1
                   else
                      icomb = icomb + 1
                      icode = lcomb(1,icomb)
                   end if
                end do
                if( icode /= -1 ) then
                   ncomb                = ncomb + 1
                   lcomb(1:mcodm,ncomb) = lcomb_gat(1:mcodm,icomb_par)
                end if
             end do
          end if
       end if

       !----------------------------------------------------------------------
       !
       ! Output bombinations LCOMB
       !
       !----------------------------------------------------------------------
       if( INOTSLAVE .and. ncomb /= 0 ) then
          call maths_heap_sort(2_ip,lcomb,NROWS=ncomb)
          coutp(1) = trim(on_what)
          call outfor(42_ip,lun_outpu,' ')
          do icomb = 1,ncomb
             ioutp(1) = lcomb(    1,icomb)
             ioutp(2) = lcomb(    2,icomb)
             ioutp(3) = lcomb(mcodm,icomb)
             ioutp(4) = mcodm
             call outfor(43_ip,lun_outpu,' ')
          end do
       end if
       call memory_deallo(memor_dom,'NCOMB_GAT4','chkcod',ncomb_gat4)
       call memory_deallo(memor_dom,'LCOMB_GAT' ,'chkcod',lcomb_gat)
       call memory_deallo(memor_dom,'LCOMB'     ,'chkcod',lcomb)
    end if

  end subroutine boundary_conditions_check_code

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-11-09
  !> @brief   Number fo boundary codes
  !> @details Number fo boundary codes
  !> 
  !-----------------------------------------------------------------------
  
  integer(ip) function boundary_conditions_number_boundary_codes()

    integer(ip), pointer :: list_codes(:)
    integer(ip)          :: iboun,icode

    if( mcodb > 0 ) then

       nullify(list_codes)

       call memory_alloca(memor_dom,'LIST_CODES','mod_boundary_conditions',list_codes,mcodb+2_ip,LBOUN=0_ip)

       if( associated(kfl_codbo) ) then
          do iboun = 1,nboun
             icode             = kfl_codbo(iboun)
             list_codes(icode) = 1
          end do
       end if
       call PAR_SUM(list_codes)

       boundary_conditions_number_boundary_codes = count(list_codes/=0)

       call memory_deallo(memor_dom,'LIST_CODES','mod_boundary_conditions',list_codes)

    end if

  end function boundary_conditions_number_boundary_codes
  
end module mod_boundary_conditions
!> @}
