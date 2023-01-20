!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    mod_domain_memory.f90
!> @author  houzeaux
!> @date    2019-09-22
!> @brief   Manage domain 
!> @details Memory, functions, etc.
!-----------------------------------------------------------------------

module mod_domain

  use def_kintyp,         only : ip
  use def_master,         only : lninv_loc
  use def_master,         only : leinv_loc
  use def_master,         only : lbinv_loc
  use def_master,         only : lginv_loc
  use def_master,         only : NPOIN_TYPE
  use def_master,         only : NBOUN_TYPE
  use def_master,         only : NELEM_TYPE
  use def_master,         only : ISEQUEN
  use def_master,         only : lumma
  use def_master,         only : I_AM_IN_SUBD
  use def_master,         only : I_AM_IN_ZONE
  use def_master,         only : intost
  use def_master,         only : lnwit
  use def_elmtyp,         only : element_end
  use def_kermod,         only : kfl_noslw_ker
  use def_kermod,         only : cowit
  use def_kermod,         only : cowit_origi
  use def_kermod,         only : dewit
  use def_kermod,         only : lewit
  use def_kermod,         only : shwit
  use def_kermod,         only : gewit
  use def_kermod,         only : mwitn
  use def_kermod,         only : nwitn
  use def_kermod,         only : nwitg
  use def_kermod,         only : nwith
  use def_kermod,         only : witness_mesh
  use def_kermod,         only : ndivi
  use mod_mpio_config,    only : mpio_config
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_graphs,         only : graphs_elepoi_deallocate
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_SUM
  use def_domain

  implicit none

  private

  public :: domain_memory_allocate
  public :: domain_memory_deallocate
  public :: domain_memory_reallocate
  public :: domain_total_node_number
  public :: domain_number_gauss_points
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-09-22
  !> @brief   Reallocate array
  !> @details Reallocate a domain array
  !> 
  !-----------------------------------------------------------------------

  subroutine domain_memory_reallocate(wvari,NUMBER1,NUMBER2)
    character(len=*), intent(in)           :: wvari
    integer(ip),      intent(in), optional :: NUMBER1
    integer(ip),      intent(in), optional :: NUMBER2
    call domain_memory_deallocate(wvari,NUMBER1,NUMBER2)
    call domain_memory_allocate  (wvari,NUMBER1,NUMBER2)
  end subroutine domain_memory_reallocate
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-09-22
  !> @brief   Allocate array
  !> @details Allocate a domain array
  !> 
  !-----------------------------------------------------------------------

  subroutine domain_memory_allocate(wvari,NUMBER1,NUMBER2)
    character(len=*), intent(in)           :: wvari
    integer(ip),      intent(in), optional :: NUMBER1
    integer(ip),      intent(in), optional :: NUMBER2
    call domain_memory_allocate_deallocate(0_ip,wvari,NUMBER1,NUMBER2)
  end subroutine domain_memory_allocate
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-09-22
  !> @brief   Deallocate array
  !> @details Deallocate a domain array
  !> 
  !-----------------------------------------------------------------------

  subroutine domain_memory_deallocate(wvari,NUMBER1,NUMBER2)
    character(len=*), intent(in)           :: wvari
    integer(ip),      intent(in), optional :: NUMBER1
    integer(ip),      intent(in), optional :: NUMBER2
    call domain_memory_allocate_deallocate(1_ip,wvari,NUMBER1,NUMBER2)
  end subroutine domain_memory_deallocate
  
  subroutine domain_memory_allocate_deallocate(itask,wvari,NUMBER1,NUMBER2)

    integer(ip),      intent(in)           :: itask
    character(len=*), intent(in)           :: wvari
    integer(ip),      intent(in), optional :: NUMBER1
    integer(ip),      intent(in), optional :: NUMBER2
    integer(ip)                            :: nbou1,ipoin,ielem,nsteps
    integer(ip)                            :: ifiel,ifiel1,ifiel2,ii

    nbou1 = nboun !max(1_ip,nboun)

    select case(trim(wvari))

    case ( 'GEOMETRY' )
       !
       ! Arrays read in reageo 
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'LTYPE'    ,'mod_domain_memory' , ltype     , nelem   )
          call memory_alloca(memor_dom,'LELCH'    ,'mod_domain_memory' , lelch     , nelem   )
          call memory_alloca(memor_dom,'LNODS'    ,'mod_domain_memory' , lnods     , mnode   , nelem )
          call memory_alloca(memor_dom,'LESUB'    ,'mod_domain_memory' , lesub     , nelem   )
          call memory_alloca(memor_dom,'LMATE'    ,'mod_domain_memory' , lmate     , nelem )
          !call memory_alloca(memor_dom,'LEINV_LOC','mod_domain_memory' , leinv_loc , nelem, 'IDENTITY') 

          call memory_alloca(memor_dom,'COORD'    ,'mod_domain_memory' , coord     , ndime   , npoin )
          call memory_alloca(memor_dom,'LNOCH'    ,'mod_domain_memory' , lnoch     , npoin   )
          call memory_alloca(memor_dom,'LMAST'    ,'mod_domain_memory' , lmast     , npoin   )
          !call memory_alloca(memor_dom,'LNINV_LOC','mod_domain_memory' , lninv_loc , npoin , 'IDENTITY')

          call memory_alloca(memor_dom,'LNODB'    ,'mod_domain_memory' , lnodb     , mnodb   , nbou1 )
          call memory_alloca(memor_dom,'LTYPB'    ,'mod_domain_memory' , ltypb     , nbou1   )
          call memory_alloca(memor_dom,'LBOCH'    ,'mod_domain_memory' , lboch     , nbou1   )
          call memory_alloca(memor_dom,'LELBO'    ,'mod_domain_memory' , lelbo     , nbou1   )
          !call memory_alloca(memor_dom,'LBINV_LOC','mod_domain_memory' , lbinv_loc , nbou1 , 'IDENTITY')

          if( ISEQUEN ) then
             call memory_alloca(memor_dom,'LNINV_LOC','mod_domain_memory',lninv_loc,npoin,'IDENTITY')
             call memory_alloca(memor_dom,'LEINV_LOC','mod_domain_memory',leinv_loc,nelem,'IDENTITY')
             call memory_alloca(memor_dom,'LBINV_LOC','mod_domain_memory',lbinv_loc,nboun,'IDENTITY')
          else
             !
             ! In parallel, these arrays are computed in the associated modules
             !
          end if
          do ielem = 1,nelem
             lesub(ielem) = 1
             lmate(ielem) = 1
          end do
       else
          call memory_deallo(memor_dom,'LTYPE'    ,'mod_domain_memory' , ltype )
          call memory_deallo(memor_dom,'LELCH'    ,'mod_domain_memory' , lelch )
          call memory_deallo(memor_dom,'LNODS'    ,'mod_domain_memory' , lnods )
          call memory_deallo(memor_dom,'LESUB'    ,'mod_domain_memory' , lesub )
          call memory_deallo(memor_dom,'LMATE'    ,'mod_domain_memory' , lmate )
          !call memory_deallo(memor_dom,'LEINV_LOC','mod_domain_memory' , leinv_loc )

          call memory_deallo(memor_dom,'COORD'    ,'mod_domain_memory' , coord )
          call memory_deallo(memor_dom,'LNOCH'    ,'mod_domain_memory' , lnoch )
          call memory_deallo(memor_dom,'LMAST'    ,'mod_domain_memory' , lmast )
          !call memory_deallo(memor_dom,'LNINV_LOC','mod_domain_memory' , lninv_loc )

          call memory_deallo(memor_dom,'LNODB'    ,'mod_domain_memory' , lnodb )
          call memory_deallo(memor_dom,'LBOEL'    ,'mod_domain_memory' , lboel )
          call memory_deallo(memor_dom,'LTYPB'    ,'mod_domain_memory' , ltypb )
          call memory_deallo(memor_dom,'LBOCH'    ,'mod_domain_memory' , lboch )
          call memory_deallo(memor_dom,'LELBO'    ,'mod_domain_memory' , lelbo )
          !call memory_deallo(memor_dom,'LBINV_LOC','mod_domain_memory' , lbinv_loc )

          if( ISEQUEN ) then
             call memory_deallo(memor_dom,'LNINV_LOC','mod_domain_memory',lninv_loc)
             call memory_deallo(memor_dom,'LEINV_LOC','mod_domain_memory',leinv_loc)
             call memory_deallo(memor_dom,'LBINV_LOC','mod_domain_memory',lbinv_loc)
          else
             !
             !
             !
          end if

       end if

    case ( 'ALL MESH' )
       !
       !
       if( itask == 0 ) then
          call runend('MOD_DOMAIN_MEMORY: DO NOT KNOW WHAT TO DO')
       else
          !
          ! Deallocate memory: Mesh arrays
          !
          call memory_deallo(memor_dom,'LTYPE'    ,'mod_domain_memory' , ltype )
          call memory_deallo(memor_dom,'LNNOD'    ,'mod_domain_memory' , lnnod )
          call memory_deallo(memor_dom,'LELCH'    ,'mod_domain_memory' , lelch )
          call memory_deallo(memor_dom,'LNODS'    ,'mod_domain_memory' , lnods )
          call memory_deallo(memor_dom,'LESUB'    ,'mod_domain_memory' , lesub )
          call memory_deallo(memor_dom,'LMATE'    ,'mod_domain_memory' , lmate )
          call memory_deallo(memor_dom,'LEINV_LOC','mod_domain_memory' , leinv_loc )

          call memory_deallo(memor_dom,'COORD'    ,'mod_domain_memory' , coord )
          call memory_deallo(memor_dom,'LNOCH'    ,'mod_domain_memory' , lnoch )
          call memory_deallo(memor_dom,'LMAST'    ,'mod_domain_memory' , lmast )
          call memory_deallo(memor_dom,'LNINV_LOC','mod_domain_memory' , lninv_loc )

          call memory_deallo(memor_dom,'LNODB'    ,'mod_domain_memory' , lnodb )
          call memory_deallo(memor_dom,'LTYPB'    ,'mod_domain_memory' , ltypb )
          call memory_deallo(memor_dom,'LBOCH'    ,'mod_domain_memory' , lboch )
          call memory_deallo(memor_dom,'LELBO'    ,'mod_domain_memory' , lelbo )
          call memory_deallo(memor_dom,'LBINV_LOC','mod_domain_memory' , lbinv_loc )
          !
          ! Dealloctae R_DOM and C_DOM to minimum size (used in deflated CG)
          !
          call memory_deallo(memor_dom,'R_DOM','mod_domain_memory' , r_dom )
          call memory_deallo(memor_dom,'C_DOM','mod_domain_memory' , c_dom )
          call memory_alloca(memor_dom,'R_DOM','mod_domain_memory' , r_dom , 1_ip )
          call memory_alloca(memor_dom,'C_DOM','mod_domain_memory' , c_dom , 1_ip )
          !
          ! Deallocate LELPO and PELPO (used in parall)
          !
          call graphs_elepoi_deallocate(pelpo,lelpo,memor=memor_dom)

       end if

    case( 'LNOCH' )
       !
       ! LNOCH
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'LNOCH'    ,'mod_domain_memory' , lnoch     , npoin)
       else
          call memory_deallo(memor_dom,'LNOCH'    ,'mod_domain_memory' , lnoch     )
       end if

    case( 'LBOEL' )
       !
       ! LBOEL
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'LBOEL'    ,'mod_domain_memory' , lboel     , mnodb   , nbou1 )
       else
          call memory_deallo(memor_dom,'LBOEL'    ,'mod_domain_memory' , lboel     )
       end if

    case( 'LNINV_LOC' )
       !
       ! LNINV_LOC
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'LNINV_LOC'    ,'mod_domain_memory' , lninv_loc     , npoin )
       else
          call memory_deallo(memor_dom,'LNINV_LOC'    ,'mod_domain_memory' , lninv_loc     )
       end if

    case( 'LEINV_LOC' )
       !
       ! LEINV_LOC
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'LEINV_LOC'    ,'mod_domain_memory' , leinv_loc     , nelem )
       else
          call memory_deallo(memor_dom,'LEINV_LOC'    ,'mod_domain_memory' , leinv_loc     )
       end if

    case( 'LBINV_LOC' )
       !
       ! LBINV_LOC
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'LBINV_LOC'    ,'mod_domain_memory' , lbinv_loc     , nboun )
       else
          call memory_deallo(memor_dom,'LBINV_LOC'    ,'mod_domain_memory' , lbinv_loc     )
       end if

    case( 'LGINV_LOC' )
       !
       ! LGINV_LOC
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'LGINV_LOC'    ,'mod_domain_memory' , lginv_loc     , nedge )
       else
          call memory_deallo(memor_dom,'LGINV_LOC'    ,'mod_domain_memory' , lginv_loc     )
       end if

    case( 'SKCOS' )
       !
       ! Allocate SKCOS
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'SKCOS','mod_domain_memory' , skcos , ndime , ndime , nbopo )
       else
          call memory_deallo(memor_dom,'SKCOS','mod_domain_memory' , skcos )
       end if

    case( 'EXNOR' )
       !
       ! Allocate EXNOR
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'EXNOR','mod_domain_memory' , exnor , ndime , ndime , nbopo )
       else
          call memory_deallo(memor_dom,'EXNOR','mod_domain_memory' , exnor )
       end if

    case( 'LPOTY' )
       !
       ! Allocate LPOTY
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'LPOTY','mod_domain_memory' , lpoty , npoin )
       else
          call memory_deallo(memor_dom,'LPOTY','mod_domain_memory' , lpoty )
       end if

    case( 'LNLEV' )
       !
       ! Allocate LNLEV
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'LNLEV','mod_domain_memory' , lnlev , npoin )
       else
          call memory_deallo(memor_dom,'LNLEV','mod_domain_memory' , lnlev )
       end if

    case( 'BOUNDARY GRAPH' , 'R_BOU AND C_BOU' )
       !
       ! Allocate  boundary graph: R_BOU and C_BOU
       !
       if( itask == 0 ) then
          if( kfl_crbou == 0 ) then
             call memory_alloca(memor_dom,'R_BOU','mod_domain_memory' , r_bou , npoin+1 )
             call memory_alloca(memor_dom,'C_BOU','mod_domain_memory' , c_bou , nzbou   )
             kfl_crbou = 1
          end if
       else
          if( kfl_crbou == 1 ) then
             call memory_deallo(memor_dom,'R_BOU','mod_domain_memory' , r_bou )
             call memory_deallo(memor_dom,'C_BOU','mod_domain_memory' , c_bou )         
             kfl_crbou = 0
          end if
       end if

    case( 'YWALB, YWALP AND YWALE' , 'VARIABLE WALL DISTANCE' )
       !
       ! Arrays for variable wall distance on boundaries and boundary nodes
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'YWALB','mod_domain_memory' , ywalb , nboun )
          call memory_alloca(memor_dom,'YWALP','mod_domain_memory' , ywalp , nbopo )
          if ( kfl_noslw_ker /= 0_ip) call memory_alloca(memor_dom,'YWALE','mod_domain_memory' , ywale , nelem )
       else
          call memory_deallo(memor_dom,'YWALB','mod_domain_memory' , ywalb )
          call memory_deallo(memor_dom,'YWALP','mod_domain_memory' , ywalp )
          call memory_deallo(memor_dom,'YWALE','mod_domain_memory' , ywale )          
       end if

    case( 'YSCAB, YSCAP' , 'VARIABLE SCALE DISTANCE' )
       !
       ! Arrays for variable wall distance for machine learning
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'YSCAB','mod_domain_memory' , yscab , nboun )
          call memory_alloca(memor_dom,'YSCAP','mod_domain_memory' , yscap , nbopo )
       else
          call memory_deallo(memor_dom,'YSCAB','mod_domain_memory' , yscab )
          call memory_deallo(memor_dom,'YSCAP','mod_domain_memory' , yscap )
       end if

    case ( 'LGROU_DOM' )
       !
       ! Groups for deflated CG
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'LGROU_DOM','mod_domain_memory' , lgrou_dom , npoin )
       else
          call memory_deallo(memor_dom,'LGROU_DOM','mod_domain_memory' , lgrou_dom )          
       end if

    case ( 'XFIEL % A' , 'FIELDS' )
       !
       ! Fields values: allocate
       !../../Sources/kernel/mpio/def_mpio.f90
       ! Master should not allocate these.
       !
       if( present(NUMBER1) ) then
          if( NUMBER1 == 0 ) then
             ifiel1 = 1
             ifiel2 = nfiel
          else
             ifiel1 = NUMBER1
             ifiel2 = NUMBER1
          end if
       else
          ifiel1 = 1
          ifiel2 = nfiel
       end if
       
       do ifiel = ifiel1,ifiel2                
          if( itask == 0 ) then                
             if ( (kfl_field(6,ifiel) == 1) .AND. (.not. mpio_config%output%post_process%export_only) ) then !ondemand 
                nsteps = nsteps_fiel_ondemand
             else
                nsteps = kfl_field(4,ifiel)
             end if
             
             if(      kfl_field(2,ifiel) == NELEM_TYPE ) then
                call memory_alloca(memor_dom,'XFIEL % A','mod_domain_memory' , xfiel(ifiel) % a , kfl_field(1,ifiel) , nelem , nsteps )
             else if( kfl_field(2,ifiel) == NPOIN_TYPE ) then  
                call memory_alloca(memor_dom,'XFIEL % A','mod_domain_memory' , xfiel(ifiel) % a , kfl_field(1,ifiel) , npoin , nsteps )
             else if( kfl_field(2,ifiel) == NBOUN_TYPE ) then                                                                        
                call memory_alloca(memor_dom,'XFIEL % A','mod_domain_memory' , xfiel(ifiel) % a , kfl_field(1,ifiel) , nbou1 , nsteps )
             end if
          else
             call memory_deallo(memor_dom,'XFIEL % A','mod_domain_memory'  , xfiel(ifiel) % a )
          end if
       end do

    case ( 'VMASS' )
       !
       ! VMASS
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'VMASS','mod_domain_memory' , vmass , npoin )
       else
          call memory_deallo(memor_dom,'VMASS','mod_domain_memory' , vmass )
       end if

    case ( 'VMASC' )
       !
       ! VMASC
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'VMASC','mod_domain_memory' , vmasc , npoin )
       else
          call memory_deallo(memor_dom,'VMASC','mod_domain_memory' , vmasc )
       end if

    case ( 'DMASS' )
       !
       ! DMASS
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'DMASS','mod_domain_memory' , dmass , npoin )
       else
          call memory_deallo(memor_dom,'DMASS','mod_domain_memory' , dmass )
       end if

    case ( 'CMASS' )
       !
       ! CMASS
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'CMASS','mod_domain_memory' , cmass , nzdom )
       else
          call memory_deallo(memor_dom,'CMASS','mod_domain_memory' , cmass )
       end if

    case ( 'CMASS_WEIGHTED' )
       !
       ! CMASS_WEIGHTED
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'CMASS_WEIGHTED','mod_domain_memory' , cmass_weighted , nzdom )
       else
          call memory_deallo(memor_dom,'CMASS_WEIGHTED','mod_domain_memory' , cmass_weighted )
       end if

    case ( 'MASS_RIGHT_FACT' )
       !
       ! MASS_RIGHT_FACT
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'MASS_RIGHT_FACT','mod_domain_memory', mass_right_fact , nzdom )
       else
          call memory_deallo(memor_dom,'MASS_RIGHT_FACT','mod_domain_memory', mass_right_fact )
       end if

    case( 'LMATE' )
       !
       ! Allocate materials 
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'LMATE','mod_domain_memory' , lmate , nelem )
       else
          call memory_deallo(memor_dom,'LMATE','mod_domain_memory' , lmate )
       end if

    case( 'MATERIALS' )
       !
       ! Allocate materials: NODEMAT, LMATN and NMATN
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'LMATN'  ,'mod_domain_memory' , lmatn , nmate )
          call memory_alloca(memor_dom,'NMATN'  ,'mod_domain_memory' , nmatn , nmate )
          call memory_alloca(memor_dom,'NODEMAT','mod_domain_memory' , nodemat,npoin )
          do ipoin = 1,npoin
             nodemat(ipoin) = -1_ip
          end do
       else
          call memory_deallo(memor_dom,'LMATN'  ,'mod_domain_memory' , lmatn   )
          call memory_deallo(memor_dom,'NMATN'  ,'mod_domain_memory' , nmatn   )
          call memory_deallo(memor_dom,'NODEMAT','mod_domain_memory' , nodemat )          
       end if

    case( 'LMATN % L' )
       !
       ! LMATN(IMATE) % L
       !
       if( itask == 0 ) then
          if( present(NUMBER1) .and. present(NUMBER2) ) then 
             call memory_alloca(memor_dom,'LMATN % L','mod_domain_memory' , lmatn(NUMBER1)%l , NUMBER2 )
          else
             call runend('MOD_DOMAIN_MEMORY: MISSING ARGUMENTS')
          end if
       else
          if( present(NUMBER1) ) then 
             call memory_deallo(memor_dom,'LMATN % L','mod_domain_memory' , lmatn(NUMBER1)%l )
          else
             call runend('MOD_DOMAIN_MEMORY: MISSING ARGUMENT')
          end if
       end if

    case( 'PELEL' )
       !
       ! PELEL
       !
       if( itask == 0 ) then
          kfl_pelel = 1
          call memory_alloca(memor_dom,'PELEL','mod_domain_memory' , pelel , nelem+1_ip )
       else
          kfl_pelel = 0
          call memory_deallo(memor_dom,'PELEL','mod_domain_memory' , pelel )
       end if

    case( 'LELEL' )
       !
       ! LELEL
       !
       if( itask == 0 ) then
          if( present(NUMBER1) ) then             
             call memory_alloca(memor_dom,'LELEL','mod_domain_memory' , lelel , NUMBER1 )
          else
             call memory_alloca(memor_dom,'LELEL','mod_domain_memory' , lelel , nedge )
          end if
       else
          call memory_deallo(memor_dom,'LELEL','mod_domain_memory' , lelel )
       end if

    case( 'LPERI' )
       !
       ! LPERI: periodicity
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'LPERI','mod_domain_memory' , lperi , 2_ip , nperi )
       else
          call memory_deallo(memor_dom,'LPERI','mod_domain_memory' , lperi )
       end if

    case( 'LEZDO AND LBZDO' )
       !
       ! LEZDO and LBZDO
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'LEZDO','mod_domain_memory' , lezdo , mnode , mnode , nelem )
          call memory_alloca(memor_dom,'LBZDO','mod_domain_memory' , lbzdo , mnodb , mnodb , nbou1 )
       else
          call memory_deallo(memor_dom,'LEZDO','mod_domain_memory' , lezdo )
          call memory_deallo(memor_dom,'LBZDO','mod_domain_memory' , lbzdo )          
       end if

    case ( 'I_AM_IN_ZONE AND I_AM_IN_SUBD' )
       !
       ! Which zones and subdomain I'm in   
       !     
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'I_AM_IN_ZONE','mod_domain_memory',I_AM_IN_ZONE,nzone+1_ip,'INITIALIZE',lboun=0_ip)
          call memory_alloca(memor_dom,'I_AM_IN_SUBD','mod_domain_memory',I_AM_IN_SUBD,nsubd+1_ip,'INITIALIZE',lboun=0_ip)
       else
          call memory_deallo(memor_dom,'I_AM_IN_ZONE','mod_domain_memory',I_AM_IN_ZONE)
          call memory_deallo(memor_dom,'I_AM_IN_SUBD','mod_domain_memory',I_AM_IN_SUBD)
       end if

    case ( 'LBONO' )
       !
       ! List of boundary nodes
       !     
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'LBONO','mod_domain_memory',lbono,nbono)
       else
          call memory_deallo(memor_dom,'LBONO','mod_domain_memory',lbono)
       end if

    case ( 'LNNOB' )
       !
       ! Number of nodes per boundary
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'LNNOB','mod_domain_memory' , lnnob ,  nboun )
       else
          call memory_deallo(memor_dom,'LNNOB','mod_domain_memory' , lnnob )
       end if

    case ( 'LMAST' )
       !
       ! Allocate list of masters
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'LMAST','mod_domain_memory' , lmast , npoin )
       else
          call memory_deallo(memor_dom,'LMAST','mod_domain_memory' , lmast )
       end if

    case ( 'LUMMA' )
       !
       ! lumma - for use in dual time step preconditioner
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'LUMMA','mod_domain_memory' , lumma , npoin )      !nzdom * ndime * ndime    ! If I want full matrix
       else
          call memory_deallo(memor_dom,'LUMMA','mod_domain_memory' , lumma )     
       end if

    case ( 'LGAUS' )
       !
       ! Number of Gauss points per elements
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'LGAUS','mod_domain_memory' , lgaus ,  nelem )
       else
          call memory_deallo(memor_dom,'LGAUS','mod_domain_memory' , lgaus )
       end if

    case ( 'LNNOD' )
       !
       ! lnnod
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'LNNOD','mod_domain_memory' , lnnod ,  nelem   )
       else
          call memory_deallo(memor_dom,'LNNOD','mod_domain_memory' , lnnod )
       end if

    case ( 'TIME_FIELD % A' )
       !
       ! TIME FIELD
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'TIME_FIELD % A','mod_domain_memory' , time_field(NUMBER1) % a , kfl_field(4,NUMBER1) )
       else
          call memory_deallo(memor_dom,'TIME_FIELD % A','mod_domain_memory' , time_field(NUMBER1) % a )
       end if

    case ( 'KFL_CODNO' )
       !
       ! KFL_CODNO
       !
       if( itask == 0 ) then
          if( kfl_icodn > 0 ) then 
             call memory_alloca(memor_dom,'KFL_CODNO','mod_domain_memory',kfl_codno,mcono,npoin,INIT_VALUE=mcodb+1_ip)
          end if
       else
          if( kfl_icodn > 0 ) then
             call memory_deallo(memor_dom,'KFL_CODNO','mod_domain_memory',kfl_codno)
          end if
       end if

    case ( 'KFL_CODBO' )
       !
       ! KFL_CODBO
       !
       if( itask == 0 ) then
          if( kfl_icodb > 0 ) then
             call memory_alloca(memor_dom,'KFL_CODBO','mod_domain_memory',kfl_codbo,max(1_ip,nboun))
          end if
       else          
          if( kfl_icodb > 0 ) then
             call memory_deallo(memor_dom,'KFL_CODBO','mod_domain_memory',kfl_codbo)
          end if
       end if

    case ( 'KFL_GEOBO' )
       !
       ! KFL_GEOBO
       !
       if( itask == 0 ) then
          if( kfl_icodb > 0 ) then
             call memory_alloca(memor_dom,'KFL_GEOBO','mod_domain_memory',kfl_geobo,nboun_2)
          end if
       else
          if( kfl_icodb > 0 ) then
             call memory_deallo(memor_dom,'KFL_GEOBO','mod_domain_memory',kfl_geobo)
          end if
       end if

    case ( 'KFL_GEONO' )
       !
       ! KFL_GEONO
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'KFL_GEONO','mod_domain_memory',kfl_geono,nbopo)
       else
          call memory_deallo(memor_dom,'KFL_GEONO','mod_domain_memory',kfl_geono)
       end if

    case ( 'KFL_CODED' )
       !
       ! KFL_CODED
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'KFL_CODED','mod_domain_memory',kfl_coded,2_ip,meshe(ndivi) % nedge)
       else
          call memory_deallo(memor_dom,'KFL_CODED','mod_domain_memory',kfl_coded)
       end if

    case ( 'LESET' )
       !
       ! LESET
       !
       if( itask == 0 ) then
          if( neset > 0 ) then
             call memory_alloca(memor_dom,'LESET','mod_domain_memory',leset,nelem)
          end if
       else
          if( neset > 0 ) then
             call memory_deallo(memor_dom,'LESET','mod_domain_memory',leset)
          end if
       end if

    case ( 'LBSET' )
       !
       ! LBSET
       !
       if( itask == 0 ) then
          if( nbset > 0 ) then
             call memory_alloca(memor_dom,'LBSET','mod_domain_memory',lbset,nboun)
          end if
       else
          if( nbset > 0 ) then
             call memory_deallo(memor_dom,'LBSET','mod_domain_memory',lbset)
          end if
       end if
       
    case ( 'LNSET' )
       !
       ! LNSET
       !
       if( itask == 0 ) then
          if( nnset > 0 ) then
             call memory_alloca(memor_dom,'LNSET','mod_domain_memory',lnset,npoin)
          end if
       else
          if( associated(lnset) ) then
             call memory_deallo(memor_dom,'LNSET','mod_domain_memory',lnset)
          end if
       end if

    case ( 'LNSEC' )
       !
       ! LNSEC
       !
       if( itask == 0 ) then
          if( nnset > 0 ) then
             call memory_alloca(memor_dom,'LNSEC','mod_domain_memory',lnsec,nnset)
          end if
       else
          if( nnset > 0 ) then
             call memory_deallo(memor_dom,'LNSEC','mod_domain_memory',lnsec)
          end if
       end if

    case ( 'LBSEC' )
       !
       ! LBSEC
       !
       if( itask == 0 ) then
          if( nbset > 0 ) then
             call memory_alloca(memor_dom,'LBSEC','mod_domain_memory',lbsec,nbset)
          end if
       else
          if( nbset > 0 ) then
             call memory_deallo(memor_dom,'LBSEC','mod_domain_memory',lbsec)
          end if
       end if

    case ( 'LESEC' )
       !
       ! LESEC
       !
       if( itask == 0 ) then
          if( neset > 0 ) then
             call memory_alloca(memor_dom,'LESEC','mod_domain_memory',lesec,neset)
          end if
       else
          if( neset > 0 ) then
             call memory_deallo(memor_dom,'LESEC','mod_domain_memory',lesec)
          end if
       end if

    case ( 'WITNESS POINTS' , 'SHWIT, DEWIT AND LEWIT' )
       !
       ! SHWIT, DEWIT, LEWIT
       !
       if( itask == 0 ) then
          if( mwitn > 0 ) then
             call memory_alloca(memor_dom,'SHWIT','mod_domain_memory',shwit,mnode,mwitn)
             call memory_alloca(memor_dom,'DEWIT','mod_domain_memory',dewit,ndime,mnode,mwitn)
             call memory_alloca(memor_dom,'LEWIT','mod_domain_memory',lewit,mwitn)
          end if
       else
          if( mwitn > 0 ) then
             call memory_deallo(memor_dom,'SHWIT','mod_domain_memory',shwit)
             call memory_deallo(memor_dom,'DEWIT','mod_domain_memory',dewit)
             call memory_deallo(memor_dom,'LEWIT','mod_domain_memory',lewit)
          end if
       end if

    case ( 'COWIT' )
       !
       ! COWIT
       !
       if( itask == 0 ) then
          if( mwitn > 0 ) then
             call memory_alloca(memor_dom,'COWIT','mod_domain_memory',cowit,3_ip,mwitn)
          end if
       else
          if( mwitn > 0 ) then
             call memory_deallo(memor_dom,'COWIT','mod_domain_memory',cowit)
          end if
       end if
       
    case ( 'GEWIT' )
       !
       ! GEWIT
       !
       if( itask == 0 ) then
          if( nwitg > 0 ) then
             allocate(gewit(nwitg))
             do ii = 1,nwitg
                gewit(ii) % param = 0.0_rp
             end do
          end if
       else
          if( nwitg > 0 ) then
             deallocate(gewit)
          end if
       end if
 
    case ( 'WITNESS_MESH' )
       !
       ! WITNESS_MESH
       !
       if( itask == 0 ) then
          if( nwith > 0 ) then
             allocate(witness_mesh(nwith))             
             do ii = 1,nwith
                call witness_mesh(ii) % init()
             end do
          end if
       else
          if( nwith > 0 ) then
             deallocate(witness_mesh)
          end if
       end if
 
    case ( 'COWIT_ORIGI' )
       !
       ! COWIT_ORIGI
       !
       if( itask == 0 ) then
          if( mwitn > 0 ) then
             call memory_alloca(memor_dom,'COWIT_ORIGI','mod_domain_memory',cowit_origi,3_ip,mwitn)
          end if
       else
          if( mwitn > 0 ) then
             call memory_deallo(memor_dom,'COWIT_ORIGI','mod_domain_memory',cowit_origi)
          end if
       end if

    case ( 'LNWIT' )
       !
       ! LNWIT
       !
       if( itask == 0 ) then
          call memory_alloca(memor_dom,'LNWIT','mod_domain_memory',lnwit,nwitn)
       else
          call memory_deallo(memor_dom,'LNWIT','mod_domain_memory',lnwit)
       end if

    case ( 'LOBAS' )
       !
       ! Coordinate systems
       !
       if( itask == 0 ) then
          if( present(NUMBER1) ) then
             allocate(lobas(NUMBER1))
             do ii = 1,NUMBER1
                nullify(lobas(ii) % local_basis)
                lobas(ii) % type  = 0_ip
                lobas(ii) % param = 0.0_rp
             end do
          end if
       else
          deallocate(lobas)
       end if
       
    case ( 'LOBAS % LOCAL_BASIS' )
       !
       ! LOBAS: Coordinate system
       !
       if( present(NUMBER1) .and. associated(lobas) ) then
          if( itask == 0 ) then
             call memory_alloca(memor_dom,'LOBAS % LOCAL_BASIS','mod_domain_memory',lobas(NUMBER1) % local_basis,ndime,ndime,nbopo)
          else
             call memory_deallo(memor_dom,'LOBAS % LOCAL_BASIS','mod_domain_memory',lobas(NUMBER1) % local_basis)
          end if
       else
          call runend('MOD_DOMAIN: INTEGER ARGUMENT IS REQUIRED')
       end if
       
     case default

       call runend('DOMAIN_MEMORY: UNKNWON CASE= '//trim(wvari))

    end select

  end subroutine domain_memory_allocate_deallocate

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    13/12/2019
  !> @brief   Max global numbering
  !> @details Max global numbering
  !>
  !----------------------------------------------------------------------

  function domain_total_node_number() result(npoin_max)

    integer(ip) :: npoin_max

    if( nperi <= -1 ) then
       if( ISEQUEN ) then
          npoin_max = npoin
       else
          if( npoin > 0 ) then
             npoin_max = maxval(lninv_loc)
          else
             npoin_max = 0
          end if
          call PAR_MAX(npoin_max)
       end if
    else
       npoin_max = npoin_own
       call PAR_SUM(npoin_max)
    end if
    
  end function domain_total_node_number

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    13/12/2019
  !> @brief   Number of Gauss points
  !> @details Number of Gauss points
  !>
  !----------------------------------------------------------------------

  pure function domain_number_gauss_points(pelty) result(pgaus)

    integer(ip), intent(in) :: pelty
    integer(ip)             :: pgaus
    integer(ip)             :: qelty

    qelty = abs(pelty)
    if( qelty <= element_end ) then
       pgaus = ngaus(qelty)
    else
       pgaus = 1_ip             
    end if

  end function domain_number_gauss_points
  
end module mod_domain
!> @}
