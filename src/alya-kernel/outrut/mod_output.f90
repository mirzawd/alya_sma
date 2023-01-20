!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @defgroup Output_Toolbox
!> Toolbox for output of meshes, element matrices
!> @{
!> @name    ToolBox for output
!> @file    mod_output.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for output
!> @details ToolBox for output, mainly for debugging
!
!-----------------------------------------------------------------------

module mod_output

  use def_kintyp,                only : ip,rp,lg
  use def_master,                only : fil_pos00
  use def_master,                only : fil_pos01
  use def_master,                only : fil_pos02
  use def_master,                only : lun_pos00
  use def_master,                only : lun_pos01
  use def_master,                only : lun_pos02
  use def_master,                only : fil_pos00_save
  use def_master,                only : fil_pos01_save
  use def_master,                only : fil_pos02_save
  use def_master,                only : INOTSLAVE
  use def_master,                only : intost
  use def_domain,                only : nnode
  use def_domain,                only : memor_dom,nnode
  use def_domain,                only : mesh_type
  use def_elmgeo,                only : element_type
  use mod_memory,                only : memory_alloca
  use mod_memory,                only : memory_deallo
  use mod_memory,                only : memory_size
  use mod_messages,              only : messages_live
  use def_parall,                only : num_repart_par
  use def_parall,                only : kfl_repart_par
  use def_parall,                only : kfl_repart_post_par
  use def_AMR,                   only : kfl_amr
  use def_AMR,                   only : kfl_amr_post
  use def_AMR,                   only : num_amr
  use mod_iofile,                only : iofile_close_unit
  use mod_iofile,                only : iofile_open_unit
  use mod_iofile_basic,          only : iofile_create_directory
  use mod_communications_tools,  only : PAR_BARRIER
  use mod_strings,               only : integer_to_string
  use mod_std
#ifdef __PGI
#define MEMPGI )
#else
#define MEMPGI ,memor=memor_dom)
#endif
  implicit none

  private

  interface output_result_gid_format
     module procedure &
          output_result_gid_format_RP_1,&
          output_result_gid_format_IP_1
  end interface output_result_gid_format

  public :: output_mesh_gid_format
  public :: output_element_system
  public :: output_domain
  public :: output_file_names
  public :: output_open_files
  public :: output_result_gid_format
  public :: output_domain_alya_format
  
contains

  !-----------------------------------------------------------------------
  !> @author  Guillaume Houzeaux
  !> @date    18/09/2012
  !> @brief   Output boudnary mesh
  !> @details Output boudnary mesh in STL format
  !-----------------------------------------------------------------------

  subroutine output_stl()

    use def_kintyp_basic,   only :  ip,rp,i1p
    use def_kintyp_comm,    only :  comm_data_par
    use def_master,         only :  INOTMASTER
    use def_master,         only :  INOTSLAVE
    use def_master,         only :  lun_tempo
    use def_master,         only :  namda
    use def_master,         only :  zeror
    use def_master,         only :  kfl_paral
    use def_master,         only :  intost
    use def_kermod,         only :  kfl_oustl
    use def_domain,         only :  ndime,nboun
    use def_domain,         only :  ltypb,lnnob,coord,npoin
    use def_domain,         only :  lnodb

    use def_domain,         only :  nelem,mnode,mnodb,nelty
    use def_domain,         only :  lnods,lnnod,ltype
    use def_domain,         only :  pelpo,lelpo
    use def_domain,         only :  mface

    use def_elmtyp,         only :  TRI03,QUA04
    use mod_communications, only :  PAR_DEFINE_COMMUNICATOR
    use mod_communications, only :  PAR_GATHER
    use mod_communications, only :  PAR_GATHERV
    use mod_parall,         only :  PAR_CODE_SIZE
    use mod_iofile,         only :  iofile
    use mod_graphs,         only :  graphs_list_faces
    use mod_graphs,         only :  graphs_deallocate_list_faces
    use mod_boundary_coarsening, only: boundary_coarsening,boundary_coarsening_graph
    use mod_messages,       only : livinf
    use mod_iofile,         only : iofile_flush_unit
    implicit none
!    type(comm_data_par), pointer :: commu
    integer(ip)                  :: idime,istl
    integer(ip)                  :: pblty,iboun,inodb,ipoin
    integer(ip)                  :: ipart,ipoi1,ipoi2,ipoi3
    integer(ip)                  :: ifacg,ielem,iface,ielty
    integer(4)                   :: nstl4,nstl4_tot
    real(rp)                     :: vec(3,3),vnor,rboun
    integer(4),          pointer :: nstl4_gat(:)
    real(rp),            pointer :: xstl_gat(:)
    real(rp),            pointer :: xstl(:,:)
    character(150)               :: fil_tempo

    integer(ip)                  :: nfacg
    integer(ip),         pointer :: lfacg(:,:)
    type(i1p),           pointer :: lelfa(:)

    integer(ip)                  :: nboun_coarse          !< Number of boundaries
    integer(ip)                  :: npoin_coarse          !< Number of nodes
    integer(ip),         pointer :: lnodb_coarse(:,:)     !< Boundary connectivity
    real(rp),            pointer :: coord_coarse(:,:)     !< Boundary connectivity

    integer(ip)                  :: nbounThreshold
    integer(ip)                  :: coarseningPoints
!    integer(ip)                  :: maxBoundElems

    select case ( kfl_oustl )

    case ( 1_ip )

       !-------------------------------------------------------------------
       !
       ! Master outputs the geometrical STL using lnodb
       !
       !-------------------------------------------------------------------

       if( ndime == 3 ) then

          call livinf(0_ip,'OUTPUT BOUNDARY MESH IN STL FORMAT',0_ip)
          nullify(nstl4_gat)
          nullify(xstl_gat)
          nullify(xstl)
          nstl4 = 0
          !call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
          if( INOTSLAVE ) then
             !
             ! Open file
             !
             fil_tempo = adjustl(trim(namda))//'.stl'
             call iofile(0_ip,lun_tempo,fil_tempo,'SYSTEM INFO')
          end if
          !
          ! Gather all STL
          !
          allocate( nstl4_gat(0:PAR_CODE_SIZE-1) )
          !
          ! Workers convert Alya format to STL
          !
          if( INOTMASTER ) then
             do iboun = 1,nboun
                pblty = ltypb(iboun)
                if( pblty == TRI03 ) then
                   nstl4  = nstl4 + 1
                else if( pblty == QUA04 ) then
                   nstl4  = nstl4 + 2
                end if
             end do
             allocate( xstl(9,nstl4) )
             nstl4 = 0
             do iboun = 1,nboun
                pblty = ltypb(iboun)
                if( pblty == TRI03 ) then
                   nstl4  = nstl4 + 1
                   idime = 1
                   do inodb = 1,lnnob(iboun)
                      ipoin = lnodb(inodb,iboun)
                      xstl(idime:idime+2,nstl4) = coord(1:ndime,ipoin)
                      idime = idime + 3
                   end do
                else if( pblty == QUA04 ) then
                   nstl4  = nstl4 + 1
                   idime = 1
                   ipoi1 = lnodb(1,iboun)
                   ipoi2 = lnodb(2,iboun)
                   ipoi3 = lnodb(3,iboun)
                   xstl(idime:idime+2,nstl4) = coord(1:ndime,ipoi1) ; idime = idime + 3
                   xstl(idime:idime+2,nstl4) = coord(1:ndime,ipoi2) ; idime = idime + 3
                   xstl(idime:idime+2,nstl4) = coord(1:ndime,ipoi3)
                   nstl4  = nstl4 + 1
                   idime = 1
                   ipoi1 = lnodb(1,iboun)
                   ipoi2 = lnodb(3,iboun)
                   ipoi3 = lnodb(4,iboun)
                   xstl(idime:idime+2,nstl4) = coord(1:ndime,ipoi1) ; idime = idime + 3
                   xstl(idime:idime+2,nstl4) = coord(1:ndime,ipoi2) ; idime = idime + 3
                   xstl(idime:idime+2,nstl4) = coord(1:ndime,ipoi3)
                else
                   call runend('OUTSTL: CANNOT GENERATE STL FOR THIS BOUNDARY TYPE '//element_type(pblty) % name)
                end if
             end do
          end if
          !
          ! Gather all STL
          !
          nstl4 = 9_4*nstl4
          call PAR_GATHER(nstl4,nstl4_gat,'IN MY CODE')

          if( INOTSLAVE ) then
             nstl4_tot = 0
             do ipart = 0,PAR_CODE_SIZE-1
                nstl4_tot = nstl4_tot + nstl4_gat(ipart)
             end do
             allocate( xstl_gat(nstl4_tot) )
          end if

          call PAR_GATHERV(xstl,xstl_gat,nstl4_gat,'IN MY CODE')
          !
          ! Master outputs STL
          !
          if( INOTSLAVE ) then
             write(lun_tempo,'(a)') 'solid mesh_boundary'
             nstl4 = 0
             do ipart = 0,PAR_CODE_SIZE-1
                do istl = 1,nstl4_gat(ipart)/9
                   vec(1,1) = xstl_gat(nstl4+4) - xstl_gat(nstl4+1)
                   vec(2,1) = xstl_gat(nstl4+5) - xstl_gat(nstl4+2)
                   vec(3,1) = xstl_gat(nstl4+6) - xstl_gat(nstl4+3)
                   vec(1,2) = xstl_gat(nstl4+7) - xstl_gat(nstl4+1)
                   vec(2,2) = xstl_gat(nstl4+8) - xstl_gat(nstl4+2)
                   vec(3,2) = xstl_gat(nstl4+9) - xstl_gat(nstl4+3)
                   call vecpro(vec(1,1),vec(1,2),vec(1,3),ndime)
                   call vecnor(vec(1,3),ndime,vnor,2_ip)
                   vec(1:3,3) = vec(1:3,3) / ( vnor + zeror )
                   write(lun_tempo,1)     'facet normal',vec(1:3,3)
                   write(lun_tempo,'(a)') 'outer loop'
                   write(lun_tempo,1)     'vertex',xstl_gat(nstl4+1:nstl4+3)
                   write(lun_tempo,1)     'vertex',xstl_gat(nstl4+4:nstl4+6)
                   write(lun_tempo,1)     'vertex',xstl_gat(nstl4+7:nstl4+9)
                   write(lun_tempo,'(a)') 'endloop'
                   write(lun_tempo,'(a)') 'endfacet'
                   nstl4 = nstl4 + 9
                end do
             end do
             write(lun_tempo,'(a)') 'endsolid mesh_boundary'
             deallocate( nstl4_gat )
             deallocate( xstl_gat  )
             call iofile_flush_unit(lun_tempo)
             close(lun_tempo)
          else
             deallocate( xstl )
          end if

       end if

    case ( 2_ip )

       !-------------------------------------------------------------------
       !
       ! Each slaves output its STL including subdomain interfaces
       !
       !-------------------------------------------------------------------

       if( INOTMASTER ) then
          nullify(lfacg)
          nullify(lelfa)
          call graphs_list_faces(&
               nelem,mnode,mnodb,nelty,mface,lnods,ltype,&
               pelpo,lelpo,nfacg,lfacg,lelfa MEMPGI

          nbounThreshold = 1000
          rboun = real(nboun,rp)
          !maxBoundElems = 100000
          !coarseningPoints = maxBoundElems*(tanh( nboun/dfloat(maxBoundElems) ) +1)/2.0
          if(nboun>nbounThreshold) then!!!!
             !
             if(nboun<10000) then
                coarseningPoints = 2 + 1*(int(tanh(rboun/1000.0_rp-1.0_rp),ip)  +1)/2
             else if(nboun<100000) then
                coarseningPoints = 4 + 2*(int(tanh(rboun/10000.0_rp-1.0_rp),ip) +1)/2
             else
                coarseningPoints = 6 + 4*(int(tanh(rboun/100000.0_rp-1.0_rp),ip)+1)/2
             end if

             !-***default call
             !            call boundary_coarsening_graph(&
             !                nelem,mnode,mnodb,nelty,mface,lnods,lnnod,ltype,&
             !                nnodf,lface,nface,pelpo,lelpo,nfacg,lfacg,lelfa,&
             !                ndime,npoin,coord,&
             !                nboun_coarse,npoin_coarse,lnodb_coarse,coord_coarse)!Default:'divide',10
             !-*** call with optional arguments
             !-** 'divide' requires the factor to compute npoin_coarse=npoin/factor
             call boundary_coarsening_graph(&
                  nelem,mnode,mnodb,nelty,mface,lnods,lnnod,ltype,&
                  pelpo,lelpo,nfacg,lfacg,lelfa,&
                  ndime,npoin,coord,&
                  nboun_coarse,npoin_coarse,lnodb_coarse,coord_coarse,&
                  'divide',coarseningPoints)
             !-** 'npoin_coarse' provides npoin_coarse (approximately)
             !            call boundary_coarsening_graph(&
             !                nelem,mnode,mnodb,nelty,mface,lnods,lnnod,ltype,&
             !                nnodf,lface,nface,pelpo,lelpo,nfacg,lfacg,lelfa,&
             !                ndime,npoin,coord,&
             !                nboun_coarse,npoin_coarse,lnodb_coarse,coord_coarse,&
             !                'npoin_coarse',coarseningPoints)

             fil_tempo = adjustl(trim(namda))//'-'//trim(intost(kfl_paral))//'.stl'
             call iofile(0_ip,lun_tempo,fil_tempo,'SYSTEM INFO')
             write(lun_tempo,'(a)') 'solid subdomain_boundary'
             do iboun = 1,nboun_coarse
                ipoi1    = lnodb_coarse(1,iboun)
                ipoi2    = lnodb_coarse(2,iboun)
                ipoi3    = lnodb_coarse(3,iboun)
                vec(1,1) = coord_coarse(1,ipoi2) - coord_coarse(1,ipoi1)
                vec(2,1) = coord_coarse(2,ipoi2) - coord_coarse(2,ipoi1)
                vec(3,1) = coord_coarse(3,ipoi2) - coord_coarse(3,ipoi1)
                vec(1,2) = coord_coarse(1,ipoi3) - coord_coarse(1,ipoi1)
                vec(2,2) = coord_coarse(2,ipoi3) - coord_coarse(2,ipoi1)
                vec(3,2) = coord_coarse(3,ipoi3) - coord_coarse(3,ipoi1)
                call vecpro(vec(1,1),vec(1,2),vec(1,3),ndime)
                call vecnor(vec(1,3),ndime,vnor,2_ip)
                vec(1:3,3) = vec(1:3,3) / ( vnor + zeror )
                write(lun_tempo,1)     'facet normal',vec(1:3,3)
                write(lun_tempo,'(a)') 'outer loop'
                write(lun_tempo,1)     'vertex',coord_coarse(1:3,ipoi1)
                write(lun_tempo,1)     'vertex',coord_coarse(1:3,ipoi2)
                write(lun_tempo,1)     'vertex',coord_coarse(1:3,ipoi3)
                write(lun_tempo,'(a)') 'endloop'
                write(lun_tempo,'(a)') 'endfacet'
             end do
             deallocate( coord_coarse )
             deallocate( lnodb_coarse )

          else!!!!

             fil_tempo = adjustl(trim(namda))//'-'//trim(intost(kfl_paral))//'.stl'
             call iofile(0_ip,lun_tempo,fil_tempo,'SYSTEM INFO')
             write(lun_tempo,'(a)') 'solid subdomain_boundary'
             do ifacg = 1,nfacg
                if( lfacg(2,ifacg) == 0 ) then
                   iface = lfacg(3,ifacg)
                   ielem = lfacg(1,ifacg)
                   ielty = ltype(ielem)
                   pblty = element_type(ielty) % type_faces(iface) 
                   if( pblty == TRI03 ) then
                      ipoi1    = lnods(element_type(ielty) % list_faces(1,iface),ielem)   
                      ipoi2    = lnods(element_type(ielty) % list_faces(2,iface),ielem)
                      ipoi3    = lnods(element_type(ielty) % list_faces(3,iface),ielem)
                      vec(1,1) = coord(1,ipoi2) - coord(1,ipoi1)
                      vec(2,1) = coord(2,ipoi2) - coord(2,ipoi1)
                      vec(3,1) = coord(3,ipoi2) - coord(3,ipoi1)
                      vec(1,2) = coord(1,ipoi3) - coord(1,ipoi1)
                      vec(2,2) = coord(2,ipoi3) - coord(2,ipoi1)
                      vec(3,2) = coord(3,ipoi3) - coord(3,ipoi1)
                      call vecpro(vec(1,1),vec(1,2),vec(1,3),ndime)
                      call vecnor(vec(1,3),ndime,vnor,2_ip)
                      vec(1:3,3) = vec(1:3,3) / ( vnor + zeror )
                      write(lun_tempo,1)     'facet normal',vec(1:3,3)
                      write(lun_tempo,'(a)') 'outer loop'
                      write(lun_tempo,1)     'vertex',coord(1:3,ipoi1)
                      write(lun_tempo,1)     'vertex',coord(1:3,ipoi2)
                      write(lun_tempo,1)     'vertex',coord(1:3,ipoi3)
                      write(lun_tempo,'(a)') 'endloop'
                      write(lun_tempo,'(a)') 'endfacet'
                   else if( pblty == QUA04 ) then
                      ipoi1    = lnods(element_type(ielty) % list_faces(1,iface),ielem)
                      ipoi2    = lnods(element_type(ielty) % list_faces(2,iface),ielem)
                      ipoi3    = lnods(element_type(ielty) % list_faces(3,iface),ielem)
                      vec(1,1) = coord(1,ipoi2) - coord(1,ipoi1)
                      vec(2,1) = coord(2,ipoi2) - coord(2,ipoi1)
                      vec(3,1) = coord(3,ipoi2) - coord(3,ipoi1)
                      vec(1,2) = coord(1,ipoi3) - coord(1,ipoi1)
                      vec(2,2) = coord(2,ipoi3) - coord(2,ipoi1)
                      vec(3,2) = coord(3,ipoi3) - coord(3,ipoi1)
                      call vecpro(vec(1,1),vec(1,2),vec(1,3),ndime)
                      call vecnor(vec(1,3),ndime,vnor,2_ip)
                      vec(1:3,3) = vec(1:3,3) / ( vnor + zeror )
                      write(lun_tempo,1)     'facet normal',vec(1:3,3)
                      write(lun_tempo,'(a)') 'outer loop'
                      write(lun_tempo,1)     'vertex',coord(1:3,ipoi1)
                      write(lun_tempo,1)     'vertex',coord(1:3,ipoi2)
                      write(lun_tempo,1)     'vertex',coord(1:3,ipoi3)
                      write(lun_tempo,'(a)') 'endloop'
                      write(lun_tempo,'(a)') 'endfacet'
                      ipoi1    = lnods(element_type(ielty) % list_faces(1,iface),ielem)
                      ipoi2    = lnods(element_type(ielty) % list_faces(3,iface),ielem)
                      ipoi3    = lnods(element_type(ielty) % list_faces(4,iface),ielem)
                      vec(1,1) = coord(1,ipoi2) - coord(1,ipoi1)
                      vec(2,1) = coord(2,ipoi2) - coord(2,ipoi1)
                      vec(3,1) = coord(3,ipoi2) - coord(3,ipoi1)
                      vec(1,2) = coord(1,ipoi3) - coord(1,ipoi1)
                      vec(2,2) = coord(2,ipoi3) - coord(2,ipoi1)
                      vec(3,2) = coord(3,ipoi3) - coord(3,ipoi1)
                      call vecpro(vec(1,1),vec(1,2),vec(1,3),ndime)
                      call vecnor(vec(1,3),ndime,vnor,2_ip)
                      vec(1:3,3) = vec(1:3,3) / ( vnor + zeror )
                      write(lun_tempo,1)     'facet normal',vec(1:3,3)
                      write(lun_tempo,'(a)') 'outer loop'
                      write(lun_tempo,1)     'vertex',coord(1:3,ipoi1)
                      write(lun_tempo,1)     'vertex',coord(1:3,ipoi2)
                      write(lun_tempo,1)     'vertex',coord(1:3,ipoi3)
                      write(lun_tempo,'(a)') 'endloop'
                      write(lun_tempo,'(a)') 'endfacet'
                   else
                      call runend('OUTSTL: CANNOT GENERATE STL FOR THIS BOUNDARY TYPE '//element_type(pblty) % name)
                   end if
                end if
             end do

          end if!!!!!

          write(lun_tempo,'(a)') 'endsolid subdomain_boundary'
          flush(lun_tempo)
          close(lun_tempo)
          call graphs_deallocate_list_faces(lfacg,lelfa MEMPGI
       end if

    end select

1   format(a,3(1x,e13.6),1x,i7)

  end subroutine output_stl


  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    05/01/2018
  !> @brief   Mesh in Gid format
  !> @details Output a mesh or a submesh in GiD format
  !>          1. Offsets can be used to call this subroutine recursively,
  !>             in order to append meshes: NODE_OFFSET, ELEMENT_OFFSET
  !>          2. To output a submesh, you can mark nodes with
  !>             MARK_NPOIN_OPT or elements with mark_nelem_opt
  !>          3. Boundary output can be disabled using OUTPUT_BOUNDARY='OFF'
  !
  !-----------------------------------------------------------------------

  subroutine output_mesh_gid_format(&
       meshe,title,lunit,mark_npoin_opt,mark_nelem_opt,&
       RESULT_INT,NODE_OFFSET,ELEMENT_OFFSET,OUTPUT_BOUNDARY,&
       MARK_ELEMENTS,RESULT_ELEM_INT,RESULT_UNIT)

    type(mesh_type), intent(inout)                 :: meshe               !< Mesh type
    character(*),    intent(in)                    :: title               !< Title of the mesh
    integer(ip),     intent(in)                    :: lunit               !< Output unit
    logical(lg),     intent(in), pointer, optional :: mark_npoin_opt(:)   !< If a submesh is required
    logical(lg),     intent(in), pointer, optional :: mark_nelem_opt(:)   !< If a submesh is required
    integer(ip),     intent(in), pointer, optional :: RESULT_INT(:)       !< Integer results
    integer(ip),     intent(in),          optional :: NODE_OFFSET         !< Offset for node numbering
    integer(ip),     intent(in),          optional :: ELEMENT_OFFSET      !< Offset for element numbering
    character(*),    intent(in),          optional :: OUTPUT_BOUNDARY     !< Output boundary mesh
    integer(ip),     intent(in), pointer, optional :: MARK_ELEMENTS(:)    !< Mark elements
    integer(ip),     intent(in), pointer, optional :: RESULT_ELEM_INT(:)  !< Integer results
    integer(ip),     intent(in),          optional :: RESULT_UNIT         !< Result unit

    integer(ip)                                    :: iesto,iesta,ibsta
    integer(ip)                                    :: ifirs,inode
    integer(ip)                                    :: ielty,pnodb
    integer(ip)                                    :: mnode,npoin,nelem
    integer(ip)                                    :: nboun,iboun
    integer(ip)                                    :: mnodb
    integer(ip)                                    :: ndime,ipoin,ielem
    integer(ip)                                    :: ipoin_offset
    integer(ip)                                    :: ielem_offset
    integer(ip)                                    :: lun_result_unit
    character(150)                                 :: dumml
    integer(ip),   pointer                         :: lnods(:,:)
    integer(ip),   pointer                         :: ltype(:)
    integer(ip),   pointer                         :: lnodb(:,:)
    integer(ip),   pointer                         :: ltypb(:)
    real(rp),      pointer                         :: coord(:,:)
    logical(lg),   pointer                         :: mark_nelem(:)
    logical(lg),   pointer                         :: mark_nboun(:)
    logical(lg)                                    :: output_type
    logical(lg)                                    :: boundary
    integer(ip),   pointer                         :: permu_npoin(:)
    !
    ! Nullify
    !
    if( meshe % nelem <= 0 ) return
    nullify(lnods)
    nullify(ltype)
    nullify(lnodb)
    nullify(ltypb)
    nullify(coord)
    nullify(mark_nelem)
    nullify(mark_nboun)
    nullify(permu_npoin)
    !
    ! Output boundary
    !
    boundary = .true.
    if( present(OUTPUT_BOUNDARY) ) then
       if( trim(OUTPUT_BOUNDARY) == 'OFF' .or. trim(OUTPUT_BOUNDARY) == 'NO' ) then
          boundary = .false.
       end if
    end if
    !
    ! Submesh dimensions
    !
    mnode = meshe % mnode
    mnodb = meshe % mnodb
    ndime = meshe % ndime
    !
    ! Offsets
    !
    if( present(NODE_OFFSET) ) then
       ipoin_offset = NODE_OFFSET
    else
       ipoin_offset = 0
    end if
    if( present(ELEMENT_OFFSET) ) then
       ielem_offset = ELEMENT_OFFSET
    else
       ielem_offset = 0
    end if

    !--------------------------------------------------------------------
    !
    ! Compute submeshes
    !
    !--------------------------------------------------------------------

    if( present(mark_npoin_opt) ) then
       !
       ! Submesh where nodes are marked: mark associated elements and boundaries
       !
       call memory_alloca(memor_dom,'MARK_NELEM','output_mesh_gid_format',mark_nelem,meshe % nelem)
       call memory_alloca(memor_dom,'MARK_NBOUN','output_mesh_gid_format',mark_nboun,meshe % nboun)
       do ielem = 1,meshe % nelem
          loop_inode: do inode = 1,meshe % lnnod(ielem)
             ipoin = meshe % lnods(inode,ielem)
             if( mark_npoin_opt(ipoin) ) then
                mark_nelem(ielem) = .true.
                exit loop_inode
             end if
          end do loop_inode
       end do
       do iboun = 1,meshe % nboun
          pnodb = meshe % lnnob(iboun)
          ielem = meshe % lelbo(iboun)
          if( mark_nelem(ielem) ) mark_nboun(iboun) = .true.
       end do
    else if( present(mark_nelem_opt) ) then
       mark_nelem => mark_nelem_opt
    end if

    if( present(mark_npoin_opt) .or. present(mark_nelem_opt) ) then
       !
       ! Submesh where elements are marked
       !
       call memory_alloca(memor_dom,'PERMU_NPOIN','output_mesh_gid_format',permu_npoin,meshe % npoin)

       nelem = count(mark_nelem(1:meshe % nelem),KIND=ip)
       nboun = count(mark_nboun(1:meshe % nboun),KIND=ip)
       call memory_alloca(memor_dom,'LNODS','output_mesh_gid_format',lnods,mnode,nelem)
       call memory_alloca(memor_dom,'LTYPE','output_mesh_gid_format',ltype,nelem)
       call memory_alloca(memor_dom,'LNODB','output_mesh_gid_format',lnodb,mnodb,nboun)
       call memory_alloca(memor_dom,'LTYPB','output_mesh_gid_format',ltypb,nboun)
       nelem = 0
       do ielem = 1,meshe % nelem
          if( mark_nelem(ielem) ) then
             nelem = nelem + 1
             lnods(:,nelem) = meshe % lnods(:,ielem)
             ltype(nelem)   = meshe % ltype(ielem)
          end if
       end do
       nboun = 0
       do iboun = 1,meshe % nboun
          if( mark_nboun(iboun) ) then
             nboun = nboun + 1
             lnodb(:,nboun) = meshe % lnodb(:,iboun)
             ltypb(nboun)   = meshe % ltypb(iboun)
          end if
       end do

       npoin = 0
       do ielem = 1,meshe % nelem
          if( mark_nelem(ielem) ) then
             do inode = 1,meshe % lnnod(ielem)
                ipoin = meshe % lnods(inode,ielem)
                permu_npoin(ipoin) = 1
             end do
          end if
       end do
       npoin = count(permu_npoin==1,KIND=ip)
       call memory_alloca(memor_dom,'COORD','output_mesh_gid_format',coord,ndime,npoin)

       npoin = 0
       do ipoin = 1,meshe % npoin
          if( permu_npoin(ipoin) == 1 ) then
             npoin = npoin + 1
             permu_npoin(ipoin) = npoin
             coord(1:ndime,npoin) = meshe % coord(1:ndime,ipoin)
          end if
       end do
       do ielem = 1,nelem
          do inode = 1,mnode
             ipoin = lnods(inode,ielem)
             if( ipoin > 0 ) lnods(inode,ielem) = permu_npoin(ipoin)
          end do
       end do
       do iboun = 1,nboun
          do inode = 1,mnodb
             ipoin = lnodb(inode,iboun)
             if( ipoin > 0 ) lnodb(inode,iboun) = permu_npoin(ipoin)
          end do
       end do

    else
       !
       ! Whole mesh
       !
       npoin =  meshe % npoin
       nelem =  meshe % nelem
       nboun =  meshe % nboun
       lnods => meshe % lnods
       ltypb => meshe % ltypb
       lnodb => meshe % lnodb
       ltype => meshe % ltype
       coord => meshe % coord

    end if
    !
    ! Element range
    !
    if(      ndime == 1 ) then
       ibsta =  1
       iesta =  2
       iesto =  9
    else if( ndime == 2 ) then
       ibsta =  2
       iesta = 10
       iesto = 29
    else if( ndime == 3 ) then
       ibsta = 10
       iesta = 30
       iesto = 50
    end if

    !--------------------------------------------------------------------
    !
    ! Output mesh
    !
    !--------------------------------------------------------------------

    ifirs = 0
    do ielty = ibsta,iesto

       output_type = .false.
       if( any(abs(ltype)==ielty) ) output_type = .true.
       if( nboun > 0 ) then
          if( any(abs(ltypb)==ielty) ) output_type = .true.
       end if

       if( output_type ) then

          dumml = adjustl(trim(title)//'_'//element_type(ielty) % name)
          !
          ! Header
          !
          write(lunit,1)&
               adjustl(trim(dumml)),ndime,&
               adjustl(trim(element_type(ielty) % nametopo)),nnode(ielty)
          !
          ! Coordinates
          !
          if( ifirs == 0 .and. npoin > 0 ) then
             ifirs = 1
             write(lunit,2) 'coordinates'
             if( ndime == 1 ) then
                do ipoin = 1,npoin
                   write(lunit,3) ipoin+ipoin_offset,coord(1,ipoin),0.0_rp
                end do
             else
                do ipoin = 1,npoin
                   write(lunit,3) ipoin+ipoin_offset,coord(1:ndime,ipoin)
                end do
             end if
             write(lunit,2) 'end coordinates'
          end if

          write(lunit,2) 'elements'

          if( ielty >= iesta ) then
             !
             ! Element connectivity
             !
             if( present(MARK_ELEMENTS) ) then
                do ielem = 1,nelem
                   if( abs(ltype(ielem)) == ielty ) then
                      write(lunit,4) ielem+ielem_offset,(lnods(inode,ielem),inode=1,nnode(ielty),MARK_ELEMENTS(ielem))
                   end if
                end do
             else
                do ielem = 1,nelem
                   if( abs(ltype(ielem)) == ielty ) then
                      write(lunit,4) ielem+ielem_offset,(lnods(inode,ielem),inode=1,nnode(ielty))
                   end if
                end do
             end if

          else if( boundary ) then
             !
             ! Boundary connectivity
             !
             do iboun = 1,nboun
                if( abs(ltypb(iboun)) == ielty ) then
                   write(lunit,4) iboun+ielem_offset,(lnodb(inode,iboun),inode=1,nnode(ielty))
                end if
             end do
          end if
          write(lunit,2) 'end elements'
          write(lunit,2) ''

       end if

    end do
    !
    ! Results
    !
    if( present(RESULT_UNIT) ) then
       lun_result_unit = RESULT_UNIT
    else
       lun_result_unit = lunit
    end if

    if( present(RESULT_INT) ) then
       if( present(mark_npoin_opt) .or. present(mark_nelem_opt) ) then
          call output_result_gid_format(lun_result_unit,meshe % npoin,RESULT_INT,'GENERIC',permu_npoin)
       else
          call output_result_gid_format(lun_result_unit,meshe % npoin,RESULT_INT,'GENERIC')
       end if
    end if
    if( present(RESULT_ELEM_INT) ) then
       call output_result_gid_format(lun_result_unit,meshe % nelem,RESULT_ELEM_INT,'GENERIC',wherein='ON ELEMENTS')
    end if
    !
    ! Deallocate
    !
    call memory_deallo(memor_dom,'PERMU_NPOIN','output_mesh_gid_format',permu_npoin)

    if( present(mark_npoin_opt) ) then

       call memory_deallo(memor_dom,'MARK_NELEM','output_mesh_gid_format',mark_nelem)
       call memory_deallo(memor_dom,'MARK_NBOUN','output_mesh_gid_format',mark_nboun)

    end if

    if( present(mark_npoin_opt) .or. present(mark_nelem_opt) ) then

       call memory_deallo(memor_dom,'LNODS','output_mesh_gid_format',lnods)
       call memory_deallo(memor_dom,'LTYPE','output_mesh_gid_format',ltype)
       call memory_deallo(memor_dom,'LNODB','output_mesh_gid_format',lnodb)
       call memory_deallo(memor_dom,'LTYPB','output_mesh_gid_format',ltypb)
       call memory_deallo(memor_dom,'COORD','output_mesh_gid_format',coord)

    end if

1   format('MESH ',a,' dimension ',i1,' Elemtype ',a,' Nnode ',i2)
2   format(a)
3   format(i9, 3(1x,e16.8e3))
4   format(i9,50(1x,i9))

  end subroutine output_mesh_gid_format

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    05/01/2018
  !> @brief   Result in Gid format
  !> @details Output a result with possible permutation
  !
  !-----------------------------------------------------------------------

  subroutine output_result_gid_format_RP_1(lunit,nenti,RESULT_INT,NAME,permu_nenti,wherein,TIME_STEP,RESET_POSTPROCESS)

    integer(ip),     intent(in)                    :: lunit
    integer(ip),     intent(in)                    :: nenti
    real(rp),        intent(in), pointer           :: RESULT_INT(:)       !< Integer results
    character(*),    intent(in),          optional :: NAME
    integer(ip),     intent(in), pointer, optional :: permu_nenti(:)
    character(len=*),intent(in),          optional :: wherein
    real(rp),        intent(in),          optional :: TIME_STEP
    logical(lg),     intent(in),          optional :: RESET_POSTPROCESS
    integer(ip)                                    :: ipoin,kpoin,lmax
    character(10)                                  :: my_name
    integer(ip),                          save     :: lunit_old=0
    integer(ip)                                    :: entity_type
    real(rp)                                       :: rstep
    logical(lg)                                    :: if_reset_postprocess

    if( present(wherein) ) then
       if( trim(wherein) == 'ON NODES' ) then
          entity_type = 0
       else if( trim(wherein) == 'ON ELEMENTS' ) then
          entity_type = 1
       else
          call runend('OUTPUT_RESULT_GID_FORMAT: UNKNOWN ENTITY')
       end if
    else
       entity_type = 0
    end if

    if( present(RESET_POSTPROCESS) ) then
       if_reset_postprocess = RESET_POSTPROCESS
    else
       if_reset_postprocess = .false.
    end if

    if( lunit /= lunit_old .or. if_reset_postprocess ) then
       write(lunit,1) 'GiD Post Results File 1.0'
       write(lunit,1) ' '
    end if

    if( present(NAME) ) then
       lmax = min(10,len(NAME))
       my_name = trim(NAME(1:lmax))
    else
       my_name = 'GENERIC'
    end if

    if( present(TIME_STEP) ) then
       rstep = TIME_STEP
    else
       rstep = 1.0_rp
    end if

    if(      entity_type == 0 ) then
       write(lunit,2) trim(my_name),'ALYA',rstep,'Scalar','OnNodes'
    else if( entity_type == 1 ) then
       write(lunit,'(a)') 'GaussPoints GP Elemtype Triangle'
       write(lunit,'(a)') 'Number of Gauss Points:   1'
       write(lunit,'(a)') 'Natural Coordinates: Internal'
       write(lunit,'(a)') 'End GaussPoints'
       write(lunit,'(a)') 'GaussPoints GP Elemtype Quadrilateral'
       write(lunit,'(a)') 'Number of Gauss Points:   1'
       write(lunit,'(a)') 'Natural Coordinates: Internal'
       write(lunit,'(a)') 'End GaussPoints'
       write(lunit,2) trim(my_name),'ALYA',rstep,'Scalar','OnGaussPoints GP'
    end if

    write(lunit,3) trim(my_name)


    write(lunit,1) 'Values'

    if( present(permu_nenti) ) then
       do ipoin = 1,nenti
          if( permu_nenti(ipoin) > 0 ) then
             kpoin = permu_nenti(ipoin)
             write(lunit,6) kpoin,RESULT_INT(kpoin)
          end if
       end do
    else
       do ipoin = 1,nenti
          write(lunit,6) ipoin,RESULT_INT(ipoin)
       end do
    end if

    write(lunit,1) 'End values'

    lunit_old = lunit

1   format(a)
2   format('Result ',a,' ',a,' ',e15.8,' ',a,' ',a)
3   format('ComponentNames ',a)
6   format(i9, 3(1x,e13.6))

  end subroutine output_result_gid_format_RP_1

  subroutine output_result_gid_format_IP_1(lunit,nenti,RESULT_INT,NAME,permu_nenti,wherein,TIME_STEP,RESET_POSTPROCESS)

    integer(ip),     intent(in)                    :: lunit
    integer(ip),     intent(in)                    :: nenti
    integer(ip),     intent(in), pointer           :: RESULT_INT(:)       !< Integer results
    character(*),    intent(in),          optional :: NAME
    integer(ip),     intent(in), pointer, optional :: permu_nenti(:)
    character(len=*),intent(in),          optional :: wherein
    real(rp),        intent(in),          optional :: TIME_STEP
    logical(lg),     intent(in),          optional :: RESET_POSTPROCESS
    real(rp),                    pointer           :: result_rea(:)
    integer(ip)                                    :: ipoin

    allocate(result_rea(nenti))
    do ipoin = 1,nenti
       result_rea(ipoin) = real(RESULT_INT(ipoin),rp)
    end do
    call output_result_gid_format_RP_1(lunit,nenti,RESULT_rea,NAME,permu_nenti,wherein,TIME_STEP,RESET_POSTPROCESS)
    deallocate(result_rea)

  end subroutine output_result_gid_format_IP_1

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux version of nastin (adapted by ECR)
  !> @brief   Output system
  !> @details Output elemental system
  !>
  !-----------------------------------------------------------------------

  subroutine output_element_system(kdime,pnode,elmat,elrhs)

    use def_master, only :  intost

    implicit none

    integer(ip), intent(in)           :: kdime
    integer(ip), intent(in)           :: pnode
    ! Element matrices
    real(rp),    intent(in)           :: elmat(kdime*pnode,kdime*pnode)
    real(rp),    intent(in)           :: elrhs(kdime,pnode)

    integer(ip),   parameter          :: lreal=9
    character(13)                     :: FMT1
    integer(ip)                       :: inode,idime,idofn,jnode,jdime,jdofn
    character(3)                      :: vanam

    !
    ! Format
    !
    FMT1 = '(1x,e'//trim(intost(lreal))//'.'//trim(intost(lreal-7))//')'
    !
    ! First line
    !
    write(99,'(a)',advance='no') '+-'
    do idofn = 1,kdime*pnode
       do idime = 1,lreal+1
          write(99,'(a)',advance='no') ' '
       end do
    end do
    write(99,'(a)',advance='no') ' '
    do idofn = 1,pnode
       do idime = 1,lreal+1
          write(99,'(a)',advance='no') ' '
       end do
    end do

    write(99,'(a)',advance='no') '-+'
    write(99,'(a)',advance='no') '  +-   -+'

    write(99,'(a)',advance='no') '     +-'
    do idime = 1,lreal
       write(99,'(a)',advance='no') ' '
    end do
    write(99,'(a)',advance='no') '-+'

    write(99,*)

    do inode = 1,pnode
       do idime = 1,kdime
          !
          ! Auu
          !
          write(99,'(a)',advance='no') '|'
          idofn = (inode-1)*kdime+idime
          do jnode = 1,pnode
             do jdime = 1,kdime
                jdofn = (jnode-1)*kdime+jdime
                write(99,FMT1,advance='no') elmat(idofn,jdofn)
             end do
          end do
          !
          ! ui
          !
          write(99,'(a)',advance='no') ' |'
          if( idime == 1 ) then
             vanam = 'u'
          else if( idime == 2 ) then
             vanam = 'v'
          else if( idime == 3 ) then
             vanam = 'w'
          end if
          vanam = trim(vanam) // trim(intost(inode))
          write(99,'(a)',advance='no') '  | '//vanam//' |'
          !
          ! bu
          !
          write(99,'(a)',advance='no') '     |'
          write(99,FMT1,advance='no') elrhs(idime,inode)
          write(99,'(a)',advance='no') ' |'

          write(99,*)
       end do
    end do
    !
    ! Last line
    !
    write(99,'(a)',advance='no') '+-'
    do idofn = 1,pnode*kdime
       do idime = 1,lreal+1
          write(99,'(a)',advance='no') ' '
       end do
    end do

    write(99,'(a)',advance='no') ' '
    do idofn = 1,pnode
       do idime = 1,lreal+1
          write(99,'(a)',advance='no') ' '
       end do
    end do


    write(99,'(a)',advance='no') '-+'
    write(99,'(a)',advance='no') '  +-   -+'

    write(99,'(a)',advance='no') '     +-'
    do idime = 1,lreal
       write(99,'(a)',advance='no') ' '
    end do
    write(99,'(a)',advance='no') '-+'

    write(99,*)

  end subroutine output_element_system

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2019-05-02
  !> @brief   Output domain
  !> @details Output the mesh
  !>
  !-----------------------------------------------------------------------

  subroutine output_domain(CURRENT_MESH,ONLY_MESH)

    use def_master
    use def_kermod
    use def_domain
    use def_mpio
    use mod_postpr
    use mod_mpio_config, only : mpio_config
    use mod_iofile,      only : iofile_flush_unit
    use mod_iofile,      only : iofile_open_unit
    use mod_iofile,      only : iofile_close_unit
    use mod_memory
    implicit none

    integer(ip), optional, intent(in) :: CURRENT_MESH
    logical(lg), optional, intent(in) :: ONLY_MESH
    integer(ip)                       :: ipart,ifiel,kdime,istep,idivi
    logical(lg)                       :: output_boundaries
    real(rp),    pointer              :: dumm2(:,:),dumm1(:)
    logical(lg)                       :: if_only_mesh
    !
    ! Options
    !
    if( present(ONLY_MESH) ) then
       if_only_mesh = ONLY_MESH
    else
       if_only_mesh = .false.
    end if
    if( present(CURRENT_MESH) ) then
       idivi = CURRENT_MESH
    else
       idivi = kfl_posdi
    end if

    nullify(dumm2,dumm1)

    if( maxval(kfl_oumes) == 1 .and. output_check_mesh_postprocess() ) then

       !-----------------------------------------------------------------
       !
       ! Output mesh information
       !
       !-----------------------------------------------------------------

       call messages_live('OUTPUT MESH','START SECTION')

       if( .not. if_only_mesh .and. INOTSLAVE ) then
          call iofile_open_unit(lun_pos00,fil_pos00,'POSTPROCESS INFO')
          if( IMASTER ) then
             if( mpio_config%output%merge ) then
                !
                ! Merged parallel
                !
                write(lun_pos00,1) 1_ip
                write(lun_pos00,2) 1_ip,&
                     &             meshe(idivi) % nelem_total,&
                     &             meshe(idivi) % npoin_origi,&
                     &             meshe(idivi) % nboun_total
             else
                !
                ! Distributed parallel
                !
                write(lun_pos00,1) npart
                do ipart = 1,npart
                   write(lun_pos00,2) ipart,&
                        &             meshe(idivi) % nelem_par(ipart),&
                        &             meshe(idivi) % npoin_par(ipart),&
                        &             meshe(idivi) % nboun_par(ipart)
                end do
             end if
          else if( ISEQUEN ) then
             !
             ! Sequential
             !
             write(lun_pos00,1) 1_ip
             write(lun_pos00,2) 1_ip,&
                  &          meshe(idivi) % nelem,&
                  &          meshe(idivi) % npoin,&
                  &          meshe(idivi) % nboun
          end if
          call iofile_flush_unit(lun_pos00)
          call iofile_close_unit(lun_pos00,fil_pos00,'POSTPROCESS INFO')
       end if

       !-----------------------------------------------------------------
       !
       ! Output basic element, boundary and node arrays
       !
       !-----------------------------------------------------------------
       
       if( kfl_oumes(1) == 1  .or. mpio_config%output%post_process%export_only ) then

          call postpr(meshe(idivi) % coord,    postp(1) % wopos(:,16),ittim,cutim,ndime,'ORIGI')    ! COORD
          call postpr(meshe(idivi) % lnods,    postp(1) % wopos(:,15),ittim,cutim,mnode,'ORIGI')    ! LNODS
          call postpr(meshe(idivi) % ltype,    postp(1) % wopos(:,17),ittim,cutim,      'ORIGI')    ! LTYPE
          call postpr(meshe(idivi) % lninv_loc,postp(1) % wopos(:,18),ittim,cutim,      'ORIGI')    ! LNINV_LOC
          call postpr(meshe(idivi) % leinv_loc,postp(1) % wopos(:,23),ittim,cutim,      'ORIGI')    ! LEINV_LOC

          if ( .not. mpio_config%output%post_process%light ) then
            call postpr(meshe(idivi) % lelch,    postp(1) % wopos(:,19),ittim,cutim,      'ORIGI')  ! LELCH
            call postpr(meshe(idivi) % lnoch,    postp(1) % wopos(:,27),ittim,cutim,      'ORIGI')  ! LNOCH
            call postpr(meshe(idivi) % lesub,    postp(1) % wopos(:,22),ittim,cutim,      'ORIGI')  ! LESUB
            call postpr(meshe(idivi) % lmate,    postp(1) % wopos(:,24),ittim,cutim,      'ORIGI')  ! LMATE
            call postpr(meshe(idivi) % lmast,    postp(1) % wopos(:,75),ittim,cutim,      'ORIGI')  ! LMAST
          end if

          output_boundaries = .false.
          if( ( nboun_total > 0  .or.  (mpio_config%output%post_process%enabled .and. .not. mpio_config%output%post_process%light) ) ) output_boundaries = .true.

          if( output_boundaries ) then
             call postpr(meshe(idivi) % lnodb,    postp(1) % wopos(:,20),ittim,cutim,mnodb,'ORIGI') ! LNODB
             call postpr(meshe(idivi) % ltypb,    postp(1) % wopos(:,21),ittim,cutim,      'ORIGI') ! LTYPB
             call postpr(meshe(idivi) % lboch,    postp(1) % wopos(:,26),ittim,cutim,      'ORIGI') ! LBOCH
             call postpr(meshe(idivi) % lelbo,    postp(1) % wopos(:,33),ittim,cutim,      'ORIGI') ! LELBO
             call postpr(meshe(idivi) % lbinv_loc,postp(1) % wopos(:,25),ittim,cutim,      'ORIGI') ! LBINV_LOC
          end if
       end if

       !-----------------------------------------------------------------
       !
       ! Output SETS field
       !
       !-----------------------------------------------------------------

       if(  mpio_config%output%post_process%export_only .or. kfl_oumes(2) == 1 .and. .not. mpio_config%output%post_process%light ) then
          if( nnset_origi > 0 ) call postpr(lnset,postp(1) % wopos(:,34),ittim,cutim,'ORIGI')                          ! LNSET
          if( neset_origi > 0 ) call postpr(leset,postp(1) % wopos(:,30),ittim,cutim,'ORIGI')                          ! LESET
          if( nbset_origi > 0 .and. output_boundaries ) call postpr(lbset,postp(1) % wopos(:,31),ittim,cutim,'ORIGI')  ! LBSET
       end if

       !-----------------------------------------------------------------
       !
       ! Output BOUNDARY_CONDITIONS
       !
       !-----------------------------------------------------------------

       if( mpio_config%output%post_process%export_only .or.  kfl_oumes(3) == 1 .and. .not. mpio_config%output%post_process%light ) then
          if( kfl_icodb > 0 .and. output_boundaries ) call postpr(kfl_codbo,postp(1) % wopos(:,28),ittim,cutim,      'ORIGI') ! CODBO
          if( kfl_icodn > 0 )                         call postpr(kfl_codno,postp(1) % wopos(:,29),ittim,cutim,mcono,'ORIGI') ! CODNO
       end if

       !-----------------------------------------------------------------
       !
       ! Output FIELDS
       !
       !-----------------------------------------------------------------

       if( mpio_config%output%post_process%export_only .or. kfl_oumes(4) == 1 .and. .not. mpio_config%output%post_process%light ) then

          do ifiel = 1,nfiel

             if ( (kfl_field(6,ifiel) /= 1 ) .OR. mpio_config%output%post_process%export_only ) then !not on demand or export

                if(      kfl_field(2,ifiel) == NELEM_TYPE ) then
                   postp(1) % wopos(3,32) = 'NELEM'
                else if( kfl_field(2,ifiel) == NPOIN_TYPE ) then
                   postp(1) % wopos(3,32) = 'NPOIN'
                else if( kfl_field(2,ifiel) == NBOUN_TYPE ) then
                   postp(1) % wopos(3,32) = 'NBOUN'
                else
                   call runend('OUTDOM: UNDEFINED FIELD TYPE')
                end if
                kdime = kfl_field(1,ifiel)


                do istep = 1,kfl_field(4,ifiel)
                   if( kdime == 1 ) then
                      postp(1) % wopos(2,32) = 'SCALA'
                      if( associated(xfiel(ifiel) % a) ) then
                         dumm1 => xfiel(ifiel) % a(1,:,istep)
                         call postpr(dumm1,postp(1) % wopos(:,32),ittim,time_field(ifiel) % a(istep),'ORIGI',TAG1=ifiel,TAG2=istep)          ! XFIEL
                      else
                         call postpr(dumm1,postp(1) % wopos(:,32),ittim,cutim,'ORIGI',TAG1=ifiel,TAG2=istep)
                      end if
                   else
                      postp(1) % wopos(2,32) = 'VECTO'
                      if( associated(xfiel(ifiel) % a) ) then
                         dumm2 => xfiel(ifiel) % a(:,:,istep)
                         call postpr(dumm2,postp(1) % wopos(:,32),ittim,time_field(ifiel) % a(istep),kdime,'ORIGI',TAG1=ifiel,TAG2=istep)    ! XFIEL
                      else
                         call postpr(dumm2,postp(1) % wopos(:,32),ittim,cutim,kdime,'ORIGI',TAG1=ifiel,TAG2=istep)
                      end if
                   end if
                end do
             end if !not on demand
          end do
       end if

       call messages_live('OUTPUT MESH','END SECTION')

    end if

    !-----------------------------------------------------------------
    !
    ! Output boundary mesh in STL format
    !
    !-----------------------------------------------------------------

    if( .not. if_only_mesh .and. kfl_oustl /= 0 ) call output_stl()

1   format(i9)
2   format(10(1x,i9))

  end subroutine output_domain

  !-----------------------------------------------------------------
  !
  ! Change file name in case of repartitioning
  !
  !-----------------------------------------------------------------

  subroutine output_file_names(REPARTITIONING,AMR)

    logical(lg), optional, intent(in) :: REPARTITIONING
    logical(lg), optional, intent(in) :: AMR
    logical(lg)                       :: if_repartitioning
    logical(lg)                       :: if_amr
    integer(ip)                       :: ii,jj

    if( present(REPARTITIONING) ) then
       if_repartitioning = REPARTITIONING
    else
       if_repartitioning = .false.
    end if
    if( present(AMR) ) then
       if_amr = AMR
    else
       if_amr = .false.
    end if

    if( kfl_repart_par /= 0 .and. INOTSLAVE .and. if_repartitioning ) then
       if( kfl_repart_post_par == 1 ) then
          ii        = index(fil_pos00_save,'.post')
          jj        = len(trim(fil_pos00_save))
          fil_pos00 = trim(fil_pos00_save(1:ii-1))//'-repartition'//trim(intost(num_repart_par))//trim(fil_pos00_save(ii:jj))

          ii        = index(fil_pos01_save,'.post')
          jj        = len(trim(fil_pos01_save))
          fil_pos01 = trim(fil_pos01_save(1:ii-1))//'-repartition'//trim(intost(num_repart_par))//trim(fil_pos01_save(ii:jj))

          ii        = index(fil_pos02_save,'.post')
          jj        = len(trim(fil_pos02_save))
          fil_pos02 = trim(fil_pos02_save(1:ii-1))//'-repartition'//trim(intost(num_repart_par))//trim(fil_pos02_save(ii:jj))
       else if( kfl_repart_post_par == 2 ) then
          fil_pos00 = 'repartition'//trim(intost(num_repart_par))//'/'//trim(fil_pos00_save)
          fil_pos01 = 'repartition'//trim(intost(num_repart_par))//'/'//trim(fil_pos01_save)
          fil_pos02 = 'repartition'//trim(intost(num_repart_par))//'/'//trim(fil_pos02_save)
       end if
    end if

    if( kfl_amr /= 0 .and. INOTSLAVE .and. if_amr ) then
       if( kfl_amr_post == 1 ) then
          ii        = index(fil_pos00_save,'.post')
          jj        = len(trim(fil_pos00_save))
          fil_pos00 = trim(fil_pos00_save(1:ii-1))//'-amr'//trim(intost(num_amr))//trim(fil_pos00_save(ii:jj))

          ii        = index(fil_pos01_save,'.post')
          jj        = len(trim(fil_pos01_save))
          fil_pos01 = trim(fil_pos01_save(1:ii-1))//'-amr'//trim(intost(num_amr))//trim(fil_pos01_save(ii:jj))

          ii        = index(fil_pos02_save,'.post')
          jj        = len(trim(fil_pos02_save))
          fil_pos02 = trim(fil_pos02_save(1:ii-1))//'-amr'//trim(intost(num_amr))//trim(fil_pos02_save(ii:jj))
       else if( kfl_amr_post == 2 ) then
          fil_pos00 = 'amr'//trim(intost(num_amr))//'/'//trim(fil_pos00_save)
          fil_pos01 = 'amr'//trim(intost(num_amr))//'/'//trim(fil_pos01_save)
          fil_pos02 = 'amr'//trim(intost(num_amr))//'/'//trim(fil_pos02_save)
       end if
    end if

  end subroutine output_file_names

  subroutine output_open_files(REPARTITIONING,AMR)

    logical(lg), optional, intent(in) :: REPARTITIONING
    logical(lg), optional, intent(in) :: AMR
    logical(lg)                       :: if_repartitioning
    logical(lg)                       :: if_amr

    if( present(REPARTITIONING) ) then
       if_repartitioning = REPARTITIONING
    else
       if_repartitioning = .false.
    end if
    if( present(AMR) ) then
       if_amr = AMR
    else
       if_amr = .false.
    end if
    !
    ! Close some files. LUN_POS00 is opened and closed on the fly in mod_output
    !
    if( INOTSLAVE ) then
       call output_file_names(REPARTITIONING,AMR)
       if( if_repartitioning ) then
          if( kfl_repart_par /= 0 .and. kfl_repart_post_par /= 0 ) then
             if( kfl_repart_post_par == 2 ) then
                call iofile_create_directory('repartition'//integer_to_string(num_repart_par))
             end if
             call iofile_close_unit(lun_pos01)
             call iofile_close_unit(lun_pos02)
             call iofile_open_unit(lun_pos01,fil_pos01,'POSTPROCESS FILE NAMES')
             call iofile_open_unit(lun_pos02,fil_pos02,'POSTPROCESS LOG')
          end if
       end if
       if( if_amr ) then
          if( kfl_amr /= 0 .and. kfl_amr_post /= 0 ) then
             if( kfl_amr_post == 2 ) then
                call iofile_create_directory('amr'//integer_to_string(num_amr))
             end if
             call iofile_close_unit(lun_pos01)
             call iofile_close_unit(lun_pos02)
             call iofile_open_unit(lun_pos01,fil_pos01,'POSTPROCESS FILE NAMES')
             call iofile_open_unit(lun_pos02,fil_pos02,'POSTPROCESS LOG')
          end if
       end if
    end if

    call PAR_BARRIER()

  end subroutine output_open_files

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2019-09-19
  !> @brief   Check if mesh should be postprocess
  !> @details Check if mesh should be postprocessed
  !>
  !-----------------------------------------------------------------------

  logical(lg) function output_check_mesh_postprocess()

    use def_parall, only : kfl_repart_par
    use def_parall, only : kfl_repart_post_par

    integer(ip), save :: ipass = 0

    ipass = ipass + 1
    !
    ! Reaprtitioning: postprocess mesh all the time if files are tagged
    !
    if( ipass == 1 ) then
       output_check_mesh_postprocess = .true.
    else
       output_check_mesh_postprocess = .false.
       if( kfl_repart_par /= 0 .and. kfl_repart_post_par > 0 ) output_check_mesh_postprocess =.true.
       if( kfl_amr        /= 0 .and. kfl_amr_post        > 0 ) output_check_mesh_postprocess =.true.
       if( kfl_amr        /= 0                               ) output_check_mesh_postprocess =.true.
    end if

  end function output_check_mesh_postprocess

  subroutine output_domain_alya_format(meshe,lunit)

    use mod_elmgeo, only : element_type
    use def_domain, only : nelty,lquad,ngaus,kfl_extra
    use def_master, only : NELEM_TYPE
    use def_master, only : NBOUN_TYPE
    use def_master, only : NPOIN_TYPE

    type(mesh_type), intent(in) :: meshe
    integer(ip),     intent(in) :: lunit
    integer(4)                  :: lunit4
    integer(ip)                 :: lexis_loc(nelty)
    integer(ip)                 :: lexib_loc(nelty)
    integer(ip)                 :: ielty,nelty_loc,ielem,kelty
    integer(ip)                 :: ifiel,iboun,nblty_loc,ipoin
    integer(ip)                 :: ienti,nenti

    lunit4 = int(lunit,4)
    !
    ! Dimensions
    !
    lexis_loc = 0
    do ielem = 1,meshe % nelem
       ielty = meshe % ltype(ielem)
       lexis_loc(ielty) =  1
    end do
    lexib_loc = 0
    do iboun = 1,meshe % nboun
       ielty = meshe % ltypb(iboun)
       lexib_loc(ielty) =  1
    end do
    nelty_loc = 0
    do ielty = 1,nelty
       if( lexis_loc(ielty) /= 0 ) then
          nelty_loc = nelty_loc + 1
          lexis_loc(nelty_loc) = ielty
       end if
    end do
    nblty_loc = 0
    do ielty = 1,nelty
       if( lexib_loc(ielty) /= 0 ) then
          nblty_loc = nblty_loc + 1
          lexib_loc(nblty_loc) = ielty
       end if
    end do

    write(lunit4,1)      &
         meshe % npoin,  &
         meshe % nelem,  &
         meshe % ndime,  &
         meshe % nboun
    write(lunit4,2) &
         lexis_loc(1:nelty_loc)
    if( meshe % nfiel > 0 ) then
       write(lunit4,3) meshe % nfiel
       do ifiel = 1,meshe % nfiel
          if(      meshe % kfl_field(2,ifiel) == NELEM_TYPE ) then
             write(lunit4,4) ifiel,meshe % kfl_field(1,ifiel),'ELEMENT'
          else if( meshe % kfl_field(2,ifiel) == NPOIN_TYPE ) then
             write(lunit4,4) ifiel,meshe % kfl_field(1,ifiel),'NODE'
          else if( meshe % kfl_field(2,ifiel) == NBOUN_TYPE ) then
             write(lunit4,4) ifiel,meshe % kfl_field(1,ifiel),'BOUNDARY'
          end if
       end do
       write(lunit4,5)
    end if
    write(lunit4,6)
    !
    ! Strategy
    !
    write(lunit4,10)
    if( maxval(lquad) == 1 ) then
       write(lunit4,11) 'CLOSE'
    else
       write(lunit4,11) 'OPEN'
    end if
    write(lunit4,12)
    do kelty = 1,nblty_loc
       ielty = lexib_loc(kelty)
       write(lunit4,13) element_type(ielty) % name,ngaus(ielty)
    end do
    do kelty = 1,nelty_loc
       ielty = lexis_loc(kelty)
       write(lunit4,13) element_type(ielty) % name,ngaus(ielty)
    end do
    write(lunit4,14)
    if( kfl_extra == 1 ) then
       write(lunit4,15) 'ON'
    else
       write(lunit4,15) 'OFF'
    end if
    write(lunit4,16)
    if( meshe % kfl_ngrou > 0 ) then
       write(lunit4,17) meshe % kfl_ngrou
    end if
    write(lunit4,19)
    !
    ! Geometry
    !
    write(lunit4,100)
    write(lunit4,101)
    do ielem = 1,meshe % nelem
       write(lunit4,*) ielem,meshe % ltype(ielem)
    end do
    write(lunit4,102)
    write(lunit4,103)
    do ielem = 1,meshe % nelem
       write(lunit4,*) ielem,meshe % lnods(1:meshe%lnnod(ielem),ielem)
    end do
    write(lunit4,104)
    if( associated(meshe % lmate) ) then
       write(lunit4,111)
       do ielem = 1,meshe % nelem
          write(lunit4,*) ielem,meshe % lmate(ielem)
       end do
       write(lunit4,112)
    end if
    if( associated(meshe % lesub) ) then
       write(lunit4,113)
       do ielem = 1,meshe % nelem
          write(lunit4,*) ielem,meshe % lesub(ielem)
       end do
       write(lunit4,114)
    end if

    write(lunit4,105)
    do ipoin = 1,meshe % npoin
       write(lunit4,'(i6,3(1x,e13.6))') ipoin,meshe % coord(1:meshe % ndime,ipoin)
    end do
    write(lunit4,106)
    write(lunit4,107)
    do iboun = 1,meshe % nboun
       write(lunit4,*) iboun,meshe % lelbo(iboun)
    end do
    write(lunit4,108)
    write(lunit4,109)
    do iboun = 1,meshe % nboun
       write(lunit4,*) iboun,meshe % lnodb(1:meshe % lnnob(iboun),iboun)
    end do
    write(lunit4,110)
    write(lunit4,199)
    !
    ! Sets
    !
    write(lunit4,200)
    if( associated(meshe % leset) ) then
       write(lunit4,203)
       do ielem = 1,meshe % nelem
          write(lunit4,*) ielem,meshe % leset(ielem)
       end do
       write(lunit4,204)
    end if
    if( associated(meshe % lbset) ) then
       write(lunit4,205)
       do iboun = 1,meshe % nboun
          write(lunit4,*) iboun,meshe % lbset(iboun)
       end do
       write(lunit4,206)
    end if
    if( associated(meshe % lnset) ) then
       write(lunit4,207)
       do ipoin = 1,meshe % npoin
          write(lunit4,*) ipoin,meshe % lnset(ipoin)
       end do
       write(lunit4,208)
    end if
    write(lunit4,299)
    !
    ! Boundary conditions
    !
    write(lunit4,300)
    if( associated(meshe % kfl_codbo) ) then
       write(lunit4,301)
       do iboun = 1,meshe % nboun
          write(lunit4,*) iboun,meshe % kfl_codbo(iboun)
       end do
       write(lunit4,302)
    end if
    write(lunit4,399)
    !
    ! Fields
    !
    write(lunit4,400)
    do ifiel = 1,meshe % nfiel
       if( associated(meshe % xfiel(ifiel) % a) ) then
          if(      meshe % kfl_field(2,ifiel) == NELEM_TYPE ) then
             nenti = meshe % nelem
          else if( meshe % kfl_field(2,ifiel) == NPOIN_TYPE ) then
             nenti = meshe % npoin
          else if( meshe % kfl_field(2,ifiel) == NBOUN_TYPE ) then
             nenti = meshe % nboun
          end if
          write(lunit4,401) ifiel
          do ienti = 1,nenti
             write(lunit4,'(i6,10(1x,e13.6))') ienti,meshe % xfiel(ifiel) % a(:,ienti,1)
          end do
          write(lunit4,402)
       end if
    end do
    write(lunit4,499)

1   format(&
         & '$------------------------------------------------------------',/,&
         & 'DIMENSIONS                                                   ',/,&
         & '  NODAL_POINTS=      ',i7,/,&
         & '  ELEMENTS=          ',i7,/,&
         & '  SPACE_DIMENSIONS=  ',i7,/,&
         & '  BOUNDARIES=        ',i7)
2   format(&
         & '  TYPES_OF_ELEMENTS= ',10(1x,i3))
3   format(&
         & '  FIELDS=            ',i7)
4   format(&
         & '    FIELD=',i4,' DIMENSION =',i4,' ,',a)
5   format(&
         & '  END_FIELDS')
6   format(&
         & 'END_DIMENSIONS                                               ',/,&
         & '$------------------------------------------------------------')

10  format(&
         & 'STRATEGY')
11  format(&
         & '  INTEGRATION_RULE:               ',a)
12  format(&
         & '  GAUSS_POINTS')
13  format(&
         & '    ',a5,1x,i3)
14  format(&
         & '  END_GAUSS_POINTS')
15  format(&
         & '  EXTRAPOLATE_BOUNDARY_CONDITIONS: ',a)
16  format(&
         & '  BOUNDARY_ELEMENT:                On')
17  format(&
         & '  GROUPS: FIELD= ',i5)
19 format(&
         & 'END_STRATEGY                                                 ',/,&
         & '$------------------------------------------------------------')
100 format(&
         & 'GEOMETRY')
101 format('  TYPES')
102 format('  END_TYPES')
103 format('  ELEMENTS')
104 format('  END_ELEMENTS')
105 format('  COORDINATES')
106 format('  END_COORDINATES')
107 format('  LELBO')
108 format('  END_LELBO')
109 format('  BOUNDARIES')
110 format('  END_BOUNDARIES')
111 format('  MATERIALS')
112 format('  END_MATERIALS')
113 format('  SUBDOMAINS')
114 format('  END_SUBDOMAINS')
199 format(&
         & 'END_GEOMETRY                                                 ',/,&
         & '$------------------------------------------------------------')
200 format(&
         & 'SETS')
203 format(&
         & '  ELEMENTS')
204 format(&
         & '  END_ELEMENTS')
205 format(&
         & '  BOUNDARIES')
206 format(&
         & '  END_BOUNDARIES')
207 format(&
         & '  NODES')
208 format(&
         & '  END_NODES')
299 format(&
         & 'END_SETS                                                     ',/,&
         & '$------------------------------------------------------------')
300 format(&
         & 'BOUNDARY_CONDITIONS')
301 format(&
         & '  ON_BOUNDARIES')
302 format(&
         & '  END_ON_BOUNDARIES')
399 format(&
         & 'END_BOUNDARY_CONDITIONS                                      ',/,&
         & '$------------------------------------------------------------')
400 format(&
         & 'FIELDS')
401 format(&
         & '  FIELD: ',i5)
402 format(&
         & '  END_FIELD')
499 format(&
         & 'END_FIELDS                                      ',/,&
         & '$------------------------------------------------------------')

  end subroutine output_domain_alya_format

end module mod_output
!> @}
