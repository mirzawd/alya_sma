!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    mod_par_output_partition.f90
!> @author  Guillaume Houzeaux
!> @date    20/04/2014
!> @brief   Output statistics
!> @details Output some information about the partition, the mesh,
!>          the hybrid parallelization and the coupling.
!>          This subroutine is called twice: at the beginning of the
!>          run and at the end to output timings
!>          The files can be read with GiD.
!>
!>          To visualize any distribution with gnuplot. For example,
!>          the histogram of number of elements;
!>
!>          1. Copy NUM_ELEMENTS field (2nd column) in a separate file caca.txt
!>          2. Execute the following in gnuplot (change binwidth if necessary).
!>          
!>          set xlabel 'Number of elements'
!>          set ylabel 'Number of subdomains' 
!>          binwidth=5
!>          bin(x,width)=width*floor(x/width)
!>          set style fill solid border -1
!>          plot 'caca.txt' using (bin($1,binwidth)):(1.0) not smooth freq with boxes
!----------------------------------------------------------------------

module mod_par_output_partition

  use def_kintyp
  use def_parame
  use mod_parall
  use def_domain
  use def_master
  use def_coupli
  use def_kermod,                    only : kfl_edge_elements
  use def_parall,                    only : kfl_matri_par
  use def_parall,                    only : kfl_connectivity_par
  use def_parall,                    only : fil_conne_par
  use mod_parall,                    only : num_subd_par
  use mod_parall,                    only : num_pack_par
  use mod_parall,                    only : list_elements_par
  use mod_elmgeo,                    only : element_type
  use mod_maths_sort,                only : maths_heap_sort
  use mod_iofile,                    only : iofile
  use mod_iofile,                    only : iofile_flush_unit
  use mod_iofile,                    only : iofile_open_unit
  use mod_iofile,                    only : iofile_close_unit
  use mod_iofile,                    only : iofile_append_tag
  use mod_memory,                    only : memory_alloca
  use mod_memory,                    only : memory_deallo
  use mod_memory,                    only : memory_resize 
  use mod_memory,                    only : memory_size 
  use mod_outfor,                    only : outfor
  use mod_messages,                  only : livinf
  use mod_messages,                  only : messages_live
  use mod_par_parallel_partitioning, only : par_partition_total_weight
  use mod_par_affinity,              only : par_affinity
  use mod_matrix_market,             only : matrix_market_matrix
  use mod_matrix_market,             only : matrix_market_vector
  use def_parall,                    only : kfl_output_partition
  use def_parall,                    only : kfl_output_node_comm_arrays
  use def_parall,                    only : kfl_output_edge_comm_arrays
  use mod_strings,                   only : integer_to_string
  use def_AMR,                       only : kfl_amr
  use def_AMR,                       only : num_amr
  use def_search_method,             only : search_method
  use def_maths_bin,                 only : maths_bin
  use def_maths_tree,                only : maths_octree
  use def_maths_tree,                only : maths_kdtree
  use def_mpi
#include "def_mpi.inc"
  use mod_std
  use mod_communications
  implicit none  
  private

  integer(8) :: output_memor(2)

  interface par_output_value
     module procedure par_output_value_I,&
          &           par_output_value_R
  end interface par_output_value

  interface par_output_coupling_timings
     module procedure par_output_coupling_timings_one,&
          &           par_output_coupling_timings_all
  end interface par_output_coupling_timings

  public :: par_output_partition
  public :: par_output_maphys
  public :: par_output_global_matrix
  public :: par_output_solvers
  public :: par_output_coupling_timings
  public :: par_output_search_timings
  
contains

  !-----------------------------------------------------------------------
  !
  !> @date    19/12/2016
  !> @author  Guillaume Houzeaux
  !> @brief   Output of partitions and timings
  !> @details Output statistics on the parallelization
  !
  !-----------------------------------------------------------------------

  subroutine par_output_partition(itask)

    integer(ip), intent(in)       :: itask                  !< Task to do
    integer(ip)                   :: ielem,ipoin,imodu
    integer(ip)                   :: ipart,idime,dummi
    integer(ip)                   :: dummm,imate,nmate_total
    integer(ip)                   :: icode,ipart_world
    integer(ip)                   :: kelem,ii,isubd,ipack
    integer(ip)                   :: number_wet_points
    integer(ip)                   :: npoin_wet,ineig,jpart       
    integer(ip)                   :: nboun_wet,jelem,ielty     
    integer(ip)                   :: nneig_source,icoup    
    integer(ip)                   :: iposi(1),jpart_last
    type(comm_data_par), pointer  :: commu
    MY_MPI_COMM                   :: PAR_COMM_TO_USE
    real(rp)                      :: dummr,xaver,xx
    real(rp)                      :: xmini(3)
    real(rp)                      :: xmaxi(3)
    real(rp)                      :: xmini_send(3)
    real(rp)                      :: xmaxi_send(3)
    real(rp)                      :: xcoor_send(3)
    real(rp)                      :: cpu_ave_assembly(mmodu)
    real(rp)                      :: cpu_ave_boundary(mmodu)
    integer(ip)                   :: nneig
    integer(ip), pointer          :: lneig(:)
    integer(ip), pointer          :: sneig(:)
    real(rp),    pointer          :: xcoor(:)
    real(rp),    pointer          :: xcoor_gat(:)
    integer(ip), pointer          :: nneig_gat(:)
    integer(4),  pointer          :: nneig4_gat(:)
    integer(ip), pointer          :: lneig_gat(:)
    integer(ip), pointer          :: sneig_gat(:)
    integer(ip), pointer          :: lresu_gat(:)
    real(rp),    pointer          :: rresu_gat(:)
    integer(ip), pointer          :: lcode_gat(:)
    type(i1p),   pointer          :: ledge(:)
    type(i1p),   pointer          :: sedge(:)
    integer(ip), pointer          :: ja(:)

    integer(ip)                   :: lexis_gat(nelty)
    integer(ip)                   :: element_stat(21)
    integer(ip)                   :: boundary_stat(21)
    integer(ip)                   :: neighbor_stat(21)

    integer(ip), pointer          :: kfl_modul_tmp(:)

    integer(ip)                   :: dom_i
    
    integer(ip)                   :: par_omp_num_threads_gat
    integer(ip)                   :: par_hybrid_ompss_gat
    integer(ip)                   :: par_hybrid_openmp_gat
    integer(ip)                   :: num_cores_max
    integer(ip)                   :: host_num,num_hosts
    integer(ip),      pointer     :: core_num(:)
    character(len=:), pointer     :: fil_conne_coupling

    if( ISEQUEN ) return

    nullify( commu      )
    nullify( lneig      )
    nullify( sneig      )
    nullify( xcoor      )
    nullify( xcoor_gat  )
    nullify( nneig_gat  )
    nullify( nneig4_gat )
    nullify( lneig_gat  )
    nullify( sneig_gat  )
    nullify( lresu_gat  )
    nullify( rresu_gat  )
    nullify( lcode_gat  )
    nullify( ledge      )
    nullify( sedge      )
    nullify( core_num   )
    nullify( fil_conne_coupling )
    nullify( ja )

    select case ( itask )

    case ( 1_ip ) 
       !
       ! Node communication arrays
       !
       if( kfl_output_node_comm_arrays /= 0 .and. INOTMASTER ) then
          allocate(ja(commd % bound_dim))
          do ineig = 1,commd % nneig
             dom_i = commd % neights(ineig)
             do ii = commd % bound_size(ineig),commd % bound_size(ineig+1)-1
                ja(ii) = dom_i
             end do
          end do
          call matrix_market_matrix(1_ip,1_ip,commd % nneig,IA=commd % bound_size,JA=ja,A=commd % bound_perm,&
               FILENAME=trim(namda)//'-bound_perm_glo-'//integer_to_string(kfl_paral)//'.mtx',PERMA=lninv_loc)
          call matrix_market_matrix(1_ip,1_ip,commd % nneig,IA=commd % bound_size,JA=ja,A=commd % bound_perm,&
               FILENAME=trim(namda)//'-bound_perm_loc-'//integer_to_string(kfl_paral)//'.mtx')

          call matrix_market_matrix(1_ip,1_ip,commd % nneig,IA=commd % bound_size,JA=ja,A=commd % bound_invp,&
               FILENAME=trim(namda)//'-bound_invp_glo-'//integer_to_string(kfl_paral)//'.mtx',PERMA=lninv_loc)
          call matrix_market_matrix(1_ip,1_ip,commd % nneig,IA=commd % bound_size,JA=ja,A=commd % bound_invp,&
               FILENAME=trim(namda)//'-bound_invp_loc-'//integer_to_string(kfl_paral)//'.mtx')

          call matrix_market_vector(commd % bound_size,&
               FILENAME=trim(namda)//'-bound_size-'//integer_to_string(kfl_paral)//'.mtx')
          call matrix_market_vector(commd % neights,&
               FILENAME=trim(namda)//'-neights-'//integer_to_string(kfl_paral)//'.mtx')
          deallocate(ja)          
       end if
       !
       ! Edge communication arrays
       !
       if( kfl_output_edge_comm_arrays /= 0 .and. kfl_edge_elements == 1 .and. INOTMASTER ) then
          allocate(ja(commd % bedge_dim))
          do ineig = 1,commd % nneig
             dom_i = commd % neights(ineig)
             do ii = commd % bedge_size(ineig),commd % bedge_size(ineig+1)-1
                ja(ii) = dom_i
             end do
          end do
          call matrix_market_matrix(1_ip,1_ip,commd % nneig,IA=commd % bedge_size,JA=ja,A=commd % bedge_perm,&
               FILENAME=trim(namda)//'-bedge_perm_glo-'//integer_to_string(kfl_paral)//'.mtx',PERMA=lginv_loc)
          call matrix_market_matrix(1_ip,1_ip,commd % nneig,IA=commd % bedge_size,JA=ja,A=commd % bedge_perm,&
               FILENAME=trim(namda)//'-bedge_perm_loc-'//integer_to_string(kfl_paral)//'.mtx')
          call matrix_market_vector(commd % bedge_size,&
               FILENAME=trim(namda)//'-bedge_size-'//integer_to_string(kfl_paral)//'.mtx')
          call matrix_market_vector(commd % neights,&
               FILENAME=trim(namda)//'-neights-'//integer_to_string(kfl_paral)//'.mtx')
           deallocate(ja)          
       end if
       !
       ! Partition info in gid format
       !
       if( kfl_output_partition /= 0 ) then

          if( PAR_MY_WORLD_RANK == 0 ) then
             call par_output_partition_open(MESH=.true.,RESULTS=.true.)
             !if( kfl_amr /= 0 ) then
             !   ii = index(fil_parti_msh,'.post')
             !   jj = len(trim(fil_parti_msh))
             !   call iofile_open_unit(lun_parti_msh,fil_parti_msh(1:ii-1)//'-amr'//integer_to_string(num_amr)//trim(fil_parti_msh(ii:jj),'PARALL PARTITION MESH')                
             !   ii = index(fil_parti_res,'.post')
             !   jj = len(trim(fil_parti_res))
             !   call iofile_open_unit(lun_parti_res,fil_parti_res(1:ii-1)//'-amr'//integer_to_string(num_amr)//trim(fil_parti_res(ii:jj),'PARALL PARTITION RESULT')
             !else
             !   call iofile_open_unit(lun_parti_msh,fil_parti_msh,'PARALL PARTITION MESH')
             !   call iofile_open_unit(lun_parti_res,fil_parti_res,'PARALL PARTITION RESULT')
             !end if 
          end if

          allocate( xcoor(3) )
          call livinf(0_ip,'PARALL: OUTPUT PARTITION',0_ip)

          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)

          element_stat  =  0 
          boundary_stat =  0
          neighbor_stat =  0
          xmini         =  huge(1.0_rp)*0.1_rp
          xmaxi         = -huge(1.0_rp)*0.1_rp
          xmini_send    =  huge(1.0_rp)*0.1_rp
          xmaxi_send    = -huge(1.0_rp)*0.1_rp
          xcoor         =  0.0_rp
          nneig         =  0

          if( PAR_MY_WORLD_RANK == 0 ) then
             allocate( xcoor_gat (  3*PAR_WORLD_SIZE  ) )
             allocate( nneig_gat (  0:PAR_WORLD_SIZE-1) )
             allocate( nneig4_gat(  0:PAR_WORLD_SIZE-1) )
             allocate( lresu_gat (  0:PAR_WORLD_SIZE-1) )
             allocate( rresu_gat (  0:PAR_WORLD_SIZE-1) )
             allocate( lcode_gat (  0:PAR_WORLD_SIZE-1) )
             allocate( ledge     (  0:PAR_WORLD_SIZE-1) )
             allocate( sedge     (  0:PAR_WORLD_SIZE-1) )
             do ipart = 0,PAR_WORLD_SIZE-1 
                nullify(ledge(ipart) % l)
                nullify(sedge(ipart) % l)
             end do
          end if
          !
          ! Partition coordinates
          !
          if( INOTMASTER ) then
             do ipoin = 1,npoin
                xcoor(1:ndime) = xcoor(1:ndime) + coord(1:ndime,ipoin)
                xmini(1:ndime) = min(xmini(1:ndime),coord(1:ndime,ipoin))
                xmaxi(1:ndime) = max(xmaxi(1:ndime),coord(1:ndime,ipoin))
             end do
             if( npoin > 0 ) xcoor(1:ndime) = xcoor(1:ndime) / real(npoin,rp)
             nneig = commu % nneig 
          end if
          !
          ! Master's coordinates
          ! 
          do icode = 1,mcode
             if( current_code == icode ) then
                xcoor_send = xcoor
             else       
                xcoor_send = 0.0_rp
             end if
             call PAR_SUM(3_ip,xcoor_send,'IN THE WORLD')
             if( IMASTER .and. icode == current_code ) then
                xcoor = xcoor_send / real(PAR_CODE_SIZE,rp)
             end if
          end do
          call PAR_GATHER(xcoor,xcoor_gat,'IN THE WORLD')
          !
          ! Number and list of neighbors
          !
          call PAR_GATHER(nneig,nneig_gat,'IN THE WORLD')

          if( PAR_MY_WORLD_RANK == 0 ) then
             dummi = 0
             do ipart = 0,PAR_WORLD_SIZE-1 
                dummi = dummi + nneig_gat(ipart)
                nneig4_gat(ipart) = int(nneig_gat(ipart),4)
             end do
             allocate( lneig_gat(dummi) )
             allocate( sneig_gat(dummi) )
          else
             if( INOTMASTER ) then
                if( commu % nneig > 0 ) allocate( lneig(commu % nneig) )
                if( commu % nneig > 0 ) allocate( sneig(commu % nneig) )
                do ineig = 1,commu % nneig   
                   ipart        = commd % neights(ineig)
                   ipart_world  = par_world_rank_of_a_code_neighbor(ipart,current_code)
                   lneig(ineig) = ipart_world
                   sneig(ineig) = commu % bound_size(ineig+1)-commu % bound_size(ineig)
                end do
             end if
          end if
          call PAR_GATHERV(lneig,lneig_gat,nneig4_gat,'IN THE WORLD')
          call PAR_GATHERV(sneig,sneig_gat,nneig4_gat,'IN THE WORLD')
          !
          ! Code number
          !
          call PAR_GATHER(current_code,lcode_gat,'IN THE WORLD') 

          !----------------------------------------------------------------------
          !
          ! Geometry
          !
          !----------------------------------------------------------------------

          if( PAR_MY_WORLD_RANK == 0 ) then
             !
             ! Partition graph LEDGE
             !
             ielem = 0
             do ipart = 0,PAR_WORLD_SIZE-1
                if( nneig_gat(ipart) > 0 ) then
                   allocate( ledge(ipart) % l(nneig_gat(ipart)) )
                   allocate( sedge(ipart) % l(nneig_gat(ipart)) )
                   do ineig = 1,nneig_gat(ipart)
                      ielem = ielem + 1
                      ledge(ipart) % l(ineig) = lneig_gat(ielem)
                      sedge(ipart) % l(ineig) = sneig_gat(ielem)
                   end do
                end if
             end do
             !
             ! Coordinates (Master are also present!)
             !
             ipoin = 0
             icode = 1
             write(lun_parti_msh,'(a,i1,a)') 'MESH CODE_'//trim(intost(icode))//' dimension ',ndime,' Elemtype Linear Nnode 2'     
             write(lun_parti_msh,'(a)') 'coordinates'  
             do ipart = 0,PAR_WORLD_SIZE-1 
                ipoin = ipoin + 1
                write(lun_parti_msh,'(i5,3(1x,e12.6))') ipoin,(xcoor_gat(ipart*3+idime),idime=1,ndime)
             end do
             write(lun_parti_msh,'(a)') 'end coordinates'
             !
             ! Elements: edges of the communication
             !
             ielem = 0
             do icode = 1,mcode
                if( icode > 1 ) &
                     write(lun_parti_msh,'(a,i1,a)') 'MESH CODE_'//trim(intost(icode))//' dimension ',ndime,' Elemtype Linear Nnode 2'
                write(lun_parti_msh,'(a)') 'elements'
                do ipart = 0,PAR_WORLD_SIZE-1 
                   if( lcode_gat(ipart) == icode ) then
                      do ineig = 1,nneig_gat(ipart)
                         ielem = ielem + 1
                         write(lun_parti_msh,'(3(1x,i7))') ielem,ipart+1,ledge(ipart) % l(ineig)+1
                      end do
                   end if
                end do
                write(lun_parti_msh,'(a)') 'end elements'
             end do
             !
             ! Connectivity file for gnuplot
             !
             ! set title  "Interface size"
             ! set xlabel "MPI rank"
             ! set ylabel "MPI rank"
             ! set xrange
             ! set yrange reverse
             ! plot '*-connectivity.par.res' matrix with image not
             !
             if( kfl_connectivity_par /= 0 ) then
                call iofile_open_unit(lun_conne_par,fil_conne_par,'MPI CONNECTIVITY FILE')             
                do ipart = 0,PAR_WORLD_SIZE-1
                   if( nneig_gat(ipart) > 0 ) then
                      !call maths_heap_sort(2_ip,nneig_gat(ipart),ledge(ipart) % l)
                      do jpart = 0,PAR_WORLD_SIZE-1
                         iposi = maxloc(ledge(ipart) % l(:),ledge(ipart) % l(:)==jpart)
                         if( iposi(1) <= 0 ) then
                            write(lun_conne_par,fmt='(i6,1x)',advance='no') 0
                         else
                            write(lun_conne_par,fmt='(i6,1x)',advance='no') sedge(ipart)%l(iposi(1))
                         end if
                      end do
                   else
                      do jpart = 0,PAR_WORLD_SIZE-1
                         write(lun_conne_par,fmt='(i6,1x)',advance='no') 0
                      end do
                   end if
                   write(lun_conne_par,*)
                end do
                call iofile_close_unit(lun_conne_par)
             end if
             !
             ! Deallocate
             !
             do ipart = 0,PAR_WORLD_SIZE-1
                if( nneig_gat(ipart) > 0 ) then
                   deallocate( ledge(ipart) % l )
                   deallocate( sedge(ipart) % l )
                end if
             end do
          end if
          !
          ! Code bounding boxes
          !
          ipoin = PAR_WORLD_SIZE
          do icode = 1,mcode
             if( current_code == icode ) then
                xmini_send = xmini
                xmaxi_send = xmaxi
             else       
                xmini_send =  huge(1.0_rp)
                xmaxi_send = -huge(1.0_rp)
             end if
             call PAR_MIN(3_ip,xmini_send,'IN THE WORLD')
             call PAR_MAX(3_ip,xmaxi_send,'IN THE WORLD')
             if( PAR_MY_WORLD_RANK == 0 ) then
                ielem = ielem + 1
                if( ndime == 2 ) then
                   write(lun_parti_msh,'(a)') 'MESH CODE'//trim(intost(icode))//'_BOUNDING_BOX dimension 2 Elemtype Quadrilateral Nnode 4'       
                   write(lun_parti_msh,'(a)') 'coordinates'  
                   ipoin = ipoin + 1 ; write(lun_parti_msh,'(i5,3(1x,e12.6))') ipoin,xmini_send(1),xmini_send(2)
                   ipoin = ipoin + 1 ; write(lun_parti_msh,'(i5,3(1x,e12.6))') ipoin,xmaxi_send(1),xmini_send(2)
                   ipoin = ipoin + 1 ; write(lun_parti_msh,'(i5,3(1x,e12.6))') ipoin,xmaxi_send(1),xmaxi_send(2)
                   ipoin = ipoin + 1 ; write(lun_parti_msh,'(i5,3(1x,e12.6))') ipoin,xmini_send(1),xmaxi_send(2)
                   write(lun_parti_msh,'(a)') 'end coordinates' 
                   write(lun_parti_msh,'(a)') 'elements'  
                   write(lun_parti_msh,'(5(1x,i7))') ielem,ipoin-3,ipoin-2,ipoin-1,ipoin
                   write(lun_parti_msh,'(a)') 'end elements'  
                else
                   write(lun_parti_msh,'(a)') 'MESH CODE'//trim(intost(icode))//'_BOUNDING_BOX dimension 3 Elemtype Hexahedra Nnode 8'       
                   write(lun_parti_msh,'(a)') 'coordinates'  
                   ipoin = ipoin + 1 ; write(lun_parti_msh,'(i5,3(1x,e12.6))') ipoin,xmini_send(1),xmini_send(2),xmini_send(3)
                   ipoin = ipoin + 1 ; write(lun_parti_msh,'(i5,3(1x,e12.6))') ipoin,xmaxi_send(1),xmini_send(2),xmini_send(3)
                   ipoin = ipoin + 1 ; write(lun_parti_msh,'(i5,3(1x,e12.6))') ipoin,xmaxi_send(1),xmaxi_send(2),xmini_send(3)
                   ipoin = ipoin + 1 ; write(lun_parti_msh,'(i5,3(1x,e12.6))') ipoin,xmini_send(1),xmaxi_send(2),xmini_send(3)
                   ipoin = ipoin + 1 ; write(lun_parti_msh,'(i5,3(1x,e12.6))') ipoin,xmini_send(1),xmini_send(2),xmaxi_send(3)
                   ipoin = ipoin + 1 ; write(lun_parti_msh,'(i5,3(1x,e12.6))') ipoin,xmaxi_send(1),xmini_send(2),xmaxi_send(3)
                   ipoin = ipoin + 1 ; write(lun_parti_msh,'(i5,3(1x,e12.6))') ipoin,xmaxi_send(1),xmaxi_send(2),xmaxi_send(3)
                   ipoin = ipoin + 1 ; write(lun_parti_msh,'(i5,3(1x,e12.6))') ipoin,xmini_send(1),xmaxi_send(2),xmaxi_send(3)
                   write(lun_parti_msh,'(a)') 'end coordinates' 
                   write(lun_parti_msh,'(a)') 'elements'  
                   write(lun_parti_msh,'(9(1x,i7))') ielem,ipoin-3,ipoin-2,ipoin-1,ipoin-0,&
                        &                               ipoin-7,ipoin-6,ipoin-5,ipoin-4
                   write(lun_parti_msh,'(a)') 'end elements'  
                end if
             end if
          end do

          !----------------------------------------------------------------------
          !
          ! Results: partition
          !
          !---------------------------------------------------------------------- 

          if( PAR_MY_WORLD_RANK == 0 ) then
             write(lun_parti_res,'(a)') 'GiD Post Results File 1.0'
             write(lun_parti_res,'(a)') ' '
             call iofile_flush_unit(lun_parti_res)
          end if
          !
          ! Number of neighbors
          !
          call par_output_value(0_ip,commu % nneig,'NUM_NEIGHBORS','PARTITION')
          !
          ! Compute neighbors stats
          ! 
          if( PAR_MY_WORLD_RANK == 0 ) then
             !
             ! xxxx--+--+--+--xxxx--+--+--+--xxxx--+--+--+--xxxx--+--+--+--xxxx
             !  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21
             ! xxxx--+--+--+--xxxx--+--+--+--xxxx--+--+--+--xxxx--+--+--+--xxxx
             ! av/2          av/1.5           av           1.5xav          2xav
             !
             xaver = 0.0_rp 
             do ipart = 0,PAR_WORLD_SIZE-1  
                xaver = xaver + real(nneig_gat(ipart),rp)
             end do
             xaver = xaver / real(PAR_WORLD_SIZE-mcode,rp) ! Remove master's contributions
             do ipart = 0,PAR_WORLD_SIZE-1  
                if( nneig_gat(ipart) > 0 ) then
                   xx    = real(nneig_gat(ipart),rp)
                   dummr = max(xx/xaver,xaver/xx)
                   if( xx >= xaver ) then
                      ii =   1 + int(10.0_rp*dummr,ip)
                   else
                      ii =  21 - int(10.0_rp*dummr,ip) 
                   end if
                   ii = min(max( 1_ip,ii),21_ip)
                   neighbor_stat(ii) = neighbor_stat(ii) + 1
                end if
             end do
          end if
          !
          ! My rank in my code
          !
          call par_output_value(0_ip,kfl_paral,'CODE_RANK','PARTITION','DO NOT AVERAGE')
          !
          ! Number of interior nodes per subdomain
          !     
          call par_output_value(0_ip,npoi1,'NUM_INTERIOR_NODES','PARTITION')
          !
          ! Number of boundary nodes per subdomain
          !
          call par_output_value(0_ip,npoin-npoi2+1_ip,'NUM_BOUNDARY_NODES','PARTITION','AVERAGE',lresu_gat)
          !
          ! Number of own nodes per subdomain
          !     
          call par_output_value(0_ip,npoin_own,'NUM_OWN_NODES','PARTITION')
          !
          ! Number of computational halo nodes per subdomain
          !     
          call par_output_value(0_ip,npoin_halo-npoin_own,'NUM_COMPUTATIONAL_HALO_NODES','PARTITION')
          !
          ! Compute boundary stats
          ! 
          if( PAR_MY_WORLD_RANK == 0 ) then
             !
             ! xxxx--+--+--+--xxxx--+--+--+--xxxx--+--+--+--xxxx--+--+--+--xxxx
             !  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21
             ! xxxx--+--+--+--xxxx--+--+--+--xxxx--+--+--+--xxxx--+--+--+--xxxx
             ! av/2          av/1.5           av           1.5xav          2xav
             !
             xaver = 0.0_rp 
             do ipart = 0,PAR_WORLD_SIZE-1  
                xaver = xaver + real(lresu_gat(ipart),rp)
             end do
             xaver = xaver / real(PAR_WORLD_SIZE-mcode,rp) ! Remove master's contributions
             do ipart = 0,PAR_WORLD_SIZE-1  
                if( lresu_gat(ipart) > 0 ) then
                   xx    = real(lresu_gat(ipart),rp)
                   dummr = max(xx/xaver,xaver/xx)
                   if( xx >= xaver ) then
                      ii =   1 + int(10.0_rp*dummr,ip)
                   else
                      ii =  21 - int(10.0_rp*dummr,ip) 
                   end if
                   ii = min(max( 1_ip,ii),21_ip)
                   boundary_stat(ii) = boundary_stat(ii) + 1
                end if
             end do
          end if
          !
          ! Number of elements per subdomain
          !
          call par_output_value(0_ip,nelem,'NUM_ELEMENTS','PARTITION')
          !
          ! Number of weighted elements per subdomain
          !
          !dummi = 0
          !do jelem = 1,nelem
          !   dummi = dummi + ngaus(abs(ltype(jelem)))
          !end do
          call par_partition_total_weight(dummi)
          call par_output_value(0_ip,dummi,'WEIGHTED_ELEMENTS','PARTITION','AVERAGE',lresu_gat)
          !
          ! Compute element stats
          ! 
          if( PAR_MY_WORLD_RANK == 0 ) then
             !
             ! xxxx--+--+--+--xxxx--+--+--+--xxxx--+--+--+--xxxx--+--+--+--xxxx
             !  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21
             ! xxxx--+--+--+--xxxx--+--+--+--xxxx--+--+--+--xxxx--+--+--+--xxxx
             ! av/2          av/1.5           av           1.5xav          2xav
             !
             xaver = 0.0_rp 
             do ipart = 0,PAR_WORLD_SIZE-1  
                xaver = xaver + real(lresu_gat(ipart),rp)
             end do
             xaver = xaver / real(PAR_WORLD_SIZE-mcode,rp) ! Remove master's contributions
             do ipart = 0,PAR_WORLD_SIZE-1  
                if( lresu_gat(ipart) > 0 ) then
                   xx    = real(lresu_gat(ipart),rp)
                   dummr = max(xx/xaver,xaver/xx)
                   if( xx >= xaver ) then
                      ii =   1 + int(10.0_rp*dummr,ip)
                   else
                      ii =  21 - int(10.0_rp*dummr,ip) 
                   end if
                   ii = min(max( 1_ip,ii),21_ip)
                   element_stat(ii) = element_stat(ii) + 1
                end if
             end do
          end if
          !
          ! Size of conenctivity
          !
          if( PAR_MY_WORLD_RANK == 0 ) then        
             write(lun_parti_res,'(a)')  'GaussPoints GP Elemtype Linear'
             write(lun_parti_res,'(a)')  'Number of Gauss Points: 1'
             write(lun_parti_res,'(a)')  'Natural Coordinates: Internal'
             write(lun_parti_res,'(a)')  'End GaussPoints'
             write(lun_parti_res,'(a)')  'Result INTERFACE_SIZE PARTITION 0 Scalar OnGaussPoints GP'
             write(lun_parti_res,'(a)')  'ComponentNames INTERFACE_SIZE'
             write(lun_parti_res,'(a)')  'Values'
             kelem = 0
             do ipart = 0,PAR_WORLD_SIZE-1
                do ineig = 1,nneig_gat(ipart)
                   kelem = kelem + 1
                   write(lun_parti_res,'(3(1x,i7))') kelem,sneig_gat(kelem)
                end do
             end do
             write(lun_parti_res,'(a)')  'End Values'  
          end if

          !----------------------------------------------------------------------
          !
          ! Host information... core_num is not postprocessed for the moment
          !
          !----------------------------------------------------------------------

          call memory_alloca(par_memor,'CORE_NUM','mod_par_output_partition',core_num,max(1_ip,par_omp_num_threads))
          call par_affinity(core_num,host_num,num_hosts)
          call par_output_value(0_ip,host_num,'HOST_NUMBER','HOST_INFORMATION')

          do ii = 1,memory_size(core_num)
             dummi = core_num(ii)
             if( dummi > 0 ) then
                if( count(core_num==dummi) > 1 ) &
                     call messages_live('MOD_PAR_OUTPUT_PARTITION: CHECK THE BINDING','WARNING')
             end if
          end do
          call par_output_value(0_ip,core_num(1),'CORE_NUMBER_LOCAL' ,'HOST_INFORMATION')
          num_cores_max = core_num(1)
          call PAR_MAX(num_cores_max)
          dummi = (host_num-1)*num_cores_max + core_num(1)
          call par_output_value(0_ip,dummi   ,'CORE_NUMBER_GLOBAL','HOST_INFORMATION')

          call memory_deallo(par_memor,'CORE_NUM','mod_par_output_partition',core_num)

          !----------------------------------------------------------------------
          !
          ! Hybrid parallelization
          ! The following treatment works for co-execution where the different
          ! Alya have a different hybrid strategy
          !
          !----------------------------------------------------------------------

          par_omp_num_threads_gat = par_omp_num_threads
          if( par_hybrid == PAR_OMPSS ) then
             par_hybrid_ompss_gat  = 1
             par_hybrid_openmp_gat = 0
          else
             par_hybrid_ompss_gat  = 0
             par_hybrid_openmp_gat = 1
          end if
          call PAR_MAX(par_omp_num_threads_gat,'IN THE WORLD')
          call PAR_MAX(par_hybrid_ompss_gat   ,'IN THE WORLD')
          call PAR_MAX(par_hybrid_openmp_gat  ,'IN THE WORLD')
          dummi = 0

          if( par_omp_num_threads_gat > 0 ) then

             if( par_hybrid_ompss_gat == 1 ) then
                !
                ! Number of OMPSS subdomains
                !
                if( par_hybrid == PAR_OMPSS ) then
                   if( INOTMASTER ) then
                      if( associated(ompss_domains) ) then
                         dummi = size(ompss_domains,KIND=ip)
                      else
                         dummi = 0
                      end if
                   end if
                else
                   dummi = 0
                end if
                call par_output_value(0_ip,dummi,'OMPSS_SUBDOMAINS','HYBRID_PARALLELIZATION')

                if( par_hybrid == PAR_OMPSS ) then
                   if( INOTMASTER ) then
                      if( associated(ompss_boundaries) ) then
                         dummi = size(ompss_boundaries,KIND=ip)
                      else
                         dummi = 0
                      end if
                   end if
                else
                   dummi = 0
                end if
                call par_output_value(0_ip,dummi,'OMPSS_BOUNDARIES','HYBRID_PARALLELIZATION')
             end if
             if( par_hybrid_openmp_gat == 1 ) then
                !
                ! Number of OpenMP colors
                !
                if( par_hybrid /= PAR_OMPSS ) then             
                   if( par_omp_num_threads > 0 ) then
                      dummi = par_omp_num_colors
                   else
                      dummi = 0
                   end if
                else
                   dummi = 0
                end if
                call par_output_value(0_ip,dummi,'OMP_COLOR','HYBRID_PARALLELIZATION')
             end if
             !
             ! Element OMP chunk size
             !
             if( par_omp_num_threads > 0 ) then
                dummi = par_omp_nelem_chunk
             else
                dummi = 0
             end if
             call par_output_value(0_ip,dummi,'OMP_CHUNK_ELEMENT','HYBRID_PARALLELIZATION')
             !
             ! Node OMP chunk size
             !
             if( par_omp_num_threads > 0 ) then
                dummi = par_omp_npoin_chunk
             else
                dummi = 0
             end if
             call par_output_value(0_ip,par_omp_npoin_chunk,'OMP_CHUNK_NODE','HYBRID_PARALLELIZATION')

          end if

          !----------------------------------------------------------------------
          !
          ! Co-execution
          !
          !----------------------------------------------------------------------

          dummi = 0
#ifdef OPENACCHHH
          dummi = 1
#endif
          dummm = dummi
          call PAR_MAX(dummm)
          if( dummm > 0 ) then
             call par_output_value(0_ip,dummi,'CPU_0_GPU_1','CO_EXECUTION')
             dummi = 0
             ii    = 0
             do isubd = 1,num_subd_par
                do ipack = 1,num_pack_par(isubd)
                   dummi = dummi + memory_size(list_elements_par(isubd) % packs(ipack) % l)
                   ii    = ii + 1
                end do
             end do
             dummr = real(dummi,rp)/real(max(ii,1_ip),rp) 
             call par_output_value(0_ip,dummr,'AVERAGE_VECTOR_SIZE','CO_EXECUTION')
          end if

          !----------------------------------------------------------------------
          !
          ! Mesh
          !
          !----------------------------------------------------------------------
          !
          ! Number elements per type
          !
          lexis_gat(1:nelty) = lexis(1:nelty)
          call PAR_SUM(nelty,lexis_gat,'IN THE WORLD') 
          do ielty = 1_ip,iesto_dom
             if( lexis_gat(ielty) > 0 ) then
                dummi = 0
                if( INOTMASTER ) then
                   if( ielty >= iesta_dom ) then
                      if( nelem > 0 ) dummi = count(ltype(1:nelem) == ielty,KIND=ip)
                   else
                      if( nboun > 0 ) dummi = count(ltypb(1:nboun) == ielty,KIND=ip)
                   end if
                end if
                call par_output_value(0_ip,dummi,'NUM_'//trim(element_type(ielty) % name),'MESH')
             end if
          end do
          !
          ! Non-zeros in the graph
          !
          call par_output_value(0_ip,nzdom,'NUMBER_NODE_GRAPH_ENTRIES','MESH')

          !----------------------------------------------------------------------
          !
          ! Materials
          !
          !----------------------------------------------------------------------

          nmate_total = nmate
          call PAR_MAX(nmate_total,'IN THE WORLD') 
          if( nmate_total > 1 ) then
             do imate = 1,nmate_total
                dummi = 0
                do ielem = 1,nelem
                   if( lmate(ielem) == imate ) dummi = dummi + 1
                end do
                call par_output_value(0_ip,dummi,'MAT_NUM_'//trim(intost(imate)),'MATERIALS')          
             end do
          end if

          !----------------------------------------------------------------------
          !
          ! Couplings
          !
          !----------------------------------------------------------------------

          if( mcoup > 0 ) then
             ! 
             ! My code number
             !
             call par_output_value(0_ip,current_code,'CODE_NUMBER','PARTITION','DO NOT AVERAGE')
             !
             ! My world rank 
             !
             call par_output_value(0_ip,PAR_MY_WORLD_RANK,'WORLD_RANK','PARTITION','DO NOT AVERAGE')
          end if

          do icoup = 1,mcoup
             color_target      = coupling_type(icoup) % color_target
             color_source      = coupling_type(icoup) % color_source
             number_wet_points = coupling_type(icoup) % wet   % number_wet_points
             npoin_wet         = coupling_type(icoup) % wet   % npoin_wet
             nboun_wet         = coupling_type(icoup) % wet   % nboun_wet
             nneig_source      = coupling_type(icoup) % commd % nneig
             !
             ! Geometry
             !
             call PAR_GATHER(nneig_source,nneig_gat,'IN THE WORLD') 

             if( PAR_MY_WORLD_RANK == 0 ) then        
                dummi = 0
                do ipart = 0,PAR_WORLD_SIZE-1 
                   dummi = dummi + nneig_gat(ipart)
                   nneig4_gat(ipart) = int(nneig_gat(ipart),4)
                end do
                if( associated(lneig_gat) ) deallocate( lneig_gat )
                allocate( lneig_gat(dummi) )
                if( associated(sneig_gat) ) deallocate( sneig_gat )
                allocate( sneig_gat(dummi) )
             else
                if( INOTMASTER ) then
                   if( associated(lneig) ) deallocate( lneig )
                   allocate( lneig(coupling_type(icoup) % commd % nneig) )
                   if( associated(sneig) ) deallocate( sneig )
                   allocate( sneig(coupling_type(icoup) % commd % nneig) )
                   do ineig = 1,coupling_type(icoup) % commd % nneig   
                      ipart        = coupling_type(icoup) % commd % neights(ineig)
                      ipart_world  = par_color_coupling_rank_to_world(ipart,color_target,color_source)
                      lneig(ineig) = ipart_world
                   end do
                end if
             end if

             call PAR_GATHERV(lneig,lneig_gat,nneig4_gat,'IN THE WORLD')

             if( PAR_MY_WORLD_RANK == 0 ) then        
                write(lun_parti_msh,'(a,i1,a)') 'MESH COUPLING'//trim(intost(icoup))//' dimension ',ndime,' Elemtype Linear Nnode 2'
                write(lun_parti_msh,'(a)') 'elements'
                kelem = 0
                jelem = ielem
                do ipart = 0,PAR_WORLD_SIZE-1
                   do ineig = 1,nneig_gat(ipart)
                      kelem = kelem + 1
                      ielem = ielem + 1
                      write(lun_parti_msh,'(3(1x,i7))') ielem,ipart+1,lneig_gat(kelem)+1
                   end do
                end do
                write(lun_parti_msh,'(a)') 'end elements' 
             end if
             !
             ! Weight
             !
             if( INOTMASTER ) then
                do ineig = 1,coupling_type(icoup) % commd % nneig   
                   sneig(ineig) =   coupling_type(icoup) % commd % lrecv_size(ineig+1) &
                        &         - coupling_type(icoup) % commd % lrecv_size(ineig) 
                end do
             end if

             call PAR_GATHERV(sneig,sneig_gat,nneig4_gat,'IN THE WORLD')

             if( PAR_MY_WORLD_RANK == 0 ) then
                write(lun_parti_res,'(a)')  'GaussPoints GP Elemtype Linear'
                write(lun_parti_res,'(a)')  'Number of Gauss Points: 1'
                write(lun_parti_res,'(a)')  'Natural Coordinates: Internal'
                write(lun_parti_res,'(a)')  'End GaussPoints'
                write(lun_parti_res,'(a)')  'Result COUPLING_'//trim(intost(icoup))//'_SIZE COUPLING 0 Scalar OnGaussPoints GP'
                write(lun_parti_res,'(a)')  'ComponentNames COUPLING_'//trim(intost(icoup))//'_SIZE'
                write(lun_parti_res,'(a)')  'Values'
                kelem = 0
                do ipart = 0,PAR_WORLD_SIZE-1
                   do ineig = 1,nneig_gat(ipart)
                      kelem = kelem + 1
                      jelem = jelem + 1
                      write(lun_parti_res,'(3(1x,i7))') jelem,sneig_gat(kelem)
                   end do
                end do
                write(lun_parti_res,'(a)')  'End Values'
                !
                ! Connectivity file for gnuplot:
                !
                ! set title  "Interface size"
                ! set xlabel "MPI rank"
                ! set ylabel "MPI rank"
                ! set xrange
                ! set yrange reverse
                ! plot '*-connectivity-*.par.res' matrix with image not
                !
                if( kfl_connectivity_par /= 0 ) then
                   call memory_alloca(output_memor,'FIL_CONNE_PAR','membcs',fil_conne_coupling,len(fil_conne_par,kind=ip))
                   call iofile_append_tag(fil_conne_par,icoup,fil_conne_coupling)
                   call iofile_open_unit(lun_conne_par,fil_conne_coupling,'MPI CONNECTIVITY FILE FOR COUPLING')
                   kelem = 0
                   do ipart = 0,PAR_WORLD_SIZE-1
                      if( nneig_gat(ipart) == 0 ) then
                         do jpart = 0,PAR_WORLD_SIZE-1
                            write(lun_conne_par,fmt='(i6,1x)',advance='no') 0
                         end do
                      else
                         jpart_last = -1
                         call maths_heap_sort(2_ip,nneig_gat(ipart),lneig_gat(kelem+1:))
                         do ineig = 1,nneig_gat(ipart)
                            kelem = kelem + 1
                            do jpart = jpart_last+1,lneig_gat(kelem)-1
                               write(lun_conne_par,fmt='(i6,1x)',advance='no') 0
                            end do
                            write(lun_conne_par,fmt='(i6,1x)',advance='no') sneig_gat(kelem)
                            jpart_last = lneig_gat(kelem)
                         end do
                         do jpart = jpart_last+1,PAR_WORLD_SIZE-1
                            write(lun_conne_par,fmt='(i6,1x)',advance='no') 0
                         end do
                      end if
                      write(lun_conne_par,*)
                   end do
                   call iofile_close_unit(lun_conne_par)
                   call memory_deallo(output_memor,'FIL_CONNE_PAR','membcs',fil_conne_coupling)
                end if

             end if
             !
             ! Results
             !
             call par_output_value(0_ip,number_wet_points,                       'COUPLI'//trim(intost(icoup))//'_WET_POINTS'     ,'COUPLING')
             call par_output_value(0_ip,npoin_wet,                               'COUPLI'//trim(intost(icoup))//'_WET_NODES'      ,'COUPLING')
             call par_output_value(0_ip,nboun_wet,                               'COUPLI'//trim(intost(icoup))//'_WET_BOUNDARIES' ,'COUPLING')
             call par_output_value(0_ip,nneig_source,                            'COUPLI'//trim(intost(icoup))//'_NEIGHBORS'      ,'COUPLING')
             call par_output_value(0_ip,coupling_type(icoup) % commd % lrecv_dim,'COUPLI'//trim(intost(icoup))//'_RECEIVE_SIZE'   ,'COUPLING')
             call par_output_value(0_ip,coupling_type(icoup) % commd % lsend_dim,'COUPLI'//trim(intost(icoup))//'_SEND_SIZE'      ,'COUPLING')

          end do
          !
          ! Deallocate
          !
          if( associated(lneig)      ) deallocate( lneig      )
          if( associated(sneig)      ) deallocate( sneig      )
          if( associated(xcoor)      ) deallocate( xcoor      )
          if( associated(xcoor_gat)  ) deallocate( xcoor_gat  )
          if( associated(nneig_gat)  ) deallocate( nneig_gat  )
          if( associated(nneig4_gat) ) deallocate( nneig4_gat )
          if( associated(lneig_gat)  ) deallocate( lneig_gat  )
          if( associated(sneig_gat)  ) deallocate( sneig_gat  )
          if( associated(lresu_gat)  ) deallocate( lresu_gat  )
          if( associated(rresu_gat)  ) deallocate( rresu_gat  )
          if( associated(lcode_gat)  ) deallocate( lcode_gat  )
          if( associated(ledge)      ) deallocate( ledge      )
          if( associated(sedge)      ) deallocate( sedge      )

          if( PAR_MY_WORLD_RANK == 0 ) then
             call iofile_close_unit(lun_parti_msh)
             call iofile_close_unit(lun_parti_res)
             !close(lun_parti_res)     
          end if

          !----------------------------------------------------------------------
          !
          ! Output statistics: See using gnuplot as follows:
          !
          ! reset
          ! set xtics nomirror; set ytics nomirror; set border front;
          ! set xrange [-1:21]; set yrange [0:100]
          ! set xlabel 'Dispersion'; set ylabel 'Percentage'
          ! set style data histogram; set style histogram cluster gap 0; set style fill solid
          ! set xtics ("<0.5*average" 1, "average" 11, ">2*average" 21)
          ! set ytics ("0" 0, "20" 20, "40" 40, "60" 60, "80" 80, "100" 100)
          ! plot 'caca.txt' u 1 t 'Elements' lc rgb "red", '' u 2 t 'Boundaries' lc rgb "green", '' u 3 t 'Neighbors' lc rgb "blue"
          !
          !----------------------------------------------------------------------

          call PAR_BROADCAST(21_ip,element_stat, 'IN THE WORLD')
          call PAR_BROADCAST(21_ip,boundary_stat,'IN THE WORLD')
          call PAR_BROADCAST(21_ip,neighbor_stat,'IN THE WORLD')
          xaver = real(PAR_WORLD_SIZE-mcode,rp) 
          call outfor(70_ip,lun_outpu_par,' ')
          do ii = 1,21
             routp(1) = real(element_stat(ii), rp) / xaver * 100.0_rp
             routp(2) = real(boundary_stat(ii),rp) / xaver * 100.0_rp
             routp(3) = real(neighbor_stat(ii),rp) / xaver * 100.0_rp
             call outfor(71_ip,lun_outpu_par,' ')
          end do

       end if

    case ( 2_ip )

       if( kfl_output_partition /= 0 ) then

          if( PAR_MY_WORLD_RANK == 0 ) call par_output_partition_open(RESULTS=.true.,APPEND=.true.)
          
          !----------------------------------------------------------------------
          !
          ! Timings
          !
          !----------------------------------------------------------------------

          !if( INOTSLAVE ) & 
          !     call iofile(zero,lun_parti_res,fil_parti_res,'PARALL PARTITION RESULT','old','formatted','append')
          if( PAR_MY_WORLD_RANK == 0 ) then
             allocate( lresu_gat (  0:PAR_WORLD_SIZE-1) )
             allocate( rresu_gat (  0:PAR_WORLD_SIZE-1) )
          end if
          !
          ! Agree on modules that are solver
          !
          nullify( kfl_modul_tmp )
          allocate( kfl_modul_tmp(lbound(kfl_modul,1):ubound(kfl_modul,1)) )
          kfl_modul_tmp = abs(kfl_modul)
          call PAR_MAX(kfl_modul_tmp,'IN THE WORLD')

          do imodu = 1,mmodu-1

             if( kfl_modul_tmp(imodu) /= 0 ) then
                !
                ! my_assembly_time 
                !
                cpu_ave_assembly(imodu) = cpu_modul(CPU_ASSEMBLY,imodu)
                call par_output_value(0_ip,cpu_ave_assembly(imodu),'ELEMENT_ASSEMBLY_TIME_'//trim(namod(imodu)),'TIMING')
                call PAR_AVERAGE(cpu_ave_assembly(imodu),'IN MY CODE')
                !
                ! my_assembly_time / average_assembly_time
                !
                dummr = cpu_modul(CPU_ASSEMBLY,imodu) / (cpu_ave_assembly(imodu)+zeror)
                call par_output_value(0_ip,dummr,'ELEMENT_ASSEMBLY_RATIO_'//trim(namod(imodu)),'TIMING')
                !
                ! my_boundary_time 
                !
                cpu_ave_boundary(imodu) = cpu_modul(CPU_ASSEMBLY_BOUNDARY,imodu)
                call par_output_value(0_ip,cpu_ave_boundary(imodu),'BOUNDARY_ASSEMBLY_TIME_'//trim(namod(imodu)),'TIMING')
                call PAR_AVERAGE(cpu_ave_boundary(imodu),'IN MY CODE')
                !
                ! my_boundary_time / average_boundary_time
                !
                dummr = cpu_modul(CPU_ASSEMBLY_BOUNDARY,imodu) / (cpu_ave_boundary(imodu)+zeror)
                call par_output_value(0_ip,dummr,'BOUNDARY_ASSEMBLY_RATIO_'//trim(namod(imodu)),'TIMING')
             end if
          end do

          !close(lun_parti_res)     
          deallocate(kfl_modul_tmp)

          if( associated(lresu_gat)  ) deallocate( lresu_gat  )
          if( associated(rresu_gat)  ) deallocate( rresu_gat  )
          !
          ! Coupling
          !
          if( mcoup > 0 ) call par_output_coupling_timings()

          if( PAR_MY_WORLD_RANK == 0 ) call iofile_close_unit(lun_parti_res)
          
       end if

    end select

  end subroutine par_output_partition

  !-----------------------------------------------------------------------
  !>
  !> @date    08/01/2018
  !> @author  Guillaume Houzeaux
  !> @brief   Output solvers statistics
  !> @details Output solvers statistics
  !>          1. Chunk size when OpenMP is used:
  !>             0 ... OPENMP OFF
  !>            -2 ... STATIC SCHEDULING
  !>            -3 ... GUIDED SCHEDULING
  !>            >0 ... CHUNK SIZE FOR DYNAMIC SCHEDULING
  !>          2. Matrix size: number of coefficients
  !>          3. Number of groups
  !>
  !-----------------------------------------------------------------------

  subroutine par_output_solvers()

    use def_master, only : solve_sol
    use def_master, only : mmodu,momod,modul,kfl_modul

    integer(ip)                        :: ivari
    integer(ip)                        :: tmp_size_solve(mmodu)
    integer(ip)                        :: tmp_ivari_solve

    if( IPARALL ) then

       if( PAR_MY_WORLD_RANK == 0 ) call par_output_partition_open(RESULTS=.true.,APPEND=.true.)
       !
       ! Check which modules are solved
       !
       tmp_size_solve = 0
       do modul = 1,mmodu 
          if( kfl_modul(modul) == 1 ) then 
             if( associated(momod(modul) % solve) ) then
                tmp_size_solve(modul) = size(momod(modul) % solve)
             end if
          end if
       end do
       call PAR_MAX(mmodu,tmp_size_solve,'IN THE WORLD')

       do modul = 1,mmodu 
          !
          ! This module is solved, check which solver should be postprocessed
          !
          do ivari = 1,tmp_size_solve(modul)
             nullify(solve_sol)
             tmp_ivari_solve = 0
             if( kfl_modul(modul) == 1 ) then 
                if( associated(momod(modul) % solve) ) then
                   if(    momod(modul) % solve(ivari) % kfl_solve /= 0 .and. &
                        & momod(modul) % solve(ivari) % kfl_algso >  0 .and. &
                        & momod(modul) % solve(ivari) % kfl_algso /= 9 ) then
                      tmp_ivari_solve =  1
                      solve_sol       => momod(modul) % solve(ivari:)
                   end if
                end if
             end if
             call PAR_MAX(tmp_ivari_solve,'IN THE WORLD')
             if( tmp_ivari_solve == 1 ) then
                call par_output_single_solver(modul,ivari,solve_sol)
             end if
          end do
       end do

       if( PAR_MY_WORLD_RANK == 0 ) call iofile_close_unit(lun_parti_res)

    end if

  end subroutine par_output_solvers

  !-----------------------------------------------------------------------
  !
  !> @date    19/12/2016
  !> @author  Guillaume Houzeaux
  !> @brief   Output solvers statistics
  !> @details Output solvers statistics
  !
  !-----------------------------------------------------------------------

  subroutine par_output_single_solver(imodu,ivari,solve)

    use def_solver
    use mod_maths, only : maths_list_different_elements

    integer(ip),           intent(in) :: imodu
    integer(ip),           intent(in) :: ivari
    type(soltyp), pointer, intent(in) :: solve(:)
!    integer(ip)                       :: omp_chunk_size,my_ngrou
    integer(ip)                       :: my_nzmat
!    integer(ip)                       :: par_hybrid_max
    real(rp)                          :: spmv_total,spmv_ratio(2)
    real(rp)                          :: cpu_spmv(2) 
    !
    ! Matrix size and SpMV timings
    !
    if( associated(solve) ) then
       my_nzmat      = solve(1) % nzmat
       cpu_spmv(1:2) = solve(1) % cpu_spmv(1:2)
       spmv_total    = solve(1) % cpu_spmv(1) + solve(1) % cpu_spmv(2)
       spmv_ratio(1) = solve(1) % cpu_spmv(1)
       spmv_ratio(2) = spmv_total
       call PAR_AVERAGE(2_ip,spmv_ratio,'IN MY CODE')
       spmv_ratio(1) = solve(1) % cpu_spmv(1) / (spmv_ratio(1)+zeror)
       spmv_ratio(2) = spmv_total             / (spmv_ratio(2)+zeror)
    else
       my_nzmat      = 0
       cpu_spmv      = 0.0_rp
       spmv_total    = 0.0_rp
       spmv_ratio(1) = 0.0_rp
       spmv_ratio(2) = 0.0_rp
    end if
       
    call par_output_value(0_ip,my_nzmat     ,'MATRIX_NUM_COEFF_'      //namod(imodu)//'_PROB_'//integer_to_string(ivari),'SOLVERS') 
    call par_output_value(0_ip,cpu_spmv(1)  ,'SPMV_TIME_COMPUTATION_' //namod(imodu)//'_PROB_'//integer_to_string(ivari),'SOLVERS')
    call par_output_value(0_ip,spmv_total   ,'SPMV_TIME_TOTAL_'       //namod(imodu)//'_PROB_'//integer_to_string(ivari),'SOLVERS')
    call par_output_value(0_ip,spmv_ratio(1),'SPMV_RATIO_COMPUTATION_'//namod(imodu)//'_PROB_'//integer_to_string(ivari),'SOLVERS')
    call par_output_value(0_ip,spmv_ratio(2),'SPMV_RATIO_TOTAL_'      //namod(imodu)//'_PROB_'//integer_to_string(ivari),'SOLVERS')
    !
    ! Number of groups
    !
    !if(    solve % kfl_algso == SOL_SOLVER_DEFLATED_CG           .or. &
    !     & solve % kfl_algso == SOL_SOLVER_PIPELINED_DEFLATED_CG ) then
    !   if( INOTMASTER ) call maths_list_different_elements(solve % lgrou,my_ngrou,OPTION='STRICTLY POSITIVE')
    !   call par_output_value(0_ip,my_ngrou,'NUMBER_OF_GROUPS_'//trim(solve % wprob),'SOLVERS')
    !else
    !   call par_output_value(0_ip,0_ip,'NUMBER_OF_GROUPS_'//trim(solve % wprob),'SOLVERS')
    !end if
    !
    ! Chunk size
    !
    !par_hybrid_max = par_hybrid
    !call PAR_MAX(par_hybrid_max)
    !if( par_hybrid_max /= 0 ) then
    !   if(      solve % omp_schedule == SOL_OMP_OFF     ) then
    !      omp_chunk_size =  0
    !   else if( solve % omp_schedule == SOL_OMP_STATIC  ) then
    !      omp_chunk_size = -2
    !   else if( solve % omp_schedule == SOL_OMP_GUIDED  ) then
    !      omp_chunk_size = -3
    !   else if( solve % omp_schedule == SOL_OMP_DYNAMIC ) then           
    !      omp_chunk_size = solve % omp_chunk_size
    !   end if
    !   !call par_output_value(0_ip,omp_chunk_size,'OMP_CHUNK_SIZE_'//trim(solve % wprob),'SOLVERS',wherein='IN MY CODE')
    !end if    

  end subroutine par_output_single_solver

  !-----------------------------------------------------------------------
  !
  !> @date    19/12/2016
  !> @author  Guillaume Houzeaux
  !> @brief   Output MAPHYS statistics
  !> @details Output some MAPHYS statistics
  !
  !-----------------------------------------------------------------------

  subroutine par_output_maphys(&
       mem_maphys_facto, mem_maphys_schur, mem_maphys_precond ,mem_maphys_solver,&
       time_maphys_facto,time_maphys_analysis,time_maphys_precond,time_maphys_solver)

    integer(ip), intent(in) :: mem_maphys_facto
    integer(ip), intent(in) :: mem_maphys_schur
    integer(ip), intent(in) :: mem_maphys_precond
    integer(ip), intent(in) :: mem_maphys_solver
    real(rp),    intent(in) :: time_maphys_facto 
    real(rp),    intent(in) :: time_maphys_analysis
    real(rp),    intent(in) :: time_maphys_precond
    real(rp),    intent(in) :: time_maphys_solver

    integer(ip)             :: jtime

    !----------------------------------------------------------------------
    !
    ! Timing MAPHYS
    !
    !----------------------------------------------------------------------

    !if( INOTSLAVE ) & 
    !     call iofile(zero,lun_parti_res,fil_parti_res,'PARALL PARTITION RESULT','old','formatted','append')
    jtime = solve_sol(1) % nsolv + 1
    !
    ! Facto memory
    !
    call par_output_value(jtime,mem_maphys_facto    ,'MEM_FACTO_' //trim(solve_sol(1) % wprob),'MAPHYS')
    call par_output_value(jtime,mem_maphys_schur    ,'MEM_SCHUR_' //trim(solve_sol(1) % wprob),'MAPHYS')
    call par_output_value(jtime,mem_maphys_precond  ,'MEM_PRECO_' //trim(solve_sol(1) % wprob),'MAPHYS')
    call par_output_value(jtime,mem_maphys_solver   ,'MEM_SOLVE_' //trim(solve_sol(1) % wprob),'MAPHYS')

    call par_output_value(jtime,time_maphys_facto   ,'TIME_FACTO_'//trim(solve_sol(1) % wprob),'MAPHYS')
    call par_output_value(jtime,time_maphys_analysis,'TIME_ANALY_'//trim(solve_sol(1) % wprob),'MAPHYS')
    call par_output_value(jtime,time_maphys_precond ,'TIME_PRECO_'//trim(solve_sol(1) % wprob),'MAPHYS')
    call par_output_value(jtime,time_maphys_solver  ,'TIME_SOLVE_'//trim(solve_sol(1) % wprob),'MAPHYS')

    !close(lun_parti_res)     

  end subroutine par_output_maphys

  !-----------------------------------------------------------------------
  !
  !> @date    19/12/2016
  !> @author  Guillaume Houzeaux
  !> @brief   Output coupling timings
  !> @details Output coupling timings
  !
  !-----------------------------------------------------------------------

  subroutine par_output_coupling_timings_all()

    integer(ip)   :: jtime,icoup,jcoup
    character(20) :: wcoup

    !if( INOTSLAVE ) & 
    !     call iofile(zero,lun_parti_res,fil_parti_res,'PARALL PARTITION RESULT','old','formatted','append')
    jtime = 0_ip

    do icoup = 1,mcoup
       jcoup = coupling_type(icoup) % number
       wcoup = 'COUPLI'//trim(intost(icoup))//'_'
       call par_output_value(jtime,coupling_type(icoup) % cputim( 1),trim(wcoup)//'TIME_DEFINE_WET_NODES'      ,'COUPLING')
       call par_output_value(jtime,coupling_type(icoup) % cputim(10),trim(wcoup)//'TIME_HOLCUT'                ,'COUPLING')
       call par_output_value(jtime,coupling_type(icoup) % cputim( 2),trim(wcoup)//'TIME_COMPUTE_KD_TREES'      ,'COUPLING')
       call par_output_value(jtime,coupling_type(icoup) % cputim( 6),trim(wcoup)//'TIME_DISTRIBUTE_WET_POINTS' ,'COUPLING')
       call par_output_value(jtime,coupling_type(icoup) % cputim( 3),trim(wcoup)//'TIME_SEND_WET_COORD_TO_CPUS','COUPLING')
       call par_output_value(jtime,coupling_type(icoup) % cputim( 4),trim(wcoup)//'TIME_CPUS_CHECK_OWNING'     ,'COUPLING')
       call par_output_value(jtime,coupling_type(icoup) % cputim( 5),trim(wcoup)//'TIME_CREATE_COUPLING_COMM'  ,'COUPLING')
       call par_output_value(jtime,coupling_type(icoup) % cputim( 7),trim(wcoup)//'TIME_GENERATE_TRANSM_MATRIX','COUPLING')
       call par_output_value(jtime,coupling_type(icoup) % cputim( 8),trim(wcoup)//'TIME_SEND_RECV_COUPLING'    ,'COUPLING')
       call par_output_value(jtime,coupling_type(icoup) % cputim( 9),trim(wcoup)//'TIME_MULT_TRANSM_MATRIX'    ,'COUPLING')
    end do

    !close(lun_parti_res)     

  end subroutine par_output_coupling_timings_all

  subroutine par_output_coupling_timings_one(coupling,MY_TIME)

    type(typ_color_coupling), intent(inout)          :: coupling
    integer(ip),              intent(in),   optional :: MY_TIME
    integer(ip)                                      :: jtime
    integer(ip),              save                   :: jcoup
    character(20)                                    :: wcoup

    !if( INOTSLAVE ) & 
    !     call iofile(zero,lun_parti_res,fil_parti_res,'PARALL PARTITION RESULT','old','formatted','append')
    if( present(MY_TIME) ) then
       jtime = MY_TIME
    else
       jtime = 0_ip
    end if
    jcoup = coupling % number

    wcoup = 'COUPLI'//trim(intost(jcoup))//'_'
    call par_output_value(jtime,coupling % cputim( 1),trim(wcoup)//'TIME_DEFINE_WET_NODES'      ,'COUPLING')
    call par_output_value(jtime,coupling % cputim(10),trim(wcoup)//'TIME_HOLCUT'                ,'COUPLING')
    call par_output_value(jtime,coupling % cputim( 2),trim(wcoup)//'TIME_COMPUTE_KD_TREES'      ,'COUPLING')
    call par_output_value(jtime,coupling % cputim( 6),trim(wcoup)//'TIME_DISTRIBUTE_WET_POINTS' ,'COUPLING')
    call par_output_value(jtime,coupling % cputim( 3),trim(wcoup)//'TIME_SEND_WET_COORD_TO_CPUS','COUPLING')
    call par_output_value(jtime,coupling % cputim( 4),trim(wcoup)//'TIME_CPUS_CHECK_OWNING'     ,'COUPLING')
    call par_output_value(jtime,coupling % cputim( 5),trim(wcoup)//'TIME_CREATE_COUPLING_COMM'  ,'COUPLING')
    call par_output_value(jtime,coupling % cputim( 7),trim(wcoup)//'TIME_GENERATE_TRANSM_MATRIX','COUPLING')
    call par_output_value(jtime,coupling % cputim( 8),trim(wcoup)//'TIME_SEND_RECV_COUPLING'    ,'COUPLING')
    call par_output_value(jtime,coupling % cputim( 9),trim(wcoup)//'TIME_MULT_TRANSM_MATRIX'    ,'COUPLING')

    !close(lun_parti_res)     

  end subroutine par_output_coupling_timings_one

  subroutine par_output_search_timings(search,MY_TIME)

    class(search_method),     intent(inout)          :: search
    integer(ip),              intent(in),   optional :: MY_TIME
    integer(ip)                                      :: jtime
    character(20)                                    :: wcoup

    if( kfl_output_partition /= 0 .and. IPARALL ) then

       if( INOTSLAVE ) call par_output_partition_open(RESULTS=.true.,APPEND=.true.)

       if( present(MY_TIME) ) then
          jtime = MY_TIME
       else
          jtime = 0_ip
       end if
       
       wcoup = 'SEARCH_'//trim(search % name)//'_'
       call par_output_value(jtime,search % times(1),trim(wcoup)//'FILL'      ,'SEARCH',wherein='IN MY CODE')
       call par_output_value(jtime,search % times(2),trim(wcoup)//'CANDIDATE' ,'SEARCH',wherein='IN MY CODE')
       
       select type ( search )
       class is ( maths_bin ) 
          call par_output_value(jtime,search % stats(1),trim(wcoup)//'AVE_ELEMENTS','SEARCH',wherein='IN MY CODE')
          call par_output_value(jtime,search % stats(2),trim(wcoup)//'MAX_ELEMENTS','SEARCH',wherein='IN MY CODE')
       class is ( maths_octree ) 
          call par_output_value(jtime,search % stats(1),trim(wcoup)//'NUM_LEAVES'  ,'SEARCH',wherein='IN MY CODE')
          call par_output_value(jtime,search % stats(3),trim(wcoup)//'NUM_LEVELS'  ,'SEARCH',wherein='IN MY CODE')
          call par_output_value(jtime,search % stats(4),trim(wcoup)//'MAX_ELEMENTS','SEARCH',wherein='IN MY CODE')
       end select

       if( INOTSLAVE ) call iofile_close_unit(lun_parti_res)

    end if
    
  end subroutine par_output_search_timings

  !-----------------------------------------------------------------------
  !
  !> @date    19/12/2016
  !> @author  Guillaume Houzeaux
  !> @brief   Output values in Gid format
  !> @details Perform a gather of a value to be output by the master of 
  !>          the universe. The values of the code masters is the average
  !>          in the code
  !
  !-----------------------------------------------------------------------

  subroutine par_output_value_I(jtime,my_ivalue,message1,message2,message3,ivalue_gat,wherein)

    integer(ip),  intent(in)                    :: jtime
    integer(ip),  intent(in)                    :: my_ivalue
    character(*), intent(in)                    :: message1
    character(*), intent(in)                    :: message2
    character(*), intent(in), optional          :: message3
    integer(ip),  intent(in), optional, pointer :: ivalue_gat(:)
    character(*), intent(in), optional          :: wherein
    integer(ip),  pointer                       :: ivalue(:)
    integer(ip)                                 :: ave_ivalue,dummi,ipart
    logical(lg)                                 :: laverage
    logical(lg)                                 :: in_my_code

    in_my_code = .false.
    if( present(wherein) ) then
       if(trim(wherein) == 'IN MY CODE' ) in_my_code = .true.
    end if

    laverage = .true.
    if( present(message3) ) then
       if( trim(message3) == 'DO NOT AVERAGE') laverage = .false.
    end if

    nullify(ivalue)

    if( in_my_code .and. PAR_MY_CODE_RANK == 0 ) then
       allocate( ivalue(0:PAR_WORLD_SIZE-1) )
       ivalue = 0_ip       
    else if( .not. in_my_code .and. PAR_MY_WORLD_RANK == 0 ) then
       allocate( ivalue(0:PAR_WORLD_SIZE-1) )
       ivalue = 0_ip
    end if

    if( laverage ) then
       ave_ivalue = my_ivalue
       call PAR_AVERAGE(ave_ivalue,'IN MY CODE')
       if( IMASTER ) then
          dummi = ave_ivalue
       else
          dummi = my_ivalue
       end if
    else
       dummi = my_ivalue
    end if

    if( in_my_code ) then
       call PAR_GATHER(dummi,ivalue,'IN MY CODE')
    else
       call PAR_GATHER(dummi,ivalue,'IN THE WORLD') 
    end if

    if( PAR_MY_WORLD_RANK == 0 ) then
       write(lun_parti_res,'(a,i4,a)') 'Result '//trim(message1)//' '//trim(message2)//' ',jtime,' Scalar OnNodes'
       write(lun_parti_res,'(a)') 'ComponentNames '//trim(message1)
       write(lun_parti_res,'(a)') 'Values'
       do ipart = 0,PAR_WORLD_SIZE-1  
          write(lun_parti_res,'(1x,i7,1x,e13.6)') ipart+1,real(ivalue(ipart),rp)
          !write(lun_parti_res,'(1x,i7,1x,i12)') ipart+1,ivalue(ipart)
       end do
       write(lun_parti_res,'(a)') 'End Values'  
    end if

    if( laverage .and. PAR_MY_WORLD_RANK == 0 ) then
       ivalue(PAR_MY_WORLD_RANK) = 0
    end if

    if( PAR_MY_WORLD_RANK == 0 .and. present(ivalue_gat) ) ivalue_gat(0:PAR_WORLD_SIZE-1) = ivalue(0:PAR_WORLD_SIZE-1)
    if( associated(ivalue)  ) deallocate( ivalue  )

  end subroutine par_output_value_I

  subroutine par_output_value_R(jtime,my_rvalue,message1,message2,message3,wherein)

    integer(ip),  intent(in)           :: jtime
    real(rp),     intent(in)           :: my_rvalue
    character(*), intent(in)           :: message1
    character(*), intent(in)           :: message2
    character(*), intent(in), optional :: message3
    real(rp),     pointer              :: rvalue(:)
    character(*), intent(in), optional :: wherein
    integer(ip)                        :: ipart
    real(rp)                           :: ave_rvalue,dummr
    logical(lg)                        :: laverage
    logical(lg)                        :: in_my_code

    in_my_code = .false.
    if( present(wherein) ) then
       if(trim(wherein) == 'IN MY CODE' ) in_my_code = .true.
    end if

    laverage = .true.
    if( present(message3) ) then
       if( trim(message3) == 'DO NOT AVERAGE') laverage = .false.
    end if

    nullify(rvalue)

    if( in_my_code .and. PAR_MY_CODE_RANK == 0 ) then
       allocate( rvalue(0:PAR_WORLD_SIZE-1) )
       rvalue = 0_ip       
    else if( .not. in_my_code .and. PAR_MY_WORLD_RANK == 0 ) then
       allocate( rvalue(0:PAR_WORLD_SIZE-1) )
       rvalue = 0_ip
    end if

    if( laverage ) then
       ave_rvalue = my_rvalue
       call PAR_AVERAGE(ave_rvalue,'IN MY CODE')
       if( IMASTER ) then
          dummr = ave_rvalue
       else
          dummr = my_rvalue
       end if
    else
       dummr = my_rvalue
    end if

    if( in_my_code ) then
       call PAR_GATHER(dummr,rvalue,'IN MY CODE')
    else
       call PAR_GATHER(dummr,rvalue,'IN THE WORLD') 
    end if

    if( PAR_MY_WORLD_RANK == 0 ) then
       write(lun_parti_res,'(a,i4,a)') 'Result '//trim(message1)//' '//trim(message2)//' ',jtime,' Scalar OnNodes'
       write(lun_parti_res,'(a)') 'ComponentNames '//trim(message1)
       write(lun_parti_res,'(a)') 'Values'
       do ipart = 0,PAR_WORLD_SIZE-1  
          write(lun_parti_res,'(1x,i7,1x,e13.6)') ipart+1,rvalue(ipart)
       end do
       write(lun_parti_res,'(a)') 'End Values'    
    end if

    if( associated(rvalue)  ) deallocate( rvalue  )

  end subroutine par_output_value_R

  !----------------------------------------------------------------------
  !> @addtogroup Parall
  !> @{
  !> @file    par_output_global_matrix.f90
  !> @author  Guillaume Houzeaux
  !> @date    20/04/2014
  !> @brief   Output global matrix
  !> @details Output global matrix in GiD format
  !> @} 
  !----------------------------------------------------------------------

  subroutine par_output_global_matrix()


    integer(ip), parameter       :: BLOCK_AII=1
    integer(ip), parameter       :: BLOCK_AIB=2
    integer(ip), parameter       :: BLOCK_ABI=3
    integer(ip), parameter       :: BLOCK_ABB=4
    integer(ip)                  :: ielem,ipoin,ipart
    integer(ip)                  :: nboxes,isize,iboxe,jboxe
    integer(ip)                  :: ipoin1,ipoin2,ipoin3,ipoin4
    integer(ip)                  :: current_block,jboun,jpart
    integer(ip)                  :: ncomb,icomb,icomb_max
    integer(ip)                  :: jj,jneig,jpoin,iboun
    integer(ip)                  :: icode,ineig,ii
    integer(ip)                  :: icombi1,icombi2
    integer(ip)                  :: ncomb_gat,icomb_gat
    MY_MPI_COMM                  :: PAR_COMM_TO_USE
    real(rp)                     :: offsetx,offsety,xx,yy

    type(comm_data_par), pointer :: commu
    integer(ip),         pointer :: int_nodes_gat(:)
    logical(lg),         pointer :: gcombi(:,:)
    integer(ip),         pointer :: lcombi(:,:)
    integer(ip),         pointer :: scombi(:)
    integer(ip),         pointer :: lcombi_gat(:,:,:)
    integer(ip),         pointer :: scombi_gat(:,:)
    integer(ip),         pointer :: num_neigh(:)
    integer(ip),         pointer :: lcombi_tmp(:,:)
    integer(ip),         pointer :: scombi_tmp(:,:)


    if( ISEQUEN ) return

    icode = -1_ip 

    if( kfl_matri_par == 1 ) then

       call livinf(0_ip,'PARALL: OUTPUT GLOBAL MATRIX',0_ip)

       nullify( int_nodes_gat  )

       nullify(gcombi)
       nullify(lcombi)
       nullify(scombi)
       nullify(num_neigh)
       nullify(lcombi_tmp)

       nullify(lcombi_gat)
       nullify(scombi_gat)

       if( PAR_MY_WORLD_RANK == 0 ) then
          call memory_alloca(output_memor,'int_nodes_gat','par_output_global_matrix',int_nodes_gat,PAR_WORLD_SIZE,'INITIALIZE',0_ip)

          !allocate( int_nodes_gat (0:PAR_WORLD_SIZE-1) )
       end if
       !
       ! Number of combinations
       !
       if( INOTMASTER ) then

          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)

          call memory_alloca(output_memor,'gcombi','par_output_global_matrix',gcombi,commd % bound_dim,commu % nneig)
          call memory_alloca(output_memor,'scombi','par_output_global_matrix',scombi,commd % bound_dim)

          do ineig = 1,commu % nneig
             do ii = commd % bound_size(ineig),commd % bound_size(ineig+1)-1
                gcombi(ii,ineig) = .true.
             end do
          end do

          do ineig = 1,commu % nneig
             do ii = commd % bound_size(ineig),commd % bound_size(ineig+1)-1
                if( gcombi(ii,ineig) ) then
                   ipoin = commd % bound_perm(ii)
                   do jneig = ineig+1,commu % nneig
                      loop_jj: do jj = commd % bound_size(jneig),commd % bound_size(jneig+1)-1
                         jpoin = commd % bound_perm(jj)
                         if( ipoin == jpoin ) then
                            gcombi(ii,jneig) = .true.
                            gcombi(jj,jneig) = .false.
                            exit loop_jj
                         end if
                      end do loop_jj
                   end do
                end if
             end do
          end do

          ncomb = 0
          do ii = 1,commd % bound_dim
             if( any(gcombi(ii,:)) ) then
                icomb = 1          
                do jj = ii+1,commd % bound_dim
                   if( all(gcombi(ii,:) .eqv. gcombi(jj,:) ) ) then
                      icomb = icomb + 1
                      gcombi(jj,:) = .false.
                   end if
                end do
                ncomb           = ncomb + min(1_ip,icomb)
                gcombi(ncomb,:) = gcombi(ii,:)
                scombi(ncomb)   = scombi(ncomb) + icomb
             end if
          end do

          icomb_max = 0
          do icomb = 1,ncomb
             icomb_max = max(icomb_max,count( gcombi(icomb,:) .eqv. .true. ,KIND=ip))
          end do
          icomb_max = icomb_max + 1

          call memory_alloca(output_memor,'lcombi','par_output_global_matrix',lcombi,ncomb,icomb_max)

          do icomb = 1,ncomb
             ii = 1
             lcombi(icomb,ii) = kfl_paral
             do ineig = 1,commu % nneig
                if( gcombi(icomb,ineig) ) then
                   ii = ii + 1
                   lcombi(icomb,ii) = commu % neights(ineig)
                end if
             end do
          end do

          do icomb = 1,ncomb
             if( kfl_paral < maxval(lcombi(icomb,:)) ) then
                lcombi(icomb,:) = 0
                scombi(icomb)   = 0
             end if
          end do

          !do icomb = 1,ncomb
          !   !write(6,'(a,10(1x,i3))') 'a=',kfl_paral,lcombi(icomb,1:icomb_max)
          !   write(6,'(a,10(1x,i5))') 'a=',kfl_paral,lcombi(icomb,1:icomb_max),scombi(icomb)
          !end do
       end if

       ncomb_gat = ncomb
       icomb_gat = icomb_max
       call PAR_MAX(ncomb_gat)
       call PAR_MAX(icomb_gat)

       call memory_resize(output_memor,'lcombi'    ,'par_output_global_matrix',lcombi,ncomb_gat,icomb_gat)
       call memory_resize(output_memor,'scombi'    ,'par_output_global_matrix',scombi,ncomb_gat)
       if( IMASTER ) then
          call memory_alloca(output_memor,'lcombi_gat','par_output_global_matrix',lcombi_gat,ncomb_gat,icomb_gat,npart+1_ip,'INITIALIZE',1_ip,1_ip,0_ip)
          call memory_alloca(output_memor,'scombi_gat','par_output_global_matrix',scombi_gat,ncomb_gat,npart+1_ip,'INITIALIZE',1_ip,0_ip)
       end if

       call PAR_GATHER(lcombi,lcombi_gat,'IN THE WORLD')       
       call PAR_GATHER(scombi,scombi_gat,'IN THE WORLD')       

       if( IMASTER ) then 
          nboxes = PAR_WORLD_SIZE-1
          ncomb  = 0
          do ipart = 1,PAR_WORLD_SIZE-1
             do icomb = 1,ncomb_gat
                isize = 0
                if( scombi_gat(icomb,ipart) /= 0 ) then 

                   loop_isize: do while( isize < icomb_gat )
                      isize = isize + 1
                      if( lcombi_gat(icomb,isize,ipart) == 0 ) then
                         isize = isize - 1
                         exit loop_isize                     
                      end if
                   end do loop_isize

                   nboxes = nboxes + 1      
                   ncomb  = ncomb + 1
                   call maths_heap_sort(2_ip,isize,lcombi_gat(icomb,1:isize,ipart))
                end if
             end do
          end do

          call memory_resize(output_memor,'lcombi','par_output_global_matrix',lcombi,ncomb,icomb_gat)
          call memory_resize(output_memor,'scombi','par_output_global_matrix',scombi,ncomb)
          lcombi = 0
          ncomb = 0
          do ipart = 1,PAR_WORLD_SIZE-1
             do icomb = 1,ncomb_gat
                if( scombi_gat(icomb,ipart) /= 0 ) then 
                   ncomb = ncomb + 1
                   lcombi(ncomb,:) = lcombi_gat(icomb,:,ipart)
                   scombi(ncomb)   = scombi_gat(icomb,ipart)
                end if
             end do
          end do

       end if
       !
       ! Order
       !
       if( IMASTER .and. ncomb > 1 ) then

          call maths_heap_sort(1_ip,ncomb,scombi,lcombi)

          call memory_alloca(output_memor,'num_neigh','par_output_global_matrix',num_neigh,ncomb)

          do icomb = 1,ncomb
             num_neigh(icomb) = count( lcombi(icomb,:) /= 0 ,KIND=ip)
          end do

          icombi1 = 1
          icomb   = 0
          ineig   = scombi(1)
          do while( icomb < ncomb )
             icomb = icomb + 1
             if( scombi(icomb) /= ineig .or. icomb == ncomb ) then
                icombi2 = icomb - 1
                ii      = icombi2-icombi1+1
                allocate( lcombi_tmp(ii,icomb_gat) )
                allocate( scombi_tmp(ii,1) )
                lcombi_tmp(1:ii,1:icomb_gat) = lcombi(icombi1:icombi2,1:icomb_gat)
                scombi_tmp(1:ii,1)           = scombi(icombi1:icombi2)
                call maths_heap_sort(1_ip,ii,num_neigh(icombi1:),lcombi_tmp,scombi_tmp)
                lcombi(icombi1:icombi2,1:icomb_gat) = lcombi_tmp(1:ii,1:icomb_gat) 
                scombi(icombi1:icombi2)             = scombi_tmp(1:ii,1) 
                icombi1 = icomb
                ineig   = scombi(icomb)
                deallocate( lcombi_tmp , scombi_tmp )
             end if
          end do
       end if
       !
       ! Interior nodes
       !
       call PAR_GATHER(npoi1,int_nodes_gat,'IN THE WORLD') 
       !
       ! Coordinates (Master are also present!)
       ! 
       if( IMASTER ) then
          write(lun_matri_msh,'(a,i1,a)') 'MESH CODE_'//trim(intost(icode))//' dimension ',2_ip,' Elemtype Quadrilateral Nnode 4'   
          write(lun_matri_msh,'(a)') 'coordinates'  

          ipoin   = 0
          offsety = 0.0_rp
          do jboxe = 1,nboxes+1

             yy = - offsety
             if( jboxe <= nboxes ) then
                if( jboxe <= npart ) then
                   offsety = offsety + real(int_nodes_gat(jboxe),rp)
                else
                   offsety = offsety + real(scombi(jboxe-npart),rp)
                end if
             end if

             offsetx = 0.0_rp
             xx = 0.0_rp
             do iboxe = 1,nboxes+1 

                xx = offsetx

                if( iboxe <= nboxes ) then
                   if( iboxe <= npart ) then
                      offsetx = offsetx + real(int_nodes_gat(iboxe),rp)
                   else
                      offsetx = offsetx + real(scombi(iboxe-npart),rp)
                   end if
                end if

                ipoin = ipoin + 1
                write(lun_matri_msh,'(i9,3(1x,e12.6))') ipoin,xx,yy
             end do

          end do

          write(lun_matri_msh,'(a)') 'end coordinates'
          !
          ! Elements: edges of the communication
          !
          write(lun_matri_msh,'(a)') 'elements'
          ielem = 0 
          do jboxe = 1,nboxes
             do iboxe = 1,nboxes
                ipoin1 = jboxe*(nboxes+1)+iboxe
                ipoin2 = ipoin1 + 1
                ipoin4 = ipoin1 - (nboxes+1) 
                ipoin3 = ipoin4 + 1
                ielem  = ielem + 1

                if(      iboxe <= npart .and. jboxe <= npart ) then
                   current_block = BLOCK_AII
                else if( iboxe <= npart .and. jboxe > npart  ) then
                   current_block = BLOCK_AIB
                else if( iboxe > npart  .and. jboxe > npart  ) then
                   current_block = BLOCK_ABB
                else
                   current_block = BLOCK_ABI
                end if

                ii = 0
                if(      current_block == BLOCK_AII ) then
                   if( iboxe == jboxe ) ii = 1

                else if( current_block == BLOCK_AIB ) then
                   jboun = jboxe - npart 
                   ipart = iboxe
                   if( any(lcombi(jboun,:)==ipart) ) ii = 2

                else if( current_block == BLOCK_ABI ) then
                   iboun = iboxe - npart 
                   jpart = jboxe
                   if( any(lcombi(iboun,:)==jpart) ) ii = 2

                else if( current_block == BLOCK_ABB ) then
                   iboun = iboxe - npart  
                   jboun = jboxe - npart 
                   if( par_any_included(icomb_gat,icomb_gat,lcombi(iboun,:),lcombi(jboun,:)) ) ii = 3

                end if

                write(lun_matri_msh,'(10(1x,i9))') ielem,ipoin1,ipoin2,ipoin3,ipoin4,ii
             end do
          end do
          write(lun_matri_msh,'(a)') 'end elements'

          write(lun_matri_res,*) 'GiD Post Results File 1.0'
          write(lun_matri_res,*) ' '

          write(lun_matri_res,*) 'GaussPoints GP_QUAD4 Elemtype Quadrilateral'
          write(lun_matri_res,*) 'Number of Gauss Points: 1'
          write(lun_matri_res,*) 'Natural Coordinates: Internal'
          write(lun_matri_res,*) 'End GaussPoints'

          write(lun_matri_res,*) 'Result 1D_COORD 1D_COORD 0 Scalar OnGaussPoints GP_QUAD4'
          write(lun_matri_res,*) 'ComponentNames 1D_COORD'
          write(lun_matri_res,*) 'values'
          ielem = 0
          do jboxe = 1,nboxes
             do iboxe = 1,nboxes
                ielem = ielem + 1

                if(      iboxe <= npart .and. jboxe <= npart ) then
                   current_block = BLOCK_AII
                else if( iboxe <= npart .and. jboxe > npart  ) then
                   current_block = BLOCK_AIB
                else if( iboxe > npart  .and. jboxe > npart  ) then
                   current_block = BLOCK_ABB
                else
                   current_block = BLOCK_ABI
                end if

                ii = 0
                if(      current_block == BLOCK_AII ) then
                   if( iboxe == jboxe ) ii = int_nodes_gat(iboxe) ** 2

                else if( current_block == BLOCK_AIB ) then
                   jboun = jboxe - npart 
                   ipart = iboxe
                   if( any(lcombi(jboun,:)==ipart) ) ii = scombi(jboun) * int_nodes_gat(ipart) 

                else if( current_block == BLOCK_ABI ) then
                   iboun = iboxe - npart 
                   jpart = jboxe
                   if( any(lcombi(iboun,:)==jpart) ) ii = scombi(iboun) * int_nodes_gat(jpart) 

                else if( current_block == BLOCK_ABB ) then
                   iboun = iboxe - npart  
                   jboun = jboxe - npart 
                   if( par_any_included(icomb_gat,icomb_gat,lcombi(iboun,:),lcombi(jboun,:)) ) ii = scombi(iboun) * scombi(jboun) 

                end if

                write(lun_matri_res,'(10(1x,i9))') ielem,ii
             end do
          end do
          write(lun_matri_res,*) 'end values' 

          call iofile_flush_unit(lun_matri_msh)
          call iofile_flush_unit(lun_matri_res)

       end if
       !
       ! Deallocate
       !
       call memory_deallo(output_memor,'int_nodes_gat','par_output_global_matrix',int_nodes_gat)
       call memory_deallo(output_memor,'gcombi',       'par_output_global_matrix',gcombi)
       call memory_deallo(output_memor,'scombi',       'par_output_global_matrix',scombi)
       call memory_deallo(output_memor,'lcombi',       'par_output_global_matrix',lcombi)
       call memory_deallo(output_memor,'num_neigh',    'par_output_global_matrix',num_neigh)
       call memory_deallo(output_memor,'scombi_gat',   'par_output_global_matrix',scombi_gat)
       call memory_deallo(output_memor,'lcombi_gat',   'par_output_global_matrix',lcombi_gat)


    end if

  end subroutine par_output_global_matrix

  function par_any_nonordered(array1,array2,my_mask)

    integer(ip), intent(in), pointer  :: array1(:)
    integer(ip), intent(in), pointer  :: array2(:)
    integer(ip), intent(in), optional :: my_mask
    logical(lg)                       :: par_any_nonordered
    integer(ip)                       :: ii,jj,nn1,nn2

    par_any_nonordered = .false.
    nn1 = memory_size(array1)
    nn2 = memory_size(array2)

    if( present(my_mask) ) then
       do ii = 1,nn1
          if( array1(ii) /= 0 ) then
             do jj = 1,nn2
                if( array1(ii) == array2(jj) ) then
                   par_any_nonordered = .true.
                   return
                end if
             end do
          end if
       end do
    else
       do ii = 1,nn1
          do jj = 1,nn2
             if( array1(ii) == array2(jj) ) then
                par_any_nonordered = .true.
                return
             end if
          end do
       end do
    end if

  end function par_any_nonordered

  function par_any_included(nn1,nn2,array1,array2)

    integer(ip), intent(in)  :: nn1
    integer(ip), intent(in)  :: nn2
    integer(ip), intent(in)  :: array1(nn1)
    integer(ip), intent(in)  :: array2(nn2)
    logical(lg)              :: par_any_included
    integer(ip)              :: ii,jj,kk,ii0,jj0

    par_any_included = .false.

    ii0 = nn1
    loop_ii: do ii = 1,nn1
       if( array1(ii) == 0 ) then
          ii0 = ii - 1
          exit loop_ii
       end if
    end do loop_ii

    jj0 = nn2
    loop_jj: do jj = 1,nn2
       if( array2(jj) == 0 ) then
          jj0 = jj - 1
          exit loop_jj
       end if
    end do loop_jj

    if( ii0 <= jj0  ) then
       ii  = 0
       kk  = 0
       do while( ii < nn1 )
          ii = ii + 1
          if( array1(ii) /= 0 ) then
             jj = 0
             do while( jj < nn2 )
                jj = jj + 1
                if( array1(ii) == array2(jj) ) then
                   kk = kk + 1
                   jj = nn2
                end if
             end do
          else
             ii  = nn1
          end if
       end do
    else
       ii  = 0
       kk  = 0
       do while( ii < nn2 )
          ii = ii + 1
          if( array2(ii) /= 0 ) then
             jj = 0
             do while( jj < nn1 )
                jj = jj + 1
                if( array2(ii) == array1(jj) ) then
                   kk = kk + 1
                   jj = nn1
                end if
             end do
          else
             ii  = nn2
          end if
       end do
    end if

    if( kk == min(ii0,jj0)) par_any_included = .true.

  end function par_any_included

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-10-07
  !> @brief   Open result file
  !> @details Open result file and append
  !> 
  !-----------------------------------------------------------------------

  subroutine par_output_partition_open(MESH,RESULTS,APPEND)

    logical(lg), optional, intent(in) :: MESH
    logical(lg), optional, intent(in) :: RESULTS
    logical(lg), optional, intent(in) :: APPEND
    integer(ip)                       :: ii,jj
    logical(lg)                       :: if_mesh
    logical(lg)                       :: if_results
    logical(lg)                       :: if_append

    if_mesh    = optional_argument(.false.,MESH)
    if_results = optional_argument(.false.,RESULTS)
    if_append  = optional_argument(.false.,APPEND)

    if( PAR_MY_WORLD_RANK == 0 ) then
       if( kfl_amr /= 0 ) then
          if( if_mesh ) then
             ii = index(fil_parti_msh,'.post')
             jj = len(trim(fil_parti_msh))
             call iofile_open_unit(lun_parti_msh,fil_parti_msh(1:ii-1)//'-amr'//integer_to_string(num_amr)&
                  //trim(fil_parti_msh(ii:jj)),'PARALL PARTITION MESH')
          end if
          if( if_results ) then
             ii = index(fil_parti_res,'.post')
             jj = len(trim(fil_parti_res))
             if( if_append ) then
                call iofile_open_unit(lun_parti_res,fil_parti_res(1:ii-1)//'-amr'//integer_to_string(num_amr)&
                     //trim(fil_parti_res(ii:jj)),'PARALL PARTITION RESULT',POSIO='APPEND')
             else
                call iofile_open_unit(lun_parti_res,fil_parti_res(1:ii-1)//'-amr'//integer_to_string(num_amr)&
                     //trim(fil_parti_res(ii:jj)),'PARALL PARTITION RESULT')
             end if
          end if
       else
          if( if_mesh ) then
             call iofile_open_unit(lun_parti_msh,fil_parti_msh,'PARALL PARTITION MESH')
          end if
          if( if_results ) then
             if( if_append ) then
                call iofile_open_unit(lun_parti_res,fil_parti_res,'PARALL PARTITION RESULT',POSIO='APPEND')
             else
                call iofile_open_unit(lun_parti_res,fil_parti_res,'PARALL PARTITION RESULT')
             end if
          end if
       end if
    end if

  end subroutine par_output_partition_open
    
end module mod_par_output_partition
!> @}
