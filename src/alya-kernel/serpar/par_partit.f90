!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_partit()
  !-------------------------------------------------------------------------------
  !****f* Parall/par_arrays
  ! NAME
  !    par_arrays
  ! DESCRIPTION
  !
  ! INPUT
  !    Element graph
  ! OUTPUT
  !    Partition of the graph
  !    lepar_par
  !    leper_par
  !    leinv_par
  !    lnpar_par
  !    lnper_par
  !    lninv_par
  !    lneig_par
  !    ginde_par 
  !    lcomm_par
  ! USED BY
  !    par_partit
  !***
  !-------------------------------------------------------------------------------
  use def_parame 
  use def_elmtyp
  use def_domain 
  use def_parall
  use def_master
  use mod_memory
  use mod_parall
  use def_kintyp,                    only : ip, rp
  use mod_partition_sfc,             only : partition_sfc
  use mod_redistribute,              only : par_redistribute
  use mod_redistribute,              only : gather_to_master
  use mod_partitioning,              only : partitioning_oriented_bin
  use mod_partitioning,              only : partitioning_metis
  use def_domain,                    only : ndime
  use mod_parall,                    only : PAR_COMM_MY_CODE_WM
  use mod_par_parallel_partitioning, only : par_partition_weights
  use mod_alya2metis,                only : alya2metis_METIS_PartGraph
  implicit none
  integer(ip), pointer    :: wvert_par(:)
  integer(ip)             :: wedge_par(1)   
  integer(ip)             :: idime,pnode,npart_sfc
  integer(ip)             :: dummi,ielem
  integer(ip)             :: kfl_weigh_metis_par
  integer(ip), pointer    :: ladja_tmp(:)
  integer(ip), pointer    :: padja_tmp(:)
  integer(ip), pointer    :: wvert_tmp(:)
  integer(ip), pointer    :: permr_tmp(:)
  integer(ip), pointer    :: invpr_tmp(:)
  integer(ip), pointer    :: npart_tmp(:)
  integer(ip), pointer    :: lepar_tmp(:)
  real(rp),    pointer    :: coord_elem(:,:)
  real(rp)                :: time0,time1,time2
  !
  ! Intermediate variables for partitioning
  !
  integer(ip)             :: nenti_
  real(rp), pointer       :: lenti_(:,:) 
  integer(ip), pointer    :: lweig_(:)   
  integer(4),  pointer    :: nelem_part_master4_(:) 
  integer(4),  pointer    :: npoin_part_master4_(:) 

  nullify(lenti_,lweig_,nelem_part_master4_,npoin_part_master4_)
  !
  ! Nullify pointers
  !
  nullify(ladja_tmp)
  nullify(padja_tmp)
  nullify(wvert_tmp)
  nullify(permr_tmp)
  nullify(invpr_tmp)
  nullify(npart_tmp) 
  nullify(lepar_tmp)
  nullify(wvert_par)
  nullify(coord_elem)

  !----------------------------------------------------------------------
  !
  ! Compute Vertex weights WVERT_PAR = # of Gauss points of the element
  !
  !----------------------------------------------------------------------
  !kfl_weigh_par = 5
  !print *, '## -->', kfl_weigh_par
  call par_livinf(20_ip,'  ALLOCATE MEMORY...',dummi) 
  call memory_alloca(par_memor,'WVERT_PAR','par_prepro',wvert_par,nelem)
  !
  ! Weights for partitioning...
  ! If hole elements are taken off, one subdomain can end up with all the hole elements
  ! as they have null weight!
  ! 
  call par_livinf(20_ip,'  SET WEIGHTS... ',dummi)

  call par_partition_weights(wvert_par)
  !
  ! Define flag for METIS weights
  ! 0: No weights, 1: edges, 2: vertices, 3: both
  !
  if( kfl_weigh_par == -1 ) then
     kfl_weigh_metis_par = 0 
  else
     kfl_weigh_metis_par = 2
  end if

  wedge_par = 0

  !----------------------------------------------------------------------
  !
  ! Compute element graph: PADJA_PAR (PELEL) and LADJA_PAR (LELEL)
  ! using node or face connectivity
  !
  !----------------------------------------------------------------------

  call cputim(time0)
  call par_memory(1_ip)
  call par_memory(8_ip)
  call par_livinf(2_ip,' ',dummi)

  call par_elmgra()
  call par_livinf(7_ip,' (# EDGES= '//trim(intost(nedge))//', MAX # EDGES/ELEMENT= '//trim(intost(medge))//')',0_ip)

  call cputim(time1)
  cpu_paral(3) = time1 - time0

  !----------------------------------------------------------------------
  !
  ! Partition mesh with METIS
  !
  !----------------------------------------------------------------------

  call par_livinf(3_ip,' ',dummi)

  if( kfl_partition_par == PAR_METIS4 ) then
     !
     ! METIS
     !
     call alya2metis_METIS_PartGraph(npart_par,nelem,padja_par,ladja_par,wvert_par,lepar_par)
     
  else if( kfl_partition_par < 0 ) then
     !
     ! Read from a field
     !
     do ielem = 1,nelem 
        lepar_par(ielem) = int(xfiel(-kfl_partition_par) % a(1,ielem,1),ip)
     end do

  else if( kfl_partition_par == PAR_ORIENTED_BIN ) then
     !
     ! Oriented bin
     !
     call memory_alloca(par_memor,'COORD_ELEM','par_prepro',coord_elem,ndime,nelem)
     !$OMP  PARALLEL DO SCHEDULE (STATIC)            & 
     !$OMP  DEFAULT (NONE)                           &
     !$OMP  PRIVATE ( ielem,idime,pnode )            &
     !$OMP  SHARED  ( nelem,nnode,ltype,coord,lnods, &
#ifndef NDIMEPAR
     !$OMP            ndime,                         &
#endif
     !$OMP            coord_elem )

     do ielem = 1,nelem
        pnode = nnode(abs(ltype(ielem)))
        do idime = 1,ndime 
           coord_elem(idime,ielem) = sum(coord(idime,lnods(1:pnode,ielem)))
        end do
        coord_elem(1:ndime,ielem) = coord_elem(1:ndime,ielem) / real(pnode,rp)
     end do
     !$OMP END PARALLEL DO

     call partitioning_oriented_bin(&
          ndime,npart_par,nelem,boxes_fine_par,vect_partition_par, &
          coord_elem,lepar_par,'LOAD BALANCE WITH OPENMP',wvert_par)

     call memory_deallo(par_memor,'COORD_ELEM','par_prepro',coord_elem)

  else if( kfl_partition_par == PAR_SFC ) then
     !
     ! SFC
     !
     npart_sfc = PAR_WORLD_SIZE - 1_ip 

     call par_redistribute(npart_sfc,lenti_,lweig_,nenti_,nelem_part_master4_,npoin_part_master4_)
     !call partition_sfc(lepar_par,npart_sfc,lenti_,lweig_,nenti_,ndime,PAR_COMM_MY_CODE_WM)
     call gather_to_master(lepar_par,npart)

     if(associated(lenti_)) deallocate(lenti_)
     if(associated(lweig_)) deallocate(lweig_)
     if(associated(nelem_part_master4_)) deallocate(nelem_part_master4_)
     if(associated(npoin_part_master4_)) deallocate(npoin_part_master4_)

  end if

  call cputim(time2)
  cpu_paral(5) = time2 - time1

  !----------------------------------------------------------------------
  !
  ! Deallocate memory
  !
  !----------------------------------------------------------------------

  call memory_deallo(par_memor,'WVERT_PAR','par_prepro_D',wvert_par)

end subroutine par_partit
