!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine inidom()
  !-----------------------------------------------------------------------
  !****f* domain/inidom
  ! NAME
  !    inidom
  ! DESCRIPTION
  !    This routines initialize domain variables.
  ! USED BY
  !    domain
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_meshin, only : memor_msh
  use def_inpout
  use mod_iofile
  use mod_memchk
  use mod_ecoute, only : ecoute_set_read_unit
  use mod_ecoute, only : ecoute_set_write_unit

  implicit none
  integer(ip) :: ielty,ii
  !
  ! Nullify pointers from reageo
  !
  nullify(lnods)
  nullify(ltype)
  nullify(lnodb)
  nullify(ltypb)
  nullify(lboch)
  nullify(lelbo)
  nullify(lnnod)
  nullify(lelch)
  nullify(lmate)
  nullify(lesub)
  nullify(lnoch)
  nullify(lgrou_dom)
  nullify(lperi)
  nullify(coord)
  nullify(cooin)
  nullify(skcos)
  nullify(lhang)
  !
  ! Nullify pointers from readim
  !
  do ii = 1,mfiel
     nullify(time_field(ii) % a)
     nullify(xfiel(ii) % a)
     k_tran_fiel(ii) = 0_ip
     k_tran_fiel_real(ii) = 0_ip
  end do
  !
  ! Nullify pointers from reaset
  !
  nullify(leset)
  nullify(lbset)
  nullify(lnset)
  nullify(lesec)
  nullify(lbsec)
  nullify(lnsec)
  !
  ! Nullify pointers from reabcs
  !
  nullify(kfl_codno)
  nullify(kfl_coded)
  nullify(kfl_codbo)
  !
  ! Nullify zones and nodal materials
  !
  nullify(lpoib)
  nullify(nmatn)
  nullify(lbono)
  nullify(lmatn)
  nullify(lnnob)
  !
  ! Nullify other pointers
  !
  nullify(lgaus)
  nullify(kfl_coded)
  nullify(nepoi)
  nullify(pelpo)
  nullify(lelpo)
  nullify(pelpo_2)
  nullify(lelpo_2)
  nullify(pelel)
  nullify(lelel)
  nullify(pelel_2)
  nullify(lelel_2)
  nullify(lezdo)
  nullify(lbzdo)
  nullify(r_sol)
  nullify(c_sol)
  nullify(r_dom)
  nullify(c_dom)
  nullify(r_bou)
  nullify(c_bou)
  nullify(r_dom_aii)
  nullify(permr_aii)
  nullify(invpr_aii)
  nullify(c_dom_aii)
  nullify(r_dom_aib)
  nullify(c_dom_aib)
  nullify(r_dom_abi)
  nullify(c_dom_abi)
  nullify(r_dom_abb)
  nullify(c_dom_abb)
  nullify(permr_abb)
  nullify(invpr_abb)
  nullify(r_dom_prec)
  nullify(c_dom_prec)
  nullify(permr_prec)
  nullify(invpr_prec)
  nullify(r_dom_own)
  nullify(c_dom_own)
  nullify(r_sym)
  nullify(c_sym)
  nullify(ompss_domains)
  nullify(ompss_boundaries)
  nullify(lpoty)
  nullify(lfacs)
  nullify(lboel)
  nullify(lgaib)
  nullify(lrenn)
  nullify(lfcnt)
  nullify(lncnt)
  nullify(lessl)
  nullify(kfl_fixno)
  nullify(kfl_fixbo)
  nullify(kfl_funno)
  nullify(kfl_funbo)
  nullify(kfl_funtn)
  nullify(kfl_funtb)
  nullify(kfl_fixrs)
  nullify(kfl_geobo)
  nullify(kfl_geono)
  nullify(lpoin)
  nullify(bvess)
  nullify(bvnat)
  nullify(tncod)
  nullify(tgcod)
  nullify(exnor)
  nullify(lobas)
  nullify(vmass)
  nullify(vmasc)
  nullify(cmass)
  nullify(cmass_weighted)
  nullify(dmass)
  nullify(mass_right_fact)
  nullify(walld)
  nullify(wallo)
  nullify(wallcoor)
  nullify(rough)
  nullify(canhe)
  nullify(canla)
  nullify(heiov)
  nullify(ywalb)
  nullify(ywalp)
  nullify(yscab)
  nullify(yscap)
  nullify(elmar)
  nullify(pefpo)
  nullify(lefpo)
  nullify(lnuew)
  nullify(leldo)
  nullify(facel)
  nullify(lnlev)
  nullify(lelev)
  nullify(lblev)
  nullify(lmast)
  nullify(meshe)
  nullify(lpmsh)
  nullify(lemsh)
  nullify(lbmsh)
  nullify(cutel)
  nullify(element_bin)
  nullify(vomat)

  nullify(xscal_fields)
  
  nullify(materials_nlaye)            ! Automatic generation of materials
  nullify(materials_icode)            ! Automatic generation of materials
  nullify(materials_imate)            ! Automatic generation of materials

  !
  ! Edges
  !
  nullify(edge_to_node)
  nullify(ledgs)
  nullify(ledgb)
  nullify(lnned)
  nullify(lnneb)
  nullify(c_edg)
  nullify(r_edg)
  !
  ! Finite volume
  !
  nullify(fv_center_coord)
  nullify(fv_cell_volume)
  nullify(fv_face_area)
  nullify(fv_face_normal)
  nullify(fv_face_orientation)
  nullify(fv_center_distance)
  nullify(fv_center_vector)
  nullify(fv_center_face)
  nullify(fv_face_boundary)
  nullify(fv_face_graph)
  nullify(fv_graph_diag)
  !
  ! General
  !
  call ecoute_set_read_unit (lun_pdata_dom) ! Reading file
  call ecoute_set_write_unit(lun_outpu_dom) ! Writing file
  nzsky         = 0              ! Components of the matrix
  memor_dom     = 0              ! Memory counters
  nperi         = 0              ! # of periodic nodes
  necnt         = 0              ! Number of contact elements
  nncnt         = 0              ! Number of contact elements
  kfl_parex     = 0              ! No automatic parallelization
  nbono         = 0              ! Number of boundary nodes
  kfl_crbou     = 0              ! Boundary graph has not been computed
  kfl_pelel     = 0              ! No element graph
  kfl_domar     = 0              ! Do not recompute mesh arrays (do not go to domarr)
  kfl_domar_world = 0            ! Do not recompute añya world mesh arrays (like bin structure)
  bandw_dom     = 0              ! Bandwidth
  profi_dom     = 0.0_rp         ! Profile
  kfl_horde     = 0              ! High order element do not exist
  kfl_elcoh     = 0              ! No cohesive element
  kfl_elint     = 0              ! No interface element
  nfacg         = 0              ! Number of faces
  mface         = 0              ! Max number of faces
  num_lobas     = 0              ! Number of local basis
  ncodb         = 0              ! Number of boundary codes
  do ielty = 1,nelty
     lnuty(ielty) = 0
     lnuib(ielty) = 0
  end do
  kexist_tran_fiel = 0
  !
  ! Immersed boundary
  !
  mnoib         = 0              ! IB boundary nodes
  mnodi         = 0              ! IB element nodes
  mgaib         = 0              ! IB Gauss points
  mgaui         = 0              ! IB Gauss points
  kfl_savda     = 0              ! Do not save element data base
  !
  ! Mesh generator
  !
  memor_msh     = 0              !  Memory counters
  !
  ! Variables read in readim
  !
  kfl_autbo =  0                 ! Automatic boundaries off
  npoin     = -1                 ! Obligatory
  nelem     = -1                 ! Obligatory
#ifndef NDIMEPAR
  ndime     = -1                 ! Obligatory
#endif
  nboun     = -1                 ! Optional
  if( ISEQUEN ) then
     allocate( npoin_par(1) )
     allocate( nelem_par(1) )
     allocate( nboun_par(1) )
  end if
  npoin_origi = 0
  npoin_total = 0
  nelem_total = 0
  nboun_total = 0
  nedge     =  0
  medge     =  0
  nperi     =  0                 ! Optional
  nboib     =  0                 ! Optional: # IB
  npoib     =  0                 ! Optional: # IB points
  nimbo     =  0                 ! Optional
  nrbod     =  0                 ! Optional
  nmate     =  1                 ! No material
  nsubd     =  1                 ! No material
  nfiel     =  0                 ! Number of fields
  nzone     =  1                 ! Only one zone
  kfl_field =  0                 ! Field dimensions
  ngrou_dom =  0                 ! Number of groups
  !
  ! Variables read in reastr
  !
  lquad      = 0                 ! Open integration rule
  kfl_xscal_fields = 0           ! No field scaling
  !
  ! Variables set in reageo
  !
  kfl_chege = 0                  ! Don't check geometry
  kfl_naxis = 0                  ! Cartesian coordinate system
  kfl_spher = 0                  ! Cartesian coordinate system
  kfl_bouel = 0                  ! Element # connected to boundary unknown
  !
  ! Variable read in reabcs
  !
  mcono      = 3                 ! Max number of conditions per node
  kfl_icodn  = 0                 ! There is no node code array
  kfl_icode  = 0                 ! There is no edge code array
  kfl_icodb  = 0                 ! There is no boundary code array
  kfl_crbou  = 0                 ! Boundary graph allocated?
  !
  ! Element bin
  !
  element_bin_boxes(1:3) = 4     ! Read by kermod
  !
  ! Graph
  !
  nzdom      = 0                 ! Node graph
  nzdom_own  = 0                 ! Own node graph
  !
  ! Element volumes: negative values means not initialized: keep this!
  !
  voave = -1.0_rp
  vomin = -1.0_rp
  vomax = -1.0_rp
  !
  ! Element array
  !
  call cshder(1_ip)

end subroutine inidom

