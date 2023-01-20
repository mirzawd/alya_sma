!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_memall()
  !-----------------------------------------------------------------------
  !****f*  magnet/mag_memall
  ! NAME
  !    mag_memall
  ! DESCRIPTION
  !    This routine allocates memory for arrays needed for the module
  ! USES
  ! USED BY
  !    mag_turnon
  !***
  !-----------------------------------------------------------------------

  use def_master, only: mem_modul, modul, INOTMASTER
  use def_domain, only: ndime, meshe, ltype, ngaus
  use def_magnet
  use def_kermod, only: ndivi
  use mod_memory

  use mod_communications, only: PAR_SUM, PAR_ALLGATHER

  implicit none

  integer(ip) :: iconstr, iedge, jedge, pedge, ielem, pelty, pgaus
  integer(ip) :: nleleb, nledgb, ngedgb
  !
  call mag_memsol()
  !
  call mag_memcnt()
  !
  !-----------------------------------------------------------------------------
  ! Boundary data for self field calculations
  !-----------------------------------------------------------------------------
  !
  ! proc_nbedg array stores the number of local edges on the boundary for each processor
  !
  allocate(proc_nbedg(nproc_par))
  !
  ! Boolean marking Dirichlet edges: .true. if Dirichlet, otherwise .false.
  !
  call memory_alloca(mem_modul(1:2,modul), 'EDGDIR', 'mag_memall', edgdir_mag, meshe(ndivi) % nedge)
  !
  ! Boolean marking boundary edges: .true. if boundary edge, otherwise .false.
  !
  call memory_alloca(mem_modul(1:2,modul), 'EDGBOU', 'mag_memall', edgbou_mag, meshe(ndivi) % nedge)
  !
  ! nleleb - number of local elements on boundary
  ! nledgb - number of local edges on boundary
  ! ngedgb - number of global edges on boundary
  !
  nleleb = memory_size(meshe(ndivi) % ledgb, 2_ip)
  !
  do ielem = 1, nleleb
    ! Number of edges in current boundary element
    pedge = meshe(ndivi) % lnneb(ielem)
    do iedge = 1, pedge
      edgdir_mag(meshe(ndivi) % ledgb(iedge, ielem)) = .true.
      edgbou_mag(meshe(ndivi) % ledgb(iedge, ielem)) = .true.
    end do
  end do
  !
  if( nleleb > 0 ) then
     nledgb = count(edgbou_mag)
  else
     nledgb = 0
  end if
  ngedgb = nledgb
  !
  call PAR_ALLGATHER(int(nledgb * ndime, 4), proc_nbedg)
  !
  call PAR_SUM(ngedgb)
  !
  ! List of Dirichlet edges
  ! In 3-D some boundary edges will be repeated across processors
  ! This needs to be somehow taken into account in self field calculations????
  !
  call memory_alloca(mem_modul(1:2,modul), 'DIREDG', 'mag_memall', diredg_mag, nedgdir_mag)
  !
  call memory_alloca(mem_modul(1:2,modul), 'BOUEDG', 'mag_memall', bouedg_mag, nledgb)
  !
  jedge = 0
  do iedge = 1, meshe(ndivi) % nedge
    if (edgbou_mag(iedge)) then
      jedge = jedge + 1_ip
      bouedg_mag(jedge) = iedge
    end if
  end do

  call memory_alloca(mem_modul(1:2,modul), 'HBS', 'mag_memall', Hbs_mag, meshe(ndivi) % ndime, ngedgb)

  call memory_alloca(mem_modul(1:2,modul), 'HXTL', 'mag_memall', Hxtl_mag, nledgb)
  call memory_alloca(mem_modul(1:2,modul), 'HSFL', 'mag_memall', Hsfl_mag, nledgb)
  call memory_alloca(mem_modul(1:2,modul), 'HSFG', 'mag_memall', Hsfg_mag, ngedgb)

  call memory_alloca(mem_modul(1:2,modul), 'LOCBOUCEN', 'mag_memall', locboucen_mag, meshe(ndivi) % ndime, nledgb)
  call memory_alloca(mem_modul(1:2,modul), 'GLOBOUCEN', 'mag_memall', globoucen_mag, meshe(ndivi) % ndime, ngedgb)

  call memory_alloca(mem_modul(1:2,modul), 'LOCBOUTAN', 'mag_memall', locboutan_mag, meshe(ndivi) % ndime, nledgb)
  call memory_alloca(mem_modul(1:2,modul), 'GLOBOUTAN', 'mag_memall', globoutan_mag, meshe(ndivi) % ndime, ngedgb)
  !-----------------------------------------------------------------------------
  !
  if (INOTMASTER) then
    !-----------------------------------------------------------------------------
    ! Memory Allocation 
    !-----------------------------------------------------------------------------
    !
    ! Edge sign data
    !
    call memory_alloca(mem_modul(1:2,modul), 'ELESIG', 'mag_memall', elesig_mag, meshe(ndivi) % medge, meshe(ndivi) % nelem)
    !
    ! Edge global signs
    !
    call memory_alloca(mem_modul(1:2,modul), 'EDGSIG', 'mag_memall', edgsig_mag, meshe(ndivi) % nedge)
    !
    ! Edge length data
    !
    call memory_alloca(mem_modul(1:2,modul), 'EDGLEN', 'mag_memall', edglen_mag, meshe(ndivi) % nedge)
    !
    ! Edge connectivity matrix
    !
    call memory_alloca(mem_modul(1:2,modul), 'ELEEDG', 'mag_memall', eleedg_mag, meshe(ndivi) % medge, meshe(ndivi) % nelem)
    !
    ! Element area
    !
    call memory_alloca(mem_modul(1:2,modul), 'ELEVOL', 'mag_memall', elevol_mag, meshe(ndivi) % nelem)
    !
    ! Edge - Node connectivity (global)
    !
    call memory_alloca(mem_modul(1:2,modul), 'EDGNOD', 'mag_memall', edgnod_mag, 2_ip, meshe(ndivi) % nedge)
    !
    ! Solution field at time step n
    !
    call memory_alloca(mem_modul(1:2,modul), 'HE', 'mag_memall', he_mag, meshe(ndivi) % nedge)
    !
    ! Solution field at intermediate time step: Crank-Nicholson
    !
    call memory_alloca(mem_modul(1:2,modul), 'HA', 'mag_memall', ha_mag, meshe(ndivi) % nedge)
    !
    ! Solution field at time step n-1
    !
    call memory_alloca(mem_modul(1:2,modul), 'HP', 'mag_memall', hp_mag, meshe(ndivi) % nedge)
    !
    ! Solution update
    !
    call memory_alloca(mem_modul(1:2,modul), 'DH', 'mag_memall', dh_mag, meshe(ndivi) % nedge)
    !
    ! Constraints
    !
    do iconstr = 1_ip, constr_total
      call memory_alloca(mem_modul(1:2,modul), 'CONSTR_CNSTR', 'mag_memall', constr_mag(iconstr) % cnstr, meshe(ndivi) % nedge)
      call memory_alloca(mem_modul(1:2,modul), 'CONSTR_UNKNO', 'mag_memall', constr_mag(iconstr) % unkno, meshe(ndivi) % nedge)
    end do
    !
    ! bhsid & fhsid
    ! 
    call memory_alloca(mem_modul(1:2,modul), 'BHSID', 'mag_memall', bhsid, meshe(ndivi) % nedge)
    call memory_alloca(mem_modul(1:2,modul), 'FHSID', 'mag_memall', fhsid, meshe(ndivi) % nedge)
    !
    ! hhsid
    !
    call memory_alloca(mem_modul(1:2,modul), 'HHSID', 'mag_memall', hhsid, meshe(ndivi) % npoin)
    !
    ! Post-processing
    !
    call memory_alloca(mem_modul(1:2,modul), 'HC', 'mag_memall', Hc_mag, meshe(ndivi) % ndime, meshe(ndivi) % nelem)
    call memory_alloca(mem_modul(1:2,modul), 'HN', 'mag_memall', Hn_mag, meshe(ndivi) % ndime, meshe(ndivi) % npoin)
    call memory_alloca(mem_modul(1:2,modul), 'BC', 'mag_memall', Bc_mag, meshe(ndivi) % ndime, meshe(ndivi) % nelem)
    call memory_alloca(mem_modul(1:2,modul), 'BN', 'mag_memall', Bn_mag, meshe(ndivi) % ndime, meshe(ndivi) % npoin)
    call memory_alloca(mem_modul(1:2,modul), 'JC', 'mag_memall', Jc_mag, meshe(ndivi) % ndime, meshe(ndivi) % nelem)
    call memory_alloca(mem_modul(1:2,modul), 'JN', 'mag_memall', Jn_mag, meshe(ndivi) % ndime, meshe(ndivi) % npoin)
    call memory_alloca(mem_modul(1:2,modul), 'JCZ', 'mag_memall', Jcz_mag, meshe(ndivi) % nelem)
    call memory_alloca(mem_modul(1:2,modul), 'JNZ', 'mag_memall', Jnz_mag, meshe(ndivi) % npoin)
    call memory_alloca(mem_modul(1:2,modul), 'FC', 'mag_memall', Fc_mag, meshe(ndivi) % ndime, meshe(ndivi) % nelem)
    call memory_alloca(mem_modul(1:2,modul), 'FN', 'mag_memall', Fn_mag, meshe(ndivi) % ndime, meshe(ndivi) % npoin)
    call memory_alloca(mem_modul(1:2,modul), 'GC', 'mag_memall', Gc_mag, meshe(ndivi) % nelem)
    call memory_alloca(mem_modul(1:2,modul), 'GN', 'mag_memall', Gn_mag, meshe(ndivi) % npoin)
    !
    call memory_alloca(mem_modul(1:2,modul), 'HGP', 'mag_memall', Hgp_mag, meshe(ndivi) % nelem)
    call memory_alloca(mem_modul(1:2,modul), 'JGP', 'mag_memall', Jgp_mag, meshe(ndivi) % nelem)
    call memory_alloca(mem_modul(1:2,modul), 'BGP', 'mag_memall', Bgp_mag, meshe(ndivi) % nelem)
    call memory_alloca(mem_modul(1:2,modul), 'FGP', 'mag_memall', Fgp_mag, meshe(ndivi) % nelem)
    call memory_alloca(mem_modul(1:2,modul), 'GGP', 'mag_memall', Ggp_mag, meshe(ndivi) % nelem)
    !
    do ielem = 1, meshe(ndivi) % nelem
      pelty = ltype(ielem)
      pgaus = ngaus(pelty)
      call memory_alloca(mem_modul(1:2,modul), 'HGP_MAG(IELEM)', 'mag_memall', Hgp_mag(ielem)%a, meshe(ndivi) % ndime, pgaus, 1_ip)
      call memory_alloca(mem_modul(1:2,modul), 'BGP_MAG(IELEM)', 'mag_memall', Bgp_mag(ielem)%a, meshe(ndivi) % ndime, pgaus, 1_ip)
      call memory_alloca(mem_modul(1:2,modul), 'FGP_MAG(IELEM)', 'mag_memall', Fgp_mag(ielem)%a, meshe(ndivi) % ndime, pgaus, 1_ip)
      call memory_alloca(mem_modul(1:2,modul), 'GGP_MAG(IELEM)', 'mag_memall', Ggp_mag(ielem)%a, 1_ip, pgaus, 1_ip)

      if (ndime == 2_ip) then 
        call memory_alloca(mem_modul(1:2,modul), 'JGP_MAG(IELEM)', 'mag_memall', Jgp_mag(ielem)%a, 1_ip, pgaus, 1_ip)
      else
        call memory_alloca(mem_modul(1:2,modul), 'JGP_MAG(IELEM)', 'mag_memall', Jgp_mag(ielem)%a, meshe(ndivi) % ndime, pgaus, 1_ip)
      end if
    end do    
    !
    ! BDF
    !
    if (bdfode_mag % s >= 1_ip) call memory_alloca(mem_modul(1:2,modul), 'BDFODE_Y1', 'mag_memall', bdfode_mag % y1, meshe(ndivi) % nedge)
    if (bdfode_mag % s >= 2_ip) call memory_alloca(mem_modul(1:2,modul), 'BDFODE_Y2', 'mag_memall', bdfode_mag % y2, meshe(ndivi) % nedge)
    if (bdfode_mag % s >= 3_ip) call memory_alloca(mem_modul(1:2,modul), 'BDFODE_Y3', 'mag_memall', bdfode_mag % y3, meshe(ndivi) % nedge)
    if (bdfode_mag % s >= 4_ip) call memory_alloca(mem_modul(1:2,modul), 'BDFODE_Y4', 'mag_memall', bdfode_mag % y4, meshe(ndivi) % nedge)
    if (bdfode_mag % s >= 5_ip) call memory_alloca(mem_modul(1:2,modul), 'BDFODE_Y5', 'mag_memall', bdfode_mag % y5, meshe(ndivi) % nedge)
    if (bdfode_mag % s >= 6_ip) call memory_alloca(mem_modul(1:2,modul), 'BDFODE_Y6', 'mag_memall', bdfode_mag % y6, meshe(ndivi) % nedge)
    call memory_alloca(mem_modul(1:2,modul), 'BDFODE_YPI', 'mag_memall', bdfode_mag % ypi, meshe(ndivi) % nedge)
    call memory_alloca(mem_modul(1:2,modul), 'BDFODE_YPE', 'mag_memall', bdfode_mag % ype, meshe(ndivi) % nedge)
    !
    !------------------------------------------------------------------------------
  else
    call memory_alloca(mem_modul(1:2,modul), 'HE', 'mag_memall', he_mag, 1_ip)

    call memory_alloca(mem_modul(1:2,modul), 'HA', 'mag_memall', ha_mag, 1_ip)

    call memory_alloca(mem_modul(1:2,modul), 'HP', 'mag_memall', hp_mag, 1_ip)

    call memory_alloca(mem_modul(1:2,modul), 'DH', 'mag_memall', dh_mag, 1_ip)

    call memory_alloca(mem_modul(1:2,modul), 'BHSID', 'mag_memall', bhsid, 1_ip)

    call memory_alloca(mem_modul(1:2,modul), 'FHSID', 'mag_memall', fhsid, 1_ip)

    call memory_alloca(mem_modul(1:2,modul), 'HHSID', 'mag_memall', hhsid, 1_ip)

    do iconstr = 1, constr_total
      call memory_alloca(mem_modul(1:2,modul), 'CONSTR_CNSTR', 'mag_memall', constr_mag(iconstr) % cnstr, 1_ip)
      call memory_alloca(mem_modul(1:2,modul), 'CONSTR_UNKNO', 'mag_memall', constr_mag(iconstr) % unkno, 1_ip)
    end do
    !
    call memory_alloca(mem_modul(1:2,modul), 'HC', 'mag_memall', Hc_mag, 1_ip, 1_ip)
    call memory_alloca(mem_modul(1:2,modul), 'HN', 'mag_memall', Hn_mag, 1_ip, 1_ip)
    call memory_alloca(mem_modul(1:2,modul), 'BC', 'mag_memall', Bc_mag, 1_ip, 1_ip)
    call memory_alloca(mem_modul(1:2,modul), 'BN', 'mag_memall', Bn_mag, 1_ip, 1_ip)
    call memory_alloca(mem_modul(1:2,modul), 'JC', 'mag_memall', Jc_mag, 1_ip, 1_ip)
    call memory_alloca(mem_modul(1:2,modul), 'JN', 'mag_memall', Jn_mag, 1_ip, 1_ip)
    call memory_alloca(mem_modul(1:2,modul), 'JCZ', 'mag_memall', Jcz_mag, 1_ip)
    call memory_alloca(mem_modul(1:2,modul), 'JNZ', 'mag_memall', Jnz_mag, 1_ip)
    call memory_alloca(mem_modul(1:2,modul), 'FC', 'mag_memall', Fc_mag, 1_ip, 1_ip)
    call memory_alloca(mem_modul(1:2,modul), 'FN', 'mag_memall', Fn_mag, 1_ip, 1_ip)
    call memory_alloca(mem_modul(1:2,modul), 'GC', 'mag_memall', Gc_mag, 1_ip)
    call memory_alloca(mem_modul(1:2,modul), 'GN', 'mag_memall', Gn_mag, 1_ip)
    !
    ! BDF
    !
    if (bdfode_mag % s >= 1_ip) call memory_alloca(mem_modul(1:2,modul), 'BDFODE_Y1', 'mag_memall', bdfode_mag % y1, 1_ip)
    if (bdfode_mag % s >= 2_ip) call memory_alloca(mem_modul(1:2,modul), 'BDFODE_Y2', 'mag_memall', bdfode_mag % y2, 1_ip)
    if (bdfode_mag % s >= 3_ip) call memory_alloca(mem_modul(1:2,modul), 'BDFODE_Y3', 'mag_memall', bdfode_mag % y3, 1_ip)
    if (bdfode_mag % s >= 4_ip) call memory_alloca(mem_modul(1:2,modul), 'BDFODE_Y4', 'mag_memall', bdfode_mag % y4, 1_ip)
    if (bdfode_mag % s >= 5_ip) call memory_alloca(mem_modul(1:2,modul), 'BDFODE_Y5', 'mag_memall', bdfode_mag % y5, 1_ip)
    if (bdfode_mag % s >= 6_ip) call memory_alloca(mem_modul(1:2,modul), 'BDFODE_Y6', 'mag_memall', bdfode_mag % y6, 1_ip) 
    call memory_alloca(mem_modul(1:2,modul), 'BDFODE_YPI', 'mag_memall', bdfode_mag % ypi, 1_ip)
    call memory_alloca(mem_modul(1:2,modul), 'BDFODE_YPE', 'mag_memall', bdfode_mag % ype, 1_ip)
    !
  end if
  !
end subroutine mag_memall
