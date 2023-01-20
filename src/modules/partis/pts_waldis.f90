!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @file    pts_waldis.f90
!> @author  Guillaume Houzeaux
!> @date    16/11/1966   
!> @brief   Compute some useful variables
!> @details Variables computed:
!>          HLENG_PTS:               Element length                     -
!>          LBOUE_PTS:               Boundary condition on boundaries   -
!>          LELEBOUN_PTS(IELEM) % L: List of boundariesof element IELEM x
!>          BOUNO_PTS:               Normals at boudnaries              -
!>          WALLD_SLIP_PTS:          Distance to slip boundary          -
!>          FRICTION_PTS:            Friction on nodes                  -
!>          
!> @} 
!-----------------------------------------------------------------------
subroutine pts_waldis()
  use def_kintyp,        only : ip,rp,i1p
  use def_master
  use def_domain
  use def_kermod
  use def_partis
  use mod_kdtree
  use mod_memory
  use mod_elmgeo,         only : elmgeo_element_characteristic_length
  use mod_maths,          only : maths_heap_sort
  use mod_graphs,         only : graphs_elepoi
  use mod_graphs,         only : graphs_dealep
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications, only : PAR_GHOST_BOUNDARY_EXCHANGE
  use mod_communications, only : PAR_SUM,PAR_MAX
  use mod_ADR,            only : ADR_assemble_laplacian
  use mod_solver,         only : solver_solve
  use mod_solver,         only : solver_initialize_matrix_and_rhs
  use mod_gradie,         only : gradie
  use mod_projec,         only : projec_boundaries_to_nodes
  use mod_elmgeo,         only : element_type
  use mod_outfor,         only : outfor
  use mod_strings,        only : integer_to_string 
  use mod_pts_arrays
  use mod_bouder
  implicit none
  integer(ip)             :: iboun,pblty,pnodb,jnode,itype,knodb
  integer(ip)             :: ielem,inode,ipoin,ierro,pelty,ptopo
  integer(ip)             :: jelem,iz,kboun,pnode,jpoin,iedge
  integer(ip)             :: lelem_face(100),inodb
  real(rp)                :: baloc(ndime,ndime),eucta,dimax
  real(rp)                :: bocod(ndime,mnodb),fact1,dista
  real(rp)                :: elcod(ndime,mnode),fact2,hleng(3)
  type(i1p),    pointer   :: lelem_boun_2(:)
  integer(ip),  pointer   :: nelem_boun_2(:)
  logical(lg),  pointer   :: boundary_mask(:)
  
  nullify(nelem_boun_2)
  nullify(lelem_boun_2)
  nullify(boundary_mask)
  !
  ! Redimension arrays to account for halos
  !
  !write(90+kfl_paral,*) mem_modul(1:2,modul)
  call memory_resize(mem_modul(1:2,modul),'KFL_FIXBO_PTS','pts_inibcs',kfl_fixbo_pts,nboun_2)
  call memory_resize(mem_modul(1:2,modul),'BVNAT_PTS'    ,'pts_inibcs',bvnat_pts,1_ip,nboun_2)

  call pts_arrays('RESIZE')
  !call memory_resize(mem_modul(1:2,modul),'MOMEN'        ,'pts_inibcs',momentum_sink,ndime,npoin_2)
  !call memory_resize(mem_modul(1:2,modul),'HEATS'        ,'pts_inibcs',heat_sink,npoin_2)
  !call memory_resize(mem_modul(1:2,modul),'MASSK'        ,'pts_inibcs',mass_sink,npoin_2)
  !call memory_resize(mem_modul(1:2,modul),'RESID'        ,'pts_inibcs',resid_pts,ntyla_pts,nelem_2)
  !call memory_resize(mem_modul(1:2,modul),'DEPOE'        ,'pts_inibcs',depoe_pts,ntyla_pts,nelem_2)
  !call memory_resize(mem_modul(1:2,modul),'DEPOB'        ,'pts_inibcs',depob_pts,ntyla_pts,nboun_2)
    
  if( INOTMASTER ) then
     !
     ! ELement charcateristic minimum length
     !
     do ielem = 1,nelem
        pelty = abs(ltype(ielem))
        pnode = nnode(pelty)
        ptopo = ltopo(pelty)
        do inode = 1,pnode
           ipoin = lnods(inode,ielem)
           elcod(1:ndime,inode) = coord(1:ndime,ipoin)
        end do
        call elmgeo_element_characteristic_length(& 
             ndime,pnode,elmar(pelty) % dercg,elcod,hleng)
        hleng_pts(ielem) = hleng(ndime)
     end do
     !
     ! Mark node in wall elements
     !
     do iboun = 1,nboun
        ielem = lelbo(iboun)
        lboue_pts(ielem) = 1
     end do
     do iboun = 1,nboun
        if( kfl_fixbo_pts(iboun) == PTS_SLIP_CONDITION ) then
           ielem = lelbo(iboun)
           lboue_pts(ielem) = PTS_SLIP_CONDITION
        else if( kfl_fixbo_pts(iboun) == PTS_BOUNCING_CONDITION ) then
           ielem = lelbo(iboun)
           lboue_pts(ielem) = PTS_BOUNCING_CONDITION
        end if
     end do
     !
     ! List of boundaries connected to boundary elements (LBOUE_PTS)
     !    
     call memory_alloca(mem_modul(1:2,modul),'LELEBOUN_PTS','pts_waldis',leleboun_pts,nelem)     
     call memory_alloca(mem_modul(1:2,modul),'NELEM_BOUN_2','pts_waldis',nelem_boun_2,nelem_2)
     call memory_alloca(mem_modul(1:2,modul),'LELEM_BOUN_2','pts_waldis',lelem_boun_2,nelem_2)
     !
     ! Boundaries connected to neigboring elements using LBOEL
     !
     do iboun = 1,nboun_2
        ielem = lelbo(iboun)
        nelem_boun_2(ielem) = nelem_boun_2(ielem) + 1
     end do
     do ielem = 1,nelem_2
        call memory_alloca(mem_modul(1:2,modul),'LELEM_BOUN_2 % L','pts_waldis',lelem_boun_2(ielem) % l,nelem_boun_2(ielem))
        nelem_boun_2(ielem) = 0
     end do
     do iboun = 1,nboun_2
        ielem = lelbo(iboun)
        nelem_boun_2(ielem) = nelem_boun_2(ielem) + 1
        lelem_boun_2(ielem) % l(nelem_boun_2(ielem)) = iboun 
     end do
     
     do ielem = 1,nelem
        if( lboue_pts(ielem) > 0 ) then
           iboun = 0
           iz    = pelel_2(ielem)-1
           do while( iz <= pelel_2(ielem+1)-1 )
              if( iz == pelel_2(ielem)-1 ) then
                 jelem = ielem
              else
                 jelem = lelel_2(iz)          
              end if
              do kboun = 1,nelem_boun_2(jelem)
                 iboun             = iboun + 1
                 lelem_face(iboun) = lelem_boun_2(jelem) % l(kboun)
              end do
              iz = iz + 1
           end do
           if( iboun > 0 ) then
              call maths_heap_sort(2_ip,iboun,lelem_face,ELIMINATE_REPLICATES=.true.)
              call memory_alloca(mem_modul(1:2,modul),'LELEBOUN_PTS % L','pts_waldis',leleboun_pts(ielem) % l,iboun)
              do kboun = 1,iboun
                 leleboun_pts(ielem) % l(kboun) = lelem_face(kboun)
              end do
           end if
        end if
     end do

     call memory_deallo(mem_modul(1:2,modul),'NELEM_BOUN_2','pts_waldis',nelem_boun_2)
     call memory_deallo(mem_modul(1:2,modul),'LELEM_BOUN_2','pts_waldis',lelem_boun_2)
     !
     ! Boundary normals
     !
     call memory_alloca(mem_modul(1:2,modul),'BOUNO_PTS','pts_waldis',bouno_pts,ndime,nboun_2)               
     do iboun = 1,nboun
        pnodb = lnnob(iboun)
        pblty = abs(ltypb(iboun))
        ielem = lelbo(iboun)
        pnode = lnnod(ielem)
        bocod(1:ndime,1:pnodb) = coord(1:ndime,lnodb(1:pnodb,iboun))
        elcod(1:ndime,1:pnode) = coord(1:ndime,lnods(1:pnode,ielem))
        call bouder(&
             pnodb,ndime,ndimb,elmar(pblty) % dercg,&
             bocod,baloc,eucta)
        call chenor(pnode,baloc,bocod,elcod) 
        bouno_pts(1:ndime,iboun) = baloc(1:ndime,ndime)
     end do
     !
     ! Exchange boundary codes and normals
     !
     call PAR_GHOST_BOUNDARY_EXCHANGE(kfl_fixbo_pts,'SUBSTITUTE','IN MY CODE')
     call PAR_GHOST_BOUNDARY_EXCHANGE(bvnat_pts    ,'SUBSTITUTE','IN MY CODE')
     call PAR_GHOST_BOUNDARY_EXCHANGE(bouno_pts,    'SUBSTITUTE','IN MY CODE')

  end if
  
  !----------------------------------------------------------------------
  !
  ! Slip wall distance
  !
  !    Compute the generalized distance to the wall via a Poisson equation:
  !    1. Solve Lapl(f) = -1, with f = 0 on wall
  !    2. d = sqrt[ grad(f)^2 +2*f ] - sqrt[grad(f)^2]
  !    See the following references:
  !    P.G. Tucker, Differential equation-based wall distance computation for
  !         DES and RANS, J. Comp. Phys. 190 (2003) 229-248.
  !    P.G. Tucker, Int. J. Numer. Fluids 33 (2000) 869.
  !    P.G. Tucker, Appl. Math. Model. 22 (1998) 293.
  !
  !----------------------------------------------------------------------

  solve_sol => solve(1:)
  if( kfl_slip_wall_pts > 0 ) then
     do ipoin = 1,npoin
        walld_slip_pts(ipoin) =  1.0_rp
     end do
     call solver_initialize_matrix_and_rhs(solve_sol,amatr,rhsid)     
     call ADR_assemble_laplacian(meshe(ndivi),elmar,amatr,walld_slip_pts,rhsid)
     do ipoin = 1,npoin
        walld_slip_pts(ipoin) =  0.0_rp
        unkno(ipoin)          =  0.0_rp
     end do
     call solver_solve(solve_sol,amatr,rhsid,unkno)    
     if( INOTMASTER ) then
        call memgen(0_ip,ndime,npoin)
        call gradie(unkno,gevec)
        do ipoin = 1,npoin
           fact1 = dot_product(gevec(1:ndime,ipoin),gevec(1:ndime,ipoin))
           fact2 = fact1 + 2.0_rp * max(unkno(ipoin),0.0_rp)
           walld_slip_pts(ipoin) = sqrt(fact2) - sqrt(fact1)  
        end do
        call memgen(2_ip,ndime,npoin)
     end if
  end if

  solve_sol => solve(2:)
  if( kfl_bouncing_wall_pts > 0 ) then
     do ipoin = 1,npoin
        walld_bouncing_pts(ipoin) =  1.0_rp
     end do
     call solver_initialize_matrix_and_rhs(solve_sol,amatr,rhsid)     
     call ADR_assemble_laplacian(meshe(ndivi),elmar,amatr,walld_bouncing_pts,rhsid)
     do ipoin = 1,npoin
        walld_bouncing_pts(ipoin) =  0.0_rp
        unkno(ipoin)              =  0.0_rp
     end do
     call solver_solve(solve_sol,amatr,rhsid,unkno)    
     if( INOTMASTER ) then
        call memgen(0_ip,ndime,npoin)
        call gradie(unkno,gevec)
        do ipoin = 1,npoin
           fact1 = dot_product(gevec(1:ndime,ipoin),gevec(1:ndime,ipoin))
           fact2 = fact1 + 2.0_rp * max(unkno(ipoin),0.0_rp)
           walld_bouncing_pts(ipoin) = sqrt(fact2) - sqrt(fact1)  
        end do
        call memgen(2_ip,ndime,npoin)
     end if
  end if

  !----------------------------------------------------------------------
  !
  ! Compute friction, use same solver as previous one
  !
  !----------------------------------------------------------------------

  solve_sol => solve(1:)
  if( kfl_slip_wall_pts > 0 ) then
     if( INOTMASTER ) then
        call memory_alloca(mem_modul(1:2,modul),'BOUNDARY_MASK','pts_waldis',boundary_mask,nboun)
        call solver_initialize_matrix_and_rhs(solve_sol,amatr,rhsid)             
        do iboun = 1,nboun
           if( kfl_fixbo_pts(iboun) == PTS_SLIP_CONDITION ) boundary_mask(iboun) = .true.        
        end do
        call projec_boundaries_to_nodes(bvnat_pts,meshe(ndivi),elmar,friction_pts,boundary_mask)   
        do ipoin = 1,npoin
           rhsid(ipoin) = 0.0_rp
        end do
        call ADR_assemble_laplacian(meshe(ndivi),elmar,amatr) 
     end if
     call solver_solve(momod(modul) % solve,amatr,rhsid,friction_pts)     
     call memory_deallo(mem_modul(1:2,modul),'BOUNDARY_MASK','pts_waldis',boundary_mask)
  end if

  !----------------------------------------------------------------------
  !
  ! Check smallest distance and hmin
  !
  !----------------------------------------------------------------------

  if( kfl_walld_pts == 0 ) then
     !
     ! Maximum particle diameter
     !
     dimax = 0.0_rp
     do itype = 1,ntyla_pts
        if( parttyp(itype) % kfl_exist /= 0 ) then
           dimax = max(dimax,parttyp(itype) % diame)
        end if
     end do

     ierro = 0
     if( INOTMASTER ) then
        iboun_loop: do iboun = 1,nboun
           ielem = lelbo(iboun)
           pelty = abs(ltype(ielem))
           do iedge = 1,element_type(pelty) % number_edges 
              inode = element_type(pelty) % list_edges(1,iedge) 
              jnode = element_type(pelty) % list_edges(2,iedge)
              ipoin = lnods(inode,ielem)
              jpoin = lnods(jnode,ielem)
              dista = sqrt(dot_product(coord(:,ipoin)-coord(:,jpoin),coord(:,ipoin)-coord(:,jpoin)))
              if( dimax > 0.9_rp * dista ) then
                 ierro = 1
                 exit iboun_loop
              end if
           end do
        end do iboun_loop
     end if

     call PAR_SUM(ierro)
     if( ierro > 0 ) call outfor(2_ip,momod(modul)%lun_outpu,&
          'SOME PARTICLES ARE LARGER THAN SMALLEST WALL ELEMENT: ACTIVATE WALLD: ON OPTION')

  else
     !
     ! Check WALLD has been calculated
     !
     if( kfl_walld == 0 ) then
        call runend('PTS_WALDIS: ACTIVATE WALLD IN KERMOD')
     end if
     
  end if
  
  !----------------------------------------------------------------------
  !
  ! Check if slip walls are not no-slip walls.
  ! At this stage, kfl_fixbo_walld_ker is not available
  !
  !----------------------------------------------------------------------

  if( kfl_walld_pts /= 0 .and. kfl_slip_wall_pts /= 0 ) then
     ierro = 0
     kboun = 0
     do iboun = 1,nboun
        if( kfl_fixbo_pts(iboun) == PTS_SLIP_CONDITION ) then
           pnodb = element_type(ltypb(iboun)) % number_nodes
           knodb = 0
           do inodb = 1,pnodb
              ipoin = lnodb(inodb,iboun)
              if( kfl_fixno_walld_ker(1,ipoin) == 1 ) knodb = knodb + 1
           end do
           if( knodb == pnodb ) then
              ierro = ierro + 1
              kboun = max(kboun,lbinv_loc(iboun))
           end if
        end if
     end do
     call PAR_SUM(ierro)
     if( ierro /= 0 ) then
        call PAR_MAX(kboun)
        call runend('PTS_WALDIS: SOME WALLS ARE BOTH SLIP AND NO-SLIP WALL. E.G., CHECK BOUNDARY '//integer_to_string(kboun))
     end if
  end if
  
  call memory_deallo(mem_modul(1:2,modul),'KFL_FIXNO_WALLD_SLIP_PTS','pts_waldis',    kfl_fixno_walld_slip_pts)
  call memory_deallo(mem_modul(1:2,modul),'KFL_FIXNO_WALLD_BOUNCING_PTS','pts_waldis',kfl_fixno_walld_bouncing_pts)

end subroutine pts_waldis
