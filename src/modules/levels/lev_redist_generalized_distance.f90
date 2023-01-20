!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Levels
!> @{
!> @file    lev_redist_generalized_distance.f90
!> @author  houzeaux
!> @date    2020-09-07
!> @brief   Redistantiation
!> @details Redistantiation based on the generlized distance
!> @} 
!-----------------------------------------------------------------------

subroutine lev_redist_generalized_distance()

  use def_kintyp_mesh_basic, only : mesh_type_basic
  use def_elmtyp
  use def_master
  use def_domain 
  use def_kermod
  use def_levels
  use mod_elmgeo
  use mod_memory
  use mod_solver
  use mod_elsest,               only : elsest_host_element
  use mod_generalized_distance, only : generalized_distance_assembly
  use mod_generalized_distance, only : generalized_distance_solution
  use mod_generalized_distance, only : generalized_distance_update
  use mod_messages,             only : messages_live
  use mod_communications,       only : PAR_SUM
  use mod_mesh_type_basic,      only : mesh_type_basic_parallel
  use mod_mesh_type_basic,      only : mesh_type_basic_output
  use mod_parall,               only : PAR_COMM_MY_CODE
  use mod_bouder
  
  implicit none

  type(mesh_type_basic) :: mesh
  integer(ip)           :: ipoin,ielem,igaub
  integer(ip)           :: pelty,pnode,pgaub
  integer(ip)           :: inode,jnode,inodb
  integer(ip)           :: iboun,pblty,pnodb
  integer(ip)           :: nelem_sum
  real(rp)              :: baloc(ndime,ndime)
  real(rp)              :: eucta,xx(3)
!  real(rp)              :: err,den
  real(rp)              :: coloc(3),dNds(ndime,mnode)
  real(rp)              :: N(mnode),toler,bocod(ndime,mnode)
  real(rp)              :: elmat(mnode,mnode)
  real(rp),     pointer :: dista(:)
  !
  ! Create mesh
  !
  !do ipoin = 1,npoin
  !   fleve(ipoin,1)=cos(coord(1,ipoin))+coord(2,ipoin)-2.0_rp
  !end do
  dista => fleve(:,1)
  call mesh % init()
  call mesh % cut_level(meshe(ndivi),dista,mem_modul(1:2,modul))
  nelem_sum = mesh % nelem
  call PAR_SUM(nelem_sum)
  
  if( 1 == 2 ) then
     mesh % name = 'LEVEL'
     call mesh_type_basic_parallel(mesh,PAR_COMM_MY_CODE,kfl_paral)
     call mesh_type_basic_output(mesh)
  end if
  
  if( nelem_sum > 0 ) then
     
     call messages_live('REDISTANCIATION USING GENERALIZED DISTANCE')
     !
     ! Assemble
     !
     solve_sol => solve(5:)
     call generalized_distance_assembly(solve_sol,amatr,rhsid)
     !
     ! Boundary condition
     !
     do iboun = 1,mesh % nelem
        pblty = mesh % ltype(iboun)
        pgaub = ngaus(pblty)
        pnodb = nnode(pblty)
        do inodb = 1,pnodb
           ipoin = mesh % lnods(inodb,iboun)
           bocod(1:ndime,inodb) = mesh % coord(1:ndime,ipoin)
        end do
        do igaub = 1,pgaub
           call bouder(&
                pnodb,ndime,ndimb,elmar(pblty)%deriv(:,:,igaub),& 
                bocod,baloc,eucta)               
           eucta = elmar(pblty) % weigp(igaub) * eucta
           xx    = 0.0_rp
           do inodb = 1,pnodb
              xx(1:ndime) = xx(1:ndime) + elmar(pblty) % shape(inodb,igaub) * bocod(1:ndime,inodb)
           end do
           call elsest_host_element(&
                ielse,relse,1_ip,meshe(ndivi),xx,ielem,&
                N,dNds,coloc,toler)
           if( ielem > 0 ) then
              pelty = ltype(ielem)
              pnode = element_type(pelty) % number_nodes
              do inode = 1,pnode
                 do jnode = 1,pnode
                    elmat(inode,jnode) = 10.0_rp**9 * N(inode) * N(jnode) * eucta
                 end do
              end do
              call solver_assemble_element_matrix_scalar(&
                   solve_sol(1),1_ip,pnode,pnode,ielem,lnods(:,ielem),elmat,amatr)
           end if
        end do
     end do
     !
     ! Solve Lu = -1
     !
     do ipoin = 1,npoin
        unkno(ipoin) = abs(fleve(ipoin,1))
     end do
     call generalized_distance_solution(solve_sol,amatr,rhsid,unkno)
     !
     ! Compute distance
     !
     call generalized_distance_update(unkno,dista_lev)

     !call kdtree_construct(&
     !     mesh % nelem,mesh % npoin,mesh % lnods,mesh % ltype,coord,coupling_type(icoup) % geome % kdtree,&
     !     coupling_type(jcoup) % wet % lboun_wet)

     !
     ! Correct level set
     !
     do ipoin = 1,npoin
        if( kfl_fixno_lev(1,ipoin) < 1 ) then
           if( fleve(ipoin,1) >= 0.0_rp ) then
              fleve(ipoin,1) = dista_lev(ipoin)
           else
              fleve(ipoin,1) =-dista_lev(ipoin)
           endif
        end if
     end do

     !err = 0.0_rp
     !den = 0.0_rp
     !do ipoin = 1,npoin
     !   err = err + (fleve(ipoin,1)-(cos(coord(1,ipoin))+coord(2,ipoin)-2.0_rp))**2
     !   den = den + ((cos(coord(1,ipoin))+coord(2,ipoin)-2.0_rp))**2
     !end do
     !call PAR_SUM(err)
     !call PAR_SUM(den)
     !if( INOTSLAVE ) print*,'Error= ',sqrt(err/den)
  end if
  !
  ! Deallocate mesh
  !
  call mesh % deallo()

end subroutine lev_redist_generalized_distance
