!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_zero_divergence_velocity.f90
!> @author  Guillaume Houzeaux
!> @brief   Zero divergence velocity
!> @details Correct initial velocity so that div(u)=0
!> @} 
!------------------------------------------------------------------------

subroutine nsi_zero_divergence_velocity()

  use def_parame 
  use def_master
  use def_elmtyp
  use def_kermod
  use def_domain
  use def_nastin 
  use mod_nsi_schur_operations
  use mod_matrix,         only : matrix_assemble_element_matrix_to_CSR
  use mod_matrix,         only : matrix_assemble_element_RHS
  use mod_communications, only : PAR_MAX
  use mod_gradie,         only : gradie
  use mod_solver,         only : solver_solve
  use mod_solver,         only : solver_initialize_matrix_and_rhs
  use mod_local_basis,    only : local_basis_global_to_local
  use mod_local_basis,    only : local_basis_local_to_global
  implicit none

  integer(ip) :: ielem,pelty,pnode,idime
  integer(ip) :: pgaus,plapl,ipoin,ibopo
  real(rp)    :: elmat(mnode,mnode)
  real(rp)    :: elrhs(mnode) 
  real(rp)    :: elvel(ndime,mnode)
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: gpsha(mnode,mgaus)                    ! N
  real(rp)    :: gpder(ndime,mnode,mgaus)              ! dN/dsi
  real(rp)    :: gpcar(ndime,mnode,mgaus)              ! dN/dxi
  real(rp)    :: gphes(ntens,mnode,mgaus)              ! d2N/dxidxj
  real(rp)    :: gpvol(mgaus)                          ! w*|J|, |J|
  real(rp)    :: alpha(3)
  real(rp)    :: udotn
  !
  ! Initialize solver
  !
  solve_sol => solve(6:)
  call solver_initialize_matrix_and_rhs(solve_sol,amatr,rhsid)
  !solve(6) % kfl_assem = SOL_MATRIX_HAS_CHANGED
  !
  ! Adaptive boundary conditions
  !
  do ipoin = 1,npoin
     if( kfl_fixno_div_nsi(1,ipoin) == 8 ) then
        ibopo = lpoty(ipoin)
        if( ibopo /= 0 ) then
           udotn = dot_product(veloc(1:ndime,ipoin,nprev_nsi),exnor(1:ndime,1,ibopo))
           if( udotn <= 0.0_rp ) kfl_fixno_div_nsi(1,ipoin) = 0
        end if
     end if
  end do
  !
  ! Assemble system A,b
  !
  if( INOTMASTER ) then

     unkno(1:npoin) = 0.0_rp

     elements: do ielem = 1,nelem
        pelty = ltype(ielem)
        if( pelty > 0 ) then
           pnode = nnode(pelty)
           pgaus = ngaus(pelty)
           plapl = 0
           !
           ! Gather
           ! 
           elvel(1:ndime,1:pnode) = veloc(1:ndime,lnods(1:pnode,ielem),nprev_nsi)
           elcod(1:ndime,1:pnode) = coord(1:ndime,lnods(1:pnode,ielem))
           !
           ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
           !
           call elmca2(&
                pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
                elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpsha,&
                gpder,gpcar,gphes,ielem)
           ! 
           ! Assemble matrix
           !
           call nsi_zero_divergence_velocity_assembly(&
                ndime,pnode,pgaus,gpvol,gpsha,gpcar,elvel,elmat,elrhs)
           !
           ! Assembly
           !
           call matrix_assemble_element_RHS(&
                1_ip,1_ip,pnode,lnods(:,ielem),elrhs,rhsid)
           call matrix_assemble_element_matrix_to_CSR(&
                kfl_element_to_csr,1_ip,pnode,pnode,&
                ielem,lnods(:,ielem),elmat,r_dom,c_dom,amatr,lezdo)
        end if
     end do elements 
           
  end if 
  !
  ! Solver problem Ax=b
  !   
  call solver_solve(solve_sol,amatr,rhsid,unkno)
  !
  ! Correction velocity
  !
  if( INOTMASTER ) then
     alpha(1:2) = 0.5_rp
     alpha(3)   = 0.5_rp * divcorrec_nsi**2
     call memgen(0_ip,ndime,npoin) 
     call gradie(unkno,gevec)     
     !
     ! Correction: u = u + alpha^2/2 grad(phi)
     !
     do ipoin = 1,npoin
        gevec(1:ndime,ipoin) = gevec(1:ndime,ipoin) * alpha(1:ndime)
     end do
     !
     ! Boundary conditions: impose grad(phi)=0 manually where velocity is prescribed
     !
     call local_basis_global_to_local(kfl_fixrs_nsi,gevec) ! Gobal to local
     do ipoin = 1,npoin 
        do idime = 1,ndime
           if( kfl_fixno_nsi(idime,ipoin) > 0 ) gevec(idime,ipoin) = 0.0_rp
        end do
     end do
     call local_basis_local_to_global(kfl_fixrs_nsi,gevec) ! Gobal to local
     do ipoin = 1,npoin
        veloc(1:ndime,ipoin,nprev_nsi) = veloc(1:ndime,ipoin,nprev_nsi) + gevec(1:ndime,ipoin)
     end do
     call memgen(2_ip,ndime,npoin)     
  end if

end subroutine nsi_zero_divergence_velocity

!----------------------------------------------------------------------
!
!> @author  Guillaume Houzeaux
!> @date    13/10/2014
!> @brief   Impose div(u)=0 to a given velocity field u0
!> @details Solve the following equation:
!>          div(k*grad(phi)) = -2 div(u0) 
!>          with k_ij = delta_ij (1,1,alpha^2)
!>          and then compute u = u0 + alpha^2/2 grad(phi)
!
!----------------------------------------------------------------------

subroutine nsi_zero_divergence_velocity_assembly(&
     ndime,pnode,pgaus,gpvol,gpsha,gpcar,elvel,elmat,elrhs)

  use def_kintyp, only : ip,rp
  use def_domain, only : mnode
  use def_nastin, only : divcorrec_nsi
  implicit none
  integer(ip), intent(in)  :: ndime
  integer(ip), intent(in)  :: pnode
  integer(ip), intent(in)  :: pgaus
  real(rp),    intent(in)  :: gpvol(pgaus)
  real(rp),    intent(in)  :: gpsha(pnode,pgaus)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)  :: elvel(ndime,pnode)
  real(rp),    intent(out) :: elmat(pnode,pnode)
  real(rp),    intent(out) :: elrhs(pnode)

  integer(ip)              :: inode,jnode,igaus
  real(rp)                 :: divu,gpdif(3)

  elmat      = 0.0_rp
  elrhs      = 0.0_rp
  gpdif(1:2) = 1.0_rp
  gpdif(3)   = divcorrec_nsi**2

  do igaus = 1,pgaus
     do jnode = 1,pnode
        divu = 0.0_rp
        do inode = 1,pnode
           elmat(inode,jnode) = elmat(inode,jnode) &
                + gpvol(igaus) * dot_product(gpdif(1:ndime)*gpcar(1:ndime,inode,igaus),gpcar(1:ndime,jnode,igaus))
           divu = divu + dot_product(gpcar(1:ndime,inode,igaus),elvel(1:ndime,inode))
        end do
        elrhs(jnode) = elrhs(jnode) + 2.0_rp * gpvol(igaus) * gpsha(jnode,igaus) * divu
     end do

  end do

end subroutine nsi_zero_divergence_velocity_assembly
