!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup NastinInput
!> @{
!> @file    nsi_elmresdiff_all.f90
!> @author  Guillaume Houzeaux
!> @brief   Compute the derivation of GPRES r. t. design variables
!
!-------------------------------------------------------------------
!        +-           -+    +-                      -+   +- -+   +-  -+
!        | elmresudiff  |   | elauudiff    elaupdiff |   | u |   | elrbudiff |
!        |              | = |                        | * |   | - |           |
!        | elmresbdiff  |   | elapudiff    elappdiff |   | p |   | elrbpdiff |
!        +-           -+    +-                      -+   +- -+   +-  -+
!-------------------------------------------------------------------
!
!> @} 
!-----------------------------------------------------------------------
subroutine nsi_elmresdiff_all(pnode,pgaus,pevat,gpvol,gpvol_der,gpsha,lnods, &
                                    p1vec,p2vec,p2sca,elvel,elpre,wgrgr,h,h_der, &
                                    gpsp1,gpsp2,rmomu,gpden,gpadv,gpvis,gpcar,gpcar_der,rcont,sens_mesh,ielem)

  use def_kintyp, only     :  ip,rp
  use def_kermod, only     :  kfl_dvar_type
!  use def_kermod, only     :  kfl_ndvars_opt
  use def_domain, only     :  ndime,mnode
  use def_parame, only     :  pi
  use def_master

  
  implicit none
  integer(ip), intent(in)    :: pnode,pgaus,pevat,ielem
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(in)    :: p1vec(pnode,pgaus)
  real(rp),    intent(in)    :: elvel(ndime, pnode,*)
  real(rp),    intent(in)    :: elpre(pnode)
  real(rp),    intent(in)    :: p2vec(ndime,pnode,pgaus)
  real(rp),    intent(in)    :: p2sca(pnode,pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus),gpsha(pnode,pgaus),gpvol_der(ndime,pnode,pgaus)
  real(rp),    intent(in)    :: wgrgr(pnode,pnode,pgaus)
  real(rp),    intent(in)    :: h
  real(rp),    intent(in)    :: h_der(ndime,pnode)
  real(rp),    intent(in)    :: gpsp1(pgaus)                      ! tau1'
  real(rp),    intent(in)    :: gpsp2(pgaus)                      ! tau2'
  real(rp),    intent(in)    :: rmomu(pnode,pgaus)
  real(rp),    intent(in)    :: gpden(pgaus)
  real(rp),    intent(in)    :: gpadv(ndime,pgaus)
  real(rp),    intent(in)    :: gpvis(pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gpcar_der(ndime,mnode,ndime,mnode,pgaus)
  real(rp),    intent(in)    :: rcont(ndime,pnode,pgaus)
!   real(rp),    intent(inout) :: resdiff_nsi(kfl_ndvars_opt,*)
  real(rp),    intent(inout) :: sens_mesh(ndime,*)

  integer(ip)                :: idime,inode,ipoin
  real(rp)	             :: elresudiff(ndime,pnode)
  real(rp)	             :: elresucoodiff(ndime,pnode,ndime,pnode)
  real(rp)	             :: elrespdiff(pnode)
  real(rp)	             :: elrespcoodiff(pnode,ndime,pnode)
  
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!        design variables will be viscousity                 !!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (kfl_dvar_type == 2) then 
  
    call nsi_elmresdiff_visco(pnode,pgaus,pevat,gpvol,gpsha,lnods, &
			  p1vec,p2vec,elvel,wgrgr,h, &
			  gpsp1,gpsp2,rmomu,gpden,gpadv,gpcar,rcont,elresudiff,elrespdiff)

    ! 
    !    Effect of Dirichlet boundary conditions (elresudiff = 0.0)
    !
    if( solve(1) % kfl_iffix == 0 ) then
	call nsi_elmdirdiff(pnode,pevat,lnods,elresudiff)
    end if
    !
    ! Assembling for RESDIFF
    !
!     call nsi_assresdiff(1_ip,pnode,pevat,lnods,elresudiff,elrespdiff,resdiff_nsi,1_ip,elvel,elpre)
  
  endif	

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!        design variables will be geometry                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!   calculate dR/dX*Landa and send to alefor RHD where X is the coordinates of the domain points  !!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  if (kfl_dvar_type == 5) then 
  
    call nsi_elmresdiff_coo(pnode,pgaus,pevat,gpvol,gpvol_der,gpsha,lnods, &
			  p1vec,p2vec,p2sca,elvel,wgrgr,h,h_der,&
			  gpsp1,gpsp2,rmomu,gpden,gpadv,gpvis,gpcar,gpcar_der,rcont,elresucoodiff,elrespcoodiff,ielem)

    ! 
    !    Effect of Dirichlet boundary conditions (elresucoodiff = 0.0)
    !
!     if( solve(1) % kfl_iffix == 0 ) then
! 	call nsi_elmdirdiff(pnode,pevat,lnods,elresucoodiff)
!     end if
    !
    ! Assembling for RESDIFF
    !
    do idime = 1, ndime
      do inode = 1, pnode
        ipoin = lnods(inode)
        call nsi_assresdiff(1_ip,pnode,pevat,lnods,elresucoodiff(:,:,idime,inode),elrespcoodiff(:,idime,inode),sens_mesh,idime,ipoin,elvel,elpre)
      enddo
    enddo
  
  endif	
  
  
  
end subroutine nsi_elmresdiff_all


subroutine nsi_elmdirdiff(&
     pnode,pevat,lnods,elrbudiff)
  !------------------------------------------------------------------------
  !****f* Nastin/nsi_elmdirdiff
  ! NAME 
  !    nsi_elmdirdiff
  ! DESCRIPTION
  !    This routine prescribes the Dirichlet conditions to the elresdiff_des
  ! USES
  ! USED BY
  !    nsi_bouope
  !------------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_nastin, only       :  kfl_fixno_nsi
  use def_domain, only       :  ndime
  use def_elmtyp, only       :  ELHAN
  implicit none
  integer(ip), intent(in)    :: pnode,pevat
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(inout) :: elrbudiff(pevat)
  integer(ip)                :: inode,ipoin,idime,ievav
  
  do inode = 1,pnode
     ievav = (inode-1) * ndime
     ipoin = lnods(inode)
     do idime = 1,ndime
        ievav = ievav+1
        if(      kfl_fixno_nsi(idime,ipoin) ==  1 &
             .or.kfl_fixno_nsi(idime,ipoin) ==  8 &
             .or.kfl_fixno_nsi(idime,ipoin) ==  9 &
             .or.kfl_fixno_nsi(idime,ipoin) == 11 ) then
	   elrbudiff(ievav) = 0.0_rp
        end if
     enddo
  end do
  
  
end subroutine nsi_elmdirdiff

