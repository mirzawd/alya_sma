!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_resdiff_all(&
     kfl_ortho, kfl_shock,shock_tur, pnode,plapl,pgaus,gpvel,gpdif,&
     gprea,gptur,gpgrd,gprhs,gpden,gpsha,gpcar,gpvol,chale, &
     elunk, gphes, sreac, itask, gppro, gpprr, gppgr,gpmut,lnods,eltur,elwal,gpvis,&
     gprea_diff,gprhs_diff, gpreawal_diff,gprhswal_diff, gpvol_der, gpcar_der,h_der,elwao,sens_mesh,sens_wall,ielem)
  !------------------------------------------------------------------------
  !****f* Turbul/tur_resdiff_all
  ! NAME 
  !    tur_resdiff_all
  ! DESCRIPTION
  !    
  !      This calculates the partial derivative of the residual w.r.t. design variables (dR/dD)
  !      Depondency to the design variables and couplement 
  ! USES
  ! USED BY
  !    
  !------------------------------------------------------------------------

  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,ntens,mnode
  use def_turbul, only       :  nturb_tur
!  use def_turbul, only       :  iunkn_tur,resdiff_tur
  use def_master
  use def_kermod, only       :  kfl_dvar_type,kfl_walld
  !use mod_matrix
  use mod_assresdiff

  
  
  implicit none
  integer(ip), intent(in)    :: pnode,plapl,pgaus, kfl_ortho, kfl_shock, itask,ielem
  real(rp),    intent(in)    :: gpvel(ndime,pgaus), gptur(nturb_tur,3,pgaus)
  real(rp),    intent(in)    :: gpdif(pgaus),gprea(pgaus)
  real(rp),    intent(in)    :: gprhs(pgaus), shock_tur
  real(rp),    intent(in)    :: gpden(pgaus)
  real(rp),    intent(in)    :: gpgrd(ndime,pgaus)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus), chale(2), elunk(pnode, 2),sreac(pgaus)
  real(rp),    intent(in)    :: gphes(ntens,mnode,pgaus)             ! dNk/dxidxj
  real(rp),    intent(in)    :: gppro(pgaus), gpprr(pgaus), gppgr(ndime, pgaus)
  real(rp),    intent(in)    :: gpmut(pgaus)
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(in)    :: eltur(nturb_tur,pnode,*)
  real(rp),    intent(in)    :: elwal(pnode)
  real(rp),    intent(in)    :: gpvis(pgaus)
  real(rp),    intent(inout) :: sens_mesh(ndime,*)
  real(rp),    intent(inout) :: sens_wall(ndime,*)
  integer(ip), intent(in)    :: elwao(pnode)

  real(rp),    intent(in)    :: gpvol_der(ndime,pnode,pgaus),gpcar_der(ndime,mnode,ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gprea_diff(ndime,pnode,pgaus),gprhs_diff(ndime,pnode,pgaus)
  real(rp),    intent(in)    :: gprhswal_diff(ndime,pnode,pgaus),gpreawal_diff(ndime,pnode,pgaus)
  real(rp),    intent(in)    :: h_der(ndime,pnode)
  
  integer(ip)                :: inode,idime,ipoin,ipoinwal
!  real(rp)                   :: elresdiff(pnode)
  real(rp)                   :: elres_diff(pnode,ndime,pnode)
  real(rp)                   :: elreswal_diff(pnode,ndime,pnode)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!        design variables will be constants of turbulent visco   !!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (kfl_dvar_type == 6 .and. itask == 30) then 
!   if (kfl_dvar_type == 6) then 
      
!         call tur_elmmsudiff_tvcon(&
!                    kfl_ortho,kfl_shock,shock_tur,pnode,plapl,pgaus,gpvel,gpdif,&
!                    gprea,gptur,gpgrd,gprhs,gpden,gpsha,gpcar,gpvol,chale,&
!                    elunk,gphes,sreac,itask,gppro,gpprr, gppgr,gpmut,eltur,elwal,gpvis,elresdiff)
!         !
!         ! Effect of Dirichlet boundary conditions (elresdiff = 0.0)
!         !
!         if( solve(iunkn_tur) % kfl_iffix == 0 ) then
! 	    call tur_elmdirdiff(pnode,lnods,elresdiff)
!         end if
!         !
!         ! Assembly
!         !
!         call assresdiff(solve(1) % ndofn, 1_ip, iunkn_tur, pnode, lnods, elresdiff, resdiff_tur)      

        call tur_elmmsudiff_spacon(&
                   kfl_ortho,kfl_shock,shock_tur,pnode,plapl,pgaus,gpvel,gpdif,&
                   gprea,gptur,gpgrd,gprhs,gpden,gpsha,gpcar,gpvol,chale,&
                   elunk,gphes,sreac,gppro,gpprr, gppgr,gpmut,eltur,elwal,gpvis,elres_diff,gprhs_diff,gprea_diff)

        call tur_assresdiff(1_ip,pnode,elres_diff(:,1,1),sens_mesh,1,1,eltur)
                
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!        design variables will be constants of turbulent visco   !!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  elseif (kfl_dvar_type == 5 .and. itask == 30) then 
      
        call tur_elmmsudiff_coo(&
                   kfl_ortho,kfl_shock,shock_tur,pnode,plapl,pgaus,gpvel,gpdif,&
                   gprea,gptur,gpgrd,gprhs,gpden,gpsha,gpcar,gpvol,chale,&
                   elunk,gphes,sreac,gppro,gpprr, gppgr,gpmut,eltur,elwal,gpvis,elres_diff,elreswal_diff,&
                   gprea_diff,gprhs_diff, gpreawal_diff,gprhswal_diff, gpvol_der, gpcar_der,h_der,ielem)
        !
        ! Assembling for RESDIFF
        !
        do idime = 1, ndime
          do inode = 1, pnode
            ipoin = lnods(inode)
            call tur_assresdiff(1_ip,pnode,elres_diff(:,idime,inode),sens_mesh,idime,ipoin,eltur)
            if (kfl_walld == 2 .or. kfl_walld == 3) then
!               ipoinwal = elwai(inode)
              ipoinwal = elwao(inode)
              call tur_assresdiff(1_ip,pnode,elreswal_diff(:,idime,inode),sens_wall,idime,ipoinwal,eltur)  
            endif
          enddo
        enddo
        
  endif
  

end subroutine tur_resdiff_all
