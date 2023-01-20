!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_der_aux_all(&
     ielem,pgaus,pnode,gpsha,gpcar,dummi,dsou_dtem,dsou_dcon,gpdde)
  !-----------------------------------------------------------------------
  !****f* Temper/tem_der_aux_all
  ! NAME
  !   tem_der_aux_all
  ! DESCRIPTION
  !    calculate the derivatives of the some variables w.r.t. (u,p,T,Y)
  !    until now 
  !              d( SUM[hk*Wk] + SUM[Cpk*Difk*dYk/dxi]*dT/dxi )/dT
  !              d( SUM[hk*Wk] + SUM[Cpk*Difk*dYk/dxi]*dT/dxi )/dYs 
  ! USES
  ! USED BY
  !    tem_elmope_new 
  !***
  !-----------------------------------------------------------------------
  use mod_ker_proper
  use def_kintyp, only       :  ip,rp
  use def_master
  use def_temper, only       :  kfl_regim_tem
  use def_domain, only       :  ndime,mnode

  implicit none 
  integer(ip), intent(in)    :: ielem,pgaus,pnode
  integer(ip), intent(in)    :: dummi
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(inout) :: gpdde(pnode,pgaus)                    ! Density derivatives w.r.t nodal temperature
  real(rp),    intent(inout) :: dsou_dtem(pnode,pgaus)                ! d(hW)/dT
  real(rp),    intent(inout) :: dsou_dcon(nspec,pnode,pgaus)          ! d(hW)/dYk
  integer(ip)                :: igaus,inode,ispec
  
  
  !
  ! low mach regime: calculate d(rho)/dT 
  !
  if ( kfl_regim_tem==3 ) then
    call ker_proper('DRDEN','PGAUS',dummi,ielem,gpdde,pnode,pgaus,gpsha,gpcar)
  endif
  !
  ! coupling with chemic: calculate d(wh)/dT and d(wh)/dYk
  !
  if (kfl_coupl(ID_TEMPER,ID_CHEMIC) >= 1 ) then
      ! d( SUM[hk*Wk] )/dT ... the derivatives of the remaining term SUM[Cpk*Difk*dYk/dxi]*dT/dxi w.r.t. T is almost zero
      do inode=1,pnode
	  do igaus=1,pgaus
	      dsou_dtem(inode,igaus) = dchemicalHeat_dtem(ielem)%a(inode,igaus) !+ ddiv_enthalpy_transport_dtem(ielem)%a(inode,igaus)
	  end do
      enddo
      ! d( SUM[hk*Wk] )/dYs ... the derivatives of the remaining term SUM[Cpk*Difk*dYk/dxi]*dT/dxi w.r.t Y is almost zero
      do ispec = 1, nspec
	  do inode=1,pnode
	      do igaus=1,pgaus
		  dsou_dcon(ispec,inode,igaus) = dchemicalHeat(ielem)%a(inode,igaus,ispec)
	      end do
	  enddo
      enddo
  endif
  
  
end subroutine tem_der_aux_all
