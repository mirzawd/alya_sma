!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_updeda()
  !------------------------------------------------------------------------
  !****f* Turbul/tur_updeda
  ! NAME 
  !    tur_updeda
  ! DESCRIPTION
  !    This subroutine updates turb viscosity at gauss points
  !    Only valid for standard K eps for the moment
  ! USES
  !    ker_proper for density
  ! USED BY
  !    tur_updeda
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_turbul  , only     :  turvi_tur, nturb_tur, param_tur, TUR_K_EPS_STD
  use def_kermod
  use mod_ker_proper
  

  implicit none
  real(rp)      :: gpsha(mnode,mgaus), gpmut , gptur( nturb_tur)                  
  real(rp)      :: cmu, relat, gpden(1)
  integer(ip)   :: ielem, pnode, pgaus, pelty, ipoin
  integer(ip)   :: igaus,inode

  if (.not.TUR_K_EPS_STD) return
  
  relat =0.5_rp ! relaxation for turbulent viscosity
  Cmu  = param_tur(6)

  elements: do ielem = 1,nelem
     !
     ! Initialize
     !
     pelty = ltype(ielem)
     pnode = nnode(pelty)
     pgaus = ngaus(pelty)

     do igaus = 1,pgaus         
        gptur = 0.0_rp
        do inode = 1,pnode
           ipoin = lnods(inode, ielem)
        
           gpsha(inode,igaus) =  elmar(pelty)%shape(inode,igaus)
           gptur(1:nturb_tur) =  gptur(1:nturb_tur) &
                + gpsha(inode, igaus)*untur(1:nturb_tur, ipoin, 1)      
        end do
        
        call ker_proper('DENSI','IGAUS',igaus,ielem,gpden,pnode,pgaus,gpsha)
        
        if (kfl_logva==0) then
           gpmut = gpden(1)*cmu * gptur(1)*gptur(1)/gptur(2)
        else  if (kfl_logva==1) then
           gpmut =  gpden(1)*cmu * exp(2.0_rp*gptur(1)- gptur(2))
        end if
        turvi_tur(1,igaus, ielem) = max(relat*gpmut + (1.0_rp-relat)*turvi_tur(1,igaus, ielem), 5.0e-5_rp)
        turvi_tur(1,igaus, ielem) = min(  turvi_tur(1,igaus, ielem), 1.0e12_rp)
     end do

  end do elements
  turvi_tur(2,1:mgaus, 1:nelem)= turvi_tur(1,1:mgaus, 1:nelem)
end subroutine tur_updeda
