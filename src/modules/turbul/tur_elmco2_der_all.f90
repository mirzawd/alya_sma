!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_elmco2_der_all(&
     itask,pnode,gpvis,gpden,gpsha,gpcar,elvel,elmsh,eltur,eledd,&
     elwal,elust,eltem,elrg2,elsqk,elgrp,elpro,elprd,gpvol,&
     chale,gpdif,gpvel,gpgrv,gptur,&
     hleng,gpcan, elprr, pmate, gpcar_der,&
     drea_dtur,dsou_dtur, &
     dsou_dvel,gpmut,dgpmut,dgpmut_dvel,gprea_diff,gprhs_diff,gpreawal_diff,gprhswal_diff,&
     elcod,elwai_coo)

     
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_elmco2_der_all
  ! NAME
  !   tur_elmco2_der_all
  ! DESCRIPTION
  !    Compute derivative of coefficient of the equation of SST K-W
  ! OUTPUT 
  ! USES
  ! USED BY
  !    tur_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_master, only     :  ID_TURBUL,ID_ALEFOR
  use def_domain, only     :  ndime
  use def_kermod, only     :  kfl_walld
  use def_turbul, only     :  nturb_tur,TUR_SPALART_ALLMARAS,TUR_SST_K_OMEGA
  implicit none
  integer(ip), intent(in)  :: itask,pnode
  real(rp),    intent(in)  :: gpvis,gpden
  real(rp),    intent(in)  :: gpsha(pnode),gpcar(ndime,pnode)
  real(rp),    intent(in)  :: elvel(ndime,pnode)
  real(rp),    intent(in)  :: elmsh(ndime,pnode)
  real(rp),    intent(in)  :: eltur(nturb_tur,pnode)
  real(rp),    intent(in)  :: eledd(pnode),elwal(pnode)
  real(rp),    intent(in)  :: elust(pnode),eltem(pnode)
  real(rp),    intent(in)  :: gpvol,chale(2)
  real(rp),    intent(in)  :: elrg2(pnode),elsqk(pnode)
  real(rp),    intent(in)  :: elgrp(ndime,pnode)
  real(rp),    intent(in)  :: elpro(pnode)
  real(rp),    intent(in)  :: elprd(pnode)
  real(rp),    intent(in)  :: elprr(pnode)
  real(rp),    intent(in)  :: gptur(nturb_tur)
  real(rp),    intent(in)  :: gpgrv(ndime,ndime)
  real(rp),    intent(in)  :: gpdif
  real(rp),    intent(in)  :: gpvel(ndime)
  real(rp),    intent(in)  :: gpmut
  real(rp),    intent(in)  :: dgpmut(nturb_tur,pnode)
  real(rp),    intent(in)  :: dgpmut_dvel(ndime,pnode)
  
  real(rp),    intent(out) :: drea_dtur(nturb_tur,pnode)
  real(rp),    intent(out) :: dsou_dtur(nturb_tur,pnode)
  real(rp),    intent(out) :: dsou_dvel(ndime,pnode)
  real(rp),    intent(out) :: gprea_diff(ndime,pnode)
  real(rp),    intent(out) :: gprhs_diff(ndime,pnode)
  real(rp),    intent(out) :: gpreawal_diff(ndime,pnode)
  real(rp),    intent(out) :: gprhswal_diff(ndime,pnode)
  
  real(rp),    intent(in)  :: hleng(ndime)
  real(rp),    intent(in)  :: gpcan
  integer(ip), intent(in)  :: pmate
  real(rp),    intent(in)  :: gpcar_der(ndime,pnode,ndime,pnode)
  
  real(rp),    intent(in)  :: elcod(ndime,pnode)
  real(rp),    intent(in)  :: elwai_coo(ndime,pnode)
  
  
  integer(ip)              :: idime,inode
  real(rp)                 :: gpwal,gpwal_der(ndime,pnode)
    
  gpwal = 0.0_rp
  do inode = 1,pnode
     gpwal = gpwal + gpsha(inode) *  elwal(inode)
  end do
  
  gpwal_der = 0.0_rp
  if (kfl_walld == 2 .or. kfl_walld == 3) then
    do inode = 1,pnode
      do idime = 1,ndime
        if (elwal(inode) /= 0.0_rp) then 
          gpwal_der(idime,inode) = gpsha(inode) * (elwai_coo(idime,inode) - elcod(idime,inode))/elwal(inode)
        else
          gpwal_der(idime,inode) = 0.0_rp
        endif
      enddo
    end do
  endif
    
  
  !------------------------------------------------------------------------
  !
  ! Compute coefficients for the different models
  !
  !------------------------------------------------------------------------
 
  if( TUR_SPALART_ALLMARAS ) then
     ! 
     ! Spalart-Allmaras
     !
     call tur_spaalm_der(&
                         itask,pnode,gpvis,gpden,gpsha,gpcar,elvel,elmsh,eltur,eledd,&
                         gpwal,elust,eltem,elrg2,elsqk,elgrp,elpro,elprd,gpvol,&
                         chale,gpdif,gpvel,gpgrv,gptur,gpcar_der,gpwal_der,&
                         drea_dtur,dsou_dtur,dsou_dvel,gprea_diff,gprhs_diff,gpreawal_diff,gprhswal_diff)

  else if( TUR_SST_K_OMEGA ) then
     !
     ! Menter SST k-w
     !
     call tur_sstkom_der(& 
                         pnode,gpvis,gpden,gpsha,gpcar,elvel,elmsh,eltur,eledd,&
                         gpwal,elust,eltem,elrg2,elsqk,elgrp,elpro,elprd,gpvol,&
                         chale,gpdif,gpvel,gpgrv,gptur,&
                         hleng,gpcan, elprr, pmate, &
                         drea_dtur,dsou_dtur, &
                         dsou_dvel,gpmut,dgpmut,dgpmut_dvel)
          
  endif
       
end subroutine tur_elmco2_der_all
