!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_elmmat_der_all(                                &
	      pnode,pgaus,pevat,gpden,gpvis,gppor,gpgvi,gpsp1, &
	      gptt1,gpsp2,gptt2,gpvol,gpsha,gpcar,gplap,gphes, &
	      gpadv,gpvep,gpprp,gpgrp,elauu,elaup,elapp,elapu,gpst1,&
	      gpstrm,gpstrc,elrbu,elrbp,elvel,ielem,gptem,gpgde,&
	      gpdvi,elaut,elapt,elpre,gpgve,p1vec,gpgdv,gpdde,gpgdd,chale,&
	      dgpmut_dvel,dgpmut_dtur)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmmat_der_all
  ! NAME 
  !    nsi_elmmat_der_all
  ! DESCRIPTION
  !    Compute element matrix and RHS related to the exact linealization
  ! USES
  ! USED BY
  !    nsi_elmope_omp
  !***
  !----------------------------------------------------------------------- 
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,mnode,ntens
  use def_nastin, only       :  nbdfp_nsi
  use def_nastin, only       :  ncomp_nsi
  use def_master, only       :  ID_TEMPER,ID_CHEMIC,ID_NASTIN,ID_TURBUL,kfl_coupl
  use mod_ker_proper
  
  implicit none
  integer(ip), intent(in)    :: pnode,pgaus,pevat,ielem
  real(rp),    intent(in)    :: gpden(pgaus),gpvis(pgaus)
  real(rp),    intent(in)    :: gppor(pgaus),gpgvi(ndime,pgaus),gpdvi(pnode,pgaus),gpgdv(ndime,pnode,pgaus)
  real(rp),    intent(inout) :: gpsp1(pgaus),gptt1(pgaus)
  real(rp),    intent(inout) :: gpsp2(pgaus),gptt2(pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus),gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gpadv(ndime,pgaus)
  real(rp),    intent(inout) :: gpvep(ndime,pgaus)
  real(rp),    intent(in)    :: gpprp(pgaus)      
  real(rp),    intent(in)    :: gpgrp(ndime,pgaus)
  real(rp),    intent(in)    :: gplap(pnode,pgaus)
  real(rp),    intent(in)    :: gphes(ntens,mnode,pgaus)
  real(rp),    intent(inout) :: elaut(pnode*ndime,pnode)
  real(rp),    intent(inout) :: elauu(pnode*ndime,pnode*ndime)
  real(rp),    intent(inout) :: elaup(pnode*ndime,pnode)
  real(rp),    intent(inout) :: elapp(pnode,pnode)
  real(rp),    intent(inout) :: elapu(pnode,pnode*ndime)
  real(rp),    intent(inout) :: elrbu(ndime,pnode)
  real(rp),    intent(inout) :: elrbp(pnode)
  real(rp),    intent(in)    :: gptem(pgaus, nbdfp_nsi)

  real(rp),    intent(inout) :: elapt(pnode,pnode)
  real(rp),    intent(in)    :: gpst1(pgaus)
  real(rp),    intent(in)    :: gpstrm(ndime,pgaus)
  real(rp),    intent(in)    :: gpstrc(pgaus)
  real(rp),    intent(in)    :: chale(2)
  real(rp),    intent(in)    :: elvel(ndime,pnode,ncomp_nsi),elpre(pnode)
  real(rp),    intent(in)    :: gpgde(ndime,pgaus)                    ! grad(den)
  real(rp),    intent(in)    :: dgpmut_dvel(ndime,pnode,pgaus)        ! Turbulence viscosity derivatives w.r.t velocity
  real(rp),    intent(in)    :: dgpmut_dtur(1,pnode,pgaus)            ! Turbulence viscosity derivatives w.r.t turbulence unk
  real(rp),    intent(in)    :: gpdde(pnode,pgaus)                    ! Density derivatives w.r.t nodal temperature
  real(rp),    intent(in)    :: gpgdd(ndime,pnode,pgaus)              ! Density derivatives w.r.t nodal temperature and coordinates
  real(rp),    intent(in)    :: p1vec(pnode,pgaus)
  real(rp),    intent(in)    :: gpgve(ndime,ndime,pgaus)

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!        ELAUU   ---->  ELAUU  +  dR_m/du               !!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!        ELAUP   ---->  ELAUP  +  dR_m/dp               !!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!        ELAPU   ---->  ELAPU  +  dR_c/du               !!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!        ELAPP   ---->  ELAPP  +  dR_c/dp               !!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  call nsi_elmmat_st(                                    &
		    pnode,pgaus,pevat,gpden,gpvis,gppor,gpgvi,gpsp1, &
		    gptt1,gpsp2,gptt2,gpvol,gpsha,gpcar,gplap,gphes, &
		    gpadv,gpvep,gpprp,gpgrp,elauu,elaup,elapp,elapu,gpst1,&
		    gpstrm,gpstrc,chale,elrbu,elrbp,elvel,ielem,gptem,gpgde,dgpmut_dvel,gpgve)
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!  RhsadjTem_nsi (to be sent to adjoint temper)   <------  Trans[dRm/dT] [Lambda_u]  +   Trans[dRc/dT] [Lambda_p]   !!!!!!!
  !!!!!!!!!!!!!!!              elrbu   <-------   elrbu   -    Trans[dRt/du] [Lambda_t]       (sent by temper)               !!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Coupling with TEMPER
  !    
  if (kfl_coupl(ID_NASTIN,ID_TEMPER) >= 1 ) then

    call nsi_elmmat_tem(                                &
	pnode,pgaus,pevat,gpden,gpvis,gppor,gpgvi,gpdvi,gpsp1, &
	gptt1,gpsp2,gptt2,gpvol,gpsha,gpcar,gplap,gphes, &
	gpadv,gpvep,gpprp,gpgrp,elaut,elapt,gpst1,&
	gpstrm,gpstrc,chale,elvel,elpre,gpgve,ielem,p1vec,gpgdv,&
	gpgde,gpdde,gpgdd,elrbu)
	
  endif
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!      elrbu   <-------   elrbu   -   Trans[dRs/du] [Lambda_s]  (sent by chemic)         !!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Coupling with CHEMIC
  ! 
  if (kfl_coupl(ID_CHEMIC,ID_NASTIN) >= 1 ) then
  
    call nsi_elmmat_chm(pnode,ielem,elrbu)
    
  endif
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!  RhsadjTur_nsi(1) (to be sent to adjoint turbul)   <------  Trans[dRm/dk] [Lambda_u]  +   Trans[dRc/dk] [Lambda_p] !!!!!!!
  !!!!!!  RhsadjTur_nsi(2) (to be sent to adjoint turbul)   <------  Trans[dRm/dw] [Lambda_u]  +   Trans[dRc/dw] [Lambda_p]  !!!!!!
  !!!!!!!!!!!!!!!        elrbu   <-------   elrbu   -    Trans[dRk/du] [Lambda_k] - Trans[dRw/du] [Lambda_w] (sent by turbul)!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Coupling with turbul
  !    
  if (kfl_coupl(ID_NASTIN,ID_TURBUL) >= 1 ) then

    call nsi_elmmat_tur(                                &
	pnode,pgaus,pevat,gpden,gpvis,gppor,gpgvi,gpdvi,gpsp1, &
	gptt1,gpsp2,gptt2,gpvol,gpsha,gpcar,gplap,gphes, &
	gpadv,gpvep,gpprp,gpgrp,gpst1,&
	gpstrm,gpstrc,chale,elvel,elpre,gpgve,ielem,p1vec,gpgdv,&
	gpgde,gpdde,gpgdd,elrbu,dgpmut_dtur)
	
  endif
  
end subroutine nsi_elmmat_der_all
