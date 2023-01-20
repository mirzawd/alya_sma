!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_calc_scalar_dissip_fast(&
     VECTOR_DIM,pnode,pgaus,elcon,elvel,gpsha,gpcar,&
     hleng,gptur,gpdtr,gpden,res_Chi_Yc,sgs_Chi_Yc,res_Chi_z,&
     sgs_Chi_z)

  !-----------------------------------------------------------------------
  !****f* Chemic/chm_calc_scalar_dissip
  ! NAME
  !    chm_calc_scalar_dissip_fast
  ! DESCRIPTION
  !    Compute scalar dissipation rate at Gauss points for flamelet combustion model using vectorization
  ! USES
  ! USED BY
  !    chm_post_scalar_dissipation_rate_fast
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      :  ip,rp
  use def_domain, only      :  ndime,mnode,mgaus
  use def_chemic, only      :  nclas_chm,mixedEq_eqs_chm
  use def_chemic, only      :  kfl_izmean_chm,kfl_izvar_chm,kfl_icmean_chm,kfl_icvar_chm
  use mod_chm_mixedEq, only : CHM_EQ_ZVAR
  use mod_chm_mixedEq, only : CHM_EQ_ZZ
  use mod_chm_mixedEq, only : CHM_EQ_YCVAR
  use mod_chm_mixedEq, only : CHM_EQ_YCYC

  implicit none
  integer(ip),  intent(in)  :: VECTOR_DIM
  integer(ip),  intent(in)  :: pnode,pgaus
  real(rp),     intent(in)  :: elcon(VECTOR_DIM,pnode,nclas_chm,*)
  real(rp),     intent(in)  :: elvel(VECTOR_DIM,ndime,pnode)
  real(rp),     intent(in)  :: gpsha(pnode,pgaus)
  real(rp),     intent(in)  :: gpcar(VECTOR_DIM,ndime,mnode,pgaus)
  real(rp),     intent(in)  :: hleng(VECTOR_DIM,3)
  real(rp),     intent(in)  :: gptur(VECTOR_DIM,pgaus)                  ! Turbulence viscosity at gauss points
  real(rp),     intent(in)  :: gpdtr(VECTOR_DIM,mgaus)                  ! Dt*rho at gauss points
  real(rp),     intent(in)  :: gpden(VECTOR_DIM,mgaus)                  ! Density at gauss points
  real(rp),     intent(out) :: res_Chi_Yc(VECTOR_DIM,pgaus)
  real(rp),     intent(out) :: sgs_Chi_Yc(VECTOR_DIM,pgaus)
  real(rp),     intent(out) :: res_Chi_z(VECTOR_DIM,pgaus)
  real(rp),     intent(out) :: sgs_Chi_z(VECTOR_DIM,pgaus)
  real(rp)                  :: gpcon(VECTOR_DIM,pgaus,nclas_chm,4)
  real(rp)                  :: seci4(VECTOR_DIM),delta(VECTOR_DIM),delta2(VECTOR_DIM),rdelta2(VECTOR_DIM)
  real(rp)                  :: gpgve(VECTOR_DIM,ndime,ndime,pgaus)
  real(rp)                  :: rtau_turb(VECTOR_DIM)                     ! Inverse of turbulent time scale
  real(rp)                  :: z_var(VECTOR_DIM),Yc_var(VECTOR_DIM)
  real(rp)                  :: grad_Z(VECTOR_DIM,pgaus,ndime),mod_grad_Z(VECTOR_DIM,pgaus)
  real(rp)                  :: grad_Yc(VECTOR_DIM,pgaus,ndime),mod_grad_Yc(VECTOR_DIM,pgaus)
  real(rp)                  :: diff_therm(VECTOR_DIM)        ! diffusivity assuming Le=1
  integer(ip)               :: igaus,iclas,inode,idime,jdime,ivect

  !
  ! Initialization
  !
  res_Chi_z  = 0.0_rp
  res_Chi_Yc = 0.0_rp
  sgs_Chi_z  = 0.0_rp
  sgs_Chi_Yc = 0.0_rp

  if (kfl_icmean_chm > 0 .or. kfl_izmean_chm > 0) then
     !
     ! Concentration
     !
     do iclas = 1,nclas_chm
        do igaus = 1,pgaus
           gpcon(:,igaus,iclas,1) = 0.0_rp
           do inode = 1,pnode
              gpcon(:,igaus,iclas,1) = gpcon(:,igaus,iclas,1)&
                                   + gpsha(inode,igaus) * elcon(:,inode,iclas,1)
           end do
        end do
     end do

     !
     ! Scalar dissipation rates of the reaction progress variable Yc and
     !   mixture fraction Z: X_Yc & X_Z
     !
     grad_Yc     = 0.0_rp
     mod_grad_Yc = 0.0_rp
     if (kfl_icmean_chm > 0) then
        do igaus = 1,pgaus
           do inode = 1,pnode
              do idime = 1,ndime
                 grad_Yc(:,igaus,idime) = grad_Yc(:,igaus,idime) + gpcar(:,idime,inode,igaus) * elcon(:,inode,kfl_icmean_chm,1)
              end do
           end do
        end do
     endif

     grad_Z     = 0.0_rp
     mod_grad_Z = 0.0_rp
     if (kfl_izmean_chm > 0) then
        do igaus = 1,pgaus
           do inode = 1,pnode
              do idime = 1,ndime
                 grad_Z(:,igaus,idime) = grad_Z(:,igaus,idime) + gpcar(:,idime,inode,igaus) * elcon(:,inode,kfl_izmean_chm,1)
              end do
           end do
        end do
     end if

     !
     ! Get scalar dissipation rate on each GP
     !
     do igaus = 1,pgaus

        diff_therm = gpdtr(:,igaus) / gpden(:,igaus)

        !
        ! Compute scalar dissipation rate (resolved) X_Yc = 2 * D * |grad Yc|^2
        !
        if (kfl_icmean_chm > 0) then
           do ivect=1, VECTOR_DIM
              mod_grad_Yc(ivect,igaus) = dot_product(grad_Yc(ivect,igaus,:),grad_Yc(ivect,igaus,:))
           end do
           res_Chi_Yc(:,igaus) = 2.0_rp * diff_therm / mixedEq_eqs_chm(kfl_icmean_chm) % Lewis * mod_grad_Yc(:,igaus)
        endif

        !
        ! Compute scalar dissipation rate (resolved) X_Z = 2 * D * |grad Z|^2
        !
        if (kfl_izmean_chm > 0) then
           do ivect=1, VECTOR_DIM
              mod_grad_Z(ivect,igaus) = dot_product(grad_Z(ivect,igaus,:),grad_Z(ivect,igaus,:))
           end do
           res_Chi_z(:,igaus) = 2.0_rp * diff_therm / mixedEq_eqs_chm(kfl_icmean_chm) % Lewis * mod_grad_Z(:,igaus)
        end if

        !
        ! Compute subgrid scalar dissipation rates: X_Yc^sgs & X_Z^sgs
        !
        if (kfl_icvar_chm /= 0 .or. kfl_izvar_chm /= 0 ) then

           delta(:)   = ( hleng(:,1) * hleng(:,2) * hleng(:,3) )** 0.3333333_rp

           do ivect=1, VECTOR_DIM
               delta2(ivect)  = delta(ivect)*delta(ivect)
           end do

           rdelta2(:) = 1.0_rp / delta2(:)


           !
           ! Compute strain rate square |S_ij|*|S_ij|
           !
           gpgve = 0.0_rp
           do inode = 1,pnode
              do idime = 1,ndime
                 do jdime = 1,ndime
                    gpgve(:,jdime,idime,igaus) = gpgve(:,jdime,idime,igaus) + elvel(:,idime,inode) * gpcar(:,jdime,inode,igaus)
                 end do
              end do
           end do

           seci4 = 0.0_rp
           do idime = 1,ndime                     ! |S_ij|*|S_ij|
              do jdime = 1,ndime
                 seci4(:) = seci4(:) + 0.5_rp * gpgve(:,idime,jdime,igaus) * (gpgve(:,idime,jdime,igaus)&
                     + gpgve(:,jdime,idime,igaus))
              end do
           end do

           !
           ! Compute inverse of turbulent time scale: 1 / tau = ( C_eps **2 * nu_t * |S|**2 / Delta**2 ) **(1/3)
           !
           rtau_turb(:) = ( 3.24_rp * gptur(:,igaus) * seci4(:) * rdelta2(:) ) ** (0.33333_rp)

           !
           ! Compute X_Yc (subgrid) from transport Y_c^2 equation
           !
           if (kfl_icvar_chm > 0 ) then
              if ( mixedEq_eqs_chm(kfl_icvar_chm) % kfl_eqtype == CHM_EQ_YCYC ) then
                Yc_var(:) = gpcon(:,igaus,kfl_icvar_chm,1) - gpcon(:,igaus,kfl_icmean_chm,1) * gpcon(:,igaus,kfl_icmean_chm,1)
                sgs_Chi_Yc(:,igaus) = 2.0_rp * rtau_turb(:) * Yc_var(:)

              !
              ! Compute X_Yc (subgrid) from transport Y_c variance equation
              !
              else if ( mixedEq_eqs_chm(kfl_icvar_chm) % kfl_eqtype == CHM_EQ_YCVAR ) then
                sgs_Chi_Yc(:,igaus) = 2.0_rp * rtau_turb(:) * gpcon(:,igaus,kfl_icvar_chm,1)

              end if
           end if

           !
           ! Compute X_Z (subgrid) from transport Z variance equation
           !
           if (kfl_izvar_chm > 0 ) then
              if (mixedEq_eqs_chm(kfl_izvar_chm) % kfl_eqtype == CHM_EQ_ZVAR) then
                sgs_Chi_z(:,igaus) = 2.0_rp * rtau_turb(:) * gpcon(:,igaus,kfl_izvar_chm,1)

              !
              ! Compute X_Z (subgrid) from transport Z^2 equation
              !
              elseif (mixedEq_eqs_chm(kfl_izvar_chm) % kfl_eqtype == CHM_EQ_ZZ) then
                z_var(:) = gpcon(:,igaus,kfl_izvar_chm,1) - gpcon(:,igaus,kfl_izmean_chm,1) * gpcon(:,igaus,kfl_izmean_chm,1)
                sgs_Chi_z(:,igaus) = 2.0_rp * rtau_turb(:) * z_var(:)

              endif
           endif

        end if

     end do
  endif

end subroutine chm_calc_scalar_dissip_fast
