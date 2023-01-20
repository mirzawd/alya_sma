!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_turbul(&
     itask,jtask,pnode,pgaus,igaui,igauf,&
     gpsha,gpcar,elvel,gpden,gpvis,gpmut, &
     gpgvi,grvis,gpgve,ielem,kfl_kxmod)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_turbul
  ! NAME 
  !    nsi_turbul
  ! DESCRIPTION
  !    Compute viscosity and its gradient due to turbulence
  !    GBMUT ....... mut
  !    JTASK = 1 ... Compute gradient GRVIS = grad(mut)
  !            0 ... Do not compute gradient GRVIS = grad(mut)
  !    ITASK = 1 ... GPVIS <= GPVIS + GPMUT
  !            0 ... Do not change GPVIS
  !    GPGVI = Viscosity gradient grad(mu+mut)
  ! USES
  !
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp,nvart  
  use def_domain, only       :  mnode,ndime
  use def_master, only       :  zeror
  use def_nastin, only       :  zensi,turmu_nsi
  use def_kermod, only       :  turmu_ker
  use mod_ker_regularization, only : regul_k, regul_e
  use mod_output_postprocess, only : output_postprocess_check_variable_postprocess
  use mod_arrays,             only : arrays_number
  implicit none 
  integer(ip),  intent(in)   :: itask,jtask,pnode,pgaus,igaui,igauf
  integer(ip),  intent(in)   :: kfl_kxmod, ielem
  real(rp),     intent(in)   :: gpsha(pnode,pgaus)
  real(rp),     intent(in)   :: gpcar(ndime,mnode,pgaus)
  real(rp),     intent(in)   :: elvel(ndime,pnode)
  real(rp),     intent(in)   :: gpden(pgaus)
  real(rp),     intent(out)  :: gpvis(pgaus)
  real(rp),     intent(inout):: gpmut(pgaus)
  real(rp),     intent(inout):: gpgvi(ndime,pgaus)
  real(rp),     intent(out)  :: grvis(ndime,pgaus)
  real(rp),     intent(out)  :: gpgve(ndime,ndime,pgaus)
  integer(ip)                :: igaus,idime


  if(turmu_ker % kfl_exist /= 0_ip) then ! matias suggested I use this - seems good idea 
     
     !
     ! Compute mu_t = rho * nu_t for LES
     !
     gpmut(igaui:igauf) = gpden(igaui:igauf) * gpmut(igaui:igauf)    ! DEJARLO CONSISTENTE   en element operations esta fuera!!!!!!
     !
     ! Compute effective viscosity & gradient of viscosity: 
     !    mu_eff        = mu + mu_t 
     !    grad (mu_eff) = grad(mu) + grad(mu_t)
     !
     if (itask == 1) then
        do igaus = igaui,igauf
           gpvis(igaus) = gpvis(igaus) + gpmut(igaus)
           do idime = 1,ndime
              gpgvi(idime,igaus) = gpgvi(idime,igaus) + grvis(idime,igaus)
           end do
        end do
     end if

     !----------------------------------------------------------------------
     !
     ! Needed for postprocessing turbulent viscosity
     !
     !----------------------------------------------------------------------
     if ( jtask == 1 ) then
        if(         output_postprocess_check_variable_postprocess(arrays_number('TURMU'))  &
             & .or. output_postprocess_check_variable_postprocess(arrays_number('AVMUT'))  &
             & .or. output_postprocess_check_variable_postprocess(arrays_number('AVSTX'))  &
             & .or. output_postprocess_check_variable_postprocess(arrays_number('AVSTY'))  &
             & .or. output_postprocess_check_variable_postprocess(arrays_number('AVSTZ'))  &
             & .or. output_postprocess_check_variable_postprocess(arrays_number('ENMUT'))  &
             & .or. output_postprocess_check_variable_postprocess(arrays_number('ENSTX'))  &
             & .or. output_postprocess_check_variable_postprocess(arrays_number('ENSTY'))  &
             & .or. output_postprocess_check_variable_postprocess(arrays_number('ENSTZ'))  ) then
           do igaus = igaui,igauf
              turmu_nsi(ielem)%a(1,igaus,1) = gpmut(igaus)
           end do
        end if
     end if

  end if

end subroutine nsi_turbul
