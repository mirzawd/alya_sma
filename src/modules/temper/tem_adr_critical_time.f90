!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_adr_critical_time(pgaus,chale,gp_velocity,gp_diffusion,gp_reaction,gp_rhs,gp_density,gp_specific_heat,gp_temperature,crit_time)
  !-----------------------------------------------------------------------
  !    This routine computes the element time step
  ! USED BY
  !    tem_elmope
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime
  use def_temper, only       :  staco_tem,source_safet_tem
  use mod_tauadr, only       :  tauadr
  implicit none
  integer(ip), intent(in)  :: pgaus
  real(rp),    intent(in)  :: chale(3)
  real(rp),    intent(in)  :: gp_velocity(ndime,pgaus)
  real(rp),    intent(in)  :: gp_diffusion(pgaus)
  real(rp),    intent(in)  :: gp_reaction(pgaus)
  real(rp),    intent(in)  :: gp_rhs(pgaus)
  real(rp),    intent(in)  :: gp_density(pgaus)
  real(rp),    intent(in)  :: gp_specific_heat(pgaus)
  real(rp),    intent(in)  :: gp_temperature(pgaus)
  real(rp),    intent(out) :: crit_time

  integer(ip)                :: idime,igaus
  real(rp)                   :: adv,dif,rea,sou,vel(pgaus)
  real(rp)                   :: source_time,rhocp
  
  adv=0.0_rp
  dif=0.0_rp
  rea=0.0_rp
  sou=0.0_rp
  vel=0.0_rp
  do igaus=1,pgaus
     do idime=1,ndime
        vel(igaus)=vel(igaus)+gp_velocity(idime,igaus)**2
     end do
     vel(igaus) = sqrt(vel(igaus))
  enddo
  do igaus=1,pgaus
     if( gp_density(igaus)*gp_specific_heat(igaus)*vel(igaus) >= adv ) then
        adv = gp_density(igaus)*gp_specific_heat(igaus)*vel(igaus)
        rhocp = gp_density(igaus)*gp_specific_heat(igaus)
     endif
     dif=max(dif, gp_diffusion(igaus))
     rea=max(rea, gp_reaction(igaus))
     sou=max(sou, abs(gp_rhs(igaus)) / (gp_density(igaus)*gp_specific_heat(igaus)*gp_temperature(igaus)) )
  end do

  call tauadr(&
       1_ip,staco_tem,adv,dif,rea,&
       chale(1),chale(2),crit_time)

  source_time = 1.e6_rp
  if (sou /= 0.0_rp) source_time = 1.0_rp / sou

  crit_time = min ( rhocp*crit_time, source_time*source_safet_tem )

end subroutine tem_adr_critical_time
