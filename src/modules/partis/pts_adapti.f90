!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @file    pts_adapti.f90
!> @author  Guillaume Houzeaux
!> @date    16/11/1966   
!> @brief   Adaptive time step strategy
!> @details Compute the stretching factor 
!>          d=1.0
!>          s=1.2
!>          f(x) = x>= 1 ? (s-1.0) * tanh( (x-1.0)/(d*(s-1.0)) ) + 1 : (1.0-1.0/s) * tanh( (x-1.0)/(d*(1.0-1.0/s)) ) + 1.0
!>          set xrange[0.01:100]
!>          set log x
!>          plot f(x)
!>
!>          Normalization ....... h
!>          Stretching factor ... d,s
!>          Truncation error .... EPS_TRN
!>          Error estimate ...... EPS_ERR
!>          Error method ........ Estimate, truncation, both
!> @} 
!-----------------------------------------------------------------------
subroutine pts_adapti(&
     kfl_tstep_pts,h,x_new,x_old,u_new,u_old,u_fluid,a_new,&
     a_old,tau,dt,dto,safet_pts,alpha_str,accept_time_step,&
     kfl_modla)
 
  use def_kintyp, only : ip,rp,lg
  use def_master, only : zeror
  use def_domain, only : ndime
  use def_partis, only : beta_pts,gamma_pts
  implicit none
  integer(ip), intent(in)  :: kfl_tstep_pts
  real(rp),    intent(in)  :: h
  real(rp),    intent(in)  :: x_new(ndime)
  real(rp),    intent(in)  :: x_old(ndime)
  real(rp),    intent(in)  :: u_new(ndime)
  real(rp),    intent(in)  :: u_old(ndime)
  real(rp),    intent(in)  :: u_fluid(ndime)
  real(rp),    intent(in)  :: a_new(ndime)
  real(rp),    intent(in)  :: a_old(ndime)
  real(rp),    intent(in)  :: tau
  real(rp),    intent(in)  :: dt
  real(rp),    intent(in)  :: dto
  real(rp),    intent(in)  :: safet_pts
  real(rp),    intent(out) :: alpha_str
  logical(lg), intent(out) :: accept_time_step
  integer(ip), intent(in)  :: kfl_modla
  integer(ip)              :: idime
  real(rp)                 :: aa,dadt,aa_dadt,eps_trn,eps_err
  real(rp)                 :: dadt_idime,aa_idime,aa_dadt_idime
  real(rp)                 :: dt_err,s,d,alpha,uu,uu_idime
!  real(rp)                 :: dt_trn
  real(rp)                 :: uu_dadt
  real(rp)                 :: diffu_idime,diffu,diffa_idime,diffa

  accept_time_step = .true.
  s       = 1.2_rp
  d       = 2.5_rp
  eps_err = 1.0e-3_rp
  eps_trn = 1.0_rp

  if( kfl_tstep_pts == 1 ) then

     !-------------------------------------------------------------------
     !
     ! Error estimate based on position
     !
     !-------------------------------------------------------------------
     !
     ! alpha = dt_err / dt
     ! if alpha > 1, means that the time step can be greater
     !
     dadt    = 0.0_rp
     aa      = 0.0_rp
     aa_dadt = 0.0_rp
     uu      = 0.0_rp
     diffa   = 0.0_rp
     diffu   = 0.0_rp

     do idime = 1,ndime
        diffa_idime   = abs(a_new(idime)-a_old(idime))
        dadt_idime    = diffa_idime/dt + zeror
        aa_idime      = abs(a_new(idime)) + zeror
        aa_dadt_idime = aa_idime / dadt_idime
        diffu_idime   = abs(u_new(idime)-u_old(idime))
        aa            = max(aa     ,aa_idime)
        dadt          = max(dadt   ,dadt_idime)
        aa_dadt       = max(aa_dadt,aa_dadt_idime)
        diffu         = max(diffu,diffu_idime)
        diffa         = max(diffa,diffa_idime)
     end do 
     if( kfl_modla == 2 ) then
        !
        ! NOSOTROS
        !
        !dt_err = ( h * eps_err / ( dadt * max(0.1_rp,abs(1.0_rp/6.0_rp - beta_pts)) ) ) ** (1.0_rp/3.0_rp)
        !dt_err = ( 6.0_rp * h * eps_err / dadt ) ** (1.0_rp/3.0_rp)
        dt_err = sqrt( eps_err * h / abs(1.0_rp/6.0_rp - beta_pts) / ( diffa + zeror ) )
        !dt_trn = eps_trn * 2.0_rp * (diffu+zeror)/(diffa+zeror)
        !dt_err = min(dt_err,dt_trn)
        !
        ! RICKELT
        !
        !dt_err = ( eps_err * gamma_pts * h / beta_pts / (uu+zeror) )
        !
        ! ZIENKIEWICZ & XIE 
        !
        !dt_err = ( eps_err * h / (uu+zeror) * 2.0_rp * gamma_pts / abs( 2.0_rp * beta_pts - gamma_pts ) )
     else
        dt_err = max( eps_err * h * 4.0_rp / ( dto * ( diffa + zeror ) ) , sqrt( eps_err * h * 6.0_rp / ( diffa + zeror )))
     end if

     alpha  = dt_err/ dt

  else if( kfl_tstep_pts == 3 ) then

     !-------------------------------------------------------------------
     !
     ! Safety factor
     !
     !-------------------------------------------------------------------

     alpha = safet_pts * tau / dt

  else if( kfl_tstep_pts == 4 ) then

     !-------------------------------------------------------------------
     !
     ! Error estimate based on velocity
     !
     !-------------------------------------------------------------------
     !
     ! alpha = dt_err / dt
     ! if alpha > 1, means that the time step can be greater
     !
     uu_dadt = 0.0_rp
     uu = 0.0_rp
     do idime = 1,ndime
        uu = uu + u_fluid(idime)**2
     end do
     uu = sqrt(uu)

     do idime = 1,ndime
        dadt_idime = abs(a_new(idime)-a_old(idime))/dt + zeror
        uu_idime   = abs(uu/dadt_idime)
        uu_dadt    = max(uu_dadt,uu_idime)
     end do 

     dt_err  = sqrt( uu_dadt * eps_err / abs(0.5_rp-gamma_pts) ) 
     alpha   = dt_err / dt

  end if
  !
  ! Accept or not?
  ! In principle, if alpha is near one, we could accept the solution and be less restrictive
  !
  ! alpha > 1 ... dt is incerased
  ! alpha < 1 ... dt is decreased
  !
  if( alpha >= 1.0_rp ) then
     alpha_str = (s-1.0_rp) * tanh( (alpha-1.0_rp)/(d*(s-1.0_rp)) ) + 1.0_rp
  else
     alpha_str = (1.0_rp-1.0_rp/s) * tanh( (alpha-1.0_rp)/(d*(1.0_rp-1.0_rp/s)) ) + 1.0_rp
  end if
  if( alpha >= 0.9_rp ) then
     accept_time_step = .true.
  else
     accept_time_step = .false.
  end if

end subroutine pts_adapti
