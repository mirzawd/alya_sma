!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine pts_rungk4(  &
     ndime,kfl_drafo,kfl_extfo,grafo,buofo,gravi,veloc_fluid,  &
     visco_fluid,densi_fluid,v,x,densi_parti,                  &
     diame_parti,spher_parti,dista,t,dt,vf,xf)
  !
  ! Returns final (position, velocity) tuple after time dt has passed
  !
  use def_kintyp, only : ip,rp
  implicit none
  integer(ip), intent(in)  :: ndime
  integer(ip), intent(in)  :: kfl_drafo
  integer(ip), intent(in)  :: kfl_extfo
  real(rp),    intent(in)  :: grafo
  real(rp),    intent(in)  :: buofo
  real(rp),    intent(in)  :: gravi(*)
  real(rp),    intent(in)  :: veloc_fluid(*)
  real(rp),    intent(in)  :: visco_fluid
  real(rp),    intent(in)  :: densi_fluid
  real(rp),    intent(in)  :: v(*)             !< Initial velocity
  real(rp),    intent(in)  :: x(*)             !< Initial position
  real(rp),    intent(in)  :: densi_parti
  real(rp),    intent(in)  :: diame_parti
  real(rp),    intent(in)  :: spher_parti
  real(rp),    intent(in)  :: dista            !< Time 
  real(rp),    intent(in)  :: t                !< Time 
  real(rp),    intent(in)  :: dt               !< Time step
  real(rp),    intent(out) :: vf(*)            !< Final velocity
  real(rp),    intent(out) :: xf(*)            !< Final position

  integer(ip)              :: idime
  real(rp)                 :: x1(3),v1(3),a1(3)
  real(rp)                 :: x2(3),v2(3),a2(3)
  real(rp)                 :: x3(3),v3(3),a3(3)
  real(rp)                 :: x4(3),v4(3),a4(3)

  do idime = 1,ndime
     x1(idime) = x(idime)
     v1(idime) = v(idime)
  end do

  call pts_accele(                                               &
       ndime,kfl_drafo,kfl_extfo,grafo,buofo,gravi,veloc_fluid,  &
       visco_fluid,densi_fluid,v1,x1,densi_parti,                &
       diame_parti,spher_parti,dista,t,a1)

  do idime = 1,ndime
     x2(idime) = x(idime) + 0.5_rp * v1(idime) * dt
     v2(idime) = v(idime) + 0.5_rp * a1(idime) * dt
  end do

  call pts_accele(                                               &
       ndime,kfl_drafo,kfl_extfo,grafo,buofo,gravi,veloc_fluid,  &
       visco_fluid,densi_fluid,v2,x2,densi_parti,                &
       diame_parti,spher_parti,dista,t,a2)

  do idime = 1,ndime
     x3(idime) = x(idime) + 0.5_rp * v2(idime) * dt
     v3(idime) = v(idime) + 0.5_rp * a2(idime) * dt
   end do

  call pts_accele(                                               &
       ndime,kfl_drafo,kfl_extfo,grafo,buofo,gravi,veloc_fluid,  &
       visco_fluid,densi_fluid,v3,x3,densi_parti,                &
       diame_parti,spher_parti,dista,t,a3)
   
  do idime = 1,ndime
     x4(idime) = x(idime) + v3(idime) * dt
     v4(idime) = v(idime) + a3(idime) * dt
  end do
 
  call pts_accele(                                               &
       ndime,kfl_drafo,kfl_extfo,grafo,buofo,gravi,veloc_fluid,  &
       visco_fluid,densi_fluid,v4,x4,densi_parti,                &
       diame_parti,spher_parti,dista,t,a4)

  do idime = 1,ndime
     xf(idime) = x(idime) + (dt/6.0_rp) * (v1(idime) + 2.0_rp * v2(idime) + 2.0_rp * v3(idime) + v4(idime) )
     vf(idime) = v(idime) + (dt/6.0_rp) * (a1(idime) + 2.0_rp * a2(idime) + 2.0_rp * a3(idime) + a4(idime) )
  end do
  
end subroutine pts_rungk4
