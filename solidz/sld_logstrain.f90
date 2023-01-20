!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_logstrain(gpgdi,gpgre,logstrain)

  use def_kintyp, only      : ip,rp,lg
  use def_domain, only      : ndime
  use mod_sld_stress_model_comput, only : SM_polar_decomposition

  implicit none
  real(rp), intent(in)       :: gpgdi(ndime,ndime)
  real(rp), intent(in)       :: gpgre(ndime,ndime)
  real(rp), intent(out)      :: logstrain(ndime,ndime)

  real(rp)                   :: rstr1(ndime,ndime), udefo(ndime,ndime)  ! F = R*U
  real(rp)                   :: gpfin(ndime,ndime)   !b = F F^t
  real(rp)                   :: gpeva(ndime), gpeve(ndime,ndime)
  real(rp)                   :: bidon
  integer(ip)                :: idime,jdime,kdime

   !Get R and U (F=R*U : gpgdi=rstr*udefo)
   call SM_polar_decomposition(1_ip,ndime,gpgdi(1,1),rstr1(1,1),udefo(1,1))

   !Calculate b=F F^T
   do kdime=1,ndime
     do jdime=1,ndime
       gpfin(jdime,kdime) = 0.0_rp
       do idime=1,ndime
         gpfin(jdime,kdime) = gpfin(jdime,kdime) + gpgdi(jdime,idime)*gpgdi(kdime,idime)
       end do
     end do
   end do

   !Find eigenvalues and eigenvectors of b
   call spcdec(gpfin(1,1),gpeva(1),gpeve(1,1),bidon,1_ip,'SLD_LOGSTRAIN')

   !Calculate epsilon_log
   logstrain = 0.0_rp
   do kdime=1,ndime
     do idime=1,ndime
       do jdime=1,ndime
         logstrain(idime,jdime) = logstrain(idime,jdime) + &
              log(sqrt(gpeva(kdime)))*gpeve(idime,kdime)*gpeve(jdime,kdime)
       end do
     end do
   end do

end subroutine
