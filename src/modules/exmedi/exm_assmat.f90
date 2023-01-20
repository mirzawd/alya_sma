!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!subroutine exm_assmat(&
!     elmat,amatr,pnode,ndime,ndofn,mevat,nequa,lnods,lpexm,algso) 
!-----------------------------------------------------------------------
!
! This routine performs the assembly of the elementary matrices 
! according to the selected storage method
!
!-----------------------------------------------------------------------
!  use      def_kintyp
!  use      def_solver
!
!  use      Mod_skyase
!
!  implicit none
!  integer(ip), intent(in)    :: pnode,ndime,ndofn,nequa,mevat
!  integer(ip), intent(in)    :: algso
!  integer(ip), intent(in)    :: lpexm(nequa)
!  integer(ip), intent(in)    :: lnods(pnode)
!  real(rp),    intent(in)    :: elmat(mevat,mevat)
!  real(rp),    intent(inout) :: amatr(*)
!
!  if(algso==0) then
!     call skyase(elmat,amatr,pnode,ndime,ndofn,mevat,nequa,lnods,lpexm,2_ip)
!  else                                         
!     !call csrase(iposi,jposi,r_sol,c_sol,amatr,acont,1)
!  end if
!  
!end subroutine exm_assmat
