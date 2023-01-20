!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_assresdiff(&
                         itask,pnode,elres_diff,sens_mesh,idime_dof,ipoin_dof,eltur)
  !-----------------------------------------------------------------------
  !****f* mathru/tur_assresdiff
  ! NAME 
  !    tur_assresdiff
  ! DESCRIPTION
  !    
  ! USES
  ! USED BY
  !    ***_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime
  use def_turbul, only       :  nturb_tur
  
  implicit none
  integer(ip), intent(in)    :: itask,pnode
  real(rp),    intent(in)    :: elres_diff(pnode)
  real(rp),    intent(inout) :: sens_mesh(ndime,*)
  real(rp),    intent(in)    :: eltur(nturb_tur,pnode,*)
  integer(ip), intent(in)    :: idime_dof,ipoin_dof
  integer(ip)                :: inode

  select case(itask)
        
  case(1_ip)
     
     do inode = 1,pnode
       ! dR_t/dX Lambda_t
       sens_mesh(idime_dof,ipoin_dof) = sens_mesh(idime_dof,ipoin_dof) + elres_diff(inode)*eltur(1,inode,3)
     end do
     
  end select

end subroutine tur_assresdiff
