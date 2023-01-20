!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_assresdiff(&
                         itask,pnode,pevat,lnods,elrbu,elrbp,sens_mesh,idime_dof,ipoin_dof,elvel,elpre)
  !-----------------------------------------------------------------------
  !****f* mathru/nsi_assresdiff
  ! NAME 
  !    nsi_assresdiff
  ! DESCRIPTION
  !    Assembly an elemental matrix ELMAT in global matrix AMATR
  ! USES
  ! USED BY
  !    ***_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime
!  use def_domain, only       :  npoin
!  use def_kermod, only       :  kfl_ndvars_opt
  implicit none
  integer(ip), intent(in)    :: itask,pnode,pevat
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(in)    :: elrbu(ndime,pnode)
  real(rp),    intent(in)    :: elrbp(pnode)
!   real(rp),    intent(inout) :: resdiff(kfl_ndvars_opt,*)
  real(rp),    intent(inout) :: sens_mesh(ndime,*)
  real(rp),    intent(in)    :: elvel(ndime, pnode,*)
  real(rp),    intent(in)    :: elpre(pnode)
!   integer(ip), intent(in)    :: idesvar
  integer(ip), intent(in)    :: idime_dof,ipoin_dof
  integer(ip)                :: inode
  integer(ip)                :: idime
!  integer(ip)                :: ipoin,ind

  select case(itask)
  
!   case(1_ip)
!      !
!      ! Momentum and continuity RHS only
!      !
!      do inode = 1,pnode
!         ipoin = lnods(inode)
!         do idime = 1,ndime
!            ind = (ipoin-1)*ndime + idime
!            !$OMP ATOMIC
!            resdiff(idesvar,ind) = resdiff(idesvar,ind) + elrbu(idime,inode)
!         end do
!         !$OMP ATOMIC
!         resdiff(idesvar,npoin*ndime + ipoin) = resdiff(idesvar,npoin*ndime + ipoin) + elrbp(inode)
!      end do
     
  case(1_ip)
     !
     ! Momentum and continuity RHS only
     !
     do inode = 1,pnode
       do idime = 1,ndime
         ! dR_m/dX Lambda_m
         !$OMP ATOMIC
         sens_mesh(idime_dof,ipoin_dof) = sens_mesh(idime_dof,ipoin_dof) + elrbu(idime,inode)*elvel(idime,inode,2)
       enddo
       ! dR_c/dX Lambda_c
       !$OMP ATOMIC
       sens_mesh(idime_dof,ipoin_dof) = sens_mesh(idime_dof,ipoin_dof) + elrbp(inode)*elpre(inode)
     end do

  end select

end subroutine nsi_assresdiff
