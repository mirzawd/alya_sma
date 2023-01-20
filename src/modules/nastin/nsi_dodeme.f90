!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_dodeme(itask,Auu,Aup,Apu,App,Q,bu,bp)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_dodeme
  ! NAME 
  !    nsi_dodeme
  ! DESCRIPTION
  !    Interfaces Nastin with Dodeme
  !    ITASK = 1 ... Change pressure boundary conditions
  !          = 2 ... Impose Dirichlet b.c.
  !                  Impose hole nodes velocity and pressure to zero
  ! USES
  ! USED BY
  !    nsi_dodeme
  !***
  !----------------------------------------------------------------------- 
  use def_master
  use def_elmtyp
  use def_domain
  use def_solver
  use def_nastin
  implicit none
  integer(ip), intent(in)    :: itask
  real(rp),    intent(inout) :: Auu(ndime,ndime,nzdom) 
  real(rp),    intent(inout) :: Aup(ndime,nzdom)
  real(rp),    intent(inout) :: Apu(ndime,nzdom)
  real(rp),    intent(inout) :: App(*) 
  real(rp),    intent(inout) :: Q(*)
  real(rp),    intent(inout) :: bu(ndime,*),bp(*) 
!  real(rp)                   :: qdiag,appdiag,auudiag(3),time1,p,u
!  integer(ip)                :: ipoin,ii,idime,ibopo,intij,iboun,jzdom,jpoin,ppoin
!  integer(ip)                :: izdom,kpoin,jdime,ielem,inode,imeth,inodb
!  integer(ip)                :: isubd,jsubd,itype,pnode,pnodb,pblty


end subroutine nsi_dodeme
