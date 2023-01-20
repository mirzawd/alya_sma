!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine elmcar_bub(&
     pnode,pgaus,plapl,deriv,shapf_bub,deriv_bub,elcod,&
     gpsha,gpcar,gpsha_bub,gpcar_bub,ielem)
  !-----------------------------------------------------------------------
  !****f* domain/elmca2
  ! NAME
  !    elmca2
  ! DESCRIPTION
  !    This routine calculates:
  !    GPCAR: Cartesian derivatives
  !    GPVOL: Unit volume
  !    GPHES: Hessian matrix
  ! USES
  !    invmtx
  ! USED BY
  !    ***_elmope
  !    extnor
  ! SOURCE
  !-----------------------------------------------------------------------
  use def_parame, only       :  twopi
  use def_domain, only       :  ndime,mnode
  use def_elmtyp, only       :  ELCUT
  use def_kintyp, only       :  ip,rp
  use mod_cutele, only       :  elmcar_cut
  use mod_elmgeo, only       :  elmgeo_jacobian_matrix
  implicit none
  integer(ip), intent(in)    :: pnode
  integer(ip), intent(in)    :: plapl
  integer(ip), intent(in)    :: ielem
  integer(ip), intent(inout) :: pgaus
  real(rp),    intent(in)    :: deriv(ndime,pnode,*)
  real(rp),    intent(in)    :: elcod(ndime,pnode)
  real(rp),    intent(in)    :: shapf_bub(*)
  real(rp),    intent(in)    :: deriv_bub(ndime,*)
  real(rp),    intent(out)   :: gpsha(pnode,*)
  real(rp),    intent(out)   :: gpcar(ndime,mnode,*)
  real(rp),    intent(out)   :: gpsha_bub(*)
  real(rp),    intent(out)   :: gpcar_bub(ndime,*)

  integer(ip)                :: igaus,idime,jdime
  real(rp)                   :: gpdet
  real(rp)                   :: xjaci(ndime,ndime)

  do igaus = 1,pgaus
     gpsha_bub(igaus) = shapf_bub(igaus)
     call elmgeo_jacobian_matrix(&
          ndime,pnode,elcod,deriv(1,1,igaus),&
          gpdet,xjaci)
     do idime = 1,ndime
        gpcar_bub(idime,igaus) = 0.0_rp
        do jdime = 1,ndime
           gpcar_bub(idime,igaus) = gpcar_bub(idime,igaus) + xjaci(jdime,idime) * deriv_bub(jdime,igaus)
        end do
     end do
  end do

end subroutine elmcar_bub
