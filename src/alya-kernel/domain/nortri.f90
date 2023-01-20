!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nortri(p1,p2,p3,coord,vec,ndime)
!-----------------------------------------------------------------------
!****f* Domain/nortri
! NAME
!    nortri
! DESCRIPTION
!    This routine computes the boundary normals
! USES
! USED BY
!    bounor
!***
!-----------------------------------------------------------------------
  use      def_master
  use      def_domain, only : npoin
  use      def_kintyp, only : ip,rp
  implicit none
  integer(ip), intent(in)  :: p1,p2,p3,ndime
  real(rp),    intent(in)  :: coord(ndime,*)
  real(rp),    intent(out) :: vec(3,3)
  !if (p1 < 0 .or. p2 < 0 .or. p3 < 0) print *,'nortri', npoin, p1, p2, p3

  vec(1,1) = coord(1,p2) - coord(1,p1)
  vec(2,1) = coord(2,p2) - coord(2,p1)
  vec(3,1) = coord(3,p2) - coord(3,p1)
  vec(1,2) = coord(1,p3) - coord(1,p1)
  vec(2,2) = coord(2,p3) - coord(2,p1)
  vec(3,2) = coord(3,p3) - coord(3,p1)
  call vecpro(vec(1,1),vec(1,2),vec(1,3),ndime)

end subroutine nortri
