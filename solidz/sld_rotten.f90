!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_rotten(irotype,tenold,baslo,tennew,ndime)
  !-----------------------------------------------------------------------
  !****f* Solidz/sld_stress_model_151
  ! NAME
  !    sld_stress_model_151
  ! DESCRIPTION
  !This subroutine perdoms rotations and inverse rotations on tenold
  !
  !tennew(i,j)= baslo(i,p)^(irotype) * baslo(j,q)^(irotype) * tenold(p,q)
  !
  ! IMPLEMENTED BY
  !    Adrià Quintanas (u1902492@campus.udg.edu)
  ! VERSION
  !    October 2014 
  !***
  !-----------------------------------------------------------------------

  use      def_kintyp  

  implicit none

  integer(ip),    intent(in)     ::    irotype,ndime
  real(rp),       intent(in)     ::    baslo(ndime,ndime), tenold(ndime,ndime)
  real(rp),       intent(out)    ::    tennew(ndime,ndime)
  integer(ip)                    ::    i,j,p,q

  ! -------------------------------
  ! INITIALIZE VARIABLES
  !
  tennew = 0.0_rp


  ! -------------------------------
  ! ROTATING THE OLD TENSOR IN FUNCTION OF THE BASIS
  !

  !
  ! From LCSYS to GCSYS
  ! T_{ij} = Q_{ip}^(-1) * Q_{jq}^(-1) * T_{pq}
  !
  if (irotype.lt.0) then

     do i = 1, ndime
        do j = 1, ndime
           do p = 1, ndime
              do q = 1, ndime
                 tennew(i, j) = tennew(i, j) + baslo(p, i)*baslo(q, j)*tenold(p, q)
              end do
           end do
        end do
     end do

     !
     ! From GCSYS to LGCSYS
     !	 T_{ij} = Q_{ip} * Q_{jq} * T_{pq}
     !
  else if (irotype.gt.0) then

     do i = 1, ndime
        do j = 1, ndime
           do p = 1, ndime
              do q = 1, ndime
                 tennew(i,j) = tennew(i, j) + baslo(i, p)*baslo(j, q)*tenold(p, q)
              end do
           end do
        end do
     end do

  end if

end subroutine sld_rotten
