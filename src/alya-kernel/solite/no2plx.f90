!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine no2plx(nbnodes,nbvar,xx,sumxx)

  !---------------------------------------------------------------------------------
  ! Sources/kernel/solite/no2plx.f90
  ! NAME 
  !    no2plx
  ! DESCRIPTION
  !    This routine computes the 2-norm (Eucledian norm) of a complex vector XX:
  !    SUMXX = ||XX||_2 = ( sum_i |XX_i|^2 )^{1/2}
  ! INPUT ARGUMENTS
  !    NBNODES .................. Number of nodes
  !    NBVAR .................... Number of variables per node
  !    XX(NBVAR*NBNODES) ........ Vector
  ! OUTPUT ARGUMENTS
  !    SUMXX .................... Eucledian norm
  ! USES
  ! USED BY 
  !***
  !---------------------------------------------------------------------------------


  use def_kintyp, only : ip,rp
  use def_master, only : kfl_paral,npoi1,npoi2,npoi3
  use mod_communications, only : PAR_SUM
  !Declaration statements
  implicit none

  !Dummy arguments
  integer(ip), intent(in)  :: nbnodes,nbvar
  complex(rp), intent(in)  :: xx(*)
  real(rp),    intent(out) :: sumxx

  !Local variables
  integer(ip)              :: ii,temp

  !--------------------------------------------------------------------------------------------------------
  ! Euclidean norm: sumxx = ||x||
  !--------------------------------------------------------------------------------------------------------
  sumxx = 0.0_rp

  if (kfl_paral == -1_ip) then
    !Sequential
! !$omp parallel do                 &
! !$omp            default(shared)  &
! !$omp            private(ii)      &
! !$omp            schedule(static) &
! !$omp            reduction(+:sumxx)
    do ii = 1,nbvar*nbnodes
      sumxx = sumxx + real(conjg(xx(ii)) * xx(ii),rp)       !Inner product: sumxx = <xx, xx> = ||x||^2
    enddo
! !$omp end parallel do
  else if (kfl_paral >= 1_ip) then
    !Parallel: slaves
! !$omp parallel do                 &
! !$omp            default(shared)  &
! !$omp            private(ii)      &
! !$omp            schedule(static) &
! !$omp            reduction(+:sumxx)
    do ii = 1,nbvar*npoi1
      sumxx = sumxx + real(conjg(xx(ii)) * xx(ii),rp)
    enddo
! !$omp end parallel do
    temp = nbvar*(npoi2 - 1_ip) + 1_ip
! !$omp parallel do                 &
! !$omp            default(shared)  &
! !$omp            private(ii)      &
! !$omp            schedule(static) &
! !$omp            reduction(+:sumxx)
    do ii = temp,nbvar*npoi3
      sumxx = sumxx + real(conjg(xx(ii)) * xx(ii),rp)
    enddo
! !$omp end parallel do
  endif
  
  !
  ! Parallel: reduce sum
  !
  if( kfl_paral >= 0 ) then
     call PAR_SUM(sumxx)
  end if

  sumxx = sqrt(sumxx)

end subroutine no2plx

