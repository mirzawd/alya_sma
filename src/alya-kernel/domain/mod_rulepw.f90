!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    mod_rulepw.f90
!> @author  houzeaux
!> @date    2020-04-04
!> @brief   Iso-parametric arrays
!> @details An n-point Gaussian quadrature rule is a quadrature rule constructed 
!>          to yield an exact result for polynomials of degree 2n-1
!>          POSGP(NDIME,NGAUS) ... Local coordinates of Gauss points
!>          WEIGP(NGAUS) ......... Weights of Gauss points
!-----------------------------------------------------------------------

module mod_rulepw

  use def_kintyp_basic,           only : ip,rp
  use def_quadrature,             only : set_quadrature

  implicit none
  private

  interface rulepw
     module procedure &
          rulepw_scalar,&
          rulepw_vector
  end interface rulepw
  
  public :: rulepw
    
contains

  subroutine rulepw_scalar(ndime,ngaus,topo,type,posgp,weigp,ierro)

    integer(ip),           intent(in)    :: ndime
    integer(ip),           intent(in)    :: ngaus
    integer(ip),           intent(in)    :: topo
    integer(ip),           intent(in)    :: type
    real(rp),              intent(out)   :: posgp(max(1_ip,ndime))
    real(rp),              intent(out)   :: weigp
    integer(ip),           intent(inout) :: ierro
    real(rp)                             :: weigp_vec(1)
    real(rp)                             :: posgp_vec(max(1_ip,ndime),1)

    call rulepw_vector(ndime,ngaus,topo,type,posgp_vec,weigp_vec,ierro)

    posgp = posgp_vec(:,1)
    weigp = weigp_vec(1)
    
  end subroutine rulepw_scalar
  
  subroutine rulepw_vector(ndime,ngaus,topo,type,posgp,weigp,ierro)

    integer(ip),           intent(in)    :: ndime
    integer(ip),           intent(in)    :: ngaus
    integer(ip),           intent(in)    :: topo
    integer(ip),           intent(in)    :: type
    real(rp),              intent(out)   :: posgp(max(1_ip,ndime),ngaus)
    real(rp),              intent(out)   :: weigp(ngaus)
    integer(ip),           intent(inout) :: ierro

    ierro = 0
    posgp = 0.0_rp
    weigp = 0.0_rp

    call set_quadrature(ndime,ngaus,topo,type,posgp,weigp,ierro)

  end subroutine rulepw_vector

end module mod_rulepw
!> @}
