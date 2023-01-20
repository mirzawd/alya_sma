!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @{
!> @addtogroup Mathru
!> @file    mod_shape_descriptors.f90
!> @author  margarida moragues
!> @date    2020-02-06
!> @brief   Shape descriptors
!> @details Compute shape descriptors
!> @}
!-----------------------------------------------------------------------

module mod_shape_descriptors

  use def_kintyp,      only: ip, rp
  use def_parame
  
  use def_domain
  use def_master
  use mod_func
  use mod_integrals
  
  implicit none

  public :: compute_compactness_2
  public :: compute_moments  

contains

  !-----------------------------------------------------------------------
  !> @{ 
  !> @author  margarida moragues
  !> @date    2020-02-06
  !> @brief   Compactness_2
  !> @details Compactness_2:
  !> This program computes the Compactness_2 for a group clusters (or subdomains).
  !> Note that for any element, ielem, legro(ielem) is the cluster (or subdomain) 
  !> number to which ielem belongs. 
  !> The Compactness_2 or Circularity is a measure of how much a shape resembles 
  !> a sphere in 3D or a circle in 2D. It is called Compactness_2 because is a
  !> modification of the standard compactness.
  !>
  !> The standard compactness of a given shape S is defined as:
  !>
  !>    Compactness_Std = 4 * Pi * Area(S) / Perimeter(S).
  !>
  !> The Compactness_2 in 3D is defined as [Martínez-Ortiz 2010]:
  !>
  !>    Compactness_2 = (3 * Volume(S))^(5/3) / 
  !>                       [5 * (4*Pi)^(2/3) * (Mu_2_0_0(S) + Mu_0_2_0(S) + Mu_0_0_2(S))]
  !>
  !> and in 2D is defined as [Zunic et al. 2010]:
  !>
  !>    Compactness_2 (or Circularity) = Volume(S)^2 / [(2*Pi) * (Mu_2_0(S) + Mu_0_2(S))],
  !>
  !> where Mu are the central moments of inertia (described below).
  !>
  !> References:
  !> C. Martínez-Ortiz (2010) 2D and 3D shape descriptors. PhD Thesis. 
  !> J. Zunic et al. (2010) A hu moment invariant as a shape circularity measure. 
  !>        Pattern Recognition, 43(1):47–57. 21, 25
  !> @}
  !-----------------------------------------------------------------------

  subroutine compute_compactness_2(legro,nclus,concentration_measure,volume,centroid,compactness_2)

    use def_domain
    use def_master
    use mod_func
    use mod_integrals

    implicit none
    integer(ip),          pointer,           intent(in)      :: legro(:)
    integer(ip),                             intent(in)      :: nclus
    real(rp),             pointer,           intent(in)      :: concentration_measure(:) 
    real(rp),             pointer,           intent(inout)   :: volume(:)
    real(rp),             pointer,           intent(inout)   :: centroid(:,:)
    real(rp),             pointer,           intent(inout)   :: compactness_2(:)
    integer(ip)                                              :: idime, iclus
    real(rp),             allocatable                        :: moments_1(:,:), moments_2(:,:), central_moments_2(:,:)
       
    allocate(moments_1(ndime,nclus))
    allocate(moments_2(ndime,nclus))
    allocate(central_moments_2(ndime,nclus))
    
    call compute_moments(legro,nclus,concentration_measure,volume,centroid,moments_1,moments_2,central_moments_2) 
    
    do iclus = 1, nclus
       compactness_2(iclus) = 0.0_rp
       do idime = 1, ndime
          compactness_2(iclus) = compactness_2(iclus) + central_moments_2(idime,iclus)
       end do
       if (ndime == 2) then
          compactness_2(iclus) = ((volume(iclus)**2.0_rp) / (2.0_rp * pi)) / compactness_2(iclus)
       else if (ndime == 3) then
          compactness_2(iclus) = (((3.0_rp * volume(iclus))**(5.0_rp/3.0_rp)) / (5.0_rp * (4.0_rp*pi)**(2.0_rp/3.0_rp))) / compactness_2(iclus)
       end if
    end do
    
    deallocate(moments_1)
    deallocate(moments_2)
    deallocate(central_moments_2)
    
  end subroutine compute_compactness_2


  !-----------------------------------------------------------------------
  !> @{
  !> @author  margarida moragues
  !> @date    2020-03-02
  !> @brief   Moments of inertia
  !> @details Moments of inertia:
  !> This program computes the moments of inertia (or geometrical moments) of a group of clusters (or subdomains).
  !> Note that for any element, ielem, legro(ielem) is the cluster (or subdomain) 
  !> number to which ielem belongs. 
  !>
  !> The moments of inertia with respect to the origin of coordinates of a given shape/subdomain S are:
  !>
  !>    m_p_q_r = integral_S (x^p * y^q * z^r) dS.
  !> 
  !> The centroid of S is:
  !>
  !>    C = (m_1_0_0(S), m_0_1_0(S), m_0_0_1(S)) / Volume(S).
  !>
  !> The central moments of inertia are:
  !>
  !>    Mu_p_q_r = integral_S ((x-bar(x))^p * (y-bar(y))^q * (z-bar(z))^r) dS,
  !>
  !> where C=(bar(x),bar(y),bar(z)) is the centroid of S.
  !> After the Parallel Axis Theorem the central moments can be computed as:
  !>
  !>    Mu_2_0_0(S) = m_2_0_0(S) - Volume(S) * bar(x)^2
  !>    Mu_0_2_0(S) = m_0_2_0(S) - Volume(S) * bar(y)^2
  !>    Mu_0_0_2(S) = m_0_0_2(S) - Volume(S) * bar(z)^2.
  !>
  !> References:
  !> C. Martínez-Ortiz (2010) 2D and 3D shape descriptors. PhD Thesis. 
  !> J. Zunic et al. (2010) A hu moment invariant as a shape circularity measure. 
  !>        Pattern Recognition, 43(1):47?~@~S57. 21, 25
  !> 
  !> @}
  !-----------------------------------------------------------------------
  
  subroutine compute_moments(legro,nclus,concentration_measure,volume,centroid,moments_1,moments_2,central_moments_2)

    use def_domain
    use def_master
    use mod_func
    use mod_integrals

    implicit none
    integer(ip),          pointer,           intent(in)      :: legro(:)
    integer(ip),                             intent(in)      :: nclus
    real(rp),             pointer,           intent(in)      :: concentration_measure(:)
    real(rp),             pointer,           intent(inout)   :: volume(:)
    real(rp),             pointer,           intent(inout)   :: centroid(:,:) 
    real(rp),             allocatable,       intent(inout)   :: moments_1(:,:), moments_2(:,:), central_moments_2(:,:)
    integer(ip)                                              :: n_fields, n_integrals
    integer(ip)                                              :: idime, i_integral, ipoin, iclus
    type(func_ptr),       allocatable                        :: array_func(:,:)
    type(field_arrays),   allocatable                        :: array_fields(:)
    real(rp),             allocatable                        :: integrals(:,:)

    n_fields = ndime+1
    n_integrals = 2*ndime+1
        
    allocate(integrals(n_integrals,nclus))
    allocate(array_fields(n_fields))
    allocate(array_func(n_integrals,n_fields))

    call func_initialization(array_func)
    do i_integral = 1, n_integrals
       array_func(i_integral,1) % f => func_identity
    end do
    do idime = 1, ndime
       array_func(1+idime,1+idime) % f => func_identity
       array_func(1+ndime+idime,1+idime) % f => func_square
    end do

    do idime = 1, n_fields
       allocate(array_fields(idime)%a(1,npoin))
    end do
    do ipoin = 1,npoin
       array_fields(1) % a(1,ipoin) = concentration_measure(ipoin)
       do idime = 1, ndime
          array_fields(1+idime) % a(1,ipoin) = coord(idime,ipoin)
       end do
    end do

    call integrals_volume(array_fields,array_func,legro,integrals)

    do iclus = 1, nclus
       volume(iclus) = integrals(1,iclus)
       do idime = 1, ndime
          moments_1(idime,iclus) = integrals(1+idime,iclus)
          centroid(idime,iclus)  = moments_1(idime,iclus) / volume(iclus)
          moments_2(idime,iclus) = integrals(1+ndime+idime,iclus)
          central_moments_2(idime,iclus) = moments_2(idime,iclus) - volume(iclus) * centroid(idime,iclus)**2
       end do
    end do

    do idime = 1, n_fields
       deallocate(array_fields(idime)%a)
    end do
     
    deallocate(array_fields)
    deallocate(integrals)
    deallocate(array_func)

  end subroutine compute_moments

end module mod_shape_descriptors
