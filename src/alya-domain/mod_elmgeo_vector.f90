!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!>
!> @defgroup Elemental_Geometric_Toolbox
!> ToolBox for elemental and general geometrical operations
!> @{
!> @file    mod_elmgeo.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for elements
!> @details Different functions useful in finite element implementations
!>
!------------------------------------------------------------------------

module mod_elmgeo_vector

#include "def_vector_size.inc"
  use def_vectorization
  use def_kintyp_basic, only : ip,rp
  use mod_elmgeo,       only : elem_typ
  use mod_elmgeo,       only : elmgeo_cartesian_derivatives_jacobian
  use mod_elmgeo,       only : elmgeo_element_characteristic_length
  implicit none
  private

  real(rp), parameter :: epsilgeo_div = epsilon(1.0_rp) !epsilgeo used to avoid divisions by zero

  type, extends(elem_typ) :: elem_typ_vector
  end type elem_typ_vector

  interface elmgeo_cartesian_derivatives_jacobian
     module procedure elmgeo_cartesian_derivatives_jacobian_vectorized
  end interface elmgeo_cartesian_derivatives_jacobian

  interface elmgeo_element_characteristic_length
     module procedure elmgeo_element_characteristic_length_vectorized
  end interface elmgeo_element_characteristic_length
  
  public :: elmgeo_cartesian_derivatives_jacobian
  public :: elmgeo_cartesian_derivatives_jacobian_vectorized_cpu
  public :: elmgeo_element_characteristic_length
  
#ifdef OPENACCHHH
#define DEF_VECT ivect
#else
#define DEF_VECT 1:VECTOR_SIZE
#endif
#define DEF_VECT_CPU 1:VECTOR_SIZE_CPU

contains

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    15/09/2016
  !> @brief   Computes the Jacobian, Jacobian determinant and
  !>          Cartesian derivatives of shape function of an element
  !> @details The jacobian is
  !>                             _           _
  !>                            | dx/ds dx/dt |                t
  !>          Jacobian: XJACM = |             | = ELCOD * DERIV
  !>                            |_dy/ds dy/dt_|
  !>          with
  !>                   _        _             _                    _
  !>                  | x1 x2 x3 |           | dN1/ds dN2/ds dN3/ds |
  !>          ELCOD = |          |,  DERIV = |                      |
  !>                  |_y1 y2 y3_|           |_dN1/dt dN2/dt dN3/dt_|
  !>
  !>          => Jacobian determinant: GPDET = det(XJACM)
  !>
  !>          P1 Element in 2D
  !>          ----------------
  !>
  !>          gpdet=(-x1+x2)*(-y1+y3)-(-y1+y2)*(-x1+x3)
  !>                         _                       _
  !>                        | -y3+y2  -y1+y3   y1-y2  |
  !>          gpcar=1/gpdet |                         |
  !>                        |_ x3-x2   x1-x3  -x1+x2 _|
  !>
  !>          P1 Element in 3D
  !>          ----------------
  !>                 _                     _
  !>                |  x2-x1  x3-x1  x4-x1  |
  !>          xjacm=|  y2-y1  y3-y1  y4-y1  |
  !>                |_ z2-z1  z3-z1  z4-z1 _|
  !>
  !>                     -1
  !>          xjaci=xjacm
  !>                 _                             _     _          _
  !>                |  dN1/ds dN2/ds dN3/ds dN4/ds  |   |  -1 1 0 0  |
  !>          deriv=|  dN1/dt dN2/dt dN3/dt_dN4/dt  | = |  -1 0 1 0  |
  !>                |_ dN1/dz dN2/dz dN3/dz_dN4/dz _|   |_ -1 0 0 1 _|
  !>
  !>                     t
  !>          cartd=xjaci *deriv
  !>
  !-----------------------------------------------------------------------  

  pure subroutine elmgeo_cartesian_derivatives_jacobian_vectorized(VECTOR_DIM,ndime,mnode,pnode,pgaus,elcod,deriv,xjaci,gpcar,gpdet)

    integer(ip),           intent(in)  :: VECTOR_DIM
    integer(ip),           intent(in)  :: ndime
    integer(ip),           intent(in)  :: mnode
    integer(ip),           intent(in)  :: pnode
    integer(ip),           intent(in)  :: pgaus
    real(rp),              intent(in)  :: elcod(VECTOR_SIZE,ndime,pnode)
    real(rp),              intent(in)  :: deriv(ndime,pnode,pgaus)
    real(rp),              intent(out) :: xjaci(VECTOR_SIZE,ndime,ndime,pgaus)
    real(rp),              intent(out) :: gpcar(VECTOR_SIZE,ndime,mnode,pgaus)
    real(rp),              intent(out) :: gpdet(VECTOR_SIZE,pgaus)
#ifdef OPENACCHHH
    integer(ip)                        :: ivect
#endif
    integer(ip)                        :: j,k,idime,inode,igaus,jdime
    real(rp)                           :: xjacm(VECTOR_SIZE,ndime,ndime)
    real(rp)                           :: t1(VECTOR_SIZE)
    real(rp)                           :: t2(VECTOR_SIZE)
    real(rp)                           :: t3(VECTOR_SIZE)
    real(rp)                           :: denom(VECTOR_SIZE)

    if( ndime == 2 .and. pnode == 3 ) then
       !
       ! 2D P1 element
       !
#ifdef OPENACCHHH
       
       !$acc enter data create(denom)

       !$acc parallel loop gang vector default(present)
       do ivect = 1,VECTOR_SIZE
#endif       
          gpdet(DEF_VECT,1)     = (    -elcod(DEF_VECT,1,1) + elcod(DEF_VECT,1,2)) * (-elcod(DEF_VECT,2,1) + elcod(DEF_VECT,2,3) ) &
               &                    -( -elcod(DEF_VECT,2,1) + elcod(DEF_VECT,2,2)) * (-elcod(DEF_VECT,1,1) + elcod(DEF_VECT,1,3) )

          denom(DEF_VECT)       = 1.0_rp / (sign(1.0_rp,gpdet(DEF_VECT,1))*max(abs(gpdet(DEF_VECT,1)),epsilgeo_div))

          gpcar(DEF_VECT,1,1,1) = ( -elcod(DEF_VECT,2,3) + elcod(DEF_VECT,2,2) ) * denom(DEF_VECT)
          gpcar(DEF_VECT,1,2,1) = ( -elcod(DEF_VECT,2,1) + elcod(DEF_VECT,2,3) ) * denom(DEF_VECT)
          gpcar(DEF_VECT,1,3,1) = (  elcod(DEF_VECT,2,1) - elcod(DEF_VECT,2,2) ) * denom(DEF_VECT)
          gpcar(DEF_VECT,2,1,1) = (  elcod(DEF_VECT,1,3) - elcod(DEF_VECT,1,2) ) * denom(DEF_VECT)
          gpcar(DEF_VECT,2,2,1) = (  elcod(DEF_VECT,1,1) - elcod(DEF_VECT,1,3) ) * denom(DEF_VECT)
          gpcar(DEF_VECT,2,3,1) = ( -elcod(DEF_VECT,1,1) + elcod(DEF_VECT,1,2) ) * denom(DEF_VECT)

          do igaus = 2,pgaus
             gpdet(DEF_VECT,igaus) = gpdet(DEF_VECT,1)
             do jdime = 1,2
                do idime = 1,2
                   xjaci(DEF_VECT,idime,jdime,igaus) = xjaci(DEF_VECT,idime,jdime,1)
                end do
             end do
             do inode = 1,3
                do idime = 1,2
                   gpcar(DEF_VECT,idime,inode,igaus) = gpcar(DEF_VECT,idime,inode,1)
                end do
             end do
          end do
#ifdef OPENACCHHH
       end do
      !$acc end parallel loop

      !$acc exit data delete(denom) 

#endif       

    else if( ndime == 3 .and. pnode == 4 ) then
       !
       ! 3D P1 element
       !
#ifdef OPENACCHHH
       !$acc enter data create(t1, t2, t3, denom) 

       !$acc parallel loop gang vector default(present)
       do ivect = 1,VECTOR_SIZE
#endif       
          gpcar(DEF_VECT,1,1,1) =  elcod(DEF_VECT,1,2)   - elcod(DEF_VECT,1,1)
          gpcar(DEF_VECT,1,2,1) =  elcod(DEF_VECT,1,3)   - elcod(DEF_VECT,1,1)
          gpcar(DEF_VECT,1,3,1) =  elcod(DEF_VECT,1,4)   - elcod(DEF_VECT,1,1)
          gpcar(DEF_VECT,2,1,1) =  elcod(DEF_VECT,2,2)   - elcod(DEF_VECT,2,1)
          gpcar(DEF_VECT,2,2,1) =  elcod(DEF_VECT,2,3)   - elcod(DEF_VECT,2,1)
          gpcar(DEF_VECT,2,3,1) =  elcod(DEF_VECT,2,4)   - elcod(DEF_VECT,2,1)
          gpcar(DEF_VECT,3,1,1) =  elcod(DEF_VECT,3,2)   - elcod(DEF_VECT,3,1)
          gpcar(DEF_VECT,3,2,1) =  elcod(DEF_VECT,3,3)   - elcod(DEF_VECT,3,1)
          gpcar(DEF_VECT,3,3,1) =  elcod(DEF_VECT,3,4)   - elcod(DEF_VECT,3,1)
          t1(DEF_VECT)          =  gpcar(DEF_VECT,2,2,1) * gpcar(DEF_VECT,3,3,1) - gpcar(DEF_VECT,3,2,1) * gpcar(DEF_VECT,2,3,1)
          t2(DEF_VECT)          = -gpcar(DEF_VECT,2,1,1) * gpcar(DEF_VECT,3,3,1) + gpcar(DEF_VECT,3,1,1) * gpcar(DEF_VECT,2,3,1)
          t3(DEF_VECT)          =  gpcar(DEF_VECT,2,1,1) * gpcar(DEF_VECT,3,2,1) - gpcar(DEF_VECT,3,1,1) * gpcar(DEF_VECT,2,2,1)
          gpdet(DEF_VECT,1)     =  gpcar(DEF_VECT,1,1,1) * t1(DEF_VECT) + gpcar(DEF_VECT,1,2,1) * t2(DEF_VECT) + gpcar(DEF_VECT,1,3,1) * t3(DEF_VECT)

          denom(DEF_VECT)       =  1.0_rp / (sign(1.0_rp,gpdet(DEF_VECT,1))*max(abs(gpdet(DEF_VECT,1)),epsilgeo_div))

          xjaci(DEF_VECT,1,1,1) =  t1(DEF_VECT) * denom(DEF_VECT)
          xjaci(DEF_VECT,2,1,1) =  t2(DEF_VECT) * denom(DEF_VECT)
          xjaci(DEF_VECT,3,1,1) =  t3(DEF_VECT) * denom(DEF_VECT)
          xjaci(DEF_VECT,2,2,1) = ( gpcar(DEF_VECT,1,1,1) * gpcar(DEF_VECT,3,3,1) - gpcar(DEF_VECT,3,1,1) * gpcar(DEF_VECT,1,3,1)) * denom(DEF_VECT)
          xjaci(DEF_VECT,3,2,1) = (-gpcar(DEF_VECT,1,1,1) * gpcar(DEF_VECT,3,2,1) + gpcar(DEF_VECT,1,2,1) * gpcar(DEF_VECT,3,1,1)) * denom(DEF_VECT)
          xjaci(DEF_VECT,3,3,1) = ( gpcar(DEF_VECT,1,1,1) * gpcar(DEF_VECT,2,2,1) - gpcar(DEF_VECT,2,1,1) * gpcar(DEF_VECT,1,2,1)) * denom(DEF_VECT)
          xjaci(DEF_VECT,1,2,1) = (-gpcar(DEF_VECT,1,2,1) * gpcar(DEF_VECT,3,3,1) + gpcar(DEF_VECT,3,2,1) * gpcar(DEF_VECT,1,3,1)) * denom(DEF_VECT)
          xjaci(DEF_VECT,1,3,1) = ( gpcar(DEF_VECT,1,2,1) * gpcar(DEF_VECT,2,3,1) - gpcar(DEF_VECT,2,2,1) * gpcar(DEF_VECT,1,3,1)) * denom(DEF_VECT)
          xjaci(DEF_VECT,2,3,1) = (-gpcar(DEF_VECT,1,1,1) * gpcar(DEF_VECT,2,3,1) + gpcar(DEF_VECT,2,1,1) * gpcar(DEF_VECT,1,3,1)) * denom(DEF_VECT)

          gpcar(DEF_VECT,1,1,1) = -xjaci(DEF_VECT,1,1,1) - xjaci(DEF_VECT,2,1,1) - xjaci(DEF_VECT,3,1,1)
          gpcar(DEF_VECT,1,2,1) =  xjaci(DEF_VECT,1,1,1)
          gpcar(DEF_VECT,1,3,1) =  xjaci(DEF_VECT,2,1,1)
          gpcar(DEF_VECT,1,4,1) =  xjaci(DEF_VECT,3,1,1)
          gpcar(DEF_VECT,2,1,1) = -xjaci(DEF_VECT,1,2,1) - xjaci(DEF_VECT,2,2,1) - xjaci(DEF_VECT,3,2,1)
          gpcar(DEF_VECT,2,2,1) =  xjaci(DEF_VECT,1,2,1)
          gpcar(DEF_VECT,2,3,1) =  xjaci(DEF_VECT,2,2,1)
          gpcar(DEF_VECT,2,4,1) =  xjaci(DEF_VECT,3,2,1)
          gpcar(DEF_VECT,3,1,1) = -xjaci(DEF_VECT,1,3,1) - xjaci(DEF_VECT,2,3,1) - xjaci(DEF_VECT,3,3,1)
          gpcar(DEF_VECT,3,2,1) =  xjaci(DEF_VECT,1,3,1)
          gpcar(DEF_VECT,3,3,1) =  xjaci(DEF_VECT,2,3,1)
          gpcar(DEF_VECT,3,4,1) =  xjaci(DEF_VECT,3,3,1)

          do igaus = 2,pgaus
             gpdet(DEF_VECT,igaus) = gpdet(DEF_VECT,1)
             do jdime = 1,3
                do idime = 1,3
                   xjaci(DEF_VECT,idime,jdime,igaus) = xjaci(DEF_VECT,idime,jdime,1)
                end do
             end do
             do inode = 1,4
                do idime = 1,3
                   gpcar(DEF_VECT,idime,inode,igaus) = gpcar(DEF_VECT,idime,inode,1)
                end do
             end do
          end do
#ifdef OPENACCHHH
       end do
      !$acc end parallel loop

      !$acc exit data delete(denom, t1, t2 ,t3 ) 

#endif       

    else if ( ndime == 1 ) then
       !
       ! 1D
       !
#ifdef OPENACCHHH
       !$acc enter data create( xjacm)

       !$acc parallel loop gang vector default(present)
       do ivect = 1,VECTOR_SIZE
#endif       
          do igaus = 1,pgaus
             xjacm(DEF_VECT,1,1) = 0.0_rp
             do k = 1,pnode
                xjacm(DEF_VECT,1,1) = xjacm(DEF_VECT,1,1) + elcod(DEF_VECT,1,k) * deriv(1,k,igaus)
             end do
             gpdet(DEF_VECT,igaus) = xjacm(DEF_VECT,1,1)
             xjaci(DEF_VECT,1,1,igaus) = 1.0_rp / xjacm(DEF_VECT,1,1)
             do j = 1,pnode
                gpcar(DEF_VECT,1,j,igaus) = xjaci(DEF_VECT,1,1,igaus) * deriv(1,j,igaus)
             end do
          end do
#ifdef OPENACCHHH
       end do
      !$acc end parallel loop

      !$acc exit data delete(xjacm)
      
#endif       

    else if ( ndime == 2 ) then
       !
       ! 2D
       !
#ifdef OPENACCHHH
       !$acc enter data create(xjacm, denom)

       !$acc parallel loop gang vector default(present)
       do ivect = 1,VECTOR_SIZE
#endif       
          do igaus = 1,pgaus

             xjacm(DEF_VECT,1,1) = 0.0_rp
             xjacm(DEF_VECT,1,2) = 0.0_rp
             xjacm(DEF_VECT,2,1) = 0.0_rp
             xjacm(DEF_VECT,2,2) = 0.0_rp
             do k = 1,pnode
                xjacm(DEF_VECT,1,1) = xjacm(DEF_VECT,1,1) + elcod(DEF_VECT,1,k) * deriv(1,k,igaus)
                xjacm(DEF_VECT,1,2) = xjacm(DEF_VECT,1,2) + elcod(DEF_VECT,1,k) * deriv(2,k,igaus)
                xjacm(DEF_VECT,2,1) = xjacm(DEF_VECT,2,1) + elcod(DEF_VECT,2,k) * deriv(1,k,igaus)
                xjacm(DEF_VECT,2,2) = xjacm(DEF_VECT,2,2) + elcod(DEF_VECT,2,k) * deriv(2,k,igaus)
             end do

             gpdet(DEF_VECT,igaus) = xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,2,2) - xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,1,2)

             denom(DEF_VECT)       = 1.0_rp / (sign(1.0_rp,gpdet(DEF_VECT,igaus))*max(abs(gpdet(DEF_VECT,igaus)),epsilgeo_div))

             xjaci(DEF_VECT,1,1,igaus) =  xjacm(DEF_VECT,2,2) * denom(DEF_VECT)
             xjaci(DEF_VECT,2,2,igaus) =  xjacm(DEF_VECT,1,1) * denom(DEF_VECT)
             xjaci(DEF_VECT,2,1,igaus) = -xjacm(DEF_VECT,2,1) * denom(DEF_VECT)
             xjaci(DEF_VECT,1,2,igaus) = -xjacm(DEF_VECT,1,2) * denom(DEF_VECT)

             do j = 1, pnode
                gpcar(DEF_VECT,1,j,igaus) =   xjaci(DEF_VECT,1,1,igaus) * deriv(1,j,igaus) &
                     &                           + xjaci(DEF_VECT,2,1,igaus) * deriv(2,j,igaus)

                gpcar(DEF_VECT,2,j,igaus) =   xjaci(DEF_VECT,1,2,igaus) * deriv(1,j,igaus) &
                     &                           + xjaci(DEF_VECT,2,2,igaus) * deriv(2,j,igaus)
             end do

          end do
#ifdef OPENACCHHH
       end do
      !$acc end parallel loop

      !$acc exit data delete(xjacm, denom ) 

#endif       

    else if ( ndime == 3 ) then
       !
       ! 3D
       !
       ! xjacm = elcod * deriv^t
       ! xjaci = xjacm^-1
       ! gpcar = xjaci^t * deriv
       !
#ifdef OPENACCHHH
       !$acc enter data create(xjacm, denom, t1, t2, t3)

       !$acc parallel loop gang vector default(present)
       do ivect = 1,VECTOR_SIZE
#endif       
          do igaus = 1,pgaus

             xjacm(DEF_VECT,1:3,1:3) = 0.0_rp
             do k = 1,pnode
                xjacm(DEF_VECT,1,1) = xjacm(DEF_VECT,1,1) + elcod(DEF_VECT,1,k) * deriv(1,k,igaus)
                xjacm(DEF_VECT,1,2) = xjacm(DEF_VECT,1,2) + elcod(DEF_VECT,1,k) * deriv(2,k,igaus)
                xjacm(DEF_VECT,1,3) = xjacm(DEF_VECT,1,3) + elcod(DEF_VECT,1,k) * deriv(3,k,igaus)
                xjacm(DEF_VECT,2,1) = xjacm(DEF_VECT,2,1) + elcod(DEF_VECT,2,k) * deriv(1,k,igaus)
                xjacm(DEF_VECT,2,2) = xjacm(DEF_VECT,2,2) + elcod(DEF_VECT,2,k) * deriv(2,k,igaus)
                xjacm(DEF_VECT,2,3) = xjacm(DEF_VECT,2,3) + elcod(DEF_VECT,2,k) * deriv(3,k,igaus)
                xjacm(DEF_VECT,3,1) = xjacm(DEF_VECT,3,1) + elcod(DEF_VECT,3,k) * deriv(1,k,igaus)
                xjacm(DEF_VECT,3,2) = xjacm(DEF_VECT,3,2) + elcod(DEF_VECT,3,k) * deriv(2,k,igaus)
                xjacm(DEF_VECT,3,3) = xjacm(DEF_VECT,3,3) + elcod(DEF_VECT,3,k) * deriv(3,k,igaus)
             end do

             t1(DEF_VECT)              =  xjacm(DEF_VECT,2,2) * xjacm(DEF_VECT,3,3) - xjacm(DEF_VECT,3,2) * xjacm(DEF_VECT,2,3)
             t2(DEF_VECT)              = -xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,3,3) + xjacm(DEF_VECT,3,1) * xjacm(DEF_VECT,2,3)
             t3(DEF_VECT)              =  xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,3,2) - xjacm(DEF_VECT,3,1) * xjacm(DEF_VECT,2,2)
             gpdet(DEF_VECT,igaus)     =  xjacm(DEF_VECT,1,1) * t1(DEF_VECT) + xjacm(DEF_VECT,1,2) * t2(DEF_VECT) + xjacm(DEF_VECT,1,3) * t3(DEF_VECT)

             denom(DEF_VECT)           = 1.0_rp / (sign(1.0_rp,gpdet(DEF_VECT,igaus))*max(abs(gpdet(DEF_VECT,igaus)),epsilgeo_div))
             xjaci(DEF_VECT,1,1,igaus) = t1(DEF_VECT) * denom(DEF_VECT)
             xjaci(DEF_VECT,2,1,igaus) = t2(DEF_VECT) * denom(DEF_VECT)
             xjaci(DEF_VECT,3,1,igaus) = t3(DEF_VECT) * denom(DEF_VECT)
             xjaci(DEF_VECT,2,2,igaus) = ( xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,3,3) - xjacm(DEF_VECT,3,1) * xjacm(DEF_VECT,1,3)) * denom(DEF_VECT)
             xjaci(DEF_VECT,3,2,igaus) = (-xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,3,2) + xjacm(DEF_VECT,1,2) * xjacm(DEF_VECT,3,1)) * denom(DEF_VECT)
             xjaci(DEF_VECT,3,3,igaus) = ( xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,2,2) - xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,1,2)) * denom(DEF_VECT)
             xjaci(DEF_VECT,1,2,igaus) = (-xjacm(DEF_VECT,1,2) * xjacm(DEF_VECT,3,3) + xjacm(DEF_VECT,3,2) * xjacm(DEF_VECT,1,3)) * denom(DEF_VECT)
             xjaci(DEF_VECT,1,3,igaus) = ( xjacm(DEF_VECT,1,2) * xjacm(DEF_VECT,2,3) - xjacm(DEF_VECT,2,2) * xjacm(DEF_VECT,1,3)) * denom(DEF_VECT)
             xjaci(DEF_VECT,2,3,igaus) = (-xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,2,3) + xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,1,3)) * denom(DEF_VECT)

             do j = 1, pnode
                gpcar(DEF_VECT,1,j,igaus) =   xjaci(DEF_VECT,1,1,igaus) * deriv(1,j,igaus) &
                     &                           + xjaci(DEF_VECT,2,1,igaus) * deriv(2,j,igaus) &
                     &                           + xjaci(DEF_VECT,3,1,igaus) * deriv(3,j,igaus)

                gpcar(DEF_VECT,2,j,igaus) =   xjaci(DEF_VECT,1,2,igaus) * deriv(1,j,igaus) &
                     &                           + xjaci(DEF_VECT,2,2,igaus) * deriv(2,j,igaus) &
                     &                           + xjaci(DEF_VECT,3,2,igaus) * deriv(3,j,igaus)

                gpcar(DEF_VECT,3,j,igaus) =   xjaci(DEF_VECT,1,3,igaus) * deriv(1,j,igaus) &
                     &                           + xjaci(DEF_VECT,2,3,igaus) * deriv(2,j,igaus) &
                     &                           + xjaci(DEF_VECT,3,3,igaus) * deriv(3,j,igaus)
             end do

          end do
#ifdef OPENACCHHH
       end do
      !$acc end parallel loop


      !$acc exit data delete(xjacm, denom, &
      !$acc                  t1, t2, t3 ) 


#endif       
    end if

  end subroutine elmgeo_cartesian_derivatives_jacobian_vectorized

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    15/09/2016
  !> @brief   Computes the Jacobian, Jacobian determinant and
  !>          Cartesian derivatives of shape function of an element
  !> @details The jacobian is
  !>                             _           _
  !>                            | dx/ds dx/dt |                t
  !>          Jacobian: XJACM = |             | = ELCOD * DERIV
  !>                            |_dy/ds dy/dt_|
  !>          with
  !>                   _        _             _                    _
  !>                  | x1 x2 x3 |           | dN1/ds dN2/ds dN3/ds |
  !>          ELCOD = |          |,  DERIV = |                      |
  !>                  |_y1 y2 y3_|           |_dN1/dt dN2/dt dN3/dt_|
  !>
  !>          => Jacobian determinant: GPDET = det(XJACM)
  !>
  !>          P1 Element in 2D
  !>          ----------------
  !>
  !>          gpdet=(-x1+x2)*(-y1+y3)-(-y1+y2)*(-x1+x3)
  !>                         _                       _
  !>                        | -y3+y2  -y1+y3   y1-y2  |
  !>          gpcar=1/gpdet |                         |
  !>                        |_ x3-x2   x1-x3  -x1+x2 _|
  !>
  !>          P1 Element in 3D
  !>          ----------------
  !>                 _                     _
  !>                |  x2-x1  x3-x1  x4-x1  |
  !>          xjacm=|  y2-y1  y3-y1  y4-y1  |
  !>                |_ z2-z1  z3-z1  z4-z1 _|
  !>
  !>                     -1
  !>          xjaci=xjacm
  !>                 _                             _     _          _
  !>                |  dN1/ds dN2/ds dN3/ds dN4/ds  |   |  -1 1 0 0  |
  !>          deriv=|  dN1/dt dN2/dt dN3/dt_dN4/dt  | = |  -1 0 1 0  |
  !>                |_ dN1/dz dN2/dz dN3/dz_dN4/dz _|   |_ -1 0 0 1 _|
  !>
  !>                     t
  !>          cartd=xjaci *deriv
  !>
  !-----------------------------------------------------------------------  

  pure subroutine elmgeo_cartesian_derivatives_jacobian_vectorized_cpu(ndime,mnode,pnode,pgaus,elcod,deriv,xjaci,gpcar,gpdet)

    integer(ip),           intent(in)    :: ndime
    integer(ip),           intent(in)    :: mnode
    integer(ip),           intent(in)    :: pnode
    integer(ip),           intent(in)    :: pgaus
    real(rp),              intent(in)    :: elcod(VECTOR_SIZE_CPU,ndime,pnode)
    real(rp),              intent(in)    :: deriv(ndime,pnode,pgaus)
    real(rp),              intent(inout) :: xjaci(VECTOR_SIZE_CPU,ndime,ndime,pgaus)
    real(rp),              intent(out)   :: gpcar(VECTOR_SIZE_CPU,ndime,mnode,pgaus)
    real(rp),              intent(out)   :: gpdet(VECTOR_SIZE_CPU,pgaus)
    integer(ip)                          :: j,k,idime,inode,igaus,jdime
    real(rp)                             :: xjacm(VECTOR_SIZE_CPU,ndime,ndime)
    real(rp)                             :: t1(VECTOR_SIZE_CPU)
    real(rp)                             :: t2(VECTOR_SIZE_CPU)
    real(rp)                             :: t3(VECTOR_SIZE_CPU)
    real(rp)                             :: denom(VECTOR_SIZE_CPU)

    if( ndime == 2 .and. pnode == 3 ) then
       !
       ! 2D P1 element
       !
       gpdet(DEF_VECT_CPU,1)     = (    -elcod(DEF_VECT_CPU,1,1) + elcod(DEF_VECT_CPU,1,2)) * (-elcod(DEF_VECT_CPU,2,1) + elcod(DEF_VECT_CPU,2,3) ) &
            &                    -( -elcod(DEF_VECT_CPU,2,1) + elcod(DEF_VECT_CPU,2,2)) * (-elcod(DEF_VECT_CPU,1,1) + elcod(DEF_VECT_CPU,1,3) )

       denom(DEF_VECT_CPU)       = 1.0_rp / (sign(1.0_rp,gpdet(DEF_VECT_CPU,1))*max(abs(gpdet(DEF_VECT_CPU,1)),epsilgeo_div))

       gpcar(DEF_VECT_CPU,1,1,1) = ( -elcod(DEF_VECT_CPU,2,3) + elcod(DEF_VECT_CPU,2,2) ) * denom(DEF_VECT_CPU)
       gpcar(DEF_VECT_CPU,1,2,1) = ( -elcod(DEF_VECT_CPU,2,1) + elcod(DEF_VECT_CPU,2,3) ) * denom(DEF_VECT_CPU)
       gpcar(DEF_VECT_CPU,1,3,1) = (  elcod(DEF_VECT_CPU,2,1) - elcod(DEF_VECT_CPU,2,2) ) * denom(DEF_VECT_CPU)
       gpcar(DEF_VECT_CPU,2,1,1) = (  elcod(DEF_VECT_CPU,1,3) - elcod(DEF_VECT_CPU,1,2) ) * denom(DEF_VECT_CPU)
       gpcar(DEF_VECT_CPU,2,2,1) = (  elcod(DEF_VECT_CPU,1,1) - elcod(DEF_VECT_CPU,1,3) ) * denom(DEF_VECT_CPU)
       gpcar(DEF_VECT_CPU,2,3,1) = ( -elcod(DEF_VECT_CPU,1,1) + elcod(DEF_VECT_CPU,1,2) ) * denom(DEF_VECT_CPU)

       !TODO: here xjaci is used uninitialized potentially
       do igaus = 2,pgaus
          gpdet(DEF_VECT_CPU,igaus) = gpdet(DEF_VECT_CPU,1)
          do jdime = 1,2
             do idime = 1,2
                xjaci(DEF_VECT_CPU,idime,jdime,igaus) = xjaci(DEF_VECT_CPU,idime,jdime,1)
             end do
          end do
          do inode = 1,3
             do idime = 1,2
                gpcar(DEF_VECT_CPU,idime,inode,igaus) = gpcar(DEF_VECT_CPU,idime,inode,1)
             end do
          end do
       end do

    else if( ndime == 3 .and. pnode == 4 ) then
       !
       ! 3D P1 element
       !
       gpcar(DEF_VECT_CPU,1,1,1) =  elcod(DEF_VECT_CPU,1,2)   - elcod(DEF_VECT_CPU,1,1)
       gpcar(DEF_VECT_CPU,1,2,1) =  elcod(DEF_VECT_CPU,1,3)   - elcod(DEF_VECT_CPU,1,1)
       gpcar(DEF_VECT_CPU,1,3,1) =  elcod(DEF_VECT_CPU,1,4)   - elcod(DEF_VECT_CPU,1,1)
       gpcar(DEF_VECT_CPU,2,1,1) =  elcod(DEF_VECT_CPU,2,2)   - elcod(DEF_VECT_CPU,2,1)
       gpcar(DEF_VECT_CPU,2,2,1) =  elcod(DEF_VECT_CPU,2,3)   - elcod(DEF_VECT_CPU,2,1)
       gpcar(DEF_VECT_CPU,2,3,1) =  elcod(DEF_VECT_CPU,2,4)   - elcod(DEF_VECT_CPU,2,1)
       gpcar(DEF_VECT_CPU,3,1,1) =  elcod(DEF_VECT_CPU,3,2)   - elcod(DEF_VECT_CPU,3,1)
       gpcar(DEF_VECT_CPU,3,2,1) =  elcod(DEF_VECT_CPU,3,3)   - elcod(DEF_VECT_CPU,3,1)
       gpcar(DEF_VECT_CPU,3,3,1) =  elcod(DEF_VECT_CPU,3,4)   - elcod(DEF_VECT_CPU,3,1)
       t1(DEF_VECT_CPU)          =  gpcar(DEF_VECT_CPU,2,2,1) * gpcar(DEF_VECT_CPU,3,3,1) - gpcar(DEF_VECT_CPU,3,2,1) * gpcar(DEF_VECT_CPU,2,3,1)
       t2(DEF_VECT_CPU)          = -gpcar(DEF_VECT_CPU,2,1,1) * gpcar(DEF_VECT_CPU,3,3,1) + gpcar(DEF_VECT_CPU,3,1,1) * gpcar(DEF_VECT_CPU,2,3,1)
       t3(DEF_VECT_CPU)          =  gpcar(DEF_VECT_CPU,2,1,1) * gpcar(DEF_VECT_CPU,3,2,1) - gpcar(DEF_VECT_CPU,3,1,1) * gpcar(DEF_VECT_CPU,2,2,1)
       gpdet(DEF_VECT_CPU,1)     =  gpcar(DEF_VECT_CPU,1,1,1) * t1(DEF_VECT_CPU) + gpcar(DEF_VECT_CPU,1,2,1) * t2(DEF_VECT_CPU) + gpcar(DEF_VECT_CPU,1,3,1) * t3(DEF_VECT_CPU)

       denom(DEF_VECT_CPU)       =  1.0_rp / (sign(1.0_rp,gpdet(DEF_VECT_CPU,1))*max(abs(gpdet(DEF_VECT_CPU,1)),epsilgeo_div))

       xjaci(DEF_VECT_CPU,1,1,1) =  t1(DEF_VECT_CPU) * denom(DEF_VECT_CPU)
       xjaci(DEF_VECT_CPU,2,1,1) =  t2(DEF_VECT_CPU) * denom(DEF_VECT_CPU)
       xjaci(DEF_VECT_CPU,3,1,1) =  t3(DEF_VECT_CPU) * denom(DEF_VECT_CPU)
       xjaci(DEF_VECT_CPU,2,2,1) = ( gpcar(DEF_VECT_CPU,1,1,1) * gpcar(DEF_VECT_CPU,3,3,1) - gpcar(DEF_VECT_CPU,3,1,1) * gpcar(DEF_VECT_CPU,1,3,1)) * denom(DEF_VECT_CPU)
       xjaci(DEF_VECT_CPU,3,2,1) = (-gpcar(DEF_VECT_CPU,1,1,1) * gpcar(DEF_VECT_CPU,3,2,1) + gpcar(DEF_VECT_CPU,1,2,1) * gpcar(DEF_VECT_CPU,3,1,1)) * denom(DEF_VECT_CPU)
       xjaci(DEF_VECT_CPU,3,3,1) = ( gpcar(DEF_VECT_CPU,1,1,1) * gpcar(DEF_VECT_CPU,2,2,1) - gpcar(DEF_VECT_CPU,2,1,1) * gpcar(DEF_VECT_CPU,1,2,1)) * denom(DEF_VECT_CPU)
       xjaci(DEF_VECT_CPU,1,2,1) = (-gpcar(DEF_VECT_CPU,1,2,1) * gpcar(DEF_VECT_CPU,3,3,1) + gpcar(DEF_VECT_CPU,3,2,1) * gpcar(DEF_VECT_CPU,1,3,1)) * denom(DEF_VECT_CPU)
       xjaci(DEF_VECT_CPU,1,3,1) = ( gpcar(DEF_VECT_CPU,1,2,1) * gpcar(DEF_VECT_CPU,2,3,1) - gpcar(DEF_VECT_CPU,2,2,1) * gpcar(DEF_VECT_CPU,1,3,1)) * denom(DEF_VECT_CPU)
       xjaci(DEF_VECT_CPU,2,3,1) = (-gpcar(DEF_VECT_CPU,1,1,1) * gpcar(DEF_VECT_CPU,2,3,1) + gpcar(DEF_VECT_CPU,2,1,1) * gpcar(DEF_VECT_CPU,1,3,1)) * denom(DEF_VECT_CPU)

       gpcar(DEF_VECT_CPU,1,1,1) = -xjaci(DEF_VECT_CPU,1,1,1) - xjaci(DEF_VECT_CPU,2,1,1) - xjaci(DEF_VECT_CPU,3,1,1)
       gpcar(DEF_VECT_CPU,1,2,1) =  xjaci(DEF_VECT_CPU,1,1,1)
       gpcar(DEF_VECT_CPU,1,3,1) =  xjaci(DEF_VECT_CPU,2,1,1)
       gpcar(DEF_VECT_CPU,1,4,1) =  xjaci(DEF_VECT_CPU,3,1,1)
       gpcar(DEF_VECT_CPU,2,1,1) = -xjaci(DEF_VECT_CPU,1,2,1) - xjaci(DEF_VECT_CPU,2,2,1) - xjaci(DEF_VECT_CPU,3,2,1)
       gpcar(DEF_VECT_CPU,2,2,1) =  xjaci(DEF_VECT_CPU,1,2,1)
       gpcar(DEF_VECT_CPU,2,3,1) =  xjaci(DEF_VECT_CPU,2,2,1)
       gpcar(DEF_VECT_CPU,2,4,1) =  xjaci(DEF_VECT_CPU,3,2,1)
       gpcar(DEF_VECT_CPU,3,1,1) = -xjaci(DEF_VECT_CPU,1,3,1) - xjaci(DEF_VECT_CPU,2,3,1) - xjaci(DEF_VECT_CPU,3,3,1)
       gpcar(DEF_VECT_CPU,3,2,1) =  xjaci(DEF_VECT_CPU,1,3,1)
       gpcar(DEF_VECT_CPU,3,3,1) =  xjaci(DEF_VECT_CPU,2,3,1)
       gpcar(DEF_VECT_CPU,3,4,1) =  xjaci(DEF_VECT_CPU,3,3,1)

       do igaus = 2,pgaus
          gpdet(DEF_VECT_CPU,igaus) = gpdet(DEF_VECT_CPU,1)
          do jdime = 1,3
             do idime = 1,3
                xjaci(DEF_VECT_CPU,idime,jdime,igaus) = xjaci(DEF_VECT_CPU,idime,jdime,1)
             end do
          end do
          do inode = 1,4
             do idime = 1,3
                gpcar(DEF_VECT_CPU,idime,inode,igaus) = gpcar(DEF_VECT_CPU,idime,inode,1)
             end do
          end do
       end do

    else if ( ndime == 1 ) then
       !
       ! 1D
       !
       do igaus = 1,pgaus
          xjacm(DEF_VECT_CPU,1,1) = 0.0_rp
          do k = 1,pnode
             xjacm(DEF_VECT_CPU,1,1) = xjacm(DEF_VECT_CPU,1,1) + elcod(DEF_VECT_CPU,1,k) * deriv(1,k,igaus)
          end do
          gpdet(DEF_VECT_CPU,igaus) = xjacm(DEF_VECT_CPU,1,1)
          xjaci(DEF_VECT_CPU,1,1,igaus) = 1.0_rp / xjacm(DEF_VECT_CPU,1,1)
          do j = 1,pnode
             gpcar(DEF_VECT_CPU,1,j,igaus) = xjaci(DEF_VECT_CPU,1,1,igaus) * deriv(1,j,igaus)
          end do
       end do

    else if ( ndime == 2 ) then
       !
       ! 2D
       !
       do igaus = 1,pgaus

          xjacm(DEF_VECT_CPU,1,1) = 0.0_rp
          xjacm(DEF_VECT_CPU,1,2) = 0.0_rp
          xjacm(DEF_VECT_CPU,2,1) = 0.0_rp
          xjacm(DEF_VECT_CPU,2,2) = 0.0_rp
          do k = 1,pnode
             xjacm(DEF_VECT_CPU,1,1) = xjacm(DEF_VECT_CPU,1,1) + elcod(DEF_VECT_CPU,1,k) * deriv(1,k,igaus)
             xjacm(DEF_VECT_CPU,1,2) = xjacm(DEF_VECT_CPU,1,2) + elcod(DEF_VECT_CPU,1,k) * deriv(2,k,igaus)
             xjacm(DEF_VECT_CPU,2,1) = xjacm(DEF_VECT_CPU,2,1) + elcod(DEF_VECT_CPU,2,k) * deriv(1,k,igaus)
             xjacm(DEF_VECT_CPU,2,2) = xjacm(DEF_VECT_CPU,2,2) + elcod(DEF_VECT_CPU,2,k) * deriv(2,k,igaus)
          end do

          gpdet(DEF_VECT_CPU,igaus) = xjacm(DEF_VECT_CPU,1,1) * xjacm(DEF_VECT_CPU,2,2) - xjacm(DEF_VECT_CPU,2,1) * xjacm(DEF_VECT_CPU,1,2)

          denom(DEF_VECT_CPU)       = 1.0_rp / (sign(1.0_rp,gpdet(DEF_VECT_CPU,igaus))*max(abs(gpdet(DEF_VECT_CPU,igaus)),epsilgeo_div))

          xjaci(DEF_VECT_CPU,1,1,igaus) =  xjacm(DEF_VECT_CPU,2,2) * denom(DEF_VECT_CPU)
          xjaci(DEF_VECT_CPU,2,2,igaus) =  xjacm(DEF_VECT_CPU,1,1) * denom(DEF_VECT_CPU)
          xjaci(DEF_VECT_CPU,2,1,igaus) = -xjacm(DEF_VECT_CPU,2,1) * denom(DEF_VECT_CPU)
          xjaci(DEF_VECT_CPU,1,2,igaus) = -xjacm(DEF_VECT_CPU,1,2) * denom(DEF_VECT_CPU)

          do j = 1, pnode
             gpcar(DEF_VECT_CPU,1,j,igaus) =   xjaci(DEF_VECT_CPU,1,1,igaus) * deriv(1,j,igaus) &
                  &                           + xjaci(DEF_VECT_CPU,2,1,igaus) * deriv(2,j,igaus)

             gpcar(DEF_VECT_CPU,2,j,igaus) =   xjaci(DEF_VECT_CPU,1,2,igaus) * deriv(1,j,igaus) &
                  &                           + xjaci(DEF_VECT_CPU,2,2,igaus) * deriv(2,j,igaus)
          end do

       end do

    else if ( ndime == 3 ) then
       !
       ! 3D
       !
       ! xjacm = elcod * deriv^t
       ! xjaci = xjacm^-1
       ! gpcar = xjaci^t * deriv
       ! 
       do igaus = 1,pgaus

          xjacm(DEF_VECT_CPU,1:3,1:3) = 0.0_rp
          do k = 1,pnode
             xjacm(DEF_VECT_CPU,1,1) = xjacm(DEF_VECT_CPU,1,1) + elcod(DEF_VECT_CPU,1,k) * deriv(1,k,igaus)
             xjacm(DEF_VECT_CPU,1,2) = xjacm(DEF_VECT_CPU,1,2) + elcod(DEF_VECT_CPU,1,k) * deriv(2,k,igaus)
             xjacm(DEF_VECT_CPU,1,3) = xjacm(DEF_VECT_CPU,1,3) + elcod(DEF_VECT_CPU,1,k) * deriv(3,k,igaus)
             xjacm(DEF_VECT_CPU,2,1) = xjacm(DEF_VECT_CPU,2,1) + elcod(DEF_VECT_CPU,2,k) * deriv(1,k,igaus)
             xjacm(DEF_VECT_CPU,2,2) = xjacm(DEF_VECT_CPU,2,2) + elcod(DEF_VECT_CPU,2,k) * deriv(2,k,igaus)
             xjacm(DEF_VECT_CPU,2,3) = xjacm(DEF_VECT_CPU,2,3) + elcod(DEF_VECT_CPU,2,k) * deriv(3,k,igaus)
             xjacm(DEF_VECT_CPU,3,1) = xjacm(DEF_VECT_CPU,3,1) + elcod(DEF_VECT_CPU,3,k) * deriv(1,k,igaus)
             xjacm(DEF_VECT_CPU,3,2) = xjacm(DEF_VECT_CPU,3,2) + elcod(DEF_VECT_CPU,3,k) * deriv(2,k,igaus)
             xjacm(DEF_VECT_CPU,3,3) = xjacm(DEF_VECT_CPU,3,3) + elcod(DEF_VECT_CPU,3,k) * deriv(3,k,igaus)
          end do

          t1(DEF_VECT_CPU)              =  xjacm(DEF_VECT_CPU,2,2) * xjacm(DEF_VECT_CPU,3,3) - xjacm(DEF_VECT_CPU,3,2) * xjacm(DEF_VECT_CPU,2,3)
          t2(DEF_VECT_CPU)              = -xjacm(DEF_VECT_CPU,2,1) * xjacm(DEF_VECT_CPU,3,3) + xjacm(DEF_VECT_CPU,3,1) * xjacm(DEF_VECT_CPU,2,3)
          t3(DEF_VECT_CPU)              =  xjacm(DEF_VECT_CPU,2,1) * xjacm(DEF_VECT_CPU,3,2) - xjacm(DEF_VECT_CPU,3,1) * xjacm(DEF_VECT_CPU,2,2)
          gpdet(DEF_VECT_CPU,igaus)     =  xjacm(DEF_VECT_CPU,1,1) * t1(DEF_VECT_CPU) + xjacm(DEF_VECT_CPU,1,2) * t2(DEF_VECT_CPU) + xjacm(DEF_VECT_CPU,1,3) * t3(DEF_VECT_CPU)

          denom(DEF_VECT_CPU)           = 1.0_rp / (sign(1.0_rp,gpdet(DEF_VECT_CPU,igaus))*max(abs(gpdet(DEF_VECT_CPU,igaus)),epsilgeo_div))
          xjaci(DEF_VECT_CPU,1,1,igaus) = t1(DEF_VECT_CPU) * denom(DEF_VECT_CPU)
          xjaci(DEF_VECT_CPU,2,1,igaus) = t2(DEF_VECT_CPU) * denom(DEF_VECT_CPU)
          xjaci(DEF_VECT_CPU,3,1,igaus) = t3(DEF_VECT_CPU) * denom(DEF_VECT_CPU)
          xjaci(DEF_VECT_CPU,2,2,igaus) = ( xjacm(DEF_VECT_CPU,1,1) * xjacm(DEF_VECT_CPU,3,3) - xjacm(DEF_VECT_CPU,3,1) * xjacm(DEF_VECT_CPU,1,3)) * denom(DEF_VECT_CPU)
          xjaci(DEF_VECT_CPU,3,2,igaus) = (-xjacm(DEF_VECT_CPU,1,1) * xjacm(DEF_VECT_CPU,3,2) + xjacm(DEF_VECT_CPU,1,2) * xjacm(DEF_VECT_CPU,3,1)) * denom(DEF_VECT_CPU)
          xjaci(DEF_VECT_CPU,3,3,igaus) = ( xjacm(DEF_VECT_CPU,1,1) * xjacm(DEF_VECT_CPU,2,2) - xjacm(DEF_VECT_CPU,2,1) * xjacm(DEF_VECT_CPU,1,2)) * denom(DEF_VECT_CPU)
          xjaci(DEF_VECT_CPU,1,2,igaus) = (-xjacm(DEF_VECT_CPU,1,2) * xjacm(DEF_VECT_CPU,3,3) + xjacm(DEF_VECT_CPU,3,2) * xjacm(DEF_VECT_CPU,1,3)) * denom(DEF_VECT_CPU)
          xjaci(DEF_VECT_CPU,1,3,igaus) = ( xjacm(DEF_VECT_CPU,1,2) * xjacm(DEF_VECT_CPU,2,3) - xjacm(DEF_VECT_CPU,2,2) * xjacm(DEF_VECT_CPU,1,3)) * denom(DEF_VECT_CPU)
          xjaci(DEF_VECT_CPU,2,3,igaus) = (-xjacm(DEF_VECT_CPU,1,1) * xjacm(DEF_VECT_CPU,2,3) + xjacm(DEF_VECT_CPU,2,1) * xjacm(DEF_VECT_CPU,1,3)) * denom(DEF_VECT_CPU)

          do j = 1, pnode
             gpcar(DEF_VECT_CPU,1,j,igaus) =   xjaci(DEF_VECT_CPU,1,1,igaus) * deriv(1,j,igaus) &
                  &                           + xjaci(DEF_VECT_CPU,2,1,igaus) * deriv(2,j,igaus) &
                  &                           + xjaci(DEF_VECT_CPU,3,1,igaus) * deriv(3,j,igaus)

             gpcar(DEF_VECT_CPU,2,j,igaus) =   xjaci(DEF_VECT_CPU,1,2,igaus) * deriv(1,j,igaus) &
                  &                           + xjaci(DEF_VECT_CPU,2,2,igaus) * deriv(2,j,igaus) &
                  &                           + xjaci(DEF_VECT_CPU,3,2,igaus) * deriv(3,j,igaus)

             gpcar(DEF_VECT_CPU,3,j,igaus) =   xjaci(DEF_VECT_CPU,1,3,igaus) * deriv(1,j,igaus) &
                  &                           + xjaci(DEF_VECT_CPU,2,3,igaus) * deriv(2,j,igaus) &
                  &                           + xjaci(DEF_VECT_CPU,3,3,igaus) * deriv(3,j,igaus)
          end do

       end do
    end if

  end subroutine elmgeo_cartesian_derivatives_jacobian_vectorized_cpu

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   Element characteristic length
  !> @details \verbatim
  !>          HLENG ... Element length with:
  !>                    HLENG(1)     = Max length
  !>                    HLENG(NDIME) = Min length
  !>          \endverbatim
  !>
  !-----------------------------------------------------------------------

  pure subroutine elmgeo_element_characteristic_length_vectorized(&
       VECTOR_DIM,ndime,pnode,deriv,elcod,hleng,hnatu_opt,tragl_opt)

    integer(ip), intent(in)           :: VECTOR_DIM
    integer(ip), intent(in)           :: ndime
    integer(ip), intent(in)           :: pnode
    real(rp),    intent(in)           :: deriv(ndime,pnode)
    real(rp),    intent(in)           :: elcod(VECTOR_DIM,ndime,pnode)
    real(rp),    intent(out)          :: hleng(VECTOR_DIM,ndime)
    real(rp),    intent(in), optional :: hnatu_opt
    real(rp),    intent(out),optional :: tragl_opt(VECTOR_DIM,ndime,ndime)
    integer(ip)                       :: ivect
    integer(ip)                       :: inode
    real(rp)                          :: enor0(VECTOR_DIM)
    real(rp)                          :: tragl(VECTOR_DIM,ndime,ndime)
    real(rp)                          :: xjacm(VECTOR_DIM,3,3)
    real(rp)                          :: gpdet(VECTOR_DIM)
    real(rp)                          :: denom(VECTOR_DIM)
    real(rp)                          :: t1(VECTOR_DIM)
    real(rp)                          :: t2(VECTOR_DIM)
    real(rp)                          :: t3(VECTOR_DIM)
    real(rp)                          :: h_tem,hnatu

    if( present(hnatu_opt) ) then
       hnatu = hnatu_opt
    else
       hnatu = 1.0_rp
    end if

    if( ndime == 1 ) then

#ifdef OPENACCHHH
       !$acc enter data create(xjacm, tragl, enor0)

       !$acc parallel loop gang vector default(present)
       do ivect = 1,VECTOR_DIM
#endif       
          xjacm(DEF_VECT,:,:) = 0.0_rp
          do inode = 1,pnode
             xjacm(DEF_VECT,1,1) = xjacm(DEF_VECT,1,1) + elcod(DEF_VECT,1,inode) * deriv(1,inode)
          end do
          tragl(DEF_VECT,1,1) = 1.0_rp / (sign(1.0_rp,xjacm(DEF_VECT,1,1))*max(abs(xjacm(DEF_VECT,1,1)),epsilgeo_div))
          enor0(DEF_VECT)     = tragl(DEF_VECT,1,1) * tragl(DEF_VECT,1,1)
          hleng(DEF_VECT,1)   = hnatu/sqrt(enor0(DEF_VECT))
#ifdef OPENACCHHH
       end do
      !$acc end parallel loop

!@      !$acc update self (hleng)

      !$acc exit data delete(xjacm, tragl, enor0)


#endif       

    else if( ndime == 2 ) then

#ifdef OPENACCHHH
       !$acc enter data create(xjacm, gpdet, denom, tragl, enor0) 

       !$acc parallel loop gang vector default(present)
       do ivect = 1,VECTOR_DIM
#endif

          xjacm(DEF_VECT,:,:) = 0.0_rp
          do inode = 1,pnode
             xjacm(DEF_VECT,1,1) = xjacm(DEF_VECT,1,1) + elcod(DEF_VECT,1,inode) * deriv(1,inode)
             xjacm(DEF_VECT,1,2) = xjacm(DEF_VECT,1,2) + elcod(DEF_VECT,1,inode) * deriv(2,inode)
             xjacm(DEF_VECT,2,1) = xjacm(DEF_VECT,2,1) + elcod(DEF_VECT,2,inode) * deriv(1,inode)
             xjacm(DEF_VECT,2,2) = xjacm(DEF_VECT,2,2) + elcod(DEF_VECT,2,inode) * deriv(2,inode)
          end do

          gpdet(DEF_VECT)     =  xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,2,2) - xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,1,2)
          denom(DEF_VECT)     =  1.0_rp / (sign(1.0_rp,gpdet(DEF_VECT))*max(abs(gpdet(DEF_VECT)),epsilgeo_div))
          !denom(DEF_VECT)     =  1.0_rp / gpdet(DEF_VECT)

          tragl(DEF_VECT,1,1) =  xjacm(DEF_VECT,2,2) * denom(DEF_VECT)
          tragl(DEF_VECT,2,2) =  xjacm(DEF_VECT,1,1) * denom(DEF_VECT)
          tragl(DEF_VECT,2,1) = -xjacm(DEF_VECT,2,1) * denom(DEF_VECT)
          tragl(DEF_VECT,1,2) = -xjacm(DEF_VECT,1,2) * denom (DEF_VECT)

          enor0(DEF_VECT)     =  tragl(DEF_VECT,1,1) * tragl(DEF_VECT,1,1) + tragl(DEF_VECT,1,2) * tragl(DEF_VECT,1,2)
          denom(DEF_VECT)    =   1.0_rp / max(sqrt(enor0(DEF_VECT)),epsilgeo_div)
          hleng(DEF_VECT,1)  =   hnatu * denom(DEF_VECT)

          enor0(DEF_VECT)     =  tragl(DEF_VECT,2,1) * tragl(DEF_VECT,2,1) + tragl(DEF_VECT,2,2) * tragl(DEF_VECT,2,2)
          denom(DEF_VECT)    =   1.0_rp / max(sqrt(enor0(DEF_VECT)),epsilgeo_div)
          hleng(DEF_VECT,2)  =   hnatu * denom(DEF_VECT)

#ifdef OPENACCHHH
       end do
      !$acc end parallel loop

!@      !$acc update self (hleng)

      !$acc exit data delete(xjacm, gpdet, denom, tragl, enor0)


#endif       

       do ivect = 1,VECTOR_DIM
          if( hleng(ivect,2) > hleng(ivect,1) ) then
             h_tem          = hleng(ivect,2)
             hleng(ivect,2) = hleng(ivect,1)
             hleng(ivect,1) = h_tem
          end if
       end do

    else if( ndime == 3 ) then
       !
       ! tragl = xjacm^-1
       ! xjacm = elcod * deriv^t
       !
#ifdef OPENACCHHH
       !$acc enter data create(xjacm, gpdet, denom, tragl, enor0, t1, &
       !$acc                   t2, t3) 

       !$acc parallel loop gang vector default(present)
       do ivect = 1,VECTOR_DIM
#endif
          xjacm(DEF_VECT,:,:) = 0.0_rp
          do inode = 1,pnode
             xjacm(DEF_VECT,1,1) = xjacm(DEF_VECT,1,1) + elcod(DEF_VECT,1,inode) * deriv(1,inode)
             xjacm(DEF_VECT,1,2) = xjacm(DEF_VECT,1,2) + elcod(DEF_VECT,1,inode) * deriv(2,inode)
             xjacm(DEF_VECT,1,3) = xjacm(DEF_VECT,1,3) + elcod(DEF_VECT,1,inode) * deriv(3,inode)
             xjacm(DEF_VECT,2,1) = xjacm(DEF_VECT,2,1) + elcod(DEF_VECT,2,inode) * deriv(1,inode)
             xjacm(DEF_VECT,2,2) = xjacm(DEF_VECT,2,2) + elcod(DEF_VECT,2,inode) * deriv(2,inode)
             xjacm(DEF_VECT,2,3) = xjacm(DEF_VECT,2,3) + elcod(DEF_VECT,2,inode) * deriv(3,inode)
             xjacm(DEF_VECT,3,1) = xjacm(DEF_VECT,3,1) + elcod(DEF_VECT,3,inode) * deriv(1,inode)
             xjacm(DEF_VECT,3,2) = xjacm(DEF_VECT,3,2) + elcod(DEF_VECT,3,inode) * deriv(2,inode)
             xjacm(DEF_VECT,3,3) = xjacm(DEF_VECT,3,3) + elcod(DEF_VECT,3,inode) * deriv(3,inode)
          end do

          t1(DEF_VECT)        =  xjacm(DEF_VECT,2,2) * xjacm(DEF_VECT,3,3) - xjacm(DEF_VECT,3,2) * xjacm(DEF_VECT,2,3)
          t2(DEF_VECT)        = -xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,3,3) + xjacm(DEF_VECT,3,1) * xjacm(DEF_VECT,2,3)
          t3(DEF_VECT)        =  xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,3,2) - xjacm(DEF_VECT,3,1) * xjacm(DEF_VECT,2,2)

          gpdet(DEF_VECT)     =  xjacm(DEF_VECT,1,1) * t1(DEF_VECT) + xjacm(DEF_VECT,1,2) * t2 (DEF_VECT)+ xjacm(DEF_VECT,1,3) * t3(DEF_VECT)
          denom(DEF_VECT)     =  1.0_rp / (sign(1.0_rp,gpdet(DEF_VECT))*max(abs(gpdet(DEF_VECT)),epsilgeo_div))
          !denom(DEF_VECT)     =  1.0_rp / gpdet(DEF_VECT)

          tragl(DEF_VECT,1,1) =  t1(DEF_VECT) * denom(DEF_VECT)
          tragl(DEF_VECT,2,1) =  t2(DEF_VECT) * denom(DEF_VECT)
          tragl(DEF_VECT,3,1) =  t3(DEF_VECT) * denom(DEF_VECT)
          tragl(DEF_VECT,2,2) =  ( xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,3,3) - xjacm(DEF_VECT,3,1) * xjacm(DEF_VECT,1,3)) * denom(DEF_VECT)
          tragl(DEF_VECT,3,2) =  (-xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,3,2) + xjacm(DEF_VECT,1,2) * xjacm(DEF_VECT,3,1)) * denom(DEF_VECT)
          tragl(DEF_VECT,3,3) =  ( xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,2,2) - xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,1,2)) * denom(DEF_VECT)
          tragl(DEF_VECT,1,2) =  (-xjacm(DEF_VECT,1,2) * xjacm(DEF_VECT,3,3) + xjacm(DEF_VECT,3,2) * xjacm(DEF_VECT,1,3)) * denom(DEF_VECT)
          tragl(DEF_VECT,1,3) =  ( xjacm(DEF_VECT,1,2) * xjacm(DEF_VECT,2,3) - xjacm(DEF_VECT,2,2) * xjacm(DEF_VECT,1,3)) * denom(DEF_VECT)
          tragl(DEF_VECT,2,3) =  (-xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,2,3) + xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,1,3)) * denom(DEF_VECT)
          !
          ! Element length HLENG
          !
          enor0(DEF_VECT)    =   tragl(DEF_VECT,1,1) * tragl(DEF_VECT,1,1) &
               &                    + tragl(DEF_VECT,1,2) * tragl(DEF_VECT,1,2) &
               &                    + tragl(DEF_VECT,1,3) * tragl(DEF_VECT,1,3)

          denom(DEF_VECT)    =   1.0_rp / max(sqrt(enor0(DEF_VECT)),epsilgeo_div)
          hleng(DEF_VECT,1)  =   hnatu * denom(DEF_VECT)

          enor0(DEF_VECT)    =   tragl(DEF_VECT,2,1) * tragl(DEF_VECT,2,1) &
               &                    + tragl(DEF_VECT,2,2) * tragl(DEF_VECT,2,2) &
               &                    + tragl(DEF_VECT,2,3) * tragl(DEF_VECT,2,3)
          denom(DEF_VECT)    =   1.0_rp / max(sqrt(enor0(DEF_VECT)),epsilgeo_div)
          hleng(DEF_VECT,2)  =   hnatu * denom(DEF_VECT)

          enor0(DEF_VECT)    =   tragl(DEF_VECT,3,1) * tragl(DEF_VECT,3,1) &
               &                    + tragl(DEF_VECT,3,2) * tragl(DEF_VECT,3,2) &
               &                    + tragl(DEF_VECT,3,3) * tragl(DEF_VECT,3,3)
          denom(DEF_VECT)    =   1.0_rp / max(sqrt(enor0(DEF_VECT)),epsilgeo_div)
          hleng(DEF_VECT,3)  =   hnatu * denom(DEF_VECT)

#ifdef OPENACCHHH
       end do
      !$acc end parallel loop

!@      !$acc update self (hleng)

      !$acc exit data delete(xjacm, gpdet, denom, tragl, enor0,  &
      !$acc                  t1, t2, t3)


#endif       
       !
       ! Sort hleng: hleng(1)=max; hleng(ndime)=min
       !
       !$acc parallel loop gang vector default(present)
       do ivect = 1,VECTOR_DIM
          if( hleng(ivect,2) > hleng(ivect,1) ) then
             h_tem          = hleng(ivect,2)
             hleng(ivect,2) = hleng(ivect,1)
             hleng(ivect,1) = h_tem
          end if
          if( hleng(ivect,3) > hleng(ivect,1) ) then
             h_tem          = hleng(ivect,3)
             hleng(ivect,3) = hleng(ivect,1)
             hleng(ivect,1) = h_tem
          end if
          if( hleng(ivect,3) > hleng(ivect,2) ) then
             h_tem          = hleng(ivect,3)
             hleng(ivect,3) = hleng(ivect,2)
             hleng(ivect,2) = h_tem
          end if
       end do
      !$acc end parallel loop


    end if

    if( present(tragl_opt) ) then
       tragl_opt = tragl
    end if

  end subroutine elmgeo_element_characteristic_length_vectorized
  
end module mod_elmgeo_vector
!> @}
