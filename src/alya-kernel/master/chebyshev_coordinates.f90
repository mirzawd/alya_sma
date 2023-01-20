!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    chebyshev_coordinates.f90
!> @author  houzeaux and gasparino
!> @date    2022-07-26
!> @brief   Move mesh ndoes
!> @details Place interior nodes of fintie element on their
!>          Chebyshev coordinates
!> @} 
!-----------------------------------------------------------------------

subroutine chebyshev_coordinates()

  use def_kintyp,        only : ip,rp
  use def_isoparametric, only : CHEBYSHEV_INTERPOLATION
  use def_chebyshev,     only : lagrange_roots
  use def_chebyshev,     only : tripleTensorProduct
  use def_chebyshev,     only : var_interpolate
  use def_elmtyp,        only : element_num_ini
  use def_elmtyp,        only : element_num_end 
  use def_elmgeo,        only : element_type
  use def_domain,        only : ndime
  use def_domain,        only : nelem
  use def_domain,        only : mnode
  use def_domain,        only : coord
  use def_domain,        only : linte
  use def_domain,        only : elmar
  use def_domain,        only : coord
  use def_domain,        only : lnods
  use def_domain,        only : ltype
  use def_domain,        only : lexis
  use def_domain,        only : memor_dom
  use mod_memory,        only : memory_copy
  use mod_memory,        only : memory_deallo
  use mod_messages,      only : messages_live

  implicit none

  integer(ip)               :: pnode,idime,igaus,porde
  integer(ip)               :: pgaus,ielem,pelty
  integer(ip), parameter    :: morde=11
  real(rp),    allocatable  :: N(:,:)
  real(rp)                  :: elcod(ndime,mnode)
  real(rp),    pointer      :: coord_cpy(:,:)

  if( any(linte==CHEBYSHEV_INTERPOLATION) ) then

     call messages_live('TRANSFORM COORDINATES FROM LAGRANGE TO CHEBYSHEV')
     nullify(coord_cpy)
     call memory_copy(memor_dom,'COORD_CPY','chebyshev_coordinates',coord,coord_cpy,'DO_NOT_DEALLOCATE')

     do pelty = element_num_ini(ndime),element_num_end(ndime)
        if( linte(pelty) == CHEBYSHEV_INTERPOLATION .and. lexis(pelty) == 1 ) then
           porde = element_type(pelty) % order
           pnode = element_type(pelty) % number_nodes
           pgaus = pnode
           allocate(N(pnode,pgaus))
           do igaus = 1,pgaus
              call chebyshev_lagrangian_shape(pnode,porde,N(:,igaus),elmar(pelty) % posgp(:,igaus))
           end do
           do ielem = 1,nelem    
              if( ltype(ielem) == pelty ) then
                 elcod(1:ndime,1:pnode) = coord_cpy(1:ndime,lnods(1:pnode,ielem))
                 do igaus = 1,pnode
                    do idime = 1,ndime
                       call var_interpolate(pnode,elcod(idime,1:pnode),N(:,igaus),coord(idime,lnods(igaus,ielem)))
                    end do
                 end do
              end if

           end do
           deallocate(N)
        end if
     end do
     call memory_deallo(memor_dom,'COORD_CPY','chebyshev_coordinates',coord_cpy)

  end if

contains

  subroutine chebyshev_lagrangian_shape(pnode,porde,N,x)
    implicit none

    integer(ip),           intent(in)  :: pnode
    integer(ip),           intent(in)  :: porde
    real(rp),              intent(in)  :: x(3)
    real(rp),    optional, intent(out) :: N(pnode)
    real(rp)                           :: dN(3,pnode)
    real(rp)                           :: xi_grid(porde+1)

    call lagrange_roots     (porde,xi_grid)
    call tripleTensorProduct(pnode,porde,xi_grid,x(1),x(2),x(3),N,dN)

  end subroutine chebyshev_lagrangian_shape

end subroutine chebyshev_coordinates
