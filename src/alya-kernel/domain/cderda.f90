!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine cderda()
  !-----------------------------------------------------------------------
  !****f* Domain/cderda
  ! NAME
  !    cderda
  ! DESCRIPTION
  !    This routine defines the derivated parameters of the element
  ! OUTPUT
  !    NDIMB ... Boundary dimension
  !    NTENS ... # of independent Hessian matrix values
  !    NINER ... # of independent tensor of inertia values
  !    LRULE ... 1 Quad/Hexa:  open   
  !          ... 2 Quad/Hexa:  closed 
  !          ... 3 Tria/Tetra: open
  !          ... 4 Tria/Tetra: closed
  !          ... 5   - /Penta: open
  !          ... 6   - /Penta: closed
  !          ... 7   - /Pyram: open
  !          ... 8   - /Pyram: closed
  !    HNATU ... 2 Quad/Hexa
  !          ... 1 Others
  !    MNODE ... Maximum # of element nodes in the mesh
  !    MGAUS ... Maximum # of element Gauss points in the mesh
  !    MNODB ... Maximum # of boundary nodes in the mesh
  !    MGAUB ... Maximum # of boundary Gaus points in the mesh
  !    MLAPL ... 1 if at least one element needs Laplacian
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
  use def_elmtyp
  use def_master
  use def_domain
  use def_elmtyp,        only : element_end
  use mod_memory,        only : memory_alloca
  use def_elmgeo,        only : element_type
  use def_quadrature,    only : GAUSS_LEGENDRE_RULE   
  use def_quadrature,    only : TRAPEZOIDAL_RULE      
  use def_quadrature,    only : CLOSED_RULE           
  use def_quadrature,    only : CHEBYSHEV_RULE           
  implicit none
  integer(ip) :: ielty,pgaus,pquad,ndife
  !
  ! Where element and boundary element types start and stop
  ! We put IBSTA_DOM=1 to treat BARD3D elements, which boundary
  ! is the 1D POINT element
  !
  iesta_dom = element_num_ini(ndime)
  iesto_dom = element_num_end(ndime)
  ibsta_dom = element_num_ini(ndime-1)
  ibsto_dom = element_num_end(ndime-1)
  if( lexis(BAR3D) /= 0 ) ibsta_dom = 1
  !
  ! Gobal dimensions
  !
  ndimb = ndime-1
  select case ( ndime )
  case ( 1 )   ; ntens = 1
  case default ; ntens = 3 * ndime - 3
  end select
  !
  ! Compute LRULE and natural length HNATU
  ! 
  do ielty = 1,nelty

     select case ( lquad(ielty) ) 
     case ( GAUSS_LEGENDRE_RULE , CLOSED_RULE ) ; lrule(ielty) =  2*ltopo(ielty)+lquad(ielty)+1
     case ( TRAPEZOIDAL_RULE                  ) ; lrule(ielty) = -3
     case ( CHEBYSHEV_RULE                    ) ; lrule(ielty) = -4
     end select

     hnatu(ielty) = element_type(ielty) % natural_length 
     
  end do
  !
  ! If Volume Gauss points have not been assigned, put default option
  !
  do ielty = 1,nelty     
     if( lexis(ielty) /= 0 .and. ngaus(ielty) <= 0 ) then
        if( lquad(ielty) == 0 ) then
           !
           ! Open rule
           !
           select case ( ielty ) 
           case ( POINT ) ; ngaus(ielty) =  1
           case ( BAR04 ) ; ngaus(ielty) = 11
           case ( QUA08 ) ; ngaus(ielty) =  9
           case ( QUA09 ) ; ngaus(ielty) =  9
           case ( TRI10 ) ; ngaus(ielty) = 13
           case ( TET20 ) ; ngaus(ielty) = 29
           case ( HEX27 ) ; ngaus(ielty) = 27
           case ( PEN18 ) ; ngaus(ielty) = 21
           case ( TET10 ) ; ngaus(ielty) = 14
           case default   ; ngaus(ielty) = nnode(ielty) 
           end select
           
        else
           !
           ! Close rule
           !
           ngaus(ielty) = nnode(ielty)
        end if
     end if
  end do
  !
  ! Treat elements of dimension ndime-1: NGAUS, LQUAD and LEXIS
  !
  do ielty = iesta_dom,iesto_dom 
     if( lexis(ielty) == 1 ) then
        pgaus = ngaus(ielty)
        pquad = lquad(ielty)
        call bouele(nelty,pgaus,pquad,ielty,ngaus,lquad,lexis)
     end if
  end do
  do ielty = 1,nelty
     if( lexis(ielty) == 1 .and. ngaus(ielty) == 0 ) then
        call runend('CDERDA: INTEGRATION RULE NOT DEFINED FOR ELEMENT '//element_type(ielty) % name)
     end if
  end do
  !
  ! Maximum values
  ! MNODE can be optionally given in dimension field
  !
#ifndef PNODE_VALUE
  mnode = max(mnode,-1_ip)
#endif
  mgaus = -1
  mnodb = -1
  mgaub = -1
  mlapl = -1
  do ielty = iesta_dom,iesto_dom
     if( lexis(ielty) == 1 ) then
#ifndef PNODE_VALUE
        mnode = max(mnode,nnode(ielty))
#endif
        mgaus = max(mgaus,ngaus(ielty))
        mlapl = max(mlapl,llapl(ielty))
     end if
  end do
  do ielty = ibsta_dom,ibsto_dom
     if( lexis(ielty) == 1 ) then
        mnodb = max(mnodb,nnode(ielty))
        mgaub = max(mgaub,ngaus(ielty))
     end if
  end do
  if( lexis(SHELL) == 1 ) mnodb = max(2_ip,mnodb)
  if( lexis(BAR3D) == 1 ) mnodb = max(2_ip,mnodb)
  mnoga = max(mnode,mgaus)
  !
  ! Faces
  !
  mface = maxval(element_type(:) % number_faces)
  !
  ! Check errors
  !  
  if( mnode == -1 ) call runend('DOMAIN: NO ELEMENT TYPE HAS BEEN DECLARED')
  if( mgaus == -1 ) call runend('DOMAIN: NO NUMERICAL INTEGRATION HAS BEEN DEFINED')
  
end subroutine cderda
