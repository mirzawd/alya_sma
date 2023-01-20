!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Domain 
!> @{
!> @file    cshder.f90
!> @author  Guillaume Houzeaux
!> @date    10/10/1972
!> @brief   Shape functions, derivatives and Hessians
!> @details Shape functions, derivatives and Hessians:
!>          \verbatim
!>          - For each element type, using user integration rule:
!>            WEIGP(ngaus)
!>            SHAPE(nnode,ngaus)
!>            DERIV(ndime,nnode,ngaus)
!>            HESLO(ntens,nnode,ngaus)
!>          - For each element type, the bubble:
!>            SHAPE_BUB(ngaus)
!>            DERIV_BUB(ndime,ngaus)
!>            HESLO_BUB(ntens,ngaus)
!>          - For each element type, using a closed integration rule:
!>            WEIGC(nnode)
!>            SHAPC(nnode,nnode)
!>            DERIC(ndime,nnode,nnode)
!>            HESLC(ntens,nnode,nnode)
!>          - Center of gravity:
!>            SHACG(nnode)
!>            DERCG(ndime,nnode)
!>            WEICG
!>          - Element Gauss points to nodes: 
!>            SHAGA(ngaus,nnode) 
!>          \endverbatim
!> @} 
!------------------------------------------------------------------------

subroutine cshder(itask)
  use def_kintyp,        only : ip,rp
  use def_domain,        only : elmar
  use def_domain,        only : nnode
  use def_domain,        only : ngaus
  use def_domain,        only : ldime
  use def_domain,        only : lrule
  use def_domain,        only : linte
  use def_domain,        only : lquad
  use def_domain,        only : ltopo
  use def_domain,        only : lexis
  use def_domain,        only : mnode
  use def_domain,        only : nelty
  use def_domain,        only : ntens
  use def_domain,        only : ndime
  use def_domain,        only : memor_dom
  use mod_memory,        only : memory_alloca
  use mod_memory,        only : memory_deallo
  use def_elmgeo,        only : element_type
  use mod_elmgeo,        only : elmgeo_shapf_deriv_heslo_bubble
  use mod_elmgeo,        only : elmgeo_shapf_deriv_heslo
  use def_isoparametric, only : LAGRANGE_INTERPOLATION
  use def_isoparametric, only : CHEBYSHEV_INTERPOLATION
  use mod_rulepw,        only : rulepw
  use def_quadrature,    only : CLOSED_RULE
  use def_quadrature,    only : CHEBYSHEV_RULE
  use def_quadrature,    only : GAUSS_LEGENDRE_RULE   
  use def_elmtyp
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: igaus,ielty,inode,mdime
  integer(ip)             :: pnode,pgaus,pdime,prule,ptopo,pinte
  integer(ip)             :: pquad
  real(rp)                :: poscg(ndime) 
  real(rp)                :: posgw(ndime,mnode) 
  real(rp)                :: weigw(mnode)
  real(rp)                :: posgc(ndime,mnode)
  integer(ip)             :: ierro=0

  if( itask == 1 ) then

     !----------------------------------------------------------------------
     !
     ! Allocate memory for structure ELMAR
     !
     !----------------------------------------------------------------------

     allocate( elmar(nelty) )     
     do ielty = 1,nelty
        nullify( elmar(ielty) % shape     )
        nullify( elmar(ielty) % deriv     )
        nullify( elmar(ielty) % heslo     )
        nullify( elmar(ielty) % posgp     )    
        nullify( elmar(ielty) % weigp     )
        nullify( elmar(ielty) % shape_bub )
        nullify( elmar(ielty) % deriv_bub )
        nullify( elmar(ielty) % heslo_bub )
        nullify( elmar(ielty) % shacg     )
        nullify( elmar(ielty) % dercg     )
        nullify( elmar(ielty) % hescg     )   
        nullify( elmar(ielty) % shaga     )
        nullify( elmar(ielty) % shapc     )
        nullify( elmar(ielty) % deric     )
        nullify( elmar(ielty) % heslc     )
        nullify( elmar(ielty) % weigc     )    
        nullify( elmar(ielty) % shaib     )
        nullify( elmar(ielty) % derib     )
        nullify( elmar(ielty) % weiib     )
     end do

  else if( itask == 2 ) then 

     !----------------------------------------------------------------------
     !
     ! Deallocate memory for structure ELMAR
     !
     !----------------------------------------------------------------------

     do ielty = 1,nelty
        call memory_deallo(memor_dom,'SHAPE'    ,'cshder',elmar(ielty) % shape)
        call memory_deallo(memor_dom,'DERIV'    ,'cshder',elmar(ielty) % deriv)
        call memory_deallo(memor_dom,'HESLO'    ,'cshder',elmar(ielty) % heslo)
        call memory_deallo(memor_dom,'POSGP'    ,'cshder',elmar(ielty) % posgp)
        call memory_deallo(memor_dom,'WEIGP'    ,'cshder',elmar(ielty) % weigp)

        call memory_deallo(memor_dom,'SHAPE_BUB','cshder',elmar(ielty) % shape_bub)
        call memory_deallo(memor_dom,'DERIV_BUB','cshder',elmar(ielty) % deriv_bub)
        call memory_deallo(memor_dom,'HESLO_BUB','cshder',elmar(ielty) % heslo_bub)

        call memory_deallo(memor_dom,'SHACG'    ,'cshder',elmar(ielty) % shacg)
        call memory_deallo(memor_dom,'DERCG'    ,'cshder',elmar(ielty) % dercg)
        call memory_deallo(memor_dom,'HESCG'    ,'cshder',elmar(ielty) % hescg)

        call memory_deallo(memor_dom,'SHAGA'    ,'cshder',elmar(ielty) % shaga)

        call memory_deallo(memor_dom,'SHAPC'    ,'cshder',elmar(ielty) % shapc)
        call memory_deallo(memor_dom,'DERIC'    ,'cshder',elmar(ielty) % deric)
        call memory_deallo(memor_dom,'HESLC'    ,'cshder',elmar(ielty) % heslc)
        call memory_deallo(memor_dom,'WEIGC'    ,'cshder',elmar(ielty) % weigc)
     end do

     deallocate(elmar)

  else if( itask == 3 ) then

     !----------------------------------------------------------------------
     !
     ! Allocate: 3. Deallocate: 4
     !
     !----------------------------------------------------------------------

     ierro = 0
     poscg = 0.0_rp
     posgw = 0.0_rp
     weigw = 0.0_rp
     posgc = 0.0_rp
     !
     ! Element shape function and derivatives: SHAPE,DERIV,HESLO,WEIGP 
     !
     do ielty = 1,nelty
        if( lexis(ielty) == 1 ) then 
           
           pnode = nnode(ielty) 
           pgaus = ngaus(ielty)
           pdime = ldime(ielty)
           prule = lrule(ielty)
           pinte = linte(ielty)
           ptopo = ltopo(ielty)
           pquad = lquad(ielty)
           mdime = max(1_ip,pdime)
           elmar(ielty) % pgaus = pgaus
           
           call memory_alloca(memor_dom,'SHAPE'    ,'cshder',elmar(ielty) % shape,      pnode,pgaus)
           call memory_alloca(memor_dom,'DERIV'    ,'cshder',elmar(ielty) % deriv,mdime,pnode,pgaus)
           call memory_alloca(memor_dom,'HESLO'    ,'cshder',elmar(ielty) % heslo,ntens,pnode,pgaus)
           call memory_alloca(memor_dom,'POSGP'    ,'cshder',elmar(ielty) % posgp,ndime,      pgaus)
           call memory_alloca(memor_dom,'WEIGP'    ,'cshder',elmar(ielty) % weigp,            pgaus)

           call memory_alloca(memor_dom,'SHAPE_BUB','cshder',elmar(ielty) % shape_bub,      1_ip,pgaus)
           call memory_alloca(memor_dom,'DERIV_BUB','cshder',elmar(ielty) % deriv_bub,mdime,1_ip,pgaus)
           call memory_alloca(memor_dom,'HESLO_BUB','cshder',elmar(ielty) % heslo_bub,ntens,1_ip,pgaus)
        
           call rulepw(pdime,pgaus,ptopo,pquad,elmar(ielty) % posgp,elmar(ielty) % weigp,ierro)

           if( ierro /= 0 ) then
              call runend('CSHDER: WRONG INTEGRATION RULE FOR ELEMENT '//element_type(ielty) % name)
           else
              call elmgeo_shapf_deriv_heslo(&
                   pdime,pnode,pgaus,elmar(ielty) % posgp,elmar(ielty) % shape,&
                   elmar(ielty) % deriv,elmar(ielty) % heslo,ierro,pinte)
              call elmgeo_shapf_deriv_heslo_bubble(&
                   mdime,1_ip,pgaus,elmar(ielty) % posgp,elmar(ielty) % shape_bub,&
                   elmar(ielty) % deriv_bub,elmar(ielty) % heslo_bub)
           end if
        end if
     end do
     !
     ! Element Center of gravity shape function and derivatives SHACG,DERCG,HESCG,WEICG 
     !
     do ielty = 1,nelty
        if( lexis(ielty) == 1 ) then

           pnode = nnode(ielty)
           pgaus = ngaus(ielty)
           pdime = ldime(ielty)
           mdime = max(1_ip,pdime)
           pquad = GAUSS_LEGENDRE_RULE   
           ptopo = ltopo(ielty)
           pinte = LAGRANGE_INTERPOLATION
           
           call memory_alloca(memor_dom,'SHACG','cshder',elmar(ielty) % shacg,pnode)
           call memory_alloca(memor_dom,'DERCG','cshder',elmar(ielty) % dercg,mdime,pnode)
           call memory_alloca(memor_dom,'HESCG','cshder',elmar(ielty) % hescg,ntens,pnode)

           call rulepw(pdime,1_ip,ptopo,pquad,poscg,elmar(ielty)%weicg,ierro)
           
           call elmgeo_shapf_deriv_heslo(&
                pdime,pnode,poscg,elmar(ielty)%shacg,elmar(ielty)%dercg,elmar(ielty)%hescg,&
                ierro,pinte)
           
        end if
     end do
     !
     ! Computes the interpolation functions associated to the integration
     ! points in order to extrapolate values from these integration points
     ! to the nodes
     !
     do ielty = 1,nelty
        if( lexis(ielty) == 1 ) then

           pnode = nnode(ielty)
           pgaus = ngaus(ielty)
           pdime = ldime(ielty)
           pquad = lquad(ielty)
           mdime = max(1_ip,pdime)
           pquad = CLOSED_RULE
           ptopo = ltopo(ielty)
           pinte = linte(ielty)

           call memory_alloca(memor_dom,'SHAGA','cshder',elmar(ielty) % shaga,pgaus,pnode)
           
           if( lquad(ielty) == 1 .or. pinte == CHEBYSHEV_INTERPOLATION ) then

              do inode = 1,pnode
                 do igaus = 1,pgaus
                    elmar(ielty) % shaga(igaus,inode) = 0.0_rp
                 end do
                 elmar(ielty) % shaga(inode,inode) = 1.0_rp
              end do
              
           else
              
              call rulepw(pdime,pnode,ptopo,pquad,posgw,weigw,ierro)
              call shafga(&
                   posgw,pdime,ptopo,pgaus,pnode,&
                   elmar(ielty)%shaga,ierro)
           end if

        end if
     end do
     !
     ! Element shape function and derivatives SHAPC,DERIC,HESLC,WEIGC 
     ! for a close rule
     !
     do ielty = 1,nelty
        if( lexis(ielty) == 1 ) then

           pnode = nnode(ielty)
           pdime = ldime(ielty)
           mdime = max(1_ip,pdime)
           ptopo = ltopo(ielty)
           pinte = linte(ielty)
           if( pinte == LAGRANGE_INTERPOLATION ) then
              pquad = CLOSED_RULE              
           else
              pquad = CHEBYSHEV_RULE              
           end if

           call memory_alloca(memor_dom,'SHAPC','cshder',elmar(ielty) % shapc,pnode,pnode)
           call memory_alloca(memor_dom,'DERIC','cshder',elmar(ielty) % deric,mdime,pnode,pnode)
           call memory_alloca(memor_dom,'HESLC','cshder',elmar(ielty) % heslc,ntens,pnode,pnode)
           call memory_alloca(memor_dom,'WEIGC','cshder',elmar(ielty) % weigc,pnode)

           call rulepw(pdime,pnode,ptopo,pquad,posgc,elmar(ielty)%weigc,ierro)
           
           call elmgeo_shapf_deriv_heslo(&
                pdime,pnode,pnode,posgc,elmar(ielty) % shapc,&
                elmar(ielty) % deric,elmar(ielty) % heslc,ierro,pinte)
           
        end if
     end do

  end if

end subroutine cshder
