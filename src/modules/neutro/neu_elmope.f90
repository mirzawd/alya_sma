!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Neutro
!> @{
!> @file    neu_elmope.f90
!> @author  Guillaume Houzeaux
!> @brief   Element assembly
!> @details We assemble the elementary matrix of each element into the general matrix of the system
!> @} 
!-----------------------------------------------------------------------

subroutine neu_elmope(order)

  use def_parame
  use def_elmtyp
  use def_master
  use def_domain
  use def_kermod
  use mod_ker_proper 
  use def_neutro
  use mod_ADR,    only : ADR_element_assembly
  use mod_ADR,    only : mreac_adr
  use mod_ADR,    only : ELEMENT_ASSEMBLY
  use mod_matrix, only : matrix_assemble_element_matrix_to_CSR
  use mod_matrix, only : matrix_assemble_element_RHS

  implicit none

  integer(ip), intent(in) :: order                      ! =2: compute SGS only

  real(rp)    :: elmat(mnode,mnode)                     !> Ae, elementary matrix (whose dimensions are 
                                                        !>     the maximum number of nodes per element)
  real(rp)    :: elrhs(mnode)                           !> be, elementary rhs (with the same length as one dimension of the matrix)
  integer(ip) :: ielem,pelty,pmate,pnode                !> Indices and dimensions: element identifier, element types, number of 
                                                        !>                         materials and nodes in the element
  integer(ip) :: pgaus,plapl,porde,ptopo

  real(rp)    :: elrad(num_energies_neu,num_directions_neu,mnode) !> Gather: here we save the elementary radiation for 
                                                                  !>         each energy, direction and node 
  real(rp)    :: elunk(mnode,ADR_neu % ntime)           ! Gather 
  real(rp)    :: elcod(ndime,mnode)               !> Element code

  real(rp)    :: tragl(ndime,ndime)                     ! Stabilization
  real(rp)    :: chale(2),hleng(ndime),dummr2(1)!,dummr

  real(rp)    :: gpvol(mgaus)                           !> |J|*w volume
  real(rp)    :: gprea(mgaus,mreac_adr)                 !> This actually accounts for the lost neutrons due to total interaction
  real(rp)    :: gpadv(ndime,mgaus)                     !> a: advection term, that is, the velocity field in 
                                                        !>    each dimension and Gauss point
  real(rp)    :: gpdif(mgaus),gpgrd(ndime,mgaus)        !> k+kt, grad(k+kt)
  real(rp)    :: gprhs(mgaus)                           !> f (all terms), right-hand side of the equation at Gauss point mgaus
  real(rp)    :: gpden(mgaus)                           !> rho and then rho*cp, density parameter (which I think 
                                                        !>     we don't need for the moment)
  real(rp)    :: gpunk(mgaus,ADR_neu % ntime)           !> u (unkown=neutron flux, at a certain time)
  real(rp)    :: gpsou(mgaus)                           !> Q, external source term (depending on the maximum number of Gauss points)
  real(rp)    :: gpsha(mnode,mgaus)                     !> N, that is, the shape functions we use to approximate the 
                                                        !>    solution (depends on the maximum number of nodes and Gauss points)
  real(rp)    :: gpder(ndime,mnode,mgaus)               !> dNk/dsj, derivatives of the shape function with respect to the local 
                                                        !>          coordinates of the element
  real(rp)    :: gpcar(ndime,mnode,mgaus)               !> dNk/dxj, derivative of shape function k with respect to 
                                                        !>          coordinate x_j for each node and Gauss point
  real(rp)    :: gphes(ntens,mnode,mgaus)               !> dNk/dxidxj, Hessian, derivative of shape function k with respect to 
                                                        !>             coordinates x_i and x_j
  real(rp)    :: gpabs(mgaus)                           !> Absorption coefficient (in the future this must be 
                                                        !>                         converted into a matrix)
  real(rp)    :: gpsca(mgaus)                           !> Scattering coefficient (in the future this must be 
                                                        !>                         converted into a matrix)

!   integer(ip) :: inode,jnode,ipoin 
!   real(rp)    :: adiag,xvalu

  external :: neu_elmgat, elmlen, elmchl, elmca2, neu_coefficients

  gpden = 1.0_rp !> We set the density=1 Does this have any effect on our problem? 
  gpdif = 0.0_rp !> Diffusion coefficient set to zero since there is no diffusion term in our equation
  gpgrd = 0.0_rp !> Gradient of the diffusion coefficient (also zero since no diffusion)
  gphes = 0.0_rp !> We initialize the Hessian
  gprea = 0.0_rp !> The total XS for the moment is initialized at zero
  gpsou = 0.0_rp !> The source term is also initialized at zero
  !
  ! Loop over elements   ! aca hay que incluir 
  !                      /1   TotalXS,    **Listo
  !                      /2   Scatt matrix, **Listo
  !                      /3   armonicos para la scatt matrix, **Listo
  !                      4   fision matrix   **Falta a futuro
  !  
  elements: do ielem = 1,nelem
     !
     ! Element dimensions
     !
     pelty = ltype(ielem)

     if( pelty > 0 ) then !> We check that the element is at least 1D
        pnode = lnnod(ielem) !> number of nodes that element ielem has
        pgaus = ngaus(pelty) !> number of Gauss points there are
        plapl = 0 !> What are these variables??
        porde = lorde(pelty)
        ptopo = ltopo(pelty)
        pmate = 1 !> We set that the current material is of type 1 
        if( nmate > 1 ) pmate = lmate(ielem) !> If there is more than one type of material,
                                             !> we identify the type of material of the current element
        !
        ! Gather operations 
        !
        call neu_elmgat(&
             pnode,lnods(1:pnode,ielem),elcod,elrad,elunk)
        !
        ! hleng and tragl at center of gravity
        !
        call elmlen(&
             ndime,pnode,elmar(pelty) % dercg,tragl,elcod,hnatu(pelty),hleng)
        !
        ! Compute the characteristic length CHALE
        !
        call elmchl(&
             tragl,hleng,elcod,dummr2,dummr2,chale,pelty,pnode,&
             porde,hnatu(pelty),0_ip,0_ip) 
        ! 
        ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, GPVOL
        !
        call elmca2(&
             pnode,pgaus,plapl,elmar(pelty) % weigp,elmar(pelty) % shape,&
             elmar(pelty) % deriv,elmar(pelty) % heslo,elcod,gpvol,gpsha,&
             gpder,gpcar,gphes,ielem)
        !
        ! Properties: GPABS (absorption) and GPSCA (scattering) read for each material
        !
       
       ! call ker_proper_scalar('ABSOR','PGAUS',0_ip,ielem,gpabs,pnode,pgaus,gpsha,gpcar) 
       ! call ker_proper_scalar('SCATT','PGAUS',0_ip,ielem,gpsca,pnode,pgaus,gpsha,gpcar)

        gpabs=absor_neu(pmate,current_energy_neu)  ! OJO absorbcion vector en funcion de grupo
       
       !
       ! Equation coefficients: GPADV, GPREA, GPRHS, GPUNK
       !
        !> We determine the total XS (gprea), the RHS (gprhs), the current direction vector (gpadv)
        call neu_coefficients(pnode,pgaus,pmate,gpsha,gpabs,gpsca,elrad,elcod,elunk,gpadv,gprea,gprhs,gpunk) 


        ! "DIRICHLET" 
        ! inside direction in vacuum, 0.0 
     !        do inode = 1,pnode 
     !              ipoin = lnods(inode,ielem) !lnods(inode) 
     !             if( kfl_fixno_neu(current_energy_neu,current_direction_neu)% l(1,ipoin) == 1 )  then
     !                 adiag = elmat(inode,inode) 
     !                 xvalu = 0.0 
     !                do jnode = 1,pnode 
     !                     elmat(inode,jnode) = 0.0_rp 
     !                     elrhs(jnode) = elrhs(jnode) - elmat(jnode,inode) * xvalu 
     !                     elmat(jnode,inode) = 0.0_rp
     !                end do 
     !                if( abs(adiag) > 0.0_rp ) then 
     !                        elmat(inode,inode) = adiag 
     !                     elrhs(inode) = adiag * xvalu 
     !                else 
     !                 elmat(inode,inode) = 1.0_rp 
     !                   elrhs(inode) = xvalu 
     !               end if 
     !           endif
     !    end do


       
        if( order == ELEMENT_ASSEMBLY ) then !> (This always happens when the function is called from neu_solite.f90)
           ! 
           ! Compute and assemble element matrix and RHS 
           ! 
           call ADR_element_assembly(&
                ielem,pnode,pgaus,elcod,gpsha,gpcar,elmar(pelty) % deriv,gphes,gpvol,chale,&
                elmar(pelty) % shape_bub,elmar(pelty) % deriv_bub,ADR_neu,&
                cutim,gpden,gpadv,gpdif,gpgrd,gprea,gprhs,&
                gpunk,elunk,elmat,elrhs) !> Elemental assembly of the ADR equation (see kernel/mathru/mod_ADR.f90)
           call matrix_assemble_element_RHS(&
                1_ip,1_ip,pnode,lnods(:,ielem),elrhs,rhsid) !> We assemble the elementary rhs into the system RHS
           call matrix_assemble_element_matrix_to_CSR(&
                kfl_element_to_csr,1_ip,pnode,pnode,&
                ielem,lnods(:,ielem),elmat,r_dom,c_dom,amatr,lezdo) !> We assemble the elementary matrix into the system 
                                                                    ! matrix (see kernel/domain/mod_matrix.f90) 

        else
           !
           ! Others
           !
           call runend('NEU_ELMOPE: DONT KNOW WHAT TO DO!')

        end if

     end if

  end do elements

end subroutine neu_elmope
