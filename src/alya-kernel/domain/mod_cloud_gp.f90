!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
! Set of subroutines to manage a cloud of Gaussian points
!------------------------------------------------------------------------

module mod_cloud_gp
  use def_kintyp,              only   : ip,rp
  use def_domain,              only   : ndime
  use mod_rulepw,              only   : rulepw
  use def_quadrature,          only   : GAUSS_LEGENDRE_RULE
  use mod_elmgeo,              only   : elmgeo_shapf_deriv_heslo
  implicit none

  private

  public :: generate_cloud_gp

 contains

  subroutine generate_cloud_gp
     !-----------------------------------------------------------------------
     ! Do actions to generate and define a cloud of Gaussian points
     !-----------------------------------------------------------------------

     implicit none

     call int_rule_cloud_gp
     call allocate_memory_cloud_gp
     call initialize_var_cloud_gp
     call calc_coord_cloud_gp
  end subroutine generate_cloud_gp


  subroutine int_rule_cloud_gp
     !-----------------------------------------------------------------------
     ! Define integration rule
     !-----------------------------------------------------------------------
     use def_domain,    only :  nelem, nelty
     use def_domain,    only :  ltype, lrule, ngaus
     use def_domain,    only :  gp_total_cloud_gp, lrule_cloud_gp, &
                                ngaus_cloud_gp, mgaus_cloud_gp
     use def_coupli,    only :  ngaus_elem

     implicit none
     integer(ip)             :: ielty, pelty, ielem

     gp_total_cloud_gp = 0_ip
     mgaus_cloud_gp    = 1_ip

     do ielty = 1, nelty
        lrule_cloud_gp(ielty) = lrule(ielty)
        if (ngaus_elem(1) == 0_ip) then
           ngaus_cloud_gp(ielty) = ngaus(ielty)
        else
           ngaus_cloud_gp(ielty) = ngaus_elem(ielty)
        end if
        if (mgaus_cloud_gp < ngaus_cloud_gp(ielty))   mgaus_cloud_gp = ngaus_cloud_gp(ielty)
     end do


     do ielem = 1, nelem
        pelty = ltype(ielem)
        gp_total_cloud_gp = gp_total_cloud_gp + ngaus_cloud_gp(pelty)
     end do

  end subroutine int_rule_cloud_gp


  subroutine allocate_memory_cloud_gp
     !-----------------------------------------------------------------------
     ! Allocate memory for the different components of the structure elmar_cloud_gp
     !-----------------------------------------------------------------------
     use def_master,    only :  INOTEMPTY
     use mod_memory
     use def_domain,    only :  memor_dom
     use def_domain,    only :  elmar_cloud_gp, coord_cloud_gp, lrule_cloud_gp
     use def_domain,    only :  gp_total_cloud_gp
     use def_domain,    only :  nelty, lexis, nnode, ngaus_cloud_gp, ldime, ntens

     implicit none
     integer(ip)             :: ielty
     integer(ip)             :: pnode, pgaus, pdime, mdime, prule

     allocate(elmar_cloud_gp(nelty))
     do ielty = 1, nelty
        nullify(elmar_cloud_gp(ielty) % shape)
        nullify(elmar_cloud_gp(ielty) % deriv)
        nullify(elmar_cloud_gp(ielty) % heslo)
        nullify(elmar_cloud_gp(ielty) % weigp)
        nullify(elmar_cloud_gp(ielty) % posgp)
     end do
     nullify(coord_cloud_gp)

     do ielty = 1, nelty
        if( lexis(ielty) == 1_ip ) then
           pnode = nnode(ielty)
           pgaus = ngaus_cloud_gp(ielty)
           pdime = ldime(ielty)
           mdime = max(1_ip,pdime)
           prule = lrule_cloud_gp(ielty)
           elmar_cloud_gp(ielty) % pgaus = pgaus

           call memory_alloca(memor_dom,'SHAPE','mod_cloud_gp',elmar_cloud_gp(ielty) % shape,pnode,pgaus)
           call memory_alloca(memor_dom,'DERIV','mod_cloud_gp',elmar_cloud_gp(ielty) % deriv,mdime,pnode,pgaus)
           call memory_alloca(memor_dom,'HESLO','mod_cloud_gp',elmar_cloud_gp(ielty) % heslo,ntens,pnode,pgaus)
           call memory_alloca(memor_dom,'POSGP','mod_cloud_gp',elmar_cloud_gp(ielty) % posgp,ndime,pgaus)
           call memory_alloca(memor_dom,'WEIGP','mod_cloud_gp',elmar_cloud_gp(ielty) % weigp,pgaus)
         end if
     end do

     if (INOTEMPTY) then
        call memory_alloca(memor_dom,'GAUSS_COORD','mod_cloud_gp',coord_cloud_gp,ndime,gp_total_cloud_gp)
     end if

  end subroutine allocate_memory_cloud_gp

  
  subroutine initialize_var_cloud_gp
     !-----------------------------------------------------------------------
     ! This routine fills structure elmar_cloud_gp which contains the
     ! characteristics of the cloud_gp mesh for interpolation.
     !-----------------------------------------------------------------------
     use def_domain,    only :  nelty, lexis, ldime, nnode, ltopo
     use def_domain,    only :  elmar_cloud_gp, lrule_cloud_gp

     implicit none
     integer(ip)             :: ielty, ierro
     integer(ip)             :: pdime, pgaus, prule, pnode, pquad, ptopo

     ierro = 0_ip

     do ielty = 1,nelty
        if( lexis(ielty) == 1 ) then

           pdime = ldime(ielty)
           prule = lrule_cloud_gp(ielty)
           pnode = nnode(ielty)
           pgaus = elmar_cloud_gp(ielty) % pgaus
           ptopo = ltopo(ielty)
           pquad = GAUSS_LEGENDRE_RULE
           
           ! Position in the canonical element and weights for Gaussian points
           call rulepw(pdime,pgaus,ptopo,pquad,elmar_cloud_gp(ielty)%posgp, &
                   elmar_cloud_gp(ielty)%weigp,ierro)

           if (ierro == 1_ip)  call runend('DOMAIN CLOUD POINTS: not valid number of Gauss points')

           ! Shape functions, derivatives and hessian in Gaussian points
           call elmgeo_shapf_deriv_heslo(&
                pdime,pnode,pgaus,elmar_cloud_gp(ielty)%posgp,&
                elmar_cloud_gp(ielty)%shape,elmar_cloud_gp(ielty)%deriv,&
                elmar_cloud_gp(ielty)%heslo,ierro)
        end if
     end do
  end subroutine initialize_var_cloud_gp


  subroutine calc_coord_cloud_gp
     !-----------------------------------------------------------------------
     ! Find the global coordinates for the Gauss points of the cloud_gp mesh for
     ! interpolation from the shape functions and the nodal coordinates.
     !-----------------------------------------------------------------------
     use def_master,    only :  INOTEMPTY
     use def_domain,    only :  nelem, ltype, nnode, lnods
     use def_domain,    only :  elmar_cloud_gp, coord_cloud_gp, &
                                ngaus_cloud_gp
     use def_domain,    only :  coord

     implicit none
     integer(ip)             :: ielem, inode, ipoin, igaus, kgaus, count_gaus
     integer(ip)             :: pnode, pelty, pgaus
     real(rp)                :: elcod(ndime)

     if (INOTEMPTY) then
        count_gaus = 0_ip
        coord_cloud_gp = 0.0_rp

        ! Loop over all the elements
        elements: do ielem = 1, nelem
           pelty = ltype(ielem)
           pnode = nnode(pelty)
           pgaus = ngaus_cloud_gp(pelty)

           do inode = 1,pnode
              ipoin          = lnods(inode,ielem)
              elcod(1:ndime) = coord(1:ndime,ipoin)
              do igaus = 1, pgaus
                 kgaus = count_gaus + igaus
                 coord_cloud_gp(1:ndime,kgaus) = &
                    coord_cloud_gp(1:ndime,kgaus) + &
                    elmar_cloud_gp(pelty) % shape(inode,igaus) * elcod
              end do
           end do

           count_gaus = count_gaus + pgaus
        end do elements
     end if

  end subroutine calc_coord_cloud_gp

end module mod_cloud_gp
