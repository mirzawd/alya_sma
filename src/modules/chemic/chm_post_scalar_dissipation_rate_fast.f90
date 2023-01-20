!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_post_scalar_dissipation_rate_fast(VECTOR_DIM,pnode,pgaus,list_elements,ivari)
   !------------------------------------------------------------------------
   ! NAME
   !    chm_post_scalar_dissipation_rate
   ! DESCRIPTION
   !    Elemental operations for post-processing scalar dissipation rates
   !    in flamelet combustion model
   ! USES
   ! USED BY
   !    chm_outvar
   !***
   !------------------------------------------------------------------------
   use def_kintyp, only      : ip,rp
   use def_domain, only      : ndime,mnode,mgaus,ltype,   &
                               llapl,lorde,ltopo,   &
                               lnods,elmar,hnatu,ntens

   use def_chemic, only      : kfl_advec_chm,kfl_ellen_chm,     &
                               ADR_chm,nclas_chm, &
                               DtRho_gp,kfl_DtRho_tab_index_chm
   use mod_ker_proper, only  : ker_proper
   use mod_chm_element_operations_fast, only :chm_elmlen_fast, chm_elmchl_fast
   use mod_chm_element_operations_fast, only :chm_elmgac_flamLet_fast
   use mod_element_integration, only : element_shape_function_derivatives_jacobian


   implicit none
   integer(ip), intent(in)          :: VECTOR_DIM
   integer(ip), intent(in)          :: pnode                                    !< Number of nodes
   integer(ip), intent(in)          :: pgaus                                    !< Number of Gauss points
   integer(ip), intent(in)          :: list_elements(VECTOR_DIM)                !< List of elements
   integer(ip),  intent(in) :: ivari

   integer(ip)               :: list_elements_p(VECTOR_DIM)           ! List of elements (always positive)

   integer(ip) :: ielem,ivect,ielemone
   integer(ip) :: pelty
   integer(ip) :: plapl,porde,ptopo
   integer(ip) :: dummi

   real(rp)    :: elcon(VECTOR_DIM,pnode,nclas_chm,ADR_chm(1) % ntime)
   real(rp)    :: elcod(VECTOR_DIM,ndime,pnode)
   real(rp)    :: elvel(VECTOR_DIM,ndime,pnode)
   real(rp)    :: elmas(VECTOR_DIM,pnode,nclas_chm)
   real(rp)    :: gpvol(VECTOR_DIM,pgaus)
   real(rp)    :: gphco(VECTOR_DIM,pgaus)
   real(rp)    :: gpsph(VECTOR_DIM,pgaus)
   real(rp)    :: gpdtr(VECTOR_DIM,pgaus)
   real(rp)    :: gpden(VECTOR_DIM,pgaus)
   real(rp)    :: gptur(VECTOR_DIM,pgaus)
   real(rp)    :: gpcar(VECTOR_DIM,ndime,mnode,mgaus)
   real(rp)    :: gphes(VECTOR_DIM,ntens,mnode,mgaus)
   real(rp)    :: gp_res_Chi_Yc(VECTOR_DIM,pgaus)
   real(rp)    :: gp_sgs_Chi_Yc(VECTOR_DIM,pgaus)
   real(rp)    :: gp_res_Chi_z(VECTOR_DIM,pgaus)
   real(rp)    :: gp_sgs_Chi_z(VECTOR_DIM,pgaus)
   real(rp)    :: gp_Chi_out(VECTOR_DIM,pgaus)

   real(rp)    :: gpder(VECTOR_DIM,ndime,mnode,pgaus)              ! dNk/dsj
   real(rp)    :: gpsha(VECTOR_DIM,pnode,pgaus)                    ! N

   real(rp)    :: dummr(VECTOR_DIM,mgaus*ndime)
   real(rp)    :: chale(VECTOR_DIM,3),chave(VECTOR_DIM,3),hleng(VECTOR_DIM,3),tragl(VECTOR_DIM,9)

   integer(ip) :: lnods_loc(VECTOR_DIM,pnode)

   external    :: elmlen
   external    :: elmchl
   external    :: elmcar
   external    :: chm_calc_scalar_dissip
   external    :: chm_calc_scalar_dissip_fast
   external    :: chm_post_projection
   external    :: chm_post_gather

   !
   ! Initialization
   !
   do ivect = 1,VECTOR_DIM
       ielem = abs(list_elements(ivect))
       if( ielem /= 0 ) then
           list_elements_p(ivect)   = list_elements(ivect)
           lnods_loc(ivect,1:pnode) = lnods(1:pnode,ielem)
       else
           list_elements_p(ivect)   = list_elements(1)
           lnods_loc(ivect,1:pnode) = lnods(1:pnode,list_elements(1))
       end if
   end do

   !
   ! Element dimensions
   !

   ielemone = list_elements_p(1)

   pelty = ltype(ielemone)
   plapl = llapl(pelty)
   porde = lorde(pelty)
   ptopo = ltopo(pelty)

   !
   ! Initialization variables
   !
   gptur = 0.0_rp

   gp_res_Chi_Yc = 0.0_rp
   gp_sgs_Chi_Yc = 0.0_rp
   gp_res_Chi_z  = 0.0_rp
   gp_sgs_Chi_z  = 0.0_rp


    !
    ! Gather
    !
    call chm_elmgac_flamLet_fast(&
              VECTOR_DIM,pnode,lnods_loc,elcod,elcon,elvel,elmas)

    !
    ! CHALE, HLENG and TRAGL
    !
    call chm_elmlen_fast(&
        VECTOR_DIM,ndime,pnode,elmar(pelty)%dercg,tragl,elcod,hnatu(pelty),hleng)

    call chm_elmchl_fast(&
           VECTOR_DIM,tragl,hleng,elcod,dummr,chave,chale,pnode,&
           porde,hnatu(pelty),kfl_advec_chm,kfl_ellen_chm)


    !
    ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, GPVOL
    !
     call element_shape_function_derivatives_jacobian(&
    pnode,pgaus,plapl,elmar(pelty) % weigp,elmar(pelty) % shape,&
    elmar(pelty) % deriv,elmar(pelty) % heslo,&
    elcod,gpvol,gpsha,gpder,gpcar,gphes,list_elements=list_elements_p)


    call ker_proper('DENSI','PGAUS',dummi,list_elements_p,gpden,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM)
    call ker_proper('TURBU','PGAUS',dummi,list_elements_p,gptur,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM)

    !
    ! Dt*Rho on Gauss points:
    !
    gpdtr = 1.0_rp
    if (kfl_DtRho_tab_index_chm > 0) then
       do ivect = 1,VECTOR_DIM
          ielem = list_elements_p(ivect)
          if(ivect== 1 .or. (ivect > 1 .and. ielem /=ielemone)) then
             gpdtr(ivect,1:pgaus) = DtRho_gp(ielem) % a(1:pgaus,1,1)
          end if
       end do
    else
       call ker_proper('CONDU','PGAUS',dummi,list_elements_p,gphco,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM)
       call ker_proper('SPHEA','PGAUS',dummi,list_elements_p,gpsph,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM)
       gpdtr(:,1:pgaus) = gphco(:,1:pgaus) / gpsph(:,1:pgaus)
    endif

    call chm_calc_scalar_dissip_fast(&
           VECTOR_DIM,pnode,pgaus,elcon,elvel,elmar(pelty)%shape,gpcar,hleng,gptur,&
           gpdtr,gpden,gp_res_Chi_Yc,gp_sgs_Chi_Yc,gp_res_Chi_z,gp_sgs_Chi_z)



    select case (ivari)

      case (23_ip)
         gp_Chi_out(:,:) = gp_res_Chi_Yc(:,:)

      case (24_ip)
         gp_Chi_out(:,:) = gp_res_Chi_Z(:,:)

      case (25_ip)
         gp_Chi_out(:,:) = gp_sgs_Chi_Yc(:,:)

      case (26_ip)
         gp_Chi_out(:,:) = gp_sgs_Chi_Z(:,:)

    end select

    elements: do ivect = 1,VECTOR_DIM
     ielem = list_elements_p(ivect)
     if(ivect== 1 .or. (ivect > 1 .and. ielem /=ielemone)) then

      !
      ! Projection
      !
      call chm_post_projection(ivari,ielem,pnode,pgaus,elmar(pelty)%shape,gp_Chi_out(ivect,:),gpvol(ivect,:))
    end if
   end do elements

end subroutine chm_post_scalar_dissipation_rate_fast
