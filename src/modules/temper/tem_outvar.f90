!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_outvar(ivari,imesh)
  !------------------------------------------------------------------------
  !****f* Temper/tem_output
  ! NAME 
  !    tem_output
  ! DESCRIPTION
  !    Output a postprocess variable
  ! USES
  !    postpr
  !    memgen
  ! USED BY
  !    tem_output
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_temper
  use mod_gradie
  use mod_lodi_tem
  use def_kermod,          only : kfl_waexl_ker,kfl_vefun
  use mod_communications,  only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_finite_volume,   only : finite_volume_element_to_nodes
  use mod_ADR,             only : ADR_manufactured_nodal_error
  use mod_projec,          only : projec_elements_to_nodes
  use mod_tem_entropy,     only : tem_entropy_postprocess
  use mod_memory,          only : memory_alloca, memory_deallo
  use mod_solver,          only : solver_lumped_mass_system 
  use mod_outvar,          only : outvar
  use mod_arrays,          only : arrays_name
  implicit none
  integer(ip), intent(in) :: ivari
  integer(ip), intent(in) :: imesh   !< Mesh to postprocess on
  integer(ip)             :: ipoin,iline,jpoin,icont,idime
  real(rp)                :: dummr,rutim
  real(rp), pointer       :: tempeAux(:) 
  real(rp), pointer       :: grad_tempe(:,:) 

  rutim = cutim

  select case ( arrays_name(ivari) )  

  case ( 'TEMPE' )
     !
     ! Temperature
     !
     if( kfl_discr_tem == NODAL_SCHEME ) then
        gesca => tempe(:,1)
     else        
        if( INOTMASTER ) then
          call memgen(zero,npoin,zero)
         !  call finite_volume_element_to_nodes(meshe(ndivi),1_ip,tempe,gesca)
          tempeAux => tempe(:,1)
          call projec_elements_to_nodes(tempeAux,gesca)
        end if
     end if

  case ( 'HEATF' )
     !
     ! Heat flux
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call tem_outhfl()
     end if

  case ( 'TESGS' )
     !
     ! Tesgs
     !
     !if( INOTMASTER ) ger3p => ADR_tem % sgs 
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call projec_elements_to_nodes(ADR_tem % sgs,gesca)
     end if

  case ( 'ERROR' )
     !
     ! Error w/r manufactured solution
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call ADR_manufactured_nodal_error(ADR_tem,cutim,tempe,gesca)
     end if

  case ( 'AVTEM' )
     !
     ! Average temperature
     !
     if( rutim > avtim_tem ) then
        dummr = rutim - avtim_tem
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin = 1,npoin
              gesca(ipoin) =  avtem_tem(ipoin)/dummr
              avtem_tem(ipoin)=0.0_rp
           end do
        end if
     end if

  case ( 'VELOC' )
     !
     ! Velocity
     !
     if( kfl_advec_tem /= 0 ) then
        if( INOTEMPTY ) then        
           if( kfl_vefun == 0 ) then
              gevec => veloc(:,:,1)
           else 
              gevec => advec(:,:,1)
           end if
        end if
     else
        return
     end if

  case ( 'TURVI' )
     !
     ! Turbulent viscosity
     !
     if( INOTMASTER ) then
        if(size(turmu,1)>1) then
           gesca => turmu
        end if
     end if

  case ( 'RESID' )
     !
     ! Residual
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           gesca(ipoin) = therm(ipoin,1)-teold_tem(ipoin)
        end do
     end if
     rutim = real(ittot_tem,rp)

  case ( 'GROUP' )
     !
     ! GROUPS FOR DEFLATED CG
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,0_ip)
        do ipoin=1,npoin
           gesca(ipoin)=real(solve(1)%lgrou(ipoin),rp)
        end do
     end if

  case ( 'PROJ1' )
     !
     ! Projection
     !
     gesca => ADR_tem % proje1

  case ( 'LIMIT' )
     !
     ! Limiter
     !
     if( kfl_limit_tem /= 0 ) then
       call runend('TEM_OUTVAR: LIMITER NOT CODED') 
     end if

  case ( 'LINTE' )
     !
     ! LINTE: Linelets of preconditioner 
     !
     if( INOTMASTER ) then
        icont=0
        do ipoin=1,npoin
           rhsid(ipoin)=0.0_rp
        end do
        do iline=1,solve(1)%nline
           icont=icont+1
           do ipoin=solve(1)%lline(iline),solve(1)%lline(iline+1)-1
              jpoin=solve(1)%lrenup(ipoin)
              rhsid(jpoin)=real(icont,rp)
           end do
        end do
        gesca => rhsid
     end if

  case ( 'TESTS' )
     !
     ! TEST
     !
     if( INOTMASTER ) then
        call memgen(0_ip,ndime,npoin)
        call memgen(0_ip,npoin,0_ip)
        do ipoin = 1,npoin
           gesca(ipoin) = 2.0_rp*coord(1,ipoin)+3.0_rp*coord(2,ipoin)
        end do
        call grasca(gesca,gevec)
        call memgen(2_ip,npoin,0_ip)
     end if

  case ( 'WATVA' )
     !
     ! WAT_VAPOR: Water vapor
     !
     if( INOTMASTER ) then
        call memgen(0_ip,npoin,0_ip)
        call tem_poswat()
     end if

  case ( 'FIXTE' )
     !
     ! KFL_FIXNO_TEM
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           gesca(ipoin) = real(kfl_fixno_tem(1,ipoin),rp)
        end do
     end if

  case ( 'AVTE2' )
     !
     ! Average tempe*tempe
     !
     if( rutim > avtim_tem ) then
        dummr = rutim - avtim_tem
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin = 1,npoin
              gesca(ipoin) =  avte2_tem(ipoin)/dummr
              avte2_tem(ipoin)=0.0_rp
           end do
        end if
     end if

  case ( 'AVTEV' )
     !
     ! Average veloc*tempe
     !
     if( rutim > avtim_tem ) then
        dummr = rutim - avtim_tem
        if( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin = 1,npoin
              do idime=1, ndime
                 gevec(idime,ipoin) =  avtev_tem(idime,ipoin)/dummr
                 avtev_tem(idime,ipoin)=0.0_rp
              end do
           end do
        end if
     end if
     
  case ( 'AVDEN' )
     !
     ! Average temperature
     !
     if( rutim > avtim_tem ) then
        dummr = rutim - avtim_tem
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin = 1,npoin
              gesca(ipoin) =  avden_tem(ipoin)/dummr
              avden_tem(ipoin)=0.0_rp
           end do
        end if
     end if
     
  case ( 'FVVEL' )
     !
     ! Favre average velocity rho*veloc
     !
     if( rutim > avtim_tem ) then
        dummr = rutim - avtim_tem
        if( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin = 1,npoin
              do idime=1, ndime
                 gevec(idime,ipoin) =  fvvel_tem(idime,ipoin)/dummr
                 fvvel_tem(idime,ipoin)=0.0_rp
              end do
           end do
        end if
     end if

  case ( 'HEATN' )
     !
     ! Heat flux computed from matrix RHS
     ! 
     if( INOTMASTER ) then
        call memgen(zero,npoin,0_ip)
        if ( kfl_waexl_ker == 0) then 
           do ipoin = 1,npoin
              gesca(ipoin) = solve_sol(1) % reaction(1,ipoin) 
           end do
        !call PAR_INTERFACE_NODE_EXCHANGE(gesca,'SUM','IN MY CODE')
        else if ( kfl_waexl_ker /= 0) then ! MATIAS temporal condition
           call memory_alloca(mem_modul(1:2,modul),'HEATF_TEM','tem_outvar',heatf_tem,npoin)
           call memory_alloca(mem_modul(1:2,modul),'MASSB_TEM','tem_outvar',massb_tem,npoin)
           
           call tem_boundary_heat_flux() ! loads variational heat flux

           gesca(1:npoin) = heatf_tem(1:npoin) 

           call memory_deallo(mem_modul(1:2,modul),'MASSB_TEM','tem_outvar',massb_tem)
           call memory_deallo(mem_modul(1:2,modul),'HEATF_TEM','tem_outvar',massb_tem)
        end if
     end if 

  case ( 'GRATE' )
     !
     ! grad(T)
     !
     if( INOTMASTER ) then
        call memgen(zero,ndime,npoin)
        if( kfl_discr_tem == NODAL_SCHEME ) then
           call grasca(tempe,gevec)
        else
           nullify(grad_tempe)
           call memory_alloca(mem_modul(1:2,modul),'GRAD_TEMPE','tem_outvar',grad_tempe,ndime,nelem)
           !select case ( my_grad ) 
           !case ( fv_grad_method_gauss )
              call tem_gradient_gauss(tempe,grad_tempe)
           !case (fv_grad_method_ls)
           !   call tem_gradient_ls(tempe,grad_tempe)
           !case default
            !  call tem_gradient_gauss(tempe,grad_tempe) ! is the one more robust
           !end select
           call projec_elements_to_nodes(grad_tempe,gevec)
           call memory_deallo(mem_modul(1:2,modul),'GRAD_TEMPE','tem_outvar',grad_tempe)
        end if
     end if

  case ( 'CHARA' )
     
     if( INOTMASTER ) then
        !call memgen(zero,npoin,zero)
        !call ker_proper('DENSI','NPOIN',dummi,dummi,gesca)
        !gesca(1:npoin) = prthe(1)/(gasco*tempe(1:npoin,1))
        !CHARAs%detem => gesca 

       !CHARAs%gamme = gasco
        call lodi_tem_allocate( CHARAs )
        call lodi_tem_get_characteristics( CHARAs )
        !
        call memgen(zero,ndime,npoin) 
        do idime = 1,ndime
          do ipoin = 1,npoin
            gevec(idime,ipoin) = CHARAs%chrc(CHARAs%idofn,ipoin,idime) ! chrc_tem(ndofn_tem,npoin,ndime) 
          enddo
        enddo
        !
        call lodi_tem_deallocate( CHARAs )
     end if

  case ( 'BVNAT' )
     !
     ! Heat flux computed from matrix RHS
     ! 
     if( INOTMASTER ) then
        call memgen(zero,npoin,0_ip)
        do ipoin = 1,npoin
           gesca(ipoin) = solve_sol(1) % bvnat(1,ipoin) 
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(gesca,'SUM','IN MY CODE')
     end if 

  case ( 'ENTHA' )
     !
     ! Ethalpy
     !
     gesca => therm(:,1)


  case ( 'REACT' )
     !
     ! 'REACT'
     !
     if( INOTMASTER ) then
       if(.not.associated(solve_sol(1)%lpoin_reaction) )  call runend('ERROR: -->POSTPROCESS REACT<-- ')
       call memgen(zero,npoin,0_ip)
       gesca(1:npoin) = -1.0_rp
       where( solve_sol(1) % lpoin_reaction(1:npoin) ) gesca(1:npoin) = 1.0_rp
     endif 

  case ( 'PROJ2' )
     !
     ! Projection
     !
     if( associated(ADR_tem % proje2) ) then
        gesca => ADR_tem % proje2
     else
        return
     end if

  case ( 'RESHE' )
     !
     ! 'RESHE' Heat flux interpolated at the nodes computed from Residuals
     !
     if( INOTMASTER ) then
       call memgen(zero,npoin,zero)
       gesca(1:npoin) =  solve(1)%reaction(1,1:npoin)
     endif

  case ( 'AVRES' )
     !
     ! Average heat flux interpolated at the nodes computed from Residuals
     !
     if( rutim > avtim_tem ) then
        dummr = rutim - avtim_tem
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin = 1,npoin
              gesca(ipoin) =  avres_tem(ipoin)/dummr
              avres_tem(ipoin)=0.0_rp
           end do
        end if
     end if

  case ( 'FVTEM' )
     !
     ! FVTEM 
     !
     if( rutim > avtim_tem ) then
        dummr = rutim - avtim_tem
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin = 1,npoin
              gesca(ipoin) =  fvtem_tem(ipoin)/dummr
              fvtem_tem(ipoin)=0.0_rp
           end do
        end if
     end if
     
  case ( 'FVTE2' )
     !
     ! FVTE2
     !
     if( rutim > avtim_tem ) then
        dummr = rutim - avtim_tem
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin = 1,npoin
              gesca(ipoin) =  fvte2_tem(ipoin)/dummr
              fvte2_tem(ipoin)=0.0_rp
           end do
        end if
     end if
     
  case ( 'ENVIT' )
     !
     ! Entropy viscosity for the Temperature
     !
     call memgen(zero,max(npoin,1_ip),zero)
     call tem_entropy_postprocess(gesca)


  case ( 'HEASK' )
     !
     ! heat source from partis 
     !
     call memgen(zero,max(1_ip,npoin),zero)
     if (associated(heat_sink)) then
         do ipoin=1,npoin
            gesca(ipoin) = heat_sink(ipoin)
         enddo
        call solver_lumped_mass_system(1_ip,gesca,EXCHANGE=.false.)
     else
        do ipoin=1,npoin
            gesca(ipoin) = 0.0_rp 
        end do
     endif

  case ( 'AVHSK' )
     !
     ! Average heat source from partis
     !
     if( rutim > avtim_tem ) then
        dummr = rutim - avtim_tem
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin = 1,npoin
              gesca(ipoin) =  avhsk_tem(ipoin)/dummr
              avhsk_tem(ipoin)=0.0_rp
           end do
        end if
     end if

  case ( 'HEATS' )
     !
     ! Heat nodal source
     !
     call memgen(zero,npoin,zero)
     do ipoin = 1,npoin
        gesca(ipoin) = heat_source(1,ipoin)
     end do 
     call PAR_INTERFACE_NODE_EXCHANGE(gesca,'SUM','IN MY CODE')

  end select


  call outvar(&
       ivari,&
       ittim,rutim,postp(1) % wopos(:,ivari),MESH_ID=imesh)

end subroutine tem_outvar
