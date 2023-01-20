!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_outvar.f90
!> @author  Guillaume Houzeaux
!> @date    20/02/2013
!> @brief   Output a postprocess variable
!> @details Output a postprocess variable for nastin.
!> @} 
!-----------------------------------------------------------------------
subroutine nsi_outvar(jvari,imesh)

  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_nastin
  use mod_memchk
  use mod_lodi_nsi 
  use mod_ker_proper
  use mod_communications,      only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_projec,              only : projec_elements_to_nodes
  use mod_solver,              only : solver_solve
  use mod_solver,              only : solver_postprocess
  use mod_solver,              only : solver_lumped_mass_system
  use mod_ADR,                 only : ADR_assemble_extension
  use mod_ADR,                 only : ADR_assemble_laplacian
  use mod_memory,              only : memory_alloca
  use mod_memory,              only : memory_deallo
  use mod_memory,              only : memory_size
  use mod_projec,              only : projec_boundaries_to_nodes
  use mod_arrays,              only : arrays
  use mod_outvar,              only : outvar 
  use mod_maths_arrays,        only : maths_findloc
  use mod_frivel,              only : frivel
  use def_coupli,              only : coupling_type
  use def_coupli,              only : kfl_efect
  use mod_arrays,              only : arrays_number
  use mod_arrays,              only : arrays_name
  use mod_nsi_eigen_time_step, only : nsi_eigen_time_step_all
  use mod_streamlines,         only : strfun
  implicit none
  integer(ip), intent(in) :: jvari   !< 1 for veloc , 2 press , etc...
  integer(ip), intent(in) :: imesh   !< Mesh to postprocess on
  integer(ip)             :: ibopo,ivari
  integer(ip)             :: idime,ipoin,iline,icont,jpoin
  integer(ip)             :: izdom,jzdom,jdime,iboun
!  integer(4)              :: istat
  real(rp)                :: vmodu,xsoun,sofac,xfact,h,rho(2),mu(2),u,auxi
  real(rp)                :: tveno,delta_aux,vikin,gesma(3)
  real(rp),    pointer    :: geve2(:,:)
  real(rp)                :: rot_matrix(ndime,ndime)

  integer(ip)             :: dummi
  integer(ip), pointer    :: kfl_fixno_tmp(:,:)
  logical(lg)             :: kfl_intfo
  real(rp),    pointer    :: bvess_tmp(:,:)
  real(rp),    pointer    :: normal_tmp(:,:)
  
  nullify(geve2)
  !
  ! Define postprocess variable
  !
  ibopo = 0
  ivari = abs(jvari)
 
  select case ( arrays_name(ivari) )  

  case ( 'VELOC' )
     !
     ! VELOC: Velocity
     !     
     if( INOTMASTER ) then        
        if( kfl_inert_nsi == 0 ) then
           gevec => veloc(:,:,1)
        else
           !
           ! u <= u + w x r
           !
           call memgen(zero,ndime,npoin)
           call rotmat(ndime,cutim*fvnoa_nsi,fvdia_nsi,rot_matrix)
           do ipoin = 1,npoin
              if( ndime == 2 ) then
                 gesma(1) = veloc(1,ipoin,1)&
                      - fvela_nsi(3) * coord(2,ipoin)            ! u-wz*y
                 gesma(2) = veloc(2,ipoin,1)&  
                      + fvela_nsi(3) * coord(1,ipoin)            ! v+wz*x
              else
                 gesma(1) = veloc(1,ipoin,1)&
                      - fvela_nsi(3) * coord(2,ipoin)&           ! u-wz*y
                      + fvela_nsi(2) * coord(3,ipoin)            !  +wy*z
                 gesma(2) = veloc(2,ipoin,1)&
                      + fvela_nsi(3) * coord(1,ipoin)&           ! v+wz*x
                      - fvela_nsi(1) * coord(3,ipoin)            !  -wx*z
                 gesma(3) = veloc(3,ipoin,1)&
                      + fvela_nsi(1) * coord(2,ipoin)&           ! w+wx*y
                      - fvela_nsi(2) * coord(1,ipoin)            !  -wy*x
              end if
              do idime=1,ndime
                 gevec(idime,ipoin) = 0.0_rp              
                 do jdime=1,ndime
                    gevec(idime,ipoin) = gevec(idime,ipoin) + rot_matrix(idime,jdime)*gesma(jdime)
                 end do
              end do
           end do
        end if
     end if

  case ( 'PRESS' )
     !
     ! PRESS: Pressure
     !
     if( kfl_regim_nsi == 2 ) then
        if( INOTMASTER ) then
           do ipoin = 1,npoin
              rhsid(ipoin) = densi(ipoin,1)*gasco*tempe(ipoin,1)
           end do
           gesca => rhsid
        end if
     else
        if( abs(press_lev_nsi) > 0.0_rp ) then
           call memgen(zero,npoin,zero)
           do ipoin = 1,npoin
              gesca(ipoin) = press(ipoin,1) + press_lev_nsi
           end do
        else
           gesca => press(:,1) 
        end if
     end if

  case ( 'STREA' )
     !
     ! STREA: Streamlines 
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call strfun(veloc,gesca)
     end if

  case ( 'RESID' )
     !
     ! RESID: Velocity residual
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
     end if
     call nsi_outres()

  case ( 'VESGS' )
     !
     ! VESGS: Velocity subgrid scale
     !
     if( INOTMASTER ) ger3p => vesgs(:) 

  case ( 'BOUND' )
     !
     ! BOUND: Boundary conditions
     !
     if( INOTMASTER ) then 
        call memgen(zero,ndime,npoin)
        do ipoin=1,npoin
           do idime=1,ndime
              gevec(idime,ipoin)=unkno((ipoin-1)*ndime+idime)
           end do
        end do
     end if

  case ( 'DENSI' )
     !
     ! DENSI: Density
     !
     if( INOTMASTER ) then
        if ( kfl_prope /= 0_ip ) call runend ('DENSITY must be output by kernel not nsi')
        call memgen(zero,npoin,zero)
        if (kfl_regim_nsi == 3) then !Low Mach regime computes its own density from temperature
           do ipoin = 1,npoin
              gesca(ipoin) = prthe(1)/(gasco*tempe(ipoin,1))
           end do
        else
           do ipoin = 1,npoin
              gesca(ipoin) = prope_nsi(1,ipoin)
           end do
        endif
     end if

  case ( 'PMV  ' )
     !
     ! PMV
     !
     if( INOTMASTER ) then
        call nsi_outpmv(1_ip,1_ip,npoin,rhsid)
        gesca => rhsid
     end if

  case ( 'PPD  ' )
     !
     ! PPD
     !
     if( INOTMASTER ) then
        call nsi_outpmv(2_ip,1_ip,npoin,rhsid)
        gesca => rhsid
     end if

  case ( 'MACHN' )
     !
     ! MACHN: Mach number
     !
     if( INOTMASTER ) then
        sofac=sqrt(gamth_nsi*gasco)
        do ipoin=1,npoin
           xsoun = sofac*sqrt(abs(tempe(ipoin,1)))
           if (kfl_coupl(ID_NASTIN,ID_CHEMIC) > 0 ) xsoun = xsoun / sqrt(abs(wmean(ipoin,1)))
           vmodu = 0.0_rp
           do idime=1,ndime
              vmodu=vmodu+veloc(idime,ipoin,1)*veloc(idime,ipoin,1)
           end do
           vmodu = sqrt(vmodu)
           rhsid(ipoin) = vmodu/xsoun
        end do
        gesca => rhsid
     end if

  case ( 'TANGE' )
     !
     ! TANGE: Tangential stress
     !
     if( INOTMASTER ) then
        call memgen(zero,ndime,npoin)
        call nsi_outtan(0_ip)    
     end if
     !ibopo=1

     !!case (12_ip)
     !!   !
     !!   ! VORTI
     !!   !
     !!   if(ndime==2) then
     !!      if( INOTMASTER ) then
     !!         call memgen(zero,ndime+1_ip,npoin)
     !!         vorti => gevec
     !!         call vortic(1_ip)
     !!         do ipoin=1,npoin
     !!            rhsid(ipoin)=sqrt(&
     !!                 &  vorti(1,ipoin)*vorti(1,ipoin)&
     !!                 & +vorti(2,ipoin)*vorti(2,ipoin))
     !!         end do
     !!         gesca => rhsid           
     !!      end if
     !!   else
     !!      print*, 'call case 12'
     !!      if( INOTMASTER ) then
     !!         call memgen(zero,ndime,npoin)
     !!         allocate(geve2(ndime+1,npoin),stat=istat)
     !!         call memchk(zero,istat,mem_modul(1:2,modul),'GEVE2','nsi_outvar',geve2)
     !!         nullify(geve2)
     !!         vorti => geve2
     !!         call vortic(1_ip)
     !!         do ipoin = 1,npoin
     !!            do idime = 1,ndime
     !!               gevec(idime,ipoin) = geve2(idime,ipoin)
     !!            end do
     !!         end do
     !!         deallocate(geve2,stat=istat)
     !!      end if
     !!   end if
     !!   nullify(vorti)

  case ( 'VORTI' )
     !
     ! VORTI
     !
     if( INOTMASTER ) then
        call memgen(zero,ndime,npoin)
        nullify(geve2)
        call memory_alloca(mem_modul(1:2,modul),'GEVE2','nsi_outvar',geve2,ndime+1_ip,npoin)
        vorti => geve2
        call vortic(1_ip)
        do ipoin = 1,npoin
           gevec(1:ndime,ipoin) = geve2(1:ndime,ipoin)
        end do
        call memory_deallo(mem_modul(1:2,modul),'GEVE2','nsi_outvar',geve2)
     end if
     nullify(vorti)

  case ( 'MODVO' )
     !
     ! MODVO
     !
     if( INOTMASTER ) then
        call memgen(zero,ndime+1_ip,npoin)
        vorti => gevec
        call vortic(1_ip)
        do ipoin=1,npoin
           rhsid(ipoin)=0.0_rp
           do idime=1,ndime
              rhsid(ipoin)=rhsid(ipoin)+vorti(idime,ipoin)*vorti(idime,ipoin)
           end do
           rhsid(ipoin)=sqrt(rhsid(ipoin))
        end do
        gesca => rhsid
     end if
     nullify(vorti)

  case ( 'LAMB2' )
     !
     ! LAMB2
     !
     if( INOTMASTER ) then
        call memgen(zero,ndime+1_ip,npoin)
        call memgen(zero,npoin,0_ip)
        vorti => gevec
        call vortic(4_ip)
        do ipoin = 1,npoin
           gesca(ipoin) = gevec(ndime+1_ip,ipoin)
        end do
     end if
     nullify(vorti)

  case ( 'GROUP' )
     !
     ! GROUP: GROUPS FOR DEFLATED CG
     !
     if( INOTMASTER ) then
        do ipoin=1,npoin
           rhsid(ipoin)=real(solve(2)%lgrou(ipoin),rp)
        end do
        gesca => rhsid
     end if

  case ( 'LINEL' )
     !
     ! LINEL: Linelets of preconditioner CG
     !
     if( INOTMASTER ) then
        do ipoin=1,npoin
           rhsid(ipoin)=0.0_rp
        end do
        do iline=1,solve(2)%nline
           do ipoin=solve(2)%lline(iline),solve(2)%lline(iline+1)-1
              jpoin=solve(2)%lrenup(ipoin)
              rhsid(jpoin)=real(iline,rp)
           end do
        end do
        gesca => rhsid
     end if

  case ( 'VISCO' )
     !
     ! VISCO: Viscosity
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           gesca(ipoin) = prope_nsi(2,ipoin)
        end do
     end if

  case ( 'FIXPR' )
     ! 
     ! FIXPR: KFL_FIXPR_NSI
     !
     if( NSI_SCHUR_COMPLEMENT .or. NSI_FRACTIONAL_STEP ) then
        if( INOTMASTER ) then
           do ipoin=1,npoin
              if(lpoty(ipoin)/=0) then
                 if( kfl_fixpr_nsi(1,ipoin) > 0 ) then
                    rhsid(ipoin)=1.0_rp
                 else
                    rhsid(ipoin)=0.0_rp
                 end if
              else
                 rhsid(ipoin)=0.0_rp                 
              end if
           end do
           gesca => rhsid
        end if
     end if

  case ( 'FIXNO' )
     !
     ! FIXNO: KFL_FIXNO_NSI
     !
     if( INOTMASTER ) then
        call memgen(zero,ndime,npoin)
        do ipoin=1,npoin
           do idime=1,ndime
              gevec(idime,ipoin)=real(kfl_fixno_nsi(idime,ipoin),rp)
           end do
        end do
     end if

  case ( 'TAU  ' )
     !
     ! TAU
     !
     call runend('NSI_OUTVAR: NOT CODED')

  case ( 'AVVEL' )
     !
     ! AVVEL_NSI: averaged velocity
     !
     if(cutim>avtim_nsi) then
        auxi = cutim - avtim_nsi
        if( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin=1,npoin
              do idime=1,ndime
                 gevec(idime,ipoin) = avvel_nsi(idime,ipoin) / auxi
              end do
           end do
        end if
     else
        return
     end if

  case ( 'SCHUR' )
     !
     ! SCHUR: Schur complement residual
     !
     if( INOTMASTER ) gesca => resch_nsi

  case ( 'VEPRO' )
     !
     ! VEPRO: Velocity projection
     !
     gevec => vepro_nsi 

  case ( 'PRPRO' )
     !
     ! PRPRO: Pressure projection
     !
     gesca => prpro_nsi 

  case ( 'YPLUS' )
     ! 
     ! YPLUS: Dimensionless distance to the wall y+
     !  
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call nsi_wyplus()
!!$        ! Extend exterior normals
!!$        call memgen(one, npoin,zero)
!!$        do iboun = 1,nboun
!!$           if( kfl_fixbo_nsi(iboun) == 4 ) then
!!$              gisca(lnodb(1:lnnob(iboun),iboun)) = 1
!!$           end if
!!$        end do
!!$        call PAR_INTERFACE_NODE_EXCHANGE(gisca,'MAX')
!!$        allocate( kfl_fixno_tmp(1,npoin) )
!!$        allocate( bvess_tmp(1,npoin) )
!!$         do ipoin = 1,npoin
!!$            if( gisca(ipoin) == 1 ) then
!!$               kfl_fixno_tmp(1,ipoin) = 1
!!$               gesca(ipoin)           = coord(1,ipoin)
!!$               bvess_tmp(1,ipoin)     = gesca(ipoin)
!!$               unkno(ipoin)           = gesca(ipoin)
!!$            else
!!$               kfl_fixno_tmp(1,ipoin) = 0
!!$               gesca(ipoin)           = 0.0_rp
!!$               bvess_tmp(1,ipoin)     = 0.0_rp 
!!$               unkno(ipoin)           = 0.0_rp 
!!$            end if
!!$        end do
!!$        solve_sol                => momod(modul) % solve(7:)
!!$        solve_sol(1) % kfl_fixno => kfl_fixno_tmp 
!!$        solve_sol(1) % bvess     => bvess_tmp   
!!$        solve_sol(1) % kfl_iffix =  1
!!$        call ADR_assemble_laplacian(meshe(ndivi),elmar,amatr,unkno,rhsid) 
     end if
!!$     call solver_solve(momod(modul) % solve,amatr,rhsid,unkno)    
!!$     if( INOTMASTER ) then
!!$        gesca(1:npoin) = unkno(1:npoin)
!!$        call memgen(three,npoin,zero)
!!$     end if

  case ( 'LIMIT' )
     !
     ! LIMIT: Limiter 
     !
     if( kfl_stabi_nsi == 2 ) then
        if( INOTMASTER ) then
           do ipoin = 1,npoin
              rhsid(ipoin) = 0.0_rp
           end do
           call nsi_elmope_omp(5_ip)
           call rhsmod(1_ip,rhsid)
           do ipoin = 1,npoin
              rhsid(ipoin) = rhsid(ipoin) / vmass(ipoin)
           end do
           gesca => rhsid
        end if
     end if

  case ( 'GRPRO' )
     !
     ! GRPRO: Pressure gradient projection
     !
     gevec => grpro_nsi 

  case ( 'AVPRE' )
     !
     ! AVPRE_NSI: averaged pressure
     !
     if(cutim>avtim_nsi) then
        auxi = cutim - avtim_nsi
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin=1,npoin
              gesca(ipoin) = avpre_nsi(ipoin) / auxi
           end do
        end if
     else
        return
     end if

  case ( 'VEERR' )
     !
     ! Velocity error w/r exact solution
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call nsi_exaerr(3_ip)
     end if

  case ( 'PECLE' )
     !
     ! PECLET: Pe = rho*u*h/(2*mu)
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        xfact = 1.0_rp / real(ndime,rp)
        do ipoin = 1,npoin
           call ker_proper('DENSI','IPOIN',ipoin,dummi,rho)
           call ker_proper('VISCO','IPOIN',ipoin,dummi,mu)
           h     = vmass(ipoin) ** xfact 
           u     = 0.0_rp
           do idime = 1,ndime
              u = u + veloc(idime,ipoin,1)*veloc(idime,ipoin,1)
           end do
           u = sqrt(u)
           gesca(ipoin) = 0.5_rp * rho(1) * u * h / mu(1)
        end do
     end if

  case ( 'HYDRO' )
     !
     ! Hydrostatic pressure
     !
     if( INOTMASTER ) then        
        gesca => bpess_nsi(1,:,1)
     end if

  case ( 'DUDT ' )
     !
     ! du/dt
     !
     if( INOTMASTER ) then
        call memgen(zero,ndime,npoin)
        call nsi_outdt2(1_ip,ndime,veloc)
     end if

  case ( 'D2UDT' )
     !
     ! d^2u/dt^2
     !
     if( INOTMASTER ) then
        call memgen(zero,ndime,npoin)
        call nsi_outdt2(2_ip,ndime,veloc)
     end if

  case ( 'DT1  ' )
     !
     ! Dt criterion 1
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call nsi_outdt2(3_ip,ndime,veloc)
     end if

  case ( 'DT2  ' )
     !
     ! Dt criterion 2
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call nsi_outdt2(4_ip,ndime,veloc)
     end if

  case ( 'PRERR' )
     !
     ! Pressure error w/r exact solution
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call nsi_exaerr(4_ip)
     end if

  case ( 'LINVE' )
     !
     ! LINEL: Linelets of preconditioner CG
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

  case ( 'DEFOR' )
     !
     ! DEFOR: deformaiton of the mesh
     !
     if( INOTMASTER ) then
        call memgen(zero,ndime,npoin)
        do ipoin = 1,npoin
           do idime = 1,ndime
              gevec(idime,ipoin) = 0.0_rp
           end do
        end do
        do ipoin = 1,npoin
           do idime = 1,ndime
              gevec(idime,ipoin) = dispm(idime,ipoin,1) 
           end do
        end do
     end if

  case ( 'NODPR' )
     !
     ! NODPR
     !
     if( INOTMASTER ) then
        do ipoin = 1,npoin
           rhsid(ipoin) = 0.0_rp
        end do
        if( nodpr_nsi > 0 ) rhsid( nodpr_nsi ) = 1.0_rp
        gesca => rhsid
     end if

  case ( 'MODVE' )
     !
     ! MODVE: module of velocity vector
     !
     if( INOTMASTER ) then
        do ipoin=1,npoin
           rhsid(ipoin)=0.0_rp
           do idime=1,ndime
              rhsid(ipoin)=rhsid(ipoin)+veloc(idime,ipoin,1)*veloc(idime,ipoin,1)
           end do
           rhsid(ipoin)=sqrt(rhsid(ipoin))
        end do
        gesca => rhsid
     end if

  case ( 'AVTAN' )
     !
     ! AVTAN_NSI: averaged tangential force
     !
     if(cutim>avtim_nsi) then
        auxi = cutim - avtim_nsi
        if( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin=1,npoin
              do idime=1,ndime
                 gevec(idime,ipoin) = avtan_nsi(idime,ipoin) / auxi
              end do
           end do
        end if
     else
        return
     end if

  case ( 'VELOM' )
     !
     ! VELOM
     !
     if( INOTMASTER ) gevec => velom(:,:)

  case ( 'FORCF' )
     !
     ! VELOM
     !
     if( INOTMASTER ) gevec => forcf

  case ( 'INTFO' )
     !
     ! INTFO: internal force
     !
     if( kfl_intfo_nsi /= 2 ) then !return  !< no info in time step==0, !< 2016ABRIL04   
        call memgen(zero,ndime,npoin)
     else   
        if( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin = 1,npoin
              do idime = 1,ndime
                 kfl_intfo = .false.
                 if( lpoty(ipoin) > 0 ) then ! If it's a boundary node
                    kfl_intfo = .true.
                 else if( kfl_efect ) then   ! Or if using Embedded Finite Element Coupling Technique
                    kfl_intfo = (coupling_type(2) % wet % kfl_fringe_wetnodes(ipoin) == 1) .and. &
                                (kfl_fixno_nsi(idime,ipoin) > 0)
                 end if
                 ! Output INTFO
                 if( kfl_intfo ) then
                    gevec(idime,ipoin) = intfo_nsi(ipoin) % bu(idime)
                    jzdom = 0
                    do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                       jpoin = c_dom(izdom)
                       jzdom = jzdom + 1
                       do jdime = 1,ndime
                          gevec(idime,ipoin) = gevec(idime,ipoin) - intfo_nsi(ipoin) % Auu(jdime,idime,jzdom) * veloc(jdime,jpoin,1)
                       end do
                       gevec(idime,ipoin) = gevec(idime,ipoin) - intfo_nsi(ipoin) % Aup(idime,jzdom) * press(jpoin,1)
                    end do
                 end if
              end do
           end do
           call rhsmod(ndime,gevec)
        end if
     end if ! kfl_intfo_nsi /= 2
     
  case ( 'USTAR' )
     !
     ! USTAR
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,0_ip)
        if( kfl_rough > 0 ) & ! variable roughness
             rough_dom = rough(ipoin)
        call runend('NSI_OUTVAR: RECODE WITH PROPERTIES')
        do ipoin = 1,npoin
           if( lpoty(ipoin) > 0 ) then
              if( kfl_delta==1) then
                 delta_aux = ywalp(lpoty(ipoin))
              else
                 delta_aux = delta_nsi  
              end if
              !vikin = visco_nsi(1,1)/densi_nsi(1,1)
              tveno = 0.0_rp
              do idime = 1,ndime
                 tveno = tveno + veloc(idime,ipoin,1) * veloc(idime,ipoin,1) 
              end do
              tveno = sqrt(tveno)
              call frivel(kfl_ustar,delta_aux,rough_dom,tveno,vikin,gesca(ipoin))
           end if
        end do
     end if

  case ( 'TRACT' )
     !
     ! TRACT: Traction ! the same as TTRAC
     !
     if( INOTMASTER ) then
        call memgen(zero,ndime,npoin)
        call nsi_outtra(1_ip)
     end if

  case ( 'FIXRS' )
     !
     ! KFL_FIXRS_NSI
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,0_ip)
        do ipoin = 1,npoin
           if( lpoty(ipoin) > 0 ) then
              gesca(ipoin) = real(kfl_fixrs_nsi(ipoin),rp)
           else
              gesca(ipoin) = 0.0_rp 
           end if
        end do
     end if

  case ( 'FVELO' )
     !
     ! FVELO: Velocity for filter
     !
     if( INOTMASTER ) then
        if( kfl_inert_nsi == 0 ) then
           gevec => veloc(:,:,1) 
        else
           call memgen(zero,ndime,npoin)
           if( ndime == 2 ) then
              do ipoin = 1,npoin
                 gevec(1,ipoin) = veloc(1,ipoin,1)&
                      - fvela_nsi(3) * coord(2,ipoin)            ! u-wz*y
                 gevec(2,ipoin) = veloc(2,ipoin,1)&  
                      + fvela_nsi(3) * coord(1,ipoin)            ! v+wz*x
              end do
           else
              do ipoin = 1,npoin
                 gevec(1,ipoin) = veloc(1,ipoin,1)&
                      - fvela_nsi(3) * coord(2,ipoin)&           ! u-wz*y
                      + fvela_nsi(2) * coord(3,ipoin)            !  +wy*z
                 gevec(2,ipoin) = veloc(2,ipoin,1)&
                      + fvela_nsi(3) * coord(1,ipoin)&           ! v+wz*x
                      - fvela_nsi(1) * coord(3,ipoin)            !  -wx*z
                 gevec(3,ipoin) = veloc(3,ipoin,1)&
                      + fvela_nsi(1) * coord(2,ipoin)&           ! w+wx*y
                      - fvela_nsi(2) * coord(1,ipoin)            !  -wy*x
              end do
           end if
        end if
     end if

  case ( 'FPRES' )
     !
     ! FPRES: Pressure for filter 
     !
     if( kfl_regim_nsi == 2 ) then
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin = 1,npoin
              rhsid(ipoin) = densi(ipoin,1)*gasco*tempe(ipoin,1)
           end do
           gesca => rhsid
        end if

     else 
        gesca => press(:,1) 
     end if

  case ( 'FTANG' )
     !
     ! FTANG: Tangential stress for filter
     !
     if( INOTMASTER ) then
        call memgen(zero,ndime,npoin)
        call nsi_outtan(1_ip)  ! for no slip wall law, variational forces are obtained - comes in gevec
     end if

  case ( 'AVVE2' )
     !
     ! AVVE2_NSI: averaged velocity**2
     !
     if(cutim>avtim_nsi) then
        auxi = cutim - avtim_nsi
        if( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin=1,npoin
              do idime=1,ndime
                 gevec(idime,ipoin) = avve2_nsi(idime,ipoin) / auxi
              end do
           end do
        end if
     else
        return
     end if

  case ( 'AVVXY' )
     !
     ! AVVXY_NSI: averaged vx * vy
     !
     if(cutim>avtim_nsi) then
        auxi = cutim - avtim_nsi
        if( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin=1,npoin
              do idime=1,ndime
                 gevec(idime,ipoin) = avvxy_nsi(idime,ipoin) / auxi
              end do
           end do
        end if
     else
        return
     end if

  case ( 'AVPR2' )
     !
     ! AVPR2_NSI: averaged pressure**2
     !
     if(cutim>avtim_nsi) then
        auxi = cutim - avtim_nsi
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin=1,npoin
              gesca(ipoin) = avpr2_nsi(ipoin) / auxi 
           end do
        end if
     else
        return
     end if

  case ( 'MESHR' )
     !
     ! Mesh rotation
     !
     if( INOTMASTER ) then
        call memgen(zero,ndime,npoin)           
        if( kfl_inert_nsi == 0 ) then
           gevec = 0.0_rp
        else
           !
           ! Compute rotation matrix
           !
           call rotmat(ndime,cutim*fvnoa_nsi,fvdia_nsi,rot_matrix)
           do ipoin = 1,npoin
              do idime=1,ndime
                 gevec(idime,ipoin)=0.0_rp              
                 do jdime=1,ndime
                    gevec(idime,ipoin) = gevec(idime,ipoin) + rot_matrix(idime,jdime)*coord(jdime,ipoin) ! Rotated coordinate
                 enddo
                 gevec(idime,ipoin) =  gevec(idime,ipoin) - coord(idime,ipoin) ! But we must save the displacement
              end do
           enddo
        end if
     end if

  case ( 'NODEF' )
     !
     ! NODE_FORCE
     !
     if( INOTMASTER ) then
        call memgen(zero,ndime,npoin)   
        do ipoin = 1,npoin
           gesca(ipoin) = real(intfo_nsi(ipoin) % kfl_exist,rp)              
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(gesca,'MAX')
     end if

  case ( 'MOMEN' )
     !
     ! Momentum residual
     !
     if( INOTMASTER ) gevec => remom_nsi

  case ( 'FIXPP' )
     ! 
     ! FIXPP_NSI
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)   
        do ipoin = 1,npoin
           gesca(ipoin) = real(kfl_fixpp_nsi(1,ipoin),rp)
        end do
     end if

  case ( 'BVNAT' )
     ! 
     ! Neumann condition at algebraic level
     !
     if( INOTMASTER ) then
        !call memgen(zero,ndime,npoin)
        !do ipoin = 1,npoin
        !   gevec(1:ndime,ipoin) = solve_sol(1) % block_array(1) % bvnat(1:ndime,ipoin) 
        !end do
        !call PAR_INTERFACE_NODE_EXCHANGE(gevec,'SUM','IN MY ZONE')
        !
        !        call memgen(zero,npoin,0_ip)
        !        print*,associated(solve_sol(1) % block_array(2) % bvnat),size(solve_sol(1) % block_array(2) % bvnat,KIND=ip),npoin
        !        do ipoin = 1,npoin
        !           gesca(ipoin) = solve_sol(1) % block_array(2) % bvnat(1,ipoin) 
        !        end do
        !        call PAR_INTERFACE_NODE_EXCHANGE(gesca,'SUM','IN MY ZONE')
        !
        call memgen(zero,npoin,0_ip)
        solve_sol => momod(modul) % solve(1_ip:)                          !< 2016ABRIL04 
        if(solve_sol(1)%kfl_bvnat /= 1 .and.(memory_size(solve_sol(1)%block_array(2)%bvnat)==npoin) ) &  
             gesca(1:npoin) = solve_sol(1) % block_array(2) % bvnat(1,1:npoin)  
        call PAR_INTERFACE_NODE_EXCHANGE(gesca,'SUM','IN MY ZONE')

     end if

  case ( 'WETNO' )
     if(inotmaster)then
        call memgen(zero,npoin,zero)
        call nsi_wetno(gesca(1:npoin))
     endif

  case ( 'TURMU' )
     !
     ! TURMU: LES turbulent viscosity
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call projec_elements_to_nodes(turmu_nsi,gesca)
     end if

  case ( 'AVMUT' )
     !
     ! AVMUT_NSI: averaged turbulent viscosity
     !
     if(cutim>avtim_nsi) then
        auxi = cutim - avtim_nsi
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin=1,npoin
              gesca(ipoin) = avmut_nsi(ipoin) / auxi
           end do
        end if
     else
        return
     end if

  case ( 'AVSTX' )
     !
     ! AVSTX_NSI: average stress mu_t grad(u)
     !
     if(cutim>avtim_nsi) then
        auxi = cutim - avtim_nsi
        if( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin=1,npoin
              do idime=1,ndime
                 gevec(idime,ipoin) = avstx_nsi(idime,ipoin) / auxi
              end do
           end do
        end if
     else
        return
     end if

  case ( 'AVSTY' )
     !
     ! AVSTY_NSI: average stress mu_t grad(v)
     !
     if(cutim>avtim_nsi) then
        auxi = cutim - avtim_nsi
        if( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin=1,npoin
              do idime=1,ndime
                 gevec(idime,ipoin) = avsty_nsi(idime,ipoin) / auxi
              end do
           end do
        end if
     else
        return
     end if

  case ( 'AVSTZ' )
     !
     ! AVSTZ_NSI: average stress mu_t grad(w)
     !
     if(cutim>avtim_nsi) then
        auxi = cutim - avtim_nsi
        if( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin=1,npoin
              do idime=1,ndime
                 gevec(idime,ipoin) = avstz_nsi(idime,ipoin) / auxi
              end do
           end do
        end if
     else
        return
     end if

  case ( 'BUBBL' )
     !
     ! BUBBLE_NSI: pressure bubble
     !
     if( kfl_bubbl_nsi == 0 ) return
     gesca => bubble_nsi

  case ( 'NORMA' )
     ! 
     ! NORMA: exterior normals
     !  
     call memory_alloca(mem_modul(1:2,modul),'KFL_FIXNO_TMP','nsi_outvar',kfl_fixno_tmp,1_ip ,max(1_ip,npoin))
     call memory_alloca(mem_modul(1:2,modul),'BVESS_TMP',    'nsi_outvar',bvess_tmp    ,1_ip ,max(1_ip,npoin))
     call memory_alloca(mem_modul(1:2,modul),'NORMAL_TMP',   'nsi_outvar',normal_tmp   ,ndime,max(1_ip,npoin))
     !call memgen(zero,ndime,npoin)      
     call memgen(zero,npoin,0_ip)

     solve_sol                => momod(modul) % solve(7:)
     solve_sol(1) % kfl_fixno => kfl_fixno_tmp 
     solve_sol(1) % bvess     => bvess_tmp
     solve_sol(1) % kfl_iffix =  1

     if( INOTMASTER ) then
        call memgen(one,npoin,0_ip)      
        do iboun = 1,nboun
           !if( kfl_fixbo_nsi(iboun) == 5 .or. kfl_fixbo_nsi(iboun) == 5 ) then
           if( kfl_codbo(iboun) == 1 ) then
              gisca(lnodb(1:lnnob(iboun),iboun)) = 1
           end if
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(gisca,'MAX')
        !
        ! Test =>
        !
        rhsid = 0.0_rp
        amatr = 0.0_rp
        !do idime = 1,ndime
        do ipoin = 1,npoin
           ibopo = lpoty(ipoin)
           unkno(ipoin)              = 0.0_rp
           normal_tmp(1:ndime,ipoin) = 0.0_rp
           normal_tmp(ndime,ipoin)   = 1.0_rp           
           rhsid(ipoin)              = 0.0_rp
           if( gisca(ipoin) == 1 .and. ibopo > 0 ) then
              kfl_fixno_tmp(1,ipoin)  =  1
              unkno(ipoin)            =  2.0_rp ! coord(1,ipoin)
              bvess_tmp(1,ipoin)      =  unkno(ipoin)
           end if
        end do
        !call ADR_assemble_laplacian(meshe(ndivi),elmar,amatr)
        call ADR_assemble_extension(1_ip,meshe(ndivi),elmar,normal_tmp,amatr,rhsid)
        !end do
     end if
     call solver_solve(solve_sol,amatr,rhsid,unkno)
     if( INOTMASTER ) then
        gesca(1:npoin) = unkno(1:npoin)
     end if

     goto 10
     !
     ! <= Test
     !

     !
     ! Smooth interior normal
     !
     do idime = 1,ndime
        do ipoin = 1,npoin
           ibopo = lpoty(ipoin)
           unkno(ipoin) = normal_tmp(idime,ipoin)
           rhsid(ipoin) = 0.0_rp
           if( gisca(ipoin) == 1 .and. ibopo > 0 ) then
              kfl_fixno_tmp(1,ipoin)  =  1
              unkno(ipoin)            = -exnor(idime,1,ibopo)
              normal_tmp(idime,ipoin) =  unkno(ipoin)
              bvess_tmp(1,ipoin)      =  unkno(ipoin)
           end if
        end do
        call ADR_assemble_laplacian(meshe(ndivi),elmar,amatr)
        !call ADR_assemble_extension(1_ip,meshe(ndivi),elmar,normal_tmp,amatr,rhsid) 
        call solver_solve(solve_sol,amatr,rhsid,unkno)   
        normal_tmp(idime,1:npoin) = unkno(1:npoin)
     end do
     !do ipoin = 1,npoin
     !   do idime = 1,ndime
     !      gevec(idime,ipoin) = normal_tmp(idime,ipoin)
     !   end do
     !end do

     do ipoin = 1,npoin
        ibopo = lpoty(ipoin)
        unkno(ipoin) = 0.0_rp
        rhsid(ipoin) = 0.0_rp
        if( gisca(ipoin) == 1 .and. ibopo > 0 ) then
           kfl_fixno_tmp(1,ipoin) =  1
           unkno(ipoin)           =  coord(1,ipoin)
           bvess_tmp(1,ipoin)     =  unkno(ipoin)
        end if
     end do
     call ADR_assemble_extension(1_ip,meshe(ndivi),elmar,normal_tmp,amatr,rhsid)
     call solver_solve(solve_sol,amatr,rhsid,unkno)   

     gesca(1:npoin) = unkno(1:npoin)

     if( INOTMASTER ) call memgen(three,npoin,zero)

10   continue
     call memory_deallo(mem_modul(1:2,modul),'KFL_FIXNO_TMP','nsi_outvar',kfl_fixno_tmp)
     call memory_deallo(mem_modul(1:2,modul),'BVESS_TMP',    'nsi_outvar',bvess_tmp    )
     call memory_deallo(mem_modul(1:2,modul),'NORMAL_TMP',   'nsi_outvar',normal_tmp   )

  case ( 'SENSM' )
     !
     ! SENSM: mesh sensitivities
     !
     if( INOTMASTER ) then
        call memgen(zero,ndime,npoin)
        do ipoin = 1,npoin
           do idime = 1,ndime
              gevec(idime,ipoin) = 0.0_rp
           end do
        end do
        do ipoin = 1,npoin
           do idime = 1,ndime
              gevec(idime,ipoin) = sens_mesh(idime,ipoin) 
           end do
        end do
     end if

  case ( 'LAGRA' )
     !
     ! LAGRA_NSI: Lagrange multiplied
     !
     if( kfl_immer_nsi == 1 ) then
        if( INOTMASTER ) gevec => lagra_nsi(:,:,1) 
     else
        return
     end if


  case ( 'ENVEL' )
     !
     ! ENVEL_NSI: ensemble velocity
     !
     if(cutim>avtim_nsi) then
        auxi = entim_nsi - avtim_nsi
        if( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin=1,npoin
              do idime=1,ndime
                 gevec(idime,ipoin) = envel_nsi(idime,ipoin) / auxi
              end do
           end do
        end if
     else
        return
     end if

  case ( 'ENPRE' )
     !
     ! ENPRE_NSI: ensemble pressure
     !
     if(cutim>avtim_nsi) then
        auxi = entim_nsi - avtim_nsi
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin=1,npoin
              gesca(ipoin) = enpre_nsi(ipoin) / auxi
           end do
        end if
     else
        return
     end if

  case ( 'ENVE2' )
     !
     ! ENVE2_NSI: ensemble velocity**2
     !
     if(cutim>avtim_nsi) then
        auxi = entim_nsi - avtim_nsi
        if( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin=1,npoin
              do idime=1,ndime
                 gevec(idime,ipoin) = enve2_nsi(idime,ipoin) / auxi
              end do
           end do
        end if
     else
        return
     end if

  case ( 'ENVXY' )
     !
     ! ENVXY_NSI: ensemble vx * vy
     !
     if(cutim>avtim_nsi) then
        auxi = entim_nsi - avtim_nsi
        if( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin=1,npoin
              do idime=1,ndime
                 gevec(idime,ipoin) = envxy_nsi(idime,ipoin) / auxi
              end do
           end do
        end if
     else
        return
     end if

  case ( 'ENPR2' )
     !
     ! ENPR2_NSI: ensemble pressure**2
     !
     if(cutim>avtim_nsi) then
        auxi = entim_nsi - avtim_nsi
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin=1,npoin
              gesca(ipoin) = enpr2_nsi(ipoin) / auxi 
           end do
        end if
     else
        return
     end if

  case ( 'ENTAN' )
     !
     ! ENTAN_NSI: ensemble tangential force
     !
     if(cutim>avtim_nsi) then
        auxi = entim_nsi - avtim_nsi
        if( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin=1,npoin
              do idime=1,ndime
                 gevec(idime,ipoin) = entan_nsi(idime,ipoin) / auxi
              end do
           end do
        end if
     else
        return
     end if

  case ( 'ENMUT' )
     !
     ! ENMUT_NSI: ensemble turbulent viscosity
     !
     if(cutim>avtim_nsi) then
        auxi = entim_nsi - avtim_nsi
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin=1,npoin
              gesca(ipoin) = enmut_nsi(ipoin) / auxi
           end do
        end if
     else
        return
     end if

  case ( 'ENSTX' )
     !
     ! ENSTX_NSI: ensemble stress mu_t grad(u)
     !
     if(cutim>avtim_nsi) then
        auxi = entim_nsi - avtim_nsi
        if( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin=1,npoin
              do idime=1,ndime
                 gevec(idime,ipoin) = enstx_nsi(idime,ipoin) / auxi
              end do
           end do
        end if
     else
        return
     end if

  case ( 'ENSTY' )
     !
     ! ENSTY_NSI: ensemble stress mu_t grad(v)
     !
     if(cutim>avtim_nsi) then
        auxi = entim_nsi - avtim_nsi
        if( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin=1,npoin
              do idime=1,ndime
                 gevec(idime,ipoin) = ensty_nsi(idime,ipoin) / auxi
              end do
           end do
        end if
     else
        return
     end if

  case ( 'ENSTZ' )
     !
     ! ENSTZ_NSI: ensemble stress mu_t grad(w)
     !
     if(cutim>avtim_nsi) then
        auxi = entim_nsi - avtim_nsi
        if( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin=1,npoin
              do idime=1,ndime
                 gevec(idime,ipoin) = enstz_nsi(idime,ipoin) / auxi
              end do
           end do
        end if
     else
        return
     end if

  case ( 'BTRAC' )
     !
     ! VELAV_KER
     !
     gesca => dt_rho_nsi 


  case ( 'VAFOR' )
     !
     ! VAFOR_NSI : variational forces
     !
     if( INOTMASTER ) then
        gevec => vafor_nsi
     end if

  case ( 'AVVAF' )
     !
     ! AVVAF_NSI: averaged variational forces
     !
     if( INOTMASTER ) then
        gevec => avvaf_nsi
     end if

  case ( 'NOTRA' )
     !
     ! NOTRA : VARIATIONAL Tangential traction
     !
     if( INOTMASTER ) then
        gevec => notra_nsi
     end if

  case ( 'AVNTR' )
     !
     ! AVNTR : Average VARIATIONAL Tangential traction
     !
     if( INOTMASTER ) then
        gevec => avntr_nsi
     end if

  case ( 'AVGTR' )
     !
     ! AVGTR : Average gradient based Tangential traction
     !
     if( INOTMASTER ) then
        gevec => avgtr_nsi
     end if

  case ( 'PORFO' )
     !
     ! PORFO : Porous force
     !
     if( INOTMASTER ) then
        gevec => bupor_nsi
     end if

  case ( 'MASRH' )
     !
     ! drho/dt of Low-Mach restarts without lumping
     !
     if( INOTMASTER ) then
        gesca => mass_rho_nsi(:,1) 
     endif

  case ( 'MOMSK' )
     !
     ! MOMSK : Momentum source from Partis
     !
     call memgen(zero,ndime,max(1_ip,npoin))
     if (associated(momentum_sink)) then
        do ipoin=1,npoin
           gevec(1:ndime,ipoin) =  momentum_sink(1:ndime,ipoin)
        enddo
        call solver_lumped_mass_system(ndime,gevec,EXCHANGE=.false.)
     else
        do ipoin=1,npoin
           gevec(1:ndime,ipoin) = 0.0_rp 
        end do
     endif

  case ( 'DRHOD' )
     !
     ! drho/dt of Low-Mach
     !
     if( INOTEMPTY ) then
        call memgen(zero,npoin,zero)
        do ipoin=1,npoin
           gesca(ipoin) = drhodt_nsi(ipoin)
        end do
        call solver_lumped_mass_system(1_ip,gesca,EXCHANGE=.false.)  
     end if

  case ( 'VEPAV' )
     !
     ! VEPAV: VELAV_KER on nodes
     !
     call memgen(0_ip,ndime,npoin)
     call projec_boundaries_to_nodes(velav_ker,meshe(ndivi),gevec)

  case ( 'QCRIT' )
     !
     ! QCRIT
     !
     if( INOTMASTER ) then
        call memgen(zero,ndime+1_ip,npoin)
        call memgen(zero,npoin,0_ip)
        vorti => gevec
        call vortic(1_ip)
        do ipoin = 1,npoin
           gesca(ipoin) = gevec(ndime+1_ip,ipoin)
        end do
     end if
     nullify(vorti)

  case ( 'AVMOS' )
     !
     ! AVMOS_NSI: averaged momentum source from spray
     !
     if(cutim>avtim_nsi) then
        auxi = cutim - avtim_nsi
        if( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin=1,npoin
              do idime=1,ndime
                 gevec(idime,ipoin) = avmos_nsi(idime,ipoin) / auxi
                 avmos_nsi(idime,ipoin) = 0.0_rp
              end do
           end do
        end if
     else
        return
     end if

  case ( 'EIGEN' )
     !
     ! EIGENVALUE_DT
     !
     if( .not. NSI_FRACTIONAL_STEP ) return
     gesca => tau_nsi
     call nsi_eigen_time_step_all(dtp=tau_nsi)

  case ( 'AVMFL' )
     !
     ! AV_MASS_FLUX_NSI: average mass flux
     !
     if(cutim>avtim_nsi) then
        auxi = cutim - avtim_nsi
        if( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin=1,npoin
              do idime=1,ndime
                 gevec(idime,ipoin) = av_mass_flux_nsi(idime,ipoin) / auxi
                 av_mass_flux_nsi(idime,ipoin) = 0.0_rp
              end do
           end do
        end if
     else
        return
     end if

  case ( 'AVRUU' )
     !
     ! AV_MOM_FLUX_DIAG_NSI: average momentum flux diagonal terms
     !
     if(cutim>avtim_nsi) then
        auxi = cutim - avtim_nsi
        if( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin=1,npoin
              do idime=1,ndime
                 gevec(idime,ipoin) = av_mom_flux_diag_nsi(idime,ipoin) / auxi
                 av_mom_flux_diag_nsi(idime,ipoin) = 0.0_rp
              end do
           end do
        end if
     else
        return
     end if

  case ( 'AVRUV' )
     !
     ! AV_MOM_FLUX_OFF_NSI: average momentum flux off-diagonal terms
     !
     if(cutim>avtim_nsi) then
        auxi = cutim - avtim_nsi
        if( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin=1,npoin
              do idime=1,ndime
                 gevec(idime,ipoin) = av_mom_flux_off_nsi(idime,ipoin) / auxi
                 av_mom_flux_off_nsi(idime,ipoin) = 0.0_rp
              end do
           end do
        end if
     else
        return
     end if

  case ( 'DISTO' )
     !
     ! DISTO
     !
     if( INOTMASTER ) then
        call memgen(zero,ndime+1_ip,npoin)
        call memgen(zero,npoin,zero)
        vorti => gevec
        call vortic(5_ip)
        do ipoin = 1,npoin
           gesca(ipoin) = vorti(1_ip,ipoin)
        end do
     end if
     nullify(vorti)
     
  case ( 'VMAXP' )
     !
     ! VMAXP
     !
     if(cutim>avtim_nsi) then
        if( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin=1,npoin
              do idime=1,ndime
                 gevec(idime,ipoin) = vmaxp_nsi(idime,ipoin)
              end do
           end do
        end if
     endif
     
  case ( 'VMINP' )
     !
     ! VMINP
     !
     if(cutim>avtim_nsi) then
        if( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin=1,npoin
              do idime=1,ndime
                 gevec(idime,ipoin) = vminp_nsi(idime,ipoin)
              end do
           end do
        end if
     endif
     
  case ( 'VAVEP' )
     !
     ! VAVEP
     !
     if(cutim>avtim_nsi) then
        auxi = cutim - avtim_nsi
        if( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin=1,npoin
              do idime=1,ndime
                 gevec(idime,ipoin) = vavep_nsi(idime,ipoin) / auxi
              end do
           end do
        end if
     endif
     
  case ( 'PINDE' )
     !
     ! PINDE
     !
     if(cutim>avtim_nsi) then
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin=1,npoin
              do idime=1,ndime
                 gesca(ipoin) = pinde_nsi(1,ipoin)
              end do
           end do
        end if
     endif

  case ( 'SBVNA' )
     !
     ! SBVNA: algebraic Neumann (natural) boundary conditions
     !
     if( INOTMASTER ) then
        call memgen(zero,ndime,npoin)
        gevec = solve_sol(1) % bvnat
        call PAR_INTERFACE_NODE_EXCHANGE(gevec,'SUM','IN MY ZONE')
     end if

  case ( 'FSIFO' )
     !
     ! FSIFO: FSI body force received from solid in IB coupling for deformable bodies
     !
     if( INOTMASTER ) then
        call memgen(zero,ndime,npoin)
        gevec = fsifo_nsi 
        call PAR_INTERFACE_NODE_EXCHANGE(gevec,'SUM','IN MY ZONE')
     end if

  case ( 'DISTA' )
     !
     ! DISTA:
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           gesca(ipoin) = spare_meshes(1) % dista(ipoin)
        end do
     end if

  case ( 'VELDI' )
     !
     ! VELdi: veloc-velom
     !     
     if( associated(velom) ) then
        call memgen(zero,ndime,npoin)
        do ipoin = 1,npoin
           gevec(1:ndime,ipoin) = veloc(1:ndime,ipoin,1)-velom(1:ndime,ipoin)
        end do
     else
        gevec => veloc(:,:,1)
     end if
     
  case default

     return

  end select
  !
  ! Postprocess
  !
  call outvar(jvari,ittim,cutim,postp(1) % wopos(:,ivari),MESH_ID=imesh)

end subroutine nsi_outvar

subroutine nsi_reset_averages(jvari)
  
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_nastin
  use mod_memchk
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_size
  use mod_arrays,         only : arrays
  use mod_outvar,         only : outvar
  use mod_arrays,         only : arrays_number
  use mod_arrays,         only : arrays_name
  
  implicit none
  
  integer(ip), intent(in) :: jvari   !< 1 for veloc , 2 press , etc...
  integer(ip)             :: ivari
  integer(ip)             :: idime,ipoin

  !
  ! Define postprocess variable
  !
  
  ivari = abs(jvari)

  select case ( arrays_name(ivari) )  

  case ( 'AVVEL' )
     !
     ! AVVEL_NSI: averaged velocity
     !
     if(cutim>avtim_nsi) then
        if( INOTMASTER ) then
           do ipoin=1,npoin
              do idime=1,ndime
                 avvel_nsi(idime,ipoin) = 0.0_rp
              end do
           end do
        end if
     else
        return
     end if

  case ( 'AVPRE' )
     !
     ! AVPRE_NSI: averaged pressure
     !
     if(cutim>avtim_nsi) then
        if( INOTMASTER ) then
           do ipoin=1,npoin
              avpre_nsi(ipoin) = 0.0_rp
           end do
        end if
     else
        return
     end if

  case ( 'AVTAN' )
     !
     ! AVTAN_NSI: averaged tangential force
     !
     if(cutim>avtim_nsi) then
        if( INOTMASTER ) then
           do ipoin=1,npoin
              do idime=1,ndime
                 avtan_nsi(idime,ipoin) = 0.0_rp
              end do
           end do
        end if
     else
        return
     end if

  case ( 'AVVE2' )
     !
     ! AVVE2_NSI: averaged velocity**2
     !
     if(cutim>avtim_nsi) then
        if( INOTMASTER ) then
           do ipoin=1,npoin
              do idime=1,ndime
                 avve2_nsi(idime,ipoin) = 0.0_rp
              end do
           end do
        end if
     else
        return
     end if

  case ( 'AVVXY' )
     !
     ! AVVXY_NSI: averaged vx * vy
     !
     if(cutim>avtim_nsi) then
        if( INOTMASTER ) then
           do ipoin=1,npoin
              do idime=1,ndime
                 avvxy_nsi(idime,ipoin) = 0.0_rp
              end do
           end do
        end if
     else
        return
     end if

  case ( 'AVPR2' )
     !
     ! AVPR2_NSI: averaged pressure**2
     !
     if(cutim>avtim_nsi) then
        if( INOTMASTER ) then
           do ipoin=1,npoin
              avpr2_nsi(ipoin) = 0.0_rp
           end do
        end if
     else
        return
     end if

  case ( 'AVMUT' )
     !
     ! AVMUT_NSI: averaged turbulent viscosity
     !
     if(cutim>avtim_nsi) then
        if( INOTMASTER ) then
           do ipoin=1,npoin
              avmut_nsi(ipoin) = 0.0_rp
           end do
        end if
     else
        return
     end if

  case ( 'AVSTX' )
     !
     ! AVSTX_NSI: average stress mu_t grad(u)
     !
     if(cutim>avtim_nsi) then
        if( INOTMASTER ) then
           do ipoin=1,npoin
              do idime=1,ndime
                 avstx_nsi(idime,ipoin) = 0.0_rp
              end do
           end do
        end if
     else
        return
     end if

  case ( 'AVSTY' )
     !
     ! AVSTY_NSI: average stress mu_t grad(v)
     !
     if(cutim>avtim_nsi) then
        if( INOTMASTER ) then
           do ipoin=1,npoin
              do idime=1,ndime
                 avsty_nsi(idime,ipoin) = 0.0_rp
              end do
           end do
        end if
     else
        return
     end if

  case ( 'AVSTZ' )
     !
     ! AVSTZ_NSI: average stress mu_t grad(w)
     !
     if(cutim>avtim_nsi) then
        if( INOTMASTER ) then
           do ipoin=1,npoin
              do idime=1,ndime
                 avstz_nsi(idime,ipoin) = 0.0_rp
              end do
           end do
        end if
     else
        return
     end if

  case ( 'AVMOS' )
     !
     ! AVMOS_NSI: averaged momentum source from spray
     !
     if(cutim>avtim_nsi) then
        if( INOTMASTER ) then
           do ipoin=1,npoin
              do idime=1,ndime
                 avmos_nsi(idime,ipoin) = 0.0_rp
              end do
           end do
        end if
     else
        return
     end if

  case ( 'AVMFL' )
     !
     ! AV_MASS_FLUX_NSI: average mass flux
     !
     if(cutim>avtim_nsi) then
        if( INOTMASTER ) then
           do ipoin=1,npoin
              do idime=1,ndime
                 av_mass_flux_nsi(idime,ipoin) = 0.0_rp
              end do
           end do
        end if
     else
        return
     end if
    
  case ( 'AVRUU' )
     !
     ! AV_MOM_FLUX_DIAG_NSI: average moemntum flux, diagonal terms
     !
     if(cutim>avtim_nsi) then
        if( INOTMASTER ) then
           do ipoin=1,npoin
              do idime=1,ndime
                 av_mom_flux_diag_nsi(idime,ipoin) = 0.0_rp
              end do
           end do
        end if
     else
        return
     end if
    
  case ( 'AVRUV' )
     !
     ! AV_MOM_FLUX_OFF_NSI: average moemntum flux, off-diagonal terms
     !
     if(cutim>avtim_nsi) then
        if( INOTMASTER ) then
           do ipoin=1,npoin
              do idime=1,ndime
                 av_mom_flux_off_nsi(idime,ipoin) = 0.0_rp
              end do
           end do
        end if
     else
        return
     end if
     
  case ( 'VMAXP' )
     !
     ! VMAXP
     !
     if(cutim>avtim_nsi .and. INOTMASTER ) vmaxp_nsi = 0.0_rp
     
  case ( 'VMINP' )
     !
     ! VMINP
     !
     if(cutim>avtim_nsi .and. INOTMASTER ) vminp_nsi = huge(1.0_rp)
     
  case ( 'VAVEP' )
     !
     ! VAVEP
     !
     if(cutim>avtim_nsi .and. INOTMASTER ) vavep_nsi = 0.0_rp
     
  case ( 'PINDE' )
     !
     ! PINDE
     !
     if(cutim>avtim_nsi .and. INOTMASTER ) pinde_nsi = 0.0_rp
    
  case default

     return

  end select

end subroutine nsi_reset_averages
