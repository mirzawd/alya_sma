!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_iniunk.f90
!> @author  Guillaume Houzeaux
!> @brief   Initial value of velocity and pressurfe
!> @details Set up the initial condition for velocity and pressure.
!>          If this is a restart, initial condition are read from files.
!>          The different options for the velocity are:
!>          \verbatim
!>          Velocity
!>          --------
!>          KFL_INITI_NSI = 0 ... Values from CODES fields
!>                        = 1 ... Constant value
!>                        = 2 ... Constant value in inertial frame of reference
!>                        = 3 ... Stokes flow (Solve a PDE)
!>                        = 4 ... Potential flow (Solve a PDE)
!>                        = 5 ... Laplacian (Solve a PDE)
!>                        = 6 ... Bulk Poiseuille distribution
!>                        = 7 ... Boundary Poiseuille distribution
!>                        < 0 ... From a value function 
!>
!>          Pressure
!>          --------
!>          KFL_INIPR_NSI = 0 ... Hydrostatic pressure if coupling with LEVELS
!>                        = 1 ... Solve Laplace equation (Solve a PDE)
!>                        < 0 ... From a value function 
!>
!>          Other
!>          -----
!>          KFL_INICO_NSI = 1 ... After imposing initial condition,
!>                                solve for one Richardson iteration
!>                                preconditioned by a coarse problem
!>          \endverbatim
!>          Values are stored in position:
!>          \verbatim
!>          VELOC(1:NDIME,1:NPOIN,NPREV_NSI) 
!>          PRESS(1:NPOIN,NPREV_NSI) 
!>          \endverbatim
!> @} 
!------------------------------------------------------------------------
subroutine nsi_iniunk()
  use def_parame 
  use def_master
  use def_elmtyp
  use def_kermod
  use def_domain
  use def_nastin 
  use def_coupli,           only : kfl_efect
  use mod_chktyp,           only : check_type
  use mod_communications,   only : PAR_MAX, PAR_SUM,PAR_INTERFACE_NODE_EXCHANGE
  use mod_array_operations, only : array_operations_copy
  use mod_nsi_hydrostatic,  only : nsi_hydrostatic_pressure
!  use def_kermod,           only : kfl_noslw_ker,ndivi,kfl_fixbo_nsw_ker
  use mod_projec,           only : calc_massb
  use mod_memory,           only : memory_alloca,memory_deallo
  use mod_ker_tendencies,   only : ker_tendencies_initialize_u, kfl_tendencies_ker
  use mod_nsi_efect,        only : nsi_set_fringe_node_bcs

  use mod_ker_space_time_function

  use mod_vortex,           only : VORTEX_ID, vortex_iniunk 

  implicit none
  integer(ip)              :: ipoin,idime,inodb,iboun,ibopo,izone
  integer(ip)              :: kfl_value,ifunc
  real(rp)                 :: Rdist,dummr(3),rot_matrix(ndime,ndime)
  real(rp)                 :: vefun(3)
  real(rp)                 :: aux,avgu(ndime), avgu2(ndime)
  integer(ip)              :: totalp 
  integer(4)               :: size_seed
  integer(4),  allocatable :: seed(:)
  integer(ip)              :: aux_size_seed
  integer(ip), allocatable :: aux_seed(:)
  logical(lg), pointer     :: boundary_mask(:)

  
  izone = lzone(ID_NASTIN)

  if( kfl_rstar == 0 ) then  
     !
     ! Info and latex files
     !
     call nsi_outinf(1_ip)

     !-------------------------------------------------------------------
     !
     ! Non-PDE-based initial solution in global system
     !
     !-------------------------------------------------------------------

     if( INOTMASTER ) then

        if( kfl_exacs_nsi /= 0 ) then
           !
           ! Exact solution
           !
           call nsi_exaerr(5_ip)

        else if( kfl_initi_nsi < 0 ) then
           !
           ! Take initial condition from a value function field
           ! 
           kfl_value = -kfl_initi_nsi
           call check_type(xfiel,kfl_value,ndime,npoin) ! Check if value function exist                      
           do ipoin = 1,npoin
              veloc(1:ndime,ipoin,1) = xfiel(kfl_value) % a(1:ndime,ipoin,1)
           end do

        else if( kfl_initi_nsi == 0 ) then
         
           if (kfl_tendencies_ker) then
              !
              !  Initial conditions from WRF in tendencies file (when ini cond are not prescribed)
              !
              PRInt *, 'INICIALIZANDO'
              print *, 'inicializando nastin con tendencias'
              call ker_tendencies_initialize_u(veloc, walld)
           else
              !
              ! Take initial conditions as prescribed in codes fields
              !
              do ipoin = 1,npoin
                 veloc(1:ndime,ipoin,1) = bvess_nsi(1:ndime,ipoin,1)
              end do
           end if
        else if( kfl_initi_nsi == 1 ) then
           !
           ! Constant value prescribed
           !
           do ipoin = 1,npoin
              veloc(1:ndime,ipoin,1) = velin_nsi(1:ndime)
           end do
          
        else if( kfl_initi_nsi == 2 ) then
           !
           ! Constant value prescribed in inertial frame of reference
           !
           call nsi_updfor()
           dummr = velin_nsi
           call rotmat(ndime,-cutim*fvnoa_nsi,fvdia_nsi,rot_matrix)
           call mbvab0(velin_nsi,rot_matrix,dummr,ndime,ndime)
           if( ndime == 2 ) then
              do ipoin = 1,npoin
                 veloc(1,ipoin,1) = velin_nsi(1) + fvela_nsi(3) * coord(2,ipoin)
                 veloc(2,ipoin,1) = velin_nsi(2) - fvela_nsi(3) * coord(1,ipoin)
              end do
           else
              do ipoin = 1,npoin
                 veloc(1,ipoin,1) = velin_nsi(1) + fvela_nsi(3) * coord(2,ipoin) &
                      &                                  - fvela_nsi(2) * coord(3,ipoin) 
                 veloc(2,ipoin,1) = velin_nsi(2) - fvela_nsi(3) * coord(1,ipoin) &
                      &                                  + fvela_nsi(1) * coord(3,ipoin)
                 veloc(3,ipoin,1) = velin_nsi(3) - fvela_nsi(1) * coord(2,ipoin) &
                      &                                  + fvela_nsi(2) * coord(1,ipoin)       
              end do
           end if
           velin_nsi = dummr  
           
        else if(kfl_initi_nsi==6) then
           !
           ! Bulk Poiseuille distribution
           !
           do ipoin = 1,npoin
              Rdist = (coord(1,ipoin)-poise_nsi(4))**2 + (coord(2,ipoin)-poise_nsi(5))**2 
              if (ndime == 3) Rdist = Rdist + (coord(3,ipoin)-poise_nsi(6))**2
              do idime = 1,ndime
                 if (int(poise_nsi(1),ip)==idime) then  ! This dim is the axis of flow
                    Rdist = Rdist - (coord(idime,ipoin)-poise_nsi(3+idime))**2  ! Only use distance to the axis of flow
                    Rdist = min( Rdist, poise_nsi(2)**2 ) ! Zero out veloxity outside the radius of influence
                    veloc(idime,ipoin,1) = poise_nsi(3)*(1.0_rp - Rdist/poise_nsi(2)**2 )
                 else
                    veloc(idime,ipoin,1) = 0.0_rp
                 endif
                 if(kfl_fixno_nsi(idime,ipoin) == 16 ) then ! We modify the inlet boundary condition to the same distribution
                    bvess_nsi(idime,ipoin,1) = veloc(idime,ipoin,1)
                    kfl_fixno_nsi(idime,ipoin) = 1 ! And we make it a Dirichlet condition from now on
                 end if
              end do              
           enddo
           
        else if(kfl_initi_nsi==7) then
           !
           ! Boundary Poiseuille distribution, Bulk is constant, hopefully loaded in velin_nsi
           !
           do ipoin = 1,npoin
              do idime = 1,ndime
                 if(kfl_fixno_nsi(idime,ipoin) == 16 ) then ! We modify the inlet boundary condition to the Poiseuille distribution
                    Rdist = (coord(1,ipoin)-poise_nsi(4))**2 + (coord(2,ipoin)-poise_nsi(5))**2 
                    if (ndime == 3) Rdist = Rdist + (coord(3,ipoin)-poise_nsi(6))**2
                    if (int(poise_nsi(1),ip)==idime) then  ! This dim is the axis of flow
                       Rdist = Rdist - (coord(idime,ipoin)-poise_nsi(3+idime))**2  ! Only use distance to the axis of flow
                       Rdist = min( Rdist, poise_nsi(2)**2 ) ! Zero out veloxity outside the radius of influence
                       veloc(idime,ipoin,1) = poise_nsi(3)*(1.0_rp - Rdist/poise_nsi(2)**2 )
                    else
                       veloc(idime,ipoin,1) = 0.0_rp
                    endif
                    bvess_nsi(idime,ipoin,1) = veloc(idime,ipoin,1)
                    kfl_fixno_nsi(idime,ipoin) = 1 ! And we make it a Dirichlet condition from now on                
                 else
                    veloc(idime,ipoin,1) = velin_nsi(idime)
                 end if
              end do
           end do


        else if(kfl_initi_nsi==VORTEX_ID) then
           call vortex_iniunk() 

        else if(kfl_initi_nsi>100) then
           !
           ! Space time function
           !
           ifunc = kfl_initi_nsi - 100
           do ipoin = 1,npoin
              call ker_space_time_function(&
                   ifunc,coord(1,ipoin),coord(2,ipoin),coord(ndime,ipoin),cutim,vefun(1:ndime))
              do idime = 1,ndime
                 veloc(idime,ipoin,1) = vefun(idime) 
              end do
           end do     
           
        end if
        !
        ! Functions and specific boundary conditions
        !
     end if
     !
     ! Activate flag for Low Mach in combustion (properties on GP or Node)
     !
     kfl_lookg_nsi=0
     if (kfl_regim_nsi==3 .and. associated(tempe_gp)) kfl_lookg_nsi = 1
     
     if( INOTEMPTY ) then
        !
        ! Bemol
        !
        if( bemol_nsi > 0.0_rp ) then
           if( kfl_local_nsi == 0 ) kfl_local_nsi = -1
           do iboun = 1,nboun
              if(kfl_fixbo_nsi(iboun) == 13 ) then
                 do inodb = 1,nnode(ltypb(iboun))
                    ipoin = lnodb(inodb,iboun)
                    ibopo = lpoty(ipoin)
                    if( ibopo > 0 ) then
                       if( kfl_fixrs_nsi(ipoin) == 0 ) then 
                          kfl_fixrs_nsi(ipoin)   = -4
                          kfl_fixno_nsi(1,ipoin) =  1
                       end if
                    end if
                 end do
              end if
           end do
        end if
        !
        ! Impose boundary/initial condition values
        !
        call nsi_rotunk(3_ip,veloc(1,1,1)) ! Global to local
        do ipoin = 1,npoin
           do idime = 1,ndime
              if(       kfl_fixno_nsi(idime,ipoin) ==  1 &
                   .or. kfl_fixno_nsi(idime,ipoin) ==  8 &
                   .or. kfl_fixno_nsi(idime,ipoin) ==  9 &
                   .or. kfl_fixno_nsi(idime,ipoin) ==  5 &
                   .or. kfl_fixno_nsi(idime,ipoin) ==  6 &
                   .or. kfl_fixno_nsi(idime,ipoin) == 10 ) then
                 veloc(idime,ipoin,1) = bvess_nsi(idime,ipoin,1)
              end if
           end do
        end do
        call nsi_rotunk(4_ip,veloc(1,1,1)) ! Local to global
        !
        ! Pressure initial condition
        !
        if( kfl_exacs_nsi == 0 ) then
           if( kfl_regim_nsi == 0 .or. kfl_regim_nsi == 3 ) then
              do ipoin = 1,npoin
                 press(ipoin,1) = valpr_nsi
              end do
           else if( kfl_regim_nsi == 1  ) then
              do ipoin = 1,npoin
                 press(ipoin,1) = prthe(1)
              end do
           else if( kfl_regim_nsi == 2 ) then
              call runend('NOT CODED')
              !do ipoin = 1,npoin
              !   densi(ipoin,1) = 
              !end do
           end if
        end if
        !
        ! Bemol
        !
        if( bemol_nsi > 0.0_rp ) then
           if( kfl_local_nsi == -1 ) kfl_local_nsi = 0
           do iboun = 1,nboun
              if(kfl_fixbo_nsi(iboun) == 13 ) then
                 do inodb = 1,nnode(ltypb(iboun))
                    ipoin = lnodb(inodb,iboun)
                    ibopo = lpoty(ipoin)
                    if( ibopo > 0 ) then
                       if( kfl_fixrs_nsi(ipoin) == -4 ) then 
                          kfl_fixrs_nsi(ipoin)   = 0
                          kfl_fixno_nsi(1,ipoin) = 0
                       end if
                    end if
                 end do
              end if
           end do
        end if

     end if

     if(kfl_initi_nsi==8) then
        ! Random field here because we need comunications
        ! from Kempf et al. (2005)
        size_seed = 1_4
        call random_seed(size = size_seed)
        allocate(seed(size_seed))
        allocate(aux_seed(size_seed))
        call random_seed(get=seed)
        aux_size_seed = INT(size_seed, ip )
        aux_seed = INT(seed, ip )
        call PAR_MAX(aux_size_seed,aux_seed)

        if( INOTMASTER ) then
           call random_seed(put=seed)
           avgu = 0.0_rp
           totalp = 0_ip
           veloc(:,:,1) = 0_rp 
           do ipoin = 1,npoi1
              do idime = 1,ndime
                 call random_number(aux)
                 veloc(idime,ipoin,1) = aux
                 avgu(idime) = avgu(idime) +  veloc(idime,ipoin,1)
              end do
              totalp = totalp + 1
           end do
           do ipoin = npoi2,npoi3
              do idime = 1,ndime
                 call random_number(aux)
                 veloc(idime,ipoin,1) = aux
                 avgu(idime) = avgu(idime) +  veloc(idime,ipoin,1)
              end do
              totalp = totalp + 1
           end do
        end if
        deallocate (seed)
        call PAR_SUM(ndime,avgu)
        call PAR_SUM(totalp)
        ! fix  <u> = 0
        if( INOTMASTER ) then
           avgu2 = 0.0_rp
           do ipoin = 1,npoi1
              do idime = 1,ndime
                 veloc(idime,ipoin,1) = veloc(idime,ipoin,1) - (avgu(idime)/(real(totalp,rp)*1.0_rp))
                 avgu2(idime) = avgu2(idime) +  veloc(idime,ipoin,1)**2.0_rp
              end do
           end do
           do ipoin = npoi2,npoi3
              do idime = 1,ndime
                 veloc(idime,ipoin,1) = veloc(idime,ipoin,1) - (avgu(idime)/(real(totalp,rp)*1.0_rp))
                 avgu2(idime) = avgu2(idime) +  veloc(idime,ipoin,1)**2.0_rp
              end do
           end do
        end if
        call PAR_SUM(ndime,avgu2)
        !fix  <uu> = 1
        !fix  ui = ui/sqrt(Vol)
        if( INOTMASTER ) then
           do idime = 1,ndime
              avgu2(idime) = avgu2(idime)/(real(totalp,rp)*1.0_rp)
              avgu2(idime) = avgu2(idime)**(0.5_rp)
           end do
           do ipoin = 1,npoi1
              do idime = 1,ndime
                 veloc(idime,ipoin,1) = veloc(idime,ipoin,1) / (avgu2(idime))
                 veloc(idime,ipoin,1) = veloc(idime,ipoin,1) / (vmass(ipoin)**(0.5_rp))
              end do
           end do
           do ipoin = npoi2,npoi3
              do idime = 1,ndime
                 veloc(idime,ipoin,1) = veloc(idime,ipoin,1) / (avgu2(idime))
                 veloc(idime,ipoin,1) = veloc(idime,ipoin,1) / (vmass(ipoin)**(0.5_rp))
              end do
           end do
           call PAR_INTERFACE_NODE_EXCHANGE(ndime,veloc(:,:,1),'SUM','IN MY CODE')
        end if
     end if
     !
     ! Set zero-velocity Dirichlet BCs on fringe nodes for Embedded Finite
     ! Element Coupling Technique
     !
     if( kfl_efect ) call nsi_set_fringe_node_bcs(2_ip)

     !-------------------------------------------------------------------
     !
     ! Apply boundary conditions
     !
     !-------------------------------------------------------------------
     
     call nsi_updbcs(ITASK_INIUNK)

     !-------------------------------------------------------------------
     !
     ! Initial solution: solve coarse grid
     !
     !-------------------------------------------------------------------

     if( kfl_inico_nsi == 1 ) then
        call nsi_inicoa() 
     end if

     !-------------------------------------------------------------------
     !
     ! Initial pressure
     !
     !-------------------------------------------------------------------

     if( kfl_inipr_nsi == 1 ) then
        !
        ! Solve a Laplacian
        !
        call nsi_inipre()

     else if ( kfl_inipr_nsi < 0 .and. INOTMASTER ) then
        !
        ! Take initial condition from a value function field
        ! 
        kfl_value = -kfl_inipr_nsi
        call check_type(xfiel,kfl_value,1_ip,npoin) ! Check if value function exist                      
        press(1:npoin,1) = xfiel(kfl_value) % a(1,1:npoin,1)

     end if

     !-------------------------------------------------------------------
     !
     ! Compute hydrostatic pressure if required
     !
     !-------------------------------------------------------------------

     call nsi_hydrostatic_pressure(ITASK_INIUNK)

     !-------------------------------------------------------------------
     !
     ! First guess based on PDE: Stokes, potential, Laplacian
     !
     !------------------------------------------------------------------- 

     if( kfl_initi_nsi == 3 ) then
        call nsi_stokes()                 ! Initial solution is Stokes flow
     else if( kfl_initi_nsi == 4 ) then
        call nsi_potent()                 ! Initial solution is potential flow
     else if( kfl_initi_nsi == 5 ) then
        call runend('NSI_INIUNK: OBSOLETE OPTION')
     end if

     !-------------------------------------------------------------------
     !
     ! Divergence correction
     !
     !-------------------------------------------------------------------
     
     if( kfl_divcorrec_nsi /= 0 ) call nsi_zero_divergence_velocity()
     
     !-------------------------------------------------------------------
     !
     ! Interpolation from coarse to fine mesh
     !
     !-------------------------------------------------------------------
     
     if(kfl_meshi_nsi /= 0_ip) call nsi_coarfine(1_ip)

     !-------------------------------------------------------------------
     !
     ! (:,3) <= (:,1)
     !
     !-------------------------------------------------------------------
     
     call nsi_updunk( ITASK_INIUNK )

     !-------------------------------------------------------------------
     !
     ! Initial value for time-averaged velocity
     !
     !-------------------------------------------------------------------

     if( (kfl_wlaav_ker == 1_ip) ) then
        call nsi_wallav(3_ip)
        call nsi_wallav(4_ip)  ! beware tell georgios that I have added this that affects his stuff. I belive it is correct as I have it now
     end if

  end if


  !-------------------------------------------------------------------
  !
  ! Initial flags
  !
  !-------------------------------------------------------------------

  call nsi_inivar(2_ip)
  call nsi_coupli(ITASK_INIUNK) 

  !-------------------------------------------------------------------
  !
  ! Read forward values from fields for adjoint
  !
  !-------------------------------------------------------------------
  
  if( INOTMASTER .and. kfl_adj_prob == 1 ) then
    if( nfiel_nsi(1) == 0 .or. nfiel_nsi(2) == 0 ) call runend('NSI_INIUNK: NO FIELDS WERE GIVEN FOR FORWARD VALUES')
    do ipoin=1,npoin
      press_forw(ipoin,:) = xfiel(-nfiel_nsi(1))%a(1,ipoin,1)
      do idime=1,ndime
        veloc_forw(idime,ipoin,:) = xfiel(-nfiel_nsi(2))%a(idime,ipoin,1)
      enddo
      if( kfl_rstar == 0 ) then 
        veloc(:,ipoin,:) = 0.0_rp
        press(ipoin,:)   = 0.0_rp
      endif  
    end do
  end if

  !----------------------------------------------------------------------
  !
  ! these 9 lines where previouly in nsi_turnon - I pass them here because nsi_turnon is called before ker_turnon   -- Strange
  ! Therefore kfl_fixbo_nsw_ker is not allocated  -- The problem arose when I changed kfl_fixbo_nsi to  kfl_fixbo_nsw_ker
  ! NO ES BUENA IDEA QUE ESTEN ACA sino tendre que cambiar massb-nsi a massb-ker --  ver porqeu se llama antes a nsiturnon que a ker_turnon
  ! yo pensarái que debrái ser al reves -- hablar con guillaume
  ! => GOTO BEGRUN and allocate in nsi_solmem
  !----------------------------------------------------------------------

  if( INOTEMPTY .and. kfl_noslw_ker /= 0 ) then
     nullify(boundary_mask)
     call memory_alloca(mem_modul(1:2,modul),'BOUNDARY_MASK','nsi_iniunk',boundary_mask,max(meshe(ndivi)%nboun,1_ip))
     boundary_mask = .false.
     do iboun =1,meshe(ndivi)%nboun
        if ( kfl_fixbo_nsw_ker(iboun) == 1_ip ) boundary_mask(iboun) = .true.
     end do
     call calc_massb(meshe(ndivi),elmar,massb_nsi,boundary_mask)
     call memory_deallo(mem_modul(1:2,modul),'BOUNDARY_MASK','nsi_iniunk',boundary_mask)
  end if

  !-------------------------------------------------------------------
  !
  ! For time & space BC from file: get boundary nodes
  !
  !-------------------------------------------------------------------
  
!  if (kfl_bnods_nsi == 1) then
!     ibopo = 0
!     do ipoin = 1,npoin
!        iboun_nsi(ipoin) = 0
!        do iboun = 1,nbnod_nsi
!           if(lninv_loc(ipoin) == int(bntab_nsi(iboun,1),ip)) then
!              iboun_nsi(ipoin) = iboun
!              ibopo = ibopo + 1
!           end if
!        end do
!     end do
!     if (ibopo /= 0) write(*,*) 'total number of boundary nodes: ', nbnod_nsi, '. computed by processor ', kfl_paral,': ', ibopo
!  end if

  !-------------------------------------------------------------------
  !
  ! Identify which points belong in a wall law boundary for postprocessing
  ! => GOTO BEGRUN and allocate in nsi_solmem
  !
  !-------------------------------------------------------------------

  if (INOTMASTER) then
     do ipoin = 1,npoin
        kfl_wlawf_nsi(ipoin) = 0_ip
     end do
     do iboun = 1, nboun
        if (   kfl_fixbo_nsi(iboun) ==  3  .or. &    ! Wall law
             & kfl_fixbo_nsi(iboun) == 13  .or. &    ! Wall law + open pressure
             & kfl_fixbo_nsi(iboun) == 18) then      ! u.n in weak form
!             & kfl_fixbo_nsi(iboun) == 22 ) then    ! Boundary traction is imposed from an auxiliary RANS simulation (Two-layer wall model)
           do inodb = 1,nnode(ltypb(iboun))
              ipoin = lnodb(inodb,iboun)
              kfl_wlawf_nsi(ipoin) = 1
           end do
        end if
     end do
     call PAR_INTERFACE_NODE_EXCHANGE(kfl_wlawf_nsi,'SUM','IN MY CODE')
  end if

end subroutine nsi_iniunk

