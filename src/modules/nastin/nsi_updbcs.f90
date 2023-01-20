!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup NastinInput
!> @{
!> @file    nsi_updbcs.f90
!> @date    29/01/2018
!> @author  Guillaume Houzeaux
!> @brief   Update boundary conditions
!> @details This routine updates the velocity boundary conditions:
!>          1. Before a time step begins
!>          2. Before a global iteration begins
!>          3. Before an inner iteration begins
!>          4. Force boundary conditions in UNKNO
!> @} 
!-----------------------------------------------------------------------
subroutine nsi_updbcs(itask)

  use def_parame,                  only :  pi
  use def_elmtyp
  use def_master
  use def_kermod
  use def_domain
  use def_nastin
  use mod_ker_space_time_function
  use mod_communications,          only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications,          only : PAR_MIN
  use mod_couplings,               only : COU_PUT_VALUE_ON_TARGET
  use mod_projec,                  only : projec_mass_conservation
  use mod_memory,                  only : memory_alloca
  use mod_memory,                  only : memory_alloca_min
  use mod_memory,                  only : memory_deallo
  use mod_memory,                  only : memory_size
  use mod_output,                  only : output_mesh_gid_format
  use mod_couplings,               only : couplings_impose_dirichlet
  use mod_std
  use mod_ker_subdomain,           only : ker_subdomain_motion_exists
  use mod_windk,                   only : mod_windk_interact
  use mod_pump,                    only : pump_curve_flow
  use mod_ker_functions,           only : ker_functions
  use mod_local_basis,             only : local_basis_global_to_local
  use mod_local_basis,             only : local_basis_local_to_global
  use mod_maths_arrays,            only : maths_findloc
  use mod_messages,                only : messages_live

  implicit none
  integer(ip), intent(in) :: itask
  integer(ip), save       :: ipass = 0
  integer(ip)             :: ipoin,idime,iboun,dummi,itotv,izdom,jpoin
  integer(ip)             :: inodb,inode,ifunc,ibopo,jdime,kpoin,ii
  integer(ip)             :: iflow,itype,itable,ibtim
  integer(ip)             :: i_left,i_right,ind_left,ind_right
  real(rp)                :: venew,veold,vefun(3),xxmin(3)
  real(rp)                :: rot_matrix(ndime,ndime),udotn
  real(rp)                :: bvess_new(3),xnorm,bvnew,bvold
  real(rp)                :: alpha,center_gravity(3),p_in,dummr
  real(rp),    external   :: funcre
  logical(lg), pointer    :: lboun_nsi(:)
  logical(lg), pointer    :: lpoin_nsi(:)
  real(rp),    pointer    :: bvess_loc(:,:)
  
  nullify(lboun_nsi)
  nullify(lpoin_nsi)
  nullify(bvess_loc)

  select case( itask )

  case( ITASK_TURNON )

     !-------------------------------------------------------------------
     !
     ! NODPR_NSI: At the beginning of the run
     !
     !-------------------------------------------------------------------

     nodpr_nsi = 0

     if( kfl_confi_nsi == 0 ) then
        !
        ! Prescribe pressure: select automatically the node with lowest global numbering
        !    
        inode     = huge(ip)
        nodpr_nsi = 0 
        do ipoin = 1,npoin
           if( lpoty(ipoin) /= 0 .and. lmast(ipoin) == 0 ) then
              if( lninv_loc(ipoin) < inode ) inode = lninv_loc(ipoin)
           end if
        end do
        call PAR_MIN(inode)
        do ipoin = 1,npoin
           if( lninv_loc(ipoin) == inode ) then
              nodpr_nsi = ipoin
           end if
        end do

     else if( kfl_confi_nsi > 0 ) then 
        !
        ! Prescribe pressure: On node
        ! If MM is used, the numbering refers to the MM(0) level, that is the original one
        !     
        if( INOTMASTER ) then 
           inode     = abs(nodpr_global_nsi)
           nodpr_nsi = 0
           kpoin     = 0
           do while( kpoin < npoin )
              kpoin = kpoin + 1
              ipoin = kpoin
              if( lninv_loc(ipoin) == inode ) then
                 nodpr_nsi = ipoin
                 kpoin = npoin
                 if( lmast(ipoin) /= 0 ) call runend('NSI_UPDBCS: PRESSURE CANNOT BE PRESCRIBED ON A PERDIODIC NODE')
              end if
           end do
        end if

     end if

     !-------------------------------------------------------------------
     !
     ! KFL_FIXPR_NSI
     !
     !-------------------------------------------------------------------

     if( kfl_confi_nsi >= 0 .and. nodpr_nsi > 0 .and. INOTMASTER ) then
        !
        ! Confied flow
        !
        !if( nodpr_global_nsi > 0 ) then           
        !   ibopo = lpoty(nodpr_nsi)
        !   kfl_fixpr_nsi(1,nodpr_nsi) = 1
        !else
        ibopo = lpoty(nodpr_nsi)
        kfl_fixpr_nsi(1,nodpr_nsi) = 1           
        !end if

     else if( kfl_confi_nsi == -1 ) then 
        !
        ! Automatic bc for pressure: Prescribe when velocity is free
        !     
        if( INOTMASTER  ) then
           do ipoin = 1,npoin
              ibopo = lpoty(ipoin)
              if( ibopo > 0 ) then
                 if( kfl_fixpr_nsi(1,ipoin) == 0 ) then
                    dummi = 0
                    do idime = 1,ndime
                       if( kfl_fixno_nsi(idime,ipoin) == 0 ) dummi = dummi + 1
                    end do
                    if( dummi == ndime ) kfl_fixpr_nsi(1,ipoin)=1  
                    !! You should also check above that the node is not periodic
                 end if
              end if
           end do
        end if
     end if
     !
     ! Free pressure on internal boundary nodes
     ! Eliminate internal boundary nodes
     !
     call memgen(1_ip,npoin,0_ip)
     do iboun = 1,nboun
        if( lboch(iboun) == BOINT ) then
           do inodb = 1,nnode(ltypb(iboun))
              ipoin = lnodb(inodb,iboun)
              if( gisca(ipoin) /= 2 ) gisca(ipoin) = 1
           end do
        else
           do inodb = 1,nnode(ltypb(iboun))
              ipoin = lnodb(inodb,iboun)
              gisca(ipoin) = 2
           end do
        end if
     end do
     call PAR_INTERFACE_NODE_EXCHANGE(gisca,'MAX')
     do ipoin = 1,npoin
        if( gisca(ipoin) == 1 ) then
           kfl_fixpr_nsi(1,ipoin) = 0
        end if
     end do
     call memgen(3_ip,npoin,0_ip)     
     ! 
     ! Bemol: modify if using Parall
     !
     if( bemol_nsi > 0.0_rp .and. INOTMASTER ) then
        do iboun = 1,nboun
           if( kfl_fixbo_nsi(iboun) == 13 ) then
              do inodb = 1,nnode(ltypb(iboun))
                 ipoin = lnodb(inodb,iboun)
                 kfl_fixpr_nsi(1,ipoin) = 1
              end do
           end if
        end do
     end if
     !
     ! No-penetration in weak form
     !
     if( INOTEMPTY ) then
        call memgen(1_ip,npoin,0_ip)
        do iboun = 1,nboun
           if( kfl_fixbo_nsi(iboun) == 18 ) then
              do inodb = 1,nnode(abs(ltypb(iboun)))
                 ipoin = lnodb(inodb,iboun)
                 gisca(ipoin) = 1
              end do
           end if
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(gisca,'SUM')
        do ipoin = 1,npoin
           if( gisca(ipoin) /= 0 ) then
              kfl_fixpr_nsi(1,ipoin) = 0
           end if
        end do
        call memgen(3_ip,npoin,0_ip)
     end if

     if( INOTEMPTY .and. exfpr_nsi == 1 ) then
        write (*,*) 'CONDITION EXFPR, EXTENDING FIXPR NODES '
        call memgen(1_ip,npoin,0_ip)
        do ipoin = 1,npoin
           if( kfl_fixpr_nsi(1,ipoin) == 1 ) then 
              do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                 jpoin = c_dom(izdom)
                 if( kfl_fixpr_nsi(1,jpoin) == 0 ) gisca(jpoin) = 1
              end do
           end if
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(gisca,'SUM')
        do ipoin = 1,npoin
           if( gisca(ipoin) /= 0 .and. lpoty(ipoin)/=0) then
              kfl_fixpr_nsi(1,ipoin) = 1
           end if
        end do
        call memgen(3_ip,npoin,0_ip)
     end if

     !-------------------------------------------------------------------
     !
     ! BPESS_NSI: Pressure Dirichlet
     !
     !-------------------------------------------------------------------

     if( nodpr_nsi > 0 .and. INOTMASTER ) then
        bpess_nsi(1,nodpr_nsi,1)   = valpr_nsi
        kfl_fixpp_nsi(1,nodpr_nsi) = 1
     end if
     !
     ! Neumann BC in Fractional-Step scheme
     !
     call nsi_bougra()

     !-------------------------------------------------------------------
     !
     ! Other stuffs
     !
     !-------------------------------------------------------------------
     !
     ! Modify manually the boundary conditions
     !
     call nsi_modbcs()
     !
     ! Check outflows with too few nodes
     !
     call nsi_detect_outflows()

  case( ITASK_BEGSTE , ITASK_INIUNK )

     !-------------------------------------------------------------------
     !  
     ! Before a time step
     !     
     !-------------------------------------------------------------------
     !  
     ! Fix pressure on lower-left node
     !          
     if( kfl_confi_nsi == -3 ) then
        
        inode     = huge(ip)
        xxmin     = huge(1.0_rp)
        nodpr_nsi = 0

        call memory_alloca(mem_modul(1:2,modul),'LPOIN_NSI','nsi_updbcs',lpoin_nsi,npoin)
        do ipoin = 1,npoin
           if( lpoty(ipoin) /= 0 .and. lmast(ipoin) == 0 ) then
              if( coord(1,ipoin) <= xxmin(1) ) then
                 lpoin_nsi(ipoin) = .true.
                 xxmin(1)         = coord(1,ipoin)
              end if
           end if
        end do
        do idime = 2,ndime
           do ipoin = 1,npoin
              if( lpoin_nsi(ipoin) ) then
                 lpoin_nsi(ipoin) = .false.
                 if( abs(coord(idime-1,ipoin) - xxmin(idime-1)) <= zeror ) then
                    if( coord(idime,ipoin) <= xxmin(idime) ) then
                       xxmin(idime)     = coord(idime,ipoin)
                       lpoin_nsi(ipoin) = .true.
                       inode            = ipoin
                    end if
                 end if
              end if
           end do
        end do        
        call PAR_MIN(ndime,xxmin)
        do ipoin = 1,npoin
           if( lpoin_nsi(ipoin) ) then
              idime = 0
              do while( idime < ndime )
                 idime = idime + 1
                 if( abs(coord(idime,ipoin)-xxmin(idime)) > 1.0e-16_rp ) then
                    idime = ndime+1
                 else if( idime == ndime ) then
                    nodpr_nsi = ipoin
                 end if
              end do
           end if
        end do
        call memory_deallo(mem_modul(1:2,modul),'LPOIN_NSI','nsi_updbcs',lpoin_nsi)

        do ipoin = 1,npoin
           kfl_fixpr_nsi(1,ipoin) = 0
        end do
        if( nodpr_nsi > 0 .and. INOTMASTER ) then
           bpess_nsi(1,nodpr_nsi,1)   = valpr_nsi
           kfl_fixpp_nsi(1,nodpr_nsi) = 1
           kfl_fixpr_nsi(1,nodpr_nsi) = 1
        end if
        
     end if
     
     if( INOTEMPTY .and. kfl_conbc_nsi == 0 ) then
        !
        ! Space time function, fields and coded functions
        !
        do ipoin = 1,npoin

           ifunc = kfl_funno_nsi(ipoin)
           itype = kfl_funtn_nsi(ipoin)

           call ker_functions(ipoin,ifunc,itype,bvess_nsi(:,ipoin,2),vefun)
           !
           ! Velocity update according to time discretization
           !
           if( ifunc /= 0 ) then
              if( kfl_timei_nsi /= 0 .and. kfl_tiacc_nsi == 2 .and. kfl_tisch_nsi == 1 ) then
                 do idime = 1,ndime
                    veold = veloc(idime,ipoin,ncomp_nsi)
                    bvess_nsi(idime,ipoin,1) = 0.50_rp*(vefun(idime)+veold)
                 end do
              else
                 do idime = 1,ndime
                    bvess_nsi(idime,ipoin,1) = vefun(idime)
                 end do
              end if
           end if
        end do
        !
        ! Nodes with fixity 7 (open or close uses) bvess_nsi which has already been updated if it comes from a transient field
        !
        if( kfl_exist_fixi7_nsi==1 .and. itask==ITASK_BEGSTE ) then    
           if ( ipass==0) then
              ipass = 1
              do ipoin = 1,npoin
                 if( abs(kfl_fixno_nsi(1,ipoin)) == 7.or. abs(kfl_fixno_nsi(2,ipoin)) == 7.or. &
                      abs(kfl_fixno_nsi(ndime,ipoin)) == 7 ) then
                    if ( kfl_fixrs_nsi(ipoin) /= 0_ip ) then
                       print*,'ipoin with skew',lninv_loc(ipoin)
                       call runend('nsi_updbcs:ipoin with skew')
                    end if
                 end if
              end do
           end if

           call memgen(1_ip,npoin,0_ip)   !allocate gisca
           call open_close(1_ip,kfl_fixno_nsi,bvess_nsi,ndime)  !returns gisca(ipoin) >0 inflow or =< 0 
           do ipoin = 1,npoin
              !
              ! Open or closed point - the decision has already been taken in open_close and is stored in gisca
              ! gisca =1 for inflow                   
              if( gisca(ipoin) > 0 ) then ! inflow boundary
                 do idime =1, ndime
                    if( abs(kfl_fixno_nsi(idime,ipoin)) == 7_ip ) then ! if any
                       kfl_fixno_nsi(idime,ipoin) = 7                  ! fix
                       ! correct initial guess for outer iterations - this will be used by tem & tur_updbcs
                       veloc(idime,ipoin,1) = bvess_nsi(idime,ipoin,1)  
                       kfl_fixpr_nsi(1,ipoin) = 0 
                    end if
                 end do
                 ! for fixpr I don not care to use a 7 as it is what happens with velocity that decides
              else ! outflow boundary
                 do idime =1, ndime                    
                    if( abs(kfl_fixno_nsi(idime,ipoin)) == 7_ip ) then
                       kfl_fixno_nsi(idime,ipoin) = -7 !free
                       kfl_fixpr_nsi(1,ipoin) = 1   
                    end if
                 end do
              end if
           end do

           call memgen(3_ip,npoin,0_ip)   !deallocate gisca

        end if
        !
        ! BVESS: Dirichlet with specific code 9,8,6,5
        !
        do ipoin = 1,npoin

           if( kfl_fixno_nsi(1,ipoin) == 9 ) then
              !
              ! U(t) = U_inf - w x r
              !
              bvess_nsi(:,ipoin,1)=bvess_nsi(:,ipoin,2)              
              if(ndime==2) then
                 bvess_nsi(1,ipoin,1) = bvess_nsi(1,ipoin,1) &
                      + fvela_nsi(3) * coord(2,ipoin)
                 bvess_nsi(2,ipoin,1) = bvess_nsi(2,ipoin,1) &
                      - fvela_nsi(3) * coord(1,ipoin)
              else
                 bvess_nsi(1,ipoin,1) = bvess_nsi(1,ipoin,1) &
                      + fvela_nsi(3) * coord(2,ipoin) &
                      - fvela_nsi(2) * coord(3,ipoin)
                 bvess_nsi(2,ipoin,1) = bvess_nsi(2,ipoin,1) &
                      - fvela_nsi(3) * coord(1,ipoin) &
                      + fvela_nsi(1) * coord(3,ipoin)
                 bvess_nsi(3,ipoin,1) = bvess_nsi(3,ipoin,1) &
                      - fvela_nsi(1) * coord(2,ipoin) &
                      + fvela_nsi(2) * coord(1,ipoin)        
              end if

           else if( abs(kfl_fixno_nsi(1,ipoin)) == 8 ) then
              !
              !             R(theta)
              ! No inertial ========> Lab
              !
              !             R(-theta)
              !         Lab ========> No inertial
              !
              ! u(t) = R(-theta).U_inf - w x r
              ! U(t) = R(theta).( u + w x r ) (used in INERTIAL postprocess, to be used with postprocess of Mesh rotation)
              !
              call rotmat(ndime,-cutim*fvnoa_nsi,fvdia_nsi,rot_matrix)
              bvess_new = 0.0_rp
              do jdime = 1,ndime
                 bvess_new(1:ndime) = bvess_new(1:ndime) + rot_matrix(1:ndime,jdime) * bvess_nsi(jdime,ipoin,1)
              end do
              if( ndime == 2 ) then
                 bvess_new(1) = bvess_new(1) + fvela_nsi(3) * coord(2,ipoin)
                 bvess_new(2) = bvess_new(2) - fvela_nsi(3) * coord(1,ipoin)
              else
                 bvess_new(1) = bvess_new(1) + fvela_nsi(3) * coord(2,ipoin) - fvela_nsi(2) * coord(3,ipoin)
                 bvess_new(2) = bvess_new(2) - fvela_nsi(3) * coord(1,ipoin) + fvela_nsi(1) * coord(3,ipoin)
                 bvess_new(3) = bvess_new(3) - fvela_nsi(1) * coord(2,ipoin) + fvela_nsi(2) * coord(1,ipoin)        
              end if
              !
              !
              ! Detrermine inflow/outflow
              ibopo = lpoty(ipoin)
              if( ibopo > 0 ) then
                 xnorm = sqrt(dot_product(bvess_new(1:ndime),bvess_new(1:ndime)))
                 if( xnorm > zeror ) xnorm = 1.0_rp / xnorm
                 udotn = dot_product(exnor(1:ndime,1,ibopo),bvess_new(1:ndime))
                 udotn = udotn * xnorm                
                 if( udotn <= inflow_cosine_nsi ) then  !  Inflow
                    kfl_fixno_nsi(1:ndime,ipoin) =  8
                    kfl_fixpr_nsi(1,ipoin)       =  0
                    bvess_nsi(1:ndime,ipoin,1)   =  bvess_new(1:ndime)
                 else                                   !  Outflow
                    kfl_fixno_nsi(1:ndime,ipoin) = -8
                    kfl_fixpr_nsi(1,ipoin)       =  1
                 end if
              end if

           else if( kfl_fixno_nsi(1,ipoin) == 5 ) then
              !
              ! SEM: synthetic eddy method
              !
              !    bvess_nsi(1,ipoin,1) = 3.7_rp
              !    do idime = 2,ndime
              !       bvess_nsi(idime,ipoin,1) = 0.0_rp                 
              !    end do             

           else if( kfl_fixno_nsi(1,ipoin) == 6 ) then
              !
              ! SPACE-TIME boundary condition read from file
              !
              !if (ittim > 0 ) then 

              !   ipass = int( mod(cutim,nbtdt_nsi*(nbtim_nsi-1_ip)) / nbtdt_nsi ,ip)   ! Location of the plane
              !   alpha = mod(cutim,nbtdt_nsi) / nbtdt_nsi                              ! Interpolation between planes

              !   do idime = 1,ndime
              !      kfl_fixno_nsi(idime,ipoin) = 6
              !      bvess_nsi(idime,ipoin,1) =  alpha * bnval_nsi(ipass*nbval_nsi+iboun_nsi(ipoin),idime) &
              !           + (1.0_rp - alpha) * bnval_nsi((ipass + 1_ip)*nbval_nsi+iboun_nsi(ipoin),idime)
              !   end do
              !else
              !   do idime = 1,ndime
              !      kfl_fixno_nsi(idime,ipoin) = 6
              !      bvess_nsi(idime,ipoin,1) = bnval_nsi(iboun_nsi(ipoin),idime)
              !   end do
              !end if
              !!!if (ittim > 0 ) then

              itable = iboun_nsi(ipoin,2)
              alpha = mod(cutim, nbtdt_nsi(itable)) / nbtdt_nsi(itable)    ! weight of right value  (mod(t,dt)/dt) 
              ibtim = nbtim_pos_nsi(itable+1) - nbtim_pos_nsi(itable)
              i_left = mod( int( cutim / nbtdt_nsi(itable), ip), ibtim)    ! Index of left value (0..nbtim_nsi-1)
              if (i_left /= (ibtim-1_ip) ) then                            ! Index of right value (0..nbtim_nsi-1)
                 i_right = i_left + 1_ip
              else
                 i_right = 0_ip
              end if

              !
              ! indeces in table e.g.:
              ! 
              ind_left  = nbtim_nod_pos_nsi(itable) + i_left  * &
                             (nbnod_pos_nsi(itable+1)-nbnod_pos_nsi(itable)) + iboun_nsi(ipoin,1)
              ind_right = nbtim_nod_pos_nsi(itable) + i_right * &
                             (nbnod_pos_nsi(itable+1)-nbnod_pos_nsi(itable)) + iboun_nsi(ipoin,1)

              do idime = 1,ndime
                 kfl_fixno_nsi(idime,ipoin) = 6
                 bvess_nsi(idime,ipoin,1) =  (1.0_rp-alpha) * bnval_nsi(ind_left ,idime) &
                      +  alpha          * bnval_nsi(ind_right,idime)
              end do

              !!!else
              !!!   do idime = 1,ndime
              !!!      kfl_fixno_nsi(idime,ipoin) = 6
              !!!      bvess_nsi(idime,ipoin,1) = bnval_nsi(iboun_nsi(ipoin)  ,idime)
              !!!   end do
              !!!end if

           end if
        end do
        !
        ! BVNAT: Neumann space/time or time functions
        !
        do iboun = 1,nboun

           if( kfl_fixbo_nsi(iboun) /= 0 ) then

              center_gravity = 0.0_rp
              do inodb= 1,nnode(ltypb(iboun))
                 center_gravity(1:ndime) = center_gravity(1:ndime) + coord(1:ndime,lnodb(inodb,iboun))/real(nnode(ltypb(iboun)),rp)
              end do

              ifunc = kfl_funbo_nsi(iboun)
              itype = kfl_funtb_nsi(iboun)

              if( itype == FUNCTION_TIME ) then
                 !
                 ! Function: U(t)=f(t)*U(0)
                 !                 
                 venew = bvnat_nsi(1,iboun,2) *funcre(   &
                      time_function(ifunc) % parameters, & 
                      time_function(ifunc) % npara,      &
                      time_function(ifunc) % kfl_type,   &
                      cutim)
                 if( kfl_timei_nsi /= 0 .and. kfl_tiacc_nsi == 2 .and. kfl_tisch_nsi == 1 ) then
                    bvnat_nsi(1,iboun,1) = 0.50_rp * ( venew + bvnat_nsi(1,iboun,1) )
                 else
                    bvnat_nsi(1,iboun,1) = venew                 
                 end if

              else if( itype == FUNCTION_SPACE_TIME ) then
                 !
                 ! Field
                 !                 
                 call ker_space_time_function(&
                      ifunc,center_gravity(1),center_gravity(2),center_gravity(ndime),cutim,bvnew)
                 bvnew = bvnew * bvnat_nsi(1,iboun,2)
                 if( kfl_timei_nsi /= 0) then
                    if ( (kfl_tiacc_nsi == 2 .and. kfl_tisch_nsi == 1) .or. (kfl_tiacc_nsi == 4 .and. kfl_tisch_nsi == 4) ) then
                       bvold = bvnat_nsi(1,iboun,1)
                       bvnat_nsi(1,iboun,1) = 0.50_rp*(bvnew+bvold)
                    else
                       bvnat_nsi(1,iboun,1) = bvnew
                    end if
                 end if

              else if( itype == FUNCTION_WINDKESSEL ) then
                 !
                 ! Windkessel
                 !                 
                 p_in  = 0.0_rp                 
                 dummr = 0.0_rp
                 call  mod_windk_interact(                                        &
                      & dummr, p_in,'NASTIN', kfl_codbo(iboun), &
                      & IS_EXPLICIT=NSI_FRACTIONAL_STEP, OPT_direction='OU')
                 bvnat_nsi(5,iboun,1) = p_in 

              else if (kfl_cadan_nsi == 1) then

                 if (kfl_fixbo_nsi(iboun) == 2) bvnat_nsi(1,iboun,1) = press_cadan_nsi

              end if
           end if
        end do

     end if
     !
     ! Flow rate constraint
     !
     call nsi_impose_flow_rates()
     !
     ! Impose boundary conditions on VELOC
     !
     call local_basis_global_to_local(kfl_fixrs_nsi,veloc,LAST_COMPONENT=1_ip) ! Global to local
     do ipoin = 1,npoin
        do idime = 1,ndime
           if( kfl_fixno_nsi(idime,ipoin) > 0 ) then
              veloc(idime,ipoin,1) = bvess_nsi(idime,ipoin,1)
           end if
        end do
     end do
     call local_basis_local_to_global(kfl_fixrs_nsi,veloc,LAST_COMPONENT=1_ip) ! Local to global
     !
     ! Exact solution
     !
     if( kfl_exacs_nsi /= 0 .and. kfl_conbc_nsi == 0 ) call nsi_exaerr(1_ip)
     !
     ! Dirichlet coupling
     !
     if( INOTEMPTY ) then
        call couplings_impose_dirichlet(solve(1),veloc(:,:,1))
        call couplings_impose_dirichlet(solve(2),press(:,1))
     end if

  case( ITASK_BEGITE )

     !-------------------------------------------------------------------
     !
     ! Before a global iteration
     !  
     !-------------------------------------------------------------------
     !
     ! Add mesh velocity to prescribed velocity
     !
     if( associated(velom) .and. kfl_velom_nsi /= 0 ) then
        call local_basis_global_to_local(kfl_fixrs_nsi,velom)           ! Global to local     
        select case ( kfl_conbc_nsi )
        case ( 1_ip )
           !
           ! Constant boundary conditions
           !
           do ipoin = 1,npoin
              bvess_nsi(:,ipoin,1) = velom(:,ipoin)
           end do
        case ( 0_ip )
           !
           ! Not to overwrite the bc from a function
           !
           do ipoin = 1,npoin
              if( kfl_funno_nsi(ipoin) == 0 ) then
                 bvess_nsi(:,ipoin,1) = bvess_nsi(:,ipoin,2) + velom(:,ipoin)
              end if
           end do
        end select
        call local_basis_local_to_global(kfl_fixrs_nsi,velom)           ! Local to global
     end if

  case( ITASK_BEGINN )

     !-------------------------------------------------------------------
     !
     ! Before an inner iteration
     !  
     !-------------------------------------------------------------------
     !
     ! Pressure level mainly for fluid structure interactions
     !
     if( abs(val_press_lev_nsi) > zeror ) then
        press_lev_nsi = val_press_lev_nsi
     end if
     if( kfl_press_lev_nsi /= 0 ) then
        ii = maths_findloc(lbsec,kfl_press_lev_nsi)
        if( ii /= 0 ) then
           press_lev_nsi = val_press_lev_nsi - vbset(1,ii)
        else
           call runend('NSI_UPDBCS: WRONG SET TO COMPUTE PRESSURE LEVEL')
        end if
     end if
     if( abs(val_press_lev_nsi) > zeror .or. kfl_press_lev_nsi /= 0 ) then
        solve(2) % reaction_level = press_lev_nsi
     end if
     !
     ! If Dynamic pressure cond. exist, compute average dynamic pressure
     !
     if( ittot_nsi >= itebc_nsi ) call nsi_bouave()
     !
     ! Basilev condition
     !
     if( kfl_conbc_nsi == 0 ) then 

        if( kfl_exist_fib20_nsi == 1 ) then
           call nsi_outflow_mass()
           do ifunc = 1,max_windk_systems
              if( ( windk_systems(ifunc) % sysid /= 0_ip ) .and. ( itinn(modul)==1_ip ) ) then
                 call  mod_windk_interact(                                                &
                      & outflow_mass(windk_systems(ifunc) % iflow_nsi), &
                      & dummr, 'NASTIN', windk_systems(ifunc) % tag_in, &
                      & IS_EXPLICIT=NSI_FRACTIONAL_STEP, OPT_direction='BI')
              else
                 exit
              end if
           end do
           do iboun = 1,nboun           
              if( kfl_fixbo_nsi(iboun) == 20 ) then

                 ifunc = kfl_funbo_nsi(iboun)
                 itype = kfl_funtb_nsi(iboun)
                 
                 if( itype == FUNCTION_WINDKESSEL ) then
                    !
                    ! Windkessel
                    !                 
                    p_in  = 0.0_rp                 
                    dummr = 0.0_rp
                    call  mod_windk_interact(                                        &
                         & dummr, p_in,'NASTIN', kfl_codbo(iboun), &
                         & IS_EXPLICIT=NSI_FRACTIONAL_STEP, OPT_direction='OU')
                    bvnat_nsi(5,iboun,1) = p_in
                 else
                    !                       
                    ! Normal outflow
                    !
                    iflow = int(bvnat_nsi(4,iboun,1),ip)
                    bvnat_nsi(5,iboun,1) = bvnat_nsi(2,iboun,1) * outflow_mass(iflow)
                 end if
                 
              end if
           end do
        end if
     end if
     !
     ! Update pressure Dirichlet condition for fractional step
     !
     call nsi_bougra()
     !
     ! Impose boundary conditions on PRESS
     !
     if( ( kfl_exist_fib02_nsi == 1 .or. kfl_exist_fib20_nsi == 1 ) .and. NSI_FRACTIONAL_STEP ) then
        do ipoin = 1,npoin
           if( kfl_fixpp_nsi(1,ipoin) > 0 ) then
              press(ipoin,1) = bpess_nsi(1,ipoin,1)
           end if
        end do
     end if

  case(4_ip)

     !-------------------------------------------------------------------
     !
     ! Enforce boundary conditions
     !
     !-------------------------------------------------------------------

     if( INOTEMPTY ) then
        call local_basis_global_to_local(kfl_fixrs_nsi,unkno) ! Global to local
        if( NSI_MONOLITHIC ) then
           do ipoin = 1,npoin
              itotv = (ipoin-1)*(ndime+1)
              do idime = 1,ndime
                 itotv = itotv + 1
                 if(  kfl_fixno_nsi(idime,ipoin) == 1 .or. &
                      kfl_fixno_nsi(idime,ipoin) == 8 .or. &
                      kfl_fixno_nsi(idime,ipoin) == 9 .or. &
                      kfl_fixno_nsi(idime,ipoin) == 6   ) then
                    unkno(itotv) = bvess_nsi(idime,ipoin,1)
                 end if
              end do
           end do
        else
           itotv = 0
           do ipoin = 1,npoin
              do idime = 1,ndime
                 itotv = itotv + 1
                 if(  kfl_fixno_nsi(idime,ipoin) == 1 .or. &
                      kfl_fixno_nsi(idime,ipoin) == 8 .or. &
                      kfl_fixno_nsi(idime,ipoin) == 9 .or. &
                      kfl_fixno_nsi(idime,ipoin) == 6   ) then
                    unkno(itotv) = bvess_nsi(idime,ipoin,1)
                 end if
              end do
           end do
        end if
        call local_basis_local_to_global(kfl_fixrs_nsi,unkno) ! Local to global 
     end if

  end select

end subroutine nsi_updbcs



