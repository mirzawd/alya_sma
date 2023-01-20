!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_inibcs.f90
!> @author  Guillaume Houzeaux
!> @brief   Impose boundary conditions
!> @details Impose boundary conditions\n
!>    For conditions on nodes, bvess_nsi(idime,ipoin,1) contains the\n
!>    velocity value. \n
!>    The different codes for kfl_fixno_nsi(ipoin) are:\n
!>    =  1 ... Dirichlet\n
!>    =  0 ... Free or initial\n
!>    =  2 ... Adaptive Inflow (Dirichlet) or outflow (Neumann)\n
!>    =  3 ... Inertial Dirichlet (same as 1 but + mesh velocity)\n
!>    =  5 ... Synthetic eddy method (SEM) \n
!>    =  6 ... u(y)=0\n
!>    =  7 ... u(z)=0\n
!>    =  8 ... Interpolated Dirichlet\n
!>    =  9 ... \n
!>
!>    For conditions on boundaries, bvnat_nsi(iboun) is allocated ONLY if\n
!>    a condition is imposed on it. Overmore, its length depends on the\n
!>    condition type. The different codes for kfl_fixbo_nsi(iboun) are:\n
!>    =  1 ... Dirichlet ................... u\n
!>    =  2 ... Pressure imposed  ........... sig.n=-p n\n
!>    =  3 ... Wall law .................... sig.n.t=-rho*U_*^2*(u_tan-u_fix_tan)/|u_tan-u_fix_tan|\n
!>    =  4 ... Symmetry/Slip wall .......... sig.n.t=0, u.n=0\n
!>    =  5 ... Dynamic pressure ............ sig.n=-1/2*rho*u^2\n
!>    =  6 ... Open flow ................... Assemble 2*mu*Sym(grad(u).n and prescribe p\n
!>    =  7 ... No slip wall ................ u=0\n
!>    =  8 ... Interpolated fixed .......... u=I*u_background\n
!>    =  9 ... Non-inertial ................ u-wxr\n
!>    = 10 ... Dyn. press.+open flow ....... idem\n
!>    = 11 ... Inflow depending on angle ... u\n
!>    = 20 ... Stable outflow .............. -p0 n - C * \int_S u.n ds + rho*(u.n)_ u * beta \n
!>    At the end of the subroutine conditions on boundaries are\n
!>    transfered to conditions on nodes.!!\n
!> @} 
!-----------------------------------------------------------------------
subroutine nsi_inibcs()
  use def_parame
  use def_inpout
  use def_master
  use def_kermod
  use def_domain
  use def_nastin
  use mod_nastin,         only : nastin_solution_strategy
  use mod_communications, only : PAR_MAX, PAR_INTERFACE_NODE_EXCHANGE
  implicit none
  integer(ip)  :: ipoin,pnodb,iboun,inodb,idime,ibopo,pblty
  !
  ! Define parameters according to the solution strategy
  !
  call nastin_solution_strategy()

  !-------------------------------------------------------------
  !
  ! Allocate memory for the vectors needed to define the BC's 
  !
  !-------------------------------------------------------------

  call nsi_membcs(1_ip) 

  if( INOTMASTER ) then

     !-------------------------------------------------------------
     !
     ! Geometrical codes
     !
     !-------------------------------------------------------------

     call nsi_geobcs()

     !-------------------------------------------------------------
     !
     ! Node codes
     !
     !-------------------------------------------------------------
 
     if( kfl_icodn > 0 ) then
        !
        ! Velocity
        !
        if( kfl_conbc_nsi == 0 ) then
           iffun      =  1
           kfl_funno  => kfl_funno_nsi
           kfl_funtn  => kfl_funtn_nsi
        else
           iffun      =  0
        end if
        ifbop     =  0
        ifloc     =  1
        skcos_nsi => skcos
        kfl_fixrs => kfl_fixrs_nsi
        kfl_fixno => kfl_fixno_nsi
        bvess     => bvess_nsi(:,:,1)
        tncod     => tncod_nsi(1:) 
        call reacod(IMPOSE_NODE_CODES)
        !
        ! Pressure
        !
        if( kfl_conbc_nsi == 0 ) then
           iffun  =  1
        else
           iffun  =  0
        end if
        ifbes     =  1
        ifbop     =  0
        ifloc     =  0
        kfl_fixno => kfl_fixpp_nsi
        bvess     => bpess_nsi(:,:,1)
        tncod     => tncod_nsi(2:)
        call reacod(IMPOSE_NODE_CODES)
        !
        ! Divergence correction
        !
        if( kfl_divcorrec_nsi /= 0 ) then
           ifbes     =  0
           kfl_fixno => kfl_fixno_div_nsi
           tncod     => tncod_nsi(3:)
           call reacod(IMPOSE_NODE_CODES)
        end if

     end if

     !-------------------------------------------------------------
     !
     ! Boundary codes
     !
     !-------------------------------------------------------------

     if( kfl_icodb > 0 ) then
        if( kfl_conbc_nsi == 0 ) then
           iffun      =  1
           kfl_funbo  => kfl_funbo_nsi
           kfl_funtb  => kfl_funtb_nsi
        end if
        kfl_fixbo => kfl_fixbo_nsi
        bvnat     => bvnat_nsi(:,:,1)
        tbcod     => tbcod_nsi(1:)
        !tncod     => momod(ID_NASTIN) % tncod(1:)    
        call reacod(IMPOSE_BOUNDARY_CODES)
     end if

     !-------------------------------------------------------------
     !
     ! Exact solution
     !
     !-------------------------------------------------------------

     call nsi_exaerr(1_ip)

     !-------------------------------------------------------------
     !
     ! Check if there is a local prescription
     !
     !-------------------------------------------------------------

     nodes: do ipoin = 1,npoin
        ibopo = lpoty(ipoin)
        if( ibopo > 0 ) then
           if( kfl_fixrs_nsi(ipoin) /= 0 ) then
              kfl_local_nsi = 1
              exit nodes
           end if
        end if
     end do nodes

     !------------------------------------------------------------------
     !
     ! Put no slip to wall law boundaries if delta_nsi is negative
     !
     !------------------------------------------------------------------

     if( delta_dom > zensi .or. delta_dom < -zensi ) then
        delta_nsi = delta_dom
     end if

     if(  (abs(delta_nsi) <= zensi) .and. (kfl_delta /= 1) ) then

        do iboun = 1,nboun
           if( kfl_fixbo_nsi(iboun) == 3 ) kfl_fixbo_nsi(iboun) = 4
        end do
     end if

     if( ( delta_nsi <- zensi ) .or. (kfl_delta == -1) )  then

        do iboun = 1,nboun
           if( kfl_fixbo_nsi(iboun) == 3 ) then
              kfl_fixbo_nsi(iboun) = 7
              pnodb=nnode(ltypb(iboun))
              do inodb=1,pnodb
                 ipoin=lnodb(inodb,iboun) 
                 if(  kfl_fixno_nsi(    1,ipoin)==-1.and.&
                      kfl_fixno_nsi(    2,ipoin)==-1.and.&
                      kfl_fixno_nsi(ndime,ipoin)==-1) then
                    ibopo=lpoty(ipoin)
                    kfl_fixno_nsi(1:ndime,ipoin)=1
                    bvess_nsi(1:ndime,ipoin,1)=0.0_rp
                    kfl_fixrs_nsi(ipoin)=0
                    if(kfl_conbc_nsi==0) kfl_funno_nsi(ipoin)=0
                 end if
              end do
           end if
        end do
     end if
     delta_nsi = max(delta_nsi,0.0_rp)
     if (kfl_delta == -1)  delta_nsi = 0.0_rp
     !
     !Prepare Poiseuille bcs
     !
     if (kfl_initi_nsi == 6 .or. kfl_initi_nsi == 7 ) then
        call runend('NSI_INIBCS: THIS DOES NOT WORK IN PARALLEL')
        do iboun = 1,nboun
           if( kfl_fixbo_nsi(iboun) == 16 ) then
              kfl_fixbo_nsi(iboun) = 0
              pnodb=nnode(ltypb(iboun))
              do inodb=1,pnodb
                 ipoin=lnodb(inodb,iboun) 
                 ibopo=lpoty(ipoin)
                 kfl_fixno_nsi(1:ndime,ipoin)=16
                 kfl_fixrs_nsi(ipoin)=0
                 if(kfl_conbc_nsi==0) kfl_funno_nsi(ipoin)=0
              end do
           endif
        end do
     endif
     !
     ! Put -1 flag to 0
     !
     do ipoin = 1,npoin
        do idime = 1,ndime
           if( kfl_fixno_nsi(idime,ipoin) == -1 ) kfl_fixno_nsi(idime,ipoin) = 0
        end do
     end do

     !------------------------------------------------------------------
     !
     ! Non-inertial boundary condition: u - w x r
     !
     !------------------------------------------------------------------

     if( kfl_conbc_nsi == 1 ) then
        call nsi_updfor()
        do ipoin = 1,npoin
           if( kfl_fixno_nsi(1,ipoin) == 9 ) then
              if( ndime == 2 ) then
                 bvess_nsi(1,ipoin,1) = bvess_nsi(1,ipoin,1)&
                      + fvela_nsi(3) * coord(2,ipoin)            ! u+wz*y
                 bvess_nsi(2,ipoin,1) = bvess_nsi(2,ipoin,1)&  
                      - fvela_nsi(3) * coord(1,ipoin)            ! v-wz*x
              else
                 bvess_nsi(1,ipoin,1) = bvess_nsi(1,ipoin,1)&
                      + fvela_nsi(3) * coord(2,ipoin)&           ! u+wz*y
                      - fvela_nsi(2) * coord(3,ipoin)            !  -wy*z
                 bvess_nsi(2,ipoin,1) = bvess_nsi(2,ipoin,1)&
                      - fvela_nsi(3) * coord(1,ipoin)&           ! v-wz*x
                      + fvela_nsi(1) * coord(3,ipoin)            !  +wx*z
                 bvess_nsi(3,ipoin,1) = bvess_nsi(3,ipoin,1)&
                      - fvela_nsi(1) * coord(2,ipoin)&           ! w-wx*y
                      + fvela_nsi(2) * coord(1,ipoin)            !  +wy*x
              end if
           end if
        end do
     end if

     !------------------------------------------------------------------
     !
     ! Condition with angle: modify pressure b.c.
     !
     !------------------------------------------------------------------

     if( NSI_SCHUR_COMPLEMENT ) then
        do ipoin = 1,npoin
           ibopo = lpoty(ipoin)
           if(  ibopo > 0                         .and. &
                kfl_fixno_nsi(    1,ipoin) == 999 .and. &
                kfl_fixno_nsi(    2,ipoin) == 999 .and. &
                kfl_fixno_nsi(ndime,ipoin) == 999 ) then
              kfl_fixpr_nsi(1,ipoin) = 1
           end if
        end do
        do iboun = 1,nboun
           inodb = 0
           pnodb = nnode(ltypb(iboun))
           do while( inodb < pnodb )
              inodb = inodb+1
              ipoin = lnodb(inodb,iboun)
              if(  kfl_fixno_nsi(    1,ipoin) == 999 .and. &
                   kfl_fixno_nsi(    2,ipoin) == 999 .and. &
                   kfl_fixno_nsi(ndime,ipoin) == 999 ) then
                 do inodb = 1,pnodb
                    ipoin = lnodb(inodb,iboun)
                    ibopo = lpoty(ipoin)
                    !kfl_fixpr_nsi(1,ipoin) = 1
                 end do
              end if
           end do
        end do
     end if

     !------------------------------------------------------------------
     !
     ! Recover Neumann b.c. for velocity
     !
     !------------------------------------------------------------------

     do ipoin = 1,npoin
        do idime = 1,ndime
           if( kfl_fixno_nsi(idime,ipoin) == 999 ) then
              kfl_fixno_nsi(idime,ipoin) = 0
           end if
        end do
     end do

     !------------------------------------------------------------------
     !
     ! Non-constant boundary conditions
     !
     !------------------------------------------------------------------

     if( kfl_conbc_nsi == 0 ) then
        do ipoin = 1,npoin
           do idime = 1,ndime
              bvess_nsi(idime,ipoin,2) = bvess_nsi(idime,ipoin,1)
           end do
           bpess_nsi(1,ipoin,2) = bpess_nsi(1,ipoin,1)
        end do
        do iboun = 1,nboun
           bvnat_nsi(:,iboun,2) = bvnat_nsi(:,iboun,1)
        end do
     end if
     
     !------------------------------------------------------------------
     !
     ! Check if there is a stable outflow condition and check data
     !
     !------------------------------------------------------------------

     iboun = 0
     do while( iboun < nboun )
        iboun = iboun + 1
        if( kfl_fixbo_nsi(iboun) == 20 ) then
           kfl_exist_fib20_nsi = 1
           if( int(bvnat_nsi(4,iboun,1),ip) < 1_ip ) then
              call runend('NSI_INIBCS: WRONG STABLE OUTFLOW CONDITION DATA')
           end if
        end if
     end do

     !------------------------------------------------------------------
     !
     ! Check if pressure is imposed at values /= 0
     ! Used for initial pressure profile in nsi_iniunk
     !
     !------------------------------------------------------------------

     iboun = 0
     do while( iboun < nboun )
        iboun = iboun + 1
        if( kfl_fixbo_nsi(iboun) == 2 .or. kfl_fixbo_nsi(iboun) == 6 .or. kfl_fixbo_nsi(iboun) == 20 ) then
           if( bvnat_nsi(1,iboun,1) /= 0.0_rp ) then
              kfl_inipr_nsi = 1
              iboun = nboun
           end if
        end if
     end do
     !
     ! Obtain lpset(nbopo) when INTERNAL_FORCES: RESID is used
     !     
     if ( kfl_intfo_nsi >= 1 ) then
        allocate(lbpse(npoin))
        lbpse = 0_ip
        if ( associated(lbset) ) then
           do iboun = 1,nboun
              pblty = ltypb(iboun) 
              pnodb = nnode(pblty)
              do inodb = 1,pnodb
                 ipoin = lnodb(inodb,iboun)
                 lbpse(ipoin) = max (lbpse(ipoin),lbset(iboun)) 
              end do
           end do
        end if        
        call PAR_INTERFACE_NODE_EXCHANGE(lbpse,'MAX','IN MY CODE')
     end if
     
  end if
  !
  ! Do fixity 7 exist ?
  !
  if( kfl_conbc_nsi == 0 ) then
     kfl_exist_fixi7_nsi = 0
     do ipoin = 1,npoin
        if( kfl_fixno_nsi(1,ipoin) == 7 ) kfl_exist_fixi7_nsi = 1
     end do
  end if
  !
  ! Do fixbo = 2 exist ?
  !
  kfl_exist_fib02_nsi = 0
  do iboun = 1,nboun
     if( kfl_fixbo_nsi(iboun) == 2 ) kfl_exist_fib02_nsi = 1
  end do
  !
  ! Should be broadcast because subdomain can modify it while others not
  !
  call PAR_MAX(kfl_inipr_nsi      ,'IN MY CODE')
  call PAR_MAX(kfl_exist_fixi7_nsi,'IN MY CODE')
  call PAR_MAX(kfl_exist_fib20_nsi,'IN MY CODE')
  call PAR_MAX(kfl_exist_fib02_nsi,'IN MY CODE')

end subroutine nsi_inibcs
