!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_iniunk()
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_iniunk
  ! NAME 
  !    tur_iniunk
  ! DESCRIPTION
  !    This routine sets up the initial condition for the turbulence
  !    variables
  ! USED BY
  !    tur_begste
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_kermod
  use mod_ker_proper
  use def_domain
  use def_turbul
  use mod_chktyp,           only : check_type
  use mod_ker_space_time_function
  use mod_ker_regularization, only : inv_regul_k, inv_regul_e, kfl_regularization
  use mod_communications, only : PAR_MAX
  implicit none
  integer(ip)          :: ipoin,iturb,kfl_value,dummi,ifunc
  real(rp)             :: maxtu,rho(1),mu(1),nu,nup,nut
  real(rp),    target  :: turbm(nturb_tur)
  
  
  avtim_tur = 0.0_rp  ! Accumulated time for time-averaging variables
!  avkey_tur = 0.0_rp
!  avome_tur = 0.0_rp
!  avtvi_tur = 0.0_rp

  if( kfl_rstar == 0 ) then
        !
        ! TURMU: Initial turbulent viscosity
        !
        if( INOTMASTER ) then
           call ker_proper('VISCO','NPOIN',dummi,dummi,turmu)
           if( nutnu_tur >= 0.0_rp ) then
              do ipoin = 1,npoin 
                 turmu(ipoin) = turmu(ipoin) * nutnu_tur
              end do
           end if
           
        end if
        !
        ! Treat special cases
        !
        call tur_vortic()                            ! Vorticity magnitude
        call tur_adapti()                            ! Adaptive b.c. according to velocity
        call tur_frivel()                            ! Friction velocity
        call tur_bouave(1_ip)                        ! Average velocity
        call tur_updbcs(TUR_BEFORE_GLOBAL_ITERATION) ! Boundary conditions
        !
        ! Look for maximum values
        !    
        if( TUR_SPALART_ALLMARAS ) then

           if( INOTMASTER ) then 
              do ipoin = 1,npoin
                 if(  kfl_fixno_tur(1,ipoin,1) <= 0 .or.&
                      kfl_fixno_tur(1,ipoin,1) == 3 .or.&
                      kfl_fixno_tur(1,ipoin,1) == 4 ) then                 
                    call ker_proper('DENSI','IPOIN',ipoin,dummi,rho(1:))
                    call ker_proper('VISCO','IPOIN',ipoin,dummi,mu(1:))
                    nu  = mu(1)/rho(1) 
                    nut = turmu(ipoin)/rho(1)
                    call tur_spanut(1_ip,nu,nup,nut)
                    untur(1,ipoin,1) = nup
                 else
                    untur(1,ipoin,1) = bvess_tur(1,ipoin,1)
                 end if
              end do
           end if

        else

           if( INOTMASTER ) then
              turbm = 0.0_rp
              do ipoin = 1,npoin
                 do iturb = 1,nturb_tur
                    if( kfl_fixno_tur(1,ipoin,iturb) /= 3 .and. kfl_fixno_tur(1,ipoin,iturb) /=4 &
                         .and. bvess_tur(1,ipoin,iturb) > turbm(iturb) ) &
                         turbm(iturb) = bvess_tur(1,ipoin,iturb)
                 end do
              end do
           end if
           !
           ! Look for maximum over whole mesh
           !
           call PAR_MAX(nturb_tur,turbm)
           !
           ! Initial values for turbulence variables
           !
           if( INOTMASTER ) then

              if( kfl_inifi_tur(1) > 0 ) then
                 !
                 ! initialization from fields
                 !
                 do iturb = 1,nturb_tur
                    do ipoin = 1,npoin
                       untur(iturb,ipoin,1) = xfiel(-nfiel_tur(iturb))%a(1,ipoin,1)
                    end do
                 end do
              else if( kfl_initi_tur == 1 ) then 
                 !
                 ! All from uniform value
                 !
                 do ipoin = 1,npoin
                    do iturb = 1,nturb_tur
                       untur(iturb,ipoin,1) = xinit_tur(iturb)
                    end do
                 end do 
                  
              else if( kfl_initi_tur == 2 ) then 
                 !
                 ! All from bvess 
                 !
                 do ipoin = 1,npoin
                    do iturb = 1,nturb_tur
                       untur(iturb,ipoin,1) = bvess_tur(1,ipoin,iturb)
                    end do
                 end do 
                  
              else if( kfl_initi_tur == -1 ) then
                 !
                 ! Value from a function
                 ! 
                 do iturb = 1,nturb_tur
                    kfl_value = kfl_valbc_tur(iturb)
                    call check_type(xfiel,kfl_value,1_ip,npoin) ! Check if value function exist
                    do ipoin = 1,npoin
                       untur(iturb,ipoin,1) = xfiel(kfl_value) % a(1,ipoin,1)
                    end do
                 end do
                  
              else if( kfl_initi_tur == 3 ) then 
                 ! 
                 ! Space time function
                 !
                do iturb = 1,nturb_tur
                   ifunc = int(xinit_tur(iturb),ip)
                   if( ifunc > 0 ) then
                      do ipoin = 1,npoin
                         call ker_space_time_function(&
                              ifunc,coord(1,ipoin),coord(2,ipoin),coord(ndime,ipoin),cutim, untur(iturb,ipoin,1))
                      end do
                   end if
                end do
                
              else 
                 !
                 ! Based on maximum values
                 !
                 do ipoin = 1,npoin
                    do iturb = 1,nturb_tur
                       if(  kfl_fixno_tur(1,ipoin,iturb) <= 0 .or.&
                            kfl_fixno_tur(1,ipoin,iturb) == 3 .or.&
                            kfl_fixno_tur(1,ipoin,iturb) == 4 ) then
                          untur(iturb,ipoin,1) = turbm(iturb)
                       else 
                          untur(iturb,ipoin,1) = bvess_tur(1,ipoin,iturb)
                       end if
                    end do
                 end do

              end if

           end if

        end if
        !
        ! When solving a manufactured solution, impose exact Dirichlet bc
        !
        !if( kfl_exacs_tur /= 0 ) then
        !   call tur_manufactured_Dirichlet_condition()
        !end if

        if( INOTMASTER ) then
           !
           ! Initialize turbulent unknows ( nprev_tur -> 1_ip)
           !
           call tur_updunk(ITASK_INIUNK)
           !
           ! If initial values are zero
           !
           maxtu = 0.0_rp
           do iturb = 1,nturb_tur
              if( turbm(iturb) > maxtu ) maxtu = turbm(iturb)
           end do
           if( maxtu == 0.0_rp ) then
              do ipoin = 1,npoin 
                 call ker_proper('DENSI','IPOIN',ipoin,dummi,rho(1:))
                 call ker_proper('VISCO','IPOIN',ipoin,dummi,mu(1:))
                 turmu(ipoin) = 10.0_rp*mu(1)
              end do
           end if

        end if
        
  else
     !
     ! Restart
     !
     avtim_tur = cutim
    
  end if
  
  if (kfl_lmaxi_tur==1_ip)   call tur_maxlen()
  !-------------------------------------------------------------------
  !                                                 
  ! Interpolation from coarse to fine mesh      
  ! 
  !-------------------------------------------------------------------      

  if(kfl_meshi_tur /= 0_ip) call tur_coarfine(1_ip)

  if (INOTMASTER.and. kfl_rstar == 0) then
    
     if (kfl_logva==1) then
        do ipoin = 1,npoin
           untur(1:nturb_tur,ipoin,1) = log (untur(1:nturb_tur,ipoin,1))
        end do
        call tur_updunk(ITASK_INIUNK)
      
     else if (kfl_regularization) then
     
        do ipoin = 1, npoin 
           untur(1,ipoin,1) = inv_regul_k(untur(1,ipoin,1))
           untur(2,ipoin,1) = inv_regul_e(untur(2,ipoin,1))
        end do
        call tur_updunk(ITASK_INIUNK)
     
     end if
  end if
  
  
  !-------------------------------------------------------------------
  !
  ! Read forward values from fields for adjoint
  !
  !-------------------------------------------------------------------
  if( INOTMASTER .and. kfl_adj_prob == 1 ) then
    if( nfiel_tur(1) == 0 ) call runend('TUR_INIUNK: NO FIELDS WERE GIVEN FOR FORWARD VALUES')
    do ipoin=1,npoin
      untur_forw(1,ipoin,1) = xfiel(-nfiel_tur(1))%a(1,ipoin,1)
      untur_forw(1,ipoin,2) = xfiel(-nfiel_tur(1))%a(1,ipoin,1)
      if( kfl_rstar == 0 ) untur(:,ipoin,:) = 0.0_rp
    end do
  endif
  
end subroutine tur_iniunk
