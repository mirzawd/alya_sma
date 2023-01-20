!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_outinf(itask)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_outinf
  ! NAME 
  !    nsi_outinf
  ! DESCRIPTION
  !    This routine writes on the incompressible Navier-Stokes files
  ! USED BY
  !    nsi_turnon
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_solver
  use def_nastin
  use def_kermod, only : turmu_ker, poros_ker
!  use def_kermod, only : gasco
  use mod_outfor, only : outfor
  use mod_iofile, only : iofile_flush_unit
  implicit none
  integer(ip), intent(in) :: itask
  character(100)          :: equat
  character(100)          :: equa2
  character(100)          :: equa3
  character(10)           :: lvisc
  integer(ip)             :: ierhs
!  character(10)           :: Auu,Aup,Apu,App
!  character(20)           :: bu,bp
 
  if( INOTSLAVE ) then

     select case(itask)

     case(1)
        !
        ! Write information in Result file
        !
        if(kfl_rstar/=2) then   
           equat=''
           if(kfl_timei_nsi==1) equat=trim(equat)//'rho*du/dt'
           if(kfl_advec_nsi>=1) then
              !
              ! Convection term
              !
              if(     kfl_convection_type_nsi == NSI_CONVECTION_SKEW ) then                 
                 equat=trim(equat)//' + u.grad(u) + 1/2 (div u) u'
              else if( kfl_convection_type_nsi == NSI_CONVECTION_CONSERVATIVE ) then
                 equat=trim(equat)//' + u.grad(u) + (div u) u'
              else if( kfl_convection_type_nsi == NSI_CONVECTION_EMAC ) then
                 equat=trim(equat)//' + u.grad(u) + u.grad(u)^t + (div u) u'
              else
                 equat=trim(equat)//' + rho*(u.grad)u'
              end if
           end if
           if(kfl_visco_nsi/=0) then
              !
              ! Diffusion term
              !
              lvisc='mu*'
              if( turmu_ker % kfl_exist /= 0_ip ) lvisc='(mu+mu_t)*'
              if(int(fvins_nsi)==0) then
                 equat=trim(equat)//' - div['//trim(lvisc)//'grad(u)] '
              else if(int(fvins_nsi)==1) then
                 equat=trim(equat)//' - div[2*'//trim(lvisc)//'*eps(u)] '
              else
                 equat=trim(equat)//' - div[2*'//trim(lvisc)//'*(eps(u)-1/3*div(u) I)] '
              end if
           end if
           
           if(poros_ker % kfl_exist==1) equat=trim(equat)//'+ rho* drag* u'
           
           if(fvnoa_nsi>zensi) then
              !
              ! Rotation term
              !
              equat=trim(equat)//' + 2*rho*w x u '
           end if

           equat=trim(equat)//' + grad(p)'
           equat=trim(equat)//' = '
           ierhs=0
           if(grnor_nsi>zensi)  then
              ierhs=1
              equat=trim(equat)//' + rho*g'
           end if
           if(fvnoa_nsi>zensi) then
              equat=trim(equat)//' - rho*w x w x (r-r0) -rho*dw/dt x (r-r0) '
           end if
           if(fvnol_nsi>zensi) then
              equat=trim(equat)//' - rho*a' 
           end if
           if(kfl_cotem_nsi/=0) then
              ierhs=1
              equat=trim(equat)//' - rho*g*beta*(T-Tref)'
           end if
           if(ierhs==0) equat=trim(equat)//'0 '

           call outfor(25_ip,momod(modul) % lun_outpu,'DIFFERENTIAL EQUATION')
           
           write(momod(modul) % lun_outpu,110) trim(equat)
           if(kfl_regim_nsi==1) then
              write(momod(modul) % lun_outpu,113) '1/(RT)*dp/dt + 1/(RT)*u.grad(p) + rho*div(u) = rho/T*[ dT/dt + u.grad(T) ]'
           else if(kfl_regim_nsi==3) then
              write(momod(modul) % lun_outpu,113) 'rho*div(u) = rho/T*[ dT/dt + u.grad(T) ] - rho/p0*dp0/dt'
              if(kfl_prthe_nsi==1) then
                 write(momod(modul) % lun_outpu,113) 'p0(t) = p0(0)*[int_V 1/T(0) dV]/[int_V 1/T(t) dV]'
              else
                 write(momod(modul) % lun_outpu,113) &
                      'V/gam*dp0/dt + p0*[int u.n ds] = (gam-1)/gam*[ (int_S q.n ds) + (int_V Q dV) ]'
              end if
           else
              write(momod(modul) % lun_outpu,113) 'div(u)=0'
           end if
           !write(momod(modul) % lun_outpu,112)
           !write(momod(modul) % lun_outpu,111) 'rho [ M / L^3 ]     = ',densi_nsi(1,1)
           !write(momod(modul) % lun_outpu,111) 'mu  [ M / L T ]     = ',visco_nsi(1,1)
           !if(kfl_regim_nsi>=1) then
           !   write(momod(modul) % lun_outpu,111) 'R   [ L^2 / T^2.K ] = ',gasco
           !   write(momod(modul) % lun_outpu,111) 'Cp  [ L^2 / T^2 K ] = ',sphea_nsi
           !   write(momod(modul) % lun_outpu,111) 'gam=Cp/(Cp-R)       = ',gamth_nsi
           !end if
           write(momod(modul) % lun_outpu,*)
           !
           ! Others
           !
           call outfor(25_ip,momod(modul) % lun_outpu,'NUMERICAL PARAMETERS')

           coutp = 'NOT DEFINED'
           if(      NSI_MONOLITHIC ) then
              coutp(1) = 'MONOLITHC'
           else if( NSI_SCHUR_COMPLEMENT ) then
              coutp(1) = 'SCHUR COMPLEMENT'
           else
              coutp(1) = 'FRACTIONAL STEP'
           end if           
           if(       kfl_matdi_nsi == NSI_DIRICHLET_MATRIX ) then
              coutp(2) = 'ON MATRIX'
           else if( kfl_matdi_nsi == NSI_DIRICHLET_ELEMENT ) then
              coutp(2) = 'ON ELEMENT'
           else if( kfl_matdi_nsi == NSI_DIRICHLET_ALGORITHM ) then
              coutp(2) = 'ON ALGORITHM'
           end if
           if( kfl_tisch_nsi == 1 ) then
              coutp(3) = 'TRAPEZOIDAL'
           else if( kfl_tisch_nsi == 2 ) then
              coutp(3) = 'BDF'
           else if( kfl_tisch_nsi == 3 ) then
              coutp(3) = 'ADAMS-BASHFORTH'
           else if( kfl_tisch_nsi == 4 ) then
              coutp(3) = 'RUNGE-KUTTA'
           end if
           ioutp(1) = kfl_tiacc_nsi
           if( kfl_stabi_nsi == NSI_GALERKIN ) then
              if( NSI_FRACTIONAL_STEP ) then
                 coutp(4) = 'GALERKIN FOR VELOCITY, ALGEBRAIC SPLIT OSS FOR PRESSURE'
              else
                 coutp(4) = 'GALERKIN'
              end if
           else if( kfl_stabi_nsi == NSI_ALGEBRAIC_SPLIT_OSS ) then
              coutp(4) = 'GALERKIN FOR VELOCITY, ALGEBRAIC SPLIT OSS FOR PRESSURE'
           else if( kfl_stabi_nsi == NSI_ASGS ) then
              coutp(4) = 'ASGS'
           else if( kfl_stabi_nsi == NSI_OSS ) then
              coutp(4) = 'FULL OSS'
           else if( kfl_stabi_nsi == NSI_SPLIT_OSS ) then
              coutp(4) = 'SPLIT OSS'
           end if
           if(      kfl_sgsco_nsi == 1 .and. kfl_sgsti_nsi == 0 ) then       
              coutp(5) = 'CONVECTION'
           else if( kfl_sgsco_nsi == 1 .and. kfl_sgsti_nsi == 1 ) then       
              coutp(5) = 'CONVECTION AND TIME'
           else if( kfl_sgsco_nsi == 0 .and. kfl_sgsti_nsi == 1 ) then       
              coutp(5) = 'TIME'
           else
              coutp(5) = 'NONE'
           end if           
           write(momod(modul) % lun_outpu,106) &
                trim(coutp(1)),trim(coutp(2)),&
                trim(coutp(3)),ioutp(1),&
                trim(coutp(4)),trim(coutp(5))
           !
           ! Specific parameters
           !
           if( NSI_FRACTIONAL_STEP ) then
              if( gamma_nsi == 1.0_rp ) then
                 coutp(3) = 'INCREMENTAL'
              else
                 coutp(3) = 'NON-INCREMENTAL'
              end if
              if( kfl_press_stab_nsi == 0 ) then
                 coutp(5) = 'Q = -Tim L'
              else
                 coutp(5) = 'Q = -Tau L'
              end if
              if( kfl_grad_div_nsi == 1 ) then
                 coutp(1) = 'CONSTANT = G'
                 coutp(2) = 'CONSTANT = D=G^t'
                 coutp(4) = 'CONSTANT = L'
              else
                 coutp(1) = 'RECOMPUTED'
                 coutp(2) = 'RECOMPUTED'                 
                 coutp(4) = 'RECOMPUTED'                 
              end if
              if( kfl_press_stab_nsi == 0 ) then
                 if( kfl_grad_div_nsi == 1 ) then
                    equat = 'S   =  D Tim M^{-1} G - Tim L'                    
                 else
                    equat = 'S   =  Apu Tim M^{-1} Aup - Tim L'
                 end if
                 equa2 = 'Tim := dt/rho'
                 equa3 = 'L   := - int_V grad(p).grad(q) dV'
              else
                 if( kfl_grad_div_nsi == 1 ) then
                    equat = 'S   =  Apu Tau M^{-1} Aup - Tau L'
                 else
                    equat = 'S   =  D Tau M^{-1} G - Tau L' 
                 end if
                 equa2 = 'Tau := tau/rho'
                 equa3 = 'L   := - int_V grad(p).grad(q) dV'               
              end if
              write(momod(modul) % lun_outpu,108) &
                   trim(coutp(3)),gamma_nsi,&
                   trim(coutp(1)),&
                   trim(coutp(2)),&
                   trim(coutp(4)),trim(coutp(5)),&
                   trim(equat),&
                   trim(equa2),&
                   trim(equa3)
           end if

!!$           if( NSI_FRACTION_STEP ) then
!!$              Auu = 'Tim^{-1} M'
!!$              if( kfl_grad_div_nsi /= 0 ) then
!!$                 bu  = 'bu - Auu u^n + Tim^{-1} M u^n - gamma Aup p^n'
!!$              else
!!$                 bu  = 'bu - Auu u^n + Tim^{-1} M u^n - gamma G p^n'
!!$              end if
!!$           else
!!$              Auu = 'Auu'
!!$              bu  = 'bu'
!!$           end if
!!$           if( kfl_grad_div_nsi /= 0 ) then
!!$              Aup = 'G'
!!$              Apu = 'D'
!!$           else
!!$              Aup = 'Aup'
!!$              Apu = 'Apu'
!!$           end if
!!$
!!$           if( kfl_stabi_nsi == NSI_ALGEBRAIC_SPLIT_OSS ) then
!!$              App = 'S'
!!$              bp  = 'bp + gamma * S * p^n'
!!$           else if( kfl_stabi_nsi == NSI_ASGS ) then
!!$              if( kfl_sgsti_nsi == 1 ) then
!!$                 App = 'S1'
!!$                 bp  = 'bp'
!!$              else
!!$                 App = 'S2'
!!$                 bp  = 'bp'
!!$              end if
!!$           else if( kfl_stabi_nsi == NSI_OSS ) then
!!$              App = 'S1'
!!$              bp  = 'bp - P1 p^{i-1}'
!!$           else if( kfl_stabi_nsi == NSI_SPLIT_OSS ) then
!!$              App = ''
!!$           end if
           call iofile_flush_unit(momod(modul) % lun_outpu)
        end if

     case(2)

     end select
  end if
  !
  ! Formats
  !
110 format(/,10x,a)
112 format(//,&
       & 10x,'Physical properties: ',/,&
       & 10x,'------------------------------------------------------------ ')
111 format(&
       & 10x,a,10(1x,e12.6))
106 format(&
         & 5x,'   SOLUTION METHOD:         ',a,/,&
         & 5x,'   DIRICHLET CONDITIONS:    ',a,/,&
         & 5x,'   TIME INTEGRATION SCHEME: ',a,', ORDER=',i2,/,&
         & 5x,'   STABILIZATION METHOD:    ',a,/,&
         & 5x,'   SUBGRID SCALE TRACKING:  ',a)         
108 format(/,&
         & 5x,'   FRACTIONAL STEP OPTIONS:',/,&
         & 5x,'   ------------------------',//,&
         & 5x,'   PRESSURE SCHEME=        ',a,', GAMMA= ',f3.1,/,&
         & 5x,'   MATRIX Aup=             ',a,/,&
         & 5x,'   MATRIX Apu=             ',a,/,&
         & 5x,'   MATRIX L=               ',a,', ',a,/,&
         & 5x,'   PRESSURE STABILIZATION: ',a,/,&
         & 5x,'   SCALING FACTOR:         ',a,/,&
         & 5x,'   LAPLACIAN MATRIX:       ',a)
107 format(&
         & 5x,'+-          -+ +-        -+   +-  -+',/,&
         & 5x,'! Auu  Aup   | | u^{n+1}  |   | bu |',/,&
         & 5x,'!            | |          | = |    |',/,&
         & 5x,'! Apu  App   | | p^{n+1}  |   | bp |',/,&
         & 5x,'+-          -+ +-        -+   +-  -+')

113 format(&
       & 10x,a)
104 format(&
       '#',/,&
       '# CONTINUITY EQ:',/,&
       '# --------------')
105 format(&
       '#',/,&
       '# MOMENTUM EQS:',/,&
       '# -------------')

end subroutine nsi_outinf
      
