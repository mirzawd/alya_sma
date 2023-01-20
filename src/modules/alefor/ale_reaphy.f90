!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_reaphy()
  !------------------------------------------------------------------------
  !****f* Temper/ale_reaphy
  ! NAME 
  !    ale_reaphy
  ! DESCRIPTION
  !    This routine reads the physical problem definition for the
  !    temperature equation.
  ! USES
  ! USED BY
  !    ale_turnon
  !------------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_alefor
  use def_domain
  use mod_memory,  only :  memory_alloca
  use mod_ecoute,  only :  ecoute
  use mod_alefor
  implicit none
  integer(ip) :: ipara,iimbo

  if( INOTSLAVE ) then
     !
     ! Initializations (defaults)
     !
     kfl_rigid_ale =  0             ! no RB 
     nstro_ale =  0                 ! RB start rotation from step 0
     nstli_ale =  0                 ! RB start linear motion from step 0
     kfl_mvext_ale =  0             ! RB do not move exterior
     kfl_ralei_ale =  0             ! RB no Raleigh damping
     nstra_ale =  1000000           ! RB start eliminating Raleigh damping at a very high step
     nenra_ale =  2000000           ! RB end eliminating Raleigh damping at a very high step

     kfl_grafo_ale =  0_ip          ! Gravity force
     kfl_catfo_ale =  0_ip          ! Catamaran force
     kfl_sprin_ale =  0_ip          ! Linear Spring forces kx-nv
     kfl_topmo_ale =  0_ip          ! Non-spinning top model
     kfl_pertu_ale =  0_ip          ! Use perturbation for gneralized forces 
     kfl_genco_ale =  0_ip          ! Solve the rigid body equations in generalized coordinates
     xline_ale     =  0.0_rp        ! Linear motion - now no longer read but calculated fron nstli_ibm
     xrota_ale     =  0.0_rp        ! Rotation motion - now no longer read but calculated fron nstro_ibm

     iimbo = 1_ip   ! for the moment only one rigid body
     !
     ! Reach the section
     !
     call ecoute('ale_reaphy')
     do while( words(1) /= 'PHYSI' )
        call ecoute('ale_reaphy')
     end do
     !
     ! Begin to read data
     !
     do while( words(1) /= 'ENDPH' )
        call ecoute('ale_reaphy')


        if( words(1) == 'RIGID' ) then
           kfl_rigid_ale = 1
           nrbod = 1   ! for the moment we only allow one rigid body
                       ! then I should create nrbse and rbset for each rigid body

           call alefor_memory_allocate('RBBOU')

           call ecoute('ale_reaphy')
           do while(words(1)/='ENDRI')
! esto es lo que iba en dimensions
              if( words(1) == 'SETS ' ) then
                 rbbou(iimbo) % nrbse = 0
                 do ipara = 1,nnpar
                    rbbou(iimbo) % lrbse(ipara) = nint(param(ipara))
                    if ( rbbou(iimbo) % lrbse(ipara) /= 0 ) rbbou(iimbo) % nrbse = rbbou(iimbo) % nrbse + 1
                 end do
              else if( words(1) == 'MASS ' ) then
                 rbbou(iimbo) % massa = param(1)
                 if( param(1) <= 0.0_rp ) call runend('RIGID_BODY: WRONG PARTICLE MASS')

              else if( words(1) == 'DENSI' ) then
                 rbbou(iimbo) % densi = param(1)
                 if( param(1) <= 0.0_rp ) call runend('RIGID_BODY: WRONG PARTICLE DENSITY')

              else if( words(1) == 'VOLUM' ) then
                 rbbou(iimbo) % volum = param(1)

              else if( words(1) == 'MOMEN' ) then
                 if( param(1) < 0.0_rp ) call runend('RIGID_BODY: WRONG PARTICLE MOMENTUM OF INERTIA')
                 rbbou(iimbo) % momin      = 0.0_rp  ! initialize to zero in case not all are defined
                 do ipara = 1,max(6_ip,nnpar)
                    rbbou(iimbo) % momin(ipara) = param(ipara)
                 end do

              else if( words(1) == 'POSGR' ) then
                 if( param(1) < -0.5e12_rp ) call runend('RIGID_BODY: WRONG PARTICLE CENTER OF GRAVITY')
                 rbbou(iimbo) % posgr      = 0.0_rp  ! initialize to zero in case not all are defined
                 do ipara = 1,max(3_ip,nnpar)
                    rbbou(iimbo) % posgr(ipara) = param(ipara)
                 end do
                 !
                 ! Generalized coordinates
                 !
              else if( words(1) == 'GENER' ) then
                 kfl_genco_ale = 1_ip
                 !
                 ! Non spinning top model (developed for vortex-bladeless)
                 !
                 if( words(2) == 'TOPMO' ) then
                    if( exists('PERTU') ) kfl_pertu_ale = 1_ip
                    kfl_topmo_ale = 1_ip
                    sprin_ale     = 0_rp
                    !
                    ! Stiffness at the fixed point (hookean law)
                    !
                    sprin_ale(1_ip) = param(2)
                    !
                    ! damping at the fixed point (hookean law)
                    !
                    sprin_ale(2_ip) = param(3)
                    !
                    ! distance between the center of mass and the fixed point (hookean law)
                    !
                    sprin_ale(3_ip) = param(4)


                 else if( words(2) == 'SPRIN' ) then
                    kfl_sprin_ale   = 1
                    sprin_ale       = 0_rp
                    if( words(3) == 'RESID' )then
                       kfl_forca_res = 1_ip
                       ! k
                       sprin_ale(1_ip) = param(3)
                       ! damping
                       sprin_ale(2_ip) = param(4)
                       ! force scaling (for non-dimensional variables)
                       sprin_ale(3_ip) = param(5)
                    else
                       kfl_forca_res = 0_ip
                       ! k
                       sprin_ale(1_ip) = param(2)
                       ! damping
                       sprin_ale(2_ip) = param(3)
                       ! force scaling (for non-dimensional variables)
                       sprin_ale(3_ip) = param(4)                       
                    end if
                 end if
! estos ver si realmente valen la pena
              else if( words(1) == 'POSIL' ) then
                 do ipara = 1,3
                    rbbou(iimbo) % posil(ipara,1) = param(ipara)
                 end do
              else if( words(1) == 'VELOL' ) then
                 do ipara = 1,3
                    rbbou(iimbo) % velol(ipara,1) = param(ipara)
                 end do
              else if( words(1) == 'ACCEL' ) then
                 do ipara = 1,3
                    rbbou(iimbo) % accel(ipara,1) = param(ipara)
                 end do
              else if( words(1) == 'POSIA' ) then
                 do ipara = 1,3
                    rbbou(iimbo) % posia(ipara,1) = param(ipara)
                 end do
              else if( words(1) == 'VELOA' ) then
                 do ipara = 1,3
                    rbbou(iimbo) % veloa(ipara,1) = param(ipara)
                 end do
              else if( words(1) == 'ACCEA' ) then
                 do ipara = 1,3
                    rbbou(iimbo) % accea(ipara,1) = param(ipara)
                 end do
!
              else if(words(1)=='ROTAT' ) then
                 if(words(2)=='START')   nstro_ale(1:3) = int(param(2:4),ip)

              else if(words(1)=='TRANS' ) then
                 if(words(2)=='START')   nstli_ale(1:3) = int(param(2:4),ip)

              else if(words(1)=='MOVEE' ) then
                 kfl_mvext_ale = 1

              else if(words(1)=='RALEI' ) then
                 kfl_ralei_ale = 1
                 ralei_ale     = 1.0_rp
                 if(words(2)=='VALUE')   ralei_ale = param(2)
                 if(words(3)=='START')   nstra_ale = int(param(3),ip)
                 if(words(4)=='END  ')   nenra_ale = int(param(4),ip)


              else if( words(1) == 'FORCE' ) then
                 !
                 ! Forces
                 !
                 if( exists('CATAM') ) kfl_catfo_ale = 1
                 if( exists('GRAVI') ) kfl_grafo_ale = 1

              end if

              call ecoute('ale_reaphy')
           end do    ! END_RIGID_BODY
        else if(words(1)=='SENSI') then
           if(words(2)=='FIELD') then
              kfl_sensi_ale = - getint('FIELD',-1_ip,'#Field number for the mesh sensitivities')
           end if
        else if (words(1)=='INITI') then              
           if(words(2)=='DISPL') then
              moddi_ale = -getint('FIELD',1_ip,'#Field Number for displacements')
           else if(words(2)=='VELOC') then
              modvi_ale = -getint('FIELD',1_ip,'#Field Number for velocities')              
           end if
        end if  
     end do

  end if

end subroutine ale_reaphy

