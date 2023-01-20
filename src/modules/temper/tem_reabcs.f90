!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_reabcs()
  !-----------------------------------------------------------------------
  !****f* Temper/tem_reabcs
  ! NAME
  !    tem_reabcs
  ! DESCRIPTION
  !    This routine reads the boundary conditions for the temperature
  !    equation.
  !
  !    * For conditions on boundaries, bvnat(iboun)
  !
  !         Neumann ... fixbo(iboun)=2 -> bvnat(1,iboun)
  !         Robin ..... fixbo(iboun)=3 -> bvnat(3,iboun)
  !
  !    * For conditions on nodes, Neumann and Robin conditions are stored
  !      temporarily in bvess_tem(ipoin) and then changed to conditions on
  !      boundaries  in the routine tem_bcntoe (called at the end)
  !
  !    * Conditions on nodes have priority over conditions on boundaries.
  !      This is done and explained in tem_bcntoe. CHANGE?
  ! OUTPUT
  ! USES
  ! USED BY
  !    tem_turnon
  !***
  !> @addtogroup NastinInput
  !> @{
  !> @file    tem_reabcs.f90
  !> @author  Matias
  !> @brief   Read boundary conditions
  !> @details Read boundary conditions, initial conditions and parameters
  !> @}
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_temper
  use mod_opebcs
  use mod_ker_space_time_function
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip)  :: ifunc,ncodf,nbcod,dummi
  character(5) :: wfname
  !
  ! Allocate memory
  !
  if( kfl_icodn > 0  ) then
     call opnbcs(1_ip,1_ip,dummi,dummi,tncod_tem) ! Memory for structure
     call opnbcs(2_ip,1_ip, 1_ip, 0_ip,tncod_tem) ! Memory for variable

  end if
  npnat_tem = 4                                   ! # parameters natural bc
  if( kfl_icodb > 0  ) then
     call opbbcs(0_ip,1_ip,npnat_tem,tbcod_tem)
  end if
  call tem_membcs(4_ip)

  if( INOTSLAVE  ) then
     !
     ! Initializations
     !
     kfl_conbc_tem = 1             ! Constant boundary conditions
     kfl_inidi_tem = 0             ! If initial problem is pure diffusion = 1
     kfl_inico_tem = 0             ! No initial condition function
     kfl_intbc_tem = 0             ! Do not interpolate temperature
     delta_tem     = 0.0_rp        ! Distance to the wall
     initial_tem   = 0.0_rp        ! Initial constant temperature
     !
     ! Read flags
     !
     iknbo_tem     =  0            ! Boundary condition is known
     ncodf         =  1            ! Temper has 1 degree of freedom
     nbcod         = -1

     call ecoute('tem_reabcs')
     do while(words(1)/='BOUND')
        call ecoute('tem_reabcs')
     end do
     !--><group>
     !-->       <groupName>BOUNDARY_CONDITIONS</groupName>
     !-->       <groupType>bouco</groupType>
     !-->       <groupTypeDefinition><![CDATA[
     !-->       CODES, BOUNDARIES
     !-->         ...
     !-->         int1 int2 real1                                                              $  int1=code, int2=type of condition, real1=value
     !-->         ...
     !-->       END_CODES]]></groupTypeDefinition>
     !--></group>
     if(exists('WALLD') ) then
        delta_tem=getrea('WALLD',0.0_rp,'#Distance to the wall')
     end if
     if(exists('DIFFU') ) then
        kfl_inidi_tem=1
     end if
     if(exists('NONCO') ) then
        kfl_conbc_tem=0
     else
        kfl_conbc_tem=1
     end if
     if(exists('UNKNO') ) then
        iknbo_tem=1
     end if
     if(exists('TIMEI') ) then
        kfl_intbc_tem=2
     end if
     !
     ! Read data
     !
     call ecoute('tem_reabcs')
     do while( words(1) /= 'ENDBO' )

        if( words(1) == 'PARAM' ) then

           call ecoute('tem_reabcs')
           do while( words(1) /= 'ENDPA' )
              !
              ! Initial conditions
              !
              if( words(1) == 'INITI' ) then
                 if(      exists('CONST') ) then
                    kfl_inico_tem = 1
                    initial_tem = getrea('CONST',0.0_rp,'#Constant initial temperature')
                    if( exists('VALUE') ) then
                       initial_tem = getrea('VALUE',0.0_rp,'#Constant initial temperature')
                    end if
                 else if( exists('WARMB') ) then
                    kfl_inico_tem = 2
                 else if( exists('DIFFU') ) then
                    kfl_inico_tem = 3
                 else if( exists('SPACE') ) then
                    wfname        = getcha('SPACE','NULL ','#Space/time Function name')
                    kfl_inico_tem = -space_time_function_number(wfname)
                 else if( exists('FIELD') ) then
                    kfl_inico_tem = -100-getint('FIELD',1_ip,'#FIELD NUMBER')
                 end if
              end if
              call ecoute('tem_reabcs')
           end do

        else if( words(1) == 'CODES'.and.exists('NODES') ) then
           !
           ! User-defined codes on nodes
           !
           tncod => tncod_tem
           call boundary_conditions_read_node_codes('TEMPERATURE')

        else if( words(1) == 'CODES'.and.exists('BOUND') ) then
           !
           ! User-defined codes on boundaries
           !
           kfl_fixbo => kfl_fixbo_tem
           bvnat     => bvnat_tem(:,:,1)
           tbcod     => tbcod_tem(1:)
           call boundary_conditions_read_boundary_codes('TEMPERATURE')

        else if( words(1) == 'CODES' ) then
           call runend('TEM_REABCS: SPECIFY IF CODES ARE APPLIED TO NODES OR BOUNDARIES')

        else if( words(1) == 'INITI' ) then
           !
           ! Initial condition
           !
           if( words(2) == 'CONST' ) then
              kfl_inico_tem = 1
              ! GGU: It seems that this is for retrocompatability but initial_tem is not defined!
           else if( words(2) == 'WARMB' ) then
              kfl_inico_tem = 2
           else if( words(2) == 'SPACE' ) then
              wfname        = getcha('SPACE','NULL ','#Space/time Function name')
              kfl_inico_tem = -space_time_function_number(wfname)
           else if( words(2) == 'VALUE' ) then
              kfl_inico_tem = -100 &
                         -getint('VALUE',1_ip,'#Initial condition is from value function')
           end if

        else if( words(1) == 'FUNCT' ) then
           !
           ! Functions
           !
           call ecoute('tem_reabcs')
           do while(words(1)/='ENDFU')
              if(kfl_conbc_tem==0 ) then
                 ifunc=getint('FUNCT',1_ip,'#FUNCTION NUMBER')
                 if(ifunc<0.or.ifunc>10 ) then
                    call runend('tem_reabcs: WRONG FUNCION NUMBER')
                 else
                    if( words(2) == 'PARAB' ) then
                       kfl_funty_tem(ifunc)=1
                    else if( words(2) == 'PERIO' ) then
                       kfl_funty_tem(ifunc)=2
                    else if( words(2) == 'MAREK' ) then
                       kfl_funty_tem(ifunc)=5
                    end if
                    if(kfl_funty_tem(ifunc)>0) funpa_tem(1:6,ifunc)=param(3:8)
                 end if
              end if
              call ecoute('tem_reabcs')
           end do

        end if

        call ecoute('tem_reabcs')

     end do

  end if

end subroutine tem_reabcs
