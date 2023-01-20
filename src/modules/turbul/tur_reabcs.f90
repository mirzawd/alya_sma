!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_reabcs()
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_reabcs
  ! NAME 
  !    tur_reabcs
  ! DESCRIPTION
  !    This routine reads TURBUL boundary conditions.
  !
  !    The different codes for KFL_FIXNO_TUR(1,IPOIN) are:
  !    = 0 ... Free or initial
  !    = 1 ... Dirichlet: fixed value
  !    = 2 ... Nothing
  !    = 3 ... Wall law 
  !    = 4 ... Wall
  !    = 5 ... Presctibe eps knowing the turbulent length scale
  !    = 6 ... Inflow condition
  !    = 7 ... Adaptive depending on the angle between normal to the boundary and velocity. 
  !    = 8 ... Adaptive inflow condition
  ! USES
  ! USED BY
  !    tur_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_turbul
  use mod_opebcs
  use mod_ker_space_time_function
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip)  :: dummi,iturb
  character(5) :: wfname
  !
  ! Allocate memory
  !
  if( kfl_icodn > 0 ) then
     call opnbcs(1_ip,nturb_tur,dummi,dummi,tncod_tur)       ! Memory for structure
     do iturb = 1,nturb_tur
        call opnbcs(2_ip,iturb,1_ip, 0_ip,tncod_tur)         ! Memory for variable
     end do
  end if
   
  if( kfl_icodb > 0 ) then
     call opebcs_initialization_structure(nturb_tur,tbcod_tur)     ! Allocates for nturb
     do iturb=1, nturb_tur
        call opebcs_initialization_variable (npnat_tur,tbcod_tur(iturb))
     end do
!     call opbbcs(0_ip,nturb_tur,1_ip,tbcod_tur)      
  end if
  if( kfl_geome > 0 ) then
     call opnbcs(0_ip,nturb_tur,1_ip, 0_ip,tgcod_tur)
  end if

  if( INOTSLAVE ) then
     !
     ! Dimensions
     !
     kfl_inidi_tur =  0         ! No initial diffusion problem
     kfl_wallw_tur =  0         ! Classical b.c.
     kfl_infl1_tur =  0         ! Inflow type for 1st turbulence variable (cst)
     kfl_infl2_tur =  0         ! Inflow type for 2nd turbulence variable (cst)
     kfl_usrbc_tur =  0         ! User boundary condition
     kfl_initi_tur =  0         ! Initial solution
     kfl_valbc_tur =  0         ! Value function for initial condition
     delta_tur     =  0.0_rp    ! Delta for the law of the wall
     turin_tur     = -0.01_rp   ! Intensity of turbulence
     turle_tur     =  1.0_rp    ! Length scale of turbulence (for eps b.c.)
     nutnu_tur     = -1.0_rp    ! Ratio nut/nu
     hdiam_tur     =  0.0_rp    ! Hydraulic diameter
     rebcs_tur     =  1.0_rp    ! Realxation of boundary conditions
     xinit_tur     =  0.0_rp    ! Initial value
     !
     ! Reach the boundary condition section
     !
     call ecoute('tur_reabcs')
     do while(words(1)/='BOUND')
        call ecoute('tur_reabcs')
     end do
     if(exists('WALLD')) then
        delta_tur=getrea('WALLD',0.0_rp,'#Distance to the wall')
     end if
     if(exists('INTEN')) then
        turin_tur=getrea('INTEN',0.01_rp,'#Intensity of turbulence')
     end if
     if(exists('LENGT')) then
        turle_tur=getrea('LENGT',1.0_rp,'#Length scale of turbulence')
     end if
     if(exists('HYDRA')) then
        hdiam_tur=getrea('HYDRA',0.0_rp,'#Hydraulic diameter')
     end if
     if(exists('DIFFU')) then
        kfl_inidi_tur=1
     else
        kfl_inidi_tur=0
     end if
     if(exists('OMEGA')) then
        if(exists('BREDB')) kfl_wallw_tur=1
     end if
     !
     ! Read data
     !
     call ecoute('tur_reabcs')
     do while(words(1)/='ENDBO')

        if(words(1)=='CODES'.and.exists('NODES')) then 

           if(exists('GEOME')) then
              !
              ! Geometrical b.c.
              !                        
              if( exists('VARIA') ) then
                 iturb = getint('VARIA',1_ip,'*USER BOUNDARY CONDITIONS')
                 tgcod => tgcod_tur(iturb:)
                 call reacod(4_ip)
              else
                 tgcod => tgcod_tur(1:)
                 call reacod(4_ip)
                 do iturb = 2,nturb_tur
                    call cpybcs(1_ip,iturb,tgcod_tur) 
                 end do
              end if

           else
              !
              ! User-defined codes on nodes
              !
              if( exists('VARIA') ) then
                 iturb = getint('VARIA',1_ip,'*USER BOUNDARY CONDITIONS')
                 tncod => tncod_tur(iturb:)
                 call reacod(1_ip)
              else
                 tncod => tncod_tur(1:)
                 call reacod(1_ip)
                 do iturb = 2,nturb_tur
                    call cpybcs(1_ip,iturb,tncod_tur) 
                 end do
              end if

           end if

        else if(words(1)=='CODES'.and.exists('BOUND')) then 
           !
           !  User defined codes on boundaries
           !
           if( exists('VARIA') ) then
              iturb = getint('VARIA',1_ip,'*USER BOUNDARY CONDITIONS')
              kfl_fixbo => kfl_fixbo_tur
              bvnat     => bvnat_tur(:,:,iturb)
              tbcod     => tbcod_tur(iturb:)
              call boundary_conditions_read_boundary_codes('TURBULENCE')
           else
              kfl_fixbo => kfl_fixbo_tur
              bvnat     => bvnat_tur(:,:,1)
              tbcod     => tbcod_tur(1:)
              call boundary_conditions_read_boundary_codes('TURBULENCE')     
              do iturb =2, nturb_tur
                 call cpybcs_boundaries(1_ip,iturb,tbcod_tur) 
              end do
           end if

        else if(words(1)=='USEUS') then
           !
           ! Use user boundary conditions 
           !
           kfl_usrbc_tur = getint('USEUS',1_ip,'*USER BOUNDARY CONDITIONS')

        else if(words(1)=='INITI') then
           !
           ! Initial
           !
           if(exists('BVESS')) then   ! uses the values in bvess_tur for the initial condition
              kfl_initi_tur = 2 
           else if( exists('VALUE') ) then    
              !              
              ! Value function
              !
              kfl_initi_tur = -1
              if( .not. exists('VARIA') ) then
                 call runend('TUR_REABCS: TURBULENCE VARIABLE MUST BE SPECIFIED')
              end if
              iturb                =  getint('VARIA',1_ip,'#Turbulence variable')
              if( iturb < 1 .or. iturb > 4 ) then
                 call runend('TUR_REABCS: WRONG TURBULENCE VARIABLE')
              end if
              kfl_valbc_tur(iturb) = getint('VALUE',1_ip,'#Initial condition is from value function') 
           else if( exists('SPACE') ) then                                           ! >100: Space time function
              wfname           = getcha('SPACE','NULL ','#Space/time Function name')
              iturb            = getint('VARIA',1_ip,'#Turbulence variable')
              kfl_initi_tur    = 3 
              xinit_tur(iturb) = real(space_time_function_number(wfname),rp)              
           else         
              kfl_initi_tur = 1
              xinit_tur(1:nturb_tur) = param(1:nturb_tur)
           end if

        else if(words(1)=='PARAM') then
           !
           ! Parameters
           !
           call ecoute('tur_reabcs')
           do while(words(1)/='ENDPA')

              if(words(1)=='INLET'.or.words(1)=='INFLO') then ! Inflow (codes 6 and 8)
                 ! 1st variable
                 if(exists('CONST')) kfl_infl1_tur = 0
                 if(exists('INTEN')) kfl_infl1_tur = 1
                 if(exists('CHANN')) kfl_infl1_tur = 2
                 if(exists('SMOOT')) kfl_infl1_tur = 3
                 ! 2nd variable
                 if(exists('CONST')) kfl_infl2_tur = 0
                 if(exists('LENGT')) kfl_infl2_tur = 1
                 if(exists('CHANN')) kfl_infl2_tur = 2
                 if(exists('RATIO')) kfl_infl2_tur = 3
                 ! ABL model (atmospheric boundary layer)
                 if(exists('ABL  ')) kfl_infl1_tur = 4
                 if(exists('ABL  ')) kfl_infl2_tur = 4

              else if(words(1)=='INTEN') then
                 turin_tur=getrea('INTEN',0.01_rp,'#Intensity of turbulence')

              else if(words(1)=='RATIO') then
                 nutnu_tur=getrea('RATIO',1.0_rp,'#nut over nu')

              else if(words(1)=='VISCO') then
                 nutnu_tur=getrea('VISCO',1.0_rp,'#nut over nu')

              else if(words(1)=='LENGT') then
                 turle_tur=getrea('LENGT',1.0_rp,'#Length scale of turbulence')

              else if(words(1)=='HYDRA') then
                 hdiam_tur=getrea('HYDRA',0.0_rp,'#Hydraulic diameter')

              else if(words(1)=='OMEGA') then
                 if(words(2)=='BREDB') kfl_wallw_tur=1

              else if(words(1)=='RELAX') then
                 rebcs_tur=getrea('RELAX',1.0_rp,'#Relaxation of boundary conditions')

              end if

              call ecoute('tur_reabcs')
           end do
        end if

        call ecoute('tur_reabcs')
     end do

  end if

end subroutine tur_reabcs
