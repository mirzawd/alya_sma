!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_reabcs()
  !-----------------------------------------------------------------------
  !****f* wavequ/lev_reabcs
  ! NAME 
  !    lev_reabcs
  ! DESCRIPTION
  !    This routine reads the boundary conditions
  ! USES
  !    lev_membc
  ! USED BY
  !    lev_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_levels
  use mod_opebcs
  use mod_ker_space_time_function
  use mod_ecoute, only : ecoute
  use mod_messages, only : livinf
  implicit none
  integer(ip)  :: dummi
  integer(ip)  :: icobd, ivcod, ivcob
!  integer(ip)  :: icond, idime, ipoin, kdime
  integer(ip)  :: iread
  character(5) :: wfname

  if( kfl_icodn > 0 ) then
     call opnbcs(1_ip,1_ip,dummi,dummi,tncod_lev) ! Memory for structure
     call opnbcs(2_ip,1_ip, 1_ip, 0_ip,tncod_lev) ! Memory for variable
  end if
  if( kfl_geome > 0 ) then
     call opnbcs(0_ip,1_ip,1_ip,0_ip,tgcod_lev)
  end if

  if( INOTSLAVE ) then 
     !
     ! Initialization
     !
     kfl_inlev_lev = 0_ip
     kfl_conbc_lev = 1_ip
     height_lev    = 0.0_rp

     !call livinf(0_ip,'READ BOUNDARY CONDITIONS',0_ip)
     !
     ! Initializations and allocations
     !
     icobd     = -99999
     !     kfl_geome = 0             ! Do not compute geometrical normals -- esto esta mal lo metio bsc21235 ??? jode mod_opebcs L472!!!
     kfl_convx = 1             ! Use exnor for convex angle nodes
     kfl_frees = 0             ! Freestream criteria. 0: use wind angle to decide inflow/outflow
     npbcs     = 0             ! Max # of possibe geometrical conditions
     lsbcs     = 0             ! List of boundaries for a given geometrical code
     awind     =  0.0_rp       ! Wind angle for freestream
     tolan     =  0.0_rp       ! Tolerance used to define inflow from freestream 
     geoan     = 45.0_rp       ! Angle for geometrical normals

     ivcod = 0
     ivcob = 0
     iread = 0

     !
     ! Reach the nodal-wise section
     !
     call ecoute('lev_reabcs')
     do while(words(1)/='BOUND')
        call ecoute('lev_reabcs')
     end do
     !
     ! Loop over nodes and or boundaries
     !
     do while(words(1)/='ENDBO')
        call ecoute('lev_reabcs')

        if( words(1) == 'PARAM' ) then

           !-------------------------------------------------------------
           !
           ! Parameters
           !
           !-------------------------------------------------------------

           call ecoute('nsi_reabcs')
           do while(words(1)/='ENDPA') 
              if(words(1)=='INITI') then                 
                 !
                 ! Initialization of level set  
                 !
                 if(exists('HEIGH')) then                 
                    kfl_inlev_lev = -1
                    if( exists('ONX  ')) kfl_inlev_lev = -3
                    if( exists('ONY  ')) kfl_inlev_lev = -2
                    if( exists('ONZ  ')) kfl_inlev_lev = -1
                    if( exists('ONXIN')) kfl_inlev_lev = -6
                    if( exists('ONYIN')) kfl_inlev_lev = -5
                    if( exists('ONZIN')) kfl_inlev_lev = -4
                    height_lev    = getrea('HEIGH',0.0_rp,  '#Height of the level set')
                 else if(exists('FUNCT')) then  
                    kfl_inlev_lev = getint('FUNCT',1_ip,  '#Initialization function')
                 else if( exists('SPACE') ) then                                           ! >100: Space time function
                    wfname = getcha('SPACE','NULL ','#Space/time Function name')
                    kfl_initi_lev = 100 + space_time_function_number(wfname)
                 else if( exists('VALUE') ) then                                           ! Guillaume para Simone
                    kfl_inlev_lev = &
                         -1000-getint('VALUE',1_ip,'#Initial condition is from value function')

                 end if
              else if( words(1) == 'VARIA' ) then
                 !
                 ! ADOC[2]> VARIATION:         CONSTANT | NON_CONSTANT                                     $ Boundary conditions constant/variable in time
                 ! ADOC[d]> VARIATION:
                 ! ADOC[d]> If a time function is used, NON_CONSTANT option should be selected.
                 ! ADOC[d]> Alya does not identify non-constant boundary conditions automatically.
                 !
                 ! if( words(2) == 'NONCO' .and. kfl_conbc_nsi /= 0 ) then
                 if( words(2) == 'NONCO' ) then
                    kfl_conbc_lev = 0
                 else
                    kfl_conbc_lev = 1
                 end if

              end if
              call ecoute('lev_reabcs')


!!$              !Start SM
!!$              if(words(1)=='VALUE') then
!!$                 kdime = 0
!!$                 if( exists('FUNCT') ) then
!!$                    nvcod = getint('FUNCT',1_ip,'*FUNCTION NUMBER')
!!$                    print*,' NVCOD lev: ', nvcod
!!$                    if( nvcod > mvcod .or. nvcod <= 0 ) &
!!$                         call runend('REABCS: WRONG DEFAULT VALUE FUNCTION')
!!$                 end if
!!$                 if( exists('DIMEN') ) then
!!$                    kdime = getint('DIMEN',1_ip,'*DIMENSION')
!!$                 end if
!!$                 icond = mcodb+1
!!$                 if( exists('CODE ') ) then
!!$                    icond = getint('CODE ',0_ip,'*DEFAULT CODE ON NODES')
!!$                 end if
!!$                 if( nvcod < 1 ) call runend('REABCS: WRONG CODE NUMBER')
!!$                 call ecoute('lev_reabcs')
!!$                 if( kdime == 0 ) kdime = nnpar - 1
!!$                 ivcod = ivcod + 1
!!$                 call lev_memall(9_ip)                
!!$                 if( icond /= mcodb + 1 ) then
!!$                    do while( words(1) /= 'ENDVA' )
!!$                       ipoin = int(param(1))
!!$                       !kfl_codno(1,ipoin) = icond
!!$                       !kfl_codno(2:mcono,ipoin) = mcodb+1
!!$                       do idime = 1,kdime
!!$                          bvcod_lev(ipoin) = param(idime+1)
!!$                          !bvcod(nvcod)%a(idime,ipoin) = param(idime+1)
!!$                       end do
!!$                       call ecoute('lev_reabcs')
!!$                    end do
!!$                 else
!!$                    do while( words(1) /= 'ENDVA' )
!!$                       ipoin = int(param(1))
!!$                       do idime = 1,kdime
!!$                          bvcod_lev(ipoin) = param(idime+1)
!!$                          !print*, ipoin, bvcod_lev(ipoin)
!!$                          !bvcod(nvcod)%a(idime,ipoin) = param(idime+1)
!!$                       end do
!!$                       call ecoute('lev_reabcs')
!!$
!!$                    end do
!!$                 end if                 
!!$              end if
!!$              !tncod => tncod_lev
!!$              !call reacod(10_ip)
!!$              ! !!
!!$              !END SM
!!$              ! !!
           end do

        else if( words(1) == 'CODES' .and. exists('NODES') ) then 

           !-------------------------------------------------------------
           !
           ! User-defined codes on nodes
           !           
           !-------------------------------------------------------------

           if(exists('GEOME')) then
              !
              ! Geometrical b.c.
              !              
              tgcod => tgcod_lev
              call reacod(4_ip) 

           else
              !
              ! Node boundary condition
              !
              tncod => tncod_lev
              call reacod(1_ip)
           end if

        else if(words(1)=='ONNOD') then

           call runend('LEV_REABCS: THIS OPTION DOES NOT EXIST ANY MORE')

        end if

     end do

  end if

end subroutine lev_reabcs
