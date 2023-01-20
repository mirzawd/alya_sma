!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup DomainInput
!> @{
!> @file    readim.f90
!> @author  Guillaume Houzeaux
!> @date    18/09/2012
!> @brief   Read mesh dimensions
!> @details Read mesh dimensions:
!>          \verbatim
!>          - NPOIN ...................... Number of nodes
!>          - NELEM ...................... Number of elements
!>          - NDIME ...................... Number of space dimensions (1,2,3)
!>          - NBOUN ...................... Number of boundary elements
!>          - NHANG....... ............... Number of hanging nodes
!>          - LEXIS(NELTY) ............... List of existing element types
!>          - NFIEL ...................... Number of fields
!>          - KFL_FIELD(5,NFIEL) ......... Dimension/type of fields,steps
!>          - KFL_FIELD(1,II) ... Dimensions (1, NDIMNE, NTENS, etc,)
!>          - KFL_FIELD(2,II) ... Type (on elements, boundaries, nodes)
!>          - KFL_FIELD(3,II) ... Negative values excluded (=1)
!>          - KFL_FIELD(4,II) ... Steps to be saved
!>          - KFL_FIELD(5,II) ... Size (=nelem,nboun,npoin)
!>          - KFL_FIELD(6,II) ... 0: read at beginning, 1=read when needed
!>          \endverbatim
!> @}
!-----------------------------------------------------------------------

subroutine readim()

  use def_parame
  use def_master
  use def_domain
  use def_inpout
  use def_elmtyp
  use mod_memchk
  use mod_ecoute,        only : ecoute
  use mod_reafie,        only : reafie_read_header
  use mod_elmgeo,        only : elmgeo_element_name_to_type
  use mod_domain,        only : domain_memory_allocate
  use mod_ecoute,        only : ecoute_set_read_unit
  use mod_ecoute,        only : ecoute_set_write_unit

  implicit none

  integer(ip) :: ielty,ipara,inode,ifiel
#ifdef NDIMEPAR
  integer(ip) :: ndime_aux
#endif

  if( INOTSLAVE ) then

     call ecoute_set_read_unit (lun_pdata_dom) ! Reading file
     call ecoute_set_write_unit(lun_outpu_dom) ! Writing file
     !
     ! Initializations
     !
#ifdef NDIMEPAR
     continue                   ! NDIME is a parameter
#else
     ndime     = -1             ! Obligatory
#endif
     kfl_autbo =  0             ! Automatic boundaries off
     npoin     = -1             ! Obligatory
     nelem     = -1             ! Obligatory
     nboun     = -1             ! Optional
     nboib     =  0             ! Optional
     npoib     =  0             ! Optional
     nhang     =  0             ! Optional
     nfiel     =  0             ! Number of fields
     kfl_field =  0             ! Fields dimensions
     mcodb     =  99            ! Optional, max boundary codes
#ifndef PNODE_VALUE
     mnode     = -1             ! Optional (useful for virtual elements)
#endif
     !
     ! Reach the section
     !
     call ecoute('readim')
     do while( words(1) /= 'DIMEN' )
        call ecoute('readim')
     end do

     call ecoute('readim')
     !
     !
     !.md<module>kernel
     !.md<input>case.dom.dat
     !.md<pos>0
     !.md<sec>
     !.md<0># Dimensions of the problem
     !.md<>
     !.md<code>
     !.md<0><b>DIMENSIONS</b>
     do while( words(1) /= 'ENDDI')

        if( words(1) == 'NODAL' ) then
           !
           !.md<1><b>NODAL_POINTS</b>   = int                                                          $ Number of nodes
           !.md<field>NODAL_POINTS
           !.md<com>Number of nodes
           !
           npoin = getint('NODAL',0_ip,'*NUMBER OF NODAL POINTS')

        else if( words(1) == 'ELEME' ) then
           !
           !.md<1><b>ELEMENTS</b>       = int                                                          $ Number of elements
           !.md<field>ELEMENTS
           !.md<com>Number of elements
           !
           nelem = getint('ELEME',0_ip,'*NUMBER OF ELEMENTS')
           nelwh = nelem

        else if( words(1) == 'MAXIM' ) then
           !
           !.md<1>MAXIMUM        = int                                                                 $ Max number of nodes per element
           !.md<field>MAXIMUM
           !.md<com>Maximum number of nodes per element
           !
#ifndef PNODE_VALUE
           mnode = getint('MAXIM',0_ip,'*MAXIMUM NUMBER OF NODES PER ELEMENT')
#endif
        else if( words(1) == 'SPACE' ) then
           !
           !.md<1><b>SPACE</b>          = int                                                          $ Space dimension (1,2,3)
           !.md<field>SPACE
           !.md<com>Space dimension
           !
#ifdef NDIMEPAR
           ndime_aux = getint('SPACE',2_ip,'*SPACE DIMENSION')
           if (ndime_aux/=ndime) call runend('READIM: running with ifdef NDIMEPAR and ndime_aux/=ndime')
#else
           ndime = getint('SPACE',2_ip,'*SPACE DIMENSION')
#endif
        else if( words(1) == 'FIELD' ) then
           !
           !.md<1>FIELDS         = int, PERIODIC                                                $ Number of fields (0 by default)
           !.md<field>FIELDS
           !.md<com>Number of fields. `PERIODIC`: optional setting. If specified in `FIELDS`, all the fields with timesteps will be treated as periodic.
           !.md<com>If specified in `FIELD`, only that field will be treated as periodic.
           !
           ! PERIODIC - optional setting. If specified in FIELDS, all the fields with timesteps will be treated as periodic. If specified in FIELD, only that field will be treated as periodic.
           ! As a consequence, it is expected that the last timestep is equal to the first timestep.
           ! If the field is not periodic and the simulation tries to continue beyond the last proided time step - Alya will crash
           ! For example, a field with period 0.6s:
           ! DIMENSIONS
           !   FIELDS = 1, PERIODIC
           !       FIELD: 1, DIMENSIONS=3, NODES, STEPS=11, PERIODIC, ONDEMAND
           !          TIMES
           !          0.0
           !          0.2
           !          0.4
           !          0.6
           !          END_TIMES
           !       END_FIELD
           !   END_FIELDS
           ! END_DIMENSIONS
           ! FIELDS
           !    FIELD: 1
           !       STEP 1
           !          INCLUDE veloc.0.in
           !       END_STEP
           !       STEP 2
           !          INCLUDE veloc.1.in
           !       END_STEP
           !       STEP 3
           !          INCLUDE veloc.2.in
           !       END_STEP
           !       STEP 4
           !          INCLUDE veloc.0.in
           !       END_STEP
           !    END_FIELD
           ! END_FIELDS


           if( exists('NUMBE' ) ) then
              nfiel = getint('NUMBE',1_ip,'#NUMBER OF FIELDS')
           else
              nfiel = getint('FIELD',1_ip,'#NUMBER OF FIELDS')
           end if
           !.md<2>FIELD:int, DIMENSIONS=int, NODES | ELEMENTS | BOUNDARIES, <i>EXCLUDED</i>, <i>STEPS=int</i>
           !.md<field>FIELD
           !.md<com>Field metadata
           !.md<3><i>TIME
           !.md<4>real
           !.md<4>real
           !.md<4>...
           !.md<3>END_TIMES</i>
           !.md<2>END_FIELD
           !.md<1>END_FIELDS

           if( exists('PERIO') ) then
              !print *,'All fields will be treated as periodic'
              kfl_field(7,:) = 1_ip
           end if

           call ecoute('readim')
           do while( words(1) /= 'ENDFI')
              if( words(1) == 'FIELD') then
                 call reafie_read_header(ifiel)
                 call domain_memory_allocate('TIME_FIELD % A',NUMBER1=ifiel)
                 if( exists('STEPS') ) then
                    ipara = getint('STEPS',1_ip,'#NUMBER OF STEPS')
                    if( ipara > 1 ) then
                       call ecoute('readim')
                       do while( words(1) /= 'ENDFI' )
                          if( words(1) == 'TIMES' ) then
                             ipara = 0
                             call ecoute('readim')
                             do while( words(1) /= 'ENDTI' )
                                ipara = ipara + 1
                                if( ipara > kfl_field(4,ifiel) ) call runend('READIM: WRONG NUMBER OF STEPS WHILE READING FIELD')
                                time_field(ifiel) % a(ipara) = param(1)
                                call ecoute('readim')
                             end do
                          end if
                          call ecoute('readim')
                       end do
                    end if
                 end if
              end if
              call ecoute('readim')
           end do

        else if( words(1) == 'BOUND' ) then
           !
           !.md<1>BOUNDARIES     = AUTOMATIC/int                                                $ Number of boundary elements
           !.md<field>BOUNDARIES
           !.md<com>This field is not compulsory. Put `AUTOMATIC` to generate the boundary automatically.
           !.md<com>In this case no reference to the boundary is possible and it is only useful for postprocess purpose.
           !
           if(words(2) == 'AUTOM' ) then
              kfl_autbo = 1
           else
              nboun = getint('BOUND',0_ip,'*NUMBER OF BOUNDARY ELEMENTS')
           end if

        else if( words(1) == 'MCODB' ) then
           !
           !.md<1>MCODB          = int                                                          $ Maximum number of code
           !.md<field>MCODB
           !.md<com>Maximum number of codes
           !
           mcodb = getint('MCODB',0_ip,'*MCODB')

        else if( words(1) == 'TYPES' ) then
           !
           !.md<1>TYPES          = char, char...                                                $ List of volume elements in the mesh
           !.md<field>TYPES
           !.md<com>Elements available:
           !.md<com>    -  1D: BAR02,BAR03,BAR04
           !.md<com>    -  2D: TRI03,TRI06,QUA04,QUA08,QUA09,QUA16
           !.md<com>    -  3D: TET04,TET10,PYR05,PYR14,PEN06,PEN15,PEN18,HEX08,HEX20,HEX27,HEX64
           !
           if( nnpar == 0 ) then
              ipara = 2
              do while( trim(words(ipara)) /= '' )
                 call elmgeo_element_name_to_type(words(ipara),ielty)
                 if( ielty >= 2 .and. ielty <= nelty ) lexis(ielty) = 1
                 ipara = ipara + 1
              end do
           else
              do ipara = 1,nnpar
                 ielty = int(param(ipara))
                 if( ielty >= 2 .and. ielty <= nelty ) lexis(ielty) = 1
              end do
           end if

        else if( words(1) == 'NODES' ) then
           do ipara = 1,nnpar
              inode = int(param(ipara))
              if( ndime ==0 ) &
                   call runend('READIM: SPACE DIMENSION SHOULD BE DEFINED FIRST')
              if( inode >= 1 ) then
                 call fintyp(ndime,inode,ielty)
                 lexis(ielty) = 1
              end if
           end do

        else if( words(1) == 'HANGI' ) then
           !
           !.md<1>HANGING_NODES  = int                                                          $ Number of hanging nodes
           !.md<field>HANGING_NODES
           !.md<com>Number of hanging nodes
           !
           nhang = getint('HANGI',0_ip,'*NUMBER OF HANGING NODES')

        end if
        call ecoute('readim')
     end do
     !
     !.md<0><b>END_DIMENSIONS</b>
     !.md<0><b>FIELDS</b>
     !.md<1>FIELD: int              $ Place this group at the end of the file
     !.md<2>STEP 1                  $ STEP/END_STEP only if TIME vector is present in FIELD definition
     !.md<3>1 float float ...
     !.md<3>2 float float ...
     !.md<3>...
     !.md<2>END_STEP
     !.md<2>STEP 2
     !.md<3>1 float float ...
     !.md<3>2 float float ...
     !.md<3>...
     !.md<2>END_STEP
     !.md<2>...
     !.md<1>END_FIELD
     !.md<0><b>END_FIELDS</b>     
     !.md</code>
     !.md<>
     !
     ! Check errors
     !
     if( npoin == -1 ) call runend('READIM: WRONG NUMBER OF NODES')
     if( nelem == -1 ) call runend('READIM: WRONG NUMBER OF ELEMENTS')
     if( ndime == -1 ) call runend('READIM: WRONG NUMBER OF SPACE DIMENSIONS')

     nboun   = max( nboun , 0_ip )
     npoin_2 = npoin
     nelem_2 = nelem
     nboun_2 = nboun
     !
     ! Manage unique type of elements
     !
     nexis = 0
     nexis = sum(lexis)
     utype = -1
     if( nexis == 1 ) then
        do ielty = 1,nelty
            if( lexis(ielty) == 1 ) then
                utype = ielty
            end if
        end do
     end if

  end if

end subroutine readim
