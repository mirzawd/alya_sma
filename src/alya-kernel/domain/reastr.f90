!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup DomainInput
!> @{
!> @file    reastr.f90
!> @author  Guillaume Houzeaux
!> @date    18/09/2012
!> @brief   Read mesh strategy
!> @details Read mesh strategy:       
!>          \verbatim
!>          - LQUAD(NELTY) ... Type of integration rule (1 = closed, 0 = open)
!>          - NGAUS(NELTY) ... Number of domain integration points
!>          - NGROU_DOM ...... Number of groups for deflation based solvers
!>          - SCALE(3) ....... Geometrical scale factor
!>          - TRANS(3) ....... Geometrical translation vector
!>          - KFL_GEOME ...... If geometrical conditions exist
!>          - KFL_CONVX ...... What to do with convex node (if geometrical conditions are applied)
!>          - KFL_ICODN ...... If there are codes on nodes
!>          - KFL_ICODB ...... If the are codes on boundaries
!>          - MCONO .......... Max # codes per nodes
!>          - KFL_EXTRA ...... If boundary codes are extrapolated to nodes
!>          - NPBCS .......... Max # of possibe geometrical conditions
!>          - GEOAN .......... Geometrical angle to decide if we have corners (if geometrical conditions are applied)
!>          - LSBCS(8,100) ... List of boundaries for a given geometrical code (first argument is condition)
!>          \endverbatim
!> @} 
!-----------------------------------------------------------------------
subroutine reastr()
  use def_parame
  use def_master 
  use def_domain
  use def_inpout
  use def_elmtyp
  use mod_memchk
  use mod_iofile
  use mod_ecoute,        only : ecoute
  use mod_reabcs,        only : reabcs_geometrical
  use mod_elmgeo,        only : elmgeo_element_name_to_type
  use mod_memory,        only : memory_alloca, memory_size
  use mod_messages,      only : messages_live
  use def_isoparametric, only : LAGRANGE_INTERPOLATION 
  use def_isoparametric, only : CHEBYSHEV_INTERPOLATION
  use def_quadrature,    only : GAUSS_LEGENDRE_RULE   
  use def_quadrature,    only : TRAPEZOIDAL_RULE      
  use def_quadrature,    only : CLOSED_RULE           
  use def_quadrature,    only : CHEBYSHEV_RULE           

  implicit none
  integer(ip) :: ielty,jelty,imate,jmate,icode,ifiel
  integer(ip) :: ncodes
  !
  ! Allocate memory
  ! 
  if( INOTSLAVE ) then

     lquad              =  0          ! Open integration rule
     ngaus              =  0          ! No Gauss points
     linte              =  0          ! Interpolation function (Lagrange by default)

     kfl_ngrou          =  0          ! Strategy to construct groups
     ngrou_dom          =  0          ! Number of groups (for deflated): -2 (slaves), -1 (automatic, master), >0 (master)
     ngrou_dom_target   =  0          ! Number of groups (for deflated): -2 (slaves), -1 (automatic, master), >0 (master)
     ngrou_boxes_coarse =  1          ! Number of coarse boxes for SFC
     ngrou_boxes_fine   =  128        ! Number of fines boxes for SFC

     xscal(1)           =  1.0_rp     ! X scale factor
     xscal(2)           =  1.0_rp     ! Y scale factor
     xscal(3)           =  1.0_rp     ! Z scale factor

     trans(1)           =  0.0_rp     ! X translation
     trans(2)           =  0.0_rp     ! Y translation 
     trans(3)           =  0.0_rp     ! Z translation

     kfl_extra          =  0          ! Extrapolate from boundary to node

     kfl_geome          =  0          ! Do not compute geometrical normals
     kfl_convx          =  1          ! Use exnor for convex angle nodes
     kfl_frees          =  0          ! Freestream criteria. 0: use wind angle to decide inflow/outflow
     npbcs              =  0          ! Max # of possibe geometrical conditions
     lsbcs              =  0          ! List of boundaries for a given geometrical code
     awind              =  0.0_rp     ! Wind angle for freestream
     tolan              =  5.0_rp*pi/180.0_rp ! Tolerance used to define inflow from freestream
     geoan              = 45.0_rp     ! Angle for geometrical normals

     kfl_chege          =  0          ! Don't check geometry
     kfl_naxis          =  0          ! Cartesian coordinate system
     kfl_spher          =  0          ! Cartesian coordinate system
     kfl_bouel          =  1          ! Element # connected to boundary is known
     kfl_divid          =  0          ! Divide elements into tetra
     curvatureDataField =  0          ! curvature data field id
     curvatureField     =  0          ! curvature data field id

     ncodes             = -1_ip
     !
     ! Reach the section 
     !
     imate = 1
     jmate = 1
     call ecoute('REASTR')
     do while( words(1) /= 'STRAT' )
        call ecoute('REASTR')
     end do
     !
     !
     !.md<module>kernel
     !.md<input>case.dom.dat
     !.md<pos>1
     !.md<sec>
     !.md<0># Strategy for the mesh related operations
     !.md<>
     !.md<code>
     !.md<0><b>STRATEGY</b>
     do while(words(1) /= 'ENDST')
        call ecoute('REASTR')

        if( words(1) == 'INTEG' ) then
           !
           !.md<1>INTEGRATION_RULE = Open/Close, Open/Close...              $ Integration rule (Open by default)
           !.md<field>INTEGRATION_RULE
           !.md<com>List of intgeration rules for each of the element type declared in DIMENSION field.
           !.md<com>The list should be written in the following ordering:
           !.md<com>    - 1D: BAR02,BAR03,BAR04
           !.md<com>    - 2D: TRI03,TRI06,QUA04,QUA08,QUA09,QUA16
           !.md<com>    - 3D: TET04,TET10,PYR05,PYR14,PEN06,PEN15,PEN18,HEX08,HEX20,HEX27,HEX64
           !.md<com>
           !.md<com>If only one option is given, it is applied to all element types.
           !
           if( words(3) == '' ) then
              select case ( words(2) )
              case ( 'OPEN ' , 'GAUSS' ) ; lquad = GAUSS_LEGENDRE_RULE
              case ( 'CLOSE'           ) ; lquad = CLOSED_RULE
              case ( 'TRAPE'           ) ; lquad = TRAPEZOIDAL_RULE
              case ( 'CHEBY'           ) ; lquad = CHEBYSHEV_RULE
              case default               ; lquad = GAUSS_LEGENDRE_RULE
              end select
           else              
              do ielty = 1,min(nelty,size(words,KIND=ip)-1)
                 select case ( words(ielty+1) )
                 case ( 'OPEN ' , 'GAUSS' ) ; lquad(ielty) = GAUSS_LEGENDRE_RULE
                 case ( 'CLOSE'           ) ; lquad(ielty) = CLOSED_RULE
                 case ( 'TRAPE'           ) ; lquad(ielty) = TRAPEZOIDAL_RULE
                 case ( 'CHEBY'           ) ; lquad(ielty) = CHEBYSHEV_RULE
                 case default               ; lquad(ielty) = GAUSS_LEGENDRE_RULE
                 end select
              end do
           end if

        else if( words(1) == 'INTER' ) then
           !
           !.md<1>INTERPOLATION_FUNCTION = Lagrange/Chebyshev...              $ Interpolation function
           !.md<field>INTEGRATION_FUNCTION
           !.md<com>Interpolaiton method:
           !.md<com>    - Lagrange for Lagrange polynomial interpolation functions
           !.md<com>    - Chebyshev for Chebyshev polynomial interpolation functions
           !
           if( words(3) == '' ) then
              select case ( words(2) )
              case ( 'LAGRA' ) ; linte = LAGRANGE_INTERPOLATION
              case ( 'CHEBY' ) ; linte = CHEBYSHEV_INTERPOLATION
              case default     ; call runend('REASTR: UNKNOWN INTERPOLATION FUNCTION')
              end select
           else              
              do ielty = 1,min(nelty,size(words,KIND=ip)-1)
                 select case ( words(ielty+1) )
                 case ( 'LAGRA' ) ; linte(ielty) = LAGRANGE_INTERPOLATION
                 case ( 'CHEBY' ) ; linte(ielty) = CHEBYSHEV_INTERPOLATION
                 case default     ; call runend('REASTR: UNKNOWN INTERPOLATION FUNCTION')
                 end select
              end do
           end if

        else if( words(1) == 'DOMAI' ) then 
           !
           !.md<1>DOMAIN_INTEGRATION_POINTS = int, int...                   $ Number of Integration points
           !.md<field>DOMAIN_INTEGRATION_POINTS
           !.md<com>Number of integration points (0 for automatic) for each of the element type declared in DIMENSION field.
           !.md<com>When automatic, the number of integration points COINCIDES with the number of element nodes.
           !.md<com>To explicitly give a value, the list should be written using the following element ordering:
           !.md<com>    - 1D: BAR02,BAR03,BAR04
           !.md<com>    - 2D: TRI03,TRI06,QUA04,QUA08,QUA09,QUA16
           !.md<com>    - 3D: TET04,TET10,PYR05,PYR14,PEN06,PEN15,PEN18,HEX08,HEX20,HEX27,HEX64
           !
           jelty = 0
           do ielty = 1,nelty
              if( lexis(ielty) > 0 ) then
                 jelty = jelty + 1
                 ngaus(ielty) = int(param(jelty),ip)
                 if( ngaus(ielty) == 1 ) lquad(ielty) = 0      ! One gauss point forced to open rule 
              end if
           end do

        else if( words(1) == 'GAUSS' ) then 
           !
           ! Gauss points
           !
           call ecoute('REAST')
           do while( words(1) /= 'ENDGA' )
              call elmgeo_element_name_to_type(words(1),ielty)
              if( ielty < 1 .or. ielty > nelty ) then
                 call runend('REASTR: WRONG ELEMENT TYPE')
              else
                 ngaus(ielty) = int(param(1),ip)
              end if
              call ecoute('REAST')
           end do

        else if( words(1) == 'GROUP' ) then
           !
           !.md<1>GROUPS = int, SEQUENTIAL_FRONTAL | PARALLEL_FRONTAL | PARTITION | SFC | FLOOR | FIELD $ Number of groups (used for deflation)
           !.md<field>GROUPS
           !.md<com>This option is the number of groups to be constructed for deflation based solvers.
           !.md<com>    - GROUPS = -1 ... The number of groups is computed in an automatic way.
           !.md<com>    - GROUPS = -2 ... In parallel, one group per subdomain.
           !.md<com>    - GROUPS >  0 ... Number of groups specified by the user (highly recommanded).
           !.md<com>    - Options:
           !.md<com>        - `SEQUENTIAL_FRONTAL` (default)
           !.md<com>        - `PARALLEL_FRONTAL` (whole parallel execution)
           !.md<com>        - `PARTITION`
           !.md<com>        - `FLOOR`
           !.md<com>        - `SFC, COARSE=int | FINE=int`
           !.md<com>        - `FIELD`

           ngrou_dom = getint('GROUP',0_ip,'#NUMBER OF GROUPS')
           if(      exists('SEQUE') .or. exists('FRONT' ) ) then
              kfl_ngrou = -1
           else if( exists('PARAL') ) then
              kfl_ngrou = -2
           else if( exists('PARTI') ) then
              kfl_ngrou = -3
           else if( exists('SFC  ') ) then
              kfl_ngrou = -4
              if( exists('COARS') ) ngrou_boxes_coarse = getint('COARS',1_ip  ,'#Number of coarse bins')
              if( exists('FINE ') ) ngrou_boxes_fine   = getint('FINE ',128_ip,'#Number of fine bins')
           else if( exists('FIELD') ) then
              kfl_ngrou = getint('FIELD',1_ip,'#groups from field')
           else if( exists('FLOOR') ) then
              kfl_ngrou = -5
           else
              kfl_ngrou = -1
           end if
           ngrou_dom_target = ngrou_dom
           
        else if(words(1)=='SCALE') then 
           !
           !.md<1>SCALE : XSCALE = real, YSCALE = real, ZSCALE = real       $ Scaling factors for the geometry
           !.md<field>SCALE
           !.md<com>Multiply each coordinates by its corresponding scaling factor.
           !.md<com>In postprocess, the geometry is written using this scaling.
           !
           xscal(1) = getrea('XSCAL',1.0_rp,'#x-factor')
           xscal(2) = getrea('YSCAL',1.0_rp,'#y-factor')
           xscal(3) = getrea('ZSCAL',1.0_rp,'#z-factor')

        else if( words(1) == 'TRANS' ) then 
           !
           !.md<1>TRANSLATION = XTRANS = real, YTRANS = real, ZTRANS = real $ Translation of the geometry
           !.md<field>TRANSLATION
           !.md<com>Translate the geometry.
           !
           trans(1) = getrea('XTRAN',0.0_rp,'#x-translation')
           trans(2) = getrea('YTRAN',0.0_rp,'#y-translation')
           trans(3) = getrea('ZTRAN',0.0_rp,'#z-translation')

        else if( words(1) == 'GEOME' ) then
           !
           ! Read geometrical boundary conditions
           !
           call reabcs_geometrical()

        else if( words(1) == 'BOUND' ) then
           !
           ! Compute LELBO because it is unknown
           !
           if( words(2) == 'KNOWN' ) then
              kfl_bouel = 1
           else if( words(2) == 'UNKNO' ) then
              kfl_bouel = 0
           else if( option('BOUND') ) then
              kfl_bouel = 1
           else
              kfl_bouel = 0
           end if

        else if( words(1) == 'CHECK' ) then
           !
           ! Check geometry
           !
           if( option('CHECK') ) kfl_chege = 1

        else if( words(1) == 'CRVDA' ) then
           !
           ! Reading curvatureDataField
           !
           curvatureDataField = getint('CRVDA',1_ip,'#Curve dat field id')

        else if( words(1) == 'CRVGE' ) then
           !
           ! Reading curvatureField
           !
           curvatureField = getint('CRVGE',1_ip,'#Curved geometry field id')

        else if( words(1) == 'AXISY' ) then
           !
           ! Cylindrical coordinates
           !
           if( option('AXISY') ) kfl_naxis = 1

        else if( words(1) == 'SPHER' ) then
           !
           ! Spherical coordinates
           !
           if( option('SPHER') ) kfl_spher = 1

        else if( words(1) == 'DIVID' ) then
           !
           ! Divide mesh into tetras
           !
           if( words(2) == 'ON   ' .or. words(2) == 'YES  ' ) kfl_divid = 1

        else if( exists('EXTRA') ) then
           !
           ! Extrapolate boundary conditions
           !
           if( option('EXTRA') ) kfl_extra = 1

        else if( exists('FIELD') ) then
           !
           ! Extrapolate boundary conditions
           !
           call memory_alloca(memor_dom,'XSCAL_FIELDS','reastr',xscal_fields,nfiel,INIT_VALUE=1.0_rp)
           call ecoute('reastr')
           kfl_xscal_fields = 1
           do while( words(1) /= 'ENDFI' )
              ifiel               = int(param(1),ip)
              xscal_fields(ifiel) = param(2)              
              call ecoute('reastr')
           end do
           
        else if( words(1) == 'MATER' ) then

           !-------------------------------------------------------------
           !
           ! Materials depending on boundary codes
           !
           !-------------------------------------------------------------
           !
           !.md<1>MATERIALS, BOUNDARIES, NCODES=int, DEFAULT= int
           !.md<2>CODE =     int,int...                                                        $ From which codes material should be generated
           !.md<2>MATERIAL = int                                                               $ Material number
           !.md<2>LAYERS =   int                                                               $ Number of element layers
           !.md<2>CODE =     int,int...                                                        $ From which codes material should be generated
           !.md<2>MATERIAL = int                                                               $ Material number
           !.md<2>LAYERS =   int                                                               $ Number of element layers
           !.md<2>...
           !.md<1>END_MATERIALS
           !.md<field>MATERIALS, BOUNDARIES
           !.md<com>This option enables to generate automatically a certain number of element layers (LAYERS= int) of one material (MATERIAL= int)
           !.md<com>starting from boundaries with a given codes (CODE= int, ...). NCODES - is the number of codes to expect in total for all the materials. 
           !.md<com>The option DEFAULT= int assign automatically this material to all the element before applying the layers. 
           !.md<com>It can be useful for example in CFD to define a different material at outflows on which the viscosity will be higher.
           !
           !TODO: in the future it would be better to throw an error instead of defaulting to 5 which makes no sense 
           !if ( exists('NCODE') ) then 
               ncodes = getint('NCODE',5_ip,'#Number of boundary codes for the layers') !for backwards compatibilty default is 5
               
               if (ncodes < 1_ip) then
                  call runend("NCODES<1 in MATERIAL, BOUNDARY in STRATEGY group od dom.dat")
               end if

               call memory_alloca( memor_dom, 'materials_nlaye', 'reastr', materials_nlaye, ncodes)
               call memory_alloca( memor_dom, 'materials_icode', 'reastr', materials_icode, ncodes)
               call memory_alloca( memor_dom, 'materials_imate', 'reastr', materials_imate, ncodes)
               materials_nlaye    = 0_ip           ! Automatic generation of materials
               materials_icode    = 0_ip           ! Automatic generation of materials
               materials_imate    = 1_ip           ! Automatic generation of materials
           !else
           !    call runend("Missing NCODES in MATERIAL, BOUNDARY in STRATEGY group od dom.dat")
           !end if

           if( words(2) == 'BOUND' ) then

              icode = 0_ip
              jmate = 0_ip
              call ecoute('reastr')
              do while(words(1)/='ENDMA')
                 if( words(1)== 'CODE ' ) then
                    imate = jmate + 1_ip
                    jmate = imate + nnpar - 1_ip
                    icode = nnpar
                    if( jmate > memory_size(materials_icode) ) call runend('REASTR: MATERIAL, BOUNDARY(.dom.dat) HAS MORE CODES ('//trim(intost(nnpar))//') THAN ALLOCATED (NCODES='//trim(intost(ncodes))//')')                    
                    materials_icode(imate:jmate) = int(param(1:nnpar),ip)
                 else
                    if( icode == 0 ) then
                       call runend('REASTR: IN MATERIAL, BOUNDARY(.dom.dat), CODE FIELD SHOULD BE DEFINED FIRST')
                    else if(  words(1)== 'MATER' ) then
                       materials_imate(imate:jmate) = int(param(1),ip)
                    else if( words(1) == 'LAYER' ) then
                       materials_nlaye(imate:jmate) = int(param(1),ip)
                    end if
                 end if
                 call ecoute('reastr')
               end do

           else
              call messages_live("MATERIAL in STRATEGY group of dom.dat is without BOUNDARY layers will not be generatet","WARNING")
           end if

           !---------------------------------------------------
           !
           ! Save important info to log
           !
           if (associated(materials_icode)) then
               write( momod(ID_KERMOD)%lun_outpu, * ) "AUTOMATICALLY GENERATED MATERIAL LAYERS"
               write( momod(ID_KERMOD)%lun_outpu, * ) "    CODE MATERIAL NLAYERS"
               do icode = 1, ncodes
                  write( momod(ID_KERMOD)%lun_outpu, '(2x, I7, 2x, I7, 1x, I7)' )  materials_icode(icode), materials_imate(icode), materials_nlaye(icode)
               end do
               flush( momod(ID_KERMOD)%lun_outpu )
           end if
           !
           ! Save important info to log
           !
           !---------------------------------------------------


        end if

     end do
     !.md<0><b>END_STRATEGY</b>
     !.md</code>
     !.md<>

  end if

end subroutine reastr
