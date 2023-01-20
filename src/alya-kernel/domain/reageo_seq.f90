!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup DomainInput
!> @{
!> @file    reageo.f90
!> @author  Guillaume Houzeaux
!> @date    18/09/2012
!> @brief   Read mesh arrays
!> @details Read mesh arrays:
!>          \verbatim
!>          - Element arrays:
!>            - LEXIS(NELTY) ............... Wether an element type exists or not
!>            - LNODS(MNODE,NELEM) ......... Element connectivity
!>            - LTYPE(NELEM) ............... Types of element
!>            - LESUB(NELEM) ............... Subdomains
!>            - LMATE(NELEM) ............... Materials
!>            - LELCH(NELEM) ............... Element characteristic: 0,1,2 (finite element, extension, hole)
!>         - Nodal arrays:
!>            - COORD(NDIME,NPOIN) ......... Node coordinates
!>            - LNOCH(NPOIN) ............... Node characteristic: 0,2 (finite element, hole)
!>            - LMAST(NPOIN) ............... List of master nodes for periodicity
!>          - Boundary arrays:
!>            - LNODB(MNODB,NBOUN) ......... Boundary connectivity
!>            - LELBO(NBOUN) ............... List of elements connected to boundaries
!>            - LTYPB(NBOUN) ............... Types of boundary
!>            - LBOCH(NBOUN) ............... Boundary characteristic: 0,1 (finite boundary, extension)
!>          \endverbatim
!> @}
!-----------------------------------------------------------------------
subroutine reageo_seq()
  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use def_inpout
  use mod_iofile
  use mod_memory
  use def_elmtyp
  use mod_mpio_seq_log
  use mod_reafie
  use mod_elmgeo,             only : elmgeo_element_type
  use mod_ecoute,             only : ecoute
  use mod_messages,           only : messages_live
  use mod_read_domain_arrays, only : read_domain_arrays_types
  use mod_domain,             only : domain_memory_allocate
  implicit none
  integer(ip)           :: ielem,ipoin,jpoin,inode,idime
  integer(ip)           :: iboun,inodb,ielty,ktype,dummi
  integer(ip)           :: imate,ktypb,klelb,iblty,knode,ii
  integer(ip)           :: kfl_icgns,ipoin_master
  integer(ip)           :: isubd,ipoin_slave

  if( ISEQUEN .or. ( IMASTER .and. kfl_ptask /= 2 ) ) then
     !
     ! Local variables
     !
     kfl_icgns = 0                   ! New Alya format for types of elements
     ktype     = 0                   ! Element # of nodes is not calculated
     ktypb     = 0                   ! Boundary types have been read
     klelb     = 0                   ! If boundary elements have been read
     !
     ! Memory allocation
     !
     call domain_memory_allocate('GEOMETRY')
     !
     ! Unique type of element
     !
     if( utype > 0 ) then
        utype = 0
        !ltype = utype
        !ktype = nelem
     end if
     !
     ! Read options and arrays
     !
     call ecoute('REAGEO')
     do while( words(1) /= 'GEOME' )
        call ecoute('REAGEO')
     end do
     !.md<module>kernel
     !.md<input>case.dom.dat
     !.md<pos>2
     !.md<sec>
     !.md<0># Geometry
     !.md<>
     !.md<code>
     !.md<0><b>GEOMETRY</b>
     do while(words(1)/='ENDGE')
        call ecoute('reageo')

        if( words(1) == 'NODES' ) then
           !
           ! LTYPE: Element types
           !
           call start_timer('LTYPE', 'READ ', 'ASCII', 1_ip, 1_ip)
           call messages_live('READ ELEMENT  TYPE')
           if( utype > 0 ) then
              do while( words(1) /= 'ENDNO' )
                 call ecoute('reageo')
              end do
           else
              ktype = nelem
              do ielem = 1,nelem
                 read(nunit,*,err=1) dummi,knode
                 ielty        = elmgeo_element_type(ndime,knode)
                 lexis(ielty) = 1
                 ltype(ielem) = ielty
              end do
              call ecoute('reageo')
              if( words(1) /= 'ENDNO' ) &
                   call runend('REAGEO: ERROR READING NODES PER ELEMENT, POSSIBLY MISMATCH IN NUMBER OF NODES IN DOM.DAT')
           end if
           call end_timer()

        else if( words(1) == 'TYPES' ) then
           !
           !.md<1><b>TYPES</b> [, ALL= char] [ELEMENTS, BOUNDARIES]
           !.md<2>int int                                                                     $ Element, type number
           !.md<2>...
           !.md<1><b>END_TYPES</b>
           !.md<field>TYPES
           !.md<com>This field contains the element types. At each line, the first figure is the element number and the second
           !.md<com>one the element type. If all the elements are the same (say TET04), the following option can be added to the header:
           !.md<com>`TYPES OF ELEMENTS, ALL=TET04` and then the list should be empty.
           !.md<com>The correspondence between element type and element number is the following:
           !.md<com>    -  1D elements:
           !.md<com>        -  BAR02 =    2
           !.md<com>        -  BAR03 =    3
           !.md<com>        -  BAR04 =    4
           !.md<com>    -  2D elements:
           !.md<com>        -  TRI03 =   10
           !.md<com>        -  TRI06 =   11
           !.md<com>        -  QUA04 =   12
           !.md<com>        -  QUA08 =   13
           !.md<com>        -  QUA09 =   14
           !.md<com>        -  QUA16 =   15
           !.md<com>    -  3D elements:
           !.md<com>        -  TET04 =   30
           !.md<com>        -  TET10 =   31
           !.md<com>        -  PYR05 =   32
           !.md<com>        -  PYR14 =   33
           !.md<com>        -  PEN06 =   34
           !.md<com>        -  PEN15 =   35
           !.md<com>        -  PEN18 =   36
           !.md<com>        -  HEX08 =   37
           !.md<com>        -  HEX20 =   38
           !.md<com>        -  HEX27 =   39
           !.md<com>        -  HEX64 =   40
           !.md<com>    -  3D 2D-elements:
           !.md<com>        -  SHELL =   50

           call messages_live('READ ELEMENT  TYPE')
           if( utype > 0 ) then
              do while( words(1) /= 'ENDTY' )
                 call ecoute('reageo')
              end do
           else
              if( exists('BOUND') ) then
                 call read_domain_arrays_types(kfl_icgns,nboun,ktypb,ltypb,lexib)
              else
                 call start_timer('LTYPE', 'READ ', 'ASCII', 1_ip, 1_ip)
                 call read_domain_arrays_types(kfl_icgns,nelem,ktype,ltype,lexis)
                 call end_timer()
              end if
           end if

        else if(words(1) == 'ELEME') then
           !
           !.md<1><b>ELEMENTS</b>
           !.md<2>int int int ...                                                             $ Element, node1, node2, node3 ...
           !.md<2>...
           !.md<1><b>END_ELEMENTS</b>
           !.md<field>ELEMENTS
           !.md<com>This field contains the element connectivity. The first figure is the element number,
           !.md<com>followed by the list of it nodes.
           !
           call start_timer('LNODS', 'READ ', 'ASCII', 1_ip, 1_ip)
           call messages_live('READ ELEMENT  CONNECTIVITY')

           if( ktype == 0 .or. exists('NONOR') ) then
              call ecoute('reageo')
              do while( words(1) /= 'ENDEL' )
                 ielem = int(param(1_ip))
                 if( ielem < 0 .or. ielem > nelem ) &
                      call runend('REAGEO: WRONG NUMBER OF ELEMENT '&
                      //adjustl(trim(intost(ielem))))
                 knode        = nnpar-1
                 ielty        = elmgeo_element_type(ndime,knode)
                 ltype(ielem) = ielty
                 do inode = 1,nnode(ielty)
                    lnods(inode,ielem) = int(param(inode+1))
                 end do
                 call ecoute('reageo')
              end do
           else
              if( associated(lelch) ) then
                 do ielem = 1,nelem
                    ielty = ltype(ielem)
                    if( ielty > 0 .and. lelch(ielem) == ELFEM ) then
                       read(nunit,*,err=2) dummi,(lnods(inode,ielem),inode=1,nnode(ielty))
                    else
                       read(nunit,*,err=2) dummi,(lnods(inode,ielem),inode=1,mnode)
                    end if
                 end do
                 call ecoute('reageo')                 
              else
                 do ielem = 1,nelem
                    ielty = ltype(ielem)
                    if( ielty > 0 ) then
                       read(nunit,*,err=2) dummi,(lnods(inode,ielem),inode=1,nnode(ielty))
                    else
                       read(nunit,*,err=2) dummi,(lnods(inode,ielem),inode=1,mnode)
                    end if
                 end do
                 call ecoute('reageo')
              end if
           end if
           if( words(1) /= 'ENDEL' )&
                call runend('REAGEO: WRONG ELEMENT FIELD')
           call end_timer()

        else if( words(1) == 'CHARA' ) then
           !
           !.md<1>CHARACTERISTICS [,ELEMENTS/BOUNDARIES/NODES]
           !.md<2>int int                                                                     $ Element/Boundary/Nodes, characteristic number
           !.md<2>...
           !.md<1>END_CHARACTERISTICS
           !.md<field>CHARACTERISTICS
           !.md<com>This field contains the element or boundary charcteristics. The first figure is the element/boundary number followed by
           !.md<com>the element/boundary characteristic=0,1,2 for finite element/boundary/node, extension element/boundary/node (Chimera),
           !.md<com>hole element/node (Chimera, Embedded mesh).
           !
           if( exists('BOUND') ) then
              call messages_live('READ BOUNDARY CHARACTERISTIC')
              call ecoute('reageo')
              do while( words(1) /= 'ENDCH' )
                 iboun        = int(param(1))
                 lboch(iboun) = int(param(2))
                 call ecoute('reageo')
              end do
           else if( exists('NODES') ) then
              call messages_live('READ NODE CHARACTERISTIC')
              call ecoute('reageo')
              do while( words(1) /= 'ENDCH' )
                 ipoin        = int(param(1))
                 lnoch(ipoin) = int(param(2))
                 call ecoute('reageo')
              end do
           else
              call messages_live('READ ELEMENT CHARACTERISTIC')
              call ecoute('reageo')
              do while( words(1) /= 'ENDCH' )
                 ielem        = int(param(1))
                 lelch(ielem) = int(param(2))
                 call ecoute('reageo')
              end do
           end if

        else if(words(1) == 'LELBO' ) then
           !
           !.md<1>LELBO                 $ Boundary, element number
           !.md<2>int int
           !.md<2>...
           !.md<1>END_LELBO
           !.md<field>LELBO
           !.md<com>List of boundary elements
           !
           call messages_live('READ BOUNDARY ELEMENT')
           call ecoute('reageo')
           iboun = 0
           do while( words(1) /= 'ENDLE' )
              iboun        = int(param(1))
              lelbo(iboun) = int(param(2))
              call ecoute('reageo')
           end do
           if( iboun == nboun ) klelb = 1

        else if(words(1) == 'BOUND' .and. kfl_autbo /= 1 ) then
           !
           !.md<1>BOUNDARIES [, ELEMENT]
           !.md<2>int int... [int]                                                            $ Boundary, node1, node2, node3 [element number]
           !.md<2>...
           !.md<1>END_BOUNDARIES
           !.md<field>BOUNDARIES
           !.md<com>This field contains the boundary connectivity. The first figure is the boundary number followed by its nodes.
           !.md<com>If option `ELEMENT` is present in the header, the last parameter should be the element to which the boundary belongs.
           !.md<com>_Note: this option is now deprecated and `LELBO` should be prefered_.
           !.md<com>It should be given if possible to avoid Alya compute it.
           !
           call start_timer('LNODB', 'READ ', 'ASCII', 1_ip, 1_ip)
           if(exists('ELEME')) kfl_bouel=1
           ii = 1
           if(exists('NOELE')) ii = 2
           if(exists('ELEME')) ii = 2

           call ecoute('reageo')
           if( kfl_bouel == 0 ) then
              call messages_live('READ BOUNDARY CONNECTIVITY AND TYPE')
              do while( words(1) /= 'ENDBO' )
                 iboun = int(param(1))
                 if( iboun < 0 .or. iboun > nboun ) call runend('REAGEO: WRONG NUMBER OF BOUNDARIES')
                 if( ktypb == 0 ) then
                    knode = nnpar-ii
                    iblty = elmgeo_element_type(ndimb,knode)
                    ltypb(iboun) = iblty
                 else
                    iblty = ltypb(iboun)
                 end if
                 do inodb = 1,nnode(iblty)
                    lnodb(inodb,iboun) = int(param(inodb+1))
                 end do
                 call ecoute('reageo')
              end do

              call end_timer()
              
           else
              if(      klelb == 0 .and. ktypb == 0 ) then
                 call messages_live('READ BOUNDARY CONNECTIVITY, TYPE AND ELEMENT')
              else if( klelb == 1 .and. ktypb == 0 ) then
                 call messages_live('READ BOUNDARY CONNECTIVITY AND ELEMENT')
              else if( klelb == 0 .and. ktypb == 1 ) then
                 call messages_live('READ BOUNDARY CONNECTIVITY AND TYPE')
              else if( klelb == 1 .and. ktypb == 1 ) then
                 call messages_live('READ BOUNDARY CONNECTIVITY')
              end if
              do while( words(1) /= 'ENDBO' )
                 iboun = int(param(1))
                 if( iboun < 0 .or. iboun > nboun ) call runend('REAGEO: WRONG NUMBER OF BOUNDARIES')
                 if( ktypb == 0 ) then
                    knode = nnpar-ii
                    iblty = elmgeo_element_type(ndimb,knode)
                    ltypb(iboun) = iblty
                 else
                    iblty = ltypb(iboun)
                 end if
                 do inodb = 1,nnode(iblty)
                    lnodb(inodb,iboun) = int(param(inodb+1))
                 end do
                 if( klelb == 0 ) then
                    lelbo(iboun) = int(param(nnpar))
                 end if
                 call ecoute('reageo')
              end do
              call end_timer()
              
           end if

        else if( words(1) == 'SUBDO' ) then
           !
           !.md<1>SUBDOMAIN, [DEFAULT= int]
           !.md<2>int int                                                                     $ Element, subdomain
           !.md<2>...
           !.md<1>END_SUBDOMAIN
           !.md<field>SUBDOMAIN
           !.md<com>List of element subdomain used for the HERMESH method.
           !
           isubd = 1
           if( exists('DEFAU') ) isubd = getint('DEFAU',1_ip,'#DEFAULT SUBDOMAIN')
           do ielem = 1,nelem
              lesub(ielem) = isubd
           end do
           call messages_live('READ ELEMENT SUBDOMAIN')

           call ecoute('reageo')
           do while( words(1) /= 'ENDSU' )
              ielem = int(param(1_ip))
              lesub(ielem) = int(param(2))
              call ecoute('reageo')
           end do

        else if(words(1) == 'COORD' ) then

           call start_timer('COORD', 'READ ', 'ASCII', 1_ip, 1_ip)
           !
           !.md<1><b>COORDINATES</b> [NOT_SORTED]
           !.md<2>int real real real      $ id_node, coor_x, coor_y, coor_z
           !.md<2>...
           !.md<1><b>END_COORDINATES</b>
           !.md<field>COORDINATES
           !.md<com>This field contains the coordinates of the nodes. The first figure is the node number
           !.md<com>followed by its coordinates.
           !.md<com>"id_node" is not taken into account, unless the option `NOT_SORTED` is present, which means
           !.md<com>that nodes are not sorted by id_node.
           !.md<com>_Note that this option is strongly discouraged since it could prevent external tools from working efficiently._
           !
           call messages_live('READ COORDINATE')
           if (words(2) == 'NOTSO') then
              do ipoin=1,npoin
                 read(nunit,*,err=3) jpoin,(coord(idime,jpoin),idime=1,ndime)
              end do
           else
              do ipoin=1,npoin
                 read(nunit,*,err=3) dummi,(coord(idime,ipoin),idime=1,ndime)
              end do
           end if
           call ecoute('reageo')
           if(words(1)/='ENDCO')&
                call runend('REAGEO: WRONG COORDINATE FIELD')
           call end_timer()

        else if( words(1) == 'PERIO' ) then
           !
           call runend('REAGEO: PERIODIC NO LONGRER ACCEPTED')

        else if( words(1) == 'LMAST' ) then
           !
           !.md<1>LMAST
           !.md<2>int int                                                                     $ Slave, master
           !.md<2>...
           !.md<1>END_LMAST
           !.md<field>LMAST
           !.md<com>This field contains the list of slave (first figure) and master (second figure) periodic nodes.
           !.md<com>It should be noted that the results will not be exactly the same in sequential and in parallel
           !.md<com>due to the different treaments of periodicity in both cases. The difference is in the way scalar
           !.md<com>product are computed: they involve both the slaves and master in sequential but only the master
           !.md<com>in parallel.
           !
           call messages_live('READ PERIODICITY (LMAST)')
           call ecoute('reageo')
           do while( words(1) /= 'ENDLM' )
              if( nnpar > 0 ) then
                 ipoin_slave        = int(param(1),ip)
                 ipoin_master       = int(param(2),ip)
                 lmast(ipoin_slave) = ipoin_master
              end if
              call ecoute('reageo')
           end do

        else if( words(1) == 'MATER' .and. words(2) /= 'OFF  ' ) then
           !
           !.md<1>MATERIALS, NUMBER= int, [DEFAULT= int]
           !.md<2>int int                                                                     $ Element, material number
           !.md<2>...
           !.md<1>END_MATERIALS
           !.md<field>MATERIALS
           !.md<com>This field defines elemental materials. The first figure is the element and the second the material mumber.
           !.md<com>The header option NUMBER= int is the total number of materials. The header option DEFAULT= imate is the default
           !.md<com>material assigned to elements when they do not appear in the list.
           !
           imate = 0
           if( exists('DEFAU') ) imate = getint('DEFAU',1_ip,'#DEFAULT MATERIAL')
           call messages_live('READ MATERIAL')
           call ecoute('reageo')
           if( imate > 0 ) then
              do ielem = 1,nelem
                 lmate(ielem) = imate
              end do
           end if
           do while( words(1) /= 'ENDMA' )
              ielem = int(param(1_ip))
              lmate(ielem) = int(param(2))
              call ecoute('reageo')
           end do

        else if( words(1) == 'FIELD' ) then

           call runend('REAGEO: FIELD ARE NO LONGER DECLARED HERE')

        else if( words(1) == 'GROUP' ) then

           call runend('REAGEO_SEQ; GROUPS SHOULD BE DECLARED IN STRATEGY FIELD')

        else if( words(1) == 'CRVDA' ) then

           call runend('REAGEO_SEQ: CURVATURE DATA SHOULD BE DECLARED IN STRATEGY FIELD')

        else if( words(1) == 'CRVGE' ) then

           call runend('REAGEO_SEQ: CURVATURE DATA SHOULD BE DECLARED IN STRATEGY FIELD')

        end if

     end do
     !
     !.md<0><b>END_GEOMETRY</b>
     !.md</code>
     !.md<>
  end if

  return

1 call runend('REAGEO: WRONG NUMBER OF NODES FOR ELEMENT '//trim(intost(ielem)))
2 call runend('REAGEO: WRONG ELEMENT CONNECTIVITY FOR ELEMENT '//trim(intost(ielem)))
3 call runend('REAGEO: WRONG COORDINATES FOR NODE '//trim(intost(ipoin)))
4 call runend('REAGEO: WRONG NUMBER OF FIELD COMPONENTS '//trim(intost(ipoin)))
5 call runend('REAGEO: WRONG NUMBER OF FIELD COMPONENTS when reading curvatureField '//trim(intost(ipoin)))
6 call runend('REAGEO: WRONG NUMBER OF FIELD COMPONENTS when reading curvatureField  '//trim(intost(ipoin)))
7 call runend('REAGEO: WRONG NUMBER OF FIELD COMPONENTS when reading a field with many components '//trim(intost(ipoin)))

end subroutine reageo_seq
