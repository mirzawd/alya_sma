!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup DomainInput
!> @ingroup    Domain
!> @{
!> @file    reabcs.f90
!> @author  Guillaume Houzeaux
!> @date    18/09/2012
!> @brief   Read boundary conditions
!> @details Read boundary conditions:
!!          \verbatim
!!          - Scalars:
!!            - KFL_ICODN ................. # codes on nodes (computed later on; here 1 or 0)
!!            - KFL_ICODB ................. # codes on boundaries (computed later on; here 1 or 0)
!!          - Arrays:
!!            - KFL_CODNO(MCONO,NPOIN) .... Codes on nodes (can be multiple if extrapolated from KFL_CODBO)
!!            - KFL_CODBO(NBOUN) .......... Codes on boundaries
!!          \endverbatim
!> @}
!-----------------------------------------------------------------------
module mod_reabcs

  use mod_messages,      only : messages_live
  use mod_domain,        only : domain_memory_allocate

  implicit none
  private

  public :: reabcs_seq
  public :: reabcs_init
  public :: reabcs_geometrical

contains

  subroutine reabcs_init()
    use def_kintyp
    use def_parame
    use def_master
    use def_domain
    use def_inpout
    implicit none

    if( ISEQUEN .or. ( IMASTER .and. kfl_ptask /= 2 ) ) then
       !
       ! Initializations and allocations
       !
       kfl_icodn =  0
       kfl_icodb =  0
    end if

  end subroutine reabcs_init

  subroutine reabcs_seq()
    use def_kintyp
    use def_parame
    use def_master
    use def_domain
    use def_inpout
    use mod_ecoute, only :  ecoute
    implicit none
    integer(ip)          :: ipoin,iboun,icode,pnodb,knodb(60)
    integer(ip)          :: icond,icobd
    integer(ip)          :: icono,ivcod,ivcob,iposi
    integer(ip)          :: iread
    integer(ip), pointer :: lpbcs(:)
    integer(ip), pointer :: lnodb_tmp(:,:)
    character(300)       :: messa
    logical              :: boundary_processed

    nullify(lpbcs)
    nullify(lnodb_tmp)
    boundary_processed = .FALSE.

    call reabcs_init()

    if( ISEQUEN .or. ( IMASTER .and. kfl_ptask /= 2 ) ) then

       icobd     = -99999
       !
       ! Reach section
       !

       call ecoute('REABCS')
       do while( words(1) /= 'BOUND' )
          call ecoute('REABCS')
       end do
       !
       ! Check if boundary codes are extrapolated to nodes
       !
       if( exists('EXTRA') ) then
          call runend("REABCS: YOU ARE STILL USING THE OLD FORMAT! EXTRAPOLATE_BOUNDARY_CONDITIONS MUST BE DEFINED IN THE STRATEGY SECTION")
       end if

       if (kfl_extra == 0) then
          messa = 'KERNEL WILL NOT EXTRAPOLATE FROM BOUNDARY CODES TO NODE CODES'
       else if (kfl_extra == 1) then
          messa = 'KERNEL WILL EXTRAPOLATE FROM BOUNDARY CODES TO NODE CODES'
       end if
       call messages_live(messa,'WARNING')
       !
       ! Initialization of local variables
       !
       ivcod = 0
       ivcob = 0
       iread = 0
       !
       ! Code permutation
       !
       allocate( lpbcs(-mcodb:mcodb) )
       do icode = -mcodb,mcodb
          lpbcs(icode) = icode
       end do
       !.md<module>kernel
       !.md<input>case.dom.dat
       !.md<pos>4
       !.md<sec>
       !.md<0># Boundary conditions
       !.md<>
       !.md<code>
       !.md<0><b>BOUNDARY_CONDITIONS [, EXTRAPOLATE]</b>
       !.md<field>BOUNDARY_CONDITIONS
       !.md<com>In this field, codes are assigned to nodes or boundaries. These codes
       !.md<com>are generic, just numbers, and do not have any physical interpretation.
       !.md<com>Each module should interpret these codes in the BOUNDARY_CONDITIONS field
       !.md<com>of its data file in terms of Dirichlet or Neumann conditions.<p>
       !.md<com>The EXTRAPOLATE option should be used to extrapolate boundary codes imposed on boundaries to nodes.
       !.md<com>If this option is used, a node can have multiple codes. For example, if a corner node belongs
       !.md<com>to 3 boundary elements with code 4,7 and 3, respectively, then the node inherits the three codes.
       !.md<com>Extrapolated boundary conditions from boundaries to nodes should be preferred if one plans
       !.md<com>to use the mesh multiplication algorithm: this is in fact the only way to uniquely define the
       !.md<com>boundary conditions on the refined meshes.
       !.md<com>Note that the ON_NODES conditions overwrites the extrapolated ones.
       !
       do while( words(1) /= 'ENDBO' )
          call ecoute('reabcs')

          if( words(1) == 'PERMU' ) then

             call runend('PERMUTATION OF BOUNDARY CODE NO LONGER EXISTS')

          else if( words(1) == 'ONNOD' ) then

             !-------------------------------------------------------------
             !
             ! Node codes
             !
             !-------------------------------------------------------------
             !
             !.md<1>ON_NODES
             !.md<2>...
             !.md<2>int1 int2                                                                   $ Node, code
             !.md<2>...
             !.md<1>END_ON_NODES
             !.md<field>ON_NODES
             !.md<com>This field contains the list of node codes. int1 is the node and int2 the
             !.md<com>associated code.
             !
             call messages_live('READ BOUNDARY CODES ON NODES')

             boundary_processed = .TRUE.

             if( exists('POSIT') ) then
                iposi = 1
             else
                iposi = 0
             end if

             if( exists('DEFAU') ) then
                icond = getint('DEFAU',0_ip,'*DEFAULT CODE ON NODES')
                if(icond>mcodb.or.icond<-mcodb)&
                     call runend('REABCS: WRONG DEFAULT NODE CODE')
             else
                icond = -99999
             end if

             call ecoute('reabcs')
             kfl_icodn = 1
             if( iread == 0 ) then ! Could have been allocated by mesh generator
                iread = 1
                call domain_memory_allocate('KFL_CODNO') ! KFL_CODNO

                if(icond==-99999) then
                   do ipoin = 1,npoin
                      do icono = 1,mcono
                         kfl_codno(icono,ipoin) = mcodb+1
                      end do
                   end do
                else
                   do ipoin = 1,npoin
                      kfl_codno(1,ipoin) = icond
                   end do
                end if
             else
                do ipoin = 1,npoin
                   icode = kfl_codno(1,ipoin)
                end do
             end if

             if( iposi == 1 ) then
                !
                ! Positional
                !
                call poscod(1_ip)

             else
                !
                ! Given code
                !
                do while( words(1) /= 'ENDON' )
                   ipoin = int(param(1),ip)
                   icode = int(param(2),ip)
                   kfl_codno(1,ipoin) = icode
                   call ecoute('reabcs')
                end do

             end if

          else if( words(1) == 'ONBOU' ) then

             !-------------------------------------------------------------
             !
             ! Boundary codes
             !
             !-------------------------------------------------------------
             !
             !.md<1>ON_BOUNDARIES [, DEFAULT= int3]
             !.md<2>...
             !.md<2>int1 int2                                                                   $ Node, code
             !.md<2>...
             !.md<1>END_ON_BOUNDARIES
             !.md<field>BOUNDARIES
             !.md<com>This field contains the list of boundary codes. int1 is the boundary number and int2 the associated code.
             !.md<com>The default option enables to impose automatically the code int3 to all the boundaries.
             !
             call messages_live('READ BOUNDARY CODES ON BOUNDARIES')

             boundary_processed = .TRUE.

             if( exists('POSIT') ) then
                iposi = 1
             else
                iposi = 0
             end if

             if( exists('DEFAU') ) then
                icond = getint('DEFAU',0_ip,'*DEFAULT CODE ON BOUNDARIES')
                if( icond > mcodb .or. icond < -mcodb ) &
                     call runend('REABCS: WRONG DEFAULT BOUNDARIES CODE')
             else
                icond = mcodb + 1
             end if

             if( words(2) == 'UNKNO' ) then
                call ecoute('reabcs')

                allocate(lnodb_tmp(mnodb,nboun))
                do iboun = 1,nboun
                   pnodb = nnode(ltypb(iboun))
                   lnodb_tmp(1:mnodb,iboun) = lnodb(1:mnodb,iboun)
                   call heapsorti1(2_ip,pnodb,lnodb_tmp(:,iboun))
                end do
                if( words(1) /= 'ENDON' ) then
                   kfl_icodb = 1
                   if( .not. associated(kfl_codbo) ) then
                      call domain_memory_allocate('KFL_CODBO')
                      do iboun = 1,nboun
                         kfl_codbo(iboun) = icond
                      end do
                   end if
                   if( iposi == 1 ) then
                      !
                      ! Positional
                      !
                      !call poscod(2_ip)  !falta crear
                      call runend('POSCOD: POSITIONAL makes no sense for  UNKNO')

                   else
                      !
                      ! Given code
                      !
                      do while( words(1) /= 'ENDON')
                         pnodb = int(param(2))
                         knodb(1:pnodb) = int(param(3:2+pnodb),ip)
                         call heapsorti1(2_ip,pnodb,knodb)
                         call finbou_fast(pnodb,knodb,iboun,lnodb_tmp)
                         !call finbou(pnodb,knodb,iboun)
                         if( iboun == 0 .or. iboun > nboun ) then
                            icode = int(param(3+pnodb),ip)
                            print*,'merde=',pnodb,param
                            print*,'BOUNDARY= ',knodb(1:pnodb), 'LAST CODE=',icode
                            call runend('REABCS: BOUNDARY CODES HAVE BEEN IMPOSED ON A NON-EXISTING BOUNDARY')
                         else
                            icode = int(param(3+pnodb),ip)
                            kfl_codbo(iboun) = icode
                         end if
                         call ecoute('reabcs')
                      end do
                   end if
                end if
                deallocate(lnodb_tmp)

             else
                call ecoute('reabcs')
                if( words(1) /= 'ENDON' ) then
                   kfl_icodb = 1
                   if( .not. associated(kfl_codbo) ) then
                      call domain_memory_allocate('KFL_CODBO')
                      do iboun = 1,nboun
                         kfl_codbo(iboun) = mcodb + 1
                      end do
                   end if
                   if( iposi == 1 ) then
                      !
                      ! Positional
                      !
                      call poscod(2_ip)
                   else
                      !
                      ! Given code
                      !
                      do while( words(1) /= 'ENDON' )
                         iboun = int(param(1),ip)
                         icode = int(param(2),ip)
                         if (iboun == 0) then
                            messa = 'KERNEL FOUND NO DIRECT CONDITIONS ON_NODES'
                            call messages_live(messa,'WARNING')
                         else
                            kfl_codbo(iboun) = icode
                         end if
                         call ecoute('reabcs')
                      end do
                   end if
                end if
             end if

          else if( words(1) == 'MATER' ) then

             call runend('MATERIALS FROM BOUNDARIES ARE NOW IN STRATEGY FIELD')

          else if( words(1) == 'GEOME' ) then

             call runend('GEOMETRICAL BOUNDARY CONDITIONS ARE NOW IN STRATEGY FIELD')
          end if

       end do

       if( .NOT. boundary_processed ) then
         call messages_live('Empty BOUNDARIES...END_BOUNDARIES in dom.dat. Likely missing ONNOD or ONBOU.','WARNING')
       end if
     !
     !.md<0><b>END_BOUNDARY_CONDITIONS</b>
     !.md</code>
     !.md<>

    end if

  end subroutine reabcs_seq

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-10-24
  !> @brief   Read geometrical bcs
  !> @details Read data for geometrical boundary conditions
  !>
  !-----------------------------------------------------------------------

  subroutine reabcs_geometrical()

    use def_kintyp
    use def_parame
    use def_master
    use def_inpout
    use def_domain
    use mod_ecoute

    integer(ip) :: icond,ipara

    kfl_geome = 1

    do while( words(1) /= 'ENDGE' )

       icond = 0
       if(      words(1) == 'PRESC' .or. words(1) == 'INFLO' ) then
          icond = 1
       else if( words(1) == 'FREES' ) then
          icond = 2
       else if( words(1) == 'WALLL' ) then
          icond = 3
       else if( words(1) == 'SYMME' .or.  words(1) == 'SLIPW' ) then
          icond = 4
       else if( words(1) == 'NOSLI' .or. words(1) == 'WALL ' ) then
          icond = 5
       else if( words(1) == 'FREE ' .or. words(1) == 'OUTFL' ) then
          icond = 6

       else if( words(1) == 'ANGLE' ) then
          geoan = getrea('GEOAN',45.0_rp,'#Geometric angle')

       else if( words(1) == 'CONVE' ) then
          if( words(2) == 'FREE ') then
             kfl_convx = 0
          else if(words(2) == 'EXNOR') then
             kfl_convx = 1
          else if(words(2) == 'FIXED') then
             kfl_convx = 2
          else if(words(2) == 'GEOME') then
             kfl_convx = 3
          end if

       else if( words(1) == 'CRITE' ) then
          if( words(2) == 'WINDA') then
             kfl_frees = 0
             awind = getrea('WINDA',0.0_rp,'#Wind angle')
             awind = awind / 180.0_rp * pi

          else if(words(2) == 'VALUE') then
             kfl_frees = getint('VALUE',1_ip,'*FUNCTION NUMBER')
          end if

          if( exists('TOLER') ) then
             tolan = getrea('TOLER',5.0_rp,'#Tolerance for determining inflow')
             tolan = tolan / 180.0_rp * pi
          end if

       end if
       if( icond /= 0 ) then
          npbcs(icond) = nnpar
          if( nnpar > size(lsbcs,1,KIND=ip) ) call runend('REABCS: TOO MANY PARAMETERS')
          do ipara = 1, nnpar
             lsbcs(ipara,icond) = int(param(ipara),ip)
          end do
       end if
       call ecoute('reabcs')

    end do

  end subroutine reabcs_geometrical

end module mod_reabcs
