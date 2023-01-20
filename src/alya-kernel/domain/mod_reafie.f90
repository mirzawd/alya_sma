!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup DomainInput
!> @{
!> @file    reafie.f90
!> @author  Guillaume Houzeaux
!> @date    15/02/2018
!> @brief   Read fields
!> @details Read fields
!>          ITASK == 0 ........ Old way, FIELDS are in GEOMETRY section
!>                == 1 ........ New way, FIELDS are in a separate section
!> @}
!-----------------------------------------------------------------------
module mod_reafie

  use mod_messages,      only : livinf
  use mod_domain,        only : domain_memory_allocate
  use mod_mpio_config,   only : mpio_config
  implicit none
  private

  public :: reafie_seq
  public :: reafie_read_header

contains

  subroutine reafie_init()

    use def_kintyp
    use def_domain
    implicit none

    call livinf(0_ip,'READ FIELDS',0_ip)
    !
    ! Initializations and allocations
    !
    curvatureDataField = 0                   ! curvature data field id
    curvatureField     = 0                   ! curvature data field id

  end subroutine reafie_init

  !-----------------------------------------------------------------------
  !>
  !> @author  Deamien Dosimont
  !> @date    2018-10-24
  !> @brief   Read field dimensions
  !> @details Read field dimensions
  !>
  !-----------------------------------------------------------------------

  subroutine reafie_read_dimensions(ifiel,kfl_defau,defau,ONLY_HEADER)

    use def_inpout
    use mod_ecoute,   only : ecoute
    integer(ip), intent(out)           :: ifiel
    integer(ip), intent(out), optional :: kfl_defau
    real(rp),    intent(out), optional :: defau
    logical(lg), intent(in),  optional :: ONLY_HEADER
    logical(lg)                        :: if_only_header

    if( present(ONLY_HEADER) ) then
       if_only_header = ONLY_HEADER
    else
       if_only_header = .false.
    end if

    call ecoute('reafi')
    do while( words(1) /= 'ENDDI' )
       if( words(1) == 'FIELD' ) then
          call reafie_read_header(ifiel,kfl_defau,defau)
          if( .not. if_only_header ) call reafie_allocate(ifiel,kfl_defau,defau)
       end if
       call ecoute('reafi')
    end do

  end subroutine reafie_read_dimensions

  !-----------------------------------------------------------------------
  !>
  !> @author  Deamien Dosimont
  !> @date    2018-10-24
  !> @brief   Read field
  !> @details Read field
  !>
  !-----------------------------------------------------------------------

  subroutine reafie_read_field(ifiel,kfl_defau,defau)
    use def_kintyp
    use def_parame
    use def_master
    use def_domain
    use def_inpout
    use mod_ecoute,   only : ecoute
    use mod_opfpos,   only : postpr_intto8

    implicit none

    integer(ip), intent(in)    :: ifiel,kfl_defau
    real(rp),    intent(inout) :: defau
    integer(ip)                :: ipoin,idime,jstep
    !
    ! Read field
    !
    !-If there are more than XX numbers in this field read different (ecoute only reads a maximum of maxwp characters per line)
    if(kfl_field(1,ifiel) < 6_ip) then

       call ecoute('reafi')
       if( kfl_field(4,ifiel) == 1_ip ) then ! 1 Step - default behaviour
          do while( words(1) /= 'ENDFI' )
             ipoin = int(param(1))
             if (ipoin==0) then
                 call runend( "Point/node id==0 in field "//postpr_intto8(ifiel) )
             end if

             do jstep = 1,kfl_field(4,ifiel)
                do idime = 1,kfl_field(1,ifiel)
                   xfiel(ifiel) % a(idime,ipoin,jstep) = param(idime+1)
                end do
             end do
             call ecoute('reafi')
          end do

       else

          do while( words(1) /= 'ENDFI' )
             !
             ! If there is more than 1 step the syntax is
             ! FIELD
             !  STEP 1 $, TIME=0.0, $time is defined now in the DIMENSIONS section of dom.dat
             !   7   3.5
             !   9   4.0
             !  END_STEP
             !  STEP 2 $, TIME=3600.0
             !   7   5.0
             !   9   6.6
             !  END_STEP
             ! END_FIELD
             !
             if( words(1) == 'STEP ' ) then
                jstep = nint(param(1))
                if (jstep > kfl_field(4,ifiel)) call runend('REAFIE: jstep > kfl_field(4,ifiel)  - Check field steps')
                call ecoute('reafi')
                do while( words(1) /= 'ENDST' )
                   if( words(1) == 'STEP ' ) then 
                      call runend('Step '//postpr_intto8(ifiel)//' is missing END_STEP.')
                   end if

                   ipoin = int(param(1))
                   do idime = 1,kfl_field(1,ifiel)
                      xfiel(ifiel) % a(idime,ipoin,jstep) = param(idime+1)
                   end do
               
                   call ecoute('reafi')
                end do

             end if
             call ecoute('reafi')
          end do
       end if    ! 1 or more steps
    else
       !-Reading fields with many components (that do not fit in 60 characters)
       !if((multiply_with_curvature.eq.1).and.(ifiel.eq.curvatureField)) then
       !if( (ifiel.gt.0).and.(ifiel.eq.curvatureField) ) then
       if( (ifiel == curvatureField) ) then
          do ipoin=1,size(xfiel(ifiel) % a,2,KIND=ip)
             if(xfiel(curvatureDataField) % a(1,ipoin,1) == 0) then
                xfiel(ifiel) % a(:,ipoin,1) = 0_rp
                read(nunit,*,err=5) defau,xfiel(ifiel) % a(1,ipoin,1)
             else
                read(nunit,*,err=6) defau,xfiel(ifiel) % a(:,ipoin,1)
             end if
          end do
       else! if(ifiel.gt.0) then
          do ipoin=1,size(xfiel(ifiel) % a,2,KIND=ip)
             read(nunit,*,err=7) defau,xfiel(ifiel) % a(:,ipoin,1)
          end do
       end if
       call ecoute('reafi')
       if(words(1)/='ENDFI') then
          call runend('REAFI: WRONG FIELD')
       end if
    end if

    return

5   call runend('REAFI: WRONG NUMBER OF FIELD COMPONENTS when reading curvatureField '//trim(intost(ipoin)))
6   call runend('REAFI: WRONG NUMBER OF FIELD COMPONENTS when reading curvatureField  '//trim(intost(ipoin)))
7   call runend('REAFI: WRONG NUMBER OF FIELD COMPONENTS  when reading a field with many components '//trim(intost(ipoin)))

  end subroutine reafie_read_field

  subroutine reafie_seq()
    use def_kintyp
    use def_parame
    use def_master
    use def_domain
    use def_inpout
    use mod_messages, only : messages_live
    use mod_ecoute,   only : ecoute
    implicit none
    integer(ip)             :: kfl_defau,ifiel
    real(rp)                :: defau
    logical(lg)             :: already_read

    if( ISEQUEN .or. ( IMASTER .and. kfl_ptask /= 2 ) ) then
       !
       ! Reach section
       !

       call ecoute('reafie',STOP_END_OF_FILE=.false., DO_NOT_READ_INCLUDE=.true.)
       if( words(1) == 'EOF  ' ) return
       do while( words(1) /= 'FIELD' )
          call ecoute('reafie',STOP_END_OF_FILE=.false., DO_NOT_READ_INCLUDE=.true.)
          if( words(1) == 'EOF  ' ) return
       end do

       call reafie_init()

       already_read = .false.

       if( words(1) == 'FIELD' ) then
          !
          ! ADOC[1]> FIELDS, NUMBER= int
          ! ADOC[2]> FIELD= int, DIMENSION= int, NODE/ELEMENT/BOUNDARY, DEFAULT=real, EXCLUDE_NEGATIVE_VALUES , STEPS= int
          ! ADOC[2]>   int real...                                                              $ Element/node/boundary, value1, value2...
          ! ADOC[2]> END_FIELD
          ! ADOC[1]> END_FIELDS
          ! ADOC[d]> The field FIELDS deserves special attention. It declares fields of 1 or more dimensions defined on nodes or elements.
          ! ADOC[d]> These fields can be used in the modules to assign nodeal/elemental values to some quantities.
          ! ADOC[d]> For example, in the case of SOLIDZ module, a field can represent some material directions.
          ! ADOC[d]> FIELDS, NUMBER= int, int is the total number of fields to be declared.
          ! ADOC[d]> FIELD= ifiel, DIMENSION= ndimf, NODE/ELEMENT/BOUNDARY. Field ifiel is now going to be defined. The values are ndimf dimensional and the
          ! ADOC[d]> fields is defined on node, element or boundary. The list of values has the form: ii, value(1,ii)...value(ndimf,ii) where
          ! ADOC[d]> ii= ipoin for values on nodes and ii= ielem for values on elements.
          ! ADOC[d]> The following section is outdated. Multiple steps are specified as above in the source code using STEP/END_STEP.
          ! ADOC[d]> Now it can also deal with several steps. If there are more than 1 step they are stored they are stored :
          ! ADOC[d]> value(1 of step 1,ii)...value(ndimf of step 1,ii),value(1 of step 2,ii)...value(ndimf of step 2,ii),............
          ! ADOC[d]> DEFAULT: Impose a default value of real. EXCLUDE_NEGATIVE_VALUES can be useful when using the mesh multiplciation
          ! ADOC[d]> not to take these values nto account when interpolating
          !
          ! FIELDS: NFIEL, KFL_FIELD, XFIEL
          !
          if( words(2) == 'NUMBE' ) nfiel = getint('NUMBE',0_ip,'#NUMBER OF FIELDS')
          call messages_live('READ '//trim(intost(nfiel))//' FIELDS')
          !
          ! Allocate field type XFIEL(NFIEL)
          !
          call ecoute('reafi')
          do while( words(1) /= 'ENDFI' )

             if(      words(1) == 'DIMEN' ) then
                call reafie_read_dimensions(ifiel, kfl_defau, defau)
                already_read = .true.
             else if(words(1) == 'FIELD' ) then
                !
                ! Read header
                !
                if( .not. already_read ) then
                   call reafie_read_header(ifiel,kfl_defau,defau)
                   call reafie_allocate(ifiel,kfl_defau,defau)
                else
                   ifiel = getint('FIELD',0_ip,'#FIELD NUMBER')
                end if
                call reafie_read_field(ifiel, kfl_defau, defau)
             end if   !words(1) == 'FIELD'
             call ecoute('reafi')
          end do
       end if

    end if

  end subroutine reafie_seq

  subroutine reafie_allocate(ifiel,kfl_defau,defau)

    use def_kintyp
    use def_parame
    use def_master
    use def_domain
    use def_inpout

    implicit none

    integer(ip), intent(in)           :: ifiel
    integer(ip), intent(in), optional :: kfl_defau
    real(rp),    intent(in), optional :: defau
    !
    ! Allocate memory
    !
    call domain_memory_allocate('XFIEL % A',NUMBER1=ifiel)
    !
    ! Put default value
    !
    if( present(kfl_defau) .and. present(defau) ) then
       if( kfl_defau == 1 ) then
          xfiel(ifiel) % a = defau
       end if
    end if

  end subroutine reafie_allocate

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-10-24
  !> @brief   Read header of fields
  !> @details Read the header of fields
  !>
  !-----------------------------------------------------------------------

  subroutine reafie_read_header(ifiel,kfl_defau,defau)

    use def_kintyp
    use def_parame
    use def_master
    use def_domain
    use def_inpout

    implicit none
    integer(ip), intent(out)           :: ifiel
    integer(ip), intent(out), optional :: kfl_defau
    real(rp),    intent(out), optional :: defau

    ifiel = getint('FIELD',1_ip,'#FIELD NUMBER')

    if( ifiel < 0 .or. ifiel > nfiel ) call runend('REAFI: WRONG FIELD NUMBER')
    !
    ! KFL_FIELD(1,:): Dimension
    !
    if( exists('DIMEN') ) then
       kfl_field(1,ifiel) = getint('DIMEN',1_ip,'#FIELD DIMENSION')

       !
       ! KFL_FIELD(2,:): Type
       !
       if( exists('ELEME') ) then
          kfl_field(2,ifiel) = NELEM_TYPE
       else if( exists('NODES') .or. exists('NODE ') ) then
          kfl_field(2,ifiel) = NPOIN_TYPE
       else if( exists('BOUND') ) then
          kfl_field(2,ifiel) = NBOUN_TYPE
       else
          call runend('REAFI: DEFINE IF FIELD IS DEFINED ON NODES OR ELEMENTS')
       end if
       !
       ! KFL_FIELD(3,:): Negative value excluded?
       !
       if( exists('EXCLU') ) then
          kfl_field(3,ifiel) =  getint('EXCLU',1_ip,'#FIELD STEPS EXCLUDED')
       else
          kfl_field(3,ifiel) = 0
       end if
       !
       ! KFL_FIELD(4,:): Steps
       !
       if( exists('STEPS') ) then
          kfl_field(4,ifiel) = getint('STEPS',1_ip,'#FIELD DIMENSION')
       else
          kfl_field(4,ifiel) = 1_ip
       end if
       !
       ! KFL_FIELD(6,:): When field should be read
       !
       if( exists('BEGIN') ) then
          kfl_field(6,ifiel) = 0_ip    ! At runtime
       else if( exists('WHENN') .or. exists('ONDEM') ) then
          if (.NOT. mpio_config%enabled) then
             call runend('ONDEMAND fields work only with MPIO ON')
          end if
          kfl_field(6,ifiel) = 1_ip    ! When needed
       end if
       !
       ! KFL_FIELD(7,:): Periodic field
       !
       if( exists('NONPE') ) then
          kfl_field(7,ifiel) = 0_ip    ! Non-periodic
       else if( exists('PERIO') ) then
          kfl_field(7,ifiel) = 1_ip    ! Periodic
       end if
    end if
    !
    ! Default value
    !
    if( present(kfl_defau) .and. present(defau) ) then
       if( exists('DEFAU') ) then
          kfl_defau = 1
          defau     = getrea('DEFAU',1.0_rp,'#Default value')
       else
          kfl_defau = 0
       end if
    end if
    
  end subroutine reafie_read_header

end module mod_reafie
