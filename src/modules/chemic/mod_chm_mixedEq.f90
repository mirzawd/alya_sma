!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_chm_mixedEq

  use def_kintyp_basic,   only      : ip,rp,lg
  use def_chemic,         only      : typ_chm_mixedEq_group, typ_chm_mixedEq_equa
  implicit none

  !
  ! GROUP TYPES
  !
  integer(ip),   parameter :: &
       CHM_GR_PASSIVE                  = 1_ip, &
       CHM_GR_CONTROL                  = 2_ip, &
       CHM_GR_REACTIVE                 = 3_ip, &
       CHM_GR_SECTIONAL                = 4_ip

  !
  ! EQUATION TYPES
  ! The intention here is that it
  ! starts with the digit of the group,
  ! As of now, I do not expect more than 9
  ! different equation types in a group type
  !
  integer(ip),   parameter :: &
       CHM_EQ_PASSIVE                  = 11_ip, &
       CHM_EQ_Z                        = 21_ip, &
       CHM_EQ_YC                       = 22_ip, &
       CHM_EQ_ZVAR                     = 23_ip, &
       CHM_EQ_YCVAR                    = 24_ip, &
       CHM_EQ_ZZ                       = 25_ip, &
       CHM_EQ_YCYC                     = 26_ip, &
       CHM_EQ_REACTIVE                 = 31_ip, &
       CHM_EQ_SECTIONAL                = 41_ip

  !
  ! VOLUMETRIC SOURCE TERM TYPES
  !
  integer(ip),   parameter :: &
       CHM_SRC_OFF                     = 0_ip, &
       CHM_SRC_TABLE                   = 1_ip, &
       CHM_SRC_DSM                     = 2_ip, &
       CHM_SRC_ANN                     = 3_ip

  !
  ! INITIALIZATION TYPES
  !
  integer(ip),   parameter :: &
       CHM_INI_OLDWAY                  = -666_ip, &
       CHM_INI_OFF                     = 0_ip, &
       CHM_INI_FIELD                   = 1_ip, &
       CHM_INI_CONST                   = 2_ip


  public CHM_GR_PASSIVE
  public CHM_GR_CONTROL
  public CHM_GR_REACTIVE
  public CHM_GR_SECTIONAL

  public CHM_EQ_PASSIVE
  public CHM_EQ_Z
  public CHM_EQ_YC
  public CHM_EQ_ZVAR
  public CHM_EQ_YCVAR
  public CHM_EQ_ZZ
  public CHM_EQ_YCYC
  public CHM_EQ_REACTIVE
  public CHM_EQ_SECTIONAL

  public CHM_SRC_OFF
  public CHM_SRC_TABLE
  public CHM_SRC_DSM
  public CHM_SRC_ANN

  public CHM_INI_OLDWAY
  public CHM_INI_OFF
  public CHM_INI_FIELD
  public CHM_INI_CONST

  public chm_mixedEq_getGroupType
  public chm_mixedEq_getEqType
  public chm_mixedEq_getSrcType
  public chm_mixedEq_getIniType

  public chm_mixedEq_setGroupType
  public chm_mixedEq_EqTyp2kfl
  public chm_mixedEq_setEqType
  public chm_mixedEq_setSrcType
  public chm_mixedEq_setIniType

  public chm_mixedEq_initialize
  public chm_mixedEq_par_exchange
  private

contains
  !
  ! Get type of group
  !
  function chm_mixedEq_getGroupType(group) result(typ)
    type(typ_chm_mixedEq_group), intent(inout)   :: group
    character(10)                                :: typ

    select case(group % kfl_grtype)
    case(CHM_GR_PASSIVE)
       typ = 'PASSIVE'
    case(CHM_GR_CONTROL)
       typ = 'CONTROL'
    case(CHM_GR_REACTIVE)
       typ = 'REACTIVE'
    case(CHM_GR_SECTIONAL)
       typ = 'SECTIONAL'
    case default
       typ = 'NONE'
    end select
  end function chm_mixedEq_getGroupType

  !
  ! Get type of equation
  !
  function chm_mixedEq_getEqType(equa,isshort) result(typ)
    type(typ_chm_mixedEq_equa), intent(inout)   :: equa
    logical(lg), optional,      intent(in)      :: isshort
    character(30)                               :: typ
    logical(lg)                                 :: isshort_loc

    !
    ! Dafault lenght
    !
    isshort_loc = .false.
    if (present(isshort)) isshort_loc = isshort

    select case(equa % kfl_eqtype)
    case(CHM_EQ_PASSIVE)
       typ = 'PASSIVE'
       if (isshort_loc) typ = 'PASSI'
    case(CHM_EQ_Z)
       typ = 'MIXTURE FRACION'
       if (isshort_loc) typ = 'ZMEAN'
    case(CHM_EQ_YC)
       typ = 'PROGRESS VARIABLE'
       if (isshort_loc) typ = 'CMEAN'
    case(CHM_EQ_ZVAR)
       typ = 'MIXTURE FRACION  VARIANCE'
       if (isshort_loc) typ = 'ZVAR '
    case(CHM_EQ_YCVAR)
       typ = 'PROGRESS VARIABLE VARIANCE'
       if (isshort_loc) typ = 'CVAR '
    case(CHM_EQ_ZZ)
       typ = 'MIXTURE FRACION  SQUARE'
       if (isshort_loc) typ = 'ZZ   '
    case(CHM_EQ_YCYC)
       typ = 'PROGRESS VARIABLE SQUARE'
       if (isshort_loc) typ = 'CC   '
    case(CHM_EQ_REACTIVE)
       typ = 'REACTIVE'
       if (isshort_loc) typ = 'REACT'
    case(CHM_EQ_SECTIONAL)
       typ = 'SECTIONAL'
       if (isshort_loc) typ = 'SECTI'
    case default
       typ = 'NONE '
    end select
  end function chm_mixedEq_getEqType

  !
  ! Get type of source
  !
  function chm_mixedEq_getSrcType(equa) result(typ)
    type(typ_chm_mixedEq_equa), intent(inout)   :: equa
    character(10)                               :: typ

    select case(equa % kfl_source_type)
    case(CHM_SRC_OFF)
       typ = 'NONE '
    case(CHM_SRC_TABLE)
       typ = 'TABLE'
    case(CHM_SRC_DSM)
       typ = 'DSM'
    case(CHM_SRC_ANN)
       typ = 'ANN'
    case default
       typ = 'NONE '
    end select
  end function chm_mixedEq_getSrcType

  !
  ! Get type of initial condition
  !
  function chm_mixedEq_getIniType(equa) result(typ)
    type(typ_chm_mixedEq_equa), intent(inout)   :: equa
    character(10)                               :: typ

    select case(equa % kfl_ini_type)
    case(CHM_INI_OLDWAY)
       typ = 'OLDWA'
    case(CHM_INI_OFF)
       typ = 'NONE '
    case(CHM_INI_FIELD)
       typ = 'FIELD'
    case(CHM_INI_CONST)
       typ = 'CONST'
    case default
       typ = 'NONE '
    end select
  end function chm_mixedEq_getIniType



  !
  ! Get equation type flag from type
  !
  function chm_mixedEq_EqTyp2kfl(typ) result(kfl_eqtype)
    character(5), intent(in)        :: typ
    integer(ip)                     :: kfl_eqtype

    select case(typ)
    case('PASSI')
       kfl_eqtype = CHM_EQ_PASSIVE
    case('Z    ','ZMEAN','MIXFR')
       kfl_eqtype = CHM_EQ_Z
    case('YC   ','C    ','YCMEA','CMEAN','PROGV')
       kfl_eqtype = CHM_EQ_YC
    case('ZVAR ')
       kfl_eqtype = CHM_EQ_ZVAR
    case('YCVAR','CVAR ')
       kfl_eqtype = CHM_EQ_YCVAR
    case('ZZ   ')
       kfl_eqtype = CHM_EQ_ZZ
    case('YCYC ','CC   ')
       kfl_eqtype = CHM_EQ_YCYC
    case('REACT')
       kfl_eqtype = CHM_EQ_REACTIVE
    case('SECTI')
       kfl_eqtype = CHM_EQ_SECTIONAL
    case default
       kfl_eqtype = 0_ip
    end select
  end function chm_mixedEq_EqTyp2kfl

  !
  ! Get source type flag from type
  !
  function chm_mixedEq_ScrTyp2kfl(typ) result(kfl_source_type)
    character(5), intent(in)        :: typ
    integer(ip)                     :: kfl_source_type

    select case(typ)
    case('PASSI','OFF  ','NONE ')
       kfl_source_type = CHM_SRC_OFF
    case('TABLE')
       kfl_source_type = CHM_SRC_TABLE
    case('DSM')
       kfl_source_type = CHM_SRC_DSM
    case('ANN')
       kfl_source_type = CHM_SRC_ANN
    case default
       kfl_source_type = CHM_SRC_OFF
    end select
  end function chm_mixedEq_ScrTyp2kfl

  !
  ! Get initial condition type flag from type
  !
  function chm_mixedEq_IniTyp2kfl(typ) result(kfl_ini_type)
    character(5), intent(in)        :: typ
    integer(ip)                     :: kfl_ini_type

    select case(typ)
    case('OLDWA')
       kfl_ini_type = CHM_INI_OLDWAY
    case('OFF  ','NONE ','NULL ')
       kfl_ini_type = CHM_INI_OFF
    case('FIELD')
       kfl_ini_type = CHM_INI_FIELD
    case('CONST')
       kfl_ini_type = CHM_INI_CONST
    case default
       kfl_ini_type = CHM_INI_OFF
    end select
  end function chm_mixedEq_IniTyp2kfl

  !
  ! Set group type from input
  !
  subroutine chm_mixedEq_setGroupType(group,typ)
    type(typ_chm_mixedEq_group), intent(inout)   :: group
    character(5),                intent(in)      :: typ

    select case(typ)
    case('PASSI')
       group % kfl_grtype = CHM_GR_PASSIVE
    case('CONTR')
       group % kfl_grtype = CHM_GR_CONTROL
    case('REACT')
       group % kfl_grtype = CHM_GR_REACTIVE
    case('SECTI')
       group % kfl_grtype = CHM_GR_SECTIONAL
    case default
       group % kfl_grtype = 0_ip
    end select
    if (group % name == '') group % name = trim(typ(1:5))
  end subroutine chm_mixedEq_setGroupType


  !
  ! Set equation type from input
  !
  subroutine chm_mixedEq_setEqType(equa,typ)
    type(typ_chm_mixedEq_equa), intent(inout)   :: equa
    character(5),               intent(in)      :: typ

    equa % kfl_eqtype = chm_mixedEq_EqTyp2kfl(typ)
    if (equa % name == '') equa % name = trim(chm_mixedEq_getEqType(equa,isshort=.true.))
  end subroutine chm_mixedEq_setEqType

  !
  ! Set equation source type from input
  !
  subroutine chm_mixedEq_setSrcType(equa,srctyp)
    type(typ_chm_mixedEq_equa), intent(inout)   :: equa
    character(5),               intent(in)      :: srctyp

    equa % kfl_source_type = chm_mixedEq_ScrTyp2kfl(srctyp)
  end subroutine chm_mixedEq_setSrcType

  !
  ! Set equation initial condition type from input
  !
  subroutine chm_mixedEq_setIniType(equa,inityp)
    type(typ_chm_mixedEq_equa), intent(inout)   :: equa
    character(5),               intent(in)      :: inityp

    equa % kfl_ini_type = chm_mixedEq_IniTyp2kfl(inityp)
  end subroutine chm_mixedEq_setIniType




  !
  ! Initialize a group
  !
  subroutine chm_mixedEq_init_group(group)
    type(typ_chm_mixedEq_group), intent(inout)   :: group
    group % name            = ''
    group % kfl_grtype      = 0_ip
    group % nequa           = 0_ip
    group % i_start         = 0_ip
    group % i_end           = 0_ip
    group % kfl_therm_phor  = 0_ip
  end subroutine chm_mixedEq_init_group




  !
  ! Initialize an equation
  !
  subroutine chm_mixedEq_init_equa(equa)
    type(typ_chm_mixedEq_equa), intent(inout)   :: equa
    equa % name               = ''
    equa % kfl_eqtype         = 0_ip
    equa % kfl_source_type    = 0_ip
    equa % kfl_source_fw      = 0_ip
    equa % kfl_source_ann     = 0_ip
    equa % kfl_source_col     = 0_ip
    equa % kfl_consum_col     = 1_ip
    equa % kfl_ini_type       = 0_ip
    equa % kfl_ini_field      = 0_ip
    equa % ini_value          = 0.0_rp
    equa % kfl_ieq_mean       = 0_ip
    equa % kfl_do_post        = 0_ip
    equa % Lewis              = 1.0_rp    ! Default is 1
    equa % kfl_source_split   = 0_ip
    equa % kfl_fix_diffusion  = 0_ip
    equa % diffusivity        = 0.0_rp
    equa % kfl_diffsource_fw  = 0_ip
    equa % kfl_premsource_fw  = 0_ip
    equa % kfl_diffsource_col = 0_ip
    equa % kfl_premsource_col = 0_ip
  end subroutine chm_mixedEq_init_equa




  !
  ! Initialize groups and equations
  !
  subroutine chm_mixedEq_initialize(ngrou,nequa,groups,equas)
    integer(ip), intent(in)                               :: ngrou
    integer(ip), intent(in)                               :: nequa
    type(typ_chm_mixedEq_group), pointer, intent(inout)   :: groups(:)
    type(typ_chm_mixedEq_equa),  pointer, intent(inout)   :: equas(:)

    integer(ip) ii

    do ii = 1,ngrou
       call chm_mixedEq_init_group(groups(ii))
    enddo

    do ii = 1,nequa
       call chm_mixedEq_init_equa(equas(ii))
    enddo
  end subroutine chm_mixedEq_initialize




  !
  ! Communicate groups and equations
  !
  subroutine chm_mixedEq_par_exchange(ngrou,nequa,groups,equas)
    use mod_exchange,       only : exchange_init
    use mod_exchange,       only : exchange_end
    integer(ip), intent(in)                               :: ngrou
    integer(ip), intent(in)                               :: nequa
    type(typ_chm_mixedEq_group), pointer, intent(inout)   :: groups(:)
    type(typ_chm_mixedEq_equa),  pointer, intent(inout)   :: equas(:)

    integer(ip) ii

    call exchange_init()
    do ii = 1,ngrou
       call chm_mixedEq_par_ex_group(groups(ii))
    enddo

    do ii = 1,nequa
       call chm_mixedEq_par_ex_equa(equas(ii))
    enddo
    call exchange_end()

  end subroutine chm_mixedEq_par_exchange






  !
  ! Communicate a group
  !
  subroutine chm_mixedEq_par_ex_group(group)

    use mod_memory,         only : memory_alloca
    use mod_memory,         only : memory_deallo
    use mod_exchange,       only : exchange_init
    use mod_exchange,       only : exchange_add


    type(typ_chm_mixedEq_group),  intent(inout)   :: group

    call exchange_add(group % name          )
    call exchange_add(group % kfl_grtype    )
    call exchange_add(group % nequa         )
    call exchange_add(group % i_start       )
    call exchange_add(group % i_end         )
    call exchange_add(group % kfl_therm_phor)


  end subroutine chm_mixedEq_par_ex_group





  !
  ! Communicate an equation
  !
  subroutine chm_mixedEq_par_ex_equa(equa)

    use mod_exchange,       only : exchange_init
    use mod_exchange,       only : exchange_add
    use mod_exchange,       only : exchange_end

    type(typ_chm_mixedEq_equa), intent(inout)   :: equa


    call exchange_add(equa % name           )
    call exchange_add(equa % kfl_eqtype     )
    call exchange_add(equa % kfl_source_type)
    call exchange_add(equa % kfl_source_fw  )
    call exchange_add(equa % kfl_source_ann )
    call exchange_add(equa % kfl_source_col )
    call exchange_add(equa % kfl_consum_col )
    call exchange_add(equa % kfl_ini_type   )
    call exchange_add(equa % kfl_ini_field  )
    call exchange_add(equa % ini_value      )
    call exchange_add(equa % kfl_ieq_mean   )
    call exchange_add(equa % kfl_do_post    )
    call exchange_add(equa % Lewis          )
    call exchange_add(equa % kfl_source_split)
    call exchange_add(equa % kfl_fix_diffusion)
    call exchange_add(equa % diffusivity)
    call exchange_add(equa % kfl_diffsource_fw )           
    call exchange_add(equa % kfl_premsource_fw )           
    call exchange_add(equa % kfl_diffsource_col)
    call exchange_add(equa % kfl_premsource_col)

  end subroutine chm_mixedEq_par_ex_equa


end module mod_chm_mixedEq
