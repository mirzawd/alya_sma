!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_sld_atm.f90
!> @author  Adria Quintanas-Corominas
!> @author  Gerard Guillamet
!> @author  Daniel Mira
!> @date    2019
!> @brief   Thermomechanical framework
!> @details Thermomechanical framework    
!> @}
!------------------------------------------------------------------------

module mod_sld_atm

  use def_kintyp_basic, only                  :  ip, rp, lg

  implicit none

  integer(ip), parameter                     :: ntmp_atm = 4_ip
  integer(ip), public                        :: kfl_intem_sld           !< Initial condition: temperature
  logical(lg), public                        :: kfl_coupl_atm = .False. !< Flag for multi-module coupling
  logical(lg), public                        :: kfl_therm_sld = .False. !< Flag for temperature
  real(rp), dimension(:), pointer, protected :: tmp0_atm => null()
  real(rp)                                   :: intem_sld               !< Initial value for temperature

  public                                     :: sld_atm_set_flags
  public                                     :: sld_atm_allocate_memmory
  public                                     :: sld_atm_parall
  public                                     :: sld_atm_get_temperature_at_GP
  public                                     :: sld_atm_set_initial_temperature_at_POINTS
  public                                     :: sld_atm_iniunk
  public                                     :: sld_atm_reabcs

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  aquintan
  !> @date    2022-10-07
  !> @brief   Initializations thermomechanical framework
  !> @details Initializations thermomechanical framework
  !> 
  !-----------------------------------------------------------------------
  
  subroutine sld_atm_set_flags()
    use def_master,  only : coupling
    implicit none    
    ! Multi-module coupling
    if( coupling('SOLIDZ','TEMPER') >= 1_ip .or. coupling('TEMPER','SOLIDZ') >= 1_ip )then
       kfl_coupl_atm = .True.
    else
       kfl_coupl_atm = .False.
    endif
    ! One-module coupling
    if( kfl_therm_sld )then
       kfl_coupl_atm = .True.
    endif
  end subroutine sld_atm_set_flags
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  aquintan
  !> @date    2022-10-07
  !> @brief   Memmory for thermomechanics arrays
  !> @details Memmory for thermomechanics arrays
  !> 
  !-----------------------------------------------------------------------
  
  subroutine sld_atm_allocate_memmory()
    use def_domain, only : npoin
    use def_master, only : mem_modul, modul, INOTMASTER
    use mod_memory, only : memory_alloca
    implicit none
    call memory_alloca(mem_modul(1:2,modul), 'TMP0_ATM', 'mod_sld_atm', tmp0_atm, npoin)
    if( INOTMASTER ) tmp0_atm(:) = 0.0_rp
  end subroutine sld_atm_allocate_memmory

  !-----------------------------------------------------------------------
  !> 
  !> @author  aquintan
  !> @date    2022-10-07
  !> @brief   Send data related to thermomechanics
  !> @details Send data related to thermomechanics
  !> 
  !-----------------------------------------------------------------------
  
  subroutine sld_atm_parall()
    use mod_exchange, only : exchange_init
    use mod_exchange, only : exchange_add
    use mod_exchange, only : exchange_end
    implicit none
    call exchange_init()
    call exchange_add(kfl_coupl_atm)
    call exchange_add(kfl_therm_sld)
    call exchange_add(kfl_intem_sld)
    call exchange_add(intem_sld)
    call exchange_end()
  end subroutine sld_atm_parall

  !-----------------------------------------------------------------------
  !> 
  !> @author  aquintan
  !> @date    2022-10-07
  !> @brief   Initial temperature
  !> @details Initial temperature
  !> 
  !-----------------------------------------------------------------------
  
  subroutine sld_atm_set_initial_temperature_at_POINTS()
    use def_master, only : IMASTER, therm
    use def_domain, only : npoin
    implicit none
    if( IMASTER ) return
    tmp0_atm(1:npoin) = therm(1:npoin,3)
  end subroutine sld_atm_set_initial_temperature_at_POINTS
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  aquintan
  !> @date    2022-10-07
  !> @brief   Temperature variable at gauss point level
  !> @details Temperature variable at gauss point level
  !> 
  !-----------------------------------------------------------------------
  
  subroutine sld_atm_get_temperature_at_GP(pnode,pgaus,gpsha,lnods,gp_temp)
    use def_master, only : therm
    implicit none
    integer(ip), intent(in)  :: pnode, pgaus, lnods(:)
    real(rp),    intent(in)  :: gpsha(:,:)
    real(rp),    intent(out) :: gp_temp(:,:)
    integer(ip)              :: inode, igaus, ipoin
    gp_temp = 0.0_rp
    do inode = 1, pnode
       ipoin = lnods(inode)
       do igaus = 1, pgaus
          gp_temp(4,igaus) = gp_temp(4,igaus) + gpsha(inode,igaus)*tmp0_atm(ipoin) ! Initial-Simualtion-value
          gp_temp(3,igaus) = gp_temp(3,igaus) + gpsha(inode,igaus)*therm(ipoin,3)  ! Past-Time-value
          gp_temp(2,igaus) = gp_temp(2,igaus) + gpsha(inode,igaus)*therm(ipoin,2)  ! Past-Iter-value
          gp_temp(1,igaus) = gp_temp(1,igaus) + gpsha(inode,igaus)*therm(ipoin,1)  ! Current-Iter-value
       enddo
    enddo
  end subroutine sld_atm_get_temperature_at_GP
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillam
  !> @date    2022-10-07
  !> @brief   Initial conditions for temperature
  !> @details Initial conditions for temperature
  !> 
  !-----------------------------------------------------------------------
  
  subroutine sld_atm_iniunk()
    use def_master, only : therm
    use def_master, only : ITER_K
    use def_domain, only : xfiel
    use def_domain, only : npoin
    implicit none
    integer(ip)          :: ipoin
    if( kfl_intem_sld > 0 ) then
       do ipoin = 1,npoin
          therm(ipoin,ITER_K) = xfiel(kfl_intem_sld) % a(1,ipoin,1)
       end do
    else
       do ipoin = 1,npoin
          therm(ipoin,ITER_K) = intem_sld
       end do
    end if
  end subroutine sld_atm_iniunk
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillam
  !> @date    2022-10-07
  !> @brief   Read initial conditions for temperature
  !> @details Read initial conditions for temperature
  !> 
  !-----------------------------------------------------------------------
  
  subroutine sld_atm_reabcs()
    use def_inpout, only : words
    use def_inpout, only : getint
    use def_inpout, only : getrea
    implicit none
    !
    ! Initializations
    !
    kfl_intem_sld  = 0_ip                 
    intem_sld = 0.0_rp
    !
    ! Read data
    !
    if(      words(3) == 'FIELD' ) then
       kfl_intem_sld = getint('FIELD',1_ip,'#Field Number for initial temperature variable')
    else if( words(3) == 'VALUE' ) then
       kfl_intem_sld = -1_ip
       intem_sld = getrea('VALUE',0.0_rp,'#Value for initial temperature')
    end if
  end subroutine sld_atm_reabcs
  
end module mod_sld_atm

