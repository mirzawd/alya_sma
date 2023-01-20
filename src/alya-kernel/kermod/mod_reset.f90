!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Reset
!> @{
!> @file    mod_reset.f90
!> @author  guillaume
!> @date    2021-03-24
!> @brief   Manage the reset
!> @details Manage the reset
!-----------------------------------------------------------------------

module mod_reset

  use def_kintyp_basic, only : ip,rp
  use def_master,       only : idmod
  use def_master,       only : mmodu
  use def_master,       only : momod
  use def_master,       only : kfl_modul
  use def_master,       only : ittim
  use def_master,       only : namod
  use def_master,       only : ISLAVE
  use def_kermod,       only : kfl_reset
  use def_kermod,       only : reset_factor
  use def_kermod,       only : kfl_reset_steps
  use mod_messages,     only : messages_live
  use mod_strings,      only : integer_to_string
  use mod_ecoute,       only : ecoute
  use mod_ecoute,       only : ecoute_reach_section
  use def_inpout,       only : param
  use def_inpout,       only : getint
  use def_inpout,       only : getrea
  use def_inpout,       only : getcha
  use def_inpout,       only : words
  use mod_exchange,     only : exchange_init
  use mod_exchange,     only : exchange_add
  use mod_exchange,     only : exchange_end
  use def_kintyp_reset
  
  implicit none
  private

  public :: reset_initialization
  public :: reset_read
  public :: reset_parall
  public :: reset_check
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-03-24
  !> @brief   Reset read
  !> @details Reset read
  !> 
  !-----------------------------------------------------------------------
  
  subroutine reset_parall()

    integer(ip) :: imodu

    if( kfl_reset == 0 ) then

       call exchange_init()
       if( ISLAVE ) call reset_allocate()
       call exchange_add(kfl_reset_steps)
       call exchange_add(reset_factor)
       do imodu = 1,mmodu
          if( kfl_modul(imodu) /= 0 ) then
             call exchange_add(momod(imodu) % reset % num_criteria)
             call exchange_add(momod(imodu) % reset % kfl_criterion)
             call exchange_add(momod(imodu) % reset % param) 
          end if
       end do
       call exchange_end()

    end if

  end subroutine reset_parall
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-03-24
  !> @brief   Reset read
  !> @details Reset read
  !> 
  !-----------------------------------------------------------------------
  
  subroutine reset_read()

    integer(ip) :: imodu,icrit


    if( words(2) /= 'OFF  ') then
       kfl_reset = 0
       call reset_allocate()
       
       do while( words(1) /= 'ENDRE' )
         
          if(      words(1) == 'TIMES' ) then
             kfl_reset_steps = getint('TIMES',1_ip  ,'#RESET NUMBER OF TIME STEP TO GO BACK')
          else if ( words(1) == 'FACTO' ) then
             reset_factor    = getrea('FACTO',1.0_rp,'#TIME STEP FACTOR FOR NEW TIME STEP')
          else if(  words(1) == 'METHO' ) then
             imodu = 0
             icrit = 0             
             call ecoute('reset_read')
             do while ( words(1) /= 'ENDME' ) 
                if(      words(1) == 'MODUL' ) then
                   imodu = idmod(getcha('MODUL','NULL ','#moodule name'))
                   momod(imodu) % reset % num_criteria = momod(imodu) % reset % num_criteria + 1
                   icrit = momod(imodu) % reset % num_criteria
                else if( words(1) == 'CRITE' ) then
                   if( imodu == 0 ) call runend('MOD_REASET: MODULE HAS NOT BEEN DEFINED')
                   if(      words(2) == 'SOLVE' ) then
                      momod(imodu) % reset % kfl_criterion(icrit) = RESET_SOLVER_NOT_CONVERGED 
                   else if( words(2) == 'INNER' ) then
                      momod(imodu) % reset % kfl_criterion(icrit) = RESET_INNER_NOT_CONVERGED  
                   else if( words(2) == 'TIME ' ) then
                      momod(imodu) % reset % kfl_criterion(icrit) = RESET_TIME
                   else
                      if( imodu == 0 ) call runend('MOD_REASET: UNKNOWN CRITERION')            
                   end if                   
                else if( words(1) == 'PARAM' ) then
                   if( imodu == 0 ) call runend('MOD_REASET: MODULE    HAS NOT BEEN DEFINED')
                   if( icrit == 0 ) call runend('MOD_REASET: CRITERION HAS NOT BEEN DEFINED')
                   momod(imodu) % reset % param(1:10,icrit) = param(1:10)
                end if
                call ecoute('reset_read')
             end do
          
          end if
          call ecoute('reset_read')
       end do
       
    else

       call ecoute_reach_section('ENDRE')
       
    end if
    
  end subroutine reset_read
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-03-24
  !> @brief   Reset initialization
  !> @details Reset initialization
  !> 
  !-----------------------------------------------------------------------
  
  subroutine reset_initialization()

    integer(ip) :: imodu
    
    do imodu = 1,mmodu
       nullify(momod(imodu) % reset)
    end do
    
  end subroutine reset_initialization
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-03-24
  !> @brief   Reset allocate
  !> @details Reset allocate
  !> 
  !-----------------------------------------------------------------------
  
  subroutine reset_allocate()

    integer(ip) :: imodu
    
    if( kfl_reset /= -1 ) then
       do imodu = 1,mmodu
          allocate(momod(imodu) % reset)
          call momod(imodu) % reset % init()
       end do
    end if
    
  end subroutine reset_allocate
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-03-24
  !> @brief   Reset check
  !> @details Reset check...
  !> 
  !-----------------------------------------------------------------------
  
  subroutine reset_check()

    integer(ip) :: imodu,icrit
    
    if_reset: if( kfl_reset == 0 ) then
       do imodu = 1,mmodu
          if( kfl_modul(imodu) /= 0 ) then
             do icrit = 1,momod(imodu) % reset % num_criteria
                if( momod(imodu) % reset % kfl_criterion(icrit) /= RESET_NULL ) then
                   call reset_criterion(imodu,icrit)
                   if( kfl_reset == 1 ) exit if_reset
                end if
             end do
          end if
       end do
    end if if_reset

  end subroutine reset_check

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-03-24
  !> @brief   Reset check a specific criterion
  !> @details Reset check a specific criterion
  !> 
  !-----------------------------------------------------------------------
  
  subroutine reset_criterion(imodu,icrit)

    integer(ip), intent(in)  :: imodu
    integer(ip), intent(in)  :: icrit
    integer(ip), save        :: ipass = 0
    integer(ip)              :: ivari
    
    if( ipass /= -1 ) ipass = ittim
    
    select case ( momod(imodu) % reset % kfl_criterion(icrit) )

    case ( RESET_SOLVER_NOT_CONVERGED )
       !
       ! Check all the solver convergences
       !
       if( associated(momod(imodu) % solve) ) then
          do ivari = 1,size(momod(imodu) % solve,KIND=ip)
             if( momod(imodu) % solve(ivari) % iters >= momod(imodu) % solve(ivari) % miter ) then
                kfl_reset = 1
                call messages_live('RESET: SOLVER FOR PROBLEM '//trim(momod(imodu) % solve(ivari) % wprob)&
                     //' OF MODULE '//namod(imodu)&
                     //' HAS NOT CONVERGED IN '//integer_to_string(momod(imodu) % solve(ivari) % miter)//' ITERATIONS')
                return
             end if
          end do
       end if
       
    case ( RESET_TIME )
       !
       ! Time has been reached
       !
       if( ipass == int(momod(imodu) % reset % param(1,icrit),ip) ) then
          ipass     = -1
          kfl_reset =  1
       end if
    
    end select
    
  end subroutine reset_criterion
  
end module mod_reset
!> @}

