!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Temper
!> @{
!> @file    mod_temper.f90
!> @author  houzeaux
!> @date    2020-04-02
!> @brief   Temper main module
!> @details Manage memory of temper, etc.
!-----------------------------------------------------------------------

module mod_temper

  use def_master
  use def_domain
  use mod_memory, only : memory_alloca
  use mod_memory, only : memory_deallo
  use def_temper

  private

  public :: temper_memory_allocate
  public :: temper_memory_deallocate
  public :: temper_memory_reallocate
  public :: temper_initialization
  public :: temper_main
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-02
  !> @brief   IMain subroutine
  !> @details Main subroutine
  !> 
  !-----------------------------------------------------------------------

  subroutine temper_main(order)
    implicit none
    integer(ip), intent(in) :: order
    call Temper(order)
  end subroutine temper_main

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-02
  !> @brief   Initialization
  !> @details Initialization
  !> 
  !-----------------------------------------------------------------------

  subroutine temper_initialization()

    nullify(kfl_funty_tem)   
    nullify(funpa_tem)     
    nullify(tncod_tem)       
    nullify(tgcod_tem)       
    nullify(tbcod_tem)       
    nullify(kfl_fixno_tem)
    nullify(kfl_fixbo_tem)  
    nullify(kfl_funno_tem)  
    nullify(kfl_funbo_tem)   
    nullify(kfl_funtn_tem)  
    nullify(kfl_funtb_tem)   
    nullify(bvess_tem)  
    nullify(bvnat_tem)   
    nullify(fixnb_tem)        
    nullify(lmatb_tem)      
    nullify(gradc_tem)
    nullify(power_tem)      
    nullify(avtem_tem)      
    nullify(fvtem_tem)      
    nullify(fvte2_tem)      
    nullify(avres_tem)
    nullify(grtem_tem)    
    nullify(teold_tem)     
    nullify(avte2_tem) 
    nullify(avtev_tem)     
    nullify(avden_tem)      
    nullify(fvvel_tem)     
    nullify(dt_rho_cp_tem)
    nullify(massb_tem)     
    nullify(heatf_tem)     
    nullify(lsour_material_tem)
    nullify(xsour_material_tem)
    nullify(avhsk_tem)      

  end subroutine temper_initialization

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-09-22
  !> @brief   Reallocate array
  !> @details Reallocate a domain array
  !> 
  !-----------------------------------------------------------------------

  subroutine temper_memory_reallocate(wvari,NUMBER1,NUMBER2)
    character(len=*), intent(in)           :: wvari
    integer(ip),      intent(in), optional :: NUMBER1
    integer(ip),      intent(in), optional :: NUMBER2
    call temper_memory_deallocate(wvari,NUMBER1,NUMBER2)
    call temper_memory_allocate  (wvari,NUMBER1,NUMBER2)
  end subroutine temper_memory_reallocate

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-09-22
  !> @brief   Allocate array
  !> @details Allocate a domain array
  !> 
  !-----------------------------------------------------------------------

  subroutine temper_memory_allocate(wvari,NUMBER1,NUMBER2)
    character(len=*), intent(in)           :: wvari
    integer(ip),      intent(in), optional :: NUMBER1
    integer(ip),      intent(in), optional :: NUMBER2
    call temper_memory_allocate_deallocate(0_ip,wvari,NUMBER1,NUMBER2)
  end subroutine temper_memory_allocate

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-09-22
  !> @brief   Deallocate array
  !> @details Deallocate a domain array
  !> 
  !-----------------------------------------------------------------------

  subroutine temper_memory_deallocate(wvari,NUMBER1,NUMBER2)
    character(len=*), intent(in)           :: wvari
    integer(ip),      intent(in), optional :: NUMBER1
    integer(ip),      intent(in), optional :: NUMBER2
    call temper_memory_allocate_deallocate(1_ip,wvari,NUMBER1,NUMBER2)
  end subroutine temper_memory_deallocate

  subroutine temper_memory_allocate_deallocate(itask,wvari,NUMBER1,NUMBER2)

    integer(ip),      intent(in)           :: itask
    character(len=*), intent(in)           :: wvari
    integer(ip),      intent(in), optional :: NUMBER1
    integer(ip),      intent(in), optional :: NUMBER2
    integer(ip)                            :: nbou1

    nbou1 = max(1_ip,nboun)

    select case(trim(wvari))

    case ( 'LSOUR_MATERIAL_TEM' )

       if( itask == 0 ) then
          call memory_alloca(mem_modul(1:2,modul),'LFORC_MATERIAL_NSI','tem_memphy',lsour_material_tem,nmate)
       else
          call memory_deallo(mem_modul(1:2,modul),'LFORC_MATERIAL_NSI','tem_memphy',lsour_material_tem)
       end if

    case ( 'XSOUR_MATERIAL_TEM' )

       if( itask == 0 ) then
          call memory_alloca(mem_modul(1:2,modul),'XFORC_MATERIAL_NSI','tem_memphy',xsour_material_tem,msour_material_tem,nmate)
       else
          call memory_deallo(mem_modul(1:2,modul),'XFORC_MATERIAL_NSI','tem_memphy',xsour_material_tem)
       end if

    end select

  end subroutine temper_memory_allocate_deallocate

end module mod_temper
!> @}
