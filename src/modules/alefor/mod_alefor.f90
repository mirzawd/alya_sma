!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Alefor
!> @{
!> @file    mod_alefor.f90
!> @author  houzeaux
!> @date    2020-04-10
!> @brief   Alefor main subroutine
!> @details Mainly for memory management
!-----------------------------------------------------------------------

module mod_alefor
  
  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use def_inpout
  use def_alefor
  use mod_memory

  implicit none
  
  public :: alefor_memory_allocate
  public :: alefor_memory_deallocate
  public :: alefor_initialization
  public :: alefor_main
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-02
  !> @brief   IMain subroutine
  !> @details Main subroutine
  !> 
  !-----------------------------------------------------------------------

  subroutine alefor_main(order)
    implicit none
    integer(ip), intent(in) :: order
    call Alefor(order)
  end subroutine alefor_main

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-02
  !> @brief   Initialization
  !> @details Initialization
  !> 
  !-----------------------------------------------------------------------

  subroutine alefor_initialization()
    
     nullify(lnodb_ad)      
     nullify(ltypb_ad)       
     nullify(coord_ad)       
     nullify(tncod_ale)      
     nullify(tgcod_ale)      
     nullify(tbcod_ale)      
     nullify(kfl_funno_ale)  
     nullify(kfl_funbo_ale)  
     nullify(kfl_funtn_ale)  
     nullify(kfl_funtb_ale)  
     nullify(kfl_fixbo_ale)  
     nullify(kfl_fixrs_ale)  
     nullify(coord_ale)      
     nullify(coord_ori)      
     nullify(rbbou)
     
  end subroutine alefor_initialization

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-09-22
  !> @brief   Allocate array
  !> @details Allocate a domain array
  !> 
  !-----------------------------------------------------------------------

  subroutine alefor_memory_allocate(wvari,NUMBER1,NUMBER2)
    character(len=*), intent(in)           :: wvari
    integer(ip),      intent(in), optional :: NUMBER1
    integer(ip),      intent(in), optional :: NUMBER2
    call alefor_memory_allocate_deallocate(0_ip,wvari,NUMBER1,NUMBER2)
  end subroutine alefor_memory_allocate

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-09-22
  !> @brief   Deallocate array
  !> @details Deallocate a domain array
  !> 
  !-----------------------------------------------------------------------

  subroutine alefor_memory_deallocate(wvari,NUMBER1,NUMBER2)
    character(len=*), intent(in)           :: wvari
    integer(ip),      intent(in), optional :: NUMBER1
    integer(ip),      intent(in), optional :: NUMBER2
    call alefor_memory_allocate_deallocate(1_ip,wvari,NUMBER1,NUMBER2)
  end subroutine alefor_memory_deallocate

  subroutine alefor_memory_allocate_deallocate(itask,wvari,NUMBER1,NUMBER2)

    integer(ip),      intent(in)           :: itask
    character(len=*), intent(in)           :: wvari
    integer(ip),      intent(in), optional :: NUMBER1
    integer(ip),      intent(in), optional :: NUMBER2
    integer(ip)                            :: iimbo,ppoib,pboib

    select case(trim(wvari))

    case ( 'BOUNDARY CONDITIONS')
       !
       ! Boundary conditions
       !
       call memory_alloca(mem_modul(1:2,modul),'KFL_FIXNO_ALE','mod_alefor',kfl_fixno_ale,ndime,npoin)
       call memory_alloca(mem_modul(1:2,modul),'KFL_FIXBO_ALE','mod_alefor',kfl_fixbo_ale,nboun)
       call memory_alloca(mem_modul(1:2,modul),'KFL_FIXRS_ALE','mod_alefor',kfl_fixrs_ale,npoin)
       call memory_alloca(mem_modul(1:2,modul),'BVESS_ALE'    ,'mod_alefor',bvess_ale,ndime,npoin,2_ip)
       
       call memory_alloca(mem_modul(1:2,modul),'KFL_FUNNO_ALE','mod_alefor',kfl_funno_ale,npoin)
       call memory_alloca(mem_modul(1:2,modul),'KFL_FUNTN_ALE','mod_alefor',kfl_funtn_ale,npoin)
       call memory_alloca(mem_modul(1:2,modul),'KFL_FUNBO_ALE','mod_alefor',kfl_funbo_ale,nboun)
       call memory_alloca(mem_modul(1:2,modul),'KFL_FUNTB_ALE','mod_alefor',kfl_funtb_ale,nboun)
       
    case ( 'RBBOU' )
       !
       ! Rigid body type
       !
       if( nrbod > 0 .and. kfl_rigid_ale ==1 ) then
          if( itask == 0 ) then
             call memory_alloca(mem_modul(1:2,modul),'RBBOU','mod_alefor',rbbou,nrbod)
          else
             call memory_deallo(mem_modul(1:2,modul),'RBBOU','mod_alefor',rbbou)
          end if          
       end if
       
    case ( 'RBBOU TOPOLOGY' )
       !
       ! Immersed boundaries (IB) of body fitted type
       !
       if( kfl_rigid_ale == 1 ) then
          if( present(NUMBER1) ) then
             iimbo = NUMBER1
          else
             call runend('MOD_ALEFOR: NUMBER1 IS NECESSARY')
          end if
          if( itask == 0 ) then
             ppoib = rbbou(iimbo) % npoib
             pboib = rbbou(iimbo) % nboib
             call memory_alloca(mem_modul(1:2,modul),'RBBOU','mod_alefor',rbbou(iimbo),ndime,ppoib,pboib,mnoib)
          end if
       end if

    case default

       call runend('MOD_ALEFOR: UNKNOW VARIABLE')
       
    end select

  end subroutine alefor_memory_allocate_deallocate

end module mod_alefor
!> @}
