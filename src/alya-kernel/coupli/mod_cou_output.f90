!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Coupling
!> @{
!> @file    mod_cou_output.f90
!> @author  houzeaux
!> @date    2018-06-19
!> @brief   Output of coupling
!> @details Different output of the couplings
!-----------------------------------------------------------------------

module mod_cou_output

  use def_kintyp,         only : ip,rp
  use def_master,         only : routp,coutp,IMASTER,zeror,intost
  use def_master,         only : lun_outpu
  use def_coupli,         only : mcoup,lun_coupl_res,coupling_type
  use def_coupli,         only : typ_color_coupling
  use def_coupli,         only : lun_coupl_res
  use mod_communications, only : PAR_MAX
  use mod_outfor,         only : outfor
  implicit none
  private

  interface cou_output_timings
     module procedure &
          & cou_output_timings_one,&
          & cou_output_timings_all
  end interface cou_output_timings

  public :: cou_output_timings
 
contains
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-06-19
  !> @brief   CPU time of coupling
  !> @details CPU time of coupling
  !> 
  !-----------------------------------------------------------------------
  
  subroutine cou_output_timings_one(coupling)
    type(typ_color_coupling), intent(inout) :: coupling

    call cou_output_timings_single(coupling,lun_outpu)
    
  end subroutine cou_output_timings_one

  subroutine cou_output_timings_all() 

    integer(ip) :: icoup
   
   do icoup = 1,mcoup
      call cou_output_timings_single(coupling_type(icoup),lun_coupl_res)
    end do
    
  end subroutine cou_output_timings_all
  
  subroutine cou_output_timings_single(coupling,lunit)
    
    type(typ_color_coupling), intent(inout) :: coupling
    integer(ip),              intent(in)    :: lunit
    integer(ip)                             :: nsize,icoup
    real(rp)                                :: xfact

     nsize = size(coupling % cputim)
     call PAR_MAX(nsize,coupling % cputim,'IN THE WORLD')

     if( IMASTER ) then

        icoup     = coupling % number
        coutp(1)  = 'COUPLING '//trim(intost(icoup))
        routp(1)  = max(sum(coupling % cputim),zeror)
        call outfor(90_ip,lunit,' ')
        xfact     = 100.0_rp / routp(1)

        coutp(1)  = 'DEFINE WET NODES' 
        routp(1)  = coupling % cputim(1)
        routp(2)  = xfact * routp(1)
        call outfor(30_ip,lunit,' ')

        coutp(1)  = 'HOLCUT' 
        routp(1)  = coupling % cputim(10)
        routp(2)  = xfact * routp(1)
        call outfor(30_ip,lunit,' ')

        coutp(1)  = 'COMPUTE KD-TREES' 
        routp(1)  = coupling % cputim(2)
        routp(2)  = xfact * routp(1)
        call outfor(30_ip,lunit,' ')

        coutp(1)  = 'DISTRIBUTE WET POINTS' 
        routp(1)  = coupling % cputim(6)
        routp(2)  = xfact * routp(1)
        call outfor(30_ip,lunit,' ')

        coutp(1)  = 'SEND WET COORD TO CPUS' 
        routp(1)  = coupling % cputim(3)
        routp(2)  = xfact * routp(1)
        call outfor(30_ip,lunit,' ')

        coutp(1)  = 'CPUS CHECK OWNING' 
        routp(1)  = coupling % cputim(4)
        routp(2)  = xfact * routp(1)
        call outfor(30_ip,lunit,' ')

        coutp(1)  = 'CREATE COUPLING COMM.' 
        routp(1)  = coupling % cputim(5)
        routp(2)  = xfact * routp(1)
        call outfor(30_ip,lunit,' ')

        coutp(1)  = 'GENERATE TRANSM. MATRIX' 
        routp(1)  = coupling % cputim(7)
        routp(2)  = xfact * routp(1)
        call outfor(30_ip,lunit,' ')

        coutp(1)  = 'SEND/RECV COUPLING' 
        routp(1)  = coupling % cputim(8)
        routp(2)  = xfact * routp(1)
        call outfor(30_ip,lunit,' ')

        coutp(1)  = 'MULT. TRANSM. MATRIX' 
        routp(1)  = coupling % cputim(9)
        routp(2)  = xfact * routp(1)
        call outfor(30_ip,lunit,' ')

     end if

   end subroutine cou_output_timings_single

end module mod_cou_output
!> @}
