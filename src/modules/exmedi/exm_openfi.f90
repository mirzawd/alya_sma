!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_openfi.f90
!> @author  Mariano Vazquez
!> @date    22/11/2016
!> @brief   Open files
!> @details Open files
!> @} 
!-----------------------------------------------------------------------
subroutine exm_openfi(itask)
!-----------------------------------------------------------------------
!****f* Exmedi/exm_openfi
! NAME 
!    exm_openfi
! DESCRIPTION
!    This subroutine gets ALL the file names and open them to be used by 
!    the module in two possible ways:
! 
!    1. Recalling them from the environment, when Zephyr is launched
!       encapsulated in a shell script, or
! 
!    2. Composing the names out of the problem name which is given as
!       argument when the binary file Zephyr is launched "naked".
! USES
!    iofile
! USED BY
!    exm_turnon
!***
!-----------------------------------------------------------------------
  use      def_exmedi
  use      def_master
  use      def_parame
  use      mod_iofile

  implicit none
  integer(ip), intent(in)  :: itask !> itask
  character(150) :: fil_maxmi_exm
  character(150) :: fil_vinte_exm

  if( INOTSLAVE ) then
     
     if (itask == 1) then
        
        !
        ! kfl_naked is set in the kernel subroutine getnam
        !
        if (kfl_naked==0) then
           !  encapsulated, then get names from the environment  
           call GET_ENVIRONMENT_VARIABLE('FOR510',fil_maxmi_exm)
           call GET_ENVIRONMENT_VARIABLE('FOR511',fil_vinte_exm)
        else if (kfl_naked==1) then
           !  naked, then compose the names     
           fil_maxmi_exm = adjustl(trim(namda))//'.'//exmod(modul)//'.mxm'
           fil_vinte_exm = adjustl(trim(namda))//'.'//exmod(modul)//'.vin'
        end if

        !
        ! Open files
        !
        if(kfl_rstar == 2) then
            call iofile(zero,lun_vinte_exm,fil_vinte_exm,'EXMEDI VOLUME INTEGRALS','old','formatted','append')
        else
            call iofile(zero,lun_vinte_exm,fil_vinte_exm,'EXMEDI VOLUME INTEGRALS')
        end if                

     else if (itask == 2) then
        

     else if (itask== 3) then  
        
     else if (itask == 12) then
        
        
     end if
     
  end if

end subroutine exm_openfi

