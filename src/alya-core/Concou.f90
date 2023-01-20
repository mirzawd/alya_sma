!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine Concou()
  
  !-----------------------------------------------------------------------
  !****f* master/Concou
  ! NAME
  !    Concou
  ! DESCRIPTION
  !    This routine checks the general convergence of the run.
  ! USES
  !    Nastin
  !    Temper
  !    Codire
  !    Turbul
  !    Exmedi
  !    Nastal
  !    Alefor
  !    Latbol
  ! USED BY
  !    Alya
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_coupli
  use mod_parall
  use mod_outfor,   only : outfor
  use mod_messages, only : livinf
  use mod_iofile,   only : iofile_flush_unit
  use mod_moduls,   only : moduls
  implicit none

  integer(ip)             :: imodu
  real(rp)                :: cpu_refer
  integer(ip)             :: iorde
  integer(ip), save       :: itert = 0

  !
  ! Live information
  ! 
  call livinf(7_ip,' ',zero)
  ! 
  ! Initializations
  !
  do iorde = 1,mmodu-1
     if( lmord(iorde,iblok) > 0 ) glres(lmord(iorde,iblok)) = 0.0_rp
  end do
  kfl_gocou = 0
! *******************
! Added for FSI tests
!  kfl_gocou = 1
! *******************
  call cputim(cpu_refer)
  cpu_refer = cpu_refer - cpu_initi
  routp(1)  = cpu_refer 
  call outfor(10_ip,lun_outpu,' ')
  !
  ! Write convergence in the order of module iterations
  !
  call moduls(ITASK_CONCOU)
  !
  ! Coupling convergence
  !
  !call cou_cvgunk()  
  !
  ! Write to convergence file and keep iterating...
  !
  itert = itert + 1
  if( INOTSLAVE ) then
     if( itert == 1 .and. kfl_rstar /= 2 ) then
        write(lun_conve,11)
        write(lun_conve,12,advance='no')
        do imodu = 1,mmodu-1
           write(lun_conve,13,advance='no') namod(imodu),imodu+3
        end do
        write(lun_conve,*)
     end if
     write(lun_conve,14,advance='no') ittim,itcou,cutim
     do imodu = 1,mmodu-1
        if( kfl_modul(imodu) /= 0 ) then
           write(lun_conve,15,advance='no') glres(imodu)
        else
           write(lun_conve,16,advance='no') zero        
        end if
     end do
     write(lun_conve,*)
     call iofile_flush_unit(lun_conve) 
  end if
  !
  ! False convergence
  !
  itcou = itcou + 1
  if( itcou > micou(iblok) ) kfl_gocou = 0  
  !
  ! Formats
  !
11 format('# ','       Time','     Global','       Current','  Module residuals -->')
12 format('# ','       step','  iteration','          time')
13 format('   ',a6,' (',i2,')')
14 format(4x,i9,2x,i9,2x,e12.6)
15 format(2x,e12.6)
16 format(2x,i12)

end subroutine Concou
      
