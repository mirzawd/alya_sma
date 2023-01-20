!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_outcpu()
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_outcpu
  ! NAME
  !    chm_outcpu
  ! DESCRIPTION
  !    This routine writes a summary of spent computing time. The
  !    quantities computed by the code are:
  !
  !    CPUTI_CHM(1) .... PDE Element assembly
  !    CPUTI_CHM(2) .... PDE Boundary assembly
  !    CPUTI_CHM(3) .... PDE Solver
  !    CPUTI_CHM(4) .... PDE Element assembly
  !    CPUTI_CHM(5) .... PDE Boundary assembly
  !    CPUTI_CHM(6) .... PDE Solver
  !
  ! USES
  !    outfor
  ! USED BY
  !    chm_turnof
  !***
  !-----------------------------------------------------------------------
  use def_master, only : momod, INOTSLAVE, modul, routp, coutp, cpu_modul
  use def_kintyp, only : ip, rp
  use def_chemic, only : cputi_chm
  use mod_outfor, only : outfor
  implicit none
  real(rp) :: xfact,cpuot,cputo

  if( INOTSLAVE ) then
     !
     ! Initializations
     !
     cputo = cpu_modul(30,modul)
     !cputo    = cputi_chm(1)+cputi_chm(2)+cputi_chm(3)+cputi_chm(4)+cputi_chm(5)
     routp(1) = cputo
     call outfor(29_ip,momod(modul)%lun_outpu,' ')

     if( routp(1) > 0.0_rp ) then
        xfact = 100.0_rp / routp(1)
     else
        xfact = 1.0_rp
     end if
     cpuot = 0.0_rp

     coutp(1)  = 'PDE MATRIX AND RHS ELEMENT ASSEMBLIES'
     routp(1)  = cputi_chm(1)
     routp(2)  = xfact * routp(1)
     cpuot     = cpuot + routp(1)
     call outfor(30_ip,momod(modul)%lun_outpu,' ')

     coutp(1)  = 'PDE MATRIX AND RHS BOUNDARY ASSEMBLIES'
     routp(1)  = cputi_chm(2)
     routp(2)  = xfact * routp(1)
     cpuot     = cpuot + routp(1)
     call outfor(30_ip,momod(modul)%lun_outpu,' ')

     coutp(1)  = 'PDE SOLVER'
     routp(1)  = cputi_chm(3)
     routp(2)  = xfact * routp(1)
     cpuot     = cpuot + routp(1)
     call outfor(30_ip,momod(modul)%lun_outpu,' ')

     coutp(1)  = 'ODE ASSEMBLIES'
     routp(1)  = cputi_chm(4)
     routp(2)  = xfact * routp(1)
     cpuot     = cpuot + routp(1)
     call outfor(30_ip,momod(modul)%lun_outpu,' ')

     coutp(1)  = 'ODE SOLVER'
     routp(1)  = cputi_chm(5)
     routp(2)  = xfact * routp(1)
     cpuot     = cpuot + routp(1)
     call outfor(30_ip,momod(modul)%lun_outpu,' ')
     !
     ! Others
     !
     cpuot     = cputo - cpuot
     coutp(1)  = 'OTHERS'
     routp(1)  = cpuot
     routp(2)  = xfact * routp(1)
     call outfor(30_ip,momod(modul)%lun_outpu,' ')

  end if

end subroutine chm_outcpu
