!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup CPU_Time
!> @{
!> @file    par_outcpu.f90
!> @author  houzeaux
!> @date    2018-12-30
!> @brief   Parallelization CPU time
!> @details This routine outputs Parall CPU time table
!>          cpu_paral( 1) ... Initialization
!>          cpu_paral( 2) ... Reading
!>          cpu_paral( 3) ... par_partit: Element graphs
!>          cpu_paral( 4) ... par_partit: Exchange Master-slaves
!>          par_arrays:
!>          cpu_paral( 5) ... METIS
!>          Communication strategy:
!>          cpu_paral( 6) ... par_domgra: Subdomain graph
!>          cpu_paral( 7) ... par_duagra: Dual graph: communication graph
!>          cpu_paral( 8) ... par_colgra: Color dual graph
!>          cpu_paral( 9) ... par_commun: Communication table
!>          Permutation arrays:
!>          cpu_paral(11) ... Split interior/boundary nodes and renumber 
!>          cpu_paral(12) ... Permutation arrays
!> @} 
!-----------------------------------------------------------------------

subroutine par_outcpu()
  use      def_master
  use      def_parall
  use mod_outfor, only : outfor
  use mod_parall
  implicit none
  real(rp)    :: cpu_parto,cpu_parpe,cpu_parco,cpu_parmp
  real(rp)    :: onvto,onvco,onvpe
  integer(ip) :: ii

  return 

  if(IMASTER) then

     cpu_parto = sum(cpu_paral)
     cpu_parco = cpu_paral(6)+cpu_paral(7)+cpu_paral(8)+cpu_paral(9)
     cpu_parpe = cpu_paral(11)+cpu_paral(12)
     cpu_parmp = 0.0_rp
     do ii = 21,23
        cpu_parmp = cpu_parmp+cpu_paral(ii)
     end do
     !
     ! Title
     !
     routp(1)=cpu_parto
     call outfor(29_ip,lun_outpu_par,' ')
     !
     ! Avoid dividing by zero
     !
     if(cpu_parto == 0.0_rp) then
        onvto=1.0_rp
     else
        onvto=1.0_rp/cpu_parto
     end if
     if(cpu_parco==0.0_rp) then
        onvco=1.0_rp
     else
        onvco=1.0_rp/cpu_parco
     end if
     if(cpu_parpe==0.0_rp) then
        onvpe=1.0_rp
     else
        onvpe=1.0_rp/cpu_parpe
     end if
     !
     ! Starting operations
     !
     coutp(1) = 'INITIALIZE MPI:'
     routp(1) = cpu_paral(1)
     routp(2) = 100.0_rp*routp(1)*onvto
     call outfor(30_ip,lun_outpu_par,' ')

     coutp(1) = 'READ DATA:'
     routp(1) = cpu_paral(2)
     routp(2) = 100.0_rp*routp(1)*onvto
     call outfor(30_ip,lun_outpu_par,' ')

     coutp(1) = 'ELEMENT GRAPH:'
     routp(1) = cpu_paral(3)
     routp(2) = 100.0_rp*routp(1)*onvto
     call outfor(30_ip,lun_outpu_par,' ')

     coutp(1) = 'METIS - PARTITION GRAPH:'
     routp(1) = cpu_paral(5)
     routp(2) = 100.0_rp*routp(1)*onvto
     call outfor(30_ip,lun_outpu_par,' ')

     !-----------------------------------------------

     coutp(1) = 'COMMUNICATION STRATEGY:'
     routp(1) = cpu_parco
     routp(2) = 100.0_rp*routp(1)*onvto
     call outfor(30_ip,lun_outpu_par,' ')

     coutp(1) = 'Subdomain graph:'
     routp(1) = cpu_paral(6)
     routp(2) = 100.0_rp*routp(1)*onvco
     call outfor(31_ip,lun_outpu_par,' ')

     coutp(1) = 'Dual graph (commun.):'
     routp(1) = cpu_paral(7)
     routp(2) = 100.0_rp*routp(1)*onvco
     call outfor(31_ip,lun_outpu_par,' ')

     coutp(1) = 'Color dual graph:'
     routp(1) = cpu_paral(8)
     routp(2) = 100.0_rp*routp(1)*onvco
     call outfor(31_ip,lun_outpu_par,' ')

     coutp(1) = 'Communcation table:'
     routp(1) = cpu_paral(9)
     routp(2) = 100.0_rp*routp(1)*onvco
     call outfor(31_ip,lun_outpu_par,' ')

     !-----------------------------------------------

     coutp(1) = 'PERMUTATION ARRAYS:'
     routp(1) = cpu_parpe
     routp(2) = 100.0_rp*routp(1)*onvto
     call outfor(30_ip,lun_outpu_par,' ')

     coutp(1) = 'Split interior/boundary:'
     routp(1) = cpu_paral(11)
     routp(2) = 100.0_rp*routp(1)*onvpe
     call outfor(31_ip,lun_outpu_par,' ')

     coutp(1) = 'Compute arrays:'
     routp(1) = cpu_paral(12)
     routp(2) = 100.0_rp*routp(1)*onvpe
     call outfor(31_ip,lun_outpu_par,' ')

     !-----------------------------------------------
     !
     ! Communication times
     !
     coutp(1) = 'TOTAL COMMUNICATIONS:'
     routp(1) = cpu_parmp
     routp(2) = 100.0_rp*routp(1)*onvto
     call outfor(30_ip,lun_outpu_par,' ')  

     coutp(1) = 'Send:'
     routp(1) = cpu_paral(21)
     routp(2) = 100.0_rp*routp(1)*onvto
     call outfor(31_ip,lun_outpu_par,' ')

     coutp(1) = 'Receive:'
     routp(1) = cpu_paral(22)
     routp(2) = 100.0_rp*routp(1)*onvto
     call outfor(31_ip,lun_outpu_par,' ')

     coutp(1) = 'Broadcast:'
     routp(1) = cpu_paral(23)
     routp(2) = 100.0_rp*routp(1)*onvto
     call outfor(31_ip,lun_outpu_par,' ')

     coutp(1) = 'Operations:'
     routp(1) = cpu_paral(24)
     routp(2) = 100.0_rp*routp(1)*onvto
     call outfor(31_ip,lun_outpu_par,' ')

     coutp(1) = 'Send_Receive:'
     routp(1) = cpu_paral(25)
     routp(2) = 100.0_rp*routp(1)*onvto
     call outfor(31_ip,lun_outpu_par,' ')

  end if

end subroutine par_outcpu
