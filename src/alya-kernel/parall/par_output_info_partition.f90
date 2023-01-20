!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_output_info_partition()
  !------------------------------------------------------------------------
  !****f* Parall/par_outinf
  ! NAME
  !    par_outinf
  ! DESCRIPTION
  !    This routine output info of the parallelization
  ! USED BY
  !    Domain
  !***
  !------------------------------------------------------------------------
  use def_master
  use def_domain
  use def_parall
  use mod_parall
  use mod_communications,            only : PAR_SUM
  use mod_communications,            only : PAR_MIN
  use mod_communications,            only : PAR_MAX
  use mod_communications,            only : PAR_AVERAGE
  use mod_par_parallel_partitioning, only : par_partition_weights
  use mod_outfor,                    only : outfor
  implicit none  
  integer(ip) :: koutp 
  integer(ip) :: npoin_sum,nelem_sum,nboun_sum,nzdom_sum
  integer(ip) :: npoin_min,nelem_min,nboun_min,nzdom_min,nneig_min,bdim_min
  integer(ip) :: npoin_max,nelem_max,nboun_max,nzdom_max,nneig_max,bdim_max
  integer(ip) :: npoin_ave,nelem_ave,nboun_ave,nzdom_ave,nneig_ave,bdim_ave
  integer(ip) :: npoin_min_subd,nelem_min_subd,nboun_min_subd,nzdom_min_subd,nneig_min_subd,bdim_min_subd
  integer(ip) :: npoin_max_subd,nelem_max_subd,nboun_max_subd,nzdom_max_subd,nneig_max_subd,bdim_max_subd
  real(rp)    :: relem_ave,rpoin_ave,rboun_ave,rzdom_ave

  if( IPARALL ) then

     npoin_sum = npoin
     nelem_sum = nelem
     nboun_sum = nboun
     nzdom_sum = nzdom
     call PAR_SUM(npoin_sum)
     call PAR_SUM(nelem_sum)
     call PAR_SUM(nboun_sum)
     call PAR_SUM(nzdom_sum)

     npoin_min = npoin
     nelem_min = nelem
     nboun_min = nboun
     nzdom_min = nzdom
     call PAR_MIN(npoin_min,'IN MY CODE',npoin_min_subd)
     call PAR_MIN(nelem_min,'IN MY CODE',nelem_min_subd)
     call PAR_MIN(nboun_min,'IN MY CODE',nboun_min_subd)
     call PAR_MIN(nzdom_min,'IN MY CODE',nzdom_min_subd)

     npoin_max = npoin
     nelem_max = nelem
     nboun_max = nboun
     nzdom_max = nzdom
     call PAR_MAX(npoin_max,'IN MY CODE',npoin_max_subd)
     call PAR_MAX(nelem_max,'IN MY CODE',nelem_max_subd)
     call PAR_MAX(nboun_max,'IN MY CODE',nboun_max_subd)
     call PAR_MAX(nzdom_max,'IN MY CODE',nzdom_max_subd)

     rpoin_ave = real(npoin,rp)
     relem_ave = real(nelem,rp)
     rboun_ave = real(nboun,rp)
     rzdom_ave = real(nzdom,rp)
     call PAR_AVERAGE(rpoin_ave)
     call PAR_AVERAGE(relem_ave)
     call PAR_AVERAGE(rboun_ave)
     call PAR_AVERAGE(rzdom_ave)
     npoin_ave = int(rpoin_ave,ip)
     nelem_ave = int(relem_ave,ip)
     nboun_ave = int(rboun_ave,ip)
     nzdom_ave = int(rzdom_ave,ip)

     koutp =  0
     koutp = koutp + 1; ioutp(koutp) = npart
     koutp = koutp + 1; ioutp(koutp) = npoin_sum
     !
     ! NPOIN
     ! 
     koutp = koutp+1 ; ioutp(koutp) = npoin_min
     koutp = koutp+1 ; ioutp(koutp) = npoin_min_subd
     koutp = koutp+1 ; ioutp(koutp) = npoin_max
     koutp = koutp+1 ; ioutp(koutp) = npoin_max_subd
     koutp = koutp+1 ; ioutp(koutp) = npoin_ave
     !
     ! NELEM
     ! 
     koutp = koutp+1 ; ioutp(koutp) = nelem_min
     koutp = koutp+1 ; ioutp(koutp) = nelem_min_subd
     koutp = koutp+1 ; ioutp(koutp) = nelem_max
     koutp = koutp+1 ; ioutp(koutp) = nelem_max_subd
     koutp = koutp+1 ; ioutp(koutp) = nelem_ave
     routp(1) = real(nelem_ave,rp) / real(max(1_ip,nelem_max),rp)
     !
     ! Weights
     !
     call par_partition_weights(TOTAL_WEIGHT=relem_ave)

     nelem_min = int(relem_ave,ip)
     nelem_max = nelem_min
     call PAR_MIN(nelem_min,'IN MY CODE',nelem_min_subd)
     call PAR_MAX(nelem_max,'IN MY CODE',nelem_max_subd)
     
     call PAR_AVERAGE(relem_ave)
     nelem_ave = int(relem_ave,ip)
    
     koutp = koutp+1 ; ioutp(koutp) = nelem_min
     koutp = koutp+1 ; ioutp(koutp) = nelem_min_subd
     koutp = koutp+1 ; ioutp(koutp) = nelem_max
     koutp = koutp+1 ; ioutp(koutp) = nelem_max_subd
     koutp = koutp+1 ; ioutp(koutp) = nelem_ave
     routp(2) = real(nelem_ave,rp) / real(max(1_ip,nelem_max),rp)
     !
     ! NBOUN
     !
     koutp = koutp+1 ; ioutp(koutp) = nboun_min
     koutp = koutp+1 ; ioutp(koutp) = nboun_min_subd
     koutp = koutp+1 ; ioutp(koutp) = nboun_max
     koutp = koutp+1 ; ioutp(koutp) = nboun_max_subd
     koutp = koutp+1 ; ioutp(koutp) = nboun_ave
     !
     ! NZDOM
     !
     koutp = koutp+1 ; ioutp(koutp) = nzdom_min
     koutp = koutp+1 ; ioutp(koutp) = nzdom_min_subd
     koutp = koutp+1 ; ioutp(koutp) = nzdom_max
     koutp = koutp+1 ; ioutp(koutp) = nzdom_max_subd
     koutp = koutp+1 ; ioutp(koutp) = nzdom_ave
     !
     ! NNEIG: number of neighbors
     !
     nneig_min = PAR_COMM_MY_CODE_ARRAY(1) % nneig
     nneig_max = PAR_COMM_MY_CODE_ARRAY(1) % nneig
     nneig_ave = PAR_COMM_MY_CODE_ARRAY(1) % nneig

     call PAR_MIN(nneig_min,'IN MY CODE',nneig_min_subd)
     call PAR_MAX(nneig_max,'IN MY CODE',nneig_max_subd)
     call PAR_AVERAGE(nneig_ave)

     koutp = koutp+1 ; ioutp(koutp) = nneig_min
     koutp = koutp+1 ; ioutp(koutp) = nneig_min_subd
     koutp = koutp+1 ; ioutp(koutp) = nneig_max
     koutp = koutp+1 ; ioutp(koutp) = nneig_max_subd
     koutp = koutp+1 ; ioutp(koutp) = nneig_ave
     !
     ! NBBOU: Number of boundary nodes 
     ! 
     bdim_min = npoin - npoi1
     bdim_max = npoin - npoi1
     bdim_ave = npoin - npoi1

     call PAR_MIN(bdim_min,'IN MY CODE',bdim_min_subd)
     call PAR_MAX(bdim_max,'IN MY CODE',bdim_max_subd)
     call PAR_AVERAGE(bdim_ave)

     koutp = koutp+1 ; ioutp(koutp) = bdim_min
     koutp = koutp+1 ; ioutp(koutp) = bdim_min_subd
     koutp = koutp+1 ; ioutp(koutp) = bdim_max
     koutp = koutp+1 ; ioutp(koutp) = bdim_max_subd
     koutp = koutp+1 ; ioutp(koutp) = bdim_ave
     !
     ! Partitioning method
     !
     call outfor(38_ip,lun_outpu_par)
     call outfor(38_ip,lun_outpu)

  end if

end subroutine par_output_info_partition
