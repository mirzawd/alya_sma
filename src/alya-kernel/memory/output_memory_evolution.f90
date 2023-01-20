!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Memory
!> @{
!> @file    output_memory_evolution.f90
!> @author  houzeaux
!> @date    2018-09-04
!> @brief   Memory
!> @details Output memory evolution. For the slave count, the maximum
!>          is taken. Output current memory and maximum memory.
!>          Also check for a possible memory leak. Result is in Mb.
!> @} 
!-----------------------------------------------------------------------

subroutine output_memory_evolution()
  
  use def_master
  use mod_memory,         only : mem_curre
  use mod_memory,         only : mem_maxim
  use mod_memory,         only : Mbytes
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_AVERAGE
  use mod_messages,       only : messages_live
  use mod_iofile,         only : iofile_flush_unit
  
  implicit none

  integer(ip)            :: rank_curre
  integer(ip)            :: rank_maxim
  real(rp)               :: mem_curre_rp
  real(rp)               :: mem_maxim_rp
  real(rp)               :: mem_avera_rp
  real(rp)               :: mem_curre_master_rp
  real(rp)               :: mem_maxim_master_rp
  integer(8),  save      :: mem_history(10)
  integer(8)             :: mem_history_cpy(10)
  integer(ip)            :: kk,ii
  integer(ip), save      :: ipass=0
  
  mem_curre_rp = real(mem_curre,rp)
  mem_avera_rp = real(mem_curre,rp)
  mem_maxim_rp = real(mem_maxim,rp)

  call PAR_MAX(mem_curre_rp,rank_max_owner=rank_curre)
  call PAR_MAX(mem_maxim_rp,rank_max_owner=rank_maxim)
  call PAR_AVERAGE(mem_avera_rp)
  !
  ! Output info
  !
  if( INOTSLAVE ) then
     if( ittim==1 ) write(int(lun_memory,4),1)
     
     mem_curre_master_rp = real(mem_curre,rp)
     mem_maxim_master_rp = real(mem_maxim,rp)
     write(int(lun_memory,4),2)&
          ittim,&
          mem_maxim_master_rp/Mbytes,&
          mem_curre_master_rp/Mbytes,&
          mem_avera_rp/Mbytes,&
          mem_maxim_rp/Mbytes,&
          rank_maxim,&
          mem_curre_rp/Mbytes,&
          rank_curre
     call iofile_flush_unit(lun_memory)
  end if
  !
  ! Check possible memory leak
  !
  ipass = ipass + 1
  if( ipass > 10 ) then
     mem_history_cpy = mem_history
     mem_history(1:9) = mem_history_cpy(2:10)
  end if

  mem_history(min(ipass,10_ip)) = mem_curre
  kk = 0
  do ii = 2,10
     if( mem_history(ii) > mem_history(ii-1) ) kk = kk + 1
  end do
  call PAR_MAX(kk)
  if( kk == 9 ) then
     call messages_live('WE HAVE DETECTED A POSSIBLE MEMORY LEAK: CHECK MEMORY EVOLUTION FILE','WARNING')
  end if
  
1 format(&
       &    '#',/,&
       &    '#     Memory units are Mbytes',/,&
       &    '#',/,&
       &    '# ',&
       &    '    Time step',1x,&
       &    '  Max. Master',1x,&
       &    '  Cur. Master',1x,&
       &    '  Ave. Slaves',1x,&
       &    '  Max. Slaves',1x,&
       &    '    Max. Rank',1x,&
       &    '  Cur. Slaves',1x,&
       &    '    Cur. Rank',1x&
       &    )
2 format(2x,(1x,i13),4(1x,e13.6),(1x,i13),(1x,e13.6),(1x,i13))
  
end subroutine output_memory_evolution
