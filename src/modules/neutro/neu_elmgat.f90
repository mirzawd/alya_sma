!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Neutro
!> @{
!> @file    neu_elmgat.f90
!> @date    31/03/2016
!> @author  Guillaume Houzeaux
!> @brief   Element gather
!> @details Element gather
!> @} 
!-----------------------------------------------------------------------

subroutine neu_elmgat(pnode,lnods,elcod,elrad,elunk)

  use def_kintyp, only     :  ip,rp 
  use def_domain, only     :  ndime,coord!,npoin
  use def_master, only     :  neutr
  use def_neutro, only     :  ADR_neu,num_energies_neu,num_directions_neu
  use def_neutro, only     :  current_energy_neu,current_direction_neu
  use mod_ADR,    only     :  BDF
  implicit none

  integer(ip), intent(in)  :: pnode 
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(out) :: elcod(ndime,pnode)
  real(rp),    intent(out) :: elrad(num_energies_neu,num_directions_neu,pnode) 
  real(rp),    intent(out) :: elunk(pnode,ADR_neu % ntime)  
  integer(ip)              :: inode,ipoin,itime!,idime,dummi
  !
  ! Current neutroature and coordinates
  !
  do inode = 1,pnode
     ipoin = lnods(inode)
     elrad(1:num_energies_neu,1:num_directions_neu,inode) &
          = neutr(1:num_energies_neu,1:num_directions_neu,ipoin,1)
     elcod(1:ndime,inode) = coord(1:ndime,ipoin)
     elunk(inode,1)       = neutr(current_energy_neu,current_direction_neu,ipoin,1)
  end do
  !
  ! Time integration
  !
  if( ADR_neu % kfl_time_integration /= 0 ) then
     do inode = 1,pnode
        ipoin = lnods(inode)
        elunk(inode,2) = neutr(current_energy_neu,current_direction_neu,ipoin,3)       
     end do
     if( ADR_neu % kfl_time_scheme == BDF ) then
        do itime = 3,ADR_neu % kfl_time_order + 1
           do inode = 1,pnode 
              ipoin = lnods(inode)
              elunk(inode,itime) = neutr(current_energy_neu,current_direction_neu,ipoin,itime+1) 
           end do
        end do
     end if
  end if 

end subroutine neu_elmgat
