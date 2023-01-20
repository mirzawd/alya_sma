!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Alefor
!> @{
!> @file    ale_restar.f90
!> @author  Guillaume Houzeaux
!> @date    16/11/1966
!> @brief   Restarting.
!> @details Restarting.
!>          ITASK = 1 ... Reads the initial values from the restart file
!>                  2 ... Writes restart file
!> @} 
!-----------------------------------------------------------------------

subroutine ale_restar(itask)

  use def_master
  use def_domain
  use def_alefor
  use mod_ale_arrays
  use mod_communications, only : PAR_BROADCAST

  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,idime,iimbo
  integer(ip)             :: kpoin

  !----------------------------------------------------------------------
  !
  ! Primary arrays
  !
  !----------------------------------------------------------------------

  if( itask == ITASK_READ_RESTART ) then

     call ale_arrays('READ RESTART')
     do ipoin = 1,npoin
        coord(:,ipoin) = coord_ale(:,ipoin,1)
     end do
     kfl_domar = 1

  else if( itask == ITASK_WRITE_RESTART ) then

     call ale_arrays('WRITE RESTART')

  end if

  !----------------------------------------------------------------------
  !
  ! Rigid Body 
  !
  !----------------------------------------------------------------------

  if( kfl_rigid_ale == 1 ) then

     if( itask == ITASK_READ_RESTART ) then

        if( INOTSLAVE ) then

           do iimbo = 1,nrbod

              read(momod(modul) % lun_rstar) rbbou(iimbo) % massa
              read(momod(modul) % lun_rstar) rbbou(iimbo) % densi
              read(momod(modul) % lun_rstar) rbbou(iimbo) % volum
              read(momod(modul) % lun_rstar) rbbou(iimbo) % momin
              read(momod(modul) % lun_rstar) rbbou(iimbo) % posgr

              read(momod(modul) % lun_rstar) rbbou(iimbo) % posil  
              read(momod(modul) % lun_rstar) rbbou(iimbo) % velol  
              read(momod(modul) % lun_rstar) rbbou(iimbo) % accel  
              read(momod(modul) % lun_rstar) rbbou(iimbo) % force 
              read(momod(modul) % lun_rstar) rbbou(iimbo) % vpfor 
              read(momod(modul) % lun_rstar) rbbou(iimbo) % pforce
              read(momod(modul) % lun_rstar) rbbou(iimbo) % vforce

              read(momod(modul) % lun_rstar) rbbou(iimbo) % posia 
              read(momod(modul) % lun_rstar) rbbou(iimbo) % veloa 
              read(momod(modul) % lun_rstar) rbbou(iimbo) % accea 
              read(momod(modul) % lun_rstar) rbbou(iimbo) % rotac 
              read(momod(modul) % lun_rstar) rbbou(iimbo) % torqu 
              read(momod(modul) % lun_rstar) rbbou(iimbo) % vptor 
              read(momod(modul) % lun_rstar) rbbou(iimbo) % ptorqu
              read(momod(modul) % lun_rstar) rbbou(iimbo) % vtorqu
              read(momod(modul) % lun_rstar) rbbou(iimbo) % quate 
              read(momod(modul) % lun_rstar) rbbou(iimbo) % q_dot 
           end do
        end if

        do iimbo = 1,nrbod
           call PAR_BROADCAST(          rbbou(iimbo) % massa )
           call PAR_BROADCAST(          rbbou(iimbo) % densi )
           call PAR_BROADCAST(          rbbou(iimbo) % volum )
           call PAR_BROADCAST(6_ip,     rbbou(iimbo) % momin )
           call PAR_BROADCAST(3_ip,     rbbou(iimbo) % posgr )

           call PAR_BROADCAST(3_ip,4_ip,rbbou(iimbo) % posil )
           call PAR_BROADCAST(3_ip,4_ip,rbbou(iimbo) % velol )
           call PAR_BROADCAST(3_ip,4_ip,rbbou(iimbo) % accel )
           call PAR_BROADCAST(3_ip,4_ip,rbbou(iimbo) % force )
           call PAR_BROADCAST(3_ip,4_ip,rbbou(iimbo) % vpfor )
           call PAR_BROADCAST(3_ip,     rbbou(iimbo) % pforce)
           call PAR_BROADCAST(3_ip,     rbbou(iimbo) % vforce)

           call PAR_BROADCAST(3_ip,4_ip,rbbou(iimbo) % posia )
           call PAR_BROADCAST(3_ip,4_ip,rbbou(iimbo) % veloa )
           call PAR_BROADCAST(3_ip,4_ip,rbbou(iimbo) % accea )
           call PAR_BROADCAST(3_ip,4_ip,rbbou(iimbo) % rotac )
           call PAR_BROADCAST(3_ip,4_ip,rbbou(iimbo) % torqu )
           call PAR_BROADCAST(3_ip,4_ip,rbbou(iimbo) % vptor )
           call PAR_BROADCAST(3_ip,     rbbou(iimbo) % ptorqu)
           call PAR_BROADCAST(3_ip,     rbbou(iimbo) % vtorqu)
           call PAR_BROADCAST(4_ip,4_ip,rbbou(iimbo) % quate )
           call PAR_BROADCAST(4_ip,4_ip,rbbou(iimbo) % q_dot )
        end do

     else if( itask == ITASK_WRITE_RESTART ) then

        if( INOTSLAVE ) then           
           do iimbo = 1,nrbod
              
              write(momod(modul) % lun_rstar) rbbou(iimbo) % massa
              write(momod(modul) % lun_rstar) rbbou(iimbo) % densi
              write(momod(modul) % lun_rstar) rbbou(iimbo) % volum
              write(momod(modul) % lun_rstar) rbbou(iimbo) % momin
              write(momod(modul) % lun_rstar) rbbou(iimbo) % posgr

              write(momod(modul) % lun_rstar) rbbou(iimbo) % posil   
              write(momod(modul) % lun_rstar) rbbou(iimbo) % velol  
              write(momod(modul) % lun_rstar) rbbou(iimbo) % accel  
              write(momod(modul) % lun_rstar) rbbou(iimbo) % force 
              write(momod(modul) % lun_rstar) rbbou(iimbo) % vpfor 
              write(momod(modul) % lun_rstar) rbbou(iimbo) % pforce
              write(momod(modul) % lun_rstar) rbbou(iimbo) % vforce

              write(momod(modul) % lun_rstar) rbbou(iimbo) % posia 
              write(momod(modul) % lun_rstar) rbbou(iimbo) % veloa 
              write(momod(modul) % lun_rstar) rbbou(iimbo) % accea 
              write(momod(modul) % lun_rstar) rbbou(iimbo) % rotac 
              write(momod(modul) % lun_rstar) rbbou(iimbo) % torqu 
              write(momod(modul) % lun_rstar) rbbou(iimbo) % vptor 
              write(momod(modul) % lun_rstar) rbbou(iimbo) % ptorqu
              write(momod(modul) % lun_rstar) rbbou(iimbo) % vtorqu
              write(momod(modul) % lun_rstar) rbbou(iimbo) % quate  
              write(momod(modul) % lun_rstar) rbbou(iimbo) % q_dot 
           end do
        end if

     end if

  end if

  if( kfl_rigid_ale == 1 .and. itask == ITASK_READ_RESTART ) then
     do iimbo = 1,nrbod
        do kpoin = 1,rbbou(iimbo) % npoib
           ipoin = rbbou(iimbo) % lninv(kpoin) 
           do idime = 1,ndime
              rbbou(iimbo) % cooin(idime,kpoin) = coord_ori(idime,ipoin)
              rbbou(iimbo) % cooib(idime,kpoin) = coord(idime,ipoin)
           end do
        end do
     end do
  end if
  
end subroutine ale_restar
 
