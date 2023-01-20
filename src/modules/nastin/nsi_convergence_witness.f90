!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_convergence_witness()
  !------------------------------------------------------------------------
  !****f* Nastin/nsi_convergence_witness
  ! NAME 
  !    nsi_convergence_witness
  ! DESCRIPTION
  !    Convergence based on witness points - Done by the master since he has all the values
  !    Beware it will not work with restart -- THINK
  !
  ! For Iberdrola in .dat file put
  !  TURBUL_MODULE          On
  !    CONVERGENCE:         OFF
  !  END_TURBUL_MODULE
  !
  ! set a logical number of time steps ej 3000
  !
  ! In nsi.dat put STOP_BY_WITNESS
  !
  ! RESTART is not ready - one easy option would be to read the values directly from the nsi.wit file
  !
  !
  ! USES
  ! USED BY
  !    nsi_outwit
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_nastin
  use mod_outfor, only : outfor
  implicit none
  integer(ip) :: iwitn
  
  integer(ip)          :: iaux,jaux 
  integer(ip)                    :: interval
  integer(ip),parameter          :: interval_size = 100   ! OJO VOLVER A 100
  real(rp)                       :: aver(3,2)
  real(rp),parameter             :: toler = 0.01_rp   ! OJO volver a 0.01
  integer(ip)                    :: kount_converged
  integer(ip)                    :: ninterval

!  real(rp),allocatable           :: avwit(:,:)
  integer(ip),parameter          :: max_interval = 500 
  integer(ip),parameter          :: max_wit = 2000 
  real(rp),save                  :: avwit(max_wit,max_interval)=0.0_rp  ! Temporary solution  !otherwise I would have to initialize and alloacte outside time loop

  
  ninterval = 1 + (mitim - 1) / interval_size    ! Beware in .dat file use a logical max time steps

  interval = 1 + (ittim - 1) / interval_size 
  

  if(ninterval > max_interval) call runend('nsi_convergence_witness: ninterval > max_interval')
  if(nwitn > max_wit) call runend('nsi_convergence_witness: nwitn > max_iwit')
!  allocate(avwit(nwitn,ninterval))
!  avwit = 0.0_rp
  !
  ! Convergence based on witness points 
  !
  if( IMASTER .and. ittim >=1 ) then

     iaux = postp(1) % npp_witne(1) + postp(1) % npp_witne(2)
     if ( ndime == 3_ip) iaux = iaux + postp(1) % npp_witne(3)

     if( iaux /= ndime ) call runend('nsi_convergence_witness: witness velocity needed') 
     
     !witne(1,iwitn),witne(2,iwitn),witne(3,iwitn)  - vx , vy , vz

     if( ndime == 3 ) then
        do iwitn = 1,nwitn
           avwit(iwitn,interval) = avwit(iwitn,interval) + sqrt( witne(1,iwitn)**2 + witne(2,iwitn)**2 + witne(3,iwitn)**2 ) 
        end do
     else
        do iwitn = 1,nwitn
           avwit(iwitn,interval) = avwit(iwitn,interval) + sqrt( witne(1,iwitn)**2 + witne(2,iwitn)**2 ) 
        end do
     end if

     if( mod(ittim, interval_size ) == 0 ) then
        kount_converged = 0
        do iwitn = 1,nwitn
           !
           ! last 2
           !
           if (interval >= 2 ) then
              aver(1,1) = avwit(iwitn,interval)
              aver(1,2) = avwit(iwitn,interval-1)
           else  ! so that it does not converge
              aver(1,1) =  1000.0_rp
              aver(1,2) = -1000.0_rp
           end if
           !
           ! last 4 (2 & 2)
           !
           if (interval >= 4 ) then
              aver(2,1) = ( avwit(iwitn,interval  ) + avwit(iwitn,interval-1) ) / 2.0_rp
              aver(2,2) = ( avwit(iwitn,interval-2) + avwit(iwitn,interval-3) ) / 2.0_rp
           else  ! so that it does not converge
              aver(2,1) =  1000.0_rp
              aver(2,2) = -1000.0_rp
           end if
           !
           ! last 8 (4 & 4)
           !
           if (interval >= 8 ) then
              aver(3,1) = ( avwit(iwitn,interval  ) + avwit(iwitn,interval-1) + avwit(iwitn,interval-2) + &
                   avwit(iwitn,interval-3)  ) / 4.0_rp
              aver(3,2) = ( avwit(iwitn,interval-4) + avwit(iwitn,interval-5) + avwit(iwitn,interval-6) + &
                   avwit(iwitn,interval-7)  ) / 4.0_rp
           else  ! so that it does not converge
              aver(3,1) =  1000.0_rp
              aver(3,2) = -1000.0_rp
           end if
           !
           ! Check if converged by any of the 3 criteria
           !
           iaux = 0
           do jaux = 1,3
              if ( abs ( 2.0_rp * ( aver(jaux,1) - aver(jaux,2) )  / ( aver(jaux,1) + aver(jaux,2) + toler*0.00001_rp ) ) < toler )  iaux = 1
           end do
           if ( iaux == 1 ) kount_converged = kount_converged + 1
        end do
        if ( kount_converged == nwitn ) then
           write(*,*) 'all witness have converged at step', ittim
           !
           ! Perform at least two time steps
           ! Example: coupling with temperature where initial temperature is constant
           !          and flow confined: it starts to move at second iteration.
           !
           if( ittim < 2 ) call runend('nsi_convergence_witness should not have converged that fast')
           kfl_stead_nsi = 1
           kfl_gotim = 0     ! I had to add this here because Output is called after Endste that is the one that controls gotim
           call outfor(28_ip,momod(modul) % lun_outpu,' ')
        else
           write(*,*) kount_converged, ' witness have converged at step', ittim
        end if
     end if

  end if
  !
  ! Sends to slaves so that everybody stops correctly  - had to send it to nsi_outwit so that everybody enters
  !

end subroutine nsi_convergence_witness

