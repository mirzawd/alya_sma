!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Neutro
!> @{
!> @file    neu_iniunk.f90
!> @date    31/03/2016
!> @author  Guillaume Houzeaux
!> @brief   Actualize unknown
!> @details Actualize unknown
!> @} 
!-----------------------------------------------------------------------

subroutine neu_updunk(itask)

  use def_kintyp
  use def_master
  use def_domain
  use def_neutro
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,idirection,ienergy
  real(rp)                :: rela1

  if( INOTMASTER ) then

     select case ( itask )

     case ( 1_ip ) 
        !
        ! Assign u(n,0,*) <-- u(n-1,*,*), initial guess for outer iterations
        !
        do idirection = 1,num_directions_neu
           do ienergy = 1,num_energies_neu
              do ipoin = 1,npoin
                 neutr(ienergy,idirection,ipoin,2) = neutr(ienergy,idirection,ipoin,nprev_neu) 
!                 write(6,*) ienergy,idirection, ipoin, neutr(ienergy,idirection,ipoin,2)
              end do
           end do
        end do


     case ( 2_ip )
        !
        ! Assign u(n,i,0) <-- u(n,i-1,*), initial guess for inner iterations
        !
        do idirection = 1,num_directions_neu
           do ienergy = 1,num_energies_neu
              do ipoin = 1,npoin
                 neutr(ienergy,idirection,ipoin,1) = neutr(ienergy,idirection,ipoin,2) 
              end do
           end do
        end do

     case ( 3_ip )
        !
        ! Assign u(n,i,j-1) <-- u(n,i,j), update of unknown after solver
        !
!        write(6,*) relax_neu,current_energy_neu,current_direction_neu
        rela1 = (1.0_rp-relax_neu)
        do ipoin = 1,npoin
           neutr(current_energy_neu,current_direction_neu,ipoin,1) = &
                &   rela1     * neutr(current_energy_neu,current_direction_neu,ipoin,1) &
                & + relax_neu * unkno(ipoin)
          
!           write(6,*) ipoin, neutr(ienergy,idirection,ipoin,1),unkno(ipoin)
        end do

     case ( 4_ip )
        !
        ! Assign u(n,i-1,*) <-- u(n,i,*)
        !        
        !
        do idirection = 1,num_directions_neu
           do ienergy = 1,num_energies_neu
              do ipoin = 1,npoin
                 neutr(ienergy,idirection,ipoin,2) = neutr(ienergy,idirection,ipoin,1) 
              end do
           end do
        end do
        
     case ( 5_ip )
        !
        ! u(n-1,*,*) <-- u(n,*,*)
        !
        !if( ADR_neu % kfl_time_integration /= 0 ) then
        !  
        !end if

     case ( 6_ip )
        !
        ! Initial guess for solver
        !
        do ipoin = 1,npoin
           unkno(ipoin) = neutr(current_energy_neu,current_direction_neu,ipoin,1) 
        end do

     case ( 11_ip )
        !
        ! Assign u(n,i,*)  <-- u(n-1,*,*), initial guess after iniunk or reading restart
        !        u'(n,i,*) <-- u'(n-1,*,*)
        !
        do idirection = 1,num_directions_neu
           do ienergy = 1,num_energies_neu
              do ipoin = 1,npoin
                 neutr(ienergy,idirection,ipoin,1) = neutr(ienergy,idirection,ipoin,nprev_neu)
                 neutr(ienergy,idirection,ipoin,2) = neutr(ienergy,idirection,ipoin,nprev_neu) 
              end do
           end do
        end do

     end select

  end if

end subroutine neu_updunk

