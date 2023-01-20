!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Neutro 
!> @{
!> @file    neu_doiter.f90
!> @date    29/03/2016
!> @author  Guillaume Houzeaux
!> @brief   Solve inner iterations
!> @details Solve inner iterations
!> @}
!------------------------------------------------------------------------

subroutine neu_doiter()

  use def_parame
  use def_master
  use def_solver
  use def_domain
  use def_neutro

  implicit none
  external :: neu_begite, neu_solite, neu_endite, neu_updheat

  if ( kfl_stead_neu == 0 ) then
      call neu_begite()

      if( INOTSLAVE .AND. ittim == 1 .AND. output_level_neu > 0) print *,''
      if( INOTSLAVE ) write(momod(modul)%lun_outpu,*) ''

      do while( kfl_goite_neu == 1 )

         itinn(modul) = itinn(modul) + 1

         if( INOTSLAVE .AND. output_level_neu == 1 ) &
            print '(A,I6,A,I5,A)',' Solving= iter ',itinn(modul),', ',num_energies_neu,' energy groups'
         
         residual_neu = 0.0_rp
         do current_energy_neu = num_energies_neu,1,-1 !num_energies_neu
            if( INOTSLAVE .AND. output_level_neu >= 2 ) &
               print '(A,I6,A,I5)',' Solving= iter ',itinn(modul),' energy group ',current_energy_neu
            
            ! hay que mantener siempre el grupo de mas alta energia
            ! if (current_energy_neu==num_energies_neu .OR. resid_energy_group_neu(current_energy_neu)/=0.0_rp) then
            if (resid_energy_group_neu(current_energy_neu)/=0.0_rp) then

               do current_direction_neu = 1,num_directions_neu
         !           if( INOTSLAVE ) print*,'Solving= ',itinn(modul),current_energy_neu,current_direction_neu
                  call neu_solite()
                  call neu_endite(1_ip)
               end do
            endif
            if( INOTSLAVE ) write(momod(modul)%lun_outpu,*) '              Grupo: ',&
                                 current_energy_neu,'      error= ',resid_energy_group_neu(current_energy_neu)
         end do
         call neu_endite(3_ip) ! Check global convergence
         if( INOTSLAVE .AND. output_level_neu >= 1 ) &
            print '(A,I6,A,ES15.6)',' Solving= iter ',itinn(modul),'           residual ',residual_neu
         if( INOTSLAVE ) write(momod(modul)%lun_outpu,*) '              Paso:  ',itinn(modul),'      error= ',residual_neu
         if( INOTSLAVE ) write(momod(modul)%lun_outpu,*) ''

         !   call neu_outvar(10_ip,itinn(modul))

      end do

      call neu_endite(2_ip)
      call neu_updheat()
   endif

end subroutine neu_doiter

