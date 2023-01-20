!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_store_tempe(itask)
  !-----------------------------------------------------------------------
  !****f* Temper/tem_store_tempe
  ! NAME 
  !    tem_store_tempe
  ! DESCRIPTION
  !    This routine performs stores the updated values of temperature
  !    after computed from the enthalpy.
  ! USED BY
  !    tem_endite (itask=1, inner loop) 
  !    tem_endste (itask=2)
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_temper
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: itime,ipoin,ielem,igaus,pelty,pgaus

  if( INOTMASTER ) then

     select case (itask)

     case(1_ip)

       if (kfl_lookg_tem >0) then
          ! 
          ! Store current Gauss point temperature (:,2) <=  (:,1)
          !
          do ielem = 1,nelem
             !
             ! Element dimensions
             !
             pelty = ltype(ielem)
             pgaus = ngaus(pelty)
      
             do igaus = 1,pgaus
               tempe_gp(ielem) % a(igaus,2,1) = tempe_gp(ielem) % a(igaus,1,1)
             enddo
          end do

          if(kfl_tisch_tem==3) &
              call runend('Temper store_tempe: lookup by Gauss points is not implemented for Adams-Bashforth (tisch=3)')
       else
          ! 
          ! Store current temperature (:,2) <=  (:,1)
          !
          !do ipoin = 1,npoin
          !   tempe(ipoin,3) = tempe(ipoin,2)
          !end do 
          do ipoin = 1,npoin
             tempe(ipoin,2) = tempe(ipoin,1)
          end do 
          if(kfl_tisch_tem==3) then
             do ipoin=1,npoin
                tempe(ipoin,4) = tempe(ipoin,3)
             end do 
             do ipoin=1,npoin
                tempe(ipoin,3) = tempe(ipoin,1)
             end do 
          end if
       endif

     case(2_ip)
       ! 
       ! High-order temporal schemes 
       !
       if(kfl_tisch_tem==2) then
          !
          ! BDF scheme
          !
          do itime=2+kfl_tiaor_tem,4,-1
             do ipoin=1,npoin
                tempe(ipoin,itime) = tempe(ipoin,itime-1)
             end do
          end do
       end if
       ! 
       ! Store current temperature (:,3) <=  (:,1)
       !
       if(kfl_tisch_tem/=3) then
          do ipoin=1,npoin
             tempe(ipoin,3) = tempe(ipoin,1)
          end do 
       end if
       if( kfl_regim_tem == 4_ip .and. kfl_lookg_tem > 0 ) then
          do ielem = 1,nelem
             pelty = ltype(ielem)
             pgaus = ngaus(pelty)
             do igaus=1,pgaus
                tempe_gp(ielem) % a(igaus,3,1) = tempe_gp(ielem) % a(igaus,1,1)
             end do
          end do
       end if
       
     end select

  end if

end subroutine tem_store_tempe
