!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



  !----------------------------------------------------------------------
  !> @addtogroup Nastin
  !> @{
  !> @file    nsi_velobl.f90
  !> @author  Guillaume Houzeaux
  !> @date    01/02/2014
  !> @brief   Interpolate velocity in boundary layer
  !> @details Interpolate velocity in boundary layer
  !> @} 
  !----------------------------------------------------------------------

subroutine nsi_test(itask)
  use def_kintyp,        only :  ip,rp
  use def_master,        only :  current_code,current_zone
  use def_master,        only :  veloc,INOTMASTER
  use def_domain,        only :  ndime
  use def_domain,        only :  npoin
  use def_domain,        only :  coord
  use mod_couplings, only :  COU_INIT_INTERPOLATE_POINTS_VALUES
  use mod_interpolation, only :  COU_GET_INTERPOLATE_POINTS_VALUES
  use def_coupli,        only :  typ_color_coupling
  use mod_parall,        only :  PAR_MY_CODE_RANK
  use mod_parall,        only :  par_code_zone_subd_to_color
  implicit none
  integer(ip), intent(in)     :: itask
  integer(ip)                 :: icolo,jcolo,pp,ii
  real(rp),    pointer        :: xcoor(:,:),xvalu(:,:)
  type(typ_color_coupling), save    :: coupling

  nullify(xcoor) 
  nullify(xvalu) 
  icolo            = par_code_zone_subd_to_color(2_ip,current_zone,0_ip)
  jcolo            = par_code_zone_subd_to_color(1_ip,current_zone,0_ip)
  coupling % itype = 1
  if( current_code == 2 .and. INOTMASTER ) then
     allocate(xcoor(ndime,npoin))
     xcoor(1:ndime,1:npoin) = coord(1:ndime,1:npoin)
     allocate(xvalu(ndime,npoin))
  end if

!!!!  if( itask == 1 ) call COU_INIT_INTERPOLATE_POINTS_VALUES(xcoor,icolo,jcolo,coupling)    ! OJO

  if( itask == 2 .and. current_code == 1 ) then
     call COU_GET_INTERPOLATE_POINTS_VALUES(veloc,xvalu,coupling)
  else if( itask == 1 .and. current_code == 2 ) then
     call COU_GET_INTERPOLATE_POINTS_VALUES(veloc,xvalu,coupling)     
     if( associated(xvalu) ) veloc(:,:,1) = xvalu
  end if

  if( associated(xvalu) ) deallocate( xvalu )
  return

  if( PAR_MY_CODE_RANK /= 0 ) then
     veloc(1,1:npoin,1) = coord(1,1:npoin)
     veloc(2,1:npoin,1) = coord(2,1:npoin)
  end if

  nullify(xcoor)
  nullify(xvalu)
  pp = 0

  if( PAR_MY_CODE_RANK == 1 ) then
     pp = 3
     allocate( xcoor(ndime,pp) )
     xcoor(1,1) = 0.7_rp
     xcoor(2,1) = 0.7_rp
     xcoor(1,2) = 2.1_rp
     xcoor(2,2) = 2.9_rp
     xcoor(1,3) = 0.5_rp
     xcoor(2,3) = 0.5_rp
  end if
  if( PAR_MY_CODE_RANK == 2 ) then
     pp = 2
     allocate( xcoor(ndime,pp) )
     xcoor(1,1) = 0.3_rp
     xcoor(2,1) = 0.7_rp
     xcoor(1,2) = 1.5_rp
     xcoor(2,2) = 1.5_rp
  end if
  if( pp > 1 ) then
     allocate(xvalu(ndime,pp))
     xvalu  = 0.0_rp
  end if
  !
  ! Intialize interpolation
  !
  icolo = par_code_zone_subd_to_color(current_code,current_zone,0_ip)
  jcolo = par_code_zone_subd_to_color(current_code,current_zone,0_ip)
  coupling % itype = 1
!!!  call COU_INIT_INTERPOLATE_POINTS_VALUES(xcoor,icolo,jcolo,coupling)    ! OJO
  !
  ! Interpolate
  !
  call COU_GET_INTERPOLATE_POINTS_VALUES(veloc,xvalu,coupling)

  do ii = 1,pp
     if( coupling % geome % status(ii) /= 0 ) then
        print*,ii,PAR_MY_CODE_RANK,xcoor(1,ii)-xvalu(1,ii),xcoor(2,ii)-xvalu(2,ii)
     else
        print*,ii,PAR_MY_CODE_RANK,' is lost'
     end if
  end do

  if( associated(xcoor) ) deallocate( xcoor )
  if( associated(xvalu) ) deallocate( xvalu )

  call runend('O.K.!')

end subroutine nsi_test
