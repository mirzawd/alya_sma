!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_detjac.f90
!> @author  Guillaume Houzeaux
!> @brief   Negative Jacobian detection 
!> @details Negative Jacobian detection during assembly
!>
!> @} 
!------------------------------------------------------------------------

subroutine sld_detjac(ielem_min,gpdet_min)

  use def_kintyp,            only : ip,rp
  use def_master,            only : displ
  use def_master,            only : INOTMASTER,kfl_paral
  use def_kermod,            only : kfl_detection               ! Event detection flag
  use mod_ker_detection,     only : ker_detection_min_max_value ! Event detection output
  use def_domain,            only : lnnod,lnods
  use def_domain,            only : ndime
  use def_domain,            only : coord
  use mod_communications,    only : PAR_MAX
  use def_solidz,            only : ledam_sld
  use def_solidz,            only : kfl_xfeme_sld
  use mod_messages,          only : livinf

  implicit none
  integer(ip), intent(in)    :: ielem_min
  real(rp),    intent(inout) :: gpdet_min

  integer(ip)                :: kelem,inode,ipoin
  real(rp)                   :: xcoor_min(3)

  kelem = ielem_min
  !
  ! When using event detection, all partitions should be aware something's happening
  !
  if( kfl_detection /= 0 ) then
     call PAR_MAX(kelem,'IN MY ZONE')      
  end if 

  if(( kelem > 0 ).and.( kfl_detection /= 0 )) then
     !
     ! A partition has detected a negative Jacobian:
     ! look for element IELEM_MIN with lowest Jacobian GPDET_MIN
     !
     if( kfl_detection /= 0 ) then
        !
        ! Event detection: it occurs at the center of gravity of the element
        !
        xcoor_min = 0.0_rp
        if( INOTMASTER .and. ielem_min > 0 ) then
           do inode = 1,lnnod(ielem_min)
              xcoor_min(1:ndime) = xcoor_min(1:ndime) + coord(1:ndime,lnods(inode,ielem_min))
           end do
           xcoor_min(1:ndime) = xcoor_min(1:ndime) / real(lnnod(ielem_min),rp)
        end if
        call ker_detection_min_max_value(&
             'MIN',gpdet_min,ndime,displ,xcoor_min,&
             'NEGATIVE_JACOBIAN','DISPL','AN_ELEMENT_HAS_NEGATIVE_JACOBIAN')
        call livinf(10000_ip,'A NEGATIVE JACOBIAN HAS BEEN DETECTED AS AN EVENT',0_ip)
        call runend('O.K.!')

     elseif ( kfl_detection == 0 .and. ledam_sld(kelem) == 0 ) then 
        ! 
        ! Print the information on screen
        !
        if( INOTMASTER .and. gpdet_min < 0.0_rp ) then
           write(*,*) '----> SLD_ELMope:  **ERROR** NEGATIVE GPDET. Check element ',ielem_min,' subdomain ',kfl_paral
           write(*,*) '      Value: ',gpdet_min
           write(*,*) '      Nodal coordinates: '
           if( kfl_xfeme_sld == 1 ) write(*,*) 'The element is enriched'
           do inode = 1,lnnod(ielem_min)
              ipoin = lnods(inode,ielem_min)
              write(*,*) '      Node=',inode,'  Position=',coord(1:ndime,ipoin)
           end do

           print*, "    ALYA:"
           print*, " Element:", ielem_min 
           print*, "   Nodes:", lnods(:,ielem_min)
           !             print *, " Coords:", coord(1:ndime,lnods(:,ielem_min))

        end if
        call runend('SLD_ELMOPE: NEGATIVE JACOBIAN(S) HAS(HAVE) BEEN FOUND')

     elseif ( kfl_detection == 0 .and. ledam_sld(kelem) > 0 ) then
        ! 
        ! Warning that there is a negative jacobian but coincides with a damaged element
        !
        if( INOTMASTER .and. gpdet_min < 0.0_rp ) then
           write(*,*) '----> SLD_ELMope:  **WARNING** NEGATIVE GPDET. Check damaged element ',ielem_min,' subdomain ',kfl_paral
           write(*,*) '      Value: ',gpdet_min
           write(*,*) '      Nodal coordinates: '

           print*, "    ALYA:"
           print*, " Element:", ielem_min
           print*, "   Nodes:", lnods(:,ielem_min)
           !             print *, " Coords:", coord(1:ndime,lnods(:,ielem_min))

        end if

     end if

  end if

end subroutine sld_detjac
