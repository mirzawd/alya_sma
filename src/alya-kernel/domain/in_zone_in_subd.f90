!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    in_zine_in_subd.f90
!> @author  Guillaume Houzeaux
!> @brief   In zone or in subd
!> @details Define if i am in a zone and in a subdomain
!>          I am in a zone or a subdomain if I have at least one
!>          node in the zone or in the subdomain
!>
!> @} 
!-----------------------------------------------------------------------
subroutine in_zone_in_subd()
  use def_kintyp,         only :  ip,rp,lg
  use def_master,         only :  I_AM_IN_ZONE
  use def_master,         only :  I_AM_IN_SUBD
  use def_master,         only :  IMASTER,kfl_paral,leinv_loc
  use def_domain,         only :  nsubd,nzone 
  use def_domain,         only :  npoin
  use def_domain,         only :  nelem
  use def_domain,         only :  lesub,ltype
  use def_domain,         only :  lnods
  use def_domain,         only :  lnnod
  use mod_communications, only :  PAR_INTERFACE_NODE_EXCHANGE
  use mod_domain,         only :  domain_memory_allocate

  implicit none
  integer(ip)                  :: ipoin,izone,isubd,ielem,inode
  integer(ip), pointer         :: in_out(:)
  !
  ! Allocate memory
  !
  call domain_memory_allocate('I_AM_IN_ZONE AND I_AM_IN_SUBD')

  if( IMASTER ) then
     !
     ! Master is always in all zones and subdomains
     !
     I_AM_IN_ZONE(0:nzone) = .true.
     I_AM_IN_SUBD(0:nsubd) = .true.
  else 

     nullify( in_out )
     allocate( in_out(npoin) )
     do izone = 1,nzone
        I_AM_IN_ZONE(izone) = .false.
     end do
     do isubd = 1,nsubd
        I_AM_IN_SUBD(isubd) = .false.
     end do
     I_AM_IN_ZONE(0) = .true.
     I_AM_IN_SUBD(0) = .true.
     !
     ! In zone
     !
     izone = 1
     I_AM_IN_ZONE(izone) = .true.
     in_out = 0
     !
     ! In subdomains
     !    
     do isubd = 1,nsubd
        do ielem = 1,nelem
           if( lesub(ielem) == isubd ) then
              if(lnnod(ielem)==0) then
                 print*,'error in in_zone_in_subd.f90=',ielem,kfl_paral,leinv_loc(ielem),ltype(ielem),lnnod(ielem)
              end if
              do inode = 1,lnnod(ielem)
                 ipoin = lnods(inode,ielem)
                 in_out(ipoin) = 1
              end do
           end if
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(in_out,'SUM','IN MY CODE','SYNCHRONOUS')      
        ipoin = 0
        do while( ipoin < npoin )
           ipoin = ipoin + 1
           if( in_out(ipoin) > 0 ) then
              I_AM_IN_SUBD(isubd) = .true.
              in_out = 0
              ipoin  = npoin
           end if
        end do
        
     end do
     deallocate( in_out )

  end if

end subroutine in_zone_in_subd
