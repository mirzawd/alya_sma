!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine soldod(nbvar,an,bb,xx)
  !-----------------------------------------------------------------------
  !****f* solite/soldod
  ! NAME 
  !    soldod
  ! DESCRIPTION
  !    Take off nodes out of my zone
  !    Take off hole nodes
  ! USES
  ! USED BY
  !    nsi_dodem1
  !***
  !----------------------------------------------------------------------- 
  use def_kintyp, only       :  ip,rp,lg
  use def_master, only       :  INOTMASTER
  use def_domain, only       :  npoin,c_dom,r_dom,nzone
  use def_elmtyp, only       :  NOHOL
  implicit none
  integer(ip), intent(in)    :: nbvar
  real(rp),    intent(inout) :: an(nbvar,nbvar,*)
  real(rp),    intent(inout) :: bb(nbvar,*)
  real(rp),    intent(inout) :: xx(nbvar,*)
  integer(ip)                :: ipoin,izdom,jpoin,kk,ll
  logical(lg), pointer       :: not_in_my_zone(:)
  !
  ! Cancel lines for nodes not in current zone
  !
  if( nzone > 1 .and. INOTMASTER ) then

     allocate( not_in_my_zone(npoin) )
     do ipoin = 1,npoin
        not_in_my_zone(ipoin) = .true.
     end do
     do ipoin = 1,npoin
        not_in_my_zone(ipoin) = .false.
     end do

     do ipoin = 1,npoin
        if( not_in_my_zone(ipoin) ) then
           do kk = 1,nbvar
              bb(kk,ipoin) = 0.0_rp
              xx(kk,ipoin) = 0.0_rp
           end do
           do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
              jpoin = c_dom(izdom)
              do kk = 1,nbvar
                 do ll = 1,nbvar
                    an(ll,kk,izdom) = 0.0_rp
                 end do
              end do
              if( ipoin == jpoin ) then                 
                 do kk = 1,nbvar
                    an(kk,kk,izdom) = 0.0_rp
                 end do
              end if
           end do
        end if
     end do

     deallocate( not_in_my_zone )

  end if

end subroutine soldod
