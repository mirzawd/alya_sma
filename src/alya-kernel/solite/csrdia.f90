!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine csrdia(ipoin,kfl_symme,izmat)
  !-----------------------------------------------------------------------
  !****f* master/csrdia
  ! NAME 
  !    csrdia
  ! DESCRIPTION
  !    Diagonal
  ! USES
  !    memchk
  !    mediso
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp,lg
  use def_domain, only       :  r_sym,r_dom,c_dom
  use def_master, only       :  INOTMASTER
  implicit none
  integer(ip), intent(in)    :: ipoin,kfl_symme
  integer(ip), intent(out)   :: izmat
  logical(lg)                :: ifoun

  if( INOTMASTER ) then

     if( kfl_symme == 1 ) then
       izmat = r_sym(ipoin+1)-1

     else
        izmat = r_dom(ipoin)-1
        ifoun = .true.
        do while( ifoun )
           izmat = izmat + 1
           if( c_dom(izmat) == ipoin ) ifoun = .false.
        end do
        
     end if

  end if
 
end subroutine csrdia
