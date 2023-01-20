!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    par_domlis.f90
!> @author  Guillaume Houzeaux
!> @brief   List of subdomain connected to a node
!> @details gives the list of subdomains DOMLI to which 
!!          node ipoin belongs. The number of subdomains is NDOMI.
!> @} 
!------------------------------------------------------------------------
subroutine par_domlis(pelpo,lelpo,ipoin,domain,ndomi,domli)
  use def_kintyp, only     :  ip,lg
  implicit none 
  integer(ip), intent(in)  :: ipoin                   !< Node to inquire
  integer(ip), intent(in)  :: pelpo(*)                !< Node to element graph: pointer
  integer(ip), intent(in)  :: lelpo(*)                !< Node to element graph: list
  integer(ip), intent(in)  :: domain(*)               !< List of element subdomains
  integer(ip), intent(out) :: ndomi                   !< Number of subdomains node IPOIN belongs to
  integer(ip), intent(out) :: domli(*)                !< List of subdomain IPOIN belongs to
  integer(ip)              :: ii,jj,ielem,domin
  logical(lg)              :: liste

  ndomi = 0
  do ii = pelpo(ipoin),pelpo(ipoin+1)-1
     ielem = lelpo(ii)
     domin = domain(ielem)
     jj    = 1
     liste = .true.
     do while( jj <= ndomi .and. liste )
        if( domli(jj) == domin ) then
           liste = .false.
        end if
        jj = jj + 1
     end do
     if( liste ) then
        ndomi        = ndomi + 1
        domli(ndomi) = domin
     end if
  end do

end subroutine par_domlis
