!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine elecsr()
  !-----------------------------------------------------------------------
  !****f* Domain/elecsr
  ! NAME
  !    elecsr
  ! DESCRIPTION
  !    This routine creates the arry from element to csr.
  !    Given ELMAT(INODE,JNODE,IELEM) => AMATR(IZDOM)
  !    with IZDOM = LZDOM(INODE,JNODE,IELEM)
  !    What to do when matrix is not referenced from position 1? 
  ! OUTPUT
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
  use def_kintyp,         only : ip
  use def_domain,         only : nelem,nboun
  use def_domain,         only : lnnod,lnodb
  use def_domain,         only : lnods,lnnob
  use def_domain,         only : lezdo,lbzdo
  use def_domain,         only : c_sol,r_sol
  use def_kermod,         only : kfl_element_to_csr
  use def_master,         only : INOTMASTER
  use mod_domain,         only : domain_memory_allocate
  implicit none
  integer(ip) :: ielem,pnode,ipoin,jnode,izsol,jcolu
  integer(ip) :: iboun,inodb,jnodb,pnodb,inode,jpoin

  if( INOTMASTER .and. kfl_element_to_csr == 1 ) then
     !
     ! Allocate memory
     !
     call domain_memory_allocate('LEZDO AND LBZDO')     
     !
     ! Element to CSR
     !
     do ielem = 1,nelem
        pnode = lnnod(ielem)
        do inode = 1,pnode
           ipoin = lnods(inode,ielem)
           do jnode = 1,pnode
              jpoin = lnods(jnode,ielem)
              izsol = r_sol(ipoin)
              jcolu = c_sol(izsol)
              do while( jcolu /= jpoin )
                 izsol = izsol + 1
                 jcolu = c_sol(izsol)
              end do
              lezdo(inode,jnode,ielem) = izsol
           end do
        end do
     end do
     !
     ! Boundary to CSR
     !
     do iboun = 1,nboun
        pnodb = lnnob(iboun)
        do inodb = 1,pnodb
           ipoin = lnodb(inodb,iboun)
           do jnodb = 1,pnodb
              jpoin = lnodb(jnodb,iboun)
              izsol = r_sol(ipoin)
              jcolu = c_sol(izsol)
              do while( jcolu /= jpoin )
                 izsol = izsol + 1
                 jcolu = c_sol(izsol)
              end do
              lbzdo(inodb,jnodb,iboun) = izsol
           end do
        end do
     end do

  end if

end subroutine elecsr
