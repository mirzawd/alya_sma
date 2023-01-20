!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mergli( lista, lsize, nodes, nnode, me )
!-----------------------------------------------------------------------
!****f* Domain/mergli
! NAME
!    mergli
! DESCRIPTION
!    This routine merges to list of nodes
! OUTPUT
!    LISTA
!    LSIZE
! USED BY
!    domgra
!***
!-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,lg
  implicit none
  integer(ip), intent(inout) :: lsize,lista(*)
  integer(ip), intent(in)    :: nnode,me
  integer(ip), intent(in)    :: nodes(nnode)
  integer(ip)                :: ii, jj, n1, n2
  logical(lg)                :: noEncontrado

  if( me < 0 ) then
    do ii = 1,nnode
      n1 = nodes(ii)
      jj = 1
      noEncontrado = .true.
      do while( jj <= lsize .and. noEncontrado )
        n2 = lista(jj)
        if ( n1 == n2 ) then
          noEncontrado = .false.
        end if
        jj = jj + 1
      end do
      if ( noEncontrado ) then
        lsize = lsize + 1
        lista(lsize) = n1
      end if
    end do

  else

    do ii = 1, nnode
      n1 = nodes(ii)
      if( n1 /= me ) then
        jj = 1
        noEncontrado = .true.
        do while( jj <= lsize .and. noEncontrado )
          n2 = lista(jj)
          if( n1 == n2 ) then
            noEncontrado = .false.
          end if
          jj = jj + 1
        end do
        if( noEncontrado ) then
          lsize = lsize + 1
          lista(lsize) = n1
        end if
      end if
    end do
  end if

end subroutine mergli
