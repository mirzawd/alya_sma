!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mergl3( ipoin, lista, lsize, lnods, pnode )
  !-----------------------------------------------------------------------
  !****f* Domain/mergl2
  ! NAME
  !    mergl2
  ! DESCRIPTION
  !    This routine merges to list of lnods
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
  integer(ip), intent(in)    :: ipoin,pnode
  integer(ip), intent(in)    :: lnods(pnode)
  integer(ip)                :: inode, jj, n1, n2
  logical(lg)                :: noEncontrado

  do inode = 1, pnode
     n1 = lnods(inode)
     if( n1 < ipoin ) then
        jj = 1
        noEncontrado = .true.
        do while( jj <= lsize .and. noEncontrado )
           n2 = lista(jj)
           if( n1 == n2 ) then
              noEncontrado = .false.
           end if
           jj = jj + 1
        end do
        if (noEncontrado) then
           lsize = lsize + 1
           lista(lsize) = n1
        end if
     end if
  end do

end subroutine mergl3

subroutine mergl4( ipoin, lista, lsize, lnods, pnode )
  !-----------------------------------------------------------------------
  !****f* Domain/mergl2
  ! NAME
  !    mergl2
  ! DESCRIPTION
  !    This routine merges to list of lnods using only edges
  ! OUTPUT
  !    LISTA
  !    LSIZE
  ! USED BY
  !    domgra
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,lg
  use def_domain, only       :  ndime
  implicit none
  integer(ip), intent(inout) :: lsize,lista(*)
  integer(ip), intent(in)    :: ipoin,pnode
  integer(ip), intent(in)    :: lnods(pnode)
  integer(ip)                :: inode, jj, n1, n2
  integer(ip)                :: lnod1(10), knode
  logical(lg)                :: noEncontrado

  if( ndime == 2 .and. pnode == 4 ) then
     if( lnods(1) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(2)
        lnod1(3) = lnods(4)
     else if( lnods(2) == ipoin ) then
        lnod1(1) = lnods(1)
        lnod1(2) = ipoin
        lnod1(3) = lnods(3)
     else if( lnods(3) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(2)
        lnod1(3) = lnods(4)
     else if( lnods(4) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(1)
        lnod1(3) = lnods(3)
     end if
     do inode = 1,3
        n1 = lnod1(inode)
        if( n1 < ipoin ) then
           jj = 1
           noEncontrado = .true.
           do while( jj <= lsize .and. noEncontrado )
              n2 = lista(jj)
              if( n1 == n2 ) then
                 noEncontrado = .false.
              end if
              jj = jj + 1
           end do
           if (noEncontrado) then
              lsize = lsize + 1
              lista(lsize) = n1
           end if
        end if
     end do
     
  else if( ndime == 3 .and. pnode == 8 ) then
     
     if( lnods(1) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(2)
        lnod1(3) = lnods(4)
        lnod1(4) = lnods(5)
     else if( lnods(2) == ipoin ) then
        lnod1(1) = lnods(1)
        lnod1(2) = ipoin
        lnod1(3) = lnods(3)
        lnod1(4) = lnods(6)
     else if( lnods(3) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(2)
        lnod1(3) = lnods(4)
        lnod1(4) = lnods(7)
     else if( lnods(4) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(1)
        lnod1(3) = lnods(3)
        lnod1(4) = lnods(8)
     else if( lnods(5) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(1)
        lnod1(3) = lnods(6)
        lnod1(4) = lnods(8)
     else if( lnods(6) == ipoin ) then
        lnod1(1) = lnods(5)
        lnod1(2) = ipoin
        lnod1(3) = lnods(7)
        lnod1(4) = lnods(2)
     else if( lnods(7) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(6)
        lnod1(3) = lnods(8)
        lnod1(4) = lnods(3)
     else if( lnods(8) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(5)
        lnod1(3) = lnods(7)
        lnod1(4) = lnods(4)
     end if
     do inode = 1,4
        n1 = lnod1(inode)
        if( n1 < ipoin ) then
           jj = 1
           noEncontrado = .true.
           do while( jj <= lsize .and. noEncontrado )
              n2 = lista(jj)
              if( n1 == n2 ) then
                 noEncontrado = .false.
              end if
              jj = jj + 1
           end do
           if (noEncontrado) then
              lsize = lsize + 1
              lista(lsize) = n1
           end if
        end if
     end do 
     
  else if( ndime == 3 .and. pnode == 6 ) then
     
     if( lnods(1) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(2)
        lnod1(3) = lnods(3)
        lnod1(4) = lnods(4)
     else if( lnods(2) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(1)
        lnod1(3) = lnods(3)
        lnod1(4) = lnods(5)
     else if( lnods(3) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(1)
        lnod1(3) = lnods(2)
        lnod1(4) = lnods(6)
     else if( lnods(4) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(1)
        lnod1(3) = lnods(5)
        lnod1(4) = lnods(6)
     else if( lnods(5) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(2)
        lnod1(3) = lnods(4)
        lnod1(4) = lnods(6)
     else if( lnods(6) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(3)
        lnod1(3) = lnods(4)
        lnod1(4) = lnods(5)
     end if
     do inode = 1,4
        n1 = lnod1(inode)
        if( n1 < ipoin ) then
           jj = 1
           noEncontrado = .true.
           do while( jj <= lsize .and. noEncontrado )
              n2 = lista(jj)
              if( n1 == n2 ) then
                 noEncontrado = .false.
              end if
              jj = jj + 1
           end do
           if (noEncontrado) then
              lsize = lsize + 1
              lista(lsize) = n1
           end if
        end if
     end do 
     
  else if( ndime == 3 .and. pnode == 5 ) then
     
     knode = 4
     if( lnods(1) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(2)
        lnod1(3) = lnods(4)
        lnod1(4) = lnods(5)
     else if( lnods(2) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(1)
        lnod1(3) = lnods(3)
        lnod1(4) = lnods(5)
     else if( lnods(3) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(2)
        lnod1(3) = lnods(4)
        lnod1(4) = lnods(5)
     else if( lnods(4) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(1)
        lnod1(3) = lnods(3)
        lnod1(4) = lnods(5)
     else if( lnods(5) == ipoin ) then
        knode    = 5
        lnod1(1) = ipoin
        lnod1(2) = lnods(1)
        lnod1(3) = lnods(2)
        lnod1(4) = lnods(3)
        lnod1(5) = lnods(4)
     end if
     do inode = 1,knode
        n1 = lnod1(inode)
        if( n1 < ipoin ) then
           jj = 1
           noEncontrado = .true.
           do while( jj <= lsize .and. noEncontrado )
              n2 = lista(jj)
              if( n1 == n2 ) then
                 noEncontrado = .false.
              end if
              jj = jj + 1
           end do
           if (noEncontrado) then
              lsize = lsize + 1
              lista(lsize) = n1
           end if
        end if
     end do 
     
  else

     do inode = 1, pnode
        n1 = lnods(inode)
        if( n1 < ipoin ) then
           jj = 1
           noEncontrado = .true.
           do while( jj <= lsize .and. noEncontrado )
              n2 = lista(jj)
              if( n1 == n2 ) then
                 noEncontrado = .false.
              end if
              jj = jj + 1
           end do
           if (noEncontrado) then
              lsize = lsize + 1
              lista(lsize) = n1
           end if
        end if
     end do
  end if

end subroutine mergl4

subroutine mergl5( ipoin, lista, lsize, lnods, pnode )
  !-----------------------------------------------------------------------
  !****f* Domain/mergl2
  ! NAME
  !    mergl2
  ! DESCRIPTION
  !    This routine merges to list of lnods
  ! OUTPUT
  !    LISTA
  !    LSIZE
  ! USED BY
  !    domgra
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,lg
  use def_domain, only       :  ndime
  implicit none
  integer(ip), intent(inout) :: lsize,lista(*)
  integer(ip), intent(in)    :: ipoin,pnode
  integer(ip), intent(in)    :: lnods(pnode)
  integer(ip)                :: inode, jj, n1, n2
  integer(ip)                :: lnod1(10), knode
  logical(lg)                :: noEncontrado

  if( ndime == 2 .and. pnode == 4 ) then
     if( lnods(1) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(2)
        lnod1(3) = lnods(4)
     else if( lnods(2) == ipoin ) then
        lnod1(1) = lnods(1)
        lnod1(2) = ipoin
        lnod1(3) = lnods(3)
     else if( lnods(3) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(2)
        lnod1(3) = lnods(4)
     else if( lnods(4) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(1)
        lnod1(3) = lnods(3)
     end if
     do inode = 1,3
        n1 = lnod1(inode)
        if( n1 /= ipoin ) then
           jj = 1
           noEncontrado = .true.
           do while( jj <= lsize .and. noEncontrado )
              n2 = lista(jj)
              if( n1 == n2 ) then
                 noEncontrado = .false.
              end if
              jj = jj + 1
           end do
           if (noEncontrado) then
              lsize = lsize + 1
              lista(lsize) = n1
           end if
        end if
     end do
     
  else if( ndime == 3 .and. pnode == 8 ) then
     
     if( lnods(1) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(2)
        lnod1(3) = lnods(4)
        lnod1(4) = lnods(5)
     else if( lnods(2) == ipoin ) then
        lnod1(1) = lnods(1)
        lnod1(2) = ipoin
        lnod1(3) = lnods(3)
        lnod1(4) = lnods(6)
     else if( lnods(3) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(2)
        lnod1(3) = lnods(4)
        lnod1(4) = lnods(7)
     else if( lnods(4) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(1)
        lnod1(3) = lnods(3)
        lnod1(4) = lnods(8)
     else if( lnods(5) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(1)
        lnod1(3) = lnods(6)
        lnod1(4) = lnods(8)
     else if( lnods(6) == ipoin ) then
        lnod1(1) = lnods(5)
        lnod1(2) = ipoin
        lnod1(3) = lnods(7)
        lnod1(4) = lnods(2)
     else if( lnods(7) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(6)
        lnod1(3) = lnods(8)
        lnod1(4) = lnods(3)
     else if( lnods(8) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(5)
        lnod1(3) = lnods(7)
        lnod1(4) = lnods(4)
     end if
     do inode = 1,4
        n1 = lnod1(inode)
        if( n1 /= ipoin ) then
           jj = 1
           noEncontrado = .true.
           do while( jj <= lsize .and. noEncontrado )
              n2 = lista(jj)
              if( n1 == n2 ) then
                 noEncontrado = .false.
              end if
              jj = jj + 1
           end do
           if (noEncontrado) then
              lsize = lsize + 1
              lista(lsize) = n1
           end if
        end if
     end do 
     
  else if( ndime == 3 .and. pnode == 6 ) then
     
     if( lnods(1) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(2)
        lnod1(3) = lnods(3)
        lnod1(4) = lnods(4)
     else if( lnods(2) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(1)
        lnod1(3) = lnods(3)
        lnod1(4) = lnods(5)
     else if( lnods(3) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(1)
        lnod1(3) = lnods(2)
        lnod1(4) = lnods(6)
     else if( lnods(4) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(1)
        lnod1(3) = lnods(5)
        lnod1(4) = lnods(6)
     else if( lnods(5) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(2)
        lnod1(3) = lnods(4)
        lnod1(4) = lnods(6)
     else if( lnods(6) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(3)
        lnod1(3) = lnods(4)
        lnod1(4) = lnods(5)
     end if
     do inode = 1,4
        n1 = lnod1(inode)
        if( n1 /= ipoin ) then
           jj = 1
           noEncontrado = .true.
           do while( jj <= lsize .and. noEncontrado )
              n2 = lista(jj)
              if( n1 == n2 ) then
                 noEncontrado = .false.
              end if
              jj = jj + 1
           end do
           if (noEncontrado) then
              lsize = lsize + 1
              lista(lsize) = n1
           end if
        end if
     end do 
     
  else if( ndime == 3 .and. pnode == 5 ) then
     
     knode = 4
     if( lnods(1) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(2)
        lnod1(3) = lnods(4)
        lnod1(4) = lnods(5)
     else if( lnods(2) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(1)
        lnod1(3) = lnods(3)
        lnod1(4) = lnods(5)
     else if( lnods(3) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(2)
        lnod1(3) = lnods(4)
        lnod1(4) = lnods(5)
     else if( lnods(4) == ipoin ) then
        lnod1(1) = ipoin
        lnod1(2) = lnods(1)
        lnod1(3) = lnods(3)
        lnod1(4) = lnods(5)
     else if( lnods(5) == ipoin ) then
        knode    = 5
        lnod1(1) = ipoin
        lnod1(2) = lnods(1)
        lnod1(3) = lnods(2)
        lnod1(4) = lnods(3)
        lnod1(5) = lnods(4)
     end if
     do inode = 1,knode
        n1 = lnod1(inode)
        if( n1 /= ipoin ) then
           jj = 1
           noEncontrado = .true.
           do while( jj <= lsize .and. noEncontrado )
              n2 = lista(jj)
              if( n1 == n2 ) then
                 noEncontrado = .false.
              end if
              jj = jj + 1
           end do
           if (noEncontrado) then
              lsize = lsize + 1
              lista(lsize) = n1
           end if
        end if
     end do 
     
  else

     do inode = 1, pnode
        n1 = lnods(inode)
        if( n1 /= ipoin ) then
           jj = 1
           noEncontrado = .true.
           do while( jj <= lsize .and. noEncontrado )
              n2 = lista(jj)
              if( n1 == n2 ) then
                 noEncontrado = .false.
              end if
              jj = jj + 1
           end do
           if (noEncontrado) then
              lsize = lsize + 1
              lista(lsize) = n1
           end if
        end if
     end do
  end if

end subroutine mergl5
