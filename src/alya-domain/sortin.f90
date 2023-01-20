!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sortin(n,a)

!-----------------------------------------------------------------------
!
! Sort a vector
!
!-----------------------------------------------------------------------
  use       def_kintyp
  implicit none
  integer(ip) :: n,i,j,t
  integer(ip) :: a(n)
  
  if( n <= 1 ) return

  do i=1,n-1
     do j=i+1,n
        if (a(i)>a(j)) then
           t   =a(i)
           a(i)=a(j) 
           a(j)=t
        end if
     end do
  end do
  
end subroutine sortin

subroutine sorti3(n,a)

  !-----------------------------------------------------------------------
  !
  ! Sort a vector c according to a and then b
  !
  !-----------------------------------------------------------------------
  use       def_kintyp
  implicit none
  integer(ip) :: n,i,j,t1,t2,t3
  integer(ip) :: a(3,n)

  if( n <= 1 ) return
  
  do i = 1,n-1
     do j = i+1,n
        if( a(1,i) > a(1,j) ) then
           t1     = a(1,i)
           a(1,i) = a(1,j) 
           a(1,j) = t1
           t2     = a(2,i)
           a(2,i) = a(2,j) 
           a(2,j) = t2
           t3     = a(3,i)
           a(3,i) = a(3,j) 
           a(3,j) = t3
        end if
     end do
  end do

  do i = 1,n-1
     do j = i+1,n
        if( a(1,i) == a(1,j) ) then
           if( a(2,i) > a(2,j) ) then
              t1     = a(1,i)
              a(1,i) = a(1,j) 
              a(1,j) = t1
              t2     = a(2,i)
              a(2,i) = a(2,j) 
              a(2,j) = t2
              t3     = a(3,i)
              a(3,i) = a(3,j) 
              a(3,j) = t3
           end if
        end if
     end do
  end do

  do i = 1,n-1
     do j = i+1,n
        if( a(1,i) == a(1,j) .and. a(2,i) == a(2,j) ) then
           if( a(3,i) > a(3,j) ) then
              t1     = a(1,i)
              a(1,i) = a(1,j) 
              a(1,j) = t1
              t2     = a(2,i)
              a(2,i) = a(2,j) 
              a(2,j) = t2
              t3     = a(3,i)
              a(3,i) = a(3,j) 
              a(3,j) = t3
           end if
        end if
     end do
  end do

end subroutine sorti3
