!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine gather(itask,pgaus,pnode,ndime,lnods,gpsha,unkno,gpunk)
  !-----------------------------------------------------------------------
  !****f* mathru/gather
  ! NAME 
  !    gather
  ! DESCRIPTION
  !    Gather at Gauss point
  !    ITASK=1 ... From global array
  !         =2 ... From elemental array
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  implicit none
  integer(ip), intent(in)  :: itask,pgaus,pnode,ndime
  integer(ip), intent(in)  :: lnods(*)
  real(rp),    intent(in)  :: gpsha(pnode,pgaus),unkno(ndime,*)
  real(rp),    intent(out) :: gpunk(ndime,pgaus)
  integer(ip)              :: inode,idime,igaus,ipoin

  select case(itask)

  case(1)
     !
     ! Global to Gauss point
     !
     if(ndime==1) then
        !
        ! One-dimension
        !
        do igaus=1,pgaus
           gpunk(1,igaus)=0.0_rp
        end do
        do inode=1,pnode
           ipoin=lnods(inode)
           do igaus=1,pgaus
              gpunk(1,igaus)=gpunk(1,igaus)&
                   +unkno(1,ipoin)*gpsha(inode,igaus)
           end do
        end do
     else
        !
        ! Multidimension
        !
        do igaus=1,pgaus
           do idime=1,ndime
              gpunk(idime,igaus)=0.0_rp
           end do
        end do
        do inode=1,pnode
           ipoin=lnods(inode)
           do igaus=1,pgaus
              do idime=1,ndime
                 gpunk(idime,igaus)=gpunk(idime,igaus)&
                      +unkno(idime,ipoin)*gpsha(inode,igaus)
              end do
           end do
        end do
     end if

  case(2)
     !
     ! Elemental to Gauss point
     !
     if(ndime==1) then
        !
        ! One-dimension
        !
        do igaus=1,pgaus
           gpunk(1,igaus)=0.0_rp
        end do
        do inode=1,pnode
           do igaus=1,pgaus
              gpunk(1,igaus)=gpunk(1,igaus)&
                   +unkno(1,inode)*gpsha(inode,igaus)
           end do
        end do
     else
        !
        ! Multidimension
        !
        do igaus=1,pgaus
           do idime=1,ndime
              gpunk(idime,igaus)=0.0_rp
           end do
        end do
        do inode=1,pnode
           do igaus=1,pgaus
              do idime=1,ndime
                 gpunk(idime,igaus)=gpunk(idime,igaus)&
                      +unkno(idime,inode)*gpsha(inode,igaus)
              end do
           end do
        end do
     end if

  end select

end subroutine gather

