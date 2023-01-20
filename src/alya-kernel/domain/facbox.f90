!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine facbox(sabox,blink,bobox,nboun,candi,nlist)
  !-----------------------------------------------------------------------
  ! NAME
  !    ibm_fabbib
  ! DESCRIPTION
  !    This routines find the faces of a particle that intersect a bounding box
  !    INPUT 
  !       bobox: boundary box
  !       nboun: number of faces of the particle
  !    OUTPUT
  !       candi: face candidate list that intersect the boundary box
  !       nlist: number of candidate faces  
  ! USED BY
  !    contib 
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime
  implicit none  
  integer(ip), intent(in)  :: nboun
  real(rp),    intent(in)  :: sabox(2,ndime,2*nboun-1)
  integer(ip), intent(in)  :: blink(2*nboun-1)
  real(rp),    intent(in)  :: bobox(3,2)
  integer(ip), intent(out) :: candi(nboun),nlist
  integer(ip)              :: idime,cond,ilist
  integer(ip)              :: indst,curr
  integer(4)               :: istat
  integer(ip), pointer     :: struc(:) => null()

  !
  ! Tree structures of the particle
  !
  allocate( struc(2*nboun-1), stat=istat)

  struc(1) = 1_ip
  indst    = 1_ip
  ilist    = 0_ip

  do while( indst > 0 )
     curr  = struc(indst)
     indst = indst - 1_ip
     !
     ! Check if the current safety box intersect with the bounding box
     !
     cond = 0
     do idime = 1,ndime
         if ( bobox(idime,1) < sabox(2,idime,curr) .and. bobox(idime,2) > sabox(1,idime,curr) ) then
            cond = cond + 1_ip
         end if
     end do

     if( cond == ndime ) then
        !
        ! If ths current safety box is a face bounding box save as an intesection
        !  
        if( blink(curr) < 0 ) then
           ilist        =  ilist + 1_ip
           candi(ilist) = -blink(curr)           
        else
           indst        = indst + 1_ip
           struc(indst) = blink(curr)
           indst        = indst + 1_ip
           struc(indst) = blink(curr) + 1_ip
        end if

     end if
  end do

  nlist = ilist
  deallocate( struc, stat=istat)

end subroutine facbox
