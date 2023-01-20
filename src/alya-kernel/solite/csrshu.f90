!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




subroutine csrshu(elmat,Aii,Aib,Abi,Abb,ndofn,pnode,mevat,lnode,iprob)
  !-----------------------------------------------------------------------
  !****f* solite/csrshu
  ! NAME
  !    csrshu
  ! DESCRIPTION
  !    This routine assemble the matrix using the bcsr format
  ! INPUT
  !    ELMAT
  !    NDOFN
  !    PNODE
  !    MEVAT
  !    LNODE
  !    IPROB
  ! OUTPUT
  !    AMATR
  ! USED BY
  !    ***_assmat
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_master, only       :  npoi1
  use def_domain, only       :  r_dom_aii,r_dom_aib,r_dom_abi,r_dom_abb
  use def_domain, only       :  c_dom_aii,c_dom_aib,c_dom_abi,c_dom_abb
  implicit none 
  integer(ip), intent(in)    :: iprob
  integer(ip), intent(in)    :: ndofn,pnode,mevat
  integer(ip), intent(in)    :: lnode(pnode)
  real(rp),    intent(in)    :: elmat(mevat,mevat)
  real(rp),    intent(inout) :: Aii(ndofn*ndofn,*)
  real(rp),    intent(inout) :: Aib(ndofn*ndofn,*)
  real(rp),    intent(inout) :: Abi(ndofn*ndofn,*) 
  real(rp),    intent(inout) :: Abb(ndofn*ndofn,*)
  integer(ip)                :: inode,jnode,kpoin,qpoin,ipoin2,jpoin2
  integer(ip)                :: ipoin,jpoin,izsol,jcolu

  if( iprob == 2 ) then

     if( ndofn == 1 ) then
        !
        ! General case: 1 unknown
        !   
        do inode = 1,pnode
           ipoin = lnode(inode)
           if( ipoin <= npoi1 ) then
              !ipoin2 = permr_Aii(ipoin)
              ipoin2 = ipoin
              do jnode = 1,pnode
                 jpoin = lnode(jnode)
                 if( jpoin <= npoi1 ) then
                    !jpoin2 = permr_Aii(jpoin)
                    jpoin2 = jpoin
                    izsol  = r_dom_Aii(ipoin2)
                    jcolu  = c_dom_Aii(izsol) 
                    do while( jcolu /= jpoin2 )
                       izsol = izsol + 1
                       jcolu = c_dom_aii(izsol) 
                       jcolu = c_dom_Aii(izsol)
                    end do

                    !$OMP ATOMIC
                    Aii(1,izsol) = Aii(1,izsol) + elmat(inode,jnode)

                 else
                    qpoin = jpoin - npoi1
                    izsol = r_dom_Aib(ipoin) 
                    jcolu = c_dom_Aib(izsol)
                    do while( jcolu /= qpoin )
                       izsol = izsol + 1
                       jcolu = c_dom_aib(izsol)
                    end do

                    !$OMP ATOMIC
                    Aib(1,izsol) = Aib(1,izsol) + elmat(inode,jnode)                    
                 end if
              end do
           else
              kpoin = ipoin - npoi1
              do jnode = 1,pnode
                 jpoin = lnode(jnode)
                 if( jpoin <= npoi1 ) then
                    izsol = r_dom_Abi(kpoin)
                    jcolu = c_dom_Abi(izsol)
                    do while( jcolu /= jpoin )
                       izsol = izsol + 1
                       jcolu = c_dom_abi(izsol)
                    end do

                    !$OMP ATOMIC
                    Abi(1,izsol) = Abi(1,izsol) + elmat(inode,jnode)
                 else
                    qpoin = jpoin - npoi1
                    izsol = r_dom_Abb(kpoin)
                    jcolu = c_dom_Abb(izsol)
                    do while( jcolu /= qpoin )
                       izsol = izsol + 1
                       jcolu = c_dom_abb(izsol)
                    end do

                    !$OMP ATOMIC
                    Abb(1,izsol) = Abb(1,izsol) + elmat(inode,jnode)                    
                 end if
              end do
           end if
        end do

     end if
  end if

end subroutine csrshu
