!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_elmdir(&
     pnode,lnods,elmat,elrhs,ielem)
  !------------------------------------------------------------------------
  !****f* Temper/tem_elmdir
  ! NAME 
  !    tem_elmdir
  ! DESCRIPTION
  !    This routine prescribes the boundary conditions for the 
  !    temperature equations. 
  ! USES
  ! USED BY
  !    tem_elmope
  !    tem_bouope
  !------------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_temper, only       :  bvess_tem,kfl_fixno_tem
  use def_domain, only       :  nhang,lhang,lelch
  use def_elmtyp, only       :  ELHAN
  implicit none
  integer(ip), intent(in)    :: pnode,ielem
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(inout) :: elmat(pnode,pnode),elrhs(pnode)
  integer(ip)                :: inode,ipoin,jnode,knode,jpoin,ihang
  real(rp)                   :: adiag,xvalu
  
  if( lelch(ielem) == ELHAN ) then
     do inode = 1,pnode
        ipoin = lnods(inode)
        do ihang = 1,nhang
           if( lhang(ihang) % l(1) == ipoin ) then
              do jnode = 1,pnode
                 jpoin = lnods(jnode)
                 if( lhang(ihang) % l(2) == jpoin .or. lhang(ihang) % l(3) == jpoin ) then
                    elrhs(jnode) = elrhs(jnode) + elrhs(inode)
                    elrhs(inode) = 0.0_rp
                    do knode = 1,pnode
                       elmat(jnode,knode) = elmat(jnode,knode) + elmat(inode,knode)
                       elmat(inode,knode) = 0.0_rp
                    end do
                    do knode = 1,pnode
                       elmat(knode,jnode) = elmat(knode,jnode) + elmat(knode,inode)
                       elmat(knode,inode) = 0.0_rp
                    end do
                    elmat(inode,jnode) = -1.0_rp
                    elmat(inode,inode) =  1.0_rp
                 end if
              end do
           end if
        end do
     end do
  end if
  !if(ielem==2) then
  !   inode = 1
  !   ipoin = lnods(inode)
  !   do while( ipoin /= 6 )
  !      inode = inode + 1
  !      ipoin = lnods(inode)
  !   end do
  !   jnode = 1
  !   jpoin = lnods(jnode)
  !   do while( jpoin /= 3 )
  !      jnode = jnode + 1
  !      jpoin = lnods(jnode)
  !   end do
  !end if
  !if(ielem==3) then
  !   inode = 1
  !   ipoin = lnods(inode)
  !   do while( ipoin /= 6 )
  !      inode = inode + 1
  !      ipoin = lnods(inode)
  !   end do
  !   jnode = 1
  !   jpoin = lnods(jnode)
  !  do while( jpoin /= 9 )
  !      jnode = jnode + 1
  !      jpoin = lnods(jnode)
  !   end do
  !end if
  !if( ielem == 2 .or. ielem == 3 ) then
  !   elrhs(jnode) = elrhs(jnode) + elrhs(inode)
  !   elrhs(inode) = 0.0_rp
  !   do knode = 1,pnode
  !      elmat(jnode,knode) = elmat(jnode,knode) + elmat(inode,knode)
  !      elmat(inode,knode) = 0.0_rp
  !   end do
  !   do knode = 1,pnode
  !      elmat(knode,jnode) = elmat(knode,jnode) + elmat(knode,inode)
  !      elmat(knode,inode) = 0.0_rp
  !   end do
  !   elmat(inode,jnode) = -1.0_rp
  !   elmat(inode,inode) =  1.0_rp
  !end if

  do inode = 1,pnode
     ipoin = lnods(inode)
     if(  kfl_fixno_tem(1,ipoin) == 1 .or.&
          kfl_fixno_tem(1,ipoin) == 4 .or.&
          kfl_fixno_tem(1,ipoin) == 5 .or.&
          kfl_fixno_tem(1,ipoin) == 7 ) then
        adiag = elmat(inode,inode)
        xvalu = bvess_tem(1,ipoin,1)
        do jnode = 1,pnode
           elmat(inode,jnode) = 0.0_rp
           elrhs(jnode)       = elrhs(jnode) - elmat(jnode,inode) * xvalu
           elmat(jnode,inode) = 0.0_rp
        end do
        if( abs(adiag) > 0.0_rp ) then
           elmat(inode,inode) = adiag
           elrhs(inode)       = adiag * xvalu
        else
           elmat(inode,inode) = 1.0_rp
           elrhs(inode)       = xvalu           
        end if
     end if
  end do

end subroutine tem_elmdir
