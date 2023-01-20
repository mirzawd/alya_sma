!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_elmdir(itask,pnode,pevat,lnods,elmat,elrhs)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_elmdir
  ! NAME
  !   tur_elmdir
  ! DESCRIPTION
  ! This routine prescribes the boundary conditions for the 
  ! turbulence equations. 
  ! USES
  ! USED BY
  !    tur_elmop1
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
!  use def_domain, only       :  lmatn
  use def_kermod, only       :  kfl_adj_prob
  use def_turbul, only       :  nturb_tur,kfl_fixno_tur,bvess_tur,&
       &                        iunkn_tur,kfl_algor_tur
  implicit none
  integer(ip), intent(in)    :: itask,pnode,pevat
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(inout) :: elmat(pevat,pevat),elrhs(pevat)
  integer(ip)                :: inode,ipoin,jnode,junkn_tur,ievat,jevat
  real(rp)                   :: adiag,turdi

  if( kfl_algor_tur == 1 ) then

     if( itask /= 3 ) then
        !
        ! Fluid element
        !
        do inode = 1,pnode
           ipoin = lnods(inode)
           if(  kfl_fixno_tur(1,ipoin,iunkn_tur) > 0 ) then
              turdi = bvess_tur(1,ipoin,iunkn_tur)
              adiag = elmat(inode,inode)
              do jnode = 1,pnode
                 elmat(inode,jnode) = 0.0_rp
                 !elrhs(jnode)=elrhs(jnode)-elmat(jnode,inode)*turdi
!                  if (kfl_adj_prob == 1) elrhs(jnode)=elrhs(jnode)-elmat(jnode,inode)*turdi
                 if (kfl_adj_prob == 1) elmat(jnode,inode)=0.0_rp
              end do
              elmat(inode,inode) = adiag
              elrhs(inode)       = adiag*turdi
           end if
        end do

     else
        !
        ! Solid element
        !
        do inode=1,pnode
           do jnode=1,pnode
              elmat(inode,jnode)=0.0_rp
              !elmat(jnode,inode)=0.0_rp
           end do
           elrhs(inode)=0.0_rp
        end do
        do inode=1,pnode
           ipoin=lnods(inode)
           !if(lmatn(ipoin)==-1) elmat(inode,inode)=1.0_rp
        end do

     end if

  else if(kfl_algor_tur==2) then

     if(itask/=3) then
        !
        ! Fluid element
        !
        do inode=1,pnode
           ipoin=lnods(inode)
           do iunkn_tur=1,nturb_tur
              if(  kfl_fixno_tur(1,ipoin,iunkn_tur)>=1) then
                 turdi=bvess_tur(1,ipoin,iunkn_tur)
                 ievat=(inode-1)*nturb_tur+iunkn_tur
                 adiag=elmat(ievat,ievat)
                 do jnode=1,pnode
                    do junkn_tur=1,nturb_tur
                       jevat=(jnode-1)*nturb_tur+junkn_tur
                       elmat(ievat,jevat)=0.0_rp
                       !elrhs(jevat)=elrhs(jevat)-elmat(jevat,ievat)*turdi
                       !elmat(jevat,ievat)=0.0_rp
                    end do
                 end do
                 elmat(ievat,ievat)=adiag
                 elrhs(ievat)=adiag*turdi
              end if
           end do
        end do

     else
        !
        ! Solid element
        !
        do ievat=1,pevat
           do jevat=1,pevat
              elmat(ievat,jevat)=0.0_rp
              elmat(jevat,ievat)=0.0_rp
           end do
           elrhs(ievat)=0.0_rp
        end do
        do inode=1,pnode
           ipoin=lnods(inode)
           !if(lmatn(ipoin)==-1) then
           !   ievat=(inode-1)*nturb_tur
           !   do iunkn_tur=1,nturb_tur
           !      ievat=ievat+1
           !      elmat(ievat,ievat)=1.0_rp
           !   end do
           !end if
        end do

     end if

  end if

end subroutine tur_elmdir
