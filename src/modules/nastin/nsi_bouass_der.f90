!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_bouass_der(&
     pevat,ndofn,pnodb,pnode,lboel,setfp_derp,&
     setfp_derx,elrhs,eldcost_dx)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_bouass_der
  ! NAME 
  !    nsi_bouass_der
  ! DESCRIPTION
  ! USES
  ! USED BY
  !    nsi_matrix
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime
  use def_master, only     :  vbset
  implicit none
  integer(ip), intent(in)  :: pevat,ndofn,pnodb,pnode
  integer(ip), intent(in)  :: lboel(pnodb) 
  real(rp),    intent(in)  :: setfp_derp(ndime,pnodb),setfp_derx(ndime,ndime,pnodb)
  real(rp),    intent(out) :: elrhs(pevat),eldcost_dx(pnode*ndime)
  integer(ip)              :: inodb,idime,ievab
  real(rp)                 :: cd,cl
  
!   cd = (vbset(6,1) + vbset(6,2) + vbset(6,3) + vbset(6,4))
!   cl = (vbset(8,1) + vbset(8,2) + vbset(8,3) + vbset(8,4))

     cl = (0.96_rp*(vbset(7,2)+vbset(7,3)+vbset(7,4)) - 0.28_rp*(vbset(6,2)+vbset(6,3)+vbset(6,4)))
     cd = (0.96_rp*(vbset(6,2)+vbset(6,3)+vbset(6,4)) + 0.28_rp*(vbset(7,2)+vbset(7,3)+vbset(7,4)))
     
!      print *, cd, cl
     
!      -10*cd*cd + (cl+2495.0_rp)*(cl+2495.0_rp)

     
  do inodb = 1,pnodb
     ievab = lboel(inodb) * (ndime + 1) 
     elrhs(ievab) = elrhs(ievab) + setfp_derp(1,inodb)
!      elrhs(ievab) = elrhs(ievab) - (0.96_rp*setfp_derp(1,inodb) + 0.28_rp*setfp_derp(2,inodb)) 
!      elrhs(ievab) = elrhs(ievab) + 20*cd*(0.96_rp*setfp_derp(1,inodb) + 0.28_rp*setfp_derp(2,inodb)) &
!                                  + 2*(cl+2491.9730_rp)*(0.96_rp*setfp_derp(2,inodb) - 0.28_rp*setfp_derp(1,inodb))
     ievab = (lboel(inodb)-1) * ndime
     do idime = 1,ndime
        ievab = ievab + 1
        eldcost_dx(ievab) = eldcost_dx(ievab) + setfp_derx(1,idime,inodb)
!         eldcost_dx(ievab) = eldcost_dx(ievab) - (0.96_rp*setfp_derx(1,idime,inodb) + 0.28_rp*setfp_derx(2,idime,inodb))
!         eldcost_dx(ievab) = eldcost_dx(ievab) + 20*cd*(0.96_rp*setfp_derx(1,idime,inodb) + 0.28_rp*setfp_derx(2,idime,inodb)) + &
!                                                 2*(cl+2491.9730_rp)*(0.96_rp*setfp_derx(2,idime,inodb) - 0.28_rp*setfp_derx(1,idime,inodb))
     end do
  end do

  
!   do inode = 1,pnode
!      do idime = 1,ndime
!         ievae = (inode-1) * (ndime) + idime
!         elrhs(ievae) = elrhs(ievae) - setfv_derv(3,idime,inode)
!      end do     
! !      ievab = (lboel(inodb)-1) * ndime
! !      do idime = 1,ndime
! !         ievab = ievab + 1
! !         eldcost_dx(ievab) = eldcost_dx(ievab) - setfp_derx(1,idime,inodb)*100
! !      end do
!   end do
  
  
end subroutine nsi_bouass_der

