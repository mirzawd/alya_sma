!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine gpcabo(&
     pnode,pgaus,pnodb,lboel,shaga,gpcar,gbsha,gbcar)
  !-----------------------------------------------------------------------
  !****f* Domain/gpcabo
  ! NAME 
  !    gpcabo
  ! DESCRIPTION
  !    This routine computes the cartesian derivates on the 
  !    boundary Gauss points 
  ! USES
  ! USED BY
  !    nsi_bouset
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime,mnode
  implicit none
  integer(ip), intent(in)  :: pnode,pgaus,pnodb 
  integer(ip), intent(in)  :: lboel(pnodb)
  real(rp),    intent(in)  :: shaga(pgaus,pnode)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)  :: gbsha(pnodb)
  real(rp),    intent(out) :: gbcar(ndime,pnode)
  integer(ip)              :: inode,idime,inodb,igaus

  do inode=1,pnode                            
     do idime=1,ndime                          
        gbcar(idime,inode)=0.0_rp                
        do inodb=1,pnodb                         
           do igaus=1,pgaus
              gbcar(idime,inode)=gbcar(idime,inode)&
                   +shaga(igaus,lboel(inodb))&
                   *gbsha(inodb)*gpcar(idime,inode,igaus)
           end do
        end do
     end do
  end do

end subroutine gpcabo
