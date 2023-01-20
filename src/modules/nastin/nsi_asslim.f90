!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_asslim(&
     pnode,pgaus,lnods,gpden,gpsp1,gpadv,gpvep,gpvol,&
     elvel,gpsha,gpcar,agrau,elrhs,rhsid)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_asslim
  ! NAME 
  !    nsi_asslim
  ! DESCRIPTION
  !    Compute Limiter
  ! USES
  ! USED BY
  !    nsi_elmop3
  !***
  !----------------------------------------------------------------------- 
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime,mnode
  use def_nastin, only     :  kfl_stabi_nsi,kfl_limit_nsi
  implicit none
  integer(ip), intent(in)  :: pnode,pgaus
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(in)  :: gpden(pgaus)
  real(rp),    intent(in)  :: gpsp1(pgaus)
  real(rp),    intent(in)  :: gpadv(ndime,pgaus)
  real(rp),    intent(in)  :: gpvep(ndime,pgaus)
  real(rp),    intent(in)  :: gpvol(pgaus)
  real(rp),    intent(in)  :: elvel(ndime,pnode)
  real(rp),    intent(in)  :: gpsha(pnode,pgaus)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(out) :: agrau(pnode,pgaus)
  real(rp),    intent(out) :: elrhs(pnode)
  real(rp),    intent(out) :: rhsid(*)
  integer(ip)              :: igaus,inode,idime
  real(rp)                 :: c1,c2,c3,c4,alpha,beta

  if( kfl_stabi_nsi == 2 .and. kfl_limit_nsi /= 0 ) then

     do inode = 1,pnode
        elrhs(inode) = 0.0_rp
     end do

     if( kfl_limit_nsi > 0 ) then

        if( ndime == 2 ) then
           do igaus = 1,pgaus
              do inode = 1,pnode
                 agrau(inode,igaus) =  gpden(igaus) * (                    &
                      &                gpadv(1,igaus)*gpcar(1,inode,igaus) &
                      &              + gpadv(2,igaus)*gpcar(2,inode,igaus) )
              end do
           end do
        else
           do igaus = 1,pgaus
              do inode = 1,pnode
                 agrau(inode,igaus) =  gpden(igaus) * (                    &
                      &                gpadv(1,igaus)*gpcar(1,inode,igaus) &
                      &              + gpadv(2,igaus)*gpcar(2,inode,igaus) &
                      &              + gpadv(3,igaus)*gpcar(3,inode,igaus) )
              end do
           end do
        end if

        do igaus = 1,pgaus
           c1 = 0.0_rp
           c2 = 0.0_rp
           c3 = 0.0_rp
           do idime = 1,ndime
              c4 = 0.0_rp
              do inode = 1,pnode
                 c4 = c4 + agrau(inode,igaus) * elvel(idime,inode)
              end do
              c4 = gpsp1(igaus) * c4
              c1 = c1 + ( gpvep(idime,igaus) - c4 )**2
              c3 = c3 + gpvep(idime,igaus) * gpvep(idime,igaus)
              c2 = c2 + c4 * c4
           end do
           c3 = sqrt( c2 ) + sqrt( c3 )
           c1 = sqrt( c1 )
           if( c3 /= 0.0_rp ) then
              beta  = c1 / c3
           else
              beta  = 0.0_rp
           end if
           if( kfl_limit_nsi == 1 ) then
              alpha = min(1.0_rp,2.0_rp*(1.0_rp-beta))
           else if( kfl_limit_nsi == 2 ) then
              alpha = 0.5_rp*(tanh(20.0_rp*(beta-0.8_rp))+1.0_rp)
           end if
           do inode = 1,pnode
              elrhs(inode) = elrhs(inode) + alpha * gpvol(igaus) * gpsha(inode,igaus)
           end do
        end do

     end if

     call assrhs(1_ip,pnode,lnods,elrhs,rhsid)

  end if

end subroutine nsi_asslim
