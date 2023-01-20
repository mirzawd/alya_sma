!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine elmhes(&
     heslo,gphes,ndime,pnode,ntens,xjaci,d2sdx,deriv,elcod)
  !-----------------------------------------------------------------------
  !****f* Domain/elmhes
  ! NAME
  !    elmhes
  ! DESCRIPTION
  !    This routine computes the Gphesan matrix in the Cartesian system
  !    of coordinates according to the rule
  !    d^2 N / d x_i d x_j
  !       = (d^2 N / d s_k d s_l) (d s_k/ d x_i)(d s_l/ d x_j)
  !       + (d N / d s_k) (d^2 s_k / d x_i d x_j) 
  ! USES
  !    vetoma
  !    btdbma
  !    matove
  ! USED BY
  !    nsm_elmope
  !    tem_elmope
  !    cdr_elmope
  ! SOURCE
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  implicit none
  integer(ip), intent(in)    :: ndime,pnode,ntens
  real(rp),    intent(in)    :: xjaci(ndime,ndime)
  real(rp),    intent(in)    :: deriv(ndime,pnode),elcod(ndime,pnode)
  real(rp),    intent(in)    :: heslo(ntens,pnode)
  real(rp),    intent(inout) :: gphes(ntens,pnode) 
  real(rp),    intent(out)   :: d2sdx(ndime,ndime,ndime)
  integer(ip)                :: inode
  real(rp)                   :: heslo1, heslo2, heslo3, heslo4, heslo5, heslo6
  real(rp)                   :: xjaci11, xjaci12, xjaci13
  real(rp)                   :: xjaci21, xjaci22, xjaci23
  real(rp)                   :: xjaci31, xjaci32, xjaci33

  if(ndime==1) then
     !
     ! 1D
     !
     do inode=1,pnode
        gphes(1,inode)=xjaci(1,1)*heslo(1,inode)*xjaci(1,1)&
             -deriv(1,inode)*xjaci(1,1)*gphes(1,inode)*elcod(1,inode)
     end do

  else if(ndime==2) then
     !
     ! 2D
     !
     xjaci11 = xjaci(1,1)
     xjaci12 = xjaci(1,2)
     xjaci21 = xjaci(2,1)
     xjaci22 = xjaci(2,2)

     do inode=1,pnode

        gphes(1,inode)   = 0.0_rp
        gphes(2,inode)   = 0.0_rp
        gphes(3,inode)   = 0.0_rp

        heslo1           = heslo(1,inode)
        heslo2           = heslo(2,inode)
        heslo3           = heslo(3,inode)
        
        gphes(1,inode)=gphes(1,inode)+xjaci11*heslo1*xjaci11
        gphes(1,inode)=gphes(1,inode)+xjaci11*heslo3*xjaci21
        gphes(1,inode)=gphes(1,inode)+xjaci21*heslo3*xjaci11
        gphes(1,inode)=gphes(1,inode)+xjaci21*heslo2*xjaci21

        gphes(3,inode)=gphes(3,inode)+xjaci11*heslo1*xjaci12
        gphes(3,inode)=gphes(3,inode)+xjaci11*heslo3*xjaci22
        gphes(3,inode)=gphes(3,inode)+xjaci21*heslo3*xjaci12
        gphes(3,inode)=gphes(3,inode)+xjaci21*heslo2*xjaci22

        gphes(2,inode)=gphes(2,inode)+xjaci12*heslo1*xjaci12
        gphes(2,inode)=gphes(2,inode)+xjaci12*heslo3*xjaci22
        gphes(2,inode)=gphes(2,inode)+xjaci22*heslo3*xjaci12
        gphes(2,inode)=gphes(2,inode)+xjaci22*heslo2*xjaci22

     end do

     d2sdx(1,1,1)=0.0_rp
     d2sdx(2,1,1)=0.0_rp
     d2sdx(1,1,2)=0.0_rp
     d2sdx(2,1,2)=0.0_rp
     d2sdx(1,2,2)=0.0_rp
     d2sdx(2,2,2)=0.0_rp

     do inode=1,pnode
        d2sdx(1,1,1)=d2sdx(1,1,1) - xjaci11*gphes(1,inode)*elcod(1,inode)
        d2sdx(1,1,1)=d2sdx(1,1,1) - xjaci12*gphes(1,inode)*elcod(2,inode)
        d2sdx(2,1,1)=d2sdx(2,1,1) - xjaci21*gphes(1,inode)*elcod(1,inode)
        d2sdx(2,1,1)=d2sdx(2,1,1) - xjaci22*gphes(1,inode)*elcod(2,inode)
        d2sdx(1,1,2)=d2sdx(1,1,2) - xjaci11*gphes(3,inode)*elcod(1,inode)
        d2sdx(1,1,2)=d2sdx(1,1,2) - xjaci12*gphes(3,inode)*elcod(2,inode)
        d2sdx(2,1,2)=d2sdx(2,1,2) - xjaci21*gphes(3,inode)*elcod(1,inode)
        d2sdx(2,1,2)=d2sdx(2,1,2) - xjaci22*gphes(3,inode)*elcod(2,inode)
        d2sdx(1,2,2)=d2sdx(1,2,2) - xjaci11*gphes(2,inode)*elcod(1,inode)
        d2sdx(1,2,2)=d2sdx(1,2,2) - xjaci12*gphes(2,inode)*elcod(2,inode)
        d2sdx(2,2,2)=d2sdx(2,2,2) - xjaci21*gphes(2,inode)*elcod(1,inode)
        d2sdx(2,2,2)=d2sdx(2,2,2) - xjaci22*gphes(2,inode)*elcod(2,inode)
     end do
     !
     ! Computes the second Cartesian derivatives of the shape functions
     !
     do inode=1,pnode

        gphes(1,inode)=gphes(1,inode)+deriv(1,inode)*d2sdx(1,1,1)
        gphes(1,inode)=gphes(1,inode)+deriv(2,inode)*d2sdx(2,1,1)

        gphes(3,inode)=gphes(3,inode)+deriv(1,inode)*d2sdx(1,1,2)
        gphes(3,inode)=gphes(3,inode)+deriv(2,inode)*d2sdx(2,1,2)

        gphes(2,inode)=gphes(2,inode)+deriv(1,inode)*d2sdx(1,2,2)
        gphes(2,inode)=gphes(2,inode)+deriv(2,inode)*d2sdx(2,2,2)

     end do

  else
     !
     ! 3D
     !
     xjaci11 = xjaci(1,1)
     xjaci12 = xjaci(1,2)
     xjaci13 = xjaci(1,3)
     xjaci21 = xjaci(2,1)
     xjaci22 = xjaci(2,2)
     xjaci23 = xjaci(2,3)
     xjaci31 = xjaci(3,1)
     xjaci32 = xjaci(3,2)
     xjaci33 = xjaci(3,3)

     do inode=1,pnode

        gphes(1,inode)   = 0.0_rp
        gphes(4,inode)   = 0.0_rp
        gphes(5,inode)   = 0.0_rp
        gphes(2,inode)   = 0.0_rp
        gphes(6,inode)   = 0.0_rp
        gphes(3,inode)   = 0.0_rp

        heslo1           = heslo(1,inode)
        heslo2           = heslo(2,inode)
        heslo3           = heslo(3,inode)
        heslo4           = heslo(4,inode)
        heslo5           = heslo(5,inode)
        heslo6           = heslo(6,inode)
        
        gphes(1,inode)=gphes(1,inode)+xjaci11*heslo1*xjaci11
        gphes(1,inode)=gphes(1,inode)+xjaci11*heslo4*xjaci21
        gphes(1,inode)=gphes(1,inode)+xjaci11*heslo5*xjaci31
        gphes(1,inode)=gphes(1,inode)+xjaci21*heslo4*xjaci11
        gphes(1,inode)=gphes(1,inode)+xjaci21*heslo2*xjaci21
        gphes(1,inode)=gphes(1,inode)+xjaci21*heslo6*xjaci31
        gphes(1,inode)=gphes(1,inode)+xjaci31*heslo5*xjaci11
        gphes(1,inode)=gphes(1,inode)+xjaci31*heslo6*xjaci21
        gphes(1,inode)=gphes(1,inode)+xjaci31*heslo3*xjaci31
        gphes(4,inode)=gphes(4,inode)+xjaci11*heslo1*xjaci12
        gphes(4,inode)=gphes(4,inode)+xjaci11*heslo4*xjaci22
        gphes(4,inode)=gphes(4,inode)+xjaci11*heslo5*xjaci32
        gphes(4,inode)=gphes(4,inode)+xjaci21*heslo4*xjaci12
        gphes(4,inode)=gphes(4,inode)+xjaci21*heslo2*xjaci22
        gphes(4,inode)=gphes(4,inode)+xjaci21*heslo6*xjaci32
        gphes(4,inode)=gphes(4,inode)+xjaci31*heslo5*xjaci12
        gphes(4,inode)=gphes(4,inode)+xjaci31*heslo6*xjaci22
        gphes(4,inode)=gphes(4,inode)+xjaci31*heslo3*xjaci32
        gphes(5,inode)=gphes(5,inode)+xjaci11*heslo1*xjaci13
        gphes(5,inode)=gphes(5,inode)+xjaci11*heslo4*xjaci23
        gphes(5,inode)=gphes(5,inode)+xjaci11*heslo5*xjaci33
        gphes(5,inode)=gphes(5,inode)+xjaci21*heslo4*xjaci13
        gphes(5,inode)=gphes(5,inode)+xjaci21*heslo2*xjaci23
        gphes(5,inode)=gphes(5,inode)+xjaci21*heslo6*xjaci33
        gphes(5,inode)=gphes(5,inode)+xjaci31*heslo5*xjaci13
        gphes(5,inode)=gphes(5,inode)+xjaci31*heslo6*xjaci23
        gphes(5,inode)=gphes(5,inode)+xjaci31*heslo3*xjaci33

        gphes(2,inode)=gphes(2,inode)+xjaci12*heslo1*xjaci12
        gphes(2,inode)=gphes(2,inode)+xjaci12*heslo4*xjaci22
        gphes(2,inode)=gphes(2,inode)+xjaci12*heslo5*xjaci32
        gphes(2,inode)=gphes(2,inode)+xjaci22*heslo4*xjaci12
        gphes(2,inode)=gphes(2,inode)+xjaci22*heslo2*xjaci22
        gphes(2,inode)=gphes(2,inode)+xjaci22*heslo6*xjaci32
        gphes(2,inode)=gphes(2,inode)+xjaci32*heslo5*xjaci12
        gphes(2,inode)=gphes(2,inode)+xjaci32*heslo6*xjaci22
        gphes(2,inode)=gphes(2,inode)+xjaci32*heslo3*xjaci32
        gphes(6,inode)=gphes(6,inode)+xjaci12*heslo1*xjaci13
        gphes(6,inode)=gphes(6,inode)+xjaci12*heslo4*xjaci23
        gphes(6,inode)=gphes(6,inode)+xjaci12*heslo5*xjaci33
        gphes(6,inode)=gphes(6,inode)+xjaci22*heslo4*xjaci13
        gphes(6,inode)=gphes(6,inode)+xjaci22*heslo2*xjaci23
        gphes(6,inode)=gphes(6,inode)+xjaci22*heslo6*xjaci33
        gphes(6,inode)=gphes(6,inode)+xjaci32*heslo5*xjaci13
        gphes(6,inode)=gphes(6,inode)+xjaci32*heslo6*xjaci23
        gphes(6,inode)=gphes(6,inode)+xjaci32*heslo3*xjaci33

        gphes(3,inode)=gphes(3,inode)+xjaci13*heslo1*xjaci13
        gphes(3,inode)=gphes(3,inode)+xjaci13*heslo4*xjaci23
        gphes(3,inode)=gphes(3,inode)+xjaci13*heslo5*xjaci33
        gphes(3,inode)=gphes(3,inode)+xjaci23*heslo4*xjaci13
        gphes(3,inode)=gphes(3,inode)+xjaci23*heslo2*xjaci23
        gphes(3,inode)=gphes(3,inode)+xjaci23*heslo6*xjaci33
        gphes(3,inode)=gphes(3,inode)+xjaci33*heslo5*xjaci13
        gphes(3,inode)=gphes(3,inode)+xjaci33*heslo6*xjaci23
        gphes(3,inode)=gphes(3,inode)+xjaci33*heslo3*xjaci33

     end do

     d2sdx(1,1,1)=0.0_rp
     d2sdx(2,1,1)=0.0_rp
     d2sdx(3,1,1)=0.0_rp
     d2sdx(1,1,2)=0.0_rp
     d2sdx(2,1,2)=0.0_rp
     d2sdx(3,1,2)=0.0_rp
     d2sdx(1,2,2)=0.0_rp
     d2sdx(2,2,2)=0.0_rp
     d2sdx(3,2,2)=0.0_rp
     d2sdx(1,1,3)=0.0_rp
     d2sdx(2,1,3)=0.0_rp
     d2sdx(3,1,3)=0.0_rp
     d2sdx(1,2,3)=0.0_rp
     d2sdx(2,2,3)=0.0_rp
     d2sdx(3,2,3)=0.0_rp
     d2sdx(1,3,3)=0.0_rp
     d2sdx(2,3,3)=0.0_rp
     d2sdx(3,3,3)=0.0_rp

     do inode=1,pnode
        d2sdx(1,1,1)=d2sdx(1,1,1) - xjaci11*gphes(1,inode)*elcod(1,inode)
        d2sdx(1,1,1)=d2sdx(1,1,1) - xjaci12*gphes(1,inode)*elcod(2,inode)
        d2sdx(1,1,1)=d2sdx(1,1,1) - xjaci13*gphes(1,inode)*elcod(3,inode)
        d2sdx(2,1,1)=d2sdx(2,1,1) - xjaci21*gphes(1,inode)*elcod(1,inode)
        d2sdx(2,1,1)=d2sdx(2,1,1) - xjaci22*gphes(1,inode)*elcod(2,inode)
        d2sdx(2,1,1)=d2sdx(2,1,1) - xjaci23*gphes(1,inode)*elcod(3,inode)
        d2sdx(3,1,1)=d2sdx(3,1,1) - xjaci31*gphes(1,inode)*elcod(1,inode)
        d2sdx(3,1,1)=d2sdx(3,1,1) - xjaci32*gphes(1,inode)*elcod(2,inode)
        d2sdx(3,1,1)=d2sdx(3,1,1) - xjaci33*gphes(1,inode)*elcod(3,inode)
        d2sdx(1,1,2)=d2sdx(1,1,2) - xjaci11*gphes(4,inode)*elcod(1,inode)
        d2sdx(1,1,2)=d2sdx(1,1,2) - xjaci12*gphes(4,inode)*elcod(2,inode)
        d2sdx(1,1,2)=d2sdx(1,1,2) - xjaci13*gphes(4,inode)*elcod(3,inode)
        d2sdx(2,1,2)=d2sdx(2,1,2) - xjaci21*gphes(4,inode)*elcod(1,inode)
        d2sdx(2,1,2)=d2sdx(2,1,2) - xjaci22*gphes(4,inode)*elcod(2,inode)
        d2sdx(2,1,2)=d2sdx(2,1,2) - xjaci23*gphes(4,inode)*elcod(3,inode)
        d2sdx(3,1,2)=d2sdx(3,1,2) - xjaci31*gphes(4,inode)*elcod(1,inode)
        d2sdx(3,1,2)=d2sdx(3,1,2) - xjaci32*gphes(4,inode)*elcod(2,inode)
        d2sdx(3,1,2)=d2sdx(3,1,2) - xjaci33*gphes(4,inode)*elcod(3,inode)
        d2sdx(1,2,2)=d2sdx(1,2,2) - xjaci11*gphes(2,inode)*elcod(1,inode)
        d2sdx(1,2,2)=d2sdx(1,2,2) - xjaci12*gphes(2,inode)*elcod(2,inode)
        d2sdx(1,2,2)=d2sdx(1,2,2) - xjaci13*gphes(2,inode)*elcod(3,inode)
        d2sdx(2,2,2)=d2sdx(2,2,2) - xjaci21*gphes(2,inode)*elcod(1,inode)
        d2sdx(2,2,2)=d2sdx(2,2,2) - xjaci22*gphes(2,inode)*elcod(2,inode)
        d2sdx(2,2,2)=d2sdx(2,2,2) - xjaci23*gphes(2,inode)*elcod(3,inode)
        d2sdx(3,2,2)=d2sdx(3,2,2) - xjaci31*gphes(2,inode)*elcod(1,inode)
        d2sdx(3,2,2)=d2sdx(3,2,2) - xjaci32*gphes(2,inode)*elcod(2,inode)
        d2sdx(3,2,2)=d2sdx(3,2,2) - xjaci33*gphes(2,inode)*elcod(3,inode)
        d2sdx(1,1,3)=d2sdx(1,1,3) - xjaci11*gphes(5,inode)*elcod(1,inode)
        d2sdx(1,1,3)=d2sdx(1,1,3) - xjaci12*gphes(5,inode)*elcod(2,inode)
        d2sdx(1,1,3)=d2sdx(1,1,3) - xjaci13*gphes(5,inode)*elcod(3,inode)
        d2sdx(2,1,3)=d2sdx(2,1,3) - xjaci21*gphes(5,inode)*elcod(1,inode)
        d2sdx(2,1,3)=d2sdx(2,1,3) - xjaci22*gphes(5,inode)*elcod(2,inode)
        d2sdx(2,1,3)=d2sdx(2,1,3) - xjaci23*gphes(5,inode)*elcod(3,inode)
        d2sdx(3,1,3)=d2sdx(3,1,3) - xjaci31*gphes(5,inode)*elcod(1,inode)
        d2sdx(3,1,3)=d2sdx(3,1,3) - xjaci32*gphes(5,inode)*elcod(2,inode)
        d2sdx(3,1,3)=d2sdx(3,1,3) - xjaci33*gphes(5,inode)*elcod(3,inode)
        d2sdx(1,2,3)=d2sdx(1,2,3) - xjaci11*gphes(6,inode)*elcod(1,inode)
        d2sdx(1,2,3)=d2sdx(1,2,3) - xjaci12*gphes(6,inode)*elcod(2,inode)
        d2sdx(1,2,3)=d2sdx(1,2,3) - xjaci13*gphes(6,inode)*elcod(3,inode)
        d2sdx(2,2,3)=d2sdx(2,2,3) - xjaci21*gphes(6,inode)*elcod(1,inode)
        d2sdx(2,2,3)=d2sdx(2,2,3) - xjaci22*gphes(6,inode)*elcod(2,inode)
        d2sdx(2,2,3)=d2sdx(2,2,3) - xjaci23*gphes(6,inode)*elcod(3,inode)
        d2sdx(3,2,3)=d2sdx(3,2,3) - xjaci31*gphes(6,inode)*elcod(1,inode)
        d2sdx(3,2,3)=d2sdx(3,2,3) - xjaci32*gphes(6,inode)*elcod(2,inode)
        d2sdx(3,2,3)=d2sdx(3,2,3) - xjaci33*gphes(6,inode)*elcod(3,inode)
        d2sdx(1,3,3)=d2sdx(1,3,3) - xjaci11*gphes(3,inode)*elcod(1,inode)
        d2sdx(1,3,3)=d2sdx(1,3,3) - xjaci12*gphes(3,inode)*elcod(2,inode)
        d2sdx(1,3,3)=d2sdx(1,3,3) - xjaci13*gphes(3,inode)*elcod(3,inode)
        d2sdx(2,3,3)=d2sdx(2,3,3) - xjaci21*gphes(3,inode)*elcod(1,inode)
        d2sdx(2,3,3)=d2sdx(2,3,3) - xjaci22*gphes(3,inode)*elcod(2,inode)
        d2sdx(2,3,3)=d2sdx(2,3,3) - xjaci23*gphes(3,inode)*elcod(3,inode)
        d2sdx(3,3,3)=d2sdx(3,3,3) - xjaci31*gphes(3,inode)*elcod(1,inode)
        d2sdx(3,3,3)=d2sdx(3,3,3) - xjaci32*gphes(3,inode)*elcod(2,inode)
        d2sdx(3,3,3)=d2sdx(3,3,3) - xjaci33*gphes(3,inode)*elcod(3,inode)
     end do
     !
     ! Computes the second Cartesian derivatives of the shape functions
     !
     do inode=1,pnode

        gphes(1,inode)=gphes(1,inode)+deriv(1,inode)*d2sdx(1,1,1)
        gphes(1,inode)=gphes(1,inode)+deriv(2,inode)*d2sdx(2,1,1)
        gphes(1,inode)=gphes(1,inode)+deriv(3,inode)*d2sdx(3,1,1)
        gphes(4,inode)=gphes(4,inode)+deriv(1,inode)*d2sdx(1,1,2)
        gphes(4,inode)=gphes(4,inode)+deriv(2,inode)*d2sdx(2,1,2)
        gphes(4,inode)=gphes(4,inode)+deriv(3,inode)*d2sdx(3,1,2)
        gphes(5,inode)=gphes(5,inode)+deriv(1,inode)*d2sdx(1,1,3)
        gphes(5,inode)=gphes(5,inode)+deriv(2,inode)*d2sdx(2,1,3)
        gphes(5,inode)=gphes(5,inode)+deriv(3,inode)*d2sdx(3,1,3)

        gphes(2,inode)=gphes(2,inode)+deriv(1,inode)*d2sdx(1,2,2)
        gphes(2,inode)=gphes(2,inode)+deriv(2,inode)*d2sdx(2,2,2)
        gphes(2,inode)=gphes(2,inode)+deriv(3,inode)*d2sdx(3,2,2)
        gphes(6,inode)=gphes(6,inode)+deriv(1,inode)*d2sdx(1,2,3)
        gphes(6,inode)=gphes(6,inode)+deriv(2,inode)*d2sdx(2,2,3)
        gphes(6,inode)=gphes(6,inode)+deriv(3,inode)*d2sdx(3,2,3)

        gphes(3,inode)=gphes(3,inode)+deriv(1,inode)*d2sdx(1,3,3)
        gphes(3,inode)=gphes(3,inode)+deriv(2,inode)*d2sdx(2,3,3)
        gphes(3,inode)=gphes(3,inode)+deriv(3,inode)*d2sdx(3,3,3)

     end do

  end if

end subroutine elmhes
