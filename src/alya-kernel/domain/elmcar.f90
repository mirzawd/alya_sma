!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine elmcar(&
     pnode,pgaus,plapl,weigp,gpsha,deriv,heslo,&
     elcod,gpvol,gpcar,gphes,ielem)
  !-----------------------------------------------------------------------
  !****f* domain/elmcar
  ! NAME
  !    elmcar
  ! DESCRIPTION
  !    This routine calculates:
  !    GPCAR: Cartesian derivatives
  !    GPVOL: Unit volume
  !    GPHES: Hessian matrix
  ! USES
  !    invmtx
  ! USED BY
  !    ***_elmope
  !    extnor
  ! SOURCE
  !-----------------------------------------------------------------------
  use def_parame, only     :  twopi
  use def_domain, only     :  ndime,ntens,kfl_naxis,kfl_spher,mnode,&
       &                      kfl_savda,elmda
  use def_kintyp, only     :  ip,rp
  implicit none
  integer(ip), intent(in)  :: pnode,pgaus,plapl,ielem
  real(rp),    intent(in)  :: weigp(pgaus)
  real(rp),    intent(in)  :: gpsha(pnode,pgaus)
  real(rp),    intent(in)  :: deriv(ndime,pnode,pgaus)
  real(rp),    intent(in)  :: heslo(ntens,pnode,pgaus)
  real(rp),    intent(in)  :: elcod(ndime,pnode)
  real(rp),    intent(out) :: gpvol(pgaus)
  real(rp),    intent(out) :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(out) :: gphes(ntens,mnode,pgaus)
  integer(ip)              :: igaus,inode,idime,itens
  real(rp)                 :: gpcod,gpdet,d2sdx(27)
  real(rp)                 :: xjacm(9),xjaci(9)

  if( kfl_savda == 2 .and. plapl == 1 ) then

     do igaus = 1,pgaus
        gpvol(igaus) = elmda(ielem)%gpvol(igaus) 
        do inode = 1,pnode
           do itens = 1,ntens
              gphes(itens,inode,igaus) = elmda(ielem)%gphes(itens,inode,igaus) 
           end do
           do idime = 1,ndime
              gpcar(idime,inode,igaus) = elmda(ielem)%gpcar(idime,inode,igaus) 
           end do
        end do
     end do

  else if( kfl_savda == 2 .and. plapl == 0 ) then

     do igaus = 1,pgaus
        gpvol(igaus) = elmda(ielem)%gpvol(igaus) 
        do inode = 1,pnode
           do idime = 1,ndime
              gpcar(idime,inode,igaus) = elmda(ielem)%gpcar(idime,inode,igaus) 
           end do
        end do
     end do

  else

     if( plapl == 0 ) then

        !-------------------------------------------------------------------
        !
        ! GPCAR and GPVOL: Hessian not needed
        !
        !-------------------------------------------------------------------

        if( ( ndime == 2 .and. pnode == 3 .and. pgaus <= 3 ) .or. ( ndime == 3 .and. pnode == 4 .and. pgaus <= 4 ) ) then
           !
           ! 2D-3D P1 element (linear elements)
           !
           call elmdel(pnode,ndime,elcod,gpcar,gpdet)
           gpvol(1)=weigp(1)*gpdet
           do igaus=2,pgaus
              do inode=1,pnode
                 do idime=1,ndime
                    gpcar(idime,inode,igaus)=gpcar(idime,inode,1)
                 end do
              end do
              gpvol(igaus)=gpvol(1)
           end do

        else
           !
           ! Other elements
           !
           do igaus=1,pgaus
              call jacobi(&
                   ndime,pnode,elcod,deriv(1,1,igaus),&
                   xjacm,xjaci,gpcar(1,1,igaus),gpdet)  
              gpvol(igaus)=weigp(igaus)*gpdet
           end do
        end if

     else

        !-------------------------------------------------------------------
        !
        ! GPCAR, GPHES and GPVOL: Hessian needed
        !
        !-------------------------------------------------------------------

        do igaus=1,pgaus
           call jacobi(&
                ndime,pnode,elcod,deriv(1,1,igaus),&
                xjacm,xjaci,gpcar(1,1,igaus),gpdet)
           gpvol(igaus)=weigp(igaus)*gpdet
           call elmhes(&
                heslo(1,1,igaus),gphes(1,1,igaus),ndime,pnode,ntens,&
                xjaci,d2sdx,deriv(1,1,igaus),elcod)
        end do
     end if

     if( kfl_naxis == 1 ) then
        !
        ! Axi-symmetric coordinates
        !
        do igaus = 1,pgaus
           gpcod = 0.0_rp
           do inode = 1,pnode
              gpcod = gpcod + elcod(1,inode) * gpsha(inode,igaus)
           end do
           gpvol(igaus) = twopi * gpcod * gpvol(igaus)
        end do

     else if( kfl_spher == 1 ) then
        !
        ! Spherical coordinates
        ! 
        do igaus = 1,pgaus
           gpcod = 0.0_rp
           do inode = 1,pnode
              gpcod = gpcod + elcod(1,inode) * elcod(1,inode) * gpsha(inode,igaus)
           end do
           gpvol(igaus) = 2.0_rp * twopi * gpcod * gpvol(igaus)
        end do

     end if

  end if

end subroutine elmcar
