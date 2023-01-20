!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine elmca2(&
     pnode,pgaus,plapl,weigp,shapf,deriv,heslo,&
     elcod,gpvol,gpsha,gpder,gpcar,gphes,ielem)
  !-----------------------------------------------------------------------
  !****f* domain/elmca2
  ! NAME
  !    elmca2
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
  use def_parame, only       :  twopi
  use def_domain, only       :  ndime,ntens,kfl_naxis,kfl_spher,mnode,&
       &                        kfl_savda,elmda,lelch,ltype
  use def_elmtyp, only       :  ELCUT
  use def_kintyp, only       :  ip,rp
  use mod_cutele, only       :  elmcar_cut
  implicit none
  integer(ip), intent(in)    :: pnode,plapl,ielem
  integer(ip), intent(inout) :: pgaus
  real(rp),    intent(in)    :: weigp(*)
  real(rp),    intent(in)    :: shapf(pnode,*)
  real(rp),    intent(in)    :: deriv(ndime,pnode,*)
  real(rp),    intent(in)    :: heslo(ntens,pnode,*)
  real(rp),    intent(in)    :: elcod(ndime,pnode)
  real(rp),    intent(out)   :: gpvol(*)
  real(rp),    intent(out)   :: gpsha(pnode,*)
  real(rp),    intent(out)   :: gpder(ndime,pnode,*)
  real(rp),    intent(out)   :: gpcar(ndime,mnode,*)
  real(rp),    intent(out)   :: gphes(ntens,mnode,*)
  integer(ip)                :: igaus,inode
  real(rp)                   :: gpcod,gpdet,d2sdx(27)
  real(rp)                   :: xjacm(9),xjaci(9)
  
  if( lelch(ielem) == ELCUT ) then
     !
     ! Cut elements
     !
     call elmcar_cut(&
          ltype(ielem),pnode,plapl,ielem,elcod,pgaus,gpvol,gpsha,gpder,gpcar,gphes)
  else

     gpsha(1:pnode,1:pgaus)         = shapf(1:pnode,1:pgaus)
     gpder(1:ndime,1:pnode,1:pgaus) = deriv(1:ndime,1:pnode,1:pgaus) 
     
     if( kfl_savda == 2 .and. plapl == 1 ) then
        !
        ! Element data base with Hessian
        !
        gpvol(1:pgaus)                 = elmda(ielem) % gpvol(1:pgaus) 
        gpcar(1:ndime,1:pnode,1:pgaus) = elmda(ielem) % gpcar(1:ndime,1:pnode,1:pgaus) 
        gphes(1:ntens,1:pnode,1:pgaus) = elmda(ielem) % gphes(1:ntens,1:pnode,1:pgaus) 

     else if( kfl_savda == 2 .and. plapl == 0 ) then
        !
        ! Element data base without Hessian
        !
        gpvol(1:pgaus)                 = elmda(ielem) % gpvol(1:pgaus) 
        gpcar(1:ndime,1:pnode,1:pgaus) = elmda(ielem) % gpcar(1:ndime,1:pnode,1:pgaus) 
        gphes(1:ntens,1:pnode,1:pgaus) = 0.0_rp

     else

        if( plapl == 0 ) then

           !-------------------------------------------------------------------
           !
           ! GPCAR and GPVOL, Hessian GPHES not needed
           !
           !-------------------------------------------------------------------

           gphes(1:ntens,1:pnode,1:pgaus) = 0.0_rp

           if( ( ndime == 2 .and. pnode == 3 .and. pgaus <= 3 ) .or. ( ndime == 3 .and. pnode == 4 .and. pgaus <= 4 ) ) then
              !
              ! 2D-3D P1 element (linear elements)
              !
              call elmdel(pnode,ndime,elcod,gpcar,gpdet)
              gpvol(1)       = weigp(1) * gpdet
              gpvol(2:pgaus) = gpvol(1)

              do igaus = 2,pgaus
                 gpcar(1:ndime,1:pnode,igaus) = gpcar(1:ndime,1:pnode,1)
              end do

           else
              !
              ! Other elements
              !
              do igaus = 1,pgaus
                 call jacobi(&
                      ndime,pnode,elcod,deriv(1,1,igaus),&
                      xjacm,xjaci,gpcar(1,1,igaus),gpdet)  
                 gpvol(igaus) = weigp(igaus) * gpdet
              end do
           end if

        else

           !-------------------------------------------------------------------
           !
           ! GPCAR, GPHES and GPVOL: Hessian needed
           !
           !-------------------------------------------------------------------

           do igaus = 1,pgaus
              call jacobi(&
                   ndime,pnode,elcod,deriv(1,1,igaus),&
                   xjacm,xjaci,gpcar(1,1,igaus),gpdet)
              gpvol(igaus) = weigp(igaus) * gpdet
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
                 gpcod = gpcod + elcod(1,inode) * shapf(inode,igaus)
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
                 gpcod = gpcod + elcod(1,inode) * elcod(1,inode) * shapf(inode,igaus)
              end do
              gpvol(igaus) = 2.0_rp * twopi * gpcod * gpvol(igaus)
           end do

        end if

     end if

  end if

end subroutine elmca2
