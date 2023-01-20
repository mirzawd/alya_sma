!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine elmca2_der(&
     pnode,pgaus,plapl,weigp,deriv,elcod,gpvol_der,gpcar_der,ielem)
  !-----------------------------------------------------------------------
  !****f* domain/elmca2_der
  ! NAME
  !    elmca2_der
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
  use def_domain, only       :  ndime,kfl_naxis,kfl_spher,kfl_savda
  use def_elmtyp, only       :  ELCUT
  use def_kintyp, only       :  ip,rp
  use mod_cutele, only       :  elmcar_cut
  use def_kermod, only       :  kfl_adj_prob
  implicit none
  integer(ip), intent(in)    :: pnode,plapl,ielem
  integer(ip), intent(inout) :: pgaus
  real(rp),    intent(in)    :: weigp(*)
  real(rp),    intent(in)    :: deriv(ndime,pnode,*)
  real(rp),    intent(in)    :: elcod(ndime,pnode)
  real(rp),    intent(out)   :: gpvol_der(ndime,pnode,*)
  real(rp),    intent(out)   :: gpcar_der(ndime,pnode,ndime,pnode,*)
  integer(ip)                :: igaus,inode,idime,jdime,jnode
  real(rp)                   :: gpdet_der(ndime,pnode)

  if( kfl_adj_prob == 1 ) then
     
     if( kfl_savda == 2 .and. plapl == 1 ) then

     else if( kfl_savda == 2 .and. plapl == 0 ) then

     else

        if( plapl == 0 ) then

           !-------------------------------------------------------------------
           !
           ! GPVOL_DER: Hessian not needed
           !
           !-------------------------------------------------------------------


           if( ( ndime == 2 .and. pnode == 3 ) .or. ( ndime == 3 .and. pnode == 4 ) ) then
              !
              ! 2D-3D P1 element (linear elements)
              !
              call elmdel_der(pnode,ndime,elcod,gpdet_der,gpcar_der)
              do inode = 1,pnode
                do idime = 1,ndime
                  gpvol_der(idime,inode,1) = weigp(1) * gpdet_der(idime,inode)
                enddo
              enddo
              do igaus = 2,pgaus
                 do inode = 1,pnode
                   do idime = 1,ndime
                     gpvol_der(idime,inode,igaus) = gpvol_der(idime,inode,1)
                     do jnode = 1,pnode
                       do jdime = 1,ndime
                         gpcar_der(idime,inode,jdime,jnode,igaus) = gpcar_der(idime,inode,jdime,jnode,1)
                       end do
                     end do
                     
                   enddo
                 enddo
              end do

           else
              !
              ! Other elements
              !
              do igaus = 1,pgaus
                 call jacobi_der(ndime,pnode,elcod,deriv(1,1,igaus),gpdet_der,gpcar_der(1,1,1,1,igaus))
                 do inode = 1,pnode
                    do idime = 1,ndime
                      gpvol_der(idime,inode,igaus) = weigp(igaus) * gpdet_der(idime,inode)
                    enddo
                 enddo
                 
              end do
           end if

        else

           !-------------------------------------------------------------------
           !
           ! GPCAR, GPHES and GPVOL: Hessian needed
           !
           !-------------------------------------------------------------------
           do igaus = 1,pgaus
              call jacobi_der(ndime,pnode,elcod,deriv(1,1,igaus),gpdet_der,gpcar_der(1,1,1,1,igaus))
              do inode = 1,pnode
                 do idime = 1,ndime
                   gpvol_der(idime,inode,igaus) = weigp(igaus) * gpdet_der(idime,inode)
                 enddo
              enddo              
           end do

        end if

        if( kfl_naxis == 1 ) then
           !
           ! Axi-symmetric coordinates
           !
           call runend('No esta implementado') 

        else if( kfl_spher == 1 ) then
           !
           ! Spherical coordinates
           !
           call runend('No esta implementado') 

        end if

     end if

  end if

end subroutine elmca2_der
