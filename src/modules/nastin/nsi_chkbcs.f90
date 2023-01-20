!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_chkbcs()
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_chkbcs
  ! NAME 
  !    nsi_chkbcs
  ! DESCRIPTION
  !    This routine check boundary conditions
  ! USED BY
  !    nsi_iniunk
  !***
  !-----------------------------------------------------------------------
  use def_parame 
  use def_master
  use def_domain
  use def_nastin
  use mod_local_basis, only : local_basis_global_to_local
  use mod_local_basis, only : local_basis_local_to_global
  implicit none
  real(rp)     :: elmat(nevat_nsi,nevat_nsi)
  real(rp)     :: elrhs(nevat_nsi)
  real(rp)     :: elcod(ndime,mnode)
  real(rp)     :: xjaci(9),xjacm(9) 
  integer(ip)  :: ielem,igaus,idime,icomp,inode,ipoin
  integer(ip)  :: pnode,pgaus,pevat,pelty
  real(rp)     :: gpcar(ndime,mnode,mgaus)
  real(rp)     :: gpvol(mgaus),gpdet

  if( kfl_rstar == 0 ) then
     !
     ! Initialization
     !
     call nsi_inisol(4_ip)

     if( INOTMASTER ) then

        do ipoin=1,npoin
           icomp=(ipoin-1)*ndime
           do idime=1,ndime
              icomp=icomp+1
              unkno(icomp)=veloc(idime,ipoin,nprev_nsi)
           end do
        end do
        do icomp=1,ndime*npoin
           rhsid(icomp)=0.0_rp
        end do
        do icomp=1,solve(3)%nzmat
           amatr(icomp)=0.0_rp
        end do

        elements: do ielem = 1,nelem
           pelty = ltype(ielem)

           if( pelty > 0 ) then
              !
              ! Element properties and dimensions
              !
              pnode=nnode(pelty)
              pgaus=ngaus(pelty)
              pevat=ndime*pnode
              !
              !
              ! 1st and 2nd order Cartesian derivatives GPCAR, GPHES, and GPVOL=dV=|J|*wg
              !
              do inode=1,pnode
                 ipoin=lnods(inode,ielem)
                 do idime=1,ndime
                    elcod(idime,inode)=coord(idime,ipoin)
                 end do
              end do
              do igaus=1,pgaus     
                 call elmder(&
                      pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&        ! Cartesian derivative
                      elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)          ! and Jacobian
                 gpvol(igaus)=elmar(pelty)%weigp(igaus)*gpdet            ! |J|*wg
              end do
              !
              ! Compute Laplacian
              !
              call nsi_elmlbc(&
                   pnode,pgaus,pevat,elmar(pelty)%shape,gpcar,&
                   gpvol,elmat,elrhs)
              !
              ! Prescribe Dirichlet boundary conditions
              !
              call nsi_elmdir(&
                   1_ip,2_ip,pnode,pevat,ndime,lnods(1,ielem),&
                   elmat,elrhs)
              !
              ! Assemble matrix and RHS
              !
              call assrhs(&
                   ndime,pnode,lnods(1,ielem),elrhs,rhsid)
              call assmat(&
                   ndime,pnode,pevat,ndime*npoin,&
                   solve(3)%kfl_algso,ielem,lnods(1,ielem),elmat,amatr)

           end if

        end do elements

     end if
     !
     ! Solve the algebraic system
     !
     call local_basis_global_to_local(kfl_fixrs_nsi,unkno)           ! Global to local     
     call solver(rhsid,unkno,amatr,pmatr)
     call local_basis_local_to_global(kfl_fixrs_nsi,unkno)           ! Global to local     
     !
     ! Postprocess bcs
     !
     do ipoin=1,npoin
        icomp=(ipoin-1)*ndime
        do idime=1,ndime
           icomp=icomp+1
           veloc(idime,ipoin,nprev_nsi)=unkno(icomp)
        end do
     end do

  end if

104 format(/,&
       '#    LAPLACE EQ, BOUNDARY CONDITIONS',/,&
       '#    -------------------------------')

end subroutine nsi_chkbcs
