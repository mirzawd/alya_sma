!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine extens(ndofn,kfl_fixno_tmp,xdire,unkno_tmp)
  !-----------------------------------------------------------------------
  !****f* mathru/extens
  ! NAME 
  !    extens
  ! DESCRIPTION
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_elmtyp
  use def_domain
  implicit none
  integer(ip), intent(in)    :: ndofn
  integer(ip), intent(in)    :: kfl_fixno_tmp(ndofn,npoin)
  real(rp),    intent(in)    :: xdire(ndime) 
  real(rp),    intent(inout) :: unkno_tmp(ndofn,*)
  integer(ip)                :: kfl_fixel(ndofn,mnode)
  integer(ip)                :: ielem,pnode,pgaus,idime,inode,ipoin
  integer(ip)                :: pelty,pmate,plapl,porde,ptopo,idofn
  real(rp)                   :: elmat(mnode,mnode)
  real(rp)                   :: elrhs(mnode)
  real(rp)                   :: elcod(ndime,mnode)
  real(rp)                   :: gpcar(ndime,mnode,mgaus)  ! dN/dxi
  real(rp)                   :: gphes(ntens,mnode,mgaus)  ! d2N/dxidxj
  real(rp)                   :: gpvol(mgaus)              ! w*|J|, |J|
  real(rp)                   :: tragl(9)                  ! Stabilization
  real(rp)                   :: hleng(3)                  ! Element length
  real(rp)                   :: bvele(ndofn,mnode)

  if( INOTMASTER ) then

     elements: do ielem = 1,nelem

        if( lelch(ielem) /= ELHOL ) then
           !
           ! Element properties and dimensions
           !
           pelty = ltype(ielem)
           pnode = nnode(pelty)
           pgaus = ngaus(pelty)
           plapl = 0
           porde = lorde(pelty)
           ptopo = ltopo(pelty)
           pmate = 1
           !
           ! Check if element is a solid
           !
           if( pmate /= -1 ) then
              !
              ! Gather operations: ELVEL, ELCOD
              !
              do inode = 1,pnode
                 ipoin = lnods(inode,ielem)
                 do idofn = 1,ndofn
                    kfl_fixel(idofn,inode) = kfl_fixno_tmp(idofn,ipoin)
                    bvele(idofn,inode)     = unkno_tmp(idofn,ipoin)
                 end do
                 do idime = 1,ndime
                    elcod(idime,inode) = coord(idime,ipoin)
                 end do
              end do
              !
              ! Cartesian derivatives, Hessian matrix and volume: GPCAR, PGVOL
              !
              call elmcar(&
                   pnode,pgaus,plapl,elmar(pelty) % weigp,elmar(pelty) % shape,&
                   elmar(pelty) % deriv,elmar(pelty) % heslo,elcod,gpvol,gpcar,&
                   gphes,ielem)
              !
              ! HLENG and TRAGL at center of gravity 
              !
              call elmlen(&
                   ndime,pnode,elmar(pelty) % dercg,tragl,elcod,&
                   hnatu(pelty),hleng)
              !
              ! Assemble elemental matrix
              !
              call elexte(&
                   ndofn,pgaus,pnode,lnods(:,ielem),kfl_fixel,gpvol,gpcar,&
                   elmar(pelty) % shape,hleng,bvele,xdire,elmat,&
                   elrhs,ielem)
              !
              ! Assembly: AMATR and RHSID
              !
              call assrhs(&
                   1_ip,pnode,lnods(1,ielem),elrhs,rhsid)
              call assmat(&
                   solve_sol(1) % ndofn,pnode,pnode,solve_sol(1) % nunkn,&
                   solve_sol(1) % kfl_algso,ielem,lnods(:,ielem),elmat,amatr)

           end if

        end if

     end do elements

  end if
  !
  ! Solve system
  !
  call solver(rhsid,unkno_tmp,amatr,pmatr) 
   print*,'extens.f90'
end subroutine extens

subroutine elexte(&
     ndofn,pgaus,pnode,lnods,kfl_fixel,gpvol,gpcar,gpsha,&
     hleng,bvele,xdire,elmat,elrhs,ielem)

  !-----------------------------------------------------------------------
  !****f* mathru/elexte
  ! NAME 
  !    elexte
  ! DESCRIPTION
  !    Assemble the elemental matrix from Gauss point contributions
  !    Solve the following advection problem with source term:
  !
  !         grad(p) = rho g    =>
  !    g  . grad(p) = rho g^2  =>
  !    g' . grad(p) = rho g
  !
  !    where g' is the unit gravity vector.
  !
  !    ( rho g - g'.grad(p) , v + tau*g'.grad(v) ) = 0
  !   
  ! USES
  ! USED BY
  !    nsi_elmhyd
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_elmtyp, only     :  ELEXT
  use def_domain, only     :  mnode,ndime,lelch
  use def_master, only     :  zeror
  implicit none
  integer(ip), intent(in)  :: ndofn
  integer(ip), intent(in)  :: pgaus,pnode,ielem
  integer(ip), intent(in)  :: lnods(pnode)
  integer(ip), intent(in)  :: kfl_fixel(ndofn,pnode)
  real(rp),    intent(in)  :: gpvol(pgaus)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)  :: gpsha(pnode,pgaus)
  real(rp),    intent(in)  :: hleng(ndime)
  real(rp),    intent(in)  :: bvele(ndofn,pnode)
  real(rp),    intent(in)  :: xdire(ndime)
  real(rp),    intent(out) :: elmat(pnode,pnode)
  real(rp),    intent(out) :: elrhs(pnode)
  integer(ip)              :: igaus,inode,jnode,ipoin,idime
  real(rp)                 :: fact1,tau,ggrap(mnode)
  real(rp)                 :: vtest,dummr,u,h

  !----------------------------------------------------------------------
  !
  ! Compute LHS and RHS
  !
  !----------------------------------------------------------------------

  do inode = 1,pnode  
     elrhs(inode) = 0.0_rp
     do jnode = 1,pnode
        elmat(jnode,inode) = 0.0_rp
     end do
  end do

  u = 0.0_rp
  do idime = 1,ndime
     u = u + xdire(idime) * xdire(idime)
  end do
  h   = hleng(ndime)
  tau = h/(2.0_rp*sqrt(u))

  do igaus = 1,pgaus
     !
     ! a . grad(Ni)
     !
     do inode = 1,pnode
        ggrap(inode)  = 0.0_rp
        do idime = 1,ndime
           ggrap(inode) = ggrap(inode) + xdire(idime) * gpcar(idime,inode,igaus)
        end do
     end do

     do inode = 1,pnode 
        vtest = ( gpsha(inode,igaus) + tau * ggrap(inode) ) * gpvol(igaus)
        do jnode = 1,pnode
           elmat(inode,jnode) = elmat(inode,jnode) + vtest * ggrap(jnode)
        end do
     end do

  end do

  !----------------------------------------------------------------------
  !
  ! Extension elements
  !
  !----------------------------------------------------------------------

  if( lelch(ielem) == ELEXT ) then
     call elmext(&
          4_ip,1_ip,pnode,dummr,dummr,dummr,dummr,elmat,&
          dummr,elrhs,dummr)
  end if

  !----------------------------------------------------------------------
  !
  ! Prescribe boundary conditions at z-plane z=hydro_nsi
  !
  !----------------------------------------------------------------------

  do inode = 1,pnode
     ipoin = lnods(inode)
     if( kfl_fixel(1,inode) /= 0 ) then
        fact1 = max(zeror,elmat(inode,inode))
        do jnode = 1,pnode
           elrhs(jnode)       = elrhs(jnode) - elmat(jnode,inode) * bvele(1,inode)
           elmat(jnode,inode) = 0.0_rp
           elmat(inode,jnode) = 0.0_rp              
        end do
        elmat(inode,inode) = fact1
        elrhs(inode)       = fact1 * bvele(1,inode)
     end if
  end do

end subroutine elexte
