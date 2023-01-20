!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!==============================================================================!
!!
!!
!! Dont put all your eggs in one basket !!
!!
!
!  Patm  = 103400.0            !> [Nm2/s2/m3]
!  Eint  = 5.0/2.0*287.0*300.0 !> [J/kg] = [J/kg/K][K], Eint = 215250.0 (gas ideal)
!
!==============================================================================!
module mod_lodi_nsi
  use def_kermod, only: kfl_prope
  use def_parame, only: ip, rp  
  use def_domain, only: ltype, lnods, coord, hnatu
  use def_domain, only: nnode, ngaus, mnode, mgaus, nelem, ndime, npoin
  use def_master, only: inotmaster !, veloc, umome, densi, energ, densi, press
  use def_domain, only: vmasc, vmass, elmar
  implicit none
  !
  integer(ip)           :: ielem, inode, pelty, pgaus
  integer(ip)           :: pnode, igaus, idofn, idime, ipoin
  !
  integer(ip)           :: ndofn
!  real(rp), pointer     :: my_gamma_nsi(:)        !npoin
  !
  type LODI
    integer(ip)           :: idofn           =  1 
    integer(ip)           :: ndofn           =  0 
    real(rp)              :: gamme           =  1.2_rp
    real(rp), pointer     :: xdvol(:,:)      => null()
    real(rp), pointer     :: xchrc(:,:,:,:)  => null()
    real(rp), pointer     ::  chrc(:,:,:)    => null()
    !
    real(rp), pointer     ::  densi(:)      => null()
    real(rp), pointer     ::  veloc(:,:)    => null()
    real(rp), pointer     ::  press(:)      => null()
  end type 
  type(LODI),save :: CHARAs
  !
  private
  
  public :: CHARAs
  public :: lodi_nsi_allocate
  public :: lodi_nsi_deallocate
  public :: lodi_nsi_get_characteristics

contains 


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine lodi_nsi_allocate( STRUCT ) !, dof, den, vel, pre )
    use def_master, only: veloc, press, tempe
!    use def_master, only: densi
    use def_kermod, only: gasco       ! Gas constant R [J/K Kg] 
!    use def_nastin, only: prthe_nsi   ! Thermodynamic pressure = 101325.0_rp
    use def_master, only:  prthe 
    implicit none 

    type(LODI),  intent(inout) :: STRUCT
    !integer(ip), intent(in   ) :: dof
    !real(rp), pointer     ::  den(:)  
    !real(rp), pointer     ::  vel(:,:)
    !real(rp), pointer     ::  pre(:)     
    !---------------------------------------------------------------------||---!
    if(INOTMASTER) then
      ndofn        = ndime+2
      STRUCT%ndofn = ndofn
!      STRUCT%densi => densi(         1:npoin, 1)
      allocate( STRUCT%densi(npoin) )
      STRUCT%densi(1:npoin) = prthe(1)/(gasco*tempe(1:npoin,1))
      STRUCT%veloc => veloc(1:ndime, 1:npoin, 1)
      STRUCT%press => press(         1:npoin, 1)
      !
      allocate( STRUCT%xchrc(ndofn, nelem, mgaus, ndime) )
      !call memchk(zero,istat,mem_modul(1:2,modul),'LODI_NSA','nsi_membcs', xchrc_nsi)
      !
      allocate( STRUCT%chrc(ndofn, npoin, ndime) )
      !call memchk(zero,istat,mem_modul(1:2,modul),'LODI_NSA','nsi_membcs', chrc_nsi)
      !
      allocate( STRUCT%xdvol(nelem, mgaus) )

      STRUCT%xchrc = 0.0_rp
      STRUCT% chrc = 0.0_rp
      STRUCT%xdvol = 0.0_rp
      !STRUCT%idofn = 1
!
!print *, shape(STRUCT%densi), shape(STRUCT%veloc), shape(STRUCT%press) 
!print *, STRUCT%idofn
!
    endif 
    !---------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!
  
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine lodi_nsi_deallocate( STRUCT )
    implicit none 
    type(LODI),  intent(inout) :: STRUCT
    !
    if(INOTMASTER) then
      deallocate( STRUCT%xchrc )
      deallocate( STRUCT% chrc )
      deallocate( STRUCT%xdvol )
      deallocate( STRUCT%densi )
    endif 
    !---------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine lodi_nsi_get_characteristics( STRUCT )
    use def_master, only:  prthe
  implicit none 
  type(LODI), intent(inout) :: STRUCT
  !-----------------------------------------------------------------------||---!
  !
  !real(rp), intent(inout) :: xdvol(nelem, mgaus)
  !real(rp), intent(inout) :: xchrc(ndofn, nelem, mgaus, ndime)
  !real(rp), intent(inout) ::  chrc(ndofn, npoin, ndime)
  !
  real(rp)         :: hmini
  real(rp)         :: elcod(ndime, mnode)
  real(rp)         :: detjm
  real(rp)         :: xjaci(ndime, ndime), xjacm(ndime, ndime), cartd(ndime, mnode, mgaus)
  real(rp)         :: tragl(ndime, ndime), hleng(ndime) 
  !
  real(rp)    :: xdphit(ndofn, ndime, mgaus)  
  real(rp)    ::  phit(ndofn, mnode), xphit(ndofn, mgaus)
  real(rp)    ::  xvelmo(mgaus), xsound(mgaus) 
  real(rp)    :: hinv
  real(rp)    :: Sxinv(1:ndofn,1:ndofn)
  real(rp)    :: elgamme(mnode)
  real(rp)    ::  xgamme(mgaus)
  !-----------------------------------------------------------------------||---!

  if(INOTMASTER) then
  !-----------------------------------------------------------------------||---!
  elementary_loop: & 
  do ielem = 1,nelem 
    pelty = ltype(ielem)
    pnode = nnode(pelty)
    pgaus = ngaus(pelty)
    !---------------------------------------------------------------------||---!
    do inode = 1,pnode
      ipoin = lnods(inode, ielem) 
      !
      elcod(1:ndime,inode) = coord(1:ndime,ipoin)
      elgamme(inode)       = STRUCT%gamme !my_gamma_nsi(ipoin) 
      !
      ! PHIt -> [rho, vel, p] 
      phit(1,         inode) = STRUCT%densi(         ipoin)  
      phit(2:ndime+1, inode) = STRUCT%veloc(1:ndime, ipoin)
!      phit(  ndime+2, inode) = prthe(1) !STRUCT%press(         ipoin) 
      phit(  ndime+2, inode) = STRUCT%press(         ipoin) + prthe(1)
    enddo
    !---------------------------------------------------------------------||---!
    !---------------------------------------------------------------------||---!
    call elmlen(ndime, pnode, elmar(pelty)%dercg, tragl, elcod, hnatu(pelty), hleng)
    hmini = hleng(ndime) !> hleng(ndime) is the smallest
    hinv  = 1.0_rp/hleng(ndime)
    !
    gauss_loop01: &
    do igaus = 1,pgaus
      !-------------------------------------------------------------------||---!
      call elmder(pnode, ndime, elmar(pelty)%deriv(1,1,igaus), elcod, cartd(1,1,igaus), detjm, xjacm, xjaci)
      STRUCT%xdvol(ielem,igaus) = elmar(pelty)%weigp(igaus) * detjm 
      !-------------------------------------------------------------------||---!
      !-------------------------------------------------------------------||---!
      xdphit(1:ndofn,1:ndime,igaus) = 0.0_rp
      do idofn = 1,ndofn
        do idime = 1,ndime
          xdphit(idofn,idime,igaus) = dot_product( phit(idofn,1:pnode), cartd(idime,1:pnode,igaus) )
        enddo
      enddo
      !STRUCT%xchrc(1,ielem,igaus,1) = sqrt( dot_product( xdphit(2,1:ndime,igaus), xdphit(2,1:ndime,igaus) )  ) 
      !STRUCT%xchrc(1,ielem,igaus,2) = sqrt( dot_product( xdphit(3,1:ndime,igaus), xdphit(3,1:ndime,igaus) )  )            
      !STRUCT%xchrc(1,ielem,igaus,3) = sqrt( dot_product( xdphit(4,1:ndime,igaus), xdphit(4,1:ndime,igaus) )  )            
      !
      do idofn = 1,ndofn
        xphit(idofn,igaus) = dot_product( phit(idofn,1:pnode), elmar(pelty)%shape(1:pnode,igaus) )
      enddo
      !-------------------------------------------------------------------||---!
      !-------------------------------------------------------------------||---!
      xgamme(igaus) = dot_product( elgamme(1:pnode), elmar(pelty)%shape(1:pnode,igaus) )
      call psound(xphit(1:ndofn,igaus), xgamme(igaus), xsound(igaus), xvelmo(igaus))
!print *, "xgamme, xsound",  xgamme(igaus), xsound(igaus)  
      !>  L = \Lambda^k * S_k^{-1} * dU/dr_k
      do idime = 1,ndime
       call SINVMTRX( Sxinv(1:ndofn,1:ndofn), xphit(1:ndofn,igaus), idime, xsound(igaus) )
       STRUCT%xchrc(1:ndofn,ielem,igaus,idime) = matmul( Sxinv(1:ndofn,1:ndofn), xdphit(1:ndofn,idime,igaus) )
      enddo
      !-------------------------------------------------------------------||---!
    enddo gauss_loop01
    !---------------------------------------------------------------------||---!
  end do elementary_loop
  !-----------------------------------------------------------------------||---!
  !
  !---------------------------------------------------------| gauss2nodes |---<! 
  do idofn = 1,ndofn
    do idime = 1,ndime
      call gauss2nodes(STRUCT%xchrc(idofn,1:nelem,1:mgaus,idime), STRUCT%xdvol(1:nelem,1:mgaus), STRUCT%chrc(idofn,1:npoin,idime)) 
    enddo
  enddo
  !-----------------------------------------------------------------------||---!
  !
  endif
  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !--------------------------------------------------------------| PRIVATE |---!
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine SMTRX(Ak, Wk, kk, cvel)
    implicit none
    real(rp), intent(inout) :: Ak(ndofn, ndofn)
    real(rp), intent(in)    :: Wk(ndofn), cvel
    integer(ip), intent(in) :: kk
    real(rp)    :: inv_c, inv_c2, inv_c_rho 
    integer(ip) :: ii
    !> W  = [rho, u1, u2, u3, p]^T 
    inv_c     = 1.0_rp/cvel
    inv_c2    = inv_c*inv_c
    inv_c_rho = 1.0_rp/Wk(1)*inv_c
    Ak(1:ndofn,1:ndofn)   = 0.0_rp
    !> |
    Ak(1,(/1_ip,ndime+2/)) = inv_c2*0.5_rp 
    Ak(1,         kk+1) = inv_c2
    Ak(ndime+2, (/1_ip,ndime+2/)) = 0.5_rp
    !> _
    Ak(kk+1,      1) = -0.5_rp*inv_c_rho
    Ak(kk+1,ndime+2) =  0.5_rp*inv_c_rho 
    !> |_ 
    do ii = 2,ndime+1
      Ak(ii,ii) = 1
    enddo    
    Ak(kk+1,kk+1) = 0.0_rp 
  end subroutine 
!!
 !subroutine SINVMTRX(Wk, kk, cvel)
 subroutine SINVMTRX(Ak, Wk, kk, cvel)
    implicit none
    real(rp), intent(out)   :: Ak(ndofn, ndofn)
    real(rp), intent(in)    :: Wk(ndofn), cvel
    integer(ip), intent(in) :: kk
    real(rp)    :: c2, c_rho, rho, vk 
    real(rp)    :: lambda(ndofn)  
    integer(ip) :: ii
    !> W  = [rho, u1, u2, u3, p]^T 
    !> S^{-1}_k * F^{k} S_k = \Lambda_k
    !> F^k = [rho*u_i*u_k+p, rho*u_k, (rho*e+p)*u_k]^T  
    rho   = Wk(1)
    vk    = Wk(kk+1)
    c2    = cvel*cvel
    c_rho = cvel*rho
    Ak(:,:) = 0.0_rp
    !> |
    Ak(1,ndime+2) =  1.0_rp 
    Ak(1,   kk+1) = -c_rho 
    Ak(ndime+2,   kk+1) = c_rho 
    Ak(ndime+2,ndime+2) = 1.0_rp 
    !> _
    Ak(kk+1,      1) =  c2
    Ak(kk+1,ndime+2) = -1.0_rp
    !> |_ 
    do ii = 2,ndofn
      Ak(ii,ii) = 1
    enddo    
    Ak(kk+1,kk+1) = 0.0_rp

    lambda(2:ndime+1) = vk
    lambda(        1) = vk - cvel 
    lambda(ndofn) = vk + cvel

    do ii = 1,ndofn
      Ak(1:ndofn,ii) = Ak(1:ndofn,ii)*lambda(1:ndofn)
    enddo
  end subroutine SINVMTRX
!!
  subroutine PMTRX(Ak, Wk, gamme, vkvk) 
    implicit none
    real(rp), intent(out) :: Ak(ndofn, ndofn)
    real(rp), intent(in)  :: Wk(ndofn), gamme, vkvk
    real(rp)    :: rho 
    integer(ip) :: ii

    rho = Wk(1)

    Ak(:,:) = 0.0_rp
    !> |_ 
    do ii = 2,ndime+1 
      Ak(     ii,ii) = 0.5_rp*rho
      Ak(     ii, 1) = 0.5_rp*Wk(ii)
      Ak(ndime+2,ii) = Wk(ii)*rho 
    enddo
    Ak(      1,      1) =  1.0_rp
    Ak(ndime+2,ndime+2) =  1.0_rp/(gamme-1.0_rp)
    Ak(ndime+2,      1) = -0.5_rp*vkvk*vkvk 
  end subroutine 
!!
  subroutine cons2phys(Uk, gamme, Wk)
    implicit none
    real(rp), intent(in)  :: Uk(ndofn), gamme  
    real(rp), intent(out) :: Wk(ndofn) 
    real(rp) :: aux(ndofn) 
    !> P = (gamma-1.0) * rho * (Et - 1/2 |V|) 
    !>   = (gamme-1.0) * (ENE - 0.5*sum(MOM(1:ndime)*MOM(1:ndime))/RHO)
    !Wk(1:ndime) = Uk(1:ndime)/Uk(ndime+1) 
    !Wk(ndime+1) = Uk(ndime+1) 
    !Wk(ndime+2) = (gamme-1.0)*(Uk(ndime+2) - 0.5/Uk(ndime+1)*sum(Uk(1:ndime)*Uk(1:ndime)))  
    !> [mom, rho, et] -> [rho, vel, p] 
    aux(      1)   = Uk(ndime+1) 
    aux(2:ndime+1) = Uk(1:ndime)/Uk(ndime+1) 
    aux(ndime+2)   = (gamme-1.0_rp)*(Uk(ndime+2) - 0.5_rp/Uk(ndime+1)*sum(Uk(1:ndime)*Uk(1:ndime))) 
    Wk = aux 
  end subroutine cons2phys
!!
  subroutine phys2cons(Wk, gamme, Uk)
    implicit none
    real(rp), intent(in)  :: Wk(ndofn), gamme
    real(rp), intent(out) :: Uk(ndofn)
    real(rp) :: aux(ndofn)
    !> rho*Et = P/(gamma-1) + 1/2*rho*|V| 
    !>        =   rho*Cv*T  + 1/2*rho*|V|
    !Uk(1:ndime) = Wk(1:ndime)*Wk(ndime+1) 
    !Uk(ndime+1) = Wk(ndime+1)
    !Uk(ndime+2) = Wk(ndime+2)/(gamme-1.0) + 0.5*Wk(ndime+1)*sum(Wk(1:ndime)*Wk(1:ndime))
    !> [rho, vel, p] -> [mom, rho, et]
    aux(1:ndime) = Wk(2:ndime+1)*Wk(1)
    aux(ndime+1) = Wk(1)
    aux(ndime+2) = Wk(ndime+2)/(gamme-1.0_rp) + 0.5_rp*Wk(1)*sum(Wk(2:ndime+1)*Wk(2:ndime+1))
    Uk = aux 
  end subroutine phys2cons  
!!
  subroutine csound(Uk, gamme, c, vmod)
    implicit none
    real(rp), intent(in)  :: Uk(ndofn), gamme
    real(rp), intent(out) :: c, vmod  
    !> c2 = gamma*P/rho = gamma(gamma-1)(Et-1/2|V|) 
    vmod = sum(Uk(1:ndime)/Uk(ndime+1)*Uk(1:ndime)/Uk(ndime+1)) 
    c    = Uk(ndime+2)/Uk(ndime+1) - 0.5_rp*vmod
    vmod = sqrt(vmod) 
    c    = sqrt(gamme*(gamme-1.0_rp)*c)
  end subroutine csound
!!
  subroutine psound(Wk, gamme, c, vmod)
    implicit none
    real(rp), intent(in)  :: Wk(ndofn), gamme
    real(rp), intent(out) :: c, vmod 
    !> c^2 = gamma*p/rho = gamma*(gamma-1)*cv*T = gamma*(gamma-1)*e
    !> W = [rho, vel, p] 
    c    = Wk(ndime+2)/Wk(1) 
    c    = sqrt( gamme*c )
    vmod = sqrt( dot_product(Wk(2:ndime+1), Wk(2:ndime+1)) )
  end subroutine psound
!!
  subroutine gauss2nodes(gprop, dvol, nprop) 
!    use def_domain, only: mgaus  
    implicit none

    real(rp), intent(in ) :: gprop(nelem, mgaus) 
    real(rp), intent(in ) ::  dvol(nelem, mgaus)
    real(rp), intent(out) :: nprop(npoin)
    real(rp)    :: gshape(mgaus)
    real(rp)    ::   gfac(mgaus)
    integer(ip) :: col, dimsize  

    col     = 1
    dimsize = 1 
    nprop(1:npoin) = 0.0_rp

    do ielem = 1,nelem
      pelty = ltype(ielem)

      pgaus = ngaus(pelty)
      gfac(1:pgaus) = dvol(ielem,1:pgaus) * gprop(ielem,1:pgaus)

      do inode = 1,nnode(pelty)
        gshape(1:pgaus) = elmar(pelty)%shape(inode,1:pgaus)

        ipoin = lnods(inode,ielem)
        nprop(ipoin) = nprop(ipoin) + dot_product(gshape(1:pgaus),gfac(1:pgaus))
      enddo
    enddo

    call rhsmod(dimsize, nprop(1:npoin))
    nprop(1:npoin) = nprop(1:npoin)/vmass(1:npoin)
  end subroutine 
!!
  subroutine getgammaker(Runiv, shecp, ggamme)
    use def_master,     only: ID_NASTAL, ID_CHEMIC, kfl_coupl
    use def_master,     only: wmean
    use mod_ker_proper, only: ker_proper
    implicit none
    real(rp), intent(in   ) :: Runiv
    real(rp), intent(inout) :: shecp(npoin)
    real(rp), intent(out)   :: ggamme(npoin)
!    real(rp)    :: dummy(ndime,ndime)
    integer(ip) :: dummi

    !> Specificic heat Cp
    if(kfl_prope/=0) then
      !call ker_proper('SPHEA','NPOIN',dummi,dummi,shecp(1:npoin),dummi,dummi,dummy(:,1),dummy(:,1))
      call ker_proper('SPHEA','NPOIN',dummi,dummi,shecp(1:npoin))
      !> Mean molecular weight (wmean) calcuted from chemic
      if(kfl_coupl(ID_NASTAL,ID_CHEMIC)==1) ggamme(1:npoin) = 1.0_rp/(1.0_rp - runiv/wmean(1:npoin,1)/shecp(1:npoin))

      !> thermal conductivity, kappa, used by diffussive matrix 
      !call ker_proper('CONDU','IGAUS',1_ip,ielem,xdith,pnode,1_ip)
    endif

    !> Visco
    !if(kfl_prope/=0) call ker_proper('VISCO','NPOIN',dummi,dummi,visco(:,icomp),dummi,dummi,xrano,xrano)
  end subroutine
!!
  subroutine getmaxmach(mom, rho, ene, gamme, mach)
    implicit none
    real(rp), intent(in ) :: mom(ndime, npoin), rho(npoin), ene(npoin), gamme(npoin)
    real(rp), intent(out) :: mach(npoin)
    real(rp) :: Uk(ndofn), c, vmod

    !if(INOTMASTER) then
      mach(1:npoin) = -666.666_rp

      do ipoin = 1,npoin
        Uk(1:ndofn) = (/ mom(1:ndime,ipoin), rho(ipoin), ene(ipoin) /)
        call csound(Uk, gamme(ipoin), c, vmod)
        if(c/=0.0_rp) mach(ipoin) = vmod/c
      end do
    !endif
  end subroutine
!!
  subroutine getvorticity(dvel, vort)
    implicit none
    real(rp), intent(in ) :: dvel(ndime,ndime)
    real(rp), intent(out) :: vort 
    real(rp) :: omega(3)

    omega(1:3) = 0.0_rp 
    omega(3) = dvel(2,1) - dvel(1,2) !> dv_y/r_x - dv_x/r_y 
    if(ndime==3) then 
    omega(2) = dvel(1,3) - dvel(3,1) !> dv_x/r_z - dv_z/r_x
    omega(1) = dvel(3,2) - dvel(2,3) !> dv_z/r_y - dv_y/r_z
    endif
 
    vort = sqrt( dot_product(omega(1:3), omega(1:3)) )
  end subroutine
!!
!==============================================================================! 
end module 
