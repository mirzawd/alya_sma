!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_elmmsudiff_coo(&
     kfl_ortho, kfl_shock,shock_tur, pnode,plapl,pgaus,gpvel,gpdif,&
     gprea,gptur,gpgrd, gprhs,gpden,gpsha,gpcar,gpvol,h, elunk, gphes, sreac, &
     gppro, gpprr, gppgr,gpmut,eltur,elwal,gpvis,elres_diff,elreswal_diff,&
     gprea_diff,gprhs_diff, gpreawal_diff,gprhswal_diff,gpvol_der, gpcar_der,h_der,ielem)
  !------------------------------------------------------------------------
  !****f* Turbul/tur_elmmsudiff_coo
  ! NAME 
  !    tur_elmmsudiff_coo
  ! DESCRIPTION
  !    
  !      This calculates the partial derivative of the residual w.r.t. constants of turbulent visco (dR/dD)
  ! USES
  ! USED BY
  !    
  !------------------------------------------------------------------------

  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,ntens,mnode
  use def_turbul, only       :  nturb_tur,&
       &                        kfl_taust_tur, staco_tur
  use mod_matrix
  use mod_tauadr, only       :  tauadr
  
  implicit none
  integer(ip), intent(in)    :: pnode,plapl,pgaus, kfl_ortho, kfl_shock,ielem
  real(rp),    intent(in)    :: gpvel(ndime,pgaus), gptur(nturb_tur,3,pgaus)
  real(rp),    intent(in)    :: gpdif(pgaus),gprea(pgaus)
  real(rp),    intent(in)    :: gprhs(pgaus), shock_tur
  real(rp),    intent(in)    :: gpden(pgaus)
  real(rp),    intent(in)    :: gpgrd(ndime,pgaus)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus), h, elunk(pnode, 2),sreac(pgaus)
  real(rp),    intent(in)    :: gphes(ntens,mnode,pgaus)             ! dNk/dxidxj
  real(rp),    intent(in)    :: gppro(pgaus), gpprr(pgaus), gppgr(ndime, pgaus)
  real(rp),    intent(in)    :: gpmut(pgaus)
  real(rp),    intent(in)    :: eltur(nturb_tur,pnode)
  real(rp),    intent(in)    :: elwal(pnode)
  real(rp),    intent(in)    :: gpvis(pgaus)
  real(rp),    intent(out)   :: elres_diff(pnode,ndime,pnode)
  real(rp),    intent(out)   :: elreswal_diff(pnode,ndime,pnode)
  
  real(rp),    intent(in)    :: gpvol_der(ndime,pnode,pgaus),gpcar_der(ndime,mnode,ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gprea_diff(ndime,pnode,pgaus),gprhs_diff(ndime,pnode,pgaus)
  real(rp),    intent(in)    :: gpreawal_diff(ndime,pnode,pgaus),gprhswal_diff(ndime,pnode,pgaus)
  real(rp),    intent(in)    :: h_der(ndime,pnode)
  
  integer(ip)                :: inode,jnode,idime,igaus,knode,kdime
  real(rp)                   :: fact2, resid(pnode), gpper(pnode), rhsit
  real(rp)                   :: tau, gpadv(pnode), gpnve
  real(rp)                   :: rhnve, gpad1(pnode)
  
  real(rp)                   :: elmat_diff(pnode,pnode,ndime,pnode),elrhs_diff(pnode,ndime,pnode)
  real(rp)                   :: tau_diff(ndime,pnode)
  real(rp)                   :: resid_diff(pnode,ndime,pnode),gpad1_diff(pnode,ndime,pnode)
  real(rp)                   :: gpadv_diff(pnode,ndime,pnode),gpper_diff(pnode,ndime,pnode)
  real(rp)                   :: rhsit_diff(ndime,pnode),xmuit_diff,fact2_diff
  real(rp)                   :: rhsitwal_diff(ndime,pnode),residwal_diff(pnode,ndime,pnode)
  real(rp)                   :: elmatwal_diff(pnode,pnode,ndime,pnode),elrhswal_diff(pnode,ndime,pnode)
  
    
  tau_diff = 0.0_rp
  
  elres_diff = 0.0_rp
  elreswal_diff = 0.0_rp
  
    elmat_diff = 0.0_rp
    elrhs_diff = 0.0_rp
    elmatwal_diff = 0.0_rp
    elrhswal_diff = 0.0_rp
  !-------------------------------------------------------------------
  !
  ! Assembly
  !
  !-------------------------------------------------------------------
  do igaus = 1,pgaus
        
    !
    ! Stabilization parameter without reactive (nonlinear) term :
    ! tau = 1.0_rp /( 4.0_rp*gpdif(igaus)/chale(2)/chale(2) + 2.0_rp*rhnve/chale(1) ) 
    !
    call vecnor(gpvel(1,igaus),ndime,gpnve,2_ip)      
    rhnve = gpden(igaus) * gpnve 
    call tauadr(&
         kfl_taust_tur,staco_tur,rhnve,gpdif(igaus),sreac(igaus),&
         h,h,tau) 
             
    do kdime = 1,ndime
      do knode = 1,pnode
        tau_diff(kdime,knode) = h_der(kdime,knode)/h*( 4.0_rp*gpdif(igaus)/(h*h) + 1/tau )*tau*tau
      enddo
    enddo
             
    !
    ! Calculus of residual resid and perturbation function Pi=gppre
    !
    ! Rj  = rho*Nj/dt + rho*a.grad(Nj) + s*Nj
    ! R2j = -grad(k).grad(Nj) - k*Lap(Nj)
    ! Pi  = Ni*(1-tau*s) + tau*a.grad(Ni)
    !
    do inode = 1,pnode
        
      resid(inode) = 0.0_rp
      gpad1(inode) = 0.0_rp
      do idime = 1,ndime
        gpad1(inode) = gpad1(inode) + gpvel(idime,igaus) * gpcar(idime,inode,igaus)   
      end do
      gpadv(inode) = gpad1(inode) * gpden(igaus)
      resid(inode) = resid(inode) + gpden(igaus)*gpad1(inode) + gprea(igaus) * gpsha(inode,igaus)
      gpper(inode) = ( gpsha(inode,igaus)*(1.0_rp-tau*sreac(igaus)) + tau*gpadv(inode)) * gpvol(igaus)        
      !gpper(inode) = (gpsha(inode,igaus) + tau*gpadv(inode)) * gpvol(igaus)
        
      do kdime = 1,ndime
        do knode = 1,pnode
          resid_diff(inode,kdime,knode) = 0.0_rp
          residwal_diff(inode,kdime,knode) = 0.0_rp
          gpad1_diff(inode,kdime,knode) = 0.0_rp
          do idime = 1,ndime
            gpad1_diff(inode,kdime,knode) = gpad1_diff(inode,kdime,knode) + gpvel(idime,igaus) * gpcar_der(idime,inode,kdime,knode,igaus)   
          end do
          gpadv_diff(inode,kdime,knode) = gpad1_diff(inode,kdime,knode) * gpden(igaus)
          resid_diff(inode,kdime,knode) = resid_diff(inode,kdime,knode) + gpden(igaus)*gpad1_diff(inode,kdime,knode) + gprea_diff(kdime,knode,igaus) * gpsha(inode,igaus)
          residwal_diff(inode,kdime,knode) = residwal_diff(inode,kdime,knode) + gpreawal_diff(kdime,knode,igaus) * gpsha(inode,igaus)
          gpper_diff(inode,kdime,knode) = (gpsha(inode,igaus) + tau*gpadv(inode)) * gpvol_der(kdime,knode,igaus) + &
                                          (tau_diff(kdime,knode)*gpadv(inode) + tau*gpadv_diff(inode,kdime,knode)) * gpvol(igaus)
                                              
        enddo !knode
      enddo !kdime
          
    end do !inode
        
    rhsit = gprhs(igaus)     
    do kdime = 1,ndime
      do knode = 1,pnode
        rhsit_diff(kdime,knode) = gprhs_diff(kdime,knode,igaus)
        rhsitwal_diff(kdime,knode) = gprhswal_diff(kdime,knode,igaus)
      enddo
    enddo
    ! 
    ! Diffusion Term
    !
    fact2 = gpvol(igaus) * gpdif(igaus)  ! Diffusion
    do jnode = 1,pnode
      do inode = 1,pnode
          do kdime = 1,ndime
            do knode = 1,pnode
              xmuit_diff = 0.0_rp
              fact2_diff = gpvol_der(kdime,knode,igaus) * gpdif(igaus)  ! Diffusion
              do idime = 1,ndime
                xmuit_diff = xmuit_diff + (gpcar_der(idime,jnode,kdime,knode,igaus) * gpcar(idime,inode,igaus) &
                                        + gpcar(idime,jnode,igaus) * gpcar_der(idime,inode,kdime,knode,igaus)) * fact2 &
                                        + gpcar(idime,jnode,igaus) * gpcar(idime,inode,igaus) * fact2_diff
              enddo
              elmat_diff(inode,jnode,kdime,knode) = elmat_diff(inode,jnode,kdime,knode) + xmuit_diff
            enddo !knode
          enddo !kdime
       end do !inode
    end do !jnode
    !
    ! Assembly of the matrix and rhs for ASGS, SUPG or Full OSS 
    !
    if( kfl_ortho /= 2 ) then 
       do kdime = 1,ndime
          do knode = 1,pnode
            do inode = 1,pnode
              do jnode = 1,pnode
                elmat_diff(inode,jnode,kdime,knode) = elmat_diff(inode,jnode,kdime,knode) + resid_diff(jnode,kdime,knode) * gpper(inode) + &
                                                      resid(jnode) * gpper_diff(inode,kdime,knode)
                elmatwal_diff(inode,jnode,kdime,knode) = elmatwal_diff(inode,jnode,kdime,knode) + residwal_diff(jnode,kdime,knode) * gpper(inode)
              enddo !jnode
              elrhs_diff(inode,kdime,knode) = elrhs_diff(inode,kdime,knode) + rhsit_diff(kdime,knode) * gpper(inode) + &
                                                                              rhsit * gpper_diff(inode,kdime,knode)
              elrhswal_diff(inode,kdime,knode) = elrhswal_diff(inode,kdime,knode) + rhsitwal_diff(kdime,knode) * gpper(inode)
            enddo!inode
                
          end do !knode
       end do !kdime
    end if !kfl_ortho
    
                
  end do !igaus

  
  
    do kdime = 1,ndime
      do knode = 1,pnode
        do inode = 1,pnode
          do jnode = 1,pnode
            elres_diff(inode,kdime,knode) = elres_diff(inode,kdime,knode) + elmat_diff(inode,jnode,kdime,knode)*eltur(1,jnode)
            elreswal_diff(inode,kdime,knode) = elreswal_diff(inode,kdime,knode) + elmatwal_diff(inode,jnode,kdime,knode)*eltur(1,jnode)
          enddo
          elres_diff(inode,kdime,knode) = elres_diff(inode,kdime,knode) - elrhs_diff(inode,kdime,knode)
          elreswal_diff(inode,kdime,knode) = elreswal_diff(inode,kdime,knode) - elrhswal_diff(inode,kdime,knode)
        enddo
      enddo !knode
    enddo !kdime
    
!     if (ielem == 5439) print *, "111xxx", elrhswal_diff(1,1,:)
!     if (ielem == 5439) print *, "111yyy", elrhswal_diff(1,2,:)
!     if (ielem == 5439) print *, "222xxx", elrhswal_diff(2,1,:)
!     if (ielem == 5439) print *, "222yyy", elrhswal_diff(2,2,:)
!     if (ielem == 5440) print *, "xxx", elmatwal_diff(4,3,1,:)
!     if (ielem == 5440) print *, "yyy", elmatwal_diff(4,3,2,:)
  
  

end subroutine tur_elmmsudiff_coo
