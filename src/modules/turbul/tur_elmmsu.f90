!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_elmmsu(&
     kfl_ortho, kfl_shock,shock_tur, pnode,plapl,pgaus,gpvel,gpdif,&
     gprea,gptur,gpgrd, gprhs,gpden,gpsha,gpcar,gpvol,&
     elmat,elrhs, chale, elunk, gphes, sreac, itask, gpres, gprec,&
     gpcon, gppro, gpprr, gppgr, grtup, conve)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_elmmsu
  ! NAME
  !   tur_elmsu
  ! DESCRIPTION
  !   This routine is an alternative to elmadr routine, 
  !   with the advantage that permits to compute SUPG stabilization
  !    Compute elemental matrix and rhs  with supg stabilization
  !    It is also capable of adding 
  !    Add ASGS stabilization, can give better results
  !    Add shock capturing techniques
  ! OUTPUT
  !    ELMAT ... LHS matrix for current Gauss point
  !    ELRHS ... RHS vector for current Gauss point
  ! USES
  ! USED BY
  !    tur_elmop2
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,ntens,mnode
  use def_turbul, only       :  nturb_tur,iunkn_tur,&
       &                        dtinv_tur, &
       &                        kfl_taust_tur, staco_tur,nbdfp_tur,pabdf_tur
  use def_master, only       :  kfl_lumped
  use mod_tauadr, only       :  tauadr
  use mod_ker_regularization, only :  dregularization, kfl_second, kfl_regularization

  implicit none
  integer(ip), intent(in)    :: pnode,plapl,pgaus, kfl_ortho, kfl_shock, itask
  real(rp),    intent(in)    :: gpvel(ndime,pgaus)
  real(rp),    intent(in)    :: gptur(nturb_tur,nbdfp_tur+1,pgaus)
  real(rp),    intent(in)    :: conve(ndime,pgaus)
  real(rp),    intent(in)    :: gpdif(pgaus),gprea(pgaus)
  real(rp),    intent(in)    :: gprhs(pgaus), shock_tur
  real(rp),    intent(in)    :: gpden(pgaus)
  real(rp),    intent(in)    :: gpgrd(ndime,pgaus)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus), chale(2), elunk(pnode, 2),sreac(pgaus)
  real(rp),    intent(in)    :: gphes(ntens,mnode,pgaus)             ! dNk/dxidxj
  real(rp),    intent(in)    :: gppro(pgaus), gpprr(pgaus), gppgr(ndime, pgaus)
  real(rp),    intent(out)   :: elmat(pnode,pnode),elrhs(pnode)
  real(rp),    intent(out)   :: gprec(pgaus), gpcon(pgaus), gpres(pgaus), grtup(ndime, pgaus) ! Projections
  integer(ip)                :: inode,jnode,idime,igaus,itime
  real(rp)                   :: fact1,fact2,fact3, resid(pnode), gpper(pnode), rhsit , uscoe
  real(rp)                   :: xmuit, tau, gpadv(pnode), gpnve, grvgr, resi2(pnode)
  real(rp)                   :: rhnve, gpad1(pnode),gpad2(pnode), rresi, ugrau, grtur(ndime), grnor,  xmui3, switc
  real(rp)                   :: Dsupg,umbra, CD, SD, factt, gplap, gppe2(pnode)
  real(rp)                   :: dtinv_res, facto_time, fact2_time
  !
  ! Initialization
  !
  factt = 0.75_rp
  umbra = 1.0e-8_rp
  CD    = 0.0_rp
  SD    = 0.0_rp
  do inode = 1,pnode
     elrhs(inode) = 0.0_rp
     do jnode = 1,pnode
        elmat(jnode,inode) = 0.0_rp
     end do
  end do
  !----------------------------------------------------------------------
  !
  ! Time step
  !
  !----------------------------------------------------------------------

  if( kfl_lumped == 2 ) then
     dtinv_res = 0.0_rp
  else
     dtinv_res = dtinv_tur
  end if

  if( itask /= 4 ) then 

     !-------------------------------------------------------------------
     !
     ! Assembly
     !
     !-------------------------------------------------------------------
     facto_time = 1.0_rp
     fact2_time = 1.0_rp
        
     Gauss: do igaus = 1,pgaus
        
        if (kfl_regularization.and.kfl_second) then
           facto_time = dregularization ( gptur(iunkn_tur,1,igaus))
           fact2_time = 0.0_rp ! facto_time
        end if
        !
        ! Stabilization parameter without reactive (nonlinear) term :
        ! tau = 1.0_rp /( 4.0_rp*gpdif(igaus)/chale(2)/chale(2) + 2.0_rp*rhnve/chale(1) ) 
        !
!         call vecnor(conve(1,igaus),ndime,gpnve,2_ip)
        call vecnor(gpvel(1,igaus),ndime,gpnve,2_ip)
        rhnve = gpden(igaus) * gpnve 
        call tauadr(&
             kfl_taust_tur,staco_tur,rhnve,gpdif(igaus),sreac(igaus),&
             chale(1),chale(2),tau, dtinv_res * pabdf_tur(1))   
        !
        ! Calculus of residual resid and perturbation function Pi=gppre
        !
        ! Rj  = rho*Nj/dt + rho*a.grad(Nj) + s*Nj
        ! R2j = -grad(k).grad(Nj) - k*Lap(Nj)
        ! Pi  = Ni*(1-tau*s) + tau*a.grad(Ni)
        !
        do inode = 1,pnode
           resid(inode) = facto_time*gpden(igaus) * gpsha(inode,igaus) * dtinv_res * pabdf_tur(1)
           gpad1(inode) = 0.0_rp
           gpad2(inode) = 0.0_rp
           grvgr        = 0.0_rp
           gplap        = 0.0_rp
           do idime = 1,ndime
!               gpad1(inode) = gpad1(inode) + conve(idime,igaus) * gpcar(idime,inode,igaus)
              gpad1(inode) = gpad1(inode) + gpvel(idime,igaus) * gpcar(idime,inode,igaus) ! perturbation function
              gpad2(inode) = gpad2(inode) + conve(idime,igaus) * gpcar(idime,inode,igaus) ! residual
              grvgr        = grvgr + gpgrd(idime,igaus) * gpcar(idime,inode,igaus)          
              gplap        = gplap + gphes(idime,inode,igaus)
           end do
           gpadv(inode) = gpad1(inode) * gpden(igaus)
           resid(inode) = resid(inode) + gpden(igaus)*gpad2(inode) + gprea(igaus) * gpsha(inode,igaus)
           resi2(inode) = - grvgr - gplap * gpdif(igaus)
           gpper(inode) = ( gpsha(inode,igaus)*(1.0_rp-tau*sreac(igaus)) + tau*gpadv(inode)) * gpvol(igaus)
           !gpper(inode) = (gpsha(inode,igaus) + tau*gpadv(inode)) * gpvol(igaus)        
        end do
        rhsit =  gprhs(igaus)
        do itime = 3,nbdfp_tur+1
           rhsit = rhsit - fact2_time*dtinv_res * pabdf_tur(itime-1) * gpden(igaus) * gptur(iunkn_tur,itime,igaus)
        end do
        
        !
        ! Shock capturing
        !
        if( kfl_shock /= 0 ) then
           rresi = rhsit - (dtinv_res * pabdf_tur(1) * gpden(igaus) + gprea(igaus)) * gptur(iunkn_tur,1,igaus) 
           ugrau = 0.0_rp
           do idime = 1,ndime
              grtur(idime) = 0.0_rp
              do inode = 1, pnode
                 grtur(idime) = grtur(idime) + elunk(inode, 1)*gpcar(idime,inode, igaus)
              end do
              ugrau = ugrau + ( gpden(igaus)*gpvel(idime, igaus) &
                   -gpgrd(idime, igaus))* grtur(idime)                      
           end do
           rresi =rresi - ugrau
           if( plapl > 0 ) then ! laplacian
              do inode = 1, pnode
                 gplap =  0.0_rp
                 do idime=1,ndime           
                    gplap = gplap + gphes(idime,inode,igaus)
                 end do
                 rresi = rresi + elunk(inode, 1)*gplap*gpdif(igaus)
              end do
           end if
           grnor = 0.0_rp
           do idime =1, ndime
              grnor = grnor + grtur(idime)*grtur(idime)
           end do
           grnor = sqrt(grnor)
           if(grnor > umbra ) then
              uscoe = abs(rresi/grnor)
              ! cross diffusion 
              CD = max(factt*0.5_rp*shock_tur*chale(2)*uscoe - gpdif(igaus),0.0_rp)              
              ! Orthogonal gradient
              switc = 0.0_rp 
              do idime =1, ndime
                 switc = switc + (gppgr(idime, igaus) - grtur(idime))*(gppgr(idime, igaus) - grtur(idime))
              end do
              switc = sqrt(switc)/grnor
              !! uncomment the following line for shoch with  Orthogonal gradient
!              CD = max(factt*0.5_rp*shock_tur*chale(2)*(0.8_rp*rhnve + 0.0*chale(2)*gprea(igaus))*switc - gpdif(igaus),0.0_rp)              
              ! streamline diffusion tau rho^2*u^2
              if ( gpnve > umbra.and.kfl_shock==2 ) then ! anisotropic
                 Dsupg = tau*rhnve*rhnve !supg 
!                 SD    = max(CD-Dsupg,0.0_rp)  -CD ! streamline diffusion
                 SD    = -CD  !only crosswind
                 SD    =  SD/gpnve/gpnve   
                 !              SD=  0.0_rp
              else
                 SD    = 0.0_rp
              end if
           else 
              CD = 0.0_rp
              SD = 0.0_rp
           end if
        end if
        ! 
        ! Diffusion Term
        !
        fact2 = gpvol(igaus) * (gpdif(igaus)+CD)  ! Diffusion
        fact3 = SD * gpvol(igaus)                 ! Streamline negative diffusion
!yor! optim [originally loop on line 158] simplified loop and memory access :
        do jnode = 1,pnode
           do inode = 1,pnode
              xmuit = 0.0_rp
              do idime = 1,ndime
                 xmuit = xmuit + gpcar(idime,jnode,igaus) * gpcar(idime,inode,igaus)
              end do
              xmuit = xmuit * fact2
              xmui3 = gpad1(inode) * gpad1(jnode) * fact3
              elmat(inode,jnode) = elmat(inode,jnode) + xmuit + xmui3
           end do
        end do

        !
        ! Assembly of the matrix and rhs
        !
        if( kfl_ortho /= 2 ) then 
           !
           ! ASGS, SUPG or Full OSS 
           !
           do inode = 1,pnode
              gppe2(inode) = gpper(inode) - gpsha(inode,igaus) * gpvol(igaus)
              do jnode = 1,pnode
                 elmat(inode,jnode) = elmat(inode,jnode) &
                      + resid(jnode) * gpper(inode) &        
                      + resi2(jnode) * gppe2(inode)
              end do
              elrhs(inode) = elrhs(inode) + rhsit * gpper(inode) &
                   + gppro(igaus) * gppe2(inode)

           end do
        else
           ! 
           ! Split OSS
           !
           do inode = 1,pnode
              do jnode = 1,pnode
                 elmat(inode, jnode) = elmat(inode, jnode) + resid(jnode)* gpsha(inode,igaus)*gpvol(igaus) &
                      + gpden(igaus)*gpad2(jnode)*tau*gpadv(inode) * gpvol(igaus)  &
                      - gprea(igaus)*gpsha(jnode,igaus)* gpsha(inode,igaus)*tau*sreac(igaus)*gpvol(igaus)                       
              end do
              elrhs(inode)=elrhs(inode) +rhsit*gpsha(inode,igaus)*gpvol(igaus)*&
                   (1.0_rp-tau*sreac(igaus)) +&
                   gppro(igaus)*gpadv(inode)*gpvol(igaus) + &
                   gpprr(igaus)*tau*sreac(igaus)*gpsha(inode,igaus)*gpvol(igaus) 

           end do
        end if
        
        if( kfl_lumped == 1.and..not.kfl_second ) then
           ! 
           ! Lumped mass evolution matrix
           !
           do inode = 1,pnode
              fact1 = gpvol(igaus) * gpden(igaus) * gpsha(inode,igaus) * dtinv_res
              ! add lumped
              elmat(inode, inode) = elmat(inode, inode) + fact1
              do itime = 3,nbdfp_tur+1
                 elrhs(inode) = elrhs(inode) + fact1 * gptur(iunkn_tur, itime, igaus) * pabdf_tur(itime-1)
              end do
              ! substract consistent
              !
              ! The next line is only valid for CN or BDF1 - ELSE elunk would need to include previous values or use eltur
              ! In any case, Lumped mass evolution matrix only makes sense with local time step and in tur_outerr local is only
              ! allowed with 1st order time accuracy 
              elrhs(inode) = elrhs(inode) - fact1 * elunk(inode, 2) * pabdf_tur(2)
              do jnode =1, pnode
                 elmat(inode, jnode) = elmat(inode, jnode) - fact1 * gpsha(jnode, igaus) * pabdf_tur(1)
              end do
           end do
        end if !lumped
        
     end do Gauss

  else if( itask == 4 ) then 

     !-------------------------------------------------------------------
     !
     ! Compute residuals for orthogonal projection
     !
     !-------------------------------------------------------------------

     do igaus = 1,pgaus
        !    reactive term
        gprec(igaus) = gprea(igaus)*gptur(iunkn_tur,1, igaus) ! linearized

        ugrau =0.0_rp
        do idime =1, ndime
           grtur(idime) = 0.0_rp
           do inode =1, pnode
              grtur(idime)=grtur(idime) + elunk(inode, 1)*gpcar(idime,inode, igaus)
           end do
           ugrau = ugrau + ( gpden(igaus)*gpvel(idime, igaus) &
                -gpgrd(idime, igaus))* grtur(idime)                      
        end do
        ! convective term (plus a diffusive part )
        gpcon(igaus) =  ugrau 
        ! residual
        gpres(igaus) = gprhs(igaus) + dtinv_res*gpden(igaus)* &
             (gptur(iunkn_tur,3,igaus)-gptur(iunkn_tur,1, igaus) ) &
             - gprec(igaus) -  gpcon(igaus)

        gprec(igaus) = gprhs(igaus) - gprec(igaus)
        ! proyection of turbulence gradient
        grtup(1:ndime, igaus) = grtur(1:ndime)
     end do

  end if

end subroutine tur_elmmsu

