!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_elmmsudiff_tvcon(&
     kfl_ortho, kfl_shock,shock_tur, pnode,plapl,pgaus,gpvel,gpdif,&
     gprea,gptur,gpgrd, gprhs,gpden,gpsha,gpcar,gpvol,&
     chale, elunk, gphes, sreac, itask, &
     gppro, gpprr, gppgr,gpmut,eltur,elwal,gpvis,elresdiff)
  !------------------------------------------------------------------------
  !****f* Turbul/tur_elmmsudiff_tvcon
  ! NAME 
  !    tur_elmmsudiff_tvcon
  ! DESCRIPTION
  !    
  !      This calculates the partial derivative of the residual w.r.t. constants of turbulent visco (dR/dD)
  ! USES
  ! USED BY
  !    
  !------------------------------------------------------------------------

  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,ntens,mnode
  use def_turbul, only       :  nturb_tur,iunkn_tur, kfl_taust_tur, staco_tur,kfl_clipp_tur
  use mod_matrix
  use mod_tauadr, only       :  tauadr
  
  implicit none
  integer(ip), intent(in)    :: pnode,plapl,pgaus, kfl_ortho, kfl_shock, itask
  real(rp),    intent(in)    :: gpvel(ndime,pgaus), gptur(nturb_tur,3,pgaus)
  real(rp),    intent(in)    :: gpdif(pgaus),gprea(pgaus)
  real(rp),    intent(in)    :: gprhs(pgaus), shock_tur
  real(rp),    intent(in)    :: gpden(pgaus)
  real(rp),    intent(in)    :: gpgrd(ndime,pgaus)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus), chale(2), elunk(pnode, 2),sreac(pgaus)
  real(rp),    intent(in)    :: gphes(ntens,mnode,pgaus)             ! dNk/dxidxj
  real(rp),    intent(in)    :: gppro(pgaus), gpprr(pgaus), gppgr(ndime, pgaus)
  real(rp),    intent(in)    :: gpmut(pgaus)
  real(rp),    intent(in)    :: eltur(nturb_tur,pnode)
  real(rp),    intent(in)    :: elwal(pnode)
  real(rp),    intent(in)    :: gpvis(pgaus)
  real(rp),    intent(out)   :: elresdiff(pnode)
  integer(ip)                :: inode,idime,igaus
  real(rp)                   :: gpper(pnode)
  real(rp)                   :: tau, gpadv(pnode), gpnve
  real(rp)                   :: rhnve, gpad1(pnode)
  real(rp)                   :: CD
  real(rp)                   :: grtur_forw(pgaus,ndime),gpstadiff,resid_st
  real(rp)                   :: F1,gradk(ndime),gradw(ndime),grakw,sigw2,gpkin,gpome,xsmal,gpwal
  real(rp)                   :: arg1,c1,c2,c3,betas
  
  !
  ! Initialization
  !
  CD    = 0.0_rp
  sigw2 = 0.856_rp 
  xsmal = 1.0e-10_rp
  betas = 0.09_rp
  do inode = 1,pnode
     elresdiff(inode) = 0.0_rp
  end do
  grtur_forw = 0.0_ip

  !
  ! Calculate some variables
  !
  do igaus = 1, pgaus
    do idime=1,ndime
      do inode=1,pnode
	grtur_forw(igaus,idime) = grtur_forw(igaus,idime) + gpcar(idime,inode,igaus)*elunk(inode, 1)
      enddo
    enddo
  enddo
  
  !-------------------------------------------------------------------
  !
  ! Assembly
  !
  !-------------------------------------------------------------------
 
     do igaus = 1,pgaus
        
        gpkin = max(0.0_rp, gptur(1,1,igaus))
        gpome = gptur(2,1,igaus)
        gpwal = 0.0_rp
        do inode = 1,pnode
          gpwal = gpwal + gpsha(inode,igaus) * elwal(inode)
        end do
        
        gradk = 0.0_rp
        gradw = 0.0_rp
        grakw = 0.0_rp
        do inode = 1,pnode
          do idime = 1,ndime
            gradk(idime) = gradk(idime) + eltur(1,inode) * gpcar(idime,inode,igaus)
            gradw(idime) = gradw(idime) + eltur(2,inode) * gpcar(idime,inode,igaus)
          end do
        end do
        do idime = 1,ndime
          grakw = grakw + gradk(idime) * gradw(idime)
        end do
        if( ( gpome > xsmal ) .or. ( kfl_clipp_tur > 1 ) ) then
          CD = max( 1.0e-10_rp , 2.0_rp*gpden(igaus)*sigw2*grakw/gpome )
        else
          CD = max( 1.0e-10_rp , 2.0_rp*gpden(igaus)*sigw2*grakw/xsmal )
        end if
        !
        ! F1 function. k-w model corresponds to F1 = 1
        !
        if( ( gpome > xsmal ) .or. ( kfl_clipp_tur > 1 ) ) then
          c1   = sqrt( gpkin ) / ( betas * gpome * gpwal )
          c2   = 500.0_rp * gpvis(igaus) / ( gpden(igaus) * gpwal * gpwal * gpome )
          c3   = 4.0_rp * gpden(igaus) * sigw2 * gpkin / ( CD * gpwal * gpwal )
          arg1 = min ( max (  c1 , c2 ) , c3 )
        else
          arg1 = 4.0_rp * gpden(igaus) * sigw2 * gpkin / ( CD * gpwal * gpwal )
        end if
        F1 = tanh( arg1**4.0_rp )
     
        !
        ! Stabilization parameter without reactive (nonlinear) term :
        ! tau = 1.0_rp /( 4.0_rp*gpdif(igaus)/chale(2)/chale(2) + 2.0_rp*rhnve/chale(1) ) 
        !
        call vecnor(gpvel(1,igaus),ndime,gpnve,2_ip)      
        rhnve = gpden(igaus) * gpnve 
        call tauadr(&
             kfl_taust_tur,staco_tur,rhnve,gpdif(igaus),sreac(igaus),&
             chale(1),chale(2),tau)   
        !
        !  gpstadiff
        !
        gpstadiff = -4.0_rp*tau*tau*F1*gpmut(igaus)/(chale(2)*chale(2))
        !
        ! Calculus of residual_st
        !
        resid_st =  gprea(igaus)*gptur(iunkn_tur,1,igaus) !- gprhs(igaus)
        do idime=1,ndime
          resid_st = resid_st + gpden(igaus)*gpvel(idime,igaus)*grtur_forw(igaus,idime)
        end do
        !
        ! calculate d(perturbation)/d(des_var)
        !
        do inode=1,pnode
          gpad1(inode) =  0.0_rp
          do idime=1,ndime
            gpad1(inode) = gpad1(inode) + gpvel(idime,igaus)*gpcar(idime,inode,igaus)           
          end do
          gpadv(inode)=gpad1(inode)*gpden(igaus)
          gpper(inode)  = ( gpstadiff*gpadv(inode) ) * gpvol(igaus)        
        end do           
        !    
        ! Assembly: term related to d( mu + mut)/dD dN/dxi * dK/dxi
        !
        do inode = 1, pnode
          do idime= 1, ndime
            elresdiff(inode) = elresdiff(inode) + F1*gpmut(igaus) * grtur_forw(igaus,idime)*gpcar(idime, inode,igaus)*gpvol(igaus)
          end do
        enddo  
        !
        ! Assembly: term related to d(tau)/dk * vi * dN/dxi * ( T/dt + vi * dT/dxi)
        ! 
        do inode =1, pnode
          elresdiff(inode) = elresdiff(inode) + resid_st*gpper(inode)
        end do
        
     end do !igaus
  

end subroutine tur_elmmsudiff_tvcon
