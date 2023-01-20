!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_elmmsu_der_all(&
     kfl_ortho, kfl_shock,shock_tur, pnode,plapl,pgaus,gpvel,gpdif,gprea,&
     drea_dtur,dsou_dtur,dsou_dvel,gpmut,dgpmut,dgpmut_dvel,&
     gptur,gpgrd, gprhs,gpden,gpsha,gpcar,gpvol,&
     elmat,elrhs, chale, gphes, sreac,eltur,ielem)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_elmmsu_der_all
  ! NAME
  !   tur_elmmsu_der_all
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
!  use def_domain, only       :  lnods
  use def_turbul, only       :  nturb_tur,iunkn_tur, dtinv_tur, kfl_taust_tur,&
                                staco_tur,param_tur,TUR_SPALART_ALLMARAS
!  use def_turbul, only       :  kfl_fixno_tur, Rhsadjtur_tur
  use def_master, only       :  kfl_lumped,RhsadjNas_tur,RhsadjTur_nsi,kfl_coupl,ID_NASTIN,ID_TURBUL
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
  real(rp),    intent(in)    :: gpvol(pgaus), chale(2),sreac(pgaus)
  real(rp),    intent(in)    :: gphes(ntens,mnode,pgaus)             ! dNk/dxidxj

  real(rp),    intent(in)    :: gpmut(pgaus)
  real(rp),    intent(in)    :: dgpmut(nturb_tur,pnode,pgaus)
  real(rp),    intent(in)    :: dgpmut_dvel(ndime,pnode,pgaus)  
  real(rp),    intent(in)    :: drea_dtur(nturb_tur,pnode,pgaus)
  real(rp),    intent(in)    :: dsou_dtur(nturb_tur,pnode,pgaus)
  real(rp),    intent(in)    :: dsou_dvel(ndime,pnode,pgaus)
  real(rp),    intent(in)    :: eltur(nturb_tur,pnode,3)
  
  
  real(rp),    intent(inout) :: elmat(pnode,pnode),elrhs(pnode)
  integer(ip)                :: inode,jnode,idime,igaus,jdime,iturb_tur
!  integer(ip)                :: inode,jnode,idime,igaus,jdime,ipoin,iturb_tur
  real(rp)                   :: fact1,fact2, resid(pnode), gpper(pnode)
  real(rp)                   :: xmuit, tau, gpadv(pnode), gpnve
  real(rp)                   :: rhnve, gpad1(pnode), grtur(ndime)
  real(rp)                   :: dtau_dtur(nturb_tur,pnode),dtau_dvel(pnode)
  real(rp)                   :: dgpper_dtur(pnode,pnode),resid_st,elmat_aux(pnode,pnode),dgpper_du(pnode,pnode)
  real(rp)                   :: dgpmut_dtur(nturb_tur,pnode,pgaus),elmat_sou(pnode,pnode),diff_param_tur(nturb_tur)
  
    
  if (TUR_SPALART_ALLMARAS) then !!!!! For spalart the kernel values have to be modified
    diff_param_tur(1) = param_tur(3)
  else
    do iturb_tur = 1,nturb_tur
      diff_param_tur(iturb_tur) = param_tur(iturb_tur)
    enddo
  endif
  !
  ! Initialization
  !
  do inode = 1,pnode
    do igaus = 1,pgaus
      if (TUR_SPALART_ALLMARAS) then !!!!! For spalart the kernel values have to be modified
        dgpmut_dtur(1,inode,igaus) = gpsha(inode,igaus)
      else
        do iturb_tur = 1,nturb_tur
          dgpmut_dtur(iturb_tur,inode,igaus) = dgpmut(iturb_tur,inode,igaus)
        enddo
      endif
    enddo
  enddo
  
  if (iunkn_tur == 1) then 
!     Rhsadjtur_tur(ielem)%a(1,:) = 0.0_rp
    RhsadjNas_tur(ielem)%a = 0.0_rp
  else
!     Rhsadjtur_tur(ielem)%a(2,:) = 0.0_rp
  endif
  !
  ! For the case of adjoint 
  !
  elrhs = 0.0_rp ! elrhs starts from zero for adjoint but elmat is modified
   
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
	  chale(1),chale(2),tau)   
    !
    ! Derivatives of stabilization parameter w.r.t the k and w
    !
    do iturb_tur = 1,nturb_tur
      do inode =1, pnode
        dtau_dtur(iturb_tur,inode) = -4.0_rp*gpden(igaus)*dgpmut_dtur(iturb_tur,inode,igaus)*tau*tau/(chale(2)*chale(2)*diff_param_tur(iturb_tur))
      enddo
    enddo
    !
    ! Calculute gradient of k and w
    !
    do idime = 1, ndime
      grtur(idime) = 0.0_rp
      do inode =1, pnode
        grtur(idime) = grtur(idime) +  eltur(iunkn_tur,inode,1)*gpcar(idime,inode,igaus)
      end do                   
    end do
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!        ELMAT   ---->   dRk/dK  or  dRw/dW           !!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
    elmat_sou = 0.0_rp
    !
    ! Calculus of d(residual_st)/d(unk) AND perturbation 
    ! unk = k for kin equation and unk = w for omega equation
    !    
    do inode = 1,pnode
      gpad1(inode) = 0.0_rp
      do idime = 1,ndime
        gpad1(inode) = gpad1(inode) + gpvel(idime,igaus) * gpcar(idime,inode,igaus)           
      end do
      gpadv(inode) = gpad1(inode) * gpden(igaus)
      resid(inode) = drea_dtur(iunkn_tur,inode,igaus) + dsou_dtur(iunkn_tur,inode,igaus)
      gpper(inode) = ( gpsha(inode,igaus) + tau*gpadv(inode) ) * gpvol(igaus)        
    end do
    !
    ! Calculus of residual_st AND d(perturbation)/d(unk)
    !
    resid_st = gprea(igaus)*gptur(iunkn_tur,1,igaus) - gprhs(igaus)
    do idime=1,ndime
      resid_st = resid_st + gpden(igaus)*gpvel(idime,igaus)*grtur(idime)         
    end do
    do inode = 1,pnode
      gpad1(inode) = 0.0_rp
      do idime = 1,ndime
        gpad1(inode) = gpad1(inode) + gpvel(idime,igaus) * gpcar(idime,inode,igaus)           
      end do
      gpadv(inode) = gpad1(inode) * gpden(igaus)
      do jnode=1,pnode
        dgpper_dtur(inode,jnode) = dtau_dtur(iunkn_tur,jnode)*gpadv(inode) * gpvol(igaus)        
      enddo !jnode
    end do
    ! 
    ! Diffusion Term
    !        
    do inode = 1,pnode
      do jnode = 1,pnode
        fact2 = gpvol(igaus) * gpden(igaus)*dgpmut_dtur(iunkn_tur,jnode,igaus)/diff_param_tur(iunkn_tur)  ! turbulent visc derivatives
        xmuit = 0.0_rp
        do idime = 1,ndime
          xmuit = xmuit + grtur(idime) * gpcar(idime,inode,igaus)
        end do
!         elmat(inode,jnode) = elmat(inode,jnode) + xmuit * fact2 
        elmat_sou(inode,jnode) = elmat_sou(inode,jnode) + xmuit * fact2
      end do
    end do
    !
    ! Assembly of the matrix for ASGS, SUPG or Full OSS
    ! 
    if( kfl_ortho /= 2 ) then 
      do inode = 1,pnode
        do jnode = 1,pnode
!           elmat(inode,jnode) = elmat(inode,jnode) + resid(jnode) * gpper(inode) + dgpper_dtur(inode,jnode) * resid_st
          elmat_sou(inode,jnode) = elmat_sou(inode,jnode) + resid(jnode) * gpper(inode) + dgpper_dtur(inode,jnode)*resid_st
        end do
      end do
    end if
    !
    ! No time term have been added up to now: add Galerkin term
    !
    if ( kfl_lumped == 2 ) then
        do inode = 1,pnode
            fact1 = gpvol(igaus) * gpden(igaus) * dtinv_tur
            elmat(inode, inode) = elmat(inode, inode) + fact1 * gpsha(inode,igaus)
            elrhs(inode) = elrhs(inode) + fact1 * gpsha(inode,igaus) * eltur(iunkn_tur,inode,3)  !elunk(inode, 2)
        end do
    endif
    !
    ! Additional terms to the RHSID
    !          
!     do inode = 1, pnode
!       do jnode = 1, pnode
!         elrhs(inode) = elrhs(inode) - elmat_sou(jnode,inode)*eltur(iunkn_tur,jnode,3)
!       enddo
!     enddo  
    
    
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!  Rhsadjtur_tur(1)   <------  Transpose[dRk/dw] [Lambda_k]     to be sent to omega      !!!!!!!
  !!!!!!!!  Rhsadjtur_tur(2)   <------  Transpose[dRw/dk] [Lambda_w]     to be sent to kin        !!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
!     elmat_aux = 0.0_rp   
!     !
!     ! calculus of d(residual_st)/d(unk1) and perturbation function gppre 
!     ! unk1 = w for k equation and unk1 = k for omega equation
!     !
!     do inode = 1,pnode
!       gpad1(inode) = 0.0_rp
!       do idime = 1,ndime
! 	gpad1(inode) = gpad1(inode) + gpvel(idime,igaus) * gpcar(idime,inode,igaus)           
!       end do
!       gpadv(inode) = gpad1(inode) * gpden(igaus)
!       if (iunkn_tur == 1) then 
! 	resid(inode) = drea_dtur(2,inode,igaus) - dsou_dtur(2,inode,igaus)
!       else
! 	resid(inode) = drea_dtur(1,inode,igaus) - dsou_dtur(1,inode,igaus)
!       endif
!       gpper(inode) = ( gpsha(inode,igaus) + tau*gpadv(inode) ) * gpvol(igaus)        
!     end do
!     !
!     ! Calculus of residual_st AND d(perturbation)/d(unk1)
!     !
!     resid_st = gprea(igaus)*gptur(iunkn_tur,1,igaus) - gprhs(igaus)
!     do idime=1,ndime
! 	resid_st = resid_st + gpden(igaus)*gpvel(idime,igaus)*grtur(idime)
!     end do
!     do inode = 1,pnode
!       gpad1(inode) = 0.0_rp
!       do idime = 1,ndime
! 	gpad1(inode) = gpad1(inode) + gpvel(idime,igaus) * gpcar(idime,inode,igaus)           
!       end do
!       gpadv(inode) = gpad1(inode) * gpden(igaus)
!       do jnode=1,pnode
! 	    dgpper_dtur(inode,jnode) = dtau_dtur(iunkn_tur,jnode)*gpadv(inode) * gpvol(igaus)
!       enddo !jnode
!     end do
!     ! 
!     ! Diffusion Term
!     !
!     do inode = 1,pnode
!       do jnode = 1,pnode
! 	    fact2 = gpvol(igaus) * dgpmut_dtur(iunkn_tur,jnode,igaus)/diff_param_tur(iunkn_tur)  ! turbulent visc derivatives
! 	    xmuit = 0.0_rp
! 	    do idime = 1,ndime
! 	      xmuit = xmuit + grtur(idime) * gpcar(idime,inode,igaus)
! 	    end do
! 	    elmat_aux(inode,jnode) = elmat_aux(inode,jnode) + xmuit * fact2
!       end do
!     end do
!     !
!     ! Assembly of the matrix for ASGS, SUPG or Full OSS
!     !
!     if( kfl_ortho /= 2 ) then 
!       do inode = 1,pnode
! 	do jnode = 1,pnode
! 	  elmat_aux(inode,jnode) = elmat_aux(inode,jnode) + resid(jnode) * gpper(inode) + dgpper_dtur(inode,jnode)*resid_st
! 	end do
!       end do
!     end if
!     !
!     ! production of the Transpose[dRk/dw] [Lambda_k] to be sent to omega
!     ! production of the Transpose[dRw/dk] [Lambda_w] to be sent to kin
!     !
!     do inode = 1, pnode
!       do jnode = 1, pnode
! 	Rhsadjtur_tur(ielem)%a(iunkn_tur,inode) = Rhsadjtur_tur(ielem)%a(iunkn_tur,inode) + elmat_aux(jnode,inode)*eltur(iunkn_tur,jnode,3)
!       enddo
!     enddo  
                 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!  RhsadjNas_tur   <-------  Transpose[dRk/du] [Lambda_k]    to be sent to nastin   !!!!!!!
    !!!!!!!!  RhsadjNas_tur   <-------  Transpose[dRw/du] [Lambda_w]    to be sent to nastin   !!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Coupling with nastin
    !
    if (kfl_coupl(ID_TURBUL,ID_NASTIN) >= 1 ) then      
      !
      ! loop over ndime to calculate RhsadjNas_tur
      !
      do idime=1,ndime
        elmat_aux = 0.0_rp
        !
        ! stabilization parameter derivatives w.r.t velocity
        !
        call vecnor(gpvel(1,igaus),ndime,gpnve,2_ip)
        if (gpnve /= 0.0_rp) then
          do inode = 1,pnode
            dtau_dvel(inode) = -2.0_rp*gpden(igaus)*tau*tau*gpsha(inode,igaus)*gpvel(idime,igaus)/(gpnve*chale(1)) &
                               -4.0_rp*dgpmut_dvel(idime,inode,igaus)*tau*tau/(chale(2)*chale(2)*diff_param_tur(iunkn_tur))
          enddo
        else
          dtau_dvel(:) = 0.0_rp
        endif
        !
        ! calculus of perturbation * d(residual_st)/du
        !
        do inode=1,pnode
          gpad1(inode) =  0.0_rp
          do jdime=1,ndime
            gpad1(inode) = gpad1(inode) + gpvel(jdime,igaus) * gpcar(jdime,inode,igaus)
          end do
          gpadv(inode)=gpad1(inode)*gpden(igaus)
          resid(inode) = gpden(igaus)*gpsha(inode,igaus)*grtur(idime) + dsou_dvel(idime,inode,igaus)
          gpper(inode)  = ( gpsha(inode,igaus) + tau*gpadv(inode) ) * gpvol(igaus)
        end do
        !
        ! Calculus of residual_st AND d(perturbation)/d(vel)
        !
        resid_st = gprea(igaus)*gptur(iunkn_tur,1,igaus) - gprhs(igaus)
        do jdime=1,ndime
          resid_st = resid_st + gpden(igaus)*gpvel(jdime,igaus)*grtur(jdime)
        end do
        do inode = 1,pnode
          gpad1(inode) = 0.0_rp
            do jdime = 1,ndime
              gpad1(inode) = gpad1(inode) + gpvel(jdime,igaus) * gpcar(jdime,inode,igaus)
            end do
          gpadv(inode) = gpad1(inode) * gpden(igaus)
          do jnode=1,pnode
            dgpper_du(inode,jnode)  = (dtau_dvel(jnode)*gpadv(inode) + tau*gpden(igaus)*gpsha(jnode,igaus)*gpcar(idime,inode,igaus))*gpvol(igaus)
          enddo
        end do !inode
        !
        ! Diffusion Term
        !
!         do inode = 1,pnode
!           do jnode = 1,pnode
!             fact2 = gpvol(igaus) * dgpmut_dvel(idime,jnode,igaus)/diff_param_tur(iunkn_tur)  ! turbulent visc derivatives
!             xmuit = 0.0_rp
!             do jdime = 1,ndime
!               xmuit = xmuit + grtur(jdime) * gpcar(jdime,inode,igaus)
!             end do
!             elmat_aux(inode,jnode) = elmat_aux(inode,jnode) + xmuit * fact2
!           end do
!         end do
        !
        ! Assembly of the matrix for ASGS, SUPG or Full OSS
        !
        if( kfl_ortho /= 2 ) then
          do inode = 1,pnode
            do jnode = 1,pnode
              elmat_aux(inode,jnode) = elmat_aux(inode,jnode) + resid(jnode) * gpper(inode) + dgpper_du(inode,jnode) * resid_st
            end do
          end do
        end if
        !
        ! Apply B. C. 
        !
!         do inode = 1,pnode
!             ipoin = lnods(inode,ielem)
!             if(  kfl_fixno_tur(1,ipoin,iunkn_tur) > 0 ) elmat_aux(inode,:) = 0.0_rp
!         end do
        !
        ! production of the Transpose[dRk/du]*[Lambda_k]
        ! production of the Transpose[dRw/du]*[Lambda_w]
        !
        do inode = 1, pnode
          do jnode = 1, pnode
            RhsadjNas_tur(ielem)%a(idime, inode) = RhsadjNas_tur(ielem)%a(idime, inode) + elmat_aux(jnode,inode)*eltur(iunkn_tur,jnode,3)
          enddo
        enddo
	
      enddo ! idime        
    endif ! coupling 
       
    
  end do ! igaus

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!           elrhs(k)   <------- Rhsadjtur_tur(2)                   sent from omega         !!!!!!!
  !!!!!!!           elrhs(w)   <------- Rhsadjtur_tur(1)                   sent from kin           !!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   do inode = 1,pnode
!     if (iunkn_tur == 1) then
!       elrhs(inode) = elrhs(inode) - Rhsadjtur_tur(ielem)%a(2,inode)
!     else
!       elrhs(inode) = elrhs(inode) - Rhsadjtur_tur(ielem)%a(1,inode)
!     endif
!   enddo
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!           elrhs(k)     <-------   RhsadjTur_nsi(1)             sent from nastin       !!!!!!!
  !!!!!!!           elrhs(w)     <-------   RhsadjTur_nsi(2)             sent from nastin       !!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  if (kfl_coupl(ID_TURBUL,ID_NASTIN) >= 1 ) then
    do inode = 1,pnode
      elrhs(inode) = elrhs(inode) - RhsadjTur_nsi(ielem)%a(iunkn_tur,inode)
    enddo
  endif  

  
end subroutine tur_elmmsu_der_all

