!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_elmmul_der_all(&
     pnode     ,pgaus    ,gpvel ,gpdif    ,gprea, &
     gptem, gpgrd    ,gpden ,gpsha    ,gpcar,gpvol,&
     elmat,elrhs, chale    ,gphes ,elunk    ,&
     gpsph, ielem,lnods,gpdde,dsou_dtem,dsou_dcon,eltem_forw)
  !-----------------------------------------------------------------------
  !****f* Temper/tem_elmmul_der_all
  ! NAME
  !   tem_elmmul_der_all
  ! DESCRIPTION
  !   This routine is an alternative to elmadr routine, 
  ! OUTPUT
  !    ELMAT ... LHS matrix for current Gauss point
  !    ELRHS ... RHS vector for current Gauss point
  ! USES
  ! USED BY
  !    tem_elmope
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,ntens,mnode
  use def_temper, only       :  kfl_regim_tem,kfl_taust_tem,kfl_sgsti_tem,kfl_shock_tem,shock_tem,& 
                                staco_tem,kfl_ortho_tem
  use mod_tauadr, only       :  tauadr
  
  implicit none 
  integer(ip), intent(in)    :: pnode,pgaus, ielem
  real(rp),    intent(in)    :: gpvel(ndime,pgaus), gptem(pgaus,3)
  real(rp),    intent(in)    :: gpdif(pgaus),gprea(pgaus)
  real(rp),    intent(in)    :: gpden(pgaus)
  real(rp),    intent(in)    :: gpdde(pnode,pgaus)                    ! Density derivatives w.r.t nodal temperature
  real(rp),    intent(in)    :: gpgrd(ndime,pgaus)                    ! Density gradients
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus), chale(2)
  real(rp),    intent(in)    :: gphes(ntens,mnode,pgaus)              ! dNk/dxidxj
  real(rp),    intent(in)    :: elunk(pnode)
  real(rp),    intent(in)    :: gpsph(pgaus)
  real(rp),    intent(inout) :: elmat(pnode,pnode)
  real(rp),    intent(out)   :: elrhs(pnode)                          ! nothing from ADR_element_assembly
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(in)    :: dsou_dtem(pnode,pgaus)
  real(rp),    intent(in)    :: dsou_dcon(nspec,pnode,pgaus)
  real(rp),    intent(inout) :: eltem_forw(pnode)                     ! tempe_forw at elemental nodes
  
  integer(ip)                :: inode,jnode,idime,igaus,jdime,ipoin,ispec
  real(rp)                   :: resid(pnode, pgaus), gpper(pnode)
  real(rp)                   :: gpadv, gpnve, grvgr, resi2(pnode), sreac(pgaus)
  real(rp)                   :: rhnve, gpad1(pnode), gplap
  real(rp)                   :: gpsta(pgaus,2), gpstt
  real(rp)                   :: dresid_du(pnode, pgaus),dresid_dT(pnode, pgaus)
  real(rp)                   :: gpgrt(ndime,pgaus) 
  real(rp)                   :: gpstader(pnode),resid_st(pgaus),dgpper_du(pnode,pnode),dgpper_dT(pnode,pnode)
!  real(rp)                   :: elmder(pnode,pgaus)              ! elemental derivatives of the variables for the adjoint
  real(rp)                   :: elmat_aux(pnode,pnode)           ! auxiliary matrix for local derivatives
  real(rp)                   :: shock,staco(3)
  integer(ip)                :: kfl_taust,kfl_ortho,kfl_shock,kfl_sgsti
  
  !
  ! Initialization
  !
  kfl_ortho = kfl_ortho_tem 
  kfl_shock = kfl_shock_tem
  shock = shock_tem
  kfl_taust = kfl_taust_tem
  kfl_sgsti = kfl_sgsti_tem
  staco = staco_tem
  

  elrhs = 0.0_rp       ! elrhs starts from zero for adjoint but elmat is modified
  gpstader = 0.0_rp


  if (kfl_ortho ==0 ) then ! ASGS
     do igaus = 1, pgaus
        sreac (igaus) =  gprea (igaus)
     end do
  else                     ! SUPG, OSS
     do igaus = 1, pgaus
        sreac (igaus) = 0.0_rp
     end do
  end if
  ! 
  ! calculate eltem_forw
  !
  do inode=1,pnode
    ipoin=lnods(inode)
    eltem_forw(inode) = tempe_forw(ipoin,1)
  end do
  !
  !  calculate gpgrt using tempe_forw 
  !
  if (kfl_coupl(ID_TEMPER,ID_NASTIN) >= 1) then
    ! calculate gpgrt    
    do igaus = 1,pgaus
      do idime = 1,ndime
          gpgrt(idime,igaus) = 0.0_rp
	  do inode = 1,pnode
	    gpgrt(idime,igaus) = gpgrt(idime,igaus) + gpcar(idime,inode,igaus) * eltem_forw(inode)
	  end do
      end do
    end do
  endif
  
  !
  ! stabilization parameter
  !
  do igaus=1,pgaus
    call vecnor(gpvel(1,igaus),ndime,gpnve,2_ip)      
    rhnve = gpden(igaus)*gpnve
    call tauadr(&
	  kfl_taust, staco,rhnve,gpdif(igaus),sreac(igaus),&
	  chale(1),chale(2),gpsta(igaus,2))      
    gpsta(igaus,1) = gpsta(igaus,2)
    gpstt = 1.0_rp    
  enddo  
    
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!        ELMAT   ---->   dR_t/dT                  !!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  !
  ! Coupling with CHEMIC: term d( SUM[hk*Wk] )/dT 
  !  
  if (kfl_coupl(ID_TEMPER,ID_CHEMIC) >= 1) then
     
     do igaus=1,pgaus
        !
        ! calculus of residual and perturbation function gppre
        !        
        do inode=1,pnode
!            resid(inode, igaus) = gpden(igaus)*gpsha(inode,igaus)*dtinv*pabdf(1)
           gpad1(inode) =  0.0_rp
           grvgr =  0.0_rp
           gplap =  0.0_rp
           do idime=1,ndime
              gpad1(inode) = gpad1(inode) + gpvel(idime,igaus)*gpcar(idime,inode,igaus)           
              grvgr = grvgr + gpgrd(idime,igaus)*gpcar(idime,inode,igaus)       
              gplap = gplap + gphes(idime,inode,igaus)
           end do
           gpadv=gpad1(inode)*gpden(igaus)
!            resid (inode, igaus) = resid(inode, igaus) + gpadv + gprea(igaus)*gpsha(inode,igaus) - grvgr -gplap &
!                                                       + elmder(inode,igaus)*gpsha(inode,igaus)
           resid (inode, igaus) = dsou_dtem(inode,igaus)*gpsha(inode,igaus)
           resi2 (inode) = 0.0_rp !- grvgr -gplap
           gpper(inode)  = ( gpsha(inode,igaus)*(gpstt-gpsta(igaus,1)*sreac(igaus))  + gpsta(igaus,1)*gpadv ) * gpvol(igaus)        
        end do
        !
        ! Assembly of the matrix and rhs
        !        
        do inode =1, pnode
           do jnode =1, pnode
              elmat(inode, jnode) = elmat(inode, jnode) + resid(jnode, igaus)*gpper(inode) &        
                   -resi2(jnode)*gpsha(inode,igaus)*gpvol(igaus)
           end do
        end do
        
     end do !igaus
  end if
    
  !
  ! for low mach: the effecto of d(rho)/dT and d(tau)/dT
  !
  if ( kfl_regim_tem==3 ) then 
               
    do igaus=1,pgaus
      !
      ! stabilization parameter derivatives w.r.t. temperature
      !
      call vecnor(gpvel(1,igaus),ndime,gpnve,2_ip)
      if (gpnve /= 0.0_rp) then
	do inode = 1,pnode
	    gpstader(inode) = - 2.0_rp*gpnve*gpsta(igaus,1)**2/chale(1) * gpsph(igaus) * gpdde(inode,igaus)
	enddo
      endif
      !
      ! calculus of perturbation * d(residual)/dT = perturbation * d(rho)/dT*Cp*ui*dT/dxi
      !  
      do inode=1,pnode
	gpad1(inode) =  0.0_rp
	dresid_dT(inode, igaus) = 0.0_rp
	do jdime=1,ndime
	  gpad1(inode) = gpad1(inode) + gpvel(jdime,igaus)*gpcar(jdime,inode,igaus)           
	  dresid_dT(inode, igaus) = dresid_dT(inode, igaus) + gpsph(igaus)*gpdde(inode,igaus)*gpvel(jdime,igaus)*gpgrt(jdime,igaus)
	end do
	gpadv=gpad1(inode)*gpden(igaus)
	gpper(inode)  = ( gpsha(inode,igaus)*(gpstt-gpsta(igaus,1)*sreac(igaus))  + gpsta(igaus,1)*gpadv ) * gpvol(igaus)
      end do
      !
      ! calculus elmat = perturbation * d(residual)/dT
      !
      do inode =1, pnode
	do jnode =1, pnode
	    elmat(inode, jnode) = elmat(inode, jnode) + dresid_dT(jnode, igaus)*gpper(inode)
	end do
      end do
      !
      ! calculus of d(perturbation)/dT * residual_st = ( Cp*d(rho)/dT*tau*ui*dN/dxi + Cp*rho*d(tau)/dT*ui*dN/dxi ) * residual_st
      !
      
      ! residual_st term  --->   rho*ui*dT/dxi
      resid_st(igaus) = 0.0_rp 
      do jdime=1,ndime
	  resid_st(igaus) = resid_st(igaus) + gpden(igaus)*gpvel(jdime,igaus)*gpgrt(jdime,igaus)          
      end do
      
      ! d(perturbation)/dT term ---> ( Cp*rho*d(tau)/dT*ui*dN/dxi + Cp*d(rho)/dT*tau*ui*dN/dxi )
      do inode=1,pnode
	gpad1(inode) =  0.0_rp
	do jdime=1,ndime
	  gpad1(inode) = gpad1(inode) + gpvel(jdime,igaus)*gpcar(jdime,inode,igaus)           
	end do
	gpadv=gpad1(inode)*gpden(igaus)
	do jnode=1,pnode   
	  dgpper_dT(inode,jnode)  = ( gpstader(jnode)*gpadv + gpsta(igaus,1)*gpsph(igaus)*gpdde(jnode,igaus)*gpad1(inode) ) * gpvol(igaus)        
	enddo
      end do
      !
      ! calculus elmat = d(perturbation)/du * residual
      !       
      do inode =1, pnode
	do jnode =1, pnode
	    elmat(inode, jnode) = elmat(inode, jnode) + dgpper_dT(inode,jnode)*resid_st(igaus)
	end do
      end do

    end do !gauss  
                  
  endif
      
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!  RhsadjChe_tem   <------  Transpose[dRt/dYs] [Lambda_t]     to be sent to chemic    !!!!!!!
  !!!!!!!           elrhs   <------- RhsadjTem_chm                     sent from chemic        !!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Coupling with CHEMIC 
  !    
  if (kfl_coupl(ID_TEMPER,ID_CHEMIC) >= 1 ) then
    
    ! put RhsadjChe_tem to zero
    RhsadjChe_tem(ielem)%a = 0.0_rp
    
    !
    ! loop over the species to calculate RhsadjChe_tem
    !
    do ispec = 1, nspec 
      
      ! put zero for iteration on ispec
      elmat_aux = 0.0_rp
        
      do igaus=1,pgaus
        !
        ! calculus of residual and perturbation function gppre
        !        
        do inode=1,pnode
!            resid(inode, igaus) = gpden(igaus)*gpsha(inode,igaus)*dtinv*pabdf(1)
           gpad1(inode) =  0.0_rp
           grvgr =  0.0_rp
           gplap =  0.0_rp
           do idime=1,ndime
              gpad1(inode) = gpad1(inode) + gpvel(idime,igaus)*gpcar(idime,inode,igaus)           
              grvgr = grvgr + gpgrd(idime,igaus)*gpcar(idime,inode,igaus)       
              gplap = gplap + gphes(idime,inode,igaus)
           end do
           gpadv=gpad1(inode)*gpden(igaus)
!            resid (inode, igaus) = resid(inode, igaus) + gpadv + gprea(igaus)*gpsha(inode,igaus) - grvgr -gplap &
!                                                       + elmder(inode,igaus)*gpsha(inode,igaus)
           resid (inode, igaus) = dsou_dcon(ispec,inode,igaus)*gpsha(inode,igaus)
           resi2 (inode) = 0.0_rp !- grvgr -gplap
           gpper(inode)  = ( gpsha(inode,igaus)*(gpstt-gpsta(igaus,1)*sreac(igaus))  + gpsta(igaus,1)*gpadv ) * gpvol(igaus)        
        end do
        !
        ! Assembly of the matrix and rhs
        !        
        do inode =1, pnode
           do jnode =1, pnode
              elmat_aux(inode, jnode) = elmat_aux(inode, jnode) + resid(jnode, igaus)*gpper(inode) &        
                   -resi2(jnode)*gpsha(inode,igaus)*gpvol(igaus)
           end do
        end do

      end do ! igaus
      !
      ! production of the Transpose[dRt/dYs] [Lambda_t] to be sent to chemic
      !
      do inode = 1, pnode
	  do jnode = 1, pnode
	    RhsadjChe_tem(ielem)%a(inode, ispec) = RhsadjChe_tem(ielem)%a(inode, ispec) + elmat_aux(jnode,inode)*elunk(jnode)
	  enddo
      enddo
           
    enddo ! ispec

    !
    ! add RhsadjTem_chm to elrhs 
    !
    do inode = 1,pnode
      elrhs(inode) = elrhs(inode) - RhsadjTem_chm(ielem)%a(inode)
    enddo
    
  endif
    
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!  RhsadjNas_tem   <-------  Transpose[dRt/du] [Lambda_t]    to be sent to nastin   !!!!!!!
  !!!!!!!           elrhs      <-------   RhsadjTem_nsi             sent from nastin         !!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Coupling with nastin
  !
  if (kfl_coupl(ID_TEMPER,ID_NASTIN) >= 1 ) then
     
     ! put RhsadjNas_tem to zero
     RhsadjNas_tem(ielem)%a = 0.0_rp
     
     !
     ! loop over ndime to calculate RhsadjNas_tem
     !
     do idime=1,ndime
     
       ! put zero for iteration on idime
       elmat_aux = 0.0_rp
     
       do igaus=1,pgaus
	  !
	  ! stabilization parameter derivatives
	  !
	  call vecnor(gpvel(1,igaus),ndime,gpnve,2_ip)
	  if (gpnve /= 0.0) then
	    do inode = 1,pnode
               gpstader(inode) = -2.0_rp*gpden(igaus)*gpsta(igaus,1)*gpsta(igaus,1)*gpsha(inode,igaus)*gpvel(idime,igaus)/(gpnve*chale(1))
	    enddo
	  endif
	  !
	  ! calculus of perturbation * d(residual)/du 
	  !      
	  do inode=1,pnode    
	    gpad1(inode) =  0.0_rp
	    do jdime=1,ndime
	      gpad1(inode) = gpad1(inode) + gpvel(jdime,igaus)*gpcar(jdime,inode,igaus)           
	    end do
	    gpadv=gpad1(inode)*gpden(igaus)
	    dresid_du(inode, igaus) = gpden(igaus)*gpsha(inode,igaus)*gpgrt(idime,igaus)
	    gpper(inode)  = ( gpsha(inode,igaus)*(gpstt-gpsta(igaus,1)*sreac(igaus))  + gpsta(igaus,1)*gpadv ) * gpvol(igaus)        
	  end do
	  !
	  ! calculus elmat_aux = perturbation * d(residual)/du
	  !       
	  do inode =1, pnode
	    do jnode =1, pnode
		elmat_aux(inode, jnode) = elmat_aux(inode, jnode) + dresid_du(jnode, igaus)*gpper(inode)
	    end do
	  end do
	  !
	  ! calculus of d(perturbation)/du * residual
	  !
	  ! strong residual
	  resid_st(igaus) = 0.0_rp ! rho*ui*dT/dxi
	  do jdime=1,ndime
	      resid_st(igaus) = resid_st(igaus) + gpden(igaus)*gpvel(jdime,igaus)*gpgrt(jdime,igaus)          
	  end do
	  do inode=1,pnode
	    gpad1(inode) =  0.0_rp
	    do jdime=1,ndime
	      gpad1(inode) = gpad1(inode) + gpvel(jdime,igaus)*gpcar(jdime,inode,igaus)           
	    end do
	    gpadv=gpad1(inode)*gpden(igaus)
	    do jnode=1,pnode   
	      dgpper_du(inode,jnode)  = ( gpstader(jnode)*gpadv + gpsta(igaus,1)*gpden(igaus)*gpsha(jnode,igaus)*gpcar(idime,inode,igaus) ) * gpvol(igaus)        
	    enddo
	  end do
	  !
	  ! calculus elmat_aux = d(perturbation)/du * residual
	  !       
	  do inode =1, pnode
	    do jnode =1, pnode
		elmat_aux(inode, jnode) = elmat_aux(inode, jnode) + dgpper_du(inode,jnode)*resid_st(igaus)
	    end do
	  end do

       end do !gauss
       !
       ! production of the Transpose[dRt/du]*[Lambda_t]
       !
       do inode = 1, pnode
 	  do jnode = 1, pnode
	    RhsadjNas_tem(ielem)%a(idime, inode) = RhsadjNas_tem(ielem)%a(idime, inode) + elmat_aux(jnode,inode)*elunk(jnode)
	  enddo
       enddo
              
    enddo !idime
    
     !
     ! add RhsadjTem_nsi to elrhs 
     !
     do inode = 1,pnode
       elrhs(inode) = elrhs(inode) - RhsadjTem_nsi(ielem)%a(inode)
     enddo


  endif

end subroutine tem_elmmul_der_all

