!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_resdiff_all(&
			  pnode,pgaus,gpden,gpsha,gpdif,gpcar,lnods,ielem,gprea, &
			  chale,gpvel,gpcod,gpvol,gptem,gprhs,gphes,eltem_forw)
  !------------------------------------------------------------------------
  !****f* Temper/tem_resdiff_all
  ! NAME 
  !    tem_resdiff_all
  ! DESCRIPTION
  !    
  !      This calculates the partial derivative of the residual w.r.t. design variables (dR/dD)
  !      Depondency to the design variables and couplement 
  ! USES
  ! USED BY
  !    tem_elmope
  !------------------------------------------------------------------------

  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,ntens,mnode
  use def_kermod, only       :  kfl_dvar_type
  use def_temper, only       :  kfl_taust_tem,kfl_sgsti_tem,kfl_shock_tem,shock_tem,& 
				  staco_tem,kfl_ortho_tem,resdiff_tem
  use def_master
  use def_kermod, only       : kfl_ndvars_opt
!  use mod_matrix
  use mod_assresdiff

  
  implicit none

  integer(ip), intent(in)    :: pnode,pgaus
  
  integer(ip),    intent(in)        :: lnods(pnode),ielem
  real(rp),       intent(in)        :: gpsha(pnode,pgaus)
  real(rp),       intent(in)        :: gpden(pgaus),gpdif(pgaus)
  real(rp),       intent(in)        :: gpcar(ndime,pnode),gprea(pgaus)
  real(rp),       intent(in)        :: chale(2),gpvel(ndime,pgaus),gpcod(ndime,pgaus)
  real(rp),       intent(in)        :: gpvol(pgaus), gptem(pgaus)
  real(rp),       intent(inout)     :: gprhs(pgaus)
  real(rp),       intent(in)        :: gphes(ntens,mnode,pgaus)              ! dNk/dxidxj
  real(rp),       intent(in)        :: eltem_forw(pnode)
  
  integer(ip) :: idesvar
  real(rp)    :: elresdiff(pnode),chemicalHeatmassdiff_des(pgaus),& 
                       chemicalHeatmassdiff(pgaus, kfl_ndvars_opt)
    
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!        design variables will be thermal conductivity       !!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (kfl_dvar_type == 1) then 
    
      call tem_elmmuldiff_cond(&
	      kfl_ortho_tem,kfl_shock_tem,shock_tem,kfl_taust_tem,kfl_sgsti_tem, &
	      staco_tem ,pnode     ,pgaus ,gpvel, gpdif,gprea       ,&
	      gptem     ,gprhs ,gpden, gpsha, gpcar     ,gpvol,&
	      elresdiff     ,chale ,gphes, &
	      ielem,lnods,eltem_forw)

      !
      ! Effect of Dirichlet boundary conditions (elresdiff = 0.0)
      !
      if( solve(1) % kfl_iffix == 0 ) then
	  call tem_elmdirdiff(pnode,lnods,elresdiff)
      end if
      !
      ! Assembly
      !
      call assresdiff(solve(1) % ndofn, 1_ip, 1_ip, pnode, lnods, elresdiff, resdiff_tem)
  
  endif
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!        design variables will be Af, Ab             !!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (kfl_dvar_type == 3) then 
      !
      ! read chemicalmassdiff = d( SUM[hk*Wk )/d(des_var) sent from chemic 
      !
      call tem_chemic_diff(1_ip,ielem, pgaus, pnode, chemicalHeatmassdiff)
      
      do idesvar = 1,kfl_ndvars_opt
    
	  chemicalHeatmassdiff_des = chemicalHeatmassdiff(:,idesvar) 
	  !
	  ! calculate elresdiff
	  !
	  call tem_elmmuldiff(&
		  ielem,idesvar, kfl_taust_tem, kfl_ortho_tem,&
		  staco_tem ,pnode     ,pgaus ,gpvel, gpdif,gprea,&
		  gpden, gpcod, gpsha, gpcar ,gpvol, chemicalHeatmassdiff_des,&
		  elresdiff ,chale)
	  !
	  ! Effect of Dirichlet boundary conditions (elresdiff = 0.0)
	  !
	  if( solve(1) % kfl_iffix == 0 ) then
	      call tem_elmdirdiff(pnode,lnods,elresdiff)
	  end if
	  !
	  ! Assembly
	  !
	  call assresdiff(solve(1) % ndofn,1_ip,idesvar,pnode,lnods,elresdiff,resdiff_tem)
      end do 
  endif

end subroutine tem_resdiff_all


subroutine tem_elmmuldiff(&
     ielem,idesvar, kfl_taust, kfl_ortho, &
     staco, pnode ,pgaus ,gpvel ,gpdif ,gprea, &
     gpden ,gpcod, gpsha ,gpcar,gpvol,chemicalHeatmassdiff_des,&
     elresdiff ,chale)
  !-----------------------------------------------------------------------
  !****f* Temper/tem_elmmuldiff
  ! NAME
  !   tem_elmmuldiff
  ! DESCRIPTION
  !   This routine calculates the derivative of the residual w.r.t Af and Ab preexponencial factors 
  !
  ! USES
  ! USED BY
  !    tem_elmmuldiff
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,mnode
  use mod_tauadr, only       :  tauadr

!   use def_temper, only       :  kfl_regim_tem
  implicit none
  integer(ip), intent(in)    :: pnode,pgaus,idesvar,kfl_ortho,ielem
  integer(ip), intent(in)    :: kfl_taust
  real(rp),    intent(in)    :: gpvel(ndime,pgaus)
  real(rp),    intent(in)    :: gpdif(pgaus),gprea(pgaus),staco(3)
  real(rp),    intent(in)    :: gpden(pgaus),gpcod(ndime,pgaus)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus), chale(2)
  real(rp),    intent(in)    :: chemicalHeatmassdiff_des(pgaus)
  real(rp),    intent(out)   :: elresdiff(pnode) 
  integer(ip)                :: inode,jnode,idime,igaus
  real(rp)                   :: gpper(pnode)
  real(rp)                   :: gpadv, gpnve, sreac(pgaus)
  real(rp)                   :: rhnve, gpad1(pnode)
  real(rp)                   :: gpsta(2), gpstt
  real(rp)                   :: elmatdiff(pnode,pnode),elrhsdiff(pnode),gprhsdiff(pgaus)
  !
  ! Initialization
  !

  do inode=1,pnode
     elrhsdiff(inode)=0.0_rp
     elresdiff(inode)=0.0_rp
     do jnode=1,pnode
        elmatdiff(jnode,inode)=0.0_rp
     end do
  end do

  if (kfl_ortho ==0 ) then ! ASGS
    do igaus = 1, pgaus
	sreac (igaus) =  gprea (igaus)
    end do
  else                     ! SUPG, OSS
    do igaus = 1, pgaus
	sreac (igaus) = 0.0_rp
    end do
  end if
  
  do igaus = 1,pgaus
      !
      ! stabilization parameter
      !
      call vecnor(gpvel(1,igaus),ndime,gpnve,2_ip)      
      rhnve = gpden(igaus)*gpnve
      call tauadr(&
	    kfl_taust, staco,rhnve,gpdif(igaus),sreac(igaus),&
	    chale(1),chale(2),gpsta(2))      
      gpsta(1) = gpsta(2)
      gpstt = 1.0_rp    
      !
      ! Calculate gprhsdiff
      !
      gprhsdiff(igaus) =  chemicalHeatmassdiff_des(igaus)
      !
      ! Calculate elrhsdiff
      !
      do inode=1,pnode
	  gpad1(inode) =  0.0_rp
	  do idime=1,ndime
	    gpad1(inode) = gpad1(inode) + gpvel(idime,igaus)*gpcar(idime,inode,igaus)           
	  end do
	  gpadv=gpad1(inode)*gpden(igaus)
	  gpper(inode)=(gpsha(inode,igaus)*(gpstt-gpsta(1)*sreac(igaus))+gpsta(1)*gpadv)*&	
			gpvol(igaus)        
      end do
      
      do inode =1, pnode
	elrhsdiff(inode) = elrhsdiff(inode) + gprhsdiff(igaus)*gpper(inode)
      end do
      
  end do
  !
  ! Calculate elresdiff (No contribution from elmatdiff until now)
  !
  do inode =1, pnode
    elresdiff(inode) = elrhsdiff(inode) !+ elmatdiff terms
  end do
  
end subroutine tem_elmmuldiff

subroutine tem_elmdirdiff(&
                           pnode,lnods,elresdiff)
  !------------------------------------------------------------------------
  !****f* Temper/tem_elmdirdiff
  ! NAME 
  !    tem_elmdirdiff
  ! DESCRIPTION
  !    This routine prescribes the Dirichlet conditions to the elresdiff
  ! USES
  ! USED BY
  !    tem_elmdirdiff
  !    
  !------------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_temper, only       :  kfl_fixno_tem
  use def_elmtyp, only       :  ELHAN
  implicit none
  integer(ip), intent(in)    :: pnode
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(inout) :: elresdiff(pnode)
  integer(ip)                :: inode,ipoin
  
  do inode = 1,pnode
    ipoin = lnods(inode)
    if(  kfl_fixno_tem(1,ipoin) == 1 .or. kfl_fixno_tem(1,ipoin) == 4 .or. kfl_fixno_tem(1,ipoin) == 5 ) then
	elresdiff(inode) = 0.0_rp
    end if
  end do

end subroutine tem_elmdirdiff

subroutine tem_elmmuldiff_cond(&
		    kfl_ortho, kfl_shock,shock ,kfl_taust,kfl_sgsti, &
		    staco, pnode     ,pgaus    ,gpvel ,gpdif    ,gprea, &
		    gptem,gprhs    ,gpden ,gpsha    ,gpcar,gpvol,&
		    elresdiff     ,chale    ,gphes, &
		    ielem,lnods,eltem_forw)
  !-----------------------------------------------------------------------
  !****f* Temper/tem_elmmuldiff_cond
  ! NAME
  !   ten_elmmul
  ! DESCRIPTION
  !   This routine calculates the derivative of the residual w.r.t thermal conductivity
  ! USES
  ! USED BY
  !    tem_elmmuldiff_cond
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,ntens,mnode
  use mod_tauadr, only       :  tauadr

  implicit none
  integer(ip), intent(in)    :: pnode,pgaus, kfl_sgsti, ielem
  integer(ip), intent(in)    :: lnods(pnode)
  integer(ip), intent(in)    :: kfl_ortho, kfl_shock, kfl_taust
  real(rp),    intent(in)    :: gpvel(ndime,pgaus), gptem(pgaus,3)
  real(rp),    intent(in)    :: gpdif(pgaus),gprea(pgaus),staco(3)
  real(rp),    intent(inout) :: gprhs(pgaus)
  real(rp),    intent(in)    :: gpden(pgaus)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus), chale(2)
  real(rp),    intent(in)    :: gphes(ntens,mnode,pgaus)              ! dNk/dxidxj
  real(rp),    intent(in)    :: shock
  real(rp),    intent(in)    :: eltem_forw(pnode)
  real(rp),    intent(out)   :: elresdiff(pnode)
  integer(ip)                :: inode,idime,igaus
  real(rp)                   :: gpper(pnode)
  real(rp)                   :: gpadv, gpnve, sreac(pgaus)
  real(rp)                   :: rhnve, gpad1(pnode)
  real(rp)                   :: umbra, CD, SD, factt            ! Shock capturring variables
  real(rp)                   :: gpsta(2), gpstt
  real(rp)                   :: gpstadiff(pgaus),gpgte_forw(pgaus,ndime)
  
  !
  ! Initialization
  !
  factt = 0.75_rp
  umbra = 1.0e-6_rp
  CD=0.0_rp
  SD=0.0_rp
  do inode=1,pnode
     elresdiff(inode)=0.0_rp
  end do
  
  if (kfl_ortho ==0 ) then ! ASGS
     do igaus = 1, pgaus
        sreac (igaus) =  gprea (igaus)
     end do
  else                     ! SUPG, OSS
     do igaus = 1, pgaus
        sreac (igaus) = 0.0_rp
     end do
  end if
  
  do igaus = 1, pgaus
    do idime = 1, ndime
      gpgte_forw(igaus,idime) = 0.0_rp
    enddo
  enddo
  
  !
  ! Calculate some variables
  !
  do igaus = 1, pgaus
    do idime=1,ndime
      do inode=1,pnode
	gpgte_forw(igaus,idime) = gpgte_forw(igaus,idime) + gpcar(idime,inode,igaus)*eltem_forw(inode)
      enddo
    enddo
  enddo
  
   
  do igaus=1,pgaus
    !
    ! stabilization parameter
    !
    call vecnor(gpvel(1,igaus),ndime,gpnve,2_ip)      
    rhnve = gpden(igaus)*gpnve
    call tauadr(&
	  kfl_taust, staco,rhnve,gpdif(igaus),sreac(igaus),&
	  chale(1),chale(2),gpsta(2))      
    gpsta(1) = gpsta(2)
    gpstt = 1.0_rp    
    !
    !  gpstadiff
    !
    gpstadiff(igaus) = -4.0_rp*gpsta(1)*gpsta(1)/(chale(2)*chale(2))      
    !
    ! calculus of residual and perturbation function gppre
    !
    do inode=1,pnode
	gpad1(inode) =  0.0_rp
	do idime=1,ndime
	  gpad1(inode) = gpad1(inode) + gpvel(idime,igaus)*gpcar(idime,inode,igaus)           
	end do
	gpadv=gpad1(inode)*gpden(igaus)
	gpper(inode)  = ( gpstadiff(igaus)*gpadv ) * gpvol(igaus)        
    end do
    
    do idime=1,ndime
      gprhs(igaus) = gprhs(igaus) + gpden(igaus)*gpvel(idime,igaus)*gpgte_forw(igaus,idime)
    end do
    !
    ! Assembly of the matrix and rhs
    !  
    
    ! term related to dN/dxi * dT/dxi
    do inode = 1, pnode
      do idime= 1, ndime
	elresdiff(inode) = elresdiff(inode) + gpgte_forw(igaus,idime)*gpcar(idime, inode,igaus)*gpvol(igaus)
      end do
    enddo  
    
    !       term related to d(tau)/dk * vi * dN/dxi * ( T/dt + vi * dT/dxi)
    do inode =1, pnode
	elresdiff(inode) = elresdiff(inode) + gprhs(igaus)*gpper(inode)
    end do
	    
  end do !igaus
  

  
end subroutine tem_elmmuldiff_cond

