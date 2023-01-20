!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Tempre
!> @{
!> @file    tem_elmpre.f90
!> @author  houzeaux
!> @date    2020-04-02
!> @brief   Compute some Gauss values
!> @details The variables are:
!>          GPGRD(NDIME) ... grad(k) coefficient
!>          GPTEM .......... Temperature of previous iterations and time steps
!>          GPVEL .......... Advection a
!>          GPADV .......... Advection term a.grad(Ni)
!>          GPRGD .......... Thermal conductivity gradient grad(k+kt)
!> @} 
!-----------------------------------------------------------------------

subroutine tem_elmpre(&
     ielem,pnode,pgaus,pmate,densi,gpden,gpsph,gpsgv,&
     gpsha,gpcar,gphes,elvel,eltem,elcod,elmsh,&
     gpvel,gptem,gprhs,gpcod,gpgrt,lnods,gpmsh,&
     gprea, gptke)
 
  use def_kintyp,                  only : ip,rp
  use def_master,                  only : dpthe,cutim,vesgs,velom,FUNCTION_DISCRETE, untur
!  use def_master,                  only : dpthe,cutim,vesgs,velom,FUNCTION_DISCRETE, kfl_paral, untur
  use def_domain,                  only : mnode,ndime,ntens, walld
  use def_domain,                  only : xfiel,canhe,canla,heiov
  use def_temper,                  only : kfl_advec_tem,ADR_tem
  use def_temper,                  only : kfl_sourc_tem,kfl_sgsve_tem
  use def_temper,                  only : kfl_regim_tem
  use def_temper,                  only : sourc_tem,kfl_rhs_scal_tem
  use def_temper,                  only : kfl_ellen_tem
  use def_temper,                  only : kfl_sonum_tem
  use def_temper,                  only : lsour_material_tem
  use def_temper,                  only : xsour_material_tem
  use def_temper,                  only : SOURCE_TERM_SPACE_TIME
  use def_temper,                  only : SOURCE_TERM_FIELD
  use def_temper,                  only : SOURCE_TERM_MATERIAL
  use def_temper,                  only : kfl_nudgi_tem, nudgi_tem, kfl_forest_tem
  use mod_ADR,                     only : BDF
  use def_kermod,                  only : gasco , kfl_prtur_abl_ker
  use mod_ker_space_time_function, only : ker_space_time_function
  use mod_ker_tendencies,          only : kfl_tendencies_ker, get_tendencies_th
  use mod_ker_functions         ,  only : ker_functions
  implicit none
  
  integer(ip), intent(in)    :: ielem
  integer(ip), intent(in)    :: pnode
  integer(ip), intent(in)    :: pgaus
  integer(ip), intent(in)    :: pmate
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(in)    :: gpsph(pgaus),gpsgv(ndime,pgaus)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gphes(ntens,mnode,pgaus)
  real(rp),    intent(in)    :: elvel(ndime,pnode)
  real(rp),    intent(in)    :: eltem(pnode,*)
  real(rp),    intent(in)    :: elcod(ndime,pnode)
  real(rp),    intent(in)    :: elmsh(ndime,pnode)
  real(rp),    intent(in)    :: densi(pgaus) ! density
  real(rp),    intent(out)   :: gpden(pgaus) ! coefficient multiplying the equation, can be rho,rhocp or rhocv,...  
  real(rp),    intent(out)   :: gpvel(ndime,pgaus),gpcod(ndime,pgaus)
  real(rp),    intent(out)   :: gprhs(pgaus)
  real(rp),    intent(out)   :: gptem(pgaus,*)
  real(rp),    intent(out)   :: gpgrt(ndime,pgaus)
  real(rp),    intent(out)   :: gpmsh(ndime,pgaus)
  real(rp),    intent(out)   :: gptke(pgaus)
  real(rp),    intent(inout) :: gprea(pgaus,1)
  integer(ip)                :: idime,inode,igaus,itime,mxdim
  real(rp)                   :: dummr, gphei, elhei(pnode), ellad(pnode), elhot(pnode), gptem_ref
  real(rp)                   :: gphot, qradf, gplad, eltke(pnode)
  !
  ! Temperature: GPTEM
  !
  if( ADR_tem % kfl_time_integration /= 0 ) then
     do igaus = 1,pgaus
        gptem(igaus,2) = 0.0_rp
        do inode = 1,pnode
           gptem(igaus,2) = gptem(igaus,2) + gpsha(inode,igaus) * eltem(inode,2)
        end do
     end do
     if( ADR_tem % kfl_time_scheme == BDF ) then
        do itime = 3,ADR_tem % kfl_time_order + 1
           do igaus = 1,pgaus
              gptem(igaus,itime) = 0.0_rp
              do inode = 1,pnode
                 gptem(igaus,itime) = gptem(igaus,itime) &
                      + eltem(inode,itime) * gpsha(inode,igaus)
              end do
           end do
        end do
     end if
  end if
  !
  ! Density: GPDEN=rho*Cp
  !
  if( kfl_regim_tem == 1 ) then !Este se puede comentar se usan 3y4
     do igaus = 1,pgaus
        gpden(igaus) = densi(igaus) * (gpsph(igaus)-gasco)     ! rho*Cv=rho*(Cp-R)
     end do
  else if (kfl_regim_tem /= 4) then
     do igaus = 1,pgaus
        gpden(igaus) = densi(igaus) * gpsph(igaus)             ! rho*Cp
     end do
  else
     gpden(1:pgaus) = densi(1:pgaus)
  end if
  
  if ( kfl_rhs_scal_tem > 0 )  gpden(1:pgaus) = 1.0_rp
  !
  ! Coordinates
  !
  do igaus = 1,pgaus
     do idime = 1,ndime
        gpcod(idime,igaus) = 0.0_rp
     end do
     do inode = 1,pnode
        do idime = 1,ndime
           gpcod(idime,igaus) = gpcod(idime,igaus)&
                + gpsha(inode,igaus) * elcod(idime,inode)
        end do
     end do
  end do
  !
  ! Velocity GPVEL=a and advection term GPADV=a.grad(Ni)-(k/r,0).grad(Ni)
  !
  do igaus = 1,pgaus
     do idime = 1,ndime
        gpvel(idime,igaus) = 0.0_rp
     end do
  end do

  if( kfl_advec_tem /= 0 ) then

     if( kfl_advec_tem == 1 ) then
        do igaus = 1,pgaus
           do inode = 1,pnode
              do idime = 1,ndime
                 gpvel(idime,igaus) = gpvel(idime,igaus) &
                      + gpsha(inode,igaus) * elvel(idime,inode)
              end do
           end do
        end do
        if( kfl_sgsve_tem == 1 ) then 
           do igaus = 1,pgaus
              do idime = 1,ndime
                 gpvel(idime,igaus) = gpvel(idime,igaus) +vesgs(ielem) % a(idime,igaus,1)
              end do
           end do
        end if

     else if( kfl_advec_tem >= 2 ) then
        call tem_velfun(pgaus,gpcod,gpvel)

     else if( kfl_advec_tem == -1 ) then
        do idime = 1,ndime
           gpgrt(idime,igaus) = 0.0_rp
        end do
        do inode = 1,pnode
           do idime = 1,ndime
              gpgrt(idime,igaus) = gpgrt(idime,igaus)&
                   + gpcar(idime,inode,igaus)*eltem(inode,1)
           end do
        end do
        dummr = 0.0_rp
        do idime = 1,ndime
           dummr = dummr +gpgrt(idime,igaus) * gpgrt(idime,igaus)
        end do
        dummr = sqrt(dummr)
        if( dummr /= 0.0_rp ) then
           do idime = 1,ndime
              gpvel(idime,igaus) = -gpgrt(idime,igaus) / dummr
           end do
        end if
     end if

     if( associated(velom) ) then
        gpmsh = 0.0_rp
        do igaus = 1,pgaus
           do inode = 1,pnode
              do idime = 1,ndime
                 gpmsh(idime,igaus) = gpmsh(idime,igaus)&
                      + gpsha(inode,igaus) * elmsh(idime,inode)
              end do
           end do
           gpvel(1:ndime,igaus) = gpvel(1:ndime,igaus) - gpmsh(1:ndime,igaus)
        end do        
     end if
  end if
  !
  ! Source term: GPRHS
  !
  if( kfl_sourc_tem == SOURCE_TERM_MATERIAL ) then
     !
     ! Material
     !
     if( lsour_material_tem(pmate) == 2 ) then
        do igaus = 1,pgaus
           gprhs(igaus) = xsour_material_tem(1,pmate) * sourc_tem
        end do
     end if
     
  else if( kfl_sourc_tem == SOURCE_TERM_SPACE_TIME ) then
     !
     ! Space time
     !
     mxdim = min(2_ip,ndime)
     do igaus = 1,pgaus
        call ker_space_time_function(&
             kfl_sonum_tem,gpcod(1,igaus),gpcod(mxdim,igaus),gpcod(ndime,igaus),cutim,gprhs(igaus))
        gprhs(igaus) = gprhs(igaus) * sourc_tem
     end do
     
  else if( kfl_sourc_tem == SOURCE_TERM_FIELD ) then
     !
     ! Field
     !
     do igaus = 1,pgaus
        gprhs(igaus) = xfiel(kfl_sonum_tem) % a(1,ielem,1) * sourc_tem
     end do    
 
  else
     gprhs(1:pgaus) = 0.0_rp     
  end if

  ! nudging term  + nudg(tempe - teref)
  if (kfl_nudgi_tem==1_ip) then
     if (.not.kfl_tendencies_ker)  call runend('tem_elmpre: nudging term only works with tendencies')
     ! reactive term ! left hand side
     gprea(1:pgaus,1) =  gprea(1:pgaus,1) + gpden(1:pgaus)*nudgi_tem
     ! height of gauss points to interpolate temperature
     ! cast elemet height from wall distante
     elhei(1:pnode)   =  walld(lnods(1:pnode)) 
     do igaus = 1,pgaus
        gphei =0.0_rp
        do inode =1, pnode
          gphei = gphei + gpsha(inode, igaus)*elhei(inode)          
        end do
        call get_tendencies_th(gphei, cutim, gptem_ref)       
        gprhs(igaus) =  gprhs(igaus) + gpden(igaus)*nudgi_tem*gptem_ref
     end do
  end if 
  if( kfl_forest_tem==1_ip) then ! thermal forest, adds radiative heat transfer 

!     interpolate qradf(t), radiative heat at canopy top     
     call ker_functions(1_ip, 2_ip, FUNCTION_DISCRETE, 1.0_rp, qradf)

!     if (kfl_paral.lt.10)  write (18, *) cutim, qradf
     elhei(1:pnode) = canhe( lnods(1:pnode))
     ellad(1:pnode) = canla( lnods(1:pnode))
     elhot(1:pnode) = heiov( lnods(1:pnode))
     do igaus =1, pgaus
        gphei =0.0_rp
        gplad =0.0_rp
        gphot =0.0_rp  ! height over terrain
        do inode =1, pnode 
           gphei =   gphei + gpsha(inode, igaus)*elhei(inode) ! top of canopy
           gplad =   gplad + gpsha(inode, igaus)*ellad(inode) ! lad
           gphot =   gphot + gpsha(inode, igaus)*elhot(inode) ! height over terrain
        end do
        ! div ofradiation heat 
        if (gphot.lt.gphei.and.gphot.gt.0.4_rp*gphei) then
!           gprhs(igaus) =  gprhs(igaus) + qradf * 0.6_rp* gplad* exp(-0.6_rp*gplad *(gphei-gphot))
           gprhs(igaus) =  gprhs(igaus) +   1.6666_rp*qradf * 0.6_rp* gplad* exp(-0.6_rp*1.6666_rp*gplad *(gphei-gphot))
        end if
         
     end do
     
  end if
  !
  ! Low-Mach: dp0/dt
  ! prefactor alpha * T ~ 1 
  !
  if( kfl_regim_tem >= 3 ) then
     do igaus = 1,pgaus
        gprhs(igaus) = gprhs(igaus) + dpthe 
     end do

  end if
  !
  ! GPGRT: grad(T)
  !
  if( kfl_ellen_tem == -1.or.kfl_prtur_abl_ker==1_ip ) then
     do igaus = 1,pgaus
        do idime = 1,ndime
           gpgrt(idime,igaus) = 0.0_rp
        end do
        do inode = 1,pnode
           do idime = 1,ndime
              gpgrt(idime,igaus) = gpgrt(idime,igaus) &
                   + gpcar(idime,inode,igaus) * eltem(inode,1)
           end do
        end do
     end do
  end if
  
  if( kfl_prtur_abl_ker==1_ip ) then
     eltke(1:pnode) = untur( 1,lnods(1:pnode),1)
     do igaus = 1,pgaus
        gptke(igaus) = 0.0_rp
        do inode = 1,pnode
           gptke(igaus) = gptke(igaus) + gpsha(inode,igaus) * eltke(inode)
        end do
     end do    
  end if

end subroutine tem_elmpre

!-----------------------------------------------------------------------
!> 
!> @author  guillaume
!> @date    2021-06-14
!> @brief   Source term
!> @details Compute source term
!> 
!-----------------------------------------------------------------------

subroutine tem_source(ielem,pmate,gpcod,gpsou)

  use def_kintyp,                  only : ip,rp
  use def_master,                  only : cutim
  use def_domain,                  only : xfiel
  use def_domain,                  only : ndime
  use mod_ker_space_time_function, only : ker_space_time_function
  use def_temper,                  only : kfl_sourc_tem
  use def_temper,                  only : xsour_material_tem
  use def_temper,                  only : lsour_material_tem
  use def_temper,                  only : kfl_sonum_tem
  use def_temper,                  only : sourc_tem
  use def_temper,                  only : SOURCE_TERM_SPACE_TIME
  use def_temper,                  only : SOURCE_TERM_FIELD
  use def_temper,                  only : SOURCE_TERM_MATERIAL
  implicit none
  integer(ip), intent(in)    :: ielem
  integer(ip), intent(in)    :: pmate
  real(rp),    intent(out)   :: gpcod(ndime)
  real(rp),    intent(out)   :: gpsou
  integer(ip)                :: mxdim
  !
  ! Source term: GPRHS
  !
  if( kfl_sourc_tem == SOURCE_TERM_MATERIAL ) then
     !
     ! Material
     !
     if( lsour_material_tem(pmate) == 2 ) then
        gpsou = xsour_material_tem(1,pmate) * sourc_tem
     end if

  else if( kfl_sourc_tem == SOURCE_TERM_SPACE_TIME ) then
     !
     ! Space time
     !
     mxdim = min(2_ip,ndime)
     call ker_space_time_function(&
          kfl_sonum_tem,gpcod(1),gpcod(mxdim),gpcod(ndime),cutim,gpsou)
     gpsou = gpsou * sourc_tem

  else if( kfl_sourc_tem == SOURCE_TERM_FIELD ) then
     !
     ! Field
     !
     gpsou = xfiel(kfl_sonum_tem) % a(1,ielem,1) * sourc_tem

  else
     gpsou = 0.0_rp     
  end if

end subroutine tem_source
