!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_elmtss(&
     pelty,pmate,pnode,lnods,ielem,elcod,elvel,&
     gpcar,chale,hleng,dtcri)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmtss
  ! NAME
  !    nsi_elmtss
  ! DESCRIPTION
  !    This routine computes the element time step
  ! USED BY
  !    nsi_updtss
  !***
  !-----------------------------------------------------------------------
  use def_kintyp,     only   :  ip,rp
  use def_master,     only   :  ID_NASTIN,ID_CHEMIC
  use def_domain,     only   :  mnode,ndime,elmar
  use def_nastin,     only   :  ncoef_nsi,kfl_taust_nsi,&
       &                        staco_nsi,&
       &                        kfl_advec_nsi,&
       &                        corio_nsi
  use def_nastin,     only   :  kfl_surte_nsi
  use def_nastin,     only   :  surte_nsi
  use def_kermod,     only   :  turmu_ker
  use def_kermod,     only   :  densi_ker
  use mod_ker_proper, only   :  ker_proper
  use def_kermod,     only   :  kfl_kxmod_ker
  use mod_tauadr,     only   :  tauadr
  use def_parame,     only   :  pi

  implicit none
  integer(ip), intent(in)    :: pelty,pmate,pnode
  integer(ip), intent(in)    :: lnods(pnode),ielem
  real(rp),    intent(in)    :: elcod(ndime,mnode)
  real(rp),    intent(in)    :: elvel(ndime,mnode)
  real(rp),    intent(in)    :: chale(ndime),hleng(ndime)
  real(rp),    intent(inout) :: dtcri
  real(rp),    intent(out)   :: gpcar(ndime,mnode)
  integer(ip)                :: inode, taust
  real(rp)                   :: xjaci(9),xjacm(9)
  real(rp)                   :: gpdet,gpvis(1),gpden(1),gppor(1),gpgvi(3)
  real(rp)                   :: gpvel(3),gplev,gpvno,adv,dif,rea
  real(rp)                   :: gpgve(ndime,ndime),gppre,gptem,rnode
  real(rp)                   :: gpnut(1),grvis(ndime,ndime)
  real(rp)                   :: rho_L,rho_g,vol,dtCap,uCap
  !
  ! Velocity
  !
  gplev = 0.0_rp
  gppre = 0.0_rp
  gpvel = 0.0_rp
  gpgve = 0.0_rp
  gptem = 0.0_rp
  gpnut = 0.0_rp
  gpvno = 0.0_rp
  dtcri = 0.0_rp
  rnode = 1.0_rp/real(pnode,rp)
  !
  ! GPCAR: Cartesian derivative
  !
  call elmder(&
       pnode,ndime,elmar(pelty)%dercg,elcod,&
       gpcar,gpdet,xjacm,xjaci)
  !
  ! Values at center of gravity
  !
  if( kfl_advec_nsi /= 0 ) then                        ! GPVEL: Velocity
     do inode = 1,pnode
        gpvel(1:ndime) = gpvel(1:ndime) + elvel(1:ndime,inode)
     end do
     gpvel = rnode*gpvel
     gpvno = sqrt(dot_product(gpvel(1:ndime),gpvel(1:ndime)))
  end if
  !
  ! GPVIS, GPDEN, GPPOR: Compute properties
  !
  call ker_proper('DENSI','COG  ',1_ip,ielem,gpden)
  call ker_proper('VISCO','COG  ',1_ip,ielem,gpvis)
  call ker_proper('POROS','COG  ',1_ip,ielem,gppor)
  call ker_proper('TURBU','COG  ',1_ip,ielem,gpnut)

 
  if(turmu_ker % kfl_exist /= 0_ip) then ! matias suggested I use this - seems good idea 

     gpgvi = 0.0_rp
     grvis = 0.0_rp
     !adds gpmut to turbul
     call nsi_turbul(&
          1_ip,0_ip,pnode,1_ip,1_ip,1_ip,&
          elmar(pelty)%shacg,gpcar,elvel,gpden(1),gpvis(1),gpnut,&
          gpgvi,grvis,gpgve,ielem,kfl_kxmod_ker)
  end if
  !
  ! DTCRI: Critical time step
  !
  !! CORIOLIS IS ZERO IN THE FIRST TIME STEP BECAUSE IT IS ACTUALIZED IN NSI_UPDFOR, CALLED BY NSI_BEGITE
  adv        = gpden(1) * gpvno                      ! Convective term
  dif        = gpvis(1)                              ! Viscous term (accounts for gptur) 
  rea        = gpden(1) * corio_nsi + abs(gppor(1))  ! Porosity + Coriolis term
  !! until this branch the diffusive term was not accounting for turbulent viscosity in the critical time step. I believe is an important modification. Now, gpvis = mu + mu_t

  ! When tau depends on time step. Let us use codina
  if (kfl_taust_nsi.gt.4) then
     taust=1
  else
     taust = kfl_taust_nsi
  endif
     
  call tauadr(&
     taust,staco_nsi,adv,dif,rea,&
     chale(1),chale(2),dtcri) 

  if( adv == 0.0_rp .and. dif == 0.0_rp .and. rea == 0.0_rp ) dtcri = huge(1.0_rp)

  dtcri = gpden(1)*dtcri

  !
  ! Capilarity effects in time-step
  ! Brackville et al., JCP (1992)
  ! 
  ! dt_cap <= min{0.5*(rho_1 + rho_2) / (4 * pi() * sigma) * Delta^(3/2)} 
  ! 
  if (kfl_surte_nsi == 2) then 

     rho_L = densi_ker % rlaws(1,1)
     rho_g = densi_ker % rlaws(2,1)
     if (ndime == 3) then
        vol = ( hleng(1) * hleng(2) * hleng(3) )
     else
        vol = ( hleng(1) * hleng(2) )
     end if

     uCap = ( rho_L + rho_g ) / (4.0_rp * pi * surte_nsi) 

     dtCap = ( uCap * vol )**(0.5_rp)

     dtcri = min(dtcri,dtCap)

  end if


end subroutine nsi_elmtss
