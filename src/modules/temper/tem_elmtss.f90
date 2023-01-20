!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_elmtss(&
     ielem,pelty,pnode,pmate,elcod,elvel,&
     eltem,gpcar,chale,dtcri)
  !-----------------------------------------------------------------------
  !****f* Temper/tem_elmtss
  ! NAME 
  !    tur_elmtss
  ! DESCRIPTION
  !    This routine computes the element time step
  ! USED BY
  !    tem_updtss
  !    tem_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  mnode,ndime,elmar
  use def_master, only       :  kfl_coupl,div_enthalpy_transport,&
       &                        ID_TEMPER,ID_CHEMIC
  use def_temper, only       :  kfl_advec_tem,prtur_tem,&
       &                        staco_tem,kfl_regim_tem,&
       &                        source_safet_tem
  use def_kermod
  use mod_ker_proper
  use mod_tauadr, only       :  tauadr

  implicit none
  integer(ip), intent(in)    :: ielem,pelty,pnode,pmate
  real(rp),    intent(in)    :: elcod(ndime,pnode)
  real(rp),    intent(in)    :: elvel(ndime,pnode)
  real(rp),    intent(in)    :: eltem(pnode),chale(2)
  real(rp),    intent(inout) :: dtcri
  real(rp),    intent(out)   :: gpcar(ndime,mnode)
  integer(ip)                :: idime,inode
  real(rp)                   :: gpcon(1),gpden(1),gpsph(1),gprcp,gprea,gpdif(1)
  real(rp)                   :: gptem,rnode,gpvno,adv,dif,rea
  real(rp)                   :: gptur(1)
  real(rp)                   :: gpvol,gphes,gpvel(3),chtim,rocpt
  !
  ! Initialization
  !
  gptem = 0.0_rp
  gpvno = 0.0_rp
  gpvel = 0.0_rp
  gptur = 0.0_rp 
  dtcri = 0.0_rp
  rnode = 1.0_rp/real(pnode,rp)
  !
  ! GPTEM: Temperature used when k=k(T), Cp=Cp(T)
  !
  do inode=1,pnode
     gptem=gptem+eltem(inode)
  end do
  gptem=rnode*gptem
  !
  ! GPVNO: Velocity norm
  !
  if(kfl_advec_tem/=0) then 
     do inode=1,pnode
        do idime=1,ndime
           gpvel(idime)=gpvel(idime)+elvel(idime,inode)
        end do
     end do
     gpvel=rnode*gpvel
     do idime=1,ndime
        gpvno=gpvno+gpvel(idime)*gpvel(idime)
     end do
     gpvno=sqrt(gpvno)
  end if
  !
  ! GPDEN, GPSPH, GPDIF, GPREA: Properties
  !
  if(kfl_regim_tem==1.or.kfl_regim_tem==2.or.kfl_regim_tem>=3) then
     call elmcar(&
          pnode,1_ip,0_ip,elmar(pelty)%weicg,elmar(pelty)%shacg,&
          elmar(pelty)%dercg,elmar(pelty)%hescg,elcod,gpvol,gpcar,&
          gphes,ielem)
  end if

  call ker_proper('DENSI','IGAUS',1_ip,ielem,gpden,pnode,1_ip,elmar(pelty)%shacg,gpcar)
  call ker_proper('CONDU','IGAUS',1_ip,ielem,gpcon,pnode,1_ip,elmar(pelty)%shacg,gpcar)
  call ker_proper('SPHEA','IGAUS',1_ip,ielem,gpsph,pnode,1_ip,elmar(pelty)%shacg,gpcar)
  !
  ! Turbulent viscosity mut 
  !
  call ker_proper('TURBU','IGAUS',1_ip,ielem,gptur,pnode,1_ip,elmar(pelty)%shacg,gpcar)
  gptur(1) = gpden(1)*gptur(1)
  ! 
  ! Coupling with turbul
  !
  if(  turmu_ker % kfl_exist /= 0_ip ) then
     gpdif(1) = gpcon(1) + gpsph(1)*gptur(1)/prtur_tem     
  else
     gpdif(1) = gpcon(1)
  end if
  !
  ! Reaction term
  !
  call tem_elmrea( &
       1_ip,pnode,1_ip,1_ip,1_ip,elvel,gpden,gpcar,&
       gprea)
 
  gprcp=gpden(1)*gpsph(1)
  if (kfl_regim_tem == 4) then
    gpdif(1) = gpdif(1)/gpsph(1)
    gprcp = gpden(1) 
  endif
  !
  ! DTCRI: critical time step
  !
  adv = gprcp*gpvno
  dif = gpdif(1)
  rea = gprea

  call tauadr(&
       1_ip,staco_tem,adv,dif,rea,&
       chale(1),chale(2),dtcri)

  dtcri = gprcp*dtcri

  if(kfl_coupl(ID_TEMPER,ID_CHEMIC)==1) then
     ! Compute chemical heat release critical time
     ! The inverse critical time step of the chemical reaction is
     !  ( Chemical Heat Release Rate ) / ( rho * Cp * T)
     !
     ! We add a safety factor to make sure source is not large
     !
     rocpT = gpden(1) * gptem * gpsph(1)
     if (kfl_regim_tem == 4_ip) rocpT = gpden(1) * gptem
     if (div_enthalpy_transport(ielem)%a(1,1,1) /= 0.0_rp) then
        chtim = source_safet_tem*rocpT/(abs(div_enthalpy_transport(ielem)%a(1,1,1)))
     else
        chtim=1.0e12_rp
     endif
     dtcri = min(dtcri,chtim)
  endif

end subroutine tem_elmtss
