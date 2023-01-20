!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_elmtss(&
     pelty,pnode,ielem,elvel,eledd,dtcri,eltur,elfle,&
     lnods,chale,hleng)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_elmtss
  ! NAME 
  !    tur_elmtss
  ! DESCRIPTION
  !    This routine computes the element time step 
  !                     1
  !    dt = -------------------------------
  !          4*(mu*smole+mut*sturb)     2*u
  !          ----------------------  +  --- 
  !                 rho*h^2              h 
  ! USED BY
  !    tur_updtss
  !    tur_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime
  use def_turbul, only       :  nturb_tur,iunkn_tur,&
       &                        param_tur,TUR_SPALART_ALLMARAS,&
       &                        kfl_advec_tur,staco_tur
  use def_kermod
  use mod_ker_proper, only   :  ker_proper
  use mod_tauadr, only       :  tauadr

  implicit none
  integer(ip), intent(in)    :: pelty,pnode
  integer(ip), intent(in)    :: lnods(pnode),ielem
  real(rp),    intent(in)    :: elvel(ndime,pnode)
  real(rp),    intent(in)    :: eledd(pnode)
  real(rp),    intent(in)    :: eltur(nturb_tur,pnode)  
  real(rp),    intent(in)    :: elfle(pnode)  
  real(rp),    intent(in)    :: chale(2),hleng(*)
  real(rp),    intent(inout) :: dtcri
  integer(ip)                :: idime,inode,dummi
  real(rp)                   :: gpvis(1),gpden(1),gpvno,gpnut, gptvi(1)
  real(rp)                   :: smole,sturb,rnode,gpvel(3),adv,dif,rea
  !
  ! Initialization
  !
  gpvno = 0.0_rp
  gpvel = 0.0_rp
  rnode = 1.0_rp/real(pnode,rp)
  !
  ! GPDEN and GPVIS: Properties
  !
  call ker_proper('DENSI','COG  ',dummi,ielem,gpden)
  call ker_proper('VISCO','COG  ',dummi,ielem,gpvis)
  call ker_proper('TURBU','COG  ',dummi,ielem,gptvi)
  gpnut = gptvi(1)

  !
  ! GPNUT: Model dependent variables
  !
 
  if(TUR_SPALART_ALLMARAS) then 
     smole = 1.0_rp/param_tur(3)
     sturb = smole
  else 
     smole = 1.0_rp
     sturb = 1.0_rp/param_tur(iunkn_tur)
  end if

  !
  ! GPVNO: velocity norm
  !
  if(kfl_advec_tur/=0) then
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
  ! DTCRI: Critical time step
  ! 
  adv = gpvno
  dif = gpvis(1)/gpden(1)*smole+gpnut*sturb
  rea = 0.0_rp    
  call tauadr(&
       1_ip,staco_tur,adv,dif,rea,&
       chale(1),chale(2),dtcri)

end subroutine tur_elmtss
