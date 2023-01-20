!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_elmsta(&
     pgaus,gpvel,gpdif,gprea,gpden,chale,gpsta)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_elmsta
  ! NAME
  !   tur_elmsta
  ! DESCRIPTION
  !    Compute the stabilization parameter
  ! OUTPUT
  !    GPSTA ... Stabilization parameter with unit [T/rho]=[T*L^3/M]
  ! USES
  ! USED BY
  !    tur_matrix
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime
  use def_turbul, only     :  staco_tur,kfl_advec_tur,kfl_taust_tur,&
       &                      kfl_algor_tur,iunkn_tur,TUR_K_EPS_V2_F,&
       &                      TUR_K_EPS_PHI_F
  use mod_tauadr, only     :  tauadr
  implicit none
  integer(ip), intent(in)  :: pgaus
  real(rp),    intent(in)  :: gpvel(ndime,pgaus)
  real(rp),    intent(in)  :: gpdif(kfl_algor_tur,pgaus)
  real(rp),    intent(in)  :: gprea(kfl_algor_tur,kfl_algor_tur,pgaus)
  real(rp),    intent(in)  :: gpden(pgaus)
  real(rp),    intent(in)  :: chale(2)
  real(rp),    intent(out) :: gpsta(kfl_algor_tur,pgaus)
  integer(ip)              :: igaus
  real(rp)                 :: gpnve,stac2_tur(3)

  if( (TUR_K_EPS_V2_F.or.TUR_K_EPS_PHI_F).and.iunkn_tur==4.and.kfl_algor_tur==1) then
     gpsta=0.0_rp
  else
     if(kfl_algor_tur==1) then
        if(kfl_advec_tur==0) then
           stac2_tur(1)=0.0_rp       ! A
           stac2_tur(2)=staco_tur(2) ! D
           stac2_tur(3)=staco_tur(3) ! R
           do igaus=1,pgaus
              call tauadr(&
                   kfl_taust_tur,staco_tur,gpnve,gpdif(1,igaus),gprea(1,1,igaus),&
                   chale(1),chale(2),gpsta(1,igaus))
           end do
        else
           do igaus=1,pgaus
              call vecnor(gpvel(1,igaus),ndime,gpnve,2_ip) 
              gpnve=gpden(igaus)*gpnve
              call tauadr(&
                   kfl_taust_tur,staco_tur,gpnve,gpdif(1,igaus),gprea(1,1,igaus),&
                   chale(1),chale(2),gpsta(1,igaus))
           end do
        end if
     else
        if(kfl_advec_tur==0) then
           stac2_tur(1)=0.0_rp       ! A
           stac2_tur(2)=staco_tur(2) ! D
           stac2_tur(3)=staco_tur(3) ! R
           do igaus=1,pgaus
              call tauadr(&
                   kfl_taust_tur,staco_tur,gpnve,gpdif(1,igaus),gprea(1,1,igaus),&
                   chale(1),chale(2),gpsta(1,igaus))
              call tauadr(&
                   kfl_taust_tur,staco_tur,gpnve,gpdif(2,igaus),gprea(2,2,igaus),&
                   chale(1),chale(2),gpsta(2,igaus))
           end do
        else
           do igaus=1,pgaus
              call vecnor(gpvel(1,igaus),ndime,gpnve,2_ip)        
              gpnve=gpden(igaus)*gpnve
              call tauadr(&
                   kfl_taust_tur,staco_tur,gpnve,gpdif(1,igaus),gprea(1,1,igaus),&
                   chale(1),chale(2),gpsta(1,igaus))
              call tauadr(&
                   kfl_taust_tur,staco_tur,gpnve,gpdif(2,igaus),gprea(2,2,igaus),&
                   chale(1),chale(2),gpsta(2,igaus))
           end do
        end if
     end if
  end if

end subroutine tur_elmsta
