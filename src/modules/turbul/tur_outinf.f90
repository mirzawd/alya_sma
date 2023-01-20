!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_outinf(itask)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_outinf
  ! NAME 
  !    tur_outinf
  ! DESCRIPTION
  !    This routine writes some info in output files 
  ! USES
  ! USED BY
  !    tur_outinf
  !***
  !-----------------------------------------------------------------------
  use def_master
  use mod_outfor, only : outfor
  use def_turbul
  implicit none
  integer(ip), intent(in) :: itask
  character(100)          :: equat

  if(kfl_paral<=0) then

     select case(itask)

     case(1_ip)
        !
        ! Write information in Result file
        !
        if(TUR_SPALART_ALLMARAS) then

           call outfor(25_ip,momod(modul)%lun_outpu,'DIFFERENTIAL EQUATION')
           equat='dnu''/dt+u.grad(nu'')-1/sig*div.[(nu+nu'')*grad(nu'')]'
           write(momod(modul)%lun_outpu,110) equat
           equat='   =cb1*(1-ft2)*S''*nu''+cb2/sig*grad(nu'')^2-[cw1*fw-cb1*ft2/k^2]*nu''^2/d^2'     
           write(momod(modul)%lun_outpu,111) equat

           write(momod(modul)%lun_outpu,111) ''
           write(momod(modul)%lun_outpu,111) 'Let X=nu''/nu.'
           write(momod(modul)%lun_outpu,111) ''
           write(momod(modul)%lun_outpu,111) 'The destruction function fw is given by'
           write(momod(modul)%lun_outpu,111) 'r   = nu''/(S''*k^2*d^2)'
           write(momod(modul)%lun_outpu,111) 'g   = r+cw2*(r^6-r)'
           write(momod(modul)%lun_outpu,111) 'fw  = g*[(1+cw3^6)/(g^6+cw3^6)]^(1/6)'
           write(momod(modul)%lun_outpu,111) ''
           write(momod(modul)%lun_outpu,111) 'S'' is a function of the vorticity S such that'
           write(momod(modul)%lun_outpu,111) 'S'' = S*fv3+nu''/(k^2*d^2)*fv2'
           write(momod(modul)%lun_outpu,111) ''
           if(ipara_tur(1)==0) then
              write(momod(modul)%lun_outpu,111) 'and the functions are given by SA original model:'
              write(momod(modul)%lun_outpu,111) 'fv1 = X^3/(X^3+cv1^3)'
              write(momod(modul)%lun_outpu,111) 'fv2 = 1-X/(1+fv1*X)'
              write(momod(modul)%lun_outpu,111) 'fv3 = 1'
              write(momod(modul)%lun_outpu,111) 'ft2 = ct3*exp(-ct4*X^2)'
           else if(ipara_tur(1)==1) then
              write(momod(modul)%lun_outpu,111) 'and the functions are given by SA original model'
              write(momod(modul)%lun_outpu,111) 'together with Geuzaine, Delanaye and Liu corrections'
              write(momod(modul)%lun_outpu,111) 'fv1 = X^3/(X^3+cv1^3)'
              write(momod(modul)%lun_outpu,111) 'fv2 = (1+X/cv2)^-3'
              write(momod(modul)%lun_outpu,111) 'fv3 = (1+fv1*X)*(1-fv2)/X'
              write(momod(modul)%lun_outpu,111) 'ft2 = ct3*exp(-ct4*X^2)'
           else
              write(momod(modul)%lun_outpu,111) 'and the functions are given by SA original high Reynolds number model'
              write(momod(modul)%lun_outpu,111) 'fv1 = 1'
              write(momod(modul)%lun_outpu,111) 'fv2 = 0'
              write(momod(modul)%lun_outpu,111) 'fv3 = 1'
              write(momod(modul)%lun_outpu,111) 'ft2 = 0'
           end if

           write(momod(modul)%lun_outpu,113) 
           write(momod(modul)%lun_outpu,111) 'cb1 = ',param_tur(1)
           write(momod(modul)%lun_outpu,111) 'cb2 = ',param_tur(2)
           write(momod(modul)%lun_outpu,111) 'sig = ',param_tur(3)
           write(momod(modul)%lun_outpu,111) 'cv1 = ',param_tur(4)
           write(momod(modul)%lun_outpu,111) 'cw2 = ',param_tur(5)
           write(momod(modul)%lun_outpu,111) 'cw3 = ',param_tur(6)
           write(momod(modul)%lun_outpu,111) 'k   = ',param_tur(7)
           write(momod(modul)%lun_outpu,112) 'cw1 = ',param_tur(8),' = cb1/k^2+(1+cb2)/sig'
           if(ipara_tur(1)==0) then
              write(momod(modul)%lun_outpu,111) 'ct3 = ',1.1_rp
              write(momod(modul)%lun_outpu,111) 'ct4 = ',2.0_rp
           else if(ipara_tur(1)==1) then
              write(momod(modul)%lun_outpu,111) 'cv2 = ',5.0_rp
              write(momod(modul)%lun_outpu,111) 'ct3 = ',1.1_rp
              write(momod(modul)%lun_outpu,111) 'ct4 = ',2.0_rp
           end if
           if (kfl_rotat_tur==1) write(momod(modul)%lun_outpu,111) 'Rotating reference frame (simplified SST-RC-Hellsten model see http://turbmodels.larc.nasa.gov/sst.html)'
           write(momod(modul)%lun_outpu,*) 
        else if(TUR_K_XU_CHIEN) then
           write(momod(modul)%lun_outpu,111) 'Chien k-l model'
        else if(TUR_K_EPS_STD) then
           write(momod(modul)%lun_outpu,111) 'Standard high-Re k-eps model'
           if (kfl_kxmod_tur==1) write(momod(modul)%lun_outpu,111) 'RNG model'
           if (kfl_kxmod_tur==2) write(momod(modul)%lun_outpu,111) 'Realizable model'
           if (kfl_kxmod_tur==3) write(momod(modul)%lun_outpu,111) 'ke-f-p model'
        else if(TUR_K_EPS_LAUNDER_SHARMA) then
           write(momod(modul)%lun_outpu,111) 'Launder-Sharma k-eps model'
        else if(TUR_K_EPS_CHIEN) then
           write(momod(modul)%lun_outpu,111) 'Chien k-eps model'
        else if(TUR_K_EPS_NAGANO_TAGAWA) then
           write(momod(modul)%lun_outpu,111) 'Nagano-Tagawa k-eps model'
        else if(TUR_K_EPS_LAM_BREMHORST) then
           write(momod(modul)%lun_outpu,111) 'Lam-Bremhorst k-eps model'
        else if(TUR_K_EPS_JAW_HWANG) then
           write(momod(modul)%lun_outpu,111) 'Jaw-Hwang k-eps model'
        else if(TUR_K_EPS_V2_F) then
           write(momod(modul)%lun_outpu,111) 'k-eps-v2-f model'
        else if(TUR_K_EPS_PHI_F) then
           write(momod(modul)%lun_outpu,111) 'k-eps-phi-f model'
        else if(TUR_TWO_LAYER_RODI) then
           write(momod(modul)%lun_outpu,111) 'Two-layer k model'
        else if(TUR_TWO_LAYER_XU_CHEN) then
           write(momod(modul)%lun_outpu,111) 'Two-layer Xu-Chen k-eps/k-l model'
        else if(TUR_K_OMEGA) then
           write(momod(modul)%lun_outpu,111) 'Wilcox k-w model'
        else if(TUR_K_OMEGA_BREDBERG) then
           write(momod(modul)%lun_outpu,111) 'Bredberg-Peng-Davidson k-w model'
        end if

     end select

  end if
  !
  ! Formats
  !
110 format(/,&
       & 10x,a)
113 format(/,&
       & 10x,'Constants',/,&
       & 10x,'-----------------------------------------')
111 format(&
       & 10x,a,e12.6)
112 format(&
       & 10x,a,e12.6,a)

end subroutine tur_outinf

