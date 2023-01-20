!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_updrel()
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_updrel
  ! NAME 
  !    tur_updrel
  ! DESCRIPTION
  !    This routine updates the relaxation factor
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_turbul
  implicit none
  integer(ip) :: ipoin,idofn

  if(kfl_relax_tur==2) then
     !
     ! Aitken's method
     !
     if(ittot_tur/=1) then 
        call aitken(&
             0_ip,nturb_tur,1_ip,1_ip,untur,unkno,dunkn_tur,&
             one,one,one,one,relpa_tur)
     else
        if(kfl_paral/=0) then
           idofn=(iunkn_tur-1)*npoin
           do ipoin=1,npoin
              idofn=idofn+1
              dunkn_tur(idofn)=untur(iunkn_tur,ipoin,1)-unkno(ipoin)
           end do
        end if
        relpa_tur(1) = 0.0_rp
     end if
     relax_tur = 1.0_rp-relpa_tur(1)
  end if

end subroutine tur_updrel
