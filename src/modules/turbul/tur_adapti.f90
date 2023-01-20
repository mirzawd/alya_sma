!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_adapti()
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_adapti
  ! NAME 
  !    tur_updbcs
  ! DESCRIPTION
  !    This routine deals with adaptive boundary conditions
  ! USED BY
  !    tur_begite
  !    tur_iniunk
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_turbul
  implicit none
  integer(ip)   :: ipoin,ibopo,idime,iunkn
  real(rp)      :: udotn

  if( kfl_adapt_tur==1 .and. INOTMASTER ) then
     do ipoin=1,npoin
        ibopo = lpoty(ipoin)
        if(ibopo/=0) then
           if(abs(kfl_fixno_tur(1,ipoin,1))==8) then
              udotn=0.0_rp
              do idime=1,ndime
                 udotn=udotn+advec(idime,ipoin,1)*exnor(idime,1,ibopo)                    
              end do
              do iunkn=1,nturb_tur
                 if(udotn>=0.0_rp) then
                    kfl_fixno_tur(1,ipoin,iunkn)=-8              ! Outflow
                 else
                    kfl_fixno_tur(1,ipoin,iunkn)= 8              ! Inflow
                 end if
              end do
           end if
        end if
     end do
  end if

end subroutine tur_adapti
