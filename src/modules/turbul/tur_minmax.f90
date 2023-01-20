!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_minmax(turmi, turma, nutmi, nutma)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_minmax
  ! NAME 
  !    tur_minmax
  ! DESCRIPTION
  !    This routine computes some min and max
  ! USES
  ! USED BY
  !    nsi_cvgunk
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_turbul
  use mod_communications
  implicit none	
  real(rp), intent(out)   :: turmi(2), turma(2), nutmi, nutma
  integer(ip), parameter  :: nmima	=3
  integer(ip)             :: ipoin, iunk2
  real(rp),    target     :: rmini(nmima),rmaxi(nmima)  
  
  turmi(1) =  1.0e25_rp
  turma(1) = -1.0e25_rp
  turmi(2) = turmi(1)
  turma(2) = turma(1)
  nutmi    =  1.0e25_rp
  nutma    = - 1.0e25_rp
  !
  ! TURMI and TURMA: Min and max turbul unknowns

  !
  if( INOTMASTER ) then
     iunk2 = 1+ nturb_tur-iunkn_tur ! is the other unknow that iunkn_tur
     do ipoin=1, npoin
        turmi(iunkn_tur) = min(unkno(ipoin),turmi(iunkn_tur))
        turma(iunkn_tur) = max(unkno(ipoin),turma(iunkn_tur))
        turmi(iunk2) = min(untur(iunk2,ipoin,1),turmi(iunk2))
        turma(iunk2) = max(untur(iunk2,ipoin,1),turma(iunk2))
        nutmi        = min(turmu(ipoin),nutmi)  
        nutma        = max(turmu(ipoin),nutma)
     end do
     
  end if
  !
  ! Call to parall service
  !
  if( IPARALL ) then
     !Minimum
     rmini(1)  =  turmi(1)
     rmini(2)  =  turmi(2)
     rmini(3)  =  nutmi
     call PAR_MIN(3_ip,rmini)
     turmi(1) =  rmini(1)
     turmi(2) =  rmini(2)
     nutmi    =  rmini(3)
     
     ! maximum
     rmaxi(1)  =  turma(1)
     rmaxi(2)  =  turma(2)
     rmaxi(3)  =  nutma
     call PAR_MAX(3_ip,rmaxi)
     turma(1) =  rmaxi(1)
     turma(2) =  rmaxi(2) 
     nutma    =  rmaxi(3) 

  end if

end subroutine tur_minmax

