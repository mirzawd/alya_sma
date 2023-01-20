!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_refere()
  !------------------------------------------------------------------------
  !****f* Nastin/nsi_refere
  ! NAME 
  !    nsi_refere
  ! DESCRIPTION
  !    This routine computes the error with respect to a reference solution
  ! USES
  ! USED BY
  !    nsi_output
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastin
  use def_inpout
  use mod_communications, only : PAR_SUM
  implicit none
  integer(ip) :: idime,ipoin,itotn,ntotn
  real(rp)    :: residual(4)

  if( kfl_refer_nsi > 0 ) then 

     if( INOTMASTER ) then
        
        residual = 0.0_rp
        ntotn    = npoin*ndime
        
        do ipoin = 1,npoin
           itotn = (ipoin-1)*ndime + 1
           do idime = 1,ndime
              itotn = itotn + 1
              residual(1) = residual(1) + ( xfiel(kfl_refer_nsi) % a(idime,ipoin,1) - unkno(itotn) ) ** 2
              residual(2) = residual(2) + ( xfiel(kfl_refer_nsi) % a(idime,ipoin,1) ) ** 2
           end do
           residual(3) = residual(3) + ( xfiel(kfl_refer_nsi+1) % a(1,ipoin,1) - unkno(ntotn+ipoin) ) ** 2
           residual(4) = residual(4) + ( xfiel(kfl_refer_nsi+1) % a(1,ipoin,1) ) ** 2
        end do
     end if

     call PAR_SUM(4_ip,residual,'IN THE WORLD')
     
     difve_nsi = sqrt( residual(1) / (residual(2)+zeror) )
     difpr_nsi = sqrt( residual(3) / (residual(4)+zeror) )

  end if

end subroutine nsi_refere
