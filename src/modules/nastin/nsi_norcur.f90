!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_norcur
  !------------------------------------------------------------------------
  !****f* Nastin/nsi_norcur
  ! NAME 
  !    nsi_norcur
  ! DESCRIPTION
  !    This routine computes normal and curvature for surface tension
  ! USES
  ! USED BY
  !    nsi_norcur
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastin
  use mod_memchk
  use mod_postpr
  use mod_gradie
  use mod_memory, only: memory_alloca, memory_deallo
  implicit none
  integer(ip)             :: ipoin,idime,iclas_phi
  real(rp)                :: normn
  real(rp), pointer       :: mask_phi(:)

  select case (kfl_surte_nsi)

  !
  ! Classical FEM level set momentum conservation
  !
  case(1_ip)

     if ( kfl_modul(ID_LEVELS) == 0 ) call runend(' NSI_NORCUR: only makes sense with LEVEL SET coupling')

     call gradie(fleve(1:npoin,1),norle_nsi)    ! grad LS
     
     do ipoin = 1,npoin  ! Normal = grad LS / | grad LS |
        normn = sqrt (dot_product ( norle_nsi(1:ndime,ipoin) , norle_nsi(1:ndime,ipoin) ) )
        if (normn /= 0.0_rp) then
           do idime = 1,ndime
              norle_nsi(idime,ipoin) = norle_nsi(idime,ipoin) / normn
           end do
        end if
     end do

     call divvec(norle_nsi,curle_nsi)    ! curvature = div NORMAL

  !
  ! ELSA model momentum conservation
  !
  case(2_ip)

     if (INOTMASTER) then
        !
        ! Identify iclas for volume fraction phi_L
        !
        iclas_phi = 3

        nullify(mask_phi) 
        call memory_alloca(mem_modul(1:2,modul),'mask_phi','nsi_norcur',mask_phi,npoin)

        do ipoin = 1,npoin
           mask_phi(ipoin) = 1.0_rp 
           if ( (conce(ipoin,iclas_phi,1)<0.005_rp) .or. (conce(ipoin,iclas_phi,1)>0.995_rp) ) then
              !
              ! TODO This has to be revised
              !
              mask_phi(ipoin) = 1.0_rp
           endif
        enddo
        
        !
        ! Compute normal vector n
        !
        call gradie(conce(1:npoin,iclas_phi,1),norle_nsi)

        do ipoin = 1,npoin  ! Normal = grad LS / | grad LS |
           normn = sqrt (dot_product ( norle_nsi(1:ndime,ipoin) , norle_nsi(1:ndime,ipoin) ) )
           if (normn /= 0.0_rp) then
              do idime = 1,ndime
                 norle_nsi(idime,ipoin) = norle_nsi(idime,ipoin) / normn
              end do
           end if
        end do

        !
        ! Compute curvature K(phi_L) = div NORMAL
        !
        call divvec(norle_nsi,curle_nsi)

        !
        ! Re-compute gradient volume fraction phi_L (unscaled)
        !
        call gradie(conce(1:npoin,iclas_phi,1),norle_nsi)


        do ipoin= 1,npoin
           curle_nsi(ipoin) = curle_nsi(ipoin) * mask_phi(ipoin)
           norle_nsi(:,ipoin) = norle_nsi(:,ipoin) * mask_phi(ipoin)
        enddo

        
        call memory_deallo(mem_modul(1:2,modul),'mask_phi','nsi_norcur',mask_phi)

     end if

  end select

end subroutine nsi_norcur
