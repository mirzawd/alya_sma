!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_begrun.f90
!> @date    14/06/2019
!> @author  Guillaume Houzeaux
!> @brief   Beginning the run... 
!> @details Beginning the run... we can compute matrices!
!> @}
!-----------------------------------------------------------------------

subroutine nsi_begrun()

  use def_master
  use def_nastin
  use def_domain
  use mod_nsi_fractional_step,  only : nsi_fractional_step_matrices
  use mod_nsi_multi_step_fs,    only : nsi_multi_step_fs_matrices
  use mod_nsi_semi_implicit,    only : nsi_semi_implicit_matrices
  use mod_nsi_schur_complement, only : nsi_schur_complement_matrices
  implicit none
  integer(ip) :: ipoin,iboun,ibopo,ii,found
  
  !-------------------------------------------------------------------
  !
  ! Solution strategies
  !
  !-------------------------------------------------------------------
  
  if( NSI_FRACTIONAL_STEP ) then
     !
     ! Graction astep
     !
     if(kfl_tisch_nsi == 4) then
        call nsi_multi_step_fs_matrices()
     else
        call nsi_fractional_step_matrices()       
     end if
     
  else if(NSI_SEMI_IMPLICIT ) then
     call nsi_semi_implicit_matrices()       
     
  else if( NSI_SCHUR_COMPLEMENT ) then
     !
     ! Schur complement
     !
     call nsi_schur_complement_matrices()       
  end if
  !
  ! If pressure matrix comes from a Schur complement
  ! Must be here because amatr must be allocated
  !
  if( NSI_SCHUR_COMPLEMENT ) then
     solve(2) % A1       => amatr(poapp_nsi:)
     solve(2) % A2       => amatr(poapu_nsi:)
     solve(2) % A3       => amatr(poauu_nsi:)
     solve(2) % A4       => amatr(poaup_nsi:)
     solve(2) % ndofn_A3 =  ndime
     nullify(solve(2) % invA3)
  end if
  
  !-------------------------------------------------------------------
  !
  ! For time & space BC from file: get boundary nodes
  !
  !-------------------------------------------------------------------
  
  if (kfl_bnods_nsi == 1) then
     ibopo = 0
     do ipoin = 1,npoin
        iboun_nsi(ipoin,1) = 0_ip
        do iboun = 1,nbnod_nsi
           if(lninv_loc(ipoin) == int(bntab_nsi(iboun,1),ip)) then
              ii = 1_ip
              found = 0_ip
              do while (found == 0_ip)
                 if (nbnod_pos_nsi(ii) < iboun .and. iboun <= nbnod_pos_nsi(ii+1)) then
                    found = 1_ip
                    iboun_nsi(ipoin,2) = ii
                 else
                    ii = ii + 1_ip
                 end if
              end do
              iboun_nsi(ipoin,1) = iboun - nbnod_pos_nsi(ii)
              ibopo = ibopo + 1_ip
           end if
        end do
     end do
     !if (ibopo /= 0_ip) write(*,*) 'total number of boundary nodes: ', nbnod_nsi, '. computed by processor ', kfl_paral,': ', ibopo
  end if

end subroutine nsi_begrun
