!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_coupli.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Nastal coupling 
!> @details Nastal coupling
!> @} 
!-----------------------------------------------------------------------
subroutine nsi_coupli(itask)
  use def_master
  use def_domain
  use def_elmtyp
  use def_nastin
  use def_coupli,      only : FIXED_UNKNOWN
  use mod_local_basis, only : local_basis_global_to_local
  use mod_local_basis, only : local_basis_local_to_global
  implicit none
  integer(ip), intent(in) :: itask

  select case ( itask )

  case ( ITASK_CONCOU )


  case ( ITASK_BEGSTE )


  case ( ITASK_INIUNK )

     !-------------------------------------------------------------------
     !
     ! INIUNK
     !
     !-------------------------------------------------------------------

  case ( ITASK_BEGITE )

     !-------------------------------------------------------------------
     !
     ! BEGITE
     !
     !-------------------------------------------------------------------
     
  case ( ITASK_ENDITE )

     !-------------------------------------------------------------------
     !
     ! ENDITE
     !
     !-------------------------------------------------------------------
     !
     ! Coupling with Immbou: compute force
     !
     if( coupling('IMMBOU','NASTIN') >= 1 ) then
        call nsi_immbou()
     end if
     !
     ! Coupling with Alefor: compute force for rigid body
     !
     if( ( coupling('ALEFOR','NASTIN') >= 1 ) .and. ( nrbod > 0_ip ) ) then
        call nsi_rbobou()
     end if

  case ( ITASK_MATRIX )

     !-------------------------------------------------------------------
     !
     ! MATRIX: After assembly
     !
     !-------------------------------------------------------------------

     if( coupling('NASTIN','IMMBOU') == 1 ) then
        !
        ! Coupling with Immbou: mass matrix conservation variables
        !                
        call nsi_massma()
        if( INOTMASTER ) then
           !
           ! Coupling with Immbou: impose force FORCF: interpolate b.c.
           !     
           call nsi_embedd(&
                amatr(poauu_nsi),amatr(poaup_nsi),amatr(poapu_nsi),amatr(poapp_nsi),&
                lapla_nsi,rhsid,rhsid(ndbgs_nsi+1),unkno,unkno(ndbgs_nsi+1))           
        end if
     end if
     !
     ! Coupling with Solidz in IBM: add force from solidz 
     !
     call nsi_dimbou()
     !
     ! Coupling with Partis: take off momentum from momentum equations
     !    
     call nsi_partis()

  end select

end subroutine nsi_coupli
