!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_inibcs.f90
!> @author  Guillaume Hozeaux
!> @date    January, 2011
!>          - Subroutine creation
!> @author  Gerard Guillamet
!> @date    January, 2019
!>          - Local axes and non-constant bcs
!> @brief   This routine reads the boundary conditions
!> @details
!>          For conditions on nodes, bvess_sld(idime,ipoin,1) contains the\n
!>          displacement value.\n
!>          The different codes for kfl_fixno_sld(ipoin) are:\n
!>          =  1 ... Dirichlet\n
!>          =  0 ... Free or initial\n
!>          =  2 ... Pericardium\n
!>          =  3 ... Contact nodes\n
!>
!>          For conditions on boundaries, bvnat_sld(idime,iboun,1) is\n
!>          allocated ONLY if a condition is imposed on it. Moreover,\n
!>          its length depends on the condition type.\n
!>          The different codes for kfl_fixbo_sld(iboun) are:\n
!>          =  1 ... Dirichlet ................... u\n
!>          =  2 ... \n
!>          =  3 ... Normal Pressure\n
!>          =  4 ... \n
!>          =  5 ... \n
!>          =  6 ... \n
!>          =  7 ... \n
!>          =  8 ... \n
!>
!>          Local boundary conditions, kfl_local_sld\n
!>          =  0 ... Global axes\n
!>          =  1 ... Local axes\n
!>          The different codes for kfl_fixrs_sld(ipoin) are:\n
!>          =  0 ... Exterior normal\n
!>          = -1 ... User Local axes\n
!>
!>          At the end of the subroutine conditions on boundaries are\n
!>          transfered to conditions on nodes.\n
!> @}
!-----------------------------------------------------------------------

subroutine sld_inibcs()

  use def_kintyp,  only : ip, rp
  use def_master,  only : INOTEMPTY, INOTMASTER
  use def_master,  only : ITER_K, ITER_AUX
  use def_master,  only : IMPOSE_NODE_CODES, IMPOSE_BOUNDARY_CODES
  use def_domain,  only : lpoty
  use def_domain,  only : npoin, ndime, nboun
  use def_domain,  only : bvess, bvnat
  use def_domain,  only : ifbop, iffun, ifloc
  use def_domain,  only : kfl_fixno, kfl_fixbo, kfl_fixrs
  use def_domain,  only : kfl_funno, kfl_funbo, kfl_funtb, kfl_funtn
  use def_domain,  only : kfl_icodn, kfl_icodb
  use def_domain,  only : tncod, tbcod
  use mod_communications, only : PAR_MAX
  use def_solidz,  only : bvess_sld, bvnat_sld
  use def_solidz,  only : pressEndoIntegralNew_sld,pushForwardIntegralNew_sld
  use def_solidz,  only : pressEndoIntegralOld_sld,pushForwardIntegralOld_sld
  use def_solidz,  only : kfl_fixno_sld, kfl_fixbo_sld, kfl_fixrs_sld
  use def_solidz,  only : kfl_funno_sld, kfl_funbo_sld, kfl_funtn_sld
  use def_solidz,  only : kfl_funtb_sld
  use def_solidz,  only : kfl_conbc_sld, kfl_local_sld
  use def_solidz,  only : tncod_sld, tbcod_sld
  use def_solidz,  only : kfl_rigid_sld, kfl_rbfix_sld

  implicit none

  integer(ip)          :: ipoin,iboun,ibopo,idime

  !-------------------------------------------------------------
  !
  ! Some boundary integrals initializations
  !
  !-------------------------------------------------------------
  
  pressEndoIntegralNew_sld   = 0.0_rp
  pushForwardIntegralNew_sld = 0.0_rp
  pressEndoIntegralOld_sld   = 0.0_rp
  pushForwardIntegralOld_sld = 0.0_rp
  
  !-------------------------------------------------------------
  !
  ! Allocate memory for the vectors needed to define the BC's
  !
  !-------------------------------------------------------------

  call sld_membcs(1_ip)

  if ( INOTEMPTY ) then

     !-------------------------------------------------------------
     !
     ! Node codes
     !
     !-------------------------------------------------------------

     if ( kfl_icodn > 0 ) then
        if ( kfl_conbc_sld == 0 ) then
           iffun      =  1
           kfl_funno  => kfl_funno_sld
           kfl_funtn  => kfl_funtn_sld
        else
           iffun      =  0
        end if
        ifbop     =  0
        ifloc     =  1
        kfl_fixrs => kfl_fixrs_sld
        kfl_fixno => kfl_fixno_sld
        bvess     => bvess_sld(:,:,ITER_K)
        tncod     => tncod_sld(:)
        call reacod(IMPOSE_NODE_CODES)

     end if

     !-------------------------------------------------------------
     !
     ! Boundary codes
     !
     !-------------------------------------------------------------

     if ( kfl_icodb > 0 ) then
        if ( kfl_conbc_sld == 0 ) then
           iffun      =  1
           kfl_funbo   => kfl_funbo_sld
           kfl_funtb   => kfl_funtb_sld
        end if
        kfl_fixbo => kfl_fixbo_sld
        bvnat     => bvnat_sld(:,:,ITER_K)
        tbcod     => tbcod_sld(1:)
        !tncod     => momod(ID_NASTIN) % tncod(1:)
        call reacod(IMPOSE_BOUNDARY_CODES)
     end if

     !-------------------------------------------------------------
     !
     ! Exact solution
     !
     !-------------------------------------------------------------

     call sld_exaerr(1_ip)

     !-------------------------------------------------------------
     !
     ! Check if there is a local prescription
     !
     !-------------------------------------------------------------

     nodes: do ipoin = 1,npoin
        ibopo = lpoty(ipoin)
        if ( ibopo > 0 ) then
           if ( kfl_fixrs_sld(ipoin) /= 0 ) then
              kfl_local_sld = 1
              exit nodes
           end if
        end if
     end do nodes

     !-------------------------------------------------------------
     !
     ! Non-constant boundary conditions
     !
     !-------------------------------------------------------------

     if ( kfl_conbc_sld == 0 ) then
        do ipoin = 1,npoin
           bvess_sld(1:ndime,ipoin,ITER_AUX) = bvess_sld(1:ndime,ipoin,ITER_K)
        end do
        do iboun = 1,nboun
           bvnat_sld(1:ndime,iboun,ITER_AUX) = bvnat_sld(1:ndime,iboun,ITER_K)
        end do
     end if

  end if

  !-------------------------------------------------------------
  !
  ! Rigid body fixity
  !
  !-------------------------------------------------------------

  if( kfl_rigid_sld == 1 ) then

     kfl_rbfix_sld(:) = 0
     if( INOTMASTER ) then
        do idime = 1,ndime
           kfl_rbfix_sld(idime) = maxval(kfl_fixno_sld(idime,:))
        end do
     end if
     !
     ! Should be broadcast because subdomain can modify it while others not
     !
     call PAR_MAX(3_ip,kfl_rbfix_sld,'IN MY CODE')

  end if

end subroutine sld_inibcs
