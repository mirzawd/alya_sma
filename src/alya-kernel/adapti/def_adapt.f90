!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Adaptivity
!> @{
!> @file    def_adapt.f90
!> @author  abel.gargallo
!> @date    2021-04-20
!> @brief   def_adapt
!> @details def_adapt
!>
!>          To add further details
!>
!>
!-----------------------------------------------------------------------

MODULE def_adapt
!***************************************************************
!*
!*  Module for performing adaptive topological operation
!*
!***************************************************************
use def_kintyp_basic, only : ip,rp,lg

implicit none
private
public :: memor_adapt

integer(8)               :: memor_adapt(2)

END MODULE def_adapt


!> @}