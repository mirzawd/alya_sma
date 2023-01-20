!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_sld_bclocal.f90
!> @author  Gerard Guillamet (gerard.guillamet@bsc.es)
!> @date    April, 2019
!>          - Subroutine written
!> @brief   Prescribe BCs using local axes
!>
!> @details Prescribe BCs using local axes, this includes:
!>
!>            - General Dirichlet type
!>            - PDN-contact Dirichlet type
!>            - Rotation RHS and ELMAT
!>
!> @}
!------------------------------------------------------------------------

module mod_sld_bclocal

  use def_kintyp,   only : ip, rp
  use def_domain,   only : lpoty
  use def_solidz,   only : kfl_fixrs_sld, kfl_fixno_sld
  use mod_sld_csys, only : sld_csys_rotuni

  public :: sld_bclocal_displ
  public :: sld_bclocal_rhsid_fixity
  public :: sld_bclocal_roback
  public :: sld_bclocal_elmat_and_elrhs

contains

  !-----------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    April, 2019
  !> @brief   Displacement and unkno rotations due to local axes
  !> @details Displacement and unkno rotations due to local axes
  !>          For the PDN-contact the displacement in the normal
  !>          direction vector is fixed.
  !>          Rotations: LOCAL -> GLOBAL and GLOBAL -> LOCAL
  !-----------------------------------------------------------------------

  subroutine sld_bclocal_displ(ndofn,ipoin)

    use def_master,          only : ITER_K
    use def_master,          only : displ, unkno
    use def_solidz,          only : kfl_timet_sld
    use def_solidz,          only : SLD_IMPLICIT_SCHEME 

    implicit none

    integer(ip), intent(in) :: ipoin               !< Point number
    integer(ip), intent(in) :: ndofn               !< Dimensions
    integer(ip)             :: ibopo, itott

    ibopo = lpoty(ipoin)
    itott = (ipoin-1)*ndofn

    if ( ibopo > 0 ) then
       !
       ! Local Axes
       !
       if ( kfl_fixrs_sld(ipoin) /= 0_ip   .and. &
            kfl_fixno_sld(1,ipoin) /= 3_ip .and. &
            kfl_timet_sld == SLD_IMPLICIT_SCHEME ) then
          !
          ! User local axes (Local --> Global)
          !
          call sld_csys_rotuni(2_ip,ndofn,ipoin,unkno(itott+1:itott+ndofn))
          call sld_csys_rotuni(2_ip,ndofn,ipoin,displ(1:ndofn,ipoin,ITER_K))

       end if
    end if

  end subroutine sld_bclocal_displ

  !-----------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    April, 2019
  !> @brief   RHS rotation and fixity due to local axes
  !> @details
  !>          RHS is fixed when local and then is rotated back.
  !> @note    Rotation is performed in sld_elmmat subroutine
  !>          Rotations: GLOBAL -> LOCAL
  !-----------------------------------------------------------------------

  subroutine sld_bclocal_rhsid_fixity()

    use def_master, only : rhsid
    use def_domain, only : npoin, ndime

    implicit none

    integer(ip) :: ipoin, ibopo, itott

    do ipoin = 1,npoin
       itott = (ipoin-1) * ndime
       ibopo = lpoty(ipoin)

       if ( ibopo > 0 ) then

          if ( kfl_fixno_sld(1,ipoin) == 2_ip .or. kfl_fixno_sld(1,ipoin) == 3_ip ) then
             !
             ! PDN Contact: Fix normal in local CSYS and rotate back
             !
             rhsid(itott+1) = 0.0_rp                                            ! fix normal value
             call sld_csys_rotuni(2_ip,ndime,ipoin,rhsid(itott+1:itott+ndime))  ! rotate back

          end if

       end if
    end do

  end subroutine sld_bclocal_rhsid_fixity

  !-----------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    April, 2019
  !> @brief   Rotation back unknown after solve equation
  !> @details Rotation back unknown after solve equation
  !>          Rotations: LOCAL -> GLOBAL
  !-----------------------------------------------------------------------

  subroutine sld_bclocal_roback(ndofn,ipoin)

    use def_solidz, only : dunkn_sld

    implicit none

    integer(ip), intent(in) :: ipoin         !< Point number
    integer(ip), intent(in) :: ndofn         !< Dimensions
    integer(ip)             :: itott,ibopo

    itott = (ipoin-1) * ndofn
    ibopo = lpoty(ipoin)

    if ( ibopo > 0 ) then
       !
       ! Local Axes (Local --> Global)
       !
       if ( kfl_fixno_sld(1,ipoin) == 2_ip .or. &
            kfl_fixno_sld(1,ipoin) == 3_ip .or. kfl_fixrs_sld(ipoin) /= 0_ip ) then

          call sld_csys_rotuni(2_ip,ndofn,ipoin,dunkn_sld(itott+1:itott+ndofn))  ! rotate back

       end if

    end if

  end subroutine sld_bclocal_roback

  !-----------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    April, 2019
  !> @brief   Rotation elemental rhs and matrix
  !> @details Rotation elemental rhs and matrix
  !>          Rotations: GLOBAL -> LOCAL
  !-----------------------------------------------------------------------

  subroutine sld_bclocal_elmat_and_elrhs(ndofn,pnode,lnods,elrhs,elmat)

    implicit none

    integer(ip), intent(in)    :: ndofn                          !< Dimensions
    integer(ip), intent(in)    :: pnode                          !< Number of nodes for element
    integer(ip), intent(in)    :: lnods(pnode)                   !< List nodes element
    real(rp),    intent(inout) :: elrhs(pnode*ndofn)             !< Elemental RHS
    real(rp),    optional      :: elmat(pnode*ndofn,pnode*ndofn) !< Elemental matrix

    integer(ip)                :: inode,ipoin,ibopo,itott
    !
    ! Local bases: Global --> Local
    !
    do inode = 1, pnode
       ipoin = lnods(inode)
       ibopo = lpoty(ipoin)
       itott = (inode - 1)*ndofn
       if (ibopo > 0) then
          if ( kfl_fixno_sld(1,ipoin) == 2_ip .or. &
               kfl_fixno_sld(1,ipoin) == 3_ip .or. kfl_fixrs_sld(ipoin) /= 0_ip ) then
             !
             ! Rotate ELRHS and ELMAT (if necessary)
             !
             call sld_csys_rotuni(1_ip,ndofn,ipoin,elrhs(itott+1:itott+ndofn))

             if ( present(elmat) ) then

                call runend('MOD_SLD_BCLOCAL: ROTATION ELMAT NOT IMPLEMENTED')

             end if

          end if
       end if
    end do

  end subroutine sld_bclocal_elmat_and_elrhs

end module mod_sld_bclocal
