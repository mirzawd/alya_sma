!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_stabil.f90
!> @author  Gerard Guillamet
!> @date    September, 2018
!> @brief
!> @details
!>          \verbatim
!>
!>          1. Material dependent Rayleigh Damping
!>          2. Automatic stabilization with constant damping for static problems
!>
!>          \endverbatim
!>
!> @}
!-----------------------------------------------------------------------

subroutine sld_stabil(itask,pnode,pevat,pgaus,pmate,gpvol,gpsha,&
     eldis,elvel,elmuu,elsti,elfda,elrhs)

  use def_kintyp, only : ip, rp
  use def_master, only : ITER_K, TIME_N, dtime
  use def_domain, only : ndime
  use def_solidz, only : kfl_timei_sld
  use def_solidz, only : SLD_IMPLICIT_SCHEME, SLD_EXPLICIT_SCHEME
  use def_solidz, only : SLD_STATIC_PROBLEM, SLD_DYNAMIC_PROBLEM
  use def_solidz, only : ncomp_sld
  use def_solidz, only : kfl_dampi_sld, dampi_sld
  use def_solidz, only : kfl_stabi_sld, dafac_sld

  implicit none
  !
  ! Input and output variables
  !
  integer(ip), intent(in)    :: itask                        !< What to do
  integer(ip), intent(in)    :: pnode                        !< Number of element nodes
  integer(ip), intent(in)    :: pevat                        !< Number of element dof
  integer(ip), intent(in)    :: pgaus                        !< Number of Gauss points
  integer(ip), intent(in)    :: pmate                        !< Material code
  real(rp),    intent(in)    :: gpvol(pgaus)                 !< w*|J|, |J|
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)           !< Nk
  real(rp),    intent(in)    :: eldis(ndime,pnode,ncomp_sld) !< Element node displacement
  real(rp),    intent(in)    :: elvel(ndime,pnode,ncomp_sld) !< Element node displacement
  real(rp),    intent(in)    :: elmuu(pnode)                 !< Element lumped mass matrix
  real(rp),    intent(in)    :: elsti(pevat,pevat)           !< Element stiffness matrix
  real(rp),    intent(inout) :: elrhs(pevat)                 !< Elemental total RHS
  real(rp),    intent(out)   :: elfda(pevat)                 !< Elemental damping force
  !
  ! Indices and dimensions
  !
  integer(ip)    :: inode,jnode
  integer(ip)    :: idime,jdime
  integer(ip)    :: igaus
  integer(ip)    :: ievat,jevat
  !
  ! Element matrices and vectors
  !
  real(rp)       :: elral(pevat)
  real(rp)       :: elrbe(pevat)
  real(rp)       :: eluma(pnode,pnode)
  !
  ! Gather
  !
  real(rp)       :: elvel_static(ndime,pnode)
  !
  ! Other
  !
  real(rp)       :: volux
  !
  ! Intitializations
  !
  elfda(:) = 0.0_rp

  !------------------------------------------------------------------------
  !
  ! Material dependent Rayleigh damping
  !
  ! [C] = Alpha*[M] + Beta*[K]
  !------------------------------------------------------------------------

  if ( kfl_dampi_sld(pmate) == 1_ip .and. kfl_timei_sld == SLD_DYNAMIC_PROBLEM ) then
     !
     ! Initializations
     !
     elral(:) = 0.0_rp
     elrbe(:) = 0.0_rp
     !
     ! Mass matrix term Alpha*[M]
     !
     do inode = 1, pnode
        do idime = 1, ndime
           ievat = ( inode - 1 ) * ndime + idime
           elral(ievat) = elral(ievat) &
                + dampi_sld(1,pmate) &
                * elmuu(inode) &
                * elvel(idime,inode,ITER_K)
        end do
     end do

     if ( itask == SLD_IMPLICIT_SCHEME ) then
        !
        ! Stiffness matrix term Beta*[K]
        !
        do inode = 1, pnode
           do idime = 1, ndime
              ievat = ( inode - 1 ) * ndime + idime
              do jnode = 1, pnode
                 do jdime = 1, ndime
                    jevat = ( inode - 1 ) * ndime + jdime
                    elrbe(ievat) = elrbe(ievat) &
                         + dampi_sld(2,pmate) &
                         * elsti(ievat,jevat) &
                         * elvel(idime,inode,ITER_K)
                 end do
              end do
           end do
        end do
     end if
     !
     ! [RHS] = [RHS] - [Fdam]
     !
     elfda(1:pevat) = elral(1:pevat) + elrbe(1:pevat)
     elrhs(1:pevat) = elrhs(1:pevat) - elfda(1:pevat)

  end if

  !------------------------------------------------------------------------
  !
  ! Automatic stabilization with constant damping factor
  !
  ! [RHS] = [Fext] - [Fint] - [Fart]
  ! where,
  !  [Fda] = c*[M]*[v] v = \delta{u} / \delta{t}
  !------------------------------------------------------------------------

  if ( kfl_stabi_sld == 1_ip .and. itask == SLD_IMPLICIT_SCHEME &
       .and. kfl_timei_sld == SLD_STATIC_PROBLEM ) then

     !
     ! Elemental vector of nodal velocities
     !
     elvel_static(1:ndime,1:pnode) = (eldis(1:ndime,1:pnode,ITER_K) &
          - eldis(1:ndime,1:pnode,TIME_N))/dtime
     !
     ! Artificial mass matrix with unity density
     !
     eluma(:,:) = 0.0_rp
     do igaus = 1, pgaus
        volux = gpvol(igaus)
        do inode = 1, pnode
           do jnode = 1, pnode
              !
              ! Unity mass matrix
              !   M_ab = int_\Omega * Na * Nb d\Omega
              eluma(inode,jnode) = eluma(inode,jnode) &
                   + volux*gpsha(inode,igaus)*gpsha(jnode,igaus)
           end do
        end do
     end do
     !
     ! Elemental viscous force vector
     !
     do inode = 1, pnode
        do jnode = 1, pnode
           do idime = 1, ndime
              ievat = ( inode - 1 ) * ndime + idime
              elfda(ievat) = elfda(ievat) &
                   + dafac_sld*eluma(inode,jnode)*elvel_static(idime,jnode)
           end do
        end do
     end do
     !
     ! [RHS] = [RHS] - [Fdam]
     !
     elrhs(1:pevat) = elrhs(1:pevat) - elfda(1:pevat)

  end if

end subroutine sld_stabil
