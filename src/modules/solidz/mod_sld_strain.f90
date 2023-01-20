!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-------------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @authors Constantine Butakoff      : cbutakoff@elem.bio
!> @date    Febreuary, 2021
!> @}
!------------------------------------------------------------------------------
module mod_sld_strain
  use def_kintyp,   only                   :  ip, rp, lg
  use def_domain,   only                   :  ndime
  ! ----------------------------------------
  implicit none

  type strain_directions
    logical(lg)                          :: on
    integer(ip)                          :: field_lng
    integer(ip)                          :: field_rad
    integer(ip)                          :: field_cir
    real(rp),   pointer                  :: lng(:,:)
    real(rp),   pointer                  :: rad(:,:)
    real(rp),   pointer                  :: cir(:,:)
    contains
      procedure,   pass                 :: init => init_strain_basis
      procedure,   pass                 :: send_data => strains_send_data
      procedure,   pass                 :: assign_fields => strain_point_fields

  end type

  type(strain_directions), public         :: strain_basis


contains

  ! ----------------------------------------
  !
  ! Public methods
  !
  ! ----------------------------------------
  subroutine init_strain_basis(self)
    implicit none
    class(strain_directions), intent(inout) :: self

    self % on = .false.
    self % field_lng = 0_ip
    self % field_rad = 0_ip
    self % field_cir = 0_ip   
    self % lng => null()
    self % rad => null()
    self % cir => null() 
    ! -------------------------------
  end subroutine init_strain_basis

  subroutine strains_send_data( self, end_exchange )
    ! -------------------------------
    use def_master,            only : ISEQUEN
    use mod_exchange,          only : exchange_end, exchange_add
    ! -------------------------------
    implicit none
    ! -------------------------------
    class(strain_directions), intent(inout) :: self
    logical(lg), intent(in), optional       :: end_exchange
    ! -------------------------------
    ! Only continue if it is parall
    if( ISEQUEN ) return

    ! Interchange initialisation 
    call exchange_add( self % on )
    call exchange_add( self % field_lng )
    call exchange_add( self % field_rad )
    call exchange_add( self % field_cir )
    if (present(end_exchange)) then
      if (end_exchange) then
        call exchange_end()
      end if
    end if
    ! -------------------------------
  end subroutine strains_send_data


  subroutine strain_point_fields( self )
    use def_domain, only : xfiel
    implicit none
    ! -------------------------------
    class(strain_directions), intent(inout) :: self
    ! -------------------------------

    if( self % on )then
      ! Associate fiber direction to the corresponding field
      self % lng => xfiel(self % field_lng) % a(:,:,1)
      self % rad => xfiel(self % field_rad) % a(:,:,1)
      self % cir => xfiel(self % field_cir) % a(:,:,1)
    endif
  end subroutine strain_point_fields


end module 
