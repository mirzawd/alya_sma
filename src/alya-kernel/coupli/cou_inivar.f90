!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Coupling
!> @{
!> @file    ker_inivar.f90
!> @author  Guillaume Houzeaux
!> @date    08/10/2014
!> @brief   Initialize coupling variables
!> @details Initialize coupling variables
!> @} 
!-----------------------------------------------------------------------

subroutine cou_inivar(itask)
  use def_kintyp
  use def_master
  use def_kermod
  use def_inpout
  use def_domain
  use def_coupli
  use mod_parall,         only : I_AM_IN_COLOR
  use mod_parall,         only : PAR_COMM_COLOR
  use mod_couplings,      only : MIRROR_COUPLING
  use mod_couplings,      only : I_AM_IN_COUPLING
  implicit none
  integer(ip), intent(in) :: itask 
  integer(ip)             :: icoup,ipoin,kpoin,ipass
  integer(ip)             :: color_target,color_source
  real(rp)                :: xvalu
  logical(lg)             :: ifimplicit

  select case ( itask )

  case ( 1_ip )

     !-------------------------------------------------------------------
     !
     ! Initialize variables before reading data
     !
     !-------------------------------------------------------------------

     kfl_gozon       = 1
     nscal_cou       = 0               ! Number of scalar couplings
     nullify(scala_cou)                ! Scalar values to broadcast
     !
     ! Boundaries (useful if we have holes)
     !
     nboun_cou       =  nboun
     lnodb_cou       => lnodb
     ltypb_cou       => ltypb
     lboch_cou       => lboch
     lnnob_cou       => lnnob
     lboel_cou       => lboel
     lelbo_cou       => lelbo

     !-------------------------------------------------------------------
     !
     ! Initialize variables after reading data
     !
     !-------------------------------------------------------------------

     do icoup = 1,mcoup

        color_target = coupling_type(icoup) % color_target
        color_source = coupling_type(icoup) % color_source        
        !
        ! Define kind of coupling
        !
        if( coupling_type(icoup) % zone_target + coupling_type(icoup) % zone_source == 0 ) then
           coupling_type(icoup) % kind = BETWEEN_SUBDOMAINS
        else           
           coupling_type(icoup) % kind = BETWEEN_ZONES
        end if
        !
        ! Define mirror coupling
        !
        coupling_type(icoup) % mirror_coupling = MIRROR_COUPLING(icoup) 
        !
        ! Subdomain coupling
        !
        if( coupling_type(icoup) % kind == BETWEEN_SUBDOMAINS ) then
           if( I_AM_IN_COLOR(color_target) .or. I_AM_IN_COLOR(color_source) ) then
              ncoup_implicit = ncoup_implicit + 1
           end if
        end if
        !
        ! Scalar couplings
        !
        if( coupling_type(icoup) % what == SCALAR ) then
           nscal_cou = nscal_cou + 1
        end if
        !
        ! Communicators
        !
        coupling_type(icoup) % commd % PAR_COMM_WORLD = PAR_COMM_COLOR(color_target,color_source)
     end do

  case ( 2_ip ) 

     !-------------------------------------------------------------------
     !
     ! Initialize variables after initializing coupling
     !
     !-------------------------------------------------------------------

     !
     ! Implicit couplings
     !
     if( mcoup > 0 ) then
        call cou_memory(2_ip) ! MASK_COU
        if( INOTMASTER ) mask_cou = 1.0_rp

        do ipass = 1,2
           ncoup_implicit_d = 0             
           ncoup_implicit_n = 0             
           do icoup = 1,mcoup
              ifimplicit = .false.
              if(    I_AM_IN_COUPLING(icoup)                            .and. &
                   & coupling_type(icoup) % kind == BETWEEN_SUBDOMAINS  .and. &
                   & coupling_type(icoup) % what == UNKNOWN ) then
                 ncoup_implicit_d = ncoup_implicit_d + 1
                 ifimplicit       = .true.
                 xvalu            = 0.0_rp
                 if( ipass == 2) lcoup_implicit_d(ncoup_implicit_d) = icoup
              end if
              if(    I_AM_IN_COUPLING(icoup)                            .and. &
                   & coupling_type(icoup) % kind == BETWEEN_SUBDOMAINS  .and. &
                   & coupling_type(icoup) % what == RESIDUAL ) then
                 ncoup_implicit_n = ncoup_implicit_n + 1
                 ifimplicit       = .true.
                 xvalu            = 1.0_rp
                 if( ipass == 2) lcoup_implicit_n(ncoup_implicit_n) = icoup
              end if
              if( ifimplicit ) then
                 if( ipass == 1 ) then
                    do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
                       ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                       mask_cou(ipoin) = xvalu
                    end do
                 end if
              end if
           end do
           if( ipass == 1 ) call cou_memory(3_ip) 
        end do
     end if
     !
     ! Allocate array for scalar coupling
     !
     call cou_memory(4_ip)
     
  end select

end subroutine cou_inivar
