!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_updbcs.f90
!> @author  Solidz
!> @date
!> @brief   Boundary conditions update
!> @details This routine updates the boundary conditions:
!>          0. At the beginning of the run
!>          1. Before a time step begins
!>          2. Before a global iteration begins
!>          3. Before an inner iteration begins
!>
!> @}
!-----------------------------------------------------------------------

subroutine sld_updbcs(itask)

  use def_kintyp,                  only : ip, rp
  use def_master,                  only : INOTEMPTY
  use def_master,                  only : FUNCTION_TIME
  use def_master,                  only : FUNCTION_MODULE
  use def_master,                  only : FUNCTION_SPACE_TIME
  use def_master,                  only : FUNCTION_DISCRETE
  use def_master,                  only : cutim
  use def_master,                  only : ITER_K, ITER_AUX, TIME_N
  use def_master,                  only : displ
  use def_domain,                  only : ndime, npoin, nboun
  use def_domain,                  only : lpoty, coord
  use mod_ker_space_time_function, only : ker_space_time_function
  use def_solidz,                  only : bvess_sld, bvnat_sld
  use def_solidz,                  only : kfl_conbc_sld, kfl_funno_sld, kfl_funbo_sld
  use def_solidz,                  only : kfl_funtn_sld, kfl_funtb_sld
  use def_solidz,                  only : kfl_conta_stent
  use def_solidz,                  only : bodyf_sld, kfl_bodyf_sld
  use mod_sld_stent,               only : sld_calculation_contact_stent
  use def_master,                  only : ITASK_TURNON, ITASK_BEGSTE
  use def_master,                  only : ITASK_BEGITE, ITASK_INNITE
  use def_solidz,                  only : kfl_fixno_sld, kfl_fixbo_sld
  use def_solidz,                  only : kfl_contn_stent
  use mod_ker_functions,           only : ker_functions
  use mod_sld_stent,               only : sld_set_boundaries

  implicit none

  external                :: sld_funbou

  integer(ip), intent(in) :: itask      !< What to do
  integer(ip)             :: ipoin,ibopo,ifunc,itype,idime,iboun
  real(rp)                :: dinew(ndime),fonew(ndime),bvnew(ndime)

  select case (itask)

  case( ITASK_TURNON )

  case( ITASK_BEGSTE )

     !-------------------------------------------------------------------
     !
     ! Before a time step
     !
     !-------------------------------------------------------------------

 
     if ( INOTEMPTY ) then
        !
        ! Update boundary conditions (Transient)
        !
        if ( kfl_conbc_sld == 0 ) then ! Non-constant (Transient)

           do ipoin = 1, npoin
              ibopo = lpoty(ipoin)

              dinew(:) = 0.0_rp
              fonew(:) = 0.0_rp
              ifunc    = kfl_funno_sld(ipoin)
              itype    = kfl_funtn_sld(ipoin)

              if ( itype == FUNCTION_MODULE ) then
                 !
                 ! Solidz Functions
                 !
                 ! Displacements
                 call sld_funbou( 2_ip, ipoin, dinew)
                 !
                 ! Forces
                 if ( kfl_bodyf_sld > 0 ) then
                    call sld_funbou( 3_ip, ipoin,fonew)
                    bodyf_sld(1:ndime) = fonew(1:ndime)
                 end if

                 else if( itype == FUNCTION_SPACE_TIME .or. itype == FUNCTION_DISCRETE ) then
                    !
                    ! Space/time or discrete functions
                    !
                    call ker_functions(ipoin,ifunc,itype,bvess_sld(:,ipoin,2),dinew,WHEREIN='ON NODES')
                       
              end if
              !
              ! Update Dircihlet displacement
              !
              do idime = 1,ndime
                 if( kfl_fixno_sld(idime,ipoin) == 1_ip ) then
                    bvess_sld(idime,ipoin,ITER_K) = dinew(idime)
                 end if
              end do
              !
              ! STent computations
              !
 !             if (kfl_conta_stent /= 0_ip) then
                   if (kfl_funno_sld(ipoin) < 0 .and. kfl_contn_stent(ipoin)==1_ip) then
                      ifunc = -kfl_funno_sld(ipoin)
                          
                          call ker_space_time_function(&
                              ifunc,coord(1,ipoin),coord(2,ipoin)+displ(2,ipoin,TIME_N),coord(ndime,ipoin),cutim,dinew(1:ndime))
                          call sld_calculation_contact_stent(ipoin,ibopo,ifunc,dinew(:))

                       if (kfl_conta_stent == 3_ip) then 
                         if ( (coord(2,ipoin)+displ(2,ipoin,TIME_N) <= 0.00_rp) .and. (coord(2,ipoin)+displ(2,ipoin,TIME_N))>-4.0_rp ) then
                            call ker_space_time_function(&
                                ifunc,coord(1,ipoin),coord(2,ipoin)+displ(2,ipoin,TIME_N),coord(ndime,ipoin),cutim,dinew(1:ndime))
                            call sld_calculation_contact_stent(ipoin,ibopo,ifunc,dinew(:))
                         else if ( (coord(2,ipoin)+displ(2,ipoin,TIME_N))<=-4.0_rp ) then
                            call sld_calculation_contact_stent(ipoin,ibopo,ifunc,dinew(:))
                         end if
!                    else if (kfl_conta_stent == 4_ip) then
!                        if ( ((coord(2,ipoin) - cutim*12.5_rp) <= 0.0_rp) .and.  ((coord(2,ipoin) - cutim*12.5_rp) > -4.0_rp) ) then
!                           call ker_space_time_function(&
!                                ifunc,coord(1,ipoin),coord(2,ipoin),coord(ndime,ipoin),cutim,dinew(1:ndime))
!                           call sld_calculation_contact_stent(ipoin,ibopo,ifunc,dinew(:))
!                        else if ( ((coord(2,ipoin) - cutim*12.5_rp) <= -4.0_rp) ) then 
!                           call sld_calculation_contact_stent(ipoin,ibopo,ifunc,dinew(1:ndime)) 
!                        end if
                    end if
                   !  if ( ibopo > 0 ) then
                   !    call sld_calculation_contact_stent(ipoin,ibopo,ifunc,dinew(:))
                   !  else
                   !    call runend('WE ARE TRYING CONTACT WITHOUT A BOUNDARY IBOPO')
                   !  end if
 !               end if
             end if    
           end do
           !
           ! Update Neumann boundary conditions
           !
           do iboun = 1,nboun
              bvnew(:) = 0.0_rp
              if( kfl_fixbo_sld(iboun) /= 0 ) then

                 ifunc = kfl_funbo_sld(iboun)
                 itype = kfl_funtb_sld(iboun)

                 if(      itype == FUNCTION_MODULE ) then
                    !
                    ! Solidz Functions
                    !
                    call sld_funbou(1_ip,iboun,bvnew(:))

                 else if( itype == FUNCTION_SPACE_TIME .or. itype == FUNCTION_DISCRETE ) then
                    !
                    ! Space/time or discrete functions
                    !
                    call ker_functions(iboun,ifunc,itype,bvnat_sld(1:ndime,iboun,ITER_AUX),bvnew,WHEREIN='ON BOUNDARIES')

                 end if
                 
              end if
              !
              ! Update Neumaan values
              !
              bvnat_sld(:,iboun,ITER_K) = bvnew(:)
              
           end do
           
        end if

     end if

  case( ITASK_BEGITE )

     !-------------------------------------------------------------------
     !
     ! Before a global iteration
     !
     !-------------------------------------------------------------------

  case( ITASK_INNITE )

     !-------------------------------------------------------------------
     !
     ! Before an inner iteration
     !
     !-------------------------------------------------------------------

  end select

end subroutine sld_updbcs
