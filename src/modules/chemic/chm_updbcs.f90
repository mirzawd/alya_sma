!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_updbcs(itask)
  !-----------------------------------------------------------------------
  !****f* partis/chm_updbcs
  ! NAME
  !    chm_updbcs
  ! DESCRIPTION
  !    This routine prepares a new time step
  ! USES
  !    chm_updunk
  ! USED BY
  !    partis
  !***
  !-----------------------------------------------------------------------
  use def_master,                  only : IMASTER, INOTEMPTY, ITASK_BEGITE, ITASK_INIUNK, FUNCTION_SPACE_TIME, cutim, conce
  use def_domain,                  only : ndime, npoin, coord
  use def_kintyp,                  only : ip, rp
  use def_chemic,                  only : nclas_chm, imixf_rk, kfl_conbc_chm, kfl_tiacc_chm, kfl_timei_chm, kfl_tisch_chm,&
                                          ncomp_chm, nvar_CMC_chm, kfl_funno_chm, kfl_funtn_chm, bvess_chm, kfl_fixno_chm,&
                                          bvess_CMC_chm
  use mod_ker_space_time_function, only: ker_space_time_function
  use mod_chm_operations_CMC,      only: chm_get_cond_fields_bc_Dir_CMC

  implicit none
  integer(ip), intent(in) :: itask
  real(rp)                :: valold, funval(nclas_chm)
  integer(ip)             :: ipoin,iclas,ifunc,itype

  if( IMASTER ) return

  select case ( itask )

  case( ITASK_INIUNK , ITASK_BEGITE )

     !-------------------------------------------------------------------
     !
     ! Before a global iteration
     !
     !-------------------------------------------------------------------

     if( INOTEMPTY .and. kfl_conbc_chm == 0 ) then
        !
        ! Space/Time functions
        !
        do iclas = 1,nclas_chm
           do ipoin = 1,npoin

              ifunc = kfl_funno_chm(ipoin,iclas)
              itype = kfl_funtn_chm(ipoin,iclas)

              if( itype == FUNCTION_SPACE_TIME ) then
                 !
                 ! Space time function
                 !
                 call ker_space_time_function(&
                      ifunc,coord(1,ipoin),coord(2,ipoin),coord(ndime,ipoin),cutim,funval)
                 if( kfl_timei_chm /= 0 .and. kfl_tiacc_chm == 2 .and. kfl_tisch_chm == 1 ) then
                    valold = conce(iclas,ipoin,ncomp_chm)
                    bvess_chm(iclas,ipoin) = 0.50_rp*(funval(iclas)+valold)
                 else
                    bvess_chm(iclas,ipoin) = funval(iclas)
                 end if
              end if
           end do
        end do
     endif
     !
     ! Impose Dirichlet bc
     !
     do iclas = 1,nclas_chm
        do ipoin = 1,npoin
           if( kfl_fixno_chm(iclas,ipoin) > 0 ) then
              conce(ipoin,iclas,1) = bvess_chm(iclas,ipoin)
           end if
        end do
     end do


  case( 1_ip )
     !-------------------------------------------------------------------
     !
     ! Compute current Dirichlet boundary conditions (CMC model)
     !
     !-------------------------------------------------------------------

     call chm_get_cond_fields_bc_Dir_CMC


  case( 2_ip )
     !-------------------------------------------------------------------
     !
     ! Update bvess_chm for the current mixture fraction level (CMC model)
     !
     !-------------------------------------------------------------------

     do iclas = 1,nvar_CMC_chm
        do ipoin = 1,npoin
           bvess_chm(iclas,ipoin) = bvess_CMC_chm(imixf_rk,ipoin,iclas)
        end do
     end do
  end select

end subroutine chm_updbcs
