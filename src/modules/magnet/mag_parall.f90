!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_parall(itask)
  !-----------------------------------------------------------------------
  !****f* Magnet/mag_parall
  ! NAME
  !    mag_parall
  ! DESCRIPTION
  !    This routine is a bridge to Parall service  
  !    All processes run this routine
  ! USED BY
  !    mag_turnon
  !***
  !-----------------------------------------------------------------------

  use def_master
  use def_magnet 
  use def_domain,             only : ndime
  use mod_opebcs,             only : boundary_conditions_exchange
  use mod_exchange,           only : exchange_init
  use mod_exchange,           only : exchange_add
  use mod_exchange,           only : exchange_end
  use mod_solver,             only : solver_parall
  use mod_output_postprocess, only : output_postprocess_parall   
  implicit none

  integer(ip), intent(in) :: itask

  integer(ip) :: pmate, icomp, idime

  if( ISEQUEN ) return

  select case ( itask )

  case(1_ip)
     !
     ! Exchange data read in mag_reaphy, mag_reanut and mag_reaous
     ! always using MPI, even if this is a partition restart
     !
     call exchange_init()
     !
     ! Read in mag_reanut
     !
     call exchange_add(dtmin_mag)
     call exchange_add(dtmax_mag)
     call exchange_add(theta_mag)
     call exchange_add(nltol_mag)
     call exchange_add(reltol_mag)
     call exchange_add(abstol_mag)

     call exchange_add(nlite_mag)
     call exchange_add(nlide_mag)
     call exchange_add(gslin_mag)
     call exchange_add(gstri_mag)
     call exchange_add(gsqua_mag)
     call exchange_add(gstet_mag)
     call exchange_add(gshex_mag)
     call exchange_add(bdfode_mag % s)
     call exchange_add(struct_mag)
     call exchange_add(postev_mag)
     call exchange_add(postce_mag)
     !
     ! Read in mag_reaous
     !
     call exchange_add(kfl_nrj_mag)
     call exchange_add(kfl_mtz_mag)
     call exchange_add(kfl_crn_mag)
     call exchange_add(kfl_vlm_mag)
     !
     ! Read in mag_reaphy
     !
     call exchange_add(kfl_lagr_mag)
     call exchange_add(kfl_self_mag)
     call exchange_add(kfl_prev_mag)

     call exchange_add(nintp1_mag)



     call exchange_add(kfl_axsym_mag)
     call exchange_add(constr_total)

     do pmate = 1_ip, maxmat_mag
        call exchange_add(Ec0_mag(pmate))
        call exchange_add(Jc0_mag(pmate))
        call exchange_add(nc0_mag(pmate))
        call exchange_add(mur_mag(pmate))
        call exchange_add(rho_mag(pmate))
        call exchange_add(B0k_mag(pmate))

        call exchange_add(constrlist_mag(pmate))
        call exchange_add(resistOpt_mag(pmate))
        call exchange_add(scalinOpt_mag(pmate))

        call exchange_add(selfList_mag(pmate))

        call exchange_add(kfl_ncriso_mag(pmate))
        call exchange_add(kfl_rhoiso_mag(pmate))
        call exchange_add(kfl_Ecriso_mag(pmate))
        call exchange_add(kfl_Jc0iso_mag(pmate))
        call exchange_add(kfl_muriso_mag(pmate))
        call exchange_add(kfl_Bc0iso_mag(pmate))
        call exchange_add(kfl_Tc0iso_mag(pmate))

        do icomp = 1_ip, ncomp_mag
           call exchange_add(resmat_mag(icomp, pmate))
           call exchange_add(Jcrmat_mag(icomp, pmate))

           call exchange_add(rhomat_mag(icomp, pmate)) 
           call exchange_add(Ecrmat_mag(icomp, pmate))
           call exchange_add(ncrmat_mag(icomp, pmate))
           call exchange_add(murmat_mag(icomp, pmate))
           call exchange_add(Jc0mat_mag(icomp, pmate))
           call exchange_add(Bc0mat_mag(icomp, pmate))
           call exchange_add(Tc0mat_mag(icomp, pmate))
        end do

        do idime = 1_ip, ndime
           call exchange_add(momori_mag(idime, pmate))
        end do
     end do
     !
     ! Solver and postprocess
     !
     call solver_parall()
     call output_postprocess_parall()

     call exchange_end()
     
  case(2_ip)

     call boundary_conditions_exchange(tecod_mag)

  end select


end subroutine mag_parall
