!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_iniunk.f90
!> @author  Mariano Vazquez
!> @author  Eva Casoni
!> @todo    To finish X-FEM
!> @todo    Revise initial conditions for displacement and stress
!> @author  Gerard Guillamet
!> @date    July, 2017-Added velocities
!> @date    January, 2022-Added states
!> @date    November, 2022-Added temperature
!> @brief   Initial conditions for solidz module
!> @details Set up the initial conditions for displacement and velocity.
!>          If this is a restart, initial conditions are read from files.
!>          The different options without a restart are:
!>          \verbatim
!>
!>          Velocity
!>          --------
!>          KFL_INVEL_SLD = 0 ... Deactivated
!>                        > 0 ... Values from CODES fields
!>                        < 0 ... User values from input file
!>
!>          State dependent variables
!>          --------
!>          KFL_INSDV_SLD = 0 ... Deactivated
!>                        > 0 ... Values from CODES fields
!>
!>          Temperature
!>          --------
!>          KFL_INTEM_SLD = 0 ... Deactivated
!>                        > 0 ... Values from CODES fields
!>                        < 0 ... User values from input file
!>          \endverbatim
!>
!>          Values are stored in position:
!>          \verbatim
!>          VELOC_SLD(1:NDIME,1:NPOIN,NCOMP_SLD)
!>          SVEGM_SLD(NELEM) % a(NSVAR_SLD,PGAUS,2)
!>          THERM(1:NPOIN,ITER_K) 
!>          \endverbatim
!> @}
!-----------------------------------------------------------------------

subroutine sld_iniunk()

  use def_kintyp,  only : ip,rp
  use def_master,  only : INOTMASTER,ITASK_INIUNK
  use def_master,  only : ITER_K,ITER_K_STATE
  use def_master,  only : mem_modul,modul
  use def_master,  only : kfl_rstar,nfacg
  use def_master,  only : displ
  use def_domain,  only : ndime,npoin,nelem,xfiel
  use def_domain,  only : ltype,lmate,ngaus
  use def_domain,  only : kfl_elcoh
  use mod_cutele,  only : cutele
  use mod_memory,  only : memory_alloca
  use mod_chktyp,  only : check_type
  use def_solidz,  only : veloc_sld,accel_sld,svegm_sld
  use def_solidz,  only : kfl_sdvar_sld
  use def_solidz,  only : kfl_invel_sld, invel_sld
  use def_solidz,  only : kfl_indis_sld
  use def_solidz,  only : kfl_insdv_sld
  use def_solidz,  only : cockf_sld,lcrkf_sld,crtip_sld,cranx_sld,crapx_sld
  use def_solidz,  only : bvess_sld
  use def_solidz,  only : lawch_sld, lawst_sld
  use def_solidz,  only : kfl_xfeme_sld
  use def_solidz,  only : kfl_rigid_sld, kfl_timei_sld, SLD_DYNAMIC_PROBLEM
  use mod_sld_rbo, only : sld_rbo_iniunk, sld_rbo_updunk
  use mod_sld_atm, only : sld_atm_iniunk
  use mod_sld_atm, only : kfl_intem_sld
  
  implicit none

  integer(ip)          :: ipoin,ielem,idime,isdva
  integer(ip)          :: pgaus,pelty

  if ( kfl_rstar == 0 ) then

     !-------------------------------------------------------------------
     !
     ! Normal run
     !
     !-------------------------------------------------------------------

     if ( kfl_rigid_sld == 0_ip ) then

        !----------------------------------------------------------------
        !
        ! Deformable body
        !
        !----------------------------------------------------------------

        if( INOTMASTER ) then

           if( kfl_indis_sld(1) < 0 ) then
              !
              ! Displacement for Pre-stress
              !
              do ipoin = 1,npoin
                 do idime=1,ndime
                    If (kfl_indis_sld(1) < 0) then ! an initial displacement is given as a field
                       displ(idime,ipoin,1)= xfiel(-kfl_indis_sld(1)) % a(idime,ipoin,1)
                    else if (kfl_indis_sld(1) == 0) then
                       ! 
                       displ(idime,ipoin,1)=bvess_sld(idime,ipoin,1)
                    end if
                    veloc_sld(idime,ipoin,1) = 0.0_rp
                    accel_sld(idime,ipoin,1) = 0.0_rp
                 end do
              end do
           end if

           if( kfl_invel_sld /= 0_ip .and. kfl_timei_sld == SLD_DYNAMIC_PROBLEM ) then
              !
              ! Velocity
              !
              if( kfl_invel_sld > 0 ) then
                 do ipoin = 1,npoin
                    veloc_sld(1:ndime,ipoin,ITER_K) = xfiel(kfl_invel_sld) % a(1:ndime,ipoin,1)
                 end do
              else
                 do ipoin = 1,npoin
                    veloc_sld(1:ndime,ipoin,ITER_K) = invel_sld(1:ndime)
                 end do
              end if

           end if

           if( kfl_insdv_sld /= 0_ip .and. kfl_sdvar_sld == 1_ip ) then
              !
              ! State dependent variables
              !
              !-----------------------------------------------------
              !$OMP PARALLEL DO                                    &
              !$OMP SCHEDULE (STATIC)                              &
              !$OMP DEFAULT  ( NONE )                              &
              !$OMP PRIVATE  ( ielem, pelty, pgaus, isdva )        &
              !$OMP SHARED   ( nelem, ltype, ngaus, lmate, xfiel,  &
              !$OMP            lawch_sld, lawst_sld, svegm_sld,    &
              !$OMP            kfl_insdv_sld                       )
              !-----------------------------------------------------
              do ielem = 1,nelem
                 pelty = ltype(ielem)
                 pgaus = ngaus(pelty)
                 if(      lawch_sld(lmate(ielem)) == 904_ip ) then
                    ! coh904 (d)
                    svegm_sld(ielem)%a(1,1:pgaus,ITER_K_STATE) = xfiel(kfl_insdv_sld) % a(1,ielem,1)
                 else if( lawch_sld(lmate(ielem)) == 905_ip ) then
                    ! coh905 (d,r)
                    do isdva = 1,2
                       svegm_sld(ielem)%a(isdva,1:pgaus,ITER_K_STATE) = xfiel(kfl_insdv_sld) % a(isdva,ielem,1)
                    end do
                 else if( lawst_sld(lmate(ielem)) == 154_ip ) then
                    ! sm154 (d1,d2,d3,d4,d5,d6)
                    do isdva = 1,6
                       svegm_sld(ielem)%a(isdva+4,1:pgaus,ITER_K_STATE) = xfiel(kfl_insdv_sld) % a(isdva,ielem,1)
                    end do
                 end if
              end do
              !$OMP END PARALLEL DO
              !-----------------------------------------------------

           end if

           if( kfl_intem_sld /= 0_ip ) then
              !
              ! Temperature
              !
              call sld_atm_iniunk()
            
           end if
           
           if ( kfl_xfeme_sld > 0 ) then
              !
              ! X-FEM
              !
              ! Allocate memory
              call memory_alloca(mem_modul(1:2,modul),'LCRKF_SLD','sld_iniunk',lcrkf_sld,nfacg)
              call memory_alloca(mem_modul(1:2,modul),'COCKF_SLD','sld_iniunk',cockf_sld,ndime,nfacg)
              call memory_alloca(mem_modul(1:2,modul),'CRTIP_SLD','sld_iniunk',crtip_sld,ndime,nfacg)
              !
              ! XFEM: initialize crapx: element center and cranx: preexisting crack
              !
              call sld_craini(1_ip)
              call sld_enrich(1_ip)
              !
              ! If the element is CUT
              !
              call cutele(1_ip,cranx_sld,crapx_sld,lcrkf_sld)
              !
              ! Initialize internal variables of cohesive laws and contact/friction
              !
              call sld_updcoh(0_ip)
           end if

           if ( kfl_elcoh > 0 ) then
              !
              ! Initialize vector of nodes where the cohesive law is activated
              !
              call sld_updcoh(0_ip)
           end if

        end if

        call sld_updunk(ITASK_INIUNK)

     else if ( kfl_rigid_sld == 1_ip ) then

        !----------------------------------------------------------------
        !
        ! Rigid body
        !
        !----------------------------------------------------------------

        call sld_rbo_iniunk()

        call sld_rbo_updunk(ITASK_INIUNK)

     end if

  end if
  !
  ! Coupling initializations
  !
  call sld_coupli(ITASK_INIUNK)

end subroutine sld_iniunk
