!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_reabcs()
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_reabcs
  ! NAME
  !    chm_reabcs
  ! DESCRIPTION
  !    This routine reads the boundary conditions
  ! OUTPUT
  ! USES
  ! USED BY
  !    chm_turnon
  !***
  !-----------------------------------------------------------------------
  use def_inpout,              only : words, exists, getint, getrea
  use def_master,              only : INOTSLAVE, intost, speci
  use def_domain,              only : kfl_icodb, kfl_icodn, tbcod, tncod
  use def_kintyp,              only : ip, rp
  use def_chemic,              only : kfl_initi_chm, xinit_chm, tncod_chm, kfl_allcl_chm, kfl_conbc_chm, kfl_fields_scale_chm,&
                                      kfl_model_chm, kfl_solve_cond_CMC_chm, nclas_chm, size_tncod_CMC_chm, tbcod_chm
  use mod_opebcs,              only : opebcs_initialization_structure, opebcs_initialization_variable,&
                                      boundary_conditions_read_node_codes
  use mod_ecoute,              only : ecoute
  use mod_chm_operations_CMC,  only : chm_initial_actions_reabcs_CMC
  implicit none
  integer(ip)  :: iclas,ipara

  external     :: chm_membcs
  external     :: runend
  external     :: reacod

  !
  ! Initialization
  !
  kfl_fields_scale_chm  = 0      ! Default use scale values


  !
  ! Allocate memory
  !
  if( kfl_icodn > 0 ) then
    if (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1) then
       call opebcs_initialization_structure(size_tncod_CMC_chm,tncod_chm)
       do iclas=1,size_tncod_CMC_chm
          call opebcs_initialization_variable ( 1_ip,tncod_chm(iclas))
       end do
    else
       call opebcs_initialization_structure(nclas_chm,tncod_chm)
       do iclas=1,nclas_chm
          call opebcs_initialization_variable ( 1_ip,tncod_chm(iclas))
       end do
    end if
  end if

  if( kfl_icodb > 0 ) then
    call opebcs_initialization_structure(nclas_chm,tbcod_chm)
    do iclas=1,nclas_chm
       call opebcs_initialization_variable ( 1_ip,tbcod_chm(iclas))
    end do
     !!AB call opbbcs(0_ip,nclas_chm,1_ip,tbcod_chm)
  end if
  call chm_membcs(3_ip)

  if( INOTSLAVE ) then

     !
     ! Initializations
     !
     kfl_conbc_chm = 1 ! Constant boundary conditions

     kfl_allcl_chm = 0
     do iclas = 1,nclas_chm
        kfl_initi_chm(iclas) = 0
        xinit_chm(iclas,:)     = 0.0_rp
     end do

     call ecoute('chm_reabcs')
     do while( words(1) /= 'BOUND')
        call ecoute('chm_reabcs')
     end do

     if(exists('NONCO') ) then
        kfl_conbc_chm=0
     else
        kfl_conbc_chm=1
     end if

     !
     ! Read data
     !
     call ecoute('chm_reabcs')
     do while( words(1) /= 'ENDBO')

        if( words(1) == 'CODES' .and. exists('NODES') ) then
           !
           ! User-defined codes on nodes
           !
           if( exists('CLASS') ) then
              if( exists('ALL  ') ) then
                 iclas = 1
                 kfl_allcl_chm = 1
              else
                 iclas = getint('CLASS',1_ip,  '#Class number')
                 kfl_allcl_chm = 0
              end if
           else
              call runend('BOUNDARY CODES NEEDS CLASS NUMBER')
           end if

           if (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1) then
              if( iclas < 1 .or. iclas > size_tncod_CMC_chm ) then
                 call runend('Chemic boundary conditions: class number is wrong')
              end if
           else
              if( iclas < 1 .or. iclas > nclas_chm ) then
                 call runend('Chemic boundary conditions: class number '//trim(intost(iclas))//' is wrong')
              end if
           end if

           tncod => tncod_chm(iclas:)
           call boundary_conditions_read_node_codes('CHEMIC')


        else if( words(1) == 'CODES' .and. exists('BOUND') ) then
           !
           ! User-defined codes on boundaries
           !
           if( exists('CLASS') ) then
              iclas = getint('CLASS',1_ip,  '#Class number')
              if( iclas >= 1 .and. iclas <= nclas_chm ) then
                 tbcod => tbcod_chm(iclas:)
                 call reacod(2_ip)
              end if
           else
              call runend('BOUNDARY CODES NEEDS CLASS NUMBER')
           end if

        else if( words(1) == 'FIELD' ) then
           !
           ! Fields provided for unscale Yc (NO corrections)
           !
           if( exists('UNSCA') ) kfl_fields_scale_chm = 0
           !
           ! Fields provided for scale c
           !
           if( exists('START') ) kfl_fields_scale_chm = 1
           !
           ! Non-uniform fields for scale c
           !
           if( exists('SCALE') ) kfl_fields_scale_chm = 2

        else if( words(1) == 'INITI' ) then
           !
           ! Initial conditions
           !
           call ecoute('chm_reabcs')
           do while( words(1) /= 'ENDIN')
              if( words(1) == 'CLASS' ) then

                 if( exists('CONST') ) then

                    if( exists('CLASS') ) then
                       iclas = getint('CLASS',1_ip,  '#Class number')
                       if( exists('EQUIL') ) then
                          kfl_initi_chm(iclas) = 1
                       else
                          kfl_initi_chm(iclas) = 2
                          xinit_chm(iclas,1) = getrea('CONST',0.0_rp,'#Class constant initial value')
                       end if
                       do while( words(1) /= 'ENDCL')
                          call ecoute('chm_reabcs')
                       end do
                    end if

                 end if

              else if( words(1) == 'SPECY' ) then

                 if( exists('CONST') ) then

                    if( exists('SPECY') ) then
                       iclas = getint('SPECY',1_ip,  '#Specy number')
                       kfl_initi_chm(iclas) = 2
                       xinit_chm(iclas,1) = getrea('CONST',0.0_rp,'#Class constant initial value')
                       do while( words(1) /= 'ENDSP')
                          call ecoute('chm_reabcs')
                       end do
                    end if
                 end if

              else if( words(1) == 'CONST' ) then

                 if( exists('ONSET') ) then
                    do iclas = 1,nclas_chm
                       if( exists(speci(iclas)%name(1:5)) ) then
                          kfl_initi_chm(iclas) = -getint('ONSET',-1_ip,'#prescribe value on set only')
                          xinit_chm(iclas,1) = getrea(speci(iclas)%name,0.0_rp,'#Species constant initial value')
                       end if
                    enddo
                 else  if( exists('BIPHA') .or. exists('LEVEL') ) then
                    if (exists('BIPHA')) then
                       do iclas = 1,nclas_chm
                          kfl_initi_chm(iclas) = 4
                       enddo
                    else if (exists('LEVEL')) then
                       do iclas = 1,nclas_chm
                          kfl_initi_chm(iclas) = 3
                       enddo
                    endif
                    do while( words(1) /= 'ENDPH')
                       if (words(1)=='PHASE') then
                          ipara = getint('PHASE',1_ip,'Phase for initial value')
                          do iclas = 1,nclas_chm
                             if( exists(speci(iclas)%name(1:5)) ) then
                                xinit_chm(iclas,ipara) = getrea(speci(iclas)%name,0.0_rp,'#Species constant initial value')
                             end if
                          enddo
                       endif
                       call ecoute('chm_reabcs')
                    end do
                 else
                    do iclas = 1,nclas_chm
                       if( exists(speci(iclas)%name(1:5)) ) then
                          kfl_initi_chm(iclas) = 2
                          xinit_chm(iclas,1) = getrea(speci(iclas)%name,0.0_rp,'#Species constant initial value')
                       end if
                    enddo
                 end if

              end if

              call ecoute('chm_reabcs')
           end do

        else if( words(1) == 'NATUR' ) then

           call runend('CHM_REABCS: NATURAL BOUNDARY CONDITIONS DEPRECTAED, CHECK REVISION BEFORE DEC 2018')

        end if
        call ecoute('chm_reabcs')
     end do

     if (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1) then
        call chm_initial_actions_reabcs_CMC
     end if

  end if

end subroutine chm_reabcs
