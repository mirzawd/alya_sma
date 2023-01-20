!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_reaphy()

  use def_magnet
  use def_master
  use def_inpout
  use mod_ecoute, only : ecoute

  implicit none

  integer(ip) :: pmate, group

  logical :: group_flg

  selfList_mag = .false.

  constrlist_mag = -1_ip
  constr_total = 0_ip

  group = 0_ip

  if (INOTSLAVE) then

    kfl_edge_mag = 0_ip
    kfl_axsym_mag = .false.

    call ecoute('mag_reaphy')
    do while (words(1) /= 'PHYSI')
      call ecoute('mag_reaphy')
    end do
    call ecoute('mag_reaphy')

    do while (words(1) /= 'ENDPH')

      if (words(1) == 'PROPE') then

        do while (words(1) /= 'ENDPR')

          if (words(1) == 'MATER') then

            pmate = int(param(1),ip)
            if (pmate > maxmat_mag) call runend('mag_reaphy: unknown material')

            do while (words(1) /= 'ENDMA')
              !
              ! Original material properties w/o tensors
              !
              if (words(1) == 'RESIS') then
                !
                ! Resistivity Library
                !
                ! Opt 1 - Constant [CONST]
                !         Parameters: RHO00
                ! Opt 2 - Power Law [POWER]
                !         Parameters: ECRIT, NCRIT, JCRIT 
                !
                if (words(2) == 'CONST') then
                  resistOpt_mag(pmate) = 0_ip
                else if (words(2) == 'POWER') then
                  resistOpt_mag(pmate) = 1_ip
                else 
                  call runend('mag_reaphy: unknown resistivity function (Available: CONST, POWER)')
                end if
              else if (words(1) == 'JCLAW') then
                !
                ! Jc Scaling Law Library
                !
                ! Opt 1 - Constant [CONST]
                !         Parameters: JCRIT
                ! Opt 2 - Kim's Model [KIMMO]
                !         Parameters: JCRIT, B0KIM
                ! Opt 9 - User-defined Model [USERD]
                !         Parameters: JOPT1, JOPT2, JOPT3, JOPT4, JOPT5, JOPT6 (Work in Progress)
                if (words(2) == 'CONST') then
                  scalinOpt_mag(pmate) = 0_ip
                else if (words(2) == 'KIMMO') then
                  scalinOpt_mag(pmate) = 1_ip
                else if (words(2) == 'USERD') then
                  scalinOpt_mag(pmate) = 9_ip
                else
                  call runend('mag_reaphy: unknown scaling law (Available: CONST)')
                end if
!              else if (words(1) == 'RHO00') then
!                rho_mag(pmate) = param(1) 
              else if (words(1) == 'ECRIT') then
                Ec0_mag(pmate) = param(1)
              else if (words(1) == 'NCRIT') then
                nc0_mag(pmate) = param(1)
              else if (words(1) == 'JCRIT') then
                Jc0_mag(pmate) = param(1)
              else if (words(1) == 'MUREL') then
                mur_mag(pmate) = param(1)
              else if (words(1) == 'B0KIM') then
                B0k_mag(pmate) = param(1)

              else if (words(1) == 'RES00') then
                !
                ! Resistivity Library
                !
                ! Opt 1 - Constant [CONST]
                !         Parameters: RHO00
                ! Opt 2 - Power Law [POWER]
                !         Parameters: ECRIT, NCRIT, JCRIT
                !
                kfl_Resiso_mag(pmate) = .true.
                !
                if (words(2) == 'CONST') then
                  resmat_mag(1:3,pmate) = 0_ip
                else if (words(2) == 'POWER') then
                  resmat_mag(1:3,pmate) = 1_ip
                else
                  call runend('mag_reaphy: unknown resistivity function (Available: CONST, POWER)')
                end if
              else if (words(1) == 'JCR00' .and. .not. kfl_Jcriso_mag(pmate)) then
                !
                ! Jc Scaling Law Library
                !
                ! Opt 1 - Constant [CONST]
                !         Parameters: JCRIT
                ! Opt 2 - Kim's Model [KIMMO]
                !         Parameters: JCRIT, B0KIM
                !
                kfl_Jcriso_mag(pmate) = .true.
                !
                if (words(2) == 'CONST') then
                  jcrmat_mag(1:3,pmate) = 0_ip
                else if (words(2) == 'KIMMO') then
                  jcrmat_mag(1:3,pmate) = 1_ip
                else if (words(2) == 'USERD') then
                  jcrmat_mag(1:3,pmate) = 9_ip
                else
                  call runend('mag_reaphy: unknown scaling law (Available: CONST)')
                end if
              else if (words(1) == 'RHO00') then
                kfl_rhoiso_mag(pmate) = .true.
                rhomat_mag(1:3,pmate) = param(1)
              else if (words(1) == 'ECR00') then
                kfl_Ecriso_mag(pmate) = .true.
                Ecrmat_mag(1:3,pmate) = param(1)
              else if (words(1) == 'NCR00') then
                kfl_ncriso_mag(pmate) = .true.
                ncrmat_mag(1:3,pmate) = param(1)
              else if (words(1) == 'JC000') then
                kfl_Jc0iso_mag(pmate) = .true.
                Jc0mat_mag(1:3,pmate) = param(1)
              else if (words(1) == 'MUR00') then
                kfl_muriso_mag(pmate) = .true.
                murmat_mag(1:3,pmate) = param(1)
              else if (words(1) == 'BC000') then
                kfl_Bc0iso_mag(pmate) = .true.
                Bc0mat_mag(1:3,pmate) = param(1)
              else if (words(1) == 'TC000') then
                kfl_Tc0iso_mag(pmate) = .true.
                Tc0mat_mag(1:3,pmate) = param(1)

              !##################################################################
              !
              ! Material properties allowing tensor form
              !
              ! ZZ components
              !
              else if (words(1) == 'RESZZ' .and. .not. kfl_Resiso_mag(pmate)) then
                !
                ! Resistivity Library
                ! Opt 1 - Constant [CONST] - Parameters: RHO
                ! Opt 2 - Power Law [POWER] - Parameters: ECR, NCR, JCR
                !
                if (words(2) == 'CONST') then
                  resmat_mag(3,pmate) = 0_ip
                else if (words(2) == 'POWER') then
                  resmat_mag(3,pmate) = 1_ip
                else
                  call runend('mag_reaphy: unknown resistivity function (Available: CONST, POWER)')
                end if
              else if (words(1) == 'JCRZZ' .and. .not. kfl_Jcriso_mag(pmate)) then
                !
                ! Jc Scaling Law Library
                ! Opt 1 - Constant [CONST] - Parameters: JCR
                ! Opt 2 - Kim's Model [KIMMO] - Parameters: JC0, BC0
                !
                if (words(2) == 'CONST') then
                  jcrmat_mag(3,pmate) = 0_ip
                else if (words(2) == 'KIMMO') then
                  jcrmat_mag(3,pmate) = 1_ip
                else if (words(2) == 'USERD') then
                  jcrmat_mag(3,pmate) = 9_ip
                else
                  call runend('mag_reaphy: unknown scaling law (Available: CONST)')
                end if
              else if (words(1) == 'RHOZZ' .and. .not. kfl_rhoiso_mag(pmate)) then
                rhomat_mag(3,pmate) = param(1)
              else if (words(1) == 'ECRZZ' .and. .not. kfl_Ecriso_mag(pmate)) then
                Ecrmat_mag(3,pmate) = param(1)
              else if (words(1) == 'NCRZZ' .and. .not. kfl_ncriso_mag(pmate)) then
                ncrmat_mag(3,pmate) = param(1)
              else if (words(1) == 'JC0ZZ' .and. .not. kfl_Jc0iso_mag(pmate)) then
                Jc0mat_mag(3,pmate) = param(1)
              else if (words(1) == 'MURZZ' .and. .not. kfl_muriso_mag(pmate)) then
                murmat_mag(3,pmate) = param(1)
              else if (words(1) == 'BC0ZZ' .and. .not. kfl_Bc0iso_mag(pmate)) then
                Bc0mat_mag(3,pmate) = param(1)
              else if (words(1) == 'TC0ZZ' .and. .not. kfl_Tc0iso_mag(pmate)) then
                Tc0mat_mag(3,pmate) = param(1)
              !
              ! XX components
              !
              else if (words(1) == 'RESXX' .and. .not. kfl_Resiso_mag(pmate)) then
                !
                ! Resistivity Library
                ! Opt 1 - Constant [CONST] - Parameters: RHO
                ! Opt 2 - Power Law [POWER] - Parameters: ECR, NCR, JCR
                !
                if (words(2) == 'CONST') then
                  resmat_mag(1,pmate) = 0_ip
                else if (words(2) == 'POWER') then
                  resmat_mag(1,pmate) = 1_ip
                else
                  call runend('mag_reaphy: unknown resistivity function (Available: CONST, POWER)')
                end if
              else if (words(1) == 'JCRXX' .and. .not. kfl_Jcriso_mag(pmate)) then
                !
                ! Jc Scaling Law Library
                ! Opt 1 - Constant [CONST] - Parameters: JCR
                ! Opt 2 - Kim's Model [KIMMO] - Parameters: JC0, BC0
                !
                if (words(2) == 'CONST') then
                  jcrmat_mag(1,pmate) = 0_ip
                else if (words(2) == 'KIMMO') then
                  jcrmat_mag(1,pmate) = 1_ip
                else if (words(2) == 'USERD') then
                  jcrmat_mag(1,pmate) = 9_ip
                else
                  call runend('mag_reaphy: unknown scaling law (Available: CONST, KIMMO)')
                end if
              else if (words(1) == 'RHOXX' .and. .not. kfl_rhoiso_mag(pmate)) then
                rhomat_mag(1,pmate) = param(1)
              else if (words(1) == 'ECRXX' .and. .not. kfl_Ecriso_mag(pmate)) then
                Ecrmat_mag(1,pmate) = param(1)
              else if (words(1) == 'NCRXX' .and. .not. kfl_ncriso_mag(pmate)) then
                ncrmat_mag(1,pmate) = param(1)
              else if (words(1) == 'JC0XX' .and. .not. kfl_Jc0iso_mag(pmate)) then
                Jc0mat_mag(1,pmate) = param(1)
              else if (words(1) == 'MURXX' .and. .not. kfl_muriso_mag(pmate)) then
                murmat_mag(1,pmate) = param(1)
              else if (words(1) == 'BC0XX' .and. .not. kfl_Bc0iso_mag(pmate)) then
                Bc0mat_mag(1,pmate) = param(1)
              else if (words(1) == 'TC0XX' .and. .not. kfl_Tc0iso_mag(pmate)) then
                Tc0mat_mag(1,pmate) = param(1)
              !
              ! YY components
              !
              else if (words(1) == 'RESYY' .and. .not. kfl_Resiso_mag(pmate)) then
                !
                ! Resistivity Library
                ! Opt 1 - Constant [CONST] - Parameters: RHO
                ! Opt 2 - Power Law [POWER] - Parameters: ECR, NCR, JCR
                !
                if (words(2) == 'CONST') then
                  resmat_mag(2,pmate) = 0_ip
                else if (words(2) == 'POWER') then
                  resmat_mag(2,pmate) = 1_ip
                else
                  call runend('mag_reaphy: unknown resistivity function (Available: CONST, POWER)')
                end if
              else if (words(1) == 'JCRYY') then
                !
                ! Jc Scaling Law Library
                ! Opt 1 - Constant [CONST] - Parameters: JCR
                ! Opt 2 - Kim's Model [KIMMO] - Parameters: JC0, BC0
                !
                if (words(2) == 'CONST') then
                  jcrmat_mag(2,pmate) = 0_ip
                else if (words(2) == 'KIMMO') then
                  jcrmat_mag(2,pmate) = 1_ip
                else if (words(2) == 'USERD') then
                  jcrmat_mag(2,pmate) = 9_ip
                else
                  call runend('mag_reaphy: unknown scaling law (Available: CONST, KIMMO)')
                end if
              else if (words(1) == 'RHOYY' .and. .not. kfl_rhoiso_mag(pmate)) then
                rhomat_mag(2,pmate) = param(1)
              else if (words(1) == 'ECRYY' .and. .not. kfl_Ecriso_mag(pmate)) then
                Ecrmat_mag(2,pmate) = param(1)
              else if (words(1) == 'NCRYY' .and. .not. kfl_ncriso_mag(pmate)) then
                ncrmat_mag(2,pmate) = param(1)
              else if (words(1) == 'JC0YY' .and. .not. kfl_Jc0iso_mag(pmate)) then
                Jc0mat_mag(2,pmate) = param(1)
              else if (words(1) == 'MURYY' .and. .not. kfl_muriso_mag(pmate)) then
                murmat_mag(2,pmate) = param(1)
              else if (words(1) == 'BC0YY' .and. .not. kfl_Bc0iso_mag(pmate)) then
                Bc0mat_mag(2,pmate) = param(1)
              else if (words(1) == 'TC0YY' .and. .not. kfl_Tc0iso_mag(pmate)) then
                Tc0mat_mag(2,pmate) = param(1)
              !##################################################################

              else if (words(1) == 'MOMOX') then
                momori_mag(1,pmate) = param(1)
              else if (words(1) == 'MOMOY') then
                momori_mag(2,pmate) = param(1)
              else if (words(1) == 'MOMOZ') then
                momori_mag(3,pmate) = param(1)
              end if

              call ecoute('mag_reaphy')

            end do

          end if

          call ecoute('mag_reaphy')

        end do !'ENDPR'

      else if (words(1) == 'SELFF') then

        do while (words(1) /= 'ENDSE')

          if (words(1) == 'ENABL') then

            kfl_self_mag = option('ENABL')

          else if (words(1) == 'PREVI') then

            kfl_prev_mag = option('PREVI')

          else if (words(1) == 'MATER') then

            pmate = int(param(1),ip)

            selfList_mag(pmate) = .true.

          end if

          call ecoute('mag_reaphy')

        end do

!      else if (words(1) == 'SELFF') then
!
!        kfl_self_mag = option('SELFF')
!
!      else if (words(1) == 'SELFP') then
!
!        kfl_prev_mag = option('SELFP')
!
      else if (words(1) == 'EDGEE') then

        if (option('EDGEE')) then
          kfl_edge_mag = 1_ip
        else
          kfl_edge_mag = 0_ip
          call runend('mag_reaphy: edge elements must be used')
        end if

      else if (words(1) == 'AXSYM') then

        if (option('AXSYM')) then
          kfl_axsym_mag = .true.
        else
          kfl_axsym_mag = .false.
        end if

      else if (words(1) == 'INTP1') then

        nintp1_mag = int(param(1),ip)

      else if (words(1) == 'CONST') then

        do while (words(1) /= 'ENDCO')

          if (words(1) == 'GROUP') then

            if (int(param(1),ip) == group + 1_ip) then
              group = int(param(1),ip)
            else if (int(param(1),ip) <= group) then
              call runend('mag_reaphy: repeated groups')
            else
              print *, "mag_reaphy: enter group ", group + 1_ip, " after group ", group
              call runend('mag_reaphy: unsorted list of groups')
            end if

            group_flg = .false.

            do while (words(1) /= 'ENDGR')
              if (words(1) == 'MATER') then
                pmate = int(param(1),ip)
                if (constrlist_mag(pmate) < 0) then
                  constrlist_mag(pmate) = group
                  group_flg = .true.
                  kfl_lagr_mag = .true.
                else
                  call runend('mag_reaphy: material already assigned to existing group')
                end if
              end if
              call ecoute('mag_reaphy') 
            end do

            if (group_flg) constr_total = constr_total + 1_ip

          end if

          call ecoute('mag_reaphy')

        end do

      end if

      call ecoute('mag_reaphy')

    end do

  end if

end subroutine mag_reaphy
