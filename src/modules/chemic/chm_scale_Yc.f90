!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_scale_Yc()
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_scale_Yc
  ! NAME
  !    chm_scale_Yc
  ! DESCRIPTION
  !    Initialization of the reaction progress variable with equilibrium values
  !
  ! USES
  ! USED BY
  !    chm_iniunk
  !
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only          : ip,rp
  use def_master, only          : conce, therm
  use def_domain, only          : npoin
  use def_chemic, only          : kfl_fields_scale_chm,ncomp_chm, &
                                  xZr_chm,xZs_chm
  use mod_interp_tab, only      : fw_scale_cont_var
  use def_chemic,     only      : table_fw
  implicit none
  integer(ip)               :: ipoin,idimt,ind_cmean

  real(rp)                  :: control(table_fw % main_table % ndim)
  real(rp)                  :: scale_control(table_fw % main_table % ndim)
  real(rp)                  :: lim_control(table_fw % main_table % ndim,2_ip)
  integer(ip)               :: ind(table_fw % main_table % ndim)


  !
  ! Initialization given by geometrical fields with c = 1 or c = 0
  !
  ind = 1_ip
  do ipoin = 1,npoin

        control = 0.0_rp
        do idimt = 1, table_fw % main_table % ndim
           if (table_fw % kfl_chm_control(idimt) > 0) then
              !
              ! >0: one of the conces
              !
              control(idimt) = conce(ipoin,table_fw % kfl_chm_control(idimt),1)
           else
              if (table_fw % kfl_chm_control(idimt) == -1) then
                 !
                 ! -1: enthalpy
                 !
                 control(idimt) = therm(ipoin,1)
              elseif (table_fw % kfl_chm_control(idimt) == -2) then
                 !
                 ! -2: scalar dissipation rate
                 !
                 control(idimt) = xZr_chm(ipoin) + xZs_chm(ipoin)
              endif
           endif
           select case (table_fw % main_table % coords(idimt) % name)
           case ('CMEAN','C    ')
               ind_cmean = idimt
           end select
        enddo

        call fw_scale_cont_var( control, scale_control, lim_control, table_fw, ind)

        select case (kfl_fields_scale_chm)
            !
            ! Impose equilibrium value: y_c_eq
            !
            case (1_ip)
                if ( conce(ipoin,1,1) > 0.0_rp ) then
                    conce(ipoin,1,1:ncomp_chm) = lim_control(ind_cmean,2)
                else
                    conce(ipoin,1,1:ncomp_chm) = 0.0_rp
                end if

            !
            ! convert c to Yc
            !
            case (2_ip)
                conce(ipoin,1,1:ncomp_chm) =lim_control(ind_cmean,1) + conce(ipoin,1,1) * (lim_control(ind_cmean,2)&
                    - lim_control(ind_cmean,1))
        end select
  end do
end subroutine chm_scale_Yc
