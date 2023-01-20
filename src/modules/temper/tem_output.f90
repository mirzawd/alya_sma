!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_output()
  !------------------------------------------------------------------------
  !****f* Temper/tem_output
  ! NAME 
  !    tem_output
  ! DESCRIPTION
  !    Output and postprocess of solution
  ! USES
  ! USED BY
  !    tem_iniunk
  !    tem_endite
  !    tem_endste
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_temper
  use mod_iofile
  use def_kermod
  use mod_matrix, only : matrix_output_gid_format
  use mod_ADR,    only : ADR_manufactured_error
  use mod_output_postprocess
  implicit none
  external              :: tem_outvar
  integer(ip)           :: ivarp
  !
  ! Initial solution, end of a time step and and of run
  !
  call output_postprocess_variables(tem_outvar)
  !
  ! Update reference time for time-averaging
  ! 
  do ivarp = 1,nvarp
     if( postp(1) % npp_stepi(ivarp,0) /= 0 ) then
        if( mod(ittim, postp(1) % npp_stepi(ivarp,0) ) == 0 ) then   
           avtim_tem = cutim                                    
        endif
     endif
  end do

  if( ittyp == ITASK_INITIA .or. ittyp == ITASK_ENDTIM ) then
     !
     ! View factor, only first iteration
     !
     if( kfl_radia_tem == 1 .and. ittim == 1 ) then
        call tem_radpos()
     end if
     !
     ! Boundary conditions
     !
     if( npp_bound_tem == 1 ) then 
        npp_bound_tem = 0
        call tem_outbcs()
     end if
     !
     ! Calculations on sets
     !
     call tem_outset()
     !
     ! Calculations on witness points
     !
     call tem_outwit()
     !
     ! Error w/r exact solution
     !
     if( kfl_exacs_tem /= 0 ) then
        call ADR_manufactured_error(ADR_tem,ittim,cutim,therm)
     end if

  else if( ittyp == ITASK_ENDRUN ) then
     !
     ! End of the run
     !
     
     if( kfl_splot_tem == 1 .and. ISEQUEN ) then
        call suplot(1,therm(1:npoin,1),lun_splot_tem)
        close(lun_splot_tem)
     end if
     
     if( kfl_psmat_tem /= 0 .and. INOTMASTER ) then
        if( kfl_psmat_tem > 0 ) then
           if( kfl_discr_tem == 0 ) then
              call pspltm(&
                   npoin,npoin,1_ip,0_ip,c_dom,r_dom,amatr,&
                   trim(title)//': '//namod(modul),0_ip,18.0_rp,'cm',&
                   0_ip,0_ip,2_ip,lun_psmat_tem)
           else
              call pspltm(&
                   nunkn_tem,nunkn_tem,1_ip,0_ip,meshe(ndivi) % c_elm_2,meshe(ndivi) % r_elm_2,amatr,&
                   trim(title)//': '//namod(modul),0_ip,18.0_rp,'cm',&
                   0_ip,0_ip,2_ip,lun_psmat_tem)
           end if
           call tem_openfi(4_ip)
        !else
        !    call matrix_output_gid_format(npoin,1_ip,r_dom,c_dom,amatr,momod(modul) % solve(1) % invpr_gs)
        end if
        kfl_psmat_tem = 0
     end if

  end if

end subroutine tem_output
