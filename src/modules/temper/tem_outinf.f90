!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_outinf()

  !-----------------------------------------------------------------------
  !
  ! This routine writes on the temperature files
  !
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_temper
  use def_kermod, only : turmu_ker
  use mod_outfor, only : outfor
  implicit none
  character(60) :: equat
  integer(ip)   :: ierhs
  character(2)  :: wcpcv

  if(kfl_paral<=0) then
     !
     ! Write information in Result file
     !
     if(kfl_rstar/=2) then
        equat=''
        if(kfl_regim_tem==1) then
           wcpcv='Cv'
        else
           wcpcv='Cp'
        end if 
        if(kfl_timei_tem==1) equat=trim(equat)//'rho*'//wcpcv//'*dT/dt '
        if(kfl_advec_tem>=1) equat=trim(equat)//' + rho*'//wcpcv//'*u.grad(T) '
        if(kfl_condu_tem/=0) then
           if( turmu_ker % kfl_exist == 0_ip ) then 
              equat=trim(equat)//' - div[k*grad(T)] ' 
           else
              equat=trim(equat)//' - div[(k+kt)*grad(T)] '
           end if
           if(kfl_naxis==1) equat=trim(equat)//' - k/r*dT/dr '
        end if
        if(react_tem>zetem)  equat=trim(equat)//' + s*T'
        if(kfl_regim_tem==1) equat=trim(equat)//' + [rho*R*div(u)]*T'
        equat=trim(equat)//'= '
        ierhs=0
        if(kfl_sourc_tem==1) then
           equat=trim(equat)//' Q '
           ierhs=1
        end if
        if(kfl_exacs_tem/=0) then
           equat=trim(equat)//' Q_exact '
           ierhs=1           
        end if
        if(kfl_regim_tem>=3) then
           equat=trim(equat)//' + dp0/dt '
           ierhs=1                      
        end if
        if (kfl_nudgi_tem==1) then
           equat=trim(equat)//'+alpha(Tref-T) '
           ierhs=1
        end if
        if(ierhs==0) equat=trim(equat)//' 0 '
        call outfor(25_ip,momod(modul)%lun_outpu,'DIFFERENTIAL EQUATION')
     end if

  end if
  !
  ! Formats
  !
110 format(/,&
       & 10x,a,//,&
       & 10x,'Physical properties: ',i2,' material(s)',/,&
       & 10x,'------------------------------------------------------------ ')
111 format(&
       & 10x,a,10(e12.6,1x))
112 format(&
       & 10x,a,e12.6,a)

end subroutine tem_outinf

