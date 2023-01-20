!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine exm_outinf
!-----------------------------------------------------------------------
!****f* Exmedi/exm_outinf
! NAME 
!    exm_outinf
! DESCRIPTION
!    This routine writes on the module output files
! USED BY
!    exm_turnon
!***
!-----------------------------------------------------------------------
  use      def_master

  use      def_exmedi

  implicit none
!  character(60) :: equat
!  integer(ip)   :: ierhs,imate
!
! Headings
!
  if(kfl_paral<=0) then
  end if

!
! Write information in Result file
!


!!!!!!!! acaaaaaaaaaaaaaa ver los outputs


!  equat=''
!  equat=trim(equat)//'rho*Cp*dT/dt '
!  if(kfl_advec_exm>=1) equat=trim(equat)//'+rho*Cp*u.grad(T) '
!                       equat=trim(equat)//'-div[k*grad(T)] '
!  if(react_exm>zeexm)  equat=trim(equat)//'+s*T'
!                       equat=trim(equat)//'= '
!  ierhs=0
!  if(kfl_sourc_exm==1) then
!     equat=trim(equat)//'Q '
!     ierhs=1
!  end if
!  if(ierhs==0) equat=trim(equat)//'0 '

!  write(lun_resul_exm,110) trim(equat),nmate_exm
!  write(lun_resul_exm,111) 'rho [ M / L^3 ]     = ',(densi_exm(imate),imate=1,nmate_exm)
!  write(lun_resul_exm,111) 'Cp  [ L^2 / T^3 K ] = ',(sphea_exm(imate),imate=1,nmate_exm)
!  write(lun_resul_exm,111) 'k   [ M L / T^3 K ] = ',(tcond_exm(imate),imate=1,nmate_exm)
!  write(lun_resul_exm,*) 
!  if(react_exm>zeexm) then
!     write(lun_resul_exm,111) 's   [ M / L T^3 K ] = ',react_exm
!     write(lun_resul_exm,*) 
!  end if

!
! Formats
!
!100 format(///,&
!         5x,'****************************** ',/,&
!         5x,'ACTION POTENTIALS RESULTS FOR MODEL: ',a10,/,&
!         5x,'****************************** ',/)
!101 format(///,&
!         5x,'******************************************* ',&
!         /,&
!         5x,'EVOLUTION OF THE ACTION POTENTIALS RUN FOR MODEL: ',&
!         a10,/,&
!         5x,'******************************************* ')
!102 format(///,&
!         5x,'***************************************** ',&
!         /,&
!         5x,'ACTION POTENTIALS SOLVER INFORMATION FOR MODEL: ',&
!         a10,/,&
!         5x,'***************************************** ',/)
!110 format(/,&
!       5x, '>>>  DIFFERENTIAL EQUATION:',///,&
!       10x,a,//,&
!       10x,'Physical properties: ',i2,' material(s)',/,&
!       10x,'------------------------------------------------------------ ')
!111 format(&
!       10x,a,10e12.6)

end subroutine exm_outinf
      
