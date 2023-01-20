!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_outvar(ivari,imesh)
  !------------------------------------------------------------------------
  !****f* Levels/lev_output
  ! NAME 
  !    lev_output
  ! DESCRIPTION
  !    Output a postprocess variable
  ! USES
  !    postpr
  !    memgen
  ! USED BY
  !    lev_output
  !***
  !------------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_levels
  use      mod_postpr
  use      mod_memchk
  use mod_outvar,         only : outvar
  implicit none
  integer(ip), intent(in) :: ivari
  integer(ip), intent(in) :: imesh
  integer(ip)             :: dummi
  real(rp)                :: rutim

  rutim = cutim

  select case (ivari)  

  case(0_ip)
     !
     ! Nothing
     !
     return

  case(1_ip)
     !
     ! Level Set function
     !
     gesca => fleve(:,1) 

  case(2_ip)
     !
     ! Velocity 
     !
     if( INOTMASTER ) then
        if(kfl_advec_lev==1) then
           gevec => veloc(:,:,1)
        else
           call memgen(zero,ndime,npoin)
           call lev_velfun(kfl_advec_lev,ndime,npoin,dummi,coord,gevec)
        end if
     endif

  case(3_ip)
     !
     ! grad phi /| grad phi| 
     !
     if( INOTMASTER ) then
        gevec => norml_lev(:,:)
     endif

  case(4_ip)
     !
     ! grad phi /| grad phi| 
     !
     gesca => dista_lev

  case(5_ip)
     !
     ! Displacement
     !
     gevec => dispm(:,:,1)

  case(6_ip)
     !
     ! Mesh velocity
     !
     gevec => velom

  end select
  !
  ! Postprocess
  !
  call outvar(&
       ivari,&
       ittim,rutim,postp(1)%wopos(:,ivari),MESH_ID=imesh)

end subroutine lev_outvar
