!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_openfi(itask)
  !------------------------------------------------------------------------
  !****f* Alefor/tem_openfi
  ! NAME 
  !    tem_openfi
  ! DESCRIPTION
  !    Openf iles
  ! USES
  ! USED BY
  !    tem_turnon
  !***
  !------------------------------------------------------------------------
  use def_alefor
  use def_parame
  use def_master
  use def_domain
  use mod_iofile
  implicit none
  integer(ip), intent(in) :: itask 
  character(150)          :: fil_outpu_ale,fil_resta_ale
  character(7)            :: statu
  character(11)           :: forma
  character(6)            :: posit

  if( INOTSLAVE ) then

     if( kfl_rstar == 2 ) then
        statu='old'
        forma='formatted'
        posit='append'
     else
        statu='unknown'
        forma='formatted'
        posit='asis'
     end if

     select case (itask)

     case ( 0_ip )

        !----------------------------------------------------------------
        !
        ! Open Rigid Body result file          
        !        
        !----------------------------------------------------------------

        if( nrbod > 0 .and. kfl_rigid_ale == 1 ) then
           if (kfl_naked==0) then
              call GET_ENVIRONMENT_VARIABLE('FOR714',fil_outpu_ale)               
           else if (kfl_naked==1) then
              fil_outpu_ale = adjustl(trim(namda))//'-RB.res'
           end if
           call iofile(zero,lun_outpu_ale,fil_outpu_ale,'RB RESULTS',statu,forma,posit)
!
           if (kfl_naked==0) then
              call GET_ENVIRONMENT_VARIABLE('FOR715',fil_resta_ale)               
           else if (kfl_naked==1) then
              fil_resta_ale = adjustl(trim(namda))//'-RB.rst'
           end if
           call iofile(zero,lun_resta_ale,fil_resta_ale,'ALE RESTART','unknown','unformatted')   
!           ojo luego podria ver de meterlo en el mismo donde se escriben los DISPLM 

        end if

     case (2_ip)
        !
        ! Open files needed occasionally
        !
        if( kfl_naked == 0 ) then
        else
        end if
        
     case (3_ip)
        !
        ! Close rigid body files
        !
        if( nrbod > 0 .and. kfl_rigid_ale == 1 ) then
           if (kfl_naked==0) then
              call GET_ENVIRONMENT_VARIABLE('FOR714',fil_outpu_ale)               
           else if (kfl_naked==1) then
              fil_outpu_ale = adjustl(trim(namda))//'-RB.res'
           end if
           call iofile(2_ip,lun_outpu_ale,fil_outpu_ale,'RB RESULTS')
!
           if (kfl_naked==0) then
              call GET_ENVIRONMENT_VARIABLE('FOR715',fil_resta_ale)               
           else if (kfl_naked==1) then
              fil_resta_ale = adjustl(trim(namda))//'-RB.rst'
           end if
           call iofile(2_ip,lun_resta_ale,fil_resta_ale,'ALE RESTART')   

        end if        
     end select

  end if

end subroutine ale_openfi
