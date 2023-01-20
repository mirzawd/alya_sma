!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine pts_openfi(itask)
  !------------------------------------------------------------------------
  !****f* partis/pts_openfi
  ! NAME 
  !    pts_openfi  
  ! DESCRIPTION
  !    itask = 1 -- open all files
  !    itask = 2 -- close previously open pts.res and open pts.000<timestep>.res 
  ! USES
  ! USED BY
  !    pts_turnon
  !***
  !------------------------------------------------------------------------
  use def_partis
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use mod_iofile
  use mod_outfor,    only : outfor
  use mod_opfpos,    only : postpr_intto8
  use mod_result_io, only : pts_result_io

  implicit none
  integer(ip), intent(in) :: itask 
  character(7)            :: statu
  character(11)           :: forma
  character(6)            :: posit

  if( INOTSLAVE .or. itask == 10 ) then

     if( kfl_rstar == 2 ) then
        statu = 'old'
        forma = 'formatted'
        posit = 'append'
     else
        statu = 'unknown'
        forma = 'formatted'
        posit = 'asis'
     end if

     select case ( itask )

     case (   1_ip )
        !
        ! Open result file
        !
        if ( kfl_naked == 0 ) then
           call GET_ENVIRONMENT_VARIABLE('FOR1725',fil_resul_pts)
           call GET_ENVIRONMENT_VARIABLE('FOR1728',fil_oudep_pts)
           call GET_ENVIRONMENT_VARIABLE('FOR1729',fil_depsu_pts)
        else if ( kfl_naked == 1 ) then
           fil_resul_pts      = adjustl(trim(namda))//'.'           //exmod(modul)//'.res'
           fil_oudep_pts      = adjustl(trim(namda))//'-deposition.'//exmod(modul)//'.csv'
           fil_depsu_pts      = adjustl(trim(namda))//'-surface_depo.'//exmod(modul)//'.res'
   !!!        fil_num_partis_pts = adjustl(trim(namda))//'-particles.'//exmod(modul)//'.res' 
           !           
           !c etait mieux avant
           !
           !fil_oudep_pts = adjustl(trim(namda))//'-deposition.'//exmod(modul)//'.res' 
        end if


        ! Results file
        call pts_result_io % open_file()

        !
        ! Convergence file
        !
        if( kfl_thermo_pts /= 0 ) then
           call outfor(55_ip,momod(modul) % lun_conve,' ')
        else
           call outfor(54_ip,momod(modul) % lun_conve,' ')
        endif

        !
        ! Deposition file
        !
        if( kfl_oudep_pts /= 0 ) then
           if( kfl_rstar == 2 ) then 
              kfl_reawr = 1
              call iofile(4_ip,lun_oudep_pts,fil_oudep_pts,'LAGRANGIAN PARTICLES POSITION','old','unformatted')
              if( kfl_reawr == 1 ) then
                 call iofile(zero,lun_oudep_pts,fil_oudep_pts,'LAGRANGIAN PARTICLES POSITION','old','formatted','append')
              else
                 call iofile(zero,lun_oudep_pts,fil_oudep_pts,'LAGRANGIAN PARTICLES POSITION')
                 igene = kfl_oudep_pts
                 call outfor(57_ip,lun_oudep_pts,' ')
              end if
              kfl_reawr = 0
           else
              call iofile(zero,lun_oudep_pts,fil_oudep_pts,'LAGRANGIAN PARTICLES POSITION')
              !igene = kfl_oudep_pts
              !call outfor(57_ip,lun_oudep_pts,' ')
              call pts_deposition_header()
              
           end if
        end if
        !
        ! Surface Deposition file
        !
        if( kfl_depos_surface_pts /= 0 ) then           
           call iofile(zero,lun_depsu_pts,fil_depsu_pts,'DEPOSITION SURFACE')
           call outfor(101_ip,lun_depsu_pts,' ')
        end if

      end select

  end if

end subroutine pts_openfi

