!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_opfpos

  use def_kintyp
  use def_master
  use def_postpr,      only : postpr_intto8
  use def_postpr,      only : fsize_pos  
  use def_postpr,      only : nunam_pos  
  use def_mpio,        only : post_ext
  use def_mpio,        only : mpio_ext
  use def_parall,      only : kfl_repart_par
  use def_parall,      only : kfl_repart_post_par
  use def_parall,      only : num_repart_par
  use def_AMR,         only : kfl_amr
  use def_AMR,         only : kfl_amr_post
  use def_AMR,         only : num_amr
  use def_postpr,      only : postpr_mesh_variable
  use mod_mpio_config, only : mpio_config
  use mod_run_config,  only : run_config
  implicit none

  public :: opfpos_name, opfposvx

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-01-09
  !> @brief   Return the name of the postprocess file
  !> @details Return the name of the postprocess file
  !> 
  !-----------------------------------------------------------------------

  subroutine opfpos_name(name,extension,mesh,TAG1,TAG2,wfile)
    
    character(5),             intent(in)             :: name
    character(*),             intent(in)             :: extension
    logical(lg),              intent(out)            :: mesh
    integer(ip),              intent(in),  optional  :: TAG1
    integer(ip),              intent(in),  optional  :: TAG2
    character(len=fsize_pos), intent(out), optional  :: wfile
    character(len=fsize_pos)                         :: wtags
    character(len=fsize_pos)                         :: filsa
    character(len=fsize_pos)                         :: filen

    mesh=.false.

    !-------------------------------------------------------------------
    !
    ! Open postprocess file
    !
    !-------------------------------------------------------------------
    !
    ! Time step identifier
    !
    nunam_pos=postpr_intto8(ittim)
    !
    ! Compose file name: example.VELOC-123456.h5
    !
    filsa = adjustl(trim(namda))
    !
    ! If repartitioning is on:
    ! KFL_REPART_POST_PAR = 0 ... example-repartition#.VELOC-123456.h5
    !                     = 1 ... repartition#/example.VELOC-123456.h5
    !
    if(      kfl_repart_par /= 0 .and. kfl_repart_post_par == 1 ) then
       filsa = trim(filsa)//'-repartition'//trim(intost(num_repart_par))
    else if( kfl_repart_par /= 0 .and. kfl_repart_post_par == 2 ) then
       filsa = 'repartition'//trim(intost(num_repart_par))//'/'//trim(filsa)
    end if
    !
    ! If AMR is on:
    ! KFL_AMR_POST = 0 ... example-amr#.VELOC-123456.h5
    !              = 1 ... amr#/example.VELOC-123456.h5
    !
    if(      kfl_amr /= 0 .and. kfl_amr_post == 1 ) then
       filsa = trim(filsa)//'-amr'//trim(intost(num_amr))
    else if( kfl_amr /= 0 .and. kfl_amr_post == 2 ) then
       filsa = 'amr'//trim(intost(num_amr))//'/'//trim(filsa)
    end if
    !
    ! File name
    !    
    filsa = trim(filsa)//'-'//trim(name)
    wtags = ''
    if( present(TAG1) ) then
       if (TAG1 /= -1) then
          wtags = trim(wtags)//'.'//postpr_intto8(TAG1)
       end if
    end if
    if( present(TAG2) ) then
       if (TAG2 /= -1) then
          wtags = trim(wtags)//'.'//postpr_intto8(TAG2)
       end if
    end if
    filsa = trim(filsa)//trim(wtags)
    !
    ! Check if this is a geometric file
    !
    if( postpr_mesh_variable(name) ) then
       !
       ! Open geometry arrays file, change name
       !
       if( mpio_config%output%post_process%export_only ) then
          filen = trim(filsa)//trim(extension)
       else
          filen = trim(filsa)//trim(post_ext)//trim(extension)          
       end if
       mesh = .true.
       
    else

       if( kfl_reawr == 0 ) then
          !
          ! Open postprocess file name
          !
          filen = trim(filsa)//'-'//adjustl(trim(nunam_pos))//trim(post_ext)//trim(extension)

       else
          if( modul > 0 ) then
             filen = trim(filsa)//'.'//exmod(modul)//'.rst'
             !filen = trim(filsa)//'.rst'
          else
             filen = trim(filsa)//'.rst'             
          end if
          if (extension==mpio_ext) then
             filen = trim(filen)//trim(extension)
          end if
          if( kfl_reawr == 2) then
             !
             ! Open restart file name for writing
             !
             if( run_config%restart%append_timestep ) call appnum(ittim,filen)
          end if
       end if
    end if

    if( present(wfile) ) then
       wfile     = trim(filen)
    else
       fil_postp = trim(filen)
    end if

  end subroutine opfpos_name
  
  !-----------------------------------------------------------------------
  !****f* Domain/opfpos
  ! NAME
  !    opfpos
  ! DESCRIPTION
  !    This subroutine gets ALL the file names to be used by Alya in two
  !    possible ways and them open them:
  !
  !    1. Recalling them from the environment, when Alya is launched
  !    encapsulated in a shell script, or
  !
  !    2. Composing the names out of the problem name which is given as argument
  !    when the binary file Alya is launched "naked".
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------

  subroutine opfposvx()

    implicit none
    character(len=fsize_pos) :: filsa

    !-------------------------------------------------------------------
    !
    ! Open voxel file
    !
    !-------------------------------------------------------------------

    if( INOTSLAVE ) then
       nunam_pos=postpr_intto8(ittim)
       filsa = adjustl(trim(namda))
       fil_posvx = trim(filsa)//'-'//trim(wopos_pos(1))//'-'//adjustl(trim(nunam_pos))//'.bvox'
       open (unit=lun_posvx, FILE = trim(fil_posvx) ,ACTION='WRITE', &
            & STATUS = 'UNKNOWN', form='unformatted', access='stream' )
    end if

  end subroutine opfposvx

end module mod_opfpos
!> @}
