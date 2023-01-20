!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup mpio
!> @{
!> @file    def_mpio_file.f90
!> @author  Damien Dosimont
!> @date    08/10/2018
!> @brief   MPI-IO file data
!> @details
!> @}
!------------------------------------------------------------------------

module mod_mpio_files

  use def_kintyp,                     only : ip,rp,lg,r1p
  use def_master
  use def_mpio

  implicit none

  public

    character(150)    ::                fil_coord,&
                                        fil_ltype,&
                                        fil_lnods,&
                                        fil_lnodb,&
                                        fil_ltypb,&
                                        fil_lelbo,&
                                        fil_lmate,&
                                        fil_lmast,&
                                        fil_lesub,&
                                        fil_lbinv,&
                                        fil_leinv,&
                                        fil_lninv,&
                                        fil_lelch,&
                                        fil_lboch,&
                                        fil_lnoch,&
                                        fil_codbo,&
                                        fil_codno,&
                                        fil_leset,&
                                        fil_lbset,&
                                        fil_lnset,&
                                        fil_field

  public              ::                find_mpio_mesh_name, find_mpio_restart_mesh_name, FILE_EXIST

  contains

    subroutine find_mpio_mesh_name()
      fil_coord = adjustl(trim(namda))//trim(coord_ext)//trim(mpio_ext)
      fil_ltype = adjustl(trim(namda))//trim(ltype_ext)//trim(mpio_ext)
      fil_lnods = adjustl(trim(namda))//trim(lnods_ext)//trim(mpio_ext)
      fil_lnodb = adjustl(trim(namda))//trim(lnodb_ext)//trim(mpio_ext)
      fil_ltypb = adjustl(trim(namda))//trim(ltypb_ext)//trim(mpio_ext)
      fil_lelbo = adjustl(trim(namda))//trim(lelbo_ext)//trim(mpio_ext)
      fil_lbset = adjustl(trim(namda))//trim(lbset_ext)//trim(mpio_ext)
      fil_leset = adjustl(trim(namda))//trim(leset_ext)//trim(mpio_ext)
      fil_lnset = adjustl(trim(namda))//trim(lnset_ext)//trim(mpio_ext)
      fil_lmate = adjustl(trim(namda))//trim(lmate_ext)//trim(mpio_ext)
      fil_lmast = adjustl(trim(namda))//trim(lmast_ext)//trim(mpio_ext)
      fil_lesub = adjustl(trim(namda))//trim(lesub_ext)//trim(mpio_ext)
      fil_lelch = adjustl(trim(namda))//trim(lelch_ext)//trim(mpio_ext)
      fil_lboch = adjustl(trim(namda))//trim(lboch_ext)//trim(mpio_ext)
      fil_lnoch = adjustl(trim(namda))//trim(lnoch_ext)//trim(mpio_ext)
      fil_codno = adjustl(trim(namda))//trim(codno_ext)//trim(mpio_ext)
      fil_codbo = adjustl(trim(namda))//trim(codbo_ext)//trim(mpio_ext)
    end subroutine

    subroutine find_mpio_restart_mesh_name()
      fil_coord = adjustl(trim(namda))//trim(coord_ext)//trim(post_ext)//trim(mpio_ext)
      fil_ltype = adjustl(trim(namda))//trim(ltype_ext)//trim(post_ext)//trim(mpio_ext)
      fil_lnods = adjustl(trim(namda))//trim(lnods_ext)//trim(post_ext)//trim(mpio_ext)
      fil_lnodb = adjustl(trim(namda))//trim(lnodb_ext)//trim(post_ext)//trim(mpio_ext)
      fil_ltypb = adjustl(trim(namda))//trim(ltypb_ext)//trim(post_ext)//trim(mpio_ext)
      fil_lelbo = adjustl(trim(namda))//trim(lelbo_ext)//trim(post_ext)//trim(mpio_ext)
      fil_lbset = adjustl(trim(namda))//trim(lbset_ext)//trim(post_ext)//trim(mpio_ext)
      fil_leset = adjustl(trim(namda))//trim(leset_ext)//trim(post_ext)//trim(mpio_ext)
      fil_lnset = adjustl(trim(namda))//trim(lnset_ext)//trim(post_ext)//trim(mpio_ext)
      fil_lmate = adjustl(trim(namda))//trim(lmate_ext)//trim(post_ext)//trim(mpio_ext)
      fil_lmast = adjustl(trim(namda))//trim(lmate_ext)//trim(post_ext)//trim(mpio_ext)
      fil_lesub = adjustl(trim(namda))//trim(lesub_ext)//trim(post_ext)//trim(mpio_ext)
      fil_lelch = adjustl(trim(namda))//trim(lelch_ext)//trim(post_ext)//trim(mpio_ext)
      fil_lboch = adjustl(trim(namda))//trim(lboch_ext)//trim(post_ext)//trim(mpio_ext)
      fil_lnoch = adjustl(trim(namda))//trim(lnoch_ext)//trim(post_ext)//trim(mpio_ext)
      fil_codno = adjustl(trim(namda))//trim(codno_ext)//trim(post_ext)//trim(mpio_ext)
      fil_codbo = adjustl(trim(namda))//trim(codbo_ext)//trim(post_ext)//trim(mpio_ext)
      fil_leinv = adjustl(trim(namda))//trim(leinv_ext)//trim(post_ext)//trim(mpio_ext)
      fil_lninv = adjustl(trim(namda))//trim(lninv_ext)//trim(post_ext)//trim(mpio_ext)
      fil_lbinv = adjustl(trim(namda))//trim(lbinv_ext)//trim(post_ext)//trim(mpio_ext)
    end subroutine

    function xfiel_name(ifiel, istep) result(fil_field)
      integer(ip), intent(in)         ::  ifiel, istep
      character(300)                    ::  fil_field
      fil_field = adjustl(trim(namda))//trim(field_ext)//"."//intto8(ifiel)//"."//intto8(istep)//trim(mpio_ext)
    end function

    function FILE_EXIST(filename) result(exists)
      character(*)                    ::  filename
      logical                         ::  exists
      exists=.false.
      INQUIRE(FILE=filename, EXIST=exists)
    end function

    function intto8(v) result(v_out)
      integer(ip),  intent(in) :: v
      character(8)             :: v_out
      write(v_out,'(a)') 'XXXXXXXX'
      if(      v < 10      ) then
         write(v_out,'(a,i1)') '0000000',v
      else if( v < 100      ) then
         write(v_out,'(a,i2)') '000000',v
      else if( v < 1000     ) then
         write(v_out,'(a,i3)') '00000',v
      else if( v < 10000    ) then
         write(v_out,'(a,i4)') '0000',v
      else if( v < 100000   ) then
         write(v_out,'(a,i5)') '000',v
      else if( v < 1000000  ) then
         write(v_out,'(a,i6)') '00',v
      else if( v < 10000000 ) then
         write(v_out,'(a,i7)') '0',v
      else if( v < 100000000 ) then
         write(v_out,'(i8)') v
      end if
    end function

end module
