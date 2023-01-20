!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Postprocess
!> @{
!> @file    mod_postpr_tools.f90
!> @author  houzeaux
!> @date    2020-04-05
!> @brief   Tools for postprocess
!> @details Simple tools used mo mod_postpr 
!-----------------------------------------------------------------------

module mod_postpr_tools

  use def_kintyp_basic,      only : ip,rp,lg
  use def_master,            only : pleng_mpio
  use def_master,            only : pdime_mpio   
  use def_master,            only : ttime_mpio   
  use def_master,            only : itste_mpio   
  use def_master,            only : tag1_mpio
  use def_master,            only : tag2_mpio
  use def_master,            only : kfl_reawr
  use def_master,            only : intost
  use def_master,            only : INOTSLAVE       
  use def_master,            only : IMASTER       
  use def_master,            only : lun_pos02
  use def_master,            only : npari
  use def_master,            only : nparr  
  use def_master,            only : nparc  
  use def_master,            only : nparx  
  use def_master,            only : pardi
  use def_master,            only : npart
  use def_master,            only : kfl_filte
  use def_master,            only : lun_postp
  use def_postpr,            only : fsize_pos
  use def_postpr,            only : alyabin_ext
  use def_postpr,            only : wopos_mpio
  use def_mpio,              only : mpio_ext
  use mod_opfpos,            only : opfpos_name
  use def_postpr,            only : iiiii
  use def_postpr,            only : rrrrr
  use def_postpr,            only : wwwww
  use def_postpr,            only : wwww8
  use def_postpr,            only : nwwww
  use def_postpr,            only : ncoun_pos
  use mod_iofile,            only : iofile_flush_unit
  use mod_optional_argument, only : optional_argument
  use mod_mpio_config,       only : mpio_config
  implicit none

contains

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Compose a name with component number
  !> @details Example: CONCE <= CO001
  !>
  !----------------------------------------------------------------------

  subroutine postpr_components(icomp,wopos,ncomp)

    integer(ip),                intent(in)    :: icomp
    character(len=5),           intent(inout) :: wopos
    integer(ip),      optional, intent(in)    :: ncomp
    integer(ip)                               :: ii,ncomp_loc

    ii        = len_trim(wopos)
    ncomp_loc = optional_argument(1_ip,ncomp)

    if( ncomp_loc < 100 ) then
       if(      icomp < 10   ) then
          wopos = wopos(1:3)//'0'//trim(intost(icomp))
       else
          wopos = wopos(1:3)//trim(intost(icomp))
       end if
    else
       if(      icomp < 10   ) then
          wopos = wopos(1:2)//'00'//trim(intost(icomp))
       else if( icomp < 100  ) then
          wopos = wopos(1:2)//'0'//trim(intost(icomp))
       else if( icomp < 1000 ) then
          wopos = wopos(1:2)//trim(intost(icomp))
       end if
    end if
    
  end subroutine postpr_components
   
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-04
  !> @brief   Tags for output
  !> @details Tags for restart, postprocess... also used by MPIO
  !> 
  !-----------------------------------------------------------------------

  character(50) function postpr_tags(tag1,tag2) result(wtag)

    integer(ip),    intent(in) :: tag1
    integer(ip),    intent(in) :: tag2
    character(20)              :: wtag1
    character(20)              :: wtag2
       
    if( tag1 >= 0 ) then
       wtag1 = trim(intost(tag1))
    else
       wtag1 = ''
    end if
    if( tag2 >= 0 ) then
       wtag2 = trim(intost(tag2))
    else
       wtag2 = ''
    end if

    if( tag1 >= 0 .and. tag2 >= 0 ) then
       wtag = '(COMPONENT '//trim(wtag1)//', TIME STEP '//trim(wtag2)//')'
    else if( tag1 >= 0 ) then
       wtag = '(COMPONENT '//trim(wtag1)//')'
    else if( tag2 >= 0 ) then
       wtag = '(TIME STEP '//trim(wtag2)//')'
    else
       wtag = ''
    end if

  end function postpr_tags

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-01-09
  !> @brief   Name of the postprocess file
  !> @details Give the complete name of the postprocess file just about
  !>          to be postprocessed
  !> 
  !-----------------------------------------------------------------------

  subroutine postpr_file_name(wfile,wopos,TAG1,TAG2)
    
    character(len=fsize_pos), intent(inout)          :: wfile
    character(len=5),         intent(in)             :: wopos(*)
    integer(ip),              intent(in),   optional :: TAG1
    integer(ip),              intent(in),   optional :: TAG2
    logical(lg)                                      :: mesh

   if( postpr_is_mpio() ) then
       call opfpos_name(wopos(1),mpio_ext,    mesh, TAG1, TAG2, wfile)
    else
       call opfpos_name(wopos(1),alyabin_ext, mesh, TAG1, TAG2, wfile)
    end if
 
  end subroutine postpr_file_name

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-01-09
  !> @brief   Type of postprocess
  !> @details Guess if postprocess is with mpio
  !> 
  !-----------------------------------------------------------------------
  
  function postpr_is_mpio(wopos,itste,ttime,pdime,pleng_opt,TAG1,TAG2,fti) result(ret)

#ifdef ALYA_FTI
    use mod_fti_config, only: fti_config
#endif

    implicit none

    character(5), intent(in), optional  :: wopos(*)
    integer(ip) , intent(in), optional  :: itste
    real(rp)    , intent(in), optional  :: ttime
    integer(ip),  intent(in), optional  :: pdime
    integer(ip),  intent(in), optional  :: pleng_opt
    integer(ip),  intent(in), optional  :: TAG1
    integer(ip),  intent(in), optional  :: TAG2
    logical(lg),  intent(in), optional  :: fti
    logical(lg)                         :: ret
    logical(lg)                         :: fti_enabled


#ifdef ALYA_FTI
    fti_enabled = fti_config%enabled
    if (present(fti)) then
        fti_enabled = fti_enabled .and. fti
    end if
#else
    fti_enabled = mpio_config%output%restart%enabled
#endif
    
    ret=.false.

    if ((mpio_config%output%post_process%enabled .AND. kfl_reawr == 0) .OR. ((mpio_config%output%restart%enabled .OR. fti_enabled).AND. kfl_reawr /= 0)) then

       if( present(wopos) .and. present(itste) .and. present(ttime) .and. present(pdime) ) then

          wopos_mpio(1) =  wopos(1)
          wopos_mpio(2) =  wopos(2)
          wopos_mpio(3) =  wopos(3)

          pdime_mpio    =  pdime
          ttime_mpio    =  ttime
          itste_mpio    =  itste
          if ( present(TAG1) ) then
             tag1_mpio  = TAG1
          else
             tag1_mpio = -1
          end if
          if ( present(TAG2) ) then
             tag2_mpio  = TAG2
          else
             tag2_mpio = -1
          end if
          if( present(pleng_opt) ) then
             pleng_mpio = pleng_opt
          else
             pleng_mpio = 0
          end if
          ret=.true.
       else
          ret = .true.
       end if

    end if
    
  end function postpr_is_mpio

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-02-23
  !> @brief   Header
  !> @details Write/read header of a postprocess file
  !>
  !-----------------------------------------------------------------------

  subroutine postpr_header(nunit)

    integer(ip), optional, intent(in) :: nunit
    integer(ip)                       :: ii
    integer(4)                        :: ihead,nunit4

    nunit4 = int(optional_argument(lun_postp,nunit),4)
    npari  = 0
    nparr  = 0
    nparc  = 0
    nparx  = 0
    pardi  = int(iiiii(1),ip)

    if( INOTSLAVE ) then
       if( IMASTER ) then
          wwwww(8) = 'PARAL'
          iiiii(3) = int(npart,4)
       else
          wwwww(8) = 'SEQUE'
          iiiii(3) = int(1,4)
       end if
       if( kfl_filte == 0 ) then
          wwwww(9) = 'NOFIL'
       else
          wwwww(9) = 'FILTE'
       end if

       if( kfl_reawr == 1 ) then
          read(nunit4)  ihead
          if( ihead /= 1234 ) call runend('MOD_POSTPR: ALYA FILE HAS A WRONG BINARY FORMAT')
          read(nunit4)  wwww8(1)
          read(nunit4)  wwww8(2)
          read(nunit4)  wwww8(3)
          read(nunit4)  wwww8(4)
          read(nunit4)  wwww8(5)
          read(nunit4)  wwww8(6)
          read(nunit4)  wwww8(7)
          read(nunit4)  wwww8(8)
          read(nunit4)  wwww8(9)
          read(nunit4)  iiiii(1)
          read(nunit4)  iiiii(2)
          read(nunit4)  iiiii(3)
          read(nunit4)  iiiii(4)
          if( wwww8(2) == 'V0003' ) then
             read(nunit4)  iiiii(5)
             read(nunit4)  iiiii(6)
          end if
          read(nunit4)  rrrrr(1)
          do ii = 1,nwwww
             wwwww(ii)=wwww8(ii)(1:5)
          end do
       else
          wwww8(1) = 'ALYAPOST'
          do ii = 2,nwwww
             wwww8(ii)=wwwww(ii)//'   '
          end do
          ihead = 1234_4
          write(nunit4) ihead
          write(nunit4) wwww8(1)
          write(nunit4) wwww8(2)
          write(nunit4) wwww8(3)
          write(nunit4) wwww8(4)
          write(nunit4) wwww8(5)
          write(nunit4) wwww8(6)
          write(nunit4) wwww8(7)
          write(nunit4) wwww8(8)
          write(nunit4) wwww8(9)
          write(nunit4) iiiii(1)
          write(nunit4) iiiii(2)
          write(nunit4) iiiii(3)
          write(nunit4) iiiii(4)
          write(nunit4) iiiii(5)
          write(nunit4) iiiii(6)
          write(nunit4) rrrrr(1)
          if( kfl_reawr == 0 ) then
             if( ncoun_pos == 0 ) then
                ncoun_pos = ncoun_pos + 1
                write( lun_pos02,'(a,1x,i9,1x,e12.5)') 'START',iiiii(4),rrrrr(1)
             end if
             write( lun_pos02,'(a)') wwwww(3)(1:5)
             call iofile_flush_unit(lun_pos02)
          end if

       end if
    end if

  end subroutine postpr_header

end module mod_postpr_tools
!> @}
  
