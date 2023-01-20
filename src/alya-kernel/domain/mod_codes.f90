!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    mod_codes.f90
!> @author  houzeaux
!> @date    2020-03-30
!> @brief   Codes
!> @details Deal with codes
!-----------------------------------------------------------------------

module mod_codes

  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use def_kermod
  use def_inpout
  use mod_iofile
  use mod_ecoute,                  only : ecoute
  use mod_outfor,                  only : outfor
  use mod_memory,                  only : memory_size
  use mod_windk,                   only : mod_windk_create_system
  use mod_ker_space_time_function, only : space_time_function_number
  use mod_ker_discrete_function,   only : ker_discrete_function_number
  use mod_tubes,                   only : tubes_register
  implicit none
  private

  integer(ip)             :: icode,idofn,iauxi
  integer(ip)             :: ndofn,ncode,nnand,iword
  integer(ip)             :: jcode,nmcod,ivcod
  integer(ip)             :: kfl_funno_tmp,kfl_funbo_tmp, kfl_funtyp_tmp, ifunc
  integer(ip)             :: kfl_fixrs_tmp
  character(5)            :: wfixrs,wfname,wtag

  public :: codes_read_on_nodes
  public :: codes_read_on_boundaries

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-30
  !> @brief   Read codes
  !> @details Read codes on nodes
  !> 
  !-----------------------------------------------------------------------

  subroutine codes_read_on_nodes(itask)
    
    integer(ip), intent(in) :: itask
    integer(ip)             :: kcode(mcono)

    call ecoute('reacod')

    if ( .not. associated(tncod) ) call runend('REACOD: NODAL CODES BUT CONDITIONS ON BOUNDARIES. MISSING EXTRAPOLATE MAYBE...')

    ncode = 0
    ndofn = tncod(1) % ndofn

    if( itask == 100 ) call outfor(84_ip,momod(modul)%lun_outpu,' ')

    do while( words(1) /= 'ENDCO' )

       nnand = 0
       if(exists('&    ')) then
          !
          ! Count number of multiple codes
          !
          nnand=0
          do iword=1,maxwp
             if(words(iword)=='&    ') nnand=nnand+1
          end do
          if(mcono==1) then
             call runend('REACOD: MULTIPLE CODES FOR NODES NOT POSSIBLE')
          end if
       end if

       if( words(2) /= 'OFF  ' ) then

          ncode = ncode + 1
          nmcod = nnand+1
          if  ( ncode > size(tncod(1) % l,1,ip) ) call runend('REACOD: NUMBER OF LINES IN BOUNDARY CONDITIONS SECTION EXCEEDS THE VALUE GIVEN  BY mcodc. MODIFY reacod.f90 AND RECOMPILE.')

          !-------------------------------------------------------------
          !
          ! Axes: local or global
          !
          !-------------------------------------------------------------

          kfl_fixrs_tmp = 0
          if( exists('AXES ') ) then
             wfixrs = getcha('AXES ','     ','#Axes')
             if( wfixrs == 'LOCAL' ) then
                !
                ! Local axes: uder defined (>0) or exnor (-1)
                !
                kfl_fixrs_tmp =  getint('LOCAL',0_ip,'#Field')
                if( kfl_fixrs_tmp > 0 ) then
                   kfl_fixrs_tmp =  1000_ip + kfl_fixrs_tmp
                else
                   kfl_fixrs_tmp = -1
                end if

             else if( wfixrs == 'FROMF' ) then
                !
                ! Field
                !
                kfl_fixrs_tmp =  getint('FROMF',1_ip,'#Field')

             else 
                kfl_fixrs_tmp =  0
             end if
          end if

          !-------------------------------------------------------------
          !
          ! Functions
          !
          !-------------------------------------------------------------

          kfl_funtyp_tmp = 0
          kfl_funno_tmp  = 0
          wfname         = '     '

          if( exists('TIMEF') ) then
             !
             ! Time function
             !
             wfname         = getcha('TIMEF','     ','#Time Function name')
             kfl_funtyp_tmp = FUNCTION_TIME

             do ifunc = 1,number_time_function
                if( trim(wfname) == trim(time_function(ifunc) % name) ) then
                   kfl_funno_tmp = ifunc
                end if
             end do
             if( kfl_funno_tmp == 0 ) call runend('REACOD: TIME FUNCTION '//trim(wfname)//' IS NOT DEFINED')

          else if( exists('FIELD') ) then
             !
             ! Field
             !
             kfl_funno_tmp  = getint('FIELD',1_ip,'#Field')   ! a default value makes no sense here ! Change sign
             kfl_funtyp_tmp = FUNCTION_FIELD

          else if( exists('MODUL') ) then
             !
             ! Module function
             !
             kfl_funno_tmp  = getint('MODUL',1_ip,'#MModule')   ! a default value makes no sense here ! Change sign
             kfl_funtyp_tmp = FUNCTION_MODULE
             
          else if( exists('SPACE') ) then
             !
             ! Space time function
             !
             wfname         = getcha('SPACE','     ','#Space/time Function name')
             kfl_funno_tmp  = space_time_function_number(wfname)
             kfl_funtyp_tmp = FUNCTION_SPACE_TIME

          else if( exists('DISCR') ) then
             !
             ! Discrete function
             !
             wfname         = getcha('DISCR','     ','#Discrete Function name')
             kfl_funno_tmp  = ker_discrete_function_number(wfname)
             kfl_funtyp_tmp = FUNCTION_DISCRETE

          else if( exists('FUNCT') ) then

             call runend('REACOD: DEPRECATED KEYWORD')

          else
             !
             ! Used for retrocompatibility
             !
             if( nnand == 0 ) then
                if( param(1) < 0 ) then
                   kfl_funno_tmp = int( param(4) )
                else
                   kfl_funno_tmp = int( param(3+ndofn) )
                end if
             else
                if( param(1) < 0 ) then
                   kfl_funno_tmp = int( param(2+nmcod+1) )
                else
                   kfl_funno_tmp = int( param(2+nmcod+ndofn) )
                end if
             end if

          end if
          
          if( kfl_funno_tmp == 0 ) kfl_funtyp_tmp = 0
          
          !-------------------------------------------------------------
          !
          ! Code tag and values
          !
          !-------------------------------------------------------------

          if( exists('TAG  ') ) then
             wtag = getcha('TAG  ','     ','#Function name')
          else
             wtag = ''
          end if
          !
          ! Values on nodes function
          !
          if( exists('VALUE') ) then
             param(1) = param(1) ! -param(1) CHANGE SIGN
             param(3) = getrea('VALUE',1.0_rp,'#VALUES ON NODE FUNCTION')
             call runend('REACOD: VALUE MUST BE REPLACED BY FIELD')
          end if

          tncod(1) % l(ncode) % tag        = wtag
          tncod(1) % l(ncode) % fname      = wfname   ! Function name

          if( itask == 100 ) ioutp(1:3) = -(mcodb+1)

          !-------------------------------------------------------------
          !
          ! Impose values on code type
          !
          !-------------------------------------------------------------

          if( nnand == 0 ) then
             !
             ! Single code per node
             !
             if( param(1) < 0.0_rp ) then

                tncod(1) % l(ncode) % lcode(1)   = abs(int( param(1) ))
                tncod(1) % l(ncode) % kfl_fixno  = int( param(2) )
                tncod(1) % l(ncode) % kfl_value  = int( param(3) )
                tncod(1) % l(ncode) % kfl_funno  = kfl_funno_tmp
                tncod(1) % l(ncode) % kfl_funtyp = kfl_funtyp_tmp
                tncod(1) % l(ncode) % kfl_fixrs  = kfl_fixrs_tmp

             else

                if( exists('FIELD') ) then
                   tncod(1) % l(ncode) % lcode(1)   = int( param(1) )
                   tncod(1) % l(ncode) % kfl_fixno  = int( param(2) )
                   tncod(1) % l(ncode) % kfl_value  = getint('FIELD',1_ip,'#FIELD TO USE AS A BOUNDARY CONDITION')
                   tncod(1) % l(ncode) % kfl_funno  = kfl_funno_tmp
                   tncod(1) % l(ncode) % kfl_funtyp = kfl_funtyp_tmp
                   tncod(1) % l(ncode) % kfl_fixrs  = kfl_fixrs_tmp
                else
                   tncod(1) % l(ncode) % lcode(1)  = int( param(1) )
                   iauxi = int( param(2) )
                   tncod(1) % l(ncode) % kfl_fixno = int( param(2) )
                   tncod(1) % l(ncode) % kfl_value = 0
                   do idofn = 1,ndofn
                      tncod(1) % l(ncode) % bvess(idofn) = param(2+idofn)
                   end do
                   tncod(1) % l(ncode) % kfl_funno  = kfl_funno_tmp
                   tncod(1) % l(ncode) % kfl_funtyp = kfl_funtyp_tmp
                   tncod(1) % l(ncode) % kfl_fixrs  = kfl_fixrs_tmp
                end if
             end if

             if( itask == 100 ) then
                ioutp(1) = tncod(1) % l(ncode) % lcode(1)
             end if

          else
             !
             ! Multiple code per node
             !
             do jcode = 1,nmcod
                kcode(jcode) = abs(int(param(jcode)))
             end do
             call heapsorti1(2_ip,nmcod,kcode)

             if( param(1) < 0.0_rp ) then
                do jcode = 1,nmcod
                   tncod(1) % l(ncode) % lcode(jcode)  = kcode(jcode)
                end do
                tncod(1) % l(ncode) % kfl_fixno  = int( param(nmcod+1)   )
                tncod(1) % l(ncode) % kfl_value  = int( param(nmcod+2)   )
                tncod(1) % l(ncode) % kfl_funno  = kfl_funno_tmp
                tncod(1) % l(ncode) % kfl_funtyp = kfl_funtyp_tmp
                tncod(1) % l(ncode) % kfl_fixrs  = kfl_fixrs_tmp
             else
                if( exists('FIELD') ) then
                   do jcode = 1,nmcod
                      tncod(1) % l(ncode) % lcode(jcode)  = kcode(jcode)
                   end do
                   tncod(1) % l(ncode) % kfl_fixno  = int( param(nmcod+1)   )
                   tncod(1) % l(ncode) % kfl_value  = getint('FIELD',1_ip,'#FIELD TO USE AS A BOUNDARY CONDITION')
                   tncod(1) % l(ncode) % kfl_funno  = kfl_funno_tmp
                   tncod(1) % l(ncode) % kfl_funtyp = kfl_funtyp_tmp
                   tncod(1) % l(ncode) % kfl_fixrs  = kfl_fixrs_tmp
                else
                   do jcode = 1,nmcod
                      tncod(1) % l(ncode) % lcode(jcode)  = kcode(jcode)
                   end do
                   tncod(1) % l(ncode) % kfl_fixno = int( param(nmcod+1) )
                   tncod(1) % l(ncode) % kfl_value = 0
                   do idofn = 1,ndofn
                      tncod(1) % l(ncode) % bvess(idofn) = param(1+nmcod+idofn)
                   end do
                   tncod(1) % l(ncode) % kfl_funno  = kfl_funno_tmp
                   tncod(1) % l(ncode) % kfl_funtyp = kfl_funtyp_tmp
                   tncod(1) % l(ncode) % kfl_fixrs  = kfl_fixrs_tmp
                end if
             end if

             if( itask == 100 ) then
                do jcode = 1,nmcod
                   ioutp(jcode) = tncod(1) % l(ncode) % lcode(jcode)
                end do

             end if

          end if
          if( itask == 100 ) call outfor(43_ip,momod(modul)%lun_outpu,' ')
          !
          ! Check errors
          !
          if( param(1) < 0 .and. ( .not. READ_AND_RUN() ) ) then
             ivcod = tncod(1) % l(ncode) % kfl_value
             if( ivcod > mfiel .or. ivcod < 0 ) then
                call outfor(48_ip,0_ip,'WRONG VALUE CODE')
             end if
             if( kfl_field(1,ivcod) == 0 ) then
                call outfor(48_ip,0_ip,'VALUE CODE WAS NOT DEFINED')
             end if
          end if

       end if

       call ecoute('reacod')
    end do

    tncod(1) % ncode = ncode

  end subroutine codes_read_on_nodes

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-30
  !> @brief   Read codes
  !> @details Read codes on nodes
  !> 
  !-----------------------------------------------------------------------

  subroutine codes_read_on_boundaries(itask)

    integer(ip), intent(in) :: itask
    
    ncode = 0
    ndofn = tbcod(1) % ndofn

    call ecoute('reacod')

    if( itask == 200 ) call outfor(85_ip,momod(modul)%lun_outpu,' ')

    do while( words(1) /= 'ENDCO' )

       ncode = ncode + 1
       
       !-------------------------------------------------------------
       !
       ! Functions
       !
       !-------------------------------------------------------------

       kfl_funbo_tmp  = 0
       kfl_funtyp_tmp = 0
       wfname = '     '
       
       if( exists('TIMEF') ) then
          !
          ! Time function
          !
          wfname         = getcha('TIMEF','     ','#Time Function name')
          kfl_funtyp_tmp = FUNCTION_TIME
          do ifunc = 1,number_time_function
             if( trim(wfname) == trim(time_function(ifunc) % name) ) then
                kfl_funbo_tmp = ifunc
             end if
          end do
          if( kfl_funbo_tmp == 0 ) call runend('REACOD: TIME FUNCTION '//trim(wfname)//' IS NOT DEFINED')
          
       else if( exists('SPACE') ) then
          !
          ! Space time function
          !
          wfname         = getcha('SPACE','     ','#Space/time Function name')
          kfl_funtyp_tmp = FUNCTION_SPACE_TIME
          kfl_funbo_tmp  = space_time_function_number(wfname)              

       else if( exists('MODUL') ) then
          !
          ! Module function
          !
          kfl_funbo_tmp  = getint('MODUL',1_ip,'#MModule')   ! a default value makes no sense here ! Change sign
          kfl_funtyp_tmp = FUNCTION_MODULE
             
       else if( exists('DISCR') ) then
          !
          ! Discrete function
          !
          wfname         = getcha('DISCR','     ','#Discrete Function name')
          kfl_funbo_tmp  = ker_discrete_function_number(wfname)
          kfl_funtyp_tmp = FUNCTION_DISCRETE
          
       else if( exists('WINDK') ) then
          !
          ! Windkessel function
          !
          wfname         = getcha('WINDK','     ','#Windkessel function name')
          kfl_funtyp_tmp = FUNCTION_WINDKESSEL
          do ifunc = 1,number_windk_systems
             if( trim(wfname) == trim(windk_systems(ifunc) % name) ) then
                kfl_funbo_tmp = ifunc
             end if
          end do        
          if( kfl_funbo_tmp == 0 ) call runend('REACOD: WINDKESSEL FUNCTION '//trim(wfname)//' IS NOT DEFINED')
          
          !! Create windkessel system
          call mod_windk_create_system(wfname,modul,int(param(1),ip), OPT_iflow_nsi=int(param(6),ip))
         
       else if( exists('PUMPF') ) then
          !
          ! Pump function
          !
          wfname         = getcha('PUMPF','     ','#Pump function name')
          kfl_funtyp_tmp = FUNCTION_PUMP
          do ifunc = 1,number_pump_curve
             if( trim(wfname) == trim(pump_curve(ifunc) % name) ) then
                kfl_funbo_tmp = ifunc
             end if
          end do
          if( kfl_funbo_tmp == 0 ) call runend('REACOD: PUMP FUNCTION '//trim(wfname)//' IS NOT DEFINED') 
          
       else if( exists('TUBES') ) then
          !
          ! Tubes function
          !
          icode          = int( param(1) )
          kfl_funbo_tmp  = getint('TUBES',1_ip,'#Tubes segment connection')
          kfl_funtyp_tmp = FUNCTION_TUBES
          if( kfl_funbo_tmp == 0 ) call runend('REACOD: TUBES FUNCTION '//trim(wfname)//' IS NOT DEFINED')

          ! Register code correspondance
          call tubes_register(kfl_funbo_tmp,icode,modul)
          
       else if( exists('FUNCT') ) then
          !
          ! Wrong keyword
          !
          call runend('REACOD: DEPRECATED OPTION')
          
       else
          
          if( param(1) < 0 ) then
             kfl_funbo_tmp = int( param(4) )
          else
             kfl_funbo_tmp = int( param(3+ndofn) )
          end if
          
       end if
 
       !-------------------------------------------------------------
       !
       ! Tags and values
       !
       !-------------------------------------------------------------

       if( exists('TAG  ') ) then
          wtag = getcha('TAG  ','     ','#Function name')
       else
          wtag = ''
       end if
       if( exists('VALUE') ) then
          param(1) = param(1) ! - param(1) CHANGE SIGN
          param(3) = getrea('VALUE',1.0_rp,'#VALUES ON NODE FUNCTION')
       end if

       !-------------------------------------------------------------
       !
       ! Impose values on code type
       !
       !-------------------------------------------------------------

       if( param(1) < 0.0_rp ) then

          tbcod(1) % l(ncode) % lcode      = int( param(1) )
          tbcod(1) % l(ncode) % kfl_fixbo  = int( param(2) )
          tbcod(1) % l(ncode) % kfl_value  = int( param(3) )
          tbcod(1) % l(ncode) % kfl_funbo  = kfl_funbo_tmp
          tbcod(1) % l(ncode) % kfl_funtyp = kfl_funtyp_tmp
          tbcod(1) % l(ncode) % tag        = wtag

       else

          if( exists('FIELD') ) then
             tbcod(1) % l(ncode) % lcode      = int( param(1) )
             tbcod(1) % l(ncode) % kfl_fixbo  = int( param(2) )
             tbcod(1) % l(ncode) % kfl_value  = getint('FIELD',1_ip,'#FIELD TO USE AS A BOUNDARY CONDITION')
             tbcod(1) % l(ncode) % kfl_funbo  = kfl_funbo_tmp
             tbcod(1) % l(ncode) % kfl_funtyp = kfl_funtyp_tmp
             tbcod(1) % l(ncode) % tag        = wtag
          else
             tbcod(1) % l(ncode) % lcode     = int( param(1) )
             tbcod(1) % l(ncode) % kfl_fixbo = int( param(2) )
             do idofn = 1,ndofn
                tbcod(1) % l(ncode) % bvnat(idofn) = param(2+idofn)
             end do
             tbcod(1) % l(ncode) % kfl_funbo  = kfl_funbo_tmp
             tbcod(1) % l(ncode) % kfl_funtyp = kfl_funtyp_tmp
             tbcod(1) % l(ncode) % tag       = wtag
          end if

       end if

       ioutp(1:3) = -(mcodb+1)
       if( itask == 200 ) then
          ioutp(1) =  tbcod(1) % l(ncode) % lcode
       end if

       call ecoute('reacod')

    end do

    tbcod(1) % ncode = ncode

  end subroutine codes_read_on_boundaries

end module mod_codes
!> @}

