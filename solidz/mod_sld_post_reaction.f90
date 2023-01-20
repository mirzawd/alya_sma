!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_sld_post_reaction.f90
!> @author  gguillamet
!> @author  aquintanas
!> @date    2021-07-08
!> @brief   Postprocess of reaction forces and displacements
!> @details Postprocess of reaction forces and displacements on a set
!>          of nodes with Dirichlet boundary conditions. The following
!>          methods can be used:
!>
!>          - By plane (axis direction) located at a certain position
!>          - By bounding box
!>
!> @}
!-----------------------------------------------------------------------

module mod_sld_post_reaction

  use def_kintyp_basic, only : ip, rp, lg
  use def_domain,       only : ndime
  use def_solidz,       only : kfl_fixno_sld, frxid_sld
  use def_solidz,       only : kfl_stead_sld
  
  implicit none
  !
  ! Method for postprocess
  !  
  integer(ip), parameter          ::        &
       SLD_PREACTION_PLANE_GREATER =  0_ip, &
       SLD_PREACTION_PLANE_LESS    =  1_ip, &
       SLD_PREACTION_BOX           =  2_ip

  logical(lg),   protected  :: kfl_preac_sld = .False.
  logical(lg),   protected  :: kfl_drofo_sld = .False.
  integer(ip)               :: num_preactset
  real(rp)                  :: dro_percentage
  integer(ip)               :: dro_onset
  real(rp)                  :: ratio_dro
  
  type typ_preact
     integer(ip)         :: method   ! method to select the sets
     real(rp)            :: param(9)
     integer(ip)         :: param_int(9)
     real(rp)            :: displ_value(3)
     real(rp)            :: react_value(3)
     real(rp)            :: react_mag
     real(rp)            :: react_max
  end type typ_preact
          
  type(typ_preact), pointer :: reactset(:)

  private :: sld_post_reaction_memory_allocate
  
  public  :: sld_post_reaction
  public  :: sld_post_reaction_write_res
  public  :: sld_post_reaction_reaous
  public  :: sld_post_reaction_parall

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillamet
  !> @author  aquintanas
  !> @date    2021-07-08
  !> @brief   Memory allocate post reaction
  !> @details Memory allocate post reaction
  !> 
  !-----------------------------------------------------------------------
  
  subroutine sld_post_reaction_memory_allocate(num_sets)

    integer(ip),      intent(in) :: num_sets
    integer(ip)                  :: ii
    
    if( num_sets > 0 .and. .not. associated(reactset) ) then
       allocate(reactset(num_sets))
       do ii = 1,num_sets
          reactset(ii) % method      = 0_ip
          reactset(ii) % param       = 0.0_rp
          reactset(ii) % param_int   = 0_ip
          reactset(ii) % displ_value = 0.0_rp
          reactset(ii) % react_value = 0.0_rp
          reactset(ii) % react_max   = -huge(1.0_rp)
          reactset(ii) % react_mag   = 0.0_rp
       end do
    end if
    
  end subroutine sld_post_reaction_memory_allocate

  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillamet
  !> @author  aquintanas
  !> @date    2021-07-08
  !> @brief   Calculate reactions from a set of nodes  
  !> @details Calculate reactions from a set of nodes   
  !> 
  !-----------------------------------------------------------------------

  subroutine sld_post_reaction

    use def_master,         only : ITER_K
    use def_master,         only : displ
    use def_domain,         only : coord, npoin_own
    use mod_communications, only : PAR_SUM, PAR_MAX, PAR_MIN

    implicit none

    integer(ip)  :: ireac, ipoin, idime, idofn
    
    do ireac = 1,num_preactset
       !
       ! Initializations
       !
       reactset(ireac) % react_value(1:ndime) = 0.0_rp
       reactset(ireac) % displ_value(1:ndime) = -huge(1.0_rp)

       if( reactset(ireac) % method == SLD_PREACTION_PLANE_GREATER ) then
          !
          ! By position greater than a fixed value (one direction)
          !
          do ipoin = 1,npoin_own
             idofn = (ipoin-1)*ndime
             if( coord(reactset(ireac) % param_int(1),ipoin) > reactset(ireac) % param(1) ) then  
                do idime = 1,ndime
                   if( kfl_fixno_sld(idime,ipoin) > 0 ) then
                      reactset(ireac) % react_value(idime) = reactset(ireac) % react_value(idime) - frxid_sld(idofn+idime) 
                      reactset(ireac) % displ_value(idime) = displ(idime,ipoin,ITER_K)
                   end if
                end do
             end if
          end do

       else if( reactset(ireac) % method == SLD_PREACTION_PLANE_LESS ) then
          !
          ! By position less than a fixed value (one direction)
          !
          do ipoin = 1,npoin_own
             idofn = (ipoin-1)*ndime
             if( coord(reactset(ireac) % param_int(1),ipoin) < reactset(ireac) % param(1) ) then  
                do idime = 1,ndime
                   if( kfl_fixno_sld(idime,ipoin) > 0 ) then
                      reactset(ireac) % react_value(idime) = reactset(ireac) % react_value(idime) - frxid_sld(idofn+idime) 
                      reactset(ireac) % displ_value(idime) = displ(idime,ipoin,ITER_K)
                   end if
                end do
             end if
          end do

       else if( reactset(ireac) % method == SLD_PREACTION_BOX ) then
          !
          ! By bounding box
          !
          if( ndime == 2_ip ) then
             do ipoin = 1,npoin_own
                idofn = (ipoin-1)*ndime
                if( (coord(1,ipoin) >= reactset(ireac) % param(1) .and. coord(1,ipoin) <= reactset(ireac) % param(3)) .and. &
                    (coord(2,ipoin) >= reactset(ireac) % param(2) .and. coord(2,ipoin) <= reactset(ireac) % param(4)) ) then   
                   do idime = 1,ndime
                      if( kfl_fixno_sld(idime,ipoin) > 0 ) then
                         reactset(ireac) % react_value(idime) = reactset(ireac) % react_value(idime) - frxid_sld(idofn+idime) 
                         reactset(ireac) % displ_value(idime) = displ(idime,ipoin,ITER_K)
                      end if
                   end do
                end if
             end do
          else
             do ipoin = 1,npoin_own
                idofn = (ipoin-1)*ndime
                if( (coord(1,ipoin) >= reactset(ireac) % param(1) .and. coord(1,ipoin) <= reactset(ireac) % param(4)) .and. &
                    (coord(2,ipoin) >= reactset(ireac) % param(2) .and. coord(2,ipoin) <= reactset(ireac) % param(5)) .and. &
                    (coord(3,ipoin) >= reactset(ireac) % param(3) .and. coord(3,ipoin) <= reactset(ireac) % param(6)) ) then     
                   do idime = 1,ndime
                      if( kfl_fixno_sld(idime,ipoin) > 0 ) then
                         reactset(ireac) % react_value(idime) = reactset(ireac) % react_value(idime) - frxid_sld(idofn+idime) 
                         reactset(ireac) % displ_value(idime) = displ(idime,ipoin,ITER_K)
                      end if
                   end do
                end if
             end do
          end if
       end if
       !
       ! Parall 
       !
       call PAR_SUM(ndime,reactset(ireac) % react_value(:),'IN MY CODE')
       call PAR_MAX(ndime,reactset(ireac) % displ_value(:),'IN MY CODE')

    end do
    !
    ! End of the run due to a drop on the reaction force
    !
    if( kfl_drofo_sld ) call sld_post_reaction_reaction_drop_stop()

  end subroutine sld_post_reaction

  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillamet
  !> @author  aquintanas
  !> @date    2022-10-29
  !> @brief   End of the run when reaction drop is detected
  !> @details End of the run when reaction drop is detected
  !> 
  !-----------------------------------------------------------------------
  
  subroutine sld_post_reaction_reaction_drop_stop

    use def_master,         only : momod
    use def_master,         only : modul
    use def_master,         only : ittim
    use mod_outfor,         only : outfor
    use mod_communications, only : PAR_MAX

    implicit none

    integer(ip) :: ireac
    !
    ! End run due to a drop on the reaction force
    !
    if( ittim > 0 ) then
       do ireac = 1,num_preactset
          if( ireac == dro_onset ) then
             reactset(ireac) % react_mag = dot_product(reactset(ireac) % react_value(:),reactset(ireac) % react_value(:))
             reactset(ireac) % react_mag = sqrt(max(0.0_rp,reactset(ireac) % react_mag))
             ratio_dro = 0.0_rp
             if( reactset(ireac) % react_mag > reactset(ireac) % react_max ) then
                reactset(ireac) % react_max = reactset(ireac) % react_mag
                call PAR_MAX(reactset(ireac) % react_max,'IN MY CODE')
             else
                ratio_dro = (1.0_rp - reactset(ireac) % react_mag/(reactset(ireac) % react_max+epsilon(1.0_rp)))*100.0_rp
                if( ratio_dro > dro_percentage ) then
                   kfl_stead_sld = 1_ip
                   call outfor(28_ip,momod(modul) % lun_outpu,' ')
                end if
             end if
          end if
       end do
    end if

  end subroutine sld_post_reaction_reaction_drop_stop

  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillamet
  !> @author  aquintanas
  !> @date    2021-07-08
  !> @brief   Read data for postprocess reaction sets
  !> @details Read data for postprocess reaction sets   
  !> 
  !-----------------------------------------------------------------------

  subroutine sld_post_reaction_reaous

    use def_inpout, only : words, exists, getint, getrea, getcha
    use mod_ecoute, only : ecoute

    implicit none
     
    external :: runend

    integer(ip) :: ireac    
    !
    ! Initializations
    !
    kfl_preac_sld  = .true.
    num_preactset  = 1
    dro_percentage = 100.0_rp
    dro_onset      = 1
    !
    ! Keywords
    !
    if( exists('NUMBE') ) num_preactset  = getint('NUMBE',1_ip,'#Number of reaction sets')
    if( exists('FORCE') ) then
       kfl_drofo_sld = .True.
       dro_percentage = getrea('FORCE',100.0_rp,'#Reaction drop percentage for the end of the run')
       if( exists('ONREA') ) dro_onset = getint('ONREA',1_ip,'#Reaction set for force drop')
    end if
    call sld_post_reaction_memory_allocate(num_preactset) ! master
    call ecoute('sld_post_reaction_reaous')
    do while(words(1)/='ENDRE')
       if( words(1) == 'REACT') then
          ireac = getint('REACT',1_ip,'#Code for reaction set')
          if(      words(2) == 'PLANE' ) then

             reactset(ireac) % param_int(1) = getint('AXIS ',1_ip,'#Axis')
             if(      exists('GREAT') ) then
                reactset(ireac) % method   = SLD_PREACTION_PLANE_GREATER
                reactset(ireac) % param(1) = getrea('GREAT',0.0_rp,'#Position')
             else if( exists('LESST') ) then
                reactset(ireac) % method   = SLD_PREACTION_PLANE_LESS
                reactset(ireac) % param(1) = getrea('LESST',0.0_rp,'#Position')
             end if
          else if( words(2) == 'BOX  ' ) then

             reactset(ireac) % method = SLD_PREACTION_BOX
             reactset(ireac) % param(1)       = getrea('XMIN ',0.0_rp,'#lower limit in x')
             reactset(ireac) % param(1+ndime) = getrea('XMAX ',0.0_rp,'#upper limit in x')
             if( ndime >=2 ) then
                reactset(ireac) % param(2)       = getrea('YMIN ',0.0_rp,'#lower limit in y')
                reactset(ireac) % param(2+ndime) = getrea('YMAX ',0.0_rp,'#upper limit in y')
             end if
             if( ndime == 3 ) then
                reactset(ireac) % param(3)       = getrea('ZMIN ',0.0_rp,'#lower limit in z')
                reactset(ireac) % param(3+ndime) = getrea('ZMAX ',0.0_rp,'#upper limit in z')
             end if
          else
             call runend('sld_reabcs: Method for OUTPUT,REACTION NOT FOUND!')
          end if
       end if
       call ecoute('sld_post_reaction_reaous')
    end do

  end subroutine sld_post_reaction_reaous

  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillamet
  !> @author  aquintanas
  !> @date    2021-07-08
  !> @brief   Parall send data for postprocess reaction sets
  !> @details Parall send data for postprocess reaction sets   
  !> 
  !-----------------------------------------------------------------------

  subroutine sld_post_reaction_parall()

    use mod_exchange,       only : exchange_init
    use mod_exchange,       only : exchange_add
    use mod_exchange,       only : exchange_end

    integer(ip)              :: ii

       call exchange_init()
       call exchange_add(kfl_preac_sld)
       call exchange_add(kfl_drofo_sld)
       call exchange_add(num_preactset)
       call exchange_add(dro_percentage)
       call exchange_add(dro_onset)   
       call exchange_end()
       call exchange_init()
       if( kfl_preac_sld ) then
          call sld_post_reaction_memory_allocate(num_preactset)
          do ii = 1,num_preactset
             call exchange_add(reactset(ii) % method  )
             call exchange_add(reactset(ii) % param )
             call exchange_add(reactset(ii) % param_int )
          end do
       end if
       call exchange_end()

  end subroutine sld_post_reaction_parall

  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillamet
  !> @author  aquintanas
  !> @date    2021-07-08
  !> @brief   Write results 
  !> @details Write results 
  !> 
  !-----------------------------------------------------------------------

  subroutine sld_post_reaction_write_res(itask)

    use def_master, only : cutim
    use def_master, only : ITER_K, INOTSLAVE
    use def_domain, only : ndime
    use mod_iofile, only : iofile_flush_unit
    use def_solidz, only : lun_react_sld

    implicit none

    integer(ip), intent(in)  :: itask
    integer(ip)              :: ii, i

    select case(itask)

    case( 1_ip )
       !
       ! Write header
       !
       if( INOTSLAVE ) then
          write(lun_react_sld,300)
          if( ndime == 2_ip ) then
             i = 1
             do ii = 1,num_preactset
                write(lun_react_sld,301) ii
                write(lun_react_sld,303, advance="no")
                write(lun_react_sld,302) i+1,'. DISPX',i+2,'. DISPY',i+3,'. REACX',i+4,'. REACY'
                i = i + ndime*2 
             end do
          else
             i = 1
             do ii = 1,num_preactset
                write(lun_react_sld,301) ii
                write(lun_react_sld,303, advance="no")
                write(lun_react_sld,302) i+1,'. DISPX',i+2,'. DISPY',i+3,'. DISPZ',i+4,'. REACX',i+5,'. REACY',i+6,'. REACZ'
                i = i + ndime*2 
             end do
          end if
          write(lun_react_sld,303)
       end if

    case( 2_ip )
       !
       ! Write results
       !
       if( INOTSLAVE ) then
          if( ndime == 2_ip ) then
             write(lun_react_sld,304, advance="no" ) cutim
             do ii = 1,num_preactset
                if( ii == num_preactset ) then
                   write(lun_react_sld,305) &             
                        reactset(ii) % displ_value(1), reactset(ii) % displ_value(2), &            
                        reactset(ii) % react_value(1), reactset(ii) % react_value(2)
                else
                   write(lun_react_sld,305, advance="no") &             
                        reactset(ii) % displ_value(1), reactset(ii) % displ_value(2), &            
                        reactset(ii) % react_value(1), reactset(ii) % react_value(2)
                end if
             end do
          else
             write(lun_react_sld,304, advance="no" ) cutim
             do ii = 1,num_preactset
                if( ii == num_preactset ) then
                   write(lun_react_sld,305) &
                        reactset(ii) % displ_value(1), reactset(ii) % displ_value(2), reactset(ii) % displ_value(3), &
                        reactset(ii) % react_value(1), reactset(ii) % react_value(2), reactset(ii) % react_value(3)
                else
                   write(lun_react_sld,305, advance="no" ) &
                        reactset(ii) % displ_value(1), reactset(ii) % displ_value(2), reactset(ii) % displ_value(3), &
                        reactset(ii) % react_value(1), reactset(ii) % react_value(2), reactset(ii) % react_value(3)
                end if
             end do
          end if
          call iofile_flush_unit(lun_react_sld)
       end if
       
    end select

300 format('# --| ALYA reactions and displacement on user-sets',/,&
         & '# --| Columns displayed:'                   ,/,&
         & '# --|'                                      ,/,&
         & '# --|  1. Current time                        ')
301 format('# --| Set #',i2)
302 format(1x,i2,a7,5(2x,i2,a7))
303 format('# --|')
304 format(2x,e16.8)
305 format(1x,6(2x,e16.8))

  end subroutine sld_post_reaction_write_res

end module mod_sld_post_reaction
