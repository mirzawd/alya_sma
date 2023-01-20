!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kermod
!> @{
!> @file    mod_ker_detection.f90
!> @author  Guillaume Houzeaux
!> @date    01/09/2014
!> @brief   Automatic detection module
!> @details Automatic detection module
!> @} 
!-----------------------------------------------------------------------
module mod_ker_detection
  use def_kintyp,         only : ip,rp,lg
  use def_master,         only : INOTSLAVE
  use def_master,         only : ISEQUEN
  use def_master,         only : INOTMASTER
  use def_master,         only : ID_NASTIN
  use def_master,         only : cutim
  use def_master,         only : ittim
  use def_master,         only : mmodu
  use def_master,         only : modul
  use def_master,         only : namod,kfl_paral
  use def_master,         only : mem_modul
  use def_master,         only : lun_detec
  use def_master,         only : kfl_modul
  use def_master,         only : itinn
  use def_master,         only : momod
  use def_master,         only : zeror
  use def_master,         only : veloc
  use def_master,         only : intost
  use def_master,         only : namda
  use def_master,         only : kfl_rstar
  use def_domain,         only : ndime,coord,npoin
  use def_domain,         only : lpoty,exnor,ntens
  use def_domain,         only : nbopo,r_dom,c_dom
  use def_domain,         only : lnnod,lnods,nelem
  use def_kermod,         only : kfl_detection
  use def_kermod,         only : detection_length
  use def_kermod,         only : detection_velocity
  use def_kermod,         only : number_event 
  use def_kermod,         only : events_directory
  use mod_ker_vortex,     only : ker_vortex
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_MIN
  use mod_communications, only : PAR_COMM_RANK_AND_SIZE
  use mod_communications, only : PAR_SEND_RECEIVE
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications, only : PAR_BROADCAST
  use mod_gradie,         only : graten
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_memory,         only : memory_size
  use mod_iofile,         only : iofile
  use mod_iofile,         only : iofile_flush_unit
  implicit none
  private
  integer(ip),    pointer :: node_to_consider(:)

  interface ker_detection_doiter
     module procedure ker_detection_doiter_1,&
          &           ker_detection_doiter_2
  end interface ker_detection_doiter
 
  public :: ker_detection_output
  public :: ker_detection_maximum_residual
  public :: ker_detection_doiter
  public :: ker_detection_boundaries
  public :: ker_detection_min_max_value
  public :: ker_detection_open_file
  public :: ker_events_directory_name
  public :: ker_events_particle_not_converged

contains 

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    06/11/2014
  !> @brief   Get event directory
  !> @details Get event directory
  !>
  !----------------------------------------------------------------------

  subroutine ker_events_directory_name()
!    character(8)   :: current_date
!    character(10)  :: current_time
!    character(32)  :: login_name

    events_directory = './events'

    if( INOTSLAVE ) then  
       !call GETLOG(login_name)
       !call DATE_AND_TIME(current_date,current_time)       
       !call mkdir(adjustl(trim(events_directory)))
       !call opendir("./",adjustl(trim(events_directory)))
       !
       ! Intrinsic 2008 FORTRAN: if it does not compile, comment it and
       ! use the next option: events_directory = './events'
       !
       !events_directory = './events_'//adjustl(trim(login_name))//'_'//adjustl(trim(current_date))//'_'&
       !     //adjustl(trim(current_time(1:6)))             
       !call execute_command_line('mkdir -p '//adjustl(trim(events_directory)))
    end if

  end subroutine ker_events_directory_name

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    06/11/2014
  !> @brief   Open events file
  !> @details Open events file:
  !>          username_date_time/problemname-events.log
  !>
  !----------------------------------------------------------------------

  subroutine ker_detection_open_file()
    character(150) :: fil_detec

    if( kfl_detection /= 0 .and. INOTSLAVE ) then
       
       fil_detec = trim(events_directory)//'/'//adjustl(trim(namda))//'-events.log'    ! Unit 33
       if( kfl_rstar == 2 ) then
          call iofile(0_ip,lun_detec,fil_detec,'AUTOMATIC DETECTION','old','formatted','append')
       else
          call iofile(0_ip,lun_detec,fil_detec,'AUTOMATIC DETECTION')
       end if
    end if

  end subroutine ker_detection_open_file

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    02/09/2014
  !> @brief   Detect non-converged modules
  !> @details Checks if convergence has not been achieved
  !>
  !----------------------------------------------------------------------

  subroutine ker_detection_doiter_1(&
       nn1,nn2,v1,v2,kcom1,kcom2,kdime,&
       wevent,wvariable,wcomment)
    integer(ip), intent(in)          :: nn1,nn2,kcom1,kcom2,kdime
    real(rp),    intent(in)          :: v1(*)
    real(rp),    intent(in)          :: v2(*)
    character(*), intent(in)         :: wevent
    character(*), intent(in)         :: wvariable
    character(*), intent(in)         :: wcomment

    call ker_detection_doiter_driver(&
         1_ip,nn1,nn2,v1,v2,kcom1,kcom2,kdime,&
         wevent,wvariable,wcomment)

  end subroutine ker_detection_doiter_1

  subroutine ker_detection_doiter_2(&
       nn1,nn2,v1,v2,kcom1,kcom2,kdime,&
       wevent,wvariable,wcomment)
    integer(ip),  intent(in)          :: nn1,nn2,kcom1,kcom2,kdime
    real(rp),     intent(in), pointer :: v1(:)
    real(rp),     intent(in), pointer :: v2(:,:,:)
    character(*), intent(in)          :: wevent
    character(*), intent(in)          :: wvariable
    character(*), intent(in)          :: wcomment
    real(rp)                          :: v1_tmp(2),v2_tmp(2)

    if( memory_size(v1) > 0 .and. memory_size(v2) > 0 ) then
       call ker_detection_doiter_driver(&
            1_ip,nn1,nn2,v1,v2,kcom1,kcom2,kdime,&
            wevent,wvariable,wcomment)
    else if( memory_size(v1) > 0 ) then
       call ker_detection_doiter_driver(&
            1_ip,nn1,nn2,v1,v2_tmp,kcom1,kcom2,kdime,&
            wevent,wvariable,wcomment)
    else
       call ker_detection_doiter_driver(&
            1_ip,nn1,nn2,v1_tmp,v2_tmp,kcom1,kcom2,kdime,&
            wevent,wvariable,wcomment)       
    end if       

  end subroutine ker_detection_doiter_2

  !subroutine ker_detection_doiter_all()
  !  integer(ip) :: nn1,nn2,kcom1,kcom2,kdime
  !  real(rp)    :: v1(2)
  !  real(rp)    :: v2(2)
  !  character(*), intent(in) :: wevent
  !  character(*), intent(in) :: wvariable
  !  character(*), intent(in) :: wcomment    
  !  call ker_detection_doiter_driver(&
  !       2_ip,nn1,nn2,v1,v2,kcom1,kcom2,kdime,&
  !       wevent,wvariable,wcomment)
  !end subroutine ker_detection_doiter_all

  subroutine ker_detection_doiter_driver(&
       itask,nn1,nn2,v1,v2,kcom1,kcom2,kdime,&
       wevent,wvariable,wcomment)
    integer(ip),  intent(in) :: itask
    integer(ip),  intent(in) :: nn1,nn2,kcom1,kcom2,kdime
    real(rp),     intent(in) :: v1(*)
    real(rp),     intent(in) :: v2(*)
    character(*), intent(in) :: wevent
    character(*), intent(in) :: wvariable
    character(*), intent(in) :: wcomment
    real(rp),     pointer    :: primitive_variable2(:,:)
    real(rp),     pointer    :: primitive_variable3(:,:,:)
    integer(ip)              :: ipoin_max,ndofn,ipoin,idofn
    integer(ip)              :: idime,i1,i2
    real(rp)                 :: resid,resid_max,vo,va

    if( kfl_detection /= 0 ) then

       select case ( itask )

       case ( 1_ip )
          !
          ! Consider only current module
          !
          if( itinn(modul) >= momod(modul) % miinn .and. momod(modul) % miinn > 1 ) then

             resid_max = -1.0_rp
             ipoin_max = 0
             if( INOTMASTER ) then
                do ipoin = 1,npoin
                   i1    = (ipoin-1) * nn1 + kcom1
                   i2    = (ipoin-1) * nn2 + kcom2
                   resid = 0.0_rp
                   do idime = 0,kdime-1
                      vo    = v2(i2+idime)
                      va    = v1(i1+idime)
                      resid = resid + (vo-va)*(vo-va)
                   end do
                   if( resid > resid_max ) then
                      resid_max = resid
                      ipoin_max = ipoin
                   end if
                end do
             end if
             call ker_detection_maximum_residual(ipoin_max,resid_max,nn1,v1,wevent,wvariable,wcomment)

          end if

       case ( 2_ip )
          !
          ! Loop over all modules
          !
          do modul = 1,mmodu

             if( kfl_modul(modul) /= 0 ) then
                !
                ! Nullify pointers
                !
                nullify(primitive_variable2)
                nullify(primitive_variable3)
                !
                ! Determine primitive variable
                !
                if( modul == ID_NASTIN ) then
                   primitive_variable3 => veloc
                end if
                !
                ! For now just nastin
                !
                if( modul == ID_NASTIN ) then
                   if( itinn(modul) >= momod(modul) % miinn .and. momod(modul) % miinn > 1 ) then

                      if( INOTMASTER ) then
                         ipoin_max =  0
                         resid_max = -1.0_rp
                         if( associated(primitive_variable3) ) then
                            ndofn = size(primitive_variable3,1)
                            do ipoin = 1,npoin
                               resid = 0.0_rp
                               do idofn = 1,ndofn
                                  resid = resid + (primitive_variable3(idofn,ipoin,1)-primitive_variable3(idofn,ipoin,3))**2
                               end do
                               if( resid > resid_max ) then
                                  resid_max = resid
                                  ipoin_max = ipoin
                               end if
                            end do
                         else 
                            ndofn = 1
                            do ipoin = 1,npoin
                               resid = (primitive_variable2(ipoin,1)-primitive_variable2(ipoin,3))**2
                               if( resid > resid_max ) then
                                  resid_max = resid
                                  ipoin_max = ipoin
                               end if
                            end do
                         end if
                      end if
                      call ker_detection_maximum_residual(ipoin_max,resid_max,ndofn,primitive_variable2,wevent,wvariable,wcomment)
                   end if

                end if
             end if
          end do

       end select

    end if

  end subroutine ker_detection_doiter_driver

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    02/09/2014
  !> @brief   Output the position of a maximum residual
  !> @details Slaves have computed a maximum residual at a node. 
  !>          This routine does the following:
  !>          1. Identify the maximum over the slaves and who is the slave
  !>             that owns the maximum,
  !>          2. The responsible slave sends the coordinate of the
  !>             node where it happens to the master. 
  !>          3. The master outputs this information.
  !>
  !----------------------------------------------------------------------

  subroutine ker_detection_maximum_residual(&
       max_residual_node,max_residual,ndofn,xarray,&
       wevent,wvariable,wcomment)
    integer(ip),  intent(in)    :: max_residual_node
    real(rp),     intent(inout) :: max_residual
    integer(ip),  intent(in)    :: ndofn
    real(rp),     intent(in)    :: xarray(ndofn,*)
    character(*), intent(in)    :: wevent
    character(*), intent(in)    :: wvariable
    character(*), intent(in)    :: wcomment
    real(rp)                    :: xcoor(3)
    integer(4)                  :: rank_max_owner4,my_rank4
    integer(ip)                 :: rank_max_owner
    !
    ! Identify the maximum (max_residual) and the rank who owns it (rank_max_owner4)
    !
    call PAR_MAX(max_residual,'IN MY ZONE',rank_max_owner)
    rank_max_owner4 = int(rank_max_owner,4)
    xcoor          = 0.0_rp

    if( max_residual > 0.0_rp ) then
       !
       ! The responsible rank rank_max_owner broadcasts the coordinates
       !
       call PAR_COMM_RANK_AND_SIZE(my_rank4,wherein='IN MY ZONE')
       if( my_rank4 == rank_max_owner4 ) then
          xcoor(1:ndime) = coord(1:ndime,max_residual_node)
       end if
       call PAR_BROADCAST(3_ip,xcoor,'IN MY ZONE',rank_max_owner)
       call ker_detection_output(adjustl(trim(wevent)),rank_max_owner4,xcoor)
    end if
    !
    ! Postprocess
    !
    if( INOTMASTER ) then
       call ker_detection_submesh(number_event,xcoor,ndofn,xarray,wevent,wvariable,wcomment) 
    end if

  end subroutine ker_detection_maximum_residual

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    02/09/2014
  !> @brief   Particle not converged
  !> @details A non-converged particle has been detected
  !>
  !----------------------------------------------------------------------

  subroutine ker_events_particle_not_converged(rank4,xcoor,ndofn,xarray,xarray_particle) 
    integer(4),     intent(in)           :: rank4
    real(rp),       intent(in)           :: xcoor(*)
    integer(ip),    intent(in)           :: ndofn
    real(rp),       intent(in)           :: xarray(ndofn,*)
    real(rp),       intent(in), optional :: xarray_particle(ndofn)
    character(100)                       :: wvariable
    character(100)                       :: wcomment
    character(100)                       :: wevent
    real(rp)                             :: xx(3)

    xx(1:3)     = 0.0_rp
    xx(1:ndofn) = xcoor(1:ndofn)
    wevent      = 'PARTICLE_NOT_CONVERGED'
    wvariable   = 'ADVEC'
    wcomment    = 'NEWTON HAS NOT CONVERGED OR TIME STEP TENDS TO ZERO'

    call ker_detection_output(wevent,rank4,xcoor)
    if( present(xarray_particle) ) then
       call ker_detection_submesh(number_event,xcoor,ndofn,xarray,&
            trim(wevent),trim(wvariable),trim(wcomment),rank4,xarray_particle) 
    else
       call ker_detection_submesh(number_event,xcoor,ndofn,xarray,&
            trim(wevent),trim(wvariable),trim(wcomment))        
    end if

  end subroutine ker_events_particle_not_converged

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    02/09/2014
  !> @brief   Output the position of a maximum residual
  !> @details Slaves have computed a maximum residual at a node. 
  !>          This routine does the following:
  !>          1. Identify the maximum over the slaves and who is the slave
  !>             that owns the maximum,
  !>          2. The responsible slave sends the coordinate of the
  !>             node where it happens to the master. 
  !>          3. The master outputs this information.
  !>
  !----------------------------------------------------------------------

  subroutine ker_detection_min_max_value(&
       what,max_value,ndofn,xarray,xcoor,wevent,wvariable,wcomment)
    character(*), intent(in)    :: what
    real(rp),     intent(inout) :: max_value
    integer(ip),  intent(in)    :: ndofn
    real(rp),     intent(in)    :: xarray(ndofn,*)
    real(rp),     intent(inout) :: xcoor(3)
    character(*), intent(in)    :: wevent
    character(*), intent(in)    :: wvariable
    character(*), intent(in)    :: wcomment
    integer(4)                  :: rank_max_owner4
    integer(ip)                 :: rank_max_owner
    !
    ! Identify the maximum (max_residual) and the rank who owns it (rank_max_owner4)
    !
    if( what == 'MAX' ) then
       call PAR_MAX(max_value,'IN MY ZONE',rank_max_owner)
    else if( what == 'MIN' ) then
       call PAR_MIN(max_value,'IN MY ZONE',rank_max_owner)
    else
       call runend('DO NOT KNOW WHAT TO DO')
    end if

    rank_max_owner4 = int(rank_max_owner,4)
    !
    ! The responsible rank rank_max_owner broadcasts the coordinates
    !
    call PAR_BROADCAST(3_ip,xcoor,'IN MY ZONE',rank_max_owner)
    call ker_detection_output(adjustl(trim(wevent)),rank_max_owner4,xcoor)
    !
    ! Postprocess
    !
    if( INOTMASTER ) then
       call ker_detection_submesh(number_event,xcoor,ndofn,xarray,&
            trim(wevent),trim(wvariable),trim(wcomment)) 
    end if

  end subroutine ker_detection_min_max_value

  !----------------------------------------------------------------------
  !
  ! Output detection information
  !
  !----------------------------------------------------------------------

  subroutine ker_detection_output(wevent,rank4,xcoor)
    character(*), intent(in) :: wevent
    integer(4),   intent(in) :: rank4
    real(rp),     intent(in) :: xcoor(*)

    number_event = number_event + 1

    if( INOTSLAVE ) then
       if( number_event == 1 ) write(lun_detec,2)
       write(lun_detec,1) cutim,number_event,namod(modul),adjustl(trim(wevent)),rank4,xcoor(1:3)
       call iofile_flush_unit(lun_detec)
    end if

1   format(e13.6,2x,i6,2x,a6,2x,a20,2x,i6,3(2x,e13.6))
2   format('# Time------>',2x,'Event ',2x,'Module',2x,'Detection type----->',2x,&
         & 'Rank->',2x,'Coordinates------------------------------->')

  end subroutine ker_detection_output

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    02/10/2014
  !> @brief   Detect reattachment
  !> @details Detect reattachment
  !>
  !----------------------------------------------------------------------

  subroutine ker_detection_vortex()
    real(rp),    pointer :: coord_vort(:,:)
    logical(lg), pointer :: marked_elements(:)
    !
    ! Mark elements
    !
    nullify(coord_vort)
    nullify(marked_elements)
    call ker_vortex(coord_vort,marked_elements)
    
    !if( INOMASTER ) then
    !   do ielem = 1,nelem
    !      do zelel = pelel(ielem),pelel
    !   end do
    !end if

  end subroutine ker_detection_vortex

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    02/10/2014
  !> @brief   Detect reattachment
  !> @details Detect reattachment
  !>
  !----------------------------------------------------------------------

  subroutine ker_detection_boundaries()
    integer(ip) :: ipoin,idime,jpoin
    integer(ip) :: ibopo,jbopo,izdom
    real(rp)    :: n1,n2,n3,s11,s12,s22,s13,s23,s33
    real(rp)    :: sn1,sn2,sn3,stnn,stx,sty,stz
    real(rp)    :: characteristic_length
    real(rp)    :: characteristic_velocity

    integer(ip) :: min_wall_shear_stress_node
    real(rp)    :: min_wall_shear_stress
    real(rp)    :: max_wall_shear_stress

    integer(ip) :: max_wall_normal_stress_node
    real(rp)    :: max_wall_normal_stress
    real(rp)    :: norm_wall_normal_stress

    real(rp),    pointer :: deformation(:,:)
    real(rp),    pointer :: wall_shear_stress(:,:)
    real(rp),    pointer :: wall_normal_stress(:,:)
    real(rp)             :: t1_dot_t2,n1_dot_n2
    real(rp)             :: wall_velocity,dummr
    integer(ip), save    :: ipass = 0

    if( kfl_detection > 0 .and. associated(veloc) ) then

       nullify(deformation)
       nullify(wall_shear_stress)
       nullify(wall_normal_stress)

       if( INOTMASTER ) then 
          call memory_alloca(mem_modul(1:2,modul),'wall_shear_stress','ker_detection_boundaries',wall_shear_stress,ndime,nbopo)
          call memory_alloca(mem_modul(1:2,modul),'wall_shear_stress','ker_detection_boundaries',wall_normal_stress,ndime,nbopo)
          call memory_alloca(mem_modul(1:2,modul),'deformation',      'ker_detection_boundaries',deformation,ntens,npoin)
          !
          ! Determine which nodes should be taken into account
          ! Avoid corners and nodes with angles greater than 20 degrees
          !
          if( ipass == 0 ) then
             !ipass = 1             
             nullify(node_to_consider)
             call memory_alloca(mem_modul(1:2,modul),'node_to_consider','ker_detection_boundaries',node_to_consider,npoin)
             do ipoin = 1,npoin
                ibopo = lpoty(ipoin)
                if( ibopo > 0 ) then
                   wall_velocity = abs(dot_product(veloc(1:ndime,ipoin,1),exnor(1:ndime,1,ibopo))) 
                else
                   wall_velocity = 0.0_rp
                end if
                if( wall_velocity > 0.0_rp ) then
                   node_to_consider(ipoin) = 1
                   do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                      jpoin = c_dom(izdom)             
                      jbopo = lpoty(jpoin)
                      if( jbopo >= 1 ) then
                         node_to_consider(jpoin) = 1
                      end if
                   end do
                else
                   if( ibopo >= 1 ) then
                      do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                         jpoin = c_dom(izdom)             
                         jbopo = lpoty(jpoin)
                         if( jbopo >= 1 ) then
                            wall_velocity = abs(dot_product(veloc(1:ndime,jpoin,1),exnor(1:ndime,1,jbopo))) 
                            if( wall_velocity > 0.0_rp ) then
                               node_to_consider(ipoin) = 1
                               node_to_consider(jpoin) = 1
                            else
                               n1_dot_n2 = dot_product(exnor(1:ndime,1,ibopo),exnor(1:ndime,1,jbopo))
                               if( n1_dot_n2 < 0.94_rp ) then
                                  node_to_consider(ipoin) = 1
                                  node_to_consider(jpoin) = 1
                               end if
                            end if
                         end if
                      end do
                   end if
                end if
             end do
             call PAR_INTERFACE_NODE_EXCHANGE(node_to_consider,'MAX','IN MY CODE')
          end if
       end if
       !
       ! Maximum velocity
       !
       if( INOTMASTER ) then 
          characteristic_velocity = 0.0_rp
          do ipoin = 1,npoin
             dummr = 0.0_rp
             do idime = 1,ndime
                dummr = dummr + veloc(idime,ipoin,1) * veloc(idime,ipoin,1)
             end do
             if( dummr > characteristic_velocity ) then
                characteristic_velocity = dummr
             end if
          end do
          characteristic_velocity = sqrt(characteristic_velocity) + zeror
          characteristic_length   = 1.0_rp
       end if
       call PAR_MAX(characteristic_velocity,'IN MY CODE')

       if( INOTMASTER ) then 
          !
          ! Compute defromation tensor = ( dui/dxi + duj/dxi )
          !
          call graten(veloc,deformation)
          ! 
          ! Compute wall shear stress
          !
          do ipoin = 1,npoin
             ibopo = lpoty(ipoin)
             if( ibopo >= 1 ) then
                n1  = exnor(1,1,ibopo)
                n2  = exnor(2,1,ibopo)
                s11 = deformation(1,ipoin) 
                s22 = deformation(2,ipoin)
                s12 = deformation(3,ipoin)
                sn1 = s11 * n1 + s12 * n2
                sn2 = s12 * n1 + s22 * n2
                if( ndime == 3 ) then
                   n3   = exnor(3,1,ibopo) 
                   s33  = deformation(4,ipoin)
                   s13  = deformation(5,ipoin)
                   s23  = deformation(6,ipoin)
                   sn1  = sn1 + s13 * n3
                   sn2  = sn2 + s23 * n3
                   sn3  = s13 * n1 + s23 * n2 + s33 * n3
                   stnn = sn1 * n1 + sn2 * n2 + sn3 * n3
                   stx  = sn1 - stnn * n1
                   sty  = sn2 - stnn * n2
                   stz  = sn3 - stnn * n3
                   wall_shear_stress(1,ibopo)  = stx         
                   wall_shear_stress(2,ibopo)  = sty       
                   wall_shear_stress(3,ibopo)  = stz      
                   wall_normal_stress(1,ibopo) = sn1    
                   wall_normal_stress(2,ibopo) = sn2       
                   wall_normal_stress(3,ibopo) = sn3      
                else
                   stnn = sn1 * n1 + sn2 * n2
                   stx  = sn1 - stnn * n1
                   sty  = sn2 - stnn * n2
                   wall_shear_stress(1,ibopo)  = stx         
                   wall_shear_stress(2,ibopo)  = sty       
                   wall_normal_stress(1,ibopo) = sn1    
                   wall_normal_stress(2,ibopo) = sn2    
                end if
             end if
          end do
          !
          ! Look for reattachment
          !
          min_wall_shear_stress       =  huge(1.0_rp)
          min_wall_shear_stress_node  =  0
          max_wall_normal_stress_node =  0
          max_wall_normal_stress      = -huge(1.0_rp)
          do ipoin = 1,npoin
             ibopo = lpoty(ipoin)        
             if( ibopo >= 1 .and. node_to_consider(ipoin) == 0 ) then
                wall_velocity = abs(dot_product(veloc(1:ndime,ipoin,1),exnor(1:ndime,1,ibopo)))            
                if( wall_velocity <= zeror ) then
                   do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                      t1_dot_t2 = 0.0_rp
                      jpoin = c_dom(izdom)             
                      jbopo = lpoty(jpoin)
                      if( jbopo >= 1 .and. node_to_consider(jpoin) == 0 ) then
                         wall_velocity = abs(dot_product(veloc(1:ndime,jpoin,1),exnor(1:ndime,1,jbopo)))            
                         if( wall_velocity <= zeror ) then
                            t1_dot_t2 = dot_product(wall_shear_stress(1:ndime,ibopo),wall_shear_stress(1:ndime,jbopo))
                            if( t1_dot_t2 < min_wall_shear_stress ) then
                               min_wall_shear_stress      = t1_dot_t2
                               min_wall_shear_stress_node = ipoin
                            end if
                         end if
                      end if
                   end do
                   norm_wall_normal_stress  = dot_product(wall_normal_stress(1:ndime,ibopo),wall_normal_stress(1:ndime,ibopo)) &
                        / ( characteristic_velocity / characteristic_length )
                   if( norm_wall_normal_stress > max_wall_normal_stress ) then
                      max_wall_normal_stress      = norm_wall_normal_stress
                      max_wall_normal_stress_node = ipoin
                   end if
                end if
             end if
          end do
!print*,min_wall_shear_stress_node
          call memory_deallo(mem_modul(1:2,modul),'wall_shear_stress', 'ker_detection_boundaries',wall_shear_stress)
          call memory_deallo(mem_modul(1:2,modul),'wall_normal_stress','ker_detection_boundaries',wall_normal_stress)
          call memory_deallo(mem_modul(1:2,modul),'deformation',       'ker_detection_boundaries',deformation)
          call memory_deallo(mem_modul(1:2,modul),'node_to_consider',  'ker_detection_boundaries',node_to_consider)

       end if

       max_wall_shear_stress = -min_wall_shear_stress
       call ker_detection_maximum_residual(&
            min_wall_shear_stress_node,max_wall_shear_stress,ndime,veloc,&
            'REATTACHMENT','WALLS','NO_COMMENT')

    end if

  end subroutine ker_detection_boundaries

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    02/10/2014
  !> @brief   Detect reattachment
  !> @details Output mesh arounf the event
  !>          Event type:
  !>          SOLVER_NOT_CONVERGED => all the same
  !>          INNER_NOT_CONVERGED  => VELOC (streamline), DISPL (move mesh), TEMPE (?), etc.
  !>          NEGATIVE_JACOBIAN    => DISPL (displacement solid), DISPM (displacement mesh, usefiul for Chimera)
  !>
  !----------------------------------------------------------------------

  subroutine ker_detection_submesh(&
       event,xx,ndofn,xarray,wevent,&
       wvariable,wcomment,rank4,&
       xarray_particle) 
    integer(ip),    intent(in)           :: event
    real(rp),       intent(in)           :: xx(3)
    integer(ip),    intent(in)           :: ndofn
    real(rp),       intent(in)           :: xarray(ndofn,*)
    character(*),   intent(in)           :: wevent
    character(*),   intent(in)           :: wvariable
    character(*),   intent(in)           :: wcomment
    integer(4),     intent(in), optional :: rank4
    real(rp),       intent(in), optional :: xarray_particle(*)
    character(250)                       :: fil_submesh
    integer(ip)                          :: lun_submesh
    integer(ip)                          :: ipoin,idime,knode,ielem
    integer(ip)                          :: inode
    integer(4)                           :: kelem4,ielem4,ipoin4
    integer(4)                           :: nelem4,npoin4,kpoin4
    real(rp)                             :: diameter
    integer(4),     pointer              :: mark_node4(:)
    logical(lg),    pointer              :: mark_elem(:)
    integer(ip),    pointer              :: permu_node4(:)
    integer(ip),    pointer              :: permu_elem4(:)
    integer(4)                           :: npoin_marked4
    integer(4)                           :: nelem_marked4
    integer(4)                           :: my_rank4
    !
    ! Write event in file event_#_cpu.bin
    ! # ..... the event number
    ! cpu ... the partition number
    !    
    if( INOTMASTER ) then
       !
       ! Mark nodes
       !
       nelem4 = int(nelem,4)
       npoin4 = int(npoin,4)
       allocate(mark_node4(npoin4))
       npoin_marked4 = 0
       do ipoin = 1,npoin
          diameter = 0.0_rp
          do idime = 1,ndime
             diameter = diameter + (coord(idime,ipoin)-xx(idime))*(coord(idime,ipoin)-xx(idime))
          end do
          if( sqrt(diameter) <= detection_length ) then
             mark_node4(ipoin) = 1_4
             npoin_marked4     = 1_4
          else
             mark_node4(ipoin) = 0_4
          end if
       end do
       if( npoin_marked4 == 1 ) then
          !
          ! Mark elements
          !
          allocate(mark_elem(nelem))
          nelem_marked4 = 0_4
          do ielem = 1,nelem
             mark_elem(ielem) = .false.
             loop_inode: do inode = 1,lnnod(ielem)
                ipoin4 = int(lnods(inode,ielem),4)
                if( mark_node4(ipoin4) == 1_4 ) then
                   do knode = 1,lnnod(ielem)
                      kpoin4 = int(lnods(knode,ielem),4)
                      if( mark_node4(kpoin4) == 0_4 ) mark_node4(kpoin4) = 2_4
                   end do
                   nelem_marked4 = nelem_marked4 + 1_4
                   mark_elem(ielem) = .true.
                   exit loop_inode
                end if
             end do loop_inode
          end do
          !
          ! Renumber nodes
          !
          npoin_marked4 = 0_4
          do ipoin4 = 1_4,int(npoin,4)
             if( mark_node4(ipoin4) /= 0_4 ) then
                npoin_marked4 = npoin_marked4 + 1_4
                mark_node4(ipoin4) = npoin_marked4
             end if
          end do
          !
          ! Permutation
          !
          allocate(permu_node4(npoin_marked4))
          allocate(permu_elem4(nelem_marked4))
          npoin_marked4 = 0
          do ipoin4 = 1,int(npoin,4)
             if( mark_node4(ipoin4) /= 0_4 ) then
                npoin_marked4 = npoin_marked4 + 1_4
                permu_node4(npoin_marked4) = ipoin4
             end if
          end do
          nelem_marked4 = 0_4
          do ielem4 = 1,int(nelem,4)
             if( mark_elem(ielem4) ) then
                nelem_marked4 = nelem_marked4 + 1_4
                permu_elem4(nelem_marked4) = ielem4
             end if
          end do
          !
          ! Open file
          !
          lun_submesh = 34
          if( ISEQUEN ) then
             fil_submesh = trim(events_directory)//'/'//adjustl(trim(namda))//'_'//trim(intost(number_event))//'_1.bin'             
          else
             fil_submesh = trim(events_directory)//'/'//adjustl(trim(namda))//'_'//trim(intost(number_event))//'_'//trim(intost(kfl_paral))//'.bin'
          end if
          call iofile(0_ip,lun_submesh,trim(fil_submesh),'DETECTION SUBMESH','unknown','unformatted')
          !
          ! Output
          !
          write(lun_submesh) 3_4                                    
          write(lun_submesh) trim(wevent)             ! EVENT TYPE                       
          write(lun_submesh) trim(namod(modul))       ! PHYSICAL MODULE 
          write(lun_submesh) trim(wcomment)           ! COMMENT      
          write(lun_submesh) cutim,int(ittim,4)       ! TIME,TIME STEP NUMBER
          write(lun_submesh) 'SPHER'                                
          write(lun_submesh) xx(1),xx(2),xx(3),detection_length     
          write(lun_submesh) int(ndime,4),npoin_marked4,nelem_marked4        
          write(lun_submesh) (coord(1:ndime,permu_node4(kpoin4)),kpoin4=1,npoin_marked4)             
          write(lun_submesh) (int(lnnod(permu_elem4(kelem4)),4),mark_node4(lnods(1:lnnod(permu_elem4(kelem4)),permu_elem4(kelem4))),kelem4=1,nelem_marked4) 
          write(lun_submesh) int(ndofn,4),wvariable
          write(lun_submesh) (xarray(1:ndofn,permu_node4(kpoin4)),kpoin4=1,npoin_marked4)  
          !
          ! Close unit
          !
          deallocate(permu_node4)
          deallocate(permu_elem4)
          deallocate(mark_elem)
          call iofile(2_ip,lun_submesh,trim(fil_submesh),'DETECTION SUBMESH')

       end if

       if( present(xarray_particle) ) then
          !
          ! Output
          !
          call PAR_COMM_RANK_AND_SIZE(my_rank4,wherein='IN MY ZONE')
          if( my_rank4 == rank4 ) then
             lun_submesh = 34
             fil_submesh = trim(events_directory)//'/'//adjustl(trim(namda))//'_'//trim(intost(number_event))//'_0.bin'   
             call iofile(0_ip,lun_submesh,trim(fil_submesh),'DETECTION SUBMESH','unknown','unformatted')
             
             write(lun_submesh) 3_4                                    
             write(lun_submesh) trim(wevent)             ! EVENT TYPE                       
             write(lun_submesh) trim(namod(modul))       ! PHYSICAL MODULE 
             write(lun_submesh) trim(wcomment)           ! COMMENT      
             write(lun_submesh) cutim,int(ittim,4)       ! TIME,TIME STEP NUMBER
             write(lun_submesh) 'SPHER'                                
             write(lun_submesh) xx(1),xx(2),xx(3),detection_length     
             write(lun_submesh) int(ndime,4),1_4,1_4
             write(lun_submesh) xx(1:ndime)
             write(lun_submesh) 1_4,1_4
             write(lun_submesh) int(ndofn,4),wvariable
             write(lun_submesh) xarray_particle(1:ndofn)

             call iofile(2_ip,lun_submesh,trim(fil_submesh),'DETECTION SUBMESH')
          end if

       end if

       deallocate(mark_node4)

    end if

  end subroutine ker_detection_submesh

end module mod_ker_detection

