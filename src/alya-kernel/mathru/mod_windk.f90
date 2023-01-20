!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!> @defgroup Module for generic windkessel functions
!> Subroutine to implement and use Windkessel functions
!> the module allows to use different modules as input
!> and output
!> @{
!> @name    Module for generic windkessel functions
!> @file    mod_windk.f90
!> @author  Alfonso Santiago
!> @brief    Module for generic windkessel functions
!> @details  Module for generic windkessel functions
!
!-----------------------------------------------------------------------

module mod_windk
  use def_kintyp,                  only : lg,ip,rp, typ_windk_system
  use def_master,                  only: ID_NASTIN, ID_EXMEDI, ID_NASTAL, ID_ALEFOR, ID_SOLIDZ
  use def_master,                  only : INOTMASTER
  use def_master,                  only : cutim, dtinv, dtime
  use def_master,                  only : ittim, itinn, itcou
  use def_master,                  only : mem_modul,modul
  use def_kermod,                  only : max_windk_systems, max_windk_params, number_windk_systems, windk_systems
  use mod_memory,                  only : memory_alloca
  use mod_memory,                  only : memory_deallo
  use mod_parall,                  only : PAR_MY_CODE_RANK
  use mod_communications,          only : PAR_SUM
  use mod_communications,          only : PAR_BARRIER
  use mod_ker_discrete_function,   only : ker_discrete_function


  implicit none



  public  :: mod_windk_create_system
  public  :: mod_windk_interact
  public  :: mod_windk_systems_par_exchange
  private :: compute_function
  private :: mod2int
  private :: retrieve_sysid
  private :: sysid2index
  private :: get_windkessel_parameter

  private
    real(rp), dimension(max_windk_params) :: params_from_function        !calculated parameters form the discrete function
    real(rp)                              :: params_from_function_cutim = -1.0_rp  !time at which the parameters were calculated

  contains

  function get_windkessel_parameter(ifunc, iparam)
    integer(ip),  intent(in)                  :: ifunc
    integer(ip),  intent(in)                  :: iparam
    real(rp)                                  :: get_windkessel_parameter

    if( windk_systems(ifunc) % discrete_function > 0 ) then
        if ( abs(params_from_function_cutim - cutim) > epsilon(1.0_rp) ) then
            params_from_function = 0.0_rp
            call ker_discrete_function( windk_systems(ifunc) % discrete_function, cutim, params_from_function(1:windk_systems(ifunc) % nparam))
            params_from_function_cutim = cutim !this one is not read from restart
        end if
        get_windkessel_parameter = windk_systems(ifunc) % params(iparam) * params_from_function(iparam)
    else
        get_windkessel_parameter = windk_systems(ifunc) % params(iparam)
    end if

  end function get_windkessel_parameter


   !-----------------------------------------------------------------------
   !>
   !> @author  Alfonso Santiago
   !> @date    2019-12-17
   !> @brief   Create windkessel system
   !> @details Assign the variables to the windkessel system, for it
   !>          for it to be used in the windkessel solver
   !>
   !>           used by /domain/reacod.f90
   !>
   !-----------------------------------------------------------------------
   subroutine mod_windk_create_system(sysname, ID_MODUL_IN, tag_in, OPT_ID_MODUL_OUT, OPT_tag_out, OPT_iflow_nsi)
     implicit none

     !Mandatory arguments
     integer(ip), intent(in)           :: ID_MODUL_IN
     character(5), intent(in)          :: sysname
     integer(ip), intent(in)           :: tag_in
     
     ! Optional arguments
     integer(ip), intent(in), optional :: OPT_ID_MODUL_OUT
     integer(ip), intent(in), optional :: OPT_tag_out
     integer(ip), intent(in), optional :: OPT_iflow_nsi

     !Internal variables
     integer(ip)                       :: ID_MODUL_OUT
     integer(ip)                       :: tag_out
     integer(ip)                       :: ifunc
     integer(ip)                       :: funcn_tmp
     integer(ip)                       :: iflow_nsi
     integer(ip), save                 :: extra_systems = 1_ip

     ! Deal with optional arguments
     if(present(OPT_ID_MODUL_OUT)) then
        ID_MODUL_OUT=OPT_ID_MODUL_OUT
        if(.not.present(OPT_tag_out)) call runend('MOD_WINDK CREATE: IF AN OUTPUT MODULE IS PRESENT, AN OUTPUT TAG SHOULD BE PRESENT')
     else
        ID_MODUL_OUT=ID_MODUL_IN
     endif

     if(present(OPT_tag_out)) then
       tag_out=OPT_tag_out
     else
       tag_out=tag_in
     endif

     if(present(OPT_iflow_nsi)) then
       iflow_nsi=OPT_iflow_nsi
     else
       iflow_nsi=-1_ip
     endif


     ! Recover function index from name
     do ifunc = 1,number_windk_systems
       if( trim(sysname) == trim(windk_systems(ifunc) % name) ) then
           funcn_tmp = ifunc
       end if
     end do


     if(windk_systems(funcn_tmp) % sysid .ne. 0_ip)then
        funcn_tmp = number_windk_systems + extra_systems
        extra_systems=extra_systems+1_ip
     endif

        
     windk_systems(funcn_tmp) % sysid     = funcn_tmp
     windk_systems(funcn_tmp) % ID_IN     = ID_MODUL_IN
     windk_systems(funcn_tmp) % ID_OUT    = ID_MODUL_OUT
     windk_systems(funcn_tmp) % tag_in    = tag_in
     windk_systems(funcn_tmp) % tag_out   = tag_out
     windk_systems(funcn_tmp) % iflow_nsi = iflow_nsi

   end subroutine mod_windk_create_system

   !-----------------------------------------------------------------------
   !>
   !> @author  Alfonso Santiago
   !> @date    2019-12-17
   !> @brief   Distribute values of the windkessel system among the slaves
   !> @details Distribute values of the windkessel system among the slaves
   !>
   !>
   !-----------------------------------------------------------------------
   subroutine mod_windk_systems_par_exchange()
     implicit none
     integer(ip)                       :: ifunc


     do ifunc = 1,number_windk_systems
       call PAR_SUM( windk_systems(ifunc) % sysid     ,"IN MY CODE",INCLUDE_ROOT=.true.)
       !wdks_model already exchanged
       !nparam     already exchanged
       !params     already exchanged
       call PAR_SUM( windk_systems(ifunc) % ID_IN     ,"IN MY CODE",INCLUDE_ROOT=.true.) 
       call PAR_SUM( windk_systems(ifunc) % ID_OUT    ,"IN MY CODE",INCLUDE_ROOT=.true.)   
       call PAR_SUM( windk_systems(ifunc) % tag_in    ,"IN MY CODE",INCLUDE_ROOT=.true.) 
       call PAR_SUM( windk_systems(ifunc) % tag_out   ,"IN MY CODE",INCLUDE_ROOT=.true.) 
       !ndxs        already exchanged
       call PAR_SUM( windk_systems(ifunc) % iflow_nsi ,"IN MY CODE",INCLUDE_ROOT=.true.) 
     end do



!     do ifunc = 1,number_windk_systems
!       write(6,*) 'in wdkpar', PAR_MY_CODE_RANK, windk_systems(ifunc) % sysid, windk_systems(ifunc) % name, windk_systems(ifunc) % wdks_model, windk_systems(ifunc) % nparam, windk_systems(ifunc) % ID_IN, windk_systems(ifunc) % ID_OUT,  windk_systems(ifunc) % tag_in , windk_systems(ifunc) % tag_out,  windk_systems(ifunc) % ndxs, windk_systems(ifunc) % iflow_nsi
!     end do

   end subroutine mod_windk_systems_par_exchange


   !-----------------------------------------------------------------------
   !>
   !> @author  Alfonso Santiago
   !> @date    2019-12-17
   !> @brief   
   !> @details 
   !>
   !>
   !>
   !>
   !-----------------------------------------------------------------------
   subroutine mod_windk_interact( &
                 & x_in, y_out, mod_in, tag_in, IS_EXPLICIT, &
                 & OPT_mod_out, OPT_tag_out, OPT_dx_in, OPT_ddx_in, OPT_direction)
        implicit none

        ! Mandatory arguments
        real(rp), intent(in)                              :: x_in
        real(rp), intent(out)                             :: y_out
        integer(ip), intent(in)                           :: tag_in
        character(6), intent(in)                          :: mod_in
        logical(lg),  intent(in)                          :: IS_EXPLICIT

        ! Optional arguments
        character(6), intent(in), optional               :: OPT_mod_out
        integer(ip), intent(in), optional                :: OPT_tag_out
        real(rp), intent(in), optional                   :: OPT_dx_in, OPT_ddx_in
        character(2), intent(in), optional               :: OPT_direction

        ! Internal variables
        integer(ip)                                      :: ID_IN, ID_OUT
        integer(ip)                                      :: tag_out
        real(rp)                                         :: dx, ddx
        integer(ip)                                      :: sysid
        integer(ip)                                      :: action
        integer(ip)                                      :: findex
        character(3)                                     :: time_scheme


        ID_IN=mod2int(mod_in)

        if(present(OPT_mod_out)) then
          ID_OUT=mod2int(OPT_mod_out)
          if(.not.present(OPT_tag_out)) call runend('MOD_WINDK: IF AN OUTPUT MODULE IS PRESENT, AN OUTPUT TAG SHOULD BE PRESENT')
        else
          ID_OUT=ID_IN
        endif

        if(present(OPT_tag_out)) then
          tag_out=OPT_tag_out
        else
          tag_out=tag_in
        endif
          

        if(present(OPT_dx_in)) then
          dx=OPT_dx_in
        endif

        if(present(OPT_ddx_in)) then
          if(.not.present(OPT_dx_in)) call runend('matr/MOD_WINDK: FOR TO THE SECOND DERIVATIVE; THE FIRST SHOULD BE PRESENT')
          ddx=OPT_ddx_in
        endif

        action=0_ip
        if(.not.present(OPT_direction)) then
          action=1_ip ! Bidirectional by default
        else
          if(trim(OPT_direction).eq.'BI')then
            action=1_ip ! Bidirectional
          elseif(trim(OPT_direction).eq.'IN')then
            action=2_ip ! IN
          elseif(trim(OPT_direction).eq.'OU')then
            action=3_ip ! OUT
          else
            call runend('matru: MOD_WINDK: action not recognised. options are IN/OUT/BIDIR')
          endif
        endif

        if ( IS_EXPLICIT ) then
           time_scheme = 'EXP' ! Explicit time-stepping scheme
        else
           time_scheme = 'IMP' ! Implicit time-stepping scheme
        end if

        sysid=retrieve_sysid(ID_IN, tag_in, ID_OUT, tag_out)
        findex=sysid2index(sysid)

        select case(action)
        case(1_ip) ! BIDIR. Read and write in the same call, in the same module.
           !if(INOTMASTER) call compute_function(sysid,x_in,y_out,time_scheme)
           !need to have it on master for restart            
           call compute_function(sysid, x_in, y_out, time_scheme)
        case(2_ip) ! IN. Just receive
           !if(INOTMASTER) call compute_function(sysid, x_in, y_out, time_scheme)
           !need to have it on master for restart            
           call compute_function(sysid, x_in, y_out, time_scheme)
        case(3_ip) ! OUT. Just send
            y_out = windk_systems(findex) % y_out
        case default
            call runend('matru: MOD_WINDK: case not recognised.')
        endselect


   end subroutine mod_windk_interact

   !-----------------------------------------------------------------------
   !>
   !> @author  Alfonso Santiago
   !> @date    2019-12-17
   !> @brief   Solve the windkessel function
   !> @details   Solve the windkessel function using the input X,
   !>            returning the output Y, and using the data of the
   !>            configured function F
   !>
   !-----------------------------------------------------------------------
   subroutine compute_function(f_sysid, x, y, time_scheme)
        implicit none

        !Mandatory variables
        integer(ip),  intent(in)                        :: f_sysid
        real(rp),     intent(in)                        :: x
        real(rp),     intent(out)                       :: y
        character(3), intent(in)                        :: time_scheme

        ! Internal variables
        real(rp), dimension(:), pointer                 :: yprev, xprev
        real(rp)                                        :: r1, c1, r2, l1
        integer(ip)                                     :: findex
        integer(ip)                                     :: i

        real(rp), dimension(2)                          :: dn

       findex=sysid2index(f_sysid)

        windk_systems(findex) % w(2) = windk_systems(findex) % w(1)
        windk_systems(findex) % yrelaxed(2) = windk_systems(findex) % yrelaxed(1)
        windk_systems(findex) % yunrelaxed(2) = windk_systems(findex) % yunrelaxed(1)


        !
        ! Nullify and allocate vectors used in this subroutine
        !
        nullify(yprev)
        nullify(xprev)
        call memory_alloca(mem_modul(1:2,modul),'YPREV','mod_windks', yprev, windk_systems(findex) % ndxs)
        call memory_alloca(mem_modul(1:2,modul),'XPREV','mod_windks', xprev, windk_systems(findex) % ndxs)

        !
        ! give initial values to the previous values to avoid large derivatives
        !
        if( abs( windk_systems(findex) % stored_time_step ) < epsilon(1.0_rp) ) then
            r1 = get_windkessel_parameter( findex, 1_ip )
            windk_systems(findex) % y_out = 0.1_rp*r1*x
            do  i=1_ip, windk_systems(findex) % ndxs
                windk_systems(findex) % yprev(i)= 0.05_rp/real(i,rp) * windk_systems(findex) % y_out
                windk_systems(findex) % xprev(i)= 0.1_rp/real(i,rp) * x
            enddo
        endif
       !
       ! If new time step, save converged values
       !
        if(windk_systems(findex) % stored_time_step .lt. cutim) then
          dn=0.0_rp
          windk_systems(findex) % w=0.0_rp
          windk_systems(findex) % yrelaxed=0.0_rp
          windk_systems(findex) % yunrelaxed=0.0_rp
          
          if(windk_systems(findex) % ndxs .gt. 1_ip) then
            do  i= windk_systems(findex) % ndxs,2_ip,-1_ip
                  windk_systems(findex) % yprev(i) = windk_systems(findex) % yprev(i-1)
                  windk_systems(findex) % xprev(i) = windk_systems(findex) % xprev(i-1)
            enddo
          endif
          windk_systems(findex) % yprev(1) = windk_systems(findex) % y_out
          windk_systems(findex) % xprev(1) = windk_systems(findex) % x_in
          windk_systems(findex) % stored_time_step = cutim
        endif


        windk_systems(findex) % x_in = x

        do  i=1_ip, windk_systems(findex) % ndxs
            yprev(i) = windk_systems(findex) % yprev(i)
            xprev(i) = windk_systems(findex) % xprev(i)
        enddo

        select case ( windk_systems(findex) % wdks_model)
        case(1_ip)
        !md-TODO   ONE_ELEMENT
        !md-TODO
        !md-TODO       x->
        !md-TODO         o------------+
        !md-TODO                      |
        !md-TODO                      \
        !md-TODO                   r1 /
        !md-TODO        y             \
        !md-TODO                      /
        !md-TODO                      |
        !md-TODO         o------------+
        !md-TODO
        !md-TODO   Network's law:
        !md-TODO
        !md-TODO      X(t)=  Y(t)/r1
        !md-TODO
        !md-TODO   In finite differences:
        !md-TODO
        !md-TODO      Y(t)= r1*X(t)
        !md-TODO
            r1 = get_windkessel_parameter(findex, 1_ip)

            y = r1* x 


        case(2_ip)
        !md-TODO   TWO_ELEMENTS
        !md-TODO
        !md-TODO       x->
        !md-TODO         o----------+----------+
        !md-TODO                    |          |
        !md-TODO                    |          \
        !md-TODO                c1 ---      r1 /
        !md-TODO        y          ---         \
        !md-TODO                    |          /
        !md-TODO                    |          |
        !md-TODO         o----------+----------+
        !md-TODO
        !md-TODO   Network's law:
        !md-TODO
        !md-TODO      X(t)= c1*dY/dt + Y/r1
        !md-TODO
        !md-TODO   In finite differences:
        !md-TODO
        !md-TODO      Y^n = \delta t * r1 / (c1*r1 + \delta t)*( X + c1/delta t * Y^{n-1})
        !md-TODO
        !md-TODO   You can check the systems in this web page:
        !md-TODO
        !md-TODO          http://www.civilized.com/mlabexamples/windkesmodel.htmld/
        !md-TODO

            r1 = get_windkessel_parameter( findex, 1_ip )
            c1 = get_windkessel_parameter( findex, 2_ip )

            y= dtime*r1/(r1*c1+dtime) * (x + c1/dtime*yprev(1))

    
        case(3_ip)
        !md-TODO   THREE_ELEMENTS
        !md-TODO
        !md-TODO       x->
        !md-TODO         o--/\/\/\-------+----------+
        !md-TODO             r2          |          |
        !md-TODO                         |          \
        !md-TODO                     c1 ---      r1 /
        !md-TODO        y               ---         \
        !md-TODO                         |          /
        !md-TODO                         |          |
        !md-TODO         o---------------+----------+
        !md-TODO
        !md-TODO
        !md-TODO   Network's law:
        !md-TODO
        !md-TODO      X*(1+r2/r1) + c1r2*dX/dt = Y/r1 + c1*dY/dt
        !md-TODO
        !md-TODO   In finite differences:
        !md-TODO
        !md-TODO      Y^n = dtime*r1/(dtime+r1c1) * [ X^n*(1+r2/r1) + c1r2/dtime*(X^n-X^{n-1}) + c1/dtime*Y^{n-1} ]
        !md-TODO
        !md-TODO   You can check the systems in this web page:
        !md-TODO
        !md-TODO          http://www.civilized.com/mlabexamples/windkesmodel.htmld/
        !md-TODO
            r1 = get_windkessel_parameter( findex, 1_ip )
            c1 = get_windkessel_parameter( findex, 2_ip )
            r2 = get_windkessel_parameter( findex, 3_ip )


            y = dtime*r1/(dtime+r1*c1) * ( x*(1+r2/r1) + c1*r2/dtime*(x-xprev(1)) + c1/dtime*yprev(1) )
            

            !! Re check with equations in:
            !! http://www.civilized.com/mlabexamples/windkesmodel.htmld/

        case(4_ip)
        !md-TODO   FOUR_ELEMENTS
        !md-TODO
        !md-TODO       x->        r2
        !md-TODO         o------/\/\/\-------+----------+
        !md-TODO                             |          |
        !md-TODO                             |          \
        !md-TODO                         c1 ---      r1 /
        !md-TODO        y                   ---         \
        !md-TODO                             |          /
        !md-TODO                 l1          |          |
        !md-TODO         o-------OOOO--------+----------+
        !md-TODO
        !md-TODO
        !md-TODO   Network's law:
        !md-TODO
        !md-TODO      X*(1+r2/r1) + (c1r2+l1/r1)*dX/dt + l1c1*d^2X/dt^2 = Y/r1 + c1*dY/dt
        !md-TODO
        !md-TODO   In finite differences:
        !md-TODO
        !md-TODO      Y^n = dtime*r1/(dtime+r1c1) * [ X^n*(1+r2/r1) + (c1r2+l1/r1)/dtime*(X^n-X^{n-1}) + l1*c1/(dtime*dtime)*(X^n -2*X^{n-1} +X^{n-2})  + c1/dtime*Y^{n-1} ]
        !md-TODO
        !md-TODO   You can check the systems in this web page:
        !md-TODO
        !md-TODO          http://www.civilized.com/mlabexamples/windkesmodel.htmld/
        !md-TODO
        !md-TODO

            r1 = get_windkessel_parameter( findex, 1_ip )
            c1 = get_windkessel_parameter( findex, 2_ip )
            r2 = get_windkessel_parameter( findex, 3_ip )
            l1 = get_windkessel_parameter( findex, 4_ip )

            y = dtime*r1/(dtime+r1*c1) * (  &
                                           x*(1+r2/r1) &
                                         + (c1*r2+l1/r1)/dtime*(x-xprev(1)) &
                                         + l1*c1/(dtime*dtime)*(x-2*xprev(1)+xprev(2)) &  
                                         + c1/dtime*yprev(1)  &
                                         )

        case default
            call runend('mathru/MOD_WINDK: TYPE OF WINDKESSEL FUNCTION NOT FOUND')
        endselect

       !
       ! Aitken to relax the output in case resistance is too large
       !
       if(windk_systems(findex) % w(2).eq.0.0_rp)then
          windk_systems(findex) % w(1)=0.0001_rp
        else
          dn(1) = windk_systems(findex) % yrelaxed(1)-y
          dn(2) = windk_systems(findex) % yrelaxed(2)- windk_systems(findex) % yunrelaxed(1)
          windk_systems(findex) % w(1)= windk_systems(findex) % w(2) + (windk_systems(findex) % w(2)-1.0_rp)*(dn(2)-dn(1))*dn(1)/((dn(2)-dn(1))*(dn(2)-dn(1)))
          windk_systems(findex) % w(1)=min(windk_systems(findex) % w(1),1.0_rp)
        endif
    

        windk_systems(findex) % yunrelaxed(1) = y

        ! If time-stepping scheme is Implicit => use Aitken relaxation
        if ( trim(time_scheme) .eq. 'IMP' ) then
           y = windk_systems(findex) % w(1)*windk_systems(findex) % y_out + (1.0_rp-windk_systems(findex) % w(1))*y
        end if
      
        windk_systems(findex) % yrelaxed(1)=y


       windk_systems(findex) % y_out = y

        call memory_deallo(mem_modul(1:2,modul),'YPREV','mod_windks', yprev)
        call memory_deallo(mem_modul(1:2,modul),'XPREV','mod_windks', xprev)

   end subroutine compute_function

   !-----------------------------------------------------------------------
   !>
   !> @author  Alfonso Santiago
   !> @date    2019-12-17
   !> @brief   Obtain module ID from module name
   !> @details   Obtain module ID from module name. Private function
   !>
   !>
   !-----------------------------------------------------------------------
   integer(ip) function mod2int(mod_in) result(ID)
        implicit none

        character(6), intent(in)   :: mod_in

        select case(mod_in)
        case('NASTIN') ! ID=1
            ID=ID_NASTIN
        case('EXMEDI') ! ID=5
            ID=ID_EXMEDI
        case('NASTAL') ! ID=6
            ID=ID_NASTAL
        case('ALEFOR') ! ID=7
            ID=ID_ALEFOR
        case('SOLIDZ') ! ID=10
            ID=ID_SOLIDZ
        case default
            call runend('mathru/MOD_WINDK: WIKNDKESSEL FOR THIS PHYSICS NOT PROGRAMMED')
        endselect
    endfunction

   !-----------------------------------------------------------------------
   !>
   !> @author  Alfonso Santiago
   !> @date    2019-12-17
   !> @brief   retrieve the sysid of a function from the input data
   !>
   !>
   !-----------------------------------------------------------------------
   integer(ip)  function retrieve_sysid(mod_in, tag_in, mod_out, tag_out) result(foundf)
        implicit none

        ! Mandatory input variables
        integer(ip), intent(in)   :: mod_in, tag_in, mod_out, tag_out

        ! Local variables
        integer(ip)               :: ifunc

       foundf=0_ip
       do ifunc=1, max_windk_systems

         if(windk_systems(ifunc) % ID_IN   .ne. mod_in)  cycle
         if(windk_systems(ifunc) % ID_OUT  .ne. mod_out) cycle
         if(windk_systems(ifunc) % tag_in  .ne. tag_in)  cycle
         if(windk_systems(ifunc) % tag_out .ne. tag_out) cycle

         foundf = windk_systems(ifunc) % sysid
       enddo

       if(foundf.eq.0_ip) call runend('kermod/MOD_WINDK: FUNCTION RETRIEVING FAILED.')

    endfunction

   !-----------------------------------------------------------------------
   !>
   !> @author  Alfonso Santiago
   !> @date    2019-12-17
   !> @brief   retrieve index from sysid
   !>
   !>
   !-----------------------------------------------------------------------
   integer(ip)  function sysid2index(sysid) result(findex)
        implicit none

        integer(ip), intent(in)  :: sysid
        integer(ip)              :: ifunc

        findex=0_ip
        do ifunc=1, max_windk_systems
            if(windk_systems(ifunc) % sysid .eq. sysid) then
                findex=ifunc
                exit
            endif
        enddo
        if(findex.eq.0_ip) call runend('matru: MOD_WINDK: sysid not found!.')
   endfunction sysid2index



endmodule mod_windk
