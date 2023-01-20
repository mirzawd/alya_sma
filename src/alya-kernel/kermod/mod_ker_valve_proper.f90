!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-------------------------------------------------------------------------------
!> @addtogroup Valve through porosity
!> @{
!> @authors Ane Beatriz Eguzkitza     : beatriz.eguzkitza@bsc.es
!> @authors Alfonso Santiago          : alfonso.santiago@bsc.es
!> @date    Nov, 2020
!> @}
!------------------------------------------------------------------------------
module mod_ker_valve_proper
    ! ----------------------------------------
    use def_kintyp_basic,               only : ip, rp, lg
    ! ----------------------------------------
    implicit none
    ! ----------------------------------------
    public :: ker_valve_proper_update_porosity
    private
        ! ----------------------------------------
        integer(ip), parameter                                :: NODE_SET=1_ip
        integer(ip), parameter                                :: BOUNDARY_SET=2_ip
        integer(ip), parameter                                :: ELEMENT_SET=3_ip
        ! ----------------------------------------
        type valve_param
            logical(lg)                                       :: initialised =.false.
            integer(ip)                                       :: id_module
            integer(ip)                                       :: last_ittim
            integer(ip)                                       :: last_ittimAP
            integer(ip)                                       :: set_P1_type
            integer(ip)                                       :: set_P2_type
            integer(ip)                                       :: set_P1
            integer(ip)                                       :: set_P2
            real(rp)                                          :: mean_press_P1
            real(rp)                                          :: mean_press_P2
            real(rp)                                          :: poros_max
            real(rp)                                          :: delta_P_ref
            real(rp)                                          :: slope
            real(rp)                                          :: w_relax 
            integer(ip)                                       :: n_it_avg
            real(rp), pointer, dimension(:)                   :: set_P1_hist
            real(rp), pointer, dimension(:)                   :: set_P2_hist
            real(rp), pointer, dimension(:)                   :: set_AP_hist
        endtype valve_param
        ! ----------------------------------------
        type(valve_param),  protected                 :: valve(50)
        ! ----------------------------------------
    
    contains
        !-------------------------------------------------------------------------------
        !> @addtogroup Update porosity
        !> @{
        !> @authors Ane Beatriz Eguzkitza     : beatriz.eguzkitza@bsc.es
        !> @authors Alfonso Santiago          : alfonso.santiago@bsc.es
        !> @date    Nov, 2020
        !> @}
        !------------------------------------------------------------------------------
        subroutine ker_valve_proper_update_porosity(prope_ker, imate)
            use def_kermod,             only :  typ_valpr_ker
            use def_master,             only : ittim
            use def_domain,             only : nelem
            use def_domain,             only : ltype,lmate,ngaus,lnnod
            implicit none
            
            type(typ_valpr_ker), intent(out)        :: prope_ker
            integer(ip), intent(in)                 :: imate
            integer(ip)                             :: ielem, igaus
            integer(ip)                             :: pelty, pgaus, pnode
            real(rp)                                :: w_relax, slope, poros_max, delta_P_ref, mp_1, mp_2, AP
            

            if(.not.valve(imate) % initialised) then
                call ker_valve_proper_initialise(prope_ker, imate)
            endif

            call ker_valve_proper_update_pressure(imate)

            ! Better with AP filtering only
            !call ker_valve_proper_average_pressure(imate)

            poros_max = valve(imate) % poros_max
            slope = valve(imate) % slope
            delta_P_ref = valve(imate) % delta_P_ref
            if(ittim.le.1_ip ) then
                w_relax=1.0_rp
            else
                w_relax= valve(imate) % w_relax
            endif

            mp_1=valve(imate) % mean_press_P1
            mp_2=valve(imate) % mean_press_P2

            AP=mp_2-mp_1
            call ker_valve_proper_filter_AP(imate,AP)

            do ielem = 1,nelem
             if( lmate(ielem) == imate ) then 
                pelty = ltype(ielem)
                if( pelty > 0 ) then
                   pgaus = ngaus(pelty)
                   pnode = lnnod(ielem)
                   do igaus =1, pgaus
                     prope_ker % value_ielem(ielem) % a(igaus) = w_relax*poros_max*(1.0_rp+tanh((AP-delta_P_ref)/slope)) &
                      + (1.0_rp-w_relax)*prope_ker % value_ielem(ielem) % a(igaus)
                   end do
                end if 
             end if
            end do


        end subroutine ker_valve_proper_update_porosity
        !-------------------------------------------------------------------------------
        !> @addtogroup Update pressure
        !> @{
        !> @authors Ane Beatriz Eguzkitza     : beatriz.eguzkitza@bsc.es
        !> @authors Alfonso Santiago          : alfonso.santiago@bsc.es
        !> @date    Nov, 2020
        !> @}
        !------------------------------------------------------------------------------
        subroutine ker_valve_proper_update_pressure(imate)
            use def_master,             only : ittim, momod
            use def_domain,             only : lbsec, lesec
            implicit none

            integer(ip), intent(in)     :: imate
            integer(ip)                 :: loc_set_P1, loc_set_P2
            integer(ip)                 :: read_set_P1, read_set_P2
            integer(ip)                 :: idmod
            integer(ip)                 :: i
            real(rp)                    :: mp_1, mp_2

            read_set_P1 = valve(imate) % set_P1
            read_set_P2 = valve(imate) % set_P2
            idmod       = valve(imate) % id_module
            
            if(valve(imate) % set_P1_type .eq. NODE_SET) then
              call runend('node measuring not implemented')
            elseif(valve(imate) % set_P1_type .eq. BOUNDARY_SET) then
              loc_set_P1 = minloc(lbsec,DIM=1,mask=(lbsec==read_set_P1))
              mp_1 = momod(idmod) % postp(1) % vbset(1, loc_set_P1)
            elseif(valve(imate) % set_P1_type .eq. ELEMENT_SET) then
              loc_set_P1 = minloc(lesec,DIM=1,mask=(lesec==read_set_P1))
              mp_1 = momod(idmod) % postp(1) % vbset(8, loc_set_P1)
            else
              call runend('presure measuring option not implemented')
            endif
            
            if(valve(imate) % set_P2_type .eq. NODE_SET) then
              call runend('node measuring not implemented')
            elseif(valve(imate) % set_P2_type .eq. BOUNDARY_SET) then
              loc_set_P2 = minloc(lbsec,DIM=1,mask=(lbsec==read_set_P2))
              mp_2 = momod(idmod) % postp(1) % vbset(1, loc_set_P2)
            elseif(valve(imate) % set_P2_type .eq. ELEMENT_SET) then
              ! NO ESTA MEANP EN VESET
              loc_set_P2 = minloc(lesec,DIM=1,mask=(lesec==read_set_P2))
              mp_2 = momod(idmod) % postp(1) % veset(8, loc_set_P2)
            else
              call runend('presure measuring option not implemented')
            endif


            valve(imate) % mean_press_P1= mp_1
            valve(imate) % mean_press_P2= mp_2


            ! Store new mean pressure values
            if(ittim .gt. valve(imate) % last_ittim) then
                valve(imate) % last_ittim=ittim
                do i=valve(imate) % n_it_avg-1,1,-1
                    valve(imate) % set_P1_hist(i+1) = valve(imate) % set_P1_hist(i)
                    valve(imate) % set_P2_hist(i+1) = valve(imate) % set_P2_hist(i)
                enddo
            endif

            valve(imate) % set_P1_hist(1) = mp_1
            valve(imate) % set_P2_hist(1) = mp_2

        end subroutine ker_valve_proper_update_pressure
        !-------------------------------------------------------------------------------
        !> @addtogroup time averaging
        !> @{
        !> @authors Ane Beatriz Eguzkitza     : beatriz.eguzkitza@bsc.es
        !> @authors Alfonso Santiago          : alfonso.santiago@bsc.es
        !> @date    Nov, 2020
        !> @}
        !------------------------------------------------------------------------------
        subroutine ker_valve_proper_filter_AP(imate,AP)
            use def_master,             only : ittim
            implicit none
            integer(ip), intent(in)     :: imate
            real(rp), intent(inout)     :: AP
            real(rp)                    :: APf
            integer(ip)                 :: i, it_avg_max, imeana, imeanb
            integer(ip)                 :: filt_type
            real(rp), allocatable       :: sortvalsap(:)

            ! Store new mean pressure values
            if(ittim .gt. valve(imate) % last_ittimAP) then
                valve(imate) % last_ittimAP=ittim
                do i=valve(imate) % n_it_avg-1,1,-1
                    valve(imate) % set_AP_hist(i+1) = valve(imate) % set_AP_hist(i)
                enddo
            endif
            valve(imate) % set_AP_hist(1) = AP !Store unfiltered value

            it_avg_max = min(ittim, valve(imate) % n_it_avg)
            allocate(sortvalsap(it_avg_max))


            filt_type=2_ip
            APf=0.0_rp
            select case(filt_type)
            case(1_ip) ! averaging filter
              do  i=1,it_avg_max
                  APf=APf+ valve(imate) % set_AP_hist(i)/real(it_avg_max,rp)
              enddo
            case(2_ip) ! median filter
                  if(it_avg_max.gt.1) then
                    call sorting(valve(imate) % set_AP_hist,sortvalsAP,it_avg_max)

                    if(mod(it_avg_max,2_ip) == 0_ip) then
                        imeana=ceiling(real(it_avg_max,rp)/2.0_rp-0.5_rp,ip)
                        imeanb=ceiling(real(it_avg_max,rp)/2.0_rp+0.5_rp,ip)
                        APf=(sortvalsap(imeana)+sortvalsap(imeanb))/2.0_rp
                    elseif(mod(it_avg_max,2_ip) == 1_ip) then
                        imeana=ceiling(real(it_avg_max,rp)/2.0_rp,ip)
                        APf= sortvalsap(imeana)
                    else
                        call runend('modulus is weird')
                    endif
                  else
                    APf=valve(imate) % set_AP_hist(1)
                  endif
            endselect

            AP=APf

            deallocate(sortvalsap)


            contains
                subroutine sorting(unsortin, sorted, len)
                    implicit none
                        real(rp), dimension(:), pointer, intent(in)  :: unsortin
                        real(rp), dimension(:), allocatable          :: unsort
                        real(rp), dimension(:), intent(out)          :: sorted
                        integer(ip), intent(in)                      :: len
                        real(rp)                                     :: minv
                        integer(ip)                                  :: i, mini

                       allocate(unsort(len))
                       do i=1,len
                            unsort(i)=unsortin(i)
                       enddo
                       sorted = 0.0_rp
                       do i=1, len
                            minv=minval(unsort,1)
                            mini=minloc(unsort,1)!,kind=ip)

                            sorted(i)=minv
                            unsort(mini)=huge(minv)
                       enddo

                end subroutine sorting
            
        end subroutine ker_valve_proper_filter_AP
        !-------------------------------------------------------------------------------
        !> @addtogroup time averaging
        !> @{
        !> @authors Ane Beatriz Eguzkitza     : beatriz.eguzkitza@bsc.es
        !> @authors Alfonso Santiago          : alfonso.santiago@bsc.es
        !> @date    Nov, 2020
        !> @}
        !------------------------------------------------------------------------------
        subroutine ker_valve_proper_average_pressure(imate)
            use def_master,             only : ittim
            implicit none
            integer(ip), intent(in)     :: imate
            real(rp)                    :: mp_1, mp_2
            integer(ip)                 :: i, it_avg_max, imeana, imeanb
            integer(ip)                 :: filt_type
            real(rp), allocatable       :: sortvalsp1(:), sortvalsp2(:)

            it_avg_max = min(ittim, valve(imate) % n_it_avg)
            allocate(sortvalsp1(it_avg_max))
            allocate(sortvalsp2(it_avg_max))

            mp_1=0.0_rp
            mp_2=0.0_rp

            filt_type=2_ip
            select case(filt_type)
            case(1_ip) ! averaging filter
              do  i=1,it_avg_max
                  mp_1=mp_1+ valve(imate) % set_P1_hist(i)/real(it_avg_max,rp)
                  mp_2=mp_2+ valve(imate) % set_P2_hist(i)/real(it_avg_max,rp)
              enddo
            case(2_ip) ! mean filter
                  if(it_avg_max.gt.1) then
                    call sorting(valve(imate) % set_P1_hist,sortvalsp1,it_avg_max)
                    call sorting(valve(imate) % set_P2_hist,sortvalsp2,it_avg_max)
                    
                    if(mod(it_avg_max,2_ip) == 0_ip) then
                        imeana=ceiling(real(it_avg_max,rp)/2.0_rp-0.5_rp,ip)
                        imeanb=ceiling(real(it_avg_max,rp)/2.0_rp+0.5_rp,ip)
                        mp_1=(sortvalsp1(imeana)+sortvalsp1(imeanb))/2.0_rp
                        mp_2=(sortvalsp2(imeana)+sortvalsp2(imeanb))/2.0_rp
                    elseif(mod(it_avg_max,2_ip) == 1_ip) then
                        imeana=ceiling(real(it_avg_max,rp)/2.0_rp,ip)
                        mp_1= sortvalsp1(imeana)
                        mp_2= sortvalsp2(imeana)
                    else
                        call runend('modulus is weird')
                    endif
                  else
                    mp_1=valve(imate) % set_P1_hist(1)
                    mp_2=valve(imate) % set_P2_hist(1)
                  endif
            endselect

            valve(imate) % mean_press_P1= mp_1
            valve(imate) % mean_press_P2= mp_2
            deallocate(sortvalsp1)
            deallocate(sortvalsp2)

            contains
                subroutine sorting(unsortin, sorted, len)
                    implicit none
                        real(rp), dimension(:), pointer, intent(in)  :: unsortin
                        real(rp), dimension(:), allocatable          :: unsort
                        real(rp), dimension(:), intent(out)          :: sorted
                        integer(ip), intent(in)                      :: len
                        real(rp)                                     :: minv
                        integer(ip)                                  :: i, mini

                       allocate(unsort(len))
                       do i=1,len
                            unsort(i)=unsortin(i)
                       enddo
                       sorted = 0.0_rp
                       do i=1, len
                            minv=minval(unsort,1)
                            mini=minloc(unsort,1)!,kind=ip)

                            sorted(i)=minv
                            unsort(mini)=huge(minv)
                       enddo

                end subroutine sorting
            
        end subroutine ker_valve_proper_average_pressure
    
        !-------------------------------------------------------------------------------
        !> @addtogroup Initialise valve through porosity
        !> @{
        !> @authors Ane Beatriz Eguzkitza     : beatriz.eguzkitza@bsc.es
        !> @authors Alfonso Santiago          : alfonso.santiago@bsc.es
        !> @date    Nov, 2020
        !> @}
        !------------------------------------------------------------------------------
        subroutine ker_valve_proper_initialise(prope_ker, imate)
            use def_master,         only : mem_modul, modul
            use def_kermod,         only : typ_valpr_ker
            use mod_memory,         only : memory_alloca
            implicit none
            type(typ_valpr_ker), intent(in)     :: prope_ker
            integer(ip), intent(in)             :: imate
            integer(ip)                         :: ilaws

            if(imate.gt.50_ip) call runend('valve cant accept more than 50 materials. change the vector size or allocate dinamically')
            if(valve(imate) % initialised) return

           ilaws = prope_ker % ilaws(imate)

            valve(imate) % id_module   = prope_ker % llaws(ilaws) % lresp(1)
            valve(imate) % last_ittim  = -1_ip
            valve(imate) % last_ittimAP  = -1_ip
            valve(imate) % set_P1      = int(prope_ker % rlaws(1,imate),ip)
            valve(imate) % set_P2      = int(prope_ker % rlaws(2,imate),ip)
            valve(imate) % poros_max   =     prope_ker % rlaws(3,imate) 
            valve(imate) % slope       =     prope_ker % rlaws(4,imate)
            valve(imate) % delta_P_ref =     prope_ker % rlaws(5,imate)
            valve(imate) % w_relax     =     prope_ker % rlaws(6,imate)
            valve(imate) % n_it_avg    = int(prope_ker % rlaws(7,imate),ip)
            valve(imate) % set_P1_type = int(prope_ker % rlaws(8,imate),ip)
            valve(imate) % set_P2_type = int(prope_ker % rlaws(9,imate),ip)

            if (valve(imate) % set_P1_type .eq. 0_ip) then
                valve(imate) % set_P1_type = 2_ip
            elseif (valve(imate) % set_P1_type .gt. 3_ip .or. valve(imate) % set_P1_type .lt. 0_ip) then
                call runend('VALVE MEASURING SET NOT PROGRAMMED FOR P1')
            endif
            if (valve(imate) % set_P2_type .eq. 0_ip) then
                valve(imate) % set_P2_type = 2_ip
            elseif (valve(imate) % set_P2_type .gt. 3_ip .or. valve(imate) % set_P2_type .lt. 0_ip) then
                call runend('VALVE MEASURING SET NOT PROGRAMMED FOR P1')
            endif

            if( valve(imate) % w_relax .eq. 0.0_rp) then
                valve(imate) % w_relax =1.0_rp
            endif

            valve(imate) % n_it_avg = max( valve(imate) % n_it_avg, 1_ip)

            call memory_alloca(mem_modul(1:2,modul),'P1_hist','history_P1',valve(imate) % set_P1_hist, valve(imate) % n_it_avg)
            call memory_alloca(mem_modul(1:2,modul),'P2_hist','history_P2',valve(imate) % set_P2_hist, valve(imate) % n_it_avg)
            call memory_alloca(mem_modul(1:2,modul),'AP_hist','history_AP',valve(imate) % set_AP_hist, valve(imate) % n_it_avg)


            valve(imate) % initialised = .true.



        end subroutine ker_valve_proper_initialise

endmodule mod_ker_valve_proper
