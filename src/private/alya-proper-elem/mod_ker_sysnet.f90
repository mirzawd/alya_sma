!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_ker_sysnet

!-----------------------------------------------------------------!
! Copyright 2021 - 2041 ELEM Biotech.                             !
! Confidential                                                    !
!-----------------------------------------------------------------!

  use def_kintyp_basic, only : ip,rp,lg
  implicit none

    ! unksys (systemic unknowns):
    !  1: V_lv, <----- from left ventricle 
    !  2: V_rv, <----- from right ventricle
    !  3: V_la, 
    !  4: V_ra, 
    !  5: P^system_arterial, 
    !  6: P^system_venous, 
    !  7: P^pulmonary_arterial, 
    !  8: P^pulmonary_venous, 
    !  9: Q^system_arterial, 
    ! 10: Q^system_ven, 
    ! 11: Q^pulmonar_arterial, 
    ! 12: Q^pulmonar_venous, 

    ! funcar (cardiac functions):
    !  1: P_lv, 
    !  2: P_rv, 
    !  3: P_la, 
    !  4: P_ra, 
    !  5: Q_mv, (mitral)
    !  6: Q_av, (aortic)
    !  7: Q_tv, (tricuspid)
    !  8: Q_pv, (pulmonic)

    ! valve resitances
    !      resvalve(1)  <-  mitral
    !      resvalve(2)  <-  aortic 
    !      resvalve(3)  <-  tricuspid 
    !      resvalve(4)  <-  pulmonic 

  
  
  ! imput parameters: in_*

  real(rp) :: &
       in_rmin(4),&
       in_rmax(4),&
       in_clo_factor(4),&
       in_ope_factor(4),&
       in_volv,&
       in_vola,&
       in_vorv,&
       in_vora,&
       in_maxvolven_synth(2),&
       in_minvolven_synth(2),&
       in_elapass,&
       in_erapass,&
       in_elamaxact,&
       in_eramaxact,&
       in_tchamber(4,4),&
       in_csysar,&
       in_csysven,&
       in_cpular,&
       in_cpulven,&
       in_lsysar,&
       in_rsysar,&
       in_lsysven,&
       in_rsysven,&
       in_lpular,&
       in_rpular,&
       in_lpulven,&
       in_rpulven,&
       in_elvpass,&
       in_ervpass,&
       in_elvmaxact,&
       in_ervmaxact,&
       in_lvtimact,&
       in_rvtimact,&
       in_startup_deltime, in_startup_tpreload, in_startup_deltpreload,in_startup_fintime,in_startup_ejection,&
       in_tcycle, in_tisosys, in_tsystole, in_tdiastole, in_tisodia

  integer(ip) :: &
       in_nline_pp, in_startup_cycles, in_startup_coun_print,in_coun_print
  

  real(rp) :: &
       pressure_system_initial(2), pressure_pulmonar_initial(2), dtime_old_sysnet
  
  logical(lg),   protected                :: kfl_sysnet= .False.  
  character(30)                           :: iphase_sysnet(2)

  integer(ip),   protected                :: ncavi = 0_ip
  character(30) :: isolo,ippqp
  character(30) :: iphase_change, ivalve_movement
    
  real(rp) :: rhs(12,4)
  real(rp) :: chron_phase(4,2)      = 0.0_rp ! (sys,isy,dia,idi) x (lv,rv)
  real(rp) :: chron_valve_move(4,2) = 0.0_rp ! (mit,aor,tri,pul) x (open,close)
  real(rp) :: chron_preload                  ! preload chron
  real(rp) :: chron_total, chron_cycle(2)    ! chron_total=total time, chron_cycle= time per cycle
  real(rp) :: resvalve(4)                    ! mit, aor, tri, pul
  real(rp) :: dvoldt(2) = 0.0_rp             ! dV/dt left and right ventricles current dvoldt
  real(rp) :: dvoldt_keep(2,10) = 0.0_rp     ! dV/dt left and right ventricles, stored for up to the last 10 steps
  real(rp) :: pressu_keep(2,10) = 0.0_rp     ! press left and right ventricles, stored for up to the last 10 steps
  real(rp) :: elastan(4)
  real(rp) :: unksys(12)
  real(rp) :: unkaux(12)
  real(rp) :: funcar( 8)
  real(rp) :: funperi

  
  real(rp) :: chamber_elastic_energy( 4 )    ! (la, lv, ra, rv)
  real(rp) :: network_elastic_energy(2,2)    ! (atr, ven) x (sys, pul)
  real(rp) :: network_kinetic_energy(2,2)    ! (atr, ven) x (sys, pul)
  real(rp) :: power_diss_valves( 4 )         ! (mit, aor, tri, pul)
  real(rp) :: power_diss_network(2,2)        ! (atr, ven) x (sys, pul)

  real(rp), pointer :: pressleft(:,:),pressright(:,:)

  ! ----------------------------------------
  type heart_cavity_sysnet
     logical(lg)                         :: initialized = .False.
     integer(ip)                         :: bouset
     integer(ip), dimension(3)           :: plane_nodset
     real(rp),    dimension(3)           :: plane_origin
     real(rp),    dimension(3)           :: plane_normal
     real(rp),    dimension(3)           :: plane_vector
     real(rp),    dimension(3,3)         :: plane_points
     real(rp),    dimension(5)           :: volume
     real(rp),    dimension(5)           :: pres = 0.0_rp
     real(rp),    dimension(3)           :: dvol
     real(rp),    dimension(3)           :: ddvol
     real(rp)                            :: ini_vol
     real(rp)                            :: end_dia_vol
     real(rp)                            :: end_sys_vol
     real(rp),    dimension(3)           :: prest = 0.0_rp
     integer(ip)                         :: phase = 0_ip
  end type heart_cavity_sysnet
  ! ----------------------------------------
  integer(ip),   parameter                :: max_cavities = 4_ip
  type(heart_cavity_sysnet)               :: cavities_sysnet(max_cavities)
  ! ----------------------------------------

contains

  !
  ! Read Alya input files
  !
  subroutine sysnet_read_data()
    use def_master,           only :  intost
    use def_domain,           only :  ndime
    use def_master,           only :  tcardiac_cycle, dtime
    use def_inpout,           only :  words, param, exists, getint, getrea, getcha, nnpar
    use mod_ecoute,           only :  ecoute
    use mod_messages,         only :  livinf, messages_live
    ! -------------------------------
    implicit none
    ! -------------------------------
    integer(ip)                    :: icavi, ii
    real(rp)                       :: p1(3), p2(3), v1(3), nrm
    ! -------------------------------
    
    kfl_sysnet = .True.
    in_startup_coun_print = 1_ip   ! startup: print results every time step
    ivalve_movement= 'sharp'

    ivalve_movement= 'smooth'
    
    call ecoute('sysnet_read_data')

    !CARDIAC_CYCLE section
    do while(words(1)/='ENDSY')
       if( words(1) == 'CAVIT' )then

          ncavi=ncavi+1

          icavi=ncavi
          cavities_sysnet(icavi) % initialized  = .True.
          cavities_sysnet(icavi) % bouset       = 0_ip
          cavities_sysnet(icavi) % plane_nodset = 0_ip
          cavities_sysnet(icavi) % plane_origin = 0.0_rp
          cavities_sysnet(icavi) % plane_normal = 0.0_rp
          cavities_sysnet(icavi) % plane_vector = 0.0_rp
          cavities_sysnet(icavi) % volume       = 0.0_rp

          in_startup_fintime=   10000_rp                  ! default: very long time, but...
          in_startup_deltime=   0.0_rp          
          in_startup_cycles =   -1_ip                     ! ... no startup cycles
          in_startup_tpreload=  -1.0_rp                   ! no preload
          in_startup_deltpreload =   dtime                ! 

          dtime_old_sysnet = dtime
          
          call ecoute('sysnet_read_data')
          do while(words(1)/='ENDCA')

             !
             ! Boundary defining the cavity
             !
             if (words(1)=='BOUND') then 
                cavities_sysnet(icavi) % bouset = int(param(1))

             endif

             !
             ! Plane definition for open cavity
             !
             if(  words(1) == 'PLANE' )then

                if( exists('POINT') )then

                   ! Plane defined by 2 points
                   p1(1:ndime) = param(2:4)
                   p2(1:ndime) = param(5:7)
                   v1(:)= p1(:) - p2(:) 
                   nrm = sqrt(v1(1)**2 + v1(2)**2 + v1(3)**2)
                   if( nrm < 0.0_rp ) call runend('MOD_KER_SYSNET: ||n|| of the plane is <= 0')
                   cavities_sysnet(icavi) % plane_vector(:) = v1(:)/nrm

                elseif( exists('NODSE') .or. exists('SETNO') )then

                  ! Plane defined by a node set
                     if ( nnpar == 2 ) then
                        cavities_sysnet(icavi) % plane_nodset = (/int(param(2),kind=ip), int(param(3),kind=ip), 0_ip/)
                     elseif( nnpar==3 ) then
                        cavities_sysnet(icavi) % plane_nodset = (/int(param(2),kind=ip), int(param(3),kind=ip), int(param(4),kind=ip)/)
                     else  
                        call runend("PLANE NODESET FOR CAVITY "//trim(intost(icavi))//" HAS TO HAVE 2 OR 3 VALUES. "//trim(intost(nnpar))//" VALUES ARE PROVIDED")
                     end if
                elseif( exists('ORIGI') )then

                   ! Plane origin
                   cavities_sysnet(icavi) % plane_origin = param(2:4)

                else
                   ! Plane defined by the normal vector
                   call runend('MOD_KER_SYSNET: PLANE POINT / OIRIGIN / NODSET')

                endif
             endif

             call ecoute('sysnet_read_data')
          end do
       else if( words(1) == 'PARAM' )then
          do while( words(1) /= 'ENDPA' )
             if( words(1) == 'PRINT') then
                in_coun_print = int(param(1))
             else if( words(1) == 'VALVE') then
                do while( words(1) /= 'ENDVA' )
                   if( words(1) == 'RESMI') then
                      in_rmin(1:4) = param(1)
                   else if( words(1) == 'RESMA') then
                      in_rmax(1:4) = param(1)
                   else if( words(1) == 'CLOSI') then                      
                      do while( words(1) /= 'ENDCL' )
                         if( words(1) == 'MITRA') then
                            in_clo_factor(1) = param(1)
                         else if( words(1) == 'AORTI') then
                            in_clo_factor(2) = param(1)
                         else if( words(1) == 'TRICU') then
                            in_clo_factor(3) = param(1)
                         else if( words(1) == 'PULMO') then
                            in_clo_factor(4) = param(1)
                         end if
                         call ecoute('sysnet_read_data')
                      end do
                   else if( words(1) == 'OPENI') then                      
                      do while( words(1) /= 'ENDOP' )
                         if( words(1) == 'MITRA') then
                            in_ope_factor(1) = param(1)
                         else if( words(1) == 'AORTI') then
                            in_ope_factor(2) = param(1)
                         else if( words(1) == 'TRICU') then
                            in_ope_factor(3) = param(1)
                         else if( words(1) == 'PULMO') then
                            in_ope_factor(4) = param(1)
                         end if
                         call ecoute('sysnet_read_data')
                      end do
                   end if
                   
                   call ecoute('sysnet_read_data')
                end do
                ! recompute closed and open valves resistances
                in_rmin(1)= in_rmin(1) * in_ope_factor(1)
                in_rmin(2)= in_rmin(2) * in_ope_factor(2)
                in_rmin(3)= in_rmin(3) * in_ope_factor(3)
                in_rmin(4)= in_rmin(4) * in_ope_factor(4)

                in_rmax(1)= in_rmax(1) * in_clo_factor(1)
                in_rmax(2)= in_rmax(2) * in_clo_factor(2)
                in_rmax(3)= in_rmax(3) * in_clo_factor(3)
                in_rmax(4)= in_rmax(4) * in_clo_factor(4)                
             else if( words(1) == 'ATRIA') then
                do while( words(1) /= 'ENDAT' )
                   if( words(1) == 'BASEV') then
                      in_vola = param(1)
                      in_vora = param(2)
                   else if( words(1) == 'PASEL') then
                      in_elapass = param(1)
                      in_erapass = param(2)
                   else if( words(1) == 'ACTEL') then
                      in_elamaxact = param(1)
                      in_eramaxact = param(2)
                   else if( words(1) == 'TACTA') then
                      in_tchamber(1:4,1) = param(1:4)
                      in_tchamber(1:4,2) = param(5:8)                      
                   end if
                   call ecoute('sysnet_read_data')
                end do
             else if( words(1) == 'CIRCU') then
                do while( words(1) /= 'ENDCI' )
                   if( words(1) == 'CSYSA') then
                      in_csysar = param(1)
                   else if( words(1) == 'LSYSA') then
                      in_lsysar = param(1)
                   else if( words(1) == 'RSYSA') then
                      in_rsysar = param(1)
                   else if( words(1) == 'CSYSV') then
                      in_csysven = param(1)
                   else if( words(1) == 'LSYSV') then
                      in_lsysven = param(1)
                   else if( words(1) == 'RSYSV') then
                      in_rsysven = param(1)
                   else if( words(1) == 'CPULA') then
                      in_cpular = param(1)
                   else if( words(1) == 'LPULA') then
                      in_lpular = param(1)
                   else if( words(1) == 'RPULA') then
                      in_rpular = param(1)
                   else if( words(1) == 'CPULV') then
                      in_cpulven = param(1)
                   else if( words(1) == 'LPULV') then
                      in_lpulven = param(1)
                   else if( words(1) == 'RPULV') then
                      in_rpulven = param(1)

                   end if
                   call ecoute('sysnet_read_data')
                end do
             else if( words(1) == 'TIMES') then
                do while( words(1) /= 'ENDTI' )
                   call ecoute('sysnet_read_data')
                end do
             else if( words(1) == 'INITI') then
                do while( words(1) /= 'ENDIN' )
                   if( words(1) == 'SYSAR') then
                      pressure_system_initial(1) = param(1)
                   else if( words(1) == 'SYSVE') then
                      pressure_system_initial(2) = param(1)
                   else if( words(1) == 'PULAR') then
                      pressure_pulmonar_initial(1) = param(1)
                   else if( words(1) == 'PULVE') then
                      pressure_pulmonar_initial(2) = param(1)
                   else if( words(1) == 'PERIC') then
                      funperi = param(1)
                   end if
                   call ecoute('sysnet_read_data')
                end do

             else if( words(1) == 'START') then
                do while( words(1) /= 'ENDST' )
                   if( words(1) == 'PRINT') then
                      in_startup_coun_print = int(param(1))
                   else if( words(1) == 'DELTI') then
                      in_startup_deltime = param(1)
!                   else if( words(1) == 'FINTI') then
!                      in_startup_fintime = param(1)
                   else if( words(1) == 'CYCLE') then
                      in_startup_cycles = int(param(1))
                   else if( words(1) == 'EJECT') then
                      in_startup_ejection = param(1)
                   else if( words(1) == 'TISOS') then
                      in_tisosys = param(1)
                   else if( words(1) == 'TSYST') then
                      in_tsystole = param(1)
                   else if( words(1) == 'TISOD') then
                      in_tisodia = param(1)
                   else if( words(1) == 'TDIAS') then
                      in_tdiastole = param(1)
                   else if( words(1) == 'TIMEP') then
                      in_startup_tpreload = param(1)
                   else if( words(1) == 'DELTP') then
                      in_startup_deltpreload = param(1)
                   end if
                   call ecoute('sysnet_read_data')
                end do

             end if

             call ecoute('sysnet_read_data')

          end do
       end if
       call ecoute('sysnet_read_data')

    end do
    !
    ! Check if the nodeset is valid
    !
    do icavi= 1, max_cavities
         if( cavities_sysnet(icavi) % initialized )then
            do ii = 1, 2
               if ( cavities_sysnet(icavi) % plane_nodset(ii) == 0_ip ) then
                  call runend("VOLUME CALCULATION FOR CAVITY "//trim(intost(icavi))//" NODE NUMBER "//trim(intost(ii))//" IS UNDEFINED. 2 OR 3 NODES ARE NEEDED")
               end if
            enddo
            if (cavities_sysnet(icavi) % plane_nodset(3) == 0_ip) then
               call messages_live("VOLUME CALCULATION FOR CAVITY "//trim(intost(icavi))//" USES 2 NODES")
            else
               call messages_live("VOLUME CALCULATION FOR CAVITY "//trim(intost(icavi))//" USES 3 NODES")
            end if
         end if   
    enddo
 


    !
    ! Derived parameters
    !
    in_tcycle= tcardiac_cycle
    
    in_tsystole    =  in_tsystole  * in_tcycle
    in_tisosys     =  in_tisosys   * in_tcycle
    in_tdiastole   =  in_tdiastole * in_tcycle
    in_tisodia     =  in_tisodia   * in_tcycle
    
    ! from Blanco 2010, which uses a in_tcycle=1.00
    ! left atrium
    in_tchamber(1,1) = in_tchamber(1,1) * in_tcycle      ! Tac
    in_tchamber(2,1) = in_tchamber(2,1) * in_tcycle      ! Tar
    in_tchamber(3,1) = in_tchamber(3,1) * in_tcycle      ! tac
    in_tchamber(4,1) = in_tchamber(4,1) * in_tcycle      ! tar
    ! right atrium   
    in_tchamber(1,2) = in_tchamber(1,2) * in_tcycle      ! Tac
    in_tchamber(2,2) = in_tchamber(2,2) * in_tcycle      ! Tar
    in_tchamber(3,2) = in_tchamber(3,2) * in_tcycle      ! tac
    in_tchamber(4,2) = in_tchamber(4,2) * in_tcycle      ! tar

    ! Initial atrial elastances
    elastan(3) = in_elapass
    elastan(4) = in_erapass
    
  end subroutine sysnet_read_data

  !
  ! Write startup result files: only for the startup phase
  !
  subroutine sysnet_write_startup_res(&
       iwrite, coun_iter)

    use def_master , only : INOTSLAVE, cutim
    implicit none
    character(30) :: iwrite
    integer(ip)   :: coun_iter, icolu,valve_status(4,2)
    real(rp)      :: rauxi(2)

300 format(30(2x,e12.5))
400 format((a),2x,(a),2x,30(2x,e12.5))

    if (INOTSLAVE) then

       if (adjustl(trim(isolo)) == 'solo') then
          ! ventricular pressures, to be send to alya
          if (adjustl(trim(iwrite)) == 'start') then
             open(10030,file='sysnet-startup-unksys.res',status='unknown')
             open(10040,file='sysnet-startup-funcar.res',status='unknown')
             open(10050,file='sysnet-startup-valves.res',status='unknown')
             open(10060,file='sysnet-startup-elastances.res',status='unknown')
             open(10070,file='sysnet-startup-volumes.res',status='unknown')       
             open(10080,file='sysnet-startup-presdif.res',status='unknown')       
          else

             valve_status = 0
             if ((funcar(1)-funcar(3)) .ge. 0) valve_status(1,1)= 1 ! mitral closed
             if ((unksys(5)-funcar(1)) .ge. 0) valve_status(2,1)= 1 ! aortic closed
             if ((funcar(2)-funcar(4)) .ge. 0) valve_status(3,1)= 1 ! tricuspid closed
             if ((unksys(7)-funcar(2)) .ge. 0) valve_status(4,1)= 1 ! pulmonic closed

             if (adjustl(trim(iphase_sysnet(1))) == 'sys') rauxi(1)=1.0_rp
             if (adjustl(trim(iphase_sysnet(1))) == 'isy') rauxi(1)=2.0_rp
             if (adjustl(trim(iphase_sysnet(1))) == 'dia') rauxi(1)=3.0_rp
             if (adjustl(trim(iphase_sysnet(1))) == 'idi') rauxi(1)=4.0_rp

             if (adjustl(trim(iphase_sysnet(2))) == 'sys') rauxi(2)=1.0_rp
             if (adjustl(trim(iphase_sysnet(2))) == 'isy') rauxi(2)=2.0_rp
             if (adjustl(trim(iphase_sysnet(2))) == 'dia') rauxi(2)=3.0_rp
             if (adjustl(trim(iphase_sysnet(2))) == 'idi') rauxi(2)=4.0_rp

             if (mod(coun_iter,in_startup_coun_print) == 0) then
                write(10030,300) &
                     (unksys(icolu),icolu=1,12), chron_total,chron_cycle,real(coun_iter),cutim
                write(10040,300) &
                     (funcar(icolu),icolu=1,8), chron_total,chron_cycle,real(coun_iter),rauxi,cutim
                write(10050,400) &
                     adjustl(trim(iphase_sysnet(1))),&
                     adjustl(trim(iphase_sysnet(2))),&
                     funcar(1),unksys(1),funcar(2),unksys(2),&
                     ((resvalve(icolu)/in_rmax(icolu)),icolu=1,4), &
                     ((resvalve(icolu)),icolu=1,4), &
                     real(valve_status(1:4,1)), &
                     chron_total,chron_cycle(1),chron_cycle(2),real(coun_iter),cutim
                write(10060,300) &
                     (elastan(icolu),icolu=1,4), chron_total,chron_cycle,real(coun_iter),cutim
                write(10070,300) funcar(1),unksys(1),funcar(2),unksys(2),&
                     dvoldt(1),dvoldt(2), &
                     chron_total,chron_cycle(1),chron_cycle(2),cutim
                write(10080,300) &
                     (funcar(1)-unksys(5)),(funcar(2)-unksys(7)),&
                     (funcar(3)-funcar(1)),(funcar(4)-funcar(2)),funcar(3),funcar(4),&
                     chron_total,chron_cycle(1),chron_cycle(2),cutim
             end if

          end if

       else  if (adjustl(trim(isolo)) == 'alya') then

       end if

    end if
  end subroutine sysnet_write_startup_res
  !
  ! Write alya result files for all the run
  !
  subroutine sysnet_write_alya_res(itask,lun_res)
    use def_master, only            :  ittim, cutim
    use def_master, only            :  TIME_N, INOTSLAVE
    implicit none
    ! -------------------------------
    integer(ip), intent(in)        :: itask, lun_res
    real(rp)                       :: rauxi(2)
    ! -------------------------------

    if (INOTSLAVE) then

       select case(itask)

       case( 1_ip )
          !
          ! header for file: sysnetHeart.sld.res
          !
          write(lun_res,*) '# --| ALYA SYSNET System Network Heart results file'
          write(lun_res,*) '# --| '
          write(lun_res,*) '# --|  1.  Current Iteration'
          write(lun_res,*) '# --|  2.  Current Time'
          write(lun_res,*) '# --|  3.  Phase'
          write(lun_res,*) '# --|  4.  Phase time'
          write(lun_res,*) '# --|  5.  Left  Ventricle Volume'
          write(lun_res,*) '# --|  6.  Left  Atrium Volume'
          write(lun_res,*) '# --|  7.  Left  Ventricle Pressure [mmHg]'
          write(lun_res,*) '# --|  8.  Left  Atrium Pressure [mmHg]'
          write(lun_res,*) '# --|  9.  Left  Mitral valve resistance'
          write(lun_res,*) '# --|  10. Left  Aortic valve resistance'
          write(lun_res,*) '# --|  11. Left  Mitral valve flowrate'
          write(lun_res,*) '# --|  12. Left  Aortic valve flowrate'
          write(lun_res,*) '# --|  13. Right Ventricle Volume'
          write(lun_res,*) '# --|  14. Right Atrium Volume'
          write(lun_res,*) '# --|  15. Right Ventricle Pressure [mmHg]'
          write(lun_res,*) '# --|  16. Right Atrium Pressure [mmHg]'
          write(lun_res,*) '# --|  17. Right Tricuspid valve resistance'
          write(lun_res,*) '# --|  18. Right Pulmonic valve resistance'
          write(lun_res,*) '# --|  19. Right Tricuspid valve flowrate'
          write(lun_res,*) '# --|  20. Right Pulmonic valve flowrate'

       case( 2_ip )
          !
          ! header for file: sysnetSystem.sld.res
          !
          write(lun_res,*) '# --| ALYA SYSNET System Network System results file'
          write(lun_res,*) '# --| '
          write(lun_res,*) '# --|  1.  Current Iteration'
          write(lun_res,*) '# --|  2.  P^system_arterial'
          write(lun_res,*) '# --|  3.  P^system_venous     '
          write(lun_res,*) '# --|  4.  P^pulmonary_arterial'
          write(lun_res,*) '# --|  5.  P^pulmonary_venous  '
          write(lun_res,*) '# --|  6.  Q^system_arterial   '
          write(lun_res,*) '# --|  7.  Q^system_ven        '
          write(lun_res,*) '# --|  8.  Q^pulmonar_arterial '
          write(lun_res,*) '# --|  9.  Q^pulmonar_venous   '

       case( 3_ip )
          !
          ! contents for file: sysnetHeart.sld.res
          !

          if (adjustl(trim(iphase_sysnet(1))) == 'sys') rauxi(1)=1.0_rp
          if (adjustl(trim(iphase_sysnet(1))) == 'isy') rauxi(1)=2.0_rp
          if (adjustl(trim(iphase_sysnet(1))) == 'dia') rauxi(1)=3.0_rp
          if (adjustl(trim(iphase_sysnet(1))) == 'idi') rauxi(1)=4.0_rp

          if (adjustl(trim(iphase_sysnet(2))) == 'sys') rauxi(2)=1.0_rp
          if (adjustl(trim(iphase_sysnet(2))) == 'isy') rauxi(2)=2.0_rp
          if (adjustl(trim(iphase_sysnet(2))) == 'dia') rauxi(2)=3.0_rp
          if (adjustl(trim(iphase_sysnet(2))) == 'idi') rauxi(2)=4.0_rp

          if (mod(ittim,in_coun_print) == 0) then
             write(lun_res, 101)  &
                  ittim,&          ! 1. current iteration
                  cutim,&          ! 2. current time
                  rauxi(1),&       ! 3. phase
                  chron_cycle(1),& ! 4. phase time
                  unksys(1),&      ! 5. lv volume
                  unksys(3),&      ! 6. la volume
                  funcar(1),&      ! 7. lv pressure
                  funcar(3),&      ! 8. la pressure
                  resvalve(1),&    ! 9. mitral resistance
                  resvalve(2),&    ! 10. aortic resistance
                  funcar(5),&      ! 11. mitral flowrate
                  funcar(6),&      ! 12. aortic florate
                  unksys(2),&      ! 13. rv volume
                  unksys(4),&      ! 14. ra volume
                  funcar(2),&      ! 15. rv pressure
                  funcar(4),&      ! 16. ra pressure
                  resvalve(3),&    ! 17. tricuspid resistance 
                  resvalve(4),&    ! 18. pulmonic resistance
                  funcar(7),&      ! 19. tricuspid flowrate
                  funcar(8),&      ! 20. pulmonic flowrate
                  sqrt(dvoldt(1)*dvoldt(1))/unksys(1),&
                  sqrt(dvoldt(2)*dvoldt(2))/unksys(2)
             flush( lun_res)
          end if
       case( 4_ip )
          !
          ! contents for file: sysnetSystem.sld.res
          !
          if (mod(ittim,in_coun_print) == 0) then
             write(lun_res, 101)  &
                  ittim,&          ! 1. current iteration
                  unksys(5:12)     ! system pressures and flow rates
             flush( lun_res)
          end if

       end select

    end if
    ! -------------------------------
101 format(2x, i9, 30(2x,e16.8e3))
    ! -------------------------------
  end subroutine sysnet_write_alya_res

  subroutine sysnet_manage_restart()
    use mod_strings,  only : integer_to_string
    ! -------------------------------
    implicit none
    ! -------------------------------
    ! -------------------------------

  end subroutine sysnet_manage_restart
  

  subroutine sysnet_sendat(itask)
    ! -------------------------------
    use mod_exchange,           only : exchange_add
    implicit none
    integer(ip),  intent(in)       :: itask
    integer(ip)                    :: icavi, i, j
    ! ------------------------------- 

    select case(itask)

    case (1_ip)
       call exchange_add(kfl_sysnet)
       call exchange_add(ncavi)

    case (2_ip)

       ! cavities
       
       do icavi = 1,size(cavities_sysnet)
          call exchange_add(cavities_sysnet(icavi) % initialized )
          call exchange_add(cavities_sysnet(icavi) % phase )
          call exchange_add(cavities_sysnet(icavi) % bouset )
          call exchange_add(cavities_sysnet(icavi) % plane_nodset)
          call exchange_add(cavities_sysnet(icavi) % plane_vector)
          call exchange_add(cavities_sysnet(icavi) % plane_origin)
          call exchange_add(cavities_sysnet(icavi) % plane_normal)
          call exchange_add(cavities_sysnet(icavi) % plane_points )
       end do


       ! input data
       call exchange_add(dtime_old_sysnet)
       call exchange_add(in_volv)
       call exchange_add(in_vola)
       call exchange_add(in_vorv)
       call exchange_add(in_vora)
       call exchange_add(in_maxvolven_synth(1))
       call exchange_add(in_maxvolven_synth(2))
       call exchange_add(in_minvolven_synth(1))
       call exchange_add(in_minvolven_synth(2))
       call exchange_add(in_elapass)
       call exchange_add(in_erapass)
       call exchange_add(in_elamaxact)
       call exchange_add(in_eramaxact)
       call exchange_add(in_csysar)
       call exchange_add(in_csysven)
       call exchange_add(in_cpular)
       call exchange_add(in_cpulven)
       call exchange_add(in_lsysar)
       call exchange_add(in_rsysar)
       call exchange_add(in_lsysven)
       call exchange_add(in_rsysven)
       call exchange_add(in_lpular)
       call exchange_add(in_rpular)
       call exchange_add(in_lpulven)
       call exchange_add(in_rpulven)
       call exchange_add(in_elvpass)
       call exchange_add(in_ervpass)
       call exchange_add(in_elvmaxact)
       call exchange_add(in_ervmaxact)
       call exchange_add(in_lvtimact)
       call exchange_add(in_rvtimact)
       call exchange_add(in_startup_deltime)
       call exchange_add(in_startup_tpreload)
       call exchange_add(in_startup_deltpreload)
       call exchange_add(in_startup_fintime)
       call exchange_add(in_startup_ejection)
       call exchange_add(in_tcycle)
       call exchange_add(in_tisosys)
       call exchange_add(in_tsystole)
       call exchange_add(in_tdiastole)
       call exchange_add(in_tisodia)

       do i=1,4          
          call exchange_add(in_clo_factor(i))
          call exchange_add(in_ope_factor(i))
          call exchange_add(      in_rmin(i))
          call exchange_add(      in_rmax(i))
          do j=1,4          
             call exchange_add(in_tchamber(i,j))
          end do
       end do

       call exchange_add(in_nline_pp)
       call exchange_add(in_startup_cycles)


       call exchange_add(pressure_system_initial(1))
       call exchange_add(pressure_system_initial(2))
       call exchange_add(pressure_pulmonar_initial(1))
       call exchange_add(pressure_pulmonar_initial(2))

    end select

  end subroutine sysnet_sendat
  !
  ! Sysnet - Alya coupling through volumes and pressures
  !
  subroutine sysnet_update_phase(itask)
    use def_master,           only : kfl_paral,ITASK_INITIA,ITASK_ENDITE
    use def_master,           only : dtime
    implicit none
    ! -------------------------------
    integer(ip),  intent(in)       :: itask
    ! -------------------------------

    select case (itask)
       
    case(ITASK_INITIA)
       !
       ! Startup loops: run solo, using a synthetic volume loop with initial volumes computed from the real geometry
       !

       isolo         = "solo"          ! solo startup initialization run       
       iphase_sysnet(1:2)   = "sys"           ! always startup phase starts in sys phase

       iphase_change = "pressure"      ! phase change strategy by pressure differences
       ippqp         = "qp"            ! q-p method: fixing ventricles flow rates

       
    case (ITASK_ENDITE)
       !
       ! Coupled loops, with two phases:
       !    1. pre: preload, alya runs alone not coupled to sysnet
       !    2. idi,sys,isy,dia: normal sysnet-alya coupled cycles   

       isolo = "alya"                  ! coupled       

       if (chron_preload < in_startup_tpreload) then
          chron_preload= chron_preload + dtime
          ! only update the chrono and return
          iphase_sysnet(1:2) = "pre"           ! always alya coupled starts in preload phase
          dtime = in_startup_deltpreload
       else
          if (adjustl(trim(iphase_sysnet(2))) == 'pre') then
             !
             ! The first phase after pre is idi, dtime is back to the good one and the first dvoldt is set to zero
             !
             iphase_sysnet(1:2) = "idi"
             dtime = dtime_old_sysnet
!!             dvoldt= 0.0_rp

             if (kfl_paral == 0) write(6,*) 'dvoldt', dvoldt
             
          end if
       end if
       
    end select
  end subroutine sysnet_update_phase

  !
  ! Sysnet - Alya coupling through volumes and pressures
  !
  subroutine sysnet_alya_couple(itask)
    use def_master,           only :  ittim, ITER_K, TIME_N, TIME_N_MINUS_1, TIME_N_MINUS_2
    use def_master,           only :  dtime
    implicit none
    ! -------------------------------
    integer(ip),  intent(in)       :: itask
    ! -------------------------------
    integer(ip)                    :: icavi,jttim,kttim,mttim_dvoldt = 1, mttim_pressu = 1, max_average, mini
    real(rp)                       :: conversion_factor, minv, v_ascendente(2,10), v_desordenado(2,10)
!    logical(lg)                    :: mean_mask(2,10)
    
    select case (itask)

    case(1_ip)

       max_average= 4_ip
       
       ! update ventricles sysnet volumes coming from alya
       
       if (ittim == 1) then !initially V^n-1 = V^n = V^i 
          cavities_sysnet(1) % volume(TIME_N_MINUS_2) = cavities_sysnet(1) % volume (TIME_N)
          cavities_sysnet(1) % volume(TIME_N_MINUS_1) = cavities_sysnet(1) % volume (TIME_N)   
          cavities_sysnet(2) % volume(TIME_N_MINUS_2) = cavities_sysnet(2) % volume (TIME_N)
          cavities_sysnet(2) % volume(TIME_N_MINUS_1) = cavities_sysnet(2) % volume (TIME_N)   

       else

          ! bdf1: dv/dt = (vol^n - vol^n-1) / dt
!!$             dvoldt(1) = (cavities_sysnet(1) % volume(TIME_N)-cavities_sysnet(1) % volume(TIME_N_MINUS_1))/dtime
!!$             dvoldt(2) = (cavities_sysnet(2) % volume(TIME_N)-cavities_sysnet(2) % volume(TIME_N_MINUS_1))/dtime

          ! average, or mean, or nothing
          dvoldt     = 0.0_rp
          mttim_dvoldt = min(mttim_dvoldt,ittim)             
          if (mttim_dvoldt >= 2) then
             do jttim= mttim_dvoldt,2,-1
                dvoldt_keep(1,jttim)= dvoldt_keep(1,jttim-1)
                dvoldt_keep(2,jttim)= dvoldt_keep(2,jttim-1)
                dvoldt(1) = dvoldt(1) + dvoldt_keep(1,jttim) / real(mttim_dvoldt)
                dvoldt(2) = dvoldt(2) + dvoldt_keep(2,jttim) / real(mttim_dvoldt)
             end do
          end if
          dvoldt_keep(1,1)= (cavities_sysnet(1) % volume(TIME_N)-cavities_sysnet(1) % volume(TIME_N_MINUS_1))/dtime
          dvoldt_keep(2,1)= (cavities_sysnet(2) % volume(TIME_N)-cavities_sysnet(2) % volume(TIME_N_MINUS_1))/dtime
          dvoldt(1) = dvoldt(1) + dvoldt_keep(1,1) / real(mttim_dvoldt)
          dvoldt(1) = dvoldt(1) + dvoldt_keep(1,1) / real(mttim_dvoldt)
          
          ! median

          v_desordenado = 0.0_rp
          v_ascendente  = 0.0_rp
          v_desordenado(1,1:mttim_dvoldt)= dvoldt_keep(1,1:mttim_dvoldt)
          do kttim= 1,mttim_dvoldt
             minv=minval(v_desordenado(1,1:mttim_dvoldt),1)
             mini=minloc(v_desordenado(1,1:mttim_dvoldt),1)!,kind=ip)             
             v_ascendente(1,kttim)=minv
             v_desordenado(1,mini)=huge(minv)             
          end do
          
          
!!$          mean_mask  = .false.
!!$          mean_mask(1:mttim_dvoldt,1) = .true.
!!$          mean_mask(1:mttim_dvoldt,2) = .true.
!!$
!!$          do jttim = 1, mttim_dvoldt
!!$             v_ascendente(1,jttim) = minval(dvoldt_keep(1:mttim_dvoldt,1),mean_mask(1:mttim_dvoldt,1))
!!$             mean_mask(minloc(dvoldt_keep(1:mttim_dvoldt,1),mean_mask(1:mttim_dvoldt,1)),1) = .false.
!!$             v_ascendente(2,jttim) = minval(dvoldt_keep(1:mttim_dvoldt,2),mean_mask(1:mttim_dvoldt,2))
!!$             mean_mask(minloc(dvoldt_keep(1:mttim_dvoldt,2),mean_mask(1:mttim_dvoldt,1)),2) = .false.
!!$          end do

          if (mttim_dvoldt < max_average) mttim_dvoldt = mttim_dvoldt + 1_ip

          
          ! bdf2: dv/dt = 3/2 (vol^n - 4/3 vol^n-1 + 1/3 vol^n-2) / dt
!!$             dvoldt(1) = 1.5_rp * (&
!!$                  cavities_sysnet(1) % volume(TIME_N) &
!!$                  - 1.33333_rp * cavities_sysnet(1) % volume(TIME_N_MINUS_1)&
!!$                  + 0.33333_rp * cavities_sysnet(1) % volume(TIME_N_MINUS_2)&
!!$                  )/dtime
!!$             dvoldt(2) = 1.5_rp * (&
!!$                  cavities_sysnet(2) % volume(TIME_N) &
!!$                  - 1.33333_rp * cavities_sysnet(2) % volume(TIME_N_MINUS_1)&
!!$                  + 0.33333_rp * cavities_sysnet(2) % volume(TIME_N_MINUS_2)&
!!$                  )/dtime


          ! update volumes in unksys, to be used by sysnet
          unksys(1) =  cavities_sysnet(1) % volume(TIME_N)         
          unksys(2) =  cavities_sysnet(2) % volume(TIME_N)         

       end if

    case (2_ip)
       ! update ventricles alya pressures coming from sysnet, convert from mm/hg to CGS

       !         the conversion_factor is 1013000.0_rp / 760.0_rp  ! pressures from mmHg to barye 
       conversion_factor= 1333.22_rp  ! pressures from mmHg to barye 
       mttim_pressu = min(mttim_pressu,ittim+1)             

       max_average= 4_ip ! this is the maximum number of pressure time points for the average

       do icavi = 1, 2
          cavities_sysnet(icavi) % pres(TIME_N_MINUS_2) = cavities_sysnet(icavi) % pres(TIME_N_MINUS_1)  !V^n-2
          cavities_sysnet(icavi) % pres(TIME_N_MINUS_1) = cavities_sysnet(icavi) % pres(TIME_N)          !V^n-1

!!$             cavities_sysnet(icavi) % pres(TIME_N)         = funcar(icavi) * conversion_factor             
!!$             if (ittim == 0) then !initially V^n-1 = V^n = V^i 
!!$                cavities_sysnet(icavi) % pres(TIME_N_MINUS_2) = cavities_sysnet(icavi) % pres (TIME_N)
!!$                cavities_sysnet(icavi) % pres(TIME_N_MINUS_1) = cavities_sysnet(icavi) % pres (TIME_N) 
!!$             end if

          cavities_sysnet(icavi) % pres(TIME_N) = 0.0_rp
          if (mttim_pressu >= 2) then
             ! compute the pressure as the mean value of the last 
             do jttim= mttim_pressu,2,-1
                pressu_keep(icavi,jttim)= pressu_keep(icavi,jttim-1)
                cavities_sysnet(icavi) % pres(TIME_N) = &
                     cavities_sysnet(icavi) % pres(TIME_N) + pressu_keep(icavi,jttim) / real(mttim_pressu)
             end do
          end if
          pressu_keep(icavi,1)=  funcar(icavi) * conversion_factor             
          cavities_sysnet(icavi) % pres(TIME_N) = cavities_sysnet(icavi) % pres(TIME_N) &
               + pressu_keep(icavi,1) / real(mttim_pressu)

          if (ittim == 0) then !initially V^n-1 = V^n = V^i 
             cavities_sysnet(icavi) % pres(TIME_N_MINUS_2) = cavities_sysnet(icavi) % pres (TIME_N)
             cavities_sysnet(icavi) % pres(TIME_N_MINUS_1) = cavities_sysnet(icavi) % pres (TIME_N) 
          end if

       end do

       if (mttim_pressu < max_average) mttim_pressu = mttim_pressu + 1_ip


    end select

  end subroutine sysnet_alya_couple
  
  !
  ! Compute cavity volumes
  !
  subroutine sysnet_compute_volumes(itask)
  ! ---------------------------------------------------------------------
  !> @author   Adria Quintanas-Corominas (adria.quintanas@bsc.es)
  !> @author   Alfonso Santiago (alfonso.santiago@bsc.es)
  !> @details  Calculation of the volume of a cavity cutted by the vasal plane:
  !>           V = Int_{dOmega_0} (x \cdot e_1) (e_1 \cdot n) d dOmega_0    
  !>           where, x : current coordinates
  !>                  e : vector in the basal plane
  !>                  n : vector normal to the cavity surface
  !>           Refs: Levrero-Florencio et al (2020) DOI: 10.1016/j.cma.2019.112762
  !>                 Quarteroni et al. (2017)       DOI: 10.1016/j.cma.2016.05.031 
  ! ---------------------------------------------------------------------
    use def_domain,           only :  mnode, mnodb, ltypb, nnode, ngaus, lelbo, ltype, nboun, npoin_own
    use def_domain,           only :  elmar, ndimb, lbset, lnodb, coord, lnods
    use def_domain,           only :  htable_lninv_loc
    use def_master,           only :  ittim, displ, intost
    use def_master,           only :  INOTMASTER, ITER_K, TIME_N, TIME_N_MINUS_1, TIME_N_MINUS_2, ITASK_ENDITE
    use mod_htable,           only :  htalid
    use mod_communications,   only :  PAR_SUM
    use mod_maths_basic,      only :  maths_cross_product, maths_normalize_vector
    use mod_bouder
    use mod_messages,         only :  messages_live

    !-----------------------------
    implicit none
    !-----------------------------
    ! -------------------------------
    integer(ip),  intent(in)       :: itask
    ! -------------------------------
    integer(ip)                    :: icavi, ielem, ipoin, idime, iboun, igaub, inodb, ii
    integer(ip)                    :: pelty, pblty, pnodb, pgaub, pnode, inode
    real(rp)                       :: gpel(3), v1(3), v2(3)
    real(rp)                       :: baloc(3,3), bocod(3,mnodb), elcod(3,mnode), eucta
    logical(lg)                    :: add_vol_contribution

    ! -------------------------------
    
    !
    ! Initialise  
    !
    do icavi= 1, max_cavities
       cavities_sysnet(icavi) % volume(ITER_K) = 0.0_rp
       cavities_sysnet(icavi) % plane_points(:,:) = 0.0_rp
    enddo

    !
    ! Locate the nodes of the plane in the domain
    !
    if( INOTMASTER ) then

      do icavi= 1, max_cavities
         if( cavities_sysnet(icavi) % initialized )then
           do ii = 1, size( cavities_sysnet(icavi) % plane_nodset, 1_ip, kind=ip )
              if ( cavities_sysnet(icavi) % plane_nodset(ii) > 0_ip ) then
                 ipoin = htalid(htable_lninv_loc,cavities_sysnet(icavi) % plane_nodset(ii))
                 if( ipoin > 0_ip .and. ipoin <= npoin_own ) then
                     cavities_sysnet(icavi) % plane_points(:,ii) = coord(:,ipoin) + displ(:,ipoin,1)
                 end if
              end if
           enddo
         end if
      end do


    end if

    ! mod_parall global_to_local
    do icavi = 1, max_cavities
       do ii = 1, size( cavities_sysnet(icavi) % plane_points, 2_ip, kind=ip )
          call PAR_SUM( 3_ip, cavities_sysnet(icavi) % plane_points(:,ii) )
       enddo
    enddo

    !
    ! Compute the new origin and normal vector of the plane
    !
    do icavi= 1, max_cavities
      if( cavities_sysnet(icavi) % initialized .and. all(cavities_sysnet(icavi) % plane_nodset(1:2) > 0_ip)  )then
           v1(:) = cavities_sysnet(icavi) % plane_points(:,1) - cavities_sysnet(icavi) % plane_points(:,2) 
           call maths_normalize_vector(3_ip, v1)
           cavities_sysnet(icavi) % plane_vector(:) = v1(:)

           if ( cavities_sysnet(icavi) % plane_nodset(3) > 0_ip ) then ! 3 nodes are provided, caculate plane normal
              ! Cavity view from top at the hole. Points clockwise on the boundary plane
              ! 
              !         v1
              !  p1 <--------- p2
              !   \           /
              !    \         / v2
              !     \       /
              !      \    |/_
              !         p3
              !
              ! v1 = p1-p2
              ! v2 = p3-p2
              ! normal = v2 x v1, pointing inward.
              v2(:) =  cavities_sysnet(icavi) % plane_points(:,3) - cavities_sysnet(icavi) % plane_points(:,2)
              call maths_normalize_vector(3_ip, v2)
              cavities_sysnet(icavi) % plane_normal(:) = maths_cross_product( v2(:), v1(:), 3_ip )
              cavities_sysnet(icavi) % plane_origin(:) = cavities_sysnet(icavi) % plane_points(:,1)
           end if
      endif
   enddo

    !
    ! Compute the volume of the cavity
    !
    if( INOTMASTER ) then

       do iboun = 1, nboun
          do icavi= 1, max_cavities
             if( cavities_sysnet(icavi) % initialized )then
                if( lbset(iboun) == cavities_sysnet(icavi) % bouset )then

                   pblty=ltypb(iboun) 
                   pnodb=nnode(pblty)
                   pgaub=ngaus(pblty)
                   ielem=lelbo(iboun)
                   pelty=ltype(ielem)
                   pnode=nnode(pelty)

                   ! Gather bondary and element coordiantes
                   do inodb=1,pnodb
                      ipoin=lnodb(inodb,iboun)
                      do idime=1, 3
                         bocod(idime,inodb) = coord(idime,ipoin) + displ(idime,ipoin,1)
                      end do
                   end do

                   do inode=1,pnode
                      ipoin=lnods(inode,ielem)
                      do idime=1, 3
                         elcod(idime,inode) = coord(idime,ipoin) + displ(idime,ipoin,1)
                      end do
                   end do

                   ! Divergence theorem
                   do igaub = 1, pgaub

                      call bouder(pnodb,3_ip,ndimb,elmar(pblty)%deriv(:,:,igaub),bocod,baloc,eucta)    
                      call chenor(pnode,baloc,bocod,elcod)

                      gpel = 0.0_rp
                      do inodb = 1_ip,pnodb
                         do idime= 1_ip, 3_ip
                            gpel(idime) = gpel(idime) + elmar(pblty)%shape(inodb,igaub) * &
                                 bocod(idime,inodb)
                         end do
                      end do


                      add_vol_contribution = .TRUE.

                      ! if the usr specified a 3 point plane, exclude the faces in the exterior of the cavity
                      if ( cavities_sysnet(icavi) % plane_nodset(3) > 0_ip ) then ! 3 nodes are provided, check on which side of the plane is the point
                        if ( dot_product( gpel(:) - cavities_sysnet(icavi) % plane_origin(:), cavities_sysnet(icavi) % plane_normal(:) ) < 0 ) then
                           add_vol_contribution = .FALSE.
                        end if
                      end if

                      ! includes shifting the gpel to the coordinate system centered at one of the vector nodes
                      if (add_vol_contribution) then
                           cavities_sysnet(icavi) % volume(ITER_K) =  cavities_sysnet(icavi) % volume(ITER_K) &
                                    &    - dot_product(gpel(:) - cavities_sysnet(icavi) % plane_points(:,1), cavities_sysnet(icavi) % plane_vector(:)) &
                                    &    * dot_product(cavities_sysnet(icavi) % plane_vector(:), baloc(:,3)) &
                                    &    * elmar(pblty) % weigp(igaub) * eucta
                      end if
                   enddo
                end if
             endif
          end do
       end do

    endif

    do icavi = 1, max_cavities
       call PAR_SUM( cavities_sysnet(icavi) % volume(ITER_K) )
    enddo

    !
    ! From solidz to master arrays
    !
    do icavi = 1, max_cavities 

       !update volume, only at the end of the iterations loop
       if (itask == ITASK_ENDITE) then
          cavities_sysnet(icavi) % volume(TIME_N_MINUS_2) = cavities_sysnet(icavi) % volume(TIME_N_MINUS_1)  !V^n-2
          cavities_sysnet(icavi) % volume(TIME_N_MINUS_1) = cavities_sysnet(icavi) % volume(TIME_N)          !V^n-1
!          cavities_sysnet(icavi) % volume(TIME_N)         = cavities_sysnet(icavi) % volume(ITER_K)          !V^n

          cavities_sysnet(icavi) % volume(TIME_N)         = 0.33333_rp * (&
               cavities_sysnet(icavi) % volume(ITER_K)      &    !V^n
               + cavities_sysnet(icavi) % volume(TIME_N)    &      !V^n-1
               + cavities_sysnet(icavi) % volume(TIME_N_MINUS_1) )  !V^n-2

          if (ittim == 1) then
             cavities_sysnet(icavi) % volume(TIME_N_MINUS_2) = cavities_sysnet(icavi) % volume (TIME_N)
             cavities_sysnet(icavi) % volume(TIME_N_MINUS_1) = cavities_sysnet(icavi) % volume (TIME_N)      !initially V^n-1 = V^n = V^i 
          end if
       end if

    end do
    
  end subroutine sysnet_compute_volumes
  !
  ! Initialize standalone case
  !
  subroutine sysnet_initialize_parameters(ippqp)

    implicit none
    character(30) :: ippqp

    ! model input parameters: in_*

    ! synthetic ventricular model
    ! left ventricle
    in_maxvolven_synth(1)= 120.0_rp
    in_minvolven_synth(1)=  40.0_rp
    ! right ventricle
    in_maxvolven_synth(2)= 110.0_rp
    in_minvolven_synth(2)=  30.0_rp
    
    ! baseline parameters:
    in_rmin= 0.0075_rp      ! mmHg s / mL -> valor mio, si pones el de quarte, se jode
    in_rmax= 75006.2_rp     ! mmHg s / mL de quarte
    in_clo_factor= 1.0_rp   ! no units
    in_ope_factor= 1.0_rp   ! no units
    
    ! atrial 0D parameters:
    ! baseline volumes
    in_vola= 28.0_rp          ! mL
    in_vora= 28.0_rp          ! mL
    ! passive elastances
    in_elapass= 0.09_rp       ! mmHg / mL - quarte 
    in_erapass= 0.07_rp       ! mmHg / mL - quarte
    ! active maximum elastances
    in_elamaxact= 0.07_rp     ! mmHg / mL
    in_eramaxact= 0.06_rp     ! mmHg / mL
    ! times:
    if (adjustl(trim(ippqp)) == 'qp') then
       in_tcycle      = 0.75_rp  ! s    
       in_tsystole    =  0.1875_rp   * in_tcycle
       in_tisosys     =  0.0625_rp   * in_tcycle
       in_tdiastole   =  0.6875_rp   * in_tcycle
       in_tisodia     =  0.0625_rp   * in_tcycle
    else if (adjustl(trim(ippqp)) == 'pp') then
       ! when pp, these timings are set from the pressure cycle
       in_tcycle      = 0.75_rp  ! s    
       in_tsystole    =  0.1875_rp   * in_tcycle
       in_tisosys     =  0.1200_rp   * in_tcycle
       in_tdiastole   =  0.6300_rp   * in_tcycle
       in_tisodia     =  0.0625_rp   * in_tcycle
    end if       
    in_tchamber   = 0.0_rp
    ! from Blanco 2010, which uses a in_tcycle=1.00
    ! left atrium
    in_tchamber(1,1) = 0.17_rp * in_tcycle      ! Tac
    in_tchamber(2,1) = 0.17_rp * in_tcycle      ! Tar
    in_tchamber(3,1) = 0.80_rp * in_tcycle      ! tac
    in_tchamber(4,1) = 0.97_rp * in_tcycle      ! tar
    ! right atrium
    in_tchamber(1,2) = 0.17_rp * in_tcycle      ! Tac
    in_tchamber(2,2) = 0.17_rp * in_tcycle      ! Tar
    in_tchamber(3,2) = 0.80_rp * in_tcycle      ! tac
    in_tchamber(4,2) = 0.97_rp * in_tcycle      ! tar


    ! dont need them, but just to complete the description: these are ventricular 0D parameters
    !    in_vorv= 10.0_rp          ! mL
    !    in_volv= 10.0_rp          ! mL
    !    in_elvpass= 0.08_rp       ! mmHg / mL
    !    in_elvmaxact= 2.75_rp     ! mmHg / mL
    !    in_ervpass= 0.05_rp       ! mmHg / mL
    !    in_ervmaxact= 0.55_rp     ! mmHg / mL


!!    in_tcycle      = 0.800_rp  ! s
!!    in_tsystole    =  0.15_rp
!!    in_tdiastole   =  0.55_rp
!!    in_tisosys     =  0.05_rp
!!    in_tisodia     =  0.05_rp

    ! valores de quarte
    ! Nota: tuve que bajar las resistencias porque sino perdia demasiada carga por el camino
    ! y salian presiones negativas al llegar a la auricula izquierda
    in_csysar = 1.200_rp      ! mL / mmHg 
    in_lsysar = 5.0e-3_rp        ! mmHg s^2 / mL  
!    in_rsysar = 0.8000_rp    ! mmHg s / mL    
!    in_rsysar = 0.6000_rp     ! mmHg s / mL  
    in_rsysar = 0.5000_rp     ! mmHg s / mL  

    in_csysven= 6.000_rp     ! mL / mmHg   
    in_lsysven= 5.0e-4_rp        ! mmHg s^2 / mL    
    in_rsysven= 0.2600_rp     ! mmHg s / mL    

    in_cpular = 10.000_rp     ! mL / mmHg  
    in_lpular = 5.0e-4_rp        ! mmHg s^2 / mL   
!    in_rpular = 0.1625_rp     ! mmHg s / mL
!    in_rpular = 0.1225_rp     ! mmHg s / mL
    in_rpular = 0.1_rp     ! mmHg s / mL

    in_cpulven= 16.000_rp     ! mL / mmHg    
    in_lpulven= 5.0e-4_rp        ! mmHg s^2 / mL   
!    in_rpulven= 0.1625_rp     ! mmHg s / mL    
!    in_rpulven= 0.1225_rp     ! mmHg s / mL    
    in_rpulven= 0.1_rp     ! mmHg s / mL    
    
    ! valores de hirschvogel
!    in_csysar = 1.2000_rp     ! mL / mmHg 
!    in_cpular = 10.000_rp     ! mL / mmHg  
!    in_rpular = 0.1625_rp     ! mmHg s / mL
!    in_rsysar = 0.8000_rp     ! mmHg s / mL  
!    in_csysven= 30.0_rp * in_csysar     ! mL / mmHg   
!    in_cpulven= 2.5_rp * in_cpular     ! mL / mmHg    
!    in_rpulven= in_rpular     ! mmHg s / mL    
!    in_rsysven= in_rsysar / 5.0_rp     ! mmHg s / mL    
    
    
    ! initial guess
    funperi = 0.0_rp

    
  end subroutine sysnet_initialize_parameters

  !
  ! Interface to sysnet iteration step
  !
  subroutine sysnet_time_iteration_step(dts)

    implicit none
    real(rp) :: dts

    !
    ! RK4:
    !
    ! aux= unk_n  
    unkaux = unksys
    ! k_1 = rhs(aux,t)
    call compute_rhs(1_ip,unkaux)

    !
    ! aux= unk + rhs(1:12,1) * dt / 2.0_rp
    unkaux = unksys + rhs(1:12,1) * dts / 2.0_rp
    ! k_2 = rhs(aux, t+dt/2)
    call compute_rhs(2_ip,unkaux)

    !
    ! aux = unk + k_2 * dt / 2
    unkaux = unksys + rhs(1:12,2) * dts / 2.0_rp
    ! k_3 = rhs(aux, t+dt/2)
    call compute_rhs(3_ip,unkaux)

    !
    ! aux = unk + k_3 * dt 
    unkaux = unksys + rhs(1:12,3) * dts
    ! k_4 = rhs(aux, t+dt  )
    call compute_rhs(4_ip,unkaux)

    !
    ! unk = unk + ( k_1 / 6 + k_2 / 3 + k_3 / 3 + k_4 / 6 ) * dt
    call update_systemic_unknowns(dts)

  end subroutine sysnet_time_iteration_step
  
  subroutine sysnet_update_functions (dts,debuni)

    real(rp) :: dts
    integer(ip)   :: debuni

    call update_cardiac_functions(elastan,unksys,funcar,funperi,debuni)
    call update_cardiac_elastances(elastan)

!  antes eran asi, y cts era chroncycle    
!    call update_cardiac_functions(cts,elastan,unksys,funcar,resvalve,funperi,debuni)
!    call update_cardiac_elastances(cts,elastan)
    
    call update_derived_functions(unksys)


  end subroutine sysnet_update_functions
  
  subroutine sysnet_update_cycle (dts,debuni,coun_cycle)

    implicit none
    real(rp) :: dts, xlapse, xpi, xsigmoid
    integer(ip) :: icavity,ires_valve,debuni,coun_cycle(2)

    if (adjustl(trim(iphase_change)) == 'duration') then
       !
       ! This is a cardiac cycle manager based on times
       !
       if (chron_cycle(1) < in_tsystole) then
          iphase_sysnet(1:2)= "sys"

          resvalve(1) = in_rmax(1)  ! mitral closed
          resvalve(2) = in_rmin(2)  ! aortic open
          resvalve(3) = in_rmax(3)  ! tricuspid closed
          resvalve(4) = in_rmin(4)  ! pulmonic open

          chron_phase(1,1)= chron_phase(1,1) + dts
          chron_phase(1,2)= chron_phase(1,2) + dts
          chron_phase(2,1:2)= 0.0_rp
          chron_phase(3,1:2)= 0.0_rp
          chron_phase(4,1:2)= 0.0_rp       

       else if (chron_cycle(1) < (in_tisosys + in_tsystole)) then       
          iphase_sysnet(1:2)= "isy" ! isovol after systole, isovolumic contraction

          resvalve(1) = in_rmax(1)  ! mitral closed
          resvalve(2) = in_rmax(2)  ! aortic closed
          resvalve(3) = in_rmax(3)  ! tricuspid closed
          resvalve(4) = in_rmax(4)  ! pulmonic closed

          chron_phase(1,1:2)= 0.0_rp
          chron_phase(2,1)= chron_phase(2,1) + dts
          chron_phase(2,2)= chron_phase(2,2) + dts
          chron_phase(3,1:2)= 0.0_rp
          chron_phase(4,1:2)= 0.0_rp       
       else if (chron_cycle(1) < (in_tisosys + in_tsystole + in_tdiastole)) then        
          iphase_sysnet(1:2)= "dia"        

          resvalve(1) = in_rmin(1)  ! mitral open
          resvalve(2) = in_rmax(2)  ! aortic closed
          resvalve(3) = in_rmin(3)  ! tricuspid open
          resvalve(4) = in_rmax(4)  ! pulmonic closed

          chron_phase(1,1:2)= 0.0_rp
          chron_phase(2,1:2)= 0.0_rp
          chron_phase(3,1)= chron_phase(3,1) + dts
          chron_phase(3,2)= chron_phase(3,2) + dts
          chron_phase(4,1:2)= 0.0_rp       
       else if (chron_cycle(1) < (in_tisosys + in_tsystole + in_tdiastole + in_tisodia)) then        
          iphase_sysnet(1:2)= "idi" ! isovol after diastole, isovolumic expansion
          resvalve(1) = in_rmax(1)  ! mitral closed
          resvalve(2) = in_rmax(2)  ! aortic closed
          resvalve(3) = in_rmax(3)  ! tricuspid closed
          resvalve(4) = in_rmax(4)  ! pulmonic closed

          chron_phase(1,1:2)= 0.0_rp       
          chron_phase(2,1:2)= 0.0_rp
          chron_phase(3,1:2)= 0.0_rp
          chron_phase(4,1)= chron_phase(4,1) + dts
          chron_phase(4,2)= chron_phase(4,2) + dts
       end if

       !!       write(6,*) funcar(1),unksys(5), iphase_sysnet(1)

       ! Reset cycle time chrono
       if (chron_cycle(1)  > in_tcycle) chron_cycle(1) = 0.0_rp 
       if (chron_cycle(2)  > in_tcycle) chron_cycle(2) = 0.0_rp 

    else if (adjustl(trim(iphase_change)) == 'pressure') then
       !
       ! This is a cardiac cycle manager based on pressures on both sides of a valve
       !

       if (adjustl(trim(isolo)) == 'solo') then

          resvalve(1) = in_rmax(1)  ! mitral closed
          resvalve(2) = in_rmax(2)  ! aortic closed
          resvalve(3) = in_rmax(3)  ! tricuspid closed
          resvalve(4) = in_rmax(4)  ! pulmonic closed       

          if (funcar(3) > funcar(1))  resvalve(1)=in_rmin(1) ! open mitral    if p_la > p_lv
          if (funcar(1) >= unksys(5)) resvalve(2)=in_rmin(2) ! open aortic    if p_lv > p_sysar
          if (funcar(4) > funcar(2))  resvalve(3)=in_rmin(3) ! open tricuspid if p_ra > p_rv
          if (funcar(2) >= unksys(7)) resvalve(4)=in_rmin(4) ! open pulmonic  if p_rv > p_pular

       else 

          resvalve(1) = in_rmin(1)  ! mitral open
          resvalve(2) = in_rmin(2)  ! aortic open
          resvalve(3) = in_rmin(3)  ! tricuspid open
          resvalve(4) = in_rmin(4)  ! pulmonic open
          
          xpi= 3.1416_rp
          if (adjustl(trim(ivalve_movement)) == 'smooth') then
             ! open mitral    if p_la > p_lv
             write(6,*) 
             if (funcar(1) > funcar(3)) then
                xlapse= 20.0_rp   ! close time will be xlapse
                chron_valve_move(1,1) = chron_valve_move(1,1) + dts
                ! xsigmoid must be in miliseconds and starting at -xlapse/2
                xsigmoid = chron_valve_move(1,1) * 1000.0_rp - xlapse / 2.0_rp
                if (xsigmoid < xlapse/2.0_rp) then
                   resvalve(1)=(in_rmax(1) + in_rmin(1)) * 0.5_rp * &
                        ( &
                        1.0_rp &
                        + 2.0_rp * xsigmoid / xlapse &
                        + 1.0_rp/xpi * sin( xpi * 2.0_rp * xsigmoid / xlapse) &
                        ) + in_rmin(1)

                else                   
                   resvalve(1)= in_rmax(1) 
                end if

             end if

             if (unksys(5) >= funcar(1)) then ! close aortic    if p_lv > p_sysar
                resvalve(2)= in_rmax(2)
             end if

             ! open tricuspic if p_ra > p_rv
             if (funcar(2) > funcar(4)) then 
                xlapse= 20.0_rp   ! closing time will be xlapse
                chron_valve_move(3,1) = chron_valve_move(3,1) + dts
                ! xsigmoid must be in miliseconds and starting at -xlapse/2
                xsigmoid = chron_valve_move(1,1) * 1000.0_rp - xlapse / 2.0_rp
                if (xsigmoid < xlapse/2.0_rp) then
                   resvalve(3)=(in_rmax(3) + in_rmin(3)) * 0.5_rp * &
                        ( &
                        1.0_rp &
                        + 2.0_rp * xsigmoid / xlapse &
                        + 1.0_rp/xpi * sin( xpi * 2.0_rp * xsigmoid / xlapse) &
                        ) + in_rmin(3)
                else
                   resvalve(3)= in_rmax(3) 
                end if
             end if

             if (unksys(7) >= funcar(2)) then
                resvalve(4)=in_rmax(4) ! close pulmonic  if p_rv > p_pular
             end if

                write(6,*) 
                write(6,*) 'pipiiii valves updated',resvalve(1:4)

                
          else

             if (funcar(3) > funcar(1))  resvalve(1)=in_rmin(1) ! open mitral    if p_la > p_lv
             if (funcar(1) >= unksys(5)) resvalve(2)=in_rmin(2) ! open aortic    if p_lv > p_sysar
             if (funcar(4) > funcar(2))  resvalve(3)=in_rmin(3) ! open tricuspic if p_ra > p_rv
             if (funcar(2) >= unksys(7)) resvalve(4)=in_rmin(4) ! open pulmonic  if p_rv > p_pular

          end if

       end if

       if (adjustl(trim(isolo)) == 'solo') then
          !
          ! This is done only for the startup solo phase
          !
          if (adjustl(trim(iphase_sysnet(1)))== "dia") then
             if (chron_phase(3,1) > in_tdiastole) then
                iphase_sysnet(1)='idi'
                coun_cycle(1)= coun_cycle(1) + 1_ip
             end if
          end if
          if (adjustl(trim(iphase_sysnet(2)))== "dia") then
             if (chron_phase(3,2) > in_tdiastole) then
                iphase_sysnet(2)='idi'
                coun_cycle(2)= coun_cycle(2) + 1_ip                
             end if
          end if

          if (adjustl(trim(iphase_sysnet(1)))== "idi") then
             if (resvalve(1) == in_rmin(1)) resvalve(1) = in_rmax(1) !force closing av valves if open
          end if
          if (adjustl(trim(iphase_sysnet(2)))== "idi") then
             if (resvalve(3) == in_rmin(3)) resvalve(3) = in_rmax(3) !force closing av valves if open
          end if

       end if

       do icavity= 1,2
          ires_valve= (icavity-1)*2_ip

          if (adjustl(trim(iphase_sysnet(icavity))) == 'sys') then             
             if (resvalve(1+ires_valve) > in_rmin(1+ires_valve)) then       ! mitral / triscupid closed
                if(resvalve(2+ires_valve) < in_rmax(2+ires_valve)) then     ! aortic / pulmonic open
                   iphase_sysnet(icavity)= "sys" 
                   chron_phase(1,icavity)= chron_phase(1,icavity) + dts
                   chron_phase(2,icavity)= 0.0_rp
                   chron_phase(3,icavity)= 0.0_rp
                   chron_phase(4,icavity)= 0.0_rp             
                else                              ! aortic / pulmonic closed
                   iphase_sysnet(icavity)= "isy"              
                   chron_phase(1,icavity)= 0.0_rp
                   chron_phase(2,icavity)= chron_phase(2,icavity) + dts
                   chron_phase(3,icavity)= 0.0_rp
                   chron_phase(4,icavity)= 0.0_rp                    
                   chron_valve_move(1+ires_valve,1:2) = 0.0_rp   ! valves ready to be opened or closed
                   chron_valve_move(2+ires_valve,1:2) = 0.0_rp   ! valves ready to be opened or closed
                end if
             end if

          else if (adjustl(trim(iphase_sysnet(icavity))) == 'isy') then

             !                write(6,*) ires_valve, icavity, iphase_sysnet(1:2)
             !                write(6,*) funcar(1), unksys(5), resvalve(2)
             !                write(6,*) funcar(2), unksys(7), resvalve(4) ! open pulmonic  if p_rv > p_pular
             !                write(6,*) funcar(3), funcar(1), resvalve(1) ! open mitral    if p_la > p_lv
             !                write(6,*) funcar(4), funcar(2), resvalve(3) ! open tricuspic if p_ra > p_rv


             if (resvalve(2+ires_valve) > in_rmin(2+ires_valve)) then       ! aortic / pulmonic closed
                if(resvalve(1+ires_valve) > in_rmin(1+ires_valve)) then     ! mitral / triscupid closed
                   iphase_sysnet(icavity)= "isy"              
                   chron_phase(1,icavity)= 0.0_rp
                   chron_phase(2,icavity)= chron_phase(2,icavity) + dts
                   chron_phase(3,icavity)= 0.0_rp
                   chron_phase(4,icavity)= 0.0_rp                    
                else                              ! mitral / triscupid open
                   iphase_sysnet(icavity)= "dia" 
                   chron_phase(1,icavity)= 0.0_rp
                   chron_phase(2,icavity)= 0.0_rp
                   chron_phase(3,icavity)= chron_phase(3,icavity) + dts
                   chron_phase(4,icavity)= 0.0_rp       
                end if
             end if

          else if (adjustl(trim(iphase_sysnet(icavity))) == 'dia') then

             if (resvalve(2+ires_valve) > in_rmin(2+ires_valve)) then       ! aortic / pulmonic closed
                if(resvalve(1+ires_valve) < in_rmax(1+ires_valve)) then     ! mitral / triscupid open
                   iphase_sysnet(icavity)= "dia" 
                   chron_phase(1,icavity)= 0.0_rp
                   chron_phase(2,icavity)= 0.0_rp
                   chron_phase(3,icavity)= chron_phase(3,icavity) + dts
                   chron_phase(4,icavity)= 0.0_rp                
                else                              ! mitral / triscupid closed                
                   iphase_sysnet(icavity)= "idi"              
                   chron_phase(1,icavity)= 0.0_rp       
                   chron_phase(2,icavity)= 0.0_rp
                   chron_phase(3,icavity)= 0.0_rp
                   chron_phase(4,icavity)= chron_phase(4,icavity) + dts
                   chron_valve_move(1+ires_valve,1:2) = 0.0_rp   ! valves ready to be opened or closed
                   chron_valve_move(2+ires_valve,1:2) = 0.0_rp   ! valves ready to be opened or closed
                end if
             end if

          else if (adjustl(trim(iphase_sysnet(icavity))) == 'idi') then          

             if (resvalve(1+ires_valve) > in_rmin(1+ires_valve)) then       ! mitral / triscupid closed
                if(resvalve(2+ires_valve) > in_rmin(2+ires_valve)) then     ! aortic / pulmonic closed
                   iphase_sysnet(icavity)= "idi"              
                   chron_phase(1,icavity)= 0.0_rp       
                   chron_phase(2,icavity)= 0.0_rp
                   chron_phase(3,icavity)= 0.0_rp
                   chron_phase(4,icavity)= chron_phase(4,icavity) + dts
                else                              ! aortic / pulmonic open
                   iphase_sysnet(icavity)= "sys" 
                   chron_phase(1,icavity)= chron_phase(1,icavity) + dts
                   chron_phase(2,icavity)= 0.0_rp
                   chron_phase(3,icavity)= 0.0_rp
                   chron_phase(4,icavity)= 0.0_rp             
                end if
             end if
          end if

       end do !!! icavity


       ! Reset cycle time chrono 
       if (chron_cycle(1)  > in_tcycle) then
          chron_cycle(1) = 0.0_rp
       end if
       if (chron_cycle(2)  > in_tcycle) then
          chron_cycle(2) = 0.0_rp 
       end if
    end if


  end subroutine sysnet_update_cycle

  !
  ! Initialize some parameters of the run
  !
  subroutine sysnet_initialize_run(isolo_in) 

    implicit none
    character(4)  :: isolo_in
    
    ! initialize sysnet
    chron_phase   = 0.0_rp        ! phase chronometer
    chron_cycle   = 0.0_rp        ! cycle chronometer
    chron_total   = 0.0_rp        ! total chronometer
    isolo  = adjustl(trim(isolo_in))        ! "solo" or "alya" initialization run
    iphase_sysnet(1:2) = "initial"            ! initialization phase 
    
    iphase_change = "pressure"    ! phase change strategy by pressure differences
    
    ippqp  = "qp"                 ! q-p method: fixing ventricles flow rates

  end subroutine sysnet_initialize_run
  
  !
  ! Initialize variables and parameters
  !
  subroutine sysnet_initialize_guess

    implicit none
    
    ! unksys (systemic unknowns):
    !  1: V_lv, (2) <----- from left ventricle 
    !  2: V_rv, (4) <----- from right ventricle
    !  3: V_la, (1)
    !  4: V_ra, (3)
    !  5: P^system_arterial, 
    !  6: P^system_venous, 
    !  7: P^pulmonary_arterial, 
    !  8: P^pulmonary_venous, 
    !  9: Q^system_arterial, 
    ! 10: Q^system_ven, 
    ! 11: Q^pulmonar_arterial, 
    ! 12: Q^pulmonar_venous, 

    ! funcar (cardiac functions):
    !  1: P_lv, (1)
    !  2: P_rv, (3)
    !  3: P_la, (2)
    !  4: P_ra, (4)
    !  5: Q_mv, (mitral)
    !  6: Q_av, (aortic)
    !  7: Q_tv, (tricuspid)
    !  8: Q_pv, (pulmonic)

    ! atrial initial volumes
    unksys(3)= in_vola + 2.0_rp
    unksys(4)= in_vora + 2.0_rp

    !  5: P^system_arterial, 
    !  6: P^system_venous, 
    !  7: P^pulmonary_arterial, 
    !  8: P^pulmonary_venous, 

!!    THESE ALL ARE GOOD VALUES 
!!    pressure_system_initial(1)   = 80.0_rp 
!!    pressure_system_initial(2)   = 10.0_rp
!!    pressure_pulmonar_initial(1) = 20.0_rp
!!    pressure_pulmonar_initial(2) = 05.0_rp

!!    pressure_system_initial(1)   = 150.0_rp 
!!    pressure_system_initial(2)   = 30.0_rp
!!    pressure_pulmonar_initial(1) = 30.0_rp
!!    pressure_pulmonar_initial(2) = 10.0_rp  
    
!!! initial pressures "pre-inflates" the full system
    unksys(5)= pressure_system_initial(1)   ! this looks like a good initialization...
    unksys(6)= pressure_system_initial(2) 
    unksys(7)= pressure_pulmonar_initial(1) 
    unksys(8)= pressure_pulmonar_initial(2) 

    ! ventricular and atrial initial pressures
    funcar(1)= unksys(5) 
    funcar(2)= unksys(7) 
    funcar(3)= unksys(6)
    funcar(4)= unksys(8)
    
  end subroutine sysnet_initialize_guess

  subroutine sysnet_update_ventricular_synthetic_volumes(dt_in)
    real(rp) :: v_new(2)
    real(rp) :: dt_in, cvol, avol, bvol

    if (adjustl(trim(iphase_sysnet(1))) == 'initial') then
       ! for the initial, both iphase_sysnet are the same, then, check only the first one 
       unksys(1) = in_maxvolven_synth(1) !<----- left ventricle, cycle starting point 
       unksys(2) = in_maxvolven_synth(2) !<----- right ventricle, cycle starting point 
       return
    end if

    ! isolo = solo means that sysnet uses an imposed volume change function 
    ! otherwise, volume changed is computed by the alya cardias model 

    ! ippqp = pp o qp method

    if (adjustl(trim(isolo)) == 'solo') then
       if (adjustl(trim(ippqp)) == 'qp') then

          !
          ! left ventricle system
          !
          if (adjustl(trim(iphase_sysnet(1)))== "sys") then
             ! pressure will be computed from volume variations in update_cardiac_functions
             ! compute volume from an imposed plot
             cvol= in_maxvolven_synth(1)
             avol= (in_maxvolven_synth(1) - in_minvolven_synth(1)) / in_tsystole / in_tsystole
             bvol= - avol * 2.0_rp * in_tsystole

             v_new(1) = avol * chron_phase(1,1) * chron_phase(1,1) + bvol * chron_phase(1,1) + cvol

          else if (adjustl(trim(iphase_sysnet(1)))== "isy") then ! isovol after systole
             ! pressure will be computed from an imposed plot in update_cardiac_functions
             ! keep v_new as it was
             v_new(1) = unksys(1)
          else if (adjustl(trim(iphase_sysnet(1)))== "dia") then
             ! pressure will be computed from volume variations in update_cardiac_functions
             ! compute volume from an imposed plot
             v_new(1) = in_minvolven_synth(1) &
                  + (in_maxvolven_synth(1) &
                  - in_minvolven_synth(1)) * chron_phase(3,1) / in_tdiastole
          else if (adjustl(trim(iphase_sysnet(1)))== "idi") then ! isovol after diastole
             ! pressure will be computed from an imposed plot in update_cardiac_functions
             ! keep v_new as it was
             v_new(1) = unksys(1)
          end if

          !
          ! right ventricle system
          !
          if (adjustl(trim(iphase_sysnet(2)))== "sys") then
             ! pressure will be computed from volume variations in update_cardiac_functions
             ! compute volume from an imposed plot
             cvol= in_maxvolven_synth(2)
             avol= (in_maxvolven_synth(2) - in_minvolven_synth(2)) / in_tsystole / in_tsystole
             bvol= - avol * 2.0_rp * in_tsystole

             v_new(2) = avol * chron_phase(1,2) * chron_phase(1,2) + bvol * chron_phase(1,2) + cvol

          else if (adjustl(trim(iphase_sysnet(2)))== "isy") then ! isovol after systole
             ! pressure will be computed from an imposed plot in update_cardiac_functions
             ! keep v_new as it was
             v_new(2) = unksys(2)
          else if (adjustl(trim(iphase_sysnet(2)))== "dia") then
             ! pressure will be computed from volume variations in update_cardiac_functions
             ! compute volume from an imposed plot
             v_new(2) = in_minvolven_synth(2) &
                  + (in_maxvolven_synth(2) &
                  - in_minvolven_synth(2)) * chron_phase(3,2) / in_tdiastole          
          else if (adjustl(trim(iphase_sysnet(2)))== "idi") then ! isovol after diastole
             ! pressure will be computed from an imposed plot in update_cardiac_functions
             ! keep v_new as it was
             v_new(2) = unksys(2)
          end if


          dvoldt(1) = (v_new(1)-unksys(1))/dt_in
          dvoldt(2) = (v_new(2)-unksys(2))/dt_in

          unksys(1) = v_new(1) !<----- from left ventricle 
          unksys(2) = v_new(2) !<----- from right ventricle

       end if
    end if


  end subroutine sysnet_update_ventricular_synthetic_volumes


  !
  ! unk = unk + ( k_1 / 6 + k_2 / 3 + k_3 / 3 + k_4 / 6 ) * dt
  !
  subroutine update_systemic_unknowns(dt_in)
    implicit none
    real(rp) :: dt_in
    integer(ip) :: iunki

    do iunki= 1,12
       unksys(iunki) = unksys(iunki) &
            + ( &
            rhs(iunki,1) / 6.0_rp + &
            rhs(iunki,2) / 3.0_rp + &
            rhs(iunki,3) / 3.0_rp + &
            rhs(iunki,4) / 6.0_rp  &
            ) * dt_in
    end do

  end subroutine update_systemic_unknowns

  subroutine update_derived_functions(unk_in)
    implicit none
    real(rp) :: unk_in(12)

    ! chamber_elastic_energy( 4 )    ! (la, lv, ra, rv)
    ! network_elastic_energy(2,2)    ! (arterial, venous) x (systemic, pulmonar)
    ! network_kinetic_energy(2,2)    ! (arterial, venous) x (systemic, pulmonar)
    ! power_diss_network(2,2)        ! (arterial, venous) x (systemic, pulmonar)
    ! power_diss_valves( 4 )         ! (mit, aor, tri, pul)

    power_diss_valves(1) = &
         -(funcar(3) - funcar(1))**2.0_rp  / resvalve(1) ! mit
    power_diss_valves(2) = &
         -(funcar(1) - unksys(5))**2.0_rp  / resvalve(2) ! aor
    power_diss_valves(3) = &
         -(funcar(4) - funcar(2))**2.0_rp  / resvalve(3) ! tri
    power_diss_valves(4) = &
         -(funcar(2) - unksys(7))**2.0_rp  / resvalve(4) ! pul
         
    power_diss_network(1,1) = &
         - in_rsysar  * unk_in( 9) * unk_in( 9)     ! (arterial, systemic)
    power_diss_network(2,1) = &
         - in_rsysven * unk_in(10) * unk_in(10)     ! (venous, systemic)
    power_diss_network(1,2) = &
         - in_rpular  * unk_in(11) * unk_in(11)     ! (arterial, pulmonar)
    power_diss_network(2,2) = &
         - in_rpulven * unk_in(12) * unk_in(12)     ! (venous, pulmonar)

    chamber_elastic_energy(1) =  &
         in_elapass * (unk_in(1) - in_vola) * (unk_in(1) - in_vola) / 2.0_rp  ! la
!    chamber_elastic_energy(2) =  &
!         in_elvpass * (unk_in(2) - in_volv) * (unk_in(2) - in_volv) / 2.0_rp  ! lv
    chamber_elastic_energy(3) =  &
         in_erapass * (unk_in(3) - in_vora) * (unk_in(3) - in_vora) / 2.0_rp  ! ra
!    chamber_elastic_energy(4) =  &
!         in_ervpass * (unk_in(4) - in_vorv) * (unk_in(4) - in_vorv) / 2.0_rp  ! rv

    network_elastic_energy(1,1) = &
         in_csysar  * unk_in(5) * unk_in(5) / 2.0_rp    ! (arterial, systemic)
    network_elastic_energy(2,1) = &
         in_csysven * unk_in(6) * unk_in(6) / 2.0_rp    ! (venous, systemic)
    network_elastic_energy(1,2) = &
         in_cpular  * unk_in(7) * unk_in(7) / 2.0_rp    ! (arterial, pulmonar)
    network_elastic_energy(2,2) = &
         in_cpulven * unk_in(8) * unk_in(8) / 2.0_rp    ! (venous, pulmonar)

    network_kinetic_energy(1,1) = &
         in_lsysar  * unk_in( 9) * unk_in( 9) / 2.0_rp    ! (arterial, systemic)
    network_kinetic_energy(2,1) = &
         in_lsysven * unk_in(10) * unk_in(10) / 2.0_rp    ! (venous, systemic)
    network_kinetic_energy(1,2) = &
         in_lpular  * unk_in(11) * unk_in(11) / 2.0_rp    ! (arterial, pulmonar)
    network_kinetic_energy(2,2) = &
         in_lpulven * unk_in(12) * unk_in(12) / 2.0_rp    ! (venous, pulmonar)
    
  end subroutine update_derived_functions
  
  !
  ! Compute cardiac functions for the next time step
  !

  subroutine update_cardiac_functions(ela_in,unk_in,funcar_out,funperi_out,debuni)
!    use def_master,           only :  kfl_paral
    implicit none
    real(rp) :: ela_in(4)
    real(rp) :: unk_in(12)
    real(rp) :: funcar_out(8)
    real(rp) :: funperi_out
    real(rp) :: resratio(2)
    real(rp) :: xpres

    integer(ip) :: iline_pp, icavity,ires_valve,debuni

    real(rp), save :: isovol_constant(2,2)

    ! unksys (systemic unknowns):
    !  1: V_lv, <----- from left ventricle 
    !  2: V_rv, <----- from right ventricle
    !  3: V_la, 
    !  4: V_ra, 
    !  5: P^system_arterial, 
    !  6: P^system_venous, 
    !  7: P^pulmonary_arterial, 
    !  8: P^pulmonary_venous, 
    !  9: Q^system_arterial, 
    ! 10: Q^system_venous, 
    ! 11: Q^pulmonar_arterial, 
    ! 12: Q^pulmonar_venous, 

    ! funcar (cardiac functions):
    !  1: P_lv, (1)
    !  2: P_rv, (3)
    !  3: P_la, (2)
    !  4: P_ra, (4)
    !  5: Q_mv, (mitral)
    !  6: Q_av, (aortic)
    !  7: Q_tv, (tricuspid)
    !  8: Q_pv, (pulmonic)

    ! external pericardial pressure, could be constant or varying with respiratory cycle:
    !     funperi

    ! elastan (chambers elastances)
    !  1: E_lv (1)
    !  2: E_rv (3)
    !  3: E_la (2)
    !  4: E_ra (4)

    ! resvalve (cardiac valves flow resistances)
    !  1: R_mv, (mitral)
    !  2: R_av, (aortic)
    !  3: R_tv, (tricuspid)
    !  4: R_pv, (pulmonic)

    ! pericardiac pressure (the reference value)
    funperi_out = 0.0_rp

    ! atrial pressures, from model
    funcar_out(3) = funperi + 10.0_rp + ela_in(3) * (unk_in(3) - in_vola)
    funcar_out(4) = funperi +  5.0_rp + ela_in(4) * (unk_in(4) - in_vora)

    !!    funcar_out(3) = 0.0_rp   ! to open in pla (left atria)
    !!    funcar_out(4) = 0.0_rp   ! to open in pra (right atria)

    if (adjustl(trim(ippqp)) == 'qp') then

       do icavity= 1,2

          ires_valve= (icavity-1)*2_ip

          !   resvalve(1+ires_valve)     ! mitral / triscupid closed
          !   resvalve(2+ires_valve)     ! aortic / pulmonic closed

          if ((adjustl(trim(iphase_sysnet(icavity)))== "sys"))then
             isovol_constant(1,icavity) = funcar(icavity)
             isovol_constant(2,icavity) = (funcar(icavity) - funcar(2+icavity))/in_tisosys
             ! ventricular pressures, to be send to alya
             resratio(icavity) = (1.0_rp / resvalve(1+ires_valve) + 1.0_rp / resvalve(2+ires_valve))
             resratio(icavity) = 1.0_rp / resratio(icavity)
             funcar_out(icavity) = resratio(icavity) &
                  * (funcar_out(2+icavity)/resvalve(1+ires_valve) &
                  + unksys(5+ires_valve)/resvalve(2+ires_valve) - dvoldt(icavity) )

!             write (6,*) 
!             write (6,*) 'ooooj ', kfl_paral,icavity
!             write (6,*) 'ooooj ', funcar_out(icavity), funcar_out(2+icavity)
!             write (6,*) 

                 
          else if ((adjustl(trim(iphase_sysnet(icavity)))== "isy")) then
             if (adjustl(trim(isolo)) == 'solo') then
                ! startup phase
                funcar_out(icavity) = &
                     + (isovol_constant(1,icavity) - isovol_constant(2,icavity) * chron_phase(2,icavity))
             else if (adjustl(trim(isolo)) == 'alya') then
                ! ventricular pressures, to be send to alya
                resratio(icavity) = (1.0_rp / resvalve(1+ires_valve) + 1.0_rp / resvalve(2+ires_valve))
                resratio(icavity) = 1.0_rp / resratio(icavity)
                funcar_out(icavity) = resratio(icavity) &
                     * (funcar_out(2+icavity)/resvalve(1+ires_valve) &
                     + unksys(5+ires_valve)/resvalve(2+ires_valve) - dvoldt(icavity) )                
             end if
          else if ((adjustl(trim(iphase_sysnet(icavity)))== "dia"))then
             isovol_constant(1,icavity) = funcar(icavity)
             isovol_constant(2,icavity) = (funcar(icavity) - unksys(5+ires_valve))/in_tisodia
             ! ventricular pressures, to be send to alya
             resratio(icavity) = (1.0_rp / resvalve(1+ires_valve) + 1.0_rp / resvalve(2+ires_valve))
             resratio(icavity) = 1.0_rp / resratio(icavity)
             funcar_out(icavity) = resratio(icavity) &
                  * (funcar_out(2+icavity)/resvalve(1+ires_valve) &
                  + unksys(5+ires_valve)/resvalve(2+ires_valve) - dvoldt(icavity) )

          else if ((adjustl(trim(iphase_sysnet(icavity)))== "idi")) then
             if (adjustl(trim(isolo)) == 'solo') then
                ! startup phase
                funcar_out(icavity) = &
                     + (isovol_constant(1,icavity) - isovol_constant(2,icavity) * chron_phase(4,icavity))
             else if (adjustl(trim(isolo)) == 'alya') then
                ! ventricular pressures, to be send to alya
                resratio(icavity) = (1.0_rp / resvalve(1+ires_valve) + 1.0_rp / resvalve(2+ires_valve))
                resratio(icavity) = 1.0_rp / resratio(icavity)
                if (icavity == 1) then
                   write (6,*)
                   write (6,*) 'iphase ',iphase_sysnet(1:2)
                end if
                
                funcar_out(icavity) = resratio(icavity) &
                     * (funcar_out(2+icavity)/resvalve(1+ires_valve) &
                     + unksys(5+ires_valve)/resvalve(2+ires_valve) - dvoldt(icavity) )
                
             end if
                
          end if
       end do
    else if (adjustl(trim(ippqp)) == 'pp') then

       iline_pp = int(real(in_nline_pp) * chron_cycle(1) / in_tcycle)

       xpres= pressleft(2,iline_pp+1) &
            + (chron_cycle(1)-pressleft(1,iline_pp)) * (pressleft(2,iline_pp+1)-pressleft(2,iline_pp)) &
            / (pressleft(1,iline_pp+1)-pressleft(1,iline_pp))
       funcar_out(1) = xpres

       xpres= pressright(2,iline_pp+1) &
            + (chron_cycle(1)-pressright(1,iline_pp)) * (pressright(2,iline_pp+1)-pressright(2,iline_pp)) &
            / (pressright(1,iline_pp+1)-pressright(1,iline_pp))
       funcar_out(2) = xpres


    end if

!                   write (6,*)
!                   write (6,*) 'reses',resvalve(1:4)

        

    ! Compute the flowrates Q
    funcar_out(5) = (funcar_out(3) - funcar_out(1)) / resvalve(1)
    funcar_out(6) = (funcar_out(1) - unk_in(5)) / resvalve(2)
    funcar_out(7) = (funcar_out(4) - funcar_out(2)) / resvalve(3)
    funcar_out(8) = (funcar_out(2) - unk_in(7)) / resvalve(4)


  end subroutine update_cardiac_functions

  subroutine update_cardiac_elastances(ela_out)
    implicit none
    real(rp) :: ela_out(4)
    real(rp) :: t_atrium
    real(rp) :: yfun(4),Tinitial_relax,Tfinal_contraction
    
    ! elastan (chambers elastances)
    !  1: E_lv (1)
    !  2: E_rv (3)
    !  3: E_la (2)
    !  4: E_ra (4)
    !    in_lvtimact
    !    in_rvtimact
    !    in_tchamber(1,1) = Tac
    !    in_tchamber(2,1) = Tar
    !    in_tchamber(3,1) = tac
    !    in_tchamber(4,1) = tar

    ! atria (both right and left are the same, using left atrium parameters)

    yfun(3)= 0.0_rp
    yfun(4)= 0.0_rp


    Tinitial_relax= in_tchamber(4,1) + in_tchamber(2,1) - in_tcycle
    Tfinal_contraction= in_tchamber(3,1) + in_tchamber(1,1)
    
    if (chron_cycle(1) < Tinitial_relax) then       
       t_atrium = (chron_cycle(1) + in_tcycle - in_tchamber(4,1)) / in_tchamber(2,1)
       yfun(3) = (1.0_rp + cos(3.1415926_rp * t_atrium)) / 2.0_rp              

       yfun(4) = yfun(3) ! function is the same for both atria
       
!       write(150,*) chron_cycle(1), '  1.0  ', yfun(3), ela_out(3)

    else if ((Tinitial_relax < chron_cycle(1)) .and. (chron_cycle(1) < in_tchamber(3,1))) then
       yfun(3) = 0.0_rp

       yfun(4) = yfun(3) ! function is the same for both atria

!       write(150,*) chron_cycle(1), '  2.0  ', yfun(3), ela_out(3)

    else if ((in_tchamber(3,1) < chron_cycle(1)) .and. (chron_cycle(1) < Tfinal_contraction)) then
       t_atrium = (chron_cycle(1) - in_tchamber(3,1)) / in_tchamber(1,1)
       yfun(3) = (1.0_rp - cos(3.1415926_rp * t_atrium)) / 2.0_rp       

       yfun(4) = yfun(3) ! function is the same for both atria

!       write(150,*) chron_cycle(1), '  3.0  ', yfun(3), ela_out(3)
       
    else if ((Tfinal_contraction < chron_cycle(1)) .and. (chron_cycle(1) < in_tcycle)) then       
       t_atrium = (chron_cycle(1) - in_tchamber(4,1)) / in_tchamber(2,1)
       yfun(3) = (1.0_rp + cos(3.1415926_rp * t_atrium)) / 2.0_rp

       yfun(4) = yfun(3) ! function is the same for both atria

!       write(150,*) chron_cycle(1), '  4.0  ', yfun(3), ela_out(3)

    end if

    ela_out(3) = in_elapass + in_elamaxact * yfun(3) 
    ela_out(4) = in_erapass + in_eramaxact * yfun(4) 
        
  end subroutine update_cardiac_elastances

  
  !
  ! Compute rhs for the next time step
  !

  subroutine compute_rhs(irk,unk_in)
    implicit none
    integer(ip) :: irk
    real(rp) :: unk_in(12)
    
    ! unksys (systemic unknowns):
    !  1: V_lv, (2) <----- from left ventricle 
    !  2: V_rv, (4) <----- from right ventricle
    !  3: V_la, (1)
    !  4: V_ra, (3)
    !  5: P^system_arterial, 
    !  6: P^system_venous, 
    !  7: P^pulmonary_arterial, 
    !  8: P^pulmonary_venous, 
    !  9: Q^system_arterial, 
    ! 10: Q^system_ven, 
    ! 11: Q^pulmonar_arterial, 
    ! 12: Q^pulmonar_venous, 

    ! funcar (cardiac functions):
    !  1: P_lv, (1)
    !  2: P_rv, (3)
    !  3: P_la, (2)
    !  4: P_ra, (4)
    !  5: Q_mv, (mitral)
    !  6: Q_av, (aortic)
    !  7: Q_tv, (tricuspid)
    !  8: Q_pv, (pulmonic)


    if (adjustl(trim(ippqp)) == 'qp') then
       rhs( 1,irk)= 0.0_rp ! <----- from left ventricle, do not update
       rhs( 2,irk)= 0.0_rp ! <----- from right ventricle, do not update       
    else if (adjustl(trim(ippqp)) == 'pp') then
       rhs( 1,irk) = funcar(5) - funcar(6)    ! <----- update left ventricle
       rhs( 2,irk) = funcar(7) - funcar(8)    ! <----- update right ventricle
    end if
       
    !    rhs( 3) = (unk_in%qpulven )- funcar%qmv
    rhs( 3,irk) = &
         unk_in(12) - funcar(5)

    !    rhs( 4) = unk_in%qsysven - funcar%qtv
    rhs( 4,irk) = &
         unk_in(10) - funcar(7)    

    !    rhs( 5) = (funcar%qav - unk_in%qsysar) / in_csysar
    rhs( 5,irk) = &
         (funcar(6) - unk_in(9)) / in_csysar    

    !    rhs( 6) = (unk_in%qsysar - unk_in%qsysven) / in_csysven
    rhs( 6,irk) = &
         (unk_in(9) - unk_in(10)) / in_csysven    

    !    rhs( 7) = (funcar%qpv - unk_in%qpular) / in_cpular
    rhs( 7,irk) = &
         (funcar(8) - unk_in(11)) / in_cpular    

    !    rhs( 8) = (unk_in%qpular - unk_in%qpulven) / in_cpulven
    rhs( 8,irk) = &
         (unk_in(11) - unk_in(12)) / in_cpulven    

    !    rhs( 9) = - (unk_in%qsysar + (unk_in%psysven-unk_in%psysar)/in_rsysar  ) * in_rsysar / in_lsysar
    rhs( 9,irk) = &
         - (unk_in( 9) + (unk_in(6) - unk_in( 5)) / in_rsysar ) * in_rsysar / in_lsysar   

    !    rhs(10) = - (unk_in%qsysven + (funcar%pra-unk_in%psysven)/in_rsysven ) * in_rsysven / in_lsysven
    rhs(10,irk) = &
         - (unk_in(10) + (funcar(4) - unk_in( 6)) / in_rsysven ) * in_rsysven / in_lsysven   

    !    rhs(11) = - (unk_in%qpular + (unk_in%ppulven-unk_in%ppular)/in_rpular  ) * in_rpular / in_lpular
    rhs(11,irk) = &
         - (unk_in(11) + (unk_in(8) - unk_in( 7)) / in_rpular ) * in_rpular / in_lpular   

    !    rhs(12) = - (unk_in%qpulven + (funcar%pla-unk_in%ppulven)/in_rpulven ) * in_rpulven / in_lpulven
    rhs(12,irk) = &
         - (unk_in(12) + (funcar(3) - unk_in( 8)) / in_rpulven ) * in_rpulven / in_lpulven   
    
  end subroutine compute_rhs

  !
  ! Compute pvloop, either as an initialization loop or at each time step
  !
  subroutine sysnet_compute_pvloop(itask)

    use def_master,           only : ITASK_INITIA,ITASK_ENDITE, TIME_N
    use def_master,           only : dtime
    implicit none
    integer(ip) :: itask
    integer(ip) :: coun_iter,coun_cycle(2),ifloop,debuni

    real(rp)    :: del_t, time_total, xauxi

    character(30):: iwrite

    select case(itask)       

    case (ITASK_INITIA)

       ! THIS IS THE STARTUP PHASE
       ! SYSNET is called in sld_inivar
       
       
       del_t         = in_startup_deltime        ! dt for the initialization run
       time_total    = in_startup_fintime        ! final time of the initialization run

       ifloop        = 1_ip            ! the loop is on
       coun_iter     = 0_ip            ! iterations counter
       debuni        = 2_ip            ! debugger level

       chron_phase   = 0.0_rp        ! phase chronometer
       chron_cycle   = 0.0_rp        ! cycle chronometer
       chron_total   = 0.0_rp        ! total chronometer

       unksys(1)= cavities_sysnet(1)%volume(1) ! left ventricle initial volume
       unksys(2)= cavities_sysnet(2)%volume(1) ! right ventricle initial volume
       in_maxvolven_synth(1)= unksys(1)
       in_maxvolven_synth(2)= unksys(2)

       if( unksys(1) == 0.0_rp ) then
          call runend('MOD_KER_SYSNET: Initial volme zero, the nodes defining the cavity 1 plane are wrongly ordered.')
       end if
       if( unksys(2) == 0.0_rp ) then
          call runend('MOD_KER_SYSNET: Initial volme zero, the nodes defining the cavity 2 plane are wrongly ordered.')
       end if
       
       xauxi= (unksys(1) + unksys(2)) * in_startup_ejection

       in_minvolven_synth(1)= in_maxvolven_synth(1) - xauxi / 2.0_rp
       in_minvolven_synth(2)= in_maxvolven_synth(2) - xauxi / 2.0_rp

       ! initialize counters and chrons
       coun_iter        = coun_iter      + 1_ip
       chron_total      = chron_total    + del_t
       chron_cycle(1)   = chron_cycle(1) + del_t
       chron_cycle(2)   = chron_cycle(2) + del_t       
       chron_preload    = 0.0_rp
       coun_cycle       = 0_ip

       ! No startup required
       if ((coun_cycle(1) .ge. in_startup_cycles) .or. (chron_total  > time_total)) then
          coun_iter = 0_ip
          chron_total      = chron_total    - del_t
          chron_cycle(1)   = chron_cycle(1) - del_t
          chron_cycle(2)   = chron_cycle(2) - del_t       
          return
       end if

       call sysnet_update_cycle(del_t,debuni,coun_cycle)  
       if (chron_total  > time_total) ifloop=0_ip

       iwrite='start'
       call sysnet_write_startup_res(&
            iwrite,coun_iter)       

       do while (ifloop == 1)

          ! compute synthetic ventricular volumes
          call sysnet_update_ventricular_synthetic_volumes(del_t)

          ! update functions
          call sysnet_update_functions(del_t,debuni)     

          ! RK update unksys
          call sysnet_time_iteration_step(del_t)

          iwrite='write'
          call sysnet_write_startup_res(&
               iwrite,coun_iter)

          ! update counters and chrons
          coun_iter        = coun_iter      + 1_ip
          chron_total      = chron_total    + del_t
          chron_cycle(1)   = chron_cycle(1) + del_t
          chron_cycle(2)   = chron_cycle(2) + del_t
          call sysnet_update_cycle(del_t,debuni,coun_cycle)          
          if (chron_total  > time_total) ifloop=0_ip
          ! the LV controls the coun_cycle
          if ((coun_cycle(1) == in_startup_cycles)) ifloop=0_ip   
          
       end do

    case (ITASK_ENDITE)

       ! THIS IS THE REGULAR PHASE
       ! SYSNET is called in sld_endite
       
       ! !!!!! OJO, SOLO PARA DEBUGUIAR!!
       !       cavities_sysnet(1) % volume(TIME_N)= unksys(1)* 0.000001_rp
       !       cavities_sysnet(2) % volume(TIME_N)= unksys(2)* 0.000001_rp       
       ! !!!

!!       write(6,*) 'holaaaaa', resvalve(1:2)
       call sysnet_update_cycle(dtime,debuni,coun_cycle)  
!!       write(6,*) 'holaaaaa', resvalve(1:2)

       ! update functions
       call sysnet_update_functions(dtime,debuni)     

       ! RK update unksys
       call sysnet_time_iteration_step(dtime)

    end select


  end subroutine sysnet_compute_pvloop

end module mod_ker_sysnet

