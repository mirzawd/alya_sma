!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_ker_tendencies
  use def_kintyp, only : rp, ip,lg

  implicit none
  save
  logical(lg)    , public               :: kfl_tendencies_ker

  ! pointers
  real(rp),    public, pointer      :: ten_ugeos(:,:,:)  ! 2,nz,nsteps
  real(rp),    public, pointer      :: ten_uadve(:,:,:)  ! 2,nz,nsteps
  real(rp),    public, pointer      :: ten_thadv(:,:)    ! nz, nsteps
  real(rp),    public, pointer      :: ten_thmeso(:,:)    ! nz, nsteps
  real(rp),    public, pointer      :: ten_umeso(:,:,:)    ! nz, nsteps
  real(rp),    public, pointer      :: ten_heigh(:,:)    ! nz, nsteps
  real(rp),    public, pointer      :: ten_time(:)       ! time in secs, nsteps
  integer(ip), public, pointer      :: ten_nzlev(:)      ! number of zlevels nsteps

  ! module variables
  integer(ip), public               :: ten_nstep        ! number of steps
  integer(ip), public               :: ten_nzlev_max    ! number of vertical levels

  integer(ip), public               :: up_istep_last = 1_ip ! initialization

  ! functions
  public  :: read_tendencies           ! reads tendencies and allocates structures
  public  :: get_tendencies_u          ! interpolate tendency value to zcoor height and ctime
  public  :: get_tendencies_th         ! interpolate T WRF    value to zcoor height and ctime
  public  :: get_tendencies_uwrf       ! interpolate U WRF    value to zcoor height and ctime
  public  :: ker_tendencies_parall     ! interpolate tendency value to zcoor height and ctime
  private :: get_tendencies_u_istep    ! interpolate tendency value to zcoor height, for step istep
  private :: tendencies_allocate       ! allocates tendencies (each processor)
contains
  !
  ! IMPORTANTE
  ! para no depender tanto de AlyaFix, podria interpolarse WRF a la malla de Alya (solo dependiente de altura)
  ! bastaria con leer desde el mismo archivo, ademas de tendencias y Th, a U y W. Pero habria que saber tratar a U y W.
  ! hay que mirar tambien si AlyaFix lee las condiciones iniciales.
  !
  subroutine  read_tendencies()
    ! reads tendencies and allocates structures
    use def_inpout,              only : words, param, exists, getint, getrea
    use mod_ecoute,              only : ecoute
   
    implicit none

    integer(ip)    :: istep, izlev

    ten_nzlev_max = 500_ip ! initialization
    
    do while( words(1) /= 'FIELD' )
       call ecoute('ker_readat')  ! STEPS
    end do
    ! number of steps
    ten_nstep =  getint('STEPS',0_ip,'#Number of time steps tendencies')
    ten_nzlev_max = getint('NZLEV',0_ip,'#Number of vertical levels tendencies')
    if (ten_nstep == 0_ip.or.ten_nzlev_max==0_ip) then
      ! write(lulog, *) 'nsteps=' , ten_nstep,'nzlevels=', ten_nzlev_max
       write(    *, *) 'nsteps=' , ten_nstep,'nzlevels=', ten_nzlev_max
       call runend('mod_ker_tendencies: number of vertical levels and time steps must be larger than zero')
    end if
    ! 
    ! allocate structures
    !
    call tendencies_allocate()
   

! read data from file (/gpfs/scratch/bsc21/bsc21811/Alya-Runs/Hornamossen> more Hornamossen_tendencies/Hornamossen_tendencies)
    istep = 0    
    call ecoute('ker_readat')   ! STEP istep TIME itime
    do while(words(1) /= 'ENDFI' )
       if(words(2) == 'STEP ' ) then
          istep = istep +1_ip       
 !         print *, 'reading step', istep
 !         if (getint('STEP ',0_ip,'#step number tendencies')/=istep) call runend('mod_ker_tendencies: error reading, step number does not match')
          ten_time(istep) = getrea('TIME ', 0.0_rp,'#time value')  
          
          izlev = 0_ip
          call ecoute('ker_readat')   ! labels
          call ecoute('ker_readat')   ! data
          do while(words(2) /= 'ENDST' )  ! we need also the inlet profile! but it should be read in other routine
             izlev = izlev +1_ip 
             ten_heigh(izlev, istep)     = param(1) ! should be given in increasing order
             ten_ugeos(1:2,izlev, istep) = param(2:3)
             ten_uadve(1:2,izlev, istep) = param(4:5)
             ten_thadv(izlev,istep)      = param(6)  ! not used
             ten_thmeso(izlev,istep)      = param(6) ! Mesoscalar Temper meso
             ten_umeso(1,izlev,istep)    = param(7) ! Mesoscalar Veloc x 
             ten_umeso(2,izlev,istep)    = param(8) ! Mesoscalar Veloc y
             call ecoute('ker_readat') ! header
          end do
          ten_nzlev(istep) = izlev ! quantity of zlevels
         
          if (ten_nzlev_max.ne.izlev) then
             print *, 'nzlev, max_nzlevs', izlev, ten_nzlev_max
             call runend('mod_ker_tendencies_ bad number of levels')
          end if
          call ecoute('ker_readat')
        end if
    end do      
    return
    
  end subroutine read_tendencies
!
!***********************************************************************************************************
!
  subroutine  get_tendencies_u_istep(zcoor,istep,gpten_geo, gpten_uadv, velmeso, thmeso)
    ! interpolate tendencies value to zcoor height, for step istep
    !
    !  for interpolation, height values are supposed to be given in increasing order
    !  Two options : 1) interpolates tendencies: present gpten_geo and gpten_uadv and not velmeso thmeso
    !                2) interpolates wrf       : present velmeso and/or thmeso and not gpen_geo and gpten _uadv
    implicit none
    real(rp), intent(in)              :: zcoor   ! z coord
    real(rp), intent(out) ,optional   :: gpten_geo(2), gpten_uadv(2)
    real(rp), intent(out), optional   :: velmeso(2), thmeso
    integer(ip), intent(in)           :: istep
    ! local variables
    integer(ip)                       ::  iz, jz, kz
    real(rp)                          ::  facto

    

    ! for a given time step
    if (zcoor.lt.ten_heigh(1, istep)) then ! coorditate less than minimum 
       gpten_geo =0.0_rp
       gpten_uadv = 0.0_rp
!    else if(zcoor.gt.z_top )  then 
    else
       iz =1_ip
       jz =ten_nzlev(istep)
       kz =jz/2_ip
       do while ((jz-iz).gt.1_ip)  ! bisection method
          if (zcoor.lt.ten_heigh(kz,istep)) then
             jz = kz                  
          else
             iz = kz
          end if
          kz = (iz+jz)/2_ip
       end do
       !interpolation factor
       facto = (zcoor-ten_heigh(iz, istep))/(ten_heigh(jz, istep) -  ten_heigh(iz, istep))
       if (present(gpten_geo)) then
          gpten_geo(1:2)  = (1.0_rp - facto)* ten_ugeos(1:2,iz, istep) + facto* ten_ugeos(1:2,jz, istep)  
          gpten_uadv(1:2) = (1.0_rp - facto)* ten_uadve(1:2,iz, istep) + facto* ten_uadve(1:2,jz, istep)
          return
       end if
       if (present(velmeso)) velmeso(1:2) =  (1.0_rp - facto)* ten_umeso(1:2,iz, istep)  + facto* ten_umeso(1:2,jz, istep)
       if (present(thmeso))  thmeso  =  (1.0_rp - facto)* ten_thmeso(iz, istep) + facto* ten_thmeso(jz, istep)
    end if
    
  end subroutine get_tendencies_u_istep
  
!
!***********************************************************************************************************
!
 subroutine  get_tendencies_u(zcoor,ctime,gpten_geo, gpten_uadv)
    ! interpolate tendency value to zcoor height and ctime
    !
    !
    implicit none
    real(rp), intent(in)    :: zcoor, ctime   ! z coord
    real(rp), intent(out)   :: gpten_geo(2), gpten_uadv(2)

    real(rp)                :: geo_old(2), uadv_old(2), geo_new(2), uadv_new(2)
    real(rp)                :: facto
    !
    ! time interpolation
    !
    do while (ctime.gt.ten_time(up_istep_last))  ! check if new time limits
         up_istep_last= up_istep_last +1_ip
    end do
   
    facto = (ctime - ten_time(up_istep_last-1_ip))/  (ten_time(up_istep_last) - ten_time(up_istep_last-1_ip))

    ! interpolate tendencies at a given step istep
    call get_tendencies_u_istep(zcoor, up_istep_last -1_ip, geo_old, uadv_old)
    call get_tendencies_u_istep(zcoor, up_istep_last   , geo_new, uadv_new)

    ! time interpolation
    gpten_geo = facto*geo_new  + (1.0_rp-facto)*geo_old
    gpten_uadv = facto*uadv_new  + (1.0_rp-facto)*uadv_old
    
  end subroutine get_tendencies_u
!
!***********************************************************************************************************
!
 subroutine  get_tendencies_th(zcoor,ctime, gptem_wrf)
    ! interpolate tendency value to zcoor height and ctime
    !
    !
    implicit none
    real(rp), intent(in)    :: zcoor, ctime   ! z coord
    real(rp), intent(out)   :: gptem_wrf

    real(rp)                :: tewrf_old,  tewrf_new
    real(rp)                :: facto
    !
    ! time interpolation
    !
    do while (ctime.gt.ten_time(up_istep_last))  ! check if new time limits
         up_istep_last= up_istep_last +1_ip
    end do
   
    facto = (ctime - ten_time(up_istep_last-1_ip))/  (ten_time(up_istep_last) - ten_time(up_istep_last-1_ip))

    ! interpolate tendencies at a given step istep
    call get_tendencies_u_istep(zcoor, up_istep_last -1_ip, thmeso=tewrf_old)
    call get_tendencies_u_istep(zcoor, up_istep_last   , thmeso=tewrf_new)

    ! time interpolation
    gptem_wrf  = facto*tewrf_new   + (1.0_rp-facto)*tewrf_old 
    
  end subroutine get_tendencies_th
!
!***********************************************************************************************************
!
 subroutine  get_tendencies_uwrf(zcoor,ctime, gpvel_wrf)
    ! interpolate tendency value to zcoor height and ctime
    !
    !
    implicit none
    real(rp), intent(in)    :: zcoor, ctime   ! z coord
    real(rp), intent(out)   :: gpvel_wrf(2)

    real(rp)                :: vewrf_old(2),  vewrf_new(2)
    real(rp)                :: facto
    !
    ! time interpolation
    !
    do while (ctime.gt.ten_time(up_istep_last))  ! check if new time limits
         up_istep_last= up_istep_last +1_ip
    end do
   
    facto = (ctime - ten_time(up_istep_last-1_ip))/  (ten_time(up_istep_last) - ten_time(up_istep_last-1_ip))

    ! interpolate tendencies at a given step istep
    call get_tendencies_u_istep(zcoor, up_istep_last -1_ip, velmeso=vewrf_old)
    call get_tendencies_u_istep(zcoor, up_istep_last   , velmeso=vewrf_new)

    ! time interpolation
    gpvel_wrf(1:2)  = facto*vewrf_new(1:2)   + (1.0_rp-facto)*vewrf_old(1:2) 
    
  end subroutine get_tendencies_uwrf
     
!
!***********************************************************************************************************
!
  subroutine  ker_tendencies_initialize_u(veloc, height)
    !initialize velocity with mesoscale values
    !
    !
    use def_domain, only : ndime, npoin
    implicit none
    real(rp), intent(out)    :: veloc(ndime, npoin)   ! z coord
    real(rp), intent(in)    :: height(npoin)

    integer(ip)   :: ipoin
    
    do ipoin=1, npoin
       call get_tendencies_u_istep(height(ipoin), 1_ip, velmeso=veloc(1:2,ipoin))
    end do   
   
    
  end subroutine ker_tendencies_initialize_u
!  
!***********************************************************************************************************
!
  subroutine  ker_tendencies_initialize_th(temper, height)
    !initialize temperature with mesoscale values
    !
    !
    use def_domain, only : npoin
    implicit none
    real(rp), intent(out)    :: temper(npoin)   ! z coord
    real(rp), intent(in)    :: height(npoin)

    integer(ip)   :: ipoin
    
    do ipoin=1, npoin
       call get_tendencies_u_istep(height(ipoin), 1_ip, thmeso=temper(ipoin))
    end do   
   
    
  end subroutine ker_tendencies_initialize_th
!
!***********************************************************************************************************
!
  
  
  subroutine tendencies_allocate()
    use def_master,              only :mem_modul, modul
    use mod_memory
    
    call memory_alloca(mem_modul(1:2,modul),'TEN_UGEOS','read_tendencies',ten_ugeos, 2_ip, ten_nzlev_max,ten_nstep) ! geost tend
    call memory_alloca(mem_modul(1:2,modul),'TEN_UADVE','read_tendencies',ten_uadve, 2_ip, ten_nzlev_max,ten_nstep) ! advective tendencies
    call memory_alloca(mem_modul(1:2,modul),'TEN_THADV','read_tendencies',ten_thadv, ten_nzlev_max,ten_nstep)       ! energy advection
    call memory_alloca(mem_modul(1:2,modul),'TEN_THMESO','read_tendencies',ten_thmeso, ten_nzlev_max,ten_nstep)     ! temper from wrf
    call memory_alloca(mem_modul(1:2,modul),'TEN_UMESO','read_tendencies',ten_umeso, 2_ip, ten_nzlev_max,ten_nstep) ! VELOC from wrf
    call memory_alloca(mem_modul(1:2,modul),'TEN_HEIGH','read_tendencies',ten_heigh, ten_nzlev_max,ten_nstep)       ! height of each vertical level (at each time step)
    call memory_alloca(mem_modul(1:2,modul),'TEN_TIME' ,'read_tendencies',ten_time , ten_nstep)                     ! time in seconds of each step  
    call memory_alloca(mem_modul(1:2,modul),'TEN_NZLEV','read_tendencies',ten_nzlev, ten_nstep)                     ! number of vertical levels at each time step) 
  end subroutine tendencies_allocate

  
  subroutine ker_tendencies_parall()
    use def_master,   only : ISLAVE
    use mod_exchange, only : exchange_init
    use mod_exchange, only : exchange_add
    use mod_exchange, only : exchange_end

    call exchange_init()
    call exchange_add(ten_nzlev_max)
    call exchange_add(ten_nstep)
    call exchange_end()

    if( ten_nzlev_max > 0 .and. ten_nstep > 0 ) then 
       if( ISLAVE) call tendencies_allocate()
       call exchange_init()
       call exchange_add(ten_time)
       call exchange_add(ten_heigh)
       call exchange_add(ten_ugeos)
       call exchange_add(ten_uadve)
       call exchange_add(ten_thadv)
       call exchange_add(ten_thmeso)
       call exchange_add(ten_umeso)
       call exchange_add(ten_nzlev)
       call exchange_end()
    end if
    
  end subroutine ker_tendencies_parall
end module mod_ker_tendencies
