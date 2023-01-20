!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_interp_tab
  !-----------------------------------------------------------------------
  !****f* mathru/interp_tab
  ! NAME 
  !    interp_tab
  ! DESCRIPTION
  !    Interpolate in arbitrary table
  ! USES
  ! USED BY
  !    chm_reatab: get Yc scaling values, get chi shape function
  !    pts_thermodynamic: get properties for evaporation model
  !***
  !-----------------------------------------------------------------------

  use def_kintyp_basic,   only : ip,rp
  use mod_messages,       only : messages_live
  use mod_strings,        only : integer_to_string
  use mod_memory_basic,   only : memory_alloca
  use mod_memory_basic,   only : memory_deallo
  use def_master,         only : mem_modul,modul
  use mod_communications, only : PAR_MAX
  use def_master,         only : INOTSLAVE
  use def_master,         only : mem_modul, modul
  use def_inpout,         only : param,words,nunit
  use def_kintyp,         only : ip,rp
  use mod_ecoute,         only : ecoute
  use mod_memory,         only : memory_alloca, memory_deallo
  use mod_strings,        only : upper_case
  use iso_fortran_env,    only : iostat_eor,iostat_end     
  use mod_exchange,       only : exchange_init
  use mod_exchange,       only : exchange_add
  use mod_exchange,       only : exchange_end
  
  implicit none
  !----------------------------------------------------------------------
  ! Lookup table type:  
  !----------------------------------------------------------------------
  type typ_tab_coord
      character(5)              :: name             ! Name of coordinate 
      integer(ip)               :: leng             ! Number of elements along coordinate
      real(rp) ,       pointer  :: x(:)             ! Discrete coordinates
  end type typ_tab_coord

  type typ_lookup_table
      integer(ip)                   :: ndim, &          ! Number of independent coordinates
                                       nvar, &          ! Number of variables in table
                                       order            ! Order of interpolation
      character(5),        pointer  :: varname(:)       ! Name of variables 
      character(30),       pointer  :: longvarname(:)   ! Full variable name in case it is relevant
      type(typ_tab_coord), pointer  :: coords(:)        ! Vector of coordinates
      real(rp),            pointer  :: tab(:,:)         ! Table (len(1)*len(2)*... times nVar array)
      integer(ip),         pointer  :: shap(:)          ! Shape of table with better locality
  end type typ_lookup_table

  type typ_lookup_fw_scaling
      type(typ_lookup_table), pointer :: tab            ! scaling tables
      real(rp)                        :: lims(2_ip)     ! scaling limits 
      integer(ip),            pointer :: i_main(:)      ! indeces of scaling table in main table
      integer(ip)                     :: kfl_tab        ! external index of table
  end type typ_lookup_fw_scaling

  type typ_lookup_framework
      type(typ_lookup_table),  pointer :: main_table     ! Main table to use 
      integer(ip),             pointer :: kfl_scale(:)   ! Scalig methods:
                                                         !-2:        un-initialized
                                                         !-1:        leave as 0
                                                         ! 0:        no scaling
                                                         ! 1:        scaling with table
                                                         ! 2:        scaling with limits
                                                         ! 3:        constant scaled property 
                                                         ! >100:     variance scaling with the (n-100)^th variable
                                                         ! <-100:    square scaling with the -(n+100)^th variable
      type(typ_lookup_fw_scaling), pointer :: scaling(:)
      integer(ip),                 pointer :: scale_order(:) ! Scalig order
      integer(ip)                          :: kfl_tab_main   ! external index of table
      !
      ! Flags for coupling with chemic and partis
      !
      integer(ip)                          :: kfl_needs_enthalpy 
      integer(ip),             pointer     :: kfl_chm_control(:)   ! Control variable indexes from chemic
                                                                   ! filled by chemic or partis
                                                                   ! >0: specific index of conce
                                                                   ! -1: enthalpy
                                                                   ! -2: scalar dissipation rate of mixture fraction
  end type typ_lookup_framework
 
 
  integer(ip), parameter            :: max_lookup_dim = 6                                  ! Maximum number of independent coordinates 
  integer(ip), save                 :: line_length = 15_ip                           ! Data per line in files 
  integer(ip), save                 :: precomp_i2inds(max_lookup_dim,2**max_lookup_dim,max_lookup_dim) ! Precomputed single index to multidimensional index for two layers (for hypercube), indeces: 1st: dimesnion, 2nd: i, 3rd: inds
  logical, save                     :: initialized = .false.
  real(rp)                          :: time_int1, time_int2, time_int3, time_int4, time_int5, time_look1, time_look2

  public typ_tab_coord
  public typ_lookup_table
  public typ_lookup_framework
  public tab_init_fw      
  public tab_init_tab      
  public tab_init_coord 
  public fw_allocate           
  public tab_allocate
  public fw_par_exchange
  public tab_par_exchange
  public tab_load_file
  public tab_interp
  public tab_interpt
  public fw_scale_cont_var
  public fw_lookup
  public i2inds
  public max_lookup_dim
  public print_tab_interp_times
  private

contains

    subroutine tab_interp_initialize()
        integer(ip) :: twos(max_lookup_dim), icube(max_lookup_dim), ii, jj
        !
        ! Initialize auxiliary variable for dealing with the hipercube
        !
        twos = 2_ip
        precomp_i2inds = 0_ip
        time_int1 = 0_rp
        time_int2 = 0_rp
        time_int3 = 0_rp
        time_int4 = 0_rp
        time_int5 = 0_rp
        time_look1 = 0_rp
        time_look2 = 0_rp

        do jj = 1, max_lookup_dim
            do ii = 1, 2**max_lookup_dim
                icube = 0_ip
                call i2inds(ii,jj,twos(1:jj),icube(1:jj)) 
                precomp_i2inds(:,ii,jj) = icube
            enddo
        enddo
        initialized = .true.
    end subroutine tab_interp_initialize

    subroutine tab_init_coord(coord)
       type(typ_tab_coord), intent(inout)   :: coord
       coord % name   = ''
       coord % leng   = 0 
       nullify(coord % x)
    end subroutine tab_init_coord

    subroutine tab_init_tab(tab)
       type(typ_lookup_table), intent(inout)   :: tab
       tab % ndim  = 0
       tab % nvar  = 0
       tab % order = 0
       nullify(tab % varname)
       nullify(tab % longvarname)
       nullify(tab % coords)
       nullify(tab % tab)
       nullify(tab % shap)
    end subroutine tab_init_tab

    subroutine tab_init_scale(sca)
       type(typ_lookup_fw_scaling), intent(inout)   :: sca
       nullify(sca % tab)
       nullify(sca % i_main)
    end subroutine tab_init_scale

    subroutine tab_init_fw(fw)
       type(typ_lookup_framework), intent(inout)   :: fw
       fw % kfl_needs_enthalpy = 0
       nullify(fw % main_table)
       nullify(fw % kfl_scale)
       nullify(fw % scaling)
       nullify(fw % scale_order)
       nullify(fw % kfl_chm_control)
    end subroutine tab_init_fw




    subroutine tab_allocate(itask, coords, tab, order)
        !
        ! Allocate the fields of the tables
        !
        integer(ip),                 intent(in)          :: itask
        type(typ_tab_coord),    pointer, intent(inout)   :: coords(:)
        type(typ_lookup_table), pointer, intent(inout)   :: tab
        integer(ip),        optional,intent(in)          :: order

        integer(ip)      :: ii, nrow
        integer(4)       :: istat
        integer(ip),save :: itab = 0

        select case(itask)
            case(1_ip)
                !
                ! Initialize outter structure
                !
                allocate(tab, stat=istat)
                call tab_init_tab(tab)
                if (present(order)) tab % order = order

                allocate(coords(max_lookup_dim), stat=istat)
                do ii=1,max_lookup_dim
                   call tab_init_coord(coords(ii))
                end do 

            case(2_ip)
                !
                ! Initialize inner structure
                !
                itab = itab + 1
                nrow = 1_ip
                do ii = 1_ip, tab%ndim
                    call memory_alloca(mem_modul(1:2,modul),'COORDS % X','tab_allocate', coords(ii) % x, coords(ii) % leng)
                    nrow = nrow * coords(ii) % leng
                enddo
                allocate(tab % coords(tab % ndim),stat=istat)
                call memory_alloca(mem_modul(1:2,modul),'TAB'//trim(integer_to_string(itab))//' % TAB','tab_allocate', tab % tab, tab % nvar, nrow)
                call memory_alloca(mem_modul(1:2,modul),'TAB'//trim(integer_to_string(itab))//' % VARNAME','tab_allocate', tab % varname, tab % nvar)
                call memory_alloca(mem_modul(1:2,modul),'TAB'//trim(integer_to_string(itab))//' % LONGVARNAME','tab_allocate', tab % longvarname, tab % nvar)
                call memory_alloca(mem_modul(1:2,modul),'TAB'//trim(integer_to_string(itab))//' % SHAP','tab_allocate', tab % shap, tab % ndim)
                if (tab%ndim > 0) &
                   tab % shap(1:tab%ndim) = coords(1:tab%ndim)%leng
        end select
    end subroutine tab_allocate



    subroutine fw_allocate(itask, fw, iscal)
        !
        ! Allocate the fields of typ_lookup_framework
        !
        integer(ip),                intent(in)      :: itask
        type(typ_lookup_framework), intent(inout)   :: fw
        integer(ip),  optional,     intent(in)      :: iscal

        integer(4)       :: istat
        integer(ip)      :: ii, jj, kk, jmean,ismin,ismax

        select case(itask)
            case(1_ip)
                !
                ! Initialize parts dependent on main_table
                !
                call memory_alloca(mem_modul(1:2,modul),'FW % KFL_SCALE','fw_allocate', fw % kfl_scale,fw % main_table % ndim, INIT_VALUE=-2_ip)
                call memory_alloca(mem_modul(1:2,modul),'FW % SCALE_ORDER','fw_allocate', fw % scale_order,fw % main_table % ndim)
                allocate( fw % scaling( fw % main_table % ndim ),stat=istat)
                do ii = 1_ip, fw % main_table % ndim
                    call tab_init_scale(fw % scaling(ii))
                enddo
                call memory_alloca(mem_modul(1:2,modul),'FW % KFL_CHM_CONTROL','fw_allocate', fw % kfl_chm_control,fw % main_table % ndim)

            case(2_ip)
                !
                ! Initialize parts dependent on scaling tables 
                !
                if (fw % main_table % ndim > 0) then
                    !
                    ! Indeces corresponding to main table
                    !
                    if (present(iscal)) then
                        if (iscal <= 0) then
                            ismin = iscal 
                            ismax = iscal - 1
                        else
                            ismin = iscal 
                            ismax = iscal
                        endif
                    else
                        ismin = 1_ip
                        ismax = fw % main_table % ndim
                    endif
                    do ii = ismin, ismax
                        if (fw % kfl_scale(ii) == 1) then
                            !nullify(fw % scaling(ii) % i_main)
                            call memory_alloca(mem_modul(1:2,modul),'FW % SCALING % I_MAIN','fw_allocate', fw % scaling(ii) % i_main,fw % scaling(ii) % tab % ndim )

                            do kk = 1_ip, fw % scaling(ii) % tab % ndim
                                do jj = 1_ip, fw % main_table % ndim
                                    if ( fw % scaling(ii) % tab % coords(kk) % name == fw % main_table  % coords(jj) % name)   fw % scaling(ii) % i_main(kk) = jj
                                enddo
                            enddo

                        endif 
                    enddo
                    !
                    ! Scaling order: Decide which scaling to be executed first
                    !
                    fw % scale_order = 0_ip
                    kk = 1
                    !
                    ! Codes -1, 0, 2, 3 are first priority
                    !
                    do ii = 1, fw % main_table % ndim
                        if ( fw % kfl_scale(ii) ==-1_ip .or. &
                             fw % kfl_scale(ii) == 0_ip .or. &
                             fw % kfl_scale(ii) == 2_ip .or. &
                             fw % kfl_scale(ii) == 3_ip) then
                            fw % scale_order(kk) = ii
                            kk = kk + 1
                        endif
                    enddo
                    !
                    ! variances with parents 0 or 2 are second priority 
                    !
                    do ii = 1, fw % main_table % ndim
                        if ( abs(fw % kfl_scale(ii)) > 100_ip) then
                            jmean = abs(fw % kfl_scale(ii)) - 100
                            if (fw % kfl_scale(jmean) == 0_ip .or. &
                                fw % kfl_scale(jmean) == 2_ip) then
                                fw % scale_order(kk) = ii
                                kk = kk + 1
                            endif
                        endif
                    enddo
                    !
                    ! tabulated scalings are third priority 
                    !
                    do ii = 1, fw % main_table % ndim
                        if ( fw % kfl_scale(ii) == 1_ip) then
                            fw % scale_order(kk) = ii
                            kk = kk + 1
                        endif
                    enddo
                    !
                    ! variances with parents 1 are fourth priority 
                    !
                    do ii = 1, fw % main_table % ndim
                        if ( abs(fw % kfl_scale(ii)) > 100_ip) then
                            jmean = abs(fw % kfl_scale(ii)) - 100
                            if (fw % kfl_scale(jmean) == 1_ip) then
                                fw % scale_order(kk) = ii
                                kk = kk + 1
                            endif
                        endif
                    enddo
                endif

        end select

    end subroutine fw_allocate



    subroutine tab_par_exchange(coords, tab)
        !
        ! Exchange tables
        !
        use def_master,         only : parii, npari, nparr, &
                                       nparc, parin, parre, &
                                       parch, nparc, mem_modul, & 
                                       modul, ISLAVE, IMASTER
        
        type(typ_tab_coord),    pointer, intent(inout)   :: coords(:)
        type(typ_lookup_table), pointer, intent(inout)   :: tab

        integer(ip)      :: ii, jj, kk
       
        !
        ! Allocate pointers
        !
        if ( ISLAVE ) then 
            call tab_allocate(1_ip, coords, tab)
        else
            if (.not. associated(tab)) call tab_allocate(1_ip, coords, tab)
        endif
        
        !
        ! Exchange dimensions
        !
        call exchange_init()
        call exchange_add(tab % ndim)           
        call exchange_add(tab % nvar)           
        call exchange_add(tab % order)           
        do ii = 1_ip, size(coords, KIND=ip)
           call exchange_add(coords(ii) % leng)           
           call exchange_add(coords(ii) % name)           
        enddo
        call exchange_end()   
    


        !
        ! Once dimensions are known, allocate
        !
        if ( ISLAVE ) then 
            call tab_allocate(2_ip, coords, tab)
        endif

        !
        ! Exchange coordinates and variable names
        !
        call exchange_init()
        kk = 1_ip ! lines in table
        do ii = 1_ip, tab % ndim
           kk = kk * coords(ii) % leng
           do jj = 1_ip, coords(ii) % leng
              call exchange_add(coords(ii) % x(jj)) 
           enddo
        enddo
        
        tab % coords => coords(1:tab % ndim)
        
        do jj = 1_ip, tab%nvar
           call exchange_add(tab % varname(jj)) 
           call exchange_add(tab % longvarname(jj)) 
        enddo
        call exchange_end()   
        !
        ! Exchange table body
        !
        call exchange_init()
        do jj = 1_ip, tab%nvar
           do ii = 1_ip, kk
              call exchange_add(tab % tab(jj,ii)) 
           enddo
        enddo
        call exchange_end()   
  
    end subroutine tab_par_exchange


    subroutine fw_par_exchange(itask,fw)
        !
        ! Exchange tables
        !
        use def_master,         only : parii, npari, nparr, &
                                       nparc, parin, parre, &
                                       nparc, mem_modul, & 
                                       modul, ISLAVE, IMASTER
        
        integer(ip),                intent(in)    :: itask
        type(typ_lookup_framework), intent(inout) :: fw
        integer(ip)                               :: ii



        select case(itask)
        case(1_ip)
            !
            ! Exchange kfl_tab_main
           !
           call exchange_init()
           call exchange_add(fw % kfl_tab_main) 
           call exchange_end()
           
        case(2_ip)
            !
            ! Once main table is exchanged, allocate
            !
            if ( ISLAVE ) then 
                call fw_allocate(1_ip, fw)
            endif

            !
            ! Exchange kfl_scale
            !
           call exchange_init()
           do ii = 1_ip, fw % main_table % ndim
              call exchange_add(fw % kfl_scale(ii)) 
              call exchange_add (fw % kfl_chm_control(ii)) 
           enddo
           call exchange_add (fw % kfl_needs_enthalpy) 
           call exchange_add(fw % kfl_tab_main) 
           call exchange_end()
 
           !
           ! Based on kfl_scale exchnge scaling
           !
           call exchange_init()
           do ii = 1_ip, fw % main_table % ndim
              call exchange_add(fw % scaling(ii) % lims(1)) 
              call exchange_add(fw % scaling(ii) % lims(2)) 
           enddo
           do ii = 1_ip, fw % main_table % ndim
              call exchange_add(fw % scaling(ii) % kfl_tab) 
           enddo
           call exchange_end()
              
        end select
    
    end subroutine fw_par_exchange


    subroutine tab_load_file(coords, tab, write_table, dont_read_data, order)
      !
      ! Load table from file
      ! 

      type(typ_tab_coord),    pointer, intent(inout)   :: coords(:)
      type(typ_lookup_table), pointer, intent(inout)   :: tab
      logical,            optional,intent(in)      :: write_table
      logical,            optional,intent(in)      :: dont_read_data
      integer(ip),        optional,intent(in)      :: order

      integer(ip)      :: ii, jj, kk, nrow, inds(max_lookup_dim), ioerror
      integer(ip)      :: nchars 
      real(rp)         :: test_output(max_lookup_dim)
      character(len=:), allocatable ::  char_read
      character(30)    :: header_element
      real(rp),pointer :: real_read(:)
      logical          :: read_data, found_break

      !
      ! Allocate pointers
      !
      call tab_allocate(1_ip, coords, tab, order)

      !
      ! Read number of discrete values of each control variable
      !
      do ii=1, max_lookup_dim
         coords(ii) % leng   = 0_ip
      end do

      call ecoute('tab_load_file')
      if (words(1) == 'NVALU') then
         do ii=1, max_lookup_dim
            jj = int(param(ii),ip)
            if (jj > 0_ip) then 
               coords(ii) % leng = jj
            else
               exit
            endif
         end do
         !
         ! Number of dimensions is where it stops
         !
         tab % ndim = ii-1_ip 
      end if


      !
      ! Read number of variables 
      !
      call ecoute('tab_load_file')
      if (words(1) == 'NVARI') then
         tab % nvar = int(param(1_ip),ip) 
      end if

      !
      ! Allocate table
      !
      call tab_allocate(2_ip, coords, tab)


      !
      ! Read line length and coordinates
      !
      call ecoute('tab_load_file')
      if (words(1) == 'NCHUN') then
         line_length = int(param(1_ip),ip) 
         call ecoute('tab_load_file')
      end if


      do ii = 1, tab % ndim
         if ( coords(ii) % leng == 1_ip ) then 
            !
            ! Coordinate length == 1
            !  
            coords(ii) % name   = 'DUMMY'
            coords(ii) % x      = -1.0e12_rp
         else
            !
            ! Read table
            !
            coords(ii) % name = words(1)
            read(nunit, *, iostat=ioerror) coords(ii) % x(1:coords(ii) % leng)

            if ((ioerror == iostat_eor).or.(ioerror == iostat_end)) ioerror = 0
            if (ioerror /= 0 ) call runend('tab_load_file: cannot read table coordinates.')
            call ecoute('tab_load_file')

         endif
      enddo

      !
      ! Assign coordinates to table
      !
      tab % coords => coords(1:tab % ndim)


      !
      ! Count rows in table
      !
      nrow = 1_ip
      do ii = 1, tab % ndim
         nrow = nrow * coords(ii) % leng
      enddo

      !
      ! Output properties
      !
      call messages_live('TABLE','START SECTION')
      !shap = 0_ip 
      !do ii = 1, tab % ndim
      !   shap(ii) = coords(ii) % leng
      !enddo
      call messages_live('tab % ndim=                    ',INT_NUMBER = tab % ndim,        FMT='(a,1x,i5)')
      call messages_live('tab % nvar=                    ',INT_NUMBER = tab % nvar,        FMT='(a,1x,i5)')
      call messages_live('shap=                          ',INT_ARRAY  = tab % shap, FMT='(a,300(1x,i5))')
      call messages_live('nrow=                       ',   INT_NUMBER = nrow,              FMT='(a,1x,i8)')

      call messages_live('size(Mbytes)=                ',   REAL_NUMBER = real(nrow*tab % nvar*kind(tab%tab),rp)/1024.0_rp**2,FMT='(a,1x,f10.2)')

      !
      ! Option for omit reading table
      ! Useful for cases where the shape is known before
      ! but the table is filled runtime 
      !
      read_data = .true.
      if (present(dont_read_data)) then
         if (dont_read_data) then
            read_data = .false.
         endif
      endif

      if (read_data) then
         !
         ! Read header
         !
         nchars = 30_ip * (tab % ndim + tab % nvar)
         allocate(character(len=nchars) :: char_read)
         backspace(nunit)
         read(nunit, "(A)", iostat=ioerror, advance="NO", size=nchars) char_read 
         !
         ! Jump back if new-line chracters were present in this file
         !
         do kk = 1,len(char_read)
            if (iachar(char_read(kk:kk)) == 10) then ! newline
               backspace(nunit)
            endif
         enddo


         kk = 1
         do ii = 1, tab % ndim + tab % nvar
            !
            ! Process header element 
            !
            found_break = .false.
            header_element = ""
            jj = 1
            do while((.not. found_break) .and. (kk < len(char_read)))
                !
                ! Skip underscore
                !
                if (char_read(kk:kk) == "_") then
                    kk = kk + 1
                    continue
                endif

                !
                ! Fill header element and incrememnt counter
                !
                header_element(jj:jj) = char_read(kk:kk)
                kk = kk + 1
                jj = jj + 1 

                !
                ! Find end of word
                !
                if (char_read(kk:kk) == " " .or. iachar(char_read(kk:kk)) == 9) then ! space or tab
                    found_break = .true.
                    do while((char_read(kk:kk) == " " .or. iachar(char_read(kk:kk)) == 9) .and. (kk < len(char_read)))
                       kk = kk + 1
                    enddo
                endif
            enddo

            !
            ! If there was a dummy coordinate, get it's name:
            !
            if (ii <= tab % ndim) then
               if ( coords(ii) % name  == 'DUMMY' ) then  
                  coords(ii) % name  = upper_case(trim(header_element))
               endif
            else
               tab % varname(ii-tab%ndim) = upper_case(trim(header_element))
               tab % longvarname(ii-tab%ndim) = upper_case(trim(header_element))
            endif
         enddo

         deallocate(char_read)

         !
         ! Read table in chunks of ndim + nvar 
         !
         nullify(real_read)
         call memory_alloca(mem_modul(1:2,modul),'REALREAD','tab_load_file',real_read, tab % ndim + tab % nvar)
         do ii = 1, nrow
            read(nunit, *, iostat=ioerror) real_read

            if ((ioerror == iostat_eor).or.(ioerror == iostat_end)) ioerror = 0
            if (ioerror /= 0 ) call runend('tab_load_file: cannot read table body.')

            tab % tab(1:tab %nvar,ii) = real_read(tab % ndim+1 : tab % ndim+tab % nvar) 
         enddo
         call memory_deallo(mem_modul(1:2,modul),'REALREAD','tab_load_file',real_read)
      else
         print'(A)','Not reading table data.'
      endif


      !
      ! Output coordinates and variable names
      !
      do ii = 1, tab % ndim
         call messages_live('tab % coords('//trim(integer_to_string(ii))//'): '//tab % coords(ii) % name)
      enddo
      if (read_data) then               
          do ii = 1,size(tab % varname, KIND=ip),5
             if( ii == 1 ) then
                call messages_live('tab % varname=                 ',CHAR_ARRAY = tab % varname(ii:min(ii+4,size(tab % varname, KIND=ip))),FMT='(a,300(1x,a))')
             else
                call messages_live('                               ',CHAR_ARRAY = tab % varname(ii:min(ii+4,size(tab % varname, KIND=ip))),FMT='(a,300(1x,a))')
             end if
          end do

          !
          ! Long variable names
          !
          do ii = 1,size(tab % varname, KIND=ip)
             if( ii == 1 ) then
                call messages_live('tab % longvarname=             ')
             end if
             call messages_live(trim(integer_to_string(ii))//" "//tab % varname(ii)//": "//tab % longvarname(ii))
          end do
      endif
      call messages_live('TABLE','END SECTION')           

      !
      ! Write table as a check 
      !
      if (present(write_table)) then
         if (write_table) then
            if (read_data) then               
               do ii = 1, nrow
                  call i2inds(ii, tab % ndim ,tab%shap,inds)
                  do jj = 1, tab%ndim
                     test_output(jj) = tab % coords(jj) % x(inds(jj))
                  enddo
                  write(660+tab%ndim,'(500(2X,E20.8))') test_output(1:tab%ndim),  tab % tab(1:tab%nvar,ii)
               enddo
            endif
         endif
      endif

    end subroutine tab_load_file




    subroutine seek_ind_seq(vec,val,leng,i,w)
        !
        ! find index of val in vec
        !
        implicit none
        integer(ip),    intent(in)     :: leng
        real(rp),       intent(in)     :: vec(leng),val
        integer(ip),    intent(inout)  :: i
        real(rp),       intent(out)    :: w(2)

        integer(ip)                    :: ii
        integer(ip)                    :: first, last, dir, i2
        real(rp)                       :: val_clip
       
        w(1)        = 1.0_rp
        w(2)        = 0.0_rp

        if (leng == 1) then
            i       = 1_ip
        else
            dir     = 1_ip
            first   = i
            last    = leng-1_ip
            ! 
            ! loop through array starting form the initial value of ind 
            ! 
            val_clip    = max(min(val,vec(leng)),vec(1_ip))   ! bound 

            !
            ! If clipped value is lower than initial guess, 
            ! and there are lower values in vec, 
            ! then go backwards
            !
            if ( val_clip < vec(i) .and. i > 1_ip ) then 
                dir    = -1_ip
                first  = i
                last   = 1_ip
            endif

            do ii = first,last,dir 
              if (vec(ii + 1_ip ) >= val_clip .and. &
                  vec(ii        ) <= val_clip) then
      
                i    = ii
                i2   = ii + 1_ip
                w(1) = ( vec(i2) - val_clip ) / ( vec(i2) - vec(i) )
                w(2) = 1.0_rp-w(1)
                exit
                 
              end if
            end do 
        endif

    end subroutine seek_ind_seq



    subroutine i2inds(i,dime,shap,inds)
        !
        ! get indeces from line number
        !
        implicit none
        integer(ip), intent(in)  :: i,dime,shap(dime)
        integer(ip), intent(out) :: inds(dime)

        integer(ip)              :: ii, itemp

        itemp = i - 1_ip

        do ii = dime,1_ip,-1_ip
            inds(ii) = mod(itemp,(shap(ii))) + 1_ip
            itemp = itemp / (shap(ii))
        enddo 
    end subroutine i2inds



    subroutine inds2i(inds,dime,shap,i)
        !
        ! get line number from indeces
        !
        implicit none
        integer(ip), intent(in)  :: dime,shap(dime),inds(dime)
        integer(ip), intent(out) :: i
        integer(ip)              :: ii

        i = inds(1_ip)-1_ip

        do ii = 2_ip,dime
            i = i * shap(ii) + (inds(ii)-1_ip)
        enddo
        i = i + 1_ip
    end subroutine inds2i

   subroutine lagrangeInterp1D(np,nvar,xc,yc,x,y)
      !
      ! Arbitrary order interpolation in 1D
      !
      implicit none
      integer(ip), intent(in)  :: np
      integer(ip), intent(in)  :: nvar
      real(rp),    intent(in)  :: xc(:)
      real(rp),    intent(in)  :: yc(:,:)
      real(rp),    intent(in)  :: x
      real(rp),    intent(out) :: y(:)

      real(rp)                 :: L(np)

      integer(ip)              :: ii,jj

      !
      ! Initialize aux variable
      !
      L = 1.0_rp
      do ii = 1,np
         do jj = 1,np
            if (ii /= jj) then
               L(ii) = L(ii) * ( (x - xc(jj))/(xc(ii)-xc(jj)) )
            endif
         enddo
      enddo

      !
      ! Interpolate
      !
      y = 0.0_rp
      do ii = 1,np
         y(:) = y(:) + L(ii) * yc(ii,:)
      enddo

   end subroutine lagrangeInterp1D

   subroutine tab_ho_interp(tab,val,res,ind)
      !
      ! high order interpolation in table 
      !
      !                                                                 
      !      x--x----x----x                                                             
      !     /  /    /    /|                                              
      !    /  /    /    / x         x--x----x----x                      
      !   /  /    /    / /|         |  |    |    |                      
      !  x--x----x----x / |  step1  x--x----x----x  step2                 step3
      !  |  |    |    |/  |  --->   |  |    |    |  --->  x--x----x----x  --->  x
      !  x--x----x----x   x         |  |    |    |                      
      !  |  |    |    |  /          |  |    |    |                      
      !  |  |    |    | /           x--x----x----x                      
      !  |  |    |    |/                                                
      !  x--x----x----x                                                 
      !                                                                 
      implicit none
      type(typ_lookup_table), intent(in)      :: tab
      real(rp),    intent(in)             :: val(:)
      real(rp),    intent(out)            :: res(:)
      integer(ip),optional,intent(inout)  :: ind(:)

      integer(ip)                         :: inbour(2_ip),                  &
                                             inear(tab%ndim,tab%order+1_ip), &
                                             ii, jj, kk, itab(tab%ndim),&
                                             nres, ord(tab%ndim), nst, shast(tab%ndim),&
                                             maxjj
      real(rp)                            :: hyper( (tab%order+1_ip)**tab%ndim, tab%nvar )
      real(rp)                            :: hyper_prev( (tab%order+1_ip)**tab%ndim, tab%nvar )
      real(rp)                            :: wnbour(2_ip) 
      
      real(rp)                            :: xc(tab%order+1_ip) 
      real(rp)                            :: yc(tab%order+1_ip,tab%nvar) 

      logical                             :: ord_parity 
                                             
      
      res         = 0.0_rp
      inear       = 0_ip
      ord         = 0_ip
      nst         = 1_ip
      shast       = 1_ip
      nres        = min( size(res, KIND=ip), tab%nvar )

      !
      ! Order parity
      !
      ord_parity = .false.
      if (mod(tab%order,2_ip)==0_ip) ord_parity = .true.

      !
      ! find neighbouring indeces along dimensions
      !
      do ii = 1,tab%ndim
         if (present(ind)) then
             inbour(1)  = ind(ii)
             inbour(2)  = ind(ii)+1_ip
         else
             inbour(1)  = 1_ip 
             inbour(2)  = 2_ip 
         endif

         call seek_ind_seq( tab % coords(ii) % x, &
              val(ii), &
              tab % coords(ii) % leng, &
              inbour(1), &
              wnbour(:))
         inbour(2) = inbour(1) + 1_ip
         !
         ! Determine real order in each dimension, and choose indeces
         ! 
         if ((inbour(1) == inbour(2))) then
            !
            ! Dummy dimension, or outside of table
            !
            ord(ii) = 0_ip
            inear(ii,1) = inbour(1)
         elseif(tab%order == 0_ip) then
            !
            ! 0 order = Nearest neigbour
            !
            ord(ii) = 0_ip
            if (wnbour(1) > 0.5_rp) then
               inear(ii,1) = inbour(1)
            else
               inear(ii,1) = inbour(2)
            endif
         elseif (abs(wnbour(1)-1.0_rp)<1e-6_rp) then
            !
            ! Near left data point, use nearest neigbour
            !
            inear(ii,1) = inbour(1)
            inear(ii,2) = inbour(1)
            ord(ii)     = 0_ip
         elseif (abs(wnbour(2)-1.0_rp)<1e-6_rp) then
            !
            ! Near right data point, use nearest neigbour
            !
            inear(ii,1) = inbour(2)
            inear(ii,2) = inbour(2)
            ord(ii)     = 0_ip
         else
            !
            ! Add neighbours:
            ! 
            inear(ii,1) = inbour(1)
            inear(ii,2) = inbour(2)
            ord(ii)     = 1_ip

            !
            ! Number of additional steps in each direction:
            ! int behaves like floor
            ! o=1 -> 0
            ! o=2 -> 0
            ! o=3 -> 1
            ! o=4 -> 1 ...
            !
            maxjj = int(real(tab%order-1_ip,rp)/2.0_rp,ip)
            do jj = 1,maxjj 
               !
               ! Increase order and add index to list
               !
               if (inbour(1)-jj > 0_ip) then
                  ord(ii) = ord(ii) + 1_ip
                  inear(ii,ord(ii)+1_ip) = inbour(1)-jj
               endif

               if (inbour(2)+jj < tab % coords(ii) % leng ) then
                  ord(ii) = ord(ii) + 1_ip
                  inear(ii,ord(ii)+1_ip) = inbour(2)+jj
               endif
            enddo

            !
            ! If order is even, make asymmetric stencil
            !
            if (ord_parity) then
               if (wnbour(1) > 0.5_rp) then
                  !
                  ! Try left side
                  !
                  if ( inbour(1)-maxjj > 1_ip ) then
                     ord(ii) = ord(ii) + 1_ip
                     inear(ii,ord(ii)+1_ip) = inbour(1)-(maxjj+1)
                  elseif( inbour(2)+maxjj < tab % coords(ii) % leng-1 ) then
                     !
                     ! Resort to right side
                     !
                     ord(ii) = ord(ii) + 1_ip
                     inear(ii,ord(ii)+1_ip) = inbour(2)+(maxjj+1)
                  endif
               else
                  !
                  ! Try right side
                  !
                  if( inbour(2)+maxjj < tab % coords(ii) % leng-1 ) then        
                     ord(ii) = ord(ii) + 1_ip
                     inear(ii,ord(ii)+1_ip) = inbour(2)+(maxjj+1)
                  elseif ( inbour(1)-maxjj > 1_ip ) then
                     !
                     ! Resort to left side
                     !
                     ord(ii) = ord(ii) + 1_ip
                     inear(ii,ord(ii)+1_ip) = inbour(1)-(maxjj+1)
                  endif
               endif
            endif
         endif
         nst        = nst * (ord(ii) + 1_ip)
         shast(ii)  = ord(ii) + 1_ip
      end do

      
      !
      ! Look up all points
      !
      hyper(1:nst,:) = 0.0_rp
      do jj = 1,nst ! any combination of the individual indeces 
         !
         ! Get index in inear
         !
         call i2inds(jj, tab % ndim ,shast, itab)
         do ii = 1,tab % ndim
            itab(ii) = inear(ii,itab(ii))
         enddo

         !
         ! Get index in table
         !
         call inds2i(itab,tab%ndim,tab%shap,kk)
         hyper(jj,:) = tab%tab(:,kk)
      enddo


      !
      ! recursive interpolation along all dimensions:
      !
      do ii = 1,tab%ndim
         hyper_prev(1:nst,:) = hyper(1:nst,:)
         !
         ! New stencil size:
         !
         nst = 1_ip
         do jj = (ii+1),tab%ndim
            nst = nst * (ord(jj)+1_ip)
         enddo
         
         !
         ! Coordinates in current dimension:
         !
         do kk = 1,(ord(ii)+1_ip)
            xc(kk) = tab % coords(ii) % x(inear(ii,kk))
         enddo

         !
         ! Interpolate on new stencil:
         !
         hyper = 0.0_rp
         do jj = 1,nst 
            do kk = 0,(ord(ii))
               yc(kk+1,:) = hyper_prev(jj+kk*nst,:)
            enddo
            call lagrangeInterp1D(ord(ii)+1_ip,tab%nvar,xc,yc,val(ii),hyper(jj,:))
         enddo 

      enddo

      res(1:nres) = hyper(1,1:nres)
         
      if (present(ind)) then
         ind(1:tab%ndim) = inear(1:tab%ndim, 1)
      endif

   end subroutine tab_ho_interp


    subroutine tab_interp(tab,val,res,ind,snap_on_nearest,irow)
      !
      ! interpolate in table 
      !
      implicit none
      type(typ_lookup_table), intent(in)  :: tab
      real(rp),    intent(in)             :: val(:)
      real(rp),    intent(out)            :: res(:)
      integer(ip), intent(inout)          :: ind(:)
      logical,    optional,intent(in)     :: snap_on_nearest
      integer(ip),optional,intent(out)    :: irow

      integer(ip)                         :: ii, itab(tab%ndim), &
                                             nres, nvertex, irows(2_ip**tab%ndim)
      real(rp)                            :: wnear(tab%ndim,2), hyper(tab%nvar,2_ip**tab%ndim),&
                                             weight(2_ip**tab%ndim)
      
      logical                             :: snap

      snap = .false.
      if ( present(snap_on_nearest) ) then
          if (snap_on_nearest) then
              snap = .true.
          endif
      endif


      if ( (tab%order /= 1) .and. (.not. snap)) then !
         !
         ! Do high order interpolation
         !
         call tab_ho_interp(tab,val,res,ind)
      else
         !
         ! Do original method
         !

         if (.not. initialized) then
            call tab_interp_initialize()
         endif

         nres                 = min( size(res, KIND=ip), tab%nvar )
         nvertex              = 2_ip**tab%ndim
         weight(1:nvertex)    = 1.0_rp


         !
         ! find indeces along dimensions, and compute weights
         !
         do ii = 1,tab%ndim
            call seek_ind_seq( tab % coords(ii) % x, &
                 val(ii), &
                 tab % coords(ii) % leng, &
                 ind(ii), &
                 wnear(ii,:) )
            weight(1:nvertex) = weight(1:nvertex) * wnear(ii,precomp_i2inds(ii,1:nvertex,tab%ndim))
         end do

         !
         ! loop through vertices of hypercube
         !
         do ii = 1,nvertex
            !
            ! find index along coordinates
            !
            itab(1:tab%ndim) = ind(1:tab%ndim) + (precomp_i2inds(1:tab%ndim,ii,tab%ndim)-1_ip)

            !
            ! find index in table 
            !
            call inds2i(itab,tab%ndim,tab%shap,irows(ii))
         end do

         !
         ! get values on vertices
         !
         hyper(1:tab%nvar,1:nvertex) = tab%tab(1:tab%nvar,irows(1:nvertex))
         if (snap) then
             !
             ! Result is the closest vertex, a.k.a. the vertex with the highest
             ! weight
             !
             ii          = maxloc(weight, DIM=1)
             irow        = irows(ii) 
             res(1:nres) = hyper(1:nres,ii)
         else
             !
             ! result is weighted average of vertex values
             !
             res(1:nres) = matmul(hyper(1:nres,:),weight(:))
         endif
             
      endif

    end subroutine tab_interp

    subroutine tab_interpt(tab,val,res,ind,snap_on_nearest,irow)
      !
      ! interpolate in table 
      !
      implicit none
      type(typ_lookup_table), intent(in)  :: tab
      real(rp),    intent(in)             :: val(:)
      real(rp),    intent(out)            :: res(:)
      integer(ip), intent(inout)          :: ind(:)
      logical,    optional,intent(in)     :: snap_on_nearest
      integer(ip),optional,intent(out)    :: irow

      integer(ip)                         :: ii, itab(tab%ndim), &
                                             nres, nvertex, irows(2_ip**tab%ndim)
      real(rp)                            :: wnear(tab%ndim,2), hyper(tab%nvar,2_ip**tab%ndim),&
                                             weight(2_ip**tab%ndim)
      
      logical                             :: snap
      real(rp)                            :: time_ini, time_end

      snap = .false.
      if ( present(snap_on_nearest) ) then
          if (snap_on_nearest) then
              snap = .true.
          endif
      endif

      
      if ( (tab%order /= 1) .and. (.not. snap)) then !
         !
         ! Do high order interpolation
         !
         call tab_ho_interp(tab,val,res,ind)
      else
         !
         ! Do original method
         !
         call cputim(time_ini)

         if (.not. initialized) then
            call tab_interp_initialize()
         endif

         nres                 = min( size(res, KIND=ip), tab%nvar )
         nvertex              = 2_ip**tab%ndim
         weight(1:nvertex)    = 1.0_rp
        
         call cputim(time_end)
         time_int1 = time_int1 + (time_end - time_ini)
         call cputim(time_ini)   

         !
         ! find indeces along dimensions, and compute weights
         !
         do ii = 1,tab%ndim
            call seek_ind_seq( tab % coords(ii) % x, &
                 val(ii), &
                 tab % coords(ii) % leng, &
                 ind(ii), &
                 wnear(ii,:) )
            weight(1:nvertex) = weight(1:nvertex) * wnear(ii,precomp_i2inds(ii,1:nvertex,tab%ndim))
         end do

         call cputim(time_end)
         time_int2 = time_int2 + (time_end - time_ini)
         call cputim(time_ini)

         !
         ! loop through vertices of hypercube
         !
         do ii = 1,nvertex
            !
            ! find index along coordinates
            !
            itab(1:tab%ndim) = ind(1:tab%ndim) + (precomp_i2inds(1:tab%ndim,ii,tab%ndim)-1_ip)

            !
            ! find index in table 
            !
            call inds2i(itab,tab%ndim,tab%shap,irows(ii))
         end do

         !
         ! get values on vertices
         !
         hyper(1:tab%nvar,1:nvertex) = tab%tab(1:tab%nvar,irows(1:nvertex))


         call cputim(time_end)
         time_int3 = time_int3 + (time_end - time_ini)
         call cputim(time_ini)


         if (snap) then
             !
             ! Result is the closest vertex, a.k.a. the vertex with the highest
             ! weight
             !
             ii          = maxloc(weight, DIM=1)
             irow        = irows(ii) 
             res(1:nres) = hyper(1:nres,ii)
         else
             !
             ! result is weighted average of vertex values
             !
             res(1:nres) = matmul(hyper(1:nres,:),weight(:))
         endif
             
         call cputim(time_end)
         time_int4 = time_int4 + (time_end - time_ini)
         call cputim(time_ini)

      endif

    end subroutine tab_interpt


    subroutine print_tab_interp_times()
        real(rp) :: tloc(9)
        real(rp) :: perc(7)

        tloc(1)   = time_int1
        tloc(2)   = time_int2
        tloc(3)   = time_int3
        tloc(4)   = time_int4
        tloc(5)   = time_int5
        tloc(6)   = time_look1
        tloc(7)   = time_look2
        tloc(8)   = sum(tloc(1:5))
        tloc(9)   = sum(tloc(6:7))
        perc(1:5) = tloc(1:5) / tloc(8)
        perc(6:7) = tloc(6:7) / tloc(9)
        call PAR_MAX(9_ip,tloc)
        call PAR_MAX(7_ip,perc)
        

        
        if ( INOTSLAVE ) then 
           !print'(A21,1X,F18.14,2X,F6.3)', 'Do original method', tloc(1),tloc(1)/tloc(8)*100.0_rp
           !print'(A21,1X,F18.14,2X,F6.3)', 'Find indices      ', tloc(2),tloc(2)/tloc(8)*100.0_rp
           !print'(A21,1X,F18.14,2X,F6.3)', 'Loop hypercube    ', tloc(3),tloc(3)/tloc(8)*100.0_rp
           !print'(A21,1X,F18.14,2X,F6.3)', 'If snap or not    ', tloc(4),tloc(4)/tloc(8)*100.0_rp
           !print'(A21,1X,F18.14,2X,F6.3)', 'If ind or not     ', tloc(5),tloc(5)/tloc(8)*100.0_rp
           !print'(A21,1X,F18.14)',         'Total time Interp ', tloc(8)

           !print'(A21,1X,F18.14)',         'Time outside scaling', tloc(6)
           !print'(A21,1X,F18.14)',         'Time outside   intep', tloc(7)

           print'(A21,1X,F18.14,2X,F6.3,"%")', 'Use ind if present  ', tloc(1),perc(1)*100.0_rp
           print'(A21,1X,F18.14,2X,F6.3,"%")', 'Find indices        ', tloc(2),perc(2)*100.0_rp
           print'(A21,1X,F18.14,2X,F6.3,"%")', 'Loop hypercube      ', tloc(3),perc(3)*100.0_rp
           print'(A21,1X,F18.14,2X,F6.3,"%")', 'Do weighted average ', tloc(4),perc(4)*100.0_rp
           print'(A21,1X,F18.14,2X,F6.3,"%")', 'Out ind if present  ', tloc(5),perc(5)*100.0_rp
           print'(A21,1X,F18.14)',         'Total time Interp   ', tloc(8)

           print'(A21,1X,F18.14,2X,F6.3,"%")', 'Time outside scaling', tloc(6),perc(6)*100.0_rp
           print'(A21,1X,F18.14,2X,F6.3,"%")', 'Time outside   intep', tloc(7),perc(7)*100.0_rp
        endif

    end subroutine print_tab_interp_times

    !
    ! Clip between 0 and 1
    !
    real(rp) function clip_cont_var(y,ymin,ymax)
        real(rp),  intent(in)           :: y       
        real(rp),  intent(in), optional :: ymin       
        real(rp),  intent(in), optional :: ymax

        real(rp) :: ymin_loc, ymax_loc

        ymin_loc = 0.0_rp
        ymax_loc = 1.0_rp
        if (present(ymin)) ymin_loc = ymin
        if (present(ymax)) ymax_loc = ymax

        clip_cont_var = min(ymax_loc,max(ymin_loc,(y)))
    end function clip_cont_var




    !
    ! Scale control variables with framework
    !
    subroutine fw_scale_cont_var( control, scale_control, lim_control, fw, ind)
        real(rp),                         intent(in)    :: control(:)
        real(rp),                         intent(out)   :: scale_control(:)
        real(rp),                         intent(out)   :: lim_control(:,:)
        type(typ_lookup_framework),       intent(in)    :: fw
        integer(ip),                      intent(inout) :: ind(:)

        integer(ip)               :: ii, jj, kk, jmean, imain 
        real(rp)                  :: xtab_loc(fw % main_table % ndim)
        real(rp)                  :: xtemp(fw % main_table % ndim)
        integer(ip)               :: indtemp(fw % main_table % ndim) 
        real(rp)                  :: res_tab(5,fw % main_table % ndim)
        real(rp)                  :: var_max, denom, numer, max2, minmax, min2, var 

        scale_control = 0.0_rp
        lim_control   = 0.0_rp
        xtab_loc      = 0.0_rp
        res_tab       = 0.0_rp

        do kk = 1, fw % main_table % ndim
            ii = fw % scale_order(kk)
            if ( fw % kfl_scale(ii) >= 0_ip .and. fw % kfl_scale(ii) <= 100_ip ) then
                 !
                 ! Simple scalings
                 !
                 select case( fw % kfl_scale(ii) )
                     case (0_ip)
                         !
                         ! No scaling required: e.g.: zmean
                         !
                         scale_control(ii) = control(ii)
                         lim_control(ii,1) = fw % main_table % coords(ii) % x(1)
                         lim_control(ii,2) = fw % main_table % coords(ii) % x(fw % main_table % coords(ii) % leng)
                     case (1_ip)
                         !
                         ! Use a table for scaling: e.g.: cmean, imean, chist
                         !
                         do jj = 1, fw % scaling(ii) % tab % ndim
                             imain = fw % scaling(ii) % i_main(jj)
                             xtemp(jj) = xtab_loc(imain)
                             indtemp(jj) = ind(imain)
                         enddo
                         call tab_interp( fw % scaling(ii) % tab, &
                                        & xtemp(1:fw % scaling(ii) % tab % ndim), &
                                        & res_tab(1:fw % scaling(ii) % tab % nvar,ii), &
                                        & indtemp(1:fw % scaling(ii) % tab % ndim))

                         do jj = 1, fw % scaling(ii) % tab % ndim
                             imain = fw % scaling(ii) % i_main(jj)
                             ind(imain) = indtemp(jj) 
                         enddo

                         scale_control(ii) = (control(ii) - res_tab(2,ii)) / max(( res_tab(1,ii) - res_tab(2,ii) ), 1.0e-10_rp)
                         lim_control(ii,1) = res_tab(2,ii)
                         lim_control(ii,2) = res_tab(1,ii)

                     case (2_ip)
                         !
                         ! Use simple limits for scaling: e.g.: imean for premixed
                         !
                         scale_control(ii) = ( control(ii) - fw % scaling(ii) % lims(1) ) /  max(( fw % scaling(ii) % lims(2) - fw % scaling(ii) % lims(1) ), 1.0e-10_rp)
                         lim_control(ii,:) = fw % scaling(ii) % lims(:)
                     case (3_ip)
                         !
                         ! Force constant scaled control variable: e.g.: imean=1.0 if we want to use a table with heat loss in the
                         ! adiabatic limit
                         !
                         scale_control(ii) = fw % scaling(ii) % lims(1)
                         lim_control(ii,1) = -1.0e12_rp
                         lim_control(ii,2) = 1.0e12_rp
                 end select
            elseif ( abs(fw % kfl_scale(ii)) >= 100_ip ) then
                !
                ! Variance or square scaling
                !
                jmean   = abs(fw % kfl_scale(ii)) - 100
                var_max = clip_cont_var(scale_control(jmean)) * (1.0_rp  - clip_cont_var(scale_control(jmean)))

                select case( fw % kfl_scale(jmean) )
                    case (0_ip) 
                        !
                        ! The mean did not require scaling: e.g.: zvar
                        !  
                        if (fw % kfl_scale(ii) < 0) then
                            !
                            ! Square
                            !
                            lim_control(ii,1) = clip_cont_var(scale_control(jmean)) * clip_cont_var(scale_control(jmean))
                            lim_control(ii,2) = 1.0_rp 
                            if (var_max /= 0.0_rp) then
                            !if (var_max /= 0.0_rp .and. control(jmean) > 0.0_rp .and. (control(ii) - lim_control(ii,1)) > 1.0e-10_rp) then
                                scale_control(ii) = (control(ii) - lim_control(ii,1)) / var_max
                            else
                                scale_control(ii) = 0.0_rp
                                !scale_control(ii) = (control(ii) - lim_control(ii,1))
                            endif
                        else
                            !
                            ! Variance
                            !
                            lim_control(ii,1) = 0.0_rp 
                            lim_control(ii,2) = var_max
                            if (var_max /= 0.0_rp) then
                                scale_control(ii) = control(ii) / lim_control(ii,2)
                            else
                                scale_control(ii) = 0.0_rp
                            endif
                        endif
                    case (1_ip,2_ip)
                        !
                        ! The mean required a tabulated scaling 
                        !
                        ! 
                        ! y        : variable
                        ! yv       : varaince of variable 
                        ! s        : scaled variable
                        ! sv       : variance of scaled variable
                        ! \tilde{} : averaging
                        ! min      : minimum of vriable
                        ! max      : maximum of vriable
                        ! min2     : square of minimum of vriable averaged
                        ! max2     : square of maximum of vriable averaged
                        ! minmax   : minimum times maximum of vriable averaged
                        !
                        ! sv = \tilde{s s} - \tilde{s} \tilde{s}
                        !
                        ! \tilde{s s} = \tilde{(y-min)/(max-min) (y-min)/(max-min)}
                        !
                        !                 \tilde{(y-min)^2}       numer   
                        ! \tilde{s s} = --------------------- = ----------
                        !                \tilde{(max-min)^2}      denom     
                        !
                        ! denom = \tilde{(max-min)^2} = max2 - 2minmax + min2
                        ! 
                        ! numer = \tilde{y^2} - \tilde{2 y min - min min}
                        !
                        ! \tilde{y^2} = \tilde{y}^2 + yv
                        ! 
                        ! - \tilde{2 y min - min min} = - \tilde{2 ((max-min)s + min) min - min min}
                        !                             = - \tilde{2 min max s - 2 min min s + 2 min min - min min}
                        !                             = - \tilde{2 min max s - 2 min min s + min min}
                        !                             = -min2 -2\tilde{s} ( minmax - min2 ) 
                        ! 
                        ! 
                        ! numer =  \tilde{y^2} -min2 -2\tilde{s} ( minmax - min2 ) 
                        ! OR
                        ! numer =  \tilde{y}^2 + yv -min2 -2\tilde{s} ( minmax - min2 ) 
                        !
                        !        numer   
                        ! sv = ---------- - \tilde{s} \tilde{s}
                        !        denom   
                        !
                        if ( fw % kfl_scale(jmean) == 1 ) then
                           max2   = res_tab(3,jmean)
                           min2   = res_tab(4,jmean)
                           minmax = res_tab(5,jmean)
                        else
                           max2   = fw % scaling(ii) % lims(2)**2
                           min2   = fw % scaling(ii) % lims(1)**2
                           minmax = fw % scaling(ii) % lims(1) * fw % scaling(ii) % lims(2)
                        endif
                        denom  = max(max2 - 2.0_rp * minmax + min2,1.0e-10_rp)
                        if (fw % kfl_scale(ii) < 0) then
                            !
                            ! Square
                            !
                            numer = control(ii)
                        else
                            !
                            ! Variance
                            !
                            numer = control(ii) +  control(jmean) * control(jmean)
                        endif
                        
                        numer = numer - min2
                        numer = numer - 2.0_rp * clip_cont_var(scale_control(jmean)) * ( minmax - min2 ) 

                        !
                        ! Scaled control varaible
                        !
                        var = numer / denom - clip_cont_var(scale_control(jmean))**2
                        if (var_max /= 0.0_rp) then
                            scale_control(ii) = var / var_max
                        else
                            scale_control(ii) = 0.0_rp
                        endif

                        !
                        ! Calculate limits
                        !
                        !   0 < sv < \tilde{s}(1-\tilde{s})
                        !
                        !   \tilde{s}^2 denom < numer < \tilde{s} denom
                        !  
                        !   \tilde{s}^2 denom + min2 + 2 \tilde{s} ( minmax - min2 ) < \tilde{y^2} < \tilde{s} denom + min2 + 2 \tilde{s} ( minmax - min2 )
                        !   OR
                        !   \tilde{s}^2 denom + min2 + 2 \tilde{s} ( minmax - min2 ) - \tilde{y}^2 < yv < \tilde{s} denom + min2 + 2 \tilde{s} ( minmax - min2 ) - \tilde{y}^2
                        !
                        lim_control(ii,1) = clip_cont_var(scale_control(jmean))**2 * denom + min2 +  2.0_rp * clip_cont_var(scale_control(jmean)) * ( minmax - min2 )
                        lim_control(ii,2) = clip_cont_var(scale_control(jmean))    * denom + min2 +  2.0_rp * clip_cont_var(scale_control(jmean)) * ( minmax - min2 )
                        if (fw % kfl_scale(ii) > 0) then
                            lim_control(ii,1) = lim_control(ii,1) - clip_cont_var( control(jmean), ymin= lim_control(jmean,1),ymax=lim_control(jmean,2))**2
                            lim_control(ii,2) = lim_control(ii,2) - clip_cont_var( control(jmean), ymin= lim_control(jmean,1),ymax=lim_control(jmean,2))**2
                        endif

                end select

            endif

            !
            ! Fill xtab 
            !
            xtab_loc(ii) = clip_cont_var( scale_control(ii),                                        &
                                      ymin= fw % main_table % coords(ii) % x(1),                    &
                                      ymax= fw % main_table % coords(ii) % x( fw % main_table % coords(ii) % leng))            

        enddo

    end subroutine fw_scale_cont_var



    subroutine fw_lookup( control, scale_control, fw, retva, ind )
        real(rp),                         intent(in)    :: control(:)
        real(rp),                         intent(out)   :: scale_control(:)
        type(typ_lookup_framework),       intent(in)    :: fw
        real(rp),                         intent(out)   :: retva(:)
        integer(ip),                      intent(inout) :: ind(:)

        real(rp)    :: lim_control(fw % main_table % ndim,2)
        !!AB real(rp)    :: time_ini, time_end
        !
        ! Sale control variables
        !
        !!AB call cputim(time_ini)
        call fw_scale_cont_var( control, scale_control, lim_control, fw, ind )
        !!AB call cputim(time_end)
        !!AB time_look1 = time_look1 + (time_end -time_ini)

        !!AB call cputim(time_ini)
        !
        ! Interpolate properties from database
        !
        call tab_interp(fw % main_table,scale_control,retva, ind) 

        !!AB call cputim(time_end)
        !!AB time_look2 = time_look2 + (time_end -time_ini)


    end subroutine fw_lookup

end module mod_interp_tab
