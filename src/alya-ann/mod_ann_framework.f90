!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_ann_framework 

  use def_kintyp_basic,   only : ip,rp,lg
  use mod_ann_scaling,    only : ann_scaling
  use mod_ann_scaling,    only : ANN_SCALING_UNSCALED
  use mod_ann_scaling,    only : ANN_SCALING_LINEAR

  implicit none
  private

  !
  ! Constant parameters
  !
  integer(ip),   parameter :: ANN_BACKEND_TORCH  = 1
  !integer(ip),   parameter :: ANN_BACKEND_???  = 2

  integer(ip),   parameter :: max_ann_stage         = 4
  integer(ip),   parameter :: ANN_STAGE_READING     = 1
  integer(ip),   parameter :: ANN_STAGE_INPUTSCALE  = 2
  integer(ip),   parameter :: ANN_STAGE_FORWARDPASS = 3
  integer(ip),   parameter :: ANN_STAGE_OUTPUTSCALE = 4
 
  integer(ip),   parameter :: ANN_PRECISION_SINGLE  = 1
  integer(ip),   parameter :: ANN_PRECISION_DOUBLE  = 2

  !
  ! Timing
  !
  real(rp)                 :: time1,time2

  !
  ! Class to manage artificial neural networks
  !
  type ann_framework
      !
      ! General attributes
      !
      integer(ip)                               :: index                 ! Index of this framework
      integer(ip)                               :: backend               ! Flag indicating the library  
      integer(ip)                               :: precision             ! Flag indicating the precision set for the model 
      character(500)                            :: path                  ! Path of the saved ann file 
      logical(lg)                               :: is_read               ! Whether the path is processed already
      real(rp)                                  :: times(max_ann_stage)  ! Timings

      !
      ! Attributes related to input and output layers
      !
      type(ann_scaling)                         :: scaling_in            ! Scaling parameters of input
      type(ann_scaling)                         :: scaling_out           ! Scaling parameters of output

    contains
      procedure,                 pass           :: init                  ! Initialize with default values 
      procedure,                 pass           :: set                   ! Set properties 
      procedure,                 pass           :: set_model_precision   ! Set floating point precision of model
      procedure,                 pass           :: par_comm              ! Parallel communication of attributes 
      procedure,                 pass           :: read_file             ! Read artificial neural network file 
      procedure,                 pass           :: print                 ! Print data to standard output 
      procedure,                 pass           :: tim_ini               ! Init time counters
      procedure,                 pass           :: tim_end               ! End time counters
      procedure,                 pass           :: tim_statistics        ! Calculate parallel statistics of timings 

      procedure,                 pass           :: eval_r1               ! Evaluate based on 1D array  
      procedure,                 pass           :: eval_r2               ! Evaluate based on 2D array (First dimension is the vector size, and second is the flattened input and output shape) 
      generic                                   :: eval   =>  &          ! Generic evaluation subroutine
           &                                       eval_r1,   &
           &                                       eval_r2     
  end type ann_framework

  !
  ! Make specific procedures and variables public
  !
  public :: ann_framework

  public :: ANN_BACKEND_TORCH

  public :: max_ann_stage
  public :: ANN_STAGE_READING
  public :: ANN_STAGE_INPUTSCALE
  public :: ANN_STAGE_FORWARDPASS
  public :: ANN_STAGE_OUTPUTSCALE

  public :: ANN_PRECISION_SINGLE
  public :: ANN_PRECISION_DOUBLE

contains

  !
  ! INITIALIZATION 
  !
  subroutine init(self)

    class(ann_framework), intent(inout) :: self
    
    self % index     = 0
    self % backend   = 0
    self % precision = ANN_PRECISION_SINGLE
    self % path      = ""
    self % is_read   = .false.
    self % times     = 0.0 

    call self % scaling_in % init()
    call self % scaling_out % init()

  end subroutine init

  !
  ! SET VALUES 
  !
  subroutine set(self, index, backend, path, precision)

    class(ann_framework), intent(inout)     :: self
    integer(4),                    optional :: index
    integer(ip),                   optional :: backend
    character(len=:), allocatable, optional :: path
    integer(ip),                   optional :: precision

    if (present(index))   self % index   = index
    if (present(backend)) self % backend = backend
    call self % set_model_precision(precision)
    if (present(path))    self % path    = path

  end subroutine set

  !
  ! PARALLEL COMMUNICATON 
  !
  subroutine par_comm(self)
    use def_master,    only : parii, npari, nparr,     &
                              nparc, parin, parre,     &
                              parch, nparc, mem_modul, & 
                              modul, ISLAVE, IMASTER
    use mod_communications, only : PAR_BROADCAST
    use mod_communications, only : PAR_EXCHANGE
    use mod_memory_basic,   only : memory_alloca
    use mod_memory_basic,   only : memory_deallo

    class(ann_framework), intent(inout)     :: self


    !
    ! Exchange constants
    !  
    do parii=1,2 
       npari=0
       nparr=0
       nparc=0

       call PAR_EXCHANGE(self % index,     parin,npari,parii)           
       call PAR_EXCHANGE(self % backend,   parin,npari,parii)           
       call PAR_EXCHANGE(self % precision, parin,npari,parii)           

       if( parii == 1 ) then
          call memory_alloca(mem_modul(1:2,modul),'PARIN','ann_framework % par_comm',parin,max(1_ip,npari))
          call memory_alloca(mem_modul(1:2,modul),'PARRE','ann_framework % par_comm',parre,max(1_ip,nparr))

          if( ISLAVE  ) call PAR_BROADCAST(parin,      'IN MY CODE')
          if( ISLAVE  ) call PAR_BROADCAST(parre,      'IN MY CODE')
          if( ISLAVE  ) call PAR_BROADCAST(nparc,parch,'IN MY CODE')
       else
          if( IMASTER ) call PAR_BROADCAST(parin,      'IN MY CODE')
          if( IMASTER ) call PAR_BROADCAST(parre,      'IN MY CODE')
          if( IMASTER ) call PAR_BROADCAST(nparc,parch,'IN MY CODE')
       end if
    enddo
    call memory_deallo(mem_modul(1:2,modul),'PARIN','ann_framework % par_comm',parin)
    call memory_deallo(mem_modul(1:2,modul),'PARRE','ann_framework % par_comm',parre)

    !
    ! EXIT if index is zero
    !
    if ( self % index == 0_ip) return

    !
    ! Exchange path
    !
    do parii=1,2 
       npari=0
       nparr=0
       nparc=0

       call PAR_EXCHANGE(self % path,    parch,nparc,parii)           

       if( parii == 1 ) then
          call memory_alloca(mem_modul(1:2,modul),'PARIN','ann_framework % par_comm',parin,max(1_ip,npari))
          call memory_alloca(mem_modul(1:2,modul),'PARRE','ann_framework % par_comm',parre,max(1_ip,nparr))

          if( ISLAVE  ) call PAR_BROADCAST(parin,      'IN MY CODE')
          if( ISLAVE  ) call PAR_BROADCAST(parre,      'IN MY CODE')
          if( ISLAVE  ) call PAR_BROADCAST(nparc,parch,'IN MY CODE')
       else
          if( IMASTER ) call PAR_BROADCAST(parin,      'IN MY CODE')
          if( IMASTER ) call PAR_BROADCAST(parre,      'IN MY CODE')
          if( IMASTER ) call PAR_BROADCAST(nparc,parch,'IN MY CODE')
       end if
    enddo
    call memory_deallo(mem_modul(1:2,modul),'PARIN','ann_framework % par_comm',parin)
    call memory_deallo(mem_modul(1:2,modul),'PARRE','ann_framework % par_comm',parre)

    !
    ! EXCHANGE SCALING DATA
    !
    call self % scaling_in % par_comm()
    call self % scaling_out % par_comm()

  end subroutine par_comm

  !
  ! READ NETWORK FILE
  !
  subroutine read_file(self)

    class(ann_framework), intent(inout)     :: self
    
    !
    ! Exit if index is zero
    !
    if ( self % index == 0_ip) return

    !
    ! Start timing
    !
    call self % tim_ini()
    
    !
    ! Read network
    !
    if (.not. self % is_read) then     
       if (self % backend == ANN_BACKEND_TORCH ) then     
#ifdef TORCH
          call read_neural_network(self % index, self % path, len(trim(self % path)))
#endif
       endif
    else
        print*, "Not reading ANN file: "//trim(self % path)//" because it was already read."
    endif
    self % is_read = .true.

    !
    ! Set precision of model 
    !
    call self % set_model_precision()
    
    !
    ! End timing
    !
    call self % tim_end(ANN_STAGE_READING)

  end subroutine read_file

  subroutine set_model_precision(self, precision)
    class(ann_framework), intent(inout)        :: self
    integer(ip),          intent(in), optional :: precision

    !
    ! Save internal variable indicating precision
    !
    if (present(precision)) then
       self % precision = precision
    endif

    !
    ! Set model precision
    !
    if (self % is_read) then
       select case( self % precision )
       case (ANN_PRECISION_SINGLE)
          if (self % backend == ANN_BACKEND_TORCH ) then     
#ifdef TORCH
             call set_neural_network_precision_to_single(self % index)
#endif
          endif
       case (ANN_PRECISION_DOUBLE)
          if (self % backend == ANN_BACKEND_TORCH ) then     
#ifdef TORCH
             call set_neural_network_precision_to_double(self % index)
#endif
          endif
       end select
    endif

  end subroutine set_model_precision

  !
  ! PRINT TO STANDARD OUTPUT
  !
  subroutine print(self)
    use mod_messages,       only : messages_live
    use mod_strings,        only : integer_to_string
    use def_master,         only : retost

    class(ann_framework), intent(in)  :: self
    integer(ip)                       :: ii

    if (self % backend == ANN_BACKEND_TORCH) then
       call messages_live('ann % backend =                '//'Torch')
    endif

    if (self % precision == ANN_PRECISION_SINGLE) then
       call messages_live('ann % precision =              '//'Single')
    elseif (self % precision == ANN_PRECISION_DOUBLE) then
        call messages_live('ann % precision =              '//'Double')
    endif
    
    call messages_live('ann % path =                   '//trim(self % path))

    !
    ! Print input scalign data
    !
    call messages_live('SCALIN IN','START SECTION')
    call messages_live('n_dim =                     '//trim(integer_to_string(self % scaling_in % n_dim)))
    call messages_live('shape =                     ',INT_ARRAY  = self % scaling_in % shape, FMT='(a,300(1x,i5))'  )
    call messages_live('n_prod_shape =              '//trim(integer_to_string(self % scaling_in % n_prod_shape)))
    do ii = 1, self % scaling_in % n_prod_shape 
       if (self % scaling_in % names(ii) /= "     ") &
          call messages_live('    '//trim(integer_to_string(ii))//' name ='//trim(self % scaling_in % names(ii)))
    enddo
    select case(self % scaling_in % all_types)
    case (ANN_SCALING_UNSCALED)
        call messages_live('All elelments are unscaled')
    case (ANN_SCALING_LINEAR)
        call messages_live('All elelments are scaled lineraly')
        do ii = 1, self % scaling_in % n_prod_shape 
           call messages_live('    '//trim(integer_to_string(ii))//' coeff='//retost(self % scaling_in % coeffs(ii))//' shift='//retost(self % scaling_in % shifts(ii)) )
        enddo
    end select
    call messages_live('SCALIN IN','END SECTION')

    !
    ! Print output scalign data
    !
    call messages_live('SCALIN OUT','START SECTION')
    call messages_live('n_dim =                     '//trim(integer_to_string(self % scaling_out % n_dim)))
    call messages_live('shape =                     ',INT_ARRAY  = self % scaling_out % shape, FMT='(a,300(1x,i5))'  )
    call messages_live('n_prod_shape =              '//trim(integer_to_string(self % scaling_out % n_prod_shape)))
    do ii = 1, self % scaling_out % n_prod_shape 
       if (self % scaling_out % names(ii) /= "     ") &
          call messages_live('    '//trim(integer_to_string(ii))//' name ='//trim(self % scaling_out % names(ii)))
    enddo
    select case(self % scaling_out % all_types)
    case (ANN_SCALING_UNSCALED)
        call messages_live('All elelments are unscaled')
    case (ANN_SCALING_LINEAR)
        call messages_live('All elelments are scaled lineraly')
        do ii = 1, self % scaling_out % n_prod_shape 
           call messages_live('    '//trim(integer_to_string(ii))//' coeff='//retost(self % scaling_out % coeffs(ii))//' shift='//retost(self % scaling_out % shifts(ii)) )
        enddo
    end select
    call messages_live('SCALIN OUT','END SECTION')
    
  end subroutine print


  !
  ! PARALLEL STATISTICS OF TIMINGS
  !
  subroutine tim_statistics(self, maximum, average, total)
    use mod_communications,    only : PAR_AVERAGE
    use mod_communications,    only : PAR_MAX
    use mod_communications,    only : PAR_SUM

    class(ann_framework), intent(in)  :: self
    real(rp), intent(out)             :: maximum(max_ann_stage) 
    real(rp), intent(out)             :: average(max_ann_stage) 
    real(rp), intent(out)             :: total(max_ann_stage)
        

    maximum = 0.0_rp
    average = 0.0_rp
    total   = 0.0_rp

    !
    ! Calculate parallel statistics
    !
    maximum = self % times(1:max_ann_stage)
    average = maximum
    total   = maximum

    call PAR_MAX    (max_ann_stage, maximum)
    call PAR_AVERAGE(max_ann_stage, average)
    call PAR_SUM    (max_ann_stage, total)

  end subroutine tim_statistics

  !
  ! Start timer
  !
  subroutine tim_ini(self)

    class(ann_framework), intent(in)  :: self
    
    call cputim(time1)

  end subroutine tim_ini

  !
  ! Accumulate time in specific bin
  !
  subroutine tim_end(self,i_stage)

    class(ann_framework), intent(inout) :: self
    integer(ip),          intent(in)    :: i_stage
    
    call cputim(time2)
    self % times(i_stage) = self % times(i_stage) + (time2 - time1)
    
  end subroutine tim_end


  !
  ! EVALUATING FUNCTIONS 
  !
  subroutine eval_r1(self, input, output)

    class(ann_framework), intent(inout)     :: self
    real(rp),   intent(inout)               :: input(self % scaling_in % n_prod_shape)
    real(rp),   intent(out)                 :: output(self % scaling_out % n_prod_shape)
    integer(4)                              :: shape_input(100)
    real(4)                                 :: output_loc4(self % scaling_out % n_prod_shape)
    real(8)                                 :: output_loc8(self % scaling_out % n_prod_shape)

    output = 0.0_rp

    !
    ! Scale inputs
    !
    call self % tim_ini()
    if (self % scaling_in % all_types /= ANN_SCALING_UNSCALED) then 
        call self % scaling_in % scale( input, input )
    endif
    call self % tim_end(ANN_STAGE_INPUTSCALE)

    !
    ! Evaluate ANN
    !
    call self % tim_ini()
    shape_input(1:self % scaling_in % n_dim) = self % scaling_in % shape 
    if (self % backend == ANN_BACKEND_TORCH ) then     
#ifdef TORCH
        !
        ! Use Torch to do a forward pass of the network
        ! 
        if (self % precision == ANN_PRECISION_SINGLE) then
           call forward_pass_neural_network(int(self % index, kind=4),                &
               &                            int(self % precision, kind=4),            &
               &                            real(input,kind=4),                       &
               &                            output_loc4,                              &
               &                            int(self % scaling_in % n_dim, kind=4),   &
               &                            shape_input(1:self % scaling_in % n_dim), & 
               &                            int(self % scaling_out % n_dim, kind=4))

        elseif (self % precision == ANN_PRECISION_DOUBLE) then
           call forward_pass_neural_network(int(self % index, kind=4),                &
               &                            int(self % precision, kind=4),            &
               &                            real(input,kind=8),                       &
               &                            output_loc8,                              &
               &                            int(self % scaling_in % n_dim, kind=4),   &
               &                            shape_input(1:self % scaling_in % n_dim), & 
               &                            int(self % scaling_out % n_dim, kind=4))
        endif
#endif
    endif
    call self % tim_end(ANN_STAGE_FORWARDPASS)

    !
    ! Unscale outputs
    !
    call self % tim_ini()
    if (self % precision == ANN_PRECISION_SINGLE) then
        output = real(output_loc4, kind=rp)
    elseif (self % precision == ANN_PRECISION_DOUBLE) then
        output = real(output_loc8, kind=rp)
    endif
    if (self % scaling_out % all_types /= ANN_SCALING_UNSCALED) then 
        call self % scaling_out % unscale( output, output )
    endif
    call self % tim_end(ANN_STAGE_OUTPUTSCALE)

  end subroutine eval_r1

  subroutine eval_r2(self, nvec, input, output )

    class(ann_framework), intent(inout)     :: self
    integer(ip),intent(in)                  :: nvec 
    real(rp),   intent(inout)               :: input(nvec, self % scaling_in % n_prod_shape)
    real(rp),   intent(out)                 :: output(nvec, self % scaling_out % n_prod_shape)
    integer(4)                              :: shape_input(100)
    real(4)                                 :: output_loc4(nvec, self % scaling_out % n_prod_shape)
    real(8)                                 :: output_loc8(nvec, self % scaling_out % n_prod_shape)

    output = 0.0_rp

    !
    ! Scale inputs
    !
    call self % tim_ini()
    if (self % scaling_in % all_types /= ANN_SCALING_UNSCALED) then 
        call self % scaling_in % scale( nvec, input, input )
    endif
    call self % tim_end(ANN_STAGE_INPUTSCALE)

    !
    ! Evaluate ANN
    !
    call self % tim_ini()
    shape_input(1) = nvec
    shape_input(2:(self % scaling_in % n_dim+1)) = self % scaling_in % shape 
    if (self % backend == ANN_BACKEND_TORCH ) then     
#ifdef TORCH
        !
        ! Use Torch to do a forward pass of the network
        ! 
        if (self % precision == ANN_PRECISION_SINGLE) then
           call forward_pass_neural_network(int(self % index, kind=4),                    &
               &                            int(self % precision, kind=4),                &
               &                            real(input,kind=4),                           &
               &                            output_loc4,                                  &
               &                            int(self % scaling_in % n_dim+1, kind=4),     &
               &                            shape_input(1:(self % scaling_in % n_dim+1)), & 
               &                            int(self % scaling_out % n_dim+1, kind=4))
        elseif (self % precision == ANN_PRECISION_DOUBLE) then
           call forward_pass_neural_network(int(self % index, kind=4),                    &
               &                            int(self % precision, kind=4),                &
               &                            real(input,kind=8),                           &
               &                            output_loc8,                                  &
               &                            int(self % scaling_in % n_dim+1, kind=4),     &
               &                            shape_input(1:(self % scaling_in % n_dim+1)), & 
               &                            int(self % scaling_out % n_dim+1, kind=4))
        endif
#endif
    endif
    call self % tim_end(ANN_STAGE_FORWARDPASS)

    !
    ! Unscale outputs
    !
    call self % tim_ini()
    if (self % precision == ANN_PRECISION_SINGLE) then
        output = real(output_loc4, kind=rp)
    elseif (self % precision == ANN_PRECISION_DOUBLE) then
        output = real(output_loc8, kind=rp)
    endif
    if (self % scaling_out % all_types /= ANN_SCALING_UNSCALED) then 
        call self % scaling_out % unscale( nvec, output, output )
    endif
    call self % tim_end(ANN_STAGE_OUTPUTSCALE)

  end subroutine eval_r2



end module mod_ann_framework
