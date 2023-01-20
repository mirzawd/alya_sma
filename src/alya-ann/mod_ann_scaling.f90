!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_ann_scaling

  use def_kintyp_basic,   only : ip,rp

  implicit none
  private

  !
  ! Constant parameters
  !
  integer(ip),   parameter :: ANN_SCALING_UNSCALED = 1
  integer(ip),   parameter :: ANN_SCALING_LINEAR   = 2

  !
  ! Class to hold scaling parameters and execute the scaling on arbitrary inputs and outputs 
  !
  type ann_scaling
      !
      ! Attributes related to input and output layers
      !
      integer(ip)                               :: n_dim                 ! Number of dimensions for layer
      integer(ip), pointer                      :: shape(:)              ! Shape 
      integer(ip)                               :: n_prod_shape          ! Product of shape 
      integer(ip), pointer                      :: types(:)              ! Types of scaling 
      integer(ip)                               :: all_types             ! Set if all scaling types coincide 
      real(rp), pointer                         :: coeffs(:)             ! Coefficients for linear scaling 
      real(rp), pointer                         :: shifts(:)             ! Shifts for linear scaling
      character(5), pointer                     :: names(:)              ! Name of each variable

      integer(ip), pointer                      :: iaux1(:)              ! Auxiliary integer variable 

    contains
      procedure,                 pass           :: init                  ! Initialize with default values 
      procedure,                 pass           :: set                   ! Set properties 
      procedure,                 pass           :: set_elem              ! Set element 
      procedure,                 pass           :: alloc_shape           ! Allocate shape 
      procedure,                 pass           :: alloc                 ! Allocate other variables depending on the shape 
      procedure,                 pass           :: par_comm              ! Parallel communication of attributes 

      procedure,                 pass           :: calc_prod_shape       ! Auxiliary function to calculate product of shape 
      procedure,                 pass           :: calc_all_types        ! Auxiliary function to check if all types coincide and set self%all_types accordingly
      procedure,                 pass           :: inds2i                ! Auxiliary function to calculate flattened index from indeces

      procedure,                 pass           :: scale_r1               ! Scale 1D array  
      procedure,                 pass           :: scale_r2               ! Scale 2D array 
      generic                                   :: scale   =>  &          ! Generic scaling subroutine
           &                                       scale_r1,   &
           &                                       scale_r2     
      procedure,                 pass           :: unscale_r1             ! Unscale 1D array  
      procedure,                 pass           :: unscale_r2             ! Unscale 2D array 
      generic                                   :: unscale   =>  &        ! Generic Unscaling subroutine
           &                                       unscale_r1,   &
           &                                       unscale_r2     
  end type ann_scaling

  !
  ! Make specific procedures public
  !
  public :: ann_scaling

  public :: ANN_SCALING_UNSCALED
  public :: ANN_SCALING_LINEAR

contains

  !
  ! INITIALIZATION 
  !
  subroutine init(self)

    class(ann_scaling), intent(inout) :: self
    
    self % n_dim          = 0
    nullify(self % shape)
    self % n_prod_shape   = 0
    nullify(self % types)
    self % all_types      = 0
    nullify(self % coeffs)
    nullify(self % shifts)
    nullify(self % names)
    nullify(self % iaux1)

  end subroutine init


  !
  ! SET VALUES 
  !
  subroutine set(self, n_dim, i_shape, n_shape, inds, typ, coeff, shift, name)

    class(ann_scaling), intent(inout)       :: self
    integer(ip),  intent(in),      optional :: n_dim
    integer(ip),  intent(in),      optional :: i_shape
    integer(ip),  intent(in),      optional :: n_shape
    integer(ip),  intent(in),      optional :: inds(self % n_dim)
    integer(ip),  intent(in),      optional :: typ
    real(rp),     intent(in),      optional :: coeff
    real(rp),     intent(in),      optional :: shift
    character(5), intent(in),      optional :: name 
    
    if (present(n_dim))   self % n_dim   = n_dim

    !
    ! Set shape: i_shape^th dimension is set to n_shape
    !
    if (present(i_shape) .and. present(n_shape)) then
        if (associated(self % shape)) then
            if (self % n_dim >= i_shape) then
                self % shape(i_shape) = n_shape
                call self % calc_prod_shape()
            endif
        endif
    endif

    !
    ! Set element defined by inds
    !
    call self % set_elem(inds, typ, coeff, shift, name)

  end subroutine set

  subroutine set_elem(self, inds, typ, coeff, shift, name)

    class(ann_scaling), intent(inout)       :: self
    integer(ip),  intent(in),      optional :: inds(self % n_dim)
    integer(ip),  intent(in),      optional :: typ
    real(rp),     intent(in),      optional :: coeff
    real(rp),     intent(in),      optional :: shift
    character(5), intent(in),      optional :: name 
    integer(ip)                             :: ii  
    
    !
    ! Set element defined by inds
    !
    if (present(inds)) then
        ii = self % inds2i(inds)
        !
        ! Set type and check if all types coincide
        !
        if (present(typ)) then
            self % types(ii) = typ 
            call self % calc_all_types()
        endif

        !
        ! Set coeff
        !
        if (present(coeff)) then
            self % coeffs(ii) = coeff
        endif

        !
        ! Set shift
        !
        if (present(shift)) then
            self % shifts(ii) = shift
        endif

        !
        ! Set name  
        !
        if (present(name)) then
            self % names(ii) = name
        endif
    endif

  end subroutine set_elem


  !
  ! ALLOCATE VARIABLES 
  !
  subroutine alloc_shape(self)

    class(ann_scaling), intent(inout)       :: self
        
    if (self % n_dim > 0 .and. (.not. associated(self % shape))) then
        allocate(self % shape(self % n_dim))    
        self % shape = 0
    endif

  end subroutine alloc_shape


  subroutine alloc(self)

    class(ann_scaling), intent(inout)       :: self
    integer(ip)                             :: ii  
        
    if (self % n_prod_shape > 0) then 
        allocate(self % types(self % n_prod_shape))    
        self % types = 0
        allocate(self % coeffs(self % n_prod_shape))    
        self % coeffs = 0.0
        allocate(self % shifts(self % n_prod_shape))    
        self % shifts = 0.0
        allocate(self % names(self % n_prod_shape))    
        do ii = 1,self % n_prod_shape
           self % names = "     "
        enddo
        allocate(self % iaux1(self % n_prod_shape))    
        self % iaux1 = 0
    endif 

  end subroutine alloc

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

    class(ann_scaling), intent(inout)     :: self
    integer(ip)                           :: ii   


    !
    ! Exchange constants
    !  
    do parii=1,2 
       npari=0
       nparr=0
       nparc=0

       call PAR_EXCHANGE(self % n_dim,   parin,npari,parii)           

       if( parii == 1 ) then
          call memory_alloca(mem_modul(1:2,modul),'PARIN','ann_scaling % par_comm',parin,max(1_ip,npari))
          call memory_alloca(mem_modul(1:2,modul),'PARRE','ann_scaling % par_comm',parre,max(1_ip,nparr))

          if( ISLAVE  ) call PAR_BROADCAST(parin,      'IN MY CODE')
          if( ISLAVE  ) call PAR_BROADCAST(parre,      'IN MY CODE')
          if( ISLAVE  ) call PAR_BROADCAST(nparc,parch,'IN MY CODE')
       else
          if( IMASTER ) call PAR_BROADCAST(parin,      'IN MY CODE')
          if( IMASTER ) call PAR_BROADCAST(parre,      'IN MY CODE')
          if( IMASTER ) call PAR_BROADCAST(nparc,parch,'IN MY CODE')
       end if
    enddo
    call memory_deallo(mem_modul(1:2,modul),'PARIN','ann_scaling % par_comm',parin)
    call memory_deallo(mem_modul(1:2,modul),'PARRE','ann_scaling % par_comm',parre)

    !
    ! EXIT if index is zero
    !
    if ( self % n_dim == 0_ip) return

    !
    ! Allocate shape 
    !
    call self % alloc_shape()

    !
    ! Exchange shape
    !
    do parii=1,2 
       npari=0
       nparr=0
       nparc=0

       call PAR_EXCHANGE(self % n_prod_shape,    parin,npari,parii)           
       call PAR_EXCHANGE(self % all_types,    parin,npari,parii)           
       do ii = 1, self % n_dim
           call PAR_EXCHANGE(self % shape(ii),    parin,npari,parii)           
       enddo

       if( parii == 1 ) then
          call memory_alloca(mem_modul(1:2,modul),'PARIN','ann_scaling % par_comm',parin,max(1_ip,npari))
          call memory_alloca(mem_modul(1:2,modul),'PARRE','ann_scaling % par_comm',parre,max(1_ip,nparr))

          if( ISLAVE  ) call PAR_BROADCAST(parin,      'IN MY CODE')
          if( ISLAVE  ) call PAR_BROADCAST(parre,      'IN MY CODE')
          if( ISLAVE  ) call PAR_BROADCAST(nparc,parch,'IN MY CODE')
       else
          if( IMASTER ) call PAR_BROADCAST(parin,      'IN MY CODE')
          if( IMASTER ) call PAR_BROADCAST(parre,      'IN MY CODE')
          if( IMASTER ) call PAR_BROADCAST(nparc,parch,'IN MY CODE')
       end if
    enddo
    call memory_deallo(mem_modul(1:2,modul),'PARIN','ann_scaling % par_comm',parin)
    call memory_deallo(mem_modul(1:2,modul),'PARRE','ann_scaling % par_comm',parre)

    !
    ! Allocate memory for slaves
    !
    if( ISLAVE )  call self % alloc()

    !
    ! Exchange data
    !
    do parii=1,2 
       npari=0
       nparr=0
       nparc=0

       do ii = 1, self % n_prod_shape
          call PAR_EXCHANGE(self % types(ii),    parin,npari,parii)           
          call PAR_EXCHANGE(self % coeffs(ii),   parre,nparr,parii)           
          call PAR_EXCHANGE(self % shifts(ii),   parre,nparr,parii)           
          call PAR_EXCHANGE(self % names(ii),    parch,nparc,parii)           
          call PAR_EXCHANGE(self % iaux1(ii),    parin,npari,parii)           
       enddo

       if( parii == 1 ) then
          call memory_alloca(mem_modul(1:2,modul),'PARIN','ann_scaling % par_comm',parin,max(1_ip,npari))
          call memory_alloca(mem_modul(1:2,modul),'PARRE','ann_scaling % par_comm',parre,max(1_ip,nparr))

          if( ISLAVE  ) call PAR_BROADCAST(parin,      'IN MY CODE')
          if( ISLAVE  ) call PAR_BROADCAST(parre,      'IN MY CODE')
          if( ISLAVE  ) call PAR_BROADCAST(nparc,parch,'IN MY CODE')
       else
          if( IMASTER ) call PAR_BROADCAST(parin,      'IN MY CODE')
          if( IMASTER ) call PAR_BROADCAST(parre,      'IN MY CODE')
          if( IMASTER ) call PAR_BROADCAST(nparc,parch,'IN MY CODE')
       end if
    enddo
    call memory_deallo(mem_modul(1:2,modul),'PARIN','ann_scaling % par_comm',parin)
    call memory_deallo(mem_modul(1:2,modul),'PARRE','ann_scaling % par_comm',parre)


  end subroutine par_comm


  !
  ! Calculate product of shape
  !
  subroutine calc_prod_shape(self)

    class(ann_scaling), intent(inout)     :: self
    integer(ip)                           :: ii   
        
    if (self % n_dim > 0) then 
       self % n_prod_shape = 1
       do ii = 1,self % n_dim
          self % n_prod_shape = self % n_prod_shape * self % shape(ii) 
       enddo 
    endif 

  end subroutine calc_prod_shape

  !
  ! Check if all scaling types coincide
  !
  subroutine calc_all_types(self)

    class(ann_scaling), intent(inout)     :: self
    integer(ip)                           :: ii   
        
    if (self % n_prod_shape > 0) then 
       !
       ! Loop through all types and check if all of them are the same
       !
       self % all_types = self % types(1)
       do ii = 2,self % n_prod_shape
          if (self % all_types /= self % types(ii)) then
              !
              ! Set to zero and quit if difference is found
              !
              self % all_types = 0
              return
          endif
       enddo 
    endif 

  end subroutine calc_all_types

  !
  ! Calculate flattened index from indeces
  !
  integer(ip) function inds2i(self, inds)
    class(ann_scaling), intent(inout)       :: self
    integer(ip),  intent(in),      optional :: inds(self % n_dim)
    integer(ip)                             :: ii   

    inds2i = inds(self % n_dim)-1

    do ii = self % n_dim-1, 1, -1
        inds2i = inds2i * self % shape(ii) + (inds(ii)-1)
    enddo
    inds2i = inds2i + 1

  end function inds2i


  !
  ! SCALING FUNCTIONS 
  !
  subroutine scale_r1(self, unscaled, scaled )

    class(ann_scaling), intent(in)        :: self
    real(rp),   intent(in)                :: unscaled(self % n_prod_shape)
    real(rp),   intent(out)               :: scaled(self % n_prod_shape)

    if (self % all_types == ANN_SCALING_UNSCALED) then 
        !
        ! Do not scale
        !
        scaled = unscaled

    elseif (self % all_types == ANN_SCALING_LINEAR) then 
        !
        ! Scale linearly
        !
        scaled = unscaled * self % coeffs + self % shifts
    endif  

  end subroutine scale_r1

  subroutine scale_r2(self, nvec, unscaled, scaled )

    class(ann_scaling), intent(in)        :: self
    integer(ip),intent(in)                :: nvec 
    real(rp),   intent(in)                :: unscaled(nvec, self % n_prod_shape)
    real(rp),   intent(out)               :: scaled(nvec, self % n_prod_shape)
    integer(ip)                           :: ii   

    if (self % all_types == ANN_SCALING_UNSCALED) then 
        !
        ! Do not scale
        !
        scaled = unscaled

    elseif (self % all_types == ANN_SCALING_LINEAR) then 
        !
        ! Scale linearly
        !
        do ii = 1,nvec
           scaled(ii,:) = unscaled(ii,:) * self % coeffs + self % shifts
        enddo
    endif  

  end subroutine scale_r2

  !
  ! UNSCALING FUNCTIONS 
  !
  subroutine unscale_r1(self, scaled, unscaled )

    class(ann_scaling), intent(in)        :: self
    real(rp),   intent(in)                :: scaled(self % n_prod_shape)
    real(rp),   intent(out)               :: unscaled(self % n_prod_shape)

    if (self % all_types == ANN_SCALING_UNSCALED) then 
        !
        ! Do not unscale
        !
        unscaled = scaled

    elseif (self % all_types == ANN_SCALING_LINEAR) then 
        !
        ! Unscale linearly
        !
        unscaled = (scaled - self % shifts) / self % coeffs
    endif  

  end subroutine unscale_r1

  subroutine unscale_r2(self, nvec, scaled, unscaled )

    class(ann_scaling), intent(in)        :: self
    integer(ip),intent(in)                :: nvec 
    real(rp),   intent(in)                :: scaled(nvec, self % n_prod_shape)
    real(rp),   intent(out)               :: unscaled(nvec, self % n_prod_shape)
    integer(ip)                           :: ii   

    if (self % all_types == ANN_SCALING_UNSCALED) then 
        !
        ! Do not unscale
        !
        unscaled = scaled

    elseif (self % all_types == ANN_SCALING_LINEAR) then 
        !
        ! Unscale linearly
        !
        do ii = 1,nvec
           unscaled(ii,:) = (scaled(ii,:) - self % shifts) / self % coeffs 
        enddo
    endif  

  end subroutine unscale_r2


end module mod_ann_scaling
