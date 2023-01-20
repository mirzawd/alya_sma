!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!>
!> @defgroup Memory_Toolbox
!> @{
!> @name    ToolBox for memory management
!> @file    mod_memory_physics.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for memory management
!> @details Tools for memory of variables defined in def_kintyp_physics
!>
!------------------------------------------------------------------------

module mod_memory_physics

  use def_kintyp,         only : ip,rp
  use def_kintyp_physics, only : rbtyp
  
  use mod_memory_tools
  use mod_memory_basic

  implicit none

  private

  interface memory_alloca
     module procedure &
          &           memory_alloca_rbtyp,memory_alloca_rbtyp_s
  end interface memory_alloca

  interface memory_deallo
     module procedure &
          &           memory_deallo_rbtyp
  end interface memory_deallo

  public :: memory_alloca
  public :: memory_deallo

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-10
  !> @brief   Rigid body
  !> @details Allotate rigid body type
  !> 
  !-----------------------------------------------------------------------

  subroutine memory_alloca_rbtyp(memor,vanam,vacal,varia,ndim1)

    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    type(rbtyp),  intent(inout), pointer  :: varia(:)      !< Variable
    integer(ip),  intent(in)              :: ndim1         !< Dimension
    integer(4)                            :: istat
    integer(ip)                           :: idim1,ii

    if( ndim1 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)

       allocate( varia(ndim1) , stat = istat )

       if( istat == 0 ) then
          lbytm = storage_size(varia)*ndim1/8
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type')

       do idim1 = 1,ndim1

          varia(idim1) % npoib  =  0         ! Number of boundary nodes
          varia(idim1) % nboib  =  0         ! Number of boundary elements
          varia(idim1) % nrbse  =  0         ! Number of sets the rigid body is formed by
          varia(idim1) % lrbse  =  0         ! List of sets the rigid body is formed by

          varia(idim1) % massa  = -1.0_rp    ! Mass
          varia(idim1) % densi  = -1.0_rp    ! Density
          varia(idim1) % volum  =  0.0_rp    ! Volume
          varia(idim1) % momin  = -1.0_rp    ! Momentum of inertia tensor
          varia(idim1) % posgr  = -1.0e12_rp ! Initial position of center of gravity

          varia(idim1) % posil  = -1.0e12_rp ! Linear position
          varia(idim1) % velol  =  0.0_rp    ! Linear velocity
          varia(idim1) % accel  =  0.0_rp    ! Linear accelaration
          varia(idim1) % force  =  0.0_rp    ! Force
          varia(idim1) % vpfor  =  0.0_rp    ! Viscous + pressure Force
          varia(idim1) % pforce =  0.0_rp    ! Pressure Force
          varia(idim1) % vforce =  0.0_rp    ! Viscous Force
          varia(idim1) % pp_pf  =  0.0_rp    ! Pressure Force for post process
          varia(idim1) % pp_vf  =  0.0_rp    ! Viscous Force for post process

          varia(idim1) % posia  =  0.0_rp    ! Angular position
          varia(idim1) % veloa  =  0.0_rp    ! Angular velocity
          varia(idim1) % accea  =  0.0_rp    ! Angular accelaration
          varia(idim1) % rotac  =  0.0_rp    ! Rotation matrix
          varia(idim1) % torqu  =  0.0_rp    ! Torque
          varia(idim1) % vptor  =  0.0_rp    ! Viscous + pressure Torque
          varia(idim1) % ptorqu =  0.0_rp    ! Pressure Torque
          varia(idim1) % vtorqu =  0.0_rp    ! Viscous Torque
          varia(idim1) % pp_pt  =  0.0_rp    ! Pressure Torque for post process
          varia(idim1) % pp_vt  =  0.0_rp    ! Viscous Torque for post process
          varia(idim1) % quate  =  0.0_rp    ! quaternion (an alternative representation of the rotation matrix)
          varia(idim1) % q_dot  =  0.0_rp    ! time derivative of the quaternion

          varia(idim1) % quate(1,:) = 1.0_rp
          do ii =  1,3
             varia(idim1) % rotac(ii,ii) = 1.0_rp
          end do
          
          nullify(varia(idim1) % cooin)      ! Initial coordinates
          nullify(varia(idim1) % cooib)      ! Coordinates
          nullify(varia(idim1) % lnoib)      ! Boundaries
          nullify(varia(idim1) % ltyib)      ! Boundary types
          nullify(varia(idim1) % lninv)      ! Permutation for nodes
          nullify(varia(idim1) % lbinv)      ! Permutation for boundaries

       end do

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_rbtyp

  subroutine memory_alloca_rbtyp_s(memor,vanam,vacal,varia,ndim1,ndim2,ndim3,ndim4)

    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    type(rbtyp),  intent(inout)           :: varia         !< Variable
    integer(ip),  intent(in)              :: ndim1         !< Dimension 1 
    integer(ip),  intent(in)              :: ndim2         !< Dimension 2
    integer(ip),  intent(in)              :: ndim3         !< Dimension 3
    integer(ip),  intent(in)              :: ndim4         !< Dimension 4

    call memory_alloca(memor,trim(vanam)//' % COOIN',vacal,varia % cooin,ndim1,ndim2)
    call memory_alloca(memor,trim(vanam)//' % COOIB',vacal,varia % cooib,ndim1,ndim2)
    call memory_alloca(memor,trim(vanam)//' % LNOIB',vacal,varia % lnoib,ndim4,ndim3)
    call memory_alloca(memor,trim(vanam)//' % LTYIB',vacal,varia % ltyib,ndim3)
    call memory_alloca(memor,trim(vanam)//' % LNINV',vacal,varia % lninv,ndim2)
    call memory_alloca(memor,trim(vanam)//' % LBINV',vacal,varia % lbinv,ndim3)

  end subroutine memory_alloca_rbtyp_s
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-10
  !> @brief   Rigid body
  !> @details Allotate rigid body type
  !> 
  !-----------------------------------------------------------------------

  subroutine memory_deallo_rbtyp(memor,vanam,vacal,varia)

    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    type(rbtyp),  intent(inout), pointer  :: varia(:)      !< Variable
    integer(4)                            :: istat
    integer(ip)                           :: idim1
    
    if( associated(varia) ) then

       do idim1 = 1,size(varia,KIND=ip)
          call memory_deallo(memor,trim(vanam)//' % COOIN',vacal,varia(idim1) % cooin)
          call memory_deallo(memor,trim(vanam)//' % COOIB',vacal,varia(idim1) % cooib)
          call memory_deallo(memor,trim(vanam)//' % LNOIB',vacal,varia(idim1) % lnoib)
          call memory_deallo(memor,trim(vanam)//' % LTYIB',vacal,varia(idim1) % ltyib)
          call memory_deallo(memor,trim(vanam)//' % LNINV',vacal,varia(idim1) % lninv)
          call memory_deallo(memor,trim(vanam)//' % LBINV',vacal,varia(idim1) % lbinv)
       end do
       
       lbytm = -storage_size(varia)/8*size(varia)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type')

    else

       lbytm = 0
      
    end if

  end subroutine memory_deallo_rbtyp
    
end module mod_memory_physics
!> @}
