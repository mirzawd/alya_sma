!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!>
!> @addtogroup Exmedi
!> @{
!> @file    mod_exm_memory.f90
!> @author  mixed
!> @date    2020-09-16
!> @brief   Memory management for exmedi types
!> @details 
!> @}
!>
!------------------------------------------------------------------------

module mod_exm_memory

   use def_kintyp,           only : ip,rp
   use mod_exm_cellmodel,    only : CELL_MODEL_OUTPUTS
   use mod_memory,           only : lbytm, kfl_alloc
   use def_exmedi,           only : ecg_coords
   use mod_exm_oharaprecalc, only : ohara_precalc

   use mod_memory_tools
   use mod_memory_basic

   implicit none

   private

   interface memory_alloca
      module procedure &
            &           memory_alloca_land, & ! OHARA_LAND_PARAMETERS
            &           memory_alloca_ecg, & ! ecg_coords
            &           memory_alloca_ohara_precalc
   end interface memory_alloca

   interface memory_deallo
      module procedure &
            &           memory_deallo_land, & ! OHARA_LAND_PARAMETERS
            &           memory_deallo_ecg, & ! ecg_coords
            &           memory_deallo_ohara_precalc
   end interface memory_deallo

   public :: memory_alloca
   public :: memory_deallo

contains




subroutine memory_alloca_ohara_precalc(memor,vanam,vacal,varia,ndim1)
    implicit none
    character(*), intent(in)                              :: vanam         !< Variable name
    character(*), intent(in)                              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)                           :: memor(2)      !< Memory counter
    type(ohara_precalc),  intent(inout), pointer          :: varia(:)      !< Variable
    integer(ip),  intent(in)                              :: ndim1         !< Dimension
    integer(4)                                            :: istat

    if( ndim1 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)

       allocate( varia(ndim1) , stat = istat )

       if( istat == 0 ) then
          lbytm = storage_size(varia)*ndim1/8
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type')

    else

       nullify(varia)

    end if

 end subroutine memory_alloca_ohara_precalc



 subroutine memory_deallo_ohara_precalc(memor,vanam,vacal,varia)
    implicit none
    character(*), intent(in)                              :: vanam         !< Variable name
    character(*), intent(in)                              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)                           :: memor(2)      !< Memory counter
    type(ohara_precalc),  intent(inout), pointer          :: varia(:)      !< Variable
    integer(4)                                            :: istat
    
    if( associated(varia) ) then
       
       lbytm = -storage_size(varia)/8*size(varia,1,kind=ip)
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

 end subroutine memory_deallo_ohara_precalc


   !-----------------------------------------------------------------------
   !> 
   !> @author  Constantine Butakoff
   !> @date    2020-09-16
   !> @brief   Allotate OHARA_LAND_PARAMETERS
   !> @details Allotate OHARA_LAND_PARAMETERS
   !> 
   !-----------------------------------------------------------------------

   subroutine memory_alloca_land(memor,vanam,vacal,varia,ndim1)
      implicit none
      character(*), intent(in)                              :: vanam         !< Variable name
      character(*), intent(in)                              :: vacal         !< Calling subroutine name
      integer(8),   intent(inout)                           :: memor(2)      !< Memory counter
      type(CELL_MODEL_OUTPUTS),  intent(inout), pointer          :: varia(:)      !< Variable
      integer(ip),  intent(in)                              :: ndim1         !< Dimension
      integer(4)                                            :: istat
      integer(ip)                                           :: idim1

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
            varia(idim1) %  S =         0.0_rp
            varia(idim1) %  W =         0.0_rp
            varia(idim1) %  CaTRPN =    0.0_rp
            varia(idim1) %  B =         0.0_rp
            varia(idim1) %  zeta_s =    0.0_rp
            varia(idim1) %  zeta_w =    0.0_rp
            varia(idim1) %  Ca50 =      0.805_rp
            varia(idim1) %  Lambda =    1.0_rp
            varia(idim1) %  Vinit  =      0.0_rp     ! initial voltage
            varia(idim1) %  nbeats =     -1_ip       ! number of taken to converge
            varia(idim1) %  toler =      -1.0_rp     ! tolerance to decide if the cell model converged
            varia(idim1) %  rmse =       -1.0_rp     ! rmse between the last two beats of the cell model use to decide convergence
         end do

      else

         nullify(varia)

      end if

   end subroutine memory_alloca_land


   !-----------------------------------------------------------------------
   !> 
   !> @author  Constantine Butakoff
   !> @date    2020-09-16
   !> @brief   Allotate OHARA_LAND_PARAMETERS
   !> @details Allotate OHARA_LAND_PARAMETERS
   !> 
   !-----------------------------------------------------------------------

   subroutine memory_deallo_land(memor,vanam,vacal,varia)
      implicit none
      character(*), intent(in)                              :: vanam         !< Variable name
      character(*), intent(in)                              :: vacal         !< Calling subroutine name
      integer(8),   intent(inout)                           :: memor(2)      !< Memory counter
      type(CELL_MODEL_OUTPUTS),  intent(inout), pointer          :: varia(:)      !< Variable
      integer(4)                                            :: istat
      
      if( associated(varia) ) then
         
         lbytm = -storage_size(varia)/8*size(varia,1,kind=ip)
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

   end subroutine memory_deallo_land

   

  !-----------------------------------------------------------------------
   !> 
   !> @author  Constantine Butakoff
   !> @date    2020-09-16
   !> @brief   Allotate ecg_coords
   !> @details Allotate ecg_coords
   !> 
   !-----------------------------------------------------------------------

   subroutine memory_alloca_ecg(memor,vanam,vacal,varia,ndim1)
      implicit none
      character(*), intent(in)                              :: vanam         !< Variable name
      character(*), intent(in)                              :: vacal         !< Calling subroutine name
      integer(8),   intent(inout)                           :: memor(2)      !< Memory counter
      type(ecg_coords),  intent(inout), pointer             :: varia(:)      !< Variable
      integer(ip),  intent(in)                              :: ndim1         !< Dimension
      integer(4)                                            :: istat
      integer(ip)                                           :: idim1

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
            varia(idim1) % coords = 0.0_rp
            varia(idim1) % label   = "     "  
         end do

      else

         nullify(varia)

      end if

   end subroutine memory_alloca_ecg

   subroutine memory_deallo_ecg(memor,vanam,vacal,varia)
      implicit none
      character(*), intent(in)                              :: vanam         !< Variable name
      character(*), intent(in)                              :: vacal         !< Calling subroutine name
      integer(8),   intent(inout)                           :: memor(2)      !< Memory counter
      type(ecg_coords),  intent(inout), pointer          :: varia(:)      !< Variable
      integer(4)                                            :: istat
      
      if( associated(varia) ) then
         
         lbytm = -storage_size(varia)/8*size(varia,1,kind=ip)
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

   end subroutine memory_deallo_ecg



end module mod_exm_memory
!> @}
