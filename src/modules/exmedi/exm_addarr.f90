!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_addarr.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   This routine computes some additional arrays
!> @details This routine computes some additional arrays
!> @} 
!-----------------------------------------------------------------------
subroutine exm_addarr()
  !-----------------------------------------------------------------------
  !****f* Exmedi/exm_addarr
  ! NAME 
  !    exm_addarr
  ! DESCRIPTION
  !   subroutine that computes the conductivity tensor for 2D and 3D 
  ! USES
  !    exm_...
  ! USED BY
  !    exm_doiter
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_exmedi
  use mod_maths,           only : maths_normalize_vector
  use mod_eccoupling,      only : kfl_exmsld_ecc, eccou_set_cell_type
  use def_kermod,          only : kfl_celltype_fun, kfl_fiber_tang_fun, kfl_fiber_norm_fun
!  use def_kermod,          only : kfl_fiber_long_fun
  use mod_messages,        only : messages_live 
  use mod_exm_diffusivity, only : exm_diffusivity_init
  use mod_biofibers,       only : kfl_biofibers, biofib_point_nodal_fibers_all


  implicit none
  logical(lg) :: ortho_model
  real(rp)    :: fib_aux(ndime), fib_aux_ortho2(ndime) !! Auxiliary fiber fields
!  real(rp)    :: fib_aux_ortho1(ndime)                 !! Auxiliary fiber fields

  ! Initialize arrays
  fib_aux(:)= 0.0_rp
  fib_aux_ortho2(:)= 0.0_rp

   ! Evaluation of orthotropic model
   ortho_model=.false.
   if( (kfl_fiber_tang_fun.ne.0_ip) .and. (kfl_fiber_norm_fun .ne. 0_ip) ) then
      ortho_model=.true.
   else if( (kfl_fiber_tang_fun.ne.0_ip) .or. (kfl_fiber_norm_fun .ne. 0_ip) ) then
      call runend('EXM_ADARR: ORTHO MODEL DEFINED BUT ONE FIBER FIELD IS MISSING')
   endif


  if( INOTEMPTY ) then
     !--------------------------------------------------------!  
     !--------------< INITIALISE FIBER FIELDS >---------------!
     !
     if( kfl_biofibers )then
       call biofib_point_nodal_fibers_all( fiber_exm, 'FIBER' ) 
       call biofib_point_nodal_fibers_all( sheet_exm, 'SHEET' ) 
       call biofib_point_nodal_fibers_all( normal_exm, 'NORMA'  ) 
     else
       call runend('EXM_ADDARR: FIBERS REQUIRED FOR ELECTROPHYSIOLOGY')
     endif
     !    
     ! For the isotropic part
     !    
!     if( kfl_fiber_long_fun < 0 ) then
!      call messages_live("EXMEDI IS USING FIELD "//trim(intost(-kfl_fiber_long_fun))//" TO READ FIBERS")
!      fiber_exm => xfiel(-kfl_fiber_long_fun) % a(:,:,1)
!     else if ( kfl_fiber_long_fun > 0 ) then
!        fib_aux(:)= 0.0_rp
!        if(kfl_fiber_long_fun.eq.1_ip) then
!          ! X ALIGNED
!          call messages_live("EXMEDI USING X ALIGNED FIBERS")
!          fib_aux(1) = 1.0_rp
!        elseif(kfl_fiber_long_fun.eq.2_ip) then
!          ! Y ALIGNED
!         call messages_live("EXMEDI USING Y ALIGNED FIBERS")
!         fib_aux(2) = 1.0_rp
!        elseif(kfl_fiber_long_fun.eq.3_ip) then
!          ! Z ALIGNED
!         call messages_live("EXMEDI USING Z ALIGNED FIBERS")
!          fib_aux(3) = 1.0_rp
!        else
!          call runend('EXM_ADDARR: FIBER PREDEFINED OPTION NOT RECOGNISED')
!        endif
!     else
!       call runend('EXM_ADDARR: FIBERS REQUIRED FOR ELECTROPHYSIOLOGY')
!     end if
!     !
!     ! Sheet and normal vectors
!     !
!     if ( ortho_model ) then
!        ! Sheet fibre field
!        if(kfl_fiber_tang_fun<0_ip)then
!          call messages_live("EXMEDI USING FIBER SHEET FIELD "//trim(intost(-kfl_fiber_tang_fun)))
!          sheet_exm => xfiel(-kfl_fiber_tang_fun) % a(:,:,1)
!        elseif(kfl_fiber_tang_fun>0_ip)then
!          fib_aux_ortho1(:)= 0.0_rp
!          if(kfl_fiber_tang_fun.eq.1_ip) then
!            ! X ALIGNED
!            call messages_live("EXMEDI USING X ALIGNED FIBER SHEET")
!            fib_aux_ortho1(1) = 1.0_rp
!          elseif(kfl_fiber_tang_fun.eq.2_ip) then
!            ! Y ALIGNED
!            call messages_live("EXMEDI USING Y ALIGNED FIBER SHEET")
!            fib_aux_ortho1(2) = 1.0_rp
!          elseif(kfl_fiber_tang_fun.eq.3_ip) then
!            ! Z ALIGNED
!            call messages_live("EXMEDI USING Z ALIGNED FIBER SHEET")
!            fib_aux_ortho1(3) = 1.0_rp
!          else
!            call runend('EXM_ADDARR: FIBER PREDEFINED OPTION NOT RECOGNISED')
!          endif
!        endif
!
!        ! Normal fibre field
!        if(kfl_fiber_norm_fun<0_ip)then
!          call messages_live("EXMEDI USING NORMAL FIBER SHEET FIELD "//trim(intost(-kfl_fiber_norm_fun)))
!          normal_exm => xfiel(-kfl_fiber_norm_fun) % a(:,:,1)
!        elseif(kfl_fiber_norm_fun>0_ip)then
!          fib_aux_ortho2(:)= 0.0_rp
!          if(kfl_fiber_norm_fun.eq.1_ip) then
!            ! X ALIGNED
!            call messages_live("EXMEDI USING X ALIGNED NORMAL FIBER SHEET")
!            fib_aux_ortho2(1) = 1.0_rp
!          elseif(kfl_fiber_norm_fun.eq.2_ip) then
!            ! Y ALIGNED
!            call messages_live("EXMEDI USING Y ALIGNED NORMAL FIBER SHEET")
!            fib_aux_ortho2(2) = 1.0_rp
!          elseif(kfl_fiber_norm_fun.eq.3_ip) then
!            ! Z ALIGNED
!            call messages_live("EXMEDI USING Z ALIGNED NORMAL FIBER SHEET")
!            fib_aux_ortho2(3) = 1.0_rp
!          else
!            call runend('EXM_ADDARR: FIBER PREDEFINED OPTION NOT RECOGNISED')
!          endif
!        endif
!     end if
     !----------///////////////////////////////---------------!  
     !----------< END INITIALISE FIBER FIELDS >---------------!

     !--------------------------------------------------------!  
     !--------------< INITIALISE CELL FIELDS >----------------!
     !
     ! CEDIF_EXM: Cell types
     !
     if( kfl_celltype_fun .ne. 0_ip ) then ! there is cell type field, otherwise celltype=1_ip everywhere
         celty_exm => xfiel(-kfl_celltype_fun) % a(1,:,1)

         call messages_live("EXMEDI IS USING FIELD "//trim(intost(-kfl_celltype_fun))//" TO READ CELLTYPE")

         if ( kfl_exmsld_ecc ) then
            call eccou_set_cell_type(celty_exm)
            call messages_live("EXMEDI PASSING CELLTYPE FIELD "//trim(intost(-kfl_celltype_fun))//" TO ECCOUPLING")
         endif
     else
         call messages_live("EXMEDI CELLTYPE FIELD UNDEFINED","WARNING")
     end if


     if( modab_exm .ne. 0_ip ) then
         atbhe_exm => xfiel(-modab_exm) % a(:,:,1)
         call messages_live("EXMEDI IS USING FIELD "//trim(intost(-modab_exm))//" TO READ APEX-BASE HETEROGENEITY")
         if (.not. associated( xfiel(-modab_exm) % a )) call runend("EXMEDI: APEX TO BASE HETER FIELD "//trim(intost(-modab_exm))//" DOES NOT EXIST")
     else
         call messages_live("EXMEDI APEX-BASE FIELD UNDEFINED","WARNING")
     end if

     if (modst_exm < 0_ip) then
         if (.not. associated( xfiel(-modst_exm) % a )) call runend("EXMEDI: STIMULI FIELD "//trim(intost(-modst_exm))//" DOES NOT EXIST")
     end if

     !----------///////////////////////////////---------------!  
     !----------< END INITIALISE CELL FIELDS >----------------!
     call exm_diffusivity_init()

     !--------------------------------------------------------!  
     !--------------< DIFFUSION TENSOR CEDIF >----------------!
  end if  

  
end subroutine exm_addarr
