!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> 
!> @author  houzeaux
!> @date    2021-11-11
!> @brief   Get properties
!> @details Get properties during an assembly
!>
!-----------------------------------------------------------------------

module mod_getpro

#include "def_vector_size.inc"
  use def_kintyp_basic,       only : ip,rp,lg
  use def_ker_proper,         only : typ_valpr_ker
  use def_ker_proper,         only : SCALAR_PROPERTY
  use def_ker_proper,         only : xxx_operations
  use def_ker_proper,         only : xxx_init
  use def_kermod,             only : densi_ker
  use def_kermod,             only : visco_ker
  use def_kermod,             only : poros_ker
  use def_kermod,             only : condu_ker
  use def_kermod,             only : sphea_ker
  use def_kermod,             only : dummy_ker
  use def_kermod,             only : turmu_ker
  use def_kermod,             only : absor_ker
  use def_kermod,             only : scatt_ker
  use def_kermod,             only : mixin_ker
  use def_kermod,             only : anipo_ker
  use def_kermod,             only : walvi_ker
  use mod_ker_proper_generic, only : ker_element_operations_value

#ifdef OPENACCHHH
#define DEF_VECT ivect
#else
#define DEF_VECT 1:VECTOR_SIZE
#endif
  !
  ! Private properties
  ! 
#ifndef PROPER_PRIVATE_OFF
  use mod_ker_proper_turmu
#endif
  
  implicit none
  private

  interface getpro_val
     module procedure getpro_val_2,&
          &           getpro_val_4
  end interface getpro_val
  
  public :: getpro_val
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-11-12
  !> @brief   Value of a property
  !> @details Compute the GP value of a property
  !> 
  !-----------------------------------------------------------------------
  
  subroutine getpro_val_2(wname,xvalu,pelty,pnode,pgaus,pmate,list_elements,lnods_loc,gpcar)

#if defined OPENACCHHH && defined _OPENACC 
    use openacc
#endif

    character(5),                  intent(in)    :: wname
    real(rp),                      intent(inout) :: xvalu(:,:)
    integer(ip),                   intent(in)    :: pelty
    integer(ip),                   intent(in)    :: pnode
    integer(ip),                   intent(in)    :: pgaus
    integer(ip),                   intent(in)    :: pmate
    integer(ip),                   intent(in)    :: list_elements(:)
    real(rp),            optional, intent(in)    :: gpcar(:,:,:,:)
    integer(ip),                   intent(in)    :: lnods_loc(:,:)

    type(typ_valpr_ker), pointer                 :: prope_ker
    integer(ip)                                  :: kk,ivect,igaus,idime,ilaws
#if defined(OPENACCHHH) && defined(_OPENACC)
    logical(lg)                                  :: lpalloc
#endif

#ifndef PROPER_PRIVATE_OFF
  procedure(xxx_operations), pointer             :: subru_operations
  procedure(xxx_init),       pointer             :: subru_init
#endif

    select case ( wname )
    case ( 'DENSI' )  ; prope_ker => densi_ker
#ifdef OPENACCHHH
       !$acc enter data copyin(densi_ker, densi_ker % value_const)
#endif
    case ( 'VISCO' )  ; prope_ker => visco_ker
#ifdef OPENACCHHH
       !$acc enter data copyin(visco_ker, visco_ker % value_const)
#endif 
    case ( 'MIXIN' )  ; prope_ker => mixin_ker
    case ( 'POROS' )  ; prope_ker => poros_ker
    case ( 'CONDU' )  ; prope_ker => condu_ker
    case ( 'SPHEA' )  ; prope_ker => sphea_ker
    case ( 'DUMMY' )  ; prope_ker => dummy_ker
    case ( 'TURBU' )  ; prope_ker => turmu_ker
    case ( 'ABSOR' )  ; prope_ker => absor_ker
    case ( 'SCATT' )  ; prope_ker => scatt_ker
    case ( 'ANIPO' )  ; prope_ker => anipo_ker
    case ( 'WALLV' )  ; prope_ker => walvi_ker         
    end select

    ilaws = prope_ker % ilaws(pmate)
    
    if( prope_ker % kfl_exist == 0 ) then
       !
       ! Put value to zero
       !
       xvalu = 0.0_rp
       
    else if ( prope_ker % llaws(ilaws) % where == 'CONST' ) then
       !
       ! Property is constant 
       !
       if( prope_ker % kfl_type == SCALAR_PROPERTY ) then
          
#ifdef OPENACCHHH
#ifdef _OPENACC
          lpalloc = acc_is_present(xvalu) 
          if ( .not. lpalloc) then
             !$acc enter data copyin (xvalu(1:VECTOR_SIZE,1:pgaus))
          end if
#endif
          !$acc parallel loop gang vector default(present) 
          do ivect = 1, VECTOR_SIZE
#endif                        
             do igaus = 1,pgaus
                xvalu(DEF_VECT,igaus) = prope_ker % value_const(1,pmate)
             end do
             
#ifdef OPENACCHHH                      
          end do
          !$acc end parallel loop
          
#ifdef _OPENACC
          if ( .not. lpalloc) then
             !$acc update self      (xvalu(1:VECTOR_SIZE,1:pgaus))
             !$acc exit data delete (xvalu(1:VECTOR_SIZE,1:pgaus)) 
          end if
#endif
#endif                 
          
       else
          
          kk = 0
          do igaus = 1,pgaus
             do idime = 1,prope_ker % dim
                kk = kk + 1
                xvalu(DEF_VECT,kk) = prope_ker % value_const(idime,pmate)
             end do
          end do
          
       end if
       
       
    else
       !
       ! Non constant
       !
       select case ( prope_ker % on_the_fly(pmate) )
          
       case ( 0_ip )
          !
          ! Properties copied from memory
          !
          ilaws = prope_ker % ilaws(pmate)
          
          select case ( prope_ker % llaws(ilaws) % where )
             
          case ( 'IELEM' ) 
             !
             ! Property is computed element-wise
             !  
             if( prope_ker % kfl_type == SCALAR_PROPERTY ) then
                
                do ivect = 1,VECTOR_SIZE
                   do igaus = 1,pgaus * prope_ker % dim
                      xvalu(ivect,igaus) = prope_ker % value_ielem(list_elements(ivect)) % a(igaus)
                   end do
                end do
                
             else 

                call runend('GETPRO: NOT CODED')
 
             end if
             
          end select
          
       case ( 1_ip )
          !
          ! On-the-fly
          !
#ifndef PROPER_PRIVATE_OFF
          
          select case ( prope_ker % wlaws(pmate) )
          case ( 'VRMAN' ) ; subru_operations => ker_proper_turmu_vreman_operations      ; subru_init => ker_proper_turmu_vreman_init
          case ( 'WALE ' ) ; subru_operations => ker_proper_turmu_wale_operations        ; subru_init => ker_proper_turmu_wale_init
          case ( 'ILSA ' ) ; subru_operations => ker_proper_turmu_ilsa_operations        ; subru_init => ker_proper_turmu_ilsa_init
          case ( 'SMAGO' ) ; subru_operations => ker_proper_turmu_smagorinsky_operations ; subru_init => ker_proper_turmu_vreman_init
          end select

          call ker_element_operations_value(&
               int(VECTOR_SIZE,ip),pelty,pnode,pgaus,pmate,list_elements,&
               lnods_loc,gpcar,prope_ker,subru_init,subru_operations,xvalu)
#endif
          
       end select

    end if

    !--------------------------------------------------------------------
    !
    ! Deallocate when copied from memory
    !
    !--------------------------------------------------------------------
    
    if( prope_ker % on_the_fly(pmate) == 0 ) then
       !
       ! This part deallocates the prope_ker from the GPU
       !
       select case ( wname )
       case ( 'DENSI' ) ; prope_ker => densi_ker
#ifdef OPENACCHHH
          !$acc exit data delete(densi_ker, densi_ker % value_const)
#endif 
       case ( 'VISCO' ) ; prope_ker => visco_ker
#ifdef OPENACCHHH
          !$acc exit data delete(visco_ker, visco_ker % value_const)
#endif 
       end select

    end if
    
  end subroutine getpro_val_2

  subroutine getpro_val_4(wname,xvalu,pelty,pnode,pgaus,pmate,list_elements,lnods_loc,gpcar)

#if defined OPENACCHHH && defined _OPENACC 
    use openacc
#endif

    character(5),                  intent(in)    :: wname
    real(rp),                      intent(inout) :: xvalu(:,:,:,:)
    integer(ip),                   intent(in)    :: pelty
    integer(ip),                   intent(in)    :: pnode
    integer(ip),                   intent(in)    :: pgaus
    integer(ip),                   intent(in)    :: pmate
    integer(ip),                   intent(in)    :: list_elements(:)
    real(rp),            optional, intent(in)    :: gpcar(:,:,:,:)
    integer(ip),                   intent(in)    :: lnods_loc(:,:)

    type(typ_valpr_ker), pointer                 :: prope_ker
    integer(ip)                                  :: kk,ivect,igaus,idime,jdime,ilaws,nn
#if defined(OPENACCHHH) && defined(_OPENACC)
    logical(lg)                                  :: lpalloc
#endif

#ifndef PROPER_PRIVATE_OFF
  procedure(xxx_operations), pointer             :: subru_operations
  procedure(xxx_init),       pointer             :: subru_init
#endif

    select case ( wname )
    case ( 'ANIPO' )  ; prope_ker => anipo_ker
    case default      ; call runend('GETPRO_VAL: PROPERTY SHOULD BE ANIPO')
    end select

    ilaws = prope_ker % ilaws(pmate)
    
    if( prope_ker % kfl_exist == 0 ) then
       !
       ! Put value to zero
       !
       xvalu = 0.0_rp
       
    else if ( prope_ker % llaws(ilaws) % where == 'CONST' ) then
       !
       ! Property is constant 
       !
       nn = int(sqrt(real(prope_ker % dim,rp)),ip)
       do igaus = 1,pgaus
          kk = 0
          do idime = 1,nn
             do jdime = 1,nn
                kk = kk + 1
                xvalu(DEF_VECT,jdime,idime,igaus) = prope_ker % value_const(kk,pmate)
             end do
          end do
       end do
       
    else if ( prope_ker % llaws(ilaws) % where == 'IELEM' ) then
       !
       ! Element wise
       !
       nn = int(sqrt(real(prope_ker % dim,rp)),ip)
       do ivect = 1,VECTOR_SIZE
          do igaus = 1,pgaus
             kk = 0
             do idime = 1,nn
                do jdime = 1,nn
                   kk = kk + 1
                   xvalu(ivect,jdime,idime,igaus) = prope_ker % value_ielem(list_elements(ivect)) % a((igaus-1)*nn*nn+kk)
                end do
             end do
          end do
       end do

    end if

  end subroutine getpro_val_4
  
end module mod_getpro
!> @}
