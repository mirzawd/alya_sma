!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    mod_par_element_loop.f90
!> @author  guillaume
!> @date    2021-10-20
!> @brief   Define arrays for hybrid parallelization
!> @details Define arrays for hybrid parallelization: includes OpenMP, 
!>          OmpSs, vectorization and CUDA.
!>          For element loops, elements are ordered according to
!>          pnode, pgaus, pmate
!>
!>          1. Without race condition:
!>          --------------------------
!>
!>             1.1 Classical element loop
!>             1.2 Otherwise:
!>                 do isubd = 1,num_subd_norace_par
!>                   do ipack = 1,num_pack_norace_par(isubd)
!>                     list_elements_norace_par
!>
!>          2. With race condition
!>          ----------------------
!>             Choose a hybrid method with flag 
!>             par_hybrid = PAR_OPENMP_COLORING
!>             par_hybrid = PAR_OPENMP_NO_COLORING
!>             par_hybrid = PAR_OMPSS
!>
!>             2.1 PAR_OPENMP_NO_COLORING:
!>                 define NO_COLORING macro and do not forget to put
!>                 ATOMIC pragma or define CRITICAL sections
!>                 using #ifdef NO_COLORING. Otherwise, DO NOT define
!>                 this macro.
!>             2.2 PAR_OMPSS:
!>                 Define ALYA_OMPSS macro
!>                 DO NOT define NO_COLORING macro!
!>             2.3 PAR_OPENMP_COLORING:
!>                 Do not do anything...
!>                 
!>          If no vectorization is used, add the following inner loop:
!>          do kelem = 1,VECTOR_SIZE
!>             ielem = list_elements_par(isubd) % packs(ipack) % l(kelem)
!>             if( ielem > 0 ) call element_assembly(ielem)
!>
!>          Important note
!>          --------------
!>          When running on GPUs, VECTOR_SIZE may be tool large for loops
!>          not ported to GPU. Then we should use
!>          (num_subd_cpu,num_pack_cpu,list_elements_cpu) and
!>          VECTOR_SIZE_CPU
!>
!>
!------------------------------------------------------------------------

module mod_par_element_loop

#include "def_vector_size.inc"
  use def_kintyp,         only : ip,rp,lg,i1p
  use def_master,         only : INOTMASTER,kfl_paral,intost
  use def_domain,         only : nelem,lnnod,lmate,lgaus
  use def_domain,         only : mgaus,mnode,nmate,ltype
  use def_domain,         only : ompss_domains
  use def_kermod,         only : kfl_vector_size
  use mod_parall,         only : typ_list_elements_par
  use mod_parall,         only : par_memor
  use mod_parall,         only : num_subd_par
  use mod_parall,         only : num_pack_par
  use mod_parall,         only : list_elements_par
  use mod_parall,         only : num_subd_cpu
  use mod_parall,         only : num_pack_cpu
  use mod_parall,         only : list_elements_cpu
  use mod_parall,         only : num_subd_norace_par
  use mod_parall,         only : num_pack_norace_par
  use mod_parall,         only : list_elements_norace_par
  use mod_parall,         only : num_subd_norace_cpu
  use mod_parall,         only : num_pack_norace_cpu
  use mod_parall,         only : list_elements_norace_cpu
  use mod_parall,         only : par_omp_num_threads
  use mod_parall,         only : par_omp_num_colors
  use mod_parall,         only : par_omp_ia_colors 
  use mod_parall,         only : par_omp_ja_colors
  use mod_parall,         only : par_hybrid
  use mod_parall,         only : PAR_OPENMP_COLORING
  use mod_parall,         only : PAR_OPENMP_NO_COLORING
  use mod_parall,         only : PAR_OMPSS
  use mod_parall,         only : PAR_HYBRID_OFF
  use mod_maths,          only : maths_mapping_1d_to_3d_x
  use mod_maths,          only : maths_mapping_1d_to_3d_y
  use mod_maths,          only : maths_mapping_1d_to_3d_z
  use mod_maths,          only : maths_mapping_3d_to_1d
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_SUM
  use mod_messages,       only : messages_live

  implicit none
  private

  integer(ip)             :: num_types
  integer(ip), pointer    :: num_element_types(:)
  integer(ip), pointer    :: element_types(:)

  public :: par_element_loop
  public :: par_boundary_loop

contains

  subroutine par_element_loop()

    integer(ip)             :: isubd,ipack,itype
    integer(ip)             :: kelem,ielem
    integer(ip)             :: pnode,pgaus,pmate
    integer(ip)             :: jtype
    integer(ip)             :: VECTOR_SIZE_LOC
    integer(ip)             :: VECTOR_SIZE_CPU_LOC

    integer(ip)             :: num_holes
    logical(lg), pointer    :: consider_element(:) 
    logical(lg)             :: if_cpu
    !
    ! Deallocate memory if necessary
    !
    call memory_deallo(par_memor,'NUM_PACK_PAR','par_element_loop',num_pack_par)
    if( associated(list_elements_par) ) then
       do isubd = 1,size(list_elements_par) 
          call memory_deallo(par_memor,'LIST_ELEMENTS_PAR % PACKS','par_element_loop',list_elements_par(isubd) % packs)
       end do
       deallocate(list_elements_par)
    end if

    call memory_deallo(par_memor,'NUM_PACK_CPU','par_element_loop',num_pack_cpu)
    if( associated(list_elements_cpu) ) then
       do isubd = 1,size(list_elements_cpu) 
          call memory_deallo(par_memor,'LIST_ELEMENTS_CPU % PACKS','par_element_loop',list_elements_cpu(isubd) % packs)
       end do
       deallocate(list_elements_cpu)
    end if

    call memory_deallo(par_memor,'NUM_PACK_NORACE_PAR','par_element_loop',num_pack_norace_par)
    if( associated(list_elements_norace_par) ) then
       do isubd = 1,size(list_elements_norace_par) 
          call memory_deallo(par_memor,'LIST_ELEMENTS_NORACE_PAR % PACKS','par_element_loop',list_elements_norace_par(isubd) % packs)
       end do
       deallocate(list_elements_norace_par)
    end if

    call memory_deallo(par_memor,'NUM_PACK_NORACE_CPU','par_element_loop',num_pack_norace_cpu)
    if( associated(list_elements_norace_cpu) ) then
       do isubd = 1,size(list_elements_norace_cpu) 
          call memory_deallo(par_memor,'LIST_ELEMENTS_NORACE_CPU % PACKS','par_element_loop',list_elements_norace_cpu(isubd) % packs)
       end do
       deallocate(list_elements_norace_cpu)
    end if
    !
    ! Errors
    !
    if( par_hybrid == PAR_HYBRID_OFF .and. par_omp_num_threads /= 0 ) then
       call runend('PARALL: CHOOSE A HYBRID PARALLELIZATION STRATEGY')
    end if
    if( par_hybrid /= PAR_HYBRID_OFF .and. par_omp_num_threads == 0 ) then
       call runend('PARALL: SET ENVIRONMENT VARIABLE OMP_NUM_THREAD')
    end if
    !
    ! Vector size
    !
    VECTOR_SIZE_LOC     = max(1_ip,int(VECTOR_SIZE,ip))
    VECTOR_SIZE_CPU_LOC = max(1_ip,int(VECTOR_SIZE_CPU,ip))

    !
    ! Check if a CPU version is required
    !
    if( VECTOR_SIZE_LOC /= VECTOR_SIZE_CPU_LOC ) then
       if_cpu = .true.
    else
       if_cpu = .false.
    end if
    
    call messages_live('HYBRID PARALLELIZATION & VECTORIZATION (ELEMENTS)','START SECTION')
    if( VECTOR_SIZE > 1 ) then
       call messages_live('VECTORIZATION WITH VECTOR SIZE= '//trim(intost(int(VECTOR_SIZE_LOC,ip))))
    else     
       call messages_live('NO VECTORIZATION')
    end if
    if( if_cpu ) then
       call messages_live('VECTORIZATION FOR CPU WITH VECTOR SIZE= '//trim(intost(int(VECTOR_SIZE_CPU_LOC,ip))))
    end if
    
    num_types = 0
    num_holes = 0
    jtype     = 0
    
    if( INOTMASTER ) then
       !
       ! Initialization
       !
       num_subd_cpu        = 0
       num_subd_par        = 0
       num_subd_norace_par = 0
       num_subd_norace_cpu = 0
       nullify(num_pack_par)
       nullify(list_elements_par)
       nullify(num_pack_cpu)
       nullify(list_elements_cpu)
       nullify(num_pack_norace_par)
       nullify(list_elements_norace_par)
       nullify(num_pack_norace_cpu)
       nullify(list_elements_norace_cpu)
       !
       ! Local variables
       !
       nullify(num_element_types)
       nullify(consider_element)
       nullify(element_types)
       !
       ! Number and list of element types: NUM_ELEMENT_TYPES, ELEMENT_TYPES
       !
       if( VECTOR_SIZE_LOC == 1 ) then
          num_types = 1
          num_holes = 0
       else
          num_types = nmate * mnode * mgaus
          num_holes = 0
       end if

       call memory_alloca(par_memor,'NUM_ELEMENT_TYPES','par_element_loop',num_element_types,num_types)
       call memory_alloca(par_memor,'ELEMENT_TYPES'    ,'par_element_loop',element_types,nelem)

       if( VECTOR_SIZE_LOC == 1 ) then
          num_element_types(1) = 0
          do ielem = 1,nelem
             if( ltype(ielem) > 0 ) then
                element_types(ielem) = 1
                num_element_types(1) = num_element_types(1) + 1
             else
                num_holes = num_holes + 1
             end if
          end do
          jtype = 1
       else
          do ielem = 1,nelem
             if( ltype(ielem) > 0 ) then
                pmate                    = lmate(ielem)
                pnode                    = lnnod(ielem)
                pgaus                    = lgaus(ielem)
                itype                    = maths_mapping_3d_to_1d(nmate,mnode,mgaus,pmate,pnode,pgaus)
                element_types(ielem)     = itype
                num_element_types(itype) = num_element_types(itype) + 1
             else
                num_holes = num_holes + 1
             end if
          end do
          jtype = 0
          do itype = 1,num_types
             jtype = jtype + min(1_ip,num_element_types(itype))
          end do
       end if

    end if

    call PAR_MAX(jtype)
    call PAR_SUM(num_holes)
    call messages_live('MAX NUMBER OF ELEMENT TYPES DETECTED= '//trim(intost(jtype)))
    if( num_holes > 0 ) then
       call messages_live('NUMBER OF ELEMENT OF NON-FINITE ELEMENT TYPE= '//trim(intost(num_holes)))     
    end if
    
    !----------------------------------------------------------------------
    !
    ! Numboer of subdomains and allocation
    !
    !----------------------------------------------------------------------

    if( INOTMASTER ) then
       !
       ! No race: NUM_SUBD_NORACE_PAR, LIST_ELEMENTS_NORACE, NUM_PACK_NORACE_PAR
       !
       num_subd_norace_par = 1
       allocate( list_elements_norace_par(num_subd_norace_par) )
       do isubd = 1,num_subd_norace_par
          nullify(list_elements_norace_par(isubd) % packs)
       end do
       call memory_alloca(par_memor,'NUM_PACK_NORACE_PAR','par_element_loop',num_pack_norace_par,num_subd_norace_par)
       
       if( if_cpu ) then
          num_subd_norace_cpu = 1
          allocate( list_elements_norace_cpu(num_subd_norace_cpu) )
          do isubd = 1,num_subd_norace_cpu
             nullify(list_elements_norace_cpu(isubd) % packs)
          end do
          call memory_alloca(par_memor,'NUM_PACK_NORACE_CPU','par_element_loop',num_pack_norace_cpu,num_subd_norace_cpu)          
       end if
       !
       ! With race: NUM_SUBD_PAR, LIST_ELEMENTS_PAR, NUM_PACK_PAR
       !
       if(      par_hybrid == PAR_HYBRID_OFF ) then
          num_subd_par = 1
       else if( par_hybrid == PAR_OPENMP_NO_COLORING ) then
          num_subd_par = 1
       else if( par_hybrid == PAR_OPENMP_COLORING ) then
          num_subd_par = par_omp_num_colors
       else if( par_hybrid == PAR_OMPSS ) then
          num_subd_par = size(ompss_domains,KIND=ip)
       end if
       num_subd_cpu = num_subd_par

       allocate( list_elements_par(num_subd_par) )
       do isubd = 1,num_subd_par
          nullify(list_elements_par(isubd) % packs)
       end do
       call memory_alloca(par_memor,'NUM_PACK_PAR','par_element_loop',num_pack_par,num_subd_par)

       if( if_cpu ) then
          allocate( list_elements_cpu(num_subd_cpu) )
          do isubd = 1,num_subd_cpu
             nullify(list_elements_cpu(isubd) % packs)
          end do
          call memory_alloca(par_memor,'NUM_PACK_CPU','par_element_loop',num_pack_cpu,num_subd_cpu)
       end if
       
    end if
    
    !----------------------------------------------------------------------
    !
    ! Loops without race condition: classical OpenMP
    !
    ! NUM_SUBD_NORACE_PAR = 1
    ! NUM_PACK_NORACE_PAR
    ! LIST_ELEMENTS_NORACE_PAR
    !
    !----------------------------------------------------------------------

    if( par_hybrid /= PAR_HYBRID_OFF ) then
       call messages_live('NUMBER OF THREADS= '//trim(intost(par_omp_num_threads)))
       call messages_live('LOOP WITHOUT RACE CONDITION: OPENMP WITHOUT COLORING')
    else
       call messages_live('OPENMP/OMPSS DESACTIVATED IN ALL LOOPS')     
    end if

    call par_element_loop_omp_norace(num_subd_norace_par,num_pack_norace_par,list_elements_norace_par,VECTOR_SIZE_LOC,'NORACE_PAR')
    if( if_cpu ) then
       call par_element_loop_omp_norace(num_subd_norace_cpu,num_pack_norace_cpu,list_elements_norace_cpu,VECTOR_SIZE_CPU_LOC,'NORACE_CPU')       
    end if

    !----------------------------------------------------------------------
    !
    ! Loops with race condition
    !
    ! NUM_SUBD_NORACE_PAR = 1
    ! NUM_PACK_NORACE_PAR
    ! LIST_ELEMENTS_NORACE_PAR
    !
    !----------------------------------------------------------------------

    if( par_hybrid == PAR_HYBRID_OFF .or. par_hybrid == PAR_OPENMP_NO_COLORING ) then
       !
       ! No hybrid parallelization 
       !
       if( par_hybrid /= PAR_HYBRID_OFF ) then
          call messages_live('LOOP WITH RACE CONDITION: OPENMP WITHOUT COLORING')
       end if

       call par_element_loop_omp_no_coloring(num_subd_par,num_pack_par,list_elements_par,VECTOR_SIZE_LOC,'PAR')
       if( if_cpu ) then          
          call par_element_loop_omp_no_coloring(num_subd_cpu,num_pack_cpu,list_elements_cpu,VECTOR_SIZE_CPU_LOC,'CPU')
       end if

    else if( par_hybrid == PAR_OPENMP_COLORING  ) then
       !
       ! OpenMP with coloring
       !
       call messages_live('LOOP WITH RACE CONDITION: OPENMP WITH COLORING')

       call par_element_loop_omp_coloring(num_subd_par,num_pack_par,list_elements_par,VECTOR_SIZE_LOC,'PAR')
       if( if_cpu ) then          
          call par_element_loop_omp_coloring(num_subd_cpu,num_pack_cpu,list_elements_cpu,VECTOR_SIZE_CPU_LOC,'CPU')          
       end if

    else if( par_hybrid == PAR_OMPSS ) then
       !
       ! OmpSs
       !
       call messages_live('LOOP WITH RACE CONDITION: OMPSS WITH MULDIDEPENDENCIES')

       call par_element_loop_ompss(num_subd_par,num_pack_par,list_elements_par,VECTOR_SIZE_LOC,'PAR')
       if( if_cpu ) then          
          call par_element_loop_omp_coloring(num_subd_cpu,num_pack_cpu,list_elements_cpu,VECTOR_SIZE_CPU_LOC,'CPU')          
       end if
       
    end if

    !-------------------------------------------------------------------
    !
    ! Some checks
    !
    !-------------------------------------------------------------------
    
    if( INOTMASTER .and. nelem > 0 ) then

       if( 1 == 1 ) then
          call memory_alloca(par_memor,'CONSIDER_ELEMENT','par_element_loop',consider_element  ,nelem)
          do isubd = 1,num_subd_par
             do ipack = 1,num_pack_par(isubd)
                do kelem = 1,size(list_elements_par(isubd) % packs(ipack) % l)
                   ielem = list_elements_par(isubd) % packs(ipack) % l(kelem)
                   if( ielem > 0 ) then
                      if( ltype(ielem) > 0 ) then
                         if( consider_element(ielem)) call runend('ELEMENT ALREADY TOUCHED')
                         consider_element(ielem) = .true.
                      end if
                   end if
                end do
             end do
          end do
          kelem = 0
          do ielem = 1,nelem
             if( .not. consider_element(ielem) .and. ltype(ielem) > 0 ) kelem=kelem+1 
          end do
          if(kelem>0) call runend('ERROR: '//trim(intost(kelem))//' ELEMENTS NOT TOUCHED')
          do isubd = 1,num_subd_par
             do ipack = 1,num_pack_par(isubd)
                ielem = list_elements_par(isubd) % packs(ipack) % l(1)
                itype = element_types(ielem)
                do kelem = 2,size(list_elements_par(isubd) % packs(ipack) % l)
                   ielem = list_elements_par(isubd) % packs(ipack) % l(kelem)
                   if( ielem > 0 ) then
                      if( element_types(ielem) /= itype ) then
                         print*,'elem 1= ',list_elements_par(isubd) % packs(ipack) % l(1),itype
                         print*,'elem 2= ',ielem,element_types(ielem)
                         call runend('ERROR: ELEMENT OF DIFFERENT TYPES IN THE SAME LIST')
                      end if
                   end if
                end do
             end do
          end do
          call memory_deallo(par_memor,'CONSIDER_ELEMENT','par_element_loop',consider_element)
       end if
       !
       ! Deallocate memory
       !
       call memory_deallo(par_memor,'ELEMENT_TYPES'    ,'par_element_loop',element_types)
       call memory_deallo(par_memor,'NUM_ELEMENT_TYPES','par_element_loop',num_element_types)

    end if

    call messages_live('HYBRID PARALLELIZATION & VECTORIZATION (ELEMENTS)','END SECTION')

  end subroutine par_element_loop

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-10-20
  !> @brief   List for OpenMP with coloring
  !> @details List for OpenMP with coloring
  !> 
  !-----------------------------------------------------------------------

  subroutine par_element_loop_omp_coloring(num_subd,num_pack,list_elements,VECTOR_SIZE_LOC,ACRONYM)

    integer(ip),                          intent(in)    :: num_subd
    integer(ip),                 pointer, intent(inout) :: num_pack(:)
    type(typ_list_elements_par), pointer, intent(inout) :: list_elements(:)
    integer(ip),                          intent(in)    :: VECTOR_SIZE_LOC
    character(LEN=*),                     intent(in)    :: ACRONYM

    integer(ip)                                         :: isubd,ipack,itype
    integer(ip)                                         :: ivect,kelem,ielem
    integer(ip)                                         :: jtype

    if( INOTMASTER ) then

       do isubd = 1,num_subd
          num_element_types = 0
          do kelem = par_omp_ia_colors(isubd),par_omp_ia_colors(isubd+1)-1
             ielem = par_omp_ja_colors(kelem)
             if( ltype(ielem) > 0 ) then
                itype = element_types(ielem)
                num_element_types(itype) = num_element_types(itype) + 1
             end if
          end do

          num_pack(isubd) = 0
          do itype = 1,num_types
             if( num_element_types(itype) > 0 ) then
                num_pack(isubd) = num_pack(isubd) + ceiling(real(num_element_types(itype),rp) / real(VECTOR_SIZE_LOC,rp), ip )
             end if
          end do

          call memory_alloca(par_memor,'LIST_ELEMENTS_'//ACRONYM//' % PACKS','par_element_loop',list_elements(isubd) % packs,num_pack(isubd))
          do ipack = 1,num_pack(isubd)
             call memory_alloca(par_memor,'LIST_ELEMENTS_'//ACRONYM//' % PACKS % L','par_element_loop',list_elements(isubd) % packs(ipack) % l,VECTOR_SIZE_LOC) 
          end do
          ivect = 0
          ipack = 0 
          jtype = 0
          do itype = 1,num_types
             if( num_element_types(itype) > 0 ) then
                do kelem = par_omp_ia_colors(isubd),par_omp_ia_colors(isubd+1)-1
                   ielem = par_omp_ja_colors(kelem)
                   if( ltype(ielem) > 0 .and. element_types(ielem) == itype ) then
                      ivect = ivect + 1
                      if( ivect > VECTOR_SIZE_LOC .or. itype /= jtype ) then                    
                         ivect = 1
                         ipack = ipack + 1
                      end if
                      if( ipack > size(list_elements(isubd) % packs) ) call runend('WRONG PACK NUMBER')
                      list_elements(isubd) % packs(ipack) % l(ivect) = ielem
                      jtype = itype
                   end if
                end do
             end if
          end do
       end do

    end if

  end subroutine par_element_loop_omp_coloring

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-10-20
  !> @brief   List for OpenMP with coloring
  !> @details List for OpenMP with coloring
  !> 
  !-----------------------------------------------------------------------

  subroutine par_element_loop_omp_no_coloring(num_subd,num_pack,list_elements,VECTOR_SIZE_LOC,ACRONYM)

    integer(ip),                          intent(in)    :: num_subd
    integer(ip),                 pointer, intent(inout) :: num_pack(:)
    type(typ_list_elements_par), pointer, intent(inout) :: list_elements(:)
    integer(ip),                          intent(in)    :: VECTOR_SIZE_LOC
    character(LEN=*),                     intent(in)    :: ACRONYM

    integer(ip)                                         :: isubd,ipack,itype
    integer(ip)                                         :: ivect,ielem
    integer(ip)                                         :: jtype

    if( INOTMASTER ) then

       isubd = 1
       if( num_types > 0 ) then
          num_element_types = 0
          do ielem = 1,nelem
             if( ltype(ielem) > 0 ) then
                itype = element_types(ielem)
                num_element_types(itype) = num_element_types(itype) + 1
             end if
          end do
       end if

       if( kfl_vector_size == -2 ) then
          num_pack(isubd) = 0 
          do itype = 1,num_types
             if( num_element_types(itype) > 0 ) then
                num_pack(isubd) = num_pack(isubd) + 1
             end if
          end do
       else
          num_pack(isubd) = 0 
          do itype = 1,num_types
             if( num_element_types(itype) > 0 ) then
                num_pack(isubd) = num_pack(isubd) + ceiling(real(num_element_types(itype),rp) / real(VECTOR_SIZE_LOC,rp), ip )
             end if
          end do
       end if

       call memory_alloca(par_memor,'LIST_ELEMENTS_'//ACRONYM//' % PACKS','par_element_loop',list_elements(isubd) % packs,num_pack(isubd))
       if( kfl_vector_size == -2 ) then
          do itype = 1,num_types
             call memory_alloca(par_memor,'LIST_ELEMENTS_'//ACRONYM//' % PACKS % L','par_element_loop',list_elements(isubd) % packs(itype) % l,num_element_types(itype)) 
          end do
       else
          do ipack = 1,num_pack(isubd)
             call memory_alloca(par_memor,'LIST_ELEMENTS_'//ACRONYM//' % PACKS % L','par_element_loop',list_elements(isubd) % packs(ipack) % l,VECTOR_SIZE_LOC) 
          end do
       end if
       
       ivect = 0
       ipack = 0 
       jtype = 0
       if( kfl_vector_size == -2 ) then
          do itype = 1,num_types
             if( num_element_types(itype) > 0 ) then
                do ielem = 1,nelem
                   if( ltype(ielem) > 0 .and. element_types(ielem) == itype ) then
                      ivect = ivect + 1
                      if( ivect > num_element_types(itype) .or. itype /= jtype ) then                    
                         ivect = 1
                         ipack = ipack + 1
                      end if
                      if( ipack > size(list_elements(isubd) % packs) ) call runend('WRONG PACK NUMBER')
                      list_elements(isubd) % packs(ipack) % l(ivect) = ielem
                      jtype = itype
                   end if
                end do
             end if
          end do
       else
          do itype = 1,num_types
             if( num_element_types(itype) > 0 ) then
                do ielem = 1,nelem
                   if( ltype(ielem) > 0 .and. element_types(ielem) == itype ) then
                      ivect = ivect + 1
                      if( ivect > VECTOR_SIZE_LOC .or. itype /= jtype ) then                    
                         ivect = 1
                         ipack = ipack + 1
                      end if
                      if( ipack > size(list_elements(isubd) % packs) ) call runend('WRONG PACK NUMBER')
                      list_elements(isubd) % packs(ipack) % l(ivect) = ielem
                      jtype = itype
                   end if
                end do
             end if
          end do
       end if
       
    end if

  end subroutine par_element_loop_omp_no_coloring

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-10-20
  !> @brief   List for OmpSs
  !> @details List for OmpSs
  !> 
  !-----------------------------------------------------------------------

  subroutine par_element_loop_ompss(num_subd,num_pack,list_elements,VECTOR_SIZE_LOC,ACRONYM)

    integer(ip),                          intent(in)    :: num_subd
    integer(ip),                 pointer, intent(inout) :: num_pack(:)
    type(typ_list_elements_par), pointer, intent(inout) :: list_elements(:)
    integer(ip),                          intent(in)    :: VECTOR_SIZE_LOC
    character(LEN=*),                     intent(in)    :: ACRONYM

    integer(ip)                                         :: isubd,ipack,itype
    integer(ip)                                         :: ivect,kelem,ielem
    integer(ip)                                         :: jtype

    if( INOTMASTER ) then
       
       if( num_types > 0 ) then
          
          do isubd = 1,num_subd
             num_element_types = 0
             do kelem = 1,size(ompss_domains(isubd) % elements)
                ielem = ompss_domains(isubd) % elements(kelem)
                if( ltype(ielem) > 0 ) then
                   itype = element_types(ielem)
                   num_element_types(itype) = num_element_types(itype) + 1
                end if
             end do

             num_pack(isubd) = 0 
             do itype = 1,num_types
                if( num_element_types(itype) > 0 ) then
                   num_pack(isubd) = num_pack(isubd) + ceiling(real(num_element_types(itype),rp) / real(VECTOR_SIZE_LOC,rp), ip )
                end if
             end do
             call memory_alloca(par_memor,'LIST_ELEMENTS_'//ACRONYM//' % PACKS','par_element_loop',list_elements(isubd) % packs,num_pack(isubd))
             do ipack = 1,num_pack(isubd)
                call memory_alloca(par_memor,'LIST_ELEMENTS_'//ACRONYM//' % PACKS % L','par_element_loop',list_elements(isubd) % packs(ipack) % l,VECTOR_SIZE_LOC) 
             end do
             ivect = 0
             ipack = 0 
             jtype = 0
             do itype = 1,num_types
                if( num_element_types(itype) > 0 ) then
                   do kelem = 1,size(ompss_domains(isubd) % elements)
                      ielem = ompss_domains(isubd) % elements(kelem)
                      if( ltype(ielem) > 0 .and. element_types(ielem) == itype ) then
                         ivect = ivect + 1
                         if( ivect > VECTOR_SIZE_LOC .or. itype /= jtype ) then                    
                            ivect = 1
                            ipack = ipack + 1
                         end if
                         if( ipack > size(list_elements(isubd) % packs) ) call runend('WRONG PACK NUMBER')
                         list_elements(isubd) % packs(ipack) % l(ivect) = ielem
                         jtype = itype
                      end if
                   end do
                end if
             end do
          end do
       end if

    end if

  end subroutine par_element_loop_ompss

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-10-20
  !> @brief   List for OMP without race condition
  !> @details List for OMP without race condition
  !> 
  !-----------------------------------------------------------------------

  subroutine par_element_loop_omp_norace(num_subd,num_pack,list_elements,VECTOR_SIZE_LOC,ACRONYM)

    integer(ip),                          intent(in)    :: num_subd
    integer(ip),                 pointer, intent(inout) :: num_pack(:)
    type(typ_list_elements_par), pointer, intent(inout) :: list_elements(:)
    integer(ip),                          intent(in)    :: VECTOR_SIZE_LOC
    character(LEN=*),                     intent(in)    :: ACRONYM

    integer(ip)                                         :: isubd,ipack,itype
    integer(ip)                                         :: ivect,ielem
    integer(ip)                                         :: jtype

    if( INOTMASTER ) then

       isubd = 1
       if( num_types > 0 ) then
          num_element_types = 0
          do ielem = 1,nelem
             if( ltype(ielem) > 0 ) then
                itype = element_types(ielem)
                num_element_types(itype) = num_element_types(itype) + 1
             end if
          end do
       end if

       num_pack(isubd) = 0 
       do itype = 1,num_types
          if( num_element_types(itype) > 0 ) then
             num_pack(isubd) = num_pack(isubd) + ceiling(real(num_element_types(itype),rp) / real(VECTOR_SIZE_LOC,rp), ip )
          end if
       end do
       call memory_alloca(par_memor,'LIST_ELEMENTS_'//ACRONYM//' % PACKS','par_element_loop',list_elements(isubd) % packs,num_pack(isubd))

       do ipack = 1,num_pack(isubd)
          call memory_alloca(par_memor,'LIST_ELEMENTS_'//ACRONYM//' % PACKS % L','par_element_loop',list_elements(isubd) % packs(ipack) % l,VECTOR_SIZE_LOC) 
       end do

       ivect = 0
       ipack = 0 
       jtype = 0
       do itype = 1,num_types
          if( num_element_types(itype) > 0 ) then
             do ielem = 1,nelem
                if( ltype(ielem) > 0 .and. element_types(ielem) == itype ) then
                   ivect = ivect + 1
                   if( ivect > VECTOR_SIZE_LOC .or. itype /= jtype ) then                    
                      ivect = 1
                      ipack = ipack + 1
                   end if
                   if( ipack > size(list_elements(isubd) % packs) ) &
                        call runend('WRONG PACK NUMBER: CALL GUILLAUME TO FIX IT')
                   list_elements(isubd) % packs(ipack) % l(ivect) = ielem
                   jtype = itype
                end if
             end do
          end if
       end do

    end if

  end subroutine par_element_loop_omp_norace

subroutine par_boundary_loop()
  use def_kintyp,         only : ip,rp,lg,i1p
  use def_master,         only : INOTMASTER,kfl_paral,intost
  use def_domain,         only : nboun,lnnob,lmate,ngaus,lelbo
  use def_domain,         only : mgaub,mnodb,nmate,ltypb,lnnod
  use def_domain,         only : ompss_boundaries
  use mod_parall,         only : par_memor
  use mod_parall,         only : num_subd_nboun_par
  use mod_parall,         only : num_pack_nboun_par
  use mod_parall,         only : list_boundaries_par
  use mod_parall,         only : num_subd_norace_nboun_par
  use mod_parall,         only : num_pack_norace_nboun_par
  use mod_parall,         only : list_boundaries_norace_par
  use mod_parall,         only : par_omp_num_threads
  use mod_parall,         only : par_omp_nboun_num_colors
  use mod_parall,         only : par_omp_nboun_ia_colors 
  use mod_parall,         only : par_omp_nboun_ja_colors
  use mod_parall,         only : par_hybrid
  use mod_parall,         only : PAR_OPENMP_COLORING 
  use mod_parall,         only : PAR_OPENMP_NO_COLORING
  use mod_parall,         only : PAR_OMPSS
  use mod_parall,         only : PAR_HYBRID_OFF
  use mod_maths,          only : maths_mapping_1d_to_3d_x
  use mod_maths,          only : maths_mapping_1d_to_3d_y
  use mod_maths,          only : maths_mapping_1d_to_3d_z
  use mod_maths,          only : maths_mapping_3d_to_1d
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_memory,         only : memory_size
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_SUM
  use mod_messages,       only : messages_live
  implicit none

  integer(ip)             :: isubd,ipack,itype
  integer(ip)             :: ivect,kelem,iboun
  integer(ip)             :: pnode,pgaub,pmate
  integer(ip)             :: inode
  integer(ip)             :: jtype,ielem
  integer(ip)             :: VECTOR_SIZE_LOC
  integer(ip)             :: pnodb

  integer(ip)             :: num_types
  integer(ip)             :: num_holes
  integer(ip), pointer    :: num_element_types(:)
  logical(lg), pointer    :: consider_element(:) 
  integer(ip), pointer    :: element_types(:)
  logical(lg)             :: adjust_vector_size
  logical(lg)             :: if_vector_size

  if_vector_size = .true.
  !
  ! Deallocate memory if necessary
  !
  call memory_deallo(par_memor,'NUM_PACK_NBOUN_PAR','par_boundary_loop',num_pack_nboun_par)
  if( associated(list_boundaries_par) ) then
     do isubd = 1,size(list_boundaries_par) 
        call memory_deallo(par_memor,'LIST_BOUNDARIES_PAR % PACKS','par_boundary_loop',list_boundaries_par(isubd) % packs)
     end do
     deallocate(list_boundaries_par)
  end if
  call memory_deallo(par_memor,'NUM_PACK_NORACE_NBOUN_PAR','par_boundary_loop',num_pack_norace_nboun_par)
  if( associated(list_boundaries_norace_par) ) then
     do isubd = 1,size(list_boundaries_norace_par) 
        call memory_deallo(par_memor,'LIST_BOUNDARIES_NORACE_PAR % PACKS','par_boundary_loop',list_boundaries_norace_par(isubd) % packs)
     end do
     deallocate(list_boundaries_norace_par)
  end if
  !
  ! Errors
  !
  if( par_hybrid == PAR_HYBRID_OFF .and. par_omp_num_threads /= 0 ) then
     call runend('PARALL: CHOOSE A HYBRID PARALLELIZATION STRATEGY')
  end if
  if( par_hybrid /= PAR_HYBRID_OFF .and. par_omp_num_threads == 0 ) then
     call runend('PARALL: SET ENVIRONMENT VARIABLE OMP_NUM_THREAD')
  end if
  !
  ! Vector size
  !
  adjust_vector_size = .false.
  if( if_vector_size ) then
     if(      VECTOR_SIZE == -1 ) then
        adjust_vector_size = .true.
     else
        VECTOR_SIZE_LOC = max(1_ip,int(VECTOR_SIZE,ip))
     end if
  else
     VECTOR_SIZE_LOC = 1
  end if

  call messages_live('HYBRID PARALLELIZATION & VECTORIZATION (BOUNDARIES)','START SECTION')
  if( if_vector_size ) then
     if( VECTOR_SIZE > 1 ) then
        call messages_live('VECTORIZATION WITH VECTOR SIZE= '//trim(intost(int(VECTOR_SIZE_LOC,ip))))
     else if( VECTOR_SIZE == -1 ) then
        call messages_live('VECTORIZATION AUTOMATICALLY ADJUSTED')
     else     
        call messages_live('NO VECTORIZATION')
     end if
  end if

  num_types = 0
  num_holes = 0
  jtype     = 0
  
  if( INOTMASTER ) then
     !
     ! Initialization
     !
     num_subd_nboun_par = 0
     num_subd_norace_nboun_par = 0
     nullify(num_pack_nboun_par)
     nullify(list_boundaries_par)
     nullify(num_pack_norace_nboun_par)
     nullify(list_boundaries_norace_par)
     !
     ! Local variables
     !
     nullify(num_element_types)
     nullify(consider_element)
     nullify(element_types)
     !
     ! Number and list of element types: NUM_ELEMENT_TYPES, ELEMENT_TYPES
     !
     if( VECTOR_SIZE_LOC == 1 .and. .not. adjust_vector_size ) then
        num_types = 1
        num_holes = 0
     else
        num_types = nmate * mnodb * mgaub
        num_holes = 0
     end if

     call memory_alloca(par_memor,'NUM_ELEMENT_TYPES','par_boundary_loop',num_element_types,num_types)
     call memory_alloca(par_memor,'ELEMENT_TYPES'    ,'par_boundary_loop',element_types,nboun)
     
     if( VECTOR_SIZE_LOC == 1 ) then
        num_element_types(1) = 0
        do iboun = 1,nboun
           if( ltypb(iboun) > 0 ) then
              element_types(iboun) = 1
              num_element_types(1) = num_element_types(1) + 1
           else
              num_holes = num_holes + 1
           end if
        end do
        jtype = 1
     else
        !
        ! This could fail because the same boudnary type could belong to different element types
        !
        do iboun = 1,nboun
           if( ltypb(iboun) > 0 ) then
              ielem                    = lelbo(iboun)
              pmate                    = lmate(ielem)
              pnode                    = lnnod(ielem)
              pnodb                    = lnnob(iboun)
              pgaub                    = ngaus(ltypb(iboun))
              itype                    = maths_mapping_3d_to_1d(nmate,mnodb,mgaub,pmate,pnodb,pgaub)
              element_types(iboun)     = itype
              num_element_types(itype) = num_element_types(itype) + 1
           else
              num_holes = num_holes + 1
           end if
        end do
        jtype = 0
        do itype = 1,num_types
           jtype = jtype + min(1_ip,num_element_types(itype))
        end do
     end if

  end if

  call PAR_MAX(jtype)
  call PAR_SUM(num_holes)
  call messages_live('MAX NUMBER OF BOUNDARY TYPES DETECTED= '//trim(intost(jtype)))
  if( num_holes > 0 ) then
     call messages_live('NUMBER OF HOLE BOUNDARIES= '//trim(intost(num_holes)))     
  end if

  !----------------------------------------------------------------------
  !
  ! Loops without race condition: classical OpenMP
  !
  ! NUM_SUBD_NORACE_NBOUN_PAR = 1
  ! NUM_PACK_NORACE_NBOUN_PAR
  ! LIST_BOUNDARIES_NORACE_PAR
  !
  !----------------------------------------------------------------------
  
  if( par_hybrid /= PAR_HYBRID_OFF ) then
     call messages_live('NUMBER OF THREADS= '//trim(intost(par_omp_num_threads)))
     call messages_live('LOOP WITHOUT RACE CONDITION: OPENMP WITHOUT COLORING')
  else
     call messages_live('OPENMP/OMPSS DESACTIVATED IN ALL LOOPS')     
  end if
  
  if( INOTMASTER ) then

     num_subd_norace_nboun_par = 1
     allocate( list_boundaries_norace_par(num_subd_norace_nboun_par) )
     do isubd = 1,num_subd_norace_nboun_par
        nullify(list_boundaries_norace_par(isubd) % packs)
     end do
     call memory_alloca(par_memor,'NUM_PACK_NORACE_NBOUN_PAR','par_boundary_loop',num_pack_norace_nboun_par,num_subd_norace_nboun_par)

     isubd = 1
     num_element_types = 0
     do iboun = 1,nboun
        if( ltypb(iboun) > 0 ) then
           itype = element_types(iboun)
           num_element_types(itype) = num_element_types(itype) + 1
        end if
     end do

     if( adjust_vector_size ) then
        num_pack_norace_nboun_par(isubd) = 0 
        do itype = 1,num_types
           if( num_element_types(itype) > 0 ) then
              num_pack_norace_nboun_par(isubd) = num_pack_norace_nboun_par(isubd) + 1
           end if
        end do
     else
        num_pack_norace_nboun_par(isubd) = 0 
        do itype = 1,num_types
           if( num_element_types(itype) > 0 ) then
              num_pack_norace_nboun_par(isubd) = num_pack_norace_nboun_par(isubd) + ceiling(real(num_element_types(itype),rp) / real(VECTOR_SIZE_LOC,rp), ip )
           end if
        end do
     end if     
     call memory_alloca(par_memor,'LIST_BOUNDARIES_NORACE_PAR % PACKS','par_boundary_loop',list_boundaries_norace_par(isubd) % packs,num_pack_norace_nboun_par(isubd))
     
     if( adjust_vector_size ) then
        ipack = 0
        do itype = 1,num_types
           if( num_element_types(itype) > 0 ) then
              ipack = ipack + 1
              call memory_alloca(par_memor,'LIST_BOUNDARIES_NORACE_PAR % PACKS % L','par_boundary_loop',list_boundaries_norace_par(isubd) % packs(ipack) % l,num_element_types(itype))
           end if
        end do
     else
        do ipack = 1,num_pack_norace_nboun_par(isubd)
           call memory_alloca(par_memor,'LIST_BOUNDARIES_NORACE_PAR % PACKS % L','par_boundary_loop',list_boundaries_norace_par(isubd) % packs(ipack) % l,VECTOR_SIZE_LOC) 
        end do
     end if
          
     ivect = 0
     ipack = 0 
     jtype = 0
     do itype = 1,num_types
        if( num_element_types(itype) > 0 ) then
           do iboun = 1,nboun
              if( ltypb(iboun) > 0 .and. element_types(iboun) == itype ) then
                 ivect = ivect + 1
                 if( adjust_vector_size ) then
                    if( ivect > num_element_types(itype) .or. itype /= jtype ) then                    
                       ivect = 1
                       ipack = ipack + 1
                    end if
                 else
                    if( ivect > VECTOR_SIZE_LOC .or. itype /= jtype ) then                    
                       ivect = 1
                       ipack = ipack + 1
                    end if
                 end if
                 if( ipack > size(list_boundaries_norace_par(isubd) % packs) ) &
                      call runend('WRONG PACK NUMBER: CALL GUILLAUME TO FIX IT')
                 list_boundaries_norace_par(isubd) % packs(ipack) % l(ivect) = iboun
                 jtype = itype
              end if
           end do
        end if
     end do
     !
     ! Number of subdomains NUM_SUBD_NBOUN_PAR
     !
     if(      par_hybrid == PAR_HYBRID_OFF ) then
        num_subd_nboun_par = 1
     else if( par_hybrid == PAR_OPENMP_NO_COLORING ) then
        num_subd_nboun_par = 1
     else if( par_hybrid == PAR_OPENMP_COLORING ) then
        num_subd_nboun_par = par_omp_nboun_num_colors
     else if( par_hybrid == PAR_OMPSS ) then
        if( associated(ompss_boundaries) ) num_subd_nboun_par = size(ompss_boundaries)           
     end if

     allocate( list_boundaries_par(num_subd_nboun_par) )
     do isubd = 1,num_subd_nboun_par
        nullify(list_boundaries_par(isubd) % packs)
     end do
     call memory_alloca(par_memor,'NUM_PACK_NBOUN_PAR','par_boundary_loop',num_pack_nboun_par,num_subd_nboun_par)

  end if

  if( par_hybrid == PAR_HYBRID_OFF .or. par_hybrid == PAR_OPENMP_NO_COLORING ) then

     !-------------------------------------------------------------------
     !
     ! No hybrid parallelization 
     !
     !-------------------------------------------------------------------

     if( par_hybrid /= PAR_HYBRID_OFF ) then
        call messages_live('LOOP WITH RACE CONDITION: OPENMP WITHOUT COLORING')
     end if

     if( INOTMASTER ) then

        isubd = 1
        num_element_types = 0
        do iboun = 1,nboun
           if( ltypb(iboun) > 0 ) then
              itype = element_types(iboun)
              num_element_types(itype) = num_element_types(itype) + 1
           end if
        end do
        
        if( adjust_vector_size ) then
           num_pack_nboun_par(isubd) = 0 
           do itype = 1,num_types
              if( num_element_types(itype) > 0 ) then
                 num_pack_nboun_par(isubd) = num_pack_nboun_par(isubd) + 1
              end if
           end do
        else
           num_pack_nboun_par(isubd) = 0 
           do itype = 1,num_types
              if( num_element_types(itype) > 0 ) then
                 num_pack_nboun_par(isubd) = num_pack_nboun_par(isubd) + ceiling(real(num_element_types(itype),rp) / real(VECTOR_SIZE_LOC,rp), ip )
              end if
           end do
        end if
        
        call memory_alloca(par_memor,'LIST_BOUNDARIES_PAR % PACKS','par_boundary_loop',list_boundaries_par(isubd) % packs,num_pack_nboun_par(isubd))
        if( adjust_vector_size ) then
           ipack = 0
           do itype = 1,num_types
              if( num_element_types(itype) > 0 ) then
                 ipack = ipack + 1
                 call memory_alloca(par_memor,'LIST_BOUNDARIES_PAR % PACKS % L','par_boundary_loop',list_boundaries_par(isubd) % packs(ipack) % l,num_element_types(itype))
              end if
           end do
        else
           do ipack = 1,num_pack_nboun_par(isubd)
              call memory_alloca(par_memor,'LIST_BOUNDARIES_PAR % PACKS % L','par_boundary_loop',list_boundaries_par(isubd) % packs(ipack) % l,VECTOR_SIZE_LOC) 
           end do
        end if
        
        ivect = 0
        ipack = 0 
        jtype = 0
        do itype = 1,num_types
           if( num_element_types(itype) > 0 ) then
              do iboun = 1,nboun
                 if( ltypb(iboun) > 0 .and. element_types(iboun) == itype ) then
                    ivect = ivect + 1
                    if( adjust_vector_size ) then
                       if( ivect > num_element_types(itype) .or. itype /= jtype ) then                    
                          ivect = 1
                          ipack = ipack + 1
                       end if
                    else
                       if( ivect > VECTOR_SIZE_LOC .or. itype /= jtype ) then                    
                          ivect = 1
                          ipack = ipack + 1
                       end if
                    end if
                    if( ipack > size(list_boundaries_par(isubd) % packs) ) call runend('WRONG PACK NUMBER')
                    list_boundaries_par(isubd) % packs(ipack) % l(ivect) = iboun
                    jtype = itype
                 end if
              end do
           end if
        end do

     end if

  else if( par_hybrid == PAR_OPENMP_COLORING  ) then

     !-------------------------------------------------------------------
     !
     ! OpenMP with coloring
     !
     !-------------------------------------------------------------------

     call messages_live('LOOP WITH RACE CONDITION: OPENMP WITH COLORING')

     if( INOTMASTER ) then

        do isubd = 1,num_subd_nboun_par
           num_element_types = 0
           do kelem = par_omp_nboun_ia_colors(isubd),par_omp_nboun_ia_colors(isubd+1)-1
              iboun = par_omp_nboun_ja_colors(kelem)
              if( ltypb(iboun) > 0 ) then
                 itype = element_types(iboun)
                 num_element_types(itype) = num_element_types(itype) + 1
              end if
           end do

           if( adjust_vector_size ) then
              num_pack_nboun_par(isubd) = num_types
           else
              num_pack_nboun_par(isubd) = 0
              do itype = 1,num_types
                 if( num_element_types(itype) > 0 ) then
                    num_pack_nboun_par(isubd) = num_pack_nboun_par(isubd) + ceiling(real(num_element_types(itype),rp) / real(VECTOR_SIZE_LOC,rp), ip )
                 end if
              end do
           end if
           
           call memory_alloca(par_memor,'LIST_BOUNDARIES_PAR % PACKS','par_boundary_loop',list_boundaries_par(isubd) % packs,num_pack_nboun_par(isubd))
           if( adjust_vector_size ) then
              do ipack = 1,num_pack_nboun_par(isubd)
                 itype = ipack
                 call memory_alloca(par_memor,'LIST_BOUNDARIES_PAR % PACKS % L','par_boundary_loop',list_boundaries_par(isubd) % packs(ipack) % l,num_element_types(itype))
              end do
           else
              do ipack = 1,num_pack_nboun_par(isubd)
                 call memory_alloca(par_memor,'LIST_BOUNDARIES_PAR % PACKS % L','par_boundary_loop',list_boundaries_par(isubd) % packs(ipack) % l,VECTOR_SIZE_LOC) 
              end do
           end if
           ivect = 0
           ipack = 0 
           jtype = 0

           do itype = 1,num_types
              if( num_element_types(itype) > 0 ) then
                 do kelem = par_omp_nboun_ia_colors(isubd),par_omp_nboun_ia_colors(isubd+1)-1
                    iboun = par_omp_nboun_ja_colors(kelem)
                    if( ltypb(iboun) > 0 .and. element_types(iboun) == itype ) then
                       ivect = ivect + 1
                       if( adjust_vector_size ) then
                          if( ivect > num_element_types(itype) .or. itype /= jtype ) then                    
                             ivect = 1
                             ipack = ipack + 1
                          end if
                       else
                          if( ivect > VECTOR_SIZE_LOC .or. itype /= jtype ) then                    
                             ivect = 1
                             ipack = ipack + 1
                          end if
                       end if
                       if( ipack > size(list_boundaries_par(isubd) % packs) ) call runend('WRONG PACK NUMBER')
                       list_boundaries_par(isubd) % packs(ipack) % l(ivect) = iboun
                       jtype = itype
                    end if
                 end do
              end if
           end do
        end do
     end if

  else if( par_hybrid == PAR_OMPSS ) then

     !-------------------------------------------------------------------
     !
     ! OmpSs
     !
     !-------------------------------------------------------------------

     call messages_live('LOOP WITH RACE CONDITION: OMPSS WITH MULDIDEPENDENCIES')

     if( INOTMASTER .and. nboun > 0 ) then
        do isubd = 1,num_subd_nboun_par
           num_element_types = 0
           do kelem = 1,size(ompss_boundaries(isubd) % elements)
              iboun = ompss_boundaries(isubd) % elements(kelem)
              if( ltypb(iboun) > 0 ) then
                 itype = element_types(iboun)
                 num_element_types(itype) = num_element_types(itype) + 1
              end if
           end do
           num_pack_nboun_par(isubd) = 0 
           do itype = 1,num_types
              if( num_element_types(itype) > 0 ) then
                 num_pack_nboun_par(isubd) = num_pack_nboun_par(isubd) + ceiling(real(num_element_types(itype),rp) / real(VECTOR_SIZE_LOC,rp), ip )
              end if
           end do
           call memory_alloca(par_memor,'LIST_BOUNDARIES_PAR % PACKS','par_boundary_loop',list_boundaries_par(isubd) % packs,num_pack_nboun_par(isubd))
           do ipack = 1,num_pack_nboun_par(isubd)
              call memory_alloca(par_memor,'LIST_BOUNDARIES_PAR % PACKS','par_boundary_loop',list_boundaries_par(isubd) % packs(ipack) % l,VECTOR_SIZE_LOC) 
           end do

           ivect = 0
           ipack = 0 
           jtype = 0
           do itype = 1,num_types
              if( num_element_types(itype) > 0 ) then
                 do kelem = 1,size(ompss_boundaries(isubd) % elements)
                    iboun = ompss_boundaries(isubd) % elements(kelem)
                    if( ltypb(iboun) > 0 .and. element_types(iboun) == itype ) then
                       ivect = ivect + 1
                       if( ivect > VECTOR_SIZE_LOC .or. itype /= jtype ) then                    
                          ivect = 1
                          ipack = ipack + 1
                       end if
                       if( ipack > size(list_boundaries_par(isubd) % packs) ) call runend('WRONG PACK NUMBER')
                       list_boundaries_par(isubd) % packs(ipack) % l(ivect) = iboun
                       jtype = itype
                    end if
                 end do
              end if
           end do
        end do

     end if

  end if
  !
  ! Check
  !
  if( INOTMASTER .and. nboun > 0 ) then

     if( 1 == 1 ) then
        call memory_alloca(par_memor,'CONSIDER_ELEMENT','par_boundary_loop',consider_element  ,nboun)
        do isubd = 1,num_subd_nboun_par
           do ipack = 1,num_pack_nboun_par(isubd)
              do kelem = 1,memory_size(list_boundaries_par(isubd) % packs(ipack) % l)
                 iboun = list_boundaries_par(isubd) % packs(ipack) % l(kelem)
                 if( iboun > 0 ) then
                    if( ltypb(iboun) > 0 ) then
                       if( consider_element(iboun)) call runend('ELEMENT ALREADY TOUCHED')
                       consider_element(iboun) = .true.
                    end if
                 end if
              end do
           end do
        end do
        kelem = 0
        inode = 0
        do iboun = 1,nboun
           if( .not. consider_element(iboun) .and. ltypb(iboun) > 0 ) kelem=kelem+1
           if( consider_element(iboun) ) inode = inode + 1
        end do
        if( inode /= nboun ) then
           print*,'error=',kfl_paral,nboun,inode
           call runend('WE ARE IN TROUBLE')
        end if
        if(kelem>0) call runend('ERROR: '//trim(intost(kelem))//' ELEMENTS NOT TOUCHED')
        do isubd = 1,num_subd_nboun_par
           do ipack = 1,num_pack_nboun_par(isubd)
              iboun = list_boundaries_par(isubd) % packs(ipack) % l(1)
              itype = element_types(iboun)
              do kelem = 2,memory_size(list_boundaries_par(isubd) % packs(ipack) % l)
                 iboun = list_boundaries_par(isubd) % packs(ipack) % l(kelem)
                 if( iboun > 0 ) then
                    if( element_types(iboun) /= itype ) then
                       print*,'elem 1= ',list_boundaries_par(isubd) % packs(ipack) % l(1),itype
                       print*,'elem 2= ',iboun,element_types(iboun)
                       call runend('ERROR: ELEMENT OF DIFFERENT TYPES IN THE SAME LIST')
                    end if
                 end if
              end do
           end do
        end do
        call memory_deallo(par_memor,'CONSIDER_ELEMENT','par_boundary_loop',consider_element)
     end if
     !
     ! Deallocate memory
     !
     call memory_deallo(par_memor,'ELEMENT_TYPES'    ,'par_boundary_loop',element_types)
     call memory_deallo(par_memor,'NUM_ELEMENT_TYPES','par_boundary_loop',num_element_types)

  end if

  call messages_live('HYBRID PARALLELIZATION & VECTORIZATION (BOUNDARIES)','END SECTION')

end subroutine par_boundary_loop
  
end module mod_par_element_loop
!> @} 
