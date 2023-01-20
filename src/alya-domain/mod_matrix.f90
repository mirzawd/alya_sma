!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!>
!> @defgroup Matrix_Toolbox
!> Toolbox for matrix operations
!> @{
!> @file    mod_matrix.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!> @brief   Definition and subroutines of MatrixToolBox
!>
!------------------------------------------------------------------------

module mod_matrix
  
  use def_kintyp,            only : ip,rp,lg,i1p,r1p
  use def_spmat,             only : spmat
  use mod_memory,            only : memory_alloca,memory_deallo
  use mod_memory,            only : memory_resize,memory_copy
  use def_master,            only : INOTMASTER,NPOIN_TYPE,ISLAVE,IPARALL,ISEQUEN
  use def_master,            only : kfl_paral,intost
  use mod_graphs,            only : graphs_find_edge
  use mod_graphs,            only : graphs_number_to_linked_list
  use mod_optional_argument, only : optional_argument
  use mod_communications,    only : PAR_MIN,PAR_MAX,PAR_SUM

  implicit none
  !
  ! CHUNK size for matrix vector product
  !
  integer(ip)              :: spmv_chunk       = 1000
  integer(ip),   parameter :: SOL_OMP_STATIC   = 2
  integer(ip),   parameter :: SOL_OMP_GUIDED   = 3
  integer(ip),   parameter :: SOL_OMP_DYNAMIC  = 4

  private

  interface matrix_initialize
     module procedure matrix_initialize_rp_1,&
          &           matrix_initialize_rp_2,&
          &           matrix_initialize_rp_3
  end interface matrix_initialize

  interface matrix_COO_SpMV_coupling
     module procedure matrix_COO_SpMV_BASE,&
          &           matrix_COO_SpMV_BLOCK,&
          &           matrix_COO_SpMV_SPMAT
  end interface matrix_COO_SpMV_coupling

  public :: matrix_initialize
  public :: matrix_asscsr
  public :: matrix_asscsr_all
  public :: matrix_assrhs
  public :: matrix_assexp
  public :: matrix_blokgs
  public :: matrix_invdia

  public :: matrix_diagonal_CSR                               ! Diagonal of a matrix in CSR format
  public :: matrix_diagonal_COO                               ! Diagonal of a matrix in COO format
  public :: matrix_diagonal_ELL                               ! Diagonal of a matrix in ELL format

  public :: matrix_assemble_element_matrix_to_CSR             ! Assemble element matrix in CSR format
  public :: matrix_assemble_element_matrix_to_COO             ! Assemble element matrix in COO format
  public :: matrix_assemble_element_matrix_to_ELL             ! Assemble element matrix in ELL format
  public :: matrix_LHS_TO_RHS                                 ! cast LHS to RHS: RHS <= RHS - LHS * U 
  public :: matrix_assemble_2by2_block_element_matrix_to_CSR
  public :: matrix_assemble_element_RHS
  public :: matrix_assemble_2by2_block_element_RHS
  
  public :: matrix_scaling_CSR                                ! Scale a matrix from left or right in CSR format
  public :: matrix_add_diagonal_CSR                           ! Add a diagonal matrix to a CSR matrix
  public :: matrix_add_CSR                                    ! Add a CSR matrix to a CSR matrix

  public :: matrix_matgr2
  public :: matrix_output_gid_format                          ! Output a matrix in GiD format
  public :: matrix_output_dense_format                        ! Output a matrix in dense format
  public :: matrix_output_pajek_net_format                    ! Output a matrix in Pajek net format
  public :: matrix_permute_and_copy_matrix
  !public :: assresdiff
  public :: matrix_copy_matrix
  public :: matrix_copy_matrix_block
  public :: matrix_copy_matrix_to_block
  public :: matrix_CSR_SpMV
  public :: matrix_COO_SpMV
  public :: matrix_ELL_SpMV
  public :: matrix_CSR_parallel_SpMV
  public :: matrix_COO_SPGEMM
  public :: matrix_COO_aggregate
  public :: matrix_full_to_half
  public :: matrix_csr_to_coo
  public :: matrix_chksym
  public :: nullify_spmat
  public :: matrix_transpose
  public :: matrix_partial_to_full_row_matrix
  public :: matrix_copy
  public :: matrix_heap_sort
  public :: matrix_remove_null_coefficients
  public :: matrix_sparsify                       ! Sparsify a matrix
  public :: matrix_rotate_system                  ! Rotate the row of a matrix and rhs
  public :: matrix_CSR_rotate
  
  public :: matrix_COO_SpMV_coupling              ! To be unifified with matrix_COO_SpMV
  
contains

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    26/09/2014
  !> @brief   Initialize a matrix (or array)
  !> @details Initialize a matrix (or array)
  !
  !----------------------------------------------------------------------

  subroutine matrix_initialize_rp_1(aa)
    real(rp),   pointer :: aa(:)
    integer(ip)         :: i1
    do i1 = 1,size(aa,KIND=ip)
       aa(i1) = 0.0_rp
    end do
  end subroutine matrix_initialize_rp_1
  subroutine matrix_initialize_rp_2(aa)
    real(rp),   pointer :: aa(:,:)
    integer(ip)         :: i1,i2
    do i2 = 1,size(aa,2,KIND=ip)
       do i1 = 1,size(aa,1,KIND=ip)
          aa(i1,i2) = 0.0_rp
       end do
    end do
  end subroutine matrix_initialize_rp_2
  subroutine matrix_initialize_rp_3(aa)
    real(rp),   pointer :: aa(:,:,:)
    integer(ip)         :: i1,i2,i3
    do i3 = 1,size(aa,3,KIND=ip)
       do i2 = 1,size(aa,2,KIND=ip)
          do i1 = 1,size(aa,1,KIND=ip)
             aa(i1,i2,i3) = 0.0_rp
          end do
       end do
    end do
  end subroutine matrix_initialize_rp_3
  !----------------------------------------------------------------------
  !
  !> @author  Ricard Borrell
  !> @date    15/07/2016
  !> @brief   Nullify pointers inside spmat
  !> @details Nullify pointers inside spmat
  !
  !----------------------------------------------------------------------
  subroutine nullify_spmat(A)

    implicit none
    type(spmat), intent(inout) :: A

    A % ndof1 = 0
    A % ndof2 = 0
    A % nrows = 0
    A % ncols = 0
    nullify(A % iA)
    nullify(A % jA)
    nullify(A % vA)

  end subroutine nullify_spmat

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Agglomerate a matrix
  !> @details Agglomerate a matrix
  !
  !----------------------------------------------------------------------

  subroutine matrix_matgr2(&
       npoin,ngrou,ndofn,lgrou,ia,ja,iagroup,jagroup,aa,aagroup,memit)
    implicit none
    integer(ip),           intent(in)    :: npoin               !< Number of nodes
    integer(ip),           intent(in)    :: ngrou               !< Number of groups
    integer(ip),           intent(in)    :: ndofn               !< Nb of dof
    integer(ip), pointer,  intent(in)    :: lgrou(:)            !< List of groups
    integer(ip),           intent(in)    :: ia(*)
    integer(ip),           intent(in)    :: ja(*)
    integer(ip), pointer,  intent(in)    :: iagroup(:)
    integer(ip), pointer,  intent(in)    :: jagroup(:)
    real(rp),              intent(in)    :: aa(ndofn,ndofn,*)
    real(rp),    pointer,  intent(inout) :: aagroup(:,:,:)
    integer(8),  optional, intent(inout) :: memit(2)            !< Memory counter
    integer(ip)                          :: idofn,izdom,izgro
    integer(ip)                          :: nzgro,igrou,jgrou
    integer(ip)                          :: jdofn,ipoin,jpoin
    integer(8)                           :: memor(2)
    logical(lg), pointer                 :: lchek(:) => null()
    !
    ! Allocate matrix if necessary
    !
    if( .not. associated(aagroup) ) then
       nzgro = size(jagroup,KIND=ip)
       if( present(memit) ) then
          call memory_alloca(memit,'AAGROUP','matrix_matgro',aagroup,ndofn,ndofn,nzgro)
       else
          memor = 0
          call memory_alloca(memor,'AAGROUP','matrix_matgro',aagroup,ndofn,ndofn,nzgro)
       end if
    end if
    !
    ! Check array
    !
    allocate( lchek(ngrou) )
    do igrou = 1,ngrou
       lchek(igrou) = .true.
    end do
    !
    ! Agglomerate matrix AAGROUP <= AA
    !
    if( ndofn == 1 ) then

       do ipoin = 1,npoin
          if( lgrou(ipoin) > 0 ) then
             igrou = lgrou(ipoin)
             do izdom = ia(ipoin),ia(ipoin+1)-1
                jpoin = ja(izdom)
                if( lgrou(jpoin) > 0 ) then
                   jgrou = lgrou(jpoin)
                   izgro = iagroup(igrou)
                   iifzgro1: do while( izgro <= iagroup(igrou+1)-1 )
                      if( jagroup(izgro) == jgrou ) exit iifzgro1
                      izgro = izgro + 1
                   end do iifzgro1
                   lchek(igrou) = .false.
                   aagroup(1,1,izgro) = aagroup(1,1,izgro) + aa(1,1,izdom)
                end if
             end do
          end if
       end do

    else

       do ipoin = 1,npoin
          if( lgrou(ipoin) > 0 ) then
             igrou = lgrou(ipoin)
             do izdom = ia(ipoin),ia(ipoin+1)-1
                jpoin = ja(izdom)
                if( lgrou(jpoin) > 0 ) then
                   jgrou = lgrou(jpoin)
                   izgro = iagroup(igrou)
                   iifzgro2: do while( izgro <= iagroup(igrou+1)-1 )
                      if( jagroup(izgro) == jgrou ) exit iifzgro2
                      izgro = izgro + 1
                   end do iifzgro2
                   lchek(igrou) = .false.
                   do jdofn = 1,ndofn
                      do idofn = 1,ndofn
                         aagroup(idofn,jdofn,izgro) = aagroup(idofn,jdofn,izgro) &
                              + aa(idofn,jdofn,izdom)
                      end do
                   end do
                end if
             end do
          end if
       end do

    end if
    !
    ! Check if group has a non-zero on its row
    !
    do igrou = 1,ngrou
       if( lchek(igrou) ) then
          izgro = iagroup(igrou)
          iifzgro3: do while( izgro <= iagroup(igrou)+1-1 )
             if( jagroup(izgro) == igrou ) exit iifzgro3
             izgro = izgro + 1
          end do iifzgro3
          do idofn = 1,ndofn
             aagroup(idofn,idofn,izgro) = 1.0_rp
          end do
       end if
    end do
    deallocate(lchek)
    !
    ! Output matrix in ps format
    !
    !call pspltm(&
    !     ngrou,ngrou,1_ip,0_ip,jagroup,iagroup,aagroup,&
    !     'caca',0_ip,18.0_rp,'cm',&
    !     0_ip,0_ip,2_ip,99)

  end subroutine matrix_matgr2


  !-----------------------------------------------------------------------
  !
  !> @brief   Check the symmetry of a matrix
  !> @details Check the symmetry of a matrix
  !> @date    04/07/2012
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine matrix_chksym(nn,nbvar,ia,ja,aa,xsymm,memit)

    integer(ip),           intent(in)    :: nn
    integer(ip),           intent(in)    :: nbvar
    integer(ip),           intent(in)    :: ia(*)
    integer(ip),           intent(in)    :: ja(*)
    real(rp),              intent(in)    :: aa(nbvar,nbvar,*)
    real(rp),              intent(out)   :: xsymm
    integer(8),  optional, intent(inout) :: memit(2)            !< Memory counter
    integer(ip)                          :: nz,ii,kk,ll,jj,qq,mm
    real(rp),    pointer                 :: bb(:,:)
    real(rp)                             :: xx(nbvar)
    integer(8)                           :: memor(2)

    nullify(bb)
    nz = (ia(nn+1)-1)*nbvar*nbvar
    if( present(memit) ) then
       call memory_alloca(memit,'BB','matrix_chksym',bb,nbvar,nz)
    else
       memor = 0
       call memory_alloca(memor,'BB','matrix_chksym',bb,nbvar,nz)
    end if

    do ii = 1,nn
       do kk = ia(ii),ia(ii+1)-1
          jj = ja(kk)
          if( jj /= ii ) then
             loop_ll: do ll = ia(jj),ia(jj+1)-1
                if( ja(ll) == ii ) then
                   if( abs(aa(1,1,kk)-aa(1,1,ll)) > 1.0e-12_rp ) then
                      bb(1,kk) = abs(aa(1,1,kk)-aa(1,1,ll))&
                           &        /(0.5_rp*abs(aa(1,1,kk))+abs(aa(1,1,ll)))
                      bb(1,ll) = bb(1,kk)
                      if( nbvar > 1 ) then
                         do qq = 2,nbvar
                            if( abs(aa(qq,qq,kk))+abs(aa(qq,qq,ll)) > 1.0e-12_rp ) then
                               bb(qq,kk) = abs(aa(qq,qq,kk)-aa(qq,qq,ll))&
                                    &         /(0.5_rp*abs(aa(qq,qq,kk))+abs(aa(qq,qq,ll)))
                               bb(qq,ll) = bb(qq,kk)
                            end if
                         end do
                      end if
                   end if
                   exit loop_ll
                end if
             end do loop_ll
          end if
       end do
    end do
    !
    ! Average symmetry
    !
    mm = 0
    xx = 0.0_rp
    do ii = 1,nn
       do kk = ia(ii),ia(ii+1)-1
          jj = ja(kk)
          if( jj /= ii ) then
             mm = mm + 1
             do qq = 1,nbvar
                if( kk> size(bb,KIND=ip))stop
                xx(qq) = xx(qq) + bb(qq,kk)
             end do
          end if
       end do
    end do
    xsymm = 0.0_rp
    do qq = 1,nbvar
       xsymm  = xsymm + xx(qq) / real(mm,rp)
    end do
    xsymm = xsymm / real(nbvar,rp)
    if( present(memit) ) then
       call memory_deallo(memit,'BB','matrix_chksym',bb)
    else
       call memory_deallo(memor,'BB','matrix_chksym',bb)
    end if

  end subroutine matrix_chksym

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Block Gauss-Seidel
  !> @details Block Gauss-Seidel
  !
  !> To come back to the previous version (nblok=ndofn) just put ndonf_n_per_block = 1
  !
  !----------------------------------------------------------------------

  subroutine matrix_blokgs(&
       itask,nequa,ndofn_per_block,nblok,kblok,ia,ja,aa,bb,xx,&
       aa_aux,bb_idofn,xx_idofn,memit,lperm)
    implicit none
    integer(ip),          intent(in)    :: itask                                   !< Small system or solution update
    integer(ip),          intent(in)    :: nequa                                   !< Number of nodes
    integer(ip),          intent(in)    :: ndofn_per_block                         !< Number of dof per node in the block
    integer(ip),          intent(in)    :: nblok                                   !< Nb of dof
    integer(ip),          intent(in)    :: kblok                                   !< Dof to be solved
    integer(ip),          intent(in)    :: ia(*)                                   !< CSR format
    integer(ip),          intent(in)    :: ja(*)                                   !< CSR format
    real(rp),             intent(in)    :: aa(ndofn_per_block*nblok,ndofn_per_block*nblok,*)   !< Original matrix
    real(rp),             intent(in)    :: bb(ndofn_per_block*nblok,*)                   !< Original RHS
    real(rp),             intent(inout) :: xx(ndofn_per_block*nblok,*)                   !< Original solution
    real(rp),    pointer                :: aa_idofn(:,:,:)                             !< IDOFN Block system matrix
    real(rp),    pointer, intent(inout)   :: aa_aux(:)
    real(rp),    pointer, intent(inout)   :: bb_idofn(:)                             !< IDOFN Block system RHS
    real(rp),    pointer, intent(inout)   :: xx_idofn(:)                             !< IDOFN Block system solution
    integer(8),  optional,intent(inout) :: memit(2)                                !< Memory counter
    integer(ip)                         :: iequa,iz,jequa,cont,jcont
    integer(ip)                         :: nz,jdofn,icont
    integer(8)                          :: memor(2)
    integer(ip), pointer, optional      :: lperm(:,:)                              !< permutation vector

    nullify(aa_idofn)
    
    select case ( itask )

    case ( 1_ip )
       !
       ! Allocate if necessary
       !
       nz = ia(nequa+1)-1
       if( present( memit ) ) then
          if( .not. associated(aa_idofn) ) call memory_alloca(memit,'AA_IDOFN','matrix_blokgs',aa_idofn,ndofn_per_block,ndofn_per_block,nz)
          if( .not. associated(xx_idofn) ) call memory_alloca(memit,'XX_IDOFN','matrix_blokgs',xx_idofn,ndofn_per_block*nequa)
          if( .not. associated(bb_idofn) ) call memory_alloca(memit,'BB_IDOFN','matrix_blokgs',bb_idofn,ndofn_per_block*nequa)
       else
          memor = 0
          if( .not. associated(aa_idofn) ) call memory_alloca(memor,'AA_IDOFN','matrix_blokgs',aa_idofn,ndofn_per_block,ndofn_per_block,nz)
          if( .not. associated(xx_idofn) ) call memory_alloca(memor,'XX_IDOFN','matrix_blokgs',xx_idofn,ndofn_per_block*nequa)
          if( .not. associated(bb_idofn) ) call memory_alloca(memor,'BB_IDOFN','matrix_blokgs',bb_idofn,ndofn_per_block*nequa)
       end if
       !
       ! Construct IDOFN system: AA_IDOFN(:) XX_IDOFN(:) = BB_IDOFN(:)
       !


       if( present( lperm ) ) then

       else
          do iequa = 1,nequa  !loop on points (nodes)
             do jdofn = 1,ndofn_per_block
                xx_idofn((iequa-1)*ndofn_per_block + jdofn) = xx((kblok-1)*ndofn_per_block+jdofn, iequa)
                bb_idofn((iequa-1)*ndofn_per_block + jdofn) = bb((kblok-1)*ndofn_per_block+jdofn, iequa)
             end do
             do iz = ia(iequa), ia(iequa+1)-1
                jequa = ja(iz)
                do icont=1,ndofn_per_block
                   do jdofn=1,ndofn_per_block
                      aa_idofn(icont,jdofn,iz) = aa((kblok-1)*ndofn_per_block + icont,(kblok-1)*ndofn_per_block+jdofn,iz)
                   end do
                end do
                do jdofn = 1,nblok
                   if( jdofn /= kblok) then
                      do icont = 1,ndofn_per_block
                         do jcont = 1,ndofn_per_block
                            bb_idofn((iequa-1)*ndofn_per_block+icont) = bb_idofn((iequa-1)*ndofn_per_block+icont) &
                                 - aa((jdofn-1)*ndofn_per_block+icont,(kblok-1)*ndofn_per_block+jcont,iz) * xx((jdofn-1)*ndofn_per_block+jcont,jequa)
                         end do
                      end do
                   end if
                end do
             end do
          end do
       end if


       ! reshape matrix aa_idofn into a vector
       call memory_alloca(memit,'AA_AUX','matrix_blokgs',aa_aux,ndofn_per_block*ndofn_per_block*nz)
       do iz=1,nz
          cont = 1
          do iequa=1,ndofn_per_block
             do jdofn=1,ndofn_per_block
                aa_aux(cont) = aa_idofn(iequa,jdofn,iz)
                cont = cont+1
             end do
          end do
       end do


    case ( 2_ip )
       !
       ! Update unknown in global array XX(IDOFN,:) <= XX(:)
       !
       if( present( lperm ) ) then
       else
          do iequa=1,nequa
             do icont = 1,ndofn_per_block
                xx(icont,iequa) = xx_idofn((iequa-1)*ndofn_per_block+icont)
             end do
          end do
       end if


    case ( 3_ip )
       !
       ! Deallocate memory
       !
       if( present( memit ) ) then
          call memory_deallo(memit,'AA_IDOFN','matrix_blokgs',aa_idofn)
          call memory_deallo(memit,'XX_IDOFN','matrix_blokgs',xx_idofn)
          call memory_deallo(memit,'BB_IDOFN','matrix_blokgs',bb_idofn)
       else
          call memory_deallo(memor,'AA_IDOFN','matrix_blokgs',aa_idofn)
          call memory_deallo(memor,'XX_IDOFN','matrix_blokgs',xx_idofn)
          call memory_deallo(memor,'BB_IDOFN','matrix_blokgs',bb_idofn)
       end if

    end select

  end subroutine matrix_blokgs

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Block Gauss-Seidel
  !> @details Block Gauss-Seidel
  !
  !----------------------------------------------------------------------

  !subroutine matrix_blokgs(&
  !     itask,nequa,ndofn,idofn,ia,ja,aa,bb,xx,&
  !     aa_idofn,bb_idofn,xx_idofn,memit)
  !  implicit none
  !  integer(ip),          intent(in)    :: itask               !< Small system or solution update
  !  integer(ip),          intent(in)    :: nequa               !< Number of nodes
  !  integer(ip),          intent(in)    :: ndofn               !< Nb of dof
  !  integer(ip),          intent(in)    :: idofn               !< Dof to be solved
  !  integer(ip),          intent(in)    :: ia(*)               !< CSR format
  !  integer(ip),          intent(in)    :: ja(*)               !< CSR format
  !  real(rp),             intent(in)    :: aa(ndofn,ndofn,*)   !< Original matrix
  !  real(rp),             intent(in)    :: bb(ndofn,*)         !< Original RHS
  !  real(rp),             intent(inout) :: xx(ndofn,*)         !< Original solution
  !  real(rp),    pointer, intent(inout)   :: aa_idofn(:)         !< IDOFN Block system matrix
  !  real(rp),    pointer, intent(inout)   :: bb_idofn(:)         !< IDOFN Block system RHS
  !  real(rp),    pointer, intent(inout)   :: xx_idofn(:)         !< IDOFN Block system solution
  !  integer(8),  optional,intent(inout) :: memit(2)            !< Memory counter
  !  integer(ip)                         :: iequa,iz,jequa
  !  integer(ip)                         :: nz,jdofn
  !  integer(8)                          :: memor(2)

  !  select case ( itask )

  !  case ( 1_ip )
  !
  ! Allocate if necessary
  !
  !     nz = ia(nequa+1)-1
  !     if( present( memit ) ) then
  !        if( .not. associated(aa_idofn) ) call memory_alloca(memit,'AA_IDOFN','matrix_blokgs',aa_idofn,nz)
  !        if( .not. associated(xx_idofn) ) call memory_alloca(memit,'XX_IDOFN','matrix_blokgs',xx_idofn,nequa)
  !        if( .not. associated(bb_idofn) ) call memory_alloca(memit,'BB_IDOFN','matrix_blokgs',bb_idofn,nequa)
  !     else
  !        memor = 0
  !        if( .not. associated(aa_idofn) ) call memory_alloca(memor,'AA_IDOFN','matrix_blokgs',aa_idofn,nz)
  !        if( .not. associated(xx_idofn) ) call memory_alloca(memor,'XX_IDOFN','matrix_blokgs',xx_idofn,nequa)
  !        if( .not. associated(bb_idofn) ) call memory_alloca(memor,'BB_IDOFN','matrix_blokgs',bb_idofn,nequa)
  !     end if
  !     !
  !     ! Construct IDOFN system: AA_IDOFN(:) XX_IDOFN(:) = BB_IDOFN(:)
  !     !
  !     do iequa = 1,nequa
  !        xx_idofn(iequa) = xx(idofn,iequa)
  !        bb_idofn(iequa) = bb(idofn,iequa)
  !        do iz = ia(iequa),ia(iequa+1)-1
  !           jequa = ja(iz)
  !           aa_idofn(iz) = aa(idofn,idofn,iz)
  !           do jdofn = 1,ndofn
  !              if( jdofn /= idofn ) &
  !                   bb_idofn(iequa) = bb_idofn(iequa) - aa(jdofn,idofn,iz) * xx(jdofn,jequa)
  !           end do
  !        end do
  !     end do

  !  case ( 2_ip )
  !     !
  !     ! Update unknown in global array XX(IDOFN,:) <= XX(:)
  !     !
  !     do iequa = 1,nequa
  !        xx(idofn,iequa) = xx_idofn(iequa)
  !     end do

  !  case ( 3_ip )
  !
  ! Deallocate memory
  !
  !     if( present( memit ) ) then
  !        call memory_deallo(memit,'AA_IDOFN','matrix_blokgs',aa_idofn)
  !        call memory_deallo(memit,'XX_IDOFN','matrix_blokgs',xx_idofn)
  !        call memory_deallo(memit,'BB_IDOFN','matrix_blokgs',bb_idofn)
  !     else
  !        call memory_deallo(memor,'AA_IDOFN','matrix_blokgs',aa_idofn)
  !        call memory_deallo(memor,'XX_IDOFN','matrix_blokgs',xx_idofn)
  !        call memory_deallo(memor,'BB_IDOFN','matrix_blokgs',bb_idofn)
  !     end if

  !  end select

  !end subroutine matrix_blokgs_ECR

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Agglomerate a matrix
  !> @details Agglomerate a matrix
  !
  !----------------------------------------------------------------------


  subroutine matrix_parsum(ndimr,rvarr)
    implicit none
    integer(ip),  intent(in) :: ndimr
    real(rp),     target     :: rvarr(ndimr)
    if(IPARALL) then
      call PAR_SUM(ndimr,rvarr)
    end if
  end subroutine matrix_parsum

  subroutine matrix_parslx(ndimr,rvarr)
    use mod_periodicity_sequential
    use def_domain, only: nperi, npoin
    use def_master, only: npari, nparc, pard1
    implicit none
    integer(ip),  intent(in) :: ndimr
    real(rp),     target     :: rvarr(ndimr)
    npari = 0
    nparc = 0
    if( ISEQUEN .and. nperi /= 0) then
       pard1 = ndimr/npoin
       call periodicity_sequential(pard1,rvarr)
    end if
  end subroutine matrix_parslx

  subroutine matrix_matgro(&
       itask,npoin,ngrou,ndofn,kfl_assem_group,kfl_symme,lgrou,&
       ia,ja,iagroup,jagroup,nsizegroup,aa,aagroup)
    implicit none
    integer(ip),          intent(in)    :: itask                  !< Perform an all reduce
    integer(ip),          intent(in)    :: npoin                  !< Number of nodes
    integer(ip),          intent(in)    :: ngrou                  !< Number of groups
    integer(ip),          intent(in)    :: ndofn                  !< Nb of dof
    integer(ip),          intent(in)    :: kfl_assem_group        !< Assembly type of coarse matrix (0: skyline, 1=CSR)
    integer(ip),          intent(in)    :: kfl_symme              !< Symmetric assembly (1) or not (0)
    integer(ip), pointer, intent(in)    :: lgrou(:)               !< List of groups
    integer(ip),          intent(in)    :: ia(*)                  !< Matrix graph
    integer(ip),          intent(in)    :: ja(*)                  !< Matrix graph
    integer(ip), pointer, intent(in)    :: iagroup(:)             !< Skyline: iskyl, CSR: iagrou
    integer(ip), pointer, intent(in)    :: jagroup(:)             !< Skyline: null, CSR: jagrou
    integer(ip),          intent(in)    :: nsizegroup             !< Size of the coarse matrix
    real(rp),             intent(in)    :: aa(ndofn,ndofn,*)      !< Fine matrix
    real(rp),             intent(inout) :: aagroup(ndofn,ndofn,*) !< Coarse matrix
    integer(ip)                         :: idofn,izdom,izgro
    integer(ip)                         :: igrou,jgrou
    integer(ip)                         :: jdofn,ipoin,jpoin
    integer(ip)                         :: kskyl
    logical(lg), pointer                :: lchek(:)
    integer(ip), pointer                :: iskyl(:)

    nullify(iskyl)
    nullify(lchek)
    !
    ! Check array
    !
    allocate( lchek(ngrou) )
    do igrou = 1,ngrou
       lchek(igrou) = .true.
    end do
    !
    ! Agglomerate matrix AAGROUP <= AA
    !
    if( kfl_assem_group ==  0 ) then
       !
       ! Skyline format
       !
       iskyl => iagroup

       if( ndofn == 1 ) then

          if( kfl_symme == 0 ) then
             !
             ! Skyline format: A(1,1,:) is not symmetric
             !
             do ipoin = 1,npoin
                if( lgrou(ipoin) > 0 ) then
                   igrou = lgrou(ipoin)
                   do izdom = ia(ipoin),ia(ipoin+1)-1
                      jpoin = ja(izdom)
                      if( jpoin < ipoin ) then
                         if( lgrou(jpoin) > 0 ) then
                            jgrou = lgrou(jpoin)
                            if( igrou > jgrou ) then
                               kskyl              = iskyl(igrou+1) - 1 - (igrou-jgrou)
                               aagroup(1,1,kskyl) = aagroup(1,1,kskyl) + aa(1,1,izdom)
                            else if( igrou < jgrou ) then
                               kskyl              = iskyl(jgrou+1) - 1 - (jgrou-igrou)
                               aagroup(1,1,kskyl) = aagroup(1,1,kskyl) + aa(1,1,izdom)
                            else
                               kskyl              = iskyl(igrou+1) - 1
                               aagroup(1,1,kskyl) = aagroup(1,1,kskyl) + 2.0_rp*aa(1,1,izdom)
                            end if
                         end if
                      else if( ipoin == jpoin ) then
                         kskyl              = iskyl(igrou+1) - 1
                         aagroup(1,1,kskyl) = aagroup(1,1,kskyl) + aa(1,1,izdom)
                      end if
                   end do
                end if
             end do

          else
             !
             ! Skyline format: A(1,1,:) is symmetric
             !
             do ipoin = 1,npoin
                if( lgrou(ipoin) > 0 ) then
                   igrou = lgrou(ipoin)
                   do izdom = ia(ipoin),ia(ipoin+1)-2
                      jpoin = ja(izdom)
                      if( lgrou(jpoin) > 0 ) then
                         jgrou = lgrou(jpoin)
                         if( igrou > jgrou ) then
                            kskyl              = iskyl(igrou+1)-1-(igrou-jgrou)
                            aagroup(1,1,kskyl) = aagroup(1,1,kskyl)+aa(1,1,izdom)
                         else if( igrou < jgrou ) then
                            kskyl              = iskyl(jgrou+1)-1-(jgrou-igrou)
                            aagroup(1,1,kskyl) = aagroup(1,1,kskyl)+aa(1,1,izdom)
                         else
                            kskyl              = iskyl(igrou+1)-1
                            aagroup(1,1,kskyl) = aagroup(1,1,kskyl)+2.0_rp*aa(1,1,izdom)
                         end if
                      end if
                   end do
                   izdom              = ia(ipoin+1)-1
                   kskyl              = iskyl(igrou+1)-1
                   aagroup(1,1,kskyl) = aagroup(1,1,kskyl) + aa(1,1,izdom)
                end if
             end do


          end if

       else

          call runend('MATRIX_MATGRO: FORMAT NOT CODED')

       end if

    else if( kfl_assem_group ==  1 ) then
       !
       ! CSR format
       !
       if( ndofn == 1 ) then

          if( kfl_symme == 0 ) then
             !
             ! CSR format: A(1,1,:) is not symmetric
             !
             do ipoin = 1,npoin
                if( lgrou(ipoin) > 0 ) then
                   igrou = lgrou(ipoin)
                   do izdom = ia(ipoin),ia(ipoin+1)-1
                      jpoin = ja(izdom)
                      if( lgrou(jpoin) > 0 ) then
                         jgrou = lgrou(jpoin)
                         izgro = iagroup(igrou)
                         iifzgro1: do while( izgro <= iagroup(igrou+1)-1 )
                            if( jagroup(izgro) == jgrou ) exit iifzgro1
                            izgro = izgro + 1
                         end do iifzgro1
                         lchek(igrou) = .false.
                         aagroup(1,1,izgro) = aagroup(1,1,izgro) + aa(1,1,izdom)
                      end if
                   end do
                end if
             end do

          else
             !
             ! CSR format: A(1,1,:) is symmetric
             !
             call runend('MATRIX_MATGRO: NOT CODED')

          end if

       else

          if( kfl_symme == 0 ) then
             !
             ! CSR format: A(NDOFN,NDOFN,1,:) is not symmetric
             !
             do ipoin = 1,npoin
                if( lgrou(ipoin) > 0 ) then
                   igrou = lgrou(ipoin)
                   do izdom = ia(ipoin),ia(ipoin+1)-1
                      jpoin = ja(izdom)
                      if( lgrou(jpoin) > 0 ) then
                         jgrou = lgrou(jpoin)
                         izgro = iagroup(igrou)
                         iifzgro3: do while( izgro <= iagroup(igrou+1)-1 )
                            if( jagroup(izgro) == jgrou ) exit iifzgro3
                            izgro = izgro + 1
                         end do iifzgro3
                         lchek(igrou) = .false.
                         do jdofn = 1,ndofn
                            do idofn = 1,ndofn
                               aagroup(idofn,jdofn,izgro) = aagroup(idofn,jdofn,izgro) &
                                    + aa(idofn,jdofn,izdom)
                            end do
                         end do
                      end if
                   end do
                end if
             end do

          else
             !
             ! CSR format: A(NDOFN,NDOFN,1,:) is symmetric
             !
             call runend('MATRIX_MATGRO: NOT CODED')

          end if

       end if

    end if
    !
    ! Check if group has a non-zero on its row
    !
    do igrou = 1,ngrou
       if( lchek(igrou) ) then
          izgro = iagroup(igrou)
          iifzgro4: do while( izgro <= iagroup(igrou)+1-1 )
             if( jagroup(izgro) == igrou ) exit iifzgro4
             izgro = izgro + 1
          end do iifzgro4
          do idofn = 1,ndofn
             aagroup(idofn,idofn,izgro) = 1.0_rp
          end do
       end if
    end do
    deallocate(lchek)
    !
    ! Perform an all reduce
    !
    if( itask == 1 ) then
      call matrix_parsum(nsizegroup,aagroup)
    end if
    !
    ! Output matrix in ps format
    !
    !call pspltm(&
    !     ngrou,ngrou,1_ip,0_ip,jagroup,iagroup,aagroup,&
    !     'caca',0_ip,18.0_rp,'cm',&
    !     0_ip,0_ip,2_ip,99)

  end subroutine matrix_matgro

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Assemble a RHS
  !> @details Assemble a RHS
  !
  !----------------------------------------------------------------------

  subroutine matrix_assrhs(ndofn_glo,ndofn_ele,pnode,npoin,lnods,elrhs,rhsid,kdofn)
    implicit none
    integer(ip),  intent(in)           :: ndofn_glo
    integer(ip),  intent(in)           :: ndofn_ele
    integer(ip),  intent(in)           :: pnode
    integer(ip),  intent(in)           :: npoin
    integer(ip),  intent(in)           :: lnods(pnode)
    real(rp),     intent(in)           :: elrhs(*)
    real(rp),     intent(inout)        :: rhsid(*)
    integer(ip),  intent(in), optional :: kdofn
    integer(ip)                        :: inode,ipoin,idofl,idofg,idofn,ndof2

    if( present(kdofn) ) then
       !
       ! NDOFN_GLO DOF's but only KDOFN is assembled
       !
       if( ndofn_ele == 1 ) then
          do inode = 1,pnode
             ipoin = lnods(inode)
             idofg = (ipoin-1) * ndofn_glo + kdofn
             rhsid(idofg) = rhsid(idofg) + elrhs(inode)
          end do
       else
          call runend('MATRIX_ASSRHS: VERY STRANGE RHS ASSEMBLY INDEED...')
       end if

    else

       if(ndofn_glo==1) then
          !
          ! 1 DOF
          !
          do inode=1,pnode
             ipoin=lnods(inode)
             rhsid(ipoin)=rhsid(ipoin)+elrhs(inode)
          end do

       else if(ndofn_glo==2) then
          !
          ! 2 DOF's
          !
          do inode=1,pnode
             ipoin=lnods(inode)
             idofg=2*ipoin-1
             idofl=2*inode-1
             rhsid(idofg)=rhsid(idofg)+elrhs(idofl)
             idofg=idofg+1
             idofl=idofl+1
             rhsid(idofg)=rhsid(idofg)+elrhs(idofl)
          end do

       else if(ndofn_glo>2) then
          !
          ! >2 DOF's
          !
          do inode=1,pnode
             ipoin=lnods(inode)
             idofg=(ipoin-1)*ndofn_glo
             idofl=(inode-1)*ndofn_glo
             do idofn=1,ndofn_glo
                idofg=idofg+1
                idofl=idofl+1
                rhsid(idofg)=rhsid(idofg)+elrhs(idofl)
             end do
          end do

       else if(ndofn_glo<0) then
          !
          ! >2 DOF's
          !
          ndof2=abs(ndofn_glo)
          do inode=1,pnode
             ipoin=lnods(inode)
             idofg=(ipoin-1)*ndof2
             idofl=(inode-1)*ndof2
             do idofn=1,ndof2
                idofg=(idofn-1)*npoin+ipoin
                idofl=idofl+1
                rhsid(idofg)=rhsid(idofg)+elrhs(idofl)
             end do
          end do

       end if
    end if

  end subroutine matrix_assrhs


  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Assemble a RHS
  !> @details Assemble a RHS
  !
  !----------------------------------------------------------------------

  subroutine matrix_assexp(ndofn_glo,ndofn_ele,pnode,npoin,lnods,elrhs,elmat,elunk,rhsid,kdofn)
    implicit none
    integer(ip),  intent(in)           :: ndofn_glo                              !< Number of dof of global RHS
    integer(ip),  intent(in)           :: ndofn_ele                              !< Number of dof of element matrix and RHS
    integer(ip),  intent(in)           :: pnode                                  !< Number of nodes of element
    integer(ip),  intent(in)           :: npoin                                  !< Number of nodes of the mesh
    integer(ip),  intent(in)           :: lnods(pnode)                           !< Element connectivity
    real(rp),     intent(in)           :: elrhs(*)                               !< Element RHS
    real(rp),     intent(in)           :: elmat(pnode*ndofn_ele,pnode*ndofn_ele) !< Element matrix
    real(rp),     intent(in)           :: elunk(pnode*ndofn_ele)                 !< Element unknwon
    real(rp),     intent(inout)        :: rhsid(*)                               !< Global RHS
    integer(ip),  intent(in), optional :: kdofn                                  !< Local dof to assemble
    integer(ip)                        :: inode,ipoin,idofg
    integer(ip)                        :: jnode

    if( ndofn_glo > 1 .and. ndofn_ele == 1 .and. ( .not. present(kdofn) ) ) then
       !
       ! NDOFN_GLO DOF's but only KDOFN is assembled
       !
       call runend('MATRIX_ASSEXP: KDOFN IS MISSING')

    else

       if( ndofn_glo == 1 ) then
          !
          ! 1 DOF 
          !
          do inode = 1,pnode
             ipoin = lnods(inode)
             rhsid(ipoin) = rhsid(ipoin) + elrhs(inode)
             do jnode = 1,pnode
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
                rhsid(ipoin) = rhsid(ipoin) - elmat(inode,jnode) * elunk(jnode)
             end do
          end do

       else if( ndofn_glo > 1 ) then
          !
          ! DOF > 1
          !
          if( ndofn_ele == 1 ) then
             do inode = 1,pnode
                ipoin = lnods(inode)
                idofg = (ipoin-1)*ndofn_glo+kdofn
                rhsid(idofg) = rhsid(idofg) + elrhs(inode)
                do jnode = 1,pnode
#ifdef NO_COLORING
                   !$OMP ATOMIC
#endif
                   rhsid(idofg) = rhsid(idofg) - elmat(inode,jnode) * elunk(jnode)
                end do
             end do
          else
             call runend('NOT CODED, ESPECIE DE VAGO')
          end if

       end if
    end if

  end subroutine matrix_assexp

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Assemble a BCSR matrix
  !> @details Assemble a BCSR matrix
  !
  !----------------------------------------------------------------------

  subroutine matrix_asscsr(ndofn_glo,ndofn_ele,ia,ja,pnode,lnods,elmat,amatr,kdofn)
    implicit none
    integer(ip), intent(in)           :: ndofn_glo
    integer(ip), intent(in)           :: ndofn_ele
    integer(ip), intent(in)           :: ia(*)
    integer(ip), intent(in)           :: ja(*)
    integer(ip), intent(in)           :: pnode
    integer(ip), intent(in)           :: lnods(pnode)
    real(rp),    intent(in)           :: elmat(ndofn_ele*pnode,ndofn_ele*pnode)
    real(rp),    intent(inout)        :: amatr(ndofn_glo,ndofn_glo,*)
    integer(ip), intent(in), optional :: kdofn
    integer(ip)                       :: ievat,jevat,idofn,jdofn,inode,jnode
    integer(ip)                       :: ipoin,jpoin,izsol,jcolu

    if( present(kdofn) ) then
       !
       ! NDOFN unknowns but assemble only KDOFN
       !
       if( ndofn_ele == 1 ) then

          do inode = 1,pnode
             ipoin = lnods(inode)
             do jnode = 1,pnode
                jpoin = lnods(jnode)
                izsol = ia(ipoin)
                jcolu = ja(izsol)
                do while( jcolu /= jpoin .and. izsol < ia(ipoin+1)-1 )
                   izsol = izsol + 1
                   jcolu = ja(izsol)
                end do
                if( jpoin == jcolu ) then
                   idofn = kdofn
                   jdofn = kdofn
                   ievat = (inode-1) * ndofn_glo + idofn
                   jevat = (jnode-1) * ndofn_glo + jdofn
#ifdef NO_COLORING
                   !$OMP ATOMIC
#endif
                   amatr(jdofn,idofn,izsol) = amatr(jdofn,idofn,izsol) + elmat(inode,jnode)
                end if
             end do
          end do

       else

          call runend('MATRIX_ASSCSR: VERY STRANGE CASE INDEED...')

       end if

    else

       if( ndofn_glo == 1 ) then
          !
          ! 1 unknown
          !
          do inode = 1,pnode
             ipoin = lnods(inode)
             do jnode = 1,pnode
                jpoin = lnods(jnode)
                izsol = ia(ipoin)
                jcolu = ja(izsol)
                do while( jcolu /= jpoin .and. izsol < ia(ipoin+1)-1)
                   izsol = izsol + 1
                   jcolu = ja(izsol)
                end do
                if( jcolu == jpoin ) then
#ifdef NO_COLORING
                   !$OMP ATOMIC
#endif
                   amatr(1,1,izsol) = amatr(1,1,izsol) + elmat(inode,jnode)
                end if
             end do
          end do

       else
          !
          ! NDOFN unknown
          !
          do inode = 1,pnode
             ipoin = lnods(inode)
             do jnode=1,pnode
                jpoin = lnods(jnode)
                izsol = ia(ipoin)
                jcolu = ja(izsol)
                do while(jcolu/=jpoin .and. izsol < ia(ipoin+1)-1 )
                   izsol = izsol + 1
                   jcolu = ja(izsol)
                end do
                if( jpoin == jcolu ) then
                   do idofn = 1,ndofn_glo
                      ievat = (inode-1) * ndofn_glo + idofn
                      do jdofn = 1,ndofn_glo
                         jevat = (jnode-1) * ndofn_glo + jdofn
#ifdef NO_COLORING
                         !$OMP ATOMIC
#endif
                         amatr(jdofn,idofn,izsol) = amatr(jdofn,idofn,izsol) + elmat(ievat,jevat)
                      end do
                   end do
                end if
             end do
          end do
       end if

    end if

  end subroutine matrix_asscsr

  !----------------------------------------------------------------------
  !
  !> @author  Mohammad Kouhi
  !> @date    01/01/2016
  !> @brief   Assemble a BCSR matrix with two dofns
  !> @details Assemble a BCSR matrix with two dofns
  !
  !----------------------------------------------------------------------

  subroutine matrix_asscsr_all(ndofn_glo,ndofn_ele,ia,ja,pnode,lnods,elmat,amatr,kdofn,ldofn)
    implicit none
    integer(ip), intent(in)           :: ndofn_glo
    integer(ip), intent(in)           :: ndofn_ele
    integer(ip), intent(in)           :: ia(*)
    integer(ip), intent(in)           :: ja(*)
    integer(ip), intent(in)           :: pnode
    integer(ip), intent(in)           :: lnods(pnode)
    real(rp),    intent(in)           :: elmat(ndofn_ele*pnode,ndofn_ele*pnode)
    real(rp),    intent(inout)        :: amatr(ndofn_glo,ndofn_glo,*)
    integer(ip), intent(in), optional :: kdofn
    integer(ip), intent(in), optional :: ldofn
    integer(ip)                       :: ievat,jevat,idofn,jdofn,inode,jnode
    integer(ip)                       :: ipoin,jpoin,izsol,jcolu

    if( present(kdofn) ) then
       if( present(ldofn) ) then
          !
          ! NDOFN unknowns but assemble KDOFN and LDOFN
          !
          if( ndofn_ele == 1 ) then

             do inode = 1,pnode
                ipoin = lnods(inode)
                do jnode = 1,pnode
                   jpoin = lnods(jnode)
                   izsol = ia(ipoin)
                   jcolu = ja(izsol)
                   do while( jcolu /= jpoin .and. izsol < ia(ipoin+1)-1 )
                      izsol = izsol + 1
                      jcolu = ja(izsol)
                   end do
                   if( jpoin == jcolu ) then
                      idofn = kdofn
                      jdofn = ldofn
                      amatr(jdofn,idofn,izsol) = amatr(jdofn,idofn,izsol) + elmat(inode,jnode) ! elmat(ievat,jevat) WHY THIS IEVAT IF THE MATRIX ONLY BRINGS IPOINS
                   end if
                end do
             end do
          else
             call runend('MATRIX_ASSCSR: VERY STRANGE CASE INDEED...')
          end if
       else
          !
          ! NDOFN unknowns but assemble only KDOFN
          !
          if( ndofn_ele == 1 ) then

             do inode = 1,pnode
                ipoin = lnods(inode)
                do jnode = 1,pnode
                   jpoin = lnods(jnode)
                   izsol = ia(ipoin)
                   jcolu = ja(izsol)
                   do while( jcolu /= jpoin .and. izsol < ia(ipoin+1)-1 )
                      izsol = izsol + 1
                      jcolu = ja(izsol)
                   end do
                   if( jpoin == jcolu ) then
                      idofn = kdofn
                      jdofn = kdofn
                      ievat = (inode-1) * ndofn_glo + idofn
                      jevat = (jnode-1) * ndofn_glo + jdofn
                      amatr(jdofn,idofn,izsol) = amatr(jdofn,idofn,izsol) + elmat(inode,jnode) ! elmat(ievat,jevat) WHY THIS IEVAT IF THE MATRIX ONLY BRINGS IPOINS
                   end if
                end do
             end do

          else

             call runend('MATRIX_ASSCSR: VERY STRANGE CASE INDEED...')

          end if
       endif
    else

       if( ndofn_glo == 1 ) then
          !
          ! 1 unknown
          !
          do inode = 1,pnode
             ipoin = lnods(inode)
             do jnode = 1,pnode
                jpoin = lnods(jnode)
                izsol = ia(ipoin)
                jcolu = ja(izsol)
                do while( jcolu /= jpoin .and. izsol < ia(ipoin+1)-1)
                   izsol = izsol + 1
                   jcolu = ja(izsol)
                end do
                if( jcolu == jpoin ) then
                   amatr(1,1,izsol) = amatr(1,1,izsol) + elmat(inode,jnode)
                end if
             end do
          end do

       else
          !
          ! NDOFN unknown
          !
          do inode = 1,pnode
             ipoin = lnods(inode)
             do jnode=1,pnode
                jpoin = lnods(jnode)
                izsol = ia(ipoin)
                jcolu = ja(izsol)
                do while(jcolu/=jpoin .and. izsol < ia(ipoin+1)-1 )
                   izsol = izsol + 1
                   jcolu = ja(izsol)
                end do
                if( jpoin == jcolu ) then
                   do idofn = 1,ndofn_glo
                      ievat = (inode-1) * ndofn_glo + idofn
                      do jdofn = 1,ndofn_glo
                         jevat = (jnode-1) * ndofn_glo + jdofn
                         amatr(jdofn,idofn,izsol) = amatr(jdofn,idofn,izsol) + elmat(ievat,jevat)
                      end do
                   end do
                end if
             end do
          end do
       end if

    end if

  end subroutine matrix_asscsr_all

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Compute the diagonal of a matrix
  !> @details Computer the diagonal of a matrix inCSR format
  !
  !----------------------------------------------------------------------

  subroutine matrix_diagonal_CSR(nn,ndof,kfl_symme,ia,ja,an,diagonal)

    integer(ip),  intent(in)          :: nn
    integer(ip),  intent(in)          :: ndof
    integer(ip),  intent(in)          :: kfl_symme
    integer(ip),  intent(in), pointer :: ia(:)
    integer(ip),  intent(in), pointer :: ja(:)
    real(rp),     intent(in)          :: an(ndof,ndof,*)
    real(rp),     intent(out)         :: diagonal(*)
    integer(ip)                       :: ii,kk,ll,jj

    !--------------------------------------------------------------
    !
    ! CSR format
    !
    !--------------------------------------------------------------

    if( kfl_symme == 1 ) then
       !
       ! Symmetric graph
       !
       if( ndof == 1 ) then
          do ii = 1, nn
             ll = ia(ii+1)-1
             diagonal(ii) = an(1,1,ll)
          end do
       else
          do ii= 1, nn
             ll = ia(ii+1)-1
             jj = (ii-1) * ndof
             do kk= 1, ndof
                diagonal(jj+kk) = an(kk,kk,ll)
             end do
          end do
       end if

    else
       !
       ! Unsymmetric graph
       !
       if( ndof == 1 ) then
          do ii= 1, nn
             jj = ia(ii)
             ll = -1
             do while( jj < ia(ii+1) .and. ll ==-1 )
                if( ja(jj) == ii ) then
                   ll = jj
                end if
                jj = jj+1
             end do
             if( ll /= -1 ) then
                diagonal(ii) = an(1,1,ll)
             else
                diagonal(ii) = 0.0_rp
             end if
          end do
       else
          do ii= 1, nn
             jj = ia(ii)
             ll = -1
             do while( jj< ia(ii+1) .and. ll ==-1 )
                if( ja(jj) == ii ) then
                   ll = jj
                end if
                jj = jj + 1
             end do
             if( ll /= -1 ) then
                jj = (ii-1) * ndof
                do kk= 1, ndof
                   diagonal(jj+kk) = an(kk,kk,ll)
                end do
             else
                jj = (ii-1) * ndof
                do kk= 1, ndof
                   diagonal(jj+kk) = 0.0_rp
                end do
             end if
          end do
       end if

    end if

  end subroutine matrix_diagonal_CSR

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Compute the diagonal of a matrix in COO format
  !> @details Compute the diagonal of a matrix in COO format
  !
  !----------------------------------------------------------------------

  subroutine matrix_diagonal_COO(nn,ndof,rows,cols,an,diagonal)

    integer(ip),  intent(in)          :: nn
    integer(ip),  intent(in)          :: ndof
    integer(ip),  intent(in), pointer :: cols(:)
    integer(ip),  intent(in), pointer :: rows(:)
    real(rp),     intent(in)          :: an(ndof,ndof,*)
    real(rp),     intent(out)         :: diagonal(*)
    integer(ip)                       :: ii,kk,iz,nz

    if( size(cols,KIND=ip) /= size(rows,KIND=ip) ) &
         call runend('MATRIX_DIAGONAL_COO: SOMETHING IS STRANGE ABOUT COO FORMAT')
    nz = size(cols,KIND=ip)

    do iz = 1,nz
       if( rows(iz) == cols(iz) ) then
          ii = (rows(iz)-1) * ndof
          do kk = 1,ndof
             diagonal(ii+kk) = an(kk,kk,iz)
          end do
       end if
    end do

  end subroutine matrix_diagonal_COO

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Compute the diagonal of a matrix in ELL format
  !> @details Compute the diagonal of a matrix in ELL format
  !
  !----------------------------------------------------------------------

  subroutine matrix_diagonal_ELL(nn,ndof,cols_ell,an,diagonal)

    integer(ip),  intent(in)          :: nn
    integer(ip),  intent(in)          :: ndof
    integer(ip),  intent(in), pointer :: cols_ell(:,:)
    real(rp),     intent(in)          :: an(ndof,ndof,*)
    real(rp),     intent(out)         :: diagonal(*)
    integer(ip)                       :: ii,iz,ncols

    ncols = size(cols_ell,1,KIND=ip)

    if( ndof == 1 ) then
       do ii = 1,nn
          iz = 1
          do while( cols_ell(iz,ii) /= ii )
             iz = iz + 1
          end do
          iz = (ii-1)*ncols+iz
          diagonal(ii) = an(1,1,iz)
       end do
    else
       call runend('DIAGONAL FOR ELL FORMAT NOT COMPUTED')
    end if

  end subroutine matrix_diagonal_ELL

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Compute the inverse diagonal of a matrix
  !> @details Assemble a BCSR matrix
  !
  !----------------------------------------------------------------------

  subroutine matrix_invdia(npoin,nbvar,kfl_symme,ia,ja,an,invdiag,memit,what)
    integer(ip), intent(in)             :: npoin,nbvar,kfl_symme
    integer(ip), intent(in)             :: ia(*),ja(*)
    real(rp),    intent(in)             :: an(nbvar,nbvar,*)
    real(rp),    intent(inout), pointer :: invdiag(:)
    integer(8),  intent(inout)          :: memit(2)            !< Memory counter
    character(*),intent(in),   optional :: what                !< Options
    integer(ip)                         :: ii,jj,kk,ll
    logical(lg)                         :: exchange

    exchange = .true.
    if( present(what) ) then
       if( trim(what) == 'DO NOT EXCHANGE' .or. trim(what) == 'DO NOT ASSEMBLE' ) then
          exchange = .false.
       end if
    end if

    if( INOTMASTER ) then

       if( .not. associated(invdiag) ) then
          call memory_alloca(memit,'INVDIAG','matrix_invdia',invdiag,npoin*nbvar)
       end if

       if( kfl_symme == 1 ) then
          !
          ! Symmetric graph
          !
          if( nbvar == 1 ) then
             do ii= 1, npoin
                ll = ia(ii+1)-1
                invdiag(ii) = an(1,1,ll)
             end do
          else
             do ii= 1, npoin
                ll = ia(ii+1)-1
                jj = (ii-1) * nbvar
                do kk= 1, nbvar
                   invdiag(jj+kk) = an(kk,kk,ll)
                end do
             end do
          end if

       else
          !
          ! Unsymmetric graph
          !
          if( nbvar == 1 ) then
             do ii= 1, npoin
                jj = ia(ii)
                ll = -1
                do while (jj< ia(ii+1) .and. ll ==-1)
                   if(ja(jj)==ii) then
                      ll = jj
                   end if
                   jj = jj+1
                end do
                if(ll/=-1) then
                   invdiag(ii)=an(1,1,ll)
                else
                   invdiag(ii)=0.0_rp
                end if
             end do
          else
             do ii= 1, npoin
                jj = ia(ii)
                ll = -1
                do while (jj< ia(ii+1) .and. ll ==-1)
                   if(ja(jj)==ii) then
                      ll = jj
                   end if
                   jj = jj+1
                end do
                if(ll/=-1) then
                   jj = (ii-1) * nbvar
                   do kk= 1, nbvar
                      invdiag(jj+kk)=an(kk,kk,ll)
                   end do
                else
                   jj = (ii-1) * nbvar
                   do kk= 1, nbvar
                      invdiag(jj+kk)=0.0_rp
                   end do
                end if
             end do
          end if

       end if
       !
       ! Periodicity and Parallelization
       !
       if( exchange ) call matrix_parslx(nbvar*npoin,invdiag)
       !
       ! Inverse diagonal
       !
       do ii = 1,npoin*nbvar
          if( invdiag(ii) /= 0.0_rp ) then
             invdiag(ii) = 1.0_rp / invdiag(ii)
          end if
       end do

    end if

  end subroutine matrix_invdia

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Diagonal of a Schur complement
  !> @details Diagonal of a Schur complement
  !
  !----------------------------------------------------------------------

  subroutine matrix_diagonal_schur(nbnodes,nbvar,ndofn_A3,A1,A2,invA3,A4,ia,ja,diag)
    implicit none
    integer(ip), intent(in)          :: nbnodes,nbvar,ndofn_A3
    integer(ip), intent(in)          :: ia(*)
    integer(ip), intent(in)          :: ja(*)
    real(rp),    intent(in)          :: A1(nbvar,nbvar,*)
    real(rp),    intent(in)          :: A2(ndofn_A3,*)
    real(rp),    intent(in)          :: invA3(ndofn_A3,nbnodes)
    real(rp),    intent(in)          :: A4(ndofn_A3,*)
    real(rp),    intent(out)         :: diag(*)
    integer(ip)                      :: izdom,kzdom,lzdom,ii
    integer(ip)                      :: jj,kk,mm,idofn

    if( INOTMASTER ) then

       if( nbvar == 1 ) then

          do ii = 1,nbnodes
             do izdom = ia(ii),ia(ii+1)-1
                jj = ja(izdom)
                if( ii == jj ) then
                   diag(ii) = A1(1,1,izdom)
                   do kzdom = ia(ii),ia(ii+1)-1
                      kk    = ja(kzdom)
                      lzdom = ia(kk)
                      do while( lzdom <= ia(kk+1)-1 )
                         mm = ja(lzdom)
                         if( mm == jj ) then
                            do idofn = 1,ndofn_A3
                               diag(ii) = diag(ii) - A2(idofn,kzdom) * A4(idofn,lzdom) * invA3(idofn,kk)
                            end do
                            lzdom = ia(kk+1)
                         end if
                         lzdom = lzdom + 1
                      end do
                   end do
                end if
             end do
          end do

       end if

    end if

  end subroutine matrix_diagonal_schur

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    11/10/2017
  !> @brief   Output a matrix in dense format
  !> @details Output a matrix in a dense format
  !
  !----------------------------------------------------------------------

  subroutine matrix_output_dense_format(npoin,nbvar,ia,ja,an)
    integer(ip), intent(in)           :: npoin             !< Number of nodes
    integer(ip), intent(in)           :: nbvar             !< Number of degrees of freedom per node
    integer(ip), intent(in)           :: ia(*)             !< CSR graph
    integer(ip), intent(in)           :: ja(*)             !< CSR graph
    real(rp),    intent(in)           :: an(nbvar,nbvar,*) !< Matrix
    integer(ip)                       :: ii,jj,kk,ll
    integer(ip)                       :: idofn,jdofn
    integer(4)                        :: iunit

    iunit = 98
    open( unit = iunit , file = 'matrix_'//trim(intost(kfl_paral))//'.dat' , status = 'unknown' )

    do ii = 1,npoin
       do idofn = 1,nbvar
          do jj = 1,npoin
             kk = 0
             ll = ia(ii)-1
             do while( ll < ia(ii+1)-1 .and. kk == 0 )
                ll = ll + 1
                if( ja(ll) == jj ) kk = ll
             end do
             do jdofn = 1,nbvar
                if( kk /= 0 ) then
                   write(iunit,'(1x,e12.6,a)',advance='no') an(jdofn,idofn,kk),' ,'
                else
                   write(iunit,'(1x,e12.6,a)',advance='no') 0.0_rp,' ,'
                end if
             end do
          end do
          write(iunit,*)
       end do
    end do

    close(iunit)

  end subroutine matrix_output_dense_format

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    11/10/2017
  !> @brief   Output a matrix in PAJEK NET format
  !> @details Output a matrix in PAJEK NET format
  !
  !----------------------------------------------------------------------

  subroutine matrix_output_pajek_net_format(npoin,nbvar,ia,ja,an,ABSOLUTE_VALUE)
    
    integer(ip), intent(in)           :: npoin             !< Number of nodes
    integer(ip), intent(in)           :: nbvar             !< Number of degrees of freedom per node
    integer(ip), intent(in)           :: ia(*)             !< CSR graph
    integer(ip), intent(in)           :: ja(*)             !< CSR graph
    real(rp),    intent(in)           :: an(nbvar,nbvar,*) !< Matrix
    logical(lg), intent(in), optional :: ABSOLUTE_VALUE    !< Use absolute value of matrix
    integer(ip)                       :: ii,iz
    integer(4)                        :: iunit
    logical(lg)                       :: use_absolute_value

    if( nbvar > 1 ) call runend('MATRIX_OUTPUT_NET_FORMAT: NOT CODED')
    iunit = 98
    open( unit = iunit , file = 'matrix_'//trim(intost(abs(kfl_paral)))//'.net' , status = 'unknown' )

    use_absolute_value = .false.
    if( present(ABSOLUTE_VALUE) ) use_absolute_value = ABSOLUTE_VALUE
    
    write(iunit,'(a,1x,i7)') '*Vertices ',npoin
    write(iunit,'(a)')       '*arcs'
    
    do ii = 1,npoin
       do iz = ia(ii),ia(ii+1)-1
          if( abs(an(1,1,iz)) > 1.0e-12_rp ) then
             if( use_absolute_value ) then
                write(iunit,'(2(1x,i6),1x,e12.6,a)',advance='no') ii,ja(iz),abs(an(1,1,iz))
             else
                write(iunit,'(2(1x,i6),1x,e12.6,a)',advance='no') ii,ja(iz),an(1,1,iz)
             end if
          end if
       end do
    end do

    close(iunit)

  end subroutine matrix_output_pajek_net_format

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    29/06/2015
  !> @brief   Output a matrix in GiD format
  !> @details Output a matrix in GiD format
  !
  !----------------------------------------------------------------------

  subroutine matrix_output_gid_format(npoin,nbvar,ia,ja,an,invpr)
    integer(ip), intent(in)           :: npoin
    integer(ip), intent(in)           :: nbvar
    integer(ip), intent(in)           :: ia(*)
    integer(ip), intent(in)           :: ja(*)
    real(rp),    intent(in)           :: an(nbvar,nbvar,*)
    integer(ip), intent(in), optional :: invpr(*)
    !integer(ip), intent(in), optional :: lgrou(*)
    integer(ip)                       :: ii,jj,kk,ke,ie,je
    integer(ip)                       :: nx,ny,ex,ey,izdom,ipoin
    integer(ip)                       :: ii_old,jj_old
    real(rp)                          :: xx,yy
    logical(lg)                       :: notfound

    open( unit = 98 , file = 'matrix_'//trim(intost(kfl_paral))//'.post.msh' , status = 'unknown' )
    open( unit = 99 , file = 'matrix_'//trim(intost(kfl_paral))//'.post.res' , status = 'unknown' )

    nx = npoin + 1
    ny = npoin + 1
    ex = npoin
    ey = npoin
    !
    ! Mesh
    !
    write(98,'(a)') 'MESH MATRIX dimension 2 Elemtype Quadrilateral Nnode 4'
    kk = 0
    write(98,*) 'coordinates'
    do jj = 1,nx
       yy = real(ny-(jj-1)*nbvar,rp)
       do ii = 1,ny
          xx = real((ii-1)*nbvar,rp)
          kk = kk + 1
          write(98,*) kk,xx,yy
       end do
    end do
    write(98,*) 'end coordinates'

    ke = 0
    write(98,*) 'elements'
    do je = 1,ey
       do ie = 1,ex
          ke = ke + 1
          ii = nx*(je-1) + ie
          write(98,*) ke,ii,ii+1,ii+nx+1,ii+nx
       end do
    end do
    write(98,*) 'end elements'
    !
    ! Matrix coefficients
    !
    write(99,*) 'GiD Post Results File 1.0'
    write(99,*) ' '
    write(99,*) 'GaussPoints GP Elemtype Quadrilateral'
    write(99,*) 'Number of Gauss Points: 1'
    write(99,*) 'Natural Coordinates: Internal'
    write(99,*) 'End GaussPoints'
    write(99,*) 'Result MATRIX ALYA  0.00000000E+00 Scalar OnGaussPoints GP'
    write(99,*) 'ComponentNames MATRIX'
    write(99,*) 'Values'

    if( .not. present(invpr) ) then
       do ii = 1,npoin
          do jj = 1,npoin
             izdom    = ia(ii)-1
             notfound = .true.
             do while( notfound .and. izdom < ia(ii+1)-1 )
                izdom = izdom + 1
                if( ja(izdom) == jj ) notfound = .false.
             end do
             ipoin = (ii-1)*ex + jj
             if( notfound ) then
                write(99,*) ipoin,0.0_rp
             else
                write(99,*) ipoin,abs(an(1,1,izdom))
             end if
          end do
       end do
    else
       do ii = 1,npoin
          ii_old = invpr(ii)
          do jj = 1,npoin
             jj_old   = invpr(jj)
             izdom    = ia(ii_old)-1
             notfound = .true.
             do while( notfound .and. izdom < ia(ii_old+1)-1 )
                izdom = izdom + 1
                if( ja(izdom) == jj_old ) notfound = .false.
             end do
             ipoin = (ii-1)*ex + jj
             if( .not. notfound ) then
                write(99,*) ipoin,abs(an(1,1,izdom))
             else
                write(99,*) ipoin,0.0_rp
             end if
          end do
       end do
    end if
    write(99,*) 'End Values'

!!$    if( present(lgrou) ) then
!!$       !
!!$       ! Groups
!!$       !
!!$       write(99,*) 'Result GROUPS ALYA  0.00000000E+00 Scalar OnGaussPoints GP'
!!$       write(99,*) 'ComponentNames GROUPS'
!!$       write(99,*) 'Values'
!!$
!!$       if( .not. present(invpr) ) then
!!$          do ii = 1,npoin
!!$             do jj = 1,npoin
!!$                ipoin = (ii-1)*ex + jj
!!$                if( notfound ) then
!!$                   write(99,*) ipoin,0.0_rp
!!$                else
!!$                   write(99,*) ipoin,lgrou(ipoin)
!!$                end if
!!$             end do
!!$          end do
!!$       else
!!$          do ii = 1,npoin
!!$             ii_old = invpr(ii)
!!$             do jj = 1,npoin
!!$                jj_old = invpr(jj)
!!$                ipoin = (ii-1)*ex + jj
!!$                if( lgrou(ii_old) == lgrou(jj_old) ) then
!!$                   write(99,*) ipoin,lgrou(ii_old)
!!$                else
!!$                   write(99,*) ipoin,0.0_rp
!!$                end if
!!$             end do
!!$          end do
!!$       end if
!!$       write(99,*) 'End Values'
!!$    end if

    close( unit = 98 )
    close( unit = 99 )

  end subroutine matrix_output_gid_format


  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    04/02/2016
  !> @brief   Assemble element matrix
  !> @details Assemble an element matrix ELMAT into the global matrix AN
  !>
  !----------------------------------------------------------------------

#ifndef _OPENMP  
  pure subroutine matrix_assemble_element_matrix_to_ELL(&
       ndof,pnode,pevat,lnods,elmat,ia,an)
#else
  subroutine matrix_assemble_element_matrix_to_ELL(&
       ndof,pnode,pevat,lnods,elmat,ia,an)
#endif  
    integer(ip), intent(in)             :: ndof
    integer(ip), intent(in)             :: pnode
    integer(ip), intent(in)             :: pevat
    integer(ip), intent(in)             :: lnods(pnode)
    integer(ip), intent(in),   pointer  :: ia(:,:)
    real(rp),    intent(in)             :: elmat(pevat,pevat)
    real(rp),    intent(inout)          :: an(ndof,ndof,*)
    integer(ip)                         :: ii,jj,iz,ncols
    integer(ip)                         :: inode,jnode,jz
    integer(ip)                         :: idof,jdof
    integer(ip)                         :: idol,jdol

    ncols = size(ia,1,KIND=ip)

    if( ndof == 1 ) then
       do inode = 1,pnode
          ii = lnods(inode)
          do jnode  = 1,pnode
             jj     = lnods(jnode)
             iz     = 1
             do while( ia(iz,ii) /= jj )
                iz = iz + 1
             end do
             jz         = (ii-1)*ncols + iz
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             an(1,1,jz) = an(1,1,jz) + elmat(inode,jnode)
          end do
       end do
    else
       do inode = 1,pnode
          ii = lnods(inode)
          do jnode  = 1,pnode
             jj     = lnods(jnode)
             iz     = 1
             do while( ia(iz,ii) /= jj )
                iz = iz + 1
             end do
             jz = (ii-1)*ncols + iz
             do idof = 1,ndof
                idol = (idof-1)*ndof+inode
                do jdof = 1,ndof
                   jdol = (jdof-1)*ndof+jnode
#ifdef NO_COLORING
                   !$OMP ATOMIC
#endif
                   an(jdof,idof,jz) = an(jdof,idof,jz) + elmat(idol,jdol)
                end do
             end do
          end do
       end do
    end if

  end subroutine matrix_assemble_element_matrix_to_ELL

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    04/02/2016
  !> @brief   Assemble element matrix
  !> @details Assemble an element matrix ELMAT into the global matrices
  !>          in CSR format. Additional block matrices can also be asembled
  !>          at the same time
  !>
  !----------------------------------------------------------------------
  
  subroutine matrix_assemble_2by2_block_element_matrix_to_CSR(&
       kfl_element_to_csr,nbva1,nbva2,pnode,pevat,&
       ielem,lnods,elmat,ia,ja,a11,a12,a21,a22,&
       element_to_csr,elb11,elb22,b11,b22)

    integer(ip), intent(in)                      :: kfl_element_to_csr
    integer(ip), intent(in)                      :: nbva1
    integer(ip), intent(in)                      :: nbva2
    integer(ip), intent(in)                      :: pnode
    integer(ip), intent(in)                      :: pevat
    integer(ip), intent(in)                      :: ielem
    integer(ip), intent(in)                      :: lnods(pnode)
    integer(ip), intent(in)                      :: ia(*)
    integer(ip), intent(in)                      :: ja(*)
    real(rp),    intent(in)                      :: elmat(pevat,pevat)
    real(rp),    intent(inout)                   :: a11(nbva1,nbva1,*)
    real(rp),    intent(inout)                   :: a12(nbva2,nbva1,*)
    real(rp),    intent(inout)                   :: a21(nbva1,nbva2,*)
    real(rp),    intent(inout)                   :: a22(nbva2,nbva2,*)
    integer(ip), intent(in),   pointer, optional :: element_to_csr(:,:,:)
    
    real(rp),    intent(in),            optional :: elb11(pnode*nbva1,pnode*nbva1)
    real(rp),    intent(in),            optional :: elb22(pnode*nbva2,pnode*nbva2)
    real(rp),    intent(inout),         optional :: b11(nbva1,nbva1,*)
    real(rp),    intent(inout),         optional :: b22(nbva2,nbva2,*)
    integer(ip)                                  :: inode,jnode
    integer(ip)                                  :: ipoin,jpoin,izsol,jcolu
    logical(lg)                                  :: use_element_to_csr

    if( kfl_element_to_csr == 1.and. present(element_to_csr) ) then
       use_element_to_csr = .true.
    else
       use_element_to_csr = .false.
    end if
    if( pevat/pnode /= nbva1+nbva2 ) call runend('MOD_MATRIX: WRONG MATRIX FORMAT FOR 2 BY 2 ASSEMBLY')

    if( nbva1 == 1 .and. nbva2 == 1 ) then
       !
       ! 1 unknown in each block
       !
       if( use_element_to_csr ) then
          do inode = 1,pnode
             do jnode = 1,pnode
                izsol = element_to_csr(inode,jnode,ielem)
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
                a11(1,1,izsol) = a11(1,1,izsol) + elmat((inode-1)*2+1,(jnode-1)*2+1)
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
                a12(1,1,izsol) = a12(1,1,izsol) + elmat((inode-1)*2+1,(jnode-1)*2+2)
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
                a21(1,1,izsol) = a21(1,1,izsol) + elmat((inode-1)*2+2,(jnode-1)*2+1)
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
                a22(1,1,izsol) = a22(1,1,izsol) + elmat((inode-1)*2+2,(jnode-1)*2+2)
                if( present(elb11) .and. present(b11) ) then
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif                   
                   b11(1,1,izsol) = b11(1,1,izsol) + elb11(inode,jnode)
                end if
                if( present(elb22) .and. present(b22) ) then
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif                   
                   b22(1,1,izsol) = b22(1,1,izsol) + elb22(inode,jnode)
                end if
             end do
          end do
       else
          do inode = 1,pnode
             ipoin = lnods(inode)
             do jnode = 1,pnode
                jpoin = lnods(jnode)
                izsol = ia(ipoin)
                jcolu = ja(izsol)
                do while( jcolu /= jpoin .and. izsol < ia(ipoin+1)-1)
                   izsol = izsol + 1
                   jcolu = ja(izsol)
                end do
                if( jcolu == jpoin ) then
#ifdef NO_COLORING
                   !$OMP ATOMIC
#endif
                   a11(1,1,izsol) = a11(1,1,izsol) + elmat((inode-1)*2+1,(jnode-1)*2+1)
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
                   a12(1,1,izsol) = a12(1,1,izsol) + elmat((inode-1)*2+1,(jnode-1)*2+2)
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
                   a21(1,1,izsol) = a21(1,1,izsol) + elmat((inode-1)*2+2,(jnode-1)*2+1)
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
                   a22(1,1,izsol) = a22(1,1,izsol) + elmat((inode-1)*2+2,(jnode-1)*2+2)

                   if( present(elb11) .and. present(b11) ) then
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
                      b11(1,1,izsol) = b11(1,1,izsol) + elb11(inode,jnode)
                   end if
                   if( present(elb22) .and. present(b22) ) then
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
                      b22(1,1,izsol) = b22(1,1,izsol) + elb22(inode,jnode)
                   end if
                end if
             end do
          end do
       end if

    else
       
       call runend('MOD_MATRIX: 2 BY 2 ASSEMBLY NOT CODED')

    end if

  end subroutine matrix_assemble_2by2_block_element_matrix_to_CSR
     
#ifndef _OPENMP  
  pure subroutine matrix_assemble_element_matrix_to_CSR(&
       kfl_element_to_csr,nbvar,pnode,pevat,&
       ielem,lnods,elmat,ia,ja,an,element_to_csr)
#else
  subroutine matrix_assemble_element_matrix_to_CSR(&
       kfl_element_to_csr,nbvar,pnode,pevat,&
       ielem,lnods,elmat,ia,ja,an,element_to_csr)    
#endif  
    integer(ip), intent(in)                      :: kfl_element_to_csr
    integer(ip), intent(in)                      :: nbvar
    integer(ip), intent(in)                      :: pnode
    integer(ip), intent(in)                      :: pevat
    integer(ip), intent(in)                      :: ielem
    integer(ip), intent(in)                      :: lnods(pnode)
    integer(ip), intent(in)                      :: ia(*)
    integer(ip), intent(in)                      :: ja(*)
    real(rp),    intent(in)                      :: elmat(pevat,pevat)
    real(rp),    intent(inout)                   :: an(nbvar,nbvar,*)
    integer(ip), intent(in),   pointer, optional :: element_to_csr(:,:,:)
    integer(ip)                                  :: ievat,jevat,idofn,jdofn
    integer(ip)                                  :: inode,jnode
    integer(ip)                                  :: ipoin,jpoin,izsol,jcolu
    logical(lg)                                  :: use_element_to_csr

    if( kfl_element_to_csr == 1 .and. present(element_to_csr) ) then
       use_element_to_csr = .true.
    else
       use_element_to_csr = .false.
    end if

    if( nbvar == 1 ) then
       !
       ! 1 unknown
       !
       if( use_element_to_csr ) then
          do inode = 1,pnode
             do jnode = 1,pnode
                izsol = element_to_csr(inode,jnode,ielem)
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
                an(1,1,izsol) = an(1,1,izsol) + elmat(inode,jnode)
             end do
          end do
       else
          do inode = 1,pnode
             ipoin = lnods(inode)
             do jnode = 1,pnode
                jpoin = lnods(jnode)
                izsol = ia(ipoin)
                jcolu = ja(izsol)
                do while( jcolu /= jpoin .and. izsol < ia(ipoin+1)-1)
                   izsol = izsol + 1
                   jcolu = ja(izsol)
                end do
                if( jcolu == jpoin ) then
#ifdef NO_COLORING
                   !$OMP ATOMIC
#endif
                   an(1,1,izsol) = an(1,1,izsol) + elmat(inode,jnode)
                end if
             end do
          end do
       end if

    else if( nbvar == 2 ) then
       !
       ! General case: 2 unknowns
       !
       if( use_element_to_csr ) then
          do inode = 1,pnode
             do jnode = 1,pnode
                izsol = element_to_csr(inode,jnode,ielem)
                do idofn = 1,2
                   ievat = (inode-1) * 2 + idofn
                   do jdofn = 1,2
                      jevat = (jnode-1) * 2 + jdofn
#ifdef NO_COLORING
                      !$OMP ATOMIC
#endif
                      an(jdofn,idofn,izsol) = an(jdofn,idofn,izsol) + elmat(ievat,jevat)
                   end do
                end do
             end do
          end do
       else
          do inode = 1,pnode
             ipoin = lnods(inode)
             do jnode = 1,pnode
                jpoin = lnods(jnode)
                izsol = ia(ipoin)
                jcolu = ja(izsol)
                do while( jcolu /= jpoin .and. izsol < ia(ipoin+1)-1 )
                   izsol = izsol + 1
                   jcolu = ja(izsol)
                end do
                if( jpoin == jcolu ) then
                   do idofn = 1,2
                      ievat = (inode-1) * 2 + idofn
                      do jdofn = 1,2
                         jevat = (jnode-1) * 2 + jdofn
#ifdef NO_COLORING
                         !$OMP ATOMIC
#endif
                         an(jdofn,idofn,izsol) = an(jdofn,idofn,izsol) + elmat(ievat,jevat)
                      end do
                   end do
                end if
             end do
          end do
       end if

    else if( nbvar == 3 ) then
       !
       ! General case: 3 unknowns
       !
       if( use_element_to_csr ) then
          do inode = 1,pnode
             do jnode = 1,pnode
                izsol = element_to_csr(inode,jnode,ielem)
                do idofn = 1,3
                   ievat = (inode-1) * 3 + idofn
                   do jdofn = 1,3
                      jevat = (jnode-1) * 3 + jdofn
#ifdef NO_COLORING
                      !$OMP ATOMIC
#endif
                      an(jdofn,idofn,izsol) = an(jdofn,idofn,izsol) + elmat(ievat,jevat)
                   end do
                end do
             end do
          end do
       else
          do inode = 1,pnode
             ipoin = lnods(inode)
             do jnode = 1,pnode
                jpoin = lnods(jnode)
                izsol = ia(ipoin)
                jcolu = ja(izsol)
                do while( jcolu /= jpoin .and. izsol < ia(ipoin+1)-1 )
                   izsol = izsol + 1
                   jcolu = ja(izsol)
                end do
                if( jpoin == jcolu ) then
                   do idofn = 1,3
                      ievat = (inode-1) * 3 + idofn
                      do jdofn = 1,3
                         jevat = (jnode-1) * 3 + jdofn
#ifdef NO_COLORING
                         !$OMP ATOMIC
#endif
                         an(jdofn,idofn,izsol) = an(jdofn,idofn,izsol) + elmat(ievat,jevat)
                      end do
                   end do
                end if
             end do
          end do
       end if

    else
       !
       ! General case: NBVAR unknowns
       !
       if( use_element_to_csr ) then
          do inode = 1,pnode
             do jnode = 1,pnode
                izsol = element_to_csr(inode,jnode,ielem)
                do idofn = 1,nbvar
                   ievat = (inode-1) * nbvar + idofn
                   do jdofn = 1,nbvar
                      jevat = (jnode-1) * nbvar + jdofn
#ifdef NO_COLORING
                      !$OMP ATOMIC
#endif
                      an(jdofn,idofn,izsol) = an(jdofn,idofn,izsol) + elmat(ievat,jevat)
                   end do
                end do
             end do
          end do
       else
          do inode = 1,pnode
             ipoin = lnods(inode)
             do jnode = 1,pnode
                jpoin = lnods(jnode)
                izsol = ia(ipoin)
                jcolu = ja(izsol)
                do while( jcolu /= jpoin .and. izsol < ia(ipoin+1)-1 )
                   izsol = izsol + 1
                   jcolu = ja(izsol)
                end do
                if( jpoin == jcolu ) then
                   do idofn = 1,nbvar
                      ievat = (inode-1) * nbvar + idofn
                      do jdofn = 1,nbvar
                         jevat = (jnode-1) * nbvar + jdofn
#ifdef NO_COLORING
                         !$OMP ATOMIC
#endif
                         an(jdofn,idofn,izsol) = an(jdofn,idofn,izsol) + elmat(ievat,jevat)
                      end do
                   end do
                end if
             end do
          end do
       end if

    end if

  end subroutine matrix_assemble_element_matrix_to_CSR

#ifndef _OPENMP  
  pure subroutine matrix_assemble_element_matrix_to_CSR_explicit(&
       nbvar_glo,nbvar_ele,pnode,lnods,elrhs,rhsid,elmat,elunk)
#else
  subroutine matrix_assemble_element_matrix_to_CSR_explicit(&
       nbvar_glo,nbvar_ele,pnode,lnods,elrhs,rhsid,elmat,elunk)    
#endif  
    integer(ip),  intent(in)           :: nbvar_glo
    integer(ip),  intent(in)           :: nbvar_ele
    integer(ip),  intent(in)           :: pnode
    integer(ip),  intent(in)           :: lnods(pnode)
    real(rp),     intent(in)           :: elrhs(nbvar_ele,pnode)
    real(rp),     intent(inout)        :: rhsid(nbvar_glo,*)
    real(rp),    intent(in)             :: elmat(pnode,pnode)
    real(rp),     intent(in)           :: elunk(pnode)
    integer(ip)                         :: inode,jnode
    integer(ip)                         :: ipoin

    if( nbvar_glo == 1 ) then
       !
       ! 1 DOF
       !
       do inode = 1,pnode
          ipoin = lnods(inode)
          do jnode = 1,pnode
             rhsid(1,ipoin) = rhsid(1,ipoin) - elmat(inode,jnode)*elunk(jnode)
          end do
       end do
    end if

  end subroutine matrix_assemble_element_matrix_to_CSR_explicit

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Assemble a RHS
  !> @details Assemble a RHS
  !
  !----------------------------------------------------------------------

#ifndef _OPENMP  
  pure subroutine matrix_assemble_element_RHS(&
       nbvar_glo,nbvar_ele,pnode,lnods,elrhs,rhsid,kdofn)
#else
  subroutine matrix_assemble_element_RHS(&
       nbvar_glo,nbvar_ele,pnode,lnods,elrhs,rhsid,kdofn)    
#endif  
    integer(ip),  intent(in)           :: nbvar_glo
    integer(ip),  intent(in)           :: nbvar_ele
    integer(ip),  intent(in)           :: pnode
    integer(ip),  intent(in)           :: lnods(pnode)
    real(rp),     intent(in)           :: elrhs(nbvar_ele,pnode)
    real(rp),     intent(inout)        :: rhsid(nbvar_glo,*)
    integer(ip),  intent(in), optional :: kdofn
    integer(ip)                        :: inode,ipoin,idofn

    if( present(kdofn) ) then
       !
       ! NBVAR_GLO DOF's but only KDOFN is assembled
       !
       if( nbvar_ele == 1 ) then
          do inode = 1,pnode
             ipoin = lnods(inode)
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             rhsid(kdofn,ipoin) = rhsid(kdofn,ipoin) + elrhs(1,inode)
          end do
       else
          !call runend('MATRIX_ASSRHS: VERY STRANGE RHS ASSEMBLY INDEED...')
       end if

    else

       if( nbvar_glo == 1 ) then
          !
          ! 1 DOF
          !
          do inode = 1,pnode
             ipoin = lnods(inode)
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             rhsid(1,ipoin) = rhsid(1,ipoin) + elrhs(1,inode)
          end do

       else if( nbvar_glo == 2 ) then
          !
          ! 2 DOF's
          !
          do inode = 1,pnode
             ipoin = lnods(inode)
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             rhsid(1,ipoin) = rhsid(1,ipoin) + elrhs(1,inode)
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             rhsid(2,ipoin) = rhsid(2,ipoin) + elrhs(2,inode)
          end do

       else if( nbvar_glo == 3 ) then
          !
          ! 3 DOF's
          !
          do inode = 1,pnode
             ipoin = lnods(inode)
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             rhsid(1,ipoin) = rhsid(1,ipoin) + elrhs(1,inode)
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             rhsid(2,ipoin) = rhsid(2,ipoin) + elrhs(2,inode)
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             rhsid(3,ipoin) = rhsid(3,ipoin) + elrhs(3,inode)
          end do

       else if( nbvar_glo > 2 ) then
          !
          ! >2 DOF's
          !
          do inode = 1,pnode
             ipoin = lnods(inode)
             do idofn = 1,nbvar_glo
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
                rhsid(idofn,ipoin) = rhsid(idofn,ipoin) + elrhs(idofn,inode)
             end do
          end do

       else if(nbvar_glo < 0 ) then
          !
          ! >2 DOF's
          !
          !call runend('MATRIX_ASSEMBLE_ELEMENT_RHS: NOT CODED')

       end if
    end if

  end subroutine matrix_assemble_element_RHS

  subroutine matrix_assemble_2by2_block_element_RHS(&
       nbvar_glo,nbvar_ele,pnode,lnods,elrhs,rhsid,nn,kdofn)
    
    integer(ip),  intent(in)           :: nbvar_glo
    integer(ip),  intent(in)           :: nbvar_ele
    integer(ip),  intent(in)           :: pnode
    integer(ip),  intent(in)           :: lnods(pnode)
    real(rp),     intent(in)           :: elrhs(nbvar_ele,pnode)
    real(rp),     intent(inout)        :: rhsid(*)
    integer(ip),  intent(in)           :: nn
    integer(ip),  intent(in), optional :: kdofn
    integer(ip)                        :: inode,ipoin,idofn

    if( present(kdofn) ) then
       !
       ! NBVAR_GLO DOF's but only KDOFN is assembled
       !
       if( nbvar_ele == 1 ) then
          do inode = 1,pnode
             ipoin = lnods(inode)
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             rhsid(nn*(kdofn-1)+ipoin) = rhsid(nn*(kdofn-1)+ipoin) + elrhs(1,inode)
          end do
       else
          call runend('MATRIX_ASSRHS: VERY STRANGE RHS ASSEMBLY INDEED...')
       end if

    else

       if( nbvar_glo == 1 ) then
          !
          ! 1 DOF
          !
          do inode = 1,pnode
             ipoin = lnods(inode)
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             rhsid(ipoin) = rhsid(ipoin) + elrhs(1,inode)
          end do

       else if( nbvar_glo == 2 ) then
          !
          ! 2 DOF's
          !
          do inode = 1,pnode
             ipoin = lnods(inode)
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             rhsid(ipoin) = rhsid(ipoin) + elrhs(1,inode)
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             rhsid(nn+ipoin) = rhsid(nn+ipoin) + elrhs(2,inode)
             
          end do

       else if( nbvar_glo == 3 ) then
          !
          ! 3 DOF's
          !
          do inode = 1,pnode
             ipoin = lnods(inode)
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             rhsid(ipoin) = rhsid(ipoin) + elrhs(1,inode)
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             rhsid(nn+ipoin) = rhsid(nn+ipoin) + elrhs(2,inode)
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             rhsid(2*nn+ipoin) = rhsid(2*nn+ipoin) + elrhs(3,inode)
          end do

       else if( nbvar_glo > 2 ) then
          !
          ! >2 DOF's
          !
          do inode = 1,pnode
             ipoin = lnods(inode)
             do idofn = 1,nbvar_glo
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
                rhsid((idofn-1)*nn+ipoin) = rhsid((idofn-1)*nn+ipoin) + elrhs(idofn,inode)
             end do
          end do

       else if(nbvar_glo < 0 ) then
          !
          ! >2 DOF's
          !
          call runend('MATRIX_ASSEMBLE_ELEMENT_RHS: NOT CODED')

       end if
    end if

  end subroutine matrix_assemble_2by2_block_element_RHS

  subroutine matrix_copy_matrix(nbrows,ndof,iA,jA,Ain,Aout,iAout,jAout,invpr)

    !----------------------------------------------------------------------
    !>
    !> @author  Guillaume Houzeaux
    !> @date    19/07/2017
    !> @brief   Copy a matrix
    !> @details Copy a matrix Ain to a matrix Aout. If the graph of the
    !>          output matrix is present, then look for corresponding
    !>          positions of Ain
    !>
    !----------------------------------------------------------------------

    integer(ip), intent(in)                    :: nbrows            !< Number of rows
    integer(ip), intent(in)                    :: ndof              !< Number of dof per row
    integer(ip), intent(in), pointer           :: iA(:)             !< Matrix graph
    integer(ip), intent(in), pointer           :: jA(:)             !< Matrix graph
    real(rp)                                   :: Ain(ndof,ndof,*)  !< Input matrix
    real(rp),    intent(out)                   :: Aout(ndof,ndof,*) !< Output matrix
    integer(ip), intent(in), pointer, optional :: iAout(:)          !< Output matrix graph
    integer(ip), intent(in), pointer, optional :: jAout(:)          !< Output matrix graph
    integer(ip), intent(in), pointer, optional :: invpr(:)          !< Permutation
    integer(ip)                                :: ii,jj,iz,jz
    integer(ip)                                :: ii_loc,jj_loc
    !
    ! Allocate if necessary
    !
    !if( .not. associated(Aout) ) then
    !   iz = iA(nbrows+1)-1
    !   call memory_alloca(memor,'Aout','matrix_copy_matrix',Aout)
    !end if
    if( present(invpr) ) then
       if( present(iAout) .and. present(jAout) ) then
          do ii = 1,nbrows
             ii_loc = invpr(ii)
             if( ii_loc > 0 ) then
                do iz = iA(ii),iA(ii+1)-1
                   jj = jA(iz)
                   jj_loc = invpr(jj)
                   if( jj_loc > 0 ) then
                      call graphs_find_edge(ii_loc,jj_loc,iAout,jAout,jz)
                      if( jz /= 0 ) then
                         Aout(1:ndof,1:ndof,jz) = Ain(1:ndof,1:ndof,iz)                         
                      end if
                   end if
                end do
             end if
          end do
       end if
    else
       if( present(iAout) .and. present(jAout) ) then
          do ii = 1,nbrows
             do iz = iA(ii),iA(ii+1)-1
                jj = jA(iz)
                call graphs_find_edge(ii,jj,iAout,jAout,jz)
                if( jz /= 0 ) &
                     Aout(1:ndof,1:ndof,jz) = Ain(1:ndof,1:ndof,iz)
             end do
          end do
          
       else
          iz = iA(nbrows+1)-1
          Aout(1:ndof,1:ndof,1:iz) = Ain(1:ndof,1:ndof,1:iz)
       end if
    end if
  end subroutine matrix_copy_matrix

  subroutine matrix_copy_matrix_block(nbrows,ndof,iblock,iA,jA,Ain,Aout,iAout,jAout,invpr)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    matrix_copy_matrix_block
    ! DESCRIPTION
    !    This routine copy a block of a matrix to another
    ! INPUT ARGUMENTS
    !    NBROWS .... Dimension of the sysytem
    !    NDOF ...... Number of degrees of freedom in each node
    !    INVPR ..... Permutation array: NEW = INVPR(OLD)
    !     B ........ RHS of the system
    ! OUTPUT ARGUMENTS
    !     X ........ Solution vector of the system
    !---------------------------------------------------------------------------------------------------------------------------

    integer(ip), intent(in)             :: nbrows            !< Number of rows
    integer(ip), intent(in)             :: ndof              !< Number of dof per row
    integer(ip), intent(in)             :: iblock            !< Block to copy
    integer(ip), intent(in)             :: iA(*)             !< Matrix graph
    integer(ip), intent(in)             :: jA(*)             !< Matrix graph
    real(rp)                            :: Ain(ndof,ndof,*)  !< Input matrix
    real(rp),    intent(out)            :: Aout(*)           !< Output matrix
    integer(ip), intent(in), pointer, optional :: iAout(:)          !< Output matrix graph
    integer(ip), intent(in), pointer, optional :: jAout(:)          !< Output matrix graph
    integer(ip), intent(in),   pointer, optional :: invpr(:)          !< Permutation
    integer(ip)                         :: ii,jj,iz
    integer(ip)                         :: ii_loc,jj_loc,jz

    if( iblock > ndof ) call runend('MATRIX_COPY_MATRIX_BLOCK: WRING BLOCK NUMBER')
    
    if( present(invpr) .and. present(iAout) .and. present(jAout) ) then
       do ii = 1,nbrows
          ii_loc = invpr(ii)
          if( ii_loc > 0 ) then
             do iz = iA(ii),iA(ii+1)-1
                jj = jA(iz)
                jj_loc = invpr(jj)
                if( jj_loc > 0 ) then
                   call graphs_find_edge(ii_loc,jj_loc,iAout,jAout,jz)
                   if( jz /= 0 ) then
                      Aout(jz) = Ain(iblock,iblock,iz)
                   end if
                end if
             end do
          end if
       end do
    else
       do ii = 1,nbrows
          do iz = iA(ii),iA(ii+1)-1
             Aout(iz) = Ain(iblock,iblock,iz)
          end do
       end do
    end if
    
  end subroutine matrix_copy_matrix_block

  subroutine matrix_copy_matrix_to_block(nbrows,ndof,iA,jA,Ain,Aout)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    matrix_copy_matrix_block
    ! DESCRIPTION
    !    This routine copy a block of a matrix to another
    ! INPUT ARGUMENTS
    !    NBROWS .... Dimension of the sysytem
    !    NDOF ...... Number of degrees of freedom in each node
    !    INVPR ..... Permutation array: NEW = INVPR(OLD)
    !     B ........ RHS of the system
    ! OUTPUT ARGUMENTS
    !     X ........ Solution vector of the system
    !---------------------------------------------------------------------------------------------------------------------------

    integer(ip), intent(in)    :: nbrows            !< Number of rows
    integer(ip), intent(in)    :: ndof              !< Number of dof per row
    integer(ip), intent(in)    :: iA(*)             !< Matrix graph
    integer(ip), intent(in)    :: jA(*)             !< Matrix graph
    real(rp),    intent(in)    :: Ain(*)            !< Input matrix
    real(rp),    intent(inout) :: Aout(ndof,ndof,*) !< Output matrix
    integer(ip)                :: ii,iz
    integer(ip)                :: idof

    do ii = 1,nbrows
       do iz = iA(ii),iA(ii+1)-1
          Aout(1:ndof,1:ndof,iz) = 0.0_rp
          do idof = 1,ndof
             Aout(idof,idof,iz) = Ain(iz)
          end do
       end do
    end do

  end subroutine matrix_copy_matrix_to_block

  subroutine matrix_permute_and_copy_matrix(nbrows,ndof,iAin,jAin,iAout,jAout,invpR,Ain,Aout)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    matrix_permute_and_copy_matrix
    ! DESCRIPTION
    !    This routine copy a matrix to another
    ! INPUT ARGUMENTS
    !    NBROWS .... Dimension of the sysytem
    !    NDOF ...... Number of degrees of freedom in each node
    !    INVPR ..... Permutation array: NEW = INVPR(OLD)
    !     B ........ RHS of the system
    ! OUTPUT ARGUMENTS
    !     X ........ Solution vector of the system
    !---------------------------------------------------------------------------------------------------------------------------

    integer(ip), intent(in)             :: nbrows            !< Number of rows
    integer(ip), intent(in)             :: ndof              !< Number of dof per row
    integer(ip), intent(in)             :: iAin(*)           !< Matrix graph
    integer(ip), intent(in)             :: jAin(*)           !< Matrix graph
    integer(ip), intent(in)             :: iAout(*)          !< Matrix graph
    integer(ip), intent(in)             :: jAout(*)          !< Matrix graph
    integer(ip), pointer                :: invpR(:)          !< Inverse permutation
    real(rp)                            :: Ain(ndof,ndof,*)  !< Input matrix
    real(rp),    intent(out),  optional :: Aout(ndof,ndof,*) !< Output matrix
    integer(ip)                         :: ii,jj,iz
    integer(ip)                         :: jjold,iiold,kz
    real(rp),    pointer                :: Acpy(:,:,:)

    if( present(Aout) ) then

       if( .not. associated(invpR) ) then
          do ii = 1,nbrows
             do iz = iAin(ii),iAin(ii+1)-1
                Aout(1:ndof,1:ndof,iz) = Ain(1:ndof,1:ndof,iz)
             end do
          end do
       else
          do ii = 1,nbrows
             iiold = invpR(ii)
             do iz = iAout(ii),iAout(ii+1)-1
                jj    = jAout(iz)
                jjold = invpR(jj)
                kz    = iAin(iiold)
                do while( jAin(kz) /= jjold )
                   kz = kz + 1
                end do
                Aout(1:ndof,1:ndof,iz) = Ain(1:ndof,1:ndof,kz)
             end do
          end do
       end if

    else
       if( .not. associated(invpR) ) then
          return
       else
          nullify(Acpy)
          kz = iAin(nbrows+1)-1
          allocate(Acpy(ndof,ndof,kz))
          do ii = 1,nbrows
             do iz = iAin(ii),iAin(ii+1)-1
                Acpy(1:ndof,1:ndof,iz) = Ain(1:ndof,1:ndof,iz)
             end do
          end do
          do ii = 1,nbrows
             iiold = invpR(ii)
             do iz = iAout(ii),iAout(ii+1)-1
                jj    = jAout(iz)
                jjold = invpR(jj)
                kz    = iAin(iiold)
                do while( jAin(kz) /= jjold )
                   kz = kz + 1
                end do
                Ain(1:ndof,1:ndof,iz) = Acpy(1:ndof,1:ndof,kz)
             end do
          end do
          deallocate(Acpy)
       end if
    end if

  end subroutine matrix_permute_and_copy_matrix

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   SpMV
  !> @details Sparse matrix vector product in BCSR format
  !>          y = A x for rows n1:n2
  !>          Initial and end graphs can be different too
  !>          OpenMP can be activated through argument
  !>          YY can be initialized or not through argument
  !>
  !----------------------------------------------------------------------

  subroutine matrix_CSR_SpMV(n1,n2,ndofr,ndofc,ndofr_an,ndofc_an,ia,ja,an,xx,yy,OPENMP,INITIALIZATION,CHUNK,invpr,ia_opt,SCHEDULE)
    !
    ! Dummy arguments
    !
    integer(ip),  intent(in)                    :: n1                      !< Starting node
    integer(ip),  intent(in)                    :: n2                      !< Final node
    integer(ip),  intent(in)                    :: ndofr                   !< Number of dof per node (cols)
    integer(ip),  intent(in)                    :: ndofc                   !< Number of dof per node (rows)
    integer(ip),  intent(in)                    :: ndofr_an                !< Number of dof per node (cols)
    integer(ip),  intent(in)                    :: ndofc_an                !< Number of dof per node (rows)
    integer(ip),  intent(in), pointer           :: ia(:)                   !< Matrix graph ia
    integer(ip),  intent(in), pointer           :: ja(:)                   !< Matrix graph ja
    real(rp),     intent(in)                    :: an(ndofc_an,ndofr_an,*) !< Matrix
    real(rp),     intent(in)                    :: xx(ndofc,*)             !< Input vector
    real(rp),     intent(inout)                 :: yy(ndofr,*)             !< Output vector
    integer(ip),  intent(in), pointer, optional :: invpr(:)                !< Permutation
    logical(lg),  intent(in),          optional :: OPENMP                  !< If OpenMP should be used
    logical(lg),  intent(in),          optional :: INITIALIZATION          !< If array should be initialized
    integer(ip),  intent(in),          optional :: CHUNK                   !< Chunks for dynamic scheduling of OpenMP
    integer(ip),  intent(in), pointer, optional :: ia_opt(:)               !< Matrix graph end ia
    character(*), intent(in),          optional :: SCHEDULE                !< OpenMP schedule
    integer(ip)                                 :: ii,jj,iz,ki,kj
    integer(ip)                                 :: jjold,iiold
    integer(ip)                                 :: col,ll
    real(rp)                                    :: raux1,raux
    real(rp)                                    :: raux2,raux3
    logical(lg)                                 :: use_openmp
    logical(lg)                                 :: do_initialize
    integer(ip)                                 :: my_schedule
    integer(ip),              pointer           :: ia1(:)
    integer(ip),              pointer           :: ia2(:)
    integer(ip)                                 :: my_chunk
#if defined(ALYA_OMPSS) || defined(_OPENMP)
    integer(ip)                                 :: kk
#endif
    !
    ! If linked list should be different
    !
    if( present(ia_opt) ) then
       ia1 => ia
       ia2 => ia_opt
    else
       ia1 => ia
       ia2 => ia
    end if
    use_openmp    = .false.
    do_initialize = .true.

    if( present(OPENMP)         ) use_openmp = OPENMP
    if( present(INITIALIZATION) ) do_initialize = INITIALIZATION
    !
    ! OpenMP options
    !    
    if( use_openmp ) then

       my_chunk    = spmv_chunk
       my_schedule = SOL_OMP_STATIC

       if( present(CHUNK) ) my_chunk = CHUNK
       if( present(SCHEDULE) ) then
          if(      trim(SCHEDULE) == 'STATIC' ) then
             my_schedule = SOL_OMP_STATIC
          else if( trim(SCHEDULE) == 'DYNAMIC' ) then
             if( my_chunk <= 1 ) then
                my_schedule = SOL_OMP_STATIC
             else
                my_schedule = SOL_OMP_DYNAMIC
             end if
          else if( trim(SCHEDULE) == 'GUIDED' ) then
             my_schedule = SOL_OMP_GUIDED
          end if
       end if

    end if
    !
    ! Initialize solution
    !
    if( do_initialize ) then
       if( present(invpr) ) then
          do ii = n1,n2
             iiold = invpr(ii)
             yy(1:ndofr,iiold) = 0.0_rp
          end do
       else
          do ii = n1,n2
             yy(1:ndofr,ii) = 0.0_rp
          end do
       end if
    end if

    if( present(invpr) ) then
       !
       ! Input  x is in old numbering
       ! Output y is in old numbering (old=invpr(new))
       ! Matrix and graphs are in new numbering
       !
       do ii = n1,n2
          iiold = invpr(ii)
          do iz = ia1(ii),ia2(ii+1)-1
             jj = ja(iz)
             jjold = invpr(jj)
             do kj = 1,ndofc
                do ki = 1,ndofr
                   yy(ki,iiold) = yy(ki,iiold) + an(kj,ki,iz) * xx(kj,jjold)
                end do
             end do
          end do
       end do

    else

       
       if( ndofr_an /= ndofr .and. ndofr_an == 1 ) then
          !
          ! NDOF = whatever
          !
          if( use_openmp ) then
             if( my_schedule == SOL_OMP_STATIC ) then
                !$OMP PARALLEL   DO                                               &
                !$OMP SCHEDULE ( STATIC )                                         &
                !$OMP DEFAULT  ( NONE )                                           &
                !$OMP SHARED   ( an, ia1, ia2, ja, n1, n2, ndofr, xx, yy )        &
                !$OMP PRIVATE  ( col, ii, jj, kk, ll )
                do ii = n1,n2
                   do jj  = ia1(ii),ia2(ii+1)-1
                      col = ja(jj)
                      do ll = 1,ndofr
                         yy(ll,ii) = yy(ll,ii) + an(1,1,jj) * xx(ll,col)
                      end do
                   end do
                end do
                !$OMP END PARALLEL DO
             else if( my_schedule == SOL_OMP_GUIDED ) then
                !$OMP PARALLEL   DO                                               &
                !$OMP SCHEDULE ( GUIDED )                                         &
                !$OMP DEFAULT  ( NONE )                                           &
                !$OMP SHARED   ( an, ia1, ia2, ja, n1, n2, ndofr, xx, yy )        &
                !$OMP PRIVATE  ( col, ii, jj, kk, ll )
                do ii = n1,n2
                   do jj  = ia1(ii),ia2(ii+1)-1
                      col = ja(jj)
                      do ll = 1,ndofr
                         yy(ll,ii) = yy(ll,ii) + an(1,1,jj) * xx(ll,col)
                      end do
                   end do
                end do
                !$OMP END PARALLEL DO
             else if( my_schedule == SOL_OMP_DYNAMIC ) then
                !$OMP PARALLEL   DO                                               &
                !$OMP SCHEDULE ( DYNAMIC , my_chunk )                             &
                !$OMP DEFAULT  ( NONE )                                           &
                !$OMP SHARED   ( an, ia1, ia2, ja, n1, n2, ndofr, xx, yy )        &
                !$OMP PRIVATE  ( col, ii, jj, kk, ll )
                do ii = n1,n2
                   do jj  = ia1(ii),ia2(ii+1)-1
                      col = ja(jj)
                      do ll = 1,ndofr
                         yy(ll,ii) = yy(ll,ii) + an(1,1,jj) * xx(ll,col)
                      end do
                   end do
                end do
                !$OMP END PARALLEL DO
             end if
          else
             do ii = n1,n2
                do jj  = ia1(ii),ia2(ii+1)-1
                   col = ja(jj)
                   do ll = 1,ndofr
                      raux = xx(ll,col)
                      yy(ll,ii) = yy(ll,ii) + an(1,1,jj) * xx(ll,col)
                   end do
                end do
             end do
          end if

       else if( ndofr == 1 .and. ndofc == 1 ) then
          !
          ! NDOF=1
          !
          if( use_openmp ) then
             if( my_schedule == SOL_OMP_STATIC ) then
                !$OMP PARALLEL   DO                                 &
                !$OMP SCHEDULE ( STATIC )                           &
                !$OMP DEFAULT  ( NONE )                             &
                !$OMP SHARED   ( an, ia1, ia2, ja, n1, n2, xx, yy ) &
                !$OMP PRIVATE  ( col, ii, jj )
                do ii = n1,n2
                   do jj = ia1(ii),ia2(ii+1)-1
                      col      = ja(jj)
                      yy(1,ii) = yy(1,ii) + an(1,1,jj) * xx(1,col)
                   end do
                end do
                !$OMP END PARALLEL DO
             else if( my_schedule == SOL_OMP_GUIDED ) then
                !$OMP PARALLEL   DO                                 &
                !$OMP SCHEDULE ( GUIDED )                           &
                !$OMP DEFAULT  ( NONE )                             &
                !$OMP SHARED   ( an, ia1, ia2, ja, n1, n2, xx, yy ) &
                !$OMP PRIVATE  ( col, ii, jj )
                do ii = n1,n2
                   do jj = ia1(ii),ia2(ii+1)-1
                      col      = ja(jj)
                      yy(1,ii) = yy(1,ii) + an(1,1,jj) * xx(1,col)
                   end do
                end do
                !$OMP END PARALLEL DO
             else if( my_schedule == SOL_OMP_DYNAMIC ) then
                !$OMP PARALLEL   DO                                 &
                !$OMP SCHEDULE ( DYNAMIC, my_chunk )                &
                !$OMP DEFAULT  ( NONE )                             &
                !$OMP SHARED   ( an, ia1, ia2, ja, n1, n2, xx, yy ) &
                !$OMP PRIVATE  ( col, ii, jj )
                do ii = n1,n2
                   do jj = ia1(ii),ia2(ii+1)-1
                      col      = ja(jj)
                      yy(1,ii) = yy(1,ii) + an(1,1,jj) * xx(1,col)
                   end do
                end do
                !$OMP END PARALLEL DO
             end if
          else
             do ii = n1,n2
                do jj = ia1(ii),ia2(ii+1)-1
                   col      = ja(jj)
                   yy(1,ii) = yy(1,ii) + an(1,1,jj) * xx(1,col)
                end do
             end do
          end if

       else if( ndofr == 2 .and. ndofc == 2 ) then
          !
          ! NDOF=2
          !
          if( use_openmp ) then
             if( my_schedule == SOL_OMP_STATIC ) then
                !$OMP PARALLEL  DO                                    &
                !$OMP SCHEDULE ( STATIC )                             &
                !$OMP DEFAULT  ( NONE )                               &
                !$OMP SHARED   ( an, ia1, ia2, ja, n1, n2, xx, yy )   &
                !$OMP PRIVATE  ( col, ii, jj, raux1, raux2 )
                do ii = n1,n2
                   do jj       = ia1(ii),ia2(ii+1)-1
                      col      = ja(jj)
                      raux1    = xx(1,col)
                      raux2    = xx(2,col)
                      yy(1,ii) = yy(1,ii) + an(1,1,jj) * raux1 + an(2,1,jj) * raux2
                      yy(2,ii) = yy(2,ii) + an(1,2,jj) * raux1 + an(2,2,jj) * raux2
                   end do
                end do
                !$OMP END PARALLEL DO
             else if( my_schedule == SOL_OMP_GUIDED ) then
                !$OMP PARALLEL   DO                                   &
                !$OMP SCHEDULE ( GUIDED )                             &
                !$OMP DEFAULT  ( NONE )                               &
                !$OMP SHARED   ( an, ia1, ia2, ja, n1, n2, xx, yy )   &
                !$OMP PRIVATE  ( col, ii, jj, raux1, raux2 )
                do ii = n1,n2
                   do jj       = ia1(ii),ia2(ii+1)-1
                      col      = ja(jj)
                      raux1    = xx(1,col)
                      raux2    = xx(2,col)
                      yy(1,ii) = yy(1,ii) + an(1,1,jj) * raux1 + an(2,1,jj) * raux2
                      yy(2,ii) = yy(2,ii) + an(1,2,jj) * raux1 + an(2,2,jj) * raux2
                   end do
                end do
                !$OMP END PARALLEL DO
             else if( my_schedule == SOL_OMP_DYNAMIC ) then
                !$OMP PARALLEL   DO                                   &
                !$OMP SCHEDULE ( DYNAMIC , my_chunk )                 &
                !$OMP DEFAULT  ( NONE )                               &
                !$OMP SHARED   ( an, ia1, ia2, ja, n1, n2, xx, yy )   &
                !$OMP PRIVATE  ( col, ii, jj, raux1, raux2 )
                do ii = n1,n2
                   do jj       = ia1(ii),ia2(ii+1)-1
                      col      = ja(jj)
                      raux1    = xx(1,col)
                      raux2    = xx(2,col)
                      yy(1,ii) = yy(1,ii) + an(1,1,jj) * raux1 + an(2,1,jj) * raux2
                      yy(2,ii) = yy(2,ii) + an(1,2,jj) * raux1 + an(2,2,jj) * raux2
                   end do
                end do
                !$OMP END PARALLEL DO
             end if
          else
             do ii = n1,n2
                do jj       = ia1(ii),ia2(ii+1)-1
                   col      = ja(jj)
                   raux1    = xx(1,col)
                   raux2    = xx(2,col)
                   yy(1,ii) = yy(1,ii) + an(1,1,jj) * raux1 + an(2,1,jj) * raux2
                   yy(2,ii) = yy(2,ii) + an(1,2,jj) * raux1 + an(2,2,jj) * raux2
                end do
             end do
          end if

       else if( ndofr == 3 .and. ndofc == 3 ) then
          !
          ! NDOF=3
          !
          if( use_openmp ) then
             if( my_schedule == SOL_OMP_STATIC ) then
                !$OMP PARALLEL   DO                                    &
                !$OMP SCHEDULE ( STATIC )                              &
                !$OMP DEFAULT  ( NONE )                                &
                !$OMP SHARED   ( an, ia1, ia2, ja, n1, n2, xx, yy )    &
                !$OMP PRIVATE  ( col, ii, jj, raux1, raux2, raux3 )
                do ii = n1,n2
                   do jj       = ia1(ii),ia2(ii+1)-1
                      col      = ja(jj)
                      raux1    = xx(1,col)
                      raux2    = xx(2,col)
                      raux3    = xx(3,col)
                      yy(1,ii) = yy(1,ii) + an(1,1,jj) * raux1 + an(2,1,jj) * raux2 + an(3,1,jj) * raux3
                      yy(2,ii) = yy(2,ii) + an(1,2,jj) * raux1 + an(2,2,jj) * raux2 + an(3,2,jj) * raux3
                      yy(3,ii) = yy(3,ii) + an(1,3,jj) * raux1 + an(2,3,jj) * raux2 + an(3,3,jj) * raux3
                   end do
                end do
                !$OMP END PARALLEL DO
             else if( my_schedule == SOL_OMP_GUIDED ) then
                !$OMP PARALLEL   DO                                    &
                !$OMP SCHEDULE ( GUIDED )                              &
                !$OMP DEFAULT  ( NONE )                                &
                !$OMP SHARED   ( an, ia1, ia2, ja, n1, n2, xx, yy )    &
                !$OMP PRIVATE  ( col, ii, jj, raux1, raux2, raux3 )
                do ii = n1,n2
                   do jj       = ia1(ii),ia2(ii+1)-1
                      col      = ja(jj)
                      raux1    = xx(1,col)
                      raux2    = xx(2,col)
                      raux3    = xx(3,col)
                      yy(1,ii) = yy(1,ii) + an(1,1,jj) * raux1 + an(2,1,jj) * raux2 + an(3,1,jj) * raux3
                      yy(2,ii) = yy(2,ii) + an(1,2,jj) * raux1 + an(2,2,jj) * raux2 + an(3,2,jj) * raux3
                      yy(3,ii) = yy(3,ii) + an(1,3,jj) * raux1 + an(2,3,jj) * raux2 + an(3,3,jj) * raux3
                   end do
                end do
                !$OMP END PARALLEL DO
             else if( my_schedule == SOL_OMP_DYNAMIC ) then
                !$OMP PARALLEL   DO                                    &
                !$OMP SCHEDULE ( DYNAMIC , my_chunk )                  &
                !$OMP DEFAULT  ( NONE )                                &
                !$OMP SHARED   ( an, ia1, ia2, ja, n1, n2, xx, yy )    &
                !$OMP PRIVATE  ( col, ii, jj, raux1, raux2, raux3 )
                do ii = n1,n2
                   do jj       = ia1(ii),ia2(ii+1)-1
                      col      = ja(jj)
                      raux1    = xx(1,col)
                      raux2    = xx(2,col)
                      raux3    = xx(3,col)
                      yy(1,ii) = yy(1,ii) + an(1,1,jj) * raux1 + an(2,1,jj) * raux2 + an(3,1,jj) * raux3
                      yy(2,ii) = yy(2,ii) + an(1,2,jj) * raux1 + an(2,2,jj) * raux2 + an(3,2,jj) * raux3
                      yy(3,ii) = yy(3,ii) + an(1,3,jj) * raux1 + an(2,3,jj) * raux2 + an(3,3,jj) * raux3
                   end do
                end do
                !$OMP END PARALLEL DO
             end if
          else
             do ii = n1,n2
                do jj       = ia1(ii),ia2(ii+1)-1
                   col      = ja(jj)
                   raux1    = xx(1,col)
                   raux2    = xx(2,col)
                   raux3    = xx(3,col)
                   yy(1,ii) = yy(1,ii) + an(1,1,jj) * raux1 + an(2,1,jj) * raux2 + an(3,1,jj) * raux3
                   yy(2,ii) = yy(2,ii) + an(1,2,jj) * raux1 + an(2,2,jj) * raux2 + an(3,2,jj) * raux3
                   yy(3,ii) = yy(3,ii) + an(1,3,jj) * raux1 + an(2,3,jj) * raux2 + an(3,3,jj) * raux3
                end do
             end do
          end if

       else if( ndofr == 4 .and. ndofc == 4 ) then
          !
          ! NDOF=4
          !
          if( use_openmp ) then
             if( my_schedule == SOL_OMP_STATIC ) then
                !$OMP PARALLEL   DO                                   &
                !$OMP SCHEDULE ( STATIC )                             &
                !$OMP DEFAULT  ( NONE )                               &
                !$OMP SHARED   ( an, ia1, ia2, ja, n1, n2, xx, yy )   &
                !$OMP PRIVATE  ( col, ii, jj, raux )
                do ii = n1,n2
                   do jj       = ia1(ii),ia2(ii+1)-1
                      col      = ja(jj)
                      raux     = xx(1,col)
                      yy(1,ii) = yy(1,ii) + an(1,1,jj) * raux
                      yy(2,ii) = yy(2,ii) + an(1,2,jj) * raux
                      yy(3,ii) = yy(3,ii) + an(1,3,jj) * raux
                      yy(4,ii) = yy(4,ii) + an(1,4,jj) * raux
                      raux     = xx(2,col)
                      yy(1,ii) = yy(1,ii) + an(2,1,jj) * raux
                      yy(2,ii) = yy(2,ii) + an(2,2,jj) * raux
                      yy(3,ii) = yy(3,ii) + an(2,3,jj) * raux
                      yy(4,ii) = yy(4,ii) + an(2,4,jj) * raux
                      raux     = xx(3,col)
                      yy(1,ii) = yy(1,ii) + an(3,1,jj) * raux
                      yy(2,ii) = yy(2,ii) + an(3,2,jj) * raux
                      yy(3,ii) = yy(3,ii) + an(3,3,jj) * raux
                      yy(4,ii) = yy(4,ii) + an(3,4,jj) * raux
                      raux     = xx(4,col)
                      yy(1,ii) = yy(1,ii) + an(4,1,jj) * raux
                      yy(2,ii) = yy(2,ii) + an(4,2,jj) * raux
                      yy(3,ii) = yy(3,ii) + an(4,3,jj) * raux
                      yy(4,ii) = yy(4,ii) + an(4,4,jj) * raux
                   end do
                end do
                !$OMP END PARALLEL DO
             else if( my_schedule == SOL_OMP_GUIDED ) then
                !$OMP PARALLEL   DO                                   &
                !$OMP SCHEDULE ( GUIDED )                             &
                !$OMP DEFAULT  ( NONE )                               &
                !$OMP SHARED   ( an, ia1, ia2, ja, n1, n2, xx, yy )   &
                !$OMP PRIVATE  ( col, ii, jj, raux )
                do ii = n1,n2
                   do jj       = ia1(ii),ia2(ii+1)-1
                      col      = ja(jj)
                      raux     = xx(1,col)
                      yy(1,ii) = yy(1,ii) + an(1,1,jj) * raux
                      yy(2,ii) = yy(2,ii) + an(1,2,jj) * raux
                      yy(3,ii) = yy(3,ii) + an(1,3,jj) * raux
                      yy(4,ii) = yy(4,ii) + an(1,4,jj) * raux
                      raux     = xx(2,col)
                      yy(1,ii) = yy(1,ii) + an(2,1,jj) * raux
                      yy(2,ii) = yy(2,ii) + an(2,2,jj) * raux
                      yy(3,ii) = yy(3,ii) + an(2,3,jj) * raux
                      yy(4,ii) = yy(4,ii) + an(2,4,jj) * raux
                      raux     = xx(3,col)
                      yy(1,ii) = yy(1,ii) + an(3,1,jj) * raux
                      yy(2,ii) = yy(2,ii) + an(3,2,jj) * raux
                      yy(3,ii) = yy(3,ii) + an(3,3,jj) * raux
                      yy(4,ii) = yy(4,ii) + an(3,4,jj) * raux
                      raux     = xx(4,col)
                      yy(1,ii) = yy(1,ii) + an(4,1,jj) * raux
                      yy(2,ii) = yy(2,ii) + an(4,2,jj) * raux
                      yy(3,ii) = yy(3,ii) + an(4,3,jj) * raux
                      yy(4,ii) = yy(4,ii) + an(4,4,jj) * raux
                   end do
                end do
                !$OMP END PARALLEL DO
             else if( my_schedule == SOL_OMP_DYNAMIC ) then
                !$OMP PARALLEL   DO                                   &
                !$OMP SCHEDULE ( DYNAMIC , my_chunk )                 &
                !$OMP DEFAULT  ( NONE )                               &
                !$OMP SHARED   ( an, ia1, ia2, ja, n1, n2, xx, yy )   &
                !$OMP PRIVATE  ( col, ii, jj, raux )
                do ii = n1,n2
                   do jj       = ia1(ii),ia2(ii+1)-1
                      col      = ja(jj)
                      raux     = xx(1,col)
                      yy(1,ii) = yy(1,ii) + an(1,1,jj) * raux
                      yy(2,ii) = yy(2,ii) + an(1,2,jj) * raux
                      yy(3,ii) = yy(3,ii) + an(1,3,jj) * raux
                      yy(4,ii) = yy(4,ii) + an(1,4,jj) * raux
                      raux     = xx(2,col)
                      yy(1,ii) = yy(1,ii) + an(2,1,jj) * raux
                      yy(2,ii) = yy(2,ii) + an(2,2,jj) * raux
                      yy(3,ii) = yy(3,ii) + an(2,3,jj) * raux
                      yy(4,ii) = yy(4,ii) + an(2,4,jj) * raux
                      raux     = xx(3,col)
                      yy(1,ii) = yy(1,ii) + an(3,1,jj) * raux
                      yy(2,ii) = yy(2,ii) + an(3,2,jj) * raux
                      yy(3,ii) = yy(3,ii) + an(3,3,jj) * raux
                      yy(4,ii) = yy(4,ii) + an(3,4,jj) * raux
                      raux     = xx(4,col)
                      yy(1,ii) = yy(1,ii) + an(4,1,jj) * raux
                      yy(2,ii) = yy(2,ii) + an(4,2,jj) * raux
                      yy(3,ii) = yy(3,ii) + an(4,3,jj) * raux
                      yy(4,ii) = yy(4,ii) + an(4,4,jj) * raux
                   end do
                end do
                !$OMP END PARALLEL DO
             end if
          else
             do ii = n1,n2
                do jj       = ia1(ii),ia2(ii+1)-1
                   col      = ja(jj)
                   raux     = xx(1,col)
                   yy(1,ii) = yy(1,ii) + an(1,1,jj) * raux
                   yy(2,ii) = yy(2,ii) + an(1,2,jj) * raux
                   yy(3,ii) = yy(3,ii) + an(1,3,jj) * raux
                   yy(4,ii) = yy(4,ii) + an(1,4,jj) * raux
                   raux     = xx(2,col)
                   yy(1,ii) = yy(1,ii) + an(2,1,jj) * raux
                   yy(2,ii) = yy(2,ii) + an(2,2,jj) * raux
                   yy(3,ii) = yy(3,ii) + an(2,3,jj) * raux
                   yy(4,ii) = yy(4,ii) + an(2,4,jj) * raux
                   raux     = xx(3,col)
                   yy(1,ii) = yy(1,ii) + an(3,1,jj) * raux
                   yy(2,ii) = yy(2,ii) + an(3,2,jj) * raux
                   yy(3,ii) = yy(3,ii) + an(3,3,jj) * raux
                   yy(4,ii) = yy(4,ii) + an(3,4,jj) * raux
                   raux     = xx(4,col)
                   yy(1,ii) = yy(1,ii) + an(4,1,jj) * raux
                   yy(2,ii) = yy(2,ii) + an(4,2,jj) * raux
                   yy(3,ii) = yy(3,ii) + an(4,3,jj) * raux
                   yy(4,ii) = yy(4,ii) + an(4,4,jj) * raux
                end do
             end do
          end if

       else
          !
          ! NDOF = whatever
          !
          if( use_openmp ) then
             if( my_schedule == SOL_OMP_STATIC ) then
                !$OMP PARALLEL   DO                                               &
                !$OMP SCHEDULE ( STATIC )                                         &
                !$OMP DEFAULT  ( NONE )                                           &
                !$OMP SHARED   ( an, ia1, ia2, ja, n1, n2, ndofr, ndofc, xx, yy ) &
                !$OMP PRIVATE  ( col, ii, jj, kk, ll, raux )
                do ii = n1,n2
                   do jj  = ia1(ii),ia2(ii+1)-1
                      col = ja(jj)
                      do ll = 1,ndofc
                         raux = xx(ll,col)
                         yy(1:ndofr,ii) = yy(1:ndofr,ii) + an(ll,1:ndofr,jj) * raux
                      end do
                   end do
                end do
                !$OMP END PARALLEL DO
             else if( my_schedule == SOL_OMP_GUIDED ) then
                !$OMP PARALLEL   DO                                               &
                !$OMP SCHEDULE ( GUIDED )                                         &
                !$OMP DEFAULT  ( NONE )                                           &
                !$OMP SHARED   ( an, ia1, ia2, ja, n1, n2, ndofr, ndofc, xx, yy ) &
                !$OMP PRIVATE  ( col, ii, jj, kk, ll, raux )
                do ii = n1,n2
                   do jj  = ia1(ii),ia2(ii+1)-1
                      col = ja(jj)
                      do ll = 1,ndofc
                         raux = xx(ll,col)
                         yy(1:ndofr,ii) = yy(1:ndofr,ii) + an(ll,1:ndofr,jj) * raux
                      end do
                   end do
                end do
                !$OMP END PARALLEL DO
             else if( my_schedule == SOL_OMP_DYNAMIC ) then
                !$OMP PARALLEL   DO                                               &
                !$OMP SCHEDULE ( DYNAMIC , my_chunk )                             &
                !$OMP DEFAULT  ( NONE )                                           &
                !$OMP SHARED   ( an, ia1, ia2, ja, n1, n2, ndofr, ndofc, xx, yy ) &
                !$OMP PRIVATE  ( col, ii, jj, kk, ll, raux )
                do ii = n1,n2
                   do jj  = ia1(ii),ia2(ii+1)-1
                      col = ja(jj)
                      do ll = 1,ndofc
                         raux = xx(ll,col)
                         yy(1:ndofr,ii) = yy(1:ndofr,ii) + an(ll,1:ndofr,jj) * raux
                      end do
                   end do
                end do
                !$OMP END PARALLEL DO
             end if
          else
             do ii = n1,n2
                do jj  = ia1(ii),ia2(ii+1)-1
                   col = ja(jj)
                   do ll = 1,ndofc
                      raux = xx(ll,col)
                      yy(1:ndofr,ii) = yy(1:ndofr,ii) + an(ll,1:ndofr,jj) * raux
                   end do
                end do
             end do
          end if
       end if
    end if

  end subroutine matrix_CSR_SpMV

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    09/05/2017
  !> @brief   SpMV
  !> @details Copy matrices
  !
  !----------------------------------------------------------------------

  subroutine matrix_copy(&
       ndofn_rows,ndofn_cols,nzdom,aa_in,aa_out)

    integer(ip), intent(in)  :: ndofn_rows
    integer(ip), intent(in)  :: ndofn_cols
    integer(ip), intent(in)  :: nzdom
    real(rp),    intent(in)  :: aa_in(ndofn_cols,ndofn_rows,nzdom)
    real(rp),    intent(out) :: aa_out(ndofn_cols,ndofn_rows,nzdom)

    aa_out(1:ndofn_cols,1:ndofn_rows,1:nzdom) = aa_in(1:ndofn_cols,1:ndofn_rows,1:nzdom)

  end subroutine matrix_copy


  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    09/05/2017
  !> @brief   SpMV
  !> @details Parallel Sparse matrix vector product in BCSR format
  !>          y = A x using asynchronous communications
  !
  !----------------------------------------------------------------------

  subroutine matrix_CSR_parallel_asynchronous_SpMV(nrows,ncols,ndofn_rows,ndofn_cols,iA,jA,An,x,y,times,hybrid_opt)

    integer(ip), intent(in)            :: nrows                       !< Number of rows
    integer(ip), intent(in)            :: ncols                       !< Number of columns
    integer(ip), intent(in)            :: ndofn_rows                  !< Number of dof per node
    integer(ip), intent(in)            :: ndofn_cols                  !< Number of dof per node
    integer(ip), intent(in)            :: iA(*)                       !< Matrix graph ia
    integer(ip), intent(in)            :: jA(*)                       !< Matrix graph ja
    real(rp),    intent(in)            :: An(ndofn_cols,ndofn_rows,*) !< Matrix
    real(rp),    intent(in)            :: x(ndofn_cols,*)             !< Input vector
    real(rp),    intent(out)           :: y(ndofn_rows,*)             !< Output vector
    real(rp),    intent(out), optional :: times(3)                    !< Timings
    logical(lg), intent(in),  optional :: hybrid_opt                  !< Hybrid parallelism



  end subroutine matrix_CSR_parallel_asynchronous_SpMV

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    09/05/2017
  !> @brief   SpMV
  !> @details Parallel Sparse matrix vector product in BCSR format
  !>          y = A x using synchronous communications
  !
  !----------------------------------------------------------------------

  subroutine matrix_CSR_parallel_SpMV(nbnodes,ndofn_rows,ndofn_cols,iA,jA,An,x,y,times,hybrid_opt,what)

    integer(ip),  intent(in)            :: nbnodes                     !< Number of nodes
    integer(ip),  intent(in)            :: ndofn_rows                  !< Number of dof per node
    integer(ip),  intent(in)            :: ndofn_cols                  !< Number of dof per node
    integer(ip),  intent(in)            :: iA(*)                       !< Matrix graph ia
    integer(ip),  intent(in)            :: jA(*)                       !< Matrix graph ja
    real(rp),     intent(in)            :: An(ndofn_cols,ndofn_rows,*) !< Matrix
    real(rp),     intent(in)            :: x(ndofn_cols,*)             !< Input vector
    real(rp),     intent(out)           :: y(ndofn_rows,*)             !< Output vector
    real(rp),     intent(out), optional :: times(3)                    !< Timings
    character(*), intent(in),  optional :: hybrid_opt                  !< Kind of arallelism
    character(*), intent(in),  optional :: what                        !< What to do
    integer(ip)                         :: ii,jj
    integer(ip)                         :: kk
    integer(ip)                         :: ll,col
    real(rp)                            :: raux,raux1,raux2,raux3
    logical(lg)                         :: use_openmp,use_mpi

    if( INOTMASTER ) then

       use_mpi    = .true.
       use_openmp = .true.

       if( present(hybrid_opt) ) then
          if(      trim(hybrid_opt) == 'MPI/OPENMP' ) then
             use_mpi    = .true.
             use_openmp = .true.
          else if( trim(hybrid_opt) == 'MPI/NOT OPENMP' ) then
             use_mpi    = .true.
             use_openmp = .false.
          else if( trim(hybrid_opt) == 'NOT MPI/OPENMP' ) then
             use_mpi    = .false.
             use_openmp = .true.
          else if( trim(hybrid_opt) == 'NOT MPI/NOT OPENMP' ) then
             use_mpi    = .false.
             use_openmp = .false.
          else
             call runend('MATRIX_CSR_PARALLEL_SPMV: NON-EXISTING OPTION')
          end if
       end if

       if( present(times) ) call cputim(times(1))

       if( present(what) ) then

       end if

       if( ndofn_rows == 1 .and. ndofn_cols == 1 ) then
          !
          ! NDOF=1
          !
          if( use_openmp ) then
             !$OMP PARALLEL  DO                            &
             !$OMP SCHEDULE ( STATIC )                     &
             !$OMP DEFAULT  ( NONE )                       &
             !$OMP SHARED   ( an, ia, ja, nbnodes, x, y )  &
             !$OMP PRIVATE  ( col, ii, jj, raux )
             do ii = 1,nbnodes
                y(1,ii) = 0.0_rp
                do jj   = ia(ii),ia(ii+1)-1
                   col     = ja(jj)
                   raux    = x(1,col)
                   y(1,ii) = y(1,ii) + an(1,1,jj) * raux
                end do
             end do
             !$OMP END PARALLEL DO
          else
             do ii = 1,nbnodes
                y(1,ii) = 0.0_rp
                do jj   = ia(ii),ia(ii+1)-1
                   col     = ja(jj)
                   raux    = x(1,col)
                   y(1,ii) = y(1,ii) + an(1,1,jj) * raux
                end do
             end do
          end if

       else if( ndofn_rows == 2 .and. ndofn_cols == 2 ) then
          !
          ! NDOF=2
          !
          if( use_openmp ) then
             !$OMP PARALLEL  DO                            &
             !$OMP SCHEDULE ( STATIC )                     &
             !$OMP DEFAULT  ( NONE )                       &
             !$OMP SHARED   ( an, ia, ja, nbnodes, x, y )  &
             !$OMP PRIVATE  ( col, ii, jj, raux1, raux2 )
             do ii = 1,nbnodes
                y(1,ii) = 0.0_rp
                y(2,ii) = 0.0_rp
                do jj      = ia(ii),ia(ii+1)-1
                   col     = ja(jj)
                   raux1   = x(1,col)
                   raux2   = x(2,col)
                   y(1,ii) = y(1,ii) + an(1,1,jj) * raux1 + an(2,1,jj) * raux2
                   y(2,ii) = y(2,ii) + an(1,2,jj) * raux1 + an(2,2,jj) * raux2

                end do
             end do
             !$OMP END PARALLEL DO
          else
             do ii = 1,nbnodes
                y(1,ii) = 0.0_rp
                y(2,ii) = 0.0_rp
                do jj      = ia(ii),ia(ii+1)-1
                   col     = ja(jj)
                   raux1   = x(1,col)
                   raux2   = x(2,col)
                   y(1,ii) = y(1,ii) + an(1,1,jj) * raux1 + an(2,1,jj) * raux2
                   y(2,ii) = y(2,ii) + an(1,2,jj) * raux1 + an(2,2,jj) * raux2

                end do
             end do
          end if

       else if( ndofn_rows == 3 .and. ndofn_cols == 3 ) then
          !
          ! NDOF=3
          !
          if( use_openmp ) then
             !$OMP PARALLEL  DO                                   &
             !$OMP SCHEDULE ( STATIC )                            &
             !$OMP DEFAULT  ( NONE )                              &
             !$OMP SHARED   ( an, ia, ja, nbnodes, x, y )         &
             !$OMP PRIVATE  ( col, ii, jj, raux1, raux2, raux3 )
             do ii = 1,nbnodes
                y(1,ii) = 0.0_rp
                y(2,ii) = 0.0_rp
                y(3,ii) = 0.0_rp
                do jj      = ia(ii),ia(ii+1)-1
                   col     = ja(jj)
                   raux1   = x(1,col)
                   raux2   = x(2,col)
                   raux3   = x(3,col)
                   y(1,ii) = y(1,ii) + an(1,1,jj) * raux1 + an(2,1,jj) * raux2 + an(3,1,jj) * raux3
                   y(2,ii) = y(2,ii) + an(1,2,jj) * raux1 + an(2,2,jj) * raux2 + an(3,2,jj) * raux3
                   y(3,ii) = y(3,ii) + an(1,3,jj) * raux1 + an(2,3,jj) * raux2 + an(3,3,jj) * raux3
                end do
             end do
             !$OMP END PARALLEL DO
          else
             do ii = 1,nbnodes
                y(1,ii) = 0.0_rp
                y(2,ii) = 0.0_rp
                y(3,ii) = 0.0_rp
                do jj      = ia(ii),ia(ii+1)-1
                   col     = ja(jj)
                   raux1   = x(1,col)
                   raux2   = x(2,col)
                   raux3   = x(3,col)
                   y(1,ii) = y(1,ii) + an(1,1,jj) * raux1 + an(2,1,jj) * raux2 + an(3,1,jj) * raux3
                   y(2,ii) = y(2,ii) + an(1,2,jj) * raux1 + an(2,2,jj) * raux2 + an(3,2,jj) * raux3
                   y(3,ii) = y(3,ii) + an(1,3,jj) * raux1 + an(2,3,jj) * raux2 + an(3,3,jj) * raux3
                end do
             end do
          end if

       else if( ndofn_rows == 4 .and. ndofn_cols == 4 ) then
          !
          ! NDOF=4
          !
          if( use_openmp ) then
             !$OMP PARALLEL  DO                            &
             !$OMP SCHEDULE ( STATIC )                     &
             !$OMP DEFAULT  ( NONE )                       &
             !$OMP SHARED   ( an, ia, ja, nbnodes, x, y )  &
             !$OMP PRIVATE  ( col, ii, jj, raux )
             do ii = 1,nbnodes
                y(1,ii) = 0.0_rp
                y(2,ii) = 0.0_rp
                y(3,ii) = 0.0_rp
                y(4,ii) = 0.0_rp
                do jj = ia(ii),ia(ii+1)-1
                   col     = ja(jj)
                   raux    = x(1,col)
                   y(1,ii) = y(1,ii) + an(1,1,jj) * raux
                   y(2,ii) = y(2,ii) + an(1,2,jj) * raux
                   y(3,ii) = y(3,ii) + an(1,3,jj) * raux
                   y(4,ii) = y(4,ii) + an(1,4,jj) * raux
                   raux    = x(2,col)
                   y(1,ii) = y(1,ii) + an(2,1,jj) * raux
                   y(2,ii) = y(2,ii) + an(2,2,jj) * raux
                   y(3,ii) = y(3,ii) + an(2,3,jj) * raux
                   y(4,ii) = y(4,ii) + an(2,4,jj) * raux
                   raux    = x(3,col)
                   y(1,ii) = y(1,ii) + an(3,1,jj) * raux
                   y(2,ii) = y(2,ii) + an(3,2,jj) * raux
                   y(3,ii) = y(3,ii) + an(3,3,jj) * raux
                   y(4,ii) = y(4,ii) + an(3,4,jj) * raux
                   raux    = x(4,col)
                   y(1,ii) = y(1,ii) + an(4,1,jj) * raux
                   y(2,ii) = y(2,ii) + an(4,2,jj) * raux
                   y(3,ii) = y(3,ii) + an(4,3,jj) * raux
                   y(4,ii) = y(4,ii) + an(4,4,jj) * raux
                end do
             end do
             !$OMP END PARALLEL DO
          else
             do ii = 1,nbnodes
                y(1,ii) = 0.0_rp
                y(2,ii) = 0.0_rp
                y(3,ii) = 0.0_rp
                y(4,ii) = 0.0_rp
                do jj = ia(ii),ia(ii+1)-1
                   col     = ja(jj)
                   raux    = x(1,col)
                   y(1,ii) = y(1,ii) + an(1,1,jj) * raux
                   y(2,ii) = y(2,ii) + an(1,2,jj) * raux
                   y(3,ii) = y(3,ii) + an(1,3,jj) * raux
                   y(4,ii) = y(4,ii) + an(1,4,jj) * raux
                   raux    = x(2,col)
                   y(1,ii) = y(1,ii) + an(2,1,jj) * raux
                   y(2,ii) = y(2,ii) + an(2,2,jj) * raux
                   y(3,ii) = y(3,ii) + an(2,3,jj) * raux
                   y(4,ii) = y(4,ii) + an(2,4,jj) * raux
                   raux    = x(3,col)
                   y(1,ii) = y(1,ii) + an(3,1,jj) * raux
                   y(2,ii) = y(2,ii) + an(3,2,jj) * raux
                   y(3,ii) = y(3,ii) + an(3,3,jj) * raux
                   y(4,ii) = y(4,ii) + an(3,4,jj) * raux
                   raux    = x(4,col)
                   y(1,ii) = y(1,ii) + an(4,1,jj) * raux
                   y(2,ii) = y(2,ii) + an(4,2,jj) * raux
                   y(3,ii) = y(3,ii) + an(4,3,jj) * raux
                   y(4,ii) = y(4,ii) + an(4,4,jj) * raux
                end do
             end do

          end if

       else
          !
          ! NDOF = whatever
          !
          if( use_openmp ) then
             !$OMP PARALLEL  DO                                                    &
             !$OMP SCHEDULE ( STATIC )                                             &
             !$OMP DEFAULT  ( NONE )                                               &
             !$OMP SHARED   ( an, ia, ja, nbnodes, ndofn_rows, ndofn_cols, x, y )  &
             !$OMP PRIVATE  ( col, ii, jj, kk, ll, raux )
             do ii = 1,nbnodes
                do kk = 1,ndofn_rows
                   y(kk,ii) = 0.0_rp
                end do
                do jj  = ia(ii),ia(ii+1)-1
                   col = ja(jj)
                   do ll = 1,ndofn_cols
                      raux = x(ll,col)
                      do kk = 1,ndofn_rows
                         y(kk,ii) = y(kk,ii) + an(ll,kk,jj) * raux
                      end do
                   end do
                end do
             end do
             !$OMP END PARALLEL DO
          else
             do ii = 1,nbnodes
                do kk = 1,ndofn_rows
                   y(kk,ii) = 0.0_rp
                end do
                do jj  = ia(ii),ia(ii+1)-1
                   col = ja(jj)
                   do ll = 1,ndofn_cols
                      raux = x(ll,col)
                      do kk = 1,ndofn_rows
                         y(kk,ii) = y(kk,ii) + an(ll,kk,jj) * raux
                      end do
                   end do
                end do
             end do
          end if

       end if

       if( present(times) ) call cputim(times(2))
       !if( ISLAVE ) call rhsmod(ndofn,y)
       if( use_mpi .and. ISLAVE ) call matrix_parslx(ndofn_rows*nbnodes,y)
       if( present(times) ) call cputim(times(3))

    end if

  end subroutine matrix_CSR_parallel_SpMV

  !----------------------------------------------------------------------
  !
  !> @author  Ricard Borrell
  !> @date    15/07/2016
  !> @brief   SpMV
  !> @details Sparse matrix vector product in COO format
  !
  !----------------------------------------------------------------------

  subroutine matrix_COO_SpMV_BASE(nrows,nnz,ndofm,ndofv,iA,jA,vA,x,y,invpr,noini)

    implicit none
    integer(ip), intent(in)                    :: nrows               !< Number of rows
    integer(ip), intent(in)                    :: nnz                 !< Number of entries
    integer(ip), intent(in)                    :: ndofm               !< Number of dof per matrix entry
    integer(ip), intent(in)                    :: ndofv               !< Number of dof per vector entry
    integer(ip), intent(in)                    :: iA(nnz)             !< Matrix indices iA
    integer(ip), intent(in)                    :: jA(nnz)             !< Matrix indices jA
    real(rp),    intent(in)                    :: vA(ndofm,ndofm,nnz) !< Matrix values vA
    real(rp),    intent(in)                    :: x(ndofv,*)          !< Input vector
    real(rp),    intent(out)                   :: y(ndofv,*)          !< Output vector
    integer(ip), intent(in), optional, pointer :: invpr(:)            !< Permutation
    logical(lg), intent(in), optional          :: noini               !< No initialize
    integer(ip)                                :: ii,jj,ki,kj,kk
    !
    !initialize the solution vetor to 0
    !
    if((.not.present(noini)) .or. (noini .eqv. .false.))&
         y(1:ndofv,1:nrows)=0.0_rp
    if(ndofm /= ndofv .and. ndofm /= 1) call runend("Bad choice ndofm in COO_SpMV")
    if( present(invpr) ) then
       !
       ! Input  x is in old numbering
       ! Output y is in old numbering (old=invpr(new))
       ! Matrix and graphs are in new numbering
       !
       call runend(" matrix_COO_SPMV, not redy with invpr")
    else
       if(ndofm==ndofv) then
          do kk = 1,nnz
             ii=iA(kk)
             jj=jA(kk)
             do ki = 1,ndofv
                do kj = 1,ndofv
                   y(ki,ii) = y(ki,ii) + vA(kj,ki,kk) * x(kj,jj)
                end do
             end do
          end do
       else
          do kk = 1,nnz
             ii=iA(kk)
             jj=jA(kk)
             do ki = 1,ndofv
                y(ki,ii) = y(ki,ii) + vA(1,1,kk) * x(ki,jj)
             end do
          end do
       end if
    end if

  end subroutine matrix_COO_SpMV_BASE

  !----------------------------------------------------------------------
  !
  !> @author  Ricard Borrell
  !> @date    15/07/2016
  !> @brief   SpMV
  !> @details Sparse matrix vector product in COO format
  !
  !----------------------------------------------------------------------
  subroutine matrix_COO_SpMV_BLOCK(nrows,nnz,ndof,iA,jA,vA,x,y,invpr,noini)

    implicit none
    integer(ip), intent(in)                    :: nrows               !< Number of rows
    integer(ip), intent(in)                    :: nnz                 !< Number of entries
    integer(ip), intent(in)                    :: ndof                !< Number of dof
    integer(ip), intent(in)                    :: iA(nnz)             !< Matrix graph ia
    integer(ip), intent(in)                    :: jA(nnz)             !< Matrix graph ja
    real(rp),    intent(in)                    :: vA(ndof,ndof,nnz)   !< Matrix
    real(rp),    intent(in)                    :: x(ndof,*)           !< Input vector
    real(rp),    intent(out)                   :: y(ndof,*)           !< Output vector
    integer(ip), intent(in), optional, pointer :: invpr(:)            !< Permutation
    logical(lg), intent(in), optional          :: noini               !< No initialize
    logical(lg)                                :: auxlg

    auxlg = .false.
    if(present(noini)) auxlg = noini

    if(present(invpr)) then
       call matrix_COO_SpMV_BASE(nrows,nnz,ndof,ndof,iA,jA,vA,x,y,invpr,auxlg)
    else
       call matrix_COO_SpMV_BASE(nrows,nnz,ndof,ndof,iA,jA,vA,x,y,noini=auxlg)
    endif

  end subroutine matrix_COO_SpMV_BLOCK

  !----------------------------------------------------------------------
  !
  !> @author  Ricard Borrell
  !> @date    15/07/2016
  !> @brief   SpMV
  !> @details Sparse matrix vector product in COO format
  !
  !----------------------------------------------------------------------
  subroutine matrix_COO_SpMV_SPMAT(mat,ndof,x,y,invpr,noini)

    implicit none
    type(spmat), intent(in)                    :: mat                 !< Sparse matrix
    integer(ip), intent(in)                    :: ndof                !< Number of dof
    real(rp),    intent(in)                    :: x(ndof,*)           !< Input vector
    real(rp),    intent(out)                   :: y(ndof,*)           !< Output vector
    integer(ip), intent(in), optional, pointer :: invpr(:)            !< Permutation
    logical(lg), intent(in), optional          :: noini               !< No initialize
    logical(lg)                                :: auxlg

    auxlg = .false.
    if(present(noini)) auxlg = noini

    if(present(invpr)) then
       call matrix_COO_SpMV_BASE(mat % nrows,size(mat % iA,1,KIND=ip),mat % ndof1,ndof,&
            & mat % iA, mat % jA,mat % vA,x,y,invpr,auxlg)
    else
       call matrix_COO_SpMV_BASE(mat % nrows,size(mat % iA,1,KIND=ip),mat % ndof1,ndof,&
            & mat % iA, mat % jA,mat % vA,x,y,noini= auxlg)
    endif

  end subroutine matrix_COO_SpMV_SPMAT

  !----------------------------------------------------------------------
  !
  !> @author  Ricard Borrell
  !> @date    11/08/2016
  !> @brief   SPGEMM
  !> @details Sparse matrix matrix product in COO format
  !
  !----------------------------------------------------------------------
  subroutine   matrix_COO_SPGEMM(A,B,C,memit)

    implicit none
    type(spmat), intent(in)              :: A                 !< Sparse matrix
    type(spmat), intent(in)              :: B                 !< Sparse matrix
    type(spmat), intent(inout)           :: C                 !< Sparse matrix
    integer(8),  optional, intent(inout) :: memit(2)          !< Memory counter
    integer(8)                           :: memor(2)
    integer(ip), pointer                 :: entcolA(:)
    integer(ip), pointer                 :: entrowB(:)
    integer(ip), pointer                 :: entrowC(:)
    type(i1p), pointer                   :: adjA(:)
    type(i1p), pointer                   :: adjB(:)
    type(i1p), pointer                   :: adjC(:)
    type(r1p), pointer                   :: valA(:)
    type(r1p), pointer                   :: valB(:)
    type(r1p), pointer                   :: valC(:)
    integer(ip)                          :: ii,jj,kk,ll
    integer(ip)                          :: nentrC
    logical(lg)                          :: itsin
    !
    ! Check arguments
    !
    if(A % ndof1 /= 1 .or. B % ndof1 /= 1) &
         call runend("COO_SPGEMM only implemented for matrices with ndof==1")
    if(A % nrows <= 0 .or. A % ncols <= 0 .or. B % nrows <= 0 .or. B % ncols <= 0) &
         call runend("COO_SPGEMM: A or B matrix with nrows or ncols <= 0")
    if(A % ncols /= B % nrows) &
         call runend("COO_SPGEMM: A % ncols different than B % nrows")
    if(present(memit)) then
       memor = memit
    else
       memor = 0
    endif
    !
    ! Evaluate entries per row in B
    !
    nullify (entrowB)
    call memory_alloca(memor,'entrowB','matrix_COO_SPGEMM',entrowB, B % nrows)
    entrowB(:)=0
    do ii =1, size(B % iA,KIND=ip)
       entrowB(B % iA(ii)) = entrowB(B % iA(ii)) + 1
    end do
    !
    ! Evaluate entries per col in A
    !
    nullify (entcolA)
    call memory_alloca(memor,'entcolA','matrix_COO_SPGEMM',entcolA, A % ncols)
    entcolA(:)=0
    do ii =1, size(A % iA,KIND=ip)
       entcolA(A % jA(ii)) = entcolA(A % jA(ii)) + 1
    end do
    !
    ! Evaluate adjacency graph of B and valB
    !
    nullify(adjB,valB)
    call memory_alloca(memor,'adjB','matrix_COO_SPGEMM',adjB, B % nrows)
    call memory_alloca(memor,'valB','matrix_COO_SPGEMM',valB, B % nrows)
    do ii = 1, B % nrows
       call memory_alloca(memor,'adjB(ii) % l','matrix_COO_SPGEMM',adjB(ii) % l, entrowB(ii))
       call memory_alloca(memor,'valB(ii) % a','matrix_COO_SPGEMM',valB(ii) % a, entrowB(ii))
    end do
    entrowB(:)=0
    do ii =1, size(B % iA,KIND=ip)
       entrowB(B % iA(ii)) = entrowB(B % iA(ii)) + 1
       adjB(B % iA(ii)) % l(entrowB(B % iA(ii))) = B % jA(ii)
       valB(B % iA(ii)) % a(entrowB(B % iA(ii))) = B % vA(1,1,ii)
    end do
    !
    ! Evaluate adjacency graph of A and valA
    !
    nullify(adjA,valA)
    call memory_alloca(memor,'adjA','matrix_COO_SPGEMM',adjA, A % ncols)
    call memory_alloca(memor,'valA','matrix_COO_SPGEMM',valA, A % ncols)
    do ii = 1, A % ncols
       call memory_alloca(memor,'adjA(ii) % l','matrix_COO_SPGEMM',adjA(ii) % l, entcolA(ii))
       call memory_alloca(memor,'valA(ii) % a','matrix_COO_SPGEMM',valA(ii) % a, entcolA(ii))
    end do
    entcolA(:)=0
    do ii =1, size(A % iA,KIND=ip)
       entcolA(A % jA(ii)) = entcolA(A % jA(ii)) + 1
       adjA(A % jA(ii)) % l(entcolA(A % jA(ii))) = A % iA(ii)
       valA(A % jA(ii)) % a(entcolA(A % jA(ii))) = A % vA(1,1,ii)
    end do
    !
    ! Evaluate maximum entries per row in C
    !
    nullify (entrowC)
    call memory_alloca(memor,'entrowC','matrix_COO_SPGEMM',entrowC, A % nrows)
    entrowC(:) = 0
    do ii = 1, A % ncols
       do jj = 1, entcolA(ii)
          entrowC(adjA(ii) % l(jj)) = entrowC(adjA(ii) % l(jj)) + entrowB(ii)
       end do
    end do
    !
    ! Evaluate adjacency graph of C
    !
    nullify(adjC,valC)
    call memory_alloca(memor,'adjC','matrix_COO_SPGEMM',adjC, A % nrows)
    call memory_alloca(memor,'valC','matrix_COO_SPGEMM',valC, A % nrows)
    do ii = 1, A % nrows
       call memory_alloca(memor,'adjC(ii) % l','matrix_COO_SPGEMM',adjC(ii) % l, entrowC(ii))
       call memory_alloca(memor,'valC(ii) % a','matrix_COO_SPGEMM',valC(ii) % a, entrowC(ii))
       if(entrowC(ii) > 0) valC(ii) % a(:) = 0
    end do
    entrowC(:) = 0
    nentrC = 0
    do ii = 1, A % ncols
       do jj = 1, entcolA(ii)
          do kk = 1, entrowB(ii)
             itsin = .false.
             do ll = 1, entrowC(adjA(ii) % l(jj))
                if(adjC(adjA(ii) % l(jj)) % l(ll) == adjB(ii) % l(kk) ) then
                   itsin = .true.
                   valC(adjA(ii) % l(jj)) % a(ll) = valC(adjA(ii) % l(jj)) % a(ll) + valA(ii) % a(jj) * valB(ii) % a(kk)
                end if
             end do
             if( itsin .eqv. .false.) then
                entrowC(adjA(ii) % l(jj)) = entrowC(adjA(ii) % l(jj)) + 1
                adjC(adjA(ii) % l(jj)) % l(entrowC(adjA(ii) % l(jj))) = adjB(ii) % l(kk)
                valC(adjA(ii) % l(jj)) % a(entrowC(adjA(ii) % l(jj))) = valA(ii) % a(jj) * valB(ii) % a(kk)
                nentrC = nentrC + 1
             end if
          end do
       end do
    end do
    !
    ! Allocate and store entries in C
    !
    call nullify_spmat(C)
    call memory_alloca(memor,'C','matrix_COO_SPGEMM',C,1_ip,nentrC)
    C % ndof1 = 1
    C % ndof2 = 1
    C % nrows = A % nrows
    C % ncols = B % ncols
    kk = 0
    do ii= 1, C % nrows
       do jj = 1, entrowC(ii)
          kk = kk + 1
          C % iA(kk) = ii
          C % jA(kk) = adjC(ii) % l(jj)
          C % vA(1,1,kk) = valC(ii) % a(jj)
       end do
    end do
    !
    ! Deallocate pointers
    !
    call memory_deallo(memor,'entrowA','matrix_COO_SPGEMM',entcolA)
    call memory_deallo(memor,'entrowB','matrix_COO_SPGEMM',entrowB)
    call memory_deallo(memor,'entrowC','matrix_COO_SPGEMM',entrowC)
    call memory_deallo(memor,'adjA'   ,'matrix_COO_SPGEMM',adjA)
    call memory_deallo(memor,'adjB'   ,'matrix_COO_SPGEMM',adjB)
    call memory_deallo(memor,'adjC'   ,'matrix_COO_SPGEMM',adjC)
    call memory_deallo(memor,'valA'   ,'matrix_COO_SPGEMM',valA)
    call memory_deallo(memor,'valB'   ,'matrix_COO_SPGEMM',valB)
    call memory_deallo(memor,'valC'   ,'matrix_COO_SPGEMM',valC)
    !
    ! Restore memit
    !
    if( present(memit) ) memit = memor

  end subroutine matrix_COO_SPGEMM
  !----------------------------------------------------------------------
  !
  !> @author  Ricard Borrell
  !> @date    11/08/2016
  !> @brief   COO_aggregate_old
  !> @details Sum-up entries of the COO format with same row and col, 
  !           this is slow O(n^2) but does not consume additional mem
  !
  !----------------------------------------------------------------------
  subroutine   matrix_COO_aggregate_old(A,memit)

    implicit none
    type(spmat), intent(inout)              :: A                 !< Sparse matrix
    integer(8),  optional, intent(inout)    :: memit(2)          !< Memory counter
    integer(8)                              :: memor(2)
    type(spmat)                             :: B
    integer(ip)                             :: nentr,icont,jcont,kcont,nentr2
    integer(ip)                             :: lcont,reps,ndof,nrows,ncols

    character(100), PARAMETER :: vacal="matrix_COO_aggregate_old"
    !
    ! Check arguments
    !
    if(present(memit)) then
       memor = memit
    else
       memor = 0
    endif
    !
    ! Check for repeated entries
    !
    if(associated(A % iA)) then
       
       nentr = size(A % iA,KIND=ip)
       ndof  = A % ndof1
       nrows = A % nrows
       ncols = A % ncols
       reps = 0
       do icont = 1,nentr
          do jcont = 1,icont-1
             if(A % iA(icont) == A % iA(jcont) .and. A % jA(icont) == A % jA(jcont)) then
                do kcont = 1,ndof
                   do lcont = 1,ndof
                      A % vA(kcont,lcont,jcont) = A % vA(kcont,lcont,jcont) + A % vA(kcont,lcont,icont)
                   enddo
                enddo
                A % iA(icont) = -1_ip
                A % jA(icont) = -1_ip
                reps = reps + 1
             endif
          enddo
       enddo
       if(reps > 0) then
          nentr2=nentr-reps
          nullify(B%iA,B%jA,B%vA)
          call memory_alloca(memor,'B',vacal,B,ndof,nentr2)
          jcont = 1
          do icont = 1,nentr
             if(A % iA(icont) /= -1_ip .and. A % jA(icont) /= -1_ip) then
                B % iA(jcont) = A % iA(icont)
                B % jA(jcont) = A % jA(icont)
                B % vA(1:ndof,1:ndof,jcont) = A % vA(1:ndof,1:ndof,icont)
                jcont = jcont + 1
             endif
          enddo
          call memory_deallo(memor,'A',vacal,A)
          call memory_alloca(memor,'A',vacal,A,ndof,nentr2)
          A % ndof1 = ndof
          A % ndof2 = ndof
          A % nrows = nrows
          A % ncols = ncols
          A % iA(1:nentr2) = B % iA(1:nentr2)
          A % jA(1:nentr2) = B % jA(1:nentr2)
          A % vA(1:ndof,1:ndof,1:nentr2) = B % vA(1:ndof,1:ndof,1:nentr2)
          call memory_deallo(memor,'B',vacal,B)
       endif
    endif
    !
    ! Restore memit
    !
    if( present(memit) ) memit = memor

  end subroutine matrix_COO_aggregate_old
  
  !----------------------------------------------------------------------
  !
  !> @author  Ricard Borrell
  !> @date    11/08/2016
  !> @brief   COO_aggregate
  !> @details Sum-up entries of the COO format with same row and col.
  !           This option id faster than the previos but consumes
  !           more main memory (needs store the matrix two times).  
  !
  !----------------------------------------------------------------------
  subroutine   matrix_COO_aggregate(A,memit)

    implicit none
    type(spmat), intent(inout)              :: A                 !< Sparse matrix
    integer(8),  optional, intent(inout)    :: memit(2)          !< Memory counter
    integer(8)                              :: memor(2)
    type(spmat)                             :: B
    integer(ip)                             :: nentr,icont,jcont,kcont,nentr2
    integer(ip)                             :: lcont,reps,ndof,nrows,ncols,ientr
    integer(ip)                             :: icol,irow
    integer(ip), pointer                    :: nrowent(:)
    integer(ip), pointer                    :: controw(:)
    integer(ip), pointer                    :: ncolent(:)
    type(i1p),   pointer                    :: rowent(:)
   

    character(100), PARAMETER :: vacal="matrix_COO_aggregate"
    !
    ! Check arguments
    !
    if(present(memit)) then
       memor = memit
    else
       memor = 0
    endif
    !
    ! Check for repeated entries
    !
    if(associated(A % iA)) then
    
   
       nentr = size(A % iA,KIND=ip)
       ndof  = A % ndof1
       nrows = A % nrows
       ncols = A % ncols
       nullify(nrowent,ncolent,rowent,controw)
       call memory_alloca(memor,'NROWENT',vacal,nrowent,nrows)
       call memory_alloca(memor,'CONTROW',vacal,controw,nrows)
       call memory_alloca(memor,'ROWENT' ,vacal,rowent ,nrows)
       call memory_alloca(memor,'NCOLENT',vacal,ncolent,ncols)
       nrowent = 0_ip
       ncolent = 0_ip
       controw = 0_ip
       !
       ! Count entries per row and entries per column
       !
       do icont = 1,nentr
          nrowent(A % iA(icont)) = nrowent(A % iA(icont)) + 1_ip
          ncolent(A % jA(icont)) = ncolent(A % jA(icont)) + 1_ip
       enddo
       !
       ! The rows with more than one entry are stored in rowent
       !
       do icont = 1,nrows
          if(nrowent(icont) > 1_ip) then
             nullify(rowent(icont) % l)
             call memory_alloca(memor,'ROWENT % L',vacal,rowent(icont) % l,nrowent(icont))
          endif
       enddo
       do icont = 1,nentr
          if(nrowent(A % iA(icont)) > 1_ip) then
             controw(A % iA(icont)) = controw(A % iA(icont)) + 1_ip
             rowent(A % iA(icont)) % l(controw(A % iA(icont))) = icont
          endif
       enddo
       !
       ! Check reps
       !
       reps = 0_ip
       do icont = 1,nentr
          irow = A % iA(icont)
          icol = A % jA(icont)
          if(irow /= -1) then
             if(nrowent(irow) > 1_ip .and. ncolent(icol) > 1_ip ) then
                do jcont = 1,nrowent(irow)
                   ientr = rowent(irow) % l(jcont)
                   !
                   ! Aculamos sobre la primera entrada
                   !
                   if(A % iA(icont) == A % iA(ientr) .and. A % jA(icont) == A % jA(ientr) .and. icont/=ientr) then
                      do kcont = 1,ndof
                         do lcont = 1,ndof
                            A % vA(kcont,lcont,icont) = A % vA(kcont,lcont,icont) + A % vA(kcont,lcont,ientr) 
                         enddo
                      enddo
                      A % iA(ientr) = -1_ip
                      A % jA(ientr) = -1_ip
                      reps = reps + 1
                      nrowent(irow) = nrowent(irow) - 1_ip
                      ncolent(icol) = ncolent(icol) - 1_ip
                   endif
                enddo
             endif
          endif
       enddo
       !
       ! Elminate reps
       !
       if(reps > 0_ip) then
          nentr2=nentr-reps
          nullify(B%iA,B%jA,B%vA)
          call memory_alloca(memor,'B',vacal,B,ndof,nentr2)
          jcont = 1
          do icont = 1,nentr
             if(A % iA(icont) /= -1_ip .and. A % jA(icont) /= -1_ip) then
                B % iA(jcont) = A % iA(icont)
                B % jA(jcont) = A % jA(icont)
                B % vA(1:ndof,1:ndof,jcont) = A % vA(1:ndof,1:ndof,icont)
                jcont = jcont + 1
             endif
          enddo 
          call memory_deallo(memor,'A',vacal,A)
          call memory_alloca(memor,'A',vacal,A,ndof,nentr2)
          A % ndof1 = ndof
          A % ndof2 = ndof
          A % nrows = nrows
          A % ncols = ncols
          A % iA(1:nentr2) = B % iA(1:nentr2)
          A % jA(1:nentr2) = B % jA(1:nentr2)
          A % vA(1:ndof,1:ndof,1:nentr2) = B % vA(1:ndof,1:ndof,1:nentr2)
          call memory_deallo(memor,'B',vacal,B)

       endif
       call memory_deallo(memor,'NROWENT',vacal,nrowent)
       call memory_deallo(memor,'CONTROW',vacal,controw)
       call memory_deallo(memor,'ROWENT' ,vacal,rowent)
       call memory_deallo(memor,'NCOLENT',vacal,ncolent)
    endif
    !
    ! Restore memit
    !
    if( present(memit) ) memit = memor

  end subroutine matrix_COO_aggregate
  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    01/03/2016
  !> @brief   CSR to skyline
  !> @details Convert a matrix from CSR to skyline format
  !
  !-----------------------------------------------------------------------

  subroutine matrix_csr_to_skyline(nbnodes,ndofn,ia,ja,iskyl,idiag,a_csr,a_sky)

    integer(ip), intent(in)  :: nbnodes
    integer(ip), intent(in)  :: ndofn
    integer(ip), intent(in)  :: ia(*)
    integer(ip), intent(in)  :: ja(*)
    integer(ip), intent(in)  :: iskyl(*)
    integer(ip), intent(in)  :: idiag(*)
    real(rp),    intent(in)  :: a_csr(ndofn,ndofn,*)
    real(rp),    intent(out) :: a_sky(*)
    integer(ip)              :: ii,izdom,jj,idofn,jdofn
    integer(ip)              :: ii1,jj1,kskyl

    do ii = 1,nbnodes

       do izdom = ia(ii),ia(ii+1)-1
          jj = ja(izdom)

          do idofn = 1,ndofn
             do jdofn = 1,ndofn

                ii1 = (ii-1) * ndofn + idofn
                jj1 = (jj-1) * ndofn + jdofn

                if( ii1 < jj1 ) then
                   kskyl        = iskyl(jj1+1)-(jj1-ii1)
                   a_sky(kskyl) = a_csr(idofn,jdofn,izdom)
                else
                   kskyl        = idiag(ii1)-(ii1-jj1)
                   a_sky(kskyl) = a_csr(idofn,jdofn,izdom)
                endif

             end do
          end do
       end do
    end do

  end subroutine matrix_csr_to_skyline

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    01/03/2016
  !> @brief   CSR to skyline
  !> @details Convert a matrix from CSR to skyline format
  !
  !-----------------------------------------------------------------------

  subroutine matrix_transpose(nbnodes,ndof1,ndof2,ia,ja,aa_in,aa_out,xscal_opt)

    integer(ip), intent(in)              :: nbnodes
    integer(ip), intent(in)              :: ndof1
    integer(ip), intent(in)              :: ndof2
    integer(ip), intent(in)              :: ia(*)
    integer(ip), intent(in)              :: ja(*)
    real(rp),    intent(inout)           :: aa_in(ndof1,ndof2,*)
    real(rp),    intent(out),   optional :: aa_out(ndof1,ndof2,*)
    real(rp),    intent(in),    optional :: xscal_opt
    integer(ip)                          :: ii,jj,iz,jz
    real(rp)                             :: xscal

    if( present(xscal_opt) ) then
       xscal = xscal_opt
    else
       xscal = 1.0_rp
    end if

    do ii = 1,nbnodes
       do iz = ia(ii),ia(ii+1)-1
          jj = ja(iz)
          loop_jz: do jz = ia(jj),ia(jj+1)-1
             if( ja(jz) == ii ) then
                aa_out(:,:,jz) = xscal * aa_in(:,:,iz)
                exit loop_jz
             end if
          end do loop_jz
       end do
    end do

  end subroutine matrix_transpose

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    01/03/2016
  !> @brief   CSR to skyline
  !> @details Convert a matrix from CSR to skyline format
  !
  !-----------------------------------------------------------------------

  subroutine matrix_full_to_half(nbnodes,ndofn,ia,ja,aa_in,aa_out)

    integer(ip), intent(in)              :: nbnodes
    integer(ip), intent(in)              :: ndofn
    integer(ip), intent(in)              :: ia(*)
    integer(ip), intent(in)              :: ja(*)
    real(rp)                             :: aa_in(ndofn,ndofn,*)
    real(rp),    intent(inout), optional :: aa_out(ndofn,ndofn,*)
    integer(ip)                          :: ii,jj,iz,kz

    if( present(aa_out) ) then
       kz = 0
       do ii = 1,nbnodes
          do iz = ia(ii),ia(ii+1)-1
             jj = ja(iz)
             if( jj >= ii ) then
                kz = kz + 1
                aa_out(1:ndofn,1:ndofn,kz) = aa_in(1:ndofn,1:ndofn,iz)
             end if
          end do
       end do
    else
       kz = 0
       do ii = 1,nbnodes
          do iz = ia(ii),ia(ii+1)-1
             jj = ja(iz)
             if( jj >= ii ) then
                kz = kz + 1
                aa_in(1:ndofn,1:ndofn,kz) = aa_in(1:ndofn,1:ndofn,iz)
             end if
          end do
       end do
    end if

  end subroutine matrix_full_to_half

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    01/03/2016
  !> @brief   CSR to skyline
  !> @details Convert a matrix from BCSR to COO format
  !
  !-----------------------------------------------------------------------

  subroutine matrix_csr_to_coo(nbnodes,ndofn,ia,ja,aa_in,aa_out,message)

    integer(ip), intent(in)              :: nbnodes
    integer(ip), intent(in)              :: ndofn
    integer(ip), intent(in)              :: ia(*)
    integer(ip), intent(in)              :: ja(*)
    real(rp)                             :: aa_in(ndofn,ndofn,*)
    real(rp),    intent(inout), optional :: aa_out(*)
    character(*),intent(in),    optional :: message
    real(rp),    allocatable             :: aa_tmp(:,:,:)
    integer(ip)                          :: ii,jj,iz,kz,idofn,jdofn
    integer(ip)                          :: ii_idofn,jj_jdofn
    logical(lg)                          :: symmetry
    !
    ! Symmetry
    !
    symmetry = .false.
    if( present(message) ) then
       if( trim(message) == 'SYMMETRIC' ) symmetry = .true.
    end if

    if( present(aa_out) ) then

       if( symmetry ) then
          if( ndofn == 1 ) then
             kz = 0
             do ii = 1,nbnodes
                do iz = ia(ii),ia(ii+1)-1
                   jj = ja(iz)
                   if( jj >= ii ) then
                      kz         = kz + 1
                      aa_out(kz) = aa_in(1,1,iz)
                   end if
                end do
             end do
          else
             kz = 0
             do ii = 1,nbnodes
                do idofn = 1,ndofn
                   ii_idofn = (ii-1)*ndofn+idofn
                   do iz = ia(ii),ia(ii+1)-1
                      jj = ja(iz)
                      do jdofn = 1,ndofn
                         jj_jdofn   = (jj-1)*ndofn+jdofn
                         if( jj_jdofn >= ii_idofn ) then
                            kz         = kz + 1
                            aa_out(kz) = aa_in(jdofn,idofn,iz)
                         end if
                      end do
                   end do
                end do
             end do
          end if
       else
          if( ndofn == 1 ) then
             kz = 0
             do ii = 1,nbnodes
                do iz = ia(ii),ia(ii+1)-1
                   kz         = kz + 1
                   aa_out(kz) = aa_in(1,1,iz)
                end do
             end do
          else
             kz = 0
             do ii = 1,nbnodes
                do idofn = 1,ndofn
                   do iz = ia(ii),ia(ii+1)-1
                      jj = ja(iz)
                      do jdofn = 1,ndofn
                         kz         = kz + 1
                         aa_out(kz) = aa_in(jdofn,idofn,iz)
                      end do
                   end do
                end do
             end do
          end if
       end if

    else

       kz = ia(nbnodes+1)-1
       allocate(aa_tmp(ndofn,ndofn,kz))
       aa_tmp(1:ndofn,1:ndofn,1:kz) = aa_in(1:ndofn,1:ndofn,1:kz)
       kz = 0
       do ii = 1,nbnodes
          do idofn = 1,ndofn
             do iz = ia(ii),ia(ii+1)-1
                jj = ja(iz)
                do jdofn = 1,ndofn
                   kz            = kz + 1
                   aa_in(1,1,kz) = aa_tmp(jdofn,idofn,iz)
                end do
             end do
          end do
       end do
       deallocate(aa_tmp)

    end if

  end subroutine matrix_csr_to_coo

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    04/04/2017
  !> @brief   Partial to full row
  !> @details Copy a partial (square) matrix AA_PARTIAL
  !>          to a full row (rectangular) matrix AA_FULL. It is assumed
  !>          that the rows of the full row matrix coincide with that
  !>          of the partial one. This is ensured by numbering first
  !>          the own nodes in Alya.
  !
  !-----------------------------------------------------------------------

  subroutine matrix_partial_to_full_row_matrix(&
       nbnodes_full,ndofn,&
       ia_partial,ja_partial,aa_partial,&
       ia_full,ja_full,aa_full)

    integer(ip), intent(in)              :: nbnodes_full
    integer(ip), intent(in)              :: ndofn
    integer(ip), intent(in)              :: ia_partial(*)
    integer(ip), intent(in)              :: ja_partial(*)
    real(rp),    intent(in)              :: aa_partial(ndofn,ndofn,*)
    integer(ip), intent(in)              :: ia_full(*)
    integer(ip), intent(in)              :: ja_full(*)
    real(rp),    intent(out)             :: aa_full(ndofn,ndofn,*)

    integer(ip)                          :: ii,jj,iz,jz,nz_full

    nz_full    = ia_full   (nbnodes_full+1)    - 1
    aa_full(1:ndofn,1:ndofn,1:nz_full) = 0.0_rp

    do ii = 1,nbnodes_full
       do iz = ia_full(ii),ia_full(ii+1)-1
          jj = ja_full(iz)
          jz = ia_partial(ii)
          do while( jz <= ia_partial(ii+1)-1 )
             if( ja_partial(jz) == jj ) then
                aa_full(1:ndofn,1:ndofn,iz) = aa_partial(1:ndofn,1:ndofn,jz)
                jz = ia_partial(ii+1)
             else
                jz = jz + 1
             end if
          end do
       end do
    end do

  end subroutine matrix_partial_to_full_row_matrix

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/10/2014
  !> @brief   heap sort of the rows of a matrix
  !> @details Quick sorting of each row of a matrix using values of JA
  !>          ITASK = 1 ... Decreasing value, i.e., JA(1) > JA(2) > ...
  !>          ITASK = 2 ... Increasing value, i.e., JA(1) < JA(2) < ...
  !>          Then, the matrix in CSR format
  !
  !----------------------------------------------------------------------

  subroutine matrix_heap_sort(itask,ndofn,npoin,ia,ja,aa)
    integer(ip),  intent(in)     :: itask
    integer(ip),  intent(in)     :: ndofn
    integer(ip),  intent(in)     :: npoin
    integer(ip),  intent(in)     :: ia(*)
    integer(ip),  intent(inout)  :: ja(*)
    real(rp),     intent(inout)  :: aa(ndofn,ndofn,*)
    integer(ip)                  :: nrows,ipoin,iz

    do ipoin = 1,npoin
       nrows = ia(ipoin+1)-ia(ipoin)
       iz    = ia(ipoin)
       call matrix_heap_sort_single_row(itask,ndofn,nrows,ja(iz),aa(:,:,iz))
    end do

  end subroutine matrix_heap_sort

  subroutine matrix_heap_sort_single_row(itask,ndofn,nrows,ja,aa)

    integer(ip),  intent(in)    :: itask
    integer(ip),  intent(in)    :: ndofn
    integer(ip),  intent(inout) :: nrows
    integer(ip),  intent(inout) :: ja(*)
    real(rp),     intent(inout) :: aa(ndofn,ndofn,*)
    integer(ip)                 :: leni,ir,ii,jj,iaux
    real(rp)                    :: aa_aux(ndofn,ndofn)

    select case ( itask )

    case ( 1_ip )
       !
       ! Decreasing order
       !
       call runend('matrix_heap_sort: NOT CODED')

    case ( 2_ip )
       !
       ! Increasing order
       !
       if( nrows < 2 ) then
          goto 500
       end if

       leni = (nrows/2) + 1
       ir  = nrows

300    continue

       if( leni > 1 ) then
          leni         = leni - 1
          iaux        = ja(leni)
          aa_aux(:,:) = aa(:,:,leni)
       else
          iaux        = ja(ir)
          ja(ir)      = ja(1)
          aa_aux(:,:) = aa(:,:,ir)
          aa(:,:,ir)  = aa(:,:,1)

          ir = ir - 1

          if( ir == 1 ) then
             ja(1)     = iaux
             aa(:,:,1) = aa_aux(:,:)
             goto 500
          end if
       end if

       ii = leni
       jj = leni + leni

400    if( jj <= ir ) then
          if( jj < ir ) then
             if ( ja(jj) < ja(jj+1) ) then
                jj = jj + 1
             end if
          end if

          if( iaux < ja(jj) ) then
             ja(ii)     = ja(jj)
             aa(:,:,ii) = aa(:,:,jj)
             ii = jj
             jj = jj + jj
          else
             jj = ir + 1
          end if

          goto 400
       end if

       ja(ii)     = iaux
       aa(:,:,ii) = aa_aux(:,:)

       goto 300

    end select

500 continue

  end subroutine matrix_heap_sort_single_row

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    19/07/2017
  !> @brief   Remove zero from a matrix
  !> @details Remove the null coefficients from a matrix and renumber
  !>          the graph
  !
  !----------------------------------------------------------------------

  subroutine matrix_remove_null_coefficients(ndofn,nrows,nnz,ia,ja,aa)

    integer(ip),  intent(in)     :: ndofn
    integer(ip),  intent(in)     :: nrows
    integer(ip),  intent(out)    :: nnz
    integer(ip),  intent(inout)  :: ia(*)
    integer(ip),  intent(inout)  :: ja(*)
    real(rp),     intent(inout)  :: aa(ndofn,ndofn,*)
    integer(ip)                  :: ipoin,iz,kk,ii,ll,iz1,iz2
    real(rp),     parameter      :: epsil = epsilon(1.0_rp)

    nnz = 0
    do ipoin = 1,nrows
       iz1       = ia(ipoin)
       iz2       = ia(ipoin+1)-1
       ia(ipoin) = 0
       do iz = iz1,iz2
          if( maxval(abs(aa(:,:,iz))) > epsil ) then
             nnz         = nnz + 1
             ia(ipoin)   = ia(ipoin) + 1
             ja(nnz)     = ja(iz)
             aa(:,:,nnz) = aa(:,:,iz)
          end if
       end do
    end do
    !
    ! Number to linked list
    !
    kk    = ia(1)
    ia(1) = 1
    do ii = 2,nrows+1
       ll     = ia(ii)
       ia(ii) = ia(ii-1) + kk
       kk     = ll
    end do

  end subroutine matrix_remove_null_coefficients

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    19/07/2017
  !> @brief   Sparisfy a matrix
  !> @details sparsify a matrix by removing coefficients Aij such that
  !>          |Aij| < eps * 1/2 * ( |Aii| + |Ajj| )
  !>          If eps is huge ...... only the diagonal is assembled
  !>          If eps <= 1.0e-16 ... zero coefficients are removed
  !>
  !----------------------------------------------------------------------

  subroutine matrix_sparsify(eps,nbrows,ndof,nzin,iAin,jAin,Ain,nzout,iAout,jAout,Aout,memor_opt)

    real(rp),    intent(in)                       :: eps
    integer(ip), intent(in)                       :: nbrows            !< Number of rows
    integer(ip), intent(in)                       :: ndof              !< Number of dof per row
    integer(ip), intent(inout)                    :: nzin              !< Size of graph
    integer(ip), intent(inout), pointer           :: iAin(:)           !< Matrix graph
    integer(ip), intent(inout), pointer           :: jAin(:)           !< Matrix graph
    real(rp),    intent(inout), pointer           :: Ain(:,:,:)        !< Input matrix
    integer(ip), intent(out),            optional :: nzout             !< Output matrix size
    real(rp),    intent(inout), pointer, optional :: Aout(:,:,:)       !< Output matrix
    integer(ip), intent(inout), pointer, optional :: iAout(:)          !< Output graph
    integer(ip), intent(inout), pointer, optional :: jAout(:)          !< Output graph
    integer(8),  intent(inout),          optional :: memor_opt(2)
    integer(ip)                                   :: ii,jj,iz,jz,itask
    integer(8)                                    :: memor(2)
    real(rp)                                      :: Aii(ndof,ndof),eps2

    integer(ip)                                   :: nz
    real(rp),    pointer                          :: A(:,:,:)
    integer(ip), pointer                          :: iA(:)
    integer(ip), pointer                          :: jA(:)

    nullify(A)
    nullify(iA)
    nullify(jA)
    if( present(memor_opt) ) memor = memor_opt
    eps2 = 0.5_rp * eps

    if( ndof == 1 ) then
       do itask = 1,2
          nz = 0
          do ii = 1,nbrows
             !
             ! Diagonal block
             !
             call graphs_find_edge(ii,ii,iAin,jAin,iz)
             if( iz == 0 ) then
                print*,'IZ IS ZERO!!!!!!!!!!!!!!!!!!!!!!'
                stop
             end if
             Aii(1:ndof,1:ndof) = Ain(1:ndof,1:ndof,iz)

             do iz = iAin(ii),iAin(ii+1)-1
                jj = jAin(iz)
                if( ii /= jj ) then
                   call graphs_find_edge(jj,jj,iAin,jAin,jz)
                   if( abs(Ain(1,1,iz)) <= eps2 * ( abs(Aii(1,1)) + abs(Ain(1,1,jz)) ) ) then
                      continue
                   else
                      nz  = nz  + 1
                      if( itask == 2 ) then
                         if( present(nzout) ) then
                            iAout(ii)    = iAout(ii) + 1
                            jAout(nz)    = jj
                            Aout(1,1,nz) = Ain(1,1,iz)
                         else
                            iA(ii)       = iA(ii) + 1
                            jA(nz)       = jj
                            A(1,1,nz)    = Ain(1,1,iz)
                         end if
                      end if
                   end if
                else
                   nz  = nz  + 1
                   if( itask == 2 ) then
                      if( present(nzout) ) then
                         iAout(ii)    = iAout(ii) + 1
                         jAout(nz)    = jj
                         Aout(1,1,nz) = Ain(1,1,iz)
                      else
                         iA(ii)       = iA(ii) + 1
                         jA(nz)       = jj
                         A(1,1,nz)    = Ain(1,1,iz)
                      end if
                   end if
                end if
             end do
          end do
          if( itask == 1 ) then
             if( present(nzout) ) then
                call memory_alloca(memor,'AOUT' ,'matrix_sparsify', Aout,ndof,ndof,nz)
                call memory_alloca(memor,'IAOUT','matrix_sparsify',iAout,nbrows+1)
                call memory_alloca(memor,'JAOUT','matrix_sparsify',jAout,nz)
             else
                call memory_alloca(memor,'A'    ,'matrix_sparsify', A,ndof,ndof,nz)
                call memory_alloca(memor,'IA'   ,'matrix_sparsify',iA,nbrows+1)
                call memory_alloca(memor,'JA'   ,'matrix_sparsify',jA,nz)
             end if
          end if
       end do
    else
       call runend('MATRX_SPARSIFY: NODFN > 1 NOT CODED')
    end if

    if( present(nzout) ) then
       nzout = nz
       call graphs_number_to_linked_list(nbrows,iAout)
    else
       nzin = nz
       call graphs_number_to_linked_list(nbrows,iA)
       call memory_copy  (memor,'iAin','matrix_sparsify',iA,iAin)
       call memory_resize(memor,'jAin','matrix_sparsify',jAin,nzin)
       call memory_copy  (memor,'Ain' ,'matrix_sparsify',jA,jAin)
       call memory_resize(memor,'Ain' ,'matrix_sparsify',Ain,ndof,ndof,nzin)
       call memory_copy  (memor,'Ain' ,'matrix_sparsify',A,Ain)
    end if

    if( present(memor_opt) ) memor_opt = memor

  end subroutine matrix_sparsify

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    20/11/2017
  !> @brief   SpMV in ELL format
  !> @details Sparse matrix vector product in ELL format
  !
  !----------------------------------------------------------------------

  subroutine matrix_ELL_SpMV(n1,n2,ndofr,ndofc,iA,An,xx,yy,OPENMP,INITIALIZATION,CHUNK,SCHEDULE)
    !
    ! Dummy arguments
    !
    integer(ip),  intent(in)           :: n1                !< Initial node
    integer(ip),  intent(in)           :: n2                !< Final node
    integer(ip),  intent(in)           :: ndofr             !< Number of dof per node
    integer(ip),  intent(in)           :: ndofc             !< Number of dof per node
    integer(ip),  intent(in), pointer  :: iA(:,:)           !< Matrix graph ia
    real(rp),     intent(in)           :: An(ndofc,ndofr,*) !< Matrix
    real(rp),     intent(in)           :: xx(ndofc,*)       !< Input vector
    real(rp),     intent(out)          :: yy(ndofr,*)       !< Output vector
    logical(lg),  intent(in), optional :: OPENMP            !< If openmp chould be used
    logical(lg),  intent(in), optional :: INITIALIZATION    !< If results should be initialized
    integer(ip),  intent(in), optional :: CHUNK             !< Chunks for dynamic scheduling of OpenMP
    character(*), intent(in), optional :: SCHEDULE          !< OpenMP schedule
    integer(ip)                        :: ii,jj,iz,ncols
    integer(ip)                        :: icol,idof,jdof
    logical(lg)                        :: use_openmp
    logical(lg)                        :: do_initialize
    integer(ip)                        :: my_schedule
    integer(ip)                        :: my_chunk

    use_openmp    = .false.
    do_initialize = .true.
    
    if( present(OPENMP)         ) use_openmp    = OPENMP
    if( present(INITIALIZATION) ) do_initialize = INITIALIZATION
    !
    ! OpenMP options
    !        
    if( use_openmp ) then

       my_chunk    = spmv_chunk
       my_schedule = SOL_OMP_STATIC
       
       if( present(CHUNK) ) my_chunk = CHUNK
       if( present(SCHEDULE) ) then
          if(      trim(SCHEDULE) == 'STATIC' ) then
             my_schedule = SOL_OMP_STATIC
          else if( trim(SCHEDULE) == 'DYNAMIC' ) then
             if( my_chunk <= 1 ) then
                my_schedule = SOL_OMP_STATIC
             else
                my_schedule = SOL_OMP_DYNAMIC
             end if
          else if( trim(SCHEDULE) == 'GUIDED' ) then
             my_schedule = SOL_OMP_GUIDED
          end if
       end if

    end if
    !
    ! Initialize solution
    !
    if( do_initialize ) then
       yy(1:ndofr,n1:n2) = 0.0_rp
    end if

    ncols = size(iA,1,KIND=ip)

    if( ndofr == 1 .and. ndofc == 1 ) then
       if( use_openmp ) then
          if( my_schedule == SOL_OMP_STATIC ) then
             !$OMP PARALLEL   DO                              &
             !$OMP SCHEDULE ( STATIC )                        &
             !$OMP DEFAULT  ( NONE )                          &
             !$OMP SHARED   ( an, ia, ncols, n1, n2, xx, yy ) &
             !$OMP PRIVATE  ( icol, ii, jj, iz )
             do ii = n1,n2
                do icol = 1,ncols
                   jj       = iA(icol,ii)
                   iz       = (ii-1)*ncols+icol
                   yy(1,ii) = yy(1,ii) + An(1,1,iz) * xx(1,jj)
                end do
             end do
             !$OMP END PARALLEL DO
          else if( my_schedule == SOL_OMP_GUIDED ) then
             !$OMP PARALLEL   DO                              &
             !$OMP SCHEDULE ( GUIDED )                        &
             !$OMP DEFAULT  ( NONE )                          &
             !$OMP SHARED   ( an, ia, ncols, n1, n2, xx, yy ) &
             !$OMP PRIVATE  ( icol, ii, jj, iz )
             do ii = n1,n2
                do icol = 1,ncols
                   jj       = iA(icol,ii)
                   iz       = (ii-1)*ncols+icol
                   yy(1,ii) = yy(1,ii) + An(1,1,iz) * xx(1,jj)
                end do
             end do
             !$OMP END PARALLEL DO
          else if( my_schedule == SOL_OMP_DYNAMIC ) then
             !$OMP PARALLEL   DO                              &
             !$OMP SCHEDULE ( DYNAMIC , my_chunk )            &
             !$OMP DEFAULT  ( NONE )                          &
             !$OMP SHARED   ( an, ia, ncols, n1, n2, xx, yy ) &
             !$OMP PRIVATE  ( icol, ii, jj, iz )
             do ii = n1,n2
                do icol = 1,ncols
                   jj       = iA(icol,ii)
                   iz       = (ii-1)*ncols+icol
                   yy(1,ii) = yy(1,ii) + An(1,1,iz) * xx(1,jj)
                end do
             end do
             !$OMP END PARALLEL DO
          end if
       else
          do ii = n1,n2
             do icol = 1,ncols
                jj       = iA(icol,ii)
                iz       = (ii-1)*ncols+icol
                yy(1,ii) = yy(1,ii) + An(1,1,iz) * xx(1,jj)
             end do
          end do
       end if
    else
       if( use_openmp ) then
          if( my_schedule == SOL_OMP_STATIC ) then
             !$OMP PARALLEL   DO                                            &
             !$OMP SCHEDULE ( STATIC )                                      &
             !$OMP DEFAULT  ( NONE )                                        &
             !$OMP SHARED   ( an, ia, ncols, ndofr, ndofc, n1, n2, xx, yy ) &
             !$OMP PRIVATE  ( icol, ii, jj, iz, idof, jdof )
             do ii = n1,n2
                do idof = 1,ndofr
                   do icol = 1,ncols
                      jj       = iA(icol,ii)
                      iz       = (ii-1)*ncols+icol
                      do jdof = 1,ndofc
                         yy(idof,ii) = yy(idof,ii) + An(jdof,idof,iz) * xx(jdof,jj)
                      end do
                   end do
                end do
             end do
             !$OMP END PARALLEL DO
          else if( my_schedule == SOL_OMP_GUIDED ) then
             !$OMP PARALLEL   DO                                            &
             !$OMP SCHEDULE ( GUIDED )                                      &
             !$OMP DEFAULT  ( NONE )                                        &
             !$OMP SHARED   ( an, ia, ncols, ndofr, ndofc, n1, n2, xx, yy ) &
             !$OMP PRIVATE  ( icol, ii, jj, iz, idof, jdof )
             do ii = n1,n2
                do idof = 1,ndofr
                   do icol = 1,ncols
                      jj       = iA(icol,ii)
                      iz       = (ii-1)*ncols+icol
                      do jdof = 1,ndofc
                         yy(idof,ii) = yy(idof,ii) + An(jdof,idof,iz) * xx(jdof,jj)
                      end do
                   end do
                end do
             end do
             !$OMP END PARALLEL DO
          else if( my_schedule == SOL_OMP_DYNAMIC ) then
             !$OMP PARALLEL   DO                                            &
             !$OMP SCHEDULE ( DYNAMIC , my_chunk )                          &
             !$OMP DEFAULT  ( NONE )                                        &
             !$OMP SHARED   ( an, ia, ncols, ndofr, ndofc, n1, n2, xx, yy ) &
             !$OMP PRIVATE  ( icol, ii, jj, iz, idof, jdof )
             do ii = n1,n2
                do idof = 1,ndofr
                   do icol = 1,ncols
                      jj       = iA(icol,ii)
                      iz       = (ii-1)*ncols+icol
                      do jdof = 1,ndofc
                         yy(idof,ii) = yy(idof,ii) + An(jdof,idof,iz) * xx(jdof,jj)
                      end do
                   end do
                end do
             end do
             !$OMP END PARALLEL DO
          end if
       else
          do ii = n1,n2
             do idof = 1,ndofr
                do icol = 1,ncols
                   jj       = iA(icol,ii)
                   iz       = (ii-1)*ncols+icol
                   do jdof = 1,ndofc
                      yy(idof,ii) = yy(idof,ii) + An(jdof,idof,iz) * xx(jdof,jj)
                   end do
                end do
             end do
          end do
       end if
    end if

  end subroutine matrix_ELL_SpMV

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    04/02/2016
  !> @brief   Assemble element matrix
  !> @details Assemble an element matrix ELMAT into the global matrix
  !>          AN with CSR format
  !>
  !----------------------------------------------------------------------

#ifndef _OPENMP  
  pure subroutine matrix_assemble_element_matrix_to_COO(&
       kfl_element_to_csr,nbvar,pnode,pevat,&
       ielem,lnods,elmat,rows,cols,an,element_to_csr)
#else    
  subroutine matrix_assemble_element_matrix_to_COO(&
       kfl_element_to_csr,nbvar,pnode,pevat,&
       ielem,lnods,elmat,rows,cols,an,element_to_csr)
#endif  

    integer(ip), intent(in)                      :: kfl_element_to_csr
    integer(ip), intent(in)                      :: nbvar
    integer(ip), intent(in)                      :: pnode
    integer(ip), intent(in)                      :: pevat
    integer(ip), intent(in)                      :: ielem
    integer(ip), intent(in)                      :: lnods(pnode)
    integer(ip), intent(in),   pointer           :: rows(:)
    integer(ip), intent(in),   pointer           :: cols(:)
    real(rp),    intent(in)                      :: elmat(pevat,pevat)
    real(rp),    intent(inout)                   :: an(nbvar,nbvar,*)
    integer(ip), intent(in),   pointer, optional :: element_to_csr(:,:,:)
    integer(ip)                                  :: ievat,jevat,idofn,jdofn
    integer(ip)                                  :: inode,jnode
    integer(ip)                                  :: ipoin,jpoin,izsol,izsol_init
    logical(lg)                                  :: use_element_to_csr

    !if( size(rows,KIND=ip) /= size(cols,KIND=ip) ) call runend('matrix_assemble_element_matrix_to_COO: WRONG COO FORMAT')
    if( kfl_element_to_csr == 1.and. present(element_to_csr) ) then
       use_element_to_csr = .true.
    else
       use_element_to_csr = .false.
    end if

    if( use_element_to_csr ) then

       do inode = 1,pnode
          do jnode = 1,pnode
             izsol = element_to_csr(inode,jnode,ielem)
             do idofn = 1,nbvar
                ievat = (inode-1) * nbvar + idofn
                do jdofn = 1,nbvar
                   jevat = (jnode-1) * nbvar + jdofn
#ifdef NO_COLORING
                   !$OMP ATOMIC
#endif
                   an(jdofn,idofn,izsol) = an(jdofn,idofn,izsol) + elmat(ievat,jevat)
                end do
             end do
          end do
       end do

    else

       do inode = 1,pnode
          ipoin = lnods(inode)
          izsol_init = 1
          do while( ipoin /= rows(izsol_init) )
             izsol_init = izsol_init + 1
             !if( izsol_init > size(rows,KIND=ip) ) call runend('COO FORMAT IN TROUBLE')
          end do
          do jnode = 1,pnode
             jpoin = lnods(jnode)
             izsol = izsol_init
             do while( jpoin /= cols(izsol) )
                izsol = izsol + 1
                !if( izsol > size(cols,KIND=ip) ) call runend('COO FORMAT IN TROUBLE')
             end do

             do idofn = 1,nbvar
                ievat = (inode-1) * nbvar + idofn
                do jdofn = 1,nbvar
                   jevat = (jnode-1) * nbvar + jdofn
#ifdef NO_COLORING
                   !$OMP ATOMIC
#endif
                   an(jdofn,idofn,izsol) = an(jdofn,idofn,izsol) + elmat(ievat,jevat)
                end do
             end do
          end do
       end do
    end if

  end subroutine matrix_assemble_element_matrix_to_COO

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   SpMV
  !> @details Sparse matrix vector product in BCSR format
  !>          y = A x for rows n1:n2
  !>          Initial and end graphs can be different too
  !>          OpenMP can be activated through argument
  !>          YY can be initialized or not through argument
  !>
  !----------------------------------------------------------------------

  subroutine matrix_COO_SpMV(n1,n2,ndofr,ndofc,rows,cols,an,xx,yy,OPENMP,INITIALIZATION,CHUNK,SCHEDULE)
    !
    ! Dummy arguments
    !
    integer(ip),  intent(in)                    :: n1                !< Starting node
    integer(ip),  intent(in)                    :: n2                !< Final node
    integer(ip),  intent(in)                    :: ndofr             !< Number of dof per node (cols)
    integer(ip),  intent(in)                    :: ndofc             !< Number of dof per node (rows)
    integer(ip),  intent(in), pointer           :: rows(:)           !< Matrix graph ia
    integer(ip),  intent(in), pointer           :: cols(:)           !< Matrix graph ja
    real(rp),     intent(in)                    :: an(ndofc,ndofr,*) !< Matrix
    real(rp),     intent(in)                    :: xx(ndofc,*)       !< Input vector
    real(rp),     intent(out)                   :: yy(ndofr,*)       !< Output vector
    logical(lg),  intent(in),          optional :: OPENMP            !< If OpenMP should be used
    logical(lg),  intent(in),          optional :: INITIALIZATION    !< If array should be initialized
    integer(ip),  intent(in),          optional :: CHUNK             !< Chunks for dynamic scheduling of OpenMP
    character(*), intent(in),          optional :: SCHEDULE          !< OpenMP schedule
    integer(ip)                                 :: ii,iz
    integer(ip)                                 :: col,kk,ll
    real(rp)                                    :: raux
    logical(lg)                                 :: use_openmp
    logical(lg)                                 :: do_initialize
    integer(ip)                                 :: my_schedule
    integer(ip)                                 :: my_chunk,nz
#if defined(ALYA_OMPSS) || defined(_OPENMP)
    integer(ip)                                 :: jj
#endif

    nz = size(rows,KIND=ip)
    !
    ! Check if OpenMP should be used
    !
    use_openmp    = .false.
    do_initialize = .true.
    
    if( present(OPENMP)         ) use_openmp    = OPENMP
    if( present(INITIALIZATION) ) do_initialize = INITIALIZATION
    !
    ! OpenMP options
    !    
    if( use_openmp ) then
       
       my_chunk    = spmv_chunk
       my_schedule = SOL_OMP_STATIC
       
       if( present(CHUNK) ) my_chunk = CHUNK
       if( present(SCHEDULE) ) then
          if(      trim(SCHEDULE) == 'STATIC' ) then
             my_schedule = SOL_OMP_STATIC
          else if( trim(SCHEDULE) == 'DYNAMIC' ) then
             if( my_chunk <= 1 ) then
                my_schedule = SOL_OMP_STATIC
             else
                my_schedule = SOL_OMP_DYNAMIC
             end if
          else if( trim(SCHEDULE) == 'GUIDED' ) then
             my_schedule = SOL_OMP_GUIDED
          end if
       end if
       
    end if
    !
    ! Initialize solution
    !
    if( do_initialize ) then
       yy(1:ndofr,n1:n2) = 0.0_rp
    end if
    !
    ! NDOF = whatever
    !
    if( use_openmp ) then

       if( my_schedule == SOL_OMP_STATIC ) then
          !$OMP PARALLEL   DO                                                 &
          !$OMP SCHEDULE ( STATIC )                                           &
          !$OMP DEFAULT  ( NONE )                                             &
          !$OMP SHARED   ( an, rows, cols, n1, n2, ndofr, ndofc, xx, yy, nz ) &
          !$OMP PRIVATE  ( col, ii, jj, kk, ll, raux )
          do iz = 1,nz
             ii = rows(iz)
             if( ii >= n1 .and. ii <= n2 ) then
                col = cols(iz)
                do ll = 1,ndofc
                   raux = xx(ll,col)
                   do kk = 1,ndofr
                      !$OMP ATOMIC
                      yy(kk,ii) = yy(kk,ii) + an(ll,kk,iz) * raux
                   end do
                end do
             end if
          end do
          !$OMP END PARALLEL DO

       else if( my_schedule == SOL_OMP_GUIDED ) then
          !$OMP PARALLEL   DO                                                 &
          !$OMP SCHEDULE ( GUIDED )                                           &
          !$OMP DEFAULT  ( NONE )                                             &
          !$OMP SHARED   ( an, rows, cols, n1, n2, ndofr, ndofc, xx, yy, nz ) &
          !$OMP PRIVATE  ( col, ii, jj, kk, ll, raux )
          do iz = 1,nz
             ii = rows(iz)
             if( ii >= n1 .and. ii <= n2 ) then
                col = cols(iz)
                do ll = 1,ndofc
                   raux = xx(ll,col)
                   do kk = 1,ndofr
                      !$OMP ATOMIC
                      yy(kk,ii) = yy(kk,ii) + an(ll,kk,iz) * raux
                   end do
                end do
             end if
          end do
          !$OMP END PARALLEL DO

       else if( my_schedule == SOL_OMP_DYNAMIC ) then
          !$OMP PARALLEL   DO                                                 &
          !$OMP SCHEDULE ( DYNAMIC , my_chunk )                               &
          !$OMP DEFAULT  ( NONE )                                             &
          !$OMP SHARED   ( an, rows, cols, n1, n2, ndofr, ndofc, xx, yy, nz ) &
          !$OMP PRIVATE  ( col, ii, jj, kk, ll, raux )
          do iz = 1,nz
             ii = rows(iz)
             if( ii >= n1 .and. ii <= n2 ) then
                col = cols(iz)
                do ll = 1,ndofc
                   raux = xx(ll,col)
                   do kk = 1,ndofr
                      !$OMP ATOMIC
                      yy(kk,ii) = yy(kk,ii) + an(ll,kk,iz) * raux
                   end do
                end do
             end if
          end do
          !$OMP END PARALLEL DO
       end if

    else
       do iz = 1,nz
          ii = rows(iz)
          if( ii >= n1 .and. ii <= n2 ) then
             col = cols(iz)
             do ll = 1,ndofc
                raux = xx(ll,col)
                yy(1:ndofr,ii) = yy(1:ndofr,ii) + an(ll,1:ndofr,iz) * raux
             end do
          end if
       end do

    end if

  end subroutine matrix_COO_SpMV

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-12-04
  !> @brief   Scaling of a matrix
  !> @details Scaling from left and/or right of a matrix in CSR format
  !>
  !>          Unity test for matrix scaling:
  !>         
  !>                            +-          -+
  !>                            |  1   2   3 |
  !>          Original matrix = |  4   5   6 |
  !>                            |  7   8   9 |
  !>                            +-          -+
  !>         
  !>                            +-          -+
  !>                            |  2   0   0 |
  !>          Scaling =         |  0   4   0 |
  !>                            |  0   0   8 |
  !>                            +-          -+
  !>         
  !>                            +-          -+
  !>                            |  2   4   6 |
  !>          Left scaling =    | 16  20  24 |
  !>                            | 56  64  72 |
  !>                            +-          -+
  !>         
  !>                            +-          -+
  !>                            |  2   8  24 |
  !>          Right scaling =   |  8  20  48 |
  !>                            | 14  32  72 |
  !>                            +-          -+
  !>         
  !>         integer(ip)          :: ii,nn,nz,iboun
  !>         integer(ip), pointer :: ia(:)
  !>         integer(ip), pointer :: ja(:)
  !>         real(rp),    pointer :: aa(:)
  !>         real(rp),    pointer :: diag(:)
  !>                    
  !>         nn = 3
  !>         nz = 9
  !>         allocate(ia(nn+1),ja(nz))
  !>         allocate(aa(nz))
  !>         allocate(diag(nn))
  !>         ia   = (/ 1,4,7,10 /)
  !>         ja   = (/ 1,2,3,1,2,3,1,2,3 /)
  !>         aa   = (/ 1.0_rp,2.0_rp,3.0_rp,4.0_rp,5.0_rp,6.0_rp,7.0_rp,8.0_rp,9.0_rp /)
  !>         diag = (/ 2.0_rp,4.0_rp,8.0_rp /)
  !>         call matrix_scaling_CSR(1_ip,nn,1_ip,ia,ja,aa,diag,LEFT_SCALING=.true.,RIGHT_SCALING=.false.)
  !>         do ii = 1,nn
  !>            write(*,*) aa(ia(ii):ia(ii+1)-1)
  !>         end do
  !> 
  !-----------------------------------------------------------------------

  subroutine matrix_scaling_CSR(n1,n2,ndofn,ia,ja,an,diag,LEFT_SCALING,RIGHT_SCALING)

    integer(ip),  intent(in)                    :: n1                !< Starting node
    integer(ip),  intent(in)                    :: n2                !< Final node
    integer(ip),  intent(in)                    :: ndofn             !< Number of dof per node
    integer(ip),  intent(in), pointer           :: ia(:)             !< Matrix graph ia
    integer(ip),  intent(in), pointer           :: ja(:)             !< Matrix graph ja
    real(rp),     intent(inout)                 :: an(ndofn,ndofn,*) !< Matrix
    real(rp),     intent(in)                    :: diag(*)           !< Scaling diagonal matrix
    logical(lg),  intent(in),          optional :: LEFT_SCALING      !< Left scaling
    logical(lg),  intent(in),          optional :: RIGHT_SCALING     !< Right scaling
    integer(ip)                                 :: idofc
    integer(ip)                                 :: idofr
    integer(ip)                                 :: ii,jj,iz
    logical(lg)                                 :: if_left_scaling
    logical(lg)                                 :: if_right_scaling
    
    if_left_scaling  = .true.
    if_right_scaling = .false.
    if( present(LEFT_SCALING) )  if_left_scaling  = LEFT_SCALING
    if( present(RIGHT_SCALING) ) if_right_scaling = RIGHT_SCALING

    if( if_left_scaling ) then
       
       do ii = n1,n2
          do iz = ia(ii),ia(ii+1)-1
             jj = ja(iz)
             do idofr = 1,ndofn
                do idofc = 1,ndofn
                   an(idofc,idofr,iz) = an(idofc,idofr,iz) * diag((ii-1)*ndofn+idofr)
                end do
             end do
          end do
       end do

    end if

    if( if_right_scaling ) then
       
       do ii = n1,n2
          do iz = ia(ii),ia(ii+1)-1
             jj = ja(iz)
             do idofr = 1,ndofn
                do idofc = 1,ndofn
                   print*,jj,(jj-1)*ndofn+idofc
                   an(idofc,idofr,iz) = an(idofc,idofr,iz) * diag((jj-1)*ndofn+idofc)
                end do
             end do
          end do
       end do

    end if

  end subroutine matrix_scaling_CSR

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-20
  !> @brief   Rotate a matrix
  !> @details Rotate a matrix in CSR format
  !> 
  !-----------------------------------------------------------------------
  
  subroutine matrix_CSR_rotate(nn,ndofr,ndofc,iA,jA,kfl_fixrs,permr,skews,aa)

    integer(ip),          intent(in)    :: nn
    integer(ip),          intent(in)    :: ndofr
    integer(ip),          intent(in)    :: ndofc
    integer(ip), pointer, intent(in)    :: iA(:)
    integer(ip), pointer, intent(in)    :: jA(:)
    integer(ip), pointer, intent(in)    :: kfl_fixrs(:)
    integer(ip), pointer, intent(in)    :: permr(:)
    real(rp),    pointer, intent(in)    :: skews(:,:,:)
    real(rp),             intent(inout) :: aa(ndofc,ndofr,*)   
    integer(ip)                         :: ii,jj,kk,idime,jdime
    integer(ip)                         :: ibopo,kdime,iz,jz
    real(rp),    pointer                :: rotma(:,:)
    real(rp)                            :: worma(3,3)
    logical(lg)                         :: if_row
    logical(lg)                         :: if_col
 
    if( associated(skews) ) then

       kdime = size(skews,DIM=1,KIND=ip)
       if( kdime == ndofc ) then
          if_col = .true.
       else
          if_col = .false.
       end if
       if( kdime == ndofr ) then
          if_row = .true.
       else
          if_row = .false.
       end if

       do ii = 1,nn
          if( kfl_fixrs(ii) /= 0 ) then
             ibopo =  permr(ii)
             rotma => skews(:,:,ibopo)
             !
             ! Modifies column number II of AA ( AA_j,imodi <-- AA_j,imodi R )
             !
             if( if_col ) then
                do jz = ia(ii),ia(ii+1)-1
                   jj = jA(jz)
                   do iz = iA(jj),iA(jj+1)-1
                      kk = jA(iz)
                      if( kk == ii ) then                   
                         do idime = 1,ndofr
                            do jdime = 1,ndofc
                               worma(idime,jdime) = 0.0_rp
                               do kdime = 1,ndofc
                                  worma(idime,jdime) = worma(idime,jdime) + aa(kdime,idime,iz) * rotma(kdime,jdime)
                               end do
                            end do
                         end do
                         do idime = 1,ndofr
                            do jdime = 1,ndofc
                               aa(jdime,idime,iz) = worma(idime,jdime)
                            end do
                         end do
                      end if
                   end do
                end do
             end if
             !
             ! Modifies row number II of AA ( AA_imodi,j <-- R^t AA_imodi,j )
             !
             if( if_row ) then
                do iz = iA(ii),iA(ii+1)-1
                   jj = jA(iz)
                   do idime = 1,ndofr
                      do jdime = 1,ndofc
                         worma(idime,jdime) = 0.0_rp
                         do kdime = 1,ndofr
                            worma(idime,jdime) = worma(idime,jdime) &
                                 + aa(jdime,kdime,iz) * rotma(kdime,idime)
                         end do
                      end do
                   end do
                   do idime = 1,ndofr
                      do jdime = 1,ndofc
                         aa(jdime,idime,iz) = worma(idime,jdime)
                      end do
                   end do
                end do
             end if
          end if
       end do
       
    end if

  end subroutine matrix_CSR_rotate
     
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-04-17
  !> @brief   Rotate a matrix and a RHS
  !> @details Rotate a matrix AA and a RHS, BB
  !> 
  !-----------------------------------------------------------------------

  subroutine matrix_rotate_system(nn,ndofn,irow,iA,jA,rotma,bb,aa)

    integer(ip),           intent(in)    :: nn
    integer(ip),           intent(in)    :: ndofn
    integer(ip),           intent(in)    :: irow
    integer(ip),           intent(in)    :: ia(*)
    integer(ip),           intent(in)    :: ja(*)
    real(rp),              intent(in)    :: rotma(ndofn,ndofn)
    real(rp),    optional, intent(inout) :: bb(ndofn,*)
    real(rp),    optional, intent(inout) :: aa(ndofn,ndofn,*)
    integer(ip)                          :: jj,kk,idofn,jdofn,iz
    real(rp)                             :: worma(ndofn,ndofn)
    !
    ! Modifies column number II of AA ( A_j,imodi <-- A_j,imodi R )
    !
    do jj = 1,nn
       do iz = iA(jj),iA(jj+1)-1
          kk = jA(iz)
          if( kk == irow ) then
             if( present(aa) ) then
                worma = 0.0_rp
                do jdofn = 1,ndofn
                   do idofn = 1,ndofn
                      do kk = 1,ndofn
                         worma(idofn,jdofn) = worma(idofn,jdofn) &
                              + aa(kk,idofn,iz) * rotma(kk,jdofn)
                      end do
                   end do
                end do
                do idofn = 1,ndofn
                   do jdofn = 1,ndofn
                      aa(jdofn,idofn,iz) = worma(idofn,jdofn)
                   end do
                end do
             end if             
             
          end if
       end do
    end do
    !
    ! Modifies row number II of AA ( A_imodi,j <-- R^t A_imodi,j )
    !    
    do iz = iA(irow),iA(irow+1)-1
       if( present(aa) ) then
          jj = jA(iz)
          worma = 0.0_rp
          do jdofn = 1,ndofn
             do idofn = 1,ndofn
                do kk = 1,ndofn
                   worma(idofn,jdofn) = worma(idofn,jdofn) &
                        + aa(jdofn,kk,iz) * rotma(kk,idofn)
                end do
             end do
          end do
          do idofn = 1,ndofn
             do jdofn = 1,ndofn
                aa(jdofn,idofn,iz) = worma(idofn,jdofn)
             end do
          end do
       end if
       
    end do
    !
    ! Modify RHS BB
    !
    if( present(bb) ) then
       do idofn = 1,ndofn
          worma(idofn,1) = 0.0_rp
          do kk = 1,ndofn
             worma(idofn,1) = worma(idofn,1) &
                  + rotma(kk,idofn) * bb(kk,irow)
          end do
       end do
       do idofn = 1,ndofn
          bb(idofn,irow) = worma(idofn,1)
       end do
    end if
    
  end subroutine matrix_rotate_system

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    01/03/2016
  !> @brief   Add diagonal
  !> @details Add diagonal to a CSR matrix
  !>          A <= alpha * A + beta * D
  !
  !-----------------------------------------------------------------------

  subroutine matrix_add_diagonal_CSR(nbnodes,ndofn,ia,ja,diag,aa,MATRIX_FACTOR,DIAG_FACTOR)

    integer(ip),           intent(in)    :: nbnodes
    integer(ip),           intent(in)    :: ndofn
    integer(ip),           intent(in)    :: ia(*)
    integer(ip),           intent(in)    :: ja(*)
    real(rp),              intent(in)    :: diag(ndofn,*)
    real(rp),              intent(inout) :: aa(ndofn,ndofn,*)
    real(rp),    optional, intent(in)    :: MATRIX_FACTOR
    real(rp),    optional, intent(in)    :: DIAG_FACTOR
    integer(ip)                          :: ii,jj,iz,kk
    real(rp)                             :: xfact_matrix
    real(rp)                             :: xfact_diag

    xfact_matrix = optional_argument(1.0_rp,MATRIX_FACTOR)
    xfact_diag   = optional_argument(1.0_rp,DIAG_FACTOR)

    do ii = 1,nbnodes
       loop_iz: do iz = ia(ii),ia(ii+1)-1
          jj = ja(iz)
          if( jj == ii ) then
             do kk = 1,ndofn
                aa(kk,kk,iz) = xfact_matrix * aa(kk,kk,iz) + xfact_diag * diag(kk,ii)
             end do
             exit loop_iz
          end if
       end do loop_iz
    end do    

  end subroutine matrix_add_diagonal_CSR

  ! here I asuming same sparcesibity
  subroutine matrix_add_CSR(nbnodes,ndofn,ia,ja,add,aa,MATRIX_FACTOR,ADD_FACTOR)

    integer(ip),           intent(in)    :: nbnodes
    integer(ip),           intent(in)    :: ndofn
    integer(ip),           intent(in)    :: ia(*)
    integer(ip),           intent(in)    :: ja(*)
    real(rp),              intent(in)    :: add(ndofn,ndofn,*)
    real(rp),              intent(inout) :: aa(ndofn,ndofn,*)
    real(rp),    optional, intent(in)    :: MATRIX_FACTOR
    real(rp),    optional, intent(in)    :: ADD_FACTOR
    integer(ip)                          :: ii,iz,idofr,idofc
    real(rp)                             :: xfact_matrix
    real(rp)                             :: xfact_add

    xfact_matrix = optional_argument(1.0_rp,MATRIX_FACTOR)
    xfact_add    = optional_argument(1.0_rp,ADD_FACTOR)

    do ii = 1,nbnodes
       do iz = ia(ii),ia(ii+1)-1
          do idofr = 1,ndofn
             do idofc = 1,ndofn
                aa(idofc,idofr,iz) = xfact_matrix * aa(idofc,idofr,iz) + xfact_add * add(idofc,idofr,iz)
             end do
          end do
       end do 
    end do    

  end subroutine matrix_add_CSR

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-28
  !> @brief   Cast a LHS to RHS
  !> @details Send a LHS contribution to a LHS. This is usfeul
  !>          for explicit schemes
  !> 
  !-----------------------------------------------------------------------

  pure subroutine matrix_LHS_TO_RHS(nn,uu,elmat,elrhs)
    
    integer(ip), intent(in)    :: nn
    real(rp),    intent(in)    :: uu(nn)
    real(rp),    intent(inout) :: elmat(nn,nn)
    real(rp),    intent(inout) :: elrhs(nn)

    elrhs = elrhs - matmul(elmat,uu)
    elmat = 0.0_rp
    
  end subroutine matrix_LHS_TO_RHS
  
end module mod_matrix
!> @}


