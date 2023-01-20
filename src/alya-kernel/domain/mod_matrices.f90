!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Matrix_Toolbox
!> @{
!> @name    ToolBox for matrix operations
!> @file    mod_matrix.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!> @brief   ToolBox for matrix operations
!> @details ToolBox for matrix operations: fill in, etc.
!------------------------------------------------------------------------

module mod_matrices
#include "def_vector_size.inc"
  use def_kintyp,              only : ip,rp,lg,i1p,r1p
  use mod_memory,              only : memory_alloca,memory_deallo
  use def_master,              only : INOTMASTER,kfl_paral
  use def_kermod,              only : kfl_element_to_csr
  use def_domain,              only : elmar
  use def_domain,              only : mnode
  use def_domain,              only : mgaus
  use def_domain,              only : ndime
  use def_domain,              only : nelem
  use def_domain,              only : nzdom
  use def_domain,              only : r_dom
  use def_domain,              only : c_dom
  use def_domain,              only : ltype
  use def_domain,              only : lgaus
  use def_domain,              only : lezdo
  use def_domain,              only : ntens
  use def_domain,              only : lnods
  use def_domain,              only : lnnod
  use def_domain,              only : coord
  use def_domain,              only : ngaus
  use def_domain,              only : kfl_savda
  use def_domain,              only : elmda_gpvol
  use def_domain,              only : elmda_gpcar
  use def_domain,              only : memor_dom
  use def_elmgeo,              only : element_type
  use mod_element_integration, only : element_shape_function_derivatives_jacobian
  use mod_optional_argument,   only : optional_argument
  use mod_parall,              only : num_subd_par
  use mod_parall,              only : num_pack_par
  use mod_parall,              only : list_elements_par
  use mod_elmgeo_vector,       only : elmgeo_cartesian_derivatives_jacobian
  implicit none

  private

  public :: matrices_gradient_divergence
  public :: matrices_laplacian
  public :: matrices_all

contains

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    02/05/2017
  !> @brief   Compute the gradient matrix
  !> @details Compute the gradient matrix
  !>
  !----------------------------------------------------------------------

  subroutine matrices_gradient_divergence(Grad,Div,kdime,DEALLOCATE_MATRICES) 

    real(rp),              pointer, intent(inout) :: Grad(:,:)
    real(rp),    optional, pointer, intent(inout) :: Div(:,:)
    integer(ip), optional,          intent(in)    :: kdime
    logical(lg), optional,          intent(in)    :: DEALLOCATE_MATRICES
    integer(ip)                                   :: ielem
    integer(ip)                                   :: pnode
    integer(ip)                                   :: plapl
    integer(ip)                                   :: pgaus
    integer(ip)                                   :: pelty
    integer(ip)                                   :: inode

    real(rp)                                      :: elcod(ndime,mnode)
    real(rp)                                      :: gpvol(mgaus)
    real(rp)                                      :: gpsha(mnode,mgaus)
    real(rp)                                      :: gpder(ndime,mnode,mgaus)
    real(rp)                                      :: gpcar(ndime,mnode,mgaus)
    real(rp)                                      :: gphes(ntens,mnode,mgaus)
    logical(lg)                                   :: if_deallocate_matrices

    
    if( INOTMASTER ) then

       if( present(DEALLOCATE_MATRICES) ) then
          if_deallocate_matrices = DEALLOCATE_MATRICES
       else
          if_deallocate_matrices = .false.
       end if

       if( if_deallocate_matrices ) then
          call memory_deallo(memor_dom,'GRAD','matrices_gradient',Grad)
          if( present(Div) ) then
             call memory_deallo(memor_dom,'Div','matrices_gradient',Div)
          end if
       end if 
       
       if( .not. associated(Grad) ) then
          call memory_alloca(memor_dom,'GRAD','matrices_gradient',Grad,ndime,nzdom)
       else
          Grad = 0.0_rp
       end if
       if( present(Div) ) then
          if( .not. associated(Div) ) then
             call memory_alloca(memor_dom,'Div','matrices_gradient',Div,ndime,nzdom)
          else
             Div = 0.0_rp
          end if
       end if
       
       do ielem = 1,nelem
          pelty = ltype(ielem)
          if( pelty > 0 ) then
             pgaus = lgaus(ielem)
             pnode = lnnod(ielem)
             plapl = 0
             do inode = 1,pnode
                elcod(1:ndime,inode) = coord(1:ndime,lnods(inode,ielem))
             end do
             call element_shape_function_derivatives_jacobian(&
                  pnode,pgaus,plapl,elmar(pelty) % weigp,elmar(pelty) % shape,&
                  elmar(pelty) % deriv,elmar(pelty) % heslo,&
                  elcod,gpvol,gpsha,gpder,gpcar,gphes,ielem)
             if( present(Div) ) then
                call matrices_element_gradient(&
                     pgaus,pnode,ielem,lnods(:,ielem),gpvol,gpsha,gpcar,&
                     Grad,Div)
             else
                call matrices_element_gradient(&
                     pgaus,pnode,ielem,lnods(:,ielem),gpvol,gpsha,gpcar,&
                     Grad)
             end if
          end if
       end do


    else

       !
       ! Allocate of size 1 otherwise I got : Attempt to use pointer GRAD when it is not associated with a target
       !
       if( .not. associated(Grad) ) then
          call memory_alloca(memor_dom,'GRAD','matrices_gradient',Grad,ndime,1_ip)
       else
          Grad = 0.0_rp
       end if
       if( present(Div) ) then
          if( .not. associated(Div) ) then
             call memory_alloca(memor_dom,'Div','matrices_gradient',Div,ndime,1_ip)
          else
             Div = 0.0_rp
          end if
       end if
       
    end if
    
  end subroutine matrices_gradient_divergence

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    02/05/2017
  !> @brief   Element and global gradient matrices
  !> @details Element and global gradient matrices
  !>
  !----------------------------------------------------------------------

  subroutine matrices_element_gradient(&
       pgaus,pnode,ielem,lnods,gpvol,gpsha,gpcar,Grad,Div)

    integer(ip),                    intent(in)    :: pgaus
    integer(ip),                    intent(in)    :: pnode
    integer(ip),                    intent(in)    :: ielem
    integer(ip),                    intent(in)    :: lnods(pnode)
    real(rp),                       intent(in)    :: gpvol(pgaus)
    real(rp),                       intent(in)    :: gpsha(pnode,pgaus)
    real(rp),                       intent(in)    :: gpcar(ndime,mnode,pgaus)
    real(rp),    pointer,           intent(inout) :: Grad(:,:)
    real(rp),    pointer, optional, intent(inout) :: Div(:,:)
    real(rp)                                      :: elgra(ndime*pnode,pnode) ! Aup
    real(rp)                                      :: eldiv(pnode,ndime*pnode) ! Apu
    integer(ip)                                   :: igaus,inode,idime,idofn,iz
    integer(ip)                                   :: jcolu,jnode,ipoin,jpoin
    integer(ip)                                   :: jdofn
    !
    ! Element matrix
    !
    elgra = 0.0_rp
    eldiv = 0.0_rp
    do igaus = 1,pgaus
       do inode = 1,pnode
          do idime = 1,ndime
             idofn = (inode-1)*ndime+idime
             do jnode = 1,pnode
                elgra(idofn,jnode) = elgra(idofn,jnode) &
                     - gpvol(igaus) * gpsha(jnode,igaus) * gpcar(idime,inode,igaus)
                eldiv(jnode,idofn) = eldiv(jnode,idofn) &
                     + gpvol(igaus) * gpsha(jnode,igaus) * gpcar(idime,inode,igaus)
             end do
          end do
       end do
    end do

    !do inode = 1,pnode*ndime
    !   write(*,'(10(1x,e13.6))') elgra(idime,inode,1:pnode)
    !end do
    !call runend('O.K.!')
    !
    ! Assemble global system
    !
    if( kfl_element_to_csr == 1 ) then
      
       do inode = 1,pnode
          do jnode = 1,pnode
             iz = lezdo(inode,jnode,ielem)
             do idime = 1,ndime
                idofn = (inode-1) * ndime + idime
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
                Grad(idime,iz) = Grad(idime,iz) + elgra(idofn,jnode) ! G
             end do
             if( present(Div) ) then
                do idime = 1,ndime
                   jdofn = (jnode-1) * ndime + idime
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
                   Div(idime,iz) = Div(idime,iz) + eldiv(inode,jdofn) ! D                      
                end do
             end if
          end do
       end do
       
    else

          do inode = 1,pnode
             ipoin = lnods(inode)
             do jnode = 1,pnode
                jpoin = lnods(jnode)
                iz    = r_dom(ipoin)
                jcolu = c_dom(iz)
                do while( jcolu /= jpoin .and. iz < r_dom(ipoin+1)-1 )
                   iz    = iz + 1
                   jcolu = c_dom(iz)
                end do
                if( jcolu == jpoin ) then
                   do idime = 1,ndime
                      idofn = (inode-1) * ndime + idime
#ifdef NO_COLORING
                      !$OMP ATOMIC
#endif
                      Grad(idime,iz) = Grad(idime,iz) + elgra(idofn,jnode)  ! G
                   end do
                   if( present(Div) ) then
                      do idime = 1,ndime
                         jdofn = (jnode-1) * ndime + idime
#ifdef NO_COLORING
                         !$OMP ATOMIC
#endif
                         Div(idime,iz) = Div(idime,iz) + eldiv(inode,jdofn) ! D
                      end do
                   end if
                end if
             end do
          end do
       
    end if

  end subroutine matrices_element_gradient


  !----------------------------------------------------------------------
  !>
  !> @author  Herbert Owen
  !> @date    27/07/2017
  !> @brief   Compute the laplacian matrix
  !> @details Compute the laplacian matrix
  !>
  !----------------------------------------------------------------------

  subroutine matrices_laplacian(Lapl,DEALLOCATE_MATRIX,NAME) 

    real(rp),          pointer, intent(inout) :: Lapl(:)
    logical(lg),      optional, intent(in)    :: DEALLOCATE_MATRIX
    character(LEN=*), optional, intent(in)    :: NAME
    integer(ip)                               :: ielem
    integer(ip)                               :: pnode
    integer(ip)                               :: plapl
    integer(ip)                               :: pgaus
    integer(ip)                               :: pelty
    integer(ip)                               :: inode

    logical(lg)                               :: if_deallocate_matrices
    integer(ip)                               :: ipack,isubd
    character(LEN=:), allocatable             :: my_name
    
    real(rp)                                  :: elcod(VECTOR_SIZE,ndime,mnode)
    real(rp)                                  :: gpvol(VECTOR_SIZE,mgaus)
    real(rp)                                  :: gpsha(VECTOR_SIZE,mnode,mgaus)
    real(rp)                                  :: gpcar(VECTOR_SIZE,ndime,mnode,mgaus)

    my_name                = optional_argument('LAPL',NAME)
    if_deallocate_matrices = optional_argument(.false.,DEALLOCATE_MATRIX)
    
    if( INOTMASTER ) then

       if( if_deallocate_matrices ) then
          call memory_deallo(memor_dom,trim(my_name),'matrices_gradient',Lapl)
       end if
       
       if( .not. associated(Lapl) ) then
          call memory_alloca(memor_dom,trim(my_name),'matrices_gradient',Lapl,nzdom)
       else
          Lapl = 0.0_rp
       end if

       do isubd = 1,num_subd_par 
          do ipack = 1,num_pack_par(isubd)
             
             ielem = list_elements_par(isubd) % packs(ipack) % l(1)  ! Select first element
             pelty = abs(ltype(ielem))                               ! Element type
             pnode = element_type(pelty) % number_nodes              ! Number of nodes
             pgaus = ngaus(pelty)                                    ! Number of Gauss points
            
             call matrices_element_laplacian(&
                  size(list_elements_par(isubd) % packs(ipack) % l,kind=ip),&
                  pnode,pgaus,pelty,list_elements_par(isubd) % packs(ipack) % l,&
                  elcod,gpvol,gpsha,gpcar,Lapl)
          end do
          
       end do

    else

       !
       ! Allocate of size 1 otherwise I got : Attempt to use pointer L when it is not associated with a target
       !
       if( .not. associated(Lapl) ) then
          call memory_alloca(memor_dom,'LAPL','matrices_gradient',Lapl,1_ip)
       else
          Lapl = 0.0_rp
       end if
       
    end if
    
  end subroutine matrices_laplacian

  !----------------------------------------------------------------------
  !>
  !> @author  Herbert Owen
  !> @date    27/07/2017
  !> @brief   Compute the laplacian matrix
  !> @details Compute the laplacian matrix
  !>
  !----------------------------------------------------------------------

  subroutine matrices_all(Lapl,Grad,Div,DEALLOCATE_MATRIX,NAME) 

    real(rp),         optional, pointer, intent(inout) :: Lapl(:)
    real(rp),         optional, pointer, intent(inout) :: Grad(:,:)
    real(rp),         optional, pointer, intent(inout) :: Div(:,:)
    logical(lg),      optional,          intent(in)    :: DEALLOCATE_MATRIX
    character(LEN=*), optional,          intent(in)    :: NAME
    integer(ip)                                        :: ielem
    integer(ip)                                        :: pnode
    integer(ip)                                        :: plapl
    integer(ip)                                        :: pgaus
    integer(ip)                                        :: pelty
    integer(ip)                                        :: inode
    logical(lg)                                        :: if_deallocate_matrices
    integer(ip)                                        :: ipack,isubd
    character(LEN=:), allocatable                      :: my_name   
    real(rp)                                           :: elcod(VECTOR_SIZE,ndime,mnode)
    real(rp)                                           :: gpvol(VECTOR_SIZE,mgaus)
    real(rp)                                           :: gpsha(VECTOR_SIZE,mnode,mgaus)
    real(rp)                                           :: gpcar(VECTOR_SIZE,ndime,mnode,mgaus)

    my_name                = optional_argument('LAPL-GRAD-DIV',NAME)
    if_deallocate_matrices = optional_argument(.false.,DEALLOCATE_MATRIX)
    
    if( INOTMASTER ) then
       !
       ! Allocate
       !
       if( present(Lapl) ) then
          if( if_deallocate_matrices ) then
             call memory_deallo(memor_dom,'LAPL','matrices_gradient',Lapl)
          end if
          if( .not. associated(Lapl) ) then
             call memory_alloca(memor_dom,'LAPL','matrices_gradient',Lapl,nzdom)
          else
             Lapl = 0.0_rp
          end if
       end if
       if( present(Grad) ) then
          if( if_deallocate_matrices ) then
             call memory_deallo(memor_dom,'GRAD','matrices_gradient',Grad)
          end if
          if( .not. associated(Grad) ) then
             call memory_alloca(memor_dom,'GRAD','matrices_gradient',Grad,ndime,nzdom)
          else
             Grad = 0.0_rp
          end if
       end if
       if( present(Div) ) then
          if( if_deallocate_matrices ) then
             call memory_deallo(memor_dom,'GRAD','matrices_gradient',Div)
          end if
          if( .not. associated(Div) ) then
             call memory_alloca(memor_dom,'GRAD','matrices_gradient',Div,ndime,nzdom)
          else
             Div = 0.0_rp
          end if
       end if
       
       do isubd = 1,num_subd_par 
          do ipack = 1,num_pack_par(isubd)
             
             ielem = list_elements_par(isubd) % packs(ipack) % l(1)  ! Select first element
             pelty = abs(ltype(ielem))                               ! Element type
             pnode = element_type(pelty) % number_nodes              ! Number of nodes
             pgaus = ngaus(pelty)                                    ! Number of Gauss points
            
             call matrices_element_laplacian(&
                  size(list_elements_par(isubd) % packs(ipack) % l,kind=ip),&
                  pnode,pgaus,pelty,list_elements_par(isubd) % packs(ipack) % l,&
                  elcod,gpvol,gpsha,gpcar,Lapl,Grad,Div)
          end do
          
       end do

    else

       !
       ! Allocate of size 1 otherwise I got : Attempt to use pointer L when it is not associated with a target
       !
       if( present(Lapl) ) then
          if( .not. associated(Lapl) ) then
             call memory_alloca(memor_dom,'LAPL','matrices_gradient',Lapl,1_ip)
          else
             Lapl = 0.0_rp
          end if
       end if
       if( present(Grad) ) then
          if( .not. associated(Grad) ) then
             call memory_alloca(memor_dom,'GRAD','matrices_gradient',Grad,1_ip,1_ip)
          else
             Grad = 0.0_rp
          end if
       end if
       if( present(Div) ) then
          if( .not. associated(Div) ) then
             call memory_alloca(memor_dom,'GRAD','matrices_gradient',Div,1_ip,1_ip)
          else
             Div = 0.0_rp
          end if
       end if
       
    end if
    
  end subroutine matrices_all

  !----------------------------------------------------------------------
  !>
  !> @author  Herbert Owen
  !> @date    27/07/2017
  !> @brief   Element and global Laplacian matrix ( grad p , grad q )
  !> @details Element and global Laplacian matrix ( grad p , grad q )
  !>
  !----------------------------------------------------------------------
  
  subroutine matrices_element_laplacian(&
       VECTOR_DIM,pnode,pgaus,pelty,list_elements,elcod,gpvol,gpsha,gpcar,Lapl,Grad,Div)

    integer(ip),                    intent(in)    :: VECTOR_DIM                   !< Number of nodes
    integer(ip),                    intent(in)    :: pnode                        !< Number of nodes
    integer(ip),                    intent(in)    :: pgaus                        !< Number of Gauss points
    integer(ip),                    intent(in)    :: pelty                        !< Number of Gauss points
    integer(ip),                    intent(in)    :: list_elements(VECTOR_SIZE)    !< List of elements
    real(rp),                       intent(inout) :: elcod(VECTOR_SIZE,ndime,mnode)
    real(rp),                       intent(inout) :: gpvol(VECTOR_SIZE,mgaus)
    real(rp),                       intent(inout) :: gpsha(VECTOR_SIZE,mnode,mgaus)
    real(rp),                       intent(inout) :: gpcar(VECTOR_SIZE,ndime,mnode,mgaus)
    real(rp),    optional, pointer, intent(inout) :: Lapl(:)                      !< Laplacian matrix
    real(rp),    optional, pointer, intent(inout) :: Grad(:,:)
    real(rp),    optional, pointer, intent(inout) :: Div(:,:)    
    real(rp)                                      :: ellap(VECTOR_SIZE,pnode,pnode) 
    real(rp)                                      :: elgra(VECTOR_SIZE,pnode*ndime,pnode) 
    real(rp)                                      :: xjaci(VECTOR_SIZE,ndime,ndime,pgaus)    
    integer(ip)                                   :: igaus,inode,idime,iz,ivect,idofn
    integer(ip)                                   :: jcolu,jnode,ipoin,jpoin,ielem
#define DEF_VECT 1:VECTOR_SIZE
    !
    ! Gather
    !
    do ivect = 1,VECTOR_SIZE
       ielem = abs(list_elements(ivect))
       if( ielem > 0 ) then
          do inode = 1,pnode
             elcod(ivect,1:ndime,inode) = coord(1:ndime,lnods(inode,ielem))
          end do
       else
          do inode = 1,pnode
             elcod(ivect,1:ndime,inode) = coord(1:ndime,lnods(inode,list_elements(1)))
          end do
       end if
    end do
    !
    ! Cartesian derivatives and Jacobian
    !
    if( kfl_savda == 2 ) then
       do ivect = 1,VECTOR_SIZE  
          ielem = abs(list_elements(ivect))
          if( ielem > 0 ) then
             gpvol(ivect,1:pgaus)                 = elmda_gpvol(1:pgaus,ielem)
             gpcar(ivect,1:ndime,1:mnode,1:pgaus) = elmda_gpcar(1:ndime,1:mnode,1:pgaus,ielem)
          else
             gpvol(ivect,1:pgaus)                 = 0.0_rp
             gpcar(ivect,1:ndime,1:mnode,1:pgaus) = 0.0_rp
          end if
       end do
    else
       call elmgeo_cartesian_derivatives_jacobian(VECTOR_DIM,ndime,mnode,pnode,pgaus,elcod,elmar(pelty) % deriv,xjaci,gpcar,gpvol)
       do igaus = 1,pgaus
          gpvol(DEF_VECT,igaus) = gpvol(DEF_VECT,igaus) * elmar(pelty) % weigp(igaus)
       end do
    end if    
    !
    ! Element matrix
    ! We use a + sign as in nsi_element_schur
    ! In JCP 2001 - Codina - Pressuer Stab..   it is defined with a - sign.
    !
    if( present(Lapl) )                   ellap = 0.0_rp
    if( present(Grad) .or. present(Div) ) elgra = 0.0_rp
    do igaus = 1,pgaus
       do inode = 1,pnode
          if( present(Lapl) ) then
             do idime = 1,ndime
                ellap(DEF_VECT,inode,inode) = ellap(DEF_VECT,inode,inode) &
                     + gpvol(DEF_VECT,igaus) * gpcar(DEF_VECT,idime,inode,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
             end do
             do jnode = inode+1,pnode
                do idime = 1,ndime
                   ellap(DEF_VECT,inode,jnode) = ellap(DEF_VECT,inode,jnode) &
                        + gpvol(DEF_VECT,igaus) * gpcar(DEF_VECT,idime,inode,igaus) * gpcar(DEF_VECT,idime,jnode,igaus)
                   ellap(DEF_VECT,jnode,inode) = ellap(DEF_VECT,inode,jnode) 
                end do
             end do
          end if
          if( present(Grad) .or. present(Div) ) then
             do idime = 1,ndime
                idofn = (inode-1)*ndime+idime
                do jnode = 1,pnode
                   elgra(DEF_VECT,idofn,jnode) = elgra(DEF_VECT,idofn,jnode) &
                        - gpvol(DEF_VECT,igaus) * elmar(pelty) % shape(jnode,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
                end do
             end do      
          end if
       end do
    end do
    !
    ! Assemble global system
    !
    if( kfl_element_to_csr == 1 ) then

       do ivect = 1,VECTOR_SIZE
          ielem = list_elements(ivect)
          if( ielem > 0 ) then
             do inode = 1,pnode
                do jnode = 1,pnode
                   iz = lezdo(inode,jnode,ielem)
                   if( present(Lapl) ) then
#ifdef NO_COLORING
                      !$OMP ATOMIC
#endif
                      Lapl(iz) = Lapl(iz) + ellap(ivect,inode,jnode)                ! L
                   end if
                   if( present(Grad) ) then
                      do idime = 1,ndime
                         idofn = (inode-1) * ndime + idime
#ifdef NO_COLORING
                      !$OMP ATOMIC
#endif
                         Grad(idime,iz) = Grad(idime,iz) + elgra(ivect,idofn,jnode) ! G
                      end do
                   end if
                   if( present(Div) ) then
                      do idime = 1,ndime
                         idofn = (jnode-1) * ndime + idime
#ifdef NO_COLORING
                      !$OMP ATOMIC
#endif
                         Div(idime,iz) = Div(idime,iz) - elgra(ivect,idofn,inode)   ! D
                      end do
                   end if
                end do
             end do
          end if
       end do
       
    else

       do ivect = 1,VECTOR_SIZE
          ielem = list_elements(ivect)
          if( ielem > 0 ) then
             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                do jnode = 1,pnode
                   jpoin = lnods(jnode,ielem)
                   iz    = r_dom(ipoin)
                   jcolu = c_dom(iz)
                   do while( jcolu /= jpoin .and. iz < r_dom(ipoin+1)-1 )
                      iz    = iz + 1
                      jcolu = c_dom(iz)
                   end do
                   if( jcolu == jpoin ) then
                      if( present(Lapl) ) then
#ifdef NO_COLORING
                         !$OMP ATOMIC
#endif
                         Lapl(iz) = Lapl(iz) + ellap(ivect,inode,jnode)                 ! L
                      end if
                      if( present(Grad) ) then
                         do idime = 1,ndime
                            idofn = (inode-1) * ndime + idime
#ifdef NO_COLORING
                            !$OMP ATOMIC
#endif
                            Grad(idime,iz) = Grad(idime,iz) + elgra(ivect,idofn,jnode)  ! G
                         end do
                      end if
                      if( present(Div) ) then
                         do idime = 1,ndime
                            idofn = (jnode-1) * ndime + idime
#ifdef NO_COLORING
                            !$OMP ATOMIC
#endif
                            Div(idime,iz) = Div(idime,iz) - elgra(ivect,idofn,inode)    ! D
                         end do
                      end if
                   end if
                end do
             end do
          end if
       end do
    end if

  end subroutine matrices_element_laplacian

  !
  ! Ojo hay una parecida en nsi_element_laplacian
  !  
  
end module mod_matrices
!> @}
