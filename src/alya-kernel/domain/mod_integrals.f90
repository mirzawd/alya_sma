!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    mod_integrals.f90
!> @author  houzeaux
!> @date    2020-02-03
!> @brief   Integrals
!> @details Compute generic integrals
!-----------------------------------------------------------------------

module mod_integrals

  use def_domain
  use def_master
  use mod_func
  use mod_communications, only : PAR_SUM
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_bouder,         only : bouder

  implicit none
  private
  
  type field_arrays
     real(rp), pointer :: a(:,:)
  end type field_arrays

  public :: integrals_volume
  public :: integrals_boundary
  public :: field_arrays
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux and eduardo perez
  !> @date    2020-02-03
  !> @brief   Volume integrals
  !> @details Volume integrals:
  !> This program evaluates a set of volumetric integrals I1, I2, ..., Im 
  !> over a group of clusters (or subdomains) in the form:
  !> 
  !> I1 = integral( g11(f1) * g12(f2) * ... * g1n(fn) * dV )
  !> I2 = integral( g21(f1) * g22(f2) * ... * g2n(fn) * dV )
  !> .........................
  !> Im = integral( gm1(f1) * gm2(f2) * ... * gmn(fn) * dV )
  !>
  !> where the fields are f1, f2, ..., fn contained in array_fields and over the fields a set
  !> of functions g11, g12, ..., g1n (for I1), g21, g22, ..., g2n (for I2), etc. are applied
  !> which are given by the matrix of pointers array_func.
  !>
  !> array_fields: matrix with the fields to be integrated. It has dimensions (npoin, n_field).
  !> array_func: array of pointers to the functions used to evaluate the integrands. It has
  !> dimensions (n_integrals, n_field).
  !> n_integrals: number of integrals to be performed.
  !> n_field: number of fields.
  !>
  !> Note that for any element, ielem, legro(ielem) is the cluster (or subdomain) number to which ielem belongs. 
  !>
  !> Example of use
  !>   implicit none
  !>   integer(ip), parameter :: nintegrals=3
  !>   integer(ip), parameter :: nfields=3
  !>   type(func_ptr)         :: f_ptr(nintegrals,3) 
  !>   integer(ip)            :: ipoin
  !>   integer(ip), pointer   :: legro(:)
  !>   type(field_arrays)     :: xx(nfields) 
  !>   real(rp)               :: results(nintegrals,1)
  !> 
  !>   nullify(legro)
  !>   call func_initialization(f_ptr)
  !> 
  !>   f_ptr(1,1) % f => func_div
  !>   f_ptr(1,2) % f => func_unity
  !>   f_ptr(1,3) % f => func_unity
  !>   
  !>   f_ptr(2,1) % f => func_unity
  !>   f_ptr(2,2) % f => func_identity
  !>   f_ptr(2,3) % f => func_identity
  !> 
  !>   f_ptr(3,1) % f => func_rotational_norm
  !>   f_ptr(3,2) % f => func_unity 
  !>   f_ptr(3,3) % f => func_unity 
  !> 
  !>   allocate(xx(1)%a(ndime,npoin))
  !>   allocate(xx(2)%a(1,npoin))
  !>   allocate(xx(3)%a(1,npoin))
  !> 
  !>   do ipoin = 1,npoin
  !>      xx(1) % a(1,ipoin) = 2.0_rp*coord(1,ipoin)+4.0_rp*coord(2,ipoin)
  !>      xx(1) % a(2,ipoin) = 3.0_rp*coord(2,ipoin)+7.0_rp*coord(1,ipoin)
  !>      xx(2) % a(1,ipoin) = 6.0_rp
  !>      xx(3) % a(1,ipoin) = 7.0_rp
  !>   end do
  !>   call integrals_volume(xx,f_ptr,legro,results)
  !>   if(inotslave) print*,results
  !>   call runend('O.K.!')
  !> 
  !-----------------------------------------------------------------------

  subroutine integrals_volume(array_fields,array_func,legro,fields_integrals)

    type(field_arrays),   intent(in)    :: array_fields(:)
    type(func_ptr),       intent(in)    :: array_func(:,:)
    integer(ip), pointer, intent(in)    :: legro(:)
    real(rp),             intent(inout) :: fields_integrals(:,:)
    integer(ip)                         :: ipoin,ielem,igaus,inode,jdime
    integer(ip)                         :: i_field,i_integrals
    integer(ip)                         :: n_field,n_integrals,ngrou
    integer(ip)                         :: pnode,pgaus,pelty,idime,kdime
    real(rp)                            :: gpvol,gpdet
    real(rp)                            :: elcod(ndime,mnode)
    real(rp)                            :: gpcar(ndime,mnode)
    real(rp)                            :: xjaci(ndime,ndime),xjacm(ndime,ndime)
    real(rp),             allocatable   :: gpfie(:)
    real(rp),             allocatable   :: gradf(:,:)
    real(rp),             allocatable   :: elfie(:,:,:) 
    real(rp),             allocatable   :: integrands(:)

    n_integrals = size(array_func,1)
    n_field     = size(array_func,2)
    ngrou       = size(fields_integrals,2)
    kdime       = 0
    
    if( size(array_fields)     < n_field )      call runend('MOD_INTEGRALS: WRONG DIMENSION')
    if( size(fields_integrals,1) < n_integrals) call runend('MOD_INTEGRALS: WRONG DIMENSION')
    
    allocate(integrands(n_integrals))

    if( INOTEMPTY ) then
       do i_field = 1,n_field      
          kdime = max(kdime,size(array_fields(i_field)%a,DIM=1,KIND=ip))
          if( size(array_fields(i_field)%a,DIM=2,KIND=ip) < npoin ) then
             call runend('MOD_INTEGRALS: WRONG FIELD DIMENSION')
          end if      
       end do      
       allocate(elfie(kdime,mnode,n_field))
       allocate(gpfie(kdime))
       allocate(gradf(ndime,kdime))
    end if

    fields_integrals(1:n_integrals,1:ngrou) = 0.0_rp
    !
    ! Compute integral over the whole domain
    !
    !$OMP  PARALLEL DO SCHEDULE (STATIC)                           &
    !$OMP  DEFAULT (NONE)                                          &
    !$OMP  FIRSTPRIVATE ( elfie,gpfie,gradf,integrands)            &
    !$OMP  PRIVATE ( ielem, pelty, pnode, inode, ipoin, igaus,     &
    !$OMP            gpcar, gpdet, xjacm, xjaci, gpvol, elcod,     &
    !$OMP            i_field, i_integrals,                         &
    !$OMP            idime, jdime, kdime , pgaus)                  &
    !$OMP  SHARED  ( ltype, nnode, lnods, coord, nelem, ngaus,     &
    !$OMP            elmar, n_integrals,  n_field, array_func,     &
    !$OMP            legro, ngrou,                                 &
#ifndef NDIMEPAR
    !$OMP            ndime,                                        &
#endif
    !$OMP            array_fields) &
    !$OMP REDUCTION ( +:fields_integrals)                      
    !
    elements: do ielem = 1, nelem
       if( legro(ielem) > ngrou ) call runend('MOD_INTEGRALS: WRONG GROUP NUMBER')
       if( legro(ielem) /= 0 ) then
          pelty = ltype(ielem)
          if (pelty > 0) then
             pnode = nnode(pelty)
             pgaus = ngaus(pelty)
             !
             ! Gather operations
             !
             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                do i_field = 1, n_field
                   kdime = size(array_fields(i_field) % a,1)
                   do idime = 1,kdime
                      elfie(idime,inode,i_field) = array_fields(i_field) % a(idime,ipoin)
                   end do
                end do
                elcod(1:ndime,inode) = coord(1:ndime,ipoin)
             end do
             !
             ! 1st order Cartesian derivatives, and dV:=GPVOL=|J|*wg
             !
             do igaus = 1, pgaus
                call elmder(&
                     pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&        ! Cartesian derivative
                     elcod,gpcar,gpdet,xjacm,xjaci)                     ! and Jacobian
                gpvol      = elmar(pelty)%weigp(igaus) * gpdet          ! |J|*wg
                integrands = 1.0_rp
   
                do i_field = 1, n_field
                   gpfie      = 0.0_rp
                   gradf      = 0.0_rp
                   kdime      = size(array_fields(i_field) % a,1)
   
                   do inode = 1, pnode
                      do idime = 1,kdime
                         gpfie(idime) = gpfie(idime) + elmar(pelty) % shape(inode,igaus) * elfie(idime,inode,i_field)
                         do jdime = 1,ndime
                            gradf(jdime,idime) = gradf(jdime,idime) + gpcar(jdime,inode) * elfie(idime,inode,i_field)
                         end do
                      end do
                   end do
                   
                   do i_integrals = 1,n_integrals
                      integrands(i_integrals) = integrands(i_integrals) * array_func(i_integrals,i_field) % f(gpfie(1:kdime),gradf)
                   end do
                                   
                end do
                
                do i_integrals = 1,n_integrals
                   fields_integrals(i_integrals,legro(ielem)) = fields_integrals(i_integrals,legro(ielem)) + gpvol * integrands(i_integrals)
                end do
   
             end do
          end if
       end if
    end do elements
    !$OMP END PARALLEL DO

    call PAR_SUM(n_integrals,ngrou,fields_integrals)

    if(allocated(integrands)) deallocate(integrands)
    if(allocated(elfie))      deallocate(elfie)
    if(allocated(gpfie))      deallocate(gpfie)
    if(allocated(gradf))      deallocate(gradf)

    !do i_field = 1, n_field
    !   if(inotslave) print*,i_field,fields_integrals(i_field)
    !end do
    
  end subroutine integrals_volume
  
  subroutine integrals_boundary(array_fields,array_func,legro,fields_integrals)

    type(field_arrays),   intent(in)    :: array_fields(:)
    type(func_ptr),       intent(in)    :: array_func(:,:)
    integer(ip), pointer, intent(in)    :: legro(:)
    real(rp),             intent(inout) :: fields_integrals(:,:)
    integer(ip)                         :: ipoin,iboun,igaub,inode
    integer(ip)                         :: i_field,i_integrals
    integer(ip)                         :: n_field,n_integrals,ngrou
    integer(ip)                         :: pnode,pgaub,pblty,idime,kdime
    integer(ip)                         :: inodb,ielem,pnodb
    real(rp)                            :: gpvol,gpdet
    real(rp)                            :: elcod(ndime,mnode)
    real(rp)                            :: bocod(ndime,mnodb)
    real(rp)                            :: baloc(ndime,ndime)
    real(rp),             allocatable   :: gpfie(:)
    real(rp),             allocatable   :: gradf(:,:)
    real(rp),             allocatable   :: elfie(:,:,:) 
    real(rp),             allocatable   :: integrands(:)

    n_integrals = size(array_func,1)
    n_field     = size(array_func,2)
    ngrou       = size(fields_integrals,2)
    kdime       = 0
    
    if( size(array_fields)       < n_field    ) call runend('MOD_INTEGRALS: WRONG DIMENSION')
    if( size(fields_integrals,1) < n_integrals) call runend('MOD_INTEGRALS: WRONG DIMENSION')
    
    allocate(integrands(n_integrals))

    if( INOTEMPTY ) then
       do i_field = 1,n_field      
          kdime = max(kdime,size(array_fields(i_field)%a,DIM=1,KIND=ip))
          if( size(array_fields(i_field)%a,DIM=2,KIND=ip) < npoin ) then
             call runend('MOD_INTEGRALS: WRONG FIELD DIMENSION')
          end if      
       end do      
       allocate(elfie(kdime,mnodb,n_field))
       allocate(gpfie(kdime))
       allocate(gradf(ndime,kdime))
    end if

    fields_integrals(1:n_integrals,1:ngrou) = 0.0_rp
    !
    ! Compute integral over the whole domain
    !
    elements: do iboun = 1,nboun
       if( legro(iboun) > ngrou ) call runend('MOD_INTEGRALS: WRONG GROUP NUMBER')
       if( legro(iboun) /= 0 ) then
          pblty = ltypb(iboun)
          ielem = lelbo(iboun)
          if (pblty > 0) then
             pnodb = nnode(pblty)
             pgaub = ngaus(pblty)
             pnode = nnode(abs(ltype(ielem)))
             !
             ! Gather operations
             !
             do inodb = 1,pnodb
                ipoin = lnodb(inodb,iboun)
                do i_field = 1, n_field
                   kdime = size(array_fields(i_field) % a,1)
                   elfie(1:kdime,inodb,i_field) = array_fields(i_field) % a(1:kdime,ipoin)
                end do
                bocod(1:ndime,inodb) = coord(1:ndime,ipoin)
             end do
             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                elcod(1:ndime,inode) = coord(1:ndime,ipoin)               
             end do
             !
             ! 1st order Cartesian derivatives, and dV:=GPVOL=|J|*wg
             !
             do igaub = 1, pgaub
                 call bouder(&
                      pnodb,ndime,ndimb,elmar(pblty)%deriv(:,:,igaub),&    ! Cartesian derivative
                      bocod,baloc,gpdet)                                   ! and Jacobian
                gpvol = elmar(pblty)%weigp(igaub) * gpdet                  ! |J|*wg
                call chenor(pnode,baloc,bocod,elcod)                       ! Check normal
                
                integrands = 1.0_rp
                do i_field = 1, n_field
                   gpfie      = 0.0_rp
                   gradf      = 0.0_rp
                   kdime      = size(array_fields(i_field) % a,1)
   
                   do inodb = 1, pnodb
                      do idime = 1,kdime
                         gpfie(idime) = gpfie(idime) + elmar(pblty) % shape(inodb,igaub) * elfie(idime,inodb,i_field)
                      end do
                   end do 
                   
                   do i_integrals = 1,n_integrals
                      integrands(i_integrals) = integrands(i_integrals) * array_func(i_integrals,i_field) % f(gpfie(1:kdime),baloc)
                   end do
                                   
                end do
                
                do i_integrals = 1,n_integrals
                   fields_integrals(i_integrals,legro(iboun)) = fields_integrals(i_integrals,legro(iboun)) + gpvol * integrands(i_integrals)
                end do
   
             end do
          end if
       end if
    end do elements

    call PAR_SUM(n_integrals,ngrou,fields_integrals)

    if(allocated(integrands)) deallocate(integrands)
    if(allocated(elfie))      deallocate(elfie)
    if(allocated(gpfie))      deallocate(gpfie)
    if(allocated(gradf))      deallocate(gradf)
    
  end subroutine integrals_boundary
  
end module mod_integrals
!> @}
