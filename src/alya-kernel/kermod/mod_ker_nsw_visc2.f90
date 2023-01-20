!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kernel
!> @{
!> @file    mod_ker_nsw_visc2
!> @author  Herbert Owen
!> @brief   Obtains nodal viscosity for no slip wall law  for use in nsi_outtan & nsi_outset
!> @details - 
!>          - 
!> @} 
!-----------------------------------------------------------------------
module mod_ker_nsw_visc2
#include "def_vector_size.inc"
  use def_kintyp,                   only : ip,rp
  use def_kermod,                   only : kfl_nswel_ker
  use mod_ker_nsw_visc  
  
  implicit none
  private
  public :: ker_nod_nsw_visc_0

  contains

    subroutine ker_nod_nsw_visc(pnode,pgaus,list_elements)   ! borrowed from nsi_element_operations
      !
      ! obtains el_nsw_visc
      !
    use def_kintyp,            only : ip,rp
    use def_domain,            only : ltype
    use def_domain,            only : lorde,ltopo,ndime
    use def_domain,            only : lnods
    use def_domain,            only : lelch,elmar,mnode
    use def_domain,            only : lmate,ntens
    use mod_ker_proper,        only : ker_proper
    use def_elmtyp,            only : ELFEM
    use def_kermod,            only : el_nsw_visc

    use mod_element_integration,      only : element_shape_function_derivatives_jacobian 
  
    implicit none
    !
    ! Input and output variables
    !
    integer(ip), intent(in)          :: pnode                        !< Number of nodes
    integer(ip), intent(in)          :: pgaus                        !< Number of Gauss points
    integer(ip), intent(in)          :: list_elements(VECTOR_SIZE)   !< List of elements
    !
    ! Gather No slip wall law
    !
    real(rp)    :: elcod(VECTOR_SIZE,ndime,pnode)                    ! x
    real(rp)    :: elibopo(VECTOR_SIZE,pnode)
    real(rp)    :: elnnsw(VECTOR_SIZE,ndime)
    real(rp)    :: elavv(VECTOR_SIZE,ndime,pnode)
    real(rp)    :: elywal(VECTOR_SIZE)
    !
    ! Indices and dimensions
    !
    integer(ip) :: ielem,ivect
    integer(ip) :: pevat,dummi,ptopo
    integer(ip) :: pelty,plapl,porde
    integer(ip) :: lmate_loc(VECTOR_SIZE)
    integer(ip) :: lnods_loc(VECTOR_SIZE,pnode)
    integer(ip) :: lelch_loc(VECTOR_SIZE)
    integer(ip) :: list_elements_p(VECTOR_SIZE)                      ! List of elements (always positive)
    !
    ! Gauss point values
    !
    real(rp)    :: gpsha(VECTOR_SIZE,pnode,pgaus)                    ! N
    real(rp)    :: gpder(VECTOR_SIZE,ndime,pnode,pgaus)              ! dN/dsi
    real(rp)    :: gpcar(VECTOR_SIZE,ndime,mnode,pgaus)              ! dN/dxi
    real(rp)    :: gphes(VECTOR_SIZE,ntens,mnode,pgaus)              ! d2N/dxidxj
    real(rp)    :: gpvol(VECTOR_SIZE,pgaus)                          ! w*|J|, |J|
    real(rp)    :: gpvis(VECTOR_SIZE,pgaus)                          ! Viscosity
    real(rp)    :: gpvis_nsw(VECTOR_SIZE,pgaus)                      ! Viscosity for no slip wall
    real(rp)    :: gpden(VECTOR_SIZE,pgaus)                          ! Density
    real(rp)    :: gpmut(VECTOR_SIZE,pgaus)                          ! turbulent viscosity

    integer(ip) :: nsize
    real(rp)    :: timea 

    !--------------------------------------------------------------------
    !
    ! Gather: global to local
    !
    !--------------------------------------------------------------------

    ielem = list_elements(1)
    pelty = abs(ltype(ielem))
    plapl = 0
    porde = lorde(pelty)
    ptopo = ltopo(pelty)
    pevat = ndime * pnode
    nsize = VECTOR_SIZE

    plapl = 0

    call cputim(timea)

    !--------------------------------------------------------------------
    !
    ! List of elements
    !
    ! +----+----+----+----+
    ! | 23 | 24 | 25 | 0  | <= list_elements
    ! +----+----+----+----+
    !
    ! +----+----+----+----+
    ! | 23 | 24 | 25 | 23 | <= list_elements_p
    ! +----+----+----+----+
    !
    !--------------------------------------------------------------------

    list_elements_p = list_elements
    do ivect = 1,VECTOR_SIZE
       ielem = abs(list_elements(ivect))
       if( ielem /= 0 ) then
          lnods_loc(ivect,1:pnode) = lnods(1:pnode,ielem)
          lelch_loc(ivect)         = lelch(ielem)
          lmate_loc(ivect)         = lmate(ielem)
       else
          list_elements_p(ivect)   = list_elements(1)
          lnods_loc(ivect,1:pnode) = 0
          lelch_loc(ivect)         = ELFEM
          lmate_loc(ivect)         = 1
       end if
    end do

    call aux_gather_vector(&
         pnode,pgaus,list_elements,lnods_loc,elcod,elavv,elibopo,elnnsw,elywal)


    !--------------------------------------------------------------------
    !
    ! Element shape functions and derivatives
    !
    !--------------------------------------------------------------------

    call element_shape_function_derivatives_jacobian(&
         pnode,pgaus,plapl,elmar(pelty) % weigp,elmar(pelty) % shape,&
         elmar(pelty) % deriv,elmar(pelty) % heslo,&
         elcod,gpvol,gpsha,gpder,gpcar,gphes,&
         list_elements)

    !--------------------------------------------------------------------
    !
    ! Properties and local time step
    !
    !--------------------------------------------------------------------

    call ker_proper('DENSI','PGAUS',dummi,list_elements_p,gpden,pnode,pgaus,gpsha,gpcar)     ! rho
    call ker_proper('VISCO','PGAUS',dummi,list_elements_p,gpvis,pnode,pgaus,gpsha,gpcar)     ! mu
    call ker_proper('TURBU','PGAUS',dummi,list_elements_p,gpmut,pnode,pgaus,gpsha,gpcar)     ! mut

    !
    ! the value we obtain is gpvis_nsw - but actually it is the same for all gp I take the value of gp=1 as the elemental value
    !
    call ker_nsw_visc(ndime,pnode,pgaus,list_elements,elavv,gpcar,elnnsw,gpvis,&
         gpden,gpmut,elibopo,elywal,elcod,gpvis_nsw)

    do ivect = 1,VECTOR_SIZE
       ielem = abs(list_elements(ivect))
       if( ielem /= 0 ) then
          if (kfl_nswel_ker(ielem) /= 0_ip) el_nsw_visc(ielem) = gpvis_nsw(ivect,1)
       end if
    end do       
    
  end subroutine ker_nod_nsw_visc



  subroutine  ker_nod_nsw_visc_0  ! borrowed from nsi_elmope_all(itask)
    
    !
    ! obtains el_nsw_visc
    !
    
    use def_kintyp,                 only : ip,rp
    use def_domain,                 only : ompss_domain   ! Required by OMPSS
    use def_domain,                 only : lnnod
    use def_domain,                 only : lgaus


    use mod_parall,                 only : list_elements_norace_par
    use mod_parall,                 only : num_subd_norace_par
    use mod_parall,                 only : num_pack_norace_par
!    use mod_parall,                 only : list_elements_par
    use mod_parall,                 only : typ_list_elements_par

    use mod_communications,         only : PAR_BARRIER
#ifdef EXTRAE_SFC
    use extrae_module
#endif
    implicit none

    integer(ip)                          :: isubd,ipack,ielem
    integer(ip)                          :: pnode,pgaus
    integer(ip)                          :: num_subd
    integer(ip),                 pointer :: num_pack(:)
    type(typ_list_elements_par), pointer :: list_elements(:)
    !
    ! Creo que esto es lo que toca aca 
    !
    num_subd      =  num_subd_norace_par
    num_pack      => num_pack_norace_par
    list_elements => list_elements_norace_par
    !-------------------------------------------------------------------
    !
    ! Subgrid scale: no race condition
    !
    !-------------------------------------------------------------------

    do isubd = 1,num_subd    !colors
       do ipack = 1,num_pack(isubd)   !pack = vector

          ielem = list_elements(isubd) % packs(ipack) % l(1)  ! Select first element
          pnode = lnnod(ielem)                                ! Number of nodes
          pgaus = lgaus(ielem)                                ! Number of Gauss points

          call ker_nod_nsw_visc(pnode,pgaus,list_elements(isubd) % packs(ipack) % l)  

       end do
    end do



  end subroutine ker_nod_nsw_visc_0

  subroutine aux_gather_vector(&             ! borrowed from nsi_element_operations_gather_vector
       pnode,pgaus,list_elements,lnods,elcod,elavv,elibopo,elnnsw,elywal)

    use def_kintyp, only       :  ip,rp
    use def_domain, only       :  ndime,coord,ywale,lpoty
    use def_kermod, only       :  kfl_noslw_ker,normal_nsw_ker,avupo_ker, kfl_delta
    use def_kermod, only       :  delta_dom
    
    implicit none
    integer(ip), intent(in)    :: pnode
    integer(ip), intent(in)    :: pgaus
    integer(ip), intent(in)    :: list_elements(VECTOR_SIZE)
    integer(ip), intent(in)    :: lnods(VECTOR_SIZE,pnode)
    real(rp),    intent(out)   :: elcod(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(out)   :: elavv(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(out)   :: elibopo(VECTOR_SIZE,pnode)
    real(rp),    intent(out)   :: elnnsw(VECTOR_SIZE,ndime)
    real(rp),    intent(out)   :: elywal(VECTOR_SIZE)
    
    integer(ip)                :: inode,idime,ipoin,ivect,ielem

    do ivect = 1,VECTOR_SIZE
       ielem = list_elements(ivect)
       if( ielem > 0 ) then
          do inode = 1,pnode
             ipoin = lnods(ivect,inode)
             do idime = 1,ndime
                elcod(ivect,idime,inode)   = coord(idime,ipoin)
             end do
          end do
          !
          ! no slip wall law - I could add flag so that i does it only on those elements where it is needed kfl_nswel_ker(ielem)
          !
          if ( kfl_noslw_ker /= 0_ip ) then
             do inode = 1,pnode
                ipoin = lnods(ivect,inode)
                if (lpoty(ipoin) > 0 ) then
                   elibopo(ivect,inode) = 1.0_rp
                else
                   elibopo(ivect,inode) = 0.0_rp
                end if
                do idime = 1,ndime
!                   if (list_elements(1) == 32 .and. inode == 8) print*,'ipoin,avupo_ker(:,ipoin)',ipoin,avupo_ker(:,ipoin)
                   elavv(ivect,idime,inode) = avupo_ker(idime,ipoin)
                end do
             end do
             if (kfl_nswel_ker(ielem) > 0_ip ) then
                do idime = 1,ndime
                   elnnsw(ivect,idime) = normal_nsw_ker(idime,kfl_nswel_ker(ielem))
                end do
             else
                elnnsw(ivect,:) = 0.0_rp
             end if
             if (kfl_delta == 0) then
                 elywal(ivect) = delta_dom
             else
                 elywal(ivect) = ywale(ielem)
             end if
          end if
       else
          !
          ! Element number is null
          !
          elavv(ivect,:,:)   = 0.0_rp
          elnnsw(ivect,:)    = 0.0_rp
          elibopo(ivect,:)   = 0.0_rp
          elywal(ivect)      = 0.0_rp

       end if

    end do

  end subroutine aux_gather_vector
  
end module mod_ker_nsw_visc2
