!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    mod_nsi_corio_stab.f90
!> @author  Herbert Owen
!> @date    10/01/2020
!> @brief   Orthogonal Divergence Coriolis Stabilisation
!> @details Orthogonal Divergence Coriolis Stabilisation
!>          
!-----------------------------------------------------------------------

module mod_nsi_corio_stab

#include "def_vector_size.inc"
  use def_parame
  use def_master
  use def_domain
  use def_nastin
  use mod_memory,              only : memory_alloca,memory_deallo
  use mod_element_integration, only : element_shape_function_derivatives_jacobian
  implicit none

  private

  integer(ip), private                   :: idime


  public :: nsi_corio_stab_element_operations
  public :: nsi_corio_stab_matrix
  public :: element_length_corio
  public :: nsi_element_stabilization_corio
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  Herbert Owen
  !> @date    2020-01-10
  !> @brief   Contribution to elauu -  Coriolis galerkin & stab  term
  !> @details Contribution to elauu -  Coriolis galerkin & stab  term
  !> 
  !-----------------------------------------------------------------------

  subroutine nsi_corio_stab_element_operations(pnode,pgaus,gpvol,gpden,gpsha,gpcar,gpstcor,elprdivcor,elrbu,elauu)
    !
    ! Calculating taucor inside here does not seem a good option. This subrotine is called from nsi_element_assembly_split_oss_default & nsi_element_operations_fast5
    ! in the first subroutine all of the other stabilisation paarmeteres come as input thus bringing taucor would imply altering the structure by at least two levels
    ! one option could be to create a specific subroutine calc_taucor that is called earlier in both cases
    !
    !
    implicit none

    integer(ip),intent(in)                                 :: pnode,pgaus
    real(rp),intent(in)                                    :: gpvol(VECTOR_SIZE,pgaus)
    real(rp),intent(in)                                    :: gpden(VECTOR_SIZE,pgaus)
    
    real(rp),intent(in)                                    :: gpsha(VECTOR_SIZE,pnode,pgaus)
    real(rp),intent(in)                                    :: gpcar(VECTOR_SIZE,ndime,mnode,pgaus)
    real(rp),intent(in)                                    :: gpstcor(VECTOR_SIZE,pgaus)

    real(rp),intent(in)                                    :: elprdivcor(VECTOR_SIZE,pnode)

    real(rp),intent(inout)                                 :: elrbu(VECTOR_SIZE,ndime,pnode)
    real(rp),intent(inout)                                 :: elauu(VECTOR_SIZE,ndime*pnode,ndime*pnode)
    !yor! for memory padding purpose : 4*((pnode+3)/4) will give the closest number >= pnode which is a multiple of 4.
    real(rp), dimension(VECTOR_SIZE,4*((pnode+3)/4),pnode) :: elauu11,elauu22,elauu33
    real(rp), dimension(VECTOR_SIZE,4*((pnode+3)/4),pnode) :: elauu12,elauu23,elauu31
    real(rp), dimension(VECTOR_SIZE,4*((pnode+3)/4),pnode) :: elauu12g,elauu23g,elauu31g   ! I will treat the galerkin part separatelly so that I can take advange that the stab is symmetric
    real(rp), dimension(VECTOR_SIZE,4*((pnode+3)/4),pnode) :: elauu21g,elauu32g,elauu13g   ! that I can take advange that the stabilization is symmetric - rethink if is best for speed
    real(rp)                                               :: pdivcor(VECTOR_SIZE,ndime,pnode)   ! here I may put 4 insted of ndime if it is faster -- see yor 
    real(rp)                                               :: gpprdivcor(VECTOR_SIZE) 
    real(rp)                                               :: fact0(VECTOR_SIZE),fact1(VECTOR_SIZE),fact2(VECTOR_SIZE),fact3(VECTOR_SIZE)
    integer(ip)                                            :: igaus,inode,jnode

    real(rp)                                               :: rmom2(VECTOR_SIZE,ndime,ndime,pnode)
    real(rp)                                               :: factvec1(VECTOR_SIZE,4*((pnode+3)/4))
    
#ifdef OPENACC
#define DEF_VECT ivect
#else
#define DEF_VECT 1:VECTOR_SIZE
#endif
    !
    ! Coriolis stabilization --  ( tau_cor * div(omega X u) , div(omega X v) ) -- also Coriolis galerkin term  
    !
    if( corio_nsi > 1.0e-12_rp .and. ndime == 3_ip ) then ! 2d not implemented -- beware I should put if staco_corio_nsi further down - but I do not want to put for optim reasons decide later
       ! in any case with stac_nsi(5) = 0 in the end it adds nothing
       !it is a balance between haveing an if -- adds time when teher is stab
       ! not having an if and doing thing that are multiplied by 0
       elauu11(DEF_VECT,:,:) = 0.0_rp
       elauu22(DEF_VECT,:,:) = 0.0_rp
       elauu33(DEF_VECT,:,:) = 0.0_rp
       elauu12(DEF_VECT,:,:) = 0.0_rp
       elauu23(DEF_VECT,:,:) = 0.0_rp
       elauu31(DEF_VECT,:,:) = 0.0_rp
       elauu12g(DEF_VECT,:,:) = 0.0_rp
       elauu23g(DEF_VECT,:,:) = 0.0_rp
       elauu31g(DEF_VECT,:,:) = 0.0_rp
       elauu21g(DEF_VECT,:,:) = 0.0_rp
       elauu32g(DEF_VECT,:,:) = 0.0_rp
       elauu13g(DEF_VECT,:,:) = 0.0_rp
       
       do igaus = 1,pgaus
          !
          ! Coriolis force - Galerkin term 
          !
          ! 2*rho*(w x u)
          ! x-equation: w x u = wy*uz - wz*uy
          ! y-equation: w x u = wz*ux - wx*uz
          ! z-equation: w x u = wx*uy - wy*ux
          ! Borrowed from nsi_element_operations - here recalculated for each igaus
          !
          rmom2 = 0.0_rp
          fact0(DEF_VECT) = 2.0_rp * gpden(DEF_VECT,igaus)  * gpvol(DEF_VECT,igaus)   ! este 2 tengo que revisarlo !!! en otras partes creo no esta 
          fact1(DEF_VECT) = fact0(DEF_VECT) * fvela_nsi(1)
          fact2(DEF_VECT) = fact0(DEF_VECT) * fvela_nsi(2)
          fact3(DEF_VECT) = fact0(DEF_VECT) * fvela_nsi(3)
          do inode = 1,pnode
             rmom2(DEF_VECT,1,2,inode) = - fact3(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! -wz*uy
             rmom2(DEF_VECT,1,3,inode) =   fact2(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  !  wy*uz
             rmom2(DEF_VECT,2,1,inode) =   fact3(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  !  wz*ux
             rmom2(DEF_VECT,2,3,inode) = - fact1(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! -wx*uz
             rmom2(DEF_VECT,3,1,inode) = - fact2(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! -wy*ux
             rmom2(DEF_VECT,3,2,inode) =   fact1(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  !  wx*uy
          end do
          !
          ! Borrowed from nsi_element_assembly 
          !
          factvec1(DEF_VECT,1:pnode) = gpsha(DEF_VECT,1:pnode,igaus)
          do jnode = 1,pnode
             do inode = 1,pnode
                elauu21g(DEF_VECT,inode,jnode) = elauu21g(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,2,1,jnode)
                elauu31g(DEF_VECT,inode,jnode) = elauu31g(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,3,1,jnode)
                elauu12g(DEF_VECT,inode,jnode) = elauu12g(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,1,2,jnode)
                elauu32g(DEF_VECT,inode,jnode) = elauu32g(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,3,2,jnode)
                elauu13g(DEF_VECT,inode,jnode) = elauu13g(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,1,3,jnode)
                elauu23g(DEF_VECT,inode,jnode) = elauu23g(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,2,3,jnode)
             end do
          end do
          !
          ! pdivcor is recalculated for each igaus
          !
          fact0(DEF_VECT) = gpden(DEF_VECT,igaus)  !* gpvol(DEF_VECT,igaus)    ! ojo este gpvol estaba paar el culo 
          fact1(DEF_VECT) = fact0(DEF_VECT) * fvela_nsi(1)
          fact2(DEF_VECT) = fact0(DEF_VECT) * fvela_nsi(2)
          fact3(DEF_VECT) = fact0(DEF_VECT) * fvela_nsi(3)
          fact0(DEF_VECT) = gpstcor(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)   ! here I reusa it 

          gpprdivcor(DEF_VECT) = 0.0_rp
          do inode = 1,pnode 
             !
             pdivcor(DEF_VECT,1,inode) = fact3(DEF_VECT) * gpcar(DEF_VECT,2,inode,igaus) - fact2(DEF_VECT) * gpcar(DEF_VECT,3,inode,igaus) ! rho * (wz*d_y - wy* d_z) ux
             pdivcor(DEF_VECT,2,inode) = fact1(DEF_VECT) * gpcar(DEF_VECT,3,inode,igaus) - fact3(DEF_VECT) * gpcar(DEF_VECT,1,inode,igaus) ! rho * (wx*d_z - wz* d_x) uy
             pdivcor(DEF_VECT,3,inode) = fact2(DEF_VECT) * gpcar(DEF_VECT,1,inode,igaus) - fact1(DEF_VECT) * gpcar(DEF_VECT,2,inode,igaus) ! rho * (wy*d_x - wx* d_y) uz

             gpprdivcor(DEF_VECT) = gpprdivcor(DEF_VECT) + gpsha(DEF_VECT,inode,igaus) * elprdivcor(DEF_VECT,inode)  ! gauss point value of the projection of the coriolis term

          end do
          !
          ! Coriolis stabilization
          
          !
          do inode = 1,pnode
             do jnode = 1,pnode
                elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode)  + & ! Auu_xx
                     fact0(DEF_VECT) * pdivcor(DEF_VECT,1,inode) * pdivcor(DEF_VECT,1,jnode)   
                elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode)  + & ! Auu_yy
                     fact0(DEF_VECT) * pdivcor(DEF_VECT,2,inode) * pdivcor(DEF_VECT,2,jnode)   
                elauu33(DEF_VECT,inode,jnode) = elauu33(DEF_VECT,inode,jnode)  + & ! Auu_zz
                     fact0(DEF_VECT) * pdivcor(DEF_VECT,1,inode) * pdivcor(DEF_VECT,1,jnode)   
                elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode)  + & ! Auu_xy = Auu_yx
                     fact0(DEF_VECT) * pdivcor(DEF_VECT,1,inode) * pdivcor(DEF_VECT,2,jnode)   
                elauu23(DEF_VECT,inode,jnode) = elauu23(DEF_VECT,inode,jnode)   + & ! Auu_yz = Auu_zy
                     fact0(DEF_VECT) * pdivcor(DEF_VECT,2,inode) * pdivcor(DEF_VECT,3,jnode)   
                elauu31(DEF_VECT,inode,jnode) = elauu31(DEF_VECT,inode,jnode)   + & ! Auu_zx = Auu_xz
                     fact0(DEF_VECT) * pdivcor(DEF_VECT,3,inode) * pdivcor(DEF_VECT,1,jnode)
             end do   ! jnode

             do idime=1,ndime
                elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime, inode)     + &                                                      
                     fact0(DEF_VECT) * pdivcor(DEF_VECT,idime,inode) * gpprdivcor(DEF_VECT) 
             end do
          end do !inode
       end do  !igauss

       do jnode = 1,pnode
          do inode = 1,pnode
             elauu(DEF_VECT,3*inode-2,3*jnode-2) = elauu(DEF_VECT,3*inode-2,3*jnode-2) + elauu11(DEF_VECT,inode,jnode)
             elauu(DEF_VECT,3*inode-1,3*jnode-2) = elauu(DEF_VECT,3*inode-1,3*jnode-2) + elauu12(DEF_VECT,inode,jnode) + &
                  elauu21g(DEF_VECT,inode,jnode) !21 taking advantage of symmetry 
             elauu(DEF_VECT,3*inode  ,3*jnode-2) = elauu(DEF_VECT,3*inode  ,3*jnode-2) + elauu31(DEF_VECT,inode,jnode) + &
                  elauu31g(DEF_VECT,inode,jnode)
             elauu(DEF_VECT,3*inode-2,3*jnode-1) = elauu(DEF_VECT,3*inode-2,3*jnode-1) + elauu12(DEF_VECT,inode,jnode) + &
                  elauu12g(DEF_VECT,inode,jnode)
             elauu(DEF_VECT,3*inode-1,3*jnode-1) = elauu(DEF_VECT,3*inode-1,3*jnode-1) + elauu22(DEF_VECT,inode,jnode)
             elauu(DEF_VECT,3*inode  ,3*jnode-1) = elauu(DEF_VECT,3*inode  ,3*jnode-1) + elauu23(DEF_VECT,inode,jnode) + &
                  elauu32g(DEF_VECT,inode,jnode)  !32
             elauu(DEF_VECT,3*inode-2,3*jnode  ) = elauu(DEF_VECT,3*inode-2,3*jnode  ) + elauu31(DEF_VECT,inode,jnode) + &
                  elauu13g(DEF_VECT,inode,jnode)  !13
             elauu(DEF_VECT,3*inode-1,3*jnode  ) = elauu(DEF_VECT,3*inode-1,3*jnode  ) + elauu23(DEF_VECT,inode,jnode) + &
                  elauu23g(DEF_VECT,inode,jnode)
             elauu(DEF_VECT,3*inode  ,3*jnode  ) = elauu(DEF_VECT,3*inode  ,3*jnode  ) + elauu33(DEF_VECT,inode,jnode)             
          end do
       end do
    end if

  end subroutine nsi_corio_stab_element_operations

  
  !----------------------------------------------------------------------
  !>
  !> @author  Herbert Owen
  !> @date    10/01/2020
  !> @brief   Coriolis stabilisation matrix for Coriols term projection
  !> @details Coriolis stabilisation matrix for Coriols term projection
  !>
  !----------------------------------------------------------------------

  subroutine nsi_corio_stab_matrix(Scorio_nsi,DEALLOCATE_MATRICES) 
    use mod_ker_proper,        only : ker_proper

    use def_domain,            only : elmar  

    
    implicit none

    real(rp),pointer,               intent(inout) :: Scorio_nsi(:,:)
    logical(lg), optional,          intent(in)    :: DEALLOCATE_MATRICES
    integer(ip)                                   :: ielem
    integer(ip)                                   :: pnode
    integer(ip)                                   :: plapl
    integer(ip)                                   :: pgaus
    integer(ip)                                   :: pelty
    integer(ip)                                   :: inode
    integer(ip)                                   :: ipoin

    real(rp)                                      :: elcod(ndime,mnode)

    real(rp)                                      :: gpvol(mgaus)
    real(rp)                                      :: gpsha(mnode,mgaus)
    real(rp)                                      :: gpder(ndime,mnode,mgaus)
    real(rp)                                      :: gpcar(ndime,mnode,mgaus)
    real(rp)                                      :: gphes(ntens,mnode,mgaus)
    logical(lg)                                   :: if_deallocate_matrices
    real(rp)                                      :: gpden(mgaus)                          ! Density    
    
!    integer(ip)                                   :: dummi

    if ( ndime == 2_ip ) call runend('ndime == 2_ip not ready for corils stab')
    
    if( INOTMASTER ) then

       if( present(DEALLOCATE_MATRICES) ) then
          if_deallocate_matrices = DEALLOCATE_MATRICES
       else
          if_deallocate_matrices = .false.
       end if

       if( if_deallocate_matrices ) then
          call memory_deallo(mem_modul(1:2,modul),'SCORIO','nsi_corio_stab_matrix',Scorio_nsi)
       end if

       if( .not. associated(Scorio_nsi) ) then
          call memory_alloca(mem_modul(1:2,modul),'SCORIO','nsi_corio_stab_matrix',Scorio_nsi,ndime,nzdom)
       end if
       Scorio_nsi = 0.0_rp

       do ielem = 1,nelem
          pelty = ltype(ielem)
          if( pelty > 0 ) then
             pgaus = lgaus(ielem)
             pnode = lnnod(ielem)
             plapl = 0
             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                elcod(1:ndime,inode) = coord(1:ndime,ipoin)
             end do
             
             call element_shape_function_derivatives_jacobian(&
                  pnode,pgaus,plapl,elmar(pelty) % weigp,elmar(pelty) % shape,&
                  elmar(pelty) % deriv,elmar(pelty) % heslo,&
                  elcod,gpvol,gpsha,gpder,gpcar,gphes,ielem)

             !             call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,gpsha,gpcar)     ! rho
             gpden = 1.0_rp   ! truchada temporal pra gabls no hay prob porque densi=1.0

             call nsi_corio_stab_element_matrix(&
                 pgaus,pnode,ielem,lnods(:,ielem),gpvol,gpsha,gpcar,gpden,Scorio_nsi)
!                  pgaus,pnode,ielem,lnods(:,ielem),gpvol,gpsha,gpcar,gpden,gpvis,gppor,gpadv,chale,Scorio_nsi)
          end if
       end do

    else
       !
       ! Allocate of size 1 otherwise I got : Attempt to usar pointer when it is not associated with a target
       !
       if( .not. associated(Scorio_nsi) ) then
          call memory_alloca(memor_dom,'SCORIO','nsi_corio_stab_matrix',Scorio_nsi,ndime,1_ip)
       else
          Scorio_nsi = 0.0_rp
       end if

    end if

  end subroutine nsi_corio_stab_matrix




  
  !----------------------------------------------------------------------
  !>
  !> @author  Herbert Owen
  !> @date    10/01/2020
  !> @brief   Obtain Chale for Coriolis stabilization 
  !> @details Obtain Chale for Coriolis stabilization
  !>
  !----------------------------------------------------------------------


  subroutine element_length_corio(ndime,pnode,pelty,porde,elcod,chale)
    use def_domain,            only : elmar
    use mod_elmgeo,            only : element_type
    use mod_elmgeo,            only : elmgeo_element_characteristic_length

    implicit none
    integer(ip),    intent(in)  :: ndime,pnode,pelty,porde
    real(rp),       intent(in)  :: elcod(VECTOR_SIZE,ndime,pnode)
    real(rp),       intent(out) :: chale(VECTOR_SIZE,2)

    real(rp)  :: hleng(VECTOR_SIZE,ndime)   ! I leave it internal -- no need for it to go out

    call elmgeo_element_characteristic_length(&
         ndime,pnode,elmar(pelty) % dercg(:,:),elcod,hleng,element_type(pelty) % natural_length)   
    !
    ! Minimum element length  -- here I give no options  - instead in nsi_element_operations  -- split_oss all options are possible
    !
    chale(1:VECTOR_SIZE,1) = hleng(1:VECTOR_SIZE,ndime)
    chale(1:VECTOR_SIZE,2) = chale(1:VECTOR_SIZE,1)
    !
    ! Divide h by 2 for quadratic elements and 3 for cubic elements
    !
    chale(1:VECTOR_SIZE,1) = chale(1:VECTOR_SIZE,1) / real(porde,rp)
    chale(1:VECTOR_SIZE,2) = chale(1:VECTOR_SIZE,2) / real(porde,rp)
    
  end subroutine element_length_corio



  subroutine nsi_element_stabilization_corio(&
       pgaus,pnode,chale,gpadv,gpvis,gpden,  &
       gppor,gpstcor)
    
    use def_kintyp, only       :  ip,rp
    use def_domain, only       :  ndime
    use def_nastin, only       :  staco_nsi,&
         &                        kfl_taust_nsi,&
         &                        corio_nsi
!    use def_nastin, only       :  kfl_stabi_nsi
    use mod_tauadr, only       :  tauadr
    implicit none
    integer(ip), intent(in)    :: pgaus,pnode
    real(rp),    intent(in)    :: chale(VECTOR_SIZE,2)
    real(rp),    intent(in)    :: gpadv(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(in)    :: gpvis(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpden(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gppor(VECTOR_SIZE,pgaus)
    real(rp),    intent(out)   :: gpstcor(VECTOR_SIZE,pgaus)

    real(rp)                   :: gpst1

    
    integer(ip)                :: igaus,ivect
    real(rp)                   :: adv,dif,rea,h2,gpvno
    
    !    if( kfl_stabi_nsi == NSI_GALERKIN ) then      ! for the moment coriolis Stab  is always introduced 
    
    !       gpstcor = 0.0_rp
    
    !    else
    
    !----------------------------------------------------------------------
    !
    ! TAU1 and TAU2
    !
    !----------------------------------------------------------------------
    
    do igaus = 1,pgaus
       do ivect = 1,VECTOR_SIZE
          !
          ! tau1 = 1 / ( 4*mu/h^2 + 2*rho*uc/h + rho*|w| + sig )
          ! tau2 = 1 / [ (eps + 1/(h^2*4*mu/h^2 + 2*rho*uc/h + rho*|w|)) ]    ! esto esta mal  escrito!!!!!
          !
          gpvno = sqrt(dot_product(gpadv(ivect,1:ndime,igaus),gpadv(ivect,1:ndime,igaus)))
          adv   = gpden(ivect,igaus)*gpvno                                         ! Convective term: rho*|u+u'|
          dif   = gpvis(ivect,igaus)                                               ! Viscous term:    mu
          rea   = gpden(ivect,igaus)*corio_nsi + abs(gppor(ivect,igaus))           ! Coriolis: w + Porosity: sig
          h2    = chale(ivect,2) * chale(ivect,2)
          call tauadr(&
               kfl_taust_nsi,staco_nsi,adv,dif,rea,&
               chale(ivect,1),chale(ivect,2),gpst1)   
          gpstcor(ivect,igaus) = staco_corio_nsi * h2 * gpst1

       end do
    end do
    !    end if

  end subroutine nsi_element_stabilization_corio













  
  !----------------------------------------------------------------------
  !>
  !> @author  Herbert Owen
  !> @date    10/01/2020
  !> @brief   Coriolis stabilisation elemental matrix for Coriols term projection 
  !> @details Coriolis stabilisation elemental matrix for Coriols term projection
  !>
  !----------------------------------------------------------------------
  
  subroutine nsi_corio_stab_element_matrix(&
       pgaus,pnode,ielem,lnods,gpvol,gpsha,gpcar,gpden,Scorio_nsi)
    use def_kermod, only                          :  kfl_element_to_csr
    implicit none

    integer(ip),                    intent(in)    :: pgaus
    integer(ip),                    intent(in)    :: pnode
    integer(ip),                    intent(in)    :: ielem
    integer(ip),                    intent(in)    :: lnods(pnode)
    real(rp),                       intent(in)    :: gpden(pgaus)
    real(rp),                       intent(in)    :: gpvol(pgaus)
    real(rp),                       intent(in)    :: gpsha(pnode,pgaus)
    real(rp),                       intent(in)    :: gpcar(ndime,mnode,pgaus)

    !    real(rp),                       intent(in)    :: gpadv(ndime,pgaus)
    !    real(rp),                       intent(in)    :: gpvis(pgaus)
    !    real(rp),                       intent(in)    :: gppor(pgaus)
    !    real(rp),                       intent(in)    :: chale(2)

    real(rp),    pointer,           intent(inout) :: Scorio_nsi(:,:)
    real(rp)                                      :: pdivcor(ndime,pnode)   ! here I may put 4 insted of ndime if it is faster -- see yor 
    real(rp)                                      :: elscor(pnode,ndime*pnode) ! matrix for obtaining projection of coriolis term 
    real(rp)                                      :: fact0,fact1,fact2,fact3
    integer(ip)                                   :: igaus,inode,idime,iz,jdofn
    integer(ip)                                   :: jcolu,jnode,ipoin,jpoin,jdime
    !
    ! Element matrix
    !
    elscor = 0.0_rp
    do igaus = 1,pgaus

       fact0 = gpden(igaus)
       fact1 = fact0 * fvela_nsi(1)
       fact2 = fact0 * fvela_nsi(2)
       fact3 = fact0 * fvela_nsi(3)

       do inode = 1,pnode 
          !
          pdivcor(1,inode) = fact3 * gpcar(2,inode,igaus) - fact2 * gpcar(3,inode,igaus) ! rho * (wz*d_y - wy* d_z) ux
          pdivcor(2,inode) = fact1 * gpcar(3,inode,igaus) - fact3 * gpcar(1,inode,igaus) ! rho * (wx*d_z - wz* d_x) uy
          pdivcor(3,inode) = fact2 * gpcar(1,inode,igaus) - fact1 * gpcar(2,inode,igaus) ! rho * (wy*d_x - wx* d_y) uz

       end do

       do jnode = 1,pnode
          do jdime = 1,ndime
             jdofn = (jnode-1)*ndime+jdime
             do inode = 1,pnode
                elscor(inode,jdofn) = elscor(inode,jdofn) &
                     - gpvol(igaus) * gpsha(inode,igaus) * pdivcor(jdime,jnode)
             end do
          end do
       end do
    end do
    !
    ! Assemble global system
    !
    if( kfl_element_to_csr == 1 ) then

       do inode = 1,pnode
          do jnode = 1,pnode
             iz = lezdo(inode,jnode,ielem)
             do idime = 1,ndime
                jdofn = (jnode-1) * ndime + idime
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
                Scorio_nsi(idime,iz) = Scorio_nsi(idime,iz) + elscor(inode,jdofn)   !D
             end do
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
                   jdofn = (jnode-1) * ndime + idime
#ifdef NO_COLORING
                   !$OMP ATOMIC
#endif
                   Scorio_nsi(idime,iz) = Scorio_nsi(idime,iz) + elscor(inode,jdofn)  ! D
                end do
             end if
          end do
       end do

    end if

  end subroutine nsi_corio_stab_element_matrix


end module mod_nsi_corio_stab
!> @}
