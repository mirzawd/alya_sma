!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




!------------------------------------------------------------------------
!> @addtogroup NastinMatrixAssembly
!> @{
!> @file    mod_nsi_element_operations_fast.f90
!> @author  Guillaume Houzeaux
!> @brief   Navier-Stokes system element assembly and other element
!>          calculations
!> @details Elemental operations, according to ITASK using vectorized subroutines:
!>
!------------------------------------------------------------------------

module mod_nsi_boundary_operations_fast

#include "def_vector_size.inc"
  use def_kintyp,             only : ip,rp
  use mod_frivel,             only : frivel
  use def_kermod,             only : kfl_ustar
  use mod_bouder
  implicit none
  private

  public :: nsi_boundary_operations_fast
  public :: nsi_boundary_operations_fast5
  public :: nsi_boundary_operations_fast_bck1
  public :: nsi_boundary_operations_fast_bck2

#ifdef OPENACCHHH
#define DEF_VECT ivect
#else
#define DEF_VECT 1:VECTOR_SIZE
#endif 
contains

  subroutine nsi_boundary_operations_fast(&
       itask,VECTOR_DIM,pnodb,pgaub,list_boundaries)
!$acc routine(frivel) seq

    use def_parame
    use def_elmtyp
    use def_master
    use def_kermod
    use def_domain
    use def_nastin
    use mod_ker_proper
    use mod_nsi_assembly_global_system, only : nsi_assembly_fractional_step_boundary_scalar

    implicit none
    !
    ! Input and output variables
    !
    integer(ip), intent(in) :: itask                                         !< What to do
    integer(ip), intent(in) :: VECTOR_DIM                                    !< Pack size
    integer(ip), intent(in) :: pnodb                                         !< Number of boundary nodes
    integer(ip), intent(in) :: pgaub                                         !< Number of Gauss points
    integer(ip), intent(in) :: list_boundaries(VECTOR_DIM)                   !< List of boundaries

    real(rp)                :: elmat(VECTOR_DIM,mnode*ndime,mnode*ndime)     ! Element matrices
    real(rp)                :: elrhs(VECTOR_DIM,ndime,mnode)

    real(rp)                :: baloc(VECTOR_DIM,ndime,ndime)                 ! Gather
    real(rp)                :: bovel(VECTOR_DIM,ndime,pnodb)
    real(rp)                :: bocod(VECTOR_DIM,ndime,pnodb)
    real(rp)                :: elcod(VECTOR_DIM,ndime,mnode)
    real(rp)                :: elvel(VECTOR_DIM,ndime,mnode)

    real(rp)                :: gbsur(VECTOR_DIM)
    real(rp)                :: eucta(VECTOR_DIM)                             ! Values at Gauss points

    real(rp)                :: gbden(VECTOR_DIM,pgaub)
    real(rp)                :: gbvis(VECTOR_DIM,pgaub)

    real(rp)                :: gbden_tmp(pgaub)
    real(rp)                :: gbvis_tmp(pgaub)

    real(rp)                :: roughness
    integer(ip)             :: ielem,ipoin,inode,idime,jnode,jdime           ! Indices and dimensions
    integer(ip)             :: iboun,igaub,inodb,kdime,idofn
    integer(ip)             :: pblty,ievab,jevab,dummi,jnodb

    real(rp)                :: velsh(VECTOR_DIM,ndime,pnodb*ndime)
    real(rp)                :: tvelo(VECTOR_DIM,ndime)
    real(rp)                :: tvefi(VECTOR_DIM,ndime)
    real(rp)                :: tvedi(VECTOR_DIM,ndime)
    real(rp)                :: vikin(VECTOR_DIM)
    real(rp)                :: tveno(VECTOR_DIM)
    real(rp)                :: velfr(VECTOR_DIM)
    real(rp)                :: avvfr(VECTOR_DIM)
    real(rp)                :: tven2(VECTOR_DIM)
    real(rp)                :: yplus(VECTOR_DIM)
    real(rp)                :: fact1(VECTOR_DIM)
    real(rp)                :: fact2(VECTOR_DIM)
    real(rp)                :: velfi(VECTOR_DIM,ndime)
    real(rp)                :: vewal(VECTOR_DIM,ndime)
    real(rp)                :: tract(VECTOR_DIM,3)    
    real(rp)                :: delta_aux(VECTOR_DIM)
    real(rp)                :: produ(VECTOR_DIM)

    integer(ip)             :: list_boundaries_p(VECTOR_DIM)                 ! List of elements (always positive)
    integer(ip)             :: pnode(VECTOR_DIM)                             ! Number of element nodes
    integer(ip)             :: ivect,ii,jj,kk,iboun0(1)

    integer(ip)             :: pnode_max


    if( itask == 1 ) call runend('MOD_NSI_BOUNDARY_OPERATIONS_FATS: NOT CODED')

    gbden     = 0.0_rp
    gbvis     = 0.0_rp
    elmat     = 0.0_rp
    elrhs     = 0.0_rp
    elvel     = 0.0_rp
    roughness = 0.0_rp

    !--------------------------------------------------------------------
    !
    ! List of elements
    !
    ! +----+----+----+----+
    ! | 23 | 24 | 25 | 0  | <= list_boundaries
    ! +----+----+----+----+
    !
    ! +----+----+----+----+
    ! | 23 | 24 | 25 | 25 | <= list_boundaries_p
    ! +----+----+----+----+
    !
    !--------------------------------------------------------------------

    iboun = list_boundaries(1)
    pblty = abs(ltypb(iboun))

    list_boundaries_p = list_boundaries
    iboun0 = minloc(list_boundaries,list_boundaries==0)
    if(iboun0(1) > 1 ) list_boundaries_p(iboun0(1):VECTOR_DIM) = list_boundaries(iboun0(1)-1)

    do ivect = 1,VECTOR_DIM
       iboun        = list_boundaries_p(ivect)
       ielem        = lelbo(iboun)
       pnode(ivect) = lnnod(ielem)    
    end do
    pnode_max = maxval(pnode)

    if( maxval( kfl_fixbo_nsi(list_boundaries_p(:))) == 0 ) return

    !--------------------------------------------------------------------
    !
    ! Properties
    !
    !--------------------------------------------------------------------

    do ivect = 1,VECTOR_DIM
       iboun = list_boundaries_p(ivect)    
       call ker_proper('DENSI','PGAUB',dummi,iboun,gbden_tmp)
       call ker_proper('VISCO','PGAUB',dummi,iboun,gbvis_tmp)
       do igaub = 1,pgaub
          gbden(ivect,igaub) = gbden_tmp(igaub)
          gbvis(ivect,igaub) = gbvis_tmp(igaub)
       end do
    end do

    !--------------------------------------------------------------------
    !
    ! Gather operations: ELVEL, ELCOD, BOVEL, DELTA_AUX
    !
    !--------------------------------------------------------------------

    do ivect = 1,VECTOR_DIM
       iboun = list_boundaries_p(ivect)
       ielem = lelbo(iboun)
       do inode = 1,pnode(ivect)
          ipoin = lnods(inode,ielem)
          do idime = 1,ndime
             elcod(ivect,idime,inode) = coord(idime,ipoin)
             elvel(ivect,idime,inode) = veloc(idime,ipoin,1)
          end do
       end do
       do inodb = 1,pnodb
          ipoin = lnodb(inodb,iboun)
          do idime = 1,ndime
             bocod(ivect,idime,inodb) = coord(idime,ipoin)
             bovel(ivect,idime,inodb) = veloc(idime,ipoin,1)
          end do
       end do
       if ( kfl_waexl_ker == 0_ip ) then                            ! normal behaviour
          if( kfl_delta == 1 ) then
             delta_aux(ivect) = ywalb(iboun)                        ! variable wall distance
          else
             delta_aux(ivect) = delta_nsi                           ! fixed wall distance
          end if
       end if
    end do

    !--------------------------------------------------------------------
    !
    ! Shape functions
    !
    !--------------------------------------------------------------------
    !do ivect = 1,VECTOR_DIM
    !   iboun = list_boundaries_p(ivect)
    !   do igaub = 1,pgaub
    !      call bouder(&
    !         pnodb,ndime,ndimb,elmar(pblty)%deriv(1,1,igaub),&    ! Cartesian derivative
    !         bocod(ivect,:,:),baloc_tmp,eucta)                    ! and Jacobian
    !      call chenor(&
    !         pnode,baloc_tmp,bocod(ivect,:,:),elcod(ivect,:,:))   ! Normalize normal
    !      gbsur(ivect,igaub)     = elmar(pblty)%weigp(igaub)*eucta
    !      baloc(ivect,:,:,igaub) = baloc_tmp(:,:)
    !   end do
    !end do

#ifdef OPENACCHHH
!!$acc enter data create(agrau , lnods_loc , gpvol , elauu ,      &
!!$acc    eldtrho , elvel , gpmut     , gpadv , gpvel ,      &
!!$acc    gpgve   , elcod , gpsha     , gpcar ,              &
!!$acc    gpden   , xjaci , wgrgr     ,                      &
!!$acc    gpdet   , gprhs , elmurho   , elrbu , gpvis     )  &
!!$acc    copyin( weigp  , list_elements , sha_aux)
    !
    !$acc parallel loop gang vector default(present)
    !
    do ivect = 1,VECTOR_DIM
#endif       

       !--------------------------------------------------------------------
       !
       ! Element matrix
       !
       !--------------------------------------------------------------------

       do igaub = 1,pgaub
          !
          ! Shape functions
          !
          baloc(DEF_VECT,:,:) = 0.0_rp
          do ii=1,ndime
             do jj=1,ndimb
                do kk=1,pnodb
                   baloc(DEF_VECT,ii,jj)=baloc(DEF_VECT,ii,jj)&
                        +bocod(DEF_VECT,ii,kk)*elmar(pblty)%deriv(jj,kk,igaub)
                end do
             end do
          end do
          if( ndime == 3 ) then
             baloc(DEF_VECT,1,3) =   baloc(DEF_VECT,2,1) * baloc(DEF_VECT,3,2) &
                  &                - baloc(DEF_VECT,3,1) * baloc(DEF_VECT,2,2)
             baloc(DEF_VECT,2,3) =   baloc(DEF_VECT,3,1) * baloc(DEF_VECT,1,2) &
                  &                - baloc(DEF_VECT,1,1) * baloc(DEF_VECT,3,2)
             baloc(DEF_VECT,3,3) =   baloc(DEF_VECT,1,1) * baloc(DEF_VECT,2,2) &
                  &                - baloc(DEF_VECT,2,1) * baloc(DEF_VECT,1,2)
             eucta(DEF_VECT)     =   baloc(DEF_VECT,1,3) * baloc(DEF_VECT,1,3) &
                  &                + baloc(DEF_VECT,2,3) * baloc(DEF_VECT,2,3) &
                  &                + baloc(DEF_VECT,3,3) * baloc(DEF_VECT,3,3)
             eucta(DEF_VECT)     =   sqrt(eucta(DEF_VECT)) 
             baloc(DEF_VECT,1,1) =   baloc(DEF_VECT,2,3) * baloc(DEF_VECT,3,2) &
                  &                - baloc(DEF_VECT,3,3) * baloc(DEF_VECT,2,2)
             baloc(DEF_VECT,2,1) =   baloc(DEF_VECT,3,3) * baloc(DEF_VECT,1,2) &
                  &                - baloc(DEF_VECT,1,3) * baloc(DEF_VECT,3,2)
             baloc(DEF_VECT,3,1) =   baloc(DEF_VECT,1,3) * baloc(DEF_VECT,2,2) &
                  &                - baloc(DEF_VECT,2,3) * baloc(DEF_VECT,1,2)
          else
             baloc(DEF_VECT,1,2) =   baloc(DEF_VECT,2,1)
             baloc(DEF_VECT,2,2) = - baloc(DEF_VECT,1,1)
             eucta(DEF_VECT)     =   baloc(DEF_VECT,1,2) * baloc(DEF_VECT,1,2) &
                  &                + baloc(DEF_VECT,2,2) * baloc(DEF_VECT,2,2)
             eucta(DEF_VECT)     =   sqrt(eucta(DEF_VECT))             
          end if
          do ii = 1,ndime
             produ(DEF_VECT) = 0.0_rp
             do jj = 1,ndime
                produ(DEF_VECT) = produ(DEF_VECT) + baloc(DEF_VECT,jj,ii)*baloc(DEF_VECT,jj,ii)
             end do
             produ(DEF_VECT) = sqrt(produ(DEF_VECT))
             do jj = 1,ndime
                baloc(DEF_VECT,jj,ii) = produ(DEF_VECT) * baloc(DEF_VECT,jj,ii)
             end do
          end do
          gbsur(DEF_VECT) = elmar(pblty)%weigp(igaub)*eucta(DEF_VECT)

          tract(DEF_VECT,:)   = 0.0_rp
          vewal(DEF_VECT,:)   = 0.0_rp
          velfi(DEF_VECT,:)   = 0.0_rp
          velsh(DEF_VECT,:,:) = 0.0_rp
          tveno(DEF_VECT)     = 0.0_rp

          if ( kfl_waexl_ker == 0_ip ) then  ! normal behaviour
             do inodb = 1,pnodb
                do idime = 1,ndime
                   vewal(DEF_VECT,idime) = vewal(DEF_VECT,idime) &
                        + bovel(DEF_VECT,idime,inodb) * elmar(pblty) % shape(inodb,igaub)
                end do
             end do
          end if

          tvelo(DEF_VECT,:)   = vewal(DEF_VECT,:)
          tvefi(DEF_VECT,:)   = velfi(DEF_VECT,:)
          do idime = 1,ndime                                     ! Tangent veloc.
             do inodb = 1,pnodb
                idofn = (inodb-1)*ndime+idime
                velsh(DEF_VECT,idime,idofn) = velsh(DEF_VECT,idime,idofn) &
                     + elmar(pblty) % shape(inodb,igaub)
             end do
             do kdime = 1,ndime
                tvelo(DEF_VECT,idime) = tvelo(DEF_VECT,idime) &
                     - baloc(DEF_VECT,idime,ndime)   &
                     * baloc(DEF_VECT,kdime,ndime) * vewal(DEF_VECT,kdime)
                tvefi(DEF_VECT,idime) = tvefi(DEF_VECT,idime)   &
                     - baloc(DEF_VECT,idime,ndime) &
                     * baloc(DEF_VECT,kdime,ndime) * velfi(DEF_VECT,kdime)
                do inodb = 1,pnodb
                   idofn = (inodb-1)*ndime+kdime
                   velsh(DEF_VECT,idime,idofn) = velsh(DEF_VECT,idime,idofn) &
                        - baloc(DEF_VECT,idime,ndime)                  &
                        * baloc(DEF_VECT,kdime,ndime) * elmar(pblty) % shape(inodb,igaub)
                end do
             end do
             tvedi(DEF_VECT,idime) = tvelo(DEF_VECT,idime) - tvefi(DEF_VECT,idime)
          end do
          !
          ! Compute U*
          ! 
          vikin(DEF_VECT) = gbvis(DEF_VECT,igaub) / gbden(DEF_VECT,igaub)     ! nu

          do idime = 1,ndime
             tveno(DEF_VECT) = tveno(DEF_VECT) + tvedi(DEF_VECT,idime)**2
          end do
          tveno(DEF_VECT) = sqrt(tveno(DEF_VECT))

#ifndef OPENACCHHH
          do ivect = 1,VECTOR_DIM
#endif             
             iboun = list_boundaries_p(ivect)
!             include './my_frivel_include.txt'
             call frivel(kfl_ustar,delta_aux(ivect),roughness,tveno(ivect),vikin(ivect),velfr(ivect)) 
             avvfr(ivect) = velfr(ivect)
             tven2(ivect) = tveno(ivect)
             !
             ! Compute prescribed traction
             !
             yplus(ivect) = delta_aux(ivect) * velfr(ivect) / vikin(ivect)
             if( yplus(ivect) < 5.0_rp  .and. kfl_ustar == 0 ) then
                fact1(ivect) = gbden(ivect,igaub) * vikin(ivect) / delta_aux(ivect)
             else
                fact1(ivect) = gbden(ivect,igaub) * avvfr(ivect) * velfr(ivect) / tveno(ivect)
             end if
             do idime = 1,ndime
                tract(ivect,idime) = tract(ivect,idime) + fact1(ivect) * tvefi(ivect,idime)
             end do

             if( kfl_waexl_ker == 0 ) then
                do inodb = 1,pnodb
                   fact2(ivect) = fact1(ivect) * elmar(pblty) % shape(inodb,igaub)
                   inode = lboel(inodb,iboun)
                   ievab = (inode-1)*ndime
                   do idime = 1,ndime
                      ievab = ievab+1
                      do jnodb = 1,pnodb
                         jnode = lboel(jnodb,iboun)
                         jevab = (jnode-1)*ndime
                         do jdime = 1,ndime                          
                            jevab = jevab+1
                            elmat(ivect,ievab,jevab) = elmat(ivect,ievab,jevab)     &
                                 + fact2(ivect) * velsh(ivect,idime,(jnodb-1)*ndime+jdime) &
                                 * gbsur(ivect)
                         end do
                      end do
                      elrhs(ivect,idime,inode) = elrhs(ivect,idime,inode)           &
                           + tract(ivect,idime) * elmar(pblty) % shape(inodb,igaub) &
                           * gbsur(ivect) 
                   end do
                end do
             end if
          end do

#ifndef OPENACCHHH
       end do
#endif
       !--------------------------------------------------------------------
       !
       ! Send matrix to RHS
       !
       !--------------------------------------------------------------------

       do jnode = 1,pnode_max
          do jdime = 1,ndime
             jevab = (jnode-1)*ndime+jdime
             do inode = 1,pnode_max
                do idime = 1,ndime
                   ievab = (inode-1)*ndime+idime
                   elrhs(DEF_VECT,idime,inode) = elrhs(DEF_VECT,idime,inode) &
                        - elmat(DEF_VECT,ievab,jevab) * elvel(DEF_VECT,jdime,jnode) 
                end do
             end do
          end do
       end do

       !--------------------------------------------------------------------
       !
       ! Scatter element matrix to global one 
       !
       !--------------------------------------------------------------------

#ifndef OPENACCHHH
       do ivect = 1,VECTOR_DIM
#endif          
          iboun = list_boundaries(ivect)    
          if( iboun > 0 ) then
             if( kfl_fixbo_nsi(iboun) /= 0 ) then
                ielem = lelbo(iboun)  
                do inode = 1,pnode(ivect)
                   ipoin = lnods(inode,ielem)
                   do idime = 1,ndime
                      ievab = idime + (ipoin-1) * ndime 
                      !$acc atomic update
#ifdef NO_COLORING
                      !$OMP ATOMIC
#endif
                      rhsid(ievab) = rhsid(ievab) + elrhs(ivect,idime,inode)

                   end do
                end do
             end if
          end if
#ifndef OPENACCHHH
       end do
#endif           
#ifdef OPENACCHHH
    end do
    !!$acc end parallel loop
    !!$acc end data
#endif       

  end subroutine nsi_boundary_operations_fast

  subroutine nsi_boundary_operations_fast5(&
       itask,VECTOR_DIM,pnodb,pgaub,list_boundaries)
!$acc routine(frivel) seq
    
    use def_parame
    use def_elmtyp
    use def_master
    use def_kermod
    use def_domain
    use def_nastin
    use mod_ker_proper
    use mod_nsi_assembly_global_system, only : nsi_assembly_fractional_step_boundary_scalar
    use def_kermod,            only : densi_ker
    use def_kermod,            only : visco_ker

    implicit none
    !
    ! Input and output variables
    !
    integer(ip), intent(in) :: itask                                         !< What to do
    integer(ip), intent(in) :: VECTOR_DIM                                    !< Pack size
    integer(ip), intent(in) :: pnodb                                         !< Number of boundary nodes
    integer(ip), intent(in) :: pgaub                                         !< Number of Gauss points
    integer(ip), intent(in) :: list_boundaries(VECTOR_DIM)                   !< List of boundaries

    real(rp)                :: elmat(VECTOR_DIM,mnode*ndime,mnode*ndime)     ! Element matrices
    real(rp)                :: elrhs(VECTOR_DIM,ndime,mnode)

    real(rp)                :: baloc(VECTOR_DIM,ndime,ndime)                 ! Gather
    real(rp)                :: bovel(VECTOR_DIM,ndime,pnodb)
    real(rp)                :: bocod(VECTOR_DIM,ndime,pnodb)
    real(rp)                :: elcod(VECTOR_DIM,ndime,mnode)
    real(rp)                :: elvel(VECTOR_DIM,ndime,mnode)

    real(rp)                :: gbsur(VECTOR_DIM)
    real(rp)                :: eucta(VECTOR_DIM)                     ! Values at Gauss points

    real(rp)                :: gbden(VECTOR_DIM,pgaub)
    real(rp)                :: gbvis(VECTOR_DIM,pgaub)
    real(rp)                :: produ(VECTOR_DIM)

!    real(rp)                :: gbden_tmp(pgaub)
!    real(rp)                :: gbvis_tmp(pgaub)

    real(rp)                :: roughness
    integer(ip)             :: ielem,ipoin,inode,idime,jnode,jdime     ! Indices and dimensions
    integer(ip)             :: iboun,igaub,inodb,kdime,idofn
    integer(ip)             :: pblty,ievab,jevab
    integer(ip)             :: jnodb
!    integer(ip)             :: dummi

    real(rp)                :: velsh(VECTOR_DIM,ndime,pnodb*ndime)
    real(rp)                :: tvelo(VECTOR_DIM,ndime)
    real(rp)                :: tvefi(VECTOR_DIM,ndime)
    real(rp)                :: tvedi(VECTOR_DIM,ndime)
    real(rp)                :: vikin(VECTOR_DIM)
    real(rp)                :: tveno(VECTOR_DIM)
    real(rp)                :: velfr(VECTOR_DIM)
    real(rp)                :: avvfr(VECTOR_DIM)
    real(rp)                :: tven2(VECTOR_DIM)
    real(rp)                :: yplus(VECTOR_DIM)
    real(rp)                :: fact1(VECTOR_DIM)
    real(rp)                :: fact2(VECTOR_DIM)
    real(rp)                :: velfi(VECTOR_DIM,ndime)
    real(rp)                :: vewal(VECTOR_DIM,ndime)
    real(rp)                :: tract(VECTOR_DIM,3)    
    real(rp)                :: delta_aux(VECTOR_DIM)

    integer(ip)             :: pnode(VECTOR_DIM)
    integer(ip)             :: list_boundaries_p(VECTOR_DIM)                      ! List of elements (always positive)
    integer(ip)             :: ivect,ii,jj,kk

    integer(ip)             :: pnode_max,iboun0(1)
    real(rp)                :: ldensi,lvisco

    !!    gbden     = 0.0_rp
    !!    gbvis     = 0.0_rp
    !    elmat     = 0.0_rp
    !    elrhs     = 0.0_rp
    !    roughness = 0.0_rp


    ldensi = densi_ker % value_const(1,1)
    lvisco = visco_ker % value_const(1,1)

    !--------------------------------------------------------------------
    !
    ! List of elements
    !
    ! +----+----+----+----+
    ! | 23 | 24 | 25 | 0  | <= list_boundaries
    ! +----+----+----+----+
    !
    ! +----+----+----+----+
    ! | 23 | 24 | 25 | 23 | <= list_boundaries_p
    ! +----+----+----+----+
    !
    !--------------------------------------------------------------------

    iboun = list_boundaries(1)
    pblty = abs(ltypb(iboun))

    list_boundaries_p = list_boundaries
    iboun0 = minloc(list_boundaries,list_boundaries==0)
    if(iboun0(1) > 1 ) list_boundaries_p(iboun0(1):VECTOR_DIM) = list_boundaries(iboun0(1)-1)

    if( maxval( kfl_fixbo_nsi(list_boundaries_p(:))) == 0 ) return

    do ivect = 1,VECTOR_DIM
       iboun        = list_boundaries_p(ivect)
       ielem        = lelbo(iboun)
       pnode(ivect) = lnnod(ielem)    
    end do
    pnode_max = maxval(pnode)
 
    !--------------------------------------------------------------------
    !
    ! Properties
    !
    !--------------------------------------------------------------------

    !!    do ivect = 1,VECTOR_DIM
    !!       iboun = list_boundaries_p(ivect)    
    !!       call ker_proper('DENSI','PGAUB',dummi,iboun,gbden_tmp)
    !!       call ker_proper('VISCO','PGAUB',dummi,iboun,gbvis_tmp)
    !!       do igaub = 1,pgaub
    !!          gbden(ivect,igaub) = gbden_tmp(igaub)
    !!          gbvis(ivect,igaub) = gbvis_tmp(igaub)
    !!       end do
    !!    end do
    !!
    !--------------------------------------------------------------------
    !
    ! Gather operations: ELVEL, ELCOD, BOVEL, DELTA_AUX
    !
    !--------------------------------------------------------------------
!!!    
!!!    do ivect = 1,VECTOR_DIM
!!!       iboun = list_boundaries_p(ivect)
!!!       ielem = lelbo(iboun)
!!!       do inode = 1,pnode
!!!          ipoin = lnods(inode,ielem)
!!!          do idime = 1,ndime
!!!             elcod(ivect,idime,inode) = coord(idime,ipoin)
!!!             elvel(ivect,idime,inode) = veloc(idime,ipoin,1)
!!!          end do
!!!       end do
!!!       do inodb = 1,pnodb
!!!          ipoin = lnodb(inodb,iboun)
!!!          do idime = 1,ndime
!!!             bocod(ivect,idime,inodb) = coord(idime,ipoin)
!!!             bovel(ivect,idime,inodb) = veloc(idime,ipoin,1)
!!!          end do
!!!       end do
!!!       if ( kfl_waexl_ker == 0_ip ) then                            ! normal behaviour
!!!          if( kfl_delta == 1 ) then
!!!             delta_aux(ivect) = ywalb(iboun)                        ! variable wall distance
!!!          else
!!!             delta_aux(ivect) = delta_nsi                           ! fixed wall distance
!!!          end if
!!!       end if
!!!    end do
!!!
    !--------------------------------------------------------------------
    !
    ! Shape functions
    !
    !--------------------------------------------------------------------
    !do ivect = 1,VECTOR_DIM
    !   iboun = list_boundaries_p(ivect)
    !   do igaub = 1,pgaub
    !      call bouder(&
    !         pnodb,ndime,ndimb,elmar(pblty)%deriv(1,1,igaub),&    ! Cartesian derivative
    !         bocod(ivect,:,:),baloc_tmp,eucta)                    ! and Jacobian
    !      call chenor(&
    !         pnode,baloc_tmp,bocod(ivect,:,:),elcod(ivect,:,:))   ! Normalize normal
    !      gbsur(ivect,igaub)     = elmar(pblty)%weigp(igaub)*eucta
    !      baloc(ivect,:,:,igaub) = baloc_tmp(:,:)
    !   end do
    !end do

#ifdef OPENACCHHH

    !$acc enter data create(tract, vewal, velfi , velsh , tveno,    &
    !$acc tvelo,tvefi,tvedi,elcod,elvel,bocod,bovel,delta_aux,      &
    !$acc vikin,velfr,yplus,fact1,fact2,                            &
    !$acc avvfr, tven2,produ,gbsur,eucta,                           &
    !$acc  baloc,elmat,elrhs,gbvis,gbden )                          &
    !$acc copyin(                                                   & 
    !$acc elmar,                                                    &
    !$acc elmar(pblty)%shape,                                       &
    !$acc elmar(pblty)%deriv,                                       &
    !$acc elmar(pblty)%weigp,                                       &
    !$acc list_boundaries_p,list_boundaries,lelbo,        &
    !$acc lboel,kfl_fixbo_nsi,ywalb,lnodb,pnode )  

    !
    !$acc parallel loop gang vector default(present)
    !
    do ivect = 1,VECTOR_DIM
#endif       
       gbden(DEF_VECT,:) = ldensi !densi_ker % value_const(1,1)
       gbvis(DEF_VECT,:) = lvisco !visco_ker % value_const(1,1)

       !--------------------------------------------------------------------
       !
       ! Gather operations: ELVEL, ELCOD, BOVEL, DELTA_AUX
       !
       !--------------------------------------------------------------------
#ifndef OPENACCHHH   
       do ivect = 1,VECTOR_DIM
#endif

          iboun = list_boundaries_p(ivect)
          ielem = lelbo(iboun)
          do inode = 1,pnode(ivect)
             ipoin = lnods(inode,ielem)
             do idime = 1,ndime
                elcod(ivect,idime,inode) = coord(idime,ipoin)
                elvel(ivect,idime,inode) = veloc(idime,ipoin,1)
             end do
          end do
          do inodb = 1,pnodb
             ipoin = lnodb(inodb,iboun)
             do idime = 1,ndime
                bocod(ivect,idime,inodb) = coord(idime,ipoin)
                bovel(ivect,idime,inodb) = veloc(idime,ipoin,1)
             end do
          end do

          !!CHECK THIS IF, MAY CORRUPT MEMORY
          if ( kfl_waexl_ker == 0_ip ) then                            ! normal behaviour
             if( kfl_delta == 1 ) then
                delta_aux(ivect) = ywalb(iboun)                        ! variable wall distance
             else
                delta_aux(ivect) = delta_nsi                           ! fixed wall distance
             end if
          end if
#ifndef OPENACCHHH   
       end do
#endif

       !--------------------------------------------------------------------
       !
       ! Element matrix
       !
       !--------------------------------------------------------------------

       elrhs(DEF_VECT,:,:) = 0.0_rp 
       elmat(DEF_VECT,:,:) = 0.0_rp

       do igaub = 1,pgaub
          baloc(DEF_VECT,:,:) = 0.0_rp 
          !
          ! Shape functions
          !
          do ii=1,ndime
             do jj=1,ndimb
                do kk=1,pnodb
                   baloc(DEF_VECT,ii,jj)=baloc(DEF_VECT,ii,jj)&
                        +bocod(DEF_VECT,ii,kk)*elmar(pblty)%deriv(jj,kk,igaub)
                end do
             end do
          end do
          baloc(DEF_VECT,1,3) =   baloc(DEF_VECT,2,1) * baloc(DEF_VECT,3,2) &
               &                - baloc(DEF_VECT,3,1) * baloc(DEF_VECT,2,2)
          baloc(DEF_VECT,2,3) =   baloc(DEF_VECT,3,1) * baloc(DEF_VECT,1,2) &
               &                - baloc(DEF_VECT,1,1) * baloc(DEF_VECT,3,2)
          baloc(DEF_VECT,3,3) =   baloc(DEF_VECT,1,1) * baloc(DEF_VECT,2,2) &
               &                - baloc(DEF_VECT,2,1) * baloc(DEF_VECT,1,2)
          eucta(DEF_VECT)     =   baloc(DEF_VECT,1,3) * baloc(DEF_VECT,1,3) &
               &                + baloc(DEF_VECT,2,3) * baloc(DEF_VECT,2,3) &
               &                + baloc(DEF_VECT,3,3) * baloc(DEF_VECT,3,3)
          eucta(DEF_VECT)     =   sqrt(eucta(DEF_VECT)) 
          baloc(DEF_VECT,1,1) =   baloc(DEF_VECT,2,3) * baloc(DEF_VECT,3,2) &
               &                - baloc(DEF_VECT,3,3) * baloc(DEF_VECT,2,2)
          baloc(DEF_VECT,2,1) =   baloc(DEF_VECT,3,3) * baloc(DEF_VECT,1,2) &
               &                - baloc(DEF_VECT,1,3) * baloc(DEF_VECT,3,2)
          baloc(DEF_VECT,3,1) =   baloc(DEF_VECT,1,3) * baloc(DEF_VECT,2,2) &
               &                - baloc(DEF_VECT,2,3) * baloc(DEF_VECT,1,2)
          do ii = 1,ndime
             produ(DEF_VECT) = 0.0_rp
             do jj = 1,ndime
                produ(DEF_VECT) = produ(DEF_VECT) + baloc(DEF_VECT,jj,ii)*baloc(DEF_VECT,jj,ii)
             end do
             produ(DEF_VECT) = sqrt(produ(DEF_VECT))
             do jj = 1,ndime
                baloc(DEF_VECT,jj,ii) = produ(DEF_VECT) * baloc(DEF_VECT,jj,ii)
             end do
          end do
          gbsur(DEF_VECT) = elmar(pblty)%weigp(igaub)*eucta(DEF_VECT)

          tract(DEF_VECT,:)   = 0.0_rp
          vewal(DEF_VECT,:)   = 0.0_rp
          velfi(DEF_VECT,:)   = 0.0_rp
          velsh(DEF_VECT,:,:) = 0.0_rp
          tveno(DEF_VECT)     = 0.0_rp

          if ( kfl_waexl_ker == 0_ip ) then  ! normal behaviour
             do inodb = 1,pnodb
                do idime = 1,ndime
                   vewal(DEF_VECT,idime) = vewal(DEF_VECT,idime) &
                        + bovel(DEF_VECT,idime,inodb) * elmar(pblty) % shape(inodb,igaub)
                end do
             end do
          end if

          tvelo(DEF_VECT,:)   = vewal(DEF_VECT,:)
          tvefi(DEF_VECT,:)   = velfi(DEF_VECT,:)
          do idime = 1,ndime                                     ! Tangent veloc.
             do inodb = 1,pnodb
                idofn = (inodb-1)*ndime+idime
                velsh(DEF_VECT,idime,idofn) = velsh(DEF_VECT,idime,idofn) &
                     + elmar(pblty) % shape(inodb,igaub)
             end do
             do kdime = 1,ndime
                tvelo(DEF_VECT,idime) = tvelo(DEF_VECT,idime) &
                     - baloc(DEF_VECT,idime,ndime)   &
                     * baloc(DEF_VECT,kdime,ndime) * vewal(DEF_VECT,kdime)
                tvefi(DEF_VECT,idime) = tvefi(DEF_VECT,idime)   &
                     - baloc(DEF_VECT,idime,ndime) &
                     * baloc(DEF_VECT,kdime,ndime) * velfi(DEF_VECT,kdime)
                do inodb = 1,pnodb
                   idofn = (inodb-1)*ndime+kdime
                   velsh(DEF_VECT,idime,idofn) = velsh(DEF_VECT,idime,idofn) &
                        - baloc(DEF_VECT,idime,ndime)                  &
                        * baloc(DEF_VECT,kdime,ndime) * elmar(pblty) % shape(inodb,igaub)
                end do
             end do
             tvedi(DEF_VECT,idime) = tvelo(DEF_VECT,idime) - tvefi(DEF_VECT,idime)
          end do
          !
          ! Compute U*
          ! 
          vikin(DEF_VECT) = gbvis(DEF_VECT,igaub) / gbden(DEF_VECT,igaub)     ! nu

          do idime = 1,ndime
             tveno(DEF_VECT) = tveno(DEF_VECT) + tvedi(DEF_VECT,idime)**2
          end do
          tveno(DEF_VECT) = sqrt(tveno(DEF_VECT))

#ifndef OPENACCHHH
          do ivect = 1,VECTOR_DIM
#endif             
             iboun = list_boundaries_p(ivect)
!             include './my_frivel_include.txt'
             call frivel(kfl_ustar,delta_aux(ivect),roughness,tveno(ivect),vikin(ivect),velfr(ivect)) 
             avvfr(ivect) = velfr(ivect)
             tven2(ivect) = tveno(ivect)
             !
             ! Compute prescribed traction
             !
             yplus(ivect) = delta_aux(ivect) * velfr(ivect) / vikin(ivect)
             if( yplus(ivect) < 5.0_rp  .and. kfl_ustar == 0 ) then
                fact1(ivect) = gbden(ivect,igaub) * vikin(ivect) / delta_aux(ivect)
             else
                fact1(ivect) = gbden(ivect,igaub) * avvfr(ivect) * velfr(ivect) / tveno(ivect)
             end if
             do idime = 1,ndime
                tract(ivect,idime) = tract(ivect,idime) + fact1(ivect) * tvefi(ivect,idime)
             end do

             if( kfl_waexl_ker == 0 ) then
                do inodb = 1,pnodb
                   fact2(ivect) = fact1(ivect) * elmar(pblty) % shape(inodb,igaub)
                   inode = lboel(inodb,iboun)
                   ievab = (inode-1)*ndime
                   do idime = 1,ndime
                      ievab = ievab+1
                      do jnodb = 1,pnodb
                         jnode = lboel(jnodb,iboun)
                         jevab = (jnode-1)*ndime
                         do jdime = 1,ndime                          
                            jevab = jevab+1
                            elmat(ivect,ievab,jevab) = elmat(ivect,ievab,jevab)     &
                                 + fact2(ivect) * velsh(ivect,idime,(jnodb-1)*ndime+jdime) &
                                 * gbsur(ivect)
                         end do
                      end do
                      elrhs(ivect,idime,inode) = elrhs(ivect,idime,inode)           &
                           + tract(ivect,idime) * elmar(pblty) % shape(inodb,igaub) &
                           * gbsur(ivect) 
                   end do
                end do
             end if
          end do

#ifndef OPENACCHHH
       end do
#endif
       !--------------------------------------------------------------------
       !
       ! Send matrix to RHS
       !
       !--------------------------------------------------------------------

       do jnode = 1,pnode_max
          do jdime = 1,ndime
             jevab = (jnode-1)*ndime+jdime
             do inode = 1,pnode_max
                do idime = 1,ndime
                   ievab = (inode-1)*ndime+idime
                   elrhs(DEF_VECT,idime,inode) = elrhs(DEF_VECT,idime,inode) &
                        - elmat(DEF_VECT,ievab,jevab) * elvel(DEF_VECT,jdime,jnode) 
                end do
             end do
          end do
       end do

       !--------------------------------------------------------------------
       !
       ! Scatter element matrix to global one 
       !
       !--------------------------------------------------------------------

#ifndef OPENACCHHH
       do ivect = 1,VECTOR_DIM
#endif          
          iboun = list_boundaries(ivect)    
          if( iboun > 0 ) then
             if( kfl_fixbo_nsi(iboun) /= 0 ) then
                ielem = lelbo(iboun)  
                do inode = 1,pnode(ivect)
                   ipoin = lnods(inode,ielem)
                   do idime = 1,ndime
                      ievab = idime + (ipoin-1) * ndime 
                      !$acc atomic update
#ifdef NO_COLORING
                      !$OMP ATOMIC
#endif
                      rhsid(ievab) = rhsid(ievab) + elrhs(ivect,idime,inode)

                   end do
                end do
             end if
          end if
#ifndef OPENACCHHH
       end do
#endif           

#ifdef OPENACCHHH
    end do
    !$acc end parallel loop

    !$acc exit data delete(tract, vewal, velfi , velsh , tveno,    &
    !$acc tvelo,tvefi,tvedi,elcod,elvel,bocod,bovel,delta_aux,      &
    !$acc vikin,velfr,yplus,fact1,fact2,                            &
    !$acc avvfr, tven2,produ,gbsur,eucta,                           &
    !$acc  baloc,elmat,elrhs,gbvis,gbden,                           &
    !$acc elmar,                                                    &
    !$acc elmar(pblty)%shape,                                       &
    !$acc elmar(pblty)%deriv,                                       &
    !$acc elmar(pblty)%weigp,                                       &
    !$acc list_boundaries_p,list_boundaries,lelbo,        &
    !$acc lboel,kfl_fixbo_nsi,ywalb,pnode )  


#endif       

  end subroutine nsi_boundary_operations_fast5

  subroutine nsi_boundary_operations_fast_bck1(&
       itask,VECTOR_DIM,pnodb,pgaub,list_boundaries)
!$acc routine(frivel) seq

    use def_parame
    use def_elmtyp
    use def_master
    use def_kermod
    use def_domain
    use def_nastin
    use mod_ker_proper
    use mod_nsi_assembly_global_system, only : nsi_assembly_fractional_step_boundary_scalar

    implicit none
    !
    ! Input and output variables
    !
    integer(ip), intent(in) :: itask                                         !< What to do
    integer(ip), intent(in) :: VECTOR_DIM                                    !< Pack size
    integer(ip), intent(in) :: pnodb                                         !< Number of boundary nodes
    integer(ip), intent(in) :: pgaub                                         !< Number of Gauss points
    integer(ip), intent(in) :: list_boundaries(VECTOR_DIM)                   !< List of boundaries

    real(rp)                :: elmat(VECTOR_DIM,mnode*ndime,mnode*ndime)     ! Element matrices
    real(rp)                :: elrhs(VECTOR_DIM,ndime,mnode)

    real(rp)                :: baloc(VECTOR_DIM,ndime,ndime,pgaub)           ! Gather
    real(rp)                :: bovel(VECTOR_DIM,ndime,pnodb)
    real(rp)                :: bocod(VECTOR_DIM,ndime,pnodb)
    real(rp)                :: elcod(VECTOR_DIM,ndime,mnode)
    real(rp)                :: elvel(VECTOR_DIM,ndime,mnode)

    real(rp)                :: gbsur(VECTOR_DIM,pgaub)
    real(rp)                :: eucta                           ! Value at Gauss points

    real(rp)                :: gbden(VECTOR_DIM,pgaub)
    real(rp)                :: gbvis(VECTOR_DIM,pgaub)
    real(rp)                :: gbmut(VECTOR_DIM,pgaub)

    real(rp)                :: gbden_tmp(pgaub)
    real(rp)                :: gbvis_tmp(pgaub)
    real(rp)                :: gbmut_tmp(pgaub)

    real(rp)                :: baloc_tmp(ndime,ndime)
    real(rp)                :: roughness
    integer(ip)             :: ielem,ipoin,inode,idime,jnode,jdime     ! Indices and dimensions
    integer(ip)             :: iboun,igaub,inodb,kdime,idofn
    integer(ip)             :: pblty,ievab,jevab
    integer(ip)             :: dummi,jnodb

    real(rp)                :: velsh(VECTOR_DIM,ndime,pnodb*ndime)
    real(rp)                :: tvelo(VECTOR_DIM,ndime)
    real(rp)                :: tvefi(VECTOR_DIM,ndime)
    real(rp)                :: tvedi(VECTOR_DIM,ndime)
    real(rp)                :: vikin(VECTOR_DIM)
    real(rp)                :: tveno(VECTOR_DIM)
    real(rp)                :: velfr(VECTOR_DIM)
    real(rp)                :: avvfr(VECTOR_DIM)
    real(rp)                :: tven2(VECTOR_DIM)
    real(rp)                :: yplus(VECTOR_DIM)
    real(rp)                :: fact1(VECTOR_DIM)
    real(rp)                :: fact2(VECTOR_DIM)
    real(rp)                :: velfi(VECTOR_DIM,ndime)
    real(rp)                :: vewal(VECTOR_DIM,ndime)
    real(rp)                :: tract(VECTOR_DIM,3)    
    real(rp)                :: delta_aux(VECTOR_DIM)

    integer(ip)             :: pnode(VECTOR_DIM)
    integer(ip)             :: list_boundaries_p(VECTOR_DIM)                      ! List of elements (always positive)
    integer(ip)             :: ivect

    gbden     = 0.0_rp
    gbvis     = 0.0_rp
    gbmut     = 0.0_rp
    elmat     = 0.0_rp
    elrhs     = 0.0_rp
    roughness = 0.0_rp

    !--------------------------------------------------------------------
    !
    ! List of elements
    !
    ! +----+----+----+----+
    ! | 23 | 24 | 25 | 0  | <= list_boundaries
    ! +----+----+----+----+
    !
    ! +----+----+----+----+
    ! | 23 | 24 | 25 | 23 | <= list_boundaries_p
    ! +----+----+----+----+
    !
    !--------------------------------------------------------------------

    iboun = list_boundaries(1)
    pblty = abs(ltypb(iboun))

    do ivect = 1,VECTOR_DIM                      
       iboun = abs(list_boundaries(ivect))
       if( iboun /= 0 ) then
          list_boundaries_p(ivect) = list_boundaries(ivect)
       else
          list_boundaries_p(ivect) = list_boundaries(1)
       end if
    end do

    if( maxval( kfl_fixbo_nsi(list_boundaries_p(:))) == 0 ) return

    !--------------------------------------------------------------------
    !
    ! Properties
    !
    !--------------------------------------------------------------------

    do ivect = 1,VECTOR_DIM
       iboun = list_boundaries_p(ivect)    
       call ker_proper('DENSI','PGAUB',dummi,iboun,gbden_tmp)
       call ker_proper('VISCO','PGAUB',dummi,iboun,gbvis_tmp)
       call ker_proper('TURBU','PGAUB',dummi,iboun,gbmut_tmp)
       do igaub = 1,pgaub
          gbden(ivect,igaub) = gbden_tmp(igaub)
          gbvis(ivect,igaub) = gbvis_tmp(igaub)
          gbmut(ivect,igaub) = gbmut_tmp(igaub)
       end do
    end do

    !--------------------------------------------------------------------
    !
    ! Gather operations: ELVEL, ELCOD, BOVEL, DELTA_AUX
    !
    !--------------------------------------------------------------------
    
    do ivect = 1,VECTOR_DIM
       iboun = list_boundaries_p(ivect)
       ielem = lelbo(iboun)
       pnode(ivect) = lnnod(ielem)
       do inode = 1,pnode(ivect)
          ipoin = lnods(inode,ielem)
          do idime = 1,ndime
             elcod(ivect,idime,inode) = coord(idime,ipoin)
             elvel(ivect,idime,inode) = veloc(idime,ipoin,1)
          end do
       end do
       do inodb = 1,pnodb
          ipoin = lnodb(inodb,iboun)
          do idime = 1,ndime
             bocod(ivect,idime,inodb) = coord(idime,ipoin)
             bovel(ivect,idime,inodb) = veloc(idime,ipoin,1)
          end do
       end do
       if ( kfl_waexl_ker == 0_ip ) then                            ! normal behaviour
          if( kfl_delta == 1 ) then
             delta_aux(ivect) = ywalb(iboun)                        ! variable wall distance
          else
             delta_aux(ivect) = delta_nsi                           ! fixed wall distance
          end if
       end if
    end do

    !--------------------------------------------------------------------
    !
    ! Shape functions
    !
    !--------------------------------------------------------------------

    do ivect = 1,VECTOR_DIM
       iboun = list_boundaries_p(ivect)
       do igaub = 1,pgaub
          call bouder(&
             pnodb,ndime,ndimb,elmar(pblty)%deriv(:,:,igaub),&    ! Cartesian derivative
             bocod(ivect,:,:),baloc_tmp,eucta)                    ! and Jacobian
          call chenor(&
             pnode(ivect),baloc_tmp,bocod(ivect,:,:),elcod(ivect,:,:))   ! Normalize normal
          gbsur(ivect,igaub)     = elmar(pblty)%weigp(igaub)*eucta
          baloc(ivect,:,:,igaub) = baloc_tmp(:,:)
       end do
    end do

#ifdef OPENACCHHH

    !$acc enter data create(tract, vewal, velfi , velsh , tveno,    &
    !$acc tvelo,tvefi,tvedi,                                  &
    !$acc vikin,velfr,yplus,fact1,fact2,                            &
    !$acc avvfr, tven2)                                             &
    !$acc copyin(                                                   & 
    !$acc elmar,                                                    &
    !$acc elmar(pblty)%shape,                                       &
    !$acc baloc,gbvis, gbden,gbsur,                                 &
    !$acc list_boundaries_p,list_boundaries,delta_aux,lelbo,        &
    !$acc bovel,lboel,elmat,elrhs,elvel,kfl_fixbo_nsi )  

    !
    !$acc parallel loop gang vector default(present)
    !
    do ivect = 1,VECTOR_DIM
#endif       


       !--------------------------------------------------------------------
       !
       ! Element matrix
       !
       !--------------------------------------------------------------------

       do igaub = 1,pgaub

          tract(DEF_VECT,:)   = 0.0_rp
          vewal(DEF_VECT,:)   = 0.0_rp
          velfi(DEF_VECT,:)   = 0.0_rp
          velsh(DEF_VECT,:,:) = 0.0_rp
          tveno(DEF_VECT)     = 0.0_rp
 

          if ( kfl_waexl_ker == 0_ip ) then  ! normal behaviour
             do inodb = 1,pnodb
                do idime = 1,ndime
                   vewal(DEF_VECT,idime) = vewal(DEF_VECT,idime) &
                        + bovel(DEF_VECT,idime,inodb) * elmar(pblty) % shape(inodb,igaub)
                end do
             end do
          end if

          tvelo(DEF_VECT,:)   = vewal(DEF_VECT,:)
          tvefi(DEF_VECT,:)   = velfi(DEF_VECT,:)
          do idime = 1,ndime                                     ! Tangent veloc.
             do inodb = 1,pnodb
                idofn = (inodb-1)*ndime+idime
                velsh(DEF_VECT,idime,idofn) = velsh(DEF_VECT,idime,idofn) &
                     + elmar(pblty) % shape(inodb,igaub)
             end do
             do kdime = 1,ndime
                tvelo(DEF_VECT,idime) = tvelo(DEF_VECT,idime) &
                     - baloc(DEF_VECT,idime,ndime,igaub)   &
                     * baloc(DEF_VECT,kdime,ndime,igaub) * vewal(DEF_VECT,kdime)
                tvefi(DEF_VECT,idime) = tvefi(DEF_VECT,idime)   &
                     - baloc(DEF_VECT,idime,ndime,igaub) &
                     * baloc(DEF_VECT,kdime,ndime,igaub) * velfi(DEF_VECT,kdime)
                do inodb = 1,pnodb
                   idofn = (inodb-1)*ndime+kdime
                   velsh(DEF_VECT,idime,idofn) = velsh(DEF_VECT,idime,idofn) &
                        - baloc(DEF_VECT,idime,ndime,igaub)                  &
                        * baloc(DEF_VECT,kdime,ndime,igaub) * elmar(pblty) % shape(inodb,igaub)
                end do
             end do
             tvedi(DEF_VECT,idime) = tvelo(DEF_VECT,idime) - tvefi(DEF_VECT,idime)
          end do
          !
          ! Compute U*
          ! 
          vikin(DEF_VECT) = gbvis(DEF_VECT,igaub) / gbden(DEF_VECT,igaub)     ! nu

          do idime = 1,ndime
             tveno(DEF_VECT) = tveno(DEF_VECT) + tvedi(DEF_VECT,idime)**2
          end do
          tveno(DEF_VECT) = sqrt(tveno(DEF_VECT))

#ifndef OPENACCHHH
          do ivect = 1,VECTOR_DIM
#endif             
             iboun = list_boundaries_p(ivect)
!             include './my_frivel_include.txt'
             call frivel(kfl_ustar,delta_aux(ivect),roughness,tveno(ivect),vikin(ivect),velfr(ivect)) 
             avvfr(ivect) = velfr(ivect)
             tven2(ivect) = tveno(ivect)
             !
             ! Compute prescribed traction
             !
             yplus(ivect) = delta_aux(ivect) * velfr(ivect) / vikin(ivect)
             if( yplus(ivect) < 5.0_rp  .and. kfl_ustar == 0 ) then
                fact1(ivect) = gbden(ivect,igaub) * vikin(ivect) / delta_aux(ivect)
             else
                fact1(ivect) = gbden(ivect,igaub) * avvfr(ivect) * velfr(ivect) / tveno(ivect)
             end if
             do idime = 1,ndime
                tract(ivect,idime) = tract(ivect,idime) + fact1(ivect) * tvefi(ivect,idime)
             end do

             if( kfl_waexl_ker == 0 ) then
                do inodb = 1,pnodb
                   fact2(ivect) = fact1(ivect) * elmar(pblty) % shape(inodb,igaub)
                   inode = lboel(inodb,iboun)
                   ievab = (inode-1)*ndime
                   do idime = 1,ndime
                      ievab = ievab+1
                      do jnodb = 1,pnodb
                         jnode = lboel(jnodb,iboun)
                         jevab = (jnode-1)*ndime
                         do jdime = 1,ndime                          
                            jevab = jevab+1
                            elmat(ivect,ievab,jevab) = elmat(ivect,ievab,jevab)     &
                                 + fact2(ivect) * velsh(ivect,idime,(jnodb-1)*ndime+jdime) &
                                 * gbsur(ivect,igaub)
                         end do
                      end do
                      elrhs(ivect,idime,inode) = elrhs(ivect,idime,inode)           &
                           + tract(ivect,idime) * elmar(pblty) % shape(inodb,igaub) &
                           * gbsur(ivect,igaub) 
                   end do
                end do
             end if
          end do

#ifndef OPENACCHHH
       end do
#endif
       !--------------------------------------------------------------------
       !
       ! Send matrix to RHS
       !
       !--------------------------------------------------------------------

       do jnode = 1,pnode(ivect)
          do jdime = 1,ndime
             jevab = (jnode-1)*ndime+jdime
             do inode = 1,pnode(ivect)
                do idime = 1,ndime
                   ievab = (inode-1)*ndime+idime
                   elrhs(DEF_VECT,idime,inode) = elrhs(DEF_VECT,idime,inode) &
                        - elmat(DEF_VECT,ievab,jevab) * elvel(DEF_VECT,jdime,jnode) 
                end do
             end do
          end do
       end do

       !--------------------------------------------------------------------
       !
       ! Scatter element matrix to global one 
       !
       !--------------------------------------------------------------------

#ifndef OPENACCHHH
       do ivect = 1,VECTOR_DIM
#endif          
          iboun = list_boundaries(ivect)    
          if( iboun > 0 ) then
             if( kfl_fixbo_nsi(iboun) /= 0 ) then
                ielem = lelbo(iboun)  
                do inode = 1,pnode(ivect)
                   ipoin = lnods(inode,ielem)
                   do idime = 1,ndime
                      ievab = idime + (ipoin-1) * ndime 
                      !$acc atomic update
#ifdef NO_COLORING
                      !$OMP ATOMIC
#endif
                      rhsid(ievab) = rhsid(ievab) + elrhs(ivect,idime,inode)

                   end do
                end do
             end if
          end if
#ifndef OPENACCHHH
       end do
#endif           

#ifdef OPENACCHHH
    end do
    !$acc end parallel loop

    !$acc exit data delete (tract, vewal, velfi , velsh , tveno,    &
    !$acc bovel,tvelo,tvefi,tvedi,                                  &
    !$acc delta_aux,                                          &
    !$acc vikin,velfr,yplus,fact1,fact2,                            &
    !$acc avvfr, tven2,                                             &
    !$acc elmar,                                                    &
    !$acc elmar(pblty)%shape,                                       &
    !$acc baloc,gbvis, gbden,gbsur,                                 &
    !$acc list_boundaries_p,bovel,lboel,elmat,elrhs,elvel )  



#endif       

  end subroutine nsi_boundary_operations_fast_bck1

  subroutine nsi_boundary_operations_fast_bck2(&
       itask,VECTOR_DIM,pnodb,pgaub,list_boundaries)
!$acc routine(frivel) seq

    use def_parame
    use def_elmtyp
    use def_master
    use def_kermod
    use def_domain
    use def_nastin
    use mod_ker_proper
    use mod_nsi_assembly_global_system, only : nsi_assembly_fractional_step_boundary_scalar

    implicit none
    !
    ! Input and output variables
    !
    integer(ip), intent(in) :: itask                                         !< What to do
    integer(ip), intent(in) :: VECTOR_DIM                                    !< Pack size
    integer(ip), intent(in) :: pnodb                                         !< Number of boundary nodes
    integer(ip), intent(in) :: pgaub                                         !< Number of Gauss points
    integer(ip), intent(in) :: list_boundaries(VECTOR_DIM)                   !< List of boundaries

    real(rp)                :: elmat(VECTOR_DIM,mnode*ndime,mnode*ndime)     ! Element matrices
    real(rp)                :: elrhs(VECTOR_DIM,ndime,mnode)

    real(rp)                :: baloc(VECTOR_DIM,ndime,ndime,pgaub)           ! Gather
    real(rp)                :: bovel(VECTOR_DIM,ndime,pnodb)
    real(rp)                :: bocod(VECTOR_DIM,ndime,pnodb)
    real(rp)                :: elcod(VECTOR_DIM,ndime,mnode)
    real(rp)                :: elvel(VECTOR_DIM,ndime,mnode)

    real(rp)                :: gbsur(VECTOR_DIM,pgaub)
    real(rp)                :: eucta(VECTOR_DIM)                     ! Values at Gauss points

    real(rp)                :: gbden(VECTOR_DIM,pgaub)
    real(rp)                :: gbvis(VECTOR_DIM,pgaub)
    real(rp)                :: produ(VECTOR_DIM)
    
    real(rp)                :: gbden_tmp(pgaub)
    real(rp)                :: gbvis_tmp(pgaub)

!    real(rp)                :: baloc_tmp(ndime,ndime)
    real(rp)                :: roughness
    integer(ip)             :: ielem,ipoin,inode,idime,jnode,jdime     ! Indices and dimensions
    integer(ip)             :: iboun,igaub,inodb,kdime,idofn
    integer(ip)             :: pblty,ievab,jevab
    integer(ip)             :: dummi,jnodb

    real(rp)                :: velsh(VECTOR_DIM,ndime,pnodb*ndime)
    real(rp)                :: tvelo(VECTOR_DIM,ndime)
    real(rp)                :: tvefi(VECTOR_DIM,ndime)
    real(rp)                :: tvedi(VECTOR_DIM,ndime)
    real(rp)                :: vikin(VECTOR_DIM)
    real(rp)                :: tveno(VECTOR_DIM)
    real(rp)                :: velfr(VECTOR_DIM)
    real(rp)                :: avvfr(VECTOR_DIM)
    real(rp)                :: tven2(VECTOR_DIM)
    real(rp)                :: yplus(VECTOR_DIM)
    real(rp)                :: fact1(VECTOR_DIM)
    real(rp)                :: fact2(VECTOR_DIM)
    real(rp)                :: velfi(VECTOR_DIM,ndime)
    real(rp)                :: vewal(VECTOR_DIM,ndime)
    real(rp)                :: tract(VECTOR_DIM,3)    
    real(rp)                :: delta_aux(VECTOR_DIM)

    integer(ip)             :: pnode(VECTOR_DIM)
    integer(ip)             :: list_boundaries_p(VECTOR_DIM)                      ! List of elements (always positive)
    integer(ip)             :: ivect,ii,jj,kk

    gbden     = 0.0_rp
    gbvis     = 0.0_rp
    elmat     = 0.0_rp
    elrhs     = 0.0_rp
    roughness = 0.0_rp
    baloc     = 0.0_rp
    
    !--------------------------------------------------------------------
    !
    ! List of elements
    !
    ! +----+----+----+----+
    ! | 23 | 24 | 25 | 0  | <= list_boundaries
    ! +----+----+----+----+
    !
    ! +----+----+----+----+
    ! | 23 | 24 | 25 | 23 | <= list_boundaries_p
    ! +----+----+----+----+
    !
    !--------------------------------------------------------------------

    iboun = list_boundaries(1)
    pblty = abs(ltypb(iboun))

    do ivect = 1,VECTOR_DIM                      
       iboun = abs(list_boundaries(ivect))
       if( iboun /= 0 ) then
          list_boundaries_p(ivect) = list_boundaries(ivect)
       else
          list_boundaries_p(ivect) = list_boundaries(1)
       end if
    end do

    if( maxval( kfl_fixbo_nsi(list_boundaries_p(:))) == 0 ) return

    !--------------------------------------------------------------------
    !
    ! Properties
    !
    !--------------------------------------------------------------------

    do ivect = 1,VECTOR_DIM
       iboun = list_boundaries_p(ivect)    
       call ker_proper('DENSI','PGAUB',dummi,iboun,gbden_tmp)
       call ker_proper('VISCO','PGAUB',dummi,iboun,gbvis_tmp)
       do igaub = 1,pgaub
          gbden(ivect,igaub) = gbden_tmp(igaub)
          gbvis(ivect,igaub) = gbvis_tmp(igaub)
       end do
    end do

    !--------------------------------------------------------------------
    !
    ! Gather operations: ELVEL, ELCOD, BOVEL, DELTA_AUX
    !
    !--------------------------------------------------------------------
    
    do ivect = 1,VECTOR_DIM
       iboun = list_boundaries_p(ivect)
       ielem = lelbo(iboun)
       pnode(ivect) = lnnod(ielem)
       do inode = 1,pnode(ivect)
          ipoin = lnods(inode,ielem)
          do idime = 1,ndime
             elcod(ivect,idime,inode) = coord(idime,ipoin)
             elvel(ivect,idime,inode) = veloc(idime,ipoin,1)
          end do
       end do
       do inodb = 1,pnodb
          ipoin = lnodb(inodb,iboun)
          do idime = 1,ndime
             bocod(ivect,idime,inodb) = coord(idime,ipoin)
             bovel(ivect,idime,inodb) = veloc(idime,ipoin,1)
          end do
       end do
       if ( kfl_waexl_ker == 0_ip ) then                            ! normal behaviour
          if( kfl_delta == 1 ) then
             delta_aux(ivect) = ywalb(iboun)                        ! variable wall distance
          else
             delta_aux(ivect) = delta_nsi                           ! fixed wall distance
          end if
       end if
    end do

    !--------------------------------------------------------------------
    !
    ! Shape functions
    !
    !--------------------------------------------------------------------
    !do ivect = 1,VECTOR_DIM
    !   iboun = list_boundaries_p(ivect)
    !   do igaub = 1,pgaub
    !      call bouder(&
    !         pnodb,ndime,ndimb,elmar(pblty)%deriv(1,1,igaub),&    ! Cartesian derivative
    !         bocod(ivect,:,:),baloc_tmp,eucta)                    ! and Jacobian
    !      call chenor(&
    !         pnode,baloc_tmp,bocod(ivect,:,:),elcod(ivect,:,:))   ! Normalize normal
    !      gbsur(ivect,igaub)     = elmar(pblty)%weigp(igaub)*eucta
    !      baloc(ivect,:,:,igaub) = baloc_tmp(:,:)
    !   end do
    !end do

#ifdef OPENACCHHH

    !$acc enter data create(tract, vewal, velfi , velsh , tveno,    &
    !$acc tvelo,tvefi,tvedi,                                  &
    !$acc vikin,velfr,yplus,fact1,fact2,                            &
    !$acc avvfr, tven2,produ,gbsur,eucta)                           &
    !$acc copyin(                                                   & 
    !$acc elmar,                                                    &
    !$acc elmar(pblty)%shape,                                       &
    !$acc elmar(pblty)%deriv,                                       &
    !$acc elmar(pblty)%weigp,                                       &
    !$acc gbvis, gbden,                                 &
    !$acc list_boundaries_p,list_boundaries,delta_aux,lelbo,        &
    !$acc bovel,lboel,elmat,elrhs,elvel,kfl_fixbo_nsi,bocod,baloc )  

    !
    !$acc parallel loop gang vector default(present)
    !
    do ivect = 1,VECTOR_DIM
#endif       


       !--------------------------------------------------------------------
       !
       ! Element matrix
       !
       !--------------------------------------------------------------------

       do igaub = 1,pgaub
          !
          ! Shape functions
          !
          do ii=1,ndime
             do jj=1,ndimb
                do kk=1,pnodb
                   baloc(DEF_VECT,ii,jj,igaub)=baloc(DEF_VECT,ii,jj,igaub)&
                        +bocod(DEF_VECT,ii,kk)*elmar(pblty)%deriv(jj,kk,igaub)
                end do
             end do
          end do
          baloc(DEF_VECT,1,3,igaub) =   baloc(DEF_VECT,2,1,igaub) * baloc(DEF_VECT,3,2,igaub) &
               &                      - baloc(DEF_VECT,3,1,igaub) * baloc(DEF_VECT,2,2,igaub)
          baloc(DEF_VECT,2,3,igaub) =   baloc(DEF_VECT,3,1,igaub) * baloc(DEF_VECT,1,2,igaub) &
               &                      - baloc(DEF_VECT,1,1,igaub) * baloc(DEF_VECT,3,2,igaub)
          baloc(DEF_VECT,3,3,igaub) =   baloc(DEF_VECT,1,1,igaub) * baloc(DEF_VECT,2,2,igaub) &
               &                      - baloc(DEF_VECT,2,1,igaub) * baloc(DEF_VECT,1,2,igaub)
          eucta(DEF_VECT)           =   baloc(DEF_VECT,1,3,igaub) * baloc(DEF_VECT,1,3,igaub) &
               &                      + baloc(DEF_VECT,2,3,igaub) * baloc(DEF_VECT,2,3,igaub) &
               &                      + baloc(DEF_VECT,3,3,igaub) * baloc(DEF_VECT,3,3,igaub)
          eucta(DEF_VECT)           =   sqrt(eucta(DEF_VECT)) 
          baloc(DEF_VECT,1,1,igaub) =   baloc(DEF_VECT,2,3,igaub) * baloc(DEF_VECT,3,2,igaub) &
               &                      - baloc(DEF_VECT,3,3,igaub) * baloc(DEF_VECT,2,2,igaub)
          baloc(DEF_VECT,2,1,igaub) =   baloc(DEF_VECT,3,3,igaub) * baloc(DEF_VECT,1,2,igaub) &
               &                      - baloc(DEF_VECT,1,3,igaub) * baloc(DEF_VECT,3,2,igaub)
          baloc(DEF_VECT,3,1,igaub) =   baloc(DEF_VECT,1,3,igaub) * baloc(DEF_VECT,2,2,igaub) &
               &                      - baloc(DEF_VECT,2,3,igaub) * baloc(DEF_VECT,1,2,igaub)
          do ii = 1,ndime
             produ(DEF_VECT) = 0.0_rp
             do jj = 1,ndime
                produ(DEF_VECT) = produ(DEF_VECT) + baloc(DEF_VECT,jj,ii,igaub)*baloc(DEF_VECT,jj,ii,igaub)
             end do
             produ(DEF_VECT) = sqrt(produ(DEF_VECT))
             do jj = 1,ndime
                baloc(DEF_VECT,jj,ii,igaub) = produ(DEF_VECT) * baloc(DEF_VECT,jj,ii,igaub)
             end do
          end do
          gbsur(DEF_VECT,igaub) = elmar(pblty)%weigp(igaub)*eucta(DEF_VECT)



          !
          ! Traction due to law of the wall
          !
          tract(DEF_VECT,:)   = 0.0_rp
          vewal(DEF_VECT,:)   = 0.0_rp
          velfi(DEF_VECT,:)   = 0.0_rp
          velsh(DEF_VECT,:,:) = 0.0_rp
          tveno(DEF_VECT)     = 0.0_rp
 

          if ( kfl_waexl_ker == 0_ip ) then  ! normal behaviour
             do inodb = 1,pnodb
                do idime = 1,ndime
                   vewal(DEF_VECT,idime) = vewal(DEF_VECT,idime) &
                        + bovel(DEF_VECT,idime,inodb) * elmar(pblty) % shape(inodb,igaub)
                end do
             end do
          end if

          tvelo(DEF_VECT,:)   = vewal(DEF_VECT,:)
          tvefi(DEF_VECT,:)   = velfi(DEF_VECT,:)
          do idime = 1,ndime                                     ! Tangent veloc.
             do inodb = 1,pnodb
                idofn = (inodb-1)*ndime+idime
                velsh(DEF_VECT,idime,idofn) = velsh(DEF_VECT,idime,idofn) &
                     + elmar(pblty) % shape(inodb,igaub)
             end do
             do kdime = 1,ndime
                tvelo(DEF_VECT,idime) = tvelo(DEF_VECT,idime) &
                     - baloc(DEF_VECT,idime,ndime,igaub)   &
                     * baloc(DEF_VECT,kdime,ndime,igaub) * vewal(DEF_VECT,kdime)
                tvefi(DEF_VECT,idime) = tvefi(DEF_VECT,idime)   &
                     - baloc(DEF_VECT,idime,ndime,igaub) &
                     * baloc(DEF_VECT,kdime,ndime,igaub) * velfi(DEF_VECT,kdime)
                do inodb = 1,pnodb
                   idofn = (inodb-1)*ndime+kdime
                   velsh(DEF_VECT,idime,idofn) = velsh(DEF_VECT,idime,idofn) &
                        - baloc(DEF_VECT,idime,ndime,igaub)                  &
                        * baloc(DEF_VECT,kdime,ndime,igaub) * elmar(pblty) % shape(inodb,igaub)
                end do
             end do
             tvedi(DEF_VECT,idime) = tvelo(DEF_VECT,idime) - tvefi(DEF_VECT,idime)
          end do
          !
          ! Compute U*
          ! 
          vikin(DEF_VECT) = gbvis(DEF_VECT,igaub) / gbden(DEF_VECT,igaub)     ! nu

          do idime = 1,ndime
             tveno(DEF_VECT) = tveno(DEF_VECT) + tvedi(DEF_VECT,idime)**2
          end do
          tveno(DEF_VECT) = sqrt(tveno(DEF_VECT))

#ifndef OPENACCHHH
          do ivect = 1,VECTOR_DIM
#endif             
             iboun = list_boundaries_p(ivect)
!             include './my_frivel_include.txt'
             call frivel(kfl_ustar,delta_aux(ivect),roughness,tveno(ivect),vikin(ivect),velfr(ivect)) 
             avvfr(ivect) = velfr(ivect)
             tven2(ivect) = tveno(ivect)
             !
             ! Compute prescribed traction
             !
             yplus(ivect) = delta_aux(ivect) * velfr(ivect) / vikin(ivect)
             if( yplus(ivect) < 5.0_rp  .and. kfl_ustar == 0 ) then
                fact1(ivect) = gbden(ivect,igaub) * vikin(ivect) / delta_aux(ivect)
             else
                fact1(ivect) = gbden(ivect,igaub) * avvfr(ivect) * velfr(ivect) / tveno(ivect)
             end if
             do idime = 1,ndime
                tract(ivect,idime) = tract(ivect,idime) + fact1(ivect) * tvefi(ivect,idime)
             end do

             if( kfl_waexl_ker == 0 ) then
                do inodb = 1,pnodb
                   fact2(ivect) = fact1(ivect) * elmar(pblty) % shape(inodb,igaub)
                   inode = lboel(inodb,iboun)
                   ievab = (inode-1)*ndime
                   do idime = 1,ndime
                      ievab = ievab+1
                      do jnodb = 1,pnodb
                         jnode = lboel(jnodb,iboun)
                         jevab = (jnode-1)*ndime
                         do jdime = 1,ndime                          
                            jevab = jevab+1
                            elmat(ivect,ievab,jevab) = elmat(ivect,ievab,jevab)     &
                                 + fact2(ivect) * velsh(ivect,idime,(jnodb-1)*ndime+jdime) &
                                 * gbsur(ivect,igaub)
                         end do
                      end do
                      elrhs(ivect,idime,inode) = elrhs(ivect,idime,inode)           &
                           + tract(ivect,idime) * elmar(pblty) % shape(inodb,igaub) &
                           * gbsur(ivect,igaub) 
                   end do
                end do
             end if
          end do

#ifndef OPENACCHHH
       end do
#endif
       !--------------------------------------------------------------------
       !
       ! Send matrix to RHS
       !
       !--------------------------------------------------------------------

       do jnode = 1,pnode(ivect)
          do jdime = 1,ndime
             jevab = (jnode-1)*ndime+jdime
             do inode = 1,pnode(ivect)
                do idime = 1,ndime
                   ievab = (inode-1)*ndime+idime
                   elrhs(DEF_VECT,idime,inode) = elrhs(DEF_VECT,idime,inode) &
                        - elmat(DEF_VECT,ievab,jevab) * elvel(DEF_VECT,jdime,jnode) 
                end do
             end do
          end do
       end do

       !--------------------------------------------------------------------
       !
       ! Scatter element matrix to global one 
       !
       !--------------------------------------------------------------------

#ifndef OPENACCHHH
       do ivect = 1,VECTOR_DIM
#endif          
          iboun = list_boundaries(ivect)    
          if( iboun > 0 ) then
             if( kfl_fixbo_nsi(iboun) /= 0 ) then
                ielem = lelbo(iboun)  
                do inode = 1,pnode(ivect)
                   ipoin = lnods(inode,ielem)
                   do idime = 1,ndime
                      ievab = idime + (ipoin-1) * ndime 
                      !$acc atomic update
#ifdef NO_COLORING
                      !$OMP ATOMIC
#endif
                      rhsid(ievab) = rhsid(ievab) + elrhs(ivect,idime,inode)

                   end do
                end do
             end if
          end if
#ifndef OPENACCHHH
       end do
#endif           

#ifdef OPENACCHHH
    end do
    !$acc end parallel loop

    !$acc exit data delete(tract, vewal, velfi , velsh , tveno,    &
    !$acc tvelo,tvefi,tvedi,                                  &
    !$acc vikin,velfr,yplus,fact1,fact2,                            &
    !$acc avvfr, tven2,baloc,produ,gbsur,                           &
    !$acc elmar,                                                    &
    !$acc elmar(pblty)%shape,                                       &
    !$acc elmar(pblty)%deriv,                                       &
    !$acc elmar(pblty)%weigp,                                       &
    !$acc gbvis, gbden,eucta,                                 &
    !$acc list_boundaries_p,list_boundaries,delta_aux,lelbo,        &
    !$acc bovel,lboel,elmat,elrhs,elvel,kfl_fixbo_nsi,bocod )  



#endif       

  end subroutine nsi_boundary_operations_fast_bck2
  
end module mod_nsi_boundary_operations_fast
!> @}
