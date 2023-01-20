!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup NastinBubble
!> @ingroup    Nastin
!> @{
!> @file    mod_nsi_bubble.f90
!> @author  Guillaume Houzeaux
!> @brief   Bubble things
!> @details Bubble things
!>
!------------------------------------------------------------------------

module mod_nsi_bubble

#include "def_vector_size.inc"
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,nelem,lnods,lnnod
  use def_master, only       :  press,veloc
  use mod_maths,  only       :  maths_schur_complement
  use def_nastin, only       :  kfl_bubbl_nsi
  use def_nastin, only       :  bubble_nsi
  use def_nastin, only       :  bubble_Aqq_nsi
  use def_nastin, only       :  bubble_Aqu_nsi
  use def_nastin, only       :  bubble_Aqp_nsi
  use def_nastin, only       :  bubble_bq_nsi
  implicit none
  private

  real(rp), parameter :: zeror = epsilon(1.0_rp)
 
  interface nsi_bubble_assembly
     module procedure nsi_bubble_assembly_scalar,&
          &           nsi_bubble_assembly_vector
  end interface nsi_bubble_assembly


  public :: nsi_bubble_assembly
  public :: nsi_bubble_update
  public :: nsi_eliminate_bubble

contains

  !------------------------------------------------------------------------
  !> @addtogroup NastinMatrixAssembly
  !> @{
  !> @file    nsi_assemble_schur.f90
  !> @author  Guillaume Houzeaux
  !> @brief   Bubble
  !> @details Bubble Assembly:
  !>
  !>          Aqq <= Aqq^(e)
  !>          Aqu <= Aqu^(e)
  !>          Aqp <= Aqp^(e)
  !>          bq  <= bq^(e)
  !> @} 
  !------------------------------------------------------------------------

  subroutine nsi_bubble_assembly_scalar(&
       pnode,ielem,elaqq,elaqu,elaqp,elrbq)

    integer(ip), intent(in)    :: pnode
    integer(ip), intent(in)    :: ielem
    real(rp),    intent(in)    :: elaqq(1,1)
    real(rp),    intent(in)    :: elaqu(1,pnode*ndime)
    real(rp),    intent(in)    :: elaqp(1,pnode)
    real(rp),    intent(in)    :: elrbq(1)
    integer(ip)                :: inode,idime,idofv

    bubble_Aqq_nsi(ielem) = elaqq(1,1)
    bubble_bq_nsi(ielem)  = elrbq(1)

    do inode = 1,pnode
       do idime = 1,ndime
          idofv = (inode-1)*ndime+idime
          bubble_Aqu_nsi(idofv,ielem) = elaqu(1,idofv)
       end do
       bubble_Aqp_nsi(inode,ielem) = elaqp(1,inode)
    end do

  end subroutine nsi_bubble_assembly_scalar

  subroutine nsi_bubble_assembly_vector(&
       pnode,list_elements,elaqq,elaqu,elaqp,elrbq)

    integer(ip), intent(in)    :: pnode
    integer(ip), intent(in)    :: list_elements(VECTOR_SIZE)
    real(rp),    intent(in)    :: elaqq(VECTOR_SIZE,1,1)
    real(rp),    intent(in)    :: elaqu(VECTOR_SIZE,1,pnode*ndime)
    real(rp),    intent(in)    :: elaqp(VECTOR_SIZE,1,pnode)
    real(rp),    intent(in)    :: elrbq(VECTOR_SIZE,1)
    integer(ip)                :: inode,idime,idofv
    integer(ip)                :: ivect,ielem

    do ivect = 1,VECTOR_SIZE
       ielem = list_elements(ivect)
       if( ielem > 0 ) then    
          bubble_Aqq_nsi(ielem) = elaqq(ivect,1,1)
          bubble_bq_nsi(ielem)  = elrbq(ivect,1)

          do inode = 1,pnode
             do idime = 1,ndime
                idofv = (inode-1)*ndime+idime
                bubble_Aqp_nsi(inode,ielem) = elaqp(ivect,1,inode)
                bubble_Aqu_nsi(idofv,ielem) = elaqu(ivect,1,idofv)
             end do
          end do
       end if
    end do

  end subroutine nsi_bubble_assembly_vector

  !------------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   Bubble
  !> @details Bubble
  !>
  !------------------------------------------------------------------------

  subroutine nsi_bubble_update()

    integer(ip) :: inode,idime,idofv,ipoin,ielem


    if( kfl_bubbl_nsi /= 0 ) then

       do ielem = 1,nelem 
          bubble_nsi(ielem) = bubble_bq_nsi(ielem)

          do inode = 1,lnnod(ielem)
             ipoin = lnods(inode,ielem)
             do idime = 1,ndime
                idofv = (inode-1)*ndime+idime
                bubble_nsi(ielem) = bubble_nsi(ielem) - bubble_Aqu_nsi(idofv,ielem) * veloc(idime,ipoin,1)
             end do
             bubble_nsi(ielem) = bubble_nsi(ielem) - bubble_Aqp_nsi(inode,ielem) * press(ipoin,1)
          end do
          if( bubble_Aqq_nsi(ielem) /= 0.0_rp ) then
             bubble_nsi(ielem) = bubble_nsi(ielem) / bubble_Aqq_nsi(ielem)
          else
             bubble_nsi(ielem) = 0.0_rp
          end if
       end do

    end if

  end subroutine nsi_bubble_update

  !------------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   Eliminate bubble
  !> @details Elimiate bubble at element level
  !>          Auu <= Auu - Auq Aqq^-1 Aqu
  !>          Aup <= Aup - Auq Aqq^-1 Aqp
  !>          bu  <= bu  - Auq Aqq^-1 bq
  !>
  !>          Apu <= Apu - Apq Aqq^-1 Aqu
  !>          App <= App - Apq Aqq^-1 Aqp
  !>          bp  <= bp  - Apq Aqq^-1 bq
  !> 
  !>
  !------------------------------------------------------------------------

  subroutine nsi_eliminate_bubble(&
       pnode,pevat,pbubl,elauu,elaup,elapu,elapp,elrbu,elrbp,&
       elauq,elapq,elaqu,elaqp,elaqq,elrbq)

    use mod_nsi_element_assembly, only : nsi_element_system_output
    ! Sizes
    integer(ip), intent(in)    :: pnode
    integer(ip), intent(in)    :: pevat
    integer(ip), intent(in)    :: pbubl(VECTOR_SIZE)
    ! Element matrices
    real(rp),    intent(out)   :: elauu(VECTOR_SIZE,pnode*ndime,pnode*ndime)
    real(rp),    intent(out)   :: elaup(VECTOR_SIZE,pnode*ndime,pnode)
    real(rp),    intent(out)   :: elapp(VECTOR_SIZE,pnode,pnode)
    real(rp),    intent(out)   :: elapu(VECTOR_SIZE,pnode,pnode*ndime)
    real(rp),    intent(out)   :: elrbu(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(out)   :: elrbp(VECTOR_SIZE,pnode)
    ! Enrichement Element matrices
    real(rp),    intent(out)   :: elauq(VECTOR_SIZE,pnode*ndime,1)
    real(rp),    intent(out)   :: elapq(VECTOR_SIZE,pnode,1)
    real(rp),    intent(out)   :: elaqu(VECTOR_SIZE,1,pnode*ndime)
    real(rp),    intent(out)   :: elaqp(VECTOR_SIZE,1,pnode)
    real(rp),    intent(out)   :: elaqq(VECTOR_SIZE,1,1)
    real(rp),    intent(out)   :: elrbq(VECTOR_SIZE,1)

    integer(ip)                :: ivect


    do ivect = 1,VECTOR_SIZE
       if( pbubl(ivect) == 1 ) then
          if( abs(elapp(ivect,1,1)) < 1.0e-16_rp ) call runend('CANNOT ELIMINATE BUBBLE: NULL DIAGONAL TERM') 
          call maths_schur_complement(pevat,pevat, 1_ip,elauu(ivect,:,:),elauq(ivect,:,:),elaqq(ivect,:,:),elaqu(ivect,:,:))
          call maths_schur_complement(pevat,pnode, 1_ip,elaup(ivect,:,:),elauq(ivect,:,:),elaqq(ivect,:,:),elaqp(ivect,:,:))
          call maths_schur_complement(pevat, 1_ip, 1_ip,elrbu(ivect,:,:),elauq(ivect,:,:),elaqq(ivect,:,:),elrbq(ivect,:))
          
          call maths_schur_complement(pnode,pevat, 1_ip,elapu(ivect,:,:),elapq(ivect,:,:),elaqq(ivect,:,:),elaqu(ivect,:,:))
          call maths_schur_complement(pnode,pnode, 1_ip,elapp(ivect,:,:),elapq(ivect,:,:),elaqq(ivect,:,:),elaqp(ivect,:,:))
          call maths_schur_complement(pnode, 1_ip, 1_ip,elrbp(ivect,:)  ,elapq(ivect,:,:),elaqq(ivect,:,:),elrbq(ivect,:))
       end if
    end do


  end subroutine nsi_eliminate_bubble

end module mod_nsi_bubble
!> @}
