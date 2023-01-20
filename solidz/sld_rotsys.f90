!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_rotsys.f90
!> @date    25/04/2016
!> @author  Mariano Vazquez
!> @brief   Rotate elmat, elrhs
!> @details Rotate elmat, elrhs
!>
!> This routine modifies a matrix the components of which are
!> submatrices by rotating them using a matrix R ( = jacrot_du_dq) and its inverse T (=jacrot_dq_du) as follows:
!>
!> R: Local  -> Global
!> T: Global -> Local
!>
!>
!v                |           A_11  ....          A_1i  DU_DQ   ....          A_1n   |
!>                |   ............................................................   |
!>     ELMAT <--  |   DQ_DU_t A_i1  ....  DQ_DU_t A_ii  DU_DQ   ....  DQ_DU_t A_in   |
!>                |   ............................................................   |
!>                |           A_n1  ....          A_ni  DU_DQ   ....          A_nn   |
!>
!> where i ( = IMODI) is a given position. Also, the i-subvector of the
!> given vector BVECT is multiplied by T.
!>
!> The rotation matrix can be any kind of jacobian matrix, providing that
!> the inverse matrix is computed elsewhere.
!>
!> Let's suppose A u = b is the original system. Add R: q -> u and T: u -> q
!> its inverse. R is jacrot_du_dq and T is jacrot_dq_du.
!> Then,
!>         T   A        u  =  T b
!>         T   A (RT)   u  =  T b
!>        (T   A  R ) T u  =  T b
!>        (T   A  R )   q  =  T b
!>
!> the subroutine computes TAR and Tb.
!> Once this new system is solved, this subru also computes u = Rq
!>
!>
!> imodi: target node
!>
!>
!>
!> @}
!------------------------------------------------------------------------

subroutine sld_rotsys(&
     itask,imodi,pnode,pdofn,pevat,elmat,elrhs,jacrot_du_dq,jacrot_dq_du)

  use def_kintyp, only : ip,rp
  
  implicit none

  integer(ip), intent(in)    :: pnode             !> pnode
  integer(ip), intent(in)    :: pevat             !> total dof*pnode
  integer(ip), intent(in)    :: pdofn             !> total dof
  integer(ip), intent(in)    :: itask             !> rotate or rotate back?
  real(rp),    intent(in)    :: jacrot_du_dq(pdofn,pdofn)
  real(rp),    intent(in)    :: jacrot_dq_du(pdofn,pdofn)
  real(rp),    intent(inout) :: elmat(pevat,pevat)
  real(rp),    intent(inout) :: elrhs(pevat)
  
  real(rp)                   :: worma(pdofn,pdofn)
  integer(ip)                :: imodi,itot0,jtot0,inode,jnode
  integer(ip)                :: idofn,jdofn,kdofn,itott,jtott,ktott

  select case (itask)

  case(1)
     !
     ! Modifies column number IMODI of ELMAT ( A_j,imodi <-- A_j,imodi jacrot_du_dq )
     !
     !
     ! RHS and ELMAT: (Global --> Local)
     !
     jtot0=(imodi-1)*pdofn
     do inode=1,pnode
        itot0=(inode-1)*pdofn
        do idofn=1,pdofn
           itott=itot0+idofn
           do jdofn=1,pdofn
              worma(idofn,jdofn)=0.0_rp
              do kdofn=1,pdofn
                 ktott=jtot0+kdofn
                 worma(idofn,jdofn)=worma(idofn,jdofn)&
                      +elmat(itott,ktott)*jacrot_du_dq(kdofn,jdofn)
              end do
           end do
        end do
        do idofn=1,pdofn
           itott=itot0+idofn
           do jdofn=1,pdofn
              jtott=jtot0+jdofn
              elmat(itott,jtott)=worma(idofn,jdofn)
           end do
        end do
     end do
     !
     ! Modifies row number IMODI of ELMAT ( A_imodi,j <-- jacrot_dq_du A_imodi,j )
     !
     itot0=(imodi-1)*pdofn
     do jnode=1,pnode
        jtot0=(jnode-1)*pdofn
        do idofn=1,pdofn
           do jdofn=1,pdofn
              jtott=jtot0+jdofn
              worma(idofn,jdofn)=0.0_rp
              do kdofn=1,pdofn
                 ktott=itot0+kdofn
                 worma(idofn,jdofn)=worma(idofn,jdofn)&
                      +jacrot_dq_du(idofn,kdofn)*elmat(ktott,jtott)
              end do
           end do
        end do
        do idofn=1,pdofn
           itott=itot0+idofn
           do jdofn=1,pdofn
              jtott=jtot0+jdofn
              elmat(itott,jtott)=worma(idofn,jdofn)
           end do
        end do
     end do
     !
     ! Modifies the IMODI subvector of the vector ELRHS
     !
     itot0=(imodi-1)*pdofn
     do idofn=1,pdofn
        worma(idofn,1)=0.0_rp
        do kdofn=1,pdofn
           ktott=itot0+kdofn
           worma(idofn,1)=worma(idofn,1)&
                +jacrot_dq_du(idofn,kdofn)*elrhs(ktott)
        end do
     end do
     do idofn=1,pdofn
        itott=itot0+idofn
        elrhs(itott)=worma(idofn,1)
     end do

  case(2)
     !
     ! Rotate RHS: (Local --> Global)
     !
     itot0=(imodi-1)*pdofn
     do idofn=1,pdofn
        worma(idofn,1)=0.0_rp
        do kdofn=1,pdofn
           ktott=itot0+kdofn
           worma(idofn,1)=worma(idofn,1)&
                +jacrot_du_dq(idofn,kdofn)*elrhs(ktott)
        end do
     end do
     do idofn=1,pdofn
        itott=itot0+idofn
        elrhs(itott)=worma(idofn,1)
     end do

  case(3)
     !
     ! Rotate RHS: (Global --> Local)
     !
     itot0=(imodi-1)*pdofn
     do idofn=1,pdofn
        worma(idofn,1)=0.0_rp
        do kdofn=1,pdofn
           ktott=itot0+kdofn
           worma(idofn,1)=worma(idofn,1)&
                +jacrot_dq_du(idofn,kdofn)*elrhs(ktott)
        end do
     end do
     do idofn=1,pdofn
        itott=itot0+idofn
        elrhs(itott)=worma(idofn,1)
     end do

  end select


end subroutine sld_rotsys
