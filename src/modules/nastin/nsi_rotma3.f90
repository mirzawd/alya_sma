!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_rotma3(&
     imodi,pnode,pevat,elauu,elaup,elapu,elrbu,elcmm,rotma,worma)
  !------------------------------------------------------------------------
  !
  ! This routine modifies a matrix the components of which are
  ! submatrices by rotating them using a matrix R ( = ROTMA) as follows:
  !
  !                |       A_11  ....      A_1i  R   ....      A_1n   |
  !                |   ............................................   |
  !     ELMAT <--  |   R^t A_i1  ....  R^t A_ii  R   ....  R^t A_in   |
  !                |   ............................................   |
  !                |       A_n1  ....      A_ni  R   ....      A_nn   |
  !
  ! where i ( = IMODI) is a given position. Also, the i-subvector of the
  ! given vector BVECT is multiplied by R^t. 
  !
  !------------------------------------------------------------------------
  use def_parame
  use def_domain, only       :  ndime
  use def_nastin, only       :  kfl_corre_nsi
  implicit none
  integer(ip), intent(in)    :: pnode,pevat
  real(rp),    intent(in)    :: rotma(ndime,ndime)
  real(rp),    intent(inout) :: elauu(pevat,pevat)
  real(rp),    intent(inout) :: elaup(pevat,pnode)
  real(rp),    intent(inout) :: elapu(pnode,pevat)
  real(rp),    intent(inout) :: elrbu(pevat)
  real(rp),    intent(inout) :: elcmm(pevat,pevat)
  real(rp),    intent(inout) :: worma(ndime,ndime)
  integer(ip)                :: imodi,itot0,jtot0,inode,jnode,itotp,jtotp
  integer(ip)                :: idime,jdime,kdime,itotv,jtotv,ktotv
  !
  ! Modifies column number IMODI of ELMAT ( A_j,imodi <-- A_j,imodi R )
  !
  jtot0=(imodi-1)*ndime
  do inode=1,pnode
     itot0=(inode-1)*ndime
     do idime=1,ndime
        itotv=itot0+idime
        do jdime=1,ndime
           worma(idime,jdime)=0.0_rp
           do kdime=1,ndime
              ktotv=jtot0+kdime
              worma(idime,jdime)=worma(idime,jdime)&
                   +elauu(itotv,ktotv)*rotma(kdime,jdime)
           end do
        end do
     end do
     do idime=1,ndime
        itotv=itot0+idime
        do jdime=1,ndime
           jtotv=jtot0+jdime
           elauu(itotv,jtotv)=worma(idime,jdime)
        end do
     end do
  end do
  !
  ! Modify the part corresponding to a scalar unknown
  !      
  do itotv=1,pnode
     itotp=itotv
     do jdime=1,ndime
        worma(1,jdime)=0.0_rp
        do kdime=1,ndime
           ktotv=jtot0+kdime
           worma(1,jdime)=worma(1,jdime)&
                +elapu(itotp,ktotv)*rotma(kdime,jdime)
        end do
     end do
     do jdime=1,ndime
        jtotv=jtot0+jdime
        elapu(itotp,jtotv)=worma(1,jdime)
     end do
  end do
  !
  ! Modifies row number IMODI of ELMAT ( A_imodi,j <-- R^t A_imodi,j )
  !
  itot0=(imodi-1)*ndime
  do jnode=1,pnode
     jtot0=(jnode-1)*ndime
     do idime=1,ndime
        do jdime=1,ndime
           jtotv=jtot0+jdime
           worma(idime,jdime)=0.0_rp
           do kdime=1,ndime
              ktotv=itot0+kdime
              worma(idime,jdime)=worma(idime,jdime)&
                   +elauu(ktotv,jtotv)*rotma(kdime,idime)
           end do
        end do
     end do
     do idime=1,ndime
        itotv=itot0+idime
        do jdime=1,ndime
           jtotv=jtot0+jdime
           elauu(itotv,jtotv)=worma(idime,jdime)
        end do
     end do
  end do
  !
  ! Modify the part corresponding to a scalar unknown
  !     
  do jtotv=1,pnode
     jtotp=jtotv
     do idime=1,ndime
        worma(idime,1)=0.0_rp
        do kdime=1,ndime
           ktotv=itot0+kdime
           worma(idime,1)=worma(idime,1)&
                +rotma(kdime,idime)*elaup(ktotv,jtotp)
        end do
     end do
     do idime=1,ndime
        itotv=itot0+idime
        elaup(itotv,jtotp)=worma(idime,1)
     end do
  end do
  !
  ! Modifies the IMODI subvector of the vector BVECT
  !
  itot0=(imodi-1)*ndime
  do idime=1,ndime
     worma(idime,1)=0.0_rp
     do kdime=1,ndime
        ktotv=itot0+kdime
        worma(idime,1)=worma(idime,1)&
             +rotma(kdime,idime)*elrbu(ktotv)
     end do
  end do
  do idime=1,ndime
     itotv=itot0+idime
     elrbu(itotv)=worma(idime,1)
  end do

  !
  ! Modifies column number IMODI of ELCMM ( A_j,imodi <-- A_j,imodi R )
  !
  if( kfl_corre_nsi == 3 ) then
     jtot0=(imodi-1)*ndime
     do inode=1,pnode
        itot0=(inode-1)*ndime
        do idime=1,ndime
           itotv=itot0+idime
           do jdime=1,ndime
              worma(idime,jdime)=0.0_rp
              do kdime=1,ndime
                 ktotv=jtot0+kdime
                 worma(idime,jdime)=worma(idime,jdime)&
                      +elcmm(itotv,ktotv)*rotma(kdime,jdime)
              end do
           end do
        end do
        do idime=1,ndime
           itotv=itot0+idime
           do jdime=1,ndime
              jtotv=jtot0+jdime
              elcmm(itotv,jtotv)=worma(idime,jdime)
           end do
        end do
     end do
  end if

end subroutine nsi_rotma3
