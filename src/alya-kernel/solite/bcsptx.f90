!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine bcsptx(nbnodes,nbvar,an,ja,ia,xx,yy) 

  !----------------------------------------------------------------------------------   
  ! Sources/kernel/solite/bcsptx.f90
  ! NAME 
  !     bcsptx
  ! DESCRIPTION
  !     This routine multiplies a transposed non-symmetric complex matrix stored in 
  !     BCSR format with a complex vector:
  !                         YY = A^T XX 
  ! INPUTPUT ARGUMENTS
  !    NBNODES .... Number of nodes
  !    NBVAR ...... Number of variables in each node
  !    AN ......... Complex matrix
  !    JA ......... Compressed Sparse format: index vector for column numbers
  !    IA ......... Compressed Sparse format: index vector for beginning of a row block
  !    XX ......... Complex vector
  ! OUTPUT ARGUMENTS
  !    YY ......... Complex result vector
  ! USES
  ! USED BY
  !-------------------------------------------------------------------------------------

  use def_kintyp, only : ip,rp
  use def_master, only : INOTMASTER,NPOIN_TYPE
 !Declaration statements
  implicit none

 !Dummy arguments
  integer(ip), intent(in)          :: nbnodes,nbvar
  complex(rp), intent(in)          :: an(*)
  integer(ip), intent(in)          :: ja(*),ia(*)
  complex(rp), intent(in)          :: xx(*)
  complex(rp), intent(out), target :: yy(*)

 !Local variables
  integer(ip)                      :: ii,ih,jj,jh,col,colh
!  integer(ip)                      :: kk,kh,ll
!  complex(rp)                      :: temp

 !------------------------------------------------------------
 ! y = A^Tx
 !------------------------------------------------------------
  if (INOTMASTER) then
     if (nbvar == 4_ip) then
        do ii = 1,nbnodes
           ih = 4_ip * (ii - 1_ip)
           yy(ih+1) = (0.0_rp,0.0_rp)
           yy(ih+2) = (0.0_rp,0.0_rp)
           yy(ih+3) = (0.0_rp,0.0_rp)
           yy(ih+4) = (0.0_rp,0.0_rp)
        end do
        do jj = 1,nbnodes
           jh = 4_ip * (jj - 1_ip)
           do ii = ia(jj),ia(jj+1)-1_ip
              ih = 16_ip * (ii - 1_ip)
              col = ja(ii)
              colh = 4_ip * (col - 1_ip)
              yy(colh+1) = yy(colh+1) + (an(ih+1)  * xx(jh+1) &
              +  an(ih+5)  * xx(jh+2) +  an(ih+9)  * xx(jh+3) &
              +  an(ih+13)  * xx(jh+4))
              yy(colh+2) = yy(colh+2) + (an(ih+2)  * xx(jh+1) &
              +  an(ih+6)  * xx(jh+2) +  an(ih+10)  * xx(jh+3) &
              +  an(ih+14)  * xx(jh+4))
              yy(colh+3) = yy(colh+3) + (an(ih+3)  * xx(jh+1) &
              +  an(ih+7) * xx(jh+2) +  an(ih+11) * xx(jh+3) &
              +  an(ih+15) * xx(jh+4))
              yy(colh+4) = yy(colh+4) + (an(ih+4) * xx(jh+1) &
              +  an(ih+8) * xx(jh+2) +  an(ih+12) * xx(jh+3) &
              +  an(ih+16) * xx(jh+4))                                   
           end do
        end do
     !else
     !   temp = (0.0_rp,0.0_rp)
     !   do ii = 1,nbnodes
     !      ih = nbvar * (ii - 1_ip)
     !      do kk = 1,nbvar
     !         kh = nbvar * (kk - 1_ip)
     !         yy(ih+kk) = (0.0_rp,0.0_rp)
     !         do jj = ia(ii),ia(ii+1)-1_ip
     !            jh = nbvar * nbvar * (jj - 1_ip)
     !            col = ja(jj)
     !            colh = nbvar * (col - 1_ip)
     !            do ll = 1,nbvar
     !               temp = temp + an(jh+kh+ll) * xx(colh+ll)
     !            end do
     !            yy(ih+kk) = yy(ih+kk) + temp
     !         end do
     !      end do
     !   end do
     end if
     call pararx('SLX',NPOIN_TYPE,nbnodes*nbvar,yy)
     !call rhsmox(nbvar,yy)
  end if

end subroutine bcsptx

